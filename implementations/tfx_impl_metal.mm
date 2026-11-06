#include "tfx_impl_metal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if __has_feature(objc_arc)
#error "Build tfx_impl_metal.mm without ARC (-fno-objc-arc): the renderer is a plain struct that releases what it owns"
#endif

#define TFX_IMPL_METAL_INITIAL_INSTANCE_CAPACITY 10000
//tfx_GetGlobalColorRampLayerCount never goes past this
#define TFX_IMPL_METAL_MAX_COLOR_RAMP_LAYERS 256
//The thread group sizes in timelinefx.slang: tfx_group_size, tfx_ribbon_group_size and the playback passes'
#define TFX_IMPL_METAL_KERNEL_GROUP_SIZE 256
#define TFX_IMPL_METAL_RIBBON_GROUP_SIZE 256
#define TFX_IMPL_METAL_PLAYBACK_GROUP_SIZE 128
#define TFX_IMPL_METAL_RING_SIZE (4 * 1024 * 1024)
//Every upload starts on 256 bytes, which suits buffer offsets, texture copies and the kernel's sections
#define TFX_IMPL_METAL_UPLOAD_ALIGNMENT 256
//The [[buffer]] indexes -D TFX_METAL_SLOTS gives: vertex reads and compute writes from 0, vertex uniforms and compute reads from 10, the compute uniform at 20,
//and the lit fragments' view at fragment buffer 0
#define TFX_IMPL_METAL_VERTEX_UNIFORM_SLOT 10
#define TFX_IMPL_METAL_FRAGMENT_UNIFORM_SLOT 0
#define TFX_IMPL_METAL_COMPUTE_READ_SLOT 10
#define TFX_IMPL_METAL_COMPUTE_UNIFORM_SLOT 20
//simulate_particles, compiled without the defines, takes its push constants at [[buffer(0)]] and its gpuAddress table at [[buffer(1)]]
#define TFX_IMPL_METAL_KERNEL_PUSH_SLOT 0
#define TFX_IMPL_METAL_KERNEL_TABLE_SLOT 1
//TFX_METAL_BUFFER_COUNT in timelinefx.slang
#define TFX_IMPL_METAL_BUFFER_TABLE_COUNT 16

//The kernel's table slots, which its push constants' buffer indexes point at
enum {
	tfx_impl_metal_table_pool,
	tfx_impl_metal_table_instances,
	tfx_impl_metal_table_graphs,
	tfx_impl_metal_table_upload,
};

//timelinefx.slang mirrors these structs by hand, in its [Buffer_rows] and [Render_rows] sections
static_assert(sizeof(tfx_instance_t) == 48, "tfx_instance_t changed, update tfx_instance in timelinefx.slang");
static_assert(sizeof(tfx_gpu_image_data_t) == 48, "tfx_gpu_image_data_t changed, update tfx_image_data in timelinefx.slang");
static_assert(sizeof(tfx_gpu_particle_properties_t) == 32, "tfx_gpu_particle_properties_t changed, update tfx_emitter_properties in timelinefx.slang");
static_assert(sizeof(tfx_ribbon_vertex_t) == 32, "tfx_ribbon_vertex_t changed, update tfx_ribbon_vertex_row in timelinefx.slang");
static_assert(sizeof(tfx_gpu_graph_data_t) == 48, "tfx_gpu_graph_data_t changed, update tfx_gpu_graph_data in timelinefx.slang");
static_assert(sizeof(tfx_gpu_particle_push_t) == 64, "tfx_gpu_particle_push_t changed, update tfx_gpu_particle_push in timelinefx.slang");
static_assert(sizeof(tfx_ribbon_bucket_globals_t) == 112, "tfx_ribbon_bucket_globals_t changed, update tfx_ribbon_bucket_globals in timelinefx.slang");
static_assert(sizeof(tfx_sprite_data_push_t) == 48, "tfx_sprite_data_push_t changed, update tfx_sprite_data_push in timelinefx.slang");
static_assert(sizeof(tfx_ribbon_data_push_t) == 32, "tfx_ribbon_data_push_t changed, update tfx_ribbon_data_push in timelinefx.slang");
static_assert(sizeof(tfx_impl_uniforms_t) % 16 == 0, "Uniform blocks must be a multiple of 16 bytes");
static_assert(sizeof(tfx_impl_uniforms_t) <= 4096, "setVertexBytes takes 4 KB at most");

//Mirrors tfx_fixed_draw_constants in timelinefx.slang
typedef struct tfx_impl_metal_draw_constants_s {
	tfxU32 instance_offset;
	tfxU32 previous_index_offset;
	float ndc_offset_x;
	float ndc_offset_y;
} tfx_impl_metal_draw_constants_t;

//--- Metal_resources

static void tfx_impl_metal__release_buffer(tfx_impl_metal_buffer_t *buffer) {
	[buffer->buffer release];
	buffer->buffer = nil;
	buffer->capacity = 0;
}

//A command buffer keeps what it uses alive, so a buffer replaced while a frame is in flight is released at once
static bool tfx_impl_metal__create_buffer(tfx_impl_metal_t *renderer, size_t size_in_bytes, tfx_impl_metal_buffer_t *buffer) {
	tfx_impl_metal__release_buffer(buffer);
	size_in_bytes = (size_in_bytes + 15) & ~(size_t)15;
	buffer->buffer = [renderer->device newBufferWithLength:size_in_bytes options:MTLResourceStorageModePrivate];
	if (!buffer->buffer) {
		tfx_Print("Could not create a buffer of %zu bytes", size_in_bytes);
		return false;
	}
	buffer->capacity = size_in_bytes;
	return true;
}

//Grows by half again whenever it's too small, so the size settles after the first few frames. Never empty, so a shader always has a buffer bound
static bool tfx_impl_metal__ensure_buffer(tfx_impl_metal_t *renderer, size_t size_in_bytes, tfx_impl_metal_buffer_t *buffer) {
	size_in_bytes = size_in_bytes ? size_in_bytes : 16;
	if (buffer->buffer && buffer->capacity >= size_in_bytes) {
		return true;
	}
	return tfx_impl_metal__create_buffer(renderer, size_in_bytes + size_in_bytes / 2, buffer);
}

static bool tfx_impl_metal__retire(tfx_impl_metal_frame_slot_t *slot, id<MTLBuffer> buffer) {
	if (slot->retired_count == slot->retired_capacity) {
		tfxU32 capacity = slot->retired_capacity ? slot->retired_capacity * 2 : 8;
		id<MTLBuffer> *retired = (id<MTLBuffer> *)realloc(slot->retired, sizeof(id<MTLBuffer>) * capacity);
		if (!retired) {
			return false;
		}
		slot->retired = retired;
		slot->retired_capacity = capacity;
	}
	slot->retired[slot->retired_count++] = buffer;
	return true;
}

//Shared memory for this frame. Small uploads share the frame slot's ring, big ones get a buffer of their own
static unsigned char *tfx_impl_metal__allocate_upload(tfx_impl_metal_t *renderer, size_t size_in_bytes, tfx_impl_metal_upload_t *upload) {
	tfx_impl_metal_frame_slot_t *slot = &renderer->frame_slots[renderer->frame_index];
	size_t start = (slot->used + TFX_IMPL_METAL_UPLOAD_ALIGNMENT - 1) & ~(size_t)(TFX_IMPL_METAL_UPLOAD_ALIGNMENT - 1);
	if (!slot->ring || start + size_in_bytes > slot->capacity) {
		if (size_in_bytes > TFX_IMPL_METAL_RING_SIZE / 4) {
			id<MTLBuffer> dedicated = [renderer->device newBufferWithLength:size_in_bytes options:MTLResourceStorageModeShared];
			if (!dedicated || !tfx_impl_metal__retire(slot, dedicated)) {
				[dedicated release];
				tfx_Print("Could not create an upload buffer of %zu bytes", size_in_bytes);
				return NULL;
			}
			upload->buffer = dedicated;
			upload->offset = 0;
			return (unsigned char *)[dedicated contents];
		}
		//What this frame already wrote to the old ring stays where its copies read it
		if (slot->ring && !tfx_impl_metal__retire(slot, slot->ring)) {
			return NULL;
		}
		size_t capacity = slot->capacity ? slot->capacity * 2 : TFX_IMPL_METAL_RING_SIZE;
		slot->ring = [renderer->device newBufferWithLength:capacity options:MTLResourceStorageModeShared];
		if (!slot->ring) {
			slot->mapped = NULL;
			slot->capacity = 0;
			tfx_Print("Could not create an upload ring of %zu bytes", capacity);
			return NULL;
		}
		slot->mapped = (unsigned char *)[slot->ring contents];
		slot->capacity = capacity;
		start = 0;
	}
	slot->used = start + size_in_bytes;
	upload->buffer = slot->ring;
	upload->offset = start;
	return slot->mapped + start;
}

//Copies into the frame's upload memory, padded to 16 bytes so a blit can round its size up and a pass can read it in place
static bool tfx_impl_metal__stage_upload(tfx_impl_metal_t *renderer, const void *data, size_t size_in_bytes, tfx_impl_metal_upload_t *upload) {
	unsigned char *mapped = tfx_impl_metal__allocate_upload(renderer, size_in_bytes ? (size_in_bytes + 15) & ~(size_t)15 : 16, upload);
	if (!mapped) {
		return false;
	}
	memcpy(mapped, data, size_in_bytes);
	return true;
}

static id<MTLBlitCommandEncoder> tfx_impl_metal__blit(id<MTLCommandBuffer> command_buffer, id<MTLBlitCommandEncoder> *blit) {
	if (!*blit) {
		*blit = [command_buffer blitCommandEncoder];
	}
	return *blit;
}

static void tfx_impl_metal__end_blit(id<MTLBlitCommandEncoder> blit) {
	if (blit) {
		[blit endEncoding];
	}
}

//Grows the private buffer when it's too small and copies the data into it through the frame's upload memory
static bool tfx_impl_metal__write_buffer(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, id<MTLBlitCommandEncoder> *blit, tfx_impl_metal_buffer_t *buffer,
	const void *data, size_t size_in_bytes) {
	if (!tfx_impl_metal__ensure_buffer(renderer, size_in_bytes, buffer)) {
		return false;
	}
	if (size_in_bytes == 0) {
		return true;
	}
	tfx_impl_metal_upload_t upload;
	if (!tfx_impl_metal__stage_upload(renderer, data, size_in_bytes, &upload)) {
		return false;
	}
	[tfx_impl_metal__blit(command_buffer, blit) copyFromBuffer:upload.buffer sourceOffset:upload.offset toBuffer:buffer->buffer destinationOffset:0 size:(size_in_bytes + 3) & ~(size_t)3];
	return true;
}

static id<MTLTexture> tfx_impl_metal__create_texture_array(tfx_impl_metal_t *renderer, tfxU32 width, tfxU32 height, tfxU32 layer_count, tfxU32 mip_count, MTLPixelFormat format) {
	MTLTextureDescriptor *descriptor = [MTLTextureDescriptor new];
	descriptor.textureType = MTLTextureType2DArray;
	descriptor.pixelFormat = format;
	descriptor.width = width;
	descriptor.height = height;
	descriptor.arrayLength = layer_count;
	descriptor.mipmapLevelCount = mip_count;
	descriptor.storageMode = MTLStorageModePrivate;
	descriptor.usage = MTLTextureUsageShaderRead;
	id<MTLTexture> texture = [renderer->device newTextureWithDescriptor:descriptor];
	[descriptor release];
	if (!texture) {
		tfx_Print("Could not create a %u x %u texture array of %u layers", width, height, layer_count);
	}
	return texture;
}

static void tfx_impl_metal__copy_to_texture(id<MTLBlitCommandEncoder> blit, const tfx_impl_metal_upload_t *upload, size_t offset, id<MTLTexture> texture,
	tfxU32 layer, tfxU32 mip_level, tfxU32 width, tfxU32 height, tfxU32 bytes_per_pixel) {
	[blit copyFromBuffer:upload->buffer sourceOffset:upload->offset + offset sourceBytesPerRow:(NSUInteger)width * bytes_per_pixel sourceBytesPerImage:(NSUInteger)width * height * bytes_per_pixel
		sourceSize:MTLSizeMake(width, height, 1) toTexture:texture destinationSlice:layer destinationLevel:mip_level destinationOrigin:MTLOriginMake(0, 0, 0)];
}

static MTLPixelFormat tfx_impl_metal__color_ramp_format(tfxU32 bytes_per_pixel) {
	if (bytes_per_pixel == 4) {
		return MTLPixelFormatRGBA8Unorm;
	} else if (bytes_per_pixel == 16) {
		return MTLPixelFormatRGBA32Float;
	}
	return MTLPixelFormatRGBA16Float;
}

//A 2 x 2 box filter, the same as the mips D3D11 generates
static void tfx_impl_metal__downsample(const unsigned char *source, tfxU32 source_size, unsigned char *destination) {
	tfxU32 size = source_size > 1 ? source_size / 2 : 1;
	tfxU32 last = source_size - 1;
	for (tfxU32 y = 0; y != size; ++y) {
		const unsigned char *top = source + (size_t)(y * 2 < last ? y * 2 : last) * source_size * 4;
		const unsigned char *bottom = source + (size_t)(y * 2 + 1 < last ? y * 2 + 1 : last) * source_size * 4;
		for (tfxU32 x = 0; x != size; ++x) {
			tfxU32 left = (x * 2 < last ? x * 2 : last) * 4;
			tfxU32 right = (x * 2 + 1 < last ? x * 2 + 1 : last) * 4;
			for (tfxU32 channel = 0; channel != 4; ++channel) {
				tfxU32 sum = top[left + channel] + top[right + channel] + bottom[left + channel] + bottom[right + channel];
				*destination++ = (unsigned char)((sum + 2) / 4);
			}
		}
	}
}

//slangc compiles every entry point ahead of time to <entry point>.metal in the shader folder, and Metal compiles the text here
static id<MTLFunction> tfx_impl_metal__load_function(id<MTLDevice> device, const char *shader_folder, const char *entry_point) {
	char path[1024];
	snprintf(path, sizeof(path), "%s/%s.metal", shader_folder, entry_point);
	FILE *file = fopen(path, "rb");
	if (!file) {
		tfx_Print("Could not open the shader %s", path);
		return nil;
	}
	fseek(file, 0, SEEK_END);
	long size = ftell(file);
	fseek(file, 0, SEEK_SET);
	char *code = size > 0 ? (char *)malloc((size_t)size) : NULL;
	bool loaded = code && fread(code, 1, (size_t)size, file) == (size_t)size;
	fclose(file);
	if (!loaded) {
		tfx_Print("Could not read the shader %s", path);
		free(code);
		return nil;
	}
	NSString *source = [[NSString alloc] initWithBytes:code length:(NSUInteger)size encoding:NSUTF8StringEncoding];
	free(code);
	MTLCompileOptions *options = [MTLCompileOptions new];
	options.languageVersion = MTLLanguageVersion3_0;
	NSError *error = nil;
	id<MTLLibrary> library = [device newLibraryWithSource:source options:options error:&error];
	[options release];
	[source release];
	if (!library) {
		tfx_Print("Could not compile the shader %s: %s", path, error ? [[error localizedDescription] UTF8String] : "");
		return nil;
	}
	NSString *name = [[NSString alloc] initWithUTF8String:entry_point];
	id<MTLFunction> function = [library newFunctionWithName:name];
	[name release];
	[library release];
	if (!function) {
		tfx_Print("The shader %s has no function %s", path, entry_point);
	}
	return function;
}

//The kernel and the ribbon compute need a full thread group of 256, which a heavy kernel can be compiled below
static id<MTLComputePipelineState> tfx_impl_metal__create_compute_pipeline(id<MTLDevice> device, const char *shader_folder, const char *entry_point, tfxU32 group_size) {
	id<MTLFunction> function = tfx_impl_metal__load_function(device, shader_folder, entry_point);
	if (!function) {
		return nil;
	}
	NSError *error = nil;
	id<MTLComputePipelineState> pipeline = [device newComputePipelineStateWithFunction:function error:&error];
	[function release];
	if (!pipeline) {
		tfx_Print("Could not create the %s pipeline: %s", entry_point, error ? [[error localizedDescription] UTF8String] : "");
		return nil;
	}
	if ([pipeline maxTotalThreadsPerThreadgroup] < group_size) {
		tfx_Print("The %s pipeline takes %u threads per group at most, it needs %u", entry_point, (tfxU32)[pipeline maxTotalThreadsPerThreadgroup], group_size);
		[pipeline release];
		return nil;
	}
	return pipeline;
}

//Premultiplied blending, or none for a layer that writes depth, and no vertex descriptor: instances and ribbon vertexes are read from buffers
static id<MTLRenderPipelineState> tfx_impl_metal__create_render_pipeline(id<MTLDevice> device, id<MTLFunction> vertex_function, id<MTLFunction> fragment_function,
	MTLPixelFormat color_format, MTLPixelFormat depth_format, bool writes_depth) {
	if (!vertex_function || !fragment_function) {
		return nil;
	}
	MTLRenderPipelineDescriptor *descriptor = [MTLRenderPipelineDescriptor new];
	descriptor.vertexFunction = vertex_function;
	descriptor.fragmentFunction = fragment_function;
	descriptor.depthAttachmentPixelFormat = depth_format;
	MTLRenderPipelineColorAttachmentDescriptor *color = descriptor.colorAttachments[0];
	color.pixelFormat = color_format;
	color.blendingEnabled = writes_depth ? NO : YES;
	color.rgbBlendOperation = MTLBlendOperationAdd;
	color.alphaBlendOperation = MTLBlendOperationAdd;
	color.sourceRGBBlendFactor = MTLBlendFactorOne;
	color.sourceAlphaBlendFactor = MTLBlendFactorOne;
	color.destinationRGBBlendFactor = MTLBlendFactorOneMinusSourceAlpha;
	color.destinationAlphaBlendFactor = MTLBlendFactorOneMinusSourceAlpha;
	NSError *error = nil;
	id<MTLRenderPipelineState> pipeline = [device newRenderPipelineStateWithDescriptor:descriptor error:&error];
	if (!pipeline) {
		tfx_Print("Could not create the %s pipeline: %s", [[vertex_function name] UTF8String], error ? [[error localizedDescription] UTF8String] : "");
	}
	[descriptor release];
	return pipeline;
}

//--- Metal_functions

bool tfx_impl_metal_Init(tfx_impl_metal_t *renderer, id<MTLDevice> device, MTLPixelFormat color_format, MTLPixelFormat depth_format, const char *shader_folder, tfxU32 frames_in_flight) {
	memset(renderer, 0, sizeof(tfx_impl_metal_t));
	//The kernel reaches its buffers through gpuAddress
	if (![device supportsFamily:MTLGPUFamilyMetal3]) {
		tfx_Print("The TimelineFX Metal backend needs a Metal 3 device");
		return false;
	}
	renderer->device = [device retain];
	renderer->frames_in_flight = frames_in_flight < 1 ? 1 : frames_in_flight > TFX_IMPL_METAL_MAX_FRAMES_IN_FLIGHT ? TFX_IMPL_METAL_MAX_FRAMES_IN_FLIGHT : frames_in_flight;

	id<MTLFunction> billboard_vertex_functions[2];
	id<MTLFunction> ribbon_vertex_functions[2];
	id<MTLFunction> billboard_fragment_functions[2][2];
	id<MTLFunction> ribbon_fragment_functions[2][2];
	for (int lit = 0; lit != 2; ++lit) {
		billboard_vertex_functions[lit] = tfx_impl_metal__load_function(device, shader_folder, tfx_impl_billboard_vertex_entry_points[lit]);
		ribbon_vertex_functions[lit] = tfx_impl_metal__load_function(device, shader_folder, tfx_impl_ribbon_vertex_entry_points[lit]);
		for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
			billboard_fragment_functions[writes_depth][lit] = tfx_impl_metal__load_function(device, shader_folder, tfx_impl_billboard_fragment_entry_points[writes_depth][lit]);
			ribbon_fragment_functions[writes_depth][lit] = tfx_impl_metal__load_function(device, shader_folder, tfx_impl_ribbon_fragment_entry_points[writes_depth][lit]);
		}
	}
	id<MTLFunction> baked_billboard_vertex = tfx_impl_metal__load_function(device, shader_folder, "tfx_billboard_baked_vertex");
	bool pipelines_created = true;
	for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
		for (int lit = 0; lit != 2; ++lit) {
			renderer->billboard_pipelines[writes_depth][lit] = tfx_impl_metal__create_render_pipeline(device, billboard_vertex_functions[lit], billboard_fragment_functions[writes_depth][lit],
				color_format, depth_format, writes_depth != 0);
			renderer->ribbon_pipelines[writes_depth][lit] = tfx_impl_metal__create_render_pipeline(device, ribbon_vertex_functions[lit], ribbon_fragment_functions[writes_depth][lit],
				color_format, depth_format, writes_depth != 0);
			pipelines_created = pipelines_created && renderer->billboard_pipelines[writes_depth][lit] && renderer->ribbon_pipelines[writes_depth][lit];
		}
	}
	renderer->baked_billboard_pipeline = tfx_impl_metal__create_render_pipeline(device, baked_billboard_vertex, billboard_fragment_functions[0][0], color_format, depth_format, false);
	[baked_billboard_vertex release];
	for (int lit = 0; lit != 2; ++lit) {
		[billboard_vertex_functions[lit] release];
		[ribbon_vertex_functions[lit] release];
		for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
			[billboard_fragment_functions[writes_depth][lit] release];
			[ribbon_fragment_functions[writes_depth][lit] release];
		}
	}
	renderer->particle_kernel = tfx_impl_metal__create_compute_pipeline(device, shader_folder, "simulate_particles", TFX_IMPL_METAL_KERNEL_GROUP_SIZE);
	renderer->ribbon_compute_pipeline = tfx_impl_metal__create_compute_pipeline(device, shader_folder, "tfx_ribbon_compute", TFX_IMPL_METAL_RIBBON_GROUP_SIZE);
	renderer->sprite_playback_pipeline = tfx_impl_metal__create_compute_pipeline(device, shader_folder, "tfx_sprite_playback", TFX_IMPL_METAL_PLAYBACK_GROUP_SIZE);
	renderer->ribbon_playback_pipeline = tfx_impl_metal__create_compute_pipeline(device, shader_folder, "tfx_ribbon_playback", TFX_IMPL_METAL_PLAYBACK_GROUP_SIZE);

	MTLDepthStencilDescriptor *depth_descriptor = [MTLDepthStencilDescriptor new];
	depth_descriptor.depthCompareFunction = MTLCompareFunctionLessEqual;
	depth_descriptor.depthWriteEnabled = NO;
	renderer->depth_state = [device newDepthStencilStateWithDescriptor:depth_descriptor];
	depth_descriptor.depthWriteEnabled = YES;
	renderer->depth_write_state = [device newDepthStencilStateWithDescriptor:depth_descriptor];
	[depth_descriptor release];

	MTLSamplerDescriptor *sampler_descriptor = [MTLSamplerDescriptor new];
	sampler_descriptor.minFilter = MTLSamplerMinMagFilterLinear;
	sampler_descriptor.magFilter = MTLSamplerMinMagFilterLinear;
	sampler_descriptor.mipFilter = MTLSamplerMipFilterLinear;
	sampler_descriptor.sAddressMode = MTLSamplerAddressModeClampToEdge;
	sampler_descriptor.tAddressMode = MTLSamplerAddressModeClampToEdge;
	renderer->atlas_sampler = [device newSamplerStateWithDescriptor:sampler_descriptor];
	[sampler_descriptor release];

	if (!pipelines_created || !renderer->baked_billboard_pipeline || !renderer->particle_kernel || !renderer->ribbon_compute_pipeline
		|| !renderer->sprite_playback_pipeline || !renderer->ribbon_playback_pipeline || !renderer->depth_state || !renderer->depth_write_state || !renderer->atlas_sampler) {
		tfx_Print("Could not create the TimelineFX pipelines and samplers");
		return false;
	}
	return true;
}

void tfx_impl_metal_BeginFrame(tfx_impl_metal_t *renderer, tfxU32 frame_index) {
	renderer->frame_index = frame_index % renderer->frames_in_flight;
	tfx_impl_metal_frame_slot_t *slot = &renderer->frame_slots[renderer->frame_index];
	for (tfxU32 index = 0; index != slot->retired_count; ++index) {
		[slot->retired[index] release];
	}
	slot->retired_count = 0;
	slot->used = 0;
}

bool tfx_impl_metal_CreateAtlas(tfx_impl_metal_t *renderer, tfxU32 layer_size, tfxU32 layer_count) {
	tfxU32 mip_count = 1;
	while ((layer_size >> mip_count) > 0) {
		mip_count++;
	}
	[renderer->shape_atlas release];
	renderer->shape_atlas = tfx_impl_metal__create_texture_array(renderer, layer_size, layer_size, layer_count, mip_count, MTLPixelFormatRGBA8Unorm);
	if (!renderer->shape_atlas) {
		return false;
	}
	renderer->shape_atlas_layer_size = layer_size;
	renderer->shape_atlas_layer_count = layer_count;
	renderer->shape_atlas_mip_count = mip_count;
	return true;
}

bool tfx_impl_metal_UploadAtlas(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, const tfx_impl_atlas_t *atlas) {
	if (!renderer->shape_atlas || atlas->layer_size != renderer->shape_atlas_layer_size
		|| atlas->first_layer + atlas->layer_count > renderer->shape_atlas_layer_count) {
		tfx_Print("An atlas of %u layers from layer %u, %u pixels square, doesn't fit the %u layer, %u pixel texture array",
			atlas->layer_count, atlas->first_layer, atlas->layer_size, renderer->shape_atlas_layer_count, renderer->shape_atlas_layer_size);
		return false;
	}
	//Every level of every layer in one upload, the mips downsampled straight into it
	size_t chain_bytes = 0;
	for (tfxU32 mip_level = 0; mip_level != renderer->shape_atlas_mip_count; ++mip_level) {
		size_t size = atlas->layer_size >> mip_level;
		chain_bytes += size * size * 4;
	}
	tfx_impl_metal_upload_t upload;
	unsigned char *mapped = tfx_impl_metal__allocate_upload(renderer, chain_bytes * atlas->layer_count, &upload);
	if (!mapped) {
		return false;
	}
	id<MTLBlitCommandEncoder> blit = [command_buffer blitCommandEncoder];
	size_t layer_bytes = (size_t)atlas->layer_size * atlas->layer_size * 4;
	size_t offset = 0;
	for (tfxU32 layer = 0; layer != atlas->layer_count; ++layer) {
		memcpy(mapped + offset, atlas->pixels + layer * layer_bytes, layer_bytes);
		const unsigned char *previous_level = NULL;
		tfxU32 previous_size = 0;
		tfxU32 size = atlas->layer_size;
		for (tfxU32 mip_level = 0; mip_level != renderer->shape_atlas_mip_count; ++mip_level) {
			if (previous_level) {
				tfx_impl_metal__downsample(previous_level, previous_size, mapped + offset);
			}
			tfx_impl_metal__copy_to_texture(blit, &upload, offset, renderer->shape_atlas, atlas->first_layer + layer, mip_level, size, size, 4);
			previous_level = mapped + offset;
			previous_size = size;
			offset += (size_t)size * size * 4;
			size = size > 1 ? size / 2 : 1;
		}
	}
	[blit endEncoding];
	return true;
}

bool tfx_impl_metal_UploadSpriteData(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, tfx_animation_manager animation_manager, tfx_gpu_shapes shapes) {
	id<MTLBlitCommandEncoder> blit = nil;
	bool uploaded = true;
	tfxU32 bitmap_count = tfx_GetAnimationColorRampBitmapCount(animation_manager);
	[renderer->baked_color_ramps release];
	renderer->baked_color_ramps = nil;
	if (bitmap_count > 0) {
		tfx_bitmap_t *first_bitmap = tfx_GetAnimationColorRampBitmap(animation_manager, 0);
		tfxU32 width = (tfxU32)tfx_GetBitmapWidth(first_bitmap);
		tfxU32 height = (tfxU32)tfx_GetBitmapHeight(first_bitmap);
		tfxU32 bytes_per_pixel = (tfxU32)(tfx_GetBitmapSize(first_bitmap) / ((size_t)width * height));
		renderer->baked_color_ramps = tfx_impl_metal__create_texture_array(renderer, width, height, bitmap_count, 1, tfx_impl_metal__color_ramp_format(bytes_per_pixel));
		uploaded = renderer->baked_color_ramps != nil;
		for (tfxU32 index = 0; index != bitmap_count && uploaded; ++index) {
			tfx_impl_metal_upload_t upload;
			uploaded = tfx_impl_metal__stage_upload(renderer, tfx_GetBitmapData(tfx_GetAnimationColorRampBitmap(animation_manager, index)), (size_t)width * height * bytes_per_pixel, &upload);
			if (uploaded) {
				tfx_impl_metal__copy_to_texture(tfx_impl_metal__blit(command_buffer, &blit), &upload, 0, renderer->baked_color_ramps, index, 0, width, height, bytes_per_pixel);
			}
		}
	}

	uploaded = uploaded
		&& tfx_impl_metal__write_buffer(renderer, command_buffer, &blit, &renderer->baked_image_data, tfx_GetGPUShapesArray(shapes), sizeof(tfx_gpu_image_data_t) * tfx_GetGPUShapesCount(shapes))
		&& tfx_impl_metal__write_buffer(renderer, command_buffer, &blit, &renderer->baked_emitter_properties, tfx_GetAnimationEmitterPropertiesBufferPointer(animation_manager), tfx_GetAnimationEmitterPropertySizeInBytes(animation_manager))
		&& tfx_impl_metal__write_buffer(renderer, command_buffer, &blit, &renderer->baked_sprites, tfx_GetSpriteDataBufferPointer(animation_manager), tfx_GetSpriteDataSizeInBytes(animation_manager));
	if (uploaded && tfx_AnimationManagerHasRibbons(animation_manager)) {
		//The baked ribbons, their properties, segments, emitters and graph lookups only change when sprite data is loaded
		size_t emitters_bytes = tfx_GetAnimationRibbonEmittersSizeInBytes(animation_manager);
		void *emitters = malloc(emitters_bytes ? emitters_bytes : 1);
		uploaded = emitters != NULL;
		if (uploaded) {
			tfx_CopyAnimationRibbonEmitters(animation_manager, emitters);
			uploaded = tfx_impl_metal__write_buffer(renderer, command_buffer, &blit, &renderer->baked_ribbons.graph_lookups, tfx_GetAnimationRibbonGraphLookupsBufferPointer(animation_manager),
					tfx_GetAnimationRibbonGraphLookupsSizeInBytes(animation_manager))
				&& tfx_impl_metal__write_buffer(renderer, command_buffer, &blit, &renderer->baked_ribbon_data, tfx_GetAnimationRibbonDataBufferPointer(animation_manager),
					tfx_GetAnimationRibbonDataSizeInBytes(animation_manager))
				&& tfx_impl_metal__write_buffer(renderer, command_buffer, &blit, &renderer->baked_ribbon_properties, tfx_GetAnimationRibbonPropertiesBufferPointer(animation_manager),
					tfx_GetAnimationRibbonPropertySizeInBytes(animation_manager))
				&& tfx_impl_metal__write_buffer(renderer, command_buffer, &blit, &renderer->baked_ribbons.segments, tfx_GetAnimationRibbonSegmentsBufferPointer(animation_manager),
					tfx_GetAnimationRibbonSegmentsSizeInBytes(animation_manager))
				&& tfx_impl_metal__write_buffer(renderer, command_buffer, &blit, &renderer->baked_ribbons.emitters, emitters, emitters_bytes);
			free(emitters);
		}
	}
	tfx_impl_metal__end_blit(blit);
	return uploaded;
}

//The vertex and index buffers ribbon_compute writes, sized for this frame's buckets
static bool tfx_impl_metal__ensure_ribbon_geometry(tfx_impl_metal_t *renderer, const tfx_ribbon_batch_sizes_t *sizes, tfx_impl_metal_ribbon_set_t *set) {
	return tfx_impl_metal__ensure_buffer(renderer, sizes->vertex_buffer_size_in_bytes, &set->vertices)
		&& tfx_impl_metal__ensure_buffer(renderer, sizes->index_buffer_size_in_bytes, &set->indices);
}

static bool tfx_impl_metal__upload_baked_ribbons(tfx_impl_metal_t *renderer, const tfx_impl_baked_frame_t *baked) {
	if (baked->ribbon_count == 0 || !renderer->baked_ribbon_data.buffer) {
		return true;
	}
	tfx_ribbon_batch_sizes_t sizes = tfx_GetAnimationRibbonBufferSizes(baked->animation_manager, sizeof(tfx_ribbon_vertex_t));
	//ribbon_playback writes the ribbons, then ribbon_compute and the ribbon draw read them
	if (!tfx_impl_metal__stage_upload(renderer, baked->ribbon_offsets, baked->ribbon_offsets_bytes, &renderer->animation_ribbon_offsets)
		|| !tfx_impl_metal__ensure_buffer(renderer, (size_t)baked->ribbon_count * tfxRIBBON_SIZE, &renderer->baked_ribbons.instances)
		|| !tfx_impl_metal__ensure_ribbon_geometry(renderer, &sizes, &renderer->baked_ribbons)) {
		return false;
	}
	memset(&renderer->ribbon_data_push, 0, sizeof(tfx_ribbon_data_push_t));
	renderer->ribbon_data_push.animation_instances_total = baked->animation_instance_count;
	renderer->ribbon_data_push.ribbons_total = baked->ribbon_count;
	renderer->playback_ribbon_count = baked->ribbon_count;
	renderer->baked_ribbons.animation_manager = baked->animation_manager;
	return true;
}

//Animation instances advance every frame, so this runs every frame, not just on updates. Each pass reads them in place from the frame's upload memory
static bool tfx_impl_metal__upload_baked(tfx_impl_metal_t *renderer, const tfx_impl_baked_frame_t *baked) {
	renderer->playback_sprite_count = 0;
	renderer->playback_ribbon_count = 0;
	renderer->baked_ribbons.animation_manager = NULL;
	if (!renderer->baked_sprites.buffer || (baked->sprite_count == 0 && baked->ribbon_count == 0)) {
		return true;
	}
	if (!tfx_impl_metal__stage_upload(renderer, baked->animation_instances, baked->animation_instances_bytes, &renderer->animation_instances)) {
		return false;
	}
	if (baked->sprite_count > 0) {
		if (!tfx_impl_metal__stage_upload(renderer, baked->offsets, baked->offsets_bytes, &renderer->animation_offsets)
			|| !tfx_impl_metal__ensure_buffer(renderer, (size_t)baked->sprite_count * sizeof(tfx_instance_t), &renderer->playback_instances)) {
			return false;
		}
		memset(&renderer->sprite_data_push, 0, sizeof(tfx_sprite_data_push_t));
		renderer->sprite_data_push.animation_instances_total = baked->animation_instance_count;
		renderer->sprite_data_push.billboards_total = baked->sprite_count;
		renderer->sprite_data_push.animated_shapes = baked->animated_shapes ? 1 : 0;
		renderer->playback_sprite_count = baked->sprite_count;
	}
	return tfx_impl_metal__upload_baked_ribbons(renderer, baked);
}

static bool tfx_impl_metal__upload_ribbons(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, id<MTLBlitCommandEncoder> *blit, const tfx_impl_ribbon_frame_t *ribbons) {
	tfx_impl_metal_ribbon_set_t *set = &renderer->live_ribbons;
	set->stage = ribbons->stage;
	if (!ribbons->stage) {
		return true;
	}
	const tfx_ribbon_batch_sizes_t *sizes = &ribbons->sizes;
	if (!tfx_impl_metal__write_buffer(renderer, command_buffer, blit, &set->segments, ribbons->segments, sizes->segment_buffer_size_in_bytes)
		|| !tfx_impl_metal__write_buffer(renderer, command_buffer, blit, &set->instances, ribbons->ribbons, sizes->ribbon_buffer_size_in_bytes)
		|| !tfx_impl_metal__write_buffer(renderer, command_buffer, blit, &set->emitters, ribbons->emitters, sizes->emitter_buffer_size_in_bytes)
		|| !tfx_impl_metal__ensure_ribbon_geometry(renderer, sizes, set)) {
		set->stage = NULL;
		return false;
	}
	set->stage_offsets = ribbons->stage_offsets;
	return true;
}

//The old contents are copied across, because the previous buffer is still what captured_index points into. In a blit pass of its own, ahead of the
//frame's uploads, so the instance upload lands on top of the copy
static bool tfx_impl_metal__grow_instance_buffers(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, tfxU32 needed) {
	tfxU32 capacity = renderer->instance_capacity ? renderer->instance_capacity : TFX_IMPL_METAL_INITIAL_INSTANCE_CAPACITY;
	while (capacity < needed) {
		capacity += capacity / 2;
	}
	id<MTLBlitCommandEncoder> blit = nil;
	bool grown = true;
	for (int index = 0; index != 2 && grown; ++index) {
		tfx_impl_metal_buffer_t old = renderer->instance_buffers[index];
		renderer->instance_buffers[index].buffer = nil;
		grown = tfx_impl_metal__create_buffer(renderer, (size_t)capacity * sizeof(tfx_instance_t), &renderer->instance_buffers[index]);
		if (old.buffer) {
			if (grown) {
				[tfx_impl_metal__blit(command_buffer, &blit) copyFromBuffer:old.buffer sourceOffset:0 toBuffer:renderer->instance_buffers[index].buffer destinationOffset:0 size:old.capacity];
			}
			[old.buffer release];
		}
	}
	tfx_impl_metal__end_blit(blit);
	renderer->instance_capacity = grown ? capacity : 0;
	return grown;
}

//The kernel's upload is read where tfx_CopyGPUParticleUpload writes it, through the gpuAddress table
static bool tfx_impl_metal__upload_gpu_particles(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, id<MTLBlitCommandEncoder> *blit, const tfx_impl_frame_t *frame) {
	const tfx_gpu_particle_frame_t *gpu_frame = &frame->gpu_particles;
	//The pool holds the particles between frames, so it's cleared once when it's made and never again
	if (!renderer->particle_pool.buffer) {
		if (!tfx_impl_metal__create_buffer(renderer, frame->gpu_pool_sizes.pool_bytes, &renderer->particle_pool)) {
			return false;
		}
		[tfx_impl_metal__blit(command_buffer, blit) fillBuffer:renderer->particle_pool.buffer range:NSMakeRange(0, renderer->particle_pool.capacity) value:0];
	}
	unsigned char *mapped = tfx_impl_metal__allocate_upload(renderer, gpu_frame->upload_bytes ? gpu_frame->upload_bytes : 16, &renderer->particle_upload);
	if (!mapped) {
		return false;
	}
	tfx_CopyGPUParticleUpload(gpu_frame, mapped);

	//One table slot for the pool with the bases, one for the instances, one for the graphs and one for the whole upload with the bases
	tfx_gpu_particle_push_t *push = &renderer->particle_push_constants;
	memset(push, 0, sizeof(tfx_gpu_particle_push_t));
	push->hot_buffer_index = tfx_impl_metal_table_pool;
	push->cold_buffer_index = tfx_impl_metal_table_pool;
	push->noise_buffer_index = tfx_impl_metal_table_pool;
	push->instance_buffer_index = tfx_impl_metal_table_instances;
	push->graph_buffer_index = tfx_impl_metal_table_graphs;
	push->dispatch_buffer_index = tfx_impl_metal_table_upload;
	push->tick_buffer_index = tfx_impl_metal_table_upload;
	push->spawn_buffer_index = tfx_impl_metal_table_upload;
	push->location_buffer_index = tfx_impl_metal_table_upload;
	push->cold_base = frame->gpu_pool_sizes.cold_base;
	push->noise_base = frame->gpu_pool_sizes.noise_base;
	push->tick_base = gpu_frame->tick_base;
	push->spawn_base = gpu_frame->spawn_base;
	push->location_base = gpu_frame->location_base;
	renderer->particle_dispatches = gpu_frame->dispatches;
	renderer->particle_dispatch_count = gpu_frame->dispatch_count;
	return true;
}

//Global layers never move, so only the layers whose version changed go up, unless the texture array has to grow
static bool tfx_impl_metal__upload_color_ramps(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, id<MTLBlitCommandEncoder> *blit) {
	tfxU32 layer_count = tfx_GetGlobalColorRampLayerCount();
	tfx_bitmap_t *first_bitmap = NULL;
	for (tfxU32 layer = 0; layer != layer_count && !first_bitmap; ++layer) {
		first_bitmap = tfx_GetGlobalColorRampLayer(layer);
	}
	if (!first_bitmap) {
		return true;
	}
	tfxU32 width = (tfxU32)tfx_GetBitmapWidth(first_bitmap);
	tfxU32 height = (tfxU32)tfx_GetBitmapHeight(first_bitmap);
	tfxU32 bytes_per_pixel = (tfxU32)(tfx_GetBitmapSize(first_bitmap) / ((size_t)width * height));
	bool recreate = layer_count > renderer->color_ramp_layer_capacity;
	if (recreate) {
		//Headroom so that loading another library doesn't always recreate it
		tfxU32 capacity = layer_count + layer_count / 2 + 1;
		capacity = capacity < TFX_IMPL_METAL_MAX_COLOR_RAMP_LAYERS ? capacity : TFX_IMPL_METAL_MAX_COLOR_RAMP_LAYERS;
		tfxU32 *layer_versions = (tfxU32 *)realloc(renderer->color_ramp_layer_versions, sizeof(tfxU32) * capacity);
		if (!layer_versions) {
			return false;
		}
		renderer->color_ramp_layer_versions = layer_versions;
		renderer->color_ramp_layer_capacity = 0;
		[renderer->color_ramps release];
		renderer->color_ramps = tfx_impl_metal__create_texture_array(renderer, width, height, capacity, 1, tfx_impl_metal__color_ramp_format(bytes_per_pixel));
		if (!renderer->color_ramps) {
			return false;
		}
		renderer->color_ramp_layer_capacity = capacity;
	}
	for (tfxU32 layer = 0; layer != layer_count; ++layer) {
		tfxU32 version = tfx_GetGlobalColorRampLayerVersion(layer);
		if (!recreate && version == renderer->color_ramp_layer_versions[layer]) {
			continue;
		}
		renderer->color_ramp_layer_versions[layer] = version;
		//A free layer isn't drawn from, so whatever it holds can stay
		tfx_bitmap_t *bitmap = tfx_GetGlobalColorRampLayer(layer);
		if (bitmap) {
			tfx_impl_metal_upload_t upload;
			if (!tfx_impl_metal__stage_upload(renderer, tfx_GetBitmapData(bitmap), (size_t)width * height * bytes_per_pixel, &upload)) {
				return false;
			}
			tfx_impl_metal__copy_to_texture(tfx_impl_metal__blit(command_buffer, blit), &upload, 0, renderer->color_ramps, layer, 0, width, height, bytes_per_pixel);
		}
	}
	return true;
}

//Only on update frames: the instances drawn until then were written with the layout the last upload holds
static bool tfx_impl_metal__upload_library_data(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, id<MTLBlitCommandEncoder> *blit) {
	tfxU32 properties_version = tfx_GetGlobalParticlePropertiesVersion();
	if (!renderer->emitter_properties.buffer || properties_version != renderer->properties_version) {
		if (!tfx_impl_metal__write_buffer(renderer, command_buffer, blit, &renderer->emitter_properties, tfx_GetGlobalParticlePropertiesBuffer(), tfx_GetGlobalParticlePropertiesBufferSizeInBytes())) {
			return false;
		}
		renderer->properties_version = properties_version;
	}
	tfxU32 shapes_version = tfx_GetGlobalGPUShapesVersion();
	if (!renderer->image_data.buffer || shapes_version != renderer->shapes_version) {
		tfx_gpu_shapes shapes = tfx_GetGlobalGPUShapes();
		if (!tfx_impl_metal__write_buffer(renderer, command_buffer, blit, &renderer->image_data, tfx_GetGPUShapesArray(shapes), sizeof(tfx_gpu_image_data_t) * tfx_GetGPUShapesCount(shapes))) {
			return false;
		}
		renderer->shapes_version = shapes_version;
	}
	//Shared by every library already, read by the ribbon shaders and the GPU kernel
	tfxU32 graph_lookups_version = tfx_GetGPUGraphLookupsVersion();
	if (!renderer->graph_lookups.buffer || graph_lookups_version != renderer->graph_lookups_version) {
		if (!tfx_impl_metal__write_buffer(renderer, command_buffer, blit, &renderer->graph_lookups, tfx_GetGPUGraphLookupsBuffer(), tfx_GetGPUGraphLookupsBufferSizeInBytes())) {
			return false;
		}
		renderer->graph_lookups_version = graph_lookups_version;
	}
	tfxU32 color_ramps_version = tfx_GetGlobalColorRampsVersion();
	if (!renderer->color_ramps || color_ramps_version != renderer->color_ramps_version) {
		if (!tfx_impl_metal__upload_color_ramps(renderer, command_buffer, blit)) {
			return false;
		}
		renderer->color_ramps_version = color_ramps_version;
	}
	return true;
}

static bool tfx_impl_metal__upload_stage(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, id<MTLBlitCommandEncoder> *blit, const tfx_impl_frame_t *frame, tfxU32 needed, bool gpu_particles) {
	if (!tfx_impl_metal__upload_library_data(renderer, command_buffer, blit)) {
		return false;
	}
	//The other buffer now holds what was drawn last, which is what captured_index points into for CPU and GPU instances alike
	renderer->current_instance_buffer ^= 1;
	tfx_impl_metal_buffer_t *instances = &renderer->instance_buffers[renderer->current_instance_buffer];
	if (frame->instance_count > 0) {
		tfx_impl_metal_upload_t upload;
		if (!tfx_impl_metal__stage_upload(renderer, frame->instances, frame->instance_count * sizeof(tfx_instance_t), &upload)) {
			return false;
		}
		[tfx_impl_metal__blit(command_buffer, blit) copyFromBuffer:upload.buffer sourceOffset:upload.offset toBuffer:instances->buffer destinationOffset:0 size:frame->instance_count * sizeof(tfx_instance_t)];
	}
	renderer->instance_count = needed;
	memcpy(renderer->layer_start, frame->layer_start, sizeof(renderer->layer_start));
	memcpy(renderer->layer_count, frame->layer_count, sizeof(renderer->layer_count));
	memset(renderer->gpu_layer_count, 0, sizeof(renderer->gpu_layer_count));
	if (gpu_particles) {
		memcpy(renderer->gpu_layer_start, frame->gpu_particles.gpu_layer_start, sizeof(renderer->gpu_layer_start));
		memcpy(renderer->gpu_layer_count, frame->gpu_particles.gpu_layer_count, sizeof(renderer->gpu_layer_count));
		if (!tfx_impl_metal__upload_gpu_particles(renderer, command_buffer, blit, frame)) {
			return false;
		}
	}
	return tfx_impl_metal__upload_ribbons(renderer, command_buffer, blit, &frame->ribbons);
}

bool tfx_impl_metal_UploadFrame(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, const tfx_impl_frame_t *frame) {
	memcpy(renderer->layer_flags, frame->layer_flags, sizeof(renderer->layer_flags));
	if (!tfx_impl_metal__upload_baked(renderer, &frame->baked)) {
		return false;
	}
	if (!frame->stage_updated) {
		return true;
	}
	bool gpu_particles = frame->has_gpu_particles && renderer->particle_kernel;
	tfxU32 needed = gpu_particles ? frame->gpu_particles.gpu_instance_start + frame->gpu_particles.gpu_instance_count : frame->instance_count;
	if ((needed > renderer->instance_capacity || renderer->instance_capacity == 0) && !tfx_impl_metal__grow_instance_buffers(renderer, command_buffer, needed)) {
		renderer->instance_count = 0;
		return false;
	}
	id<MTLBlitCommandEncoder> blit = nil;
	bool uploaded = tfx_impl_metal__upload_stage(renderer, command_buffer, &blit, frame, needed, gpu_particles);
	tfx_impl_metal__end_blit(blit);
	return uploaded;
}

static bool tfx_impl_metal__ribbon_set_active(const tfx_impl_metal_ribbon_set_t *set) {
	return set->stage || set->animation_manager;
}

static tfx_ribbon_dispatch_t tfx_impl_metal__create_ribbon_dispatch(const tfx_impl_metal_ribbon_set_t *set) {
	return set->stage ? tfx_CreateRibbonBatchDispatch(&set->stage_offsets) : tfx_CreateAnimationRibbonDispatch();
}

static bool tfx_impl_metal__next_ribbon_dispatch(const tfx_impl_metal_ribbon_set_t *set, tfx_ribbon_dispatch_t *dispatch) {
	return set->stage ? tfx_NextRibbonDispatch(set->stage, dispatch) : tfx_NextAnimationRibbonDispatch(set->animation_manager, dispatch);
}

static id<MTLBuffer> tfx_impl_metal__ribbon_graph_lookups(const tfx_impl_metal_t *renderer, const tfx_impl_metal_ribbon_set_t *set) {
	return set->stage ? renderer->graph_lookups.buffer : set->graph_lookups.buffer;
}

static MTLSize tfx_impl_metal__groups(tfxU32 count, tfxU32 group_size) {
	return MTLSizeMake((count + group_size - 1) / group_size, 1, 1);
}

//Every fetch has to be dispatched, the next one carries on from the particle state this one leaves in the pool
static void tfx_impl_metal__simulate_gpu_particles(tfx_impl_metal_t *renderer, id<MTLComputeCommandEncoder> encoder) {
	id<MTLBuffer> instances = renderer->instance_buffers[renderer->current_instance_buffer].buffer;
	uint64_t table[TFX_IMPL_METAL_BUFFER_TABLE_COUNT] = { 0 };
	table[tfx_impl_metal_table_pool] = [renderer->particle_pool.buffer gpuAddress];
	table[tfx_impl_metal_table_instances] = [instances gpuAddress];
	table[tfx_impl_metal_table_graphs] = [renderer->graph_lookups.buffer gpuAddress];
	table[tfx_impl_metal_table_upload] = [renderer->particle_upload.buffer gpuAddress] + renderer->particle_upload.offset;
	//The kernel only reaches these through the table, so Metal has to be told what it reads and writes to keep them resident and ordered
	[encoder useResource:renderer->particle_pool.buffer usage:MTLResourceUsageRead | MTLResourceUsageWrite];
	[encoder useResource:instances usage:MTLResourceUsageRead | MTLResourceUsageWrite];
	[encoder useResource:renderer->graph_lookups.buffer usage:MTLResourceUsageRead];
	[encoder useResource:renderer->particle_upload.buffer usage:MTLResourceUsageRead];
	[encoder setComputePipelineState:renderer->particle_kernel];
	[encoder setBytes:table length:sizeof(table) atIndex:TFX_IMPL_METAL_KERNEL_TABLE_SLOT];
	tfx_gpu_particle_push_t *push = &renderer->particle_push_constants;
	for (tfxU32 index = 0; index != renderer->particle_dispatch_count; ++index) {
		tfxU32 slot_count = renderer->particle_dispatches[index].slot_count;
		if (slot_count == 0) {
			continue;
		}
		push->dispatch_index = index;
		[encoder setBytes:push length:sizeof(tfx_gpu_particle_push_t) atIndex:TFX_IMPL_METAL_KERNEL_PUSH_SLOT];
		[encoder dispatchThreadgroups:tfx_impl_metal__groups(slot_count, TFX_IMPL_METAL_KERNEL_GROUP_SIZE) threadsPerThreadgroup:MTLSizeMake(TFX_IMPL_METAL_KERNEL_GROUP_SIZE, 1, 1)];
	}
	renderer->particle_dispatch_count = 0;
}

//Both playback passes write slot 0 and read the offsets, the animation instances, the baked rows and their properties from 10
static void tfx_impl_metal__play_back(id<MTLComputeCommandEncoder> encoder, id<MTLComputePipelineState> pipeline, id<MTLBuffer> output, const tfx_impl_metal_upload_t *offsets,
	const tfx_impl_metal_upload_t *animation_instances, id<MTLBuffer> rows, id<MTLBuffer> properties, const void *push, size_t push_size, tfxU32 count) {
	[encoder setComputePipelineState:pipeline];
	[encoder setBuffer:output offset:0 atIndex:0];
	[encoder setBuffer:offsets->buffer offset:offsets->offset atIndex:TFX_IMPL_METAL_COMPUTE_READ_SLOT];
	[encoder setBuffer:animation_instances->buffer offset:animation_instances->offset atIndex:TFX_IMPL_METAL_COMPUTE_READ_SLOT + 1];
	[encoder setBuffer:rows offset:0 atIndex:TFX_IMPL_METAL_COMPUTE_READ_SLOT + 2];
	[encoder setBuffer:properties offset:0 atIndex:TFX_IMPL_METAL_COMPUTE_READ_SLOT + 3];
	[encoder setBytes:push length:push_size atIndex:TFX_IMPL_METAL_COMPUTE_UNIFORM_SLOT];
	[encoder dispatchThreadgroups:tfx_impl_metal__groups(count, TFX_IMPL_METAL_PLAYBACK_GROUP_SIZE) threadsPerThreadgroup:MTLSizeMake(TFX_IMPL_METAL_PLAYBACK_GROUP_SIZE, 1, 1)];
}

static void tfx_impl_metal__build_ribbon_geometry(tfx_impl_metal_t *renderer, id<MTLComputeCommandEncoder> encoder, tfx_impl_metal_ribbon_set_t *set, const tfx_impl_uniforms_t *uniforms) {
	[encoder setComputePipelineState:renderer->ribbon_compute_pipeline];
	[encoder setBuffer:set->vertices.buffer offset:0 atIndex:0];
	[encoder setBuffer:set->indices.buffer offset:0 atIndex:1];
	[encoder setBuffer:set->emitters.buffer offset:0 atIndex:TFX_IMPL_METAL_COMPUTE_READ_SLOT];
	[encoder setBuffer:tfx_impl_metal__ribbon_graph_lookups(renderer, set) offset:0 atIndex:TFX_IMPL_METAL_COMPUTE_READ_SLOT + 1];
	[encoder setBuffer:set->instances.buffer offset:0 atIndex:TFX_IMPL_METAL_COMPUTE_READ_SLOT + 2];
	[encoder setBuffer:set->segments.buffer offset:0 atIndex:TFX_IMPL_METAL_COMPUTE_READ_SLOT + 3];
	//One dispatch per bucket, each with the bucket's own globals. Buckets write their own ranges, so they can overlap
	tfx_ribbon_dispatch_t dispatch = tfx_impl_metal__create_ribbon_dispatch(set);
	while (tfx_impl_metal__next_ribbon_dispatch(set, &dispatch)) {
		tfx_ribbon_bucket_globals_t *globals = tfx_GetRibbonDispatchGlobals(&dispatch);
		if (set->stage) {
			globals->lerp = uniforms->timer_lerp;
		} else {
			//Baked ribbons are already interpolated, but they still face the camera
			globals->camera_position.x = uniforms->camera_position[0];
			globals->camera_position.y = uniforms->camera_position[1];
			globals->camera_position.z = uniforms->camera_position[2];
		}
		[encoder setBytes:globals length:sizeof(tfx_ribbon_bucket_globals_t) atIndex:TFX_IMPL_METAL_COMPUTE_UNIFORM_SLOT];
		[encoder dispatchThreadgroups:tfx_impl_metal__groups(dispatch.total_segments, TFX_IMPL_METAL_RIBBON_GROUP_SIZE) threadsPerThreadgroup:MTLSizeMake(TFX_IMPL_METAL_RIBBON_GROUP_SIZE, 1, 1)];
	}
}

void tfx_impl_metal_Compute(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, const tfx_impl_uniforms_t *uniforms) {
	bool live_ribbons = tfx_impl_metal__ribbon_set_active(&renderer->live_ribbons);
	bool baked_ribbons = tfx_impl_metal__ribbon_set_active(&renderer->baked_ribbons);
	if (renderer->particle_dispatch_count == 0 && renderer->playback_sprite_count == 0 && renderer->playback_ribbon_count == 0 && !live_ribbons && !baked_ribbons) {
		return;
	}
	//Concurrent, because the kernel's dispatches write their own rings and instance ranges and the playback passes their own outputs. Only
	//the baked ribbon geometry reads what an earlier dispatch wrote, hence the one barrier. Metal orders the passes before and after for us
	id<MTLComputeCommandEncoder> encoder = [command_buffer computeCommandEncoderWithDispatchType:MTLDispatchTypeConcurrent];
	if (renderer->particle_dispatch_count > 0) {
		tfx_impl_metal__simulate_gpu_particles(renderer, encoder);
	}
	if (renderer->playback_sprite_count > 0) {
		tfx_impl_metal__play_back(encoder, renderer->sprite_playback_pipeline, renderer->playback_instances.buffer, &renderer->animation_offsets, &renderer->animation_instances,
			renderer->baked_sprites.buffer, renderer->baked_emitter_properties.buffer, &renderer->sprite_data_push, sizeof(tfx_sprite_data_push_t), renderer->playback_sprite_count);
	}
	if (renderer->playback_ribbon_count > 0) {
		tfx_impl_metal__play_back(encoder, renderer->ribbon_playback_pipeline, renderer->baked_ribbons.instances.buffer, &renderer->animation_ribbon_offsets, &renderer->animation_instances,
			renderer->baked_ribbon_data.buffer, renderer->baked_ribbon_properties.buffer, &renderer->ribbon_data_push, sizeof(tfx_ribbon_data_push_t), renderer->playback_ribbon_count);
		[encoder memoryBarrierWithScope:MTLBarrierScopeBuffers];
	}
	if (live_ribbons) {
		tfx_impl_metal__build_ribbon_geometry(renderer, encoder, &renderer->live_ribbons, uniforms);
	}
	if (baked_ribbons) {
		tfx_impl_metal__build_ribbon_geometry(renderer, encoder, &renderer->baked_ribbons, uniforms);
	}
	[encoder endEncoding];
}

//Live and baked ribbons each have their own shape frames and color ramps, the atlas is shared
//Every bucket is walked, because the dispatch carries the offsets from one to the next, and only the ones whose layer is in the group are drawn
static void tfx_impl_metal__render_ribbons(tfx_impl_metal_t *renderer, id<MTLRenderCommandEncoder> render_encoder, tfx_impl_metal_ribbon_set_t *set, id<MTLBuffer> image_data, int writes_depth) {
	[render_encoder setVertexBuffer:set->vertices.buffer offset:0 atIndex:0];
	[render_encoder setVertexBuffer:set->instances.buffer offset:0 atIndex:1];
	[render_encoder setVertexBuffer:set->emitters.buffer offset:0 atIndex:2];
	[render_encoder setVertexBuffer:tfx_impl_metal__ribbon_graph_lookups(renderer, set) offset:0 atIndex:3];
	[render_encoder setVertexBuffer:image_data offset:0 atIndex:4];
	//The vertex shader doesn't read the lerp or the camera position, so the globals go as the library gives them
	id<MTLRenderPipelineState> bound_pipeline = nil;
	tfx_ribbon_dispatch_t dispatch = tfx_impl_metal__create_ribbon_dispatch(set);
	while (tfx_impl_metal__next_ribbon_dispatch(set, &dispatch)) {
		const tfxLayerFlags flags = renderer->layer_flags[dispatch.layer];
		if (dispatch.index_count == 0 || tfx_impl_LayerWritesDepth(flags) != writes_depth) {
			continue;
		}
		id<MTLRenderPipelineState> pipeline = renderer->ribbon_pipelines[writes_depth][tfx_impl_LayerIsLit(flags)];
		if (pipeline != bound_pipeline) {
			[render_encoder setRenderPipelineState:pipeline];
			bound_pipeline = pipeline;
		}
		[render_encoder setVertexBytes:tfx_GetRibbonDispatchGlobals(&dispatch) length:sizeof(tfx_ribbon_bucket_globals_t) atIndex:TFX_IMPL_METAL_VERTEX_UNIFORM_SLOT + 1];
		[render_encoder drawIndexedPrimitives:MTLPrimitiveTypeTriangle indexCount:dispatch.index_count indexType:MTLIndexTypeUInt32 indexBuffer:set->indices.buffer
			indexBufferOffset:(NSUInteger)dispatch.index_offset * sizeof(tfxU32)];
	}
}

//The first instance goes in the draw constants, like the other fixed binding backends, though Metal's instance_id would include a base instance
static void tfx_impl_metal__draw_instance_range(id<MTLRenderCommandEncoder> render_encoder, tfxU32 start, tfxU32 count) {
	if (count == 0) {
		return;
	}
	tfx_impl_metal_draw_constants_t draw_constants = { start, 0, 0.f, 0.f };
	[render_encoder setVertexBytes:&draw_constants length:sizeof(draw_constants) atIndex:TFX_IMPL_METAL_VERTEX_UNIFORM_SLOT + 1];
	[render_encoder drawPrimitives:MTLPrimitiveTypeTriangle vertexStart:0 vertexCount:6 instanceCount:count];
}

//Each layer in the group as its CPU instances then its GPU instances
static void tfx_impl_metal__render_billboards(tfx_impl_metal_t *renderer, id<MTLRenderCommandEncoder> render_encoder, int writes_depth) {
	tfxU32 current = renderer->current_instance_buffer;
	[render_encoder setVertexBuffer:renderer->instance_buffers[current].buffer offset:0 atIndex:0];
	[render_encoder setVertexBuffer:renderer->instance_buffers[current ^ 1].buffer offset:0 atIndex:1];
	[render_encoder setVertexBuffer:renderer->image_data.buffer offset:0 atIndex:2];
	[render_encoder setVertexBuffer:renderer->emitter_properties.buffer offset:0 atIndex:3];
	id<MTLRenderPipelineState> bound_pipeline = nil;
	for (tfxU32 layer = 0; layer != tfxLAYERS; ++layer) {
		const tfxLayerFlags flags = renderer->layer_flags[layer];
		if (tfx_impl_LayerWritesDepth(flags) != writes_depth || renderer->layer_count[layer] + renderer->gpu_layer_count[layer] == 0) {
			continue;
		}
		id<MTLRenderPipelineState> pipeline = renderer->billboard_pipelines[writes_depth][tfx_impl_LayerIsLit(flags)];
		if (pipeline != bound_pipeline) {
			[render_encoder setRenderPipelineState:pipeline];
			bound_pipeline = pipeline;
		}
		tfx_impl_metal__draw_instance_range(render_encoder, renderer->layer_start[layer], renderer->layer_count[layer]);
		tfx_impl_metal__draw_instance_range(render_encoder, renderer->gpu_layer_start[layer], renderer->gpu_layer_count[layer]);
	}
}

void tfx_impl_metal_Render(tfx_impl_metal_t *renderer, id<MTLRenderCommandEncoder> render_encoder, const tfx_impl_uniforms_t *uniforms) {
	if (!renderer->shape_atlas || (renderer->instance_count == 0 && renderer->playback_sprite_count == 0
		&& !tfx_impl_metal__ribbon_set_active(&renderer->live_ribbons) && !tfx_impl_metal__ribbon_set_active(&renderer->baked_ribbons))) {
		return;
	}
	[render_encoder setVertexBytes:uniforms length:sizeof(tfx_impl_uniforms_t) atIndex:TFX_IMPL_METAL_VERTEX_UNIFORM_SLOT];
	//The lit fragments read the light from the same uniforms
	[render_encoder setFragmentBytes:uniforms length:sizeof(tfx_impl_uniforms_t) atIndex:TFX_IMPL_METAL_FRAGMENT_UNIFORM_SLOT];
	[render_encoder setFragmentTexture:renderer->shape_atlas atIndex:0];
	[render_encoder setFragmentSamplerState:renderer->atlas_sampler atIndex:0];
	//The live draws need the library data, which goes up with the first update
	const bool live_data = renderer->color_ramps && renderer->image_data.buffer;
	//The layers that write depth draw first with blending off, so they hide the blended layers behind them
	for (int group = 0; group != 2; ++group) {
		const int writes_depth = group == 0;
		[render_encoder setDepthStencilState:writes_depth ? renderer->depth_write_state : renderer->depth_state];
		if (live_data) {
			[render_encoder setFragmentTexture:renderer->color_ramps atIndex:1];
			if (renderer->instance_count > 0) {
				tfx_impl_metal__render_billboards(renderer, render_encoder, writes_depth);
			}
			if (tfx_impl_metal__ribbon_set_active(&renderer->live_ribbons)) {
				tfx_impl_metal__render_ribbons(renderer, render_encoder, &renderer->live_ribbons, renderer->image_data.buffer, writes_depth);
			}
		}
		if (renderer->baked_color_ramps) {
			//Baked sprites come with their own shape frames, emitter properties and color ramps, the atlas is shared
			[render_encoder setFragmentTexture:renderer->baked_color_ramps atIndex:1];
			//Baked billboards come as one range across their layers, so they all draw blended and unlit
			if (!writes_depth && renderer->playback_sprite_count > 0) {
				[render_encoder setRenderPipelineState:renderer->baked_billboard_pipeline];
				[render_encoder setVertexBuffer:renderer->playback_instances.buffer offset:0 atIndex:0];
				[render_encoder setVertexBuffer:renderer->baked_image_data.buffer offset:0 atIndex:1];
				[render_encoder setVertexBuffer:renderer->baked_emitter_properties.buffer offset:0 atIndex:2];
				tfx_impl_metal__draw_instance_range(render_encoder, 0, renderer->playback_sprite_count);
			}
			if (tfx_impl_metal__ribbon_set_active(&renderer->baked_ribbons)) {
				tfx_impl_metal__render_ribbons(renderer, render_encoder, &renderer->baked_ribbons, renderer->baked_image_data.buffer, writes_depth);
			}
		}
	}
}

static void tfx_impl_metal__release_ribbon_set(tfx_impl_metal_ribbon_set_t *set) {
	tfx_impl_metal__release_buffer(&set->segments);
	tfx_impl_metal__release_buffer(&set->instances);
	tfx_impl_metal__release_buffer(&set->emitters);
	tfx_impl_metal__release_buffer(&set->vertices);
	tfx_impl_metal__release_buffer(&set->indices);
	tfx_impl_metal__release_buffer(&set->graph_lookups);
}

void tfx_impl_metal_Shutdown(tfx_impl_metal_t *renderer) {
	if (!renderer->device) {
		return;
	}
	for (tfxU32 frame = 0; frame != TFX_IMPL_METAL_MAX_FRAMES_IN_FLIGHT; ++frame) {
		tfx_impl_metal_frame_slot_t *slot = &renderer->frame_slots[frame];
		for (tfxU32 index = 0; index != slot->retired_count; ++index) {
			[slot->retired[index] release];
		}
		free(slot->retired);
		[slot->ring release];
	}
	for (int index = 0; index != 2; ++index) {
		tfx_impl_metal__release_buffer(&renderer->instance_buffers[index]);
	}
	tfx_impl_metal__release_buffer(&renderer->image_data);
	tfx_impl_metal__release_buffer(&renderer->emitter_properties);
	tfx_impl_metal__release_buffer(&renderer->graph_lookups);
	tfx_impl_metal__release_ribbon_set(&renderer->live_ribbons);
	tfx_impl_metal__release_ribbon_set(&renderer->baked_ribbons);
	tfx_impl_metal__release_buffer(&renderer->particle_pool);
	tfx_impl_metal__release_buffer(&renderer->baked_image_data);
	tfx_impl_metal__release_buffer(&renderer->baked_emitter_properties);
	tfx_impl_metal__release_buffer(&renderer->baked_sprites);
	tfx_impl_metal__release_buffer(&renderer->playback_instances);
	tfx_impl_metal__release_buffer(&renderer->baked_ribbon_data);
	tfx_impl_metal__release_buffer(&renderer->baked_ribbon_properties);
	[renderer->shape_atlas release];
	[renderer->color_ramps release];
	[renderer->baked_color_ramps release];
	[renderer->atlas_sampler release];
	[renderer->depth_state release];
	[renderer->depth_write_state release];
	for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
		for (int lit = 0; lit != 2; ++lit) {
			[renderer->billboard_pipelines[writes_depth][lit] release];
			[renderer->ribbon_pipelines[writes_depth][lit] release];
		}
	}
	[renderer->baked_billboard_pipeline release];
	[renderer->particle_kernel release];
	[renderer->ribbon_compute_pipeline release];
	[renderer->sprite_playback_pipeline release];
	[renderer->ribbon_playback_pipeline release];
	[renderer->device release];
	free(renderer->color_ramp_layer_versions);
	memset(renderer, 0, sizeof(tfx_impl_metal_t));
}
