#include "tfx_impl_webgpu.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TFX_IMPL_WEBGPU_INITIAL_INSTANCE_CAPACITY 10000
//tfx_GetGlobalColorRampLayerCount never goes past this
#define TFX_IMPL_WEBGPU_MAX_COLOR_RAMP_LAYERS 256
//Matches tfx_ribbon_group_size in timelinefx.slang
#define TFX_IMPL_WEBGPU_RIBBON_GROUP_SIZE 256
//Dynamic offsets have to be multiples of minUniformBufferOffsetAlignment, which is never more than 256
#define TFX_IMPL_WEBGPU_UNIFORM_BLOCK_SIZE 256
//Render's blocks: the view, then the draw constants of each layer's CPU and GPU range and of the baked billboards
#define TFX_IMPL_WEBGPU_RENDER_BLOCKS (2 + tfxLAYERS * 2)
#define TFX_IMPL_WEBGPU_INITIAL_UNIFORM_BLOCKS 64

//timelinefx.slang mirrors these structs by hand, in its [Render_rows] section
static_assert(sizeof(tfx_instance_t) == 48, "tfx_instance_t changed, update tfx_instance in timelinefx.slang");
static_assert(sizeof(tfx_gpu_image_data_t) == 48, "tfx_gpu_image_data_t changed, update tfx_image_data in timelinefx.slang");
static_assert(sizeof(tfx_gpu_particle_properties_t) == 32, "tfx_gpu_particle_properties_t changed, update tfx_emitter_properties in timelinefx.slang");
static_assert(sizeof(tfx_ribbon_vertex_t) == 32, "tfx_ribbon_vertex_t changed, update tfx_ribbon_vertex_row in timelinefx.slang");
static_assert(sizeof(tfx_gpu_graph_data_t) == 48, "tfx_gpu_graph_data_t changed, update tfx_gpu_graph_data in timelinefx.slang");
static_assert(sizeof(tfx_impl_uniforms_t) % 16 == 0, "Uniform buffers must be a multiple of 16 bytes");
static_assert(sizeof(tfx_impl_uniforms_t) <= TFX_IMPL_WEBGPU_UNIFORM_BLOCK_SIZE && sizeof(tfx_ribbon_bucket_globals_t) <= TFX_IMPL_WEBGPU_UNIFORM_BLOCK_SIZE,
	"Every uniform block has to fit in one step of the uniform buffer");

//Mirrors tfx_fixed_draw_constants in timelinefx.slang
typedef struct tfx_impl_webgpu_draw_constants_s {
	tfxU32 instance_offset;
	tfxU32 previous_index_offset;
	float ndc_offset_x;
	float ndc_offset_y;
} tfx_impl_webgpu_draw_constants_t;

//--- Webgpu_resources

static WGPUStringView tfx_impl_webgpu__string(const char *text) {
	WGPUStringView view = { text, WGPU_STRLEN };
	return view;
}

static void tfx_impl_webgpu__release_buffer(tfx_impl_webgpu_buffer_t *buffer) {
	if (buffer->buffer) {
		wgpuBufferRelease(buffer->buffer);
	}
	buffer->buffer = NULL;
	buffer->capacity = 0;
}

//WebGPU zeroes every new buffer. A buffer released here stays alive until the GPU and every bind group holding it are done with it
static bool tfx_impl_webgpu__create_buffer(tfx_impl_webgpu_t *renderer, tfxU32 size_in_bytes, WGPUBufferUsage usage, tfx_impl_webgpu_buffer_t *buffer) {
	tfx_impl_webgpu__release_buffer(buffer);
	renderer->resource_generation++;
	size_in_bytes = (size_in_bytes + 15) & ~15u;
	WGPUBufferDescriptor descriptor = WGPU_BUFFER_DESCRIPTOR_INIT;
	descriptor.usage = usage;
	descriptor.size = size_in_bytes;
	buffer->buffer = wgpuDeviceCreateBuffer(renderer->device, &descriptor);
	if (!buffer->buffer) {
		tfx_Print("Could not create a buffer of %u bytes", size_in_bytes);
		return false;
	}
	buffer->capacity = size_in_bytes;
	return true;
}

//Grows by half again whenever it's too small, so the size settles after the first few frames. Never empty, so a shader always has a buffer bound
static bool tfx_impl_webgpu__ensure_buffer(tfx_impl_webgpu_t *renderer, tfxU32 size_in_bytes, WGPUBufferUsage usage, tfx_impl_webgpu_buffer_t *buffer) {
	size_in_bytes = size_in_bytes ? size_in_bytes : 16;
	if (buffer->buffer && buffer->capacity >= size_in_bytes) {
		return true;
	}
	return tfx_impl_webgpu__create_buffer(renderer, size_in_bytes + size_in_bytes / 2, usage, buffer);
}

//The queue takes a copy before returning, and the write lands after every command buffer already submitted and before the next
static void tfx_impl_webgpu__write_buffer(tfx_impl_webgpu_t *renderer, WGPUBuffer buffer, const void *data, size_t size_in_bytes) {
	if (size_in_bytes > 0) {
		wgpuQueueWriteBuffer(renderer->queue, buffer, 0, data, size_in_bytes);
	}
}

static bool tfx_impl_webgpu__upload_buffer(tfx_impl_webgpu_t *renderer, const void *data, tfxU32 size_in_bytes, WGPUBufferUsage usage, tfx_impl_webgpu_buffer_t *buffer) {
	if (!tfx_impl_webgpu__ensure_buffer(renderer, size_in_bytes, usage | WGPUBufferUsage_CopyDst, buffer)) {
		return false;
	}
	tfx_impl_webgpu__write_buffer(renderer, buffer->buffer, data, size_in_bytes);
	return true;
}

static void tfx_impl_webgpu__release_texture(WGPUTexture *texture, WGPUTextureView *view) {
	if (*view) {
		wgpuTextureViewRelease(*view);
	}
	if (*texture) {
		wgpuTextureRelease(*texture);
	}
	*view = NULL;
	*texture = NULL;
}

//Viewed as an array even with one layer, where a default view would be a plain 2D texture that the shaders' texture_2d_array can't take
static bool tfx_impl_webgpu__create_texture_array(tfx_impl_webgpu_t *renderer, tfxU32 width, tfxU32 height, tfxU32 layer_count, tfxU32 mip_count, WGPUTextureFormat format,
	WGPUTexture *texture, WGPUTextureView *view) {
	tfx_impl_webgpu__release_texture(texture, view);
	renderer->resource_generation++;
	WGPUTextureDescriptor descriptor = WGPU_TEXTURE_DESCRIPTOR_INIT;
	descriptor.usage = WGPUTextureUsage_TextureBinding | WGPUTextureUsage_CopyDst;
	descriptor.dimension = WGPUTextureDimension_2D;
	descriptor.size.width = width;
	descriptor.size.height = height;
	descriptor.size.depthOrArrayLayers = layer_count;
	descriptor.format = format;
	descriptor.mipLevelCount = mip_count;
	*texture = wgpuDeviceCreateTexture(renderer->device, &descriptor);
	if (!*texture) {
		tfx_Print("Could not create a %u x %u texture array of %u layers", width, height, layer_count);
		return false;
	}
	WGPUTextureViewDescriptor view_descriptor = WGPU_TEXTURE_VIEW_DESCRIPTOR_INIT;
	view_descriptor.format = format;
	view_descriptor.dimension = WGPUTextureViewDimension_2DArray;
	view_descriptor.mipLevelCount = mip_count;
	view_descriptor.arrayLayerCount = layer_count;
	view_descriptor.aspect = WGPUTextureAspect_All;
	*view = wgpuTextureCreateView(*texture, &view_descriptor);
	return *view != NULL;
}

static void tfx_impl_webgpu__write_texture(tfx_impl_webgpu_t *renderer, WGPUTexture texture, tfxU32 layer, tfxU32 mip_level, tfxU32 width, tfxU32 height, tfxU32 bytes_per_pixel, const void *data) {
	WGPUTexelCopyTextureInfo destination = WGPU_TEXEL_COPY_TEXTURE_INFO_INIT;
	destination.texture = texture;
	destination.mipLevel = mip_level;
	destination.origin.z = layer;
	WGPUTexelCopyBufferLayout layout = WGPU_TEXEL_COPY_BUFFER_LAYOUT_INIT;
	layout.bytesPerRow = width * bytes_per_pixel;
	layout.rowsPerImage = height;
	WGPUExtent3D size = { width, height, 1 };
	wgpuQueueWriteTexture(renderer->queue, &destination, data, (size_t)width * height * bytes_per_pixel, &layout, &size);
}

static WGPUTextureFormat tfx_impl_webgpu__color_ramp_format(tfxU32 bytes_per_pixel) {
	if (bytes_per_pixel == 4) {
		return WGPUTextureFormat_RGBA8Unorm;
	} else if (bytes_per_pixel == 16) {
		return WGPUTextureFormat_RGBA32Float;
	}
	return WGPUTextureFormat_RGBA16Float;
}

//A 2 x 2 box filter, the same as the mips D3D11 generates
static void tfx_impl_webgpu__downsample(const unsigned char *source, tfxU32 source_size, unsigned char *destination) {
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

//Grown with headroom, and replacing the buffer remakes every bind group, so Render's bind groups never need it to grow
static bool tfx_impl_webgpu__ensure_uniform_blocks(tfx_impl_webgpu_t *renderer, tfxU32 block_count) {
	if (renderer->uniforms.buffer && block_count <= renderer->uniform_block_capacity) {
		return true;
	}
	tfxU32 capacity = block_count + block_count / 2;
	unsigned char *staging = (unsigned char *)realloc(renderer->uniform_staging, (size_t)capacity * TFX_IMPL_WEBGPU_UNIFORM_BLOCK_SIZE);
	if (!staging) {
		return false;
	}
	renderer->uniform_staging = staging;
	renderer->uniform_block_capacity = 0;
	if (!tfx_impl_webgpu__create_buffer(renderer, capacity * TFX_IMPL_WEBGPU_UNIFORM_BLOCK_SIZE, WGPUBufferUsage_Uniform | WGPUBufferUsage_CopyDst, &renderer->uniforms)) {
		return false;
	}
	renderer->uniform_block_capacity = capacity;
	return true;
}

static void *tfx_impl_webgpu__uniform_block(tfx_impl_webgpu_t *renderer, tfxU32 block) {
	return renderer->uniform_staging + (size_t)block * TFX_IMPL_WEBGPU_UNIFORM_BLOCK_SIZE;
}

static void tfx_impl_webgpu__write_uniform_blocks(tfx_impl_webgpu_t *renderer, tfxU32 first_block, tfxU32 end_block) {
	if (end_block > first_block) {
		size_t offset = (size_t)first_block * TFX_IMPL_WEBGPU_UNIFORM_BLOCK_SIZE;
		wgpuQueueWriteBuffer(renderer->queue, renderer->uniforms.buffer, offset, renderer->uniform_staging + offset, (size_t)(end_block - first_block) * TFX_IMPL_WEBGPU_UNIFORM_BLOCK_SIZE);
	}
}

static WGPUBindGroupLayoutEntry tfx_impl_webgpu__buffer_layout_entry(tfxU32 binding, WGPUShaderStage visibility, WGPUBufferBindingType type, bool dynamic_offset, uint64_t minimum_size) {
	WGPUBindGroupLayoutEntry entry = WGPU_BIND_GROUP_LAYOUT_ENTRY_INIT;
	entry.binding = binding;
	entry.visibility = visibility;
	entry.buffer.type = type;
	entry.buffer.hasDynamicOffset = dynamic_offset;
	entry.buffer.minBindingSize = minimum_size;
	return entry;
}

static tfxU32 tfx_impl_webgpu__add_storage_layout_entries(WGPUBindGroupLayoutEntry *entries, tfxU32 count, tfxU32 first_binding, tfxU32 binding_count, WGPUShaderStage visibility, WGPUBufferBindingType type) {
	for (tfxU32 binding = first_binding; binding != first_binding + binding_count; ++binding) {
		entries[count++] = tfx_impl_webgpu__buffer_layout_entry(binding, visibility, type, false, 0);
	}
	return count;
}

//Every render pipeline ends with the view at 5, the draw's own block at 6, then the atlas, the color ramps and the atlas's sampler at 7, 8 and 9,
//and the view again at 10, where the lit fragments read the light. The unlit pipelines share the layout and leave 10 unread
static tfxU32 tfx_impl_webgpu__add_render_layout_entries(WGPUBindGroupLayoutEntry *entries, tfxU32 count, uint64_t draw_block_size) {
	entries[count++] = tfx_impl_webgpu__buffer_layout_entry(5, WGPUShaderStage_Vertex, WGPUBufferBindingType_Uniform, false, sizeof(tfx_impl_uniforms_t));
	entries[count++] = tfx_impl_webgpu__buffer_layout_entry(10, WGPUShaderStage_Fragment, WGPUBufferBindingType_Uniform, false, sizeof(tfx_impl_uniforms_t));
	entries[count++] = tfx_impl_webgpu__buffer_layout_entry(6, WGPUShaderStage_Vertex, WGPUBufferBindingType_Uniform, true, draw_block_size);
	WGPUBindGroupLayoutEntry atlas = WGPU_BIND_GROUP_LAYOUT_ENTRY_INIT;
	atlas.binding = 7;
	atlas.visibility = WGPUShaderStage_Fragment;
	atlas.texture.sampleType = WGPUTextureSampleType_Float;
	atlas.texture.viewDimension = WGPUTextureViewDimension_2DArray;
	entries[count++] = atlas;
	//Only loaded, never sampled, so a 32 bit float ramp works without the float32-filterable feature
	WGPUBindGroupLayoutEntry color_ramps = atlas;
	color_ramps.binding = 8;
	color_ramps.texture.sampleType = WGPUTextureSampleType_UnfilterableFloat;
	entries[count++] = color_ramps;
	WGPUBindGroupLayoutEntry sampler = WGPU_BIND_GROUP_LAYOUT_ENTRY_INIT;
	sampler.binding = 9;
	sampler.visibility = WGPUShaderStage_Fragment;
	sampler.sampler.type = WGPUSamplerBindingType_Filtering;
	entries[count++] = sampler;
	return count;
}

static WGPUBindGroupLayout tfx_impl_webgpu__create_layout(WGPUDevice device, const WGPUBindGroupLayoutEntry *entries, tfxU32 entry_count) {
	WGPUBindGroupLayoutDescriptor descriptor = WGPU_BIND_GROUP_LAYOUT_DESCRIPTOR_INIT;
	descriptor.entryCount = entry_count;
	descriptor.entries = entries;
	return wgpuDeviceCreateBindGroupLayout(device, &descriptor);
}

static WGPUBindGroupEntry tfx_impl_webgpu__buffer_entry(tfxU32 binding, WGPUBuffer buffer, uint64_t size) {
	WGPUBindGroupEntry entry = WGPU_BIND_GROUP_ENTRY_INIT;
	entry.binding = binding;
	entry.buffer = buffer;
	entry.size = size;
	return entry;
}

static tfxU32 tfx_impl_webgpu__add_render_entries(const tfx_impl_webgpu_t *renderer, WGPUBindGroupEntry *entries, tfxU32 count, uint64_t draw_block_size, WGPUTextureView color_ramps_view) {
	entries[count++] = tfx_impl_webgpu__buffer_entry(5, renderer->uniforms.buffer, sizeof(tfx_impl_uniforms_t));
	entries[count++] = tfx_impl_webgpu__buffer_entry(10, renderer->uniforms.buffer, sizeof(tfx_impl_uniforms_t));
	entries[count++] = tfx_impl_webgpu__buffer_entry(6, renderer->uniforms.buffer, draw_block_size);
	WGPUBindGroupEntry atlas = WGPU_BIND_GROUP_ENTRY_INIT;
	atlas.binding = 7;
	atlas.textureView = renderer->shape_atlas_view;
	entries[count++] = atlas;
	WGPUBindGroupEntry color_ramps = WGPU_BIND_GROUP_ENTRY_INIT;
	color_ramps.binding = 8;
	color_ramps.textureView = color_ramps_view;
	entries[count++] = color_ramps;
	WGPUBindGroupEntry sampler = WGPU_BIND_GROUP_ENTRY_INIT;
	sampler.binding = 9;
	sampler.sampler = renderer->atlas_sampler;
	entries[count++] = sampler;
	return count;
}

//Made again only when a buffer or texture has been replaced since, otherwise the one made last time is still right
static WGPUBindGroup tfx_impl_webgpu__bind_group(tfx_impl_webgpu_t *renderer, tfx_impl_webgpu_bind_group_t *cached, WGPUBindGroupLayout layout, const WGPUBindGroupEntry *entries, tfxU32 entry_count) {
	if (cached->bind_group && cached->resource_generation == renderer->resource_generation) {
		return cached->bind_group;
	}
	if (cached->bind_group) {
		wgpuBindGroupRelease(cached->bind_group);
	}
	WGPUBindGroupDescriptor descriptor = WGPU_BIND_GROUP_DESCRIPTOR_INIT;
	descriptor.layout = layout;
	descriptor.entryCount = entry_count;
	descriptor.entries = entries;
	cached->bind_group = wgpuDeviceCreateBindGroup(renderer->device, &descriptor);
	cached->resource_generation = renderer->resource_generation;
	return cached->bind_group;
}

static void tfx_impl_webgpu__release_bind_group(tfx_impl_webgpu_bind_group_t *cached) {
	if (cached->bind_group) {
		wgpuBindGroupRelease(cached->bind_group);
	}
	cached->bind_group = NULL;
}

//slangc compiles every entry point ahead of time to <entry point>.wgsl in the shader folder
static WGPUShaderModule tfx_impl_webgpu__load_shader(WGPUDevice device, const char *shader_folder, const char *entry_point) {
	char path[1024];
	snprintf(path, sizeof(path), "%s/%s.wgsl", shader_folder, entry_point);
	FILE *file = fopen(path, "rb");
	if (!file) {
		tfx_Print("Could not open the shader %s", path);
		return NULL;
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
		return NULL;
	}
	WGPUShaderSourceWGSL source = WGPU_SHADER_SOURCE_WGSL_INIT;
	source.code.data = code;
	source.code.length = (size_t)size;
	WGPUShaderModuleDescriptor descriptor = WGPU_SHADER_MODULE_DESCRIPTOR_INIT;
	descriptor.nextInChain = &source.chain;
	descriptor.label = tfx_impl_webgpu__string(entry_point);
	WGPUShaderModule module = wgpuDeviceCreateShaderModule(device, &descriptor);
	free(code);
	return module;
}

static WGPUPipelineLayout tfx_impl_webgpu__create_pipeline_layout(WGPUDevice device, WGPUBindGroupLayout layout) {
	WGPUPipelineLayoutDescriptor descriptor = WGPU_PIPELINE_LAYOUT_DESCRIPTOR_INIT;
	descriptor.bindGroupLayoutCount = 1;
	descriptor.bindGroupLayouts = &layout;
	return wgpuDeviceCreatePipelineLayout(device, &descriptor);
}

static WGPUComputePipeline tfx_impl_webgpu__create_compute_pipeline(WGPUDevice device, const char *shader_folder, const char *entry_point, WGPUBindGroupLayout layout) {
	WGPUShaderModule module = tfx_impl_webgpu__load_shader(device, shader_folder, entry_point);
	if (!module || !layout) {
		if (module) wgpuShaderModuleRelease(module);
		return NULL;
	}
	WGPUPipelineLayout pipeline_layout = tfx_impl_webgpu__create_pipeline_layout(device, layout);
	WGPUComputePipelineDescriptor descriptor = WGPU_COMPUTE_PIPELINE_DESCRIPTOR_INIT;
	descriptor.label = tfx_impl_webgpu__string(entry_point);
	descriptor.layout = pipeline_layout;
	descriptor.compute.module = module;
	descriptor.compute.entryPoint = tfx_impl_webgpu__string(entry_point);
	WGPUComputePipeline pipeline = wgpuDeviceCreateComputePipeline(device, &descriptor);
	wgpuPipelineLayoutRelease(pipeline_layout);
	wgpuShaderModuleRelease(module);
	return pipeline;
}

//Premultiplied blending, depth tested but not written, and no vertex buffers: instances and ribbon vertexes are read from storage buffers
//Premultiplied blending, depth tested but not written, or for a layer that writes depth, blending off and depth written
static WGPURenderPipeline tfx_impl_webgpu__create_render_pipeline(WGPUDevice device, WGPUShaderModule vertex_module, const char *vertex_entry_point,
	WGPUShaderModule fragment_module, const char *fragment_entry_point, WGPUBindGroupLayout layout, WGPUTextureFormat color_format, WGPUTextureFormat depth_format, bool writes_depth) {
	if (!vertex_module || !fragment_module || !layout) {
		return NULL;
	}
	WGPUBlendState blend = WGPU_BLEND_STATE_INIT;
	blend.color.operation = WGPUBlendOperation_Add;
	blend.color.srcFactor = WGPUBlendFactor_One;
	blend.color.dstFactor = WGPUBlendFactor_OneMinusSrcAlpha;
	blend.alpha = blend.color;
	WGPUColorTargetState color_target = WGPU_COLOR_TARGET_STATE_INIT;
	color_target.format = color_format;
	color_target.blend = writes_depth ? NULL : &blend;
	WGPUFragmentState fragment = WGPU_FRAGMENT_STATE_INIT;
	fragment.module = fragment_module;
	fragment.entryPoint = tfx_impl_webgpu__string(fragment_entry_point);
	fragment.targetCount = 1;
	fragment.targets = &color_target;
	WGPUDepthStencilState depth_stencil = WGPU_DEPTH_STENCIL_STATE_INIT;
	depth_stencil.format = depth_format;
	depth_stencil.depthWriteEnabled = writes_depth ? WGPUOptionalBool_True : WGPUOptionalBool_False;
	depth_stencil.depthCompare = WGPUCompareFunction_LessEqual;
	WGPUPipelineLayout pipeline_layout = tfx_impl_webgpu__create_pipeline_layout(device, layout);
	WGPURenderPipelineDescriptor descriptor = WGPU_RENDER_PIPELINE_DESCRIPTOR_INIT;
	descriptor.label = tfx_impl_webgpu__string(vertex_entry_point);
	descriptor.layout = pipeline_layout;
	descriptor.vertex.module = vertex_module;
	descriptor.vertex.entryPoint = tfx_impl_webgpu__string(vertex_entry_point);
	descriptor.primitive.topology = WGPUPrimitiveTopology_TriangleList;
	descriptor.primitive.cullMode = WGPUCullMode_None;
	descriptor.depthStencil = &depth_stencil;
	descriptor.fragment = &fragment;
	WGPURenderPipeline pipeline = wgpuDeviceCreateRenderPipeline(device, &descriptor);
	wgpuPipelineLayoutRelease(pipeline_layout);
	return pipeline;
}

static void tfx_impl_webgpu__init_scope_callback(WGPUPopErrorScopeStatus status, WGPUErrorType type, WGPUStringView message, void *userdata1, void *userdata2) {
	(void)userdata2;
	tfx_impl_webgpu_t *renderer = (tfx_impl_webgpu_t *)userdata1;
	renderer->init_scope_popped = true;
	if (status == WGPUPopErrorScopeStatus_Success && type != WGPUErrorType_NoError) {
		renderer->init_scope_failed = true;
		tfx_Print("Could not create the TimelineFX pipelines: %.*s", (int)message.length, message.data ? message.data : "");
	}
}

//--- Webgpu_functions

bool tfx_impl_webgpu_Init(tfx_impl_webgpu_t *renderer, WGPUDevice device, WGPUTextureFormat color_format, WGPUTextureFormat depth_format, const char *shader_folder) {
	memset(renderer, 0, sizeof(tfx_impl_webgpu_t));
	renderer->device = device;
	renderer->queue = wgpuDeviceGetQueue(device);
	//A shader or layout mistake shows up here rather than as an invalid pipeline on the first draw
	wgpuDevicePushErrorScope(device, WGPUErrorFilter_Validation);

	//The slots of each entry point are in the Shaders readme. Vertex shaders only ever read storage buffers, which WebGPU requires
	WGPUBindGroupLayoutEntry entries[12];
	tfxU32 count = tfx_impl_webgpu__add_storage_layout_entries(entries, 0, 0, 4, WGPUShaderStage_Vertex, WGPUBufferBindingType_ReadOnlyStorage);
	count = tfx_impl_webgpu__add_render_layout_entries(entries, count, sizeof(tfx_impl_webgpu_draw_constants_t));
	renderer->billboard_layout = tfx_impl_webgpu__create_layout(device, entries, count);
	count = tfx_impl_webgpu__add_storage_layout_entries(entries, 0, 0, 3, WGPUShaderStage_Vertex, WGPUBufferBindingType_ReadOnlyStorage);
	count = tfx_impl_webgpu__add_render_layout_entries(entries, count, sizeof(tfx_impl_webgpu_draw_constants_t));
	renderer->baked_billboard_layout = tfx_impl_webgpu__create_layout(device, entries, count);
	count = tfx_impl_webgpu__add_storage_layout_entries(entries, 0, 0, 5, WGPUShaderStage_Vertex, WGPUBufferBindingType_ReadOnlyStorage);
	count = tfx_impl_webgpu__add_render_layout_entries(entries, count, sizeof(tfx_ribbon_bucket_globals_t));
	renderer->ribbon_layout = tfx_impl_webgpu__create_layout(device, entries, count);
	//Compute writes from binding 0 and reads from binding 2, with its block at 6
	count = tfx_impl_webgpu__add_storage_layout_entries(entries, 0, 0, 2, WGPUShaderStage_Compute, WGPUBufferBindingType_Storage);
	count = tfx_impl_webgpu__add_storage_layout_entries(entries, count, 2, 2, WGPUShaderStage_Compute, WGPUBufferBindingType_ReadOnlyStorage);
	entries[count++] = tfx_impl_webgpu__buffer_layout_entry(6, WGPUShaderStage_Compute, WGPUBufferBindingType_Uniform, true, sizeof(tfx_gpu_particle_push_t));
	renderer->particle_kernel_layout = tfx_impl_webgpu__create_layout(device, entries, count);
	count = tfx_impl_webgpu__add_storage_layout_entries(entries, 0, 0, 2, WGPUShaderStage_Compute, WGPUBufferBindingType_Storage);
	count = tfx_impl_webgpu__add_storage_layout_entries(entries, count, 2, 4, WGPUShaderStage_Compute, WGPUBufferBindingType_ReadOnlyStorage);
	entries[count++] = tfx_impl_webgpu__buffer_layout_entry(6, WGPUShaderStage_Compute, WGPUBufferBindingType_Uniform, true, sizeof(tfx_ribbon_bucket_globals_t));
	renderer->ribbon_compute_layout = tfx_impl_webgpu__create_layout(device, entries, count);
	count = tfx_impl_webgpu__add_storage_layout_entries(entries, 0, 0, 1, WGPUShaderStage_Compute, WGPUBufferBindingType_Storage);
	count = tfx_impl_webgpu__add_storage_layout_entries(entries, count, 2, 4, WGPUShaderStage_Compute, WGPUBufferBindingType_ReadOnlyStorage);
	//No minimum size, because the two playback blocks differ in size, so it's checked at each dispatch instead
	entries[count++] = tfx_impl_webgpu__buffer_layout_entry(6, WGPUShaderStage_Compute, WGPUBufferBindingType_Uniform, true, 0);
	renderer->playback_layout = tfx_impl_webgpu__create_layout(device, entries, count);

	WGPUShaderModule billboard_vertex_modules[2];
	WGPUShaderModule ribbon_vertex_modules[2];
	WGPUShaderModule billboard_fragment_modules[2][2];
	WGPUShaderModule ribbon_fragment_modules[2][2];
	for (int lit = 0; lit != 2; ++lit) {
		billboard_vertex_modules[lit] = tfx_impl_webgpu__load_shader(device, shader_folder, tfx_impl_billboard_vertex_entry_points[lit]);
		ribbon_vertex_modules[lit] = tfx_impl_webgpu__load_shader(device, shader_folder, tfx_impl_ribbon_vertex_entry_points[lit]);
		for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
			billboard_fragment_modules[writes_depth][lit] = tfx_impl_webgpu__load_shader(device, shader_folder, tfx_impl_billboard_fragment_entry_points[writes_depth][lit]);
			ribbon_fragment_modules[writes_depth][lit] = tfx_impl_webgpu__load_shader(device, shader_folder, tfx_impl_ribbon_fragment_entry_points[writes_depth][lit]);
		}
	}
	WGPUShaderModule baked_billboard_vertex = tfx_impl_webgpu__load_shader(device, shader_folder, "tfx_billboard_baked_vertex");
	for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
		for (int lit = 0; lit != 2; ++lit) {
			renderer->billboard_pipelines[writes_depth][lit] = tfx_impl_webgpu__create_render_pipeline(device, billboard_vertex_modules[lit], tfx_impl_billboard_vertex_entry_points[lit],
				billboard_fragment_modules[writes_depth][lit], tfx_impl_billboard_fragment_entry_points[writes_depth][lit], renderer->billboard_layout, color_format, depth_format, writes_depth != 0);
			renderer->ribbon_pipelines[writes_depth][lit] = tfx_impl_webgpu__create_render_pipeline(device, ribbon_vertex_modules[lit], tfx_impl_ribbon_vertex_entry_points[lit],
				ribbon_fragment_modules[writes_depth][lit], tfx_impl_ribbon_fragment_entry_points[writes_depth][lit], renderer->ribbon_layout, color_format, depth_format, writes_depth != 0);
		}
	}
	renderer->baked_billboard_pipeline = tfx_impl_webgpu__create_render_pipeline(device, baked_billboard_vertex, "tfx_billboard_baked_vertex", billboard_fragment_modules[0][0], "tfx_billboard_fragment",
		renderer->baked_billboard_layout, color_format, depth_format, false);
	//The pipelines keep what they need from the modules
	WGPUShaderModule modules[13] = { baked_billboard_vertex, billboard_vertex_modules[0], billboard_vertex_modules[1], ribbon_vertex_modules[0], ribbon_vertex_modules[1],
		billboard_fragment_modules[0][0], billboard_fragment_modules[0][1], billboard_fragment_modules[1][0], billboard_fragment_modules[1][1],
		ribbon_fragment_modules[0][0], ribbon_fragment_modules[0][1], ribbon_fragment_modules[1][0], ribbon_fragment_modules[1][1] };
	for (int index = 0; index != 13; ++index) {
		if (modules[index]) {
			wgpuShaderModuleRelease(modules[index]);
		}
	}
	bool pipelines_created = true;
	for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
		for (int lit = 0; lit != 2; ++lit) {
			pipelines_created = pipelines_created && renderer->billboard_pipelines[writes_depth][lit] && renderer->ribbon_pipelines[writes_depth][lit];
		}
	}
	renderer->particle_kernel = tfx_impl_webgpu__create_compute_pipeline(device, shader_folder, "simulate_particles", renderer->particle_kernel_layout);
	renderer->ribbon_compute_pipeline = tfx_impl_webgpu__create_compute_pipeline(device, shader_folder, "tfx_ribbon_compute", renderer->ribbon_compute_layout);
	renderer->sprite_playback_pipeline = tfx_impl_webgpu__create_compute_pipeline(device, shader_folder, "tfx_sprite_playback", renderer->playback_layout);
	renderer->ribbon_playback_pipeline = tfx_impl_webgpu__create_compute_pipeline(device, shader_folder, "tfx_ribbon_playback", renderer->playback_layout);

	WGPUSamplerDescriptor sampler_descriptor = WGPU_SAMPLER_DESCRIPTOR_INIT;
	sampler_descriptor.magFilter = WGPUFilterMode_Linear;
	sampler_descriptor.minFilter = WGPUFilterMode_Linear;
	sampler_descriptor.mipmapFilter = WGPUMipmapFilterMode_Linear;
	renderer->atlas_sampler = wgpuDeviceCreateSampler(device, &sampler_descriptor);

	//wgpu-native calls this before returning. A WebGPU that calls it later still prints the error, but Init can't fail on it
	WGPUPopErrorScopeCallbackInfo scope_callback = WGPU_POP_ERROR_SCOPE_CALLBACK_INFO_INIT;
	scope_callback.mode = WGPUCallbackMode_AllowSpontaneous;
	scope_callback.callback = tfx_impl_webgpu__init_scope_callback;
	scope_callback.userdata1 = renderer;
	wgpuDevicePopErrorScope(device, scope_callback);

	if (renderer->init_scope_failed || !pipelines_created || !renderer->baked_billboard_pipeline || !renderer->particle_kernel
		|| !renderer->ribbon_compute_pipeline || !renderer->sprite_playback_pipeline || !renderer->ribbon_playback_pipeline || !renderer->atlas_sampler
		|| !tfx_impl_webgpu__ensure_uniform_blocks(renderer, TFX_IMPL_WEBGPU_INITIAL_UNIFORM_BLOCKS)) {
		tfx_Print("Could not create the TimelineFX pipelines and samplers");
		return false;
	}
	return true;
}

bool tfx_impl_webgpu_CreateAtlas(tfx_impl_webgpu_t *renderer, tfxU32 layer_size, tfxU32 layer_count) {
	tfxU32 mip_count = 1;
	while ((layer_size >> mip_count) > 0) {
		mip_count++;
	}
	if (!tfx_impl_webgpu__create_texture_array(renderer, layer_size, layer_size, layer_count, mip_count, WGPUTextureFormat_RGBA8Unorm, &renderer->shape_atlas, &renderer->shape_atlas_view)) {
		return false;
	}
	renderer->shape_atlas_layer_size = layer_size;
	renderer->shape_atlas_layer_count = layer_count;
	renderer->shape_atlas_mip_count = mip_count;
	return true;
}

bool tfx_impl_webgpu_UploadAtlas(tfx_impl_webgpu_t *renderer, const tfx_impl_atlas_t *atlas) {
	if (!renderer->shape_atlas || atlas->layer_size != renderer->shape_atlas_layer_size
		|| atlas->first_layer + atlas->layer_count > renderer->shape_atlas_layer_count) {
		tfx_Print("An atlas of %u layers from layer %u, %u pixels square, doesn't fit the %u layer, %u pixel texture array",
			atlas->layer_count, atlas->first_layer, atlas->layer_size, renderer->shape_atlas_layer_count, renderer->shape_atlas_layer_size);
		return false;
	}
	//The queue copies each level as it's written, so one level's worth of scratch does for every layer
	size_t layer_bytes = (size_t)atlas->layer_size * atlas->layer_size * 4;
	unsigned char *mips = (unsigned char *)malloc(layer_bytes / 2);
	if (!mips) {
		return false;
	}
	for (tfxU32 layer = 0; layer != atlas->layer_count; ++layer) {
		const unsigned char *source = atlas->pixels + layer * layer_bytes;
		tfxU32 size = atlas->layer_size;
		tfx_impl_webgpu__write_texture(renderer, renderer->shape_atlas, atlas->first_layer + layer, 0, size, size, 4, source);
		unsigned char *destination = mips;
		for (tfxU32 mip_level = 1; mip_level < renderer->shape_atlas_mip_count; ++mip_level) {
			tfx_impl_webgpu__downsample(source, size, destination);
			size = size > 1 ? size / 2 : 1;
			tfx_impl_webgpu__write_texture(renderer, renderer->shape_atlas, atlas->first_layer + layer, mip_level, size, size, 4, destination);
			source = destination;
			destination += (size_t)size * size * 4;
		}
	}
	free(mips);
	return true;
}

bool tfx_impl_webgpu_UploadSpriteData(tfx_impl_webgpu_t *renderer, tfx_animation_manager animation_manager, tfx_gpu_shapes shapes) {
	tfxU32 bitmap_count = tfx_GetAnimationColorRampBitmapCount(animation_manager);
	tfx_impl_webgpu__release_texture(&renderer->baked_color_ramps, &renderer->baked_color_ramps_view);
	if (bitmap_count > 0) {
		tfx_bitmap_t *first_bitmap = tfx_GetAnimationColorRampBitmap(animation_manager, 0);
		tfxU32 width = (tfxU32)tfx_GetBitmapWidth(first_bitmap);
		tfxU32 height = (tfxU32)tfx_GetBitmapHeight(first_bitmap);
		tfxU32 bytes_per_pixel = (tfxU32)(tfx_GetBitmapSize(first_bitmap) / ((size_t)width * height));
		if (!tfx_impl_webgpu__create_texture_array(renderer, width, height, bitmap_count, 1, tfx_impl_webgpu__color_ramp_format(bytes_per_pixel),
			&renderer->baked_color_ramps, &renderer->baked_color_ramps_view)) {
			return false;
		}
		for (tfxU32 index = 0; index != bitmap_count; ++index) {
			tfx_bitmap_t *bitmap = tfx_GetAnimationColorRampBitmap(animation_manager, index);
			tfx_impl_webgpu__write_texture(renderer, renderer->baked_color_ramps, index, 0, width, height, bytes_per_pixel, tfx_GetBitmapData(bitmap));
		}
	}

	tfxU32 shape_bytes = (tfxU32)sizeof(tfx_gpu_image_data_t) * tfx_GetGPUShapesCount(shapes);
	tfxU32 property_bytes = (tfxU32)tfx_GetAnimationEmitterPropertySizeInBytes(animation_manager);
	tfxU32 sprite_bytes = (tfxU32)tfx_GetSpriteDataSizeInBytes(animation_manager);
	if (!tfx_impl_webgpu__upload_buffer(renderer, tfx_GetGPUShapesArray(shapes), shape_bytes, WGPUBufferUsage_Storage, &renderer->baked_image_data)
		|| !tfx_impl_webgpu__upload_buffer(renderer, tfx_GetAnimationEmitterPropertiesBufferPointer(animation_manager), property_bytes, WGPUBufferUsage_Storage, &renderer->baked_emitter_properties)
		|| !tfx_impl_webgpu__upload_buffer(renderer, tfx_GetSpriteDataBufferPointer(animation_manager), sprite_bytes, WGPUBufferUsage_Storage, &renderer->baked_sprites)) {
		return false;
	}
	if (!tfx_AnimationManagerHasRibbons(animation_manager)) {
		return true;
	}

	//The baked ribbons, their properties, segments, emitters and graph lookups only change when sprite data is loaded
	tfxU32 emitters_bytes = (tfxU32)tfx_GetAnimationRibbonEmittersSizeInBytes(animation_manager);
	void *emitters = malloc(emitters_bytes ? emitters_bytes : 1);
	if (!emitters) {
		return false;
	}
	tfx_CopyAnimationRibbonEmitters(animation_manager, emitters);
	bool uploaded = tfx_impl_webgpu__upload_buffer(renderer, tfx_GetAnimationRibbonGraphLookupsBufferPointer(animation_manager), (tfxU32)tfx_GetAnimationRibbonGraphLookupsSizeInBytes(animation_manager),
			WGPUBufferUsage_Storage, &renderer->baked_ribbons.graph_lookups)
		&& tfx_impl_webgpu__upload_buffer(renderer, tfx_GetAnimationRibbonDataBufferPointer(animation_manager), (tfxU32)tfx_GetAnimationRibbonDataSizeInBytes(animation_manager),
			WGPUBufferUsage_Storage, &renderer->baked_ribbon_data)
		&& tfx_impl_webgpu__upload_buffer(renderer, tfx_GetAnimationRibbonPropertiesBufferPointer(animation_manager), (tfxU32)tfx_GetAnimationRibbonPropertySizeInBytes(animation_manager),
			WGPUBufferUsage_Storage, &renderer->baked_ribbon_properties)
		&& tfx_impl_webgpu__upload_buffer(renderer, tfx_GetAnimationRibbonSegmentsBufferPointer(animation_manager), (tfxU32)tfx_GetAnimationRibbonSegmentsSizeInBytes(animation_manager),
			WGPUBufferUsage_Storage, &renderer->baked_ribbons.segments)
		&& tfx_impl_webgpu__upload_buffer(renderer, emitters, emitters_bytes, WGPUBufferUsage_Storage, &renderer->baked_ribbons.emitters);
	free(emitters);
	return uploaded;
}

//The vertex and index buffers ribbon_compute writes, sized for this frame's buckets
static bool tfx_impl_webgpu__ensure_ribbon_geometry(tfx_impl_webgpu_t *renderer, const tfx_ribbon_batch_sizes_t *sizes, tfx_impl_webgpu_ribbon_set_t *set) {
	return tfx_impl_webgpu__ensure_buffer(renderer, (tfxU32)sizes->vertex_buffer_size_in_bytes, WGPUBufferUsage_Storage, &set->vertices)
		&& tfx_impl_webgpu__ensure_buffer(renderer, (tfxU32)sizes->index_buffer_size_in_bytes, WGPUBufferUsage_Storage | WGPUBufferUsage_Index, &set->indices);
}

static bool tfx_impl_webgpu__upload_baked_ribbons(tfx_impl_webgpu_t *renderer, const tfx_impl_baked_frame_t *baked) {
	if (baked->ribbon_count == 0 || !renderer->baked_ribbon_data.buffer) {
		return true;
	}
	tfx_ribbon_batch_sizes_t sizes = tfx_GetAnimationRibbonBufferSizes(baked->animation_manager, sizeof(tfx_ribbon_vertex_t));
	//ribbon_playback writes the ribbons, then ribbon_compute and the ribbon draw read them
	if (!tfx_impl_webgpu__upload_buffer(renderer, baked->ribbon_offsets, (tfxU32)baked->ribbon_offsets_bytes, WGPUBufferUsage_Storage, &renderer->animation_ribbon_offsets)
		|| !tfx_impl_webgpu__ensure_buffer(renderer, baked->ribbon_count * tfxRIBBON_SIZE, WGPUBufferUsage_Storage, &renderer->baked_ribbons.instances)
		|| !tfx_impl_webgpu__ensure_ribbon_geometry(renderer, &sizes, &renderer->baked_ribbons)) {
		return false;
	}
	memset(&renderer->ribbon_data_push, 0, sizeof(tfx_ribbon_data_push_t));
	renderer->ribbon_data_push.animation_instances_total = baked->animation_instance_count;
	renderer->ribbon_data_push.ribbons_total = baked->ribbon_count;
	renderer->playback_ribbon_count = baked->ribbon_count;
	renderer->baked_ribbons.animation_manager = baked->animation_manager;
	return true;
}

//Animation instances advance every frame, so this runs every frame, not just on updates
static bool tfx_impl_webgpu__upload_baked(tfx_impl_webgpu_t *renderer, const tfx_impl_baked_frame_t *baked) {
	renderer->playback_sprite_count = 0;
	renderer->playback_ribbon_count = 0;
	renderer->baked_ribbons.animation_manager = NULL;
	if (!renderer->baked_sprites.buffer || (baked->sprite_count == 0 && baked->ribbon_count == 0)) {
		return true;
	}
	if (!tfx_impl_webgpu__upload_buffer(renderer, baked->animation_instances, (tfxU32)baked->animation_instances_bytes, WGPUBufferUsage_Storage, &renderer->animation_instances)) {
		return false;
	}
	if (baked->sprite_count > 0) {
		if (!tfx_impl_webgpu__upload_buffer(renderer, baked->offsets, (tfxU32)baked->offsets_bytes, WGPUBufferUsage_Storage, &renderer->animation_offsets)
			|| !tfx_impl_webgpu__ensure_buffer(renderer, baked->sprite_count * (tfxU32)sizeof(tfx_instance_t), WGPUBufferUsage_Storage, &renderer->playback_instances)) {
			return false;
		}
		memset(&renderer->sprite_data_push, 0, sizeof(tfx_sprite_data_push_t));
		renderer->sprite_data_push.animation_instances_total = baked->animation_instance_count;
		renderer->sprite_data_push.billboards_total = baked->sprite_count;
		renderer->sprite_data_push.animated_shapes = baked->animated_shapes ? 1 : 0;
		renderer->playback_sprite_count = baked->sprite_count;
	}
	return tfx_impl_webgpu__upload_baked_ribbons(renderer, baked);
}

static bool tfx_impl_webgpu__upload_ribbons(tfx_impl_webgpu_t *renderer, const tfx_impl_ribbon_frame_t *ribbons) {
	tfx_impl_webgpu_ribbon_set_t *set = &renderer->live_ribbons;
	set->stage = ribbons->stage;
	if (!ribbons->stage) {
		return true;
	}
	const tfx_ribbon_batch_sizes_t *sizes = &ribbons->sizes;
	if (!tfx_impl_webgpu__upload_buffer(renderer, ribbons->segments, (tfxU32)sizes->segment_buffer_size_in_bytes, WGPUBufferUsage_Storage, &set->segments)
		|| !tfx_impl_webgpu__upload_buffer(renderer, ribbons->ribbons, (tfxU32)sizes->ribbon_buffer_size_in_bytes, WGPUBufferUsage_Storage, &set->instances)
		|| !tfx_impl_webgpu__upload_buffer(renderer, ribbons->emitters, (tfxU32)sizes->emitter_buffer_size_in_bytes, WGPUBufferUsage_Storage, &set->emitters)
		|| !tfx_impl_webgpu__ensure_ribbon_geometry(renderer, sizes, set)) {
		set->stage = NULL;
		return false;
	}
	set->stage_offsets = ribbons->stage_offsets;
	return true;
}

//The old contents are copied across, because the previous buffer is still what captured_index points into. The copy is submitted
//on its own at once: queue writes run before the command buffers of the submit they precede, so recorded with the frame it would
//land on top of this frame's instance upload
static bool tfx_impl_webgpu__grow_instance_buffers(tfx_impl_webgpu_t *renderer, tfxU32 needed) {
	tfxU32 capacity = renderer->instance_capacity ? renderer->instance_capacity : TFX_IMPL_WEBGPU_INITIAL_INSTANCE_CAPACITY;
	while (capacity < needed) {
		capacity += capacity / 2;
	}
	//The GPU particle kernel writes its instances after the CPU ones
	WGPUBufferUsage usage = WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst | WGPUBufferUsage_CopySrc;
	WGPUCommandEncoder encoder = NULL;
	bool grown = true;
	for (int index = 0; index != 2 && grown; ++index) {
		tfx_impl_webgpu_buffer_t old = renderer->instance_buffers[index];
		renderer->instance_buffers[index].buffer = NULL;
		grown = tfx_impl_webgpu__create_buffer(renderer, capacity * (tfxU32)sizeof(tfx_instance_t), usage, &renderer->instance_buffers[index]);
		if (old.buffer) {
			if (grown) {
				encoder = encoder ? encoder : wgpuDeviceCreateCommandEncoder(renderer->device, NULL);
				wgpuCommandEncoderCopyBufferToBuffer(encoder, old.buffer, 0, renderer->instance_buffers[index].buffer, 0, old.capacity);
			}
			wgpuBufferRelease(old.buffer);
		}
	}
	if (encoder) {
		WGPUCommandBuffer commands = wgpuCommandEncoderFinish(encoder, NULL);
		wgpuQueueSubmit(renderer->queue, 1, &commands);
		wgpuCommandBufferRelease(commands);
		wgpuCommandEncoderRelease(encoder);
	}
	renderer->instance_capacity = grown ? capacity : 0;
	return grown;
}

static bool tfx_impl_webgpu__upload_gpu_particles(tfx_impl_webgpu_t *renderer, const tfx_impl_frame_t *frame) {
	const tfx_gpu_particle_frame_t *gpu_frame = &frame->gpu_particles;
	//The pool holds the particles between frames. It's sized once and starts zeroed, as every new WebGPU buffer does
	if (!renderer->particle_pool.buffer && !tfx_impl_webgpu__create_buffer(renderer, (tfxU32)frame->gpu_pool_sizes.pool_bytes, WGPUBufferUsage_Storage, &renderer->particle_pool)) {
		return false;
	}
	//tfx_CopyGPUParticleUpload fills memory in place, so it goes through a CPU copy that the queue copies again
	if (gpu_frame->upload_bytes > renderer->particle_upload_staging_capacity) {
		size_t capacity = gpu_frame->upload_bytes + gpu_frame->upload_bytes / 2;
		void *staging = realloc(renderer->particle_upload_staging, capacity);
		if (!staging) {
			return false;
		}
		renderer->particle_upload_staging = staging;
		renderer->particle_upload_staging_capacity = capacity;
	}
	tfx_CopyGPUParticleUpload(gpu_frame, renderer->particle_upload_staging);
	if (!tfx_impl_webgpu__upload_buffer(renderer, renderer->particle_upload_staging, (tfxU32)gpu_frame->upload_bytes, WGPUBufferUsage_Storage, &renderer->particle_upload)) {
		return false;
	}

	//With fixed bindings the buffer indexes are unused, the bases say where each section starts in the pool and the upload
	tfx_gpu_particle_push_t *push = &renderer->particle_push_constants;
	memset(push, 0, sizeof(tfx_gpu_particle_push_t));
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
static bool tfx_impl_webgpu__upload_color_ramps(tfx_impl_webgpu_t *renderer) {
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
		capacity = capacity < TFX_IMPL_WEBGPU_MAX_COLOR_RAMP_LAYERS ? capacity : TFX_IMPL_WEBGPU_MAX_COLOR_RAMP_LAYERS;
		tfxU32 *layer_versions = (tfxU32 *)realloc(renderer->color_ramp_layer_versions, sizeof(tfxU32) * capacity);
		if (!layer_versions) {
			return false;
		}
		renderer->color_ramp_layer_versions = layer_versions;
		renderer->color_ramp_layer_capacity = 0;
		if (!tfx_impl_webgpu__create_texture_array(renderer, width, height, capacity, 1, tfx_impl_webgpu__color_ramp_format(bytes_per_pixel), &renderer->color_ramps, &renderer->color_ramps_view)) {
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
			tfx_impl_webgpu__write_texture(renderer, renderer->color_ramps, layer, 0, width, height, bytes_per_pixel, tfx_GetBitmapData(bitmap));
		}
	}
	return true;
}

//Only on update frames: the instances drawn until then were written with the layout the last upload holds
static bool tfx_impl_webgpu__upload_library_data(tfx_impl_webgpu_t *renderer) {
	tfxU32 properties_version = tfx_GetGlobalParticlePropertiesVersion();
	if (!renderer->emitter_properties.buffer || properties_version != renderer->properties_version) {
		if (!tfx_impl_webgpu__upload_buffer(renderer, tfx_GetGlobalParticlePropertiesBuffer(), (tfxU32)tfx_GetGlobalParticlePropertiesBufferSizeInBytes(),
			WGPUBufferUsage_Storage, &renderer->emitter_properties)) {
			return false;
		}
		renderer->properties_version = properties_version;
	}
	tfxU32 shapes_version = tfx_GetGlobalGPUShapesVersion();
	if (!renderer->image_data.buffer || shapes_version != renderer->shapes_version) {
		tfx_gpu_shapes shapes = tfx_GetGlobalGPUShapes();
		if (!tfx_impl_webgpu__upload_buffer(renderer, tfx_GetGPUShapesArray(shapes), (tfxU32)sizeof(tfx_gpu_image_data_t) * tfx_GetGPUShapesCount(shapes),
			WGPUBufferUsage_Storage, &renderer->image_data)) {
			return false;
		}
		renderer->shapes_version = shapes_version;
	}
	//Shared by every library already, read by the ribbon shaders and the GPU kernel
	tfxU32 graph_lookups_version = tfx_GetGPUGraphLookupsVersion();
	if (!renderer->graph_lookups.buffer || graph_lookups_version != renderer->graph_lookups_version) {
		if (!tfx_impl_webgpu__upload_buffer(renderer, tfx_GetGPUGraphLookupsBuffer(), tfx_GetGPUGraphLookupsBufferSizeInBytes(), WGPUBufferUsage_Storage, &renderer->graph_lookups)) {
			return false;
		}
		renderer->graph_lookups_version = graph_lookups_version;
	}
	tfxU32 color_ramps_version = tfx_GetGlobalColorRampsVersion();
	if (!renderer->color_ramps || color_ramps_version != renderer->color_ramps_version) {
		if (!tfx_impl_webgpu__upload_color_ramps(renderer)) {
			return false;
		}
		renderer->color_ramps_version = color_ramps_version;
	}
	return true;
}

bool tfx_impl_webgpu_UploadFrame(tfx_impl_webgpu_t *renderer, const tfx_impl_frame_t *frame) {
	memcpy(renderer->layer_flags, frame->layer_flags, sizeof(renderer->layer_flags));
	if (!tfx_impl_webgpu__upload_baked(renderer, &frame->baked)) {
		return false;
	}
	if (!frame->stage_updated) {
		return true;
	}
	if (!tfx_impl_webgpu__upload_library_data(renderer)) {
		return false;
	}
	bool gpu_particles = frame->has_gpu_particles && renderer->particle_kernel;
	tfxU32 needed = gpu_particles ? frame->gpu_particles.gpu_instance_start + frame->gpu_particles.gpu_instance_count : frame->instance_count;
	if ((needed > renderer->instance_capacity || renderer->instance_capacity == 0) && !tfx_impl_webgpu__grow_instance_buffers(renderer, needed)) {
		renderer->instance_count = 0;
		return false;
	}

	//The other buffer now holds what was drawn last, which is what captured_index points into for CPU and GPU instances alike
	renderer->current_instance_buffer ^= 1;
	tfx_impl_webgpu__write_buffer(renderer, renderer->instance_buffers[renderer->current_instance_buffer].buffer, frame->instances, frame->instance_count * sizeof(tfx_instance_t));
	renderer->instance_count = needed;
	memcpy(renderer->layer_start, frame->layer_start, sizeof(renderer->layer_start));
	memcpy(renderer->layer_count, frame->layer_count, sizeof(renderer->layer_count));
	memset(renderer->gpu_layer_count, 0, sizeof(renderer->gpu_layer_count));
	if (gpu_particles) {
		memcpy(renderer->gpu_layer_start, frame->gpu_particles.gpu_layer_start, sizeof(renderer->gpu_layer_start));
		memcpy(renderer->gpu_layer_count, frame->gpu_particles.gpu_layer_count, sizeof(renderer->gpu_layer_count));
		if (!tfx_impl_webgpu__upload_gpu_particles(renderer, frame)) {
			return false;
		}
	}
	return tfx_impl_webgpu__upload_ribbons(renderer, &frame->ribbons);
}

static bool tfx_impl_webgpu__ribbon_set_active(const tfx_impl_webgpu_ribbon_set_t *set) {
	return set->stage || set->animation_manager;
}

static tfx_ribbon_dispatch_t tfx_impl_webgpu__create_ribbon_dispatch(const tfx_impl_webgpu_ribbon_set_t *set) {
	return set->stage ? tfx_CreateRibbonBatchDispatch(&set->stage_offsets) : tfx_CreateAnimationRibbonDispatch();
}

static bool tfx_impl_webgpu__next_ribbon_dispatch(const tfx_impl_webgpu_ribbon_set_t *set, tfx_ribbon_dispatch_t *dispatch) {
	return set->stage ? tfx_NextRibbonDispatch(set->stage, dispatch) : tfx_NextAnimationRibbonDispatch(set->animation_manager, dispatch);
}

static tfxU32 tfx_impl_webgpu__count_ribbon_buckets(const tfx_impl_webgpu_ribbon_set_t *set) {
	if (!tfx_impl_webgpu__ribbon_set_active(set)) {
		return 0;
	}
	tfxU32 bucket_count = 0;
	tfx_ribbon_dispatch_t dispatch = tfx_impl_webgpu__create_ribbon_dispatch(set);
	while (tfx_impl_webgpu__next_ribbon_dispatch(set, &dispatch)) {
		bucket_count++;
	}
	return bucket_count;
}

static WGPUBuffer tfx_impl_webgpu__ribbon_graph_lookups(const tfx_impl_webgpu_t *renderer, const tfx_impl_webgpu_ribbon_set_t *set) {
	return set->stage ? renderer->graph_lookups.buffer : set->graph_lookups.buffer;
}

//Every fetch has to be dispatched, the next one carries on from the particle state this one leaves in the pool
static tfxU32 tfx_impl_webgpu__simulate_gpu_particles(tfx_impl_webgpu_t *renderer, WGPUComputePassEncoder compute_pass, tfxU32 block) {
	tfxU32 current = renderer->current_instance_buffer;
	WGPUBindGroupEntry entries[5] = {
		tfx_impl_webgpu__buffer_entry(0, renderer->particle_pool.buffer, WGPU_WHOLE_SIZE),
		tfx_impl_webgpu__buffer_entry(1, renderer->instance_buffers[current].buffer, WGPU_WHOLE_SIZE),
		tfx_impl_webgpu__buffer_entry(2, renderer->graph_lookups.buffer, WGPU_WHOLE_SIZE),
		tfx_impl_webgpu__buffer_entry(3, renderer->particle_upload.buffer, WGPU_WHOLE_SIZE),
		tfx_impl_webgpu__buffer_entry(6, renderer->uniforms.buffer, sizeof(tfx_gpu_particle_push_t)),
	};
	WGPUBindGroup bind_group = tfx_impl_webgpu__bind_group(renderer, &renderer->particle_kernel_bind_groups[current], renderer->particle_kernel_layout, entries, 5);
	wgpuComputePassEncoderSetPipeline(compute_pass, renderer->particle_kernel);
	tfx_gpu_particle_push_t *push = &renderer->particle_push_constants;
	for (tfxU32 index = 0; index != renderer->particle_dispatch_count; ++index) {
		tfxU32 slot_count = renderer->particle_dispatches[index].slot_count;
		if (slot_count == 0) {
			continue;
		}
		//No push constants in WebGPU, so each dispatch reads its own block for its dispatch_index
		push->dispatch_index = index;
		memcpy(tfx_impl_webgpu__uniform_block(renderer, block), push, sizeof(tfx_gpu_particle_push_t));
		tfxU32 offset = block++ * TFX_IMPL_WEBGPU_UNIFORM_BLOCK_SIZE;
		wgpuComputePassEncoderSetBindGroup(compute_pass, 0, bind_group, 1, &offset);
		wgpuComputePassEncoderDispatchWorkgroups(compute_pass, (slot_count + 255) / 256, 1, 1);
	}
	renderer->particle_dispatch_count = 0;
	return block;
}

static tfxU32 tfx_impl_webgpu__play_back(tfx_impl_webgpu_t *renderer, WGPUComputePassEncoder compute_pass, WGPUComputePipeline pipeline, tfx_impl_webgpu_bind_group_t *cached,
	const WGPUBindGroupEntry *entries, const void *push, tfxU32 push_size, tfxU32 count, tfxU32 block) {
	WGPUBindGroup bind_group = tfx_impl_webgpu__bind_group(renderer, cached, renderer->playback_layout, entries, 6);
	memcpy(tfx_impl_webgpu__uniform_block(renderer, block), push, push_size);
	tfxU32 offset = block * TFX_IMPL_WEBGPU_UNIFORM_BLOCK_SIZE;
	wgpuComputePassEncoderSetPipeline(compute_pass, pipeline);
	wgpuComputePassEncoderSetBindGroup(compute_pass, 0, bind_group, 1, &offset);
	wgpuComputePassEncoderDispatchWorkgroups(compute_pass, (count + 127) / 128, 1, 1);
	return block + 1;
}

static tfxU32 tfx_impl_webgpu__build_ribbon_geometry(tfx_impl_webgpu_t *renderer, WGPUComputePassEncoder compute_pass, tfx_impl_webgpu_ribbon_set_t *set,
	const tfx_impl_uniforms_t *uniforms, tfxU32 block) {
	WGPUBindGroupEntry entries[7] = {
		tfx_impl_webgpu__buffer_entry(0, set->vertices.buffer, WGPU_WHOLE_SIZE),
		tfx_impl_webgpu__buffer_entry(1, set->indices.buffer, WGPU_WHOLE_SIZE),
		tfx_impl_webgpu__buffer_entry(2, set->emitters.buffer, WGPU_WHOLE_SIZE),
		tfx_impl_webgpu__buffer_entry(3, tfx_impl_webgpu__ribbon_graph_lookups(renderer, set), WGPU_WHOLE_SIZE),
		tfx_impl_webgpu__buffer_entry(4, set->instances.buffer, WGPU_WHOLE_SIZE),
		tfx_impl_webgpu__buffer_entry(5, set->segments.buffer, WGPU_WHOLE_SIZE),
		tfx_impl_webgpu__buffer_entry(6, renderer->uniforms.buffer, sizeof(tfx_ribbon_bucket_globals_t)),
	};
	WGPUBindGroup bind_group = tfx_impl_webgpu__bind_group(renderer, &set->compute_bind_group, renderer->ribbon_compute_layout, entries, 7);
	wgpuComputePassEncoderSetPipeline(compute_pass, renderer->ribbon_compute_pipeline);
	set->first_globals_block = block;

	//One dispatch per bucket, each with the bucket's own globals, which the ribbon draw reads from the same block
	tfx_ribbon_dispatch_t dispatch = tfx_impl_webgpu__create_ribbon_dispatch(set);
	while (tfx_impl_webgpu__next_ribbon_dispatch(set, &dispatch)) {
		tfx_ribbon_bucket_globals_t *globals = tfx_GetRibbonDispatchGlobals(&dispatch);
		if (set->stage) {
			globals->lerp = uniforms->timer_lerp;
		} else {
			//Baked ribbons are already interpolated, but they still face the camera
			globals->camera_position.x = uniforms->camera_position[0];
			globals->camera_position.y = uniforms->camera_position[1];
			globals->camera_position.z = uniforms->camera_position[2];
		}
		memcpy(tfx_impl_webgpu__uniform_block(renderer, block), globals, sizeof(tfx_ribbon_bucket_globals_t));
		tfxU32 offset = block++ * TFX_IMPL_WEBGPU_UNIFORM_BLOCK_SIZE;
		wgpuComputePassEncoderSetBindGroup(compute_pass, 0, bind_group, 1, &offset);
		wgpuComputePassEncoderDispatchWorkgroups(compute_pass, (dispatch.total_segments + TFX_IMPL_WEBGPU_RIBBON_GROUP_SIZE - 1) / TFX_IMPL_WEBGPU_RIBBON_GROUP_SIZE, 1, 1);
	}
	return block;
}

void tfx_impl_webgpu_Compute(tfx_impl_webgpu_t *renderer, WGPUCommandEncoder encoder, const tfx_impl_uniforms_t *uniforms) {
	//Every block has to be in place before a bind group can take the uniform buffer, so they're counted first
	tfxU32 block_count = 0;
	for (tfxU32 index = 0; index != renderer->particle_dispatch_count; ++index) {
		block_count += renderer->particle_dispatches[index].slot_count > 0 ? 1 : 0;
	}
	block_count += renderer->playback_sprite_count > 0 ? 1 : 0;
	block_count += renderer->playback_ribbon_count > 0 ? 1 : 0;
	block_count += tfx_impl_webgpu__count_ribbon_buckets(&renderer->live_ribbons);
	block_count += tfx_impl_webgpu__count_ribbon_buckets(&renderer->baked_ribbons);
	if (block_count == 0 || !tfx_impl_webgpu__ensure_uniform_blocks(renderer, TFX_IMPL_WEBGPU_RENDER_BLOCKS + block_count)) {
		renderer->particle_dispatch_count = 0;
		return;
	}

	//One pass: WebGPU puts a barrier between dispatches that write what a later one reads
	WGPUComputePassEncoder compute_pass = wgpuCommandEncoderBeginComputePass(encoder, NULL);
	tfxU32 block = TFX_IMPL_WEBGPU_RENDER_BLOCKS;
	if (renderer->particle_dispatch_count > 0) {
		block = tfx_impl_webgpu__simulate_gpu_particles(renderer, compute_pass, block);
	}
	if (renderer->playback_sprite_count > 0) {
		WGPUBindGroupEntry entries[6] = {
			tfx_impl_webgpu__buffer_entry(0, renderer->playback_instances.buffer, WGPU_WHOLE_SIZE),
			tfx_impl_webgpu__buffer_entry(2, renderer->animation_offsets.buffer, WGPU_WHOLE_SIZE),
			tfx_impl_webgpu__buffer_entry(3, renderer->animation_instances.buffer, WGPU_WHOLE_SIZE),
			tfx_impl_webgpu__buffer_entry(4, renderer->baked_sprites.buffer, WGPU_WHOLE_SIZE),
			tfx_impl_webgpu__buffer_entry(5, renderer->baked_emitter_properties.buffer, WGPU_WHOLE_SIZE),
			tfx_impl_webgpu__buffer_entry(6, renderer->uniforms.buffer, sizeof(tfx_sprite_data_push_t)),
		};
		block = tfx_impl_webgpu__play_back(renderer, compute_pass, renderer->sprite_playback_pipeline, &renderer->sprite_playback_bind_group, entries,
			&renderer->sprite_data_push, sizeof(tfx_sprite_data_push_t), renderer->playback_sprite_count, block);
	}
	if (renderer->playback_ribbon_count > 0) {
		WGPUBindGroupEntry entries[6] = {
			tfx_impl_webgpu__buffer_entry(0, renderer->baked_ribbons.instances.buffer, WGPU_WHOLE_SIZE),
			tfx_impl_webgpu__buffer_entry(2, renderer->animation_ribbon_offsets.buffer, WGPU_WHOLE_SIZE),
			tfx_impl_webgpu__buffer_entry(3, renderer->animation_instances.buffer, WGPU_WHOLE_SIZE),
			tfx_impl_webgpu__buffer_entry(4, renderer->baked_ribbon_data.buffer, WGPU_WHOLE_SIZE),
			tfx_impl_webgpu__buffer_entry(5, renderer->baked_ribbon_properties.buffer, WGPU_WHOLE_SIZE),
			tfx_impl_webgpu__buffer_entry(6, renderer->uniforms.buffer, sizeof(tfx_ribbon_data_push_t)),
		};
		block = tfx_impl_webgpu__play_back(renderer, compute_pass, renderer->ribbon_playback_pipeline, &renderer->ribbon_playback_bind_group, entries,
			&renderer->ribbon_data_push, sizeof(tfx_ribbon_data_push_t), renderer->playback_ribbon_count, block);
	}
	if (tfx_impl_webgpu__ribbon_set_active(&renderer->live_ribbons)) {
		block = tfx_impl_webgpu__build_ribbon_geometry(renderer, compute_pass, &renderer->live_ribbons, uniforms, block);
	}
	if (tfx_impl_webgpu__ribbon_set_active(&renderer->baked_ribbons)) {
		block = tfx_impl_webgpu__build_ribbon_geometry(renderer, compute_pass, &renderer->baked_ribbons, uniforms, block);
	}
	wgpuComputePassEncoderEnd(compute_pass);
	wgpuComputePassEncoderRelease(compute_pass);
	//Before the submit, so the queue writes the blocks ahead of the pass that reads them
	tfx_impl_webgpu__write_uniform_blocks(renderer, TFX_IMPL_WEBGPU_RENDER_BLOCKS, block);
}

//Live and baked ribbons each have their own shape frames and color ramps, the atlas is shared
//Only the buckets whose layer is in the group are drawn
static void tfx_impl_webgpu__render_ribbons(tfx_impl_webgpu_t *renderer, WGPURenderPassEncoder render_pass, tfx_impl_webgpu_ribbon_set_t *set,
	WGPUBuffer image_data, WGPUTextureView color_ramps_view, int writes_depth) {
	WGPUBindGroupEntry entries[11] = {
		tfx_impl_webgpu__buffer_entry(0, set->vertices.buffer, WGPU_WHOLE_SIZE),
		tfx_impl_webgpu__buffer_entry(1, set->instances.buffer, WGPU_WHOLE_SIZE),
		tfx_impl_webgpu__buffer_entry(2, set->emitters.buffer, WGPU_WHOLE_SIZE),
		tfx_impl_webgpu__buffer_entry(3, tfx_impl_webgpu__ribbon_graph_lookups(renderer, set), WGPU_WHOLE_SIZE),
		tfx_impl_webgpu__buffer_entry(4, image_data, WGPU_WHOLE_SIZE),
	};
	tfxU32 entry_count = tfx_impl_webgpu__add_render_entries(renderer, entries, 5, sizeof(tfx_ribbon_bucket_globals_t), color_ramps_view);
	WGPUBindGroup bind_group = tfx_impl_webgpu__bind_group(renderer, &set->render_bind_group, renderer->ribbon_layout, entries, entry_count);
	wgpuRenderPassEncoderSetIndexBuffer(render_pass, set->indices.buffer, WGPUIndexFormat_Uint32, 0, WGPU_WHOLE_SIZE);

	//The same walk as the compute, so each bucket finds the globals it wrote, skipped ones included
	WGPURenderPipeline bound_pipeline = NULL;
	tfxU32 block = set->first_globals_block;
	tfx_ribbon_dispatch_t dispatch = tfx_impl_webgpu__create_ribbon_dispatch(set);
	while (tfx_impl_webgpu__next_ribbon_dispatch(set, &dispatch)) {
		tfxU32 offset = block++ * TFX_IMPL_WEBGPU_UNIFORM_BLOCK_SIZE;
		const tfxLayerFlags flags = renderer->layer_flags[dispatch.layer];
		if (tfx_impl_LayerWritesDepth(flags) != writes_depth) {
			continue;
		}
		WGPURenderPipeline pipeline = renderer->ribbon_pipelines[writes_depth][tfx_impl_LayerIsLit(flags)];
		if (pipeline != bound_pipeline) {
			wgpuRenderPassEncoderSetPipeline(render_pass, pipeline);
			bound_pipeline = pipeline;
		}
		wgpuRenderPassEncoderSetBindGroup(render_pass, 0, bind_group, 1, &offset);
		wgpuRenderPassEncoderDrawIndexed(render_pass, dispatch.index_count, 1, dispatch.index_offset, 0, 0);
	}
}

static tfxU32 tfx_impl_webgpu__draw_instance_range(tfx_impl_webgpu_t *renderer, WGPURenderPassEncoder render_pass, WGPUBindGroup bind_group, tfxU32 start, tfxU32 count, tfxU32 block) {
	if (count == 0) {
		return block;
	}
	//The first instance goes in the draw constants, like the other fixed binding backends, though WGSL's instance_index would include it
	tfx_impl_webgpu_draw_constants_t *draw_constants = (tfx_impl_webgpu_draw_constants_t *)tfx_impl_webgpu__uniform_block(renderer, block);
	memset(draw_constants, 0, sizeof(tfx_impl_webgpu_draw_constants_t));
	draw_constants->instance_offset = start;
	tfxU32 offset = block * TFX_IMPL_WEBGPU_UNIFORM_BLOCK_SIZE;
	wgpuRenderPassEncoderSetBindGroup(render_pass, 0, bind_group, 1, &offset);
	wgpuRenderPassEncoderDraw(render_pass, 6, count, 0, 0);
	return block + 1;
}

//Each layer in the group as its CPU instances then its GPU instances
static tfxU32 tfx_impl_webgpu__render_billboards(tfx_impl_webgpu_t *renderer, WGPURenderPassEncoder render_pass, WGPUBindGroup bind_group, int writes_depth, tfxU32 block) {
	WGPURenderPipeline bound_pipeline = NULL;
	for (tfxU32 layer = 0; layer != tfxLAYERS; ++layer) {
		const tfxLayerFlags flags = renderer->layer_flags[layer];
		if (tfx_impl_LayerWritesDepth(flags) != writes_depth || renderer->layer_count[layer] + renderer->gpu_layer_count[layer] == 0) {
			continue;
		}
		WGPURenderPipeline pipeline = renderer->billboard_pipelines[writes_depth][tfx_impl_LayerIsLit(flags)];
		if (pipeline != bound_pipeline) {
			wgpuRenderPassEncoderSetPipeline(render_pass, pipeline);
			bound_pipeline = pipeline;
		}
		block = tfx_impl_webgpu__draw_instance_range(renderer, render_pass, bind_group, renderer->layer_start[layer], renderer->layer_count[layer], block);
		block = tfx_impl_webgpu__draw_instance_range(renderer, render_pass, bind_group, renderer->gpu_layer_start[layer], renderer->gpu_layer_count[layer], block);
	}
	return block;
}

void tfx_impl_webgpu_Render(tfx_impl_webgpu_t *renderer, WGPURenderPassEncoder render_pass, const tfx_impl_uniforms_t *uniforms) {
	if (!renderer->shape_atlas || (renderer->instance_count == 0 && renderer->playback_sprite_count == 0
		&& !tfx_impl_webgpu__ribbon_set_active(&renderer->live_ribbons) && !tfx_impl_webgpu__ribbon_set_active(&renderer->baked_ribbons))) {
		return;
	}
	//Every vertex shader reads the frame's uniforms from block 0
	memcpy(tfx_impl_webgpu__uniform_block(renderer, 0), uniforms, sizeof(tfx_impl_uniforms_t));
	tfxU32 block = 1;

	//The live draws need the library data, which goes up with the first update
	const bool live_data = renderer->color_ramps && renderer->image_data.buffer;
	WGPUBindGroup billboard_bind_group = NULL;
	if (live_data && renderer->instance_count > 0) {
		tfxU32 current = renderer->current_instance_buffer;
		WGPUBindGroupEntry entries[10] = {
			tfx_impl_webgpu__buffer_entry(0, renderer->instance_buffers[current].buffer, WGPU_WHOLE_SIZE),
			tfx_impl_webgpu__buffer_entry(1, renderer->instance_buffers[current ^ 1].buffer, WGPU_WHOLE_SIZE),
			tfx_impl_webgpu__buffer_entry(2, renderer->image_data.buffer, WGPU_WHOLE_SIZE),
			tfx_impl_webgpu__buffer_entry(3, renderer->emitter_properties.buffer, WGPU_WHOLE_SIZE),
		};
		tfxU32 entry_count = tfx_impl_webgpu__add_render_entries(renderer, entries, 4, sizeof(tfx_impl_webgpu_draw_constants_t), renderer->color_ramps_view);
		billboard_bind_group = tfx_impl_webgpu__bind_group(renderer, &renderer->billboard_bind_groups[current], renderer->billboard_layout, entries, entry_count);
	}
	//The layers that write depth draw first with blending off, so they hide the blended layers behind them
	for (int group = 0; group != 2; ++group) {
		const int writes_depth = group == 0;
		if (billboard_bind_group) {
			block = tfx_impl_webgpu__render_billboards(renderer, render_pass, billboard_bind_group, writes_depth, block);
		}
		if (live_data && tfx_impl_webgpu__ribbon_set_active(&renderer->live_ribbons)) {
			tfx_impl_webgpu__render_ribbons(renderer, render_pass, &renderer->live_ribbons, renderer->image_data.buffer, renderer->color_ramps_view, writes_depth);
		}
		if (!renderer->baked_color_ramps) {
			continue;
		}
		//Baked billboards come as one range across their layers, so they all draw blended and unlit
		if (!writes_depth && renderer->playback_sprite_count > 0) {
			//Baked sprites come with their own shape frames, emitter properties and color ramps, the atlas is shared
			WGPUBindGroupEntry entries[9] = {
				tfx_impl_webgpu__buffer_entry(0, renderer->playback_instances.buffer, WGPU_WHOLE_SIZE),
				tfx_impl_webgpu__buffer_entry(1, renderer->baked_image_data.buffer, WGPU_WHOLE_SIZE),
				tfx_impl_webgpu__buffer_entry(2, renderer->baked_emitter_properties.buffer, WGPU_WHOLE_SIZE),
			};
			tfxU32 entry_count = tfx_impl_webgpu__add_render_entries(renderer, entries, 3, sizeof(tfx_impl_webgpu_draw_constants_t), renderer->baked_color_ramps_view);
			WGPUBindGroup bind_group = tfx_impl_webgpu__bind_group(renderer, &renderer->baked_billboard_bind_group, renderer->baked_billboard_layout, entries, entry_count);
			wgpuRenderPassEncoderSetPipeline(render_pass, renderer->baked_billboard_pipeline);
			block = tfx_impl_webgpu__draw_instance_range(renderer, render_pass, bind_group, 0, renderer->playback_sprite_count, block);
		}
		if (tfx_impl_webgpu__ribbon_set_active(&renderer->baked_ribbons)) {
			tfx_impl_webgpu__render_ribbons(renderer, render_pass, &renderer->baked_ribbons, renderer->baked_image_data.buffer, renderer->baked_color_ramps_view, writes_depth);
		}
	}
	tfx_impl_webgpu__write_uniform_blocks(renderer, 0, block);
}

static void tfx_impl_webgpu__release_ribbon_set(tfx_impl_webgpu_ribbon_set_t *set) {
	tfx_impl_webgpu__release_buffer(&set->segments);
	tfx_impl_webgpu__release_buffer(&set->instances);
	tfx_impl_webgpu__release_buffer(&set->emitters);
	tfx_impl_webgpu__release_buffer(&set->vertices);
	tfx_impl_webgpu__release_buffer(&set->indices);
	tfx_impl_webgpu__release_buffer(&set->graph_lookups);
	tfx_impl_webgpu__release_bind_group(&set->compute_bind_group);
	tfx_impl_webgpu__release_bind_group(&set->render_bind_group);
}

void tfx_impl_webgpu_Shutdown(tfx_impl_webgpu_t *renderer) {
	if (!renderer->device) {
		return;
	}
	for (int index = 0; index != 2; ++index) {
		tfx_impl_webgpu__release_buffer(&renderer->instance_buffers[index]);
		tfx_impl_webgpu__release_bind_group(&renderer->billboard_bind_groups[index]);
		tfx_impl_webgpu__release_bind_group(&renderer->particle_kernel_bind_groups[index]);
	}
	tfx_impl_webgpu__release_bind_group(&renderer->sprite_playback_bind_group);
	tfx_impl_webgpu__release_bind_group(&renderer->baked_billboard_bind_group);
	tfx_impl_webgpu__release_bind_group(&renderer->ribbon_playback_bind_group);
	tfx_impl_webgpu__release_buffer(&renderer->uniforms);
	tfx_impl_webgpu__release_buffer(&renderer->image_data);
	tfx_impl_webgpu__release_buffer(&renderer->emitter_properties);
	tfx_impl_webgpu__release_buffer(&renderer->graph_lookups);
	tfx_impl_webgpu__release_ribbon_set(&renderer->live_ribbons);
	tfx_impl_webgpu__release_ribbon_set(&renderer->baked_ribbons);
	tfx_impl_webgpu__release_buffer(&renderer->particle_pool);
	tfx_impl_webgpu__release_buffer(&renderer->particle_upload);
	tfx_impl_webgpu__release_buffer(&renderer->baked_image_data);
	tfx_impl_webgpu__release_buffer(&renderer->baked_emitter_properties);
	tfx_impl_webgpu__release_buffer(&renderer->baked_sprites);
	tfx_impl_webgpu__release_buffer(&renderer->animation_offsets);
	tfx_impl_webgpu__release_buffer(&renderer->animation_instances);
	tfx_impl_webgpu__release_buffer(&renderer->playback_instances);
	tfx_impl_webgpu__release_buffer(&renderer->baked_ribbon_data);
	tfx_impl_webgpu__release_buffer(&renderer->baked_ribbon_properties);
	tfx_impl_webgpu__release_buffer(&renderer->animation_ribbon_offsets);
	tfx_impl_webgpu__release_texture(&renderer->shape_atlas, &renderer->shape_atlas_view);
	tfx_impl_webgpu__release_texture(&renderer->color_ramps, &renderer->color_ramps_view);
	tfx_impl_webgpu__release_texture(&renderer->baked_color_ramps, &renderer->baked_color_ramps_view);
	if (renderer->atlas_sampler) wgpuSamplerRelease(renderer->atlas_sampler);
	for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
		for (int lit = 0; lit != 2; ++lit) {
			if (renderer->billboard_pipelines[writes_depth][lit]) wgpuRenderPipelineRelease(renderer->billboard_pipelines[writes_depth][lit]);
			if (renderer->ribbon_pipelines[writes_depth][lit]) wgpuRenderPipelineRelease(renderer->ribbon_pipelines[writes_depth][lit]);
		}
	}
	if (renderer->baked_billboard_pipeline) wgpuRenderPipelineRelease(renderer->baked_billboard_pipeline);
	if (renderer->particle_kernel) wgpuComputePipelineRelease(renderer->particle_kernel);
	if (renderer->ribbon_compute_pipeline) wgpuComputePipelineRelease(renderer->ribbon_compute_pipeline);
	if (renderer->sprite_playback_pipeline) wgpuComputePipelineRelease(renderer->sprite_playback_pipeline);
	if (renderer->ribbon_playback_pipeline) wgpuComputePipelineRelease(renderer->ribbon_playback_pipeline);
	WGPUBindGroupLayout layouts[6] = { renderer->billboard_layout, renderer->baked_billboard_layout, renderer->ribbon_layout,
		renderer->particle_kernel_layout, renderer->ribbon_compute_layout, renderer->playback_layout };
	for (int index = 0; index != 6; ++index) {
		if (layouts[index]) {
			wgpuBindGroupLayoutRelease(layouts[index]);
		}
	}
	if (renderer->queue) wgpuQueueRelease(renderer->queue);
	free(renderer->color_ramp_layer_versions);
	free(renderer->uniform_staging);
	free(renderer->particle_upload_staging);
	memset(renderer, 0, sizeof(tfx_impl_webgpu_t));
}
