#include "tfx_impl_sdlgpu.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TFX_IMPL_SDLGPU_INITIAL_INSTANCE_CAPACITY 10000
//tfx_GetGlobalColorRampLayerCount never goes past this
#define TFX_IMPL_SDLGPU_MAX_COLOR_RAMP_LAYERS 256
//Matches tfx_ribbon_group_size in timelinefx.slang
#define TFX_IMPL_SDLGPU_RIBBON_GROUP_SIZE 256
//D3D12 places texture copies at 512 byte steps in the transfer buffer
#define TFX_IMPL_SDLGPU_TEXTURE_UPLOAD_ALIGNMENT 512

//timelinefx.slang mirrors these structs by hand, in its [Render_rows] section
static_assert(sizeof(tfx_instance_t) == 48, "tfx_instance_t changed, update tfx_instance in timelinefx.slang");
static_assert(sizeof(tfx_gpu_image_data_t) == 48, "tfx_gpu_image_data_t changed, update tfx_image_data in timelinefx.slang");
static_assert(sizeof(tfx_gpu_particle_properties_t) == 32, "tfx_gpu_particle_properties_t changed, update tfx_emitter_properties in timelinefx.slang");
static_assert(sizeof(tfx_ribbon_vertex_t) == 32, "tfx_ribbon_vertex_t changed, update tfx_ribbon_vertex_row in timelinefx.slang");
static_assert(sizeof(tfx_gpu_graph_data_t) == 48, "tfx_gpu_graph_data_t changed, update tfx_gpu_graph_data in timelinefx.slang");
static_assert(sizeof(tfx_impl_uniforms_t) % 16 == 0, "Uniform buffers must be a multiple of 16 bytes");

//Mirrors tfx_fixed_draw_constants in timelinefx.slang
typedef struct tfx_impl_sdlgpu_draw_constants_s {
	tfxU32 instance_offset;
	tfxU32 previous_index_offset;
	float ndc_offset_x;
	float ndc_offset_y;
} tfx_impl_sdlgpu_draw_constants_t;

//--- Sdlgpu_resources

static void tfx_impl_sdlgpu__release_buffer(SDL_GPUDevice *device, tfx_impl_sdlgpu_buffer_t *buffer) {
	if (buffer->buffer) {
		SDL_ReleaseGPUBuffer(device, buffer->buffer);
	}
	buffer->buffer = NULL;
	buffer->capacity = 0;
}

static bool tfx_impl_sdlgpu__create_buffer(SDL_GPUDevice *device, tfxU32 size_in_bytes, SDL_GPUBufferUsageFlags usage, tfx_impl_sdlgpu_buffer_t *buffer) {
	tfx_impl_sdlgpu__release_buffer(device, buffer);
	size_in_bytes = (size_in_bytes + 15) & ~15u;
	SDL_GPUBufferCreateInfo create_info = {};
	create_info.usage = usage;
	create_info.size = size_in_bytes;
	buffer->buffer = SDL_CreateGPUBuffer(device, &create_info);
	if (!buffer->buffer) {
		tfx_Print("Could not create a buffer of %u bytes: %s", size_in_bytes, SDL_GetError());
		return false;
	}
	buffer->capacity = size_in_bytes;
	return true;
}

//Grows by half again whenever it's too small, so the size settles after the first few frames. Never empty, so a shader always has a buffer bound
static bool tfx_impl_sdlgpu__ensure_buffer(SDL_GPUDevice *device, tfxU32 size_in_bytes, SDL_GPUBufferUsageFlags usage, tfx_impl_sdlgpu_buffer_t *buffer) {
	size_in_bytes = size_in_bytes ? size_in_bytes : 16;
	if (buffer->buffer && buffer->capacity >= size_in_bytes) {
		return true;
	}
	return tfx_impl_sdlgpu__create_buffer(device, size_in_bytes + size_in_bytes / 2, usage, buffer);
}

static tfx_impl_sdlgpu_upload_t *tfx_impl_sdlgpu__add_upload(tfx_impl_sdlgpu_t *renderer, const void *data, tfxU32 size) {
	if (renderer->upload_count == renderer->upload_capacity) {
		tfxU32 capacity = renderer->upload_capacity ? renderer->upload_capacity * 2 : 64;
		tfx_impl_sdlgpu_upload_t *uploads = (tfx_impl_sdlgpu_upload_t *)realloc(renderer->uploads, sizeof(tfx_impl_sdlgpu_upload_t) * capacity);
		if (!uploads) {
			return NULL;
		}
		renderer->uploads = uploads;
		renderer->upload_capacity = capacity;
	}
	tfx_impl_sdlgpu_upload_t *upload = &renderer->uploads[renderer->upload_count++];
	memset(upload, 0, sizeof(tfx_impl_sdlgpu_upload_t));
	upload->data = data;
	upload->size = size;
	return upload;
}

//The data has to stay valid until the next tfx_impl_sdlgpu__flush_uploads
static bool tfx_impl_sdlgpu__upload_buffer(tfx_impl_sdlgpu_t *renderer, const void *data, tfxU32 size_in_bytes, SDL_GPUBuffer *buffer) {
	if (size_in_bytes == 0) {
		return true;
	}
	tfx_impl_sdlgpu_upload_t *upload = tfx_impl_sdlgpu__add_upload(renderer, data, size_in_bytes);
	if (!upload) {
		return false;
	}
	upload->buffer = buffer;
	return true;
}

static bool tfx_impl_sdlgpu__upload_texture(tfx_impl_sdlgpu_t *renderer, const void *data, tfxU32 size_in_bytes, SDL_GPUTexture *texture, tfxU32 layer, tfxU32 mip_level, tfxU32 width, tfxU32 height) {
	tfx_impl_sdlgpu_upload_t *upload = tfx_impl_sdlgpu__add_upload(renderer, data, size_in_bytes);
	if (!upload) {
		return false;
	}
	upload->texture = texture;
	upload->layer = layer;
	upload->mip_level = mip_level;
	upload->width = width;
	upload->height = height;
	return true;
}

//Copies every queued upload into the transfer buffer and records them in one copy pass, grow copies first
static bool tfx_impl_sdlgpu__flush_uploads(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer) {
	SDL_GPUDevice *device = renderer->device;
	bool flushed = true;
	if (renderer->upload_count == 0 && renderer->copy_count == 0) {
		return true;
	}
	tfxU32 transfer_size = 0;
	for (tfxU32 index = 0; index != renderer->upload_count; ++index) {
		tfx_impl_sdlgpu_upload_t *upload = &renderer->uploads[index];
		tfxU32 alignment = upload->texture ? TFX_IMPL_SDLGPU_TEXTURE_UPLOAD_ALIGNMENT : 16;
		transfer_size = (transfer_size + alignment - 1) & ~(alignment - 1);
		upload->transfer_offset = transfer_size;
		transfer_size += upload->size;
	}
	if (transfer_size > renderer->transfer_capacity) {
		if (renderer->transfer_buffer) {
			SDL_ReleaseGPUTransferBuffer(device, renderer->transfer_buffer);
		}
		SDL_GPUTransferBufferCreateInfo create_info = {};
		create_info.usage = SDL_GPU_TRANSFERBUFFERUSAGE_UPLOAD;
		create_info.size = transfer_size + transfer_size / 2;
		renderer->transfer_buffer = SDL_CreateGPUTransferBuffer(device, &create_info);
		renderer->transfer_capacity = renderer->transfer_buffer ? create_info.size : 0;
		if (!renderer->transfer_buffer) {
			tfx_Print("Could not create a transfer buffer of %u bytes: %s", create_info.size, SDL_GetError());
			flushed = false;
		}
	}
	//Cycled, so a transfer buffer the GPU is still reading from an earlier frame isn't overwritten
	unsigned char *mapped = flushed && transfer_size ? (unsigned char *)SDL_MapGPUTransferBuffer(device, renderer->transfer_buffer, true) : NULL;
	if (flushed && transfer_size && !mapped) {
		tfx_Print("Could not map the transfer buffer: %s", SDL_GetError());
		flushed = false;
	}
	if (flushed) {
		for (tfxU32 index = 0; index != renderer->upload_count; ++index) {
			const tfx_impl_sdlgpu_upload_t *upload = &renderer->uploads[index];
			if (upload->data) {
				memcpy(mapped + upload->transfer_offset, upload->data, upload->size);
			} else {
				tfx_CopyGPUParticleUpload(renderer->particle_frame, mapped + upload->transfer_offset);
			}
		}
		if (mapped) {
			SDL_UnmapGPUTransferBuffer(device, renderer->transfer_buffer);
		}
		SDL_GPUCopyPass *copy_pass = SDL_BeginGPUCopyPass(command_buffer);
		for (tfxU32 index = 0; index != renderer->copy_count; ++index) {
			const tfx_impl_sdlgpu_copy_t *copy = &renderer->copies[index];
			SDL_GPUBufferLocation source = { copy->source, 0 };
			SDL_GPUBufferLocation destination = { copy->destination, 0 };
			SDL_CopyGPUBufferToBuffer(copy_pass, &source, &destination, copy->size, false);
		}
		for (tfxU32 index = 0; index != renderer->upload_count; ++index) {
			const tfx_impl_sdlgpu_upload_t *upload = &renderer->uploads[index];
			if (upload->texture) {
				SDL_GPUTextureTransferInfo source = { renderer->transfer_buffer, upload->transfer_offset, upload->width, upload->height };
				SDL_GPUTextureRegion destination = {};
				destination.texture = upload->texture;
				destination.mip_level = upload->mip_level;
				destination.layer = upload->layer;
				destination.w = upload->width;
				destination.h = upload->height;
				destination.d = 1;
				SDL_UploadToGPUTexture(copy_pass, &source, &destination, false);
			} else {
				SDL_GPUTransferBufferLocation source = { renderer->transfer_buffer, upload->transfer_offset };
				SDL_GPUBufferRegion destination = { upload->buffer, upload->buffer_offset, upload->size };
				//Never cycled: the previous instance buffer and the particle pool have to keep what's in them
				SDL_UploadToGPUBuffer(copy_pass, &source, &destination, false);
			}
		}
		SDL_EndGPUCopyPass(copy_pass);
	}
	for (tfxU32 index = 0; index != renderer->copy_count; ++index) {
		SDL_ReleaseGPUBuffer(device, renderer->copies[index].source);
	}
	renderer->copy_count = 0;
	renderer->upload_count = 0;
	free(renderer->upload_scratch);
	renderer->upload_scratch = NULL;
	return flushed;
}

//slangc compiles every shader ahead of time to <entry point>.spv or .dxil in the shader folder. The caller frees the code
static void *tfx_impl_sdlgpu__load_shader_code(const char *shader_folder, const char *entry_point, SDL_GPUShaderFormat format, size_t *size) {
	char path[1024];
	snprintf(path, sizeof(path), "%s/%s.%s", shader_folder, entry_point, format == SDL_GPU_SHADERFORMAT_DXIL ? "dxil" : "spv");
	void *code = SDL_LoadFile(path, size);
	if (!code) {
		tfx_Print("Could not load the shader %s: %s", path, SDL_GetError());
	}
	return code;
}

static SDL_GPUShaderFormat tfx_impl_sdlgpu__shader_format(SDL_GPUDevice *device) {
	SDL_GPUShaderFormat formats = SDL_GetGPUShaderFormats(device);
	return (formats & SDL_GPU_SHADERFORMAT_SPIRV) ? SDL_GPU_SHADERFORMAT_SPIRV : SDL_GPU_SHADERFORMAT_DXIL;
}

//The counts are each entry point's slots in the Shaders readme. SPIR-V keeps the entry point's name because slangc is given -fvk-use-entrypoint-name
static SDL_GPUShader *tfx_impl_sdlgpu__create_shader(SDL_GPUDevice *device, const char *shader_folder, const char *entry_point, SDL_GPUShaderStage stage,
	tfxU32 sampler_count, tfxU32 storage_buffer_count, tfxU32 uniform_buffer_count) {
	SDL_GPUShaderFormat format = tfx_impl_sdlgpu__shader_format(device);
	size_t size = 0;
	void *code = tfx_impl_sdlgpu__load_shader_code(shader_folder, entry_point, format, &size);
	if (!code) {
		return NULL;
	}
	SDL_GPUShaderCreateInfo create_info = {};
	create_info.code_size = size;
	create_info.code = (const Uint8 *)code;
	create_info.entrypoint = entry_point;
	create_info.format = format;
	create_info.stage = stage;
	create_info.num_samplers = sampler_count;
	create_info.num_storage_buffers = storage_buffer_count;
	create_info.num_uniform_buffers = uniform_buffer_count;
	SDL_GPUShader *shader = SDL_CreateGPUShader(device, &create_info);
	SDL_free(code);
	if (!shader) {
		tfx_Print("Could not create the shader %s: %s", entry_point, SDL_GetError());
	}
	return shader;
}

static SDL_GPUComputePipeline *tfx_impl_sdlgpu__create_compute_pipeline(SDL_GPUDevice *device, const char *shader_folder, const char *entry_point,
	tfxU32 read_count, tfxU32 write_count, tfxU32 thread_count) {
	SDL_GPUShaderFormat format = tfx_impl_sdlgpu__shader_format(device);
	size_t size = 0;
	void *code = tfx_impl_sdlgpu__load_shader_code(shader_folder, entry_point, format, &size);
	if (!code) {
		return NULL;
	}
	SDL_GPUComputePipelineCreateInfo create_info = {};
	create_info.code_size = size;
	create_info.code = (const Uint8 *)code;
	create_info.entrypoint = entry_point;
	create_info.format = format;
	create_info.num_readonly_storage_buffers = read_count;
	create_info.num_readwrite_storage_buffers = write_count;
	create_info.num_uniform_buffers = 1;
	create_info.threadcount_x = thread_count;
	create_info.threadcount_y = 1;
	create_info.threadcount_z = 1;
	SDL_GPUComputePipeline *pipeline = SDL_CreateGPUComputePipeline(device, &create_info);
	SDL_free(code);
	if (!pipeline) {
		tfx_Print("Could not create the compute pipeline %s: %s", entry_point, SDL_GetError());
	}
	return pipeline;
}

//Premultiplied blending, depth tested but not written, and no vertex input: instances and ribbon vertexes are read from storage buffers
//Premultiplied blending, depth tested but not written, or for a layer that writes depth, blending off and depth written
static SDL_GPUGraphicsPipeline *tfx_impl_sdlgpu__create_graphics_pipeline(SDL_GPUDevice *device, SDL_GPUShader *vertex_shader, SDL_GPUShader *fragment_shader,
	SDL_GPUTextureFormat color_format, SDL_GPUTextureFormat depth_format, bool writes_depth) {
	if (!vertex_shader || !fragment_shader) {
		return NULL;
	}
	SDL_GPUColorTargetDescription color_target = {};
	color_target.format = color_format;
	color_target.blend_state.enable_blend = !writes_depth;
	color_target.blend_state.src_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
	color_target.blend_state.dst_color_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
	color_target.blend_state.color_blend_op = SDL_GPU_BLENDOP_ADD;
	color_target.blend_state.src_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE;
	color_target.blend_state.dst_alpha_blendfactor = SDL_GPU_BLENDFACTOR_ONE_MINUS_SRC_ALPHA;
	color_target.blend_state.alpha_blend_op = SDL_GPU_BLENDOP_ADD;
	SDL_GPUGraphicsPipelineCreateInfo create_info = {};
	create_info.vertex_shader = vertex_shader;
	create_info.fragment_shader = fragment_shader;
	create_info.primitive_type = SDL_GPU_PRIMITIVETYPE_TRIANGLELIST;
	create_info.rasterizer_state.fill_mode = SDL_GPU_FILLMODE_FILL;
	create_info.rasterizer_state.cull_mode = SDL_GPU_CULLMODE_NONE;
	create_info.rasterizer_state.enable_depth_clip = true;
	create_info.depth_stencil_state.enable_depth_test = true;
	create_info.depth_stencil_state.enable_depth_write = writes_depth;
	create_info.depth_stencil_state.compare_op = SDL_GPU_COMPAREOP_LESS_OR_EQUAL;
	create_info.target_info.color_target_descriptions = &color_target;
	create_info.target_info.num_color_targets = 1;
	create_info.target_info.depth_stencil_format = depth_format;
	create_info.target_info.has_depth_stencil_target = true;
	SDL_GPUGraphicsPipeline *pipeline = SDL_CreateGPUGraphicsPipeline(device, &create_info);
	if (!pipeline) {
		tfx_Print("Could not create a graphics pipeline: %s", SDL_GetError());
	}
	return pipeline;
}

static SDL_GPUTextureFormat tfx_impl_sdlgpu__color_ramp_format(tfxU32 bytes_per_pixel) {
	if (bytes_per_pixel == 4) {
		return SDL_GPU_TEXTUREFORMAT_R8G8B8A8_UNORM;
	} else if (bytes_per_pixel == 16) {
		return SDL_GPU_TEXTUREFORMAT_R32G32B32A32_FLOAT;
	}
	return SDL_GPU_TEXTUREFORMAT_R16G16B16A16_FLOAT;
}

static SDL_GPUTexture *tfx_impl_sdlgpu__create_color_ramp_texture(SDL_GPUDevice *device, tfxU32 width, tfxU32 height, tfxU32 bytes_per_pixel, tfxU32 layer_count) {
	SDL_GPUTextureCreateInfo create_info = {};
	create_info.type = SDL_GPU_TEXTURETYPE_2D_ARRAY;
	create_info.format = tfx_impl_sdlgpu__color_ramp_format(bytes_per_pixel);
	create_info.usage = SDL_GPU_TEXTUREUSAGE_SAMPLER;
	create_info.width = width;
	create_info.height = height;
	create_info.layer_count_or_depth = layer_count;
	create_info.num_levels = 1;
	SDL_GPUTexture *texture = SDL_CreateGPUTexture(device, &create_info);
	if (!texture) {
		tfx_Print("Could not create a color ramp texture of %u layers: %s", layer_count, SDL_GetError());
	}
	return texture;
}

//A 2 x 2 box filter, the same as the mips D3D11 generates
static void tfx_impl_sdlgpu__downsample(const unsigned char *source, tfxU32 source_size, unsigned char *destination) {
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

//--- Sdlgpu_functions

bool tfx_impl_sdlgpu_Init(tfx_impl_sdlgpu_t *renderer, SDL_GPUDevice *device, SDL_GPUTextureFormat color_format, SDL_GPUTextureFormat depth_format, const char *shader_folder) {
	memset(renderer, 0, sizeof(tfx_impl_sdlgpu_t));
	renderer->device = device;

	//The lit fragments take the frame's uniforms in fragment slot 0 for the light
	SDL_GPUShader *billboard_vertex_shaders[2];
	SDL_GPUShader *ribbon_vertex_shaders[2];
	SDL_GPUShader *billboard_fragment_shaders[2][2];
	SDL_GPUShader *ribbon_fragment_shaders[2][2];
	for (int lit = 0; lit != 2; ++lit) {
		billboard_vertex_shaders[lit] = tfx_impl_sdlgpu__create_shader(device, shader_folder, tfx_impl_billboard_vertex_entry_points[lit], SDL_GPU_SHADERSTAGE_VERTEX, 0, 4, 2);
		ribbon_vertex_shaders[lit] = tfx_impl_sdlgpu__create_shader(device, shader_folder, tfx_impl_ribbon_vertex_entry_points[lit], SDL_GPU_SHADERSTAGE_VERTEX, 0, 5, 2);
		for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
			billboard_fragment_shaders[writes_depth][lit] = tfx_impl_sdlgpu__create_shader(device, shader_folder, tfx_impl_billboard_fragment_entry_points[writes_depth][lit], SDL_GPU_SHADERSTAGE_FRAGMENT, 2, 0, lit);
			ribbon_fragment_shaders[writes_depth][lit] = tfx_impl_sdlgpu__create_shader(device, shader_folder, tfx_impl_ribbon_fragment_entry_points[writes_depth][lit], SDL_GPU_SHADERSTAGE_FRAGMENT, 2, 0, lit);
		}
	}
	SDL_GPUShader *baked_billboard_vertex = tfx_impl_sdlgpu__create_shader(device, shader_folder, "tfx_billboard_baked_vertex", SDL_GPU_SHADERSTAGE_VERTEX, 0, 3, 2);
	bool pipelines_created = true;
	for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
		for (int lit = 0; lit != 2; ++lit) {
			renderer->billboard_pipelines[writes_depth][lit] = tfx_impl_sdlgpu__create_graphics_pipeline(device, billboard_vertex_shaders[lit], billboard_fragment_shaders[writes_depth][lit],
				color_format, depth_format, writes_depth != 0);
			renderer->ribbon_pipelines[writes_depth][lit] = tfx_impl_sdlgpu__create_graphics_pipeline(device, ribbon_vertex_shaders[lit], ribbon_fragment_shaders[writes_depth][lit],
				color_format, depth_format, writes_depth != 0);
			pipelines_created = pipelines_created && renderer->billboard_pipelines[writes_depth][lit] && renderer->ribbon_pipelines[writes_depth][lit];
		}
	}
	renderer->baked_billboard_pipeline = tfx_impl_sdlgpu__create_graphics_pipeline(device, baked_billboard_vertex, billboard_fragment_shaders[0][0], color_format, depth_format, false);
	//The pipelines keep what they need from the shaders
	SDL_GPUShader *shaders[13] = { baked_billboard_vertex, billboard_vertex_shaders[0], billboard_vertex_shaders[1], ribbon_vertex_shaders[0], ribbon_vertex_shaders[1],
		billboard_fragment_shaders[0][0], billboard_fragment_shaders[0][1], billboard_fragment_shaders[1][0], billboard_fragment_shaders[1][1],
		ribbon_fragment_shaders[0][0], ribbon_fragment_shaders[0][1], ribbon_fragment_shaders[1][0], ribbon_fragment_shaders[1][1] };
	for (int index = 0; index != 13; ++index) {
		if (shaders[index]) {
			SDL_ReleaseGPUShader(device, shaders[index]);
		}
	}
	renderer->particle_kernel = tfx_impl_sdlgpu__create_compute_pipeline(device, shader_folder, "simulate_particles", 2, 2, 256);
	renderer->ribbon_compute_pipeline = tfx_impl_sdlgpu__create_compute_pipeline(device, shader_folder, "tfx_ribbon_compute", 4, 2, TFX_IMPL_SDLGPU_RIBBON_GROUP_SIZE);
	renderer->sprite_playback_pipeline = tfx_impl_sdlgpu__create_compute_pipeline(device, shader_folder, "tfx_sprite_playback", 4, 1, 128);
	renderer->ribbon_playback_pipeline = tfx_impl_sdlgpu__create_compute_pipeline(device, shader_folder, "tfx_ribbon_playback", 4, 1, 128);

	SDL_GPUSamplerCreateInfo sampler_info = {};
	sampler_info.min_filter = SDL_GPU_FILTER_LINEAR;
	sampler_info.mag_filter = SDL_GPU_FILTER_LINEAR;
	sampler_info.mipmap_mode = SDL_GPU_SAMPLERMIPMAPMODE_LINEAR;
	sampler_info.address_mode_u = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
	sampler_info.address_mode_v = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
	sampler_info.address_mode_w = SDL_GPU_SAMPLERADDRESSMODE_CLAMP_TO_EDGE;
	sampler_info.max_lod = 1000.f;
	renderer->atlas_sampler = SDL_CreateGPUSampler(device, &sampler_info);
	sampler_info.min_filter = SDL_GPU_FILTER_NEAREST;
	sampler_info.mag_filter = SDL_GPU_FILTER_NEAREST;
	sampler_info.mipmap_mode = SDL_GPU_SAMPLERMIPMAPMODE_NEAREST;
	sampler_info.max_lod = 0.f;
	renderer->color_ramp_sampler = SDL_CreateGPUSampler(device, &sampler_info);

	if (!pipelines_created || !renderer->baked_billboard_pipeline || !renderer->particle_kernel || !renderer->ribbon_compute_pipeline
		|| !renderer->sprite_playback_pipeline || !renderer->ribbon_playback_pipeline || !renderer->atlas_sampler || !renderer->color_ramp_sampler) {
		tfx_Print("Could not create the TimelineFX pipelines and samplers");
		return false;
	}
	return true;
}

bool tfx_impl_sdlgpu_CreateAtlas(tfx_impl_sdlgpu_t *renderer, tfxU32 layer_size, tfxU32 layer_count) {
	if (renderer->shape_atlas) {
		SDL_ReleaseGPUTexture(renderer->device, renderer->shape_atlas);
	}
	tfxU32 mip_count = 1;
	while ((layer_size >> mip_count) > 0) {
		mip_count++;
	}
	SDL_GPUTextureCreateInfo create_info = {};
	create_info.type = SDL_GPU_TEXTURETYPE_2D_ARRAY;
	create_info.format = SDL_GPU_TEXTUREFORMAT_R8G8B8A8_UNORM;
	create_info.usage = SDL_GPU_TEXTUREUSAGE_SAMPLER;
	create_info.width = layer_size;
	create_info.height = layer_size;
	create_info.layer_count_or_depth = layer_count;
	create_info.num_levels = mip_count;
	renderer->shape_atlas = SDL_CreateGPUTexture(renderer->device, &create_info);
	if (!renderer->shape_atlas) {
		tfx_Print("Could not create the %u x %u x %u shape atlas: %s", layer_size, layer_size, layer_count, SDL_GetError());
		return false;
	}
	renderer->shape_atlas_layer_size = layer_size;
	renderer->shape_atlas_layer_count = layer_count;
	renderer->shape_atlas_mip_count = mip_count;
	return true;
}

bool tfx_impl_sdlgpu_UploadAtlas(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer, const tfx_impl_atlas_t *atlas) {
	if (!renderer->shape_atlas || atlas->layer_size != renderer->shape_atlas_layer_size
		|| atlas->first_layer + atlas->layer_count > renderer->shape_atlas_layer_count) {
		tfx_Print("An atlas of %u layers from layer %u, %u pixels square, doesn't fit the %u layer, %u pixel texture array",
			atlas->layer_count, atlas->first_layer, atlas->layer_size, renderer->shape_atlas_layer_count, renderer->shape_atlas_layer_size);
		return false;
	}
	//SDL_GenerateMipmapsForGPUTexture does every layer, and its Vulkan and D3D12 mips differ enough to show on particles a pixel or two across
	size_t layer_bytes = (size_t)atlas->layer_size * atlas->layer_size * 4;
	unsigned char *mips = (unsigned char *)malloc(layer_bytes / 2);
	if (!mips) {
		return false;
	}
	bool uploaded = true;
	for (tfxU32 layer = 0; layer != atlas->layer_count && uploaded; ++layer) {
		const unsigned char *source = atlas->pixels + layer * layer_bytes;
		tfxU32 size = atlas->layer_size;
		uploaded = tfx_impl_sdlgpu__upload_texture(renderer, source, (tfxU32)layer_bytes, renderer->shape_atlas, atlas->first_layer + layer, 0, size, size);
		unsigned char *destination = mips;
		for (tfxU32 mip_level = 1; mip_level < renderer->shape_atlas_mip_count && uploaded; ++mip_level) {
			tfx_impl_sdlgpu__downsample(source, size, destination);
			size = size > 1 ? size / 2 : 1;
			uploaded = tfx_impl_sdlgpu__upload_texture(renderer, destination, size * size * 4, renderer->shape_atlas, atlas->first_layer + layer, mip_level, size, size);
			source = destination;
			destination += (size_t)size * size * 4;
		}
		//One layer at a time, so the mips only ever need a third of a layer
		uploaded = tfx_impl_sdlgpu__flush_uploads(renderer, command_buffer) && uploaded;
	}
	free(mips);
	return uploaded;
}

bool tfx_impl_sdlgpu_UploadSpriteData(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer, tfx_animation_manager animation_manager, tfx_gpu_shapes shapes) {
	SDL_GPUDevice *device = renderer->device;
	tfxU32 bitmap_count = tfx_GetAnimationColorRampBitmapCount(animation_manager);
	if (renderer->baked_color_ramps) {
		SDL_ReleaseGPUTexture(device, renderer->baked_color_ramps);
		renderer->baked_color_ramps = NULL;
	}
	if (bitmap_count > 0) {
		tfx_bitmap_t *first_bitmap = tfx_GetAnimationColorRampBitmap(animation_manager, 0);
		tfxU32 width = (tfxU32)tfx_GetBitmapWidth(first_bitmap);
		tfxU32 height = (tfxU32)tfx_GetBitmapHeight(first_bitmap);
		tfxU32 bytes_per_pixel = (tfxU32)(tfx_GetBitmapSize(first_bitmap) / ((size_t)width * height));
		renderer->baked_color_ramps = tfx_impl_sdlgpu__create_color_ramp_texture(device, width, height, bytes_per_pixel, bitmap_count);
		if (!renderer->baked_color_ramps) {
			return false;
		}
		for (tfxU32 index = 0; index != bitmap_count; ++index) {
			tfx_bitmap_t *bitmap = tfx_GetAnimationColorRampBitmap(animation_manager, index);
			if (!tfx_impl_sdlgpu__upload_texture(renderer, tfx_GetBitmapData(bitmap), (tfxU32)tfx_GetBitmapSize(bitmap), renderer->baked_color_ramps, index, 0, width, height)) {
				return false;
			}
		}
	}

	tfxU32 shape_bytes = (tfxU32)sizeof(tfx_gpu_image_data_t) * tfx_GetGPUShapesCount(shapes);
	tfxU32 property_bytes = (tfxU32)tfx_GetAnimationEmitterPropertySizeInBytes(animation_manager);
	tfxU32 sprite_bytes = (tfxU32)tfx_GetSpriteDataSizeInBytes(animation_manager);
	if (!tfx_impl_sdlgpu__create_buffer(device, shape_bytes ? shape_bytes : 16, SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ, &renderer->baked_image_data)
		|| !tfx_impl_sdlgpu__create_buffer(device, property_bytes ? property_bytes : 16, SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ | SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_READ, &renderer->baked_emitter_properties)
		|| !tfx_impl_sdlgpu__create_buffer(device, sprite_bytes ? sprite_bytes : 16, SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_READ, &renderer->baked_sprites)
		|| !tfx_impl_sdlgpu__upload_buffer(renderer, tfx_GetGPUShapesArray(shapes), shape_bytes, renderer->baked_image_data.buffer)
		|| !tfx_impl_sdlgpu__upload_buffer(renderer, tfx_GetAnimationEmitterPropertiesBufferPointer(animation_manager), property_bytes, renderer->baked_emitter_properties.buffer)
		|| !tfx_impl_sdlgpu__upload_buffer(renderer, tfx_GetSpriteDataBufferPointer(animation_manager), sprite_bytes, renderer->baked_sprites.buffer)) {
		return false;
	}
	if (!tfx_AnimationManagerHasRibbons(animation_manager)) {
		return tfx_impl_sdlgpu__flush_uploads(renderer, command_buffer);
	}

	//The baked ribbons, their properties, segments, emitters and graph lookups only change when sprite data is loaded
	tfxU32 graph_bytes = (tfxU32)tfx_GetAnimationRibbonGraphLookupsSizeInBytes(animation_manager);
	tfxU32 emitters_bytes = (tfxU32)tfx_GetAnimationRibbonEmittersSizeInBytes(animation_manager);
	tfxU32 ribbon_data_bytes = (tfxU32)tfx_GetAnimationRibbonDataSizeInBytes(animation_manager);
	tfxU32 ribbon_property_bytes = (tfxU32)tfx_GetAnimationRibbonPropertySizeInBytes(animation_manager);
	tfxU32 segments_bytes = (tfxU32)tfx_GetAnimationRibbonSegmentsSizeInBytes(animation_manager);
	void *emitters = malloc(emitters_bytes ? emitters_bytes : 1);
	if (!emitters) {
		return false;
	}
	tfx_CopyAnimationRibbonEmitters(animation_manager, emitters);
	renderer->upload_scratch = emitters;
	SDL_GPUBufferUsageFlags ribbon_read = SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_READ | SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ;
	bool created = tfx_impl_sdlgpu__create_buffer(device, graph_bytes ? graph_bytes : 16, ribbon_read, &renderer->baked_ribbons.graph_lookups)
		&& tfx_impl_sdlgpu__create_buffer(device, ribbon_data_bytes, SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_READ, &renderer->baked_ribbon_data)
		&& tfx_impl_sdlgpu__create_buffer(device, ribbon_property_bytes, SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_READ, &renderer->baked_ribbon_properties)
		&& tfx_impl_sdlgpu__create_buffer(device, segments_bytes, SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_READ, &renderer->baked_ribbons.segments)
		&& tfx_impl_sdlgpu__create_buffer(device, emitters_bytes, ribbon_read, &renderer->baked_ribbons.emitters)
		&& tfx_impl_sdlgpu__upload_buffer(renderer, tfx_GetAnimationRibbonGraphLookupsBufferPointer(animation_manager), graph_bytes, renderer->baked_ribbons.graph_lookups.buffer)
		&& tfx_impl_sdlgpu__upload_buffer(renderer, tfx_GetAnimationRibbonDataBufferPointer(animation_manager), ribbon_data_bytes, renderer->baked_ribbon_data.buffer)
		&& tfx_impl_sdlgpu__upload_buffer(renderer, tfx_GetAnimationRibbonPropertiesBufferPointer(animation_manager), ribbon_property_bytes, renderer->baked_ribbon_properties.buffer)
		&& tfx_impl_sdlgpu__upload_buffer(renderer, tfx_GetAnimationRibbonSegmentsBufferPointer(animation_manager), segments_bytes, renderer->baked_ribbons.segments.buffer)
		&& tfx_impl_sdlgpu__upload_buffer(renderer, emitters, emitters_bytes, renderer->baked_ribbons.emitters.buffer);
	return tfx_impl_sdlgpu__flush_uploads(renderer, command_buffer) && created;
}

//The vertex and index buffers ribbon_compute writes, sized for this frame's buckets
static bool tfx_impl_sdlgpu__ensure_ribbon_geometry(SDL_GPUDevice *device, const tfx_ribbon_batch_sizes_t *sizes, tfx_impl_sdlgpu_ribbon_set_t *set) {
	return tfx_impl_sdlgpu__ensure_buffer(device, (tfxU32)sizes->vertex_buffer_size_in_bytes, SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_WRITE | SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ, &set->vertices)
		&& tfx_impl_sdlgpu__ensure_buffer(device, (tfxU32)sizes->index_buffer_size_in_bytes, SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_WRITE | SDL_GPU_BUFFERUSAGE_INDEX, &set->indices);
}

static bool tfx_impl_sdlgpu__upload_baked_ribbons(tfx_impl_sdlgpu_t *renderer, const tfx_impl_baked_frame_t *baked) {
	if (baked->ribbon_count == 0 || !renderer->baked_ribbon_data.buffer) {
		return true;
	}
	SDL_GPUDevice *device = renderer->device;
	tfx_ribbon_batch_sizes_t sizes = tfx_GetAnimationRibbonBufferSizes(baked->animation_manager, sizeof(tfx_ribbon_vertex_t));
	//ribbon_playback writes the ribbons and ribbon_compute and the ribbon draw read them
	SDL_GPUBufferUsageFlags ribbon_usage = SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_WRITE | SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_READ | SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ;
	if (!tfx_impl_sdlgpu__ensure_buffer(device, (tfxU32)baked->ribbon_offsets_bytes, SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_READ, &renderer->animation_ribbon_offsets)
		|| !tfx_impl_sdlgpu__upload_buffer(renderer, baked->ribbon_offsets, (tfxU32)baked->ribbon_offsets_bytes, renderer->animation_ribbon_offsets.buffer)
		|| !tfx_impl_sdlgpu__ensure_buffer(device, baked->ribbon_count * tfxRIBBON_SIZE, ribbon_usage, &renderer->baked_ribbons.instances)
		|| !tfx_impl_sdlgpu__ensure_ribbon_geometry(device, &sizes, &renderer->baked_ribbons)) {
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
static bool tfx_impl_sdlgpu__upload_baked(tfx_impl_sdlgpu_t *renderer, const tfx_impl_baked_frame_t *baked) {
	renderer->playback_sprite_count = 0;
	renderer->playback_ribbon_count = 0;
	renderer->baked_ribbons.animation_manager = NULL;
	if (!renderer->baked_sprites.buffer || (baked->sprite_count == 0 && baked->ribbon_count == 0)) {
		return true;
	}
	SDL_GPUDevice *device = renderer->device;
	if (!tfx_impl_sdlgpu__ensure_buffer(device, (tfxU32)baked->animation_instances_bytes, SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_READ, &renderer->animation_instances)
		|| !tfx_impl_sdlgpu__upload_buffer(renderer, baked->animation_instances, (tfxU32)baked->animation_instances_bytes, renderer->animation_instances.buffer)) {
		return false;
	}
	if (baked->sprite_count > 0) {
		if (!tfx_impl_sdlgpu__ensure_buffer(device, (tfxU32)baked->offsets_bytes, SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_READ, &renderer->animation_offsets)
			|| !tfx_impl_sdlgpu__upload_buffer(renderer, baked->offsets, (tfxU32)baked->offsets_bytes, renderer->animation_offsets.buffer)
			|| !tfx_impl_sdlgpu__ensure_buffer(device, baked->sprite_count * (tfxU32)sizeof(tfx_instance_t), SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_WRITE | SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ, &renderer->playback_instances)) {
			return false;
		}
		memset(&renderer->sprite_data_push, 0, sizeof(tfx_sprite_data_push_t));
		renderer->sprite_data_push.animation_instances_total = baked->animation_instance_count;
		renderer->sprite_data_push.billboards_total = baked->sprite_count;
		renderer->sprite_data_push.animated_shapes = baked->animated_shapes ? 1 : 0;
		renderer->playback_sprite_count = baked->sprite_count;
	}
	return tfx_impl_sdlgpu__upload_baked_ribbons(renderer, baked);
}

static bool tfx_impl_sdlgpu__upload_ribbons(tfx_impl_sdlgpu_t *renderer, const tfx_impl_ribbon_frame_t *ribbons) {
	tfx_impl_sdlgpu_ribbon_set_t *set = &renderer->live_ribbons;
	set->stage = ribbons->stage;
	if (!ribbons->stage) {
		return true;
	}
	SDL_GPUDevice *device = renderer->device;
	const tfx_ribbon_batch_sizes_t *sizes = &ribbons->sizes;
	SDL_GPUBufferUsageFlags ribbon_read = SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_READ | SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ;
	if (!tfx_impl_sdlgpu__ensure_buffer(device, (tfxU32)sizes->segment_buffer_size_in_bytes, SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_READ, &set->segments)
		|| !tfx_impl_sdlgpu__ensure_buffer(device, (tfxU32)sizes->ribbon_buffer_size_in_bytes, ribbon_read, &set->instances)
		|| !tfx_impl_sdlgpu__ensure_buffer(device, (tfxU32)sizes->emitter_buffer_size_in_bytes, ribbon_read, &set->emitters)
		|| !tfx_impl_sdlgpu__ensure_ribbon_geometry(device, sizes, set)
		|| !tfx_impl_sdlgpu__upload_buffer(renderer, ribbons->segments, (tfxU32)sizes->segment_buffer_size_in_bytes, set->segments.buffer)
		|| !tfx_impl_sdlgpu__upload_buffer(renderer, ribbons->ribbons, (tfxU32)sizes->ribbon_buffer_size_in_bytes, set->instances.buffer)
		|| !tfx_impl_sdlgpu__upload_buffer(renderer, ribbons->emitters, (tfxU32)sizes->emitter_buffer_size_in_bytes, set->emitters.buffer)) {
		set->stage = NULL;
		return false;
	}
	set->stage_offsets = ribbons->stage_offsets;
	return true;
}

//The old contents are copied across, because the previous buffer is still what captured_index points into
static bool tfx_impl_sdlgpu__grow_instance_buffers(tfx_impl_sdlgpu_t *renderer, tfxU32 needed) {
	tfxU32 capacity = renderer->instance_capacity ? renderer->instance_capacity : TFX_IMPL_SDLGPU_INITIAL_INSTANCE_CAPACITY;
	while (capacity < needed) {
		capacity += capacity / 2;
	}
	//The GPU particle kernel writes its instances after the CPU ones
	SDL_GPUBufferUsageFlags usage = SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ | SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_WRITE;
	for (int index = 0; index != 2; ++index) {
		tfx_impl_sdlgpu_buffer_t grown = {};
		if (!tfx_impl_sdlgpu__create_buffer(renderer->device, capacity * (tfxU32)sizeof(tfx_instance_t), usage, &grown)) {
			return false;
		}
		tfx_impl_sdlgpu_buffer_t *old = &renderer->instance_buffers[index];
		if (old->buffer) {
			tfx_impl_sdlgpu_copy_t *copy = &renderer->copies[renderer->copy_count++];
			copy->source = old->buffer;
			copy->destination = grown.buffer;
			copy->size = old->capacity;
		}
		*old = grown;
	}
	renderer->instance_capacity = capacity;
	return true;
}

static bool tfx_impl_sdlgpu__upload_gpu_particles(tfx_impl_sdlgpu_t *renderer, const tfx_impl_frame_t *frame) {
	const tfx_gpu_particle_frame_t *gpu_frame = &frame->gpu_particles;
	SDL_GPUDevice *device = renderer->device;
	//The pool holds the particles between frames, it's sized once and has to start zeroed, which SDL buffers don't
	if (!renderer->particle_pool.buffer) {
		tfxU32 pool_bytes = (tfxU32)frame->gpu_pool_sizes.pool_bytes;
		void *zeroes = calloc(1, pool_bytes);
		if (!zeroes || !tfx_impl_sdlgpu__create_buffer(device, pool_bytes, SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_WRITE, &renderer->particle_pool)) {
			free(zeroes);
			return false;
		}
		renderer->upload_scratch = zeroes;
		if (!tfx_impl_sdlgpu__upload_buffer(renderer, zeroes, pool_bytes, renderer->particle_pool.buffer)) {
			return false;
		}
	}
	if (!tfx_impl_sdlgpu__ensure_buffer(device, (tfxU32)gpu_frame->upload_bytes, SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_READ, &renderer->particle_upload)) {
		return false;
	}
	tfx_impl_sdlgpu_upload_t *upload = tfx_impl_sdlgpu__add_upload(renderer, NULL, (tfxU32)gpu_frame->upload_bytes);
	if (!upload) {
		return false;
	}
	upload->buffer = renderer->particle_upload.buffer;
	renderer->particle_frame = gpu_frame;

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

static bool tfx_impl_sdlgpu__upload_global_buffer(tfx_impl_sdlgpu_t *renderer, const void *data, tfxU32 size_in_bytes, SDL_GPUBufferUsageFlags usage, tfx_impl_sdlgpu_buffer_t *buffer) {
	return tfx_impl_sdlgpu__ensure_buffer(renderer->device, size_in_bytes, usage, buffer)
		&& tfx_impl_sdlgpu__upload_buffer(renderer, data, size_in_bytes, buffer->buffer);
}

//Global layers never move, so only the layers whose version changed go up, unless the texture array has to grow
static bool tfx_impl_sdlgpu__upload_color_ramps(tfx_impl_sdlgpu_t *renderer) {
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
		capacity = capacity < TFX_IMPL_SDLGPU_MAX_COLOR_RAMP_LAYERS ? capacity : TFX_IMPL_SDLGPU_MAX_COLOR_RAMP_LAYERS;
		tfxU32 *layer_versions = (tfxU32 *)realloc(renderer->color_ramp_layer_versions, sizeof(tfxU32) * capacity);
		if (!layer_versions) {
			return false;
		}
		renderer->color_ramp_layer_versions = layer_versions;
		if (renderer->color_ramps) {
			SDL_ReleaseGPUTexture(renderer->device, renderer->color_ramps);
		}
		renderer->color_ramp_layer_capacity = 0;
		renderer->color_ramps = tfx_impl_sdlgpu__create_color_ramp_texture(renderer->device, width, height, bytes_per_pixel, capacity);
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
		if (bitmap && !tfx_impl_sdlgpu__upload_texture(renderer, tfx_GetBitmapData(bitmap), (tfxU32)tfx_GetBitmapSize(bitmap), renderer->color_ramps, layer, 0, width, height)) {
			return false;
		}
	}
	return true;
}

//Only on update frames: the instances drawn until then were written with the layout the last upload holds
static bool tfx_impl_sdlgpu__upload_library_data(tfx_impl_sdlgpu_t *renderer) {
	tfxU32 properties_version = tfx_GetGlobalParticlePropertiesVersion();
	if (!renderer->emitter_properties.buffer || properties_version != renderer->properties_version) {
		if (!tfx_impl_sdlgpu__upload_global_buffer(renderer, tfx_GetGlobalParticlePropertiesBuffer(), (tfxU32)tfx_GetGlobalParticlePropertiesBufferSizeInBytes(),
			SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ, &renderer->emitter_properties)) {
			return false;
		}
		renderer->properties_version = properties_version;
	}
	tfxU32 shapes_version = tfx_GetGlobalGPUShapesVersion();
	if (!renderer->image_data.buffer || shapes_version != renderer->shapes_version) {
		tfx_gpu_shapes shapes = tfx_GetGlobalGPUShapes();
		if (!tfx_impl_sdlgpu__upload_global_buffer(renderer, tfx_GetGPUShapesArray(shapes), (tfxU32)sizeof(tfx_gpu_image_data_t) * tfx_GetGPUShapesCount(shapes),
			SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ, &renderer->image_data)) {
			return false;
		}
		renderer->shapes_version = shapes_version;
	}
	//Shared by every library already, read by the ribbon shaders and the GPU kernel
	tfxU32 graph_lookups_version = tfx_GetGPUGraphLookupsVersion();
	if (!renderer->graph_lookups.buffer || graph_lookups_version != renderer->graph_lookups_version) {
		if (!tfx_impl_sdlgpu__upload_global_buffer(renderer, tfx_GetGPUGraphLookupsBuffer(), tfx_GetGPUGraphLookupsBufferSizeInBytes(),
			SDL_GPU_BUFFERUSAGE_COMPUTE_STORAGE_READ | SDL_GPU_BUFFERUSAGE_GRAPHICS_STORAGE_READ, &renderer->graph_lookups)) {
			return false;
		}
		renderer->graph_lookups_version = graph_lookups_version;
	}
	tfxU32 color_ramps_version = tfx_GetGlobalColorRampsVersion();
	if (!renderer->color_ramps || color_ramps_version != renderer->color_ramps_version) {
		if (!tfx_impl_sdlgpu__upload_color_ramps(renderer)) {
			return false;
		}
		renderer->color_ramps_version = color_ramps_version;
	}
	return true;
}

static bool tfx_impl_sdlgpu__queue_frame(tfx_impl_sdlgpu_t *renderer, const tfx_impl_frame_t *frame) {
	memcpy(renderer->layer_flags, frame->layer_flags, sizeof(renderer->layer_flags));
	renderer->library_layer_count = frame->library_layer_count;
	if (!tfx_impl_sdlgpu__upload_baked(renderer, &frame->baked)) {
		return false;
	}
	if (!frame->stage_updated) {
		return true;
	}
	if (!tfx_impl_sdlgpu__upload_library_data(renderer)) {
		return false;
	}
	bool gpu_particles = frame->has_gpu_particles && renderer->particle_kernel;
	tfxU32 needed = gpu_particles ? frame->gpu_particles.gpu_instance_start + frame->gpu_particles.gpu_instance_count : frame->instance_count;
	if ((needed > renderer->instance_capacity || renderer->instance_capacity == 0) && !tfx_impl_sdlgpu__grow_instance_buffers(renderer, needed)) {
		renderer->instance_count = 0;
		return false;
	}

	//The other buffer now holds what was drawn last, which is what captured_index points into for CPU and GPU instances alike
	renderer->current_instance_buffer ^= 1;
	if (!tfx_impl_sdlgpu__upload_buffer(renderer, frame->instances, frame->instance_count * (tfxU32)sizeof(tfx_instance_t), renderer->instance_buffers[renderer->current_instance_buffer].buffer)) {
		return false;
	}
	renderer->instance_count = needed;
	memcpy(renderer->layer_start, frame->layer_start, sizeof(renderer->layer_start));
	memcpy(renderer->layer_count, frame->layer_count, sizeof(renderer->layer_count));
	memset(renderer->gpu_layer_count, 0, sizeof(renderer->gpu_layer_count));
	if (gpu_particles) {
		memcpy(renderer->gpu_layer_start, frame->gpu_particles.gpu_layer_start, sizeof(renderer->gpu_layer_start));
		memcpy(renderer->gpu_layer_count, frame->gpu_particles.gpu_layer_count, sizeof(renderer->gpu_layer_count));
		if (!tfx_impl_sdlgpu__upload_gpu_particles(renderer, frame)) {
			return false;
		}
	}
	return tfx_impl_sdlgpu__upload_ribbons(renderer, &frame->ribbons);
}

bool tfx_impl_sdlgpu_UploadFrame(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer, const tfx_impl_frame_t *frame) {
	//Whatever was queued before a failure still has to go, the grown instance buffers' copies included
	bool queued = tfx_impl_sdlgpu__queue_frame(renderer, frame);
	return tfx_impl_sdlgpu__flush_uploads(renderer, command_buffer) && queued;
}

//Every fetch has to be dispatched, the next one carries on from the particle state this one leaves in the pool
static void tfx_impl_sdlgpu__simulate_gpu_particles(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer) {
	//Never cycled: the pool carries the particles over and the instance buffer already holds this frame's CPU instances
	SDL_GPUStorageBufferReadWriteBinding writes[2] = {};
	writes[0].buffer = renderer->particle_pool.buffer;
	writes[1].buffer = renderer->instance_buffers[renderer->current_instance_buffer].buffer;
	SDL_GPUComputePass *compute_pass = SDL_BeginGPUComputePass(command_buffer, NULL, 0, writes, 2);
	SDL_BindGPUComputePipeline(compute_pass, renderer->particle_kernel);
	SDL_GPUBuffer *reads[2] = { renderer->graph_lookups.buffer, renderer->particle_upload.buffer };
	SDL_BindGPUComputeStorageBuffers(compute_pass, 0, reads, 2);
	tfx_gpu_particle_push_t *push = &renderer->particle_push_constants;
	for (tfxU32 index = 0; index != renderer->particle_dispatch_count; ++index) {
		tfxU32 slot_count = renderer->particle_dispatches[index].slot_count;
		if (slot_count == 0) {
			continue;
		}
		push->dispatch_index = index;
		SDL_PushGPUComputeUniformData(command_buffer, 0, push, sizeof(tfx_gpu_particle_push_t));
		SDL_DispatchGPUCompute(compute_pass, (slot_count + 255) / 256, 1, 1);
	}
	SDL_EndGPUComputePass(compute_pass);
	renderer->particle_dispatch_count = 0;
}

//Compute passes that rewrite all of what they write cycle it, so they never wait on a frame still drawing from it
static void tfx_impl_sdlgpu__play_back(SDL_GPUCommandBuffer *command_buffer, SDL_GPUComputePipeline *pipeline, SDL_GPUBuffer **reads, SDL_GPUBuffer *write,
	const void *push, Uint32 push_size, tfxU32 count) {
	SDL_GPUStorageBufferReadWriteBinding write_binding = {};
	write_binding.buffer = write;
	write_binding.cycle = true;
	SDL_GPUComputePass *compute_pass = SDL_BeginGPUComputePass(command_buffer, NULL, 0, &write_binding, 1);
	SDL_BindGPUComputePipeline(compute_pass, pipeline);
	SDL_BindGPUComputeStorageBuffers(compute_pass, 0, reads, 4);
	SDL_PushGPUComputeUniformData(command_buffer, 0, push, push_size);
	SDL_DispatchGPUCompute(compute_pass, (count + 127) / 128, 1, 1);
	SDL_EndGPUComputePass(compute_pass);
}

static bool tfx_impl_sdlgpu__ribbon_set_active(const tfx_impl_sdlgpu_ribbon_set_t *set) {
	return set->stage || set->animation_manager;
}

static tfx_ribbon_dispatch_t tfx_impl_sdlgpu__create_ribbon_dispatch(const tfx_impl_sdlgpu_ribbon_set_t *set) {
	return set->stage ? tfx_CreateRibbonBatchDispatch(&set->stage_offsets) : tfx_CreateAnimationRibbonDispatch();
}

static bool tfx_impl_sdlgpu__next_ribbon_dispatch(const tfx_impl_sdlgpu_ribbon_set_t *set, tfx_ribbon_dispatch_t *dispatch) {
	return set->stage ? tfx_NextRibbonDispatch(set->stage, dispatch) : tfx_NextAnimationRibbonDispatch(set->animation_manager, dispatch);
}

static SDL_GPUBuffer *tfx_impl_sdlgpu__ribbon_graph_lookups(const tfx_impl_sdlgpu_t *renderer, const tfx_impl_sdlgpu_ribbon_set_t *set) {
	return set->stage ? renderer->graph_lookups.buffer : set->graph_lookups.buffer;
}

static void tfx_impl_sdlgpu__build_ribbon_geometry(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer, const tfx_impl_sdlgpu_ribbon_set_t *set, const tfx_impl_uniforms_t *uniforms) {
	SDL_GPUStorageBufferReadWriteBinding writes[2] = {};
	writes[0].buffer = set->vertices.buffer;
	writes[0].cycle = true;
	writes[1].buffer = set->indices.buffer;
	writes[1].cycle = true;
	SDL_GPUComputePass *compute_pass = SDL_BeginGPUComputePass(command_buffer, NULL, 0, writes, 2);
	SDL_BindGPUComputePipeline(compute_pass, renderer->ribbon_compute_pipeline);
	SDL_GPUBuffer *reads[4] = { set->emitters.buffer, tfx_impl_sdlgpu__ribbon_graph_lookups(renderer, set), set->instances.buffer, set->segments.buffer };
	SDL_BindGPUComputeStorageBuffers(compute_pass, 0, reads, 4);

	//One dispatch per bucket, each with the bucket's own globals
	tfx_ribbon_dispatch_t dispatch = tfx_impl_sdlgpu__create_ribbon_dispatch(set);
	while (tfx_impl_sdlgpu__next_ribbon_dispatch(set, &dispatch)) {
		tfx_ribbon_bucket_globals_t *globals = tfx_GetRibbonDispatchGlobals(&dispatch);
		if (set->stage) {
			globals->lerp = uniforms->timer_lerp;
		} else {
			//Baked ribbons are already interpolated, but they still face the camera
			globals->camera_position.x = uniforms->camera_position[0];
			globals->camera_position.y = uniforms->camera_position[1];
			globals->camera_position.z = uniforms->camera_position[2];
		}
		SDL_PushGPUComputeUniformData(command_buffer, 0, globals, sizeof(tfx_ribbon_bucket_globals_t));
		SDL_DispatchGPUCompute(compute_pass, (dispatch.total_segments + TFX_IMPL_SDLGPU_RIBBON_GROUP_SIZE - 1) / TFX_IMPL_SDLGPU_RIBBON_GROUP_SIZE, 1, 1);
	}
	SDL_EndGPUComputePass(compute_pass);
}

void tfx_impl_sdlgpu_Compute(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer, const tfx_impl_uniforms_t *uniforms) {
	if (renderer->particle_dispatch_count > 0) {
		tfx_impl_sdlgpu__simulate_gpu_particles(renderer, command_buffer);
	}
	if (renderer->playback_sprite_count > 0) {
		SDL_GPUBuffer *reads[4] = { renderer->animation_offsets.buffer, renderer->animation_instances.buffer, renderer->baked_sprites.buffer, renderer->baked_emitter_properties.buffer };
		tfx_impl_sdlgpu__play_back(command_buffer, renderer->sprite_playback_pipeline, reads, renderer->playback_instances.buffer,
			&renderer->sprite_data_push, sizeof(tfx_sprite_data_push_t), renderer->playback_sprite_count);
	}
	if (renderer->playback_ribbon_count > 0) {
		SDL_GPUBuffer *reads[4] = { renderer->animation_ribbon_offsets.buffer, renderer->animation_instances.buffer, renderer->baked_ribbon_data.buffer, renderer->baked_ribbon_properties.buffer };
		tfx_impl_sdlgpu__play_back(command_buffer, renderer->ribbon_playback_pipeline, reads, renderer->baked_ribbons.instances.buffer,
			&renderer->ribbon_data_push, sizeof(tfx_ribbon_data_push_t), renderer->playback_ribbon_count);
	}
	if (tfx_impl_sdlgpu__ribbon_set_active(&renderer->live_ribbons)) {
		tfx_impl_sdlgpu__build_ribbon_geometry(renderer, command_buffer, &renderer->live_ribbons, uniforms);
	}
	if (tfx_impl_sdlgpu__ribbon_set_active(&renderer->baked_ribbons)) {
		tfx_impl_sdlgpu__build_ribbon_geometry(renderer, command_buffer, &renderer->baked_ribbons, uniforms);
	}
}

//Live and baked ribbons each have their own shape frames and color ramps, the atlas is shared. Every bucket is walked, because the
//dispatch carries the offsets from one to the next, and only the ones whose layer is in the group are drawn
static void tfx_impl_sdlgpu__render_ribbons(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer, SDL_GPURenderPass *render_pass,
	const tfx_impl_sdlgpu_ribbon_set_t *set, SDL_GPUBuffer *image_data, SDL_GPUTexture *color_ramps, int writes_depth) {
	SDL_BindGPUGraphicsPipeline(render_pass, renderer->ribbon_pipelines[writes_depth][0]);
	SDL_GPUGraphicsPipeline *bound_pipeline = renderer->ribbon_pipelines[writes_depth][0];
	SDL_GPUBufferBinding index_binding = { set->indices.buffer, 0 };
	SDL_BindGPUIndexBuffer(render_pass, &index_binding, SDL_GPU_INDEXELEMENTSIZE_32BIT);
	SDL_GPUBuffer *vertex_reads[5] = { set->vertices.buffer, set->instances.buffer, set->emitters.buffer, tfx_impl_sdlgpu__ribbon_graph_lookups(renderer, set), image_data };
	SDL_BindGPUVertexStorageBuffers(render_pass, 0, vertex_reads, 5);
	SDL_GPUTextureSamplerBinding samplers[2] = { { renderer->shape_atlas, renderer->atlas_sampler }, { color_ramps, renderer->color_ramp_sampler } };
	SDL_BindGPUFragmentSamplers(render_pass, 0, samplers, 2);

	tfx_ribbon_dispatch_t dispatch = tfx_impl_sdlgpu__create_ribbon_dispatch(set);
	while (tfx_impl_sdlgpu__next_ribbon_dispatch(set, &dispatch)) {
		const tfxLayerFlags flags = renderer->layer_flags[dispatch.layer];
		if (tfx_impl_LayerWritesDepth(flags) != writes_depth) {
			continue;
		}
		SDL_GPUGraphicsPipeline *pipeline = renderer->ribbon_pipelines[writes_depth][tfx_impl_LayerIsLit(flags)];
		if (pipeline != bound_pipeline) {
			SDL_BindGPUGraphicsPipeline(render_pass, pipeline);
			bound_pipeline = pipeline;
		}
		SDL_PushGPUVertexUniformData(command_buffer, 1, tfx_GetRibbonDispatchGlobals(&dispatch), sizeof(tfx_ribbon_bucket_globals_t));
		SDL_DrawGPUIndexedPrimitives(render_pass, dispatch.index_count, 1, dispatch.index_offset, 0, 0);
	}
}

static void tfx_impl_sdlgpu__draw_instance_range(SDL_GPUCommandBuffer *command_buffer, SDL_GPURenderPass *render_pass, tfxU32 start, tfxU32 count) {
	if (count == 0) {
		return;
	}
	//The first instance goes in the draw constants: D3D12 leaves it out of SV_InstanceID, and Slang subtracts it on Vulkan
	tfx_impl_sdlgpu_draw_constants_t draw_constants = {};
	draw_constants.instance_offset = start;
	SDL_PushGPUVertexUniformData(command_buffer, 1, &draw_constants, sizeof(draw_constants));
	SDL_DrawGPUPrimitives(render_pass, 6, count, 0, 0);
}

//Each layer in the group as its CPU instances then its GPU instances
static void tfx_impl_sdlgpu__render_billboards(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer, SDL_GPURenderPass *render_pass, int writes_depth) {
	SDL_BindGPUGraphicsPipeline(render_pass, renderer->billboard_pipelines[writes_depth][0]);
	SDL_GPUGraphicsPipeline *bound_pipeline = renderer->billboard_pipelines[writes_depth][0];
	SDL_GPUBuffer *vertex_reads[4] = {
		renderer->instance_buffers[renderer->current_instance_buffer].buffer,
		renderer->instance_buffers[renderer->current_instance_buffer ^ 1].buffer,
		renderer->image_data.buffer,
		renderer->emitter_properties.buffer,
	};
	SDL_BindGPUVertexStorageBuffers(render_pass, 0, vertex_reads, 4);
	SDL_GPUTextureSamplerBinding samplers[2] = { { renderer->shape_atlas, renderer->atlas_sampler }, { renderer->color_ramps, renderer->color_ramp_sampler } };
	SDL_BindGPUFragmentSamplers(render_pass, 0, samplers, 2);
	for (tfxU32 layer = 0; layer != renderer->library_layer_count; ++layer) {
		const tfxLayerFlags flags = renderer->layer_flags[layer];
		if (tfx_impl_LayerWritesDepth(flags) != writes_depth || renderer->layer_count[layer] + renderer->gpu_layer_count[layer] == 0) {
			continue;
		}
		SDL_GPUGraphicsPipeline *pipeline = renderer->billboard_pipelines[writes_depth][tfx_impl_LayerIsLit(flags)];
		if (pipeline != bound_pipeline) {
			SDL_BindGPUGraphicsPipeline(render_pass, pipeline);
			bound_pipeline = pipeline;
		}
		tfx_impl_sdlgpu__draw_instance_range(command_buffer, render_pass, renderer->layer_start[layer], renderer->layer_count[layer]);
		tfx_impl_sdlgpu__draw_instance_range(command_buffer, render_pass, renderer->gpu_layer_start[layer], renderer->gpu_layer_count[layer]);
	}
}

void tfx_impl_sdlgpu_Render(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer, SDL_GPURenderPass *render_pass, const tfx_impl_uniforms_t *uniforms) {
	if (!renderer->shape_atlas || (renderer->instance_count == 0 && renderer->playback_sprite_count == 0
		&& !tfx_impl_sdlgpu__ribbon_set_active(&renderer->live_ribbons) && !tfx_impl_sdlgpu__ribbon_set_active(&renderer->baked_ribbons))) {
		return;
	}
	//Every vertex shader takes the frame's uniforms at slot 0 and its draw's own at slot 1, and the lit fragments take the frame's at fragment slot 0
	SDL_PushGPUVertexUniformData(command_buffer, 0, uniforms, sizeof(tfx_impl_uniforms_t));
	SDL_PushGPUFragmentUniformData(command_buffer, 0, uniforms, sizeof(tfx_impl_uniforms_t));

	//The live draws need the library data, which goes up with the first update
	const bool live_data = renderer->color_ramps && renderer->image_data.buffer;
	//The layers that write depth draw first with blending off, so they hide the blended layers behind them
	for (int group = 0; group != 2; ++group) {
		const int writes_depth = group == 0;
		if (live_data && renderer->instance_count > 0) {
			tfx_impl_sdlgpu__render_billboards(renderer, command_buffer, render_pass, writes_depth);
		}
		if (live_data && tfx_impl_sdlgpu__ribbon_set_active(&renderer->live_ribbons)) {
			tfx_impl_sdlgpu__render_ribbons(renderer, command_buffer, render_pass, &renderer->live_ribbons, renderer->image_data.buffer, renderer->color_ramps, writes_depth);
		}
		if (!renderer->baked_color_ramps) {
			continue;
		}
		//Baked billboards come as one range across their layers, so they all draw blended and unlit
		if (!writes_depth && renderer->playback_sprite_count > 0) {
			//Baked sprites come with their own shape frames, emitter properties and color ramps, the atlas is shared
			SDL_BindGPUGraphicsPipeline(render_pass, renderer->baked_billboard_pipeline);
			SDL_GPUBuffer *vertex_reads[3] = { renderer->playback_instances.buffer, renderer->baked_image_data.buffer, renderer->baked_emitter_properties.buffer };
			SDL_BindGPUVertexStorageBuffers(render_pass, 0, vertex_reads, 3);
			SDL_GPUTextureSamplerBinding samplers[2] = { { renderer->shape_atlas, renderer->atlas_sampler }, { renderer->baked_color_ramps, renderer->color_ramp_sampler } };
			SDL_BindGPUFragmentSamplers(render_pass, 0, samplers, 2);
			tfx_impl_sdlgpu__draw_instance_range(command_buffer, render_pass, 0, renderer->playback_sprite_count);
		}
		if (tfx_impl_sdlgpu__ribbon_set_active(&renderer->baked_ribbons)) {
			tfx_impl_sdlgpu__render_ribbons(renderer, command_buffer, render_pass, &renderer->baked_ribbons, renderer->baked_image_data.buffer, renderer->baked_color_ramps, writes_depth);
		}
	}
}

static void tfx_impl_sdlgpu__release_ribbon_set(SDL_GPUDevice *device, tfx_impl_sdlgpu_ribbon_set_t *set) {
	tfx_impl_sdlgpu__release_buffer(device, &set->segments);
	tfx_impl_sdlgpu__release_buffer(device, &set->instances);
	tfx_impl_sdlgpu__release_buffer(device, &set->emitters);
	tfx_impl_sdlgpu__release_buffer(device, &set->vertices);
	tfx_impl_sdlgpu__release_buffer(device, &set->indices);
	tfx_impl_sdlgpu__release_buffer(device, &set->graph_lookups);
}

void tfx_impl_sdlgpu_Shutdown(tfx_impl_sdlgpu_t *renderer) {
	SDL_GPUDevice *device = renderer->device;
	if (!device) {
		return;
	}
	for (int index = 0; index != 2; ++index) {
		tfx_impl_sdlgpu__release_buffer(device, &renderer->instance_buffers[index]);
	}
	tfx_impl_sdlgpu__release_buffer(device, &renderer->image_data);
	tfx_impl_sdlgpu__release_buffer(device, &renderer->emitter_properties);
	tfx_impl_sdlgpu__release_buffer(device, &renderer->graph_lookups);
	tfx_impl_sdlgpu__release_ribbon_set(device, &renderer->live_ribbons);
	tfx_impl_sdlgpu__release_ribbon_set(device, &renderer->baked_ribbons);
	tfx_impl_sdlgpu__release_buffer(device, &renderer->particle_pool);
	tfx_impl_sdlgpu__release_buffer(device, &renderer->particle_upload);
	tfx_impl_sdlgpu__release_buffer(device, &renderer->baked_image_data);
	tfx_impl_sdlgpu__release_buffer(device, &renderer->baked_emitter_properties);
	tfx_impl_sdlgpu__release_buffer(device, &renderer->baked_sprites);
	tfx_impl_sdlgpu__release_buffer(device, &renderer->animation_offsets);
	tfx_impl_sdlgpu__release_buffer(device, &renderer->animation_instances);
	tfx_impl_sdlgpu__release_buffer(device, &renderer->playback_instances);
	tfx_impl_sdlgpu__release_buffer(device, &renderer->baked_ribbon_data);
	tfx_impl_sdlgpu__release_buffer(device, &renderer->baked_ribbon_properties);
	tfx_impl_sdlgpu__release_buffer(device, &renderer->animation_ribbon_offsets);
	SDL_GPUTexture *textures[3] = { renderer->shape_atlas, renderer->color_ramps, renderer->baked_color_ramps };
	for (int index = 0; index != 3; ++index) {
		if (textures[index]) {
			SDL_ReleaseGPUTexture(device, textures[index]);
		}
	}
	if (renderer->transfer_buffer) SDL_ReleaseGPUTransferBuffer(device, renderer->transfer_buffer);
	if (renderer->atlas_sampler) SDL_ReleaseGPUSampler(device, renderer->atlas_sampler);
	if (renderer->color_ramp_sampler) SDL_ReleaseGPUSampler(device, renderer->color_ramp_sampler);
	for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
		for (int lit = 0; lit != 2; ++lit) {
			if (renderer->billboard_pipelines[writes_depth][lit]) SDL_ReleaseGPUGraphicsPipeline(device, renderer->billboard_pipelines[writes_depth][lit]);
			if (renderer->ribbon_pipelines[writes_depth][lit]) SDL_ReleaseGPUGraphicsPipeline(device, renderer->ribbon_pipelines[writes_depth][lit]);
		}
	}
	if (renderer->baked_billboard_pipeline) SDL_ReleaseGPUGraphicsPipeline(device, renderer->baked_billboard_pipeline);
	if (renderer->particle_kernel) SDL_ReleaseGPUComputePipeline(device, renderer->particle_kernel);
	if (renderer->ribbon_compute_pipeline) SDL_ReleaseGPUComputePipeline(device, renderer->ribbon_compute_pipeline);
	if (renderer->sprite_playback_pipeline) SDL_ReleaseGPUComputePipeline(device, renderer->sprite_playback_pipeline);
	if (renderer->ribbon_playback_pipeline) SDL_ReleaseGPUComputePipeline(device, renderer->ribbon_playback_pipeline);
	free(renderer->color_ramp_layer_versions);
	free(renderer->uploads);
	free(renderer->upload_scratch);
	memset(renderer, 0, sizeof(tfx_impl_sdlgpu_t));
}
