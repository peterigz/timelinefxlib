#ifndef TFX_IMPL_SDLGPU_HEADER
#define TFX_IMPL_SDLGPU_HEADER

/*
	TimelineFX SDL3 GPU (SDL_gpu) backend, for its Vulkan and Direct3D 12 drivers. Needs SDL 3.4 or later and the shaders slangc
	builds from ../Shaders/timelinefx.slang with SDL_gpu's -D line, see README.md.

	Every call records into a command buffer you acquire, so uploads, compute and drawing go in the order you call them:
	tfx_impl_sdlgpu_UploadFrame, then tfx_impl_sdlgpu_Compute, then tfx_impl_sdlgpu_Render inside your render pass.

	Table of contents, search for the keyword to jump to the section:
	[Sdlgpu_types]		The renderer struct that owns every SDL_gpu object
	[Sdlgpu_resources]	Buffer, upload and shader loading helpers (in the .cpp)
	[Sdlgpu_functions]	Init, upload, compute, render and shutdown
*/

#include <SDL3/SDL_gpu.h>
#include "tfx_impl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

//--- Sdlgpu_types

typedef struct tfx_impl_sdlgpu_buffer_s {
	SDL_GPUBuffer *buffer;
	tfxU32 capacity;				//In bytes
} tfx_impl_sdlgpu_buffer_t;

//One copy from the frame's transfer buffer into a buffer, or into one mip of one texture layer
typedef struct tfx_impl_sdlgpu_upload_s {
	const void *data;				//Null for the GPU particle upload, which tfx_CopyGPUParticleUpload writes in place
	tfxU32 size;
	tfxU32 transfer_offset;
	SDL_GPUBuffer *buffer;
	tfxU32 buffer_offset;
	SDL_GPUTexture *texture;
	tfxU32 layer;
	tfxU32 mip_level;
	tfxU32 width;
	tfxU32 height;
} tfx_impl_sdlgpu_upload_t;

//Grown instance buffers keep their contents, so the old one is copied across and released once the copy is recorded
typedef struct tfx_impl_sdlgpu_copy_s {
	SDL_GPUBuffer *source;
	SDL_GPUBuffer *destination;
	tfxU32 size;
} tfx_impl_sdlgpu_copy_t;

//What ribbon_compute reads and writes for one source of ribbons: a stage's live ribbons or an animation manager's baked ones
typedef struct tfx_impl_sdlgpu_ribbon_set_s {
	tfx_impl_sdlgpu_buffer_t segments;
	tfx_impl_sdlgpu_buffer_t instances;		//tfx_ribbon_t, uploaded for live ribbons and written by ribbon_playback for baked ones
	tfx_impl_sdlgpu_buffer_t emitters;
	tfx_impl_sdlgpu_buffer_t vertices;		//Built by ribbon_compute every frame
	tfx_impl_sdlgpu_buffer_t indices;		//Written by ribbon_compute and bound as the index buffer
	tfx_impl_sdlgpu_buffer_t graph_lookups;	//Baked only: tfx_GetAnimationRibbonGraphLookupsBufferPointer. Live ribbons read the renderer's shared table
	//The buckets are walked every frame because the geometry is rebuilt with each frame's lerp. One of these is set while there's something to draw
	tfx_stage stage;
	tfx_ribbon_batch_offsets_t stage_offsets;
	tfx_animation_manager animation_manager;
} tfx_impl_sdlgpu_ribbon_set_t;

typedef struct tfx_impl_sdlgpu_s {
	SDL_GPUDevice *device;

	SDL_GPUGraphicsPipeline *billboard_pipeline;
	SDL_GPUGraphicsPipeline *baked_billboard_pipeline;
	SDL_GPUGraphicsPipeline *ribbon_pipeline;
	SDL_GPUComputePipeline *particle_kernel;		//simulate_particles from timelinefx.slang
	SDL_GPUComputePipeline *ribbon_compute_pipeline;
	SDL_GPUComputePipeline *sprite_playback_pipeline;
	SDL_GPUComputePipeline *ribbon_playback_pipeline;
	SDL_GPUSampler *atlas_sampler;
	SDL_GPUSampler *color_ramp_sampler;			//The ramps are only loaded, but a combined image sampler still needs one

	//Every upload of a frame goes through one transfer buffer and one copy pass
	SDL_GPUTransferBuffer *transfer_buffer;
	tfxU32 transfer_capacity;
	tfx_impl_sdlgpu_upload_t *uploads;
	tfxU32 upload_count;
	tfxU32 upload_capacity;
	tfx_impl_sdlgpu_copy_t copies[2];
	tfxU32 copy_count;
	void *upload_scratch;							//Freed once the uploads that read it are copied, the zeroed pool for one
	const tfx_gpu_particle_frame_t *particle_frame;	//For the upload with no data

	SDL_GPUTexture *shape_atlas;
	tfxU32 shape_atlas_layer_size;
	tfxU32 shape_atlas_layer_count;
	tfxU32 shape_atlas_mip_count;
	//Every loaded library's data, from the global getters, so one set of buffers draws effects from any library
	SDL_GPUTexture *color_ramps;
	tfxU32 color_ramp_layer_capacity;
	tfxU32 *color_ramp_layer_versions;			//The global layer versions last uploaded, one per layer of color_ramps
	tfx_impl_sdlgpu_buffer_t image_data;		//tfx_gpu_image_data_t per shape frame
	tfx_impl_sdlgpu_buffer_t emitter_properties;	//tfx_gpu_particle_properties_t per emitter
	tfx_impl_sdlgpu_buffer_t graph_lookups;		//tfx_GetGPUGraphLookupsBuffer
	tfxU32 properties_version;
	tfxU32 shapes_version;
	tfxU32 color_ramps_version;
	tfxU32 graph_lookups_version;

	//Two instance buffers that swap on every update: the latest instances and the ones captured_index points into
	tfx_impl_sdlgpu_buffer_t instance_buffers[2];
	tfxU32 current_instance_buffer;
	tfxU32 instance_capacity;
	tfxU32 instance_count;			//CPU and GPU instances
	//Each layer is drawn as two ranges, its CPU instances then its GPU instances
	tfxU32 layer_start[tfxLAYERS];
	tfxU32 layer_count[tfxLAYERS];
	tfxU32 gpu_layer_start[tfxLAYERS];
	tfxU32 gpu_layer_count[tfxLAYERS];

	tfx_impl_sdlgpu_buffer_t particle_pool;		//The hot, cold and noise sections in one buffer
	tfx_impl_sdlgpu_buffer_t particle_upload;	//The frame's tfx_CopyGPUParticleUpload
	tfx_gpu_particle_push_t particle_push_constants;
	//From the last upload and still valid until the next stage update, dispatched once by the next tfx_impl_sdlgpu_Compute
	const tfx_gpu_particle_dispatch_t *particle_dispatches;
	tfxU32 particle_dispatch_count;

	//Baked sprite data playback has its own color ramps, shape frames and emitter properties, from the animation manager
	SDL_GPUTexture *baked_color_ramps;
	tfx_impl_sdlgpu_buffer_t baked_image_data;
	tfx_impl_sdlgpu_buffer_t baked_emitter_properties;
	tfx_impl_sdlgpu_buffer_t baked_sprites;		//Every baked sprite of every effect
	tfx_impl_sdlgpu_buffer_t animation_offsets;	//Uploaded every frame
	tfx_impl_sdlgpu_buffer_t animation_instances;
	tfx_impl_sdlgpu_buffer_t playback_instances;	//The billboard instances the playback compute writes
	tfx_sprite_data_push_t sprite_data_push;
	tfxU32 playback_sprite_count;
	tfx_impl_sdlgpu_buffer_t baked_ribbon_data;	//tfx_ribbon_instance_data_t, every baked ribbon of every effect
	tfx_impl_sdlgpu_buffer_t baked_ribbon_properties;
	tfx_impl_sdlgpu_buffer_t animation_ribbon_offsets;	//Uploaded every frame
	tfx_ribbon_data_push_t ribbon_data_push;
	tfxU32 playback_ribbon_count;

	tfx_impl_sdlgpu_ribbon_set_t live_ribbons;	//Uploaded when an update ran
	tfx_impl_sdlgpu_ribbon_set_t baked_ribbons;	//Segments and emitters uploaded once, the rest every frame
} tfx_impl_sdlgpu_t;

//--- Sdlgpu_functions

//Loads the shaders and creates the pipelines for the given targets. shader_folder holds each entry point of timelinefx.slang compiled to
//<entry point>.spv and, where D3D12 is wanted, <entry point>.dxil, with TFX_FIXED_BINDINGS and SDL_gpu's -D line
bool tfx_impl_sdlgpu_Init(tfx_impl_sdlgpu_t *renderer, SDL_GPUDevice *device, SDL_GPUTextureFormat color_format, SDL_GPUTextureFormat depth_format, const char *shader_folder);
//Creates the texture array that every library's atlas and the sprite data's atlas are uploaded into, with room for all their layers
bool tfx_impl_sdlgpu_CreateAtlas(tfx_impl_sdlgpu_t *renderer, tfxU32 layer_size, tfxU32 layer_count);
//Writes an atlas into its layers of the texture array with mips built on the CPU, so Vulkan and D3D12 sample the same ones
bool tfx_impl_sdlgpu_UploadAtlas(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer, const tfx_impl_atlas_t *atlas);
//Uploads an animation manager's baked sprite data, color ramps, shape frames and emitter properties, after tfx_BuildAnimationManagerGPUShapeData
bool tfx_impl_sdlgpu_UploadSpriteData(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer, tfx_animation_manager animation_manager, tfx_gpu_shapes shapes);
//Uploads the baked playback every frame, and on update frames the instances, the ribbon data and any library data tfx_UpdateGlobalGPUData changed
bool tfx_impl_sdlgpu_UploadFrame(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer, const tfx_impl_frame_t *frame);
//Runs the GPU particle kernel, the baked playback and builds the ribbon geometry, every frame before tfx_impl_sdlgpu_Render, outside any pass
void tfx_impl_sdlgpu_Compute(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer, const tfx_impl_uniforms_t *uniforms);
//Draws inside a render pass with the color and depth formats given to tfx_impl_sdlgpu_Init, every frame, interpolating with uniforms->timer_lerp
void tfx_impl_sdlgpu_Render(tfx_impl_sdlgpu_t *renderer, SDL_GPUCommandBuffer *command_buffer, SDL_GPURenderPass *render_pass, const tfx_impl_uniforms_t *uniforms);
//Call after SDL_WaitForGPUIdle
void tfx_impl_sdlgpu_Shutdown(tfx_impl_sdlgpu_t *renderer);

#ifdef __cplusplus
}
#endif

#endif
