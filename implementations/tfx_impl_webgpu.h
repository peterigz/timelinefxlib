#ifndef TFX_IMPL_WEBGPU_HEADER
#define TFX_IMPL_WEBGPU_HEADER

/*
	TimelineFX WebGPU backend, written against the standard webgpu.h (webgpu-headers) as wgpu-native v29 ships it. Needs the
	shaders slangc builds from ../Shaders/timelinefx.slang as WGSL with -D TFX_FIXED_BINDINGS, see README.md.

	Uploads go straight through the queue with wgpuQueueWriteBuffer and wgpuQueueWriteTexture, so they take no encoder. Each
	frame call tfx_impl_webgpu_UploadFrame, then tfx_impl_webgpu_Compute with your command encoder, then tfx_impl_webgpu_Render
	inside your render pass on the same encoder, and submit it after all three.

	Table of contents, search for the keyword to jump to the section:
	[Webgpu_types]		The renderer struct that owns every WebGPU object
	[Webgpu_resources]	Buffer, bind group, upload and shader loading helpers (in the .cpp)
	[Webgpu_functions]	Init, upload, compute, render and shutdown
*/

#include <webgpu/webgpu.h>
#include "tfx_impl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

//--- Webgpu_types

typedef struct tfx_impl_webgpu_buffer_s {
	WGPUBuffer buffer;
	tfxU32 capacity;				//In bytes
} tfx_impl_webgpu_buffer_t;

//Bind groups can't be changed once made, so each one is kept until a buffer or texture it might hold is replaced
typedef struct tfx_impl_webgpu_bind_group_s {
	WGPUBindGroup bind_group;
	tfxU32 resource_generation;		//The renderer's resource_generation when it was made
} tfx_impl_webgpu_bind_group_t;

//What ribbon_compute reads and writes for one source of ribbons: a stage's live ribbons or an animation manager's baked ones
typedef struct tfx_impl_webgpu_ribbon_set_s {
	tfx_impl_webgpu_buffer_t segments;
	tfx_impl_webgpu_buffer_t instances;		//tfx_ribbon_t, uploaded for live ribbons and written by ribbon_playback for baked ones
	tfx_impl_webgpu_buffer_t emitters;
	tfx_impl_webgpu_buffer_t vertices;		//Built by ribbon_compute every frame
	tfx_impl_webgpu_buffer_t indices;		//Written by ribbon_compute and bound as the index buffer
	tfx_impl_webgpu_buffer_t graph_lookups;	//Baked only: tfx_GetAnimationRibbonGraphLookupsBufferPointer. Live ribbons read the renderer's shared table
	tfx_impl_webgpu_bind_group_t compute_bind_group;
	tfx_impl_webgpu_bind_group_t render_bind_group;
	tfxU32 first_globals_block;				//Its buckets' globals in the uniform buffer, written by Compute and drawn with by Render
	//The buckets are walked every frame because the geometry is rebuilt with each frame's lerp. One of these is set while there's something to draw
	tfx_stage stage;
	tfx_ribbon_batch_offsets_t stage_offsets;
	tfx_animation_manager animation_manager;
} tfx_impl_webgpu_ribbon_set_t;

typedef struct tfx_impl_webgpu_s {
	WGPUDevice device;
	WGPUQueue queue;

	WGPUBindGroupLayout billboard_layout;
	WGPUBindGroupLayout baked_billboard_layout;
	WGPUBindGroupLayout ribbon_layout;
	WGPUBindGroupLayout particle_kernel_layout;
	WGPUBindGroupLayout ribbon_compute_layout;
	WGPUBindGroupLayout playback_layout;			//Sprite and ribbon playback read and write the same slots
	//Each layer picks its pipelines from its tfx_GetLayerFlags, by [writes depth][lit]. Baked billboards are always blended and unlit
	WGPURenderPipeline billboard_pipelines[2][2];
	WGPURenderPipeline baked_billboard_pipeline;
	WGPURenderPipeline ribbon_pipelines[2][2];
	WGPUComputePipeline particle_kernel;			//simulate_particles from timelinefx.slang
	WGPUComputePipeline ribbon_compute_pipeline;
	WGPUComputePipeline sprite_playback_pipeline;
	WGPUComputePipeline ribbon_playback_pipeline;
	WGPUSampler atlas_sampler;
	bool init_scope_popped;							//Set by the error scope around the pipelines, see tfx_impl_webgpu_Init
	bool init_scope_failed;

	//Every uniform block of a frame at 256 byte steps, bound with dynamic offsets: Render's view and draw constants first, then Compute's
	tfx_impl_webgpu_buffer_t uniforms;
	unsigned char *uniform_staging;					//Filled on the CPU and written with one wgpuQueueWriteBuffer by Compute and one by Render
	tfxU32 uniform_block_capacity;
	tfxU32 resource_generation;						//Goes up whenever a buffer or texture is replaced, so the bind groups that might hold it are remade

	WGPUTexture shape_atlas;
	WGPUTextureView shape_atlas_view;
	tfxU32 shape_atlas_layer_size;
	tfxU32 shape_atlas_layer_count;
	tfxU32 shape_atlas_mip_count;
	//Every loaded library's data, from the global getters, so one set of buffers draws effects from any library
	WGPUTexture color_ramps;
	WGPUTextureView color_ramps_view;
	tfxU32 color_ramp_layer_capacity;
	tfxU32 *color_ramp_layer_versions;				//The global layer versions last uploaded, one per layer of color_ramps
	tfx_impl_webgpu_buffer_t image_data;			//tfx_gpu_image_data_t per shape frame
	tfx_impl_webgpu_buffer_t emitter_properties;	//tfx_gpu_particle_properties_t per emitter
	tfx_impl_webgpu_buffer_t graph_lookups;			//tfx_GetGPUGraphLookupsBuffer
	tfxU32 properties_version;
	tfxU32 shapes_version;
	tfxU32 color_ramps_version;
	tfxU32 graph_lookups_version;

	//Two instance buffers that swap on every update: the latest instances and the ones captured_index points into
	tfx_impl_webgpu_buffer_t instance_buffers[2];
	tfx_impl_webgpu_bind_group_t billboard_bind_groups[2];	//One for each instance buffer being the current one
	tfxU32 current_instance_buffer;
	tfxU32 instance_capacity;
	tfxU32 instance_count;							//CPU and GPU instances
	//Each layer is drawn as two ranges, its CPU instances then its GPU instances
	tfxU32 layer_start[tfxLAYERS];
	tfxU32 layer_count[tfxLAYERS];
	tfxU32 gpu_layer_start[tfxLAYERS];
	tfxU32 gpu_layer_count[tfxLAYERS];
	tfxLayerFlags layer_flags[tfxLAYERS];

	tfx_impl_webgpu_buffer_t particle_pool;			//The hot, cold and noise sections in one buffer
	tfx_impl_webgpu_buffer_t particle_upload;		//The frame's tfx_CopyGPUParticleUpload
	void *particle_upload_staging;					//tfx_CopyGPUParticleUpload writes here first, the queue takes a copy
	size_t particle_upload_staging_capacity;
	tfx_impl_webgpu_bind_group_t particle_kernel_bind_groups[2];
	tfx_gpu_particle_push_t particle_push_constants;
	//From the last upload and still valid until the next stage update, dispatched once by the next tfx_impl_webgpu_Compute
	const tfx_gpu_particle_dispatch_t *particle_dispatches;
	tfxU32 particle_dispatch_count;

	//Baked sprite data playback has its own color ramps, shape frames and emitter properties, from the animation manager
	WGPUTexture baked_color_ramps;
	WGPUTextureView baked_color_ramps_view;
	tfx_impl_webgpu_buffer_t baked_image_data;
	tfx_impl_webgpu_buffer_t baked_emitter_properties;
	tfx_impl_webgpu_buffer_t baked_sprites;			//Every baked sprite of every effect
	tfx_impl_webgpu_buffer_t animation_offsets;		//Uploaded every frame
	tfx_impl_webgpu_buffer_t animation_instances;
	tfx_impl_webgpu_buffer_t playback_instances;	//The billboard instances the playback compute writes
	tfx_impl_webgpu_bind_group_t sprite_playback_bind_group;
	tfx_impl_webgpu_bind_group_t baked_billboard_bind_group;
	tfx_sprite_data_push_t sprite_data_push;
	tfxU32 playback_sprite_count;
	tfx_impl_webgpu_buffer_t baked_ribbon_data;		//tfx_ribbon_instance_data_t, every baked ribbon of every effect
	tfx_impl_webgpu_buffer_t baked_ribbon_properties;
	tfx_impl_webgpu_buffer_t animation_ribbon_offsets;	//Uploaded every frame
	tfx_impl_webgpu_bind_group_t ribbon_playback_bind_group;
	tfx_ribbon_data_push_t ribbon_data_push;
	tfxU32 playback_ribbon_count;

	tfx_impl_webgpu_ribbon_set_t live_ribbons;		//Uploaded when an update ran
	tfx_impl_webgpu_ribbon_set_t baked_ribbons;		//Segments and emitters uploaded once, the rest every frame
} tfx_impl_webgpu_t;

//--- Webgpu_functions

//Loads the shaders and creates the pipelines for the color and depth formats you render to. shader_folder holds each entry point
//of timelinefx.slang compiled to <entry point>.wgsl with TFX_FIXED_BINDINGS
bool tfx_impl_webgpu_Init(tfx_impl_webgpu_t *renderer, WGPUDevice device, WGPUTextureFormat color_format, WGPUTextureFormat depth_format, const char *shader_folder);
//Creates the texture array that every library's atlas and the sprite data's atlas are uploaded into, with room for all their layers
bool tfx_impl_webgpu_CreateAtlas(tfx_impl_webgpu_t *renderer, tfxU32 layer_size, tfxU32 layer_count);
//Writes an atlas into its layers of the texture array, with mips built on the CPU because WebGPU can't generate them
bool tfx_impl_webgpu_UploadAtlas(tfx_impl_webgpu_t *renderer, const tfx_impl_atlas_t *atlas);
//Uploads an animation manager's baked sprite data, color ramps, shape frames and emitter properties, after tfx_BuildAnimationManagerGPUShapeData
bool tfx_impl_webgpu_UploadSpriteData(tfx_impl_webgpu_t *renderer, tfx_animation_manager animation_manager, tfx_gpu_shapes shapes);
//Uploads the baked playback every frame, and on update frames the instances, the ribbon data and any library data tfx_UpdateGlobalGPUData changed
bool tfx_impl_webgpu_UploadFrame(tfx_impl_webgpu_t *renderer, const tfx_impl_frame_t *frame);
//Runs the GPU particle kernel, the baked playback and builds the ribbon geometry in one compute pass, every frame before
//tfx_impl_webgpu_Render, outside any render pass
void tfx_impl_webgpu_Compute(tfx_impl_webgpu_t *renderer, WGPUCommandEncoder encoder, const tfx_impl_uniforms_t *uniforms);
//Draws inside a render pass with the color and depth formats given to tfx_impl_webgpu_Init, every frame, interpolating with uniforms->timer_lerp
void tfx_impl_webgpu_Render(tfx_impl_webgpu_t *renderer, WGPURenderPassEncoder render_pass, const tfx_impl_uniforms_t *uniforms);
void tfx_impl_webgpu_Shutdown(tfx_impl_webgpu_t *renderer);

#ifdef __cplusplus
}
#endif

#endif
