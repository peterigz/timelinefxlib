#ifndef TFX_IMPL_METAL_HEADER
#define TFX_IMPL_METAL_HEADER

/*
	TimelineFX Metal backend, Objective-C++ built without ARC (-fno-objc-arc), so the renderer is a plain struct that releases what it owns.
	Needs Metal 3 (macOS 13) for the GPU particle kernel's gpuAddress table, and the shaders slangc builds from ../Shaders/timelinefx.slang
	as Metal source: the render passes with -D TFX_FIXED_BINDINGS -D TFX_METAL_SLOTS and the kernel with no defines, see README.md.

	Each frame call tfx_impl_metal_BeginFrame once the GPU has finished the last frame that used its slot, then tfx_impl_metal_UploadFrame
	and tfx_impl_metal_Compute with your command buffer, then tfx_impl_metal_Render inside your render pass on the same command buffer.

	Table of contents, search for the keyword to jump to the section:
	[Metal_types]		The renderer struct that owns every Metal object
	[Metal_resources]	Buffers, the upload rings, textures and pipelines (in the .mm)
	[Metal_functions]	Init, upload, compute, render and shutdown
*/

#import <Metal/Metal.h>
#include "tfx_impl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

//--- Metal_types

#define TFX_IMPL_METAL_MAX_FRAMES_IN_FLIGHT 3

//A private buffer, filled by blits from the frame's upload ring or written by a compute pass
typedef struct tfx_impl_metal_buffer_s {
	id<MTLBuffer> buffer;
	size_t capacity;				//In bytes
} tfx_impl_metal_buffer_t;

//A place in shared memory the GPU reads this frame only: a slice of the frame slot's ring, or a big upload's buffer of its own
typedef struct tfx_impl_metal_upload_s {
	id<MTLBuffer> buffer;			//Owned by the frame slot
	size_t offset;
} tfx_impl_metal_upload_t;

//Each frame slot's shared upload memory. The CPU only writes it once the GPU has finished the slot's last frame
typedef struct tfx_impl_metal_frame_slot_s {
	id<MTLBuffer> ring;
	unsigned char *mapped;
	size_t capacity;
	size_t used;
	//Big uploads and outgrown rings, released when the slot comes round again because later passes of the frame may still bind them
	id<MTLBuffer> *retired;
	tfxU32 retired_count;
	tfxU32 retired_capacity;
} tfx_impl_metal_frame_slot_t;

//What ribbon_compute reads and writes for one source of ribbons: a stage's live ribbons or an animation manager's baked ones
typedef struct tfx_impl_metal_ribbon_set_s {
	tfx_impl_metal_buffer_t segments;
	tfx_impl_metal_buffer_t instances;		//tfx_ribbon_t, uploaded for live ribbons and written by ribbon_playback for baked ones
	tfx_impl_metal_buffer_t emitters;
	tfx_impl_metal_buffer_t vertices;		//Built by ribbon_compute every frame
	tfx_impl_metal_buffer_t indices;		//Written by ribbon_compute and drawn as the index buffer
	tfx_impl_metal_buffer_t graph_lookups;	//Baked only: tfx_GetAnimationRibbonGraphLookupsBufferPointer. Live ribbons read the renderer's shared table
	//The buckets are walked every frame because the geometry is rebuilt with each frame's lerp. One of these is set while there's something to draw
	tfx_stage stage;
	tfx_ribbon_batch_offsets_t stage_offsets;
	tfx_animation_manager animation_manager;
} tfx_impl_metal_ribbon_set_t;

typedef struct tfx_impl_metal_s {
	id<MTLDevice> device;
	tfxU32 frames_in_flight;
	tfxU32 frame_index;
	tfx_impl_metal_frame_slot_t frame_slots[TFX_IMPL_METAL_MAX_FRAMES_IN_FLIGHT];

	id<MTLRenderPipelineState> billboard_pipeline;
	id<MTLRenderPipelineState> baked_billboard_pipeline;
	id<MTLRenderPipelineState> ribbon_pipeline;
	id<MTLComputePipelineState> particle_kernel;			//simulate_particles from timelinefx.slang
	id<MTLComputePipelineState> ribbon_compute_pipeline;
	id<MTLComputePipelineState> sprite_playback_pipeline;
	id<MTLComputePipelineState> ribbon_playback_pipeline;
	id<MTLDepthStencilState> depth_state;					//Tested, not written
	id<MTLSamplerState> atlas_sampler;

	id<MTLTexture> shape_atlas;
	tfxU32 shape_atlas_layer_size;
	tfxU32 shape_atlas_layer_count;
	tfxU32 shape_atlas_mip_count;
	//Every loaded library's data, from the global getters, so one set of buffers draws effects from any library
	id<MTLTexture> color_ramps;
	tfxU32 color_ramp_layer_capacity;
	tfxU32 *color_ramp_layer_versions;				//The global layer versions last uploaded, one per layer of color_ramps
	tfx_impl_metal_buffer_t image_data;				//tfx_gpu_image_data_t per shape frame
	tfx_impl_metal_buffer_t emitter_properties;		//tfx_gpu_particle_properties_t per emitter
	tfx_impl_metal_buffer_t graph_lookups;			//tfx_GetGPUGraphLookupsBuffer
	tfxU32 properties_version;
	tfxU32 shapes_version;
	tfxU32 color_ramps_version;
	tfxU32 graph_lookups_version;

	//Two instance buffers that swap on every update: the latest instances and the ones captured_index points into
	tfx_impl_metal_buffer_t instance_buffers[2];
	tfxU32 current_instance_buffer;
	tfxU32 instance_capacity;
	tfxU32 instance_count;							//CPU and GPU instances
	//Each layer is drawn as two ranges, its CPU instances then its GPU instances
	tfxU32 layer_start[tfxLAYERS];
	tfxU32 layer_count[tfxLAYERS];
	tfxU32 gpu_layer_start[tfxLAYERS];
	tfxU32 gpu_layer_count[tfxLAYERS];

	//The hot, cold and noise sections in one buffer, reached through the kernel's gpuAddress table like everything else it reads
	tfx_impl_metal_buffer_t particle_pool;
	tfx_impl_metal_upload_t particle_upload;		//The frame's tfx_CopyGPUParticleUpload, read in place
	tfx_gpu_particle_push_t particle_push_constants;
	//From the last upload and still valid until the next stage update, dispatched once by the next tfx_impl_metal_Compute
	const tfx_gpu_particle_dispatch_t *particle_dispatches;
	tfxU32 particle_dispatch_count;

	//Baked sprite data playback has its own color ramps, shape frames and emitter properties, from the animation manager
	id<MTLTexture> baked_color_ramps;
	tfx_impl_metal_buffer_t baked_image_data;
	tfx_impl_metal_buffer_t baked_emitter_properties;
	tfx_impl_metal_buffer_t baked_sprites;			//Every baked sprite of every effect
	tfx_impl_metal_upload_t animation_offsets;		//Uploaded every frame and read in place
	tfx_impl_metal_upload_t animation_instances;
	tfx_impl_metal_buffer_t playback_instances;		//The billboard instances the playback compute writes
	tfx_sprite_data_push_t sprite_data_push;
	tfxU32 playback_sprite_count;
	tfx_impl_metal_buffer_t baked_ribbon_data;		//tfx_ribbon_instance_data_t, every baked ribbon of every effect
	tfx_impl_metal_buffer_t baked_ribbon_properties;
	tfx_impl_metal_upload_t animation_ribbon_offsets;	//Uploaded every frame and read in place
	tfx_ribbon_data_push_t ribbon_data_push;
	tfxU32 playback_ribbon_count;

	tfx_impl_metal_ribbon_set_t live_ribbons;		//Uploaded when an update ran
	tfx_impl_metal_ribbon_set_t baked_ribbons;		//Segments and emitters uploaded once, the rest every frame
} tfx_impl_metal_t;

//--- Metal_functions

//Compiles the shaders and creates the pipelines for the color and depth formats you render to. shader_folder holds each entry point of
//timelinefx.slang as <entry point>.metal, compiled at start up with newLibraryWithSource
bool tfx_impl_metal_Init(tfx_impl_metal_t *renderer, id<MTLDevice> device, MTLPixelFormat color_format, MTLPixelFormat depth_format, const char *shader_folder, tfxU32 frames_in_flight);
//Call first every frame, once the GPU has finished the last frame that used frame_index (0 to frames_in_flight - 1). Releases the slot's
//retired upload buffers and starts its ring again
void tfx_impl_metal_BeginFrame(tfx_impl_metal_t *renderer, tfxU32 frame_index);
//Creates the texture array that every library's atlas and the sprite data's atlas are uploaded into, with room for all their layers
bool tfx_impl_metal_CreateAtlas(tfx_impl_metal_t *renderer, tfxU32 layer_size, tfxU32 layer_count);
//Records an atlas's upload into its layers, with mips built on the CPU, through the frame slot's upload memory
bool tfx_impl_metal_UploadAtlas(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, const tfx_impl_atlas_t *atlas);
//Records the upload of an animation manager's baked sprite data, color ramps, shape frames and emitter properties, after tfx_BuildAnimationManagerGPUShapeData
bool tfx_impl_metal_UploadSpriteData(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, tfx_animation_manager animation_manager, tfx_gpu_shapes shapes);
//Records the baked playback's upload every frame, and on update frames the instances, the ribbon data and any library data tfx_UpdateGlobalGPUData changed
bool tfx_impl_metal_UploadFrame(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, const tfx_impl_frame_t *frame);
//Records the GPU particle kernel, the baked playback and the ribbon geometry in one compute pass, every frame before tfx_impl_metal_Render
void tfx_impl_metal_Compute(tfx_impl_metal_t *renderer, id<MTLCommandBuffer> command_buffer, const tfx_impl_uniforms_t *uniforms);
//Records the draws into your render pass, whose attachments have the formats given to tfx_impl_metal_Init, interpolating with uniforms->timer_lerp
void tfx_impl_metal_Render(tfx_impl_metal_t *renderer, id<MTLRenderCommandEncoder> render_encoder, const tfx_impl_uniforms_t *uniforms);
//Once the GPU is idle
void tfx_impl_metal_Shutdown(tfx_impl_metal_t *renderer);

#ifdef __cplusplus
}
#endif

#endif
