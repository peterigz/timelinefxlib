#ifndef TFX_IMPL_DX12_HEADER
#define TFX_IMPL_DX12_HEADER

/*
	TimelineFX Direct3D 12 backend, with the shaders' bindless model: every buffer and texture is a descriptor in one shader visible
	heap, indexed by the push constants, which are root constants here. Needs resource binding tier 3 and the shaders slangc builds
	from ../Shaders/timelinefx.slang as DXIL without TFX_FIXED_BINDINGS, see README.md.

	Table of contents, search for the keyword to jump to the section:
	[Dx12_types]		The renderer struct that owns every D3D12 object
	[Dx12_resources]	Descriptors, buffers, uploads, barriers and pipelines (in the .cpp)
	[Dx12_functions]	Init, upload, render and shutdown
*/

#include <d3d12.h>
#include "tfx_impl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

//--- Dx12_types

#define TFX_IMPL_DX12_MAX_FRAMES_IN_FLIGHT 3
#define TFX_IMPL_DX12_NO_DESCRIPTOR 0xFFFFFFFF

//A default heap buffer and the descriptor heap slots of its views, which are the indexes the shaders read it by
typedef struct tfx_impl_dx12_buffer_s {
	ID3D12Resource *resource;
	tfxU32 capacity;				//In bytes
	D3D12_RESOURCE_STATES state;
	tfxU32 shader_read_index;		//A structured view at the stride it's read at, or TFX_IMPL_DX12_NO_DESCRIPTOR
	tfxU32 shader_write_index;		//A structured unordered access view, for the buffers a compute shader writes
} tfx_impl_dx12_buffer_t;

//Resources and descriptor slots that a frame still in flight may use, released when that frame slot comes around again
typedef struct tfx_impl_dx12_release_list_s {
	IUnknown **objects;
	tfxU32 object_count;
	tfxU32 object_capacity;
	tfxU32 *descriptors;
	tfxU32 descriptor_count;
	tfxU32 descriptor_capacity;
} tfx_impl_dx12_release_list_t;

//A frame slot's persistently mapped upload buffer, filled from the start every time the slot comes around
typedef struct tfx_impl_dx12_upload_ring_s {
	ID3D12Resource *resource;
	unsigned char *mapped;
	size_t capacity;
	size_t used;
} tfx_impl_dx12_upload_ring_t;

//What ribbon_compute reads and writes for one source of ribbons: a stage's live ribbons or an animation manager's baked ones
typedef struct tfx_impl_dx12_ribbon_set_s {
	tfx_impl_dx12_buffer_t segments;
	tfx_impl_dx12_buffer_t instances;	//tfx_ribbon_t, uploaded for live ribbons and written by ribbon_playback for baked ones
	tfx_impl_dx12_buffer_t emitters;
	tfx_impl_dx12_buffer_t vertices;	//Built by ribbon_compute every frame, then the draw's vertex buffer
	tfx_impl_dx12_buffer_t indices;		//Built by ribbon_compute every frame, then the draw's index buffer
	tfx_impl_dx12_buffer_t graph_lookups;	//Baked only: tfx_GetAnimationRibbonGraphLookupsBufferPointer. Live ribbons read the renderer's shared table
	//The buckets are walked every frame because the geometry is rebuilt with each frame's lerp. One of these is set while there's something to draw
	tfx_stage stage;
	tfx_ribbon_batch_offsets_t stage_offsets;
	tfx_animation_manager animation_manager;
} tfx_impl_dx12_ribbon_set_t;

typedef struct tfx_impl_dx12_s {
	ID3D12Device *device;
	tfxU32 frames_in_flight;
	tfxU32 frame_index;
	tfx_impl_dx12_release_list_t release_lists[TFX_IMPL_DX12_MAX_FRAMES_IN_FLIGHT];
	tfx_impl_dx12_upload_ring_t upload_rings[TFX_IMPL_DX12_MAX_FRAMES_IN_FLIGHT];

	//Every view the shaders read, in one shader visible heap. Slot 0 is the atlas, because the uv lookup leaves the image index at 0
	ID3D12DescriptorHeap *descriptor_heap;
	ID3D12DescriptorHeap *sampler_heap;		//Slot 0, a clamped linear sampler with mips
	tfxU32 descriptor_size;
	tfxU32 descriptor_capacity;
	tfxU32 *free_descriptors;
	tfxU32 free_descriptor_count;

	//One root signature per push constant block, because together they're over the 64 DWORD limit
	ID3D12RootSignature *particle_root_signature;
	ID3D12RootSignature *billboard_root_signature;
	ID3D12RootSignature *baked_billboard_root_signature;
	ID3D12RootSignature *ribbon_root_signature;
	ID3D12RootSignature *sprite_playback_root_signature;
	ID3D12RootSignature *ribbon_playback_root_signature;
	ID3D12PipelineState *particle_kernel;		//simulate_particles from timelinefx.slang
	//Each layer picks its pipelines from its tfx_GetLayerFlags, by [writes depth][lit]. Baked billboards are always blended and unlit
	ID3D12PipelineState *billboard_pipelines[2][2];
	ID3D12PipelineState *baked_billboard_pipeline;
	ID3D12PipelineState *ribbon_pipelines[2][2];
	ID3D12PipelineState *ribbon_compute_pipeline;
	ID3D12PipelineState *sprite_playback_pipeline;
	ID3D12PipelineState *ribbon_playback_pipeline;

	tfx_impl_dx12_buffer_t view_uniforms;		//tfx_impl_uniforms_t behind a constant buffer view

	ID3D12Resource *shape_atlas;
	D3D12_RESOURCE_STATES shape_atlas_state;
	tfxU32 shape_atlas_layer_size;
	tfxU32 shape_atlas_layer_count;
	tfxU32 shape_atlas_mip_count;
	//Every loaded library's data, from the global getters, so one set of buffers draws effects from any library
	ID3D12Resource *color_ramps;
	D3D12_RESOURCE_STATES color_ramps_state;
	tfxU32 color_ramps_index;
	tfxU32 color_ramp_layer_capacity;
	tfxU32 *color_ramp_layer_versions;			//The global layer versions last uploaded, one per layer of color_ramps
	tfx_impl_dx12_buffer_t image_data;			//tfx_gpu_image_data_t per shape frame
	tfx_impl_dx12_buffer_t emitter_properties;	//tfx_gpu_particle_properties_t per emitter
	tfxU32 properties_version;
	tfxU32 shapes_version;
	tfxU32 color_ramps_version;
	tfxU32 graph_lookups_version;

	//Two instance buffers that swap on every update: the latest instances and the ones captured_index points into
	tfx_impl_dx12_buffer_t instance_buffers[2];
	tfxU32 current_instance_buffer;
	tfxU32 instance_capacity;
	tfxU32 instance_count;			//CPU and GPU instances
	//Each layer is drawn as two ranges, its CPU instances then its GPU instances
	tfxU32 layer_start[tfxLAYERS];
	tfxU32 layer_count[tfxLAYERS];
	tfxU32 gpu_layer_start[tfxLAYERS];
	tfxU32 gpu_layer_count[tfxLAYERS];
	tfxLayerFlags layer_flags[tfxLAYERS];

	//The pool's hot, cold and noise sections in one buffer, a view at each section's stride
	tfx_impl_dx12_buffer_t particle_pool;
	tfxU32 pool_cold_index;
	tfxU32 pool_noise_index;
	//The frame's tfx_CopyGPUParticleUpload, with a view per row stride: dispatches, the 64 byte tick rows, spawns and locations
	tfx_impl_dx12_buffer_t particle_upload;
	tfxU32 upload_tick_index;
	tfxU32 upload_spawn_index;
	tfxU32 upload_location_index;
	tfx_gpu_particle_push_t particle_push_constants;
	//From the last upload and still valid until the next stage update, dispatched once by the next tfx_impl_dx12_Compute
	const tfx_gpu_particle_dispatch_t *particle_dispatches;
	tfxU32 particle_dispatch_count;

	//Baked sprite data playback has its own color ramps, shape frames and emitter properties, from the animation manager
	ID3D12Resource *baked_color_ramps;
	D3D12_RESOURCE_STATES baked_color_ramps_state;
	tfxU32 baked_color_ramps_index;
	tfx_impl_dx12_buffer_t baked_image_data;
	tfx_impl_dx12_buffer_t baked_emitter_properties;
	tfx_impl_dx12_buffer_t baked_sprites;		//Every baked sprite of every effect
	tfx_impl_dx12_buffer_t animation_offsets;	//Uploaded every frame
	tfx_impl_dx12_buffer_t animation_instances;
	tfx_impl_dx12_buffer_t playback_instances;	//The billboard instances the playback compute writes, then a vertex buffer
	tfx_sprite_data_push_t sprite_data_push;
	tfxU32 playback_sprite_count;
	tfx_impl_dx12_buffer_t baked_ribbon_data;	//tfx_ribbon_instance_data_t, every baked ribbon of every effect
	tfx_impl_dx12_buffer_t baked_ribbon_properties;
	tfx_impl_dx12_buffer_t animation_ribbon_offsets;	//Uploaded every frame
	tfx_ribbon_data_push_t ribbon_data_push;
	tfxU32 playback_ribbon_count;

	tfx_impl_dx12_buffer_t graph_lookups;		//tfx_GetGPUGraphLookupsBuffer
	tfx_impl_dx12_ribbon_set_t live_ribbons;	//Uploaded when an update ran
	tfx_impl_dx12_ribbon_set_t baked_ribbons;	//Segments and emitters uploaded once, the rest every frame
} tfx_impl_dx12_t;

//--- Dx12_functions

//Loads the shaders and creates the descriptor heaps, root signatures and pipelines. shader_folder holds each entry point of timelinefx.slang as
//<entry point>.dxil, bindless. The pipelines draw into a render target of render_target_format with a depth buffer of depth_format
bool tfx_impl_dx12_Init(tfx_impl_dx12_t *renderer, ID3D12Device *device, DXGI_FORMAT render_target_format, DXGI_FORMAT depth_format, const char *shader_folder, tfxU32 frames_in_flight);
//Call first every frame, once the GPU has finished the last frame that used frame_index (0 to frames_in_flight - 1). Releases what that frame
//replaced and starts its upload ring again
void tfx_impl_dx12_BeginFrame(tfx_impl_dx12_t *renderer, tfxU32 frame_index);
//Creates the texture array that every library's atlas and the sprite data's atlas are uploaded into, with room for all their layers. Only while the GPU is idle
bool tfx_impl_dx12_CreateAtlas(tfx_impl_dx12_t *renderer, tfxU32 layer_size, tfxU32 layer_count);
//Records an atlas's upload into its layers, with mips built on the CPU. The rest of the library data goes up with each frame in tfx_impl_dx12_UploadFrame
bool tfx_impl_dx12_UploadAtlas(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, const tfx_impl_atlas_t *atlas);
//Records the upload of an animation manager's baked sprite data, color ramps, shape frames and emitter properties, after tfx_BuildAnimationManagerGPUShapeData
bool tfx_impl_dx12_UploadSpriteData(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, tfx_animation_manager animation_manager, tfx_gpu_shapes shapes);
//Records the baked playback's upload every frame, and on update frames the instances, the ribbon data and any library data tfx_UpdateGlobalGPUData changed
bool tfx_impl_dx12_UploadFrame(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, const tfx_impl_frame_t *frame);
//Records the GPU particle kernel, the baked playback and the ribbon geometry, every frame before tfx_impl_dx12_Render
void tfx_impl_dx12_Compute(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, const tfx_impl_uniforms_t *uniforms);
//Records the draws into the render target, depth buffer, viewport and scissor you've set, interpolating with uniforms->timer_lerp
void tfx_impl_dx12_Render(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, const tfx_impl_uniforms_t *uniforms);
//Once the GPU is idle
void tfx_impl_dx12_Shutdown(tfx_impl_dx12_t *renderer);

#ifdef __cplusplus
}
#endif

#endif
