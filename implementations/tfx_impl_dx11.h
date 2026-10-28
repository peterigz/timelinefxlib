#ifndef TFX_IMPL_DX11_HEADER
#define TFX_IMPL_DX11_HEADER

/*
	TimelineFX Direct3D 11 backend. Needs feature level 11.0 and shaders/tfx_dx11.hlsl, see README.md.

	Table of contents, search for the keyword to jump to the section:
	[Dx11_types]		The renderer struct that owns every D3D11 object
	[Dx11_resources]	Buffer, texture and shader compile helpers (in the .cpp)
	[Dx11_functions]	Init, upload, render and shutdown
*/

#include <d3d11.h>
#include "tfx_impl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

//--- Dx11_types

typedef struct tfx_impl_dx11_buffer_s {
	ID3D11Buffer *buffer;
	ID3D11ShaderResourceView *view;
	ID3D11UnorderedAccessView *unordered_view;	//Only for buffers a compute shader writes
	tfxU32 capacity;				//In elements, or bytes for raw buffers
} tfx_impl_dx11_buffer_t;

//What ribbon_compute reads and writes for one source of ribbons: a stage's live ribbons or an animation manager's baked ones
typedef struct tfx_impl_dx11_ribbon_set_s {
	tfx_impl_dx11_buffer_t segments;
	tfx_impl_dx11_buffer_t instances;	//tfx_ribbon_t, uploaded for live ribbons and written by ribbon_playback for baked ones
	tfx_impl_dx11_buffer_t emitters;
	tfx_impl_dx11_buffer_t vertices;	//Built by ribbon_compute every frame
	tfx_impl_dx11_buffer_t indices;		//Raw, written by ribbon_compute and bound as the index buffer
	tfx_impl_dx11_buffer_t graph_lookups;	//Raw, baked only: tfx_GetAnimationRibbonGraphLookupsBufferPointer. Live ribbons read the renderer's shared table
	//The buckets are walked every frame because the geometry is rebuilt with each frame's lerp. One of these is set while there's something to draw
	tfx_stage stage;
	tfx_ribbon_batch_offsets_t stage_offsets;
	tfx_animation_manager animation_manager;
} tfx_impl_dx11_ribbon_set_t;

typedef struct tfx_impl_dx11_s {
	ID3D11Device *device;
	ID3D11DeviceContext *device_context;

	ID3D11VertexShader *billboard_vertex_shader;
	ID3D11PixelShader *billboard_pixel_shader;
	ID3D11BlendState *blend_state;
	ID3D11DepthStencilState *depth_stencil_state;
	ID3D11RasterizerState *rasterizer_state;
	ID3D11SamplerState *sampler_state;
	ID3D11Buffer *frame_constants;	//tfx_impl_uniforms_t at b0
	ID3D11Buffer *draw_constants;	//The draw's instance offsets at b1

	ID3D11Texture2D *shape_atlas;
	ID3D11ShaderResourceView *shape_atlas_view;
	ID3D11Texture2D *color_ramps;
	ID3D11ShaderResourceView *color_ramps_view;
	tfx_impl_dx11_buffer_t image_data;			//tfx_gpu_image_data_t per shape frame
	tfx_impl_dx11_buffer_t emitter_properties;	//tfx_gpu_particle_properties_t per emitter

	//Two instance buffers that swap on every update: the latest instances and the ones captured_index points into
	tfx_impl_dx11_buffer_t instance_buffers[2];
	tfxU32 current_instance_buffer;
	tfxU32 instance_capacity;
	tfxU32 instance_count;			//CPU and GPU instances
	//Each layer is drawn as two ranges, its CPU instances then its GPU instances
	tfxU32 layer_start[tfxLAYERS];
	tfxU32 layer_count[tfxLAYERS];
	tfxU32 gpu_layer_start[tfxLAYERS];
	tfxU32 gpu_layer_count[tfxLAYERS];

	ID3D11ComputeShader *particle_kernel;		//tfx_gpu_particles.slang built by slangc, null without one
	ID3D11Buffer *particle_push;				//tfx_gpu_particle_push_t at b0
	tfx_impl_dx11_buffer_t particle_pool;		//Raw, the hot, cold and noise sections in one buffer
	tfx_impl_dx11_buffer_t particle_upload;	//Raw and dynamic, the frame's tfx_CopyGPUParticleUpload
	tfx_gpu_particle_push_t particle_push_constants;
	//From the last upload and still valid until the next stage update, dispatched once by the next tfx_impl_dx11_Compute
	const tfx_gpu_particle_dispatch_t *particle_dispatches;
	tfxU32 particle_dispatch_count;

	//Baked sprite data playback has its own color ramps, shape frames and emitter properties, from the animation manager
	ID3D11ComputeShader *sprite_playback_shader;
	ID3D11VertexShader *billboard_static_vertex_shader;
	ID3D11Buffer *sprite_data_constants;		//tfx_sprite_data_push_t at b3
	ID3D11Texture2D *baked_color_ramps;
	ID3D11ShaderResourceView *baked_color_ramps_view;
	tfx_impl_dx11_buffer_t baked_image_data;
	tfx_impl_dx11_buffer_t baked_emitter_properties;
	tfx_impl_dx11_buffer_t baked_sprites;		//Raw, every baked sprite of every effect
	tfx_impl_dx11_buffer_t animation_offsets;	//Raw and dynamic, uploaded every frame
	tfx_impl_dx11_buffer_t animation_instances;
	tfx_impl_dx11_buffer_t playback_instances;	//Raw, the billboard instances the playback compute writes
	tfx_sprite_data_push_t sprite_data_push;
	tfxU32 playback_sprite_count;
	ID3D11ComputeShader *ribbon_playback_shader;
	ID3D11Buffer *ribbon_data_constants;		//tfx_ribbon_data_push_t at b4
	tfx_impl_dx11_buffer_t baked_ribbon_data;	//tfx_ribbon_instance_data_t, every baked ribbon of every effect
	tfx_impl_dx11_buffer_t baked_ribbon_properties;
	tfx_impl_dx11_buffer_t animation_ribbon_offsets;	//Raw and dynamic, uploaded every frame
	tfx_ribbon_data_push_t ribbon_data_push;
	tfxU32 playback_ribbon_count;

	ID3D11ComputeShader *ribbon_compute_shader;
	ID3D11VertexShader *ribbon_vertex_shader;
	ID3D11PixelShader *ribbon_pixel_shader;
	ID3D11Buffer *ribbon_constants;	//A bucket's tfx_ribbon_bucket_globals_t at b2
	tfx_impl_dx11_buffer_t graph_lookups;		//Raw, tfx_GetGPUGraphLookupsBuffer
	tfx_impl_dx11_ribbon_set_t live_ribbons;	//Uploaded when an update ran
	tfx_impl_dx11_ribbon_set_t baked_ribbons;	//Segments and emitters uploaded once, the rest every frame
} tfx_impl_dx11_t;

//--- Dx11_functions

//Compiles the shaders and creates the states. shader_path is the path to tfx_dx11.hlsl, kernel_path the compiled GPU particle kernel or null to go without
bool tfx_impl_dx11_Init(tfx_impl_dx11_t *renderer, ID3D11Device *device, ID3D11DeviceContext *device_context, const char *shader_path, const char *kernel_path);
//Uploads the shape atlas, color ramps, shape frame data and emitter properties. Call after tfx_UpdateLibraryGPUImageData and creating your effect templates
bool tfx_impl_dx11_UploadLibrary(tfx_impl_dx11_t *renderer, tfx_library library, const tfx_impl_atlas_t *atlas);
//Uploads an animation manager's baked sprite data, color ramps, shape frames and emitter properties, after tfx_BuildAnimationManagerGPUShapeData
bool tfx_impl_dx11_UploadSpriteData(tfx_impl_dx11_t *renderer, tfx_animation_manager animation_manager, tfx_gpu_shapes shapes);
//Uploads the baked playback every frame, and the instances and ribbon data when the frame ran an update
bool tfx_impl_dx11_UploadFrame(tfx_impl_dx11_t *renderer, const tfx_impl_frame_t *frame);
//Runs the GPU particle kernel, the baked playback and builds the ribbon geometry, every frame before tfx_impl_dx11_Render
void tfx_impl_dx11_Compute(tfx_impl_dx11_t *renderer, const tfx_impl_uniforms_t *uniforms);
//Draws into whatever render target and viewport are bound, every frame, interpolating with uniforms->timer_lerp
void tfx_impl_dx11_Render(tfx_impl_dx11_t *renderer, const tfx_impl_uniforms_t *uniforms);
void tfx_impl_dx11_Shutdown(tfx_impl_dx11_t *renderer);

#ifdef __cplusplus
}
#endif

#endif
