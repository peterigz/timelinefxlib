#ifndef TFX_IMPL_GL_HEADER
#define TFX_IMPL_GL_HEADER

/*
	TimelineFX OpenGL backend. Needs a 4.3 core context and the shaders slangc and spirv-cross build from
	../Shaders/timelinefx.slang, see README.md. It includes no GL header and loads every function it calls through the
	loader you pass to tfx_impl_gl_Init, so it works with GLFW, SDL or any other window library.

	Table of contents, search for the keyword to jump to the section:
	[Gl_types]			The renderer struct that owns every GL object
	[Gl_resources]		Function loading, buffer, texture and shader helpers (in the .cpp)
	[Gl_functions]		Init, upload, render and shutdown
*/

#include "tfx_impl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

//--- Gl_types

//glfwGetProcAddress and SDL_GL_GetProcAddress both fit. It has to return the GL 1.1 functions too, which wglGetProcAddress alone doesn't
typedef void (*tfx_impl_gl_proc)(void);
typedef tfx_impl_gl_proc (*tfx_impl_gl_loader)(const char *name);

typedef struct tfx_impl_gl_buffer_s {
	tfxU32 buffer;
	tfxU32 capacity;				//In bytes
} tfx_impl_gl_buffer_t;

//What ribbon_compute reads and writes for one source of ribbons: a stage's live ribbons or an animation manager's baked ones
typedef struct tfx_impl_gl_ribbon_set_s {
	tfx_impl_gl_buffer_t segments;
	tfx_impl_gl_buffer_t instances;		//tfx_ribbon_t, uploaded for live ribbons and written by ribbon_playback for baked ones
	tfx_impl_gl_buffer_t emitters;
	tfx_impl_gl_buffer_t vertices;		//Built by ribbon_compute every frame
	tfx_impl_gl_buffer_t indices;		//Written by ribbon_compute and bound as the element buffer
	tfx_impl_gl_buffer_t graph_lookups;	//Baked only: tfx_GetAnimationRibbonGraphLookupsBufferPointer. Live ribbons read the renderer's shared table
	//The buckets are walked every frame because the geometry is rebuilt with each frame's lerp. One of these is set while there's something to draw
	tfx_stage stage;
	tfx_ribbon_batch_offsets_t stage_offsets;
	tfx_animation_manager animation_manager;
} tfx_impl_gl_ribbon_set_t;

typedef struct tfx_impl_gl_s {
	//Each layer picks its programs from its tfx_GetLayerFlags, by [writes depth][lit]. Baked billboards are always blended and unlit
	tfxU32 billboard_programs[2][2];
	tfxU32 billboard_baked_program;
	tfxU32 ribbon_programs[2][2];
	tfxU32 particle_kernel;			//simulate_particles from timelinefx.slang
	tfxU32 ribbon_compute_program;
	tfxU32 sprite_playback_program;
	tfxU32 ribbon_playback_program;
	tfxU32 empty_vertex_array;		//Core profile draws need one bound, though nothing reads vertex attributes
	tfxU32 max_storage_block_size;

	tfx_impl_gl_buffer_t frame_constants;	//tfx_impl_uniforms_t, uniform binding 5, and 10 for the lit fragments
	tfx_impl_gl_buffer_t draw_constants;	//The draw's instance offsets, or a ribbon bucket's globals, uniform binding 6
	tfx_impl_gl_buffer_t compute_constants;	//Each dispatch's block, uniform binding 6

	tfxU32 shape_atlas;
	tfxU32 shape_atlas_layer_size;
	tfxU32 shape_atlas_layer_count;
	tfxU32 shape_atlas_mip_count;
	//Every loaded library's data, from the global getters, so one set of buffers draws effects from any library
	tfxU32 color_ramps;
	tfxU32 color_ramp_layer_capacity;
	tfxU32 *color_ramp_layer_versions;			//The global layer versions last uploaded, one per layer of color_ramps
	tfx_impl_gl_buffer_t image_data;			//tfx_gpu_image_data_t per shape frame
	tfx_impl_gl_buffer_t emitter_properties;	//tfx_gpu_particle_properties_t per emitter
	tfxU32 properties_version;
	tfxU32 shapes_version;
	tfxU32 color_ramps_version;
	tfxU32 graph_lookups_version;

	//Two instance buffers that swap on every update: the latest instances and the ones captured_index points into
	tfx_impl_gl_buffer_t instance_buffers[2];
	tfxU32 current_instance_buffer;
	tfxU32 instance_capacity;
	tfxU32 instance_count;			//CPU and GPU instances
	//Each layer is drawn as two ranges, its CPU instances then its GPU instances
	tfxU32 layer_start[tfxLAYERS];
	tfxU32 layer_count[tfxLAYERS];
	tfxU32 gpu_layer_start[tfxLAYERS];
	tfxU32 gpu_layer_count[tfxLAYERS];
	tfxLayerFlags layer_flags[tfxLAYERS];

	tfx_impl_gl_buffer_t particle_pool;		//The hot, cold and noise sections in one buffer
	tfx_impl_gl_buffer_t particle_upload;	//The frame's tfx_CopyGPUParticleUpload
	tfx_gpu_particle_push_t particle_push_constants;
	//From the last upload and still valid until the next stage update, dispatched once by the next tfx_impl_gl_Compute
	const tfx_gpu_particle_dispatch_t *particle_dispatches;
	tfxU32 particle_dispatch_count;

	//Baked sprite data playback has its own color ramps, shape frames and emitter properties, from the animation manager
	tfxU32 baked_color_ramps;
	tfx_impl_gl_buffer_t baked_image_data;
	tfx_impl_gl_buffer_t baked_emitter_properties;
	tfx_impl_gl_buffer_t baked_sprites;			//Every baked sprite of every effect
	tfx_impl_gl_buffer_t animation_offsets;		//Uploaded every frame
	tfx_impl_gl_buffer_t animation_instances;
	tfx_impl_gl_buffer_t playback_instances;	//The billboard instances the playback compute writes
	tfx_sprite_data_push_t sprite_data_push;
	tfxU32 playback_sprite_count;
	tfx_impl_gl_buffer_t baked_ribbon_data;		//tfx_ribbon_instance_data_t, every baked ribbon of every effect
	tfx_impl_gl_buffer_t baked_ribbon_properties;
	tfx_impl_gl_buffer_t animation_ribbon_offsets;	//Uploaded every frame
	tfx_ribbon_data_push_t ribbon_data_push;
	tfxU32 playback_ribbon_count;

	tfx_impl_gl_buffer_t graph_lookups;			//tfx_GetGPUGraphLookupsBuffer
	tfx_impl_gl_ribbon_set_t live_ribbons;		//Uploaded when an update ran
	tfx_impl_gl_ribbon_set_t baked_ribbons;		//Segments and emitters uploaded once, the rest every frame
} tfx_impl_gl_t;

//--- Gl_functions

//Loads the GL functions, compiles the shaders and creates the uniform buffers, with the context current. shader_folder holds each
//entry point of timelinefx.slang as <entry point>.glsl, built through SPIR-V and spirv-cross. Messages from a debug context go to tfx_Print
bool tfx_impl_gl_Init(tfx_impl_gl_t *renderer, tfx_impl_gl_loader loader, const char *shader_folder);
//Creates the texture array that every library's atlas and the sprite data's atlas are uploaded into, with room for all their layers
bool tfx_impl_gl_CreateAtlas(tfx_impl_gl_t *renderer, tfxU32 layer_size, tfxU32 layer_count);
//Writes an atlas into its layers of the texture array and generates their mips. The rest of the library data goes up with each
//frame in tfx_impl_gl_UploadFrame
bool tfx_impl_gl_UploadAtlas(tfx_impl_gl_t *renderer, const tfx_impl_atlas_t *atlas);
//Uploads an animation manager's baked sprite data, color ramps, shape frames and emitter properties, after tfx_BuildAnimationManagerGPUShapeData
bool tfx_impl_gl_UploadSpriteData(tfx_impl_gl_t *renderer, tfx_animation_manager animation_manager, tfx_gpu_shapes shapes);
//Uploads the baked playback every frame, and on update frames the instances, the ribbon data and any library data tfx_UpdateGlobalGPUData changed
bool tfx_impl_gl_UploadFrame(tfx_impl_gl_t *renderer, const tfx_impl_frame_t *frame);
//Runs the GPU particle kernel, the baked playback and builds the ribbon geometry, every frame before tfx_impl_gl_Render
void tfx_impl_gl_Compute(tfx_impl_gl_t *renderer, const tfx_impl_uniforms_t *uniforms);
//Draws into the bound framebuffer and viewport, every frame, interpolating with uniforms->timer_lerp. Leaves depth writes on
void tfx_impl_gl_Render(tfx_impl_gl_t *renderer, const tfx_impl_uniforms_t *uniforms);
//With the context still current
void tfx_impl_gl_Shutdown(tfx_impl_gl_t *renderer);

#ifdef __cplusplus
}
#endif

#endif
