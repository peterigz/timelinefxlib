#ifndef TFX_EXAMPLE_SCENE_HEADER
#define TFX_EXAMPLE_SCENE_HEADER

/*
	The demo scene every platform example plays, so each example's main.cpp only has to open a window, create a device
	and call its backend. It's example code showing how a game uses the library, see README.md for what it shows.

	Table of contents, search for the keyword to jump to the section:
	[Scene_types]		The scene struct and the features a backend supports
	[Scene_functions]	Load, start, update, input and shutdown
*/

#include "tfx_impl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

//--- Scene_types

//What the backend can draw. Effects that need something the backend doesn't have yet are left out of the scene
typedef enum example_scene_features {
	example_scene_feature_none = 0,
	example_scene_feature_ribbons = 1 << 0,
	example_scene_feature_gpu_simulation = 1 << 1,
	example_scene_feature_baked = 1 << 2,
} example_scene_features;

//An effect that fires every so often at a random spot near its center
typedef struct example_scene_burst_s {
	tfx_effect_template effect_template;
	tfxEffectID effect_id;			//The one running effect when it bursts at user spawn locations
	tfx_impl_vec3_t center;
	float spread;					//Bursts land up to this far from the center on x and z
	float interval;					//Milliseconds between bursts
	float next_burst_time;
} example_scene_burst_t;

typedef struct example_scene_s {
	tfxU32 features;
	tfx_impl_shapes_t shapes;
	tfx_impl_ribbon_staging_t ribbon_staging;
	tfx_library library;
	tfx_stage stage;
	tfx_random_t random;
	tfx_impl_timer_t timer;
	tfx_impl_camera_t camera;
	float camera_yaw;
	float camera_distance;
	bool orbit_paused;
	double time;					//Milliseconds since the scene started

	//Moved around a circle with tfx_SetEffectPosition every frame
	tfx_effect_template torch_template;
	tfxEffectID torch_id;
	//Fire and forget: every burst adds a new effect, which frees itself once it's finished
	example_scene_burst_t sparks;
	//Spawn at user locations: one effect, and every burst adds a location that removes itself once the effect has played out there
	example_scene_burst_t leaves;
	bool burst_requested;
	//Ribbons
	tfx_effect_template ribbon_template;
	tfxEffectID ribbon_id;
	//Run on GPU emitters always spawn at user locations, this one has a single location that's moved around
	tfx_effect_template gpu_template;
	tfxEffectID gpu_id;
	tfxSpawnLocationID gpu_location;
	//Baked sprite data, played back by a compute shader. Each burst is a cheap animation instance at a random spot
	tfx_animation_manager animation_manager;
	tfx_gpu_shapes baked_shapes;
	double next_baked_time;
} example_scene_t;

//--- Scene_functions

//Loads the library (and the baked sprite data with example_scene_feature_baked), creates the templates and builds the shape atlas,
//which the backend uploads before calling example_scene_Start
bool example_scene_Load(example_scene_t *scene, const char *library_path, const char *sprite_data_path, tfxU32 features, tfx_impl_atlas_t *atlas);
//Creates the stage and adds the effects
void example_scene_Start(example_scene_t *scene);
//Moves things, updates the stage and gathers what the backend needs to upload and draw this frame
void example_scene_Update(example_scene_t *scene, double elapsed_milliseconds, float screen_width, float screen_height, tfx_impl_clip_space_flags clip_space_flags, tfx_impl_frame_t *frame, tfx_impl_uniforms_t *uniforms);
void example_scene_Burst(example_scene_t *scene);
void example_scene_ToggleOrbit(example_scene_t *scene);
//Positive steps zoom in
void example_scene_Zoom(example_scene_t *scene, float steps);
void example_scene_Shutdown(example_scene_t *scene);

#ifdef __cplusplus
}
#endif

#endif
