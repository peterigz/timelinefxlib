#include "example_scene.h"
#include <math.h>
#include <string.h>

#define EXAMPLE_SCENE_UPDATE_FREQUENCY 60.0
#define EXAMPLE_SCENE_ORBIT_SPEED 0.00015f		//Radians per millisecond
#define EXAMPLE_SCENE_TORCH_SPEED 0.0012f		//Radians per millisecond
#define EXAMPLE_SCENE_GPU_SPEED 0.0005f
#define EXAMPLE_SCENE_SPARKS_INTERVAL 1200.f	//Milliseconds between bursts
#define EXAMPLE_SCENE_LEAVES_INTERVAL 3000.f
#define EXAMPLE_SCENE_BAKED_INTERVAL 250.0
#define EXAMPLE_SCENE_MAX_BAKED_INSTANCES 256
#define EXAMPLE_SCENE_SECOND_LIBRARY_LAYERS 2	//At least this many layers kept for the second library, so a reload can grow into them

//--- Scene_functions

static tfx_effect_template example_scene__create_template(tfx_library library, const char *path) {
	if (!tfx_IsValidDescriptorPath(library, path)) {
		tfx_Print("There's no effect called \"%s\" in the library", path);
		return NULL;
	}
	return tfx_CreateEffectTemplate(library, path);
}

//Every Run on GPU emitter gets a ring sized for one spawn location, and the stage's pool has to hold every ring: the torch's embers,
//the leaves and the GPU Simulated fountain
static tfxU32 example_scene__set_gpu_capacities(example_scene_t *scene) {
	//Leaf bursts overlap, so the leaves need room for every location that can be alive at once
	float leaves_lifetime = tfx_GetEffectLifetime(scene->leaves.effect_template);
	tfxU32 overlapping_locations = (tfxU32)(leaves_lifetime / EXAMPLE_SCENE_LEAVES_INTERVAL) + 1;
	return tfx_SetTemplateGPUCapacityForLocations(scene->leaves.effect_template, overlapping_locations) + tfx_GetTemplateGPUPoolSlots(scene->gpu_template)
		+ tfx_GetTemplateGPUPoolSlots(scene->torch_template);
}

//Each library loads its shapes into its own list, so they get their own layers of the atlas
static tfx_library example_scene__load_library(const char *path, tfx_impl_shapes_t *shapes) {
	tfx_library library = tfx_LoadEffectLibrary(path, tfx_impl_ShapeLoader, tfx_impl_GetUV, shapes);
	tfxErrorFlags status = library ? tfx_GetLibraryErrorStatus(library) : tfxErrorCode_library_object_could_not_be_created;
	//These two still leave a fully loaded library
	if (status & tfxErrorCode_file_version_out_of_date) {
		tfx_Print("%s was saved in file format %u, resave it in the editor to update it", path, tfx_GetLibraryFileVersion(library));
	}
	status &= ~(tfxErrorCode_file_version_out_of_date | tfxErrorCode_shapes_changed_on_disk);
	if (status != 0) {
		tfx_Print("Could not load %s (error flags 0x%x)", path, (unsigned)status);
		if (library) {
			tfx_FreeLibrary(library);
		}
		return NULL;
	}
	return library;
}

//Packs a set's shapes into layers from first_layer on. A set with no shapes takes no layers
static bool example_scene__build_atlas(example_scene_t *scene, example_scene_shape_set_t *shape_set, tfxU32 first_layer) {
	if (shape_set->shapes.count == 0) {
		memset(&shape_set->atlas, 0, sizeof(tfx_impl_atlas_t));
		shape_set->atlas.first_layer = first_layer;
		return true;
	}
	if (!tfx_impl_BuildShapeAtlas(&shape_set->shapes, &shape_set->atlas, first_layer, scene->atlas_layer_size)) {
		tfx_Print("Could not build the shape atlas");
		return false;
	}
	return true;
}

//The library goes before its shapes, because it holds pointers to them
static void example_scene__free_second_library(example_scene_t *scene) {
	example_scene_shape_set_t *shape_set = &scene->shape_sets[example_scene_shapes_second_library];
	if (scene->second_library) {
		//Also removes its effects from the stage and frees its templates
		tfx_FreeLibrary(scene->second_library);
	}
	scene->second_library = NULL;
	scene->second_library_template = NULL;
	tfx_impl_FreeShapes(&shape_set->shapes);
	tfx_impl_FreeAtlasPixels(&shape_set->atlas);
}

//Called on an update frame before the stage update, so the instances drawn next already have the new layout
static void example_scene__reload_second_library(example_scene_t *scene) {
	example_scene_shape_set_t *shape_set = &scene->shape_sets[example_scene_shapes_second_library];
	tfxU32 first_layer = shape_set->atlas.first_layer;
	example_scene__free_second_library(scene);
	scene->second_library = example_scene__load_library(scene->second_library_path, &shape_set->shapes);
	if (!scene->second_library) {
		example_scene__free_second_library(scene);
		return;
	}
	scene->second_library_template = example_scene__create_template(scene->second_library, "Fireflies");
	if (!example_scene__build_atlas(scene, shape_set, first_layer)) {
		example_scene__free_second_library(scene);
		return;
	}
	if (shape_set->atlas.layer_count > shape_set->reserved_layers) {
		tfx_Print("%s needs %u atlas layers but only %u are kept for it", scene->second_library_path, shape_set->atlas.layer_count, shape_set->reserved_layers);
		example_scene__free_second_library(scene);
		return;
	}
	tfx_UpdateLibraryGPUImageData(scene->second_library);
	if (!tfx_UpdateGlobalGPUData()) {
		tfx_Print("The libraries need more room than the global GPU buffers have");
		example_scene__free_second_library(scene);
		return;
	}
	if (scene->second_library_template) {
		scene->second_library_id = tfx_AddEffectTemplateToStage(scene->stage, scene->second_library_template);
		tfx_SetEffectPosition(scene->stage, scene->second_library_id, 4.f, 1.f, -4.f);
	}
	tfx_Print("Reloaded %s into atlas layers %u to %u", scene->second_library_path, first_layer, first_layer + shape_set->reserved_layers - 1);
}

bool example_scene_Load(example_scene_t *scene, const char *library_path, const char *second_library_path, const char *sprite_data_path, tfxU32 features) {
	memset(scene, 0, sizeof(example_scene_t));
	scene->features = features;
	scene->second_library_path = second_library_path;
	scene->library = example_scene__load_library(library_path, &scene->shape_sets[example_scene_shapes_library].shapes);
	if (!scene->library) {
		return false;
	}
	//The scene still plays without it, with every effect from the one library
	scene->second_library = example_scene__load_library(second_library_path, &scene->shape_sets[example_scene_shapes_second_library].shapes);

	//Templates are made before the atlas and the GPU data so that everything they add is in the upload
	scene->torch_template = example_scene__create_template(scene->library, "Spitting Torch");
	scene->sparks.effect_template = example_scene__create_template(scene->library, "Burst Sparks");
	scene->leaves.effect_template = example_scene__create_template(scene->library, "Leaf Explosion");
	if (features & example_scene_feature_ribbons) {
		scene->ribbon_template = example_scene__create_template(scene->library, "Space Weed");
	}
	if (features & example_scene_feature_gpu_simulation) {
		scene->gpu_template = example_scene__create_template(scene->library, "GPU Simulated");
	}
	if (!scene->torch_template || !scene->sparks.effect_template || !scene->leaves.effect_template) {
		return false;
	}
	if (scene->second_library) {
		scene->second_library_template = example_scene__create_template(scene->second_library, "Fireflies");
	}

	//The sprite data carries its own shapes, which go through the same shape loader into their own layers of the atlas
	if (features & example_scene_feature_baked) {
		scene->animation_manager = tfx_CreateAnimationManager(EXAMPLE_SCENE_MAX_BAKED_INSTANCES, 0);
		tfxErrorFlags sprite_data_status = tfx_LoadSpriteData(sprite_data_path, scene->animation_manager, tfx_impl_ShapeLoader, &scene->shape_sets[example_scene_shapes_sprite_data].shapes);
		if (sprite_data_status != 0) {
			tfx_Print("Could not load %s (error flags 0x%x)", sprite_data_path, (unsigned)sprite_data_status);
			return false;
		}
	}

	//One layer size for every set so they share a texture array, then each set's layers follow the one before
	scene->atlas_layer_size = 0;
	for (tfxU32 index = 0; index != example_scene_shapes_count; ++index) {
		tfxU32 layer_size = tfx_impl_GetAtlasLayerSize(&scene->shape_sets[index].shapes);
		scene->atlas_layer_size = layer_size > scene->atlas_layer_size ? layer_size : scene->atlas_layer_size;
	}
	//The atlas only exists once every shape is loaded, so the uv lookup runs again once it's built
	scene->atlas_layer_count = 0;
	for (tfxU32 index = 0; index != example_scene_shapes_count; ++index) {
		example_scene_shape_set_t *shape_set = &scene->shape_sets[index];
		if (!example_scene__build_atlas(scene, shape_set, scene->atlas_layer_count)) {
			return false;
		}
		shape_set->reserved_layers = shape_set->atlas.layer_count;
		if (index == example_scene_shapes_second_library && shape_set->reserved_layers < EXAMPLE_SCENE_SECOND_LIBRARY_LAYERS) {
			shape_set->reserved_layers = EXAMPLE_SCENE_SECOND_LIBRARY_LAYERS;
		}
		scene->atlas_layer_count += shape_set->reserved_layers;
	}
	tfx_UpdateLibraryGPUImageData(scene->library);
	if (scene->second_library) {
		tfx_UpdateLibraryGPUImageData(scene->second_library);
	}
	if (scene->animation_manager) {
		scene->baked_shapes = tfx_CreateGPUShapesList();
		tfx_BuildAnimationManagerGPUShapeData(scene->animation_manager, scene->baked_shapes, tfx_impl_GetUV);
	}
	//Lays every library out in the global buffers the backend uploads. Effects can only be added from libraries it placed
	if (!tfx_UpdateGlobalGPUData()) {
		tfx_Print("The libraries need more room than the global GPU buffers have");
		return false;
	}
	return true;
}

void example_scene_Start(example_scene_t *scene) {
	tfx_stage_info_t stage_info = tfx_CreateStageInfo(tfxStageSetup_none);
	//Effects from both libraries share one set of buffers and every draw range
	stage_info.global_library_indexes = true;
	if (scene->gpu_template) {
		stage_info.gpu_particle_capacity = example_scene__set_gpu_capacities(scene);
	}
	scene->stage = tfx_CreateStage(stage_info);
	scene->random = tfx_CreateRandom(1);
	scene->timer = tfx_impl_CreateTimer(EXAMPLE_SCENE_UPDATE_FREQUENCY);
	scene->camera = tfx_impl_CreateCamera(tfx_impl_Vec3(0.f, 2.f, 15.f), tfx_impl_Vec3(1.f, 1.f, 0.f), 0.9f);
	scene->camera_distance = 15.f;

	scene->torch_id = tfx_AddEffectTemplateToStage(scene->stage, scene->torch_template);
	float torch_position[3] = { -3.f, 0.f, 0.f };
	scene->torch_location = tfx_AddSpawnLocation(scene->stage, scene->torch_id, torch_position, tfxSpawnLocationAdd_none);

	scene->sparks.center = tfx_impl_Vec3(0.f, 1.f, 0.f);
	scene->sparks.spread = 1.5f;
	scene->sparks.interval = EXAMPLE_SCENE_SPARKS_INTERVAL;

	scene->leaves.effect_id = tfx_AddEffectTemplateToStage(scene->stage, scene->leaves.effect_template);
	scene->leaves.center = tfx_impl_Vec3(6.f, 1.f, 0.f);
	scene->leaves.spread = 0.5f;
	scene->leaves.interval = EXAMPLE_SCENE_LEAVES_INTERVAL;

	if (scene->ribbon_template) {
		scene->ribbon_id = tfx_AddEffectTemplateToStage(scene->stage, scene->ribbon_template);
		tfx_SetEffectPosition(scene->stage, scene->ribbon_id, 0.f, 0.f, -4.f);
	}
	if (scene->gpu_template) {
		scene->gpu_id = tfx_AddEffectTemplateToStage(scene->stage, scene->gpu_template);
		float position[3] = { 0.f, 0.f, 4.f };
		scene->gpu_location = tfx_AddSpawnLocation(scene->stage, scene->gpu_id, position, tfxSpawnLocationAdd_none);
	}
	if (scene->second_library_template) {
		scene->second_library_id = tfx_AddEffectTemplateToStage(scene->stage, scene->second_library_template);
		tfx_SetEffectPosition(scene->stage, scene->second_library_id, 4.f, 1.f, -4.f);
	}
}

//Returns true when a burst is due and fills in where it goes
static bool example_scene__next_burst(example_scene_t *scene, example_scene_burst_t *burst, bool fire_now, float position[3]) {
	if (!fire_now && scene->time < burst->next_burst_time) {
		return false;
	}
	burst->next_burst_time = (float)scene->time + burst->interval;
	position[0] = burst->center.x + tfx_RandomRangeFromTo(&scene->random, -burst->spread, burst->spread);
	position[1] = burst->center.y;
	position[2] = burst->center.z + tfx_RandomRangeFromTo(&scene->random, -burst->spread, burst->spread);
	return true;
}

void example_scene_Update(example_scene_t *scene, double elapsed_milliseconds, float screen_width, float screen_height, tfx_impl_clip_space_flags clip_space_flags, tfx_impl_frame_t *frame, tfx_impl_uniforms_t *uniforms) {
	scene->time += elapsed_milliseconds;
	float time = (float)scene->time;

	if (!scene->orbit_paused) {
		scene->camera_yaw += (float)elapsed_milliseconds * EXAMPLE_SCENE_ORBIT_SPEED;
	}
	tfx_impl_OrbitCamera(&scene->camera, scene->camera_yaw, 0.25f, scene->camera_distance);
	tfx_impl_UpdateCameraMatrices(&scene->camera, screen_width / screen_height, clip_space_flags);

	//Everything that changes the stage goes in before the update, mutators aren't safe while it runs
	tfx_impl_vec3_t camera_front = tfx_impl_CameraFront(&scene->camera);
	tfx_SetStageCamera(scene->stage, &camera_front.x, &scene->camera.position.x);
	float torch_angle = time * EXAMPLE_SCENE_TORCH_SPEED;
	float torch_position[3] = { -4.f + cosf(torch_angle), 0.f, sinf(torch_angle) };
	tfx_UpdateSpawnLocation(scene->stage, scene->torch_location, torch_position);
	float burst_position[3];
	if (example_scene__next_burst(scene, &scene->sparks, scene->burst_requested, burst_position)) {
		tfxEffectID sparks_id = tfx_AddEffectTemplateToStage(scene->stage, scene->sparks.effect_template);
		tfx_SetEffectPositionVec3(scene->stage, sparks_id, burst_position);
	}
	if (example_scene__next_burst(scene, &scene->leaves, false, burst_position)) {
		tfx_AddSpawnLocation(scene->stage, scene->leaves.effect_id, burst_position, tfxSpawnLocationAdd_none);
	}
	scene->burst_requested = false;
	if (scene->gpu_location != tfxINVALID_SPAWN_LOCATION && scene->gpu_template) {
		float gpu_angle = time * EXAMPLE_SCENE_GPU_SPEED;
		float position[3] = { sinf(gpu_angle) * 1.5f, 0.f, 4.f };
		tfx_UpdateSpawnLocation(scene->stage, scene->gpu_location, position);
	}

	tfxU32 ticks = tfx_impl_AdvanceTimer(&scene->timer, elapsed_milliseconds);
	if (ticks > 0) {
		if (scene->reload_requested) {
			scene->reload_requested = false;
			example_scene__reload_second_library(scene);
		}
		//Picks up any library or template change before the update writes its indexes, and returns straight away when there's none
		tfx_UpdateGlobalGPUData();
		//One call per frame that writes the instances once, so the previous instance buffer is always the one drawn last
		tfx_UpdateStageSubSteps(scene->stage, scene->timer.update_length, ticks);
	}
	tfx_impl_GatherFrame(scene->stage, scene->library, ticks > 0, (scene->features & example_scene_feature_ribbons) ? &scene->ribbon_staging : NULL,
		(scene->features & example_scene_feature_gpu_simulation) != 0, frame);

	if (scene->animation_manager) {
		if (scene->time >= scene->next_baked_time) {
			scene->next_baked_time = scene->time + EXAMPLE_SCENE_BAKED_INTERVAL;
			tfxAnimationID animation_id = tfx_AddAnimationInstance(scene->animation_manager, "Burst", 0);
			if (tfx_AnimationIDIsValid(animation_id)) {
				float position[3] = {
					tfx_RandomRangeFromTo(&scene->random, -6.f, 6.f),
					tfx_RandomRangeFromTo(&scene->random, 3.f, 6.f),
					tfx_RandomRangeFromTo(&scene->random, -10.f, -6.f),
				};
				tfx_SetAnimationPosition(scene->animation_manager, animation_id, position);
				tfx_SetAnimationScale(scene->animation_manager, animation_id, tfx_RandomRangeFromTo(&scene->random, 0.5f, 1.2f));
			}
		}
		//Baked instances play back in real time and the compute shader interpolates between their frames, so they don't wait for the stage's ticks
		tfx_UpdateAnimationManager(scene->animation_manager, (float)elapsed_milliseconds);
		tfx_impl_GatherBakedFrame(scene->animation_manager, &frame->baked);
	}
	tfx_impl_SetUniforms(uniforms, &scene->camera, screen_width, screen_height, &scene->timer, scene->stage);
}

void example_scene_Burst(example_scene_t *scene) {
	scene->burst_requested = true;
}

void example_scene_ToggleOrbit(example_scene_t *scene) {
	scene->orbit_paused = !scene->orbit_paused;
}

void example_scene_ReloadSecondLibrary(example_scene_t *scene) {
	scene->reload_requested = true;
}

tfx_impl_atlas_t *example_scene_GetPendingAtlas(example_scene_t *scene) {
	for (tfxU32 index = 0; index != example_scene_shapes_count; ++index) {
		if (scene->shape_sets[index].atlas.pixels) {
			return &scene->shape_sets[index].atlas;
		}
	}
	return NULL;
}

void example_scene_Zoom(example_scene_t *scene, float steps) {
	scene->camera_distance *= powf(0.9f, steps);
}

//A stage goes before the libraries it plays, and the libraries (which own the templates) before the shapes they point to
void example_scene_Shutdown(example_scene_t *scene) {
	if (scene->stage) {
		tfx_FreeStage(scene->stage);
	}
	if (scene->library) {
		tfx_FreeLibrary(scene->library);
	}
	if (scene->second_library) {
		tfx_FreeLibrary(scene->second_library);
	}
	if (scene->animation_manager) {
		tfx_FreeAnimationManager(scene->animation_manager);
		tfx_FreeGPUShapesList(scene->baked_shapes);
	}
	for (tfxU32 index = 0; index != example_scene_shapes_count; ++index) {
		tfx_impl_FreeShapes(&scene->shape_sets[index].shapes);
		tfx_impl_FreeAtlasPixels(&scene->shape_sets[index].atlas);
	}
	tfx_impl_FreeRibbonStaging(&scene->ribbon_staging);
	scene->stage = NULL;
	scene->library = NULL;
	scene->second_library = NULL;
	scene->animation_manager = NULL;
}
