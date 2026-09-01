#ifndef TFX_LIBRARY_HEADER
#define TFX_LIBRARY_HEADER

#define TFX_VERSION_NAME "Alpha"
#define TFX_VERSION_MAJOR 0
#define TFX_VERSION_MINOR 33
#define TFX_VERSION_PATCH 0

#define TFX_STRINGIFY(x) #x
#if defined(__clang__)
#define TFX_DISABLE_COMPILER_WARNING(w) \
	_Pragma("clang diagnostic push") \
	_Pragma(TFX_STRINGIFY(clang diagnostic ignored w))
#define TFX_ENABLE_COMPILER_WARNING() \
	_Pragma("clang diagnostic pop")
#elif defined(__GNUC__)
#define TFX_DISABLE_COMPILER_WARNING(w)
#define TFX_ENABLE_COMPILER_WARNING()
#else
#define TFX_DISABLE_COMPILER_WARNING(w)
#define TFX_ENABLE_COMPILER_WARNING()
#endif

//Profiling is done through Tracy (https://github.com/wolfpld/tracy). Define tfxTRACY
//(along with Tracy's own TRACY_ENABLE) to turn the tfxPROFILE markers scattered through
//the library into Tracy zones. Without it they compile away to nothing, so a shipping
//build pays for no instrumentation at all.
//#define tfxTRACY
#define TFX_THREAD_SAFE
//#define TFX_EXTRA_DEBUGGING
#define SSE41		//Steam survey currently has this at 99.83% coverage 12 April 2025. I will probably make this the minimum requirement

//Override this for more layers, although currently the editor is fixed at 4
#ifndef tfxLAYERS
#define tfxLAYERS 4
#endif

//Enable this to process 8 particles at a time.
//#define tfxUSEAVX

//Enable fused multiply add in simd calculations
//#define tfxUSEFMA

/*
	Timeline FX C++ library

	This library is for implementing particle effects into your games and applications.

	This library is render agnostic, so you will have to provide your own means to render the particles. There are various API functions in the library that help you do this.

	Currently tested on Windows and MacOS, Intel and Mac based ARM processors.

	Table of contents
	Sections in this header file, you can search for the following keywords to jump to that section:

	[Header_Includes_and_Typedefs]		Base typedefs, version + platform/SIMD detection, and the tfxAPI / tfxINTERNAL markers
	[Enums]								All the definitions for enums and bit flags
	[forward_declarations]				Opaque handles, struct typedefs and the public value-type structs (tfx_random_t, tfx_stage_info_t, ...)
	[Callback_typedefs]					Callback function pointer typedefs
	[API functions]						The main functions for use by users of the library
		-[Random numbers]				Seeded random number generation
		-[Initialisation_functions]		Startup and shutdown timelinefx
		-[Global_variable_access]		Any functions that give you access to global variables relating to timelinefx
		-[Library_functions]			Functions for loading and accessing timelinefx libraries
		-[Particle_Manager_functions]	Create and update functions for effect managers where the main work is done to update particles every frame
		-[Animation_manager]			Animation manager functions for playing pre-recorded effect data
		-[Effect_templates]				Functions for working with effect templates which help modify effects in the library without actually changing the base effect in the library.
		-[Editing_graphs]				Functions to configure effect/emitter graphs
		-[General_helpers]				General math functions and other helpers.

	Everything internal - the Zest pocket allocator, containers, SIMD, work queues, hashing, profiling, file IO and all
	tfxINTERNAL / tfxAPI_EDITOR declarations - now lives in timelinefx_internal.h, included by timelinefx.cpp and the
	editor but never by a game.
*/

/*    Functions come in 3 flavours:
1) INTERNAL where they're only meant for internal use by the library and not for any use outside it. Note that these functions are declared as static.
2) API where they're meant for access within your games that you're developing. These functions are c compatible.
3) EDITOR where they can be accessed from outside the library but really they're mainly useful for editing the effects such as in in the TimelineFX Editor. These
   functions are c++ compatabile only and currently not available if you're including the library in a c project.

All functions in the library will be marked this way for clarity and naturally the API functions will all be properly documented.
*/

#ifdef __cplusplus
#define tfxAPI extern "C"
#else
#define tfxAPI 
#endif    
#if defined(__GNUC__) || defined(__clang__)
#define tfxINTERNAL static __attribute__((unused))
#else
#define tfxINTERNAL static
#endif

//----------------------------------------------------------
//Header_Includes_and_Typedefs
//----------------------------------------------------------
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64)
#define tfxX86
#include <immintrin.h>
#elif defined(__arm__) || defined(__aarch64__)
#include <arm_neon.h>
#include <mach/mach_time.h>
#define tfxARM
#endif

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#define tfxWINDOWS
#elif __APPLE__
#define tfxMAC
#elif __linux__
#define tfxLINUX
#endif

#include <stdint.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>

//type defs
typedef uint16_t tfxU16;
typedef unsigned short tfxHalf;
typedef uint32_t tfxU32;
typedef uint32_t tfxIndex;
typedef unsigned char tfxU8;
typedef unsigned int tfxEmitterID;
typedef int32_t tfxS32;
typedef uint64_t tfxU64;
typedef int64_t tfxS64;
typedef tfxU32 tfxEffectID;
typedef tfxU32 tfxAnimationID;
typedef tfxU64 tfxKey;
typedef tfxU32 tfxParticleID;
typedef short tfxShort;
typedef unsigned short tfxUShort;

typedef struct tfx_float16x4_s {
	union {
		struct {
			tfxU16 x : 16;
			tfxU16 y : 16;
			tfxU16 z : 16;
			tfxU16 w : 16;
		};
		struct {
			tfxU32 xy : 32;
			tfxU32 zw : 32;
		};
		struct { tfxU64 packed; };
	};
} tfx_float16x4_t;

typedef struct tfx_float16x2_s {
	union {
		struct {
			tfxU16 x : 16;
			tfxU16 y : 16;
		};
		struct { tfxU32 packed; };
	};
} tfx_float16x2_t;

typedef struct tfx_float8x4_s {
	union {
		struct {
			tfxU8 x : 8;
			tfxU8 y : 8;
			tfxU8 z : 8;
			tfxU8 w : 8;
		};
		struct { tfxU32 packed; };
	};
} tfx_float8x4_t;

typedef struct tfx_float32x2_s {
	float x, y;
}tfx_float32x2_t;

typedef struct tfx_float32x3_s {
	float x, y, z;
}tfx_float32x3_t;

typedef struct tfx_float32x4_s {
	float x, y, z, w;
}tfx_float32x4_t;

typedef struct tfx_instance_s {		//48 bytes
	tfx_float32x4_t position;						//The position of the billboard with stretch in w
	tfxU64 quaternion;								//Rotation of the billboard stored as a 16-bit snorm quaternion
	tfx_float16x2_t size;							//Size/Scale of the sprite
	tfx_float8x4_t alignment;						//normalised alignment vector 3 8bit floats packed into 32 bits. Free byte here.
	tfx_float16x2_t intensity_gradient_map;			//Multiplier for the color and life of particle
	tfx_float8x4_t curved_alpha_life;				//Sharpness and dissolve amount value for fading the image plus the age of the particle value packed into 3 bit unorms. Free byte here.
	tfxU32 indexes;									//[gpu property index, capture flag (1 bit << 15), image data index max 8191 images]
	tfxU32 captured_index;							//Index to the sprite in the buffer from the previous frame for interpolation
} tfx_instance_t;

typedef struct tfx_ribbon_vertex_s {
	tfx_float32x3_t position;
	tfxU32 segment_index;
	tfx_float32x2_t uv_offset_scale;
	tfxU32 ribbon_index;
	tfxU32 clipped;
} tfx_ribbon_vertex_t;

typedef struct tfx_gpu_particle_properties_s {
	tfx_float32x2_t image_handle;			
	tfxU32 color_ramp_indexes;			//[Row of color ramp bitmap, texture array]
	tfxU32 flags;						//Flags like billboard alignment type
	tfxU32 start_frame_index;
	float animation_frames;
	int padding[2];
} tfx_gpu_particle_properties_t;

typedef struct tfx_gpu_image_data_s {
	tfx_float32x4_t uv;
	tfxU64 uv_packed;
	tfx_float32x2_t image_size;
	tfxU32 texture_array_index;
	float animation_frames;
	float padding;	//Was max_radius which isn't used anymore
}tfx_gpu_image_data_t;

//------------------------------------------------------------
//Section: Enums
//------------------------------------------------------------

typedef enum tfx_color_format {
	tfx_color_format_rgba8,
	tfx_color_format_rgba16f,
	tfx_color_format_rgba32f,
} tfx_color_format;

typedef enum {
	tfxStageSetup_none,
	tfxStageSetup_group_sprites_by_effect,
} tfx_stage_setup;

typedef enum {
	tfxEffect_global_life_index,
	tfxEffect_global_amount_index,
	tfxEffect_global_velocity_index,
	tfxEffect_global_noise_index,
	tfxEffect_global_width_index,
	tfxEffect_global_height_index,
	tfxEffect_global_weight_index,
	tfxEffect_global_roll_spin_index,
	tfxEffect_global_pitch_spin_index,
	tfxEffect_global_yaw_spin_index,
	tfxEffect_global_stretch_index,
	tfxEffect_global_overall_scale_index,
	tfxEffect_global_intensity_index,
	tfxEffect_global_splatter_index,
	tfxEffect_global_emitter_width_index,
	tfxEffect_global_emitter_height_index,
	tfxEffect_global_emitter_depth_index,
	tfxEffectGraphs_max_index,
} tfx_global_graph_index;

typedef enum {
	tfxEmitter_property_emission_pitch_index,
	tfxEmitter_property_emission_yaw_index,
	tfxEmitter_property_emission_range_index,
	tfxEmitter_property_splatter_index,
	tfxEmitter_property_width_index,        //Also used for linear extrusion for paths as well
	tfxEmitter_property_height_index,
	tfxEmitter_property_depth_index,
	tfxEmitter_property_extrusion_index,
	tfxEmitter_property_arc_size_index,
	tfxEmitter_property_arc_offset_index,

	tfxEmitter_base_life_index,
	tfxEmitter_base_amount_index,
	tfxEmitter_base_velocity_index,
	tfxEmitter_base_width_index,
	tfxEmitter_base_height_index,
	tfxEmitter_base_weight_index,
	tfxEmitter_base_pitch_spin_index,
	tfxEmitter_base_yaw_spin_index,
	tfxEmitter_base_roll_spin_index,

	tfxEmitter_variation_life_index,
	tfxEmitter_variation_amount_index,
	tfxEmitter_variation_velocity_index,
	tfxEmitter_variation_width_index,
	tfxEmitter_variation_height_index,
	tfxEmitter_variation_weight_index,
	tfxEmitter_variation_path_trajectory_scale_index,
	tfxEmitter_variation_pitch_spin_index,
	tfxEmitter_variation_yaw_spin_index,
	tfxEmitter_variation_roll_spin_index,
	tfxEmitter_variation_noise_resolution_index,
	tfxEmitter_variation_motion_randomness_index,

	tfxEmitter_overtime_red_index,
	tfxEmitter_overtime_green_index,
	tfxEmitter_overtime_blue_index,
	tfxEmitter_overtime_blendfactor_index,
	tfxEmitter_overtime_velocity_adjuster_index,
	tfxEmitter_overtime_intensity_index,
	tfxEmitter_overtime_alpha_sharpness_index,
	tfxEmitter_overtime_curved_alpha_index,
	tfxEmitter_overtime_heat_response_index,
	tfxEmitter_overtime_gradient_mapper_index,
	tfxEmitter_overtime_velocity_index,
	tfxEmitter_overtime_width_index,
	tfxEmitter_overtime_height_index,
	tfxEmitter_overtime_weight_index,
	tfxEmitter_overtime_pitch_spin_index,
	tfxEmitter_overtime_yaw_spin_index,
	tfxEmitter_overtime_roll_spin_index,
	tfxEmitter_overtime_stretch_index,
	tfxEmitter_overtime_velocity_turbulance_index,
	tfxEmitter_overtime_direction_turbulance_index,
	tfxEmitter_overtime_direction_index,
	tfxEmitter_overtime_noise_resolution_index,
	tfxEmitter_overtime_motion_randomness_index,

	tfxEmitter_factor_life_index,
	tfxEmitter_factor_size_index,
	tfxEmitter_factor_velocity_index,
	tfxEmitter_factor_intensity_index,

	tfxEmitterGraphs_max_index,
} tfx_emitter_graph_index;

typedef enum {
	tfxRibbon_property_splatter_index,
	tfxRibbon_property_width_index,        //Also used for linear extrusion for paths as well
	tfxRibbon_property_height_index,
	tfxRibbon_property_depth_index,
	tfxRibbon_property_extrusion_index,
	tfxRibbon_property_arc_size_index,
	tfxRibbon_property_arc_offset_index,

	tfxRibbon_base_life_index,
	tfxRibbon_base_amount_index,
	tfxRibbon_base_width_index,

	tfxRibbon_variation_life_index,
	tfxRibbon_variation_amount_index,
	tfxRibbon_variation_width_index,

	tfxRibbon_overtime_red_index,
	tfxRibbon_overtime_green_index,
	tfxRibbon_overtime_blue_index,
	tfxRibbon_overtime_blendfactor_index,
	tfxRibbon_overtime_intensity_index,
	tfxRibbon_overtime_alpha_sharpness_index,
	tfxRibbon_overtime_curved_alpha_index,
	tfxRibbon_overtime_gradient_mapper_index,
	tfxRibbon_overtime_heat_response_index,
	tfxRibbon_overtime_width_index,
	tfxRibbon_overtime_overall_scale_index,
	tfxRibbon_overtime_uv_offset_y_index,
	tfxRibbon_overtime_uv_scale_y_index,
	tfxRibbon_overtime_clip_offset_index,
	tfxRibbon_overtime_clip_size_index,
	tfxRibbon_overtime_morph_amount_index,
	tfxRibbon_overtime_noise_amount_index,
	tfxRibbon_overtime_noise_scroll_index,
	tfxRibbon_overtime_lag_amount_index,

	tfxRibbon_overlength_intensity_index,
	tfxRibbon_overlength_alpha_sharpness_index,
	tfxRibbon_overlength_curved_alpha_index,
	tfxRibbon_overlength_gradient_map_index,
	tfxRibbon_overlength_width_index,
	tfxRibbon_overlength_fixed_angle_index,
	tfxRibbon_overlength_morph_bias_index,
	tfxRibbon_overlength_noise_envelope_index,
	tfxRibbon_overlength_lag_profile_index,


	tfxRibbonGraphs_max_index,

	tfxRibbon_property_start_index = 0,
	tfxRibbon_base_start_index = tfxRibbon_base_life_index,
	tfxRibbon_variation_start_index = tfxRibbon_variation_life_index,
	tfxRibbon_overtime_start_index = tfxRibbon_overtime_red_index,
	tfxRibbon_property_end_index = tfxRibbon_property_arc_offset_index + 1,
	tfxRibbon_base_end_index = tfxRibbon_base_width_index + 1,
	tfxRibbon_variation_end_index = tfxRibbon_variation_width_index + 1,
	tfxRibbon_overtime_end_index = tfxRibbon_overtime_lag_amount_index + 1,
	tfxRibbon_overlength_start = tfxRibbon_overlength_intensity_index,
	tfxRibbon_overlength_end = tfxRibbon_overlength_lag_profile_index + 1,
} tfx_ribbon_graph_index;

//tfx_effect_descriptor_t type - effect contains emitters, and emitters spawn particles, but they both share the same struct for simplicity
typedef enum {
	tfxEffectType,
	tfxEmitterType,
	tfxRibbonType,
	tfxFolder,
	//Not a descriptor type an effect can have - it only ever tags a tfx_graph_list_t so that the list knows
	//how many graphs it holds and which initialiser rebuilds it. Appended rather than inserted because the
	//preceding values are saved in the file as ordinals.
	tfxForceType,
	tfxMaxDescriptorTypes
} tfx_effect_descriptor_type;

typedef enum {
	tfxColorInterpolation_linear_srgb = 0,
	tfxColorInterpolation_oklch,
	tfxColorInterpolation_hsl,
	tfxColorInterpolation_linear_rgb,
	tfxColorInterpolation_max
} tfx_color_interpolation_mode;

typedef enum {
	tfxGraphEasingType_constant                                 = 0,
	tfxGraphEasingType_smoothstep                               = 17,
	tfxGraphEasingType_out_in									= 18,
	tfxGraphEasingType_in										= 4,
	tfxGraphEasingType_out										= 5,
	tfxGraphEasingType_in_out									= 6,
	tfxGraphEasingType_linear                                   = 16,
	//Unused
	tfxGraphEasingType_ease_in_quad                             = 1,
	tfxGraphEasingType_ease_out_quad                            = 2,
	tfxGraphEasingType_ease_in_out_quad                         = 3,
	tfxGraphEasingType_ease_in_circular                         = 13,
	tfxGraphEasingType_ease_out_circular                        = 14,
	tfxGraphEasingType_ease_in_out_circular                     = 15,
} tfx_graph_easing_type;

typedef enum {
	tfxErrorCode_success                                        = 0,
	tfxErrorCode_incorrect_package_format                       = 1 << 0,
	tfxErrorCode_data_could_not_be_loaded                       = 1 << 1,
	tfxErrorCode_could_not_add_shape                            = 1 << 2,
	tfxErrorCode_error_loading_shapes                           = 1 << 3,
	tfxErrorCode_some_data_not_loaded                           = 1 << 4,
	tfxErrorCode_unable_to_open_file                            = 1 << 5,
	tfxErrorCode_unable_to_read_file                            = 1 << 6,
	tfxErrorCode_wrong_file_size                                = 1 << 7,
	tfxErrorCode_invalid_format                                 = 1 << 8,
	tfxErrorCode_no_inventory                                   = 1 << 9,
	tfxErrorCode_invalid_inventory                              = 1 << 10,
	tfxErrorCode_file_version_out_of_date                       = 1 << 11,
	tfxErrorCode_library_loaded_without_shape_loader            = 1 << 13,
	tfxErrorCode_library_object_could_not_be_created            = 1 << 14,
	tfxErrorCode_some_images_loaded_without_user_ptr            = 1 << 15
} tfx_error_flag_bits;

typedef tfxU32 tfxErrorFlags;                   //tfx_error_flag_bits

//-----------------------------------------------------------
//Section: forward_declarations
//-----------------------------------------------------------
#define tfxMAKE_HANDLE(handle) typedef struct handle##_s* handle;

//For allocating a new object with handle. Only used internally.
#define tfxNEW(type) (type)tfxALLOCATE(sizeof(type##_t))
#define tfxNEW_ALIGNED(type, alignment) (type)tfxALLOCATE_ALIGNED(sizeof(type##_t), alignment)

typedef struct tfx_package_s tfx_package_t;

tfxMAKE_HANDLE(tfx_package)
tfxMAKE_HANDLE(tfx_context);
tfxMAKE_HANDLE(tfx_library);
tfxMAKE_HANDLE(tfx_stage);
tfxMAKE_HANDLE(tfx_animation_manager);
tfxMAKE_HANDLE(tfx_effect_descriptor);
tfxMAKE_HANDLE(tfx_effect_template);
tfxMAKE_HANDLE(tfx_ribbon_buffer_requirements);
tfxMAKE_HANDLE(tfx_ribbon_dispatch);
tfxMAKE_HANDLE(tfx_gpu_shapes);

typedef struct tfx_image_data_s tfx_image_data_t;
typedef struct tfx_gpu_image_data_s tfx_gpu_image_data_t;
typedef struct tfx_bitmap_s tfx_bitmap_t;
typedef struct tfx_gpu_graph_data_s tfx_gpu_graph_data_t;
typedef struct tfx_instance_s tfx_instance_t;
typedef struct tfx_ribbon_bucket_s tfx_ribbon_bucket_t;
typedef struct tfx_sprite_data_s tfx_sprite_data_t;
typedef struct tfx_effect_index_s tfx_effect_index_t;
typedef struct tfx_effect_instance_data_s tfx_effect_instance_data_t;
typedef struct tfx_animation_instance_s tfx_animation_instance_t;
typedef struct tfx_frame_meta_s tfx_frame_meta_t;
typedef struct tfx_sprite_data_settings_s tfx_sprite_data_settings_t;
typedef struct tfx_emitter_path_s tfx_emitter_path_t;

typedef void *(*tfx_allocate_callback)(void *user_data, size_t size, size_t alignment);
typedef void (*tfx_deallocate_callback)(void *user_data, void *memory, size_t size, size_t alignment);

typedef struct tfx_allocation_callbacks_s {
	void *user_data;
	tfx_allocate_callback allocate;
	tfx_deallocate_callback deallocate;
} tfx_allocation_callbacks_t;

typedef struct tfx_random_s {
	tfxU64 seeds[2];
}tfx_random_t;

typedef struct tfx_version_s {
	const char *name;
	int major;
	int minor;
	int patch;
} tfx_version_t;

typedef struct tfx_ribbon_buffer_info_s {
	tfxU32 vertices_per_segment; 
	tfxU32 triangles_per_segment; 
	tfxU32 indices_per_segment;  
	tfxU32 total_segments; 
	tfxU32 index_count; 
	tfxKey pipeline_index;
} tfx_ribbon_buffer_info_t;

//This struct is used for configuring a effect manager on creation
typedef struct tfx_stage_info_s {
	double warmup_delta_time;				//The frame length tick amount for warming up effects. Higher is more performant at the cost of accuracy.
	tfxU32 max_particles;					//The maximum number of instance_data for each layer. This setting is not relevent if dynamic_sprite_allocation is set to true or group_sprites_by_effect is true.
	tfxU32 max_effects;                     //The maximum number of effects that can be updated at the same time.
	tfxU32 max_ribbon_segments;             //All segments for ribbons are stored in a single buffer. You will need to create buffers for rendering and so whatever you decide the max segments should be your buffers
											//should be big enough to contain all ribbon segments that you might need. You can call tfx_GetSegmentBufferSizeInBytes after creating the effect manager to get the byte
											//value that you can use to create the buffers. Also note that segments are always created in multiples of 32, so whatever number you put here it will be rounded to the
											//nearest multiple of 32.
	tfxU32 max_ribbons;						//The maximum number of ribbon instances that can exist at the same time across all emitters in this effect manager.
	tfxU32 ribbon_tessellation;				//The amount of tessellation used for ribbons. Currently this is set globally. 1 is generally enough for most cases.
	tfxU32 multi_threaded_batch_size;       //The size of each batch of particles to be processed when multithreading. Must be a power of 2 and 256 or greater.
	tfxU32 sort_passes;                     //when in order by depth mode (not guaranteed order) set the number of sort passes for more accuracy. Anything above 5 and you should just be guaranteed order.
	bool double_buffer_sprites;             //Set to true to double buffer instance_data so that you can interpolate between the old and new positions for smoother animations.
	bool dynamic_sprite_allocation;         //Set to true to automatically resize the sprite buffers if they run out of space. Not applicable when grouping instance_data by effect.
	bool group_sprites_by_effect;           //Set to true to group all instance_data by effect. Effects can then be drawn in specific orders or not drawn at all on an effect by effect basis.
	bool auto_order_effects;                //When group_sprites_by_effect is true then you can set this to true to sort the effects each frame. Use tfx_SetStageCamera in 3d to set the effect depth to the distance the camera.
	void *user_data;						//User data that will get passed into the grow_staging_buffer_callback function which you can use to grow the buffer
	//If you need the staging buffer to be grown dynamically then you can use this call back to do that. It should return true if the buffer was successfully grown or false otherwise.
	bool(*grow_staging_buffer_callback)(tfxU32 new_size, tfx_stage pm, void *user_data);
} tfx_stage_info_t;

typedef struct tfx_ribbon_dispatch_s {
	tfx_ribbon_bucket_t *ribbon_data;
	tfxU32 index_offset;
	tfxU32 vertex_offset;
	tfxU32 index_count;
	tfxU32 vertex_count;
	tfxU32 ribbon_offset;
	tfxU32 segment_offset;
	tfxU32 total_segments;
	tfxU32 last_index_offset;
	tfxU32 last_vertex_offset;
	tfxU32 last_ribbon_offset;
	tfxU32 last_segment_offset;
} tfx_ribbon_dispatch_t;

typedef struct tfx_ribbon_buffer_requirements_s {
	tfxU32 segment_buffer_size_in_bytes;
	tfxU32 ribbon_buffer_size_in_bytes;
	tfxU32 emitter_buffer_size_in_bytes;
} tfx_ribbon_buffer_requirements_t;

typedef struct tfx_animation_buffer_metrics_s {
	size_t sprite_data_size;
	tfxU32 offsets_size;
	tfxU32 instances_size;
	size_t offsets_size_in_bytes;
	size_t instances_size_in_bytes;
	tfxU32 total_sprites_to_draw;
	//Ribbon Metrics
	size_t ribbon_data_size;
	tfxU32 ribbon_offsets_size;
	tfxU32 ribbon_offsets_size_in_bytes;
	tfxU32 total_ribbons_to_draw;
	size_t ribbon_segment_data_size;
}tfx_animation_buffer_metrics_t;

//This can be sent as a push constant to the gpu
typedef struct tfx_ribbon_bucket_globals_s  {
	tfx_float32x4_t camera_position;
	tfxU32 segment_count;
	tfxU32 tessellation;  
	tfxU32 index_offset;
	tfxU32 vertex_offset;
	tfxU32 ribbon_count;
	tfxU32 ribbon_offset;
	tfxU32 segment_offset;
	tfxU32 uniform_index;
	tfxU32 emitters_index;
	tfxU32 graphs_index;
	tfxU32 ribbons_index;
	tfxU32 ribbon_segments_index;
	tfxU32 vertexes_index;
	tfxU32 indexes_index;
	tfxU32 image_data_index;
	tfxU32 sampler_index;
	tfxU32 particle_texture_index;
	tfxU32 color_ramp_texture_index;
	float lerp;
	float time;
	float ndc_offset_x;
	float ndc_offset_y;
} tfx_ribbon_bucket_globals_t;

typedef struct tfx_sprite_data_push_s {
	tfxU32 animation_instances_total;
	tfxU32 billboards_total;
	tfxU32 animated_shapes;	
	tfxU32 offsets_index;
	tfxU32 animation_instances_index;
	tfxU32 billboards_index;
	tfxU32 sprite_data_index;
	tfxU32 image_data_index;
	tfxU32 emitter_properties_index;
	tfxU32 bounding_boxes_index;
} tfx_sprite_data_push_t;

typedef struct tfx_gpu_graph_data_s {
	tfx_float32x4_t node_data;
	tfx_float32x4_t oscillator;
	tfx_graph_easing_type easing_type;
	int flags;
	int padding[2];
} tfx_gpu_graph_data_t;

//------------------------------------------------------------
//Section: Callback_typedefs
//------------------------------------------------------------
typedef void(*tfx_shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data);
typedef void(*tfx_uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset);
typedef bool(*tfx_maybe_render_instance_callback)(tfx_animation_manager animation_manager, tfx_float32x3_t position, float radius, void *user_data);


tfxAPI void tfx_UpdateAnimationManagerBufferMetrics(tfx_animation_manager animation_manager);
tfxAPI float tfx_DegreesToRadians(float degrees);
tfxAPI float tfx_RadiansToDegrees(float radians);

//--------------------------------
//Random numbers
//--------------------------------
tfxAPI tfx_random_t tfx_CreateRandom(tfxU32 seed);
tfxAPI void tfx_AdvanceRandom(tfx_random_t *random);
tfxAPI void tfx_RandomReseedTime(tfx_random_t *random);
tfxAPI void tfx_RandomReseed2(tfx_random_t *random, tfxU64 seed1, tfxU64 seed2);
tfxAPI void tfx_RandomReseed(tfx_random_t *random, tfxU64 seed);
tfxAPI float tfx_GenerateRandom(tfx_random_t *random);
tfxAPI float tfx_RandomRangeZeroToMax(tfx_random_t *random, float max);
tfxAPI float tfx_RandomRangeFromTo(tfx_random_t *random, float from, float to);
tfxAPI int tfx_RandomRangeFromToInt(tfx_random_t *random, int from, int to);
tfxAPI tfxU32 tfx_RandomRangeZeroToMaxUInt(tfx_random_t *random, tfxU32 max);
tfxAPI int tfx_RandomRangeZeroToMaxInt(tfx_random_t *random, int max);
tfxAPI void tfx_AlterRandomSeedU64(tfx_random_t *random, tfxU64 amount);
tfxAPI void tfx_AlterRandomSeedU32(tfx_random_t *random, tfxU32 amount);

//[API functions]
//All the functions below represent all that you will need to call to implement TimelineFX

/*
Quickstart
----------
The minimum lifecycle to get an effect on screen. Error handling and renderer specifics are omitted; every
tfx_* function used here is documented in full below.

    // 1. Initialise the library and its memory pool (once, at startup).
    tfx_BeginTimelineFX(tfx_GetDefaultThreadCount(), 128 * 1024 * 1024);

    // 2. Load an effects library from a .tfx file. shape_loader uploads each image to your renderer;
    //    uv_lookup returns its uv rect. Then build the GPU shape/graph/color-ramp data and upload it.
    tfx_library library = tfx_LoadEffectLibrary("effects.tfx", shape_loader, uv_lookup, user_data);

    // 3. Create a stage (the runtime effect manager) to simulate effects in.
    tfx_stage stage = tfx_CreateStage(tfx_CreateStageInfo(tfxStageSetup_none));

    // 4. Create a template for the effect you want, then add an instance of it to the stage. The returned
    //    tfxEffectID is your handle for manipulating that live instance (position, scale, expiry, ...).
    tfx_effect_template explosion = tfx_CreateEffectTemplate(library, "Explosion");
    tfxEffectID effect_id = tfx_AddEffectTemplateToStage(stage, explosion);

    // 5. Per-frame update loop, ideally in a fixed timestep (see tfx_UpdateStage for the full pattern).
    tfx_UpdateStage(stage, 16.66667);   // elapsed ms since last update

    // 6. Read the results for rendering. These getters complete the stage's in-flight work for you
    //    (see the threading contract in the Particle_Manager section). Copy the instance buffer to your
    //    GPU staging buffer and draw; drive ribbons with the tfx_*RibbonDispatch iterator.
    tfxU32 count = tfx_GetInstanceCount(stage);
    tfx_instance_t *instances = tfx_GetInstanceBuffer(stage);
    // ... upload `count` instances and issue your draw call ...

    // 7. Shutdown (see "Ownership and teardown ordering" below).
    tfx_FreeEffectTemplate(explosion);
    tfx_FreeStage(stage);
    tfx_FreeLibrary(library);
    tfx_EndTimelineFX();

Ownership and teardown ordering
-------------------------------
Everything is allocated from the memory pool created by tfx_BeginTimelineFX, so tfx_EndTimelineFX must be the
very last TimelineFX call you make — it tears down the worker threads and the pool itself.

tfx_EndTimelineFX sweeps up any stages, animation managers and libraries you did not free yourself (it completes
and frees every stage, then frees every animation manager, then every library). Freeing any of them yourself
first is fine — they deregister from the sweep as they go, so there is no double free either way. It does NOT
sweep effect templates, so you must free those explicitly with tfx_FreeEffectTemplate before calling
tfx_EndTimelineFX or the built-in leak check will report leaked blocks.

If you free objects yourself rather than leaving them to the sweep, free them in dependency order: free/complete
a stage before the library whose effects it is simulating (a live stage references the library's effect data);
an effect template holds its own detached clone of the effect so it is independent of the library once created,
but freeing templates before the library keeps the ordering simple. A stage has outstanding worker-thread work
after tfx_UpdateStage, so call tfx_CompleteStageWork (or use tfx_FreeStage, which completes it for you) before
tearing anything down mid-frame.
*/

//--------------------------------
//Initialisation_functions
//--------------------------------

/*
You don't have to call this, you can just call tfx_BeginTimelineFX in order to initialise the memory, but I created this for the sake of the editor which
needs to load in an ini file before initialising timelinefx which requires the memory pool to be created before hand. Passing NULL allocation_callbacks 
uses the default internal allocator and deallocator.
* @param memory_pool_size    The size of each memory pool to contain all objects created in TimelineFX, recommended to be at least 64MB
* @param allocation_callbacks    The callbacks used to allocate and deallocate every TimelineFX backing pool, or NULL for the default internal allocator and deallocator
*/
#ifdef __cplusplus
tfxAPI tfx_context tfx_InitialiseTimelineFXMemory(size_t memory_pool_size, const tfx_allocation_callbacks_t *allocation_callbacks = nullptr);
#else
tfxAPI tfx_context tfx_InitialiseTimelineFXMemory(size_t memory_pool_size, const tfx_allocation_callbacks_t *allocation_callbacks);
#endif

/*
Initialise TimelineFX. Must be called before any functionality of TimelineFX is used.
* @param max_threads        The number of worker threads to use in addition to the main thread. Pass 0 to run single threaded.
*                            The count is clamped to the number of hardware threads available on the machine.
* @param memory_pool_size    The size of each memory pool to contain all objects created in TimelineFX, recommended to be at least 64MB
* @param allocation_callbacks    The callbacks used to allocate and deallocate every TimelineFX backing pool, or NULL for the default internal allocator and deallocator
*/
#ifdef __cplusplus
tfxAPI tfx_context tfx_BeginTimelineFX(int max_threads, size_t memory_pool_size, const tfx_allocation_callbacks_t *allocation_callbacks = nullptr);
#else
tfxAPI tfx_context tfx_BeginTimelineFX(int max_threads, size_t memory_pool_size, const tfx_allocation_callbacks_t *allocation_callbacks);
#endif

/*
Cleanup up all threads and memory used by timelinefx
*/
tfxAPI void tfx_EndTimelineFX(void);

/*
Get the context used internally for bookkeeping. If allocation callbacks are passed during initialisation, they allocate the context.
*/
tfxAPI tfx_context tfx_GetContext(void);

/*
Adopt a context obtained from tfx_GetContext into the currently running module and rebind every function pointer that
the module owns. Call this after a hot-reload of the dynamic library hosting the TimelineFX instance, before calling
tfx_ResumeTimelineFX.
Passing non-NULL allocation_callbacks replaces the stored callbacks. Passing NULL preserves the stored callbacks.
The caller must provide fresh callback addresses when the prior callback addresses are expected to be invalid.
Any update or uv lookup callbacks previously registered on effects, libraries and animation managers are cleared and
must be registered again by the caller.
* @returns    True if the context was adopted
*/
#ifdef __cplusplus
tfxAPI bool tfx_SetContext(tfx_context context, const tfx_allocation_callbacks_t *allocation_callbacks = nullptr);
#else
tfxAPI bool tfx_SetContext(tfx_context context, const tfx_allocation_callbacks_t *allocation_callbacks);
#endif

/*
Drain and stop all current TimelineFX asynchronous work without deallocating memory.
*/
tfxAPI void tfx_SuspendTimelineFX(void);

/*
Restart the worker threads after a call to tfx_SuspendTimelineFX.
*/
tfxAPI void tfx_ResumeTimelineFX(void);

//--------------------------------
//Global_variable_access
//--------------------------------
/*
Set the color format used for storing color ramps. Color ramps are generated by particle emitters and dictate how the particle colors change over the lifetime
of the particle. They can be uploaded to the GPU so you can set your preference for the color format as you need. The format should be immediately after you
call tfx_BeginTimelineFX
* @param color_format        The tfx_color_format to store color ramps in: tfx_color_format_rgba8, tfx_color_format_rgba16f or tfx_color_format_rgba32f
*/
tfxAPI void tfx_SetColorRampFormat(tfx_color_format color_format);

//--------------------------------
//Library_functions
//--------------------------------

/*
Create a new, empty library handle. Normally you will use tfx_LoadEffectLibrary to create a library
already populated from a tfx file; use this only if you want to build a library up from scratch.
Free it with tfx_FreeLibrary.
* @returns tfx_library    A handle to a new, empty library.
*/
tfxAPI tfx_library tfx_CreateLibrary(void);

/*
Count the number of particle shapes stored in a tfx file without loading the whole library. Useful
if you need to reserve storage for the shapes up front before loading.
* @param filename        The name of the tfx file to inspect
* @returns int            The number of shapes in the file.
*/
tfxAPI int tfx_GetShapeCountInLibrary(const char *filename);

/*
Validate a timelinefx tfx file to make sure that it's valid.
* @param filename        The name of the file where you want to count the number of shapes
* @returns int            Returns 0 if the file successfully validated or a tfxErrorFlags if something went wrong
*/
tfxAPI int tfx_ValidateEffectPackage(const char *filename);

/**
* Loads an effect library package from the specified filename into the provided tfx_library_t object.
*
* @param filename         A pointer to a null-terminated string that contains the path and filename of the effect library package to be loaded.
* @param shape_loader     A pointer to a function that will be used to load image data into the effect library package.
*                         The function has the following signature: void shape_loader(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data).
* @param user_data        A pointer to user-defined data that will be passed to the shape_loader function. This parameter is optional and can be set to nullptr if not needed.
* @return tfx_library	  A handle to a library object
*/
tfxAPI tfx_library tfx_LoadEffectLibrary(const char *filename, tfx_shape_loader shape_loader, tfx_uv_lookup uv_lookup, void *user_data);

/**
* Loads an effect library package from memory into the provided tfx_library_t object pointer.
*
* @param data             A pointer to a memory buffer containing the library to be loaded
* @param size             The size of the memory buffer containing the library to be loaded
* @param shape_loader     A pointer to a function that will be used to load image data into the effect library package.
*                         The function has the following signature: void shape_loader(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data).
* @param user_data        A pointer to user-defined data that will be passed to the shape_loader function. This parameter is optional and can be set to nullptr if not needed.
* @return tfx_library	  A handle to a library object
*/
tfxAPI tfx_library tfx_LoadEffectLibraryFromMemory(const void *data, tfxU32 size, tfx_shape_loader shape_loader, tfx_uv_lookup uv_lookup, void *user_data);

/*
Get the error flags from a library. When you load a library from file or memory, if something goes wrong then the error status is stored in the library object and you can retrieve it with this command.
* @param lib            A handle to a tfx_library object that will hold the loaded effect library data.
* @return A tfxErrorFlags value that indicates whether the function succeeded or failed. The possible return values are:
	tfxErrorCode_success = 0
	tfxErrorCode_incorrect_package_format
	tfxErrorCode_data_could_not_be_loaded
	tfxErrorCode_could_not_add_shape
	tfxErrorCode_error_loading_shapes
	tfxErrorCode_some_data_not_loaded
	tfxErrorCode_unable_to_open_file
	tfxErrorCode_unable_to_read_file
	tfxErrorCode_wrong_file_size
	tfxErrorCode_invalid_format
	tfxErrorCode_no_inventory
	tfxErrorCode_invalid_inventory
*/
tfxAPI tfxErrorFlags tfx_GetLibraryErrorStatus(tfx_library library);

/**
* Loads a sprite data file into an animation manager
*
* @param filename        A pointer to a null-terminated string that contains the path and filename of the effect library package to be loaded.
* @param lib            A reference to a tfx_animation_manager_t object that will hold the loaded sprite data.
* @param shape_loader    A pointer to a function that will be used to load image data into the effect library package.
*                        The function has the following signature: void shape_loader(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data).
* @param user_data        A pointer to user-defined data that will be passed to the shape_loader function. This parameter is optional and can be set to nullptr if not needed.
*
* @return A tfxErrorFlags value that indicates whether the function succeeded or failed. The possible return values are:
	tfxErrorCode_success = 0
	tfxErrorCode_incorrect_package_format
	tfxErrorCode_data_could_not_be_loaded
	tfxErrorCode_could_not_add_shape
	tfxErrorCode_error_loading_shapes
	tfxErrorCode_some_data_not_loaded
	tfxErrorCode_unable_to_open_file
	tfxErrorCode_unable_to_read_file
	tfxErrorCode_wrong_file_size
	tfxErrorCode_invalid_format
	tfxErrorCode_no_inventory
	tfxErrorCode_invalid_inventory
*/
tfxAPI tfxErrorFlags tfx_LoadSpriteData(const char *filename, tfx_animation_manager animation_manager, tfx_shape_loader shape_loader, void *user_data);

/*
* Updates all the image data in the library using the uv_lookup that you set when loading a library. This allows you to add all of the uv data for
* the shapes that are loaded into the texture. You must have set the uv_lookup callback when loading the library, otherwise you can loop over the 
* shapes in the library and update the data yourself using the tfx_GetLibraryShapeArray and related functions.
* @param tfx_library                A valid pointer to a tfx_library
*/
tfxAPI void tfx_UpdateLibraryGPUImageData(tfx_library library);

/*
Get the number of shapes stored in the library
* @param tfx_library                A valid pointer to a tfx_library
* @return tfxU32					Count of shapes
*/
tfxAPI tfxU32 tfx_GetLibraryImageCount(tfx_library library);

/*
Get a particle image from a library by it's index
* @param tfx_library                A valid pointer to a tfx_library
* @return image						A tfx_image_data_t object with all the details of the image
*/
tfxAPI tfx_image_data_t *tfx_GetLibraryImage(tfx_library library, tfxU32 index);

tfxAPI void tfx_SetImagePointer(tfx_image_data_t *image, void *pointer);
tfxAPI void tfx_SetGPUImageTextureInfo(tfx_gpu_image_data_t *image, float x, float y, float z, float w, int array_index);
tfxAPI int tfx_GetImageFrameCount(tfx_image_data_t *image);
tfxAPI int tfx_GetImageWidth(tfx_image_data_t *image);
tfxAPI int tfx_GetImageHeight(tfx_image_data_t *image);
tfxAPI void* tfx_GetBitmapData(tfx_bitmap_t *bitmap);
tfxAPI size_t tfx_GetBitmapSize(tfx_bitmap_t *bitmap);
tfxAPI int tfx_GetBitmapWidth(tfx_bitmap_t *bitmap);
tfxAPI int tfx_GetBitmapHeight(tfx_bitmap_t *bitmap);

/*
Output all the effect names in a library to the console
* @param tfx_library                A valid pointer to a tfx_library
*/
tfxAPI void tfx_ListEffectNames(tfx_library library);

/*
Get an effect in the library by it's index. If you need to get an effect in a folder or an emitter then you can use tfx_GetLibraryEffectPath instead.
* @param tfx_library                A valid pointer to a tfx_library
*/
tfxAPI tfx_effect_descriptor tfx_GetEffectByIndex(tfx_library library, int index);

/*
Get an effect in the library by it's path. So for example, if you want to get a pointer to the emitter "spark" in effect "explosion" then you could do GetEffect("explosion/spark")
You will need this function to apply user data and update callbacks to effects and emitters before adding the effect to the effect manager
* @param tfx_library_t                A valid pointer to a tfx_library_t
* @param const char *path             Path to the effect or emitter
*/
tfxAPI tfx_effect_descriptor tfx_GetLibraryEffectPath(tfx_library library, const char *path);

/*
Check whether a path resolves to an effect or emitter in the library. tfx_GetLibraryEffectPath asserts
in debug builds and returns NULL in release builds when a path is missing, so use this to probe a path
safely before looking it up.
* @param library                A valid pointer to a tfx_library_t
* @param path                   Path to the effect or emitter
* @returns                      true if the path exists in the library
*/
tfxAPI bool tfx_IsValidEffectPath(tfx_library library, const char *path);

/*
Free all the memory used by a library
* param tfx_library				A pointer to the library that you want to free
*/
tfxAPI void tfx_FreeLibrary(tfx_library library);

/*
Create the image data required for shaders from a TimelineFX library. The image data will contain data such as uv coordinates. Once you have built the data you can use GetLibraryImageData to get the buffer
and upload it to the gpu.
* @param library                  A pointer to a tfx_library_t object
* @param shapes                   A pointer to a tfx_gpu_shapes_t object which will fill a buffer with all the shapes
* @param uv_lookup                A function pointer to a function that you need to set up in order to get the uv coordinates from whatever renderer you're using
*/
tfxAPI void tfx_BuildLibraryGPUShapeData(tfx_library library, tfx_gpu_shapes shapes, tfx_uv_lookup uv_lookup);

/*
Get a pointer to the particle shapes data in a library. This can be used with tfx_BuildGPUShapeData when you want to upload the data to the GPU
* @param library        A pointer to a tfx_library_t
* @param count			A pointer to an int that will be filled with the nubmer of images in the image data array that's returned
*/
tfxAPI tfx_image_data_t *tfx_GetParticleShapesLibrary(tfx_library library, int *count);

/*
Get a count of the number of color ramp bitmaps in the library. Color ramps are used to change the color of particles over time and you will need to upload them to the GPU.
* @param library        A pointer to a tfx_library_t
*/
tfxAPI tfxU32 tfx_GetColorRampBitmapCount(tfx_library library);

/*
Get a pointer to a color ramp bitmap in a library. You can use this data to upload the bitmaps to the GPU.
* @param library        A pointer to a tfx_library_t
*/
tfxAPI tfx_bitmap_t *tfx_GetColorRampBitmap(tfx_library library, tfxU32 index);

/*
Get a count of the number of color ramp bitmaps in an animation manager.
* @param animation_manager        A handle to a tfx_animation_manager
*/
tfxAPI tfxU32 tfx_GetAnimationColorRampBitmapCount(tfx_animation_manager animation_manager);

/*
Get a pointer to a color ramp bitmap in an animation manager. You can use this data to upload the bitmaps to the GPU.
* @param animation_manager        A handle to a tfx_animation_manager
* @param index                    The index of the bitmap
*/
tfxAPI tfx_bitmap_t *tfx_GetAnimationColorRampBitmap(tfx_animation_manager animation_manager, tfxU32 index);

/*
Check to see if a library has been initialised or not
* @param library        A pointer to a tfx_library_t
*/
tfxAPI bool tfx_LibraryIsInitialised(tfx_library library);

/*
Get the gpu shapes handle in library. The gpu shapes handle can be used to upload the image data for particle shapes to the gpu using functions like tfx_GetGPUShapesArray, tfx_GetGPUShapesSizeInBytes, 
tfx_GetGPUShapesCount etc.
* @param library        A handle to a tfx_library
*/
tfxAPI tfx_gpu_shapes tfx_GetLibraryGPUShapes(tfx_library library);

/*
Get a pointer to the global gpu graph lookup data. This is a single shared buffer maintained across all
libraries in tfxStore, so it takes no library argument. Shaders use this lookup data to update attributes of
particles and ribbons (currently ribbons only). Use with tfx_GetGPUGraphLookupsBufferSizeInBytes to upload it
to your GPU buffer.
* @returns tfx_gpu_graph_data_t *   A pointer to the shared gpu graph lookup data.
*/
tfxAPI tfx_gpu_graph_data_t *tfx_GetGPUGraphLookupsBuffer(void);

/*
Get the size in bytes of the global gpu graph lookup data returned by tfx_GetGPUGraphLookupsBuffer. This is the
single shared buffer in tfxStore, so it takes no library argument.
* @returns tfxU32                   The size of the shared gpu graph lookup data in bytes.
*/
tfxAPI tfxU32 tfx_GetGPUGraphLookupsBufferSizeInBytes(void);

/*
Get buffer info for ribbons based on the tessellation value. Returns a tfx_ribbon_buffer_info_t object with
vertices per segment, triangles per segment and indices per segment requirements.
* @param tfxU32        					The number of tessellations for the ribbons
* @returns tfx_ribbon_buffer_info_t		Info containing vertices, triangles and indices per segment
*/
tfxAPI tfx_ribbon_buffer_info_t tfx_GenerateRibbonBufferInfo(tfxU32 tessellation);

//--------------------------------
//Particle_Manager_functions
//--------------------------------

/*
Threading contract
------------------
A stage updates particles across the worker threads created by tfx_BeginTimelineFX. tfx_UpdateStage kicks off
that work and returns immediately while the work is still in flight on the worker threads. This means that after
tfx_UpdateStage returns, the sprite/instance and ribbon buffers are NOT yet safe to read.

Two ways to synchronise before you read that data:

  - Call tfx_CompleteStageWork(pm) explicitly. This blocks until all of the stage's outstanding update work has
    finished, after which every buffer is safe to read.
  - Or rely on the getters that synchronise implicitly. The following complete the stage's outstanding work for
    you before returning, so you can call them directly after tfx_UpdateStage without calling
    tfx_CompleteStageWork first:
        tfx_GetInstanceBuffer, tfx_GetInstanceBufferByLayer, tfx_GetInstanceCount, tfx_GetInstanceCountByLayer,
        tfx_GetNextInstanceBuffer, tfx_GetRibbonBuffers, tfx_HasRibbonsToDraw, tfx_NextRibbonDispatch,
        tfx_GetRibbonBufferRequirements, tfx_CopyRibbonDataToStagingBuffers, tfx_ClearStage and tfx_FreeStage.

Anything that reads stage state through a path NOT in that list (for example reading the instance buffer pointer
you cached last frame, or touching the raw buffers directly) must be preceded by tfx_CompleteStageWork, otherwise
you will race the worker threads. Mutators such as tfx_SetEffectPosition are also only safe once the update work
for the frame has completed so it's a good idea to set positions before the call to tfx_UpdateStage. 
*/

/*
Create a tfx_stage_info_t object which contains configuration data that you can pass to tfx_CreateStage to setup a effect manager. You can tweak the config after calling this
function if needed to fine tune the settings.
* @param setup                    A tfx_stage_setup enum which you can use to set the info based on some commonly used templates
*/
tfxAPI tfx_stage_info_t tfx_CreateStageInfo(tfx_stage_setup setup);

/*
Create a stage (effect manager) configured with a tfx_stage_info_t object. Build the info with tfx_CreateStageInfo
and tweak it as needed before passing it in. If you later want to run the stage in a different mode then call
tfx_ReconfigureStage rather than creating a new one.
* @param info                   A tfx_stage_info_t containing the configuration for the stage.
* @returns tfx_stage            A handle to the new stage, or NULL if the allocation failed (out of memory).
*/
tfxAPI tfx_stage tfx_CreateStage(tfx_stage_info_t info);

/*
Reconfigure a effect manager to make it work in a different mode. A effect manager can only run in a single mode at time like unordered, depth ordered etc so use this to change that. Also bear
in mind that you can just use more than one effect manager and utilised different modes that way as well. The modes that you need will depend on the effects that you're adding to the effect manager.
* @param pm                       A pointer to an intialised tfx_stage_t.
* @param sort_passes              The number of sort passes if you're using depth sorted effects
*/
tfxAPI void tfx_ReconfigureStage(tfx_stage pm, tfxU32 sort_passes);

/*
Block until all of the stage's in-flight update work (started by tfx_UpdateStage) has finished. Call this before
reading the stage's sprite/instance or ribbon buffers, or mutating effects, unless you are only using the getters
that synchronise implicitly (see the threading contract above). Safe to call redundantly; it is a no-op if there
is no work outstanding.
* @param pm                       A handle to an initialised tfx_stage_t.
*/
tfxAPI void tfx_CompleteStageWork(tfx_stage pm);

/*
Turn on and off whether the effect manager should sort the effects by depth order. Use tfx_SetStageCamera to set the position of the camera that the effect manager will
use to update the depth of each effect in the scene.
* @param pm                       A pointer to an intialised tfx_stage_t.
* @param yesno                    A boolean, set to true or false if you want auto ordering on or off respectively
*/
tfxAPI void tfx_ToggleStageOrderEffects(tfx_stage pm, bool yesno);

/*
Get the billboard buffer in the effect manager containing all the sprite instances that were created in the most recent frame. You can use this to copy to a staging buffer to upload to the gpu.
* @param pm                       A pointer to an intialised tfx_stage_t.
*/
tfxAPI tfx_instance_t *tfx_GetInstanceBuffer(tfx_stage  pm);

/*
Get the billboard buffer in the effect manager containing all the sprite instances for a specific layer that were created in the most recent frame. You can use this to copy to a staging buffer to upload to the gpu.
You can then use tfx_GetInstanceCountByLayer for the draw call.
* @param pm                       A pointer to an intialised tfx_stage_t.
*/
tfxAPI tfx_instance_t *tfx_GetInstanceBufferByLayer(tfx_stage pm, tfxU32 layer);

/*
Get the number of instances within the instance buffer of a effect manager
* @param pm                       A pointer to an intialised tfx_stage_t.
*/
tfxAPI int tfx_GetInstanceCount(tfx_stage pm);

/*
Get the number of instances within the instance buffer of a effect manager for a specific layer.
* @param pm                       A pointer to an intialised tfx_stage_t.
*/
tfxAPI int tfx_GetInstanceCountByLayer(tfx_stage pm, tfxU32 layer);

/*
Get the update time being used by the effect manager.
* @param pm                       A handle to an intialised tfx_stage_t.
*/
tfxAPI double tfx_GetUpdateTime(tfx_stage pm);

/*
Get the ribbon buffer for a given segment size. This will give you all the necessary info and buffer pointers for uploading the ribbon data to the GPU for processing and converting into
a vertex buffer for rendering.
* @param pm                       A pointer to an intialised tfx_stage_t.
* @param segment_count            An unsigned int specifying the ribbon length that you want the rendering info for.
* @returns						  A pointer to a tfx_ribbon_bucket_t
*/
tfxAPI tfx_ribbon_bucket_t *tfx_GetRibbonBuffers(tfx_stage pm, tfxKey bucket_id);

/*
Call this to determine whether or not any effect manager has ribbon_emitters to draw this frame.
* @returns						  True or false
*/
tfxAPI bool tfx_HasRibbonsToDraw(tfx_stage pm);

/*
Get a struct containing the info you need to compute and render ribbon_emitters of a specific length. 
* @param pm                       A pointer to an intialised tfx_stage_t.
* @param segment_count            An unsigned int specifying the ribbon length that you want the rendering info for.
* @returns						  A tfx_ribbon_buffer_info_t struct
*/
tfxAPI tfx_ribbon_buffer_info_t tfx_GetRibbonBufferInfo(tfx_stage pm, tfxKey bucket_id);

/*
--------------------------------
Ribbon rendering / GPU dispatch
--------------------------------
Ribbons are turned into geometry on the GPU by a compute shader and then drawn. Unlike billboard
sprites, which you upload wholesale from a single instance buffer, ribbons are batched into "buckets" (one per
distinct segment length / tessellation) and you drive the upload and dispatch by iterating those buckets.

Per-frame flow, after tfx_UpdateStage:

  1. Size and (re)allocate your GPU buffers. Use tfx_GetRibbonBufferRequirements() for the exact bytes needed
     this frame across every stage, or the tfx_Get*MaxSizeInBytes / tfx_GetTotal*MaxSizeInBytes helpers to size
     for the worst case up front so you never resize.
  2. Copy the CPU ribbon data into your mapped staging buffers with tfx_CopyRibbonDataToStagingBuffers, then
     transfer them to your device-local segment / ribbon / emitter buffers.
  3. Iterate the buckets to record compute dispatches: create a tfx_ribbon_dispatch_t with
     tfx_CreateRibbonDispatch(), then loop while tfx_NextRibbonDispatch(pm, &dispatch) returns true. Each call
     fills in the per-bucket offsets/counts and populates dispatch.ribbon_data->globals with the push-constant
     block; dispatch total_segments tells you how many threads to launch.
  4. Iterate the buckets again the same way to record the draw calls, using dispatch.index_count and
     dispatch.index_offset for the indexed draw.

The dispatch iterator advances an internal cursor on the stage, so if you need to walk the buckets more than
once in a frame (compute pass, then render pass) call tfx_ResetRibbonDispatchIterator(pm) between passes.
tfx_NextRibbonDispatch / tfx_GetRibbonBufferRequirements / tfx_CopyRibbonDataToStagingBuffers all implicitly
complete any in-flight stage update work first (see the threading contract above), so they are safe to call
straight after tfx_UpdateStage. See the RibbonComputeFunction / RenderRibbons functions in the editor for a
complete worked example.
*/

/*
Create a zero-initialised tfx_ribbon_dispatch_t to drive the ribbon bucket iteration. Pass its address to
tfx_NextRibbonDispatch. The struct carries the running offsets between buckets, so create a fresh one at the
start of each iteration pass (or reset the fields to zero).
* @returns tfx_ribbon_dispatch_t   A zeroed dispatch struct ready to pass to tfx_NextRibbonDispatch.
*/
tfxAPI tfx_ribbon_dispatch_t tfx_CreateRibbonDispatch(void);

/*
Advance to the next ribbon bucket that has ribbons to draw. Call in a while loop until it returns false. On each
true return, ribbon_dispatch is filled with this bucket's vertex/index offsets and counts, total_segments (the
work item count for the compute dispatch) and a pointer to the bucket in ribbon_data whose globals member is the
push-constant block for both the compute and render passes. Implicitly completes any in-flight stage update.
* @param pm                        A handle to an initialised tfx_stage_t.
* @param ribbon_dispatch           A pointer to a tfx_ribbon_dispatch_t (from tfx_CreateRibbonDispatch) updated in place.
* @returns bool                    true while there is another bucket to process, false when iteration is finished.
*/
tfxAPI bool tfx_NextRibbonDispatch(tfx_stage pm, tfx_ribbon_dispatch_t *ribbon_dispatch);

/*
Reset the stage's internal ribbon bucket iterator back to the start. Call this between two passes over the
buckets in the same frame (for example after the compute pass, before the render pass) so tfx_NextRibbonDispatch
walks them again from the beginning.
* @param pm                        A handle to an initialised tfx_stage_t.
*/
tfxAPI void tfx_ResetRibbonDispatchIterator(tfx_stage pm);

tfxAPI tfx_ribbon_bucket_globals_t *tfx_GetRibbonDispatchGlobals(tfx_ribbon_dispatch_t *ribbon_dispatch);

/*
Get the exact ribbon buffer sizes needed this frame, summed across every registered stage, for allocating a
single shared set of GPU buffers. Implicitly completes any in-flight stage update on each stage. Use this when
you would rather allocate to the exact per-frame need than to the worst case (see tfx_GetTotal*MaxSizeInBytes).
* @returns tfx_ribbon_buffer_requirements_t   Segment, ribbon and emitter buffer sizes in bytes for this frame.
*/
tfxAPI tfx_ribbon_buffer_requirements_t tfx_GetRibbonBufferRequirements(void);

/*
Copy this stage's ribbon segment, ribbon instance and emitter data into your mapped staging buffers, ready to
transfer to the GPU. Size the destinations with tfx_GetRibbonBufferRequirements or the max-size helpers.
Implicitly completes any in-flight stage update.
* @param stage                     A handle to an initialised tfx_stage_t.
* @param segments_dst              Destination for the ribbon segment data (must be non-null).
* @param ribbons_dst               Destination for the ribbon instance data (must be non-null).
* @param emitters_dst              Destination for the ribbon emitter data (must be non-null).
*/
tfxAPI void tfx_CopyRibbonDataToStagingBuffers(tfx_stage stage, void *segments_dst, void *ribbons_dst, void *emitters_dst);

/*
Worst-case size in bytes of the ribbon segment buffer for a single stage, based on its configured max_ribbon_segments.
Use to size a buffer once at startup instead of resizing per frame.
* @param pm                        A handle to an initialised tfx_stage_t.
*/
tfxAPI size_t tfx_GetSegmentBufferMaxSizeInBytes(tfx_stage pm);

/*
Worst-case size in bytes of the ribbon vertex buffer for a single stage. Pass the size of your own vertex struct,
or 0 to use the library's tfx_ribbon_vertex_t size.
* @param pm                        A handle to an initialised tfx_stage_t.
* @param vertex_size               The size of a single vertex in bytes, or 0 for the default tfx_ribbon_vertex_t.
*/
tfxAPI size_t tfx_GetSegmentVertexBufferMaxSizeInBytes(tfx_stage pm, tfxU32 vertex_size);

/*
Worst-case size in bytes of the ribbon index buffer for a single stage, based on its configured max_ribbon_segments
and tessellation.
* @param pm                        A handle to an initialised tfx_stage_t.
*/
tfxAPI size_t tfx_GetSegmentIndexBufferMaxSizeInBytes(tfx_stage pm);

/*
Worst-case size in bytes of the ribbon instance buffer for a single stage, based on its configured max_ribbons.
* @param pm                        A handle to an initialised tfx_stage_t.
*/
tfxAPI size_t tfx_GetRibbonBufferMaxSizeInBytes(tfx_stage pm);

/*
Worst-case size in bytes of the ribbon emitter buffer for a single stage.
* @param pm                        A handle to an initialised tfx_stage_t.
*/
tfxAPI size_t tfx_GetEmitterBufferMaxSizeInBytes(tfx_stage pm);

/*
Get the size in bytes of the per-emitter particle properties buffer for a library. These properties are looked up
by the particle and ribbon shaders; use with tfx_GetParticlePropertiesBuffer to upload them to the GPU.
* @param library                  A handle to a tfx_library.
*/
tfxAPI size_t tfx_GetParticlePropertiesBufferSizeInBytes(tfx_library library);

/*
Get a pointer to the per-emitter particle properties buffer for a library, for uploading to the GPU. Size it with
tfx_GetParticlePropertiesBufferSizeInBytes.
* @param library                  A handle to a tfx_library.
*/
tfxAPI void *tfx_GetParticlePropertiesBuffer(tfx_library library);

/*
Get the total worst-case buffer sizes across all registered stages, for creating a single shared set of ribbon GPU
buffers. These iterate every registered stage and sum the individual tfx_Get*MaxSizeInBytes requirements. Use these
to allocate for the worst case once; use tfx_GetRibbonBufferRequirements instead to size to the exact per-frame need.
*/
tfxAPI size_t tfx_GetTotalSegmentBufferMaxSizeInBytes(void);

tfxAPI size_t tfx_GetTotalSegmentVertexBufferMaxSizeInBytes(tfxU32 vertex_size);

tfxAPI size_t tfx_GetTotalSegmentIndexBufferMaxSizeInBytes(void);

tfxAPI size_t tfx_GetTotalRibbonBufferMaxSizeInBytes(void);

tfxAPI size_t tfx_GetTotalEmitterBufferMaxSizeInBytes(void);

/*
When a effect manager updates particles it creates work queues to handle the work. By default these each have a maximum amount of 1000 entries which should be
more than enough for most situations. However you can increase the sizes here if needed. You only need to set this manually if you hit one of the asserts when these
run out of space or you anticipate a huge amount of emitters and particles to be used (> million). On the other hand, you might be tight on memory in which case you
could reduce the numbers as well if needed (they don't take a lot of space though)
* @param pm                        A pointer to an intialised tfx_stage_t.
* @param spawn_work_max            The maximum amount of spawn work entries
* @param control_work_max          The maximum amount of control work entries
* @param age_work_max              The maximum amount of age_work work entries
*/
tfxAPI void tfx_SetStageWorkQueueSizes(tfx_stage pm, tfxU32 spawn_work_max, tfxU32 control_work_max, tfxU32 age_work_max);

/*Free the memory for a specific emitter type. When an emitter is created it creates memory to store all of the particles that it updates each frame. If you have
multiple emitters of the same type then their particle lists are resused rather then freed as they expire. When they're freed then the unused list is added to a list
of free particle banks for that emitter type so that they can then be recycled if another emitter of the same type is created. If you want to free the memory for a
specific emitter then you can call this function to do that.
NOTE: No emitters of the type passed to the function must be in use in the effect manager.
* @param pm                        A pointer to an intialised tfx_stage_t.
* @param emitter                   A pointer to a valid tfx_effect_descriptor_t of type tfxEmitterType
*/
tfxAPI void tfx_FreeParticleListsMemory(tfx_stage pm, tfx_effect_descriptor emitter);

/*
Free all the memory that is associated with an effect. Depending on the configuration of the effect manager this might be instance_data, particle lists and spawn location lists.
* @param pm                        A pointer to an intialised tfx_stage_t.
* @param emitter                   A pointer to a valid tfx_effect_descriptor_t of type tfxEffectType
*/
tfxAPI void tfx_FreeEffectListsMemory(tfx_stage pm, tfx_effect_descriptor effect);

/*
Get the current particle count for a effect manager
* @param pm                        A pointer to an tfx_stage_t
* @returns tfxU32                  The total number of particles currently being updated
*/
tfxAPI tfxU32 tfx_GetParticleCount(tfx_stage pm);

/*
Get the current ribbon count for a effect manager
* @param pm                        A pointer to an tfx_stage_t
* @returns tfxU32                  The total number of particles currently being updated
*/
tfxAPI tfxU32 tfx_GetRibbonCount(tfx_stage pm);

/*
Get the current number of effects that are currently being updated by a effect manager
* @param pm                        A pointer to an tfx_stage_t
* @returns tfxU32                  The total number of effects currently being updated
*/
tfxAPI tfxU32 tfx_GetEffectCount(tfx_stage pm);

/*
Get the current number of emitters that are currently being updated by a effect manager
* @param pm                        A pointer to an tfx_stage_t
* @returns tfxU32                  The total number of emitters currently being updated
*/
tfxAPI tfxU32 tfx_GetEmitterCount(tfx_stage pm);

/*
Get the current number of free effects in the stage that are ready to be re-used
* @param pm                        A pointer to an tfx_stage_t
* @returns tfxU32                  The total number of free effects in the stage
*/
tfxAPI tfxU32 tfx_GetFreeEffectCount(tfx_stage pm);

/*
Set the seed for the effect manager for random number generation. Setting the seed can determine how an emitters spawns particles, so if you set the seed before adding an effect to the effect manager
then the effect will look the same each time. Note that seed of 0 is invalid, it must be 1 or greater.
* @param pm                        A pointer to an initialised tfx_stage_t. 
* @param seed                      An unsigned int representing the seed (Any value other then 0)
*/
tfxAPI void tfx_SetStageSeed(tfx_stage pm, tfxU64 seed);

/*
Add an effect to a tfx_stage_t from an effect template
* @param pm                         A pointer to an initialised tfx_stage_t.
* @param effect_template			The tfx_effect_template_t object that you want to add to the effect manager. It must have already been prepared by calling tfx_CreateEffectTemplate
* @param effect_id					pointer to a tfxEffectID of the effect which will be set after it's been added to the effect manager. This index can then be used to manipulate the effect in the effect manager as it's update
									For example by calling tfx_SetEffectPosition. This will be set to tfxINVALID if the function is unable to add the effect to the effect manager if it's out of space and reached it's effect limit.
  @returns							True if the effect was succesfully added.
*/
tfxAPI tfxEffectID tfx_AddEffectTemplateToStage(tfx_stage pm, tfx_effect_template effect);

/*
Add an effect to a tfx_stage_t. Generally you should always call tfx_AddEffectTemplateToStage and use templates to organise your effects but if you want to just
test things out you can add an effect direct from a library using this command.
* @param pm							A pointer to an initialised tfx_stage_t. 
* @param effect						tfx_effect_descriptor_t object that you want to add to the effect manager.
* @param effect_id					pointer to a tfxEffectID of the effect which will be set after it's been added to the effect manager. This index can then be used to manipulate the effect in the effect manager as it's update
									For example by calling tfx_SetEffectPosition. This will be set to tfxINVALID if the function is unable to add the effect to the effect manager if it's out of space and reached it's effect limit.
  @returns							True if the effect was succesfully added.
*/
tfxAPI tfxEffectID tfx_AddRawEffectToStage(tfx_stage pm, tfx_effect_descriptor effect);

tfxAPI bool tfx_EffectIDIsValid(tfxEffectID id);
tfxAPI bool tfx_AnimationIDIsValid(tfxAnimationID id);

/*
Set the warmup time for an effect template. When the effect template is added to the stage it will be warmed up
first so that effectively it begins it's simulation from a set point in time. 
* @param tfx_effect_template		A handle to the effect template
* @param float						The number of milliseconds to warmup for
*/
tfxAPI void tfx_SetEffectTemplateWarmupTime(tfx_effect_template effect, float millisecs);

/*
Set the delta time used when warming up effects. Either match the fixed time step you're using like 1000/60 for 60fps 
or used mulitples of that value for higher performance at the cost of precision
* @param pm							A pointer to an initialised tfx_stage_t. 
* @param double						The delta time measured in milliseconds
*/
tfxAPI void tfx_SetWarmUpDeltaTime(tfx_stage pm, double delta_time);

/*
If an effect is already in the stage and you want to advance it be a given amount of time
without rendering those frame then you can use this function.
* @param pm							A pointer to an initialised tfx_stage_t. 
* @param tfxEffectID				The effect id that is in the stage
* @param float						The time to advance by measured in milliseconds
*/
tfxAPI void tfx_AdvanceEffectTime(tfx_stage pm, tfxEffectID effect_id, float time);

/*
Update a stage which manages effects. If you are interpolating particles in the vertex shader then it's important to only call this function once per frame only and idealy in a fixed step loop.
That means that if your fixed loop has to run twice to catch up (because of low frame rates) then you should still only call this function once but you can multiply the elapsed time
by the number of ticks. The ellapsed time should be the amount of time that has passed since the last frame so in a fixed step loop this will simply be the update rate in millisecs.
For example if you're updating 60 frames per second then elapsed time would be 16.666667. Psuedo code would look something like this:

	TimerAccumulate(game->timer);
	int pending_ticks = TimerPendingTicks(game->timer);	//The number of times the update loop will run this frame.

	while (tfx_TimerDoUpdate(game->timer)) {
		if (pending_ticks > 0) {
			tfx_UpdateStage(&game->pm, FrameLength * pending_ticks);
			//Set the pending ticks to 0 so we don't run the update again this frame.
			pending_ticks = 0;
		}

		TimerUnAccumulate(game->timer);
	}
	TimerSet(game->timer);	//Set the timer and calculate the interpolation value. You can pass that to a uniform or push constant for the shader

	//Only upload the sprite/billboard buffer to the gpu if the effect manager was updated.
	if (TimerUpdateWasRun(game->timer)) {
		RenderParticles(game->pm, game);
	}

* @param pm                    A pointer to an initialised tfx_stage_t.
* @param double                the amount of time that elapsed since the last frame
*/
tfxAPI void tfx_UpdateStage(tfx_stage pm, double elapsed);

/*
Get the image pointer for a sprite. Use this when rendering particles in your renderer. The pointer that is returned will be the pointer that you set in your shape loader function
used when loading an effect library. Generally you shouldn't need to use this function, simply copy the whole instance buffer in the effect manager to your staging buffer to be
copied to the gpu in one go.
* @param pm                    A pointer to an initialised tfx_stage_t. 
* @param property_indexes    The value in the instance_data->property_indexs[i] when iterating over the instance_data in your render function
  @returns                    void* pointer to the image
*/
tfxAPI void *tfx_GetSpriteImagePointer(tfx_stage pm, tfxU32 property_indexes);

/*
Get the total number of instances ready for rendering in the effect manager.
* @param pm                    A pointer to an initialised tfx_stage_t.
*/
tfxAPI tfxU32 tfx_TotalSpriteCount(tfx_stage pm);

/*
Clear all particles, instance_data and effects in a effect manager. If you don't need to use the effect manager again then call tfx_FreeStage to also
free all the memory associated with the effect manager.
* @param pm                        A pointer to an initialised tfx_stage_t.
* @param free_particle_banks    Set to true if you want to free the memory associated with the particle banks and release back to the memory pool
*/
tfxAPI void tfx_ClearStage(tfx_stage pm, bool free_particle_banks, bool free_sprite_buffers);

/*
Free all the memory used in the effect manager. You don't have to call this - tfx_EndTimelineFX frees any stages
that are still alive - but you can free them earlier if you want to reclaim the memory. Calling it is safe either
way, the stage is deregistered from the shutdown sweep so it won't be freed twice.
* @param pm                        A pointer to an initialised tfx_stage_t.
*/
tfxAPI void tfx_FreeStage(tfx_stage pm);

//[Effects functions for altering effects that are currently playing out in a effect manager]

/*
Expire an effect by telling it to stop spawning particles. This means that the effect will eventually be removed from the effect manager after all of it's remaining particles have expired.
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index    The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToStage
*/
tfxAPI void tfx_SoftExpireEffect(tfx_stage pm, tfxEffectID effect_index);

/*
Soft expire all the effects in a effect manager so that the particles complete their animation first
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
*/
tfxAPI void tfx_SoftExpireAll(tfx_stage pm);

/*
Expire an effect by telling it to stop spawning particles and remove all associated particles immediately.
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index    The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToStage
*/
tfxAPI void tfx_HardExpireEffect(tfx_stage pm, tfxEffectID effect_index);

/*
Get effect user data
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index    The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToStage
* @returns                void* pointing to the user data set in the effect. See tfx_effect_template_t::SetUserData() and tfx__set_effect_user_data()
*/
tfxAPI void *tfx_GetEffectUserData(tfx_stage pm, tfxEffectID effect_index);

/*
More for use in the editor, this function updates emitter base values for any effects that are currently running after their graph values have been changed.
*/
tfxAPI void tfx_UpdateStageBaseValues(tfx_stage pm);

/*
Set the effect manager camera. This is used to calculate particle depth if you're using depth ordered particles so it needs to be updated each frame.
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
* @param front            An array of 3 floats representing a normalised 3d vector describing the direction that the camera is pointing
* @param position        An array of 3 floats representing the position of the camera in 3d space
*/
tfxAPI void tfx_SetStageCamera(tfx_stage pm, float front[3], float position[3]);

/*
Each effect in the effect manager can have bounding box which you can decide to keep updated or not if you wanted to do any offscreen culling of effects. Theres some
extra overhead to keep the bounding boxes updated but that can be made back if you have a number of effect particles offscreen that don't need to be drawn.
* @param pm					A pointer to a tfx_stage_t where the effect is being managed
* @param yesno				Set to true or false if you want the bounding boxes to be udpated.
*/
tfxAPI void tfx_KeepBoundingBoxesUpdated(tfx_stage pm, bool yesno);

/*
Set the effect user data for an effect already added to a effect manager
* @param pm					A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index		The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param user_data			A void* pointing to the user_data that you want to store in the effect
*/
tfxAPI void tfx_SetEffectUserData(tfx_stage pm, tfxU32 effect_index, void *data);

/*
Force a effect manager to only run in single threaded mode. In other words, only use the main thread to update particles
* @param pm                A pointer to a tfx_stage_t.
* @param switch_on        true or false to use a single thread or not
*/
tfxAPI void tfx_ForceStageSingleThreaded(tfx_stage pm, bool switch_on);

/*
Get the transform vectors for a sprite's previous position so that you can use that to interpolate between that and the current sprite position
* @param pm                A pointer to a tfx_stage_t.
* @param layer            The index of the sprite layer
* @param index            The sprite index of the sprite that you want the captured sprite for.
* @param position         This should be a pointer to a vec3 that you pass in that will get loaded with the position of the instance
*/
tfxAPI void tfx_GetCapturedInstanceTransform(tfx_stage pm, tfxU32 layer, tfxU32 index, float out_position[3]);

/*
Get the end index offset into the sprite memory for sprite data containing a pre recorded effect animation. 
a for loop to iterate over the instance_data in a pre-recorded effect
* @param sprite_data    A pointer to tfx_sprite_data_t containing all the instance_data and frame data
* @param frame            The index of the frame you want the end index for
* @param layer            The sprite layer
* @returns                tfxU32 containing the end offset
*/
tfxAPI tfxU32 tfx_SpriteDataEndIndex(tfx_sprite_data_t *sprite_data, tfxU32 frame, tfxU32 layer);

/*
Get the end index offset into the ribbon memory for ribbon data containing a pre recorded effect animation. 
a for loop to iterate over the instance_data in a pre-recorded effect
* @param sprite_data    A pointer to tfx_sprite_data_t containing all the instance_data and frame data
* @param frame            The index of the frame you want the end index for
* @returns                tfxU32 containing the end offset
*/
tfxAPI tfxU32 tfx_RibbonDataEndIndex(tfx_sprite_data_t *sprite_data, tfxU32 frame);

/*
Make a effect manager stop spawning. This will mean that all emitters in the effect manager will no longer spawn any particles so all currently running effects will expire
as the remaining particles come to the end of their life. Any single particles will also get flagged to expire
* @param pm                A pointer to a tfx_stage_t.
* @param yesno            True = disable spawning, false = enable spawning
*/
tfxAPI void tfx_DisableStageSpawning(tfx_stage pm, bool yesno);

/*
Get the buffer of effect indexes in the effect manager.
* @param pm               A pointer to a tfx_stage_t.
* @param depth            The depth of the list that you want. 0 are top level effects and anything higher are sub effects within those effects
* @param count			  A pointer to an int that you can pass in that will be filled with the count of effects in the array
* @returns                Pointer to the array of effect indexes
*/
tfxAPI tfx_effect_index_t *tfx_GetStageEffectBuffer(tfx_stage pm, int *count);

/*
Get the buffer of emitter indexes in the effect manager.
* @param pm                A pointer to a tfx_stage_t.
* @param depth            The depth of the list that you want. 0 are top level emitters and anything higher are sub emitters within those effects
* @param count			  A pointer to an int that you can pass in that will be filled with the count of emitters in the array
* @returns                Pointer to the tfxvec of effect indexes
*/
tfxAPI tfxU32 *tfx_GetStageEmitterBuffer(tfx_stage pm, int *count);

/*
Error-handling contract for the effect manager functions that take a tfxEffectID (below):
Every such function validates the id. In debug builds an invalid id trips an assert at the call site.
In release builds the call is a documented no-op instead of undefined behaviour: mutators return
without touching any state, functions that return a pointer return NULL, and functions that return a
value return 0. This means a stale effect id (for example one that expired a frame earlier than the
caller expected) can never write through a bad index and corrupt the effect manager. Effect ids are
produced by tfx_AddEffectTemplateToStage and should be treated as invalid once the effect has
expired or been removed. The same contract applies to the animation manager functions that take a
tfxAnimationID.
*/

/*
Set the position of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param x                The x value of the position
* @param y                The y value of the position
* @param z                The z value of the position
*/
tfxAPI void tfx_SetEffectPosition(tfx_stage pm, tfxEffectID effect_index, float x, float y, float z);

/*
Set the position of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param position        A float[3] array containing the x, y and z coordinates
*/
tfxAPI void tfx_SetEffectPositionVec3(tfx_stage pm, tfxEffectID effect_index, float position[3]);

/*
Move an Effect by a specified amount relative to the effect's current position
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param amount            A float[3] array containing the amount to move in the x, y and z planes
*/
tfxAPI void tfx_MoveEffectVec3(tfx_stage pm, tfxEffectID effect_index, float amount[3]);

/*
Move an Effect by a specified amount relative to the effect's current position
* @param pm               A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index     The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param x                The amount to move in the x plane
* @param y                The amount to move in the y plane
* @param z                The amount to move in the z plane
*/
tfxAPI void tfx_MoveEffect(tfx_stage pm, tfxEffectID effect_index, float x, float y, float z);

/*
Get the current position of an effect
* @param pm              A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @return                float[3] array containing the effect position
*/
tfxAPI void tfx_GetEffectPositionVec3(tfx_stage pm, tfxEffectID effect_index, float out_position[3]);

/*
You can use this function to get the billboard buffer of a specific effect. 
* @param pm						A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index			The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param tfxU32					Pass in a pointer to a tfxU32 which will be set to the number of instance_data in the buffer.
* @return						tfx_instance_t pointer to the buffer
*/
tfxAPI tfx_instance_t *tfx_GetEffectInstanceBuffer(tfx_stage pm, tfxEffectID effect_index, tfxU32 *sprite_count);

/*
You can use this function to get each billboard buffer for every effect that is currently active in the effect manager. Generally you would call this inside a for loop for each layer.
* @param pm						A pointer to a tfx_stage_t where the effect is being managed
* @param tfx_sprite_billboard_t	Pass in a pointer which will be set to the current sprite buffer containing all of the sprite data for this frame.
* @param tfx_effect_instance_data_t   Pass in a second pointer which will be set to the tfx_effect_instance_data_t containing all of the sprite buffer data. This can be used to gain access to all the sprite data if using double buffered instance_data (to interpolated with the previous frame).
* @param tfxU32					Pass in a pointer to a tfxU32 which will be set to the number of instance_data in the buffer.
* @return						true or false if the next billboard buffer was found. False will be returned once there are no more effect sprite buffers in the effect manager
*/
tfxAPI bool tfx_GetNextInstanceBuffer(tfx_stage pm, tfx_instance_t **sprites_soa, tfx_effect_instance_data_t **effect_sprites, tfxU32 *sprite_count);

/*After calling GetNextBillboard/SpriteBuffer in a while loop you can call this to reset the index for the next frame
* @param pm						A pointer to a tfx_stage_t
*/
tfxAPI void tfx_ResetInstanceBufferLoopIndex(tfx_stage pm);

/*
Set the roll of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current roll of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param roll            A float of the amount that you want to set the roll too
*/
tfxAPI void tfx_SetEffectRoll(tfx_stage pm, tfxEffectID effect_index, float roll);

/*
Set the pitch of a effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current pitch of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param pitch            A float of the amount that you want to set the pitch too
*/
tfxAPI void tfx_SetEffectPitch(tfx_stage pm, tfxEffectID effect_index, float pitch);

/*
Set the yaw of a effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current yaw of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param yaw            A float of the amount that you want to set the yaw too
*/
tfxAPI void tfx_SetEffectYaw(tfx_stage pm, tfxEffectID effect_index, float yaw);

/*
Set the width of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current width of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param width            A float of the amount that you want to set the width multiplier too. The width multiplier will multiply all widths of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void tfx_SetEffectWidthMultiplier(tfx_stage pm, tfxEffectID effect_index, float width);

/*
Set the height of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current height of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param height            A float of the amount that you want to set the height multiplier too. The height multiplier will multiply all heights of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void tfx_SetEffectHeightMultiplier(tfx_stage pm, tfxEffectID effect_index, float height);

/*
Set the depth of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current depth of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param depth            A float of the amount that you want to set the depth multiplier too. The depth multiplier will multiply all heights of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void tfx_SetEffectDepthMultiplier(tfx_stage pm, tfxEffectID effect_index, float depth);

/*
Set the life multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current life of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param life            A float of the amount that you want to set the life multiplier too. The life mulitplier will affect how long all particles emitted within the effect will last before expiring.
*/
tfxAPI void tfx_SetEffectLifeMultiplier(tfx_stage pm, tfxEffectID effect_index, float life);

/*
Set the particle width multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current particle width of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param width            A float of the amount that you want to set the particle width multiplier too. The particle width mulitplier will affect the width of each particle if the emitter has a non uniform particle size, otherwise
						it will uniformly size the particle
*/
tfxAPI void tfx_SetEffectParticleWidthMultiplier(tfx_stage pm, tfxEffectID effect_index, float width);

/*
Set the particle height multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current particle width of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param height            A float of the amount that you want to set the particle height multiplier too. The particle height mulitplier will affect the height of each particle if the emitter has a non uniform particle size, otherwise
						this function will have no effect.
*/
tfxAPI void tfx_SetEffectParticleHeightMultiplier(tfx_stage pm, tfxEffectID effect_index, float height);

/*
Set the velocity multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current velocity of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param velocity        A float of the amount that you want to set the particle velocity multiplier too. The particle velocity mulitplier will affect the base velocity of a particle at spawn time.
*/
tfxAPI void tfx_SetEffectVelocityMultiplier(tfx_stage pm, tfxEffectID effect_index, float velocity);

/*
Set the spin multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current spin of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param spin            A float of the amount that you want to set the particle spin multiplier too. The particle spin mulitplier will affect the base spin of a particle at spawn time.
*/
tfxAPI void tfx_SetEffectSpinMultiplier(tfx_stage pm, tfxEffectID effect_index, float spin);

/*
Set the intensity multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current intensity of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param intensity        A float of the amount that you want to set the particle intensity multiplier too. The particle intensity mulitplier will instantly affect the opacity of all particles currently emitted by the effect.
*/
tfxAPI void tfx_SetEffectIntensityMultiplier(tfx_stage pm, tfxEffectID effect_index, float intensity);

/*
Set the splatter multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current splatter of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param splatter        A float of the amount that you want to set the particle splatter multiplier too. The particle splatter mulitplier will change the amount of random offset all particles emitted in the effect will have.
*/
tfxAPI void tfx_SetEffectSplatterMultiplier(tfx_stage pm, tfxEffectID effect_index, float splatter);

/*
Set the weight multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current weight of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param weight            A float of the amount that you want to set the particle weight multiplier too. The particle weight mulitplier will change the weight applied to particles in the effect at spawn time.
*/
tfxAPI void tfx_SetEffectWeightMultiplier(tfx_stage pm, tfxEffectID effect_index, float weight);

/*
Set the overal scale of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current weight of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param overall_scale    A float of the amount that you want to set the overal scale to. The overal scale is an simply way to change the size of an effect
*/
tfxAPI void tfx_SetEffectOverallScale(tfx_stage pm, tfxEffectID effect_index, float overall_scale);

/*
Set the base noise offset for an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param noise_offset    A float of the amount that you want to set the effect noise offset to. By default when an effect is added to a effect manager a random noise offset will be set based on the Base Noise Offset Range property. Here you can override that
						value by setting it here. The most ideal time to set this would be immediately after you have added the effect to the effect manager, but you could call it any time you wanted for a constantly changing noise offset.
*/
tfxAPI void tfx_SetEffectBaseNoiseOffset(tfx_stage pm, tfxEffectID effect_index, float noise_offset);

/*
Get the name of an effect
* @param pm                A pointer to the effect
* @returns                const char * name
*/
tfxAPI const char *tfx_GetEffectName(tfx_effect_descriptor effect);


//--------------------------------
//Animation_manager
//--------------------------------

/*
Set the position of an animation instance
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index            The index of the effect. This is the index returned when calling tfx_AddAnimationInstance
* @param position                A float[3] array containing the x, y and z coordinates
*/
tfxAPI void tfx_SetAnimationPosition(tfx_animation_manager animation_manager, tfxAnimationID animation_id, float position[3]);

/*
Set the scale of an animation instance
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index            The index of the effect. This is the index returned when calling tfx_AddAnimationInstance
* @param scale                    A multiplier that will determine the overal size/scale of the effect
*/
tfxAPI void tfx_SetAnimationScale(tfx_animation_manager animation_manager, tfxAnimationID animation_id, float scale);

/*
Get an animation instance from an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param tfxAnimationID            The index of the effect. This is the index returned when calling tfx_AddAnimationInstance
* @returns pointer to instance    Pointer to a tfx_animation_instance_t
*/
tfxAPI tfx_animation_instance_t *tfx_GetAnimationInstance(tfx_animation_manager animation_manager, tfxAnimationID animation_id);

/*
Initialise an Animation Manager for use with instance data. This must be run before using an animation manager. An animation manager is used
to playback pre recorded particle effects as opposed to using a effect manager that simulates the particles in
real time. This pre-recorded data can be uploaded to the gpu for a compute shader to do all the interpolation work
to calculate the state of particles between frames for smooth animation.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param max_instances            The maximum number of animation instances that you want to be able to play at one time.
* @param initial_capacity        Optionally, you can set an initial capacity for the sprite data. The data will grow if you add
								beyond this amount but it gives you a chance to reserve a decent amount to start with to
								save too much mem copies as the data grows
*/
tfxAPI tfx_animation_manager tfx_CreateAnimationManager(tfxU32 max_instances, tfxU32 initial_sprite_data_capacity);

/*
Set the callback that you can use to determine whether or not a tfx_animation_instance_t should be added to the next frame's render queue. You can use this
to cull instances that are outside of the view frustum for example
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param callback                Pointer to the callback you want to use. It must have the following signature:
								bool(*maybe_render_instance_callback(tfx_animation_manager animation_manager, tfx_animation_instance_t *instance, tfx_frame_meta_t *meta, void *user_data))
								Values passed into the callback function are a pointer to the animation manager, a pointer to the instance being processed, a pointer to
								the frame meta of the instance, this will contain the bounding box and radius of the instance from the current frame of the instance and a pointer
								to any user data that you set that might contain the camera frustum that you want to check against.
*/
tfxAPI void tfx_SetAnimationManagerInstanceCallback(tfx_animation_manager animation_manager, tfx_maybe_render_instance_callback maybe_render_instance_callback);

/*
Get the sprite data settings for an effect in a library. Sprite data settings are the settings for an effect in the editor relating to setting up pre-baked effects
* @param library				Pointer to the tfx_library_t where the effect is stored
* @param effect					Pointer the the effect that you want the sprite settings for.
* @returns						Pointer to the tfx_sprite_data_settings
*/
tfxAPI tfx_sprite_data_settings_t *tfx_GetEffectSpriteDataSettings(tfx_library library, tfx_effect_descriptor effect);

/*
Get the sprite data settings for an effect in a library by it's path. Sprite data settings are the settings for an effect in the editor relating to setting up pre-baked effects
* @param library				Pointer to the tfx_library_t where the effect is stored
* @param path					const char* string of the path to the effect. Must be the path to a root effect.
* @returns						Pointer to the tfx_sprite_data_settings
*/
tfxAPI tfx_sprite_data_settings_t *tfx_GetEffectSpriteDataSettingsByPath(tfx_library library, const char *path);

/*
Get the index offset into the sprite memory for sprite data containing a pre recorded effect animation. Can be used along side tfx_SpriteDataEndIndex to create
a for loop to iterate over the instance_data in a pre-recorded effect
* @param sprite_data    A pointer to tfx_sprite_data_t containing all the instance_data and frame data
* @param frame            The index of the frame you want the offset for
* @param layer            The sprite layer
* @returns                tfxU32 containing the index offset
*/
tfxAPI tfxU32 tfx_SpriteDataIndexOffset(tfx_sprite_data_t *sprite_data, tfxU32 frame, tfxU32 layer);

/*
Get the index offset into the ribbon memory for ribbon data containing a pre recorded effect animation. Can be used along side tfx_RibbonDataEndIndex to create
a for loop to iterate over the instance_data in a pre-recorded effect
* @param sprite_data      A pointer to tfx_sprite_data_t containing all the instance_data and frame data
* @param frame            The index of the frame you want the offset for
* @returns                tfxU32 containing the index offset
*/
tfxAPI tfxU32 tfx_RibbonDataIndexOffset(tfx_sprite_data_t *sprite_data, tfxU32 frame);

/*
Set the user data in a tfx_animation_manager_t which can get passed through to callback functions when updated the animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param user_data                void* pointer to the data that you want to set
*/
tfxAPI void tfx_SetAnimationManagerUserData(tfx_animation_manager animation_manager, void *user_data);

/*
Add sprite data to an animation manager sprite data buffer from an effect. This will record the
animation if necessary and then convert the sprite data to tfx_sprite_data_t ready for uploading
to the GPU
* @param animation_manager       A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index            The index of the effect. This is the index returned when calling tfx_AddAnimationInstance
* @param position                A float[3] array containing the x, y and z coordinates
* @param sprite_data             Optional advanced usage if you want to use sprite data you recorded in another process. Pass 0 to just use the record or lookup pre-recorded data for the effect.
*/
tfxAPI void tfx_AddSpriteData(tfx_animation_manager animation_manager, tfx_effect_descriptor effect, tfx_stage pm, float camera_position[3], tfx_sprite_data_t *sprite_data);

/*
Add an animation instance to the animation manager.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param path                    tfxKey path hash of the effect name and path: effect.path_hash
* @param start_frame            Starting frame of the animation
* @returns                        The index id of the animation instance. You can use this to reference the animation when changing position, scale etc
								Return tfxINVALID if there is no room in the animation manager
*/
tfxAPI tfxAnimationID tfx_AddAnimationInstanceByKey(tfx_animation_manager animation_manager, tfxKey path, tfxU32 start_frame);

/*
Add an animation instance to the animation manager.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param path                    const char * name of the effect. If the effect was in a folder then specify the whole path
* @param start_frame            Starting frame of the animation
* @returns                        The index id of the animation instance. You can use this to reference the animation when changing position, scale etc
								Return tfxINVALID if there is no room in the animation manager
*/
tfxAPI tfxAnimationID tfx_AddAnimationInstance(tfx_animation_manager animation_manager, const char *path, tfxU32 start_frame);

/*
Update an animation manager to advance the time and frames of all instances currently playing.
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to update
* @param start_frame            Starting frame of the animation
*/
tfxAPI void tfx_UpdateAnimationManager(tfx_animation_manager animation_manager, float elapsed);

/*
Add an effect's shapes to an animation manager. You can use this function if you're manually recording particle effects and adding them to an animation
manager rather then just using the editor.
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to update
* @param effect                    A pointer to the effect whose shapes you want to add
*/
tfxAPI void tfx_AddEffectShapes(tfx_animation_manager animation_manager, tfx_effect_descriptor effect);

/*
Update an animation manager so that the effects do not expire they just loop forever instead regardless of whether they're a looped effect or not.
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to update
*/
tfxAPI void tfx_CycleAnimationManager(tfx_animation_manager animation_manager);

/*
Clears all animation instances currently in play in an animation manager, resulting in all currently running animations
from being drawn
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to clear
*/
tfxAPI void tfx_ClearAllAnimationInstances(tfx_animation_manager animation_manager);

/*
Clears all data from the animation manager including sprite data, metrics and instances. Essentially resetting everything back to
it's initialisation point
from being drawn
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to reset
*/
tfxAPI void tfx_ResetAnimationManager(tfx_animation_manager animation_manager);

/*
Frees all data from the animation manager including sprite data, metrics and instances and also the handle itself.
You don't have to call this - tfx_EndTimelineFX frees any animation managers that are still alive - but you can
free them earlier if you want to reclaim the memory. Calling it is safe either way, the animation manager is
deregistered from the shutdown sweep so it won't be freed twice.
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to free
*/
tfxAPI void tfx_FreeAnimationManager(tfx_animation_manager animation_manager);

/*
Get the tfx_animation_buffer_metrics_t from an animation manager. This will contain the info you need to upload the sprite data,
offsets and animation instances to the GPU. Only offsets and animation instances need to be uploaded to the GPU each frame. Sprite
data can be done ahead of time.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @returns                        tfx_animation_buffer_metrics_t containing buffer sizes
*/
tfxAPI tfx_animation_buffer_metrics_t tfx_GetAnimationBufferMetrics(tfx_animation_manager animation_manager);

/*
Get the total number of instance_data that need to be drawn by an animation manager this frame. You can use this in your renderer
to draw your sprite instances
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @returns                        tfxU32 of the number of instance_data
*/
tfxAPI tfxU32 tfx_GetTotalSpritesThatNeedDrawing(tfx_animation_manager animation_manager);

/*
Get the total number of ribbons that need to be drawn by an animation manager this frame. You can use this in your renderer
to draw your sprite instances
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @returns                        tfxU32 of the number of ribbons to draw
*/
tfxAPI tfxU32 tfx_GetTotalRibbonsThatNeedDrawing(tfx_animation_manager animation_manager);

/*
Get the total number of instances being processed by an animation manager. This will not necessarily be the same number as
the instances being rendered if some are being culled in your custom callback if your using one.
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to clear
* @returns int                    The number of instances being updated
*/
tfxAPI tfxU32 tfx_GetTotalInstancesBeingUpdated(tfx_animation_manager animation_manager);

/*
Create the image data required for shaders from a TimelineFX library. The image data will contain data such as uv coordinates. Once you have built the data you can use GetLibraryImageData to get the buffer
and upload it to the gpu.
* @param animation_manager		  A pointer to an tfx_animation_manager_t object
* @param shapes                   A pointer to a tfx_gpu_shapes_t object which will fill a buffer with all the shapes
* @param uv_lookup                A function pointer to a function that you need to set up in order to get the uv coordinates from whatever renderer you're using
*/
tfxAPI void tfx_BuildAnimationManagerGPUShapeData(tfx_animation_manager animation_manager, tfx_gpu_shapes shapes, tfx_uv_lookup uv_lookup);

/*
Get a pointer to the particle shapes data in the animation manager. This can be used with tfx_BuildGPUShapeData when you want to upload the data to the GPU
* @param animation_manager        A pointer the tfx_animation_manager_t
*/
tfxAPI tfx_image_data_t *tfx_GetParticleShapesAnimationManager(tfx_animation_manager animation_manager, int *count);

/*
Get the total number of instance_data in an animation manger's sprite data buffer
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns tfxU32                The number of instance_data in the buffer
*/
tfxAPI tfxU32 tfx_GetTotalSpriteDataCount(tfx_animation_manager animation_manager);

/*
Get the total byte size of instance_data in an animation manger's sprite data buffer
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t                 The number of instance_data in the buffer
*/
tfxAPI size_t tfx_GetSpriteDataSizeInBytes(tfx_animation_manager animation_manager);

/*
Get the total byte size of ribbon_data in an animation manger's ribbon data buffer
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns tfxU32                The number of instance_data in the buffer
*/
tfxAPI size_t tfx_GetRibbonDataSizeInBytes(tfx_animation_manager animation_manager);

/*
Get the buffer memory address for the sprite data in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                A pointer to the sprite data memory
*/
tfxAPI void *tfx_GetSpriteDataBufferPointer(tfx_animation_manager animation_manager);

/*
Get the buffer memory address for the ribbon data in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                  A pointer to the sprite data memory
*/
tfxAPI void *tfx_GetRibbonDataBufferPointer(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the offsets buffer in an animation manager for sprite data
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t                Size in bytes of the offsets buffer
*/
tfxAPI size_t tfx_GetOffsetsSizeInBytes(tfx_animation_manager animation_manager);

/*
Get the buffer pointer for the offsets buffer in an animation manager for sprite data
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t                Size in bytes of the offsets buffer
*/
tfxAPI void *tfx_GetOffsetsBufferPointer(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the offsets buffer in an animation manager for ribbons
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @returns size_t                Size in bytes of the offsets buffer
*/
tfxAPI size_t tfx_GetRibbonOffsetsSizeInBytes(tfx_animation_manager animation_manager);

/*
Get the buffer pointer for the offsets buffer in an animation manager for ribbon data
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @returns size_t                Size in bytes of the offsets buffer
*/
tfxAPI void *tfx_GetRibbonOffsetsBufferPointer(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the render queue of animation instances buffer in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t                Size in bytes of the instances buffer
*/
tfxAPI size_t tfx_GetAnimationInstancesSizeInBytes(tfx_animation_manager animation_manager);

/*
Get the count of of animation instances in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns tfxU32                 Count of the instances buffer
*/
tfxAPI tfxU32 tfx_GetAnimationInstancesCount(tfx_animation_manager animation_manager);

/*
Get the buffer pointer of of animation instances in an animation manager which you'll need to copy
the data to the GPU
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                  Pointer to the instances buffer
*/
tfxAPI void *tfx_GetAnimationInstancesBufferPointer(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the animation emitter properties list
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t                Size in bytes of the properties bufffer
*/
tfxAPI size_t tfx_GetAnimationEmitterPropertySizeInBytes(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the animation ribbon properties list
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @returns size_t                Size in bytes of the properties bufffer
*/
tfxAPI size_t tfx_GetAnimationRibbonPropertySizeInBytes(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the animation ribbon segments list
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon segment data from
* @returns size_t                Size in bytes of the ribbon segments bufffer
*/
tfxAPI size_t tfx_GetAnimationRibbonSegmentsSizeInBytes(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the animation ribbon data list
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @returns size_t                Size in bytes of the ribbon segments bufffer
*/
tfxAPI size_t tfx_GetAnimationRibbonDataSizeInBytes(tfx_animation_manager animation_manager);

/*
Get the number of emitter properties being using by the animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns tfxU32                Number of emitter properties
*/
tfxAPI tfxU32 tfx_GetAnimationEmitterPropertyCount(tfx_animation_manager animation_manager);

/*
Get the buffer memory address for the sprite emitter properties in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                A pointer to the sprite emitter properties data memory
*/
tfxAPI void *tfx_GetAnimationEmitterPropertiesBufferPointer(tfx_animation_manager animation_manager);

/*
Get the buffer memory address for the ribbon property data in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                A pointer to the ribbon properties data memory
*/
tfxAPI void *tfx_GetAnimationRibbonPropertiesBufferPointer(tfx_animation_manager animation_manager);

/*
Get the buffer memory address for the ribbon segments data in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                  A pointer to the ribbon segments data memory
*/
tfxAPI void *tfx_GetAnimationRibbonSegmentsBufferPointer(tfx_animation_manager animation_manager);

/*
Get the buffer memory address for the ribbon segments data in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                  A pointer to the ribbon segments data memory
*/
tfxAPI void *tfx_GetAnimationRibbonDataBufferPointer(tfx_animation_manager animation_manager);

/*
Returns true or false if the animation manager contains effects with ribbons
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @returns bool 	              True if the animation manager has ribbons
*/
tfxAPI bool tfx_AnimationManagerHasRibbons(tfx_animation_manager animation_manager);

/*
Returns true or false if the animation manager contains emitters that use animated shapes
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @returns bool 	              True if the animation manager does contain animated shapes
*/
tfxAPI bool tfx_HasAnimatedShapes(tfx_animation_manager animation_manager);

/*
Set the callback used by the animation manager to determin if an animation should be rendered or not
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @param callback 	              The function callback
*/
tfxAPI void tfx_SetAnimationManagerCullCallback(tfx_animation_manager animation_manager, tfx_maybe_render_instance_callback callback);

tfxAPI size_t tfx_CalculateAnimationInstanceBufferSize(size_t instance_count);
tfxAPI size_t tfx_CalculateAnimationOffsetsBufferSize(size_t instance_count);

//--------------------------------
//Effect_templates
//--------------------------------
 
/*
Prepare a tfx_effect_template_t that you can use to customise effects in the library in various ways before adding them into a effect manager for updating and rendering. Using a template like this
means that you can tweak an effect without editing the base effect in the library.
* @param library                    A reference to a tfx_library_t that should be loaded with tfx_LoadEffectLibrary
* @param name                       The name of the effect in the library that you want to use for the template. If the effect is in a folder then use normal pathing: "My Folder/My effect"
//Returns handle					Handle to the newly created effect template or nullptr if the effect couldn't be found in the library
*/
tfxAPI tfx_effect_template tfx_CreateEffectTemplate(tfx_library library, const char *name);
 
/*
Delete an effect template and free all memory associated with it
* @param effect_template            A handle to the effect template to be deleted
//Returns handle					Handle to the newly created effect template or nullptr if the effect couldn't be found in the library
*/
tfxAPI void tfx_FreeEffectTemplate(tfx_effect_template effect_template);

/*
Reset an effect template and make it empty so you can use it to store another effect.
* @param t                        A pointer to a tfx_effect_template_t
*/
tfxAPI void tfx_ResetTemplate(tfx_effect_template t);

/*
Get the root effect from the template
* @param t                        A pointer to a tfx_effect_template_t
* @returns                        A pointer to the root effect
*/
tfxAPI tfx_effect_descriptor tfx_GetEffectFromTemplate(tfx_effect_template t);

/*
Get an emitter or sub effect from an effect template.
* @param t                        A pointer to a tfx_effect_template_t
* @param path                     A path to the emitter or sub effect that you want to retrieve. Must be a valid path. Example path might be: "Explosion/Smoke"
* @returns                        A pointer to the root effect
*/
tfxAPI tfx_effect_descriptor tfx_GetEmitterFromTemplate(tfx_effect_template t, const char *path);

/*
Get an emitter path that an emitter is using. The emitter must have the path emission type set or nullptr will be returned
* @param t                        A pointer to a tfx_effect_descriptor_t
* @param path                     A path to the emitter or sub effect that you want to retrieve. Must be a valid path. Example path might be: "Explosion/Smoke"
* @returns                        A pointer to the root effect
*/
tfxAPI tfx_emitter_path_t *tfx_GetEmitterPath(tfx_effect_descriptor e);

/*
Set the user data for any effect or emitter in the effect template. This user data will get passed through to any update callback functions
* @param t                        A pointer to a tfx_effect_template_t
* @param path                     A path to the effect or emitter in the effect template
* @param data                     A pointer to the user data
*/
tfxAPI void tfx_SetTemplateUserData(tfx_effect_template t, const char *path, void *data);

/*
Set the user data for the root effect in an effect template
* @param t                        A pointer to a tfx_effect_template_t
* @param data                     A pointer to the user data
*/
tfxAPI void tfx_SetTemplateEffectUserData(tfx_effect_template t, void *data);

/*
Set the same user data for all effects and emitters/sub effects in the effect template
* @param t                        A pointer to a tfx_effect_template_t
* @param data                     A pointer to the user data that will be set to all effects and emitters in the template
*/
tfxAPI void tfx_SetTemplateUserDataAll(tfx_effect_template t, void *data);

/*
Set an update callback for the root effect in the effect template.
* @param t                        A pointer to a tfx_effect_template_t
* @param update_callback          A pointer to the call back function
*/
tfxAPI void tfx_SetTemplateEffectUpdateCallback(tfx_effect_template t, void(*update_callback)(tfx_stage pm, tfxEffectID effect_index));

/*
Pre-record this effect into a sprite cache so that you can play the effect back without the need to actually caclulate particles in realtime.
	* @param pm					  Reference to a pm that will be used to run the particle simulation and record the sprite data
	* @param path				  const *char of a path to the emitter in the effect.Must be a valid path, for example: "My Effect/My Emitter"
	* @param camera				  Array of 3 floats with the camera position (only needed for effects that are sorted by depth)
*/
tfxAPI void tfx_RecordTemplateEffect(tfx_effect_template t, tfx_stage pm, float update_frequency, float camera_position[3]);

/*
Pre-record this effect into a sprite cache so that you can play the effect back without the need to actually caclulate particles in realtime. This version
of the function allows you to pass in specific settings
	* @param pm					  Reference to a pm that will be used to run the particle simulation and record the sprite data
	* @param path				  const *char of a path to the emitter in the effect.Must be a valid path, for example: "My Effect/My Emitter"
	* @param camera				  Array of 3 floats with the camera position (only needed for effects that are sorted by depth)
*/
tfxAPI void tfx_RecordEffect(tfx_effect_descriptor e, tfx_sprite_data_settings_t *settings, tfx_sprite_data_t *sprite_data, tfx_stage pm, float update_frequency, float camera_position[3]);

/*
Disable an emitter within an effect. Disabling an emitter will stop it being added to the effect manager when calling tfx_AddEffectTemplateToStage
* @param path					  const *char of a path to the emitter in the effect. Must be a valid path, for example: "My Effect/My Emitter"
*/
tfxAPI void tfx_DisableTemplateEmitter(tfx_effect_template t, const char *path);

/*
Enable an emitter within an effect so that it is added to the effect manager when calling tfx_AddEffectTemplateToStage. Emitters are enabled by default.
* @param path					  const *char of a path to the emitter in the effect. Must be a valid path, for example: "My Effect/My Emitter"
*/
tfxAPI void tfx_EnableTemplateEmitter(tfx_effect_template t, const char *path);

/*
Scale all nodes on a global graph graph of the effect
* @param tfx_effect_template	  Handle to a valid effect template
* @param graph_index			  tfx_global_graph_index enum of the global graph that you want to scale. 
* @param amount					  A float of the amount that you want to scale the multiplier by.
*/
tfxAPI void tfx_ScaleTemplateGlobalMultiplier(tfx_effect_template t, tfx_global_graph_index graph_index, float scale_amount);

/*
Set the single spawn amount for an emitter. Only affects emitters that have the single spawn flag set.
* @param emitter_path			 const *char of the emitter path
* @param amount					 A float of the amount that you want to set the single spawn amount to.
*/
tfxAPI void tfx_SetTemplateSingleSpawnAmount(tfx_effect_template t, const char *emitter_path, tfxU32 amount);

/*
Scale all nodes on an emitter graph
* @param tfx_effect_template	  Handle to a valid effect template
* @param emitter_path			  const *char of the emitter path
* @param global_type			  tfx_graph_type of the emitter graph that you want to scale. Must be an emitter graph or an assert will be called
* @param scale amount             A float of the amount that you want to scale the graph by.
*/
tfxAPI void tfx_ScaleTemplateEmitterGraph(tfx_effect_template t, const char *emitter_path, tfx_emitter_graph_index graph_index, float scale_amount);

//--------------------------------
//Editing_graphs
//--------------------------------
tfxAPI void tfx_ClearBaseLifetimeGraph(tfx_effect_descriptor emitter, float v );
tfxAPI void tfx_ClearVariationLifetimeGraph(tfx_effect_descriptor emitter, float v);

//--------------------------------
//General_helpers
//--------------------------------

tfxAPI void tfx_GetSpriteScale(void *instance, float out_scale[2]);

tfxAPI inline float tfx_GetDistance(float fromx, float fromy, float tox, float toy) {
	float w = tox - fromx;
	float h = toy - fromy;
	return sqrtf(w * w + h * h);
}

tfxAPI tfx_effect_descriptor tfx_CreateEffectDescriptor(tfx_effect_descriptor_type type);

/*
Create a new list for containing gpu shapes. You can use this list to upload to the GPU so that shaders can have access to particle image data like UV coordinates
* @returns tfx_gpu_shapes                A handle to the shapes list
*/
tfxAPI tfx_gpu_shapes tfx_CreateGPUShapesList(void);

/*
Delete and free all the memory for a tfx_gpu_shapes list.
* @param shapes                 A tfx_gpu_shapes handle
*/
tfxAPI void tfx_FreeGPUShapesList(tfx_gpu_shapes shapes);

/*
Clear the tfx_gpu_shapes list but keep the memory associated with it.
* @param shapes                 A tfx_gpu_shapes handle
*/
tfxAPI void tfx_ClearGPUShapesList(tfx_gpu_shapes shapes);

/*
Get a pointer to the GPU shapes array which you can use in a memcpy to a staging buffer for uploading to a GPU
* @param particle_shapes        A pointer the tfx_gpu_shapes_t
*/
tfxAPI tfx_gpu_image_data_t *tfx_GetGPUShapesArray(tfx_gpu_shapes shapes);

/*
Get the number of shapes in the GPU Shape Data buffer contained within a library. Make sure you call tfx_BuildGPUShapeData first or they'll be nothing to return
* @param shapes                 A tfx_gpu_shapes handle
* @returns tfxU32               The number of shapes in the buffer
*/
tfxAPI tfxU32 tfx_GetGPUShapesCount(tfx_gpu_shapes shapes);

/*
Get the size in bytes of shapes from a tfx_gpu_shapes handle. You can use this when uploading the shape data to the GPU along with tfx_GetGPUShapesArray
* @param shapes                 A tfx_gpu_shapes handle
* @returns tfxU32               The number of shapes in the buffer
*/
tfxAPI size_t tfx_GetGPUShapesSizeInBytes(tfx_gpu_shapes shapes);

/*
Add a new tfx_gpu_image_data_t object onto a list of gpu shapes
* @param shapes                 A tfx_gpu_shapes handle
* @param image data             The tfx_gpu_image_data_t object that you want to add to the shapes list
* @returns index				An index into the array where the new image data was added
*/
tfxAPI tfxU32 tfx_AddGPUShape(tfx_gpu_shapes shapes, tfx_gpu_image_data_t image_data);

/*
Retrieve image data point from a tfx_gpu_shapes handle containing a list of tfx_gpu_image_data_t objects
* @param shapes                 A tfx_gpu_shapes handle
* @param index					
* @returns index				An index into the array where the new image data was added
*/
tfxAPI tfx_gpu_image_data_t *tfx_GetGPUShape(tfx_gpu_shapes shapes, tfxU32 index);

tfxAPI tfx_version_t tfx_GetVersion(void);

// Helper function to get a good default thread count for thread pools
// Usually hardware threads - 1 to leave a core for the OS/main thread
tfxAPI unsigned int tfx_GetDefaultThreadCount(void);

#endif
