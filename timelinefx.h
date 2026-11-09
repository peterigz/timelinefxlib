#ifndef TFX_LIBRARY_HEADER
#define TFX_LIBRARY_HEADER

#define TFX_VERSION_NAME "Alpha"
#define TFX_VERSION_MAJOR 0
#define TFX_VERSION_MINOR 39
#define TFX_VERSION_PATCH 0

#define TFX_STRINGIFY(x) #x

//Bookmarks are user time marker in an effect that are setup in the editor.
#ifndef tfxMAX_BOOKMARKS
#define tfxMAX_BOOKMARKS 8
#endif
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
//Define tfxTRACY_TRANSIENT instead of tfxTRACY when the library lives in a module that can be unloaded (e.g. a hot reloaded DLL)
//#define tfxTRACY_TRANSIENT
#if defined(tfxTRACY_TRANSIENT) && !defined(tfxTRACY)
#define tfxTRACY
#endif
#define TFX_THREAD_SAFE
//#define TFX_EXTRA_DEBUGGING
#define SSE41		//Steam survey currently has this at 99.83% coverage 12 April 2025. I will probably make this the minimum requirement

//The most layers a library can define. Define it lower to shrink a stage's per layer arrays, but never above 64
#ifndef tfxLAYERS
#define tfxLAYERS 64
#endif
#if tfxLAYERS > 64
#error tfxLAYERS must be 64 or less
#endif

//Enable this to process 8 particles at a time.
//#define tfxUSEAVX

//Enable fused multiply add in simd calculations
//#define tfxUSEFMA

/*
	Timeline FX C++ library

	This library is for implementing particle effects into your games and applications.

	This library is render agnostic, so you will have to provide your own means to render the particles. There are various API functions in the library that help you do this.

	Coding agents integrating the library: read the agents folder next to this header first (gpu-particles.md,
	user-spawn-locations.md, integration-checks.md), and route the library's messages with tfx_SetPrintCallback so you can
	read its "TimelineFX tip:" lines.

	Scale: effects are authored in the editor's world units, +Y up. The editor's default preview camera has a 60 degree
	vertical field of view, sits at (0, 3.5, 7.5) pitched down 30 degrees, and plays the effect 5 units in front of it, but
	artists move it, so treat that as a rough guide only. The right size depends on your game: scale your camera distance,
	or the effect with tfx_SetEffectOverallScale, to fit your world.

	Currently tested on Windows, Linux and MacOS, Intel and Mac based ARM processors.

	Table of contents
	Sections in this header file, you can search for the following keywords to jump to that section:

	[Header_Includes_and_Typedefs]		Base typedefs, version + platform/SIMD detection, and the tfxAPI / tfxINTERNAL markers
	[Enums]								All the definitions for enums and bit flags
	[forward_declarations]				Opaque handles, struct typedefs and the public value-type structs (tfx_random_t, tfx_stage_info_t, ...)
	[Callback_typedefs]					Callback function pointer typedefs
	[API functions]						The main functions for use by users of the library
		-[Random numbers]				Seeded random number generation
		-[Initialisation_functions]		Startup and shutdown timelinefx
		-[Library_functions]			Functions for loading and accessing timelinefx libraries
		-[Particle_Manager_functions]	Create and update functions for stages where the main work is done to update particles every frame
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
#define tfxFORCE_INLINE inline __attribute__((always_inline))
#else
#define tfxINTERNAL static
#define tfxFORCE_INLINE __forceinline
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
typedef tfxU64 tfxSpawnLocationID;
#define tfxINVALID_SPAWN_LOCATION 0xFFFFFFFFFFFFFFFFull
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

//Two signed 16 bit values, which tfx_instance_t scales into its size and intensity fields
typedef struct tfx_int16x2_s {
	union {
		struct {
			int16_t x;
			int16_t y;
		};
		struct { tfxU32 packed; };
	};
} tfx_int16x2_t;

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
	tfx_float32x4_t position;				//The position of the billboard with stretch in w
	tfxU64 quaternion;						//Rotation of the billboard stored as a 16-bit snorm quaternion
	tfx_int16x2_t size;						//Width and height scaled by 32767/64: multiply by 64/32767 to unpack, or call tfx_GetSpriteScale
	tfx_float8x4_t alignment;				//Normalised alignment vector as snorm8 in xyz, and the hue in w as snorm8 where 1 is 180 degrees
	tfx_int16x2_t intensity_gradient_map;	//Intensity in x and gradient map in y, scaled by 32767/128: multiply by 128/32767 to unpack
	tfx_float8x4_t curved_alpha_life;		//Sharpness and dissolve amount value for fading the image plus the age of the particle value packed into 3 bit unorms. The w byte holds the lerp offset in baked sprite data
	tfxU32 indexes;							//GPU property index << 16 | capture flag << 15 | animation frame in the low 13 bits (max 8191), which the shader adds to the emitter's start_frame_index. Bits 13 and 14 are unused. The property index is into tfx_GetGlobalParticlePropertiesBuffer on a stage with global_library_indexes, otherwise into the effect's library
	tfxU32 captured_index;					//Only meaningful when the capture flag in indexes is set: mask with 0x0FFFFFFF for this instance's index in the previous instance buffer. Bit 30 is which buffer, bit 31 is internal. Without the capture flag the instance is a spawn from this update or a particle that died but still holds its slot, and the reference shaders draw it at zero size that frame
} tfx_instance_t;

typedef struct tfx_ribbon_vertex_s {
	tfx_float32x3_t position;
	tfxU32 segment_index;
	tfx_float32x2_t uv_offset_scale;
	tfxU32 ribbon_index;
	tfxU32 normal;							//Octahedral unit normal as two snorm16s, which tfx_ribbon_lit_vertex reads. 0 on a layer that isn't lit
} tfx_ribbon_vertex_t;

typedef struct tfx_gpu_particle_properties_s {
	tfx_float32x2_t image_handle;			
	tfxU32 color_ramp_indexes;			//Color ramp layer << 8 | row. The layer is a global layer in tfx_GetGlobalParticlePropertiesBuffer
	tfxU32 flags;						//Flags like billboard alignment type
	tfxU32 start_frame_index;			//First shape frame, an index into tfx_GetGlobalGPUShapes in tfx_GetGlobalParticlePropertiesBuffer
	float animation_frames;
	float alpha_threshold;				//Pixels below this alpha are discarded on a layer that writes depth
	tfxU32 lighting;					//[wrap unorm8 << 16 | roundness unorm8 << 8], read by the lit entry points on a lit layer
} tfx_gpu_particle_properties_t;

//The light the lit fragments read from the view uniforms at tfxVIEW_LIGHT_OFFSET. Every member is 16 bytes so std140, D3D and WGSL agree
typedef struct tfx_view_light_s {
	tfx_float32x4_t direction;			//Unit length, pointing towards the light. w is unused
	tfx_float32x4_t color;				//The light's color times its intensity. w is unused
	tfx_float32x4_t sky_color;			//Ambient light on normals facing up (+Y)
	tfx_float32x4_t ground_color;		//Ambient light on normals facing down
} tfx_view_light_t;

//Where tfx_view_light_t starts in the view uniforms, after view, projection, screen_size, milliseconds, timer_lerp, update_time and 12 bytes of padding
#define tfxVIEW_LIGHT_OFFSET 160
#define tfxVIEW_UNIFORMS_SIZE (tfxVIEW_LIGHT_OFFSET + 64)

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

//Image formats, can be used in a shape loader to determine how the image is loaded and used
typedef enum tfx_image_format {
	//Any effect library created before file format 6 (editor version alpha 35 or before)
	//will have unknown format and will be auto determined.
	tfx_image_format_unknown = 0,
	tfx_image_format_rgba8_png = 1, 	//4 channel rgba .png
	tfx_image_format_la8_png = 2,		//2 channel luminance and alpha .png
	tfx_image_format_l8_png = 3,		//luminance only .png
	tfx_image_format_a8_png = 4,		//alpha only .png
	tfx_image_format_rgba8_raw = 5,		//4 channel rgba raw bitmap
	tfx_image_format_a_bc4_ktx2 = 6,	//alpha only BC4 .ktx2 with its mip chain. Decodes to (1, 1, 1, a), KTXswizzle 111r
	tfx_image_format_l_bc4_ktx2 = 7,	//luminance only BC4 .ktx2 with its mip chain. Decodes to (l, l, l, 1), KTXswizzle rrr1
	tfx_image_format_la_bc5_ktx2 = 8,	//luminance in red and alpha in green, BC5 .ktx2 with its mip chain. Decodes to (l, l, l, a), KTXswizzle rrrg
	tfx_image_format_rgba_bc7_ktx2 = 9,	//4 channel rgba BC7 .ktx2 with its mip chain
} tfx_image_format;

typedef enum {
	tfxStageSetup_none,
	tfxStageSetup_auto_order_effects,
} tfx_stage_setup;

//Which local axis of an effect is treated as the direction it faces, used by tfx_PointEffectAt. TimelineFX is
//a +Y up library: an emitter with no rotation emits along +Y, so tfxEffectFace_up is the default face and the
//one to use unless the effect was deliberately authored to fire along another axis.
typedef enum {
	tfxEffectFace_up,			//+Y, the axis emitters fire along by default
	tfxEffectFace_down,			//-Y
	tfxEffectFace_right,		//+X
	tfxEffectFace_left,			//-X
	tfxEffectFace_forwards,		//+Z
	tfxEffectFace_backwards,	//-Z
} tfx_effect_face;

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
	tfxEffect_global_hue_index,
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

	//Last so that the GPU graph slots before it don't move
	tfxEmitter_overtime_hue_index,

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

	//Last so that the GPU graph slots before them don't move
	tfxRibbon_overtime_hue_index,
	tfxRibbon_overlength_hue_index,

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
	//Force type and transform type are only used in tfx_graph_list_t. We should properly have a separate
	//enum for graph_list types, but they're here for now.
	tfxForceType,
	tfxTransformType,
	tfxMaxDescriptorTypes
} tfx_effect_descriptor_type;

//Used in file loading - for loading effects library
typedef enum {
	tfxString,
	tfxSInt,
	tfxUInt,
	tfxFloat,
	tfxDouble,
	tfxBool,
	tfxColor,
	tfxUInt64,
	tfxFloat3,
	tfxFloat2,
	tfxAttributeGraph,
	tfxTransformGraph,
	tfxGraphProperty,
} tfx_data_type;

//A property read with tfx_GetDescriptorProperty. type is one of tfxUInt, tfxSInt, tfxFloat, tfxBool, tfxUInt64 or tfxString
typedef struct tfx_property_value_s {
	tfx_data_type type;
	union {
		tfxU32 uint_value;
		int int_value;
		float float_value;
		bool bool_value;
		tfxU64 uint64_value;
		const char *string_value;
	};
} tfx_property_value_t;

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
	tfxErrorCode_some_images_loaded_without_user_ptr            = 1 << 15,
	tfxErrorCode_folder_effect_data_not_found     				= 1 << 16,
	tfxErrorCode_could_not_find_valid_effect_data_in_folder     = 1 << 17,
	tfxErrorCode_shapes_changed_on_disk                         = 1 << 18	//Not a failure: an image was edited outside the editor and its emitters were remapped to it
} tfx_error_flag_bits;

//Errors that mean nothing at all could be read, as opposed to a partially loaded library
#define tfxErrorCode_package_unreadable (tfxErrorCode_incorrect_package_format | tfxErrorCode_unable_to_open_file | tfxErrorCode_unable_to_read_file | tfxErrorCode_wrong_file_size | tfxErrorCode_invalid_format | tfxErrorCode_no_inventory | tfxErrorCode_invalid_inventory)

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
typedef struct tfx_ribbon_bucket_globals_s tfx_ribbon_bucket_globals_t;
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

//Configures a stage on creation. The defaults below are what tfx_CreateStageInfo sets
typedef struct tfx_stage_info_s {
	double warmup_delta_time;				//The frame length tick amount for warming up effects. Higher is more performant at the cost of accuracy. Default 1000/60 ms.
	tfxU32 max_particles;					//The starting number of instances for each layer, which grows when dynamic_sprite_allocation is true. Default 5000.
	tfxU32 max_effects;                     //The maximum number of effects that can be updated at the same time. Must be less than 65535. Default 1000.
	tfxU32 max_emitters;                    //The maximum number of particle emitters across all effects, 0 for 3 times max_effects. Adding an effect whose emitters won't fit fails, see tfx_EffectIDIsValid. Default 0.
	tfxU32 max_ribbon_segments;             //All segments for ribbons are stored in a single buffer. You will need to create buffers for rendering and so whatever you decide the max segments should be your buffers
											//should be big enough to contain all ribbon segments that you might need. You can call tfx_GetSegmentBufferMaxSizeInBytes after creating the stage to get the byte
											//value that you can use to create the buffers. Also note that segments are always created in multiples of 32, so whatever number you put here it will be rounded to the
											//nearest multiple of 32. Must be less than 16777215. Default 131072, about 25MB of ribbon vertex, index and segment buffers if you size them all from it at tessellation 1.
	tfxU32 max_ribbons;						//The maximum number of ribbon instances that can exist at the same time across all emitters in this stage. Ribbons stop spawning once it is reached. Default 4000.
	tfxU32 max_ribbon_emitters;				//The maximum number of ribbon emitters across all effects in this stage, 0 to match max_effects. Must be 65535 or less. Adding an effect whose ribbon emitters won't fit fails, see tfx_EffectIDIsValid. Default 0.
	tfxU32 ribbon_tessellation;				//The tessellation tfx_GetSegmentVertexBufferMaxSizeInBytes and tfx_GetSegmentIndexBufferMaxSizeInBytes size for. A half pipe ribbon emitter tessellated more than this
											//spends whole multiples of its segment count from max_ribbon_segments, so the buffers never overflow. Clamped to at least 1. Default 1.
	tfxU32 multi_threaded_batch_size;       //The size of each batch of particles to be processed when multithreading. Must be a power of 2 and 256 or greater. Default 4096.
	tfxU32 gpu_particle_capacity;			//Pool slots across every run_on_gpu emitter on the stage, fixed for the stage's life. 0, the default, runs those emitters on the CPU. Add up tfx_GetTemplateGPUPoolSlots or tfx_SetTemplateGPUCapacityForLocations for each effect instance that plays at once. Spawn locations share their effect's rings, so size those with tfx_SetTemplateGPUCapacityForLocations rather than by adding more. That exact sum fits at a steady 60 updates a second, slower rates and single call catch-up updates need tfx_SetGPUSizingUpdateRate. It fits when every ring is allocated once, but a stage that keeps adding and expiring effects should leave room for about one more of its largest ring, since freed rings can leave the space in pieces. Each slot is 112 bytes of GPU memory, and tfx_GetGPUParticleBufferSizes gives the exact pool size. See tfx_GetGPUParticleFrame.
	bool double_buffer_sprites;             //Set to true to double buffer instance_data so that you can interpolate between the old and new positions for smoother animations. Default true.
	bool dynamic_sprite_allocation;         //Set to true to automatically resize the sprite buffers if they run out of space. Default true.
	bool auto_order_effects;                //Set to true to sort the effects by depth each frame, so within each layer effects are drawn back to front. Use tfx_SetStageCamera in 3d to set the effect depth to the distance the camera. Default false, true with tfxStageSetup_auto_order_effects.
	bool global_library_indexes;			//Set to true to write indexes into the global buffers that tfx_UpdateGlobalGPUData builds, so one stage can draw effects from any number of libraries with one set of GPU resources. Call tfx_UpdateGlobalGPUData before updating the stage, and only add effects from libraries it placed. See "Drawing effects from more than one library". Default false.
	void *user_data;						//User data that will get passed into the grow_staging_buffer_callback function which you can use to grow the buffer. Default NULL.
	//If you need the staging buffer to be grown dynamically then you can use this call back to do that. It should return true if the buffer was successfully grown or false otherwise. Default NULL.
	bool(*grow_staging_buffer_callback)(tfxU32 new_size, tfx_stage stage, void *user_data);
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
	tfxU32 sort_key;						//The sort key of every effect whose ribbons are in this bucket, see tfx_SetEffectSortKey
	tfxU32 layer;							//The layer of every ribbon emitter in this bucket. Draw it with your depth pipeline when tfx_GetLayerFlags says the layer writes depth
	tfx_ribbon_bucket_globals_t *globals;	//This bucket's push constants, which tfx_GetRibbonDispatchGlobals returns
	tfxU32 next_bucket;						//Where tfx_NextAnimationRibbonDispatch carries on from, unused for a stage
} tfx_ribbon_dispatch_t;

//Every instance on one layer with one sort key, see tfx_GetNextInstanceKeyRange. Both starts index the instance buffer you draw from
typedef struct tfx_instance_key_range_s {
	tfxU32 sort_key;
	tfxU32 instance_start;					//The CPU instances, an index into tfx_GetInstanceBuffer
	tfxU32 instance_count;
	tfxU32 gpu_instance_start;				//The Run on GPU instances, which the kernel writes after the CPU ones
	tfxU32 gpu_instance_count;
} tfx_instance_key_range_t;

//The state of one GPU simulated particle that changes every update, one per pool slot
typedef struct tfx_gpu_particle_hot_s {		//32 bytes
	tfx_float32x3_t position;
	tfxU32 flags;							//tfxGPUParticleFlags
	tfx_float32x3_t velocity;
	float age;								//Milliseconds
} tfx_gpu_particle_hot_t;

//The state of one GPU simulated particle that is fixed when it spawns, one per pool slot
typedef struct tfx_gpu_particle_cold_s {	//64 bytes
	tfxU32 uid;
	float inv_max_age;
	tfxU32 velocity_normal;					//Packed 10 bit xyz
	float base_velocity;
	float base_weight;
	tfx_float32x2_t base_size;
	float intensity_factor;
	float random_color;
	tfxU32 rotation_offsets;				//Packed 10 bit pitch, yaw, roll, or the roll angle's float bits when the emitter can't spin pitch and yaw
	float base_roll_spin;
	float base_pitch_spin;
	float base_yaw_spin;
	float image_start_frame;				//The frame at age 0, so the current frame is image_start_frame + image_frame_rate * age
	tfxU32 spawn_location;
	tfxU32 user_value;						//The spawn location's tfx_SetSpawnLocationUserData, unless a host's after_spawn changed it. 0 without a location
} tfx_gpu_particle_cold_t;

//A GPU simulated particle's noise state, fixed when it spawns, one per pool slot. Only rings with a noise force read it
typedef struct tfx_gpu_particle_noise_s {	//16 bytes
	tfx_float32x3_t offset;					//Noise field sample offset in authored units, zero for white noise
	float resolution;
} tfx_gpu_particle_noise_t;

//One particle spawned since the last fetch, which the kernel copies into its pool slot on the tick it arrives
typedef struct tfx_gpu_particle_spawn_s {	//112 bytes
	tfx_gpu_particle_hot_t hot;
	tfx_gpu_particle_cold_t cold;
	tfx_gpu_particle_noise_t noise;			//Zero for an emitter without noise
} tfx_gpu_particle_spawn_t;

//One emitter instance's ring at one update tick. A dispatch runs every tick record it points to in order
typedef struct tfx_gpu_particle_tick_s {	//64 bytes
	tfxU32 head;							//Ring slot of the oldest particle after this tick's bump
	tfxU32 count;							//Particles in the ring after this tick's spawns
	tfxU32 spawn_count;						//The last spawn_count of those were spawned this tick
	tfxU32 upload_offset;					//Into tfx_gpu_particle_frame_t spawns for this tick's spawns
	float frame_length;						//Milliseconds
	float update_time;						//Seconds
	float overall_scale;
	float velocity_adjuster;
	float drag_alpha;						//Fraction of the velocity drag removes this tick
	float acceleration_scale;				//Seconds of acceleration that survive the drag
	float drag_alpha_variation;
	float acceleration_scale_variation;
	tfxU32 force_offset;					//Row of this tick's first force in the same array as the ticks, see tfx_gpu_particle_frame_t force_rows
	tfxU32 force_count;						//Rows, so a fluctuating vortex counts twice
	tfxU32 flags;							//tfxGPUParticleTickFlags
	float location_age_offset;				//Milliseconds the ticks after this one in the fetch added to an unpaused location's age
} tfx_gpu_particle_tick_t;

typedef enum {
	tfxGPUParticleTickFlags_none = 0,
	tfxGPUParticleTickFlags_generated = 1 << 0,							//upload_offset is the row of this tick's tfx_gpu_spawn_block_t in the same array as the ticks, not spawn rows
} tfx_gpu_particle_tick_flag_bits;

typedef tfxU32 tfxGPUParticleTickFlags;			//tfx_gpu_particle_tick_flag_bits

//tfx_gpu_particle_force_t type_flags & tfxGPU_FORCE_TYPE_MASK
typedef enum {
	tfxGPUForceType_wind,					//Every wind force of the emitter summed into one flow in axis
	tfxGPUForceType_attract,
	tfxGPUForceType_vortex,
	tfxGPUForceType_shockwave,
	tfxGPUForceType_noise,
} tfx_gpu_particle_force_type;

typedef enum {
	tfxGPUForceFlags_none = 0,
	tfxGPUForceFlags_profile_bezier = 1 << 8,							//profile is a bezier segment, otherwise points on the straight line
	tfxGPUForceFlags_bounded = 1 << 9,									//Noise only: faded by the profile out to the falloff distance from origin
	tfxGPUForceFlags_fluctuating = 1 << 10,								//Vortex only: the next row is its tfx_gpu_particle_vortex_fluctuation_t
	tfxGPUForceFlags_at_location = 1 << 11,								//Evaluated in the particle's spawn location frame, origin already less the anchor. A shockwave is timed by the location's age
	tfxGPUForceFlags_location_delay = 1 << 12,							//Wind and noise: front holds a delay each particle waits out on its location's age. A delayed wind row adds to the flow
} tfx_gpu_particle_force_flag_bits;

#define tfxGPU_FORCE_TYPE_MASK 0xFF

//One emitter force resolved for one ring at one tick, in the particles' space. The rows follow the tick records in the same
//array, so a force row and a tick record must stay the same size
typedef struct tfx_gpu_particle_force_s {	//64 bytes
	tfxU32 type_flags;						//tfx_gpu_particle_force_type | tfx_gpu_particle_force_flag_bits
	float strength;							//Times the effect's noise multiplier, and the motion randomness base for white noise
	float falloff_scale;					//Reciprocal of the distance the profile spans, in the particles' units
	float front;							//Shockwave: how far its leading edge has travelled, at_location: where it is at location time 0. Other at_location rows: the delay
	tfx_float32x3_t origin;
	float radial_ratio;						//Vortex. at_location shockwave: the front's speed per millisecond of location time
	tfx_float32x3_t axis;					//Vortex: the unit axis. Wind: the summed flow
	float axial_ratio;						//Vortex. at_location shockwave: the radius it dies at
	tfx_float32x4_t profile;				//The profile graph's segment: from, curve1, curve2, to
} tfx_gpu_particle_force_t;

//The row after a fluctuating vortex's force row; after an at_location row harmonic_cos/harmonic_sin/turn_offset hold weight, turns per ms and turn at location time 0
typedef struct tfx_gpu_particle_vortex_fluctuation_s {	//64 bytes
	tfx_float32x3_t reference;				//Angle zero about the axis
	tfxU32 lobes;
	tfx_float32x3_t reference_turn;			//axis x reference
	float turn_offset_0;
	float harmonic_cos[3];					//Each harmonic's weight, current turn and frame average folded together
	float turn_offset_1;
	float harmonic_sin[3];
	float turn_offset_2;
} tfx_gpu_particle_vortex_fluctuation_t;

typedef enum {
	tfxGPUParticleDispatchFlags_none = 0,
	tfxGPUParticleDispatchFlags_capture = 1 << 0,						//The ring was written into the previous instance buffer, so its particles can interpolate
	tfxGPUParticleDispatchFlags_align_emission = 1 << 1,
	tfxGPUParticleDispatchFlags_align_emitter = 1 << 2,
	tfxGPUParticleDispatchFlags_random_color = 1 << 3,
	tfxGPUParticleDispatchFlags_lifetime_uniform_size = 1 << 4,
	tfxGPUParticleDispatchFlags_base_uniform_size = 1 << 5,
	tfxGPUParticleDispatchFlags_do_not_render = 1 << 6,
	tfxGPUParticleDispatchFlags_spin_roll = 1 << 7,
	tfxGPUParticleDispatchFlags_spin_3d = 1 << 8,
	tfxGPUParticleDispatchFlags_spin_pitch_and_yaw = 1 << 9,			//rotation_offsets holds packed pitch, yaw and roll
	tfxGPUParticleDispatchFlags_relative_angle = 1 << 10,
	tfxGPUParticleDispatchFlags_play_once = 1 << 11,
	tfxGPUParticleDispatchFlags_reverse_animation = 1 << 12,
	tfxGPUParticleDispatchFlags_drag_variation = 1 << 13,
	tfxGPUParticleDispatchFlags_relative_location = 1 << 14,			//Positions are local to the spawn location in the particle's cold spawn_location
	tfxGPUParticleDispatchFlags_location_rotation = 1 << 16,			//A spawn location has a rotation, applied after the emitter's
	tfxGPUParticleDispatchFlags_align_rotated_emission = 1 << 17,		//Emission alignment of a relative emitter, rotated by the emitter
	tfxGPUParticleDispatchFlags_noise = 1 << 18,						//The emitter has a noise force, its algorithm is in the noise_algorithm bits
	tfxGPUParticleDispatchFlags_noise_algorithm = 3 << 19,				//(flags >> tfxGPU_NOISE_ALGORITHM_SHIFT) & 3: white 0, simplex 1, curl 2, value 3
	tfxGPUParticleDispatchFlags_noise_anisotropic = 1 << 22,			//Split the noise against the particle's heading with the two noise gains
	tfxGPUParticleDispatchFlags_location_forces = 1 << 23,				//Some force rows are at_location, read the location table by the cold spawn_location. A relative ring's position is already in the frame
	tfxGPUParticleDispatchFlags_location_forces_rotation = 1 << 24,		//Those locations can have a rotation
	tfxGPUParticleDispatchFlags_hue = 1 << 25,							//The emitter has a hue graph or an effect hue, so the kernel samples it
} tfx_gpu_particle_dispatch_flag_bits;

#define tfxGPU_NOISE_ALGORITHM_SHIFT 19

typedef tfxU32 tfxGPUParticleDispatchFlags;      //tfx_gpu_particle_dispatch_flag_bits

typedef enum {
	tfxGPUParticleFlags_none = 0,
	tfxGPUParticleFlags_remove = 1 << 0,								//Expired, drawn with zero size until the ring head passes it
} tfx_gpu_particle_flag_bits;

//Bits 16 to 31 of a GPU particle's flags belong to the host's kernel hooks. The library zeroes them at spawn and never reads them
#define tfxGPU_PARTICLE_USER_FLAGS_SHIFT 16
#define tfxGPU_PARTICLE_USER_FLAGS_MASK 0xFFFF0000u

typedef tfxU32 tfxGPUParticleFlags;              //tfx_gpu_particle_flag_bits

//One emitter instance's ring for one tfx_GetGPUParticleFrame. Upload them verbatim and dispatch once per record
typedef struct tfx_gpu_particle_dispatch_s {	//128 bytes
	tfxU32 base;							//First pool slot of the ring
	tfxU32 capacity;						//Ring slots
	tfxU32 first_tick;						//Into tfx_gpu_particle_frame_t ticks
	tfxU32 tick_count;
	tfxU32 slot_start;						//Threads cover ring slots [slot_start, slot_start + slot_count) wrapped by capacity
	tfxU32 slot_count;
	tfxU32 flags;							//tfxGPUParticleDispatchFlags
	tfxU32 graph_offset;					//Into the graph lookups from tfx_GetGPUGraphLookupsBuffer
	tfxU32 property_index;					//gpu property index << 16, as tfx_instance_t indexes expects it, so global on a stage with global_library_indexes. The low 16 bits are the effect's hue in turns as snorm16 where 1 is half a turn
	tfxU32 instance_start;					//Ring position i on the last tick writes instance instance_start + i
	tfxU32 previous_instance_start;			//Where ring position 0 was written in the previous instance buffer
	tfxU32 bumped_since_write;				//Particles that left the ring since the previous write
	tfxU32 spawned_since_write;				//Particles that joined the ring since the previous write, they have nothing to interpolate from
	float global_intensity;
	float global_stretch;
	float image_frame_rate;					//Frames per second
	float end_frame;
	float emitter_roll;
	tfxU32 location_offset;					//Into tfx_gpu_particle_frame_t spawn_locations, the effect's table indexed by location slot
	tfxU32 location_count;
	tfx_float32x4_t emitter_rotation;		//Quaternion xyzw
	tfx_float32x3_t handle;					//Added to a relative position before the emitter rotation
	float noise_swerve_gain_minus_one;
	float noise_speed_gain_minus_swerve_gain;
	tfxU32 key;								//The host's key for choosing a kernel, see tfx_SetTemplateEmitterGPUKey and tfx_SetEffectGPUKey
	float location_delay;					//The emitter's spawn delay, which each location waits out from its own age
	float location_loop_length;				//The emitter's loop length, 0 if it doesn't loop
} tfx_gpu_particle_dispatch_t;

//A user spawn location as relative particles see it on the last tick of a fetch, one per location slot
typedef struct tfx_gpu_spawn_location_s {	//32 bytes
	tfx_float32x3_t position;
	tfxU32 local_id;						//slot | generation, matched against the cold spawn_location. 0xFFFFFFFF when the slot is empty
	tfxU64 packed_rotation;					//16 bit snorm quaternion, x|y in the low word and z|w in the high one
	float age;								//Milliseconds, which times a location's own forces
	tfxU32 paused;							//1 when the location's age stopped for the fetch
} tfx_gpu_spawn_location_t;

typedef enum {
	tfxGPUSpawnBlockFlags_none = 0,
	tfxGPUSpawnBlockFlags_relative = 1 << 0,							//Positions stay local to the spawn location, which the particle follows
	tfxGPUSpawnBlockFlags_random_color = 1 << 1,
	tfxGPUSpawnBlockFlags_random_start_frame = 1 << 2,
	tfxGPUSpawnBlockFlags_base_uniform_size = 1 << 3,					//One size draw for both width and height
	tfxGPUSpawnBlockFlags_spin_pitch_and_yaw = 1 << 4,					//rotation_offsets holds packed pitch, yaw and roll, otherwise the roll as a float
	tfxGPUSpawnBlockFlags_specify_pitch = 1 << 5,						//The angle offset is used as is rather than as the top of a random range
	tfxGPUSpawnBlockFlags_specify_yaw = 1 << 6,
	tfxGPUSpawnBlockFlags_specify_roll = 1 << 7,
	tfxGPUSpawnBlockFlags_relative_angle = 1 << 8,						//A rotated location's rotation is baked into the angle offsets of non relative particles
	tfxGPUSpawnBlockFlags_emission_cone = 1 << 9,
	tfxGPUSpawnBlockFlags_noise_offsets = 1 << 10,						//A field noise algorithm reads a random offset per particle
	tfxGPUSpawnBlockFlags_noise_resolution = 1 << 11,					//Any noise, motion randomness included
	tfxGPUSpawnBlockFlags_user = 1 << 12,								//A User emitter, whose spawns go through the host's ITfxParticleSpawn
	tfxGPUSpawnBlockFlags_fill_area = 1 << 13,
	tfxGPUSpawnBlockFlags_open_ends = 1 << 14,							//A hollow area without its two end faces
	tfxGPUSpawnBlockFlags_uniform_distribution = 1 << 15,				//Ellipse and disc spawns follow a low discrepancy sequence
	tfxGPUSpawnBlockFlags_shape = 1 << 16,								//Spawns in a shape rather than at a point, so the block's shape row is read
	tfxGPUSpawnBlockFlags_grid = 1 << 17,								//Spawns on the shape's grid points, so the block's step row is read
	tfxGPUSpawnBlockFlags_grid_random = 1 << 18,						//A random grid point rather than the next one along the march
	tfxGPUSpawnBlockFlags_grid_backwards = 1 << 19,						//The march runs down from grid_start instead of up
	tfxGPUSpawnBlockFlags_stepped_emission = 1 << 20,					//The emission direction steps through divisions, read from the step row
	tfxGPUSpawnBlockFlags_stepped_roll = 1 << 21,						//The roll steps through divisions, read from the step row
} tfx_gpu_spawn_block_flag_bits;

typedef tfxU32 tfxGPUSpawnBlockFlags;          //tfx_gpu_spawn_block_flag_bits

//The values a ring's generated spawns share on one tick, sampled once from the emitter's graphs. Spawn index i of the tick has uid
//tfx__seedgen_u32(uid_base + i) and draws every random value from that uid
typedef struct tfx_gpu_spawn_block_s {		//256 bytes
	tfxU32 uid_base;
	tfxU32 flags;							//tfxGPUSpawnBlockFlags
	tfxU32 run_row;							//The first tfx_gpu_spawn_run_t, in rows after this block's first row
	tfxU32 run_count;
	tfx_float32x3_t emission_direction;		//Before the cone, with the emitter rotation baked in for non relative particles
	float emission_cone_cosine;				//Cosine of half the emission range
	tfx_float32x3_t handle;					//The handle turned by the emitter, added to non relative positions
	float splatter;
	tfx_float32x3_t noise_step_offset;
	float noise_resolution;
	tfx_float32x3_t angle_offsets;			//Pitch, yaw, roll
	float animation_frames;
	tfx_float32x2_t size;
	tfx_float32x2_t size_variation;
	float frame_length;						//Milliseconds
	float inverse_frame_length;
	float overall_scale;
	float velocity;
	float life;
	float life_variation;
	float weight;
	float weight_variation;
	float roll_spin;
	float pitch_spin;
	float yaw_spin;
	float roll_spin_variation;
	float pitch_spin_variation;
	float yaw_spin_variation;
	float spawn_impulse;
	float spawn_impulse_variation;
	float noise_offset_variation;
	float noise_shared_offset;
	float start_frame;
	float arc_size;							//Disc and cylinder: the sweep around the y axis. Ellipse: the cosine of the bottom of its polar band
	float arc_offset;						//Where the sweep starts. Ellipse: the cosine of the top of its polar band
	tfxU32 spawn_ordinal;					//Particles the emitter spawned before this tick, which Both Ways alternates on
	tfxU32 padding[2];
	//Shapes, unused by Point and User
	tfx_float32x3_t emitter_size;
	tfxU32 shape;							//tfx_emission_type
	tfx_float32x4_t shape_rotation;			//x, y, z, w: the emitter's rotation, which turns non relative shape offsets
	tfx_float32x4_t emission_rotation;		//x, y, z, w: turns a shape's local emission direction by the pitch and yaw, then the emitter for non relative particles
	tfx_float32x3_t local_handle;			//The emitter's handle before its rotation
	tfxU32 emission_direction_mode;			//tfx_emission_direction
	//Grid points and stepped angles, only read with tfxGPUSpawnBlockFlags_grid or a stepped flag
	tfxU32 grid_start;						//Where the tick's march starts
	tfxU32 grid_points[3];					//x, y, z, each at least 1
	tfx_float32x3_t emission_fan_axis;		//The axis stepped emission turns about, turned like emission_direction
	float emission_range;					//Stepped emission's arc, and the tilt of a ring or cap
	float emission_step_jitter;				//0 to 1 of a division
	float roll_step_jitter;
	tfxU32 step_divisions;					//Emission divisions in the low 16 bits, roll divisions in the high 16
	tfxU32 step_settings;					//Emission distribution | roll distribution << 2 | emission tfxAngleStepFlags << 4 | roll tfxAngleStepFlags << 8
	tfxU32 step_seed;						//The emitter's seed index, which seeds the shuffled distribution
	tfxU32 step_padding[3];
} tfx_gpu_spawn_block_t;

typedef enum {
	tfxGPUSpawnRunFlags_none = 0,
	tfxGPUSpawnRunFlags_rotated = 1 << 0,								//The location has a rotation
	tfxGPUSpawnRunFlags_scale_life = 1 << 1,							//Its life multiplier isn't 1
} tfx_gpu_spawn_run_flag_bits;

typedef tfxU32 tfxGPUSpawnRunFlags;            //tfx_gpu_spawn_run_flag_bits

//One spawn location's share of a ring's generated spawns on one tick: spawn indexes [first_spawn, first_spawn + count)
typedef struct tfx_gpu_spawn_run_s {		//64 bytes
	tfxU32 local_id;						//slot | generation, which relative particles follow
	tfxU32 first_spawn;
	tfxU32 count;
	tfxU32 multipliers_row;					//The location's tfx_gpu_spawn_multipliers_t in rows after the block's first row, 0xFFFFFFFF when it has none
	tfx_float32x3_t position;
	tfxU32 flags;							//tfxGPUSpawnRunFlags
	tfx_float32x3_t captured_position;		//Where the location was at the start of the tick, spawns tween from it to position
	float location_roll;					//The location's turn about the roll axis, for relative angle particles that only roll
	tfxU64 packed_rotation;					//16 bit snorm quaternion, as tfx_gpu_spawn_location_t packs it
	tfxU32 user_data;						//tfx_SetSpawnLocationUserData
	tfxU32 spawn_ordinal;					//Spawns the location made for this emitter before the run, which stepped angles count from
} tfx_gpu_spawn_run_t;

//A spawn location's multipliers as its generated spawns apply them
typedef struct tfx_gpu_spawn_multipliers_s {	//64 bytes
	float particle_width;
	float particle_height;					//The width again when the effect has a uniform particle size
	float velocity;
	float inverse_life;
	float spin;
	float intensity;
	float splatter;
	float weight;
	float noise_offset;						//Added rather than multiplied
	float emitter_width;
	float emitter_height;
	float emitter_depth;
	tfxU32 padding[4];
} tfx_gpu_spawn_multipliers_t;

//Everything the stage simulated on the GPU since the last tfx_GetGPUParticleFrame
typedef struct tfx_gpu_particle_frame_s {
	const tfx_gpu_particle_dispatch_t *dispatches;
	const tfx_gpu_particle_tick_t *ticks;	//tick_count tick records then force_row_count rows of forces and generated spawn blocks: upload all of them as one array
	const tfx_gpu_particle_force_t *force_rows;		//ticks + tick_count, for a host that uploads the rows separately
	const tfx_gpu_particle_spawn_t *spawns;
	const tfx_gpu_spawn_location_t *spawn_locations;
	tfxU32 dispatch_count;
	tfxU32 tick_count;
	tfxU32 force_row_count;
	tfxU32 spawn_count;
	tfxU32 spawn_location_count;
	tfxU32 gpu_instance_start;				//The GPU instances follow the CPU ones, so this is the CPU instance count
	tfxU32 gpu_instance_count;
	tfxU32 gpu_layer_start[tfxLAYERS];
	tfxU32 gpu_layer_count[tfxLAYERS];
	//A layout for every array in one upload buffer, which tfx_CopyGPUParticleUpload fills. The dispatches start at 0, and every
	//section starts on a multiple of 256 bytes and of its row size, so it can also be bound as a view or range of its own
	size_t upload_bytes;
	size_t tick_upload_offset;				//The ticks with their force rows
	size_t location_upload_offset;
	size_t spawn_upload_offset;
	tfxU32 tick_base;						//The section offsets in rows, for tfx_gpu_particle_push_t
	tfxU32 location_base;
	tfxU32 spawn_base;
} tfx_gpu_particle_frame_t;

//Push constants for the reference kernel, simulate_particles in Shaders/timelinefx.slang. Buffers are bindless indexes, and each *_base is
//the row its section starts at in that buffer: 0 for a buffer or view of its own, the section's row when sections share one
typedef struct tfx_gpu_particle_push_s {	//64 bytes
	tfxU32 dispatch_index;					//First, so only these 4 bytes change from one dispatch to the next
	tfxU32 hot_buffer_index;				//tfx_gpu_particle_hot_t per pool slot from row 0, read and written
	tfxU32 cold_buffer_index;				//tfx_gpu_particle_cold_t per pool slot, written when a particle is ingested
	tfxU32 noise_buffer_index;				//tfx_gpu_particle_noise_t per pool slot, only touched by rings with noise
	tfxU32 instance_buffer_index;			//The instance buffer the CPU instances were uploaded to
	tfxU32 graph_buffer_index;				//tfx_GetGPUGraphLookupsBuffer
	tfxU32 dispatch_buffer_index;			//tfx_gpu_particle_frame_t dispatches from row 0
	tfxU32 tick_buffer_index;				//tfx_gpu_particle_frame_t ticks, with the force rows after them
	tfxU32 spawn_buffer_index;				//tfx_gpu_particle_frame_t spawns
	tfxU32 location_buffer_index;			//tfx_gpu_particle_frame_t spawn_locations
	tfxU32 cold_base;
	tfxU32 noise_base;
	tfxU32 tick_base;
	tfxU32 spawn_base;
	tfxU32 location_base;
	tfxU32 padding;							//Keeps a host's own push constants after these on a 16 byte offset
} tfx_gpu_particle_push_t;

typedef struct tfx_gpu_particle_stats_s {
	tfxU32 capacity;						//Pool slots, tfx_stage_info_t gpu_particle_capacity
	tfxU32 allocated;						//Pool slots held by rings
	tfxU32 slots_in_use;					//Ring slots in use. Includes particles that died but whose batch hasn't expired, which are still drawn as zero size billboards
	tfxU32 estimated_live_particles;		//Particles still alive, assuming each spawn's particles die evenly over its range of lives. Compare this with tfx_GetParticleCount, which reads higher by the dead particles still holding their place
	tfxU32 ring_count;
	tfxU32 blocked_emitters;				//GPU emitters that got no ring because the pool was full, raise gpu_particle_capacity
	tfxU64 dropped_spawns;					//Spawns refused because a ring was full, since the stage was created or last cleared. For overlapping spawn locations size the template with tfx_SetTemplateGPUCapacityForLocations, for slow or catch-up updates set tfx_SetGPUSizingUpdateRate, otherwise raise the ring with tfx_SetTemplateEmitterGPUCapacity. Blocked emitters aren't counted here
} tfx_gpu_particle_stats_t;

//A stage's health in one place, see tfx_GetStageDiagnostics and tfx_DescribeStage
typedef struct tfx_stage_diagnostics_s {
	tfxU32 effect_count;
	tfxU32 emitter_count;
	tfxU32 cpu_particles;					//tfx_GetParticleCount
	tfx_gpu_particle_stats_t gpu;
	tfxU32 effects_without_locations;		//Effects that spawn at user spawn locations and have none right now, which can be deliberate
	tfxU32 unfetched_gpu_updates;			//Updates that simulated GPU particles since the last tfx_GetGPUParticleFrame
	tfxErrorFlags library_error_flags;		//Every library a running effect came from, combined, without tfxErrorCode_shapes_changed_on_disk
	bool gpu_simulation_on;					//tfx_StageGPUSimulationIsOn
} tfx_stage_diagnostics_t;

//The particle pool's hot, cold and noise sections. Each offset is a multiple of 256 bytes, so a host can put all three in one
//buffer and bind each as its own view or range, or give each a buffer of its own
typedef struct tfx_gpu_particle_buffer_sizes_s {
	size_t pool_bytes;						//All three sections in one buffer
	size_t hot_bytes;						//Hot starts at offset 0
	size_t cold_offset;
	size_t cold_bytes;
	size_t noise_offset;
	size_t noise_bytes;
	tfxU32 cold_base;						//cold_offset in cold rows, for tfx_gpu_particle_push_t when cold shares the pool buffer
	tfxU32 noise_base;						//noise_offset in noise rows
} tfx_gpu_particle_buffer_sizes_t;

typedef struct tfx_ribbon_buffer_requirements_s {
	tfxU32 segment_buffer_size_in_bytes;
	tfxU32 ribbon_buffer_size_in_bytes;
	tfxU32 emitter_buffer_size_in_bytes;
} tfx_ribbon_buffer_requirements_t;

//Element offsets (not bytes) of one stage's data within a ribbon batch, or running totals for the whole batch
typedef struct tfx_ribbon_batch_offsets_s {
	tfxU32 segment_offset;
	tfxU32 ribbon_offset;
	tfxU32 emitter_offset;
	tfxU32 vertex_offset;
	tfxU32 index_offset;
} tfx_ribbon_batch_offsets_t;

//Packs the ribbons of several stages into one shared set of buffers. See tfx_BeginRibbonBatch.
typedef struct tfx_ribbon_batch_s {
	void *segments_dst;
	void *ribbons_dst;
	void *emitters_dst;
	size_t segments_capacity;
	size_t ribbons_capacity;
	size_t emitters_capacity;
	tfx_ribbon_batch_offsets_t totals;
	tfxU32 stage_count;
} tfx_ribbon_batch_t;

//Bytes per ribbon the baked ribbon playback writes, for sizing its output buffer from tfx_impl_baked_frame_t ribbon_count
#define tfxRIBBON_SIZE 48

typedef struct tfx_ribbon_batch_sizes_s {
	size_t segment_buffer_size_in_bytes;
	size_t ribbon_buffer_size_in_bytes;
	size_t emitter_buffer_size_in_bytes;
	size_t vertex_buffer_size_in_bytes;
	size_t index_buffer_size_in_bytes;
} tfx_ribbon_batch_sizes_t;

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
struct tfx_ribbon_bucket_globals_s  {
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
	tfxU32 emitter_offset;			//Where this stage's emitters start in a ribbon batch's emitter buffer
	tfxU32 padding;					//Uniform blocks are multiples of 16 bytes
};

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
	tfxU32 padding[2];
} tfx_sprite_data_push_t;

//Push constants for tfx_ribbon_playback
typedef struct tfx_ribbon_data_push_s {
	tfxU32 animation_instances_total;
	tfxU32 ribbons_total;
	tfxU32 offsets_index;
	tfxU32 animation_instances_index;
	tfxU32 ribbon_data_index;
	tfxU32 ribbon_output_index;
	tfxU32 ribbon_properties_index;
	tfxU32 padding;
} tfx_ribbon_data_push_t;

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
typedef void(*tfx_shape_remover)(tfx_image_data_t *image_data, void *user_data);
typedef void(*tfx_uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset);
typedef bool(*tfx_maybe_render_instance_callback)(tfx_animation_manager animation_manager, tfx_float32x3_t position, float radius, void *user_data);
typedef void(*tfx_print_callback)(const char *message, void *user_data);

//One location's share, or part of it, of the particles a User emitter spawned this tick, see tfx_SetStageSpawnCallback
typedef struct tfx_spawn_callback_info_s {
	tfxEffectID effect;
	tfx_effect_descriptor emitter;
	tfxSpawnLocationID location;			//tfxINVALID_SPAWN_LOCATION when the effect doesn't spawn at user locations
	const tfxU32 *uids;						//One per particle, a stable seed for deterministic draws
	tfxU32 user_data;						//tfx_SetSpawnLocationUserData, 0 without a location
	tfxU32 index_in_run;					//Where the first of these particles is among the location's spawns this tick
	tfxU32 run_count;						//Particles the location spawned this tick
	tfxU32 count;
	bool relative;							//Positions are local to the location rather than world space
} tfx_spawn_callback_info_t;

typedef void(*tfx_spawn_callback)(const tfx_spawn_callback_info_t *info, float *position_x, float *position_y, float *position_z, void *user_data);


tfxAPI void tfx_UpdateAnimationManagerBufferMetrics(tfx_animation_manager animation_manager);
tfxAPI float tfx_DegreesToRadians(float degrees);
tfxAPI float tfx_RadiansToDegrees(float radians);

//--------------------------------
//Random numbers
//--------------------------------
//A xorshift128+ generator seeded as tfx_RandomReseed does, so any seed, 0 included, gives a well spread stream
tfxAPI tfx_random_t tfx_CreateRandom(tfxU32 seed);
//Step the generator once and throw the number away
tfxAPI void tfx_AdvanceRandom(tfx_random_t *random);
//Reseed from the millisecond clock, as tfx_RandomReseed does
tfxAPI void tfx_RandomReseedTime(tfx_random_t *random);
//Set the two state words as they are, unmixed. They mustn't both be 0, and small values give poor first draws, so pass well mixed values or use tfx_RandomReseed
tfxAPI void tfx_RandomReseed2(tfx_random_t *random, tfxU64 seed1, tfxU64 seed2);
//Reseed from any value, 0 included: the seed is mixed into both state words with split mix 64, and the same seed always gives the same stream
tfxAPI void tfx_RandomReseed(tfx_random_t *random, tfxU64 seed);
//A float from 0 to 1, which can round up to exactly 1
tfxAPI float tfx_GenerateRandom(tfx_random_t *random);
//A float from 0 to max
tfxAPI float tfx_RandomRangeZeroToMax(tfx_random_t *random, float max);
//A float between from and to, worked out as to - (to - from) * tfx_GenerateRandom, so low draws land near to rather than from
tfxAPI float tfx_RandomRangeFromTo(tfx_random_t *random, float from, float to);
//An int from from to to, both included
tfxAPI int tfx_RandomRangeFromToInt(tfx_random_t *random, int from, int to);
//An unsigned int from 0 to max - 1
tfxAPI tfxU32 tfx_RandomRangeZeroToMaxUInt(tfx_random_t *random, tfxU32 max);
//An int from 0 to max - 1
tfxAPI int tfx_RandomRangeZeroToMaxInt(tfx_random_t *random, int max);
//Move the generator onto a different stream by amount (the first state word is multiplied by it, the second has it added), to branch a copy
tfxAPI void tfx_AlterRandomSeedU64(tfx_random_t *random, tfxU64 amount);
//As tfx_AlterRandomSeedU64
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

    // 2. Load an effects library from a .tfx file. shape_loader hands you each image; uv_lookup returns its uv rect.
    //    Then pack the shapes (tfx_impl_BuildShapeAtlas does it), call tfx_UpdateLibraryGPUImageData so uv_lookup runs
    //    again, and tfx_UpdateGlobalGPUData before drawing more than one library with one set of buffers. The order is
    //    in implementations/README.md, "Setting up a library".
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
    tfx_FreeStage(stage);
    tfx_FreeLibrary(library);           // also frees the effect templates made from it
    tfx_EndTimelineFX();

Ownership and teardown ordering
-------------------------------
Everything is allocated from the memory pool created by tfx_BeginTimelineFX, so tfx_EndTimelineFX must be the
very last TimelineFX call you make, because it tears down the worker threads and the pool itself.

tfx_EndTimelineFX sweeps up anything you did not free yourself: it completes and frees every stage, then frees every
animation manager, then every library (and with it every effect template made from it), then any streams, packages and
GPU shape lists still allocated. It prints one tip listing what it freed for you. That isn't a leak, it's there so you know. Freeing any of them yourself
first is fine, they deregister from the sweep as they go, so there is no double free either way. tfx_GetMemoryInUse
tells you how much is still allocated at any point.

If you free objects yourself rather than leaving them to the sweep, either order of stage and library is fine:
tfx_FreeLibrary removes the library's effects from every stage before it frees anything, so you can unload one library
while effects from others keep playing. An effect template holds its own detached clone of the effect, but
tfx_FreeLibrary frees every template made from the library, so drop your handles to them and don't free them again. A stage has outstanding worker-thread work
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
Cleanup up all threads and memory used by timelinefx. Anything you didn't free is freed for you, and a tip lists it
*/
tfxAPI void tfx_EndTimelineFX(void);

/*
Get how much of TimelineFX's memory pools is allocated, across every pool. Call it when no stage is updating. It won't
go back to the figure after tfx_BeginTimelineFX once a library has been loaded, because some tables shared by every library
(GPU graph lookups, graph indexes and so on) are kept until tfx_EndTimelineFX. To look for leaks, compare calls taken at
the same point in a load/free cycle: growth from one cycle to the next is a leak. tfx_EndTimelineFX lists anything you
didn't free.
* @param bytes                    Set to the bytes in use, can be NULL
* @param blocks                   Set to the number of allocations, can be NULL
*/
tfxAPI void tfx_GetMemoryInUse(size_t *bytes, tfxU32 *blocks);

/*
Send every message TimelineFX prints to your own log instead of stdout. Errors start "TimelineFX error:" and are always on:
they only print when something has failed, such as running out of memory pools or calling tfx_BeginTimelineFX twice.
Tips are one line messages starting "TimelineFX tip:" that point out a likely mistake once, at the moment the library sees
it, for example an effect that spawns at user spawn locations being added to a stage. Define TFX_DISABLE_TIPS to compile
them out. In debug builds a failed assert prints "TimelineFX assert failed:" with its condition, file and line, then
stops the program with the C runtime's assert, which in a Windows GUI program is a dialog that waits for a click. Define
your own TFX_ASSERT when you compile timelinefx.cpp to change that.
Can be set before or after tfx_BeginTimelineFX. It can be called from the stage update thread and worker threads, so it
must be thread safe.
* @param callback           Receives each message without a trailing newline, or NULL to go back to printing to stdout
* @param user_data          Passed to the callback
*/
tfxAPI void tfx_SetPrintCallback(tfx_print_callback callback, void *user_data);

/*
Print a message through the callback set with tfx_SetPrintCallback, or to stdout if there isn't one. Backends and helper
layers use it so their messages land in the same log as the library's. Messages longer than 1023 characters are cut short.
* @param format             A printf format string, without a trailing newline
*/
#if defined(__MINGW32__) && !defined(__clang__)
tfxAPI void tfx_Print(const char *format, ...) __attribute__((format(gnu_printf, 1, 2)));
#elif defined(__GNUC__) || defined(__clang__)
tfxAPI void tfx_Print(const char *format, ...) __attribute__((format(printf, 1, 2)));
#else
tfxAPI void tfx_Print(const char *format, ...);
#endif

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
must be registered again by the caller with tfx_SetTemplateEffectUpdateCallback, tfx_SetLibraryUVLookup and
tfx_SetAnimationManagerInstanceCallback.
* @returns    True if the context was adopted
*/
#ifdef __cplusplus
tfxAPI bool tfx_SetContext(tfx_context context, const tfx_allocation_callbacks_t *allocation_callbacks = nullptr);
#else
tfxAPI bool tfx_SetContext(tfx_context context, const tfx_allocation_callbacks_t *allocation_callbacks);
#endif

/*
Call this function to wait on all thread work currently in process in TimelineFX. This does not close threads.
Use tfx_SuspendTimelineFX to draid and close threads.
*/
tfxAPI void tfx_DrainAllThreadWork(void);

/*
Drain and stop all current TimelineFX asynchronous work without deallocating memory.
*/
tfxAPI void tfx_SuspendTimelineFX(void);

/*
Restart the worker threads after a call to tfx_SuspendTimelineFX.
*/
tfxAPI void tfx_ResumeTimelineFX(void);

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
* @param uv_lookup        A function that fills in a shape frame's uvs for the GPU, called for every frame of every shape when the GPU shape data is built:
*                         void uv_lookup(void *ptr, tfx_gpu_image_data_t *image_data, int offset), where ptr is what your shape_loader set with
*                         tfx_SetImagePointer and offset is the animation frame. Called again by tfx_UpdateLibraryGPUImageData.
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

//What a refresh found. Effects and shapes are reported as counts, because a host acts on the library as a whole
typedef enum {
	tfxRefreshFlags_none = 0,
	tfxRefreshFlags_effects_changed = 1 << 0,   //Effects with a newer version on disk were overwritten, and anything running from them restarted
	tfxRefreshFlags_effects_added = 1 << 1,     //Effects that were only on disk were added to the library
	tfxRefreshFlags_effects_removed = 1 << 2,   //Effects no longer on disk were marked for deletion and anything running from them expired
	tfxRefreshFlags_shapes_changed = 1 << 3,    //The shape set differs, so an atlas rebuild is needed if used
	tfxRefreshFlags_unreadable = 1 << 4,        //The file could not be re-read; nothing was applied
	tfxRefreshFlags_merged = 1 << 5,            //The changed values were applied to the library in place
	tfxRefreshFlags_needs_reload = 1 << 6,      //A new shape could not be loaded, so reload the library to get it
	//Which of the library's gpu-side buffers the refresh rebuilt. A host owns the uploaded copy of each of
	//these, so it has to re-upload the ones flagged here or it goes on drawing with the old data.
	tfxRefreshFlags_gpu_shapes_changed = 1 << 7,           //tfx_GetLibraryGPUShapes was rebuilt: re-upload it. Any cached tfx_GetGPUShapesArray pointer is now dangling
	tfxRefreshFlags_particle_properties_changed = 1 << 8,  //tfx_GetParticlePropertiesBuffer changed, and may have grown: re-upload it
	tfxRefreshFlags_color_ramps_changed = 1 << 9,          //The color ramp bitmaps were rebuilt and re-indexed: re-upload the ramp texture array
	tfxRefreshFlags_graph_lookups_changed = 1 << 10,       //tfx_GetGPUGraphLookupsBuffer changed, and may have resized: re-upload it
	//Reported by tfx_RefreshShape only
	tfxRefreshFlags_shape_not_found = 1 << 11,             //No shape has that name or file name; nothing was applied
	tfxRefreshFlags_shape_matches_another = 1 << 12,       //The new image is identical to a different shape; nothing was applied
	tfxRefreshFlags_shape_layout_mismatch = 1 << 13,       //The image no longer fits the shape's frame layout, or can't be laid out as asked; nothing was applied
	tfxRefreshFlags_shape_layout_changed = 1 << 14,        //The shape took a new frame size or frame count
} tfx_refresh_flag_bits;

typedef tfxU32 tfxRefreshFlags;                 //tfx_refresh_flag_bits

typedef enum {
	tfxRefreshShapeFlags_none = 0,
	tfxRefreshShapeFlags_allow_layout_change = 1 << 0,     //Accept an image of a different size. A sheet keeps its grid, so frames scale with it
} tfx_refresh_shape_flag_bits;

typedef tfxU32 tfxRefreshShapeFlags;            //tfx_refresh_shape_flag_bits

typedef enum {
	tfxEffectTemplateFlags_none = 0,
	tfxEffectTemplateFlags_marked_for_deletion					= 1 << 0,
	tfxEffectTemplateFlags_needs_updating						= 1 << 1,
	tfxEffectTemplateFlags_deletion_reported					= 1 << 2,		//The block on adding it to a stage has been reported once
	tfxEffectTemplateFlags_user_location_tip_checked			= 1 << 3,
} tfx_effect_template_flag_bits;

typedef tfxU32 tfxEffectTemplateFlags;          //tfx_effect_template_flag_bits

typedef enum {
	tfxSpawnLocationAdd_none = 0,
	tfxSpawnLocationAdd_transient				= 1 << 0,		//Removed automatically after the next update, for one off bursts
	tfxSpawnLocationAdd_no_auto_remove			= 1 << 1,		//Keep the location after a finite effect has played out, instead of removing it. A kept location never lets tfxSpawnLocationAdd_expire_effect fire
	tfxSpawnLocationAdd_expire_effect			= 1 << 2,		//Soft expire the effect once it has no locations left, so one effect per burst ends like any finite effect
} tfx_spawn_location_add_flag_bits;

typedef tfxU32 tfxSpawnLocationAddFlags;        //tfx_spawn_location_add_flag_bits

typedef enum {
	tfxUpdateStepFlags_none = 0,
	tfxUpdateStepFlags_write_instances			= 1 << 0,		//Write the instance buffer this step. Set it on the last step before you draw
} tfx_update_step_flag_bits;

typedef tfxU32 tfxUpdateStepFlags;              //tfx_update_step_flag_bits

typedef struct tfx_refresh_result_s {
	tfxRefreshFlags flags;               //A combination of tfxRefreshFlags
	tfxU32 library_version;              //The library_version found on disk, or the loaded one if unreadable
	tfxU32 changed_count;                //Effects overwritten from disk
	tfxU32 added_count;                  //Effects added from disk
	tfxU32 removed_count;                //Effects newly missing from disk, counted by the first refresh that finds them gone
	tfxU32 added_shape_count;            //Shapes passed to shape_loader and now in the library
	tfxU32 removed_shape_count;          //Shapes passed to shape_remover and already gone from the library
} tfx_refresh_result_t;

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
	tfxErrorCode_some_images_loaded_without_user_ptr
	tfxErrorCode_shapes_changed_on_disk				Not a failure: an image was edited outside the editor and its emitters were remapped to it
	tfxErrorCode_file_version_out_of_date			Not a failure: the library loaded fine, re-save it in the editor to update its format
*/
tfxAPI tfxErrorFlags tfx_GetLibraryErrorStatus(tfx_library library);

/*
What a library contains, as text: its error flags, then a line for each effect (how long it plays and whether it spawns at
user spawn locations) followed by a line for each of its emitters (layer, emission type, which emitter an Other Emitter
spawns from, how many particles a single shot emitter spawns, whether its particles are relative, whether it draws nothing
because its width or intensity is always 0, usually a carrier other emitters spawn from, and whether it runs on the GPU,
with its ring capacity, or on the CPU and why). Emitters are named by the path the tfx_SetTemplate functions take.
* @param library                   A handle to a loaded library
* @param buffer                    Filled with as many whole lines as fit, always null terminated. Can be NULL with a size of 0
* @param buffer_size               The size of buffer in bytes
* @returns tfxU32                  The buffer size the whole text needs, terminator included
*/
tfxAPI tfxU32 tfx_DescribeLibrary(tfx_library library, char *buffer, tfxU32 buffer_size);

/*
tfx_DescribeLibrary printed a line at a time, through tfx_SetPrintCallback if one is set.
* @param library                   A handle to a loaded library
*/
tfxAPI void tfx_PrintLibrary(tfx_library library);

/*
Get the file format version that the library was loaded from. A library that was created rather than
loaded, or one whose file had no version line, reports the current format version and 0 respectively.
* @param library		A handle to the library
* @returns tfxU32		The file version found in the file
*/
tfxAPI tfxU32 tfx_GetLibraryFileVersion(tfx_library library);

typedef enum {
	tfxLayerFlags_none = 0,
	tfxLayerFlags_writes_depth = 1 << 0,		//Draw with tfx_billboard_cutout_fragment and tfx_ribbon_cutout_fragment, blending off and depth write on, before the layers that don't
	tfxLayerFlags_lit = 1 << 1,					//Draw with the lit vertex and fragment entry points, which light the particles with the view uniforms' tfx_view_light_t
} tfx_layer_flag_bits;

typedef tfxU32 tfxLayerFlags;					//tfx_layer_flag_bits

/*
Get the number of layers a library defines, which is at most tfxLAYERS. Loop your layers up to this count. Deleting a layer in the editor
leaves a gap rather than moving the layers after it, so a gap has an empty name, no flags and no instances, and it still counts.
* @param library		A handle to the library
* @returns tfxU32		The library's layer count, at least 2
*/
tfxAPI tfxU32 tfx_GetLayerCount(tfx_library library);

/*
Get the flags of one of a library's layers. Each emitter draws on one of the library's layers, and the instance buffer is grouped by layer
(tfx_GetInstanceBufferByLayer). Ribbon emitters have a layer too, which each ribbon dispatch reports in its layer field. Layer 0, Default,
blends and layer 1, Depth, writes depth in every library. The other layers are set per library in the editor. Draw the layers that write
depth first, in index order, billboards then ribbons, then the rest, and give those draws your scene's depth buffer so the scene hides the
particles behind it. When a stage holds effects from more than one library, give them the same layer table (the editor's Copy Layers does
that), because the stage has one range per layer index whichever library an effect came from. Layer names are unique within a library, so a
host that gives a custom layer its own pipeline can find its index once at load by looping over tfx_GetLayerName, rather than hard coding it.
A refresh can change the flags, so read them when you draw. A layer can also have a particle limit, set in the editor: once that many CPU
particles are alive on the layer across the stage, nothing new spawns on it until some die.

A layer with tfxLayerFlags_lit draws with tfx_billboard_lit_vertex, tfx_billboard_baked_lit_vertex and tfx_ribbon_lit_vertex, and with
tfx_billboard_lit_fragment and tfx_ribbon_lit_fragment, or their cutouts when the layer writes depth as well. Pipeline state is the same as
the unlit entry points'. The lit fragments light the color with the tfx_view_light_t your view uniforms hold at tfxVIEW_LIGHT_OFFSET, so
fill it every frame: a zeroed light draws lit particles black. A host that draws a lit layer with the unlit entry points gets it unlit.
* @param library		A handle to the library
* @param layer			The layer, below tfx_GetLayerCount
* @returns tfxLayerFlags	The layer's tfx_layer_flag_bits
*/
tfxAPI tfxLayerFlags tfx_GetLayerFlags(tfx_library library, tfxU32 layer);

/*
Get the name a library's layer was given in the editor.
* @param library		A handle to the library
* @param layer			The layer, below tfx_GetLayerCount
* @returns const char*	The layer's name, empty for a deleted layer
*/
tfxAPI const char *tfx_GetLayerName(tfx_library library, tfxU32 layer);

/*
Checks the library on disk to see if it's been updated since it was loaded.

Effects with a newer version on disk are overwritten in place and anything running from them in a stage is
restarted. Effects that are new on disk are added, and effects that have gone from disk are marked for deletion,
with anything running from them expired. Templates of a changed effect are updated and templates of a removed one
are marked for deletion, so tfx_AddEffectTemplateToStage refuses them. Re-resolve any effect descriptor you hold
after a refresh that reports effects removed.

Shapes that are added are passed to shape_loader for you to load in the renderer, and shapes that are removed
are passed to shape_remover.

The library's own memory is all this updates. Anything you uploaded to the gpu is your copy and goes stale,
so check the tfxRefreshFlags_*_changed bits in the result and re-upload what they name: the gpu shape data,
the particle properties buffer, the color ramp bitmaps and the global graph lookup buffer.

The library only has to be reloaded if the file didn't parse or tfxRefreshFlags_needs_reload is set.

It's safe to call this function regularily as it only loads the initial bytes of the file (whether it's
loaded from a folder or a .tfx package) to assertain whether the library is a newer version.

* @param library                 A handle to the loaded library to compare against its file
* @param shape_loader            Called for each shape that is new on disk, the same as when the library is loaded
* @param uv_lookup               Replaces the library's uv lookup if not null, used to rebuild the gpu shape data
* @param shape_remover           Called with a copy of each shape that is gone from disk, after it is removed from the library. Can be null
* @param user_data               Passed to shape_loader and shape_remover
* @param result                  Filled in with the flags and counts described above
*/
tfxAPI void tfx_RefreshLibrary(tfx_library library, tfx_shape_loader shape_loader, tfx_uv_lookup uv_lookup, tfx_shape_remover shape_remover, void *user_data, tfx_refresh_result_t *result);

/*
Reload one shape's image from the library's own folder or package, for when the image file was changed
directly and nothing else in the library was. Only this shape is touched; the rest of the library and
anything running in a stage carry on as they were.

shape_loader is called with the shape's image data, its ptr field cleared, the same as when the library was
loaded, and then shape_remover with a copy of the image data as it was, so you can free the old texture.
The shape keeps its place in the library, so image pointers you hold stay valid, but its image_hash
changes to the hash of the new image.

An image of a different size is refused unless tfxRefreshShapeFlags_allow_layout_change is set. With it,
a single image takes the new size and an animation keeps its grid of frames, so halving the sheet halves
each frame. Pass frame_count to change how many frames the sheet holds: the frame size is kept if the new
sheet still fits that many, otherwise the number of columns is. A frame count change restarts anything
running that uses the shape.

Check the returned flags: tfxRefreshFlags_gpu_shapes_changed and tfxRefreshFlags_particle_properties_changed
name what to re-upload, the same as tfx_RefreshLibrary. When the shape is not refreshed the flags say why:
tfxRefreshFlags_shape_not_found, tfxRefreshFlags_shape_matches_another, tfxRefreshFlags_shape_layout_mismatch
or tfxRefreshFlags_unreadable, which a library loaded from memory always reports. tfxRefreshFlags_none means
the image on disk is the one already loaded.

Animation managers keep their own copy of each shape, so sprite data recorded with the old image is not updated.

* @param library                 A handle to a library loaded from a file or folder
* @param shape_name              The shape's name or its file name. Any directory in front of it is ignored
* @param flags                   A combination of tfxRefreshShapeFlags
* @param frame_count             The number of frames the image now holds, or 0 to keep the current count
* @param shape_loader            Called to load the new image, the same as when the library is loaded
* @param shape_remover           Called with a copy of the shape as it was, after the new image is loaded. Can be null
* @param user_data               Passed to shape_loader and shape_remover
* @return tfxRefreshFlags        tfxRefreshFlags_shapes_changed when the shape was refreshed
*/
tfxAPI tfxRefreshFlags tfx_RefreshShape(tfx_library library, const char *shape_name, tfxRefreshShapeFlags flags, tfxU32 frame_count, tfx_shape_loader shape_loader, tfx_shape_remover shape_remover, void *user_data);

/*
The same as tfx_RefreshShape but with the image supplied by you, for a library loaded from memory or an image
that comes from your own asset pipeline.
* @param data                    The encoded image, a .png or .ktx2 the same as the library holds
* @param size                    The size of data in bytes
*/
tfxAPI tfxRefreshFlags tfx_RefreshShapeFromMemory(tfx_library library, const char *shape_name, const void *data, tfxU32 size, tfxRefreshShapeFlags flags, tfxU32 frame_count, tfx_shape_loader shape_loader, tfx_shape_remover shape_remover, void *user_data);

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
* Set the uv lookup callback that tfx_UpdateLibraryGPUImageData uses. Normally this is set for you when the library is loaded, but it is cleared by
* tfx_SetContext so it must be set again after hot reloading the module that owns the callback.
* @param library                   A valid pointer to a tfx_library
* @param uv_lookup                 A function pointer that fills in the uv coordinates from whatever renderer you're using
*/
tfxAPI void tfx_SetLibraryUVLookup(tfx_library library, tfx_uv_lookup uv_lookup);

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
tfxAPI void *tfx_GetImagePointer(tfx_image_data_t *image);
tfxAPI void tfx_SetGPUImageTextureInfo(tfx_gpu_image_data_t *image, float x, float y, float z, float w, int array_index);
tfxAPI int tfx_GetImageFrameCount(tfx_image_data_t *image);
tfxAPI int tfx_GetImageWidth(tfx_image_data_t *image);
tfxAPI int tfx_GetImageHeight(tfx_image_data_t *image);
tfxAPI tfx_image_format tfx_GetImageFormat(tfx_image_data_t *image);
tfxAPI void* tfx_GetBitmapData(tfx_bitmap_t *bitmap);
tfxAPI size_t tfx_GetBitmapSize(tfx_bitmap_t *bitmap);
tfxAPI int tfx_GetBitmapWidth(tfx_bitmap_t *bitmap);
tfxAPI int tfx_GetBitmapHeight(tfx_bitmap_t *bitmap);

/*
Print all the effect names in a library through the print callback, see tfx_SetPrintCallback
* @param tfx_library                A valid pointer to a tfx_library
*/
tfxAPI void tfx_ListEffectNames(tfx_library library);

/*
Get a top level descriptor in the library by it's index, which will be an effect or a folder (check with tfx_GetDescriptorType). Use with
tfx_GetLibraryRootDescriptorCount to loop over the library, and tfx_GetDescriptorChildCount/tfx_GetDescriptorChild to get the effects in a folder
or the emitters in an effect. To get a descriptor by it's path use tfx_GetLibraryDescriptor instead.
* @param tfx_library                A valid pointer to a tfx_library
* @param index                      Index of the descriptor, must be less than tfx_GetLibraryRootDescriptorCount
*/
tfxAPI tfx_effect_descriptor tfx_GetLibraryRootDescriptor(tfx_library library, tfxU32 index);

/*
Get a folder, effect or emitter in the library by it's path. So for example, if you want to get the emitter "spark" in effect "explosion" then you could do
tfx_GetLibraryDescriptor(library, "explosion/spark"). You will need this function to apply user data and update callbacks to effects and emitters before
adding the effect to the stage. Returns NULL if the path doesn't exist.
* @param tfx_library_t                A valid pointer to a tfx_library_t
* @param const char *path             Path to the folder, effect or emitter
*/
tfxAPI tfx_effect_descriptor tfx_GetLibraryDescriptor(tfx_library library, const char *path);

/*
Check whether a path resolves to a folder, effect or emitter in the library.
* @param library                A valid pointer to a tfx_library_t
* @param path                   Path to the folder, effect or emitter
* @returns                      true if the path exists in the library
*/
tfxAPI bool tfx_IsValidDescriptorPath(tfx_library library, const char *path);

/*
Get the number of top level descriptors in a library, which are effects and folders. Loop over them with tfx_GetLibraryRootDescriptor, for example to build a tree view:
	for (tfxU32 index = 0; index != tfx_GetLibraryRootDescriptorCount(library); ++index) {
		tfx_effect_descriptor descriptor = tfx_GetLibraryRootDescriptor(library, index);
		if (tfx_GetDescriptorType(descriptor) == tfxFolder) {
			for (tfxU32 child_index = 0; child_index != tfx_GetDescriptorChildCount(descriptor); ++child_index) {
				tfx_effect_descriptor effect = tfx_GetDescriptorChild(descriptor, child_index);
			}
		}
	}
* @param library                A valid pointer to a tfx_library_t
*/
tfxAPI tfxU32 tfx_GetLibraryRootDescriptorCount(tfx_library library);

/*
Get the type of an effect descriptor: tfxEffectType, tfxEmitterType (particle emitter), tfxRibbonType or tfxFolder
* @param descriptor             A valid effect descriptor
*/
tfxAPI tfx_effect_descriptor_type tfx_GetDescriptorType(tfx_effect_descriptor descriptor);

/*
Get the full path of a folder, effect or emitter in the library, for example "My Folder/Explosion/Sparks". You can pass this to tfx_GetLibraryDescriptor
* @param descriptor             A valid effect descriptor
*/
tfxAPI const char *tfx_GetDescriptorPath(tfx_effect_descriptor descriptor);

/*
Get the name of a folder, effect or emitter
* @param descriptor             A valid effect descriptor
*/
tfxAPI const char *tfx_GetDescriptorName(tfx_effect_descriptor descriptor);

/*
Get the number of children of a descriptor. A folder's children are effects and an effect's children are its particle and ribbon emitters,
so use this with tfx_GetDescriptorChild and tfx_GetDescriptorType to check whether an effect has any ribbons, for example. Emitters have no children.
* @param descriptor             A valid effect descriptor
*/
tfxAPI tfxU32 tfx_GetDescriptorChildCount(tfx_effect_descriptor descriptor);

/*
Get a child of a descriptor by index
* @param descriptor             A valid effect descriptor
* @param index                  Index of the child, must be less than tfx_GetDescriptorChildCount
*/
tfxAPI tfx_effect_descriptor tfx_GetDescriptorChild(tfx_effect_descriptor descriptor, tfxU32 index);

/*
Read a property of an effect, emitter or ribbon by the same name it's saved under in the library file, for example "spawn_amount",
"emission_type", "relative_position" or "warmup_time". Graphs are not properties, use the graph functions for those.
tfx_PrintDescriptorProperties lists every name that applies to a folder, effect or emitter, with its value.
	tfx_property_value_t value;
	if (tfx_GetDescriptorProperty(emitter, "spawn_amount", &value)) {
		printf("%u\n", value.uint_value);
	}
* @param descriptor             A valid effect descriptor
* @param property_name          Name of the property
* @param out_value              Filled with the value and its type. string_value points into the descriptor so it's only valid while the descriptor is
* @returns                      false if the property name is unknown or doesn't apply to this type of descriptor (an emitter only property read from
								an effect for example), in which case out_value is left untouched
*/
tfxAPI bool tfx_GetDescriptorProperty(tfx_effect_descriptor descriptor, const char *property_name, tfx_property_value_t *out_value);

/*
Print every property name tfx_GetDescriptorProperty takes for a folder, effect or emitter, with its current value and the
tfx_property_value_t member that holds it, a line at a time through tfx_SetPrintCallback if one is set.
	tfx_PrintDescriptorProperties(library, "Explosion/Sparks");
* @param library                A valid pointer to a tfx_library_t
* @param path                   Path to the folder, effect or emitter, as tfx_PrintLibrary lists them
*/
tfxAPI void tfx_PrintDescriptorProperties(tfx_library library, const char *path);

/*
Free all the memory used by a library, including every effect template made from it: drop your template handles, and don't
pass them to tfx_FreeEffectTemplate afterwards. In debug builds a freed handle usually trips the handle assert. Any of its effects still playing on a
stage are removed from that stage first, straight away and without waiting for an update, so their ids go invalid (see
tfx_IsEffectExpired) and there's no need to expire them yourself. Effects from other libraries keep playing.
The stage's instance buffer still holds the removed effects' last instances until its next tfx_UpdateStage, so update
the stages you'll draw and upload their instances before you draw them again. On a frame where the timer didn't tick,
or while paused, tfx_UpdateStage(stage, 0) is enough. Until then tfx_GetNextInstanceKeyRange and tfx_GetNextInstanceBuffer
already leave the removed effects out, so they no longer match the buffer, and the instance getters print a tip.
Don't call this while a stage is updating on another thread: it completes each stage's
in-flight work, but a new update mustn't start until it returns.
* @param library                A valid pointer to a tfx_library_t
*/
tfxAPI void tfx_FreeLibrary(tfx_library library);

/*
Build the GPU shape data, one tfx_gpu_image_data_t per shape frame with its uvs from uv_lookup. Upload it with tfx_GetGPUShapesArray and
tfx_GetGPUShapesSizeInBytes.
* @param library                  A pointer to a tfx_library_t object
* @param shapes                   A pointer to a tfx_gpu_shapes_t object which will fill a buffer with all the shapes
* @param uv_lookup                A function pointer to a function that you need to set up in order to get the uv coordinates from whatever renderer you're using
*/
tfxAPI void tfx_BuildLibraryGPUShapeData(tfx_library library, tfx_gpu_shapes shapes, tfx_uv_lookup uv_lookup);

/*
Get the particle shapes in a library, the list the GPU shape data is built from (see tfx_UpdateLibraryGPUImageData)
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
A number that changes whenever the library's color ramp bitmaps change. Keep the value from when you last uploaded the ramps
and compare it each frame: if it differs, upload them again from tfx_GetColorRampBitmapCount/tfx_GetColorRampBitmap. Loading,
refreshing and customising a template's color graphs with tfx_ScaleTemplateEmitterGraph change it. Creating templates doesn't.
* @param library        A handle to a tfx_library
*/
tfxAPI tfxU32 tfx_GetColorRampsVersion(tfx_library library);

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
libraries in tfxStore, so it takes no library argument. The ribbon shaders and the GPU particle kernel use this
lookup data to update attributes over time. Use with tfx_GetGPUGraphLookupsBufferSizeInBytes to upload it
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
A number that changes whenever the shared gpu graph lookup data changes. Keep the value from when you last uploaded it and
compare it each frame: if it differs, upload tfx_GetGPUGraphLookupsBuffer again. Loading or refreshing any library changes
it, because that rebuilds the lookups of every library, and so does tfx_ScaleTemplateEmitterGraph on a graph the GPU reads.
Creating templates doesn't.
*/
tfxAPI tfxU32 tfx_GetGPUGraphLookupsVersion(void);

/*
Drawing effects from more than one library
------------------------------------------
By default every index a stage writes for the shaders is relative to the effect's own library: the GPU property index in
tfx_instance_t indexes, the color ramp layer in a property and in a ribbon's texture_indexes, and the shape frames. So each
draw binds one library's properties, shapes and color ramps, and a range you draw can only hold effects from one library.

To draw effects from any number of libraries with one set of resources:
  1. Set global_library_indexes in the stage's tfx_stage_info_t. Stages without it are unaffected, and both kinds can run
     side by side.
  2. Call tfx_UpdateGlobalGPUData once per frame before updating your stages. It returns straight away when no library
     changed.
  3. Compare the global versions with the ones you last uploaded and upload what changed:
       - tfx_GetGlobalParticlePropertiesVersion: upload tfx_GetGlobalParticlePropertiesBuffer, sized by
         tfx_GetGlobalParticlePropertiesBufferSizeInBytes since it can grow and shrink.
       - tfx_GetGlobalGPUShapesVersion: upload tfx_GetGlobalGPUShapes with tfx_GetGPUShapesArray and
         tfx_GetGPUShapesSizeInBytes.
       - tfx_GetGlobalColorRampsVersion: check tfx_GetGlobalColorRampLayerVersion for each layer below
         tfx_GetGlobalColorRampLayerCount and upload only the layers that changed, from tfx_GetGlobalColorRampLayer. If
         the count grows past your texture array, recreate it and upload every layer.
     The graph lookups are already global, so keep uploading tfx_GetGPUGraphLookupsBuffer as before.
  4. Bind those buffers for every draw of the stage, whichever libraries its effects came from.
  5. Every library's shape texture_array_index values, whatever your uv_lookup wrote, must address the one texture
     binding your shaders sample. Use one atlas or one bindless table for every library.

Rules:
- Don't sync while another thread is loading, freeing or refreshing a library. A library is registered as soon as it's
  created, so the sync would read it half built. Sync on the thread that loads, or after the load call has returned.
- A library is only placed once a sync has run after it was loaded. Adding an effect from an unplaced library to an
  opted-in stage asserts in debug, and in release the effect isn't added (see tfx_EffectIDIsValid).
- Color ramp layers never move while their library is loaded, because live ribbons keep the layer they spawned with.
  Properties and shapes can move on any sync that changes something, so always upload on a version change.
- If you draw between updates, upload on the frames you upload the stage's instances. The instances you keep drawing
  until the next update were written with the layout you uploaded last.
- To unload a library while others keep playing: tfx_FreeLibrary (it removes that library's effects from every stage),
  then tfx_UpdateGlobalGPUData, then update each stage and upload its instances before drawing it again, since the free
  can move the other libraries' properties and shapes. A paused game or a frame where the timer didn't tick still needs
  that update: tfx_UpdateStage(stage, 0) rewrites the instances without advancing. Or free libraries inside your fixed
  step, before tfx_UpdateGlobalGPUData and the update. The GPU copies of its images (textures, atlas layers) have to
  outlive the instances drawn with them, so free or overwrite those after that update. Its CPU side pixels can go any
  time after tfx_FreeLibrary.
- Limits: 65536 GPU property slots and 256 color ramp layers across every library. If the loaded libraries need more,
  the sync asserts in debug and returns false, and the libraries that didn't fit get no placement. tfx_DescribeStage
  reports it.
- The per-library getters (tfx_GetParticlePropertiesBuffer, tfx_GetLibraryGPUShapes, tfx_GetColorRampBitmap) are
  unchanged and keep each library's own numbering.
*/

/*
Lay out every loaded library in the global GPU buffers so one stage can draw effects from any of them. Call it once per
frame before updating your stages, and after loading, freeing or refreshing a library. It returns straight away when no
library changed. When something did change it waits for any opted-in stage update still in flight before moving anything.
If it moved a library that effects were already drawing from, update the opted-in stages before you draw them with the
new data. Never call it while another thread is loading, freeing or refreshing a library. See "Drawing effects
from more than one library".
* @returns bool                     false if a library didn't fit in the global buffers and has no placement
*/
tfxAPI bool tfx_UpdateGlobalGPUData(void);

/*
Get the particle properties of every placed library, for a stage with global_library_indexes. Each library's properties
are one contiguous run, with color ramp layers and shape frames rewritten to global ones. Valid until the next
tfx_UpdateGlobalGPUData.
*/
tfxAPI void *tfx_GetGlobalParticlePropertiesBuffer(void);

/*
Get the size in bytes of tfx_GetGlobalParticlePropertiesBuffer. It can grow and shrink with each tfx_UpdateGlobalGPUData.
*/
tfxAPI size_t tfx_GetGlobalParticlePropertiesBufferSizeInBytes(void);

/*
A number that changes whenever tfx_UpdateGlobalGPUData changes the global particle properties. Keep the value from when you
last uploaded them and upload them again when it differs.
*/
tfxAPI tfxU32 tfx_GetGlobalParticlePropertiesVersion(void);

/*
Get the number of global color ramp layers, including free ones. Size your color ramp texture array to at least this.
*/
tfxAPI tfxU32 tfx_GetGlobalColorRampLayerCount(void);

/*
Get the bitmap in a global color ramp layer, or NULL if the layer is free. Nothing draws from a free layer, so you can leave
whatever your texture array holds there.
* @param index                  The global layer, below tfx_GetGlobalColorRampLayerCount
*/
tfxAPI tfx_bitmap_t *tfx_GetGlobalColorRampLayer(tfxU32 index);

/*
A number that changes whenever tfx_UpdateGlobalGPUData changes any global color ramp layer. When it differs from the value you
last saw, check tfx_GetGlobalColorRampLayerVersion for each layer.
*/
tfxAPI tfxU32 tfx_GetGlobalColorRampsVersion(void);

/*
A number that changes whenever a global color ramp layer changes, is taken by a library or is freed. Keep one per layer from
when you last uploaded it, and upload the layer again when it differs.
* @param index                  The global layer, below tfx_GetGlobalColorRampLayerCount
*/
tfxAPI tfxU32 tfx_GetGlobalColorRampLayerVersion(tfxU32 index);

/*
Get the shape frames of every placed library as one list, for a stage with global_library_indexes. Use it with
tfx_GetGPUShapesArray, tfx_GetGPUShapesCount and tfx_GetGPUShapesSizeInBytes. Valid until the next tfx_UpdateGlobalGPUData.
*/
tfxAPI tfx_gpu_shapes tfx_GetGlobalGPUShapes(void);

/*
A number that changes whenever tfx_UpdateGlobalGPUData changes the global shape list. Keep the value from when you last
uploaded it and upload it again when it differs.
*/
tfxAPI tfxU32 tfx_GetGlobalGPUShapesVersion(void);

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

  - Call tfx_CompleteStageWork(stage) explicitly. This blocks until all of the stage's outstanding update work has
    finished, after which every buffer is safe to read.
  - Or rely on the getters that synchronise implicitly. The following complete the stage's outstanding work for
    you before returning, so you can call them directly after tfx_UpdateStage without calling
    tfx_CompleteStageWork first:
        tfx_GetInstanceBuffer, tfx_GetInstanceBufferByLayer, tfx_GetInstanceCount, tfx_GetInstanceCountByLayer,
        tfx_GetNextInstanceBuffer, tfx_GetNextInstanceKeyRange, tfx_GetRibbonBuffers, tfx_HasRibbonsToDraw, tfx_NextRibbonDispatch,
        tfx_GetRibbonBufferRequirements, tfx_CopyRibbonDataToStagingBuffers, tfx_GetGPUParticleFrame,
        tfx_GetGPUParticleStats, tfx_GetEffectGPUParticleStats, tfx_SetStageGPUSimulation, tfx_ClearStage and tfx_FreeStage.

Anything that reads stage state through a path NOT in that list (for example reading the instance buffer pointer
you cached last frame, or touching the raw buffers directly) must be preceded by tfx_CompleteStageWork, otherwise
you will race the worker threads. Mutators such as tfx_SetEffectPosition are also only safe once the update work
for the frame has completed so it's a good idea to set positions before the call to tfx_UpdateStage.

Adding an effect (tfx_AddEffectTemplateToStage) is the exception: it waits for any update in flight to finish
before it adds, so it's always safe to call. The cheapest place to add effects is still before tfx_UpdateStage
though, otherwise the add blocks until the update is done. Effects can't be added from inside the stage's own
update (from an effect update callback for example) or while the stage is recording sprite data, those adds
are refused and return an id that fails tfx_EffectIDIsValid.
*/

/*
Create a tfx_stage_info_t object which contains configuration data that you can pass to tfx_CreateStage to set up a stage. You can tweak the config after calling this
function if needed to fine tune the settings. The tfx_stage_info_t field comments give each default.
* @param setup                    A tfx_stage_setup enum which you can use to set the info based on some commonly used templates
*/
tfxAPI tfx_stage_info_t tfx_CreateStageInfo(tfx_stage_setup setup);

/*
Create a stage (effect manager) configured with a tfx_stage_info_t object. Build the info with tfx_CreateStageInfo
and tweak it as needed before passing it in.
* @param info                   A tfx_stage_info_t containing the configuration for the stage.
* @returns tfx_stage            A handle to the new stage, or NULL if the allocation failed (out of memory).
*/
tfxAPI tfx_stage tfx_CreateStage(tfx_stage_info_t info);

/*
Clear the stage as tfx_ClearStage(stage, true) does, and also empty its instance buffers and free each emitter's
storage. Nothing else about the stage changes.
* @param stage                    A pointer to an initialised tfx_stage_t.
*/
tfxAPI void tfx_ReconfigureStage(tfx_stage stage);

/*
Block until all of the stage's in-flight update work (started by tfx_UpdateStage) has finished. Call this before
reading the stage's sprite/instance or ribbon buffers, or mutating effects, unless you are only using the getters
that synchronise implicitly (see the threading contract above). Safe to call redundantly; it is a no-op if there
is no work outstanding.
* @param stage                    A handle to an initialised tfx_stage_t.
*/
tfxAPI void tfx_CompleteStageWork(tfx_stage stage);

/*
Turn on and off whether the stage should sort the effects by depth order. Use tfx_SetStageCamera to set the position of the camera that the stage will
use to update the depth of each effect in the scene.
* @param stage                    A pointer to an initialised tfx_stage_t.
* @param yesno                    A boolean, set to true or false if you want auto ordering on or off respectively
*/
tfxAPI void tfx_ToggleStageOrderEffects(tfx_stage stage, bool yesno);

/*
Give an effect a sort key. Within each layer of the instance buffer effects are ordered by sort key (lowest first) and then by depth if tfx_ToggleStageOrderEffects
is on, so every effect sharing a key is one contiguous range that you can draw with its own resources, for example one range per library (or see
"Drawing effects from more than one library" to let libraries share a range). Use
tfx_GetNextInstanceKeyRange to get those ranges. Run on GPU emitters are ordered by key too, in a second range per key after the CPU instances, but
not by depth. Ribbons aren't in those ranges: they're kept in separate buckets per key, each ribbon dispatch reports its key in sort_key,
and the buckets aren't in key order. Where ribbons draw against billboards is up to you, so to layer effects by key draw each key's ribbon
buckets with that key's billboard ranges (see "Drawing by sort key" in implementations/README.md).
Set the key straight after adding the effect. Live ribbons are stored in their bucket, so the key of an effect can't change while any of its ribbons are
alive: the call asserts and leaves the key as it was. Keys take effect on the next tfx_UpdateStage and effects have a key of 0 until you set one.
This waits for any stage update in progress, so call it before tfx_UpdateStage rather than straight after.
Keep keys to a small set you reuse: every distinct key gets its own ribbon buckets with their own copies of the ribbon paths, and these stay allocated
(counting against max_ribbon_segments) until the stage is cleared. For example when swapping between an old and a reloaded library, alternate between
two keys rather than giving each new library load a new key.
* @param stage                    A pointer to an initialised tfx_stage_t.
* @param effect_id                The id of the effect, returned when you added it to the stage
* @param sort_key                 The key to sort the effect by
*/
tfxAPI void tfx_SetEffectSortKey(tfx_stage stage, tfxEffectID effect_id, tfxU32 sort_key);

/*
Get the sort key of an effect, see tfx_SetEffectSortKey.
* @param stage                    A pointer to an initialised tfx_stage_t.
* @param effect_id                The id of the effect, returned when you added it to the stage
* @returns tfxU32                 The effect's sort key
*/
tfxAPI tfxU32 tfx_GetEffectSortKey(tfx_stage stage, tfxEffectID effect_id);

/*
Give every Run on GPU emitter of one effect the same dispatch key, overriding the keys from tfx_SetTemplateEmitterGPUKey. It takes effect on the
next tfx_GetGPUParticleFrame and can be changed at any time.
* @param stage                    A pointer to an initialised tfx_stage_t.
* @param effect_id                The id of the effect, returned when you added it to the stage
* @param key                      The key every dispatch of the effect carries. tfx_ClearEffectGPUKey goes back to each emitter's template key
*/
tfxAPI void tfx_SetEffectGPUKey(tfx_stage stage, tfxEffectID effect_id, tfxU32 key);

/*
Get the dispatch key set with tfx_SetEffectGPUKey. It fails tfx_GPUKeyIsValid when the effect uses its emitters' template keys.
* @param stage                    A pointer to an initialised tfx_stage_t.
* @param effect_id                The id of the effect, returned when you added it to the stage
* @returns tfxU32                 The effect's GPU key
*/
tfxAPI tfxU32 tfx_GetEffectGPUKey(tfx_stage stage, tfxEffectID effect_id);

//Undo tfx_SetEffectGPUKey, so the effect's dispatches carry each emitter's template key again
tfxAPI void tfx_ClearEffectGPUKey(tfx_stage stage, tfxEffectID effect_id);

//Only checks for the value tfx_GetEffectGPUKey returns when the effect has no key of its own
tfxAPI bool tfx_GPUKeyIsValid(tfxU32 key);

/*
Get the billboard buffer in the stage containing all the sprite instances that were created in the most recent frame. You can use this to copy to a staging buffer to upload to the gpu.
* @param stage                    A pointer to an initialised tfx_stage_t.
*/
tfxAPI tfx_instance_t *tfx_GetInstanceBuffer(tfx_stage stage);

/*
Get the billboard buffer in the stage containing all the sprite instances for a specific layer that were created in the most recent frame. You can use this to copy to a staging buffer to upload to the gpu.
You can then use tfx_GetInstanceCountByLayer for the draw call. The instance buffer is grouped by layer, so each layer is one contiguous range across every effect in the stage and can be drawn with its own pipeline, which tfx_GetLayerFlags picks.
Within a layer, effects are in the same order as the stage's effect list (sorted by depth when auto_order_effects is on). Loop the layers up to tfx_GetLayerCount.
* @param stage                    A pointer to an initialised tfx_stage_t.
* @param layer                    The layer, below tfxLAYERS
*/
tfxAPI tfx_instance_t *tfx_GetInstanceBufferByLayer(tfx_stage stage, tfxU32 layer);

/*
Get the number of instances within the instance buffer of a stage
* @param stage                    A pointer to an initialised tfx_stage_t.
*/
tfxAPI tfxU32 tfx_GetInstanceCount(tfx_stage stage);

/*
Get the instances from the write before the one in tfx_GetInstanceBuffer. This is the buffer that each instance's
captured_index points into. When every update writes, which tfx_UpdateStage and tfx_UpdateStageSubSteps always do, it's the buffer you uploaded
last update, so you only need it if you don't keep that yourself, or you write instances on more than one tfx_UpdateStageSubStep call per frame.
It's only kept when the stage was created with double_buffer_sprites set (the default).
* @param stage                    A pointer to an initialised tfx_stage_t.
*/
tfxAPI tfx_instance_t *tfx_GetPreviousInstanceBuffer(tfx_stage stage);

/*
Get the number of instances in the buffer returned by tfx_GetPreviousInstanceBuffer.
* @param stage                    A pointer to an initialised tfx_stage_t.
*/
tfxAPI tfxU32 tfx_GetPreviousInstanceCount(tfx_stage stage);

/*
Get the number of instances within the instance buffer of a stage for a specific layer.
* @param stage                    A pointer to an initialised tfx_stage_t.
*/
tfxAPI tfxU32 tfx_GetInstanceCountByLayer(tfx_stage stage, tfxU32 layer);

/*
Get the length of the last update in seconds, or of the last step for tfx_UpdateStageSubSteps, after max_frame_length (240 ms) cut it.
* @param stage                    A handle to an initialised tfx_stage_t.
*/
tfxAPI double tfx_GetUpdateTime(tfx_stage stage);

/*
Get the ribbon buffer for a given segment size, tessellation, effect sort key and layer. This will give you all the necessary info and buffer pointers for uploading the ribbon data to the GPU for processing and converting into
a vertex buffer for rendering.
* @param stage                    A pointer to an initialised tfx_stage_t.
* @param segment_count            An unsigned int specifying the ribbon length that you want the rendering info for.
* @param tessellation             Subdivisions across the ribbon's width, 1 to 7. Only half pipe ribbon emitters build with more than 1.
* @param sort_key                 The sort key of the effects whose ribbons you want, see tfx_SetEffectSortKey. Pass 0 if you don't use sort keys.
* @param layer                    The layer of the ribbon emitters you want, below tfx_GetLayerCount.
* @returns						  A pointer to a tfx_ribbon_bucket_t, or null if the stage has no bucket for that segment count, tessellation, key and layer
*/
tfxAPI tfx_ribbon_bucket_t *tfx_GetRibbonBuffers(tfx_stage stage, tfxU32 segment_count, tfxU32 tessellation, tfxU32 sort_key, tfxU32 layer);

/*
Call this to determine whether or not any stage has ribbon_emitters to draw this frame.
* @returns						  True or false
*/
tfxAPI bool tfx_HasRibbonsToDraw(tfx_stage stage);

/*
Get a struct containing the info you need to compute and render ribbon_emitters of a specific length, tessellation, effect sort key and layer.
* @param stage                    A pointer to an initialised tfx_stage_t.
* @param segment_count            An unsigned int specifying the ribbon length that you want the rendering info for.
* @param tessellation             Subdivisions across the ribbon's width, 1 to 7. Only half pipe ribbon emitters build with more than 1.
* @param sort_key                 The sort key of the effects whose ribbons you want, see tfx_SetEffectSortKey. Pass 0 if you don't use sort keys.
* @param layer                    The layer of the ribbon emitters you want, below tfx_GetLayerCount.
* @returns						  A tfx_ribbon_buffer_info_t struct, zeroed if the stage has no bucket for that segment count, tessellation, key and layer
*/
tfxAPI tfx_ribbon_buffer_info_t tfx_GetRibbonBufferInfo(tfx_stage stage, tfxU32 segment_count, tfxU32 tessellation, tfxU32 sort_key, tfxU32 layer);

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
     tfx_CreateRibbonDispatch(), then loop while tfx_NextRibbonDispatch(stage, &dispatch) returns true. Each call
     fills in the per-bucket offsets/counts and points dispatch.globals (or tfx_GetRibbonDispatchGlobals) at the
     push-constant block; dispatch total_segments tells you how many threads to launch.
  4. Iterate the buckets again the same way to record the draw calls, using dispatch.index_count and
     dispatch.index_offset for the indexed draw.

Each bucket holds ribbons of one layer, which dispatch.layer reports. When tfx_GetLayerFlags says a ribbon's layer writes depth, draw
it with tfx_ribbon_cutout_fragment, blending off and depth write on, after the depth layers' billboards and before every blended layer.
To do that, walk the buckets once for the depth writing layers and again for the rest, with a fresh dispatch each time, and skip the
buckets that aren't in the group, so the offsets still add up over every bucket. Ribbon emitters on different layers that share a path
get a copy of the path in each layer's bucket, which counts against max_ribbon_segments.

The dispatch iterator advances an internal cursor on the stage, which goes back to the start when
tfx_NextRibbonDispatch returns false, so you can walk the buckets as many times a frame as you like (compute pass,
then render pass). Only a walk you stop early needs tfx_ResetRibbonDispatchIterator(stage) before the next one.
tfx_NextRibbonDispatch / tfx_GetRibbonBufferRequirements / tfx_CopyRibbonDataToStagingBuffers all implicitly
complete any in-flight stage update work first (see the threading contract above), so they are safe to call
straight after tfx_UpdateStage. For a complete worked example see tfx_impl_dx11__build_ribbon_geometry and
tfx_impl_dx11__render_ribbons in implementations/tfx_impl_dx11.cpp.
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
work item count for the compute dispatch) and globals, the push-constant block for both the compute and render passes.
Implicitly completes any in-flight stage update.
* @param stage                     A handle to an initialised tfx_stage_t.
* @param ribbon_dispatch           A pointer to a tfx_ribbon_dispatch_t (from tfx_CreateRibbonDispatch) updated in place.
* @returns bool                    true while there is another bucket to process, false when iteration is finished.
*/
tfxAPI bool tfx_NextRibbonDispatch(tfx_stage stage, tfx_ribbon_dispatch_t *ribbon_dispatch);

/*
Reset the stage's internal ribbon bucket iterator back to the start. A walk that runs until tfx_NextRibbonDispatch
returns false resets it already, so you only need this after stopping a walk early.
* @param stage                     A handle to an initialised tfx_stage_t.
*/
tfxAPI void tfx_ResetRibbonDispatchIterator(tfx_stage stage);

//The bucket's push-constant block for the compute and render passes, the same pointer as ribbon_dispatch->globals
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
--------------------------------
Ribbon batches (several stages)
--------------------------------
A ribbon batch packs the ribbons of any number of stages into one shared set of segment / ribbon / emitter
buffers, so a renderer uploads once and runs one compute pass and one draw pass for all of them. Which stages go
into the batch is decided every frame, independently of which stages were updated: a stage that wasn't updated
this frame (paused, say) still has to be added if you want it drawn.

Per-frame flow, after tfx_UpdateStage:

  1. tfx_BeginRibbonBatch with your mapped staging buffers and their capacities. Size them for the worst case
     with the tfx_GetTotal*MaxSizeInBytes helpers.
  2. tfx_AddStageToRibbonBatch for each stage you want drawn. Keep the tfx_ribbon_batch_offsets_t it returns.
  3. tfx_GetRibbonBatchSizes gives the exact bytes to upload and the vertex / index buffer sizes to allocate.
  4. For the compute pass and again for the draw pass, for each added stage: create the dispatch with
     tfx_CreateRibbonBatchDispatch(&stage_offsets), then loop tfx_NextRibbonDispatch(stage, &dispatch) exactly as
     for a single stage. The bucket globals already carry the stage's segment and emitter offsets.

A stage's bucket globals hold its batch offsets, so a stage can only be in one batch at a time, and only once.
*/

/*
Start a ribbon batch that writes into the given mapped staging buffers.
* @param segments_dst              Destination for ribbon segment data (tfx_ribbon_segment_t).
* @param segments_capacity         Size of segments_dst in bytes.
* @param ribbons_dst               Destination for ribbon instance data (tfx_ribbon_t).
* @param ribbons_capacity          Size of ribbons_dst in bytes.
* @param emitters_dst              Destination for ribbon emitter data (tfx_gpu_ribbon_emitter_t).
* @param emitters_capacity         Size of emitters_dst in bytes.
* @returns tfx_ribbon_batch_t      An empty batch.
*/
tfxAPI tfx_ribbon_batch_t tfx_BeginRibbonBatch(void *segments_dst, size_t segments_capacity, void *ribbons_dst, size_t ribbons_capacity, void *emitters_dst, size_t emitters_capacity);

/*
Copy a stage's ribbon data into the batch after the stages already added. Implicitly completes any in-flight
stage update.
* @param batch                     A batch from tfx_BeginRibbonBatch.
* @param stage                     A handle to an initialised tfx_stage_t.
* @param stage_offsets             Receives where this stage's data starts in the batch. Pass it to tfx_CreateRibbonBatchDispatch.
* @returns bool                    false if the stage doesn't fit in the batch's capacities, in which case nothing was copied.
*/
tfxAPI bool tfx_AddStageToRibbonBatch(tfx_ribbon_batch_t *batch, tfx_stage stage, tfx_ribbon_batch_offsets_t *stage_offsets);

/*
Byte sizes of everything added to the batch: how much of each staging buffer to upload, and how big the vertex
and index buffers written by the ribbon compute shader need to be.
* @param batch                     A batch from tfx_BeginRibbonBatch.
* @param vertex_size               The size of your vertex struct in bytes, or 0 for the default tfx_ribbon_vertex_t.
*/
tfxAPI tfx_ribbon_batch_sizes_t tfx_GetRibbonBatchSizes(const tfx_ribbon_batch_t *batch, tfxU32 vertex_size);

/*
Create a dispatch for iterating one batched stage's buckets with tfx_NextRibbonDispatch. The vertex, index and
ribbon offsets it produces start where the stage was placed in the batch.
* @param stage_offsets             The offsets returned by tfx_AddStageToRibbonBatch for the stage.
*/
tfxAPI tfx_ribbon_dispatch_t tfx_CreateRibbonBatchDispatch(const tfx_ribbon_batch_offsets_t *stage_offsets);

/*
Worst-case size in bytes of the ribbon segment buffer for a single stage, based on its configured max_ribbon_segments.
Use to size a buffer once at startup instead of resizing per frame.
* @param stage                     A handle to an initialised tfx_stage_t.
*/
tfxAPI size_t tfx_GetSegmentBufferMaxSizeInBytes(tfx_stage stage);

/*
Worst-case size in bytes of the ribbon vertex buffer for a single stage. Pass the size of your own vertex struct,
or 0 to use the library's tfx_ribbon_vertex_t size.
* @param stage                     A handle to an initialised tfx_stage_t.
* @param vertex_size               The size of a single vertex in bytes, or 0 for the default tfx_ribbon_vertex_t.
*/
tfxAPI size_t tfx_GetSegmentVertexBufferMaxSizeInBytes(tfx_stage stage, tfxU32 vertex_size);

/*
Worst-case size in bytes of the ribbon index buffer for a single stage, based on its configured max_ribbon_segments
and tessellation.
* @param stage                     A handle to an initialised tfx_stage_t.
*/
tfxAPI size_t tfx_GetSegmentIndexBufferMaxSizeInBytes(tfx_stage stage);

/*
Worst-case size in bytes of the ribbon instance buffer for a single stage, based on its configured max_ribbons.
* @param stage                     A handle to an initialised tfx_stage_t.
*/
tfxAPI size_t tfx_GetRibbonBufferMaxSizeInBytes(tfx_stage stage);

/*
Worst-case size in bytes of the ribbon emitter buffer for a single stage, based on its configured max_ribbon_emitters.
* @param stage                     A handle to an initialised tfx_stage_t.
*/
tfxAPI size_t tfx_GetEmitterBufferMaxSizeInBytes(tfx_stage stage);

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
A number that changes whenever the library's particle properties buffer changes. Keep the value from when you last uploaded
the buffer and compare it each frame: if it differs, upload tfx_GetParticlePropertiesBuffer again, sized by
tfx_GetParticlePropertiesBufferSizeInBytes since it can grow. Loading, refreshing, rebuilding the shape data and customising a
template's color graphs with tfx_ScaleTemplateEmitterGraph change it. Creating templates doesn't.
* @param library                  A handle to a tfx_library.
*/
tfxAPI tfxU32 tfx_GetParticlePropertiesVersion(tfx_library library);

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
--------------------------------
GPU particle simulation
--------------------------------
Emitters with Run on GPU set are simulated by a compute shader you dispatch, which writes their instances into the
same instance buffer as the CPU particles. Point, User, Area, Ellipse, Disc, Cylinder and Line emitters have their
spawns built in the kernel too, grid spawning and stepped angles included, from a few rows per tick and location, so a
burst of any size costs the CPU next to nothing. Icosphere, Path, Other Emitter and Spawn on Ribbon emitters spawn on the
CPU and upload a row per particle. Set tfx_stage_info_t gpu_particle_capacity to turn it on
for a stage. Without it those emitters simply run on the CPU, as they also do while a stage records sprite data, while
tfx_SetStageGPUSimulation has it off, and when they are the source another emitter spawns from. Every emitter force
runs in the kernel: the library resolves each emitter's forces every tick into force rows that follow the tick records,
so they need nothing more from you.

The reference kernel is in Shaders/timelinefx.slang, a single file with a contents list at its top, entry point
simulate_particles, with tfx_gpu_particle_push_t as its push constants. Every buffer access it makes goes through its
Buffer_access section, which handles each binding model. By default it indexes bindless arrays (Vulkan, D3D12). On Metal every buffer index is a slot in an argument buffer
of gpuAddress values (tfx_metal_buffers). WGSL, and every target compiled with -D TFX_FIXED_BINDINGS (D3D11, SDL_gpu,
OpenGL), binds the pool, instances, graph lookups and upload as four fixed raw buffers. See Shaders/readme.md for each.
To add custom forces, collision or anything else per particle without forking it, write your own
entry point in a Slang file of your own. Shaders/readme.md has a worked example under Adding your own behaviour:
  1. #include "timelinefx.slang", not import: an import fails on DXIL and DXBC, and hides the library's macros from
     your file. Declare a struct that implements ITfxParticleUpdate. Every tick, acceleration is
     added to the particle's acceleration before drag and integration, like the built-in forces, and after_move runs
     once the particle has moved, before its instance is written. hot.position is in the ring's space: build the
     transform with context.ring_transform(cold) and use to_world/to_ring, and direction_to_world/direction_to_ring
     for velocities, forces and normals. Set tfx_particle_remove in hot.flags to kill a particle. It keeps its slot,
     drawn with no size, until its batch expires, and neither method runs for it again.
     To place the particles of User emitters yourself, also declare a struct that implements ITfxParticleSpawn.
     after_spawn gets each particle the kernel builds as a Point at its location, to change as you like, usually its
     position. context.user_data is the location's tfx_SetSpawnLocationUserData. Positions are world space, or local to
     the location when context.relative. Spawns the CPU makes (warmup, no pool, the fallbacks below) go through
     tfx_SetStageSpawnCallback instead, so register one too for the same look everywhere.
     To change how particles are drawn, declare a struct that implements ITfxParticleWrite. before_write gets each live
     particle's tfx_particle_draw on the last tick of the fetch, after its world transform and before it's packed:
     world position, stretch, rotation, alignment, size, intensity, color ramp position, alpha and image frame.
     Per particle state for your hooks: bits 16 to 31 of hot.flags (tfxGPU_PARTICLE_USER_FLAGS_MASK) are yours, zero at
     spawn and kept every tick, and cold.user_value holds the spawn location's tfx_SetSpawnLocationUserData. after_spawn
     can set both on the spawns the kernel generates (spawn.user_flags, spawn.user_value).
  2. Declare your own push constants with tfx_gpu_particle_push_t's twin, tfx_gpu_particle_push, as the first member.
     It takes 64 bytes of the 128 Vulkan guarantees, so your own members start at offset 64. WGSL and TFX_FIXED_BINDINGS
     builds have no push constants, so also declare a uniform twin with TFX_UNIFORM_PUSH(your_push, name) and pick with
     __target_switch on TFX_UNIFORM_PUSH_CASES, as the library's simulate_particles does (Shaders/readme.md). Give
     the push constant block TFX_BINDLESS_REGISTER(b0, space0), so it takes simulate_particles' register on D3D12 and
     its buffer slot on Metal.
  3. In a [numthreads(tfx_group_size, 1, 1)] entry point, fill your update from your push constants and call
     tfx_simulate(thread_id.x, your_push.tfx, update), tfx_simulate(thread_id.x, your_push.tfx, update, spawn_hook)
     with a spawn hook, or tfx_simulate(thread_id.x, your_push.tfx, update, spawn_hook, write_hook) with a write hook.
     tfx_no_particle_update, tfx_no_particle_spawn and tfx_no_particle_write stand in for the hooks you don't have.
  4. Dispatch it exactly like simulate_particles. Send your whole push block once, then only the first 4 bytes,
     dispatch_index, per dispatch.
Your update only runs on the GPU: Run on GPU emitters that fall back to the CPU, for any of the reasons above,
simulate without it.

Three quantities, all in particles: the ring capacity is what one emitter instance can hold, and every spawn location of
that effect shares it. Pool slots are a ring capacity rounded up the way the pool hands it out. gpu_particle_capacity is
the pool, in slots.
Every emitter instance gets a ring of fixed size in the pool when its effect is added: tfx_GetEmitterGPUCapacity, which
is the authored estimate times any location count from tfx_SetTemplateGPUCapacityForLocations, or whatever
tfx_SetTemplateEmitterGPUCapacity set. Size the pool in this order:
  1. Create the effect templates and customise their amounts (tfx_SetTemplateSingleSpawnAmount and the rest). Every
     location spawns the effect's full amount, so to play an effect authored for one location at many, scale its amount
     graphs down with tfx_ScaleTemplateEmitterGraph. If you update slower than 60 times a second, or catch up in one
     tfx_UpdateStage call, set tfx_SetGPUSizingUpdateRate.
  2. tfx_SetTemplateGPUCapacityForLocations on templates that play at several spawn locations at once.
  3. tfx_SetTemplateEmitterGPUCapacity for any one emitter that needs a different ring.
  4. Add tfx_GetTemplateGPUPoolSlots for every template instance that plays at once into gpu_particle_capacity, then
     tfx_CreateStage.
A full ring drops spawns and a full pool leaves the emitter without a ring, both counted by tfx_GetGPUParticleStats and reported once per emitter with a tip (see
tfx_SetPrintCallback). A particle that dies keeps its slot, and is drawn as a zero size billboard, until the rest of the
particles spawned in the same tick have died. So an emitter with a wide life variation costs more vertex work than its
live particle count suggests, and its ring needs room for those slots too.

The current setup (subject to change, this is all WIP)

Setup, once:
  1. Create the particle pool from tfx_GetGPUParticleBufferSizes. It holds the particles between frames and never
     needs resizing or clearing. Its hot, cold and noise sections can share one buffer (pool_bytes) or have one each.
     Where storage buffers are untyped (Vulkan, and Metal's table slots) one index with cold_base and noise_base does
     for all three. D3D12 structured views carry their row size, so make one view per section, either over the
     section with a base of 0 or over the whole buffer with the base. WGSL and TFX_FIXED_BINDINGS builds (D3D11,
     SDL_gpu, OpenGL) always use one pool buffer with the bases.
  2. Upload the graph lookups from tfx_GetGPUGraphLookupsBuffer, as ribbons do.

Per frame, after tfx_UpdateStage:
  1. tfx_GetGPUParticleFrame. It returns false when no update ran since the last call: draw last frame's instance
     buffer again, as you would for the CPU instances.
  2. Upload the arrays. The simplest way is one buffer of upload_bytes, filled with tfx_CopyGPUParticleUpload straight
     into mapped memory, with that one buffer for the dispatch, tick, location and spawn indexes and the frame's
     tick_base, location_base and spawn_base in the push constants. The same rules as the pool apply on D3D12 and
     D3D11. To upload them yourself, every array goes verbatim: ticks is tick_count + force_row_count rows of 64 bytes,
     the tick records then the resolved forces their force_offset points at, and the kernel reads both from one
     buffer. Every offset in the records is relative to its own array, so an array in a buffer of its own has a base
     of 0.
  3. Size the instance buffer for gpu_instance_start + gpu_instance_count instances and upload the CPU instances
     from tfx_GetInstanceBuffer as normal. The GPU instances go after them.
  4. One compute dispatch per tfx_gpu_particle_dispatch_t, (slot_count + 255) / 256 groups of 256. Each loops over
     its own ticks, so there are no barriers between dispatches. Order the dispatches after the previous frame's.
  5. Draw each layer as two ranges: the CPU range as normal and then [gpu_layer_start, + gpu_layer_count). The
     captured_index of a GPU instance points into the previous instance buffer, just like a CPU one. With sort keys,
     draw the CPU then the GPU range of each key from tfx_GetNextInstanceKeyRange instead.

The frame stays valid until the next tfx_UpdateStage. Updates that run without a fetch in between are all returned
by the next fetch, so the GPU never misses a tick. The flip side is that every fetch must be dispatched: the next one's
ticks carry on from the particle state this one leaves in the pool. If you won't dispatch this frame,
don't fetch. But don't keep updating without fetching either: unfetched ticks and their spawn rows build up in memory,
and the next dispatch runs every one of them. A stage you stop drawing should stop updating too.

GPU particles missing or looking wrong? Print tfx_DescribeStage first: it names blocked emitters, emitters dropping
spawns, effects with no spawn locations and updates you haven't fetched. Then turn the GPU simulation off with
tfx_SetStageGPUSimulation(stage, false), clear the stage and add the effects again. If they're still missing on the CPU,
the effect isn't spawning, which the description usually explains. If they appear, the problem is in the host's GPU
integration: the upload, the dispatches or the instance ranges.
*/

/*
Byte sizes and offsets of the particle pool a stage's GPU particles live in. Fixed for the life of the stage.
* @param stage                     A handle to an initialised tfx_stage_t.
* @returns tfx_gpu_particle_buffer_sizes_t   Zero sizes when the stage has no gpu_particle_capacity.
*/
tfxAPI tfx_gpu_particle_buffer_sizes_t tfx_GetGPUParticleBufferSizes(tfx_stage stage);

/*
Everything to upload and dispatch for the GPU particles simulated since the last call. Implicitly completes any
in-flight stage update.
* @param stage                     A handle to an initialised tfx_stage_t.
* @param frame                     Filled with the frame. Its pointers stay valid until the next tfx_UpdateStage.
* @returns bool                    false when there's nothing to dispatch: no update since the last call, or no GPU emitters.
*/
tfxAPI bool tfx_GetGPUParticleFrame(tfx_stage stage, tfx_gpu_particle_frame_t *frame);

/*
Copies every array of a frame into one upload buffer at the frame's upload offsets, typically straight into mapped
staging memory. The gaps between sections are left as they are.
* @param frame                     A frame from tfx_GetGPUParticleFrame.
* @param destination               At least frame->upload_bytes of writable memory.
*/
tfxAPI void tfx_CopyGPUParticleUpload(const tfx_gpu_particle_frame_t *frame, void *destination);

/*
Pool use, ring slots in use, an estimate of the live particles and the counters that say a capacity needs raising. The
slots in use are always more than the live particles: a particle that dies keeps its slot until every particle spawned
in the same tick has died, because only then can the ring move past them. Implicitly completes any in-flight stage update.
* @param stage                     A handle to an initialised tfx_stage_t.
*/
tfxAPI tfx_gpu_particle_stats_t tfx_GetGPUParticleStats(tfx_stage stage);

/*
tfx_GetGPUParticleStats for one effect's rings, to see which effect is using the pool or dropping spawns. capacity is
still the whole pool. dropped_spawns only counts the rings the effect has now, so an emitter that has finished takes its
drops with it. Implicitly completes any in-flight stage update.
* @param stage                     A handle to an initialised tfx_stage_t.
* @param effect_id                 The effect from tfx_AddEffectTemplateToStage. An expired one returns zeros apart from capacity.
*/
tfxAPI tfx_gpu_particle_stats_t tfx_GetEffectGPUParticleStats(tfx_stage stage, tfxEffectID effect_id);

/*
Turn a stage's GPU simulation off or back on, for comparing it with the CPU simulation or falling back to the CPU.
It's the first thing to try when GPU particles look wrong, see the end of the GPU particle simulation notes above.
Effects added while it's off run their Run on GPU emitters on the CPU. Effects already running keep the simulation
they started with, so clear the stage (tfx_ClearStage) and add your effects again after switching. Implicitly completes
any in-flight stage update.
* @param stage                     A handle to an initialised tfx_stage_t.
* @param enabled                   false runs Run on GPU emitters of newly added effects on the CPU.
*/
tfxAPI void tfx_SetStageGPUSimulation(tfx_stage stage, bool enabled);

/*
Whether effects added now would simulate their Run on GPU emitters on the GPU: the stage has a gpu_particle_capacity and
tfx_SetStageGPUSimulation hasn't turned it off.
* @param stage                     A handle to an initialised tfx_stage_t.
*/
tfxAPI bool tfx_StageGPUSimulationIsOn(tfx_stage stage);

/*
Sets the slowest update rate GPU ring capacities are sized for, 60 updates a second by default. A ring holds each
particle's slot for its life plus up to one update, and an update's spawns all arrive at once, so slower updates need
bigger rings. A catch-up done in one call, tfx_UpdateStage(stage, tick * 3), is one long update: it spawns all three
ticks' particles at once, up to max_frame_length (240 ms) of them. tfx_UpdateStageSubSteps steps are an update each, so
catching up that way needs no extra room.
Call it before tfx_SetTemplateGPUCapacityForLocations and tfx_GetTemplateGPUPoolSlots: each records the rate on the
template and its rings are sized at that rate from then on, so changing it afterwards can't make the rings disagree with
the pool slots you added up. Templates you never sized, and tfx_GetEmitterGPUCapacity on a library's own emitters, use
the current rate. Rings already allocated keep their size, and it has no effect on capacities set with
tfx_SetTemplateEmitterGPUCapacity.
* @param updates_per_second        The slowest rate you'll update stages at, counting a single call catch-up of N ticks
                                  as your rate divided by N.
*/
tfxAPI void tfx_SetGPUSizingUpdateRate(float updates_per_second);

/*
The ring capacity an instance of a Run on GPU emitter gets, in particles. Without an override this is what one spawn
location needs, times the location count tfx_SetTemplateGPUCapacityForLocations set on its template. Every location
adds its own amount, so an effect at N locations needs about N times the one location figure. It's worked out from the
emitter's amounts and the sizing rate each call (the rate its template was sized at, or tfx_SetGPUSizingUpdateRate), so it follows tfx_SetTemplateSingleSpawnAmount, and tfx_ScaleTemplateEmitterGraph on the amount graphs, until tfx_SetTemplateEmitterGPUCapacity
fixes it.
* @param emitter                   A handle to an emitter descriptor.
* @returns tfxU32                  0 when the emitter doesn't have Run on GPU set.
*/
tfxAPI tfxU32 tfx_GetEmitterGPUCapacity(tfx_effect_descriptor emitter);

/*
Sets the ring capacities of every Run on GPU emitter in a template for a number of spawn locations playing at once, and
returns the pool slots they take. It stores the location count rather than fixed capacities, so amounts you change
later with tfx_SetTemplateSingleSpawnAmount are still followed. It replaces any capacity set before, including
tfx_SetTemplateEmitterGPUCapacity: call that after this to size one emitter differently, then take the slots from
tfx_GetTemplateGPUPoolSlots. Use tfx_GetEffectLifetime divided by how often you add a location, plus 1, to estimate how
many overlap. Add the return value of every template you'll play into gpu_particle_capacity. Only affects effects added
after the call.
* @param t                         A handle to an effect template
* @param location_count            Spawn locations the effect will play at, at the same time, at most
* @returns tfxU32                  The pool slots the template's rings take, rounded the way the pool hands them out
*/
tfxAPI tfxU32 tfx_SetTemplateGPUCapacityForLocations(tfx_effect_template t, tfxU32 location_count);

/*
The pool slots the template's Run on GPU emitters take with their current capacities, rounded the way the pool hands
them out. Add it into gpu_particle_capacity for every instance of the template that will play at once. It records the
current tfx_SetGPUSizingUpdateRate on the template, so the rings it adds match this count even if the rate changes later.
* @param t                         A handle to an effect template
*/
tfxAPI tfxU32 tfx_GetTemplateGPUPoolSlots(tfx_effect_template t);

/*
When a stage updates particles it creates work queues to handle the work. By default these each have room for max_emitters entries which should be
more than enough for most situations. However you can increase the sizes here if needed. You only need to set this manually if you hit one of the asserts when these
run out of space or you anticipate a huge amount of emitters and particles to be used (> million). On the other hand, you might be tight on memory in which case you
could reduce the numbers as well if needed (they don't take a lot of space though)
* @param stage                     A pointer to an initialised tfx_stage_t.
* @param spawn_work_max            The maximum amount of spawn work entries
* @param control_work_max          The maximum amount of control work entries
* @param age_work_max              The maximum amount of age_work work entries
*/
tfxAPI void tfx_SetStageWorkQueueSizes(tfx_stage stage, tfxU32 spawn_work_max, tfxU32 control_work_max, tfxU32 age_work_max);

/*
Every emitter on a stage stores its particles in a particle array that grows as it needs to. When the emitter ends its array isn't freed, it's kept
idle so the next emitter of the same type can reuse it without allocating or growing again. The following functions free those idle arrays. Arrays
that running emitters are using are never touched, so these are safe to call at any time apart from inside the stage's update, for example from an
effect update callback. They wait for an update that is running on another thread to finish first.
*/

/*
Free every idle particle array on the stage, and the stage's lists of them.
* @param stage                     A pointer to an initialised tfx_stage_t.
*/
tfxAPI void tfx_FreeAllParticleArrays(tfx_stage stage);

/*
Free the idle particle arrays of one emitter type.
* @param stage                     A pointer to an initialised tfx_stage_t.
* @param path                      The emitter's full path in its library, for example "Explosion/Sparks". An effect's path does nothing here, use
                                   tfx_FreeEffectParticleArrays to free every emitter in an effect.
*/
tfxAPI void tfx_FreeEmitterParticleArrays(tfx_stage stage, const char *path);

/*
Free the idle particle arrays of every emitter below an effect or folder. An emitter's descriptor frees just that emitter.
* @param stage                     A pointer to an initialised tfx_stage_t.
* @param descriptor                A valid tfx_effect_descriptor of an effect, folder or emitter
*/
tfxAPI void tfx_FreeEffectParticleArrays(tfx_stage stage, tfx_effect_descriptor descriptor);

/*
Free the particle arrays that have been idle for at least idle_milliseconds of stage updates. Idle arrays are kept oldest first, so it stops at the
first one of each emitter type that isn't due. That makes it cheap enough to call regularly, every second or so, or when the game changes level.
* @param stage                     A pointer to an initialised tfx_stage_t.
* @param idle_milliseconds         How long an array must have been idle, in milliseconds of updates, to be freed. 0 frees every idle array.
*/
tfxAPI void tfx_TrimParticleArrays(tfx_stage stage, float idle_milliseconds);

/*
Get the number of CPU particles a stage held after its last update, so 0 straight after adding an effect. Particles of Run on GPU emitters
aren't counted, see tfx_GetGPUParticleStats. Includes particles that have died but still hold their place until the oldest particle in their
emitter dies (drawn at zero size), so an emitter whose lives vary reads high by about half its life variation over its average life
(about 6% for lives of 2000 to 2250 ms).
* @param stage                     A pointer to a tfx_stage_t
* @returns tfxU32                  The number of CPU particles the last update held
*/
tfxAPI tfxU32 tfx_GetParticleCount(tfx_stage stage);

/*
Get the current ribbon count for a stage
* @param stage                     A pointer to a tfx_stage_t
* @returns tfxU32                  The total number of ribbons currently being updated
*/
tfxAPI tfxU32 tfx_GetRibbonCount(tfx_stage stage);

/*
Get the current number of effects that are currently being updated by a stage
* @param stage                     A pointer to a tfx_stage_t
* @returns tfxU32                  The total number of effects currently being updated
*/
tfxAPI tfxU32 tfx_GetEffectCount(tfx_stage stage);

/*
Get the number of emitters the CPU updated in the stage's last update, so 0 straight after adding an effect, until the next update.
Run on GPU emitters aren't counted.
* @param stage                     A pointer to a tfx_stage_t
* @returns tfxU32                  The number of emitters the CPU updated in the last update
*/
tfxAPI tfxU32 tfx_GetEmitterCount(tfx_stage stage);

/*
Everything that says whether a stage is healthy: counts, GPU pool use, blocked GPU emitters, dropped spawns, effects
waiting for spawn locations, GPU updates that haven't been fetched and the error flags of the libraries in use.
Implicitly completes any in-flight stage update.
* @param stage                     A handle to an initialised tfx_stage_t
*/
tfxAPI tfx_stage_diagnostics_t tfx_GetStageDiagnostics(tfx_stage stage);

/*
tfx_GetStageDiagnostics as text: a line of counts, a line for the GPU pool, then a line for each problem naming the
emitter or effect it's about, with emitters named by the path the tfx_SetTemplate functions take. Print it or compare it
between runs to check an integration without looking at the screen. Implicitly completes any in-flight stage update.
* @param stage                     A handle to an initialised tfx_stage_t
* @param buffer                    Filled with as many whole lines as fit, always null terminated. Can be NULL with a size of 0
* @param buffer_size               The size of buffer in bytes
* @returns tfxU32                  The buffer size the whole text needs, terminator included
*/
tfxAPI tfxU32 tfx_DescribeStage(tfx_stage stage, char *buffer, tfxU32 buffer_size);

/*
Get the current number of free effects in the stage that are ready to be re-used
* @param stage                     A pointer to a tfx_stage_t
* @returns tfxU32                  The total number of free effects in the stage
*/
tfxAPI tfxU32 tfx_GetFreeEffectCount(tfx_stage stage);

/*
Set the seed for the stage for random number generation. Setting the seed can determine how an emitters spawns particles, so if you set the seed before adding an effect to the stage
then the effect will look the same each time. Any seed, 0 included, is fine.
* @param stage                     A pointer to an initialised tfx_stage_t.
* @param seed                      Any value
*/
tfxAPI void tfx_SetStageSeed(tfx_stage stage, tfxU64 seed);

/*
Add an effect to a tfx_stage_t from an effect template. This is the only way to add an effect, so make a template with tfx_CreateEffectTemplate
even if you don't customise it. If the stage is part way through an update this waits for it to finish first,
see the threading contract above tfx_CreateStageInfo. An effect that spawns at user spawn locations (see
tfx_EffectTemplateSpawnsAtUserLocations) spawns nothing until you give it one with tfx_AddSpawnLocation.
* @param stage                      A pointer to an initialised tfx_stage_t.
* @param effect_template			The tfx_effect_template_t object that you want to add to the stage. It must have already been prepared by calling tfx_CreateEffectTemplate
* @returns tfxEffectID				The id to manipulate the effect with, for example with tfx_SetEffectPosition. Check it with tfx_EffectIDIsValid: it fails if the stage is
									out of room for it
*/
tfxAPI tfxEffectID tfx_AddEffectTemplateToStage(tfx_stage stage, tfx_effect_template effect);

//Only checks for the id a failed add returns, use tfx_IsEffectExpired to find out if the effect is still running
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
Returns true if the effect template is an effect that will spawn particles that eventually expire.
For an effect that spawns at user spawn locations, true means each location plays out. The effect itself keeps running
with no locations until you expire it, unless its locations were added with tfxSpawnLocationAdd_expire_effect.
* @param tfx_effect_template		A handle to the effect template
* @returns bool						True if effect will expire on it's own, false if it spawns particles forever.
*/
tfxAPI bool tfx_IsFiniteEffect(tfx_effect_template effect);

/*
Whether the effect spawns at the user spawn locations you add with tfx_AddSpawnLocation. Such an effect spawns nothing
until it has one.
* @param tfx_effect_template		A handle to the effect template
*/
tfxAPI bool tfx_EffectTemplateSpawnsAtUserLocations(tfx_effect_template effect);

/*
Whether an effect descriptor spawns at user spawn locations, for walking a library without making templates.
* @param effect						A handle to an effect descriptor
*/
tfxAPI bool tfx_EffectSpawnsAtUserLocations(tfx_effect_descriptor effect);

/*
Returns the maximum possible lifetime of an effect taking into account particle lifetime and lifetime variation.
If the effect is not finite then it returns 0. You should not use this effect in realtime as it can be expensive
depending on the effect, so use it ahead of time after creating the effect template.
For an effect that spawns at user spawn locations the effect itself keeps running, so this is how long one location
plays, from when it's added until its last particle dies. It's also when a location is removed automatically, unless
it was added with tfxSpawnLocationAdd_no_auto_remove. Divide it by how often you add locations to estimate how many
overlap, for tfx_SetTemplateGPUCapacityForLocations. To have the effect end after its locations, add them with
tfxSpawnLocationAdd_expire_effect.
A finite effect expires with its last particle, at most this long after it starts unless you raise its life with
tfx_SetEffectLifeMultiplier, and leaves the stage an update after that (see tfx_IsEffectExpired). If you add one effect per
burst, max_effects needs bursts per second * this / 1000, plus a couple.
* @param tfx_effect_template		A handle to the effect template
* @returns float					The time in milliseconds	
*/
tfxAPI float tfx_GetEffectLifetime(tfx_effect_template effect);

/*
Get the time of a specific bookmark that was set in the effect via the editor by index. You can use bookmarks to
see if a certain part of an effect has been reached yet. Bookmarks are in order of time, so index 0 will
always be the earliest bookmark
* @param tfx_effect_template		A handle to the effect template
* @param bookmark_index				The bookmark slot, 0 to tfxMAX_BOOKMARKS - 1
* @returns float					The time in milliseconds, or 0 if there is no bookmark in that slot
*/
tfxAPI float tfx_GetBookmarkTime(tfx_effect_template effect, tfxU32 bookmark_index);

/*
Get the number of bookmarks set in the effect. Bookmarks are packed, so slots 0 to count - 1 are all in use.
* @param tfx_effect_template		A handle to the effect template
* @returns tfxU32					The number of bookmarks, 0 to tfxMAX_BOOKMARKS
*/
tfxAPI tfxU32 tfx_GetBookmarkCount(tfx_effect_template effect);

/*
Get the name the bookmark was given in the editor. The pointer stays valid until the template is freed or the
library is refreshed.
* @param tfx_effect_template		A handle to the effect template
* @param bookmark_index				The bookmark slot, 0 to tfxMAX_BOOKMARKS - 1
* @returns const char*				The bookmark name, or an empty string if there is no bookmark in that slot
*/
tfxAPI const char *tfx_GetBookmarkName(tfx_effect_template effect, tfxU32 bookmark_index);

/*
Find the slot of a bookmark by its name. Bookmarks are ordered by time, so moving one in the editor can change
the slot it lives in. Look up the slots you need by name once after creating the template, and again after
calling tfx_RefreshLibrary, then pass them to tfx_IsBookmarkCrossed rather than hard coding slot numbers.
If more than one bookmark has the same name then the earliest one is returned.
* @param tfx_effect_template		A handle to the effect template
* @param name						The name of the bookmark, case sensitive
* @returns tfxU32					The bookmark slot. It fails tfx_BookmarkIndexIsValid if no bookmark has that name
*/
tfxAPI tfxU32 tfx_GetBookmarkIndex(tfx_effect_template effect, const char *name);

//Only checks for the value tfx_GetBookmarkIndex returns when no bookmark has the name
tfxAPI bool tfx_BookmarkIndexIsValid(tfxU32 bookmark_index);

/*
Check if a specific bookmark time has been crossed for an effect running in the stage. Use tfx_GetBookmarkIndex
to get the slot from the bookmark's name.
* @param stage							A pointer to an initialised tfx_stage_t
* @param effect_id					The effect id that is in the stage
* @param bookmark_index				The bookmark slot, 0 to tfxMAX_BOOKMARKS - 1
* @returns bool						True if the effect has passed that bookmark. Always false for an empty slot.
*/
tfxAPI bool tfx_IsBookmarkCrossed(tfx_stage stage, tfxEffectID effect_id, tfxU32 bookmark_index);

/*
Set the delta time used when warming up effects. Either match the fixed time step you're using like 1000/60 for 60fps 
or used mulitples of that value for higher performance at the cost of precision
* @param stage							A pointer to an initialised tfx_stage_t. 
* @param double						The delta time measured in milliseconds
*/
tfxAPI void tfx_SetWarmUpDeltaTime(tfx_stage stage, double delta_time);

/*
If an effect is already in the stage and you want to advance it be a given amount of time
without rendering those frame then you can use this function.
* @param stage							A pointer to an initialised tfx_stage_t. 
* @param tfxEffectID				The effect id that is in the stage
* @param float						The time to advance by measured in milliseconds
*/
tfxAPI void tfx_AdvanceEffectTime(tfx_stage stage, tfxEffectID effect_id, float time);

/*
Update a stage which manages effects, ideally from a fixed step loop. elapsed is in milliseconds: at 60 updates a second it's 16.666667.
If you interpolate particles in the vertex shader, the previous instance buffer has to be the one you drew last, so when the fixed loop has to
run more than once to catch up, use tfx_UpdateStageSubSteps, which runs every tick at its own length and writes the instances once, or call this
once with the elapsed time multiplied by the number of ticks rather than once per tick. Either way the previous buffer is then all of those ticks
behind, so your timer's alpha isn't the factor to interpolate with, tfx_GetStageInterpolation turns it into the one that is.
A single call is also one long update to the rings of Run on GPU emitters, which tfx_SetGPUSizingUpdateRate sizes for, sub-steps aren't.
An elapsed longer than max_frame_length (240 ms) is cut to it, for spawning, aging and motion alike, so a long stall doesn't burst.
In pseudo code, where the Timer functions are your own:

	TimerAccumulate(game->timer);
	int pending_ticks = TimerPendingTicks(game->timer);	//The number of times the update loop will run this frame.

	while (TimerDoUpdate(game->timer)) {
		if (pending_ticks > 0) {
			tfx_UpdateStage(game->stage, FrameLength * pending_ticks);
			//Set the pending ticks to 0 so we don't run the update again this frame.
			pending_ticks = 0;
		}

		TimerUnAccumulate(game->timer);
	}
	TimerSet(game->timer);	//Set the timer's alpha: how far its accumulator is into the next tick, 0 to 1
	//Every frame, not just update frames. Pass it to a uniform or push constant for the vertex shader
	float lerp = tfx_GetStageInterpolation(game->stage, TimerAlpha(game->timer), FrameLength);

	//Only upload the instance buffer to the gpu if the stage was updated.
	if (TimerUpdateWasRun(game->timer)) {
		RenderParticles(game->stage, game);
	}

An elapsed of 0 advances nothing and spawns nothing but still rewrites the instances. Use it to bring the instances up to
date on a frame where the timer didn't tick, for example straight after tfx_FreeLibrary.
* @param stage                 A pointer to an initialised tfx_stage_t.
* @param elapsed               The time to advance the stage by, in milliseconds. 0 or more
*/
tfxAPI void tfx_UpdateStage(tfx_stage stage, double elapsed);

/*
Update a stage in a number of fixed size steps with one call. All of the steps run on the stage's update thread, so like tfx_UpdateStage this returns
without waiting for them, whereas calling tfx_UpdateStage once per step makes each call wait for the step before it. Use this when you simulate at a
fixed rate and run more than one step in a frame. Only the last step writes instances, the steps before it only simulate, which makes them cheaper.
captured_index points into what the last update wrote, the buffer you uploaded then, so the vertex shader interpolates across all of this update's
steps. Positions you set on effects since the last update are spread evenly across the steps. tfx_UpdateStage(stage, elapsed) is the same as
tfx_UpdateStageSubSteps(stage, elapsed, 1).
* @param stage                 A pointer to an initialised tfx_stage_t.
* @param step_length           The length of each step in milliseconds
* @param step_count            The number of steps to run, must be at least 1
*/
tfxAPI void tfx_UpdateStageSubSteps(tfx_stage stage, double step_length, tfxU32 step_count);

/*
Update a stage by a single step, choosing whether the step writes instances. Use this instead of tfx_UpdateStageSubSteps when you need to change
things in the stage between steps, like moving an effect along with something in your own fixed step simulation. Each call waits for the step
before it to finish, so do your own work for the next step before calling it again and the two will overlap. Steps that don't write instances
are cheaper. Pass tfxUpdateStepFlags_write_instances on the last step before you draw, which is what tfx_UpdateStageSubSteps does: captured_index
then points into the buffer you uploaded last frame, so you interpolate across all the steps since then. To interpolate across only the last step,
write on the last two steps instead and upload tfx_GetPreviousInstanceBuffer alongside the current buffer. Don't do that with Run on GPU emitters:
their captured indexes point into the layout of the last tfx_GetGPUParticleFrame, which the previous buffer's CPU instances can overlap. Don't read
the instance buffer after a step that didn't write instances, it doesn't describe the particles until the next step that does.
* @param stage                 A pointer to an initialised tfx_stage_t.
* @param step_length           The length of the step in milliseconds
* @param flags                 tfxUpdateStepFlags_write_instances to write the instance buffer, or tfxUpdateStepFlags_none
*/
tfxAPI void tfx_UpdateStageSubStep(tfx_stage stage, double step_length, tfxUpdateStepFlags flags);

/*
Get the factor to interpolate between the previous and the current instance buffer with, from your fixed step timer's alpha. The previous buffer
is as far behind as the time the stage simulated since its last write before this one, which is several ticks after a catch-up update and less
than you asked for when max_frame_length clamped the update, so alpha on its own would make particles jump or stutter. This returns
1 - (1 - alpha) * tick_length / that time, clamped to 0 to 1, which draws the particles one tick behind the simulation however much the last
update covered. Call it every frame, from the thread that calls tfx_UpdateStage. If you update once a frame with the frame's own elapsed time
and don't interpolate, you don't need it.
* @param stage                 A pointer to an initialised tfx_stage_t.
* @param alpha                 How far your timer's accumulator is into the next tick, 0 to 1
* @param tick_length           The length of one tick of your timer in milliseconds, the same unit you pass to tfx_UpdateStage
* @returns                     The interpolation factor, 0 to 1. 1 before the first update and after an update of 0
*/
tfxAPI float tfx_GetStageInterpolation(tfx_stage stage, float alpha, float tick_length);

/*
Get the image pointer for a sprite. Use this when rendering particles in your renderer. The pointer that is returned will be the pointer that you set in your shape loader function
used when loading an effect library. Generally you shouldn't need to use this function, simply copy the whole instance buffer in the stage to your staging buffer to be
copied to the gpu in one go.
* @param stage                 A pointer to an initialised tfx_stage_t. 
* @param property_indexes    The value in the instance_data->property_indexs[i] when iterating over the instance_data in your render function
  @returns                    void* pointer to the image
*/
tfxAPI void *tfx_GetSpriteImagePointer(tfx_stage stage, tfxU32 property_indexes);

/*
Get the total number of instances ready for rendering in the stage.
* @param stage                 A pointer to an initialised tfx_stage_t.
*/
tfxAPI tfxU32 tfx_TotalSpriteCount(tfx_stage stage);

/*
Clear all particles, instance_data and effects in a stage, and reset the stage's GPU dropped spawn count. If you don't need to use the stage again then call tfx_FreeStage to also
free all the memory associated with the stage.
* @param stage                     A pointer to an initialised tfx_stage_t.
* @param free_particle_banks    Set to true if you want to free the memory associated with the particle banks and release back to the memory pool
*/
tfxAPI void tfx_ClearStage(tfx_stage stage, bool free_particle_banks);

/*
Free all the memory used in the stage. You don't have to call this - tfx_EndTimelineFX frees any stages
that are still alive - but you can free them earlier if you want to reclaim the memory. Calling it is safe either
way, the stage is deregistered from the shutdown sweep so it won't be freed twice.
* @param stage                     A pointer to an initialised tfx_stage_t.
*/
tfxAPI void tfx_FreeStage(tfx_stage stage);

//[Effects functions for altering effects that are currently playing out in a stage]

/*
Expire an effect by telling it to stop spawning particles. It expires on the update after its last remaining particle dies, see tfx_IsEffectExpired.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect that you want to expire, returned when you added it to the stage
*/
tfxAPI void tfx_SoftExpireEffect(tfx_stage stage, tfxEffectID effect_id);

/*
Soft expire all the effects in a stage so that the particles complete their animation first
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
*/
tfxAPI void tfx_SoftExpireAll(tfx_stage stage);

/*
Expire an effect by telling it to stop spawning particles and remove all associated particles immediately.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect that you want to expire, returned when you added it to the stage
*/
tfxAPI void tfx_HardExpireEffect(tfx_stage stage, tfxEffectID effect_id);

/*
Find out if an effect has finished. That's at once after tfx_HardExpireEffect, otherwise on the update after which it has
nothing alive and can never spawn again: a soft expired effect, or a finite one that has played out, ends with its last
particle. An effect that could still spawn, like a looping one in a gap between spawns or one that spawns by distance
travelled, has to have nothing alive for 1 second first, so that a gap doesn't end it.
Nothing of an expired effect draws after the next tfx_UpdateStage, and the update after that takes it off the stage (see
tfx_IsEffectOnStage). Its slot is then recycled for the next effect that you add, but the old id is never reused so it will keep
reporting as expired and any other calls you make with it are ignored. Use this to know when to let go of an id you're holding on to.
Like the other effect functions this should be called either before tfx_UpdateStage or after tfx_CompleteStageWork so that it isn't
racing the stage update.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @returns               True if the effect has expired or the id is stale or invalid
*/
tfxAPI bool tfx_IsEffectExpired(tfx_stage stage, tfxEffectID effect_id);

/*
Whether an effect still holds its slot on the stage. It turns false an update after tfx_IsEffectExpired turns true, when the slot
is freed for reuse, so count effects against max_effects with this rather than tfx_IsEffectExpired. Same threading rule as
tfx_IsEffectExpired.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id         The id of the effect, returned when you added it to the stage
* @returns                 False once the effect has left the stage, or if the id is stale or invalid
*/
tfxAPI bool tfx_IsEffectOnStage(tfx_stage stage, tfxEffectID effect_id);

/*
Get effect user data
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect that you want to expire, returned when you added it to the stage
* @returns                void* pointing to the user data set in the effect. See tfx_effect_template_t::SetUserData() and tfx__set_effect_user_data()
*/
tfxAPI void *tfx_GetEffectUserData(tfx_stage stage, tfxEffectID effect_id);

/*
More for use in the editor, this function updates emitter base values for any effects that are currently running after their graph values have been changed.
*/
tfxAPI void tfx_UpdateStageBaseValues(tfx_stage stage);

/*
Set the stage camera. This is used to calculate particle depth if you're using depth ordered particles so it needs to be updated each frame.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param front            An array of 3 floats representing a normalised 3d vector describing the direction that the camera is pointing
* @param position        An array of 3 floats representing the position of the camera in 3d space
*/
tfxAPI void tfx_SetStageCamera(tfx_stage stage, float front[3], float position[3]);

/*
Set the function that places the particles of User emitters, which spawn like a Point and leave the positions to you, whenever the CPU spawns
them: CPU emitters, and GPU emitters during warmup, on a stage without a GPU pool or when the GPU can't generate their spawns. Change the
positions in place. Called from worker threads during tfx_UpdateStage, once per location share of a tick's spawns, or twice when the share wraps
the particle buffer. Only read data that stays put until the update completes, and never call the stage api.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param callback          The function to call, or NULL to stop calling one
* @param user_data         Passed to every call
*/
tfxAPI void tfx_SetStageSpawnCallback(tfx_stage stage, tfx_spawn_callback callback, void *user_data);

/*
Each effect in the stage can have bounding box which you can decide to keep updated or not if you wanted to do any offscreen culling of effects. Theres some
extra overhead to keep the bounding boxes updated but that can be made back if you have a number of effect particles offscreen that don't need to be drawn.
* @param stage					A pointer to a tfx_stage_t where the effect is being managed
* @param yesno				Set to true or false if you want the bounding boxes to be udpated.
*/
tfxAPI void tfx_KeepBoundingBoxesUpdated(tfx_stage stage, bool yesno);

/*
Set the effect user data for an effect already added to a stage
* @param stage					A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id		The id of the effect that you want to expire, returned when you added it to the stage
* @param user_data			A void* pointing to the user_data that you want to store in the effect
*/
tfxAPI void tfx_SetEffectUserData(tfx_stage stage, tfxEffectID effect_id, void *data);

/*
Force a stage to only run in single threaded mode. In other words, only use the main thread to update particles
* @param stage             A pointer to a tfx_stage_t.
* @param switch_on        true or false to use a single thread or not
*/
tfxAPI void tfx_ForceStageSingleThreaded(tfx_stage stage, bool switch_on);

/*
Get the transform vectors for a sprite's previous position so that you can use that to interpolate between that and the current sprite position
* @param stage             A pointer to a tfx_stage_t.
* @param layer            The index of the sprite layer
* @param index            The sprite index of the sprite that you want the captured sprite for.
* @param position         This should be a pointer to a vec3 that you pass in that will get loaded with the position of the instance
*/
tfxAPI void tfx_GetCapturedInstanceTransform(tfx_stage stage, tfxU32 layer, tfxU32 index, float out_position[3]);

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
Make a stage stop spawning. This will mean that all emitters in the stage will no longer spawn any particles so all currently running effects will expire
as the remaining particles come to the end of their life. Any single particles will also get flagged to expire
* @param stage             A pointer to a tfx_stage_t.
* @param yesno            True = disable spawning, false = enable spawning
*/
tfxAPI void tfx_DisableStageSpawning(tfx_stage stage, bool yesno);

/*
Error-handling contract for the stage functions that take a tfxEffectID (below):
Every such function validates the id. An id that fails validation makes the call a no-op: mutators return
without touching any state, functions that return a pointer return NULL, and functions that return a
value return 0.
Effect ids are produced by tfx_AddEffectTemplateToStage and carry a generation for their slot in the stage.
When an effect expires and is removed, or the stage is cleared, its id goes stale and is never handed out
again even though the slot is recycled. So holding on to an id after the effect has gone is safe: calls made
with it are silently skipped rather than landing on whatever effect took the slot. Use tfx_IsEffectExpired
to find out when an id has gone stale. Malformed ids, such as one from a failed add (see tfx_EffectIDIsValid), still trip an assert in debug builds.
The animation manager functions that take a tfxAnimationID are validated the same way but those ids have no
generation, so they assert when invalid and should be treated as invalid once the instance has expired.
*/

/*
Set the position of an effect
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param x                The x value of the position
* @param y                The y value of the position
* @param z                The z value of the position
*/
tfxAPI void tfx_SetEffectPosition(tfx_stage stage, tfxEffectID effect_id, float x, float y, float z);

/*
Set the position of an effect
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param position        A float[3] array containing the x, y and z coordinates
*/
tfxAPI void tfx_SetEffectPositionVec3(tfx_stage stage, tfxEffectID effect_id, float position[3]);

/*
Move an effect to a new position without tweening, so nothing is spawned along the line between the old and new position and
no extra particles are added by spawn distance spacing. Use this rather than tfx_SetEffectPosition when the effect jumps somewhere
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param x                The x value of the position
* @param y                The y value of the position
* @param z                The z value of the position
*/
tfxAPI void tfx_TeleportEffect(tfx_stage stage, tfxEffectID effect_id, float x, float y, float z);

/*
Move an effect to a new position without tweening, so nothing is spawned along the line between the old and new position and
no extra particles are added by spawn distance spacing. Use this rather than tfx_SetEffectPositionVec3 when the effect jumps somewhere
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param position        A float[3] array containing the x, y and z coordinates
*/
tfxAPI void tfx_TeleportEffectVec3(tfx_stage stage, tfxEffectID effect_id, float position[3]);

/*
Add a location for an effect to spawn particles at. The effect must have "Spawn at User Locations" set in the editor. Every emitter in the
effect spawns at every location the effect has, with its shape placed around each one, so an area emitter gives you its whole area at each
location rather than a single point, and a path emitter gives you its whole path at each location. Line and path emitters that traverse their
edge or use their path as a trajectory need relative position set, since that is what anchors them to their location (edge traversal sets it
for you). Other Emitter and Spawn on Ribbon emitters carry on spawning from
their source emitter, which puts them at the locations too. Spawn on Ribbon particles are world space here, like the ribbons they spawn on,
so relative position is ignored for them. Use this to have one effect handle lots of sources, for example the
smoke trails of every ship in a fleet, which is much cheaper than an effect per ship.
Locations persist until you remove them. Changes are applied the next time tfx_UpdateStage is called, so it's safe to call this while the
previous update is still running, but only from the thread that calls tfx_UpdateStage.
Locations are in world space. Particles from emitters set to relative position follow the location they spawned at (position only, the
effect rotation still applies) and are removed along with it, which suits things like a glow on each ship. Relative particles at a transient
location only last one update.
One thing is shared rather than per location: an emitter set to spawn on a grid marches a single grid position for the whole emitter, so the
locations take it in turns through the grid instead of each getting a grid of their own. Path emitters are the exception: each location gets
its own rotated paths, rotation cycles, stagger, maximum paths and grid march, as if the effect were playing on its own there.
Ribbons are always world space here: a ribbon is anchored where its location was when it spawned and stays there for its life, so relative
position is ignored for a ribbon emitter spawning at user locations, and so is ribbon lag, which needs it.
The spawn amount of an emitter is per location and is taken from its amount graph at that location's age, so a burst that fades out does so
at every location in turn. An emitter's delay also counts from when the location was added, and single shot emitters fire once per location.
Every other graph is sampled at the effect's or emitter's own age as it would be for any effect, so it affects the particles of all the
locations equally. Bear in mind that the effect keeps ageing for as long as it's on the stage, so graphs meant to shape a single burst are
best left flat here.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id         The id of the effect, returned when you added it to the stage
* @param position          A float[3] array containing the x, y and z coordinates
* @param flags             A combination of tfxSpawnLocationAddFlags, or tfxSpawnLocationAdd_none.
                           tfxSpawnLocationAdd_transient removes the location automatically after the next update so you don't need to keep
                           the id at all, which is what you want for a one off burst like an impact.
                           If the effect is finite, meaning every emitter stops spawning of its own accord, the location removes itself once
                           the effect has played out there, so you don't have to tidy it up yourself. Pass
                           tfxSpawnLocationAdd_no_auto_remove to keep it instead. You need that if you change the effect's life with
                           tfx_SetEffectLifeMultiplier after adding the location, because the time is worked out from the effect as authored
                           and a longer life would outlive it. An effect that never ends is never auto removed either way.
                           The effect itself keeps running when it has no locations, so one shared effect can spawn burst after burst.
                           tfxSpawnLocationAdd_expire_effect soft expires the effect once it has no locations left instead, for when
                           you add one effect per burst: tfx_IsEffectExpired then turns true after the last location's particles have
                           gone, as it would for any finite effect. A location kept with tfxSpawnLocationAdd_no_auto_remove never goes,
                           so it never lets the effect expire.
* @returns                 A tfxSpawnLocationID to use with tfx_UpdateSpawnLocation and tfx_RemoveSpawnLocation. It includes the effect so
                           you don't need to keep the effect index with it. tfxINVALID_SPAWN_LOCATION if the effect doesn't spawn at user locations,
                           with a tip the first time on each stage
*/
tfxAPI tfxSpawnLocationID tfx_AddSpawnLocation(tfx_stage stage, tfxEffectID effect_id, float position[3], tfxSpawnLocationAddFlags flags);

/*
Move a spawn location. Particles spawned this update are spread along the line between the old and new position.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param location_id       The id returned by tfx_AddSpawnLocation
* @param position          A float[3] array containing the x, y and z coordinates
*/
tfxAPI void tfx_UpdateSpawnLocation(tfx_stage stage, tfxSpawnLocationID location_id, float position[3]);

/*
Move a spawn location without tweening, so nothing is spawned along the line between the old and new position and no extra
particles are added by spawn distance spacing. The last position given before the next update wins, and if any of the calls since the last
update was a teleport, none of the move is tweened.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param location_id       The id returned by tfx_AddSpawnLocation
* @param position          A float[3] array containing the x, y and z coordinates
*/
tfxAPI void tfx_TeleportSpawnLocation(tfx_stage stage, tfxSpawnLocationID location_id, float position[3]);

/*
Orient the emitter at a spawn location: the shape it spawns in and the direction its particles are emitted, so an area emitter's box turns
with it and a point emitter's emission direction turns with it. A location has an identity rotation until you call this. Line emitters that
traverse their edge or use their path as a trajectory take the rotation but not the direction, theirs is fixed along the line.
When the rotation takes effect depends on whether the emitter is set to relative position, matching what relative already means for the
location's position. Relative emitters are attached, so they turn with the location and keep turning as it does, which is what lets you aim
a beam at something. Emitters that aren't relative take the rotation as each particle spawns, so turning a location aims what it spawns from
then on and leaves what it already spawned where it was, which is what you want for a trail.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param location_id       The id returned by tfx_AddSpawnLocation
* @param rotation          A float[4] quaternion as x, y, z, w
*/
tfxAPI void tfx_SetSpawnLocationRotation(tfx_stage stage, tfxSpawnLocationID location_id, float rotation[4]);

/*
Per location multipliers, the location's own version of the tfx_SetEffect*Multiplier functions, so one effect can drive sources that differ
from each other. They multiply on top of the effect's, and they all default to 1. Every one of them is applied as each particle or ribbon
spawns, so changing one never disturbs anything already alive. That includes intensity, which on the effect also changes the particles
already out.
Life: how long particles and ribbons spawned here last. A location that's counting down to removing itself, whether auto removing or soft
expiring with a time worked out for you, has its countdown lengthened to cover a larger multiplier, but a smaller one never shortens it.
Particle width and height: the base size of particles spawned here, and width is also the width of ribbons. If the effect is set to a uniform
particle size then width scales both and height does nothing, as it does for the effect.
Velocity, spin and weight: the base velocity, spin (all three axes) and weight of particles spawned here. Velocity never scales the
emitter's spawn impulse, so it barely changes a burst that gets its speed from the impulse.
Intensity: the intensity of particles spawned here.
Splatter: the random offset given to particles and ribbons spawned here.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param location_id       The id returned by tfx_AddSpawnLocation
*/
tfxAPI void tfx_SetSpawnLocationLifeMultiplier(tfx_stage stage, tfxSpawnLocationID location_id, float life);
tfxAPI void tfx_SetSpawnLocationParticleWidthMultiplier(tfx_stage stage, tfxSpawnLocationID location_id, float width);
tfxAPI void tfx_SetSpawnLocationParticleHeightMultiplier(tfx_stage stage, tfxSpawnLocationID location_id, float height);
tfxAPI void tfx_SetSpawnLocationVelocityMultiplier(tfx_stage stage, tfxSpawnLocationID location_id, float velocity);
tfxAPI void tfx_SetSpawnLocationSpinMultiplier(tfx_stage stage, tfxSpawnLocationID location_id, float spin);
tfxAPI void tfx_SetSpawnLocationIntensityMultiplier(tfx_stage stage, tfxSpawnLocationID location_id, float intensity);
tfxAPI void tfx_SetSpawnLocationSplatterMultiplier(tfx_stage stage, tfxSpawnLocationID location_id, float splatter);
tfxAPI void tfx_SetSpawnLocationWeightMultiplier(tfx_stage stage, tfxSpawnLocationID location_id, float weight);

/*
The size of the shape an emitter spawns in at a location, the location's version of tfx_SetEffectWidthMultiplier, tfx_SetEffectHeightMultiplier
and tfx_SetEffectDepthMultiplier, and multiplied on top of them. Width, height and depth are the emitter's own x, y and z, so a line's length
is its height. All default to 1 and take effect for particles and ribbons spawned from now on. If the emitter sets its amount by its area,
a resized location's amount changes with it so the density stays the same, as it does for the effect.
Trajectory emitters, which follow their line or path every frame, only take the height, as a uniform scale of the whole trajectory.
Emitters that traverse their line or path edge ignore all three.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param location_id       The id returned by tfx_AddSpawnLocation
*/
tfxAPI void tfx_SetSpawnLocationWidthMultiplier(tfx_stage stage, tfxSpawnLocationID location_id, float width);
tfxAPI void tfx_SetSpawnLocationHeightMultiplier(tfx_stage stage, tfxSpawnLocationID location_id, float height);
tfxAPI void tfx_SetSpawnLocationDepthMultiplier(tfx_stage stage, tfxSpawnLocationID location_id, float depth);

/*
Offset the noise field sampled by particles spawned at a location, the location's version of tfx_SetEffectBaseNoiseOffset. Locations
otherwise share the effect's noise offset, so giving each its own stops them all swirling the same way. Added to the effect's offset,
0 by default, and only read by emitters using a noise field.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param location_id       The id returned by tfx_AddSpawnLocation
* @param noise_offset      The offset to add for particles spawned from now on
*/
tfxAPI void tfx_SetSpawnLocationNoiseOffset(tfx_stage stage, tfxSpawnLocationID location_id, float noise_offset);

/*
Give a spawn location a value of your own for the spawn hooks of User emitters, such as an index into your mesh, bone or trail data. 0 by default.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param location_id       The id returned by tfx_AddSpawnLocation
* @param user_data         The value to pass
*/
tfxAPI void tfx_SetSpawnLocationUserData(tfx_stage stage, tfxSpawnLocationID location_id, tfxU32 user_data);

/*
Stop a spawn location spawning without removing it. The location still exists and can still be moved, and relative particles already
spawned there stay attached to it and live out their life as normal, so pausing a location leaves what it has already made alone.
The location also stops ageing while it's paused, so its amount graph carries on from where it left off when you resume it, and a single
shot emitter that hasn't reached its delay yet still fires once it does.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param location_id       The id returned by tfx_AddSpawnLocation
*/
tfxAPI void tfx_PauseSpawnLocation(tfx_stage stage, tfxSpawnLocationID location_id);

/*
Start a paused spawn location spawning again, carrying on from the age it was paused at.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param location_id       The id returned by tfx_AddSpawnLocation
*/
tfxAPI void tfx_ResumeSpawnLocation(tfx_stage stage, tfxSpawnLocationID location_id);

/*
Retire a spawn location gently. It stops spawning straight away like tfx_PauseSpawnLocation, but rather than staying put it removes itself
once the particles already spawned there have had time to live out their life, so relative particles fade out naturally instead of vanishing
the moment the location goes, which is what tfx_RemoveSpawnLocation does to them.
The id stays valid while the location drains and you can carry on moving it with tfx_UpdateSpawnLocation, so something that's been destroyed
can keep drifting while its particles follow it out. Use tfx_SpawnLocationIsValid to find out when it has gone.
Sub effects spawned by the particles at the location aren't counted, so pass a time of your own if you have a long lived one.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param location_id       The id returned by tfx_AddSpawnLocation
* @param milliseconds      How long to wait before removing the location. Pass 0 to have it worked out from the longest lived emitter in the
                           effect, which is the worst case of its base life plus variation so it always outlasts its particles
*/
tfxAPI void tfx_SoftExpireSpawnLocation(tfx_stage stage, tfxSpawnLocationID location_id, float milliseconds);

/*
Check whether a spawn location is currently spawning. False if it has been paused. Note that this doesn't tell you the emitters are
actually putting particles out, only that the location isn't holding them back: they might be waiting out a delay or have finished.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param location_id       The id returned by tfx_AddSpawnLocation
*/
tfxAPI bool tfx_SpawnLocationIsSpawning(tfx_stage stage, tfxSpawnLocationID location_id);

/*
Remove a spawn location. Particles already spawned there live out their life, except relative ones which are removed with it. The id is
invalid straight away.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param location_id       The id returned by tfx_AddSpawnLocation
*/
tfxAPI void tfx_RemoveSpawnLocation(tfx_stage stage, tfxSpawnLocationID location_id);

/*
Remove all the spawn locations from an effect.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id         The id of the effect, returned when you added it to the stage
*/
tfxAPI void tfx_ClearSpawnLocations(tfx_stage stage, tfxEffectID effect_id);

/*
Check whether a spawn location id still refers to a location in the effect.
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param location_id       The id returned by tfx_AddSpawnLocation
*/
tfxAPI bool tfx_SpawnLocationIsValid(tfx_stage stage, tfxSpawnLocationID location_id);

/*
Move an Effect by a specified amount relative to the effect's current position
* @param stage             A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param amount            A float[3] array containing the amount to move in the x, y and z planes
*/
tfxAPI void tfx_MoveEffectVec3(tfx_stage stage, tfxEffectID effect_id, float amount[3]);

/*
Move an Effect by a specified amount relative to the effect's current position
* @param stage            A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id        The id of the effect, returned when you added it to the stage
* @param x                The amount to move in the x plane
* @param y                The amount to move in the y plane
* @param z                The amount to move in the z plane
*/
tfxAPI void tfx_MoveEffect(tfx_stage stage, tfxEffectID effect_id, float x, float y, float z);

/*
Get the current position of an effect
* @param stage           A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @return                float[3] array containing the effect position
*/
tfxAPI void tfx_GetEffectPositionVec3(tfx_stage stage, tfxEffectID effect_id, float out_position[3]);

/*
You can use this function to get the instances of a specific effect on one layer. The instance buffer is grouped by layer, so an effect's instances are only
contiguous within each layer; call this once per layer to get all of them.
* @param stage						A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id			The id of the effect, returned when you added it to the stage
* @param layer					The layer to get the instances for, below tfx_GetLayerCount
* @param tfxU32					Pass in a pointer to a tfxU32 which will be set to the number of instance_data in the buffer, 0 on a layer the effect doesn't use.
* @return						tfx_instance_t pointer to the buffer
*/
tfxAPI tfx_instance_t *tfx_GetEffectInstanceBuffer(tfx_stage stage, tfxEffectID effect_id, tfxU32 layer, tfxU32 *sprite_count);

/*
You can use this function to get each effect's instances on one layer for every effect that is currently active in the stage. Generally you would call this
in a while loop inside a for loop over each layer up to tfx_GetLayerCount, calling tfx_ResetInstanceBufferLoopIndex before each layer. It returns false
straight away on a layer with no instances.
* @param stage						A pointer to a tfx_stage_t where the effect is being managed
* @param layer					The layer to get the instances for
* @param tfx_sprite_billboard_t	Pass in a pointer which will be set to the current sprite buffer containing all of the sprite data for this frame.
* @param tfx_effect_instance_data_t   Pass in a second pointer which will be set to the tfx_effect_instance_data_t containing all of the sprite buffer data. This can be used to gain access to all the sprite data if using double buffered instance_data (to interpolated with the previous frame).
* @param tfxU32					Pass in a pointer to a tfxU32 which will be set to the number of instance_data in the buffer.
* @return						true or false if the next billboard buffer was found. False will be returned once there are no more effect sprite buffers in the stage
*/
tfxAPI bool tfx_GetNextInstanceBuffer(tfx_stage stage, tfxU32 layer, tfx_instance_t **sprites_soa, tfx_effect_instance_data_t **effect_sprites, tfxU32 *sprite_count);

/*After calling GetNextBillboard/SpriteBuffer in a while loop you can call this to reset the index for the next frame
* @param stage						A pointer to a tfx_stage_t
*/
tfxAPI void tfx_ResetInstanceBufferLoopIndex(tfx_stage stage);

/*
Get the instances on one layer one sort key at a time (see tfx_SetEffectSortKey). Call it in a while loop until it returns false, then call
tfx_ResetInstanceBufferLoopIndex before the next layer. Ranges come in key order, and each one holds two contiguous ranges of every effect sharing the key:
its CPU instances, and the instances its Run on GPU emitters write after the CPU ones (see tfx_GetGPUParticleFrame). Bind the key's resources and draw the
CPU range then the GPU range, so within a key and layer an effect's GPU particles draw over its CPU ones. Either range can be empty.
This shares its loop index with tfx_GetNextInstanceBuffer, so don't interleave the two.
* @param stage                    A pointer to an initialised tfx_stage_t.
* @param layer                    The layer to get the ranges for, below tfx_GetLayerCount
* @param range                    Set to the key and its two ranges
* @returns bool                   true while there is another range, false when there are no more on this layer
*/
tfxAPI bool tfx_GetNextInstanceKeyRange(tfx_stage stage, tfxU32 layer, tfx_instance_key_range_t *range);

/*
Set the roll of an effect in radians. Setting any of roll, pitch or yaw replaces the effect's three rotation graphs from the next tfx_UpdateStage
for the rest of its life, and the two angles you didn't set keep the value they had (0 before the first update). Like the other effect setters,
call it before tfx_UpdateStage, or once the update has completed (see the threading contract).
* @param stage           A pointer to a tfx_stage_t where the effect is being managed.
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param roll            The roll in radians
*/
tfxAPI void tfx_SetEffectRoll(tfx_stage stage, tfxEffectID effect_id, float roll);

/*
Set the pitch of an effect in radians. See tfx_SetEffectRoll, the same rules apply.
* @param stage           A pointer to a tfx_stage_t where the effect is being managed.
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param pitch           The pitch in radians
*/
tfxAPI void tfx_SetEffectPitch(tfx_stage stage, tfxEffectID effect_id, float pitch);

/*
Set the yaw of an effect in radians. See tfx_SetEffectRoll, the same rules apply.
* @param stage           A pointer to a tfx_stage_t where the effect is being managed.
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param yaw             The yaw in radians
*/
tfxAPI void tfx_SetEffectYaw(tfx_stage stage, tfxEffectID effect_id, float yaw);

/*
Set the orientation of an effect with a quaternion, overriding the rotation graphs set in the TimelineFX editor. tfx_SetEffectRoll/Pitch/Yaw
and tfx_PointEffectAt replace it rather than adjust it, they work from their own angles.
* @param stage           A pointer to a tfx_stage_t where the effect is being managed.
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param rotation        A float[4] quaternion as x, y, z, w. It should be normalized before you pass it in.
*/
tfxAPI void tfx_SetEffectRotation(tfx_stage stage, tfxEffectID effect_id, float rotation[4]);

/*
Orient an effect towards a specific point. Call it before tfx_UpdateStage, or once the update has completed (see the threading contract).
* @param stage           A pointer to a tfx_stage_t where the effect is being managed.
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param x, y, z         The coordinates that you want to point the effect at
* @param face            The local axis of the effect that should end up aiming at those coordinates. Use tfxEffectFace_up
*                        unless the effect was authored to fire along a different axis.
*                        Aiming only fixes two of the three angles, so the third is left as it is for you to drive: the
*                        up/down/left/right faces preserve the effect's roll, the forwards/backwards faces preserve its yaw.
*                        Note that this overrides the effect's authored pitch, yaw and roll graphs together, so once you
*                        call this you own all three angles for the lifetime of the effect.
*                        Coordinates are in the library's own +Y up space. If your renderer is Z up then convert the target
*                        point into that space first, the same way you convert the effect's position going in.
*/
tfxAPI void tfx_PointEffectAt(tfx_stage stage, tfxEffectID effect_id, float x, float y, float z, tfx_effect_face face);

/*
Set the emitter width multiplier of an effect. It multiplies the width of every emitter in the effect, an easy way to resize area, line, ellipse
etc. emitters. It replaces the effect's own emitter width graph from the next tfx_UpdateStage for the rest of the effect's life, and leaves height
and depth to their graphs. Like the other effect setters, call it before tfx_UpdateStage, or once the update has completed (see the threading contract).
* @param stage           A pointer to a tfx_stage_t where the effect is being managed.
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param width           The width multiplier, in place of the graph's value (1 unless the effect animates it)
*/
tfxAPI void tfx_SetEffectWidthMultiplier(tfx_stage stage, tfxEffectID effect_id, float width);

/*
Set the emitter height multiplier of an effect. The same as tfx_SetEffectWidthMultiplier for heights.
* @param stage           A pointer to a tfx_stage_t where the effect is being managed.
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param height          The height multiplier, in place of the graph's value (1 unless the effect animates it)
*/
tfxAPI void tfx_SetEffectHeightMultiplier(tfx_stage stage, tfxEffectID effect_id, float height);

/*
Set the emitter depth multiplier of an effect. The same as tfx_SetEffectWidthMultiplier for depths.
* @param stage           A pointer to a tfx_stage_t where the effect is being managed.
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param depth           The depth multiplier, in place of the graph's value (1 unless the effect animates it)
*/
tfxAPI void tfx_SetEffectDepthMultiplier(tfx_stage stage, tfxEffectID effect_id, float depth);

/*
Set the life multiplier of an effect. Multiplies on top of the effect's own global graph for it, so the authored animation still plays, and
stays until you change it. 1 is no change.
* @param stage           A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param life            How long particles and ribbons spawned from now on last before expiring.
*/
tfxAPI void tfx_SetEffectLifeMultiplier(tfx_stage stage, tfxEffectID effect_id, float life);

/*
Set the particle width multiplier of an effect. Multiplies on top of the effect's own global graph for it, so the authored animation still plays, and
stays until you change it. 1 is no change.
* @param stage           A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param width           The width of particles and ribbons spawned from now on. If the effect is set to a uniform particle size it scales the height as well.
*/
tfxAPI void tfx_SetEffectParticleWidthMultiplier(tfx_stage stage, tfxEffectID effect_id, float width);

/*
Set the particle height multiplier of an effect. Multiplies on top of the effect's own global graph for it, so the authored animation still plays, and
stays until you change it. 1 is no change.
* @param stage           A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param height          The height of particles spawned from now on. Does nothing if the effect is set to a uniform particle size, the width multiplier covers both.
*/
tfxAPI void tfx_SetEffectParticleHeightMultiplier(tfx_stage stage, tfxEffectID effect_id, float height);

/*
Set the velocity multiplier of an effect. Multiplies on top of the effect's own global graph for it, so the authored animation still plays, and
stays until you change it. 1 is no change.
* @param stage           A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param velocity        The base velocity of particles spawned from now on. Spawn impulse isn't scaled.
*/
tfxAPI void tfx_SetEffectVelocityMultiplier(tfx_stage stage, tfxEffectID effect_id, float velocity);

/*
Set the spin multiplier of an effect. Multiplies on top of the effect's own global graph for it, so the authored animation still plays, and
stays until you change it. 1 is no change.
* @param stage           A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param spin            The base roll, pitch and yaw spin of particles spawned from now on.
*/
tfxAPI void tfx_SetEffectSpinMultiplier(tfx_stage stage, tfxEffectID effect_id, float spin);

/*
Set the intensity multiplier of an effect. Multiplies on top of the effect's own global graph for it, so the authored animation still plays, and
stays until you change it. 1 is no change.
* @param stage           A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param intensity       The intensity of every particle in the effect, including the ones already alive.
*/
tfxAPI void tfx_SetEffectIntensityMultiplier(tfx_stage stage, tfxEffectID effect_id, float intensity);

/*
Rotate the hue of every emitter in an effect. Adds to the effect's own global hue graph, so the authored animation still plays, and stays
until you change it. Emitters set to exclude themselves from global hue are left alone. 0 is no change.
* @param stage           A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param degrees         The hue rotation of every particle and ribbon in the effect, including the ones already alive.
*/
tfxAPI void tfx_SetEffectHue(tfx_stage stage, tfxEffectID effect_id, float degrees);

/*
Set the splatter multiplier of an effect. Multiplies on top of the effect's own global graph for it, so the authored animation still plays, and
stays until you change it. 1 is no change.
* @param stage           A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param splatter        The random offset given to particles and ribbons spawned from now on.
*/
tfxAPI void tfx_SetEffectSplatterMultiplier(tfx_stage stage, tfxEffectID effect_id, float splatter);

/*
Set the weight multiplier of an effect. Multiplies on top of the effect's own global graph for it, so the authored animation still plays, and
stays until you change it. 1 is no change.
* @param stage           A pointer to a tfx_stage_t where the effect is being managed
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param weight          The weight of particles spawned from now on.
*/
tfxAPI void tfx_SetEffectWeightMultiplier(tfx_stage stage, tfxEffectID effect_id, float weight);

/*
Set the overall scale of an effect, a simple way to change its size. It replaces the effect's own overall scale graph from the next
tfx_UpdateStage for the rest of the effect's life. Call it before tfx_UpdateStage, or once the update has completed (see the threading contract).
* @param stage           A pointer to a tfx_stage_t where the effect is being managed.
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param overall_scale   The overall scale, in place of the graph's value (1 unless the effect animates it)
*/
tfxAPI void tfx_SetEffectOverallScale(tfx_stage stage, tfxEffectID effect_id, float overall_scale);

/*
Set the base noise offset for an effect
* @param stage             A pointer to a tfx_stage_t where the effect is being managed.
* @param effect_id       The id of the effect, returned when you added it to the stage
* @param noise_offset    A float of the amount that you want to set the effect noise offset to. By default when an effect is added to a stage a random noise offset will be set based on the Base Noise Offset Range property. Here you can override that
						value by setting it here. The most ideal time to set this would be immediately after you have added the effect to the stage, but you could call it any time you wanted for a constantly changing noise offset.
*/
tfxAPI void tfx_SetEffectBaseNoiseOffset(tfx_stage stage, tfxEffectID effect_id, float noise_offset);


//--------------------------------
//Animation_manager
//--------------------------------

/*
Set the position of an animation instance
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param animation_id            The id returned by tfx_AddAnimationInstance
* @param position                A float[3] array containing the x, y and z coordinates
*/
tfxAPI void tfx_SetAnimationPosition(tfx_animation_manager animation_manager, tfxAnimationID animation_id, float position[3]);

/*
Set the scale of an animation instance
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param animation_id            The id returned by tfx_AddAnimationInstance
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
to playback pre recorded particle effects as opposed to using a stage that simulates the particles in
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
Add an animation instance to the animation manager.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param path                    tfxKey path hash of the effect name and path: effect.path_hash
* @param start_frame            Starting frame of the animation
* @returns                        The index id of the animation instance. You can use this to reference the animation when changing position, scale etc
								It fails tfx_AnimationIDIsValid if there is no room in the animation manager or the effect recorded nothing
*/
tfxAPI tfxAnimationID tfx_AddAnimationInstanceByKey(tfx_animation_manager animation_manager, tfxKey path, tfxU32 start_frame);

/*
Add an animation instance to the animation manager.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param path                    const char * name of the effect. If the effect was in a folder then specify the whole path
* @param start_frame            Starting frame of the animation
* @returns                        The index id of the animation instance. You can use this to reference the animation when changing position, scale etc
								It fails tfx_AnimationIDIsValid if there is no room in the animation manager or the effect recorded nothing
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
Build the GPU shape data for an animation manager's sprite data, one tfx_gpu_image_data_t per shape frame with its uvs from uv_lookup. Upload it
with tfx_GetGPUShapesArray and tfx_GetGPUShapesSizeInBytes.
* @param animation_manager		  A pointer to an tfx_animation_manager_t object
* @param shapes                   A pointer to a tfx_gpu_shapes_t object which will fill a buffer with all the shapes
* @param uv_lookup                A function pointer to a function that you need to set up in order to get the uv coordinates from whatever renderer you're using
*/
tfxAPI void tfx_BuildAnimationManagerGPUShapeData(tfx_animation_manager animation_manager, tfx_gpu_shapes shapes, tfx_uv_lookup uv_lookup);

/*
Get the particle shapes in an animation manager, the list tfx_BuildAnimationManagerGPUShapeData builds from
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
Get the buffer memory address for the graph lookups of the animation manager's ribbon emitters. These are separate from
tfx_GetGPUGraphLookupsBuffer because sprite data can be played back without a library: bind this table wherever the ribbon shaders
expect the graph lookups when drawing baked ribbons. It's saved in the sprite data file and only changes when sprite data is added or
loaded.
* @param animation_manager        A handle to a tfx_animation_manager
* @returns void*                  A pointer to tfx_gpu_graph_data_t rows, tfx_GetAnimationRibbonGraphLookupsSizeInBytes long
*/
tfxAPI void *tfx_GetAnimationRibbonGraphLookupsBufferPointer(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the animation manager's ribbon graph lookups from tfx_GetAnimationRibbonGraphLookupsBufferPointer
* @param animation_manager        A handle to a tfx_animation_manager
* @returns size_t                 The size in bytes, 0 when there are no ribbons
*/
tfxAPI size_t tfx_GetAnimationRibbonGraphLookupsSizeInBytes(tfx_animation_manager animation_manager);

/*
Returns true or false if the animation manager contains effects with ribbons
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @returns bool 	              True if the animation manager has ribbons
*/
tfxAPI bool tfx_AnimationManagerHasRibbons(tfx_animation_manager animation_manager);

/*
Baked ribbons go through the same tfx_ribbon_compute and ribbon shaders as live ones. Every frame, after tfx_UpdateAnimationManager:
  1. Dispatch tfx_ribbon_playback with tfx_ribbon_data_push_t, (tfx_GetTotalRibbonsThatNeedDrawing + 127) / 128 groups. It
     reads the ribbon data, the ribbon offsets and the animation instances and writes a tfx_ribbon_t per ribbon to draw.
  2. Walk the buckets with tfx_CreateAnimationRibbonDispatch and tfx_NextAnimationRibbonDispatch, exactly as you walk a stage's with
     tfx_NextRibbonDispatch: one tfx_ribbon_compute dispatch per bucket with tfx_GetRibbonDispatchGlobals as its push constants, then
     walk them again for one indexed draw each. Set the globals' camera_position yourself, an animation manager has no camera.
tfx_ribbon_compute reads the emitter buffer from tfx_CopyAnimationRibbonEmitters, the segments from
tfx_GetAnimationRibbonSegmentsBufferPointer and the graph lookups from tfx_GetAnimationRibbonGraphLookupsBufferPointer (in place of the
library's tfx_GetGPUGraphLookupsBuffer, tfx_ribbon_vertex too), which only change when sprite data is added or loaded.
*/

/*
Get the buffer sizes baked ribbon playback needs this frame. Call after tfx_UpdateAnimationManager.
* @param animation_manager        A handle to a tfx_animation_manager
* @param vertex_size              The size of the vertex tfx_ribbon_compute writes, sizeof(tfx_ribbon_vertex_t) for the reference shaders
* @returns tfx_ribbon_batch_sizes_t   segment_buffer_size_in_bytes is the shared segments, ribbon_buffer_size_in_bytes the tfx_ribbon_t
                                  rows tfx_ribbon_playback writes, emitter_buffer_size_in_bytes the emitter rows, and the vertex and
                                  index sizes what tfx_ribbon_compute writes for this frame's buckets
*/
tfxAPI tfx_ribbon_batch_sizes_t tfx_GetAnimationRibbonBufferSizes(tfx_animation_manager animation_manager, tfxU32 vertex_size);

/*
Get the size of the emitter buffer that tfx_ribbon_compute reads for baked ribbons, one row per ribbon emitter in the animation manager
* @param animation_manager        A handle to a tfx_animation_manager
* @returns size_t                 The size in bytes
*/
tfxAPI size_t tfx_GetAnimationRibbonEmittersSizeInBytes(tfx_animation_manager animation_manager);

/*
Write the emitter rows that tfx_ribbon_compute reads for baked ribbons. Baked ribbon positions are already in world space, so every emitter has an
identity transform and only carries what the shader samples: its graph lookups, angle type, path samples, morph target and noise settings.
The graph lookup offsets index tfx_GetAnimationRibbonGraphLookupsBufferPointer, not the library's table.
* @param animation_manager        A handle to a tfx_animation_manager
* @param destination              At least tfx_GetAnimationRibbonEmittersSizeInBytes of writable memory
*/
tfxAPI void tfx_CopyAnimationRibbonEmitters(tfx_animation_manager animation_manager, void *destination);

/*
Create a zeroed tfx_ribbon_dispatch_t to walk an animation manager's ribbon buckets with tfx_NextAnimationRibbonDispatch. Create a fresh one
for each walk.
* @returns tfx_ribbon_dispatch_t  A zeroed dispatch
*/
tfxAPI tfx_ribbon_dispatch_t tfx_CreateAnimationRibbonDispatch(void);

/*
Advance to the next of this frame's baked ribbon buckets that has ribbons to draw. Call in a while loop until it returns false. Fills the
dispatch the same way tfx_NextRibbonDispatch does: total_segments for the tfx_ribbon_compute dispatch, index_count and index_offset for the draw,
and the bucket's push constants through tfx_GetRibbonDispatchGlobals. ribbon_data is null, there's no stage bucket behind it.
* @param animation_manager        A handle to a tfx_animation_manager, updated this frame
* @param ribbon_dispatch          A dispatch from tfx_CreateAnimationRibbonDispatch, updated in place
* @returns bool                   true while there's another bucket
*/
tfxAPI bool tfx_NextAnimationRibbonDispatch(tfx_animation_manager animation_manager, tfx_ribbon_dispatch_t *ribbon_dispatch);

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
Prepare a tfx_effect_template_t that you can use to customise effects in the library in various ways before adding them into a stage for updating and rendering. Using a template like this
means that you can tweak an effect without editing the base effect in the library.
Changes made through the template setters (tfx_SetTemplate*, tfx_ScaleTemplate*, tfx_Enable/DisableTemplateEmitter, tfx_SetEffectTemplateWarmupTime)
are re-applied when tfx_RefreshLibrary rebuilds the template. Edits made directly on descriptors fetched from the template are not.
* @param library                    A reference to a tfx_library_t that should be loaded with tfx_LoadEffectLibrary
* @param name                       The name of the effect in the library that you want to use for the template. If the effect is in a folder then use normal pathing: "My Folder/My effect"
//Returns handle					Handle to the newly created effect template or nullptr if the effect couldn't be found in the library, which prints a tip
*/
tfxAPI tfx_effect_template tfx_CreateEffectTemplate(tfx_library library, const char *name);
 
/*
Delete an effect template and free all memory associated with it
* @param effect_template            A handle to the effect template to be deleted
*/
tfxAPI void tfx_FreeEffectTemplate(tfx_effect_template effect_template);

/*
Check whether the effect a template was cloned from has gone from the library file, as found by tfx_RefreshLibrary.
While this is true, tfx_AddEffectTemplateToStage refuses the template and returns an id that fails tfx_EffectIDIsValid, and the refresh that found
the effect missing expires anything already running from it. The template is not freed: if a later refresh finds the
effect in the file again, this goes back to false and the template can be used again. Call tfx_FreeEffectTemplate if you
no longer need it.
* @param effect_template    A handle to the effect template
* @returns bool             True if the effect it was cloned from is no longer in the library file
*/
tfxAPI bool tfx_EffectTemplateIsMarkedForDeletion(tfx_effect_template effect_template);

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
* @param path                     A path to the emitter: "Explosion/Smoke", or the emitter's tfx_GetDescriptorPath
* @returns                        The template's copy of the emitter, or nullptr when the path isn't in the template
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
* @param path                     A path to the effect or emitter in the effect template: "Effect/Emitter", or its tfx_GetDescriptorPath
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
Set an update callback for the root effect in the effect template. The callback runs during the stage update, which may be on the
stage's update thread, so it can't add effects to the same stage.
* @param t                        A pointer to a tfx_effect_template_t
* @param update_callback          A pointer to the call back function
*/
tfxAPI void tfx_SetTemplateEffectUpdateCallback(tfx_effect_template t, void(*update_callback)(tfx_stage stage, tfxEffectID effect_id));

/*
Disable an emitter within an effect. Disabling an emitter will stop it being added to the stage when calling tfx_AddEffectTemplateToStage
* @param path					  const *char of a path to the emitter in the effect, "My Effect/My Emitter" or the emitter's tfx_GetDescriptorPath
*/
tfxAPI void tfx_DisableTemplateEmitter(tfx_effect_template t, const char *path);

/*
Enable an emitter within an effect so that it is added to the stage when calling tfx_AddEffectTemplateToStage. Emitters are enabled by default.
* @param path					  const *char of a path to the emitter in the effect, "My Effect/My Emitter" or the emitter's tfx_GetDescriptorPath
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
Set how many particles a single shot emitter spawns. Single shot emitters ignore the amount graphs, so this is how you change
how many they spawn. Only affects emitters that have the single spawn flag set. A Run on GPU emitter's ring follows the
new amount, times any location count from tfx_SetTemplateGPUCapacityForLocations, unless tfx_SetTemplateEmitterGPUCapacity
fixed its capacity. Rings are sized when an effect is added, so set it before adding the template and before summing
gpu_particle_capacity.
* @param emitter_path			 "Effect/Emitter", or the emitter's tfx_GetDescriptorPath
* @param amount					 The number of particles to spawn.
*/
tfxAPI void tfx_SetTemplateSingleSpawnAmount(tfx_effect_template t, const char *emitter_path, tfxU32 amount);

/*
Set the ring capacity of a Run on GPU emitter, replacing the authored one from tfx_GetEmitterGPUCapacity. Use it to
make room for many spawn locations or a large spawn multiplier. Only affects effects added after the call. Ribbons are ignored,
so every child of tfx_GetEffectFromTemplate can be passed without checking its type. tfx_SetTemplateGPUCapacityForLocations replaces
what this set on every emitter, so to size one emitter differently from the rest call this after it, then take the template's pool
slots from tfx_GetTemplateGPUPoolSlots rather than the earlier return value.
* @param emitter_path			 "Effect/Emitter", or the emitter's tfx_GetDescriptorPath
* @param capacity				 Particles the ring can hold, 0 to go back to the authored capacity times any location count.
*/
tfxAPI void tfx_SetTemplateEmitterGPUCapacity(tfx_effect_template t, const char *emitter_path, tfxU32 capacity);

/*
Give a Run on GPU emitter a key that its dispatches carry in tfx_gpu_particle_dispatch_t key, so a host can choose which kernel runs them,
for example to give one emitter a custom update or spawn hook and leave the rest on the reference kernel. The library never reads it.
Two templates made from the same effect can key the same emitter differently. tfx_SetEffectGPUKey overrides it for one effect.
Only affects effects added after the call. Ribbons are ignored.
* @param emitter_path			 "Effect/Emitter", or the emitter's tfx_GetDescriptorPath
* @param key					 Any value the host chooses, 0 by default.
*/
tfxAPI void tfx_SetTemplateEmitterGPUKey(tfx_effect_template t, const char *emitter_path, tfxU32 key);

/*
Scale all nodes on an emitter graph. Scaling a color graph gives the emitter a color ramp row and particle properties entry of its
own, and scaling an overtime graph the GPU reads gives it its own graph lookups, so tfx_GetColorRampsVersion,
tfx_GetParticlePropertiesVersion or tfx_GetGPUGraphLookupsVersion will tell you to upload them again. Single shot emitters
ignore the amount graphs, use tfx_SetTemplateSingleSpawnAmount for those.
A particle's spawn size is a base plus a random amount up to a variation, each from its own graph, so to make particles
twice the size scale tfxEmitter_base_width_index, tfxEmitter_variation_width_index, tfxEmitter_base_height_index and
tfxEmitter_variation_height_index by 2. An emitter with uniform size takes its height from the two width graphs, so scale
just those. Scaling the base alone leaves the variation where it was.
* @param tfx_effect_template	  Handle to a valid effect template
* @param emitter_path			  "Effect/Emitter", or the emitter's tfx_GetDescriptorPath
* @param graph_index			  tfx_emitter_graph_index of the graph that you want to scale. Must be an emitter graph or an assert will be called
* @param scale_amount             The amount to scale the graph by.
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

//Unpack the width and height of a tfx_instance_t
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
Get the number of shapes in the GPU Shape Data buffer contained within a library. It's empty until the shape data is built (tfx_UpdateLibraryGPUImageData or tfx_BuildAnimationManagerGPUShapeData)
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
