//#define tfxENABLE_PROFILING
#define tfxTRACK_MEMORY
/*
	Timeline FX C++ library

	This library is for implementing particle effects into your games and applications.

	This library is render agnostic, so you will have to provide your own means to render the particles. You will use ParticleManager::GetParticleBuffer() to get all of the active particles in the particle manager
	and then use the values in Particle struct to draw a correctly scaled and rotated particle. See example below.

	//Before you load an effects file, you will need to define a ShapeLoader function that passes the following parameters:
	//const char* filename			- this will be the filename of the image being loaded from the library. You don't have to do anything with this if you don't need to.
	//ImageData	&image_data			- A struct containing data about the image. You will have to set image_data.ptr to point to the texture in your renderer for later use in the Render function that you will create to render the particles
	//void *raw_image_data			- The raw data of the image which you can use to load the image into graphics memory
	//int image_memory_size			- The size in bytes of the raw_image_data
	//void *custom_data				- This allows you to pass through an object you can use to access whatever is necessary to load the image into graphics memory, depending on the renderer that you're using
	void ShapeLoader(const char* filename, ImageData &image_data, void *raw_image_data, int image_memory_size, void *custom_data) {
		//Cast your custom data, this can be anything you want
		TfxExample *example = static_cast<TfxExample*>(custom_data);

		//This shape loader example uses the STB image library to load the raw bitmap (png usually) data
		StbImage image;
		LoadStbImageMemory(image, (unsigned char*)raw_image_data, image_memory_size);
		//Convert the image to RGBA which is necessary for this particular renderer
		ConvertToRGBA(&image);
		//The editor has the option to convert an image to an alpha map. I will probably change this so that it gets baked into the saved effect so you won't need to apply the filter here.
		//Alpha map is where all color channels are set to 255
		if (image_data.import_filter)
			ConvertToAlpha(image);

		//You'll probably need to load the image in such a way depending on whether or not it's an animation or not
		if (image_data.animation_frames > 1) {
			//Add the spritesheet to the texture in our renderer
			tfxU32 anim_index = example->particle_textures->AddAnimation(image, (tfxU32)image_data.image_size.x, (tfxU32)image_data.image_size.y, (tfxU32)image_data.animation_frames);
			//Important step: you need to point the ImageData.ptr to the appropriate handle in the renderer to point to the texture of the particle shape
			//You'll need to use this in your render function to tell your renderer which texture to use to draw the particle
			image_data.ptr = &example->particle_textures->GetAnimation(anim_index);
		}
		else {
			//Add the image to the texture in our renderer
			tfxU32 image_index = example->particle_textures->AddImage(image);
			//Important step: you need to point the ImageData.ptr to the appropriate handle in the renderer to point to the texture of the particle shape
			//You'll need to use this in your render function to tell your renderer which texture to use to draw the particle
			image_data.ptr = &example->particle_textures->GetImage(image_index);
		}
	}

	//Here's an example of a render function that you will need to write in order to integrate timeline fx with your specific renderer that you're using
	//I should think that you could quite easily multi-thread this as well, as long as your renderer is happy with that
	void TfxExample::RenderParticles(float tween) {
		//In this example, a compute shader is used to transform all the vertices into the right place by sending a batch of quads. A quad just has the size, orientation, color and UV coords, the compute
		//shader then builds the vertex buffer by doing all the transforms to save the CPU having to do it.
		render_layer->StartQuadBatch(&particle_textures->PipelineIndex(qulkan::BlendMode::Alpha, 1));

		//Loop through all the draw layers - particles can be assigned to a specific draw layer so you can draw them in a specific order if necessary. The layer is set in the editor on the properties tab.
		for (tfxU32 layer = 0; layer != tfxLAYERS; ++layer) {
			//Use GetParticleBuffer(layer) to get all of the particles in the current layer
			for (auto p : *pm.GetParticleBuffer(layer)) {
				//In order to set the correct blendmode we need to get the property from the parent emitter that emitted the particle
				//A pointer to the parent emitter is stored in the parent member
				tfx::tfxEffectEmitter &e = *p.parent;

				//Set the correct blendmode, see timelinefx::tfxBlendMode. You may have to map the blendmodes depending on the renderer you use
				render_layer->SetBlendMode(qulkan::BlendMode(e.properties.blend_mode));
				//Set the color for the quad
				render_layer->SetColor(p.color.r, p.color.g, p.color.b, p.color.a);
				//Set the value that the color will be multiplied by, this happens in your fragment shader. You can always omit this if you're not using intensity
				render_layer->SetMultiplyFactor(p.intensity);
				//You can use render tweening to smooth particle movement from frame to frame by interpolating between captured and world states
				tfx::FormState tweened = tfx::Tween(tween, p.world, p.captured);
				//Is the particle using an image with more than one frame of animation?
				if (e.properties.image->animation_frames == 1 && e.properties.start_frame == 0) {
					//One frame of animation
					//Set the image handle, this the offset that the particle is drawn at.
					//This is where you can make use of the image->ptr from the ShapeLoader function, cast it into the appropriate type for the renderer
					qulkan::SetImageHandle(*static_cast<qulkan::QulkanImage*>(e.properties.image->ptr), p.handle.x, p.handle.y);
					//Add the particle image quad to the renderer for the next render pass at the particle position/rotation/scale
					render_layer->AddQuad(*static_cast<qulkan::QulkanImage*>(e.properties.image->ptr), tweened.position.x, tweened.position.y, tweened.rotation, tweened.scale.x, tweened.scale.y);
				}
				else {
					//Multiple frames of animation
					//Get the current frame of animation that the particle is using
					uint32_t frame = uint32_t(p.image_frame);
					//frame must be within the bounds of the animation
					assert(frame >= 0 && frame < e.properties.image->animation_frames);
					//Cast the image->ptr to the appropriate type for the renderer to get at the animation frames. In this case it's an AnimationFrames struct which contains a list of indexes
					//where to lookup each texture representing each frame of animation
					qulkan::AnimationFrames* af = static_cast<qulkan::AnimationFrames*>(e.properties.image->ptr);
					//Get the image of the current frame of animation
					qulkan::QulkanImage &image = particle_textures->GetImage(af->images[frame]);

					//Set the image handle
					SetImageHandle(image, p.handle.x, p.handle.y);
					//Add the particle frame of animation quad to the renderer for the next render pass at the particle position/rotation/scale
					render_layer->AddQuad(image, tweened.position.x, tweened.position.y, tweened.rotation, tweened.scale.x, tweened.scale.y);
				}
			}
		}
	}
*/
#pragma once

#if defined(_WIN32)
#include <SDKDDKVer.h>
#define WIN_LEAN_AND_MEAN
#endif

//Might possibly replace some of these in the future
#include <stdio.h>
#include <stdarg.h>					//va_list
#include <chrono>					//std::chrono::high_resolution_clock
#include <cctype>					//std::is_digit
#include <algorithm>
#include <stdint.h>
#include <assert.h>
#include <iostream>					//temp for std::cout
#include <immintrin.h>
#include <intrin.h>
#include <mutex>

namespace tfx {

#define TWO63 0x8000000000000000u 
#define TWO64f (TWO63*2.0)
#define tfxPI 3.14159265359f
#define tfx360Radians 6.28319f
#define tfx180Radians 3.14159f
#define tfx90Radians 1.5708f

	//----------------------------------------------------------
	//Forward declarations

	struct tfxEffectEmitter;
	struct tfxEffectTemplate;
	struct tfxParticle;
	struct tfxParticleData;
	struct tfxComputeSprite;
	struct tfxComputeParticle;
	struct tfxAnimationSettings;
	struct tfxEffectLibrary;
	struct tfxStr;
	struct tfxStr16;
	struct tfxStr32;
	struct tfxStr64;
	struct tfxStr128;
	struct tfxStr256;
	struct tfxStr512;

//--------------------------------------------------------------
//macros
#define TFX_VERSION "Alpha"
#define TFX_VERSION_NUMBER 3.29.2022

#define tfxMAX_FRAME 20000.f
#define tfxNullParent 0xFFFFFFFF
#define tfxINVALID 0xFFFFFFFF
#define tfxEmitterPropertiesCount 26

#define tfxDel << "=" <<
#define tfxCom << "," <<
#define tfxEndLine << std::endl

#define tfxDelt "=" 
#define tfxComt ","
#define tfxEndLinet "\n"

typedef std::chrono::high_resolution_clock tfxClock;

//Override this for more layers, although currently the editor is fixed at 4
#ifndef tfxLAYERS
#define tfxLAYERS 4
#define tfxEachLayer int layer = 0; layer !=tfxLAYERS; ++layer
#endif 
#ifndef tfxDEFAULT_SPRITE_ALLOCATION
#define tfxDEFAULT_SPRITE_ALLOCATION 25000
#endif
//type defs
typedef unsigned int tfxU32;
typedef int tfxS32;
typedef unsigned long long tfxU64;
typedef long long tfxS64;
typedef tfxU32 tfxEffectID;
typedef unsigned long long tfxKey;
union tfxUInt10bit
{
	struct
	{
		int x : 10;
		int y : 10;
		int z : 10;
		int w : 2;
	} data;
	tfxU32 pack;
};


	//----------------------------------------------------------
	//enums/flags

	//Blend mode property of the emitter
	//It's up to whoever is implementing this library to provide a render function for the particles and make use of these blend modes
	enum tfxBlendMode : unsigned char {
		tfxNone = 0,												//Basically not used, only alpha and additive are used
		tfxAlpha = 1,												//Alpha blend the particle with what it's being drawn on 
		tfxAdditive = 2												//Add the color of the particle with what it's being drawn on
	};

	//tfxGraph presets to determine limits and scales of different graphs, mainly used for the editor
	enum tfxGraphPreset {
		tfxGlobalPercentPreset,
		tfxGlobalOpacityPreset,
		tfxGlobalPercentPresetSigned,
		tfxAnglePreset,
		tfxArcPreset,
		tfxEmissionRangePreset,
		tfxDimensionsPreset,
		tfxLifePreset,
		tfxAmountPreset,
		tfxVelocityPreset,
		tfxVelocityOvertimePreset,
		tfxWeightPreset,
		tfxWeightVariationPreset,
		tfxNoiseOffsetVariationPreset,
		tfxNoiseResolutionPreset,
		tfxWeightOvertimePreset,
		tfxSpinPreset,
		tfxSpinVariationPreset,
		tfxSpinOvertimePreset,
		tfxDirectionOvertimePreset,
		tfxDirectionVariationPreset,
		tfxFrameratePreset,
		tfxVelocityTurbulancePreset,
		tfxOpacityOvertimePreset,
		tfxColorPreset,
		tfxPercentOvertime,
		tfxIntensityOvertimePreset
	};

	enum tfxGraphCategory : unsigned int {
		tfxGraphCategory_global,
		tfxGraphCategory_property,
		tfxGraphCategory_base,
		tfxGraphCategory_variation,
		tfxGraphCategory_overtime
	};


#define tfxGlobalCount  18
#define	tfxPropertyCount  12
#define	tfxBaseCount  8
#define	tfxVariationCount  9
#define	tfxOvertimeCount  16

#define tfxGlobalStart 0
#define	tfxPropertyStart tfxGlobalCount
#define	tfxBaseStart (tfxPropertyStart + tfxPropertyCount)
#define	tfxVariationStart (tfxBaseStart + tfxBaseCount)
#define	tfxOvertimeStart (tfxVariationStart + tfxVariationCount)

	//All the different types of graphs, split into main type: global, property, base, variation and overtime
	enum tfxGraphType : unsigned char {
		tfxGlobal_life,
		tfxGlobal_amount,
		tfxGlobal_velocity,
		tfxGlobal_width,
		tfxGlobal_height,
		tfxGlobal_weight,
		tfxGlobal_spin,
		tfxGlobal_stretch,
		tfxGlobal_overal_scale,
		tfxGlobal_intensity,
		tfxGlobal_frame_rate,
		tfxGlobal_splatter,
		tfxGlobal_effect_roll,
		tfxGlobal_effect_pitch,
		tfxGlobal_effect_yaw,
		tfxGlobal_emitter_width,
		tfxGlobal_emitter_height,
		tfxGlobal_emitter_depth,

		tfxProperty_emission_pitch,
		tfxProperty_emission_yaw,
		tfxProperty_emission_range,
		tfxProperty_emitter_roll,
		tfxProperty_emitter_pitch,
		tfxProperty_emitter_yaw,
		tfxProperty_splatter,
		tfxProperty_emitter_width,
		tfxProperty_emitter_height,
		tfxProperty_emitter_depth,
		tfxProperty_arc_size,
		tfxProperty_arc_offset,

		tfxBase_life,
		tfxBase_amount,
		tfxBase_velocity,
		tfxBase_width,
		tfxBase_height,
		tfxBase_weight,
		tfxBase_spin,
		tfxBase_noise_offset,

		tfxVariation_life,
		tfxVariation_amount,
		tfxVariation_velocity,
		tfxVariation_width,
		tfxVariation_height,
		tfxVariation_weight,
		tfxVariation_spin,
		tfxVariation_noise_offset,
		tfxVariation_noise_resolution,

		tfxOvertime_velocity,
		tfxOvertime_width,
		tfxOvertime_height,
		tfxOvertime_weight,
		tfxOvertime_spin,
		tfxOvertime_stretch,
		tfxOvertime_red,
		tfxOvertime_green,
		tfxOvertime_blue,
		tfxOvertime_blendfactor,
		tfxOvertime_velocity_turbulance,
		tfxOvertime_direction_turbulance,
		tfxOvertime_velocity_adjuster,
		tfxOvertime_intensity,
		tfxOvertime_direction,
		tfxOvertime_noise_resolution,
		tfxGraphMaxIndex,
	};

	//tfxEffectEmitter type - effect contains emitters, and emitters spawn particles, but they both share the same struct for simplicity
	enum tfxEffectEmitterType : unsigned char {
		tfxEffectType,
		tfxEmitterType,
		tfxStage,
		tfxFolder
	};

	//Different ways that particles can be emitted
	enum tfxEmissionType : unsigned char {
		tfxPoint,
		tfxArea,
		tfxLine,
		tfxEllipse
	};

	//Determines how for area, line and ellipse emitters the direction that particles should travel
	enum tfxEmissionDirection : unsigned char {
		tfxInwards,
		tfxOutwards,
		tfxBothways,
		tfxSpecified
	};

	//For line effects where traverse line is switched on
	enum tfxLineTraversalEndBehaviour : unsigned char {
		tfxLoop,
		tfxKill,
		tfxLetFree
	};

	//Mainly for the editor, maybe this can just be moved there instead?
	enum tfxExportColorOptions {
		tfxFullColor,
		tfxOneColor,
		tfxGreyScale,
		tfxOneColorAlpha,
		tfxGreyScaleAlhpa
	};

	//Mainly for the editor, maybe this can just be moved there instead?
	enum tfxExportOptions {
		tfxSpriteSheet,
		tfxStrip,
		tfxSeparateImages
	};

	//tfxGraph data can be looked up in one of 2 ways, either by just using linear/bezier interpolation (slower), or look up the value in a pre-compiled look up table.
	enum tfxLookupMode {
		tfxPrecise,
		tfxFast
	};

	//Used in file loading - for loading effects library
	enum tfxDataType {
		tfxString,
		tfxSInt,
		tfxUint,
		tfxFloat,
		tfxDouble,
		tfxBool
	};
	
	//Block designators for loading effects library
	enum tfxEffectLibraryStream : uint32_t {
		tfxStartEffect = 0x00FFFF00,
		tfxEndEffect,
		tfxStartEmitter,
		tfxEndEmitter,
		tfxStartGraphs,
		tfxEndGraphs,
		tfxStartShapes,
		tfxEndShapes,
		tfxStartAnimationSettings,
		tfxEndAnimationSettings,
		tfxStartImageData,
		tfxStartEffectData,
		tfxEndOfFile,
		tfxStartFolder,
		tfxEndFolder,
		tfxStartPreviewCameraSettings,
		tfxEndPreviewCameraSettings,
	};

	typedef tfxU32 tfxEmitterPropertyFlags;
	typedef tfxU32 tfxEffectPropertyFlags;
	typedef tfxU32 tfxVectorFieldFlags;
	typedef unsigned char tfxParticleFlags;
	typedef tfxU32 tfxEmitterStateFlags;
	typedef tfxU32 tfxParticleControlFlags;
	typedef tfxU32 tfxAttributeNodeFlags;
	typedef tfxU32 tfxAngleSettingFlags;
	typedef tfxU32 tfxEffectManagerFlags;
	typedef tfxU32 tfxErrorFlags;
	typedef tfxU32 tfxEffectCloningFlags;

	enum tfxErrorFlags_ {
		tfxErrorCode_success = 0,
		tfxErrorCode_incorrect_package_format = 1 << 0,
		tfxErrorCode_data_could_not_be_loaded = 1 << 1,
		tfxErrorCode_could_not_add_shape = 1 << 2,
		tfxErrorCode_error_loading_shapes = 1 << 3,
		tfxErrorCode_some_data_not_loaded = 1 << 4,
		tfxErrorCode_unable_to_open_file = 1 << 5,
		tfxErrorCode_unable_to_read_file = 1 << 6,
		tfxErrorCode_wrong_file_size = 1 << 7,
		tfxErrorCode_invalid_format = 1 << 8,
		tfxErrorCode_no_inventory = 1 << 9,
		tfxErrorCode_invalid_inventory = 1 << 10

	};

	enum tfxEffectCloningFlags_ {
		tfxEffectCloningFlags_none = 0,
		tfxEffectCloningFlags_keep_user_data = 1 << 0,
		tfxEffectCloningFlags_force_clone_global = 1 << 1,
		tfxEffectCloningFlags_clone_graphs = 1 << 2,
		tfxEffectCloningFlags_compile_graphs = 1 << 3
	};

	enum tfxParticleManagerModes {
		tfxParticleManagerMode_unordered,
		tfxParticleManagerMode_ordered_by_age,
		tfxParticleManagerMode_ordered_by_depth,
		tfxParticleManagerMode_ordered_by_depth_guaranteed
	};

	enum tfxBillboardingOptions {
		tfxBillboarding = 0,
		tfxBillboarding_disabled = 1,
		tfxBillboarding_disabled_align = 2,
		tfxBillboarding_align = 1 << 2 
	};

	enum tfxEffectManagerFlags_ {
		tfxEffectManagerFlags_none = 0,
		tfxEffectManagerFlags_disable_spawning = 1,
		tfxEffectManagerFlags_force_capture = 2,
		tfxEffectManagerFlags_use_compute_shader = 1 << 3,
		tfxEffectManagerFlags_order_by_depth = 1 << 4,
		tfxEffectManagerFlags_guarantee_order = 1 << 5,
		tfxEffectManagerFlags_update_base_values = 1 << 6,
		tfxEffectManagerFlags_dynamic_sprite_allocation = 1 << 7,
		tfxEffectManagerFlags_3d_effects = 1 << 8,
		tfxEffectManagerFlags_unorderd = 1 << 9,
		tfxEffectManagerFlags_ordered_by_age = 1 << 10
	};

	enum tfxVectorAlignType {
		tfxVectorAlignType_motion,
		tfxVectorAlignType_emission,
		tfxVectorAlignType_emitter,
		tfxVectorAlignType_max,
		//Not in yet, need to think about methods of implementing
		tfxVectorAlignType_surface_normal,
	};

	//Particle property that defines how a particle will rotate
	enum tfxAngleSettingFlags_ {
		tfxAngleSettingFlags_none = 0,														//No flag
		tfxAngleSettingFlags_align_roll = 1 << 0,											//Align the particle with it's direction of travel in 2d
		tfxAngleSettingFlags_random_roll = 1 << 1,											//Chose a random angle at spawn time/flags
		tfxAngleSettingFlags_specify_roll = 1 << 2,											//Specify the angle at spawn time
		tfxAngleSettingFlags_align_with_emission = 1 << 3,									//Align the particle with the emission direction only
		tfxAngleSettingFlags_random_pitch = 1 << 4,											//3d mode allows for rotating pitch and yaw when not using billboarding (when particle always faces the camera)
		tfxAngleSettingFlags_random_yaw = 1 << 5,
		tfxAngleSettingFlags_specify_pitch = 1 << 6,
		tfxAngleSettingFlags_specify_yaw = 1 << 7
	};

	//All the flags needed by the ControlParticle function put into one enum to save space
	enum tfxParticleControlFlags_ {
		tfxParticleControlFlags_none = 0,
		tfxParticleControlFlags_random_color = 1 << 0,
		tfxParticleControlFlags_relative_position = 1 << 1,
		tfxParticleControlFlags_relative_angle = 1 << 2,
		tfxParticleControlFlags_point = 1 << 3,
		tfxParticleControlFlags_area = 1 << 4,
		tfxParticleControlFlags_line = 1 << 5,
		tfxParticleControlFlags_ellipse = 1 << 6,
		tfxParticleControlFlags_loop = 1 << 7,
		tfxParticleControlFlags_kill = 1 << 8,
		tfxParticleControlFlags_letFree = 1 << 9,
		tfxParticleControlFlags_edge_traversal = 1 << 10,
		tfxParticleControlFlags_remove = 1 << 11,
		tfxParticleControlFlags_base_uniform_size = 1 << 12,
		tfxParticleControlFlags_lifetime_uniform_size = 1 << 13,
		tfxParticleControlFlags_animate = 1 << 14,
		tfxParticleControlFlags_reverse_animation = 1 << 15,
		tfxParticleControlFlags_play_once = 1 << 16,
		tfxParticleControlFlags_align = 1 << 17,
		tfxParticleControlFlags_emission = 1 << 18,
		tfxParticleControlFlags_random_roll = 1 << 19,
		tfxParticleControlFlags_specify_roll = 1 << 20,
		tfxParticleControlFlags_random_pitch = 1 << 21,
		tfxParticleControlFlags_specify_pitch = 1 << 22,
		tfxParticleControlFlags_random_yaw = 1 << 23,
		tfxParticleControlFlags_specify_yaw = 1 << 24,
	};

	enum tfxEffectPropertyFlags_ {
		tfxEffectPropertyFlags_none = 0,
		tfxEffectPropertyFlags_is_3d = 1 << 0,
		tfxEffectPropertyFlags_depth_draw_order = 1 << 1,
		tfxEffectPropertyFlags_guaranteed_order = 1 << 2,
		tfxEffectPropertyFlags_age_order = 1 << 3,
	};

	enum tfxEmitterPropertyFlags_ {
		tfxEmitterPropertyFlags_none = 0,
		tfxEmitterPropertyFlags_random_color = 1 << 0,						//Pick a random color from the color overtime gradient rather then change the color over the lifetime of the particle
		tfxEmitterPropertyFlags_relative_position = 1 << 1,					//Keep the particles position relative to the current position of the emitter
		tfxEmitterPropertyFlags_relative_angle = 1 << 2,					//Keep the angle of the particles relative to the current angle of the emitter
		tfxEmitterPropertyFlags_image_handle_auto_center = 1 << 3,			//Set the offset of the particle to the center of the image
		tfxEmitterPropertyFlags_single = 1 << 4,							//Only spawn a single particle (or number of particles specified by spawn_amount) that does not expire
		tfxEmitterPropertyFlags_specific_emission_direction = 1 << 5,		//Uses a normal vector (3d) or direction (2d) to determine emission direction
		tfxEmitterPropertyFlags_spawn_on_grid = 1 << 6,						//When using an area, line or ellipse emitter, spawn along a grid
		tfxEmitterPropertyFlags_grid_spawn_clockwise = 1 << 7,				//Spawn clockwise/left to right around the area
		tfxEmitterPropertyFlags_fill_area = 1 << 8,							//Fill the area
		tfxEmitterPropertyFlags_emitter_handle_auto_center = 1 << 9,		//Center the handle of the emitter
		tfxEmitterPropertyFlags_edge_traversal = 1 << 10,					//Line emitters only: make particles traverse the line
		tfxEmitterPropertyFlags_global_uniform_size = 1 << 11,				//Keep the global particle size uniform
		tfxEmitterPropertyFlags_base_uniform_size = 1 << 12,				//Keep the base particle size uniform
		tfxEmitterPropertyFlags_lifetime_uniform_size = 1 << 13,			//Keep the size over lifetime of the particle uniform
		tfxEmitterPropertyFlags_animate = 1 << 14,							//Animate the particle shape if it has more than one frame of animation
		tfxEmitterPropertyFlags_reverse_animation = 1 << 15,				//Make the image animation go in reverse
		tfxEmitterPropertyFlags_play_once = 1 << 16,						//Play the animation once only
		tfxEmitterPropertyFlags_random_start_frame = 1 << 17,				//Start the animation of the image from a random frame
		tfxEmitterPropertyFlags_keep_alive = 1 << 18,						//Keep the effect/emitter in the particle manager, don't remove it when it has no particles
		tfxEmitterPropertyFlags_unused = 1 << 19,							//Unused
		tfxEmitterPropertyFlags_is_in_folder = 1 << 20,						//This effect is located inside a folder
		tfxEmitterPropertyFlags_is_bottom_emitter = 1 << 21,				//This emitter has no child effects, so can spawn particles that could be used in a compute shader if it's enabled
		tfxEmitterPropertyFlags_use_spawn_ratio = 1 << 22,					//Option for area emitters to multiply the amount spawned by a ration of particles per pixels squared
		tfxEmitterPropertyFlags_can_grow_particle_memory = 1 << 23,			//Allows for expanding the memory used for particle emitters if the amount spawned is changed dynamically
		tfxEmitterPropertyFlags_is_3d = 1 << 24,							//Makes the effect run in 3d mode for 3d effects
		tfxEmitterPropertyFlags_use_dynamic = 1 << 25,						//Use a dynamic particle storage rather then a fixed one
	};

	enum tfxParticleFlags_ : unsigned char {
		tfxParticleFlags_none = 0,
		tfxParticleFlags_fresh = 1 << 0,									//Particle has just spawned this frame	
		tfxParticleFlags_capture_after_transform = 1 << 3,					//Particle will be captured after a transfrom, used for traversing lines and looping back to the beginning to avoid lerping imbetween
		tfxParticleFlags_remove = 1 << 4,									//Particle will be removed this or next frame
		tfxParticleFlags_has_velocity = 1 << 5,								//Flagged if the particle is currently moving
	};

	enum tfxEmitterStateFlags_ : unsigned int {
		tfxEmitterStateFlags_none = 0,
		tfxEmitterStateFlags_random_color = 1 << 0,
		tfxEmitterStateFlags_relative_position = 1 << 1,					//Keep the particles position relative to the current position of the emitter
		tfxEmitterStateFlags_relative_angle = 1 << 2,						//Keep the angle of the particles relative to the current angle of the emitter
		tfxEmitterStateFlags_stop_spawning = 1 << 3,						//Tells the emitter to stop spawning
		tfxEmitterStateFlags_remove = 1 << 4,								//Tells the effect/emitter to remove itself from the particle manager immediately
		tfxEmitterStateFlags_enabled = 1 << 5,								//the emitter is enabled. If flag is not set then it will not be added to the particle manager with AddEffect
		tfxEmitterStateFlags_retain_matrix = 1 << 6,						//Internal flag about matrix usage
		tfxEmitterStateFlags_no_tween_this_update = 1 << 7,					//Internal flag generally, but you could use it if you want to teleport the effect to another location
		tfxEmitterStateFlags_is_single = 1 << 8,
		tfxEmitterStateFlags_not_line = 1 << 9,
		tfxEmitterStateFlags_is_line_traversal = 1 << 10,
		tfxEmitterStateFlags_can_spin = 1 << 11,
		tfxEmitterStateFlags_base_uniform_size = 1 << 12,
		tfxEmitterStateFlags_lifetime_uniform_size = 1 << 13,				//Keep the size over lifetime of the particle uniform
		tfxEmitterStateFlags_loop = 1 << 14,
		tfxEmitterStateFlags_kill = 1 << 15,
		tfxEmitterStateFlags_play_once = 1 << 16,							//Play the animation once only
		tfxEmitterStateFlags_single_shot_done = 1 << 17,
		tfxEmitterStateFlags_is_line = 1 << 18,
		tfxEmitterStateFlags_is_area = 1 << 19,
		tfxEmitterStateFlags_no_tween = 1 << 20,
		tfxEmitterStateFlags_align_with_velocity = 1 << 21,
		tfxEmitterStateFlags_is_sub_emitter = 1 << 28
	};

	enum tfxVectorFieldFlags_: unsigned char {
		tfxVectorFieldFlags_none = 0,
		tfxVectorFieldFlags_repeat_horizontal = 1 << 0,						//Field will repeat horizontally
		tfxVectorFieldFlags_repeat_vertical = 1 << 1						//Field will repeat vertically
	};

	enum tfxAttributeNodeFlags_ {
		tfxAttributeNodeFlags_none = 0,
		tfxAttributeNodeFlags_is_curve = 1 << 0,
		tfxAttributeNodeFlags_is_left_curve = 2 << 1,
		tfxAttributeNodeFlags_is_right_curve = 3 << 2,
		tfxAttributeNodeFlags_curves_initialised = 4 << 3
	};

	//-----------------------------------------------------------
	//Constants

	const float tfxMIN_FLOAT = -2147483648.f;
	const float tfxMAX_FLOAT = 2147483647.f;
	const tfxU32 tfxMAX_UINT = 4294967295;

	const float tfxLIFE_MIN = 0.f;
	const float tfxLIFE_MAX = 100000.f;
	const float tfxLIFE_STEPS = 200.f;

	const float tfxAMOUNT_MIN = 0.f;
	const float tfxAMOUNT_MAX = 5000.f;
	const float tfxAMOUNT_STEPS = 100.f;

	const float tfxGLOBAL_PERCENT_MIN = 0.f;
	const float tfxGLOBAL_PERCENT_MAX = 20.f;
	const float tfxGLOBAL_PERCENT_STEPS = 100.f;

	const float tfxGLOBAL_PERCENT_V_MIN = 0.f;
	const float tfxGLOBAL_PERCENT_V_MAX = 10.f;
	const float tfxGLOBAL_PERCENT_V_STEPS = 200.f;

	const float tfxINTENSITY_MIN = 0.f;
	const float tfxINTENSITY_MAX = 5.f;
	const float tfxINTENSITY_STEPS = 100.f;

	const float tfxANGLE_MIN = 0.f;
	const float tfxANGLE_MAX = 1080.f;
	const float tfxANGLE_STEPS = 54.f;

	const float tfxARC_MIN = 0.f;
	const float tfxARC_MAX = 6.28319f;
	const float tfxARC_STEPS = .3141595f;

	const float tfxEMISSION_RANGE_MIN = 0.f;
	const float tfxEMISSION_RANGE_MAX = 180.f;
	const float tfxEMISSION_RANGE_STEPS = 30.f;

	const float tfxDIMENSIONS_MIN = 0.f;
	const float tfxDIMENSIONS_MAX = 4000.f;
	const float tfxDIMENSIONS_STEPS = 40.f;

	const float tfxVELOCITY_MIN = 0.f;
	const float tfxVELOCITY_MAX = 10000.f;
	const float tfxVELOCITY_STEPS = 100.f;

	const float tfxVELOCITY_OVERTIME_MIN = -20.f;
	const float tfxVELOCITY_OVERTIME_MAX = 20.f;
	const float tfxVELOCITY_OVERTIME_STEPS = 200.f;

	const float tfxWEIGHT_MIN = -2500.f;
	const float tfxWEIGHT_MAX = 2500.f;
	const float tfxWEIGHT_STEPS = 200.f;

	const float tfxWEIGHT_VARIATION_MIN = 0.f;
	const float tfxWEIGHT_VARIATION_MAX = 2500.f;
	const float tfxWEIGHT_VARIATION_STEPS = 250.f;

	const float tfxSPIN_MIN = -2000.f;
	const float tfxSPIN_MAX = 2000.f;
	const float tfxSPIN_STEPS = 100.f;

	const float tfxSPIN_VARIATION_MIN = 0.f;
	const float tfxSPIN_VARIATION_MAX = 2000.f;
	const float tfxSPIN_VARIATION_STEPS = 100.f;

	const float tfxSPIN_OVERTIME_MIN = -20.f;
	const float tfxSPIN_OVERTIME_MAX = 20.f;
	const float tfxSPIN_OVERTIME_STEPS = 200.f;

	const float tfxDIRECTION_OVERTIME_MIN = 0.f;
	const float tfxDIRECTION_OVERTIME_MAX = 4320.f;
	const float tfxDIRECTION_OVERTIME_STEPS = 216.f;

	const float tfxFRAMERATE_MIN = 0.f;
	const float tfxFRAMERATE_MAX = 200.f;
	const float tfxFRAMERATE_STEPS = 100.f;

	const float tfxMAX_DIRECTION_VARIATION = 22.5f;
	const float tfxMAX_VELOCITY_VARIATION = 30.f;
	const int tfxMOTION_VARIATION_INTERVAL = 30;
	
	//these Variables determine the timing resolution that particles are updated at. So an Update frequency of 60 would mean that the particles are updated at 60 frames per second.
	extern float tfxUPDATE_FREQUENCY;
	extern float tfxUPDATE_TIME;
	extern float tfxFRAME_LENGTH;

	//Look up frequency determines the resolution of graphs that are compiled into look up arrays.
	static float tfxLOOKUP_FREQUENCY = 10.f;
	//Overtime frequency is for lookups that will vary in length depending on the lifetime of the particle. It should generally be a higher resolution than the base graphs
	static float tfxLOOKUP_FREQUENCY_OVERTIME = 1.f;

	//-----------------------------------------------------------
	//Utility things:
	
	struct tfxMemoryTrackerEntry {
		char name[64];
		tfxU64 amount_allocated;
		void *address;
		tfxU32 allocations;
		bool is_alive;
	};

	struct tfxMemoryTrackerPair {
		tfxKey key;
		tfxMemoryTrackerEntry log;
		tfxMemoryTrackerPair(tfxKey k, tfxMemoryTrackerEntry l) : key(k), log(l) {}
	};

	struct tfxLogList {
		tfxU32 current_size;
		tfxU32 capacity;
		tfxMemoryTrackerPair* data;

		inline tfxLogList() { current_size = capacity = 0; data = NULL; }
		inline ~tfxLogList() { if (data) free(data); data = NULL; current_size = capacity = 0; }

		inline bool			empty() { return current_size == 0; }
		inline tfxMemoryTrackerPair&           operator[](tfxU32 i) { return data[i]; }
		inline const tfxMemoryTrackerPair&     operator[](tfxU32 i) const { assert(i < current_size); return data[i]; }

		inline void         free_all() { if (data) { current_size = capacity = 0; free(data); data = NULL; } }
		inline void         clear() { if (data) { current_size = 0; } }
		inline tfxMemoryTrackerPair*           begin() { return data; }
		inline const tfxMemoryTrackerPair*     begin() const { return data; }
		inline tfxMemoryTrackerPair*           end() { return data + current_size; }
		inline const tfxMemoryTrackerPair*     end() const { return data + current_size; }

		inline tfxU32       _grow_capacity(tfxU32 sz) const { tfxU32 new_capacity = capacity ? (capacity + capacity / 2) : 8; return new_capacity > sz ? new_capacity : sz; }
		inline void         resize(tfxU32 new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
		inline void         reserve(tfxU32 new_capacity) { if (new_capacity <= capacity) return; tfxMemoryTrackerPair* new_data = (tfxMemoryTrackerPair*)malloc((size_t)new_capacity * sizeof(tfxMemoryTrackerPair)); if (data) { memcpy(new_data, data, (size_t)current_size * sizeof(tfxMemoryTrackerPair)); free(data); } data = new_data; capacity = new_capacity; }

		inline tfxMemoryTrackerPair*           erase(const tfxMemoryTrackerPair* it) { assert(it >= data && it < data + current_size); const ptrdiff_t off = it - data; memmove(data + off, data + off + 1, ((size_t)current_size - (size_t)off - 1) * sizeof(tfxMemoryTrackerPair)); current_size--; return data + off; }
		inline tfxMemoryTrackerPair*           erase(const tfxMemoryTrackerPair* it, const tfxMemoryTrackerPair* it_last) { assert(it >= data && it < data + current_size && it_last > it && it_last <= data + current_size); const ptrdiff_t count = it_last - it; const ptrdiff_t off = it - data; memmove(data + off, data + off + count, ((size_t)current_size - (size_t)off - count) * sizeof(tfxMemoryTrackerPair)); current_size -= (tfxU32)count; return data + off; }
		inline tfxMemoryTrackerPair*           insert(const tfxMemoryTrackerPair* it, const tfxMemoryTrackerPair& v) { assert(it >= data && it <= data + current_size); const ptrdiff_t off = it - data; if (current_size == capacity) reserve(_grow_capacity(current_size + 1)); if (off < (ptrdiff_t)current_size) memmove(data + off + 1, data + off, ((size_t)current_size - (size_t)off) * sizeof(tfxMemoryTrackerPair)); new((void*)(data + off)) tfxMemoryTrackerPair(v); current_size++; return data + off; }

	};

	struct tfxMemoryTrackerLog{

		tfxLogList log;
		std::mutex insert_mutex;

		//Insert a new T value into the storage
		void Insert(tfxKey key, const tfxMemoryTrackerEntry &value) {
			SetIndex(key, value);
		}

		void SetIndex(tfxKey key, const tfxMemoryTrackerEntry &value) {
			tfxMemoryTrackerPair* it = LowerBound(key);
			if (it == log.end() || it->key != key)
			{
				std::lock_guard<std::mutex> lock(insert_mutex);
				log.insert(it, tfxMemoryTrackerPair(key, value));
				return;
			} 
			it->log = value;
		}

		inline tfxMemoryTrackerEntry *At(tfxKey key) {
			return GetIndex(key);
		}

		tfxMemoryTrackerEntry *GetIndex(tfxKey key) {
			tfxMemoryTrackerPair* it = LowerBound(key);
			if (it == log.end() || it->key != key)
				return NULL;
			return &it->log;
		}

		tfxMemoryTrackerPair* LowerBound(tfxKey key)
		{
			tfxMemoryTrackerPair* first = log.data;
			tfxMemoryTrackerPair* last = log.data + log.current_size;
			size_t count = (size_t)(last - first);
			while (count > 0)
			{
				size_t count2 = count >> 1;
				tfxMemoryTrackerPair* mid = first + count2;
				if (mid->key < key)
				{
					first = ++mid;
					count -= count2 + 1;
				}
				else
				{
					count = count2;
				}
			}
			return first;
		}

		inline bool ValidKey(tfxKey key) {
			return GetIndex(key) != NULL;
		}
	};

	extern tfxMemoryTrackerLog tfxMEMORY_TRACKER;
	extern char tfxMEMORY_CONTEXT[64];

#ifdef tfxTRACK_MEMORY
#define tfxALLOCATE(tracker_name, dst, size) malloc(size); tfxMemoryTrackerEntry tracker; memset(&tracker, '\0', sizeof(tfxMemoryTrackerEntry)); memcpy(tracker.name, tracker_name, strlen(tracker_name)); tracker.amount_allocated = size; tracker.address = dst; tfxU32 allocations = 0; if (data &&tfxMEMORY_TRACKER.ValidKey((tfxKey)data)) { tfxMemoryTrackerEntry *entry = tfxMEMORY_TRACKER.At((tfxKey)data); if(entry->is_alive) allocations = tfxMEMORY_TRACKER.At((tfxKey)data)->allocations; } tracker.allocations = allocations + 1; tracker.is_alive = true; tfxMEMORY_TRACKER.Insert((tfxKey)dst, tracker);
#define tfxFREE(dst) free(dst); assert(tfxMEMORY_TRACKER.ValidKey((tfxKey)dst)); tfxMEMORY_TRACKER.At((tfxKey)dst)->is_alive = false;
#define tfxINIT_VEC_NAME memset(name, '\0', 64); name[0] = 'X' 
#define tfxINIT_VEC_NAME_INIT memset(name, '\0', 64); memcpy(name, name_init, strlen(name_init));
#define tfxINIT_VEC_NAME_SRC_COPY memset(name, '\0', 64); memcpy(name, src.name, strlen(src.name));
#define tfxCONSTRUCTOR_VEC_DEF const char *name_init
#define tfxCONSTRUCTOR_VEC_INIT(name) #name
#define tfxCONSTRUCTOR_VEC_INIT2(name) name
#else
#define tfxALLOCATE(tracker_name, dst, size) malloc(size);
#define tfxFREE(dst) free(dst); 
#define tfxINIT_VEC_NAME 
#define tfxINIT_VEC_NAME_INIT 
#define tfxINIT_VEC_NAME_SRC_COPY 
#define tfxCONSTRUCTOR_VEC_DEF 
#define tfxCONSTRUCTOR_VEC_INIT(name) 
#define tfxCONSTRUCTOR_VEC_INIT2(name) 
#endif

	//Intrinsics and multithreading

	inline tfxU64 AtomicExchange64(tfxU64 volatile *value, tfxU64 new_value) {
		tfxU64 result = _InterlockedExchange64((__int64*)value, new_value);
		return result;
	}
	inline tfxU64 AtomicAdd64(tfxU64 volatile *value, tfxU64 amount_to_add) {
		tfxU64 result = _InterlockedExchangeAdd64((__int64*)value, amount_to_add);
		return result;
	}

	inline tfxU32 AtomicExchange32(tfxU32 volatile *value, tfxU32 new_value) {
		tfxU32 result = _InterlockedExchange((long*)value, new_value);
		return result;
	}
	inline tfxU32 AtomicAdd32(tfxU32 volatile *value, tfxU32 amount_to_add) {
		tfxU32 result = _InterlockedExchangeAdd((long*)value, amount_to_add);
		return result;
	}

	//Credit to ocornut https://github.com/ocornut/imgui/commits?author=ocornut
	//std::vector replacement with some extra stuff and tweaks specific to Qulkan/TimelineFX
	template<typename T>
	struct tfxvec {
#ifdef tfxTRACK_MEMORY
		char name[64];
#endif
		tfxU32 current_size;
		tfxU32 capacity;
		T* data;

		// Provide standard typedefs but we don't use them ourselves.
		typedef T                   value_type;
		typedef value_type*         iterator;
		typedef const value_type*   const_iterator;

		inline tfxvec() { current_size = capacity = 0; data = NULL; tfxINIT_VEC_NAME;  }
		inline tfxvec(const char *name_init) { current_size = capacity = 0; data = NULL; tfxINIT_VEC_NAME_INIT(name_init);  }
		inline tfxvec(const tfxvec<T> &src) { current_size = capacity = 0; data = NULL; tfxINIT_VEC_NAME; resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); }
		inline tfxvec<T>& operator=( const tfxvec<T>& src) { clear(); resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); tfxINIT_VEC_NAME_SRC_COPY; return *this; }
		inline ~tfxvec() { if (data) { tfxFREE(data) }; data = NULL; current_size = capacity = 0; }

		inline bool			empty() { return current_size == 0; }
		inline tfxU32		size() { return current_size; }
		inline const tfxU32	size() const { return current_size; }
		inline T&           operator[](tfxU32 i) { return data[i]; }
		inline const T&     operator[](tfxU32 i) const { assert(i < current_size); return data[i]; }

		inline void         free_all() { if (data) { current_size = capacity = 0; tfxFREE(data); data = NULL; } }
		inline void         clear() { if (data) { current_size = 0; } }
		inline T*           begin() { return data; }
		inline const T*     begin() const { return data; }
		inline T*           end() { return data + current_size; }
		inline const T*     end() const { return data + current_size; }
		inline T*           rend() { return data; }
		inline const T*     rend() const { return data; }
		inline T*           rbegin() { return data + current_size; }
		inline const T*     rbegin() const { return data + current_size; }
		inline T&           front() { assert(current_size > 0); return data[0]; }
		inline const T&     front() const { assert(current_size > 0); return data[0]; }
		inline T&           back() { assert(current_size > 0); return data[current_size - 1]; }
		inline const T&     back() const { assert(current_size > 0); return data[current_size - 1]; }
		inline T&           parent() { assert(current_size > 1); return data[current_size - 2]; }
		inline const T&     parent() const { assert(current_size > 1); return data[current_size - 2]; }

		inline tfxU32       _grow_capacity(tfxU32 sz) const { tfxU32 new_capacity = capacity ? (capacity + capacity / 2) : 8; return new_capacity > sz ? new_capacity : sz; }
		inline void         resize(tfxU32 new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
		inline void         resize(tfxU32 new_size, const T& v) { if (new_size > capacity) reserve(_grow_capacity(new_size)); if (new_size > current_size) for (tfxU32 n = current_size; n < new_size; n++) memcpy(&data[n], &v, sizeof(v)); current_size = new_size; }
		inline void         shrink(tfxU32 new_size) { assert(new_size <= current_size); current_size = new_size; }
		inline void         reserve(tfxU32 new_capacity) { if (new_capacity <= capacity) return; 
			T* new_data = (T*)tfxALLOCATE(name, new_data, (size_t)new_capacity * sizeof(T)); if (data) { memcpy(new_data, data, (size_t)current_size * sizeof(T)); tfxFREE(data); } data = new_data; capacity = new_capacity; }

		inline T&	        grab() {
			if (current_size == capacity) reserve(_grow_capacity(current_size + 1));
			current_size++;
			return data[current_size - 1];
		}
		inline T&	        push_back(const T& v) { 
			if (current_size == capacity) 
				reserve(_grow_capacity(current_size + 1)); 
			new((void*)(data + current_size)) T(v); 
			current_size++; return data[current_size - 1]; 
		}
		inline T&	        push_back_copy(const T& v) {
			if (current_size == capacity)
				reserve(_grow_capacity(current_size + 1));
			memcpy(&data[current_size], &v, sizeof(v)); 
			current_size++; return data[current_size - 1];
		}
		inline void         pop() { assert(current_size > 0); current_size--; }
		inline T&	        pop_back() { assert(current_size > 0); current_size--; return data[current_size]; }
		inline void         push_front(const T& v) { if (current_size == 0) push_back(v); else insert(data, v); }
		inline T*           erase(const T* it) { assert(it >= data && it < data + current_size); const ptrdiff_t off = it - data; memmove(data + off, data + off + 1, ((size_t)current_size - (size_t)off - 1) * sizeof(T)); current_size--; return data + off; }
		inline T*           erase(const T* it, const T* it_last) { assert(it >= data && it < data + current_size && it_last > it && it_last <= data + current_size); const ptrdiff_t count = it_last - it; const ptrdiff_t off = it - data; memmove(data + off, data + off + count, ((size_t)current_size - (size_t)off - count) * sizeof(T)); current_size -= (tfxU32)count; return data + off; }
		inline T*           erase_unsorted(const T* it) { assert(it >= data && it < data + current_size);  const ptrdiff_t off = it - data; if (it < data + current_size - 1) memcpy(data + off, data + current_size - 1, sizeof(T)); current_size--; return data + off; }
		inline T*           insert(const T* it, const T& v) { assert(it >= data && it <= data + current_size); const ptrdiff_t off = it - data; if (current_size == capacity) reserve(_grow_capacity(current_size + 1)); if (off < (ptrdiff_t)current_size) memmove(data + off + 1, data + off, ((size_t)current_size - (size_t)off) * sizeof(T)); new((void*)(data + off)) T(v); current_size++; return data + off; }
		inline bool         contains(const T& v) const { const T* _data = data;  const T* data_end = data + current_size; while (_data < data_end) if (*_data++ == v) return true; return false; }
		inline T*           find(const T& v) { T* _data = data;  const T* data_end = data + current_size; while (_data < data_end) if (*_data == v) break; else ++_data; return _data; }
		inline const T*     find(const T& v) const { const T* _data = data;  const T* data_end = data + current_size; while (_data < data_end) if (*_data == v) break; else ++_data; return _data; }
		inline bool         find_erase(const T& v) { const T* it = find(v); if (it < data + current_size) { erase(it); return true; } return false; }
		inline bool         find_erase_unsorted(const T& v) { const T* it = find(v); if (it < data + current_size) { erase_unsorted(it); return true; } return false; }
		inline tfxU32       index_from_ptr(const T* it) const { assert(it >= data && it < data + current_size); const ptrdiff_t off = it - data; return (tfxU32)off; }

		inline void			create_pool(tfxU32 amount) { assert(current_size == 0); T base; reserve(amount); for (tfxU32 i = 0; i != capacity; ++i) { new((void*)(data + current_size)) T(base); current_size++; } }
		inline void			create_pool_with(tfxU32 amount, const T &base) { assert(current_size == 0);  reserve(amount); for (tfxU32 i = 0; i != capacity; ++i) { new((void*)(data + current_size)) T(base); current_size++; } }

	};

	//Ring/Circular buffer
	template<typename T>
	struct tfxring {
#ifdef tfxTRACK_MEMORY
		char name[64];
#endif
		T* data;
		unsigned int current_size;
		unsigned int capacity;
		unsigned int start_index;
		int last_bump;

		inline tfxring() { start_index = current_size = capacity = last_bump = 0; data = NULL; tfxINIT_VEC_NAME; }
		inline tfxring(const char *name_init) { start_index = current_size = capacity = last_bump = 0; data = NULL; tfxINIT_VEC_NAME_INIT(name_init);  }
		inline tfxring(unsigned int qty) { start_index = current_size = capacity = last_bump = 0; data = NULL; reserve(qty); tfxINIT_VEC_NAME; }
		inline void         free_all() { if (data) { current_size = capacity = 0; tfxFREE(data); data = NULL; } }

		inline bool			empty() { return current_size == 0; }
		inline bool			full() { return current_size == capacity; }
		inline unsigned int	free_space() { return capacity - current_size; }

		inline void         clear() { start_index = current_size = last_bump = 0; }
		inline unsigned int			size() { return current_size; }
		inline const unsigned int	size() const { return current_size; }
		inline T&           operator[](unsigned int i) { return data[(i + start_index) % capacity]; }
		inline const T&     operator[](unsigned int i) const { return data[(i + start_index) % capacity]; }
		inline T&           AtAbs(unsigned int i) { return data[i]; }
		inline const T&     AtAbs(unsigned int i) const { return data[i]; }

		inline unsigned int end_index() { return (start_index + current_size) % capacity; }
		inline unsigned int before_start_index() { return start_index == 0 ? capacity - 1 : start_index - 1; }
		inline unsigned int last_index() { return (start_index + current_size - 1) % capacity; }
		inline T*			eob_ptr() { return data + (start_index + current_size) % capacity; }
		inline T&           front() { assert(current_size > 0); return data[start_index]; }
		inline const T&     front() const { assert(current_size > 0); return data[start_index]; }
		inline T&           back() { assert(current_size > 0); return data[last_index()]; }
		inline const T&     back() const { assert(current_size > 0); return data[last_index()]; }

		inline void         pop() { assert(current_size > 0); current_size--; }
		inline T&	        pop_back() { assert(current_size > 0); current_size--; return data[current_size]; }

		inline T&	        emplace_back(const T& v) {
			if (current_size == capacity) reserve(_grow_capacity(current_size + 1));
			new((void*)(data + end_index())) T(v);
			current_size++; return data[current_size - 1];
		}
		inline T&	        push_back(const T& v) {
			if (current_size == capacity) reserve(_grow_capacity(current_size + 1));
			memcpy(&data[end_index()], &v, sizeof(v));
			current_size++; return data[current_size - 1];
		}
		inline T&	        grab() {
			if (current_size == capacity) reserve(_grow_capacity(current_size + 1));
			return data[(++current_size - 1 + start_index) % capacity];
		}

		inline tfxU32       _grow_capacity(tfxU32 sz) const { tfxU32 new_capacity = capacity ? (capacity + capacity / 2) : 8; return new_capacity > sz ? new_capacity : sz; }
		inline void         reserve(tfxU32 new_capacity) {
			if (new_capacity <= capacity) return;
			T* new_data = (T*)tfxALLOCATE(name, new_data, (size_t)new_capacity * sizeof(T));
			if (data) {
				if (last_index() < start_index) {
					memcpy(new_data, data + start_index, (size_t)(capacity - start_index) * sizeof(T));
					memcpy(new_data + (capacity - start_index), data, (size_t)(start_index) * sizeof(T));
				}
				else {
					memcpy(new_data, data + start_index, (size_t)current_size * sizeof(T));
				}
				free(data);
			}
			data = new_data;
			capacity = new_capacity;
			start_index = 0;
		}

		inline void			bump() { if (current_size == 0) return; start_index++; start_index %= capacity; current_size--; }
		inline void			bump(unsigned int amount) { if (current_size == 0) return; if (amount > current_size) amount = current_size; start_index += amount; start_index %= capacity; current_size -= amount; last_bump = amount; }
		inline void			shrink(unsigned int amount) { if (amount > current_size) current_size = 0; else current_size -= amount; }
	};

#define tfxKilobyte(Value) ((Value)*1024LL)
#define tfxMegabyte(Value) (tfxKilobyte(Value)*1024LL)
#define tfxGigabyte(Value) (tfxMegabyte(Value)*1024LL)

	inline tfxU32 IsPowerOf2(tfxU32 v)
	{
		return ((v & ~(v - 1)) == v);
	}

	struct tfxMemoryBucket {
		void *data = NULL;
		void *end_ptr = NULL;
		tfxU32 next_block = tfxINVALID;
		tfxU32 unit_size;
		tfxU32 capacity = 0;
		tfxU32 current_size = 0;
		tfxU32 arena_index;

		inline void*           begin() { return data; }
		inline const void*     begin() const { return data; }
		inline void*           end() { return end_ptr; }
		inline const void*     end() const { return end_ptr; }
		inline void			   clear() { current_size = 0; end_ptr = data; }
		inline void			   reset() { current_size = 0; next_block = tfxINVALID; end_ptr = data; }
		inline size_t		   capacity_in_bytes() { return (size_t)unit_size * (size_t)capacity; }

	};

	template <typename T>
	inline T* BlockBegin(tfxMemoryBucket &block) {
		return (T*)block.data;
	}

	template <typename T>
	inline T* BlockEnd(tfxMemoryBucket &block) {
		return (T*)block.data + block.current_size;
	}

	template <typename T>
	inline T& BlockBack(tfxMemoryBucket &block) {
		return *((T*)block.data + block.current_size - 1);
	}

	template <typename T>
	inline T& BlockFront(tfxMemoryBucket &block) {
		return *((T*)block.data);
	}

	template <typename T>
	inline bool PushBack(tfxMemoryBucket &block, const T &v) {
		if (block.current_size != block.capacity) {
			*(T*)((T*)block.data + block.current_size++) = v;
			block.end_ptr = (T*)block.end_ptr + 1;
			return true;
		}
		return false;
	}

	template <typename T>
	inline tfxU32 BumpBlock(tfxMemoryBucket &block) {
		assert(block.current_size != block.capacity);	//Should not be bumping the block when there's no space
		block.current_size++;
		block.end_ptr = (T*)block.data + block.current_size;
		return block.current_size;
	}

	template <typename T>
	inline tfxU32 PopBlock(tfxMemoryBucket &block) {
		assert(block.current_size);		//Nothing to Pop
		block.current_size--;
		block.end_ptr = (T*)block.data + block.current_size;
		return block.current_size;
	}

	template <typename T>
	inline T &ValueAt(tfxMemoryBucket &block, tfxU32 index) {
		assert(index < block.current_size);		//Index was out of bounds
		return *((T*)block.data + index);
	}

	struct tfxMemoryArena {
		void *data;					//big allocation for splitting up into smaller blocks
		void *end_of_allocated;
		size_t memory_remaining;
		size_t total_memory;

		tfxMemoryArena() { data = end_of_allocated = NULL; memory_remaining = 0; }

		inline void FreeAll() {
			if (data) {
				memory_remaining = 0;
				free(data);
				data = NULL;
			}
		}

		inline void Reset() {
			memory_remaining = total_memory;
			end_of_allocated = data;
			memset((void*)data, 0, total_memory);
		}

	};

	struct tfxMemoryArenaManager {
		size_t arena_size;
		size_t size_diff_threshold;
		tfxvec<tfxMemoryArena> arenas;
		tfxvec<tfxMemoryBucket> blocks;
		tfxvec<tfxU32> free_blocks;

		inline void FreeAll() {
			for (auto &arena : arenas) {
				arena.FreeAll();
			}
			arenas.free_all();
			blocks.free_all();
			free_blocks.free_all();
		}

		inline void FreeBlock(tfxMemoryBucket *block) {
			tfxU32 index = blocks.index_from_ptr(block);
			free_blocks.push_back(index);
		}

		inline tfxU32 FreeBlocks(tfxU32 block) {
			if (block >= blocks.current_size) return 0;
			free_blocks.push_back(block);
			tfxU32 freed_count = 1;
			while (blocks[block].next_block != tfxINVALID) {
				tfxU32 prev_block = block;
				block = blocks[block].next_block;
				free_blocks.push_back(block);
				freed_count++;
				blocks[prev_block].reset();
			}
			blocks[block].reset();
			return freed_count;
		}

		inline void *End(tfxU32 block) {
			while (blocks[block].next_block != tfxINVALID) {
				block = blocks[block].next_block;
			}
			return blocks[block].end();
		}

		inline tfxMemoryBucket &LastBlock(tfxU32 starting_block) {
			assert(starting_block != tfxINVALID);
			tfxU32 found_block = starting_block;
			while (blocks[found_block].next_block != tfxINVALID) {
				found_block = blocks[found_block].next_block;
			}
			return blocks[found_block];
		}

		inline tfxMemoryBucket &FirstBlockWithSpace(tfxU32 starting_block) {
			assert(starting_block != tfxINVALID);
			tfxU32 found_block = starting_block;
			while (blocks[found_block].next_block != tfxINVALID && blocks[found_block].current_size == blocks[found_block].capacity) {
				found_block = blocks[found_block].next_block;
			}
			return blocks[found_block];
		}

		inline tfxU32 FirstBlockIndexWithSpace(tfxU32 starting_block) {
			assert(starting_block != tfxINVALID);
			tfxU32 found_block = starting_block;
			while (blocks[found_block].next_block != tfxINVALID && blocks[found_block].current_size == blocks[found_block].capacity) {
				found_block = blocks[found_block].next_block;
			}
			return found_block;
		}

		inline tfxU32 FirstEmptyBlockIndex(tfxU32 starting_block) {
			assert(starting_block != tfxINVALID);
			tfxU32 found_block = starting_block;
			while (blocks[found_block].next_block != tfxINVALID) {
				if (blocks[found_block].current_size == 0) {
					return found_block;
				}
				found_block = blocks[found_block].next_block;
			}
			if (blocks[found_block].current_size == 0) {
				return found_block;
			}
			return tfxINVALID;
		}

		inline void CutOffBlock(tfxU32 starting_block, tfxU32 block) {
			assert(block != tfxINVALID);
			tfxU32 found_block = starting_block;
			while (blocks[found_block].next_block != tfxINVALID) {
				if (blocks[found_block].next_block == block) {
					blocks[found_block].next_block = tfxINVALID;
					return;
				}
				found_block = blocks[found_block].next_block;
			}
		}

		inline tfxMemoryBucket &FirstBlockWithSpace(tfxU32 current_size, tfxU32 block) {
			assert(block != tfxINVALID);
			tfxU32 block_index = current_size / blocks[block].capacity;
			for (int i = 0; i != block_index; ++i) {
				block = blocks[block].next_block;
			}
			return blocks[block];
		}

		inline bool FirstArenaWithEnoughSpace(size_t required_space_in_bytes, tfxMemoryArena **found_arena, tfxU32 &arena_index) {
			arena_index = 0;
			for (auto &arena : arenas) {
				if (arena.memory_remaining >= required_space_in_bytes) {
					*found_arena = &arena;
					return true;
				}
				++arena_index;
			}
			*found_arena = NULL;
			return false;
		}

		inline tfxMemoryBucket &BlockByIndex(tfxU32 index, tfxU32 starting_block) {
			assert(starting_block != tfxINVALID);
			tfxU32 found_block = starting_block;
			tfxU32 index_count = 0;
			while (blocks[found_block].next_block != tfxINVALID && index_count++ != index) {
				found_block = blocks[found_block].next_block;
			}
			return blocks[found_block];
		}

		inline tfxU32 BlockIndexByIndex(tfxU32 index, tfxU32 starting_block) {
			assert(starting_block != tfxINVALID);
			tfxU32 found_block = starting_block;
			tfxU32 index_count = 0;
			while (blocks[found_block].next_block != tfxINVALID && index_count++ != index) {
				found_block = blocks[found_block].next_block;
			}
			return found_block;
		}

		inline void ClearBlocks(tfxU32 block) {
			if (block == tfxINVALID) return;
			blocks[block].clear();
			while (blocks[block].next_block != tfxINVALID) {
				block = blocks[block].next_block;
				blocks[block].clear();
			}
		}

		inline tfxU32 BlockCount(tfxU32 block) {
			if (block == tfxINVALID) return 0;
			tfxU32 count = 1;
			while (blocks[block].next_block != tfxINVALID) {
				block = blocks[block].next_block;
				count++;
			}
			return count;
		}

		inline void CopyBlockToBlock(tfxU32 from, tfxU32 to) {
			assert(blocks[from].capacity && blocks[from].capacity <= blocks[to].capacity);		//must have valid capacities
			memcpy(blocks[to].data, blocks[from].data, blocks[from].capacity * blocks[from].unit_size);
			auto &src = blocks[from];
			auto &dst = blocks[to];
			dst.current_size = src.current_size;
			dst.end_ptr = (char*)dst.data + (dst.unit_size * dst.current_size);
		}

		inline void CopyBlockToBlock(tfxMemoryBucket *from, tfxMemoryBucket *to) {
			assert(from->capacity && from->capacity <= to->capacity);		//must have valid capacities
			memcpy(to->data, from->data, from->capacity * from->unit_size);
			to->current_size = from->current_size;
			to->end_ptr = (char*)to->data + (to->unit_size * to->current_size);
		}

		inline tfxMemoryArena *AddArena() {
			tfxMemoryArena arena;
			arena.total_memory = arena_size;
			arena.memory_remaining = arena_size;
			arena.data = malloc(arena_size);
			arena.end_of_allocated = arena.data;
			memset(arena.data, 0, arena_size);
			arenas.push_back(arena);
			return &arenas.back();
		}

		inline size_t TotalMemoryCapacity() {
			size_t size = 0;
			for (auto &arena : arenas) {
				size += arena.total_memory;
			}
			return size;
		}

		inline size_t TotalMemoryInUse() {
			size_t size = 0;
			for (auto &arena : arenas) {
				size += arena.total_memory - arena.memory_remaining;
			}
			return size;
		}

		void CheckForOverlappingBlocks() {
			tfxMemoryBucket *last_block = NULL;
			for (auto &block : blocks) {
				if (last_block) {
					void *last_end = (char*)last_block->data + last_block->capacity_in_bytes();
					const ptrdiff_t offset = (char*)last_end - (char*)block.data;
					if (offset < 0) {
						printf("Overlapping memory blocks found. %zi", offset);
					}
				}
				last_block = &block;
			}
		}
	};

	inline tfxU64 NearestMultiple(tfxU64 numToRound, tfxU64 multiple)
	{
		assert(multiple);
		return ((numToRound + multiple - 1) / multiple) * multiple;
	}

	static inline tfxMemoryArenaManager CreateArenaManager(size_t size_of_each_arena, tfxU32 size_diff_threshold = 8) {
		tfxMemoryArenaManager manager;
		manager.arena_size = size_of_each_arena;
		manager.size_diff_threshold = size_diff_threshold;
		return manager;
	}

	static inline tfxMemoryArena CreateMemoryArena(size_t size_in_bytes, tfxU32 size_diff_threshold = 8) {
		assert(size_in_bytes > 1024 * 1024);	//minimum 1mb allocation
		tfxMemoryArena allocator;
		void* new_data = malloc(size_in_bytes);
		allocator.data = new_data;
		allocator.end_of_allocated = allocator.data;
		allocator.memory_remaining = size_in_bytes;
		allocator.total_memory = size_in_bytes;
		memset((void*)allocator.data, 0, size_in_bytes);
		return allocator;
	}

	template <typename T>
	static bool Allocate(tfxMemoryArenaManager &allocator, tfxU32 block_size, tfxU32 &block) {
		tfxU32 size_in_bytes = block_size * sizeof(T);
		if (size_in_bytes == 0) return false;
		if (allocator.free_blocks.current_size > 0) {
			size_t size_diff = tfxMAX_UINT;
			tfxU32 best_fit = tfxINVALID;
			tfxU32 found_index = 0;
			tfxU32 i = 0;
			for (auto free_block : allocator.free_blocks) {
				if (allocator.blocks[free_block].capacity_in_bytes() >= size_in_bytes && allocator.blocks[free_block].capacity_in_bytes() - size_in_bytes < size_diff) {
					size_diff = allocator.blocks[free_block].capacity_in_bytes() - size_in_bytes;
					best_fit = free_block;
					found_index = i;
					if (size_diff == 0)
						break;
				}
				++i;
			}
			if (best_fit != tfxINVALID && size_diff != tfxMAX_UINT && (size_diff <= allocator.size_diff_threshold)) {
				allocator.free_blocks[found_index] = allocator.free_blocks.pop_back();
				allocator.FreeBlocks(block);
				block = best_fit;
				return true;
			}
		}

		//If the block exists and it's at the end of the area then don't add it to the free list, instead just remove the block and add the new on the end
		if (block != tfxINVALID) {
			tfxU32 old_block_size = allocator.blocks[block].capacity;
			void *end_ptr = (T*)allocator.blocks[block].data + old_block_size;
			if (block == allocator.blocks.current_size - 1 && allocator.arenas[allocator.blocks[block].arena_index].end_of_allocated == end_ptr) {
				allocator.arenas[allocator.blocks[block].arena_index].end_of_allocated = (T*)allocator.arenas[allocator.blocks[block].arena_index].end_of_allocated - old_block_size;
				allocator.arenas[allocator.blocks[block].arena_index].memory_remaining += old_block_size * sizeof(T);
				allocator.blocks.pop();
			}
			else {
				allocator.FreeBlocks(block);
			}
		}

		tfxMemoryArena *arena;
		tfxU32 arena_index = 0;
		if (!allocator.FirstArenaWithEnoughSpace(size_in_bytes, &arena, arena_index)) {
			arena = allocator.AddArena();
			arena_index = allocator.arenas.current_size - 1;
		}

		tfxMemoryBucket new_block;
		new_block.unit_size = sizeof(T);
		new_block.data = arena->end_of_allocated;
		new_block.end_ptr = arena->end_of_allocated;
		new_block.capacity = block_size;
		new_block.arena_index = arena_index;
		arena->end_of_allocated = (T*)arena->end_of_allocated + block_size;
		arena->memory_remaining -= size_in_bytes;
		allocator.blocks.push_back(new_block);
		block = allocator.blocks.current_size - 1;
		return true;

	}

	template <typename T>
	static bool AllocateBucket(tfxMemoryArenaManager &allocator, tfxU32 bucket_size, tfxU32 &block) {
		tfxU32 size_in_bytes = bucket_size * sizeof(T);
		assert(bucket_size > 1 && size_in_bytes > 0);		//bucket size must be greater than 1
		if (allocator.free_blocks.current_size > 0) {
			tfxU32 i = 0;
			for (auto free_block : allocator.free_blocks) {
				if (allocator.blocks[free_block].capacity * allocator.blocks[free_block].unit_size == size_in_bytes) {
					allocator.free_blocks[i] = allocator.free_blocks.pop_back();
					if (block == tfxINVALID)
						block = free_block;
					else
						allocator.LastBlock(block).next_block = free_block;
					return true;
				}
				++i;
			}
		}

		tfxMemoryArena *arena;
		tfxU32 arena_index = 0;
		if (!allocator.FirstArenaWithEnoughSpace(size_in_bytes, &arena, arena_index)) {
			arena = allocator.AddArena();
			arena_index = allocator.arenas.current_size - 1;
		}

		tfxMemoryBucket new_block;
		new_block.unit_size = sizeof(T);
		new_block.data = arena->end_of_allocated;
		new_block.end_ptr = arena->end_of_allocated;
		new_block.capacity = bucket_size;
		new_block.arena_index = arena_index;
		arena->end_of_allocated = (T*)arena->end_of_allocated + bucket_size;
		arena->memory_remaining -= size_in_bytes;
		allocator.blocks.push_back(new_block);
		if (block == tfxINVALID)
			block = allocator.blocks.current_size - 1;
		else
			allocator.LastBlock(block).next_block = allocator.blocks.current_size - 1;

		return true;
	}

	template <typename T>
	inline T &FindValueByIndex(tfxMemoryArenaManager &allocator, tfxU32 i, tfxU32 block) {
		while (i >= allocator.blocks[block].current_size && allocator.blocks[block].next_block != tfxINVALID) {
			i -= allocator.blocks[block].current_size;
			block = allocator.blocks[block].next_block;
		}
		return *((T*)allocator.blocks[block].data + i);
	}

	template <typename T>
	inline T &FindValueByIndex(tfxMemoryBucket &range, tfxU32 i) {
		return *((T*)range.data + i);
	}

#define tfxBucket(type, index, bucket_ptr) (type*)bucket_ptr + index;

	//No Destructor, so use free_all before it goes out of scope!
	template <typename T>
	struct tfxArray {
		tfxMemoryArenaManager *allocator;			//Pointer to the allocator that manages the memory and blocks of that memory
		T *block;									//Pointer to the data storing the array. This will be somewhere in a tfxMemoryArena
		tfxU32 block_index;							//This index of the block of memory referenced in allocator->blocks
		tfxU32 capacity;							//The total capacity of the array in units of T

		tfxArray() : allocator(NULL) { block = NULL; capacity = 0; block_index = tfxINVALID; }
		tfxArray(tfxMemoryArenaManager *allocator_init, tfxU32 size) : allocator(allocator_init) { block = NULL; capacity = 0; reserve(size); }

		inline tfxU32		size() { return capacity; }
		inline const tfxU32	size() const { return capacity; }
		inline T&           operator[](tfxU32 i) {
			assert(i < capacity);		//Index is out of bounds
			return block[i];
		}
		inline const T&     operator[](tfxU32 i) const {
			assert(i < capacity);		//Index is out of bounds
			return block[i];
		}
		inline tfxArray<T>&		operator=(const tfxArray<T>& src) {
			if (!allocator)
				allocator = src.allocator;
			if (!src.block) return *this;
			assert(resize(src.capacity));
			allocator->CopyBlockToBlock(src.block_index, block_index);
			return *this;
		}

		inline void         free() { if (block != NULL) { capacity = capacity = 0; allocator->FreeBlocks(block_index); block = NULL; block_index = tfxINVALID; } }
		inline T*           begin() { return block; }
		inline const T*     begin() const { return block; }
		inline T*           end() { return block + capacity; }
		inline const T*     end() const { return block + capacity; }
		inline bool			reserve(tfxU32 size) {
			assert(allocator);		//Must assign an allocator before doing anything with a tfxBucketArray. Capacity must equal 0
			assert(capacity == 0);	//Capacity must equal 0 before reserving an array
			assert(size * sizeof(T) < allocator->arena_size);	//The size of an array must fit into an arena size
			if (capacity == 0) {
				if (!Allocate<T>(*allocator, size, block_index)) {
					return false;
				}
				capacity = size;
				allocator->blocks[block_index].current_size = capacity;
				allocator->blocks[block_index].end_ptr = (T*)allocator->blocks[block_index].end_ptr + capacity;
				block = (T*)allocator->blocks[block_index].data;
			}
			return true;
		}
		inline bool			resize(tfxU32 size, bool keep_contents = false) {
			assert(allocator);		//Must assign an allocator before doing anything with a tfxBucketArray. Capacity must equal 0
			if (size == capacity) return true;
			tfxU32 current_block = block_index;
			if (!Allocate<T>(*allocator, size, block_index)) {
				return false;
			}
			if (keep_contents)
				allocator->CopyBlockToBlock(current_block, block_index);
			capacity = size;
			block = (T*)allocator->blocks[block_index].data;
			return true;
		}

	};

	template <typename T>
	struct tfxBucketArray {
		tfxMemoryArenaManager *allocator;	//Pointer to the arena manager that handles the memory and blocks of that memory
		tfxU32 block;						//The first block in the allocator
		tfxU32 current_bucket;				//the current bucket for iterating all blocks
		tfxU32 current_size;				//Current size of the bucket array. This will be the total of all buckets if there are more then one
		tfxU32 capacity;					//The total capacity of the bucket array
		tfxU32 size_of_each_bucket;			//The size of each bucket

		tfxBucketArray() { allocator = NULL; current_bucket = block = tfxINVALID; size_of_each_bucket = 64; current_size = capacity = size_of_each_bucket = 0; }
		tfxBucketArray(tfxMemoryArenaManager *allocator_init) : allocator(allocator_init) { size_of_each_bucket = 64; current_size = capacity = 0; current_bucket = block = tfxINVALID; }
		tfxBucketArray(tfxMemoryArenaManager *allocator_init, tfxU32 bucket_size) { assert(bucket_size > 1); size_of_each_bucket = bucket_size; allocator = allocator_init; current_size = capacity = 0; current_bucket = block = tfxINVALID; }

		inline bool			empty() { return current_size == 0; }
		inline tfxU32		size() { return current_size; }
		inline const tfxU32	size() const { return current_size; }
		inline T&           operator[](tfxU32 i) {
			assert(i < current_size);		//Index is out of bounds
			return FindValueByIndex<T>(*allocator, i, block);
		}
		inline const T&     operator[](tfxU32 i) const {
			assert(i < current_size);		//Index is out of bounds
			return FindValueByIndex<T>(*allocator, i, block);
		}
		inline tfxBucketArray<T>&		operator=(const tfxBucketArray<T>& src) {
			if (!allocator)
				allocator = src.allocator;
			size_of_each_bucket = src.size_of_each_bucket;
			if (src.capacity == 0) {
				return *this;
			}
			assert(reserve(src.capacity / src.size_of_each_bucket));
			current_bucket = block;
			tfxU32 src_block = src.block;
			while (current_bucket != tfxINVALID && src_block != tfxINVALID) {
				allocator->CopyBlockToBlock(src_block, current_bucket);
				current_bucket = allocator->blocks[current_bucket].next_block;
				src_block = src.allocator->blocks[src_block].next_block;
			}
			current_bucket = block;
			current_size = src.current_size;
			return *this;
		}

		inline void         free_all() { if (block != tfxINVALID) { current_size = capacity = 0; current_bucket = tfxINVALID; allocator->FreeBlocks(block); block = tfxINVALID; } }
		inline T*           begin() { return current_bucket != tfxINVALID ? (T*)allocator->blocks[current_bucket].data : NULL; }
		inline const T*     begin() const { return current_bucket != tfxINVALID ? (T*)allocator->blocks[current_bucket].data : NULL; }
		inline T*           end() { return current_bucket != tfxINVALID ? (T*)allocator->blocks[current_bucket].end_ptr : NULL; }
		inline const T*     end() const { return current_bucket != tfxINVALID ? (T*)allocator->blocks[current_bucket].end_ptr : NULL; }
		inline T*			bucket_end() { return (T*)allocator->FirstBlockWithSpace(block).end(); }
		inline T&           front() { assert(current_size > 0); return *(T*)allocator->blocks[block].data; }
		inline const T&     front() const { assert(current_size > 0); return *(T*)allocator->blocks[block].data; }
		inline T&           back() { assert(current_size > 0); return BlockBack<T>(allocator->FirstBlockWithSpace(block)); }
		inline const T&     back() const { assert(current_size > 0); return BlockBack<T>(allocator->FirstBlockWithSpace(block)); }
		inline void         clear() {
			current_size = 0;
			allocator->ClearBlocks(block);
		}
		inline tfxU32 bump() {
			assert(current_size != capacity);
			current_size++;
			return (tfxU32)BumpBlock<T>(allocator->FirstBlockWithSpace(block)) - 1;
		}
		inline tfxU32 bump(tfxU32 block_to_bump) {
			//You must ensure that this block belongs to this bucket array
			assert(current_size != capacity);
			current_size++;
			allocator->blocks[block_to_bump].current_size++;
			allocator->blocks[block_to_bump].end_ptr = (T*)allocator->blocks[block_to_bump].data + allocator->blocks[block_to_bump].current_size;
			return current_size - 1;
		}
		inline bool			reserve(int number_of_buckets) {
			assert(allocator);										//Must assign and allocator before doing anything with a tfxBucketArray
			assert(size_of_each_bucket * number_of_buckets > 1);	//Buckets must be greater than 0
			int block_count = (int)allocator->BlockCount(block);
			if (block_count > 0) {
				number_of_buckets = (int)allocator->BlockCount(block) - (int)number_of_buckets;
			}
			for (int i = 0; i < number_of_buckets; ++i) {
				assert(AllocateBucket<T>(*allocator, size_of_each_bucket, block));		//Out of memory!
				capacity += size_of_each_bucket;
			}
			return true;
		}
		inline T&	        push_back(const T& v) {
			assert(allocator);	//Must assign an allocator before doing anything with a tfxBucketArray
			if (current_size == capacity) {
				assert(AllocateBucket<T>(*allocator, size_of_each_bucket, block));		//Out of memory!
				capacity += size_of_each_bucket;
				ResetIteratorIndex();
			}
			tfxMemoryBucket &last_block = allocator->FirstBlockWithSpace(block);
			assert(last_block.current_size < last_block.capacity);
			*(T*)((T*)last_block.data + last_block.current_size++) = v;
			last_block.end_ptr = (T*)last_block.end_ptr + 1;
			current_size++;
			return BlockBack<T>(last_block);
		}
		inline T*	insert(tfxU32 insert_index, const T &v) {
			assert(insert_index < current_size);
			if (current_size == capacity) {
				assert(AllocateBucket<T>(*allocator, size_of_each_bucket, block));		//Out of memory!
				capacity += size_of_each_bucket;
				ResetIteratorIndex();
			}
			tfxU32 bucket_index = insert_index / size_of_each_bucket;
			insert_index -= bucket_index * size_of_each_bucket;
			tfxU32 index_block = allocator->BlockIndexByIndex(bucket_index, block);
			T value_to_insert = v;
			T* return_value = NULL;
			bool initial_inserted = false;
			do {
				T* insert_point = &ValueAt<T>(allocator->blocks[index_block], insert_index);
				size_t move_size = size_of_each_bucket - insert_index - (allocator->blocks[index_block].current_size == allocator->blocks[index_block].capacity ? 1 : 0);
				T value_at_back = BlockBack<T>(allocator->blocks[index_block]);
				if (move_size > 0) {
					memmove(insert_point + 1, insert_point, move_size * sizeof(T));
				}
				*insert_point = value_to_insert;
				if (!initial_inserted) {
					initial_inserted = true;
					return_value = insert_point;
				}
				if (allocator->blocks[index_block].current_size < allocator->blocks[index_block].capacity)
					BumpBlock<T>(allocator->blocks[index_block]);
				index_block = allocator->blocks[index_block].next_block;
				if (index_block != tfxINVALID && allocator->blocks[index_block].current_size == 0) {
					PushBack<T>(allocator->blocks[index_block], value_at_back);
					break;
				}
				else {
					value_to_insert = value_at_back;
				}
				insert_index = 0;
			} while (index_block != tfxINVALID);
			current_size++;
			return return_value;
		}
		inline T*	insert(const T* it, const T &v) {
			tfxU32 index;
			T* found_value = find(*it, index);
			return insert(index, v);
		}
		inline bool	erase(tfxU32 erase_index) {
			assert(erase_index < current_size);
			tfxU32 bucket_index = erase_index / size_of_each_bucket;
			erase_index -= bucket_index * size_of_each_bucket;
			tfxU32 index_block = allocator->BlockIndexByIndex(bucket_index, block);
			tfxU32 last_block = tfxINVALID;
			do {
				T* erase_point = &ValueAt<T>(allocator->blocks[index_block], erase_index);
				size_t move_size = size_of_each_bucket - erase_index - 1;
				T value_at_back = BlockBack<T>(allocator->blocks[index_block]);
				if (move_size > 0) {
					memmove(erase_point, erase_point + 1, move_size * sizeof(T));
				}
				T &end_of_block = BlockBack<T>(allocator->blocks[index_block]);
				if (allocator->blocks[index_block].next_block != tfxINVALID) {
					end_of_block = BlockFront<T>(allocator->blocks[allocator->blocks[index_block].next_block]);
				}
				last_block = index_block;
				index_block = allocator->blocks[index_block].next_block;
				erase_index = 0;
			} while (index_block != tfxINVALID && allocator->blocks[index_block].current_size > 0);
			PopBlock<T>(allocator->blocks[last_block]);
			current_size--;
			return true;
		}
		inline bool	erase(const T* it) {
			tfxU32 index;
			T* found_value = find(*it, index);
			return erase(index);
		}
		inline T* find(const T& v) {
			assert(block != tfxINVALID);		//bucket array must be initialised
			tfxU32 current_block = block;
			T *_data = (T*)allocator->blocks[block].data;
			while (current_block != tfxINVALID && allocator->blocks[current_block].current_size > 0) {
				_data = (T*)allocator->blocks[current_block].data;
				while (_data < allocator->blocks[current_block].end_ptr) {
					if (*_data == v)
						return _data;
					++_data;
				}
				current_block = allocator->blocks[current_block].next_block;
			}
			return _data;
		}
		inline T* find(const T& v, tfxU32 &index) {
			assert(block != tfxINVALID);
			index = 0;
			tfxU32 current_block = block;
			T *_data = (T*)allocator->blocks[block].data;
			while (current_block != tfxINVALID && allocator->blocks[current_block].current_size > 0) {
				_data = (T*)allocator->blocks[current_block].data;
				while (_data < allocator->blocks[current_block].end_ptr) {
					if (*_data == v)
						return _data;
					++index;
					++_data;
				}
				current_block = allocator->blocks[current_block].next_block;
			}
			return _data;
		}
		inline void ResetIteratorIndex() {
			current_bucket = block;
		}
		inline bool EndOfBuckets() {
			current_bucket = allocator->blocks[current_bucket].next_block;
			if (current_bucket == tfxINVALID) {
				current_bucket = block;
				return true;
			}
			return false;
		}
		inline tfxMemoryBucket &NextBlock() {
			assert(allocator);	//Must assign and allocator before doing anything with a tfxBucketArray
			return allocator->blocks[current_bucket];
		}
		inline tfxMemoryBucket &GetFirstBucket() {
			assert(allocator);	//Must assign an allocator before doing anything with a tfxBucketArray
			return allocator->blocks[block];
		}
		inline void TrimBuckets() {
			tfxU32 first_empty_block = allocator->FirstEmptyBlockIndex(block);
			if (first_empty_block != tfxINVALID) {
				tfxU32 freed_blocks = allocator->FreeBlocks(first_empty_block);
				capacity -= freed_blocks * size_of_each_bucket;
				allocator->CutOffBlock(block, first_empty_block);
			}
		}
	};

	template <typename T>
	static inline tfxBucketArray<T> CreateBucketArray(tfxMemoryArena *allocator, tfxU32 bucket_size) {
		tfxBucketArray bucket_array(allocator, bucket_size);
		return bucket_array;
	}

	//A char buffer you can use to load a file into and read from
	//Has no deconstructor so make sure you call FreeAll() when done
	//This is meant for limited usage in timeline fx only and not recommended for use outside!
	struct tfxstream {
		tfxU64 size = 0;
		tfxU64 position = 0;
		char* data = NULL;

		inline tfxstream() { size = position = 0; data = NULL; }
		inline tfxstream(tfxU64 qty) { size = position = 0; data = NULL; Resize(qty); }
		inline tfxstream(const tfxstream &src) { size = 0; data = NULL; Resize(src.size); memcpy(data, src.data, (tfxU64)size * sizeof(char)); }

		inline bool Read(char* dst, tfxU64 count) {
			if (count + position <= size) { 
				memcpy(dst, data + position, count); 
				position += count; 
				return true;
			}
			else { 
				assert(false); //Trying to read beyond the data boundary
			}
			return false;
		}
		inline tfxStr512 ReadLine();
		inline bool Write(void *src, tfxU64 count) {
			if (count + position <= size) {
				memcpy(data + position, src, count);
				position += count;
				return true;
			}
			else {
				assert(false); //Trying to write beyond the data boundary
			}
			return false;
		}
		inline bool EoF() { return position >= size; }
		inline void Seek(tfxU64 offset) {
			if (offset < size)
				position = offset;
			else
				position = size;
		}

		inline bool			Empty() { return size == 0; }
		inline tfxU64		Size() { return size; }
		inline const tfxU64	Size() const { return size; }

		inline void			FreeAll() { if (data) { size = size = 0; free(data); data = NULL; } }
		inline void         Clear() { if (data) { size = 0; } }

		inline void         Resize(tfxU64 new_capacity) { if (new_capacity <= size) return; char* new_data = (char*)malloc((tfxU64)new_capacity * sizeof(char)); if (data) { memcpy(new_data, data, (tfxU64)size * sizeof(char)); free(data); } data = new_data; size = new_capacity; position = 0; }
		inline void			NullTerminate() { *(data + size) = NULL; }

	};

	//Just the very basic vector types that we need
	struct tfxVec2 {
		float x, y;

		tfxVec2() { x = y = 0.f; }
		tfxVec2(float _x, float _y) : x(_x), y(_y) {}

		inline void operator=(float v) { x = v; y = v; }

		inline tfxVec2 operator+(tfxVec2 v) { return tfxVec2(x + v.x, y + v.y); }
		inline void operator+=(tfxVec2 v) { x += v.x; y += v.y; }
		inline tfxVec2 operator-(tfxVec2 v) { return tfxVec2(x - v.x, y - v.y); }
		inline tfxVec2 operator-() { return tfxVec2(-x, -y); }
		inline void operator-=(tfxVec2 v) { x -= v.x; y -= v.y; }
		inline tfxVec2 operator*(tfxVec2 v) { return tfxVec2(x * v.x, y * v.y); }
		inline void operator*=(tfxVec2 v) { x *= v.x; y *= v.y; }
		inline tfxVec2 operator/(tfxVec2 v) { return tfxVec2(x / v.x, y / v.y); }
		inline void operator/=(tfxVec2 v) { x /= v.x; y /= v.y; }
		inline tfxVec2 operator+(float v) { return tfxVec2(x + v, y + v); }
		inline tfxVec2 operator-(float v) { return tfxVec2(x - v, y - v); }
		inline tfxVec2 operator*(float v) { return tfxVec2(x * v, y * v); }
		inline void operator*=(float v) { x *= v; y *= v; }
		inline bool operator>(tfxVec2 &v) { return x + y > v.x + v.y; }
		inline tfxVec2 operator/(float v) { return tfxVec2(x / v, y / v); }

		inline tfxVec2 operator+(float v) const { return tfxVec2(x + v, y + v); }
		inline tfxVec2 operator-(float v) const { return tfxVec2(x - v, y - v); }
		inline tfxVec2 operator*(float v) const { return tfxVec2(x * v, y * v); }
		inline tfxVec2 operator/(float v) const { return tfxVec2(x / v, y / v); }

		inline tfxVec2 operator+(const tfxVec2 &v) const { return tfxVec2(x + v.x, y + v.y); }
		inline tfxVec2 operator-(const tfxVec2 &v) const { return tfxVec2(x - v.x, y - v.y); }
		inline tfxVec2 operator*(const tfxVec2 &v) const { return tfxVec2(x * v.x, y * v.y); }
		inline tfxVec2 operator/(const tfxVec2 &v) const { return tfxVec2(x / v.x, y / v.y); }

		inline float Squared() { return x * x + y * y; }
		inline bool IsNill() { return !x && !y; }
	};
	inline tfxVec2 operator*(float ls, tfxVec2 rs) { return tfxVec2(rs.x * ls, rs.y * ls); }

	struct tfxVec3 {
		union {
			struct { float x, y, z; };
			struct { float pitch, yaw, roll; };
		};

		tfxVec3() { x = y = z = 0.f; }
		tfxVec3(float v) : x(v), y(v), z(v) {}
		tfxVec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
		inline void operator=(const tfxVec2 &v) { x = v.x; y = v.y; }

		inline tfxVec2 xy() const { return tfxVec2(x, y); }

		inline bool operator==(const tfxVec3 &v) const { return x == v.x && y == v.y && z == v.z; }

		inline tfxVec3 operator+(const tfxVec3 &v) const { return tfxVec3(x + v.x, y + v.y, z + v.z); }
		inline tfxVec3 operator-(const tfxVec3 &v) const { return tfxVec3(x - v.x, y - v.y, z - v.z); }
		inline tfxVec3 operator*(const tfxVec3 &v) const { return tfxVec3(x * v.x, y * v.y, z * v.z); }
		inline tfxVec3 operator/(const tfxVec3 &v) const { return tfxVec3(x / v.x, y / v.y, z / v.z); }

		inline tfxVec3 operator-() const { return tfxVec3(-x, -y, -z); }

		inline void operator-=(const tfxVec3 &v) { x -= v.x; y -= v.y; z -= v.z; }
		inline void operator+=(const tfxVec3 &v) { x += v.x; y += v.y; z += v.z; }
		inline void operator*=(const tfxVec3 &v) { x *= v.x; y *= v.y; z *= v.z; }
		inline void operator/=(const tfxVec3 &v) { x /= v.x; y /= v.y; z /= v.z; }

		inline void operator-=(const tfxVec2 &v) { x -= v.x; y -= v.y; }
		inline void operator+=(const tfxVec2 &v) { x += v.x; y += v.y; }
		inline void operator*=(const tfxVec2 &v) { x *= v.x; y *= v.y; }
		inline void operator/=(const tfxVec2 &v) { x /= v.x; y /= v.y; }

		inline tfxVec3 operator+(float v) const { return tfxVec3(x + v, y + v, z + v); }
		inline tfxVec3 operator-(float v) const { return tfxVec3(x - v, y - v, z - v); }
		inline tfxVec3 operator*(float v) const { return tfxVec3(x * v, y * v, z * v); }
		inline tfxVec3 operator/(float v) const { return tfxVec3(x / v, y / v, z / v); }

		inline void operator*=(float v) { x *= v; y *= v; z *= v; }
		inline void operator/=(float v) { x /= v; y /= v; z /= v; }
		inline void operator+=(float v) { x += v; y += v; z += v; }
		inline void operator-=(float v) { x -= v; y -= v; z -= v; }

		inline float Squared() { return x * x + y * y + z * z; }
		inline bool IsNill() { return !x && !y && !z; }
	};

	struct tfxVec4 {
		float x, y, z, w;

		tfxVec4() { x = y = z = w = 0.f; }
		tfxVec4(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}
		tfxVec4(tfxVec2 vec1, tfxVec2 vec2) : x(vec1.x), y(vec1.y), z(vec2.x), w(vec2.y) {}
		tfxVec4(tfxVec3 vec) : x(vec.x), y(vec.y), z(vec.z), w(0.f) {}

		inline tfxVec2 xy() { return tfxVec2(x, y); }
		inline tfxVec2 zw() { return tfxVec2(z, w); }
		inline tfxVec3 xyz() { return tfxVec3(x, y, z); }

		inline tfxVec2 xy() const{ return tfxVec2(x, y); }
		inline tfxVec2 zw() const { return tfxVec2(z, w); }
		inline tfxVec3 xyz() const { return tfxVec3(x, y, z); }

		inline void operator=(const tfxVec2 &v) { x = v.x; y = v.y; }
		inline void operator=(const tfxVec3 &v) { x = v.x; y = v.y; z = v.z; }

		inline tfxVec4 operator+(const tfxVec4 &v) { return tfxVec4(x + v.x, y + v.y, z + v.z, w + v.w); }
		inline tfxVec4 operator-(const tfxVec4 &v) { return tfxVec4(x - v.x, y - v.y, z - v.z, w - v.w); }
		inline tfxVec4 operator-() { return tfxVec4(-x, -y, -z, -w); }
		inline tfxVec4 operator*(const tfxVec4 &v) { return tfxVec4(x * v.x, y * v.y, z * v.z, w * v.w); }
		inline tfxVec4 operator/(const tfxVec4 &v) { return tfxVec4(x / v.x, y / v.y, z / v.z, w / v.w); }

		inline void operator-=(const tfxVec4 &v) { x -= v.x; y -= v.y; z -= v.z; w -= v.w; }
		inline void operator+=(const tfxVec4 &v) { x += v.x; y += v.y; z += v.z; w += v.w; }
		inline void operator*=(const tfxVec4 &v) { x *= v.x; y *= v.y; z *= v.z; w *= v.w; }
		inline void operator/=(const tfxVec4 &v) { x /= v.x; y /= v.y; z /= v.z; w /= v.w; }

		inline void operator-=(const tfxVec3 &v) { x -= v.x; y -= v.y; z -= v.z; }
		inline void operator+=(const tfxVec3 &v) { x += v.x; y += v.y; z += v.z; }
		inline void operator*=(const tfxVec3 &v) { x *= v.x; y *= v.y; z *= v.z; }
		inline void operator/=(const tfxVec3 &v) { x /= v.x; y /= v.y; z /= v.z; }

		inline void operator-=(const tfxVec2 &v) { x -= v.x; y -= v.y; }
		inline void operator+=(const tfxVec2 &v) { x += v.x; y += v.y; }
		inline void operator*=(const tfxVec2 &v) { x *= v.x; y *= v.y; }
		inline void operator/=(const tfxVec2 &v) { x /= v.x; y /= v.y; }

		inline tfxVec4 operator+(float v) const { return tfxVec4(x + v, y + v, z + v, w + v); }
		inline tfxVec4 operator-(float v) const { return tfxVec4(x - v, y - v, z - v, w - v); }
		inline tfxVec4 operator*(float v) const { return tfxVec4(x * v, y * v, z * v, w * v); }
		inline tfxVec4 operator/(float v) const { return tfxVec4(x / v, y / v, z / v, w / v); }

		inline void operator*=(float v) { x *= v; y *= v; z *= v; w *= v; }
		inline void operator/=(float v) { x /= v; y /= v; z /= v; w /= v; }
		inline void operator+=(float v) { x += v; y += v; z += v; w += v; }
		inline void operator-=(float v) { x -= v; y -= v; z -= v; w -= v; }
	};

	static inline void ScaleVec4xyz(tfxVec4 &v, float scalar) {
		v.x *= scalar;
		v.y *= scalar;
		v.z *= scalar;
	}

	struct tfxRGBA8 {
		unsigned char r, g, b, a;

		tfxRGBA8() { r = g = b = a = 0; }
		tfxRGBA8(unsigned char _r, unsigned char _g, unsigned char _b, unsigned char _a) : r(_r), g(_g), b(_b), a(_a) { }
	};

	struct tfxRGB {
		float r, g, b;
		tfxRGB() { r = g = b = 0.f; }
		tfxRGB(float _r, float _g, float _b) : r(_r), g(_g), b(_b) { }
	};

	struct tfxHSV {
		float h, s, v;
		tfxHSV() { h = s = v = 0.f; }
		tfxHSV(float _h, float _s, float _v) : h(_h), s(_s), v(_v) { }
	};

	struct tfxRGBA {
		float r, g, b, a;
		tfxRGBA() { r = g = b = a = 1.f; }
		tfxRGBA(float _r, float _g, float _b, float _a) : r(_r), g(_g), b(_b), a(_a) { }
	};

	static inline float tfxRadians(float degrees) { return degrees * 0.01745329251994329576923690768489f; }
	static inline float tfxDegrees(float radians) { return radians * 57.295779513082320876798154814105f; }
	static inline void tfxBound(tfxVec2 s, tfxVec2 b) { if (s.x < 0.f) s.x = 0.f; if (s.y < 0.f) s.y = 0.f; if (s.x >= b.x) s.x = b.x - 1.f; if (s.y >= b.y) s.y = b.y - 1.f; }
	static inline void tfxBound3d(tfxVec3 s, tfxVec3 b) { if (s.x < 0.f) s.x = 0.f; if (s.y < 0.f) s.y = 0.f; if (s.z < 0.f) s.z = 0.f; if (s.x >= b.x) s.x = b.x - 1.f; if (s.y >= b.y) s.y = b.y - 1.f; if (s.z >= b.z) s.z = b.y - 1.f; }

	static inline float LengthVec3NoSqR(tfxVec3 const &v) {
		return v.x * v.x + v.y * v.y + v.z * v.z;
	}

	static inline float LengthVec(tfxVec3 const &v) {
		return sqrtf(LengthVec3NoSqR(v));
	}

	static inline float HasLength(tfxVec3 const &v) {
		return (v.x == 0 && v.y == 0 && v.z == 0) ? 0.f : 1.f;
	}

	static inline tfxVec3 NormalizeVec(tfxVec3 const &v) {
		if (v.x == 0 && v.y == 0 && v.z == 0) return tfxVec3(1.f, 0.f, 0.f);
		float length = LengthVec(v);
		return tfxVec3(v.x / length, v.y / length, v.z / length);
	}

	static inline tfxVec3 NormalizeVec(tfxVec3 const &v, float &length) {
		if (length == 0) return tfxVec3();
		return tfxVec3(v.x / length, v.y / length, v.z / length);
	}

	static inline tfxVec3 Cross(tfxVec3 a, tfxVec3 b) {
		tfxVec3 result;

		result.x = a.y*b.z - a.z*b.y;
		result.y = a.z*b.x - a.x*b.z;
		result.z = a.x*b.y - a.y*b.x;

		return(result);
	}

	static inline float DotProduct(const tfxVec3 &a, const tfxVec3 &b)
	{
		return (a.x * b.x + a.y * b.y + a.z * b.z);
	}

	//Quake 3 inverse square root
	static inline float tfxSqrt(float number)
	{
		long i;
		float x2, y;
		const float threehalfs = 1.5F;

		x2 = number * 0.5F;
		y = number;
		i = *(long *)&y;                       // evil floating point bit level hacking
		i = 0x5f3759df - (i >> 1);               // what the fuck? 
		y = *(float *)&i;
		y = y * (threehalfs - (x2 * y * y));   // 1st iteration

		return y;
	}

	inline tfxU32 GetLayerFromID(tfxU32 index) {
		return (index & 0xF0000000) >> 28;
	}

	inline tfxU32 GetIndexFromID(tfxU32 index) {
		return index & 0x0FFFFFFF;
	}

	inline tfxU32 SetNibbleID(tfxU32 nibble, tfxU32 index) {
		assert(nibble < 16);
		return (nibble << 28) + index;
	}

	static inline float FastLength(tfxVec3 const &v) {
		return 1.f / tfxSqrt(DotProduct(v, v));
	}

	static inline tfxVec3 FastNormalizeVec(tfxVec3 const &v) {
		return v * tfxSqrt(DotProduct(v, v));
	}

	struct tfxMatrix4 {

		tfxVec4 v[4];

		inline void Set2(float aa, float ab, float ba, float bb) {
			v[0].x = aa; v[0].y = ab;
			v[1].x = ba; v[1].y = bb;
		}

	};

	//Very simple 2D Matix
	struct tfxMatrix2 {

		float aa, ab, ba, bb;

		tfxMatrix2() :aa(1.f), ab(0.f), ba(0.f), bb(1.f) {}

		void Set(float _aa = 1.f, float _ab = 1.f, float _ba = 1.f, float _bb = 1.f) {
			aa = _aa;
			ab = _ab;
			ba = _ba;
			bb = _bb;
		}

		void Transpose() {
			float abt = ab;
			ab = ba;
			ba = abt;
		}

		void Scale(float s) {
			aa *= s;
			ab *= s;
			ba *= s;
			bb *= s;
		}

		tfxMatrix2 Transform(const tfxMatrix2 &m) {
			tfxMatrix2 r;
			r.aa = aa * m.aa + ab * m.ba; r.ab = aa * m.ab + ab * m.bb;
			r.ba = ba * m.aa + bb * m.ba; r.bb = ba * m.ab + bb * m.bb;
			return r;
		}

		tfxMatrix2 Transform(const tfxMatrix4 &m) {
			tfxMatrix2 r;
			r.aa = aa * m.v[0].x + ab * m.v[1].x; r.ab = aa * m.v[0].y + ab * m.v[1].y;
			r.ba = ba * m.v[0].x + bb * m.v[1].x; r.bb = ba * m.v[0].y + bb * m.v[1].y;
			return r;
		}

		tfxVec2 TransformVector(const tfxVec2 v) {
			tfxVec2 tv = tfxVec2(0.f, 0.f);
			tv.x = v.x * aa + v.y * ba;
			tv.y = v.x * ab + v.y * bb;
			return tv;
		}

	};

	static inline tfxVec2 mmTransformVector(const tfxMatrix4 &mat, const tfxVec2 v) {
		tfxVec2 tv = tfxVec2(0.f, 0.f);
		tv.x = v.x * mat.v[0].x + v.y * mat.v[1].x;
		tv.y = v.x * mat.v[0].y + v.y * mat.v[1].y;
		return tv;
	}

	inline tfxMatrix4 M4(float v = 1.f) {
		tfxMatrix4 R =
		{ {
			{v, 0, 0, 0},
			{0, v, 0, 0},
			{0, 0, v, 0},
			{0, 0, 0, v}},
		};
		return(R);
	}

	inline tfxMatrix4 M4(tfxVec4 a, tfxVec4 b, tfxVec4 c, tfxVec4 d) {
		tfxMatrix4 R =
		{ {
			{a.x, a.y, a.z, a.w},
			{b.x, b.y, b.z, b.w},
			{c.x, c.y, c.z, c.w},
			{d.x, d.y, d.z, d.w}},
		};
		return(R);
	}

	static inline tfxMatrix4 mmXRotate(float angle) {
		float c = std::cos(angle);
		float s = std::sin(angle);
		tfxMatrix4 r =
		{ {
			{1, 0, 0, 0},
			{0, c,-s, 0},
			{0, s, c, 0},
			{0, 0, 0, 1}}, 
		};
		return r;
	}

	static inline tfxMatrix4 mmYRotate(float angle) {
		float c = std::cos(angle);
		float s = std::sin(angle);
		tfxMatrix4 r =
		{ {
			{ c, 0, s, 0},
			{ 0, 1, 0, 0},
			{-s, 0, c, 0},
			{ 0, 0, 0, 1}}, 
		};
		return r;
	}

	static inline tfxMatrix4 mmZRotate(float angle) {
		float c = std::cos(angle);
		float s = std::sin(angle);
		tfxMatrix4 r =
		{ {
			{c, -s, 0, 0},
			{s,  c, 0, 0},
			{0,  0, 1, 0},
			{0,  0, 0, 1}}, 
		};
		return r;
	}

	static inline tfxMatrix4 mmTranslate(tfxMatrix4 const &m, tfxVec3 const &v) {
		tfxMatrix4 result;
		result.v[3] = m.v[0] * v.x + m.v[1] * v.y + m.v[2] * v.z + m.v[3];
		return result;
	}

	static inline tfxMatrix4 mmScale(tfxMatrix4 const &m, tfxVec3 const &v) {
		tfxMatrix4 result;
		result.v[0] = m.v[0] * v.x;
		result.v[1] = m.v[1] * v.y;
		result.v[2] = m.v[2] * v.z;
		result.v[3] = m.v[3];
		return result;
	}

	static inline tfxMatrix4 Transpose(tfxMatrix4 &mat) {
		return M4(
			tfxVec4(mat.v[0].x, mat.v[1].x, mat.v[2].x, mat.v[3].x),
			tfxVec4(mat.v[0].y, mat.v[1].y, mat.v[2].y, mat.v[3].y),
			tfxVec4(mat.v[0].z, mat.v[1].z, mat.v[2].z, mat.v[3].z),
			tfxVec4(mat.v[0].w, mat.v[1].w, mat.v[2].w, mat.v[3].w)
		);
	}

	static inline tfxMatrix4 mmTransform2(const tfxMatrix4 &in, const tfxMatrix4 &m) {
		tfxMatrix4 r;
		r.v[0].x = in.v[0].x * m.v[0].x + in.v[0].y * m.v[1].x; r.v[0].y = in.v[0].x * m.v[0].y + in.v[0].y * m.v[1].y;
		r.v[1].x = in.v[1].x * m.v[0].x + in.v[1].y * m.v[1].x; r.v[1].y = in.v[1].x * m.v[0].y + in.v[1].y * m.v[1].y;
		return r;
	}

	static inline tfxMatrix4 mmTransform2(const tfxMatrix4 &in, const tfxMatrix2 &m) {
		tfxMatrix4 r;
		r.v[0].x = in.v[0].x * m.aa + in.v[0].y * m.ba; r.v[0].y = in.v[0].x * m.ab + in.v[0].y * m.bb;
		r.v[1].x = in.v[1].x * m.aa + in.v[1].y * m.ba; r.v[1].y = in.v[1].x * m.ab + in.v[1].y * m.bb;
		return r;
	}

	static inline tfxMatrix4 mmTransform(const tfxMatrix4 &in, const tfxMatrix4 &m) {
		tfxMatrix4 res = M4(0.f);

		__m128 in_row[4];
		in_row[0] = _mm_load_ps(&in.v[0].x);
		in_row[1] = _mm_load_ps(&in.v[1].x);
		in_row[2] = _mm_load_ps(&in.v[2].x);
		in_row[3] = _mm_load_ps(&in.v[3].x);

		__m128 m_row1 = _mm_set_ps(m.v[3].x, m.v[2].x, m.v[1].x, m.v[0].x);
		__m128 m_row2 = _mm_set_ps(m.v[3].y, m.v[2].y, m.v[1].y, m.v[0].y);
		__m128 m_row3 = _mm_set_ps(m.v[3].z, m.v[2].z, m.v[1].z, m.v[0].z);
		__m128 m_row4 = _mm_set_ps(m.v[3].w, m.v[2].w, m.v[1].w, m.v[0].w);

		for (int r = 0; r <= 3; ++r)
		{

			__m128 row1result = _mm_mul_ps(in_row[r], m_row1);
			__m128 row2result = _mm_mul_ps(in_row[r], m_row2);
			__m128 row3result = _mm_mul_ps(in_row[r], m_row3);
			__m128 row4result = _mm_mul_ps(in_row[r], m_row4);

			float tmp[4];
			_mm_store_ps(tmp, row1result);
			res.v[r].x = tmp[0] + tmp[1] + tmp[2] + tmp[3];
			_mm_store_ps(tmp, row2result);
			res.v[r].y = tmp[0] + tmp[1] + tmp[2] + tmp[3];
			_mm_store_ps(tmp, row3result);
			res.v[r].z = tmp[0] + tmp[1] + tmp[2] + tmp[3];
			_mm_store_ps(tmp, row4result);
			res.v[r].w = tmp[0] + tmp[1] + tmp[2] + tmp[3];

		}
		return res;
	}

	static inline tfxVec4 mmTransformVector(const tfxMatrix4 &mat, const tfxVec4 vec) {
		tfxVec4 v;

		__m128 v4 = _mm_set_ps(vec.w, vec.z, vec.y, vec.x);

		__m128 mrow1 = _mm_load_ps(&mat.v[0].x);
		__m128 mrow2 = _mm_load_ps(&mat.v[1].x);
		__m128 mrow3 = _mm_load_ps(&mat.v[2].x);
		__m128 mrow4 = _mm_load_ps(&mat.v[3].x);

		__m128 row1result = _mm_mul_ps(v4, mrow1);
		__m128 row2result = _mm_mul_ps(v4, mrow2);
		__m128 row3result = _mm_mul_ps(v4, mrow3);
		__m128 row4result = _mm_mul_ps(v4, mrow4);

		float tmp[4];
		_mm_store_ps(tmp, row1result);
		v.x = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		_mm_store_ps(tmp, row2result);
		v.y = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		_mm_store_ps(tmp, row3result);
		v.z = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		_mm_store_ps(tmp, row4result);
		v.w = tmp[0] + tmp[1] + tmp[2] + tmp[3];

		return v;
	}

	static inline tfxVec4 mmTransformVector(const __m128 &row1, const __m128 &row2, const __m128 &row3, const __m128 &row4, const tfxVec4 vec) {
		tfxVec4 v;

		__m128 v4 = _mm_set_ps(vec.w, vec.z, vec.y, vec.x);

		__m128 row1result = _mm_mul_ps(v4, row1);
		__m128 row2result = _mm_mul_ps(v4, row2);
		__m128 row3result = _mm_mul_ps(v4, row3);
		__m128 row4result = _mm_mul_ps(v4, row4);

		float tmp[4];
		_mm_store_ps(tmp, row1result);
		v.x = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		_mm_store_ps(tmp, row2result);
		v.y = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		_mm_store_ps(tmp, row3result);
		v.z = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		_mm_store_ps(tmp, row4result);
		v.w = tmp[0] + tmp[1] + tmp[2] + tmp[3];

		return v;
	}

	static inline tfxVec3 mmTransformVector3(const tfxMatrix4 &mat, const tfxVec4 vec) {
		tfxVec3 v;

		__m128 v4 = _mm_set_ps(vec.w, vec.z, vec.y, vec.x);

		__m128 mrow1 = _mm_load_ps(&mat.v[0].x);
		__m128 mrow2 = _mm_load_ps(&mat.v[1].x);
		__m128 mrow3 = _mm_load_ps(&mat.v[2].x);
		__m128 mrow4 = _mm_load_ps(&mat.v[3].x);

		__m128 row1result = _mm_mul_ps(v4, mrow1);
		__m128 row2result = _mm_mul_ps(v4, mrow2);
		__m128 row3result = _mm_mul_ps(v4, mrow3);
		__m128 row4result = _mm_mul_ps(v4, mrow4);

		float tmp[4];
		_mm_store_ps(tmp, row1result);
		v.x = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		_mm_store_ps(tmp, row2result);
		v.y = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		_mm_store_ps(tmp, row3result);
		v.z = tmp[0] + tmp[1] + tmp[2] + tmp[3];

		return v;
	}

	static inline tfxMatrix4 mmRotate(tfxMatrix4 const &m, float r, tfxVec3 const &v) {
		float const a = r;
		float const c = cosf(a);
		float const s = sinf(a);

		tfxVec3 axis = NormalizeVec(v);
		tfxVec3 temp = axis * (1.f - c);

		tfxMatrix4 rotate;
		rotate.v[0].x = c + temp.x * axis.x;
		rotate.v[0].y = temp.x * axis.y + s * axis.z;
		rotate.v[0].z = temp.x * axis.z - s * axis.y;

		rotate.v[1].x = temp.y * axis.x - s * axis.z;
		rotate.v[1].y = c + temp.y * axis.y;
		rotate.v[1].z = temp.y * axis.z + s * axis.x;

		rotate.v[2].x = temp.z * axis.x + s * axis.y;
		rotate.v[2].y = temp.z * axis.y - s * axis.x;
		rotate.v[2].z = c + temp.z * axis.z;

		tfxMatrix4 result = M4(1.f);
		result.v[0] = m.v[0] * rotate.v[0].x + m.v[1] * rotate.v[0].y + m.v[2] * rotate.v[0].z;
		result.v[1] = m.v[0] * rotate.v[1].x + m.v[1] * rotate.v[1].y + m.v[2] * rotate.v[1].z;
		result.v[2] = m.v[0] * rotate.v[2].x + m.v[1] * rotate.v[2].y + m.v[2] * rotate.v[2].z;
		result.v[3] = m.v[3];
		return result;
	}

	static inline float Clamp(float lower, float upper, float value) {
		if (value < lower) return lower;
		if (value > upper) return upper;
		return value;
	}

	inline tfxVec3 Clamp(float lower, float upper, tfxVec3 const &v) {
		tfxVec3 result;
		result.x = Clamp(lower, upper, v.x);
		result.y = Clamp(lower, upper, v.y);
		result.z = Clamp(lower, upper, v.z);
		return result;
	}

	inline tfxU32 Pack10bit(tfxVec3 const &v, tfxU32 extra) {
		tfxVec3 converted = Clamp(-1.f, 1.f, v) * 511.f;
		tfxUInt10bit result;
		result.data.x = (tfxU32)converted.z;
		result.data.y = (tfxU32)converted.y;
		result.data.z = (tfxU32)converted.x;
		result.data.w = extra;
		return result.pack;
	}

	static inline size_t ClampStringSize(size_t compare, size_t string_size) {
		return compare < string_size ? compare : string_size;
	}

	static inline float Distance(float fromx, float fromy, float tox, float toy) {

		float w = tox - fromx;
		float h = toy - fromy;

		return std::sqrt(w * w + h * h);

	}

	inline tfxVec2 InterpolateVec2(float tween, tfxVec2 from, tfxVec2 to) {
		return to * tween + from * (1.f - tween);
	}

	inline tfxVec3 InterpolateVec3(float tween, tfxVec3 from, tfxVec3 to) {
		return to * tween + from * (1.f - tween);
	}

	inline tfxVec4 InterpolateVec4(float tween, tfxVec4 &from, tfxVec4 &to) {
		__m128 l4 = _mm_set_ps1(tween);
		__m128 l4minus1 = _mm_set_ps1(1.f - tween);
		__m128 f4 = _mm_set_ps(from.x, from.y, from.z, from.w);
		__m128 t4 = _mm_set_ps(to.x, to.y, to.z, to.w);
		__m128 from_lerp = _mm_mul_ps(f4, l4);
		__m128 to_lerp = _mm_mul_ps(f4, l4minus1);
		__m128 result = _mm_add_ps(from_lerp, to_lerp);
		tfxVec4 vec;
		_mm_store_ps(&vec.x, result);

		return vec;
	}

	inline float Interpolatef(float tween, float from, float to) {
		return from * tween + to * (1.f - tween);
	}

	typedef union {
		__m128i m;
		int a[4];
	} m128i;

	typedef union {
		__m128 m;
		float a[4];
	} m128;

	const float gradX[] =
	{
		1,-1, 1,-1,
		1,-1, 1,-1,
		0, 0, 0, 0
	};

	const float gradY[] =
	{
		1, 1,-1,-1,
		0, 0, 0, 0,
		1,-1, 1,-1
	};

	const float gradZ[] =
	{
		0, 0, 0, 0,
		1, 1,-1,-1,
		1, 1,-1,-1
	};

	const __m128 tfxF3_4 = _mm_set_ps1(1.0f / 3.0f);
	const __m128 tfxG3_4 = _mm_set_ps1(1.0f / 6.0f);
	const __m128 tfxG32_4 = _mm_set_ps1((1.0f / 6.0f) * 2.f);
	const __m128 tfxG33_4 = _mm_set_ps1((1.0f / 6.0f) * 3.f);
	const __m128i tfxONE = _mm_set1_epi32(1);
	const __m128 tfxONEF = _mm_set1_ps(1.f);
	const __m128 tfxZERO = _mm_set1_ps(0.f);
	const __m128 tfxTHIRTYTWO = _mm_set1_ps(32.f);
	const __m128i tfxFF = _mm_set1_epi32(0xFF);
	const __m128 tfxPSIX = _mm_set_ps1(0.6f);

	static inline float Dot(float x1, float y1, float z1, float x2, float y2, float z2)
	{
		return x1 * x2 + y1 * y2 + z1 * z2;
	}

	static inline float Dot(float x1, float y1, float x2, float y2)
	{
		return x1 * x2 + y1 * y2;
	}

	static inline __m128 DotProductSIMD(const __m128 &x1, const __m128 &y1, const __m128 &z1, const __m128 &x2, const __m128 &y2, const __m128 &z2)
	{
		__m128 xx = _mm_mul_ps(x1, x2);
		__m128 yy = _mm_mul_ps(y1, y2);
		__m128 zz = _mm_mul_ps(z1, z2);
		return _mm_add_ps(xx, _mm_add_ps(yy, zz));
	}

	static const int tfxPRIME_X = 501125321;
	static const int tfxPRIME_Y = 1136930381;
	static const int tfxPRIME_Z = 1720413743;

	static inline int _fnlFastRound(float f) { return (f >= 0) ? (int)(f + 0.5f) : (int)(f - 0.5f); }

	static const float tfxGRADIENTS_3D[] =
	{
		0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
		1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
		1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
		0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
		1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
		1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
		0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
		1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
		1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
		0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
		1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
		1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
		0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
		1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
		1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
		1, 1, 0, 0,  0,-1, 1, 0, -1, 1, 0, 0,  0,-1,-1, 0
	};

	static inline int _fnlHash3D(int seed, int xPrimed, int yPrimed, int zPrimed)
	{
		int hash = seed ^ xPrimed ^ yPrimed ^ zPrimed;

		hash *= 0x27d4eb2d;
		return hash;
	}

	static inline float _fnlGradCoord3D(int seed, int xPrimed, int yPrimed, int zPrimed, float xd, float yd, float zd)
	{
		int hash = _fnlHash3D(seed, xPrimed, yPrimed, zPrimed);
		hash ^= hash >> 15;
		hash &= 63 << 2;
		return xd * tfxGRADIENTS_3D[hash] + yd * tfxGRADIENTS_3D[hash | 1] + zd * tfxGRADIENTS_3D[hash | 2];
	}

	static float _fnlSingleOpenSimplex23D(int seed, float x, float y, float z)
	{
		// 3D OpenSimplex2 case uses two offset rotated cube grids.

		/*
		 * --- Rotation moved to TransformNoiseCoordinate method ---
		 * const float R3 = (float)(2.0 / 3.0);
		 * float r = (x + y + z) * R3; // Rotation, not skew
		 * x = r - x; y = r - y; z = r - z;
		 */

		int i = _fnlFastRound(x);
		int j = _fnlFastRound(y);
		int k = _fnlFastRound(z);
		float x0 = (float)(x - i);
		float y0 = (float)(y - j);
		float z0 = (float)(z - k);

		int xNSign = (int)(-1.0f - x0) | 1;
		int yNSign = (int)(-1.0f - y0) | 1;
		int zNSign = (int)(-1.0f - z0) | 1;

		float ax0 = xNSign * -x0;
		float ay0 = yNSign * -y0;
		float az0 = zNSign * -z0;

		i *= tfxPRIME_X;
		j *= tfxPRIME_Y;
		k *= tfxPRIME_Z;

		float value = 0;
		float a = (0.6f - x0 * x0) - (y0 * y0 + z0 * z0);

		for (int l = 0; ; l++)
		{
			if (a > 0)
			{
				value += (a * a) * (a * a) * _fnlGradCoord3D(seed, i, j, k, x0, y0, z0);
			}

			float b = a + 1;
			int i1 = i;
			int j1 = j;
			int k1 = k;
			float x1 = x0;
			float y1 = y0;
			float z1 = z0;
			if (ax0 >= ay0 && ax0 >= az0)
			{
				x1 += xNSign;
				b -= xNSign * 2 * x1;
				i1 -= xNSign * tfxPRIME_X;
			}
			else if (ay0 > ax0 && ay0 >= az0)
			{
				y1 += yNSign;
				b -= yNSign * 2 * y1;
				j1 -= yNSign * tfxPRIME_Y;
			}
			else
			{
				z1 += zNSign;
				b -= zNSign * 2 * z1;
				k1 -= zNSign * tfxPRIME_Z;
			}

			if (b > 0)
			{
				value += (b * b) * (b * b) * _fnlGradCoord3D(seed, i1, j1, k1, x1, y1, z1);
			}

			if (l == 1)
				break;

			ax0 = 0.5f - ax0;
			ay0 = 0.5f - ay0;
			az0 = 0.5f - az0;

			x0 = xNSign * ax0;
			y0 = yNSign * ay0;
			z0 = zNSign * az0;

			a += (0.75f - ax0) - (ay0 + az0);

			i += (xNSign >> 1) & tfxPRIME_X;
			j += (yNSign >> 1) & tfxPRIME_Y;
			k += (zNSign >> 1) & tfxPRIME_Z;

			xNSign = -xNSign;
			yNSign = -yNSign;
			zNSign = -zNSign;

			seed = ~seed;
		}

		return value * 32.69428253173828125f;
	}

	/**
	 * @file    SimplexNoise.h
	 * @brief   A Perlin Simplex Noise C++ Implementation (1D, 2D, 3D).
	 *
	 * Copyright (c) 2014-2018 Sebastien Rombauts (sebastien.rombauts@gmail.com)
	 *
	 * Distributed under the MIT License (MIT) (See accompanying file LICENSE.txt
	 * or copy at http://opensource.org/licenses/MIT)
	 */

	 /**
	  * @brief A Perlin Simplex Noise C++ Implementation (1D, 2D, 3D, 4D).
	  */
	class tfxSimplexNoise {
	public:
		// 1D Perlin simplex noise
		static float noise(float x);
		// 2D Perlin simplex noise
		static float noise(float x, float y);
		// 3D Perlin simplex noise
		static float noise(float x, float y, float z);
		// 4 noise samples using simd
		static tfxVec4 noise4(const __m128 &x4, const __m128 &y4, const __m128 &z4);

		// Fractal/Fractional Brownian Motion (fBm) noise summation
		float fractal(size_t octaves, float x) const;
		float fractal(size_t octaves, float x, float y) const;
		float fractal(size_t octaves, float x, float y, float z) const;

		/**
		 * Constructor of to initialize a fractal noise summation
		 *
		 * @param[in] frequency    Frequency ("width") of the first octave of noise (default to 1.0)
		 * @param[in] amplitude    Amplitude ("height") of the first octave of noise (default to 1.0)
		 * @param[in] lacunarity   Lacunarity specifies the frequency multiplier between successive octaves (default to 2.0).
		 * @param[in] persistence  Persistence is the loss of amplitude between successive octaves (usually 1/lacunarity)
		 */
		explicit tfxSimplexNoise(float frequency = 1.0f,
			float amplitude = 1.0f,
			float lacunarity = 2.0f,
			float persistence = 0.5f) :
			mFrequency(frequency),
			mAmplitude(amplitude),
			mLacunarity(lacunarity),
			mPersistence(persistence) {
		}

	private:
		// Parameters of Fractional Brownian Motion (fBm) : sum of N "octaves" of noise
		float mFrequency;   ///< Frequency ("width") of the first octave of noise (default to 1.0)
		float mAmplitude;   ///< Amplitude ("height") of the first octave of noise (default to 1.0)
		float mLacunarity;  ///< Lacunarity specifies the frequency multiplier between successive octaves (default to 2.0).
		float mPersistence; ///< Persistence is the loss of amplitude between successive octaves (usually 1/lacunarity)
	};

	/*
		MIT License

		Copyright (c) 2018 Stephan Brumme

		Permission is hereby granted, free of charge, to any person obtaining a copy
		of this software and associated documentation files (the "Software"),
		to deal in the Software without restriction, including without limitation
		the rights to use, copy, modify, merge, publish, distribute, sublicense,
		and/or sell copies of the Software, and to permit persons to whom the Software
		is furnished to do so, subject to the following conditions:

		The above copyright notice and this permission notice shall be included
		in all copies or substantial portions of the Software.

		THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
		INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
		PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
		HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
		OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
		SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

	/// XXHash (64 bit), based on Yann Collet's descriptions, see http://cyan4973.github.io/xxHash/
	/** How to use:
		uint64_t myseed = 0;
		XXHash64 myhash(myseed);
		myhash.add(pointerToSomeBytes,     numberOfBytes);
		myhash.add(pointerToSomeMoreBytes, numberOfMoreBytes); // call add() as often as you like to ...
		// and compute hash:
		uint64_t result = myhash.hash();
		// or all of the above in one single line:
		uint64_t result2 = XXHash64::hash(mypointer, numBytes, myseed);
		Note: my code is NOT endian-aware !
	**/
	class tfxXXHash64
	{
	public:
		/// create new XXHash (64 bit)
		/** @param seed your seed value, even zero is a valid seed **/
		explicit tfxXXHash64(uint64_t seed)
		{
			state[0] = seed + Prime1 + Prime2;
			state[1] = seed + Prime2;
			state[2] = seed;
			state[3] = seed - Prime1;
			bufferSize = 0;
			totalLength = 0;
		}

		/// add a chunk of bytes
		/** @param  input  pointer to a continuous block of data
			@param  length number of bytes
			@return false if parameters are invalid / zero **/
		bool add(const void* input, uint64_t length)
		{
			// no data ?
			if (!input || length == 0)
				return false;

			totalLength += length;
			// byte-wise access
			const unsigned char* data = (const unsigned char*)input;

			// unprocessed old data plus new data still fit in temporary buffer ?
			if (bufferSize + length < MaxBufferSize)
			{
				// just add new data
				while (length-- > 0)
					buffer[bufferSize++] = *data++;
				return true;
			}

			// point beyond last byte
			const unsigned char* stop = data + length;
			const unsigned char* stopBlock = stop - MaxBufferSize;

			// some data left from previous update ?
			if (bufferSize > 0)
			{
				// make sure temporary buffer is full (16 bytes)
				while (bufferSize < MaxBufferSize)
					buffer[bufferSize++] = *data++;

				// process these 32 bytes (4x8)
				process(buffer, state[0], state[1], state[2], state[3]);
			}

			// copying state to local variables helps optimizer A LOT
			uint64_t s0 = state[0], s1 = state[1], s2 = state[2], s3 = state[3];
			// 32 bytes at once
			while (data <= stopBlock)
			{
				// local variables s0..s3 instead of state[0]..state[3] are much faster
				process(data, s0, s1, s2, s3);
				data += 32;
			}
			// copy back
			state[0] = s0; state[1] = s1; state[2] = s2; state[3] = s3;

			// copy remainder to temporary buffer
			bufferSize = stop - data;
			for (uint64_t i = 0; i < bufferSize; i++)
				buffer[i] = data[i];

			// done
			return true;
		}

		/// get current hash
		/** @return 64 bit XXHash **/
		uint64_t hash() const
		{
			// fold 256 bit state into one single 64 bit value
			uint64_t result;
			if (totalLength >= MaxBufferSize)
			{
				result = rotateLeft(state[0], 1) +
					rotateLeft(state[1], 7) +
					rotateLeft(state[2], 12) +
					rotateLeft(state[3], 18);
				result = (result ^ processSingle(0, state[0])) * Prime1 + Prime4;
				result = (result ^ processSingle(0, state[1])) * Prime1 + Prime4;
				result = (result ^ processSingle(0, state[2])) * Prime1 + Prime4;
				result = (result ^ processSingle(0, state[3])) * Prime1 + Prime4;
			}
			else
			{
				// internal state wasn't set in add(), therefore original seed is still stored in state2
				result = state[2] + Prime5;
			}

			result += totalLength;

			// process remaining bytes in temporary buffer
			const unsigned char* data = buffer;
			// point beyond last byte
			const unsigned char* stop = data + bufferSize;

			// at least 8 bytes left ? => eat 8 bytes per step
			for (; data + 8 <= stop; data += 8)
				result = rotateLeft(result ^ processSingle(0, *(uint64_t*)data), 27) * Prime1 + Prime4;

			// 4 bytes left ? => eat those
			if (data + 4 <= stop)
			{
				result = rotateLeft(result ^ (*(uint32_t*)data) * Prime1, 23) * Prime2 + Prime3;
				data += 4;
			}

			// take care of remaining 0..3 bytes, eat 1 byte per step
			while (data != stop)
				result = rotateLeft(result ^ (*data++) * Prime5, 11) * Prime1;

			// mix bits
			result ^= result >> 33;
			result *= Prime2;
			result ^= result >> 29;
			result *= Prime3;
			result ^= result >> 32;
			return result;
		}


		/// combine constructor, add() and hash() in one static function (C style)
		/** @param  input  pointer to a continuous block of data
			@param  length number of bytes
			@param  seed your seed value, e.g. zero is a valid seed
			@return 64 bit XXHash **/
		static uint64_t hash(const void* input, uint64_t length, uint64_t seed)
		{
			tfxXXHash64 hasher(seed);
			hasher.add(input, length);
			return hasher.hash();
		}

	private:
		/// magic constants :-)
		static const uint64_t Prime1 = 11400714785074694791ULL;
		static const uint64_t Prime2 = 14029467366897019727ULL;
		static const uint64_t Prime3 = 1609587929392839161ULL;
		static const uint64_t Prime4 = 9650029242287828579ULL;
		static const uint64_t Prime5 = 2870177450012600261ULL;

		/// temporarily store up to 31 bytes between multiple add() calls
		static const uint64_t MaxBufferSize = 31 + 1;

		uint64_t      state[4];
		unsigned char buffer[MaxBufferSize];
		uint64_t      bufferSize;
		uint64_t      totalLength;

		/// rotate bits, should compile to a single CPU instruction (ROL)
		static inline uint64_t rotateLeft(uint64_t x, unsigned char bits)
		{
			return (x << bits) | (x >> (64 - bits));
		}

		/// process a single 64 bit value
		static inline uint64_t processSingle(uint64_t previous, uint64_t input)
		{
			return rotateLeft(previous + input * Prime2, 31) * Prime1;
		}

		/// process a block of 4x4 bytes, this is the main part of the XXHash32 algorithm
		static inline void process(const void* data, uint64_t& state0, uint64_t& state1, uint64_t& state2, uint64_t& state3)
		{
			const uint64_t* block = (const uint64_t*)data;
			state0 = processSingle(state0, block[0]);
			state1 = processSingle(state1, block[1]);
			state2 = processSingle(state2, block[2]);
			state3 = processSingle(state3, block[3]);
		}
	};
	//End of xxHash code

	int FormatString(char* buf, size_t buf_size, const char* fmt, va_list args);

	//Very simple string builder
	struct tfxStr {
		char *data;
		tfxU32 current_size;
		tfxU32 capacity;
		bool is_local_buffer;

		inline tfxStr() { current_size = capacity = 0; data = NULL; is_local_buffer = false; }
		inline ~tfxStr() { if (data && !is_local_buffer) { free(data); data = NULL; } current_size = capacity = 0; }

		inline bool			empty() { return current_size == 0; }
		inline char&           operator[](tfxU32 i) { return data[i]; }
		inline const char&     operator[](tfxU32 i) const { assert(i < current_size); return data[i]; }

		inline void         free_all() { if (data) { current_size = capacity = 0; free(data); data = NULL; } }
		inline void         Clear() { if (data) { current_size = 0; } }
		inline char*           begin() { return data; }
		inline const char*     begin() const { return data; }
		inline char*           end() { return data + current_size; }
		inline const char*     end() const { return data + current_size; }
		inline char&           back() { assert(current_size > 0); return data[current_size - 1]; }
		inline const char&     back() const { assert(current_size > 0); return data[current_size - 1]; }
		inline void         pop() { assert(current_size > 0); current_size--; }
		inline void	        push_back(const char v) { if (current_size == capacity) reserve(_grow_capacity(current_size + 1)); new((void*)(data + current_size)) char(v); current_size++; }

		inline tfxU32       _grow_capacity(tfxU32 sz) const { tfxU32 new_capacity = capacity ? (capacity + capacity / 2) : 8; return new_capacity > sz ? new_capacity : sz; }
		inline void         resize(tfxU32 new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
		inline void         reserve(tfxU32 new_capacity) { 
			if (new_capacity <= capacity) return; 
			char* new_data = (char*)malloc((size_t)new_capacity * sizeof(char)); 
			if (data && !is_local_buffer) { 
				memcpy(new_data, data, (size_t)current_size * sizeof(char)); 
				free(data); 
			}
			else if (is_local_buffer) {
				memcpy(new_data, strbuffer(), (size_t)current_size * sizeof(char)); 
			}
			data = new_data; 
			is_local_buffer = false;
			capacity = new_capacity; 
		}

		tfxStr(const char *text) : data(NULL), current_size(0), capacity(0), is_local_buffer(false) { size_t length = strnlen_s(text, 512); if (!length) { Clear(); return; }; if (capacity < length) reserve((tfxU32)length); memcpy(data, text, length); current_size = (tfxU32)length; NullTerminate(); }
		tfxStr(const tfxStr &src) : data(NULL), current_size(0), capacity(0), is_local_buffer(false) { size_t length = src.Length(); if (!length) { Clear(); return; }; if (capacity < length) reserve((tfxU32)length); memcpy(data, src.data, length); current_size = (tfxU32)length; NullTerminate(); }
		inline void operator=(const char *text) { size_t length = strnlen_s(text, 512); if (!length) { Clear(); return; }; if (capacity < length) reserve((tfxU32)length); memcpy(data, text, length); current_size = (tfxU32)length; NullTerminate(); }
		inline void operator=(const tfxStr& src) { Clear(); resize(src.current_size); memcpy(data, src.strbuffer(), (size_t)current_size * sizeof(char)); }
		inline bool operator==(const char *string) { return !strcmp(string, c_str()); }
		inline bool operator==(const tfxStr string) { return !strcmp(c_str(), string.c_str()); }
		inline bool operator!=(const char *string) { return strcmp(string, c_str()); }
		inline bool operator!=(const tfxStr string) { return strcmp(c_str(), string.c_str()); }
		inline const char *strbuffer() const { return is_local_buffer ? (char*)this + sizeof(tfxStr) : data; }
		inline char *strbuffer() { return is_local_buffer ? (char*)this + sizeof(tfxStr) : data; }
		inline const char *c_str() const { return current_size ? strbuffer() : ""; }
		int Find(const char *needle);
		tfxStr Lower();
		inline tfxU32 Length() const { return current_size ? current_size - 1 : 0; }
		void AddLine(const char *format, ...);
		void Appendf(const char *format, ...);
		void Appendv(const char *format, va_list args);
		inline void Append(char c) { 
			if (current_size) {
				pop();
			}
			push_back(c); 
			NullTerminate();
		}
		inline void Pop() {
			if (!Length()) return;
			if(back() == NULL)
				pop();
			pop();
			NullTerminate();
		}
		inline void Trim(char c = 32) {
			if (!Length()) return;
			if(back() == NULL)
				pop();
			while (back() == c && current_size) {
				pop();
			}
			NullTerminate();
		}
		void NullTerminate() { push_back(NULL); }
		bool SaveToFile(const char *file_name);
		const bool IsInt() const;
		const bool IsFloat() const;
	};

#define tfxStrType(type, size)		\
	struct type : public tfxStr { \
	char buffer[size]; \
	type() { data = buffer; capacity = size; current_size = 0; is_local_buffer = true; NullTerminate(); } \
	inline void operator=(const tfxStr& src) { \
		data = buffer; \
		is_local_buffer = true; \
		capacity = size; size_t length = src.Length(); \
		if (!length) { \
			Clear(); return; \
		}; \
		resize(src.current_size); \
		memcpy(data, src.strbuffer(), length); \
		current_size = (tfxU32)length; \
		NullTerminate(); \
	} \
	inline void operator=(const type& src) { \
		data = buffer; \
		is_local_buffer = true; \
		capacity = size; size_t length = src.Length(); \
		if (!length) { \
			Clear(); return; \
		}; \
		resize(src.current_size); \
		memcpy(data, src.strbuffer(), length); \
		current_size = (tfxU32)length; \
		NullTerminate(); \
	} \
	inline void operator=(const char *text) { data = buffer; is_local_buffer = true; capacity = size; size_t length = strnlen_s(text, size); if (!length) { Clear(); return; } memcpy(data, text, length); current_size = (tfxU32)length; NullTerminate(); } \
	type(const char *text) { data = buffer; is_local_buffer = true; capacity = size; size_t length = strnlen_s(text, size); if (!length) { Clear(); return; } memcpy(data, text, length); current_size = (tfxU32)length; NullTerminate(); } \
	type(const tfxStr &src) { \
		data = buffer; \
		is_local_buffer = true; \
		capacity = size; size_t length = src.Length(); \
		if (!length) { \
			Clear(); return; \
		}; \
		resize(src.current_size); \
		memcpy(data, src.strbuffer(), length); \
		current_size = (tfxU32)length; \
		NullTerminate(); \
	} \
	type(const type &src) { \
		data = buffer; \
		is_local_buffer = true; \
		capacity = size; size_t length = src.Length(); \
		if (!length) { \
			Clear(); return; \
		}; \
		resize(src.current_size); \
		memcpy(data, src.strbuffer(), length); \
		current_size = (tfxU32)length; \
		NullTerminate(); \
	} \
	inline int Find(const char *needle) { type compare = needle; type lower = Lower(); compare = compare.Lower(); if (compare.Length() > Length()) return -1; tfxU32 pos = 0; int found = 0; while (compare.Length() + pos <= Length()) { if (strncmp(lower.data + pos, compare.data, compare.Length()) == 0) { return pos; } ++pos; } return -1; } \
	inline type Lower() { type convert = *this; for (auto &c : convert) { c = tolower(c); } return convert; } \
	};

	tfxStrType(tfxStr512, 512);
	tfxStrType(tfxStr256, 256);
	tfxStrType(tfxStr128, 128);
	tfxStrType(tfxStr64, 64);
	tfxStrType(tfxStr32, 32);
	tfxStrType(tfxStr16, 16);

	/*struct tfxStr64 : public tfxStr {
		char buffer[64]; 
		tfxStr64() { data = buffer; capacity = 64; current_size = 0; is_local_buffer = true; NullTerminate(); }
		inline void operator=(const tfxStr& src) {
			data = buffer;
			is_local_buffer = true;
			capacity = 64; size_t length = src.Length();
			if (!length) {
				Clear(); return;
			};
			resize(src.current_size);
			memcpy(data, src.strbuffer(), length);
			current_size = (tfxU32)length;
			NullTerminate();
		}
		inline void operator=(const tfxStr64& src) {
			data = buffer;
			is_local_buffer = true;
			capacity = 64; size_t length = src.Length();
			if (!length) {
				Clear(); return;
			};
			resize(src.current_size);
			memcpy(data, src.strbuffer(), length);
			current_size = (tfxU32)length;
			NullTerminate();
		}
		inline void operator=(const char *text) { data = buffer; is_local_buffer = true; capacity = 64; size_t length = strnlen_s(text, 64); if (!length) { Clear(); return; } memcpy(data, text, length); current_size = (tfxU32)length; NullTerminate(); }
		tfxStr64(const char *text) { data = buffer; is_local_buffer = true; capacity = 64; size_t length = strnlen_s(text, 64); if (!length) { Clear(); return; } memcpy(data, text, length); current_size = (tfxU32)length; NullTerminate(); }
		tfxStr64(const tfxStr &src) { 
			data = buffer; 
			is_local_buffer = true; 
			capacity = 64; size_t length = src.Length(); 
			if (!length) { 
				Clear(); return; 
			}; 
			resize(src.current_size);
			memcpy(data, src.strbuffer(), length); 
			current_size = (tfxU32)length; 
			NullTerminate(); 
		}
		tfxStr64(const tfxStr64 &src) {
			data = buffer;
			is_local_buffer = true;
			capacity = 64; size_t length = src.Length();
			if (!length) {
				Clear(); return;
			};
			resize(src.current_size);
			memcpy(data, src.strbuffer(), length);
			current_size = (tfxU32)length;
			NullTerminate();
		}
		inline int Find(const char *needle) { tfxStr64 compare = needle; tfxStr64 lower = Lower(); compare = compare.Lower(); if (compare.Length() > Length()) return -1; tfxU32 pos = 0; int found = 0; while (compare.Length() + pos <= Length()) { if (strncmp(lower.data + pos, compare.data, compare.Length()) == 0) { return pos; } ++pos; } return -1; }
		inline tfxStr64 Lower() { tfxStr64 convert = *this; for (auto &c : convert) { c = tolower(c); } return convert; }
	};*/

	//tfxStr64 StrToStr64(const tfxStr &string) {
	//}

	/*struct teststr : public tfxStr {

		char buffer[64];
		teststr() { data = buffer; capacity = 64; current_size = 0; is_local_buffer = true; NullTerminate(); }
		teststr(const char *text) {
			data = buffer;
			is_local_buffer = true;
			capacity = 64;
			size_t length = strnlen_s(text, 64);
			if (!length) {
				Clear();
				return;
			}
			memcpy(data, text, length);
			current_size = (tfxU32)length;
			NullTerminate();
		}
		teststr(const tfxStr &src) { 
			data = buffer;
			is_local_buffer = true;
			capacity = 64;
			size_t length = src.Length(); 
			if (!length) { 
				Clear(); 
				return; 
			}; 
			if (capacity < length) {
				memcpy(data, src.data, capacity - 1);
			}
			else {
				memcpy(data, src.data, length);
			}
			current_size = (tfxU32)length; 
			NullTerminate(); 
		}
		int Find(const char *needle) {

			teststr compare = needle;
			teststr lower = Lower();
			compare = compare.Lower();
			if (compare.Length() > Length()) return -1;
			tfxU32 pos = 0;
			int found = 0;
			while (compare.Length() + pos <= Length()) {

				if (strncmp(lower.data + pos, compare.data, compare.Length()) == 0) {

					return pos;
				}
				++pos;
			}
			return -1;
		}
		teststr Lower() {

			teststr convert = *this;
			for (auto &c : convert) {

				c = tolower(c);
			}
			return convert;
		}
	};*/

	//Simple storage map for storing things by key/pair. The data will be in order that you add items, but the map will be in key order so just do a foreach on the data
	//and use At() to retrieve data items by name use [] overload to fetch by index if you have that.
	//Should not be used to constantly insert/remove things every frame, it's designed for setting up lists and fetching values in loops (by index preferably), and modifying based on user interaction or setting up new situation.
	//Note that if you reference things by index and you then remove something then that index may not be valid anymore so you would need to keep checks on that.
	template<typename T>
	struct tfxStorageMap {
		struct pair {
			tfxKey key;
			tfxU32 index;
			pair(tfxKey k, tfxU32 i) : key(k), index(i) {}
		};

		tfxvec<pair> map;
		tfxvec<T> data;
		void(*remove_callback)(T &item) = nullptr;

		tfxStorageMap() : map(tfxCONSTRUCTOR_VEC_INIT("Storage Map map")), data(tfxCONSTRUCTOR_VEC_INIT("Storage Map data")) {}
		tfxStorageMap(const char *map_tracker, const char *data_tracker) : map(tfxCONSTRUCTOR_VEC_INIT2(map_tracker)), data(tfxCONSTRUCTOR_VEC_INIT2(data_tracker)) {}

		//Insert a new T value into the storage
		inline void Insert(const char *name, const T &value) {
			tfxKey key = tfxXXHash64::hash(name, strlen(name), 0);
			SetIndex(key, value);
		}

		//Insert a new T value into the storage
		inline void Insert(tfxStr name, const T &value) {
			tfxKey key = tfxXXHash64::hash(name.c_str(), name.Length(), 0);
			SetIndex(key, value);
		}

		//Insert a new T value into the storage
		void Insert(tfxKey key, const T &value) {
			SetIndex(key, value);
		}

		//Insert a new T value into the storage
		inline void InsertByInt(int name, const T &value) {
			tfxKey key = name;
			SetIndex(key, value);
		}

		inline void Clear() {
			data.clear();
			map.clear();
		}

		inline void FreeAll() {
			data.free_all();
			map.free_all();
		}

		inline tfxU32 Size() {
			return data.current_size;
		}

		inline tfxU32 LastIndex() {
			return data.current_size - 1;
		}

		inline bool ValidIndex(tfxU32 index) {
			return index < data.current_size;
		}

		inline bool ValidName(const char *name) {
			return GetIndex(name) > -1;
		}

		inline bool ValidKey(tfxKey key) {
			return GetIndex(key) > -1;
		}

		inline bool ValidIntName(tfxU32 name) {
			return GetIntIndex(name) > -1;
		}

		inline bool ValidName(const tfxStr &name) {
			return GetIndex(name) > -1;
		}

		//Remove an item from the data. Slow function, 2 memmoves and then the map has to be iterated and indexes reduced by one
		//to re-align them
		inline void Remove(const char *name) {
			tfxKey key = tfxXXHash64::hash(name, strlen(name), 0);
			pair *it = LowerBound(key);
			if (remove_callback)
				remove_callback(data[it->index]);
			tfxU32 index = it->index;
			T* list_it = &data[index];
			map.erase(it);
			data.erase(list_it);
			for (auto &p : map) {
				if (p.index < index) continue;
				p.index--;
			}
		}

		//Remove an item from the data. Slow function, 2 memmoves and then the map has to be iterated and indexes reduced by one
		//to re-align them
		inline void Remove(const tfxKey &key) {
			pair *it = LowerBound(key);
			if (remove_callback)
				remove_callback(data[it->index]);
			tfxU32 index = it->index;
			T* list_it = &data[index];
			map.erase(it);
			data.erase(list_it);
			for (auto &p : map) {
				if (p.index < index) continue;
				p.index--;
			}
		}

		//Remove an item from the data. Slow function, 2 memmoves and then the map has to be iterated and indexes reduced by one
		//to re-align them
		inline void RemoveInt(int name) {
			tfxKey key = name;
			pair *it = LowerBound(key);
			if (remove_callback)
				remove_callback(data[it->index]);
			tfxU32 index = it->index;
			T* list_it = &data[index];
			map.erase(it);
			data.erase(list_it);
			for (auto &p : map) {
				if (p.index < index) continue;
				p.index--;
			}
		}

		inline T &At(const char *name) {
			int index = GetIndex(name);
			assert(index > -1);						//Key was not found
			return data[index];
		}

		inline T &At(const tfxStr &name) {
			int index = GetIndex(name.c_str());
			assert(index > -1);						//Key was not found
			return data[index];
		}

		inline T &AtInt(int name) {
			int index = GetIntIndex(name);
			assert(index > -1);						//Key was not found
			return data[index];
		}

		inline T &At(tfxKey key) {
			int index = GetIndex(key);
			assert(index > -1);						//Key was not found
			return data[index];
		}

		inline T &operator[](const tfxU32 index) {
			assert(index < data.current_size);		//Index was out of range
			return data[index];
		}

		void SetIndex(tfxKey key, const T &value) {
			pair* it = LowerBound(key);
			if (it == map.end() || it->key != key)
			{
				data.push_back(value);
				map.insert(it, pair(key, data.current_size - 1));
				return;
			}
			data[it->index] = value;
		}

		int GetIndex(const char *name) {
			tfxKey key = tfxXXHash64::hash(name, strlen(name), 0);
			pair* it = LowerBound(key);
			if (it == map.end() || it->key != key)
				return -1;
			return it->index;
		}

		int GetIntIndex(int name) {
			tfxKey key = name;
			pair* it = LowerBound(key);
			if (it == map.end() || it->key != key)
				return -1;
			return it->index;
		}

		int GetIndex(const tfxStr &name) {
			tfxKey key = tfxXXHash64::hash(name.c_str(), name.Length(), 0);
			pair* it = LowerBound(key);
			if (it == map.end() || it->key != key)
				return -1;
			return it->index;
		}

		int GetIndex(tfxKey key) {
			pair* it = LowerBound(key);
			if (it == map.end() || it->key != key)
				return -1;
			return it->index;
		}

		pair* LowerBound(tfxKey key)
		{
			tfxStorageMap::pair* first = map.data;
			tfxStorageMap::pair* last = map.data + map.current_size;
			size_t count = (size_t)(last - first);
			while (count > 0)
			{
				size_t count2 = count >> 1;
				tfxStorageMap::pair* mid = first + count2;
				if (mid->key < key)
				{
					first = ++mid;
					count -= count2 + 1;
				}
				else
				{
					count = count2;
				}
			}
			return first;
		}

	};

	static inline uint32_t Millisecs() {
		auto now = tfxClock::now().time_since_epoch();
		auto m = std::chrono::duration_cast<std::chrono::milliseconds>(now).count();
		return (uint32_t)m;
	}
	static inline uint64_t Microsecs() {
		auto now = tfxClock::now().time_since_epoch();
		auto m = std::chrono::duration_cast<std::chrono::microseconds>(now).count();
		return m;
	}

	struct tfxProfile {
		const char *name;
		tfxU32 hit_count;
		tfxU64 run_time;
		tfxU64 cycle_count;
	};

	tfxProfile tfxPROFILE_ARRAY[];

	struct tfxProfileTag {
		tfxProfile *profile;
		tfxU64 start_cycles;

		tfxProfileTag(tfxU32 id, const char *name) {
			profile = tfxPROFILE_ARRAY + id;
			profile->name = name;
			//It's surprisingly slow to use microsecs in debug mode, QueryPerformanceCounter counter is slightly faster but not enough to actually use. Release build seems fine.
			profile->run_time -= Microsecs();
			start_cycles = __rdtsc();
			AtomicAdd32(&profile->hit_count, 1);
		}

		~tfxProfileTag() {
			profile->run_time += Microsecs();
			AtomicAdd64(&profile->cycle_count, (__rdtsc() - start_cycles));
		}

	};

#ifdef tfxENABLE_PROFILING 
#define tfxPROFILE tfxProfileTag tfx_tag((tfxU32)__COUNTER__, __FUNCTION__)
#else
#define tfxPROFILE __COUNTER__
#endif

	bool EndOfProfiles();
	tfxProfile* NextProfile();

	const tfxU32 tfxMAGIC_NUMBER = '!XFT';
	const tfxU32 tfxMAGIC_NUMBER_INVENTORY = '!VNI';
	const tfxU32 tfxFILE_VERSION = 1;	//Not doing anything with this yet

	//Basic package manager used for reading/writing effects files
	struct tfxHeader {
		tfxU32 magic_number;						//Magic number to confirm file format
		tfxU32 file_version;						//The version of the file
		tfxU32 flags;								//Any flags for the file
		tfxU32 reserved0;							//Reserved for future if needed
		tfxU64 offset_to_inventory;					//Memory offset for the inventory of files
		tfxU64 reserved1;							//More reserved space
		tfxU64 reserved2;							//More reserved space
		tfxU64 reserved3;							//More reserved space
		tfxU64 reserved4;							//More reserved space
		tfxU64 reserved5;							//More reserved space
	};

	struct tfxEntryInfo {
		tfxStr file_name;							//The file name of the name stored in the package
		tfxU64 offset_from_start_of_file;			//Offset from the start of the file to where the file is located
		tfxU64 file_size;							//The size of the file
		tfxstream data;								//The file data
		
		void FreeData();
	};

	struct tfxInventory {
		tfxU32 magic_number;						//Magic number to confirm format of the Inventory
		tfxU32 entry_count;							//Number of files in the inventory
		tfxStorageMap<tfxEntryInfo> entries;		//The inventory list

		tfxInventory() :
			entries("Inventory Map", "Inventory Data")
		{}
	};

	struct tfxPackage {
		tfxStr file_path;
		tfxHeader header;
		tfxInventory inventory;
		tfxU64 file_size;							//The total file size of the package, should match file size on disk
		tfxstream file_data;						//Dump of the data from the package file on disk

		~tfxPackage();

		tfxEntryInfo *GetFile(const char *name);
		void AddFile(tfxEntryInfo file);
		void AddFile(const char *file_name, tfxstream &data);
		void Free();

	};
	
	tfxstream ReadEntireFile(const char *file_name, bool terminate = false);
	tfxErrorFlags LoadPackage(const char *file_name, tfxPackage &package);
	tfxErrorFlags LoadPackage(tfxstream &stream, tfxPackage &package);
	tfxPackage CreatePackage(const char *file_path);
	bool SavePackageDisk(tfxPackage &package);
	tfxstream SavePackageMemory(tfxPackage &package);
	tfxU64 GetPackageSize(tfxPackage &package);
	bool ValidatePackage(tfxPackage &package);

	//------------------------------------------------------------

	//Structs
	//These are mainly internal structs
	//------------------------------------------------------------

	struct tfxAttributeNode {
		float frame;
		float value;

		tfxVec2 left;
		tfxVec2 right;

		tfxAttributeNodeFlags flags;
		tfxU32 index;

		tfxAttributeNode() : frame(0.f), value(0.f), flags(0) { }
		inline bool operator==(const tfxAttributeNode& n) { return n.frame == frame && n.value == value; }

		/*
			Set the curve points for the emitterchange
			x0 and y0 are the coordinates of the point to the left of the attribute node, x1 and y1 are the coordinates to the right of the attribute node. Setting
			these will create a bezier curve.The bezier curves are restricted so that they cannot be drawn so that they loop over or behind the frame of the attribute nodes.
		*/
		void SetCurvePoints(float x0, float y0, float x1, float y1) {
			left.x = x0;
			left.y = y0;
			right.x = x1;
			right.y = y1;
			flags |= tfxAttributeNodeFlags_is_curve;
		}

		/*
			Toggle whether this attribute node is curved or linear
		*/
		void ToggleCurve() {
			flags  =~ tfxAttributeNodeFlags_is_curve;
		}

		bool IsCurve() {
			return flags & tfxAttributeNodeFlags_is_curve;
		}

		bool CurvesAreInitialised() {
			return flags & tfxAttributeNodeFlags_curves_initialised;
		}

		bool SetCurveInitialised() {
			return flags |= tfxAttributeNodeFlags_curves_initialised;
		}

		float GetX();
		float GetY();


	};

	struct tfxRandom {
		uint64_t seeds[2];
		tfxRandom();

		void ReSeed();
		void ReSeed(uint64_t seed1, uint64_t seed2);

		inline float Generate() {
			uint64_t s1 = seeds[0];
			uint64_t s0 = seeds[1];
			uint64_t result = s0 + s1;
			seeds[0] = s0;
			s1 ^= s1 << 23; // a
			seeds[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5); // b, c
			return float((double)result / TWO64f);
		}

		inline float Range(float max) {
			return Generate() * max;
		};

		inline float Range(float from, float to) {
			float a = Generate();
			float range = to - from;
			return to - range * a;
		};

		inline int RangeInt(float from, float to) {
			float a = (to - from) * Generate() + (to - from);
			return a < 0 ? int(a - 0.5f) : int(a + 0.5f);
		};

		inline tfxU32 RangeUInt(tfxU32 max) {
			float a = Generate() * (float)max;
			return tfxU32(a);
		};

	};

	static tfxRandom random_generation;

	struct tfxGraphLookup {
		tfxArray<float> values;
		tfxU32 last_frame;
		float life;

		tfxGraphLookup() : last_frame(0), life(0) {}
	};

	struct tfxGraphID {
		tfxGraphCategory category;
		tfxGraphType type = tfxGraphMaxIndex;
		tfxU32 graph_id = 0;
		tfxU32 node_id = 0;
	};

	struct tfxGraphLookupIndex {
		tfxU32 start_index;
		tfxU32 length;
		float max_life;
		float padding1;
	};
	
	//This struct is used to store indexing data in order to index into large lists containing either the node data of graphs
	//or the lookup data of compiled graphs. This is so that we can upload that data into a buffer on the GPU to get the particles
	//updating in a compute shader.
	struct tfxEffectLookUpData {
		tfxGraphLookupIndex overtime_velocity;
		tfxGraphLookupIndex overtime_width;
		tfxGraphLookupIndex overtime_height;
		tfxGraphLookupIndex overtime_weight;
		tfxGraphLookupIndex overtime_spin;
		tfxGraphLookupIndex overtime_stretch;
		tfxGraphLookupIndex overtime_red;
		tfxGraphLookupIndex overtime_green;
		tfxGraphLookupIndex overtime_blue;
		tfxGraphLookupIndex overtime_opacity;
		tfxGraphLookupIndex overtime_velocity_turbulance;
		tfxGraphLookupIndex overtime_direction_turbulance;
		tfxGraphLookupIndex overtime_velocity_adjuster;
		tfxGraphLookupIndex overtime_intensity;
		tfxGraphLookupIndex overtime_direction;
		tfxGraphLookupIndex overtime_noise_resolution;
	};

	struct tfxGraph {
		//The ratio to transalte graph frame/value to grid x/y coords on a graph editor

		tfxVec2 min;
		tfxVec2 max;
		tfxGraphPreset graph_preset;
		tfxGraphType type;
		tfxEffectEmitter *effector;
		tfxBucketArray<tfxAttributeNode> nodes;
		tfxGraphLookup lookup;
		tfxU32 index;

		tfxGraph();
		tfxGraph(tfxMemoryArenaManager *node_allocator, tfxU32 bucket_size);
		~tfxGraph();

		tfxAttributeNode* AddNode(float frame, float value, tfxAttributeNodeFlags flags = 0, float x1 = 0, float y1 = 0, float x2 = 0, float y2 = 0);
		void AddNode(tfxAttributeNode &node);
		void SetNode(uint32_t index, float frame, float value, tfxAttributeNodeFlags flags = 0, float x1 = 0, float y1 = 0, float x2 = 0, float y2 = 0);
		float GetValue(float age);
		float GetRandomValue(float age);
		float GetValue(float age, float life);
		tfxAttributeNode *GetNextNode(tfxAttributeNode &node);
		tfxAttributeNode *GetPrevNode(tfxAttributeNode &node);
		tfxAttributeNode *GetLastNode();
		float GetFirstValue();
		tfxAttributeNode* AddCoordNode(float, float);
		tfxAttributeNode* InsertCoordNode(float, float);
		tfxAttributeNode* InsertNode(float, float);
		float *LinkFirstValue();
		float GetLastValue();
		float GetMaxValue();
		float GetMinValue();
		float GetLastFrame();
		tfxBucketArray<tfxAttributeNode>& Nodes();
		tfxAttributeNode* FindNode(const tfxAttributeNode &n);
		void ValidateCurves();
		void DeleteNode(const tfxAttributeNode &n);
		void Reset(float first_node_value, tfxGraphPreset preset, bool add_node = true);
		void DragValues(tfxGraphPreset preset, float &frame, float &value);
		void Clear();
		void Free();
		void Copy(tfxGraph &to);
		bool Sort();
		void ReIndex();
		tfxVec2 GetInitialZoom();
		tfxVec2 GetInitialZoom3d();
		bool IsOvertimeGraph();
		bool IsGlobalGraph();
		bool IsAngleGraph();
		void MultiplyAllValues(float scalar);
		void CopyToNoLookups(tfxGraph *graph);

	};

	tfxVec4 GetMinMaxGraphValues(tfxGraphPreset preset);

	//todo:: Inline a lot of these.
	tfxVec2 GetQuadBezier(tfxVec2 p0, tfxVec2 p1, tfxVec2 p2, float t, float ymin, float ymax, bool clamp = true);
	tfxVec2 GetCubicBezier(tfxVec2 p0, tfxVec2 p1, tfxVec2 p2, tfxVec2 p3, float t, float ymin, float ymax, bool clamp = true);
	float GetBezierValue(const tfxAttributeNode *lastec, const tfxAttributeNode &a, float t, float ymin, float ymax);
	float GetDistance(float fromx, float fromy, float tox, float toy);
	float GetVectorAngle(float, float);
	static bool CompareNodes(tfxAttributeNode &left, tfxAttributeNode &right);
	void CompileGraph(tfxGraph &graph);
	void CompileGraphOvertime(tfxGraph &graph);
	float GetMaxLife(tfxEffectEmitter &e);
	float GetMaxAmount(tfxEffectEmitter &e);
	float LookupFastOvertime(tfxGraph &graph, float age, float lifetime);
	float LookupFast(tfxGraph &graph, float frame);
	float LookupPreciseOvertime(tfxGraph &graph, float age, float lifetime);
	float LookupPrecise(tfxGraph &graph, float frame);
	float GetRandomFast(tfxGraph &graph, float frame);
	float GetRandomPrecise(tfxGraph &graph, float frame);

	//Node Manipluation
	bool SetNode(tfxGraph &graph, tfxAttributeNode &node, float, float, tfxAttributeNodeFlags flags, float = 0, float = 0, float = 0, float = 0);
	bool SetNode(tfxGraph &graph, tfxAttributeNode &node, float &frame, float &value);
	void SetCurve(tfxGraph &graph, tfxAttributeNode &node, bool is_left_curve, float &frame, float &value);
	bool MoveNode(tfxGraph &graph, tfxAttributeNode &node, float frame, float value, bool sort = true);
	bool SetNodeFrame(tfxGraph &graph, tfxAttributeNode &node, float &frame);
	bool SetNodeValue(tfxGraph &graph, tfxAttributeNode &node, float &value);
	void ClampNode(tfxGraph &graph, tfxAttributeNode &node);
	void ClampCurve(tfxGraph &graph, tfxVec2 &curve, tfxAttributeNode &node);
	void ClampGraph(tfxGraph &graph);
	bool IsOvertimeGraph(tfxGraphType type);
	bool IsOvertimePercentageGraph(tfxGraphType type);
	bool IsGlobalGraph(tfxGraphType type);
	bool IsGlobalPercentageGraph(tfxGraphType type);
	bool IsAngleGraph(tfxGraphType type);
	bool IsAngleOvertimeGraph(tfxGraphType type);
	bool IsEverythingElseGraph(tfxGraphType type);

	struct tfxGlobalAttributes {
		tfxGraph life;
		tfxGraph amount;
		tfxGraph velocity;
		tfxGraph width;
		tfxGraph height;
		tfxGraph weight;
		tfxGraph spin;
		tfxGraph stretch;
		tfxGraph overal_scale;
		tfxGraph intensity;
		tfxGraph frame_rate;
		tfxGraph splatter;
		tfxGraph roll;
		tfxGraph pitch;
		tfxGraph yaw;
		tfxGraph emitter_width;
		tfxGraph emitter_height;
		tfxGraph emitter_depth;

		void Initialise(tfxMemoryArenaManager *allocator, tfxMemoryArenaManager *value_allocator, tfxU32 bucket_size = 8) {
			life.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			amount.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			velocity.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			width.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			height.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			weight.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			spin.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			stretch.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			overal_scale.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			intensity.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			frame_rate.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			splatter.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			roll.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			pitch.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			yaw.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			emitter_width.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			emitter_height.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			emitter_depth.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);

			life.lookup.values.allocator = value_allocator;
			amount.lookup.values.allocator = value_allocator;
			velocity.lookup.values.allocator = value_allocator;
			width.lookup.values.allocator = value_allocator;
			height.lookup.values.allocator = value_allocator;
			weight.lookup.values.allocator = value_allocator;
			spin.lookup.values.allocator = value_allocator;
			stretch.lookup.values.allocator = value_allocator;
			overal_scale.lookup.values.allocator = value_allocator;
			intensity.lookup.values.allocator = value_allocator;
			frame_rate.lookup.values.allocator = value_allocator;
			splatter.lookup.values.allocator = value_allocator;
			roll.lookup.values.allocator = value_allocator;
			pitch.lookup.values.allocator = value_allocator;
			yaw.lookup.values.allocator = value_allocator;
			emitter_width.lookup.values.allocator = value_allocator;
			emitter_height.lookup.values.allocator = value_allocator;
			emitter_depth.lookup.values.allocator = value_allocator;
		}

		void Free() {
			life.Free();
			amount.Free();
			velocity.Free();
			width.Free();
			height.Free();
			weight.Free();
			spin.Free();
			stretch.Free();
			overal_scale.Free();
			intensity.Free();
			frame_rate.Free();
			splatter.Free();
			roll.Free();
			pitch.Free();
			yaw.Free();
			emitter_width.Free();
			emitter_height.Free();
			emitter_depth.Free();
		}

		void CopyToNoLookups(tfxGlobalAttributes *dst) {
			life.CopyToNoLookups(&dst->life);
			amount.CopyToNoLookups(&dst->amount);
			velocity.CopyToNoLookups(&dst->velocity);
			width.CopyToNoLookups(&dst->width);
			height.CopyToNoLookups(&dst->height);
			weight.CopyToNoLookups(&dst->weight);
			spin.CopyToNoLookups(&dst->spin);
			stretch.CopyToNoLookups(&dst->stretch);
			overal_scale.CopyToNoLookups(&dst->overal_scale);
			intensity.CopyToNoLookups(&dst->intensity);
			frame_rate.CopyToNoLookups(&dst->frame_rate);
			splatter.CopyToNoLookups(&dst->splatter);
			roll.CopyToNoLookups(&dst->roll);
			pitch.CopyToNoLookups(&dst->pitch);
			yaw.CopyToNoLookups(&dst->yaw);
			emitter_width.CopyToNoLookups(&dst->emitter_width);
			emitter_height.CopyToNoLookups(&dst->emitter_height);
			emitter_depth.CopyToNoLookups(&dst->emitter_depth);
		}

	};

	struct tfxPropertyAttributes {
		tfxGraph emission_pitch;
		tfxGraph emission_yaw;
		tfxGraph emission_range;
		tfxGraph roll;
		tfxGraph pitch;
		tfxGraph yaw;
		tfxGraph splatter;
		tfxGraph emitter_width;
		tfxGraph emitter_height;
		tfxGraph emitter_depth;
		tfxGraph arc_size;
		tfxGraph arc_offset;

		void Initialise(tfxMemoryArenaManager *allocator, tfxMemoryArenaManager *value_allocator, tfxU32 bucket_size = 8) {
			emission_pitch.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			emission_yaw.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			emission_range.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			roll.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			pitch.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			yaw.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			splatter.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			emitter_width.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			emitter_height.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			emitter_depth.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			arc_size.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			arc_offset.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);

			emission_pitch.lookup.values.allocator = value_allocator;
			emission_yaw.lookup.values.allocator = value_allocator;
			emission_range.lookup.values.allocator = value_allocator;
			roll.lookup.values.allocator = value_allocator;
			pitch.lookup.values.allocator = value_allocator;
			yaw.lookup.values.allocator = value_allocator;
			splatter.lookup.values.allocator = value_allocator;
			emitter_width.lookup.values.allocator = value_allocator;
			emitter_height.lookup.values.allocator = value_allocator;
			emitter_depth.lookup.values.allocator = value_allocator;
			arc_size.lookup.values.allocator = value_allocator;
			arc_offset.lookup.values.allocator = value_allocator;
		}

		void Free() {
			emission_pitch.Free();
			emission_yaw.Free();
			emission_range.Free();
			roll.Free();
			pitch.Free();
			yaw.Free();
			splatter.Free();
			emitter_width.Free();
			emitter_height.Free();
			emitter_depth.Free();
			arc_size.Free();
			arc_offset.Free();
		}

		void CopyToNoLookups(tfxPropertyAttributes *dst) {
			emission_pitch.CopyToNoLookups(&dst->emission_pitch);
			emission_yaw.CopyToNoLookups(&dst->emission_yaw);
			emission_range.CopyToNoLookups(&dst->emission_range);
			roll.CopyToNoLookups(&dst->roll);
			pitch.CopyToNoLookups(&dst->pitch);
			yaw.CopyToNoLookups(&dst->yaw);
			splatter.CopyToNoLookups(&dst->splatter);
			emitter_width.CopyToNoLookups(&dst->emitter_width);
			emitter_height.CopyToNoLookups(&dst->emitter_height);
			emitter_depth.CopyToNoLookups(&dst->emitter_depth);
			arc_size.CopyToNoLookups(&dst->arc_size);
			arc_offset.CopyToNoLookups(&dst->arc_offset);
		}

	};

	struct tfxBaseAttributes {
		tfxGraph life;
		tfxGraph amount;
		tfxGraph velocity;
		tfxGraph width;
		tfxGraph height;
		tfxGraph weight;
		tfxGraph spin;
		tfxGraph noise_offset;

		void Initialise(tfxMemoryArenaManager *allocator, tfxMemoryArenaManager *value_allocator, tfxU32 bucket_size = 8) {
			life.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			amount.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			velocity.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			width.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			height.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			weight.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			spin.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			noise_offset.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);

			life.lookup.values.allocator = value_allocator;
			amount.lookup.values.allocator = value_allocator;
			velocity.lookup.values.allocator = value_allocator;
			width.lookup.values.allocator = value_allocator;
			height.lookup.values.allocator = value_allocator;
			weight.lookup.values.allocator = value_allocator;
			spin.lookup.values.allocator = value_allocator;
			noise_offset.lookup.values.allocator = value_allocator;
		}

		void Free() {
			life.Free();
			amount.Free();
			velocity.Free();
			width.Free();
			height.Free();
			weight.Free();
			spin.Free();
			noise_offset.Free();
		}

		void CopyToNoLookups(tfxBaseAttributes *dst) {
			life.CopyToNoLookups(&dst->life);
			amount.CopyToNoLookups(&dst->amount);
			velocity.CopyToNoLookups(&dst->velocity);
			width.CopyToNoLookups(&dst->width);
			height.CopyToNoLookups(&dst->height);
			weight.CopyToNoLookups(&dst->weight);
			spin.CopyToNoLookups(&dst->spin);
			noise_offset.CopyToNoLookups(&dst->noise_offset);
		}

	};

	struct tfxVariationAttributes {
		tfxGraph life;
		tfxGraph amount;
		tfxGraph velocity;
		tfxGraph width;
		tfxGraph height;
		tfxGraph weight;
		tfxGraph spin;
		tfxGraph noise_offset;
		tfxGraph noise_resolution;

		void Initialise(tfxMemoryArenaManager *allocator, tfxMemoryArenaManager *value_allocator, tfxU32 bucket_size = 8) {
			life.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			amount.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			velocity.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			width.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			height.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			weight.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			spin.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			noise_offset.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			noise_resolution.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);

			life.lookup.values.allocator = value_allocator;
			amount.lookup.values.allocator = value_allocator;
			velocity.lookup.values.allocator = value_allocator;
			width.lookup.values.allocator = value_allocator;
			height.lookup.values.allocator = value_allocator;
			weight.lookup.values.allocator = value_allocator;
			spin.lookup.values.allocator = value_allocator;
			noise_offset.lookup.values.allocator = value_allocator;
			noise_resolution.lookup.values.allocator = value_allocator;
		}

		void Free() {
			life.Free();
			amount.Free();
			velocity.Free();
			width.Free();
			height.Free();
			weight.Free();
			spin.Free();
			noise_offset.Free();
			noise_resolution.Free();
		}

		void CopyToNoLookups(tfxVariationAttributes *dst) {
			life.CopyToNoLookups(&dst->life);
			amount.CopyToNoLookups(&dst->amount);
			velocity.CopyToNoLookups(&dst->velocity);
			width.CopyToNoLookups(&dst->width);
			height.CopyToNoLookups(&dst->height);
			weight.CopyToNoLookups(&dst->weight);
			spin.CopyToNoLookups(&dst->spin);
			noise_offset.CopyToNoLookups(&dst->noise_offset);
			noise_resolution.CopyToNoLookups(&dst->noise_resolution);
		}

	};

	struct tfxOvertimeAttributes {
		tfxGraph velocity;
		tfxGraph width;
		tfxGraph height;
		tfxGraph weight;
		tfxGraph spin;
		tfxGraph stretch;
		tfxGraph red;
		tfxGraph green;
		tfxGraph blue;
		tfxGraph blendfactor;
		tfxGraph velocity_turbulance;
		tfxGraph direction_turbulance;
		tfxGraph velocity_adjuster;
		tfxGraph intensity;
		tfxGraph direction;
		tfxGraph noise_resolution;

		void Initialise(tfxMemoryArenaManager *allocator, tfxMemoryArenaManager *value_allocator, tfxU32 bucket_size = 8) {
			velocity.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			width.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			height.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			weight.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			spin.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			stretch.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			red.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			blue.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			green.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			blendfactor.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			velocity_turbulance.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			direction_turbulance.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			velocity_adjuster.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			intensity.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			direction.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			noise_resolution.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);

			velocity.lookup.values.allocator = value_allocator;
			width.lookup.values.allocator = value_allocator;
			height.lookup.values.allocator = value_allocator;
			weight.lookup.values.allocator = value_allocator;
			spin.lookup.values.allocator = value_allocator;
			stretch.lookup.values.allocator = value_allocator;
			red.lookup.values.allocator = value_allocator;
			blue.lookup.values.allocator = value_allocator;
			green.lookup.values.allocator = value_allocator;
			blendfactor.lookup.values.allocator = value_allocator;
			velocity_turbulance.lookup.values.allocator = value_allocator;
			direction_turbulance.lookup.values.allocator = value_allocator;
			velocity_adjuster.lookup.values.allocator = value_allocator;
			intensity.lookup.values.allocator = value_allocator;
			direction.lookup.values.allocator = value_allocator;
			noise_resolution.lookup.values.allocator = value_allocator;

		}

		void Free() {
			velocity.Free();
			width.Free();
			height.Free();
			weight.Free();
			spin.Free();
			stretch.Free();
			red.Free();
			green.Free();
			blue.Free();
			blendfactor.Free();
			velocity_turbulance.Free();
			direction_turbulance.Free();
			velocity_adjuster.Free();
			intensity.Free();
			direction.Free();
			noise_resolution.Free();
		}

		void CopyToNoLookups(tfxOvertimeAttributes *dst) {
			velocity.CopyToNoLookups(&dst->velocity);
			width.CopyToNoLookups(&dst->width);
			height.CopyToNoLookups(&dst->height);
			weight.CopyToNoLookups(&dst->weight);
			spin.CopyToNoLookups(&dst->spin);
			stretch.CopyToNoLookups(&dst->stretch);
			red.CopyToNoLookups(&dst->red);
			green.CopyToNoLookups(&dst->green);
			blue.CopyToNoLookups(&dst->blue);
			blendfactor.CopyToNoLookups(&dst->blendfactor);
			velocity_turbulance.CopyToNoLookups(&dst->velocity_turbulance);
			direction_turbulance.CopyToNoLookups(&dst->direction_turbulance);
			velocity_adjuster.CopyToNoLookups(&dst->velocity_adjuster);
			intensity.CopyToNoLookups(&dst->intensity);
			direction.CopyToNoLookups(&dst->direction);
			noise_resolution.CopyToNoLookups(&dst->noise_resolution);
		}

	};

	struct tfxEmitterAttributes {
		tfxPropertyAttributes properties;
		tfxBaseAttributes base;
		tfxVariationAttributes variation;
		tfxOvertimeAttributes overtime;

		void Initialise(tfxMemoryArenaManager *allocator, tfxMemoryArenaManager *value_allocator, tfxU32 bucket_size = 8) {
			properties.Initialise(allocator, value_allocator, bucket_size);
			base.Initialise(allocator, value_allocator, bucket_size);
			variation.Initialise(allocator, value_allocator, bucket_size);
			overtime.Initialise(allocator, value_allocator, bucket_size);
		}

		void Free() {
			properties.Free();
			base.Free();
			variation.Free();
			overtime.Free();
		}
	};

	static float(*lookup_overtime_callback)(tfxGraph &graph, float age, float lifetime) = LookupFastOvertime;
	static float(*lookup_callback)(tfxGraph &graph, float age) = LookupFast;
	static float(*lookup_random_callback)(tfxGraph &graph, float age) = GetRandomFast;

	struct tfxShapeData {
		char name[64];
		tfxU32 frame_count = 0;
		tfxU32 width = 0;
		tfxU32 height = 0;
		tfxU32 shape_index = 0;
		int import_filter = 0;
	};

	struct tfxBase {
		tfxVec2 size;
		float velocity;
		float spin;
		float weight;
	};

	struct tfxAnimationCameraSettings {
		tfxVec3 camera_position;
		float camera_pitch;
		float camera_yaw;
		float camera_fov;
		float camera_floor_height;
		float camera_isometric_scale;
		bool camera_isometric;
		bool camera_hide_floor;
	};
	
	struct tfxPreviewCameraSettings {
		tfxAnimationCameraSettings camera_settings;
		float effect_z_offset;
		float camera_speed;
		bool attach_effect_to_camera = false;
	};

	//this probably only needs to be in the editor, no use for it in the library? Maybe in the future as an alternative way to play back effects...
	struct tfxAnimationSettings {
		tfxVec4 bb;
		tfxVec3 position;
		tfxVec2 frame_size;
		float scale;
		float zoom;
		int frames;
		int current_frame;
		int frame_offset;
		int extra_frames_count;
		tfxU32 seed;
		bool seamless;
		bool loop;
		bool needs_recording;
		tfxU32 needs_exporting;
		float max_radius;
		tfxU32 largest_frame;
		float playback_speed;
		tfxExportColorOptions color_option;
		tfxExportOptions export_option;
		bool export_with_transparency;
		tfxAnimationCameraSettings camera_settings;
	};

	//------------------------------------------------------------

	//API structs you can access in various ways to update and render effects in realtime

	//Image data for particle shapes. This is passed into your custom ShapeLoader function for loading image textures into whatever renderer you're using
	struct tfxImageData {
		//This can be a ptr to the image texture for rendering. You must assign this in your ShapeLoader function
		void *ptr;

		//Each particle shape saved in an effect library has a unique index
		tfxU32 shape_index;
		//The size of one frame of the image
		tfxVec2 image_size;
		//Image index refers to any index that helps you look up the correct image to use. this could be an index in a texture array for example.
		tfxU32 image_index;
		//The number of frames in the image, can be one or more
		float animation_frames;
		//Maximum distance to the nearest transparent edge of the image
		float max_radius;
		int import_filter;
		tfxU32 compute_shape_index;

		//use this definition if you need more spefic data to point to the image texture in whatever renderer you're using
		//Just define tfxCUSTOM_IMAGE_DATA before you include timelinefx.h
#ifdef tfxCUSTOM_IMAGE_DATA
		tfxCUSTOM_IMAGE_DATA
#endif // tfxCUSTOM_IMAGE_DATA

		tfxImageData() :
			ptr(nullptr),
			animation_frames(1.f),
			shape_index(0),
			max_radius(0),
			import_filter(0)
		{ }
	};

	struct tfxEmitterProperties {
		//Angle added to the rotation of the particle when spawned or random angle range if angle setting is set to tfxRandom
		tfxVec3 angle_offsets;
		//When aligning the billboard along a vector, you can set the type of vector that it aligns with
		tfxVectorAlignType vector_align_type;
		//Point, area, ellipse emitter etc.
		tfxEmissionType emission_type;
		//If single shot flag is set then you can limit how many times it will loop over it's overtime graphs before expiring
		tfxU32 single_shot_limit;
		//Animation frame rate
		float frame_rate;
		//The final frame index of the animation
		float end_frame;
		//Pointer to the ImageData in the EffectLibary. 
		tfxImageData *image;
		//For 3d effects, the type of billboarding: 0 = use billboarding (always face camera), 1 = No billboarding, 2 = No billboarding and align with motion
		tfxBillboardingOptions billboard_option;

		//The number of rows/columns/ellipse/line points in the grid when spawn on grid flag is used
		tfxVec3 grid_points;
		//The rotation of particles when they spawn, or behave overtime if tfxAlign is used
		tfxAngleSettingFlags angle_settings;
		//Layer of the particle manager that the particle is added to
		tfxU32 layer;
		//Milliseconds to delay spawing
		float delay_spawning;
		//Should particles emit towards the center of the emitter or away, or in a specific direction
		tfxEmissionDirection emission_direction;

		//How particles should behave when they reach the end of the line
		tfxLineTraversalEndBehaviour end_behaviour;
		//Bit field of various boolean flags
		tfxParticleControlFlags compute_flags;
		//Offset to draw particles at
		tfxVec2 image_handle;
		//Offset of emitters
		tfxVec3 emitter_handle;
		//When single flag is set, spawn this amount of particles in one go
		tfxU32 spawn_amount;
		//The shape being used for all particles spawned from the emitter
		tfxU32 shape_index;
		//The number of millisecs before an effect or emitter will loop back round to the beginning of it's graph lookups
		float loop_length;
		//The start frame index of the animation
		float start_frame;

		tfxEmitterProperties() :
			angle_offsets(0.f, 0.f, tfx360Radians),
			image(nullptr),
			image_handle(tfxVec2()),
			spawn_amount(1),
			single_shot_limit(0),
			emission_type(tfxEmissionType::tfxPoint),
			billboard_option(tfxBillboarding),
			vector_align_type(tfxVectorAlignType_motion),
			emission_direction(tfxEmissionDirection::tfxOutwards),
			grid_points({ 10.f, 10.f, 10.f }),
			emitter_handle(),
			end_behaviour(tfxLineTraversalEndBehaviour::tfxLoop),
			loop_length(0.f),
			layer(0),
			shape_index(1),
			start_frame(0),
			end_frame(0),
			frame_rate(30.f),
			angle_settings(tfxAngleSettingFlags_random_roll | tfxAngleSettingFlags_specify_pitch | tfxAngleSettingFlags_specify_yaw),
			delay_spawning(0.f)
		{ }
	};

	struct tfxEmitterTransform {
		//Position, scale and rotation values
		tfxVec3 local_position;
		tfxVec3 world_position;
		tfxVec3 captured_position;
		tfxVec3 local_rotations;
		tfxVec3 world_rotations;
		tfxVec3 scale;
		//Todo: save space and use a quaternion here
		tfxMatrix4 matrix;
		tfxEmitterTransform() :
			matrix(M4())
		{}
	};

	struct tfxCommon {

		tfxEmitterTransform transform;
		float frame;
		float age;
		float loop_length;
		float highest_particle_age;
		float timeout_counter;
		float timeout;
		tfxU32 active_children;
		tfxVec3 handle;
		tfxEmitterStateFlags state_flags;
		tfxEmitterPropertyFlags property_flags;
		tfxEffectLibrary *library;

		tfxCommon() :
			property_flags(tfxEmitterPropertyFlags_image_handle_auto_center | 
							tfxEmitterPropertyFlags_grid_spawn_clockwise | 
							tfxEmitterPropertyFlags_emitter_handle_auto_center | 
							tfxEmitterPropertyFlags_global_uniform_size | 
							tfxEmitterPropertyFlags_base_uniform_size | 
							tfxEmitterPropertyFlags_lifetime_uniform_size),
			frame(0.f),
			age(0.f),
			state_flags(0),
			timeout_counter(0),
			timeout(100.f),
			active_children(0)
		{ }

	};

	struct tfxEmitterState {
		tfxVec3 emitter_size;
		tfxVec3 grid_coords;
		tfxVec3 grid_direction;
		tfxVec3 emission_direction_normal;	//for 2d effects, x contains the direction - not in use yet
		tfxVec2 image_handle;
		float intensity;
		float velocity_adjuster;
		float global_opacity;
		float stretch;
		float emitter_handle_y;
		float overal_scale;
		float amount_remainder;
		float emission_alternator;
		float qty;
		//The callback to transform the particles each update. This will change based on the properties of the emitter
		void(*transform_particle_callback2d)(tfxParticleData &data, tfxVec2 &world_position, float &world_rotations, const tfxCommon &common, const tfxVec3 &from_position);
		void(*transform_particle_callback3d)(tfxParticleData &data, tfxVec3 &world_position, tfxVec3 &world_rotations, const tfxCommon &common, const tfxVec3 &from_position);

		tfxEmitterState() :
			amount_remainder(0.f),
			emission_direction_normal(0.f, 1.f, 0.f)
		{}
	};

	struct tfxEmitterSpawnControls {
		float life;
		float life_variation;
		float arc_size;
		float arc_offset;
		float weight;
		float weight_variation;
		float velocity;
		float velocity_variation;
		tfxVec2 size;
		tfxVec2 size_variation;
		float spin;
		float spin_variation;
		float splatter;
		float noise_offset_variation;
		float noise_offset;
		float noise_resolution;
		tfxVec3 grid_segment_size;
	};

	//Stores the most recent parent effect (with global attributes) spawn control values to be applied to sub emitters.
	struct tfxParentSpawnControls {
		float life;
		float size_x;
		float size_y;
		float velocity;
		float spin;
		float intensity;
		float splatter;
		float weight;
	};


	float GetEmissionDirection2d(tfxCommon &common, tfxEmitterState &current, tfxEffectEmitter *library_link, tfxVec2 local_position, tfxVec2 world_position, tfxVec2 emitter_size);
	tfxVec3 GetEmissionDirection3d(tfxCommon &common, tfxEmitterState &current, tfxEffectEmitter *library_link, float emission_pitch, float emission_yaw, tfxVec3 local_position, tfxVec3 world_position, tfxVec3 emitter_size);

	struct tfxEffectEmitterInfo {
		//Name of the effect
		tfxStr64 name;				
		//Every effect and emitter in the library gets a unique id
		tfxU32 uid;
		//The max_radius of the emitter, taking into account all the particles that have spawned and active (editor only)
		float max_radius;
		//List of sub_effects ( effects contain emitters, emitters contain sub effects )
		tfxvec<tfxEffectEmitter> sub_effectors;
		//Experiment: index into the lookup index data in the effect library
		tfxU32 lookup_node_index;
		tfxU32 lookup_value_index;
		//Index to animation settings stored in the effect library. Would like to move this at some point
		tfxU32 animation_settings;
		//Index to preview camera settings stored in the effect library. Would like to move this at some point
		tfxU32 preview_camera_settings;
		//The maximum amount of life that a particle can be spawned with taking into account base + variation life values
		float max_life;
		//The estimated maximum time that the sub emitter might last for, taking into account the parent particle lifetime
		float max_sub_emitter_life;
		//The maximum amount of particles that this effect can spawn (root effects and emitters only)
		tfxU32 max_particles[tfxLAYERS];
		tfxU32 max_sub_emitters;

		tfxEffectEmitterInfo() :
			animation_settings(0),
			preview_camera_settings(0),
			max_sub_emitters(0),
			max_sub_emitter_life(0.f),
			sub_effectors(tfxCONSTRUCTOR_VEC_INIT("sub_effectors"))
		{
			for (int i = 0; i != tfxLAYERS; ++i) {
				max_particles[i] = 0;
			}
		}
	};

	//An tfxEffectEmitter can either be an effect which stores emitters and global graphs for affecting all the attributes in the emitters
	//Or it can be an emitter which spawns all of the particles. Effectors are stored in the particle manager effects list buffer.
	//This is only for library storage, when using to update each frame this is copied to tfxEffectType and tfxEmitterType, much more compact versions more
	//suited for realtime use.
	/*
	Effect Emitter realtime fields:
	common
	type
	parent
	parent_particle
	next_ptr
	properties
	flags
	current
	highest_particle_age
	spawn_controls
	global
	property
	base
	variation
	overtime
	sub_effectors
	compute_slot_id
	*/
	struct tfxEffectEmitter {
		//Required for frame by frame updating
		//The current state of the effect/emitter used in the editor only at this point
		tfxEmitterState current;
		//Common variables needed to update the effect/emitter
		tfxCommon common;
		//Is this an tfxEffectType or tfxEmitterType
		tfxEffectEmitterType type;
		//The index within the library that this exists at
		tfxU32 library_index;
		//A hash of the directory path to the effect ie Flare/spark, and also a UID for the effect/emitter
		tfxKey path_hash;
		//The current highest particle age. When using a compute buffer we don't have any reliable way of keeping track of particle counts of individual emitters, so how do we know when to remove an emitter
		//after all it's particles have expired? We set this variable to the highest particle age each time it spawns a particle and then counts it down each frame. When it's 0 then we know that there are no
		//more particles being controlled by this emitter and can therefore time it out.
		float highest_particle_age;
		//compute slot id if a compute shader is being used. Only applied to bottom emitters (emitters with no child effects)
		tfxU32 compute_slot_id;
		//All graphs that the effect uses to lookup attribute values are stored in the library. These variables here are indexes to the array where they're stored
		tfxU32 global;
		tfxU32 emitter_attributes;
		//Pointer to the immediate parent
		tfxEffectEmitter *parent;
		//Pointer to the next pointer in the particle manager buffer. 
		tfxEffectEmitter *next_ptr;
		//Pointer to the sub effect's particle that spawned it
		tfxParticle *parent_particle;
		//State flags for emitters and effects
		tfxEmitterStateFlags flags;
		tfxEffectPropertyFlags effect_flags;
		//When not using insert sort to guarantee particle order, sort passes offers a more lax way of ordering particles over a number of frames.
		//The more passes the more quickly ordered the particles will be but at a higher cost
		tfxU32 sort_passes;
		//Custom user data, can be accessed in callback functions
		void *user_data;

		tfxU32 particles_index;
		tfxU32 sprites_count;
		tfxU32 sprites_index;
		tfxU32 info_index;
		tfxU32 property_index;

		//Update callbacks that are called as the effect is updated in the particle manager. See tfxEffectTemplate
		void(*update_effect_callback)(tfxEffectEmitter &effect_emitter, tfxParentSpawnControls &spawn_controls);		//Called after the effect state has been udpated
		void(*update_emitter_callback)(tfxEffectEmitter &effect_emitter, tfxEmitterSpawnControls &spawn_controls);		//Called after the emitter state has been udpated
		void(*particle_onspawn_callback)(tfxParticle &particle, void *user_data);
		void(*particle_update_callback)(tfxParticleData &particle, void *user_data);		//Called for each particle that has been udpated, but before it's state is updated (so you can override behaviour first) 

		tfxEffectEmitter() :
			highest_particle_age(0),
			parent(nullptr),
			parent_particle(nullptr),
			user_data(nullptr),
			flags(tfxEmitterStateFlags_no_tween_this_update | tfxEmitterStateFlags_enabled),
			effect_flags(tfxEffectPropertyFlags_none),
			sort_passes(1),
			update_effect_callback(NULL),
			update_emitter_callback(NULL),
			particle_onspawn_callback(NULL),
			particle_update_callback(NULL),
			particles_index(tfxINVALID)
		{ }
		~tfxEffectEmitter();

		//API functions
		//Tell the effect to stop spawning so that eventually particles will expire and the effect will be removed from the particle manager
		inline void SoftExpire();

		void SetUserData(void *data);
		void *GetUserData();
		void SetTimeout(float frames);

		tfxEffectEmitterInfo &GetInfo();
		tfxEmitterProperties &GetProperties();

		//Override graph functions for use in update_callback
		//Some of these change the same state and property values, but they're named differently just to make it clearer as to whether you're overriding kEffect or a kEmitter.

		//Internal functions
		tfxEffectEmitter& AddEmitter(tfxEffectEmitter &e);
		tfxEffectEmitter& AddEffect(tfxEffectEmitter &e);
		tfxEffectEmitter& AddEffect();
		tfxEffectEmitter& AddEffector(tfxEffectEmitterType type = tfxEmitterType);
		tfxEffectEmitter* GetRootEffect();
		bool IsRootEffect();
		void ReIndex();
		void CountChildren(int &emitters, int &effects);
		void ResetParents();
		tfxEffectEmitter* MoveUp(tfxEffectEmitter &effect);
		tfxEffectEmitter* MoveDown(tfxEffectEmitter &effect);
		void DeleteEmitter(tfxEffectEmitter *effect);
		void CleanUp();

		void ResetGlobalGraphs(bool add_node = true, bool compile = true);
		void ResetBaseGraphs(bool add_node = true, bool compile = true);
		void ResetPropertyGraphs(bool add_node = true, bool compile = true);
		void ResetVariationGraphs(bool add_node = true, bool compile = true);
		void ResetOvertimeGraphs(bool add_node = true, bool compile = true);
		void ResetEffectGraphs(bool add_node = true, bool compile = true);
		void ResetEmitterGraphs(bool add_node = true, bool compile = true);
		void UpdateMaxLife();
		void ResetAllBufferSizes();
		void UpdateAllBufferSizes();
		void UpdateAllSpriteAmounts();
		tfxU32 GetSubEffectSpriteCounts(tfxU32 layer, tfxU32 multiplier);
		float GetSubEffectLength();
		tfxU32 GetHighestQty(float parent_age);
		tfxGraph* GetGraphByType(tfxGraphType type);
		tfxU32 GetGraphIndexByType(tfxGraphType type);
		void CompileGraphs();
		void InitialiseUninitialisedGraphs();
		void SetName(const char *n);

		void ReSeed(uint64_t seed = 0);
		bool HasSingle();
		bool RenameSubEffector(tfxEffectEmitter &effect, const char *new_name);
		bool NameExists(tfxEffectEmitter &effect, const char *name);
		void FreeGraphs();
		void NoTweenNextUpdate();

		void ClearColors();
		void AddColorOvertime(float frame, tfxRGB color);
		void Clone(tfxEffectEmitter &clone, tfxEffectEmitter *root_parent, tfxEffectLibrary *destination_library, tfxEffectCloningFlags flags = 0);
		void EnableAllEmitters();
		void EnableEmitter();
		void DisableAllEmitters();
		void DisableAllEmittersExcept(tfxEffectEmitter &emitter);
		bool IsFiniteEffect();
		void FlagAs3D(bool flag);
		bool Is3DEffect();
		tfxU32 CountAllLookupValues();
		tfxParticleManagerModes GetRequiredParticleManagerMode();

	};

	struct EffectEmitterSnapShot {
		tfxEffectEmitter effect;
		tfxU32 index;
		char description[256];
		char path[512];
		bool is_current_revision = false;
		void SetDescription(const char *format, ...);
	};

	struct tfxComputeSprite {	//64 bytes
		tfxVec4 bounds;				//the min/max x,y coordinates of the image being drawn
		tfxVec4 uv;					//The UV coords of the image in the texture
		tfxVec4 scale_rotation;		//Scale and rotation (x, y = scale, z = rotation, w = multiply blend factor)
		tfxVec2 position;			//The position of the sprite
		tfxRGBA8 color;				//The color tint of the sprite
		tfxU32 parameters;	//4 extra parameters packed into a tfxU32: blend_mode, image layer index, shader function index, blend type
	};

	struct tfxControlData {
		tfxU32 flags;
		float velocity_adjuster;
		float global_intensity;
		float image_size_y;
		float image_frame_rate;
		float stretch;
		float emitter_size_y;
		float overal_scale;
		float angle_offset;
		float end_frame;
		void(*particle_update_callback)(tfxParticleData &particle, void *user_data);
		void *user_data;
		tfxOvertimeAttributes *graphs;
	};

	//This struct is only used for 3d. It's used to get the particle's spawn position which we can then use to order by depth to the camera.
	//We can then initialise the rest of the particle data with that position knowing that it will be in order.
	struct tfxSpawnPosition {
		tfxVec3 local_position;
		tfxVec3 captured_position;
		tfxVec3 world_position;
		tfxVec4 velocity_normal;
		float distance_to_camera;
		float base_weight;
		float base_velocity;
		float weight_acceleration;
	};

	// ----------Stage 1-----------
	//Attributes
	//		--> Age	
	//		-->	Noise		-->	Velocity	-->	Location
	//		--> Spin		--> Rotation	
	//		--> Size		--> Scale
	//		--> Alignment 
	//		--> Depth
	//		--> Frame
	//		--> Color

	// ----------Stage 2-----------
	//Transform

	// ----------Stage 3-----------
	//Output
	//Discard expired and write to next buffer
	//

	//Initial particle struct, looking to optimise this and make as small as possible
	//These are spawned by effector emitter types
	//Particles are stored in the particle manager particle buffer.
	struct tfxParticleData {
		tfxVec3 local_position;			//The local position of the particle, relative to the emitter.
		tfxVec3 local_rotations;
		tfxVec3 captured_position;
		//Read only when ControlParticle is called, only written to at spawn time
		tfxBase base;						//Base values created when the particle is spawned. They can be different per particle due to variations
		tfxVec4 velocity_normal;			//Current velocity direction, with stretch factor in w
		float noise_offset;					//Higer numbers means random movement is less uniform
		float noise_resolution;				//Higer numbers means random movement is more uniform
		tfxParticleFlags flags;				//flags for different states
		//Updated everyframe
		float age;							//The age of the particle, used by the controller to look up the current state on the graphs
		float max_age;						//max age before the particle expires
		tfxU32 single_loop_count;			//The number of times a single particle has looped over
		float image_frame;					//Current frame of the image if it's an animation
		float weight_acceleration;			//The current amount of gravity applied to the y axis of the particle each frame
		float intensity;					//Color is multiplied by this value in the shader to increase the brightness of the particles
		float depth;
		tfxRGBA8 color;						//Colour of the particle
	};

	struct tfxParticle {
		tfxParticle *next_ptr;
		tfxEffectEmitter *parent;
		tfxU32 sprite_index;
		tfxU32 prev_index;
		tfxParticleData data;
	};

	struct tfxSpriteTransform2d {
		tfxVec2 position;			//The position of the sprite, x, y - world, z, w = captured for interpolating
		tfxVec2 captured_position;
		tfxVec2 scale;				//Scale
		float rotation;
	};

	struct tfxParticleSprite2d {	//56 bytes
		tfxU32 image_frame;			//The image image of animation index. Set to tfxINVALID when the particle expires
		void *image_ptr;
		tfxSpriteTransform2d transform;
		tfxVec2 handle;				//Image handle offset of the sprite
		tfxRGBA8 color;				//The color tint of the sprite and blend factor in a
		float intensity;
	};

	struct tfxSpriteTransform3d {
		tfxVec3 position;			//The position of the sprite, x, y - world, z, w = captured for interpolating
		tfxVec3 captured_position;	
		tfxVec3 rotations;			//Scale and rotation (x, y = scale, z = rotation, w = intensity)
		tfxVec2 scale;				//Scale, stretch and intensity (x, y = scale, z = stretch, w = intensity)
	};

	struct tfxParticleSprite3d {	//88 bytes
		tfxU32 image_frame_plus;	//The image frame of animation index packed with alignment option flag
		void *image_ptr;
		tfxSpriteTransform3d transform;
		tfxU32 alignment;			//normalised alignment vector 3 floats packed into 10bits each with 2 bits left over
		tfxVec2 handle;				//Image handle offset of the sprite
		tfxRGBA8 color;				//The color tint of the sprite and blend factor in a
		float stretch;
		float intensity;
	};

	/*struct BillboardInstance {
		QVec4 uv;					//The UV coords of the image in the texture
		QVec4 position;				//The position of the sprite with roll in w
		QVec4 scale_pitch_yaw;		//The image scale of the billboard and pitch/yaw in z/w
		QVec2 handle;				//The handle of the billboard
		u32 alignment;				//normalised alignment vector 3 floats packed into 10bits each with 2 bits left over
		QRGBA8 color;				//The color tint of the sprite
		float stretch;				//Amount to stretch the billboard along it's alignment vector
		u32 blend_texture_array;	//reference for the texture array (8bits) and blend factor (24bits)
	};*/

	struct tfxComputeFXGlobalState {
		tfxU32 start_index = 0;
		tfxU32 current_length = 0;
		tfxU32 max_index = 0;
		tfxU32 end_index = 0;
	};

	struct tfxComputeController {
		tfxVec2 position;
		float line_length;
		float angle_offset;
		tfxVec4 scale_rotation;				//Scale and rotation (x, y = scale, z = rotation, w = velocity_adjuster)
		float end_frame;
		tfxU32 normalised_values;		//Contains normalized values which are generally either 0 or 255, normalised in the shader to 0 and 1 (except opacity): age_rate, line_negator, spin_negator, position_negator, opacity
		tfxParticleControlFlags flags;
		tfxU32 image_data_index;		//index into the shape buffer on the gpu. CopyComputeShapeData must be called to prepare the data.
		tfxVec2 image_handle;
		tfxVec2 emitter_handle;
		float noise_offset;
		float stretch;
		float frame_rate;
		float noise_resolution;
	};

	struct tfxComputeParticle {
		tfxVec2 local_position;
		tfxVec2 base_size;

		float base_velocity = 1;
		float base_spin = 1;
		float base_weight = 1;

		float age = 1;							//The age of the particle, used by the controller to look up the current state on the graphs
		float max_age = 1;						//max age before the particle expires
		float emission_angle = 1;				//Emission angle of the particle at spawn time
		float weight_acceleration = 1;			//The current amount of gravity applied to the y axis of the particle each frame

		float noise_offset = 1;					//The random velocity added each frame
		float noise_resolution = 1;				//The random velocity added each frame
		float image_frame = 0;
		tfxU32 control_slot_and_layer;	//index to the controller, and also stores the layer in the particle manager that the particle is on (layer << 3)
		float local_rotation;
	};

	struct tfxComputeImageData {
		tfxVec4 uv;
		tfxVec2 image_size;
		tfxU32 image_index;
		float animation_frames;
		//float max_radius;
	};

	//Struct to contain a static state of a particle in a frame of animation. Used in the editor for recording frames of animation so probably not needed here really!
	struct tfxParticleFrame {
		tfxVec3 position;
		tfxVec2 scale;
		tfxVec2 handle;
		tfxVec3 rotations;
		tfxU32 alignment;
		float stretch;
		tfxU32 image_frame;
		float depth;
		void *image_ptr;
		tfxRGBA8 color;
		float intensity;
		tfxU32 alignment_type;
	};

	static inline tfxParticleFrame ConvertToParticleFrame(const tfxParticle &p, tfxEmitterProperties &properties, tfxVec2 &handle) {
		tfxParticleFrame pf;
		//pf.position = p.data.world_position;
		//pf.scale = p.data.scale;
		//pf.alignment = p.data.alignment_vector;
		pf.stretch = p.data.velocity_normal.w;
		//pf.rotations = p.data.world_rotations;
		pf.alignment_type = properties.billboard_option;
		pf.handle = handle;
		pf.color = p.data.color;
		pf.intensity = p.data.intensity;
		pf.image_ptr = properties.image->ptr;
		pf.image_frame = (tfxU32)p.data.image_frame;
		return pf;
	}

	struct tfxEffect {
		tfxEffectEmitter *effect_ptr;
	};

	//Use the particle manager to add compute effects to your scene 
	struct tfxParticleManager {
		//In ordered mode, emitters get their own list of particles to update
		tfxvec<tfxring<tfxParticle>> particle_banks;
		tfxStorageMap<tfxvec<tfxEffectEmitter>> expired_emitters;
		//Only used when using distance from camera ordering. New particles are put in this list and then merge sorted into the particles buffer
		tfxvec<tfxSpawnPosition> new_positions;
		//Effects are also stored using double buffering. Effects stored here are "fire and forget", so you won't be able to apply changes to the effect in realtime. If you want to do that then 
		//you can use an tfxEffectTemplate and use callback funcitons. 
		tfxvec<tfxEffectEmitter> effects[2];
		//Set when an effect is updated and used to pass on global attributes to child emitters
		tfxParentSpawnControls parent_spawn_controls;

		//Banks of sprites for drawing in unordered mode
		tfxring<tfxParticleSprite3d> sprites3d[tfxLAYERS];
		tfxring<tfxParticleSprite2d> sprites2d[tfxLAYERS];

		//todo: document compute controllers once we've established this is how we'll be doing it.
		void *compute_controller_ptr;
		tfxvec<unsigned int> free_compute_controllers;
		unsigned int new_compute_particle_index;
		unsigned int new_particles_count;
		void *new_compute_particle_ptr;
		//The maximum number of effects that can be updated per frame in the particle manager. If you're running effects with particles that have sub effects then this number might need 
		//to be relatively high depending on your needs. Use Init to udpate the sizes if you need to. Best to call Init at the start with the max numbers that you'll need for your application and don't adjust after.
		unsigned int max_effects;
		//The maximum number of particles that can be updated per frame per layer. #define tfxLAYERS to set the number of allowed layers. This is currently 4 by default
		unsigned int max_cpu_particles_per_layer[tfxLAYERS];
		//The maximum number of particles that can be updated per frame per layer in the compute shader. #define tfxLAYERS to set the number of allowed layers. This is currently 4 by default
		unsigned int max_new_compute_particles;
		//The current sprite buffer in use, can be either 0 or 1
		unsigned int current_pbuff;
		//The current effect buffer in use, can be either 0 or 1
		unsigned int current_ebuff;

		tfxU32 sprite_index_point[tfxLAYERS];

		unsigned int max_compute_controllers;
		unsigned int highest_compute_controller_index;
		tfxComputeFXGlobalState compute_global_state;
		tfxU32 sort_passes;
		tfxU32 new_particles_index_start[tfxLAYERS];
		tfxLookupMode lookup_mode;
		//For when particles are ordered by distance from camera (3d effects)
		tfxVec3 camera_front;
		tfxVec3 camera_position;

		//These can possibly be removed at some point, they're debugging variables
		unsigned int particle_id;
		tfxEffectManagerFlags flags;

		tfxParticleManager() :
			flags(0),
			lookup_mode(tfxFast),
			max_effects(10000),
			current_ebuff(0),
			current_pbuff(0),
			highest_compute_controller_index(0),
			new_compute_particle_ptr(nullptr),
			compute_controller_ptr(nullptr),
			max_compute_controllers(10000),
			max_new_compute_particles(10000),
			new_compute_particle_index(0),
			new_particles_count(0),
			new_positions(tfxCONSTRUCTOR_VEC_INIT("new_positions")),
			free_compute_controllers(tfxCONSTRUCTOR_VEC_INIT("free_comput_controllers"))
		{ }
		~tfxParticleManager();
		tfxEffectEmitter &operator[] (unsigned int index);

		//Initialise the particle manager with the maximum number of particles and effects that you want the manager to update per frame
		void Reconfigure(tfxParticleManagerModes mode, bool is_3d);
		void InitForBoth(tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit = 1000, tfxParticleManagerModes mode = tfxParticleManagerMode_unordered);
		void InitFor2d(tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit = 1000, tfxParticleManagerModes mode = tfxParticleManagerMode_unordered);
		void InitFor3d(tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit = 1000, tfxParticleManagerModes mode = tfxParticleManagerMode_unordered);
		void InitFor2d(unsigned int effects_limit = 1000, tfxParticleManagerModes mode = tfxParticleManagerMode_unordered);
		void InitFor3d(unsigned int effects_limit = 1000, tfxParticleManagerModes mode = tfxParticleManagerMode_unordered);
		void CreateParticleBanksForEachLayer();
		//Update the particle manager. Call this once per frame in your logic udpate.
		void Update();
		//When paused you still might want to keep the particles in order:
		void UpdateParticleOrderOnly();
		//Add an effect to the particle manager. Pass a tfxEffectEmitter pointer if you want to change the effect on the fly. Once you add the effect to the particle manager
		//then it's location in the buffer will keep changing as effects are updated and added and removed. The tracker will be updated accordingly each frame so you will always
		//have access to the effect if you need it.
		void AddEffect(tfxEffectEmitter &effect, unsigned int buffer, bool is_sub_effect = false);
		void AddEffect(tfxEffectTemplate &effect);
		//Clear all effects and particles in the particle manager
		void ClearAll(bool free_memory = false);
		void FreeParticleBanks();
		//Soft expire all the effects so that the particles complete their animation first
		inline void SetCamera(float front_x, float	front_y, float front_z, float pos_x, float pos_y, float pos_z) {
			camera_front.x = front_x;
			camera_front.y = front_y;
			camera_front.z = front_z;
			camera_position.x = pos_x;
			camera_position.y = pos_y;
			camera_position.z = pos_z;
		}
		inline void ResetFlags() { flags = 0; }

		inline void ForceCapture(bool state) {
			if (state)
				flags |= tfxEffectManagerFlags_force_capture;
			else
				flags &= ~tfxEffectManagerFlags_force_capture;
		}

		inline void DisableSpawning(bool state) {
			if (state)
				flags |= tfxEffectManagerFlags_disable_spawning;
			else
				flags &= ~tfxEffectManagerFlags_disable_spawning;
		}
		void SoftExpireAll();

		//Internal use only
		int AddComputeController();
		inline void FreeComputeSlot(unsigned int slot_id) { free_compute_controllers.push_back(slot_id); }
		void EnableCompute() { flags |= tfxEffectManagerFlags_use_compute_shader; }
		void DisableCompute() { flags &= ~tfxEffectManagerFlags_use_compute_shader; }
		tfxParticle &GrabCPUParticle(unsigned int index);
		tfxComputeParticle &GrabComputeParticle(unsigned int layer);
		void ResetParticlePtr(void *ptr);
		void ResetControllerPtr(void *ptr);
		inline unsigned int GetControllerMemoryUsage() { return highest_compute_controller_index * sizeof(tfxComputeController); }
		inline unsigned int GetParticleMemoryUsage() { return new_compute_particle_index * sizeof(tfxComputeParticle); }
		void UpdateCompute(void *sampled_particles, unsigned int sample_size = 100);
		//float Record(unsigned int frames, unsigned int start_frame, std::array<tfxvec<ParticleFrame>, 1000> &particle_frames);
		inline tfxEffectEmitter* SetNextEffect(tfxEffectEmitter &e, unsigned int buffer);
		void UpdateBaseValues();
		tfxvec<tfxEffectEmitter> *GetEffectBuffer();
		void SetLookUpMode(tfxLookupMode mode);
		inline tfxParticle *SetNextParticle(unsigned int buffer_index, tfxParticle &p) {
			unsigned int index = particle_banks[buffer_index].current_size++;
			assert(index < particle_banks[buffer_index].capacity);
			particle_banks[buffer_index][index] = p;
			return &particle_banks[buffer_index][index];
		}

		inline bool FreeCapacity2d(int index, bool compute) {
			if (!compute) {
				return sprites2d[index].current_size < max_cpu_particles_per_layer[index] || flags & tfxEffectManagerFlags_dynamic_sprite_allocation;
			}
			else
				return new_compute_particle_index < max_new_compute_particles && new_compute_particle_index < compute_global_state.end_index - compute_global_state.current_length;
		}

		inline bool FreeCapacity3d(int index, bool compute) {
			if (!compute) {
				return sprites3d[index].current_size < max_cpu_particles_per_layer[index] || flags & tfxEffectManagerFlags_dynamic_sprite_allocation;
			}
			else
				return new_compute_particle_index < max_new_compute_particles && new_compute_particle_index < compute_global_state.end_index - compute_global_state.current_length;
		}

		inline bool FreeEffectCapacity() {
			return effects[0].current_size + effects[1].current_size < max_effects;
		}
		inline tfxU32 ParticleCount() { 
			tfxU32 count = 0;
			for (tfxEachLayer) {
				count += sprites2d[layer].current_size;
				count += sprites3d[layer].current_size;
			}
			return count;
		}
	};

	tfxU32 CreateParticleBank(tfxParticleManager &pm, tfxU32 reserve_amount = 100);

	void StopSpawning(tfxParticleManager &pm);
	void RemoveAllEffects(tfxParticleManager &pm);
	void AddEffect(tfxParticleManager &pm, tfxEffectEmitter &effect, float x = 0.f, float y = 0.f);
	void AddEffect(tfxParticleManager &pm, tfxEffectTemplate &effect, float x = 0.f, float y = 0.f);

	void Rotate(tfxEffectEmitter &e, float r);
	void SetAngle(tfxEffectEmitter &e, float a);
	void Scale(tfxEffectEmitter &e, const tfxVec3& s);
	void Scale(tfxEffectEmitter &e, float x, float y, float z = 1.f);
	void Position(tfxEffectEmitter &e, const tfxVec2& p);
	void Position(tfxEffectEmitter &e, const tfxVec3& p);
	void TransformEffector(tfxEffectEmitter &e, tfxSpriteTransform2d &parent, bool relative_position = true, bool relative_angle = false);
	void TransformEffector3d(tfxEffectEmitter &e, tfxSpriteTransform3d &parent, bool relative_position = true, bool relative_angle = false);
	void UpdatePMEmitter(tfxParticleManager &pm, tfxEffectEmitter &e);
	tfxU32 NewSpritesNeeded(tfxParticleManager &pm, tfxEffectEmitter &e);
	tfxU32 SpawnParticles2d(tfxParticleManager &pm, tfxEffectEmitter &e, tfxEmitterSpawnControls &spawn_controls, tfxU32 max_spawn_amount);
	tfxU32 SpawnParticles3d(tfxParticleManager &pm, tfxEffectEmitter &e, tfxEmitterSpawnControls &spawn_controls, tfxU32 max_spawn_amount);
	void InitCPUParticle2d(tfxParticleManager &pm, tfxEffectEmitter &e, tfxParticle &p, tfxSpriteTransform2d &sprite_transform, tfxEmitterSpawnControls &spawn_controls, float tween);
	void InitCPUParticle3d(tfxParticleManager &pm, tfxEffectEmitter &e, tfxParticle &p, tfxSpriteTransform3d &sprite_transform, tfxEmitterSpawnControls &spawn_controls, float tween);
	tfxEmitterSpawnControls UpdateEmitterState(tfxEffectEmitter &e, tfxParentSpawnControls &parent_spawn_controls);
	tfxParentSpawnControls UpdateEffectState(tfxEffectEmitter &e);
	bool ControlParticle(tfxParticleManager &pm, tfxParticle &p, tfxVec2 &sprite_scale , tfxEffectEmitter &e);
	void ControlParticles2d(tfxParticleManager &pm, tfxEffectEmitter &e, tfxU32 amount_spawned);
	void ControlParticles3d(tfxParticleManager &pm, tfxEffectEmitter &e, tfxU32 amount_spawned);
	void ControlParticlesOrdered2d(tfxParticleManager &pm);
	void ControlParticlesOrdered3d(tfxParticleManager &pm);
	void ControlParticlesDepthOrdered3d(tfxParticleManager &pm);

	struct tfxEffectLibraryStats {
		tfxU32 total_effects;
		tfxU32 total_sub_effects;
		tfxU32 total_emitters;
		tfxU32 total_attribute_nodes;
		tfxU32 total_node_lookup_indexes;
		tfxU32 total_shapes;
		tfxU64 required_graph_node_memory;
		tfxU64 required_graph_lookup_memory;
		tfxU32 reserved1;
		tfxU32 reserved2;
		tfxU32 reserved3;
		tfxU32 reserved4;
		tfxU32 reserved5;
		tfxU32 reserved6;
		tfxU32 reserved7;
	};

	struct tfxEffectLibrary {
		tfxMemoryArenaManager graph_node_allocator;
		tfxMemoryArenaManager graph_lookup_allocator;
		tfxStorageMap<tfxEffectEmitter*> effect_paths;
		tfxvec<tfxEffectEmitter> effects;
		tfxStorageMap<tfxImageData> particle_shapes;
		tfxvec<tfxEffectEmitterInfo> effect_infos;
		tfxvec<tfxEmitterProperties> emitter_properties;

		tfxvec<tfxGlobalAttributes> global_graphs;
		tfxvec<tfxEmitterAttributes> emitter_attributes;
		tfxvec<tfxAnimationSettings> animation_settings;
		tfxvec<tfxPreviewCameraSettings> preview_camera_settings;
		tfxvec<tfxAttributeNode> all_nodes;
		tfxvec<tfxEffectLookUpData> node_lookup_indexes;
		tfxvec<float> compiled_lookup_values;
		tfxvec<tfxGraphLookupIndex> compiled_lookup_indexes;
		tfxvec<tfxComputeImageData> shape_data;
		//This could probably be stored globally
		tfxvec<tfxVec4> graph_min_max;

		tfxvec<tfxU32> free_global_graphs;
		tfxvec<tfxU32> free_emitter_attributes;
		tfxvec<tfxU32> free_animation_settings;
		tfxvec<tfxU32> free_preview_camera_settings;
		tfxvec<tfxU32> free_properties;
		tfxvec<tfxU32> free_infos;

		//Get an effect from the library by index
		tfxEffectEmitter& operator[] (uint32_t index);
		tfxStr64 name;
		bool open_library = false;
		bool dirty = false;
		tfxStr library_file_path;
		tfxU32 uid = 0;

		tfxEffectLibrary() :
			effect_paths("EffectLib effect paths map", "EffectLib effect paths data"),
			particle_shapes("EffectLib shapes map", "EffectLib shapes data"),
			effects(tfxCONSTRUCTOR_VEC_INIT("effects")),
			effect_infos(tfxCONSTRUCTOR_VEC_INIT("effect_infos")),
			emitter_properties(tfxCONSTRUCTOR_VEC_INIT("emitter_properties")),
			global_graphs(tfxCONSTRUCTOR_VEC_INIT("global_graphs")),
			emitter_attributes(tfxCONSTRUCTOR_VEC_INIT("emitter_attributes")),
			animation_settings(tfxCONSTRUCTOR_VEC_INIT("animation_settings")),
			preview_camera_settings(tfxCONSTRUCTOR_VEC_INIT("preview_camera_settings")),
			all_nodes(tfxCONSTRUCTOR_VEC_INIT("all_nodes")),
			node_lookup_indexes(tfxCONSTRUCTOR_VEC_INIT("nodes_lookup_indexes")),
			compiled_lookup_values(tfxCONSTRUCTOR_VEC_INIT("compiled_lookup_values")),
			compiled_lookup_indexes(tfxCONSTRUCTOR_VEC_INIT("compiled_lookup_indexes")),
			shape_data(tfxCONSTRUCTOR_VEC_INIT("shape_data")),
			graph_min_max(tfxCONSTRUCTOR_VEC_INIT("graph_min_max")),
			free_global_graphs(tfxCONSTRUCTOR_VEC_INIT("free_global_graphs")),
			free_emitter_attributes(tfxCONSTRUCTOR_VEC_INIT("free_emitter_attributes")),
			free_animation_settings(tfxCONSTRUCTOR_VEC_INIT("free_animation_settings")),
			free_preview_camera_settings(tfxCONSTRUCTOR_VEC_INIT("free_preview_camera_settings")),
			free_properties(tfxCONSTRUCTOR_VEC_INIT("free_properties")),
			free_infos(tfxCONSTRUCTOR_VEC_INIT("free_infos"))
		{}

		//Todo: Inline a lot of these
		//Free everything in the library
		void Clear();
		//Get an effect in the library by it's path. So for example, if you want to get a pointer to the emitter "spark" in effect "explosion" then you could do GetEffect("explosion/spark")
		//You will need this function to apply user data and update callbacks to effects and emitters before adding the effect to the particle manager
		tfxEffectEmitter *GetEffect(tfxStr256 &path);
		tfxEffectEmitter *GetEffect(const char *path);
		//Get an effect by it's path hash key
		tfxEffectEmitter *GetEffect(tfxKey key);
		//Get and effect by it's index
		void PrepareEffectTemplate(tfxStr256 path, tfxEffectTemplate &effect);
		//Copy the shape data to a memory location, like a staging buffer ready to be uploaded to the GPU for use in a compute shader
		void BuildComputeShapeData(void* dst, tfxVec4(uv_lookup)(void *ptr, tfxComputeImageData &image_data, int offset));
		void CopyComputeShapeData(void* dst);
		void CopyLookupIndexesData(void* dst);
		void CopyLookupValuesData(void* dst);
		tfxU32 GetComputeShapeDataSizeInBytes();
		tfxU32 GetComputeShapeCount();
		tfxU32 GetLookupIndexCount();
		tfxU32 GetLookupValueCount();
		tfxU32 GetLookupIndexesSizeInBytes();
		tfxU32 GetLookupValuesSizeInBytes();

		inline void MaybeGrowProperties() {
			if (emitter_properties.current_size >= emitter_properties.capacity - 1) {
				emitter_properties.reserve(emitter_properties._grow_capacity(emitter_properties.capacity + 1));
			}
		}

		inline void MaybeGrowInfos() {
			if (effect_infos.current_size == effect_infos.capacity) {
				effect_infos.reserve(effect_infos._grow_capacity(effect_infos.current_size + 1));
			}
		}

		inline tfxEmitterProperties &GetProperties(tfxU32 index) {
			assert(emitter_properties.size() > index);
			return emitter_properties[index];
		}

		inline tfxEffectEmitterInfo &GetInfo(tfxEffectEmitter &e) {
			assert(effect_infos.size() > e.info_index);
			return effect_infos[e.info_index];
		}

		inline const tfxEffectEmitterInfo &GetInfo(const tfxEffectEmitter &e) {
			assert(effect_infos.size() > e.info_index);
			return effect_infos[e.info_index];
		}

		//Mainly internal functions
		void RemoveShape(tfxU32 shape_index);
		tfxEffectEmitter &AddEffect(tfxEffectEmitter &effect);
		tfxEffectEmitter &AddFolder(tfxStr64 &name);
		tfxEffectEmitter &AddFolder(tfxEffectEmitter &effect);
		void UpdateEffectPaths();
		void AddPath(tfxEffectEmitter &effect_emitter, tfxStr256 &path);
		void DeleteEffect(tfxEffectEmitter *effect);
		bool RenameEffect(tfxEffectEmitter &effect, const char *new_name);
		bool NameExists(tfxEffectEmitter &effect, const char *name);
		bool NameExists2(tfxEffectEmitter &effect, const char *name);
		void ReIndex();
		void UpdateParticleShapeReferences(tfxvec<tfxEffectEmitter> &effects, tfxU32 default_index);
		tfxEffectEmitter* MoveUp(tfxEffectEmitter &effect);
		tfxEffectEmitter* MoveDown(tfxEffectEmitter &effect);
		tfxU32 AddGlobal();
		tfxU32 AddEmitterAttributes();
		void FreeGlobal(tfxU32 index);
		void FreeEmitterAttributes(tfxU32 index);
		void FreeProperties(tfxU32 index);
		void FreeInfo(tfxU32 index);
		tfxU32 CountGlobalLookUpValues(tfxU32 index);
		tfxU32 CountEmitterLookUpValues(tfxU32 index);
		tfxU32 CloneGlobal(tfxU32 source_index, tfxEffectLibrary *destination_library);
		tfxU32 CloneEmitterAttributes(tfxU32 source_index, tfxEffectLibrary *destination_library);
		tfxU32 CloneInfo(tfxU32 source_index, tfxEffectLibrary *destination_library);
		tfxU32 CloneProperties(tfxU32 source_index, tfxEffectLibrary *destination_library);
		void AddEmitterGraphs(tfxEffectEmitter& effect);
		void AddEffectGraphs(tfxEffectEmitter& effect);
		tfxU32 AddAnimationSettings(tfxEffectEmitter& effect);
		tfxU32 AddPreviewCameraSettings(tfxEffectEmitter& effect);
		tfxU32 AddPreviewCameraSettings();
		tfxU32 AddEffectEmitterInfo();
		tfxU32 AddEmitterProperties();
		void UpdateEffectParticleStorage();
		void UpdateComputeNodes();
		void CompileAllGraphs();
		void CompileGlobalGraph(tfxU32 index);
		void CompileEmitterGraphs(tfxU32 index);
		void CompilePropertyGraph(tfxU32 index);
		void CompileBaseGraph(tfxU32 index);
		void CompileVariationGraph(tfxU32 index);
		void CompileOvertimeGraph(tfxU32 index);
		void CompileColorGraphs(tfxU32 index);
		void CompileGraphsOfEffect(tfxEffectEmitter &effect, tfxU32 depth = 0);
		void SetMinMaxData();
		float LookupPreciseOvertimeNodeList(tfxGraphType graph_type, int index, float age, float life);
		float LookupPreciseNodeList(tfxGraphType graph_type, int index, float age);
		float LookupFastOvertimeValueList(tfxGraphType graph_type, int index, float age, float life);
		float LookupFastValueList(tfxGraphType graph_type, int index, float age);

		//Debug stuff, used to check graphs are being properly recycled
		tfxU32 CountOfGraphsInUse();
		tfxU32 CountOfFreeGraphs();
	};

	struct tfxEffectTemplate {
		tfxStorageMap<tfxEffectEmitter*> paths;
		tfxEffectEmitter effect;

		tfxEffectTemplate() :
			paths("Effect template paths map", "Effect template paths data")
		{}
		void AddPath(tfxEffectEmitter &effect_emitter, tfxStr256 path) {
			paths.Insert(path, &effect_emitter);
			for (auto &sub : effect_emitter.common.library->GetInfo(effect_emitter).sub_effectors) {
				tfxStr256 sub_path = path;
				sub_path.Appendf("/%s", sub.common.library->GetInfo(sub).name.c_str());
				AddPath(sub, sub_path);
			}
		}

		inline tfxEffectEmitter &Effect() { return effect; }
		inline tfxEffectEmitter *Get(tfxStr256 &path) { if (paths.ValidName(path)) return paths.At(path); return nullptr; }
		inline void SetUserData(tfxStr256 &path, void *data) { if (paths.ValidName(path)) paths.At(path)->user_data = data; }
		inline void SetUserData(void *data) { effect.user_data = data; }
		void SetUserDataAll(void *data);
		inline void SetEffectUpdateCallback(tfxStr256 path, void(*update_callback)(tfxEffectEmitter &effect_emitter, tfxParentSpawnControls &spawn_controls)) { 
			assert(paths.ValidName(path));						//Path does not exist in library
			assert(paths.At(path)->type == tfxEffectType);		//Path must be path to an effect type
			paths.At(path)->update_effect_callback = update_callback;
		}
		inline void SetEmitterUpdateCallback(tfxStr256 path, void(*update_callback)(tfxEffectEmitter &effect_emitter, tfxEmitterSpawnControls &spawn_controls)) { 
			assert(paths.ValidName(path));						//Path does not exist in library
			assert(paths.At(path)->type == tfxEmitterType);		//Path must be a path to an emitter type
			paths.At(path)->update_emitter_callback = update_callback; }
		inline void SetEffectUpdateCallback(void(*update_callback)(tfxEffectEmitter &effect_emitter, tfxParentSpawnControls &spawn_controls)) { effect.update_effect_callback = update_callback; }
		void SetParticleUpdateCallback(tfxStr256 path, void(*particle_update_callback)(tfxParticleData &particle, void *user_data));
		void SetParticleOnSpawnCallback(tfxStr256 path, void(*particle_onspawn_callback)(tfxParticle &particle, void *user_data));
	};

	void SetTimeOut(tfxEffectTemplate &effect_template, float frames);

	/*
	Notes on updating effects emitters and particles:

	We're presented with a number of constraints when updating and drawing particles based on the asthetics that we want for the particles:
	1) Draw order of particles is important, especially when emitters emit alpha blended particles
	2) Effects can have sub effects, presenting a problem from a memory management point of view especially as particles with sub effects will be expiring at different times leaving holes in memory.
	3) You may want to draw effects in different orders, grouping together the sprites generated from effects so that they can be drawn in specific orders or not at all if not on screen etc.
	4) Particles need to be updated as quickly as possible so memory layout is very important.

	In an ideal world emitters maintain their own list of particles and update them in turn. This means that all the base values can be accessed and held in local variables. 
	The problem is that other emitters in the same effect would be spawning particles too so the draw order becomes important, you want the particles from each emitter to be mixed together when drawn
	so updating them by emitter means you lose that order.

	Sub effects start to complicate things from a memory management point of view. With each particle having it's own sub effect that need their own place in memory.

	The particle manager solves a lot of these issues by having 2 lists, an effects list and a particle list. First the effects list is updated and then the particle list is updated. When the particle
	list is updated the particles need to reference their parent emitters/effects and each particle may be referencing different emitters so I'm not sure how good that is from a caching point of view.
	Not great I would assume. 
	Another problem with the particle manager is that all the particles are in one list and so you can't separate out individual effects when drawing which would pose a problem when drawing the effects
	in different orders with other non particle related drawing.
	But the particle manager does make managing the memory a lot easier as you only need the effects and particle lists.

	My current conclusion is that for use in a game where you have already defined your particle effects in the editor you can use the tfxEffect approach where draw order of effects is more flexible and
	memory access is more efficient. The draw back is you can't really make an emitter spawn more particles, but the easy work around here is to create the effect with the most amount of particles spawning
	as you need and then scale them back dynamically as you need.

	But there is still a strong enough case for the particle manager for use in the editor where it's very useful to be able to dynamically grow the number of particles and sub effects as you work on new 
	effects. Use of compute shaders is also easier at this point with the particle manager. So at this point I think it's best to keep both methods even though it's more code to maintain.
	*/

	struct tfxMockEffect {
		tfxU32 timeout = 5;
		tfxU32 timeout_counter = 0;
		tfxU32 emitter_count = 0;
		float highest_particle_age = 0;
		float frame = 0.f;
		float age = 0.f;
		float amount_remainder = 0;
		float qty = 1.f;
		tfxEffectEmitter *library_link;
		tfxEffectLibrary *library;
		bool single_shot_done = false;
		bool started_spawning = false;
		tfxvec<float> particles[2];
	};

	//This is used to figure out how much memory each effect and emitter needs to draw particles so that the correct amount of memory can be assigned as each effect is used.
	struct tfxParticleMemoryTools {
		tfxU32 sprite_count[4];
		tfxU32 sub_effect_count;
		tfxU32 initial_effect_size = 0;
		tfxU32 emitters_removed = 0;
		float max_frames;
		float max_last_life;
		tfxU32 current_buffer;
		tfxvec<float> particles[tfxLAYERS][2];
		tfxvec<tfxMockEffect> effects[2];
		tfxEffectEmitter current_effect;

		tfxParticleMemoryTools() : current_buffer(0), sub_effect_count(0) {}

		void AddEffect(tfxEffectEmitter &effect);
		void GetEffectMaxFrames(tfxEffectEmitter &effect);
		void ProcessEffect(tfxEffectEmitter &effect);
		void Process();
		void MockUpdateEmitter(tfxMockEffect &emitter);
		void MockUpdateParticles();
	};

	struct tfxDataEntry {
		tfxDataType type;
		tfxStr32 key;
		tfxStr str_value;
		int int_value;
		bool bool_value;
		float float_value;
		double double_value;
	};

	struct tfxDataTypesDictionary {
		bool initialised = false;
		tfxStorageMap<tfxDataType> names_and_types;

		tfxDataTypesDictionary() : 
			names_and_types("Data Types Storage Map", "Data Types Storage Data")
		{}
		void Init();
	};

	extern tfxDataTypesDictionary data_types;

	//Internal functions
	//Some file IO functions
	bool HasDataValue(tfxStorageMap<tfxDataEntry> &config, tfxStr32 key);
	void AddDataValue(tfxStorageMap<tfxDataEntry> &config, tfxStr32 key, const char *value);
	void AddDataValue(tfxStorageMap<tfxDataEntry> &config, tfxStr32 key, int value);
	void AddDataValue(tfxStorageMap<tfxDataEntry> &config, tfxStr32 key, bool value);
	void AddDataValue(tfxStorageMap<tfxDataEntry> &config, tfxStr32 key, double value);
	void AddDataValue(tfxStorageMap<tfxDataEntry> &config, tfxStr32 key, float value);
	tfxStr &GetDataStrValue(tfxStorageMap<tfxDataEntry> &config, const char* key);
	int& GetDataIntValue(tfxStorageMap<tfxDataEntry> &config, const char* key);
	float& GetDataFloatValue(tfxStorageMap<tfxDataEntry> &config, const char* key);
	bool SaveDataFile(tfxStorageMap<tfxDataEntry> &config, const char* path = "");
	bool LoadDataFile(tfxStorageMap<tfxDataEntry> &config, const char* path);
	void StreamProperties(tfxEmitterProperties &property, tfxEmitterPropertyFlags &flags, tfxStr &file);
	void StreamProperties(tfxEffectEmitter &effect, tfxStr &file);
	void StreamGraph(const char * name, tfxGraph &graph, tfxStr &file);
	tfxvec<tfxStr64> SplitString(const tfxStr &s, char delim = 61);
	void SplitString(const tfxStr &s, tfxvec<tfxStr64> &pair, char delim = 61);
	bool StringIsUInt(const tfxStr &s);
	int GetDataType(const tfxStr &s);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, uint32_t value);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, float value);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, bool value);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, int value);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, tfxStr &value);
	void AssignGraphData(tfxEffectEmitter &effect, tfxvec<tfxStr64> &values);
	void AssignNodeData(tfxAttributeNode &node, tfxvec<tfxStr64> &values);
	static inline void Transform(tfxEmitterTransform &out, const tfxEmitterTransform &in) {
		float s = sin(out.local_rotations.roll);
		float c = cos(out.local_rotations.roll);

		out.matrix.Set2(c, s, -s, c);
		out.scale = in.scale;

		out.world_rotations.roll = in.world_rotations.roll + out.local_rotations.roll;

		out.matrix = mmTransform2(out.matrix, in.matrix);
		tfxVec2 rotatevec = mmTransformVector(in.matrix, tfxVec2(out.local_position.x, out.local_position.y));

		out.world_position = in.world_position.xy() + rotatevec * in.scale.xy();
	}
	static inline void Transform3d(tfxEmitterTransform &out, const tfxEmitterTransform &in) {
		tfxMatrix4 roll = mmZRotate(out.local_rotations.roll);
		tfxMatrix4 pitch = mmXRotate(out.local_rotations.pitch);
		tfxMatrix4 yaw = mmYRotate(out.local_rotations.yaw);
		out.matrix = mmTransform(yaw, pitch);
		out.matrix = mmTransform(out.matrix, roll);
		out.scale = in.scale;

		out.world_rotations = in.world_rotations + out.local_rotations;

		out.matrix = mmTransform(out.matrix, in.matrix);
		tfxVec3 rotatevec = mmTransformVector3(in.matrix, out.local_position);

		out.world_position = in.world_position + rotatevec;
	}
	static inline void TransformParticle(tfxParticleData &data, tfxVec2 &world_position, float &world_rotations, const tfxCommon &common, const tfxVec3 &from_position) {
		world_position = data.local_position.xy();
		world_rotations = data.local_rotations.roll;
	}
	static inline void TransformParticleAngle(tfxParticleData &data, tfxVec2 &world_position, float &world_rotations, const tfxCommon &common, const tfxVec3 &from_position) {
		world_position = data.local_position.xy();
		world_rotations = common.transform.world_rotations.roll + data.local_rotations.roll;
	}
	static inline void TransformParticleRelative(tfxParticleData &data, tfxVec2 &world_position, float &world_rotations, const tfxCommon &common, const tfxVec3 &from_position) {
		world_rotations = data.local_rotations.roll;
		float s = sin(data.local_rotations.roll);
		float c = cos(data.local_rotations.roll);
		tfxVec2 rotatevec = mmTransformVector(common.transform.matrix, tfxVec2(data.local_position.x, data.local_position.y) + common.handle.xy());
		world_position = from_position.xy() + rotatevec * common.transform.scale.xy();
	}
	static inline void TransformParticleRelativeLine(tfxParticleData &data, tfxVec2 &world_position, float &world_rotations, const tfxCommon &common, const tfxVec3 &from_position) {
		world_rotations = common.transform.world_rotations.roll + data.local_rotations.roll;
		float s = sin(data.local_rotations.roll);
		float c = cos(data.local_rotations.roll);
		tfxVec2 rotatevec = mmTransformVector(common.transform.matrix, tfxVec2(data.local_position.x, data.local_position.y) + common.handle.xy());
		world_position = from_position.xy() + rotatevec * common.transform.scale.xy();
	}
	static inline void TransformParticle3dPositions(tfxParticleData &data, tfxVec3 &world_position, tfxVec3 &world_rotations, const tfxCommon &common, const tfxVec3 &from_position) {
		world_position = data.local_position;
	}
	static inline void TransformParticle3dPositionsRelative(tfxParticleData &data, tfxVec3 &world_position, tfxVec3 &world_rotations, const tfxCommon &common, const tfxVec3 &from_position) {
		tfxVec4 rotatevec = mmTransformVector(common.transform.matrix, data.local_position + common.handle);
		world_position = common.transform.world_position + rotatevec.xyz();
	}
	static inline void TransformParticle3d(tfxParticleData &data, tfxVec3 &world_position, tfxVec3 &world_rotations, const tfxCommon &common, const tfxVec3 &from_position) {
		world_position = data.local_position;
		world_rotations = data.local_rotations;
	}
	static inline void TransformParticle3dAngle(tfxParticleData &data, tfxVec3 &world_position, tfxVec3 &world_rotations, const tfxCommon &common, const tfxVec3 &from_position) {
		world_position = data.local_position;
		world_rotations = common.transform.world_rotations + data.local_rotations;
	}
	static inline void TransformParticle3dRelative(tfxParticleData &data, tfxVec3 &world_position, tfxVec3 &world_rotations, const tfxCommon &common, const tfxVec3 &from_position) {
		world_rotations = data.local_rotations;
		tfxVec4 rotatevec = mmTransformVector(common.transform.matrix, data.local_position + common.handle);
		world_position = from_position + rotatevec.xyz() * common.transform.scale;
	}
	//todo: redundant function? Remove if so.
	static inline void TransformParticle3dRelativeLine(tfxParticleData &data, tfxVec3 &world_position, tfxVec3 &world_rotations, const tfxCommon &common, const tfxVec3 &from_position) {
		world_rotations = data.local_rotations;
		tfxVec4 rotatevec = mmTransformVector(common.transform.matrix, data.local_position + common.handle);
		world_position = from_position + rotatevec.xyz() * common.transform.scale;
	}
	static inline int SortParticles(void const *left, void const *right) {
		float d1 = static_cast<const tfxSpawnPosition*>(left)->distance_to_camera;
		float d2 = static_cast<const tfxSpawnPosition*>(right)->distance_to_camera;
		return (d2 > d1) - (d2 < d1);
	}

	static inline void InsertionSortParticles(tfxring<tfxParticle> &particles, tfxring<tfxParticle> &current_buffer) {
		tfxPROFILE;
		for (tfxU32 i = 1; i < particles.current_size; ++i) {
			tfxParticle key = particles[i];
			int j = i - 1;
			while (j >= 0 && key.data.depth > particles[j].data.depth) {
				particles[j + 1] = particles[j];
				current_buffer[particles[j + 1].prev_index].next_ptr = &particles[j + 1];
				--j;
			}
			particles[j + 1] = key;
			current_buffer[particles[j + 1].prev_index].next_ptr = &particles[j + 1];
		}
	}

	/*static inline void InsertionSortSprites3d(tfxring<tfxParticleSprite3d> &sprites, tfxvec<tfxring<tfxParticle>> &particle_banks) {
		tfxPROFILE;
		for (tfxU32 i = 1; i < sprites.current_size; ++i) {
			tfxParticleSprite3d key = sprites[i];
			int j = i - 1;
			while (j >= 0 && key.depth > sprites[j].depth) {
				sprites[j + 1] = sprites[j];
				particle_banks[(sprites[j + 1].particle & 0xFFF00000) >> 20][sprites[j + 1].particle & 0x000FFFFF].sprite_index = j + 1;
				--j;
			}
			sprites[j + 1] = key;
			particle_banks[(sprites[j + 1].particle & 0xFFF00000) >> 20][sprites[j + 1].particle & 0x000FFFFF].sprite_index = j + 1;
		}
	}*/

	static inline void InsertionSortParticleFrame(tfxvec<tfxParticleFrame> &particles) {
		for (tfxU32 i = 1; i < particles.current_size; ++i) {
			tfxParticleFrame key = particles[i];
			int j = i - 1;

			while (j >= 0 && key.depth > particles[j].depth) {
				particles[j + 1] = particles[j];
				--j;
			}
			particles[j + 1] = key;
		}
	}
	inline tfxVec3 Tween(float tween, const tfxVec3 &world, const tfxVec3 &captured) {
		tfxVec3 tweened;
		tweened = world * tween + captured * (1.f - tween);
		return tweened;
	}
	inline tfxVec2 Tween2d(float tween, const tfxVec4 &world) {
		tfxVec2 tweened;
		tweened = world.xy() * tween + world.zw() * (1.f - tween);
		return tweened;
	}
	inline tfxVec2 Tween2d(float tween, const tfxVec2 &world, const tfxVec2 &captured) {
		tfxVec2 tweened;
		tweened = world * tween + captured * (1.f - tween);
		return tweened;
	}
	static inline tfxVec3 SetParticleAlignment(tfxParticle &p, tfxVec3 &position, tfxEmitterProperties &properties) {
		if (properties.vector_align_type == tfxVectorAlignType_motion) {
			tfxVec3 alignment_vector = position - p.data.captured_position;
			float l = FastLength(alignment_vector);
			p.data.velocity_normal.w *= l * 10.f;
			return FastNormalizeVec(alignment_vector);
		}
		else if (properties.vector_align_type == tfxVectorAlignType_emission) {
			return p.data.velocity_normal.xyz();
		}
		else if (properties.vector_align_type == tfxVectorAlignType_emitter) {
			return mmTransformVector(p.parent->common.transform.matrix, tfxVec4(0.f, 1.f, 0.f, 0.f)).xyz();
		}
		return tfxVec3(0.f, 0.002f, 0.f);
	}
	float Interpolatef(float tween, float, float);
	int ValidateEffectPackage(const char *filename);
	void ReloadBaseValues(tfxParticle &p, tfxEffectEmitter &e);

	//Particle initialisation functions, one for 2d one for 3d effects
	void InitialiseParticle2d(tfxParticleData &data, tfxSpriteTransform2d &sprite_transform, tfxEmitterState &emitter, tfxCommon &common, tfxEmitterSpawnControls &spawn_values, tfxEffectEmitter *library_link, float tween);
	tfxSpawnPosition InitialisePosition3d(tfxEmitterState &current, tfxCommon &common, tfxEmitterSpawnControls &spawn_values, tfxEffectEmitter *library_link, float tween);
	void InitialiseParticle3d(tfxParticleData &data, tfxSpriteTransform3d &sprite_transform, tfxEmitterState &current, tfxCommon &common, tfxEmitterSpawnControls &spawn_values, tfxEffectEmitter *library_link, float tween);
	void UpdateParticle2d(tfxParticleData &data, tfxVec2 &sprite_scale, tfxControlData &c);
	void UpdateParticle3d(tfxParticleData &data, tfxVec2 &sprite_scale, tfxControlData &c);

	//Get a graph by tfxGraphID
	tfxGraph &GetGraph(tfxEffectLibrary &library, tfxGraphID &graph_id);

	//Set the udpate frequency for all particle effects - There may be options in the future for individual effects to be updated at their own specific frequency.
	void SetUpdateFrequency(float fps);

	inline float GetUpdateFrequency() { return tfxUPDATE_FREQUENCY; }
	inline float GetUpdateTime() { return tfxUPDATE_TIME; }
	inline float GetFrameLength() { return tfxFRAME_LENGTH; }
	inline void SetLookUpFrequency(float frequency) {
		tfxLOOKUP_FREQUENCY = frequency;
	}
	inline void SetLookUpFrequencyOvertime(float frequency) {
		tfxLOOKUP_FREQUENCY_OVERTIME = frequency;
	}
	int GetShapesInPackage(const char *filename);
	int GetEffectLibraryStats(const char *filename, tfxEffectLibraryStats &stats);
	tfxEffectLibraryStats CreateLibraryStats(tfxEffectLibrary &lib);
	tfxErrorFlags LoadEffectLibraryPackage(const char *filename, tfxEffectLibrary &lib, void(*shape_loader)(const char *filename, tfxImageData &image_data, void *raw_image_data, int image_size, void *user_data) = nullptr, void *user_data = nullptr, bool read_only = true);

	//---
	//Prepare an effect template for setting up function call backs to customise the behaviour of the effect in realtime
	//Returns true on success.
	bool PrepareEffectTemplate(tfxEffectLibrary &library, const char *name, tfxEffectTemplate &effect_template);
}

