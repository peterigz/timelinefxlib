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
	struct tfxParticle;
	struct tfxParticleData;
	struct tfxComputeSprite;
	struct tfxComputeParticle;
	struct tfxAnimationSettings;
	struct tfxEffectLibrary;
	struct tfxText;
	struct tfxEffectTemplate;

//--------------------------------------------------------------
//macros
#define TFX_VERSION "Alpha"
#define TFX_VERSION_NUMBER 3.29.2022

#define tfxMAX_FRAME 20000.f
#define tfxNullParent 0xFFFFFFFF
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

//type defs
typedef unsigned int tfxU32;
typedef int tfxS32;
typedef unsigned long long tfxU64;
typedef long long tfxS64;
typedef tfxU32 tfxEffectID;

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
		tfxEffectManagerFlags_update_base_values = 1 << 6
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
		tfxEmitterStateFlags_is_sub_emitter = 1 << 28,
	};

	enum tfxVectorFieldFlags_: unsigned char {
		tfxVectorFieldFlags_none = 0,
		tfxVectorFieldFlags_repeat_horizontal = 1 << 0,						//Field will repeat horizontally
		tfxVectorFieldFlags_repeat_vertical = 1 << 1						//Field will repeat vertically
	};

	enum tfxPackageErrorCode : unsigned int {
		tfxPackageErrorCode_unable_to_open_file = 1,
		tfxPackageErrorCode_unable_to_read_file,
		tfxPackageErrorCode_wrong_file_size,
		tfxPackageErrorCode_invalid_format,
		tfxPackageErrorCode_no_inventory,
		tfxPackageErrorCode_invalid_inventory,
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

	//Credit to ocornut https://github.com/ocornut/imgui/commits?author=ocornut
	//std::vector replacement with some extra stuff and tweaks specific to Qulkan/TimelineFX
	template<typename T>
	struct tfxvec {
		tfxU32 current_size;
		tfxU32 capacity;
		T* data;

		// Provide standard typedefs but we don't use them ourselves.
		typedef T                   value_type;
		typedef value_type*         iterator;
		typedef const value_type*   const_iterator;

		inline tfxvec() { current_size = capacity = 0; data = NULL; }
		inline tfxvec(tfxU32 qty) { current_size = capacity = 0; data = NULL; resize(qty); }
		inline tfxvec(T* from, T* to) { current_size = capacity = 0; data = NULL; auto current = from; while (current != to + 1) { push_back(*current); ++current; } }
		inline tfxvec(std::initializer_list<T> t) { current_size = capacity = 0; data = NULL; for (T element : t) { push_back(element); } }
		inline tfxvec(const tfxvec<T> &src) { current_size = capacity = 0; data = NULL; resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); }
		inline tfxvec<T>& operator=(const tfxvec<T>& src) { clear(); resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); return *this; }
		inline ~tfxvec() { if (data) free(data); data = NULL; current_size = capacity = 0; }

		inline bool			empty() { return current_size == 0; }
		inline tfxU32			size() { return current_size; }
		inline const tfxU32	size() const { return current_size; }
		inline T&           operator[](tfxU32 i) { return data[i]; }
		inline const T&     operator[](tfxU32 i) const { assert(i < current_size); return data[i]; }

		inline void         free_all() { if (data) { current_size = capacity = 0; free(data); data = NULL; } }
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
		inline void         reserve(tfxU32 new_capacity) { if (new_capacity <= capacity) return; T* new_data = (T*)malloc((size_t)new_capacity * sizeof(T)); if (data) { memcpy(new_data, data, (size_t)current_size * sizeof(T)); free(data); } data = new_data; capacity = new_capacity; }

		inline T&	        push_back(const T& v) { if (current_size == capacity) reserve(_grow_capacity(current_size + 1)); new((void*)(data + current_size)) T(v); current_size++; return data[current_size - 1]; }
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
		inline tfxU32          index_from_ptr(const T* it) const { assert(it >= data && it < data + current_size); const ptrdiff_t off = it - data; return (tfxU32)off; }

		inline void			create_pool(tfxU32 amount) { assert(current_size == 0); T base; reserve(amount); for (tfxU32 i = 0; i != capacity; ++i) { new((void*)(data + current_size)) T(base); current_size++; } }
		inline void			create_pool_with(tfxU32 amount, const T &base) { assert(current_size == 0);  reserve(amount); for (tfxU32 i = 0; i != capacity; ++i) { new((void*)(data + current_size)) T(base); current_size++; } }

	};

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
		inline tfxText ReadLine();
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

	template <typename T>
	struct tfxMemoryRange {
		T *data = NULL;
		tfxS32 next_index = -1;
		tfxU32 capacity = 0;
		tfxU32 current_size = 0;

		inline T*           begin() { return data; }
		inline const T*     begin() const { return data; }
		inline T*           end() { return data + current_size; }
		inline const T*     end() const { return data + current_size; }
		inline T&           back() { return data[current_size - 1]; }
		inline const T&     back() const { return data[current_size - 1]; }
	};

	template <typename T>
	struct tfxRangeAllocator {
		T *data;					//big allocation for splitting up into smaller ranges
		T *end_of_allocated;
		tfxU32 memory_remaining;
		tfxU32 total_memory;
		tfxU32 size_diff_threshold;
		tfxU32 current_buffer;

		tfxvec<tfxMemoryRange<T>> ranges;
		tfxvec<tfxU32> free_ranges;

		tfxRangeAllocator() { data = end_of_allocated = NULL; memory_remaining = size_diff_threshold = current_buffer = 0; }

		inline void FreeAll() {
			if (data) free(data);
			data = end_of_allocated = NULL;
			memory_remaining = 0;
			ranges.free_all();
			free_ranges.free_all();
		}

		inline void Reset() {
			memory_remaining = total_memory;
			ranges.clear();
			free_ranges.clear();
		}

		inline bool Allocate(tfxU32 size, tfxS32 &range_index) {
			if (size == 0) return false;
			if (free_ranges.size() > 0) {
				tfxU32 size_diff = tfxMAX_UINT;
				tfxS32 best_fit_index = -1;
				tfxS32 index = 0;
				for (auto range_index : free_ranges) {
					if (size_diff > ranges[range_index].capacity && ranges[range_index].capacity > size) {
						size_diff = ranges[range_index].capacity;
						best_fit_index = index;
					}
					index++;
				}
				if (size_diff != tfxMAX_UINT && (size_diff <= size_diff_threshold || size_diff > memory_remaining)) {
					free_ranges[best_fit_index] = free_ranges.back();
					free_ranges.pop();
					range_index = best_fit_index;
					return true;
				}
			}
			else if (memory_remaining >= size) {
				tfxMemoryRange<T> range;
				range.data = end_of_allocated;
				range.capacity = size;
				ranges.push_back(range);
				end_of_allocated = end_of_allocated + size;
				memory_remaining -= size;
				range_index = ranges.size() - 1;
				return true;
			}
			return false;
		}

		inline tfxU32 AllocateMore(tfxU32 size, tfxS32 range_index) {
			tfxS32 new_range;
			tfxU32 allocatable_size = size > memory_remaining ? memory_remaining : size;
			if (Allocate(allocatable_size, new_range)) {
				ranges[FindLastRangeIndex(range_index)].next_index = new_range;
				return allocatable_size;
			} 
			return 0;
		}

		inline tfxMemoryRange<T> &FindLastRange(tfxS32 starting_range) {
			assert(starting_range >= 0 && (tfxU32)starting_range < ranges.size());
			tfxS32 found_range = starting_range;
			while (ranges[found_range].next_index != -1) {
				found_range = ranges[found_range].next_index;
			}
			return ranges[found_range];
		}

		inline tfxU32 FindLastRangeIndex(tfxS32 starting_range) {
			assert(starting_range >= 0 && (tfxU32)starting_range < ranges.size());
			tfxS32 found_range = starting_range;
			while (ranges[found_range].next_index != -1) {
				found_range = ranges[found_range].next_index;
			}
			return found_range;
		}

		inline T *FindLastRangeBegin(tfxS32 starting_range) {
			assert(starting_range >= 0 && (tfxU32)starting_range < ranges.size());
			tfxS32 found_range = starting_range;
			while (ranges[found_range].next_index != -1) {
				found_range = ranges[found_range].next_index;
			}
			return ranges[found_range].data;
		}

		inline T &FindValueByIndex(tfxU32 i, tfxS32 range_index) {
			while (ranges[range_index].next_index != -1) {
				i -= ranges[range_index].current_size;
				range_index = ranges[range_index].next_index;
			}
			return ranges[range_index].data[i];
		}

		inline T &Back(tfxS32 range_index) {
			while (ranges[range_index].next_index != -1) {
				range_index = ranges[range_index].next_index;
			}
			return ranges[range_index].back();
		}

		inline T &End(tfxS32 range_index) {
			while (ranges[range_index].next_index != -1) {
				range_index = ranges[range_index].next_index;
			}
			return ranges[range_index].end();
		}

		inline void FreeRanges(tfxS32 range_index) {
			free_ranges.push_back(range_index);
			while (ranges[range_index].next_index != -1) {
				range_index = ranges[range_index].next_index;
				free_ranges.push_back(range_index);
			}
		}

		inline void ClearRanges(tfxS32 range_index) {
			free_ranges.push_back(range_index);
			while (ranges[range_index].next_index != -1) {
				ranges[range_index].current_size = 0;
				range_index = ranges[range_index].next_index;
			}
		}

	};

	template <typename T>
	static inline tfxRangeAllocator<T> CreateRangeAllocator(tfxU32 size, tfxU32 size_diff_threshold = 8) {
		tfxRangeAllocator<T> allocator;
		T* new_data = (T*)malloc(size * sizeof(T));
		allocator.data = new_data;
		allocator.end_of_allocated = allocator.data;
		allocator.memory_remaining = size;
		allocator.size_diff_threshold = size_diff_threshold;
		allocator.total_memory = size;
		return allocator;
	}

	template <typename T>
	struct tfxRangeVec {
		tfxRangeAllocator<T> *allocator;	//Pointer to the allocator that manages the memory and ranges of that memory
		tfxS32 range_index;				//The first range in the allocator
		tfxU32 current_size;			//Current size of the range vec. This will be the total of all ranges if there are more then one
		tfxU32 capacity;				//The total capacity of the range vec
		tfxU32 current_range_index;		//the current range index for iterating all ranges

		tfxRangeVec() : allocator(NULL) { range_index = -1; current_size = capacity = current_range_index = 0; }
		tfxRangeVec(tfxRangeAllocator<T> *allocator_init) : allocator(allocator_init) { range_index = -1; current_size = capacity = current_range_index = 0; }

		inline bool			empty() { return current_size == 0; }
		inline tfxU32		size() { return current_size; }
		inline const tfxU32	size() const { return current_size; }
		inline T&           operator[](tfxU32 i) { return allocator->FindValueByIndex(i, range_index); }
		inline const T&     operator[](tfxU32 i) const { assert(i < current_size); return allocator->FindValueByIndex(); }

		inline void         free_all() { if (range_index > -1) { current_size = capacity = 0; allocator->FreeRanges(range_index); range_index = -1; } }
		inline T*           begin() { return allocator->ranges[range_index].data; }
		inline const T*     begin() const { return allocator->ranges[range_index].data; }
		inline T*           end() { return &allocator->End(range_index); }
		inline const T*     end() const { return &allocator->End(range_index); }
		inline T&           front() { assert(current_size > 0); return allocator->ranges[range_index].data[0]; }
		inline const T&     front() const { assert(current_size > 0); return allocator->ranges[range_index].data[0]; }
		inline T&           back() { assert(current_size > 0); return allocator->Back(range_index); }
		inline const T&     back() const { assert(current_size > 0); return allocator->Back(range_index); }
		inline void         clear() { 
			current_size = 0; 
			allocator->ClearRanges(range_index);
		}

		inline tfxU32       _grow_capacity(tfxU32 sz) const { tfxU32 new_capacity = capacity ? (capacity + capacity / 2) : 8; return new_capacity > sz ? new_capacity : sz; }
		inline bool			reserve(tfxU32 size) {
			assert(allocator);	//Must assign and allocator before doing anything with a tfxRangeVec
			if (capacity == 0) {
				if (!allocator->Allocate(size, range_index)) {
					return false;
				}
			}
			else {
				tfxU32 reserved_amount = allocator->AllocateMore(size, range_index);
				if(!reserved_amount) {
					return false;
				}
				capacity += reserved_amount;
			}
			return true;
		}
		inline bool	        push_back(const T& v) { 
			assert(allocator);	//Must assign and allocator before doing anything with a tfxRangeVec
			if (capacity == 0) {
				capacity = _grow_capacity(current_size + 1);
				if (!allocator->Allocate(capacity, range_index)) {
					return false;
				}
				current_range_index = range_index;
			}
			else if (current_size == capacity) {
				tfxU32 new_capacity = _grow_capacity(current_size + 1);
				tfxU32 reserved_amount = allocator->AllocateMore(new_capacity, range_index);
				if (!reserved_amount) {
					return false;
				}
				capacity += reserved_amount;
			}
			tfxMemoryRange<T> &last_range = allocator->FindLastRange(range_index);
			new((void*)(last_range.data + last_range.current_size++)) T(v); 
			current_size++; 
			return true;
		}
		inline void ResetIteratorIndex() { 
			current_range_index = range_index; 
		}
		inline bool EndOfRanges() { 
			current_range_index = allocator->ranges[current_range_index].next_index;
			if (current_range_index == -1) { 
				current_range_index = range_index; 
				return true; 
			} 
			return false;
		}
		inline tfxMemoryRange<T> &NextRange() {
			assert(allocator);	//Must assign and allocator before doing anything with a tfxRangeVec
			return allocator->ranges[current_range_index];
		}
	};

	template <typename T>
	static inline bool MergeRangeVecs(tfxRangeVec<T> &range_vec_from, tfxRangeVec<T> &range_vec_to) {
		range_vec_to.current_size = 0;
		range_vec_to.capacity = range_vec_from.current_size;
		if (range_vec_to.capacity == 0) {
			range_vec_to.range_index = -1;
			return false;
		}
		if (!range_vec_to.allocator->Allocate(range_vec_to.capacity, range_vec_to.range_index)) {
			return false;
		}
		range_vec_to.current_range_index = range_vec_to.range_index;
		return true;
	}

	//Just the very basic vector types that we need
	struct tfxVec2 {
		float x, y;

		tfxVec2() { x = y = 0.f; }
		tfxVec2(float _x, float _y) : x(_x), y(_y) {}

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
		inline tfxVec3 xyz() { return tfxVec3(x, y, z); }

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
	struct tfxText {
		tfxvec<char> string;

		tfxText() {}
		tfxText(const char *text) { size_t length = strnlen_s(text, 512); if (!length) { Clear(); return; }; if (string.capacity < length) string.reserve((tfxU32)length); memcpy(string.data, text, length); string.current_size = (tfxU32)length; NullTerminate(); }
		inline char operator[](tfxU32 i) const { assert(i < string.current_size); return string[i]; }
		inline char operator[](tfxU32 i) { assert(i < string.current_size); return string[i]; }
		inline void operator=(const char *text) { size_t length = strnlen_s(text, 512); if (!length) { Clear(); return; }; if (string.capacity < length) string.reserve((tfxU32)length); memcpy(string.data, text, length); string.current_size = (tfxU32)length; NullTerminate(); }
		inline void operator=(const tfxText &label) { string = label.string; }
		inline bool operator==(const char *string) { return !strcmp(string, c_str()); }
		inline bool operator==(const tfxText string) { return !strcmp(c_str(), string.c_str()); }
		inline bool operator!=(const char *string) { return strcmp(string, c_str()); }
		inline bool operator!=(const tfxText string) { return strcmp(c_str(), string.c_str()); }
		const char* begin() const { return string.data ? &string.front() : NULL; }
		const char* end() const { return string.data ? &string.back() : NULL; }  
		inline const char *c_str() const { return string.current_size ? string.data : ""; }
		inline void Clear() { string.clear(); }
		int Find(const char *needle);
		tfxText Lower();
		inline tfxU32 Length() const { return string.current_size ? string.current_size - 1 : 0; }
		void AddLine(const char *format, ...);
		void Appendf(const char *format, ...);
		void Appendv(const char *format, va_list args);
		inline void Append(char c) { 
			if (string.current_size) {
				string.pop();
			}
			string.push_back(c); 
			NullTerminate();
		}
		inline void Pop() {
			if (!Length()) return;
			if(string.back() == NULL)
				string.pop();
			string.pop();
			NullTerminate();
		}
		inline void Trim(char c = 32) {
			if (!Length()) return;
			if(string.back() == NULL)
				string.pop();
			while (string.back() == c && string.current_size) {
				string.pop();
			}
			NullTerminate();
		}
		void NullTerminate() { string.push_back(NULL); }
		bool SaveToFile(const char *file_name);
		const bool IsInt() const;
		const bool IsFloat() const;
	};

	typedef unsigned long long tfxKey;

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

		tfxStorageMap() {}

		//Insert a new T value into the storage
		inline void Insert(const char *name, const T &value) {
			tfxKey key = tfxXXHash64::hash(name, strlen(name), 0);
			SetIndex(key, value);
		}

		//Insert a new T value into the storage
		inline void Insert(tfxText name, const T &value) {
			tfxKey key = tfxXXHash64::hash(name.c_str(), name.Length(), 0);
			SetIndex(key, value);
		}

		//Insert a new T value into the storage
		inline void Insert(tfxKey key, const T &value) {
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

		inline bool ValidName(const tfxText &name) {
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

		inline T &At(const tfxText &name) {
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

		int GetIndex(const tfxText &name) {
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

		tfxProfileTag(tfxU32 id, const char *name) {
			profile = tfxPROFILE_ARRAY + id;
			profile->name = name;
			profile->run_time -= Microsecs();
			profile->cycle_count -= __rdtsc();
			profile->hit_count++;
		}

		~tfxProfileTag() {
			profile->run_time += Microsecs();
			profile->cycle_count += __rdtsc();
		}

	};

#define tfxPROFILE tfxProfileTag tfx_tag((tfxU32)__COUNTER__, __FUNCTION__)

	bool EndOfProfiles();
	tfxProfile* NextProfile();

#define tfxINVALID 0xFFFFFFFF

	struct tfxrange {
		tfxU32 capacity;
		tfxU32 offset_into_memory;
	};

	struct tfxmemory {
		tfxU32 current_size;
		tfxU32 capacity;
		void* data;
		tfxvec<tfxrange> ranges;
		tfxStorageMap<tfxvec<tfxU32>> free_ranges;

		inline tfxmemory() { current_size = capacity = 0; data = NULL; }
		inline tfxmemory(tfxU32 qty) { current_size = capacity = 0; data = NULL; resize(qty); }
		inline ~tfxmemory() { if (data) free(data); data = NULL; current_size = capacity = 0; }

		inline bool					empty() { return current_size == 0; }
		inline tfxU32			size() { return current_size; }
		inline const tfxU32	size() const { return current_size; }

		inline void					free_all() { if (data) { current_size = capacity = 0; free(data); data = NULL; } }
		inline void					clear() { if (data) { current_size = 0; ranges.clear(); free_ranges.Clear(); } }

		inline tfxU32			_grow_capacity(tfxU32 sz) const { tfxU32 new_capacity = capacity ? (capacity + capacity / 2) : 8; return new_capacity > sz ? new_capacity : sz; }
		inline void					resize(tfxU32 new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
		inline void					shrink(tfxU32 new_size) { assert(new_size <= current_size); current_size = new_size; }
		inline void					reserve(tfxU32 new_capacity) { if (new_capacity <= capacity) return; char* new_data = (char*)malloc((size_t)new_capacity * sizeof(char)); if (data) { memcpy(new_data, data, (size_t)current_size * sizeof(char)); free(data); } data = new_data; capacity = new_capacity; }
		inline tfxU32			free_unused_space() { return capacity - current_size; }

		inline bool has_free_range_available(tfxU32 bytes) {
			return free_ranges.ValidKey(tfxKey(bytes));
		}
		inline tfxU32			get_range(tfxU32 bytes) { 
			if (current_size + bytes >= capacity) reserve(_grow_capacity(current_size + bytes)); 
			tfxrange range;
			range.offset_into_memory = current_size;
			range.capacity = bytes;
			if (free_ranges.ValidKey(tfxKey(bytes))) {
				tfxvec<tfxU32> &free = free_ranges.At(tfxKey(bytes));
				if (!free.empty()) {
					return free.pop_back();
				}
			}
			ranges.push_back(range);
			current_size += bytes;
			return ranges.current_size - 1;
		}
		inline void					free_range(tfxU32 index) {
			assert(index < ranges.size());		//Index to free must be within size of ranges list
			int capacity = ranges[index].capacity;
			if (!free_ranges.ValidKey(tfxKey(ranges[index].capacity))) {
				tfxvec<tfxU32> new_free_ranges;
				free_ranges.Insert((tfxKey)ranges[index].capacity, new_free_ranges);
				free_ranges.At((tfxKey)ranges[index].capacity).push_back(index);
			}
			else {
				free_ranges.At((tfxKey)ranges[index].capacity).push_back(index);
			}
		}

	};

	//Fixed size ring buffer, you must reserve the amount of memory you need ahead of time with reserve or just use the constructor specifying the number of elements you need
	//Calling reserve on a buffer that already contains data will erase that data.
	//Iterate over the ring with:
	/*
	ring.reset();
	while(!ring.eob()) {
		auto &item = ring.next();
	}
	*/
	//You must call free_all when done with the buffer, there's no deconstructor
	//You can also use this in combination with tfxmemory, assign_memory to grab a range of memory from there instead of allocating some. In this case do not call free_all or you will delete
	//the memory in tfxmem instead! Basically this is only for use in TimelineFX internally.
	template<typename T>
	struct tfxring {
		tfxU32 current_size;
		tfxU32 capacity;
		tfxU32 start_index;
		tfxU32 pos;
		tfxU32 range_index;
		T* data;
		inline tfxring() { pos = start_index = current_size = capacity = 0; data = NULL; range_index = tfxINVALID; }
		inline tfxring(tfxU32 qty) { pos = start_index = current_size = capacity = 0; data = NULL; reserve(qty); }
		inline tfxring(tfxmemory &mem, tfxU32 index, tfxU32 size) { pos = start_index = current_size = 0; capacity = size; data = mem.data + mem.ranges[index].offset_into_memory; range_index = index; }

		inline bool			empty() { return current_size == 0; }
		inline bool			full() { return current_size == capacity; }
		inline void         free_all() { assert(range_index == tfxINVALID); /*Call free_range instead if using tfxmemory*/ if (data) { pos = start_index = current_size = capacity = 0; free(data); data = NULL; } }
		inline void         free_range(tfxmemory &mem) { 
		if (range_index == tfxINVALID) return; 
			if (data) { 
				pos = start_index = current_size = capacity = 0; 
				mem.free_range(range_index); 
				range_index = tfxINVALID; 
				data = NULL; 
			} 
		}
		inline void			prep_for_new() { start_index = current_size = pos = 0; data = NULL; data = (T*)malloc((size_t)capacity * sizeof(T)); }
		inline void         clear() { if (data) { start_index = pos = current_size = 0; } }
		inline tfxU32			size() { return current_size; }
		inline const tfxU32	size() const { return current_size; }
		inline T&           operator[](tfxU32 i) { return data[(i + start_index) % capacity]; }
		inline const T&     operator[](tfxU32 i) const { return data[(i + start_index) % capacity]; }
		inline T&           AtAbs(tfxU32 i) { return data[i]; }
		inline const T&     AtAbs(tfxU32 i) const { return data[i]; }

		inline tfxU32 end_index() { return (start_index + current_size) % capacity; }
		inline tfxU32 before_start_index() { return start_index == 0 ? capacity - 1 : start_index - 1; }
		inline tfxU32 last_index() { return (start_index + current_size - 1) % capacity; }
		inline T*			eob_ptr() { return data + (start_index + current_size) % capacity; }
		inline T&           front() { assert(current_size > 0); return data[start_index]; }
		inline const T&     front() const { assert(current_size > 0); return data[start_index]; }
		inline T&           back() { assert(current_size > 0); return data[last_index()]; }
		inline const T&     back() const { assert(current_size > 0); return data[last_index()]; }

		inline void         pop() { assert(current_size > 0); current_size--; }
		inline T&	        pop_back() { assert(current_size > 0); current_size--; return data[current_size]; }
		inline bool	        push_back(const T& v) { if (current_size == capacity) return false; new((void*)(data + end_index())) T(v); current_size++; return true; }
		inline bool	        push_front(const T& v) { if (current_size == capacity) return false; new((void*)(data + before_start_index())) T(v); current_size++; start_index = before_start_index(); return true; }
		inline void         reserve(tfxU32 new_capacity) { if (new_capacity <= capacity) return; free(data); data = (T*)malloc((size_t)new_capacity * sizeof(T)); capacity = new_capacity; start_index = current_size = 0; }
		inline void			assign_memory(tfxmemory &mem, tfxU32 element_size, tfxU32 element_count) { 
			assert(range_index == tfxINVALID);	//call refresh instead if a range has already been assigned, or free_all and then assign another range;
			if (!element_count) return;
			int index = mem.get_range(element_size * element_count);
			pos = start_index = current_size = 0; range_index = index; capacity = element_count; void *ptr = (char*)mem.data + mem.ranges[index].offset_into_memory; data = static_cast<T*>(ptr); }
		inline void			refresh(tfxmemory &mem) { void *ptr = (char*)mem.data + mem.ranges[range_index].offset_into_memory; data = static_cast<T*>(ptr); }
		inline void			reset_to_null() { pos = start_index = current_size = 0; range_index = tfxINVALID; data = NULL; }
		inline void			bump() { if (current_size == 0) return; start_index++; start_index %= capacity; current_size--; }
		inline void			bump(tfxU32 amount) { if (current_size == 0) return; if (amount > current_size) amount = current_size; start_index += amount; start_index %= capacity; current_size -= amount; }
		inline void			shrink(tfxU32 amount) { if (amount > current_size) current_size = 0; else current_size -= amount; }

		inline bool			reset() { if (current_size == 0) return false; assert(data); pos = 0; return true; }
		inline T&			next() { assert(current_size > 0); tfxU32 current_pos = (start_index + pos) % capacity; pos++; return data[current_pos]; }
		inline T&			next(int start) { assert(current_size > 0); tfxU32 current_pos = (start + pos) % capacity; pos++; return data[current_pos]; }
		inline bool			eob() { return pos >= current_size; }
	};

	const tfxU32 tfxMAGIC_NUMBER = '!XFT';
	const tfxU32 tfxMAGIC_NUMBER_INVENTORY = '!VNI';
	const tfxU32 tfxFILE_VERSION = 1;	//Not doing anything with this yet

	//Basic package manager used for reading/writing effects files
	struct tfxHeader {
		tfxU32 magic_number;						//Magic number to confirm file format
		tfxU32 file_version;						//The version of the file
		tfxU32 flags;								//Any flags for the file
		tfxU32 reserved0;							//Reserved for future if needed
		tfxU64 offset_to_inventory;				//Memory offset for the inventory of files
		tfxU64 reserved1;							//More reserved space
		tfxU64 reserved2;							//More reserved space
		tfxU64 reserved3;							//More reserved space
		tfxU64 reserved4;							//More reserved space
		tfxU64 reserved5;							//More reserved space
	};

	struct tfxEntryInfo {
		tfxText file_name;						//The file name of the name stored in the package
		tfxU64 offset_from_start_of_file;			//Offset from the start of the file to where the file is located
		tfxU64 file_size;							//The size of the file
		tfxstream data;							//The file data
		
		void FreeData();
	};

	struct tfxInventory {
		tfxU32 magic_number;						//Magic number to confirm format of the Inventory
		tfxU32 entry_count;						//Number of files in the inventory
		tfxStorageMap<tfxEntryInfo> entries;	//The inventory list
	};

	struct tfxPackage {
		tfxText file_path;
		tfxHeader header;
		tfxInventory inventory;
		tfxU64 file_size;							//The total file size of the package, should match file size on disk
		tfxstream file_data;					//Dump of the data from the package file on disk

		~tfxPackage();

		tfxEntryInfo *GetFile(const char *name);
		void AddFile(tfxEntryInfo file);
		void AddFile(const char *file_name, tfxstream &data);
		void Free();

	};
	
	tfxstream ReadEntireFile(const char *file_name, bool terminate = false);
	int LoadPackage(const char *file_name, tfxPackage &package);
	int LoadPackage(tfxstream &stream, tfxPackage &package);
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
		tfxvec<float> values;
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
		tfxvec<tfxAttributeNode> nodes;
		tfxGraphLookup lookup;
		tfxU32 index;

		tfxGraph();
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
		tfxvec<tfxAttributeNode>& Nodes();
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

	};

	tfxVec4 GetMinMaxGraphValues(tfxGraphPreset preset);

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
		bool attach_effect_to_camera;
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
		//Pointer to the ImageData in the EffectLibary. 
		tfxImageData *image;
		//Assigns either alpha or additive blend to particles that are spawned
		tfxEmissionType emission_type;
		//Should particles emit towards the center of the emitter or away, or in a specific direction
		tfxEmissionDirection emission_direction;
		//How particles should behave when they reach the end of the line
		tfxLineTraversalEndBehaviour end_behaviour;
		//The rotation of particles when they spawn, or behave overtime if tfxAlign is used
		tfxAngleSettingFlags angle_settings;
		//For 3d effects, the type of billboarding: 0 = use billboarding (always face camera), 1 = No billboarding, 2 = No billboarding and align with motion
		tfxBillboardingOptions billboard_option;
		//When aligning the billboard along a vector, you can set the type of vector that it aligns with
		tfxVectorAlignType vector_align_type;

		//Bit field of various boolean flags
		tfxParticleControlFlags compute_flags;

		//Offset to draw particles at
		tfxVec2 image_handle;
		//Animation frame rate
		float frame_rate;
		//Offset of emitters
		tfxVec3 emitter_handle;
		//When single flag is set, spawn this amount of particles in one go
		tfxU32 spawn_amount;
		//If single shot flag is set then you can limit how many times it will loop over it's overtime graphs before expiring
		tfxU32 single_shot_limit;
		//Layer of the particle manager that the particle is added to
		tfxU32 layer;
		//The shape being used for all particles spawned from the emitter
		tfxU32 shape_index;

		//Angle added to the rotation of the particle when spawned or random angle range if angle setting is set to tfxRandom
		tfxVec3 angle_offsets;
		//The number of rows/columns/ellipse/line points in the grid when spawn on grid flag is used
		tfxVec3 grid_points;
		//The number of millisecs before an effect or emitter will loop back round to the beginning of it's graph lookups
		float loop_length;
		//The start frame index of the animation
		float start_frame;
		//The final frame index of the animation
		float end_frame;
		//Milliseconds to delay spawing
		float delay_spawning;

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
			timeout(1000.f),
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
		float qty_step_size;
		//The callback to transform the particles each update. This will change based on the properties of the emitter
		void(*transform_particle_callback)(tfxParticleData &data, const tfxCommon &common, const tfxVec3 &from_position);

		tfxEmitterState() :
			amount_remainder(0.f),
			emission_direction_normal(0.f, 1.f, 0.f),
			qty_step_size(0.f)
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

	float GetEmissionDirection2d(tfxCommon &common, tfxEmitterState &current, tfxEffectEmitter *library_link, tfxVec2 local_position, tfxVec2 world_position, tfxVec2 emitter_size);
	tfxVec3 GetEmissionDirection3d(tfxCommon &common, tfxEmitterState &current, tfxEffectEmitter *library_link, float emission_pitch, float emission_yaw, tfxVec3 local_position, tfxVec3 world_position, tfxVec3 emitter_size);

	struct tfxEffectEmitterInfo {
		//Not required for frame by frame updating - should be moved into an info lookup in library
		//Name of the effect
		tfxText name;						//Todo: Do we need this here?
		//A hash of the directory path to the effect ie Flare/spark
		tfxKey path_hash;
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
			max_sub_emitter_life(0.f)
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
		//The current highest particle age. When using a compute buffer we don't have any reliable way of keeping track of particle counts of individual emitters, so how do we know when to remove an emitter
		//after all it's particles have expired? We set this variable to the highest particle age each time it spawns a particle and then counts it down each frame. When it's 0 then we know that there are no
		//more particles being controlled by this emitter and can therefore time it out.
		float highest_particle_age;
		//compute slot id if a compute shader is being used. Only applied to bottom emitters (emitters with no child effects)
		tfxU32 compute_slot_id;
		//All graphs that the effect uses to lookup attribute values are stored in the library. These variables here are indexes to the array where they're stored
		tfxU32 global;
		tfxU32 property;
		tfxU32 base;
		tfxU32 variation;
		tfxU32 overtime;
		//Pointer to the immediate parent
		tfxEffectEmitter *parent;
		//Pointer to the next pointer in the particle manager buffer. 
		tfxEffectEmitter *next_ptr;
		//Pointer to the sub effect's particle that spawned it
		tfxParticle *parent_particle;
		//Container for all the particles that this emitter will spawn. The memory allocator for this is in the ParticleManager
		tfxRangeVec<tfxParticle> particles;
		//State flags for emitters and effects
		tfxEmitterStateFlags flags;
		tfxEffectPropertyFlags effect_flags;
		//When not using insert sort to guarantee particle order, sort passes offers a more lax way of ordering particles over a number of frames.
		//The more passes the more quickly ordered the particles will be but at a higher cost
		tfxU32 sort_passes;
		//Custom user data, can be accessed in callback functions
		void *user_data;

		tfxU32 info_index;
		tfxU32 property_index;

		//Custom fuction pointers that you can use to override attributes and affect the effect/emitters behaviour in realtime
		//See tfxEffectTemplate for applying these callbacks

		tfxEffectEmitter() :
			highest_particle_age(0),
			parent(nullptr),
			parent_particle(nullptr),
			user_data(nullptr),
			flags(tfxEmitterStateFlags_no_tween_this_update | tfxEmitterStateFlags_enabled),
			effect_flags(tfxEffectPropertyFlags_none),
			sort_passes(1)
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

		void ResetGlobalGraphs(bool add_node = true);
		void ResetBaseGraphs(bool add_node = true);
		void ResetPropertyGraphs(bool add_node = true);
		void ResetVariationGraphs(bool add_node = true);
		void ResetOvertimeGraphs(bool add_node = true);
		void ResetEffectGraphs(bool add_node = true);
		void ResetEmitterGraphs(bool add_node = true);
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
		void Clone(tfxEffectEmitter &clone, tfxEffectEmitter *root_parent, tfxEffectLibrary *destination_library, bool keep_user_data = false, bool force_clone_global = false);
		void EnableAllEmitters();
		void EnableEmitter();
		void DisableAllEmitters();
		void DisableAllEmittersExcept(tfxEffectEmitter &emitter);
		bool IsFiniteEffect();
		void FlagAs3D(bool flag);
		bool Is3DEffect();

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
		float emitter_handle_y;
		float overal_scale;
		float angle_offset;
		tfxOvertimeAttributes *graphs;
		tfxVec2 image_handle;
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

	struct tfxParticleVelocityData {
		tfxVec3 current_velocity;
		float base_velocity;
		tfxVec4 velocity_normal;
		float weight_acceleration;
	};

	//Initial particle struct, looking to optimise this and make as small as possible
	//These are spawned by effector emitter types
	//Particles are stored in the particle manager particle buffer.
	struct tfxParticleData {
		tfxVec3 local_position;			//The local position of the particle, relative to the emitter.
		tfxVec3 world_position;			//The world position of the particle relative to the world/screen.
		tfxVec3 captured_position;		//The captured world coords for interpolating frames
		tfxVec3 local_rotations;
		tfxVec3 world_rotations;
		tfxVec2 scale;
		//Read only when ControlParticle is called, only written to at spawn time
		tfxBase base;							//Base values created when the particle is spawned. They can be different per particle due to variations
		tfxVec4 velocity_normal;			//Current velocity direction, with stretch factor in w
		tfxVec3 alignment_vector;			//Current Vector of alignment
		float depth;						//Distance to camera and can also be used for other things too, but generally a distance to something. Not necessarily actual distance, just a relative number (to avoid needing to square root)
		float noise_offset;					//Higer numbers means random movement is less uniform
		float noise_resolution;				//Higer numbers means random movement is more uniform
		tfxParticleFlags flags;				//flags for different states
		//Updated everyframe
		float age;							//The age of the particle, used by the controller to look up the current state on the graphs
		float max_age;						//max age before the particle expires
		tfxU32 single_loop_count;		//The number of times a single particle has looped over
		float image_frame;					//Current frame of the image if it's an animation
		float weight_acceleration;			//The current amount of gravity applied to the y axis of the particle each frame
		float intensity;					//Color is multiplied by this value in the shader to increase the brightness of the particles
		tfxRGBA8 color;						//Colour of the particle
	};

	//At the moment this struct is used in the editor only, looking to unify at some point
	struct tfxParticle {
		tfxParticleData data;
		tfxEffectEmitter *parent;				//pointer to the emitter that emitted the particle.
		//Internal use variables
		tfxParticle *next_ptr;
		tfxU32 prev_index;
	};

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
		tfxVec3 alignment;
		float stretch;
		float image_frame;
		float start_frame;
		float depth;
		void *image_ptr;
		tfxRGBA8 color;
		float intensity;
		tfxU32 alignment_type;
		bool has_frames;
	};

	static inline tfxParticleFrame ConvertToParticleFrame(const tfxParticle &p) {
		tfxParticleFrame pf;
		pf.position = p.data.world_position;
		pf.scale = p.data.scale;
		pf.alignment = p.data.alignment_vector;
		pf.stretch = p.data.velocity_normal.w;
		pf.rotations = p.data.world_rotations;
		pf.alignment_type = p.parent->GetProperties().billboard_option;
		pf.handle = p.parent->current.image_handle;
		pf.color = p.data.color;
		pf.intensity = p.data.intensity;
		pf.start_frame = p.parent->GetProperties().start_frame;
		pf.image_ptr = p.parent->GetProperties().image->ptr;
		pf.image_frame = p.data.image_frame;
		pf.has_frames = p.parent->GetProperties().image->animation_frames > 1;
		return pf;
	}

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

	//Use the particle manager to add compute effects to your scene 
	struct tfxParticleManager {
		//Particles that we can't send to the compute shader (because they have sub effects attached to them) are stored and processed here
		tfxvec<tfxParticle> particles[tfxLAYERS][2];
		tfxRangeAllocator<tfxParticle> particle_allocators[tfxLAYERS][2];
		//Only used when using distance from camera ordering. New particles are put in this list and then merge sorted into the particles buffer
		tfxvec<tfxSpawnPosition> new_positions;
		//Effects are also stored using double buffering. Effects stored here are "fire and forget", so you won't be able to apply changes to the effect in realtime. If you want to do that then 
		//you can use an tfxEffectTemplate and use callback funcitons. 
		tfxvec<tfxEffectEmitter> effects[2];
		//Set when an effect is updated and used to pass on global attributes to child emitters
		tfxParentSpawnControls parent_spawn_controls;

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
		unsigned int max_cpu_particles_per_layer;
		//The maximum number of particles that can be updated per frame per layer in the compute shader. #define tfxLAYERS to set the number of allowed layers. This is currently 4 by default
		unsigned int max_new_compute_particles;
		//The current particle buffer in use, can be either 0 or 1
		unsigned int current_pbuff;
		//The current effect buffer in use, can be either 0 or 1
		unsigned int current_ebuff;
		//The current number of particles in each buffer
		unsigned int particle_count;

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
			max_cpu_particles_per_layer(50000),
			current_ebuff(0),
			current_pbuff(0),
			highest_compute_controller_index(0),
			new_compute_particle_ptr(nullptr),
			compute_controller_ptr(nullptr),
			max_compute_controllers(10000),
			max_new_compute_particles(10000),
			new_compute_particle_index(0),
			new_particles_count(0)
		{ }
		~tfxParticleManager();
		tfxEffectEmitter &operator[] (unsigned int index);

		//Initialise the particle manager with the maximum number of particles and effects that you want the manager to update per frame
		void Init(unsigned int effects_limit = 10000, unsigned int particle_limit_per_layer = 50000);
		//Update the particle manager. Call this once per frame in your logic udpate.
		void Update();
		//When paused you still might want to keep the particles in order:
		void UpdateParticleOrderOnly();
		//Get the current particle buffer that contains all particles currently active. The particle manager can have layers in order to control draw order of particles.
		//Pass the layer that you want to get.
		tfxvec<tfxParticle> *GetParticleBuffer(unsigned int layer);
		//Add an effect to the particle manager. Pass a tfxEffectEmitter pointer if you want to change the effect on the fly. Once you add the effect to the particle manager
		//then it's location in the buffer will keep changing as effects are updated and added and removed. The tracker will be updated accordingly each frame so you will always
		//have access to the effect if you need it.
		void AddEffect(tfxEffectEmitter &effect, unsigned int buffer, bool is_sub_effect = false);
		void AddEffect(tfxEffectTemplate &effect, unsigned int buffer);
		//Clear all effects and particles in the particle manager
		void ClearAll();
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
		tfxParticle &GrabCPUParticle(unsigned int layer);
		tfxComputeParticle &GrabComputeParticle(unsigned int layer);
		void ResetParticlePtr(void *ptr);
		void ResetControllerPtr(void *ptr);
		inline unsigned int GetControllerMemoryUsage() { return highest_compute_controller_index * sizeof(tfxComputeController); }
		inline unsigned int GetParticleMemoryUsage() { return new_compute_particle_index * sizeof(tfxComputeParticle); }
		void UpdateCompute(void *sampled_particles, unsigned int sample_size = 100);
		unsigned int AddParticle(unsigned int layer, tfxParticle &p);
		//float Record(unsigned int frames, unsigned int start_frame, std::array<tfxvec<ParticleFrame>, 1000> &particle_frames);
		unsigned int ParticleCount();
		inline tfxParticle* SetNextParticle(unsigned int layer, tfxParticle &p, unsigned int buffer);
		inline tfxEffectEmitter* SetNextEffect(tfxEffectEmitter &e, unsigned int buffer);
		void UpdateBaseValues();
		tfxvec<tfxEffectEmitter> *GetEffectBuffer();
		void SetLookUpMode(tfxLookupMode mode);

		inline bool FreeCapacity(unsigned int layer, bool compute) {
			if (!compute)
				return particles[layer][current_pbuff].current_size < max_cpu_particles_per_layer;
			else
				return new_compute_particle_index < max_new_compute_particles && new_compute_particle_index < compute_global_state.end_index - compute_global_state.current_length;
		}
		inline bool FreeEffectCapacity() {
			return effects[0].current_size + effects[1].current_size < max_effects;
		}
	};

	void StopSpawning(tfxParticleManager &pm);
	void RemoveAllEffects(tfxParticleManager &pm);
	void InitParticleManager(tfxParticleManager &pm, unsigned int effects_limit, unsigned int particle_limit_per_layer);
	void AddEffect(tfxParticleManager &pm, tfxEffectEmitter &effect, float x = 0.f, float y = 0.f);
	void AddEffect(tfxParticleManager &pm, tfxEffectTemplate &effect, float x = 0.f, float y = 0.f);

	void Rotate(tfxEffectEmitter &e, float r);
	void SetAngle(tfxEffectEmitter &e, float a);
	void Scale(tfxEffectEmitter &e, const tfxVec3& s);
	void Scale(tfxEffectEmitter &e, float x, float y, float z = 1.f);
	void Position(tfxEffectEmitter &e, const tfxVec2& p);
	void Position(tfxEffectEmitter &e, const tfxVec3& p);
	void TransformEffector(tfxEffectEmitter &e, tfxParticle &parent, bool relative_position = true, bool relative_angle = false);
	void TransformEffector3d(tfxEffectEmitter &e, tfxParticle &parent, bool relative_position = true, bool relative_angle = false);
	void UpdatePMEmitter(tfxParticleManager &pm, tfxEffectEmitter &e);
	void SpawnParticles(tfxParticleManager &pm, tfxEffectEmitter &e, tfxEmitterSpawnControls &spawn_controls);
	void SpawnParticles3d(tfxParticleManager &pm, tfxEffectEmitter &e, tfxEmitterSpawnControls &spawn_controls);
	void InitCPUParticle(tfxParticleManager &pm, tfxEffectEmitter &e, tfxParticle &p, tfxEmitterSpawnControls &spawn_controls, float tween);
	tfxEmitterSpawnControls UpdateEmitterState(tfxEffectEmitter &e, tfxParentSpawnControls &parent_spawn_controls);
	tfxParentSpawnControls UpdateEffectState(tfxEffectEmitter &e);
	bool ControlParticle(tfxParticleManager &pm, tfxParticle &p, tfxEffectEmitter &e);
	static inline void SetParticleAlignment(tfxParticle &p, tfxEmitterProperties &properties) {
		bool line = p.parent->common.property_flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type == tfxEmissionType::tfxLine;
		if (p.data.flags & tfxParticleFlags_fresh && properties.vector_align_type == tfxVectorAlignType_motion) {
			p.data.alignment_vector = p.data.velocity_normal.xyz();
			p.data.alignment_vector = FastNormalizeVec(p.data.alignment_vector);
		}
		else if (properties.vector_align_type == tfxVectorAlignType_motion) {
			p.data.alignment_vector = p.data.world_position - p.data.captured_position;
			float l = FastLength(p.data.alignment_vector);
			p.data.velocity_normal.w *= l * 10.f;
			p.data.alignment_vector = FastNormalizeVec(p.data.alignment_vector);
		}
		else if (properties.vector_align_type == tfxVectorAlignType_emission) {
			p.data.alignment_vector = p.data.velocity_normal.xyz();
		}
		else if (properties.vector_align_type == tfxVectorAlignType_emitter) {
			p.data.alignment_vector = mmTransformVector(p.parent->common.transform.matrix, tfxVec4(0.f, 1.f, 0.f, 0.f)).xyz();
		}
		else {
			//Set at spawn time
		}
	}

	struct tfxEffectLibrary {
		tfxStorageMap<tfxEffectEmitter*> effect_paths;
		tfxvec<tfxEffectEmitter> effects;
		tfxStorageMap<tfxImageData> particle_shapes;
		tfxvec<tfxEffectEmitterInfo> effect_infos;
		tfxvec<tfxEmitterProperties> emitter_properties;

		tfxvec<tfxGlobalAttributes> global_graphs;
		tfxvec<tfxPropertyAttributes> property_graphs;
		tfxvec<tfxBaseAttributes> base_graphs;
		tfxvec<tfxVariationAttributes> variation_graphs;
		tfxvec<tfxOvertimeAttributes> overtime_graphs;
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
		tfxvec<tfxU32> free_property_graphs;
		tfxvec<tfxU32> free_base_graphs;
		tfxvec<tfxU32> free_variation_graphs;
		tfxvec<tfxU32> free_overtime_graphs;
		tfxvec<tfxU32> free_animation_settings;
		tfxvec<tfxU32> free_preview_camera_settings;
		tfxvec<tfxU32> free_properties;
		tfxvec<tfxU32> free_infos;

		//Get an effect from the library by index
		tfxEffectEmitter& operator[] (uint32_t index);
		tfxText name;
		bool open_library = false;
		bool dirty = false;
		tfxText library_file_path;
		tfxU32 uid = 0;

		//Todo: Inline a lot of these
		//Free everything in the library
		void Clear();
		//Get an effect in the library by it's path. So for example, if you want to get a pointer to the emitter "spark" in effect "explosion" then you could do GetEffect("explosion/spark")
		//You will need this function to apply user data and update callbacks to effects and emitters before adding the effect to the particle manager
		tfxEffectEmitter *GetEffect(tfxText &path);
		tfxEffectEmitter *GetEffect(const char *path);
		//Get an effect by it's path hash key
		tfxEffectEmitter *GetEffect(tfxKey key);
		//Get and effect by it's index
		void PrepareEffectTemplate(tfxText path, tfxEffectTemplate &effect);
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
		tfxEffectEmitter &AddFolder(tfxText name);
		tfxEffectEmitter &AddFolder(tfxEffectEmitter &effect);
		void UpdateEffectPaths();
		void AddPath(tfxEffectEmitter &effectemitter, tfxText path);
		void DeleteEffect(tfxEffectEmitter *effect);
		bool RenameEffect(tfxEffectEmitter &effect, const char *new_name);
		bool NameExists(tfxEffectEmitter &effect, const char *name);
		bool NameExists2(tfxEffectEmitter &effect, const char *name);
		void ReIndex();
		void UpdateParticleShapeReferences(tfxvec<tfxEffectEmitter> &effects, tfxU32 default_index);
		tfxEffectEmitter* MoveUp(tfxEffectEmitter &effect);
		tfxEffectEmitter* MoveDown(tfxEffectEmitter &effect);
		tfxU32 AddGlobal();
		tfxU32 AddProperty();
		tfxU32 AddBase();
		tfxU32 AddVariation();
		tfxU32 AddOvertime();
		void FreeGlobal(tfxU32 index);
		void FreeProperty(tfxU32 index);
		void FreeBase(tfxU32 index);
		void FreeVariation(tfxU32 index);
		void FreeOvertime(tfxU32 index);
		void FreeProperties(tfxU32 index);
		void FreeInfos(tfxEffectEmitter &e);
		tfxU32 CloneGlobal(tfxU32 source_index, tfxEffectLibrary *destination_library);
		tfxU32 CloneProperty(tfxU32 source_index, tfxEffectLibrary *destination_library);
		tfxU32 CloneBase(tfxU32 source_index, tfxEffectLibrary *destination_library);
		tfxU32 CloneVariation(tfxU32 source_index, tfxEffectLibrary *destination_library);
		tfxU32 CloneOvertime(tfxU32 source_index, tfxEffectLibrary *destination_library);
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
		void CompilePropertyGraph(tfxU32 index);
		void CompileBaseGraph(tfxU32 index);
		void CompileVariationGraph(tfxU32 index);
		void CompileOvertimeGraph(tfxU32 index);
		void CompileColorGraphs(tfxU32 index);
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
		tfxEffectEmitter effect_template;

		void AddPath(tfxEffectEmitter &effectemitter, tfxText path) {
			paths.Insert(path, &effectemitter);
			for (auto &sub : effectemitter.common.library->GetInfo(effectemitter).sub_effectors) {
				tfxText sub_path = path;
				sub_path.Appendf("/%s", sub.common.library->GetInfo(sub).name.c_str());
				AddPath(sub, sub_path);
			}
		}

		inline tfxEffectEmitter &Effect() { return effect_template; }
		inline tfxEffectEmitter *Get(tfxText path) { if (paths.ValidName(path)) return paths.At(path); return nullptr; }
		inline void SetUserData(tfxText path, void *data) { if (paths.ValidName(path)) paths.At(path)->user_data = data; }
		inline void SetUserData(void *data) { effect_template.user_data = data; }
		void SetUserDataAll(void *data);
	};

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
		tfxText key;
		tfxText str_value;
		int int_value;
		bool bool_value;
		float float_value;
		double double_value;
	};

	struct tfxDataTypesDictionary {
		tfxStorageMap<tfxDataType> names_and_types;

		tfxDataTypesDictionary();
	};

	static tfxDataTypesDictionary data_types;

	//Internal functions
	//Some file IO functions
	bool HasDataValue(tfxStorageMap<tfxDataEntry> &config, tfxText key);
	void AddDataValue(tfxStorageMap<tfxDataEntry> &config, tfxText key, const char *value);
	void AddDataValue(tfxStorageMap<tfxDataEntry> &config, tfxText key, int value);
	void AddDataValue(tfxStorageMap<tfxDataEntry> &config, tfxText key, bool value);
	void AddDataValue(tfxStorageMap<tfxDataEntry> &config, tfxText key, double value);
	void AddDataValue(tfxStorageMap<tfxDataEntry> &config, tfxText key, float value);
	tfxText &GetDataStrValue(tfxStorageMap<tfxDataEntry> &config, const char* key);
	int& GetDataIntValue(tfxStorageMap<tfxDataEntry> &config, const char* key);
	float& GetDataFloatValue(tfxStorageMap<tfxDataEntry> &config, const char* key);
	bool SaveDataFile(tfxStorageMap<tfxDataEntry> &config, const char* path = "");
	bool LoadDataFile(tfxStorageMap<tfxDataEntry> &config, const char* path);
	void StreamProperties(tfxEmitterProperties &property, tfxEmitterPropertyFlags &flags, tfxText &file);
	void StreamProperties(tfxEffectEmitter &effect, tfxText &file);
	void StreamGraph(const char * name, tfxGraph &graph, tfxText &file);
	tfxvec<tfxText> SplitString(const tfxText &s, char delim = 61);
	bool StringIsUInt(const tfxText &s);
	int GetDataType(const tfxText &s);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxText &field, uint32_t value);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxText &field, float value);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxText &field, bool value);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxText &field, int value);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxText &field, tfxText &value);
	void AssignGraphData(tfxEffectEmitter &effect, tfxvec<tfxText> &values);
	void AssignNodeData(tfxAttributeNode &node, tfxvec<tfxText> &values);
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
	static inline void TransformParticle(tfxParticleData &data, const tfxCommon &common, const tfxVec3 &from_position) {
		data.world_position = data.local_position;
		data.world_rotations.roll = data.local_rotations.roll;
	}
	static inline void TransformParticleAngle(tfxParticleData &data, const tfxCommon &common, const tfxVec3 &from_position) {
		data.world_position = data.local_position;
		data.world_rotations.roll = common.transform.world_rotations.roll + data.local_rotations.roll;
	}
	static inline void TransformParticleRelative(tfxParticleData &data, const tfxCommon &common, const tfxVec3 &from_position) {
		data.world_rotations.roll = data.local_rotations.roll;
		float s = sin(data.local_rotations.roll);
		float c = cos(data.local_rotations.roll);
		tfxMatrix2 pmat;
		pmat.Set(c, s, -s, c);
		pmat = pmat.Transform(common.transform.matrix);
		tfxVec2 rotatevec = mmTransformVector(common.transform.matrix, tfxVec2(data.local_position.x, data.local_position.y) + common.handle.xy());
		data.world_position = from_position.xy() + rotatevec * common.transform.scale.xy();
	}
	static inline void TransformParticleRelativeLine(tfxParticleData &data, const tfxCommon &common, const tfxVec3 &from_position) {
		data.world_rotations.roll = common.transform.world_rotations.roll + data.local_rotations.roll;
		float s = sin(data.local_rotations.roll);
		float c = cos(data.local_rotations.roll);
		tfxMatrix2 pmat;
		pmat.Set(c, s, -s, c);
		pmat = pmat.Transform(common.transform.matrix);
		tfxVec2 rotatevec = mmTransformVector(common.transform.matrix, tfxVec2(data.local_position.x, data.local_position.y) + common.handle.xy());
		data.world_position = from_position.xy() + rotatevec * common.transform.scale.xy();
	}
	static inline void TransformParticle3dPositions(tfxParticleData &data, const tfxCommon &common, const tfxVec3 &from_position) {
		data.world_position = data.local_position;
	}
	static inline void TransformParticle3dPositionsRelative(tfxParticleData &data, const tfxCommon &common, const tfxVec3 &from_position) {
		tfxVec4 rotatevec = mmTransformVector(common.transform.matrix, data.local_position + common.handle);
		data.world_position = common.transform.world_position + rotatevec.xyz();
	}
	static inline void TransformParticle3d(tfxParticleData &data, const tfxCommon &common, const tfxVec3 &from_position) {
		data.world_position = data.local_position;
		data.world_rotations = data.local_rotations;
	}
	static inline void TransformParticle3dAngle(tfxParticleData &data, const tfxCommon &common, const tfxVec3 &from_position) {
		data.world_position = data.local_position;
		data.world_rotations = common.transform.world_rotations + data.local_rotations;
	}
	static inline void TransformParticle3dRelative(tfxParticleData &data, const tfxCommon &common, const tfxVec3 &from_position) {
		data.world_rotations = data.local_rotations;
		float s = sin(data.local_rotations.roll);
		float c = cos(data.local_rotations.roll);
		tfxMatrix2 pmat;
		pmat.Set(c, s, -s, c);
		pmat = pmat.Transform(common.transform.matrix);
		tfxVec4 rotatevec = mmTransformVector(common.transform.matrix, data.local_position + common.handle);
		data.world_position = from_position + rotatevec.xyz() * common.transform.scale;
	}
	static inline void TransformParticle3dRelativeLine(tfxParticleData &data, const tfxCommon &common, const tfxVec3 &from_position) {
		data.world_rotations = data.local_rotations;
		float s = sin(data.local_rotations.roll);
		float c = cos(data.local_rotations.roll);
		tfxMatrix2 pmat;
		pmat.Set(c, s, -s, c);
		pmat = pmat.Transform(common.transform.matrix);
		tfxVec4 rotatevec = mmTransformVector(common.transform.matrix, data.local_position + common.handle);
		data.world_position = from_position + rotatevec.xyz() * common.transform.scale;
	}
	static inline int SortParticles(void const *left, void const *right) {
		float d1 = static_cast<const tfxSpawnPosition*>(left)->distance_to_camera;
		float d2 = static_cast<const tfxSpawnPosition*>(right)->distance_to_camera;
		return (d2 > d1) - (d2 < d1);
	}
	static inline void InsertionSortParticles(tfxvec<tfxParticle> &particles, tfxvec<tfxParticle> &current_buffer) {
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
	tfxVec3 Tween(float tween, tfxVec3 &world, tfxVec3 &captured);
	float Interpolatef(float tween, float, float);
	int ValidateEffectPackage(const char *filename);
	void ReloadBaseValues(tfxParticle &p, tfxEffectEmitter &e);

	//Particle initialisation functions, one for 2d one for 3d effects
	void InitialiseParticle2d(tfxParticleData &data, tfxEmitterState &emitter, tfxCommon &common, tfxEmitterSpawnControls &spawn_values, tfxEffectEmitter *library_link, float tween);
	tfxSpawnPosition InitialisePosition3d(tfxEmitterState &current, tfxCommon &common, tfxEmitterSpawnControls &spawn_values, tfxEffectEmitter *library_link, float tween);
	void InitialiseParticle3d(tfxParticleData &data, tfxEmitterState &current, tfxCommon &common, tfxEmitterSpawnControls &spawn_values, tfxEffectEmitter *library_link, float tween);
	//void InitialisePostion3d(tfxParticle &p, tfxEmitter &emitter, tfxEmitterSpawnControls &spawn_values);
	void UpdateParticle2d(tfxParticleData &data, tfxControlData &c, tfxEffectEmitter *library_link);
	void UpdateParticle3d(tfxParticleData &data, tfxControlData &c, tfxEffectEmitter *library_link);

	//Helper functions

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
	int LoadEffectLibraryPackage(const char *filename, tfxEffectLibrary &lib, void(*shape_loader)(const char *filename, tfxImageData &image_data, void *raw_image_data, int image_size, void *user_data) = nullptr, void *user_data = nullptr);

	//---
	//Prepare an effect template for setting up function call backs to customise the behaviour of the effect in realtime
	//Returns true on success.
	bool PrepareEffectTemplate(tfxEffectLibrary &library, const char *name, tfxEffectTemplate &effect_template);
}

