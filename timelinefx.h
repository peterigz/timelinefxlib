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
			unsigned int anim_index = example->particle_textures->AddAnimation(image, (unsigned int)image_data.image_size.x, (unsigned int)image_data.image_size.y, (unsigned int)image_data.animation_frames);
			//Important step: you need to point the ImageData.ptr to the appropriate handle in the renderer to point to the texture of the particle shape
			//You'll need to use this in your render function to tell your renderer which texture to use to draw the particle
			image_data.ptr = &example->particle_textures->GetAnimation(anim_index);
		}
		else {
			//Add the image to the texture in our renderer
			unsigned int image_index = example->particle_textures->AddImage(image);
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
		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			//Use GetParticleBuffer(layer) to get all of the particles in the current layer
			for (auto p : *pm.GetParticleBuffer(layer)) {
				//In order to set the correct blendmode we need to get the property from the parent emitter that emitted the particle
				//A pointer to the parent emitter is stored in the parent member
				tfx::EffectEmitter &e = *p.parent;

				//Set the correct blendmode, see timelinefx::BlendMode. You may have to map the blendmodes depending on the renderer you use
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
class SimplexNoise {
public:
	// 1D Perlin simplex noise
	static float noise(float x);
	// 2D Perlin simplex noise
	static float noise(float x, float y);
	// 3D Perlin simplex noise
	static float noise(float x, float y, float z);

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
	explicit SimplexNoise(float frequency = 1.0f,
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



namespace tfx {

#define TWO63 0x8000000000000000u 
#define TWO64f (TWO63*2.0)

	//----------------------------------------------------------
	//Forward declarations

	struct EffectEmitter;
	struct EffectorStore;
	struct Particle;
	struct ComputeSprite;
	struct ParticleSprite;
	struct ComputeParticle;
	struct AnimationSettings;
	struct EffectLibrary;
	struct tfxText;
	struct tfxEffect;
	struct tfxEmitter;
	struct tfxEffectPool;

	//--------------------------------------------------------------
	//macros
#define TFX_VERSION "Alpha"
#define TFX_VERSION_NUMBER 3.29.2022

#define tfxMAX_FRAME 20000.f
#define tfxNullParent 0xFFFFFFFF
#define EmitterPropertiesCount 26

#define Del << "=" <<
#define Com << "," <<
#define EndLine << std::endl

#define Delt "=" 
#define Comt ","
#define EndLinet "\n"

typedef std::chrono::high_resolution_clock Clock;

//Override this for more layers, although currently the editor is fixed at 4
#ifndef tfxLAYERS
#define tfxLAYERS 4
#define EachLayer int layer = 0; layer !=tfxLAYERS; ++layer
#endif 

//type defs
typedef unsigned int u32;
typedef int s32;
typedef unsigned long long u64;
typedef long long s64;
typedef unsigned int tfxEffectID;

	//----------------------------------------------------------
	//enums/flags

	//Particle property that defines how a particle will rotate
	enum AngleSetting : unsigned char {
		tfxAlign,												//Align the particle with it's direction of travel
		tfxRandom,												//Chose a random angle at spawn time/flags
		tfxSpecify,												//Specify the angle at spawn time
		tfxAlignWithEmission,									//Align the particle with the emission direction only
		tfxGraph												//Unused, but the idea was to allow the angle to be changed overtime using a graph
	};

	//Blend mode property of the emitter
	//It's up to whoever is implementing this library to provide a render function for the particles and make use of these blend modes
	enum BlendMode : unsigned char {
		tfxNone = 0,												//Basically not used, only alpha and additive are used
		tfxAlpha = 1,												//Alpha blend the particle with what it's being drawn on 
		tfxAdditive = 2												//Add the color of the particle with what it's being drawn on
	};

	//Graph presets to determine limits and scales of different graphs, mainly used for the editor
	enum GraphPreset {
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
		tfxOpacityOvertimePreset,
		tfxColorPreset,
		tfxPercentOvertime,
		tfxIntensityOvertimePreset
	};

	enum GraphCategory : unsigned int {
		tfxGraphCategory_global,
		tfxGraphCategory_property,
		tfxGraphCategory_base,
		tfxGraphCategory_variation,
		tfxGraphCategory_overtime
	};


#define tfxGlobalCount  13
#define	tfxPropertyCount  8
#define	tfxBaseCount  8
#define	tfxVariationCount  9
#define	tfxOvertimeCount  16

#define tfxGlobalStart 0
#define	tfxPropertyStart tfxGlobalCount
#define	tfxBaseStart (tfxPropertyStart + tfxPropertyCount)
#define	tfxVariationStart (tfxBaseStart + tfxBaseCount)
#define	tfxOvertimeStart (tfxVariationStart + tfxVariationCount)

	//All the different types of graphs, split into main type: global, property, base, variation and overtime
	enum GraphType : unsigned char {
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
		tfxGlobal_effect_angle,

		tfxProperty_emission_angle,
		tfxProperty_emission_range,
		tfxProperty_emitter_angle,
		tfxProperty_splatter,
		tfxProperty_emitter_width,
		tfxProperty_emitter_height,
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

	//EffectEmitter type - effect contains emitters, and emitters spawn particles, but they both share the same struct for simplicity
	enum EffectEmitterType : unsigned char {
		tfxEffectType,
		tfxEmitterType,
		tfxStage,
		tfxFolder
	};

	//Different ways that particles can be emitted
	enum EmissionType : unsigned char {
		tfxPoint,
		tfxArea,
		tfxLine,
		tfxEllipse
	};

	//Determines how for area, line and ellipse emitters the direction that particles should travel
	enum EmissionDirection : unsigned char {
		tfxInwards,
		tfxOutwards,
		tfxBothways,
		tfxSpecified
	};

	//For line effects where traverse line is switched on
	enum LineTraversalEndBehaviour : unsigned char {
		tfxLoop,
		tfxKill,
		tfxLetFree
	};

	//Mainly for the editor, maybe this can just be moved there instead?
	enum ExportColorOptions {
		tfxFullColor,
		tfxOneColor,
		tfxGreyScale,
		tfxOneColorAlpha,
		tfxGreyScaleAlhpa
	};

	//Mainly for the editor, maybe this can just be moved there instead?
	enum ExportOptions {
		tfxSpriteSheet,
		tfxStrip,
		tfxSeparateImages
	};

	//Graph data can be looked up in one of 2 ways, either by just using linear/bezier interpolation (slower), or look up the value in a pre-compiled look up table.
	enum LookupMode {
		tfxPrecise,
		tfxFast
	};

	//Used in file loading - for loading effects library
	enum DataType {
		tfxString,
		tfxSInt,
		tfxUint,
		tfxFloat,
		tfxDouble,
		tfxBool
	};
	
	//Block designators for loading effects library
	enum EffectLibraryStream : uint32_t {
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
		tfxEndFolder
	};

	typedef unsigned int tfxEmitterPropertyFlags;
	typedef unsigned int tfxVectorFieldFlags;
	typedef unsigned char tfxParticleFlags;
	typedef unsigned int tfxEmitterStateFlags;
	typedef unsigned int tfxParticleControlFlags;
	typedef unsigned int tfxAttributeNodeFlags;

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
		tfxParticleControlFlags_random = 1 << 19,
		tfxParticleControlFlags_specify = 1 << 20,
		tfxParticleControlFlags_alpha = 1 << 21,
		tfxParticleControlFlags_additive = 1 << 22,
	};

	enum tfxEmitterPropertyFlags_ {
		tfxEmitterPropertyFlags_none = 0,
		tfxEmitterPropertyFlags_random_color = 1 << 0,						//Pick a random color from the color overtime gradient rather then change the color over the lifetime of the particle
		tfxEmitterPropertyFlags_relative_position = 1 << 1,					//Keep the particles position relative to the current position of the emitter
		tfxEmitterPropertyFlags_relative_angle = 1 << 2,					//Keep the angle of the particles relative to the current angle of the emitter
		tfxEmitterPropertyFlags_image_handle_auto_center = 1 << 3,			//Set the offset of the particle to the center of the image
		tfxEmitterPropertyFlags_single = 1 << 4,							//Only spawn a single particle (or number of particles specified by spawn_amount) that does not expire
		tfxEmitterPropertyFlags_one_shot = 1 << 5,							//Only spawn a single particle (or number of particles specified by spawn_amount) in one go
		tfxEmitterPropertyFlags_spawn_on_grid = 1 << 6,						//When using an area, line or ellipse emitter, spawn along a grid
		tfxEmitterPropertyFlags_grid_spawn_clockwise = 1 << 7,				//Spawn clockwise/left to right around the area
		tfxEmitterPropertyFlags_fill_area = 1 << 8,							//Fill the area *not implemented yet*
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
		tfxEmitterPropertyFlags_use_vector_field = 1 << 19,					//Enable the use of a vector field to apply forces to the particles
		tfxEmitterPropertyFlags_is_in_folder = 1 << 20,						//This effect is located inside a folder
		tfxEmitterPropertyFlags_is_bottom_emitter = 1 << 21,				//This emitter has no child effects, so can spawn particles that could be used in a compute shader if it's enabled
		tfxEmitterPropertyFlags_use_spawn_ratio = 1 << 22,					//Option for area emitters to multiply the amount spawned by a ration of particles per pixels squared
		tfxEmitterPropertyFlags_can_grow_particle_memory = 1 << 23,			//Allows for expanding the memory used for particle emitters if the amount spawned is changed dynamically
	};

	enum tfxParticleFlags_ : unsigned char {
		tfxParticleFlags_none = 0,
		tfxParticleFlags_fresh = 1 << 0,									//Particle has just spawned this frame	
		tfxParticleFlags_capture_after_transform = 1 << 3,					//Particle will be captured after a transfrom, used for traversing lines and looping back to the beginning to avoid lerping imbetween
		tfxParticleFlags_remove = 1 << 4,									//Particle will be removed this or next frame
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
		tfxEmitterStateFlags_align_with_velocity = 1 << 12,
		tfxEmitterStateFlags_lifetime_uniform_size = 1 << 13,			//Keep the size over lifetime of the particle uniform
		tfxEmitterStateFlags_loop = 1 << 14,
		tfxEmitterStateFlags_kill = 1 << 15,
		tfxEmitterStateFlags_play_once = 1 << 16,						//Play the animation once only
		tfxEmitterStateFlags_single_shot_done = 1 << 17,
		tfxEmitterStateFlags_is_line = 1 << 18,
		tfxEmitterStateFlags_is_area = 1 << 19
	};

	enum tfxVectorFieldFlags_: unsigned char {
		tfxVectorFieldFlags_none = 0,
		tfxVectorFieldFlags_repeat_horizontal = 1 << 0,							//Field will repeat horizontally
		tfxVectorFieldFlags_repeat_vertical = 1 << 1								//Field will repeat vertically
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

	const float kLIFE_MIN = 0.f;
	const float kLIFE_MAX = 100000.f;
	const float kLIFE_STEPS = 200.f;

	const float kAMOUNT_MIN = 0.f;
	const float kAMOUNT_MAX = 5000.f;
	const float kAMOUNT_STEPS = 100.f;

	const float kGLOBAL_PERCENT_MIN = 0.f;
	const float kGLOBAL_PERCENT_MAX = 20.f;
	const float kGLOBAL_PERCENT_STEPS = 100.f;

	const float kGLOBAL_PERCENT_V_MIN = 0.f;
	const float kGLOBAL_PERCENT_V_MAX = 10.f;
	const float kGLOBAL_PERCENT_V_STEPS = 200.f;

	const float kINTENSITY_MIN = 0.f;
	const float kINTENSITY_MAX = 5.f;
	const float kINTENSITY_STEPS = 100.f;

	const float kANGLE_MIN = 0.f;
	const float kANGLE_MAX = 1080.f;
	const float kANGLE_STEPS = 54.f;

	const float kARC_MIN = 0.f;
	const float kARC_MAX = 6.28319f;
	const float kARC_STEPS = .3141595f;

	const float kEMISSION_RANGE_MIN = 0.f;
	const float kEMISSION_RANGE_MAX = 180.f;
	const float kEMISSION_RANGE_STEPS = 30.f;

	const float kDIMENSIONS_MIN = 0.f;
	const float kDIMENSIONS_MAX = 4000.f;
	const float kDIMENSIONS_STEPS = 40.f;

	const float kVELOCITY_MIN = 0.f;
	const float kVELOCITY_MAX = 10000.f;
	const float kVELOCITY_STEPS = 100.f;

	const float kVELOCITY_OVERTIME_MIN = -20.f;
	const float kVELOCITY_OVERTIME_MAX = 20.f;
	const float kVELOCITY_OVERTIME_STEPS = 200.f;

	const float kWEIGHT_MIN = -2500.f;
	const float kWEIGHT_MAX = 2500.f;
	const float kWEIGHT_STEPS = 200.f;

	const float kWEIGHT_VARIATION_MIN = 0.f;
	const float kWEIGHT_VARIATION_MAX = 2500.f;
	const float kWEIGHT_VARIATION_STEPS = 250.f;

	const float kSPIN_MIN = -2000.f;
	const float kSPIN_MAX = 2000.f;
	const float kSPIN_STEPS = 100.f;

	const float kSPIN_VARIATION_MIN = 0.f;
	const float kSPIN_VARIATION_MAX = 2000.f;
	const float kSPIN_VARIATION_STEPS = 100.f;

	const float kSPIN_OVERTIME_MIN = -20.f;
	const float kSPIN_OVERTIME_MAX = 20.f;
	const float kSPIN_OVERTIME_STEPS = 200.f;

	const float kDIRECTION_OVERTIME_MIN = 0.f;
	const float kDIRECTION_OVERTIME_MAX = 4320.f;
	const float kDIRECTION_OVERTIME_STEPS = 216.f;

	const float kFRAMERATE_MIN = 0.f;
	const float kFRAMERATE_MAX = 200.f;
	const float kFRAMERATE_STEPS = 100.f;

	const float kMAX_DIRECTION_VARIATION = 22.5f;
	const float kMAX_VELOCITY_VARIATION = 30.f;
	const int kMOTION_VARIATION_INTERVAL = 30;
	
	//these Variables determine the timing resolution that particles are updated at. So an Update frequency of 60 would mean that the particles are updated at 60 frames per second.
	extern float UPDATE_FREQUENCY;
	extern float UPDATE_TIME;
	extern float FRAME_LENGTH;

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
		unsigned int current_size;
		unsigned int capacity;
		T* data;

		// Provide standard typedefs but we don't use them ourselves.
		typedef T                   value_type;
		typedef value_type*         iterator;
		typedef const value_type*   const_iterator;

		inline tfxvec() { current_size = capacity = 0; data = NULL; }
		inline tfxvec(unsigned int qty) { current_size = capacity = 0; data = NULL; resize(qty); }
		inline tfxvec(T* from, T* to) { current_size = capacity = 0; data = NULL; auto current = from; while (current != to + 1) { push_back(*current); ++current; } }
		inline tfxvec(std::initializer_list<T> t) { current_size = capacity = 0; data = NULL; for (T element : t) { push_back(element); } }
		inline tfxvec(const tfxvec<T> &src) { current_size = capacity = 0; data = NULL; resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); }
		inline tfxvec<T>& operator=(const tfxvec<T>& src) { clear(); resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); return *this; }
		inline ~tfxvec() { if (data) free(data); data = NULL; current_size = capacity = 0; }

		inline bool			empty() { return current_size == 0; }
		inline unsigned int			size() { return current_size; }
		inline const unsigned int	size() const { return current_size; }
		inline T&           operator[](unsigned int i) { return data[i]; }
		inline const T&     operator[](unsigned int i) const { assert(i < current_size); return data[i]; }

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

		inline unsigned int          _grow_capacity(unsigned int sz) const { unsigned int new_capacity = capacity ? (capacity + capacity / 2) : 8; return new_capacity > sz ? new_capacity : sz; }
		inline void         resize(unsigned int new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
		inline void         resize(unsigned int new_size, const T& v) { if (new_size > capacity) reserve(_grow_capacity(new_size)); if (new_size > current_size) for (unsigned int n = current_size; n < new_size; n++) memcpy(&data[n], &v, sizeof(v)); current_size = new_size; }
		inline void         shrink(unsigned int new_size) { assert(new_size <= current_size); current_size = new_size; }
		inline void         reserve(unsigned int new_capacity) { if (new_capacity <= capacity) return; T* new_data = (T*)malloc((size_t)new_capacity * sizeof(T)); if (data) { memcpy(new_data, data, (size_t)current_size * sizeof(T)); free(data); } data = new_data; capacity = new_capacity; }

		inline T&	        push_back(const T& v) { if (current_size == capacity) reserve(_grow_capacity(current_size + 1)); new((void*)(data + current_size)) T(v); current_size++; return data[current_size - 1]; }
		inline void         pop() { assert(current_size > 0); current_size--; }
		inline T&	        pop_back() { assert(current_size > 0); current_size--; return data[current_size]; }
		inline void         push_front(const T& v) { if (current_size == 0) push_back(v); else insert(data, v); }
		inline T*           erase(const T* it) { assert(it >= data && it < data + current_size); const ptrdiff_t off = it - data; memmove(data + off, data + off + 1, ((size_t)current_size - (size_t)off - 1) * sizeof(T)); current_size--; return data + off; }
		inline T*           erase(const T* it, const T* it_last) { assert(it >= data && it < data + current_size && it_last > it && it_last <= data + current_size); const ptrdiff_t count = it_last - it; const ptrdiff_t off = it - data; memmove(data + off, data + off + count, ((size_t)current_size - (size_t)off - count) * sizeof(T)); current_size -= (unsigned int)count; return data + off; }
		inline T*           erase_unsorted(const T* it) { assert(it >= data && it < data + current_size);  const ptrdiff_t off = it - data; if (it < data + current_size - 1) memcpy(data + off, data + current_size - 1, sizeof(T)); current_size--; return data + off; }
		inline T*           insert(const T* it, const T& v) { assert(it >= data && it <= data + current_size); const ptrdiff_t off = it - data; if (current_size == capacity) reserve(_grow_capacity(current_size + 1)); if (off < (ptrdiff_t)current_size) memmove(data + off + 1, data + off, ((size_t)current_size - (size_t)off) * sizeof(T)); new((void*)(data + off)) T(v); current_size++; return data + off; }
		inline bool         contains(const T& v) const { const T* _data = data;  const T* data_end = data + current_size; while (_data < data_end) if (*_data++ == v) return true; return false; }
		inline T*           find(const T& v) { T* _data = data;  const T* data_end = data + current_size; while (_data < data_end) if (*_data == v) break; else ++_data; return _data; }
		inline const T*     find(const T& v) const { const T* _data = data;  const T* data_end = data + current_size; while (_data < data_end) if (*_data == v) break; else ++_data; return _data; }
		inline bool         find_erase(const T& v) { const T* it = find(v); if (it < data + current_size) { erase(it); return true; } return false; }
		inline bool         find_erase_unsorted(const T& v) { const T* it = find(v); if (it < data + current_size) { erase_unsorted(it); return true; } return false; }
		inline unsigned int          index_from_ptr(const T* it) const { assert(it >= data && it < data + current_size); const ptrdiff_t off = it - data; return (unsigned int)off; }

		inline void			create_pool(unsigned int amount) { assert(current_size == 0); T base; reserve(amount); for (unsigned int i = 0; i != capacity; ++i) { new((void*)(data + current_size)) T(base); current_size++; } }
		inline void			create_pool_with(unsigned int amount, const T &base) { assert(current_size == 0);  reserve(amount); for (unsigned int i = 0; i != capacity; ++i) { new((void*)(data + current_size)) T(base); current_size++; } }

	};

	//A char buffer you can use to load a file into and read from
	//Has no deconstructor so make sure you call FreeAll() when done
	//This is meant for limited usage in timeline fx only and not recommended for use outside!
	struct tfxstream {
		u64 size = 0;
		u64 position = 0;
		char* data = NULL;

		inline tfxstream() { size = position = 0; data = NULL; }
		inline tfxstream(u64 qty) { size = position = 0; data = NULL; Resize(qty); }
		inline tfxstream(const tfxstream &src) { size = 0; data = NULL; Resize(src.size); memcpy(data, src.data, (u64)size * sizeof(char)); }

		inline bool Read(char* dst, u64 count) {
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
		inline bool Write(void *src, u64 count) {
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
		inline void Seek(u64 offset) {
			if (offset < size)
				position = offset;
			else
				position = size;
		}

		inline bool			Empty() { return size == 0; }
		inline u64			Size() { return size; }
		inline const u64	Size() const { return size; }

		inline void			FreeAll() { if (data) { size = size = 0; free(data); data = NULL; } }
		inline void         Clear() { if (data) { size = 0; } }

		inline void         Resize(u64 new_capacity) { if (new_capacity <= size) return; char* new_data = (char*)malloc((u64)new_capacity * sizeof(char)); if (data) { memcpy(new_data, data, (u64)size * sizeof(char)); free(data); } data = new_data; size = new_capacity; position = 0; }
		inline void			NullTerminate() { *(data + size) = NULL; }

	};

	template <typename T> int sgn(T val) {
		return (T(0) < val) - (val < T(0));
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
		inline float Squared() { return x * x + y * y; }
		inline bool IsNill() { return !x && !y; }
	};
	inline tfxVec2 operator*(float ls, tfxVec2 rs) { return tfxVec2(rs.x * ls, rs.y * ls); }

	struct tfxVec3 {
		float x, y, z;

		tfxVec3() { x = y = z = 0.f; }
		tfxVec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}

		inline tfxVec3 operator+(tfxVec3 v) { return tfxVec3(x + v.x, y + v.y, z + v.z); }
		inline tfxVec3 operator+=(tfxVec3 v) { return tfxVec3(x + v.x, y + v.y, z + v.z); }
		inline tfxVec3 operator-(tfxVec3 v) { return tfxVec3(x - v.x, y - v.y, z - v.z); }
		inline tfxVec3 operator-() { return tfxVec3(-x, -y, -z); }
		inline tfxVec3 operator-=(tfxVec3 v) { return tfxVec3(-x, -y, -z); }
		inline tfxVec3 operator*(tfxVec3 v) { return tfxVec3(x * v.x, y * v.y, z * v.z); }
		inline tfxVec3 operator*=(tfxVec3 v) { return tfxVec3(x * v.x, y * v.y, z * v.z); }
		inline tfxVec3 operator/(tfxVec3 v) { return tfxVec3(x / v.x, y / v.y, z / v.z); }
		inline tfxVec3 operator/=(tfxVec3 v) { return tfxVec3(x / v.x, y / v.y, z / v.z); }
		inline tfxVec3 operator+(float v) { return tfxVec3(x + v, y + v, z + v); }
		inline tfxVec3 operator-(float v) { return tfxVec3(x - v, y - v, z - v); }
		inline tfxVec3 operator*(float v) { return tfxVec3(x * v, y * v, z * v); }
		inline tfxVec3 operator/(float v) { return tfxVec3(x / v, y / v, z / v); }
		inline void operator*=(float v) { x *= v; y *= v; z *= v; }
		inline float Squared() { return x * x + y * y + z * z; }
		inline bool IsNill() { return !x && !y && !z; }
	};

	struct tfxVec4 {
		float x, y, z, w;

		tfxVec4() { x = y = z = w = 0.f; }
		tfxVec4(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}
		tfxVec4(tfxVec2 vec1, tfxVec2 vec2) : x(vec1.x), y(vec1.y), z(vec2.x), w(vec2.y) {}

		inline tfxVec4 operator+(tfxVec4 v) { return tfxVec4(x + v.x, y + v.y, z + v.z, w + v.w); }
		inline tfxVec4 operator+=(tfxVec4 v) { return tfxVec4(x + v.x, y + v.y, z + v.z, w + v.w); }
		inline tfxVec4 operator-(tfxVec4 v) { return tfxVec4(x - v.x, y - v.y, z - v.z, w - v.w); }
		inline tfxVec4 operator-() { return tfxVec4(-x, -y, -z, -w); }
		inline tfxVec4 operator-=(tfxVec4 v) { return tfxVec4(-x, -y, -z, -w); }
		inline tfxVec4 operator*(tfxVec4 v) { return tfxVec4(x * v.x, y * v.y, z * v.z, w * v.w); }
		inline tfxVec4 operator*=(tfxVec4 v) { return tfxVec4(x * v.x, y * v.y, z * v.z, w * v.w); }
		inline tfxVec4 operator/(tfxVec4 v) { return tfxVec4(x / v.x, y / v.y, z / v.z, w / v.w); }
		inline tfxVec4 operator/=(tfxVec4 v) { return tfxVec4(x / v.x, y / v.y, z / v.z, w / v.w); }
		inline tfxVec4 operator+(float v) { return tfxVec4(x + v, y + v, z + v, w + v); }
		inline tfxVec4 operator-(float v) { return tfxVec4(x - v, y - v, z - v, w - v); }
		inline tfxVec4 operator*(float v) { return tfxVec4(x * v, y * v, z * v, w * v); }
		inline tfxVec4 operator/(float v) { return tfxVec4(x / v, y / v, z / v, w / v); }
	};

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

	inline float tfxRadians(float degrees) { return degrees * 0.01745329251994329576923690768489f; }
	inline float tfxDegrees(float radians) { return radians * 57.295779513082320876798154814105f; }
	inline void tfxBound(tfxVec2 &s, tfxVec2 &b) { if (s.x < 0.f) s.x = 0.f; if (s.y < 0.f) s.y = 0.f; if (s.x >= b.x) s.x = b.x - 1; if (s.y >= b.y) s.y = b.y - 1; }

	//Very simple 2D Matix
	struct Matrix2 {

		float aa, ab, ba, bb;

		Matrix2() :aa(1.f), ab(0.f), ba(0.f), bb(1.f) {}

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

		Matrix2 Transform(Matrix2 m) {
			Matrix2 r;
			r.aa = aa * m.aa + ab * m.ba; r.ab = aa * m.ab + ab * m.bb;
			r.ba = ba * m.aa + bb * m.ba; r.bb = ba * m.ab + bb * m.bb;
			return r;
		}

		tfxVec2 TransformVector(tfxVec2 v) {
			tfxVec2 tv = tfxVec2(0.f, 0.f);
			tv.x = v.x * aa + v.y * ba;
			tv.y = v.x * ab + v.y * bb;
			return tv;
		}

	};

	inline float Clamp(float lower, float upper, float value) {
		if (value < lower) return lower;
		if (value > upper) return upper;
		return value;
	}

	inline float Distance(float fromx, float fromy, float tox, float toy) {

		float w = tox - fromx;
		float h = toy - fromy;

		return std::sqrt(w * w + h * h);

	}

	inline tfxVec2 InterpolateVec2(float tween, tfxVec2 &from, tfxVec2 &to) {
		return from * tween + to * (1.f - tween);
	}

	inline float Interpolatef(float tween, float from, float to) {
		return from * tween + to * (1.f - tween);
	}

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
	class XXHash64
	{
	public:
		/// create new XXHash (64 bit)
		/** @param seed your seed value, even zero is a valid seed **/
		explicit XXHash64(uint64_t seed)
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
			XXHash64 hasher(seed);
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
		tfxText(const char *text) { size_t length = strnlen_s(text, 512); if (!length) { Clear(); return; }; if (string.capacity < length) string.reserve((unsigned int)length); memcpy(string.data, text, length); string.current_size = (unsigned int)length; NullTerminate(); }
		inline char operator[](unsigned int i) const { assert(i < string.current_size); return string[i]; }
		inline char operator[](unsigned int i) { assert(i < string.current_size); return string[i]; }
		inline void operator=(const char *text) { size_t length = strnlen_s(text, 512); if (!length) { Clear(); return; }; if (string.capacity < length) string.reserve((unsigned int)length); memcpy(string.data, text, length); string.current_size = (unsigned int)length; NullTerminate(); }
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
		inline unsigned int Length() const { return string.current_size ? string.current_size - 1 : 0; }
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
			unsigned int index;
			pair(tfxKey k, unsigned int i) : key(k), index(i) {}
		};

		tfxvec<pair> map;
		tfxvec<T> data;
		void(*remove_callback)(T &item) = nullptr;

		tfxStorageMap() {}

		//Insert a new T value into the storage
		inline void Insert(const char *name, const T &value) {
			tfxKey key = XXHash64::hash(name, strlen(name), 0);
			SetIndex(key, value);
		}

		//Insert a new T value into the storage
		inline void Insert(tfxText name, const T &value) {
			tfxKey key = XXHash64::hash(name.c_str(), name.Length(), 0);
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

		inline unsigned int Size() {
			return data.current_size;
		}

		inline unsigned int LastIndex() {
			return data.current_size - 1;
		}

		inline bool ValidIndex(unsigned int index) {
			return index < data.current_size;
		}

		inline bool ValidName(const char *name) {
			return GetIndex(name) > -1;
		}

		inline bool ValidKey(tfxKey key) {
			return GetIndex(key) > -1;
		}

		inline bool ValidIntName(unsigned int name) {
			return GetIntIndex(name) > -1;
		}

		inline bool ValidName(const tfxText &name) {
			return GetIndex(name) > -1;
		}

		//Remove an item from the data. Slow function, 2 memmoves and then the map has to be iterated and indexes reduced by one
		//to re-align them
		inline void Remove(const char *name) {
			tfxKey key = XXHash64::hash(name, strlen(name), 0);
			pair *it = LowerBound(key);
			if (remove_callback)
				remove_callback(data[it->index]);
			unsigned int index = it->index;
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
			unsigned int index = it->index;
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
			unsigned int index = it->index;
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

		inline T &operator[](const unsigned int index) {
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
			tfxKey key = XXHash64::hash(name, strlen(name), 0);
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
			tfxKey key = XXHash64::hash(name.c_str(), name.Length(), 0);
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

	struct tfxrange {
		unsigned int capacity;
		unsigned int offset_into_memory;
	};

	struct tfxmemory {
		unsigned int current_size;
		unsigned int capacity;
		void* data;
		tfxvec<tfxrange> ranges;
		tfxStorageMap<tfxvec<unsigned int>> free_ranges;

		inline tfxmemory() { current_size = capacity = 0; data = NULL; }
		inline tfxmemory(unsigned int qty) { current_size = capacity = 0; data = NULL; resize(qty); }
		inline ~tfxmemory() { if (data) free(data); data = NULL; current_size = capacity = 0; }

		inline bool					empty() { return current_size == 0; }
		inline unsigned int			size() { return current_size; }
		inline const unsigned int	size() const { return current_size; }

		inline void					free_all() { if (data) { current_size = capacity = 0; free(data); data = NULL; } }
		inline void					clear() { if (data) { current_size = 0; ranges.clear(); free_ranges.Clear(); } }

		inline unsigned int			_grow_capacity(unsigned int sz) const { unsigned int new_capacity = capacity ? (capacity + capacity / 2) : 8; return new_capacity > sz ? new_capacity : sz; }
		inline void					resize(unsigned int new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
		inline void					shrink(unsigned int new_size) { assert(new_size <= current_size); current_size = new_size; }
		inline void					reserve(unsigned int new_capacity) { if (new_capacity <= capacity) return; char* new_data = (char*)malloc((size_t)new_capacity * sizeof(char)); if (data) { memcpy(new_data, data, (size_t)current_size * sizeof(char)); free(data); } data = new_data; capacity = new_capacity; }
		inline unsigned int			free_unused_space() { return capacity - current_size; }

		inline bool has_free_range_available(unsigned int bytes) {
			return free_ranges.ValidKey(tfxKey(bytes));
		}
		inline unsigned int			get_range(unsigned int bytes) { 
			if (current_size + bytes >= capacity) reserve(_grow_capacity(current_size + bytes)); 
			tfxrange range;
			range.offset_into_memory = current_size;
			range.capacity = bytes;
			if (free_ranges.ValidKey(tfxKey(bytes))) {
				tfxvec<unsigned int> &free = free_ranges.At(tfxKey(bytes));
				if (!free.empty()) {
					return free.pop_back();
				}
			}
			ranges.push_back(range);
			current_size += bytes;
			return ranges.current_size - 1;
		}
		inline void					free_range(unsigned int index) {
			assert(index < ranges.size());		//Index to free must be within size of ranges list
			int capacity = ranges[index].capacity;
			if (!free_ranges.ValidKey(tfxKey(ranges[index].capacity))) {
				tfxvec<unsigned int> new_free_ranges;
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
#define tfxINVALID 0xFFFFFFFF
	template<typename T>
	struct tfxring {
		unsigned int current_size;
		unsigned int capacity;
		unsigned int start_index;
		unsigned int pos;
		unsigned int range_index;
		T* data;
		inline tfxring() { pos = start_index = current_size = capacity = 0; data = NULL; range_index = tfxINVALID; }
		inline tfxring(unsigned int qty) { pos = start_index = current_size = capacity = 0; data = NULL; reserve(qty); }
		inline tfxring(tfxmemory &mem, unsigned int index, unsigned int size) { pos = start_index = current_size = 0; capacity = size; data = mem.data + mem.ranges[index].offset_into_memory; range_index = index; }

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
		inline bool	        push_back(const T& v) { if (current_size == capacity) return false; new((void*)(data + end_index())) T(v); current_size++; return true; }
		inline bool	        push_front(const T& v) { if (current_size == capacity) return false; new((void*)(data + before_start_index())) T(v); current_size++; start_index = before_start_index(); return true; }
		inline void         reserve(unsigned int new_capacity) { if (new_capacity <= capacity) return; free(data); data = (T*)malloc((size_t)new_capacity * sizeof(T)); capacity = new_capacity; start_index = current_size = 0; }
		inline void			assign_memory(tfxmemory &mem, unsigned int element_size, unsigned int element_count) { 
			assert(range_index == tfxINVALID);	//call refresh instead if a range has already been assigned, or free_all and then assign another range;
			if (!element_count) return;
			int index = mem.get_range(element_size * element_count);
			pos = start_index = current_size = 0; range_index = index; capacity = element_count; void *ptr = (char*)mem.data + mem.ranges[index].offset_into_memory; data = static_cast<T*>(ptr); }
		inline void			refresh(tfxmemory &mem) { void *ptr = (char*)mem.data + mem.ranges[range_index].offset_into_memory; data = static_cast<T*>(ptr); }
		inline void			reset_to_null() { pos = start_index = current_size = 0; range_index = tfxINVALID; data = NULL; }
		inline void			bump() { if (current_size == 0) return; start_index++; start_index %= capacity; current_size--; }
		inline void			bump(unsigned int amount) { if (current_size == 0) return; if (amount > current_size) amount = current_size; start_index += amount; start_index %= capacity; current_size -= amount; }
		inline void			shrink(unsigned int amount) { if (amount > current_size) current_size = 0; else current_size -= amount; }

		inline bool			reset() { if (current_size == 0) return false; assert(data); pos = 0; return true; }
		inline T&			next() { assert(current_size > 0); unsigned int current_pos = (start_index + pos) % capacity; pos++; return data[current_pos]; }
		inline T&			next(int start) { assert(current_size > 0); unsigned int current_pos = (start + pos) % capacity; pos++; return data[current_pos]; }
		inline bool			eob() { return pos >= current_size; }
	};

	//A version of tfxvec that has a fixed capacity where you assign memory from a tfxmemory object
	template<typename T>
	struct tfxfixedvec {
		unsigned int current_size;
		unsigned int capacity;
		unsigned int range_index;
		T* data;
		inline tfxfixedvec() { current_size = capacity = 0; data = NULL; range_index = tfxINVALID; }
		inline tfxfixedvec(tfxmemory &mem, unsigned int index, unsigned int size) { current_size = 0; capacity = size; data = mem.data + mem.ranges[index].offset_into_memory; range_index = index; }

		inline bool			empty() { return current_size == 0; }
		inline bool			full() { return current_size == capacity; }
		inline void         free_range(tfxmemory &mem) {
			if (range_index == tfxINVALID) return;
			if (data) {
				current_size = capacity = 0;
				mem.free_range(range_index);
				range_index = tfxINVALID;
				data = NULL;
			}
		}
		inline void         clear() { if (data) { current_size = 0; } }
		inline unsigned int			size() { return current_size; }
		inline const unsigned int	size() const { return current_size; }
		inline T&           operator[](unsigned int i) { return data[i]; }
		inline const T&     operator[](unsigned int i) const { return data[i]; }

		inline unsigned int end_index() { return current_size; }
		inline unsigned int last_index() { return current_size - 1; }
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
		inline T&           back() { assert(current_size > 0); return data[last_index()]; }
		inline const T&     back() const { assert(current_size > 0); return data[last_index()]; }

		inline void         pop() { assert(current_size > 0); current_size--; }
		inline T&	        pop_back() { assert(current_size > 0); current_size--; return data[current_size]; }
		inline bool	        push_back(const T& v) { if (current_size == capacity) return false; new((void*)(data + end_index())) T(v); current_size++; return true; }
		inline void			assign_memory(tfxmemory &mem, unsigned int element_size, unsigned int element_count) {
			assert(range_index == tfxINVALID);	//call refresh instead if a range has already been assigned, or free_all and then assign another range;
			if (!element_count) return;
			int index = mem.get_range(element_size * element_count);
			current_size = 0; range_index = index; capacity = element_count; void *ptr = (char*)mem.data + mem.ranges[index].offset_into_memory; data = static_cast<T*>(ptr);
		}
		inline void			refresh(tfxmemory &mem) { void *ptr = (char*)mem.data + mem.ranges[range_index].offset_into_memory; data = static_cast<T*>(ptr); }
		inline void			reset_to_null() { current_size = 0; range_index = tfxINVALID; data = NULL; }
		inline void			shrink(unsigned int amount) { if (amount > current_size) current_size = 0; else current_size -= amount; }

		inline void			copyto(tfxmemory &mem, tfxfixedvec &dst) { 
			assert(dst.capacity > capacity); 
			memcpy(dst.data, data, mem.ranges[range_index].capacity); 
			dst.current_size = current_size;
		}
	};

	const u32 tfxMAGIC_NUMBER = '!XFT';
	const u32 tfxMAGIC_NUMBER_INVENTORY = '!VNI';
	const u32 tfxFILE_VERSION = 1;	//Not doing anything with this yet

	//Basic package manager used for reading/writing effects files
	struct tfxHeader {
		u32 magic_number;						//Magic number to confirm file format
		u32 file_version;						//The version of the file
		u32 flags;								//Any flags for the file
		u32 reserved0;							//Reserved for future if needed
		u64 offset_to_inventory;				//Memory offset for the inventory of files
		u64 reserved1;							//More reserved space
		u64 reserved2;							//More reserved space
		u64 reserved3;							//More reserved space
		u64 reserved4;							//More reserved space
		u64 reserved5;							//More reserved space
	};

	struct tfxEntryInfo {
		tfxText file_name;						//The file name of the name stored in the package
		u64 offset_from_start_of_file;			//Offset from the start of the file to where the file is located
		u64 file_size;							//The size of the file
		tfxstream data;							//The file data
		
		void FreeData();
	};

	struct tfxInventory {
		u32 magic_number;						//Magic number to confirm format of the Inventory
		u32 entry_count;						//Number of files in the inventory
		tfxStorageMap<tfxEntryInfo> entries;	//The inventory list
	};

	struct tfxPackage {
		tfxText file_path;
		tfxHeader header;
		tfxInventory inventory;
		u64 file_size;							//The total file size of the package, should match file size on disk
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
	u64 GetPackageSize(tfxPackage &package);
	bool ValidatePackage(tfxPackage &package);

	//------------------------------------------------------------

	//Structs
	//These are mainly internal structs
	//------------------------------------------------------------

	typedef tfxVec2 Point;

	struct NodePair {
		Point left;
		Point right;
		Point left_curve_left;
		Point left_curve_right;
		Point right_curve_left;
		Point right_curve_right;
		tfxAttributeNodeFlags flags;
	};

	struct AttributeNode {
		float frame;
		float value;

		Point left;
		Point right;

		tfxAttributeNodeFlags flags;
		unsigned int index;

		AttributeNode() : frame(0.f), value(0.f), flags(0) { }
		inline bool operator==(const AttributeNode& n) { return n.frame == frame && n.value == value; }

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

	struct Random {
		uint64_t seeds[2];
		Random();

		void ReSeed();
		void ReSeed(uint64_t seed1, uint64_t seed2);
		double Millisecs();

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

		inline unsigned int RangeUInt(unsigned int max) {
			float a = Generate() * (float)max;
			return unsigned int(a);
		};

	};

	static Random random_generation;

	struct GraphLookup {
		tfxvec<float> values;
		unsigned int last_frame;
		float life;

		GraphLookup() : last_frame(0), life(0) {}
	};

	struct GraphID {
		GraphCategory category;
		GraphType type = tfxGraphMaxIndex;
		unsigned int graph_id = 0;
		unsigned int node_id = 0;
	};

	struct GraphLookupIndex {
		unsigned int start_index;
		unsigned int length;
		float max_life;
		float padding1;
	};
	
	//This struct is used to store indexing data in order to index into large lists containing either the node data of graphs
	//or the lookup data of compiled graphs. This is so that we can upload that data into a buffer on the GPU to get the particles
	//updating in a compute shader.
	struct EffectLookUpData {
		GraphLookupIndex overtime_velocity;
		GraphLookupIndex overtime_width;
		GraphLookupIndex overtime_height;
		GraphLookupIndex overtime_weight;
		GraphLookupIndex overtime_spin;
		GraphLookupIndex overtime_stretch;
		GraphLookupIndex overtime_red;
		GraphLookupIndex overtime_green;
		GraphLookupIndex overtime_blue;
		GraphLookupIndex overtime_opacity;
		GraphLookupIndex overtime_velocity_turbulance;
		GraphLookupIndex overtime_direction_turbulance;
		GraphLookupIndex overtime_velocity_adjuster;
		GraphLookupIndex overtime_intensity;
		GraphLookupIndex overtime_direction;
		GraphLookupIndex overtime_noise_resolution;
	};

	struct Graph {
		//The ratio to transalte graph frame/value to grid x/y coords on a graph editor
		Point min;
		Point max;
		GraphPreset graph_preset;
		GraphType type;
		EffectEmitter *effector;
		tfxvec<AttributeNode> nodes;
		GraphLookup lookup;
		unsigned int index;

		Graph();
		~Graph();

		AttributeNode* AddNode(float frame, float value, tfxAttributeNodeFlags flags = 0, float x1 = 0, float y1 = 0, float x2 = 0, float y2 = 0);
		void AddNode(AttributeNode &node);
		void SetNode(uint32_t index, float frame, float value, tfxAttributeNodeFlags flags = 0, float x1 = 0, float y1 = 0, float x2 = 0, float y2 = 0);
		float GetValue(float age);
		float GetRandomValue(float age);
		float GetValue(float age, float life);
		AttributeNode *GetNextNode(AttributeNode &node);
		AttributeNode *GetPrevNode(AttributeNode &node);
		AttributeNode *GetLastNode();
		float GetFirstValue();
		AttributeNode* AddCoordNode(float, float);
		AttributeNode* InsertCoordNode(float, float);
		AttributeNode* InsertNode(float, float);
		float *LinkFirstValue();
		float GetLastValue();
		float GetMaxValue();
		float GetMinValue();
		float GetLastFrame();
		tfxvec<AttributeNode>& Nodes();
		AttributeNode* FindNode(const AttributeNode &n);
		void ValidateCurves();
		void DeleteNode(const AttributeNode &n);
		void Reset(float first_node_value, GraphPreset preset, bool add_node = true);
		void DragValues(GraphPreset preset, float &frame, float &value);
		void Clear();
		void Free();
		void Copy(Graph &to);
		bool Sort();
		void ReIndex();
		tfxVec2 GetInitialZoom();
		bool IsOvertimeGraph();
		bool IsGlobalGraph();
		bool IsAngleGraph();

	};

	tfxVec4 GetMinMaxGraphValues(GraphPreset preset);

	Point GetQuadBezier(Point p0, Point p1, Point p2, float t, float ymin, float ymax, bool clamp = true);
	Point GetCubicBezier(Point p0, Point p1, Point p2, Point p3, float t, float ymin, float ymax, bool clamp = true);
	float GetBezierValue(const AttributeNode *lastec, const AttributeNode &a, float t, float ymin, float ymax);
	float GetDistance(float fromx, float fromy, float tox, float toy);
	float GetVectorAngle(float, float);
	static bool CompareNodes(AttributeNode &left, AttributeNode &right);
	void CompileGraph(Graph &graph);
	void CompileGraphOvertime(Graph &graph);
	float GetMaxLife(EffectEmitter &e);
	float GetMaxAmount(EffectEmitter &e);
	float LookupFastOvertime(Graph &graph, float age, float lifetime);
	float LookupFast(Graph &graph, float frame);
	float LookupPreciseOvertime(Graph &graph, float age, float lifetime);
	float LookupPrecise(Graph &graph, float frame);
	float GetRandomFast(Graph &graph, float frame);
	float GetRandomPrecise(Graph &graph, float frame);

	//Node Manipluation
	bool SetNode(Graph &graph, AttributeNode &node, float, float, tfxAttributeNodeFlags flags, float = 0, float = 0, float = 0, float = 0);
	bool SetNode(Graph &graph, AttributeNode &node, float &frame, float &value);
	void SetCurve(Graph &graph, AttributeNode &node, bool is_left_curve, float &frame, float &value);
	bool MoveNode(Graph &graph, AttributeNode &node, float frame, float value, bool sort = true);
	bool SetNodeFrame(Graph &graph, AttributeNode &node, float &frame);
	bool SetNodeValue(Graph &graph, AttributeNode &node, float &value);
	void ClampNode(Graph &graph, AttributeNode &node);
	void ClampCurve(Graph &graph, Point &curve, AttributeNode &node);
	void ClampGraph(Graph &graph);
	bool IsOvertimeGraph(GraphType type);
	bool IsOvertimePercentageGraph(GraphType type);
	bool IsGlobalGraph(GraphType type);
	bool IsGlobalPercentageGraph(GraphType type);
	bool IsAngleGraph(GraphType type);
	bool IsAngleOvertimeGraph(GraphType type);
	bool IsEverythingElseGraph(GraphType type);

	struct GlobalAttributes {
		Graph life;
		Graph amount;
		Graph velocity;
		Graph width;
		Graph height;
		Graph weight;
		Graph spin;
		Graph stretch;
		Graph overal_scale;
		Graph intensity;
		Graph frame_rate;
		Graph splatter;
		Graph effect_angle;
	};

	struct PropertyAttributes {
		Graph emission_angle;
		Graph emission_range;
		Graph emitter_angle;
		Graph splatter;
		Graph emitter_width;
		Graph emitter_height;
		Graph arc_size;
		Graph arc_offset;
	};

	struct BaseAttributes {
		Graph life;
		Graph amount;
		Graph velocity;
		Graph width;
		Graph height;
		Graph weight;
		Graph spin;
		Graph noise_offset;
	};

	struct VariationAttributes {
		Graph life;
		Graph amount;
		Graph velocity;
		Graph width;
		Graph height;
		Graph weight;
		Graph spin;
		Graph noise_offset;
		Graph noise_resolution;
	};

	struct OvertimeAttributes {
		Graph velocity;
		Graph width;
		Graph height;
		Graph weight;
		Graph spin;
		Graph stretch;
		Graph red;
		Graph green;
		Graph blue;
		Graph blendfactor;
		Graph velocity_turbulance;
		Graph direction_turbulance;
		Graph velocity_adjuster;
		Graph intensity;
		Graph direction;
		Graph noise_resolution;
	};

	static float(*lookup_overtime_callback)(Graph &graph, float age, float lifetime) = LookupFastOvertime;
	static float(*lookup_callback)(Graph &graph, float age) = LookupFast;
	static float(*lookup_random_callback)(Graph &graph, float age) = GetRandomFast;

	struct ShapeData {
		char name[64];
		unsigned int frame_count = 0;
		unsigned int width = 0;
		unsigned int height = 0;
		unsigned int shape_index = 0;
		int import_filter = 0;
	};

	//Store the current state of the object in 2d space
	struct FormState {
		tfxVec2 position;
		tfxVec2 scale;
		float rotation;
	};

	//Store the current local state of the object in 2d space (doesn't require scale here so can save the 8 bytes)
	struct LocalFormState {
		tfxVec2 position;
		float rotation;
	};

	struct Base {
		tfxVec2 size;
		tfxVec2 random_size;
		float velocity;
		float spin;
		float weight;
	};

	//this probably only needs to be in the editor, no use for it in the library? Maybe in the future as an alternative way to play back effects...
	struct AnimationSettings {
		tfxVec4 bb;
		tfxVec2 position;
		tfxVec2 frame_size;
		float scale;
		float zoom;
		int frames;
		int current_frame;
		int frame_offset;
		int extra_frames_count;
		unsigned int seed;
		bool seamless;
		bool loop;
		bool needs_recording;
		unsigned int needs_exporting;
		float max_radius;
		unsigned int largest_frame;
		float playback_speed;
		ExportColorOptions color_option;
		ExportOptions export_option;
		bool export_with_transparency;
	};

	//------------------------------------------------------------

	//API structs you can access in various ways to update and render effects in realtime

	//Image data for particle shapes. This is passed into your custom ShapeLoader function for loading image textures into whatever renderer you're using
	struct ImageData {
		//This can be a ptr to the image texture for rendering. You must assign this in your ShapeLoader function
		void *ptr;

		//Each particle shape saved in an effect library has a unique index
		unsigned int shape_index;
		//The size of one frame of the image
		tfxVec2 image_size;
		//Image index refers to any index that helps you look up the correct image to use. this could be an index in a texture array for example.
		unsigned int image_index;
		//The number of frames in the image, can be one or more
		float animation_frames;
		//Maximum distance to the nearest transparent edge of the image
		float max_radius;
		int import_filter;
		unsigned int compute_shape_index;

		//use this definition if you need more spefic data to point to the image texture in whatever renderer you're using
		//Just define tfxCUSTOM_IMAGE_DATA before you include timelinefx.h
#ifdef tfxCUSTOM_IMAGE_DATA
		tfxCUSTOM_IMAGE_DATA
#endif // tfxCUSTOM_IMAGE_DATA

		ImageData() :
			ptr(nullptr),
			animation_frames(1.f),
			shape_index(0),
			max_radius(0),
			import_filter(0)
		{ }
	};

	struct EmitterProperties {
		//Pointer to the ImageData in the EffectLibary. 
		ImageData *image;
		//Assigns either alpha or additive blend to particles that are spawned
		EmissionType emission_type;
		//Should particles emit towards the center of the emitter or away, or in a specific direction
		EmissionDirection emission_direction;
		//How particles should behave when they reach the end of the line
		LineTraversalEndBehaviour end_behaviour;
		//The rotation of particles when they spawn, or behave overtime if tfxAlign is used
		AngleSetting angle_setting = AngleSetting::tfxRandom;

		//Bit field of various boolean flags
		tfxEmitterPropertyFlags flags;
		tfxParticleControlFlags compute_flags;

		//Offset to draw particles at
		tfxVec2 image_handle;
		//Animation frame rate
		float frame_rate;
		//Offset of emitters
		tfxVec2 emitter_handle;
		//When single or one shot flags are set, spawn this amount of particles in one go
		unsigned int spawn_amount;
		//Layer of the particle manager that the particle is added to
		unsigned int layer;
		//The shape being used for all particles spawned from the emitter
		unsigned int shape_index;

		//Angle added to the rotation of the particle when spawned or random angle range if angle setting is set to tfxRandom
		float angle_offset;
		//The number of rows/columns/ellipse/line points in the grid when spawn on grid flag is used
		tfxVec2 grid_points;
		//The number of millisecs before an effect or emitter will loop back round to the beginning of it's graph lookups
		float loop_length;
		//The start frame index of the animation
		float start_frame;
		//The final frame index of the animation
		float end_frame;
		//Milliseconds to delay spawing
		float delay_spawning;

		EmitterProperties() :
			angle_offset(360),
			flags(tfxEmitterPropertyFlags_image_handle_auto_center | tfxEmitterPropertyFlags_grid_spawn_clockwise | tfxEmitterPropertyFlags_emitter_handle_auto_center | tfxEmitterPropertyFlags_global_uniform_size | tfxEmitterPropertyFlags_base_uniform_size | tfxEmitterPropertyFlags_lifetime_uniform_size),
			image(nullptr),
			image_handle(tfxVec2()),
			spawn_amount(1),
			emission_type(EmissionType::tfxPoint),
			emission_direction(EmissionDirection::tfxOutwards),
			grid_points({ 10.f, 10.f }),
			emitter_handle(tfxVec2()),
			end_behaviour(LineTraversalEndBehaviour::tfxLoop),
			loop_length(0.f),
			layer(0),
			shape_index(1),
			start_frame(0),
			end_frame(0),
			frame_rate(30.f),
			angle_setting(AngleSetting::tfxRandom),
			delay_spawning(0.f)
		{ }
	};

	//Store the current state of the Effect/Emitter. All these values can change over the lifetime of the effect/emitter.
	//It's questionable that we actually need to store these, rather then just use local variables during the ControlParticles function
	struct EmitterState {
		//Base particle size
		tfxVec2 size;
		//Particle size variation
		tfxVec2 size_variation;
		//Particle color and opacity
		tfxRGBA color;
		//Size of the emitter area that particles can spawn in. X will be used for line length for line effects
		tfxVec2 emitter_size;
		//Emitter handle - the offset at which the emitter is positioned
		tfxVec2 emitter_handle;
		//Offset to draw particles at
		tfxVec2 image_handle;
		//Current base life that particles will be spawned with (milliseconds)
		float life;
		//Current number of particles that will be spawned per second
		float amount;
		//variance of the number of particles that will be spawned per second
		float amount_variation;
		//The amount of variation of life that particles are spawned with (milliseconds)
		float life_variation;
		//The base velocity of particles (pixels per second)
		float velocity;
		//The amount that velocity will vary
		float velocity_variation;
		//Global multiplier for all currently spawned particles of this emitter.
		float velocity_adjuster;
		//Base spine that a particle will rotate (radians per second)
		float spin;
		//Amount in radians that the spin will vary
		float spin_variation;
		//Direction of travel that the particles go when spawned (radians)
		float emission_angle;
		//Amount the emission angle will vary (radians)
		float emission_angle_variation;
		//Amount that particles will randomly offset from the spawn point (radius in pixels)
		float splatter;
		//For Ellipse type emitters, you can set the arc_size to only spawn a segment of the ellipse (radians)
		float arc_size;
		//Starting point in radians for the arc segment
		float arc_offset;
		//Scales all particles uniformly so that you can make the effect bigger or smaller overal
		float overal_scale;
		//Amount that particles will stretch base on their current velocity
		float stretch;
		//The base weight of spawned particles (downward accelleration in pixels per second)
		float weight;
		//The amount the the weight will vary
		float weight_variation;
		//The more motion randomness the more particles will move about erratically
		float noise_offset_variation;
		float noise_offset;
		//The more motion randomness the more particles will move about erratically
		float noise_resolution;
		//The current age of the emitter.
		float age;
		//The current frame of the effect/emitter as calculated by age / update_frequency
		float frame;
		//internal use variables
		float amount_remainder;
		bool emission_alternator;
		bool single_shot_done;

		tfxVec2 grid_coords;
		tfxVec2 grid_direction;
		tfxVec2 grid_segment_size;

		EmitterState() : grid_coords(tfxVec2()), single_shot_done(false), age(0.f), amount_remainder(0.f) {}
	};

	struct tfxEffectState {
		//All of these values are global values that affect overal relative effects of sub emitters and effects
		float life;
		float amount;
		tfxVec2 size;
		float velocity;
		float spin;
		float intensity;
		float splatter;
		float overal_scale;
		float stretch;
		float weight;
	};

	struct tfxTransform {

		//Position, scale and rotation values
		LocalFormState local;
		FormState world;
		FormState captured;
		//2d matrix for transformations
		Matrix2 matrix;

	};

	struct tfxCommon {

		float frame;
		float age;
		float loop_length;
		float highest_particle_age;
		float timeout_counter;
		float timeout;
		unsigned int active_children;
		tfxVec2 handle;
		tfxEmitterStateFlags state_flags;
		tfxEmitterPropertyFlags property_flags;
		EffectLibrary *library;
		tfxEffect *root_effect;

		tfxCommon() :
			frame(0.f),
			age(0.f),
			state_flags(0),
			property_flags(0),
			timeout_counter(0),
			timeout(1000.f),
			active_children(0),
			root_effect(nullptr)
		{ }

	};

	struct tfxEmitterState {
		tfxVec2 emitter_size;
		tfxVec2 grid_coords;
		tfxVec2 grid_direction;
		float velocity_adjuster;
		float global_opacity;
		float stretch;
		float emitter_handle_y;
		float overal_scale;
		float amount_remainder;
		float emission_alternator;
		float qty;

		tfxEmitterState() :
			amount_remainder(0.f)
		{}
	};

	struct tfxEmitterSpawnProperties {
		unsigned int layer;
		unsigned int blend_mode;
		unsigned int spawn_amount;
		EmissionType emission_type;
		AngleSetting angle_setting;
		EmissionDirection emission_direction;
		float start_frame;
		float animation_frames;
		float end_frame;
		float frame_rate;
		void *image_ptr;
		float angle_offset;
		tfxVec2 image_size;
		tfxVec2 image_handle;
		tfxVec2 emitter_handle;
		tfxVec2 grid_points;
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
		tfxVec2 grid_segment_size;
	};

	struct tfxEmitter {
		tfxTransform transform;
		EffectEmitter *library_link;
		EffectEmitterType type;
		Particle *parent_particle;
		tfxEmitter *parent;
		tfxCommon common;
		tfxEmitterState current;
		LookupMode lookup_mode;
		tfxVec2 image_handle;
		unsigned int offset;
		tfxEmitter *next_emitter;

		tfxfixedvec<Particle> particles;

		tfxEmitter() : parent(nullptr), parent_particle(nullptr) {}
		void Reset();
		void UpdateEmitter();
		void UpdateAsSubEffect();
		bool GrowParticles(unsigned int min_amount);
		void RefreshFromLibrary();
		void SpawnParticles();
		float GetEmissionDirection(tfxVec2 &local_position, tfxVec2 &world_position, tfxVec2 &emitter_size);
		void InitCPUParticle(Particle &p, tfxEmitterSpawnControls &spawn_values, float tween);
		void ControlParticles();
		bool FreeCapacity();
		void *UserData();
		Particle &GrabParticle();
	};

	struct tfxEffect {
		//todo: Put an operator overload for = with an assert, these shouldn't be copied in that way

		tfxTransform transform;
		EffectEmitter *library_link;
		tfxEffectID id;
		LookupMode lookup_mode;
		Particle *parent_particle;
		tfxCommon common;
		tfxKey path_hash;
		tfxEffectState current;
		unsigned int max_particles[tfxLAYERS];
		unsigned int max_sub_emitters;
		tfxEffectPool *storage;
		tfxfixedvec<tfxEmitter> sub_emitters;
		tfxfixedvec<tfxEmitter> sub_effects;
		tfxfixedvec<ParticleSprite> sprites[tfxLAYERS];
		void *user_data;
		void(*update_callback)(tfxEffect &effect);		//Called after the state has been udpated

		tfxEffect() :
			parent_particle(nullptr),
			path_hash(0),
			max_sub_emitters(0),
			lookup_mode(tfxFast)
		{}
		void Reset();
		void ReleaseMemory();
		void ClearSprites();
		void CompressSprites();
		void UpdateSpritePointers();
		void RefreshFromLibrary();
		inline unsigned int ParticleCount() { unsigned int count = 0; for (EachLayer) { count += sprites[layer].current_size; } return count; }
		ParticleSprite &GrabSprite(unsigned int layer);
		tfxEmitter &GrabSubEffect();
		tfxEmitter& AddSubEffect(tfxEmitter &sub_effect);
		inline void AllowMoreParticles() { common.property_flags |= tfxEmitterPropertyFlags_can_grow_particle_memory; }
		inline void Move(float x, float y) { transform.local.position.x += x; transform.local.position.y += y; }
		inline void Position(float x, float y, bool capture = true) { transform.local.position.x = x; transform.local.position.y = y; common.state_flags |= tfxEmitterStateFlags_no_tween_this_update; }
		inline void Position(tfxVec2 pos, bool capture = true) { transform.local.position = pos; common.state_flags |= tfxEmitterStateFlags_no_tween_this_update; }
		inline void SetLookupMode(LookupMode mode) { lookup_mode = mode; }
		
	};

	//An EffectEmitter can either be an effect which stores emitters and global graphs for affecting all the attributes in the emitters
	//Or it can be an emitter which spawns all of the particles. Effectors are stored in the particle manager effects list buffer.
	//This is only for library storage, when using to update each frame this is copied to tfxEffectType and tfxEmitterType, much more compact versions more
	//suited for realtime use.
	struct EffectEmitter {
		tfxTransform transform;

		//The current state of the effect/emitter
		EmitterState current;
		unsigned int emitter_state_index;
		//All of the properties of the effect/emitter
		EmitterProperties properties;

		//Name of the effect
		tfxText name;						//Todo: Do we need this here?
		//A hash of the directory path to the effect ie Flare/spark
		tfxKey path_hash;
		//Is this a tfxEffectType or tfxEmitterType
		EffectEmitterType type;
		//A pointer to the library this effect belongs
		EffectLibrary *library;
		//The index within the library that this exists at
		unsigned int library_index;
		//The current highest particle age. When using a compute buffer we don't have any reliable way of keeping track of particle counts of individual emitters, so how do we know when to remove an emitter
		//after all it's particles have expired? We set this variable to the highest particle age each time it spawns a particle and then counts it down each frame. When it's 0 then we know that there are no
		//more particles being controlled by this emitter and can therefore time it out.
		float highest_particle_age;
		//The number of frames before this is removed from the particle manager after particle count is 0
		unsigned int timeout;
		//Internal, keep track of idle frames
		unsigned int timeout_counter;
		//Every effect and emitter in the library gets a unique id
		unsigned int uid;
		//compute slot id if a compute shader is being used. Only applied to bottom emitters (emitters with no child effects)
		unsigned int compute_slot_id;
		//The max_radius of the emitter, taking into account all the particles that have spawned and active
		float max_radius;
		//List of sub_effects ( effects contain emitters, emitters contain sub effects )
		tfxvec<EffectEmitter> sub_effectors;
		//Experiment with emitters maintaining their own list of particles. Buffer info contains info for fetching the area in the buffer stored in the particle manager
		//Custom user data, can be accessed in callback functions
		void *user_data;

		//Idea for adding you're own user define struct members
#ifdef TFX_CUSTOM_EMITTER
TFX_CUSTOM_EMITTER
#endif

		//All graphs that the effect uses to lookup attribute values are stored in the library. These variables here are indexes to the array where they're stored
		unsigned int global;
		unsigned int property;
		unsigned int base;
		unsigned int variation;
		unsigned int overtime;
		//Experitment: index into the lookup index data in the effect library
		unsigned int lookup_node_index;
		unsigned int lookup_value_index;
		//effect sprite offset
		unsigned int sprite_offset;
		//Index to animation settings stored in the effect library. Would like to move this at some point
		unsigned int animation_settings;
		//The maximum amount of life that a particle can be spawned with taking into account base + variation life values
		float max_life;
		//The estimated maximum time that the sub emitter might last for, taking into account the parent particle lifetime
		float max_sub_emitter_life;
		//The maximum amount of particles that this effect can spawn (root effects and emitters only)
		u32 max_particles[tfxLAYERS];
		u32 max_sub_emitters;
		//Pointer to the immediate parent
		EffectEmitter *parent;
		//Pointer to the next pointer in the particle manager buffer. 
		EffectEmitter *next_ptr;
		//Pointer to the sub effect's particle that spawned it
		Particle *parent_particle;

		//Custom fuction pointers that you can use to override attributes and affect the effect/emitters behaviour in realtime
		//See tfxEffectTemplate for applying these callbacks
		//tfxEffectSprites *sprites_container;
		//Callbacks for effect pool effects only:
		void(*root_effect_update_callback)(tfxEffect &effect);						//Called after the root effect state has been udpated
		void(*emitter_update_callback)(tfxEmitter &emitter);						//Called after the emitter state has been udpated
		void(*spawn_update_callback)(tfxEmitterSpawnControls &spawn_controls, tfxEmitter &emitter);				//Called before the emitter spawns particles
		void(*particle_onspawn_callback)(Particle &particle);						//Called as each particle is spawned.

		tfxEmitterStateFlags flags;

		EffectEmitter() :
			highest_particle_age(0),
			parent(nullptr),
			parent_particle(nullptr),
			user_data(nullptr),
			flags(tfxEmitterStateFlags_no_tween_this_update | tfxEmitterStateFlags_enabled),
			timeout(100),
			timeout_counter(0),
			animation_settings(0),
			library(nullptr),
			root_effect_update_callback(nullptr),
			emitter_update_callback(nullptr),
			spawn_update_callback(nullptr),
			particle_onspawn_callback(nullptr),
			max_sub_emitters(0),
			max_sub_emitter_life(0.f)
		{
			for (int i = 0; i != tfxLAYERS; ++i) {
				max_particles[i] = 0;
			}
		}
		~EffectEmitter();
		bool operator < (const EffectEmitter& e) const
		{
			return (properties.layer < e.properties.layer);
		}

		//API functions
		//Tell the effect to stop spawning so that eventually particles will expire and the effect will be removed from the particle manager
		inline void SoftExpire();

		void SetUserData(void *data);
		void *GetUserData();
		void SetTimeout(unsigned int frames);

		//Override graph functions for use in update_callback
		//Some of these change the same state and property values, but they're named differently just to make it clearer as to whether you're overriding kEffect or a kEmitter.

		//Internal functions
		EffectEmitter& AddEmitter(EffectEmitter &e);
		EffectEmitter& AddEffect(EffectEmitter &e);
		EffectEmitter& AddEffect();
		EffectEmitter& AddEffector(EffectEmitterType type = tfxEmitterType);
		EffectEmitter* GetRootEffect();
		bool IsRootEffect();
		void ReIndex();
		void CountChildren(int &emitters, int &effects);
		void ResetParents();
		EffectEmitter* MoveUp(EffectEmitter &effect);
		EffectEmitter* MoveDown(EffectEmitter &effect);
		void DeleteEmitter(EffectEmitter *effect);
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
		unsigned int GetSubEffectSpriteCounts(unsigned int layer, unsigned int multiplier);
		float GetSubEffectLength();
		unsigned int GetHighestQty(float parent_age);
		Graph* GetGraphByType(GraphType type);
		unsigned int GetGraphIndexByType(GraphType type);
		void CompileGraphs();
		void InitialiseUninitialisedGraphs();
		void SetName(const char *n);

		void ReSeed(uint64_t seed = 0);
		bool HasSingle();
		bool RenameSubEffector(EffectEmitter &effect, const char *new_name);
		bool NameExists(EffectEmitter &effect, const char *name);
		void FreeGraphs();
		void NoTweenNextUpdate();

		void ClearColors();
		void AddColorOvertime(float frame, tfxRGB color);
		void Clone(EffectEmitter &clone, EffectEmitter *root_parent, EffectLibrary *destination_library, bool keep_user_data = false, bool force_clone_global = false);
		void CopyToEffect(tfxEffectID &effect_id, tfxEffectPool &storage);
		void CopyToEmitter(tfxEmitter &emitter, tfxEffectPool &storage, bool assign_memory);
		void EnableAllEmitters();
		void EnableEmitter();
		void DisableAllEmitters();
		void DisableAllEmittersExcept(EffectEmitter &emitter);
		bool IsFiniteEffect();
	};

	struct tfxEffectTemplate {
		tfxStorageMap<EffectEmitter*> paths;
		EffectEmitter effect_template;

		void AddPath(EffectEmitter &effectemitter, tfxText path) {
			paths.Insert(path, &effectemitter);
			for (auto &sub : effectemitter.sub_effectors) {
				tfxText sub_path = path;
				sub_path.Appendf("/%s", sub.name.c_str());
				AddPath(sub, sub_path);
			}
		}

		inline EffectEmitter &Effect() { return effect_template; }
		inline EffectEmitter *Get(tfxText path) { if (paths.ValidName(path)) return paths.At(path); return nullptr; }
		inline void SetUserData(tfxText path, void *data) { if(paths.ValidName(path)) paths.At(path)->user_data = data; }
		inline void SetUserData(void *data) { effect_template.user_data = data; }
		void SetUserDataAll(void *data);
		inline void SetUpdateCallback(tfxText path, void(*root_effect_update_callback)(tfxEffect &effect)) { if (paths.ValidName(path)) paths.At(path)->root_effect_update_callback = root_effect_update_callback; }
		inline void SetUpdateCallback(tfxText path, void(*emitter_update_callback)(tfxEmitter &emitter)) { if (paths.ValidName(path)) paths.At(path)->emitter_update_callback = emitter_update_callback; }
		inline void SetUpdateCallback(tfxText path, void(*spawn_update_callback)(tfxEmitterSpawnControls &spawn_controls, tfxEmitter &emitter)) { if (paths.ValidName(path)) paths.At(path)->spawn_update_callback = spawn_update_callback; }
		void SetParticleOnSpawnCallback(tfxText path, void(*particle_onspawn_callback)(Particle &particle));
	};

	struct EffectEmitterSnapShot {
		EffectEmitter effect;
		unsigned int index;
		char description[256];
		char path[512];
		bool is_current_revision = false;
		void SetDescription(const char *format, ...);
	};

	struct ComputeSprite {	//64 bytes
		tfxVec4 bounds;				//the min/max x,y coordinates of the image being drawn
		tfxVec4 uv;					//The UV coords of the image in the texture
		tfxVec4 scale_rotation;		//Scale and rotation (x, y = scale, z = rotation, w = multiply blend factor)
		tfxVec2 position;			//The position of the sprite
		tfxRGBA8 color;				//The color tint of the sprite
		unsigned int parameters;	//4 extra parameters packed into a u32: blend_mode, image layer index, shader function index, blend type
	};

	struct ParticleSprite {	//80 bytes
		void *ptr;					//Pointer to the image data
		Particle *particle;			//We need to point to the particle in order to update it's sprite index
		FormState world;
		FormState captured;
		tfxVec2 handle;
		tfxRGBA8 color;				//The color tint of the sprite
		float intensity;			
		unsigned int parameters;	//4 extra parameters packed into a u32: blend_mode, expired flag, frame
	};

	//Initial particle struct, looking to optimise this and make as small as possible
	//These are spawned by effector emitter types
	//Particles are stored in the particle manager particle buffer.
	//I really think that tweened frames should be ditched in favour of delta time so captured can be ditched
	//168 bytes
	struct Particle {
		LocalFormState local;			//The local position of the particle, relative to the emitter.
		FormState world;				//The world position of the particle relative to the world/screen.
		FormState captured;				//The captured world coords for tweening
		Matrix2 matrix;					//Simple 2d matrix for transforms (only needed for sub effects)

		//Read only when ControlParticle is called, only written to at spawn time
		Base base;						//Base values created when the particle is spawned. They can be different per particle due to variations
		float emission_angle;			//Emission angle of the particle at spawn time
		float noise_offset;				//Higer numbers means random movement is less uniform
		float noise_resolution;			//Higer numbers means random movement is more uniform
		tfxParticleFlags flags;			//flags for different states

		//Updated everyframe
		float age;						//The age of the particle, used by the controller to look up the current state on the graphs
		float max_age;					//max age before the particle expires
		float image_frame;				//Current frame of the image if it's an animation
		float weight_acceleration;		//The current amount of gravity applied to the y axis of the particle each frame
		float intensity;				//Color is multiplied by this value in the shader to increase the brightness of the particles
		tfxRGBA8 color;					//Colour of the particle

		EffectEmitter *parent;			//pointer to the emitter that emitted the particle.
		tfxEmitter *emitter;			//pointer to the emitter that emitted the particle.
		//Internal use variables
		Particle *next_ptr;
		unsigned int sprite_index;
		unsigned int offset;

		//Override functions, you can use these inside an update_callback if you want to modify the particle's behaviour
		inline void OverridePosition(float x, float y) { local.position.x = x; local.position.y = y; }
		inline void OverrideSize(float x, float y) { world.scale.x = x; world.scale.y = y; }
		//inline void OverrideVelocity(float x, float y) { velocity.x = x; velocity.y = y; }
		//inline void OverrideVelocityScale(float v) { velocity_scale = v; }
		inline void OverrideRotation(float r) { local.rotation = r; }
		inline void OverrideRotationDegrees(float d) { local.rotation = tfxDegrees(d); }
		//inline void OverrideDirection(float r) { direction = r; }
		//inline void OverrideDirectionDegrees(float d) { direction = tfxDegrees(d); }
		//inline void OverrideSpin(float r) { direction = r; }
		//inline void OverrideSpinDegrees(float d) { spin = tfxDegrees(d); }
		inline void OverrideAge(float a) { age = a; }
		//inline void OverrideImageFrameRate(float fr) { image_frame_rate = fr; }
		inline void OverrideImageFrame(float f) { image_frame = f; }
		inline void OverrideWeight(float w) { weight_acceleration = w; }
		inline void OverrideIntensity(float i) { intensity = i; }
		inline void OverrideBaseSize(float width, float height) { base.size.x = width; base.size.y = height; }
		inline void OverrideBaseVelocity(float v) { base.velocity = v; }
		inline void OverrideBaseWeight(float w) { base.weight = w; }
		inline void OverrideBaseSpin(float s) { base.spin = s; }
		inline void OverrideBaseSpinDegrees(float d) { base.spin = tfxDegrees(d); }

	};

	struct ComputeFXGlobalState {
		u32 start_index = 0;
		u32 current_length = 0;
		u32 max_index = 0;
		u32 end_index = 0;
	};

	struct ComputeController {
		tfxVec2 position;
		float line_length;
		float angle_offset;
		tfxVec4 scale_rotation;				//Scale and rotation (x, y = scale, z = rotation, w = velocity_adjuster)
		float end_frame;
		unsigned int normalised_values;		//Contains normalized values which are generally either 0 or 255, normalised in the shader to 0 and 1 (except opacity): age_rate, line_negator, spin_negator, position_negator, opacity
		tfxParticleControlFlags flags;
		unsigned int image_data_index;		//index into the shape buffer on the gpu. CopyComputeShapeData must be called to prepare the data.
		tfxVec2 image_handle;
		tfxVec2 emitter_handle;
		float noise_offset;
		float stretch;
		float frame_rate;
		float noise_resolution;
	};

	struct ComputeParticle {
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
		unsigned int control_slot_and_layer;	//index to the controller, and also stores the layer in the particle manager that the particle is on (layer << 3)
		float local_rotation;
	};

	struct ComputeImageData {
		tfxVec4 uv;
		tfxVec2 image_size;
		unsigned int image_index;
		float animation_frames;
		//float max_radius;
	};

	//Struct to contain a static state of a particle in a frame of animation. Used in the editor for recording frames of animation
	struct ParticleFrame {
		tfxVec2 position;
		tfxVec2 scale;
		tfxVec2 handle;
		float rotation;
		float image_frame;
		float start_frame;
		void *image_ptr;
		tfxRGBA8 color;
		float intensity;
		bool has_frames;
	};

	struct EffectLibrary {
		tfxStorageMap<EffectEmitter*> effect_paths;
		tfxvec<EffectEmitter> effects;
		tfxStorageMap<ImageData> particle_shapes;

		tfxvec<GlobalAttributes> global_graphs;
		tfxvec<PropertyAttributes> property_graphs;
		tfxvec<BaseAttributes> base_graphs;
		tfxvec<VariationAttributes> variation_graphs;
		tfxvec<OvertimeAttributes> overtime_graphs;
		tfxvec<AnimationSettings> animation_settings;
		//Experiment, all nodes from graphs can get added to this one list to keep them in one place in memory, this is where data is read from when
		//when updating particles and emitters. I'm starting to think towards getting this running in a compute shader, so will need data to be stored
		//in buffers in one block as much as possible
		tfxvec<AttributeNode> all_nodes;
		tfxvec<EffectLookUpData> node_lookup_indexes;
		tfxvec<float> compiled_lookup_values;
		tfxvec<GraphLookupIndex> compiled_lookup_indexes;
		tfxvec<ComputeImageData> shape_data;
		//This could probably be stored globally
		tfxvec<tfxVec4> graph_min_max;

		tfxvec<unsigned int> free_global_graphs;
		tfxvec<unsigned int> free_property_graphs;
		tfxvec<unsigned int> free_base_graphs;
		tfxvec<unsigned int> free_variation_graphs;
		tfxvec<unsigned int> free_overtime_graphs;
		tfxvec<unsigned int> free_animation_settings;

		//Get an effect from the library by index
		EffectEmitter& operator[] (uint32_t index);
		tfxText name;
		bool open_library = false;
		bool dirty = false;
		tfxText library_file_path;
		unsigned int uid = 0;

		//Free everything in the library
		void Clear();
		//Get an effect in the library by it's path. So for example, if you want to get a pointer to the emitter "spark" in effect "explosion" then you could do GetEffect("explosion/spark")
		//You will need this function to apply user data and update callbacks to effects and emitters before adding the effect to the particle manager
		EffectEmitter *GetEffect(tfxText &path);
		EffectEmitter *GetEffect(const char *path);
		//Get an effect by it's path hash key
		EffectEmitter *GetEffect(tfxKey key);
		//Get and effect by it's index
		void PrepareEffectTemplate(tfxText path, tfxEffectTemplate &effect);
		//Copy the shape data to a memory location, like a staging buffer ready to be uploaded to the GPU for use in a compute shader
		void BuildComputeShapeData(void* dst, tfxVec4(uv_lookup)(void *ptr, ComputeImageData &image_data, int offset));
		void CopyComputeShapeData(void* dst);
		void CopyLookupIndexesData(void* dst);
		void CopyLookupValuesData(void* dst);
		u32 GetComputeShapeDataSizeInBytes();
		u32 GetComputeShapeCount();
		u32 GetLookupIndexCount();
		u32 GetLookupValueCount();
		u32 GetLookupIndexesSizeInBytes();
		u32 GetLookupValuesSizeInBytes();

		//Mainly internal functions
		void RemoveShape(unsigned int shape_index);
		EffectEmitter &AddEffect(EffectEmitter &effect);
		EffectEmitter &AddFolder(tfxText name);
		EffectEmitter &AddFolder(EffectEmitter &effect);
		void UpdateEffectPaths();
		void AddPath(EffectEmitter &effectemitter, tfxText path);
		void DeleteEffect(EffectEmitter *effect);
		bool RenameEffect(EffectEmitter &effect, const char *new_name);
		bool NameExists(EffectEmitter &effect, const char *name);
		bool NameExists2(EffectEmitter &effect, const char *name);
		void ReIndex();
		void UpdateParticleShapeReferences(tfxvec<EffectEmitter> &effects, unsigned int default_index);
		EffectEmitter* MoveUp(EffectEmitter &effect);
		EffectEmitter* MoveDown(EffectEmitter &effect);
		unsigned int AddGlobal();
		unsigned int AddProperty();
		unsigned int AddBase();
		unsigned int AddVariation();
		unsigned int AddOvertime();
		void FreeGlobal(unsigned int index);
		void FreeProperty(unsigned int index);
		void FreeBase(unsigned int index);
		void FreeVariation(unsigned int index);
		void FreeOvertime(unsigned int index);
		unsigned int CloneGlobal(unsigned int source_index, EffectLibrary *destination_library);
		unsigned int CloneProperty(unsigned int source_index, EffectLibrary *destination_library);
		unsigned int CloneBase(unsigned int source_index, EffectLibrary *destination_library);
		unsigned int CloneVariation(unsigned int source_index, EffectLibrary *destination_library);
		unsigned int CloneOvertime(unsigned int source_index, EffectLibrary *destination_library);
		void AddEmitterGraphs(EffectEmitter& effect);
		void AddEffectGraphs(EffectEmitter& effect);
		unsigned int AddAnimationSettings(EffectEmitter& effect);
		void UpdateEffectParticleStorage();
		void UpdateComputeNodes();
		void CompileAllGraphs();
		void CompileGlobalGraph(unsigned int index);
		void CompilePropertyGraph(unsigned int index);
		void CompileBaseGraph(unsigned int index);
		void CompileVariationGraph(unsigned int index);
		void CompileOvertimeGraph(unsigned int index);
		void CompileColorGraphs(unsigned int index);
		void SetMinMaxData();
		float LookupPreciseOvertimeNodeList(GraphType graph_type, int index, float age, float life);
		float LookupPreciseNodeList(GraphType graph_type, int index, float age);
		float LookupFastOvertimeValueList(GraphType graph_type, int index, float age, float life);
		float LookupFastValueList(GraphType graph_type, int index, float age);

		//Debug stuff, used to check graphs are being properly recycled
		unsigned int CountOfGraphsInUse();
		unsigned int CountOfFreeGraphs();
	};

	struct BufferInfo {
		u32 start_index;
		u32 end_index;
		u32 length;
	};

	struct tfxEffectPool {
		tfxmemory effect_memory;
		tfxmemory emitter_memory;
		tfxmemory particle_memory;
		tfxmemory sprite_memory;
		
		tfxEffectPool() {}
		~tfxEffectPool() {
			effect_memory.free_all();
			emitter_memory.free_all();
			particle_memory.free_all();
			sprite_memory.free_all();
		};

		void Init(unsigned int max_effects = 1000, unsigned int max_emitters = 10000, unsigned int max_particles = 100000);
		tfxEffectID InsertEffect();
		tfxEffect &GetEffect(tfxEffectID id);
	};

	struct tfxEffectSprites {
		tfxring<ParticleSprite> sprites;
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
		unsigned int timeout = 5;
		unsigned int timeout_counter = 0;
		unsigned int emitter_count = 0;
		float highest_particle_age = 0;
		float frame = 0.f;
		float age = 0.f;
		float amount_remainder = 0;
		float qty = 1.f;
		EffectEmitter *library_link;
		EffectLibrary *library;
		bool single_shot_done = false;
		bool started_spawning = false;
		tfxvec<float> particles[2];
	};

	//This is used to figure out how much memory each effect and emitter needs to draw particles so that the correct amount of memory can be assigned as each effect is used.
	struct tfxParticleMemoryTools {
		unsigned int sprite_count[4];
		unsigned int sub_effect_count;
		unsigned int initial_effect_size = 0;
		unsigned int emitters_removed = 0;
		float max_frames;
		float max_last_life;
		unsigned int current_buffer;
		tfxvec<float> particles[tfxLAYERS][2];
		tfxvec<tfxMockEffect> effects[2];
		EffectEmitter current_effect;

		tfxParticleMemoryTools() : current_buffer(0), sub_effect_count(0) {}

		void AddEffect(EffectEmitter &effect);
		void GetEffectMaxFrames(EffectEmitter &effect);
		void ProcessEffect(EffectEmitter &effect);
		void Process();
		void MockUpdateEmitter(tfxMockEffect &emitter);
		void MockUpdateParticles();
	};

	struct DataEntry {
		DataType type;
		tfxText key;
		tfxText str_value;
		int int_value;
		bool bool_value;
		float float_value;
		double double_value;
	};

	struct DataTypesDictionary {
		tfxStorageMap<DataType> eff;

		DataTypesDictionary();
	};

	static DataTypesDictionary data_types;

	//Internal functions
	//Some file IO functions
	bool HasDataValue(tfxStorageMap<DataEntry> &config, tfxText key);
	void AddDataValue(tfxStorageMap<DataEntry> &config, tfxText key, const char *value);
	void AddDataValue(tfxStorageMap<DataEntry> &config, tfxText key, int value);
	void AddDataValue(tfxStorageMap<DataEntry> &config, tfxText key, bool value);
	void AddDataValue(tfxStorageMap<DataEntry> &config, tfxText key, double value);
	void AddDataValue(tfxStorageMap<DataEntry> &config, tfxText key, float value);
	tfxText &GetDataStrValue(tfxStorageMap<DataEntry> &config, const char* key);
	int& GetDataIntValue(tfxStorageMap<DataEntry> &config, const char* key);
	float& GetDataFloatValue(tfxStorageMap<DataEntry> &config, const char* key);
	bool SaveDataFile(tfxStorageMap<DataEntry> &config, const char* path = "");
	bool LoadDataFile(tfxStorageMap<DataEntry> &config, const char* path);
	void StreamProperties(EmitterProperties &property, tfxText &file);
	void StreamGraph(const char * name, Graph &graph, tfxText &file);
	tfxvec<tfxText> SplitString(const tfxText &s, char delim = 61);
	bool StringIsUInt(const tfxText &s);
	int GetDataType(const tfxText &s);
	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, uint32_t value);
	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, float value);
	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, bool value);
	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, int value);
	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, tfxText &value);
	void AssignGraphData(EffectEmitter &effect, tfxvec<tfxText> &values);
	void AssignNodeData(AttributeNode &node, tfxvec<tfxText> &values);
	EffectEmitter CreateEffector(float x = 0.f, float y = 0.f);
	void TransformParticle(Particle &p, EffectEmitter &e);
	void TransformParticle(Particle &p, tfxEmitter &e);
	void Transform(tfxEmitter &emitter, Particle &parent);
	void Transform(tfxEmitter &emitter, tfxEffect &parent);
	void Transform(tfxEmitter &emitter, tfxEmitter &parent);
	void Transform(FormState &local, FormState &world, EffectEmitter &e);
	FormState Tween(float tween, FormState &world, FormState &captured);
	tfxVec2 InterpolateVec2(float, const tfxVec2&, const tfxVec2&);
	float Interpolatef(float tween, float, float);
	int ValidateEffectPackage(const char *filename);
	void ReloadBaseValues(Particle &p, EffectEmitter &e);
	bool Copy(tfxEffectPool &storage, EffectEmitter &in, tfxEffectID &out);

	//Helper functions

	//Get a graph by GraphID
	Graph &GetGraph(EffectLibrary &library, GraphID &graph_id);

	//Set the udpate frequency for all particle effects - There may be options in the future for individual effects to be updated at their own specific frequency.
	void SetUpdateFrequency(float fps);

	inline float GetUpdateFrequency() { return UPDATE_FREQUENCY; }
	inline float GetUpdateTime() { return UPDATE_TIME; }
	inline float GetFrameLength() { return FRAME_LENGTH; }
	inline void SetLookUpFrequency(float frequency) {
		tfxLOOKUP_FREQUENCY = frequency;
	}
	inline void SetLookUpFrequencyOvertime(float frequency) {
		tfxLOOKUP_FREQUENCY_OVERTIME = frequency;
	}
	int GetShapesInPackage(const char *filename);
	int LoadEffectLibraryPackage(const char *filename, EffectLibrary &lib, void(*shape_loader)(const char *filename, ImageData &image_data, void *raw_image_data, int image_size, void *user_data) = nullptr, void *user_data = nullptr);

	//---
	//Effect Pool functions for managing and playing effects
	//Add an effect to an effect pool from a library by it's name. Returns true on success. the Effect ID in stored in the tfxEffectID that you pass to the function which
	//you can then use to access the effect in the pool later.
	//Overloaded functions are for retrieving effects from the library by its hash and by index
	bool PoolEffectFromLibrary(EffectLibrary &library, tfxEffectPool &storage, const char *name, tfxEffectID &out);
	bool PoolEffectFromLibrary(EffectLibrary &library, tfxEffectPool &storage, tfxKey path_hash, tfxEffectID &out);
	bool PoolEffectFromLibrary(EffectLibrary &library, tfxEffectPool &storage, unsigned int, tfxEffectID &out);
	bool PoolEffect(tfxEffectPool &storage, EffectEmitter &effect, tfxEffectID &out);
	//Prepare an effect template for setting up function call backs to customise the behaviour of the effect in realtime
	//Returns true on success.
	bool PrepareEffectTemplate(EffectLibrary &library, const char *name, tfxEffectTemplate &effect_template);
	//Add an effect to an effect pool using a template you created with PrepareEffectTemplate. Returns true on success. the Effect ID in stored in the tfxEffectID that you pass to the function which
	//you can then use to access the effect in the pool later.
	bool PoolEffectFromTemplate(tfxEffectPool &effect_pool, tfxEffectTemplate &effect_template, tfxEffectID &effect_id);
	//Retrieve an effect from an effect pool
	tfxEffect &GetEffect(tfxEffectPool &effect_pool, tfxEffectID &effect_id);
	//Check to see if an effect id is a ID that's currently in use in the effect pool
	bool ValidEffect(tfxEffectPool &effect_pool, tfxEffectID effect_id);
	//Update an effect in the effect pool. This is necessary to do on each frame update.
	void UpdateEffect(tfxEffectPool &effect_pool, tfxEffectID effect_id);
	//Initialise an effect pool to allocate all the necessary memory to hold the effects, emitters, particles and sprites
	void InitEffectPool(tfxEffectPool &effect_pool, unsigned int max_effects = 1000, unsigned int max_emitters = 10000, unsigned int max_particles = 100000);
	//Instantly stop an effect from spawning particles and remove all currently active particles
	void HardStopEffect(tfxEffectPool &effect_pool, tfxEffectID effect_id);
	//Make an effect stop spawning but let the current particles play out to the end of life.
	void SoftStopEffect(tfxEffectPool &effect_pool, tfxEffectID effect_id);
	//Restart an effect that has been stopped
	void StartEffect(tfxEffectPool &effect_pool, tfxEffectID effect_id);
	//Clear the effect pool of all effects. Any effect ids that you're currently using will not longer be valid. 
	//This function will completely start the effect pool from scratch, removing all cached memory ranges. Can be useful to defrag the memory and start it from fresh again.
	void ClearEffectPool(tfxEffectPool &effect_pool);
	//Free all the memory used in the effect pool. You would have to initialise it but trying to use it again. all Effect ids that used this effect pool
	//will no longer be valid after calling this function
	void FreeEffectPool(tfxEffectPool &effect_pool);
}

