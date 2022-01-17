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
#include <fstream>					//std::basic_ofstream
#include <sstream>					//std::basic_stringstream
#include <stdio.h>
#include <stdarg.h>					//va_list
#include <chrono>					//std::chrono::high_resolution_clock
#include <cctype>					//std::is_digit
#include <algorithm>
#include <stdint.h>
#include <assert.h>

namespace tfx {

#define TWO63 0x8000000000000000u 
#define TWO64f (TWO63*2.0)

	//----------------------------------------------------------
	//Forward declarations

	struct EffectEmitter;
	struct ParticleManager;
	struct EffectorStore;
	struct Particle;
	struct AnimationSettings;
	struct EffectLibrary;

	//--------------------------------------------------------------
	//macros
#define TFX_VERSION "Alpha 1"
#define TFX_VERSION_NUMBER 1

#define tfxMAX_FRAME 20000.f
#define tfxNullParent 0xFFFFFFFF
#define EmitterPropertiesCount 26

#define Del << "=" <<
#define Com << "," <<
#define EndLine << std::endl

typedef std::chrono::high_resolution_clock Clock;

//Override this for more layers, although currently the editor is fixed at 4
#ifndef tfxLAYERS
#define tfxLAYERS 4
#endif 

	//----------------------------------------------------------
	//enums/flags

	//Particle property that defines how a particle will rotate
	enum AngleSetting : unsigned char {
		tfxAlign,												//Align the particle with it's direction of travel
		tfxRandom,												//Chose a random angle at spawn time/flags
		tfxSpecify,												//Specify the angle at spawn time
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
#define	tfxBaseCount  7
#define	tfxVariationCount  8
#define	tfxOvertimeCount  15

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
		tfxGlobal_effect_angle,
		tfxGlobal_stretch,
		tfxGlobal_overal_scale,
		tfxGlobal_opacity,
		tfxGlobal_frame_rate,
		tfxGlobal_splatter,

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

		tfxVariation_life,
		tfxVariation_amount,
		tfxVariation_velocity,
		tfxVariation_width,
		tfxVariation_height,
		tfxVariation_weight,
		tfxVariation_spin,
		tfxVariation_motion_randomness,

		tfxOvertime_velocity,
		tfxOvertime_width,
		tfxOvertime_height,
		tfxOvertime_weight,
		tfxOvertime_spin,
		tfxOvertime_stretch,
		tfxOvertime_red,
		tfxOvertime_green,
		tfxOvertime_blue,
		tfxOvertime_opacity,
		tfxOvertime_frame_rate,
		tfxOvertime_motion_randomness,
		tfxOvertime_velocity_adjuster,
		tfxOvertime_intensity,
		tfxOvertime_direction,
		tfxGraphMaxIndex,
	};

	//EffectEmitter type - effect contains emitters, and emitters spawn particles, but they both share the same struct for simplicity
	enum EffectEmitterType : unsigned char {
		tfxEffect,
		tfxEmitter,
		tfxStage
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
		tfxEndAnimationSettings
	};

	typedef unsigned int tfxEmitterPropertyFlags;
	typedef unsigned int tfxVectorFieldFlags;
	typedef unsigned char tfxParticleFlags;
	typedef unsigned char tfxEmitterStateFlags;

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
		tfxEmitterPropertyFlags_use_vector_field = 1 << 19					//Enable the use of a vector field to apply forces to the particles
	};

	enum tfxParticleFlags_ : unsigned char {
		tfxParticleFlags_none = 0,
		tfxParticleFlags_fresh = 1 << 0,									//Particle has just spawned this frame	
		tfxParticleFlags_remove = 1 << 1,									//Particle will be removed this or next frame
	};

	enum tfxEmitterStateFlags_ : unsigned char {
		tfxEmitterStateFlags_none = 0,
		tfxEmitterStateFlags_stop_spawning = 1 << 0,							//Tells the emitter to stop spawning
		tfxEmitterStateFlags_remove = 1 << 1,									//Tells the effect/emitter to remove itself from the particle manager immediately
		tfxEmitterStateFlags_enabled = 1 << 2,									//the emitter is enabled. If flag is not set then it will not be added to the particle manager with AddEffect
		tfxEmitterStateFlags_retain_matrix = 1 << 3,							//Internal flag about matrix usage
		tfxEmitterStateFlags_no_tween_this_update = 1 << 4						//Internal flag generally, but you could use it if you want to teleport the effect to another location
	};

	enum tfxVectorFieldFlags_: unsigned char {
		tfxVectorFieldFlags_none = 0,
		tfxVectorFieldFlags_repeat_horizontal = 1 << 0,							//Field will repeat horizontally
		tfxVectorFieldFlags_repeat_vertical = 1 << 1								//Field will repeat vertically
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
	
	//these Variables determine the timing resolution that particles are updated at. So an Update frequency of 60 would mean that the particles are updated at 60 frames per secon d.
	static float UPDATE_FREQUENCY = 60.f;
	static float UPDATE_TIME = 1.f / UPDATE_FREQUENCY;
	static float FRAME_LENGTH = 1000.f / UPDATE_FREQUENCY;

	//Look up frequency determines the resolution of graphs that are compiled into look up arrays.
	static float tfxLOOKUP_FREQUENCY = 30.f;
	//Overtime frequency is for lookups that will vary in length depending on the lifetime of the particle. It should generally be a higher resolution than the base graphs
	static float tfxLOOKUP_FREQUENCY_OVERTIME = 10.f;

	//-----------------------------------------------------------
	//Utility things:

	//Credit to ocornut https://github.com/ocornut/imgui/commits?author=ocornut
	//std::vector replacement with some extra stuff and tweaks specific to Qulkan/TimelineFX
	template<typename T>
	struct tfxvec {
		unsigned int current_size;
		unsigned int capacity;
		T* data;

		inline tfxvec() { current_size = capacity = 0; data = NULL; }
		inline tfxvec(unsigned int qty) { current_size = capacity = 0; data = NULL; resize(qty); }
		inline tfxvec(T* from, T* to) { current_size = capacity = 0; data = NULL; auto current = from; while (current != to + 1) { push_back(*current); ++current; } }
		inline tfxvec(std::initializer_list<T> t) { current_size = capacity = 0; data = NULL; for (T element : t) { push_back(element); } }
		inline tfxvec(const tfxvec<T> &src) { current_size = capacity = 0; data = NULL; resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); }
		inline tfxvec<T>& operator=(const tfxvec<T>& src) { free_all(); resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); return *this; }
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
		inline tfxVec2 operator/(float v) { return tfxVec2(x / v, y / v); }
		inline float Squared() { return x * x + y * y; }
	};
	inline tfxVec2 operator*(float ls, tfxVec2 rs) { return tfxVec2(rs.x * ls, rs.y * ls); }

	struct tfxVec4 {
		float x, y, z, w;

		tfxVec4() { x = y = z = w = 0.f; }
		tfxVec4(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}

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
		inline unsigned int Length() const { return string.current_size ? string.current_size - 1 : 0; }
		void Appendf(const char *format, ...);
		inline void Append(char c) { 
			if (string.current_size) {
				string.pop();
			}
			string.push_back(c); 
			NullTerminate();
		}
		void NullTerminate() { string.push_back(NULL); }
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

	//------------------------------------------------------------

	//Structs
	//These are mainly internal structs
	//------------------------------------------------------------

	typedef tfxVec2 Point;

	struct AttributeNode {
		float frame;
		float value;

		bool is_curve;
		bool curves_initialised;

		Point left;
		Point right;

		unsigned int index;

		AttributeNode() : frame(0.f), value(0.f), is_curve(false), curves_initialised(false) { }
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
			is_curve = true;
		}

		/*
			Toggle whether this attribute node is curved or linear
		*/
		void ToggleCurve() {
			is_curve = !is_curve;
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
	};
	
	//This struct is used to store indexing data in order to index into large lists containing either the node data of graphs
	//or the lookup data of compiled graphs. This is so that we can upload that data into a buffer on the GPU to get the particles
	//updating in a compute shader.
	struct EffectLookUpData {
		GraphLookupIndex global_life;
		GraphLookupIndex global_amount;
		GraphLookupIndex global_velocity;
		GraphLookupIndex global_width;
		GraphLookupIndex global_height;
		GraphLookupIndex global_weight;
		GraphLookupIndex global_spin;
		GraphLookupIndex global_effect_angle;
		GraphLookupIndex global_stretch;
		GraphLookupIndex global_overal_scale;
		GraphLookupIndex global_opacity;
		GraphLookupIndex global_frame_rate;
		GraphLookupIndex global_splatter;

		GraphLookupIndex property_emission_angle;
		GraphLookupIndex property_emission_range;
		GraphLookupIndex property_emitter_angle;
		GraphLookupIndex property_splatter;
		GraphLookupIndex property_emitter_width;
		GraphLookupIndex property_emitter_height;
		GraphLookupIndex property_arc_size;
		GraphLookupIndex property_arc_offset;

		GraphLookupIndex base_life;
		GraphLookupIndex base_amount;
		GraphLookupIndex base_velocity;
		GraphLookupIndex base_width;
		GraphLookupIndex base_height;
		GraphLookupIndex base_weight;
		GraphLookupIndex base_spin;

		GraphLookupIndex variation_life;
		GraphLookupIndex variation_amount;
		GraphLookupIndex variation_velocity;
		GraphLookupIndex variation_width;
		GraphLookupIndex variation_height;
		GraphLookupIndex variation_weight;
		GraphLookupIndex variation_spin;
		GraphLookupIndex variation_motion_randomness;

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
		GraphLookupIndex overtime_frame_rate;
		GraphLookupIndex overtime_motion_randomness;
		GraphLookupIndex overtime_velocity_adjuster;
		GraphLookupIndex overtime_intensity;
		GraphLookupIndex overtime_direction;
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

		AttributeNode* AddNode(float frame, float value, bool is_curve = false, float x1 = 0, float y1 = 0, float x2 = 0, float y2 = 0);
		void AddNode(AttributeNode &node);
		void SetNode(uint32_t index, float frame, float value, bool is_curve = false, float x1 = 0, float y1 = 0, float x2 = 0, float y2 = 0);
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
	float LookupFastOvertime(Graph &graph, float age, float lifetime);
	float LookupFast(Graph &graph, float frame);
	float LookupPreciseOvertime(Graph &graph, float age, float lifetime);
	float LookupPrecise(Graph &graph, float frame);
	float GetRandomFast(Graph &graph, float frame);
	float GetRandomPrecise(Graph &graph, float frame);

	//Node Manipluation
	bool SetNode(Graph &graph, AttributeNode &node, float, float, bool, float = 0, float = 0, float = 0, float = 0);
	bool SetNode(Graph &graph, AttributeNode &node, float &frame, float &value);
	void SetCurve(Graph &graph, AttributeNode &node, bool is_left_curve, float &frame, float &value);
	bool MoveNode(Graph &graph, AttributeNode &node, float frame, float value, bool sort = true);
	bool SetNodeFrame(Graph &graph, AttributeNode &node, float &frame);
	bool SetNodeValue(Graph &graph, AttributeNode &node, float &value);
	void ClampNode(Graph &graph, AttributeNode &node);
	void ClampCurve(Graph &graph, Point &curve, AttributeNode &node);
	bool IsOvertimeGraph(GraphType type);
	bool IsGlobalGraph(GraphType type);

	struct GlobalAttributes {
		Graph life;
		Graph amount;
		Graph velocity;
		Graph width;
		Graph height;
		Graph weight;
		Graph spin;
		Graph effect_angle;
		Graph stretch;
		Graph overal_scale;
		Graph opacity;
		Graph frame_rate;
		Graph splatter;
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
	};

	struct VariationAttributes {
		Graph life;
		Graph amount;
		Graph velocity;
		Graph width;
		Graph height;
		Graph weight;
		Graph spin;
		Graph motion_randomness;
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
		Graph opacity;
		Graph frame_rate;
		Graph motion_randomness;
		Graph velocity_adjuster;
		Graph intensity;
		Graph direction;
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

	struct Base {
		tfxVec2 size;
		tfxVec2 random_size;
		float velocity;
		float height;
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
		unsigned int seed;
		bool seamless;
		bool loop;
		bool needs_recording;
		unsigned int needs_exporting;
		float max_radius;
		unsigned int largest_frame;
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

		//use this definition if you need more spefic data to point to the image texture in whatever renderer you're using
		//Just define tfxCUSTOM_IMAGE_DATA before you include timelinefx.h
#ifdef tfxCUSTOM_IMAGE_DATA
		tfxCUSTOM_IMAGE_DATA
#endif // tfxCUSTOM_IMAGE_DATA

		//Each particle shape saved in an effect library has a unique index
		uint32_t shape_index;
		//The size of one frame of the image
		tfxVec2 image_size;
		//The number of frames in the image, can be one or more
		float animation_frames;
		//Maximum distance to the nearest transparent edge of the image
		float max_radius;
		int import_filter;

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
		BlendMode blend_mode;
		//Currently there are 4 types of emission, point, line, area and ellipse
		EmissionType emission_type;
		//Should particles emit towards the center of the emitter or away, or in a specific direction
		EmissionDirection emission_direction;
		//How particles should behave when they reach the end of the line
		LineTraversalEndBehaviour end_behaviour;
		//The rotation of particles when they spawn, or behave overtime if tfxAlign is used
		AngleSetting angle_setting = AngleSetting::tfxRandom;

		//Bit field of various boolean flags
		tfxEmitterPropertyFlags flags;

		//Offset to draw particles at
		tfxVec2 image_handle;
		//Offset of emitters
		tfxVec2 emitter_handle;
		//When single or one shot flags are set, spawn this amount of particles in one go
		unsigned int spawn_amount;
		//Layer of the particle manager that the particle is added to
		unsigned int layer;
		//The shape being used for all particles spawned from the emitter
		unsigned int shape_index;

		//Angle added to the rotation of the particle when spawned
		float angle_offset;
		//The number of rows/columns/ellipse/line points in the grid when spawn on grid flag is used
		tfxVec2 grid_points;
		//The number of millisecs before an effect or emitter will loop back round to the beginning of it's graph lookups
		float loop_length;
		//The start frame index of the animation
		float start_frame;
		//The final frame index of the animation
		float end_frame;

		EmitterProperties() :
			angle_offset(0),
			flags(tfxEmitterPropertyFlags_image_handle_auto_center | tfxEmitterPropertyFlags_grid_spawn_clockwise | tfxEmitterPropertyFlags_emitter_handle_auto_center | tfxEmitterPropertyFlags_global_uniform_size | tfxEmitterPropertyFlags_base_uniform_size | tfxEmitterPropertyFlags_lifetime_uniform_size),
			image(nullptr),
			image_handle(tfxVec2()),
			spawn_amount(1),
			blend_mode(BlendMode::tfxAlpha),
			emission_type(EmissionType::tfxPoint),
			emission_direction(EmissionDirection::tfxOutwards),
			grid_points({ 10.f, 10.f }),
			emitter_handle(tfxVec2()),
			end_behaviour(LineTraversalEndBehaviour::tfxLoop),
			loop_length(0.f),
			layer(0),
			shape_index(1),
			start_frame(0),
			end_frame(0)
		{ }
	};

	//Store the current state of the Effect/Emitter. All these values can change over the lifetime of the effect/emitter.
	struct EffectEmitterState {
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
		float motion_randomness;
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

		EffectEmitterState() : grid_coords(tfxVec2()), single_shot_done(false), age(0.f), amount_remainder(0.f) {}
	};

	//An EffectEmitter can either be an effect which stores emitters and global graphs for affecting all the attributes in the emitters
	//Or it can be an emitter which spawns all of the particles. Effectors are stored in the particle manager effects list buffer.
	struct EffectEmitter {

		//Position, scale and rotation values
		FormState local;
		FormState world;
		FormState captured;
		//2d matrix for transformations
		Matrix2 matrix;

		//The current state of the effect/emitter
		EffectEmitterState current;
		//All of the properties of the effect/emitter
		EmitterProperties properties;

		//Name of the effect
		tfxText name;						//Todo: Do we need this here?
		//A hash of the directory path to the effect ie Flare/spark
		tfxKey path_hash;
		//Is this a tfxEffect or tfxEmitter
		EffectEmitterType type;
		//A pointer to the library this effect belongs
		EffectLibrary *library;
		//A pointer the the particle manager that this has been added
		ParticleManager *pm;
		//The index within the library that this exists at
		unsigned int library_index;
		//The number of sub_effects still in use
		unsigned int active_children;
		//The number of particles active within this emitter
		unsigned int particle_count;
		//The number of frames before this is removed from the particle manager after particle count is 0
		unsigned int timeout;
		//Internal, keep track of idle frames
		unsigned int timeout_counter;
		//Every effect and emitter in the library gets a unique id
		unsigned int uid;
		//The max_radius of the emitter, taking into account all the particles that have spawned and active
		float max_radius;
		//List of sub_effects ( effects contain emitters, emitters contain sub effects )
		tfxvec<EffectEmitter> sub_effectors;

		//Custom user data, can be accessed in callback functions
		void *user_data;

		//All graphs that the effect uses to lookup attribute values are stored in the library. These variables here are indexes to the array where they're stored
		unsigned int global;
		unsigned int property;
		unsigned int base;
		unsigned int variation;
		unsigned int overtime;
		//Experitment: index into the lookup index data in the effect library
		unsigned int lookup_node_index;
		unsigned int lookup_value_index;
		//Index to animation settings stored in the effect library. Would like to move this at some point
		unsigned int animation_settings;
		//The maximum amount of life that a particle can be spawned with taking into account base + variation life values
		float max_life;
		//Pointer to the root effect
		EffectEmitter *root_parent;
		//Pointer to the immediate parent
		EffectEmitter *parent;
		//Pointer to the next pointer in the particle manager buffer. 
		EffectEmitter *next_ptr;
		//Pointer to the sub effect's particle that spawned it
		Particle *parent_particle;

		//Custom fuction pointers that you can use to override attributes and affect the effect/emitters behaviour in realtime
		//See EffectEmitterTemplate for applying this callbacks
		void(*update_callback)(EffectEmitter &effectemitter);		//Called after the state has been udpated
		void(*particle_update_callback)(Particle &particle);		//Called for each particle that has been udpated, but before it's state is updated (so you can override behaviour first)
		void(*particle_onspawn_callback)(Particle &particle);		//Called as each particle is spawned.

		tfxEmitterStateFlags flags;

		EffectEmitter() : particle_count(0),
			root_parent(nullptr),
			parent(nullptr),
			parent_particle(nullptr),
			user_data(nullptr),
			flags(tfxEmitterStateFlags_no_tween_this_update | tfxEmitterStateFlags_enabled),
			timeout(100),
			timeout_counter(0),
			animation_settings(0),
			library(nullptr),
			update_callback(nullptr),
			particle_update_callback(nullptr),
			particle_onspawn_callback(nullptr)
		{
		}
		~EffectEmitter();
		bool operator < (const EffectEmitter& e) const
		{
			return (properties.layer < e.properties.layer);
		}

		//API functions
		//Tell the effect to stop spawning so that eventually particles will expire and the effect will be removed from the particle manager
		inline void SoftExpire();

		void Rotate(float);
		void SetAngle(float);
		void Scale(const tfxVec2&);
		void Scale(float x, float y);
		void Move(const tfxVec2&);
		void Move(float, float);
		void Position(const tfxVec2&);
		void Position(float, float);
		void SetUpdateCallback(void(*callback)(EffectEmitter &effectemitter));
		void SetUserData(void *data);
		void *GetUserData();
		void SetTimeout(unsigned int frames);

		//Override graph functions for use in update_callback
		//Some of these change the same state and property values, but they're named differently just to make it clearer as to whether you're overriding kEffect or a kEmitter.

		//Global Effect Overrides
		inline void OverrideGlobalLife(float value) { current.life; }
		inline void OverrideGlobalAmount(float value) { current.amount = value; }
		inline void OverrideGlobalVelocity(float value) { current.velocity = value; }
		inline void OverrideGlobalWidth(float value) { current.size.x = value; }
		inline void OverrideGlobalHeight(float value) { current.size.y = value; }
		inline void OverrideGlobalSpin(float value) { current.spin = value; }
		inline void OverrideGlobalEffectAngle(float value) { local.rotation = value; }
		inline void OverrideGlobalStretch(float value) { current.stretch = value; }
		inline void OverrideGlobalOveralScale(float value) { current.overal_scale = value; }
		inline void OverrideGlobalOpacity(float value) { current.color.a = value; }
		inline void OverrideGlobalSplatter(float value) { current.splatter = value; }
		//Todo: implement this
		//inline void OverrideGlobalFrameRate(float value) { properties. = value; }

		//Property Emitter Overrides
		inline void OverridePropertyEmissionAngle(float value) { current.emission_angle = value; }
		inline void OverridePropertyEmissionRange(float value) { current.emission_angle_variation = value; }
		inline void OverridePropertyEmitterAngle(float value) { local.rotation = value; }
		inline void OverridePropertySplatter(float value) { current.splatter = value; }
		inline void OverridePropertyEmitterWidth(float value) { current.emitter_size.x = value; }
		inline void OverridePropertyEmitterHeight(float value) { current.emitter_size.y = value; }
		inline void OverridePropertyArcSize(float value) { current.arc_size = value; }
		inline void OverridePropertyArcOffset(float value) { current.arc_offset = value; }

		//Base Emitter Overrides
		inline void OverrideBaseLife(float value) { current.life = value; }
		inline void OverrideBaseAmount(float value) { current.amount = value; }
		inline void OverrideBaseVelocity(float value) { current.velocity = value; }
		inline void OverrideBaseWidth(float value) { current.size.x = value; }
		inline void OverrideBaseHeight(float value) { current.size.y = value; }
		inline void OverrideBaseWeight(float value) { current.weight = value; }
		inline void OverrideBaseSpin(float value) { current.spin = value; }

		//Variation Emitter Overrides
		inline void OverrideVariationLife(float value) { current.life_variation = value; }
		inline void OverrideVariationAmount(float value) { current.amount_variation = value; }
		inline void OverrideVariationVelocity(float value) { current.velocity_variation = value; }
		inline void OverrideVariationWidth(float value) { current.size_variation.x = value; }
		inline void OverrideVariationHeight(float value) { current.size_variation.y = value; }
		inline void OverrideVariationWeight(float value) { current.weight_variation = value; }
		inline void OverrideVariationSpin(float value) { current.spin_variation = value; }
		inline void OverrideVariationMotionRandomness(float value) { current.motion_randomness = value; }

		//Property Flags
		inline void OverridePropertyFlag(tfxEmitterPropertyFlags flag, bool value) { if (value) properties.flags |= flag; else properties.flags &= ~flag; }

		//Internal functions
		EffectEmitter& AddEmitter(EffectEmitter &e);
		EffectEmitter& AddEffect(EffectEmitter &e);
		EffectEmitter& AddEffect();
		EffectEmitter& AddEffector(EffectEmitterType type = tfxEmitter);
		EffectEmitter* GetRootEffect();
		void ReIndex();
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
		Graph* GetGraphByType(GraphType type);
		unsigned int GetGraphIndexByType(GraphType type);
		void CompileGraphs();
		void InitialiseUninitialisedGraphs();
		void SetName(const char *n);

		void TransformEffector(EffectEmitter &parent, bool relative_position = true, bool relative_angle = false);
		void TransformEffector(Particle &parent, bool relative_position = true, bool relative_angle = false);
		void Update();
		void SpawnParticles();
		void UpdateEmitterState();
		void UpdateEffectState();
		float GetEmissionDirection(Particle& p);
		void ReSeed(uint64_t seed = 0);
		bool HasSingle();
		bool RenameSubEffector(EffectEmitter &effect, const char *new_name);
		bool NameExists(EffectEmitter &effect, const char *name);
		void FreeGraphs();
		void NoTweenNextUpdate();

		void ClearColors();
		void AddColorOvertime(float frame, tfxRGB color);
		void Clone(EffectEmitter &clone, EffectEmitter *root_parent, EffectLibrary *destination_library);
		void EnableAllEmitters();
		void EnableEmitter();
		void DisableAllEmitters();
		void DisableAllEmittersExcept(EffectEmitter &emitter);
	};

	struct EffectEmitterTemplate {
		tfxStorageMap<EffectEmitter*> paths;
		EffectEmitter effect_template;

		void AddPath(EffectEmitter &effectemitter, tfxText path) {
			paths.Insert(path, &effectemitter);
			for (auto &sub : effectemitter.sub_effectors) {
				tfxText sub_path = path;
				sub_path.Appendf("/%s", sub.name);
				AddPath(sub, sub_path);
			}
		}

		inline EffectEmitter &Effect() { return effect_template; }
		inline EffectEmitter *Get(tfxText path) { if (paths.ValidName(path)) return paths.At(path); return nullptr; }
		inline void SetUserData(tfxText path, void *data) { if(paths.ValidName(path)) paths.At(path)->user_data = data; }
		inline void SetUserData(void *data) { effect_template.user_data = data; }
		void SetUserDataAll(void *data);
		inline void SetUpdateCallback(tfxText path, void(*update_callback)(EffectEmitter &effectemitter)) { if (paths.ValidName(path)) paths.At(path)->update_callback = update_callback; }
		inline void SetUpdateCallback(void(*update_callback)(EffectEmitter &effectemitter)) { effect_template.update_callback = update_callback; }
		void SetUpdateCallbackAll(void(*update_callback)(EffectEmitter &effectemitter));
		void SetParticleUpdateCallback(tfxText path, void(*particle_update_callback)(Particle &particle));
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

	//Initial particle struct, looking to optimise this and make as small as possible
	//These are spawned by effector emitter types
	//Particles are stored in the particle manager particle buffer.
	//I really think that tweened frames should be ditched in favour of delta time so captured can be ditched
	//180 bytes
	struct Particle {
		FormState local;				//The local position of the particle, relative to the emitter.
		FormState world;				//The world position of the particle relative to the screen.
		FormState captured;				//The captured world coords for tweening
		Matrix2 matrix;					//Simple 2d matrix for transforms
		Base base;						//Base values created when the particle is spawned. They can be different per particle due to variations
		//tfxVec2 velocity_normal;		//stores the current direction of travel for the particle
		//tfxVec2 velocity;				//=velocity_normal * base.velocity * velocity_scale (velocity overtime) Gets added to the particle coords each frame
		//tfxVec2 handle;				//The image handle
		float age;						//The age of the particle, used by the controller to look up the current state on the graphs
		float max_age;					//max age before the particle expires
		float image_frame_rate;			//current frame rate of the image if it's an animation
		float velocity_scale;			//Current velocity overtime
		//float direction;				//Current direction of travel in radians
		float emission_angle;			//Emission angle of the particle at spawn time
		//float spin;					//Current spin overtime
		float image_frame;				//Current frame of the image if it's an animation
		float distance_travelled;		//Used in edge traversal and kLoop to make the particle start back at the beginning of the line again
		float weight_acceleration;		//The current amount of gravity applied to the y axis of the particle each frame
		float motion_randomness;		//The random velocity added each frame
		float motion_randomness_speed;
		//float motion_randomness_direction;
		//float motion_tracker;
		float intensity;				//Color is multiplied by this value in the shader to increase the brightness of the particles
		tfxRGBA8 color;					//Colour of the particle
		tfxParticleFlags flags;			//flags for different states
		EffectEmitter *parent;			//pointer to the emitter that emitted the particle.

		//Internal use variables
		Particle *next_ptr;

		//Override functions, you can use these inside an update_callback if you want to modify the particle's behaviour
		inline void OverridePosition(float x, float y) { local.position.x = x; local.position.y = y; }
		inline void OverrideSize(float x, float y) { local.scale.x = x; local.scale.y = y; }
		//inline void OverrideVelocity(float x, float y) { velocity.x = x; velocity.y = y; }
		inline void OverrideVelocityScale(float v) { velocity_scale = v; }
		inline void OverrideRotation(float r) { local.rotation = r; }
		inline void OverrideRotationDegrees(float d) { local.rotation = tfxDegrees(d); }
		//inline void OverrideDirection(float r) { direction = r; }
		//inline void OverrideDirectionDegrees(float d) { direction = tfxDegrees(d); }
		//inline void OverrideSpin(float r) { direction = r; }
		//inline void OverrideSpinDegrees(float d) { spin = tfxDegrees(d); }
		inline void OverrideAge(float a) { age = a; }
		inline void OverrideImageFrameRate(float fr) { image_frame_rate = fr; }
		inline void OverrideImageFrame(float f) { image_frame = f; }
		inline void OverrideWeight(float w) { weight_acceleration = w; }
		inline void OverrideIntensity(float i) { intensity = i; }
		inline void OverrideBaseSize(float width, float height) { base.size.x = width; base.size.y = height; }
		inline void OverrideBaseVelocity(float v) { base.velocity = v; }
		inline void OverrideBaseWeight(float w) { base.weight = w; }
		inline void OverrideBaseSpin(float s) { base.spin = s; }
		inline void OverrideBaseSpinDegrees(float d) { base.spin = tfxDegrees(d); }

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
		BlendMode blend_mode;
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
		tfxvec<EffectLookUpData> compiled_lookup_indexes;
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
		EffectEmitter *GetEffect(tfxText path);
		void PrepareEffectTemplate(tfxText path, EffectEmitterTemplate &effect);

		//Mainly internal functions
		EffectEmitter &AddEffect(EffectEmitter &effect);
		void UpdateEffectPaths();
		void UpdateEffectPaths(EffectEmitter &effect);
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
		void UpdateAllNodes();
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

	//Use the particle manager to add effects to your scene
	struct ParticleManager {

		//Particles are stored using double buffering, every update the particle is moved to the next list, so particle deletion happens by not adding to the next list. Given that particles can have varying
		//lifetimes and can expire at anytime, this seemed like the easiest way to do this, although I dare say that there's faster ways of doing it, multi-threading it might be a good start there.
		tfxvec<Particle> particles[tfxLAYERS][2];
		//Effects are also stored using double buffering. Effects stored here are "fire and forget", so you won't be able to apply changes to the effect in realtime. If you want to do that then 
		//store the effect yourself and manually update it. You just need to assign a particle manager to the effect so that it knows which particle manager to use to update particles and any sub effects
		//You can still apply changes to effect before you actually add it to the particle manager though
		tfxvec<EffectEmitter> effects[2];
		//The maximum number of effects that can be updated per frame in the particle manager. If you're running effects with particles that have sub effects then this number might need 
		//to be relatively high depending on your needs. Use Init to udpate the sizes if you need to. Best to call Init at the start with the max numbers that you'll need for your application and don't adjust after.
		unsigned int max_effects;
		//The maximum number of particles that can be updated per frame per layer. #define tfxLAYERS to set the number of allowed layers. This is currently 4 by default
		unsigned int max_particles_per_layer;
		//The current particle buffer in use, can be either 0 or 1
		unsigned int current_pbuff;
		//The current effect buffer in use, can be either 0 or 1
		unsigned int current_ebuff;
		//The current number of particles in each buffer
		unsigned int particle_count;
		//Callback to the render function that renders all the particles - is this not actually useful now, just use GetParticleBuffer() in your own render function
		void(*render_func)(float, void*, void*);
		bool disable_spawing;
		bool force_capture;
		//Used if the emitter changed attributes and spawned particles need updating their base values. Primarily used by the editor
		bool update_base_values;	
		LookupMode lookup_mode;

		//These can possibly be removed at some point, they're debugging variables
		unsigned int particle_id;

		ParticleManager() : force_capture(false), disable_spawing(false), lookup_mode(tfxFast), max_effects(10000), max_particles_per_layer(50000), update_base_values(false) { }
		~ParticleManager();
		EffectEmitter &operator[] (unsigned int index);

		//Initialise the particle manager with the maximum number of particles and effects that you want the manager to update per frame
		void Init(unsigned int effects_limit = 10000, unsigned int particle_limit_per_layer = 50000);
		//Update the particle manager. Call this once per frame in your logic udpate.
		void Update();
		//Get the current particle buffer that contains all particles currently active. The particle manager can have layers in order to control draw order of particles.
		//Pass the layer that you want to get.
		tfxvec<Particle> *GetParticleBuffer(unsigned int layer);
		//Add an effect to the particle manager. Pass an EffectEmitter pointer if you want to change the effect on the fly. Once you add the effect to the particle manager
		//then it's location in the buffer will keep changing as effects are updated and added and removed. The tracker will be updated accordingly each frame so you will always
		//have access to the effect if you need it.
		void AddEffect(EffectEmitter &effect, unsigned int buffer);
		void AddEffect(EffectEmitterTemplate &effect, unsigned int buffer);
		//Clear all effects and particles in the particle manager
		void ClearAll();
		//Soft expire all the effects so that the particles complete their animation first
		void SoftExpireAll();

		//This can be removed soon, don't think there's much need for it, use GetParticleBuffer and render that instead in your own function
		void Render(float tween, void *data);
		void SetRenderCallback(void func(float, void*, void*));

		//Internal use only
		Particle &GrabParticle(unsigned int layer);
		unsigned int AddParticle(unsigned int layer, Particle &p);
		//float Record(unsigned int frames, unsigned int start_frame, std::array<tfxvec<ParticleFrame>, 1000> &particle_frames);
		unsigned int ParticleCount();
		void ClearDepths();
		inline Particle* SetNextParticle(unsigned int layer, Particle &p, unsigned int buffer);
		inline EffectEmitter* SetNextEffect(EffectEmitter &e, unsigned int buffer);
		void UpdateBaseValues();
		tfxvec<EffectEmitter> *GetEffectBuffer();
		void SetLookUpMode(LookupMode mode);

		inline bool FreeCapacity(unsigned int layer) { return particles[layer][current_pbuff].current_size <= max_particles_per_layer; }
		inline bool FreeEffectCapacity() {
			return effects[0].current_size + effects[1].current_size < max_effects;
		}
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
	void SaveDataFile(tfxStorageMap<DataEntry> &config, const char* path = "");
	void LoadDataFile(tfxStorageMap<DataEntry> &config, const char* path);
	void StreamProperties(EmitterProperties &property, std::stringstream &file);
	void StreamGraph(const char * name, Graph &graph, std::stringstream &file);
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
	void TransformParticlePrevious(Particle &p, EffectEmitter &e);
	bool ControlParticle(Particle &p, EffectEmitter &e);
	FormState Tween(float tween, FormState &world, FormState &captured);
	tfxVec2 InterpolateVec2(float, const tfxVec2&, const tfxVec2&);
	float Interpolatef(float tween, float, float);
	int ValidateEffectLibrary(const char *filename);
	void ReloadBaseValues(Particle &p, EffectEmitter &e);

	//Helper functions

	//Get a graph by GraphID
	Graph &GetGraph(EffectLibrary &library, GraphID &graph_id);
	//Get a node by GraphID
	Graph &GetGraphNode(EffectLibrary &library, GraphID &graph_id);

	//Set the udpate frequency for all particle effects - There may be options in the future for individual effects to be updated at their own specific frequency.
	inline void SetUpdateFrequency(float fps) {
		UPDATE_FREQUENCY = fps;
		UPDATE_TIME = 1.f / UPDATE_FREQUENCY;
		FRAME_LENGTH = 1000.f / UPDATE_FREQUENCY;
	}
	inline void SetLookUpFrequency(float frequency) {
		tfxLOOKUP_FREQUENCY = frequency;
	}
	inline void SetLookUpFrequencyOvertime(float frequency) {
		tfxLOOKUP_FREQUENCY_OVERTIME = frequency;
	}
	int GetShapesInLibrary(const char *filename);
	int LoadEffectLibrary(const char *filename, EffectLibrary &lib, void(*shape_loader)(const char *filename, ImageData &image_data, void *raw_image_data, int image_size, void *user_data) = nullptr, void *user_data = nullptr);

	//Particle manager functions
	void StopSpawning(ParticleManager &pm);
	void RemoveAllEffects(ParticleManager &pm);
	void InitParticleManager(ParticleManager &pm, unsigned int effects_limit, unsigned int particle_limit_per_layer);
	void AddEffect(ParticleManager &pm, EffectEmitter &effect, float x = 0.f, float y = 0.f);
	void AddEffect(ParticleManager &pm, EffectEmitterTemplate &effect, float x = 0.f, float y = 0.f);

}

