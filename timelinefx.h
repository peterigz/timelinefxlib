#define tfxENABLE_PROFILING
#define tfxPROFILER_SAMPLES 60
//#define tfxTRACK_MEMORY
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
		//In this example, a compute shader is used to transform_3d all the vertices into the right place by sending a batch of quads. A quad just has the size, orientation, color and UV coords, the compute
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
#include <Windows.h>
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
#include <thread>					//only using this for std::thread::hardware_ concurrency()

namespace tfx {

#define TWO63 0x8000000000000000u 
#define TWO64f (TWO63*2.0)
#define tfxPI 3.14159265359f
#define tfx360Radians 6.28319f
#define tfx180Radians 3.14159f
#define tfx90Radians 1.5708f
#define tfxMAXDEPTH 3

	//----------------------------------------------------------
	//Forward declarations

	struct tfxEffectEmitter;
	struct tfxParticleManager;
	struct tfxEffectTemplate;
	struct tfxParticleData;
	struct tfxComputeSprite;
	struct tfxComputeParticle;
	struct tfxSpriteSheetSettings;
	struct tfxSpriteDataSettings;
	struct tfxLibrary;
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
#define tfxINVALID_SPRITE 0x0FFFFFFF
#define tfxEmitterPropertiesCount 26

#define tfxDel << "=" <<
#define tfxCom << "," <<
#define tfxEndLine << std::endl

#define tfxDelt "=" 
#define tfxComt ","
#define tfxEndLinet "\n"

#define tfxMin(a, b) (((a) < (b)) ? (a) : (b))
#define tfxMax(a, b) (((a) > (b)) ? (a) : (b))

	typedef std::chrono::high_resolution_clock tfxClock;

#define tfxAPI				//Function marker for any functions meant for external/api use

	//Override this for more layers, although currently the editor is fixed at 4
#ifndef tfxLAYERS
#define tfxLAYERS 4

/*
Helper macro to place inside a for loop, for example:
for(tfxEachLayer)
You can then use layer inside the loop to get the current layer
*/
#define tfxEachLayer int layer = 0; layer != tfxLAYERS; ++layer

//Internal use macro
#define tfxEachLayerDB int layer = 0; layer != tfxLAYERS * 2; ++layer
#endif 
//type defs
	typedef unsigned int tfxU32;
	typedef unsigned int tfxEmitterID;
	typedef int tfxS32;
	typedef unsigned long long tfxU64;
	typedef long long tfxS64;
	typedef tfxU32 tfxEffectID;
	typedef unsigned long long tfxKey;
	typedef tfxU32 tfxParticleID;
	typedef short tfxShort;
	typedef unsigned short tfxUShort;

	inline tfxParticleID MakeParticleID(tfxU32 bank_index, tfxU32 particle_index) {
		return ((bank_index & 0x00000FFF) << 20) + particle_index;
	}

	inline tfxU32 ParticleIndex(tfxParticleID id) {
		return id & 0x000FFFFF;
	}

	inline tfxU32 ParticleBank(tfxParticleID id) {
		return (id & 0xFFF00000) >> 20;
	}

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

	//#define tfxUSEAVX

	//Define tfxUSEAVX if you want to compile and use AVX simd operations for updating particles, otherwise SSE will be
	//used by default
#ifdef tfxUSEAVX
#define tfxDataWidth 8	
	typedef __m256 tfxWideFloat;
	typedef __m256i tfxWideInt;
#define tfxWideLoad _mm256_load_ps
#define tfxWideLoadi _mm256_load_si256
#define tfxWideSet _mm256_set_ps
#define tfxWideSetSingle _mm256_set1_ps
#define tfxWideSeti _mm256_set_epi32
#define tfxWideSetSinglei _mm256_set1_epi32
#define tfxWideAdd _mm256_add_ps
#define tfxWideSub _mm256_sub_ps
#define tfxWideMul _mm256_mul_ps
#define tfxWideDiv _mm256_div_ps
#define tfxWideAddi _mm256_add_epi32
#define tfxWideSubi _mm256_sub_epi32
#define tfxWideMuli _mm256_mul_epi32
#define tfxWideSqrt _mm256_sqrt_ps
#define tfxWideShiftRight _mm256_srli_epi32
#define tfxWideShiftLeft _mm256_slli_epi32
#define tfxWideGreaterEqual(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_GE_OS)
#define tfxWideGreater(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_GT_OS)
#define tfxWideGreateri(v1, v2) _mm256_cmpgt_epi32(v1, v2)
#define tfxWideLess(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_LT_OS)
#define tfxWideLessEqeual(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_LE_OS)
#define tfxWideEqualsi _mm256_cmpeq_epi32 
#define tfxWideStore _mm256_store_ps
#define tfxWideStorei _mm256_store_si256
#define tfxWideCasti _mm256_castps_si256
#define tfxWideCast _mm256_castsi256_ps 
#define tfxWideConverti _mm256_cvttps_epi32 
#define tfxWideConvert	_mm256_cvtepi32_ps 
#define tfxWideMin _mm256_min_ps
#define tfxWideMax _mm256_max_ps
#define tfxWideMini _mm256_min_epi32
#define tfxWideMaxi _mm256_max_epi32
#define tfxWideOr _mm256_or_ps
#define tfxWideOri _mm256_or_si256
#define tfxWideXOri _mm256_xor_si256
#define tfxWideXOr _mm256_xor_ps
#define tfxWideAnd _mm256_and_ps
#define tfxWideAndi _mm256_and_si256
#define tfxWideAndNot _mm256_andnot_ps
#define tfxWideAndNoti _mm256_andnot_si256
#define tfxWideSetZero _mm256_setzero_si256
#define tfxWideEqualsi _mm256_cmpeq_epi32 
#define tfxWideAndNot _mm256_andnot_ps
#define tfxWideEquals _mm256_cmpeq_ps
#define tfxWideLookupSet(lookup, index) tfxWideSet(lookup[index[7]], lookup[index[6]], lookup[index[5]], lookup[index[4]], lookup[index[3]], lookup[index[2]], lookup[index[1]], lookup[index[0]] )
#define tfxWideLookupSeti(lookup, index) tfxWideSeti(lookup[index[7]], lookup[index[6]], lookup[index[5]], lookup[index[4]], lookup[index[3]], lookup[index[2]], lookup[index[1]], lookup[index[0]] )

	const __m256 tfxWIDEF3_4 = _mm256_set1_ps(1.0f / 3.0f);
	const __m256 tfxWIDEG3_4 = _mm256_set1_ps(1.0f / 6.0f);
	const __m256 tfxWIDEG32_4 = _mm256_set1_ps((1.0f / 6.0f) * 2.f);
	const __m256 tfxWIDEG33_4 = _mm256_set1_ps((1.0f / 6.0f) * 3.f);
	const __m256i tfxWIDEONEi = _mm256_set1_epi32(1);
	const __m256 tfxWIDEONE = _mm256_set1_ps(1.f);
	const __m256 tfxWIDE255 = _mm256_set1_ps(255.f);
	const __m256 tfxWIDEZERO = _mm256_set1_ps(0.f);
	const __m256 tfxWIDETHIRTYTWO = _mm256_set1_ps(32.f);
	const __m256i tfxWIDEFF = _mm256_set1_epi32(0xFF);
	const __m256 tfxPWIDESIX = _mm256_set1_ps(0.6f);

	typedef union {
		__m256i m;
		int a[8];
	} tfxWideArrayi;

	typedef union {
		__m256 m;
		float a[8];
	} tfxWideArray;

#else
#define tfxDataWidth 4	
	typedef __m128 tfxWideFloat;
	typedef __m128i tfxWideInt;
#define tfxWideLoad _mm_load_ps
#define tfxWideLoadi _mm_load_si128
#define tfxWideSet _mm_set_ps
#define tfxWideSetSingle _mm_set_ps1
#define tfxWideSeti _mm_set_epi32
#define tfxWideSetSinglei _mm_set1_epi32
#define tfxWideAdd _mm_add_ps
#define tfxWideSub _mm_sub_ps
#define tfxWideMul _mm_mul_ps
#define tfxWideDiv _mm_div_ps
#define tfxWideAddi _mm_add_epi32
#define tfxWideSubi _mm_sub_epi32
#define tfxWideMuli _mm_mul_epu32
#define tfxWideSqrt _mm_sqrt_ps
#define tfxWideShiftRight _mm_srli_epi32
#define tfxWideShiftLeft _mm_slli_epi32
#define tfxWideGreaterEqual(v1, v2) _mm_cmpge_ps(v1, v2)
#define tfxWideGreater(v1, v2) _mm_cmpgt_ps(v1, v2)
#define tfxWideGreateri(v1, v2) _mm_cmpgt_epi32(v1, v2)
#define tfxWideLessEqual(v1, v2) _mm_cmple_ps(v1, v2)
#define tfxWideLess(v1, v2) _mm_cmplt_ps(v1, v2)
#define tfxWideLessi(v1, v2) _mm_cmplt_epi32(v1, v2)
#define tfxWideStore _mm_store_ps
#define tfxWideStorei _mm_store_si128
#define tfxWideCasti _mm_castps_si128 
#define tfxWideCast _mm_castsi128_ps
#define tfxWideConverti _mm_cvttps_epi32 
#define tfxWideConvert _mm_cvtepi32_ps 
#define tfxWideMin _mm_min_ps
#define tfxWideMax _mm_max_ps
#define tfxWideMini _mm_min_epi32
#define tfxWideMaxi _mm_max_epi32
#define tfxWideOr _mm_or_ps
#define tfxWideOri _mm_or_si128
#define tfxWideXOr _mm_xor_ps
#define tfxWideXOri _mm_xor_si128
#define tfxWideAnd _mm_and_ps
#define tfxWideAndi _mm_and_si128
#define tfxWideAndNot _mm_andnot_ps
#define tfxWideAndNoti _mm_andnot_si128
#define tfxWideSetZeroi _mm_setzero_si128
#define tfxWideSetZero _mm_setzero_ps
#define tfxWideEqualsi _mm_cmpeq_epi32 
#define tfxWideEquals _mm_cmpeq_ps
#define tfxWideShufflei _mm_shuffle_epi32

#define tfxWideLookupSet(lookup, index) tfxWideSet( lookup[index.a[3]], lookup[index.a[2]], lookup[index.a[1]], lookup[index.a[0]] )
#define tfxWideLookupSetMember(lookup, member, index) tfxWideSet( lookup[index.a[3]].member, lookup[index.a[2]].member, lookup[index.a[1]].member, lookup[index.a[0]].member )
#define tfxWideLookupSetMemberi(lookup, member, index) tfxWideSeti( lookup[index.a[3]].member, lookup[index.a[2]].member, lookup[index.a[1]].member, lookup[index.a[0]].member )
#define tfxWideLookupSet2(lookup1, lookup2, index1, index2) tfxWideSet( lookup1[index1.a[3]].lookup2[index2.a[3]], lookup1[index1.a[2]].lookup2[index2.a[2]], lookup1[index1.a[1]].lookup2[index2.a[1]], lookup1[index1.a[0]].lookup2[index2.a[0]] )
#define tfxWideLookupSeti(lookup, index) tfxWideSeti( lookup[index.a[3]], lookup[index.a[2]], lookup[index.a[1]], lookup[index.a[0]] )

	const __m128 tfxWIDEF3_4 = _mm_set_ps1(1.0f / 3.0f);
	const __m128 tfxWIDEG3_4 = _mm_set_ps1(1.0f / 6.0f);
	const __m128 tfxWIDEG32_4 = _mm_set_ps1((1.0f / 6.0f) * 2.f);
	const __m128 tfxWIDEG33_4 = _mm_set_ps1((1.0f / 6.0f) * 3.f);
	const __m128i tfxWIDEONEi = _mm_set1_epi32(1);
	const __m128 tfxWIDEONE = _mm_set1_ps(1.f);
	const __m128 tfxWIDE255 = _mm_set1_ps(255.f);
	const __m128 tfxWIDEZERO = _mm_set1_ps(0.f);
	const __m128 tfxWIDETHIRTYTWO = _mm_set1_ps(32.f);
	const __m128i tfxWIDEFF = _mm_set1_epi32(0xFF);
	const __m128 tfxPWIDESIX = _mm_set_ps1(0.6f);

	typedef union {
		__m128i m;
		int a[4];
	} tfxWideArrayi;

	typedef union {
		__m128 m;
		float a[4];
	} tfxWideArray;

#endif

	typedef __m128 tfx128;
	typedef __m128i tfx128i;

	typedef union {
		__m128i m;
		int a[4];
	} tfx128iArray;

	typedef union {
		__m128i m;
		tfxU64 a[2];
	} tfx128iArray64;

	typedef union {
		__m128 m;
		float a[4];
	} tfx128Array;

	//simd floor function thanks to Stephanie Rancourt: http://dss.stephanierct.com/DevBlog/?p=8
	inline tfx128 tfxFloor128(const tfx128& x) {
		//__m128i v0 = _mm_setzero_si128();
		//__m128i v1 = _mm_cmpeq_epi32(v0, v0);
		//__m128i ji = _mm_srli_epi32(v1, 25);
		//__m128 j = *(__m128*)&_mm_slli_epi32(ji, 23); //create vector 1.0f
		//I'm not entirely sure why original code had above lines to create a vector of 1.f. It seems to me that the below works fine
		//Worth noting that we only need to floor small numbers for the noise algorithm so can get away with this function.
		__m128 j = _mm_set1_ps(1.f); //create vector 1.0f
		__m128i i = _mm_cvttps_epi32(x);
		__m128 fi = _mm_cvtepi32_ps(i);
		__m128 igx = _mm_cmpgt_ps(fi, x);
		j = _mm_and_ps(igx, j);
		return _mm_sub_ps(fi, j);
	}

	inline tfxWideFloat tfxWideAbs(tfxWideFloat v) {
		return tfxWideAnd(tfxWideCast(tfxWideShiftRight(tfxWideSetSinglei(-1), 1)), v);
	}

	inline tfxWideInt tfxWideAbsi(tfxWideInt v) {
		return tfxWideAndi(tfxWideShiftRight(tfxWideSetSinglei(-1), 1), v);
	}

	//----------------------------------------------------------
	//enums/state_flags

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
		tfxTranslationPreset,
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
		tfxGraphCategory_transform,
		tfxGraphCategory_property,
		tfxGraphCategory_base,
		tfxGraphCategory_variation,
		tfxGraphCategory_overtime
	};


#define tfxGlobalCount  15
#define	tfxPropertyCount  9
#define	tfxBaseCount  8
#define	tfxVariationCount  9
#define	tfxOvertimeCount  16
#define	tfxTransformCount  6

#define tfxGlobalStart 0
#define	tfxPropertyStart tfxGlobalCount
#define	tfxBaseStart (tfxPropertyStart + tfxPropertyCount)
#define	tfxVariationStart (tfxBaseStart + tfxBaseCount)
#define	tfxOvertimeStart (tfxVariationStart + tfxVariationCount)
#define	tfxTransformStart (tfxOvertimeStart + tfxOvertimeCount)

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
		tfxGlobal_emitter_width,
		tfxGlobal_emitter_height,
		tfxGlobal_emitter_depth,

		tfxProperty_emission_pitch,
		tfxProperty_emission_yaw,
		tfxProperty_emission_range,
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

		tfxTransform_roll,
		tfxTransform_pitch,
		tfxTransform_yaw,
		tfxTransform_translate_x,
		tfxTransform_translate_y,
		tfxTransform_translate_z,
		tfxGraphMaxIndex,
	};

	inline bool IsGraphTransformRotation(tfxGraphType type) {
		return type == tfxTransform_roll || type == tfxTransform_pitch || type == tfxTransform_yaw;
	}

	inline bool IsGraphEmitterDimension(tfxGraphType type) {
		return type == tfxProperty_emitter_width || type == tfxProperty_emitter_height || type == tfxProperty_emitter_depth;
	}

	inline bool IsGraphTranslation(tfxGraphType type) {
		return type == tfxTransform_translate_x || type == tfxTransform_translate_y || type == tfxTransform_translate_z;
	}

	inline bool IsGraphEmission(tfxGraphType type) {
		return type == tfxProperty_emission_pitch || type == tfxProperty_emission_yaw;
	}

	inline bool IsGraphParticleSize(tfxGraphType type) {
		return	type == tfxBase_width || type == tfxBase_height ||
			type == tfxVariation_width || type == tfxVariation_height ||
			type == tfxOvertime_width || type == tfxOvertime_height;
	}

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
		tfxEllipse,
		tfxCylinder,
		tfxIcosphere
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
		tfxGreyScale
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
	//The values of existing enums below must never change or older files won't load anymore!
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
		tfxStartStage,
		tfxEndStage
	};

	typedef tfxU32 tfxEmitterPropertyFlags;
	typedef tfxU32 tfxEffectPropertyFlags;
	typedef tfxU32 tfxVectorFieldFlags;
	typedef tfxU32 tfxParticleFlags;
	typedef tfxU32 tfxEmitterStateFlags;
	typedef tfxU32 tfxEffectStateFlags;
	typedef tfxU32 tfxParticleControlFlags;
	typedef tfxU32 tfxAttributeNodeFlags;
	typedef tfxU32 tfxAngleSettingFlags;
	typedef tfxU32 tfxParticleManagerFlags;
	typedef tfxU32 tfxErrorFlags;
	typedef tfxU32 tfxEffectCloningFlags;
	typedef tfxU32 tfxAnimationFlags;

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
		tfxBillboarding_disabled = 1 << 0,
		tfxBillboarding_disabled_align = 1 << 1,
		tfxBillboarding_align = 1 << 2
	};

	enum tfxParticleManagerFlags_ {
		tfxEffectManagerFlags_none = 0,
		tfxEffectManagerFlags_disable_spawning = 1,
		tfxEffectManagerFlags_force_capture = 2,
		tfxEffectManagerFlags_use_compute_shader = 1 << 3,
		tfxEffectManagerFlags_order_by_depth = 1 << 4,
		tfxEffectManagerFlags_guarantee_order = 1 << 5,
		tfxEffectManagerFlags_update_base_values = 1 << 6,
		tfxEffectManagerFlags_dynamic_sprite_allocation = 1 << 7,
		tfxEffectManagerFlags_3d_effects = 1 << 8,
		tfxEffectManagerFlags_unordered = 1 << 9,
		tfxEffectManagerFlags_ordered_by_age = 1 << 10,
		tfxEffectManagerFlags_update_age_only = 1 << 11,
		tfxEffectManagerFlags_single_threaded = 1 << 12,
		tfxEffectManagerFlags_double_buffer_sprites = 1 << 13
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
		tfxAngleSettingFlags_random_roll = 1 << 1,											//Chose a random angle at spawn time/state_flags
		tfxAngleSettingFlags_specify_roll = 1 << 2,											//Specify the angle at spawn time
		tfxAngleSettingFlags_align_with_emission = 1 << 3,									//Align the particle with the emission direction only
		tfxAngleSettingFlags_random_pitch = 1 << 4,											//3d mode allows for rotating pitch and yaw when not using billboarding (when particle always faces the camera)
		tfxAngleSettingFlags_random_yaw = 1 << 5,
		tfxAngleSettingFlags_specify_pitch = 1 << 6,
		tfxAngleSettingFlags_specify_yaw = 1 << 7
	};

	//All the state_flags needed by the ControlParticle function put into one enum to save space
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
		tfxEffectPropertyFlags_use_keyframes = 1 << 4
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
		tfxEmitterPropertyFlags_grid_spawn_random = 1 << 26,				//Spawn on grid points but randomly rather then in sequence
		tfxEmitterPropertyFlags_area_open_ends = 1 << 27,					//Only sides of the area/cylinder are spawned on when fill area is not checked
		tfxEmitterPropertyFlags_exclude_from_hue_adjustments = 1 << 28,		//Emitter will be excluded from effect hue adjustments if this flag is checked
		tfxEmitterPropertyFlags_enabled = 1 << 29							//The emitter is enabled or not, meaning it will or will not be added the particle manager with AddEffect
	};

	enum tfxParticleFlags_ : unsigned char {
		tfxParticleFlags_none = 0,
		tfxParticleFlags_fresh = 1 << 0,									//Particle has just spawned this frame	
		tfxParticleFlags_capture_after_transform = 1 << 3,					//Particle will be captured after a transfrom, used for traversing lines and looping back to the beginning to avoid lerping imbetween
		tfxParticleFlags_remove = 1 << 4,									//Particle will be removed this or next frame
		tfxParticleFlags_has_velocity = 1 << 5,								//Flagged if the particle is currently moving
		tfxParticleFlags_has_sub_effects = 1 << 6,							//Flagged if the particle has sub effects
	};

	enum tfxEmitterStateFlags_ : unsigned int {
		tfxEmitterStateFlags_none = 0,
		tfxEmitterStateFlags_random_color = 1 << 0,
		tfxEmitterStateFlags_relative_position = 1 << 1,					//Keep the particles position relative to the current position of the emitter
		tfxEmitterStateFlags_relative_angle = 1 << 2,						//Keep the angle of the particles relative to the current angle of the emitter
		tfxEmitterStateFlags_stop_spawning = 1 << 3,						//Tells the emitter to stop spawning
		tfxEmitterStateFlags_remove = 1 << 4,								//Tells the effect/emitter to remove itself from the particle manager immediately
		tfxEmitterStateFlags_unused1 = 1 << 5,								//the emitter is enabled. **moved to property state_flags**
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
		tfxEmitterStateFlags_is_line_loop_or_kill = 1 << 18,
		tfxEmitterStateFlags_is_area = 1 << 19,
		tfxEmitterStateFlags_no_tween = 1 << 20,
		tfxEmitterStateFlags_align_with_velocity = 1 << 21,
		tfxEmitterStateFlags_is_sub_emitter = 1 << 22,
		tfxEmitterStateFlags_has_noise = 1 << 23
	};

	enum tfxEffectStateFlags_ : unsigned int {
		tfxEffectStateFlags_none = 0,
		tfxEffectStateFlags_stop_spawning = 1 << 3,							//Tells the emitter to stop spawning
		tfxEffectStateFlags_remove = 1 << 4,								//Tells the effect/emitter to remove itself from the particle manager immediately
		tfxEffectStateFlags_retain_matrix = 1 << 6,							//Internal flag about matrix usage
		tfxEffectStateFlags_no_tween_this_update = 1 << 7,					//Internal flag generally, but you could use it if you want to teleport the effect to another location
		tfxEffectStateFlags_override_overal_scale = 1 << 8,					//Flagged when the over scale is overridden with SetEffectOveralScale
		tfxEffectStateFlags_override_orientiation = 1 << 9,					//Flagged when any of the effect angles are overridden
		tfxEffectStateFlags_override_size_multiplier = 1 << 10,				//Flagged when any of the effect size multipliers are overridden
		tfxEffectStateFlags_no_tween = 1 << 20
	};

	enum tfxVectorFieldFlags_ : unsigned char {
		tfxVectorFieldFlags_none = 0,
		tfxVectorFieldFlags_repeat_horizontal = 1 << 0,						//Field will repeat horizontally
		tfxVectorFieldFlags_repeat_vertical = 1 << 1						//Field will repeat vertically
	};

	enum tfxAttributeNodeFlags_ {
		tfxAttributeNodeFlags_none = 0,
		tfxAttributeNodeFlags_is_curve = 1 << 0,
		tfxAttributeNodeFlags_is_left_curve = 1 << 1,
		tfxAttributeNodeFlags_is_right_curve = 1 << 2,
		tfxAttributeNodeFlags_curves_initialised = 1 << 3
	};

	enum tfxAnimationFlags_ {
		tfxAnimationFlags_none = 0,
		tfxAnimationFlags_loop = 1 << 0,
		tfxAnimationFlags_seamless = 1 << 1,
		tfxAnimationFlags_needs_recording = 1 << 2,
		tfxAnimationFlags_export_with_transparency = 1 << 3,
		tfxAnimationFlags_auto_set_length = 1 << 4,
		tfxAnimationFlags_orthographic = 1 << 5
	};

	//-----------------------------------------------------------
	//Constants

	const float tfxMIN_FLOAT = -2147483648.f;
	const float tfxMAX_FLOAT = 2147483647.f;
	const tfxU32 tfxMAX_UINT = 4294967295;
	const int tfxMAX_INT = INT_MAX;
	const int tfxMIN_INT = INT_MIN;
	const tfxS64 tfxMAX_64i = LLONG_MAX;
	const tfxS64 tfxMIN_64i = LLONG_MIN;
	const tfxU64 tfxMAX_64u = ULLONG_MAX;
#if defined(__x86_64__) || defined(_M_X64)
	typedef tfxU64 tfxAddress;
#else
	typedef tfxU32 tfxAddress;
#endif

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
	extern tfxWideFloat tfxUPDATE_TIME_WIDE;
	extern float tfxFRAME_LENGTH;
	extern tfxWideFloat tfxFRAME_LENGTH_WIDE;

	//Look up frequency determines the resolution of graphs that are compiled into look up arrays.
	static float tfxLOOKUP_FREQUENCY = 10.f;
	//Overtime frequency is for lookups that will vary in length depending on the lifetime of the particle. It should generally be a higher resolution than the base graphs
	static float tfxLOOKUP_FREQUENCY_OVERTIME = 1.f;

	//Look up frequency determines the resolution of graphs that are compiled into look up arrays.
	static tfxWideFloat tfxLOOKUP_FREQUENCY_WIDE = tfxWideSetSingle(10.f);
	//Overtime frequency is for lookups that will vary in length depending on the lifetime of the particle. It should generally be a higher resolution than the base graphs
	static tfxWideFloat tfxLOOKUP_FREQUENCY_OVERTIME_WIDE = tfxWideSetSingle(1.f);

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
		inline void         reserve(tfxU32 new_capacity) { if (new_capacity <= capacity) return; 
			tfxMemoryTrackerPair* new_data = (tfxMemoryTrackerPair*)malloc((size_t)new_capacity * sizeof(tfxMemoryTrackerPair)); 
			assert(new_data);	//Unable to allocate memory for new_data
			if (data) { 
				memcpy(new_data, data, (size_t)current_size * sizeof(tfxMemoryTrackerPair)); 
				free(data); 
			} 
			data = new_data; 
			capacity = new_capacity; 
		}

		inline tfxMemoryTrackerPair*           erase(const tfxMemoryTrackerPair* it) { assert(it >= data && it < data + current_size); const ptrdiff_t off = it - data; memmove(data + off, data + off + 1, ((size_t)current_size - (size_t)off - 1) * sizeof(tfxMemoryTrackerPair)); current_size--; return data + off; }
		inline tfxMemoryTrackerPair*           erase(const tfxMemoryTrackerPair* it, const tfxMemoryTrackerPair* it_last) { assert(it >= data && it < data + current_size && it_last > it && it_last <= data + current_size); const ptrdiff_t count = it_last - it; const ptrdiff_t off = it - data; memmove(data + off, data + off + count, ((size_t)current_size - (size_t)off - count) * sizeof(tfxMemoryTrackerPair)); current_size -= (tfxU32)count; return data + off; }
		inline tfxMemoryTrackerPair*           insert(const tfxMemoryTrackerPair* it, const tfxMemoryTrackerPair& v) { assert(it >= data && it <= data + current_size); const ptrdiff_t off = it - data; if (current_size == capacity) reserve(_grow_capacity(current_size + 1)); if (off < (ptrdiff_t)current_size) memmove(data + off + 1, data + off, ((size_t)current_size - (size_t)off) * sizeof(tfxMemoryTrackerPair)); new((void*)(data + off)) tfxMemoryTrackerPair(v); current_size++; return data + off; }

	};

	struct tfxMemoryTrackerLog {

		tfxLogList log;
		std::mutex insert_mutex;
		char last_entry[64];

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
			memcpy(last_entry, value.name, 64);
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
#define tfxINIT_VEC_NAME if(name[0] < 41 || name[0] > 90) { memset(name, '\0', 64); name[0] = 'X'; }
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
		tfxU64 result = _InterlockedExchange64((__int64 volatile*)value, new_value);
		return result;
	}

	inline tfxU64 AtomicAdd64(tfxU64 volatile *value, tfxU64 amount_to_add) {
		tfxU64 result = _InterlockedExchangeAdd64((__int64 volatile*)value, amount_to_add);
		return result;
	}

	inline tfxU64 AtomicIncrement64(tfxU64 volatile *value) {
		return InterlockedIncrement64((__int64 volatile*)value);
	}

	inline tfxU32 AtomicIncrement32(LONG volatile *value) {
		return InterlockedIncrement((LONG volatile*)value);
	}

	inline tfxU32 AtomicExchange32(tfxU32 volatile *value, tfxU32 new_value) {
		tfxU32 result = _InterlockedExchange((LONG*)value, new_value);
		return result;
	}

	inline tfxU32 AtomicAdd32(tfxU32 volatile *value, tfxU32 amount_to_add) {
		tfxU32 result = _InterlockedExchangeAdd((LONG*)value, amount_to_add);
		return result;
	}

	inline tfxU32 AtomicExchangeCompare(tfxU32 volatile *value, tfxU32 exchange, tfxU32 compare) {
		tfxU32 result = InterlockedCompareExchange((LONG volatile *)value, exchange, compare);
		return result;
	}

#define tfxArrayCount(Array) (sizeof(Array) / sizeof((Array)[0]))

	//Storage
	//Credit to ocornut https://github.com/ocornut/imgui/commits?author=ocornut for tfxvec
	//std::vector replacement with some extra stuff and tweaks specific to TimelineFX
	template<typename T>
	struct tfxvec {
#ifdef tfxTRACK_MEMORY
		char name[64];
#endif
		tfxU32 current_size;
		tfxU32 capacity;
		tfxU32 volatile locked;
		T* data;

		// Provide standard typedefs but we don't use them ourselves.
		typedef T                   value_type;
		typedef value_type*         iterator;
		typedef const value_type*   const_iterator;

		inline tfxvec() { locked = false; current_size = capacity = 0; data = NULL; tfxINIT_VEC_NAME; }
		inline tfxvec(const char *name_init) { locked = false; current_size = capacity = 0; data = NULL; tfxINIT_VEC_NAME_INIT(name_init); }
		inline tfxvec(const tfxvec<T> &src) { locked = false; current_size = capacity = 0; data = NULL; tfxINIT_VEC_NAME_SRC_COPY; resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); }
		inline tfxvec<T>& operator=(const tfxvec<T>& src) { clear(); resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); return *this; }
		inline ~tfxvec() { if (data) { tfxFREE(data) }; data = NULL; current_size = capacity = 0; }

		inline bool			empty() { return current_size == 0; }
		inline tfxU32		size() { return current_size; }
		inline const tfxU32	size() const { return current_size; }
		inline tfxU32		size_in_bytes() { return current_size * sizeof(T); }
		inline const tfxU32	size_in_bytes() const { return current_size * sizeof(T); }
		inline T&           operator[](tfxU32 i) { return data[i]; }
		inline const T&     operator[](tfxU32 i) const { assert(i < current_size); return data[i]; }
		inline T&           ts_at(tfxU32 i) { while (locked > 0); return data[i]; }

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
		inline void         reserve(tfxU32 new_capacity) {
			if (new_capacity <= capacity)
				return;
			T* new_data = (T*)tfxALLOCATE(name, new_data, (size_t)new_capacity * sizeof(T)); 
			assert(new_data);	//Unable to allocate memory. todo: better handling
			if (data) { memcpy(new_data, data, (size_t)current_size * sizeof(T)); tfxFREE(data); } data = new_data; capacity = new_capacity;
		}

		inline T&	        grab() {
			if (current_size == capacity) reserve(_grow_capacity(current_size + 1));
			current_size++;
			return data[current_size - 1];
		}
		inline tfxU32        locked_push_back(const T& v) {
			//suspect, just use a mutex instead
			while (InterlockedCompareExchange((LONG volatile*)&locked, 1, 0) > 1);
			if (current_size == capacity)
				reserve(_grow_capacity(current_size + 1));
			new((void*)(data + current_size)) T(v);
			tfxU32 index = current_size++;
			InterlockedExchange((LONG volatile*)&locked, 0);
			return index;
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
		inline void			zero() { assert(capacity > 0); memset(data, 0, capacity * sizeof(T)); }
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
		void *user_data;
		void(*resize_callback)(tfxring<T> *ring, T *new_data, void *user_data);
		unsigned int bank_index;
		tfxring *pair;

		inline tfxring() : resize_callback(nullptr), user_data(NULL), pair(nullptr) { start_index = current_size = capacity = last_bump = bank_index = 0; data = NULL; tfxINIT_VEC_NAME; }
		inline tfxring(const char *name_init) : resize_callback(nullptr), user_data(NULL), pair(nullptr) { start_index = current_size = capacity = last_bump = 0; data = NULL; tfxINIT_VEC_NAME_INIT(name_init); }
		inline tfxring(unsigned int qty) : resize_callback(nullptr), user_data(NULL), pair(nullptr) { start_index = current_size = capacity = last_bump = 0; data = NULL; reserve(qty); tfxINIT_VEC_NAME; }
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
		inline T&	        pop_back() { assert(current_size > 0); current_size--; return data[last_index()]; }
		inline T&	        pop_front() { assert(current_size > 0); tfxU32 front_index = start_index++; start_index %= capacity; current_size--; return data[front_index]; }

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
		inline void         reserve(tfxU32 new_capacity, bool use_callback = true) {
			if (new_capacity <= capacity) return;
			T* new_data = (T*)tfxALLOCATE(name, new_data, (size_t)new_capacity * sizeof(T));
			assert(new_data);	//Unable to allocate memory. todo: better handling
			if (data) {
				if (last_index() < start_index) {
					memcpy(new_data, data + start_index, (size_t)(capacity - start_index) * sizeof(T));
					memcpy(new_data + (capacity - start_index), data, (size_t)(start_index) * sizeof(T));
				}
				else {
					memcpy(new_data, data + start_index, (size_t)current_size * sizeof(T));
				}
				if (resize_callback && use_callback)
					resize_callback(this, new_data, user_data);
				tfxFREE(data);
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

#ifndef tfxSTACK_SIZE
#define tfxSTACK_SIZE tfxMegabyte(2)
#endif

#ifndef tfxMT_STACK_SIZE
#define tfxMT_STACK_SIZE tfxMegabyte(4)
#endif

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

		tfxMemoryArena() { data = end_of_allocated = NULL; memory_remaining = total_memory = 0; }

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

		tfxMemoryArenaManager() :
			arenas(tfxCONSTRUCTOR_VEC_INIT("Memory Arena arenas")),
			blocks(tfxCONSTRUCTOR_VEC_INIT("Memory Arena blocks")),
			free_blocks(tfxCONSTRUCTOR_VEC_INIT("Memory Arena free_blocks")),
			arena_size(0),
			size_diff_threshold(0)
		{}

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
			block->clear();
			free_blocks.push_back(index);
		}

		inline void FreeBlock(tfxU32 block) {
			if (block != tfxINVALID) {
				blocks[block].clear();
				if (block == blocks.current_size - 1) {
					arenas[blocks[block].arena_index].end_of_allocated = blocks[block].data;
					arenas[blocks[block].arena_index].memory_remaining += blocks[block].capacity_in_bytes();
					blocks.pop();
				}
				else {
					free_blocks.push_back(block);
				}
			}
		}

		inline tfxU32 FreeBlocks(tfxU32 block) {
			if (block >= blocks.current_size) return 0;
			FreeBlock(block);
			tfxU32 freed_count = 1;
			while (blocks[block].next_block != tfxINVALID) {
				tfxU32 prev_block = block;
				block = blocks[block].next_block;
				FreeBlock(block);
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
			if (starting_block == tfxINVALID) return tfxINVALID;
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
			memcpy(blocks[to].data, blocks[from].data, (tfxAddress)blocks[from].capacity * blocks[from].unit_size);
			auto &src = blocks[from];
			auto &dst = blocks[to];
			dst.current_size = src.current_size;
			dst.end_ptr = (char*)dst.data + ((tfxAddress)dst.unit_size * dst.current_size);
		}

		inline void CopyBlockToBlock(tfxMemoryBucket *from, tfxMemoryBucket *to) {
			assert(from->capacity && from->capacity <= to->capacity);		//must have valid capacities
			memcpy(to->data, from->data, (tfxAddress)from->capacity * from->unit_size);
			to->current_size = from->current_size;
			to->end_ptr = (char*)to->data + ((tfxAddress)to->unit_size * to->current_size);
		}

		inline tfxMemoryArena *AddArena() {
			tfxMemoryArena arena;
			arena.total_memory = arena_size;
			arena.memory_remaining = arena_size;
			arena.data = tfxALLOCATE(0, 0, arena_size);
			assert(arena.data); //Unable to allocate memory. Todo: proper handling of out of memory
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

	inline void CopyBlockToBlock(tfxMemoryArenaManager &from_allocator, tfxMemoryArenaManager &to_allocator, tfxU32 from, tfxU32 to) {
		assert(from_allocator.blocks[from].capacity && from_allocator.blocks[from].capacity <= to_allocator.blocks[to].capacity);		//must have valid capacities
		memcpy(to_allocator.blocks[to].data, from_allocator.blocks[from].data, (tfxAddress)from_allocator.blocks[from].capacity * from_allocator.blocks[from].unit_size);
		auto &src = from_allocator.blocks[from];
		auto &dst = to_allocator.blocks[to];
		dst.current_size = src.current_size;
		dst.end_ptr = (char*)dst.data + ((tfxAddress)dst.unit_size * dst.current_size);
	}

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
		void* new_data = tfxALLOCATE(0, 0, size_in_bytes);
		assert(new_data);	//Unable to allocate memory. Todo: better handling
		allocator.data = new_data;
		allocator.end_of_allocated = allocator.data;
		allocator.memory_remaining = size_in_bytes;
		allocator.total_memory = size_in_bytes;
		memset((void*)allocator.data, 0, size_in_bytes);
		return allocator;
	}

	struct tfxSoAData {
		void *ptr = NULL;
		size_t offset = 0;
		size_t unit_size = 0;
	};

	//A buffer designed to contain structs of arrays. If the arrays need to grow then a new memory block is made and all copied over
	//together. All arrays in the struct will be the same capacity but can all have different unit sizes/types.
	//In order to use this you need to first prepare the buffer by calling AddStructArray for each struct member of the SoA you're setting up. 
	//All members must be of the same struct.
	//Then call FinishSoABufferSetup to create the memory for the struct of arrays with an initial reserve amount.
	struct tfxSoABuffer {
		size_t current_arena_size = 0;		//The current size of the arena that contains all the arrays
		size_t struct_size = 0;				//The size of the struct (each member unit size added up)
		tfxU32 current_size = 0;			//current size of each array
		tfxU32 start_index = 0;				//Start index if you're using the buffer as a ring buffer
		tfxU32 last_bump = 0;				//the amount the the start index was bumped by the last time Bump was called
		tfxU32 capacity = 0;				//capacity of each array
		tfxU32 block_size = tfxDataWidth;	//Keep the capacity to the nearest block size
		tfxvec<tfxSoAData> array_ptrs;		//Container for all the pointers into the arena
		void *user_data = NULL;
		void(*resize_callback)(tfxSoABuffer *ring, tfxU32 new_index_start) = NULL;
		void *struct_of_arrays = NULL;				//Pointer to the struct of arrays. Important that this is a stable pointer! Set with FinishSoABufferSetup
		void *data = NULL;					//Pointer to the area in memory that contains all of the array data	
	};

	//Get the amount of free space in the buffer
	static inline tfxU32 FreeSpace(tfxSoABuffer *buffer) {
		return buffer->capacity - buffer->current_size;
	}

	//Get the index based on the buffer being a ring buffer
	static inline tfxU32 GetCircularIndex(tfxSoABuffer *buffer, tfxU32 index) {
		return (buffer->start_index + index) % buffer->capacity;
	}

	//Get the index based on the buffer being a ring buffer
	static inline tfxU32 GetAbsoluteIndex(tfxSoABuffer *buffer, tfxU32 circular_index) {
		return buffer->capacity - (circular_index % buffer->capacity);
	}

	//Add an array to a SoABuffer. parse in the size of the data type and the offset to the member variable within the struct.
	//You must add all the member veriables in the struct before calling FinishSoABufferSetup
	static inline void AddStructArray(tfxSoABuffer *buffer, size_t unit_size, size_t offset) {
		tfxSoAData data;
		data.unit_size = unit_size;
		data.offset = offset;
		buffer->array_ptrs.push_back(data);
	}

	//Once you have called AddStructArray for all your member variables you must call this function in order to 
	//set up the memory for all your arrays. One block of memory will be created and all your arrays will be line up
	//inside the space
	static inline void FinishSoABufferSetup(tfxSoABuffer *buffer, void *struct_of_arrays, tfxU32 reserve_amount) {
		assert(buffer->data == NULL && buffer->array_ptrs.current_size > 0);
		assert(reserve_amount > buffer->block_size);		//reserve amount must be greater than the block_size
		for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
			buffer->struct_size += buffer->array_ptrs[i].unit_size;
		}
		reserve_amount = (reserve_amount / buffer->block_size + 1) * buffer->block_size;
		buffer->current_arena_size = reserve_amount * buffer->struct_size;
		buffer->data = tfxALLOCATE(0, 0, buffer->current_arena_size);
		assert(buffer->data);	//Unable to allocate memory. Todo: better handling
		memset(buffer->data, 0, buffer->current_arena_size);
		buffer->capacity = reserve_amount;
		buffer->struct_of_arrays = struct_of_arrays;
		size_t running_offset = 0;
		for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
			buffer->array_ptrs[i].ptr = (char*)buffer->data + running_offset;
			memcpy((char*)buffer->struct_of_arrays + buffer->array_ptrs[i].offset, &buffer->array_ptrs[i].ptr, sizeof(void*));
			running_offset += buffer->array_ptrs[i].unit_size * buffer->capacity;
		}
		if (buffer->resize_callback) {
			buffer->resize_callback(buffer, 0);
		}
	}

	//Call this function to increase the capacity of all the arrays in the buffer. Data that is already in the arrays is preserved.
	static inline void GrowArrays(tfxSoABuffer *buffer, tfxU32 first_new_index, tfxU32 new_size = 0, bool keep_data = true) {
		assert(buffer->capacity);			//buffer must already have a capacity!
		tfxU32 new_capacity = 0;
		if (new_size > 0) {
			if (new_size < buffer->capacity)
				return;
			new_capacity = new_size;
		}
		else {
			new_capacity = buffer->current_size > buffer->capacity ? buffer->current_size + buffer->current_size / 2 : buffer->capacity + buffer->capacity / 2;
		}
		new_capacity = (new_capacity / buffer->block_size + 1) * buffer->block_size;
		void *new_data = tfxALLOCATE(0, 0, new_capacity * buffer->struct_size);
		assert(new_data);	//Unable to allocate memory. Todo: better handling
		memset(new_data, 0, new_capacity * buffer->struct_size);
		size_t running_offset = 0;
		if (keep_data) {
			for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
				if (((buffer->start_index + buffer->current_size - 1) % buffer->capacity) < buffer->start_index) {
					size_t start_index = buffer->start_index * buffer->array_ptrs[i].unit_size;
					size_t capacity = buffer->capacity * buffer->array_ptrs[i].unit_size;
					memcpy((char*)new_data + running_offset, (char*)buffer->array_ptrs[i].ptr + start_index, (size_t)(capacity - start_index));
					memcpy((char*)new_data + (capacity - start_index) + running_offset, (char*)buffer->array_ptrs[i].ptr, (size_t)(start_index));
				}
				else {
					memcpy((char*)new_data + running_offset, buffer->array_ptrs[i].ptr, buffer->array_ptrs[i].unit_size * buffer->capacity);
				}
				running_offset += buffer->array_ptrs[i].unit_size * new_capacity;

			}
		}
		void *old_data = buffer->data;

		buffer->data = new_data;
		buffer->capacity = new_capacity;
		buffer->current_arena_size = new_capacity * buffer->struct_size;
		buffer->start_index = 0;
		running_offset = 0;
		for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
			buffer->array_ptrs[i].ptr = (char*)buffer->data + running_offset;
			memcpy((char*)buffer->struct_of_arrays + buffer->array_ptrs[i].offset, &buffer->array_ptrs[i].ptr, sizeof(void*));
			running_offset += buffer->array_ptrs[i].unit_size * buffer->capacity;
		}
		free(old_data);

		if (buffer->resize_callback) {
			buffer->resize_callback(buffer, first_new_index);
		}
	}

	static inline void GrowArraysBySpecificAmount(tfxSoABuffer *buffer, tfxU32 extra_size) {
		assert(extra_size > 0);		//Must specify a size
		assert(buffer->data);		//Must be a valid SoA buffer
		GrowArrays(buffer, buffer->capacity, buffer->capacity + extra_size);
	}

	//Increase current size of a SoA Buffer and grow if necessary.
	static inline void Resize(tfxSoABuffer *buffer, tfxU32 new_size) {
		assert(buffer->data);			//No data allocated in buffer
		buffer->current_size = new_size;
		if (buffer->current_size >= buffer->capacity) {
			GrowArrays(buffer, buffer->capacity);
		}
	}

	//Increase current size of a SoA Buffer and grow if necessary.
	static inline void SetCapacity(tfxSoABuffer *buffer, tfxU32 new_size) {
		assert(buffer->data);			//No data allocated in buffer
		if (new_size >= buffer->capacity) {
			GrowArrays(buffer, buffer->capacity, new_size);
		}
	}

	//Increase current size of a SoA Buffer and grow if grow is true. Returns the last index.
	static inline tfxU32 AddRow(tfxSoABuffer *buffer, bool grow = false) {
		assert(buffer->data);			//No data allocated in buffer
		buffer->current_size++;
		if (grow && buffer->current_size == buffer->capacity) {
			GrowArrays(buffer, buffer->capacity);
		}
		assert(buffer->current_size <= buffer->capacity);	//Capacity of buffer is exceeded, set grow to true or don't exceed the capacity
		return buffer->current_size - 1;
	}

	//Increase current size of a SoA Buffer and grow if grow is true. Returns the last index.
	static inline tfxU32 AddRows(tfxSoABuffer *buffer, tfxU32 amount, bool grow = false) {
		assert(buffer->data);			//No data allocated in buffer
		tfxU32 first_new_index = buffer->current_size;
		buffer->current_size += amount;
		if (grow && buffer->current_size >= buffer->capacity) {
			GrowArrays(buffer, buffer->capacity);
		}
		assert(buffer->current_size < buffer->capacity);	//Capacity of buffer is exceeded, set grow to true or don't exceed the capacity
		return first_new_index;
	}

	//Decrease the current size of a SoA Buffer.
	static inline void PopRow(tfxSoABuffer *buffer, bool grow = false) {
		assert(buffer->data && buffer->current_size > 0);			//No data allocated in buffer
		buffer->current_size--;
	}

	//Bump the start index of the SoA buffer (ring buffer usage)
	static inline void Bump(tfxSoABuffer *buffer) {
		assert(buffer->data && buffer->current_size > 0);			//No data allocated in buffer
		if (buffer->current_size == 0)
			return;
		buffer->start_index++; buffer->start_index %= buffer->capacity; buffer->current_size--;
	}

	//Bump the start index of the SoA buffer (ring buffer usage)
	static inline void Bump(tfxSoABuffer *buffer, tfxU32 amount) {
		assert(buffer->data && buffer->current_size > 0);			//No data allocated in buffer
		if (buffer->current_size == 0)
			return;
		if (amount > buffer->current_size)
			amount = buffer->current_size;
		buffer->start_index += amount;
		buffer->start_index %= buffer->capacity;
		buffer->current_size -= amount;
		buffer->last_bump = amount;
	}

	//Free the SoA buffer
	static inline void FreeSoABuffer(tfxSoABuffer *buffer) {
		buffer->current_arena_size = buffer->current_size = buffer->capacity = 0;
		if (buffer->data)
			free(buffer->data);
		buffer->data = NULL;
		buffer->array_ptrs.free_all();
	}

	//Clear the SoA buffer
	static inline void ClearSoABuffer(tfxSoABuffer *buffer) {
		buffer->current_size = buffer->start_index = 0;
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
		tfxArray(tfxMemoryArenaManager *allocator_init, tfxU32 size) : allocator(allocator_init) { block = NULL; block_index = tfxINVALID; capacity = 0; reserve(size); }

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
		inline void			zero() { assert(capacity > 0); memset(block, 0, capacity * sizeof(T)); }

	};

	template <typename T>
	struct tfxBucketArray {
		tfxMemoryArenaManager *allocator;	//Pointer to the arena manager that handles the memory and blocks of that memory
		tfxU32 block;						//The first block in the allocator
		tfxU32 current_bucket;				//the current bucket for iterating all blocks
		tfxU32 current_size;				//Current size of the bucket array. This will be the total of all buckets if there are more then one
		tfxU32 capacity;					//The total capacity of the bucket array
		tfxU32 size_of_each_bucket;			//The size of each bucket
		tfxU32 volatile locked;

		tfxBucketArray() { allocator = NULL; current_bucket = block = tfxINVALID; size_of_each_bucket = 64; current_size = capacity = size_of_each_bucket = locked = 0; }
		tfxBucketArray(tfxMemoryArenaManager *allocator_init) : allocator(allocator_init) { size_of_each_bucket = 64; current_size = capacity = locked = 0; current_bucket = block = tfxINVALID; }
		tfxBucketArray(tfxMemoryArenaManager *allocator_init, tfxU32 bucket_size) { assert(bucket_size > 1); size_of_each_bucket = bucket_size; allocator = allocator_init; current_size = locked = capacity = 0; current_bucket = block = tfxINVALID; }

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
			if (&src == this) {
				return *this;
			}
			if (!allocator)
				allocator = src.allocator;
			free_all();
			size_of_each_bucket = src.size_of_each_bucket;
			if (src.capacity == 0) {
				return *this;
			}
			assert(reserve(src.capacity / src.size_of_each_bucket));
			current_bucket = block;
			tfxU32 src_block = src.block;
			while (current_bucket != tfxINVALID && src_block != tfxINVALID) {
				CopyBlockToBlock(*src.allocator, *allocator, src_block, current_bucket);
				current_bucket = allocator->blocks[current_bucket].next_block;
				src_block = src.allocator->blocks[src_block].next_block;
			}
			while (current_bucket != tfxINVALID) {
				allocator->blocks[current_bucket].clear();
				current_bucket = allocator->blocks[current_bucket].next_block;
			}
			TrimBuckets();
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
				number_of_buckets = (int)number_of_buckets - block_count;
			}
			for (int i = 0; i < number_of_buckets; ++i) {
				assert(AllocateBucket<T>(*allocator, size_of_each_bucket, block));		//Out of memory!
				capacity += size_of_each_bucket;
			}
			return true;
		}
		inline tfxU32        locked_push_back(const T& v) {
			while (InterlockedCompareExchange((LONG volatile*)&locked, 1, 0) > 1);

			if (current_size == capacity) {
				assert(AllocateBucket<T>(*allocator, size_of_each_bucket, block));		//Out of memory!
				capacity += size_of_each_bucket;
				ResetIteratorIndex();
			}
			tfxMemoryBucket &last_block = allocator->FirstBlockWithSpace(block);
			assert(last_block.current_size < last_block.capacity);
			*(T*)((T*)last_block.data + last_block.current_size++) = v;
			last_block.end_ptr = (T*)last_block.end_ptr + 1;
			tfxU32 index = current_size++;

			InterlockedExchange((LONG volatile*)&locked, 0);
			return index;
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
				size_t move_size = (tfxAddress)size_of_each_bucket - insert_index - (allocator->blocks[index_block].current_size == allocator->blocks[index_block].capacity ? 1 : 0);
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
					if (*_data == v)		//Note that you need to have an operator overload for == for your type if you don't have one or your inserts/erase may crash on compile here
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
					if (*_data == v)		//Note that you need to have an operator overload for == for your type if you don't have one or your inserts/erase may crash on compile here
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
			if (current_bucket == tfxINVALID)
				return true;
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
			if (block == first_empty_block)
				block = tfxINVALID;
		}
	};

	template <typename T>
	struct tfxStack {
		tfxMemoryArenaManager *allocator;			//Pointer to the allocator that manages the memory and blocks of that memory
		T *block;									//Pointer to the data storing the array. This will be somewhere in a tfxMemoryArena
		tfxU32 block_index;							//This index of the block of memory referenced in allocator->blocks
		tfxU32 capacity;							//The total capacity of the array in units of T
		tfxU32 current_size;

		tfxStack() : allocator(NULL) { block = NULL; capacity = current_size = 0; block_index = tfxINVALID; }
		tfxStack(tfxMemoryArenaManager *allocator_init) : allocator(allocator_init) { block = NULL; capacity = current_size = 0; block_index = tfxINVALID; }
		tfxStack(tfxMemoryArenaManager *allocator_init, tfxU32 size) : allocator(allocator_init) { block = NULL; capacity = current_size = 0; reserve(size); block_index = tfxINVALID; }
		~tfxStack() { free(); }

		inline tfxU32		size() { return current_size; }
		inline const tfxU32	size() const { return current_size; }
		inline T&           operator[](tfxU32 i) {
			assert(i < current_size);		//Index is out of bounds
			return block[i];
		}
		inline const T&     operator[](tfxU32 i) const {
			assert(i < current_size);		//Index is out of bounds
			return block[i];
		}

		inline void         free() { if (block != NULL) { capacity = capacity = 0; allocator->FreeBlock(block_index); block = NULL; block_index = tfxINVALID; } }
		inline void			clear() { current_size = 0; }
		inline T*           begin() { return block; }
		inline const T*     begin() const { return block; }
		inline T*           end() { return block + current_size; }
		inline const T*     end() const { return block + current_size; }
		inline T&           back() { return block[current_size - 1]; }
		inline const T&     back() const { return block[current_size - 1]; }
		inline T&           parent() { assert(current_size > 1); return block[current_size - 2]; }
		inline const T&     parent() const { assert(current_size > 1); return block[current_size - 2]; }
		inline bool			empty() { return current_size == 0; }
		inline tfxU32       _grow_capacity(tfxU32 sz) const { tfxU32 new_capacity = capacity ? (capacity + capacity / 2) : 8; return new_capacity > sz ? new_capacity : sz; }
		inline T&			next() {
			if (current_size == capacity)
				assert(resize(_grow_capacity(current_size + 1), true));	//Stack overflow, try increasing the stack size
			new((void*)(block + current_size)) T();
			return block[current_size++];
		}
		inline void			pop() {
			assert(current_size > 0);		//Can't pop back if the stack is empty
			current_size--;
		}
		inline T&			pop_back() {
			assert(current_size > 0);		//Can't pop back if the stack is empty
			current_size--;
			return block[current_size];
		}
		inline T&	        push_back(const T& v) {
			if (current_size == capacity)
				assert(resize(_grow_capacity(current_size + 1), true));	//Stack overflow, try increasing the stack size
			new((void*)(block + current_size)) T(v);
			return block[current_size++];
		}
		inline T&	        push_back_copy(const T& v) {
			if (current_size == capacity)
				assert(resize(_grow_capacity(current_size + 1), true));	//Stack overflow, try increasing the stack size
			memcpy(&block[current_size], &v, sizeof(v));
			return block[current_size++];
		}
		inline bool			reserve(tfxU32 size) {
			assert(allocator);		//Must assign an allocator before doing anything with a tfxStack. Capacity must equal 0
			assert(capacity == 0);	//Capacity must equal 0 before reserving
			assert(size * sizeof(T) < allocator->arena_size);	//The size must fit into an arena size
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
			if (keep_contents && current_block != tfxINVALID)
				allocator->CopyBlockToBlock(current_block, block_index);
			capacity = size;
			block = (T*)allocator->blocks[block_index].data;
			return true;
		}

	};

	//You must called InitialiseTimelineFX() before doing anything!
#define tmpStack(type, name) assert(tfxSTACK_ALLOCATOR.arena_size > 0); tfxStack<type> name(&tfxSTACK_ALLOCATOR)
		//You must called InitialiseTimelineFX() before doing anything!
#define tmpMTStack(type, name) assert(tfxMT_STACK_ALLOCATOR.arena_size > 0); tfxStack<type> name(&tfxMT_STACK_ALLOCATOR)

	template <typename T>
	static inline tfxBucketArray<T> CreateBucketArray(tfxMemoryArena *allocator, tfxU32 bucket_size) {
		tfxBucketArray bucket_array(allocator, bucket_size);
		return bucket_array;
	}

	//A char buffer you can use to load a file into and read from
	//Has no deconstructor so make sure you call FreeAll() when done
	//This is meant for limited usage in timeline fx only and not recommended for use outside!
	struct tfxStream {
		tfxU64 size = 0;
		tfxU64 position = 0;
		char* data = NULL;

		inline tfxStream() { size = position = 0; data = NULL; }
		inline tfxStream(tfxU64 qty) { size = position = 0; data = NULL; Resize(qty); }

		inline bool Read(char* dst, tfxU64 count) {
			if (count + position <= size) {
				memcpy(dst, data + position, count);
				position += count;
				return true;
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

		inline void			FreeAll() { if (data) { size = size = 0; tfxFREE(data); data = NULL; } }
		inline void         Clear() { if (data) { size = 0; } }

		inline void         Resize(tfxU64 new_capacity) {
			if (new_capacity <= size)
				return;
			char* new_data = (char*)tfxALLOCATE("Stream", new_data, (tfxU64)new_capacity * sizeof(char));
			assert(new_data);	//Unable to allocate memory. Todo: better handling
			if (data) {
				memcpy(new_data, data, (tfxU64)size * sizeof(char));
				free(data);
			}
			data = new_data; size = new_capacity;
			position = 0;
		}
		inline void			NullTerminate() { *(data + size) = NULL; }

	};

	//Some multithreading functions - TODO: currently this is windows only, needs linux/mac etc added
	//Might end up just using std::thread but will see how the Mac API is
	//There is a single thread pool created to serve multiple queues. Currently each particle manager that you create will have it's own queue.
	struct tfxWorkQueue;

#define tfxWORKQUEUECALLBACK(name) void name(tfxWorkQueue *queue, void *data)
	typedef tfxWORKQUEUECALLBACK(tfxWorkQueueCallback);

	struct tfxSignaller {
		tfxU32 signal_slots[256];
	};

	struct tfxWorkQueueEntry {
		tfxWorkQueueCallback *call_back = NULL;
		void *data = NULL;
	};

	/*
		in_signal GetSignalSlot();
		AddWorkQueueEntry(dependencies, out_signal);
	*/

	typedef tfxU32 tfxWorkQueueFlags;

	enum tfxWorkQueueFlag_ {
		tfxWorkQueueFlag_none = 0
	};

	extern HANDLE tfxThreadSemaphore;
	extern int tfxNumberOfThreadsInAdditionToMain;

	struct tfxWorkQueue {
		tfxU32 volatile entry_completion_goal = 0;
		tfxU32 volatile entry_completion_count = 0;
		tfxU32 volatile next_read_entry = 0;
		tfxU32 volatile next_write_entry = 0;
		tfxWorkQueueEntry entries[512];
	};

	struct tfxQueueProcessor {
		std::mutex mutex;
		HANDLE empty_semaphore;
		HANDLE full_semaphore;
		tfxU32 count;
		tfxWorkQueue *queues[512];
	};

	extern tfxQueueProcessor tfxThreadQueues;

	inline void InitialiseThreadQueues(tfxQueueProcessor *queues) {
		queues->count = 0;
		queues->empty_semaphore = CreateSemaphoreEx(0, 64, 64, 0, 0, SEMAPHORE_ALL_ACCESS);
		queues->full_semaphore = CreateSemaphoreEx(0, 0, 64, 0, 0, SEMAPHORE_ALL_ACCESS);
		memset(queues->queues, 0, 64 * sizeof(void*));
	}

	inline tfxWorkQueue *tfxGetQueueWithWork(tfxQueueProcessor *thread_processor) {
		WaitForSingleObject(thread_processor->full_semaphore, INFINITE);
		std::unique_lock<std::mutex> lock(thread_processor->mutex);
		tfxWorkQueue *queue = thread_processor->queues[thread_processor->count--];
		ReleaseSemaphore(thread_processor->empty_semaphore, 1, 0);
		return queue;
	}

	inline void tfxPushQueueWork(tfxQueueProcessor *thread_processor, tfxWorkQueue *queue) {
		WaitForSingleObject(thread_processor->empty_semaphore, INFINITE);
		std::unique_lock<std::mutex> lock(thread_processor->mutex);
		thread_processor->queues[thread_processor->count++] = queue;
		ReleaseSemaphore(thread_processor->full_semaphore, 1, 0);
	}

	static bool tfxDoNextWorkQueue(tfxQueueProcessor *queue_processor) {
		bool sleep = false;

		tfxWorkQueue *queue = tfxGetQueueWithWork(queue_processor);

		if (queue) {
			tfxU32 original_read_entry = queue->next_read_entry;
			tfxU32 new_original_read_entry = (original_read_entry + 1) % (tfxU32)tfxArrayCount(queue->entries);

			if (original_read_entry != queue->next_write_entry) {
				tfxU32 index = InterlockedCompareExchange((LONG volatile *)&queue->next_read_entry, new_original_read_entry, original_read_entry);
				if (index == original_read_entry) {
					tfxWorkQueueEntry entry = queue->entries[index];
					entry.call_back(queue, entry.data);
					InterlockedIncrement((LONG volatile *)&queue->entry_completion_count);
				}
			}
			else {
				sleep = true;
			}
		}
		else {
			sleep = false;
		}

		return sleep;
	}

	static bool tfxDoNextWorkQueueEntry(tfxWorkQueue *queue) {
		bool sleep = false;

		tfxU32 original_read_entry = queue->next_read_entry;
		tfxU32 new_original_read_entry = (original_read_entry + 1) % (tfxU32)tfxArrayCount(queue->entries);

		if (original_read_entry != queue->next_write_entry) {
			tfxU32 index = InterlockedCompareExchange((LONG volatile *)&queue->next_read_entry, new_original_read_entry, original_read_entry);
			if (index == original_read_entry) {
				tfxWorkQueueEntry entry = queue->entries[index];
				entry.call_back(queue, entry.data);
				InterlockedIncrement((LONG volatile *)&queue->entry_completion_count);
			}
		}
		else {
			sleep = true;
		}

		return sleep;
	}

	inline void tfxAddWorkQueueEntry(tfxWorkQueue *queue, void *data, tfxWorkQueueCallback call_back) {
		assert(tfxNumberOfThreadsInAdditionToMain > 0);

		tfxU32 new_entry_to_write = (queue->next_write_entry + 1) % (tfxU32)tfxArrayCount(queue->entries);
		while (new_entry_to_write == queue->next_read_entry) {		//Not enough room in work queue
			//We can do this because we're single producer
			tfxDoNextWorkQueueEntry(queue);
		}
		queue->entries[queue->next_write_entry].data = data;
		queue->entries[queue->next_write_entry].call_back = call_back;
		InterlockedIncrement(&queue->entry_completion_goal);

		_WriteBarrier();

		tfxPushQueueWork(&tfxThreadQueues, queue);
		queue->next_write_entry = new_entry_to_write;

		ReleaseSemaphore(tfxThreadSemaphore, 1, 0);
	}

	inline DWORD WINAPI tfxThreadProc(LPVOID lpParameter) {

		tfxQueueProcessor *thread_processor = (tfxQueueProcessor*)lpParameter;

		for (;;) {
			if (tfxDoNextWorkQueue(thread_processor)) {
				//Suspend the thread
				WaitForSingleObjectEx(tfxThreadSemaphore, INFINITE, false);
			}
		}

	}

	static void tfxCompleteAllWork(tfxWorkQueue *queue) {
		tfxWorkQueueEntry entry = {};
		while (queue->entry_completion_goal != queue->entry_completion_count) {
			tfxDoNextWorkQueueEntry(queue);
		}
		queue->entry_completion_count = 0;
		queue->entry_completion_goal = 0;
	}

	inline void tfxInitialiseWorkQueue(tfxWorkQueue *queue) {
		queue->entry_completion_count = 0;
		queue->entry_completion_goal = 0;
		queue->next_read_entry = 0;
		queue->next_write_entry = 0;
	}

	inline bool tfxInitialiseThreads(tfxQueueProcessor *thread_queues) {
		InitialiseThreadQueues(&tfxThreadQueues);

		//todo: create a function to close all the threads 

		tfxThreadSemaphore = CreateSemaphoreEx(0, 0, tfxNumberOfThreadsInAdditionToMain, 0, 0, SEMAPHORE_ALL_ACCESS);
		tfxU32 threads_initialised = 0;
		for (int thread_index = 0; thread_index < tfxNumberOfThreadsInAdditionToMain; ++thread_index) {
			DWORD thread_id;
			HANDLE thread_handle = CreateThread(0, 0, tfxThreadProc, thread_queues, 0, &thread_id);
			if (!thread_handle) {
				tfxNumberOfThreadsInAdditionToMain = threads_initialised;
				return false;
			}
		}
		return true;
	}

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
		union {
			struct { float x, y, z, w; };
			struct { float c0, c1, c2, c3; };
		};

		tfxVec4() { x = y = z = w = 0.f; }
		tfxVec4(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}
		tfxVec4(tfxVec2 vec1, tfxVec2 vec2) : x(vec1.x), y(vec1.y), z(vec2.x), w(vec2.y) {}
		tfxVec4(tfxVec3 vec) : x(vec.x), y(vec.y), z(vec.z), w(0.f) {}
		tfxVec4(tfxVec3 vec, float _w) : x(vec.x), y(vec.y), z(vec.z), w(_w) {}

		inline tfxVec2 xy() { return tfxVec2(x, y); }
		inline tfxVec2 zw() { return tfxVec2(z, w); }
		inline tfxVec3 xyz() { return tfxVec3(x, y, z); }

		inline tfxVec2 xy() const { return tfxVec2(x, y); }
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

	//Wide simd versions of tfxVec2/3
	struct tfxWideVec3 {
		union {
			struct { tfxWideFloat x, y, z; };
			struct { tfxWideFloat pitch, yaw, roll; };
		};

		tfxWideVec3() { x = tfxWideSetSingle(0.f); y = tfxWideSetSingle(0.f); z = tfxWideSetSingle(0.f); }
		tfxWideVec3(tfxWideFloat _x, tfxWideFloat _y, tfxWideFloat _z) { x = _x; y = _y; z = _z; }

		inline tfxWideVec3 operator+(const tfxWideVec3 &v) const { return tfxWideVec3(tfxWideAdd(x, v.x), tfxWideAdd(y, v.y), tfxWideAdd(z, v.z)); }
		inline tfxWideVec3 operator-(const tfxWideVec3 &v) const { return tfxWideVec3(tfxWideSub(x, v.x), tfxWideSub(y, v.y), tfxWideSub(z, v.z)); }
		inline tfxWideVec3 operator*(const tfxWideVec3 &v) const { return tfxWideVec3(tfxWideMul(x, v.x), tfxWideMul(y, v.y), tfxWideMul(z, v.z)); }
		inline tfxWideVec3 operator/(const tfxWideVec3 &v) const { return tfxWideVec3(tfxWideDiv(x, v.x), tfxWideDiv(y, v.y), tfxWideDiv(z, v.z)); }

		inline tfxWideVec3 operator-() const { return tfxWideVec3(tfxWideSub(tfxWideSetSingle(0.f), x), tfxWideSub(tfxWideSetSingle(0.f), y), tfxWideSub(tfxWideSetSingle(0.f), z)); }

		inline void operator-=(const tfxWideVec3 &v) { x = tfxWideSub(x, v.x); y = tfxWideSub(y, v.y); z = tfxWideSub(z, v.z); }
		inline void operator+=(const tfxWideVec3 &v) { x = tfxWideAdd(x, v.x); y = tfxWideAdd(y, v.y); z = tfxWideAdd(z, v.z); }
		inline void operator*=(const tfxWideVec3 &v) { x = tfxWideMul(x, v.x); y = tfxWideMul(y, v.y); z = tfxWideMul(z, v.z); }
		inline void operator/=(const tfxWideVec3 &v) { x = tfxWideDiv(x, v.x); y = tfxWideDiv(y, v.y); z = tfxWideDiv(z, v.z); }

		inline tfxWideVec3 operator+(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfxWideVec3(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v), tfxWideAdd(z, wide_v)); }
		inline tfxWideVec3 operator-(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfxWideVec3(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v), tfxWideAdd(z, wide_v)); }
		inline tfxWideVec3 operator*(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfxWideVec3(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v), tfxWideAdd(z, wide_v)); }
		inline tfxWideVec3 operator/(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfxWideVec3(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v), tfxWideAdd(z, wide_v)); }

		inline tfxWideVec3 operator+(tfxWideFloat v) const { return tfxWideVec3(tfxWideAdd(x, v), tfxWideAdd(y, v), tfxWideAdd(z, v)); }
		inline tfxWideVec3 operator-(tfxWideFloat v) const { return tfxWideVec3(tfxWideAdd(x, v), tfxWideAdd(y, v), tfxWideAdd(z, v)); }
		inline tfxWideVec3 operator*(tfxWideFloat v) const { return tfxWideVec3(tfxWideAdd(x, v), tfxWideAdd(y, v), tfxWideAdd(z, v)); }
		inline tfxWideVec3 operator/(tfxWideFloat v) const { return tfxWideVec3(tfxWideAdd(x, v), tfxWideAdd(y, v), tfxWideAdd(z, v)); }

		inline void operator*=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideMul(x, wide_v); y = tfxWideMul(y, wide_v); z = tfxWideMul(z, wide_v); }
		inline void operator/=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideDiv(x, wide_v); y = tfxWideDiv(y, wide_v); z = tfxWideDiv(z, wide_v); }
		inline void operator+=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideAdd(x, wide_v); y = tfxWideAdd(y, wide_v); z = tfxWideAdd(z, wide_v); }
		inline void operator-=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideSub(x, wide_v); y = tfxWideSub(y, wide_v); z = tfxWideSub(z, wide_v); }

		inline tfxWideFloat Squared() { return tfxWideAdd(tfxWideMul(x, x), tfxWideAdd(tfxWideMul(y, y), tfxWideMul(z, z))); }
	};

	struct tfxWideVec2 {
		tfxWideFloat x, y;

		tfxWideVec2() { x = tfxWideSetSingle(0.f); y = tfxWideSetSingle(0.f); }
		tfxWideVec2(tfxWideFloat _x, tfxWideFloat _y) { x = _x; y = _y; }

		inline tfxWideVec2 operator+(const tfxWideVec2 &v) const { return tfxWideVec2(tfxWideAdd(x, v.x), tfxWideAdd(y, v.y)); }
		inline tfxWideVec2 operator-(const tfxWideVec2 &v) const { return tfxWideVec2(tfxWideSub(x, v.x), tfxWideSub(y, v.y)); }
		inline tfxWideVec2 operator*(const tfxWideVec2 &v) const { return tfxWideVec2(tfxWideMul(x, v.x), tfxWideMul(y, v.y)); }
		inline tfxWideVec2 operator/(const tfxWideVec2 &v) const { return tfxWideVec2(tfxWideDiv(x, v.x), tfxWideDiv(y, v.y)); }

		inline tfxWideVec2 operator-() const { return tfxWideVec2(tfxWideSub(tfxWideSetSingle(0.f), x), tfxWideSub(tfxWideSetSingle(0.f), y)); }

		inline void operator-=(const tfxWideVec2 &v) { x = tfxWideSub(x, v.x); y = tfxWideSub(y, v.y); }
		inline void operator+=(const tfxWideVec2 &v) { x = tfxWideAdd(x, v.x); y = tfxWideAdd(y, v.y); }
		inline void operator*=(const tfxWideVec2 &v) { x = tfxWideMul(x, v.x); y = tfxWideMul(y, v.y); }
		inline void operator/=(const tfxWideVec2 &v) { x = tfxWideDiv(x, v.x); y = tfxWideDiv(y, v.y); }

		inline tfxWideVec2 operator+(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfxWideVec2(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v)); }
		inline tfxWideVec2 operator-(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfxWideVec2(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v)); }
		inline tfxWideVec2 operator*(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfxWideVec2(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v)); }
		inline tfxWideVec2 operator/(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfxWideVec2(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v)); }

		inline tfxWideVec2 operator+(tfxWideFloat v) const { return tfxWideVec2(tfxWideAdd(x, v), tfxWideAdd(y, v)); }
		inline tfxWideVec2 operator-(tfxWideFloat v) const { return tfxWideVec2(tfxWideAdd(x, v), tfxWideAdd(y, v)); }
		inline tfxWideVec2 operator*(tfxWideFloat v) const { return tfxWideVec2(tfxWideAdd(x, v), tfxWideAdd(y, v)); }
		inline tfxWideVec2 operator/(tfxWideFloat v) const { return tfxWideVec2(tfxWideAdd(x, v), tfxWideAdd(y, v)); }

		inline void operator*=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideMul(x, wide_v); y = tfxWideMul(y, wide_v); }
		inline void operator/=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideDiv(x, wide_v); y = tfxWideDiv(y, wide_v); }
		inline void operator+=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideAdd(x, wide_v); y = tfxWideAdd(y, wide_v); }
		inline void operator-=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideSub(x, wide_v); y = tfxWideSub(y, wide_v); }

		inline tfxWideFloat Squared() { return tfxWideAdd(tfxWideMul(x, x), tfxWideMul(y, y)); }
	};

	inline tfxWideVec3 InterpolateWideVec3(tfxWideFloat &tween, tfxWideVec3 &from, tfxWideVec3 &to) {
		return to * tween + from * (tfxWideSub(tfxWIDEONE, tween));
	}

	static inline void ScaleVec4xyz(tfxVec4 &v, float scalar) {
		v.x *= scalar;
		v.y *= scalar;
		v.z *= scalar;
	}

	struct tfxRGBA8 {
		union {
			struct { unsigned char r, g, b, a; };
			struct { tfxU32 color; };
		};

		tfxRGBA8() { r = g = b = a = 0; }
		tfxRGBA8(unsigned char _r, unsigned char _g, unsigned char _b, unsigned char _a) : r(_r), g(_g), b(_b), a(_a) { }
		tfxRGBA8(float _r, float _g, float _b, float _a) : r((char)_r), g((char)_g), b((char)_b), a((char)_a) { }
		tfxRGBA8(tfxU32 _r, tfxU32 _g, tfxU32 _b, tfxU32 _a) : r((char)_r), g((char)_g), b((char)_b), a((char)_a) { }
		tfxRGBA8(tfxRGBA8 _c, char _a) : r(_c.r), g(_c.g), b(_c.b), a((char)_a) { }
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

	const float one_div_255 = 1 / 255.f;
	const float one_div_511 = 1 / 511.f;
	const tfxWideFloat one_div_511_wide = tfxWideSetSingle(1 / 511.f);
	const tfxWideFloat one_div_32k_wide = tfxWideSetSingle(1 / 32767.f);
	#define tfxPACKED_Y_NORMAL_3D 0x1FFFF9FF
	#define tfxPACKED_Y_NORMAL_2D 32767

	struct tfxRGBA {
		float r, g, b, a;
		tfxRGBA() { r = g = b = a = 1.f; }
		tfxRGBA(float _r, float _g, float _b, float _a) : r(_r), g(_g), b(_b), a(_a) { }
		tfxRGBA(tfxRGBA8 c) : r((float)c.r * one_div_255), g((float)c.g * one_div_255), b((float)c.b * one_div_255), a((float)c.a * one_div_255) { }
	};

	inline tfxHSV RGBtoHSV(tfxRGB in)
	{
		tfxHSV      out;
		float      min, max, delta;

		min = in.r < in.g ? in.r : in.g;
		min = min < in.b ? min : in.b;

		max = in.r > in.g ? in.r : in.g;
		max = max > in.b ? max : in.b;

		out.v = max;                                // v
		delta = max - min;
		if (delta < 0.00001f)
		{
			out.s = 0;
			out.h = 0; // undefined, maybe nan?
			return out;
		}
		if (max > 0.0f) { // NOTE: if Max is == 0, this divide would cause a crash
			out.s = (delta / max);                  // s
		}
		else {
			// if max is 0, then r = g = b = 0              
			// s = 0, h is undefined
			out.s = 0.0f;
			out.h = NAN;                            // its now undefined
			return out;
		}
		if (in.r >= max)                           // > is bogus, just keeps compilor happy
			out.h = (in.g - in.b) / delta;        // between yellow & magenta
		else
			if (in.g >= max)
				out.h = 2.0f + (in.b - in.r) / delta;  // between cyan & yellow
			else
				out.h = 4.0f + (in.r - in.g) / delta;  // between magenta & cyan

		out.h *= 60.0f;                              // degrees

		if (out.h < 0.0f)
			out.h += 360.0f;

		return out;
	}


	inline tfxRGB HSVtoRGB(tfxHSV in)
	{
		float      hh, p, q, t, ff;
		long        i;
		tfxRGB      out;

		if (in.s <= 0.0f) {       // < is bogus, just shuts up warnings
			out.r = in.v;
			out.g = in.v;
			out.b = in.v;
			return out;
		}
		hh = in.h;
		if (hh >= 360.0f) hh = 0.0f;
		hh /= 60.0f;
		i = (long)hh;
		ff = hh - i;
		p = in.v * (1.0f - in.s);
		q = in.v * (1.0f - (in.s * ff));
		t = in.v * (1.0f - (in.s * (1.0f - ff)));

		switch (i) {
		case 0:
			out.r = in.v;
			out.g = t;
			out.b = p;
			break;
		case 1:
			out.r = q;
			out.g = in.v;
			out.b = p;
			break;
		case 2:
			out.r = p;
			out.g = in.v;
			out.b = t;
			break;

		case 3:
			out.r = p;
			out.g = q;
			out.b = in.v;
			break;
		case 4:
			out.r = t;
			out.g = p;
			out.b = in.v;
			break;
		case 5:
		default:
			out.r = in.v;
			out.g = p;
			out.b = q;
			break;
		}
		return out;
	}

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

	static inline float DotProduct(const tfxVec2 &a, const tfxVec2 &b)
	{
		return (a.x * b.x + a.y * b.y);
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

	static inline float FastLength(tfxVec2 const &v) {
		return 1.f / tfxSqrt(DotProduct(v, v));
	}

	static inline float FastLength(tfxVec3 const &v) {
		return 1.f / tfxSqrt(DotProduct(v, v));
	}

	static inline tfxVec3 FastNormalizeVec(tfxVec3 const &v) {
		return v * tfxSqrt(DotProduct(v, v));
	}

	static inline tfxVec2 NormalizeVec(tfxVec2 const &v) {
		float length = FastLength(v);
		return tfxVec2(v.x / length, v.y / length);
	}

	static inline tfxVec2 NormalizeVec(tfxVec2 const &v, const float length) {
		return tfxVec2(v.x / length, v.y / length);
	}

	struct tfxMatrix4 {

		tfxVec4 v[4];

		inline void Set2(float aa, float ab, float ba, float bb) {
			v[0].c0 = aa; v[0].c1 = ab;
			v[1].c0 = ba; v[1].c1 = bb;
		}

	};

	struct tfxMatrix4Wide {
		float x[4];
		float y[4];
		float z[4];
		float w[4];
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

		tfx128 in_row[4];
		in_row[0] = _mm_load_ps(&in.v[0].x);
		in_row[1] = _mm_load_ps(&in.v[1].x);
		in_row[2] = _mm_load_ps(&in.v[2].x);
		in_row[3] = _mm_load_ps(&in.v[3].x);

		tfx128 m_row1 = _mm_set_ps(m.v[3].x, m.v[2].x, m.v[1].x, m.v[0].x);
		tfx128 m_row2 = _mm_set_ps(m.v[3].y, m.v[2].y, m.v[1].y, m.v[0].y);
		tfx128 m_row3 = _mm_set_ps(m.v[3].z, m.v[2].z, m.v[1].z, m.v[0].z);
		tfx128 m_row4 = _mm_set_ps(m.v[3].w, m.v[2].w, m.v[1].w, m.v[0].w);

		for (int r = 0; r <= 3; ++r)
		{

			tfx128 row1result = _mm_mul_ps(in_row[r], m_row1);
			tfx128 row2result = _mm_mul_ps(in_row[r], m_row2);
			tfx128 row3result = _mm_mul_ps(in_row[r], m_row3);
			tfx128 row4result = _mm_mul_ps(in_row[r], m_row4);

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

	static inline void mmWideTransformVector(const tfxMatrix4 &mat, tfxWideFloat &x, tfxWideFloat &y, tfxWideFloat &z) {
		tfxWideFloat xr = tfxWideMul(x, tfxWideSetSingle(mat.v[0].c0));
		xr = tfxWideAdd(tfxWideMul(y, tfxWideSetSingle(mat.v[0].c1)), xr);
		xr = tfxWideAdd(tfxWideMul(z, tfxWideSetSingle(mat.v[0].c2)), xr);
		tfxWideFloat yr = tfxWideMul(x, tfxWideSetSingle(mat.v[1].c0));
		yr = tfxWideAdd(tfxWideMul(y, tfxWideSetSingle(mat.v[1].c1)), yr);
		yr = tfxWideAdd(tfxWideMul(z, tfxWideSetSingle(mat.v[1].c2)), yr);
		tfxWideFloat zr = tfxWideMul(x, tfxWideSetSingle(mat.v[2].c0));
		zr = tfxWideAdd(tfxWideMul(y, tfxWideSetSingle(mat.v[2].c1)), zr);
		zr = tfxWideAdd(tfxWideMul(z, tfxWideSetSingle(mat.v[2].c2)), zr);
		x = xr;
		y = yr;
		z = zr;
	}

	static inline void mmWideTransformVector(const tfxMatrix4 &mat, tfxWideFloat &x, tfxWideFloat &y) {
		tfxWideFloat xr = tfxWideMul(x, tfxWideSetSingle(mat.v[0].c0));
		xr = tfxWideAdd(tfxWideMul(y, tfxWideSetSingle(mat.v[1].c0)), xr);
		tfxWideFloat yr = tfxWideMul(x, tfxWideSetSingle(mat.v[0].c1));
		yr = tfxWideAdd(tfxWideMul(y, tfxWideSetSingle(mat.v[1].c1)), yr);
		x = xr;
		y = yr;
	}

	static inline void mmWideTransformVector(const tfxWideFloat *r0c, const tfxWideFloat *r1c, const tfxWideFloat *r2c, tfxWideFloat &x, tfxWideFloat &y, tfxWideFloat &z, tfxWideFloat &mask, tfxWideFloat &xor_mask) {
		tfxWideFloat xr = tfxWideMul(x, r0c[0]);
		xr = tfxWideAdd(tfxWideMul(y, r0c[1]), xr);
		xr = tfxWideAdd(tfxWideMul(z, r0c[2]), xr);
		tfxWideFloat yr = tfxWideMul(x, r1c[0]);
		yr = tfxWideAdd(tfxWideMul(y, r1c[1]), yr);
		yr = tfxWideAdd(tfxWideMul(z, r1c[2]), yr);
		tfxWideFloat zr = tfxWideMul(x, r2c[0]);
		zr = tfxWideAdd(tfxWideMul(y, r2c[1]), zr);
		zr = tfxWideAdd(tfxWideMul(z, r2c[2]), zr);
		x = tfxWideAdd(tfxWideAnd(x, xor_mask), tfxWideAnd(xr, mask));
		y = tfxWideAdd(tfxWideAnd(y, xor_mask), tfxWideAnd(yr, mask));
		z = tfxWideAdd(tfxWideAnd(z, xor_mask), tfxWideAnd(zr, mask));
	}

	static inline tfxVec4 mmTransformVector(const tfxMatrix4 &mat, const tfxVec4 vec) {
		tfxVec4 v;

		tfx128 v4 = _mm_set_ps(vec.w, vec.z, vec.y, vec.x);

		tfx128 mrow1 = _mm_load_ps(&mat.v[0].c0);
		tfx128 mrow2 = _mm_load_ps(&mat.v[1].c0);
		tfx128 mrow3 = _mm_load_ps(&mat.v[2].c0);
		tfx128 mrow4 = _mm_load_ps(&mat.v[3].c0);

		tfx128 row1result = _mm_mul_ps(v4, mrow1);
		tfx128 row2result = _mm_mul_ps(v4, mrow2);
		tfx128 row3result = _mm_mul_ps(v4, mrow3);
		tfx128 row4result = _mm_mul_ps(v4, mrow4);

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

	static inline tfxVec4 mmTransformVector(const tfx128 &row1, const tfx128 &row2, const tfx128 &row3, const tfx128 &row4, const tfxVec4 vec) {
		tfxVec4 v;

		tfx128 v4 = _mm_set_ps(vec.w, vec.z, vec.y, vec.x);

		tfx128 row1result = _mm_mul_ps(v4, row1);
		tfx128 row2result = _mm_mul_ps(v4, row2);
		tfx128 row3result = _mm_mul_ps(v4, row3);
		tfx128 row4result = _mm_mul_ps(v4, row4);

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

		tfx128 v4 = _mm_set_ps(vec.w, vec.z, vec.y, vec.x);

		tfx128 mrow1 = _mm_load_ps(&mat.v[0].x);
		tfx128 mrow2 = _mm_load_ps(&mat.v[1].x);
		tfx128 mrow3 = _mm_load_ps(&mat.v[2].x);
		tfx128 mrow4 = _mm_load_ps(&mat.v[3].x);

		tfx128 row1result = _mm_mul_ps(v4, mrow1);
		tfx128 row2result = _mm_mul_ps(v4, mrow2);
		tfx128 row3result = _mm_mul_ps(v4, mrow3);
		tfx128 row4result = _mm_mul_ps(v4, mrow4);

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

	static inline int Clamp(int lower, int upper, int value) {
		if (value < lower) return lower;
		if (value > upper) return upper;
		return value;
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
		tfxVec3 converted = v * 511.f;
		tfxUInt10bit result;
		result.pack = 0;
		result.data.x = (tfxU32)converted.z;
		result.data.y = (tfxU32)converted.y;
		result.data.z = (tfxU32)converted.x;
		result.data.w = extra;
		return result.pack;
	}

	inline tfxU32 Pack10bitUnsigned(tfxVec3 const &v) {
		tfxVec3 converted = v * 511.f + 511.f;
		tfxUInt10bit result;
		result.pack = 0;
		result.data.x = (tfxU32)converted.z;
		result.data.y = (tfxU32)converted.y;
		result.data.z = (tfxU32)converted.x;
		result.data.w = 0;
		return result.pack;
	}

	inline tfxU32 Pack16bit(float x, float y) {
		union
		{
			signed short in[2];
			tfxU32 out;
		} u;

		x = x * 32767.f;
		y = y * 32767.f;

		u.in[0] = (signed short)x;
		u.in[1] = (signed short)y;

		return u.out;
	}

	inline tfxU32 Pack16bitUnsigned(float x, float y) {
		union
		{
			struct {
				tfxU32 x : 16;
				tfxU32 y : 16;
			} data;
			tfxU32 out;
		} u;

		x = x * 32767.f + 32767.f;
		y = y * 32767.f + 32767.f;

		u.data.x = (tfxU32)x;
		u.data.y = (tfxU32)y;

		return u.out;
	}

	inline tfxVec2 UnPack16bit(tfxU32 in) {
		float one_div_32k = 1.f / 32767.f;

		tfxVec2 result;
		result.x = ((signed short)(in & 0x0000FFFF)) * one_div_32k;
		result.y = ((signed short)((in & 0xFFFF0000) >> 16)) * one_div_32k;

		return result;
	}

	inline tfxVec2 UnPack16bitUnsigned(tfxU32 in) {
		float one_div_32k = 1.f / 32767.f;

		tfxVec2 result;
		result.x = ((int)(in & 0x0000FFFF) - 32767) * one_div_32k;
		result.y = ((int)((in & 0xFFFF0000) >> 16) - 32767) * one_div_32k;

		return result;
	}

	inline tfxWideInt PackWide16bit(tfxWideFloat &v_x, tfxWideFloat &v_y) {
		tfxWideFloat w32k = tfxWideSetSingle(32767.f);
		tfxWideInt bits16 = tfxWideSetSinglei(0xFFFF);
		tfxWideInt converted_y = tfxWideConverti(tfxWideMul(v_y, w32k));
		converted_y = tfxWideAndi(converted_y, bits16);
		converted_y = tfxWideShiftLeft(converted_y, 16);
		tfxWideInt converted_x = tfxWideConverti(tfxWideMul(v_x, w32k));
		converted_x = tfxWideAndi(converted_x, bits16);
		return tfxWideOri(converted_x, converted_y);
	}

	inline void UnPackWide16bit(tfxWideInt in, tfxWideFloat &x, tfxWideFloat &y) {
		tfxWideInt w32k = tfxWideSetSinglei(32767);
		x = tfxWideConvert(tfxWideSubi(tfxWideAndi(in, tfxWideSetSinglei(0x0000FFFF)), w32k));
		y = tfxWideConvert(tfxWideSubi(tfxWideShiftRight(tfxWideAndi(in, tfxWideSetSinglei(0xFFFF0000)), 16), w32k));
		x = tfxWideMul(x, one_div_32k_wide);
		y = tfxWideMul(y, one_div_32k_wide);
	}

	inline tfxWideInt PackWide10bit(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z, tfxU32 extra) {
		tfxWideFloat w511 = tfxWideSetSingle(511.f);
		tfxWideInt bits10 = tfxWideSetSinglei(0x3FF);
		tfxWideInt converted_x = tfxWideConverti(tfxWideMul(v_x, w511));
		converted_x = tfxWideAndi(converted_x, bits10);
		converted_x = tfxWideShiftLeft(converted_x, 20);
		tfxWideInt converted_y = tfxWideConverti(tfxWideMul(v_y, w511));
		converted_y = tfxWideAndi(converted_y, bits10);
		converted_y = tfxWideShiftLeft(converted_y, 10);
		tfxWideInt converted_z = tfxWideConverti(tfxWideMul(v_z, w511));
		converted_z = tfxWideAndi(converted_z, bits10);
		tfxWideInt extra_bits = tfxWideShiftLeft(tfxWideSetSinglei(extra), 30);
		return tfxWideOri(tfxWideOri(tfxWideOri(converted_x, converted_y), converted_z), extra_bits);
	}

	inline tfxWideInt PackWide10bitUnsigned(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z) {
		tfxWideFloat w511 = tfxWideSetSingle(511.f);
		tfxWideInt bits10 = tfxWideSetSinglei(0x3FF);
		tfxWideInt converted_x = tfxWideConverti(tfxWideAdd(tfxWideMul(v_x, w511), w511));
		converted_x = tfxWideAndi(converted_x, bits10);
		converted_x = tfxWideShiftLeft(converted_x, 20);
		tfxWideInt converted_y = tfxWideConverti(tfxWideAdd(tfxWideMul(v_y, w511), w511));
		converted_y = tfxWideAndi(converted_y, bits10);
		converted_y = tfxWideShiftLeft(converted_y, 10);
		tfxWideInt converted_z = tfxWideConverti(tfxWideAdd(tfxWideMul(v_z, w511), w511));
		converted_z = tfxWideAndi(converted_z, bits10);
		return tfxWideOri(tfxWideOri(converted_x, converted_y), converted_z);
	}

	inline void UnPackWide10bit(tfxWideInt in, tfxWideFloat &x, tfxWideFloat &y, tfxWideFloat &z) {
		tfxWideInt w511 = tfxWideSetSinglei(511);
		x = tfxWideConvert(tfxWideSubi(tfxWideShiftRight(tfxWideAndi(in, tfxWideSetSinglei(0x3FF00000)), 20), w511));
		y = tfxWideConvert(tfxWideSubi(tfxWideShiftRight(tfxWideAndi(in, tfxWideSetSinglei(0x000FFC00)), 10), w511));
		z = tfxWideConvert(tfxWideSubi(tfxWideAndi(in, tfxWideSetSinglei(0x000003FF)), w511));
		x = tfxWideMul(x, one_div_511_wide);
		y = tfxWideMul(y, one_div_511_wide);
		z = tfxWideMul(z, one_div_511_wide);
	}

	inline tfxWideFloat UnPackWide10bitY(tfxWideInt in) {
		return tfxWideMul(tfxWideConvert(tfxWideSubi(tfxWideShiftRight(tfxWideAndi(in, tfxWideSetSinglei(0x000FFC00)), 10), tfxWideSetSinglei(511))), one_div_511_wide);
	}

	inline tfxWideInt PackWideColor(tfxWideFloat const &v_r, tfxWideFloat const &v_g, tfxWideFloat const &v_b, tfxWideFloat v_a) {
		tfxWideInt color = tfxWideShiftLeft(tfxWideConverti(v_a), 24);
		color = tfxWideAddi(color, tfxWideShiftLeft(tfxWideConverti(v_b), 16));
		color = tfxWideAddi(color, tfxWideShiftLeft(tfxWideConverti(v_g), 8));
		color = tfxWideAddi(color, tfxWideConverti(v_r));
		return color;
	}

	inline tfxWideInt PackWide10bit(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z, tfxWideInt extra) {
		tfxWideFloat w511 = tfxWideSetSingle(511.f);
		tfxWideInt bits10 = tfxWideSetSinglei(0x3FF);
		tfxWideInt converted_x = tfxWideConverti(tfxWideMul(v_x, w511));
		converted_x = tfxWideAndi(converted_x, bits10);
		converted_x = tfxWideShiftLeft(converted_x, 20);
		tfxWideInt converted_y = tfxWideConverti(tfxWideMul(v_y, w511));
		converted_y = tfxWideAndi(converted_y, bits10);
		converted_y = tfxWideShiftLeft(converted_y, 10);
		tfxWideInt converted_z = tfxWideConverti(tfxWideMul(v_z, w511));
		converted_z = tfxWideAndi(converted_z, bits10);
		tfxWideInt extra_bits = tfxWideShiftLeft(extra, 30);
		return tfxWideOri(tfxWideOri(tfxWideOri(converted_x, converted_y), converted_z), extra_bits);
	}

	inline tfxVec4 UnPack10bit(tfxU32 in) {
		tfxUInt10bit unpack;
		unpack.pack = in;
		tfxVec3 result((float)unpack.data.z, (float)unpack.data.y, (float)unpack.data.x);
		result = result * tfxVec3(1.f * one_div_511, 1.f * one_div_511, 1.f * one_div_511);
		return tfxVec4(result, (float)unpack.data.w);
	}

	inline tfxVec3 UnPack10bitVec3(tfxU32 in) {
		tfxUInt10bit unpack;
		unpack.pack = in;
		tfxVec3 result((float)unpack.data.z, (float)unpack.data.y, (float)unpack.data.x);
		return result * tfxVec3(1.f / 511.f, 1.f / 511.f, 1.f / 511.f);
	}

	inline tfxU32 Get2bitFromPacked10bit(tfxU32 in) {
		return ((in >> 30) & 0x3);
	}

	static inline size_t ClampStringSize(size_t compare, size_t string_size) {
		return compare < string_size ? compare : string_size;
	}

	static inline float Distance(float fromx, float fromy, float tox, float toy) {

		float w = tox - fromx;
		float h = toy - fromy;

		return std::sqrt(w * w + h * h);

	}

	inline tfxUInt10bit UintToPacked10bit(tfxU32 in) {
		tfxUInt10bit out;
		out.data.x = (in & 0x3FF);
		out.data.y = ((in >> 10) & 0x3FF);
		out.data.z = ((in >> 20) & 0x3FF);
		out.data.w = ((in >> 30) & 0x3);
		return out;
	}

	inline tfxVec2 InterpolateVec2(float tween, tfxVec2 from, tfxVec2 to) {
		return to * tween + from * (1.f - tween);
	}

	inline tfxVec3 InterpolateVec3(float tween, tfxVec3 from, tfxVec3 to) {
		return to * tween + from * (1.f - tween);
	}

	inline tfxRGBA8 InterpolateRGBA(float tween, tfxRGBA8 from, tfxRGBA8 to) {
		tfxRGBA8 out;
		out.r = char((float)to.r * tween + (float)from.r * (1 - tween));
		out.g = char((float)to.g * tween + (float)from.g * (1 - tween));
		out.b = char((float)to.b * tween + (float)from.b * (1 - tween));
		out.a = char((float)to.a * tween + (float)from.a * (1 - tween));
		return out;
	}

	inline tfxU32 InterpolateAlignment(float tween, tfxU32 from, tfxU32 to) {
		tfxVec3 fromf = UnPack10bitVec3(from);
		tfxVec3 tof = UnPack10bitVec3(to);
		return Pack10bit(InterpolateVec3(tween, fromf, tof), (from >> 30) & 0x3);
	}

	inline tfxVec4 InterpolateVec4(float tween, tfxVec4 &from, tfxVec4 &to) {
		tfx128 l4 = _mm_set_ps1(tween);
		tfx128 l4minus1 = _mm_set_ps1(1.f - tween);
		tfx128 f4 = _mm_set_ps(from.x, from.y, from.z, from.w);
		tfx128 t4 = _mm_set_ps(to.x, to.y, to.z, to.w);
		tfx128 from_lerp = _mm_mul_ps(f4, l4);
		tfx128 to_lerp = _mm_mul_ps(f4, l4minus1);
		tfx128 result = _mm_add_ps(from_lerp, to_lerp);
		tfxVec4 vec;
		_mm_store_ps(&vec.x, result);
		return vec;
	}

	inline tfxWideFloat WideInterpolate(tfxWideFloat tween, tfxWideFloat &from, tfxWideFloat &to) {
		tfxWideFloat one_minus_tween = tfxWideSub(tfxWIDEONE, tween);
		tfxWideFloat to_lerp = tfxWideMul(to, tween);
		tfxWideFloat from_lerp = tfxWideMul(from, one_minus_tween);
		tfxWideFloat result = tfxWideAdd(from_lerp, to_lerp);
		return result;
	}

	inline float Interpolatef(float tween, float from, float to) {
		return to * tween + from * (1.f - tween);
	}

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

	const tfx128 tfxF3_4 = _mm_set_ps1(1.0f / 3.0f);
	const tfx128 tfxF2_4 = _mm_set_ps1(.366025403f);
	const tfx128 tfxG2_4 = _mm_set_ps1(0.211324865f);
	const tfx128 tfxG2_4x2 = _mm_set_ps1(0.42264973f);
	const tfx128 tfxG3_4 = _mm_set_ps1(1.0f / 6.0f);
	const tfx128 tfxG32_4 = _mm_set_ps1((1.0f / 6.0f) * 2.f);
	const tfx128 tfxG33_4 = _mm_set_ps1((1.0f / 6.0f) * 3.f);
	const tfx128i tfxONE = _mm_set1_epi32(1);
	const tfx128 tfxONEF = _mm_set1_ps(1.f);
	const tfx128 tfxZERO = _mm_set1_ps(0.f);
	const tfx128 tfxTHIRTYTWO = _mm_set1_ps(32.f);
	const tfx128i tfxFF = _mm_set1_epi32(0xFF);
	const tfx128 tfxPSIX = _mm_set_ps1(0.6f);

	static inline float Dot(float x1, float y1, float z1, float x2, float y2, float z2)
	{
		return x1 * x2 + y1 * y2 + z1 * z2;
	}

	static inline float Dot(float x1, float y1, float x2, float y2)
	{
		return x1 * x2 + y1 * y2;
	}

	static inline tfx128 Dot128XYZ(const tfx128 &x1, const tfx128 &y1, const tfx128 &z1, const tfx128 &x2, const tfx128 &y2, const tfx128 &z2)
	{
		tfx128 xx = _mm_mul_ps(x1, x2);
		tfx128 yy = _mm_mul_ps(y1, y2);
		tfx128 zz = _mm_mul_ps(z1, z2);
		return _mm_add_ps(xx, _mm_add_ps(yy, zz));
	}

	static inline tfx128 Dot128XY(const tfx128 &x1, const tfx128 &y1, const tfx128 &x2, const tfx128 &y2)
	{
		tfx128 xx = _mm_mul_ps(x1, x2);
		tfx128 yy = _mm_mul_ps(y1, y2);
		return _mm_add_ps(xx, yy);
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
 * Permutation table. This is just a random jumble of all numbers 0-255.
 *
 * This produce a repeatable pattern of 256, but Ken Perlin stated
 * that it is not a problem for graphic texture as the noise features disappear
 * at a distance far enough to be able to see a repeatable pattern of 256.
 *
 * This needs to be exactly the same for all instances on all platforms,
 * so it's easiest to just keep it as static explicit data.
 * This also removes the need for any initialisation of this class.
 *
 * Note that making this an uint32_t[] instead of a uint8_t[] might make the
 * code run faster on platforms with a high penalty for unaligned single
 * byte addressing. Intel x86 is generally single-byte-friendly, but
 * some other CPUs are faster with 4-aligned reads.
 * However, a char[] is smaller, which avoids cache trashing, and that
 * is probably the most important aspect on most architectures.
 * This array is accessed a *lot* by the noise functions.
 * A vector-valued noise over 3D accesses it 96 times, and a
 * float-valued 4D noise 64 times. We want this to fit in the cache!
 */
	const uint8_t perm[] =
	{ 
		151,160,137,91,90,15,
		131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
		190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
		88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
		77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
		102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
		135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
		5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
		223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
		129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
		251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
		49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
		138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
		151,160,137,91,90,15,
		131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
		190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
		88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
		77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
		102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
		135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
		5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
		223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
		129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
		251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
		49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
		138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
	};

	static const uint8_t permMOD12[] =
	{
		7, 4, 5, 7, 6, 3, 11, 1, 9, 11, 0, 5, 2, 5, 7, 9, 8, 0, 7, 6, 9, 10, 8, 3,
		1, 0, 9, 10, 11, 10, 6, 4, 7, 0, 6, 3, 0, 2, 5, 2, 10, 0, 3, 11, 9, 11, 11,
		8, 9, 9, 9, 4, 9, 5, 8, 3, 6, 8, 5, 4, 3, 0, 8, 7, 2, 9, 11, 2, 7, 0, 3, 10,
		5, 2, 2, 3, 11, 3, 1, 2, 0, 7, 1, 2, 4, 9, 8, 5, 7, 10, 5, 4, 4, 6, 11, 6,
		5, 1, 3, 5, 1, 0, 8, 1, 5, 4, 0, 7, 4, 5, 6, 1, 8, 4, 3, 10, 8, 8, 3, 2, 8,
		4, 1, 6, 5, 6, 3, 4, 4, 1, 10, 10, 4, 3, 5, 10, 2, 3, 10, 6, 3, 10, 1, 8, 3,
		2, 11, 11, 11, 4, 10, 5, 2, 9, 4, 6, 7, 3, 2, 9, 11, 8, 8, 2, 8, 10, 7, 10, 5,
		9, 5, 11, 11, 7, 4, 9, 9, 10, 3, 1, 7, 2, 0, 2, 7, 5, 8, 4, 10, 5, 4, 8, 2, 6,
		1, 0, 11, 10, 2, 1, 10, 6, 0, 0, 11, 11, 6, 1, 9, 3, 1, 7, 9, 2, 11, 11, 1, 0,
		10, 7, 1, 7, 10, 1, 4, 0, 0, 8, 7, 1, 2, 9, 7, 4, 6, 2, 6, 8, 1, 9, 6, 6, 7, 5,
		0, 0, 3, 9, 8, 3, 6, 6, 11, 1, 0, 0,
		7, 4, 5, 7, 6, 3, 11, 1, 9, 11, 0, 5, 2, 5, 7, 9, 8, 0, 7, 6, 9, 10, 8, 3,
		1, 0, 9, 10, 11, 10, 6, 4, 7, 0, 6, 3, 0, 2, 5, 2, 10, 0, 3, 11, 9, 11, 11,
		8, 9, 9, 9, 4, 9, 5, 8, 3, 6, 8, 5, 4, 3, 0, 8, 7, 2, 9, 11, 2, 7, 0, 3, 10,
		5, 2, 2, 3, 11, 3, 1, 2, 0, 7, 1, 2, 4, 9, 8, 5, 7, 10, 5, 4, 4, 6, 11, 6,
		5, 1, 3, 5, 1, 0, 8, 1, 5, 4, 0, 7, 4, 5, 6, 1, 8, 4, 3, 10, 8, 8, 3, 2, 8,
		4, 1, 6, 5, 6, 3, 4, 4, 1, 10, 10, 4, 3, 5, 10, 2, 3, 10, 6, 3, 10, 1, 8, 3,
		2, 11, 11, 11, 4, 10, 5, 2, 9, 4, 6, 7, 3, 2, 9, 11, 8, 8, 2, 8, 10, 7, 10, 5,
		9, 5, 11, 11, 7, 4, 9, 9, 10, 3, 1, 7, 2, 0, 2, 7, 5, 8, 4, 10, 5, 4, 8, 2, 6,
		1, 0, 11, 10, 2, 1, 10, 6, 0, 0, 11, 11, 6, 1, 9, 3, 1, 7, 9, 2, 11, 11, 1, 0,
		10, 7, 1, 7, 10, 1, 4, 0, 0, 8, 7, 1, 2, 9, 7, 4, 6, 2, 6, 8, 1, 9, 6, 6, 7, 5,
		0, 0, 3, 9, 8, 3, 6, 6, 11, 1, 0, 0
	};

	// 2D Perlin simplex noise (here for reference only)
	float tfxNoise(float x, float y);
	// 4 noise samples using simd
	tfx128Array tfxNoise4(const tfx128 x4, const tfx128 y4);
	tfx128Array tfxNoise4(const tfx128 &x4, const tfx128 &y4, const tfx128 &z4);

	/*
		Start of xxHash code that encompasses the following license
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
		// create new XXHash (64 bit)
		/** @param seed your seed value, even zero is a valid seed **/
		explicit tfxXXHash64(uint64_t seed)
		{
			state[0] = seed + Prime1 + Prime2;
			state[1] = seed + Prime2;
			state[2] = seed;
			state[3] = seed - Prime1;
			bufferSize = 0;
			totalLength = 0;
			memset(buffer, 0, MaxBufferSize);
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
		inline void         Clear() { current_size = 0; }
		inline char*           begin() { return strbuffer(); }
		inline const char*     begin() const { return strbuffer(); }
		inline char*           end() { return strbuffer() + current_size; }
		inline const char*     end() const { return strbuffer() + current_size; }
		inline char&           back() { assert(current_size > 0); return strbuffer()[current_size - 1]; }
		inline const char&     back() const { assert(current_size > 0); return strbuffer()[current_size - 1]; }
		inline void         pop() { assert(current_size > 0); current_size--; }
		inline void	        push_back(const char v) { if (current_size == capacity) reserve(_grow_capacity(current_size + 1)); new((void*)(data + current_size)) char(v); current_size++; }

		inline tfxU32       _grow_capacity(tfxU32 sz) const { tfxU32 new_capacity = capacity ? (capacity + capacity / 2) : 8; return new_capacity > sz ? new_capacity : sz; }
		inline void         resize(tfxU32 new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
		inline void         reserve(tfxU32 new_capacity) {
			if (new_capacity <= capacity) return;
			char* new_data = (char*)tfxALLOCATE(0, 0, (size_t)new_capacity * sizeof(char));
			assert(new_data);	//unable to allocate memory. Todo: proper handling
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

		tfxStr(const char *text) : data(NULL), current_size(0), capacity(0), is_local_buffer(false) { size_t length = strnlen_s(text, 512); if (!length) { Clear(); return; }; if (capacity < length) reserve((tfxU32)length); assert(data); memcpy(data, text, length); current_size = (tfxU32)length; NullTerminate(); }
		tfxStr(const tfxStr &src) : data(NULL), current_size(0), capacity(0), is_local_buffer(false) { size_t length = src.Length(); if (!length) { Clear(); return; }; if (capacity < length) reserve((tfxU32)length); assert(data); memcpy(data, src.data, length); current_size = (tfxU32)length; NullTerminate(); }
		inline void operator=(const char *text) { size_t length = strnlen_s(text, 512); if (!length) { Clear(); return; }; if (capacity < length) reserve((tfxU32)length); assert(data); memcpy(data, text, length); current_size = (tfxU32)length; NullTerminate(); }
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
			if (back() == NULL)
				pop();
			pop();
			NullTerminate();
		}
		inline void Trim(char c = ' ') {
			if (!Length()) return;
			if (back() == NULL)
				pop();
			while (back() == c && current_size) {
				pop();
			}
			NullTerminate();
		}
		inline void TrimFront(char c = ' ') {
			if (!Length()) return;
			tfxU32 pos = 0;
			while (strbuffer()[pos] == c && pos < current_size) {
				pos++;
			}
			if (pos < current_size) {
				memcpy(strbuffer(), strbuffer() + pos, current_size - pos);
			}
			current_size -= pos;
		}
		void NullTerminate() { push_back(NULL); }
		bool SaveToFile(const char *file_name);
		const bool IsInt() const;
		const bool IsFloat() const;
	};

#define tfxStrType(type, size)		\
	struct type : public tfxStr { \
	char buffer[size]; \
	type() { memset(buffer, 0, size); data = buffer; capacity = size; current_size = 0; is_local_buffer = true; NullTerminate(); } \
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
	type(const char *text) { memset(buffer, 0, size); data = buffer; is_local_buffer = true; capacity = size; size_t length = strnlen_s(text, size); if (!length) { Clear(); return; } memcpy(data, text, length); current_size = (tfxU32)length; NullTerminate(); } \
	type(const tfxStr &src) { \
		memset(buffer, 0, size); \
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
		memset(buffer, 0, size); \
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

	/*
	//Unwrapped for convenience when debugging. Can be removed at some point
	struct tfxStr64 : public tfxStr {
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
	//Not sure how efficient a hash lookup with this is, could probably be better, but isn't used much at all in any realtime particle updating.
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


	struct tfxFace {
		int v[3];
	};
	extern tfxvec<tfxVec3> tfxIcospherePoints[6];
	void MakeIcospheres();
	int VertexForEdge(tfxStorageMap<int> &point_cache, tfxvec<tfxVec3>& vertices, int first, int second);
	tfxvec<tfxFace> SubDivideIcosphere(tfxStorageMap<int> &point_cache, tfxvec<tfxVec3>& vertices, tfxvec<tfxFace> &triangles);

	static inline tfxU32 Millisecs() {
		auto now = tfxClock::now().time_since_epoch();
		auto m = std::chrono::duration_cast<std::chrono::milliseconds>(now).count();
		return (tfxU32)m;
	}
	static inline uint64_t Microsecs() {
		auto now = tfxClock::now().time_since_epoch();
		auto m = std::chrono::duration_cast<std::chrono::microseconds>(now).count();
		return m;
	}

	struct tfxProfileStats {
		tfxU64 cycle_high;
		tfxU64 cycle_low;
		tfxU64 cycle_average;
		tfxU64 time_high;
		tfxU64 time_low;
		tfxU64 time_average;
		tfxU32 hit_count;
	};

	struct tfxProfileSnapshot {
		tfxU32 hit_count;
		tfxU64 run_time;
		tfxU64 cycle_count;
	};

	struct tfxProfile {
		const char *name;
		tfxProfileSnapshot snapshots[tfxPROFILER_SAMPLES];
	};

	extern const tfxU32 tfxPROFILE_COUNT;
	extern tfxU32 tfxCurrentSnapshot;
	tfxProfile tfxProfileArray[];

	struct tfxProfileTag {
		tfxProfile *profile;
		tfxProfileSnapshot *snapshot;
		tfxU64 start_cycles;
		tfxU64 start_time;

		tfxProfileTag(tfxU32 id, const char *name);

		~tfxProfileTag() {
			AtomicAdd64(&snapshot->run_time, Microsecs() - start_time);
			AtomicAdd64(&snapshot->cycle_count, (__rdtsc() - start_cycles));
		}

	};

	inline void GatherStats(tfxProfile *profile, tfxProfileStats *stat) {
		stat->cycle_high = 0;
		stat->cycle_low = tfxMAX_64u;
		stat->time_high = 0;
		stat->time_low = tfxMAX_64u;
		stat->hit_count = 0;
		stat->cycle_average = 0;
		stat->time_average = 0;
		for (int i = 0; i != tfxPROFILER_SAMPLES; ++i) {
			tfxProfileSnapshot *snap = profile->snapshots + i;
			stat->cycle_high = tfxMax(snap->cycle_count, stat->cycle_high);
			stat->cycle_low = tfxMin(snap->cycle_count, stat->cycle_low);
			stat->cycle_average += snap->cycle_count;
			stat->time_high = tfxMax(snap->run_time, stat->time_high);
			stat->time_low = tfxMin(snap->run_time, stat->time_low);
			stat->time_average += snap->run_time;
			stat->hit_count += snap->hit_count;
		}
		stat->cycle_average /= tfxPROFILER_SAMPLES;
		stat->time_average /= tfxPROFILER_SAMPLES;
		stat->hit_count /= tfxPROFILER_SAMPLES;
	}

	inline void ResetSnapshot(tfxProfileSnapshot *snapshot) {
		snapshot->cycle_count = 0;
		snapshot->hit_count = 0;
		snapshot->run_time = 0;
	}

#ifdef tfxENABLE_PROFILING 
#define tfxPROFILE tfxProfileTag tfx_tag((tfxU32)__COUNTER__, __FUNCTION__)
#else
#define tfxPROFILE __COUNTER__
#endif

	const tfxU32 tfxMAGIC_NUMBER = '!XFT';
	const tfxU32 tfxMAGIC_NUMBER_INVENTORY = '!VNI';
	const tfxU32 tfxFILE_VERSION = 1;	//Not doing anything with this yet

	//Basic package manager used for reading/writing effects files
	struct tfxHeader {
		tfxU32 magic_number;						//Magic number to confirm file format
		tfxU32 file_version;						//The version of the file
		tfxU32 flags;								//Any state_flags for the file
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
		tfxU64 offset_from_start_of_file = 0;		//Offset from the start of the file to where the file is located
		tfxU64 file_size = 0;						//The size of the file
		tfxStream data;								//The file data

		void FreeData();
	};

	struct tfxInventory {
		tfxU32 magic_number;						//Magic number to confirm format of the Inventory
		tfxU32 entry_count;							//Number of files in the inventory
		tfxStorageMap<tfxEntryInfo> entries;		//The inventory list

		tfxInventory() :
			entries("Inventory Map", "Inventory Data"),
			magic_number(0),
			entry_count(0)
		{}
	};

	struct tfxPackage {
		tfxStr file_path;
		tfxHeader header;
		tfxInventory inventory;
		tfxU64 file_size = 0;						//The total file size of the package, should match file size on disk
		tfxStream file_data;						//Dump of the data from the package file on disk

		~tfxPackage();

		tfxEntryInfo *GetFile(const char *name);
		void AddFile(tfxEntryInfo file);
		void AddFile(const char *file_name, tfxStream &data);
		void Free();

	};

	tfxStream ReadEntireFile(const char *file_name, bool terminate = false);
	tfxErrorFlags LoadPackage(const char *file_name, tfxPackage &package);
	tfxErrorFlags LoadPackage(tfxStream &stream, tfxPackage &package);
	tfxPackage CreatePackage(const char *file_path);
	bool SavePackageDisk(tfxPackage &package);
	tfxStream SavePackageMemory(tfxPackage &package);
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

		tfxAttributeNode() : frame(0.f), value(0.f), flags(0), index(0) { }
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
			flags = ~tfxAttributeNodeFlags_is_curve;
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
		tfxU64 seeds[2];

		tfxRandom() {
			memset(seeds, 0, sizeof(tfxU64) * 2);
			ReSeed();
		}

		tfxRandom(tfxU32 seed) {
			memset(seeds, 0, sizeof(tfxU64) * 2);
			ReSeed(seed);
		}

		void ReSeed() {
			seeds[0] = Millisecs(); seeds[1] = (tfxU64)Millisecs() * 2; 
			Advance();
		}

		void ReSeed(tfxU64 seed1, tfxU64 seed2) {
			seeds[0] = seed1;
			seeds[1] = seed2;
			Advance();
		}

		void ReSeed(tfxU64 seed) {
			seeds[0] = seed;
			seeds[1] = seed * 2;
			Advance();
		}

		inline void Advance() {
			tfxU64 s1 = seeds[0];
			tfxU64 s0 = seeds[1];
			seeds[0] = s0;
			s1 ^= s1 << 23; // a
			seeds[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5); // b, c
		}

		inline float Generate() {
			tfxU64 s1 = seeds[0];
			tfxU64 s0 = seeds[1];
			tfxU64 result = s0 + s1;
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

		inline int Range(int from, int to) {
			float a = (to - from) * Generate() - (to - from) * .5f;
			return a < 0 ? int(a - 0.5f) : int(a + 0.5f);
		};

		inline tfxU32 Range(tfxU32 max) {
			float g = Generate();
			float a = g * (float)max;
			return tfxU32(a);
		};

		inline void AlterSeed(tfxU64 amount) {
			seeds[0] *= amount;
			seeds[1] += amount;
		};

	};

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
		float GetRandomValue(float age, tfxRandom &seed);
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
		void DeleteNodeAtFrame(float frame);
		void Reset(float first_node_value, tfxGraphPreset preset, bool add_node = true);
		void DragValues(tfxGraphPreset preset, float &frame, float &value);
		void Clear();
		void Free();
		void Copy(tfxGraph &to, bool compile = true);
		bool Sort();
		void ReIndex();
		tfxVec2 GetInitialZoom();
		tfxVec2 GetInitialZoom3d();
		bool IsOvertimeGraph();
		bool IsGlobalGraph();
		bool IsAngleGraph();
		bool IsTranslationGraph();
		void MultiplyAllValues(float scalar);
		void CopyToNoLookups(tfxGraph *graph);

	};

	tfxVec4 GetMinMaxGraphValues(tfxGraphPreset preset);

	//todo:: Inline a lot of these.
	tfxVec2 GetQuadBezier(tfxVec2 p0, tfxVec2 p1, tfxVec2 p2, float t, float ymin, float ymax, bool clamp = true);
	tfxVec2 GetCubicBezier(tfxVec2 p0, tfxVec2 p1, tfxVec2 p2, tfxVec2 p3, float t, float ymin, float ymax, bool clamp = true);
	float GetBezierValue(const tfxAttributeNode *lastec, const tfxAttributeNode &a, float t, float ymin, float ymax);
	float GetDistance(float fromx, float fromy, float tox, float toy);
	inline float GetVectorAngle(float x, float y) {
		return atan2f(x, -y);
	}
	static bool CompareNodes(tfxAttributeNode &left, tfxAttributeNode &right);
	void CompileGraph(tfxGraph &graph);
	void CompileGraphOvertime(tfxGraph &graph);
	float GetMaxLife(tfxEffectEmitter &e);
	float LookupFastOvertime(tfxGraph &graph, float age, float lifetime);
	float LookupFast(tfxGraph &graph, float frame);
	float LookupPreciseOvertime(tfxGraph &graph, float age, float lifetime);
	float LookupPrecise(tfxGraph &graph, float frame);
	float GetRandomFast(tfxGraph &graph, float frame, tfxRandom &seed);
	float GetRandomPrecise(tfxGraph &graph, float frame, tfxRandom &seed);

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
	bool IsEmitterGraph(tfxGraphType type);
	bool IsTransformGraph(tfxGraphType type);
	bool IsGlobalPercentageGraph(tfxGraphType type);
	bool IsAngleGraph(tfxGraphType type);
	bool IsAngleOvertimeGraph(tfxGraphType type);
	bool IsEverythingElseGraph(tfxGraphType type);
	inline bool HasNodeAtFrame(tfxGraph &graph, float frame) {
		graph.nodes.ResetIteratorIndex();
		do {
			for (auto &node : graph.nodes) {
				if (node.frame == frame) return true;
			}
		} while (!graph.nodes.EndOfBuckets());
		return false;
	}
	bool HasKeyframes(tfxEffectEmitter &e);
	bool HasMoreThanOneKeyframe(tfxEffectEmitter &e);
	void PushTranslationPoints(tfxEffectEmitter &e, tfxStack<tfxVec3> &points, float frame);

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
			emitter_width.CopyToNoLookups(&dst->emitter_width);
			emitter_height.CopyToNoLookups(&dst->emitter_height);
			emitter_depth.CopyToNoLookups(&dst->emitter_depth);
		}

	};

	struct tfxTransformAttributes {
		tfxGraph roll;
		tfxGraph pitch;
		tfxGraph yaw;
		tfxGraph translation_x;
		tfxGraph translation_y;
		tfxGraph translation_z;

		void Initialise(tfxMemoryArenaManager *allocator, tfxMemoryArenaManager *value_allocator, tfxU32 bucket_size = 8) {
			roll.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			pitch.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			yaw.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			translation_x.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			translation_y.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			translation_z.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);

			roll.lookup.values.allocator = value_allocator;
			pitch.lookup.values.allocator = value_allocator;
			yaw.lookup.values.allocator = value_allocator;
			translation_x.lookup.values.allocator = value_allocator;
			translation_y.lookup.values.allocator = value_allocator;
			translation_z.lookup.values.allocator = value_allocator;
		}

		void Free() {
			roll.Free();
			pitch.Free();
			yaw.Free();
			translation_x.Free();
			translation_y.Free();
			translation_z.Free();
		}

		void CopyToNoLookups(tfxTransformAttributes *dst) {
			roll.CopyToNoLookups(&dst->roll);
			pitch.CopyToNoLookups(&dst->pitch);
			yaw.CopyToNoLookups(&dst->yaw);
			translation_x.CopyToNoLookups(&dst->translation_x);
			translation_y.CopyToNoLookups(&dst->translation_y);
			translation_z.CopyToNoLookups(&dst->translation_z);
		}
	};

	inline bool HasTranslationKeyframes(tfxTransformAttributes &graphs) {
		return graphs.translation_x.nodes.size() || graphs.translation_y.nodes.size() || graphs.translation_z.nodes.size();
	}

	inline void AddTranslationNodes(tfxTransformAttributes &keyframes, float frame) {
		if (keyframes.translation_x.nodes.size()) {
			if (!HasNodeAtFrame(keyframes.translation_x, frame))
				keyframes.translation_x.AddCoordNode(frame, 0.f);
			if (!HasNodeAtFrame(keyframes.translation_y, frame))
				keyframes.translation_y.AddCoordNode(frame, 0.f);
			if (!HasNodeAtFrame(keyframes.translation_z, frame))
				keyframes.translation_z.AddCoordNode(frame, 0.f);
		}
		else {
			keyframes.translation_x.AddCoordNode(0.f, 0.f);
			keyframes.translation_y.AddCoordNode(0.f, 0.f);
			keyframes.translation_z.AddCoordNode(0.f, 0.f);
			if (frame != 0) {
				keyframes.translation_x.AddCoordNode(frame, 0.f);
				keyframes.translation_y.AddCoordNode(frame, 0.f);
				keyframes.translation_z.AddCoordNode(frame, 0.f);
			}
		}
	}

	struct tfxPropertyAttributes {
		tfxGraph emission_pitch;
		tfxGraph emission_yaw;
		tfxGraph emission_range;
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
			splatter.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			emitter_width.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			emitter_height.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			emitter_depth.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			arc_size.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
			arc_offset.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);

			emission_pitch.lookup.values.allocator = value_allocator;
			emission_yaw.lookup.values.allocator = value_allocator;
			emission_range.lookup.values.allocator = value_allocator;
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

	static float(*lookup_overtime_callback)(tfxGraph &graph, float age, float lifetime);
	static float(*lookup_callback)(tfxGraph &graph, float age);
	static float(*lookup_random_callback)(tfxGraph &graph, float age, tfxRandom &seed);

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

	struct tfxCameraSettings {
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
		tfxCameraSettings camera_settings;
		float effect_z_offset;
		float camera_speed;
		bool attach_effect_to_camera = false;
	};

	//this probably only needs to be in the editor, no use for it in the library? Maybe in the future as an alternative way to play back effects...
	struct tfxSpriteSheetSettings {
		tfxVec3 position;
		tfxVec2 frame_size;
		float scale;
		float zoom;
		int frames;
		int current_frame;
		int frame_offset;
		int extra_frames_count;
		tfxU32 seed;
		tfxAnimationFlags animation_flags;
		tfxU32 needs_exporting;
		float max_radius;
		tfxU32 largest_frame;
		float playback_speed;
		float effect_z_offset = 5.f;
		tfxExportColorOptions color_option;
		tfxExportOptions export_option;
		tfxCameraSettings camera_settings;
		tfxCameraSettings camera_settings_orthographic;
	};

	//This struct has the settings for recording sprite data frames so that they can be played back as an alternative to dynamic particle updating
	struct tfxSpriteDataSettings {
		int frames;
		int current_frame;
		float current_time;
		float animation_time;
		int frame_offset;
		int extra_frames_count;
		tfxU32 seed;
		tfxAnimationFlags animation_flags;
		tfxU32 needs_exporting;
		float max_radius;
		tfxU32 largest_frame;
		float playback_speed;
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
			image_index(0),
			ptr(nullptr),
			animation_frames(1.f),
			shape_index(0),
			max_radius(0),
			import_filter(0),
			compute_shape_index(0)
		{ }
	};

	//Struct of Arrays for the emitter properties. 
	struct tfxEmitterPropertiesSoA {
		//Angle added to the rotation of the particle when spawned or random angle range if angle setting is set to tfxRandom
		tfxVec3 *angle_offsets;
		//When aligning the billboard along a vector, you can set the type of vector that it aligns with
		tfxVectorAlignType *vector_align_type;
		//Point, area, ellipse emitter etc.
		tfxEmissionType *emission_type;
		//If single shot flag is set then you can limit how many times it will loop over it's overtime graphs before expiring
		tfxU32 *single_shot_limit;
		//Animation frame rate
		float *frame_rate;
		//The final frame index of the animation
		float *end_frame;
		//Pointer to the ImageData in the EffectLibary. 
		tfxImageData **image;
		//For 3d effects, the type of billboarding: 0 = use billboarding (always face camera), 1 = No billboarding, 2 = No billboarding and align with motion
		tfxBillboardingOptions *billboard_option;

		//The number of rows/columns/ellipse/line points in the grid when spawn on grid flag is used
		tfxVec3 *grid_points;
		//The rotation of particles when they spawn, or behave overtime if tfxAlign is used
		tfxAngleSettingFlags *angle_settings;
		//Layer of the particle manager that the particle is added to
		tfxU32 *layer;
		//Milliseconds to delay spawing
		float *delay_spawning;
		//Should particles emit towards the center of the emitter or away, or in a specific direction
		tfxEmissionDirection *emission_direction;

		//How particles should behave when they reach the end of the line
		tfxLineTraversalEndBehaviour *end_behaviour;
		//Bit field of various boolean state_flags
		tfxParticleControlFlags *compute_flags;
		//Offset to draw particles at
		tfxVec2 *image_handle;
		//Offset of emitters
		tfxVec3 *emitter_handle;
		//When single flag is set, spawn this amount of particles in one go
		tfxU32 *spawn_amount;
		//The shape being used for all particles spawned from the emitter
		tfxU32 *shape_index;
		//The number of millisecs before an effect or emitter will loop back round to the beginning of it's graph lookups
		float *loop_length;
		//The start frame index of the animation
		float *start_frame;
		//Base noise offset random range so that noise patterns don't repeat so much over multiple effects
		float *noise_base_offset_range;

		tfxEmitterPropertiesSoA() { memset(this, 0, sizeof(tfxEmitterPropertiesSoA)); }
	};

	inline void InitEmitterProperites(tfxEmitterPropertiesSoA &properties, tfxU32 i) {
		properties.angle_offsets[i] = { 0.f, 0.f, tfx360Radians };
		properties.image[i] = nullptr;
		properties.image_handle[i] = tfxVec2();
		properties.spawn_amount[i] = 1;
		properties.single_shot_limit[i] = 0;
		properties.emission_type[i] = tfxEmissionType::tfxPoint;
		properties.billboard_option[i] = tfxBillboarding;
		properties.vector_align_type[i] = tfxVectorAlignType_motion;
		properties.emission_direction[i] = tfxEmissionDirection::tfxOutwards;
		properties.grid_points[i] = { 10.f, 10.f, 10.f };
		properties.emitter_handle[i] = { 0.f, 0.f, 0.f };
		properties.end_behaviour[i] = tfxLineTraversalEndBehaviour::tfxLoop;
		properties.loop_length[i] = 0.f;
		properties.layer[i] = 0;
		properties.shape_index[i] = 1;
		properties.start_frame[i] = 0;
		properties.end_frame[i] = 0;
		properties.frame_rate[i] = 30.f;
		properties.angle_settings[i] = tfxAngleSettingFlags_random_roll | tfxAngleSettingFlags_specify_pitch | tfxAngleSettingFlags_specify_yaw;
		properties.delay_spawning[i] = 0.f;
		properties.noise_base_offset_range[i] = 1000.f;
	}

	//Use with care, no checks for out of bounds
	inline void CopyEmitterProperites(tfxEmitterPropertiesSoA &from_properties, tfxU32 from_i, tfxEmitterPropertiesSoA &to_properties, tfxU32 to_i) {
		to_properties.angle_offsets[to_i] = from_properties.angle_offsets[from_i];
		to_properties.image[to_i] = from_properties.image[from_i];
		to_properties.image_handle[to_i] = from_properties.image_handle[from_i];
		to_properties.spawn_amount[to_i] = from_properties.spawn_amount[from_i];
		to_properties.single_shot_limit[to_i] = from_properties.single_shot_limit[from_i];
		to_properties.emission_type[to_i] = from_properties.emission_type[from_i];
		to_properties.billboard_option[to_i] = from_properties.billboard_option[from_i];
		to_properties.vector_align_type[to_i] = from_properties.vector_align_type[from_i];
		to_properties.emission_direction[to_i] = from_properties.emission_direction[from_i];
		to_properties.grid_points[to_i] = from_properties.grid_points[from_i];
		to_properties.emitter_handle[to_i] = from_properties.emitter_handle[from_i];
		to_properties.end_behaviour[to_i] = from_properties.end_behaviour[from_i];
		to_properties.loop_length[to_i] = from_properties.loop_length[from_i];
		to_properties.layer[to_i] = from_properties.layer[from_i];
		to_properties.shape_index[to_i] = from_properties.shape_index[from_i];
		to_properties.start_frame[to_i] = from_properties.start_frame[from_i];
		to_properties.end_frame[to_i] = from_properties.end_frame[from_i];
		to_properties.frame_rate[to_i] = from_properties.frame_rate[from_i];
		to_properties.angle_settings[to_i] = from_properties.angle_settings[from_i];
		to_properties.delay_spawning[to_i] = from_properties.delay_spawning[from_i];
		to_properties.noise_base_offset_range[to_i] = from_properties.noise_base_offset_range[from_i];
	}

	struct tfxEmitterTransform {
		//Position, scale and rotation values
		tfxVec3 translation;
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


	float GetEmissionDirection2d(tfxParticleManager &pm, tfxLibrary *library, tfxRandom &random, tfxU32 property_index, tfxU32 index, tfxVec2 local_position, tfxVec2 world_position, tfxVec2 emitter_size);
	tfxVec3 GetEmissionDirection3d(tfxParticleManager &pm, tfxLibrary *library, tfxRandom &random, tfxU32 property_index, tfxU32 index, float emission_pitch, float emission_yaw, tfxVec3 local_position, tfxVec3 world_position, tfxVec3 emitter_size);

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
		//Index to sprite sheet settings stored in the effect library. 
		tfxU32 sprite_sheet_settings_index;
		//Index to sprite data settings stored in the effect library. 
		tfxU32 sprite_data_settings_index;
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
			lookup_node_index(0),
			lookup_value_index(0),
			sprite_data_settings_index(0),
			uid(0),
			max_radius(0),
			sprite_sheet_settings_index(0),
			preview_camera_settings(0),
			max_sub_emitters(0),
			max_life(0),
			max_sub_emitter_life(0.f),
			sub_effectors(tfxCONSTRUCTOR_VEC_INIT("sub_effectors"))
		{
			for (int i = 0; i != tfxLAYERS; ++i) {
				max_particles[i] = 0;
			}
		}
	};

#define tfx2DTRANSFORMCALLBACK(name) void *name(const tfxVec2 local_position, const float roll, tfxVec2 &world_position, float &world_rotations, const tfxVec3 &parent_rotations, const tfxMatrix4 &matrix, const tfxVec3 &handle, const tfxVec3 &scale, const tfxVec3 &from_position)
	typedef tfx2DTRANSFORMCALLBACK(tfxParticleTransformCallback2d);

	struct tfxEmitterSoA {
		void(**transform_particle_callback2d)(const float local_position_x, const float local_position_y, const float roll, tfxVec2 &world_position, float &world_rotations, const tfxVec3 &parent_rotations, const tfxMatrix4 &matrix, const tfxVec3 &handle, const tfxVec3 &scale, const tfxVec3 &from_position);

		//State data
		float *frame;
		float *age;
		float *highest_particle_age;
		float *delay_spawning;
		float *timeout_counter;
		float *timeout;
		tfxVec3 *handle;
		tfxEmitterPropertyFlags *property_flags;
		float *loop_length;
		//Position, scale and rotation values
		tfxVec3 *translation;
		tfxVec3 *local_position;
		tfxVec3 *world_position;
		tfxVec3 *captured_position;
		tfxVec3 *local_rotations;
		tfxVec3 *world_rotations;
		tfxVec3 *scale;
		//Todo: save space and use a quaternion here... maybe
		tfxMatrix4 *matrix;
		tfxVec2 *image_handle;
		float *amount_remainder;
		float *spawn_quantity;
		float *qty_step_size;

		tfxU32 *emitter_attributes;
		tfxU32 *transform_attributes;
		tfxU32 *overtime_attributes;

		tfxU32 *parent_index;
		tfxU32 *properties_index;
		tfxU32 *info_index;
		tfxU32 *hierarchy_depth;
		tfxU32 *sprites_count;
		tfxU32 *sprites_index;
		tfxKey *path_hash;

		//Spawn controls
		float *life;
		float *life_variation;
		float *arc_size;
		float *arc_offset;
		float *weight;
		float *weight_variation;
		float *velocity;
		float *velocity_variation;
		float *spin;
		float *spin_variation;
		float *splatter;
		float *noise_offset_variation;
		float *noise_offset;
		float *noise_resolution;
		tfxVec2 *size;
		tfxVec2 *size_variation;
		tfxVec3 *grid_segment_size;

		//Control Data
		tfxU32 *particles_index;
		float *overal_scale;
		float *velocity_adjuster;
		float *intensity;
		float *image_frame_rate;
		float *stretch;
		float *end_frame;
		tfxVec3 *grid_coords;
		tfxVec3 *grid_direction;
		tfxVec3 *emitter_size;
		float *emission_alternator;
		tfxEmitterStateFlags *state_flags;
		tfxVec2 *image_size;
		tfxVec3 *angle_offsets;

	};

	inline void InitEmitterSoA(tfxSoABuffer *buffer, tfxEmitterSoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(void*), offsetof(tfxEmitterSoA, transform_particle_callback2d));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, frame));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, age));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, highest_particle_age));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, delay_spawning));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, timeout_counter));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, timeout));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, handle));
		AddStructArray(buffer, sizeof(tfxEmitterPropertyFlags), offsetof(tfxEmitterSoA, property_flags));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, loop_length));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, translation));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, local_position));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, world_position));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, captured_position));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, local_rotations));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, world_rotations));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, scale));

		//Todo: save space and use a quaternion here?
		AddStructArray(buffer, sizeof(tfxMatrix4), offsetof(tfxEmitterSoA, matrix));
		AddStructArray(buffer, sizeof(tfxVec2), offsetof(tfxEmitterSoA, image_handle));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, amount_remainder));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, spawn_quantity));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, qty_step_size));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, emitter_attributes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, transform_attributes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, overtime_attributes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, parent_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, sprites_count));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, sprites_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, properties_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, info_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, hierarchy_depth));
		AddStructArray(buffer, sizeof(tfxKey), offsetof(tfxEmitterSoA, path_hash));

		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, life));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, life_variation));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, arc_size));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, arc_offset));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, weight));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, weight_variation));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, velocity));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, velocity_variation));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, spin));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, spin_variation));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, splatter));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, noise_offset_variation));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, noise_offset));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, noise_resolution));
		AddStructArray(buffer, sizeof(tfxVec2), offsetof(tfxEmitterSoA, size));
		AddStructArray(buffer, sizeof(tfxVec2), offsetof(tfxEmitterSoA, size_variation));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, grid_segment_size));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, particles_index));

		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, overal_scale));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, velocity_adjuster));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, intensity));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, image_frame_rate));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, stretch));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, end_frame));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, grid_coords));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, grid_direction));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, emitter_size));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, emission_alternator));
		AddStructArray(buffer, sizeof(tfxEmitterStateFlags), offsetof(tfxEmitterSoA, state_flags));
		AddStructArray(buffer, sizeof(tfxVec2), offsetof(tfxEmitterSoA, image_size));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, angle_offsets));

		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	struct tfxEffectSoA {
		//State data
		float *frame;
		float *age;
		float *highest_particle_age;
		float *timeout_counter;
		float *timeout;
		tfxVec3 *handle;
		tfxEmitterPropertyFlags *property_flags;
		float *loop_length;
		//Position, scale and rotation values
		tfxVec3 *translation;
		tfxVec3 *local_position;
		tfxVec3 *world_position;
		tfxVec3 *captured_position;
		tfxVec3 *local_rotations;
		tfxVec3 *world_rotations;
		tfxVec3 *scale;
		//Todo: save space and use a quaternion here?
		tfxMatrix4 *matrix;
		tfxU32 *global_attributes;
		tfxU32 *transform_attributes;

		tfxU32 *properties_index;
		tfxU32 *info_index;
		tfxU32 *parent_particle_index;
		tfxLibrary **library;

		//Spawn controls
		tfxParentSpawnControls *spawn_controls;
		tfxVec3 *emitter_size;
		float *stretch;
		float *overal_scale;
		float *noise_base_offset;
		tfxEmitterStateFlags *state_flags;

		//User Data
		void **user_data;
		void(**update_callback)(tfxParticleManager *pm, tfxEffectID effect_index);
	};

	inline void InitEffectSoA(tfxSoABuffer *buffer, tfxEffectSoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, frame));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, age));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, highest_particle_age));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, timeout_counter));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, timeout));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, handle));
		AddStructArray(buffer, sizeof(tfxEmitterPropertyFlags), offsetof(tfxEffectSoA, property_flags));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, loop_length));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, translation));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, local_position));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, world_position));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, captured_position));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, local_rotations));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, world_rotations));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, scale));

		//Todo: save space and use a quaternion here?
		AddStructArray(buffer, sizeof(tfxMatrix4), offsetof(tfxEffectSoA, matrix));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEffectSoA, global_attributes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEffectSoA, transform_attributes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEffectSoA, parent_particle_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEffectSoA, properties_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEffectSoA, info_index));
		AddStructArray(buffer, sizeof(void*), offsetof(tfxEffectSoA, library));
		AddStructArray(buffer, sizeof(tfxParentSpawnControls), offsetof(tfxEffectSoA, spawn_controls));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, emitter_size));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, stretch));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, overal_scale));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, noise_base_offset));
		AddStructArray(buffer, sizeof(tfxEffectStateFlags), offsetof(tfxEffectSoA, state_flags));

		AddStructArray(buffer, sizeof(void*), offsetof(tfxEffectSoA, user_data));
		AddStructArray(buffer, sizeof(void*), offsetof(tfxEffectSoA, update_callback));

		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	//An tfxEffectEmitter can either be an effect which stores emitters and global graphs for affecting all the attributes in the emitters
	//Or it can be an emitter which spawns all of the particles. Effectors are stored in the particle manager effects list buffer.
	//This is only for library storage, when using to update each frame this is copied to tfxEffectSoA and tfxEmitterSoA for realtime updates
	//suited for realtime use.
	struct tfxEffectEmitter {
		//Required for frame by frame updating
		//The current state of the effect/emitter used in the editor only at this point
		tfxEmitterStateFlags state_flags;
		tfxEmitterPropertyFlags property_flags;
		tfxLibrary *library;
		//Is this an tfxEffectType or tfxEmitterType
		tfxEffectEmitterType type;
		//The index within the library that this exists at
		tfxU32 library_index;
		//A hash of the directory path to the effect ie Flare/spark, and also a UID for the effect/emitter
		tfxKey path_hash;
		//All graphs that the effect uses to lookup attribute values are stored in the library. These variables here are indexes to the array where they're stored
		tfxU32 global;
		tfxU32 emitter_attributes;
		tfxU32 transform_attributes;
		//Pointer to the immediate parent
		tfxEffectEmitter *parent;
		//State state_flags for emitters and effects
		tfxEffectPropertyFlags effect_flags;
		//When not using insert sort to guarantee particle order, sort passes offers a more lax way of ordering particles over a number of frames.
		//The more passes the more quickly ordered the particles will be but at a higher cost
		tfxU32 sort_passes;
		//Custom user data, can be accessed in callback functions
		void *user_data;
		void(*update_callback)(tfxParticleManager *pm, tfxEffectID effect_index);

		tfxU32 buffer_index;

		//Indexes into library storage
		tfxU32 info_index;
		tfxU32 property_index;
		tfxU32 pm_index;

		tfxEffectEmitter() :
			buffer_index(0),
			pm_index(0),
			parent(nullptr),
			user_data(nullptr),
			update_callback(nullptr),
			effect_flags(tfxEffectPropertyFlags_none),
			sort_passes(1),
			info_index(tfxINVALID),
			property_index(tfxINVALID),
			global(tfxINVALID),
			emitter_attributes(tfxINVALID),
			transform_attributes(tfxINVALID),
			property_flags(tfxEmitterPropertyFlags_image_handle_auto_center |
				tfxEmitterPropertyFlags_grid_spawn_clockwise |
				tfxEmitterPropertyFlags_emitter_handle_auto_center |
				tfxEmitterPropertyFlags_global_uniform_size |
				tfxEmitterPropertyFlags_base_uniform_size |
				tfxEmitterPropertyFlags_lifetime_uniform_size),
			state_flags(0)
		{ }
		~tfxEffectEmitter();

		//API related functions

		void SetUserData(void *data);
		void *GetUserData();

		tfxEffectEmitterInfo &GetInfo();
		tfxEmitterPropertiesSoA &GetProperties();

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
		void ResetTransformGraphs(bool add_node = true, bool compile = true);
		void ResetBaseGraphs(bool add_node = true, bool compile = true);
		void ResetPropertyGraphs(bool add_node = true, bool compile = true);
		void ResetVariationGraphs(bool add_node = true, bool compile = true);
		void ResetOvertimeGraphs(bool add_node = true, bool compile = true);
		void ResetEffectGraphs(bool add_node = true, bool compile = true);
		void ResetEmitterGraphs(bool add_node = true, bool compile = true);
		void UpdateMaxLife();
		void ResetAllBufferSizes();
		tfxGraph* GetGraphByType(tfxGraphType type);
		tfxU32 GetGraphIndexByType(tfxGraphType type);
		void CompileGraphs();
		void InitialiseUninitialisedGraphs();
		void SetName(const char *n);

		bool HasSingle();
		bool RenameSubEffector(tfxEffectEmitter &effect, const char *new_name);
		bool NameExists(tfxEffectEmitter &effect, const char *name);
		void FreeGraphs();

		void ClearColors();
		void AddColorOvertime(float frame, tfxRGB color);
		void Clone(tfxEffectEmitter &clone, tfxEffectEmitter *root_parent, tfxLibrary *destination_library, tfxEffectCloningFlags flags = 0);
		void EnableAllEmitters();
		void EnableEmitter();
		void DisableAllEmitters();
		void DisableAllEmittersExcept(tfxEffectEmitter &emitter);
		bool IsFinite();
		void FlagAs3D(bool flag);
		bool Is3DEffect();
		tfxU32 CountAllLookupValues();
		tfxParticleManagerModes GetRequiredParticleManagerMode();
		tfxPreviewCameraSettings &GetCameraSettings();
		tfxU32 GetCameraSettingsIndex();

	};

	inline int GetDepth(tfxEffectEmitter &e) {
		tfxEffectEmitter *current_parent = e.parent;
		int depth = 0;
		while (current_parent) {
			if (current_parent->type == tfxEmitterType) {
				depth++;
			}
			current_parent = current_parent->parent;
		}
		return depth;
	}

	inline tfxU32 CountAllEffects(tfxEffectEmitter &effect, tfxU32 amount = 0) {
		for (auto &sub : effect.GetInfo().sub_effectors) {
			amount = CountAllEffects(sub, amount);
		}
		return ++amount;
	}

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

	//this is just used in sorting to store a temporary copy of the particle data
	struct tfxParticleTemp {
		tfxU32 parent_index;
		tfxU32 sprite_index;
		tfxU32 particle_index;
		tfxParticleFlags flags;
		float age;
		float max_age;
		float position_x;
		float position_y;
		float position_z;
		float captured_position_x;
		float captured_position_y;
		float captured_position_z;
		float local_rotations_x;
		float local_rotations_y;
		float local_rotations_z;
		tfxU32 velocity_normal;
		float depth;
		float base_weight;
		float base_velocity;
		float base_spin;
		float noise_offset;
		float noise_resolution;
		tfxRGBA8 color;
		float alpha;
		float image_frame;
		float base_size_x;
		float base_size_y;
		tfxU32 single_loop_count;
	};

	struct tfxVec3SoA {
		float *x;
		float *y;
		float *z;
	};

	inline void InitVec3SoA(tfxSoABuffer *buffer, tfxVec3SoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(float), offsetof(tfxVec3SoA, x));
		AddStructArray(buffer, sizeof(float), offsetof(tfxVec3SoA, y));
		AddStructArray(buffer, sizeof(float), offsetof(tfxVec3SoA, z));
		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	//These all point into a tfxSoABuffer, initialised with InitParticleSoA. Current Bandwidth: 108 bytes
	struct tfxParticleSoA {
		tfxU32 *parent_index;
		tfxU32 *sprite_index;
		tfxU32 *particle_index;
		tfxParticleFlags *flags;
		float *age;
		float *max_age;
		float *position_x;
		float *position_y;
		float *position_z;
		float *captured_position_x;
		float *captured_position_y;
		float *captured_position_z;
		float *local_rotations_x;
		float *local_rotations_y;
		float *local_rotations_z;
		tfxU32 *velocity_normal;
		float *depth;
		float *base_weight;
		float *base_velocity;
		float *base_spin;
		float *base_size_x;
		float *base_size_y;
		float *noise_offset;
		float *noise_resolution;
		tfxRGBA8 *color;
		float *image_frame;
		tfxU32 *single_loop_count;
	};

	inline void InitParticleSoA(tfxSoABuffer *buffer, tfxParticleSoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxParticleSoA, parent_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxParticleSoA, sprite_index));
		AddStructArray(buffer, sizeof(tfxParticleID), offsetof(tfxParticleSoA, particle_index));
		AddStructArray(buffer, sizeof(tfxParticleFlags), offsetof(tfxParticleSoA, flags));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, age));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, max_age));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, position_x));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, position_y));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, position_z));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, captured_position_x));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, captured_position_y));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, captured_position_z));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, local_rotations_x));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, local_rotations_y));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, local_rotations_z));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxParticleSoA, velocity_normal));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, depth));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, base_weight));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, base_velocity));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, base_spin));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, noise_offset));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, noise_resolution));
		AddStructArray(buffer, sizeof(tfxRGBA8), offsetof(tfxParticleSoA, color));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, image_frame));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, base_size_x));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, base_size_y));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxParticleSoA, single_loop_count));
		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	struct tfxSpriteTransform2d {
		tfxVec2 position;					//The position of the sprite, x, y - world, z, w = captured for interpolating
		tfxVec2 scale;						//Scale
		float rotation;
	};

	struct tfxSpriteTransform3d {
		tfxVec3 position;					//The position of the sprite, x, y - world, z, w = captured for interpolating
		tfxVec3 rotations;					//Rotations of the sprite
		tfxVec2 scale;						//Scale
	};

	//When exporting effects as sprite data each frame gets frame meta containing information about the frame such as bounding box and sprite count/offset into the buffer
	struct tfxFrameMeta {
		tfxU32 frame_index;					//The index of the frame of animation
		tfxU32 index_offset[tfxLAYERS];		//All sprite data is contained in a single buffer and this is the offset to the first sprite in the range
		tfxU32 sprite_count[tfxLAYERS];		//The number of sprites in the frame
		tfxVec3 corner1;					//Bounding box corner
		tfxVec3 corner2;					//The bounding box can be used to decide if this frame needs to be drawn
	};

	struct tfxSprite3d {	//56 bytes
		tfxU32 image_frame_plus;	//The image frame of animation index packed with alignment option flag and property_index
		tfxU32 captured_index;
		tfxSpriteTransform3d transform;
		tfxU32 alignment;			//normalised alignment vector 3 floats packed into 10bits each with 2 bits left over
		tfxRGBA8 color;				//The color tint of the sprite and blend factor in a
		float stretch;
		float intensity;
	};

	//This struct of arrays is used for both 2d and 3d sprites, but obviously the transform_3d data is either 2d or 3d depending on which effects you're using in the particle manager.
	//InitSprite3dSoA is called to initialise 3d sprites and InitSprite2dArray for 2d sprites. This is all managed internally by the particle manager. It's convenient to have both 2d and
	//3d in one struct like this as it makes it a lot easier to use the same control functions where we can.
	struct tfxSpriteSoA {	//3d takes 56 bytes of bandwidth, 2d takes 40 bytes of bandwidth
		tfxU32 *image_frame_plus;				//The image frame of animation index packed with alignment option flag and property_index
		tfxU32 *captured_index;					//The index of the sprite in the previous frame so that it can be looked up and interpolated with
		tfxSpriteTransform3d *transform_3d;		//Transform data for 3d sprites
		tfxSpriteTransform2d *transform_2d;		//Transform data for 2d sprites
		tfxU32 *alignment;						//normalised alignment vector 3 floats packed into 10bits each with 2 bits left over or 2 packed 16bit floats for 2d
		tfxRGBA8 *color;						//The color tint of the sprite and blend factor in alpha channel
		float *stretch;							//Multiplier for how much the particle is stretched in the shader (3d only)	
		float *intensity;						//The multiplier for the sprite color
	};

	inline void InitSprite3dSoA(tfxSoABuffer *buffer, tfxSpriteSoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteSoA, image_frame_plus));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteSoA, captured_index));
		AddStructArray(buffer, sizeof(tfxSpriteTransform3d), offsetof(tfxSpriteSoA, transform_3d));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteSoA, alignment));
		AddStructArray(buffer, sizeof(tfxRGBA8), offsetof(tfxSpriteSoA, color));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteSoA, stretch));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteSoA, intensity));
		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	inline void InitSpriteBothSoA(tfxSoABuffer *buffer, tfxSpriteSoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteSoA, image_frame_plus));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteSoA, captured_index));
		AddStructArray(buffer, sizeof(tfxSpriteTransform2d), offsetof(tfxSpriteSoA, transform_2d));
		AddStructArray(buffer, sizeof(tfxSpriteTransform3d), offsetof(tfxSpriteSoA, transform_3d));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteSoA, alignment));
		AddStructArray(buffer, sizeof(tfxRGBA8), offsetof(tfxSpriteSoA, color));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteSoA, stretch));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteSoA, intensity));
		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	inline void InitSprite2dSoA(tfxSoABuffer *buffer, tfxSpriteSoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteSoA, image_frame_plus));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteSoA, captured_index));
		AddStructArray(buffer, sizeof(tfxSpriteTransform2d), offsetof(tfxSpriteSoA, transform_2d));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteSoA, alignment));
		AddStructArray(buffer, sizeof(tfxRGBA8), offsetof(tfxSpriteSoA, color));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteSoA, stretch));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteSoA, intensity));
		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	struct tfxSpriteData3dSoA {	//56 bytes
		tfxU32 *image_frame_plus;	//The image frame of animation index packed with alignment option flag and property_index
		tfxU32 *captured_index;
		tfxU32 *captured_index_offset;
		tfxSpriteTransform3d *transform;
		tfxU32 *alignment;			//normalised alignment vector 3 floats packed into 10bits each with 2 bits left over
		tfxRGBA8 *color;				//The color tint of the sprite and blend factor in a
		float *stretch;
		float *intensity;
	};

	inline void InitSpriteData3dSoA(tfxSoABuffer *buffer, tfxSpriteData3dSoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteData3dSoA, image_frame_plus));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteData3dSoA, captured_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteData3dSoA, captured_index_offset));
		AddStructArray(buffer, sizeof(tfxSpriteTransform3d), offsetof(tfxSpriteData3dSoA, transform));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteData3dSoA, alignment));
		AddStructArray(buffer, sizeof(tfxRGBA8), offsetof(tfxSpriteData3dSoA, color));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteData3dSoA, stretch));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteData3dSoA, intensity));
		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	struct tfxWideLerpTransformResult {
		float position[3];
		float rotations[3];
		float scale[2];
	};

	struct tfxWideLerpOtherResult {
		float stretch;
		float intensity;
		float color[4];
		float padding[2];
	};

	inline tfxWideLerpTransformResult InterpolateSpriteTransform(const tfxWideFloat &tween, const tfxSpriteTransform3d &current, const tfxSpriteTransform3d &captured) {
#ifdef tfxUSEAVX
		tfxWideFloat to = tfxWideLoad(&current.position.x);
		tfxWideFloat from = tfxWideLoad(&captured.position.x);
		tfxWideFloat one_minus_tween = tfxWideSub(tfxWIDEONE, tween);
		tfxWideFloat to_lerp = tfxWideMul(to, tween);
		tfxWideFloat from_lerp = tfxWideMul(from, one_minus_tween);
		tfxWideFloat result = tfxWideAdd(from_lerp, to_lerp);
		tfxWideLerpTransformResult out;
		tfxWideStore(out.position, result);
		return out;
#else
		tfxWideFloat to1 = tfxWideLoad(&current.position.x);
		tfxWideFloat from1 = tfxWideLoad(&captured.position.x);
		tfxWideFloat to2 = tfxWideLoad(&current.rotations.y);
		tfxWideFloat from2 = tfxWideLoad(&captured.rotations.y);
		tfxWideFloat one_minus_tween = tfxWideSub(tfxWIDEONE, tween);
		tfxWideFloat to_lerp1 = tfxWideMul(to1, tween);
		tfxWideFloat from_lerp1 = tfxWideMul(from1, one_minus_tween);
		tfxWideFloat result = tfxWideAdd(from_lerp1, to_lerp1);
		tfxWideLerpTransformResult out;
		tfxWideStore(out.position, result);
		to_lerp1 = tfxWideMul(to2, tween);
		from_lerp1 = tfxWideMul(from2, one_minus_tween);
		result = tfxWideAdd(from_lerp1, to_lerp1);
		tfxWideStore(&out.rotations[1], result);
		return out;
#endif
	}

	struct tfxSpriteData {
		tfxU32 frame_count;
		float frame_compression;
		float animation_length_in_time;		//measured in millesecs
		tfxU32 total_sprites;
		tfxU32 total_memory_for_sprites;
		tfxSoABuffer sprites_buffer;
		tfxSpriteData3dSoA sprites;
		tfxArray<tfxFrameMeta> frame_meta;
	};

	inline void FreeSpriteData(tfxSpriteData &sprite_data) {
		FreeSoABuffer(&sprite_data.sprites_buffer);
		sprite_data.frame_meta.free();
	}

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

		float noise_offset = 1;					//The random velocity added each frame
		float noise_resolution = 1;				//The random velocity added each frame
		float image_frame = 0;
		tfxU32 control_slot_and_layer;	//index to the controller, and also stores the layer in the particle manager that the particle is on (layer << 3)
		float local_rotation;
	};

	struct tfxComputeImageData {
		tfxVec4 uv;
		tfxVec2 image_size;
		tfxU32 image_index = 0;
		float animation_frames = 0;
		//float max_radius;
	};

	//Struct to contain a static state of a particle in a frame of animation. Used in the editor for recording frames of animation so probably not needed here really!
	struct tfxParticleFrame {
		tfxU32 image_frame_plus;	//The image frame of animation index packed with alignment option flag and property_index
		tfxU32 captured_index;
		tfxSpriteTransform3d transform;
		tfxU32 alignment;			//normalised alignment vector 3 floats packed into 10bits each with 2 bits left over
		tfxRGBA8 color;				//The color tint of the sprite and blend factor in a
		float stretch;
		float intensity;
		float depth;
	};

	struct tfxSpawnWorkEntry {
		tfxParticleManager *pm;
		tfxEmitterPropertiesSoA *properties;
		tfxU32 emitter_index;
		tfxParticleSoA *particle_data;
		tfxvec<tfxEffectEmitter> *sub_effects;
		tfxU32 seed;
		float tween;
		tfxU32 max_spawn_count;
		tfxU32 amount_to_spawn = 0;
		tfxU32 end_index;
		tfxU32 spawn_start_index;
		tfxU32 next_buffer;
		int depth;
		float qty_step_size;
		float highest_particle_age;
	};

	struct tfxControlWorkEntry {
		tfxU32 start_index;
		tfxU32 end_index;
		tfxU32 wide_end_index;
		tfxU32 start_diff;
		tfxU32 sprites_index;
		tfxU32 sprite_buffer_end_index;
		tfxU32 emitter_index;
		tfxParticleManager *pm;
		tfxOvertimeAttributes *graphs;
		tfxU32 layer;
		tfxEmitterPropertiesSoA *properties;
		tfxSpriteSoA *sprites;
	};

	struct tfxControlWorkEntryOrdered {
		tfxParticleManager *pm;
		tfxU32 sprite_layer;
		tfxU32 current_buffer_index;
		tfxU32 next_buffer_index;
		tfxU32 amount_to_update;
		tfxU32 start_index;
		tfxU32 sprite_start_index;
		tfxU32 end_index;
		tfxU32 wide_end_index;
		tfxU32 start_diff;
		bool calculate_depth;
		tfxSpriteSoA *sprites;
	};

	struct tfxParticleAgeWorkEntry {
		tfxU32 start_index;
		tfxU32 emitter_index;
		tfxU32 wide_end_index;
		tfxU32 start_diff;
		tfxEmitterPropertiesSoA *properties;
		tfxParticleManager *pm;
	};

	struct tfxSortWorkEntry {
		tfxParticleManager *pm;
		tfxParticleSoA *particles;
		tfxU32 current_buffer_index;
		tfxU32 size;
		int start_index;
		int end_index;
	};

	struct tfxEffectData {
		tfxU32 *global_attributes;
		tfxU32 *transform_attributes;
		float *overal_scale;
		float *life;
		float *size_x;
		float *size_y;
		float *velocity;
		float *spin;
		float *intensity;
		float *splatter;
		float *weight;
	};

	struct tfxEffectsInUseSoA {
		tfxU32 *effects_in_use[3][2];
		tfxU32 *emitters_in_use[3][2];
		tfxU32 *free_effects[3];
		tfxU32 *free_emitters[3];
	};

	inline void InitEffectsInUse(tfxSoABuffer *buffer, tfxEffectsInUseSoA *soa, tfxU32 reserve_amount) {

		for (int i = 0; i != 3; ++i) {
			for (int j = 0; j != 2; ++j) {
				AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEffectsInUseSoA, effects_in_use[i][j]));
			}
			AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEffectsInUseSoA, free_effects[i]));
		}

		for (int i = 0; i != 3; ++i) {
			for (int j = 0; j != 2; ++j) {
				AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEffectsInUseSoA, emitters_in_use[i][j]));
			}
			AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEffectsInUseSoA, free_emitters[i]));
		}

		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	inline void tfxResizeParticleSoACallback(tfxSoABuffer *buffer, tfxU32 index) {
		tfxParticleSoA *particles = static_cast<tfxParticleSoA*>(buffer->user_data);
		for (int i = index; i != buffer->capacity; ++i) {
			particles->max_age[i] = 1.f;
			particles->age[i] = 1.f;
			particles->flags[i] = 0;
		}
	}

	//Use the particle manager to add multiple effects to your scene 
	struct tfxParticleManager {
		//In unordered mode, emitters get their own list of particles to update
		tfxvec<tfxSoABuffer> particle_array_buffers;
		tfxBucketArray<tfxParticleSoA> particle_arrays;

		tfxMemoryArenaManager particle_array_allocator;
		//In unordered mode emitters that expire have their particle banks added here to be reused
		tfxStorageMap<tfxvec<tfxU32>> free_particle_lists;
		//Only used when using distance from camera ordering. New particles are put in this list and then merge sorted into the particles buffer
		tfxControlWorkEntryOrdered ordered_age_work_entry[tfxLAYERS * 2];
		tfxSortWorkEntry sorting_work_entry[tfxLAYERS];

		tfxvec<tfxParticleID> particle_indexes;
		tfxvec<tfxU32> free_particle_indexes;
		tfxvec<tfxU32> effects_in_use[tfxMAXDEPTH][2];
		tfxvec<tfxU32> emitters_in_use[tfxMAXDEPTH][2];
		tfxvec<tfxU32> free_effects;
		tfxvec<tfxU32> free_emitters;
		tfxSoABuffer effect_buffers;
		tfxEffectSoA effects;
		tfxSoABuffer emitter_buffers;
		tfxEmitterSoA emitters;
		tfxLibrary *library;

		tfxWorkQueue work_queue;

		//Banks of sprites for drawing in unordered mode
		tfxSoABuffer sprite_buffer[2][tfxLAYERS];
		tfxSpriteSoA sprites[2][tfxLAYERS];
		tfxU32 current_sprite_buffer;

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
		//The current effect buffer in use, can be either 0 or 1
		unsigned int current_ebuff;
		unsigned int next_ebuff;
		//When using depth sorting in 3d, the particles are double buffered
		unsigned int current_pbuff;
		tfxU32 effects_start_size[tfxMAXDEPTH];
		tfxU32 emitter_start_size[tfxMAXDEPTH];

		tfxU32 sprite_index_point[tfxLAYERS];
		tfxU32 new_particles_index_start[tfxLAYERS];

		int mt_batch_size;

		tfxRandom random;
		unsigned int max_compute_controllers;
		unsigned int highest_compute_controller_index;
		tfxComputeFXGlobalState compute_global_state;
		tfxU32 sort_passes;
		tfxLookupMode lookup_mode;
		//For when particles are ordered by distance from camera (3d effects)
		tfxVec3 camera_front;
		tfxVec3 camera_position;

		tfxU32 temp_count = 100;

		//These can possibly be removed at some point, they're debugging variables
		unsigned int particle_id;
		tfxParticleManagerFlags flags;

		tfxParticleManager() :
			flags(0),
			lookup_mode(tfxFast),
			max_effects(10000),
			current_ebuff(0),
			current_pbuff(0),
			highest_compute_controller_index(0),
			new_compute_particle_ptr(NULL),
			compute_controller_ptr(NULL),
			max_compute_controllers(10000),
			max_new_compute_particles(10000),
			new_compute_particle_index(0),
			new_particles_count(0),
			mt_batch_size(512),
			current_sprite_buffer(0),
			free_compute_controllers(tfxCONSTRUCTOR_VEC_INIT(pm "free_comput_controllers")),
			library(NULL),
			sort_passes(0)
		{
		}
		~tfxParticleManager();

		//Initialise the particle manager with the maximum number of particles and effects that you want the manager to update per frame
		void Reconfigure(tfxParticleManagerModes mode, tfxU32 sort_passes, bool is_3d);
		void InitForBoth(tfxLibrary *lib, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit = 1000, tfxParticleManagerModes mode = tfxParticleManagerMode_unordered, bool double_buffer_sprites = true, bool dynamic_sprite_allocation = false, tfxU32 multi_threaded_batch_size = 512);
		void InitFor2d(tfxLibrary *lib, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit = 1000, tfxParticleManagerModes mode = tfxParticleManagerMode_unordered, bool double_buffer_sprites = true, bool dynamic_sprite_allocation = false, tfxU32 multi_threaded_batch_size = 512);
		void InitFor3d(tfxLibrary *lib, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit = 1000, tfxParticleManagerModes mode = tfxParticleManagerMode_unordered, bool double_buffer_sprites = true, bool dynamic_sprite_allocation = false, tfxU32 multi_threaded_batch_size = 512);
		void InitFor2d(tfxLibrary *lib, unsigned int effects_limit = 1000, tfxParticleManagerModes mode = tfxParticleManagerMode_unordered);
		void InitFor3d(tfxLibrary *lib, unsigned int effects_limit = 1000, tfxParticleManagerModes mode = tfxParticleManagerMode_unordered);
		inline void SetLibrary(tfxLibrary *lib) {
			library = lib;
		}
		void CreateParticleBanksForEachLayer();
		//Update the particle manager. Call this once per frame in your logic udpate.
		void Update();
		//When paused you still might want to keep the particles in order:
		void UpdateParticleOrderOnly();
		//Add an effect to the particle manager. Pass a tfxEffectEmitter pointer if you want to change the effect on the fly. Once you add the effect to the particle manager
		//then it's location in the buffer will keep changing as effects are updated and added and removed. The tracker will be updated accordingly each frame so you will always
		//have access to the effect if you need it.
		tfxU32 AddEffect(tfxEffectEmitter &effect, int buffer, int depth = 0, bool is_sub_effect = false, float add_delayed_spawning = 0);
		tfxU32 AddEffect(tfxEffectTemplate &effect);
		inline void UpdateAgeOnly(bool switch_on) { if (switch_on) flags |= tfxEffectManagerFlags_update_age_only; else flags &= ~tfxEffectManagerFlags_update_age_only; }
		inline void ForceSingleThreaded(bool switch_on) { if (switch_on) flags |= tfxEffectManagerFlags_single_threaded; else flags &= ~tfxEffectManagerFlags_single_threaded; }
		inline tfxU32 GetEffectSlot() {
			if (!free_effects.empty()) {
				return free_effects.pop_back();
			}
			if (effect_buffers.current_size == effect_buffers.capacity)
				return tfxINVALID;
			AddRow(&effect_buffers);
			return effect_buffers.current_size - 1;
		}
		inline tfxU32 GetEmitterSlot() {
			if (!free_emitters.empty()) {
				return free_emitters.pop_back();
			}
			if (emitter_buffers.current_size == emitter_buffers.capacity) {
				return tfxINVALID;
			}
			AddRow(&emitter_buffers);
			return emitter_buffers.current_size - 1;
		}
		inline tfxU32 GetParticleIndexSlot(tfxParticleID particle_id) {
			if (!free_particle_indexes.empty()) {
				particle_indexes[free_particle_indexes.back()] = particle_id;
				return free_particle_indexes.pop_back();
			}
			particle_indexes.push_back(particle_id);
			return particle_indexes.current_size - 1;
		}
		inline void FreeParticleIndex(tfxU32 &index) {
			particle_indexes[index] = tfxINVALID;
			free_particle_indexes.push_back(index);
			index = tfxINVALID;
		}
		void FreeParticleList(tfxU32 index);
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
		tfxAPI inline tfxU32 PreviousSpriteBuffer() {
			assert(flags & tfxEffectManagerFlags_double_buffer_sprites);		//Particle manager must have double buffered sprites activated
			return !current_sprite_buffer;
		}
		tfxAPI inline tfxVec3 &GetCapturedSprite3dPosition(tfxU32 layer, tfxU32 index) {
			return sprites[!current_sprite_buffer][layer].transform_3d[index].position;
		}
		tfxAPI inline tfxSpriteTransform3d &GetCapturedSprite3dTransform(tfxU32 layer, tfxU32 index) {
			return sprites[(index & 0xF0000000) >> 28][layer].transform_3d[index & 0x0FFFFFFF];
		}
		tfxAPI inline tfxSpriteTransform2d &GetCapturedSprite2dTransform(tfxU32 layer, tfxU32 index) {
			return sprites[(index & 0xF0000000) >> 28][layer].transform_2d[index & 0x0FFFFFFF];
		}
		tfxAPI inline float &GetCapturedSprite3dIntensity(tfxU32 layer, tfxU32 index) {
			return sprites[(index & 0xF0000000) >> 28][layer].intensity[index & 0x0FFFFFFF];
		}

		//Internal use only
		int AddComputeController();
		inline void FreeComputeSlot(unsigned int slot_id) { free_compute_controllers.push_back(slot_id); }
		void EnableCompute() { flags |= tfxEffectManagerFlags_use_compute_shader; }
		void DisableCompute() { flags &= ~tfxEffectManagerFlags_use_compute_shader; }

		inline tfxU32 &GetParticleSpriteIndex(tfxParticleID id) { return particle_arrays[ParticleBank(id)].sprite_index[ParticleIndex(id)]; }

		tfxComputeParticle &GrabComputeParticle(unsigned int layer);
		void ResetParticlePtr(void *ptr);
		void ResetControllerPtr(void *ptr);
		inline unsigned int GetControllerMemoryUsage() { return highest_compute_controller_index * sizeof(tfxComputeController); }
		inline unsigned int GetParticleMemoryUsage() { return new_compute_particle_index * sizeof(tfxComputeParticle); }
		void UpdateCompute(void *sampled_particles, unsigned int sample_size = 100);
		//float Record(unsigned int frames, unsigned int start_frame, std::array<tfxvec<ParticleFrame>, 1000> &particle_frames);
		void UpdateBaseValues();
		tfxvec<tfxU32> *GetEffectBuffer(tfxU32 depth = 0);
		tfxvec<tfxU32> *GetEmitterBuffer(tfxU32 depth = 0);
		void SetLookUpMode(tfxLookupMode mode);
		inline tfxParticleID SetNextParticle(tfxU32 next_index, tfxU32 current_index, tfxU32 other_index) {
			tfxParticleSoA &to_bank = particle_arrays[next_index];
			tfxParticleSoA &from_bank = particle_arrays[current_index];
			tfxU32 index = particle_array_buffers[next_index].current_size++;
			assert(index < particle_array_buffers[next_index].capacity);
			to_bank.parent_index[index] = from_bank.parent_index[other_index];
			to_bank.sprite_index[index] = from_bank.sprite_index[other_index];
			to_bank.particle_index[index] = from_bank.particle_index[other_index];
			to_bank.flags[index] = from_bank.flags[other_index];
			to_bank.age[index] = from_bank.age[other_index];
			to_bank.max_age[index] = from_bank.max_age[other_index];
			to_bank.position_x[index] = from_bank.position_x[other_index];
			to_bank.position_y[index] = from_bank.position_y[other_index];
			to_bank.position_z[index] = from_bank.position_z[other_index];
			to_bank.captured_position_x[index] = from_bank.captured_position_x[other_index];
			to_bank.captured_position_y[index] = from_bank.captured_position_y[other_index];
			to_bank.captured_position_z[index] = from_bank.captured_position_z[other_index];
			to_bank.local_rotations_x[index] = from_bank.local_rotations_x[other_index];
			to_bank.local_rotations_y[index] = from_bank.local_rotations_y[other_index];
			to_bank.local_rotations_z[index] = from_bank.local_rotations_z[other_index];
			to_bank.velocity_normal[index] = from_bank.velocity_normal[other_index];
			to_bank.depth[index] = from_bank.depth[other_index];
			to_bank.base_weight[index] = from_bank.base_weight[other_index];
			to_bank.base_velocity[index] = from_bank.base_velocity[other_index];
			to_bank.base_spin[index] = from_bank.base_spin[other_index];
			to_bank.noise_offset[index] = from_bank.noise_offset[other_index];
			to_bank.noise_resolution[index] = from_bank.noise_resolution[other_index];
			to_bank.color[index] = from_bank.color[other_index];
			to_bank.image_frame[index] = from_bank.image_frame[other_index];
			to_bank.base_size_x[index] = from_bank.base_size_x[other_index];
			to_bank.base_size_y[index] = from_bank.base_size_y[other_index];
			to_bank.single_loop_count[index] = from_bank.single_loop_count[other_index];
			return MakeParticleID(next_index, index);
		}

		inline bool FreeCapacity(int index, bool compute) {
			if (!compute) {
				return sprite_buffer[current_sprite_buffer][index].current_size < max_cpu_particles_per_layer[index] || flags & tfxEffectManagerFlags_dynamic_sprite_allocation;
			}
			else
				return new_compute_particle_index < max_new_compute_particles && new_compute_particle_index < compute_global_state.end_index - compute_global_state.current_length;
		}

		inline bool FreeEffectCapacity() {
			return emitter_buffers.current_size < max_effects;
		}
		inline tfxU32 ParticleCount() {
			tfxU32 count = 0;
			for (tfxEachLayer) {
				count += sprite_buffer[current_sprite_buffer][layer].current_size;
			}
			return count;
		}
	};

	inline void DumpSprites(tfxParticleManager &pm, tfxU32 layer) {
		for (int i = 0; i != pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size; ++i) {
			printf("%i:\t%f\t%f\t%f\t%u\n",
				i,
				pm.sprites[pm.current_sprite_buffer][layer].transform_3d[i].position.x,
				pm.sprites[pm.current_sprite_buffer][layer].transform_3d[i].position.y,
				pm.sprites[pm.current_sprite_buffer][layer].transform_3d[i].position.z,
				pm.sprites[pm.current_sprite_buffer][layer].image_frame_plus[i]
			);
		}
	}

	tfxU32 GrabParticleLists(tfxParticleManager &pm, tfxKey emitter_hash, tfxU32 reserve_amount = 100);

	void StopSpawning(tfxParticleManager &pm);
	void RemoveAllEffects(tfxParticleManager &pm);
	void AddEffect(tfxParticleManager &pm, tfxEffectEmitter &effect, tfxVec3 position);
	//Set the user data of an effect in the particle Manager. Not guaranteed that your effect index is actually in use
	void SetEffectUserData(tfxParticleManager &pm, tfxU32 effect_index, void *data);

	void TransformEffector2d(tfxVec3 &world_rotations, tfxVec3 &local_rotations, tfxVec3 &world_position, tfxVec3 &local_position, tfxMatrix4 &matrix, tfxSpriteTransform2d &parent, bool relative_position = true, bool relative_angle = false);
	void TransformEffector3d(tfxVec3 &world_rotations, tfxVec3 &local_rotations, tfxVec3 &world_position, tfxVec3 &local_position, tfxMatrix4 &matrix, tfxSpriteTransform3d &parent, bool relative_position = true, bool relative_angle = false);
	void UpdatePMEffect(tfxParticleManager &pm, tfxU32 index, tfxU32 parent_index = tfxINVALID);
	void UpdatePMEmitter(tfxParticleManager &pm, tfxSpawnWorkEntry *spawn_work_entry);
	tfxU32 NewSpritesNeeded(tfxParticleManager &pm, tfxU32 index, tfxU32 parent_index, tfxEmitterPropertiesSoA &properties);
	void UpdateEmitterState(tfxParticleManager &pm, tfxU32 index, tfxU32 parent_index, const tfxParentSpawnControls &parent_spawn_controls, tfxSpawnWorkEntry *entry);
	void UpdateEffectState(tfxParticleManager &pm, tfxU32 index);

	void CompletePMWork(tfxParticleManager &pm);
	//Wide mt versions
	tfxU32 SpawnWideParticles2d(tfxParticleManager &pm, tfxSpawnWorkEntry &spawn_work_entry, tfxU32 max_spawn_count);
	void SpawnParticlePoint2d(tfxWorkQueue *queue, void *data);
	void SpawnParticleLine2d(tfxWorkQueue *queue, void *data);
	void SpawnParticleArea2d(tfxWorkQueue *queue, void *data);
	void SpawnParticleEllipse2d(tfxWorkQueue *queue, void *data);
	void SpawnParticleMicroUpdate2d(tfxWorkQueue *queue, void *data);
	void SpawnParticleNoise(tfxWorkQueue *queue, void *data);

	void SpawnParticleWeight(tfxWorkQueue *queue, void *data);
	void SpawnParticleVelocity(tfxWorkQueue *queue, void *data);
	void SpawnParticleRoll(tfxWorkQueue *queue, void *data);
	void SpawnParticleImageFrame(tfxWorkQueue *queue, void *data);
	void SpawnParticleSize2d(tfxWorkQueue *queue, void *data);
	void SpawnParticleAge(tfxWorkQueue *queue, void *data);
	void SpawnParticleSpin2d(tfxWorkQueue *queue, void *data);

	tfxU32 SpawnWideParticles3d(tfxParticleManager &pm, tfxSpawnWorkEntry &spawn_work_entry, tfxU32 max_spawn_count);
	void SpawnParticlePoint3d(tfxWorkQueue *queue, void *data);
	void SpawnParticleLine3d(tfxWorkQueue *queue, void *data);
	void SpawnParticleArea3d(tfxWorkQueue *queue, void *data);
	void SpawnParticleEllipse3d(tfxWorkQueue *queue, void *data);
	void SpawnParticleCylinder3d(tfxWorkQueue *queue, void *data);
	void SpawnParticleIcosphereRandom3d(tfxWorkQueue *queue, void *data);
	void SpawnParticleIcosphere3d(tfxWorkQueue *queue, void *data);
	void SpawnParticleMicroUpdate3d(tfxWorkQueue *queue, void *data);
	void SpawnParticleSpin3d(tfxWorkQueue *queue, void *data);
	void SpawnParticleSize3d(tfxWorkQueue *queue, void *data);

	void ControlParticles2d(tfxParticleManager &pm, tfxU32 emitter_index, tfxControlWorkEntry &work_entry);
	void ControlParticles3d(tfxParticleManager &pm, tfxU32 emitter_index, tfxControlWorkEntry &work_entry);
	void ControlParticlesOrdered2d(tfxParticleManager &pm, tfxControlWorkEntryOrdered &work_entry);
	void ControlParticlesOrdered3d(tfxParticleManager &pm, tfxControlWorkEntryOrdered &work_entry);

	void ControlParticleAge(tfxWorkQueue *queue, void *data);
	void ControlParticleImageFrame(tfxWorkQueue *queue, void *data);
	void ControlParticleColor(tfxWorkQueue *queue, void *data);
	void ControlParticleSize(tfxWorkQueue *queue, void *data);
	void ControlParticleOrderedAge(tfxWorkQueue *queue, void *data);
	void ControlParticleOrderedDepth(tfxWorkQueue *queue, void *data);
	void ControlParticleSizeOrdered(tfxWorkQueue *queue, void *data);
	void ControlParticleColorOrdered(tfxWorkQueue *queue, void *data);
	void ControlParticleImageFrameOrdered(tfxWorkQueue *queue, void *data);

	void ControlParticlePosition2d(tfxWorkQueue *queue, void *data);
	void ControlParticlePosition2dOld(tfxWorkQueue *queue, void *data);
	void ControlParticlePositionOrdered2d(tfxWorkQueue *queue, void *data);
	void ControlParticleTransform2dOld(tfxWorkQueue *queue, void *data);
	void ControlParticleTransform2d(tfxWorkQueue *queue, void *data);

	void ControlParticlePosition3d(tfxWorkQueue *queue, void *data);
	void ControlParticleTransform3d(tfxWorkQueue *queue, void *data);
	void ControlParticlePositionOrdered3d(tfxWorkQueue *queue, void *data);
	void ControlParticleTransformOrdered3d(tfxWorkQueue *queue, void *data);

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

	struct tfxEffectsStage {
		tfxvec<tfxEffectEmitter> effects;
	};

	struct tfxLibrary {
		tfxMemoryArenaManager graph_node_allocator;
		tfxMemoryArenaManager graph_lookup_allocator;
		tfxMemoryArenaManager sprite_data_allocator;
		tfxSoABuffer emitter_properties_buffer;

		tfxStorageMap<tfxEffectEmitter*> effect_paths;
		tfxvec<tfxEffectEmitter> effects;
		tfxStorageMap<tfxImageData> particle_shapes;
		tfxvec<tfxEffectEmitterInfo> effect_infos;
		tfxEmitterPropertiesSoA emitter_properties;
		tfxStorageMap<tfxSpriteData> pre_recorded_effects;

		tfxvec<tfxGlobalAttributes> global_graphs;
		tfxvec<tfxEmitterAttributes> emitter_attributes;
		tfxvec<tfxTransformAttributes> transform_attributes;
		tfxvec<tfxSpriteSheetSettings> sprite_sheet_settings;
		tfxvec<tfxSpriteDataSettings> sprite_data_settings;
		tfxvec<tfxPreviewCameraSettings> preview_camera_settings;
		tfxvec<tfxAttributeNode> all_nodes;
		tfxvec<tfxEffectLookUpData> node_lookup_indexes;
		tfxvec<float> compiled_lookup_values;
		tfxvec<tfxGraphLookupIndex> compiled_lookup_indexes;
		tfxvec<tfxComputeImageData> shape_data;
		//This could probably be stored globally
		tfxvec<tfxVec4> graph_min_max;

		tfxvec<tfxU32> free_global_graphs;
		tfxvec<tfxU32> free_keyframe_graphs;
		tfxvec<tfxU32> free_emitter_attributes;
		tfxvec<tfxU32> free_animation_settings;
		tfxvec<tfxU32> free_preview_camera_settings;
		tfxvec<tfxU32> free_properties;
		tfxvec<tfxU32> free_infos;
		tfxvec<tfxU32> free_keyframes;

		//Get an effect from the library by index
		tfxEffectEmitter& operator[] (uint32_t index);
		tfxStr64 name;
		bool open_library = false;
		bool dirty = false;
		tfxStr library_file_path;
		tfxU32 uid;

		tfxLibrary() :
			uid(0),
			effect_paths("EffectLib effect paths map", "EffectLib effect paths data"),
			particle_shapes("EffectLib shapes map", "EffectLib shapes data"),
			effects(tfxCONSTRUCTOR_VEC_INIT("effects")),
			effect_infos(tfxCONSTRUCTOR_VEC_INIT("effect_infos")),
			global_graphs(tfxCONSTRUCTOR_VEC_INIT("global_graphs")),
			emitter_attributes(tfxCONSTRUCTOR_VEC_INIT("emitter_attributes")),
			sprite_sheet_settings(tfxCONSTRUCTOR_VEC_INIT("animation_settings")),
			preview_camera_settings(tfxCONSTRUCTOR_VEC_INIT("preview_camera_settings")),
			all_nodes(tfxCONSTRUCTOR_VEC_INIT("all_nodes")),
			node_lookup_indexes(tfxCONSTRUCTOR_VEC_INIT("nodes_lookup_indexes")),
			compiled_lookup_values(tfxCONSTRUCTOR_VEC_INIT("compiled_lookup_values")),
			compiled_lookup_indexes(tfxCONSTRUCTOR_VEC_INIT("compiled_lookup_indexes")),
			shape_data(tfxCONSTRUCTOR_VEC_INIT("shape_data")),
			graph_min_max(tfxCONSTRUCTOR_VEC_INIT("graph_min_max")),
			free_global_graphs(tfxCONSTRUCTOR_VEC_INIT("free_global_graphs")),
			free_keyframe_graphs(tfxCONSTRUCTOR_VEC_INIT("free_keyframe_graphs")),
			free_emitter_attributes(tfxCONSTRUCTOR_VEC_INIT("free_emitter_attributes")),
			free_animation_settings(tfxCONSTRUCTOR_VEC_INIT("free_animation_settings")),
			free_preview_camera_settings(tfxCONSTRUCTOR_VEC_INIT("free_preview_camera_settings")),
			free_properties(tfxCONSTRUCTOR_VEC_INIT("free_properties")),
			free_infos(tfxCONSTRUCTOR_VEC_INIT("free_infos"))
		{}

		//Free everything in the library
		void Clear();
		void Init();
		void InitEmitterProperties();
		//Get an effect in the library by it's path. So for example, if you want to get a pointer to the emitter "spark" in effect "explosion" then you could do GetEffect("explosion/spark")
		//You will need this function to apply user data and update callbacks to effects and emitters before adding the effect to the particle manager
		tfxEffectEmitter *GetEffect(tfxStr256 &path);
		tfxEffectEmitter *GetEffect(const char *path);
		//Get an effect by it's path hash key
		tfxEffectEmitter *GetEffect(tfxKey key);
		//Get and effect by it's index
		void PrepareEffectTemplate(tfxStr256 path, tfxEffectTemplate &effect);
		void PrepareEffectTemplate(tfxEffectEmitter &effect, tfxEffectTemplate &effect_template);
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

		inline void MaybeGrowProperties(tfxU32 size_offset) {
			if (emitter_properties_buffer.current_size >= emitter_properties_buffer.capacity - size_offset) {
				GrowArrays(&emitter_properties_buffer, emitter_properties_buffer.capacity);
			}
		}

		inline void MaybeGrowInfos() {
			if (effect_infos.current_size >= effect_infos.capacity - 4) {
				effect_infos.reserve(effect_infos._grow_capacity(effect_infos.current_size + 1));
			}
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
		tfxEffectEmitter &AddStage(tfxStr64 &name);
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
		void FreeKeyframes(tfxU32 index);
		void FreeEmitterAttributes(tfxU32 index);
		void FreeProperties(tfxU32 index);
		void FreeInfo(tfxU32 index);
		tfxU32 CountKeyframeLookUpValues(tfxU32 index);
		tfxU32 CountGlobalLookUpValues(tfxU32 index);
		tfxU32 CountEmitterLookUpValues(tfxU32 index);
		tfxU32 CloneGlobal(tfxU32 source_index, tfxLibrary *destination_library);
		tfxU32 CloneKeyframes(tfxU32 source_index, tfxLibrary *destination_library);
		tfxU32 CloneEmitterAttributes(tfxU32 source_index, tfxLibrary *destination_library);
		tfxU32 CloneInfo(tfxU32 source_index, tfxLibrary *destination_library);
		tfxU32 CloneProperties(tfxU32 source_index, tfxLibrary *destination_library);
		void AddEmitterGraphs(tfxEffectEmitter &effect);
		void AddEffectGraphs(tfxEffectEmitter &effect);
		void AddTransformGraphs(tfxEffectEmitter &effect);
		tfxU32 AddSpriteSheetSettings(tfxEffectEmitter& effect);
		tfxU32 AddSpriteDataSettings(tfxEffectEmitter& effect);
		tfxU32 AddPreviewCameraSettings(tfxEffectEmitter& effect);
		tfxU32 AddPreviewCameraSettings();
		tfxU32 AddEffectEmitterInfo();
		tfxU32 AddEmitterProperties();
		tfxU32 AddKeyframes();
		void UpdateComputeNodes();
		void CompileAllGraphs();
		void CompileGlobalGraph(tfxU32 index);
		void CompileKeyframeGraph(tfxU32 index);
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

	void InvalidateNewSpriteCapturedIndex(tfxParticleManager &pm);
	void RecordSpriteData2d(tfxParticleManager &pm, tfxEffectEmitter &effect);
	void RecordSpriteData3d(tfxParticleManager &pm, tfxEffectEmitter &effect);

	struct tfxEffectTemplate {
		tfxStorageMap<tfxEffectEmitter*> paths;
		tfxEffectEmitter effect;
		tfxKey original_effect_hash;

		tfxEffectTemplate() :
			original_effect_hash(0),
			paths("Effect template paths map", "Effect template paths data")
		{}
		void AddPath(tfxEffectEmitter &effect_emitter, tfxStr256 path) {
			paths.Insert(path, &effect_emitter);
			for (auto &sub : effect_emitter.library->GetInfo(effect_emitter).sub_effectors) {
				tfxStr256 sub_path = path;
				sub_path.Appendf("/%s", sub.library->GetInfo(sub).name.c_str());
				AddPath(sub, sub_path);
			}
		}

		tfxAPI inline void Reset() {
			if (paths.Size()) {
				paths.Clear();
				effect.CleanUp();
			}
		}
		tfxAPI inline tfxEffectEmitter &Effect() { return effect; }
		tfxAPI inline tfxEffectEmitter *Get(tfxStr256 &path) { if (paths.ValidName(path)) return paths.At(path); return nullptr; }
		tfxAPI inline void SetUserData(tfxStr256 &path, void *data) { if (paths.ValidName(path)) paths.At(path)->user_data = data; }
		tfxAPI inline void SetUserData(void *data) { effect.user_data = data; }
		tfxAPI void SetUserDataAll(void *data);
		tfxAPI inline void SetEffectUpdateCallback(void(*update_callback)(tfxParticleManager *pm, tfxEffectID effect_index)) {
			effect.update_callback = update_callback;
		}
		tfxAPI inline void SetEffectUpdateCallback(tfxStr256 path, void(*update_callback)(tfxParticleManager *pm, tfxEffectID effect_index)) {
			assert(paths.ValidName(path));						//Path does not exist in library
			assert(paths.At(path)->type == tfxEffectType);		//Path must be path to an effect type
		}
		tfxAPI inline void SetEmitterUpdateCallback(tfxStr256 path, void(*update_callback)(tfxParticleManager *pm, tfxEffectID emitter_index)) {
			assert(paths.ValidName(path));						//Path does not exist in library
			assert(paths.At(path)->type == tfxEmitterType);		//Path must be a path to an emitter type
		}

		/*
		Pre-record this effect into a sprite cache so that you can play the effect back without the need to actually caclulate particles in realtime.
			* @param pm			Reference to a pm that will be used to run the particle simulation and record the sprite data
			* @param path		const *char of a path to the emitter in the effect.Must be a valid path, for example: "My Effect/My Emitter"
			* / void RecordSpriteData3d(tfxParticleManager &pm, u32 frames, u32 start_frame, int extra_frames, u32 &largest_frame);
		*/
		tfxAPI void RecordSpriteData(tfxParticleManager &pm);

		/*
		Disable an emitter within an effect. Disabling an emitter will stop it being added to the particle manager when calling AddEffectToParticleManager
		* @param path		const *char of a path to the emitter in the effect. Must be a valid path, for example: "My Effect/My Emitter"
		*/
		tfxAPI inline void DisableEmitter(const char *path) {
			assert(paths.ValidName(path));			//Must be a valid path to the emitter
			tfxEffectEmitter *emitter = paths.At(path);
			assert(emitter->type == tfxEmitterType);	//Must be an emitter that you're trying to remove. Use RemoveSubEffect if you're trying to remove one of those. 
			emitter->property_flags &= ~tfxEmitterPropertyFlags_enabled;
		}

		/*
		Enable an emitter within an effect so that it is added to the particle manager when calling AddEffectToParticleManager. Emitters are enabled by default.
		* @param path		const *char of a path to the emitter in the effect. Must be a valid path, for example: "My Effect/My Emitter"
		*/
		tfxAPI inline void EnableEmitter(const char *path) {
			assert(paths.ValidName(path));			//Must be a valid path to the emitter
			tfxEffectEmitter *emitter = paths.At(path);
			assert(emitter->type == tfxEmitterType);	//Must be an emitter that you're trying to remove. Use RemoveSubEffect if you're trying to remove one of those
			emitter->property_flags |= tfxEmitterPropertyFlags_enabled;
		}

		/*
		Scale all nodes on a global graph graph of the effect
		* @param global_type		tfxGraphType of the global graph that you want to scale. Must be a global graph or an assert will be called
		* @param amount				A float of the amount that you want to scale the multiplier by.
		*/
		tfxAPI inline void ScaleGlobalMultiplier(tfxGraphType global_type, float amount) {
			assert(IsGlobalGraph(global_type));
			tfxGraph *graph = effect.GetGraphByType(global_type);
			tfxEffectEmitter *original_effect = effect.library->GetEffect(original_effect_hash);
			tfxGraph *original_graph = original_effect->GetGraphByType(global_type);
			original_graph->Copy(*graph, false);
			graph->MultiplyAllValues(amount);
			CompileGraph(*graph);
		}

		/*
		Scale all nodes on an emitter graph
		* @param emitter_path		const *char of the emitter path
		* @param global_type		tfxGraphType of the emitter graph that you want to scale. Must be an emitter graph or an assert will be called
		* @param amount				A float of the amount that you want to scale the graph by.
		*/
		tfxAPI inline void ScaleEmitterGraph(const char *emitter_path, tfxGraphType graph_type, float amount) {
			assert(IsEmitterGraph(graph_type));		//Must be an emitter graph type. This is any property, base, variaion or overtime graph
			assert(paths.ValidName(emitter_path));			//Must be a valid path to the emitter
			tfxEffectEmitter *emitter = paths.At(emitter_path);
			tfxGraph *graph = emitter->GetGraphByType(graph_type);
			tfxEffectEmitter *original_emitter = effect.library->GetEffect(emitter_path);
			tfxGraph *original_graph = original_emitter->GetGraphByType(graph_type);
			original_graph->Copy(*graph, false);
			graph->MultiplyAllValues(amount);
			CompileGraph(*graph);
		}

		/*
		Set the single spawn amount for an emitter. Only affects emitters that have the single spawn flag set.
		* @param emitter_path		const *char of the emitter path
		* @param amount				A float of the amount that you want to set the single spawn amount to.
		*/
		tfxAPI inline void SetSingleSpawnAmount(const char *emitter_path, tfxU32 amount) {
			assert(amount >= 0);							//Amount must not be less than 0
			assert(paths.ValidName(emitter_path));			//Must be a valid path to the emitter
			tfxEffectEmitter *emitter = paths.At(emitter_path);
			emitter->GetProperties().spawn_amount[emitter->property_index] = amount;
		}

	};

	/*
	Notes on updating effects emitters and particles:

	Todo: rewrite now that we've converted to SoA data layouts
	*/

	struct tfxDataEntry {
		tfxDataType type = tfxSInt;
		tfxStr32 key;
		tfxStr str_value;
		int int_value = 0;
		bool bool_value = 0;
		float float_value = 0;
		double double_value = 0;
	};

	struct tfxDataTypesDictionary {
		bool initialised = false;
		tfxStorageMap<tfxDataType> names_and_types;

		tfxDataTypesDictionary() :
			names_and_types("Data Types Storage Map", "Data Types Storage Data")
		{}
		void Init();
	};

	extern tfxDataTypesDictionary tfxDataTypes;

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
	void StreamProperties(tfxEmitterPropertiesSoA &property, tfxU32 index, tfxEmitterPropertyFlags &flags, tfxStr &file);
	void StreamProperties(tfxEffectEmitter &effect, tfxStr &file);
	void StreamGraph(const char * name, tfxGraph &graph, tfxStr &file);
	void SplitStringStack(const tfxStr &s, tfxStack<tfxStr256> &pair, char delim = 61);
	void SplitStringVec(const tfxStr &s, tfxvec<tfxStr256> &pair, char delim = 61);
	bool StringIsUInt(const tfxStr &s);
	int GetDataType(const tfxStr &s);
	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, uint32_t value);
	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, float value);
	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, bool value);
	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, int value);
	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, tfxStr &value);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, uint32_t value);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, float value);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, bool value);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, int value);
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, tfxStr &value);
	void AssignGraphData(tfxEffectEmitter &effect, tfxStack<tfxStr256> &values);
	void AssignNodeData(tfxAttributeNode &node, tfxStack<tfxStr256> &values);
	static inline void Transform2d(tfxVec3 &out_rotations, tfxVec3 &out_local_rotations, tfxVec3 &out_scale, tfxVec3 &out_position, tfxVec3 &out_local_position, tfxVec3 &out_translation, tfxMatrix4 &out_matrix, const tfxVec3 &in_rotations, const tfxVec3 &in_scale, const tfxVec3 &in_position, const tfxMatrix4 &in_matrix) {
		float s = sin(out_local_rotations.roll);
		float c = cos(out_local_rotations.roll);

		out_matrix.Set2(c, s, -s, c);
		out_scale = in_scale;

		out_rotations.roll = in_rotations.roll + out_local_rotations.roll;

		out_matrix = mmTransform2(out_matrix, in_matrix);
		tfxVec2 rotatevec = mmTransformVector(in_matrix, tfxVec2(out_local_position.x + out_translation.x, out_local_position.y + out_translation.y));

		out_position = in_position.xy() + rotatevec * in_scale.xy();
	}
	static inline void Transform3d(tfxVec3 &out_rotations, tfxVec3 &out_local_rotations, tfxVec3 &out_scale, tfxVec3 &out_position, tfxVec3 &out_local_position, tfxVec3 &out_translation, tfxMatrix4 &out_matrix, const tfxVec3 &in_rotations, const tfxVec3 &in_scale, const tfxVec3 &in_position, const tfxMatrix4 &in_matrix) {
		tfxMatrix4 roll = mmZRotate(out_local_rotations.roll);
		tfxMatrix4 pitch = mmXRotate(out_local_rotations.pitch);
		tfxMatrix4 yaw = mmYRotate(out_local_rotations.yaw);
		out_matrix = mmTransform(yaw, pitch);
		out_matrix = mmTransform(out_matrix, roll);
		out_scale = in_scale;

		out_rotations = in_rotations + out_local_rotations;

		out_matrix = mmTransform(out_matrix, in_matrix);
		tfxVec3 rotatevec = mmTransformVector3(in_matrix, out_local_position + out_translation);

		out_position = in_position + rotatevec;
	}
	//-------------------------------------------------
	//--New transform_3d particle functions for SoA data--
	//--------------------------2d---------------------
	static inline void TransformParticlePosition(const float local_position_x, const float local_position_y, const float roll, tfxVec2 &world_position, float &world_rotations, const tfxVec3 &parent_rotations, const tfxMatrix4 &matrix, const tfxVec3 &handle, const tfxVec3 &scale, const tfxVec3 &from_position) {
		world_position.x = local_position_x;
		world_position.y = local_position_y;
		world_rotations = roll;
	}
	static inline void TransformParticlePositionAngle(const float local_position_x, const float local_position_y, const float roll, tfxVec2 &world_position, float &world_rotations, const tfxVec3 &parent_rotations, const tfxMatrix4 &matrix, const tfxVec3 &handle, const tfxVec3 &scale, const tfxVec3 &from_position) {
		world_position.x = local_position_x;
		world_position.y = local_position_y;
		world_rotations = parent_rotations.roll + roll;
	}
	static inline void TransformParticlePositionRelative(const float local_position_x, const float local_position_y, const float roll, tfxVec2 &world_position, float &world_rotations, const tfxVec3 &parent_rotations, const tfxMatrix4 &matrix, const tfxVec3 &handle, const tfxVec3 &scale, const tfxVec3 &from_position) {
		world_rotations = roll;
		tfxVec2 rotatevec = mmTransformVector(matrix, tfxVec2(local_position_x, local_position_y) + handle.xy());
		world_position = from_position.xy() + rotatevec * scale.xy();
	}
	static inline void TransformParticlePositionRelativeLine(const float local_position_x, const float local_position_y, const float roll, tfxVec2 &world_position, float &world_rotations, const tfxVec3 &parent_rotations, const tfxMatrix4 &matrix, const tfxVec3 &handle, const tfxVec3 &scale, const tfxVec3 &from_position) {
		world_rotations = parent_rotations.roll + roll;
		tfxVec2 rotatevec = mmTransformVector(matrix, tfxVec2(local_position_x, local_position_y) + handle.xy());
		world_position = from_position.xy() + rotatevec * scale.xy();
	}
	//-------------------------------------------------
	//--New transform_3d particle functions for SoA data--
	//--------------------------3d---------------------
	static inline void TransformParticlePosition3d(const float local_position_x, const float local_position_y, const float local_position_z, const tfxVec3 local_rotations, tfxVec3 &world_position, tfxVec3 &world_rotations, const tfxVec3 &parent_rotations, const tfxMatrix4 &matrix, const tfxVec3 &handle, const tfxVec3 &scale, const tfxVec3 &from_position) {
		world_position.x = local_position_x;
		world_position.y = local_position_y;
		world_position.z = local_position_z;
		world_rotations = local_rotations;
	}
	static inline void TransformParticlePositionAngle3d(const float local_position_x, const float local_position_y, const float local_position_z, const tfxVec3 local_rotations, tfxVec3 &world_position, tfxVec3 &world_rotations, const tfxVec3 &parent_rotations, const tfxMatrix4 &matrix, const tfxVec3 &handle, const tfxVec3 &scale, const tfxVec3 &from_position) {
		world_position.x = local_position_x;
		world_position.y = local_position_y;
		world_position.z = local_position_z;
		world_rotations = parent_rotations + local_rotations;
	}
	static inline void TransformParticlePositionRelative3d(const float local_position_x, const float local_position_y, const float local_position_z, const tfxVec3 local_rotations, tfxVec3 &world_position, tfxVec3 &world_rotations, const tfxVec3 &parent_rotations, const tfxMatrix4 &matrix, const tfxVec3 &handle, const tfxVec3 &scale, const tfxVec3 &from_position) {
		world_rotations = local_rotations;
		tfxVec4 rotatevec = mmTransformVector(matrix, tfxVec3(local_position_x, local_position_y, local_position_z) + handle);
		world_position = from_position + rotatevec.xyz() * scale;
	}
	static inline void TransformParticlePositionRelativeLine3d(const float local_position_x, const float local_position_y, const float local_position_z, const tfxVec3 local_rotations, tfxVec3 &world_position, tfxVec3 &world_rotations, const tfxVec3 &parent_rotations, const tfxMatrix4 &matrix, const tfxVec3 &handle, const tfxVec3 &scale, const tfxVec3 &from_position) {
		world_rotations = local_rotations;
		tfxVec4 rotatevec = mmTransformVector(matrix, tfxVec3(local_position_x, local_position_y, local_position_z) + handle);
		world_position = from_position + rotatevec.xyz() * scale;
	}

	static inline void TransformWideParticlePositionRelative3d(const float local_position_x, const float local_position_y, const float local_position_z, const tfxVec3 local_rotations, tfxVec3 &world_position, const tfxVec3 &parent_rotations, const tfxMatrix4 &matrix, const tfxVec3 &handle, const tfxVec3 &scale, const tfxVec3 &from_position) {
		tfxVec4 rotatevec = mmTransformVector(matrix, tfxVec3(local_position_x, local_position_y, local_position_z) + handle);
		world_position = from_position + rotatevec.xyz() * scale;
	}
	static inline void TransformWideParticlePositionRelativeLine3d(const float local_position_x, const float local_position_y, const float local_position_z, const tfxVec3 local_rotations, tfxVec3 &world_position, const tfxVec3 &parent_rotations, const tfxMatrix4 &matrix, const tfxVec3 &handle, const tfxVec3 &scale, const tfxVec3 &from_position) {
		tfxVec4 rotatevec = mmTransformVector(matrix, tfxVec3(local_position_x, local_position_y, local_position_z) + handle);
		world_position = from_position + rotatevec.xyz() * scale;
	}

	static inline int SortDepth(void const *left, void const *right) {
		float d1 = *static_cast<const float*>(left);
		float d2 = *static_cast<const float*>(right);
		return (d2 > d1) - (d2 < d1);
	}

	static inline int SortIcospherePoints(void const *left, void const *right) {
		float d1 = static_cast<const tfxVec3*>(left)->y;
		float d2 = static_cast<const tfxVec3*>(right)->y;
		return (d2 > d1) - (d2 < d1);
	}

	static inline void SwapSoAParticle(tfxParticleSoA &particles, tfxU32 from, tfxU32 to) {
		std::swap(particles.depth[from], particles.depth[to]);
		std::swap(particles.age[from], particles.age[to]);
		std::swap(particles.parent_index[from], particles.parent_index[to]);
		std::swap(particles.sprite_index[from], particles.sprite_index[to]);
		std::swap(particles.particle_index[from], particles.particle_index[to]);
		std::swap(particles.flags[from], particles.flags[to]);
		std::swap(particles.max_age[from], particles.max_age[to]);
		std::swap(particles.position_x[from], particles.position_x[to]);
		std::swap(particles.position_y[from], particles.position_y[to]);
		std::swap(particles.position_z[from], particles.position_z[to]);
		std::swap(particles.captured_position_x[from], particles.captured_position_x[to]);
		std::swap(particles.captured_position_y[from], particles.captured_position_y[to]);
		std::swap(particles.captured_position_z[from], particles.captured_position_z[to]);
		std::swap(particles.local_rotations_x[from], particles.local_rotations_x[to]);
		std::swap(particles.local_rotations_y[from], particles.local_rotations_y[to]);
		std::swap(particles.local_rotations_z[from], particles.local_rotations_z[to]);
		std::swap(particles.velocity_normal[from], particles.velocity_normal[to]);
		std::swap(particles.base_weight[from], particles.base_weight[to]);
		std::swap(particles.base_velocity[from], particles.base_velocity[to]);
		std::swap(particles.base_spin[from], particles.base_spin[to]);
		std::swap(particles.noise_offset[from], particles.noise_offset[to]);
		std::swap(particles.noise_resolution[from], particles.noise_resolution[to]);
		std::swap(particles.color[from], particles.color[to]);
		std::swap(particles.image_frame[from], particles.image_frame[to]);
		std::swap(particles.base_size_x[from], particles.base_size_x[to]);
		std::swap(particles.base_size_y[from], particles.base_size_y[to]);
		std::swap(particles.single_loop_count[from], particles.single_loop_count[to]);
	}

	static inline void StoreSoAParticle(tfxParticleSoA &particles, tfxU32 from, tfxParticleTemp &temp) {
		temp.depth = particles.depth[from];
		temp.age = particles.age[from];
		temp.parent_index = particles.parent_index[from];
		temp.sprite_index = particles.sprite_index[from];
		temp.particle_index = particles.particle_index[from];
		temp.flags = particles.flags[from];
		temp.max_age = particles.max_age[from];
		temp.position_x = particles.position_x[from];
		temp.position_y = particles.position_y[from];
		temp.position_z = particles.position_z[from];
		temp.captured_position_x = particles.captured_position_x[from];
		temp.captured_position_y = particles.captured_position_y[from];
		temp.captured_position_z = particles.captured_position_z[from];
		temp.local_rotations_x = particles.local_rotations_x[from];
		temp.local_rotations_y = particles.local_rotations_y[from];
		temp.local_rotations_z = particles.local_rotations_z[from];
		temp.velocity_normal = particles.velocity_normal[from];
		temp.base_weight = particles.base_weight[from];
		temp.base_velocity = particles.base_velocity[from];
		temp.base_spin = particles.base_spin[from];
		temp.noise_offset = particles.noise_offset[from];
		temp.noise_resolution = particles.noise_resolution[from];
		temp.color = particles.color[from];
		temp.image_frame = particles.image_frame[from];
		temp.base_size_x = particles.base_size_x[from];
		temp.base_size_y = particles.base_size_y[from];
		temp.single_loop_count = particles.single_loop_count[from];
	}

	static inline void LoadSoAParticle(tfxParticleSoA &particles, tfxU32 from, tfxParticleTemp &temp) {
		particles.depth[from] = temp.depth;
		particles.age[from] = temp.age;
		particles.parent_index[from] = temp.parent_index;
		particles.sprite_index[from] = temp.sprite_index;
		particles.particle_index[from] = temp.particle_index;
		particles.flags[from] = temp.flags;
		particles.max_age[from] = temp.max_age;
		particles.position_x[from] = temp.position_x;
		particles.position_y[from] = temp.position_y;
		particles.position_z[from] = temp.position_z;
		particles.captured_position_x[from] = temp.captured_position_x;
		particles.captured_position_y[from] = temp.captured_position_y;
		particles.captured_position_z[from] = temp.captured_position_z;
		particles.local_rotations_x[from] = temp.local_rotations_x;
		particles.local_rotations_y[from] = temp.local_rotations_y;
		particles.local_rotations_z[from] = temp.local_rotations_z;
		particles.velocity_normal[from] = temp.velocity_normal;
		particles.base_weight[from] = temp.base_weight;
		particles.base_velocity[from] = temp.base_velocity;
		particles.base_spin[from] = temp.base_spin;
		particles.noise_offset[from] = temp.noise_offset;
		particles.noise_resolution[from] = temp.noise_resolution;
		particles.color[from] = temp.color;
		particles.image_frame[from] = temp.image_frame;
		particles.base_size_x[from] = temp.base_size_x;
		particles.base_size_y[from] = temp.base_size_y;
		particles.single_loop_count[from] = temp.single_loop_count;
	}

	static inline void InsertionSortSoAParticles(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSortWorkEntry *work_entry = static_cast<tfxSortWorkEntry*>(data);
		tfxParticleManager &pm = *work_entry->pm;
		tfxU32 size = work_entry->size;
		tfxU32 current_buffer_index = work_entry->current_buffer_index;
		tfxParticleSoA &particles = *work_entry->particles;
		tfxParticleTemp key;
		for (tfxU32 i = 1; i < size; ++i) {
			StoreSoAParticle(particles, i, key);
			int j = i - 1;
			while (j >= 0 && key.depth > particles.depth[j]) {
				SwapSoAParticle(particles, j + 1, j);
				if (particles.flags[j + 1] & tfxParticleFlags_has_sub_effects) pm.particle_indexes[particles.particle_index[j + 1]] = MakeParticleID(current_buffer_index, j + 1);
				if (particles.flags[j] & tfxParticleFlags_has_sub_effects) pm.particle_indexes[particles.particle_index[j]] = MakeParticleID(current_buffer_index, j);
				--j;
			}
			LoadSoAParticle(particles, j + 1, key);
			if (particles.flags[j + 1] & tfxParticleFlags_has_sub_effects) pm.particle_indexes[particles.particle_index[j + 1]] = MakeParticleID(current_buffer_index, j + 1);
		}
	}

	static inline void InsertionSortSoAParticles2(tfxParticleSoA &particles, int start_index, int end_index) {
		tfxPROFILE;
		tfxParticleTemp key;
		int si = (int)start_index;
		for (int i = si + 1; i < end_index; ++i) {
			StoreSoAParticle(particles, i, key);
			int j = i - 1;
			while (j >= si && key.depth > particles.depth[j]) {
				SwapSoAParticle(particles, j + 1, j);
				--j;
			}
			LoadSoAParticle(particles, j + 1, key);
		}
	}

	static inline int QuickSortPartition(tfxParticleSoA &particles, int start_index, int end_index)
	{
		float depth = particles.depth[end_index];
		int i = start_index - 1;
		for (int j = start_index; j < end_index; j++)
		{
			if (particles.depth[j] > depth)
			{
				i++;
				SwapSoAParticle(particles, i, j);
			}
		}
		i++;
		SwapSoAParticle(particles, i, end_index);
		return i;
	}

	static inline void QuickSortSoAParticles(tfxParticleSoA &particles, int start_index, int end_index) {
		if (start_index >= end_index || start_index < 0)
			return;

		if (end_index - start_index <= 15)
		{
			//todo: This needs proper benchmarking
			InsertionSortSoAParticles2(particles, start_index, end_index);
			return;
		}

		int pivot = QuickSortPartition(particles, start_index, end_index);

		QuickSortSoAParticles(particles, start_index, pivot - 1);
		QuickSortSoAParticles(particles, pivot + 1, end_index);
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

	int ValidateEffectPackage(const char *filename);

	//Get a graph by tfxGraphID
	tfxGraph &GetGraph(tfxLibrary &library, tfxGraphID &graph_id);

	extern tfxMemoryArenaManager tfxSTACK_ALLOCATOR;
	extern tfxMemoryArenaManager tfxMT_STACK_ALLOCATOR;

	int GetShapesInPackage(const char *filename);
	int GetEffectLibraryStats(const char *filename, tfxEffectLibraryStats &stats);
	tfxEffectLibraryStats CreateLibraryStats(tfxLibrary &lib);
	tfxErrorFlags LoadEffectLibraryPackage(tfxPackage &package, tfxLibrary &lib, void(*shape_loader)(const char *filename, tfxImageData &image_data, void *raw_image_data, int image_size, void *user_data), void *user_data = nullptr, bool read_only = true);
	inline float GetUpdateTime() { return tfxUPDATE_TIME; }
	inline float GetFrameLength() { return tfxFRAME_LENGTH; }
	inline void SetLookUpFrequency(float frequency) {
		tfxLOOKUP_FREQUENCY = frequency;
	}
	inline void SetLookUpFrequencyOvertime(float frequency) {
		tfxLOOKUP_FREQUENCY_OVERTIME = frequency;
	}

	//[API functions]
	//All the functions below represent all that you will need to call to implement TimelineFX

	/*
	Initialise TimelineFX. Must be called before any functionality of TimelineFX is used.
	* @param max_threads	Pass the number of threads that you want to use in addition to the main thread.
	*						Example, if there are 12 logical cores available, 0.5 will use 6 threads. 0 means only single threaded will be used.
	*/
	tfxAPI void InitialiseTimelineFX(int max_threads = 0);

	/*
	Set the udpate frequency for all particle effects - There may be options in the future for individual effects to be updated at their own specific frequency.
	* @param fps	The target number of frames to be udpated per second. If this does not match the current update rate of your game then the particles may playback slower or faster then they should
	*/
	tfxAPI void SetUpdateFrequency(float fps);

	/*
	Get the current update frequency of timelineFX
	* @return float of the the current update frequency
	*/
	tfxAPI inline float GetUpdateFrequency() { return tfxUPDATE_FREQUENCY; }

	/**
	* Loads an effect library package from the specified filename into the provided tfxLibrary object.
	*
	* @param filename		A pointer to a null-terminated string that contains the path and filename of the effect library package to be loaded.
	* @param lib			A reference to a tfxLibrary object that will hold the loaded effect library data.
	* @param shape_loader	A pointer to a function that will be used to load image data into the effect library package.
	*						The function has the following signature: void shape_loader(const char *filename, tfxImageData &image_data, void *raw_image_data, int image_size, void *user_data).
	* @param user_data		A pointer to user-defined data that will be passed to the shape_loader function. This parameter is optional and can be set to nullptr if not needed.
	* @param read_only		A boolean value that determines whether the effect library data will be loaded in read-only mode. (Maybe removed in the future).
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
	tfxAPI tfxErrorFlags LoadEffectLibraryPackage(const char *filename, tfxLibrary &lib, void(*shape_loader)(const char *filename, tfxImageData &image_data, void *raw_image_data, int image_size, void *user_data), void *user_data = nullptr, bool read_only = true);

	//[Particle Manager functions]

	/*
	Initialise a tfxParticleManager for 3d usage
	* @param pm						A pointer to an unitialised tfxParticleManager. If you want to reconfigure a particle manager for a different usage then you can call ReconfigureParticleManager.
	* @param layer_max_values		An array of unsigned ints representing the maximum amount of particles you want available for each layer. This will allocate the appropriate amount of memory ahead of time.
	* @param effects_limit			The maximum amount of effects and emitters that can be updated in a single frame. This will allocate the appropriate amount of memory ahead of time. Default: 1000.
	* @param mode					The operation mode of the particle manager regarding how particles are ordered. Default value: tfxParticleManagerMode_unordered. Possible modes are:
		tfxParticleManagerMode_unordered					Particles will be updated by emitter. No ordering is maintained, each emitter will spawn and update their particles in turn and sprites will be ordered
															according to that sequence.
		tfxParticleManagerMode_ordered_by_age				Particles will be kept in age order, older particles will be drawn first and newer ones last
		tfxParticleManagerMode_ordered_by_depth				Particles will be drawn in depth order or distance from the camera. You can specify the number of sort passes when setting up the effects in TimelineFX editor
		tfxParticleManagerMode_ordered_by_depth_guaranteed	Particles will be sorted each update and kept in depth order
	* @param double_buffer_sprites	True or False, whether the last frame of sprites is kept so that you can use to do interpolations for smoother animation
	* @param dynamic_allocation		If set to true then when the layer_max_values is hit for a layer the sprite and particle memory allocation will be grown dynamically. This can be useful when you're unsure of how
									many particles you will need to display while developing you're game/app. Default is false.
	* @param mt_batch_size			When using multithreading you can alter the size of each batch of particles that each thread will update. The default is 512

	*/
	tfxAPI void InitParticleManagerFor3d(tfxParticleManager *pm, tfxLibrary *library, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit = 1000, tfxParticleManagerModes mode = tfxParticleManagerMode_unordered, bool double_buffer_sprites = true, bool dynamic_allocation = false, tfxU32 mt_batch_size = 512);

	/*
	Initialise a tfxParticleManager for 2d usage
	* @param pm						A pointer to an unitialised tfxParticleManager. If you want to reconfigure a particle manager for a different usage then you can call ReconfigureParticleManager.
	* @param layer_max_values		An array of unsigned ints representing the maximum amount of particles you want available for each layer. This will allocate the appropriate amount of memory ahead of time.
	* @param effects_limit			The maximum amount of effects and emitters that can be updated in a single frame. This will allocate the appropriate amount of memory ahead of time. Default: 1000.
	* @param mode					The operation mode of the particle manager regarding how particles are ordered. Default value: tfxParticleManagerMode_unordered. Possible modes are:
		tfxParticleManagerMode_unordered					Particles will be updated by emitter. No ordering is maintained, each emitter will spawn and update their particles in turn and sprites will be ordered
															according to that sequence.
		tfxParticleManagerMode_ordered_by_age				Particles will be kept in age order, older particles will be drawn first and newer ones last
	* @param double_buffer_sprites	True or False, whether the last frame of sprites is kept so that you can use to do interpolations for smoother animation
	* @param dynamic_allocation		If set to true then when the layer_max_values is hit for a layer the sprite and particle memory allocation will be grown dynamically. This can be useful when you're unsure of how
									many particles you will need to display while developing you're game/app. Default is false.
	* @param mt_batch_size			When using multithreading you can alter the size of each batch of particles that each thread will update. The default is 512.

	*/
	tfxAPI void InitParticleManagerFor2d(tfxParticleManager *pm, tfxLibrary *library, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit = 1000, tfxParticleManagerModes mode = tfxParticleManagerMode_unordered, bool double_buffer_sprites = true, bool dynamic_allocation = false, tfxU32 mt_batch_size = 512);

	/*
	Enable or disable double buffering of sprites in a particle manager. When double buffered, you can used the previous frame of sprites to interpolated sprites each frame for smoother animation. You can disable it to save memory.
	* @param pm							A pointer to an initialised tfxParticleManager. The particle manager must have already been initialised by calling InitFor3d or InitFor2d
	* @param double_buffer_sprites	True or False to activate or deactivate double buffering of sprites
	*/
	tfxAPI inline void DoubleBufferSprites(tfxParticleManager *pm, bool double_buffer_sprites) {
		if (!double_buffer_sprites && pm->flags & tfxEffectManagerFlags_double_buffer_sprites) {
			for (tfxEachLayer) {
				FreeSoABuffer(&pm->sprite_buffer[1][layer]);
			}
			pm->flags &= ~tfxEffectManagerFlags_double_buffer_sprites;
		}
		else if (double_buffer_sprites && !(pm->flags & tfxEffectManagerFlags_double_buffer_sprites)) {
			for (tfxEachLayer) {
				InitSprite3dSoA(&pm->sprite_buffer[1][layer], &pm->sprites[1][layer], tfxMax(pm->max_cpu_particles_per_layer[layer], 8));
			}
			pm->flags |= tfxEffectManagerFlags_double_buffer_sprites;
		}
	}

	/*
	Set the seed for the particle manager for random number generation. Setting the seed can determine how an emitters spawns particles, so if you set the seed before adding an effect to the particle manager
	then the effect will look the same each time. Note that seed of 0 is invalid, it must be 1 or greater.
	* @param pm							A pointer to an initialised tfxParticleManager. The particle manager must have already been initialised by calling InitFor3d or InitFor2d
	* @param seed						An unsigned int representing the seed (Any value other then 0)
	*/
	tfxAPI inline void SetSeed(tfxParticleManager *pm, tfxU64 seed) {
		pm->random.ReSeed(seed == 0 ? tfxMAX_UINT : seed);
	}

	/*
	Prepare a tfxEffectTemplate that you can use to customise effects in the library in various ways before adding them into a particle manager for updating and rendering. Using a template like this
	means that you can tweak an effect without editing the base effect in the library.
	* @param library					A reference to a tfxLibrary that should be loaded with LoadEffectLibraryPackage
	* @param name						The name of the effect in the library that you want to use for the template. If the effect is in a folder then use normal pathing: "My Folder/My effect"
	* @param effect_template			The empty tfxEffectTemplate object that you want the effect loading into
	//Returns true on success.
	*/
	tfxAPI bool PrepareEffectTemplate(tfxLibrary &library, const char *name, tfxEffectTemplate &effect_template);

	/*
	Add an effect to a tfxParticleManager
	* @param pm					A pointer to an initialised tfxParticleManager. The particle manager must have already been initialised by calling InitFor3d or InitFor2d
	* @param effect_template	The tfxEffectTemplate object that you want to add to the particle manager. It must have already been prepared by calling PrepareEffectTemplate
	*
	* @return					Index of the effect after it's been added to the particle manager. This index can then be used to manipulate the effect in the particle manager as it's update
								For example by calling SetEffectPosition
	*/
	tfxAPI tfxEffectID AddEffectToParticleManager(tfxParticleManager *pm, tfxEffectTemplate &effect);

	/*
	Update a particle manager. Call this function each frame in your update loop. It should be called the same number of times per second as set with SetUpdateFrequency.
	* @param pm					A pointer to an initialised tfxParticleManager. The particle manager must have already been initialised by calling InitFor3d or InitFor2d
	*/
	tfxAPI inline void UpdateParticleManager(tfxParticleManager *pm) {
		pm->Update();
	}

	/*
	Get the total number of 3d sprites within the layer of the particle manager
	* @param pm					A pointer to an initialised tfxParticleManager.
	* @param layer				The layer of the sprites to the count of
	*/
	tfxAPI inline tfxU32 SpritesInLayer3d(tfxParticleManager *pm, tfxU32 layer) {
		return pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size;
	}

	/*
	Get the total number of 3d sprites ready for rendering in the particle manager
	* @param pm					A pointer to an initialised tfxParticleManager.
	*/
	tfxAPI inline tfxU32 TotalSpriteCount3d(tfxParticleManager *pm) {
		tfxU32 count = 0;
		for (tfxEachLayer) {
			count += pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size;
		}
		return count;
	}

	/*
	Clear all particles and effects in a particle manager
	* @param pm					A pointer to an initialised tfxParticleManager.
	*/
	tfxAPI inline void ClearParticleManager(tfxParticleManager *pm) {
		pm->ClearAll();
	}

	//[Effects functions for altering effects that are currently playing out in a particle manager]

	/*
	Expire an effect by telling it to stop spawning particles. This means that the effect will eventually be removed from the particle manager after all of it's remaining particles have expired.
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed
	* @param effect_index	The index of the effect that you want to expire. This is the index returned when calling AddEffectToParticleManager
	*/
	tfxAPI inline void SoftExpireEffect(tfxParticleManager *pm, tfxEffectID effect_index) {
		pm->effects.state_flags[effect_index] |= tfxEmitterStateFlags_stop_spawning;
	}

	/*
	Expire an effect by telling it to stop spawning particles and remove all associated particles immediately.
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed
	* @param effect_index	The index of the effect that you want to expire. This is the index returned when calling AddEffectToParticleManager
	*/
	tfxAPI inline void HardExpireEffect(tfxParticleManager *pm, tfxEffectID effect_index) {
		pm->effects.state_flags[effect_index] |= tfxEmitterStateFlags_stop_spawning;
		pm->effects.state_flags[effect_index] |= tfxEmitterStateFlags_remove;
	}

	/*
	Get effect user data
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed
	* @param effect_index	The index of the effect that you want to expire. This is the index returned when calling AddEffectToParticleManager
	* @returns				void* pointing to the user data set in the effect. See tfxEffectTemplate::SetUserData() and SetEffectUserData()
	*/
	tfxAPI inline void* GetEffectUserData(tfxParticleManager *pm, tfxEffectID effect_index) {
		return pm->effects.user_data[effect_index];
	}

	/*
	Set the effect user data for an effect already added to a particle manager
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed
	* @param effect_index	The index of the effect that you want to expire. This is the index returned when calling AddEffectToParticleManager
	* @param user_data		A void* pointing to the user_data that you want to store in the effect
	*/
	tfxAPI inline void SetEffectUserData(tfxParticleManager *pm, tfxEffectID effect_index, void* user_data) {
		pm->effects.user_data[effect_index] = user_data;
	}

	/*
	Get the index offset into the sprite memory for sprite data containing a pre recorded effect animation. Can be used along side SpriteDataEndIndex to create
	a for loop to iterate over the sprites in a pre-recorded effect
	* @param sprite_data	A pointer to tfxSpriteData containing all the sprites and frame data
	* @param frame			The index of the frame you want the offset for
	* @param layer			The sprite layer
	* @returns				tfxU32 containing the index offset
	*/
	tfxAPI inline tfxU32 SpriteDataIndexOffset(tfxSpriteData *sprite_data, tfxU32 frame, tfxU32 layer) {
		assert(frame < sprite_data->frame_meta.size());			//frame is outside index range
		assert(layer < tfxLAYERS);								//layer is outside index range
		return sprite_data->frame_meta[frame].index_offset[layer];
	}

	/*
	Get the end index offset into the sprite memory for sprite data containing a pre recorded effect animation. Can be used along side SpriteDataIndexOffset to create
	a for loop to iterate over the sprites in a pre-recorded effect
	* @param sprite_data	A pointer to tfxSpriteData containing all the sprites and frame data
	* @param frame			The index of the frame you want the end index for
	* @param layer			The sprite layer
	* @returns				tfxU32 containing the end offset
	*/
	tfxAPI inline tfxU32 SpriteDataEndIndex(tfxSpriteData *sprite_data, tfxU32 frame, tfxU32 layer) {
		assert(frame < sprite_data->frame_meta.size());			//frame is outside index range
		assert(layer < tfxLAYERS);								//layer is outside index range
		return sprite_data->frame_meta[frame].index_offset[layer] + sprite_data->frame_meta[frame].sprite_count[layer];
	}

	/*
	Get the sprite by index from sprite data containing a pre-recorded effect. Can be used along side SpriteDataIndexOffset and SpriteDataEndIndex to create
	a for loop to iterate over the sprites in a pre-recorded effect
	* @param sprite_data	A pointer to tfxSpriteData containing all the sprites and frame data
	* @param index			The index of the sprite you want to retrieve
	* @returns				tfxSpriteSoA reference containing the sprite data for drawing
	*/
	tfxAPI inline tfxSpriteTransform3d &GetSpriteData3dTransform(tfxSpriteData3dSoA &sprites, tfxU32 index) {
		return sprites.transform[index];
	}

	/*
	Set the position of a 2d effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param x				The x value of the position
	* @param y				The y value of the position
	*/
	tfxAPI void SetEffectPosition(tfxParticleManager *pm, tfxEffectID effect_index, float x, float y);

	/*
	Set the position of a 3d effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param x				The x value of the position
	* @param y				The y value of the position
	* @param z				The z value of the position
	*/
	tfxAPI void SetEffectPosition(tfxParticleManager *pm, tfxEffectID effect_index, float x, float y, float z);

	/*
	Set the position of a 2d effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param position		A tfxVec2 vector object containing the x and y coordinates
	*/
	tfxAPI void SetEffectPosition(tfxParticleManager *pm, tfxEffectID effect_index, tfxVec2 position);

	/*
	Set the position of a 3d effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param position		A tfxVec3 vector object containing the x, y and z coordinates
	*/
	tfxAPI void SetEffectPosition(tfxParticleManager *pm, tfxEffectID effect_index, tfxVec3 position);

	/*
	Move an Effect by a specified amount relative to the effect's current position
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param amount			A tfxVec3 vector object containing the amount to move in the x, y and z planes
	*/
	tfxAPI void MoveEffect(tfxParticleManager *pm, tfxEffectID effect_index, tfxVec3 amount);

	/*
	Move an Effect by a specified amount relative to the effect's current position
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param x				The amount to move in the x plane
	* @param y				The amount to move in the y plane
	* @param z				The amount to move in the z plane
	*/
	tfxAPI void MoveEffect(tfxParticleManager *pm, tfxEffectID effect_index, float x, float y, float z);

	/*
	Get the current position of an effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @return				tfxVec3 containing the effect position
	*/
	tfxAPI tfxVec3 GetEffectPosition(tfxParticleManager *pm, tfxEffectID effect_index);

	/*
	Set the rotation of a 2d effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current rotation of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param rotation		A float of the amount that you want to set the rotation too
	*/
	tfxAPI void SetEffectRotation(tfxParticleManager *pm, tfxEffectID effect_index, float rotation);

	/*
	Set the roll of a 3d effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current roll of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param roll			A float of the amount that you want to set the roll too
	*/
	tfxAPI void SetEffectRoll(tfxParticleManager *pm, tfxEffectID effect_index, float roll);

	/*
	Set the pitch of a 3d effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current pitch of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param pitch			A float of the amount that you want to set the pitch too
	*/
	tfxAPI void SetEffectPitch(tfxParticleManager *pm, tfxEffectID effect_index, float pitch);

	/*
	Set the yaw of a 3d effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current yaw of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param yaw			A float of the amount that you want to set the yaw too
	*/
	tfxAPI void SetEffectYaw(tfxParticleManager *pm, tfxEffectID effect_index, float yaw);

	/*
	Set the width of an effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current width of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param width			A float of the amount that you want to set the width multiplier too. The width multiplier will multiply all widths of emitters within the effect so it can be an easy way to alter the size
							of area, line, ellipse etc., emitters.
	*/
	tfxAPI void SetEffectWidthMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float width);

	/*
	Set the height of an effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current height of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param height			A float of the amount that you want to set the height multiplier too. The height multiplier will multiply all heights of emitters within the effect so it can be an easy way to alter the size
							of area, line, ellipse etc., emitters.
	*/
	tfxAPI void SetEffectHeightMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float height);

	/*
	Set the depth of an effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current depth of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param depth			A float of the amount that you want to set the depth multiplier too. The depth multiplier will multiply all heights of emitters within the effect so it can be an easy way to alter the size
							of area, line, ellipse etc., emitters.
	*/
	tfxAPI void SetEffectDepthMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float depth);

	/*
	Set the life multiplier of an effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current life of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param life			A float of the amount that you want to set the life multiplier too. The life mulitplier will affect how long all particles emitted within the effect will last before expiring.
	*/
	tfxAPI void SetEffectLifeMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float life);

	/*
	Set the particle width multiplier of an effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current particle width of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param width			A float of the amount that you want to set the particle width multiplier too. The particle width mulitplier will affect the width of each particle if the emitter has a non uniform particle size, otherwise
							it will uniformly size the particle
	*/
	tfxAPI void SetEffectParticleWidthMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float width);

	/*
	Set the particle height multiplier of an effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current particle width of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param height			A float of the amount that you want to set the particle height multiplier too. The particle height mulitplier will affect the height of each particle if the emitter has a non uniform particle size, otherwise
							this function will have no effect.
	*/
	tfxAPI void SetEffectParticleHeightMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float height);

	/*
	Set the velocity multiplier of an effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current velocity of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param velocity		A float of the amount that you want to set the particle velocity multiplier too. The particle velocity mulitplier will affect the base velocity of a particle at spawn time.
	*/
	tfxAPI void SetEffectVelocityMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float velocity);

	/*
	Set the spin multiplier of an effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current spin of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param spin			A float of the amount that you want to set the particle spin multiplier too. The particle spin mulitplier will affect the base spin of a particle at spawn time.
	*/
	tfxAPI void SetEffectSpinMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float spin);

	/*
	Set the intensity multiplier of an effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current intensity of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param intensity		A float of the amount that you want to set the particle intensity multiplier too. The particle intensity mulitplier will instantly affect the opacity of all particles currently emitted by the effect.
	*/
	tfxAPI void SetEffectIntensityMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float intensity);

	/*
	Set the splatter multiplier of an effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current splatter of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param splatter		A float of the amount that you want to set the particle splatter multiplier too. The particle splatter mulitplier will change the amount of random offset all particles emitted in the effect will have.
	*/
	tfxAPI void SetEffectSplatterMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float splatter);

	/*
	Set the weight multiplier of an effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current weight of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param weight			A float of the amount that you want to set the particle weight multiplier too. The particle weight mulitplier will change the weight applied to particles in the effect at spawn time.
	*/
	tfxAPI void SetEffectWeightMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float weight);

	/*
	Set the overal scale of an effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current weight of the effect that was
	*						set in the TimelineFX editor.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param overal_scale	A float of the amount that you want to set the overal scale to. The overal scale is an simply way to change the size of an effect
	*/
	tfxAPI void SetEffectOveralScale(tfxParticleManager *pm, tfxEffectID effect_index, float overal_scale);

	/*
	Set the base noise offset for an effect
	* @param pm				A pointer to a tfxParticleManager where the effect is being managed.
	* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
	* @param noise_offset	A float of the amount that you want to set the effect noise offset to. By default when an effect is added to a particle manager a random noise offset will be set based on the Base Noise Offset Range property. Here you can override that
							value by setting it here. The most ideal time to set this would be immediately after you have added the effect to the particle manager, but you could call it any time you wanted for a constantly changing noise offset.
	*/
	tfxAPI void SetEffectBaseNoiseOffset(tfxParticleManager *pm, tfxEffectID effect_index, float noise_offset);

}

