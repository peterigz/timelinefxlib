#pragma once

#include "pch.h"

//Example implementation of Timelinefx. This won't compile as it requires Qulkan 2D renderer which currently isn't available.
//It should still be a useful reference point for implementing with your own render library of choice.

using namespace qulkan;

class TfxExample {

public:

	tfx::EffectLibrary library;
	tfx::ParticleManager pm;
	tfx::EffectEmitterTemplate torch;
	qulkan::Timer *timer;
	qulkan::QulkanLayer *render_layer;
	qulkan::QulkanTextureLibrary *particle_textures;

	void Init();
	void Update(float ellapsed);
	void RenderParticles(float tween);
};

void ShapeLoader(const char* filename, tfx::ImageData &image_data, void *raw_image_data, int image_memory_size, void *custom_data);
void UpdateTorchEffect(tfx::EffectEmitter &effect);
void UpdateTorchFlames(tfx::EffectEmitter &effect);
void UpdateEmbers(tfx::Particle &particle);
