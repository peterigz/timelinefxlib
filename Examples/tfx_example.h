#pragma once

#include "pch.h"

using namespace qulkan;

class TfxExample {

public:

	tfx::EffectLibrary library;
	tfx::ParticleManager pm;
	tfx::EffectEmitterTemplate torch;
	qulkan::Timer *timer;
	u32 render_layer;
	u32 base_target;
	u32 particle_textures;
	bool paused = false;

	void Init();
	void Update(float ellapsed);
	void RenderParticles(float tween);
};

void ShapeLoader(const char* filename, tfx::ImageData &image_data, void *raw_image_data, int image_memory_size, void *custom_data);
void UpdateTorchEffect(tfx::EffectEmitter &effect);
void UpdateTorchFlames(tfx::EffectEmitter &effect);
void UpdateEmbers(tfx::Particle &particle);
