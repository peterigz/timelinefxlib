#pragma once

#include "pch.h"

using namespace qulkan;

class TfxExample {

public:

	tfx::tfxEffectLibrary library;
	tfx::tfxEffectTemplate torch;
	tfx::tfxParticleManager pm;
	tfx::tfxEffectID torch_effect_id;
	qulkan::Timer *timer;
	u32 render_layer;
	u32 base_target;
	u32 particle_textures;
	bool paused = false;

	void Init();
	void Update(float ellapsed);
	void RenderParticles(tfx::tfxParticleManager &pm, float tween);
};

void ShapeLoader(const char* filename, tfx::tfxImageData &image_data, void *raw_image_data, int image_memory_size, void *custom_data);
void UpdateTorchEffect(tfx::tfxEffectEmitter &effect, tfx::tfxParentSpawnControls &spawn_controls);
void UpdateTorchEmbers(tfx::tfxEffectEmitter &emitter, tfx::tfxEmitterSpawnControls &spawn_controls);
void UpdateEmberParticles(tfx::tfxParticleData &particle, void *user_data);
