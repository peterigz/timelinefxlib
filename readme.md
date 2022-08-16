# Timeline FX C++ library

This is an early alpha version of the TimelineFX library that you can use to incorporate 2D particle effects directly into your games and applications. You can use the [TimelineFX Particle Editor](https://www.rigzsoft.co.uk/timelinefx-alpha-version/) to design your effects.

**This is still very much a work in progress but feel free to use it and give some feedback on how to improve the interface and make it easier to use and implement.**

## Usage
The main library consists of 1 header file and 1 cpp file, there are no other dependencies. Options for including the library in your project would be one of the following three:
* Copy the files (including the Libraries folder) into your project and `#include "timelinefx.h"`.
* Add the timelinefx folder to your list of Includes Directories and `#include "timelinefx.h"`.
* Add TimelineFX.lib as a dependency and `#include "timelinefx.h"`

The library is render agnostic, so you will have to provide your own means to render the particles and to load the images that particles use to render.

To load an effect library that you create with the editor, use 

	EffectLibrary library;
	//Parameters are: filename, library object to load into, function pointer to the 
	//shape loader and any pointer to some custom data that you might want to pass through 
	//to the function.
	LoadEffectLibrary("MyEffects.zip", library, ShapeLoader, this);

### Setting up a ShapeLoader function to load the particle images
ShapeLoader is a pointer to your function that the library loader will use to load the images in the library. Each renderer will have its own way of doing that so you will have to adopt the appropriate way.

Here's an example of a Shapeloader using the same render that the Editor uses.

The ShapeLoader function expects the following parameters to be passed to it:
- `const char* filename`			*this will be the filename of the image being loaded from the library. You don't have to do anything with this if you don't need to.*
- `ImageData	&image_data`		*A struct containing data about the image. You will have to set image_data.ptr to point to the texture in your renderer for later use in the Render function that you will create to render the particles*
- `void *raw_image_data`			*The raw data of the image which you can use to load the image into graphics memory*
- `int image_memory_size`			*The size in bytes of the raw_image_data*
- `void *custom_data`				*This allows you to pass through an object you can use to access whatever is necessary to load the image into graphics memory, depending on the renderer that you're using*
```cpp
void ShapeLoader(const char* filename, tfxImageData &image_data, void *raw_image_data, int image_memory_size, void *custom_data) {
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

	//Get the texture where we're storing all the particle shapes
	QulkanTexture &texture = GetTexture(example->particle_textures);

	//You'll probably need to load the image in such a way depending on whether or not it's an animation or not
	if (image_data.animation_frames > 1) {
		//Add the spritesheet to the texture in our renderer
		u32 anim_index = texture.AddAnimation(image, (u32)image_data.image_size.x, (u32)image_data.image_size.y, (u32)image_data.animation_frames);
		//Important step: you need to point the ImageData.ptr to the appropriate handle in the renderer to point to the texture of the particle shape
		//You'll need to use this in your render function to tell your renderer which texture to use to draw the particle
		image_data.ptr = &texture.GetImage(anim_index);
	}
	else {
		//Add the image to the texture in our renderer
		u32 image_index = texture.AddImage(image);
		//Important step: you need to point the ImageData.ptr to the appropriate handle in the renderer to point to the texture of the particle shape
		//You'll need to use this in your render function to tell your renderer which texture to use to draw the particle
		image_data.ptr = &texture.GetImage(image_index);
	}
}
```

### Creating a render function to render the particles
Here's an example of a render function that you will need to write in order to integrate timeline fx with your specific renderer that you're using.
```cpp
//Here's an example of a render function that you will need to write in order to integrate timeline fx with your specific renderer that you're using
void TfxExample::RenderParticles(tfx::tfxParticleManager &pm, float tween) {
	//Renderer specific, get the layer that we will draw on (there's only one layer in this example)
	QulkanLayer &render_layer = GetLayer();

	for (int i = 0; i != tfxLAYERS; ++i) {
		for (auto &p : *pm.GetParticleBuffer(i)) {
			//Get the parent emitter of the particle to access some if it's properties
			tfxEffectEmitter &e = *p.parent;
			tfxEmitterProperties &properties = e.GetProperties();

			//Set the color for the quad. Note that the alpha value should dictate the blendfactor of the partciles, so 0 would be additive blend, 1 would be alpha blend.
			//Any value in between would be a mix of the 2.
			render_layer.SetColor(p.data.color.r, p.data.color.g, p.data.color.b, p.data.color.a);
			//Set the value that the color will be multiplied by, this happens in your fragment shader. You can always omit this if you're not using intensity
			render_layer.SetMultiplyFactor(p.data.intensity);
			//You can use render tweening to smooth particle movement from frame to frame by interpolating between captured and world states
			tfxVec3 tweened = Tween(tween, p.data.world_position, p.data.captured_position);
			//Get the frame of animation if it's an animated particle
			uint32_t frame = uint32_t(p.data.image_frame);
			//Set the image handle. It will differ from renderer to renderer how you access the right frame of animation. Here the pointer always points to the first frame, and then 
			//we can just add the current frame to the pointer to get the correct frame
			SetImageHandle(*(static_cast<QulkanImage*>(properties.image->ptr) + frame), e.current.image_handle.x, e.current.image_handle.y);
			//Add the particle frame of animation quad to the renderer for the next render pass at the particle position/rotation/scale
			render_layer.DrawSprite(
				*(static_cast<QulkanImage*>(properties.image->ptr) + frame),	//pointer to the image texture to use
				tweened.x, tweened.y,											//location of the particle
				p.data.world_rotations.roll,									//Rotation of the particle
				p.data.scale.x, p.data.scale.y);								//Scale of the particle
		}
	}
}
```

### Adding and updating effects
Use a tfxParticleManager to update effects and particles. You can create multiple particle managers if needed.

```cpp
	//Make sure that timelinefx udates effects at the same frequency as your update loop.
	tfx::SetUpdateFrequency(60);

	//Initialise a particle manager. This manages effects, emitters and the particles that they spawn
	//Depending on your needs you can use as many particle managers as you need.
	pm.Init();

	//Renderer specific
	QulkanLayer &layer = GetLayer();
	layer.SetBlendMode(Blendmode_alpha);
	//Ideally your renderer can use premultiplied alpha, this allows for particles to be drawn in a single draw call
	layer.SetBlendType(BlendType_pre_multiply);

	//In order to set update callbacks in your effects so that you can udpate them in realtime, prepare an effect template first:
	library.PrepareEffectTemplate("Torch", torch);
	//Then set the udpate callbacks
	torch.SetEffectUpdateCallback(UpdateTorchEffect);
	torch.SetEmitterUpdateCallback("Torch/Embers", UpdateTorchEmbers);
	torch.SetParticleUpdateCallback("Torch/Embers", UpdateEmberParticles);

	//Add the effect to the effect pool, the id you can use to access the effect is stored in torch_effect_id
	pm.AddEffect(torch);
```

Your callbacks need to pass in either an EffectEmitter object reference for effect/emitter callbacks or a Particle reference for a particle callback:

```cpp
//Change the position of the effect each update
void UpdateTorchEffect(tfx::tfxEffectEmitter &effect, tfxParentSpawnControls &spawn_controls) {
	Position(effect, tfxVec2(fMouseX(), fMouseY()));
}

//This callback is for a specific emitter in the effect, and changes quantity of particles spawned depending on the mouse location
void UpdateTorchEmbers(tfx::tfxEffectEmitter &emitter, tfxEmitterSpawnControls &spawn_controls) {
	emitter.current.qty = fMouseX();
}

//This callback gets called by every particle spawned by the Embers emitter in the torch effect and changes the color of the particles based on the location of the mouse pointer
void UpdateEmberParticles(tfx::tfxParticleData &particle_data, void *user_data) {
	float x = (float)particle_data.world_position.x / ScreenWidth();
	float y = (float)particle_data.world_position.y / ScreenHeight();
	particle_data.color.r = uint8_t(x * 255.f);
	particle_data.color.g = uint8_t(y * 255.f);
	particle_data.color.b = uint8_t(x * 255.f);
}

```

### Update/Render loop
You will need to update and render the particles in your main loop, and you will also need to make sure that you set the update frequency to the appropriate fps to match it to your own update loop frequencey.

	tfx::SetUpdateFrequency(60);

An update loop might look something like:

```cpp
	//Application specific update loop
	timer->Accumulate();

	while (timer->DoUpdate()) {

		//Update the particle manager
		pm.Update();

		timer->UnAccumulate();
	}

	timer->Set();

	//Call your renderer function to render all the particles
	RenderParticles(pm, timer->Tween());
```

See the files in the examples folder for a full example.

*More documentation to follow soon, but check out the header file and the Particle, tfxEffectEmitter and tfxParticleManager structs for more insights into usage*

