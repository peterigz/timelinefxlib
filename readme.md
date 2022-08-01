# Timeline FX C++ library

**Quick note to say that currently this library is out of date with the latest addition of 3d in the editor. Development is on the 3d branch which is what I'm working on next so watch this space!**

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
void TfxExample::RenderEffect(tfx::tfxEffect &effect, float tween) {
	//Renderer specific, get the layer that we will draw on (there's only one layer in this example)
	QulkanLayer &render_layer = GetLayer();

	for (EachLayer) {
		for (auto &s : effect.sprites[layer]) {
			//In order to set the correct blendmode we need to get the property from the parent emitter that emitted the particle
			//A pointer to the parent emitter is stored in the parent member
			//Set the correct blendmode, see timelinefx::BlendMode. You may have to map the blendmodes depending on the renderer you use
			render_layer.SetBlendMode(qulkan::BlendMode((s.parameters & 0xF0000000) >> 28));
			//Set the color for the quad
			render_layer.SetColor(s.color.r, s.color.g, s.color.b, s.color.a);
			//Set the value that the color will be multiplied by, this happens in your fragment shader. You can always omit this if you're not using intensity
			render_layer.SetMultiplyFactor(s.intensity);
			//You can use render tweening to smooth particle movement from frame to frame by interpolating between captured and world states
			tfx::FormState tweened = tfx::Tween(tween, s.world, s.captured);
			//Is the particle using an image with more than one frame of animation?
			//Get the current frame of animation that the particle is using
			uint32_t frame = s.parameters & 0x00FFFFFF;
			//Set the image handle. It will differ from renderer to renderer how you access the right frame of animation. Here the pointer always points to the first frame, and then 
			//we can just add the current frame to the pointer to get the correct frame
			SetImageHandle(*(static_cast<QulkanImage*>(s.ptr) + frame), s.handle.x, s.handle.y);
			//Add the particle frame of animation quad to the renderer for the next render pass at the particle position/rotation/scale
			render_layer.DrawSprite(*(static_cast<QulkanImage*>(s.ptr) + frame), tweened.position.x, tweened.position.y, tweened.rotation, tweened.scale.x, tweened.scale.y);
		}
	}
}
```

### Adding and updating effects
The memory for effects, emitters and particles are managed using a tfxEffectPool. You must call Init on the effect pool before use.

	tfx::tfxEffectPool effect_pool;
	effect_pool.Init();

You can have as many effect pools as you want if needed. The simplest way to add effects is straight from the library:

	tfx::tfxEffectID effect_id;
	tfx::PoolEffect(effect_pool, *library.GetEffect("Explosion 1"), effect_id);

The id for the effect will be stored in effect_id in this case and you can use that with GetEffect to retrieve the effect later to move it or render it etc.

```cpp
	//In order to set update callbacks in your effects so that you can udpate them in realtime, 
	//prepare an effect template first:
	tfx::EffectEmitterTemplate torch;
	library.PrepareEffectTemplate("Torch", torch);	//Torch is the name of the effect in the library.
	//You can set some custom user data which you can cast in the callback if needed. (useful if 
	//attaching the effect to an object in your game)
	torch.SetUserData(this);
	//Set the udpate callback for the effect
	torch.SetUpdateCallback(UpdateTorchEffect);
	//Set a callback for a specific emitter in the effect
	torch.SetUpdateCallback("Torch/Embers", UpdateTorchFlames);
	//Finally add the effect template to the effect pool
	tfx::PoolEffectFromTemplate(effect_pool, torch, effect_id); 
```

Your callbacks need to pass in either an EffectEmitter object reference for effect/emitter callbacks or a Particle reference for a particle callback:

```cpp
//This one uses a callback to set the opacity of the effect to something else
void UpdateTorchEffect(tfx::tfxEffect &effect) {
	effect.current.opacity(0.5f);
}

//This callback is for a specific emitter in the effect, and changes the angle of the emitter
//based on the location of the mouse pointer
void UpdateTorchFlames(tfx::tfxEmitter &emitter) {
	emitter.(fMouseX() *0.01f);
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

		//You can get an effect from the pool with GetEffect, which you might do to position the effect
		GetEffect(effect_pool, torch_effect_id).Position(fMouseX(), fMouseY());
		//In order to update the effects each frame you'll need to call UpdateEffect on each Effect you're using
		UpdateEffect(effect_pool, torch_effect_id);

		timer->UnAccumulate();
	}

	timer->Set();

	//Call your renderer function to render all the particles
	RenderParticles(GetEffect(effect_pool, torch_effect_id), timer->Tween());
```

See the files in the examples folder for a full example.

*More documentation to follow soon, but check out the header file and the Particle, EffectEmitter and ParticleManager structs for more insights into usage*

