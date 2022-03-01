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
void TfxExample::RenderParticles(float tween) {
	//Renderer specific, get the layer that we will draw on (there's only one layer in this example)
	QulkanLayer &render_layer = GetLayer();

	//Loop through all the draw layers - particles can be assigned to a specific draw layer so you can draw them in a specific order if necessary. The layer is set in the editor on the properties tab.
	for (int layer = 0; layer != tfxLAYERS; ++layer) {
		//Use GetParticleBuffer(layer) to get all of the particles in the current layer
		for (auto p : *pm.GetParticleBuffer(layer)) {
			//In order to set the correct blendmode we need to get the property from the parent emitter that emitted the particle
			//A pointer to the parent emitter is stored in the parent member
			tfx::EffectEmitter &e = *p.parent;

			//Set the correct blendmode, see timelinefx::BlendMode. You may have to map the blendmodes depending on the renderer you use
			render_layer.SetBlendMode(qulkan::BlendMode(e.properties.blend_mode));
			//Set the color for the quad
			render_layer.SetColor(p.color.r, p.color.g, p.color.b, p.color.a);
			//Set the value that the color will be multiplied by, this happens in your fragment shader. You can always omit this if you're not using intensity
			render_layer.SetMultiplyFactor(p.intensity);
			//You can use render tweening to smooth particle movement from frame to frame by interpolating between captured and world states
			tfx::FormState tweened = tfx::Tween(tween, p.world, p.captured);
			//Is the particle using an image with more than one frame of animation?
			if (e.properties.image->animation_frames == 1 && e.properties.start_frame == 0) {
				//One frame of animation
				//Set the image handle, this the offset that the particle is drawn at.
				//This is where you can make use of the image->ptr from the ShapeLoader function, cast it into the appropriate type for the renderer
				qulkan::SetImageHandle(*static_cast<qulkan::QulkanImage*>(e.properties.image->ptr), e.current.image_handle.x, e.current.image_handle.y);
				//Add the particle image quad to the renderer for the next render pass at the particle position/rotation/scale
				render_layer.DrawSprite(*static_cast<qulkan::QulkanImage*>(e.properties.image->ptr), tweened.position.x, tweened.position.y, tweened.rotation, tweened.scale.x, tweened.scale.y);
			}
			else {
				//Multiple frames of animation
				//Get the current frame of animation that the particle is using
				uint32_t frame = uint32_t(p.image_frame);
				//frame must be within the bounds of the animation
				assert(frame >= 0 && frame < e.properties.image->animation_frames);

				//Set the image handle. It will differ from renderer to renderer how you access the right frame of animation. Here the pointer always points to the first frame, and then 
				//we can just add the current frame to the pointer to get the correct frame
				SetImageHandle(*(static_cast<QulkanImage*>(e.properties.image->ptr) + frame), e.current.image_handle.x, e.current.image_handle.y);
				//Add the particle frame of animation quad to the renderer for the next render pass at the particle position/rotation/scale
				render_layer.DrawSprite(*(static_cast<QulkanImage*>(e.properties.image->ptr) + frame), tweened.position.x, tweened.position.y, tweened.rotation, tweened.scale.x, tweened.scale.y);
			}
		}
	}
}
```

### Adding and updating effects
Effects and particles are updated using a `ParticleManager`. Every particle manager must be initialised before use to create the effect and particle pools. Set the maximum amount of effects and particles that you think you'll need for each manager. The default is 1000 effects and 50000 particles.

	tfx::ParticleManager pm;
	pm.Init();

You can have as many particle managers as you want if needed. The simplest way to add effects is straight from the library:

	tfx::AddEffect(pm, *library.GetEffect("Explosion 1"), MouseX(), MouseY());

This would be a simple "fire and forget" method and you would have no way of actually controlling the effect after it's added to the particle manager. In order to do that you need to create a template effect and then add that to the manager. Template effects allow you to set up different callbacks that you can use to make changes to the effect and emitters of the effect in realtime. Even individual particles can use a callback as well.

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
	//Set a callback for any particle spawned by a specific emitter
	torch.SetParticleUpdateCallback("Torch/Embers", UpdateEmbers);
	//Finally add the effect template to the particle manager
	tfx::AddEffect(pm, torch); 
```

Your callbacks need to pass in either an EffectEmitter object reference for effect/emitter callbacks or a Particle reference for a particle callback:

```cpp
//This one uses a callback to set the position of the effect to the current location of the mouse pointer
void UpdateTorchEffect(tfx::EffectEmitter &effect) {
	effect.Position(MouseX(), MouseY());
}

//This callback is for a specific emitter in the effect, and changes the amount of particles that are 
//spawned based on the location of the mouse pointer
void UpdateTorchFlames(tfx::EffectEmitter &emitter) {
	emitter.OverrideBaseAmount(MouseX());
}

//This callback gets called by every particle spawned by the Embers emitter in the torch effect and 
//changes the color of the particles based on the location of the mouse pointer
void UpdateEmbers(tfx::Particle &particle) {
	float x = (float)particle.world.position.x / ScreenWidth();
	float y = (float)particle.world.position.y / ScreenHeight();
	particle.color.r = uint8_t(x * 255.f);
	particle.color.g = uint8_t(y * 255.f);
	particle.color.b = uint8_t(x * 255.f);

	//You could also access the parent emitters user data with particle->parent->user_data if 
	//you needed to pull in more data to update the particle with
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

		//Inside your update loop, you will need to call the particle manager's update function
		pm.Update();

		timer->UnAccumulate();
	}

	timer->Set();

	//Call your renderer function to render all the particles
	RenderParticles(timer->Tween());
```

See the files in the examples folder for a full example.

*More documentation to follow soon, but check out the header file and the Particle, EffectEmitter and ParticleManager structs for more insights into usage*

