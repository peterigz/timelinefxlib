# Timeline FX C++ library

This is an early alpha version of the TimelineFX library that you can use to incorporate 2D and 3D particle effects directly into your games and applications. You can use the [TimelineFX Particle Editor](https://www.rigzsoft.co.uk/timelinefx-alpha-version/) to design your effects.

**This is still very much a work in progress but feel free to use it and give some feedback on how to improve the interface and make it easier to use and implement.**

## Usage
The main library consists of 1 header file and 1 cpp file, there are no other dependencies. Options for including the library in your project would be one of the following three:
* Copy the timelinefx.h and .cpp files into your project and `#include "timelinefx.h"`.
* Add the timelinefx folder to your list of Includes Directories and `#include "timelinefx.h"`.
* Add TimelineFX.lib as a dependency and `#include "timelinefx.h"`

The library is render agnostic, so you will have to provide your own means to render the particles and to load the images that particles use to render.

So far some of the features of the Library:
* Fast. Utilises multithreading and SIMD intrinsics to simulate particles as fast as possible. Work is always on going to make things faster.
* Can pre-bake effects so that they can be animated entirely in a compute shader for ultimate speed. The baked data can also be compressed to save up to 90% memory depending on the effect. The compute shader does the job of interpolating the sprite data to keep things smooth.
* Render both 2d or 3d particles.
* Render agnostic, just a few simple integrations to get it drawing in your own renderer. (See shaders folder for example shaders)
* Simple API which is easy to use.
* Runs on Intel/ARM CPUs, currently tested on Windows and Mac. Linux coming soon.
* Many different types of effects are possible with more coming as features are added.

## Basic Overview of Usage

To load an effect library that you create with the editor, use 

	EffectLibrary library;
	//Parameters are: filename, library object to load into, function pointer to the 
	//shape loader and any pointer to some custom data that you might want to pass through 
	//to the function.
	LoadEffectLibrary("MyEffects.zip", &library, ShapeLoader, &my_custom_data);

### Setting up a ShapeLoader function to load the particle images
ShapeLoader is a pointer to your function that the library loader will use to load the images in the library. Each renderer will have its own way of doing that so you will have to adopt the appropriate way.

Here's an example of a Shapeloader using the same renderer that the Editor uses. This is taken from an example you can [https://github.com/peterigz/zest/blob/3685db7c7e066f43e25db5015851c13af5ef89e0/examples/GLFW/zest-timelinefx/zest-timelinefx.cpp#L54](find here): 

```
/* Before you load an effects file, you will need to define a ShapeLoader function that passes the following parameters: */
//const char* filename			- this will be the filename of the image being loaded from the library. You don't have to do anything with this if you don't need to.
//ImageData	&image_data			- A struct containing data about the image. You will have to set image_data.ptr to point to the texture in your renderer for later use in the Render function that you will create to render the particles
//void *raw_image_data			- The raw data of the image which you can use to load the image into graphics memory
//int image_memory_size			- The size in bytes of the raw_image_data
//void *custom_data				- This allows you to pass through an object you can use to access whatever is necessary to load the image into graphics memory, depending on the renderer that you're using
void ShapeLoader(const char* filename, tfx_image_data_t *image_data, void *raw_image_data, int image_memory_size, void *custom_data) {
	//Cast your custom data, this can be anything you want
	VadersGame *game = static_cast<VadersGame*>(custom_data);

	//This shape loader example uses the STB image library to load the raw bitmap (png usually) data
	zest_bitmap_t bitmap = zest_NewBitmap();
	zest_LoadBitmapImageMemory(&bitmap, (unsigned char*)raw_image_data, image_memory_size, 0);
	//Convert the image to RGBA which is necessary for this particular renderer
	zest_ConvertBitmapToRGBA(&bitmap, 255);
	//The editor has the option to convert an bitmap to an alpha map. I will probably change this so that it gets baked into the saved effect so you won't need to apply the filter here.
	//Alpha map is where all color channels are set to 255
	if (image_data->import_filter)
		zest_ConvertBitmapToAlpha(&bitmap);

	//Get the texture where we're storing all the particle shapes
	//You'll probably need to load the image in such a way depending on whether or not it's an animation or not
	if (image_data->animation_frames > 1) {
		//Add the spritesheet to the texture in our renderer
		float max_radius = 0;
		image_data->ptr = zest_AddTextureAnimationBitmap(game->particle_texture, &bitmap, (u32)image_data->image_size.x, (u32)image_data->image_size.y, (u32)image_data->animation_frames, &max_radius, 1);
		//Important step: you need to point the ImageData.ptr to the appropriate handle in the renderer to point to the texture of the particle shape
		//You'll need to use this in your render function to tell your renderer which texture to use to draw the particle
	}
	else {
		//Add the image to the texture in our renderer
		image_data->ptr = zest_AddTextureImageBitmap(game->particle_texture, &bitmap);
		//Important step: you need to point the ImageData.ptr to the appropriate handle in the renderer to point to the texture of the particle shape
		//You'll need to use this in your render function to tell your renderer which texture to use to draw the particle
	}
}
```

### Creating a render function to render the particles
Here's an example of a render function that you will need to write in order to integrate timeline fx with your specific renderer that you're using. a 2d render function would be very similar to the 3d one below.
This is taken from an example you can This is taken from an example you can [https://github.com/peterigz/zest/blob/3685db7c7e066f43e25db5015851c13af5ef89e0/examples/GLFW/zest-timelinefx/zest-timelinefx.cpp#L54](find here):
```cpp
void RenderParticles3d(tfx_particle_manager_t& pm, VadersGame* game) {
	//Let our renderer know that we want to draw to the billboard layer.
	zest_SetBillboardDrawing(game->billboard_layer, game->particle_texture, game->particle_descriptor, game->billboard_pipeline);
	//Cycle through each layer
	//There is also a macro :tfxEachLayer which you could use like so:
	//for(tfxEachLayer) {
	//and that will output exactly the same code as the below line
	for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
		//tfx_sprite_soa_t is a struct of arrays so first get a pointer to it from the particle manager
		tfx_sprite_soa_t *sprites = SpritesInLayer(&pm, layer);
		//Now loop over all the sprites in the layer
		for (int i = 0; i != SpritesInLayerCount(&pm, layer); ++i) {
			//Set the color to draw the sprite using the sprite color. Note that bacuase this is a struct of arrays we access each sprite element using an array lookup
			zest_SetLayerColor(game->billboard_layer, sprites->color[i].r, sprites->color[i].g, sprites->color[i].b, sprites->color[i].a);
			//Set the instensity of the sprite.
			zest_SetLayerIntensity(game->billboard_layer, sprites->intensity[i]);
			//Grab the image pointer from the emitter properties in the library
			zest_image image = (zest_image)GetSpriteImagePointer(&pm, sprites->property_indexes[i]);
			//Grab the sprite handle which is used to offset the position of the sprite
			tfx_vec2_t handle = GetSpriteHandle(&pm, sprites->property_indexes[i]);
			//Render specific function to draw a billboard to the screen using the values in the sprite data
			zest_DrawBillboard(game->billboard_layer, image + ((sprites->property_indexes[i] & 0x00FF0000) >> 16),
				&sprites->transform_3d[i].position.x,
				sprites->alignment[i],
				&sprites->transform_3d[i].rotations.x,
				&handle.x,
				sprites->stretch[i],
				(sprites->property_indexes[i] & 0xFF000000) >> 24,
				sprites->transform_3d[i].scale.x, sprites->transform_3d[i].scale.y);
		}
	}
}
```

### Using a particle manager
All effects are updated using a particle manager so you'll need to set one up like so:

```cpp
tfx_particle_manager_t pm;

//We can define the maximum amount of particles for each layer in an array and pass that into the initialisation function
tfxU32 layer_max_values[tfxLAYERS];
layer_max_values[0] = 5000;
layer_max_values[1] = 2500;
layer_max_values[2] = 2500;
layer_max_values[3] = 2500;
InitParticleManagerFor3d(&pm, &library, layer_max_values, 1000, tfxParticleManagerMode_unordered, true, false, 2048);
/*
Parameters are:
InitParticleManagerFor3d(	[pointer to the particle manager], 
							[pointer to the library], 
							[array of maximum particles per layer], 
							[maximum number of emitters allowed], 
							[the mode of the particle manager: unordered, age ordered, depth orderd], 
							[true/false use double buffered sprites for interpolating positions for smoother updates], 
							[true/false dynamically increase particle limits], 
							[batch size for multithreading]);
*/
```

### Effect templates
To add effects to the particle manager, first you can grab them from the library and put them into an effect template. Doing this means you can make alterations without changing the base effect in the library.

```cpp
tfx_effect_template_t effect_template;
PrepareEffectTemplate(&library, "Star Burst Flash", &effect_template);
```

You can then use the effect template to add to the particle manager:

```cpp
//Each time you add an effect to the particle manager it generates an ID which you can use to modify the effect whilst it's being updated
tfxEffectID effect_id;
//Add the effect template to the particle manager
if(AddEffectToParticleManager(&pm, &effect_template, &effect_id)) {
	//Set the effect position using the effect id
	SetEffectPosition(&game->pm, effect_id, {10.f, 0.f, 0.f});
}
```

### Updating the Particle Manager in your game loop
//You might have a fixed stepped time loop something like this
```cpp
	zest_TimerAccumulate(game->timer);
	while (zest_TimerDoUpdate(game->timer)) {

		//Update the particle manager with the length of time that elapsed since the last frame. So you could also put this outside of the timer loop
		//if you're not using one, as long as the elapsed time since the last frame is passed into the function
		UpdateParticleManager(&pm, FrameLength);

		zest_TimerUnAccumulate(game->timer);
	}
	zest_TimerSet(game->timer);
```

Then after the particle manager is updated you can just render the particles calling the render fuction as exampled above:

```cpp
	//Render the particles with our custom render function
	RenderParticles3d(game->pm, game);
```

For more in depth examples, currently you can check out a simple [space invaders game here](https://github.com/peterigz/zest/blob/main/examples/GLFW/zest-vaders/zest-vaders.cpp) and [another minimal example here](https://github.com/peterigz/zest/blob/main/examples/GLFW/zest-timelinefx/zest-timelinefx.cpp)

All of the API functions are documented so check out the [header file to view them](https://github.com/peterigz/timelinefxlib/blob/a1dcce20b608c228004b466de4a4fb0df59bb609/timelinefx.h#L6202).
