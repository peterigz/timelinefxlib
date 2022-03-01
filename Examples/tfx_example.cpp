//Example of implementing TimelineFX with a renderer.
//The TimelineFX API is render agnostic, so you will have to implement your own:
//1) ShapeLoader function to load each image in the library into the renderer texture buffer
//2) Renderer function which renders each particle
//Both functions are quite small and should be easy to implement

#include "pch.h"
#include "tfx_example.h"

using namespace qulkan;
using namespace tfx;

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

//You can utilise callbacks in effects and emitters to change the behaviour of effects however you need
//This one uses a callback to set the position of the effect to the current location of the mouse pointer
void UpdateTorchEffect(tfx::EffectEmitter &effect) {
	effect.Position((float)MouseX(), (float)MouseY());
}

//This callback is for a specific emitter in the effect, and changes the amount of particles that are spawned based on the location of the mouse pointer
void UpdateTorchFlames(tfx::EffectEmitter &emitter) {
	emitter.OverrideBaseAmount((float)MouseX());
}

//This callback gets called by every particle spawned by the Embers emitter in the torch effect and changes the color of the particles based on the location of the mouse pointer
void UpdateEmbers(tfx::Particle &particle) {
	float x = (float)particle.world.position.x / ScreenWidth();
	float y = (float)particle.world.position.y / ScreenHeight();
	particle.color.r = uint8_t(x * 255.f);
	particle.color.g = uint8_t(y * 255.f);
	particle.color.b = uint8_t(x * 255.f);

	//You could also access the parent emitters user data with particle->parent->user_data if you needed to pull in more data to update the particle with
}

void TfxExample::Init() {
	//Renderer specific - initialise the texture
	particle_textures = qulkan::CreateTexture("Particle Textures");
	QulkanTexture &texture = GetTexture(particle_textures);
	texture.SetUseFiltering(false);

	int shape_count = GetShapesInPackage("LibraryExamples.tfx");
	texture.ReserveImages(shape_count);
	//Load the effects library and pass the shape loader function pointer that you created earlier. Also pass this pointer to point to this object to give the shapeloader access to the texture we're loading the particle images into
	LoadEffectLibraryPackage("LibraryExamples.tfx", library, ShapeLoader, this);
	//Renderer specific
	texture.ProcessImages();

	//Application specific, set up a timer for the update loop
	timer = qulkan::Timer::Instance();
	timer->SetUpdateFrequency(60);

	//Make sure that timelinefx udates effects at the same frequency as your update loop.
	tfx::SetUpdateFrequency(60);

	//Renderer specific
	QulkanLayer &layer = GetLayer();
	layer.SetBlendMode(Blendmode_alpha);
	//Ideally your renderer can use premultiplied alpha, this allows for particles to be drawn in a single draw call
	layer.SetBlendType(BlendType_pre_multiply);

	//You must initialise the particle manager, and set the maximum amount of effects and particles that you're likely to be updating per frame
	//You'll need a higher effects amount if you're using a lot of effects with sub effects
	pm.Init();

	//In order to set update callbacks in your effects so that you can udpate them in realtime, prepare an effect template first:
	library.PrepareEffectTemplate("Torch", torch);
	//You can set some custom user data which you can cast in the callback if needed. (useful if attaching the effect to an object in your game)
	torch.SetUserData(this);
	//Set the udpate callback for the effect
	torch.SetUpdateCallback(UpdateTorchEffect);
	//Set a callback for a specific emitter in the effect
	torch.SetUpdateCallback("Torch/Embers", UpdateTorchFlames);
	//Set a callback for any particle spawned by a specific emitter
	torch.SetParticleUpdateCallback("Torch/Embers", UpdateEmbers);
	//Finally add the effect to the particle manager
	tfx::AddEffect(pm, torch);
}

void TfxExample::Update(float ellapsed) {
	//Application specific update loop
	timer->Accumulate();

	//Renderer specific
	SetActiveRenderQueue();

	if (MouseHit(App.Mouse.LeftButton)) {
		//You don't have to prepare a template effect before adding to the effect library. You can just add it to the particle manager and forget about it.
		//Once the effect is added though, you won't be able to access it again to make any realtime changes, so use a template if you need to do that.
		AddEffect(pm, *library.GetEffect("Explosion 1"), (float)MouseX(), (float)MouseY());
	}

	if (KeyHit(GLFW_KEY_SPACE)) {
		paused = !paused;
	}

	while (timer->DoUpdate()) {

		//Inside your update loop, you will need to call the particle manager's update function
		if (!paused)
			pm.Update();

		timer->UnAccumulate();
	}

	timer->Set();

	//Call your renderer function to render all the particles
	RenderParticles(timer->Tween());
}

//Here's an example of a render function that you will need to write in order to integrate timeline fx with your specific renderer that you're using
//I should think that you could quite easily multi-thread this as well, as long as your renderer is happy with that
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

//Application specific, this just sets the function to call each render update
void UpdateTfxExample(float ellapsed, void *data) {
	TfxExample *example = static_cast<TfxExample*>(data);
	example->Update(ellapsed);
}

int main(int argc, char** argv) {
	//Application specific
	TfxExample example;

	QulkanCreateInfo info;
	info.title = "TimelineFX Example";
	info.maximised = true;
	info.show_fps_in_title = true;

	std::cout << sizeof(Particle) << std::endl;

	Initialise(info);

	example.Init();

	SetUserData(&example);
	SetUserCallback(UpdateTfxExample);

	try {
		Start();
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}