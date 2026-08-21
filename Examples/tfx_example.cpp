//Example of implementing TimelineFX with a renderer.
//The TimelineFX API is render agnostic, so you will have to implement your own:
//1) ShapeLoader function to load each image in the library into the renderer texture buffer
//2) Renderer function which renders each particle
//Both functions are quite small and should be easy to implement depending on the render library you're using

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
}

void UpdateTorchEffect(tfx::tfxEffectEmitter &effect, tfxParentSpawnControls &spawn_controls) {
	Position(effect, tfxVec2(fMouseX(), fMouseY()));
}

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

void TfxExample::Update(float ellapsed) {
	//Application specific update loop
	timer->Accumulate();

	//Renderer specific
	SetActiveRenderQueue();

	if (KeyHit(GLFW_KEY_SPACE)) {
		paused = !paused;
	}

	while (timer->DoUpdate()) {

		//Update the particle manager
		pm.Update();

		timer->UnAccumulate();
	}

	timer->Set();

	//Call your renderer function to render all the particles
	RenderParticles(pm, timer->Tween());

	tfxStr32 pcount = "Particles: ";
	pcount.Appendf("%i", pm.ParticleCount());
	GetLayer().SetColor(255, 255, 255, 255);
	GetLayer().SetMultiplyFactor(1.f);
	GetLayer().Text(pcount.c_str(), 20.f, 50.f);
}

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

	std::cout << sizeof(tfxParticle) << std::endl;

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