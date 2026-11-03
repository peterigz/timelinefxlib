//TimelineFX example: SDL3 with the SDL_gpu API, on Vulkan or Direct3D 12. The window, device and main loop, see ../README.md

#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
//Only for the feature struct SDL_GPUVulkanOptions takes. The Vulkan SDK's own include folder carries older SDL3 headers, so these come from SDL's source
#define VK_NO_PROTOTYPES
#include <vulkan/vulkan_core.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "timelinefx.h"
#include "tfx_impl_common.h"
#include "tfx_impl_sdlgpu.h"
#include "example_scene.h"

//TFX_EXAMPLE_ROOT is set by the CMake project to the implementations folder
#ifndef TFX_EXAMPLE_ROOT
#define TFX_EXAMPLE_ROOT "."
#endif
#define EXAMPLE_LIBRARY_PATH TFX_EXAMPLE_ROOT "/example.tfx"
#define EXAMPLE_SECOND_LIBRARY_PATH TFX_EXAMPLE_ROOT "/example_second.tfx"
#define EXAMPLE_SPRITE_DATA_PATH TFX_EXAMPLE_ROOT "/example.tfxsd"
//TFX_EXAMPLE_SHADER_FOLDER is where the CMake project has slangc build the shaders
#ifndef TFX_EXAMPLE_SHADER_FOLDER
#define TFX_EXAMPLE_SHADER_FOLDER "sdl_shaders"
#endif
#define EXAMPLE_FEATURES (example_scene_feature_ribbons | example_scene_feature_gpu_simulation | example_scene_feature_baked)

#define EXAMPLE_WINDOW_WIDTH 1280
#define EXAMPLE_WINDOW_HEIGHT 800
#define EXAMPLE_MEMORY_POOL_SIZE (64ull * 1024ull * 1024ull)
#define EXAMPLE_DEPTH_FORMAT SDL_GPU_TEXTUREFORMAT_D32_FLOAT

typedef struct example_s {
	SDL_Window *window;
	SDL_GPUDevice *device;
	SDL_GPUTexture *depth_texture;
	Uint32 depth_width;
	Uint32 depth_height;
	tfx_impl_sdlgpu_t renderer;
	example_scene_t scene;
} example_t;

static example_t example;

static void ExampleLog(void *, int, SDL_LogPriority, const char *message) {
	printf("SDL: %s\n", message);
}

//Vulkan 1.2 for spirv_1_5 and shaderDrawParameters for SV_VertexID. SDL only turns features on that are asked for here
static SDL_GPUDevice *CreateDevice(bool direct3d12) {
	SDL_PropertiesID properties = SDL_CreateProperties();
	SDL_SetBooleanProperty(properties, SDL_PROP_GPU_DEVICE_CREATE_DEBUGMODE_BOOLEAN, true);
	SDL_SetBooleanProperty(properties, SDL_PROP_GPU_DEVICE_CREATE_SHADERS_SPIRV_BOOLEAN, true);
	SDL_SetBooleanProperty(properties, SDL_PROP_GPU_DEVICE_CREATE_SHADERS_DXIL_BOOLEAN, true);
	SDL_SetStringProperty(properties, SDL_PROP_GPU_DEVICE_CREATE_NAME_STRING, direct3d12 ? "direct3d12" : "vulkan");
	VkPhysicalDeviceVulkan11Features vulkan_11_features = {};
	vulkan_11_features.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_1_FEATURES;
	vulkan_11_features.shaderDrawParameters = VK_TRUE;
	SDL_GPUVulkanOptions vulkan_options = {};
	vulkan_options.vulkan_api_version = VK_MAKE_API_VERSION(0, 1, 2, 0);
	vulkan_options.feature_list = &vulkan_11_features;
	SDL_SetPointerProperty(properties, SDL_PROP_GPU_DEVICE_CREATE_VULKAN_OPTIONS_POINTER, &vulkan_options);
	//SDL turns the validation layer on in debug mode but installs no messenger, so the layer is told to print to stdout. Prove it reaches you with a deliberate error
	SDL_setenv_unsafe("VK_KHRONOS_VALIDATION_DEBUG_ACTION", "VK_DBG_LAYER_ACTION_LOG_MSG", 0);
	SDL_GPUDevice *device = SDL_CreateGPUDeviceWithProperties(properties);
	SDL_DestroyProperties(properties);
	if (!device) {
		printf("Could not create a %s GPU device: %s\n", direct3d12 ? "Direct3D 12" : "Vulkan", SDL_GetError());
	}
	return device;
}

//Matches the swapchain, which can change size with the window
static bool EnsureDepthTexture(example_t *state, Uint32 width, Uint32 height) {
	if (state->depth_texture && state->depth_width == width && state->depth_height == height) {
		return true;
	}
	if (state->depth_texture) {
		SDL_ReleaseGPUTexture(state->device, state->depth_texture);
	}
	SDL_GPUTextureCreateInfo create_info = {};
	create_info.type = SDL_GPU_TEXTURETYPE_2D;
	create_info.format = EXAMPLE_DEPTH_FORMAT;
	create_info.usage = SDL_GPU_TEXTUREUSAGE_DEPTH_STENCIL_TARGET;
	create_info.width = width;
	create_info.height = height;
	create_info.layer_count_or_depth = 1;
	create_info.num_levels = 1;
	//D3D12 warns on every clear that doesn't match the value the texture was made with
	create_info.props = SDL_CreateProperties();
	SDL_SetFloatProperty(create_info.props, SDL_PROP_GPU_TEXTURE_CREATE_D3D12_CLEAR_DEPTH_FLOAT, 1.f);
	state->depth_texture = SDL_CreateGPUTexture(state->device, &create_info);
	SDL_DestroyProperties(create_info.props);
	if (!state->depth_texture) {
		printf("Could not create the depth texture: %s\n", SDL_GetError());
		return false;
	}
	state->depth_width = width;
	state->depth_height = height;
	return true;
}

//Every atlas after the load, and the second library's after a reload
static bool UploadPendingAtlases(example_t *state, SDL_GPUCommandBuffer *command_buffer) {
	bool uploaded = true;
	tfx_impl_atlas_t *atlas = example_scene_GetPendingAtlas(&state->scene);
	while (atlas) {
		uploaded = tfx_impl_sdlgpu_UploadAtlas(&state->renderer, command_buffer, atlas) && uploaded;
		tfx_impl_FreeAtlasPixels(atlas);
		atlas = example_scene_GetPendingAtlas(&state->scene);
	}
	return uploaded;
}

static bool LoadScene(example_t *state) {
	if (!example_scene_Load(&state->scene, EXAMPLE_LIBRARY_PATH, EXAMPLE_SECOND_LIBRARY_PATH, EXAMPLE_SPRITE_DATA_PATH, EXAMPLE_FEATURES)
		|| !tfx_impl_sdlgpu_CreateAtlas(&state->renderer, state->scene.atlas_layer_size, state->scene.atlas_layer_count)) {
		return false;
	}
	SDL_GPUCommandBuffer *command_buffer = SDL_AcquireGPUCommandBuffer(state->device);
	if (!command_buffer) {
		return false;
	}
	bool uploaded = UploadPendingAtlases(state, command_buffer)
		&& (!state->scene.animation_manager || tfx_impl_sdlgpu_UploadSpriteData(&state->renderer, command_buffer, state->scene.animation_manager, state->scene.baked_shapes));
	return SDL_SubmitGPUCommandBuffer(command_buffer) && uploaded;
}

static void UpdateWindowTitle(example_t *state, Uint32 frames_per_second, const char *driver) {
	tfx_version_t version = tfx_GetVersion();
	char title[256];
	tfx_gpu_particle_stats_t gpu_stats = tfx_GetGPUParticleStats(state->scene.stage);
	snprintf(title, sizeof(title), "TimelineFX SDL3 GPU (%s) Example - %s %i.%i - %u fps - %u CPU particles - %u GPU particles - %llu dropped GPU spawns",
		driver, version.name, version.major, version.minor, frames_per_second, tfx_GetParticleCount(state->scene.stage), gpu_stats.estimated_live_particles,
		(unsigned long long)gpu_stats.dropped_spawns);
	SDL_SetWindowTitle(state->window, title);
}

//Usage: example_sdl3_gpu [--d3d12] [updates per second]
int main(int argc, char **argv) {
	bool direct3d12 = false;
	double updates_per_second = 0.0;
	for (int index = 1; index < argc; ++index) {
		if (strcmp(argv[index], "--d3d12") == 0) {
			direct3d12 = true;
		} else if (atof(argv[index]) > 0.0) {
			updates_per_second = atof(argv[index]);
		}
	}

	SDL_SetLogOutputFunction(ExampleLog, NULL);
	if (!SDL_Init(SDL_INIT_VIDEO)) {
		printf("Could not start SDL: %s\n", SDL_GetError());
		return 1;
	}
	example.window = SDL_CreateWindow("TimelineFX SDL3 GPU Example", EXAMPLE_WINDOW_WIDTH, EXAMPLE_WINDOW_HEIGHT, SDL_WINDOW_RESIZABLE);
	example.device = example.window ? CreateDevice(direct3d12) : NULL;
	if (!example.device || !SDL_ClaimWindowForGPUDevice(example.device, example.window)) {
		printf("Could not set up the window for the GPU device: %s\n", SDL_GetError());
		if (example.device) {
			SDL_DestroyGPUDevice(example.device);
		}
		SDL_DestroyWindow(example.window);
		SDL_Quit();
		return 1;
	}
	const char *driver = SDL_GetGPUDeviceDriver(example.device);

	//Starts the worker threads and the memory pool that everything in the library allocates from
	tfx_BeginTimelineFX(tfx_GetDefaultThreadCount(), EXAMPLE_MEMORY_POOL_SIZE);
	bool ready = tfx_impl_sdlgpu_Init(&example.renderer, example.device, SDL_GetGPUSwapchainTextureFormat(example.device, example.window), EXAMPLE_DEPTH_FORMAT, TFX_EXAMPLE_SHADER_FOLDER)
		&& LoadScene(&example);
	if (ready) {
		//Before the start sizes the GPU rings, so they hold a slower rate's longer updates
		if (updates_per_second > 0.0) {
			tfx_SetGPUSizingUpdateRate((float)updates_per_second);
		}
		example_scene_Start(&example.scene);
		if (updates_per_second > 0.0) {
			example.scene.timer = tfx_impl_CreateTimer(updates_per_second);
		}
		printf("Space: burst  O: pause the camera  L: free and reload the second library  Mouse wheel: zoom  Escape: quit\n");
	}

	Uint64 counter_frequency = SDL_GetPerformanceFrequency();
	Uint64 last_counter = SDL_GetPerformanceCounter();
	double title_timer = 0.0;
	Uint32 frame_count = 0;

	bool running = ready;
	while (running) {
		SDL_Event event;
		while (SDL_PollEvent(&event)) {
			if (event.type == SDL_EVENT_QUIT) {
				running = false;
			} else if (event.type == SDL_EVENT_MOUSE_WHEEL) {
				example_scene_Zoom(&example.scene, event.wheel.y);
			} else if (event.type == SDL_EVENT_KEY_DOWN && !event.key.repeat) {
				if (event.key.key == SDLK_ESCAPE) {
					running = false;
				} else if (event.key.key == SDLK_SPACE) {
					example_scene_Burst(&example.scene);
				} else if (event.key.key == SDLK_O) {
					example_scene_ToggleOrbit(&example.scene);
				} else if (event.key.key == SDLK_L) {
					example_scene_ReloadSecondLibrary(&example.scene);
				}
			}
		}
		if (!running) {
			break;
		}

		Uint64 counter = SDL_GetPerformanceCounter();
		double elapsed_milliseconds = (double)(counter - last_counter) * 1000.0 / (double)counter_frequency;
		last_counter = counter;

		int window_width = 0;
		int window_height = 0;
		SDL_GetWindowSizeInPixels(example.window, &window_width, &window_height);
		SDL_GPUCommandBuffer *command_buffer = SDL_AcquireGPUCommandBuffer(example.device);
		if (!command_buffer) {
			printf("Could not acquire a command buffer: %s\n", SDL_GetError());
			break;
		}
		//SDL_gpu flips Vulkan's clip space to match D3D, so both take the D3D row of the clip space table
		tfx_impl_frame_t frame;
		tfx_impl_uniforms_t uniforms;
		example_scene_Update(&example.scene, elapsed_milliseconds, (float)window_width, (float)window_height, tfx_impl_clip_space_depth_zero_to_one, &frame, &uniforms);
		//Recorded whether or not there's a swapchain texture to draw to, because every GPU particle fetch has to be dispatched
		UploadPendingAtlases(&example, command_buffer);
		tfx_impl_sdlgpu_UploadFrame(&example.renderer, command_buffer, &frame);
		tfx_impl_sdlgpu_Compute(&example.renderer, command_buffer, &uniforms);

		SDL_GPUTexture *swapchain_texture = NULL;
		Uint32 swapchain_width = 0;
		Uint32 swapchain_height = 0;
		if (SDL_WaitAndAcquireGPUSwapchainTexture(command_buffer, example.window, &swapchain_texture, &swapchain_width, &swapchain_height)
			&& swapchain_texture && EnsureDepthTexture(&example, swapchain_width, swapchain_height)) {
			SDL_GPUColorTargetInfo color_target = {};
			color_target.texture = swapchain_texture;
			color_target.clear_color = { 0.f, 0.1f, 0.2f, 1.f };
			color_target.load_op = SDL_GPU_LOADOP_CLEAR;
			color_target.store_op = SDL_GPU_STOREOP_STORE;
			SDL_GPUDepthStencilTargetInfo depth_target = {};
			depth_target.texture = example.depth_texture;
			depth_target.clear_depth = 1.f;
			depth_target.load_op = SDL_GPU_LOADOP_CLEAR;
			depth_target.store_op = SDL_GPU_STOREOP_DONT_CARE;
			depth_target.stencil_load_op = SDL_GPU_LOADOP_DONT_CARE;
			depth_target.stencil_store_op = SDL_GPU_STOREOP_DONT_CARE;
			SDL_GPURenderPass *render_pass = SDL_BeginGPURenderPass(command_buffer, &color_target, 1, &depth_target);
			tfx_impl_sdlgpu_Render(&example.renderer, command_buffer, render_pass, &uniforms);
			SDL_EndGPURenderPass(render_pass);
		}
		SDL_SubmitGPUCommandBuffer(command_buffer);

		frame_count++;
		title_timer += elapsed_milliseconds;
		if (title_timer >= 1000.0) {
			UpdateWindowTitle(&example, frame_count, driver);
			title_timer -= 1000.0;
			frame_count = 0;
		}
	}

	SDL_WaitForGPUIdle(example.device);
	example_scene_Shutdown(&example.scene);
	tfx_EndTimelineFX();
	tfx_impl_sdlgpu_Shutdown(&example.renderer);

	if (example.depth_texture) {
		SDL_ReleaseGPUTexture(example.device, example.depth_texture);
	}
	SDL_ReleaseWindowFromGPUDevice(example.device, example.window);
	SDL_DestroyGPUDevice(example.device);
	SDL_DestroyWindow(example.window);
	SDL_Quit();
	return 0;
}
