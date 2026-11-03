//TimelineFX example: GLFW + WebGPU through wgpu-native. The window, instance, adapter, device, surface and main loop, see ../README.md

#ifndef GLFW_INCLUDE_NONE
#define GLFW_INCLUDE_NONE
#endif
#include <GLFW/glfw3.h>
#if defined(_WIN32)
#define GLFW_EXPOSE_NATIVE_WIN32
#elif defined(__linux__)
#define GLFW_EXPOSE_NATIVE_X11
#define GLFW_EXPOSE_NATIVE_WAYLAND
#endif
#include <GLFW/glfw3native.h>
#include <webgpu/webgpu.h>
#include <webgpu/wgpu.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "timelinefx.h"
#include "tfx_impl_common.h"
#include "tfx_impl_webgpu.h"
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
#define TFX_EXAMPLE_SHADER_FOLDER "webgpu_shaders"
#endif
#define EXAMPLE_FEATURES (example_scene_feature_ribbons | example_scene_feature_gpu_simulation | example_scene_feature_baked)

#define EXAMPLE_WINDOW_WIDTH 1280
#define EXAMPLE_WINDOW_HEIGHT 800
#define EXAMPLE_MEMORY_POOL_SIZE (64ull * 1024ull * 1024ull)
#define EXAMPLE_DEPTH_FORMAT WGPUTextureFormat_Depth32Float

typedef struct example_s {
	GLFWwindow *window;
	WGPUInstance instance;
	WGPUSurface surface;
	WGPUAdapter adapter;
	WGPUDevice device;
	WGPUQueue queue;
	WGPUTextureFormat color_format;
	uint32_t width;					//The size the surface was last configured at
	uint32_t height;
	WGPUTexture depth_texture;
	WGPUTextureView depth_view;
	bool surface_stale;
	const char *backend_name;

	tfx_impl_webgpu_t renderer;
	example_scene_t scene;
} example_t;

//Filled by the adapter and device callbacks
typedef struct example_request_s {
	bool done;
	WGPUAdapter adapter;
	WGPUDevice device;
} example_request_t;

static example_t example;

static void ErrorCallback(int error, const char *description) {
	printf("GLFW error %i: %s\n", error, description);
}

static void KeyCallback(GLFWwindow *window, int key, int scancode, int action, int mods) {
	(void)scancode; (void)mods;
	if (action != GLFW_PRESS) {
		return;
	}
	if (key == GLFW_KEY_ESCAPE) {
		glfwSetWindowShouldClose(window, GLFW_TRUE);
	} else if (key == GLFW_KEY_SPACE) {
		example_scene_Burst(&example.scene);
	} else if (key == GLFW_KEY_O) {
		example_scene_ToggleOrbit(&example.scene);
	} else if (key == GLFW_KEY_L) {
		example_scene_ReloadSecondLibrary(&example.scene);
	}
}

static void ScrollCallback(GLFWwindow *window, double x_offset, double y_offset) {
	(void)window; (void)x_offset;
	example_scene_Zoom(&example.scene, (float)y_offset);
}

static WGPUStringView ExampleString(const char *text) {
	WGPUStringView view = { text, WGPU_STRLEN };
	return view;
}

//wgpu-native's own log, which carries the Vulkan validation layer's and the D3D12 debug layer's messages
static void LogCallback(WGPULogLevel level, WGPUStringView message, void *userdata) {
	(void)userdata;
	printf("wgpu %s: %.*s\n", level == WGPULogLevel_Error ? "error" : "warning", (int)message.length, message.data ? message.data : "");
}

//Every WebGPU validation error that no error scope catches comes here
static void UncapturedErrorCallback(WGPUDevice const *device, WGPUErrorType type, WGPUStringView message, void *userdata1, void *userdata2) {
	(void)device; (void)userdata1; (void)userdata2;
	printf("WebGPU %s error: %.*s\n", type == WGPUErrorType_Validation ? "validation" : type == WGPUErrorType_OutOfMemory ? "out of memory" : "internal",
		(int)message.length, message.data ? message.data : "");
}

static void DeviceLostCallback(WGPUDevice const *device, WGPUDeviceLostReason reason, WGPUStringView message, void *userdata1, void *userdata2) {
	(void)device; (void)userdata1; (void)userdata2;
	if (reason != WGPUDeviceLostReason_Destroyed && reason != WGPUDeviceLostReason_CallbackCancelled) {
		printf("The WebGPU device was lost: %.*s\n", (int)message.length, message.data ? message.data : "");
	}
}

static void RequestAdapterCallback(WGPURequestAdapterStatus status, WGPUAdapter adapter, WGPUStringView message, void *userdata1, void *userdata2) {
	(void)userdata2;
	example_request_t *request = (example_request_t *)userdata1;
	request->done = true;
	if (status == WGPURequestAdapterStatus_Success) {
		request->adapter = adapter;
	} else {
		printf("No WebGPU adapter: %.*s\n", (int)message.length, message.data ? message.data : "");
	}
}

static void RequestDeviceCallback(WGPURequestDeviceStatus status, WGPUDevice device, WGPUStringView message, void *userdata1, void *userdata2) {
	(void)userdata2;
	example_request_t *request = (example_request_t *)userdata1;
	request->done = true;
	if (status == WGPURequestDeviceStatus_Success) {
		request->device = device;
	} else {
		printf("Could not create the WebGPU device: %.*s\n", (int)message.length, message.data ? message.data : "");
	}
}

//Validation turns on the Vulkan validation layer or the D3D12 debug layer under wgpu, which needs the Vulkan SDK or the Graphics Tools.
//WebGPU's own validation is always on
static bool CreateInstance(example_t *state, WGPUInstanceBackend backends, bool gpu_validation) {
	WGPUInstanceExtras extras = {};
	extras.chain.sType = (WGPUSType)WGPUSType_InstanceExtras;
	extras.backends = backends;
	extras.flags = WGPUInstanceFlag_Default;
#ifndef NDEBUG
	extras.flags = WGPUInstanceFlag_Debug | WGPUInstanceFlag_Validation;
#endif
	if (gpu_validation) {
		extras.flags |= WGPUInstanceFlag_GPUBasedValidation;
	}
	//DXC compiles the kernel in a fraction of FXC's time. The file name alone finds dxcompiler.dll on the PATH (the Vulkan SDK has one),
	//where an empty path asks for a DXC built into wgpu-native, which the release builds don't have
	extras.dx12ShaderCompiler = WGPUDx12Compiler_Dxc;
	extras.dxcPath = ExampleString("dxcompiler.dll");
	WGPUInstanceDescriptor descriptor = WGPU_INSTANCE_DESCRIPTOR_INIT;
	descriptor.nextInChain = &extras.chain;
	wgpuSetLogCallback(LogCallback, NULL);
	wgpuSetLogLevel(WGPULogLevel_Warn);
	state->instance = wgpuCreateInstance(&descriptor);
	if (!state->instance) {
		printf("Could not create the WebGPU instance\n");
	}
	return state->instance != NULL;
}

static bool CreateSurface(example_t *state) {
	WGPUSurfaceDescriptor descriptor = WGPU_SURFACE_DESCRIPTOR_INIT;
#if defined(_WIN32)
	WGPUSurfaceSourceWindowsHWND source = WGPU_SURFACE_SOURCE_WINDOWS_HWND_INIT;
	source.hinstance = GetModuleHandle(NULL);
	source.hwnd = glfwGetWin32Window(state->window);
	descriptor.nextInChain = &source.chain;
#elif defined(__linux__)
	WGPUSurfaceSourceWaylandSurface wayland_source = WGPU_SURFACE_SOURCE_WAYLAND_SURFACE_INIT;
	WGPUSurfaceSourceXlibWindow x11_source = WGPU_SURFACE_SOURCE_XLIB_WINDOW_INIT;
	if (glfwGetPlatform() == GLFW_PLATFORM_WAYLAND) {
		wayland_source.display = glfwGetWaylandDisplay();
		wayland_source.surface = glfwGetWaylandWindow(state->window);
		descriptor.nextInChain = &wayland_source.chain;
	} else {
		x11_source.display = glfwGetX11Display();
		x11_source.window = glfwGetX11Window(state->window);
		descriptor.nextInChain = &x11_source.chain;
	}
#else
	printf("This example makes WebGPU surfaces on Windows and Linux. Elsewhere give wgpuInstanceCreateSurface your platform's surface source\n");
	return false;
#endif
	state->surface = wgpuInstanceCreateSurface(state->instance, &descriptor);
	return state->surface != NULL;
}

//The adapter's limits on buffer sizes rather than WebGPU's defaults, which cap a storage buffer binding at 128 MB
static bool CreateDevice(example_t *state, WGPUBackendType backend_type) {
	example_request_t request = {};
	WGPURequestAdapterOptions adapter_options = WGPU_REQUEST_ADAPTER_OPTIONS_INIT;
	adapter_options.powerPreference = WGPUPowerPreference_HighPerformance;
	adapter_options.backendType = backend_type;
	adapter_options.compatibleSurface = state->surface;
	WGPURequestAdapterCallbackInfo adapter_callback = WGPU_REQUEST_ADAPTER_CALLBACK_INFO_INIT;
	adapter_callback.mode = WGPUCallbackMode_AllowProcessEvents;
	adapter_callback.callback = RequestAdapterCallback;
	adapter_callback.userdata1 = &request;
	wgpuInstanceRequestAdapter(state->instance, &adapter_options, adapter_callback);
	while (!request.done) {
		wgpuInstanceProcessEvents(state->instance);
	}
	state->adapter = request.adapter;
	if (!state->adapter) {
		return false;
	}
	WGPUAdapterInfo info = WGPU_ADAPTER_INFO_INIT;
	wgpuAdapterGetInfo(state->adapter, &info);
	state->backend_name = info.backendType == WGPUBackendType_D3D12 ? "D3D12" : info.backendType == WGPUBackendType_Vulkan ? "Vulkan"
		: info.backendType == WGPUBackendType_Metal ? "Metal" : "other";
	printf("%.*s on %s, wgpu-native %u.%u.%u.%u\n", (int)info.device.length, info.device.data ? info.device.data : "", state->backend_name,
		wgpuGetVersion() >> 24, (wgpuGetVersion() >> 16) & 0xFF, (wgpuGetVersion() >> 8) & 0xFF, wgpuGetVersion() & 0xFF);
	wgpuAdapterInfoFreeMembers(info);

	WGPULimits adapter_limits = WGPU_LIMITS_INIT;
	wgpuAdapterGetLimits(state->adapter, &adapter_limits);
	WGPULimits required_limits = WGPU_LIMITS_INIT;
	required_limits.maxStorageBufferBindingSize = adapter_limits.maxStorageBufferBindingSize;
	required_limits.maxBufferSize = adapter_limits.maxBufferSize;
	WGPUDeviceDescriptor device_descriptor = WGPU_DEVICE_DESCRIPTOR_INIT;
	device_descriptor.label = ExampleString("TimelineFX example");
	device_descriptor.requiredLimits = &required_limits;
	device_descriptor.uncapturedErrorCallbackInfo.callback = UncapturedErrorCallback;
	device_descriptor.deviceLostCallbackInfo.mode = WGPUCallbackMode_AllowSpontaneous;
	device_descriptor.deviceLostCallbackInfo.callback = DeviceLostCallback;
	WGPURequestDeviceCallbackInfo device_callback = WGPU_REQUEST_DEVICE_CALLBACK_INFO_INIT;
	device_callback.mode = WGPUCallbackMode_AllowProcessEvents;
	device_callback.callback = RequestDeviceCallback;
	device_callback.userdata1 = &request;
	request.done = false;
	wgpuAdapterRequestDevice(state->adapter, &device_descriptor, device_callback);
	while (!request.done) {
		wgpuInstanceProcessEvents(state->instance);
	}
	state->device = request.device;
	if (!state->device) {
		return false;
	}
	state->queue = wgpuDeviceGetQueue(state->device);

	//A non-sRGB format to match the other examples, whose shaders write the colors as they are
	WGPUSurfaceCapabilities capabilities = WGPU_SURFACE_CAPABILITIES_INIT;
	wgpuSurfaceGetCapabilities(state->surface, state->adapter, &capabilities);
	state->color_format = capabilities.formatCount > 0 ? capabilities.formats[0] : WGPUTextureFormat_BGRA8Unorm;
	for (size_t index = 0; index != capabilities.formatCount; ++index) {
		if (capabilities.formats[index] == WGPUTextureFormat_BGRA8Unorm || capabilities.formats[index] == WGPUTextureFormat_RGBA8Unorm) {
			state->color_format = capabilities.formats[index];
			break;
		}
	}
	wgpuSurfaceCapabilitiesFreeMembers(capabilities);
	return true;
}

static void ReleaseDepth(example_t *state) {
	if (state->depth_view) wgpuTextureViewRelease(state->depth_view);
	if (state->depth_texture) wgpuTextureRelease(state->depth_texture);
	state->depth_view = NULL;
	state->depth_texture = NULL;
}

static bool ConfigureSurface(example_t *state, uint32_t width, uint32_t height) {
	WGPUSurfaceConfiguration configuration = WGPU_SURFACE_CONFIGURATION_INIT;
	configuration.device = state->device;
	configuration.format = state->color_format;
	configuration.usage = WGPUTextureUsage_RenderAttachment;
	configuration.width = width;
	configuration.height = height;
	configuration.presentMode = WGPUPresentMode_Fifo;
	wgpuSurfaceConfigure(state->surface, &configuration);
	ReleaseDepth(state);
	WGPUTextureDescriptor depth_descriptor = WGPU_TEXTURE_DESCRIPTOR_INIT;
	depth_descriptor.usage = WGPUTextureUsage_RenderAttachment;
	depth_descriptor.size.width = width;
	depth_descriptor.size.height = height;
	depth_descriptor.format = EXAMPLE_DEPTH_FORMAT;
	state->depth_texture = wgpuDeviceCreateTexture(state->device, &depth_descriptor);
	state->depth_view = state->depth_texture ? wgpuTextureCreateView(state->depth_texture, NULL) : NULL;
	state->width = width;
	state->height = height;
	state->surface_stale = false;
	return state->depth_view != NULL;
}

//Every atlas after the load, and the second library's after a reload. The queue copies them, so the pixels can go at once
static bool UploadPendingAtlases(example_t *state) {
	bool uploaded = true;
	tfx_impl_atlas_t *atlas = example_scene_GetPendingAtlas(&state->scene);
	while (atlas) {
		uploaded = tfx_impl_webgpu_UploadAtlas(&state->renderer, atlas) && uploaded;
		tfx_impl_FreeAtlasPixels(atlas);
		atlas = example_scene_GetPendingAtlas(&state->scene);
	}
	return uploaded;
}

static void UpdateWindowTitle(example_t *state, unsigned int frames_per_second) {
	tfx_version_t version = tfx_GetVersion();
	char title[256];
	tfx_gpu_particle_stats_t gpu_stats = tfx_GetGPUParticleStats(state->scene.stage);
	snprintf(title, sizeof(title), "TimelineFX WebGPU (%s) Example - %s %i.%i - %u fps - %u CPU particles - %u GPU particles - %llu dropped GPU spawns",
		state->backend_name, version.name, version.major, version.minor, frames_per_second, tfx_GetParticleCount(state->scene.stage), gpu_stats.estimated_live_particles,
		(unsigned long long)gpu_stats.dropped_spawns);
	glfwSetWindowTitle(state->window, title);
}

static void RenderFrame(example_t *state, WGPUCommandEncoder encoder, WGPUTextureView target, const tfx_impl_uniforms_t *uniforms) {
	WGPURenderPassColorAttachment color_attachment = WGPU_RENDER_PASS_COLOR_ATTACHMENT_INIT;
	color_attachment.view = target;
	color_attachment.loadOp = WGPULoadOp_Clear;
	color_attachment.storeOp = WGPUStoreOp_Store;
	color_attachment.clearValue.r = 0.0;
	color_attachment.clearValue.g = 0.1;
	color_attachment.clearValue.b = 0.2;
	color_attachment.clearValue.a = 1.0;
	WGPURenderPassDepthStencilAttachment depth_attachment = WGPU_RENDER_PASS_DEPTH_STENCIL_ATTACHMENT_INIT;
	depth_attachment.view = state->depth_view;
	depth_attachment.depthLoadOp = WGPULoadOp_Clear;
	depth_attachment.depthStoreOp = WGPUStoreOp_Discard;
	depth_attachment.depthClearValue = 1.f;
	WGPURenderPassDescriptor pass_descriptor = WGPU_RENDER_PASS_DESCRIPTOR_INIT;
	pass_descriptor.colorAttachmentCount = 1;
	pass_descriptor.colorAttachments = &color_attachment;
	pass_descriptor.depthStencilAttachment = &depth_attachment;
	WGPURenderPassEncoder render_pass = wgpuCommandEncoderBeginRenderPass(encoder, &pass_descriptor);
	tfx_impl_webgpu_Render(&state->renderer, render_pass, uniforms);
	wgpuRenderPassEncoderEnd(render_pass);
	wgpuRenderPassEncoderRelease(render_pass);
}

static void ReleaseDevice(example_t *state) {
	ReleaseDepth(state);
	if (state->surface) {
		if (state->width) wgpuSurfaceUnconfigure(state->surface);
		wgpuSurfaceRelease(state->surface);
	}
	if (state->queue) wgpuQueueRelease(state->queue);
	if (state->device) wgpuDeviceRelease(state->device);
	if (state->adapter) wgpuAdapterRelease(state->adapter);
	if (state->instance) wgpuInstanceRelease(state->instance);
}

//Usage: example_glfw_webgpu [updates per second] [--frames N] [--test-debug-output] [--gpu-validation] [--vulkan | --d3d12]
int main(int argc, char **argv) {
	double updates_per_second = 0.0;
	long frame_limit = 0;
	bool test_debug_output = false;
	bool gpu_validation = false;
	WGPUInstanceBackend backends = WGPUInstanceBackend_Primary;
	WGPUBackendType backend_type = WGPUBackendType_Undefined;
	for (int index = 1; index < argc; ++index) {
		if (strcmp(argv[index], "--gpu-validation") == 0) {
			gpu_validation = true;
		} else if (strcmp(argv[index], "--vulkan") == 0) {
			backends = WGPUInstanceBackend_Vulkan;
			backend_type = WGPUBackendType_Vulkan;
		} else if (strcmp(argv[index], "--d3d12") == 0) {
			backends = WGPUInstanceBackend_DX12;
			backend_type = WGPUBackendType_D3D12;
		} else if (strcmp(argv[index], "--frames") == 0 && index + 1 < argc) {
			frame_limit = atol(argv[++index]);
		} else if (strcmp(argv[index], "--test-debug-output") == 0) {
			test_debug_output = true;
		} else if (atof(argv[index]) > 0.0) {
			updates_per_second = atof(argv[index]);
		}
	}

	glfwSetErrorCallback(ErrorCallback);
	if (!glfwInit()) {
		return 1;
	}
	glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
	example.window = glfwCreateWindow(EXAMPLE_WINDOW_WIDTH, EXAMPLE_WINDOW_HEIGHT, "TimelineFX WebGPU Example", NULL, NULL);
	int width = 0;
	int height = 0;
	if (example.window) {
		glfwSetKeyCallback(example.window, KeyCallback);
		glfwSetScrollCallback(example.window, ScrollCallback);
		glfwGetFramebufferSize(example.window, &width, &height);
	}
	bool ready = example.window && CreateInstance(&example, backends, gpu_validation) && CreateSurface(&example) && CreateDevice(&example, backend_type)
		&& ConfigureSurface(&example, (uint32_t)width, (uint32_t)height);

	//Starts the worker threads and the memory pool that everything in the library allocates from
	tfx_BeginTimelineFX(tfx_GetDefaultThreadCount(), EXAMPLE_MEMORY_POOL_SIZE);
	ready = ready && tfx_impl_webgpu_Init(&example.renderer, example.device, example.color_format, EXAMPLE_DEPTH_FORMAT, TFX_EXAMPLE_SHADER_FOLDER)
		&& example_scene_Load(&example.scene, EXAMPLE_LIBRARY_PATH, EXAMPLE_SECOND_LIBRARY_PATH, EXAMPLE_SPRITE_DATA_PATH, EXAMPLE_FEATURES)
		&& tfx_impl_webgpu_CreateAtlas(&example.renderer, example.scene.atlas_layer_size, example.scene.atlas_layer_count)
		&& UploadPendingAtlases(&example)
		&& (!example.scene.animation_manager || tfx_impl_webgpu_UploadSpriteData(&example.renderer, example.scene.animation_manager, example.scene.baked_shapes));
	if (ready) {
		//A buffer that can't be both mapped for reading and a storage buffer, to prove the errors reach the console before trusting their silence
		if (test_debug_output) {
			WGPUBufferDescriptor invalid = WGPU_BUFFER_DESCRIPTOR_INIT;
			invalid.usage = WGPUBufferUsage_MapRead | WGPUBufferUsage_Storage;
			invalid.size = 16;
			WGPUBuffer buffer = wgpuDeviceCreateBuffer(example.device, &invalid);
			if (buffer) {
				wgpuBufferRelease(buffer);
			}
		}
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

	double last_time = glfwGetTime();
	double title_timer = 0.0;
	unsigned int frame_count = 0;
	long frames_run = 0;

	bool running = ready;
	while (running && !glfwWindowShouldClose(example.window)) {
		glfwPollEvents();
		glfwGetFramebufferSize(example.window, &width, &height);
		if (width == 0 || height == 0) {
			glfwWaitEvents();
			continue;
		}
		if (example.surface_stale || (uint32_t)width != example.width || (uint32_t)height != example.height) {
			if (!ConfigureSurface(&example, (uint32_t)width, (uint32_t)height)) {
				break;
			}
		}
		double time = glfwGetTime();
		double elapsed_milliseconds = (time - last_time) * 1000.0;
		last_time = time;

		tfx_impl_frame_t frame;
		tfx_impl_uniforms_t uniforms;
		example_scene_Update(&example.scene, elapsed_milliseconds, (float)example.width, (float)example.height, tfx_impl_clip_space_depth_zero_to_one, &frame, &uniforms);
		UploadPendingAtlases(&example);
		tfx_impl_webgpu_UploadFrame(&example.renderer, &frame);

		//The compute runs whether or not there's a surface texture to draw into, because every GPU particle fetch has to be dispatched
		WGPUCommandEncoder encoder = wgpuDeviceCreateCommandEncoder(example.device, NULL);
		tfx_impl_webgpu_Compute(&example.renderer, encoder, &uniforms);
		WGPUSurfaceTexture surface_texture = WGPU_SURFACE_TEXTURE_INIT;
		wgpuSurfaceGetCurrentTexture(example.surface, &surface_texture);
		bool acquired = surface_texture.status == WGPUSurfaceGetCurrentTextureStatus_SuccessOptimal || surface_texture.status == WGPUSurfaceGetCurrentTextureStatus_SuccessSuboptimal;
		WGPUTextureView target = acquired ? wgpuTextureCreateView(surface_texture.texture, NULL) : NULL;
		if (target) {
			RenderFrame(&example, encoder, target, &uniforms);
		}
		example.surface_stale = surface_texture.status != WGPUSurfaceGetCurrentTextureStatus_SuccessOptimal;
		WGPUCommandBuffer commands = wgpuCommandEncoderFinish(encoder, NULL);
		wgpuQueueSubmit(example.queue, 1, &commands);
		wgpuCommandBufferRelease(commands);
		wgpuCommandEncoderRelease(encoder);
		if (target) {
			wgpuSurfacePresent(example.surface);
			wgpuTextureViewRelease(target);
		}
		if (surface_texture.texture) {
			wgpuTextureRelease(surface_texture.texture);
		}

		frame_count++;
		title_timer += elapsed_milliseconds;
		if (title_timer >= 1000.0) {
			UpdateWindowTitle(&example, frame_count);
			title_timer -= 1000.0;
			frame_count = 0;
		}
		if (frame_limit > 0 && ++frames_run >= frame_limit) {
			running = false;
		}
	}

	if (example.device) {
		wgpuDevicePoll(example.device, true, NULL);
	}
	example_scene_Shutdown(&example.scene);
	tfx_EndTimelineFX();
	tfx_impl_webgpu_Shutdown(&example.renderer);

	ReleaseDevice(&example);
	if (example.window) {
		glfwDestroyWindow(example.window);
	}
	glfwTerminate();
	return 0;
}
