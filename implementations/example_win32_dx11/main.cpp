//TimelineFX example: Win32 + Direct3D 11. The window, device and main loop, see ../README.md

#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <d3d11.h>
#include <d3d11sdklayers.h>
#include <stdio.h>
#include <stdlib.h>
#include "timelinefx.h"
#include "tfx_impl_common.h"
#include "tfx_impl_dx11.h"
#include "example_scene.h"

//TFX_EXAMPLE_ROOT is set by the CMake project to the implementations folder
#ifndef TFX_EXAMPLE_ROOT
#define TFX_EXAMPLE_ROOT "."
#endif
#define EXAMPLE_SHADER_PATH TFX_EXAMPLE_ROOT "/shaders/tfx_dx11.hlsl"
#define EXAMPLE_LIBRARY_PATH TFX_EXAMPLE_ROOT "/example.tfx"
#define EXAMPLE_SPRITE_DATA_PATH TFX_EXAMPLE_ROOT "/example.tfxsd"
//TFX_EXAMPLE_KERNEL_PATH is the GPU particle kernel the CMake project builds with slangc, when it finds slangc
#ifdef TFX_EXAMPLE_KERNEL_PATH
#define EXAMPLE_KERNEL_PATH TFX_EXAMPLE_KERNEL_PATH
#define EXAMPLE_FEATURES (example_scene_feature_ribbons | example_scene_feature_gpu_simulation | example_scene_feature_baked)
#else
#define EXAMPLE_KERNEL_PATH NULL
#define EXAMPLE_FEATURES (example_scene_feature_ribbons | example_scene_feature_baked)
#endif

#define EXAMPLE_WINDOW_WIDTH 1280
#define EXAMPLE_WINDOW_HEIGHT 800
#define EXAMPLE_MEMORY_POOL_SIZE (64ull * 1024ull * 1024ull)

typedef struct example_s {
	HWND window;
	ID3D11Device *device;
	ID3D11DeviceContext *device_context;
	IDXGISwapChain *swap_chain;
	ID3D11RenderTargetView *render_target_view;
	ID3D11Texture2D *depth_texture;
	ID3D11DepthStencilView *depth_stencil_view;
	UINT width;
	UINT height;
	UINT resize_width;				//Set by WM_SIZE, applied at the start of the next frame
	UINT resize_height;
	bool swap_chain_occluded;
	ID3D11InfoQueue *info_queue;	//Debug layer messages, null without the debug layer

	tfx_impl_dx11_t renderer;
	example_scene_t scene;
} example_t;

static example_t example;

static void ReleaseRenderTargets(example_t *state) {
	if (state->depth_stencil_view) { state->depth_stencil_view->Release(); state->depth_stencil_view = NULL; }
	if (state->depth_texture) { state->depth_texture->Release(); state->depth_texture = NULL; }
	if (state->render_target_view) { state->render_target_view->Release(); state->render_target_view = NULL; }
}

static bool CreateRenderTargets(example_t *state) {
	ID3D11Texture2D *back_buffer = NULL;
	if (FAILED(state->swap_chain->GetBuffer(0, IID_PPV_ARGS(&back_buffer)))) {
		printf("Could not get the swap chain back buffer\n");
		return false;
	}
	D3D11_TEXTURE2D_DESC back_buffer_description;
	back_buffer->GetDesc(&back_buffer_description);
	HRESULT result = state->device->CreateRenderTargetView(back_buffer, NULL, &state->render_target_view);
	back_buffer->Release();
	if (FAILED(result)) {
		printf("Could not create the render target view (0x%08lx)\n", result);
		return false;
	}
	state->width = back_buffer_description.Width;
	state->height = back_buffer_description.Height;

	D3D11_TEXTURE2D_DESC depth_description = {};
	depth_description.Width = state->width;
	depth_description.Height = state->height;
	depth_description.MipLevels = 1;
	depth_description.ArraySize = 1;
	depth_description.Format = DXGI_FORMAT_D32_FLOAT;
	depth_description.SampleDesc.Count = 1;
	depth_description.Usage = D3D11_USAGE_DEFAULT;
	depth_description.BindFlags = D3D11_BIND_DEPTH_STENCIL;
	result = state->device->CreateTexture2D(&depth_description, NULL, &state->depth_texture);
	if (FAILED(result)) {
		printf("Could not create the depth buffer (0x%08lx)\n", result);
		return false;
	}
	result = state->device->CreateDepthStencilView(state->depth_texture, NULL, &state->depth_stencil_view);
	if (FAILED(result)) {
		printf("Could not create the depth stencil view (0x%08lx)\n", result);
		return false;
	}
	return true;
}

static bool CreateDevice(example_t *state) {
	DXGI_SWAP_CHAIN_DESC swap_chain_description = {};
	swap_chain_description.BufferCount = 2;
	swap_chain_description.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
	swap_chain_description.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
	swap_chain_description.OutputWindow = state->window;
	swap_chain_description.SampleDesc.Count = 1;
	swap_chain_description.Windowed = TRUE;
	swap_chain_description.SwapEffect = DXGI_SWAP_EFFECT_FLIP_DISCARD;

	UINT device_flags = 0;
#ifdef _DEBUG
	device_flags |= D3D11_CREATE_DEVICE_DEBUG;
#endif
	//Feature level 11.0 is the minimum: compute shaders, raw buffers and BC7 shapes all need it
	const D3D_FEATURE_LEVEL feature_levels[] = { D3D_FEATURE_LEVEL_11_1, D3D_FEATURE_LEVEL_11_0 };
	D3D_FEATURE_LEVEL feature_level;
	HRESULT result = D3D11CreateDeviceAndSwapChain(NULL, D3D_DRIVER_TYPE_HARDWARE, NULL, device_flags, feature_levels, 2, D3D11_SDK_VERSION,
		&swap_chain_description, &state->swap_chain, &state->device, &feature_level, &state->device_context);
	if (result == DXGI_ERROR_SDK_COMPONENT_MISSING && (device_flags & D3D11_CREATE_DEVICE_DEBUG)) {
		printf("The D3D11 debug layer isn't installed (Graphics Tools optional feature), continuing without it\n");
		device_flags &= ~D3D11_CREATE_DEVICE_DEBUG;
		result = D3D11CreateDeviceAndSwapChain(NULL, D3D_DRIVER_TYPE_HARDWARE, NULL, device_flags, feature_levels, 2, D3D11_SDK_VERSION,
			&swap_chain_description, &state->swap_chain, &state->device, &feature_level, &state->device_context);
	}
	if (FAILED(result)) {
		printf("Could not create a D3D11 device with feature level 11.0 or higher (0x%08lx)\n", result);
		return false;
	}
	if (device_flags & D3D11_CREATE_DEVICE_DEBUG) {
		state->device->QueryInterface(IID_PPV_ARGS(&state->info_queue));
	}
	return CreateRenderTargets(state);
}

//The debug layer writes to the debugger output, this repeats it on the console
static void PrintDebugLayerMessages(example_t *state) {
	if (!state->info_queue) {
		return;
	}
	UINT64 message_count = state->info_queue->GetNumStoredMessages();
	for (UINT64 message_index = 0; message_index != message_count; ++message_index) {
		SIZE_T message_size = 0;
		state->info_queue->GetMessage(message_index, NULL, &message_size);
		D3D11_MESSAGE *message = (D3D11_MESSAGE *)malloc(message_size);
		if (message && SUCCEEDED(state->info_queue->GetMessage(message_index, message, &message_size))) {
			printf("D3D11: %.*s\n", (int)message->DescriptionByteLength, message->pDescription);
		}
		free(message);
	}
	state->info_queue->ClearStoredMessages();
}

static void ReleaseDevice(example_t *state) {
	ReleaseRenderTargets(state);
	if (state->info_queue) { state->info_queue->Release(); state->info_queue = NULL; }
	if (state->swap_chain) { state->swap_chain->Release(); state->swap_chain = NULL; }
	if (state->device_context) { state->device_context->Release(); state->device_context = NULL; }
	if (state->device) { state->device->Release(); state->device = NULL; }
}

static LRESULT WINAPI WindowProcedure(HWND window, UINT message, WPARAM wparam, LPARAM lparam) {
	switch (message) {
	case WM_SIZE:
		if (wparam != SIZE_MINIMIZED) {
			example.resize_width = (UINT)LOWORD(lparam);
			example.resize_height = (UINT)HIWORD(lparam);
		}
		return 0;
	case WM_MOUSEWHEEL:
		example_scene_Zoom(&example.scene, (float)GET_WHEEL_DELTA_WPARAM(wparam) / WHEEL_DELTA);
		return 0;
	case WM_KEYDOWN:
		if (wparam == VK_ESCAPE) {
			DestroyWindow(window);
		} else if (wparam == VK_SPACE) {
			example_scene_Burst(&example.scene);
		} else if (wparam == 'O') {
			example_scene_ToggleOrbit(&example.scene);
		}
		return 0;
	case WM_SYSCOMMAND:
		//Stops alt from opening the window menu and pausing the loop
		if ((wparam & 0xfff0) == SC_KEYMENU) {
			return 0;
		}
		break;
	case WM_DESTROY:
		PostQuitMessage(0);
		return 0;
	}
	return DefWindowProcA(window, message, wparam, lparam);
}

static void UpdateWindowTitle(example_t *state, UINT frames_per_second) {
	tfx_version_t version = tfx_GetVersion();
	char title[256];
	tfx_gpu_particle_stats_t gpu_stats = tfx_GetGPUParticleStats(state->scene.stage);
	snprintf(title, sizeof(title), "TimelineFX D3D11 Example - %s %i.%i - %u fps - %u CPU particles - %u GPU particles - %llu dropped GPU spawns",
		version.name, version.major, version.minor, frames_per_second, tfx_GetParticleCount(state->scene.stage), gpu_stats.live_particles,
		(unsigned long long)gpu_stats.dropped_spawns);
	SetWindowTextA(state->window, title);
}

//Usage: example_win32_dx11 [updates per second]
int main(int argc, char **argv) {
	WNDCLASSEXA window_class = {};
	window_class.cbSize = sizeof(window_class);
	window_class.style = CS_CLASSDC;
	window_class.lpfnWndProc = WindowProcedure;
	window_class.hInstance = GetModuleHandleA(NULL);
	window_class.hCursor = LoadCursor(NULL, IDC_ARROW);
	window_class.lpszClassName = "TimelineFXExample";
	RegisterClassExA(&window_class);

	RECT window_rect = { 0, 0, EXAMPLE_WINDOW_WIDTH, EXAMPLE_WINDOW_HEIGHT };
	AdjustWindowRect(&window_rect, WS_OVERLAPPEDWINDOW, FALSE);
	example.window = CreateWindowExA(0, window_class.lpszClassName, "TimelineFX D3D11 Example", WS_OVERLAPPEDWINDOW, 100, 100,
		window_rect.right - window_rect.left, window_rect.bottom - window_rect.top, NULL, NULL, window_class.hInstance, NULL);
	if (!example.window || !CreateDevice(&example)) {
		ReleaseDevice(&example);
		if (example.window) {
			DestroyWindow(example.window);
		}
		UnregisterClassA(window_class.lpszClassName, window_class.hInstance);
		return 1;
	}
	ShowWindow(example.window, SW_SHOWDEFAULT);
	UpdateWindow(example.window);

	//Starts the worker threads and the memory pool that everything in the library allocates from
	tfx_BeginTimelineFX(tfx_GetDefaultThreadCount(), EXAMPLE_MEMORY_POOL_SIZE);
	tfx_impl_atlas_t atlas = {};
	bool ready = tfx_impl_dx11_Init(&example.renderer, example.device, example.device_context, EXAMPLE_SHADER_PATH, EXAMPLE_KERNEL_PATH)
		&& example_scene_Load(&example.scene, EXAMPLE_LIBRARY_PATH, EXAMPLE_SPRITE_DATA_PATH, EXAMPLE_FEATURES, &atlas)
		&& tfx_impl_dx11_UploadLibrary(&example.renderer, example.scene.library, &atlas)
		&& (!example.scene.animation_manager || tfx_impl_dx11_UploadSpriteData(&example.renderer, example.scene.animation_manager, example.scene.baked_shapes));
	tfx_impl_FreeAtlasPixels(&atlas);
	if (ready) {
		example_scene_Start(&example.scene);
		if (argc > 1 && atof(argv[1]) > 0.0) {
			example.scene.timer = tfx_impl_CreateTimer(atof(argv[1]));
		}
		printf("Space: burst  O: pause the camera  Mouse wheel: zoom  Escape: quit\n");
	}

	LARGE_INTEGER counter_frequency;
	LARGE_INTEGER last_counter;
	QueryPerformanceFrequency(&counter_frequency);
	QueryPerformanceCounter(&last_counter);
	double title_timer = 0.0;
	UINT frame_count = 0;

	bool running = ready;
	while (running) {
		MSG message;
		while (PeekMessageA(&message, NULL, 0, 0, PM_REMOVE)) {
			TranslateMessage(&message);
			DispatchMessageA(&message);
			if (message.message == WM_QUIT) {
				running = false;
			}
		}
		if (!running) {
			break;
		}

		if (example.swap_chain_occluded && example.swap_chain->Present(0, DXGI_PRESENT_TEST) == DXGI_STATUS_OCCLUDED) {
			Sleep(10);
			continue;
		}
		example.swap_chain_occluded = false;

		if (example.resize_width != 0 && example.resize_height != 0) {
			ReleaseRenderTargets(&example);
			example.swap_chain->ResizeBuffers(0, example.resize_width, example.resize_height, DXGI_FORMAT_UNKNOWN, 0);
			example.resize_width = example.resize_height = 0;
			if (!CreateRenderTargets(&example)) {
				break;
			}
		}

		LARGE_INTEGER counter;
		QueryPerformanceCounter(&counter);
		double elapsed_milliseconds = (double)(counter.QuadPart - last_counter.QuadPart) * 1000.0 / (double)counter_frequency.QuadPart;
		last_counter = counter;

		tfx_impl_frame_t frame;
		tfx_impl_uniforms_t uniforms;
		example_scene_Update(&example.scene, elapsed_milliseconds, (float)example.width, (float)example.height, tfx_impl_clip_space_depth_zero_to_one, &frame, &uniforms);
		tfx_impl_dx11_UploadFrame(&example.renderer, &frame);
		tfx_impl_dx11_Compute(&example.renderer, &uniforms);

		D3D11_VIEWPORT viewport = { 0.f, 0.f, (float)example.width, (float)example.height, 0.f, 1.f };
		const float clear_color[4] = { 0.f, 0.1f, 0.2f, 1.f };
		example.device_context->OMSetRenderTargets(1, &example.render_target_view, example.depth_stencil_view);
		example.device_context->RSSetViewports(1, &viewport);
		example.device_context->ClearRenderTargetView(example.render_target_view, clear_color);
		example.device_context->ClearDepthStencilView(example.depth_stencil_view, D3D11_CLEAR_DEPTH, 1.f, 0);
		tfx_impl_dx11_Render(&example.renderer, &uniforms);

		HRESULT present_result = example.swap_chain->Present(1, 0);
		example.swap_chain_occluded = present_result == DXGI_STATUS_OCCLUDED;
		PrintDebugLayerMessages(&example);

		frame_count++;
		title_timer += elapsed_milliseconds;
		if (title_timer >= 1000.0) {
			UpdateWindowTitle(&example, frame_count);
			title_timer -= 1000.0;
			frame_count = 0;
		}
	}

	example_scene_Shutdown(&example.scene);
	tfx_EndTimelineFX();
	tfx_impl_dx11_Shutdown(&example.renderer);

	ReleaseDevice(&example);
	if (IsWindow(example.window)) {
		DestroyWindow(example.window);
	}
	UnregisterClassA(window_class.lpszClassName, window_class.hInstance);
	return 0;
}
