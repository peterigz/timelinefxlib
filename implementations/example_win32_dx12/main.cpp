//TimelineFX example: Win32 + Direct3D 12 with the bindless shaders. The window, device, swap chain and main loop, see ../README.md

#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <d3d12.h>
#include <dxgi1_4.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "timelinefx.h"
#include "tfx_impl_common.h"
#include "tfx_impl_dx12.h"
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
#define TFX_EXAMPLE_SHADER_FOLDER "dx12_shaders"
#endif
#define EXAMPLE_FEATURES (example_scene_feature_ribbons | example_scene_feature_gpu_simulation | example_scene_feature_baked)

#define EXAMPLE_WINDOW_WIDTH 1280
#define EXAMPLE_WINDOW_HEIGHT 800
#define EXAMPLE_MEMORY_POOL_SIZE (64ull * 1024ull * 1024ull)
//One frame slot per back buffer, so the slot is the back buffer index
#define EXAMPLE_FRAME_COUNT 2
#define EXAMPLE_COLOR_FORMAT DXGI_FORMAT_R8G8B8A8_UNORM
#define EXAMPLE_DEPTH_FORMAT DXGI_FORMAT_D32_FLOAT
#define EXAMPLE_RELEASE(object) if (object) { (object)->Release(); (object) = NULL; }

typedef struct example_s {
	HWND window;
	ID3D12Device *device;
	ID3D12CommandQueue *queue;
	IDXGISwapChain3 *swap_chain;
	ID3D12DescriptorHeap *render_target_heap;
	ID3D12DescriptorHeap *depth_heap;
	ID3D12Resource *back_buffers[EXAMPLE_FRAME_COUNT];
	ID3D12Resource *depth_buffer;
	ID3D12CommandAllocator *command_allocators[EXAMPLE_FRAME_COUNT];
	ID3D12GraphicsCommandList *command_list;
	ID3D12Fence *fence;
	HANDLE fence_event;
	UINT64 fence_value;
	UINT64 frame_fence_values[EXAMPLE_FRAME_COUNT];
	UINT render_target_descriptor_size;
	UINT width;
	UINT height;
	UINT resize_width;				//Set by WM_SIZE, applied at the start of the next frame
	UINT resize_height;
	ID3D12InfoQueue *info_queue;	//Debug layer messages, null without the debug layer

	tfx_impl_dx12_t renderer;
	example_scene_t scene;
} example_t;

static example_t example;

//The debug layer writes to the debugger output, this repeats it on the console
static void PrintDebugLayerMessages(example_t *state) {
	if (!state->info_queue) {
		return;
	}
	UINT64 message_count = state->info_queue->GetNumStoredMessages();
	for (UINT64 message_index = 0; message_index != message_count; ++message_index) {
		SIZE_T message_size = 0;
		state->info_queue->GetMessage(message_index, NULL, &message_size);
		D3D12_MESSAGE *message = (D3D12_MESSAGE *)malloc(message_size);
		if (message && SUCCEEDED(state->info_queue->GetMessage(message_index, message, &message_size))) {
			printf("D3D12: %.*s\n", (int)message->DescriptionByteLength, message->pDescription);
		}
		free(message);
	}
	state->info_queue->ClearStoredMessages();
}

static void WaitForGpu(example_t *state) {
	state->fence_value++;
	state->queue->Signal(state->fence, state->fence_value);
	if (state->fence->GetCompletedValue() < state->fence_value) {
		state->fence->SetEventOnCompletion(state->fence_value, state->fence_event);
		WaitForSingleObject(state->fence_event, INFINITE);
	}
}

static void ReleaseRenderTargets(example_t *state) {
	for (UINT index = 0; index != EXAMPLE_FRAME_COUNT; ++index) {
		EXAMPLE_RELEASE(state->back_buffers[index]);
	}
	EXAMPLE_RELEASE(state->depth_buffer);
}

static D3D12_CPU_DESCRIPTOR_HANDLE RenderTargetView(example_t *state, UINT index) {
	D3D12_CPU_DESCRIPTOR_HANDLE handle = state->render_target_heap->GetCPUDescriptorHandleForHeapStart();
	handle.ptr += (SIZE_T)index * state->render_target_descriptor_size;
	return handle;
}

static bool CreateRenderTargets(example_t *state) {
	DXGI_SWAP_CHAIN_DESC1 swap_chain_description;
	state->swap_chain->GetDesc1(&swap_chain_description);
	state->width = swap_chain_description.Width;
	state->height = swap_chain_description.Height;
	for (UINT index = 0; index != EXAMPLE_FRAME_COUNT; ++index) {
		if (FAILED(state->swap_chain->GetBuffer(index, IID_PPV_ARGS(&state->back_buffers[index])))) {
			printf("Could not get swap chain buffer %u\n", index);
			return false;
		}
		state->device->CreateRenderTargetView(state->back_buffers[index], NULL, RenderTargetView(state, index));
	}
	D3D12_HEAP_PROPERTIES heap = {};
	heap.Type = D3D12_HEAP_TYPE_DEFAULT;
	D3D12_RESOURCE_DESC depth_description = {};
	depth_description.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;
	depth_description.Width = state->width;
	depth_description.Height = state->height;
	depth_description.DepthOrArraySize = 1;
	depth_description.MipLevels = 1;
	depth_description.Format = EXAMPLE_DEPTH_FORMAT;
	depth_description.SampleDesc.Count = 1;
	depth_description.Flags = D3D12_RESOURCE_FLAG_ALLOW_DEPTH_STENCIL;
	D3D12_CLEAR_VALUE clear_value = {};
	clear_value.Format = EXAMPLE_DEPTH_FORMAT;
	clear_value.DepthStencil.Depth = 1.f;
	HRESULT result = state->device->CreateCommittedResource(&heap, D3D12_HEAP_FLAG_NONE, &depth_description, D3D12_RESOURCE_STATE_DEPTH_WRITE, &clear_value, IID_PPV_ARGS(&state->depth_buffer));
	if (FAILED(result)) {
		printf("Could not create the depth buffer (0x%08lx)\n", result);
		return false;
	}
	state->device->CreateDepthStencilView(state->depth_buffer, NULL, state->depth_heap->GetCPUDescriptorHandleForHeapStart());
	return true;
}

//GPU based validation checks what the shaders do with each descriptor and resource state, which the debug layer can't see. It's slow
static bool CreateDevice(example_t *state, bool gpu_validation) {
	UINT factory_flags = 0;
#ifdef _DEBUG
	ID3D12Debug1 *debug = NULL;
	if (SUCCEEDED(D3D12GetDebugInterface(IID_PPV_ARGS(&debug)))) {
		debug->EnableDebugLayer();
		debug->SetEnableGPUBasedValidation(gpu_validation ? TRUE : FALSE);
		debug->Release();
		factory_flags |= DXGI_CREATE_FACTORY_DEBUG;
	} else {
		printf("The D3D12 debug layer isn't installed (Graphics Tools optional feature), continuing without it\n");
	}
#endif
	IDXGIFactory4 *factory = NULL;
	HRESULT result = CreateDXGIFactory2(factory_flags, IID_PPV_ARGS(&factory));
	if (FAILED(result)) {
		printf("Could not create a DXGI factory (0x%08lx)\n", result);
		return false;
	}
	result = D3D12CreateDevice(NULL, D3D_FEATURE_LEVEL_11_0, IID_PPV_ARGS(&state->device));
	if (FAILED(result)) {
		printf("Could not create a D3D12 device (0x%08lx)\n", result);
		factory->Release();
		return false;
	}
	state->device->QueryInterface(IID_PPV_ARGS(&state->info_queue));

	D3D12_COMMAND_QUEUE_DESC queue_description = {};
	queue_description.Type = D3D12_COMMAND_LIST_TYPE_DIRECT;
	result = state->device->CreateCommandQueue(&queue_description, IID_PPV_ARGS(&state->queue));
	DXGI_SWAP_CHAIN_DESC1 swap_chain_description = {};
	swap_chain_description.Format = EXAMPLE_COLOR_FORMAT;
	swap_chain_description.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
	swap_chain_description.BufferCount = EXAMPLE_FRAME_COUNT;
	swap_chain_description.SampleDesc.Count = 1;
	swap_chain_description.SwapEffect = DXGI_SWAP_EFFECT_FLIP_DISCARD;
	IDXGISwapChain1 *swap_chain = NULL;
	if (SUCCEEDED(result)) {
		result = factory->CreateSwapChainForHwnd(state->queue, state->window, &swap_chain_description, NULL, NULL, &swap_chain);
	}
	if (SUCCEEDED(result)) {
		result = swap_chain->QueryInterface(IID_PPV_ARGS(&state->swap_chain));
		swap_chain->Release();
	}
	factory->Release();
	if (FAILED(result)) {
		printf("Could not create the command queue and swap chain (0x%08lx)\n", result);
		return false;
	}

	D3D12_DESCRIPTOR_HEAP_DESC render_target_heap = {};
	render_target_heap.Type = D3D12_DESCRIPTOR_HEAP_TYPE_RTV;
	render_target_heap.NumDescriptors = EXAMPLE_FRAME_COUNT;
	D3D12_DESCRIPTOR_HEAP_DESC depth_heap = {};
	depth_heap.Type = D3D12_DESCRIPTOR_HEAP_TYPE_DSV;
	depth_heap.NumDescriptors = 1;
	if (FAILED(state->device->CreateDescriptorHeap(&render_target_heap, IID_PPV_ARGS(&state->render_target_heap)))
		|| FAILED(state->device->CreateDescriptorHeap(&depth_heap, IID_PPV_ARGS(&state->depth_heap)))) {
		printf("Could not create the render target heaps\n");
		return false;
	}
	state->render_target_descriptor_size = state->device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_RTV);
	for (UINT index = 0; index != EXAMPLE_FRAME_COUNT; ++index) {
		if (FAILED(state->device->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_DIRECT, IID_PPV_ARGS(&state->command_allocators[index])))) {
			return false;
		}
	}
	if (FAILED(state->device->CreateCommandList(0, D3D12_COMMAND_LIST_TYPE_DIRECT, state->command_allocators[0], NULL, IID_PPV_ARGS(&state->command_list)))
		|| FAILED(state->device->CreateFence(0, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&state->fence)))) {
		printf("Could not create the command list and fence\n");
		return false;
	}
	state->command_list->Close();
	state->fence_event = CreateEventA(NULL, FALSE, FALSE, NULL);
	return state->fence_event && CreateRenderTargets(state);
}

static void ReleaseDevice(example_t *state) {
	ReleaseRenderTargets(state);
	EXAMPLE_RELEASE(state->command_list);
	for (UINT index = 0; index != EXAMPLE_FRAME_COUNT; ++index) {
		EXAMPLE_RELEASE(state->command_allocators[index]);
	}
	EXAMPLE_RELEASE(state->fence);
	EXAMPLE_RELEASE(state->render_target_heap);
	EXAMPLE_RELEASE(state->depth_heap);
	EXAMPLE_RELEASE(state->swap_chain);
	EXAMPLE_RELEASE(state->queue);
	if (state->fence_event) {
		CloseHandle(state->fence_event);
		state->fence_event = NULL;
	}
	//Anything still alive apart from the device itself is a leak
	ID3D12DebugDevice *debug_device = NULL;
	if (state->device && SUCCEEDED(state->device->QueryInterface(IID_PPV_ARGS(&debug_device)))) {
		debug_device->ReportLiveDeviceObjects(D3D12_RLDO_DETAIL | D3D12_RLDO_IGNORE_INTERNAL);
		debug_device->Release();
		PrintDebugLayerMessages(state);
	}
	EXAMPLE_RELEASE(state->info_queue);
	EXAMPLE_RELEASE(state->device);
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
		} else if (wparam == 'L') {
			example_scene_ReloadSecondLibrary(&example.scene);
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

//Every atlas after the load, and the second library's after a reload
static bool UploadPendingAtlases(example_t *state) {
	bool uploaded = true;
	tfx_impl_atlas_t *atlas = example_scene_GetPendingAtlas(&state->scene);
	while (atlas) {
		uploaded = tfx_impl_dx12_UploadAtlas(&state->renderer, state->command_list, atlas) && uploaded;
		tfx_impl_FreeAtlasPixels(atlas);
		atlas = example_scene_GetPendingAtlas(&state->scene);
	}
	return uploaded;
}

//The atlas and sprite data go up in one command list, waited on before the first frame
static bool LoadScene(example_t *state) {
	if (!example_scene_Load(&state->scene, EXAMPLE_LIBRARY_PATH, EXAMPLE_SECOND_LIBRARY_PATH, EXAMPLE_SPRITE_DATA_PATH, EXAMPLE_FEATURES)
		|| !tfx_impl_dx12_CreateAtlas(&state->renderer, state->scene.atlas_layer_size, state->scene.atlas_layer_count)) {
		return false;
	}
	state->command_allocators[0]->Reset();
	state->command_list->Reset(state->command_allocators[0], NULL);
	tfx_impl_dx12_BeginFrame(&state->renderer, 0);
	bool uploaded = UploadPendingAtlases(state)
		&& (!state->scene.animation_manager || tfx_impl_dx12_UploadSpriteData(&state->renderer, state->command_list, state->scene.animation_manager, state->scene.baked_shapes));
	state->command_list->Close();
	ID3D12CommandList *lists[] = { state->command_list };
	state->queue->ExecuteCommandLists(1, lists);
	WaitForGpu(state);
	return uploaded;
}

static void UpdateWindowTitle(example_t *state, UINT frames_per_second) {
	tfx_version_t version = tfx_GetVersion();
	char title[256];
	tfx_gpu_particle_stats_t gpu_stats = tfx_GetGPUParticleStats(state->scene.stage);
	snprintf(title, sizeof(title), "TimelineFX D3D12 Example - %s %i.%i - %u fps - %u CPU particles - %u GPU particles - %llu dropped GPU spawns",
		version.name, version.major, version.minor, frames_per_second, tfx_GetParticleCount(state->scene.stage), gpu_stats.estimated_live_particles,
		(unsigned long long)gpu_stats.dropped_spawns);
	SetWindowTextA(state->window, title);
}

static void Barrier(ID3D12GraphicsCommandList *command_list, ID3D12Resource *resource, D3D12_RESOURCE_STATES before, D3D12_RESOURCE_STATES after) {
	D3D12_RESOURCE_BARRIER barrier = {};
	barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
	barrier.Transition.pResource = resource;
	barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
	barrier.Transition.StateBefore = before;
	barrier.Transition.StateAfter = after;
	command_list->ResourceBarrier(1, &barrier);
}

//Usage: example_win32_dx12 [updates per second] [--frames N] [--test-debug-output] [--gpu-validation]
int main(int argc, char **argv) {
	double updates_per_second = 0.0;
	long frame_limit = 0;
	bool test_debug_output = false;
	bool gpu_validation = false;
	for (int index = 1; index < argc; ++index) {
		if (strcmp(argv[index], "--gpu-validation") == 0) {
			gpu_validation = true;
		} else if (strcmp(argv[index], "--frames") == 0 && index + 1 < argc) {
			frame_limit = atol(argv[++index]);
		} else if (strcmp(argv[index], "--test-debug-output") == 0) {
			test_debug_output = true;
		} else if (atof(argv[index]) > 0.0) {
			updates_per_second = atof(argv[index]);
		}
	}

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
	example.window = CreateWindowExA(0, window_class.lpszClassName, "TimelineFX D3D12 Example", WS_OVERLAPPEDWINDOW, 100, 100,
		window_rect.right - window_rect.left, window_rect.bottom - window_rect.top, NULL, NULL, window_class.hInstance, NULL);
	if (!example.window || !CreateDevice(&example, gpu_validation)) {
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
	bool ready = tfx_impl_dx12_Init(&example.renderer, example.device, EXAMPLE_COLOR_FORMAT, EXAMPLE_DEPTH_FORMAT, TFX_EXAMPLE_SHADER_FOLDER, EXAMPLE_FRAME_COUNT)
		&& LoadScene(&example);
	PrintDebugLayerMessages(&example);
	if (ready) {
		//A zero size buffer, to prove the debug layer's messages reach the console before trusting their silence
		if (test_debug_output) {
			ID3D12Resource *invalid = NULL;
			D3D12_HEAP_PROPERTIES heap = {};
			heap.Type = D3D12_HEAP_TYPE_DEFAULT;
			D3D12_RESOURCE_DESC description = {};
			description.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
			description.Height = 1;
			description.DepthOrArraySize = 1;
			description.MipLevels = 1;
			description.SampleDesc.Count = 1;
			description.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
			example.device->CreateCommittedResource(&heap, D3D12_HEAP_FLAG_NONE, &description, D3D12_RESOURCE_STATE_COMMON, NULL, IID_PPV_ARGS(&invalid));
			EXAMPLE_RELEASE(invalid);
			PrintDebugLayerMessages(&example);
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

	LARGE_INTEGER counter_frequency;
	LARGE_INTEGER last_counter;
	QueryPerformanceFrequency(&counter_frequency);
	QueryPerformanceCounter(&last_counter);
	double title_timer = 0.0;
	UINT frame_count = 0;
	long frames_run = 0;

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

		if (example.resize_width != 0 && example.resize_height != 0) {
			WaitForGpu(&example);
			ReleaseRenderTargets(&example);
			example.swap_chain->ResizeBuffers(0, example.resize_width, example.resize_height, DXGI_FORMAT_UNKNOWN, 0);
			example.resize_width = example.resize_height = 0;
			if (!CreateRenderTargets(&example)) {
				break;
			}
		}

		//Wait until the GPU has finished the last frame that used this slot's allocator and upload memory
		UINT frame_index = example.swap_chain->GetCurrentBackBufferIndex();
		if (example.fence->GetCompletedValue() < example.frame_fence_values[frame_index]) {
			example.fence->SetEventOnCompletion(example.frame_fence_values[frame_index], example.fence_event);
			WaitForSingleObject(example.fence_event, INFINITE);
		}
		example.command_allocators[frame_index]->Reset();
		example.command_list->Reset(example.command_allocators[frame_index], NULL);
		tfx_impl_dx12_BeginFrame(&example.renderer, frame_index);

		LARGE_INTEGER counter;
		QueryPerformanceCounter(&counter);
		double elapsed_milliseconds = (double)(counter.QuadPart - last_counter.QuadPart) * 1000.0 / (double)counter_frequency.QuadPart;
		last_counter = counter;

		tfx_impl_frame_t frame;
		tfx_impl_uniforms_t uniforms;
		example_scene_Update(&example.scene, elapsed_milliseconds, (float)example.width, (float)example.height, tfx_impl_clip_space_depth_zero_to_one, &frame, &uniforms);
		UploadPendingAtlases(&example);
		tfx_impl_dx12_UploadFrame(&example.renderer, example.command_list, &frame);
		tfx_impl_dx12_Compute(&example.renderer, example.command_list, &uniforms);

		Barrier(example.command_list, example.back_buffers[frame_index], D3D12_RESOURCE_STATE_PRESENT, D3D12_RESOURCE_STATE_RENDER_TARGET);
		D3D12_CPU_DESCRIPTOR_HANDLE render_target = RenderTargetView(&example, frame_index);
		D3D12_CPU_DESCRIPTOR_HANDLE depth_target = example.depth_heap->GetCPUDescriptorHandleForHeapStart();
		const float clear_color[4] = { 0.f, 0.1f, 0.2f, 1.f };
		D3D12_VIEWPORT viewport = { 0.f, 0.f, (float)example.width, (float)example.height, 0.f, 1.f };
		D3D12_RECT scissor = { 0, 0, (LONG)example.width, (LONG)example.height };
		example.command_list->OMSetRenderTargets(1, &render_target, FALSE, &depth_target);
		example.command_list->ClearRenderTargetView(render_target, clear_color, 0, NULL);
		example.command_list->ClearDepthStencilView(depth_target, D3D12_CLEAR_FLAG_DEPTH, 1.f, 0, 0, NULL);
		example.command_list->RSSetViewports(1, &viewport);
		example.command_list->RSSetScissorRects(1, &scissor);
		tfx_impl_dx12_Render(&example.renderer, example.command_list, &uniforms);
		Barrier(example.command_list, example.back_buffers[frame_index], D3D12_RESOURCE_STATE_RENDER_TARGET, D3D12_RESOURCE_STATE_PRESENT);
		example.command_list->Close();
		ID3D12CommandList *lists[] = { example.command_list };
		example.queue->ExecuteCommandLists(1, lists);
		example.swap_chain->Present(1, 0);
		example.fence_value++;
		example.queue->Signal(example.fence, example.fence_value);
		example.frame_fence_values[frame_index] = example.fence_value;
		PrintDebugLayerMessages(&example);

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

	if (example.queue) {
		WaitForGpu(&example);
	}
	example_scene_Shutdown(&example.scene);
	tfx_EndTimelineFX();
	tfx_impl_dx12_Shutdown(&example.renderer);

	ReleaseDevice(&example);
	if (IsWindow(example.window)) {
		DestroyWindow(example.window);
	}
	UnregisterClassA(window_class.lpszClassName, window_class.hInstance);
	return 0;
}
