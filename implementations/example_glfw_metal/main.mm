//TimelineFX example: GLFW + Metal. The window, device, CAMetalLayer, frames in flight and main loop, see ../README.md

#ifndef GLFW_INCLUDE_NONE
#define GLFW_INCLUDE_NONE
#endif
#include <GLFW/glfw3.h>
#define GLFW_EXPOSE_NATIVE_COCOA
#include <GLFW/glfw3native.h>
#import <Cocoa/Cocoa.h>
#import <Metal/Metal.h>
#import <QuartzCore/CAMetalLayer.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "timelinefx.h"
#include "tfx_impl_common.h"
#include "tfx_impl_metal.h"
#include "example_scene.h"

#if __has_feature(objc_arc)
#error "Build main.mm without ARC (-fno-objc-arc), like tfx_impl_metal.mm"
#endif

//TFX_EXAMPLE_ROOT is set by the CMake project to the implementations folder
#ifndef TFX_EXAMPLE_ROOT
#define TFX_EXAMPLE_ROOT "."
#endif
#define EXAMPLE_LIBRARY_PATH TFX_EXAMPLE_ROOT "/example.tfx"
#define EXAMPLE_SECOND_LIBRARY_PATH TFX_EXAMPLE_ROOT "/example_second.tfx"
#define EXAMPLE_SPRITE_DATA_PATH TFX_EXAMPLE_ROOT "/example.tfxsd"
//TFX_EXAMPLE_SHADER_FOLDER is where the CMake project has slangc build the shaders
#ifndef TFX_EXAMPLE_SHADER_FOLDER
#define TFX_EXAMPLE_SHADER_FOLDER "metal_shaders"
#endif
#define EXAMPLE_FEATURES (example_scene_feature_ribbons | example_scene_feature_gpu_simulation | example_scene_feature_baked)

#define EXAMPLE_WINDOW_WIDTH 1280
#define EXAMPLE_WINDOW_HEIGHT 800
#define EXAMPLE_MEMORY_POOL_SIZE (64ull * 1024ull * 1024ull)
#define EXAMPLE_FRAMES_IN_FLIGHT 3
//Not sRGB, to match the other examples, whose shaders write the colors as they are
#define EXAMPLE_COLOR_FORMAT MTLPixelFormatBGRA8Unorm
#define EXAMPLE_DEPTH_FORMAT MTLPixelFormatDepth32Float

typedef struct example_s {
	GLFWwindow *window;
	id<MTLDevice> device;
	id<MTLCommandQueue> queue;
	CAMetalLayer *layer;
	id<MTLTexture> depth_texture;
	uint32_t width;					//The drawable size the depth texture was made for
	uint32_t height;
	//Counts the frame slots the GPU has finished with, signalled from each command buffer's completion handler
	dispatch_semaphore_t frame_semaphore;
	tfxU32 frame_index;

	tfx_impl_metal_t renderer;
	example_scene_t scene;
} example_t;

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

//Metal reads these when the device is created, so setting them here works as well as setting them in the shell. Its messages go to NSLog,
//which prints on stderr. In nslog mode an error is printed rather than stopping the program, like the other examples
static void SetValidationEnvironment(bool gpu_validation) {
#ifndef NDEBUG
	setenv("MTL_DEBUG_LAYER", "1", 0);
	setenv("MTL_DEBUG_LAYER_ERROR_MODE", "nslog", 0);
#endif
	if (gpu_validation) {
		setenv("MTL_SHADER_VALIDATION", "1", 0);
	}
}

static bool CreateDevice(example_t *state) {
	state->device = MTLCreateSystemDefaultDevice();
	if (!state->device) {
		printf("No Metal device\n");
		return false;
	}
	printf("%s, Metal API validation %s, shader validation %s\n", [[state->device name] UTF8String],
		getenv("MTL_DEBUG_LAYER") ? "on" : "off", getenv("MTL_SHADER_VALIDATION") ? "on" : "off");
	state->queue = [state->device newCommandQueue];
	state->frame_semaphore = dispatch_semaphore_create(EXAMPLE_FRAMES_IN_FLIGHT);
	return state->queue != nil;
}

static bool CreateLayer(example_t *state) {
	NSWindow *window = glfwGetCocoaWindow(state->window);
	state->layer = [[CAMetalLayer layer] retain];
	state->layer.device = state->device;
	state->layer.pixelFormat = EXAMPLE_COLOR_FORMAT;
	state->layer.contentsScale = [window backingScaleFactor];
	window.contentView.layer = state->layer;
	window.contentView.wantsLayer = YES;
	return true;
}

//The drawable and the depth texture follow the framebuffer, which is in pixels. The depth is memoryless: it's cleared and never stored
static bool ResizeDrawable(example_t *state, uint32_t width, uint32_t height) {
	state->layer.drawableSize = CGSizeMake(width, height);
	[state->depth_texture release];
	MTLTextureDescriptor *descriptor = [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:EXAMPLE_DEPTH_FORMAT width:width height:height mipmapped:NO];
	descriptor.usage = MTLTextureUsageRenderTarget;
	descriptor.storageMode = MTLStorageModeMemoryless;
	state->depth_texture = [state->device newTextureWithDescriptor:descriptor];
	state->width = width;
	state->height = height;
	return state->depth_texture != nil;
}

static id<MTLCommandBuffer> CreateCommandBuffer(example_t *state) {
	MTLCommandBufferDescriptor *descriptor = [[MTLCommandBufferDescriptor new] autorelease];
#ifndef NDEBUG
	//Names the encoder that failed in the error the completion handler prints
	descriptor.errorOptions = MTLCommandBufferErrorOptionEncoderExecutionStatus;
#endif
	id<MTLCommandBuffer> command_buffer = [state->queue commandBufferWithDescriptor:descriptor];
	//Shader validation leaves the command buffer completed and reports in its logs
	[command_buffer addCompletedHandler:^(id<MTLCommandBuffer> completed) {
		if ([completed status] == MTLCommandBufferStatusError) {
			printf("Metal command buffer error: %s\n", [[[completed error] description] UTF8String]);
		}
		for (id<MTLFunctionLog> log in [completed logs]) {
			printf("Metal shader validation: %s\n", [[log description] UTF8String]);
		}
	}];
	return command_buffer;
}

//Every atlas after the load, and the second library's after a reload. Its pixels are copied into the frame's upload memory, so they can go at once
static bool UploadPendingAtlases(example_t *state, id<MTLCommandBuffer> command_buffer) {
	bool uploaded = true;
	tfx_impl_atlas_t *atlas = example_scene_GetPendingAtlas(&state->scene);
	while (atlas) {
		uploaded = tfx_impl_metal_UploadAtlas(&state->renderer, command_buffer, atlas) && uploaded;
		tfx_impl_FreeAtlasPixels(atlas);
		atlas = example_scene_GetPendingAtlas(&state->scene);
	}
	return uploaded;
}

//The atlases and the sprite data, before the first frame, waiting for the GPU so the first frame slot's upload memory is free again
static bool UploadScene(example_t *state) {
	id<MTLCommandBuffer> command_buffer = CreateCommandBuffer(state);
	bool uploaded = UploadPendingAtlases(state, command_buffer)
		&& (!state->scene.animation_manager || tfx_impl_metal_UploadSpriteData(&state->renderer, command_buffer, state->scene.animation_manager, state->scene.baked_shapes));
	[command_buffer commit];
	[command_buffer waitUntilCompleted];
	return uploaded && [command_buffer status] == MTLCommandBufferStatusCompleted;
}

static void UpdateWindowTitle(example_t *state, unsigned int frames_per_second) {
	tfx_version_t version = tfx_GetVersion();
	char title[256];
	tfx_gpu_particle_stats_t gpu_stats = tfx_GetGPUParticleStats(state->scene.stage);
	snprintf(title, sizeof(title), "TimelineFX Metal Example - %s %i.%i - %u fps - %u CPU particles - %u GPU particles - %llu dropped GPU spawns",
		version.name, version.major, version.minor, frames_per_second, tfx_GetParticleCount(state->scene.stage), gpu_stats.estimated_live_particles,
		(unsigned long long)gpu_stats.dropped_spawns);
	glfwSetWindowTitle(state->window, title);
}

static void RenderFrame(example_t *state, id<MTLCommandBuffer> command_buffer, id<CAMetalDrawable> drawable, const tfx_impl_uniforms_t *uniforms) {
	MTLRenderPassDescriptor *pass = [MTLRenderPassDescriptor renderPassDescriptor];
	pass.colorAttachments[0].texture = [drawable texture];
	pass.colorAttachments[0].loadAction = MTLLoadActionClear;
	pass.colorAttachments[0].storeAction = MTLStoreActionStore;
	pass.colorAttachments[0].clearColor = MTLClearColorMake(0.0, 0.1, 0.2, 1.0);
	pass.depthAttachment.texture = state->depth_texture;
	pass.depthAttachment.loadAction = MTLLoadActionClear;
	pass.depthAttachment.storeAction = MTLStoreActionDontCare;
	pass.depthAttachment.clearDepth = 1.0;
	id<MTLRenderCommandEncoder> render_encoder = [command_buffer renderCommandEncoderWithDescriptor:pass];
	tfx_impl_metal_Render(&state->renderer, render_encoder, uniforms);
	[render_encoder endEncoding];
	[command_buffer presentDrawable:drawable];
}

//Usage: example_glfw_metal [updates per second] [--frames N] [--test-debug-output] [--gpu-validation]
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

	@autoreleasepool {
		SetValidationEnvironment(gpu_validation);
		glfwSetErrorCallback(ErrorCallback);
		if (!glfwInit()) {
			return 1;
		}
		glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
		example.window = glfwCreateWindow(EXAMPLE_WINDOW_WIDTH, EXAMPLE_WINDOW_HEIGHT, "TimelineFX Metal Example", NULL, NULL);
		int width = 0;
		int height = 0;
		if (example.window) {
			glfwSetKeyCallback(example.window, KeyCallback);
			glfwSetScrollCallback(example.window, ScrollCallback);
			glfwGetFramebufferSize(example.window, &width, &height);
		}
		bool ready = example.window && CreateDevice(&example) && CreateLayer(&example) && ResizeDrawable(&example, (uint32_t)width, (uint32_t)height);

		//Starts the worker threads and the memory pool that everything in the library allocates from
		tfx_BeginTimelineFX(tfx_GetDefaultThreadCount(), EXAMPLE_MEMORY_POOL_SIZE);
		ready = ready && tfx_impl_metal_Init(&example.renderer, example.device, EXAMPLE_COLOR_FORMAT, EXAMPLE_DEPTH_FORMAT, TFX_EXAMPLE_SHADER_FOLDER, EXAMPLE_FRAMES_IN_FLIGHT)
			&& example_scene_Load(&example.scene, EXAMPLE_LIBRARY_PATH, EXAMPLE_SECOND_LIBRARY_PATH, EXAMPLE_SPRITE_DATA_PATH, EXAMPLE_FEATURES)
			&& tfx_impl_metal_CreateAtlas(&example.renderer, example.scene.atlas_layer_size, example.scene.atlas_layer_count)
			&& UploadScene(&example);
		if (ready) {
			//A buffer of zero length, to prove the validation layer's messages reach the console before trusting their silence
			if (test_debug_output) {
				[[example.device newBufferWithLength:0 options:MTLResourceStorageModeShared] release];
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
			@autoreleasepool {
				glfwPollEvents();
				glfwGetFramebufferSize(example.window, &width, &height);
				if (width == 0 || height == 0) {
					glfwWaitEvents();
					continue;
				}
				if ((uint32_t)width != example.width || (uint32_t)height != example.height) {
					if (!ResizeDrawable(&example, (uint32_t)width, (uint32_t)height)) {
						break;
					}
				}
				double time = glfwGetTime();
				double elapsed_milliseconds = (time - last_time) * 1000.0;
				last_time = time;

				//The slot's upload memory is only written once the GPU has finished the frame that last used it
				dispatch_semaphore_wait(example.frame_semaphore, DISPATCH_TIME_FOREVER);
				tfx_impl_metal_BeginFrame(&example.renderer, example.frame_index);

				tfx_impl_frame_t frame;
				tfx_impl_uniforms_t uniforms;
				example_scene_Update(&example.scene, elapsed_milliseconds, (float)example.width, (float)example.height, tfx_impl_clip_space_depth_zero_to_one, &frame, &uniforms);
				id<MTLCommandBuffer> command_buffer = CreateCommandBuffer(&example);
				dispatch_semaphore_t frame_semaphore = example.frame_semaphore;
				[command_buffer addCompletedHandler:^(id<MTLCommandBuffer> completed) {
					(void)completed;
					dispatch_semaphore_signal(frame_semaphore);
				}];
				UploadPendingAtlases(&example, command_buffer);
				tfx_impl_metal_UploadFrame(&example.renderer, command_buffer, &frame);
				//The compute runs whether or not there's a drawable to draw into, because every GPU particle fetch has to be dispatched
				tfx_impl_metal_Compute(&example.renderer, command_buffer, &uniforms);
				id<CAMetalDrawable> drawable = [example.layer nextDrawable];
				if (drawable) {
					RenderFrame(&example, command_buffer, drawable, &uniforms);
				}
				[command_buffer commit];
				example.frame_index = (example.frame_index + 1) % EXAMPLE_FRAMES_IN_FLIGHT;

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
		}

		//Every frame slot back means the GPU is idle
		if (example.frame_semaphore) {
			for (int index = 0; index != EXAMPLE_FRAMES_IN_FLIGHT; ++index) {
				dispatch_semaphore_wait(example.frame_semaphore, DISPATCH_TIME_FOREVER);
			}
		}
		example_scene_Shutdown(&example.scene);
		tfx_EndTimelineFX();
		tfx_impl_metal_Shutdown(&example.renderer);

		[example.depth_texture release];
		[example.layer release];
		[example.queue release];
		[example.device release];
		if (example.frame_semaphore) {
			//libdispatch won't release a semaphore below the value it was created with
			for (int index = 0; index != EXAMPLE_FRAMES_IN_FLIGHT; ++index) {
				dispatch_semaphore_signal(example.frame_semaphore);
			}
			dispatch_release(example.frame_semaphore);
		}
		if (example.window) {
			glfwDestroyWindow(example.window);
		}
		glfwTerminate();
	}
	return 0;
}
