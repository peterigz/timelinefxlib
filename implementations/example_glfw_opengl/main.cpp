//TimelineFX example: GLFW + OpenGL 4.3 core. The window, context and main loop, see ../README.md

#ifndef GLFW_INCLUDE_NONE
#define GLFW_INCLUDE_NONE
#endif
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "timelinefx.h"
#include "tfx_impl_common.h"
#include "tfx_impl_gl.h"
#include "example_scene.h"

//TFX_EXAMPLE_ROOT is set by the CMake project to the implementations folder
#ifndef TFX_EXAMPLE_ROOT
#define TFX_EXAMPLE_ROOT "."
#endif
#define EXAMPLE_LIBRARY_PATH TFX_EXAMPLE_ROOT "/example.tfx"
#define EXAMPLE_SECOND_LIBRARY_PATH TFX_EXAMPLE_ROOT "/example_second.tfx"
#define EXAMPLE_SPRITE_DATA_PATH TFX_EXAMPLE_ROOT "/example.tfxsd"
//TFX_EXAMPLE_SHADER_FOLDER is where the CMake project has slangc and spirv-cross build the shaders
#ifndef TFX_EXAMPLE_SHADER_FOLDER
#define TFX_EXAMPLE_SHADER_FOLDER "gl_shaders"
#endif
#define EXAMPLE_FEATURES (example_scene_feature_ribbons | example_scene_feature_gpu_simulation | example_scene_feature_baked)

#define EXAMPLE_WINDOW_WIDTH 1280
#define EXAMPLE_WINDOW_HEIGHT 800
#define EXAMPLE_MEMORY_POOL_SIZE (64ull * 1024ull * 1024ull)

//The few GL calls the example makes itself, the backend loads its own
#define EXAMPLE_GL_COLOR_BUFFER_BIT 0x00004000
#define EXAMPLE_GL_DEPTH_BUFFER_BIT 0x00000100
#define EXAMPLE_GL_RENDERER 0x1F01
#define EXAMPLE_GL_VERSION 0x1F02
#if defined(_WIN32) && !defined(_WIN64)
#define EXAMPLE_GL_API __stdcall
#else
#define EXAMPLE_GL_API
#endif

typedef struct example_gl_s {
	void (EXAMPLE_GL_API *Viewport)(int x, int y, int width, int height);
	void (EXAMPLE_GL_API *ClearColor)(float red, float green, float blue, float alpha);
	void (EXAMPLE_GL_API *Clear)(unsigned int mask);
	void (EXAMPLE_GL_API *Enable)(unsigned int capability);
	const unsigned char *(EXAMPLE_GL_API *GetString)(unsigned int name);
} example_gl_t;

typedef struct example_s {
	GLFWwindow *window;
	example_gl_t gl;
	tfx_impl_gl_t renderer;
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

static bool LoadExampleFunctions(example_gl_t *gl) {
	gl->Viewport = (void (EXAMPLE_GL_API *)(int, int, int, int))glfwGetProcAddress("glViewport");
	gl->ClearColor = (void (EXAMPLE_GL_API *)(float, float, float, float))glfwGetProcAddress("glClearColor");
	gl->Clear = (void (EXAMPLE_GL_API *)(unsigned int))glfwGetProcAddress("glClear");
	gl->Enable = (void (EXAMPLE_GL_API *)(unsigned int))glfwGetProcAddress("glEnable");
	gl->GetString = (const unsigned char *(EXAMPLE_GL_API *)(unsigned int))glfwGetProcAddress("glGetString");
	return gl->Viewport && gl->ClearColor && gl->Clear && gl->Enable && gl->GetString;
}

//Every atlas after the load, and the second library's after a reload
static bool UploadPendingAtlases(example_t *state) {
	bool uploaded = true;
	tfx_impl_atlas_t *atlas = example_scene_GetPendingAtlas(&state->scene);
	while (atlas) {
		uploaded = tfx_impl_gl_UploadAtlas(&state->renderer, atlas) && uploaded;
		tfx_impl_FreeAtlasPixels(atlas);
		atlas = example_scene_GetPendingAtlas(&state->scene);
	}
	return uploaded;
}

static void UpdateWindowTitle(example_t *state, unsigned int frames_per_second) {
	tfx_version_t version = tfx_GetVersion();
	char title[256];
	tfx_gpu_particle_stats_t gpu_stats = tfx_GetGPUParticleStats(state->scene.stage);
	snprintf(title, sizeof(title), "TimelineFX OpenGL Example - %s %i.%i - %u fps - %u CPU particles - %u GPU particles - %llu dropped GPU spawns",
		version.name, version.major, version.minor, frames_per_second, tfx_GetParticleCount(state->scene.stage), gpu_stats.estimated_live_particles,
		(unsigned long long)gpu_stats.dropped_spawns);
	glfwSetWindowTitle(state->window, title);
}

//Usage: example_glfw_opengl [updates per second] [--frames N] [--test-debug-output]
int main(int argc, char **argv) {
	double updates_per_second = 0.0;
	long frame_limit = 0;
	bool test_debug_output = false;
	for (int index = 1; index < argc; ++index) {
		if (strcmp(argv[index], "--frames") == 0 && index + 1 < argc) {
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
	//4.3 core for storage buffers in the vertex stage and compute shaders. A debug context so the driver's messages reach the log
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GLFW_TRUE);
	glfwWindowHint(GLFW_DEPTH_BITS, 24);
	example.window = glfwCreateWindow(EXAMPLE_WINDOW_WIDTH, EXAMPLE_WINDOW_HEIGHT, "TimelineFX OpenGL Example", NULL, NULL);
	if (!example.window) {
		printf("Could not create a window with an OpenGL 4.3 core context\n");
		glfwTerminate();
		return 1;
	}
	glfwMakeContextCurrent(example.window);
	glfwSwapInterval(1);
	glfwSetKeyCallback(example.window, KeyCallback);
	glfwSetScrollCallback(example.window, ScrollCallback);

	//Starts the worker threads and the memory pool that everything in the library allocates from
	tfx_BeginTimelineFX(tfx_GetDefaultThreadCount(), EXAMPLE_MEMORY_POOL_SIZE);
	bool ready = LoadExampleFunctions(&example.gl)
		&& tfx_impl_gl_Init(&example.renderer, glfwGetProcAddress, TFX_EXAMPLE_SHADER_FOLDER)
		&& example_scene_Load(&example.scene, EXAMPLE_LIBRARY_PATH, EXAMPLE_SECOND_LIBRARY_PATH, EXAMPLE_SPRITE_DATA_PATH, EXAMPLE_FEATURES)
		&& tfx_impl_gl_CreateAtlas(&example.renderer, example.scene.atlas_layer_size, example.scene.atlas_layer_count)
		&& UploadPendingAtlases(&example)
		&& (!example.scene.animation_manager || tfx_impl_gl_UploadSpriteData(&example.renderer, example.scene.animation_manager, example.scene.baked_shapes));
	if (ready) {
		printf("%s, %s\n", example.gl.GetString(EXAMPLE_GL_RENDERER), example.gl.GetString(EXAMPLE_GL_VERSION));
		//An invalid enum, to prove the debug messages reach the log before trusting their silence
		if (test_debug_output) {
			example.gl.Enable(0);
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
		int width = 0;
		int height = 0;
		glfwGetFramebufferSize(example.window, &width, &height);
		if (width == 0 || height == 0) {
			glfwWaitEvents();
			continue;
		}
		double time = glfwGetTime();
		double elapsed_milliseconds = (time - last_time) * 1000.0;
		last_time = time;

		tfx_impl_frame_t frame;
		tfx_impl_uniforms_t uniforms;
		example_scene_Update(&example.scene, elapsed_milliseconds, (float)width, (float)height, tfx_impl_clip_space_none, &frame, &uniforms);
		UploadPendingAtlases(&example);
		tfx_impl_gl_UploadFrame(&example.renderer, &frame);
		tfx_impl_gl_Compute(&example.renderer, &uniforms);

		example.gl.Viewport(0, 0, width, height);
		example.gl.ClearColor(0.f, 0.1f, 0.2f, 1.f);
		example.gl.Clear(EXAMPLE_GL_COLOR_BUFFER_BIT | EXAMPLE_GL_DEPTH_BUFFER_BIT);
		tfx_impl_gl_Render(&example.renderer, &uniforms);
		glfwSwapBuffers(example.window);

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

	example_scene_Shutdown(&example.scene);
	tfx_EndTimelineFX();
	tfx_impl_gl_Shutdown(&example.renderer);

	glfwDestroyWindow(example.window);
	glfwTerminate();
	return 0;
}
