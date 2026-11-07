#include "tfx_impl_gl.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TFX_IMPL_GL_INITIAL_INSTANCE_CAPACITY 10000
//tfx_GetGlobalColorRampLayerCount never goes past this
#define TFX_IMPL_GL_MAX_COLOR_RAMP_LAYERS 256
//Match the thread group sizes in timelinefx.slang
#define TFX_IMPL_GL_KERNEL_GROUP_SIZE 256
#define TFX_IMPL_GL_RIBBON_GROUP_SIZE 256
#define TFX_IMPL_GL_PLAYBACK_GROUP_SIZE 128

//The default fixed bindings in timelinefx.slang. SSBOs, UBOs and texture units are separate binding spaces in GL
#define TFX_IMPL_GL_VERTEX_READ_BINDING 0
#define TFX_IMPL_GL_VERTEX_UNIFORM_BINDING 5
#define TFX_IMPL_GL_FRAGMENT_TEXTURE_BINDING 7
#define TFX_IMPL_GL_FRAGMENT_UNIFORM_BINDING 10
#define TFX_IMPL_GL_COMPUTE_WRITE_BINDING 0
#define TFX_IMPL_GL_COMPUTE_READ_BINDING 2
#define TFX_IMPL_GL_COMPUTE_UNIFORM_BINDING 6

//timelinefx.slang mirrors these structs by hand, in its [Render_rows] section
static_assert(sizeof(tfx_instance_t) == 48, "tfx_instance_t changed, update tfx_instance in timelinefx.slang");
static_assert(sizeof(tfx_gpu_image_data_t) == 48, "tfx_gpu_image_data_t changed, update tfx_image_data in timelinefx.slang");
static_assert(sizeof(tfx_gpu_particle_properties_t) == 32, "tfx_gpu_particle_properties_t changed, update tfx_emitter_properties in timelinefx.slang");
static_assert(sizeof(tfx_ribbon_vertex_t) == 32, "tfx_ribbon_vertex_t changed, update tfx_ribbon_vertex_row in timelinefx.slang");
static_assert(sizeof(tfx_gpu_graph_data_t) == 48, "tfx_gpu_graph_data_t changed, update tfx_gpu_graph_data in timelinefx.slang");
static_assert(sizeof(tfx_impl_uniforms_t) % 16 == 0, "Uniform blocks must be a multiple of 16 bytes");

//Mirrors tfx_fixed_draw_constants in timelinefx.slang
typedef struct tfx_impl_gl_draw_constants_s {
	tfxU32 instance_offset;
	tfxU32 previous_index_offset;
	float ndc_offset_x;
	float ndc_offset_y;
} tfx_impl_gl_draw_constants_t;

//--- Gl_resources

#if defined(_WIN32) && !defined(_WIN64)
#define TFX_IMPL_GL_API __stdcall
#else
#define TFX_IMPL_GL_API
#endif

typedef unsigned int tfx_gl_enum;
typedef unsigned int tfx_gl_uint;
typedef int tfx_gl_int;
typedef int tfx_gl_sizei;
typedef unsigned int tfx_gl_bitfield;
typedef unsigned char tfx_gl_boolean;
typedef ptrdiff_t tfx_gl_intptr;
typedef ptrdiff_t tfx_gl_sizeiptr;
typedef void (TFX_IMPL_GL_API *tfx_gl_debug_callback)(tfx_gl_enum source, tfx_gl_enum type, tfx_gl_uint id, tfx_gl_enum severity, tfx_gl_sizei length, const char *message, const void *user_data);

#define TFX_GL_DEPTH_TEST 0x0B71
#define TFX_GL_BLEND 0x0BE2
#define TFX_GL_CULL_FACE 0x0B44
#define TFX_GL_LEQUAL 0x0203
#define TFX_GL_ONE 1
#define TFX_GL_ONE_MINUS_SRC_ALPHA 0x0303
#define TFX_GL_TRIANGLES 0x0004
#define TFX_GL_UNSIGNED_INT 0x1405
#define TFX_GL_UNSIGNED_BYTE 0x1401
#define TFX_GL_FLOAT 0x1406
#define TFX_GL_HALF_FLOAT 0x140B
#define TFX_GL_RGBA 0x1908
#define TFX_GL_RGBA8 0x8058
#define TFX_GL_RGBA16F 0x881A
#define TFX_GL_RGBA32F 0x8814
#define TFX_GL_TEXTURE_2D_ARRAY 0x8C1A
#define TFX_GL_TEXTURE0 0x84C0
#define TFX_GL_TEXTURE_MIN_FILTER 0x2801
#define TFX_GL_TEXTURE_MAG_FILTER 0x2800
#define TFX_GL_TEXTURE_WRAP_S 0x2802
#define TFX_GL_TEXTURE_WRAP_T 0x2803
#define TFX_GL_NEAREST 0x2600
#define TFX_GL_LINEAR 0x2601
#define TFX_GL_LINEAR_MIPMAP_LINEAR 0x2703
#define TFX_GL_CLAMP_TO_EDGE 0x812F
#define TFX_GL_ARRAY_BUFFER 0x8892
#define TFX_GL_ELEMENT_ARRAY_BUFFER 0x8893
#define TFX_GL_UNIFORM_BUFFER 0x8A11
#define TFX_GL_SHADER_STORAGE_BUFFER 0x90D2
#define TFX_GL_COPY_READ_BUFFER 0x8F36
#define TFX_GL_COPY_WRITE_BUFFER 0x8F37
#define TFX_GL_STATIC_DRAW 0x88E4
#define TFX_GL_DYNAMIC_DRAW 0x88E8
#define TFX_GL_MAP_WRITE_BIT 0x0002
#define TFX_GL_MAP_INVALIDATE_BUFFER_BIT 0x0008
#define TFX_GL_VERTEX_SHADER 0x8B31
#define TFX_GL_FRAGMENT_SHADER 0x8B30
#define TFX_GL_COMPUTE_SHADER 0x91B9
#define TFX_GL_COMPILE_STATUS 0x8B81
#define TFX_GL_LINK_STATUS 0x8B82
#define TFX_GL_INFO_LOG_LENGTH 0x8B84
#define TFX_GL_SHADER_STORAGE_BARRIER_BIT 0x00002000
#define TFX_GL_ELEMENT_ARRAY_BARRIER_BIT 0x00000002
#define TFX_GL_CONTEXT_FLAGS 0x821E
#define TFX_GL_CONTEXT_FLAG_DEBUG_BIT 0x00000002
#define TFX_GL_DEBUG_OUTPUT 0x92E0
#define TFX_GL_DEBUG_OUTPUT_SYNCHRONOUS 0x8242
#define TFX_GL_DEBUG_SEVERITY_NOTIFICATION 0x826B
#define TFX_GL_DONT_CARE 0x1100
#define TFX_GL_MAX_SHADER_STORAGE_BLOCK_SIZE 0x90DE
#define TFX_GL_MAJOR_VERSION 0x821B
#define TFX_GL_MINOR_VERSION 0x821C

//Every GL function the backend calls: return type, name, parameters
#define TFX_IMPL_GL_FUNCTIONS(X) \
	X(void, GetIntegerv, (tfx_gl_enum name, tfx_gl_int *data)) \
	X(void, Enable, (tfx_gl_enum capability)) \
	X(void, Disable, (tfx_gl_enum capability)) \
	X(void, BlendFunc, (tfx_gl_enum source, tfx_gl_enum destination)) \
	X(void, DepthFunc, (tfx_gl_enum function)) \
	X(void, DepthMask, (tfx_gl_boolean flag)) \
	X(void, GenBuffers, (tfx_gl_sizei count, tfx_gl_uint *buffers)) \
	X(void, DeleteBuffers, (tfx_gl_sizei count, const tfx_gl_uint *buffers)) \
	X(void, BindBuffer, (tfx_gl_enum target, tfx_gl_uint buffer)) \
	X(void, BindBufferBase, (tfx_gl_enum target, tfx_gl_uint index, tfx_gl_uint buffer)) \
	X(void, BufferData, (tfx_gl_enum target, tfx_gl_sizeiptr size, const void *data, tfx_gl_enum usage)) \
	X(void, BufferSubData, (tfx_gl_enum target, tfx_gl_intptr offset, tfx_gl_sizeiptr size, const void *data)) \
	X(void, CopyBufferSubData, (tfx_gl_enum read_target, tfx_gl_enum write_target, tfx_gl_intptr read_offset, tfx_gl_intptr write_offset, tfx_gl_sizeiptr size)) \
	X(void *, MapBufferRange, (tfx_gl_enum target, tfx_gl_intptr offset, tfx_gl_sizeiptr length, tfx_gl_bitfield access)) \
	X(tfx_gl_boolean, UnmapBuffer, (tfx_gl_enum target)) \
	X(void, GenTextures, (tfx_gl_sizei count, tfx_gl_uint *textures)) \
	X(void, DeleteTextures, (tfx_gl_sizei count, const tfx_gl_uint *textures)) \
	X(void, BindTexture, (tfx_gl_enum target, tfx_gl_uint texture)) \
	X(void, ActiveTexture, (tfx_gl_enum texture)) \
	X(void, TexParameteri, (tfx_gl_enum target, tfx_gl_enum name, tfx_gl_int value)) \
	X(void, TexStorage3D, (tfx_gl_enum target, tfx_gl_sizei levels, tfx_gl_enum internal_format, tfx_gl_sizei width, tfx_gl_sizei height, tfx_gl_sizei depth)) \
	X(void, TexSubImage3D, (tfx_gl_enum target, tfx_gl_int level, tfx_gl_int x, tfx_gl_int y, tfx_gl_int z, tfx_gl_sizei width, tfx_gl_sizei height, tfx_gl_sizei depth, tfx_gl_enum format, tfx_gl_enum type, const void *pixels)) \
	X(void, TextureView, (tfx_gl_uint texture, tfx_gl_enum target, tfx_gl_uint original, tfx_gl_enum internal_format, tfx_gl_uint min_level, tfx_gl_uint level_count, tfx_gl_uint min_layer, tfx_gl_uint layer_count)) \
	X(void, GenerateMipmap, (tfx_gl_enum target)) \
	X(tfx_gl_uint, CreateShader, (tfx_gl_enum type)) \
	X(void, ShaderSource, (tfx_gl_uint shader, tfx_gl_sizei count, const char *const *sources, const tfx_gl_int *lengths)) \
	X(void, CompileShader, (tfx_gl_uint shader)) \
	X(void, GetShaderiv, (tfx_gl_uint shader, tfx_gl_enum name, tfx_gl_int *value)) \
	X(void, GetShaderInfoLog, (tfx_gl_uint shader, tfx_gl_sizei size, tfx_gl_sizei *length, char *log)) \
	X(void, DeleteShader, (tfx_gl_uint shader)) \
	X(tfx_gl_uint, CreateProgram, (void)) \
	X(void, AttachShader, (tfx_gl_uint program, tfx_gl_uint shader)) \
	X(void, LinkProgram, (tfx_gl_uint program)) \
	X(void, GetProgramiv, (tfx_gl_uint program, tfx_gl_enum name, tfx_gl_int *value)) \
	X(void, GetProgramInfoLog, (tfx_gl_uint program, tfx_gl_sizei size, tfx_gl_sizei *length, char *log)) \
	X(void, DeleteProgram, (tfx_gl_uint program)) \
	X(void, UseProgram, (tfx_gl_uint program)) \
	X(void, GenVertexArrays, (tfx_gl_sizei count, tfx_gl_uint *arrays)) \
	X(void, DeleteVertexArrays, (tfx_gl_sizei count, const tfx_gl_uint *arrays)) \
	X(void, BindVertexArray, (tfx_gl_uint array)) \
	X(void, DrawArraysInstanced, (tfx_gl_enum mode, tfx_gl_int first, tfx_gl_sizei count, tfx_gl_sizei instance_count)) \
	X(void, DrawElements, (tfx_gl_enum mode, tfx_gl_sizei count, tfx_gl_enum type, const void *indices)) \
	X(void, DispatchCompute, (tfx_gl_uint groups_x, tfx_gl_uint groups_y, tfx_gl_uint groups_z)) \
	X(void, MemoryBarrier, (tfx_gl_bitfield barriers)) \
	X(void, DebugMessageCallback, (tfx_gl_debug_callback callback, const void *user_data)) \
	X(void, DebugMessageControl, (tfx_gl_enum source, tfx_gl_enum type, tfx_gl_enum severity, tfx_gl_sizei count, const tfx_gl_uint *ids, tfx_gl_boolean enabled))

#define TFX_IMPL_GL_DECLARE(return_type, name, parameters) return_type (TFX_IMPL_GL_API *name) parameters;
typedef struct tfx_impl_gl__functions_s {
	TFX_IMPL_GL_FUNCTIONS(TFX_IMPL_GL_DECLARE)
} tfx_impl_gl__functions_t;
#undef TFX_IMPL_GL_DECLARE

//One set for the process, because GL function pointers are the same for every context of a pixel format on every platform this runs on
static tfx_impl_gl__functions_t tfx_gl;

static bool tfx_impl_gl__load_functions(tfx_impl_gl_loader loader) {
	bool loaded = true;
#define TFX_IMPL_GL_LOAD(return_type, name, parameters) \
	tfx_gl.name = (return_type (TFX_IMPL_GL_API *) parameters)loader("gl" #name); \
	if (!tfx_gl.name) { tfx_Print("Could not load gl" #name ", the backend needs an OpenGL 4.3 core context"); loaded = false; }
	TFX_IMPL_GL_FUNCTIONS(TFX_IMPL_GL_LOAD)
#undef TFX_IMPL_GL_LOAD
	return loaded;
}

static void TFX_IMPL_GL_API tfx_impl_gl__debug_message(tfx_gl_enum source, tfx_gl_enum type, tfx_gl_uint id, tfx_gl_enum severity, tfx_gl_sizei length, const char *message, const void *user_data) {
	(void)source; (void)type; (void)severity; (void)length; (void)user_data;
	tfx_Print("OpenGL (%u): %s", id, message);
}

static void tfx_impl_gl__release_buffer(tfx_impl_gl_buffer_t *buffer) {
	if (buffer->buffer) {
		tfx_gl.DeleteBuffers(1, &buffer->buffer);
	}
	buffer->buffer = 0;
	buffer->capacity = 0;
}

//Every buffer is raw bytes in the C layout, so one kind of buffer serves as storage, uniform or element buffer
static bool tfx_impl_gl__create_buffer(const void *data, tfxU32 size_in_bytes, tfx_gl_enum usage, tfx_impl_gl_buffer_t *buffer) {
	tfx_impl_gl__release_buffer(buffer);
	size_in_bytes = (size_in_bytes + 15) & ~15u;
	tfx_gl.GenBuffers(1, &buffer->buffer);
	if (!buffer->buffer) {
		tfx_Print("Could not create a buffer of %u bytes", size_in_bytes);
		return false;
	}
	tfx_gl.BindBuffer(TFX_GL_COPY_WRITE_BUFFER, buffer->buffer);
	tfx_gl.BufferData(TFX_GL_COPY_WRITE_BUFFER, size_in_bytes, NULL, usage);
	if (data) {
		tfx_gl.BufferSubData(TFX_GL_COPY_WRITE_BUFFER, 0, size_in_bytes, data);
	}
	buffer->capacity = size_in_bytes;
	return true;
}

//Grows by half again whenever it's too small, so the size settles after the first few frames. Never empty, so a shader always has a buffer bound
static bool tfx_impl_gl__ensure_buffer(tfxU32 size_in_bytes, tfx_gl_enum usage, tfx_impl_gl_buffer_t *buffer) {
	size_in_bytes = size_in_bytes ? size_in_bytes : 16;
	if (buffer->buffer && buffer->capacity >= size_in_bytes) {
		return true;
	}
	return tfx_impl_gl__create_buffer(NULL, size_in_bytes + size_in_bytes / 2, usage, buffer);
}

static void tfx_impl_gl__update_buffer(tfx_impl_gl_buffer_t *buffer, const void *data, size_t size_in_bytes) {
	if (size_in_bytes == 0) {
		return;
	}
	tfx_gl.BindBuffer(TFX_GL_COPY_WRITE_BUFFER, buffer->buffer);
	tfx_gl.BufferSubData(TFX_GL_COPY_WRITE_BUFFER, 0, (tfx_gl_sizeiptr)size_in_bytes, data);
}

static bool tfx_impl_gl__write_buffer(tfx_impl_gl_buffer_t *buffer, const void *data, size_t size_in_bytes, tfx_gl_enum usage) {
	if (!tfx_impl_gl__ensure_buffer((tfxU32)size_in_bytes, usage, buffer)) {
		return false;
	}
	tfx_impl_gl__update_buffer(buffer, data, size_in_bytes);
	return true;
}

static void tfx_impl_gl__bind_storage(tfxU32 binding, const tfx_impl_gl_buffer_t *buffer) {
	tfx_gl.BindBufferBase(TFX_GL_SHADER_STORAGE_BUFFER, binding, buffer->buffer);
}

static void tfx_impl_gl__write_uniforms(tfx_impl_gl_buffer_t *buffer, tfxU32 binding, const void *data, size_t size) {
	tfx_impl_gl__update_buffer(buffer, data, size);
	tfx_gl.BindBufferBase(TFX_GL_UNIFORM_BUFFER, binding, buffer->buffer);
}

//The build writes each entry point to <entry point>.glsl. The caller frees the text
static char *tfx_impl_gl__load_shader_text(const char *shader_folder, const char *entry_point) {
	char path[1024];
	snprintf(path, sizeof(path), "%s/%s.glsl", shader_folder, entry_point);
	FILE *file = fopen(path, "rb");
	if (!file) {
		tfx_Print("Could not open the shader %s", path);
		return NULL;
	}
	fseek(file, 0, SEEK_END);
	long file_size = ftell(file);
	fseek(file, 0, SEEK_SET);
	char *text = file_size > 0 ? (char *)malloc((size_t)file_size + 1) : NULL;
	if (text && fread(text, 1, (size_t)file_size, file) == (size_t)file_size) {
		text[file_size] = '\0';
	} else {
		free(text);
		text = NULL;
		tfx_Print("Could not read the shader %s", path);
	}
	fclose(file);
	return text;
}

static tfx_gl_uint tfx_impl_gl__compile_shader(const char *shader_folder, const char *entry_point, tfx_gl_enum type) {
	char *text = tfx_impl_gl__load_shader_text(shader_folder, entry_point);
	if (!text) {
		return 0;
	}
	tfx_gl_uint shader = tfx_gl.CreateShader(type);
	tfx_gl.ShaderSource(shader, 1, &text, NULL);
	tfx_gl.CompileShader(shader);
	free(text);
	tfx_gl_int compiled = 0;
	tfx_gl.GetShaderiv(shader, TFX_GL_COMPILE_STATUS, &compiled);
	if (!compiled) {
		char log[2048] = {};
		tfx_gl.GetShaderInfoLog(shader, sizeof(log), NULL, log);
		tfx_Print("Could not compile %s: %s", entry_point, log);
		tfx_gl.DeleteShader(shader);
		return 0;
	}
	return shader;
}

static tfx_gl_uint tfx_impl_gl__link_program(const char *name, tfx_gl_uint first_shader, tfx_gl_uint second_shader) {
	tfx_gl_uint program = tfx_gl.CreateProgram();
	tfx_gl.AttachShader(program, first_shader);
	if (second_shader) {
		tfx_gl.AttachShader(program, second_shader);
	}
	tfx_gl.LinkProgram(program);
	tfx_gl_int linked = 0;
	tfx_gl.GetProgramiv(program, TFX_GL_LINK_STATUS, &linked);
	if (!linked) {
		char log[2048] = {};
		tfx_gl.GetProgramInfoLog(program, sizeof(log), NULL, log);
		tfx_Print("Could not link %s: %s", name, log);
		tfx_gl.DeleteProgram(program);
		return 0;
	}
	return program;
}

static tfx_gl_uint tfx_impl_gl__create_compute_program(const char *shader_folder, const char *entry_point) {
	tfx_gl_uint shader = tfx_impl_gl__compile_shader(shader_folder, entry_point, TFX_GL_COMPUTE_SHADER);
	if (!shader) {
		return 0;
	}
	tfx_gl_uint program = tfx_impl_gl__link_program(entry_point, shader, 0);
	tfx_gl.DeleteShader(shader);
	return program;
}

static tfx_gl_uint tfx_impl_gl__create_render_program(const char *shader_folder, const char *vertex_entry_point, tfx_gl_uint fragment_shader) {
	tfx_gl_uint vertex_shader = tfx_impl_gl__compile_shader(shader_folder, vertex_entry_point, TFX_GL_VERTEX_SHADER);
	if (!vertex_shader) {
		return 0;
	}
	tfx_gl_uint program = tfx_impl_gl__link_program(vertex_entry_point, vertex_shader, fragment_shader);
	tfx_gl.DeleteShader(vertex_shader);
	return program;
}

static void tfx_impl_gl__color_ramp_format(tfxU32 bytes_per_pixel, tfx_gl_enum *internal_format, tfx_gl_enum *type) {
	if (bytes_per_pixel == 4) {
		*internal_format = TFX_GL_RGBA8;
		*type = TFX_GL_UNSIGNED_BYTE;
	} else if (bytes_per_pixel == 16) {
		*internal_format = TFX_GL_RGBA32F;
		*type = TFX_GL_FLOAT;
	} else {
		*internal_format = TFX_GL_RGBA16F;
		*type = TFX_GL_HALF_FLOAT;
	}
}

//Read with texelFetch, but a texture still has to be complete for that, so the filters say there are no mips
static tfx_gl_uint tfx_impl_gl__create_color_ramp_texture(tfxU32 width, tfxU32 height, tfxU32 layer_count, tfxU32 bytes_per_pixel) {
	tfx_gl_enum internal_format, type;
	tfx_impl_gl__color_ramp_format(bytes_per_pixel, &internal_format, &type);
	tfx_gl_uint texture = 0;
	tfx_gl.GenTextures(1, &texture);
	tfx_gl.BindTexture(TFX_GL_TEXTURE_2D_ARRAY, texture);
	tfx_gl.TexStorage3D(TFX_GL_TEXTURE_2D_ARRAY, 1, internal_format, (tfx_gl_sizei)width, (tfx_gl_sizei)height, (tfx_gl_sizei)layer_count);
	tfx_gl.TexParameteri(TFX_GL_TEXTURE_2D_ARRAY, TFX_GL_TEXTURE_MIN_FILTER, TFX_GL_NEAREST);
	tfx_gl.TexParameteri(TFX_GL_TEXTURE_2D_ARRAY, TFX_GL_TEXTURE_MAG_FILTER, TFX_GL_NEAREST);
	return texture;
}

static void tfx_impl_gl__upload_color_ramp_layer(tfx_gl_uint texture, tfxU32 layer, tfx_bitmap_t *bitmap) {
	tfxU32 width = (tfxU32)tfx_GetBitmapWidth(bitmap);
	tfxU32 height = (tfxU32)tfx_GetBitmapHeight(bitmap);
	tfx_gl_enum internal_format, type;
	tfx_impl_gl__color_ramp_format((tfxU32)(tfx_GetBitmapSize(bitmap) / ((size_t)width * height)), &internal_format, &type);
	tfx_gl.BindTexture(TFX_GL_TEXTURE_2D_ARRAY, texture);
	tfx_gl.TexSubImage3D(TFX_GL_TEXTURE_2D_ARRAY, 0, 0, 0, (tfx_gl_int)layer, (tfx_gl_sizei)width, (tfx_gl_sizei)height, 1, TFX_GL_RGBA, type, tfx_GetBitmapData(bitmap));
}

//--- Gl_functions

bool tfx_impl_gl_Init(tfx_impl_gl_t *renderer, tfx_impl_gl_loader loader, const char *shader_folder) {
	memset(renderer, 0, sizeof(tfx_impl_gl_t));
	if (!tfx_impl_gl__load_functions(loader)) {
		return false;
	}
	tfx_gl_int major = 0, minor = 0, context_flags = 0, max_block_size = 0;
	tfx_gl.GetIntegerv(TFX_GL_MAJOR_VERSION, &major);
	tfx_gl.GetIntegerv(TFX_GL_MINOR_VERSION, &minor);
	if (major < 4 || (major == 4 && minor < 3)) {
		tfx_Print("The OpenGL context is %i.%i, the backend needs 4.3 for storage buffers and compute shaders", major, minor);
		return false;
	}
	//A debug context's messages go to the print callback. Notifications are left out because some drivers send one for every buffer
	tfx_gl.GetIntegerv(TFX_GL_CONTEXT_FLAGS, &context_flags);
	if (context_flags & TFX_GL_CONTEXT_FLAG_DEBUG_BIT) {
		tfx_gl.Enable(TFX_GL_DEBUG_OUTPUT);
		tfx_gl.Enable(TFX_GL_DEBUG_OUTPUT_SYNCHRONOUS);
		tfx_gl.DebugMessageCallback(tfx_impl_gl__debug_message, NULL);
		tfx_gl.DebugMessageControl(TFX_GL_DONT_CARE, TFX_GL_DONT_CARE, TFX_GL_DEBUG_SEVERITY_NOTIFICATION, 0, NULL, 0);
	}
	tfx_gl.GetIntegerv(TFX_GL_MAX_SHADER_STORAGE_BLOCK_SIZE, &max_block_size);
	renderer->max_storage_block_size = (tfxU32)max_block_size;

	//Each shader is compiled once and linked into every program that uses it, the live and baked billboards share their fragment shaders
	tfx_gl_uint billboard_vertex_shaders[2] = {};
	tfx_gl_uint ribbon_vertex_shaders[2] = {};
	tfx_gl_uint billboard_fragment_shaders[2][2] = {};
	tfx_gl_uint ribbon_fragment_shaders[2][2] = {};
	bool compiled = true;
	for (int lit = 0; lit != 2; ++lit) {
		billboard_vertex_shaders[lit] = tfx_impl_gl__compile_shader(shader_folder, tfx_impl_billboard_vertex_entry_points[lit], TFX_GL_VERTEX_SHADER);
		ribbon_vertex_shaders[lit] = tfx_impl_gl__compile_shader(shader_folder, tfx_impl_ribbon_vertex_entry_points[lit], TFX_GL_VERTEX_SHADER);
		compiled = compiled && billboard_vertex_shaders[lit] && ribbon_vertex_shaders[lit];
		for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
			billboard_fragment_shaders[writes_depth][lit] = tfx_impl_gl__compile_shader(shader_folder, tfx_impl_billboard_fragment_entry_points[writes_depth][lit], TFX_GL_FRAGMENT_SHADER);
			ribbon_fragment_shaders[writes_depth][lit] = tfx_impl_gl__compile_shader(shader_folder, tfx_impl_ribbon_fragment_entry_points[writes_depth][lit], TFX_GL_FRAGMENT_SHADER);
			compiled = compiled && billboard_fragment_shaders[writes_depth][lit] && ribbon_fragment_shaders[writes_depth][lit];
		}
	}
	bool linked = compiled;
	if (compiled) {
		for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
			for (int lit = 0; lit != 2; ++lit) {
				renderer->billboard_programs[writes_depth][lit] = tfx_impl_gl__link_program(tfx_impl_billboard_fragment_entry_points[writes_depth][lit],
					billboard_vertex_shaders[lit], billboard_fragment_shaders[writes_depth][lit]);
				renderer->ribbon_programs[writes_depth][lit] = tfx_impl_gl__link_program(tfx_impl_ribbon_fragment_entry_points[writes_depth][lit],
					ribbon_vertex_shaders[lit], ribbon_fragment_shaders[writes_depth][lit]);
				linked = linked && renderer->billboard_programs[writes_depth][lit] && renderer->ribbon_programs[writes_depth][lit];
			}
		}
		renderer->billboard_baked_program = tfx_impl_gl__create_render_program(shader_folder, "tfx_billboard_baked_vertex", billboard_fragment_shaders[0][0]);
	}
	for (int lit = 0; lit != 2; ++lit) {
		tfx_gl_uint shaders[6] = { billboard_vertex_shaders[lit], ribbon_vertex_shaders[lit], billboard_fragment_shaders[0][lit], billboard_fragment_shaders[1][lit],
			ribbon_fragment_shaders[0][lit], ribbon_fragment_shaders[1][lit] };
		for (int index = 0; index != 6; ++index) {
			if (shaders[index]) {
				tfx_gl.DeleteShader(shaders[index]);
			}
		}
	}
	renderer->particle_kernel = tfx_impl_gl__create_compute_program(shader_folder, "simulate_particles");
	renderer->ribbon_compute_program = tfx_impl_gl__create_compute_program(shader_folder, "tfx_ribbon_compute");
	renderer->sprite_playback_program = tfx_impl_gl__create_compute_program(shader_folder, "tfx_sprite_playback");
	renderer->ribbon_playback_program = tfx_impl_gl__create_compute_program(shader_folder, "tfx_ribbon_playback");
	if (!linked || !renderer->billboard_baked_program || !renderer->particle_kernel
		|| !renderer->ribbon_compute_program || !renderer->sprite_playback_program || !renderer->ribbon_playback_program) {
		return false;
	}

	size_t draw_block_size = sizeof(tfx_ribbon_bucket_globals_t) > sizeof(tfx_impl_gl_draw_constants_t) ? sizeof(tfx_ribbon_bucket_globals_t) : sizeof(tfx_impl_gl_draw_constants_t);
	size_t compute_block_size = sizeof(tfx_gpu_particle_push_t);
	compute_block_size = sizeof(tfx_ribbon_bucket_globals_t) > compute_block_size ? sizeof(tfx_ribbon_bucket_globals_t) : compute_block_size;
	compute_block_size = sizeof(tfx_sprite_data_push_t) > compute_block_size ? sizeof(tfx_sprite_data_push_t) : compute_block_size;
	compute_block_size = sizeof(tfx_ribbon_data_push_t) > compute_block_size ? sizeof(tfx_ribbon_data_push_t) : compute_block_size;
	tfx_gl.GenVertexArrays(1, &renderer->empty_vertex_array);
	return renderer->empty_vertex_array
		&& tfx_impl_gl__create_buffer(NULL, sizeof(tfx_impl_uniforms_t), TFX_GL_DYNAMIC_DRAW, &renderer->frame_constants)
		&& tfx_impl_gl__create_buffer(NULL, (tfxU32)draw_block_size, TFX_GL_DYNAMIC_DRAW, &renderer->draw_constants)
		&& tfx_impl_gl__create_buffer(NULL, (tfxU32)compute_block_size, TFX_GL_DYNAMIC_DRAW, &renderer->compute_constants);
}

bool tfx_impl_gl_CreateAtlas(tfx_impl_gl_t *renderer, tfxU32 layer_size, tfxU32 layer_count) {
	if (renderer->shape_atlas) {
		tfx_gl.DeleteTextures(1, &renderer->shape_atlas);
		renderer->shape_atlas = 0;
	}
	//Immutable storage, so tfx_impl_gl_UploadAtlas can make a view of a range of layers to generate their mips
	tfxU32 mip_count = 1;
	while ((layer_size >> mip_count) > 0) {
		mip_count++;
	}
	tfx_gl.GenTextures(1, &renderer->shape_atlas);
	tfx_gl.BindTexture(TFX_GL_TEXTURE_2D_ARRAY, renderer->shape_atlas);
	tfx_gl.TexStorage3D(TFX_GL_TEXTURE_2D_ARRAY, (tfx_gl_sizei)mip_count, TFX_GL_RGBA8, (tfx_gl_sizei)layer_size, (tfx_gl_sizei)layer_size, (tfx_gl_sizei)layer_count);
	tfx_gl.TexParameteri(TFX_GL_TEXTURE_2D_ARRAY, TFX_GL_TEXTURE_MIN_FILTER, TFX_GL_LINEAR_MIPMAP_LINEAR);
	tfx_gl.TexParameteri(TFX_GL_TEXTURE_2D_ARRAY, TFX_GL_TEXTURE_MAG_FILTER, TFX_GL_LINEAR);
	tfx_gl.TexParameteri(TFX_GL_TEXTURE_2D_ARRAY, TFX_GL_TEXTURE_WRAP_S, TFX_GL_CLAMP_TO_EDGE);
	tfx_gl.TexParameteri(TFX_GL_TEXTURE_2D_ARRAY, TFX_GL_TEXTURE_WRAP_T, TFX_GL_CLAMP_TO_EDGE);
	renderer->shape_atlas_layer_size = layer_size;
	renderer->shape_atlas_layer_count = layer_count;
	renderer->shape_atlas_mip_count = mip_count;
	return renderer->shape_atlas != 0;
}

bool tfx_impl_gl_UploadAtlas(tfx_impl_gl_t *renderer, const tfx_impl_atlas_t *atlas) {
	if (!renderer->shape_atlas || atlas->layer_size != renderer->shape_atlas_layer_size
		|| atlas->first_layer + atlas->layer_count > renderer->shape_atlas_layer_count) {
		tfx_Print("An atlas of %u layers from layer %u, %u pixels square, doesn't fit the %u layer, %u pixel texture array",
			atlas->layer_count, atlas->first_layer, atlas->layer_size, renderer->shape_atlas_layer_count, renderer->shape_atlas_layer_size);
		return false;
	}
	tfx_gl.BindTexture(TFX_GL_TEXTURE_2D_ARRAY, renderer->shape_atlas);
	tfx_gl.TexSubImage3D(TFX_GL_TEXTURE_2D_ARRAY, 0, 0, 0, (tfx_gl_int)atlas->first_layer, (tfx_gl_sizei)atlas->layer_size, (tfx_gl_sizei)atlas->layer_size,
		(tfx_gl_sizei)atlas->layer_count, TFX_GL_RGBA, TFX_GL_UNSIGNED_BYTE, atlas->pixels);
	//A view of just these layers, so the other atlases' mips aren't regenerated
	tfx_gl_uint range_view = 0;
	tfx_gl.GenTextures(1, &range_view);
	tfx_gl.TextureView(range_view, TFX_GL_TEXTURE_2D_ARRAY, renderer->shape_atlas, TFX_GL_RGBA8, 0, renderer->shape_atlas_mip_count, atlas->first_layer, atlas->layer_count);
	tfx_gl.BindTexture(TFX_GL_TEXTURE_2D_ARRAY, range_view);
	tfx_gl.GenerateMipmap(TFX_GL_TEXTURE_2D_ARRAY);
	tfx_gl.BindTexture(TFX_GL_TEXTURE_2D_ARRAY, 0);
	tfx_gl.DeleteTextures(1, &range_view);
	return true;
}

bool tfx_impl_gl_UploadSpriteData(tfx_impl_gl_t *renderer, tfx_animation_manager animation_manager, tfx_gpu_shapes shapes) {
	if (renderer->baked_color_ramps) {
		tfx_gl.DeleteTextures(1, &renderer->baked_color_ramps);
		renderer->baked_color_ramps = 0;
	}
	tfxU32 bitmap_count = tfx_GetAnimationColorRampBitmapCount(animation_manager);
	if (bitmap_count) {
		tfx_bitmap_t *first_bitmap = tfx_GetAnimationColorRampBitmap(animation_manager, 0);
		tfxU32 width = (tfxU32)tfx_GetBitmapWidth(first_bitmap);
		tfxU32 height = (tfxU32)tfx_GetBitmapHeight(first_bitmap);
		renderer->baked_color_ramps = tfx_impl_gl__create_color_ramp_texture(width, height, bitmap_count, (tfxU32)(tfx_GetBitmapSize(first_bitmap) / ((size_t)width * height)));
		for (tfxU32 index = 0; index != bitmap_count; ++index) {
			tfx_impl_gl__upload_color_ramp_layer(renderer->baked_color_ramps, index, tfx_GetAnimationColorRampBitmap(animation_manager, index));
		}
	}

	tfxU32 shape_count = tfx_GetGPUShapesCount(shapes);
	tfx_gpu_image_data_t empty_shape = {};
	tfxU32 property_count = (tfxU32)(tfx_GetAnimationEmitterPropertySizeInBytes(animation_manager) / sizeof(tfx_gpu_particle_properties_t));
	tfx_gpu_particle_properties_t empty_properties = {};
	size_t sprite_bytes = tfx_GetSpriteDataSizeInBytes(animation_manager);
	tfx_instance_t empty_sprite = {};
	if (!tfx_impl_gl__create_buffer(shape_count ? (const void *)tfx_GetGPUShapesArray(shapes) : (const void *)&empty_shape,
			(tfxU32)sizeof(tfx_gpu_image_data_t) * (shape_count ? shape_count : 1), TFX_GL_STATIC_DRAW, &renderer->baked_image_data)
		|| !tfx_impl_gl__create_buffer(property_count ? tfx_GetAnimationEmitterPropertiesBufferPointer(animation_manager) : (const void *)&empty_properties,
			(tfxU32)sizeof(tfx_gpu_particle_properties_t) * (property_count ? property_count : 1), TFX_GL_STATIC_DRAW, &renderer->baked_emitter_properties)
		|| !tfx_impl_gl__create_buffer(sprite_bytes ? tfx_GetSpriteDataBufferPointer(animation_manager) : (const void *)&empty_sprite,
			(tfxU32)(sprite_bytes ? sprite_bytes : sizeof(tfx_instance_t)), TFX_GL_STATIC_DRAW, &renderer->baked_sprites)) {
		return false;
	}
	if (!tfx_AnimationManagerHasRibbons(animation_manager)) {
		return true;
	}

	//The baked ribbons, their properties, segments, emitters and graph lookups only change when sprite data is loaded
	tfxU32 graph_bytes = (tfxU32)tfx_GetAnimationRibbonGraphLookupsSizeInBytes(animation_manager);
	tfx_gpu_graph_data_t empty_graph = {};
	if (!tfx_impl_gl__create_buffer(graph_bytes ? tfx_GetAnimationRibbonGraphLookupsBufferPointer(animation_manager) : (const void *)&empty_graph,
		graph_bytes ? graph_bytes : sizeof(tfx_gpu_graph_data_t), TFX_GL_STATIC_DRAW, &renderer->baked_ribbons.graph_lookups)) {
		return false;
	}
	size_t emitters_bytes = tfx_GetAnimationRibbonEmittersSizeInBytes(animation_manager);
	void *emitters = malloc(emitters_bytes);
	if (!emitters) {
		return false;
	}
	tfx_CopyAnimationRibbonEmitters(animation_manager, emitters);
	bool created = tfx_impl_gl__create_buffer(tfx_GetAnimationRibbonDataBufferPointer(animation_manager),
			(tfxU32)tfx_GetAnimationRibbonDataSizeInBytes(animation_manager), TFX_GL_STATIC_DRAW, &renderer->baked_ribbon_data)
		&& tfx_impl_gl__create_buffer(tfx_GetAnimationRibbonPropertiesBufferPointer(animation_manager),
			(tfxU32)tfx_GetAnimationRibbonPropertySizeInBytes(animation_manager), TFX_GL_STATIC_DRAW, &renderer->baked_ribbon_properties)
		&& tfx_impl_gl__create_buffer(tfx_GetAnimationRibbonSegmentsBufferPointer(animation_manager),
			(tfxU32)tfx_GetAnimationRibbonSegmentsSizeInBytes(animation_manager), TFX_GL_STATIC_DRAW, &renderer->baked_ribbons.segments)
		&& tfx_impl_gl__create_buffer(emitters, (tfxU32)emitters_bytes, TFX_GL_STATIC_DRAW, &renderer->baked_ribbons.emitters);
	free(emitters);
	return created;
}

//The vertex and index buffers ribbon_compute writes, sized for this frame's buckets
static bool tfx_impl_gl__ensure_ribbon_geometry(const tfx_ribbon_batch_sizes_t *sizes, tfx_impl_gl_ribbon_set_t *set) {
	return tfx_impl_gl__ensure_buffer((tfxU32)sizes->vertex_buffer_size_in_bytes, TFX_GL_DYNAMIC_DRAW, &set->vertices)
		&& tfx_impl_gl__ensure_buffer((tfxU32)sizes->index_buffer_size_in_bytes, TFX_GL_DYNAMIC_DRAW, &set->indices);
}

static bool tfx_impl_gl__upload_baked_ribbons(tfx_impl_gl_t *renderer, const tfx_impl_baked_frame_t *baked) {
	if (baked->ribbon_count == 0 || !renderer->baked_ribbon_data.buffer) {
		return true;
	}
	tfx_ribbon_batch_sizes_t sizes = tfx_GetAnimationRibbonBufferSizes(baked->animation_manager, sizeof(tfx_ribbon_vertex_t));
	if (!tfx_impl_gl__write_buffer(&renderer->animation_ribbon_offsets, baked->ribbon_offsets, baked->ribbon_offsets_bytes, TFX_GL_DYNAMIC_DRAW)
		|| !tfx_impl_gl__ensure_buffer(baked->ribbon_count * tfxRIBBON_SIZE, TFX_GL_DYNAMIC_DRAW, &renderer->baked_ribbons.instances)
		|| !tfx_impl_gl__ensure_ribbon_geometry(&sizes, &renderer->baked_ribbons)) {
		return false;
	}
	memset(&renderer->ribbon_data_push, 0, sizeof(tfx_ribbon_data_push_t));
	renderer->ribbon_data_push.animation_instances_total = baked->animation_instance_count;
	renderer->ribbon_data_push.ribbons_total = baked->ribbon_count;
	renderer->playback_ribbon_count = baked->ribbon_count;
	renderer->baked_ribbons.animation_manager = baked->animation_manager;
	return true;
}

//Animation instances advance every frame, so this runs every frame, not just on updates
static bool tfx_impl_gl__upload_baked(tfx_impl_gl_t *renderer, const tfx_impl_baked_frame_t *baked) {
	renderer->playback_sprite_count = 0;
	renderer->playback_ribbon_count = 0;
	renderer->baked_ribbons.animation_manager = NULL;
	if (!renderer->baked_sprites.buffer || (baked->sprite_count == 0 && baked->ribbon_count == 0)) {
		return true;
	}
	if (!tfx_impl_gl__write_buffer(&renderer->animation_instances, baked->animation_instances, baked->animation_instances_bytes, TFX_GL_DYNAMIC_DRAW)) {
		return false;
	}
	if (baked->sprite_count > 0) {
		if (!tfx_impl_gl__write_buffer(&renderer->animation_offsets, baked->offsets, baked->offsets_bytes, TFX_GL_DYNAMIC_DRAW)
			|| !tfx_impl_gl__ensure_buffer(baked->sprite_count * sizeof(tfx_instance_t), TFX_GL_DYNAMIC_DRAW, &renderer->playback_instances)) {
			return false;
		}
		memset(&renderer->sprite_data_push, 0, sizeof(tfx_sprite_data_push_t));
		renderer->sprite_data_push.animation_instances_total = baked->animation_instance_count;
		renderer->sprite_data_push.billboards_total = baked->sprite_count;
		renderer->sprite_data_push.animated_shapes = baked->animated_shapes ? 1 : 0;
		renderer->playback_sprite_count = baked->sprite_count;
	}
	return tfx_impl_gl__upload_baked_ribbons(renderer, baked);
}

static bool tfx_impl_gl__upload_ribbons(tfx_impl_gl_t *renderer, const tfx_impl_ribbon_frame_t *ribbons) {
	tfx_impl_gl_ribbon_set_t *set = &renderer->live_ribbons;
	set->stage = ribbons->stage;
	if (!ribbons->stage) {
		return true;
	}
	const tfx_ribbon_batch_sizes_t *sizes = &ribbons->sizes;
	if (!tfx_impl_gl__write_buffer(&set->segments, ribbons->segments, sizes->segment_buffer_size_in_bytes, TFX_GL_DYNAMIC_DRAW)
		|| !tfx_impl_gl__write_buffer(&set->instances, ribbons->ribbons, sizes->ribbon_buffer_size_in_bytes, TFX_GL_DYNAMIC_DRAW)
		|| !tfx_impl_gl__write_buffer(&set->emitters, ribbons->emitters, sizes->emitter_buffer_size_in_bytes, TFX_GL_DYNAMIC_DRAW)
		|| !tfx_impl_gl__ensure_ribbon_geometry(sizes, set)) {
		set->stage = NULL;
		return false;
	}
	set->stage_offsets = ribbons->stage_offsets;
	return true;
}

//The old contents are copied across, because the previous buffer is still what captured_index points into
static bool tfx_impl_gl__grow_instance_buffers(tfx_impl_gl_t *renderer, tfxU32 needed) {
	tfxU32 capacity = renderer->instance_capacity ? renderer->instance_capacity : TFX_IMPL_GL_INITIAL_INSTANCE_CAPACITY;
	while (capacity < needed) {
		capacity += capacity / 2;
	}
	for (int index = 0; index != 2; ++index) {
		tfx_impl_gl_buffer_t grown = {};
		if (!tfx_impl_gl__create_buffer(NULL, capacity * sizeof(tfx_instance_t), TFX_GL_DYNAMIC_DRAW, &grown)) {
			return false;
		}
		tfx_impl_gl_buffer_t *old = &renderer->instance_buffers[index];
		if (old->buffer) {
			tfx_gl.BindBuffer(TFX_GL_COPY_READ_BUFFER, old->buffer);
			tfx_gl.BindBuffer(TFX_GL_COPY_WRITE_BUFFER, grown.buffer);
			tfx_gl.CopyBufferSubData(TFX_GL_COPY_READ_BUFFER, TFX_GL_COPY_WRITE_BUFFER, 0, 0, old->capacity);
		}
		tfx_impl_gl__release_buffer(old);
		*old = grown;
	}
	renderer->instance_capacity = capacity;
	return true;
}

static bool tfx_impl_gl__upload_gpu_particles(tfx_impl_gl_t *renderer, const tfx_impl_frame_t *frame) {
	const tfx_gpu_particle_frame_t *gpu_frame = &frame->gpu_particles;
	//The pool holds the particles between frames, it's sized once and starts zeroed
	if (!renderer->particle_pool.buffer) {
		if (frame->gpu_pool_sizes.pool_bytes > renderer->max_storage_block_size) {
			tfx_Print("The GPU particle pool is %zu bytes but this driver's GL_MAX_SHADER_STORAGE_BLOCK_SIZE is %u. Lower gpu_particle_capacity",
				frame->gpu_pool_sizes.pool_bytes, renderer->max_storage_block_size);
			return false;
		}
		void *zeroes = calloc(1, frame->gpu_pool_sizes.pool_bytes);
		bool created = zeroes && tfx_impl_gl__create_buffer(zeroes, (tfxU32)frame->gpu_pool_sizes.pool_bytes, TFX_GL_STATIC_DRAW, &renderer->particle_pool);
		free(zeroes);
		if (!created) {
			return false;
		}
	}
	if (!tfx_impl_gl__ensure_buffer((tfxU32)gpu_frame->upload_bytes, TFX_GL_DYNAMIC_DRAW, &renderer->particle_upload)) {
		return false;
	}
	tfx_gl.BindBuffer(TFX_GL_COPY_WRITE_BUFFER, renderer->particle_upload.buffer);
	void *mapped = tfx_gl.MapBufferRange(TFX_GL_COPY_WRITE_BUFFER, 0, renderer->particle_upload.capacity, TFX_GL_MAP_WRITE_BIT | TFX_GL_MAP_INVALIDATE_BUFFER_BIT);
	if (!mapped) {
		return false;
	}
	tfx_CopyGPUParticleUpload(gpu_frame, mapped);
	tfx_gl.UnmapBuffer(TFX_GL_COPY_WRITE_BUFFER);

	//With fixed bindings the buffer indexes are unused, the bases say where each section starts in the pool and the upload
	tfx_gpu_particle_push_t *push = &renderer->particle_push_constants;
	memset(push, 0, sizeof(tfx_gpu_particle_push_t));
	push->cold_base = frame->gpu_pool_sizes.cold_base;
	push->noise_base = frame->gpu_pool_sizes.noise_base;
	push->tick_base = gpu_frame->tick_base;
	push->spawn_base = gpu_frame->spawn_base;
	push->location_base = gpu_frame->location_base;
	renderer->particle_dispatches = gpu_frame->dispatches;
	renderer->particle_dispatch_count = gpu_frame->dispatch_count;
	return true;
}

//Global layers never move, so only the layers whose version changed go up, unless the texture array has to grow
static bool tfx_impl_gl__upload_color_ramps(tfx_impl_gl_t *renderer) {
	tfxU32 layer_count = tfx_GetGlobalColorRampLayerCount();
	tfx_bitmap_t *first_bitmap = NULL;
	for (tfxU32 layer = 0; layer != layer_count && !first_bitmap; ++layer) {
		first_bitmap = tfx_GetGlobalColorRampLayer(layer);
	}
	if (!first_bitmap) {
		return true;
	}
	bool recreate = layer_count > renderer->color_ramp_layer_capacity;
	if (recreate) {
		//Headroom so that loading another library doesn't always recreate it
		tfxU32 capacity = layer_count + layer_count / 2 + 1;
		capacity = capacity < TFX_IMPL_GL_MAX_COLOR_RAMP_LAYERS ? capacity : TFX_IMPL_GL_MAX_COLOR_RAMP_LAYERS;
		tfxU32 *layer_versions = (tfxU32 *)realloc(renderer->color_ramp_layer_versions, sizeof(tfxU32) * capacity);
		if (!layer_versions) {
			return false;
		}
		renderer->color_ramp_layer_versions = layer_versions;
		if (renderer->color_ramps) {
			tfx_gl.DeleteTextures(1, &renderer->color_ramps);
		}
		tfxU32 width = (tfxU32)tfx_GetBitmapWidth(first_bitmap);
		tfxU32 height = (tfxU32)tfx_GetBitmapHeight(first_bitmap);
		renderer->color_ramps = tfx_impl_gl__create_color_ramp_texture(width, height, capacity, (tfxU32)(tfx_GetBitmapSize(first_bitmap) / ((size_t)width * height)));
		renderer->color_ramp_layer_capacity = capacity;
	}
	for (tfxU32 layer = 0; layer != layer_count; ++layer) {
		tfxU32 version = tfx_GetGlobalColorRampLayerVersion(layer);
		if (!recreate && version == renderer->color_ramp_layer_versions[layer]) {
			continue;
		}
		renderer->color_ramp_layer_versions[layer] = version;
		//A free layer isn't drawn from, so whatever it holds can stay
		tfx_bitmap_t *bitmap = tfx_GetGlobalColorRampLayer(layer);
		if (bitmap) {
			tfx_impl_gl__upload_color_ramp_layer(renderer->color_ramps, layer, bitmap);
		}
	}
	return true;
}

//Only on update frames: the instances drawn until then were written with the layout the last upload holds
static bool tfx_impl_gl__upload_library_data(tfx_impl_gl_t *renderer) {
	tfxU32 properties_version = tfx_GetGlobalParticlePropertiesVersion();
	if (!renderer->emitter_properties.buffer || properties_version != renderer->properties_version) {
		if (!tfx_impl_gl__write_buffer(&renderer->emitter_properties, tfx_GetGlobalParticlePropertiesBuffer(), tfx_GetGlobalParticlePropertiesBufferSizeInBytes(), TFX_GL_STATIC_DRAW)) {
			return false;
		}
		renderer->properties_version = properties_version;
	}
	tfxU32 shapes_version = tfx_GetGlobalGPUShapesVersion();
	if (!renderer->image_data.buffer || shapes_version != renderer->shapes_version) {
		tfx_gpu_shapes shapes = tfx_GetGlobalGPUShapes();
		if (!tfx_impl_gl__write_buffer(&renderer->image_data, tfx_GetGPUShapesArray(shapes), sizeof(tfx_gpu_image_data_t) * tfx_GetGPUShapesCount(shapes), TFX_GL_STATIC_DRAW)) {
			return false;
		}
		renderer->shapes_version = shapes_version;
	}
	//Shared by every library already, read by the ribbon shaders and the GPU kernel
	tfxU32 graph_lookups_version = tfx_GetGPUGraphLookupsVersion();
	if (!renderer->graph_lookups.buffer || graph_lookups_version != renderer->graph_lookups_version) {
		tfxU32 graph_bytes = tfx_GetGPUGraphLookupsBufferSizeInBytes();
		if (!tfx_impl_gl__ensure_buffer(graph_bytes ? graph_bytes : sizeof(tfx_gpu_graph_data_t), TFX_GL_STATIC_DRAW, &renderer->graph_lookups)) {
			return false;
		}
		tfx_impl_gl__update_buffer(&renderer->graph_lookups, tfx_GetGPUGraphLookupsBuffer(), graph_bytes);
		renderer->graph_lookups_version = graph_lookups_version;
	}
	tfxU32 color_ramps_version = tfx_GetGlobalColorRampsVersion();
	if (!renderer->color_ramps || color_ramps_version != renderer->color_ramps_version) {
		if (!tfx_impl_gl__upload_color_ramps(renderer)) {
			return false;
		}
		renderer->color_ramps_version = color_ramps_version;
	}
	return true;
}

bool tfx_impl_gl_UploadFrame(tfx_impl_gl_t *renderer, const tfx_impl_frame_t *frame) {
	memcpy(renderer->layer_flags, frame->layer_flags, sizeof(renderer->layer_flags));
	renderer->library_layer_count = frame->library_layer_count;
	if (!tfx_impl_gl__upload_baked(renderer, &frame->baked)) {
		return false;
	}
	if (!frame->stage_updated) {
		return true;
	}
	if (!tfx_impl_gl__upload_library_data(renderer)) {
		return false;
	}
	bool gpu_particles = frame->has_gpu_particles && renderer->particle_kernel;
	tfxU32 needed = gpu_particles ? frame->gpu_particles.gpu_instance_start + frame->gpu_particles.gpu_instance_count : frame->instance_count;
	if ((needed > renderer->instance_capacity || renderer->instance_capacity == 0) && !tfx_impl_gl__grow_instance_buffers(renderer, needed)) {
		renderer->instance_count = 0;
		return false;
	}

	//The other buffer now holds what was drawn last, which is what captured_index points into for CPU and GPU instances alike
	renderer->current_instance_buffer ^= 1;
	tfx_impl_gl__update_buffer(&renderer->instance_buffers[renderer->current_instance_buffer], frame->instances, frame->instance_count * sizeof(tfx_instance_t));
	renderer->instance_count = needed;
	memcpy(renderer->layer_start, frame->layer_start, sizeof(renderer->layer_start));
	memcpy(renderer->layer_count, frame->layer_count, sizeof(renderer->layer_count));
	memset(renderer->gpu_layer_count, 0, sizeof(renderer->gpu_layer_count));
	if (gpu_particles) {
		memcpy(renderer->gpu_layer_start, frame->gpu_particles.gpu_layer_start, sizeof(renderer->gpu_layer_start));
		memcpy(renderer->gpu_layer_count, frame->gpu_particles.gpu_layer_count, sizeof(renderer->gpu_layer_count));
		if (!tfx_impl_gl__upload_gpu_particles(renderer, frame)) {
			return false;
		}
	}
	return tfx_impl_gl__upload_ribbons(renderer, &frame->ribbons);
}

//Every fetch has to be dispatched, the next one carries on from the particle state this one leaves in the pool. The dispatches
//write separate rings and instance ranges, so they need no barrier between them
static void tfx_impl_gl__simulate_gpu_particles(tfx_impl_gl_t *renderer) {
	tfx_gl.UseProgram(renderer->particle_kernel);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_WRITE_BINDING + 0, &renderer->particle_pool);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_WRITE_BINDING + 1, &renderer->instance_buffers[renderer->current_instance_buffer]);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_READ_BINDING + 0, &renderer->graph_lookups);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_READ_BINDING + 1, &renderer->particle_upload);
	tfx_gpu_particle_push_t *push = &renderer->particle_push_constants;
	for (tfxU32 index = 0; index != renderer->particle_dispatch_count; ++index) {
		tfxU32 slot_count = renderer->particle_dispatches[index].slot_count;
		if (slot_count == 0) {
			continue;
		}
		push->dispatch_index = index;
		tfx_impl_gl__write_uniforms(&renderer->compute_constants, TFX_IMPL_GL_COMPUTE_UNIFORM_BINDING, push, sizeof(tfx_gpu_particle_push_t));
		tfx_gl.DispatchCompute((slot_count + TFX_IMPL_GL_KERNEL_GROUP_SIZE - 1) / TFX_IMPL_GL_KERNEL_GROUP_SIZE, 1, 1);
	}
	renderer->particle_dispatch_count = 0;
}

static void tfx_impl_gl__play_back_sprites(tfx_impl_gl_t *renderer) {
	tfx_gl.UseProgram(renderer->sprite_playback_program);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_READ_BINDING + 0, &renderer->animation_offsets);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_READ_BINDING + 1, &renderer->animation_instances);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_READ_BINDING + 2, &renderer->baked_sprites);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_READ_BINDING + 3, &renderer->baked_emitter_properties);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_WRITE_BINDING + 0, &renderer->playback_instances);
	tfx_impl_gl__write_uniforms(&renderer->compute_constants, TFX_IMPL_GL_COMPUTE_UNIFORM_BINDING, &renderer->sprite_data_push, sizeof(tfx_sprite_data_push_t));
	tfx_gl.DispatchCompute((renderer->playback_sprite_count + TFX_IMPL_GL_PLAYBACK_GROUP_SIZE - 1) / TFX_IMPL_GL_PLAYBACK_GROUP_SIZE, 1, 1);
}

static void tfx_impl_gl__play_back_ribbons(tfx_impl_gl_t *renderer) {
	tfx_gl.UseProgram(renderer->ribbon_playback_program);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_READ_BINDING + 0, &renderer->animation_ribbon_offsets);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_READ_BINDING + 1, &renderer->animation_instances);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_READ_BINDING + 2, &renderer->baked_ribbon_data);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_READ_BINDING + 3, &renderer->baked_ribbon_properties);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_WRITE_BINDING + 0, &renderer->baked_ribbons.instances);
	tfx_impl_gl__write_uniforms(&renderer->compute_constants, TFX_IMPL_GL_COMPUTE_UNIFORM_BINDING, &renderer->ribbon_data_push, sizeof(tfx_ribbon_data_push_t));
	tfx_gl.DispatchCompute((renderer->playback_ribbon_count + TFX_IMPL_GL_PLAYBACK_GROUP_SIZE - 1) / TFX_IMPL_GL_PLAYBACK_GROUP_SIZE, 1, 1);
	//ribbon_compute reads the ribbons this wrote
	tfx_gl.MemoryBarrier(TFX_GL_SHADER_STORAGE_BARRIER_BIT);
}

static bool tfx_impl_gl__ribbon_set_active(const tfx_impl_gl_ribbon_set_t *set) {
	return set->stage || set->animation_manager;
}

static tfx_ribbon_dispatch_t tfx_impl_gl__create_ribbon_dispatch(const tfx_impl_gl_ribbon_set_t *set) {
	return set->stage ? tfx_CreateRibbonBatchDispatch(&set->stage_offsets) : tfx_CreateAnimationRibbonDispatch();
}

static bool tfx_impl_gl__next_ribbon_dispatch(const tfx_impl_gl_ribbon_set_t *set, tfx_ribbon_dispatch_t *dispatch) {
	return set->stage ? tfx_NextRibbonDispatch(set->stage, dispatch) : tfx_NextAnimationRibbonDispatch(set->animation_manager, dispatch);
}

static const tfx_impl_gl_buffer_t *tfx_impl_gl__ribbon_graph_lookups(const tfx_impl_gl_t *renderer, const tfx_impl_gl_ribbon_set_t *set) {
	return set->stage ? &renderer->graph_lookups : &set->graph_lookups;
}

static void tfx_impl_gl__build_ribbon_geometry(tfx_impl_gl_t *renderer, const tfx_impl_gl_ribbon_set_t *set, const tfx_impl_uniforms_t *uniforms) {
	tfx_gl.UseProgram(renderer->ribbon_compute_program);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_READ_BINDING + 0, &set->emitters);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_READ_BINDING + 1, tfx_impl_gl__ribbon_graph_lookups(renderer, set));
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_READ_BINDING + 2, &set->instances);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_READ_BINDING + 3, &set->segments);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_WRITE_BINDING + 0, &set->vertices);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_COMPUTE_WRITE_BINDING + 1, &set->indices);

	//One dispatch per bucket, each with the bucket's own globals
	tfx_ribbon_dispatch_t dispatch = tfx_impl_gl__create_ribbon_dispatch(set);
	while (tfx_impl_gl__next_ribbon_dispatch(set, &dispatch)) {
		tfx_ribbon_bucket_globals_t *globals = tfx_GetRibbonDispatchGlobals(&dispatch);
		if (set->stage) {
			globals->lerp = uniforms->timer_lerp;
		} else {
			//Baked ribbons are already interpolated, but they still face the camera
			globals->camera_position.x = uniforms->camera_position[0];
			globals->camera_position.y = uniforms->camera_position[1];
			globals->camera_position.z = uniforms->camera_position[2];
		}
		tfx_impl_gl__write_uniforms(&renderer->compute_constants, TFX_IMPL_GL_COMPUTE_UNIFORM_BINDING, globals, sizeof(tfx_ribbon_bucket_globals_t));
		tfx_gl.DispatchCompute((dispatch.total_segments + TFX_IMPL_GL_RIBBON_GROUP_SIZE - 1) / TFX_IMPL_GL_RIBBON_GROUP_SIZE, 1, 1);
	}
}

void tfx_impl_gl_Compute(tfx_impl_gl_t *renderer, const tfx_impl_uniforms_t *uniforms) {
	bool dispatched = renderer->particle_dispatch_count > 0 || renderer->playback_sprite_count > 0;
	if (renderer->particle_dispatch_count > 0) {
		tfx_impl_gl__simulate_gpu_particles(renderer);
	}
	if (renderer->playback_sprite_count > 0) {
		tfx_impl_gl__play_back_sprites(renderer);
	}
	if (renderer->playback_ribbon_count > 0) {
		tfx_impl_gl__play_back_ribbons(renderer);
	}
	if (tfx_impl_gl__ribbon_set_active(&renderer->live_ribbons)) {
		tfx_impl_gl__build_ribbon_geometry(renderer, &renderer->live_ribbons, uniforms);
		dispatched = true;
	}
	if (tfx_impl_gl__ribbon_set_active(&renderer->baked_ribbons)) {
		tfx_impl_gl__build_ribbon_geometry(renderer, &renderer->baked_ribbons, uniforms);
		dispatched = true;
	}
	//Storage writes aren't seen by later draws, or by the element fetch for the ribbon indexes, without a barrier
	if (dispatched) {
		tfx_gl.MemoryBarrier(TFX_GL_SHADER_STORAGE_BARRIER_BIT | TFX_GL_ELEMENT_ARRAY_BARRIER_BIT);
	}
	tfx_gl.UseProgram(0);
}

//Live and baked ribbons each have their own shape frames and color ramps, the atlas is shared
//Every bucket is walked, because the dispatch carries the offsets from one to the next, and only the ones whose layer is in the group are drawn
static void tfx_impl_gl__render_ribbons(tfx_impl_gl_t *renderer, const tfx_impl_gl_ribbon_set_t *set, const tfx_impl_gl_buffer_t *image_data, tfx_gl_uint color_ramps, int writes_depth) {
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_VERTEX_READ_BINDING + 0, &set->vertices);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_VERTEX_READ_BINDING + 1, &set->instances);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_VERTEX_READ_BINDING + 2, &set->emitters);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_VERTEX_READ_BINDING + 3, tfx_impl_gl__ribbon_graph_lookups(renderer, set));
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_VERTEX_READ_BINDING + 4, image_data);
	tfx_gl.ActiveTexture(TFX_GL_TEXTURE0 + TFX_IMPL_GL_FRAGMENT_TEXTURE_BINDING + 1);
	tfx_gl.BindTexture(TFX_GL_TEXTURE_2D_ARRAY, color_ramps);
	//The element buffer binding belongs to the vertex array
	tfx_gl.BindBuffer(TFX_GL_ELEMENT_ARRAY_BUFFER, set->indices.buffer);

	tfx_gl_uint bound_program = 0;
	tfx_ribbon_dispatch_t dispatch = tfx_impl_gl__create_ribbon_dispatch(set);
	while (tfx_impl_gl__next_ribbon_dispatch(set, &dispatch)) {
		const tfxLayerFlags flags = renderer->layer_flags[dispatch.layer];
		if (tfx_impl_LayerWritesDepth(flags) != writes_depth) {
			continue;
		}
		tfx_gl_uint program = renderer->ribbon_programs[writes_depth][tfx_impl_LayerIsLit(flags)];
		if (program != bound_program) {
			tfx_gl.UseProgram(program);
			bound_program = program;
		}
		tfx_impl_gl__write_uniforms(&renderer->draw_constants, TFX_IMPL_GL_VERTEX_UNIFORM_BINDING + 1, tfx_GetRibbonDispatchGlobals(&dispatch), sizeof(tfx_ribbon_bucket_globals_t));
		tfx_gl.DrawElements(TFX_GL_TRIANGLES, (tfx_gl_sizei)dispatch.index_count, TFX_GL_UNSIGNED_INT, (const void *)((uintptr_t)dispatch.index_offset * sizeof(tfxU32)));
	}
	tfx_gl.BindBuffer(TFX_GL_ELEMENT_ARRAY_BUFFER, 0);
}

//Every draw starts at instance 0 and the shader adds instance_offset, see the OpenGL notes in Shaders/readme.md
static void tfx_impl_gl__draw_instance_range(tfx_impl_gl_t *renderer, tfxU32 start, tfxU32 count) {
	if (count == 0) {
		return;
	}
	tfx_impl_gl_draw_constants_t draw_constants = {};
	draw_constants.instance_offset = start;
	tfx_impl_gl__write_uniforms(&renderer->draw_constants, TFX_IMPL_GL_VERTEX_UNIFORM_BINDING + 1, &draw_constants, sizeof(draw_constants));
	tfx_gl.DrawArraysInstanced(TFX_GL_TRIANGLES, 0, 6, (tfx_gl_sizei)count);
}

//Each layer in the group as its CPU instances then its GPU instances
static void tfx_impl_gl__render_billboards(tfx_impl_gl_t *renderer, int writes_depth) {
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_VERTEX_READ_BINDING + 0, &renderer->instance_buffers[renderer->current_instance_buffer]);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_VERTEX_READ_BINDING + 1, &renderer->instance_buffers[renderer->current_instance_buffer ^ 1]);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_VERTEX_READ_BINDING + 2, &renderer->image_data);
	tfx_impl_gl__bind_storage(TFX_IMPL_GL_VERTEX_READ_BINDING + 3, &renderer->emitter_properties);
	tfx_gl.ActiveTexture(TFX_GL_TEXTURE0 + TFX_IMPL_GL_FRAGMENT_TEXTURE_BINDING + 1);
	tfx_gl.BindTexture(TFX_GL_TEXTURE_2D_ARRAY, renderer->color_ramps);
	tfx_gl_uint bound_program = 0;
	for (tfxU32 layer = 0; layer != renderer->library_layer_count; ++layer) {
		const tfxLayerFlags flags = renderer->layer_flags[layer];
		if (tfx_impl_LayerWritesDepth(flags) != writes_depth || renderer->layer_count[layer] + renderer->gpu_layer_count[layer] == 0) {
			continue;
		}
		tfx_gl_uint program = renderer->billboard_programs[writes_depth][tfx_impl_LayerIsLit(flags)];
		if (program != bound_program) {
			tfx_gl.UseProgram(program);
			bound_program = program;
		}
		tfx_impl_gl__draw_instance_range(renderer, renderer->layer_start[layer], renderer->layer_count[layer]);
		tfx_impl_gl__draw_instance_range(renderer, renderer->gpu_layer_start[layer], renderer->gpu_layer_count[layer]);
	}
}

void tfx_impl_gl_Render(tfx_impl_gl_t *renderer, const tfx_impl_uniforms_t *uniforms) {
	if (!renderer->shape_atlas || (renderer->instance_count == 0 && renderer->playback_sprite_count == 0
		&& !tfx_impl_gl__ribbon_set_active(&renderer->live_ribbons) && !tfx_impl_gl__ribbon_set_active(&renderer->baked_ribbons))) {
		return;
	}
	tfx_impl_gl__write_uniforms(&renderer->frame_constants, TFX_IMPL_GL_VERTEX_UNIFORM_BINDING, uniforms, sizeof(tfx_impl_uniforms_t));
	//The lit fragments read the light from the same block
	tfx_gl.BindBufferBase(TFX_GL_UNIFORM_BUFFER, TFX_IMPL_GL_FRAGMENT_UNIFORM_BINDING, renderer->frame_constants.buffer);

	//Quad corners come from gl_VertexID and the instances and ribbon vertices are read from buffers, so the vertex array is empty
	tfx_gl.BindVertexArray(renderer->empty_vertex_array);
	tfx_gl.BlendFunc(TFX_GL_ONE, TFX_GL_ONE_MINUS_SRC_ALPHA);
	tfx_gl.Enable(TFX_GL_DEPTH_TEST);
	tfx_gl.DepthFunc(TFX_GL_LEQUAL);
	tfx_gl.Disable(TFX_GL_CULL_FACE);
	tfx_gl.ActiveTexture(TFX_GL_TEXTURE0 + TFX_IMPL_GL_FRAGMENT_TEXTURE_BINDING);
	tfx_gl.BindTexture(TFX_GL_TEXTURE_2D_ARRAY, renderer->shape_atlas);

	//The layers that write depth draw first with blending off, so they hide the blended layers behind them
	for (int group = 0; group != 2; ++group) {
		const int writes_depth = group == 0;
		if (writes_depth) {
			tfx_gl.Disable(TFX_GL_BLEND);
			tfx_gl.DepthMask(1);
		} else {
			tfx_gl.Enable(TFX_GL_BLEND);
			tfx_gl.DepthMask(0);
		}
		if (renderer->instance_count > 0) {
			tfx_impl_gl__render_billboards(renderer, writes_depth);
		}
		if (tfx_impl_gl__ribbon_set_active(&renderer->live_ribbons)) {
			tfx_impl_gl__render_ribbons(renderer, &renderer->live_ribbons, &renderer->image_data, renderer->color_ramps, writes_depth);
		}
		//Baked billboards come as one range across their layers, so they all draw blended and unlit
		if (!writes_depth && renderer->playback_sprite_count > 0) {
			//Baked sprites come with their own shape frames, emitter properties and color ramps, the atlas is shared
			tfx_gl.UseProgram(renderer->billboard_baked_program);
			tfx_impl_gl__bind_storage(TFX_IMPL_GL_VERTEX_READ_BINDING + 0, &renderer->playback_instances);
			tfx_impl_gl__bind_storage(TFX_IMPL_GL_VERTEX_READ_BINDING + 1, &renderer->baked_image_data);
			tfx_impl_gl__bind_storage(TFX_IMPL_GL_VERTEX_READ_BINDING + 2, &renderer->baked_emitter_properties);
			tfx_gl.ActiveTexture(TFX_GL_TEXTURE0 + TFX_IMPL_GL_FRAGMENT_TEXTURE_BINDING + 1);
			tfx_gl.BindTexture(TFX_GL_TEXTURE_2D_ARRAY, renderer->baked_color_ramps);
			tfx_impl_gl__draw_instance_range(renderer, 0, renderer->playback_sprite_count);
		}
		if (tfx_impl_gl__ribbon_set_active(&renderer->baked_ribbons)) {
			tfx_impl_gl__render_ribbons(renderer, &renderer->baked_ribbons, &renderer->baked_image_data, renderer->baked_color_ramps, writes_depth);
		}
	}
	//glClear doesn't clear depth while depth writes are off
	tfx_gl.DepthMask(1);
	tfx_gl.UseProgram(0);
	tfx_gl.BindVertexArray(0);
}

static void tfx_impl_gl__release_ribbon_set(tfx_impl_gl_ribbon_set_t *set) {
	tfx_impl_gl__release_buffer(&set->segments);
	tfx_impl_gl__release_buffer(&set->instances);
	tfx_impl_gl__release_buffer(&set->emitters);
	tfx_impl_gl__release_buffer(&set->vertices);
	tfx_impl_gl__release_buffer(&set->indices);
	tfx_impl_gl__release_buffer(&set->graph_lookups);
}

static void tfx_impl_gl__delete_program(tfxU32 *program) {
	if (*program) {
		tfx_gl.DeleteProgram(*program);
		*program = 0;
	}
}

static void tfx_impl_gl__delete_texture(tfxU32 *texture) {
	if (*texture) {
		tfx_gl.DeleteTextures(1, texture);
		*texture = 0;
	}
}

void tfx_impl_gl_Shutdown(tfx_impl_gl_t *renderer) {
	if (!tfx_gl.DeleteBuffers) {
		return;
	}
	for (int index = 0; index != 2; ++index) {
		tfx_impl_gl__release_buffer(&renderer->instance_buffers[index]);
	}
	tfx_impl_gl__release_buffer(&renderer->frame_constants);
	tfx_impl_gl__release_buffer(&renderer->draw_constants);
	tfx_impl_gl__release_buffer(&renderer->compute_constants);
	tfx_impl_gl__release_buffer(&renderer->image_data);
	tfx_impl_gl__release_buffer(&renderer->emitter_properties);
	tfx_impl_gl__release_buffer(&renderer->graph_lookups);
	tfx_impl_gl__release_ribbon_set(&renderer->live_ribbons);
	tfx_impl_gl__release_ribbon_set(&renderer->baked_ribbons);
	tfx_impl_gl__release_buffer(&renderer->particle_pool);
	tfx_impl_gl__release_buffer(&renderer->particle_upload);
	tfx_impl_gl__release_buffer(&renderer->baked_image_data);
	tfx_impl_gl__release_buffer(&renderer->baked_emitter_properties);
	tfx_impl_gl__release_buffer(&renderer->baked_sprites);
	tfx_impl_gl__release_buffer(&renderer->animation_offsets);
	tfx_impl_gl__release_buffer(&renderer->animation_instances);
	tfx_impl_gl__release_buffer(&renderer->playback_instances);
	tfx_impl_gl__release_buffer(&renderer->baked_ribbon_data);
	tfx_impl_gl__release_buffer(&renderer->baked_ribbon_properties);
	tfx_impl_gl__release_buffer(&renderer->animation_ribbon_offsets);
	tfx_impl_gl__delete_texture(&renderer->shape_atlas);
	tfx_impl_gl__delete_texture(&renderer->color_ramps);
	tfx_impl_gl__delete_texture(&renderer->baked_color_ramps);
	free(renderer->color_ramp_layer_versions);
	renderer->color_ramp_layer_versions = NULL;
	renderer->color_ramp_layer_capacity = 0;
	for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
		for (int lit = 0; lit != 2; ++lit) {
			tfx_impl_gl__delete_program(&renderer->billboard_programs[writes_depth][lit]);
			tfx_impl_gl__delete_program(&renderer->ribbon_programs[writes_depth][lit]);
		}
	}
	tfx_impl_gl__delete_program(&renderer->billboard_baked_program);
	tfx_impl_gl__delete_program(&renderer->particle_kernel);
	tfx_impl_gl__delete_program(&renderer->ribbon_compute_program);
	tfx_impl_gl__delete_program(&renderer->sprite_playback_program);
	tfx_impl_gl__delete_program(&renderer->ribbon_playback_program);
	if (renderer->empty_vertex_array) {
		tfx_gl.DeleteVertexArrays(1, &renderer->empty_vertex_array);
		renderer->empty_vertex_array = 0;
	}
}
