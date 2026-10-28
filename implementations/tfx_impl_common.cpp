#include "tfx_impl_common.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//Define TFX_IMPL_NO_STB_IMPLEMENTATION if your project already compiles stb_image and stb_rect_pack
#ifndef TFX_IMPL_NO_STB_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_PNG
#define STB_RECT_PACK_IMPLEMENTATION
#endif
#include "third_party/stb_image.h"
#include "third_party/stb_rect_pack.h"
#define BCDEC_STATIC
#define BCDEC_IMPLEMENTATION
#include "third_party/bcdec.h"

//Zstd supercompressed KTX2 shapes need third_party/zstddeclib.c compiled in, the CMake project does that and defines this
#ifdef TFX_IMPL_KTX2_ZSTD
extern "C" size_t ZSTD_decompress(void *destination, size_t destination_capacity, const void *source, size_t source_size);
extern "C" unsigned ZSTD_isError(size_t code);
#endif

#define TFX_IMPL_ATLAS_PADDING 2
#define TFX_IMPL_ATLAS_MIN_LAYER_SIZE 2048
#define TFX_IMPL_ATLAS_MAX_LAYER_SIZE 16384
#define TFX_IMPL_KTX2_HEADER_SIZE 80
#define TFX_IMPL_KTX2_LEVEL_ENTRY_SIZE 24
#define TFX_IMPL_KTX2_SUPERCOMPRESSION_NONE 0
#define TFX_IMPL_KTX2_SUPERCOMPRESSION_ZSTD 2
#define TFX_IMPL_VK_FORMAT_BC4_UNORM 139
#define TFX_IMPL_VK_FORMAT_BC5_UNORM 141
#define TFX_IMPL_VK_FORMAT_BC7_UNORM 145

//--- Impl_math

tfx_impl_vec3_t tfx_impl_Vec3(float x, float y, float z) {
	tfx_impl_vec3_t result = { x, y, z };
	return result;
}

tfx_impl_vec3_t tfx_impl_AddVec3(tfx_impl_vec3_t a, tfx_impl_vec3_t b) {
	return tfx_impl_Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

tfx_impl_vec3_t tfx_impl_SubVec3(tfx_impl_vec3_t a, tfx_impl_vec3_t b) {
	return tfx_impl_Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

tfx_impl_vec3_t tfx_impl_ScaleVec3(tfx_impl_vec3_t v, float scale) {
	return tfx_impl_Vec3(v.x * scale, v.y * scale, v.z * scale);
}

tfx_impl_vec3_t tfx_impl_CrossVec3(tfx_impl_vec3_t a, tfx_impl_vec3_t b) {
	return tfx_impl_Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

float tfx_impl_DotVec3(tfx_impl_vec3_t a, tfx_impl_vec3_t b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

tfx_impl_vec3_t tfx_impl_NormalizeVec3(tfx_impl_vec3_t v) {
	float length = sqrtf(tfx_impl_DotVec3(v, v));
	if (length == 0.f) {
		return v;
	}
	return tfx_impl_ScaleVec3(v, 1.f / length);
}

tfx_impl_mat4_t tfx_impl_IdentityMat4(void) {
	tfx_impl_mat4_t result;
	memset(&result, 0, sizeof(result));
	result.columns[0][0] = 1.f;
	result.columns[1][1] = 1.f;
	result.columns[2][2] = 1.f;
	result.columns[3][3] = 1.f;
	return result;
}

tfx_impl_mat4_t tfx_impl_MultiplyMat4(const tfx_impl_mat4_t *a, const tfx_impl_mat4_t *b) {
	tfx_impl_mat4_t result;
	for (int column = 0; column != 4; ++column) {
		for (int row = 0; row != 4; ++row) {
			float sum = 0.f;
			for (int index = 0; index != 4; ++index) {
				sum += a->columns[index][row] * b->columns[column][index];
			}
			result.columns[column][row] = sum;
		}
	}
	return result;
}

tfx_impl_mat4_t tfx_impl_LookAt(tfx_impl_vec3_t eye, tfx_impl_vec3_t target, tfx_impl_vec3_t up) {
	tfx_impl_vec3_t forward = tfx_impl_NormalizeVec3(tfx_impl_SubVec3(target, eye));
	tfx_impl_vec3_t side = tfx_impl_NormalizeVec3(tfx_impl_CrossVec3(forward, up));
	tfx_impl_vec3_t camera_up = tfx_impl_CrossVec3(side, forward);
	tfx_impl_mat4_t result = tfx_impl_IdentityMat4();
	result.columns[0][0] = side.x;
	result.columns[1][0] = side.y;
	result.columns[2][0] = side.z;
	result.columns[0][1] = camera_up.x;
	result.columns[1][1] = camera_up.y;
	result.columns[2][1] = camera_up.z;
	result.columns[0][2] = -forward.x;
	result.columns[1][2] = -forward.y;
	result.columns[2][2] = -forward.z;
	result.columns[3][0] = -tfx_impl_DotVec3(side, eye);
	result.columns[3][1] = -tfx_impl_DotVec3(camera_up, eye);
	result.columns[3][2] = tfx_impl_DotVec3(forward, eye);
	return result;
}

tfx_impl_mat4_t tfx_impl_Perspective(float field_of_view, float aspect_ratio, float near_plane, float far_plane, tfxU32 clip_space_flags) {
	float focal_length = 1.f / tanf(field_of_view * 0.5f);
	tfx_impl_mat4_t result;
	memset(&result, 0, sizeof(result));
	result.columns[0][0] = focal_length / aspect_ratio;
	result.columns[1][1] = focal_length;
	result.columns[2][3] = -1.f;
	if (clip_space_flags & tfx_impl_clip_space_depth_zero_to_one) {
		result.columns[2][2] = far_plane / (near_plane - far_plane);
		result.columns[3][2] = -(far_plane * near_plane) / (far_plane - near_plane);
	} else {
		result.columns[2][2] = -(far_plane + near_plane) / (far_plane - near_plane);
		result.columns[3][2] = -(2.f * far_plane * near_plane) / (far_plane - near_plane);
	}
	if (clip_space_flags & tfx_impl_clip_space_flip_y) {
		result.columns[1][1] *= -1.f;
	}
	return result;
}

//--- Impl_timer

tfx_impl_timer_t tfx_impl_CreateTimer(double update_frequency) {
	tfx_impl_timer_t timer;
	memset(&timer, 0, sizeof(timer));
	timer.update_frequency = update_frequency;
	timer.update_length = 1000.0 / update_frequency;
	timer.max_ticks_per_frame = 8;
	return timer;
}

tfxU32 tfx_impl_AdvanceTimer(tfx_impl_timer_t *timer, double elapsed_milliseconds) {
	timer->accumulator += elapsed_milliseconds;
	tfxU32 ticks = (tfxU32)(timer->accumulator / timer->update_length);
	if (ticks > timer->max_ticks_per_frame) {
		ticks = timer->max_ticks_per_frame;
		timer->accumulator = ticks * timer->update_length;
	}
	timer->accumulator -= ticks * timer->update_length;
	timer->lerp = (float)(timer->accumulator / timer->update_length);
	timer->ticks_this_frame = ticks;
	return ticks;
}

//--- Impl_camera

tfx_impl_camera_t tfx_impl_CreateCamera(tfx_impl_vec3_t position, tfx_impl_vec3_t target, float field_of_view) {
	tfx_impl_camera_t camera;
	memset(&camera, 0, sizeof(camera));
	camera.position = position;
	camera.target = target;
	camera.up = tfx_impl_Vec3(0.f, 1.f, 0.f);
	camera.field_of_view = field_of_view;
	camera.near_plane = 0.1f;
	camera.far_plane = 1000.f;
	camera.view = tfx_impl_IdentityMat4();
	camera.projection = tfx_impl_IdentityMat4();
	return camera;
}

void tfx_impl_OrbitCamera(tfx_impl_camera_t *camera, float yaw, float pitch, float distance) {
	float horizontal = cosf(pitch) * distance;
	tfx_impl_vec3_t offset = tfx_impl_Vec3(sinf(yaw) * horizontal, sinf(pitch) * distance, cosf(yaw) * horizontal);
	camera->position = tfx_impl_AddVec3(camera->target, offset);
}

void tfx_impl_UpdateCameraMatrices(tfx_impl_camera_t *camera, float aspect_ratio, tfxU32 clip_space_flags) {
	camera->view = tfx_impl_LookAt(camera->position, camera->target, camera->up);
	camera->projection = tfx_impl_Perspective(camera->field_of_view, aspect_ratio, camera->near_plane, camera->far_plane, clip_space_flags);
}

tfx_impl_vec3_t tfx_impl_CameraFront(const tfx_impl_camera_t *camera) {
	return tfx_impl_NormalizeVec3(tfx_impl_SubVec3(camera->target, camera->position));
}

//--- Impl_shapes

static const unsigned char tfx_impl__ktx2_identifier[12] = { 0xAB, 0x4B, 0x54, 0x58, 0x20, 0x32, 0x30, 0xBB, 0x0D, 0x0A, 0x1A, 0x0A };

static bool tfx_impl__is_alpha_only(tfx_image_format format) {
	return format == tfx_image_format_a8_png || format == tfx_image_format_a_bc4_ktx2;
}

//Expands 1, 2 or 4 channels to RGBA8: alpha only is (1, 1, 1, a), luminance (l, l, l, 1), luminance alpha (l, l, l, a)
static void tfx_impl__expand_to_rgba(const unsigned char *source, tfxU32 channels, tfx_image_format format, size_t pixel_count, unsigned char *destination) {
	bool alpha_only = tfx_impl__is_alpha_only(format);
	for (size_t pixel = 0; pixel != pixel_count; ++pixel) {
		const unsigned char *in = source + pixel * channels;
		unsigned char *out = destination + pixel * 4;
		if (channels == 4) {
			out[0] = in[0]; out[1] = in[1]; out[2] = in[2]; out[3] = in[3];
		} else if (channels == 2) {
			out[0] = in[0]; out[1] = in[0]; out[2] = in[0]; out[3] = in[1];
		} else if (alpha_only) {
			out[0] = 255; out[1] = 255; out[2] = 255; out[3] = in[0];
		} else {
			out[0] = in[0]; out[1] = in[0]; out[2] = in[0]; out[3] = 255;
		}
	}
}

static tfxU32 tfx_impl__read_u32(const unsigned char *bytes) {
	return (tfxU32)bytes[0] | ((tfxU32)bytes[1] << 8) | ((tfxU32)bytes[2] << 16) | ((tfxU32)bytes[3] << 24);
}

static tfxU64 tfx_impl__read_u64(const unsigned char *bytes) {
	return (tfxU64)tfx_impl__read_u32(bytes) | ((tfxU64)tfx_impl__read_u32(bytes + 4) << 32);
}

//Returns why the shape couldn't be read, or NULL
static const char *tfx_impl__load_png_shape(tfx_impl_shape_t *shape, tfx_image_format format, tfx_image_data_t *image_data, const unsigned char *data, int size) {
	tfxU32 channels = 4;
	if (format == tfx_image_format_la8_png) {
		channels = 2;
	} else if (format == tfx_image_format_l8_png || format == tfx_image_format_a8_png) {
		channels = 1;
	}
	int width = 0;
	int height = 0;
	int file_channels = 0;
	unsigned char *decoded = format == tfx_image_format_rgba8_raw ? NULL : stbi_load_from_memory(data, size, &width, &height, &file_channels, (int)channels);
	if (decoded) {
		shape->pixels = (unsigned char *)malloc((size_t)width * height * 4);
		if (!shape->pixels) {
			stbi_image_free(decoded);
			return "there wasn't enough memory to load it";
		}
		tfx_impl__expand_to_rgba(decoded, channels, format, (size_t)width * height, shape->pixels);
		stbi_image_free(decoded);
	} else {
		//Libraries saved before the format was recorded can hold raw RGBA sheets
		if (format != tfx_image_format_rgba8_raw && format != tfx_image_format_unknown) {
			return "its png couldn't be decoded";
		}
		width = tfx_GetImageWidth(image_data);
		height = tfx_GetImageHeight(image_data);
		size_t sheet_size = (size_t)width * height * 4;
		shape->pixels = (unsigned char *)calloc(1, sheet_size);
		if (!shape->pixels) {
			return "there wasn't enough memory to load it";
		}
		memcpy(shape->pixels, data, (size_t)size < sheet_size ? (size_t)size : sheet_size);
	}
	shape->sheet_width = (tfxU32)width;
	shape->sheet_height = (tfxU32)height;
	int frames = tfx_GetImageFrameCount(image_data);
	shape->frame_count = frames > 1 ? (tfxU32)frames : 1;
	shape->frame_width = shape->frame_count > 1 ? (tfxU32)tfx_GetImageWidth(image_data) : (tfxU32)width;
	shape->frame_height = shape->frame_count > 1 ? (tfxU32)tfx_GetImageHeight(image_data) : (tfxU32)height;
	return NULL;
}

//Every layer of the KTX2 is a frame. Only the largest mip level is decoded, the atlas makes its own mips
static const char *tfx_impl__load_ktx2_shape(tfx_impl_shape_t *shape, tfx_image_format format, const unsigned char *data, size_t size) {
	tfxU32 vk_format = TFX_IMPL_VK_FORMAT_BC7_UNORM;
	tfxU32 block_bytes = 16;
	tfxU32 channels = 4;
	if (format == tfx_image_format_a_bc4_ktx2 || format == tfx_image_format_l_bc4_ktx2) {
		vk_format = TFX_IMPL_VK_FORMAT_BC4_UNORM;
		block_bytes = 8;
		channels = 1;
	} else if (format == tfx_image_format_la_bc5_ktx2) {
		vk_format = TFX_IMPL_VK_FORMAT_BC5_UNORM;
		channels = 2;
	}
	if (size < TFX_IMPL_KTX2_HEADER_SIZE + TFX_IMPL_KTX2_LEVEL_ENTRY_SIZE || memcmp(data, tfx_impl__ktx2_identifier, sizeof(tfx_impl__ktx2_identifier)) != 0) {
		return "it isn't a KTX2 file";
	}
	if (tfx_impl__read_u32(data + 12) != vk_format) {
		return "the KTX2 vkFormat doesn't match its shape format";
	}
	if (tfx_impl__read_u32(data + 28) != 0 || tfx_impl__read_u32(data + 36) != 1) {
		return "it isn't a 2D KTX2 file with a single face";
	}
	tfxU32 width = tfx_impl__read_u32(data + 20);
	tfxU32 height = tfx_impl__read_u32(data + 24);
	tfxU32 layers = tfx_impl__read_u32(data + 32);
	tfxU32 supercompression = tfx_impl__read_u32(data + 44);
	layers = layers ? layers : 1;
	if (!width || !height) {
		return "the KTX2 file has no size";
	}
	const unsigned char *level_entry = data + TFX_IMPL_KTX2_HEADER_SIZE;
	tfxU64 offset = tfx_impl__read_u64(level_entry);
	tfxU64 stored_size = tfx_impl__read_u64(level_entry + 8);
	tfxU64 level_size = tfx_impl__read_u64(level_entry + 16);
	tfxU32 blocks_wide = (width + 3) / 4;
	tfxU32 blocks_high = (height + 3) / 4;
	size_t layer_size = (size_t)blocks_wide * blocks_high * block_bytes;
	if (level_size != layer_size * layers || !stored_size || offset > size || stored_size > size - offset) {
		return "its largest mip level is the wrong size or lies outside the file";
	}
	unsigned char *blocks = NULL;
	if (supercompression == TFX_IMPL_KTX2_SUPERCOMPRESSION_NONE) {
		if (stored_size != level_size) {
			return "its largest mip level is the wrong size";
		}
		blocks = (unsigned char *)(data + offset);
	} else if (supercompression == TFX_IMPL_KTX2_SUPERCOMPRESSION_ZSTD) {
#ifdef TFX_IMPL_KTX2_ZSTD
		blocks = (unsigned char *)malloc((size_t)level_size);
		if (!blocks) {
			return "there wasn't enough memory to inflate it";
		}
		size_t inflated_size = ZSTD_decompress(blocks, (size_t)level_size, data + offset, (size_t)stored_size);
		if (ZSTD_isError(inflated_size) || inflated_size != level_size) {
			free(blocks);
			return "its zstd supercompression couldn't be inflated";
		}
#else
		return "it's zstd supercompressed and TFX_IMPL_KTX2_ZSTD isn't defined";
#endif
	} else {
		return "it uses a KTX2 supercompression scheme this loader can't inflate";
	}

	shape->frame_count = layers;
	shape->frame_width = width;
	shape->frame_height = height;
	shape->sheet_width = width;
	shape->sheet_height = height * layers;
	shape->pixels = (unsigned char *)malloc((size_t)width * height * layers * 4);
	if (shape->pixels) {
		unsigned char texels[16 * 4];
		for (tfxU32 layer = 0; layer != layers; ++layer) {
			const unsigned char *layer_blocks = blocks + layer * layer_size;
			for (tfxU32 block_y = 0; block_y != blocks_high; ++block_y) {
				for (tfxU32 block_x = 0; block_x != blocks_wide; ++block_x) {
					const unsigned char *block = layer_blocks + ((size_t)block_y * blocks_wide + block_x) * block_bytes;
					if (vk_format == TFX_IMPL_VK_FORMAT_BC4_UNORM) {
						bcdec_bc4(block, texels, 4);
					} else if (vk_format == TFX_IMPL_VK_FORMAT_BC5_UNORM) {
						bcdec_bc5(block, texels, 4 * 2);
					} else {
						bcdec_bc7(block, texels, 4 * 4);
					}
					for (tfxU32 texel_y = 0; texel_y != 4; ++texel_y) {
						tfxU32 y = block_y * 4 + texel_y;
						for (tfxU32 texel_x = 0; texel_x != 4 && y < height; ++texel_x) {
							tfxU32 x = block_x * 4 + texel_x;
							if (x < width) {
								unsigned char *destination = shape->pixels + (((size_t)layer * height + y) * width + x) * 4;
								tfx_impl__expand_to_rgba(texels + (texel_y * 4 + texel_x) * channels, channels, format, 1, destination);
							}
						}
					}
				}
			}
		}
	}
	if (blocks != data + offset) {
		free(blocks);
	}
	return shape->pixels ? NULL : "there wasn't enough memory to decode it";
}

static tfx_impl_shape_t *tfx_impl__add_shape(tfx_impl_shapes_t *shapes) {
	if (shapes->count == shapes->capacity) {
		tfxU32 capacity = shapes->capacity ? shapes->capacity * 2 : 64;
		tfx_impl_shape_t **grown = (tfx_impl_shape_t **)realloc(shapes->shapes, sizeof(tfx_impl_shape_t *) * capacity);
		if (!grown) {
			return NULL;
		}
		shapes->shapes = grown;
		shapes->capacity = capacity;
	}
	tfx_impl_shape_t *shape = (tfx_impl_shape_t *)calloc(1, sizeof(tfx_impl_shape_t));
	if (shape) {
		shapes->shapes[shapes->count++] = shape;
	}
	return shape;
}

void tfx_impl_ShapeLoader(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_memory_size, void *user_data) {
	tfx_impl_shapes_t *shapes = (tfx_impl_shapes_t *)user_data;
	tfx_impl_shape_t *shape = tfx_impl__add_shape(shapes);
	if (!shape) {
		return;
	}
	tfx_image_format format = tfx_GetImageFormat(image_data);
	const char *error = NULL;
	if (format >= tfx_image_format_a_bc4_ktx2 && format <= tfx_image_format_rgba_bc7_ktx2) {
		error = tfx_impl__load_ktx2_shape(shape, format, (const unsigned char *)raw_image_data, (size_t)image_memory_size);
	} else if (format <= tfx_image_format_rgba8_raw) {
		error = tfx_impl__load_png_shape(shape, format, image_data, (const unsigned char *)raw_image_data, image_memory_size);
	} else {
		error = "its image format isn't one this loader recognises";
	}
	if (error) {
		free(shape->pixels);
		shape->pixels = (unsigned char *)malloc(4);
		if (shape->pixels) {
			memset(shape->pixels, 255, 4);
		}
		shape->unreadable = true;
		shape->frame_count = 1;
		shape->frame_width = shape->frame_height = 1;
		shape->sheet_width = shape->sheet_height = 1;
		shapes->unreadable_count++;
		printf("Skipped particle shape %s (image format %i) because %s. It will draw as a white square.\n", filename ? filename : "", (int)format, error);
	}
	tfx_SetImagePointer(image_data, shape);
}

void tfx_impl_FreeShapes(tfx_impl_shapes_t *shapes) {
	for (tfxU32 index = 0; index != shapes->count; ++index) {
		free(shapes->shapes[index]->pixels);
		free(shapes->shapes[index]->regions);
		free(shapes->shapes[index]);
	}
	free(shapes->shapes);
	memset(shapes, 0, sizeof(tfx_impl_shapes_t));
}

//--- Impl_atlas

typedef struct tfx_impl_atlas_frame_s {
	tfx_impl_shape_t *shape;
	tfxU32 frame;
	tfxU32 layer;
	tfxU32 x;
	tfxU32 y;
} tfx_impl_atlas_frame_t;

//Repeats the frame's edge pixels into the padding so bilinear filtering doesn't pick up its neighbours
static void tfx_impl__copy_frame_to_atlas(const tfx_impl_atlas_frame_t *atlas_frame, tfx_impl_atlas_t *atlas) {
	const tfx_impl_shape_t *shape = atlas_frame->shape;
	tfxU32 columns = shape->frame_width ? shape->sheet_width / shape->frame_width : 1;
	columns = columns ? columns : 1;
	int sheet_x = (int)((atlas_frame->frame % columns) * shape->frame_width);
	int sheet_y = (int)((atlas_frame->frame / columns) * shape->frame_height);
	int frame_width = (int)shape->frame_width;
	int frame_height = (int)shape->frame_height;
	unsigned char *layer_pixels = atlas->pixels + (size_t)atlas_frame->layer * atlas->layer_size * atlas->layer_size * 4;
	for (int padded_y = -TFX_IMPL_ATLAS_PADDING; padded_y != frame_height + TFX_IMPL_ATLAS_PADDING; ++padded_y) {
		int source_y = sheet_y + (padded_y < 0 ? 0 : (padded_y >= frame_height ? frame_height - 1 : padded_y));
		size_t destination_row = (size_t)(atlas_frame->y + TFX_IMPL_ATLAS_PADDING + padded_y) * atlas->layer_size;
		for (int padded_x = -TFX_IMPL_ATLAS_PADDING; padded_x != frame_width + TFX_IMPL_ATLAS_PADDING; ++padded_x) {
			int source_x = sheet_x + (padded_x < 0 ? 0 : (padded_x >= frame_width ? frame_width - 1 : padded_x));
			unsigned char *destination = layer_pixels + (destination_row + atlas_frame->x + TFX_IMPL_ATLAS_PADDING + padded_x) * 4;
			//A sheet whose grid doesn't hold every frame it claims leaves the missing frames blank
			if (shape->pixels && source_x < (int)shape->sheet_width && source_y < (int)shape->sheet_height) {
				memcpy(destination, shape->pixels + ((size_t)source_y * shape->sheet_width + source_x) * 4, 4);
			} else {
				memset(destination, 0, 4);
			}
		}
	}
}

bool tfx_impl_BuildShapeAtlas(tfx_impl_shapes_t *shapes, tfx_impl_atlas_t *atlas) {
	memset(atlas, 0, sizeof(tfx_impl_atlas_t));
	tfxU32 frame_total = 0;
	tfxU32 largest_side = 0;
	for (tfxU32 index = 0; index != shapes->count; ++index) {
		tfx_impl_shape_t *shape = shapes->shapes[index];
		frame_total += shape->frame_count;
		largest_side = shape->frame_width > largest_side ? shape->frame_width : largest_side;
		largest_side = shape->frame_height > largest_side ? shape->frame_height : largest_side;
	}
	if (frame_total == 0) {
		return false;
	}
	atlas->layer_size = TFX_IMPL_ATLAS_MIN_LAYER_SIZE;
	while (atlas->layer_size < largest_side + TFX_IMPL_ATLAS_PADDING * 2 && atlas->layer_size < TFX_IMPL_ATLAS_MAX_LAYER_SIZE) {
		atlas->layer_size *= 2;
	}
	if (largest_side + TFX_IMPL_ATLAS_PADDING * 2 > atlas->layer_size) {
		printf("A shape frame of %u pixels is too big for a %u pixel atlas layer\n", largest_side, atlas->layer_size);
		return false;
	}

	tfx_impl_atlas_frame_t *frames = (tfx_impl_atlas_frame_t *)malloc(sizeof(tfx_impl_atlas_frame_t) * frame_total);
	stbrp_rect *rects = (stbrp_rect *)malloc(sizeof(stbrp_rect) * frame_total);
	stbrp_node *nodes = (stbrp_node *)malloc(sizeof(stbrp_node) * atlas->layer_size);
	if (!frames || !rects || !nodes) {
		free(frames);
		free(rects);
		free(nodes);
		return false;
	}
	tfxU32 frame_index = 0;
	for (tfxU32 index = 0; index != shapes->count; ++index) {
		tfx_impl_shape_t *shape = shapes->shapes[index];
		free(shape->regions);
		shape->regions = (tfx_impl_atlas_region_t *)calloc(shape->frame_count, sizeof(tfx_impl_atlas_region_t));
		for (tfxU32 frame = 0; frame != shape->frame_count; ++frame) {
			frames[frame_index].shape = shape;
			frames[frame_index].frame = frame;
			rects[frame_index].id = (int)frame_index;
			rects[frame_index].w = (stbrp_coord)(shape->frame_width + TFX_IMPL_ATLAS_PADDING * 2);
			rects[frame_index].h = (stbrp_coord)(shape->frame_height + TFX_IMPL_ATLAS_PADDING * 2);
			frame_index++;
		}
	}

	//Each pass fills one layer with whatever fits, and what's left moves to the front for the next layer
	tfxU32 remaining = frame_total;
	while (remaining > 0) {
		stbrp_context context;
		stbrp_init_target(&context, (int)atlas->layer_size, (int)atlas->layer_size, nodes, (int)atlas->layer_size);
		stbrp_pack_rects(&context, rects, (int)remaining);
		tfxU32 unpacked = 0;
		for (tfxU32 index = 0; index != remaining; ++index) {
			if (rects[index].was_packed) {
				tfx_impl_atlas_frame_t *atlas_frame = &frames[rects[index].id];
				atlas_frame->layer = atlas->layer_count;
				atlas_frame->x = (tfxU32)rects[index].x;
				atlas_frame->y = (tfxU32)rects[index].y;
			} else {
				rects[unpacked++] = rects[index];
			}
		}
		atlas->layer_count++;
		remaining = unpacked;
	}
	free(rects);
	free(nodes);

	size_t layer_bytes = (size_t)atlas->layer_size * atlas->layer_size * 4;
	atlas->pixels = (unsigned char *)calloc(atlas->layer_count, layer_bytes);
	if (!atlas->pixels) {
		free(frames);
		return false;
	}
	float inverse_size = 1.f / (float)atlas->layer_size;
	for (frame_index = 0; frame_index != frame_total; ++frame_index) {
		tfx_impl_atlas_frame_t *atlas_frame = &frames[frame_index];
		tfx_impl__copy_frame_to_atlas(atlas_frame, atlas);
		tfx_impl_atlas_region_t *region = &atlas_frame->shape->regions[atlas_frame->frame];
		float left = (float)(atlas_frame->x + TFX_IMPL_ATLAS_PADDING);
		float top = (float)(atlas_frame->y + TFX_IMPL_ATLAS_PADDING);
		region->uv[0] = left * inverse_size;
		region->uv[1] = top * inverse_size;
		region->uv[2] = (left + atlas_frame->shape->frame_width) * inverse_size;
		region->uv[3] = (top + atlas_frame->shape->frame_height) * inverse_size;
		region->layer = atlas_frame->layer;
	}
	free(frames);
	for (tfxU32 index = 0; index != shapes->count; ++index) {
		free(shapes->shapes[index]->pixels);
		shapes->shapes[index]->pixels = NULL;
	}
	return true;
}

void tfx_impl_FreeAtlasPixels(tfx_impl_atlas_t *atlas) {
	free(atlas->pixels);
	atlas->pixels = NULL;
}

void tfx_impl_GetUV(void *pointer, tfx_gpu_image_data_t *image_data, int offset) {
	tfx_impl_shape_t *shape = (tfx_impl_shape_t *)pointer;
	//The library also calls this while it loads, before the atlas exists
	if (!shape->regions) {
		return;
	}
	tfxU32 frame = (tfxU32)offset < shape->frame_count ? (tfxU32)offset : shape->frame_count - 1;
	const tfx_impl_atlas_region_t *region = &shape->regions[frame];
	image_data->uv.x = region->uv[0];
	image_data->uv.y = region->uv[1];
	image_data->uv.z = region->uv[2];
	image_data->uv.w = region->uv[3];
	image_data->texture_array_index = region->layer;
}

//--- Impl_frame

static bool tfx_impl__ensure_staging(void **memory, size_t *capacity, size_t required) {
	if (*capacity >= required && *memory) {
		return true;
	}
	free(*memory);
	*memory = malloc(required ? required : 16);
	*capacity = *memory ? required : 0;
	return *memory != NULL;
}

//Uses the ribbon batch API even for one stage, which is the general flow for drawing any number of stages with one upload
static void tfx_impl__gather_ribbons(tfx_stage stage, tfx_impl_ribbon_staging_t *staging, tfx_impl_ribbon_frame_t *ribbon_frame) {
	if (!tfx_impl__ensure_staging(&staging->segments, &staging->segments_capacity, tfx_GetTotalSegmentBufferMaxSizeInBytes())
		|| !tfx_impl__ensure_staging(&staging->ribbons, &staging->ribbons_capacity, tfx_GetTotalRibbonBufferMaxSizeInBytes())
		|| !tfx_impl__ensure_staging(&staging->emitters, &staging->emitters_capacity, tfx_GetTotalEmitterBufferMaxSizeInBytes())) {
		return;
	}
	tfx_ribbon_batch_t batch = tfx_BeginRibbonBatch(staging->segments, staging->segments_capacity, staging->ribbons, staging->ribbons_capacity,
		staging->emitters, staging->emitters_capacity);
	if (!tfx_AddStageToRibbonBatch(&batch, stage, &ribbon_frame->stage_offsets) || batch.totals.index_offset == 0) {
		return;
	}
	ribbon_frame->stage = stage;
	ribbon_frame->sizes = tfx_GetRibbonBatchSizes(&batch, sizeof(tfx_ribbon_vertex_t));
	ribbon_frame->segments = staging->segments;
	ribbon_frame->ribbons = staging->ribbons;
	ribbon_frame->emitters = staging->emitters;
}

void tfx_impl_GatherFrame(tfx_stage stage, tfxU32 ticks, tfx_impl_ribbon_staging_t *ribbon_staging, bool fetch_gpu_particles, tfx_impl_frame_t *frame) {
	memset(frame, 0, sizeof(tfx_impl_frame_t));
	frame->ticks = ticks;
	if (ticks == 0) {
		return;
	}
	frame->instances = tfx_GetInstanceBuffer(stage);
	frame->instance_count = (tfxU32)tfx_GetInstanceCount(stage);
	for (tfxU32 layer = 0; layer != tfxLAYERS; ++layer) {
		frame->layer_start[layer] = (tfxU32)(tfx_GetInstanceBufferByLayer(stage, layer) - frame->instances);
		frame->layer_count[layer] = (tfxU32)tfx_GetInstanceCountByLayer(stage, layer);
	}
	if (ribbon_staging) {
		tfx_impl__gather_ribbons(stage, ribbon_staging, &frame->ribbons);
	}
	if (fetch_gpu_particles) {
		frame->has_gpu_particles = tfx_GetGPUParticleFrame(stage, &frame->gpu_particles);
		frame->gpu_pool_sizes = tfx_GetGPUParticleBufferSizes(stage);
	}
}

void tfx_impl_GatherBakedFrame(tfx_animation_manager animation_manager, tfx_impl_baked_frame_t *baked) {
	memset(baked, 0, sizeof(tfx_impl_baked_frame_t));
	tfx_animation_buffer_metrics_t metrics = tfx_GetAnimationBufferMetrics(animation_manager);
	baked->sprite_count = metrics.total_sprites_to_draw;
	baked->animation_instance_count = metrics.instances_size;
	baked->animated_shapes = tfx_HasAnimatedShapes(animation_manager);
	baked->offsets = tfx_GetOffsetsBufferPointer(animation_manager);
	baked->offsets_bytes = tfx_GetOffsetsSizeInBytes(animation_manager);
	baked->animation_instances = tfx_GetAnimationInstancesBufferPointer(animation_manager);
	baked->animation_instances_bytes = tfx_GetAnimationInstancesSizeInBytes(animation_manager);
	baked->animation_manager = animation_manager;
	baked->ribbon_count = metrics.total_ribbons_to_draw;
	baked->ribbon_offsets = tfx_GetRibbonOffsetsBufferPointer(animation_manager);
	baked->ribbon_offsets_bytes = tfx_GetRibbonOffsetsSizeInBytes(animation_manager);
}

void tfx_impl_FreeRibbonStaging(tfx_impl_ribbon_staging_t *ribbon_staging) {
	free(ribbon_staging->segments);
	free(ribbon_staging->ribbons);
	free(ribbon_staging->emitters);
	memset(ribbon_staging, 0, sizeof(tfx_impl_ribbon_staging_t));
}

void tfx_impl_SetUniforms(tfx_impl_uniforms_t *uniforms, const tfx_impl_camera_t *camera, float screen_width, float screen_height, const tfx_impl_timer_t *timer) {
	memset(uniforms, 0, sizeof(tfx_impl_uniforms_t));
	uniforms->view = camera->view;
	uniforms->projection = camera->projection;
	uniforms->screen_size[0] = screen_width;
	uniforms->screen_size[1] = screen_height;
	uniforms->timer_lerp = timer->lerp;
	uniforms->update_time = (float)timer->update_length;
	uniforms->camera_position[0] = camera->position.x;
	uniforms->camera_position[1] = camera->position.y;
	uniforms->camera_position[2] = camera->position.z;
}
