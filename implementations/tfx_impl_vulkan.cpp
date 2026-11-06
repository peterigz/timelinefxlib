#include "tfx_impl_vulkan.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

#define TFX_IMPL_VULKAN_INITIAL_INSTANCE_CAPACITY 10000
#define TFX_IMPL_VULKAN_BUFFER_DESCRIPTOR_CAPACITY 1024
#define TFX_IMPL_VULKAN_IMAGE_DESCRIPTOR_CAPACITY 64
//Bigger uploads get a buffer of their own rather than growing a frame's ring
#define TFX_IMPL_VULKAN_RING_SIZE (4 * 1024 * 1024)
//tfx_GetGlobalColorRampLayerCount never goes past this
#define TFX_IMPL_VULKAN_MAX_COLOR_RAMP_LAYERS 256
#define TFX_IMPL_VULKAN_VIEW_UNIFORMS_SIZE 256
//Match the thread group sizes in timelinefx.slang
#define TFX_IMPL_VULKAN_KERNEL_GROUP_SIZE 256
#define TFX_IMPL_VULKAN_RIBBON_GROUP_SIZE 256
#define TFX_IMPL_VULKAN_PLAYBACK_GROUP_SIZE 128
//The bindings of the descriptor arrays in timelinefx.slang, all in set 0
#define TFX_IMPL_VULKAN_SAMPLER_BINDING 0
#define TFX_IMPL_VULKAN_IMAGE_BINDING 3
#define TFX_IMPL_VULKAN_BUFFER_BINDING 5
#define TFX_IMPL_VULKAN_UNIFORM_BINDING 7
#define TFX_IMPL_VULKAN_SHADER_STAGES (VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT | VK_SHADER_STAGE_COMPUTE_BIT)
//Every buffer is a storage buffer that can be copied to, and copied from when it grows
#define TFX_IMPL_VULKAN_BUFFER_USAGE (VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT)
#define TFX_IMPL_VULKAN_INSTANCE_STRIDE ((tfxU32)sizeof(tfx_instance_t))

//timelinefx.slang mirrors these structs by hand, in its [Render_rows] section
static_assert(sizeof(tfx_instance_t) == 48, "tfx_instance_t changed, update tfx_instance in timelinefx.slang");
static_assert(sizeof(tfx_gpu_image_data_t) == 48, "tfx_gpu_image_data_t changed, update tfx_image_data in timelinefx.slang");
static_assert(sizeof(tfx_gpu_particle_properties_t) == 32, "tfx_gpu_particle_properties_t changed, update tfx_emitter_properties in timelinefx.slang");
static_assert(sizeof(tfx_ribbon_vertex_t) == 32, "tfx_ribbon_vertex_t changed, update tfx_ribbon_vertex_row in timelinefx.slang");
static_assert(sizeof(tfx_gpu_graph_data_t) == 48, "tfx_gpu_graph_data_t changed, update tfx_gpu_graph_data in timelinefx.slang");
static_assert(sizeof(tfx_impl_uniforms_t) <= TFX_IMPL_VULKAN_VIEW_UNIFORMS_SIZE, "The view uniforms have to fit their uniform buffer");

//tfx_billboard_push in timelinefx.slang, which is the host's to declare. Every index is a slot in a descriptor array
typedef struct tfx_impl_vulkan_billboard_push_s {
	tfxU32 color_ramp_texture_index;
	tfxU32 image_data_index;
	tfxU32 properties_index;
	tfxU32 sampler_index;
	tfxU32 prev_billboards_index;
	tfxU32 index_offset;
	tfxU32 uniform_index;
} tfx_impl_vulkan_billboard_push_t;

//tfx_baked_billboard_push, the live push's fields and then the clip space offsets, so the fragment shader reads either
typedef struct tfx_impl_vulkan_baked_billboard_push_s {
	tfx_impl_vulkan_billboard_push_t billboard;
	float ndc_offset_x;
	float ndc_offset_y;
} tfx_impl_vulkan_baked_billboard_push_t;

//One push constant range serves every pipeline, as big as the largest block
#define TFX_IMPL_VULKAN_PUSH_CONSTANTS_SIZE ((tfxU32)sizeof(tfx_ribbon_bucket_globals_t))
static_assert(TFX_IMPL_VULKAN_PUSH_CONSTANTS_SIZE <= 128, "Every Vulkan device has at least 128 bytes of push constants");
static_assert(sizeof(tfx_gpu_particle_push_t) <= TFX_IMPL_VULKAN_PUSH_CONSTANTS_SIZE && sizeof(tfx_sprite_data_push_t) <= TFX_IMPL_VULKAN_PUSH_CONSTANTS_SIZE
	&& sizeof(tfx_ribbon_data_push_t) <= TFX_IMPL_VULKAN_PUSH_CONSTANTS_SIZE && sizeof(tfx_impl_vulkan_baked_billboard_push_t) <= TFX_IMPL_VULKAN_PUSH_CONSTANTS_SIZE,
	"A push constant block outgrew the push constant range");

//--- Vulkan_resources

//A frame still in flight may use it, so it goes when this frame slot comes around again
static void tfx_impl_vulkan__defer_release(tfx_impl_vulkan_t *renderer, VkBuffer buffer, VkImage image, VkImageView view, VkDeviceMemory memory, tfxU32 buffer_index, tfxU32 image_index) {
	tfx_impl_vulkan_release_list_t *list = &renderer->release_lists[renderer->frame_index];
	if (list->count == list->capacity) {
		tfxU32 capacity = list->capacity ? list->capacity * 2 : 16;
		tfx_impl_vulkan_release_t *grown = (tfx_impl_vulkan_release_t *)realloc(list->releases, sizeof(tfx_impl_vulkan_release_t) * capacity);
		if (!grown) {
			return;
		}
		list->releases = grown;
		list->capacity = capacity;
	}
	tfx_impl_vulkan_release_t *release = &list->releases[list->count++];
	release->buffer = buffer;
	release->image = image;
	release->view = view;
	release->memory = memory;
	release->buffer_index = buffer_index;
	release->image_index = image_index;
}

static void tfx_impl_vulkan__free_index(tfx_impl_vulkan_index_list_t *list, tfxU32 index) {
	if (index != TFX_IMPL_VULKAN_NO_DESCRIPTOR && list->indexes) {
		list->indexes[list->count++] = index;
	}
}

static void tfx_impl_vulkan__process_release_list(tfx_impl_vulkan_t *renderer, tfx_impl_vulkan_release_list_t *list) {
	for (tfxU32 index = 0; index != list->count; ++index) {
		tfx_impl_vulkan_release_t *release = &list->releases[index];
		vkDestroyImageView(renderer->device, release->view, NULL);
		vkDestroyImage(renderer->device, release->image, NULL);
		vkDestroyBuffer(renderer->device, release->buffer, NULL);
		vkFreeMemory(renderer->device, release->memory, NULL);
		tfx_impl_vulkan__free_index(&renderer->free_buffer_indexes, release->buffer_index);
		tfx_impl_vulkan__free_index(&renderer->free_image_indexes, release->image_index);
	}
	list->count = 0;
}

//Every index from first_free up is free, handed out lowest first
static bool tfx_impl_vulkan__create_index_list(tfx_impl_vulkan_index_list_t *list, tfxU32 capacity, tfxU32 first_free) {
	list->indexes = (tfxU32 *)malloc(sizeof(tfxU32) * capacity);
	if (!list->indexes) {
		return false;
	}
	for (tfxU32 index = capacity; index != first_free; --index) {
		list->indexes[list->count++] = index - 1;
	}
	return true;
}

static tfxU32 tfx_impl_vulkan__allocate_index(tfx_impl_vulkan_index_list_t *list, const char *array_name) {
	if (list->count == 0) {
		tfx_Print("Every slot in the %s descriptor array is in use", array_name);
		return TFX_IMPL_VULKAN_NO_DESCRIPTOR;
	}
	return list->indexes[--list->count];
}

static VkDeviceMemory tfx_impl_vulkan__allocate_memory(tfx_impl_vulkan_t *renderer, const VkMemoryRequirements *requirements, VkMemoryPropertyFlags properties) {
	for (tfxU32 type = 0; type != renderer->memory_properties.memoryTypeCount; ++type) {
		if ((requirements->memoryTypeBits & (1u << type)) && (renderer->memory_properties.memoryTypes[type].propertyFlags & properties) == properties) {
			VkMemoryAllocateInfo allocate_info = { VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO };
			allocate_info.allocationSize = requirements->size;
			allocate_info.memoryTypeIndex = type;
			VkDeviceMemory memory = VK_NULL_HANDLE;
			VkResult result = vkAllocateMemory(renderer->device, &allocate_info, NULL, &memory);
			if (result != VK_SUCCESS) {
				tfx_Print("Could not allocate %llu bytes of device memory (%i)", (unsigned long long)requirements->size, (int)result);
			}
			return memory;
		}
	}
	tfx_Print("No memory type has the properties 0x%x", (unsigned int)properties);
	return VK_NULL_HANDLE;
}

//A buffer with its own memory. An engine would sub-allocate from its own allocator or VMA instead
static bool tfx_impl_vulkan__create_raw_buffer(tfx_impl_vulkan_t *renderer, VkDeviceSize size_in_bytes, VkBufferUsageFlags usage, VkMemoryPropertyFlags properties, VkBuffer *buffer, VkDeviceMemory *memory) {
	VkBufferCreateInfo create_info = { VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO };
	create_info.size = size_in_bytes;
	create_info.usage = usage;
	create_info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
	*buffer = VK_NULL_HANDLE;
	*memory = VK_NULL_HANDLE;
	VkResult result = vkCreateBuffer(renderer->device, &create_info, NULL, buffer);
	if (result != VK_SUCCESS) {
		tfx_Print("Could not create a buffer of %llu bytes (%i)", (unsigned long long)size_in_bytes, (int)result);
		*buffer = VK_NULL_HANDLE;
		return false;
	}
	VkMemoryRequirements requirements;
	vkGetBufferMemoryRequirements(renderer->device, *buffer, &requirements);
	*memory = tfx_impl_vulkan__allocate_memory(renderer, &requirements, properties);
	if (!*memory || vkBindBufferMemory(renderer->device, *buffer, *memory, 0) != VK_SUCCESS) {
		vkDestroyBuffer(renderer->device, *buffer, NULL);
		vkFreeMemory(renderer->device, *memory, NULL);
		*buffer = VK_NULL_HANDLE;
		*memory = VK_NULL_HANDLE;
		return false;
	}
	return true;
}

static bool tfx_impl_vulkan__create_mapped_buffer(tfx_impl_vulkan_t *renderer, VkDeviceSize size_in_bytes, VkBufferUsageFlags usage, tfx_impl_vulkan_mapped_buffer_t *mapped_buffer) {
	memset(mapped_buffer, 0, sizeof(tfx_impl_vulkan_mapped_buffer_t));
	if (!tfx_impl_vulkan__create_raw_buffer(renderer, size_in_bytes, usage, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT, &mapped_buffer->buffer, &mapped_buffer->memory)) {
		return false;
	}
	if (vkMapMemory(renderer->device, mapped_buffer->memory, 0, VK_WHOLE_SIZE, 0, (void **)&mapped_buffer->mapped) != VK_SUCCESS) {
		tfx_Print("Could not map a buffer of %llu bytes", (unsigned long long)size_in_bytes);
		vkDestroyBuffer(renderer->device, mapped_buffer->buffer, NULL);
		vkFreeMemory(renderer->device, mapped_buffer->memory, NULL);
		memset(mapped_buffer, 0, sizeof(tfx_impl_vulkan_mapped_buffer_t));
		return false;
	}
	mapped_buffer->capacity = size_in_bytes;
	return true;
}

static void tfx_impl_vulkan__release_mapped_buffer(tfx_impl_vulkan_t *renderer, tfx_impl_vulkan_mapped_buffer_t *mapped_buffer) {
	if (mapped_buffer->buffer) {
		tfx_impl_vulkan__defer_release(renderer, mapped_buffer->buffer, VK_NULL_HANDLE, VK_NULL_HANDLE, mapped_buffer->memory, TFX_IMPL_VULKAN_NO_DESCRIPTOR, TFX_IMPL_VULKAN_NO_DESCRIPTOR);
	}
	memset(mapped_buffer, 0, sizeof(tfx_impl_vulkan_mapped_buffer_t));
}

static void tfx_impl_vulkan__write_buffer_descriptor(tfx_impl_vulkan_t *renderer, tfxU32 binding, VkDescriptorType type, tfxU32 index, VkBuffer buffer) {
	VkDescriptorBufferInfo buffer_info = { buffer, 0, VK_WHOLE_SIZE };
	VkWriteDescriptorSet write = { VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET };
	write.dstSet = renderer->descriptor_set;
	write.dstBinding = binding;
	write.dstArrayElement = index;
	write.descriptorCount = 1;
	write.descriptorType = type;
	write.pBufferInfo = &buffer_info;
	vkUpdateDescriptorSets(renderer->device, 1, &write, 0, NULL);
}

static void tfx_impl_vulkan__release_buffer(tfx_impl_vulkan_t *renderer, tfx_impl_vulkan_buffer_t *buffer) {
	if (buffer->buffer) {
		tfx_impl_vulkan__defer_release(renderer, buffer->buffer, VK_NULL_HANDLE, VK_NULL_HANDLE, buffer->memory, buffer->index, TFX_IMPL_VULKAN_NO_DESCRIPTOR);
	}
	buffer->buffer = VK_NULL_HANDLE;
	buffer->memory = VK_NULL_HANDLE;
	buffer->capacity = 0;
	buffer->index = TFX_IMPL_VULKAN_NO_DESCRIPTOR;
}

//A device local storage buffer and its descriptor. extra_usage makes it a vertex or index buffer as well
static bool tfx_impl_vulkan__create_buffer(tfx_impl_vulkan_t *renderer, VkDeviceSize size_in_bytes, VkBufferUsageFlags extra_usage, tfx_impl_vulkan_buffer_t *buffer) {
	tfx_impl_vulkan__release_buffer(renderer, buffer);
	size_in_bytes = ((size_in_bytes ? size_in_bytes : 16) + 255) & ~(VkDeviceSize)255;
	if (!tfx_impl_vulkan__create_raw_buffer(renderer, size_in_bytes, TFX_IMPL_VULKAN_BUFFER_USAGE | extra_usage, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, &buffer->buffer, &buffer->memory)) {
		return false;
	}
	buffer->capacity = size_in_bytes;
	buffer->index = tfx_impl_vulkan__allocate_index(&renderer->free_buffer_indexes, "storage buffer");
	if (buffer->index == TFX_IMPL_VULKAN_NO_DESCRIPTOR) {
		tfx_impl_vulkan__release_buffer(renderer, buffer);
		return false;
	}
	tfx_impl_vulkan__write_buffer_descriptor(renderer, TFX_IMPL_VULKAN_BUFFER_BINDING, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, buffer->index, buffer->buffer);
	return true;
}

//Grows by half again whenever it's too small, so the size settles after the first few frames. Never empty, so a shader always has a buffer to read
static bool tfx_impl_vulkan__ensure_buffer(tfx_impl_vulkan_t *renderer, size_t size_in_bytes, VkBufferUsageFlags extra_usage, tfx_impl_vulkan_buffer_t *buffer) {
	VkDeviceSize size = size_in_bytes ? size_in_bytes : 16;
	if (buffer->buffer && buffer->capacity >= size) {
		return true;
	}
	return tfx_impl_vulkan__create_buffer(renderer, size + size / 2, extra_usage, buffer);
}

static void tfx_impl_vulkan__memory_barrier(VkCommandBuffer command_buffer, VkPipelineStageFlags2 source_stages, VkAccessFlags2 source_access, VkPipelineStageFlags2 destination_stages, VkAccessFlags2 destination_access) {
	VkMemoryBarrier2 barrier = { VK_STRUCTURE_TYPE_MEMORY_BARRIER_2 };
	barrier.srcStageMask = source_stages;
	barrier.srcAccessMask = source_access;
	barrier.dstStageMask = destination_stages;
	barrier.dstAccessMask = destination_access;
	VkDependencyInfo dependency = { VK_STRUCTURE_TYPE_DEPENDENCY_INFO };
	dependency.memoryBarrierCount = 1;
	dependency.pMemoryBarriers = &barrier;
	vkCmdPipelineBarrier2(command_buffer, &dependency);
}

//The copies that follow overwrite buffers the last frame's passes read and wrote, so everything recorded before them finishes first
static void tfx_impl_vulkan__begin_transfers(VkCommandBuffer command_buffer) {
	tfx_impl_vulkan__memory_barrier(command_buffer, VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, VK_ACCESS_2_MEMORY_WRITE_BIT,
		VK_PIPELINE_STAGE_2_ALL_TRANSFER_BIT, VK_ACCESS_2_TRANSFER_READ_BIT | VK_ACCESS_2_TRANSFER_WRITE_BIT);
}

//Images only change layout around their uploads: into a copy destination, from there to a blit source for the next mip, and back to shader reads
static void tfx_impl_vulkan__image_barrier(VkCommandBuffer command_buffer, VkImage image, VkImageLayout old_layout, VkImageLayout new_layout,
	tfxU32 first_mip, tfxU32 mip_count, tfxU32 first_layer, tfxU32 layer_count) {
	VkImageMemoryBarrier2 barrier = { VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER_2 };
	if (new_layout == VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL) {
		barrier.srcStageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
		barrier.srcAccessMask = VK_ACCESS_2_MEMORY_WRITE_BIT;
		barrier.dstStageMask = VK_PIPELINE_STAGE_2_ALL_TRANSFER_BIT;
		barrier.dstAccessMask = VK_ACCESS_2_TRANSFER_WRITE_BIT;
	} else if (new_layout == VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL) {
		barrier.srcStageMask = VK_PIPELINE_STAGE_2_ALL_TRANSFER_BIT;
		barrier.srcAccessMask = VK_ACCESS_2_TRANSFER_WRITE_BIT;
		barrier.dstStageMask = VK_PIPELINE_STAGE_2_ALL_TRANSFER_BIT;
		barrier.dstAccessMask = VK_ACCESS_2_TRANSFER_READ_BIT;
	} else {
		barrier.srcStageMask = VK_PIPELINE_STAGE_2_ALL_TRANSFER_BIT;
		barrier.srcAccessMask = VK_ACCESS_2_TRANSFER_WRITE_BIT;
		barrier.dstStageMask = VK_PIPELINE_STAGE_2_VERTEX_SHADER_BIT | VK_PIPELINE_STAGE_2_FRAGMENT_SHADER_BIT | VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT;
		barrier.dstAccessMask = VK_ACCESS_2_SHADER_SAMPLED_READ_BIT;
	}
	barrier.oldLayout = old_layout;
	barrier.newLayout = new_layout;
	barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
	barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
	barrier.image = image;
	barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
	barrier.subresourceRange.baseMipLevel = first_mip;
	barrier.subresourceRange.levelCount = mip_count;
	barrier.subresourceRange.baseArrayLayer = first_layer;
	barrier.subresourceRange.layerCount = layer_count;
	VkDependencyInfo dependency = { VK_STRUCTURE_TYPE_DEPENDENCY_INFO };
	dependency.imageMemoryBarrierCount = 1;
	dependency.pImageMemoryBarriers = &barrier;
	vkCmdPipelineBarrier2(command_buffer, &dependency);
}

//Upload memory for this frame. Small uploads share the frame slot's ring, big ones get a buffer of their own that goes when the slot comes round
static unsigned char *tfx_impl_vulkan__allocate_upload(tfx_impl_vulkan_t *renderer, VkDeviceSize size_in_bytes, VkDeviceSize alignment, VkBuffer *source, VkDeviceSize *offset) {
	tfx_impl_vulkan_mapped_buffer_t *ring = &renderer->upload_rings[renderer->frame_index];
	VkDeviceSize start = (ring->used + alignment - 1) / alignment * alignment;
	if (!ring->buffer || start + size_in_bytes > ring->capacity) {
		if (size_in_bytes > TFX_IMPL_VULKAN_RING_SIZE / 4) {
			tfx_impl_vulkan_mapped_buffer_t dedicated;
			if (!tfx_impl_vulkan__create_mapped_buffer(renderer, size_in_bytes, VK_BUFFER_USAGE_TRANSFER_SRC_BIT, &dedicated)) {
				return NULL;
			}
			unsigned char *mapped = dedicated.mapped;
			*source = dedicated.buffer;
			*offset = 0;
			tfx_impl_vulkan__release_mapped_buffer(renderer, &dedicated);
			return mapped;
		}
		VkDeviceSize capacity = ring->capacity ? ring->capacity * 2 : TFX_IMPL_VULKAN_RING_SIZE;
		tfx_impl_vulkan__release_mapped_buffer(renderer, ring);
		if (!tfx_impl_vulkan__create_mapped_buffer(renderer, capacity, VK_BUFFER_USAGE_TRANSFER_SRC_BIT, ring)) {
			return NULL;
		}
		start = 0;
	}
	ring->used = start + size_in_bytes;
	*source = ring->buffer;
	*offset = start;
	return ring->mapped + start;
}

static bool tfx_impl_vulkan__upload_buffer(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, tfx_impl_vulkan_buffer_t *buffer, const void *data, size_t size_in_bytes) {
	if (size_in_bytes == 0) {
		return true;
	}
	VkBuffer source = VK_NULL_HANDLE;
	VkDeviceSize offset = 0;
	unsigned char *mapped = tfx_impl_vulkan__allocate_upload(renderer, size_in_bytes, 16, &source, &offset);
	if (!mapped) {
		return false;
	}
	memcpy(mapped, data, size_in_bytes);
	VkBufferCopy region = { offset, 0, size_in_bytes };
	vkCmdCopyBuffer(command_buffer, source, buffer->buffer, 1, &region);
	return true;
}

static bool tfx_impl_vulkan__write_buffer(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, tfx_impl_vulkan_buffer_t *buffer, const void *data, size_t size_in_bytes) {
	return tfx_impl_vulkan__ensure_buffer(renderer, size_in_bytes, 0, buffer) && tfx_impl_vulkan__upload_buffer(renderer, command_buffer, buffer, data, size_in_bytes);
}

static void tfx_impl_vulkan__destroy_image(tfx_impl_vulkan_t *renderer, tfx_impl_vulkan_image_t *image) {
	vkDestroyImageView(renderer->device, image->view, NULL);
	vkDestroyImage(renderer->device, image->image, NULL);
	vkFreeMemory(renderer->device, image->memory, NULL);
	memset(image, 0, sizeof(tfx_impl_vulkan_image_t));
	image->index = TFX_IMPL_VULKAN_NO_DESCRIPTOR;
}

static void tfx_impl_vulkan__release_image(tfx_impl_vulkan_t *renderer, tfx_impl_vulkan_image_t *image) {
	if (image->image) {
		tfx_impl_vulkan__defer_release(renderer, VK_NULL_HANDLE, image->image, image->view, image->memory, TFX_IMPL_VULKAN_NO_DESCRIPTOR, image->index);
	}
	memset(image, 0, sizeof(tfx_impl_vulkan_image_t));
	image->index = TFX_IMPL_VULKAN_NO_DESCRIPTOR;
}

//A sampled texture array and its view, in undefined layout until its first upload. Its descriptor is written once it has an index
static bool tfx_impl_vulkan__create_image(tfx_impl_vulkan_t *renderer, tfxU32 width, tfxU32 height, tfxU32 layer_count, tfxU32 mip_count, VkFormat format, VkImageUsageFlags extra_usage, tfx_impl_vulkan_image_t *image) {
	memset(image, 0, sizeof(tfx_impl_vulkan_image_t));
	image->index = TFX_IMPL_VULKAN_NO_DESCRIPTOR;
	VkImageCreateInfo create_info = { VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO };
	create_info.imageType = VK_IMAGE_TYPE_2D;
	create_info.format = format;
	create_info.extent.width = width;
	create_info.extent.height = height;
	create_info.extent.depth = 1;
	create_info.mipLevels = mip_count;
	create_info.arrayLayers = layer_count;
	create_info.samples = VK_SAMPLE_COUNT_1_BIT;
	create_info.tiling = VK_IMAGE_TILING_OPTIMAL;
	create_info.usage = VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT | extra_usage;
	create_info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
	create_info.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
	VkResult result = vkCreateImage(renderer->device, &create_info, NULL, &image->image);
	if (result != VK_SUCCESS) {
		tfx_Print("Could not create a %u x %u x %u texture array (%i)", width, height, layer_count, (int)result);
		image->image = VK_NULL_HANDLE;
		return false;
	}
	VkMemoryRequirements requirements;
	vkGetImageMemoryRequirements(renderer->device, image->image, &requirements);
	image->memory = tfx_impl_vulkan__allocate_memory(renderer, &requirements, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
	VkImageViewCreateInfo view_info = { VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO };
	view_info.image = image->image;
	view_info.viewType = VK_IMAGE_VIEW_TYPE_2D_ARRAY;
	view_info.format = format;
	view_info.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
	view_info.subresourceRange.levelCount = mip_count;
	view_info.subresourceRange.layerCount = layer_count;
	if (!image->memory || vkBindImageMemory(renderer->device, image->image, image->memory, 0) != VK_SUCCESS
		|| vkCreateImageView(renderer->device, &view_info, NULL, &image->view) != VK_SUCCESS) {
		tfx_Print("Could not create the memory and view of a %u x %u x %u texture array", width, height, layer_count);
		tfx_impl_vulkan__destroy_image(renderer, image);
		return false;
	}
	image->layout = VK_IMAGE_LAYOUT_UNDEFINED;
	image->width = width;
	image->height = height;
	image->layer_count = layer_count;
	image->mip_count = mip_count;
	return true;
}

static void tfx_impl_vulkan__write_image_descriptor(tfx_impl_vulkan_t *renderer, const tfx_impl_vulkan_image_t *image) {
	VkDescriptorImageInfo image_info = { VK_NULL_HANDLE, image->view, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL };
	VkWriteDescriptorSet write = { VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET };
	write.dstSet = renderer->descriptor_set;
	write.dstBinding = TFX_IMPL_VULKAN_IMAGE_BINDING;
	write.dstArrayElement = image->index;
	write.descriptorCount = 1;
	write.descriptorType = VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE;
	write.pImageInfo = &image_info;
	vkUpdateDescriptorSets(renderer->device, 1, &write, 0, NULL);
}

//Tightly packed layers of mip 0 from the frame's upload memory, into an image in the copy destination layout
static bool tfx_impl_vulkan__upload_image_layers(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_vulkan_image_t *image, tfxU32 first_layer, tfxU32 layer_count,
	tfxU32 bytes_per_pixel, const void *pixels) {
	VkDeviceSize size_in_bytes = (VkDeviceSize)image->width * image->height * bytes_per_pixel * layer_count;
	VkBuffer source = VK_NULL_HANDLE;
	VkDeviceSize offset = 0;
	unsigned char *mapped = tfx_impl_vulkan__allocate_upload(renderer, size_in_bytes, 16, &source, &offset);
	if (!mapped) {
		return false;
	}
	memcpy(mapped, pixels, size_in_bytes);
	VkBufferImageCopy region = {};
	region.bufferOffset = offset;
	region.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
	region.imageSubresource.baseArrayLayer = first_layer;
	region.imageSubresource.layerCount = layer_count;
	region.imageExtent.width = image->width;
	region.imageExtent.height = image->height;
	region.imageExtent.depth = 1;
	vkCmdCopyBufferToImage(command_buffer, source, image->image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 1, &region);
	return true;
}

static VkFormat tfx_impl_vulkan__color_ramp_format(tfxU32 bytes_per_pixel) {
	if (bytes_per_pixel == 4) {
		return VK_FORMAT_R8G8B8A8_UNORM;
	} else if (bytes_per_pixel == 16) {
		return VK_FORMAT_R32G32B32A32_SFLOAT;
	}
	return VK_FORMAT_R16G16B16A16_SFLOAT;
}

static tfxU32 tfx_impl_vulkan__bitmap_bytes_per_pixel(tfx_bitmap_t *bitmap) {
	return (tfxU32)(tfx_GetBitmapSize(bitmap) / ((size_t)tfx_GetBitmapWidth(bitmap) * tfx_GetBitmapHeight(bitmap)));
}

//A color ramp texture array and its descriptor, replacing the old ones once the frames using them are done
static bool tfx_impl_vulkan__replace_color_ramps(tfx_impl_vulkan_t *renderer, tfx_impl_vulkan_image_t *texture, tfx_bitmap_t *first_bitmap, tfxU32 layer_count) {
	tfx_impl_vulkan__release_image(renderer, texture);
	VkFormat format = tfx_impl_vulkan__color_ramp_format(tfx_impl_vulkan__bitmap_bytes_per_pixel(first_bitmap));
	if (!tfx_impl_vulkan__create_image(renderer, (tfxU32)tfx_GetBitmapWidth(first_bitmap), (tfxU32)tfx_GetBitmapHeight(first_bitmap), layer_count, 1, format, 0, texture)) {
		return false;
	}
	texture->index = tfx_impl_vulkan__allocate_index(&renderer->free_image_indexes, "image");
	if (texture->index == TFX_IMPL_VULKAN_NO_DESCRIPTOR) {
		tfx_impl_vulkan__destroy_image(renderer, texture);
		return false;
	}
	tfx_impl_vulkan__write_image_descriptor(renderer, texture);
	return true;
}

//Every layer goes to the copy destination layout and back, so the layers that aren't uploaded keep what they hold
static void tfx_impl_vulkan__begin_image_upload(VkCommandBuffer command_buffer, tfx_impl_vulkan_image_t *texture) {
	tfx_impl_vulkan__image_barrier(command_buffer, texture->image, texture->layout, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 0, texture->mip_count, 0, texture->layer_count);
}

static void tfx_impl_vulkan__end_image_upload(VkCommandBuffer command_buffer, tfx_impl_vulkan_image_t *texture) {
	tfx_impl_vulkan__image_barrier(command_buffer, texture->image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, 0, texture->mip_count, 0, texture->layer_count);
	texture->layout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
}

static bool tfx_impl_vulkan__upload_color_ramp_layer(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_vulkan_image_t *texture, tfxU32 layer, tfx_bitmap_t *bitmap) {
	return tfx_impl_vulkan__upload_image_layers(renderer, command_buffer, texture, layer, 1, tfx_impl_vulkan__bitmap_bytes_per_pixel(bitmap), tfx_GetBitmapData(bitmap));
}

//slangc compiles every shader ahead of time to <entry point>.spv in the shader folder
static VkShaderModule tfx_impl_vulkan__load_shader(tfx_impl_vulkan_t *renderer, const char *shader_folder, const char *entry_point) {
	char path[1024];
	snprintf(path, sizeof(path), "%s/%s.spv", shader_folder, entry_point);
	FILE *file = fopen(path, "rb");
	if (!file) {
		tfx_Print("Could not open the shader %s", path);
		return VK_NULL_HANDLE;
	}
	fseek(file, 0, SEEK_END);
	long file_size = ftell(file);
	fseek(file, 0, SEEK_SET);
	uint32_t *code = file_size > 0 && file_size % 4 == 0 ? (uint32_t *)malloc((size_t)file_size) : NULL;
	if (code && fread(code, 1, (size_t)file_size, file) != (size_t)file_size) {
		free(code);
		code = NULL;
	}
	fclose(file);
	if (!code) {
		tfx_Print("Could not read the shader %s", path);
		return VK_NULL_HANDLE;
	}
	VkShaderModuleCreateInfo create_info = { VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO };
	create_info.codeSize = (size_t)file_size;
	create_info.pCode = code;
	VkShaderModule module = VK_NULL_HANDLE;
	if (vkCreateShaderModule(renderer->device, &create_info, NULL, &module) != VK_SUCCESS) {
		tfx_Print("Could not create the shader module %s", path);
	}
	free(code);
	return module;
}

static VkPipeline tfx_impl_vulkan__create_compute_pipeline(tfx_impl_vulkan_t *renderer, const char *shader_folder, const char *entry_point) {
	VkShaderModule module = tfx_impl_vulkan__load_shader(renderer, shader_folder, entry_point);
	if (!module) {
		return VK_NULL_HANDLE;
	}
	VkComputePipelineCreateInfo create_info = { VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO };
	create_info.stage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
	create_info.stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
	create_info.stage.module = module;
	create_info.stage.pName = entry_point;
	create_info.layout = renderer->pipeline_layout;
	VkPipeline pipeline = VK_NULL_HANDLE;
	VkResult result = vkCreateComputePipelines(renderer->device, VK_NULL_HANDLE, 1, &create_info, NULL, &pipeline);
	if (result != VK_SUCCESS) {
		tfx_Print("Could not create the %s pipeline (%i)", entry_point, (int)result);
	}
	vkDestroyShaderModule(renderer->device, module, NULL);
	return pipeline;
}

//Premultiplied blending, depth tested but not written, or for a layer that writes depth, blending off and depth written. No culling, and
//the viewport and scissor left to the caller
static VkPipeline tfx_impl_vulkan__create_graphics_pipeline(tfx_impl_vulkan_t *renderer, const char *shader_folder, const char *vertex_entry_point, const char *fragment_entry_point,
	const VkVertexInputBindingDescription *vertex_binding, const VkVertexInputAttributeDescription *attributes, tfxU32 attribute_count, VkFormat color_format, VkFormat depth_format, bool writes_depth) {
	VkShaderModule vertex_module = tfx_impl_vulkan__load_shader(renderer, shader_folder, vertex_entry_point);
	VkShaderModule fragment_module = tfx_impl_vulkan__load_shader(renderer, shader_folder, fragment_entry_point);
	VkPipeline pipeline = VK_NULL_HANDLE;
	if (vertex_module && fragment_module) {
		VkPipelineShaderStageCreateInfo stages[2] = {};
		stages[0].sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
		stages[0].stage = VK_SHADER_STAGE_VERTEX_BIT;
		stages[0].module = vertex_module;
		stages[0].pName = vertex_entry_point;
		stages[1].sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
		stages[1].stage = VK_SHADER_STAGE_FRAGMENT_BIT;
		stages[1].module = fragment_module;
		stages[1].pName = fragment_entry_point;
		VkPipelineVertexInputStateCreateInfo vertex_input = { VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO };
		vertex_input.vertexBindingDescriptionCount = 1;
		vertex_input.pVertexBindingDescriptions = vertex_binding;
		vertex_input.vertexAttributeDescriptionCount = attribute_count;
		vertex_input.pVertexAttributeDescriptions = attributes;
		VkPipelineInputAssemblyStateCreateInfo input_assembly = { VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO };
		input_assembly.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
		VkPipelineViewportStateCreateInfo viewport = { VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO };
		viewport.viewportCount = 1;
		viewport.scissorCount = 1;
		VkPipelineRasterizationStateCreateInfo rasterization = { VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO };
		rasterization.polygonMode = VK_POLYGON_MODE_FILL;
		rasterization.cullMode = VK_CULL_MODE_NONE;
		rasterization.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;
		rasterization.lineWidth = 1.f;
		VkPipelineMultisampleStateCreateInfo multisample = { VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO };
		multisample.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;
		VkPipelineDepthStencilStateCreateInfo depth_stencil = { VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO };
		depth_stencil.depthTestEnable = VK_TRUE;
		depth_stencil.depthWriteEnable = writes_depth ? VK_TRUE : VK_FALSE;
		depth_stencil.depthCompareOp = VK_COMPARE_OP_LESS_OR_EQUAL;
		VkPipelineColorBlendAttachmentState blend = {};
		blend.blendEnable = writes_depth ? VK_FALSE : VK_TRUE;
		blend.srcColorBlendFactor = VK_BLEND_FACTOR_ONE;
		blend.dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
		blend.colorBlendOp = VK_BLEND_OP_ADD;
		blend.srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE;
		blend.dstAlphaBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
		blend.alphaBlendOp = VK_BLEND_OP_ADD;
		blend.colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT | VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
		VkPipelineColorBlendStateCreateInfo color_blend = { VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO };
		color_blend.attachmentCount = 1;
		color_blend.pAttachments = &blend;
		const VkDynamicState dynamic_states[] = { VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR };
		VkPipelineDynamicStateCreateInfo dynamic_state = { VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO };
		dynamic_state.dynamicStateCount = 2;
		dynamic_state.pDynamicStates = dynamic_states;
		VkPipelineRenderingCreateInfo rendering = { VK_STRUCTURE_TYPE_PIPELINE_RENDERING_CREATE_INFO };
		rendering.colorAttachmentCount = 1;
		rendering.pColorAttachmentFormats = &color_format;
		rendering.depthAttachmentFormat = depth_format;
		VkGraphicsPipelineCreateInfo create_info = { VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO };
		create_info.pNext = &rendering;
		create_info.stageCount = 2;
		create_info.pStages = stages;
		create_info.pVertexInputState = &vertex_input;
		create_info.pInputAssemblyState = &input_assembly;
		create_info.pViewportState = &viewport;
		create_info.pRasterizationState = &rasterization;
		create_info.pMultisampleState = &multisample;
		create_info.pDepthStencilState = &depth_stencil;
		create_info.pColorBlendState = &color_blend;
		create_info.pDynamicState = &dynamic_state;
		create_info.layout = renderer->pipeline_layout;
		VkResult result = vkCreateGraphicsPipelines(renderer->device, VK_NULL_HANDLE, 1, &create_info, NULL, &pipeline);
		if (result != VK_SUCCESS) {
			tfx_Print("Could not create the %s pipeline (%i)", vertex_entry_point, (int)result);
			pipeline = VK_NULL_HANDLE;
		}
	}
	vkDestroyShaderModule(renderer->device, vertex_module, NULL);
	vkDestroyShaderModule(renderer->device, fragment_module, NULL);
	return pipeline;
}

static bool tfx_impl_vulkan__create_descriptor_set(tfx_impl_vulkan_t *renderer) {
	const VkDescriptorBindingFlags written_at_init = VK_DESCRIPTOR_BINDING_PARTIALLY_BOUND_BIT;
	//Images and buffers are added while earlier frames are in flight, into slots those frames don't read
	const VkDescriptorBindingFlags written_any_time = VK_DESCRIPTOR_BINDING_PARTIALLY_BOUND_BIT | VK_DESCRIPTOR_BINDING_UPDATE_AFTER_BIND_BIT | VK_DESCRIPTOR_BINDING_UPDATE_UNUSED_WHILE_PENDING_BIT;
	const VkDescriptorType types[4] = { VK_DESCRIPTOR_TYPE_SAMPLER, VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER };
	const tfxU32 binding_numbers[4] = { TFX_IMPL_VULKAN_SAMPLER_BINDING, TFX_IMPL_VULKAN_IMAGE_BINDING, TFX_IMPL_VULKAN_BUFFER_BINDING, TFX_IMPL_VULKAN_UNIFORM_BINDING };
	const tfxU32 counts[4] = { 1, TFX_IMPL_VULKAN_IMAGE_DESCRIPTOR_CAPACITY, TFX_IMPL_VULKAN_BUFFER_DESCRIPTOR_CAPACITY, TFX_IMPL_VULKAN_MAX_FRAMES_IN_FLIGHT };
	const VkDescriptorBindingFlags binding_flags[4] = { written_at_init, written_any_time, written_any_time, written_at_init };
	VkDescriptorSetLayoutBinding bindings[4] = {};
	VkDescriptorPoolSize pool_sizes[4] = {};
	for (tfxU32 index = 0; index != 4; ++index) {
		bindings[index].binding = binding_numbers[index];
		bindings[index].descriptorType = types[index];
		bindings[index].descriptorCount = counts[index];
		bindings[index].stageFlags = TFX_IMPL_VULKAN_SHADER_STAGES;
		pool_sizes[index].type = types[index];
		pool_sizes[index].descriptorCount = counts[index];
	}
	VkDescriptorSetLayoutBindingFlagsCreateInfo flags_info = { VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_BINDING_FLAGS_CREATE_INFO };
	flags_info.bindingCount = 4;
	flags_info.pBindingFlags = binding_flags;
	VkDescriptorSetLayoutCreateInfo layout_info = { VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO };
	layout_info.pNext = &flags_info;
	layout_info.flags = VK_DESCRIPTOR_SET_LAYOUT_CREATE_UPDATE_AFTER_BIND_POOL_BIT;
	layout_info.bindingCount = 4;
	layout_info.pBindings = bindings;
	VkDescriptorPoolCreateInfo pool_info = { VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO };
	pool_info.flags = VK_DESCRIPTOR_POOL_CREATE_UPDATE_AFTER_BIND_BIT;
	pool_info.maxSets = 1;
	pool_info.poolSizeCount = 4;
	pool_info.pPoolSizes = pool_sizes;
	if (vkCreateDescriptorSetLayout(renderer->device, &layout_info, NULL, &renderer->descriptor_set_layout) != VK_SUCCESS
		|| vkCreateDescriptorPool(renderer->device, &pool_info, NULL, &renderer->descriptor_pool) != VK_SUCCESS) {
		tfx_Print("Could not create the descriptor set layout and pool. Check the device has the descriptor indexing features tfx_impl_vulkan_Init lists");
		return false;
	}
	VkDescriptorSetAllocateInfo allocate_info = { VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO };
	allocate_info.descriptorPool = renderer->descriptor_pool;
	allocate_info.descriptorSetCount = 1;
	allocate_info.pSetLayouts = &renderer->descriptor_set_layout;
	if (vkAllocateDescriptorSets(renderer->device, &allocate_info, &renderer->descriptor_set) != VK_SUCCESS) {
		tfx_Print("Could not allocate the descriptor set");
		return false;
	}

	VkSamplerCreateInfo sampler_info = { VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO };
	sampler_info.magFilter = VK_FILTER_LINEAR;
	sampler_info.minFilter = VK_FILTER_LINEAR;
	sampler_info.mipmapMode = VK_SAMPLER_MIPMAP_MODE_LINEAR;
	sampler_info.addressModeU = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
	sampler_info.addressModeV = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
	sampler_info.addressModeW = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
	sampler_info.maxLod = VK_LOD_CLAMP_NONE;
	if (vkCreateSampler(renderer->device, &sampler_info, NULL, &renderer->sampler) != VK_SUCCESS) {
		tfx_Print("Could not create the sampler");
		return false;
	}
	VkDescriptorImageInfo sampler_descriptor = { renderer->sampler, VK_NULL_HANDLE, VK_IMAGE_LAYOUT_UNDEFINED };
	VkWriteDescriptorSet write = { VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET };
	write.dstSet = renderer->descriptor_set;
	write.dstBinding = TFX_IMPL_VULKAN_SAMPLER_BINDING;
	write.descriptorCount = 1;
	write.descriptorType = VK_DESCRIPTOR_TYPE_SAMPLER;
	write.pImageInfo = &sampler_descriptor;
	vkUpdateDescriptorSets(renderer->device, 1, &write, 0, NULL);
	return true;
}

//--- Vulkan_functions

bool tfx_impl_vulkan_Init(tfx_impl_vulkan_t *renderer, VkPhysicalDevice physical_device, VkDevice device, VkFormat color_format, VkFormat depth_format, const char *shader_folder, tfxU32 frames_in_flight) {
	memset(renderer, 0, sizeof(tfx_impl_vulkan_t));
	renderer->device = device;
	renderer->frames_in_flight = frames_in_flight == 0 ? 1 : frames_in_flight < TFX_IMPL_VULKAN_MAX_FRAMES_IN_FLIGHT ? frames_in_flight : TFX_IMPL_VULKAN_MAX_FRAMES_IN_FLIGHT;
	vkGetPhysicalDeviceMemoryProperties(physical_device, &renderer->memory_properties);
	renderer->shape_atlas.index = TFX_IMPL_VULKAN_NO_DESCRIPTOR;
	renderer->color_ramps.index = TFX_IMPL_VULKAN_NO_DESCRIPTOR;
	renderer->baked_color_ramps.index = TFX_IMPL_VULKAN_NO_DESCRIPTOR;
	//Image slot 0 is kept for the atlas
	if (!tfx_impl_vulkan__create_index_list(&renderer->free_buffer_indexes, TFX_IMPL_VULKAN_BUFFER_DESCRIPTOR_CAPACITY, 0)
		|| !tfx_impl_vulkan__create_index_list(&renderer->free_image_indexes, TFX_IMPL_VULKAN_IMAGE_DESCRIPTOR_CAPACITY, 1)
		|| !tfx_impl_vulkan__create_descriptor_set(renderer)) {
		return false;
	}

	VkPushConstantRange push_range = { TFX_IMPL_VULKAN_SHADER_STAGES, 0, TFX_IMPL_VULKAN_PUSH_CONSTANTS_SIZE };
	VkPipelineLayoutCreateInfo layout_info = { VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO };
	layout_info.setLayoutCount = 1;
	layout_info.pSetLayouts = &renderer->descriptor_set_layout;
	layout_info.pushConstantRangeCount = 1;
	layout_info.pPushConstantRanges = &push_range;
	if (vkCreatePipelineLayout(device, &layout_info, NULL, &renderer->pipeline_layout) != VK_SUCCESS) {
		tfx_Print("Could not create the pipeline layout");
		return false;
	}

	//The instances and ribbon vertexes are vertex buffers in the bindless model, in the formats the Shaders readme lists
	const VkVertexInputBindingDescription billboard_binding = { 0, sizeof(tfx_instance_t), VK_VERTEX_INPUT_RATE_INSTANCE };
	const VkVertexInputAttributeDescription billboard_attributes[] = {
		{ 0, 0, VK_FORMAT_R32G32B32A32_SFLOAT, offsetof(tfx_instance_t, position) },
		{ 1, 0, VK_FORMAT_R16G16B16A16_SNORM, offsetof(tfx_instance_t, quaternion) },
		{ 2, 0, VK_FORMAT_R16G16_SINT, offsetof(tfx_instance_t, size) },
		{ 3, 0, VK_FORMAT_R8G8B8A8_SNORM, offsetof(tfx_instance_t, alignment) },
		{ 4, 0, VK_FORMAT_R16G16_SINT, offsetof(tfx_instance_t, intensity_gradient_map) },
		{ 5, 0, VK_FORMAT_R8G8B8A8_UNORM, offsetof(tfx_instance_t, curved_alpha_life) },
		{ 6, 0, VK_FORMAT_R32_UINT, offsetof(tfx_instance_t, indexes) },
		{ 7, 0, VK_FORMAT_R32_UINT, offsetof(tfx_instance_t, captured_index) },
	};
	const VkVertexInputBindingDescription ribbon_binding = { 0, sizeof(tfx_ribbon_vertex_t), VK_VERTEX_INPUT_RATE_VERTEX };
	const VkVertexInputAttributeDescription ribbon_attributes[] = {
		{ 0, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(tfx_ribbon_vertex_t, position) },
		{ 1, 0, VK_FORMAT_R32_UINT, offsetof(tfx_ribbon_vertex_t, segment_index) },
		{ 2, 0, VK_FORMAT_R32G32_SFLOAT, offsetof(tfx_ribbon_vertex_t, uv_offset_scale) },
		{ 3, 0, VK_FORMAT_R32_UINT, offsetof(tfx_ribbon_vertex_t, ribbon_index) },
		{ 4, 0, VK_FORMAT_R32_UINT, offsetof(tfx_ribbon_vertex_t, normal) },			//Only the lit vertex shader reads it
	};
	for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
		for (int lit = 0; lit != 2; ++lit) {
			renderer->billboard_pipelines[writes_depth][lit] = tfx_impl_vulkan__create_graphics_pipeline(renderer, shader_folder, tfx_impl_billboard_vertex_entry_points[lit],
				tfx_impl_billboard_fragment_entry_points[writes_depth][lit], &billboard_binding, billboard_attributes, 8, color_format, depth_format, writes_depth != 0);
			renderer->ribbon_pipelines[writes_depth][lit] = tfx_impl_vulkan__create_graphics_pipeline(renderer, shader_folder, tfx_impl_ribbon_vertex_entry_points[lit],
				tfx_impl_ribbon_fragment_entry_points[writes_depth][lit], &ribbon_binding, ribbon_attributes, lit ? 5 : 4, color_format, depth_format, writes_depth != 0);
			if (!renderer->billboard_pipelines[writes_depth][lit] || !renderer->ribbon_pipelines[writes_depth][lit]) {
				return false;
			}
		}
	}
	//Baked sprites are already interpolated, so the baked vertex shader has no captured_index input
	renderer->baked_billboard_pipeline = tfx_impl_vulkan__create_graphics_pipeline(renderer, shader_folder, "tfx_billboard_baked_vertex", "tfx_billboard_fragment",
		&billboard_binding, billboard_attributes, 7, color_format, depth_format, false);
	renderer->particle_kernel = tfx_impl_vulkan__create_compute_pipeline(renderer, shader_folder, "simulate_particles");
	renderer->ribbon_compute_pipeline = tfx_impl_vulkan__create_compute_pipeline(renderer, shader_folder, "tfx_ribbon_compute");
	renderer->sprite_playback_pipeline = tfx_impl_vulkan__create_compute_pipeline(renderer, shader_folder, "tfx_sprite_playback");
	renderer->ribbon_playback_pipeline = tfx_impl_vulkan__create_compute_pipeline(renderer, shader_folder, "tfx_ribbon_playback");
	if (!renderer->baked_billboard_pipeline || !renderer->particle_kernel
		|| !renderer->ribbon_compute_pipeline || !renderer->sprite_playback_pipeline || !renderer->ribbon_playback_pipeline) {
		return false;
	}

	//Each frame slot's view uniforms stay mapped, so the draws read what tfx_impl_vulkan_Render wrote without a copy
	for (tfxU32 frame = 0; frame != renderer->frames_in_flight; ++frame) {
		if (!tfx_impl_vulkan__create_mapped_buffer(renderer, TFX_IMPL_VULKAN_VIEW_UNIFORMS_SIZE, VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT, &renderer->view_uniforms[frame])) {
			return false;
		}
		tfx_impl_vulkan__write_buffer_descriptor(renderer, TFX_IMPL_VULKAN_UNIFORM_BINDING, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, frame, renderer->view_uniforms[frame].buffer);
	}
	return true;
}

void tfx_impl_vulkan_BeginFrame(tfx_impl_vulkan_t *renderer, tfxU32 frame_index) {
	renderer->frame_index = frame_index % renderer->frames_in_flight;
	tfx_impl_vulkan__process_release_list(renderer, &renderer->release_lists[renderer->frame_index]);
	renderer->upload_rings[renderer->frame_index].used = 0;
}

bool tfx_impl_vulkan_CreateAtlas(tfx_impl_vulkan_t *renderer, tfxU32 layer_size, tfxU32 layer_count) {
	tfx_impl_vulkan__destroy_image(renderer, &renderer->shape_atlas);
	tfxU32 mip_count = 1;
	while ((layer_size >> mip_count) > 0) {
		mip_count++;
	}
	//A blit source too, for its mips
	if (!tfx_impl_vulkan__create_image(renderer, layer_size, layer_size, layer_count, mip_count, VK_FORMAT_R8G8B8A8_UNORM, VK_IMAGE_USAGE_TRANSFER_SRC_BIT, &renderer->shape_atlas)) {
		return false;
	}
	//Slot 0, because tfx_impl_GetUV leaves the image index in the upper 16 bits of texture_array_index at 0
	renderer->shape_atlas.index = 0;
	tfx_impl_vulkan__write_image_descriptor(renderer, &renderer->shape_atlas);
	return true;
}

bool tfx_impl_vulkan_UploadAtlas(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_atlas_t *atlas) {
	tfx_impl_vulkan_image_t *texture = &renderer->shape_atlas;
	if (!texture->image || atlas->layer_size != texture->width || atlas->first_layer + atlas->layer_count > texture->layer_count) {
		tfx_Print("An atlas of %u layers from layer %u, %u pixels square, doesn't fit the %u layer, %u pixel texture array",
			atlas->layer_count, atlas->first_layer, atlas->layer_size, texture->layer_count, texture->width);
		return false;
	}
	//The view covers layers that no atlas has filled yet, so the whole array starts out in the layout the draws read it in
	if (texture->layout == VK_IMAGE_LAYOUT_UNDEFINED) {
		tfx_impl_vulkan__image_barrier(command_buffer, texture->image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, 0, texture->mip_count, 0, texture->layer_count);
		texture->layout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
	}
	//Everything in these layers is replaced, so what they held can be discarded
	tfx_impl_vulkan__image_barrier(command_buffer, texture->image, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 0, texture->mip_count, atlas->first_layer, atlas->layer_count);
	bool uploaded = tfx_impl_vulkan__upload_image_layers(renderer, command_buffer, texture, atlas->first_layer, atlas->layer_count, 4, atlas->pixels);
	//Each mip is a linear blit of the one above, a 2 x 2 box filter
	tfxU32 size = texture->width;
	for (tfxU32 mip_level = 1; mip_level < texture->mip_count && uploaded; ++mip_level) {
		tfx_impl_vulkan__image_barrier(command_buffer, texture->image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL, mip_level - 1, 1, atlas->first_layer, atlas->layer_count);
		tfxU32 next_size = size > 1 ? size / 2 : 1;
		VkImageBlit blit = {};
		blit.srcSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
		blit.srcSubresource.mipLevel = mip_level - 1;
		blit.srcSubresource.baseArrayLayer = atlas->first_layer;
		blit.srcSubresource.layerCount = atlas->layer_count;
		blit.srcOffsets[1].x = (int32_t)size;
		blit.srcOffsets[1].y = (int32_t)size;
		blit.srcOffsets[1].z = 1;
		blit.dstSubresource = blit.srcSubresource;
		blit.dstSubresource.mipLevel = mip_level;
		blit.dstOffsets[1].x = (int32_t)next_size;
		blit.dstOffsets[1].y = (int32_t)next_size;
		blit.dstOffsets[1].z = 1;
		vkCmdBlitImage(command_buffer, texture->image, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL, texture->image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 1, &blit, VK_FILTER_LINEAR);
		size = next_size;
	}
	//Every mip but the last was a blit source
	tfxU32 source_mips = uploaded ? texture->mip_count - 1 : 0;
	if (source_mips) {
		tfx_impl_vulkan__image_barrier(command_buffer, texture->image, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, 0, source_mips, atlas->first_layer, atlas->layer_count);
	}
	tfx_impl_vulkan__image_barrier(command_buffer, texture->image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, source_mips, texture->mip_count - source_mips,
		atlas->first_layer, atlas->layer_count);
	return uploaded;
}

bool tfx_impl_vulkan_UploadSpriteData(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, tfx_animation_manager animation_manager, tfx_gpu_shapes shapes) {
	tfx_impl_vulkan__begin_transfers(command_buffer);
	tfxU32 bitmap_count = tfx_GetAnimationColorRampBitmapCount(animation_manager);
	if (bitmap_count) {
		tfx_impl_vulkan_image_t *texture = &renderer->baked_color_ramps;
		if (!tfx_impl_vulkan__replace_color_ramps(renderer, texture, tfx_GetAnimationColorRampBitmap(animation_manager, 0), bitmap_count)) {
			return false;
		}
		tfx_impl_vulkan__begin_image_upload(command_buffer, texture);
		bool uploaded = true;
		for (tfxU32 index = 0; index != bitmap_count && uploaded; ++index) {
			uploaded = tfx_impl_vulkan__upload_color_ramp_layer(renderer, command_buffer, texture, index, tfx_GetAnimationColorRampBitmap(animation_manager, index));
		}
		tfx_impl_vulkan__end_image_upload(command_buffer, texture);
		if (!uploaded) {
			return false;
		}
	}

	tfxU32 shape_count = tfx_GetGPUShapesCount(shapes);
	if (!tfx_impl_vulkan__write_buffer(renderer, command_buffer, &renderer->baked_image_data, tfx_GetGPUShapesArray(shapes), sizeof(tfx_gpu_image_data_t) * shape_count)
		|| !tfx_impl_vulkan__write_buffer(renderer, command_buffer, &renderer->baked_emitter_properties, tfx_GetAnimationEmitterPropertiesBufferPointer(animation_manager),
			tfx_GetAnimationEmitterPropertySizeInBytes(animation_manager))
		|| !tfx_impl_vulkan__write_buffer(renderer, command_buffer, &renderer->baked_sprites, tfx_GetSpriteDataBufferPointer(animation_manager), tfx_GetSpriteDataSizeInBytes(animation_manager))) {
		return false;
	}
	if (!tfx_AnimationManagerHasRibbons(animation_manager)) {
		return true;
	}

	//The baked ribbons, their properties, segments, emitters and graph lookups only change when sprite data is loaded
	size_t emitters_bytes = tfx_GetAnimationRibbonEmittersSizeInBytes(animation_manager);
	void *emitters = malloc(emitters_bytes);
	if (!emitters) {
		return false;
	}
	tfx_CopyAnimationRibbonEmitters(animation_manager, emitters);
	bool uploaded = tfx_impl_vulkan__write_buffer(renderer, command_buffer, &renderer->baked_ribbons.graph_lookups, tfx_GetAnimationRibbonGraphLookupsBufferPointer(animation_manager),
			tfx_GetAnimationRibbonGraphLookupsSizeInBytes(animation_manager))
		&& tfx_impl_vulkan__write_buffer(renderer, command_buffer, &renderer->baked_ribbon_data, tfx_GetAnimationRibbonDataBufferPointer(animation_manager),
			tfx_GetAnimationRibbonDataSizeInBytes(animation_manager))
		&& tfx_impl_vulkan__write_buffer(renderer, command_buffer, &renderer->baked_ribbon_properties, tfx_GetAnimationRibbonPropertiesBufferPointer(animation_manager),
			tfx_GetAnimationRibbonPropertySizeInBytes(animation_manager))
		&& tfx_impl_vulkan__write_buffer(renderer, command_buffer, &renderer->baked_ribbons.segments, tfx_GetAnimationRibbonSegmentsBufferPointer(animation_manager),
			tfx_GetAnimationRibbonSegmentsSizeInBytes(animation_manager))
		&& tfx_impl_vulkan__write_buffer(renderer, command_buffer, &renderer->baked_ribbons.emitters, emitters, emitters_bytes);
	free(emitters);
	return uploaded;
}

//The vertex and index buffers ribbon_compute writes, sized for this frame's buckets
static bool tfx_impl_vulkan__ensure_ribbon_geometry(tfx_impl_vulkan_t *renderer, const tfx_ribbon_batch_sizes_t *sizes, tfx_impl_vulkan_ribbon_set_t *set) {
	return tfx_impl_vulkan__ensure_buffer(renderer, sizes->vertex_buffer_size_in_bytes, VK_BUFFER_USAGE_VERTEX_BUFFER_BIT, &set->vertices)
		&& tfx_impl_vulkan__ensure_buffer(renderer, sizes->index_buffer_size_in_bytes, VK_BUFFER_USAGE_INDEX_BUFFER_BIT, &set->indices);
}

static bool tfx_impl_vulkan__upload_baked_ribbons(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_baked_frame_t *baked) {
	if (baked->ribbon_count == 0 || !renderer->baked_ribbon_data.buffer) {
		return true;
	}
	tfx_ribbon_batch_sizes_t sizes = tfx_GetAnimationRibbonBufferSizes(baked->animation_manager, sizeof(tfx_ribbon_vertex_t));
	if (!tfx_impl_vulkan__write_buffer(renderer, command_buffer, &renderer->animation_ribbon_offsets, baked->ribbon_offsets, baked->ribbon_offsets_bytes)
		|| !tfx_impl_vulkan__ensure_buffer(renderer, (size_t)baked->ribbon_count * tfxRIBBON_SIZE, 0, &renderer->baked_ribbons.instances)
		|| !tfx_impl_vulkan__ensure_ribbon_geometry(renderer, &sizes, &renderer->baked_ribbons)) {
		return false;
	}
	tfx_ribbon_data_push_t *push = &renderer->ribbon_data_push;
	memset(push, 0, sizeof(tfx_ribbon_data_push_t));
	push->animation_instances_total = baked->animation_instance_count;
	push->ribbons_total = baked->ribbon_count;
	push->offsets_index = renderer->animation_ribbon_offsets.index;
	push->animation_instances_index = renderer->animation_instances.index;
	push->ribbon_data_index = renderer->baked_ribbon_data.index;
	push->ribbon_output_index = renderer->baked_ribbons.instances.index;
	push->ribbon_properties_index = renderer->baked_ribbon_properties.index;
	renderer->playback_ribbon_count = baked->ribbon_count;
	renderer->baked_ribbons.animation_manager = baked->animation_manager;
	return true;
}

//Animation instances advance every frame, so this runs every frame, not just on updates
static bool tfx_impl_vulkan__upload_baked(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_baked_frame_t *baked) {
	renderer->playback_sprite_count = 0;
	renderer->playback_ribbon_count = 0;
	renderer->baked_ribbons.animation_manager = NULL;
	if (!renderer->baked_sprites.buffer || (baked->sprite_count == 0 && baked->ribbon_count == 0)) {
		return true;
	}
	if (!tfx_impl_vulkan__write_buffer(renderer, command_buffer, &renderer->animation_instances, baked->animation_instances, baked->animation_instances_bytes)) {
		return false;
	}
	if (baked->sprite_count > 0) {
		if (!tfx_impl_vulkan__write_buffer(renderer, command_buffer, &renderer->animation_offsets, baked->offsets, baked->offsets_bytes)
			|| !tfx_impl_vulkan__ensure_buffer(renderer, (size_t)baked->sprite_count * sizeof(tfx_instance_t), VK_BUFFER_USAGE_VERTEX_BUFFER_BIT, &renderer->playback_instances)) {
			return false;
		}
		tfx_sprite_data_push_t *push = &renderer->sprite_data_push;
		memset(push, 0, sizeof(tfx_sprite_data_push_t));
		push->animation_instances_total = baked->animation_instance_count;
		push->billboards_total = baked->sprite_count;
		push->animated_shapes = baked->animated_shapes ? 1 : 0;
		push->offsets_index = renderer->animation_offsets.index;
		push->animation_instances_index = renderer->animation_instances.index;
		push->billboards_index = renderer->playback_instances.index;
		push->sprite_data_index = renderer->baked_sprites.index;
		push->image_data_index = renderer->baked_image_data.index;
		push->emitter_properties_index = renderer->baked_emitter_properties.index;
		renderer->playback_sprite_count = baked->sprite_count;
	}
	return tfx_impl_vulkan__upload_baked_ribbons(renderer, command_buffer, baked);
}

static bool tfx_impl_vulkan__upload_ribbons(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_ribbon_frame_t *ribbons) {
	tfx_impl_vulkan_ribbon_set_t *set = &renderer->live_ribbons;
	set->stage = ribbons->stage;
	if (!ribbons->stage) {
		return true;
	}
	const tfx_ribbon_batch_sizes_t *sizes = &ribbons->sizes;
	if (!tfx_impl_vulkan__write_buffer(renderer, command_buffer, &set->segments, ribbons->segments, sizes->segment_buffer_size_in_bytes)
		|| !tfx_impl_vulkan__write_buffer(renderer, command_buffer, &set->instances, ribbons->ribbons, sizes->ribbon_buffer_size_in_bytes)
		|| !tfx_impl_vulkan__write_buffer(renderer, command_buffer, &set->emitters, ribbons->emitters, sizes->emitter_buffer_size_in_bytes)
		|| !tfx_impl_vulkan__ensure_ribbon_geometry(renderer, sizes, set)) {
		set->stage = NULL;
		return false;
	}
	set->stage_offsets = ribbons->stage_offsets;
	return true;
}

//The old contents are copied across, because the previous buffer is still what captured_index points into
static bool tfx_impl_vulkan__grow_instance_buffers(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, tfxU32 needed) {
	tfxU32 capacity = renderer->instance_capacity ? renderer->instance_capacity : TFX_IMPL_VULKAN_INITIAL_INSTANCE_CAPACITY;
	while (capacity < needed) {
		capacity += capacity / 2;
	}
	for (int index = 0; index != 2; ++index) {
		tfx_impl_vulkan_buffer_t grown = {};
		//The GPU particle kernel writes its instances after the CPU ones
		if (!tfx_impl_vulkan__create_buffer(renderer, (VkDeviceSize)capacity * TFX_IMPL_VULKAN_INSTANCE_STRIDE, VK_BUFFER_USAGE_VERTEX_BUFFER_BIT, &grown)) {
			return false;
		}
		tfx_impl_vulkan_buffer_t *old = &renderer->instance_buffers[index];
		if (old->buffer) {
			VkBufferCopy region = { 0, 0, old->capacity };
			vkCmdCopyBuffer(command_buffer, old->buffer, grown.buffer, 1, &region);
		}
		tfx_impl_vulkan__release_buffer(renderer, old);
		*old = grown;
	}
	renderer->instance_capacity = capacity;
	//The instance upload that follows writes into one of the copies
	tfx_impl_vulkan__memory_barrier(command_buffer, VK_PIPELINE_STAGE_2_ALL_TRANSFER_BIT, VK_ACCESS_2_TRANSFER_WRITE_BIT, VK_PIPELINE_STAGE_2_ALL_TRANSFER_BIT, VK_ACCESS_2_TRANSFER_WRITE_BIT);
	return true;
}

static bool tfx_impl_vulkan__upload_gpu_particles(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_frame_t *frame) {
	const tfx_gpu_particle_frame_t *gpu_frame = &frame->gpu_particles;
	const tfx_gpu_particle_buffer_sizes_t *pool_sizes = &frame->gpu_pool_sizes;
	//The pool holds the particles between frames, it's sized once and starts zeroed
	if (!renderer->particle_pool.buffer) {
		if (!tfx_impl_vulkan__create_buffer(renderer, pool_sizes->pool_bytes, 0, &renderer->particle_pool)) {
			return false;
		}
		vkCmdFillBuffer(command_buffer, renderer->particle_pool.buffer, 0, VK_WHOLE_SIZE, 0);
	}
	if (!tfx_impl_vulkan__ensure_buffer(renderer, gpu_frame->upload_bytes, 0, &renderer->particle_upload)) {
		return false;
	}
	if (gpu_frame->upload_bytes) {
		VkBuffer source = VK_NULL_HANDLE;
		VkDeviceSize offset = 0;
		unsigned char *mapped = tfx_impl_vulkan__allocate_upload(renderer, gpu_frame->upload_bytes, 16, &source, &offset);
		if (!mapped) {
			return false;
		}
		tfx_CopyGPUParticleUpload(gpu_frame, mapped);
		VkBufferCopy region = { offset, 0, gpu_frame->upload_bytes };
		vkCmdCopyBuffer(command_buffer, source, renderer->particle_upload.buffer, 1, &region);
	}

	//The pool's sections share its index, and the upload's arrays share theirs. The bases are in rows of each section's own stride
	tfx_gpu_particle_push_t *push = &renderer->particle_push_constants;
	memset(push, 0, sizeof(tfx_gpu_particle_push_t));
	push->hot_buffer_index = renderer->particle_pool.index;
	push->cold_buffer_index = renderer->particle_pool.index;
	push->noise_buffer_index = renderer->particle_pool.index;
	push->graph_buffer_index = renderer->graph_lookups.index;
	push->dispatch_buffer_index = renderer->particle_upload.index;
	push->tick_buffer_index = renderer->particle_upload.index;
	push->spawn_buffer_index = renderer->particle_upload.index;
	push->location_buffer_index = renderer->particle_upload.index;
	push->cold_base = pool_sizes->cold_base;
	push->noise_base = pool_sizes->noise_base;
	push->tick_base = gpu_frame->tick_base;
	push->spawn_base = gpu_frame->spawn_base;
	push->location_base = gpu_frame->location_base;
	renderer->particle_dispatches = gpu_frame->dispatches;
	renderer->particle_dispatch_count = gpu_frame->dispatch_count;
	return true;
}

//Global layers never move, so only the layers whose version changed go up, unless the texture array has to grow
static bool tfx_impl_vulkan__upload_color_ramps(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer) {
	tfxU32 layer_count = tfx_GetGlobalColorRampLayerCount();
	tfx_bitmap_t *first_bitmap = NULL;
	for (tfxU32 layer = 0; layer != layer_count && !first_bitmap; ++layer) {
		first_bitmap = tfx_GetGlobalColorRampLayer(layer);
	}
	if (!first_bitmap) {
		return true;
	}
	tfx_impl_vulkan_image_t *texture = &renderer->color_ramps;
	bool recreate = layer_count > texture->layer_count;
	if (recreate) {
		//Headroom so that loading another library doesn't always recreate it
		tfxU32 capacity = layer_count + layer_count / 2 + 1;
		capacity = capacity < TFX_IMPL_VULKAN_MAX_COLOR_RAMP_LAYERS ? capacity : TFX_IMPL_VULKAN_MAX_COLOR_RAMP_LAYERS;
		tfxU32 *layer_versions = (tfxU32 *)realloc(renderer->color_ramp_layer_versions, sizeof(tfxU32) * capacity);
		if (!layer_versions) {
			return false;
		}
		renderer->color_ramp_layer_versions = layer_versions;
		if (!tfx_impl_vulkan__replace_color_ramps(renderer, texture, first_bitmap, capacity)) {
			return false;
		}
	}
	tfx_impl_vulkan__begin_image_upload(command_buffer, texture);
	bool uploaded = true;
	for (tfxU32 layer = 0; layer != layer_count && uploaded; ++layer) {
		tfxU32 version = tfx_GetGlobalColorRampLayerVersion(layer);
		if (!recreate && version == renderer->color_ramp_layer_versions[layer]) {
			continue;
		}
		renderer->color_ramp_layer_versions[layer] = version;
		//A free layer isn't drawn from, so whatever it holds can stay
		tfx_bitmap_t *bitmap = tfx_GetGlobalColorRampLayer(layer);
		if (bitmap) {
			uploaded = tfx_impl_vulkan__upload_color_ramp_layer(renderer, command_buffer, texture, layer, bitmap);
		}
	}
	tfx_impl_vulkan__end_image_upload(command_buffer, texture);
	return uploaded;
}

//Only on update frames: the instances drawn until then were written with the layout the last upload holds
static bool tfx_impl_vulkan__upload_library_data(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer) {
	tfxU32 properties_version = tfx_GetGlobalParticlePropertiesVersion();
	if (!renderer->emitter_properties.buffer || properties_version != renderer->properties_version) {
		if (!tfx_impl_vulkan__write_buffer(renderer, command_buffer, &renderer->emitter_properties, tfx_GetGlobalParticlePropertiesBuffer(), tfx_GetGlobalParticlePropertiesBufferSizeInBytes())) {
			return false;
		}
		renderer->properties_version = properties_version;
	}
	tfxU32 shapes_version = tfx_GetGlobalGPUShapesVersion();
	if (!renderer->image_data.buffer || shapes_version != renderer->shapes_version) {
		tfx_gpu_shapes shapes = tfx_GetGlobalGPUShapes();
		if (!tfx_impl_vulkan__write_buffer(renderer, command_buffer, &renderer->image_data, tfx_GetGPUShapesArray(shapes), sizeof(tfx_gpu_image_data_t) * tfx_GetGPUShapesCount(shapes))) {
			return false;
		}
		renderer->shapes_version = shapes_version;
	}
	//Shared by every library already, read by the ribbon shaders and the GPU kernel
	tfxU32 graph_lookups_version = tfx_GetGPUGraphLookupsVersion();
	if (!renderer->graph_lookups.buffer || graph_lookups_version != renderer->graph_lookups_version) {
		if (!tfx_impl_vulkan__write_buffer(renderer, command_buffer, &renderer->graph_lookups, tfx_GetGPUGraphLookupsBuffer(), tfx_GetGPUGraphLookupsBufferSizeInBytes())) {
			return false;
		}
		renderer->graph_lookups_version = graph_lookups_version;
	}
	tfxU32 color_ramps_version = tfx_GetGlobalColorRampsVersion();
	if (!renderer->color_ramps.image || color_ramps_version != renderer->color_ramps_version) {
		if (!tfx_impl_vulkan__upload_color_ramps(renderer, command_buffer)) {
			return false;
		}
		renderer->color_ramps_version = color_ramps_version;
	}
	return true;
}

bool tfx_impl_vulkan_UploadFrame(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_frame_t *frame) {
	memcpy(renderer->layer_flags, frame->layer_flags, sizeof(renderer->layer_flags));
	tfx_impl_vulkan__begin_transfers(command_buffer);
	if (!tfx_impl_vulkan__upload_baked(renderer, command_buffer, &frame->baked)) {
		return false;
	}
	if (!frame->stage_updated) {
		return true;
	}
	if (!tfx_impl_vulkan__upload_library_data(renderer, command_buffer)) {
		return false;
	}
	bool gpu_particles = frame->has_gpu_particles && renderer->particle_kernel;
	tfxU32 needed = gpu_particles ? frame->gpu_particles.gpu_instance_start + frame->gpu_particles.gpu_instance_count : frame->instance_count;
	if ((needed > renderer->instance_capacity || renderer->instance_capacity == 0) && !tfx_impl_vulkan__grow_instance_buffers(renderer, command_buffer, needed)) {
		renderer->instance_count = 0;
		return false;
	}

	//The other buffer now holds what was drawn last, which is what captured_index points into for CPU and GPU instances alike
	renderer->current_instance_buffer ^= 1;
	if (!tfx_impl_vulkan__upload_buffer(renderer, command_buffer, &renderer->instance_buffers[renderer->current_instance_buffer], frame->instances, frame->instance_count * sizeof(tfx_instance_t))) {
		return false;
	}
	renderer->instance_count = needed;
	memcpy(renderer->layer_start, frame->layer_start, sizeof(renderer->layer_start));
	memcpy(renderer->layer_count, frame->layer_count, sizeof(renderer->layer_count));
	memset(renderer->gpu_layer_count, 0, sizeof(renderer->gpu_layer_count));
	if (gpu_particles) {
		memcpy(renderer->gpu_layer_start, frame->gpu_particles.gpu_layer_start, sizeof(renderer->gpu_layer_start));
		memcpy(renderer->gpu_layer_count, frame->gpu_particles.gpu_layer_count, sizeof(renderer->gpu_layer_count));
		if (!tfx_impl_vulkan__upload_gpu_particles(renderer, command_buffer, frame)) {
			return false;
		}
	}
	return tfx_impl_vulkan__upload_ribbons(renderer, command_buffer, &frame->ribbons);
}

static void tfx_impl_vulkan__push(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const void *constants, tfxU32 size_in_bytes) {
	vkCmdPushConstants(command_buffer, renderer->pipeline_layout, TFX_IMPL_VULKAN_SHADER_STAGES, 0, size_in_bytes, constants);
}

//Every fetch has to be dispatched, the next one carries on from the particle state this one leaves in the pool. The dispatches
//write separate rings and instance ranges, so they need no barrier between them
static void tfx_impl_vulkan__simulate_gpu_particles(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer) {
	vkCmdBindPipeline(command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE, renderer->particle_kernel);
	tfx_gpu_particle_push_t *push = &renderer->particle_push_constants;
	push->instance_buffer_index = renderer->instance_buffers[renderer->current_instance_buffer].index;
	tfx_impl_vulkan__push(renderer, command_buffer, push, sizeof(tfx_gpu_particle_push_t));
	for (tfxU32 index = 0; index != renderer->particle_dispatch_count; ++index) {
		tfxU32 slot_count = renderer->particle_dispatches[index].slot_count;
		if (slot_count == 0) {
			continue;
		}
		//dispatch_index is the first push constant, so it's the only one that changes
		tfx_impl_vulkan__push(renderer, command_buffer, &index, sizeof(tfxU32));
		vkCmdDispatch(command_buffer, (slot_count + TFX_IMPL_VULKAN_KERNEL_GROUP_SIZE - 1) / TFX_IMPL_VULKAN_KERNEL_GROUP_SIZE, 1, 1);
	}
	renderer->particle_dispatch_count = 0;
}

static bool tfx_impl_vulkan__ribbon_set_active(const tfx_impl_vulkan_ribbon_set_t *set) {
	return set->stage || set->animation_manager;
}

static tfx_ribbon_dispatch_t tfx_impl_vulkan__create_ribbon_dispatch(const tfx_impl_vulkan_ribbon_set_t *set) {
	return set->stage ? tfx_CreateRibbonBatchDispatch(&set->stage_offsets) : tfx_CreateAnimationRibbonDispatch();
}

static bool tfx_impl_vulkan__next_ribbon_dispatch(const tfx_impl_vulkan_ribbon_set_t *set, tfx_ribbon_dispatch_t *dispatch) {
	return set->stage ? tfx_NextRibbonDispatch(set->stage, dispatch) : tfx_NextAnimationRibbonDispatch(set->animation_manager, dispatch);
}

//The library fills a bucket's counts and offsets, the descriptor indexes are the host's
static void tfx_impl_vulkan__set_ribbon_indexes(tfx_impl_vulkan_t *renderer, tfx_impl_vulkan_ribbon_set_t *set, tfx_ribbon_bucket_globals_t *globals) {
	globals->uniform_index = renderer->frame_index;
	globals->emitters_index = set->emitters.index;
	globals->graphs_index = set->stage ? renderer->graph_lookups.index : set->graph_lookups.index;
	globals->ribbons_index = set->instances.index;
	globals->ribbon_segments_index = set->segments.index;
	globals->vertexes_index = set->vertices.index;
	globals->indexes_index = set->indices.index;
	globals->image_data_index = set->stage ? renderer->image_data.index : renderer->baked_image_data.index;
	globals->sampler_index = 0;
	globals->particle_texture_index = renderer->shape_atlas.index;
	globals->color_ramp_texture_index = set->stage ? renderer->color_ramps.index : renderer->baked_color_ramps.index;
}

static void tfx_impl_vulkan__build_ribbon_geometry(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, tfx_impl_vulkan_ribbon_set_t *set, const tfx_impl_uniforms_t *uniforms) {
	vkCmdBindPipeline(command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE, renderer->ribbon_compute_pipeline);
	//One dispatch per bucket, each with the bucket's own globals
	tfx_ribbon_dispatch_t dispatch = tfx_impl_vulkan__create_ribbon_dispatch(set);
	while (tfx_impl_vulkan__next_ribbon_dispatch(set, &dispatch)) {
		tfx_ribbon_bucket_globals_t *globals = tfx_GetRibbonDispatchGlobals(&dispatch);
		if (set->stage) {
			globals->lerp = uniforms->timer_lerp;
		} else {
			//Baked ribbons are already interpolated, but they still face the camera
			globals->camera_position.x = uniforms->camera_position[0];
			globals->camera_position.y = uniforms->camera_position[1];
			globals->camera_position.z = uniforms->camera_position[2];
		}
		tfx_impl_vulkan__set_ribbon_indexes(renderer, set, globals);
		tfx_impl_vulkan__push(renderer, command_buffer, globals, sizeof(tfx_ribbon_bucket_globals_t));
		vkCmdDispatch(command_buffer, (dispatch.total_segments + TFX_IMPL_VULKAN_RIBBON_GROUP_SIZE - 1) / TFX_IMPL_VULKAN_RIBBON_GROUP_SIZE, 1, 1);
	}
}

void tfx_impl_vulkan_Compute(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_uniforms_t *uniforms) {
	//This frame's uploads, and the last frame's reads of anything the passes below write again, finish first
	tfx_impl_vulkan__memory_barrier(command_buffer, VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, VK_ACCESS_2_MEMORY_WRITE_BIT,
		VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT, VK_ACCESS_2_SHADER_STORAGE_READ_BIT | VK_ACCESS_2_SHADER_STORAGE_WRITE_BIT);
	vkCmdBindDescriptorSets(command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE, renderer->pipeline_layout, 0, 1, &renderer->descriptor_set, 0, NULL);
	if (renderer->particle_dispatch_count > 0) {
		tfx_impl_vulkan__simulate_gpu_particles(renderer, command_buffer);
	}
	if (renderer->playback_sprite_count > 0) {
		vkCmdBindPipeline(command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE, renderer->sprite_playback_pipeline);
		tfx_impl_vulkan__push(renderer, command_buffer, &renderer->sprite_data_push, sizeof(tfx_sprite_data_push_t));
		vkCmdDispatch(command_buffer, (renderer->playback_sprite_count + TFX_IMPL_VULKAN_PLAYBACK_GROUP_SIZE - 1) / TFX_IMPL_VULKAN_PLAYBACK_GROUP_SIZE, 1, 1);
	}
	if (renderer->playback_ribbon_count > 0) {
		vkCmdBindPipeline(command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE, renderer->ribbon_playback_pipeline);
		tfx_impl_vulkan__push(renderer, command_buffer, &renderer->ribbon_data_push, sizeof(tfx_ribbon_data_push_t));
		vkCmdDispatch(command_buffer, (renderer->playback_ribbon_count + TFX_IMPL_VULKAN_PLAYBACK_GROUP_SIZE - 1) / TFX_IMPL_VULKAN_PLAYBACK_GROUP_SIZE, 1, 1);
		//ribbon_compute reads the ribbons the playback just wrote
		tfx_impl_vulkan__memory_barrier(command_buffer, VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT, VK_ACCESS_2_SHADER_STORAGE_WRITE_BIT,
			VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT, VK_ACCESS_2_SHADER_STORAGE_READ_BIT);
	}
	if (tfx_impl_vulkan__ribbon_set_active(&renderer->live_ribbons)) {
		tfx_impl_vulkan__build_ribbon_geometry(renderer, command_buffer, &renderer->live_ribbons, uniforms);
	}
	if (tfx_impl_vulkan__ribbon_set_active(&renderer->baked_ribbons)) {
		tfx_impl_vulkan__build_ribbon_geometry(renderer, command_buffer, &renderer->baked_ribbons, uniforms);
	}
	//The draws read what the uploads and these passes wrote
	tfx_impl_vulkan__memory_barrier(command_buffer, VK_PIPELINE_STAGE_2_ALL_TRANSFER_BIT | VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT, VK_ACCESS_2_TRANSFER_WRITE_BIT | VK_ACCESS_2_SHADER_STORAGE_WRITE_BIT,
		VK_PIPELINE_STAGE_2_VERTEX_ATTRIBUTE_INPUT_BIT | VK_PIPELINE_STAGE_2_INDEX_INPUT_BIT | VK_PIPELINE_STAGE_2_VERTEX_SHADER_BIT | VK_PIPELINE_STAGE_2_FRAGMENT_SHADER_BIT,
		VK_ACCESS_2_VERTEX_ATTRIBUTE_READ_BIT | VK_ACCESS_2_INDEX_READ_BIT | VK_ACCESS_2_SHADER_STORAGE_READ_BIT);
}

//Live and baked ribbons each have their own shape frames and color ramps, the atlas is shared
//Every bucket is walked, because the dispatch carries the offsets from one to the next, and only the ones whose layer is in the group are drawn
static void tfx_impl_vulkan__render_ribbons(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, tfx_impl_vulkan_ribbon_set_t *set, int writes_depth) {
	VkDeviceSize offset = 0;
	vkCmdBindVertexBuffers(command_buffer, 0, 1, &set->vertices.buffer, &offset);
	vkCmdBindIndexBuffer(command_buffer, set->indices.buffer, 0, VK_INDEX_TYPE_UINT32);
	VkPipeline bound_pipeline = VK_NULL_HANDLE;
	tfx_ribbon_dispatch_t dispatch = tfx_impl_vulkan__create_ribbon_dispatch(set);
	while (tfx_impl_vulkan__next_ribbon_dispatch(set, &dispatch)) {
		const tfxLayerFlags flags = renderer->layer_flags[dispatch.layer];
		if (tfx_impl_LayerWritesDepth(flags) != writes_depth) {
			continue;
		}
		VkPipeline pipeline = renderer->ribbon_pipelines[writes_depth][tfx_impl_LayerIsLit(flags)];
		if (pipeline != bound_pipeline) {
			vkCmdBindPipeline(command_buffer, VK_PIPELINE_BIND_POINT_GRAPHICS, pipeline);
			bound_pipeline = pipeline;
		}
		tfx_ribbon_bucket_globals_t *globals = tfx_GetRibbonDispatchGlobals(&dispatch);
		tfx_impl_vulkan__set_ribbon_indexes(renderer, set, globals);
		tfx_impl_vulkan__push(renderer, command_buffer, globals, sizeof(tfx_ribbon_bucket_globals_t));
		vkCmdDrawIndexed(command_buffer, dispatch.index_count, 1, dispatch.index_offset, 0, 0);
	}
}

static void tfx_impl_vulkan__draw_instance_range(VkCommandBuffer command_buffer, tfxU32 start, tfxU32 count) {
	//The first instance starts the instance rate attributes at the range, so nothing extra goes in the push constants
	if (count) {
		vkCmdDraw(command_buffer, 6, count, 0, start);
	}
}

//Each layer in the group as its CPU instances then its GPU instances
static void tfx_impl_vulkan__render_billboards(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_vulkan_billboard_push_t *push, int writes_depth) {
	VkDeviceSize offset = 0;
	vkCmdBindVertexBuffers(command_buffer, 0, 1, &renderer->instance_buffers[renderer->current_instance_buffer].buffer, &offset);
	VkPipeline bound_pipeline = VK_NULL_HANDLE;
	for (tfxU32 layer = 0; layer != tfxLAYERS; ++layer) {
		const tfxLayerFlags flags = renderer->layer_flags[layer];
		if (tfx_impl_LayerWritesDepth(flags) != writes_depth || renderer->layer_count[layer] + renderer->gpu_layer_count[layer] == 0) {
			continue;
		}
		VkPipeline pipeline = renderer->billboard_pipelines[writes_depth][tfx_impl_LayerIsLit(flags)];
		if (pipeline != bound_pipeline) {
			vkCmdBindPipeline(command_buffer, VK_PIPELINE_BIND_POINT_GRAPHICS, pipeline);
			tfx_impl_vulkan__push(renderer, command_buffer, push, sizeof(tfx_impl_vulkan_billboard_push_t));
			bound_pipeline = pipeline;
		}
		tfx_impl_vulkan__draw_instance_range(command_buffer, renderer->layer_start[layer], renderer->layer_count[layer]);
		tfx_impl_vulkan__draw_instance_range(command_buffer, renderer->gpu_layer_start[layer], renderer->gpu_layer_count[layer]);
	}
}

static void tfx_impl_vulkan__render_baked_billboards(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_vulkan_billboard_push_t *push) {
	//Baked sprites come with their own shape frames, emitter properties and color ramps, the atlas is shared. The baked push starts
	//with the live push's fields, which the shared fragment shader reads
	tfx_impl_vulkan_baked_billboard_push_t baked_push = {};
	baked_push.billboard = *push;
	baked_push.billboard.color_ramp_texture_index = renderer->baked_color_ramps.index;
	baked_push.billboard.image_data_index = renderer->baked_image_data.index;
	baked_push.billboard.properties_index = renderer->baked_emitter_properties.index;
	baked_push.billboard.prev_billboards_index = 0;
	vkCmdBindPipeline(command_buffer, VK_PIPELINE_BIND_POINT_GRAPHICS, renderer->baked_billboard_pipeline);
	tfx_impl_vulkan__push(renderer, command_buffer, &baked_push, sizeof(baked_push));
	VkDeviceSize offset = 0;
	vkCmdBindVertexBuffers(command_buffer, 0, 1, &renderer->playback_instances.buffer, &offset);
	tfx_impl_vulkan__draw_instance_range(command_buffer, 0, renderer->playback_sprite_count);
}

void tfx_impl_vulkan_Render(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_uniforms_t *uniforms) {
	if (!renderer->shape_atlas.image || (renderer->instance_count == 0 && renderer->playback_sprite_count == 0
		&& !tfx_impl_vulkan__ribbon_set_active(&renderer->live_ribbons) && !tfx_impl_vulkan__ribbon_set_active(&renderer->baked_ribbons))) {
		return;
	}
	//Host coherent, and this frame slot's previous draws are done, so the GPU reads it when this frame runs
	memcpy(renderer->view_uniforms[renderer->frame_index].mapped, uniforms, sizeof(tfx_impl_uniforms_t));
	vkCmdBindDescriptorSets(command_buffer, VK_PIPELINE_BIND_POINT_GRAPHICS, renderer->pipeline_layout, 0, 1, &renderer->descriptor_set, 0, NULL);

	tfx_impl_vulkan_billboard_push_t push = {};
	push.sampler_index = 0;
	push.uniform_index = renderer->frame_index;
	push.color_ramp_texture_index = renderer->color_ramps.index;
	push.image_data_index = renderer->image_data.index;
	push.properties_index = renderer->emitter_properties.index;
	push.prev_billboards_index = renderer->instance_buffers[renderer->current_instance_buffer ^ 1].index;
	//The layers that write depth draw first with blending off, so they hide the blended layers behind them
	for (int group = 0; group != 2; ++group) {
		const int writes_depth = group == 0;
		if (renderer->instance_count > 0) {
			tfx_impl_vulkan__render_billboards(renderer, command_buffer, &push, writes_depth);
		}
		if (tfx_impl_vulkan__ribbon_set_active(&renderer->live_ribbons)) {
			tfx_impl_vulkan__render_ribbons(renderer, command_buffer, &renderer->live_ribbons, writes_depth);
		}
		//Baked billboards come as one range across their layers, so they all draw blended and unlit
		if (!writes_depth && renderer->playback_sprite_count > 0) {
			tfx_impl_vulkan__render_baked_billboards(renderer, command_buffer, &push);
		}
		if (tfx_impl_vulkan__ribbon_set_active(&renderer->baked_ribbons)) {
			tfx_impl_vulkan__render_ribbons(renderer, command_buffer, &renderer->baked_ribbons, writes_depth);
		}
	}
}

static void tfx_impl_vulkan__release_ribbon_set(tfx_impl_vulkan_t *renderer, tfx_impl_vulkan_ribbon_set_t *set) {
	tfx_impl_vulkan__release_buffer(renderer, &set->segments);
	tfx_impl_vulkan__release_buffer(renderer, &set->instances);
	tfx_impl_vulkan__release_buffer(renderer, &set->emitters);
	tfx_impl_vulkan__release_buffer(renderer, &set->vertices);
	tfx_impl_vulkan__release_buffer(renderer, &set->indices);
	tfx_impl_vulkan__release_buffer(renderer, &set->graph_lookups);
}

void tfx_impl_vulkan_Shutdown(tfx_impl_vulkan_t *renderer) {
	if (!renderer->device) {
		return;
	}
	tfx_impl_vulkan_buffer_t *buffers[] = { &renderer->image_data, &renderer->emitter_properties, &renderer->instance_buffers[0], &renderer->instance_buffers[1],
		&renderer->particle_pool, &renderer->particle_upload, &renderer->baked_image_data, &renderer->baked_emitter_properties, &renderer->baked_sprites, &renderer->animation_offsets,
		&renderer->animation_instances, &renderer->playback_instances, &renderer->baked_ribbon_data, &renderer->baked_ribbon_properties, &renderer->animation_ribbon_offsets, &renderer->graph_lookups };
	for (tfxU32 index = 0; index != sizeof(buffers) / sizeof(buffers[0]); ++index) {
		tfx_impl_vulkan__release_buffer(renderer, buffers[index]);
	}
	tfx_impl_vulkan__release_ribbon_set(renderer, &renderer->live_ribbons);
	tfx_impl_vulkan__release_ribbon_set(renderer, &renderer->baked_ribbons);
	for (tfxU32 frame = 0; frame != TFX_IMPL_VULKAN_MAX_FRAMES_IN_FLIGHT; ++frame) {
		tfx_impl_vulkan__release_mapped_buffer(renderer, &renderer->upload_rings[frame]);
		tfx_impl_vulkan__release_mapped_buffer(renderer, &renderer->view_uniforms[frame]);
	}
	tfx_impl_vulkan__destroy_image(renderer, &renderer->shape_atlas);
	tfx_impl_vulkan__destroy_image(renderer, &renderer->color_ramps);
	tfx_impl_vulkan__destroy_image(renderer, &renderer->baked_color_ramps);
	//The GPU is idle, so every frame's leftovers can go now
	for (tfxU32 frame = 0; frame != TFX_IMPL_VULKAN_MAX_FRAMES_IN_FLIGHT; ++frame) {
		tfx_impl_vulkan__process_release_list(renderer, &renderer->release_lists[frame]);
		free(renderer->release_lists[frame].releases);
	}
	memset(renderer->release_lists, 0, sizeof(renderer->release_lists));
	free(renderer->color_ramp_layer_versions);
	free(renderer->free_buffer_indexes.indexes);
	free(renderer->free_image_indexes.indexes);
	renderer->color_ramp_layer_versions = NULL;
	memset(&renderer->free_buffer_indexes, 0, sizeof(tfx_impl_vulkan_index_list_t));
	memset(&renderer->free_image_indexes, 0, sizeof(tfx_impl_vulkan_index_list_t));
	VkPipeline pipelines[] = { renderer->particle_kernel, renderer->baked_billboard_pipeline, renderer->ribbon_compute_pipeline, renderer->sprite_playback_pipeline, renderer->ribbon_playback_pipeline };
	for (tfxU32 index = 0; index != sizeof(pipelines) / sizeof(pipelines[0]); ++index) {
		vkDestroyPipeline(renderer->device, pipelines[index], NULL);
	}
	for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
		for (int lit = 0; lit != 2; ++lit) {
			vkDestroyPipeline(renderer->device, renderer->billboard_pipelines[writes_depth][lit], NULL);
			vkDestroyPipeline(renderer->device, renderer->ribbon_pipelines[writes_depth][lit], NULL);
		}
	}
	vkDestroyPipelineLayout(renderer->device, renderer->pipeline_layout, NULL);
	vkDestroyDescriptorPool(renderer->device, renderer->descriptor_pool, NULL);
	vkDestroyDescriptorSetLayout(renderer->device, renderer->descriptor_set_layout, NULL);
	vkDestroySampler(renderer->device, renderer->sampler, NULL);
	memset(renderer, 0, sizeof(tfx_impl_vulkan_t));
}
