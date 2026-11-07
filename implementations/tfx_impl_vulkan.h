#ifndef TFX_IMPL_VULKAN_HEADER
#define TFX_IMPL_VULKAN_HEADER

/*
	TimelineFX Vulkan backend, with the shaders' bindless model: every buffer and image is a descriptor in one descriptor set, indexed
	by the push constants. Needs Vulkan 1.3 and the device features tfx_impl_vulkan_Init lists, and the shaders slangc builds from
	../Shaders/timelinefx.slang as SPIR-V without TFX_FIXED_BINDINGS, see README.md.

	Table of contents, search for the keyword to jump to the section:
	[Vulkan_types]		The renderer struct that owns every Vulkan object
	[Vulkan_resources]	Memory, descriptors, buffers, images, uploads, barriers and pipelines (in the .cpp)
	[Vulkan_functions]	Init, upload, compute, render and shutdown
*/

#include <vulkan/vulkan.h>
#include "tfx_impl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

//--- Vulkan_types

#define TFX_IMPL_VULKAN_MAX_FRAMES_IN_FLIGHT 3
#define TFX_IMPL_VULKAN_NO_DESCRIPTOR 0xFFFFFFFF

//A device local buffer and its slot in the storage buffer array, which is the index the shaders read it by, whatever its row type
typedef struct tfx_impl_vulkan_buffer_s {
	VkBuffer buffer;
	VkDeviceMemory memory;
	VkDeviceSize capacity;
	tfxU32 index;
} tfx_impl_vulkan_buffer_t;

//A host visible buffer that stays mapped: a frame slot's upload ring, or its view uniforms
typedef struct tfx_impl_vulkan_mapped_buffer_s {
	VkBuffer buffer;
	VkDeviceMemory memory;
	unsigned char *mapped;
	VkDeviceSize capacity;
	VkDeviceSize used;
} tfx_impl_vulkan_mapped_buffer_t;

//A texture array, its view and its slot in the image array
typedef struct tfx_impl_vulkan_image_s {
	VkImage image;
	VkDeviceMemory memory;
	VkImageView view;
	VkImageLayout layout;			//Of every layer and mip, between uploads
	tfxU32 width;
	tfxU32 height;
	tfxU32 layer_count;
	tfxU32 mip_count;
	tfxU32 index;
} tfx_impl_vulkan_image_t;

//One object and the descriptor slots it held, released when the frame slot that let it go comes around again
typedef struct tfx_impl_vulkan_release_s {
	VkBuffer buffer;
	VkImage image;
	VkImageView view;
	VkDeviceMemory memory;
	tfxU32 buffer_index;
	tfxU32 image_index;
} tfx_impl_vulkan_release_t;

typedef struct tfx_impl_vulkan_release_list_s {
	tfx_impl_vulkan_release_t *releases;
	tfxU32 count;
	tfxU32 capacity;
} tfx_impl_vulkan_release_list_t;

typedef struct tfx_impl_vulkan_index_list_s {
	tfxU32 *indexes;
	tfxU32 count;
} tfx_impl_vulkan_index_list_t;

//What ribbon_compute reads and writes for one source of ribbons: a stage's live ribbons or an animation manager's baked ones
typedef struct tfx_impl_vulkan_ribbon_set_s {
	tfx_impl_vulkan_buffer_t segments;
	tfx_impl_vulkan_buffer_t instances;		//tfx_ribbon_t, uploaded for live ribbons and written by ribbon_playback for baked ones
	tfx_impl_vulkan_buffer_t emitters;
	tfx_impl_vulkan_buffer_t vertices;		//Built by ribbon_compute every frame, then the draw's vertex buffer
	tfx_impl_vulkan_buffer_t indices;		//Built by ribbon_compute every frame, then the draw's index buffer
	tfx_impl_vulkan_buffer_t graph_lookups;	//Baked only: tfx_GetAnimationRibbonGraphLookupsBufferPointer. Live ribbons read the renderer's shared table
	//The buckets are walked every frame because the geometry is rebuilt with each frame's lerp. One of these is set while there's something to draw
	tfx_stage stage;
	tfx_ribbon_batch_offsets_t stage_offsets;
	tfx_animation_manager animation_manager;
} tfx_impl_vulkan_ribbon_set_t;

typedef struct tfx_impl_vulkan_s {
	VkDevice device;
	VkPhysicalDeviceMemoryProperties memory_properties;
	tfxU32 frames_in_flight;
	tfxU32 frame_index;
	tfx_impl_vulkan_release_list_t release_lists[TFX_IMPL_VULKAN_MAX_FRAMES_IN_FLIGHT];
	tfx_impl_vulkan_mapped_buffer_t upload_rings[TFX_IMPL_VULKAN_MAX_FRAMES_IN_FLIGHT];
	//tfx_impl_uniforms_t for each frame slot, at the slot's index in the uniform buffer array
	tfx_impl_vulkan_mapped_buffer_t view_uniforms[TFX_IMPL_VULKAN_MAX_FRAMES_IN_FLIGHT];

	//One set holds everything: the sampler at binding 0, images at 3, storage buffers at 5 and the view uniforms at 7
	VkDescriptorSetLayout descriptor_set_layout;
	VkDescriptorPool descriptor_pool;
	VkDescriptorSet descriptor_set;
	VkSampler sampler;						//A clamped linear sampler with mips, at index 0
	tfx_impl_vulkan_index_list_t free_buffer_indexes;
	tfx_impl_vulkan_index_list_t free_image_indexes;

	//Every pipeline shares one layout, whose push constant range fits the largest block
	VkPipelineLayout pipeline_layout;
	VkPipeline particle_kernel;				//simulate_particles from timelinefx.slang
	//Each layer picks its pipelines from its tfx_GetLayerFlags, by [writes depth][lit]. Baked billboards are always blended and unlit
	VkPipeline billboard_pipelines[2][2];
	VkPipeline baked_billboard_pipeline;
	VkPipeline ribbon_pipelines[2][2];
	VkPipeline ribbon_compute_pipeline;
	VkPipeline sprite_playback_pipeline;
	VkPipeline ribbon_playback_pipeline;

	tfx_impl_vulkan_image_t shape_atlas;	//Image index 0, because the uv lookup leaves the image index at 0
	//Every loaded library's data, from the global getters, so one set of buffers draws effects from any library
	tfx_impl_vulkan_image_t color_ramps;
	tfxU32 *color_ramp_layer_versions;		//The global layer versions last uploaded, one per layer of color_ramps
	tfx_impl_vulkan_buffer_t image_data;		//tfx_gpu_image_data_t per shape frame
	tfx_impl_vulkan_buffer_t emitter_properties;	//tfx_gpu_particle_properties_t per emitter
	tfxU32 properties_version;
	tfxU32 shapes_version;
	tfxU32 color_ramps_version;
	tfxU32 graph_lookups_version;

	//Two instance buffers that swap on every update: the latest instances and the ones captured_index points into
	tfx_impl_vulkan_buffer_t instance_buffers[2];
	tfxU32 current_instance_buffer;
	tfxU32 instance_capacity;
	tfxU32 instance_count;			//CPU and GPU instances
	//Each layer is drawn as two ranges, its CPU instances then its GPU instances
	tfxU32 layer_start[tfxLAYERS];
	tfxU32 layer_count[tfxLAYERS];
	tfxU32 gpu_layer_start[tfxLAYERS];
	tfxU32 gpu_layer_count[tfxLAYERS];
	tfxLayerFlags layer_flags[tfxLAYERS];
	tfxU32 library_layer_count;

	//The pool's hot, cold and noise sections in one buffer, and the frame's tfx_CopyGPUParticleUpload in another. Storage buffers
	//are untyped, so one index each serves every row type, with the push constants' bases
	tfx_impl_vulkan_buffer_t particle_pool;
	tfx_impl_vulkan_buffer_t particle_upload;
	tfx_gpu_particle_push_t particle_push_constants;
	//From the last upload and still valid until the next stage update, dispatched once by the next tfx_impl_vulkan_Compute
	const tfx_gpu_particle_dispatch_t *particle_dispatches;
	tfxU32 particle_dispatch_count;

	//Baked sprite data playback has its own color ramps, shape frames and emitter properties, from the animation manager
	tfx_impl_vulkan_image_t baked_color_ramps;
	tfx_impl_vulkan_buffer_t baked_image_data;
	tfx_impl_vulkan_buffer_t baked_emitter_properties;
	tfx_impl_vulkan_buffer_t baked_sprites;		//Every baked sprite of every effect
	tfx_impl_vulkan_buffer_t animation_offsets;	//Uploaded every frame
	tfx_impl_vulkan_buffer_t animation_instances;
	tfx_impl_vulkan_buffer_t playback_instances;	//The billboard instances the playback compute writes, then a vertex buffer
	tfx_sprite_data_push_t sprite_data_push;
	tfxU32 playback_sprite_count;
	tfx_impl_vulkan_buffer_t baked_ribbon_data;	//tfx_ribbon_instance_data_t, every baked ribbon of every effect
	tfx_impl_vulkan_buffer_t baked_ribbon_properties;
	tfx_impl_vulkan_buffer_t animation_ribbon_offsets;	//Uploaded every frame
	tfx_ribbon_data_push_t ribbon_data_push;
	tfxU32 playback_ribbon_count;

	tfx_impl_vulkan_buffer_t graph_lookups;		//tfx_GetGPUGraphLookupsBuffer
	tfx_impl_vulkan_ribbon_set_t live_ribbons;	//Uploaded when an update ran
	tfx_impl_vulkan_ribbon_set_t baked_ribbons;	//Segments and emitters uploaded once, the rest every frame
} tfx_impl_vulkan_t;

//--- Vulkan_functions

//Loads the shaders and creates the descriptor set, pipeline layout and pipelines. shader_folder holds each entry point of timelinefx.slang
//as <entry point>.spv, bindless, compiled with -fvk-use-entrypoint-name. The pipelines draw with dynamic rendering into a color attachment
//of color_format and a depth attachment of depth_format. The device needs Vulkan 1.3 with these features on:
//	1.1: shaderDrawParameters
//	1.2: runtimeDescriptorArray, descriptorBindingPartiallyBound, descriptorBindingSampledImageUpdateAfterBind,
//	     descriptorBindingStorageBufferUpdateAfterBind, descriptorBindingUpdateUnusedWhilePending, shaderSampledImageArrayNonUniformIndexing
//	1.3: dynamicRendering, synchronization2
bool tfx_impl_vulkan_Init(tfx_impl_vulkan_t *renderer, VkPhysicalDevice physical_device, VkDevice device, VkFormat color_format, VkFormat depth_format, const char *shader_folder, tfxU32 frames_in_flight);
//Call first every frame, once the GPU has finished the last frame that used frame_index (0 to frames_in_flight - 1). Releases what that frame
//replaced and starts its upload ring again
void tfx_impl_vulkan_BeginFrame(tfx_impl_vulkan_t *renderer, tfxU32 frame_index);
//Creates the texture array that every library's atlas and the sprite data's atlas are uploaded into, with room for all their layers. Only while the GPU is idle
bool tfx_impl_vulkan_CreateAtlas(tfx_impl_vulkan_t *renderer, tfxU32 layer_size, tfxU32 layer_count);
//Records an atlas's upload into its layers, and blits its mips. The rest of the library data goes up with each frame in tfx_impl_vulkan_UploadFrame
bool tfx_impl_vulkan_UploadAtlas(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_atlas_t *atlas);
//Records the upload of an animation manager's baked sprite data, color ramps, shape frames and emitter properties, after tfx_BuildAnimationManagerGPUShapeData
bool tfx_impl_vulkan_UploadSpriteData(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, tfx_animation_manager animation_manager, tfx_gpu_shapes shapes);
//Records the baked playback's upload every frame, and on update frames the instances, the ribbon data and any library data tfx_UpdateGlobalGPUData changed
bool tfx_impl_vulkan_UploadFrame(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_frame_t *frame);
//Records the GPU particle kernel, the baked playback and the ribbon geometry, and the barriers the draws need. Call every frame, outside
//vkCmdBeginRendering and before tfx_impl_vulkan_Render, even when there's nothing to dispatch
void tfx_impl_vulkan_Compute(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_uniforms_t *uniforms);
//Records the draws inside the rendering you've begun, with the viewport and scissor you've set, interpolating with uniforms->timer_lerp
void tfx_impl_vulkan_Render(tfx_impl_vulkan_t *renderer, VkCommandBuffer command_buffer, const tfx_impl_uniforms_t *uniforms);
//Once the GPU is idle
void tfx_impl_vulkan_Shutdown(tfx_impl_vulkan_t *renderer);

#ifdef __cplusplus
}
#endif

#endif
