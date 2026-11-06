#include "tfx_impl_dx12.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TFX_IMPL_DX12_INITIAL_INSTANCE_CAPACITY 10000
#define TFX_IMPL_DX12_DESCRIPTOR_CAPACITY 4096
#define TFX_IMPL_DX12_SAMPLER_CAPACITY 16
//Bigger uploads get a buffer of their own rather than growing a frame's ring
#define TFX_IMPL_DX12_RING_SIZE (4 * 1024 * 1024)
//tfx_GetGlobalColorRampLayerCount never goes past this
#define TFX_IMPL_DX12_MAX_COLOR_RAMP_LAYERS 256
//Match the thread group sizes in timelinefx.slang
#define TFX_IMPL_DX12_KERNEL_GROUP_SIZE 256
#define TFX_IMPL_DX12_RIBBON_GROUP_SIZE 256
#define TFX_IMPL_DX12_PLAYBACK_GROUP_SIZE 128
//The spaces of the bindless descriptor arrays in timelinefx.slang: t0 in spaces 1 to 25, u0 in 1 to 8, b0 and s0 in space 1
#define TFX_IMPL_DX12_READ_SPACES 25
#define TFX_IMPL_DX12_WRITE_SPACES 8
//Every state a buffer is read in, so a buffer goes to it once and serves any draw or dispatch
#define TFX_IMPL_DX12_SHADER_READ (D3D12_RESOURCE_STATE_VERTEX_AND_CONSTANT_BUFFER | D3D12_RESOURCE_STATE_INDEX_BUFFER | D3D12_RESOURCE_STATE_NON_PIXEL_SHADER_RESOURCE | D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE)
#define TFX_IMPL_DX12_RELEASE(object) if (object) { (object)->Release(); (object) = NULL; }

//Row sizes of the structs timelinefx.slang reads, which a structured view bakes in. The ones timelinefx.h doesn't declare are mirrored here
#define TFX_IMPL_DX12_INSTANCE_STRIDE ((tfxU32)sizeof(tfx_instance_t))
#define TFX_IMPL_DX12_IMAGE_DATA_STRIDE ((tfxU32)sizeof(tfx_gpu_image_data_t))
#define TFX_IMPL_DX12_PROPERTIES_STRIDE ((tfxU32)sizeof(tfx_gpu_particle_properties_t))
#define TFX_IMPL_DX12_GRAPH_STRIDE ((tfxU32)sizeof(tfx_gpu_graph_data_t))
#define TFX_IMPL_DX12_RIBBON_VERTEX_STRIDE ((tfxU32)sizeof(tfx_ribbon_vertex_t))
#define TFX_IMPL_DX12_RIBBON_STRIDE tfxRIBBON_SIZE
#define TFX_IMPL_DX12_RIBBON_EMITTER_STRIDE 416
#define TFX_IMPL_DX12_RIBBON_SEGMENT_STRIDE 16
#define TFX_IMPL_DX12_ANIMATION_INSTANCE_STRIDE 64
#define TFX_IMPL_DX12_BAKED_RIBBON_STRIDE 64
#define TFX_IMPL_DX12_BAKED_RIBBON_PROPERTIES_STRIDE 68
#define TFX_IMPL_DX12_UINT_STRIDE 4
#define TFX_IMPL_DX12_HOT_STRIDE 32
#define TFX_IMPL_DX12_COLD_STRIDE 64
#define TFX_IMPL_DX12_NOISE_STRIDE 16
#define TFX_IMPL_DX12_DISPATCH_STRIDE 128
#define TFX_IMPL_DX12_TICK_STRIDE 64
#define TFX_IMPL_DX12_SPAWN_STRIDE 112
#define TFX_IMPL_DX12_LOCATION_STRIDE 32

//timelinefx.slang mirrors these structs by hand, in its [Render_rows] section
static_assert(sizeof(tfx_instance_t) == 48, "tfx_instance_t changed, update tfx_instance in timelinefx.slang");
static_assert(sizeof(tfx_gpu_image_data_t) == 48, "tfx_gpu_image_data_t changed, update tfx_image_data in timelinefx.slang");
static_assert(sizeof(tfx_gpu_particle_properties_t) == 32, "tfx_gpu_particle_properties_t changed, update tfx_emitter_properties in timelinefx.slang");
static_assert(sizeof(tfx_ribbon_vertex_t) == 32, "tfx_ribbon_vertex_t changed, update tfx_ribbon_vertex_row in timelinefx.slang");
static_assert(sizeof(tfx_gpu_graph_data_t) == 48, "tfx_gpu_graph_data_t changed, update tfx_gpu_graph_data in timelinefx.slang");
static_assert(sizeof(tfx_impl_uniforms_t) <= 256, "The view uniforms have to fit one constant buffer view");

//tfx_billboard_push in timelinefx.slang, which is the host's to declare. Every index is a slot in the descriptor heap
typedef struct tfx_impl_dx12_billboard_push_s {
	tfxU32 color_ramp_texture_index;
	tfxU32 image_data_index;
	tfxU32 properties_index;
	tfxU32 sampler_index;
	tfxU32 prev_billboards_index;
	tfxU32 index_offset;
	tfxU32 uniform_index;
} tfx_impl_dx12_billboard_push_t;

//tfx_baked_billboard_push, the live push's fields and then the clip space offsets
typedef struct tfx_impl_dx12_baked_billboard_push_s {
	tfx_impl_dx12_billboard_push_t billboard;
	float ndc_offset_x;
	float ndc_offset_y;
} tfx_impl_dx12_baked_billboard_push_t;

//Root constant registers of the push constant blocks, b0 to b6 in space 0
typedef enum tfx_impl_dx12_push_register {
	tfx_impl_dx12_push_particles = 0,
	tfx_impl_dx12_push_billboard = 1,
	tfx_impl_dx12_push_baked_billboard = 2,
	tfx_impl_dx12_push_ribbon = 3,
	tfx_impl_dx12_push_sprite_data = 4,
	tfx_impl_dx12_push_ribbon_data = 5,
} tfx_impl_dx12_push_register;

//The root parameters every root signature starts with
#define TFX_IMPL_DX12_DESCRIPTOR_TABLE_PARAMETER 0
#define TFX_IMPL_DX12_SAMPLER_TABLE_PARAMETER 1
#define TFX_IMPL_DX12_PUSH_PARAMETER 2

//--- Dx12_resources

static void tfx_impl_dx12__grow_array(void **array, tfxU32 *capacity, tfxU32 needed, size_t element_size) {
	if (needed <= *capacity) {
		return;
	}
	tfxU32 new_capacity = *capacity ? *capacity * 2 : 16;
	while (new_capacity < needed) {
		new_capacity *= 2;
	}
	void *grown = realloc(*array, element_size * new_capacity);
	if (grown) {
		*array = grown;
		*capacity = new_capacity;
	}
}

//A frame still in flight may use it, so it goes when this frame slot comes around again
static void tfx_impl_dx12__defer_release(tfx_impl_dx12_t *renderer, IUnknown *object) {
	if (!object) {
		return;
	}
	tfx_impl_dx12_release_list_t *list = &renderer->release_lists[renderer->frame_index];
	tfx_impl_dx12__grow_array((void **)&list->objects, &list->object_capacity, list->object_count + 1, sizeof(IUnknown *));
	if (list->object_count < list->object_capacity) {
		list->objects[list->object_count++] = object;
	}
}

static void tfx_impl_dx12__defer_free_descriptor(tfx_impl_dx12_t *renderer, tfxU32 index) {
	if (index == TFX_IMPL_DX12_NO_DESCRIPTOR) {
		return;
	}
	tfx_impl_dx12_release_list_t *list = &renderer->release_lists[renderer->frame_index];
	tfx_impl_dx12__grow_array((void **)&list->descriptors, &list->descriptor_capacity, list->descriptor_count + 1, sizeof(tfxU32));
	if (list->descriptor_count < list->descriptor_capacity) {
		list->descriptors[list->descriptor_count++] = index;
	}
}

static void tfx_impl_dx12__process_release_list(tfx_impl_dx12_t *renderer, tfx_impl_dx12_release_list_t *list) {
	for (tfxU32 index = 0; index != list->object_count; ++index) {
		list->objects[index]->Release();
	}
	for (tfxU32 index = 0; index != list->descriptor_count; ++index) {
		renderer->free_descriptors[renderer->free_descriptor_count++] = list->descriptors[index];
	}
	list->object_count = 0;
	list->descriptor_count = 0;
}

static tfxU32 tfx_impl_dx12__allocate_descriptor(tfx_impl_dx12_t *renderer) {
	if (renderer->free_descriptor_count == 0) {
		tfx_Print("The descriptor heap's %u slots are all in use", renderer->descriptor_capacity);
		return TFX_IMPL_DX12_NO_DESCRIPTOR;
	}
	return renderer->free_descriptors[--renderer->free_descriptor_count];
}

static D3D12_CPU_DESCRIPTOR_HANDLE tfx_impl_dx12__descriptor(tfx_impl_dx12_t *renderer, tfxU32 index) {
	D3D12_CPU_DESCRIPTOR_HANDLE handle = renderer->descriptor_heap->GetCPUDescriptorHandleForHeapStart();
	handle.ptr += (SIZE_T)index * renderer->descriptor_size;
	return handle;
}

static tfxU32 tfx_impl_dx12__create_read_view(tfx_impl_dx12_t *renderer, ID3D12Resource *resource, tfxU32 size_in_bytes, tfxU32 stride) {
	tfxU32 index = tfx_impl_dx12__allocate_descriptor(renderer);
	if (index != TFX_IMPL_DX12_NO_DESCRIPTOR) {
		D3D12_SHADER_RESOURCE_VIEW_DESC description = {};
		description.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
		description.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
		description.Buffer.NumElements = size_in_bytes / stride;
		description.Buffer.StructureByteStride = stride;
		renderer->device->CreateShaderResourceView(resource, &description, tfx_impl_dx12__descriptor(renderer, index));
	}
	return index;
}

static tfxU32 tfx_impl_dx12__create_write_view(tfx_impl_dx12_t *renderer, ID3D12Resource *resource, tfxU32 size_in_bytes, tfxU32 stride) {
	tfxU32 index = tfx_impl_dx12__allocate_descriptor(renderer);
	if (index != TFX_IMPL_DX12_NO_DESCRIPTOR) {
		D3D12_UNORDERED_ACCESS_VIEW_DESC description = {};
		description.ViewDimension = D3D12_UAV_DIMENSION_BUFFER;
		description.Buffer.NumElements = size_in_bytes / stride;
		description.Buffer.StructureByteStride = stride;
		renderer->device->CreateUnorderedAccessView(resource, NULL, &description, tfx_impl_dx12__descriptor(renderer, index));
	}
	return index;
}

static void tfx_impl_dx12__create_texture_view(tfx_impl_dx12_t *renderer, ID3D12Resource *texture, DXGI_FORMAT format, tfxU32 layer_count, tfxU32 index) {
	D3D12_SHADER_RESOURCE_VIEW_DESC description = {};
	description.Format = format;
	description.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2DARRAY;
	description.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
	description.Texture2DArray.MipLevels = (UINT)-1;
	description.Texture2DArray.ArraySize = layer_count;
	renderer->device->CreateShaderResourceView(texture, &description, tfx_impl_dx12__descriptor(renderer, index));
}

static void tfx_impl_dx12__release_buffer(tfx_impl_dx12_t *renderer, tfx_impl_dx12_buffer_t *buffer) {
	tfx_impl_dx12__defer_release(renderer, buffer->resource);
	tfx_impl_dx12__defer_free_descriptor(renderer, buffer->shader_read_index);
	tfx_impl_dx12__defer_free_descriptor(renderer, buffer->shader_write_index);
	buffer->resource = NULL;
	buffer->capacity = 0;
	buffer->shader_read_index = TFX_IMPL_DX12_NO_DESCRIPTOR;
	buffer->shader_write_index = TFX_IMPL_DX12_NO_DESCRIPTOR;
}

static ID3D12Resource *tfx_impl_dx12__create_resource(ID3D12Device *device, D3D12_HEAP_TYPE heap_type, UINT64 size_in_bytes, D3D12_RESOURCE_FLAGS flags, D3D12_RESOURCE_STATES state) {
	D3D12_HEAP_PROPERTIES heap = {};
	heap.Type = heap_type;
	D3D12_RESOURCE_DESC description = {};
	description.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
	description.Width = size_in_bytes;
	description.Height = 1;
	description.DepthOrArraySize = 1;
	description.MipLevels = 1;
	description.SampleDesc.Count = 1;
	description.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
	description.Flags = flags;
	ID3D12Resource *resource = NULL;
	HRESULT result = device->CreateCommittedResource(&heap, D3D12_HEAP_FLAG_NONE, &description, state, NULL, IID_PPV_ARGS(&resource));
	if (FAILED(result)) {
		tfx_Print("Could not create a buffer of %llu bytes (0x%08lx)", (unsigned long long)size_in_bytes, result);
		return NULL;
	}
	return resource;
}

//A default heap buffer with a structured view at each stride that isn't 0. A write stride makes it writable by compute shaders
static bool tfx_impl_dx12__create_buffer(tfx_impl_dx12_t *renderer, tfxU32 size_in_bytes, tfxU32 read_stride, tfxU32 write_stride, tfx_impl_dx12_buffer_t *buffer) {
	tfx_impl_dx12__release_buffer(renderer, buffer);
	size_in_bytes = ((size_in_bytes ? size_in_bytes : 16) + 255) & ~255u;
	buffer->resource = tfx_impl_dx12__create_resource(renderer->device, D3D12_HEAP_TYPE_DEFAULT, size_in_bytes,
		write_stride ? D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS : D3D12_RESOURCE_FLAG_NONE, D3D12_RESOURCE_STATE_COMMON);
	if (!buffer->resource) {
		return false;
	}
	buffer->capacity = size_in_bytes;
	buffer->state = D3D12_RESOURCE_STATE_COMMON;
	if (read_stride) {
		buffer->shader_read_index = tfx_impl_dx12__create_read_view(renderer, buffer->resource, size_in_bytes, read_stride);
	}
	if (write_stride) {
		buffer->shader_write_index = tfx_impl_dx12__create_write_view(renderer, buffer->resource, size_in_bytes, write_stride);
	}
	return (!read_stride || buffer->shader_read_index != TFX_IMPL_DX12_NO_DESCRIPTOR) && (!write_stride || buffer->shader_write_index != TFX_IMPL_DX12_NO_DESCRIPTOR);
}

//Grows by half again whenever it's too small, so the size settles after the first few frames. Never empty, so a shader always has a view to read
static bool tfx_impl_dx12__ensure_buffer(tfx_impl_dx12_t *renderer, size_t size_in_bytes, tfxU32 read_stride, tfxU32 write_stride, tfx_impl_dx12_buffer_t *buffer) {
	tfxU32 size = (tfxU32)(size_in_bytes ? size_in_bytes : 16);
	if (buffer->resource && buffer->capacity >= size) {
		return true;
	}
	return tfx_impl_dx12__create_buffer(renderer, size + size / 2, read_stride, write_stride, buffer);
}

static void tfx_impl_dx12__transition(ID3D12GraphicsCommandList *command_list, ID3D12Resource *resource, D3D12_RESOURCE_STATES *state, D3D12_RESOURCE_STATES new_state) {
	if (!resource || *state == new_state) {
		return;
	}
	D3D12_RESOURCE_BARRIER barrier = {};
	barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
	barrier.Transition.pResource = resource;
	barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
	barrier.Transition.StateBefore = *state;
	barrier.Transition.StateAfter = new_state;
	command_list->ResourceBarrier(1, &barrier);
	*state = new_state;
}

static void tfx_impl_dx12__transition_buffer(ID3D12GraphicsCommandList *command_list, tfx_impl_dx12_buffer_t *buffer, D3D12_RESOURCE_STATES new_state) {
	tfx_impl_dx12__transition(command_list, buffer->resource, &buffer->state, new_state);
}

//Into unordered access, or a barrier between the last frame's writes and these if it's there already
static void tfx_impl_dx12__begin_writes(ID3D12GraphicsCommandList *command_list, tfx_impl_dx12_buffer_t *buffer) {
	if (buffer->state == D3D12_RESOURCE_STATE_UNORDERED_ACCESS) {
		D3D12_RESOURCE_BARRIER barrier = {};
		barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_UAV;
		barrier.UAV.pResource = buffer->resource;
		command_list->ResourceBarrier(1, &barrier);
	} else {
		tfx_impl_dx12__transition_buffer(command_list, buffer, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
	}
}

static ID3D12Resource *tfx_impl_dx12__create_upload_resource(ID3D12Device *device, size_t size_in_bytes, unsigned char **mapped) {
	ID3D12Resource *resource = tfx_impl_dx12__create_resource(device, D3D12_HEAP_TYPE_UPLOAD, size_in_bytes, D3D12_RESOURCE_FLAG_NONE, D3D12_RESOURCE_STATE_GENERIC_READ);
	D3D12_RANGE no_reads = { 0, 0 };
	if (resource && FAILED(resource->Map(0, &no_reads, (void **)mapped))) {
		resource->Release();
		return NULL;
	}
	return resource;
}

//Upload memory for this frame. Small uploads share the frame slot's ring, big ones get a buffer of their own that goes when the slot comes round
static unsigned char *tfx_impl_dx12__allocate_upload(tfx_impl_dx12_t *renderer, size_t size_in_bytes, size_t alignment, ID3D12Resource **resource, UINT64 *offset) {
	tfx_impl_dx12_upload_ring_t *ring = &renderer->upload_rings[renderer->frame_index];
	size_t start = (ring->used + alignment - 1) / alignment * alignment;
	if (!ring->resource || start + size_in_bytes > ring->capacity) {
		if (size_in_bytes > TFX_IMPL_DX12_RING_SIZE / 4) {
			unsigned char *mapped = NULL;
			*resource = tfx_impl_dx12__create_upload_resource(renderer->device, size_in_bytes, &mapped);
			tfx_impl_dx12__defer_release(renderer, *resource);
			*offset = 0;
			return *resource ? mapped : NULL;
		}
		tfx_impl_dx12__defer_release(renderer, ring->resource);
		ring->capacity = ring->capacity ? ring->capacity * 2 : TFX_IMPL_DX12_RING_SIZE;
		ring->resource = tfx_impl_dx12__create_upload_resource(renderer->device, ring->capacity, &ring->mapped);
		ring->used = 0;
		start = 0;
		if (!ring->resource) {
			ring->capacity = 0;
			return NULL;
		}
	}
	ring->used = start + size_in_bytes;
	*resource = ring->resource;
	*offset = start;
	return ring->mapped + start;
}

//Leaves the buffer as a copy destination, it goes to the state it's used in when it's bound
static bool tfx_impl_dx12__upload_buffer(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, tfx_impl_dx12_buffer_t *buffer, const void *data, size_t size_in_bytes) {
	if (size_in_bytes == 0) {
		return true;
	}
	ID3D12Resource *source = NULL;
	UINT64 offset = 0;
	unsigned char *mapped = tfx_impl_dx12__allocate_upload(renderer, size_in_bytes, 16, &source, &offset);
	if (!mapped) {
		return false;
	}
	memcpy(mapped, data, size_in_bytes);
	tfx_impl_dx12__transition_buffer(command_list, buffer, D3D12_RESOURCE_STATE_COPY_DEST);
	command_list->CopyBufferRegion(buffer->resource, 0, source, offset, size_in_bytes);
	return true;
}

static bool tfx_impl_dx12__write_buffer(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, tfx_impl_dx12_buffer_t *buffer, const void *data, size_t size_in_bytes, tfxU32 read_stride) {
	return tfx_impl_dx12__ensure_buffer(renderer, size_in_bytes, read_stride, 0, buffer) && tfx_impl_dx12__upload_buffer(renderer, command_list, buffer, data, size_in_bytes);
}

//Texture rows are copied at a pitch of 256 bytes, from an offset aligned to 512
static bool tfx_impl_dx12__upload_texture(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, ID3D12Resource *texture, UINT subresource, DXGI_FORMAT format,
	tfxU32 width, tfxU32 height, tfxU32 bytes_per_pixel, const void *pixels) {
	tfxU32 row_bytes = width * bytes_per_pixel;
	tfxU32 pitch = (row_bytes + D3D12_TEXTURE_DATA_PITCH_ALIGNMENT - 1) & ~(D3D12_TEXTURE_DATA_PITCH_ALIGNMENT - 1);
	ID3D12Resource *source = NULL;
	UINT64 offset = 0;
	unsigned char *mapped = tfx_impl_dx12__allocate_upload(renderer, (size_t)pitch * height, D3D12_TEXTURE_DATA_PLACEMENT_ALIGNMENT, &source, &offset);
	if (!mapped) {
		return false;
	}
	for (tfxU32 row = 0; row != height; ++row) {
		memcpy(mapped + (size_t)row * pitch, (const unsigned char *)pixels + (size_t)row * row_bytes, row_bytes);
	}
	D3D12_TEXTURE_COPY_LOCATION destination = {};
	destination.pResource = texture;
	destination.Type = D3D12_TEXTURE_COPY_TYPE_SUBRESOURCE_INDEX;
	destination.SubresourceIndex = subresource;
	D3D12_TEXTURE_COPY_LOCATION source_location = {};
	source_location.pResource = source;
	source_location.Type = D3D12_TEXTURE_COPY_TYPE_PLACED_FOOTPRINT;
	source_location.PlacedFootprint.Offset = offset;
	source_location.PlacedFootprint.Footprint.Format = format;
	source_location.PlacedFootprint.Footprint.Width = width;
	source_location.PlacedFootprint.Footprint.Height = height;
	source_location.PlacedFootprint.Footprint.Depth = 1;
	source_location.PlacedFootprint.Footprint.RowPitch = pitch;
	command_list->CopyTextureRegion(&destination, 0, 0, 0, &source_location, NULL);
	return true;
}

static ID3D12Resource *tfx_impl_dx12__create_texture_array(ID3D12Device *device, tfxU32 width, tfxU32 height, tfxU32 layer_count, tfxU32 mip_count, DXGI_FORMAT format) {
	D3D12_HEAP_PROPERTIES heap = {};
	heap.Type = D3D12_HEAP_TYPE_DEFAULT;
	D3D12_RESOURCE_DESC description = {};
	description.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;
	description.Width = width;
	description.Height = height;
	description.DepthOrArraySize = (UINT16)layer_count;
	description.MipLevels = (UINT16)mip_count;
	description.Format = format;
	description.SampleDesc.Count = 1;
	ID3D12Resource *texture = NULL;
	HRESULT result = device->CreateCommittedResource(&heap, D3D12_HEAP_FLAG_NONE, &description, D3D12_RESOURCE_STATE_COPY_DEST, NULL, IID_PPV_ARGS(&texture));
	if (FAILED(result)) {
		tfx_Print("Could not create a %u x %u x %u texture array (0x%08lx)", width, height, layer_count, result);
		return NULL;
	}
	return texture;
}

static DXGI_FORMAT tfx_impl_dx12__color_ramp_format(tfxU32 bytes_per_pixel) {
	if (bytes_per_pixel == 4) {
		return DXGI_FORMAT_R8G8B8A8_UNORM;
	} else if (bytes_per_pixel == 16) {
		return DXGI_FORMAT_R32G32B32A32_FLOAT;
	}
	return DXGI_FORMAT_R16G16B16A16_FLOAT;
}

static tfxU32 tfx_impl_dx12__bitmap_bytes_per_pixel(tfx_bitmap_t *bitmap) {
	return (tfxU32)(tfx_GetBitmapSize(bitmap) / ((size_t)tfx_GetBitmapWidth(bitmap) * tfx_GetBitmapHeight(bitmap)));
}

static bool tfx_impl_dx12__upload_color_ramp_layer(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, ID3D12Resource *texture, tfxU32 layer, tfx_bitmap_t *bitmap) {
	tfxU32 bytes_per_pixel = tfx_impl_dx12__bitmap_bytes_per_pixel(bitmap);
	return tfx_impl_dx12__upload_texture(renderer, command_list, texture, layer, tfx_impl_dx12__color_ramp_format(bytes_per_pixel),
		(tfxU32)tfx_GetBitmapWidth(bitmap), (tfxU32)tfx_GetBitmapHeight(bitmap), bytes_per_pixel, tfx_GetBitmapData(bitmap));
}

//A 2 x 2 box filter, the same as the mips D3D11 generates. D3D12 has no GenerateMips
static void tfx_impl_dx12__downsample(const unsigned char *source, tfxU32 source_size, unsigned char *destination) {
	tfxU32 size = source_size > 1 ? source_size / 2 : 1;
	tfxU32 last = source_size - 1;
	for (tfxU32 y = 0; y != size; ++y) {
		const unsigned char *top = source + (size_t)(y * 2 < last ? y * 2 : last) * source_size * 4;
		const unsigned char *bottom = source + (size_t)(y * 2 + 1 < last ? y * 2 + 1 : last) * source_size * 4;
		for (tfxU32 x = 0; x != size; ++x) {
			tfxU32 left = (x * 2 < last ? x * 2 : last) * 4;
			tfxU32 right = (x * 2 + 1 < last ? x * 2 + 1 : last) * 4;
			for (tfxU32 channel = 0; channel != 4; ++channel) {
				tfxU32 sum = top[left + channel] + top[right + channel] + bottom[left + channel] + bottom[right + channel];
				*destination++ = (unsigned char)((sum + 2) / 4);
			}
		}
	}
}

//slangc compiles every shader ahead of time to <entry point>.dxil in the shader folder. The caller frees the code
static void *tfx_impl_dx12__load_shader_code(const char *shader_folder, const char *entry_point, size_t *size) {
	char path[MAX_PATH];
	snprintf(path, sizeof(path), "%s/%s.dxil", shader_folder, entry_point);
	FILE *file = fopen(path, "rb");
	if (!file) {
		tfx_Print("Could not open the shader %s", path);
		return NULL;
	}
	fseek(file, 0, SEEK_END);
	long file_size = ftell(file);
	fseek(file, 0, SEEK_SET);
	void *code = file_size > 0 ? malloc((size_t)file_size) : NULL;
	if (code && fread(code, 1, (size_t)file_size, file) != (size_t)file_size) {
		free(code);
		code = NULL;
	}
	fclose(file);
	if (!code) {
		tfx_Print("Could not read the shader %s", path);
		return NULL;
	}
	*size = (size_t)file_size;
	return code;
}

//The descriptor table covers the whole heap once per bindless array: every range starts at the heap's first slot, so a slot's index is
//the same whichever array a shader reads it through. Then the sampler table, then one root constant block per push constant register
static ID3D12RootSignature *tfx_impl_dx12__create_root_signature(ID3D12Device *device, const tfx_impl_dx12_push_register *registers, const tfxU32 *dword_counts, tfxU32 push_count, bool graphics) {
	D3D12_DESCRIPTOR_RANGE ranges[2 + TFX_IMPL_DX12_READ_SPACES + TFX_IMPL_DX12_WRITE_SPACES] = {};
	tfxU32 range_count = 0;
	//Space 1 is the view the vertex shaders read and space 2 the same view as the lit fragments read it, with the light
	for (tfxU32 space = 1; space <= 2; ++space) {
		ranges[range_count].RangeType = D3D12_DESCRIPTOR_RANGE_TYPE_CBV;
		ranges[range_count].NumDescriptors = UINT_MAX;
		ranges[range_count].RegisterSpace = space;
		range_count++;
	}
	for (tfxU32 space = 1; space <= TFX_IMPL_DX12_READ_SPACES; ++space) {
		ranges[range_count].RangeType = D3D12_DESCRIPTOR_RANGE_TYPE_SRV;
		ranges[range_count].NumDescriptors = UINT_MAX;
		ranges[range_count].RegisterSpace = space;
		range_count++;
	}
	for (tfxU32 space = 1; space <= TFX_IMPL_DX12_WRITE_SPACES; ++space) {
		ranges[range_count].RangeType = D3D12_DESCRIPTOR_RANGE_TYPE_UAV;
		ranges[range_count].NumDescriptors = UINT_MAX;
		ranges[range_count].RegisterSpace = space;
		range_count++;
	}
	D3D12_DESCRIPTOR_RANGE sampler_range = {};
	sampler_range.RangeType = D3D12_DESCRIPTOR_RANGE_TYPE_SAMPLER;
	sampler_range.NumDescriptors = UINT_MAX;
	sampler_range.RegisterSpace = 1;

	D3D12_ROOT_PARAMETER parameters[TFX_IMPL_DX12_PUSH_PARAMETER + 2] = {};
	parameters[TFX_IMPL_DX12_DESCRIPTOR_TABLE_PARAMETER].ParameterType = D3D12_ROOT_PARAMETER_TYPE_DESCRIPTOR_TABLE;
	parameters[TFX_IMPL_DX12_DESCRIPTOR_TABLE_PARAMETER].DescriptorTable.NumDescriptorRanges = range_count;
	parameters[TFX_IMPL_DX12_DESCRIPTOR_TABLE_PARAMETER].DescriptorTable.pDescriptorRanges = ranges;
	parameters[TFX_IMPL_DX12_SAMPLER_TABLE_PARAMETER].ParameterType = D3D12_ROOT_PARAMETER_TYPE_DESCRIPTOR_TABLE;
	parameters[TFX_IMPL_DX12_SAMPLER_TABLE_PARAMETER].DescriptorTable.NumDescriptorRanges = 1;
	parameters[TFX_IMPL_DX12_SAMPLER_TABLE_PARAMETER].DescriptorTable.pDescriptorRanges = &sampler_range;
	for (tfxU32 index = 0; index != push_count; ++index) {
		D3D12_ROOT_PARAMETER *parameter = &parameters[TFX_IMPL_DX12_PUSH_PARAMETER + index];
		parameter->ParameterType = D3D12_ROOT_PARAMETER_TYPE_32BIT_CONSTANTS;
		parameter->Constants.ShaderRegister = (UINT)registers[index];
		parameter->Constants.Num32BitValues = dword_counts[index];
	}
	D3D12_ROOT_SIGNATURE_DESC description = {};
	description.NumParameters = TFX_IMPL_DX12_PUSH_PARAMETER + push_count;
	description.pParameters = parameters;
	description.Flags = graphics ? D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT : D3D12_ROOT_SIGNATURE_FLAG_NONE;
	ID3DBlob *blob = NULL;
	ID3DBlob *errors = NULL;
	ID3D12RootSignature *root_signature = NULL;
	HRESULT result = D3D12SerializeRootSignature(&description, D3D_ROOT_SIGNATURE_VERSION_1_0, &blob, &errors);
	if (SUCCEEDED(result)) {
		result = device->CreateRootSignature(0, blob->GetBufferPointer(), blob->GetBufferSize(), IID_PPV_ARGS(&root_signature));
	}
	if (FAILED(result)) {
		tfx_Print("Could not create a root signature (0x%08lx) %s", result, errors ? (const char *)errors->GetBufferPointer() : "");
	}
	TFX_IMPL_DX12_RELEASE(blob);
	TFX_IMPL_DX12_RELEASE(errors);
	return root_signature;
}

static ID3D12PipelineState *tfx_impl_dx12__create_compute_pipeline(ID3D12Device *device, ID3D12RootSignature *root_signature, const char *shader_folder, const char *entry_point) {
	size_t size = 0;
	void *code = tfx_impl_dx12__load_shader_code(shader_folder, entry_point, &size);
	if (!code) {
		return NULL;
	}
	D3D12_COMPUTE_PIPELINE_STATE_DESC description = {};
	description.pRootSignature = root_signature;
	description.CS.pShaderBytecode = code;
	description.CS.BytecodeLength = size;
	ID3D12PipelineState *pipeline = NULL;
	HRESULT result = device->CreateComputePipelineState(&description, IID_PPV_ARGS(&pipeline));
	free(code);
	if (FAILED(result)) {
		tfx_Print("Could not create the %s pipeline (0x%08lx)", entry_point, result);
	}
	return pipeline;
}

//Premultiplied blending, depth tested but not written, no culling
//Premultiplied blending, depth tested but not written, or for a layer that writes depth, blending off and depth written
static ID3D12PipelineState *tfx_impl_dx12__create_graphics_pipeline(ID3D12Device *device, ID3D12RootSignature *root_signature, const char *shader_folder, const char *vertex_entry_point,
	const char *fragment_entry_point, const D3D12_INPUT_ELEMENT_DESC *elements, tfxU32 element_count, DXGI_FORMAT render_target_format, DXGI_FORMAT depth_format, bool writes_depth) {
	size_t vertex_size = 0;
	size_t fragment_size = 0;
	void *vertex_code = tfx_impl_dx12__load_shader_code(shader_folder, vertex_entry_point, &vertex_size);
	void *fragment_code = tfx_impl_dx12__load_shader_code(shader_folder, fragment_entry_point, &fragment_size);
	ID3D12PipelineState *pipeline = NULL;
	if (vertex_code && fragment_code) {
		D3D12_GRAPHICS_PIPELINE_STATE_DESC description = {};
		description.pRootSignature = root_signature;
		description.VS.pShaderBytecode = vertex_code;
		description.VS.BytecodeLength = vertex_size;
		description.PS.pShaderBytecode = fragment_code;
		description.PS.BytecodeLength = fragment_size;
		D3D12_RENDER_TARGET_BLEND_DESC *blend = &description.BlendState.RenderTarget[0];
		blend->BlendEnable = writes_depth ? FALSE : TRUE;
		blend->SrcBlend = D3D12_BLEND_ONE;
		blend->DestBlend = D3D12_BLEND_INV_SRC_ALPHA;
		blend->BlendOp = D3D12_BLEND_OP_ADD;
		blend->SrcBlendAlpha = D3D12_BLEND_ONE;
		blend->DestBlendAlpha = D3D12_BLEND_INV_SRC_ALPHA;
		blend->BlendOpAlpha = D3D12_BLEND_OP_ADD;
		blend->LogicOp = D3D12_LOGIC_OP_NOOP;
		blend->RenderTargetWriteMask = D3D12_COLOR_WRITE_ENABLE_ALL;
		description.SampleMask = UINT_MAX;
		description.RasterizerState.FillMode = D3D12_FILL_MODE_SOLID;
		description.RasterizerState.CullMode = D3D12_CULL_MODE_NONE;
		description.RasterizerState.DepthClipEnable = TRUE;
		description.DepthStencilState.DepthEnable = TRUE;
		description.DepthStencilState.DepthWriteMask = writes_depth ? D3D12_DEPTH_WRITE_MASK_ALL : D3D12_DEPTH_WRITE_MASK_ZERO;
		description.DepthStencilState.DepthFunc = D3D12_COMPARISON_FUNC_LESS_EQUAL;
		description.InputLayout.pInputElementDescs = elements;
		description.InputLayout.NumElements = element_count;
		description.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
		description.NumRenderTargets = 1;
		description.RTVFormats[0] = render_target_format;
		description.DSVFormat = depth_format;
		description.SampleDesc.Count = 1;
		HRESULT result = device->CreateGraphicsPipelineState(&description, IID_PPV_ARGS(&pipeline));
		if (FAILED(result)) {
			tfx_Print("Could not create the %s pipeline (0x%08lx)", vertex_entry_point, result);
		}
	}
	free(vertex_code);
	free(fragment_code);
	return pipeline;
}

//--- Dx12_functions

bool tfx_impl_dx12_Init(tfx_impl_dx12_t *renderer, ID3D12Device *device, DXGI_FORMAT render_target_format, DXGI_FORMAT depth_format, const char *shader_folder, tfxU32 frames_in_flight) {
	memset(renderer, 0, sizeof(tfx_impl_dx12_t));
	renderer->device = device;
	device->AddRef();
	renderer->frames_in_flight = frames_in_flight < TFX_IMPL_DX12_MAX_FRAMES_IN_FLIGHT ? frames_in_flight : TFX_IMPL_DX12_MAX_FRAMES_IN_FLIGHT;

	//Unbounded unordered access and constant buffer ranges need tier 3
	D3D12_FEATURE_DATA_D3D12_OPTIONS options = {};
	if (FAILED(device->CheckFeatureSupport(D3D12_FEATURE_D3D12_OPTIONS, &options, sizeof(options))) || options.ResourceBindingTier < D3D12_RESOURCE_BINDING_TIER_3) {
		tfx_Print("The bindless shaders need resource binding tier 3, this device has tier %i", (int)options.ResourceBindingTier);
		return false;
	}

	D3D12_DESCRIPTOR_HEAP_DESC heap_description = {};
	heap_description.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
	heap_description.NumDescriptors = TFX_IMPL_DX12_DESCRIPTOR_CAPACITY;
	heap_description.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
	D3D12_DESCRIPTOR_HEAP_DESC sampler_heap_description = {};
	sampler_heap_description.Type = D3D12_DESCRIPTOR_HEAP_TYPE_SAMPLER;
	sampler_heap_description.NumDescriptors = TFX_IMPL_DX12_SAMPLER_CAPACITY;
	sampler_heap_description.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
	renderer->free_descriptors = (tfxU32 *)malloc(sizeof(tfxU32) * TFX_IMPL_DX12_DESCRIPTOR_CAPACITY);
	if (!renderer->free_descriptors || FAILED(device->CreateDescriptorHeap(&heap_description, IID_PPV_ARGS(&renderer->descriptor_heap)))
		|| FAILED(device->CreateDescriptorHeap(&sampler_heap_description, IID_PPV_ARGS(&renderer->sampler_heap)))) {
		tfx_Print("Could not create the descriptor heaps");
		return false;
	}
	renderer->descriptor_size = device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);
	renderer->descriptor_capacity = TFX_IMPL_DX12_DESCRIPTOR_CAPACITY;
	//Slot 0 is kept for the atlas, and the rest are handed out lowest first
	for (tfxU32 index = TFX_IMPL_DX12_DESCRIPTOR_CAPACITY - 1; index != 0; --index) {
		renderer->free_descriptors[renderer->free_descriptor_count++] = index;
	}
	D3D12_SAMPLER_DESC sampler = {};
	sampler.Filter = D3D12_FILTER_MIN_MAG_MIP_LINEAR;
	sampler.AddressU = D3D12_TEXTURE_ADDRESS_MODE_CLAMP;
	sampler.AddressV = D3D12_TEXTURE_ADDRESS_MODE_CLAMP;
	sampler.AddressW = D3D12_TEXTURE_ADDRESS_MODE_CLAMP;
	sampler.MaxLOD = D3D12_FLOAT32_MAX;
	device->CreateSampler(&sampler, renderer->sampler_heap->GetCPUDescriptorHandleForHeapStart());

	//One root signature per push constant block. The baked billboards' vertex shader reads b2 and the fragment shader it shares b1
	const tfx_impl_dx12_push_register particle_registers[] = { tfx_impl_dx12_push_particles };
	const tfxU32 particle_dwords[] = { sizeof(tfx_gpu_particle_push_t) / 4 };
	const tfx_impl_dx12_push_register billboard_registers[] = { tfx_impl_dx12_push_billboard };
	const tfxU32 billboard_dwords[] = { sizeof(tfx_impl_dx12_billboard_push_t) / 4 };
	const tfx_impl_dx12_push_register baked_registers[] = { tfx_impl_dx12_push_billboard, tfx_impl_dx12_push_baked_billboard };
	const tfxU32 baked_dwords[] = { sizeof(tfx_impl_dx12_billboard_push_t) / 4, sizeof(tfx_impl_dx12_baked_billboard_push_t) / 4 };
	const tfx_impl_dx12_push_register ribbon_registers[] = { tfx_impl_dx12_push_ribbon };
	const tfxU32 ribbon_dwords[] = { sizeof(tfx_ribbon_bucket_globals_t) / 4 };
	const tfx_impl_dx12_push_register sprite_registers[] = { tfx_impl_dx12_push_sprite_data };
	const tfxU32 sprite_dwords[] = { sizeof(tfx_sprite_data_push_t) / 4 };
	const tfx_impl_dx12_push_register ribbon_data_registers[] = { tfx_impl_dx12_push_ribbon_data };
	const tfxU32 ribbon_data_dwords[] = { sizeof(tfx_ribbon_data_push_t) / 4 };
	renderer->particle_root_signature = tfx_impl_dx12__create_root_signature(device, particle_registers, particle_dwords, 1, false);
	renderer->billboard_root_signature = tfx_impl_dx12__create_root_signature(device, billboard_registers, billboard_dwords, 1, true);
	renderer->baked_billboard_root_signature = tfx_impl_dx12__create_root_signature(device, baked_registers, baked_dwords, 2, true);
	renderer->ribbon_root_signature = tfx_impl_dx12__create_root_signature(device, ribbon_registers, ribbon_dwords, 1, true);
	renderer->sprite_playback_root_signature = tfx_impl_dx12__create_root_signature(device, sprite_registers, sprite_dwords, 1, false);
	renderer->ribbon_playback_root_signature = tfx_impl_dx12__create_root_signature(device, ribbon_data_registers, ribbon_data_dwords, 1, false);
	if (!renderer->particle_root_signature || !renderer->billboard_root_signature || !renderer->baked_billboard_root_signature || !renderer->ribbon_root_signature
		|| !renderer->sprite_playback_root_signature || !renderer->ribbon_playback_root_signature) {
		return false;
	}

	//The instances and ribbon vertexes are vertex buffers in the bindless model, in the formats the Shaders readme lists
	const D3D12_INPUT_ELEMENT_DESC billboard_elements[] = {
		{ "POSITION", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, offsetof(tfx_instance_t, position), D3D12_INPUT_CLASSIFICATION_PER_INSTANCE_DATA, 1 },
		{ "QUATERNION", 0, DXGI_FORMAT_R16G16B16A16_SNORM, 0, offsetof(tfx_instance_t, quaternion), D3D12_INPUT_CLASSIFICATION_PER_INSTANCE_DATA, 1 },
		{ "SIZE", 0, DXGI_FORMAT_R16G16_SINT, 0, offsetof(tfx_instance_t, size), D3D12_INPUT_CLASSIFICATION_PER_INSTANCE_DATA, 1 },
		{ "ALIGNMENT", 0, DXGI_FORMAT_R8G8B8A8_SNORM, 0, offsetof(tfx_instance_t, alignment), D3D12_INPUT_CLASSIFICATION_PER_INSTANCE_DATA, 1 },
		{ "INTENSITY", 0, DXGI_FORMAT_R16G16_SINT, 0, offsetof(tfx_instance_t, intensity_gradient_map), D3D12_INPUT_CLASSIFICATION_PER_INSTANCE_DATA, 1 },
		{ "CURVED_ALPHA", 0, DXGI_FORMAT_R8G8B8A8_UNORM, 0, offsetof(tfx_instance_t, curved_alpha_life), D3D12_INPUT_CLASSIFICATION_PER_INSTANCE_DATA, 1 },
		{ "INDEXES", 0, DXGI_FORMAT_R32_UINT, 0, offsetof(tfx_instance_t, indexes), D3D12_INPUT_CLASSIFICATION_PER_INSTANCE_DATA, 1 },
		{ "CAPTURED_INDEX", 0, DXGI_FORMAT_R32_UINT, 0, offsetof(tfx_instance_t, captured_index), D3D12_INPUT_CLASSIFICATION_PER_INSTANCE_DATA, 1 },
	};
	const D3D12_INPUT_ELEMENT_DESC ribbon_elements[] = {
		{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, offsetof(tfx_ribbon_vertex_t, position), D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
		{ "SEGMENT_INDEX", 0, DXGI_FORMAT_R32_UINT, 0, offsetof(tfx_ribbon_vertex_t, segment_index), D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
		{ "UV_OFFSET_SCALE", 0, DXGI_FORMAT_R32G32_FLOAT, 0, offsetof(tfx_ribbon_vertex_t, uv_offset_scale), D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
		{ "RIBBON_INDEX", 0, DXGI_FORMAT_R32_UINT, 0, offsetof(tfx_ribbon_vertex_t, ribbon_index), D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
		{ "NORMAL", 0, DXGI_FORMAT_R32_UINT, 0, offsetof(tfx_ribbon_vertex_t, normal), D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },	//Only the lit vertex shader reads it
	};
	for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
		for (int lit = 0; lit != 2; ++lit) {
			renderer->billboard_pipelines[writes_depth][lit] = tfx_impl_dx12__create_graphics_pipeline(device, renderer->billboard_root_signature, shader_folder, tfx_impl_billboard_vertex_entry_points[lit],
				tfx_impl_billboard_fragment_entry_points[writes_depth][lit], billboard_elements, 8, render_target_format, depth_format, writes_depth != 0);
			renderer->ribbon_pipelines[writes_depth][lit] = tfx_impl_dx12__create_graphics_pipeline(device, renderer->ribbon_root_signature, shader_folder, tfx_impl_ribbon_vertex_entry_points[lit],
				tfx_impl_ribbon_fragment_entry_points[writes_depth][lit], ribbon_elements, lit ? 5 : 4, render_target_format, depth_format, writes_depth != 0);
			if (!renderer->billboard_pipelines[writes_depth][lit] || !renderer->ribbon_pipelines[writes_depth][lit]) {
				return false;
			}
		}
	}
	renderer->baked_billboard_pipeline = tfx_impl_dx12__create_graphics_pipeline(device, renderer->baked_billboard_root_signature, shader_folder, "tfx_billboard_baked_vertex", "tfx_billboard_fragment",
		billboard_elements, 8, render_target_format, depth_format, false);
	renderer->particle_kernel = tfx_impl_dx12__create_compute_pipeline(device, renderer->particle_root_signature, shader_folder, "simulate_particles");
	renderer->ribbon_compute_pipeline = tfx_impl_dx12__create_compute_pipeline(device, renderer->ribbon_root_signature, shader_folder, "tfx_ribbon_compute");
	renderer->sprite_playback_pipeline = tfx_impl_dx12__create_compute_pipeline(device, renderer->sprite_playback_root_signature, shader_folder, "tfx_sprite_playback");
	renderer->ribbon_playback_pipeline = tfx_impl_dx12__create_compute_pipeline(device, renderer->ribbon_playback_root_signature, shader_folder, "tfx_ribbon_playback");
	if (!renderer->baked_billboard_pipeline || !renderer->particle_kernel
		|| !renderer->ribbon_compute_pipeline || !renderer->sprite_playback_pipeline || !renderer->ribbon_playback_pipeline) {
		return false;
	}

	//Every buffer starts without views, so a release before the first create frees nothing
	tfx_impl_dx12_buffer_t *buffers[] = { &renderer->view_uniforms, &renderer->image_data, &renderer->emitter_properties, &renderer->instance_buffers[0], &renderer->instance_buffers[1],
		&renderer->particle_pool, &renderer->particle_upload, &renderer->baked_image_data, &renderer->baked_emitter_properties, &renderer->baked_sprites, &renderer->animation_offsets,
		&renderer->animation_instances, &renderer->playback_instances, &renderer->baked_ribbon_data, &renderer->baked_ribbon_properties, &renderer->animation_ribbon_offsets, &renderer->graph_lookups,
		&renderer->live_ribbons.segments, &renderer->live_ribbons.instances, &renderer->live_ribbons.emitters, &renderer->live_ribbons.vertices, &renderer->live_ribbons.indices, &renderer->live_ribbons.graph_lookups,
		&renderer->baked_ribbons.segments, &renderer->baked_ribbons.instances, &renderer->baked_ribbons.emitters, &renderer->baked_ribbons.vertices, &renderer->baked_ribbons.indices, &renderer->baked_ribbons.graph_lookups };
	for (tfxU32 index = 0; index != sizeof(buffers) / sizeof(buffers[0]); ++index) {
		buffers[index]->shader_read_index = TFX_IMPL_DX12_NO_DESCRIPTOR;
		buffers[index]->shader_write_index = TFX_IMPL_DX12_NO_DESCRIPTOR;
	}
	renderer->color_ramps_index = TFX_IMPL_DX12_NO_DESCRIPTOR;
	renderer->baked_color_ramps_index = TFX_IMPL_DX12_NO_DESCRIPTOR;
	renderer->pool_cold_index = TFX_IMPL_DX12_NO_DESCRIPTOR;
	renderer->pool_noise_index = TFX_IMPL_DX12_NO_DESCRIPTOR;
	renderer->upload_tick_index = TFX_IMPL_DX12_NO_DESCRIPTOR;
	renderer->upload_spawn_index = TFX_IMPL_DX12_NO_DESCRIPTOR;
	renderer->upload_location_index = TFX_IMPL_DX12_NO_DESCRIPTOR;

	//The view uniforms are a constant buffer view, which has to be 256 bytes
	if (!tfx_impl_dx12__create_buffer(renderer, 256, 0, 0, &renderer->view_uniforms)) {
		return false;
	}
	renderer->view_uniforms.shader_read_index = tfx_impl_dx12__allocate_descriptor(renderer);
	D3D12_CONSTANT_BUFFER_VIEW_DESC constant_view = {};
	constant_view.BufferLocation = renderer->view_uniforms.resource->GetGPUVirtualAddress();
	constant_view.SizeInBytes = 256;
	device->CreateConstantBufferView(&constant_view, tfx_impl_dx12__descriptor(renderer, renderer->view_uniforms.shader_read_index));
	return true;
}

void tfx_impl_dx12_BeginFrame(tfx_impl_dx12_t *renderer, tfxU32 frame_index) {
	renderer->frame_index = frame_index % renderer->frames_in_flight;
	tfx_impl_dx12__process_release_list(renderer, &renderer->release_lists[renderer->frame_index]);
	renderer->upload_rings[renderer->frame_index].used = 0;
}

bool tfx_impl_dx12_CreateAtlas(tfx_impl_dx12_t *renderer, tfxU32 layer_size, tfxU32 layer_count) {
	TFX_IMPL_DX12_RELEASE(renderer->shape_atlas);
	tfxU32 mip_count = 1;
	while ((layer_size >> mip_count) > 0) {
		mip_count++;
	}
	renderer->shape_atlas = tfx_impl_dx12__create_texture_array(renderer->device, layer_size, layer_size, layer_count, mip_count, DXGI_FORMAT_R8G8B8A8_UNORM);
	if (!renderer->shape_atlas) {
		return false;
	}
	renderer->shape_atlas_state = D3D12_RESOURCE_STATE_COPY_DEST;
	renderer->shape_atlas_layer_size = layer_size;
	renderer->shape_atlas_layer_count = layer_count;
	renderer->shape_atlas_mip_count = mip_count;
	//Slot 0, because tfx_impl_GetUV leaves the image index in the upper 16 bits of texture_array_index at 0
	tfx_impl_dx12__create_texture_view(renderer, renderer->shape_atlas, DXGI_FORMAT_R8G8B8A8_UNORM, layer_count, 0);
	return true;
}

bool tfx_impl_dx12_UploadAtlas(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, const tfx_impl_atlas_t *atlas) {
	if (!renderer->shape_atlas || atlas->layer_size != renderer->shape_atlas_layer_size
		|| atlas->first_layer + atlas->layer_count > renderer->shape_atlas_layer_count) {
		tfx_Print("An atlas of %u layers from layer %u, %u pixels square, doesn't fit the %u layer, %u pixel texture array",
			atlas->layer_count, atlas->first_layer, atlas->layer_size, renderer->shape_atlas_layer_count, renderer->shape_atlas_layer_size);
		return false;
	}
	size_t layer_bytes = (size_t)atlas->layer_size * atlas->layer_size * 4;
	unsigned char *mips = (unsigned char *)malloc(layer_bytes / 2);
	if (!mips) {
		return false;
	}
	tfx_impl_dx12__transition(command_list, renderer->shape_atlas, &renderer->shape_atlas_state, D3D12_RESOURCE_STATE_COPY_DEST);
	bool uploaded = true;
	for (tfxU32 layer = 0; layer != atlas->layer_count && uploaded; ++layer) {
		const unsigned char *source = atlas->pixels + layer * layer_bytes;
		tfxU32 size = atlas->layer_size;
		tfxU32 array_layer = atlas->first_layer + layer;
		uploaded = tfx_impl_dx12__upload_texture(renderer, command_list, renderer->shape_atlas, array_layer * renderer->shape_atlas_mip_count, DXGI_FORMAT_R8G8B8A8_UNORM, size, size, 4, source);
		unsigned char *destination = mips;
		for (tfxU32 mip_level = 1; mip_level < renderer->shape_atlas_mip_count && uploaded; ++mip_level) {
			tfx_impl_dx12__downsample(source, size, destination);
			size = size > 1 ? size / 2 : 1;
			uploaded = tfx_impl_dx12__upload_texture(renderer, command_list, renderer->shape_atlas, mip_level + array_layer * renderer->shape_atlas_mip_count, DXGI_FORMAT_R8G8B8A8_UNORM,
				size, size, 4, destination);
			source = destination;
			destination += (size_t)size * size * 4;
		}
	}
	free(mips);
	tfx_impl_dx12__transition(command_list, renderer->shape_atlas, &renderer->shape_atlas_state, D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE);
	return uploaded;
}

//A color ramp texture array and its view, replacing the old ones once the frames using them are done
static ID3D12Resource *tfx_impl_dx12__replace_color_ramps(tfx_impl_dx12_t *renderer, ID3D12Resource *old_texture, tfxU32 *view_index, tfx_bitmap_t *first_bitmap, tfxU32 layer_count) {
	tfx_impl_dx12__defer_release(renderer, old_texture);
	tfx_impl_dx12__defer_free_descriptor(renderer, *view_index);
	*view_index = TFX_IMPL_DX12_NO_DESCRIPTOR;
	DXGI_FORMAT format = tfx_impl_dx12__color_ramp_format(tfx_impl_dx12__bitmap_bytes_per_pixel(first_bitmap));
	ID3D12Resource *texture = tfx_impl_dx12__create_texture_array(renderer->device, (tfxU32)tfx_GetBitmapWidth(first_bitmap), (tfxU32)tfx_GetBitmapHeight(first_bitmap), layer_count, 1, format);
	if (texture) {
		*view_index = tfx_impl_dx12__allocate_descriptor(renderer);
		if (*view_index != TFX_IMPL_DX12_NO_DESCRIPTOR) {
			tfx_impl_dx12__create_texture_view(renderer, texture, format, layer_count, *view_index);
		}
	}
	return texture;
}

bool tfx_impl_dx12_UploadSpriteData(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, tfx_animation_manager animation_manager, tfx_gpu_shapes shapes) {
	tfxU32 bitmap_count = tfx_GetAnimationColorRampBitmapCount(animation_manager);
	if (bitmap_count) {
		renderer->baked_color_ramps = tfx_impl_dx12__replace_color_ramps(renderer, renderer->baked_color_ramps, &renderer->baked_color_ramps_index,
			tfx_GetAnimationColorRampBitmap(animation_manager, 0), bitmap_count);
		if (!renderer->baked_color_ramps) {
			return false;
		}
		renderer->baked_color_ramps_state = D3D12_RESOURCE_STATE_COPY_DEST;
		for (tfxU32 index = 0; index != bitmap_count; ++index) {
			if (!tfx_impl_dx12__upload_color_ramp_layer(renderer, command_list, renderer->baked_color_ramps, index, tfx_GetAnimationColorRampBitmap(animation_manager, index))) {
				return false;
			}
		}
		tfx_impl_dx12__transition(command_list, renderer->baked_color_ramps, &renderer->baked_color_ramps_state, D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE);
	}

	tfxU32 shape_count = tfx_GetGPUShapesCount(shapes);
	size_t property_bytes = tfx_GetAnimationEmitterPropertySizeInBytes(animation_manager);
	size_t sprite_bytes = tfx_GetSpriteDataSizeInBytes(animation_manager);
	if (!tfx_impl_dx12__write_buffer(renderer, command_list, &renderer->baked_image_data, tfx_GetGPUShapesArray(shapes), sizeof(tfx_gpu_image_data_t) * shape_count, TFX_IMPL_DX12_IMAGE_DATA_STRIDE)
		|| !tfx_impl_dx12__write_buffer(renderer, command_list, &renderer->baked_emitter_properties, tfx_GetAnimationEmitterPropertiesBufferPointer(animation_manager), property_bytes, TFX_IMPL_DX12_PROPERTIES_STRIDE)
		|| !tfx_impl_dx12__write_buffer(renderer, command_list, &renderer->baked_sprites, tfx_GetSpriteDataBufferPointer(animation_manager), sprite_bytes, TFX_IMPL_DX12_INSTANCE_STRIDE)) {
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
	bool uploaded = tfx_impl_dx12__write_buffer(renderer, command_list, &renderer->baked_ribbons.graph_lookups, tfx_GetAnimationRibbonGraphLookupsBufferPointer(animation_manager),
			tfx_GetAnimationRibbonGraphLookupsSizeInBytes(animation_manager), TFX_IMPL_DX12_GRAPH_STRIDE)
		&& tfx_impl_dx12__write_buffer(renderer, command_list, &renderer->baked_ribbon_data, tfx_GetAnimationRibbonDataBufferPointer(animation_manager),
			tfx_GetAnimationRibbonDataSizeInBytes(animation_manager), TFX_IMPL_DX12_BAKED_RIBBON_STRIDE)
		&& tfx_impl_dx12__write_buffer(renderer, command_list, &renderer->baked_ribbon_properties, tfx_GetAnimationRibbonPropertiesBufferPointer(animation_manager),
			tfx_GetAnimationRibbonPropertySizeInBytes(animation_manager), TFX_IMPL_DX12_BAKED_RIBBON_PROPERTIES_STRIDE)
		&& tfx_impl_dx12__write_buffer(renderer, command_list, &renderer->baked_ribbons.segments, tfx_GetAnimationRibbonSegmentsBufferPointer(animation_manager),
			tfx_GetAnimationRibbonSegmentsSizeInBytes(animation_manager), TFX_IMPL_DX12_RIBBON_SEGMENT_STRIDE)
		&& tfx_impl_dx12__write_buffer(renderer, command_list, &renderer->baked_ribbons.emitters, emitters, emitters_bytes, TFX_IMPL_DX12_RIBBON_EMITTER_STRIDE);
	free(emitters);
	return uploaded;
}

//The vertex and index buffers ribbon_compute writes, sized for this frame's buckets
static bool tfx_impl_dx12__ensure_ribbon_geometry(tfx_impl_dx12_t *renderer, const tfx_ribbon_batch_sizes_t *sizes, tfx_impl_dx12_ribbon_set_t *set) {
	return tfx_impl_dx12__ensure_buffer(renderer, sizes->vertex_buffer_size_in_bytes, 0, TFX_IMPL_DX12_RIBBON_VERTEX_STRIDE, &set->vertices)
		&& tfx_impl_dx12__ensure_buffer(renderer, sizes->index_buffer_size_in_bytes, 0, TFX_IMPL_DX12_UINT_STRIDE, &set->indices);
}

static bool tfx_impl_dx12__upload_baked_ribbons(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, const tfx_impl_baked_frame_t *baked) {
	if (baked->ribbon_count == 0 || !renderer->baked_ribbon_data.resource) {
		return true;
	}
	tfx_ribbon_batch_sizes_t sizes = tfx_GetAnimationRibbonBufferSizes(baked->animation_manager, sizeof(tfx_ribbon_vertex_t));
	//ribbon_playback writes the ribbons and ribbon_compute reads them, so they need both views
	if (!tfx_impl_dx12__write_buffer(renderer, command_list, &renderer->animation_ribbon_offsets, baked->ribbon_offsets, baked->ribbon_offsets_bytes, TFX_IMPL_DX12_UINT_STRIDE)
		|| !tfx_impl_dx12__ensure_buffer(renderer, (size_t)baked->ribbon_count * tfxRIBBON_SIZE, TFX_IMPL_DX12_RIBBON_STRIDE, TFX_IMPL_DX12_RIBBON_STRIDE, &renderer->baked_ribbons.instances)
		|| !tfx_impl_dx12__ensure_ribbon_geometry(renderer, &sizes, &renderer->baked_ribbons)) {
		return false;
	}
	tfx_ribbon_data_push_t *push = &renderer->ribbon_data_push;
	memset(push, 0, sizeof(tfx_ribbon_data_push_t));
	push->animation_instances_total = baked->animation_instance_count;
	push->ribbons_total = baked->ribbon_count;
	push->offsets_index = renderer->animation_ribbon_offsets.shader_read_index;
	push->animation_instances_index = renderer->animation_instances.shader_read_index;
	push->ribbon_data_index = renderer->baked_ribbon_data.shader_read_index;
	push->ribbon_output_index = renderer->baked_ribbons.instances.shader_write_index;
	push->ribbon_properties_index = renderer->baked_ribbon_properties.shader_read_index;
	renderer->playback_ribbon_count = baked->ribbon_count;
	renderer->baked_ribbons.animation_manager = baked->animation_manager;
	return true;
}

//Animation instances advance every frame, so this runs every frame, not just on updates
static bool tfx_impl_dx12__upload_baked(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, const tfx_impl_baked_frame_t *baked) {
	renderer->playback_sprite_count = 0;
	renderer->playback_ribbon_count = 0;
	renderer->baked_ribbons.animation_manager = NULL;
	if (!renderer->baked_sprites.resource || (baked->sprite_count == 0 && baked->ribbon_count == 0)) {
		return true;
	}
	if (!tfx_impl_dx12__write_buffer(renderer, command_list, &renderer->animation_instances, baked->animation_instances, baked->animation_instances_bytes, TFX_IMPL_DX12_ANIMATION_INSTANCE_STRIDE)) {
		return false;
	}
	if (baked->sprite_count > 0) {
		if (!tfx_impl_dx12__write_buffer(renderer, command_list, &renderer->animation_offsets, baked->offsets, baked->offsets_bytes, TFX_IMPL_DX12_UINT_STRIDE)
			|| !tfx_impl_dx12__ensure_buffer(renderer, (size_t)baked->sprite_count * sizeof(tfx_instance_t), 0, TFX_IMPL_DX12_INSTANCE_STRIDE, &renderer->playback_instances)) {
			return false;
		}
		tfx_sprite_data_push_t *push = &renderer->sprite_data_push;
		memset(push, 0, sizeof(tfx_sprite_data_push_t));
		push->animation_instances_total = baked->animation_instance_count;
		push->billboards_total = baked->sprite_count;
		push->animated_shapes = baked->animated_shapes ? 1 : 0;
		push->offsets_index = renderer->animation_offsets.shader_read_index;
		push->animation_instances_index = renderer->animation_instances.shader_read_index;
		push->billboards_index = renderer->playback_instances.shader_write_index;
		push->sprite_data_index = renderer->baked_sprites.shader_read_index;
		push->image_data_index = renderer->baked_image_data.shader_read_index;
		push->emitter_properties_index = renderer->baked_emitter_properties.shader_read_index;
		renderer->playback_sprite_count = baked->sprite_count;
	}
	return tfx_impl_dx12__upload_baked_ribbons(renderer, command_list, baked);
}

static bool tfx_impl_dx12__upload_ribbons(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, const tfx_impl_ribbon_frame_t *ribbons) {
	tfx_impl_dx12_ribbon_set_t *set = &renderer->live_ribbons;
	set->stage = ribbons->stage;
	if (!ribbons->stage) {
		return true;
	}
	const tfx_ribbon_batch_sizes_t *sizes = &ribbons->sizes;
	if (!tfx_impl_dx12__write_buffer(renderer, command_list, &set->segments, ribbons->segments, sizes->segment_buffer_size_in_bytes, TFX_IMPL_DX12_RIBBON_SEGMENT_STRIDE)
		|| !tfx_impl_dx12__write_buffer(renderer, command_list, &set->instances, ribbons->ribbons, sizes->ribbon_buffer_size_in_bytes, TFX_IMPL_DX12_RIBBON_STRIDE)
		|| !tfx_impl_dx12__write_buffer(renderer, command_list, &set->emitters, ribbons->emitters, sizes->emitter_buffer_size_in_bytes, TFX_IMPL_DX12_RIBBON_EMITTER_STRIDE)
		|| !tfx_impl_dx12__ensure_ribbon_geometry(renderer, sizes, set)) {
		set->stage = NULL;
		return false;
	}
	set->stage_offsets = ribbons->stage_offsets;
	return true;
}

//The old contents are copied across, because the previous buffer is still what captured_index points into
static bool tfx_impl_dx12__grow_instance_buffers(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, tfxU32 needed) {
	tfxU32 capacity = renderer->instance_capacity ? renderer->instance_capacity : TFX_IMPL_DX12_INITIAL_INSTANCE_CAPACITY;
	while (capacity < needed) {
		capacity += capacity / 2;
	}
	for (int index = 0; index != 2; ++index) {
		tfx_impl_dx12_buffer_t grown = {};
		grown.shader_read_index = TFX_IMPL_DX12_NO_DESCRIPTOR;
		grown.shader_write_index = TFX_IMPL_DX12_NO_DESCRIPTOR;
		//The GPU particle kernel writes its instances after the CPU ones
		if (!tfx_impl_dx12__create_buffer(renderer, capacity * TFX_IMPL_DX12_INSTANCE_STRIDE, TFX_IMPL_DX12_INSTANCE_STRIDE, TFX_IMPL_DX12_INSTANCE_STRIDE, &grown)) {
			return false;
		}
		tfx_impl_dx12_buffer_t *old = &renderer->instance_buffers[index];
		if (old->resource) {
			tfx_impl_dx12__transition_buffer(command_list, old, D3D12_RESOURCE_STATE_COPY_SOURCE);
			tfx_impl_dx12__transition_buffer(command_list, &grown, D3D12_RESOURCE_STATE_COPY_DEST);
			command_list->CopyBufferRegion(grown.resource, 0, old->resource, 0, old->capacity);
		}
		tfx_impl_dx12__release_buffer(renderer, old);
		*old = grown;
	}
	renderer->instance_capacity = capacity;
	return true;
}

static void tfx_impl_dx12__replace_view(tfx_impl_dx12_t *renderer, tfxU32 *index, tfxU32 new_index) {
	tfx_impl_dx12__defer_free_descriptor(renderer, *index);
	*index = new_index;
}

static bool tfx_impl_dx12__upload_gpu_particles(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, const tfx_impl_frame_t *frame) {
	const tfx_gpu_particle_frame_t *gpu_frame = &frame->gpu_particles;
	const tfx_gpu_particle_buffer_sizes_t *pool_sizes = &frame->gpu_pool_sizes;
	//The pool holds the particles between frames, it's sized once and starts zeroed. Each section's offset is a multiple of its stride
	if (!renderer->particle_pool.resource) {
		tfxU32 pool_bytes = (tfxU32)pool_sizes->pool_bytes;
		if (!tfx_impl_dx12__create_buffer(renderer, pool_bytes, 0, TFX_IMPL_DX12_HOT_STRIDE, &renderer->particle_pool)) {
			return false;
		}
		renderer->pool_cold_index = tfx_impl_dx12__create_write_view(renderer, renderer->particle_pool.resource, renderer->particle_pool.capacity, TFX_IMPL_DX12_COLD_STRIDE);
		renderer->pool_noise_index = tfx_impl_dx12__create_write_view(renderer, renderer->particle_pool.resource, renderer->particle_pool.capacity, TFX_IMPL_DX12_NOISE_STRIDE);
		ID3D12Resource *source = NULL;
		UINT64 offset = 0;
		unsigned char *zeroes = tfx_impl_dx12__allocate_upload(renderer, pool_bytes, 16, &source, &offset);
		if (!zeroes) {
			return false;
		}
		memset(zeroes, 0, pool_bytes);
		tfx_impl_dx12__transition_buffer(command_list, &renderer->particle_pool, D3D12_RESOURCE_STATE_COPY_DEST);
		command_list->CopyBufferRegion(renderer->particle_pool.resource, 0, source, offset, pool_bytes);
	}
	//A view per row stride over the whole upload, each section starting on a multiple of its stride
	if (!renderer->particle_upload.resource || renderer->particle_upload.capacity < gpu_frame->upload_bytes) {
		if (!tfx_impl_dx12__ensure_buffer(renderer, gpu_frame->upload_bytes, TFX_IMPL_DX12_DISPATCH_STRIDE, 0, &renderer->particle_upload)) {
			return false;
		}
		ID3D12Resource *upload = renderer->particle_upload.resource;
		tfxU32 capacity = renderer->particle_upload.capacity;
		tfx_impl_dx12__replace_view(renderer, &renderer->upload_tick_index, tfx_impl_dx12__create_read_view(renderer, upload, capacity, TFX_IMPL_DX12_TICK_STRIDE));
		tfx_impl_dx12__replace_view(renderer, &renderer->upload_spawn_index, tfx_impl_dx12__create_read_view(renderer, upload, capacity, TFX_IMPL_DX12_SPAWN_STRIDE));
		tfx_impl_dx12__replace_view(renderer, &renderer->upload_location_index, tfx_impl_dx12__create_read_view(renderer, upload, capacity, TFX_IMPL_DX12_LOCATION_STRIDE));
	}
	if (gpu_frame->upload_bytes) {
		ID3D12Resource *source = NULL;
		UINT64 offset = 0;
		unsigned char *mapped = tfx_impl_dx12__allocate_upload(renderer, gpu_frame->upload_bytes, 16, &source, &offset);
		if (!mapped) {
			return false;
		}
		tfx_CopyGPUParticleUpload(gpu_frame, mapped);
		tfx_impl_dx12__transition_buffer(command_list, &renderer->particle_upload, D3D12_RESOURCE_STATE_COPY_DEST);
		command_list->CopyBufferRegion(renderer->particle_upload.resource, 0, source, offset, gpu_frame->upload_bytes);
	}

	//The bases are in rows of each section's own stride, which is what the views count in
	tfx_gpu_particle_push_t *push = &renderer->particle_push_constants;
	memset(push, 0, sizeof(tfx_gpu_particle_push_t));
	push->hot_buffer_index = renderer->particle_pool.shader_write_index;
	push->cold_buffer_index = renderer->pool_cold_index;
	push->noise_buffer_index = renderer->pool_noise_index;
	push->graph_buffer_index = renderer->graph_lookups.shader_read_index;
	push->dispatch_buffer_index = renderer->particle_upload.shader_read_index;
	push->tick_buffer_index = renderer->upload_tick_index;
	push->spawn_buffer_index = renderer->upload_spawn_index;
	push->location_buffer_index = renderer->upload_location_index;
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
static bool tfx_impl_dx12__upload_color_ramps(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list) {
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
		capacity = capacity < TFX_IMPL_DX12_MAX_COLOR_RAMP_LAYERS ? capacity : TFX_IMPL_DX12_MAX_COLOR_RAMP_LAYERS;
		tfxU32 *layer_versions = (tfxU32 *)realloc(renderer->color_ramp_layer_versions, sizeof(tfxU32) * capacity);
		if (!layer_versions) {
			return false;
		}
		renderer->color_ramp_layer_versions = layer_versions;
		renderer->color_ramp_layer_capacity = 0;
		renderer->color_ramps = tfx_impl_dx12__replace_color_ramps(renderer, renderer->color_ramps, &renderer->color_ramps_index, first_bitmap, capacity);
		if (!renderer->color_ramps) {
			return false;
		}
		renderer->color_ramps_state = D3D12_RESOURCE_STATE_COPY_DEST;
		renderer->color_ramp_layer_capacity = capacity;
	}
	tfx_impl_dx12__transition(command_list, renderer->color_ramps, &renderer->color_ramps_state, D3D12_RESOURCE_STATE_COPY_DEST);
	for (tfxU32 layer = 0; layer != layer_count; ++layer) {
		tfxU32 version = tfx_GetGlobalColorRampLayerVersion(layer);
		if (!recreate && version == renderer->color_ramp_layer_versions[layer]) {
			continue;
		}
		renderer->color_ramp_layer_versions[layer] = version;
		//A free layer isn't drawn from, so whatever it holds can stay
		tfx_bitmap_t *bitmap = tfx_GetGlobalColorRampLayer(layer);
		if (bitmap && !tfx_impl_dx12__upload_color_ramp_layer(renderer, command_list, renderer->color_ramps, layer, bitmap)) {
			return false;
		}
	}
	tfx_impl_dx12__transition(command_list, renderer->color_ramps, &renderer->color_ramps_state, D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE);
	return true;
}

//Only on update frames: the instances drawn until then were written with the layout the last upload holds
static bool tfx_impl_dx12__upload_library_data(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list) {
	tfxU32 properties_version = tfx_GetGlobalParticlePropertiesVersion();
	if (!renderer->emitter_properties.resource || properties_version != renderer->properties_version) {
		if (!tfx_impl_dx12__write_buffer(renderer, command_list, &renderer->emitter_properties, tfx_GetGlobalParticlePropertiesBuffer(), tfx_GetGlobalParticlePropertiesBufferSizeInBytes(),
			TFX_IMPL_DX12_PROPERTIES_STRIDE)) {
			return false;
		}
		renderer->properties_version = properties_version;
	}
	tfxU32 shapes_version = tfx_GetGlobalGPUShapesVersion();
	if (!renderer->image_data.resource || shapes_version != renderer->shapes_version) {
		tfx_gpu_shapes shapes = tfx_GetGlobalGPUShapes();
		if (!tfx_impl_dx12__write_buffer(renderer, command_list, &renderer->image_data, tfx_GetGPUShapesArray(shapes), sizeof(tfx_gpu_image_data_t) * tfx_GetGPUShapesCount(shapes),
			TFX_IMPL_DX12_IMAGE_DATA_STRIDE)) {
			return false;
		}
		renderer->shapes_version = shapes_version;
	}
	//Shared by every library already, read by the ribbon shaders and the GPU kernel
	tfxU32 graph_lookups_version = tfx_GetGPUGraphLookupsVersion();
	if (!renderer->graph_lookups.resource || graph_lookups_version != renderer->graph_lookups_version) {
		if (!tfx_impl_dx12__write_buffer(renderer, command_list, &renderer->graph_lookups, tfx_GetGPUGraphLookupsBuffer(), tfx_GetGPUGraphLookupsBufferSizeInBytes(), TFX_IMPL_DX12_GRAPH_STRIDE)) {
			return false;
		}
		renderer->graph_lookups_version = graph_lookups_version;
	}
	tfxU32 color_ramps_version = tfx_GetGlobalColorRampsVersion();
	if (!renderer->color_ramps || color_ramps_version != renderer->color_ramps_version) {
		if (!tfx_impl_dx12__upload_color_ramps(renderer, command_list)) {
			return false;
		}
		renderer->color_ramps_version = color_ramps_version;
	}
	return true;
}

bool tfx_impl_dx12_UploadFrame(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, const tfx_impl_frame_t *frame) {
	memcpy(renderer->layer_flags, frame->layer_flags, sizeof(renderer->layer_flags));
	if (!tfx_impl_dx12__upload_baked(renderer, command_list, &frame->baked)) {
		return false;
	}
	if (!frame->stage_updated) {
		return true;
	}
	if (!tfx_impl_dx12__upload_library_data(renderer, command_list)) {
		return false;
	}
	bool gpu_particles = frame->has_gpu_particles && renderer->particle_kernel;
	tfxU32 needed = gpu_particles ? frame->gpu_particles.gpu_instance_start + frame->gpu_particles.gpu_instance_count : frame->instance_count;
	if ((needed > renderer->instance_capacity || renderer->instance_capacity == 0) && !tfx_impl_dx12__grow_instance_buffers(renderer, command_list, needed)) {
		renderer->instance_count = 0;
		return false;
	}

	//The other buffer now holds what was drawn last, which is what captured_index points into for CPU and GPU instances alike
	renderer->current_instance_buffer ^= 1;
	if (!tfx_impl_dx12__upload_buffer(renderer, command_list, &renderer->instance_buffers[renderer->current_instance_buffer], frame->instances, frame->instance_count * sizeof(tfx_instance_t))) {
		return false;
	}
	renderer->instance_count = needed;
	memcpy(renderer->layer_start, frame->layer_start, sizeof(renderer->layer_start));
	memcpy(renderer->layer_count, frame->layer_count, sizeof(renderer->layer_count));
	memset(renderer->gpu_layer_count, 0, sizeof(renderer->gpu_layer_count));
	if (gpu_particles) {
		memcpy(renderer->gpu_layer_start, frame->gpu_particles.gpu_layer_start, sizeof(renderer->gpu_layer_start));
		memcpy(renderer->gpu_layer_count, frame->gpu_particles.gpu_layer_count, sizeof(renderer->gpu_layer_count));
		if (!tfx_impl_dx12__upload_gpu_particles(renderer, command_list, frame)) {
			return false;
		}
	}
	return tfx_impl_dx12__upload_ribbons(renderer, command_list, &frame->ribbons);
}

static void tfx_impl_dx12__bind_heaps(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list) {
	ID3D12DescriptorHeap *heaps[2] = { renderer->descriptor_heap, renderer->sampler_heap };
	command_list->SetDescriptorHeaps(2, heaps);
}

static void tfx_impl_dx12__bind_compute(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, ID3D12RootSignature *root_signature, ID3D12PipelineState *pipeline) {
	command_list->SetComputeRootSignature(root_signature);
	command_list->SetPipelineState(pipeline);
	command_list->SetComputeRootDescriptorTable(TFX_IMPL_DX12_DESCRIPTOR_TABLE_PARAMETER, renderer->descriptor_heap->GetGPUDescriptorHandleForHeapStart());
	command_list->SetComputeRootDescriptorTable(TFX_IMPL_DX12_SAMPLER_TABLE_PARAMETER, renderer->sampler_heap->GetGPUDescriptorHandleForHeapStart());
}

static void tfx_impl_dx12__bind_graphics(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, ID3D12RootSignature *root_signature, ID3D12PipelineState *pipeline) {
	command_list->SetGraphicsRootSignature(root_signature);
	command_list->SetPipelineState(pipeline);
	command_list->SetGraphicsRootDescriptorTable(TFX_IMPL_DX12_DESCRIPTOR_TABLE_PARAMETER, renderer->descriptor_heap->GetGPUDescriptorHandleForHeapStart());
	command_list->SetGraphicsRootDescriptorTable(TFX_IMPL_DX12_SAMPLER_TABLE_PARAMETER, renderer->sampler_heap->GetGPUDescriptorHandleForHeapStart());
}

//Every fetch has to be dispatched, the next one carries on from the particle state this one leaves in the pool. The dispatches
//write separate rings and instance ranges, so they need no barrier between them, only one against the last frame's
static void tfx_impl_dx12__simulate_gpu_particles(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list) {
	tfx_impl_dx12_buffer_t *instances = &renderer->instance_buffers[renderer->current_instance_buffer];
	tfx_impl_dx12__begin_writes(command_list, &renderer->particle_pool);
	tfx_impl_dx12__begin_writes(command_list, instances);
	tfx_impl_dx12__transition_buffer(command_list, &renderer->particle_upload, TFX_IMPL_DX12_SHADER_READ);
	tfx_impl_dx12__transition_buffer(command_list, &renderer->graph_lookups, TFX_IMPL_DX12_SHADER_READ);
	tfx_impl_dx12__bind_compute(renderer, command_list, renderer->particle_root_signature, renderer->particle_kernel);
	tfx_gpu_particle_push_t *push = &renderer->particle_push_constants;
	push->instance_buffer_index = instances->shader_write_index;
	command_list->SetComputeRoot32BitConstants(TFX_IMPL_DX12_PUSH_PARAMETER, sizeof(tfx_gpu_particle_push_t) / 4, push, 0);
	for (tfxU32 index = 0; index != renderer->particle_dispatch_count; ++index) {
		tfxU32 slot_count = renderer->particle_dispatches[index].slot_count;
		if (slot_count == 0) {
			continue;
		}
		//dispatch_index is the first push constant, so it's the only one that changes
		command_list->SetComputeRoot32BitConstant(TFX_IMPL_DX12_PUSH_PARAMETER, index, 0);
		command_list->Dispatch((slot_count + TFX_IMPL_DX12_KERNEL_GROUP_SIZE - 1) / TFX_IMPL_DX12_KERNEL_GROUP_SIZE, 1, 1);
	}
	renderer->particle_dispatch_count = 0;
}

static void tfx_impl_dx12__play_back_sprites(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list) {
	tfx_impl_dx12__transition_buffer(command_list, &renderer->animation_offsets, TFX_IMPL_DX12_SHADER_READ);
	tfx_impl_dx12__transition_buffer(command_list, &renderer->animation_instances, TFX_IMPL_DX12_SHADER_READ);
	tfx_impl_dx12__transition_buffer(command_list, &renderer->baked_sprites, TFX_IMPL_DX12_SHADER_READ);
	tfx_impl_dx12__transition_buffer(command_list, &renderer->baked_emitter_properties, TFX_IMPL_DX12_SHADER_READ);
	tfx_impl_dx12__begin_writes(command_list, &renderer->playback_instances);
	tfx_impl_dx12__bind_compute(renderer, command_list, renderer->sprite_playback_root_signature, renderer->sprite_playback_pipeline);
	command_list->SetComputeRoot32BitConstants(TFX_IMPL_DX12_PUSH_PARAMETER, sizeof(tfx_sprite_data_push_t) / 4, &renderer->sprite_data_push, 0);
	command_list->Dispatch((renderer->playback_sprite_count + TFX_IMPL_DX12_PLAYBACK_GROUP_SIZE - 1) / TFX_IMPL_DX12_PLAYBACK_GROUP_SIZE, 1, 1);
}

static void tfx_impl_dx12__play_back_ribbons(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list) {
	tfx_impl_dx12__transition_buffer(command_list, &renderer->animation_ribbon_offsets, TFX_IMPL_DX12_SHADER_READ);
	tfx_impl_dx12__transition_buffer(command_list, &renderer->animation_instances, TFX_IMPL_DX12_SHADER_READ);
	tfx_impl_dx12__transition_buffer(command_list, &renderer->baked_ribbon_data, TFX_IMPL_DX12_SHADER_READ);
	tfx_impl_dx12__transition_buffer(command_list, &renderer->baked_ribbon_properties, TFX_IMPL_DX12_SHADER_READ);
	tfx_impl_dx12__begin_writes(command_list, &renderer->baked_ribbons.instances);
	tfx_impl_dx12__bind_compute(renderer, command_list, renderer->ribbon_playback_root_signature, renderer->ribbon_playback_pipeline);
	command_list->SetComputeRoot32BitConstants(TFX_IMPL_DX12_PUSH_PARAMETER, sizeof(tfx_ribbon_data_push_t) / 4, &renderer->ribbon_data_push, 0);
	command_list->Dispatch((renderer->playback_ribbon_count + TFX_IMPL_DX12_PLAYBACK_GROUP_SIZE - 1) / TFX_IMPL_DX12_PLAYBACK_GROUP_SIZE, 1, 1);
}

static bool tfx_impl_dx12__ribbon_set_active(const tfx_impl_dx12_ribbon_set_t *set) {
	return set->stage || set->animation_manager;
}

static tfx_ribbon_dispatch_t tfx_impl_dx12__create_ribbon_dispatch(const tfx_impl_dx12_ribbon_set_t *set) {
	return set->stage ? tfx_CreateRibbonBatchDispatch(&set->stage_offsets) : tfx_CreateAnimationRibbonDispatch();
}

static bool tfx_impl_dx12__next_ribbon_dispatch(const tfx_impl_dx12_ribbon_set_t *set, tfx_ribbon_dispatch_t *dispatch) {
	return set->stage ? tfx_NextRibbonDispatch(set->stage, dispatch) : tfx_NextAnimationRibbonDispatch(set->animation_manager, dispatch);
}

static tfx_impl_dx12_buffer_t *tfx_impl_dx12__ribbon_graph_lookups(tfx_impl_dx12_t *renderer, tfx_impl_dx12_ribbon_set_t *set) {
	return set->stage ? &renderer->graph_lookups : &set->graph_lookups;
}

//The library fills a bucket's counts and offsets, the descriptor indexes are the host's
static void tfx_impl_dx12__set_ribbon_indexes(tfx_impl_dx12_t *renderer, tfx_impl_dx12_ribbon_set_t *set, tfx_ribbon_bucket_globals_t *globals, const tfx_impl_dx12_buffer_t *image_data, tfxU32 color_ramps_index) {
	globals->uniform_index = renderer->view_uniforms.shader_read_index;
	globals->emitters_index = set->emitters.shader_read_index;
	globals->graphs_index = tfx_impl_dx12__ribbon_graph_lookups(renderer, set)->shader_read_index;
	globals->ribbons_index = set->instances.shader_read_index;
	globals->ribbon_segments_index = set->segments.shader_read_index;
	globals->vertexes_index = set->vertices.shader_write_index;
	globals->indexes_index = set->indices.shader_write_index;
	globals->image_data_index = image_data->shader_read_index;
	globals->sampler_index = 0;
	globals->particle_texture_index = 0;
	globals->color_ramp_texture_index = color_ramps_index;
}

static void tfx_impl_dx12__build_ribbon_geometry(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, tfx_impl_dx12_ribbon_set_t *set, const tfx_impl_uniforms_t *uniforms) {
	tfx_impl_dx12__transition_buffer(command_list, &set->emitters, TFX_IMPL_DX12_SHADER_READ);
	tfx_impl_dx12__transition_buffer(command_list, tfx_impl_dx12__ribbon_graph_lookups(renderer, set), TFX_IMPL_DX12_SHADER_READ);
	tfx_impl_dx12__transition_buffer(command_list, &set->instances, TFX_IMPL_DX12_SHADER_READ);
	tfx_impl_dx12__transition_buffer(command_list, &set->segments, TFX_IMPL_DX12_SHADER_READ);
	tfx_impl_dx12__begin_writes(command_list, &set->vertices);
	tfx_impl_dx12__begin_writes(command_list, &set->indices);
	tfx_impl_dx12__bind_compute(renderer, command_list, renderer->ribbon_root_signature, renderer->ribbon_compute_pipeline);
	const tfx_impl_dx12_buffer_t *image_data = set->stage ? &renderer->image_data : &renderer->baked_image_data;
	tfxU32 color_ramps_index = set->stage ? renderer->color_ramps_index : renderer->baked_color_ramps_index;

	//One dispatch per bucket, each with the bucket's own globals
	tfx_ribbon_dispatch_t dispatch = tfx_impl_dx12__create_ribbon_dispatch(set);
	while (tfx_impl_dx12__next_ribbon_dispatch(set, &dispatch)) {
		tfx_ribbon_bucket_globals_t *globals = tfx_GetRibbonDispatchGlobals(&dispatch);
		if (set->stage) {
			globals->lerp = uniforms->timer_lerp;
		} else {
			//Baked ribbons are already interpolated, but they still face the camera
			globals->camera_position.x = uniforms->camera_position[0];
			globals->camera_position.y = uniforms->camera_position[1];
			globals->camera_position.z = uniforms->camera_position[2];
		}
		tfx_impl_dx12__set_ribbon_indexes(renderer, set, globals, image_data, color_ramps_index);
		command_list->SetComputeRoot32BitConstants(TFX_IMPL_DX12_PUSH_PARAMETER, sizeof(tfx_ribbon_bucket_globals_t) / 4, globals, 0);
		command_list->Dispatch((dispatch.total_segments + TFX_IMPL_DX12_RIBBON_GROUP_SIZE - 1) / TFX_IMPL_DX12_RIBBON_GROUP_SIZE, 1, 1);
	}
}

void tfx_impl_dx12_Compute(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, const tfx_impl_uniforms_t *uniforms) {
	tfx_impl_dx12__bind_heaps(renderer, command_list);
	if (renderer->particle_dispatch_count > 0) {
		tfx_impl_dx12__simulate_gpu_particles(renderer, command_list);
	}
	if (renderer->playback_sprite_count > 0) {
		tfx_impl_dx12__play_back_sprites(renderer, command_list);
	}
	if (renderer->playback_ribbon_count > 0) {
		tfx_impl_dx12__play_back_ribbons(renderer, command_list);
	}
	if (tfx_impl_dx12__ribbon_set_active(&renderer->live_ribbons)) {
		tfx_impl_dx12__build_ribbon_geometry(renderer, command_list, &renderer->live_ribbons, uniforms);
	}
	if (tfx_impl_dx12__ribbon_set_active(&renderer->baked_ribbons)) {
		tfx_impl_dx12__build_ribbon_geometry(renderer, command_list, &renderer->baked_ribbons, uniforms);
	}
}

//Live and baked ribbons each have their own shape frames and color ramps, the atlas is shared
//Every bucket is walked, because the dispatch carries the offsets from one to the next, and only the ones whose layer is in the group are drawn
static void tfx_impl_dx12__render_ribbons(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, tfx_impl_dx12_ribbon_set_t *set, const tfx_impl_dx12_buffer_t *image_data, tfxU32 color_ramps_index, int writes_depth) {
	D3D12_VERTEX_BUFFER_VIEW vertex_view = { set->vertices.resource->GetGPUVirtualAddress(), set->vertices.capacity, TFX_IMPL_DX12_RIBBON_VERTEX_STRIDE };
	D3D12_INDEX_BUFFER_VIEW index_view = { set->indices.resource->GetGPUVirtualAddress(), set->indices.capacity, DXGI_FORMAT_R32_UINT };
	command_list->IASetVertexBuffers(0, 1, &vertex_view);
	command_list->IASetIndexBuffer(&index_view);
	ID3D12PipelineState *bound_pipeline = NULL;
	tfx_ribbon_dispatch_t dispatch = tfx_impl_dx12__create_ribbon_dispatch(set);
	while (tfx_impl_dx12__next_ribbon_dispatch(set, &dispatch)) {
		const tfxLayerFlags flags = renderer->layer_flags[dispatch.layer];
		if (tfx_impl_LayerWritesDepth(flags) != writes_depth) {
			continue;
		}
		ID3D12PipelineState *pipeline = renderer->ribbon_pipelines[writes_depth][tfx_impl_LayerIsLit(flags)];
		if (pipeline != bound_pipeline) {
			tfx_impl_dx12__bind_graphics(renderer, command_list, renderer->ribbon_root_signature, pipeline);
			bound_pipeline = pipeline;
		}
		tfx_ribbon_bucket_globals_t *globals = tfx_GetRibbonDispatchGlobals(&dispatch);
		tfx_impl_dx12__set_ribbon_indexes(renderer, set, globals, image_data, color_ramps_index);
		command_list->SetGraphicsRoot32BitConstants(TFX_IMPL_DX12_PUSH_PARAMETER, sizeof(tfx_ribbon_bucket_globals_t) / 4, globals, 0);
		command_list->DrawIndexedInstanced(dispatch.index_count, 1, dispatch.index_offset, 0, 0);
	}
}

//Every buffer a draw reads goes to the read state before the first draw, which also waits for the compute writes
static void tfx_impl_dx12__prepare_draws(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list) {
	tfx_impl_dx12_buffer_t *reads[] = { &renderer->view_uniforms, &renderer->instance_buffers[0], &renderer->instance_buffers[1], &renderer->image_data, &renderer->emitter_properties,
		&renderer->graph_lookups, &renderer->playback_instances, &renderer->baked_image_data, &renderer->baked_emitter_properties,
		&renderer->live_ribbons.vertices, &renderer->live_ribbons.indices, &renderer->live_ribbons.instances, &renderer->live_ribbons.emitters,
		&renderer->baked_ribbons.vertices, &renderer->baked_ribbons.indices, &renderer->baked_ribbons.instances, &renderer->baked_ribbons.emitters, &renderer->baked_ribbons.graph_lookups };
	for (tfxU32 index = 0; index != sizeof(reads) / sizeof(reads[0]); ++index) {
		tfx_impl_dx12__transition_buffer(command_list, reads[index], TFX_IMPL_DX12_SHADER_READ);
	}
	tfx_impl_dx12__transition(command_list, renderer->shape_atlas, &renderer->shape_atlas_state, D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE);
}

static void tfx_impl_dx12__draw_instance_range(ID3D12GraphicsCommandList *command_list, tfxU32 start, tfxU32 count) {
	//The input assembler starts the instance buffer at the first instance, so nothing extra goes in the push constants
	if (count) {
		command_list->DrawInstanced(6, count, 0, start);
	}
}

//Each layer in the group as its CPU instances then its GPU instances
static void tfx_impl_dx12__render_billboards(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, const tfx_impl_dx12_billboard_push_t *push, int writes_depth) {
	tfx_impl_dx12_buffer_t *current = &renderer->instance_buffers[renderer->current_instance_buffer];
	D3D12_VERTEX_BUFFER_VIEW vertex_view = { current->resource->GetGPUVirtualAddress(), current->capacity, TFX_IMPL_DX12_INSTANCE_STRIDE };
	command_list->IASetVertexBuffers(0, 1, &vertex_view);
	ID3D12PipelineState *bound_pipeline = NULL;
	for (tfxU32 layer = 0; layer != tfxLAYERS; ++layer) {
		const tfxLayerFlags flags = renderer->layer_flags[layer];
		if (tfx_impl_LayerWritesDepth(flags) != writes_depth || renderer->layer_count[layer] + renderer->gpu_layer_count[layer] == 0) {
			continue;
		}
		ID3D12PipelineState *pipeline = renderer->billboard_pipelines[writes_depth][tfx_impl_LayerIsLit(flags)];
		if (pipeline != bound_pipeline) {
			//Setting the root signature clears the root constants, so they go again with each pipeline
			tfx_impl_dx12__bind_graphics(renderer, command_list, renderer->billboard_root_signature, pipeline);
			command_list->SetGraphicsRoot32BitConstants(TFX_IMPL_DX12_PUSH_PARAMETER, sizeof(tfx_impl_dx12_billboard_push_t) / 4, push, 0);
			bound_pipeline = pipeline;
		}
		tfx_impl_dx12__draw_instance_range(command_list, renderer->layer_start[layer], renderer->layer_count[layer]);
		tfx_impl_dx12__draw_instance_range(command_list, renderer->gpu_layer_start[layer], renderer->gpu_layer_count[layer]);
	}
}

static void tfx_impl_dx12__render_baked_billboards(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, const tfx_impl_dx12_billboard_push_t *push) {
	//Baked sprites come with their own shape frames, emitter properties and color ramps, the atlas is shared. The fragment shader reads
	//the live push at b1 and the vertex shader the baked one at b2
	tfx_impl_dx12_baked_billboard_push_t baked_push = {};
	baked_push.billboard = *push;
	baked_push.billboard.color_ramp_texture_index = renderer->baked_color_ramps_index;
	baked_push.billboard.image_data_index = renderer->baked_image_data.shader_read_index;
	baked_push.billboard.properties_index = renderer->baked_emitter_properties.shader_read_index;
	baked_push.billboard.prev_billboards_index = 0;
	tfx_impl_dx12__bind_graphics(renderer, command_list, renderer->baked_billboard_root_signature, renderer->baked_billboard_pipeline);
	command_list->SetGraphicsRoot32BitConstants(TFX_IMPL_DX12_PUSH_PARAMETER, sizeof(tfx_impl_dx12_billboard_push_t) / 4, &baked_push.billboard, 0);
	command_list->SetGraphicsRoot32BitConstants(TFX_IMPL_DX12_PUSH_PARAMETER + 1, sizeof(baked_push) / 4, &baked_push, 0);
	D3D12_VERTEX_BUFFER_VIEW vertex_view = { renderer->playback_instances.resource->GetGPUVirtualAddress(), renderer->playback_instances.capacity, TFX_IMPL_DX12_INSTANCE_STRIDE };
	command_list->IASetVertexBuffers(0, 1, &vertex_view);
	tfx_impl_dx12__draw_instance_range(command_list, 0, renderer->playback_sprite_count);
}

void tfx_impl_dx12_Render(tfx_impl_dx12_t *renderer, ID3D12GraphicsCommandList *command_list, const tfx_impl_uniforms_t *uniforms) {
	if (!renderer->shape_atlas || (renderer->instance_count == 0 && renderer->playback_sprite_count == 0
		&& !tfx_impl_dx12__ribbon_set_active(&renderer->live_ribbons) && !tfx_impl_dx12__ribbon_set_active(&renderer->baked_ribbons))) {
		return;
	}
	if (!tfx_impl_dx12__upload_buffer(renderer, command_list, &renderer->view_uniforms, uniforms, sizeof(tfx_impl_uniforms_t))) {
		return;
	}
	tfx_impl_dx12__prepare_draws(renderer, command_list);
	tfx_impl_dx12__bind_heaps(renderer, command_list);
	command_list->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

	tfx_impl_dx12_billboard_push_t push = {};
	push.sampler_index = 0;
	push.uniform_index = renderer->view_uniforms.shader_read_index;
	push.color_ramp_texture_index = renderer->color_ramps_index;
	push.image_data_index = renderer->image_data.shader_read_index;
	push.properties_index = renderer->emitter_properties.shader_read_index;
	push.prev_billboards_index = renderer->instance_buffers[renderer->current_instance_buffer ^ 1].shader_read_index;
	//The layers that write depth draw first with blending off, so they hide the blended layers behind them
	for (int group = 0; group != 2; ++group) {
		const int writes_depth = group == 0;
		if (renderer->instance_count > 0) {
			tfx_impl_dx12__render_billboards(renderer, command_list, &push, writes_depth);
		}
		if (tfx_impl_dx12__ribbon_set_active(&renderer->live_ribbons)) {
			tfx_impl_dx12__render_ribbons(renderer, command_list, &renderer->live_ribbons, &renderer->image_data, renderer->color_ramps_index, writes_depth);
		}
		//Baked billboards come as one range across their layers, so they all draw blended and unlit
		if (!writes_depth && renderer->playback_sprite_count > 0) {
			tfx_impl_dx12__render_baked_billboards(renderer, command_list, &push);
		}
		if (tfx_impl_dx12__ribbon_set_active(&renderer->baked_ribbons)) {
			tfx_impl_dx12__render_ribbons(renderer, command_list, &renderer->baked_ribbons, &renderer->baked_image_data, renderer->baked_color_ramps_index, writes_depth);
		}
	}
}

static void tfx_impl_dx12__release_ribbon_set(tfx_impl_dx12_t *renderer, tfx_impl_dx12_ribbon_set_t *set) {
	tfx_impl_dx12__release_buffer(renderer, &set->segments);
	tfx_impl_dx12__release_buffer(renderer, &set->instances);
	tfx_impl_dx12__release_buffer(renderer, &set->emitters);
	tfx_impl_dx12__release_buffer(renderer, &set->vertices);
	tfx_impl_dx12__release_buffer(renderer, &set->indices);
	tfx_impl_dx12__release_buffer(renderer, &set->graph_lookups);
}

void tfx_impl_dx12_Shutdown(tfx_impl_dx12_t *renderer) {
	if (!renderer->device) {
		return;
	}
	tfx_impl_dx12_buffer_t *buffers[] = { &renderer->view_uniforms, &renderer->image_data, &renderer->emitter_properties, &renderer->instance_buffers[0], &renderer->instance_buffers[1],
		&renderer->particle_pool, &renderer->particle_upload, &renderer->baked_image_data, &renderer->baked_emitter_properties, &renderer->baked_sprites, &renderer->animation_offsets,
		&renderer->animation_instances, &renderer->playback_instances, &renderer->baked_ribbon_data, &renderer->baked_ribbon_properties, &renderer->animation_ribbon_offsets, &renderer->graph_lookups };
	for (tfxU32 index = 0; index != sizeof(buffers) / sizeof(buffers[0]); ++index) {
		tfx_impl_dx12__release_buffer(renderer, buffers[index]);
	}
	tfx_impl_dx12__release_ribbon_set(renderer, &renderer->live_ribbons);
	tfx_impl_dx12__release_ribbon_set(renderer, &renderer->baked_ribbons);
	tfx_impl_dx12__defer_release(renderer, renderer->shape_atlas);
	tfx_impl_dx12__defer_release(renderer, renderer->color_ramps);
	tfx_impl_dx12__defer_release(renderer, renderer->baked_color_ramps);
	renderer->shape_atlas = renderer->color_ramps = renderer->baked_color_ramps = NULL;
	//The GPU is idle, so every frame's leftovers can go now
	for (tfxU32 frame = 0; frame != TFX_IMPL_DX12_MAX_FRAMES_IN_FLIGHT; ++frame) {
		if (renderer->free_descriptors) {
			tfx_impl_dx12__process_release_list(renderer, &renderer->release_lists[frame]);
		}
		free(renderer->release_lists[frame].objects);
		free(renderer->release_lists[frame].descriptors);
		TFX_IMPL_DX12_RELEASE(renderer->upload_rings[frame].resource);
	}
	memset(renderer->release_lists, 0, sizeof(renderer->release_lists));
	free(renderer->color_ramp_layer_versions);
	free(renderer->free_descriptors);
	renderer->color_ramp_layer_versions = NULL;
	renderer->free_descriptors = NULL;
	TFX_IMPL_DX12_RELEASE(renderer->particle_kernel);
	for (int writes_depth = 0; writes_depth != 2; ++writes_depth) {
		for (int lit = 0; lit != 2; ++lit) {
			TFX_IMPL_DX12_RELEASE(renderer->billboard_pipelines[writes_depth][lit]);
			TFX_IMPL_DX12_RELEASE(renderer->ribbon_pipelines[writes_depth][lit]);
		}
	}
	TFX_IMPL_DX12_RELEASE(renderer->baked_billboard_pipeline);
	TFX_IMPL_DX12_RELEASE(renderer->ribbon_compute_pipeline);
	TFX_IMPL_DX12_RELEASE(renderer->sprite_playback_pipeline);
	TFX_IMPL_DX12_RELEASE(renderer->ribbon_playback_pipeline);
	TFX_IMPL_DX12_RELEASE(renderer->particle_root_signature);
	TFX_IMPL_DX12_RELEASE(renderer->billboard_root_signature);
	TFX_IMPL_DX12_RELEASE(renderer->baked_billboard_root_signature);
	TFX_IMPL_DX12_RELEASE(renderer->ribbon_root_signature);
	TFX_IMPL_DX12_RELEASE(renderer->sprite_playback_root_signature);
	TFX_IMPL_DX12_RELEASE(renderer->ribbon_playback_root_signature);
	TFX_IMPL_DX12_RELEASE(renderer->descriptor_heap);
	TFX_IMPL_DX12_RELEASE(renderer->sampler_heap);
	TFX_IMPL_DX12_RELEASE(renderer->device);
}
