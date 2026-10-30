#include "tfx_impl_dx11.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TFX_IMPL_DX11_INITIAL_INSTANCE_CAPACITY 10000
//Matches tfx_ribbon_group_size in timelinefx.slang
#define TFX_IMPL_DX11_RIBBON_GROUP_SIZE 256
#define TFX_IMPL_DX11_RELEASE(object) if (object) { (object)->Release(); (object) = NULL; }

//timelinefx.slang mirrors these structs by hand, in its [Render_rows] section
static_assert(sizeof(tfx_instance_t) == 48, "tfx_instance_t changed, update tfx_instance in timelinefx.slang");
static_assert(sizeof(tfx_gpu_image_data_t) == 48, "tfx_gpu_image_data_t changed, update tfx_image_data in timelinefx.slang");
static_assert(sizeof(tfx_gpu_particle_properties_t) == 32, "tfx_gpu_particle_properties_t changed, update tfx_emitter_properties in timelinefx.slang");
static_assert(sizeof(tfx_ribbon_vertex_t) == 32, "tfx_ribbon_vertex_t changed, update tfx_ribbon_vertex_row in timelinefx.slang");
static_assert(sizeof(tfx_gpu_graph_data_t) == 48, "tfx_gpu_graph_data_t changed, update tfx_gpu_graph_data in timelinefx.slang");
static_assert(sizeof(tfx_impl_uniforms_t) % 16 == 0, "Constant buffers must be a multiple of 16 bytes");

//Not in the public header, these sizes are what timelinefx.slang's structs are laid out for
#define TFX_IMPL_DX11_RIBBON_SEGMENT_SIZE 16
#define TFX_IMPL_DX11_RIBBON_SIZE 48
#define TFX_IMPL_DX11_RIBBON_EMITTER_SIZE 416
#define TFX_IMPL_DX11_RIBBON_DATA_SIZE 64
#define TFX_IMPL_DX11_RIBBON_PROPERTIES_SIZE 68

//Mirrors tfx_fixed_draw_constants in timelinefx.slang
typedef struct tfx_impl_dx11_draw_constants_s {
	tfxU32 instance_offset;
	tfxU32 previous_index_offset;
	float ndc_offset_x;
	float ndc_offset_y;
} tfx_impl_dx11_draw_constants_t;

//--- Dx11_resources

static void tfx_impl_dx11__release_buffer(tfx_impl_dx11_buffer_t *buffer) {
	TFX_IMPL_DX11_RELEASE(buffer->unordered_view);
	TFX_IMPL_DX11_RELEASE(buffer->view);
	TFX_IMPL_DX11_RELEASE(buffer->buffer);
	buffer->capacity = 0;
}

//A compute shader can write it when unordered is true
static bool tfx_impl_dx11__create_structured_buffer(ID3D11Device *device, const void *data, tfxU32 element_size, tfxU32 element_count, bool unordered, tfx_impl_dx11_buffer_t *buffer) {
	tfx_impl_dx11__release_buffer(buffer);
	D3D11_BUFFER_DESC description = {};
	description.ByteWidth = element_size * element_count;
	description.Usage = D3D11_USAGE_DEFAULT;
	description.BindFlags = D3D11_BIND_SHADER_RESOURCE | (unordered ? D3D11_BIND_UNORDERED_ACCESS : 0);
	description.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_STRUCTURED;
	description.StructureByteStride = element_size;
	D3D11_SUBRESOURCE_DATA initial_data = {};
	initial_data.pSysMem = data;
	HRESULT result = device->CreateBuffer(&description, data ? &initial_data : NULL, &buffer->buffer);
	if (FAILED(result)) {
		printf("Could not create a structured buffer of %u bytes (0x%08lx)\n", description.ByteWidth, result);
		return false;
	}
	D3D11_SHADER_RESOURCE_VIEW_DESC view_description = {};
	view_description.Format = DXGI_FORMAT_UNKNOWN;
	view_description.ViewDimension = D3D11_SRV_DIMENSION_BUFFER;
	view_description.Buffer.NumElements = element_count;
	result = device->CreateShaderResourceView(buffer->buffer, &view_description, &buffer->view);
	if (SUCCEEDED(result) && unordered) {
		D3D11_UNORDERED_ACCESS_VIEW_DESC unordered_description = {};
		unordered_description.Format = DXGI_FORMAT_UNKNOWN;
		unordered_description.ViewDimension = D3D11_UAV_DIMENSION_BUFFER;
		unordered_description.Buffer.NumElements = element_count;
		result = device->CreateUnorderedAccessView(buffer->buffer, &unordered_description, &buffer->unordered_view);
	}
	if (FAILED(result)) {
		printf("Could not create a structured buffer view (0x%08lx)\n", result);
		return false;
	}
	buffer->capacity = element_count;
	return true;
}

//Raw buffers hold rows in the C layout. bind_flags picks the views: shader resource, unordered access, and index buffer
static bool tfx_impl_dx11__create_raw_buffer(ID3D11Device *device, const void *data, tfxU32 size_in_bytes, UINT bind_flags, tfx_impl_dx11_buffer_t *buffer) {
	tfx_impl_dx11__release_buffer(buffer);
	size_in_bytes = (size_in_bytes + 15) & ~15u;
	D3D11_BUFFER_DESC description = {};
	description.ByteWidth = size_in_bytes;
	description.Usage = D3D11_USAGE_DEFAULT;
	description.BindFlags = bind_flags;
	description.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_ALLOW_RAW_VIEWS;
	D3D11_SUBRESOURCE_DATA initial_data = {};
	initial_data.pSysMem = data;
	HRESULT result = device->CreateBuffer(&description, data ? &initial_data : NULL, &buffer->buffer);
	if (FAILED(result)) {
		printf("Could not create a raw buffer of %u bytes (0x%08lx)\n", size_in_bytes, result);
		return false;
	}
	if (bind_flags & D3D11_BIND_SHADER_RESOURCE) {
		D3D11_SHADER_RESOURCE_VIEW_DESC view_description = {};
		view_description.Format = DXGI_FORMAT_R32_TYPELESS;
		view_description.ViewDimension = D3D11_SRV_DIMENSION_BUFFEREX;
		view_description.BufferEx.NumElements = size_in_bytes / 4;
		view_description.BufferEx.Flags = D3D11_BUFFEREX_SRV_FLAG_RAW;
		result = device->CreateShaderResourceView(buffer->buffer, &view_description, &buffer->view);
	}
	if (SUCCEEDED(result) && (bind_flags & D3D11_BIND_UNORDERED_ACCESS)) {
		D3D11_UNORDERED_ACCESS_VIEW_DESC unordered_description = {};
		unordered_description.Format = DXGI_FORMAT_R32_TYPELESS;
		unordered_description.ViewDimension = D3D11_UAV_DIMENSION_BUFFER;
		unordered_description.Buffer.NumElements = size_in_bytes / 4;
		unordered_description.Buffer.Flags = D3D11_BUFFER_UAV_FLAG_RAW;
		result = device->CreateUnorderedAccessView(buffer->buffer, &unordered_description, &buffer->unordered_view);
	}
	if (FAILED(result)) {
		printf("Could not create a raw buffer view (0x%08lx)\n", result);
		return false;
	}
	buffer->capacity = size_in_bytes;
	return true;
}

//Dynamic, so each frame's data is written straight into it with Map
static bool tfx_impl_dx11__create_upload_buffer(ID3D11Device *device, tfxU32 size_in_bytes, tfx_impl_dx11_buffer_t *buffer) {
	tfx_impl_dx11__release_buffer(buffer);
	size_in_bytes = (size_in_bytes + 15) & ~15u;
	D3D11_BUFFER_DESC description = {};
	description.ByteWidth = size_in_bytes;
	description.Usage = D3D11_USAGE_DYNAMIC;
	description.BindFlags = D3D11_BIND_SHADER_RESOURCE;
	description.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
	description.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_ALLOW_RAW_VIEWS;
	HRESULT result = device->CreateBuffer(&description, NULL, &buffer->buffer);
	if (SUCCEEDED(result)) {
		D3D11_SHADER_RESOURCE_VIEW_DESC view_description = {};
		view_description.Format = DXGI_FORMAT_R32_TYPELESS;
		view_description.ViewDimension = D3D11_SRV_DIMENSION_BUFFEREX;
		view_description.BufferEx.NumElements = size_in_bytes / 4;
		view_description.BufferEx.Flags = D3D11_BUFFEREX_SRV_FLAG_RAW;
		result = device->CreateShaderResourceView(buffer->buffer, &view_description, &buffer->view);
	}
	if (FAILED(result)) {
		printf("Could not create an upload buffer of %u bytes (0x%08lx)\n", size_in_bytes, result);
		return false;
	}
	buffer->capacity = size_in_bytes;
	return true;
}

//Grows by half again whenever it's too small, so the size settles after the first few frames
static bool tfx_impl_dx11__ensure_structured_buffer(ID3D11Device *device, tfxU32 element_size, tfxU32 element_count, bool unordered, tfx_impl_dx11_buffer_t *buffer) {
	if (buffer->buffer && buffer->capacity >= element_count) {
		return true;
	}
	tfxU32 capacity = element_count + element_count / 2 + 1;
	return tfx_impl_dx11__create_structured_buffer(device, NULL, element_size, capacity, unordered, buffer);
}

static bool tfx_impl_dx11__ensure_raw_buffer(ID3D11Device *device, tfxU32 size_in_bytes, UINT bind_flags, tfx_impl_dx11_buffer_t *buffer) {
	if (buffer->buffer && buffer->capacity >= size_in_bytes) {
		return true;
	}
	return tfx_impl_dx11__create_raw_buffer(device, NULL, size_in_bytes + size_in_bytes / 2, bind_flags, buffer);
}

static void tfx_impl_dx11__update_buffer_bytes(ID3D11DeviceContext *device_context, ID3D11Buffer *buffer, const void *data, tfxU32 size_in_bytes) {
	if (size_in_bytes == 0) {
		return;
	}
	D3D11_BOX box = { 0, 0, 0, size_in_bytes, 1, 1 };
	device_context->UpdateSubresource(buffer, 0, &box, data, 0, 0);
}

static void tfx_impl_dx11__write_constants(ID3D11DeviceContext *device_context, ID3D11Buffer *buffer, const void *data, size_t size) {
	D3D11_MAPPED_SUBRESOURCE mapped;
	if (SUCCEEDED(device_context->Map(buffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped))) {
		memcpy(mapped.pData, data, size);
		device_context->Unmap(buffer, 0);
	}
}

static ID3D11Buffer *tfx_impl_dx11__create_constant_buffer(ID3D11Device *device, size_t size) {
	D3D11_BUFFER_DESC description = {};
	description.Usage = D3D11_USAGE_DYNAMIC;
	description.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
	description.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
	description.ByteWidth = (UINT)((size + 15) & ~(size_t)15);
	ID3D11Buffer *buffer = NULL;
	device->CreateBuffer(&description, NULL, &buffer);
	return buffer;
}

//slangc compiles every shader ahead of time to <entry point>.cso in the shader folder. The caller frees the code
static void *tfx_impl_dx11__load_shader_code(const char *shader_folder, const char *entry_point, size_t *size) {
	char path[MAX_PATH];
	snprintf(path, sizeof(path), "%s/%s.cso", shader_folder, entry_point);
	FILE *file = fopen(path, "rb");
	if (!file) {
		printf("Could not open the shader %s\n", path);
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
		printf("Could not read the shader %s\n", path);
		return NULL;
	}
	*size = (size_t)file_size;
	return code;
}

static bool tfx_impl_dx11__create_vertex_shader(ID3D11Device *device, const char *shader_folder, const char *entry_point, ID3D11VertexShader **shader) {
	size_t size = 0;
	void *code = tfx_impl_dx11__load_shader_code(shader_folder, entry_point, &size);
	bool created = code && SUCCEEDED(device->CreateVertexShader(code, size, NULL, shader));
	free(code);
	return created;
}

static bool tfx_impl_dx11__create_pixel_shader(ID3D11Device *device, const char *shader_folder, const char *entry_point, ID3D11PixelShader **shader) {
	size_t size = 0;
	void *code = tfx_impl_dx11__load_shader_code(shader_folder, entry_point, &size);
	bool created = code && SUCCEEDED(device->CreatePixelShader(code, size, NULL, shader));
	free(code);
	return created;
}

static bool tfx_impl_dx11__create_compute_shader(ID3D11Device *device, const char *shader_folder, const char *entry_point, ID3D11ComputeShader **shader) {
	size_t size = 0;
	void *code = tfx_impl_dx11__load_shader_code(shader_folder, entry_point, &size);
	bool created = code && SUCCEEDED(device->CreateComputeShader(code, size, NULL, shader));
	free(code);
	return created;
}

static DXGI_FORMAT tfx_impl_dx11__color_ramp_format(tfxU32 bytes_per_pixel) {
	if (bytes_per_pixel == 4) {
		return DXGI_FORMAT_R8G8B8A8_UNORM;
	} else if (bytes_per_pixel == 16) {
		return DXGI_FORMAT_R32G32B32A32_FLOAT;
	}
	return DXGI_FORMAT_R16G16B16A16_FLOAT;
}

//One 256 x 256 layer per color ramp bitmap, each row is an emitter's ramp
static bool tfx_impl_dx11__create_color_ramps(ID3D11Device *device, tfx_bitmap_t **bitmaps, tfxU32 bitmap_count, ID3D11Texture2D **texture, ID3D11ShaderResourceView **view) {
	TFX_IMPL_DX11_RELEASE(*view);
	TFX_IMPL_DX11_RELEASE(*texture);
	if (bitmap_count == 0) {
		return true;
	}
	tfxU32 width = (tfxU32)tfx_GetBitmapWidth(bitmaps[0]);
	tfxU32 height = (tfxU32)tfx_GetBitmapHeight(bitmaps[0]);
	tfxU32 bytes_per_pixel = (tfxU32)(tfx_GetBitmapSize(bitmaps[0]) / ((size_t)width * height));
	D3D11_SUBRESOURCE_DATA *layers = (D3D11_SUBRESOURCE_DATA *)malloc(sizeof(D3D11_SUBRESOURCE_DATA) * bitmap_count);
	if (!layers) {
		return false;
	}
	for (tfxU32 index = 0; index != bitmap_count; ++index) {
		layers[index].pSysMem = tfx_GetBitmapData(bitmaps[index]);
		layers[index].SysMemPitch = width * bytes_per_pixel;
		layers[index].SysMemSlicePitch = 0;
	}
	D3D11_TEXTURE2D_DESC ramp_description = {};
	ramp_description.Width = width;
	ramp_description.Height = height;
	ramp_description.MipLevels = 1;
	ramp_description.ArraySize = bitmap_count;
	ramp_description.Format = tfx_impl_dx11__color_ramp_format(bytes_per_pixel);
	ramp_description.SampleDesc.Count = 1;
	ramp_description.Usage = D3D11_USAGE_IMMUTABLE;
	ramp_description.BindFlags = D3D11_BIND_SHADER_RESOURCE;
	HRESULT result = device->CreateTexture2D(&ramp_description, layers, texture);
	free(layers);
	if (FAILED(result) || FAILED(device->CreateShaderResourceView(*texture, NULL, view))) {
		printf("Could not create a color ramp texture (0x%08lx)\n", result);
		return false;
	}
	return true;
}

//--- Dx11_functions

bool tfx_impl_dx11_Init(tfx_impl_dx11_t *renderer, ID3D11Device *device, ID3D11DeviceContext *device_context, const char *shader_folder) {
	memset(renderer, 0, sizeof(tfx_impl_dx11_t));
	renderer->device = device;
	renderer->device_context = device_context;
	device->AddRef();
	device_context->AddRef();

	if (!tfx_impl_dx11__create_vertex_shader(device, shader_folder, "tfx_billboard_vertex", &renderer->billboard_vertex_shader)
		|| !tfx_impl_dx11__create_pixel_shader(device, shader_folder, "tfx_billboard_fragment", &renderer->billboard_pixel_shader)
		|| !tfx_impl_dx11__create_compute_shader(device, shader_folder, "tfx_ribbon_compute", &renderer->ribbon_compute_shader)
		|| !tfx_impl_dx11__create_vertex_shader(device, shader_folder, "tfx_ribbon_vertex", &renderer->ribbon_vertex_shader)
		|| !tfx_impl_dx11__create_pixel_shader(device, shader_folder, "tfx_ribbon_fragment", &renderer->ribbon_pixel_shader)
		|| !tfx_impl_dx11__create_compute_shader(device, shader_folder, "tfx_sprite_playback", &renderer->sprite_playback_shader)
		|| !tfx_impl_dx11__create_vertex_shader(device, shader_folder, "tfx_billboard_baked_vertex", &renderer->billboard_static_vertex_shader)
		|| !tfx_impl_dx11__create_compute_shader(device, shader_folder, "tfx_ribbon_playback", &renderer->ribbon_playback_shader)
		|| !tfx_impl_dx11__create_compute_shader(device, shader_folder, "simulate_particles", &renderer->particle_kernel)) {
		return false;
	}

	D3D11_BLEND_DESC blend_description = {};
	blend_description.RenderTarget[0].BlendEnable = TRUE;
	blend_description.RenderTarget[0].SrcBlend = D3D11_BLEND_ONE;
	blend_description.RenderTarget[0].DestBlend = D3D11_BLEND_INV_SRC_ALPHA;
	blend_description.RenderTarget[0].BlendOp = D3D11_BLEND_OP_ADD;
	blend_description.RenderTarget[0].SrcBlendAlpha = D3D11_BLEND_ONE;
	blend_description.RenderTarget[0].DestBlendAlpha = D3D11_BLEND_INV_SRC_ALPHA;
	blend_description.RenderTarget[0].BlendOpAlpha = D3D11_BLEND_OP_ADD;
	blend_description.RenderTarget[0].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;

	D3D11_DEPTH_STENCIL_DESC depth_description = {};
	depth_description.DepthEnable = TRUE;
	depth_description.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ZERO;
	depth_description.DepthFunc = D3D11_COMPARISON_LESS_EQUAL;

	D3D11_RASTERIZER_DESC rasterizer_description = {};
	rasterizer_description.FillMode = D3D11_FILL_SOLID;
	rasterizer_description.CullMode = D3D11_CULL_NONE;
	rasterizer_description.DepthClipEnable = TRUE;

	D3D11_SAMPLER_DESC sampler_description = {};
	sampler_description.Filter = D3D11_FILTER_MIN_MAG_MIP_LINEAR;
	sampler_description.AddressU = D3D11_TEXTURE_ADDRESS_CLAMP;
	sampler_description.AddressV = D3D11_TEXTURE_ADDRESS_CLAMP;
	sampler_description.AddressW = D3D11_TEXTURE_ADDRESS_CLAMP;
	sampler_description.MaxLOD = D3D11_FLOAT32_MAX;

	renderer->frame_constants = tfx_impl_dx11__create_constant_buffer(device, sizeof(tfx_impl_uniforms_t));
	renderer->draw_constants = tfx_impl_dx11__create_constant_buffer(device, sizeof(tfx_impl_dx11_draw_constants_t));
	renderer->ribbon_constants = tfx_impl_dx11__create_constant_buffer(device, sizeof(tfx_ribbon_bucket_globals_t));
	renderer->particle_push = tfx_impl_dx11__create_constant_buffer(device, sizeof(tfx_gpu_particle_push_t));
	renderer->sprite_data_constants = tfx_impl_dx11__create_constant_buffer(device, sizeof(tfx_sprite_data_push_t));
	renderer->ribbon_data_constants = tfx_impl_dx11__create_constant_buffer(device, sizeof(tfx_ribbon_data_push_t));
	if (FAILED(device->CreateBlendState(&blend_description, &renderer->blend_state))
		|| FAILED(device->CreateDepthStencilState(&depth_description, &renderer->depth_stencil_state))
		|| FAILED(device->CreateRasterizerState(&rasterizer_description, &renderer->rasterizer_state))
		|| FAILED(device->CreateSamplerState(&sampler_description, &renderer->sampler_state))
		|| !renderer->frame_constants || !renderer->draw_constants || !renderer->ribbon_constants || !renderer->particle_push || !renderer->sprite_data_constants
		|| !renderer->ribbon_data_constants) {
		printf("Could not create the TimelineFX render states\n");
		return false;
	}
	return true;
}

bool tfx_impl_dx11_UploadLibrary(tfx_impl_dx11_t *renderer, tfx_library library, const tfx_impl_atlas_t *atlas) {
	ID3D11Device *device = renderer->device;

	//The atlas makes its own mip chain, so it's created empty, its top level uploaded and the rest generated
	TFX_IMPL_DX11_RELEASE(renderer->shape_atlas_view);
	TFX_IMPL_DX11_RELEASE(renderer->shape_atlas);
	D3D11_TEXTURE2D_DESC atlas_description = {};
	atlas_description.Width = atlas->layer_size;
	atlas_description.Height = atlas->layer_size;
	atlas_description.MipLevels = 0;
	atlas_description.ArraySize = atlas->layer_count;
	atlas_description.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
	atlas_description.SampleDesc.Count = 1;
	atlas_description.Usage = D3D11_USAGE_DEFAULT;
	atlas_description.BindFlags = D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_RENDER_TARGET;
	atlas_description.MiscFlags = D3D11_RESOURCE_MISC_GENERATE_MIPS;
	HRESULT result = device->CreateTexture2D(&atlas_description, NULL, &renderer->shape_atlas);
	if (FAILED(result) || FAILED(device->CreateShaderResourceView(renderer->shape_atlas, NULL, &renderer->shape_atlas_view))) {
		printf("Could not create the %u x %u x %u shape atlas (0x%08lx)\n", atlas->layer_size, atlas->layer_size, atlas->layer_count, result);
		return false;
	}
	renderer->shape_atlas->GetDesc(&atlas_description);
	size_t layer_bytes = (size_t)atlas->layer_size * atlas->layer_size * 4;
	for (tfxU32 layer = 0; layer != atlas->layer_count; ++layer) {
		UINT subresource = D3D11CalcSubresource(0, layer, atlas_description.MipLevels);
		renderer->device_context->UpdateSubresource(renderer->shape_atlas, subresource, NULL, atlas->pixels + layer * layer_bytes, atlas->layer_size * 4, 0);
	}
	renderer->device_context->GenerateMips(renderer->shape_atlas_view);

	tfxU32 bitmap_count = tfx_GetColorRampBitmapCount(library);
	tfx_bitmap_t **bitmaps = (tfx_bitmap_t **)malloc(sizeof(tfx_bitmap_t *) * (bitmap_count ? bitmap_count : 1));
	for (tfxU32 index = 0; bitmaps && index != bitmap_count; ++index) {
		bitmaps[index] = tfx_GetColorRampBitmap(library, index);
	}
	bool ramps_created = bitmaps && tfx_impl_dx11__create_color_ramps(device, bitmaps, bitmap_count, &renderer->color_ramps, &renderer->color_ramps_view);
	free(bitmaps);
	if (!ramps_created) {
		return false;
	}

	tfx_gpu_shapes gpu_shapes = tfx_GetLibraryGPUShapes(library);
	tfxU32 shape_count = tfx_GetGPUShapesCount(gpu_shapes);
	tfx_gpu_image_data_t empty_shape = {};
	if (!tfx_impl_dx11__create_structured_buffer(device, shape_count ? (const void *)tfx_GetGPUShapesArray(gpu_shapes) : (const void *)&empty_shape,
		sizeof(tfx_gpu_image_data_t), shape_count ? shape_count : 1, false, &renderer->image_data)) {
		return false;
	}

	tfxU32 property_count = (tfxU32)(tfx_GetParticlePropertiesBufferSizeInBytes(library) / sizeof(tfx_gpu_particle_properties_t));
	tfx_gpu_particle_properties_t empty_properties = {};
	if (!tfx_impl_dx11__create_structured_buffer(device, property_count ? tfx_GetParticlePropertiesBuffer(library) : (const void *)&empty_properties,
		sizeof(tfx_gpu_particle_properties_t), property_count ? property_count : 1, false, &renderer->emitter_properties)) {
		return false;
	}

	//The graph lookups are shared by every library and read by the ribbon shaders
	tfxU32 graph_bytes = tfx_GetGPUGraphLookupsBufferSizeInBytes();
	tfx_gpu_graph_data_t empty_graph = {};
	if (!tfx_impl_dx11__create_raw_buffer(device, graph_bytes ? (const void *)tfx_GetGPUGraphLookupsBuffer() : (const void *)&empty_graph,
		graph_bytes ? graph_bytes : sizeof(tfx_gpu_graph_data_t), D3D11_BIND_SHADER_RESOURCE, &renderer->graph_lookups)) {
		return false;
	}
	return true;
}

bool tfx_impl_dx11_UploadSpriteData(tfx_impl_dx11_t *renderer, tfx_animation_manager animation_manager, tfx_gpu_shapes shapes) {
	ID3D11Device *device = renderer->device;
	tfxU32 bitmap_count = tfx_GetAnimationColorRampBitmapCount(animation_manager);
	tfx_bitmap_t **bitmaps = (tfx_bitmap_t **)malloc(sizeof(tfx_bitmap_t *) * (bitmap_count ? bitmap_count : 1));
	for (tfxU32 index = 0; bitmaps && index != bitmap_count; ++index) {
		bitmaps[index] = tfx_GetAnimationColorRampBitmap(animation_manager, index);
	}
	bool ramps_created = bitmaps && tfx_impl_dx11__create_color_ramps(device, bitmaps, bitmap_count, &renderer->baked_color_ramps, &renderer->baked_color_ramps_view);
	free(bitmaps);
	if (!ramps_created) {
		return false;
	}

	tfxU32 shape_count = tfx_GetGPUShapesCount(shapes);
	tfx_gpu_image_data_t empty_shape = {};
	tfxU32 property_count = (tfxU32)(tfx_GetAnimationEmitterPropertySizeInBytes(animation_manager) / sizeof(tfx_gpu_particle_properties_t));
	tfx_gpu_particle_properties_t empty_properties = {};
	size_t sprite_bytes = tfx_GetSpriteDataSizeInBytes(animation_manager);
	tfx_instance_t empty_sprite = {};
	if (!tfx_impl_dx11__create_structured_buffer(device, shape_count ? (const void *)tfx_GetGPUShapesArray(shapes) : (const void *)&empty_shape,
			sizeof(tfx_gpu_image_data_t), shape_count ? shape_count : 1, false, &renderer->baked_image_data)
		|| !tfx_impl_dx11__create_structured_buffer(device, property_count ? tfx_GetAnimationEmitterPropertiesBufferPointer(animation_manager) : (const void *)&empty_properties,
			sizeof(tfx_gpu_particle_properties_t), property_count ? property_count : 1, false, &renderer->baked_emitter_properties)
		|| !tfx_impl_dx11__create_raw_buffer(device, sprite_bytes ? tfx_GetSpriteDataBufferPointer(animation_manager) : (const void *)&empty_sprite,
			(tfxU32)(sprite_bytes ? sprite_bytes : sizeof(tfx_instance_t)), D3D11_BIND_SHADER_RESOURCE, &renderer->baked_sprites)) {
		return false;
	}
	if (!tfx_AnimationManagerHasRibbons(animation_manager)) {
		return true;
	}

	//The baked ribbons, their properties, segments, emitters and graph lookups only change when sprite data is loaded
	tfxU32 graph_bytes = (tfxU32)tfx_GetAnimationRibbonGraphLookupsSizeInBytes(animation_manager);
	tfx_gpu_graph_data_t empty_graph = {};
	if (!tfx_impl_dx11__create_raw_buffer(device, graph_bytes ? tfx_GetAnimationRibbonGraphLookupsBufferPointer(animation_manager) : (const void *)&empty_graph,
		graph_bytes ? graph_bytes : sizeof(tfx_gpu_graph_data_t), D3D11_BIND_SHADER_RESOURCE, &renderer->baked_ribbons.graph_lookups)) {
		return false;
	}
	size_t emitters_bytes = tfx_GetAnimationRibbonEmittersSizeInBytes(animation_manager);
	void *emitters = malloc(emitters_bytes);
	if (!emitters) {
		return false;
	}
	tfx_CopyAnimationRibbonEmitters(animation_manager, emitters);
	bool created = tfx_impl_dx11__create_structured_buffer(device, tfx_GetAnimationRibbonDataBufferPointer(animation_manager), TFX_IMPL_DX11_RIBBON_DATA_SIZE,
			(tfxU32)(tfx_GetAnimationRibbonDataSizeInBytes(animation_manager) / TFX_IMPL_DX11_RIBBON_DATA_SIZE), false, &renderer->baked_ribbon_data)
		&& tfx_impl_dx11__create_structured_buffer(device, tfx_GetAnimationRibbonPropertiesBufferPointer(animation_manager), TFX_IMPL_DX11_RIBBON_PROPERTIES_SIZE,
			(tfxU32)(tfx_GetAnimationRibbonPropertySizeInBytes(animation_manager) / TFX_IMPL_DX11_RIBBON_PROPERTIES_SIZE), false, &renderer->baked_ribbon_properties)
		&& tfx_impl_dx11__create_structured_buffer(device, tfx_GetAnimationRibbonSegmentsBufferPointer(animation_manager), TFX_IMPL_DX11_RIBBON_SEGMENT_SIZE,
			(tfxU32)(tfx_GetAnimationRibbonSegmentsSizeInBytes(animation_manager) / TFX_IMPL_DX11_RIBBON_SEGMENT_SIZE), false, &renderer->baked_ribbons.segments)
		&& tfx_impl_dx11__create_structured_buffer(device, emitters, TFX_IMPL_DX11_RIBBON_EMITTER_SIZE,
			(tfxU32)(emitters_bytes / TFX_IMPL_DX11_RIBBON_EMITTER_SIZE), false, &renderer->baked_ribbons.emitters);
	free(emitters);
	return created;
}

static bool tfx_impl_dx11__write_upload_buffer(tfx_impl_dx11_t *renderer, tfx_impl_dx11_buffer_t *buffer, const void *data, size_t size_in_bytes) {
	if (buffer->capacity < size_in_bytes && !tfx_impl_dx11__create_upload_buffer(renderer->device, (tfxU32)(size_in_bytes + size_in_bytes / 2), buffer)) {
		return false;
	}
	D3D11_MAPPED_SUBRESOURCE mapped;
	if (FAILED(renderer->device_context->Map(buffer->buffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped))) {
		return false;
	}
	memcpy(mapped.pData, data, size_in_bytes);
	renderer->device_context->Unmap(buffer->buffer, 0);
	return true;
}

//The vertex and index buffers ribbon_compute writes, sized for this frame's buckets
static bool tfx_impl_dx11__ensure_ribbon_geometry(ID3D11Device *device, const tfx_ribbon_batch_sizes_t *sizes, tfx_impl_dx11_ribbon_set_t *set) {
	return tfx_impl_dx11__ensure_structured_buffer(device, sizeof(tfx_ribbon_vertex_t), (tfxU32)(sizes->vertex_buffer_size_in_bytes / sizeof(tfx_ribbon_vertex_t)), true, &set->vertices)
		&& tfx_impl_dx11__ensure_raw_buffer(device, (tfxU32)sizes->index_buffer_size_in_bytes, D3D11_BIND_UNORDERED_ACCESS | D3D11_BIND_INDEX_BUFFER, &set->indices);
}

static bool tfx_impl_dx11__upload_baked_ribbons(tfx_impl_dx11_t *renderer, const tfx_impl_baked_frame_t *baked) {
	if (baked->ribbon_count == 0 || !renderer->baked_ribbon_data.buffer) {
		return true;
	}
	tfx_ribbon_batch_sizes_t sizes = tfx_GetAnimationRibbonBufferSizes(baked->animation_manager, sizeof(tfx_ribbon_vertex_t));
	//ribbon_playback writes the ribbons and ribbon_compute reads them, so they need both views
	if (!tfx_impl_dx11__write_upload_buffer(renderer, &renderer->animation_ribbon_offsets, baked->ribbon_offsets, baked->ribbon_offsets_bytes)
		|| !tfx_impl_dx11__ensure_structured_buffer(renderer->device, TFX_IMPL_DX11_RIBBON_SIZE, baked->ribbon_count, true, &renderer->baked_ribbons.instances)
		|| !tfx_impl_dx11__ensure_ribbon_geometry(renderer->device, &sizes, &renderer->baked_ribbons)) {
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
static bool tfx_impl_dx11__upload_baked(tfx_impl_dx11_t *renderer, const tfx_impl_baked_frame_t *baked) {
	renderer->playback_sprite_count = 0;
	renderer->playback_ribbon_count = 0;
	renderer->baked_ribbons.animation_manager = NULL;
	if (!renderer->baked_sprites.buffer || (baked->sprite_count == 0 && baked->ribbon_count == 0)) {
		return true;
	}
	if (!tfx_impl_dx11__write_upload_buffer(renderer, &renderer->animation_instances, baked->animation_instances, baked->animation_instances_bytes)) {
		return false;
	}
	if (baked->sprite_count > 0) {
		if (!tfx_impl_dx11__write_upload_buffer(renderer, &renderer->animation_offsets, baked->offsets, baked->offsets_bytes)
			|| !tfx_impl_dx11__ensure_raw_buffer(renderer->device, baked->sprite_count * sizeof(tfx_instance_t), D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS, &renderer->playback_instances)) {
			return false;
		}
		memset(&renderer->sprite_data_push, 0, sizeof(tfx_sprite_data_push_t));
		renderer->sprite_data_push.animation_instances_total = baked->animation_instance_count;
		renderer->sprite_data_push.billboards_total = baked->sprite_count;
		renderer->sprite_data_push.animated_shapes = baked->animated_shapes ? 1 : 0;
		renderer->playback_sprite_count = baked->sprite_count;
	}
	return tfx_impl_dx11__upload_baked_ribbons(renderer, baked);
}

static bool tfx_impl_dx11__upload_ribbons(tfx_impl_dx11_t *renderer, const tfx_impl_ribbon_frame_t *ribbons) {
	tfx_impl_dx11_ribbon_set_t *set = &renderer->live_ribbons;
	set->stage = ribbons->stage;
	if (!ribbons->stage) {
		return true;
	}
	ID3D11Device *device = renderer->device;
	const tfx_ribbon_batch_sizes_t *sizes = &ribbons->sizes;
	if (!tfx_impl_dx11__ensure_structured_buffer(device, TFX_IMPL_DX11_RIBBON_SEGMENT_SIZE, (tfxU32)(sizes->segment_buffer_size_in_bytes / TFX_IMPL_DX11_RIBBON_SEGMENT_SIZE), false, &set->segments)
		|| !tfx_impl_dx11__ensure_structured_buffer(device, TFX_IMPL_DX11_RIBBON_SIZE, (tfxU32)(sizes->ribbon_buffer_size_in_bytes / TFX_IMPL_DX11_RIBBON_SIZE), false, &set->instances)
		|| !tfx_impl_dx11__ensure_structured_buffer(device, TFX_IMPL_DX11_RIBBON_EMITTER_SIZE, (tfxU32)(sizes->emitter_buffer_size_in_bytes / TFX_IMPL_DX11_RIBBON_EMITTER_SIZE), false, &set->emitters)
		|| !tfx_impl_dx11__ensure_ribbon_geometry(device, sizes, set)) {
		set->stage = NULL;
		return false;
	}
	ID3D11DeviceContext *device_context = renderer->device_context;
	tfx_impl_dx11__update_buffer_bytes(device_context, set->segments.buffer, ribbons->segments, (tfxU32)sizes->segment_buffer_size_in_bytes);
	tfx_impl_dx11__update_buffer_bytes(device_context, set->instances.buffer, ribbons->ribbons, (tfxU32)sizes->ribbon_buffer_size_in_bytes);
	tfx_impl_dx11__update_buffer_bytes(device_context, set->emitters.buffer, ribbons->emitters, (tfxU32)sizes->emitter_buffer_size_in_bytes);
	set->stage_offsets = ribbons->stage_offsets;
	return true;
}

//The old contents are copied across, because the previous buffer is still what captured_index points into
static bool tfx_impl_dx11__grow_instance_buffers(tfx_impl_dx11_t *renderer, tfxU32 needed) {
	tfxU32 capacity = renderer->instance_capacity ? renderer->instance_capacity : TFX_IMPL_DX11_INITIAL_INSTANCE_CAPACITY;
	while (capacity < needed) {
		capacity += capacity / 2;
	}
	//The GPU particle kernel writes its instances after the CPU ones
	UINT bind_flags = D3D11_BIND_SHADER_RESOURCE | (renderer->particle_kernel ? D3D11_BIND_UNORDERED_ACCESS : 0);
	for (int index = 0; index != 2; ++index) {
		tfx_impl_dx11_buffer_t grown = {};
		if (!tfx_impl_dx11__create_raw_buffer(renderer->device, NULL, capacity * sizeof(tfx_instance_t), bind_flags, &grown)) {
			return false;
		}
		tfx_impl_dx11_buffer_t *old = &renderer->instance_buffers[index];
		if (old->buffer) {
			D3D11_BOX box = { 0, 0, 0, old->capacity, 1, 1 };
			renderer->device_context->CopySubresourceRegion(grown.buffer, 0, 0, 0, 0, old->buffer, 0, &box);
		}
		tfx_impl_dx11__release_buffer(old);
		*old = grown;
	}
	renderer->instance_capacity = capacity;
	return true;
}

static bool tfx_impl_dx11__upload_gpu_particles(tfx_impl_dx11_t *renderer, const tfx_impl_frame_t *frame) {
	const tfx_gpu_particle_frame_t *gpu_frame = &frame->gpu_particles;
	//The pool holds the particles between frames, it's sized once and starts zeroed
	if (!renderer->particle_pool.buffer) {
		void *zeroes = calloc(1, frame->gpu_pool_sizes.pool_bytes);
		bool created = zeroes && tfx_impl_dx11__create_raw_buffer(renderer->device, zeroes, (tfxU32)frame->gpu_pool_sizes.pool_bytes,
			D3D11_BIND_UNORDERED_ACCESS, &renderer->particle_pool);
		free(zeroes);
		if (!created) {
			return false;
		}
	}
	if (renderer->particle_upload.capacity < gpu_frame->upload_bytes
		&& !tfx_impl_dx11__create_upload_buffer(renderer->device, (tfxU32)(gpu_frame->upload_bytes + gpu_frame->upload_bytes / 2), &renderer->particle_upload)) {
		return false;
	}
	D3D11_MAPPED_SUBRESOURCE mapped;
	if (FAILED(renderer->device_context->Map(renderer->particle_upload.buffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped))) {
		return false;
	}
	tfx_CopyGPUParticleUpload(gpu_frame, mapped.pData);
	renderer->device_context->Unmap(renderer->particle_upload.buffer, 0);

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

bool tfx_impl_dx11_UploadFrame(tfx_impl_dx11_t *renderer, const tfx_impl_frame_t *frame) {
	if (!tfx_impl_dx11__upload_baked(renderer, &frame->baked)) {
		return false;
	}
	if (frame->ticks == 0) {
		return true;
	}
	bool gpu_particles = frame->has_gpu_particles && renderer->particle_kernel;
	tfxU32 needed = gpu_particles ? frame->gpu_particles.gpu_instance_start + frame->gpu_particles.gpu_instance_count : frame->instance_count;
	if ((needed > renderer->instance_capacity || renderer->instance_capacity == 0) && !tfx_impl_dx11__grow_instance_buffers(renderer, needed)) {
		renderer->instance_count = 0;
		return false;
	}

	//The other buffer now holds what was drawn last, which is what captured_index points into for CPU and GPU instances alike
	renderer->current_instance_buffer ^= 1;
	tfx_impl_dx11__update_buffer_bytes(renderer->device_context, renderer->instance_buffers[renderer->current_instance_buffer].buffer,
		frame->instances, frame->instance_count * sizeof(tfx_instance_t));
	renderer->instance_count = needed;
	memcpy(renderer->layer_start, frame->layer_start, sizeof(renderer->layer_start));
	memcpy(renderer->layer_count, frame->layer_count, sizeof(renderer->layer_count));
	memset(renderer->gpu_layer_count, 0, sizeof(renderer->gpu_layer_count));
	if (gpu_particles) {
		memcpy(renderer->gpu_layer_start, frame->gpu_particles.gpu_layer_start, sizeof(renderer->gpu_layer_start));
		memcpy(renderer->gpu_layer_count, frame->gpu_particles.gpu_layer_count, sizeof(renderer->gpu_layer_count));
		if (!tfx_impl_dx11__upload_gpu_particles(renderer, frame)) {
			return false;
		}
	}
	return tfx_impl_dx11__upload_ribbons(renderer, &frame->ribbons);
}

//Every fetch has to be dispatched, the next one carries on from the particle state this one leaves in the pool
static void tfx_impl_dx11__simulate_gpu_particles(tfx_impl_dx11_t *renderer) {
	ID3D11DeviceContext *device_context = renderer->device_context;
	ID3D11UnorderedAccessView *unordered_views[2] = { renderer->particle_pool.unordered_view, renderer->instance_buffers[renderer->current_instance_buffer].unordered_view };
	ID3D11ShaderResourceView *views[2] = { renderer->graph_lookups.view, renderer->particle_upload.view };
	device_context->CSSetShader(renderer->particle_kernel, NULL, 0);
	device_context->CSSetUnorderedAccessViews(0, 2, unordered_views, NULL);
	device_context->CSSetShaderResources(0, 2, views);
	device_context->CSSetConstantBuffers(0, 1, &renderer->particle_push);
	tfx_gpu_particle_push_t *push = &renderer->particle_push_constants;
	for (tfxU32 index = 0; index != renderer->particle_dispatch_count; ++index) {
		tfxU32 slot_count = renderer->particle_dispatches[index].slot_count;
		if (slot_count == 0) {
			continue;
		}
		push->dispatch_index = index;
		tfx_impl_dx11__write_constants(device_context, renderer->particle_push, push, sizeof(tfx_gpu_particle_push_t));
		device_context->Dispatch((slot_count + 255) / 256, 1, 1);
	}
	renderer->particle_dispatch_count = 0;
	ID3D11UnorderedAccessView *no_unordered_views[2] = {};
	ID3D11ShaderResourceView *no_views[2] = {};
	device_context->CSSetUnorderedAccessViews(0, 2, no_unordered_views, NULL);
	device_context->CSSetShaderResources(0, 2, no_views);
	device_context->CSSetShader(NULL, NULL, 0);
}

static void tfx_impl_dx11__play_back_sprites(tfx_impl_dx11_t *renderer) {
	ID3D11DeviceContext *device_context = renderer->device_context;
	ID3D11ShaderResourceView *sprite_views[3] = { renderer->baked_sprites.view, NULL, renderer->baked_emitter_properties.view };
	ID3D11ShaderResourceView *animation_views[2] = { renderer->animation_offsets.view, renderer->animation_instances.view };
	device_context->CSSetShader(renderer->sprite_playback_shader, NULL, 0);
	device_context->CSSetShaderResources(1, 3, sprite_views);
	device_context->CSSetShaderResources(11, 2, animation_views);
	device_context->CSSetUnorderedAccessViews(2, 1, &renderer->playback_instances.unordered_view, NULL);
	tfx_impl_dx11__write_constants(device_context, renderer->sprite_data_constants, &renderer->sprite_data_push, sizeof(tfx_sprite_data_push_t));
	device_context->CSSetConstantBuffers(3, 1, &renderer->sprite_data_constants);
	device_context->Dispatch((renderer->playback_sprite_count + 127) / 128, 1, 1);

	ID3D11ShaderResourceView *no_views[3] = {};
	ID3D11UnorderedAccessView *no_unordered_view = NULL;
	device_context->CSSetShaderResources(1, 3, no_views);
	device_context->CSSetShaderResources(11, 2, no_views);
	device_context->CSSetUnorderedAccessViews(2, 1, &no_unordered_view, NULL);
	device_context->CSSetShader(NULL, NULL, 0);
}

static void tfx_impl_dx11__play_back_ribbons(tfx_impl_dx11_t *renderer) {
	ID3D11DeviceContext *device_context = renderer->device_context;
	ID3D11ShaderResourceView *views[4] = { renderer->animation_ribbon_offsets.view, renderer->animation_instances.view, renderer->baked_ribbon_data.view, renderer->baked_ribbon_properties.view };
	device_context->CSSetShader(renderer->ribbon_playback_shader, NULL, 0);
	device_context->CSSetShaderResources(11, 4, views);
	device_context->CSSetUnorderedAccessViews(3, 1, &renderer->baked_ribbons.instances.unordered_view, NULL);
	tfx_impl_dx11__write_constants(device_context, renderer->ribbon_data_constants, &renderer->ribbon_data_push, sizeof(tfx_ribbon_data_push_t));
	device_context->CSSetConstantBuffers(4, 1, &renderer->ribbon_data_constants);
	device_context->Dispatch((renderer->playback_ribbon_count + 127) / 128, 1, 1);

	ID3D11ShaderResourceView *no_views[4] = {};
	ID3D11UnorderedAccessView *no_unordered_view = NULL;
	device_context->CSSetShaderResources(11, 4, no_views);
	device_context->CSSetUnorderedAccessViews(3, 1, &no_unordered_view, NULL);
	device_context->CSSetShader(NULL, NULL, 0);
}

static bool tfx_impl_dx11__ribbon_set_active(const tfx_impl_dx11_ribbon_set_t *set) {
	return set->stage || set->animation_manager;
}

static tfx_ribbon_dispatch_t tfx_impl_dx11__create_ribbon_dispatch(const tfx_impl_dx11_ribbon_set_t *set) {
	return set->stage ? tfx_CreateRibbonBatchDispatch(&set->stage_offsets) : tfx_CreateAnimationRibbonDispatch();
}

static bool tfx_impl_dx11__next_ribbon_dispatch(const tfx_impl_dx11_ribbon_set_t *set, tfx_ribbon_dispatch_t *dispatch) {
	return set->stage ? tfx_NextRibbonDispatch(set->stage, dispatch) : tfx_NextAnimationRibbonDispatch(set->animation_manager, dispatch);
}

static ID3D11ShaderResourceView *tfx_impl_dx11__ribbon_graph_lookups(const tfx_impl_dx11_t *renderer, const tfx_impl_dx11_ribbon_set_t *set) {
	return set->stage ? renderer->graph_lookups.view : set->graph_lookups.view;
}

static void tfx_impl_dx11__build_ribbon_geometry(tfx_impl_dx11_t *renderer, const tfx_impl_dx11_ribbon_set_t *set, const tfx_impl_uniforms_t *uniforms) {
	ID3D11DeviceContext *device_context = renderer->device_context;
	ID3D11ShaderResourceView *views[4] = { set->emitters.view, tfx_impl_dx11__ribbon_graph_lookups(renderer, set), set->instances.view, set->segments.view };
	ID3D11UnorderedAccessView *unordered_views[2] = { set->vertices.unordered_view, set->indices.unordered_view };
	device_context->CSSetShader(renderer->ribbon_compute_shader, NULL, 0);
	device_context->CSSetShaderResources(6, 4, views);
	device_context->CSSetUnorderedAccessViews(0, 2, unordered_views, NULL);
	device_context->CSSetConstantBuffers(2, 1, &renderer->ribbon_constants);

	//One dispatch per bucket, each with the bucket's own globals
	tfx_ribbon_dispatch_t dispatch = tfx_impl_dx11__create_ribbon_dispatch(set);
	while (tfx_impl_dx11__next_ribbon_dispatch(set, &dispatch)) {
		tfx_ribbon_bucket_globals_t *globals = tfx_GetRibbonDispatchGlobals(&dispatch);
		if (set->stage) {
			globals->lerp = uniforms->timer_lerp;
		} else {
			//Baked ribbons are already interpolated, but they still face the camera
			globals->camera_position.x = uniforms->camera_position[0];
			globals->camera_position.y = uniforms->camera_position[1];
			globals->camera_position.z = uniforms->camera_position[2];
		}
		tfx_impl_dx11__write_constants(device_context, renderer->ribbon_constants, globals, sizeof(tfx_ribbon_bucket_globals_t));
		device_context->Dispatch((dispatch.total_segments + TFX_IMPL_DX11_RIBBON_GROUP_SIZE - 1) / TFX_IMPL_DX11_RIBBON_GROUP_SIZE, 1, 1);
	}

	//The vertex and index buffers are read by the draw next, so they can't stay bound for writing
	ID3D11ShaderResourceView *no_views[4] = {};
	ID3D11UnorderedAccessView *no_unordered_views[2] = {};
	device_context->CSSetShaderResources(6, 4, no_views);
	device_context->CSSetUnorderedAccessViews(0, 2, no_unordered_views, NULL);
	device_context->CSSetShader(NULL, NULL, 0);
}

void tfx_impl_dx11_Compute(tfx_impl_dx11_t *renderer, const tfx_impl_uniforms_t *uniforms) {
	if (renderer->particle_dispatch_count > 0) {
		tfx_impl_dx11__simulate_gpu_particles(renderer);
	}
	if (renderer->playback_sprite_count > 0) {
		tfx_impl_dx11__play_back_sprites(renderer);
	}
	if (renderer->playback_ribbon_count > 0) {
		tfx_impl_dx11__play_back_ribbons(renderer);
	}
	if (tfx_impl_dx11__ribbon_set_active(&renderer->live_ribbons)) {
		tfx_impl_dx11__build_ribbon_geometry(renderer, &renderer->live_ribbons, uniforms);
	}
	if (tfx_impl_dx11__ribbon_set_active(&renderer->baked_ribbons)) {
		tfx_impl_dx11__build_ribbon_geometry(renderer, &renderer->baked_ribbons, uniforms);
	}
}

//Live and baked ribbons each have their own shape frames and color ramps, the atlas is shared
static void tfx_impl_dx11__render_ribbons(tfx_impl_dx11_t *renderer, const tfx_impl_dx11_ribbon_set_t *set, ID3D11ShaderResourceView *image_data_view, ID3D11ShaderResourceView *color_ramps_view) {
	ID3D11DeviceContext *device_context = renderer->device_context;
	device_context->IASetIndexBuffer(set->indices.buffer, DXGI_FORMAT_R32_UINT, 0);
	ID3D11ShaderResourceView *vertex_views[5] = { set->emitters.view, tfx_impl_dx11__ribbon_graph_lookups(renderer, set), set->instances.view, NULL, set->vertices.view };
	device_context->VSSetShader(renderer->ribbon_vertex_shader, NULL, 0);
	device_context->VSSetConstantBuffers(2, 1, &renderer->ribbon_constants);
	device_context->VSSetShaderResources(2, 1, &image_data_view);
	device_context->VSSetShaderResources(6, 5, vertex_views);
	device_context->PSSetShader(renderer->ribbon_pixel_shader, NULL, 0);
	device_context->PSSetShaderResources(5, 1, &color_ramps_view);

	tfx_ribbon_dispatch_t dispatch = tfx_impl_dx11__create_ribbon_dispatch(set);
	while (tfx_impl_dx11__next_ribbon_dispatch(set, &dispatch)) {
		tfx_impl_dx11__write_constants(device_context, renderer->ribbon_constants, tfx_GetRibbonDispatchGlobals(&dispatch), sizeof(tfx_ribbon_bucket_globals_t));
		device_context->DrawIndexed(dispatch.index_count, dispatch.index_offset, 0);
	}

	//Unbound so the next compute pass can write them
	ID3D11ShaderResourceView *no_views[5] = {};
	device_context->VSSetShaderResources(6, 5, no_views);
	device_context->IASetIndexBuffer(NULL, DXGI_FORMAT_R32_UINT, 0);
}

static void tfx_impl_dx11__draw_instance_range(tfx_impl_dx11_t *renderer, tfxU32 start, tfxU32 count) {
	if (count == 0) {
		return;
	}
	tfx_impl_dx11_draw_constants_t draw_constants = {};
	draw_constants.instance_offset = start;
	tfx_impl_dx11__write_constants(renderer->device_context, renderer->draw_constants, &draw_constants, sizeof(draw_constants));
	renderer->device_context->DrawInstanced(6, count, 0, 0);
}

void tfx_impl_dx11_Render(tfx_impl_dx11_t *renderer, const tfx_impl_uniforms_t *uniforms) {
	if (!renderer->color_ramps_view || (renderer->instance_count == 0 && renderer->playback_sprite_count == 0
		&& !tfx_impl_dx11__ribbon_set_active(&renderer->live_ribbons) && !tfx_impl_dx11__ribbon_set_active(&renderer->baked_ribbons))) {
		return;
	}
	ID3D11DeviceContext *device_context = renderer->device_context;
	tfx_impl_dx11__write_constants(device_context, renderer->frame_constants, uniforms, sizeof(tfx_impl_uniforms_t));

	//Quad corners come from SV_VertexID and the instances and ribbon vertices are read from buffers, so there's no input layout or vertex buffer
	device_context->IASetInputLayout(NULL);
	device_context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
	const float blend_factor[4] = { 0.f, 0.f, 0.f, 0.f };
	device_context->OMSetBlendState(renderer->blend_state, blend_factor, 0xffffffff);
	device_context->OMSetDepthStencilState(renderer->depth_stencil_state, 0);
	device_context->RSSetState(renderer->rasterizer_state);
	ID3D11Buffer *constant_buffers[2] = { renderer->frame_constants, renderer->draw_constants };
	ID3D11ShaderResourceView *pixel_views[2] = { renderer->shape_atlas_view, renderer->color_ramps_view };
	device_context->VSSetConstantBuffers(0, 2, constant_buffers);
	device_context->PSSetShaderResources(4, 2, pixel_views);
	device_context->PSSetSamplers(0, 1, &renderer->sampler_state);

	if (renderer->instance_count > 0) {
		ID3D11ShaderResourceView *vertex_views[4] = {
			renderer->instance_buffers[renderer->current_instance_buffer].view,
			renderer->instance_buffers[renderer->current_instance_buffer ^ 1].view,
			renderer->image_data.view,
			renderer->emitter_properties.view,
		};
		device_context->VSSetShader(renderer->billboard_vertex_shader, NULL, 0);
		device_context->VSSetShaderResources(0, 4, vertex_views);
		device_context->PSSetShader(renderer->billboard_pixel_shader, NULL, 0);
		for (tfxU32 layer = 0; layer != tfxLAYERS; ++layer) {
			tfx_impl_dx11__draw_instance_range(renderer, renderer->layer_start[layer], renderer->layer_count[layer]);
			tfx_impl_dx11__draw_instance_range(renderer, renderer->gpu_layer_start[layer], renderer->gpu_layer_count[layer]);
		}
		//Unbound so the GPU particle kernel can write the current buffer next frame
		ID3D11ShaderResourceView *no_views[2] = {};
		device_context->VSSetShaderResources(0, 2, no_views);
	}
	if (tfx_impl_dx11__ribbon_set_active(&renderer->live_ribbons)) {
		tfx_impl_dx11__render_ribbons(renderer, &renderer->live_ribbons, renderer->image_data.view, renderer->color_ramps_view);
	}
	if (renderer->playback_sprite_count > 0) {
		//Baked sprites come with their own shape frames, emitter properties and color ramps, the atlas is shared
		ID3D11ShaderResourceView *vertex_views[4] = { renderer->playback_instances.view, NULL, renderer->baked_image_data.view, renderer->baked_emitter_properties.view };
		device_context->VSSetShader(renderer->billboard_static_vertex_shader, NULL, 0);
		device_context->VSSetShaderResources(0, 4, vertex_views);
		device_context->PSSetShader(renderer->billboard_pixel_shader, NULL, 0);
		device_context->PSSetShaderResources(5, 1, &renderer->baked_color_ramps_view);
		tfx_impl_dx11__draw_instance_range(renderer, 0, renderer->playback_sprite_count);
		ID3D11ShaderResourceView *no_views[1] = {};
		device_context->VSSetShaderResources(0, 1, no_views);
	}
	if (tfx_impl_dx11__ribbon_set_active(&renderer->baked_ribbons)) {
		tfx_impl_dx11__render_ribbons(renderer, &renderer->baked_ribbons, renderer->baked_image_data.view, renderer->baked_color_ramps_view);
	}
}

static void tfx_impl_dx11__release_ribbon_set(tfx_impl_dx11_ribbon_set_t *set) {
	tfx_impl_dx11__release_buffer(&set->segments);
	tfx_impl_dx11__release_buffer(&set->instances);
	tfx_impl_dx11__release_buffer(&set->emitters);
	tfx_impl_dx11__release_buffer(&set->vertices);
	tfx_impl_dx11__release_buffer(&set->indices);
	tfx_impl_dx11__release_buffer(&set->graph_lookups);
}

void tfx_impl_dx11_Shutdown(tfx_impl_dx11_t *renderer) {
	for (int index = 0; index != 2; ++index) {
		tfx_impl_dx11__release_buffer(&renderer->instance_buffers[index]);
	}
	tfx_impl_dx11__release_buffer(&renderer->image_data);
	tfx_impl_dx11__release_buffer(&renderer->emitter_properties);
	tfx_impl_dx11__release_buffer(&renderer->graph_lookups);
	tfx_impl_dx11__release_ribbon_set(&renderer->live_ribbons);
	tfx_impl_dx11__release_ribbon_set(&renderer->baked_ribbons);
	tfx_impl_dx11__release_buffer(&renderer->particle_pool);
	tfx_impl_dx11__release_buffer(&renderer->particle_upload);
	TFX_IMPL_DX11_RELEASE(renderer->particle_push);
	TFX_IMPL_DX11_RELEASE(renderer->particle_kernel);
	tfx_impl_dx11__release_buffer(&renderer->baked_image_data);
	tfx_impl_dx11__release_buffer(&renderer->baked_emitter_properties);
	tfx_impl_dx11__release_buffer(&renderer->baked_sprites);
	tfx_impl_dx11__release_buffer(&renderer->animation_offsets);
	tfx_impl_dx11__release_buffer(&renderer->animation_instances);
	tfx_impl_dx11__release_buffer(&renderer->playback_instances);
	tfx_impl_dx11__release_buffer(&renderer->baked_ribbon_data);
	tfx_impl_dx11__release_buffer(&renderer->baked_ribbon_properties);
	tfx_impl_dx11__release_buffer(&renderer->animation_ribbon_offsets);
	TFX_IMPL_DX11_RELEASE(renderer->ribbon_data_constants);
	TFX_IMPL_DX11_RELEASE(renderer->ribbon_playback_shader);
	TFX_IMPL_DX11_RELEASE(renderer->baked_color_ramps_view);
	TFX_IMPL_DX11_RELEASE(renderer->baked_color_ramps);
	TFX_IMPL_DX11_RELEASE(renderer->sprite_data_constants);
	TFX_IMPL_DX11_RELEASE(renderer->billboard_static_vertex_shader);
	TFX_IMPL_DX11_RELEASE(renderer->sprite_playback_shader);
	TFX_IMPL_DX11_RELEASE(renderer->color_ramps_view);
	TFX_IMPL_DX11_RELEASE(renderer->color_ramps);
	TFX_IMPL_DX11_RELEASE(renderer->shape_atlas_view);
	TFX_IMPL_DX11_RELEASE(renderer->shape_atlas);
	TFX_IMPL_DX11_RELEASE(renderer->ribbon_constants);
	TFX_IMPL_DX11_RELEASE(renderer->draw_constants);
	TFX_IMPL_DX11_RELEASE(renderer->frame_constants);
	TFX_IMPL_DX11_RELEASE(renderer->sampler_state);
	TFX_IMPL_DX11_RELEASE(renderer->rasterizer_state);
	TFX_IMPL_DX11_RELEASE(renderer->depth_stencil_state);
	TFX_IMPL_DX11_RELEASE(renderer->blend_state);
	TFX_IMPL_DX11_RELEASE(renderer->ribbon_pixel_shader);
	TFX_IMPL_DX11_RELEASE(renderer->ribbon_vertex_shader);
	TFX_IMPL_DX11_RELEASE(renderer->ribbon_compute_shader);
	TFX_IMPL_DX11_RELEASE(renderer->billboard_pixel_shader);
	TFX_IMPL_DX11_RELEASE(renderer->billboard_vertex_shader);
	TFX_IMPL_DX11_RELEASE(renderer->device_context);
	TFX_IMPL_DX11_RELEASE(renderer->device);
}
