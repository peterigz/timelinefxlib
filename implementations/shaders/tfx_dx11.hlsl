/*
	TimelineFX Direct3D 11 shaders (shader model 5.0), ported from the GLSL reference shaders in ../../Shaders.
	Variable names follow the GLSL so the two can be compared side by side.

	Table of contents, search for the keyword to jump to the section:
	[Shader_resources]		Constant buffers, buffers, textures and the C struct mirrors
	[Instance_unpacking]	Reading tfx_instance_t out of a raw buffer
	[Billboard_vertex]		timelinefx3d.vert
	[Billboard_pixel]		timelinefx.frag
	[Ribbon_resources]		The ribbon structs, buffers and the per bucket constant buffer
	[Graph_lookup]			Sampling the graph lookups, shared by the ribbon compute and vertex shaders
	[Ribbon_noise]			Value and gradient noise for ribbon displacement
	[Ribbon_compute]		ribbons.comp, builds the ribbon vertex and index buffers
	[Ribbon_vertex]			ribbon_3d.vert
	[Ribbon_pixel]			ribbon.frag
	[Sprite_playback]		sprite_data_playback.comp, interpolates baked sprite data into billboard instances
	[Static_billboard]		static_billboard.vert, draws the playback's instances with no interpolation of its own
	[Ribbon_playback]		ribbon_data_playback.comp, interpolates baked ribbon data into the ribbons ribbon_compute builds
*/

//--- Shader_resources

//tfx_impl_uniforms_t. Matrices are uploaded column major, which is HLSL's default packing, so mul(matrix, vector) matches GLSL's matrix * vector
cbuffer tfx_frame_constants : register(b0) {
	float4x4 view;
	float4x4 projection;
	float2 screen_size;
	uint milliseconds;
	float timer_lerp;
	float update_time;
	float3 eye_position;
};

cbuffer tfx_draw_constants : register(b1) {
	uint instance_offset;			//SV_InstanceID doesn't include StartInstanceLocation in D3D11, so the draw's first instance comes in here
	uint previous_index_offset;
	uint2 draw_padding;
};

//tfx_gpu_image_data_t, 48 bytes
struct tfx_image_data {
	float4 uv;
	uint2 uv_packed;
	float2 image_size;
	uint texture_array_index;		//The atlas layer. The reference shaders keep a bindless image index in the upper 16 bits, here they're 0
	float animation_frames;
	float max_radius;
	float padding;
};

//tfx_gpu_particle_properties_t, 32 bytes
struct tfx_emitter_properties {
	float2 image_handle;
	uint color_ramp_indexes;		//[Row of color ramp bitmap, texture array]
	uint flags;						//Billboard alignment in the lowest 2 bits
	uint start_frame_index;
	float animation_frames;
	uint2 padding;
};

ByteAddressBuffer current_instances : register(t0);
ByteAddressBuffer previous_instances : register(t1);
StructuredBuffer<tfx_image_data> image_data : register(t2);
StructuredBuffer<tfx_emitter_properties> emitter_properties : register(t3);
Texture2DArray particle_images : register(t4);
Texture2DArray color_ramps : register(t5);
SamplerState particle_sampler : register(s0);

//--- Instance_unpacking

//tfx_instance_t unpacked the way the reference pipeline's vertex formats read it
struct tfx_instance {
	float4 position;				//Stretch in w
	float4 quaternion;				//r16g16b16a16_snorm
	float2 size;					//r16g16_sscaled
	float3 alignment;				//r8g8b8_snorm
	float2 intensity_gradient_map;	//r16g16_sscaled
	float3 curved_alpha_life;		//r8g8b8_unorm
	uint indexes;					//[gpu property index << 16, capture flag << 15, image frame]
	uint captured_index;			//The sprite's index in the previous buffer
};

float2 tfx_unpack_snorm16x2(uint packed) {
	int2 value = int2(int(packed << 16) >> 16, int(packed) >> 16);
	return max(float2(value) / 32767.0, -1.0);
}

float2 tfx_unpack_sscaled16x2(uint packed) {
	return float2(int(packed << 16) >> 16, int(packed) >> 16);
}

float3 tfx_unpack_snorm8x3(uint packed) {
	int3 value = int3(int(packed << 24) >> 24, int(packed << 16) >> 24, int(packed << 8) >> 24);
	return max(float3(value) / 127.0, -1.0);
}

float3 tfx_unpack_unorm8x3(uint packed) {
	return float3(packed & 0xFF, (packed >> 8) & 0xFF, (packed >> 16) & 0xFF) / 255.0;
}

tfx_instance tfx_load_instance(ByteAddressBuffer buffer, uint index) {
	uint address = index * 48;
	uint4 first = buffer.Load4(address);
	uint4 second = buffer.Load4(address + 16);
	uint4 third = buffer.Load4(address + 32);
	tfx_instance instance;
	instance.position = asfloat(first);
	instance.quaternion = float4(tfx_unpack_snorm16x2(second.x), tfx_unpack_snorm16x2(second.y));
	instance.size = tfx_unpack_sscaled16x2(second.z);
	instance.alignment = tfx_unpack_snorm8x3(second.w);
	instance.intensity_gradient_map = tfx_unpack_sscaled16x2(third.x);
	instance.curved_alpha_life = tfx_unpack_unorm8x3(third.y);
	instance.indexes = third.z;
	instance.captured_index = third.w;
	return instance;
}

//GLSL's mat3(a, b, c) takes columns, HLSL's float3x3(a, b, c) takes rows
float3x3 tfx_matrix_from_columns(float3 first, float3 second, float3 third) {
	return transpose(float3x3(first, second, third));
}

//--- Billboard_vertex

//#define NO_LERP

static const uint quad_indexes[6] = { 0, 1, 2, 2, 1, 3 };
//An epsilon avoids nans with the cross product
static const float3 up = float3(0, 1, 0.00001);
static const float size_max_value = 64.0 / 32767.0;
static const float intensity_max_value = 128.0 / 32767.0;

struct tfx_billboard_output {
	float4 position : SV_Position;
	float3 tex_coord : TEXCOORD0;
	nointerpolation int3 texture_indexes : TEXCOORD1;
	float4 intensity_curved_alpha_map : TEXCOORD2;
};

//GLSL fills m[column][row], HLSL indexes m[row][column]
float3x3 QuaternionToRotationMatrix(float4 q) {
	float xx = q.x * q.x; float xy = q.x * q.y; float xz = q.x * q.z; float xw = q.x * q.w;
	float yy = q.y * q.y; float yz = q.y * q.z; float yw = q.y * q.w;
	float zz = q.z * q.z; float zw = q.z * q.w;
	float3x3 m;
	m[0][0] = 1.0 - 2.0 * (yy + zz); m[1][0] = 2.0 * (xy + zw);       m[2][0] = 2.0 * (xz - yw);
	m[0][1] = 2.0 * (xy - zw);       m[1][1] = 1.0 - 2.0 * (xx + zz); m[2][1] = 2.0 * (yz + xw);
	m[0][2] = 2.0 * (xz + yw);       m[1][2] = 2.0 * (yz - xw);       m[2][2] = 1.0 - 2.0 * (xx + yy);
	return m;
}

tfx_billboard_output billboard_vertex(uint vertex_id : SV_VertexID, uint instance_id : SV_InstanceID) {
	tfx_instance instance = tfx_load_instance(current_instances, instance_id + instance_offset);
	uint emitter_index = (instance.indexes >> 16) & 0xFFFFu;
	tfx_emitter_properties props = emitter_properties[emitter_index];

	float2 size = instance.size * size_max_value;
	float2 handle = props.image_handle;
	float4 quaternion = normalize(instance.quaternion);

	uint prev_index = (instance.captured_index & 0x0FFFFFFF) + previous_index_offset;

	//Not interpolated on the first frame of a particle or when it loops back to the start of a path, and it's hidden that frame
	float interpolate_is_active = float((instance.indexes & 0x00008000) >> 15);

#ifdef NO_LERP
	float2 lerped_size = size;
	float3 lerped_position = instance.position.xyz;
	float4 lerped_quaternion = quaternion;
	float3 motion = float3(0, 0, 0);
#else
	tfx_instance previous = tfx_load_instance(previous_instances, prev_index);
	float2 lerped_size = lerp(previous.size * size_max_value, size, timer_lerp) * interpolate_is_active;
	float3 lerped_position = lerp(previous.position.xyz, instance.position.xyz, timer_lerp);
	float4 lerped_quaternion = lerp(normalize(previous.quaternion), quaternion, timer_lerp);
	float3 motion = instance.position.xyz - previous.position.xyz;
#endif
	motion.z += 0.000001;
	bool has_alignment = dot(instance.alignment, instance.alignment) > 0;
	float travel_distance = length(motion);
	float stretch_factor = instance.position.w * interpolate_is_active;

	motion = normalize(motion);
	float3 final_alignment = has_alignment ? instance.alignment : motion;

	//01 = camera facing off, 11 = align to the vector with no camera facing, 10 = face the camera and roll to the vector
	uint align_flags = props.flags & 0x3u;
	bool billboarding = (align_flags & 0x1u) > 0u;
	bool align_type = align_flags == 0x3u;
	bool vector_align = align_flags == 0x2u;

	//GLSL's v * inverse(mat3(view)), the view rotation is orthonormal so its inverse transpose is itself
	float3x3 view_rotation = (float3x3)view;
	float3 camera_relative_aligment = mul(view_rotation, motion);

	uint index = quad_indexes[vertex_id];

	float3x3 base_spin_mat = QuaternionToRotationMatrix(lerped_quaternion);

	float3 normal = normalize(final_alignment);
	float3 tangent = normalize(cross(up, normal));
	if (length(tangent) < 0.001) {
		tangent = normalize(cross(float3(1.0, 0.0, 0.0), normal));
	}
	float3 bitangent = normalize(cross(normal, tangent));
	float3x3 align_mat = tfx_matrix_from_columns(tangent, bitangent, normal);

	float3 camera_relative_alignment = mul(view_rotation, final_alignment);
	float2 align_xy = normalize(camera_relative_alignment.xy + float2(0.00001, 0.00001));
	float2 view_up_xy = float2(0, 1);
	float c_vec = dot(align_xy, view_up_xy);
	float s_vec = align_xy.x * view_up_xy.y - align_xy.y * view_up_xy.x;
	float final_c_roll = c_vec;
	float final_s_roll = -s_vec;
	float3x3 vector_align_roll_mat = tfx_matrix_from_columns(float3(final_c_roll, final_s_roll, 0.0), float3(-final_s_roll, final_c_roll, 0.0), float3(0.0, 0.0, 1.0));

	float3x3 final_rot_mat = base_spin_mat;
	if (align_type) {
		final_rot_mat = mul(align_mat, base_spin_mat);
	} else if (vector_align) {
		final_rot_mat = mul(vector_align_roll_mat, base_spin_mat);
	}

	float3 identity_bounds[4] = {
		float3(lerped_size.x * (0 - handle.x), -lerped_size.y * (0 - handle.y), 0),
		float3(lerped_size.x * (1 - handle.x), -lerped_size.y * (0 - handle.y), 0),
		float3(lerped_size.x * (0 - handle.x), -lerped_size.y * (1 - handle.y), 0),
		float3(lerped_size.x * (1 - handle.x), -lerped_size.y * (1 - handle.y), 0)
	};
	float3 pos = mul(final_rot_mat, identity_bounds[index]);

	//Stretch follows the direction of travel. Each vertex moves along it by its own projection over the quad's projected half extent, which keeps it continuous
	float3 stretch_direction = billboarding ? motion : normalize(camera_relative_aligment);
	float vertex_offset_magnitude = 0.5 * travel_distance * stretch_factor;
	float projected_half_extent = 0.0;
	for (int corner = 0; corner < 4; ++corner) {
		projected_half_extent = max(projected_half_extent, abs(dot(mul(final_rot_mat, identity_bounds[corner]), stretch_direction)));
	}
	float projected_position = dot(pos, stretch_direction);
	float3 potential_offset = stretch_direction * (projected_position / max(projected_half_extent, 0.00001)) * vertex_offset_magnitude;
	bool apply_stretch = travel_distance > 0.00001 && stretch_factor > 0.0;
	pos += apply_stretch ? potential_offset : float3(0.0, 0.0, 0.0);

	//GLSL replaces modelView's rotation with identity for camera facing quads, which is the corner offset added in view space
	float4 view_position = mul(view, float4(lerped_position, 1.0));
	float3 view_offset = billboarding ? mul(view_rotation, pos) : pos;

	tfx_billboard_output output;
	output.position = mul(projection, float4(view_position.xyz + view_offset, view_position.w));

	//The sprite holds the animation frame only, the emitter's properties hold where its image starts
	uint image_index = props.start_frame_index + (instance.indexes & 0x00001FFF);
	float4 uv = image_data[image_index].uv;
	float2 uvs[4] = { uv.xy, uv.zy, uv.xw, uv.zw };

	int life = int(instance.curved_alpha_life.z * 255);
	int ramp_y = int(props.color_ramp_indexes & 0xFFu);
	int ramp_array = int((props.color_ramp_indexes >> 8) & 0xFFu);
	output.texture_indexes = int3(ramp_y, ramp_array, life);
	uint packed_texture_index = image_data[image_index].texture_array_index;
	output.tex_coord = float3(uvs[index], float(packed_texture_index & 0xFFFFu));
	output.intensity_curved_alpha_map = float4(instance.intensity_gradient_map.x * intensity_max_value, instance.curved_alpha_life.x, instance.curved_alpha_life.y, instance.intensity_gradient_map.y * intensity_max_value);
	return output;
}

//--- Billboard_pixel

//Writes premultiplied alpha, so blend with ONE, INV_SRC_ALPHA
float4 billboard_pixel(tfx_billboard_output input) : SV_Target {
	float4 texel = particle_images.Sample(particle_sampler, input.tex_coord);
	float lookup = saturate(texel.r * input.intensity_curved_alpha_map.w);
	int ramp_x = int(lookup * 255);
	float4 ramp_texel = color_ramps.Load(int4(ramp_x, input.texture_indexes.x, input.texture_indexes.y, 0));
	ramp_texel *= input.intensity_curved_alpha_map.x;
	ramp_texel.a = min(1, ramp_texel.a);
	float curved_alpha = 1 - smoothstep(texel.a * input.intensity_curved_alpha_map.z, texel.a, 1 - input.intensity_curved_alpha_map.y);
	float4 out_color;
	out_color.rgb = texel.rgb * ramp_texel.rgb * curved_alpha * texel.a;
	out_color.a = texel.a * ramp_texel.a * curved_alpha;
	return out_color;
}

//--- Ribbon_resources

//tfx_ribbon_bucket_globals_t, set by the library per bucket. The bindless indexes in it are unused here
cbuffer tfx_ribbon_constants : register(b2) {
	float4 camera_position;
	uint segment_count;
	uint tessellation;
	uint index_offset;
	uint vertex_offset;
	uint ribbon_count;
	uint ribbon_offset;
	uint segment_offset;
	uint uniform_index;
	uint emitters_index;
	uint graphs_index;
	uint ribbons_index;
	uint ribbon_segments_index;
	uint vertexes_index;
	uint indexes_index;
	uint image_data_index;
	uint sampler_index;
	uint particle_texture_index;
	uint color_ramp_texture_index;
	float ribbon_lerp_time;			//pc.lerp in the GLSL
	float time;
	float ndc_offset_x;
	float ndc_offset_y;
	uint emitter_offset;
};

static const float tfx_packed_scale_unpack = 128.0 / 32767.0;

static const int tfx_intensity_offset                  = 0;
static const int tfx_alpha_sharpness_offset            = 1;
static const int tfx_curved_alpha_offset               = 2;
static const int tfx_gradient_map_offset               = 3;
static const int tfx_width_offset                      = 4;
static const int tfx_scale_offset                      = 5;
static const int tfx_uv_offset_y_offset                = 6;
static const int tfx_uv_scale_y_offset                 = 7;
static const int tfx_clip_offset_offset                = 8;
static const int tfx_clip_size_offset                  = 9;
static const int tfx_morph_amount_offset               = 10;
static const int tfx_noise_amount_offset               = 11;
static const int tfx_noise_scroll_offset               = 12;
static const int tfx_lag_amount_offset                 = 13;
static const int tfx_intensity_overlength_offset       = 14;
static const int tfx_alpha_sharpness_overlength_offset = 15;
static const int tfx_curved_alpha_overlength_offset    = 16;
static const int tfx_gradient_map_overlength_offset    = 17;
static const int tfx_width_overlength_offset           = 18;
static const int tfx_fixed_angle_offset                = 19;
static const int tfx_morph_bias_offset                 = 20;
static const int tfx_noise_envelope_offset             = 21;
static const int tfx_lag_profile_offset                = 22;

static const uint tfx_active_mask = 0x1u;
static const uint tfx_relative_mask = 0x2u;
static const uint tfx_ribbon_flags_shift = 24u;
static const uint tfx_ribbon_start_index_mask = 0x00FFFFFFu;

//Hash axes the ribbon noise phase offset is drawn on, they must differ from each other and from every other axis in the system
static const uint tfx_ribbon_noise_phase_x_axis = 0x27d4eb2du;
static const uint tfx_ribbon_noise_phase_y_axis = 0x165667b1u;
static const uint tfx_ribbon_noise_phase_z_axis = 0x9e3779b1u;

static const uint tfx_invalid_index = 0xFFFFFFFFu;

static const uint tfx_angle_type_face_camera = 0;
static const uint tfx_angle_type_fixed_angle = 1;
static const uint tfx_angle_type_frenet_serret = 2;

#define tfx_lag_spine_samples 16

//tfx_ribbon_t, 48 bytes
struct tfx_ribbon {
	float3 position;
	uint age_width;					//Normalised age as a unorm16 in the low half, width as a half float in the high half
	uint2 quaternion;				//16-bit snorm quaternion: .x = X|Y, .y = Z|W
	uint start_index_flags;			//Segment start index in the low 24 bits, flags in the top 8
	uint emitter_index_phase_seed;	//Emitter index in the low 16 bits, phase offset seed in the high 16
	uint texture_indexes;
	uint intensity_gradient_map;
	uint curved_alpha;
	uint scale;						//The parent effect's overall scale as a half float in the low 16 bits
};

//tfx_gpu_ribbon_emitter_t, 416 bytes. Structured buffers pack tightly, so the C padding is spelled out
struct tfx_ribbon_emitter {
	uint2 quaternion;
	uint lookup_offset;
	uint angle_type;
	float3 position;
	float age;
	float3 captured_position;
	uint morph_segment_start_index;
	float noise_frequency;
	float noise_speed;
	float noise_phase_range;
	float noise_lock_rate;
	float3 fixed_angle_normal;
	uint sample_count_noise;		//[sample count in the low 24 bits, noise algorithm and octaves in the top 8]
	uint start_frame_index;
	float lag_time;
	float lag_span;
	float lag_spine_position_x[tfx_lag_spine_samples];
	float lag_spine_position_y[tfx_lag_spine_samples];
	float lag_spine_position_z[tfx_lag_spine_samples];
	uint lag_spine_padding;
	uint2 lag_spine_quaternion[tfx_lag_spine_samples];
};

//tfx_ribbon_vertex_t, 32 bytes
struct tfx_ribbon_vertex {
	float3 position;
	uint segment_index;
	float2 uv_offset_scale;
	uint ribbon_index;
	uint padding;
};

//tfx_gpu_graph_data_t
struct tfx_gpu_graph_data {
	float4 node_data;
	float4 oscillator;
	int easing_type;
	int flags;
};

StructuredBuffer<tfx_ribbon_emitter> ribbon_emitters : register(t6);
//Raw so the GPU particle kernel can read the same buffer
ByteAddressBuffer graph_lookups : register(t7);
StructuredBuffer<tfx_ribbon> ribbon_instances : register(t8);
StructuredBuffer<float4> ribbon_segments : register(t9);
StructuredBuffer<tfx_ribbon_vertex> ribbon_vertices : register(t10);
RWStructuredBuffer<tfx_ribbon_vertex> ribbon_vertices_out : register(u0);
RWByteAddressBuffer ribbon_indices_out : register(u1);

//--- Graph_lookup

tfx_gpu_graph_data tfx_load_graph(uint index) {
	uint address = index * 48;
	uint4 third = graph_lookups.Load4(address + 32);
	tfx_gpu_graph_data graph;
	graph.node_data = asfloat(graph_lookups.Load4(address));
	graph.oscillator = asfloat(graph_lookups.Load4(address + 16));
	graph.easing_type = asint(third.x);
	graph.flags = asint(third.y);
	return graph;
}

float tfx_easing_in_cubic(float t) { return t * t * t; }
float tfx_easing_out_cubic(float t) { t = 1 - t; return 1 - (t * t * t); }
float tfx_easing_in_out_cubic(float t) { float t3 = -2 * t + 2; t3 = t3 * t3 * t3; return t < 0.5 ? 4 * t * t * t : 1 - t3 * 0.5; }
float tfx_easing_smoothstep(float t) { return t * t * (3.0 - 2.0 * t); }
float tfx_easing_out_in(float t) { return t < 0.5 ? t * (2.0 - 2.0 * t) : 2.0 * (t - 0.5) * (t - 0.5) + 0.5; }

float tfx_apply_easing(int easing_type, float t) {
	switch (easing_type) {
		case 0:  return 0.0;
		case 4:  return tfx_easing_in_cubic(t);
		case 5:  return tfx_easing_out_cubic(t);
		case 6:  return tfx_easing_in_out_cubic(t);
		case 16: return t;
		case 17: return tfx_easing_smoothstep(t);
		case 18: return tfx_easing_out_in(t);
		default: return t;
	}
}

float tfx_do_bezier_lerp(float t, float first_node, float curve1, float curve2, float second_node) {
	float u = 1 - t;
	float w1 = u * u * u;
	float w2 = 3 * u * u * t;
	float w3 = 3 * u * t * t;
	float w4 = t * t * t;
	return w1 * first_node + w2 * curve1 + w3 * curve2 + w4 * second_node;
}

float tfx_do_linear_lerp(float t, float first_node, float second_node) {
	return first_node + (second_node - first_node) * t;
}

static const int tfx_use_bezier_sampling_flag = 1;
static const int tfx_enable_oscillator_flag = 2;
static const float tfx_two_pi = 6.283185307;

float tfx_graph_lookup(tfx_gpu_graph_data graph, float t) {
	float eased_t = tfx_apply_easing(graph.easing_type, t);
	float result;
	if ((graph.flags & tfx_use_bezier_sampling_flag) != 0) {
		result = tfx_do_bezier_lerp(eased_t, graph.node_data.x, graph.node_data.y, graph.node_data.z, graph.node_data.w);
	} else {
		result = tfx_do_linear_lerp(eased_t, graph.node_data.x, graph.node_data.w);
	}
	if ((graph.flags & tfx_enable_oscillator_flag) != 0) {
		//oscillator: .x = frequency, .y = amplitude, .z = offset_x shifts the phase, .w = offset_y lifts the multiplier
		float oscillator_sin = 0.5 + sin((eased_t + graph.oscillator.z) * graph.oscillator.x * tfx_two_pi) * graph.oscillator.y;
		result = result * (oscillator_sin + graph.oscillator.w);
	}
	return result;
}

float tfx_emitter_graph(uint lookup_offset, int graph_offset, float t) {
	return tfx_graph_lookup(tfx_load_graph(graph_offset + lookup_offset), t);
}

//--- Ribbon_noise

//Murmur3 finalizer, a direct port of tfx__seedgen_u32 so the GPU field matches the CPU one's character
uint tfx_seedgen_u32(uint h) {
	h ^= h >> 15;
	h *= 0x85ebca6bu;
	h ^= h >> 13;
	h *= 0xc2b2ae35u;
	h ^= h >> 16;
	return h;
}

//Port of tfx__white_unit, [-1, 1]: the signed reinterpret is what centres it
float tfx_white_unit(uint seed, uint axis) {
	return (float(asint(tfx_seedgen_u32(seed ^ axis))) / 4294967295.0) * 2.0;
}

//One corner of the value noise lattice, a single hash unpacked into all three components
float3 tfx_value_noise_corner(int3 lattice) {
	uint h = tfx_seedgen_u32(uint(lattice.x + lattice.y) ^ uint(lattice.z));
	return float3(float((h >> 20) & 0x3FFu), float((h >> 10) & 0x3FFu), float(h & 0x3FFu)) * (1.0 / 511.0) - 1.0;
}

static const int3 tfx_noise_primes = int3(1619, 31337, 6971);

//Trilinear value noise, a port of tfx__wide_value_noise_3d
float3 tfx_value_noise_3d(float3 position) {
	float3 floor_position = floor(position);
	float3 fade = position - floor_position;
	fade = fade * fade * (3.0 - 2.0 * fade);

	int3 lattice0 = int3(floor_position) * tfx_noise_primes;
	int3 lattice1 = lattice0 + tfx_noise_primes;

	float3 near_edge = lerp(tfx_value_noise_corner(int3(lattice0.x, lattice0.y, lattice0.z)),
	                        tfx_value_noise_corner(int3(lattice1.x, lattice0.y, lattice0.z)), fade.x);
	float3 far_edge = lerp(tfx_value_noise_corner(int3(lattice0.x, lattice1.y, lattice0.z)),
	                       tfx_value_noise_corner(int3(lattice1.x, lattice1.y, lattice0.z)), fade.x);
	float3 near_slice = lerp(near_edge, far_edge, fade.y);

	near_edge = lerp(tfx_value_noise_corner(int3(lattice0.x, lattice0.y, lattice1.z)),
	                 tfx_value_noise_corner(int3(lattice1.x, lattice0.y, lattice1.z)), fade.x);
	far_edge = lerp(tfx_value_noise_corner(int3(lattice0.x, lattice1.y, lattice1.z)),
	                tfx_value_noise_corner(int3(lattice1.x, lattice1.y, lattice1.z)), fade.x);
	float3 far_slice = lerp(near_edge, far_edge, fade.y);

	return lerp(near_slice, far_slice, fade.z);
}

//Perlin's gradient select: one of the cube edge vectors picked by the low bits of the hash
float tfx_gradient_dot(uint hash, float3 distance_to_corner) {
	uint h = hash & 15u;
	float u = h < 8u ? distance_to_corner.x : distance_to_corner.y;
	float v = h < 4u ? distance_to_corner.y : ((h == 12u || h == 14u) ? distance_to_corner.x : distance_to_corner.z);
	return ((h & 1u) == 0u ? u : -u) + ((h & 2u) == 0u ? v : -v);
}

//One scalar field of gradient noise with a quintic fade
float tfx_gradient_noise_field(float3 position, uint seed) {
	float3 floor_position = floor(position);
	float3 relative = position - floor_position;
	float3 fade = relative * relative * relative * (relative * (relative * 6.0 - 15.0) + 10.0);

	int3 lattice0 = int3(floor_position) * tfx_noise_primes;
	int3 lattice1 = lattice0 + tfx_noise_primes;

	float near_lower = lerp(tfx_gradient_dot(tfx_seedgen_u32(uint(lattice0.x + lattice0.y) ^ uint(lattice0.z) ^ seed), relative),
	                        tfx_gradient_dot(tfx_seedgen_u32(uint(lattice1.x + lattice0.y) ^ uint(lattice0.z) ^ seed), relative - float3(1.0, 0.0, 0.0)), fade.x);
	float far_lower = lerp(tfx_gradient_dot(tfx_seedgen_u32(uint(lattice0.x + lattice1.y) ^ uint(lattice0.z) ^ seed), relative - float3(0.0, 1.0, 0.0)),
	                       tfx_gradient_dot(tfx_seedgen_u32(uint(lattice1.x + lattice1.y) ^ uint(lattice0.z) ^ seed), relative - float3(1.0, 1.0, 0.0)), fade.x);
	float lower_slice = lerp(near_lower, far_lower, fade.y);

	float near_upper = lerp(tfx_gradient_dot(tfx_seedgen_u32(uint(lattice0.x + lattice0.y) ^ uint(lattice1.z) ^ seed), relative - float3(0.0, 0.0, 1.0)),
	                        tfx_gradient_dot(tfx_seedgen_u32(uint(lattice1.x + lattice0.y) ^ uint(lattice1.z) ^ seed), relative - float3(1.0, 0.0, 1.0)), fade.x);
	float far_upper = lerp(tfx_gradient_dot(tfx_seedgen_u32(uint(lattice0.x + lattice1.y) ^ uint(lattice1.z) ^ seed), relative - float3(0.0, 1.0, 1.0)),
	                       tfx_gradient_dot(tfx_seedgen_u32(uint(lattice1.x + lattice1.y) ^ uint(lattice1.z) ^ seed), relative - float3(1.0, 1.0, 1.0)), fade.x);
	float upper_slice = lerp(near_upper, far_upper, fade.y);

	return lerp(lower_slice, upper_slice, fade.z);
}

//Measured so gradient noise lands at the same rms as the value noise
static const float tfx_gradient_noise_scale = 1.397;

float3 tfx_gradient_noise_3d(float3 position) {
	return float3(tfx_gradient_noise_field(position, 0u),
	              tfx_gradient_noise_field(position, 0x9e3779b9u),
	              tfx_gradient_noise_field(position, 0x85ebca6bu)) * tfx_gradient_noise_scale;
}

static const uint tfx_noise_algorithm_simplex = 2u;
static const uint tfx_noise_algorithm_value = 4u;

//Octaves at doubling frequency and halving amplitude, normalised by the root of the summed squares
float3 tfx_noise_octaves(float3 position, uint algorithm, uint octaves) {
	float3 sum = float3(0.0, 0.0, 0.0);
	float amplitude = 1.0;
	float frequency = 1.0;
	float normalisation = 0.0;
	for (uint octave = 0u; octave < octaves; ++octave) {
		sum += amplitude * (algorithm == tfx_noise_algorithm_value
			? tfx_value_noise_3d(position * frequency)
			: tfx_gradient_noise_3d(position * frequency));
		normalisation += amplitude * amplitude;
		amplitude *= 0.5;
		frequency *= 2.0;
	}
	return sum * rsqrt(max(normalisation, 1e-8));
}

//--- Ribbon_compute

//256 rather than the GLSL's 1024, which fxc warns is over its register budget. Any multiple of 32 keeps ribbons warp aligned
#define tfx_ribbon_group_size 256

uint pack_segment_data(uint segment_index, uint ribbon_lerped, uint tessellation_count, uint tessellation_index, uint side) {
	return (tessellation_count << 29) | (tessellation_index << 26) | (side << 25) | (ribbon_lerped << 12) | segment_index;
}

float3 rotate_with_quaternion(float3 position, float4 quaternion) {
	return position + 2.0 * cross(cross(position, quaternion.xyz) + quaternion.w * position, quaternion.xyz);
}

//With a clip size only the clip window + clip offset section of the path is sampled
float3 sample_ribbon_path(int segment, uint base_index, uint sample_count, float clip_offset, float clip_size) {
	float segment_lerp = float(segment) / float(segment_count - 1);
	float path_t = clip_offset + segment_lerp * clip_size;
	float stored = clamp(path_t, 0.0, 1.0) * float(sample_count - 1);
	uint i0 = min(uint(stored), sample_count - 2);
	float fraction = clamp(stored - float(i0), 0.0, 1.0);
	return lerp(ribbon_segments[base_index + i0].xyz, ribbon_segments[base_index + i0 + 1].xyz, fraction);
}

float3 safe_normalize(float3 v) {
	float length_squared = dot(v, v);
	return v * rsqrt(length_squared + 1e-16);
}

float4 multiply_quaternions(float4 q1, float4 q2) {
	return float4(
		q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
		q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x,
		q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w,
		q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z
	);
}

//Path, then morph, then noise, for the centre segment and both neighbours so every angle type's frame comes from displaced geometry
float3 sample_ribbon_shape(int segment, uint base_index, uint morph_base_index, uint sample_count,
                           float clip_offset, float clip_size, float morph_amount,
                           tfx_gpu_graph_data morph_bias_graph, tfx_gpu_graph_data noise_envelope_graph,
                           float noise_amplitude, float noise_frequency, float3 noise_sample_offset, uint noise_packed) {
	float3 position = sample_ribbon_path(segment, base_index, sample_count, clip_offset, clip_size);
	float segment_lerp = float(segment) / float(segment_count - 1);
	if (morph_base_index != tfx_invalid_index) {
		float morph = clamp(morph_amount + tfx_graph_lookup(morph_bias_graph, clamp(segment_lerp, 0.0, 1.0)), 0.0, 1.0);
		float3 morph_position = sample_ribbon_path(segment, morph_base_index, sample_count, clip_offset, clip_size);
		position = lerp(position, morph_position, morph);
	}

	if (noise_amplitude == 0.0 || noise_packed == 0u) {
		return position;
	}
	float envelope = tfx_graph_lookup(noise_envelope_graph, clamp(segment_lerp, 0.0, 1.0));
	float3 sample_position = noise_sample_offset + float3(segment_lerp * noise_frequency, 0.0, 0.0);
	return position + noise_amplitude * envelope * tfx_noise_octaves(sample_position, noise_packed & 0xFu, (noise_packed >> 4) & 0xFu);
}

struct tfx_lag_transform {
	float3 position;
	float4 quaternion;
};

tfx_lag_transform tfx_emitter_transform_at(uint emitter_index, float lag_time, float lag_span, float depth) {
	float frame_length = lag_span - lag_time;
	float distance_back = (1.0 - ribbon_lerp_time) * frame_length + depth * lag_time;
	float spine_position = clamp(distance_back / lag_span, 0.0, 1.0) * float(tfx_lag_spine_samples - 1);
	uint first_index = min(uint(spine_position), uint(tfx_lag_spine_samples - 2));
	uint second_index = first_index + 1u;
	float fraction = spine_position - float(first_index);

	float3 first_position = float3(ribbon_emitters[emitter_index].lag_spine_position_x[first_index],
	                               ribbon_emitters[emitter_index].lag_spine_position_y[first_index],
	                               ribbon_emitters[emitter_index].lag_spine_position_z[first_index]);
	float3 second_position = float3(ribbon_emitters[emitter_index].lag_spine_position_x[second_index],
	                                ribbon_emitters[emitter_index].lag_spine_position_y[second_index],
	                                ribbon_emitters[emitter_index].lag_spine_position_z[second_index]);

	uint2 first_packed = ribbon_emitters[emitter_index].lag_spine_quaternion[first_index];
	uint2 second_packed = ribbon_emitters[emitter_index].lag_spine_quaternion[second_index];
	float4 first_quaternion = float4(tfx_unpack_snorm16x2(first_packed.x), tfx_unpack_snorm16x2(first_packed.y));
	float4 second_quaternion = float4(tfx_unpack_snorm16x2(second_packed.x), tfx_unpack_snorm16x2(second_packed.y));
	second_quaternion *= dot(first_quaternion, second_quaternion) < 0.0 ? -1.0 : 1.0;

	tfx_lag_transform result;
	result.position = lerp(first_position, second_position, fraction);
	result.quaternion = normalize(lerp(first_quaternion, second_quaternion, fraction));
	return result;
}

float3 tfx_ribbon_world_position(float3 segment_position, float ribbon_scale, float effect_scale, float4 ribbon_quaternion, float3 ribbon_position, float4 emitter_quaternion, float3 emitter_position) {
	float3 position = rotate_with_quaternion(segment_position * ribbon_scale * effect_scale, ribbon_quaternion) + ribbon_position;
	return rotate_with_quaternion(position, emitter_quaternion) + emitter_position;
}

//Segment counts are always a multiple of 32, so a ribbon never shares a warp with another and the angle type and easing branches stay uniform
[numthreads(tfx_ribbon_group_size, 1, 1)]
void ribbon_compute(uint3 dispatch_thread_id : SV_DispatchThreadID) {
	uint global_index = dispatch_thread_id.x;
	if (global_index >= segment_count * ribbon_count) return;

	uint local_ribbon_index = global_index / segment_count;
	uint global_ribbon_index = local_ribbon_index + ribbon_offset;
	uint global_offset = segment_offset;

	tfx_ribbon ribbon = ribbon_instances[global_ribbon_index];
	float ribbon_age = float(ribbon.age_width & 0xFFFFu) * (1.0 / 65535.0);
	float ribbon_width = f16tof32(ribbon.age_width >> 16);
	float effect_scale = f16tof32(ribbon.scale & 0xFFFFu);
	uint ribbon_flags = ribbon.start_index_flags >> tfx_ribbon_flags_shift;
	uint phase_seed = ribbon.emitter_index_phase_seed >> 16;
	uint emitter_index = (ribbon.emitter_index_phase_seed & 0xFFFFu) + emitter_offset;
	//Fields are read from the buffer one at a time rather than copying the 416 byte emitter
	uint lookup_offset = ribbon_emitters[emitter_index].lookup_offset;
	uint angle_type = ribbon_emitters[emitter_index].angle_type;
	uint local_segment_index = global_index % segment_count;
	float ribbon_lerp = float(local_segment_index) / float(segment_count - 1);

	uint base_segment_index = global_offset + (ribbon.start_index_flags & tfx_ribbon_start_index_mask);
	//A unique id per segment for the packed vertex attribute
	uint global_segment_index = base_segment_index + local_segment_index;

	float clip_offset = tfx_emitter_graph(lookup_offset, tfx_clip_offset_offset, ribbon_age);
	float clip_size = tfx_emitter_graph(lookup_offset, tfx_clip_size_offset, ribbon_age);

	float morph_amount = tfx_emitter_graph(lookup_offset, tfx_morph_amount_offset, ribbon_age);
	tfx_gpu_graph_data morph_bias_graph = tfx_load_graph(tfx_morph_bias_offset + lookup_offset);
	//The morph block lives in the same bucket as the static one so it rebases with the same offset
	uint morph_segment_start_index = ribbon_emitters[emitter_index].morph_segment_start_index;
	uint morph_base_segment_index = morph_segment_start_index == tfx_invalid_index ? tfx_invalid_index : global_offset + morph_segment_start_index;

	tfx_gpu_graph_data noise_envelope_graph = tfx_load_graph(tfx_noise_envelope_offset + lookup_offset);
	float noise_amplitude = tfx_emitter_graph(lookup_offset, tfx_noise_amount_offset, ribbon_age);
	float noise_lock_rate = ribbon_emitters[emitter_index].noise_lock_rate;
	float noise_time = ribbon_age * ribbon_emitters[emitter_index].noise_speed;
	if (noise_lock_rate > 0.0) {
		noise_time = floor(noise_time * noise_lock_rate) / noise_lock_rate;
	}
	float noise_scroll = tfx_emitter_graph(lookup_offset, tfx_noise_scroll_offset, ribbon_age);
	float3 phase_offset = float3(tfx_white_unit(phase_seed, tfx_ribbon_noise_phase_x_axis),
	                             tfx_white_unit(phase_seed, tfx_ribbon_noise_phase_y_axis),
	                             tfx_white_unit(phase_seed, tfx_ribbon_noise_phase_z_axis))
	                      * ribbon_emitters[emitter_index].noise_phase_range;
	float3 noise_sample_offset = phase_offset + float3(noise_scroll, noise_time, 0.0);
	float noise_frequency = ribbon_emitters[emitter_index].noise_frequency;

	bool relative = (ribbon_flags & tfx_relative_mask) > 0;

	uint2 emitter_packed_quaternion = ribbon_emitters[emitter_index].quaternion;
	float4 ribbon_quaternion = float4(tfx_unpack_snorm16x2(ribbon.quaternion.x), tfx_unpack_snorm16x2(ribbon.quaternion.y));
	float4 emitter_unpacked_quaternion = float4(tfx_unpack_snorm16x2(emitter_packed_quaternion.x), tfx_unpack_snorm16x2(emitter_packed_quaternion.y));
	float4 emitter_quaternion = relative ? emitter_unpacked_quaternion : float4(0, 0, 0, 1.0);
	float3 emitter_position = relative ? lerp(ribbon_emitters[emitter_index].captured_position, ribbon_emitters[emitter_index].position, ribbon_lerp_time) : float3(0.0, 0.0, 0.0);
	int neighbour_segment_index = (local_segment_index < segment_count - 1) ? int(local_segment_index) + 1 : int(local_segment_index) - 1;
	float4 neighbour_emitter_quaternion = emitter_quaternion;
	float3 neighbour_emitter_position = emitter_position;
	bool lagging = false;
	float lag_time = ribbon_emitters[emitter_index].lag_time;
	if (relative && lag_time > 0.0) {
		float lag_amount = tfx_emitter_graph(lookup_offset, tfx_lag_amount_offset, ribbon_age);
		if (lag_amount > 0.0) {
			float lag_span = ribbon_emitters[emitter_index].lag_span;
			tfx_gpu_graph_data lag_profile_graph = tfx_load_graph(tfx_lag_profile_offset + lookup_offset);
			float neighbour_lerp = clamp(float(neighbour_segment_index) / float(segment_count - 1), 0.0, 1.0);
			tfx_lag_transform centre_transform = tfx_emitter_transform_at(emitter_index, lag_time, lag_span, lag_amount * tfx_graph_lookup(lag_profile_graph, ribbon_lerp));
			tfx_lag_transform neighbour_transform = tfx_emitter_transform_at(emitter_index, lag_time, lag_span, lag_amount * tfx_graph_lookup(lag_profile_graph, neighbour_lerp));
			emitter_quaternion = centre_transform.quaternion;
			emitter_position = centre_transform.position;
			neighbour_emitter_quaternion = neighbour_transform.quaternion;
			neighbour_emitter_position = neighbour_transform.position;
			lagging = true;
		}
	}
	uint sample_count_noise = ribbon_emitters[emitter_index].sample_count_noise;
	uint sample_count = sample_count_noise & 0xFFFFFFu;
	uint noise_packed = sample_count_noise >> 24;
	float3 current_segment_pos = sample_ribbon_shape(int(local_segment_index), base_segment_index, morph_base_segment_index, sample_count, clip_offset, clip_size, morph_amount, morph_bias_graph, noise_envelope_graph, noise_amplitude, noise_frequency, noise_sample_offset, noise_packed);

	float ribbon_scale = tfx_emitter_graph(lookup_offset, tfx_scale_offset, ribbon_age);

	float3 next_segment_pos = sample_ribbon_shape(int(local_segment_index) + 1, base_segment_index, morph_base_segment_index, sample_count, clip_offset, clip_size, morph_amount, morph_bias_graph, noise_envelope_graph, noise_amplitude, noise_frequency, noise_sample_offset, noise_packed);
	float3 prev_segment_pos = sample_ribbon_shape(int(local_segment_index) - 1, base_segment_index, morph_base_segment_index, sample_count, clip_offset, clip_size, morph_amount, morph_bias_graph, noise_envelope_graph, noise_amplitude, noise_frequency, noise_sample_offset, noise_packed);

	bool has_next = local_segment_index < segment_count - 1;
	float3 current_pos = tfx_ribbon_world_position(current_segment_pos, ribbon_scale, effect_scale, ribbon_quaternion, ribbon.position, emitter_quaternion, emitter_position);
	//The direction along the ribbon in world space, towards the next segment, or from the previous one at the end
	float3 neighbour_pos = tfx_ribbon_world_position(has_next ? next_segment_pos : prev_segment_pos, ribbon_scale, effect_scale, ribbon_quaternion, ribbon.position, neighbour_emitter_quaternion, neighbour_emitter_position);
	float3 world_forward = safe_normalize(has_next ? neighbour_pos - current_pos : current_pos - neighbour_pos);
	float3 right;

	if (angle_type == tfx_angle_type_face_camera) {
		float3 to_camera = normalize(camera_position.xyz - current_pos);
		float3 bisector = normalize(to_camera + world_forward);
		right = safe_normalize(cross(world_forward, bisector));

	} else if (angle_type == tfx_angle_type_fixed_angle) {
		float fixed_angle = tfx_emitter_graph(lookup_offset, tfx_fixed_angle_offset, ribbon_lerp);
		float c = cos(fixed_angle);
		float s = sin(fixed_angle);
		float3 fixed_angle_normal = ribbon_emitters[emitter_index].fixed_angle_normal;
		if (lagging) {
			float3 up_vector = rotate_with_quaternion(fixed_angle_normal, multiply_quaternions(ribbon_quaternion, emitter_quaternion));
			right = safe_normalize(cross(world_forward, up_vector));
			right = safe_normalize(right * c + cross(world_forward, right) * s + world_forward * dot(world_forward, right) * (1.0 - c));
		} else {
			//Built in the ribbon's own space and rotated into the world after
			float3 forward = safe_normalize(has_next ? next_segment_pos - current_segment_pos : current_segment_pos - prev_segment_pos);
			right = safe_normalize(cross(forward, fixed_angle_normal));
			right = right * c + cross(forward, right) * s + forward * dot(forward, right) * (1.0 - c);
			right = safe_normalize(rotate_with_quaternion(right, multiply_quaternions(ribbon_quaternion, emitter_quaternion)));
		}

	} else {
		//Frenet-serret frame
		float3 up_vector = float3(0.0, 1.0, 0.0);
		if (abs(dot(world_forward, up_vector)) > 0.999) {
			up_vector = float3(1.0, 0.0, 0.0);
		}
		float3 binormal = safe_normalize(cross(world_forward, up_vector));
		right = safe_normalize(cross(binormal, world_forward));
	}

	uint ribbon_lerp_packed = uint(ribbon_lerp * 2047);
	float ribbon_is_active = float(ribbon_flags & tfx_active_mask);

	//The ribbon's scale, which is what the positions above are already scaled by
	float width_overtime_and_length = tfx_emitter_graph(lookup_offset, tfx_width_overlength_offset, ribbon_lerp);
	width_overtime_and_length = width_overtime_and_length * effect_scale * tfx_emitter_graph(lookup_offset, tfx_width_offset, ribbon_age) * ribbon_is_active;

	float uv_offset = tfx_emitter_graph(lookup_offset, tfx_uv_offset_y_offset, ribbon_age);
	float uv_scale = tfx_emitter_graph(lookup_offset, tfx_uv_scale_y_offset, ribbon_age);

	uint vertices_per_segment = (tessellation + 1) * 2;
	uint vertices_per_ribbon = segment_count * vertices_per_segment;
	uint vertex_base_index = (local_ribbon_index * vertices_per_ribbon) + (local_segment_index * vertices_per_segment) + vertex_offset;

	//Tessellation is the same for the whole dispatch, so this never diverges
	if (tessellation > 0) {
		for (uint i = 0; i <= tessellation; i++) {
			float t = float(i) / float(tessellation);
			float width_factor = lerp(.5, 0, t) * width_overtime_and_length;

			tfx_ribbon_vertex vertex;
			vertex.position = current_pos + right * ribbon_width * width_factor;
			vertex.segment_index = pack_segment_data(global_segment_index, ribbon_lerp_packed, tessellation, i, 0);
			vertex.uv_offset_scale = float2(uv_offset, uv_scale);
			vertex.ribbon_index = global_ribbon_index;
			vertex.padding = 0;
			ribbon_vertices_out[vertex_base_index + i] = vertex;

			width_factor = lerp(0, -0.5, t) * width_overtime_and_length;
			vertex.position = current_pos + right * ribbon_width * width_factor;
			vertex.segment_index = pack_segment_data(global_segment_index, ribbon_lerp_packed, tessellation, i, 1);
			ribbon_vertices_out[vertex_base_index + (tessellation + 1) + i] = vertex;
		}

		if (has_next) {
			uint indices_per_segment = tessellation * 12;
			uint indices_per_ribbon = segment_count * indices_per_segment - indices_per_segment;
			uint index_base = index_offset + (local_ribbon_index * indices_per_ribbon) + (local_segment_index * indices_per_segment);
			uint next_base_vertex = vertex_base_index + vertices_per_segment;

			for (uint quad = 0; quad < tessellation; quad++) {
				uint current_quad_index = index_base + (quad * 12);
				uint v0 = vertex_base_index + quad;
				uint v1 = vertex_base_index + quad + 1;
				uint v2 = next_base_vertex + quad;
				uint v3 = next_base_vertex + quad + 1;
				ribbon_indices_out.Store3((current_quad_index + 0) * 4, uint3(v0, v1, v2));
				ribbon_indices_out.Store3((current_quad_index + 3) * 4, uint3(v2, v1, v3));

				uint next_t = tessellation + 1;
				v0 += next_t;
				v1 += next_t;
				v2 += next_t;
				v3 += next_t;
				ribbon_indices_out.Store3((current_quad_index + 6) * 4, uint3(v0, v1, v2));
				ribbon_indices_out.Store3((current_quad_index + 9) * 4, uint3(v2, v1, v3));
			}
		}
	} else {
		//Ported as it is: unlike the tessellated path this one leaves uv_offset_scale and ribbon_index unwritten
		float width_factor = ribbon_width * 0.5;
		ribbon_vertices_out[vertex_base_index].position = current_pos + right * width_factor;
		ribbon_vertices_out[vertex_base_index].segment_index = global_segment_index;
		ribbon_vertices_out[vertex_base_index + 1].position = current_pos - right * width_factor;
		ribbon_vertices_out[vertex_base_index + 1].segment_index = global_segment_index;

		if (has_next) {
			uint index_base = index_offset + (global_segment_index * 6);
			uint v0 = vertex_base_index;
			uint v1 = vertex_base_index + 2;
			uint v2 = vertex_base_index + 1;
			uint v3 = vertex_base_index + 3;
			ribbon_indices_out.Store3((index_base + 0) * 4, uint3(v0, v1, v2));
			ribbon_indices_out.Store3((index_base + 3) * 4, uint3(v2, v1, v3));
		}
	}
}

//--- Ribbon_vertex

struct tfx_ribbon_output {
	float4 position : SV_Position;
	float3 tex_coord : TEXCOORD0;
	nointerpolation int3 texture_indexes : TEXCOORD1;
	float4 intensity_curved_alpha_map : TEXCOORD2;
	float ribbon_across : TEXCOORD3;	//0 to 1 across the width, its screen space derivative is the ribbon's width in pixels
};

//DrawIndexed's SV_VertexID is the index buffer value, which ribbon_compute wrote as an absolute vertex index
tfx_ribbon_output ribbon_vertex(uint vertex_id : SV_VertexID) {
	tfx_ribbon_vertex vertex = ribbon_vertices[vertex_id];
	tfx_ribbon_output output;
	output.position = mul(projection, mul(view, float4(vertex.position, 1.0)));
	output.position.x += ndc_offset_x * output.position.w;
	output.position.y += ndc_offset_y * output.position.w;
	uint segment_index = vertex.segment_index;
	float ribbon_position = float((segment_index & 0x007FF000) >> 12) / 2047.0;
	tfx_ribbon ribbon = ribbon_instances[vertex.ribbon_index];
	uint emitter_index = (ribbon.emitter_index_phase_seed & 0xFFFFu) + emitter_offset;
	uint lookup_offset = ribbon_emitters[emitter_index].lookup_offset;
	//The ribbon holds the animation frame, its emitter holds where the image starts
	uint image_index = ribbon_emitters[emitter_index].start_frame_index + (ribbon.texture_indexes & 0x00001FFF);

	float tessellation_count = float((segment_index & 0xE0000000) >> 29);
	float tessellation_index = float((segment_index & 0x1C000000) >> 26);
	float side = float((segment_index & 0x02000000) >> 25);
	float t = tessellation_index / tessellation_count;
	float4 uv = image_data[image_index].uv;
	float texture_uv_width = uv.z - uv.x;

	float sampled_gradient_map = tfx_emitter_graph(lookup_offset, tfx_gradient_map_overlength_offset, ribbon_position);
	float sampled_intensity = tfx_emitter_graph(lookup_offset, tfx_intensity_overlength_offset, ribbon_position);
	float sampled_curved_alpha = tfx_emitter_graph(lookup_offset, tfx_curved_alpha_overlength_offset, ribbon_position);
	float sampled_alpha_sharpness = tfx_emitter_graph(lookup_offset, tfx_alpha_sharpness_overlength_offset, ribbon_position);

	//Ping pong the u coordinate on every other wrap, so a texture repeating along the ribbon has no flipped seam where it wraps
	float scaled_position = ribbon_position * vertex.uv_offset_scale.y + vertex.uv_offset_scale.x;
	float wrap_index = floor(scaled_position);
	float wrap_fraction = frac(scaled_position);
	//GLSL's mod, which unlike fmod is never negative
	if (wrap_index - 2.0 * floor(wrap_index * 0.5) == 1.0) {
		wrap_fraction = 1.0 - wrap_fraction;
	}

	float uv_x = wrap_fraction * texture_uv_width + uv.x;
	float ribbon_across = t * 0.5 + (side * 0.5);
	output.ribbon_across = ribbon_across;
	float uv_y = ribbon_across * (uv.w - uv.y) + uv.y;
	float2 ribbon_intensity_gradient_map = tfx_unpack_sscaled16x2(ribbon.intensity_gradient_map) * tfx_packed_scale_unpack;
	float2 ribbon_curved_alpha = tfx_unpack_sscaled16x2(ribbon.curved_alpha) * tfx_packed_scale_unpack;

	int life = int(ribbon_position * 255);
	uint packed_texture_index = image_data[image_index].texture_array_index;
	output.tex_coord = float3(uv_x, uv_y, float(packed_texture_index & 0xFFFFu));
	output.texture_indexes = int3((ribbon.texture_indexes & 0xFF000000) >> 24, (ribbon.texture_indexes & 0x00FF0000) >> 16, life);
	output.intensity_curved_alpha_map = float4(sampled_intensity * ribbon_intensity_gradient_map.x, sampled_curved_alpha * ribbon_curved_alpha.x,
		sampled_alpha_sharpness * ribbon_curved_alpha.y, sampled_gradient_map * ribbon_intensity_gradient_map.y);
	return output;
}

//--- Ribbon_pixel

//Below this many pixels across a ribbon starts dropping fragments and breaks into dashes, so it's faded out as it approaches it
static const float tfx_min_ribbon_pixels = 1.5;

float4 ribbon_pixel(tfx_ribbon_output input) : SV_Target {
	float4 texel = particle_images.Sample(particle_sampler, input.tex_coord);
	float lookup = saturate(texel.r * input.intensity_curved_alpha_map.w);
	int ramp_x = int(lookup * 255);
	float4 ramp_texel = color_ramps.Load(int4(ramp_x, input.texture_indexes.x, input.texture_indexes.y, 0));
	ramp_texel *= input.intensity_curved_alpha_map.x;
	ramp_texel.a = min(1, ramp_texel.a);
	float curved_alpha = 1 - smoothstep(texel.a * input.intensity_curved_alpha_map.z, texel.a, 1 - input.intensity_curved_alpha_map.y);
	float across_per_pixel = fwidth(input.ribbon_across);
	float ribbon_pixels = across_per_pixel > 0.0 ? 1.0 / across_per_pixel : tfx_min_ribbon_pixels;
	float thin_fade = smoothstep(0.0, tfx_min_ribbon_pixels, ribbon_pixels);
	float4 out_color;
	out_color.rgb = texel.rgb * ramp_texel.rgb * curved_alpha * texel.a * thin_fade;
	out_color.a = texel.a * ramp_texel.a * curved_alpha * thin_fade;
	return out_color;
}

//--- Sprite_playback

//tfx_sprite_data_push_t. The bindless indexes in it are unused here
cbuffer tfx_sprite_data_constants : register(b3) {
	uint animation_instances_total;
	uint billboards_total;
	uint animated_shapes;
	uint offsets_index;
	uint animation_instances_index;
	uint billboards_index;
	uint sprite_data_index;
	uint image_data_index_unused;
	uint emitter_properties_index;
	uint bounding_boxes_index;
};

ByteAddressBuffer baked_sprites : register(t11);			//tfx_GetSpriteDataBufferPointer, tfx_instance_t
ByteAddressBuffer animation_offsets : register(t12);		//tfx_GetOffsetsBufferPointer, a running sprite total per animation instance
ByteAddressBuffer animation_instances : register(t13);		//tfx_GetAnimationInstancesBufferPointer, tfx_animation_instance_t
RWByteAddressBuffer playback_instances_out : register(u2);

#define tfx_invalid 0xFFFFFFFF
#define tfx_reverse_animation_flag 8192
static const float size_max_divisor = 32767.0 / 64.0;
static const float intensity_max_divisor = 32767.0 / 128.0;

//tfx_animation_instance_t, 64 bytes: position, scale | sprite_count, ribbon_count, frame_count, offset_into_sprite_data |
//offset_into_ribbon_data, info_index, current_time, animation_length_in_time | tween, flags, padding
struct tfx_animation_instance {
	float3 position;
	float scale;
	uint frame_count;
	uint offset_into_sprite_data;
	uint offset_into_ribbon_data;
	float current_time;
	float animation_length_in_time;
};

tfx_animation_instance tfx_load_animation_instance(uint index) {
	uint address = index * 64;
	uint4 first = animation_instances.Load4(address);
	uint4 second = animation_instances.Load4(address + 16);
	uint4 third = animation_instances.Load4(address + 32);
	tfx_animation_instance instance;
	instance.position = asfloat(first.xyz);
	instance.scale = asfloat(first.w);
	instance.frame_count = second.z;
	instance.offset_into_sprite_data = second.w;
	instance.offset_into_ribbon_data = third.x;
	instance.current_time = asfloat(third.z);
	instance.animation_length_in_time = asfloat(third.w);
	return instance;
}

//The GLSL builtins the playback shader packs and unpacks with
float4 unpack_unorm4x8(uint packed) {
	return float4(packed & 0xFF, (packed >> 8) & 0xFF, (packed >> 16) & 0xFF, packed >> 24) / 255.0;
}

uint pack_unorm4x8(float4 value) {
	uint4 scaled = uint4(round(saturate(value) * 255.0));
	return scaled.x | (scaled.y << 8) | (scaled.z << 16) | (scaled.w << 24);
}

uint pack_snorm4x8(float4 value) {
	int4 scaled = int4(round(clamp(value, -1.0, 1.0) * 127.0));
	return (uint(scaled.x) & 0xFF) | ((uint(scaled.y) & 0xFF) << 8) | ((uint(scaled.z) & 0xFF) << 16) | ((uint(scaled.w) & 0xFF) << 24);
}

uint pack_snorm2x16(float2 value) {
	int2 scaled = int2(round(clamp(value, -1.0, 1.0) * 32767.0));
	return (uint(scaled.x) & 0xFFFF) | (uint(scaled.y) << 16);
}

//The GLSL reads these halves unsigned, so it's kept that way
float2 unpack_16bit_2_scale(uint packed, float max_value) {
	return float2(packed & 0xFFFF, (packed >> 16) & 0xFFFF) * max_value;
}

[numthreads(128, 1, 1)]
void sprite_playback(uint3 dispatch_thread_id : SV_DispatchThreadID) {
	uint index = dispatch_thread_id.x;
	if (index >= billboards_total) return;

	//The offsets are a running total of sprites per instance, so count the instances whose total is at or below this sprite
	uint instance_index = 0;
	for (uint i = 0; i < animation_instances_total; ++i) {
		instance_index += uint(animation_offsets.Load(instance_index * 4) <= index);
	}
	uint sprite_index_offset = instance_index == 0 ? 0 : animation_offsets.Load((instance_index - 1) * 4);
	tfx_animation_instance animation = tfx_load_animation_instance(instance_index);

	//Frame 0 is 0 to 1, frame 1 is 1 to 2 and so on, and the fraction is how far between this frame and the next
	float frame_time = (animation.current_time / animation.animation_length_in_time) * float(animation.frame_count);
	float lerp_amount = frame_time - floor(frame_time);

	uint sprite_data_index = (index - sprite_index_offset) + animation.offset_into_sprite_data;
	uint4 current_first = baked_sprites.Load4(sprite_data_index * 48);
	uint4 current_second = baked_sprites.Load4(sprite_data_index * 48 + 16);
	uint4 current_third = baked_sprites.Load4(sprite_data_index * 48 + 32);
	uint captured_index = current_third.w;
	//A sprite with no previous frame, its first, interpolates with itself
	if (captured_index == tfx_invalid) captured_index = sprite_data_index;
	uint4 captured_first = baked_sprites.Load4(captured_index * 48);
	uint4 captured_second = baked_sprites.Load4(captured_index * 48 + 16);
	uint4 captured_third = baked_sprites.Load4(captured_index * 48 + 32);

	uint indexes = current_third.z;
	uint property_index = (indexes & 0xFFFF0000) >> 16;
	//lerp_offset is how far into the frame the sprite spawns, 1 once it's lived a whole frame. Not spawned yet means no size
	float4 curved_alpha_life = unpack_unorm4x8(current_third.y);
	//The spawn record holds the lerp_offset and it's the captured sprite of the next frame, the first one that interpolates from it
	float lerp_offset = unpack_unorm4x8(captured_third.y).w;
	float spawn_scale = ((lerp_offset >= (1 - lerp_amount)) && (captured_index != sprite_data_index)) ? 1 : 0;
	//The captured state is the sprite at its spawn time, so interpolate from the spawn point rather than the start of the frame
	lerp_amount = saturate((lerp_amount - (1 - lerp_offset)) / max(lerp_offset, 1.0 / 255.0));

	float3 lerped_position = lerp(asfloat(captured_first.xyz), asfloat(current_first.xyz), lerp_amount) * animation.scale;
	float4 quaternion = float4(tfx_unpack_snorm16x2(current_second.x), tfx_unpack_snorm16x2(current_second.y));
	float4 captured_quaternion = float4(tfx_unpack_snorm16x2(captured_second.x), tfx_unpack_snorm16x2(captured_second.y));
	float4 lerped_quaternion = lerp(captured_quaternion, quaternion, lerp_amount);
	float2 scale = unpack_16bit_2_scale(current_second.z, size_max_value);
	float2 captured_scale = unpack_16bit_2_scale(captured_second.z, size_max_value);
	float2 lerped_scale = lerp(captured_scale, scale, lerp_amount) * spawn_scale * animation.scale;
	float3 lerped_alignment = lerp(tfx_unpack_snorm8x3(captured_second.w), tfx_unpack_snorm8x3(current_second.w), lerp_amount);
	float lerped_stretch = lerp(asfloat(captured_first.w), asfloat(current_first.w), lerp_amount) * spawn_scale;
	tfx_emitter_properties properties = emitter_properties[property_index];
	uint lerped_frame = properties.start_frame_index;
	uint out_indexes;
	if (animated_shapes != 0) {
		//The same for every thread, so it never diverges. Interpolates the shape's frame, allowing for it wrapping round at the end
		float reverse_direction = float((properties.flags & tfx_reverse_animation_flag) > 1);
		float current_frame = float(indexes & 0x00001FFF);
		float prev_frame = float(captured_third.z & 0x00001FFF);
		float frames = properties.animation_frames;
		float wrap = (1 - reverse_direction) * float(current_frame < prev_frame);
		float reverse_wrap = reverse_direction * float(prev_frame < current_frame);
		float frame = lerp(prev_frame + frames * reverse_wrap, current_frame + frames * wrap, lerp_amount);
		//GLSL's mod
		lerped_frame += uint(frame - frames * floor(frame / frames));
		out_indexes = (indexes & 0xFFFFE000) | lerped_frame;
	} else {
		out_indexes = indexes + lerped_frame;
	}
	float2 intensity_gradient_map = unpack_16bit_2_scale(current_third.x, intensity_max_value);
	float2 captured_intensity_gradient_map = unpack_16bit_2_scale(captured_third.x, intensity_max_value);
	float2 lerped_intensity_gradient_map = lerp(captured_intensity_gradient_map, intensity_gradient_map, lerp_amount);
	float3 captured_curved_alpha_life = unpack_unorm4x8(captured_third.y).xyz;
	float3 lerped_curved_alpha_life = lerp(captured_curved_alpha_life, curved_alpha_life.xyz, lerp_amount);

	lerped_position += animation.position;

	uint address = index * 48;
	playback_instances_out.Store4(address, uint4(asuint(lerped_position), asuint(lerped_stretch)));
	playback_instances_out.Store4(address + 16, uint4(
		pack_snorm2x16(lerped_quaternion.xy),
		pack_snorm2x16(lerped_quaternion.zw),
		(uint(lerped_scale.y * size_max_divisor) << 16) + uint(lerped_scale.x * size_max_divisor),
		pack_snorm4x8(float4(lerped_alignment, 0))));
	playback_instances_out.Store4(address + 32, uint4(
		(uint(lerped_intensity_gradient_map.y * intensity_max_divisor) << 16) + uint(lerped_intensity_gradient_map.x * intensity_max_divisor),
		pack_unorm4x8(float4(lerped_curved_alpha_life, 0)),
		out_indexes,
		index));
}

//--- Static_billboard

//The playback has already interpolated, resolved the image frame (start frame included) and baked the stretch, so nothing reads a previous buffer here
tfx_billboard_output billboard_static_vertex(uint vertex_id : SV_VertexID, uint instance_id : SV_InstanceID) {
	tfx_instance instance = tfx_load_instance(current_instances, instance_id + instance_offset);
	uint emitter_index = (instance.indexes >> 16) & 0xFFFFu;
	tfx_emitter_properties props = emitter_properties[emitter_index];

	float2 size = instance.size * size_max_value;
	float2 handle = props.image_handle;
	float4 quaternion = normalize(instance.quaternion);

	uint align_flags = props.flags & 0x3u;
	bool billboarding = (align_flags & 0x1u) > 0u;
	bool align_type = align_flags == 0x3u;
	bool vector_align = align_flags == 0x2u;

	uint index = quad_indexes[vertex_id];
	float3x3 base_spin_mat = QuaternionToRotationMatrix(quaternion);

	float3 normal = normalize(instance.alignment);
	float3 tangent = normalize(cross(up, normal));
	if (length(tangent) < 0.001) {
		tangent = normalize(cross(float3(1.0, 0.0, 0.0), normal));
	}
	float3 bitangent = normalize(cross(normal, tangent));
	float3x3 align_mat = tfx_matrix_from_columns(tangent, bitangent, normal);

	float3x3 view_rotation = (float3x3)view;
	float3 camera_relative_alignment = mul(view_rotation, instance.alignment);
	float2 align_xy = normalize(camera_relative_alignment.xy + float2(0.00001, 0.00001));
	float2 view_up_xy = float2(0, 1);
	float c_vec = dot(align_xy, view_up_xy);
	float s_vec = align_xy.x * view_up_xy.y - align_xy.y * view_up_xy.x;
	float final_c_roll = c_vec;
	float final_s_roll = -s_vec;
	float3x3 vector_align_roll_mat = tfx_matrix_from_columns(float3(final_c_roll, final_s_roll, 0.0), float3(-final_s_roll, final_c_roll, 0.0), float3(0.0, 0.0, 1.0));

	float3x3 final_rot_mat = base_spin_mat;
	if (align_type) {
		final_rot_mat = mul(align_mat, base_spin_mat);
	} else if (vector_align) {
		final_rot_mat = mul(vector_align_roll_mat, base_spin_mat);
	}

	float3 identity_bounds[4] = {
		float3(size.x * (0 - handle.x), -size.y * (0 - handle.y), 0),
		float3(size.x * (1 - handle.x), -size.y * (0 - handle.y), 0),
		float3(size.x * (0 - handle.x), -size.y * (1 - handle.y), 0),
		float3(size.x * (1 - handle.x), -size.y * (1 - handle.y), 0)
	};
	float3 pos = mul(final_rot_mat, identity_bounds[index]);

	//position.w is the baked stretch, applied forwards and backwards along the alignment
	float3 stretch_direction = billboarding ? instance.alignment : normalize(camera_relative_alignment);
	float vertex_offset_magnitude = 0.5 * instance.position.w;
	float front_back_sign = sign(dot(pos, stretch_direction));
	pos += (instance.position.w > 0.00001) ? stretch_direction * front_back_sign * vertex_offset_magnitude : float3(0.0, 0.0, 0.0);

	float4 view_position = mul(view, float4(instance.position.xyz, 1.0));
	float3 view_offset = billboarding ? mul(view_rotation, pos) : pos;

	tfx_billboard_output output;
	output.position = mul(projection, float4(view_position.xyz + view_offset, view_position.w));

	uint image_index = instance.indexes & 0x00001FFF;
	float4 uv = image_data[image_index].uv;
	float2 uvs[4] = { uv.xy, uv.zy, uv.xw, uv.zw };

	int ramp_y = int(props.color_ramp_indexes & 0xFFu);
	int ramp_array = int((props.color_ramp_indexes >> 8) & 0xFFu);
	output.texture_indexes = int3(ramp_y, ramp_array, 0);
	uint packed_texture_index = image_data[image_index].texture_array_index;
	output.tex_coord = float3(uvs[index], float(packed_texture_index & 0xFFFFu));
	output.intensity_curved_alpha_map = float4(instance.intensity_gradient_map.x * intensity_max_value, instance.curved_alpha_life.x, instance.curved_alpha_life.y, instance.intensity_gradient_map.y * intensity_max_value);
	return output;
}

//--- Ribbon_playback

//tfx_ribbon_data_push_t. The bindless indexes in it are unused here
cbuffer tfx_ribbon_data_constants : register(b4) {
	uint ribbon_animation_instances_total;
	uint ribbons_total;
	uint ribbon_offsets_index;
	uint ribbon_animation_instances_index;
	uint ribbon_data_index;
	uint ribbon_output_index;
	uint ribbon_properties_index;
};

//tfx_ribbon_instance_data_t, 64 bytes
struct tfx_ribbon_data {
	float4 position;				//xyz position, w = normalised age
	float width;
	uint start_index;				//Into the shared segments
	uint flags;
	uint captured_index;			//The same ribbon in the previous frame
	uint2 quaternion;
	float scale;
	uint texture_indexes;
	uint intensity_gradient_map;
	uint curved_alpha;
	uint property_index_lerp_offset;	//Unorm lerp offset in 0xFFFF and the property index in 0xFFFF0000
	uint phase_seed;
};

//tfx_animation_ribbon_properties_t, 68 bytes
struct tfx_ribbon_emitter_properties {
	uint segment_count;
	uint sample_count;
	uint tessellation;
	uint shader_type;
	uint color_ramp_index;
	uint segment_data_offset;
	float animation_frames;
	uint start_frame_index;
	uint flags;
	uint graph_lookup_offset;
	uint bucket_index;
	uint morph_segment_start_index;
	float noise_frequency;
	float noise_speed;
	float noise_phase_range;
	float noise_lock_rate;
	uint noise_packed;
};

StructuredBuffer<tfx_ribbon_data> baked_ribbons : register(t14);		//tfx_GetAnimationRibbonDataBufferPointer
ByteAddressBuffer animation_ribbon_offsets : register(t15);			//tfx_GetRibbonOffsetsBufferPointer, a running ribbon total per animation instance
StructuredBuffer<tfx_ribbon_emitter_properties> baked_ribbon_properties : register(t16);
RWStructuredBuffer<tfx_ribbon> playback_ribbons_out : register(u3);

[numthreads(128, 1, 1)]
void ribbon_playback(uint3 dispatch_thread_id : SV_DispatchThreadID) {
	uint index = dispatch_thread_id.x;
	if (index >= ribbons_total) return;

	uint instance_index = 0;
	for (uint i = 0; i < ribbon_animation_instances_total; ++i) {
		instance_index += uint(animation_ribbon_offsets.Load(instance_index * 4) <= index);
	}
	uint ribbon_index_offset = instance_index == 0 ? 0 : animation_ribbon_offsets.Load((instance_index - 1) * 4);
	tfx_animation_instance animation = tfx_load_animation_instance(instance_index);

	float frame_time = (animation.current_time / animation.animation_length_in_time) * float(animation.frame_count);
	float lerp_amount = frame_time - floor(frame_time);

	uint ribbon_record_index = (index - ribbon_index_offset) + animation.offset_into_ribbon_data;
	tfx_ribbon_data current = baked_ribbons[ribbon_record_index];
	uint captured_index = current.captured_index;
	//A ribbon in its first frame interpolates with itself
	if (captured_index == tfx_invalid) captured_index = ribbon_record_index;
	tfx_ribbon_data captured = baked_ribbons[captured_index];

	//The spawn record holds the lerp_offset and it's the captured ribbon of the next frame, the first one that interpolates from it
	float lerp_offset = float(captured.property_index_lerp_offset & 0xFFFF) * (1.0 / 65535.0);
	float spawn_scale = ((lerp_offset >= (1 - lerp_amount)) && (captured_index != ribbon_record_index)) ? 1 : 0;
	//The captured state is the ribbon at its spawn time, so interpolate from the spawn point rather than the start of the frame
	lerp_amount = saturate((lerp_amount - (1 - lerp_offset)) / max(lerp_offset, 1.0 / 65535.0));

	float3 lerped_position = lerp(captured.position.xyz, current.position.xyz, lerp_amount) * animation.scale;
	float lerped_age = lerp(captured.position.w, current.position.w, lerp_amount);
	//No width until it's spawned
	float lerped_width = lerp(captured.width, current.width, lerp_amount) * spawn_scale;
	float4 current_quaternion = float4(tfx_unpack_snorm16x2(current.quaternion.x), tfx_unpack_snorm16x2(current.quaternion.y));
	float4 captured_quaternion = float4(tfx_unpack_snorm16x2(captured.quaternion.x), tfx_unpack_snorm16x2(captured.quaternion.y));
	float4 lerped_quaternion = lerp(captured_quaternion, current_quaternion, lerp_amount);
	float2 lerped_intensity_gradient_map = lerp(unpack_16bit_2_scale(captured.intensity_gradient_map, intensity_max_value), unpack_16bit_2_scale(current.intensity_gradient_map, intensity_max_value), lerp_amount);
	float4 lerped_curved_alpha = lerp(unpack_unorm4x8(captured.curved_alpha), unpack_unorm4x8(current.curved_alpha), lerp_amount);

	lerped_position += animation.position;
	uint property_index = (current.property_index_lerp_offset & 0xFFFF0000) >> 16;

	tfx_ribbon ribbon;
	ribbon.position = lerped_position;
	ribbon.age_width = uint(clamp(lerped_age, 0.0, 1.0) * 65535.0 + 0.5) | (f32tof16(lerped_width) << 16);
	ribbon.quaternion = uint2(pack_snorm2x16(lerped_quaternion.xy), pack_snorm2x16(lerped_quaternion.zw));
	ribbon.start_index_flags = (current.start_index & 0x00FFFFFFu) | (current.flags << 24);
	//The baked emitter buffer has one row per ribbon property, so the property index is the emitter index
	ribbon.emitter_index_phase_seed = (property_index & 0xFFFFu) | (current.phase_seed << 16);
	ribbon.texture_indexes = (current.texture_indexes & 0xFFFFE000) | baked_ribbon_properties[property_index].start_frame_index;
	ribbon.intensity_gradient_map = (uint(lerped_intensity_gradient_map.y * intensity_max_divisor) << 16) + uint(lerped_intensity_gradient_map.x * intensity_max_divisor);
	ribbon.curved_alpha = pack_unorm4x8(lerped_curved_alpha);
	ribbon.scale = f32tof16(lerp(captured.scale, current.scale, lerp_amount));
	playback_ribbons_out[index] = ribbon;
}
