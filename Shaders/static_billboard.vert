#version 450
#extension GL_EXT_nonuniform_qualifier : require
//This is an all in one billboard vertex shader which can handle both free align and align to camera.
//If you have an effect with multiple emitters using both align methods then you can use this shader to
//draw the effect with one draw call. Otherwise if your effect only uses one or the other you can optimise by
//using the shaders that only handle whatever you need for the effect.

//Quad indexes
const int indexes[6] = int[6]( 0, 1, 2, 2, 1, 3 );

//We add an epsilon to avoid nans with the cross product
const vec3 up = vec3( 0, 1, 0.00001 );
const vec3 front = vec3( 0, 0, 1 );
const vec3 left = vec3( 1, 0, 0 );
const float size_max_value = 64.0 / 32767.0;
const float handle_max_value = 128.0 / 32767.0;
const float intensity_max_value = 128.0 / 32767.0;

struct ImageData {
	vec4 uv;
	vec2 padding;
	vec2 image_size;
	uint texture_array_index;
	float animation_frames;
	float max_radius;
};

struct Properties {
	vec2 image_handle;
	uint color_ramp_indexes;			//[Row of color ramp bitmap, texture array]
	uint flags;							//Flags like billboard alignment type
	uint start_frame_index;
	float animation_frames;
	uint padding1;
	uint padding2;
};

layout(binding = 7) uniform UboView
{
    mat4 view;
    mat4 proj;
    vec2 screen_size;
	uint millisecs;
    float timer_lerp;
    float update_time;
} ub[];

layout(std430, binding = 5) readonly buffer InImageData {
	ImageData data[];
}in_image_data[];

layout(std430, binding = 5) readonly buffer InProperties {
	Properties data[];
} in_properties[];

layout(push_constant) uniform quad_index
{
    vec4 offset;
	uint color_ramp_texture_index;
	uint image_data_index;
	uint properties_index;
	uint sampler_index;
	uint prev_billboards_index;
	uint index_offset;
	uint uniform_index;
	float ndc_offset_x;
	float ndc_offset_y;
} pc;

//Vertex
//layout(location = 0) in vec2 vertex_position;

//Instance
layout(location = 0) in vec4 position;
layout(location = 1) in vec4 in_quaternion;
layout(location = 2) in vec2 size_in;
layout(location = 3) in vec3 alignment;
layout(location = 4) in vec2 intensity_gradient_map;
layout(location = 5) in vec3 curved_alpha_life;
layout(location = 6) in uint texture_indexes;
layout(location = 7) in uint captured_index;

layout(location = 0) out vec3 out_tex_coord;
layout(location = 1) out flat ivec3 out_texture_indexes;
layout(location = 2) out vec4 out_intensity_curved_alpha_map;
//Every shape has its own image so the descriptor index has to reach the fragment shader per particle
layout(location = 3) out flat uint out_image_index;

mat3 QuaternionToRotationMatrix(vec4 q) {
    float xx = q.x * q.x; float xy = q.x * q.y; float xz = q.x * q.z; float xw = q.x * q.w;
    float yy = q.y * q.y; float yz = q.y * q.z; float yw = q.y * q.w;
    float zz = q.z * q.z; float zw = q.z * q.w;
    mat3 m;
    m[0][0] = 1.0 - 2.0 * (yy + zz); m[0][1] = 2.0 * (xy + zw);       m[0][2] = 2.0 * (xz - yw);
    m[1][0] = 2.0 * (xy - zw);       m[1][1] = 1.0 - 2.0 * (xx + zz); m[1][2] = 2.0 * (yz + xw);
    m[2][0] = 2.0 * (xz + yw);       m[2][1] = 2.0 * (yz - xw);       m[2][2] = 1.0 - 2.0 * (xx + yy);
    return m;
}

void main() {
    uint emitter_index = (texture_indexes >> 16) & 0xFFFFu;
    Properties props = in_properties[pc.properties_index].data[emitter_index];

    vec2 size = size_in * size_max_value;
    vec2 handle = props.image_handle;
    vec4 quaternion = normalize(in_quaternion);

    //Info about how to align the billboard is stored in the lowest 2 bits of props.flags.

    //Billboarding determines whether the quad will align to the camera or not. 0 means that it will
    //align to the camera. This value is determined by the first bit: 01
    uint align_flags = props.flags & 0x3u;
    bool billboarding = (align_flags & 0x1u) > 0u;

    //align_type is set to 1 when we want the quad to align to the vector stored in alignment.xyz with
    //no billboarding. Billboarding will always be set to 1 in this case, so in other words both bits
    //will be set: 11.
    bool align_type = align_flags == 0x3u;

    //vector_align is set to 1 when we want the billboard to align to the camera and the vector
    //stored in alignment.xyz. billboarding and align_type will always be 0 if this is the case. the second
    //bit is the only bit set if this is the case: 10
    bool vector_align = align_flags == 0x2u;

    int index = indexes[gl_VertexIndex];

    mat3 final_rot_mat;
    mat3 base_spin_mat = QuaternionToRotationMatrix(quaternion);

    //align_type = true
    //------------------------
	vec3 normal = normalize(alignment);
	vec3 tangent = normalize(cross(up, normal));
	if (length(tangent) < 0.001) {
		tangent = normalize(cross(vec3(1.0, 0.0, 0.0), normal));
	}
	vec3 bitangent = normalize(cross(normal, tangent));
	mat3 align_mat = mat3(tangent, bitangent, normal);
    //------------------------

	//align to camera and vector
    //------------------------
	vec3 camera_relative_alignment = alignment * inverse(mat3(ub[pc.uniform_index].view));
	vec2 align_xy = normalize(camera_relative_alignment.xy + vec2(0.00001));
	vec2 view_up_xy = vec2(0, 1);

	float c_vec = dot(align_xy, view_up_xy);
	float s_vec = align_xy.x * view_up_xy.y - align_xy.y * view_up_xy.x;

	float final_c_roll = c_vec;
	float final_s_roll = -s_vec;

	mat3 vector_align_roll_mat = mat3(final_c_roll,  final_s_roll, 0.0,
									  -final_s_roll, final_c_roll, 0.0,
									   0.0,          0.0,           1.0);

	final_rot_mat = align_type ? align_mat * base_spin_mat : (vector_align ? vector_align_roll_mat * base_spin_mat : base_spin_mat);
	//------------------------

    const vec3 identity_bounds[4] = vec3[4](
        vec3( size.x * (0 - handle.x), -size.y * (0 - handle.y), 0),
        vec3( size.x * (1 - handle.x), -size.y * (0 - handle.y), 0),
        vec3( size.x * (0 - handle.x), -size.y * (1 - handle.y), 0),
        vec3( size.x * (1 - handle.x), -size.y * (1 - handle.y), 0)
    );

    vec3 pos = final_rot_mat * identity_bounds[index];

    mat4 model = mat4(1.0);
    model[3].xyz = position.xyz + pc.offset.xyz;

	mat4 modelView = ub[pc.uniform_index].view * model;

    //Stretch effect: position.w contains the pre-baked stretch magnitude (stretch_factor * travel_distance)
    //Offset vertices forward/backward along the alignment direction
    vec3 stretch_direction = billboarding ? alignment : normalize(camera_relative_alignment);
    float vertex_offset_magnitude = 0.5 * position.w;
    float front_back_sign = sign(dot(pos, stretch_direction));
    pos += (position.w > 0.00001) ? stretch_direction * front_back_sign * vertex_offset_magnitude : vec3(0.0);

    //Billboarding. If billboarding = 0 then billboarding is active and the quad will always face the camera,
    //otherwise the modelView matrix is used as it is.
    modelView[0].xyz = !billboarding ? vec3(1, 0, 0) : modelView[0].xyz;
    modelView[1].xyz = !billboarding ? vec3(0, 1, 0) : modelView[1].xyz;
    modelView[2].xyz = !billboarding ? vec3(0, 0, 1) : modelView[2].xyz;

    vec2 uvs[4];
	uint image_index = texture_indexes & 0x00001FFF;
	vec4 uv = in_image_data[pc.image_data_index].data[image_index].uv;
    uvs[0] = uv.xy;
    uvs[1] = uv.zy;
    uvs[2] = uv.xw;
    uvs[3] = uv.zw;

    vec4 p = modelView * vec4(pos, 1.0);
	gl_Position = ub[pc.uniform_index].proj * p;
	gl_Position.x += pc.ndc_offset_x * gl_Position.w;
	gl_Position.y += pc.ndc_offset_y * gl_Position.w;

    //----------------
	int ramp_y = int(props.color_ramp_indexes & 0xFFu);
	int ramp_array = int((props.color_ramp_indexes >> 8) & 0xFFu);
	out_texture_indexes = ivec3(ramp_y, ramp_array, 0);
	uint packed_texture_index = in_image_data[pc.image_data_index].data[image_index].texture_array_index;
	//Packed as [bindless image index << 16 | animation frame]
	out_image_index = packed_texture_index >> 16;
	out_tex_coord = vec3(uvs[index], float(packed_texture_index & 0xFFFFu));
	out_intensity_curved_alpha_map = vec4(intensity_gradient_map.x * intensity_max_value, curved_alpha_life.x, curved_alpha_life.y, intensity_gradient_map.y * intensity_max_value);
}
