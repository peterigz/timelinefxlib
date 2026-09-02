#version 450
#extension GL_EXT_nonuniform_qualifier : require
//This is an all in one billboard vertex shader which can handle both free align and align to camera.
//If you have an effect with multiple emitters using both align methods then you can use this shader to 
//draw the effect with one draw call. Otherwise if your effect only uses one or the other you can optimise by
//using the shaders that only handle whatever you need for the effect.

//#define NO_LERP

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
	uint flags;						//Flags like billboard alignment type
	uint start_frame_index;
	float animation_frames;
	uint padding1;
	uint padding2;
};

struct BillboardInstance {					//48 bytes, mirrors tfx_instance_t
	vec4 position;							//The position of the sprite, x, y - world, z, w = captured for interpolating
	uvec2 quaternion;			            //Rotations of the sprite packed into a 16-bit snorm quaternion: .x = X|Y, .y = Z|W
	uint size;								//Size of the sprite in pixels packed into a u32 (2 16bit floats)
	uint alignment;	    					//normalised alignment vector 2 floats packed into 16bits or 3 8bit floats for 3d
	uint intensity_gradient_map;			//Multiplier for the color and the gradient map value packed into 16bit scaled floats
	uint curved_alpha_life;					//Sharpness and dissolve amount value for fading the image plus the life of the partile packed into 3 unorms
	uint indexes;							//[gpu property index << 16, capture flag << 15, image data index max 8191 images]
	uint captured_index;					//Index to the sprite in the buffer from the previous frame for interpolation
};

layout(push_constant) uniform quad_index
{
    uint particle_texture_index;
    uint color_ramp_texture_index;
    uint image_data_index;
    uint properties_index;
	uint sampler_index;
    uint prev_billboards_index;
    uint index_offset;
	uint uniform_index;
} pc;

layout(binding = 7) uniform UboView
{
    mat4 view;
    mat4 proj;
    vec2 screen_size;
	uint millisecs;
    float timer_lerp;
    float update_time;
} ub[];

layout (std430, set = 0, binding = 5) readonly buffer InBillboardInstances {
	BillboardInstance data[];
} in_prev_billboards[];

layout (std430, set = 0, binding = 5) readonly buffer InImageData {
	ImageData data[];
} in_image_data[];

layout (std430, set = 0, binding = 5) readonly buffer InProperties {
	Properties data[];
} in_properties[];

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

    uint prev_index = (captured_index & 0x0FFFFFFF) + pc.index_offset;  //Add on an offset if the are multiple draw instructions

    //We won't interpolate the position of the particle for a few possible reasons:
    //  It's the first frame of the particle and there's no previous frame to interpolate with
    //  We're looping back to the start of a path or line
    //  It maybe cancelled becuase too much time passed since the last frame
	float interpolate_is_active = float((texture_indexes & 0x00008000) >> 15);

	uint prev_size_packed = in_prev_billboards[pc.prev_billboards_index].data[prev_index].size;

	//vec3 lerped_position = interpolate_is_active == 1 ? mix(in_prev_billboards[pc.prev_billboards_index].data[prev_index].position.xyz, position.xyz, ub[pc.uniform_index].timer_lerp) : position.xyz;

    #ifdef NO_LERP
	vec2 lerped_size = size;
	vec3 lerped_position = position.xyz;
	vec4 lerped_quaternion = quaternion;
    vec3 motion = vec3(0);
    #else
    //Otherwise just hide the first frame of the particle which is a little more efficient:
	vec2 lerped_size = vec2(float(prev_size_packed & 0xFFFF) * size_max_value, float((prev_size_packed & 0xFFFF0000) >> 16) * size_max_value);
	lerped_size = mix(lerped_size, size, ub[pc.uniform_index].timer_lerp.x) * interpolate_is_active;
	vec3 lerped_position = mix(in_prev_billboards[pc.prev_billboards_index].data[prev_index].position.xyz, position.xyz, ub[pc.uniform_index].timer_lerp.x);
    //(peviously lerped_rotation)
	vec4 lerped_quaternion = mix(normalize(vec4(unpackSnorm2x16(in_prev_billboards[pc.prev_billboards_index].data[prev_index].quaternion.x), unpackSnorm2x16(in_prev_billboards[pc.prev_billboards_index].data[prev_index].quaternion.y))), quaternion, ub[pc.uniform_index].timer_lerp.x);
    vec3 motion = position.xyz - in_prev_billboards[pc.prev_billboards_index].data[prev_index].position.xyz;
    #endif
    motion.z += 0.000001;
	bool has_alignment = dot(alignment, alignment) > 0;
    float travel_distance = length(motion); // Actual distance travelled this frame. Stretch always scales with the real motion now.
    float stretch_factor = position.w * interpolate_is_active;

    motion = normalize(motion); // Normalize for direction
    vec3 final_alignment = has_alignment ? alignment : motion; // Vector the billboard is oriented to: motion for motion-align, otherwise the supplied alignment vector.
    //----

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

    //Stretch always follows the direction of travel (motion), never the alignment vector.
    vec3 camera_relative_aligment = motion * inverse(mat3(ub[pc.uniform_index].view));

    int index = indexes[gl_VertexIndex];

    mat3 final_rot_mat;
    mat3 base_spin_mat = QuaternionToRotationMatrix(lerped_quaternion); // Particle's intrinsic rotation

    //align_type = true
    //------------------------
	vec3 normal = normalize(final_alignment); // Z axis = alignment vector
	vec3 tangent = normalize(cross(up, normal)); // X axis
	if (length(tangent) < 0.001) { // Handle case where normal is parallel to up
		tangent = normalize(cross(vec3(1.0, 0.0, 0.0), normal));
	}
	vec3 bitangent = normalize(cross(normal, tangent)); // Y axis
	mat3 align_mat = mat3(tangent, bitangent, normal); // World axes for the aligned frame
    //------------------------

	//align to camera and vector
    //------------------------
	// Project alignment vector onto view plane (XY plane in view space)
	vec3 camera_relative_alignment = final_alignment * inverse(mat3(ub[pc.uniform_index].view));
	vec2 align_xy = normalize(camera_relative_alignment.xy + vec2(0.00001)); // Add epsilon
	vec2 view_up_xy = vec2(0, 1);

	// Calculate cos/sin of angle between projected vector and view Up
	float c_vec = dot(align_xy, view_up_xy);
	float s_vec = align_xy.x * view_up_xy.y - align_xy.y * view_up_xy.x; // Determinant for sine

	// --- Simplification: Use ONLY the vector alignment for roll ---
	// This creates a roll relative to the view plane based purely on the alignment vector.
	// It ignores the intrinsic roll potentially contained within lerped_quaternion.
	// If combined roll is needed, extract roll from lerped_quaternion (complex) & combine sin/cos.
	float final_c_roll = c_vec;
	float final_s_roll = -s_vec;
	// --- End Simplification ---

	// Construct Z-rotation matrix using this calculated roll
	mat3 vector_align_roll_mat = mat3(final_c_roll,  final_s_roll, 0.0,
									  -final_s_roll, final_c_roll, 0.0,
									   0.0,          0.0,           1.0);

	final_rot_mat = align_type ? align_mat * base_spin_mat : (vector_align ? vector_align_roll_mat * base_spin_mat : base_spin_mat);
	// --- Final Rotation ---

    const vec3 identity_bounds[4] = vec3[4](
        vec3( lerped_size.x * (0 - handle.x), -lerped_size.y * (0 - handle.y), 0),
        vec3( lerped_size.x * (1 - handle.x), -lerped_size.y * (0 - handle.y), 0),
        vec3( lerped_size.x * (0 - handle.x), -lerped_size.y * (1 - handle.y), 0),
        vec3( lerped_size.x * (1 - handle.x), -lerped_size.y * (1 - handle.y), 0)
    );
    vec3 pos = final_rot_mat * identity_bounds[index];

    mat4 model = mat4(1.0);
    model[3].xyz = lerped_position;

    mat4 modelView = ub[pc.uniform_index].view * model;

    //Calculate the amount to stretch the particles. A stretch value (set in the editor) of 1 means that the particle is stretched
    //by the amount travelled.

    //Get the direction to stretch based on the billboarding option
	vec3 stretch_direction = billboarding ? motion : normalize(camera_relative_aligment); // Stretch along the direction of travel
    //Half the amount in order to apply half of the stretch in the forward direction and the other half behind
	float vertex_offset_magnitude = 0.5 * travel_distance * stretch_factor;

	//The stretch elongates the sprite by travel_distance: "front" corners are pushed forward along stretch_direction and
	//"back" corners backward. This used to be sign(dot(pos, stretch_direction)) * magnitude, but sign() is discontinuous -
	//as stretch_direction sweeps past a quad corner (dot -> 0) that corner flips from +magnitude to -magnitude and the
	//silhouette snaps. That was hidden while the quad was rotated to line up with stretch_direction (a clean 2-front/2-back
	//split), but fixed-align particles rotate the quad to their align vector (path gradient / emission dir) while stretch
	//follows motion - a different direction - so the split is unstable. Instead scale each vertex along stretch_direction by
	//its own projection, normalised by the quad's projected half-extent. This is continuous (no snapping), keeps the total
	//elongation equal to travel_distance, and reduces exactly to the old rigid front/back push when the quad *is* aligned
	//with the stretch direction.
	float projected_half_extent = 0.0;
	for (int corner = 0; corner < 4; ++corner) {
		projected_half_extent = max(projected_half_extent, abs(dot(final_rot_mat * identity_bounds[corner], stretch_direction)));
	}
	float projected_position = dot(pos, stretch_direction);
	vec3 potential_offset = stretch_direction * (projected_position / max(projected_half_extent, 0.00001)) * vertex_offset_magnitude;

    //Apply the stretch if applicable
	bool apply_stretch = (travel_distance > 0.00001 && stretch_factor > 0.0);
	pos += apply_stretch ? potential_offset : vec3(0.0);

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

    //----------------
	int life = int(curved_alpha_life.z * 255);
	int ramp_y = int(props.color_ramp_indexes & 0xFFu);
	int ramp_array = int((props.color_ramp_indexes >> 8) & 0xFFu);
	out_texture_indexes = ivec3(ramp_y, ramp_array, life);
	out_tex_coord = vec3(uvs[index], in_image_data[pc.image_data_index].data[image_index].texture_array_index);
	out_intensity_curved_alpha_map = vec4(intensity_gradient_map.x * intensity_max_value, curved_alpha_life.x, curved_alpha_life.y, intensity_gradient_map.y * intensity_max_value);
}
