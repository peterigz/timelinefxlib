//Spawn generation for GPU particles, the mirror of tfx__generate_gpu_spawn: one particle from a tick's spawn block, its location's run and
//multipliers. Everything here is a pure function of its arguments, so a test entry compiled with slangc -target cpp can be checked
//against the library. Change the two together.

//tfx_gpu_spawn_block_flag_bits
static const uint tfx_spawn_relative = 1 << 0;
static const uint tfx_spawn_random_color = 1 << 1;
static const uint tfx_spawn_random_start_frame = 1 << 2;
static const uint tfx_spawn_base_uniform_size = 1 << 3;
static const uint tfx_spawn_spin_pitch_and_yaw = 1 << 4;
static const uint tfx_spawn_specify_pitch = 1 << 5;
static const uint tfx_spawn_specify_yaw = 1 << 6;
static const uint tfx_spawn_specify_roll = 1 << 7;
static const uint tfx_spawn_relative_angle = 1 << 8;
static const uint tfx_spawn_emission_cone = 1 << 9;
static const uint tfx_spawn_noise_offsets = 1 << 10;
static const uint tfx_spawn_noise_resolution = 1 << 11;
static const uint tfx_spawn_user = 1 << 12;
static const uint tfx_spawn_fill_area = 1 << 13;
static const uint tfx_spawn_open_ends = 1 << 14;
static const uint tfx_spawn_uniform_distribution = 1 << 15;
static const uint tfx_spawn_shape = 1 << 16;
static const uint tfx_spawn_grid = 1 << 17;
static const uint tfx_spawn_grid_random = 1 << 18;
static const uint tfx_spawn_grid_backwards = 1 << 19;
static const uint tfx_spawn_stepped_emission = 1 << 20;
static const uint tfx_spawn_stepped_roll = 1 << 21;

//tfx_angle_step_distribution and tfx_angle_step_flag_bits
static const uint tfx_step_ping_pong = 2;
static const uint tfx_step_shuffled = 3;
static const uint tfx_step_restart_each_spawn = 1 << 0;
static const uint tfx_step_ring = 1 << 1;
static const uint tfx_step_cap = 1 << 2;

//tfx_emission_type and tfx_emission_direction values the shapes use
static const uint tfx_emission_area = 1;
static const uint tfx_emission_line = 2;
static const uint tfx_emission_ellipse = 3;
static const uint tfx_emission_cylinder = 4;
static const uint tfx_emission_disc = 9;
static const uint tfx_direction_inwards = 0;
static const uint tfx_direction_outwards = 1;
static const uint tfx_direction_bothways = 2;
static const uint tfx_direction_surface = 4;

//tfx_gpu_spawn_run_flag_bits
static const uint tfx_spawn_run_rotated = 1 << 0;
static const uint tfx_spawn_run_scale_life = 1 << 1;

//tfx_spawn_hash_axis
static const uint tfx_spawn_axis_life = 0x51a70001u;
static const uint tfx_spawn_axis_random_color = 0x51a70002u;
static const uint tfx_spawn_axis_start_frame = 0x51a70003u;
static const uint tfx_spawn_axis_width = 0x51a70004u;
static const uint tfx_spawn_axis_height = 0x51a70005u;
static const uint tfx_spawn_axis_weight = 0x51a70006u;
static const uint tfx_spawn_axis_roll_spin = 0x51a70007u;
static const uint tfx_spawn_axis_pitch_spin = 0x51a70008u;
static const uint tfx_spawn_axis_yaw_spin = 0x51a70009u;
static const uint tfx_spawn_axis_pitch_offset = 0x51a7000au;
static const uint tfx_spawn_axis_yaw_offset = 0x51a7000bu;
static const uint tfx_spawn_axis_roll_offset = 0x51a7000cu;
static const uint tfx_spawn_axis_splatter_z = 0x51a7000du;
static const uint tfx_spawn_axis_splatter_phi = 0x51a7000eu;
static const uint tfx_spawn_axis_splatter_radius = 0x51a7000fu;
static const uint tfx_spawn_axis_cone_z = 0x51a70010u;
static const uint tfx_spawn_axis_cone_phi = 0x51a70011u;
static const uint tfx_spawn_axis_noise_x = 0x51a70012u;
static const uint tfx_spawn_axis_noise_y = 0x51a70013u;
static const uint tfx_spawn_axis_noise_z = 0x51a70014u;
static const uint tfx_spawn_axis_shape_side = 0x51a70015u;
static const uint tfx_spawn_axis_shape_u = 0x51a70016u;
static const uint tfx_spawn_axis_shape_v = 0x51a70017u;
static const uint tfx_spawn_axis_shape_w = 0x51a70018u;
static const uint tfx_spawn_axis_sequence_theta = 0x51a70019u;
static const uint tfx_spawn_axis_sequence_phi = 0x51a7001au;
static const uint tfx_spawn_axis_sequence_radius = 0x51a7001bu;
static const uint tfx_spawn_axis_emission_jitter = 0x51a7001cu;
static const uint tfx_spawn_axis_roll_jitter = 0x51a7001du;
//Must match tfxSPAWN_IMPULSE_HASH_AXIS, which tfx_apply_velocity relaunches looping singles on
static const uint tfx_spawn_impulse_hash_axis = 0x2545f491u;
//tfxPI2
static const float tfx_spawn_pi2 = 6.283185307;
//tfxGOLDEN_ANGLE
static const float tfx_spawn_golden_angle = 2.39996323;

//tfx_gpu_spawn_block_t is five 64 byte rows of the tick array, read through a view per row
struct tfx_gpu_spawn_block_a {
	uint uid_base;
	uint flags;
	uint run_row;						//From the block's first row
	uint run_count;
	float3 emission_direction;
	float emission_cone_cosine;
	float3 handle;
	float splatter;
	float3 noise_step_offset;
	float noise_resolution;
};

struct tfx_gpu_spawn_block_b {
	float3 angle_offsets;				//Pitch, yaw, roll
	float animation_frames;
	float2 size;
	float2 size_variation;
	float frame_length;
	float inverse_frame_length;
	float overall_scale;
	float velocity;
	float life;
	float life_variation;
	float weight;
	float weight_variation;
};

struct tfx_gpu_spawn_block_c {
	float roll_spin;
	float pitch_spin;
	float yaw_spin;
	float roll_spin_variation;
	float pitch_spin_variation;
	float yaw_spin_variation;
	float spawn_impulse;
	float spawn_impulse_variation;
	float noise_offset_variation;
	float noise_shared_offset;
	float start_frame;
	float arc_size;						//Disc and cylinder sweep. Ellipse: the cosine of the bottom of its polar band
	float arc_offset;					//Ellipse: the cosine of the top of its polar band
	uint spawn_ordinal;
	uint padding[2];
};

//Only read when the block has tfx_spawn_shape
struct tfx_gpu_spawn_block_d {
	float3 emitter_size;
	uint shape;
	float4 shape_rotation;
	float4 emission_rotation;
	float3 local_handle;
	uint emission_direction_mode;
};

//Only read when the block has tfx_spawn_grid or a stepped flag. Scalars, so the layout matches the C struct under std430
struct tfx_gpu_spawn_block_e {
	uint grid_start;
	uint grid_point_x;
	uint grid_point_y;
	uint grid_point_z;
	float3 emission_fan_axis;
	float emission_range;
	float emission_step_jitter;
	float roll_step_jitter;
	uint step_divisions;				//Emission in the low 16 bits, roll in the high 16
	uint step_settings;					//Emission distribution | roll distribution << 2 | emission step flags << 4 | roll step flags << 8
	uint step_seed;
	uint step_padding[3];
};

//Mirrors tfx_gpu_spawn_run_t
struct tfx_gpu_spawn_run {
	uint local_id;
	uint first_spawn;
	uint count;
	uint multipliers_row;				//From the block's first row, 0xFFFFFFFF when the location has none
	float3 position;
	uint flags;
	float3 captured_position;
	float location_roll;
	uint2 packed_rotation;
	uint user_data;
	uint spawn_ordinal;
};

//Mirrors tfx_gpu_spawn_multipliers_t
struct tfx_gpu_spawn_multipliers {
	float particle_width;
	float particle_height;
	float velocity;
	float inverse_life;
	float spin;
	float intensity;
	float splatter;
	float weight;
	float noise_offset;
	float emitter_width;
	float emitter_height;
	float emitter_depth;
	uint4 padding;
};

//One generated particle, as the CPU bank holds it at spawn. image_frame is the frame at spawn, not the frame at age 0
struct tfx_generated_spawn {
	float3 position;
	float age;
	float3 velocity;
	uint uid;
	float inv_max_age;
	uint velocity_normal;
	float base_velocity;
	float base_weight;
	float2 base_size;
	float intensity_factor;
	float random_color;
	uint rotation_offsets;
	float base_roll_spin;
	float base_pitch_spin;
	float base_yaw_spin;
	float image_frame;
	uint spawn_location;
	float3 noise_offset;
	float noise_resolution;
};

//tfx__seedgen_u32, a murmur3 finalizer
uint tfx_seedgen(uint hash) {
	hash ^= hash >> 15;
	hash *= 0x85ebca6bu;
	hash ^= hash >> 13;
	hash *= 0xc2b2ae35u;
	hash ^= hash >> 16;
	return hash;
}

//Port of tfx__white_unit: the signed reinterpret is what centres it on 0
float tfx_white_unit(uint seed, uint axis) {
	return (float(int(tfx_seedgen(seed ^ axis))) / 4294967295.0) * 2.0;
}

//tfx__hash_unit, exact on both sides because 24 bits convert to a float without rounding
float tfx_hash_unit(uint seed, uint axis) {
	return float(tfx_seedgen(seed ^ axis) >> 8) * (1.0 / 16777216.0);
}

//tfx__wide_unpack10bit: x in the top 10 bits, z in the bottom
float3 tfx_unpack10bit(uint packed) {
	const float one_div_511 = 1.0 / 511.0;
	return float3(float((packed >> 20) & 0x3FFu), float((packed >> 10) & 0x3FFu), float(packed & 0x3FFu)) * one_div_511 - 1.0;
}

//tfx__pack10bit_unsigned
uint tfx_pack10bit_unsigned(float3 v) {
	const uint x = uint(mad(v.x, 511.0, 511.0)) & 0x3FFu;
	const uint y = uint(mad(v.y, 511.0, 511.0)) & 0x3FFu;
	const uint z = uint(mad(v.z, 511.0, 511.0)) & 0x3FFu;
	return z | (y << 10) | (x << 20);
}

//tfx__wide_unpack16bit: each lane sign extended, x and y from .x, z and w from .y
float4 tfx_unpack16bit_quaternion(uint2 packed) {
	const int4 values = int4(int(packed.x << 16) >> 16, int(packed.x) >> 16, int(packed.y << 16) >> 16, int(packed.y) >> 16);
	return float4(values) * (1.0 / 32767.0);
}

//The CPU's rotate by quaternion, tfx__wide_transform_quaternion_vec3
float3 tfx_rotate_vector(float3 v, float4 q) {
	return v + 2.0 * cross(cross(v, q.xyz) + q.w * v, q.xyz);
}

//tfx_quaternion_t's operator*, xyzw
float4 tfx_multiply_quaternions(float4 a, float4 b) {
	return float4(
		a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
		a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
		a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
		a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z);
}

//tfx__spin_euler_to_quaternion, angles as pitch, yaw, roll
float4 tfx_spin_euler_to_quaternion(float3 angles) {
	const float sin_pitch = sin(angles.x * 0.5), cos_pitch = cos(angles.x * 0.5);
	const float sin_yaw = sin(angles.y * 0.5), cos_yaw = cos(angles.y * 0.5);
	const float sin_roll = sin(angles.z * 0.5), cos_roll = cos(angles.z * 0.5);
	return float4(
		cos_roll * sin_pitch * cos_yaw + sin_roll * cos_pitch * sin_yaw,
		cos_roll * cos_pitch * sin_yaw - sin_roll * sin_pitch * cos_yaw,
		sin_roll * cos_pitch * cos_yaw - cos_roll * sin_pitch * sin_yaw,
		cos_roll * cos_pitch * cos_yaw + sin_roll * sin_pitch * sin_yaw);
}

//tfx__quaternion_to_spin_euler
float3 tfx_quaternion_to_spin_euler(float4 q) {
	const float sin_pitch = 2.0 * (q.w * q.x - q.y * q.z);
	return float3(
		asin(clamp(sin_pitch, -1.0, 1.0)),
		atan2(2.0 * (q.w * q.y + q.z * q.x), 1.0 - 2.0 * (q.x * q.x + q.y * q.y)),
		atan2(2.0 * (q.w * q.z + q.x * q.y), 1.0 - 2.0 * (q.z * q.z + q.x * q.x)));
}

//tfx__quake_sqrt, the inverse square root with one Newton step
float tfx_quake_inverse_sqrt(float number) {
	const float half_number = number * 0.5;
	float estimate = asfloat(0x5f3759dfu - (asuint(number) >> 1));
	estimate = estimate * (1.5 - (half_number * estimate * estimate));
	return estimate;
}

//tfx__normalize_vec3_fast
float3 tfx_normalize_fast(float3 v) {
	return v * tfx_quake_inverse_sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

//tfx__vector_in_cone: z is the cosine of the angle from cone_direction and phi the turn about it
float3 tfx_vector_in_cone(float3 cone_direction, float z, float phi) {
	const float ring = sqrt(1.0 - z * z);
	const float x = ring * cos(phi);
	const float y = ring * sin(phi);
	float3 random_vector = float3(x, y, z);
	//On either pole there is no rotation axis, so mirror the sample for the south pole
	if (cone_direction.x == 0.0 && cone_direction.y == 0.0) {
		if (cone_direction.z < 0.0) {
			random_vector.z = -z;
		}
		return random_vector;
	}
	const float3 rotation_axis = tfx_normalize_fast(float3(-cone_direction.y, cone_direction.x, 0.0));
	const float rotation_angle = acos(clamp(cone_direction.z, -1.0, 1.0));
	const float cosine = cos(rotation_angle);
	const float sine = sin(rotation_angle);
	const float dot_axis = rotation_axis.x * x + rotation_axis.y * y;
	return random_vector * cosine + cross(rotation_axis, random_vector) * sine + rotation_axis * dot_axis * (1.0 - cosine);
}

//tfx__radical_inverse
float tfx_radical_inverse(uint index, uint base) {
	const float inverse_base = 1.0 / float(base);
	float digit_scale = inverse_base;
	float result = 0.0;
	while (index != 0) {
		result += float(index % base) * digit_scale;
		digit_scale *= inverse_base;
		index /= base;
	}
	return result;
}

//tfx__shifted_sequence
float tfx_shifted_sequence(uint sequence_index, uint base, uint uid_base, uint axis) {
	const float unit = tfx_radical_inverse(sequence_index, base) + tfx_hash_unit(uid_base, axis);
	return unit >= 1.0 ? unit - 1.0 : unit;
}

//tfx__grid_march_length
uint tfx_grid_march_length(uint shape, uint flags, uint3 grid_points) {
	const uint x = grid_points.x, y = grid_points.y, z = grid_points.z;
	if (shape == tfx_emission_area) {
		if ((flags & tfx_spawn_fill_area) != 0) return x * y * z;
		const uint ring = y == 1 ? x : 2 * x + 2 * (y - 2);
		if ((flags & tfx_spawn_open_ends) != 0) return z * ring;
		return z == 1 ? x * y : 2 * x * y + (z - 2) * ring;
	}
	if (shape == tfx_emission_line) return y;
	if (shape == tfx_emission_cylinder) return x * y;
	return x;
}

//tfx__hollow_area_ring_point
uint3 tfx_hollow_area_ring_point(uint ring_index, uint x, uint y) {
	const uint ring = y == 1 ? x : 2 * x + 2 * (y - 2);
	if (y == 1 || ring_index < x) {
		return uint3(ring_index, 0, 0);
	}
	if (ring_index >= ring - x) {
		return uint3(ring_index - (ring - x), y - 1, 0);
	}
	const uint side_index = ring_index - x;
	return uint3((side_index & 1) != 0 ? x - 1 : 0, 1 + side_index / 2, 0);
}

//tfx__grid_march_point
uint3 tfx_grid_march_point(uint shape, uint flags, uint3 grid_points, uint position) {
	const uint x = grid_points.x, y = grid_points.y, z = grid_points.z;
	if (shape == tfx_emission_area) {
		if ((flags & tfx_spawn_fill_area) != 0) {
			return uint3(position % x, (position / x) % y, position / (x * y));
		}
		const uint ring = y == 1 ? x : 2 * x + 2 * (y - 2);
		if ((flags & tfx_spawn_open_ends) != 0) {
			uint3 coordinates = tfx_hollow_area_ring_point(position % ring, x, y);
			coordinates.z = position / ring;
			return coordinates;
		}
		if (z == 1 || position < x * y) {
			return uint3(position % x, position / x, 0);
		}
		const uint side_position = position - x * y;
		if (side_position < (z - 2) * ring) {
			uint3 coordinates = tfx_hollow_area_ring_point(side_position % ring, x, y);
			coordinates.z = 1 + side_position / ring;
			return coordinates;
		}
		const uint end_position = side_position - (z - 2) * ring;
		return uint3(end_position % x, end_position / x, z - 1);
	}
	if (shape == tfx_emission_line) return uint3(0, position, 0);
	if (shape == tfx_emission_cylinder) return uint3(position % x, position / x, 0);
	return uint3(position, 0, 0);
}

//tfx__generate_grid_offset
float3 tfx_generate_grid_offset(tfx_gpu_spawn_block_a block_a, tfx_gpu_spawn_block_c block_c, tfx_gpu_spawn_block_d block_d, tfx_gpu_spawn_block_e block_e, uint uid, uint spawn_index, inout float3 face_normal) {
	const uint shape = block_d.shape;
	const uint3 grid_points = uint3(block_e.grid_point_x, block_e.grid_point_y, block_e.grid_point_z);
	const float3 emitter_size = block_d.emitter_size;
	const float3 half_size = emitter_size * 0.5;
	const bool hollow_area = shape == tfx_emission_area && (block_a.flags & tfx_spawn_fill_area) == 0;
	uint3 coordinates = uint3(0, 0, 0);
	if ((block_a.flags & tfx_spawn_grid_random) != 0) {
		coordinates.x = uint(tfx_hash_unit(uid, tfx_spawn_axis_shape_u) * float(grid_points.x));
		coordinates.y = uint(tfx_hash_unit(uid, tfx_spawn_axis_shape_v) * float(grid_points.y));
		coordinates.z = uint(tfx_hash_unit(uid, tfx_spawn_axis_shape_w) * float(grid_points.z));
		if (hollow_area) {
			const uint sides = (block_a.flags & tfx_spawn_open_ends) != 0 ? 4u : 6u;
			const uint side = uint(tfx_hash_unit(uid, tfx_spawn_axis_shape_side) * float(sides));
			face_normal = float3(0.0, 0.0, 0.0);
			if (side == 0) { coordinates.x = 0; face_normal.x = -1.0; }
			else if (side == 1) { coordinates.x = grid_points.x - 1; face_normal.x = 1.0; }
			else if (side == 2) { coordinates.y = 0; face_normal.y = -1.0; }
			else if (side == 3) { coordinates.y = grid_points.y - 1; face_normal.y = 1.0; }
			else if (side == 4) { coordinates.z = grid_points.z - 1; face_normal.z = 1.0; }
			else { coordinates.z = 0; face_normal.z = -1.0; }
		}
	} else {
		const uint length = tfx_grid_march_length(shape, block_a.flags, grid_points);
		const uint travel = spawn_index % length;
		const uint position = (block_a.flags & tfx_spawn_grid_backwards) != 0 ? (block_e.grid_start + length - 1 - travel) % length : (block_e.grid_start + travel) % length;
		coordinates = tfx_grid_march_point(shape, block_a.flags, grid_points, position);
		if (hollow_area) {
			face_normal.x = coordinates.x == 0 ? -1.0 : (coordinates.x == grid_points.x - 1 ? 1.0 : 0.0);
			face_normal.y = coordinates.y == 0 ? -1.0 : (coordinates.y == grid_points.y - 1 ? 1.0 : 0.0);
			face_normal.z = coordinates.z == 0 ? -1.0 : (coordinates.z == grid_points.z - 1 ? 1.0 : 0.0);
		}
	}
	if (shape == tfx_emission_area) {
		const float3 segment = emitter_size / max(float3(grid_points) - 1.0, float3(1.0, 1.0, 1.0));
		return float3(coordinates) * segment - half_size;
	}
	if (shape == tfx_emission_line) {
		return float3(0.0, float(coordinates.y) * emitter_size.y / max(float(grid_points.y) - 1.0, 1.0), 0.0);
	}
	if (shape == tfx_emission_cylinder) {
		const float theta = float(coordinates.x) * block_c.arc_size / float(grid_points.x) + block_c.arc_offset;
		return float3(cos(theta) * half_size.x, float(coordinates.y) * emitter_size.y / max(float(grid_points.y) - 1.0, 1.0), -sin(theta) * half_size.z);
	}
	const float theta = block_c.arc_offset + float(coordinates.x) * block_c.arc_size / float(grid_points.x);
	return float3(half_size.x * cos(theta), 0.0, half_size.z * -sin(theta));
}

//tfx__permute_index
uint tfx_permute_index(uint index, uint length, uint seed) {
	uint mask = length - 1;
	mask |= mask >> 1;
	mask |= mask >> 2;
	mask |= mask >> 4;
	mask |= mask >> 8;
	mask |= mask >> 16;
	do {
		index ^= seed;
		index *= 0xe170893du;
		index ^= seed >> 16;
		index ^= (index & mask) >> 4;
		index ^= seed >> 8;
		index *= 0x0929eb3fu;
		index ^= seed >> 23;
		index ^= (index & mask) >> 1;
		index *= 1u | seed >> 27;
		index *= 0x6935fa69u;
		index ^= (index & mask) >> 11;
		index *= 0x74dcb303u;
		index ^= (index & mask) >> 2;
		index *= 0x9e501cc3u;
		index ^= (index & mask) >> 2;
		index *= 0xc860a3dfu;
		index &= mask;
		index ^= index >> 5;
	} while (index >= length);
	return (index + seed) % length;
}

//tfx__angle_step_span, start in x and step in y
float2 tfx_angle_step_span(float range, float centre, uint divisions) {
	divisions = max(divisions, 1u);
	const bool full_circle = range >= tfx_spawn_pi2 - 0.001;
	const float step = full_circle ? range / float(divisions) : (divisions > 1 ? range / float(divisions - 1) : 0.0);
	const float start = full_circle || divisions == 1 ? centre : centre - range * 0.5;
	return float2(start, step);
}

//tfx__angle_step_division
uint tfx_angle_step_division(uint distribution, uint divisions, uint ordinal, uint seed) {
	divisions = max(divisions, 1u);
	if (divisions == 1) {
		return 0;
	}
	const uint period = distribution == tfx_step_ping_pong ? (divisions - 1) * 2 : divisions;
	const uint position = ordinal % period;
	if (distribution == tfx_step_ping_pong) {
		return position < divisions ? position : period - position;
	}
	if (distribution == tfx_step_shuffled) {
		return tfx_permute_index(position, divisions, seed ^ ((ordinal / period) * 0x9E3779B9u));
	}
	return position;
}

//tfx__generate_stepped_direction
float3 tfx_generate_stepped_direction(tfx_gpu_spawn_block_e block_e, tfx_gpu_spawn_run run, uint uid, uint spawn_index, float3 direction) {
	const uint step_flags = (block_e.step_settings >> 4) & 0xFu;
	const uint divisions = max(block_e.step_divisions & 0xFFFFu, 1u);
	const uint ordinal = ((step_flags & tfx_step_restart_each_spawn) != 0 ? 0 : run.spawn_ordinal) + (spawn_index - run.first_spawn);
	const uint division = tfx_angle_step_division(block_e.step_settings & 3u, divisions, ordinal, block_e.step_seed);
	const float3 axis = block_e.emission_fan_axis;
	const float half_range = block_e.emission_range * 0.5;
	if ((step_flags & tfx_step_cap) != 0) {
		const float cap_height_step = (divisions > 1 ? 1.0 - cos(half_range) : 0.0) / float(divisions);
		const float cos_tilt = 1.0 - cap_height_step * (float(division) + 0.5);
		const float sin_tilt = sqrt(max(1.0 - cos_tilt * cos_tilt, 0.0));
		const float angle = float(division) * tfx_spawn_golden_angle;
		const float3 radial = axis * cos(angle) + cross(direction, axis) * sin(angle);
		float3 stepped = direction * cos_tilt + radial * sin_tilt;
		const float jitter_cone = block_e.emission_step_jitter * 0.5 * sqrt(tfx_spawn_pi2 * cap_height_step);
		if (jitter_cone != 0.0) {
			const float z = 1.0 - (1.0 - cos(jitter_cone)) * tfx_hash_unit(uid, tfx_spawn_axis_cone_z);
			const float phi = tfx_spawn_pi2 * tfx_hash_unit(uid, tfx_spawn_axis_cone_phi);
			stepped = tfx_vector_in_cone(tfx_normalize_fast(stepped), z, phi);
		}
		return stepped;
	}
	const bool ring = (step_flags & tfx_step_ring) != 0;
	const float2 span = tfx_angle_step_span(ring ? tfx_spawn_pi2 : block_e.emission_range, 0.0, divisions);
	float angle = span.x + span.y * float(division);
	const float half_jitter = span.y * block_e.emission_step_jitter * 0.5;
	if (half_jitter > 0.0) {
		angle += (2.0 * tfx_hash_unit(uid, tfx_spawn_axis_emission_jitter) - 1.0) * half_jitter;
	}
	if (ring) {
		const float3 radial = axis * cos(angle) + cross(direction, axis) * sin(angle);
		return direction * cos(half_range) + radial * sin(half_range);
	}
	return direction * cos(angle) + cross(axis, direction) * sin(angle);
}

//tfx__generate_stepped_roll
float tfx_generate_stepped_roll(tfx_gpu_spawn_block_b block_b, tfx_gpu_spawn_block_e block_e, tfx_gpu_spawn_run run, uint uid, uint spawn_index) {
	const uint step_flags = (block_e.step_settings >> 8) & 0xFu;
	const uint divisions = max(block_e.step_divisions >> 16, 1u);
	const uint ordinal = ((step_flags & tfx_step_restart_each_spawn) != 0 ? 0 : run.spawn_ordinal) + (spawn_index - run.first_spawn);
	const uint division = tfx_angle_step_division((block_e.step_settings >> 2) & 3u, divisions, ordinal, block_e.step_seed);
	const float range = block_b.angle_offsets.z;
	const float centre = range >= tfx_spawn_pi2 - 0.001 ? 0.0 : range * 0.5;
	const float2 span = tfx_angle_step_span(range, centre, divisions);
	float angle = span.x + span.y * float(division);
	const float half_jitter = span.y * block_e.roll_step_jitter * 0.5;
	if (half_jitter > 0.0) {
		angle += (2.0 * tfx_hash_unit(uid, tfx_spawn_axis_roll_jitter) - 1.0) * half_jitter;
	}
	return angle;
}

//tfx__generate_shape_offset
float3 tfx_generate_shape_offset(tfx_gpu_spawn_block_a block_a, tfx_gpu_spawn_block_c block_c, tfx_gpu_spawn_block_d block_d, tfx_gpu_spawn_block_e block_e, uint uid, uint spawn_index, out float3 face_normal) {
	const float3 emitter_size = block_d.emitter_size;
	const float3 half_size = emitter_size * 0.5;
	const bool fill_area = (block_a.flags & tfx_spawn_fill_area) != 0;
	face_normal = float3(0.0, 0.0, 0.0);
	if ((block_a.flags & tfx_spawn_grid) != 0) {
		return tfx_generate_grid_offset(block_a, block_c, block_d, block_e, uid, spawn_index, face_normal);
	}
	float unit_u = tfx_hash_unit(uid, tfx_spawn_axis_shape_u);
	float unit_v = tfx_hash_unit(uid, tfx_spawn_axis_shape_v);
	float unit_w = tfx_hash_unit(uid, tfx_spawn_axis_shape_w);
	if (block_d.shape == tfx_emission_line) {
		return float3(0.0, unit_u * emitter_size.y, 0.0);
	}
	if (block_d.shape == tfx_emission_area) {
		const float3 inside = float3(-half_size.x + unit_u * emitter_size.x, -half_size.y + unit_v * emitter_size.y, -half_size.z + unit_w * emitter_size.z);
		if (fill_area) {
			return inside;
		}
		const uint sides = (block_a.flags & tfx_spawn_open_ends) != 0 ? 4u : 6u;
		const uint side = uint(tfx_hash_unit(uid, tfx_spawn_axis_shape_side) * float(sides));
		if (side == 0) { face_normal.x = -1.0; return float3(-half_size.x, inside.y, inside.z); }
		if (side == 1) { face_normal.x = 1.0; return float3(half_size.x, inside.y, inside.z); }
		if (side == 2) { face_normal.y = -1.0; return float3(inside.x, -half_size.y, inside.z); }
		if (side == 3) { face_normal.y = 1.0; return float3(inside.x, half_size.y, inside.z); }
		if (side == 4) { face_normal.z = 1.0; return float3(inside.x, inside.y, half_size.z); }
		face_normal.z = -1.0;
		return float3(inside.x, inside.y, -half_size.z);
	}
	if (block_d.shape == tfx_emission_ellipse) {
		if ((block_a.flags & tfx_spawn_uniform_distribution) != 0) {
			const uint sequence_index = block_a.uid_base + spawn_index + 1;
			unit_u = tfx_shifted_sequence(sequence_index, 2, block_a.uid_base, tfx_spawn_axis_sequence_theta);
			unit_v = tfx_shifted_sequence(sequence_index, 3, block_a.uid_base, tfx_spawn_axis_sequence_phi);
			unit_w = tfx_shifted_sequence(sequence_index, 5, block_a.uid_base, tfx_spawn_axis_sequence_radius);
		}
		const float theta = unit_u * tfx_spawn_pi2;
		const float cos_phi = block_c.arc_size + unit_v * (block_c.arc_offset - block_c.arc_size);
		const float radius = fill_area ? pow(unit_w, 1.0 / 3.0) : 1.0;
		const float sin_phi = sqrt(max(0.0, 1.0 - cos_phi * cos_phi));
		return float3(half_size.x * radius * sin_phi * cos(theta), half_size.y * radius * cos_phi, half_size.z * radius * sin_phi * sin(theta));
	}
	if (block_d.shape == tfx_emission_disc) {
		if ((block_a.flags & tfx_spawn_uniform_distribution) != 0) {
			const uint sequence_index = block_a.uid_base + spawn_index + 1;
			unit_u = tfx_shifted_sequence(sequence_index, 2, block_a.uid_base, tfx_spawn_axis_sequence_theta);
			unit_v = tfx_shifted_sequence(sequence_index, 3, block_a.uid_base, tfx_spawn_axis_sequence_radius);
		}
		const float theta = block_c.arc_offset + unit_u * block_c.arc_size;
		const float radius = fill_area ? sqrt(unit_v) : 1.0;
		return float3(half_size.x * radius * cos(theta), 0.0, half_size.z * radius * -sin(theta));
	}
	if (block_d.shape == tfx_emission_cylinder) {
		if (fill_area) {
			const float radius = sqrt(unit_u);
			const float angle = unit_v * tfx_spawn_pi2;
			return float3(half_size.x * radius * cos(angle), unit_w * emitter_size.y, half_size.z * radius * sin(angle));
		}
		const float theta = unit_u * block_c.arc_size + block_c.arc_offset;
		return float3(cos(theta) * half_size.x, unit_v * emitter_size.y, -sin(theta) * half_size.z);
	}
	return float3(0.0, 0.0, 0.0);
}

//tfx__normalize_vec3, exact rather than fast
float3 tfx_normalize_exact(float3 v) {
	const float length = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	return length > 0.0 ? v / length : v;
}

//tfx__ellipse_surface_normal and tfx__cylinder_surface_normal (y = 0, height = 1)
float3 tfx_ellipse_surface_normal(float3 position, float3 half_size) {
	const float3 gradient = 2.0 * position / (half_size * half_size);
	const float length = 1.0 / tfx_quake_inverse_sqrt(dot(gradient, gradient));
	return gradient / length;
}

//tfx__generate_shape_direction
float3 tfx_generate_shape_direction(tfx_gpu_spawn_block_a block_a, tfx_gpu_spawn_block_c block_c, tfx_gpu_spawn_block_d block_d, float3 spawn_offset, float3 shape_position, float3 shape_size, uint spawn_index, float3 face_normal) {
	const uint mode = block_d.emission_direction_mode;
	const float3 radial = spawn_offset.x == 0.0 && spawn_offset.y == 0.0 && spawn_offset.z == 0.0 ? shape_size : spawn_offset;
	float3 to_handle = float3(0.0, 1.0, 0.0);
	if (mode == tfx_direction_outwards) {
		to_handle = tfx_normalize_exact(radial);
	} else if (mode == tfx_direction_inwards) {
		to_handle = tfx_normalize_exact(-radial);
	} else if (mode == tfx_direction_bothways) {
		to_handle = tfx_normalize_exact(((block_c.spawn_ordinal + spawn_index) & 1) != 0 ? radial : -radial);
	} else if (mode == tfx_direction_surface) {
		if (block_d.shape == tfx_emission_ellipse) {
			to_handle = tfx_ellipse_surface_normal(shape_position, shape_size * 0.5);
		} else if (block_d.shape == tfx_emission_cylinder) {
			const float3 normal = tfx_ellipse_surface_normal(float3(shape_position.x, 0.0, shape_position.z), float3(shape_size.x * 0.5, 1.0, shape_size.z * 0.5));
			to_handle = float3(normal.x, 0.0, normal.z);
		} else if (block_d.shape == tfx_emission_area && (block_a.flags & tfx_spawn_fill_area) == 0) {
			to_handle = face_normal;
		}
	}
	return tfx_rotate_vector(to_handle, block_d.emission_rotation);
}

//tfx__generate_gpu_spawn
tfx_generated_spawn tfx_generate_spawn(tfx_gpu_spawn_block_a block_a, tfx_gpu_spawn_block_b block_b, tfx_gpu_spawn_block_c block_c, tfx_gpu_spawn_block_d block_d, tfx_gpu_spawn_block_e block_e, tfx_gpu_spawn_run run, tfx_gpu_spawn_multipliers multipliers, bool has_multipliers, uint spawn_index) {
	const uint flags = block_a.flags;
	const bool relative = (flags & tfx_spawn_relative) != 0;
	const bool rotate_from_location = !relative && (run.flags & tfx_spawn_run_rotated) != 0;
	tfx_generated_spawn spawn;
	spawn.uid = tfx_seedgen(block_a.uid_base + spawn_index);
	const uint uid = spawn.uid;
	//Spread over the whole tick per location, so a moving location doesn't leave its particles bunched together
	spawn.age = block_b.frame_length * float(spawn_index - run.first_spawn) / float(run.count);

	spawn.inv_max_age = 1.0 / max(block_b.life + tfx_hash_unit(uid, tfx_spawn_axis_life) * block_b.life_variation, 1.0);
	if (has_multipliers && (run.flags & tfx_spawn_run_scale_life) != 0) {
		spawn.inv_max_age = min(spawn.inv_max_age * multipliers.inverse_life, block_b.inverse_frame_length);
	}

	float3 anchor = float3(0.0, 0.0, 0.0);
	if (!relative) {
		const float tween = block_b.frame_length > 0.0 ? 1.0 - spawn.age / block_b.frame_length : 1.0;
		anchor = run.position * tween + run.captured_position * (1.0 - tween);
	}
	float3 splat = float3(0.0, 0.0, 0.0);
	if (block_a.splatter != 0.0) {
		//Uniform in the ball, the same distribution as the CPU's rejection sampling without the loop
		const float z = 1.0 - 2.0 * tfx_hash_unit(uid, tfx_spawn_axis_splatter_z);
		const float phi = tfx_spawn_pi2 * tfx_hash_unit(uid, tfx_spawn_axis_splatter_phi);
		const float radius = pow(tfx_hash_unit(uid, tfx_spawn_axis_splatter_radius), 1.0 / 3.0) * block_a.splatter;
		const float ring = sqrt(max(1.0 - z * z, 0.0)) * radius;
		splat = float3(ring * cos(phi), ring * sin(phi), z * radius);
		if (has_multipliers) {
			splat = splat * multipliers.splatter;
		}
	}
	float4 location_rotation = float4(0.0, 0.0, 0.0, 1.0);
	if (rotate_from_location) {
		location_rotation = tfx_unpack16bit_quaternion(run.packed_rotation);
	}

	float3 direction;
	if ((flags & tfx_spawn_shape) == 0) {
		spawn.position = relative ? float3(0.0, 0.0, 0.0) : anchor + block_a.handle;
		if (block_a.splatter != 0.0) {
			spawn.position += relative ? splat : splat * block_b.overall_scale;
		}
		direction = block_a.emission_direction;
	} else {
		float3 face_normal;
		float3 local_offset = tfx_generate_shape_offset(block_a, block_c, block_d, block_e, uid, spawn_index, face_normal);
		float3 shape_size = block_d.emitter_size;
		if (has_multipliers) {
			const float3 shape_scale = float3(multipliers.emitter_width, multipliers.emitter_height, multipliers.emitter_depth);
			local_offset = local_offset * shape_scale;
			shape_size = shape_size * shape_scale;
		}
		const float3 shape_position = local_offset + splat;
		const float3 spawn_offset = shape_position + block_d.local_handle;
		if (relative) {
			spawn.position = shape_position;
		} else {
			float4 shape_rotation = block_d.shape_rotation;
			if (rotate_from_location) {
				shape_rotation = tfx_multiply_quaternions(location_rotation, shape_rotation);
			}
			spawn.position = anchor + tfx_rotate_vector(spawn_offset, shape_rotation) * block_b.overall_scale;
		}
		direction = tfx_generate_shape_direction(block_a, block_c, block_d, spawn_offset, shape_position, shape_size, spawn_index, face_normal);
	}
	if ((flags & tfx_spawn_stepped_emission) != 0) {
		direction = tfx_generate_stepped_direction(block_e, run, uid, spawn_index, direction);
	}
	if ((flags & tfx_spawn_emission_cone) != 0) {
		const float z = 1.0 - (1.0 - block_a.emission_cone_cosine) * tfx_hash_unit(uid, tfx_spawn_axis_cone_z);
		const float phi = tfx_spawn_pi2 * tfx_hash_unit(uid, tfx_spawn_axis_cone_phi);
		direction = tfx_vector_in_cone(direction, z, phi);
	}
	direction = tfx_normalize_fast(direction);
	if (rotate_from_location) {
		direction = tfx_rotate_vector(direction, location_rotation);
	}
	spawn.velocity_normal = tfx_pack10bit_unsigned(direction);

	spawn.velocity = float3(0.0, 0.0, 0.0);
	if (block_c.spawn_impulse != 0.0 || block_c.spawn_impulse_variation != 0.0) {
		const float impulse = block_c.spawn_impulse + tfx_white_unit(uid, tfx_spawn_impulse_hash_axis) * block_c.spawn_impulse_variation;
		spawn.velocity = tfx_unpack10bit(spawn.velocity_normal) * impulse;
	}

	spawn.base_weight = 0.0;
	if (block_b.weight != 0.0) {
		spawn.base_weight = block_b.weight;
		if (block_b.weight_variation > 0.0) {
			spawn.base_weight += (2.0 * tfx_hash_unit(uid, tfx_spawn_axis_weight) - 1.0) * block_b.weight_variation;
		}
	}

	const float random_width = tfx_hash_unit(uid, tfx_spawn_axis_width) * block_b.size_variation.x;
	const float random_height = (flags & tfx_spawn_base_uniform_size) != 0 ? random_width : tfx_hash_unit(uid, tfx_spawn_axis_height) * block_b.size_variation.y;
	spawn.base_size = float2(random_width + block_b.size.x, random_height + block_b.size.y);

	spawn.base_roll_spin = (2.0 * tfx_hash_unit(uid, tfx_spawn_axis_roll_spin) - 1.0) * block_c.roll_spin_variation + block_c.roll_spin;
	spawn.base_pitch_spin = 0.0;
	spawn.base_yaw_spin = 0.0;
	if ((flags & tfx_spawn_spin_pitch_and_yaw) != 0) {
		spawn.base_yaw_spin = (2.0 * tfx_hash_unit(uid, tfx_spawn_axis_yaw_spin) - 1.0) * block_c.yaw_spin_variation + block_c.yaw_spin;
		spawn.base_pitch_spin = (2.0 * tfx_hash_unit(uid, tfx_spawn_axis_pitch_spin) - 1.0) * block_c.pitch_spin_variation + block_c.pitch_spin;
		float3 offsets;
		offsets.x = (flags & tfx_spawn_specify_pitch) != 0 ? block_b.angle_offsets.x : tfx_hash_unit(uid, tfx_spawn_axis_pitch_offset) * block_b.angle_offsets.x;
		offsets.y = (flags & tfx_spawn_specify_yaw) != 0 ? block_b.angle_offsets.y : tfx_hash_unit(uid, tfx_spawn_axis_yaw_offset) * block_b.angle_offsets.y;
		offsets.z = (flags & tfx_spawn_specify_roll) != 0 ? block_b.angle_offsets.z : tfx_hash_unit(uid, tfx_spawn_axis_roll_offset) * block_b.angle_offsets.z;
		if (rotate_from_location && (flags & tfx_spawn_relative_angle) != 0) {
			//The shader composes conj(emitter) * particle, and the location sits in front of the emitter
			const float4 inverse_rotation = float4(-location_rotation.xyz, location_rotation.w);
			offsets = tfx_quaternion_to_spin_euler(tfx_multiply_quaternions(inverse_rotation, tfx_spin_euler_to_quaternion(offsets)));
		}
		spawn.rotation_offsets = tfx_pack10bit_unsigned(offsets * (1.0 / tfx_spawn_pi2));
	} else {
		float roll = (flags & tfx_spawn_specify_roll) != 0 ? block_b.angle_offsets.z : tfx_hash_unit(uid, tfx_spawn_axis_roll_offset) * block_b.angle_offsets.z;
		if ((flags & tfx_spawn_stepped_roll) != 0) {
			roll = tfx_generate_stepped_roll(block_b, block_e, run, uid, spawn_index);
		}
		if (rotate_from_location && (flags & tfx_spawn_relative_angle) != 0) {
			roll += run.location_roll;
		}
		spawn.rotation_offsets = asuint(roll);
	}

	spawn.base_velocity = block_b.velocity;
	spawn.intensity_factor = 1.0;
	spawn.noise_offset = float3(0.0, 0.0, 0.0);
	if ((flags & tfx_spawn_noise_offsets) != 0) {
		spawn.noise_offset.x = tfx_hash_unit(uid, tfx_spawn_axis_noise_x) * block_c.noise_offset_variation + block_c.noise_shared_offset + block_a.noise_step_offset.x;
		spawn.noise_offset.y = tfx_hash_unit(uid, tfx_spawn_axis_noise_y) * block_c.noise_offset_variation + block_c.noise_shared_offset + block_a.noise_step_offset.y;
		spawn.noise_offset.z = tfx_hash_unit(uid, tfx_spawn_axis_noise_z) * block_c.noise_offset_variation + block_c.noise_shared_offset + block_a.noise_step_offset.z;
	}
	if (has_multipliers) {
		spawn.base_size.x *= multipliers.particle_width;
		spawn.base_size.y *= multipliers.particle_height;
		spawn.base_velocity *= multipliers.velocity;
		spawn.base_roll_spin *= multipliers.spin;
		spawn.base_pitch_spin *= multipliers.spin;
		spawn.base_yaw_spin *= multipliers.spin;
		spawn.intensity_factor *= multipliers.intensity;
		spawn.base_weight *= multipliers.weight;
		if ((flags & tfx_spawn_noise_offsets) != 0) {
			spawn.noise_offset += multipliers.noise_offset;
		}
	}

	spawn.image_frame = (flags & tfx_spawn_random_start_frame) != 0 ? tfx_hash_unit(uid, tfx_spawn_axis_start_frame) * block_b.animation_frames : block_c.start_frame;
	spawn.random_color = (flags & tfx_spawn_random_color) != 0 ? tfx_hash_unit(uid, tfx_spawn_axis_random_color) : 0.0;
	spawn.spawn_location = relative ? run.local_id : 0xFFFFFFFFu;
	spawn.noise_resolution = (flags & tfx_spawn_noise_resolution) != 0 ? block_a.noise_resolution : 0.0;
	return spawn;
}
