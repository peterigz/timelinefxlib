//The reference GPU particle kernel. One dispatch per tfx_gpu_particle_dispatch_t from tfx_GetGPUParticleFrame, one thread per ring
//slot. Each thread runs every tick record of the fetch in registers, mirroring the CPU's control, instance and age passes in their
//order, and writes its instance on the last tick. Ring positions count from the tick's head, so a slot that is bumped and reused
//within one fetch just ingests its new occupant at the tick it was spawned: from a spawn row, or generated from the tick's spawn block
//when the tick is flagged generated (tfx_gpu_spawn.slang).

import tfx_graphs;
import tfx_noise;
//Exported so a host's ITfxParticleSpawn sees tfx_generated_spawn and the hash draws
__exported import tfx_gpu_spawn;

static const uint tfx_group_size = 256;

//Each particle emitter's graphs sit in emitter graph list order from its graph_offset, see tfx__update_library_compute_nodes
static const uint tfx_graph_intensity = 0;
static const uint tfx_graph_alpha_sharpness = 1;
static const uint tfx_graph_curved_alpha = 2;
static const uint tfx_graph_gradient_mapper = 3;
static const uint tfx_graph_velocity = 4;
static const uint tfx_graph_width = 5;
static const uint tfx_graph_height = 6;
static const uint tfx_graph_weight = 7;
static const uint tfx_graph_pitch_spin = 8;
static const uint tfx_graph_yaw_spin = 9;
static const uint tfx_graph_roll_spin = 10;
static const uint tfx_graph_stretch = 11;
static const uint tfx_graph_velocity_turbulance = 12;
static const uint tfx_graph_noise_resolution = 15;
static const uint tfx_graph_motion_randomness = 16;

//tfx_gpu_particle_dispatch_flag_bits
static const uint tfx_dispatch_capture = 1 << 0;
static const uint tfx_dispatch_align_emission = 1 << 1;
static const uint tfx_dispatch_align_emitter = 1 << 2;
static const uint tfx_dispatch_random_color = 1 << 3;
static const uint tfx_dispatch_lifetime_uniform_size = 1 << 4;
static const uint tfx_dispatch_base_uniform_size = 1 << 5;
static const uint tfx_dispatch_do_not_render = 1 << 6;
static const uint tfx_dispatch_spin_roll = 1 << 7;
static const uint tfx_dispatch_spin_3d = 1 << 8;
static const uint tfx_dispatch_spin_pitch_and_yaw = 1 << 9;
static const uint tfx_dispatch_relative_angle = 1 << 10;
static const uint tfx_dispatch_play_once = 1 << 11;
static const uint tfx_dispatch_reverse_animation = 1 << 12;
static const uint tfx_dispatch_drag_variation = 1 << 13;
static const uint tfx_dispatch_relative_location = 1 << 14;
static const uint tfx_dispatch_location_rotation = 1 << 16;
static const uint tfx_dispatch_align_rotated_emission = 1 << 17;
static const uint tfx_dispatch_noise = 1 << 18;
static const uint tfx_dispatch_noise_algorithm_shift = 19;
static const uint tfx_dispatch_noise_anisotropic = 1 << 22;
static const uint tfx_dispatch_location_forces = 1 << 23;
static const uint tfx_dispatch_location_forces_rotation = 1 << 24;

//The noise algorithm bits, tfx_noise_type - tfxWhiteNoise
static const uint tfx_noise_white = 0;
static const uint tfx_noise_simplex = 1;
static const uint tfx_noise_curl = 2;
static const uint tfx_noise_value = 3;

//tfx_gpu_particle_flag_bits
static const uint tfx_particle_remove = 1 << 0;

//tfx_gpu_particle_tick_flag_bits
static const uint tfx_tick_generated = 1 << 0;

//tfx_gpu_particle_force_type in type_flags & tfx_force_type_mask, tfx_gpu_particle_force_flag_bits above it
static const uint tfx_force_type_mask = 0xFF;
static const uint tfx_force_wind = 0;
static const uint tfx_force_attract = 1;
static const uint tfx_force_vortex = 2;
static const uint tfx_force_shockwave = 3;
static const uint tfx_force_noise = 4;
static const uint tfx_force_profile_bezier = 1 << 8;
static const uint tfx_force_bounded = 1 << 9;
static const uint tfx_force_fluctuating = 1 << 10;
static const uint tfx_force_at_location = 1 << 11;

//Must match tfxDRAG_VARIATION_HASH_AXIS or the drag drawn here won't be the one the CPU draws for the same uid
static const uint tfx_drag_variation_hash_axis = 0x9e3779b9u;
static const float tfx_size_packed_scale = 32767.0 / 64.0;
static const float tfx_color_packed_scale = 32767.0 / 128.0;
static const uint tfx_capture_bit = 1 << 15;
//tfxSPAWN_LOCATION_SLOT_MASK
static const uint tfx_spawn_location_slot_mask = (1u << 20) - 1;

//Mirrors tfx_gpu_particle_hot_t
struct tfx_gpu_particle_hot {
	float3 position;
	uint flags;
	float3 velocity;
	float age;							//Milliseconds
};

//Mirrors tfx_gpu_particle_cold_t
struct tfx_gpu_particle_cold {
	uint uid;
	float inv_max_age;
	uint velocity_normal;
	float base_velocity;
	float base_weight;
	float base_size_x;
	float base_size_y;
	float intensity_factor;
	float random_color;
	uint rotation_offsets;
	float base_roll_spin;
	float base_pitch_spin;
	float base_yaw_spin;
	float image_start_frame;
	uint spawn_location;
	uint padding;
};

//Mirrors tfx_gpu_particle_noise_t
struct tfx_gpu_particle_noise {
	float3 offset;						//Authored units, zero for white noise
	float resolution;
};

//Mirrors tfx_gpu_particle_spawn_t
struct tfx_gpu_particle_spawn {
	tfx_gpu_particle_hot hot;
	tfx_gpu_particle_cold cold;
	tfx_gpu_particle_noise noise;
};

//Mirrors tfx_gpu_particle_tick_t
struct tfx_gpu_particle_tick {
	uint head;
	uint count;
	uint spawn_count;
	uint upload_offset;
	float frame_length;					//Milliseconds
	float update_time;					//Seconds
	float overall_scale;
	float velocity_adjuster;
	float drag_alpha;
	float acceleration_scale;
	float drag_alpha_variation;
	float acceleration_scale_variation;
	uint force_offset;					//Row of the first force row in the same buffer
	uint force_count;
	uint flags;							//tfxGPUParticleTickFlags
	uint padding;
};

//Mirrors tfx_gpu_particle_dispatch_t
struct tfx_gpu_particle_dispatch {
	uint base;
	uint capacity;
	uint first_tick;
	uint tick_count;
	uint slot_start;
	uint slot_count;
	uint flags;
	uint graph_offset;
	uint property_index;
	uint instance_start;
	uint previous_instance_start;
	uint bumped_since_write;
	uint spawned_since_write;
	float global_intensity;
	float global_stretch;
	float image_frame_rate;
	float end_frame;
	float emitter_roll;
	uint location_offset;
	uint location_count;
	float4 emitter_rotation;			//Quaternion xyzw
	float3 handle;
	float noise_swerve_gain_minus_one;
	float noise_speed_gain_minus_swerve_gain;
	uint key;							//The host's tfx_SetTemplateEmitterGPUKey / tfx_SetEffectGPUKey, which the library never reads
	uint padding_1;
	uint padding_2;
};

//Mirrors tfx_gpu_spawn_location_t
struct tfx_gpu_spawn_location {
	float3 position;
	uint local_id;						//slot | generation, tfxINVALID when the slot is empty
	uint2 packed_rotation;				//16 bit snorm: .x = X|Y, .y = Z|W
	uint2 padding;
};

//Mirrors tfx_gpu_particle_force_t. Every float3 starts on a 16 byte offset, so std430 matches the C layout
struct tfx_gpu_particle_force {
	uint type_flags;
	float strength;
	float falloff_scale;
	float front;
	float3 origin;
	float radial_ratio;
	float3 axis;						//Vortex axis, or the summed wind flow
	float axial_ratio;
	float4 profile;						//from, curve1, curve2, to
};

//Mirrors tfx_gpu_particle_vortex_fluctuation_t, the row after a fluctuating vortex
struct tfx_gpu_particle_vortex_fluctuation {
	float3 reference;
	uint lobes;
	float3 reference_turn;
	float padding_turn;
	float3 harmonic_cos;
	float padding_cos;
	float3 harmonic_sin;
	float padding_sin;
};

//Mirrors tfx_instance_t
struct tfx_instance {
	float4 position;					//Stretch in w
	uint2 quaternion;					//16 bit snorm: .x = X|Y, .y = Z|W
	uint size;
	uint alignment;
	uint intensity_gradient_map;
	uint curved_alpha_life;
	uint indexes;						//[gpu property index << 16, capture flag << 15, image frame]
	uint captured_index;
};

//Mirrors tfx_gpu_particle_push_t
struct tfx_gpu_particle_push {
	uint dispatch_index;
	uint hot_buffer_index;
	uint cold_buffer_index;
	uint noise_buffer_index;
	uint instance_buffer_index;
	uint graph_buffer_index;
	uint dispatch_buffer_index;
	uint tick_buffer_index;
	uint spawn_buffer_index;
	uint location_buffer_index;
	uint cold_base;						//Each base is the row its section starts at in its buffer
	uint noise_base;
	uint tick_base;
	uint spawn_base;
	uint location_base;
	uint padding;
};

//The pool's hot, cold and noise sections can be views of one buffer, as can the frame's arrays
[[vk::binding(5, 0)]] RWStructuredBuffer<tfx_gpu_particle_hot> hot_buffers[];
[[vk::binding(5, 0)]] RWStructuredBuffer<tfx_gpu_particle_cold> cold_buffers[];
[[vk::binding(5, 0)]] RWStructuredBuffer<tfx_gpu_particle_noise> noise_buffers[];
[[vk::binding(5, 0)]] RWStructuredBuffer<tfx_instance> instance_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_particle_spawn> spawn_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_particle_dispatch> dispatch_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_particle_tick> tick_buffers[];
//The force rows after the tick records, the same buffer as tick_buffers
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_particle_force> force_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_particle_vortex_fluctuation> fluctuation_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_graph_data> graph_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_spawn_location> location_buffers[];
//A generated tick's spawn block, runs and multipliers, 64 byte rows of the same buffer as tick_buffers
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_spawn_block_a> spawn_block_a_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_spawn_block_b> spawn_block_b_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_spawn_block_c> spawn_block_c_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_spawn_block_d> spawn_block_d_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_spawn_block_e> spawn_block_e_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_spawn_run> spawn_run_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_spawn_multipliers> spawn_multipliers_buffers[];

//The pack helpers truncate like _mm_cvttps_epi32, and mask the two's complement bits as the CPU does
uint tfx_pack16bit_2sscaled(float x, float y, float max_value) {
	return (uint(int(x * max_value)) & 0xFFFFu) | ((uint(int(y * max_value)) & 0xFFFFu) << 16);
}

uint tfx_pack8bit_xyz(float3 v) {
	return (uint(int(v.x * 127.0)) & 0xFFu) | ((uint(int(v.y * 127.0)) & 0xFFu) << 8) | ((uint(int(v.z * 127.0)) & 0xFFu) << 16);
}

uint tfx_pack8bitunorm_xyz(float3 v) {
	return (uint(int(v.x * 255.0)) & 0xFFu) | ((uint(int(v.y * 255.0)) & 0xFFu) << 8) | ((uint(int(v.z * 255.0)) & 0xFFu) << 16);
}

uint2 tfx_pack16bit_quaternion(float4 q) {
	const uint x = uint(int(q.x * 32767.0)) & 0xFFFFu;
	const uint y = uint(int(q.y * 32767.0)) & 0xFFFFu;
	const uint z = uint(int(q.z * 32767.0)) & 0xFFFFu;
	const uint w = uint(int(q.w * 32767.0)) & 0xFFFFu;
	return uint2(x | (y << 16), z | (w << 16));
}

//The euler order of tfx__block_particle_spin, xyzw
float4 tfx_euler_to_quaternion(float pitch, float yaw, float roll) {
	const float sin_pitch = sin(pitch * 0.5), cos_pitch = cos(pitch * 0.5);
	const float sin_yaw = sin(yaw * 0.5), cos_yaw = cos(yaw * 0.5);
	const float sin_roll = sin(roll * 0.5), cos_roll = cos(roll * 0.5);
	const float cos_pitch_cos_yaw = cos_pitch * cos_yaw;
	const float sin_pitch_sin_yaw = sin_pitch * sin_yaw;
	const float sin_pitch_cos_yaw = sin_pitch * cos_yaw;
	const float cos_pitch_sin_yaw = cos_pitch * sin_yaw;
	return float4(
		cos_roll * sin_pitch_cos_yaw + sin_roll * cos_pitch_sin_yaw,
		cos_roll * cos_pitch_sin_yaw - sin_roll * sin_pitch_cos_yaw,
		sin_roll * cos_pitch_cos_yaw - cos_roll * sin_pitch_sin_yaw,
		cos_roll * cos_pitch_cos_yaw + sin_roll * sin_pitch_sin_yaw);
}

tfx_gpu_graph_data tfx_graph(tfx_gpu_particle_push push, uint graph_offset, uint graph) {
	return graph_buffers[push.graph_buffer_index][graph_offset + graph];
}

//tfx__wide_sample_force_profile: one segment, no easing or oscillator, and the far edge is killed rather than clamped
float tfx_sample_force_profile(tfx_gpu_particle_force force, float profile_time) {
	const bool inside = profile_time < 1.0;
	const float time = min(profile_time, 1.0);
	float sample;
	if ((force.type_flags & tfx_force_profile_bezier) != 0) {
		const float inverse = 1.0 - time;
		const float weight1 = (inverse * inverse) * inverse;
		const float weight2 = ((3.0 * inverse) * inverse) * time;
		const float weight3 = ((3.0 * inverse) * time) * time;
		const float weight4 = (time * time) * time;
		sample = ((weight1 * force.profile.x + weight2 * force.profile.y) + weight3 * force.profile.z) + weight4 * force.profile.w;
	} else {
		sample = force.profile.x + (force.profile.w - force.profile.x) * time;
	}
	return inside ? sample : 0.0;
}

//The distance the profiles are sampled at. Floored rather than branched on, so a particle on the origin reads 0 rather than a NaN
float tfx_force_inverse_length(float length_squared) {
	return rsqrt(max(length_squared, tfx_noise_epsilon));
}

//tfx__wide_apply_vortex_force, with tfx__wide_vortex_fluctuation from the row after it when the force fluctuates
void tfx_apply_vortex_force(tfx_gpu_particle_force force, tfx_gpu_particle_vortex_fluctuation fluctuation, float3 position, inout float3 medium) {
	const float3 axis = force.axis;
	const float3 offset = position - force.origin;
	const float axial_distance = dot(offset, axis);
	float3 radial = offset - axial_distance * axis;
	const float radial_length_squared = dot(radial, radial);
	const float inverse_radial_length = tfx_force_inverse_length(radial_length_squared);
	radial *= inverse_radial_length;
	const float radial_distance = radial_length_squared * inverse_radial_length;

	float scale = force.strength * tfx_sample_force_profile(force, radial_distance * force.falloff_scale);
	if ((force.type_flags & tfx_force_fluctuating) != 0) {
		const float arm_real = dot(radial, fluctuation.reference);
		const float arm_imaginary = dot(radial, fluctuation.reference_turn);
		//The arm raised to the lobe count by squaring, the first factor taken as it is
		float power_real = arm_real;
		float power_imaginary = arm_imaginary;
		bool have_power = false;
		float base_real = arm_real;
		float base_imaginary = arm_imaginary;
		for (uint exponent = fluctuation.lobes; exponent != 0; exponent >>= 1) {
			if ((exponent & 1) != 0) {
				if (have_power) {
					const float next_real = power_real * base_real - power_imaginary * base_imaginary;
					power_imaginary = power_real * base_imaginary + power_imaginary * base_real;
					power_real = next_real;
				} else {
					power_real = base_real;
					power_imaginary = base_imaginary;
					have_power = true;
				}
			}
			if (exponent > 1) {
				const float next_real = base_real * base_real - base_imaginary * base_imaginary;
				base_imaginary = (base_real + base_real) * base_imaginary;
				base_real = next_real;
			}
		}
		//Each further harmonic has one more lobe than the last, so it is one more multiply by the arm
		float variation = 1.0;
		for (uint harmonic_index = 0; harmonic_index != 3; ++harmonic_index) {
			if (harmonic_index > 0) {
				const float next_real = power_real * arm_real - power_imaginary * arm_imaginary;
				power_imaginary = power_real * arm_imaginary + power_imaginary * arm_real;
				power_real = next_real;
			}
			variation = power_real * fluctuation.harmonic_cos[harmonic_index] + variation;
			variation = power_imaginary * fluctuation.harmonic_sin[harmonic_index] + variation;
		}
		scale *= variation;
	}

	const float3 tangential = cross(axis, radial);
	const float3 flow = axis * force.axial_ratio + (radial * force.radial_ratio + tangential);
	medium = flow * scale + medium;
}

//tfx__wide_apply_attract_force: the pull grows with distance towards the edge of the field
void tfx_apply_attract_force(tfx_gpu_particle_force force, float3 position, inout float3 medium) {
	const float3 offset = position - force.origin;
	const float length_squared = dot(offset, offset);
	const float inverse_length = tfx_force_inverse_length(length_squared);
	const float distance = length_squared * inverse_length;
	float scale = force.strength * tfx_sample_force_profile(force, distance * force.falloff_scale);
	scale *= distance;
	medium = (offset * inverse_length) * scale + medium;
}

//tfx__wide_apply_shockwave_force: the profile runs back from the leading edge, so only particles it hasn't reached need masking
void tfx_apply_shockwave_force(tfx_gpu_particle_force force, float3 position, inout float3 medium) {
	const float3 offset = position - force.origin;
	const float length_squared = dot(offset, offset);
	const float inverse_length = tfx_force_inverse_length(length_squared);
	const float distance = length_squared * inverse_length;
	const float profile_time = (force.front - distance) * force.falloff_scale;
	const float scale = profile_time >= 0.0 ? force.strength * tfx_sample_force_profile(force, profile_time) : 0.0;
	medium = (offset * inverse_length) * scale + medium;
}

//The noise force from its row. Like the CPU it samples the field at the start of tick position and splits against the previous velocity.
//A bounded field fades by its distance from bounded_position, the location frame position for an at_location row
float3 tfx_noise_medium(tfx_gpu_particle_push push, tfx_gpu_particle_dispatch dispatch, tfx_gpu_particle_force force, tfx_gpu_particle_tick tick, tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, tfx_gpu_particle_noise noise, float life, float3 bounded_position, bool bounded_alive) {
	const uint flags = dispatch.flags;
	const uint algorithm = (flags >> tfx_dispatch_noise_algorithm_shift) & 3u;
	float3 direction;
	float amplitude;
	if (algorithm == tfx_noise_white) {
		direction = tfx_white_noise_direction(cold.uid, hot.age, noise.resolution);
		amplitude = force.strength * tfx_sample_graph(tfx_graph(push, dispatch.graph_offset, tfx_graph_motion_randomness), life);
	} else {
		const float noise_resolution = (tfx_sample_graph(tfx_graph(push, dispatch.graph_offset, tfx_graph_noise_resolution), life) * noise.resolution) * tick.overall_scale;
		const float3 sample_position = hot.position / noise_resolution + noise.offset * tick.overall_scale;
		if (algorithm == tfx_noise_simplex) {
			direction = tfx_simplex_noise_direction(sample_position, life);
		} else if (algorithm == tfx_noise_curl) {
			direction = tfx_curl_noise_direction(sample_position, life);
		} else {
			direction = tfx_value_noise_direction(sample_position);
		}
		amplitude = force.strength * tfx_sample_graph(tfx_graph(push, dispatch.graph_offset, tfx_graph_velocity_turbulance), life);
	}

	if ((force.type_flags & tfx_force_bounded) != 0) {
		const float3 offset = bounded_position - force.origin;
		const float length_squared = dot(offset, offset);
		const float distance = length_squared * tfx_force_inverse_length(length_squared);
		amplitude *= bounded_alive ? tfx_sample_force_profile(force, distance * force.falloff_scale) : 0.0;
	}
	float3 medium = direction * amplitude;

	if ((flags & tfx_dispatch_noise_anisotropic) != 0) {
		medium = tfx_split_medium_velocity(hot.velocity, medium, dispatch.noise_swerve_gain_minus_one, dispatch.noise_speed_gain_minus_swerve_gain);
	}
	return medium;
}

float4 tfx_conjugate(float4 q) {
	return float4(-q.xyz, q.w);
}

//tfx__block_particle_transform's two spaces: local at a spawn location, and world
struct tfx_ring_transform {
	uint flags;
	float3 handle;
	float4 emitter_rotation;
	uint2 location_rotation;			//16 bit snorm, as tfx_gpu_spawn_location packs it
	float3 location_position;
	uint location_id;					//tfxINVALID when the particle's location slot is out of the table
	float overall_scale;

	float3 to_world(float3 position) {
		if ((flags & tfx_dispatch_relative_location) != 0) {
			float3 world = tfx_rotate_vector(position + handle, emitter_rotation);
			//After the emitter rotation, so the emitter keeps its own angle and the location turns the whole thing into its frame
			if ((flags & tfx_dispatch_location_rotation) != 0) {
				world = tfx_rotate_vector(world, tfx_unpack16bit_quaternion(location_rotation));
			}
			return world * overall_scale + location_position;
		}
		return position;
	}

	float3 to_ring(float3 world) {
		if ((flags & tfx_dispatch_relative_location) != 0) {
			float3 position = (world - location_position) / overall_scale;
			if ((flags & tfx_dispatch_location_rotation) != 0) {
				position = tfx_rotate_vector(position, tfx_conjugate(tfx_unpack16bit_quaternion(location_rotation)));
			}
			return tfx_rotate_vector(position, tfx_conjugate(emitter_rotation)) - handle;
		}
		return world;
	}

	//For velocities and normals: rotation only, because the scale is uniform and velocities are kept in authored units
	float3 direction_to_world(float3 direction) {
		if ((flags & tfx_dispatch_relative_location) == 0) {
			return direction;
		}
		float3 world = tfx_rotate_vector(direction, emitter_rotation);
		if ((flags & tfx_dispatch_location_rotation) != 0) {
			world = tfx_rotate_vector(world, tfx_unpack16bit_quaternion(location_rotation));
		}
		return world;
	}

	float3 direction_to_ring(float3 world) {
		if ((flags & tfx_dispatch_relative_location) == 0) {
			return world;
		}
		float3 direction = world;
		if ((flags & tfx_dispatch_location_rotation) != 0) {
			direction = tfx_rotate_vector(direction, tfx_conjugate(tfx_unpack16bit_quaternion(location_rotation)));
		}
		return tfx_rotate_vector(direction, tfx_conjugate(emitter_rotation));
	}
};

//Only a relative location ring reads the location table
tfx_ring_transform tfx_particle_ring_transform(tfx_gpu_particle_push push, tfx_gpu_particle_dispatch dispatch, tfx_gpu_particle_tick tick, tfx_gpu_particle_cold cold) {
	tfx_ring_transform transform;
	transform.flags = dispatch.flags;
	transform.handle = dispatch.handle;
	transform.emitter_rotation = dispatch.emitter_rotation;
	transform.location_rotation = uint2(0u, 0x7FFF0000u);
	transform.location_position = float3(0.0, 0.0, 0.0);
	transform.location_id = 0xFFFFFFFFu;
	transform.overall_scale = tick.overall_scale;
	if ((dispatch.flags & tfx_dispatch_relative_location) != 0) {
		const uint slot = cold.spawn_location & tfx_spawn_location_slot_mask;
		if (slot < dispatch.location_count) {
			const tfx_gpu_spawn_location location = location_buffers[push.location_buffer_index][push.location_base + dispatch.location_offset + slot];
			transform.location_rotation = location.packed_rotation;
			transform.location_position = location.position;
			transform.location_id = location.local_id;
		}
	}
	return transform;
}

//A relative particle whose location has gone is removed here, before the write, so it's drawn with no size and isn't captured, which
//is what the CPU's location_alive mask does in the instance pass ahead of its age pass kill
float3 tfx_particle_world_position(tfx_gpu_particle_push push, tfx_gpu_particle_dispatch dispatch, tfx_gpu_particle_tick tick, inout tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, out uint2 location_rotation) {
	const tfx_ring_transform transform = tfx_particle_ring_transform(push, dispatch, tick, cold);
	if ((dispatch.flags & tfx_dispatch_relative_location) != 0 && transform.location_id != cold.spawn_location) {
		hot.flags |= tfx_particle_remove;
	}
	location_rotation = transform.location_rotation;
	return transform.to_world(hot.position);
}

//What a particle update sees of the particle's ring and tick
struct tfx_particle_context {
	tfx_gpu_particle_push push;
	tfx_gpu_particle_dispatch dispatch;
	tfx_gpu_particle_tick tick;
	tfx_gpu_particle_noise noise;		//Only loaded for a ring with a noise force, zero otherwise
	float life;
	uint ring_position;
	bool last_tick;						//The tick the instance is written on, after after_move

	tfx_ring_transform ring_transform(tfx_gpu_particle_cold cold) {
		return tfx_particle_ring_transform(push, dispatch, tick, cold);
	}
};

//tfx__wide_apply_location_forces' frame: the particle less its live location, turned back by the location's rotation. The rows have the
//anchor taken off their origin already. A particle whose location has gone feels none of these rows
struct tfx_location_frame {
	float3 position;
	float4 rotation;
	bool alive;
};

tfx_location_frame tfx_particle_location_frame(tfx_particle_context context, tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold) {
	tfx_location_frame frame;
	frame.position = hot.position;
	frame.rotation = float4(0.0, 0.0, 0.0, 1.0);
	frame.alive = false;
	const uint slot = cold.spawn_location & tfx_spawn_location_slot_mask;
	if (slot < context.dispatch.location_count) {
		const tfx_gpu_spawn_location location = location_buffers[context.push.location_buffer_index][context.push.location_base + context.dispatch.location_offset + slot];
		frame.alive = location.local_id == cold.spawn_location;
		frame.position = hot.position - location.position;
		if ((context.dispatch.flags & tfx_dispatch_location_forces_rotation) != 0) {
			frame.rotation = tfx_unpack16bit_quaternion(location.packed_rotation);
			frame.position = tfx_rotate_vector(frame.position, tfx_conjugate(frame.rotation));
		}
	}
	return frame;
}

//tfx_apply_forces over the tick's force rows, which the CPU packed in its order: wind, vortices, attracts, shockwaves, the same at
//locations, noise
float3 tfx_force_medium(tfx_particle_context context, tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold) {
	const uint buffer_index = context.push.tick_buffer_index;
	const uint first_row = context.push.tick_base + context.tick.force_offset;
	const uint end_row = first_row + context.tick.force_count;
	bool location_forces = (context.dispatch.flags & tfx_dispatch_location_forces) != 0;
	tfx_location_frame frame = { hot.position, float4(0.0, 0.0, 0.0, 1.0), false };
	if (location_forces) {
		frame = tfx_particle_location_frame(context, hot, cold);
	}
	float3 medium = float3(0.0, 0.0, 0.0);
	float3 location_medium = float3(0.0, 0.0, 0.0);
	for (uint row = first_row; row < end_row; ++row) {
		const tfx_gpu_particle_force force = force_buffers[buffer_index][row];
		const uint type = force.type_flags & tfx_force_type_mask;
		if (type == tfx_force_wind) {
			//The CPU seeds the medium with the summed wind rather than adding it
			medium = force.axis;
			continue;
		}
		if (type == tfx_force_noise) {
			//The location flow joins before the noise, as on the CPU
			if (location_forces && frame.alive) {
				medium += tfx_rotate_vector(location_medium, frame.rotation);
			}
			const bool noise_at_location = location_forces && (force.type_flags & tfx_force_at_location) != 0;
			location_forces = false;
			medium += tfx_noise_medium(context.push, context.dispatch, force, context.tick, hot, cold, context.noise, context.life, noise_at_location ? frame.position : hot.position, !noise_at_location || frame.alive);
			continue;
		}
		const bool at_location = location_forces && (force.type_flags & tfx_force_at_location) != 0;
		const float3 position = at_location ? frame.position : hot.position;
		float3 target = at_location ? location_medium : medium;
		if (type == tfx_force_vortex) {
			tfx_gpu_particle_vortex_fluctuation fluctuation = {};
			if ((force.type_flags & tfx_force_fluctuating) != 0) {
				fluctuation = fluctuation_buffers[buffer_index][++row];
			}
			tfx_apply_vortex_force(force, fluctuation, position, target);
		} else if (type == tfx_force_attract) {
			tfx_apply_attract_force(force, position, target);
		} else if (type == tfx_force_shockwave) {
			tfx_apply_shockwave_force(force, position, target);
		}
		if (at_location) {
			location_medium = target;
		} else {
			medium = target;
		}
	}
	if (location_forces && frame.alive) {
		medium += tfx_rotate_vector(location_medium, frame.rotation);
	}
	return medium;
}

//Per particle behaviour a host adds with its own entry point over tfx_simulate, such as custom forces or collision. Both run every
//tick in ring space, and neither runs for a removed particle
interface ITfxParticleUpdate {
	//Added to the acceleration before drag and integration, like the built-in forces, from the state at the start of the tick
	float3 acceleration(tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, tfx_particle_context context);
	//After the move and before the instance write. Set tfx_particle_remove to kill the particle
	void after_move(inout tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, tfx_particle_context context);
}

struct tfx_no_particle_update : ITfxParticleUpdate {
	float3 acceleration(tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, tfx_particle_context context) {
		return float3(0.0, 0.0, 0.0);
	}
	void after_move(inout tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, tfx_particle_context context) {
	}
};

//What a spawn hook sees of the particle's ring, tick and spawn location
struct tfx_spawn_context {
	tfx_gpu_particle_push push;
	tfx_gpu_particle_dispatch dispatch;	//key and property_index tell a host's User emitters apart
	tfx_gpu_particle_tick tick;
	uint location_id;					//slot | generation, the low 32 bits of the tfxSpawnLocationID
	uint user_data;						//tfx_SetSpawnLocationUserData
	uint index_in_run;					//Where the particle is among the location's spawns this tick
	uint run_count;
	bool relative;						//Positions are local to the location rather than world space
};

//Places the particles of User emitters the kernel generates, the GPU side of tfx_SetStageSpawnCallback. Spawns the CPU made arrive as
//rows already placed by the callback, so this never sees them
interface ITfxParticleSpawn {
	//After the library has generated the particle as a Point at the location. tfx_hash_unit(spawn.uid, axis) gives deterministic draws on
	//the host's own axes. A spawn impulse has already launched velocity along velocity_normal, so change the two together
	void after_spawn(inout tfx_generated_spawn spawn, tfx_spawn_context context);
}

struct tfx_no_particle_spawn : ITfxParticleSpawn {
	void after_spawn(inout tfx_generated_spawn spawn, tfx_spawn_context context) {
	}
};

//tfx_apply_lookup_velocity, tfx_apply_lookup_weight, tfx_apply_forces, tfx_apply_velocity and tfx_apply_position
void tfx_control_particle<TUpdate : ITfxParticleUpdate>(tfx_particle_context context, inout tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, tfx_gpu_graph_data velocity_graph, tfx_gpu_graph_data weight_graph, bool drag_variation, TUpdate update) {
	const tfx_gpu_particle_push push = context.push;
	const tfx_gpu_particle_dispatch dispatch = context.dispatch;
	const tfx_gpu_particle_tick tick = context.tick;
	const float life = context.life;
	const float velocity = tfx_sample_graph(velocity_graph, life) * cold.base_velocity;
	const float weight = tfx_sample_graph(weight_graph, life) * cold.base_weight;
	float3 acceleration = tfx_unpack10bit(cold.velocity_normal) * velocity;
	//A force before its delay, or with nothing to describe this tick, has no row
	if (tick.force_count != 0) {
		acceleration += tfx_force_medium(context, hot, cold);
	}
	acceleration += update.acceleration(hot, cold, context);
	acceleration.y -= weight;

	//A particle that has only existed for part of the frame takes that part of the acceleration and drag. A 0ms tick gives the
	//CPU's _mm_min_ps 1 for NaN, where FMin with a NaN is undefined, so it's handled here
	const float age_fraction = tick.frame_length > 0.0 ? min(hot.age / tick.frame_length, 1.0) : 1.0;
	float drag_alpha = tick.drag_alpha;
	float acceleration_scale = tick.acceleration_scale;
	if (drag_variation) {
		const float draw = tfx_white_unit(cold.uid, tfx_drag_variation_hash_axis) * 0.5 + 0.5;
		drag_alpha += tick.drag_alpha_variation * draw;
		acceleration_scale += tick.acceleration_scale_variation * draw;
	}
	drag_alpha *= age_fraction;
	acceleration_scale *= age_fraction;
	hot.velocity = (1.0 - drag_alpha) * hot.velocity + acceleration_scale * acceleration;

	hot.position += hot.velocity * tick.update_time * tick.velocity_adjuster * age_fraction * tick.overall_scale;
}

//tfx__control_particle_age for a particle that isn't single
void tfx_age_particle(inout tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, tfx_gpu_particle_tick tick) {
	hot.age += tick.frame_length;
	if (hot.age * cold.inv_max_age >= 1.0) {
		hot.flags |= tfx_particle_remove;
	}
}

//tfx__block_particle_spin
//q = conj(a) * q
float4 tfx_conjugate_multiply(float4 a, float4 q) {
	return float4(
		(a.w * q.x - a.x * q.w) + (a.z * q.y - a.y * q.z),
		(a.w * q.y - a.y * q.w) + (a.x * q.z - a.z * q.x),
		(a.w * q.z - a.z * q.w) + (a.y * q.x - a.x * q.y),
		(a.w * q.w + a.x * q.x) + (a.y * q.y + a.z * q.z));
}

//tfx__wide_location_roll_to_packed_quaternion: the roll plus the location's twist about the roll axis, 2 * atan2(x, w), without the atan2
float4 tfx_location_roll_quaternion(float roll, float4 location) {
	const float sin_roll = sin(roll * 0.5);
	const float cos_roll = cos(roll * 0.5);
	const float length_squared = location.x * location.x + location.w * location.w;
	float twist_w = 1.0;
	float twist_x = 0.0;
	if (length_squared > 0.0) {
		const float inverse_length = rsqrt(max(length_squared, 1e-30)) * 0.999;
		twist_w = location.w * inverse_length;
		twist_x = location.x * inverse_length;
	}
	return float4(0.0, 0.0, sin_roll * twist_w + cos_roll * twist_x, cos_roll * twist_w - sin_roll * twist_x);
}

//location_rotation is only read when a relative location ring's locations have rotations
uint2 tfx_particle_quaternion(tfx_gpu_particle_push push, tfx_gpu_particle_dispatch dispatch, tfx_gpu_particle_cold cold, float life, uint2 location_rotation) {
	const uint flags = dispatch.flags;
	const bool rotate_from_location = (flags & tfx_dispatch_location_rotation) != 0;
	if ((flags & tfx_dispatch_spin_pitch_and_yaw) == 0) {
		float roll = 0.0;
		if ((flags & tfx_dispatch_spin_roll) != 0) {
			roll = tfx_sample_graph(tfx_graph(push, dispatch.graph_offset, tfx_graph_roll_spin), life) * cold.base_roll_spin;
		}
		if ((flags & tfx_dispatch_relative_angle) != 0) {
			roll += dispatch.emitter_roll;
		}
		roll += asfloat(cold.rotation_offsets);
		if ((flags & tfx_dispatch_relative_angle) != 0 && rotate_from_location) {
			return tfx_pack16bit_quaternion(tfx_location_roll_quaternion(roll, tfx_unpack16bit_quaternion(location_rotation)));
		}
		return tfx_pack16bit_quaternion(tfx_euler_to_quaternion(0.0, 0.0, roll));
	}

	const float3 offsets = tfx_unpack10bit(cold.rotation_offsets) * tfx_two_pi;
	float pitch = 0.0;
	float yaw = 0.0;
	float roll = 0.0;
	if ((flags & tfx_dispatch_spin_3d) != 0) {
		pitch = tfx_sample_graph(tfx_graph(push, dispatch.graph_offset, tfx_graph_pitch_spin), life) * cold.base_pitch_spin;
		yaw = tfx_sample_graph(tfx_graph(push, dispatch.graph_offset, tfx_graph_yaw_spin), life) * cold.base_yaw_spin;
		roll = tfx_sample_graph(tfx_graph(push, dispatch.graph_offset, tfx_graph_roll_spin), life) * cold.base_roll_spin;
	}
	float4 q = tfx_euler_to_quaternion(pitch + offsets.x, yaw + offsets.y, roll + offsets.z);

	//conj(emitter) * particle, because the CPU's position transform rotates the opposite way to the vertex shader
	if ((flags & tfx_dispatch_relative_angle) != 0) {
		//The location sits in front of the emitter: conj(location * emitter) = conj(emitter) * conj(location)
		if (rotate_from_location) {
			q = tfx_conjugate_multiply(tfx_unpack16bit_quaternion(location_rotation), q);
		}
		q = tfx_conjugate_multiply(dispatch.emitter_rotation, q);
	}
	return tfx_pack16bit_quaternion(q);
}

//tfx__block_particle_image_frame without the stored frame: the CPU advances it once per tick before using it, so on the tick a
//particle is aged age it reads image_start_frame + rate * (age + update_time)
uint tfx_particle_image_frame(tfx_gpu_particle_dispatch dispatch, tfx_gpu_particle_cold cold, tfx_gpu_particle_hot hot, tfx_gpu_particle_tick tick) {
	float image_frame = cold.image_start_frame + dispatch.image_frame_rate * 0.001 * hot.age + dispatch.image_frame_rate * tick.update_time;
	const float frames = dispatch.end_frame + 1.0;
	const bool play_once = (dispatch.flags & tfx_dispatch_play_once) != 0;
	const bool reverse_animation = (dispatch.flags & tfx_dispatch_reverse_animation) != 0;
	if (reverse_animation && play_once) {
		image_frame = max(dispatch.end_frame - image_frame, 0.0);
	} else if (play_once) {
		image_frame = min(image_frame, dispatch.end_frame);
	} else {
		//tfxWideMod truncates the quotient
		image_frame = image_frame - float(int(image_frame / frames)) * frames;
		if (reverse_animation) {
			image_frame = dispatch.end_frame - image_frame;
		}
	}
	return dispatch.property_index | uint(int(image_frame));
}

//A removed particle still holding its slot: uncaptured, so the vertex shader zeroes size and stretch and only needs a valid quaternion and property index
void tfx_write_empty_instance(tfx_gpu_particle_push push, tfx_gpu_particle_dispatch dispatch, float3 position, uint ring_position) {
	const uint instance_index = dispatch.instance_start + ring_position;
	tfx_instance instance;
	instance.position = float4(position, 0.0);
	instance.quaternion = uint2(0u, 0x7FFF0000u);
	instance.size = 0;
	instance.alignment = 0;
	instance.intensity_gradient_map = 0;
	instance.curved_alpha_life = 0;
	instance.indexes = dispatch.property_index;
	instance.captured_index = instance_index;
	instance_buffers[push.instance_buffer_index][instance_index] = instance;
}

//tfx__control_particle_instances for one particle, at ring position ring_position of the last tick
void tfx_write_particle_instance(tfx_gpu_particle_push push, tfx_gpu_particle_dispatch dispatch, tfx_gpu_particle_tick tick, tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, float3 position, float life, uint ring_position, uint2 location_rotation) {
	const uint flags = dispatch.flags;
	const uint graph_offset = dispatch.graph_offset;
	const bool alive = (hot.flags & tfx_particle_remove) == 0;
	if (!alive) {
		tfx_write_empty_instance(push, dispatch, position, ring_position);
		return;
	}
	tfx_instance instance;

	//----Transform
	const float stretch = tfx_sample_graph(tfx_graph(push, graph_offset, tfx_graph_stretch), life) * dispatch.global_stretch;
	instance.position = float4(position, stretch);
	instance.alignment = 0;
	if ((flags & tfx_dispatch_align_emission) != 0) {
		float3 alignment = tfx_unpack10bit(cold.velocity_normal);
		if ((flags & tfx_dispatch_align_rotated_emission) != 0) {
			alignment = tfx_rotate_vector(alignment, dispatch.emitter_rotation);
		}
		if ((flags & tfx_dispatch_location_rotation) != 0) {
			alignment = tfx_rotate_vector(alignment, tfx_unpack16bit_quaternion(location_rotation));
		}
		instance.alignment = tfx_pack8bit_xyz(alignment);
	} else if ((flags & tfx_dispatch_align_emitter) != 0) {
		float3 alignment = tfx_rotate_vector(float3(0.0, 1.0, 0.0), dispatch.emitter_rotation);
		if ((flags & tfx_dispatch_location_rotation) != 0) {
			alignment = tfx_rotate_vector(alignment, tfx_unpack16bit_quaternion(location_rotation));
		}
		instance.alignment = tfx_pack8bit_xyz(alignment);
	}

	//----Spin
	instance.quaternion = tfx_particle_quaternion(push, dispatch, cold, life, location_rotation);

	//----Size
	instance.size = 0;
	if ((flags & tfx_dispatch_do_not_render) == 0) {
		const float size_life = min(life, 1.0);
		const float lookup_width = tfx_sample_graph(tfx_graph(push, graph_offset, tfx_graph_width), size_life);
		float scale_x = cold.base_size_x * lookup_width;
		float scale_y;
		if ((flags & tfx_dispatch_lifetime_uniform_size) != 0) {
			scale_y = lookup_width * cold.base_size_y;
			if ((flags & tfx_dispatch_base_uniform_size) != 0) {
				scale_y = min(scale_x, scale_y);
			}
		} else {
			scale_y = tfx_sample_graph(tfx_graph(push, graph_offset, tfx_graph_height), size_life) * cold.base_size_y;
		}
		instance.size = tfx_pack16bit_2sscaled(scale_x * tick.overall_scale, scale_y * tick.overall_scale, tfx_size_packed_scale);
	}

	//----Color
	const float intensity = dispatch.global_intensity * tfx_sample_graph(tfx_graph(push, graph_offset, tfx_graph_intensity), life) * cold.intensity_factor;
	const float curved_alpha = tfx_sample_graph(tfx_graph(push, graph_offset, tfx_graph_curved_alpha), life);
	const float alpha_sharpness = tfx_sample_graph(tfx_graph(push, graph_offset, tfx_graph_alpha_sharpness), life);
	const float gradient_map = (flags & tfx_dispatch_random_color) != 0 ? cold.random_color : tfx_sample_graph(tfx_graph(push, graph_offset, tfx_graph_gradient_mapper), life);
	instance.intensity_gradient_map = tfx_pack16bit_2sscaled(intensity, gradient_map, tfx_color_packed_scale);
	instance.curved_alpha_life = tfx_pack8bitunorm_xyz(float3(curved_alpha, alpha_sharpness, life));

	//----Image frame and interpolation. Only particles that were in the ring at the previous write have somewhere to come from
	const uint instance_index = dispatch.instance_start + ring_position;
	const bool capture = (flags & tfx_dispatch_capture) != 0 && ring_position + dispatch.spawned_since_write < tick.count;
	instance.indexes = tfx_particle_image_frame(dispatch, cold, hot, tick) | (capture ? tfx_capture_bit : 0);
	instance.captured_index = capture ? dispatch.previous_instance_start + ring_position + dispatch.bumped_since_write : instance_index;

	instance_buffers[push.instance_buffer_index][instance_index] = instance;
}

//A spawn the kernel builds from its tick's block rather than loading a row, see tfx__record_gpu_groups
void tfx_generate_particle<TSpawn : ITfxParticleSpawn>(tfx_gpu_particle_push push, tfx_gpu_particle_dispatch dispatch, tfx_gpu_particle_tick tick, uint spawn_index, TSpawn spawn_hook, out tfx_gpu_particle_hot hot, out tfx_gpu_particle_cold cold, out tfx_gpu_particle_noise noise) {
	const uint buffer_index = push.tick_buffer_index;
	const uint block_row = push.tick_base + tick.upload_offset;
	const tfx_gpu_spawn_block_a block_a = spawn_block_a_buffers[buffer_index][block_row];
	const tfx_gpu_spawn_block_b block_b = spawn_block_b_buffers[buffer_index][block_row + 1];
	const tfx_gpu_spawn_block_c block_c = spawn_block_c_buffers[buffer_index][block_row + 2];
	tfx_gpu_spawn_block_d block_d = {};
	if ((block_a.flags & tfx_spawn_shape) != 0) {
		block_d = spawn_block_d_buffers[buffer_index][block_row + 3];
	}
	tfx_gpu_spawn_block_e block_e = {};
	if ((block_a.flags & (tfx_spawn_grid | tfx_spawn_stepped_emission | tfx_spawn_stepped_roll)) != 0) {
		block_e = spawn_block_e_buffers[buffer_index][block_row + 4];
	}
	//The last run starting at or before the spawn. Neighbouring threads are in the same run, so the search stays coherent
	const uint first_run = block_row + block_a.run_row;
	uint low = 0;
	uint high = block_a.run_count - 1;
	while (low < high) {
		const uint middle = (low + high + 1) >> 1;
		if (spawn_run_buffers[buffer_index][first_run + middle].first_spawn <= spawn_index) {
			low = middle;
		} else {
			high = middle - 1;
		}
	}
	const tfx_gpu_spawn_run run = spawn_run_buffers[buffer_index][first_run + low];
	const bool has_multipliers = run.multipliers_row != 0xFFFFFFFFu;
	tfx_gpu_spawn_multipliers multipliers = {};
	if (has_multipliers) {
		multipliers = spawn_multipliers_buffers[buffer_index][block_row + run.multipliers_row];
	}
	tfx_generated_spawn spawn = tfx_generate_spawn(block_a, block_b, block_c, block_d, block_e, run, multipliers, has_multipliers, spawn_index);
	if ((block_a.flags & tfx_spawn_user) != 0) {
		tfx_spawn_context context;
		context.push = push;
		context.dispatch = dispatch;
		context.tick = tick;
		context.location_id = run.local_id;
		context.user_data = run.user_data;
		context.index_in_run = spawn_index - run.first_spawn;
		context.run_count = run.count;
		context.relative = (block_a.flags & tfx_spawn_relative) != 0;
		spawn_hook.after_spawn(spawn, context);
	}
	hot.position = spawn.position;
	hot.flags = 0;
	hot.velocity = spawn.velocity;
	hot.age = spawn.age;
	cold.uid = spawn.uid;
	cold.inv_max_age = spawn.inv_max_age;
	cold.velocity_normal = spawn.velocity_normal;
	cold.base_velocity = spawn.base_velocity;
	cold.base_weight = spawn.base_weight;
	cold.base_size_x = spawn.base_size.x;
	cold.base_size_y = spawn.base_size.y;
	cold.intensity_factor = spawn.intensity_factor;
	cold.random_color = spawn.random_color;
	cold.rotation_offsets = spawn.rotation_offsets;
	cold.base_roll_spin = spawn.base_roll_spin;
	cold.base_pitch_spin = spawn.base_pitch_spin;
	cold.base_yaw_spin = spawn.base_yaw_spin;
	//The frame at age 0, as tfx__pack_gpu_group sets it from the bank
	const float image_frame_rate_per_millisecond = dispatch.image_frame_rate * 0.001;
	cold.image_start_frame = spawn.image_frame - image_frame_rate_per_millisecond * spawn.age;
	cold.spawn_location = spawn.spawn_location;
	cold.padding = 0;
	noise.offset = spawn.noise_offset;
	noise.resolution = spawn.noise_resolution;
}

//position % capacity for a position below 2 * capacity, which every ring index sum here is, without an integer divide
uint tfx_wrap_ring(uint position, uint capacity) {
	return position >= capacity ? position - capacity : position;
}

//The kernel for one ring, push.dispatch_index, with a thread per slot. Dispatch in groups of tfx_group_size
void tfx_simulate<TUpdate : ITfxParticleUpdate, TSpawn : ITfxParticleSpawn>(uint thread_index, tfx_gpu_particle_push push, TUpdate update, TSpawn spawn_hook) {
	const tfx_gpu_particle_dispatch dispatch = dispatch_buffers[push.dispatch_buffer_index][push.dispatch_index];
	if (thread_index >= dispatch.slot_count) {
		return;
	}
	const uint ring_slot = tfx_wrap_ring(dispatch.slot_start + thread_index, dispatch.capacity);
	const uint pool_slot = dispatch.base + ring_slot;
	const tfx_gpu_graph_data velocity_graph = tfx_graph(push, dispatch.graph_offset, tfx_graph_velocity);
	const tfx_gpu_graph_data weight_graph = tfx_graph(push, dispatch.graph_offset, tfx_graph_weight);
	const bool drag_variation = (dispatch.flags & tfx_dispatch_drag_variation) != 0;
	const bool has_noise = (dispatch.flags & tfx_dispatch_noise) != 0;

	tfx_gpu_particle_hot hot;
	tfx_gpu_particle_cold cold;
	tfx_gpu_particle_noise noise = {};
	bool live = false;
	bool ingested = false;
	//Removed before this fetch and waiting for its batch to be bumped: cold is never loaded and hot is never stored
	bool dead = false;
	const uint end_tick = dispatch.first_tick + dispatch.tick_count;
	for (uint tick_index = dispatch.first_tick; tick_index != end_tick; ++tick_index) {
		const tfx_gpu_particle_tick tick = tick_buffers[push.tick_buffer_index][push.tick_base + tick_index];
		const uint ring_position = tfx_wrap_ring(ring_slot + dispatch.capacity - tick.head, dispatch.capacity);
		if (ring_position >= tick.count) {
			live = false;
			continue;
		}
		//Older particles only ever move towards the head, so the slots past count - spawn_count can only hold this tick's spawns
		const uint first_spawn = tick.count - tick.spawn_count;
		if (ring_position >= first_spawn) {
			if ((tick.flags & tfx_tick_generated) != 0) {
				tfx_generate_particle(push, dispatch, tick, ring_position - first_spawn, spawn_hook, hot, cold, noise);
			} else {
				const uint row = push.spawn_base + tick.upload_offset + ring_position - first_spawn;
				hot = spawn_buffers[push.spawn_buffer_index][row].hot;
				cold = spawn_buffers[push.spawn_buffer_index][row].cold;
				if (has_noise) {
					noise = spawn_buffers[push.spawn_buffer_index][row].noise;
				}
			}
			ingested = true;
			dead = false;
		} else if (!live) {
			hot = hot_buffers[push.hot_buffer_index][pool_slot];
			dead = (hot.flags & tfx_particle_remove) != 0;
			if (!dead) {
				cold = cold_buffers[push.cold_buffer_index][push.cold_base + pool_slot];
				if (has_noise) {
					noise = noise_buffers[push.noise_buffer_index][push.noise_base + pool_slot];
				}
			}
		}
		live = true;

		const bool last_tick = tick_index + 1 == end_tick;
		if (dead) {
			if (last_tick) {
				tfx_write_empty_instance(push, dispatch, hot.position, ring_position);
			}
			continue;
		}
		const float life = hot.age * cold.inv_max_age;
		//An expired or killed particle holds its slot until its batch expires, and is written with no size
		if ((hot.flags & tfx_particle_remove) == 0) {
			tfx_particle_context context;
			context.push = push;
			context.dispatch = dispatch;
			context.tick = tick;
			context.noise = noise;
			context.life = life;
			context.ring_position = ring_position;
			context.last_tick = last_tick;
			tfx_control_particle(context, hot, cold, velocity_graph, weight_graph, drag_variation, update);
			update.after_move(hot, cold, context);
		}
		if (last_tick) {
			uint2 location_rotation;
			const float3 position = tfx_particle_world_position(push, dispatch, tick, hot, cold, location_rotation);
			tfx_write_particle_instance(push, dispatch, tick, hot, cold, position, life, ring_position, location_rotation);
		}
		tfx_age_particle(hot, cold, tick);
	}

	if (live && !dead) {
		hot_buffers[push.hot_buffer_index][pool_slot] = hot;
		if (ingested) {
			cold_buffers[push.cold_buffer_index][push.cold_base + pool_slot] = cold;
			if (has_noise) {
				noise_buffers[push.noise_buffer_index][push.noise_base + pool_slot] = noise;
			}
		}
	}
}

//For a host with an update but no spawn hook
void tfx_simulate<TUpdate : ITfxParticleUpdate>(uint thread_index, tfx_gpu_particle_push push, TUpdate update) {
	tfx_simulate(thread_index, push, update, tfx_no_particle_spawn());
}

//Only the library's entry point reads these, so a host's entry point can declare its own push constants
[[vk::push_constant]] ConstantBuffer<tfx_gpu_particle_push> tfx_push_constants;

[shader("compute")]
[numthreads(tfx_group_size, 1, 1)]
void simulate_particles(uint3 thread_id : SV_DispatchThreadID) {
	tfx_simulate<tfx_no_particle_update>(thread_id.x, tfx_push_constants, tfx_no_particle_update());
}
