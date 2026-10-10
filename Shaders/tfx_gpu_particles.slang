//The reference GPU particle kernel. One dispatch per tfx_gpu_particle_dispatch_t from tfx_GetGPUParticleFrame, one thread per ring
//slot. Each thread runs every tick record of the fetch in registers, mirroring the CPU's control, instance and age passes in their
//order, and writes its instance on the last tick. Ring positions count from the tick's head, so a slot that is bumped and reused
//within one fetch just ingests its new occupant at the tick it was spawned.

import tfx_graphs;

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
static const uint tfx_dispatch_relative_scale_only = 1 << 15;
static const uint tfx_dispatch_location_rotation = 1 << 16;
static const uint tfx_dispatch_align_rotated_emission = 1 << 17;

//tfx_gpu_particle_flag_bits
static const uint tfx_particle_remove = 1 << 0;

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
	uint padding;
};

//Mirrors tfx_gpu_spawn_location_t
struct tfx_gpu_spawn_location {
	float3 position;
	uint local_id;						//slot | generation, tfxINVALID when the slot is empty
	uint2 packed_rotation;				//16 bit snorm: .x = X|Y, .y = Z|W
	uint2 padding;
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
	uint hot_buffer_index;
	uint cold_buffer_index;
	uint instance_buffer_index;
	uint spawn_hot_buffer_index;
	uint spawn_cold_buffer_index;
	uint dispatch_buffer_index;
	uint tick_buffer_index;
	uint graph_buffer_index;
	uint location_buffer_index;
	uint dispatch_index;
};

[[vk::binding(5, 0)]] RWStructuredBuffer<tfx_gpu_particle_hot> hot_buffers[];
[[vk::binding(5, 0)]] RWStructuredBuffer<tfx_gpu_particle_cold> cold_buffers[];
[[vk::binding(5, 0)]] RWStructuredBuffer<tfx_instance> instance_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_particle_hot> spawn_hot_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_particle_cold> spawn_cold_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_particle_dispatch> dispatch_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_particle_tick> tick_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_graph_data> graph_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_gpu_spawn_location> location_buffers[];
[[vk::push_constant]] ConstantBuffer<tfx_gpu_particle_push> push;

//tfx__wide_unpack10bit: x in the top 10 bits, z in the bottom
float3 tfx_unpack10bit(uint packed) {
	const float one_div_511 = 1.0 / 511.0;
	return float3(float((packed >> 20) & 0x3FFu), float((packed >> 10) & 0x3FFu), float(packed & 0x3FFu)) * one_div_511 - 1.0;
}

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

//tfx__wide_unpack16bit: each lane sign extended, x and y from .x, z and w from .y
float4 tfx_unpack16bit_quaternion(uint2 packed) {
	const int4 values = int4(int(packed.x << 16) >> 16, int(packed.x) >> 16, int(packed.y << 16) >> 16, int(packed.y) >> 16);
	return float4(values) * (1.0 / 32767.0);
}

//Port of tfx__white_unit: the signed reinterpret is what centres it on 0
float tfx_white_unit(uint seed, uint axis) {
	uint hash = seed ^ axis;
	hash ^= hash >> 15;
	hash *= 0x85ebca6bu;
	hash ^= hash >> 13;
	hash *= 0xc2b2ae35u;
	hash ^= hash >> 16;
	return (float(int(hash)) / 4294967295.0) * 2.0;
}

//The CPU's rotate by quaternion, tfx__wide_transform_quaternion_vec3
float3 tfx_rotate_vector(float3 v, float4 q) {
	return v + 2.0 * cross(cross(v, q.xyz) + q.w * v, q.xyz);
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

tfx_gpu_graph_data tfx_graph(uint graph_offset, uint graph) {
	return graph_buffers[push.graph_buffer_index][graph_offset + graph];
}

//tfx_apply_lookup_velocity, tfx_apply_lookup_weight, tfx_apply_velocity and tfx_apply_position, without forces
void tfx_control_particle(inout tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, tfx_gpu_particle_tick tick, float life, tfx_gpu_graph_data velocity_graph, tfx_gpu_graph_data weight_graph, bool drag_variation) {
	const float velocity = tfx_sample_graph(velocity_graph, life) * cold.base_velocity;
	const float weight = tfx_sample_graph(weight_graph, life) * cold.base_weight;
	float3 acceleration = tfx_unpack10bit(cold.velocity_normal) * velocity;
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

//tfx__block_particle_transform. A relative particle whose location has gone is removed here, before the write, so it's drawn with no
//size and isn't captured, which is what the CPU's location_alive mask does in the instance pass ahead of its age pass kill
float3 tfx_particle_world_position(tfx_gpu_particle_dispatch dispatch, tfx_gpu_particle_tick tick, inout tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold) {
	const uint flags = dispatch.flags;
	if ((flags & tfx_dispatch_relative_location) != 0) {
		const uint slot = cold.spawn_location & tfx_spawn_location_slot_mask;
		float3 location_position = float3(0.0, 0.0, 0.0);
		uint2 location_rotation = uint2(0u, 0x7FFF0000u);
		uint location_id = 0xFFFFFFFFu;
		if (slot < dispatch.location_count) {
			const tfx_gpu_spawn_location location = location_buffers[push.location_buffer_index][dispatch.location_offset + slot];
			location_position = location.position;
			location_rotation = location.packed_rotation;
			location_id = location.local_id;
		}
		if (location_id != cold.spawn_location) {
			hot.flags |= tfx_particle_remove;
		}
		float3 position = tfx_rotate_vector(hot.position + dispatch.handle, dispatch.emitter_rotation);
		//After the emitter rotation, so the emitter keeps its own angle and the location turns the whole thing into its frame
		if ((flags & tfx_dispatch_location_rotation) != 0) {
			position = tfx_rotate_vector(position, tfx_unpack16bit_quaternion(location_rotation));
		}
		return position * tick.overall_scale + location_position;
	}
	if ((flags & tfx_dispatch_relative_scale_only) != 0) {
		return hot.position * tick.overall_scale;
	}
	return hot.position;
}

//tfx__control_particle_age for a particle that isn't single
void tfx_age_particle(inout tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, tfx_gpu_particle_tick tick) {
	hot.age += tick.frame_length;
	if (hot.age * cold.inv_max_age >= 1.0) {
		hot.flags |= tfx_particle_remove;
	}
}

//tfx__block_particle_spin
uint2 tfx_particle_quaternion(tfx_gpu_particle_dispatch dispatch, tfx_gpu_particle_cold cold, float life) {
	const uint flags = dispatch.flags;
	if ((flags & tfx_dispatch_spin_pitch_and_yaw) == 0) {
		float roll = 0.0;
		if ((flags & tfx_dispatch_spin_roll) != 0) {
			roll = tfx_sample_graph(tfx_graph(dispatch.graph_offset, tfx_graph_roll_spin), life) * cold.base_roll_spin;
		}
		if ((flags & tfx_dispatch_relative_angle) != 0) {
			roll += dispatch.emitter_roll;
		}
		roll += asfloat(cold.rotation_offsets);
		return tfx_pack16bit_quaternion(tfx_euler_to_quaternion(0.0, 0.0, roll));
	}

	const float3 offsets = tfx_unpack10bit(cold.rotation_offsets) * tfx_two_pi;
	float pitch = 0.0;
	float yaw = 0.0;
	float roll = 0.0;
	if ((flags & tfx_dispatch_spin_3d) != 0) {
		pitch = tfx_sample_graph(tfx_graph(dispatch.graph_offset, tfx_graph_pitch_spin), life) * cold.base_pitch_spin;
		yaw = tfx_sample_graph(tfx_graph(dispatch.graph_offset, tfx_graph_yaw_spin), life) * cold.base_yaw_spin;
		roll = tfx_sample_graph(tfx_graph(dispatch.graph_offset, tfx_graph_roll_spin), life) * cold.base_roll_spin;
	}
	float4 q = tfx_euler_to_quaternion(pitch + offsets.x, yaw + offsets.y, roll + offsets.z);

	//conj(emitter) * particle, because the CPU's position transform rotates the opposite way to the vertex shader
	if ((flags & tfx_dispatch_relative_angle) != 0) {
		const float4 e = dispatch.emitter_rotation;
		q = float4(
			(e.w * q.x - e.x * q.w) + (e.z * q.y - e.y * q.z),
			(e.w * q.y - e.y * q.w) + (e.x * q.z - e.z * q.x),
			(e.w * q.z - e.z * q.w) + (e.y * q.x - e.x * q.y),
			(e.w * q.w + e.x * q.x) + (e.y * q.y + e.z * q.z));
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

//tfx__control_particle_instances for one particle, at ring position ring_position of the last tick
void tfx_write_particle_instance(tfx_gpu_particle_dispatch dispatch, tfx_gpu_particle_tick tick, tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, float3 position, float life, uint ring_position) {
	const uint flags = dispatch.flags;
	const uint graph_offset = dispatch.graph_offset;
	const bool alive = (hot.flags & tfx_particle_remove) == 0;
	tfx_instance instance;

	//----Transform
	const float stretch = tfx_sample_graph(tfx_graph(graph_offset, tfx_graph_stretch), life) * dispatch.global_stretch;
	instance.position = float4(position, stretch);
	instance.alignment = 0;
	if ((flags & tfx_dispatch_align_emission) != 0) {
		float3 alignment = tfx_unpack10bit(cold.velocity_normal);
		if ((flags & tfx_dispatch_align_rotated_emission) != 0) {
			alignment = tfx_rotate_vector(alignment, dispatch.emitter_rotation);
		}
		instance.alignment = tfx_pack8bit_xyz(alignment);
	} else if ((flags & tfx_dispatch_align_emitter) != 0) {
		instance.alignment = tfx_pack8bit_xyz(tfx_rotate_vector(float3(0.0, 1.0, 0.0), dispatch.emitter_rotation));
	}

	//----Spin
	instance.quaternion = tfx_particle_quaternion(dispatch, cold, life);

	//----Size, zero for particles that expired but still hold their slot
	instance.size = 0;
	if ((flags & tfx_dispatch_do_not_render) == 0 && alive) {
		const float size_life = min(life, 1.0);
		const float lookup_width = tfx_sample_graph(tfx_graph(graph_offset, tfx_graph_width), size_life);
		float scale_x = cold.base_size_x * lookup_width;
		float scale_y;
		if ((flags & tfx_dispatch_lifetime_uniform_size) != 0) {
			scale_y = lookup_width * cold.base_size_y;
			if ((flags & tfx_dispatch_base_uniform_size) != 0) {
				scale_y = min(scale_x, scale_y);
			}
		} else {
			scale_y = tfx_sample_graph(tfx_graph(graph_offset, tfx_graph_height), size_life) * cold.base_size_y;
		}
		instance.size = tfx_pack16bit_2sscaled(scale_x * tick.overall_scale, scale_y * tick.overall_scale, tfx_size_packed_scale);
	}

	//----Color
	const float intensity = dispatch.global_intensity * tfx_sample_graph(tfx_graph(graph_offset, tfx_graph_intensity), life) * cold.intensity_factor;
	const float curved_alpha = tfx_sample_graph(tfx_graph(graph_offset, tfx_graph_curved_alpha), life);
	const float alpha_sharpness = tfx_sample_graph(tfx_graph(graph_offset, tfx_graph_alpha_sharpness), life);
	const float gradient_map = (flags & tfx_dispatch_random_color) != 0 ? cold.random_color : tfx_sample_graph(tfx_graph(graph_offset, tfx_graph_gradient_mapper), life);
	instance.intensity_gradient_map = tfx_pack16bit_2sscaled(intensity, gradient_map, tfx_color_packed_scale);
	instance.curved_alpha_life = tfx_pack8bitunorm_xyz(float3(curved_alpha, alpha_sharpness, life));

	//----Image frame and interpolation. Only particles that were in the ring at the previous write have somewhere to come from, and
	//an expired particle isn't captured, so it doesn't shrink towards zero over its last frame
	const uint instance_index = dispatch.instance_start + ring_position;
	const bool capture = (flags & tfx_dispatch_capture) != 0 && alive && ring_position + dispatch.spawned_since_write < tick.count;
	instance.indexes = tfx_particle_image_frame(dispatch, cold, hot, tick) | (capture ? tfx_capture_bit : 0);
	instance.captured_index = capture ? dispatch.previous_instance_start + ring_position + dispatch.bumped_since_write : instance_index;

	instance_buffers[push.instance_buffer_index][instance_index] = instance;
}

[shader("compute")]
[numthreads(tfx_group_size, 1, 1)]
void simulate_particles(uint3 thread_id : SV_DispatchThreadID) {
	const tfx_gpu_particle_dispatch dispatch = dispatch_buffers[push.dispatch_buffer_index][push.dispatch_index];
	if (thread_id.x >= dispatch.slot_count) {
		return;
	}
	const uint ring_slot = (dispatch.slot_start + thread_id.x) % dispatch.capacity;
	const uint pool_slot = dispatch.base + ring_slot;
	const tfx_gpu_graph_data velocity_graph = tfx_graph(dispatch.graph_offset, tfx_graph_velocity);
	const tfx_gpu_graph_data weight_graph = tfx_graph(dispatch.graph_offset, tfx_graph_weight);
	const bool drag_variation = (dispatch.flags & tfx_dispatch_drag_variation) != 0;

	tfx_gpu_particle_hot hot;
	tfx_gpu_particle_cold cold;
	bool live = false;
	bool ingested = false;
	const uint end_tick = dispatch.first_tick + dispatch.tick_count;
	for (uint tick_index = dispatch.first_tick; tick_index != end_tick; ++tick_index) {
		const tfx_gpu_particle_tick tick = tick_buffers[push.tick_buffer_index][tick_index];
		const uint ring_position = (ring_slot + dispatch.capacity - tick.head) % dispatch.capacity;
		if (ring_position >= tick.count) {
			live = false;
			continue;
		}
		//Older particles only ever move towards the head, so the slots past count - spawn_count can only hold this tick's spawns
		const uint first_spawn = tick.count - tick.spawn_count;
		if (ring_position >= first_spawn) {
			const uint row = tick.upload_offset + ring_position - first_spawn;
			hot = spawn_hot_buffers[push.spawn_hot_buffer_index][row];
			cold = spawn_cold_buffers[push.spawn_cold_buffer_index][row];
			ingested = true;
		} else if (!live) {
			hot = hot_buffers[push.hot_buffer_index][pool_slot];
			cold = cold_buffers[push.cold_buffer_index][pool_slot];
		}
		live = true;

		const float life = hot.age * cold.inv_max_age;
		tfx_control_particle(hot, cold, tick, life, velocity_graph, weight_graph, drag_variation);
		if (tick_index + 1 == end_tick) {
			const float3 position = tfx_particle_world_position(dispatch, tick, hot, cold);
			tfx_write_particle_instance(dispatch, tick, hot, cold, position, life, ring_position);
		}
		tfx_age_particle(hot, cold, tick);
	}

	if (live) {
		hot_buffers[push.hot_buffer_index][pool_slot] = hot;
		if (ingested) {
			cold_buffers[push.cold_buffer_index][pool_slot] = cold;
		}
	}
}
