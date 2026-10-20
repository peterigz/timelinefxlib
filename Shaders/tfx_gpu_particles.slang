/*
	TimelineFX GPU particle kernel

	The reference compute kernel for Run on GPU emitters, one file a host can copy into its renderer and compile with Slang for
	Vulkan, D3D12, D3D11, Metal, OpenGL or WebGPU. Search for a keyword below to jump to its section (in Vim, * on the keyword):

	[About_the_kernel]                  What the kernel does, one dispatch per ring and one thread per ring slot

	--Mirrors of the library, change these together with the C functions they name
	[Shared_hashes]                     Two pi, the seed hash and its unit draws
	[Buffer_rows]                       Row structs, byte for byte mirrors of the C structs in timelinefx.h
	[Graph_lookups]                     Graph sampling, the twin of tfx__sample_wide_graph
	[Noise_fields]                      White, simplex, curl and value noise, ports of the CPU's wide noise
	[Spawn_generation]                  One particle from a tick's spawn block, the twin of tfx__generate_gpu_spawn
		[Spawn_constants]               Block flags, hash axes, emission types
		[Spawn_rows]                    The block's five rows, runs, multipliers and the generated spawn
		[Spawn_math]                    Packing, quaternions, the fast inverse square root and cone draws
		[Grid_and_sequences]            Uniform distribution sequences and grid marching
		[Stepped_angles]                Stepped emission and roll divisions
		[Spawn_shapes]                  Area, line, ellipse, disc and cylinder positions and directions
		[Generate_spawn]                The whole particle

	--Binding
	[Buffer_access]                     Every buffer read and write, for each binding model
		[Binding_macros]                TFX_FIXED_BINDINGS and the per target cases
		[Push_constants]                The library's push constants, and their uniform buffer twin for WGSL and OpenGL
		[Bindless_arrays]               Vulkan and D3D12
		[Metal_buffer_table]            Metal's argument buffer of gpuAddress values
		[Fixed_raw_buffers]             WGSL, and every target compiled with TFX_FIXED_BINDINGS
		[Buffer_accessors]              The tfx_load_* and tfx_store_* functions

	--The kernel
	[Kernel_constants]                  Graph slots, dispatch, particle, tick and force flags
	[Instance_packing]                  Pack helpers and the spin quaternion
	[Emitter_forces]                    Wind, attract, vortex, shockwave and noise from one force row
	[Ring_space]                        Ring space to world, the particle context and the location frame
	[Force_rows]                        Location timed forces and tfx_force_medium over the tick's rows
	[Host_hooks]                        ITfxParticleUpdate, ITfxParticleSpawn and ITfxParticleWrite
	[Particle_control]                  Velocity, forces, weight, drag, integration and age
	[Instance_write]                    Spin, image frame and the instance write
	[Spawn_ingest]                      A tick's spawns, generated from its block or loaded from rows
	[Kernel_body]                       tfx_simulate and its overloads
	[Entry_point]                       Picking the push constants per target, and simulate_particles
*/

//--- About_the_kernel

//The reference GPU particle kernel. One dispatch per tfx_gpu_particle_dispatch_t from tfx_GetGPUParticleFrame, one thread per ring
//slot. Each thread runs every tick record of the fetch in registers, mirroring the CPU's control, instance and age passes in their
//order, and writes its instance on the last tick. Ring positions count from the tick's head, so a slot that is bumped and reused
//within one fetch just ingests its new occupant at the tick it was spawned: from a spawn row, or generated from the tick's spawn block
//when the tick is flagged generated (see the spawn generation section).

//--- Shared_hashes

//The constants and hashes the rest of the file shares, each the twin of the CPU function it names

//tfxPI2
static const float tfx_two_pi = 6.283185307;

//tfx__seedgen_u32, a murmur3 finalizer. tfx__wide_seedgen_base(base, h) is tfx_seedgen(base ^ h)
uint tfx_seedgen(uint hash) {
	hash ^= hash >> 15;
	hash *= 0x85ebca6bu;
	hash ^= hash >> 13;
	hash *= 0xc2b2ae35u;
	hash ^= hash >> 16;
	return hash;
}

//tfx__white_unit and tfx__wide_white_unit: the signed reinterpret is what centres it on 0
float tfx_white_unit(uint seed, uint axis) {
	return (float(int(tfx_seedgen(seed ^ axis))) / 4294967295.0) * 2.0;
}

//tfx__hash_unit, exact on both sides because 24 bits convert to a float without rounding
float tfx_hash_unit(uint seed, uint axis) {
	return float(tfx_seedgen(seed ^ axis) >> 8) * (1.0 / 16777216.0);
}

//--- Buffer_rows

//The GPU particle buffers' rows, each the byte for byte mirror of its C struct in timelinefx.h. Check the SPIR-V member offsets after changing either side

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
	uint user_value;					//The spawn location's user data unless after_spawn changed it, 0 without a location
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
	float location_age_offset;			//Milliseconds the later ticks of the fetch added to an unpaused location's age
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
	float location_delay;				//The emitter's spawn delay, waited out from each location's age
	float location_loop_length;			//0 if the emitter doesn't loop
};

//Mirrors tfx_gpu_spawn_location_t
struct tfx_gpu_spawn_location {
	float3 position;
	uint local_id;						//slot | generation, tfxINVALID when the slot is empty
	uint2 packed_rotation;				//16 bit snorm: .x = X|Y, .y = Z|W
	float age;							//Milliseconds, after the last tick of the fetch
	uint paused;
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
	float turn_offset_0;
	float3 harmonic_cos;
	float turn_offset_1;
	float3 harmonic_sin;
	float turn_offset_2;
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

//--- Graph_lookups

//Graph lookups for the GPU kernels. Mirrors tfx__sample_wide_graph: one segment, easing, then bezier or linear, then the oscillator

static const int tfx_graph_use_bezier_sampling = 1;
static const int tfx_graph_enable_oscillator = 2;

static const int tfx_easing_constant = 0;
static const int tfx_easing_in = 4;
static const int tfx_easing_out = 5;
static const int tfx_easing_in_out = 6;
static const int tfx_easing_smoothstep = 17;
static const int tfx_easing_out_in = 18;

//Mirrors tfx_gpu_graph_data_t, one per graph from tfx_GetGPUGraphLookupsBuffer
struct tfx_gpu_graph_data {
	float4 node_data;					//from, curve1, curve2, to
	float4 oscillator;					//frequency, amplitude, offset_x, offset_y
	int easing_type;
	int flags;
	int padding[2];
};

float tfx_apply_easing(int easing_type, float t) {
	switch (easing_type) {
		case tfx_easing_constant: return 0.0;
		case tfx_easing_in: return t * t * t;
		case tfx_easing_out: {
			float inverse = 1.0 - t;
			return 1.0 - inverse * inverse * inverse;
		}
		case tfx_easing_in_out: {
			float t3 = -2.0 * t + 2.0;
			t3 = t3 * t3 * t3;
			return t <= 0.5 ? 4.0 * t * t * t : 1.0 - t3 * 0.5;
		}
		case tfx_easing_smoothstep: return t * t * (3.0 - 2.0 * t);
		case tfx_easing_out_in: return t <= 0.5 ? t * (2.0 - 2.0 * t) : 2.0 * (t - 0.5) * (t - 0.5) + 0.5;
		default: return t;
	}
}

//Constant easing samples t = 0, which is why the CPU skips the bezier and the oscillator for it, and so does this
float tfx_sample_graph(tfx_gpu_graph_data graph, float t) {
	const float time = tfx_apply_easing(graph.easing_type, t);
	const bool constant = graph.easing_type == tfx_easing_constant;
	float value;
	if ((graph.flags & tfx_graph_use_bezier_sampling) != 0 && !constant) {
		const float inverse = 1.0 - time;
		const float weight1 = inverse * inverse * inverse;
		const float weight2 = 3.0 * inverse * inverse * time;
		const float weight3 = 3.0 * inverse * time * time;
		const float weight4 = time * time * time;
		value = weight1 * graph.node_data.x + weight2 * graph.node_data.y + weight3 * graph.node_data.z + weight4 * graph.node_data.w;
	} else {
		value = graph.node_data.x + (graph.node_data.w - graph.node_data.x) * time;
	}
	if ((graph.flags & tfx_graph_enable_oscillator) != 0 && !constant && graph.oscillator.x != 0.0 && graph.oscillator.y != 0.0) {
		const float oscillator = 0.5 + sin((time + graph.oscillator.z) * graph.oscillator.x * tfx_two_pi) * graph.oscillator.y;
		value *= oscillator + graph.oscillator.w;
	}
	return value;
}

//--- Noise_fields

//Noise fields for the GPU particle kernel, scalar ports of the CPU's wide noise in timelinefx_internal.h. Every hash, prime and
//operation order follows the CPU so the two only differ by rounding.

static const float tfx_noise_epsilon = 1e-12;
static const float tfx_noise_offset = 100.0;			//tfxWIDENOISEOFFSET
static const float tfx_simplex_kernel_radius_squared = 0.5;	//tfxSIMPLEX_KERNEL_RADIUS_SQUARED
static const float tfx_curl_scale = 30.857;

//---- Value noise, tfx__wide_value_noise_3d

float3 tfx_value_noise_corner(uint lattice_x, uint lattice_y, uint lattice_z) {
	const uint h = tfx_seedgen((lattice_x + lattice_y) ^ lattice_z);
	return float3(float((h >> 20) & 0x3FFu), float((h >> 10) & 0x3FFu), float(h & 0x3FFu)) * (1.0 / 511.0) - 1.0;
}

float3 tfx_value_noise_lerp(float3 a, float3 b, float t) {
	return t * (b - a) + a;
}

float3 tfx_value_noise_3d(float3 position) {
	const float3 floor_position = floor(position);
	float3 fade = position - floor_position;
	fade = fade * fade * (3.0 - 2.0 * fade);

	const uint3 primes = uint3(1619u, 31337u, 6971u);
	const uint3 lattice0 = uint3(int3(floor_position)) * primes;
	const uint3 lattice1 = lattice0 + primes;

	float3 edge = tfx_value_noise_lerp(tfx_value_noise_corner(lattice0.x, lattice0.y, lattice0.z), tfx_value_noise_corner(lattice1.x, lattice0.y, lattice0.z), fade.x);
	const float3 slice = tfx_value_noise_lerp(edge, tfx_value_noise_lerp(tfx_value_noise_corner(lattice0.x, lattice1.y, lattice0.z), tfx_value_noise_corner(lattice1.x, lattice1.y, lattice0.z), fade.x), fade.y);
	edge = tfx_value_noise_lerp(tfx_value_noise_corner(lattice0.x, lattice0.y, lattice1.z), tfx_value_noise_corner(lattice1.x, lattice0.y, lattice1.z), fade.x);
	return tfx_value_noise_lerp(slice, tfx_value_noise_lerp(edge, tfx_value_noise_lerp(tfx_value_noise_corner(lattice0.x, lattice1.y, lattice1.z), tfx_value_noise_corner(lattice1.x, lattice1.y, lattice1.z), fade.x), fade.y), fade.z);
}

//---- Simplex noise, tfx__wide_simplex_cell and friends

//Named corners rather than arrays, so nothing is left for a driver to index dynamically or spill to scratch
struct tfx_simplex_cell {
	float3 offset0, offset1, offset2, offset3;
	uint hash0, hash1, hash2, hash3;
};

tfx_simplex_cell tfx_simplex_cell_at(float3 position) {
	const float skew = (position.x + (position.y + position.z)) * (1.0 / 3.0);
	const float3 cell = floor(position + skew);
	const float unskew = ((cell.x + cell.y) + cell.z) * (1.0 / 6.0);
	const float3 offset0 = position - (cell - unskew);

	const bool x_ge_y = offset0.x >= offset0.y;
	const bool y_ge_z = offset0.y >= offset0.z;
	const bool x_ge_z = offset0.x >= offset0.z;
	const uint3 first_step = uint3(x_ge_y && x_ge_z ? 1u : 0u, !x_ge_y && y_ge_z ? 1u : 0u, !(x_ge_z || y_ge_z) ? 1u : 0u);
	const uint3 second_step = uint3(x_ge_y || x_ge_z ? 1u : 0u, !x_ge_y || y_ge_z ? 1u : 0u, !(y_ge_z && x_ge_z) ? 1u : 0u);

	tfx_simplex_cell result;
	result.offset0 = offset0;
	result.offset1 = (offset0 - float3(first_step)) + (1.0 / 6.0);
	result.offset2 = (offset0 - float3(second_step)) + (1.0 / 3.0);
	result.offset3 = offset0 - 0.5;

	//floor() is already integral, so the truncating convert is exact for negative cells too
	const uint3 primes = uint3(1619u, 31337u, 6971u);
	const uint3 lattice = uint3(int3(cell)) * primes;
	const uint3 first_corner = lattice + first_step * primes;
	const uint3 second_corner = lattice + second_step * primes;
	const uint3 last_corner = lattice + primes;
	result.hash0 = tfx_seedgen((lattice.x + lattice.y) ^ lattice.z);
	result.hash1 = tfx_seedgen((first_corner.x + first_corner.y) ^ first_corner.z);
	result.hash2 = tfx_seedgen((second_corner.x + second_corner.y) ^ second_corner.z);
	result.hash3 = tfx_seedgen((last_corner.x + last_corner.y) ^ last_corner.z);
	return result;
}

float tfx_simplex_falloff(float3 offset) {
	return max(((tfx_simplex_kernel_radius_squared - offset.x * offset.x) - offset.y * offset.y) - offset.z * offset.z, 0.0);
}

//Perlin's 12 edge gradient, picked by the low 4 hash bits as tfx__wide_simplex_corner_edge does
float3 tfx_simplex_corner_edge(uint hash) {
	const bool u_is_x = (hash & 8u) == 0u;
	const bool v_is_y = (hash & 12u) == 0u;
	const bool v_is_x = (hash & 13u) == 12u;
	const bool v_is_z = !(v_is_y || v_is_x);
	const float u_signed_one = (hash & 1u) != 0u ? -1.0 : 1.0;
	const float v_signed_one = (hash & 2u) != 0u ? -1.0 : 1.0;
	return float3(u_is_x ? u_signed_one : (v_is_x ? v_signed_one : 0.0),
				  !u_is_x ? u_signed_one : (v_is_y ? v_signed_one : 0.0),
				  v_is_z ? v_signed_one : 0.0);
}

float tfx_simplex_corner(uint hash, float3 offset) {
	const float3 edge = tfx_simplex_corner_edge(hash);
	float falloff = tfx_simplex_falloff(offset);
	falloff *= falloff;
	falloff *= falloff;
	return falloff * (edge.x * offset.x + (edge.y * offset.y + edge.z * offset.z));
}

float tfx_simplex_noise_3d(float3 position) {
	const tfx_simplex_cell cell = tfx_simplex_cell_at(position);
	float sum = tfx_simplex_corner(cell.hash0, cell.offset0);
	sum += tfx_simplex_corner(cell.hash1, cell.offset1);
	sum += tfx_simplex_corner(cell.hash2, cell.offset2);
	sum += tfx_simplex_corner(cell.hash3, cell.offset3);
	return sum * 32.0;
}

//tfx__wide_simplex_noise_3d_gradient: without the x32, which tfx_curl_scale carries
//d/dp of falloff^4 * dot(edge, offset), which is falloff^4 * edge - 8 * falloff^3 * dot * offset
float3 tfx_simplex_corner_gradient(uint hash, float3 offset) {
	const float3 edge = tfx_simplex_corner_edge(hash);
	const float falloff = tfx_simplex_falloff(offset);
	const float falloff_squared = falloff * falloff;
	const float falloff_fourth = falloff_squared * falloff_squared;
	const float dot_edge = edge.x * offset.x + (edge.y * offset.y + edge.z * offset.z);
	const float offset_scale = (-8.0 * (falloff_squared * falloff)) * dot_edge;
	return falloff_fourth * edge + offset_scale * offset;
}

float3 tfx_simplex_noise_3d_gradient(float3 position) {
	const tfx_simplex_cell cell = tfx_simplex_cell_at(position);
	float3 gradient = tfx_simplex_corner_gradient(cell.hash0, cell.offset0);
	gradient += tfx_simplex_corner_gradient(cell.hash1, cell.offset1);
	gradient += tfx_simplex_corner_gradient(cell.hash2, cell.offset2);
	gradient += tfx_simplex_corner_gradient(cell.hash3, cell.offset3);
	return gradient;
}

//---- The four algorithms as a unit direction (curl is left unnormalised, as on the CPU)

float3 tfx_noise_normalise(float3 noise) {
	return noise * rsqrt(max(dot(noise, noise), tfx_noise_epsilon));
}

//tfx__wide_apply_noise_simplex_force: three decorrelated samples of one field
float3 tfx_simplex_noise_direction(float3 sample_position, float life) {
	const float shift = tfx_noise_offset + life;
	return tfx_noise_normalise(float3(tfx_simplex_noise_3d(sample_position),
									  tfx_simplex_noise_3d(float3(sample_position.x, sample_position.y + shift, sample_position.z)),
									  tfx_simplex_noise_3d(float3(sample_position.x, sample_position.y, sample_position.z + shift))));
}

//tfx__wide_apply_noise_curl_force: Bridson curl from the analytic gradients of three potentials
float3 tfx_curl_noise_direction(float3 sample_position, float life) {
	const float psi2_shift = tfx_noise_offset + life;
	const float psi3_shift = (tfx_noise_offset + tfx_noise_offset) + life;
	const float3 psi1 = tfx_simplex_noise_3d_gradient(sample_position);
	const float3 psi2 = tfx_simplex_noise_3d_gradient(sample_position + psi2_shift);
	const float3 psi3 = tfx_simplex_noise_3d_gradient(sample_position + psi3_shift);
	return float3(psi3.y - psi2.z, psi1.z - psi3.x, psi2.x - psi1.y) * tfx_curl_scale;
}

float3 tfx_value_noise_direction(float3 sample_position) {
	return tfx_noise_normalise(tfx_value_noise_3d(sample_position));
}

//tfx__wide_apply_white_noise_force's direction, already scaled by the inverse length the CPU folds into its amplitude
float3 tfx_white_noise_direction(uint uid, float age, float noise_resolution) {
	const float inv_period_scale = 1.0 / (250.0 * 300.0);
	const uint axis_x = 0x1b56c4f9u;
	const uint axis_y = 0x68e31da4u;
	const uint axis_z = 0xb5297a4du;
	const uint axis_phase = 0x27d4eb2fu;

	const float phase_offset = (tfx_white_unit(uid, axis_phase) * 0.5 + 0.5) * 256.0;
	const float phase = phase_offset + age * (noise_resolution * inv_period_scale);
	const uint bucket = uint(int(phase));
	float fraction = phase - float(int(bucket));
	fraction = fraction * fraction * (3.0 - 2.0 * fraction);

	const uint direction_seed = tfx_seedgen(bucket ^ uid);
	const uint next_direction_seed = tfx_seedgen((bucket + 1u) ^ uid);
	const float one_minus_fraction = 1.0 - fraction;
	const float3 noise = float3(tfx_white_unit(direction_seed, axis_x) * one_minus_fraction + tfx_white_unit(next_direction_seed, axis_x) * fraction,
								tfx_white_unit(direction_seed, axis_y) * one_minus_fraction + tfx_white_unit(next_direction_seed, axis_y) * fraction,
								tfx_white_unit(direction_seed, axis_z) * one_minus_fraction + tfx_white_unit(next_direction_seed, axis_z) * fraction);
	return noise * rsqrt(max(dot(noise, noise), tfx_noise_epsilon));
}

//tfx__wide_split_medium_velocity: steer and surge against the particle's heading by the two noise gains
float3 tfx_split_medium_velocity(float3 travel, float3 medium, float swerve_gain_minus_one, float speed_gain_minus_swerve_gain) {
	const float length_squared = dot(travel, travel);
	const float inverse_length = rsqrt(max(length_squared, 1e-8));
	travel *= inverse_length;
	const float movement_fade = length_squared * (inverse_length * inverse_length);
	const float swerve_gain = movement_fade * swerve_gain_minus_one + 1.0;
	const float along_scale = dot(medium, travel) * (movement_fade * speed_gain_minus_swerve_gain);
	return along_scale * travel + swerve_gain * medium;
}

//--- Spawn_generation

//Spawn generation for GPU particles, the mirror of tfx__generate_gpu_spawn: one particle from a tick's spawn block, its location's run and
//multipliers. Everything here is a pure function of its arguments, so a test entry compiled with slangc -target cpp can be checked
//against the library. Change the two together.

//--- Spawn_constants

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
//tfxGOLDEN_ANGLE
static const float tfx_spawn_golden_angle = 2.39996323;

//--- Spawn_rows

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
	uint user_value;					//The location's user data by default, kept in the particle's cold row for the host's hooks
	uint user_flags;					//The particle's starting user flag bits, tfx_particle_user_flags_mask of its hot flags
};

//--- Spawn_math

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

//--- Grid_and_sequences

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

//--- Stepped_angles

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
	const bool full_circle = range >= tfx_two_pi - 0.001;
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
		const float jitter_cone = block_e.emission_step_jitter * 0.5 * sqrt(tfx_two_pi * cap_height_step);
		if (jitter_cone != 0.0) {
			const float z = 1.0 - (1.0 - cos(jitter_cone)) * tfx_hash_unit(uid, tfx_spawn_axis_cone_z);
			const float phi = tfx_two_pi * tfx_hash_unit(uid, tfx_spawn_axis_cone_phi);
			stepped = tfx_vector_in_cone(tfx_normalize_fast(stepped), z, phi);
		}
		return stepped;
	}
	const bool ring = (step_flags & tfx_step_ring) != 0;
	const float2 span = tfx_angle_step_span(ring ? tfx_two_pi : block_e.emission_range, 0.0, divisions);
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
	const float centre = range >= tfx_two_pi - 0.001 ? 0.0 : range * 0.5;
	const float2 span = tfx_angle_step_span(range, centre, divisions);
	float angle = span.x + span.y * float(division);
	const float half_jitter = span.y * block_e.roll_step_jitter * 0.5;
	if (half_jitter > 0.0) {
		angle += (2.0 * tfx_hash_unit(uid, tfx_spawn_axis_roll_jitter) - 1.0) * half_jitter;
	}
	return angle;
}

//--- Spawn_shapes

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
		const float theta = unit_u * tfx_two_pi;
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
			const float angle = unit_v * tfx_two_pi;
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

//--- Generate_spawn

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
		const float phi = tfx_two_pi * tfx_hash_unit(uid, tfx_spawn_axis_splatter_phi);
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
		const float phi = tfx_two_pi * tfx_hash_unit(uid, tfx_spawn_axis_cone_phi);
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
		spawn.rotation_offsets = tfx_pack10bit_unsigned(offsets * (1.0 / tfx_two_pi));
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
	spawn.spawn_location = run.local_id;
	spawn.noise_resolution = (flags & tfx_spawn_noise_resolution) != 0 ? block_a.noise_resolution : 0.0;
	spawn.user_value = run.user_data;
	spawn.user_flags = 0;
	return spawn;
}

//--- Buffer_access

//Every buffer access the kernel makes, so a port to another binding model replaces only this section
//Rows are relative to their own array and each accessor adds that array's base from the push constants, as tfx_GetGPUParticleFrame lays them out
//Bindings are shared on purpose here: the bindless arrays alias each other, and a host entry point's push constants take the library's binding
#pragma warning(disable: 39001)

//--- Binding_macros

//-D TFX_FIXED_BINDINGS gives every target the four fixed raw buffers WGSL always uses, for D3D11, OpenGL and anything else without descriptor arrays
//HLSL keeps tfx_push_constants under TFX_FIXED_BINDINGS: it's a plain cbuffer at b0, which is what D3D11 wants
#if defined(TFX_FIXED_BINDINGS)
#define TFX_FIXED_CASES case wgsl: case hlsl: case glsl: case spirv:
#define TFX_FIXED_REGISTER(slot) : register(slot)
#define TFX_BINDLESS_ACCESS
#define TFX_UNIFORM_PUSH_CASES case wgsl: case glsl:
#else
#define TFX_FIXED_CASES case wgsl:
#define TFX_FIXED_REGISTER(slot)
//GLSL needs the extension to index runtime descriptor arrays, and Slang only emits it from a function, never from the declarations
#define TFX_BINDLESS_ACCESS __glsl_extension(GL_EXT_nonuniform_qualifier)
#define TFX_UNIFORM_PUSH_CASES case wgsl:
#endif

//--- Push_constants

//Only the library's entry point reads these, so a host's entry point can declare its own push constants
[[vk::push_constant]] ConstantBuffer<tfx_gpu_particle_push> tfx_push_constants;
//WebGPU and OpenGL have no push constants, so there they're a uniform buffer moved along per dispatch. The register keeps it out of Metal's buffer slots
[[vk::binding(4, 0)]] ConstantBuffer<tfx_gpu_particle_push> tfx_uniform_push_constants : register(b15);

//--- Bindless_arrays

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

//--- Metal_buffer_table

//Metal has no descriptor arrays: there each push constant buffer index is a slot in this argument buffer of gpuAddress values
#ifndef TFX_METAL_BUFFER_COUNT
#define TFX_METAL_BUFFER_COUNT 16
#endif
struct tfx_metal_buffer_table {
	RWByteAddressBuffer buffers[TFX_METAL_BUFFER_COUNT];
};
ParameterBlock<tfx_metal_buffer_table> tfx_metal_buffers;

//--- Fixed_raw_buffers

//The four buffers of the one-buffer setup, for WGSL (no descriptor arrays, 8 storage buffers) and TFX_FIXED_BINDINGS. These ignore the indexes
[[vk::binding(0, 0)]] RWByteAddressBuffer tfx_fixed_pool TFX_FIXED_REGISTER(u0);
[[vk::binding(1, 0)]] RWByteAddressBuffer tfx_fixed_instances TFX_FIXED_REGISTER(u1);
[[vk::binding(2, 0)]] ByteAddressBuffer tfx_fixed_graphs TFX_FIXED_REGISTER(t0);
//The frame's upload: dispatches at its start, then the tick, location and spawn arrays at their bases
[[vk::binding(3, 0)]] ByteAddressBuffer tfx_fixed_upload TFX_FIXED_REGISTER(t1);

//The C sizes of the rows the raw paths read, which pack members as C does, whose Load/Store pack members as C does where a Metal float3 would take 16 bytes
static const uint tfx_hot_row_bytes = 32;
static const uint tfx_cold_row_bytes = 64;
static const uint tfx_noise_row_bytes = 16;
static const uint tfx_instance_row_bytes = 48;
static const uint tfx_spawn_row_bytes = 112;
static const uint tfx_spawn_cold_offset = 32;
static const uint tfx_spawn_noise_offset = 96;
static const uint tfx_dispatch_row_bytes = 128;
static const uint tfx_location_row_bytes = 32;
static const uint tfx_tick_row_bytes = 64;
static const uint tfx_graph_row_bytes = 48;
//offsetof(tfx_gpu_spawn_run_t, first_spawn)
static const uint tfx_spawn_run_first_spawn_offset = 4;

//--- Buffer_accessors

//---- The pool, by pool slot

TFX_BINDLESS_ACCESS tfx_gpu_particle_hot tfx_load_hot(tfx_gpu_particle_push push, uint slot) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.hot_buffer_index].Load<tfx_gpu_particle_hot>(slot * tfx_hot_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_pool.Load<tfx_gpu_particle_hot>(slot * tfx_hot_row_bytes);
	default: return hot_buffers[push.hot_buffer_index][slot];
	}
}

TFX_BINDLESS_ACCESS void tfx_store_hot(tfx_gpu_particle_push push, uint slot, tfx_gpu_particle_hot hot) {
	__target_switch {
	case metal: tfx_metal_buffers.buffers[push.hot_buffer_index].Store<tfx_gpu_particle_hot>(slot * tfx_hot_row_bytes, hot);
	TFX_FIXED_CASES tfx_fixed_pool.Store<tfx_gpu_particle_hot>(slot * tfx_hot_row_bytes, hot);
	default: hot_buffers[push.hot_buffer_index][slot] = hot;
	}
}

TFX_BINDLESS_ACCESS tfx_gpu_particle_cold tfx_load_cold(tfx_gpu_particle_push push, uint slot) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.cold_buffer_index].Load<tfx_gpu_particle_cold>((push.cold_base + slot) * tfx_cold_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_pool.Load<tfx_gpu_particle_cold>((push.cold_base + slot) * tfx_cold_row_bytes);
	default: return cold_buffers[push.cold_buffer_index][push.cold_base + slot];
	}
}

TFX_BINDLESS_ACCESS void tfx_store_cold(tfx_gpu_particle_push push, uint slot, tfx_gpu_particle_cold cold) {
	__target_switch {
	case metal: tfx_metal_buffers.buffers[push.cold_buffer_index].Store<tfx_gpu_particle_cold>((push.cold_base + slot) * tfx_cold_row_bytes, cold);
	TFX_FIXED_CASES tfx_fixed_pool.Store<tfx_gpu_particle_cold>((push.cold_base + slot) * tfx_cold_row_bytes, cold);
	default: cold_buffers[push.cold_buffer_index][push.cold_base + slot] = cold;
	}
}

TFX_BINDLESS_ACCESS tfx_gpu_particle_noise tfx_load_noise(tfx_gpu_particle_push push, uint slot) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.noise_buffer_index].Load<tfx_gpu_particle_noise>((push.noise_base + slot) * tfx_noise_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_pool.Load<tfx_gpu_particle_noise>((push.noise_base + slot) * tfx_noise_row_bytes);
	default: return noise_buffers[push.noise_buffer_index][push.noise_base + slot];
	}
}

TFX_BINDLESS_ACCESS void tfx_store_noise(tfx_gpu_particle_push push, uint slot, tfx_gpu_particle_noise noise) {
	__target_switch {
	case metal: tfx_metal_buffers.buffers[push.noise_buffer_index].Store<tfx_gpu_particle_noise>((push.noise_base + slot) * tfx_noise_row_bytes, noise);
	TFX_FIXED_CASES tfx_fixed_pool.Store<tfx_gpu_particle_noise>((push.noise_base + slot) * tfx_noise_row_bytes, noise);
	default: noise_buffers[push.noise_buffer_index][push.noise_base + slot] = noise;
	}
}

//---- The instance buffer, shared with the CPU instances

TFX_BINDLESS_ACCESS void tfx_store_instance(tfx_gpu_particle_push push, uint index, tfx_instance instance) {
	__target_switch {
	case metal: tfx_metal_buffers.buffers[push.instance_buffer_index].Store<tfx_instance>(index * tfx_instance_row_bytes, instance);
	TFX_FIXED_CASES tfx_fixed_instances.Store<tfx_instance>(index * tfx_instance_row_bytes, instance);
	default: instance_buffers[push.instance_buffer_index][index] = instance;
	}
}

//---- The frame's arrays from tfx_GetGPUParticleFrame

TFX_BINDLESS_ACCESS tfx_gpu_particle_dispatch tfx_load_dispatch(tfx_gpu_particle_push push, uint index) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.dispatch_buffer_index].Load<tfx_gpu_particle_dispatch>(index * tfx_dispatch_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_particle_dispatch>(index * tfx_dispatch_row_bytes);
	default: return dispatch_buffers[push.dispatch_buffer_index][index];
	}
}

TFX_BINDLESS_ACCESS tfx_gpu_particle_tick tfx_load_tick(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.tick_buffer_index].Load<tfx_gpu_particle_tick>((push.tick_base + row) * tfx_tick_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_particle_tick>((push.tick_base + row) * tfx_tick_row_bytes);
	default: return tick_buffers[push.tick_buffer_index][push.tick_base + row];
	}
}

TFX_BINDLESS_ACCESS tfx_gpu_particle_force tfx_load_force(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.tick_buffer_index].Load<tfx_gpu_particle_force>((push.tick_base + row) * tfx_tick_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_particle_force>((push.tick_base + row) * tfx_tick_row_bytes);
	default: return force_buffers[push.tick_buffer_index][push.tick_base + row];
	}
}

TFX_BINDLESS_ACCESS tfx_gpu_particle_vortex_fluctuation tfx_load_fluctuation(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.tick_buffer_index].Load<tfx_gpu_particle_vortex_fluctuation>((push.tick_base + row) * tfx_tick_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_particle_vortex_fluctuation>((push.tick_base + row) * tfx_tick_row_bytes);
	default: return fluctuation_buffers[push.tick_buffer_index][push.tick_base + row];
	}
}

TFX_BINDLESS_ACCESS tfx_gpu_spawn_location tfx_load_location(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.location_buffer_index].Load<tfx_gpu_spawn_location>((push.location_base + row) * tfx_location_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_spawn_location>((push.location_base + row) * tfx_location_row_bytes);
	default: return location_buffers[push.location_buffer_index][push.location_base + row];
	}
}

//A CPU spawned particle's row, loaded a part at a time so a ring without noise never reads the noise
TFX_BINDLESS_ACCESS tfx_gpu_particle_hot tfx_load_spawn_hot(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.spawn_buffer_index].Load<tfx_gpu_particle_hot>((push.spawn_base + row) * tfx_spawn_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_particle_hot>((push.spawn_base + row) * tfx_spawn_row_bytes);
	default: return spawn_buffers[push.spawn_buffer_index][push.spawn_base + row].hot;
	}
}

TFX_BINDLESS_ACCESS tfx_gpu_particle_cold tfx_load_spawn_cold(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.spawn_buffer_index].Load<tfx_gpu_particle_cold>((push.spawn_base + row) * tfx_spawn_row_bytes + tfx_spawn_cold_offset);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_particle_cold>((push.spawn_base + row) * tfx_spawn_row_bytes + tfx_spawn_cold_offset);
	default: return spawn_buffers[push.spawn_buffer_index][push.spawn_base + row].cold;
	}
}

TFX_BINDLESS_ACCESS tfx_gpu_particle_noise tfx_load_spawn_noise(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.spawn_buffer_index].Load<tfx_gpu_particle_noise>((push.spawn_base + row) * tfx_spawn_row_bytes + tfx_spawn_noise_offset);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_particle_noise>((push.spawn_base + row) * tfx_spawn_row_bytes + tfx_spawn_noise_offset);
	default: return spawn_buffers[push.spawn_buffer_index][push.spawn_base + row].noise;
	}
}

//---- A generated tick's spawn block, its runs and multipliers, rows of the tick array

TFX_BINDLESS_ACCESS tfx_gpu_spawn_block_a tfx_load_spawn_block_a(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.tick_buffer_index].Load<tfx_gpu_spawn_block_a>((push.tick_base + row) * tfx_tick_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_spawn_block_a>((push.tick_base + row) * tfx_tick_row_bytes);
	default: return spawn_block_a_buffers[push.tick_buffer_index][push.tick_base + row];
	}
}

TFX_BINDLESS_ACCESS tfx_gpu_spawn_block_b tfx_load_spawn_block_b(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.tick_buffer_index].Load<tfx_gpu_spawn_block_b>((push.tick_base + row) * tfx_tick_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_spawn_block_b>((push.tick_base + row) * tfx_tick_row_bytes);
	default: return spawn_block_b_buffers[push.tick_buffer_index][push.tick_base + row];
	}
}

TFX_BINDLESS_ACCESS tfx_gpu_spawn_block_c tfx_load_spawn_block_c(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.tick_buffer_index].Load<tfx_gpu_spawn_block_c>((push.tick_base + row) * tfx_tick_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_spawn_block_c>((push.tick_base + row) * tfx_tick_row_bytes);
	default: return spawn_block_c_buffers[push.tick_buffer_index][push.tick_base + row];
	}
}

TFX_BINDLESS_ACCESS tfx_gpu_spawn_block_d tfx_load_spawn_block_d(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.tick_buffer_index].Load<tfx_gpu_spawn_block_d>((push.tick_base + row) * tfx_tick_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_spawn_block_d>((push.tick_base + row) * tfx_tick_row_bytes);
	default: return spawn_block_d_buffers[push.tick_buffer_index][push.tick_base + row];
	}
}

TFX_BINDLESS_ACCESS tfx_gpu_spawn_block_e tfx_load_spawn_block_e(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.tick_buffer_index].Load<tfx_gpu_spawn_block_e>((push.tick_base + row) * tfx_tick_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_spawn_block_e>((push.tick_base + row) * tfx_tick_row_bytes);
	default: return spawn_block_e_buffers[push.tick_buffer_index][push.tick_base + row];
	}
}

TFX_BINDLESS_ACCESS tfx_gpu_spawn_run tfx_load_spawn_run(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.tick_buffer_index].Load<tfx_gpu_spawn_run>((push.tick_base + row) * tfx_tick_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_spawn_run>((push.tick_base + row) * tfx_tick_row_bytes);
	default: return spawn_run_buffers[push.tick_buffer_index][push.tick_base + row];
	}
}

//Just the field the run search compares, so each step loads 4 bytes rather than the whole run
TFX_BINDLESS_ACCESS uint tfx_load_spawn_run_first_spawn(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.tick_buffer_index].Load((push.tick_base + row) * tfx_tick_row_bytes + tfx_spawn_run_first_spawn_offset);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<uint>((push.tick_base + row) * tfx_tick_row_bytes + tfx_spawn_run_first_spawn_offset);
	default: return spawn_run_buffers[push.tick_buffer_index][push.tick_base + row].first_spawn;
	}
}

TFX_BINDLESS_ACCESS tfx_gpu_spawn_multipliers tfx_load_spawn_multipliers(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.tick_buffer_index].Load<tfx_gpu_spawn_multipliers>((push.tick_base + row) * tfx_tick_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_upload.Load<tfx_gpu_spawn_multipliers>((push.tick_base + row) * tfx_tick_row_bytes);
	default: return spawn_multipliers_buffers[push.tick_buffer_index][push.tick_base + row];
	}
}

//---- The graph lookups from tfx_GetGPUGraphLookupsBuffer

TFX_BINDLESS_ACCESS tfx_gpu_graph_data tfx_load_graph(tfx_gpu_particle_push push, uint row) {
	__target_switch {
	case metal: return tfx_metal_buffers.buffers[push.graph_buffer_index].Load<tfx_gpu_graph_data>(row * tfx_graph_row_bytes);
	TFX_FIXED_CASES return tfx_fixed_graphs.Load<tfx_gpu_graph_data>(row * tfx_graph_row_bytes);
	default: return graph_buffers[push.graph_buffer_index][row];
	}
}

//--- Kernel_constants

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
//tfxGPU_PARTICLE_USER_FLAGS_SHIFT and _MASK: the host's own per particle state, which the library never reads
static const uint tfx_particle_user_flags_shift = 16;
static const uint tfx_particle_user_flags_mask = 0xFFFF0000u;

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
static const uint tfx_force_location_delay = 1 << 12;

//Must match tfxDRAG_VARIATION_HASH_AXIS or the drag drawn here won't be the one the CPU draws for the same uid
static const uint tfx_drag_variation_hash_axis = 0x9e3779b9u;
static const float tfx_size_packed_scale = 32767.0 / 64.0;
static const float tfx_color_packed_scale = 32767.0 / 128.0;
static const uint tfx_capture_bit = 1 << 15;
//tfxSPAWN_LOCATION_SLOT_MASK
static const uint tfx_spawn_location_slot_mask = (1u << 20) - 1;


tfx_gpu_graph_data tfx_graph(tfx_gpu_particle_push push, uint graph_offset, uint graph) {
	return tfx_load_graph(push, graph_offset + graph);
}

//--- Instance_packing

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

//--- Emitter_forces

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
void tfx_apply_shockwave_force(tfx_gpu_particle_force force, float front, float3 position, inout float3 medium) {
	const float3 offset = position - force.origin;
	const float length_squared = dot(offset, offset);
	const float inverse_length = tfx_force_inverse_length(length_squared);
	const float distance = length_squared * inverse_length;
	const float profile_time = (front - distance) * force.falloff_scale;
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

//--- Ring_space

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
			const tfx_gpu_spawn_location location = tfx_load_location(push, dispatch.location_offset + slot);
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
//anchor taken off their origin already. A relative ring is in that frame already. A particle whose location has gone feels none of these rows
struct tfx_location_frame {
	float3 position;
	float4 rotation;
	float time;							//tfx__wide_location_time: the emitter's age as this location sees it
	bool alive;
};

tfx_location_frame tfx_particle_location_frame(tfx_particle_context context, tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold) {
	tfx_location_frame frame;
	frame.position = hot.position;
	frame.rotation = float4(0.0, 0.0, 0.0, 1.0);
	frame.time = 0.0;
	frame.alive = false;
	const uint slot = cold.spawn_location & tfx_spawn_location_slot_mask;
	if (slot < context.dispatch.location_count) {
		const tfx_gpu_spawn_location location = tfx_load_location(context.push, context.dispatch.location_offset + slot);
		frame.alive = location.local_id == cold.spawn_location;
		const float age = location.paused != 0 ? location.age : location.age - context.tick.location_age_offset;
		const float running_time = max(age - context.dispatch.location_delay, 0.0);
		const float loop_length = context.dispatch.location_loop_length;
		frame.time = loop_length > 0.0 ? running_time - floor(running_time / loop_length) * loop_length : running_time;
		if ((context.dispatch.flags & tfx_dispatch_relative_location) == 0) {
			frame.position = hot.position - location.position;
			if ((context.dispatch.flags & tfx_dispatch_location_forces_rotation) != 0) {
				frame.rotation = tfx_unpack16bit_quaternion(location.packed_rotation);
				frame.position = tfx_rotate_vector(frame.position, tfx_conjugate(frame.rotation));
			}
		}
	}
	return frame;
}

//--- Force_rows

//tfx__wide_apply_location_shockwave_force: an at_location row holds the front at location time 0, its speed in radial_ratio and its radius in axial_ratio
void tfx_apply_location_shockwave_force(tfx_gpu_particle_force force, float location_time, float3 position, inout float3 medium) {
	const float front = location_time * force.radial_ratio + force.front;
	if (front <= force.axial_ratio) {
		tfx_apply_shockwave_force(force, front, position, medium);
	}
}

//tfx__wide_location_fluctuation: an at_location row holds each harmonic's weight, turns per millisecond and turn at location time 0
tfx_gpu_particle_vortex_fluctuation tfx_location_fluctuation(tfx_gpu_particle_vortex_fluctuation row, float location_time) {
	float3 turns = row.harmonic_sin * location_time + float3(row.turn_offset_0, row.turn_offset_1, row.turn_offset_2);
	turns -= floor(turns);
	const float3 angle = turns * tfx_two_pi;
	row.harmonic_sin = row.harmonic_cos * sin(angle);
	row.harmonic_cos = row.harmonic_cos * cos(angle);
	return row;
}

//tfx_apply_forces over the tick's force rows, which the CPU packed in its order: wind, vortices, attracts, shockwaves, the same at
//locations, noise
float3 tfx_force_medium(tfx_particle_context context, tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold) {
	const uint first_row = context.tick.force_offset;
	const uint end_row = first_row + context.tick.force_count;
	bool location_forces = (context.dispatch.flags & tfx_dispatch_location_forces) != 0;
	tfx_location_frame frame = { hot.position, float4(0.0, 0.0, 0.0, 1.0), 0.0, false };
	if (location_forces) {
		frame = tfx_particle_location_frame(context, hot, cold);
	}
	float3 medium = float3(0.0, 0.0, 0.0);
	float3 location_medium = float3(0.0, 0.0, 0.0);
	for (uint row = first_row; row < end_row; ++row) {
		const tfx_gpu_particle_force force = tfx_load_force(context.push, row);
		const uint type = force.type_flags & tfx_force_type_mask;
		if (type == tfx_force_wind) {
			if ((force.type_flags & tfx_force_location_delay) == 0) {
				//The CPU seeds the medium with the summed wind rather than adding it
				medium = force.axis;
			} else if (!frame.alive || frame.time >= force.front) {
				medium += force.axis;
			}
			continue;
		}
		if (type == tfx_force_noise) {
			//The location flow joins before the noise, as on the CPU
			if (location_forces && frame.alive) {
				medium += tfx_rotate_vector(location_medium, frame.rotation);
			}
			const bool noise_at_location = location_forces && (force.type_flags & tfx_force_at_location) != 0;
			//tfx__wide_location_delay_reached: a gone location's particles keep the noise
			const bool noise_on = (force.type_flags & tfx_force_location_delay) == 0 || !frame.alive || frame.time >= force.front;
			location_forces = false;
			if (noise_on) {
				medium += tfx_noise_medium(context.push, context.dispatch, force, context.tick, hot, cold, context.noise, context.life, noise_at_location ? frame.position : hot.position, !noise_at_location || frame.alive);
			}
			continue;
		}
		const bool at_location = location_forces && (force.type_flags & tfx_force_at_location) != 0;
		const float3 position = at_location ? frame.position : hot.position;
		float3 target = at_location ? location_medium : medium;
		//An at_location attract or vortex holds its delay in front, waited out from the particle's location time
		const bool switched_on = !at_location || frame.time >= force.front;
		if (type == tfx_force_vortex) {
			tfx_gpu_particle_vortex_fluctuation fluctuation = {};
			if ((force.type_flags & tfx_force_fluctuating) != 0) {
				fluctuation = tfx_load_fluctuation(context.push, ++row);
				if (at_location) {
					fluctuation = tfx_location_fluctuation(fluctuation, frame.time);
				}
			}
			if (switched_on) {
				tfx_apply_vortex_force(force, fluctuation, position, target);
			}
		} else if (type == tfx_force_attract) {
			if (switched_on) {
				tfx_apply_attract_force(force, position, target);
			}
		} else if (type == tfx_force_shockwave) {
			if (at_location) {
				tfx_apply_location_shockwave_force(force, frame.time, position, target);
			} else {
				tfx_apply_shockwave_force(force, force.front, position, target);
			}
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

//--- Host_hooks

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

//A live particle's instance in world space before it's packed
struct tfx_particle_draw {
	float3 position;
	float stretch;						//The instance's position.w, the stretch graph
	float4 rotation;					//Quaternion xyzw
	float3 alignment;					//Unit vector the billboard aligns to, zero for none
	float intensity;					//Packed 16 bit, within +-128
	float2 size;						//Overall scale applied, packed 16 bit within +-64
	float gradient_map;					//Where the color ramp is sampled
	float curved_alpha;
	float alpha_sharpness;
	uint image_frame;					//Must stay below the emitter's frame count
};

//Changes how a live particle is drawn after its world transform, on the last tick of the fetch only. hot is in ring space
interface ITfxParticleWrite {
	void before_write(inout tfx_particle_draw draw, tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, tfx_particle_context context);
}

struct tfx_no_particle_write : ITfxParticleWrite {
	void before_write(inout tfx_particle_draw draw, tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, tfx_particle_context context) {
	}
};

//--- Particle_control

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

//--- Instance_write

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
float4 tfx_particle_quaternion(tfx_gpu_particle_push push, tfx_gpu_particle_dispatch dispatch, tfx_gpu_particle_cold cold, float life, uint2 location_rotation) {
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
			return tfx_location_roll_quaternion(roll, tfx_unpack16bit_quaternion(location_rotation));
		}
		return tfx_euler_to_quaternion(0.0, 0.0, roll);
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
	return q;
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
	return uint(int(image_frame));
}

//A removed particle still holding its slot: uncaptured, so the vertex shader zeroes size and stretch, and nothing reads its position, only a valid quaternion and property index
void tfx_write_empty_instance(tfx_gpu_particle_push push, tfx_gpu_particle_dispatch dispatch, uint ring_position) {
	const uint instance_index = dispatch.instance_start + ring_position;
	tfx_instance instance;
	instance.position = float4(0.0, 0.0, 0.0, 0.0);
	instance.quaternion = uint2(0u, 0x7FFF0000u);
	instance.size = 0;
	instance.alignment = 0;
	instance.intensity_gradient_map = 0;
	instance.curved_alpha_life = 0;
	instance.indexes = dispatch.property_index;
	instance.captured_index = instance_index;
	tfx_store_instance(push, instance_index, instance);
}

//tfx__control_particle_instances for one particle, at the context's ring position on the last tick
void tfx_write_particle_instance<TWrite : ITfxParticleWrite>(tfx_particle_context context, tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, float3 position, uint2 location_rotation, TWrite write_hook) {
	const tfx_gpu_particle_push push = context.push;
	const tfx_gpu_particle_dispatch dispatch = context.dispatch;
	const tfx_gpu_particle_tick tick = context.tick;
	const float life = context.life;
	const uint ring_position = context.ring_position;
	const uint flags = dispatch.flags;
	const uint graph_offset = dispatch.graph_offset;
	const bool alive = (hot.flags & tfx_particle_remove) == 0;
	if (!alive) {
		tfx_write_empty_instance(push, dispatch, ring_position);
		return;
	}
	tfx_particle_draw draw;
	draw.position = position;

	//----Transform
	const float stretch = tfx_sample_graph(tfx_graph(push, graph_offset, tfx_graph_stretch), life) * dispatch.global_stretch;
	draw.stretch = stretch;
	draw.alignment = float3(0.0, 0.0, 0.0);
	if ((flags & tfx_dispatch_align_emission) != 0) {
		float3 alignment = tfx_unpack10bit(cold.velocity_normal);
		if ((flags & tfx_dispatch_align_rotated_emission) != 0) {
			alignment = tfx_rotate_vector(alignment, dispatch.emitter_rotation);
		}
		if ((flags & tfx_dispatch_location_rotation) != 0) {
			alignment = tfx_rotate_vector(alignment, tfx_unpack16bit_quaternion(location_rotation));
		}
		draw.alignment = alignment;
	} else if ((flags & tfx_dispatch_align_emitter) != 0) {
		float3 alignment = tfx_rotate_vector(float3(0.0, 1.0, 0.0), dispatch.emitter_rotation);
		if ((flags & tfx_dispatch_location_rotation) != 0) {
			alignment = tfx_rotate_vector(alignment, tfx_unpack16bit_quaternion(location_rotation));
		}
		draw.alignment = alignment;
	}

	//----Spin
	draw.rotation = tfx_particle_quaternion(push, dispatch, cold, life, location_rotation);

	//----Size
	draw.size = float2(0.0, 0.0);
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
		draw.size = float2(scale_x, scale_y) * tick.overall_scale;
	}

	//----Color
	draw.intensity = dispatch.global_intensity * tfx_sample_graph(tfx_graph(push, graph_offset, tfx_graph_intensity), life) * cold.intensity_factor;
	draw.curved_alpha = tfx_sample_graph(tfx_graph(push, graph_offset, tfx_graph_curved_alpha), life);
	draw.alpha_sharpness = tfx_sample_graph(tfx_graph(push, graph_offset, tfx_graph_alpha_sharpness), life);
	draw.gradient_map = (flags & tfx_dispatch_random_color) != 0 ? cold.random_color : tfx_sample_graph(tfx_graph(push, graph_offset, tfx_graph_gradient_mapper), life);
	draw.image_frame = tfx_particle_image_frame(dispatch, cold, hot, tick);

	write_hook.before_write(draw, hot, cold, context);

	tfx_instance instance;
	instance.position = float4(draw.position, draw.stretch);
	instance.quaternion = tfx_pack16bit_quaternion(draw.rotation);
	instance.size = tfx_pack16bit_2sscaled(draw.size.x, draw.size.y, tfx_size_packed_scale);
	instance.alignment = tfx_pack8bit_xyz(draw.alignment);
	instance.intensity_gradient_map = tfx_pack16bit_2sscaled(draw.intensity, draw.gradient_map, tfx_color_packed_scale);
	instance.curved_alpha_life = tfx_pack8bitunorm_xyz(float3(draw.curved_alpha, draw.alpha_sharpness, life));

	//----Image frame and interpolation. Only particles that were in the ring at the previous write have somewhere to come from
	const uint instance_index = dispatch.instance_start + ring_position;
	const bool capture = (flags & tfx_dispatch_capture) != 0 && ring_position + dispatch.spawned_since_write < tick.count;
	instance.indexes = dispatch.property_index | draw.image_frame | (capture ? tfx_capture_bit : 0);
	instance.captured_index = capture ? dispatch.previous_instance_start + ring_position + dispatch.bumped_since_write : instance_index;

	tfx_store_instance(push, instance_index, instance);
}

//--- Spawn_ingest

//A spawn the kernel builds from its tick's block rather than loading a row, see tfx__record_gpu_groups
void tfx_generate_particle<TSpawn : ITfxParticleSpawn>(tfx_gpu_particle_push push, tfx_gpu_particle_dispatch dispatch, tfx_gpu_particle_tick tick, uint spawn_index, TSpawn spawn_hook, out tfx_gpu_particle_hot hot, out tfx_gpu_particle_cold cold, out tfx_gpu_particle_noise noise) {
	const uint block_row = tick.upload_offset;
	const tfx_gpu_spawn_block_a block_a = tfx_load_spawn_block_a(push, block_row);
	const tfx_gpu_spawn_block_b block_b = tfx_load_spawn_block_b(push, block_row + 1);
	const tfx_gpu_spawn_block_c block_c = tfx_load_spawn_block_c(push, block_row + 2);
	tfx_gpu_spawn_block_d block_d = {};
	if ((block_a.flags & tfx_spawn_shape) != 0) {
		block_d = tfx_load_spawn_block_d(push, block_row + 3);
	}
	tfx_gpu_spawn_block_e block_e = {};
	if ((block_a.flags & (tfx_spawn_grid | tfx_spawn_stepped_emission | tfx_spawn_stepped_roll)) != 0) {
		block_e = tfx_load_spawn_block_e(push, block_row + 4);
	}
	//The last run starting at or before the spawn. Neighbouring threads are in the same run, so the search stays coherent
	const uint first_run = block_row + block_a.run_row;
	uint low = 0;
	uint high = block_a.run_count - 1;
	while (low < high) {
		const uint middle = (low + high + 1) >> 1;
		if (tfx_load_spawn_run_first_spawn(push, first_run + middle) <= spawn_index) {
			low = middle;
		} else {
			high = middle - 1;
		}
	}
	const tfx_gpu_spawn_run run = tfx_load_spawn_run(push, first_run + low);
	const bool has_multipliers = run.multipliers_row != 0xFFFFFFFFu;
	tfx_gpu_spawn_multipliers multipliers = {};
	if (has_multipliers) {
		multipliers = tfx_load_spawn_multipliers(push, block_row + run.multipliers_row);
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
	hot.flags = spawn.user_flags & tfx_particle_user_flags_mask;
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
	cold.user_value = spawn.user_value;
	noise.offset = spawn.noise_offset;
	noise.resolution = spawn.noise_resolution;
}

//position % capacity for a position below 2 * capacity, which every ring index sum here is, without an integer divide
uint tfx_wrap_ring(uint position, uint capacity) {
	return position >= capacity ? position - capacity : position;
}

//--- Kernel_body

//The kernel for one ring, push.dispatch_index, with a thread per slot. Dispatch in groups of tfx_group_size
void tfx_simulate<TUpdate : ITfxParticleUpdate, TSpawn : ITfxParticleSpawn, TWrite : ITfxParticleWrite>(uint thread_index, tfx_gpu_particle_push push, TUpdate update, TSpawn spawn_hook, TWrite write_hook) {
	const tfx_gpu_particle_dispatch dispatch = tfx_load_dispatch(push, push.dispatch_index);
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
		const tfx_gpu_particle_tick tick = tfx_load_tick(push, tick_index);
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
				const uint row = tick.upload_offset + ring_position - first_spawn;
				hot = tfx_load_spawn_hot(push, row);
				cold = tfx_load_spawn_cold(push, row);
				if (has_noise) {
					noise = tfx_load_spawn_noise(push, row);
				}
			}
			ingested = true;
			dead = false;
		} else if (!live) {
			hot = tfx_load_hot(push, pool_slot);
			dead = (hot.flags & tfx_particle_remove) != 0;
			if (!dead) {
				cold = tfx_load_cold(push, pool_slot);
				if (has_noise) {
					noise = tfx_load_noise(push, pool_slot);
				}
			}
		}
		live = true;

		const bool last_tick = tick_index + 1 == end_tick;
		if (dead) {
			if (last_tick) {
				tfx_write_empty_instance(push, dispatch, ring_position);
			}
			continue;
		}
		tfx_particle_context context;
		context.push = push;
		context.dispatch = dispatch;
		context.tick = tick;
		context.noise = noise;
		context.life = hot.age * cold.inv_max_age;
		context.ring_position = ring_position;
		context.last_tick = last_tick;
		//An expired or killed particle holds its slot until its batch expires, and is written with no size
		if ((hot.flags & tfx_particle_remove) == 0) {
			tfx_control_particle(context, hot, cold, velocity_graph, weight_graph, drag_variation, update);
			update.after_move(hot, cold, context);
		}
		if (last_tick) {
			uint2 location_rotation;
			const float3 position = tfx_particle_world_position(push, dispatch, tick, hot, cold, location_rotation);
			tfx_write_particle_instance(context, hot, cold, position, location_rotation, write_hook);
		}
		tfx_age_particle(hot, cold, tick);
	}

	if (live && !dead) {
		tfx_store_hot(push, pool_slot, hot);
		if (ingested) {
			tfx_store_cold(push, pool_slot, cold);
			if (has_noise) {
				tfx_store_noise(push, pool_slot, noise);
			}
		}
	}
}

//For a host with an update and a spawn hook but no write hook
void tfx_simulate<TUpdate : ITfxParticleUpdate, TSpawn : ITfxParticleSpawn>(uint thread_index, tfx_gpu_particle_push push, TUpdate update, TSpawn spawn_hook) {
	tfx_simulate(thread_index, push, update, spawn_hook, tfx_no_particle_write());
}

//For a host with an update but no spawn or write hook
void tfx_simulate<TUpdate : ITfxParticleUpdate>(uint thread_index, tfx_gpu_particle_push push, TUpdate update) {
	tfx_simulate(thread_index, push, update, tfx_no_particle_spawn(), tfx_no_particle_write());
}

//--- Entry_point

tfx_gpu_particle_push tfx_entry_push_constants() {
	__target_switch {
	TFX_UNIFORM_PUSH_CASES return tfx_uniform_push_constants;
	default: return tfx_push_constants;
	}
}

[shader("compute")]
[numthreads(tfx_group_size, 1, 1)]
void simulate_particles(uint3 thread_id : SV_DispatchThreadID) {
	tfx_simulate<tfx_no_particle_update>(thread_id.x, tfx_entry_push_constants(), tfx_no_particle_update());
}
