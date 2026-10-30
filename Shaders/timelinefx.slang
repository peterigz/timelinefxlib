/*
	TimelineFX shaders

	The library's reference shaders, one file a host can copy into its renderer and compile with Slang, one entry point at a time: the
	compute kernel for Run on GPU emitters, which compiles for Vulkan, D3D12, D3D11, Metal, OpenGL or WebGPU, and the billboard, ribbon,
	baked playback and bounds passes, which use the bindless model. Search for a keyword below to jump to its section (in Vim, * on the
	keyword):

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

	--The render passes
	[Render_rows]                       Row and push constant structs, byte for byte mirrors of their C structs
	[Render_access]                     Push constant blocks, bindless arrays and the accessors every read and write goes through
		[Fixed_render_bindings]         The D3D11 registers every render read and write goes to under TFX_FIXED_BINDINGS
	[Render_math]                       Unpacking and packing, the quaternion matrix and the billboard rotation
	[Billboard_vertex]                  tfx_billboard_vertex, the live billboard draw
	[Billboard_baked_vertex]            tfx_billboard_baked_vertex, the baked billboard draw
	[Particle_shading]                  tfx_billboard_fragment for live and baked billboards, and tfx_ribbon_fragment
	[Ribbon_shape]                      Ribbon graph lookups, gradient noise and path sampling
	[Ribbon_build]                      tfx_ribbon_compute, a bucket's ribbon vertexes and indexes
	[Ribbon_draw]                       tfx_ribbon_vertex
	[Sprite_playback]                   tfx_sprite_playback, a baked frame's sprites into billboard instances
	[Ribbon_playback]                   tfx_ribbon_playback, a baked frame's ribbons into the ribbons tfx_ribbon_compute builds
	[Bounds_reduction]                  tfx_sprite_bounds and tfx_ribbon_bounds, a baked frame's bounding box
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

//f16tof32 and f32tof16 lower to 16 bit types on SPIR-V, which need shaderInt16 and shaderFloat16. The GLSL.std.450 pair needs neither
float tfx_half_to_float(uint value) {
	__target_switch {
	case spirv:
		return spirv_asm {
			%unpacked:$$float2 = OpExtInst glsl450 UnpackHalf2x16 $value;
			result:$$float = OpCompositeExtract %unpacked 0
		};
	default:
		return f16tof32(value);
	}
}

//The half in the low 16 bits, the high 16 bits zero
uint tfx_float_to_half(float value) {
	__target_switch {
	case spirv:
		float2 packing = float2(value, 0.0);
		return spirv_asm {
			result:$$uint = OpExtInst glsl450 PackHalf2x16 $packing
		};
	default:
		return f32tof16(value);
	}
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
[[vk::push_constant]] ConstantBuffer<tfx_gpu_particle_push> tfx_push_constants TFX_FIXED_REGISTER(b0);
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

//--- Render_rows

//The rows the draw, ribbon and baked playback passes read, each the byte for byte mirror of its C struct, padded where C pads

//Mirrors tfx_gpu_image_data_t
struct tfx_image_data {
	float4 uv;
	uint2 uv_packed;
	float2 image_size;
	uint texture_array_index;			//[bindless image index << 16 | animation frame]
	float animation_frames;
	float padding[2];
};

//Mirrors tfx_gpu_particle_properties_t
struct tfx_emitter_properties {
	float2 image_handle;
	uint color_ramp_indexes;			//[color ramp layer << 8 | color ramp row]
	uint flags;							//The billboard alignment in the lowest 2 bits
	uint start_frame_index;
	float animation_frames;
	int padding[2];
};

//Mirrors tfx_uniform_buffer_data_t. column_major because the Slang API's sessions default to row major, where slangc defaults to column major
struct tfx_view_uniforms {
	column_major float4x4 view;
	column_major float4x4 projection;
	float2 screen_size;
	uint milliseconds;
	float timer_lerp;
	float update_time;
};

//Mirrors tfx_ribbon_t
struct tfx_ribbon {
	float3 position;
	uint age_width;						//Normalised age as a unorm16 in the low half, width as a half float in the high half
	uint2 quaternion;					//16 bit snorm: .x = X|Y, .y = Z|W
	uint start_index_flags;				//Segment start index in the low 24 bits, flags in the top 8
	uint emitter_index_phase_seed;		//Emitter index in the low 16 bits, phase offset seed in the high 16
	uint texture_indexes;				//[color ramp row << 24 | color ramp layer << 16 | image frame]
	uint intensity_gradient_map;		//Two 16 bit sscaled values
	uint curved_alpha;					//Curved alpha and sharpness, two 16 bit sscaled values
	uint scale;							//The parent effect's overall scale as a half float in the low 16 bits
};

//Mirrors tfx_ribbon_vertex_t
struct tfx_ribbon_vertex_row {
	float3 position;
	uint segment_index;					//[tessellation << 29 | tessellation index << 26 | side << 25 | position along the ribbon << 12 | segment]
	float2 uv_offset_scale;
	uint ribbon_index;
	uint clipped;
};

//tfxRIBBON_LAG_SPINE_SAMPLES
static const int tfx_lag_spine_samples = 16;

//Mirrors tfx_gpu_ribbon_emitter_t, 416 bytes
struct tfx_ribbon_emitter {
	uint2 quaternion;
	uint lookup_offset;
	uint angle_type;
	float3 position;
	float age;
	float3 captured_position;
	uint morph_segment_start_index;		//tfx_invalid_index when the emitter doesn't morph
	float noise_frequency;
	float noise_speed;
	float noise_phase_range;
	float noise_lock_rate;
	float3 fixed_angle_normal;
	uint sample_count_noise;			//[noise algorithm and octaves << 24 | sample count]
	uint start_frame_index;
	float lag_time;
	float lag_span;
	float lag_spine_position_x[tfx_lag_spine_samples];
	float lag_spine_position_y[tfx_lag_spine_samples];
	float lag_spine_position_z[tfx_lag_spine_samples];
	uint lag_spine_padding;
	uint2 lag_spine_quaternion[tfx_lag_spine_samples];
};

//Mirrors tfx_animation_instance_t
struct tfx_animation_instance {
	float3 position;
	float scale;
	uint sprite_count;
	uint ribbon_count;
	uint frame_count;
	uint offset_into_sprite_data;
	uint offset_into_ribbon_data;
	uint info_index;
	float current_time;
	float animation_length_in_time;
	float tween;
	uint flags;
	uint padding0;
	uint padding1;
};

//Mirrors tfx_ribbon_instance_data_t, a baked ribbon
struct tfx_ribbon_instance_data {
	float4 position;					//w is the normalised age
	float width;
	uint start_index;
	uint flags;
	uint captured_index;				//The same ribbon in the previous frame, tfx_invalid_index in its first
	uint2 quaternion;
	float scale;
	uint texture_indexes;
	uint intensity_gradient_map;
	uint curved_alpha;
	uint additional;					//[property index << 16 | lerp offset as a unorm16]
	uint phase_seed;
};

//Mirrors tfx_animation_ribbon_properties_t
struct tfx_animation_ribbon_properties {
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
	uint morph_segment_start_index;		//Absolute, tfx_invalid_index when the emitter doesn't morph
	float noise_frequency;
	float noise_speed;
	float noise_phase_range;
	float noise_lock_rate;
	uint noise_packed;
};

//Mirrors the editor's sprite_data_bb, one box per bounds group
struct tfx_bounding_box {
	float3 min_corner;
	float max_stretch;
	float3 max_corner;
	float padding;
};

//Mirrors tfx_push_constants_t in zest's impl_timelinefx.h, the live billboard draw's push constants
struct tfx_billboard_push {
	uint color_ramp_texture_index;
	uint image_data_index;
	uint properties_index;
	uint sampler_index;
	uint prev_billboards_index;
	uint index_offset;					//Added to every captured index, for a draw whose previous instances don't start at 0
	uint uniform_index;
};

//Mirrors the editor's tfx_sprite_sheet_push_t. It starts with tfx_billboard_push's fields so tfx_billboard_fragment serves both draws
struct tfx_baked_billboard_push {
	uint color_ramp_texture_index;
	uint image_data_index;
	uint properties_index;
	uint sampler_index;
	uint prev_billboards_index;
	uint index_offset;
	uint uniform_index;
	float ndc_offset_x;					//Clip space offsets for drawing a seamless sprite sheet's neighbouring tiles
	float ndc_offset_y;
};

//Mirrors tfx_ribbon_bucket_globals_t, the push constants of all three ribbon passes, set per bucket by the library
struct tfx_ribbon_bucket_globals {
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
	float lerp;
	float time;
	float ndc_offset_x;
	float ndc_offset_y;
	uint emitter_offset;
};

//Mirrors tfx_sprite_data_push_t, for the sprite playback and the sprite bounds
struct tfx_sprite_data_push {
	uint animation_instances_total;
	uint billboards_total;
	uint animated_shapes;
	uint offsets_index;
	uint animation_instances_index;
	uint billboards_index;
	uint sprite_data_index;
	uint image_data_index;
	uint emitter_properties_index;
	uint bounding_boxes_index;
};

//Mirrors tfx_ribbon_data_push_t
struct tfx_ribbon_data_push {
	uint animation_instances_total;
	uint ribbons_total;
	uint offsets_index;
	uint animation_instances_index;
	uint ribbon_data_index;
	uint ribbon_output_index;
	uint ribbon_properties_index;
};

//Mirrors the editor's ribbon_bb_push_t
struct tfx_ribbon_bounds_push {
	uint ribbons_total;
	uint animation_instances_index;
	uint ribbon_data_index;
	uint ribbon_segments_index;
	uint ribbon_properties_index;
	uint bounding_boxes_index;
	uint graphs_index;
};

//--- Render_access

//Every buffer, image and push constant block the draw and playback passes use. The descriptor layout is the editor's: samplers at
//binding 0, images at 3, storage buffers at 5 and uniform buffers at 7, each an array indexed from the push constants

//---- Push constants, one block per pass

//Under TFX_FIXED_BINDINGS the billboard draws take tfx_fixed_draw in place of theirs, which are all indexes but for the offsets
[[vk::push_constant]] ConstantBuffer<tfx_billboard_push> tfx_billboard_push_constants;
[[vk::push_constant]] ConstantBuffer<tfx_baked_billboard_push> tfx_baked_billboard_push_constants;
[[vk::push_constant]] ConstantBuffer<tfx_ribbon_bucket_globals> tfx_ribbon_push_constants TFX_FIXED_REGISTER(b2);
[[vk::push_constant]] ConstantBuffer<tfx_sprite_data_push> tfx_sprite_data_push_constants TFX_FIXED_REGISTER(b3);
[[vk::push_constant]] ConstantBuffer<tfx_ribbon_data_push> tfx_ribbon_data_push_constants TFX_FIXED_REGISTER(b4);
[[vk::push_constant]] ConstantBuffer<tfx_ribbon_bounds_push> tfx_ribbon_bounds_push_constants;

//---- Bindless arrays

[[vk::binding(0, 0)]] SamplerState tfx_samplers[];
[[vk::binding(3, 0)]] Texture2DArray tfx_images[];
[[vk::binding(7, 0)]] ConstantBuffer<tfx_view_uniforms> tfx_view_uniform_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_instance> tfx_instance_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_image_data> tfx_image_data_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_emitter_properties> tfx_emitter_properties_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_ribbon> tfx_ribbon_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_ribbon_emitter> tfx_ribbon_emitter_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<float4> tfx_ribbon_segment_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<uint> tfx_offset_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_animation_instance> tfx_animation_instance_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_ribbon_instance_data> tfx_baked_ribbon_buffers[];
[[vk::binding(5, 0)]] StructuredBuffer<tfx_animation_ribbon_properties> tfx_ribbon_properties_buffers[];
[[vk::binding(5, 0)]] RWStructuredBuffer<tfx_ribbon> tfx_ribbon_output_buffers[];
[[vk::binding(5, 0)]] RWStructuredBuffer<tfx_ribbon_vertex_row> tfx_ribbon_vertex_buffers[];
[[vk::binding(5, 0)]] RWStructuredBuffer<uint> tfx_ribbon_index_buffers[];
[[vk::binding(5, 0)]] RWStructuredBuffer<tfx_bounding_box> tfx_bounding_box_buffers[];

//--- Fixed_render_bindings

//One register per kind of row, the layout the D3D11 example binds. Every shape is in one atlas, so the bindless image index is ignored
#if defined(TFX_FIXED_BINDINGS)
//SV_InstanceID leaves out the draw's first instance on D3D11, so it comes in here with the billboard push constants that aren't indexes
struct tfx_fixed_draw_constants {
	uint instance_offset;
	uint index_offset;
	float ndc_offset_x;
	float ndc_offset_y;
};

ConstantBuffer<tfx_view_uniforms> tfx_fixed_view : register(b0);
ConstantBuffer<tfx_fixed_draw_constants> tfx_fixed_draw : register(b1);
ByteAddressBuffer tfx_fixed_draw_instances : register(t0);										//The instances a billboard draw draws
ByteAddressBuffer tfx_fixed_instance_rows : register(t1);										//The previous instances, or the baked sprites
StructuredBuffer<tfx_image_data> tfx_fixed_image_data : register(t2);
StructuredBuffer<tfx_emitter_properties> tfx_fixed_emitter_properties : register(t3);
Texture2DArray tfx_fixed_particle_images : register(t4);										//The shape atlas
Texture2DArray tfx_fixed_color_ramps : register(t5);
StructuredBuffer<tfx_ribbon_emitter> tfx_fixed_ribbon_emitters : register(t6);
ByteAddressBuffer tfx_fixed_graph_rows : register(t7);
StructuredBuffer<tfx_ribbon> tfx_fixed_ribbons : register(t8);
StructuredBuffer<float4> tfx_fixed_ribbon_segments : register(t9);
StructuredBuffer<tfx_ribbon_vertex_row> tfx_fixed_ribbon_vertices : register(t10);
ByteAddressBuffer tfx_fixed_offsets : register(t11);
ByteAddressBuffer tfx_fixed_animation_instances : register(t12);
StructuredBuffer<tfx_ribbon_instance_data> tfx_fixed_baked_ribbons : register(t13);
StructuredBuffer<tfx_animation_ribbon_properties> tfx_fixed_ribbon_properties : register(t14);
RWStructuredBuffer<tfx_ribbon_vertex_row> tfx_fixed_ribbon_vertex_output : register(u0);
RWByteAddressBuffer tfx_fixed_ribbon_index_output : register(u1);
RWByteAddressBuffer tfx_fixed_instance_output : register(u2);
RWStructuredBuffer<tfx_ribbon> tfx_fixed_ribbon_output : register(u3);
SamplerState tfx_fixed_sampler : register(s0);

static const uint tfx_animation_instance_row_bytes = 64;
#define TFX_RIBBON_EMITTER_ROW(buffer_index, row) tfx_fixed_ribbon_emitters[row]
#else
#define TFX_RIBBON_EMITTER_ROW(buffer_index, row) tfx_ribbon_emitter_buffers[buffer_index][row]
#endif

//---- Accessors

#if defined(TFX_FIXED_BINDINGS)
tfx_view_uniforms tfx_load_view(uint buffer_index) {
	return tfx_fixed_view;
}

tfx_instance tfx_load_draw_instance(uint instance_id) {
	return tfx_fixed_draw_instances.Load<tfx_instance>((instance_id + tfx_fixed_draw.instance_offset) * tfx_instance_row_bytes);
}

tfx_ribbon_vertex_row tfx_load_draw_ribbon_vertex(uint vertex_id) {
	return tfx_fixed_ribbon_vertices[vertex_id];
}

tfx_instance tfx_load_instance_row(uint buffer_index, uint row) {
	return tfx_fixed_instance_rows.Load<tfx_instance>(row * tfx_instance_row_bytes);
}

void tfx_store_instance_row(uint buffer_index, uint row, tfx_instance instance) {
	tfx_fixed_instance_output.Store<tfx_instance>(row * tfx_instance_row_bytes, instance);
}

tfx_image_data tfx_load_image_data(uint buffer_index, uint row) {
	return tfx_fixed_image_data[row];
}

tfx_emitter_properties tfx_load_emitter_properties(uint buffer_index, uint row) {
	return tfx_fixed_emitter_properties[row];
}

tfx_gpu_graph_data tfx_load_graph_row(uint buffer_index, uint row) {
	return tfx_fixed_graph_rows.Load<tfx_gpu_graph_data>(row * tfx_graph_row_bytes);
}

tfx_ribbon tfx_load_ribbon(uint buffer_index, uint row) {
	return tfx_fixed_ribbons[row];
}

void tfx_store_ribbon(uint buffer_index, uint row, tfx_ribbon ribbon) {
	tfx_fixed_ribbon_output[row] = ribbon;
}

float3 tfx_load_ribbon_segment(uint buffer_index, uint row) {
	return tfx_fixed_ribbon_segments[row].xyz;
}

void tfx_store_ribbon_vertex(uint buffer_index, uint row, tfx_ribbon_vertex_row vertex) {
	tfx_fixed_ribbon_vertex_output[row] = vertex;
}

void tfx_store_ribbon_triangle(uint buffer_index, uint row, uint3 triangle) {
	tfx_fixed_ribbon_index_output.Store3(row * 4, triangle);
}

uint tfx_load_offset(uint buffer_index, uint row) {
	return tfx_fixed_offsets.Load(row * 4);
}

tfx_animation_instance tfx_load_animation_instance(uint buffer_index, uint row) {
	return tfx_fixed_animation_instances.Load<tfx_animation_instance>(row * tfx_animation_instance_row_bytes);
}

tfx_ribbon_instance_data tfx_load_baked_ribbon(uint buffer_index, uint row) {
	return tfx_fixed_baked_ribbons[row];
}

tfx_animation_ribbon_properties tfx_load_ribbon_properties(uint buffer_index, uint row) {
	return tfx_fixed_ribbon_properties[row];
}

float4 tfx_sample_particle_image(uint image_index, uint sampler_index, float3 tex_coord) {
	return tfx_fixed_particle_images.Sample(tfx_fixed_sampler, tex_coord);
}

float4 tfx_load_color_ramp(uint texture_index, int3 texel) {
	return tfx_fixed_color_ramps.Load(int4(texel, 0));
}

tfx_billboard_push tfx_billboard_entry_push() {
	tfx_billboard_push push = {};
	push.index_offset = tfx_fixed_draw.index_offset;
	return push;
}

tfx_baked_billboard_push tfx_baked_billboard_entry_push() {
	tfx_baked_billboard_push push = {};
	push.ndc_offset_x = tfx_fixed_draw.ndc_offset_x;
	push.ndc_offset_y = tfx_fixed_draw.ndc_offset_y;
	return push;
}
#else
TFX_BINDLESS_ACCESS tfx_view_uniforms tfx_load_view(uint buffer_index) {
	return tfx_view_uniform_buffers[buffer_index];
}

TFX_BINDLESS_ACCESS tfx_instance tfx_load_instance_row(uint buffer_index, uint row) {
	return tfx_instance_buffers[buffer_index][row];
}

TFX_BINDLESS_ACCESS void tfx_store_instance_row(uint buffer_index, uint row, tfx_instance instance) {
	instance_buffers[buffer_index][row] = instance;
}

TFX_BINDLESS_ACCESS tfx_image_data tfx_load_image_data(uint buffer_index, uint row) {
	return tfx_image_data_buffers[buffer_index][row];
}

TFX_BINDLESS_ACCESS tfx_emitter_properties tfx_load_emitter_properties(uint buffer_index, uint row) {
	return tfx_emitter_properties_buffers[buffer_index][row];
}

TFX_BINDLESS_ACCESS tfx_gpu_graph_data tfx_load_graph_row(uint buffer_index, uint row) {
	return graph_buffers[buffer_index][row];
}

TFX_BINDLESS_ACCESS tfx_ribbon tfx_load_ribbon(uint buffer_index, uint row) {
	return tfx_ribbon_buffers[buffer_index][row];
}

TFX_BINDLESS_ACCESS void tfx_store_ribbon(uint buffer_index, uint row, tfx_ribbon ribbon) {
	tfx_ribbon_output_buffers[buffer_index][row] = ribbon;
}

TFX_BINDLESS_ACCESS float3 tfx_load_ribbon_segment(uint buffer_index, uint row) {
	return tfx_ribbon_segment_buffers[buffer_index][row].xyz;
}

TFX_BINDLESS_ACCESS void tfx_store_ribbon_vertex(uint buffer_index, uint row, tfx_ribbon_vertex_row vertex) {
	tfx_ribbon_vertex_buffers[buffer_index][row] = vertex;
}

TFX_BINDLESS_ACCESS void tfx_store_ribbon_triangle(uint buffer_index, uint row, uint3 triangle) {
	tfx_ribbon_index_buffers[buffer_index][row] = triangle.x;
	tfx_ribbon_index_buffers[buffer_index][row + 1] = triangle.y;
	tfx_ribbon_index_buffers[buffer_index][row + 2] = triangle.z;
}

TFX_BINDLESS_ACCESS uint tfx_load_offset(uint buffer_index, uint row) {
	return tfx_offset_buffers[buffer_index][row];
}

TFX_BINDLESS_ACCESS tfx_animation_instance tfx_load_animation_instance(uint buffer_index, uint row) {
	return tfx_animation_instance_buffers[buffer_index][row];
}

TFX_BINDLESS_ACCESS tfx_ribbon_instance_data tfx_load_baked_ribbon(uint buffer_index, uint row) {
	return tfx_baked_ribbon_buffers[buffer_index][row];
}

TFX_BINDLESS_ACCESS tfx_animation_ribbon_properties tfx_load_ribbon_properties(uint buffer_index, uint row) {
	return tfx_ribbon_properties_buffers[buffer_index][row];
}

//Every shape has its own image, so the index varies per particle rather than per draw and the descriptor read is non uniform
TFX_BINDLESS_ACCESS float4 tfx_sample_particle_image(uint image_index, uint sampler_index, float3 tex_coord) {
	return tfx_images[NonUniformResourceIndex(image_index)].Sample(tfx_samplers[sampler_index], tex_coord);
}

TFX_BINDLESS_ACCESS float4 tfx_load_color_ramp(uint texture_index, int3 texel) {
	return tfx_images[texture_index].Load(int4(texel, 0));
}

tfx_billboard_push tfx_billboard_entry_push() {
	return tfx_billboard_push_constants;
}

tfx_baked_billboard_push tfx_baked_billboard_entry_push() {
	return tfx_baked_billboard_push_constants;
}
#endif

//The ribbon emitter up to its lag spine. The ribbon passes read the 416 byte row field by field rather than copying it
struct tfx_ribbon_emitter_frame {
	uint2 quaternion;
	uint lookup_offset;
	uint angle_type;
	float3 position;
	float3 captured_position;
	uint morph_segment_start_index;
	float noise_frequency;
	float noise_speed;
	float noise_phase_range;
	float noise_lock_rate;
	float3 fixed_angle_normal;
	uint sample_count_noise;
	float lag_time;
	float lag_span;
};

TFX_BINDLESS_ACCESS tfx_ribbon_emitter_frame tfx_load_ribbon_emitter_frame(uint buffer_index, uint row) {
	tfx_ribbon_emitter_frame frame;
	frame.quaternion = TFX_RIBBON_EMITTER_ROW(buffer_index, row).quaternion;
	frame.lookup_offset = TFX_RIBBON_EMITTER_ROW(buffer_index, row).lookup_offset;
	frame.angle_type = TFX_RIBBON_EMITTER_ROW(buffer_index, row).angle_type;
	frame.position = TFX_RIBBON_EMITTER_ROW(buffer_index, row).position;
	frame.captured_position = TFX_RIBBON_EMITTER_ROW(buffer_index, row).captured_position;
	frame.morph_segment_start_index = TFX_RIBBON_EMITTER_ROW(buffer_index, row).morph_segment_start_index;
	frame.noise_frequency = TFX_RIBBON_EMITTER_ROW(buffer_index, row).noise_frequency;
	frame.noise_speed = TFX_RIBBON_EMITTER_ROW(buffer_index, row).noise_speed;
	frame.noise_phase_range = TFX_RIBBON_EMITTER_ROW(buffer_index, row).noise_phase_range;
	frame.noise_lock_rate = TFX_RIBBON_EMITTER_ROW(buffer_index, row).noise_lock_rate;
	frame.fixed_angle_normal = TFX_RIBBON_EMITTER_ROW(buffer_index, row).fixed_angle_normal;
	frame.sample_count_noise = TFX_RIBBON_EMITTER_ROW(buffer_index, row).sample_count_noise;
	frame.lag_time = TFX_RIBBON_EMITTER_ROW(buffer_index, row).lag_time;
	frame.lag_span = TFX_RIBBON_EMITTER_ROW(buffer_index, row).lag_span;
	return frame;
}

//What the ribbon vertex shader needs of the emitter: .x = lookup_offset, .y = start_frame_index
TFX_BINDLESS_ACCESS uint2 tfx_load_ribbon_emitter_draw(uint buffer_index, uint row) {
	return uint2(TFX_RIBBON_EMITTER_ROW(buffer_index, row).lookup_offset, TFX_RIBBON_EMITTER_ROW(buffer_index, row).start_frame_index);
}

TFX_BINDLESS_ACCESS float3 tfx_load_lag_spine_position(uint buffer_index, uint row, uint sample_index) {
	return float3(TFX_RIBBON_EMITTER_ROW(buffer_index, row).lag_spine_position_x[sample_index],
	              TFX_RIBBON_EMITTER_ROW(buffer_index, row).lag_spine_position_y[sample_index],
	              TFX_RIBBON_EMITTER_ROW(buffer_index, row).lag_spine_position_z[sample_index]);
}

TFX_BINDLESS_ACCESS float4 tfx_load_lag_spine_quaternion(uint buffer_index, uint row, uint sample_index) {
	return tfx_unpack16bit_quaternion(TFX_RIBBON_EMITTER_ROW(buffer_index, row).lag_spine_quaternion[sample_index]);
}

//The bounds passes are editor only, so they have no fixed binding
TFX_BINDLESS_ACCESS void tfx_store_bounding_box(uint buffer_index, uint row, tfx_bounding_box box) {
	tfx_bounding_box_buffers[buffer_index][row] = box;
}

//--- Render_math

static const uint tfx_invalid_index = 0xFFFFFFFFu;
static const uint tfx_quad_indexes[6] = { 0, 1, 2, 2, 1, 3 };
//The epsilon keeps cross(up, alignment) from reaching zero for an alignment straight up or down, which is why the alignment frame needs no fallback
static const float3 tfx_billboard_up = float3(0, 1, 0.00001);
static const float tfx_size_unpack_scale = 64.0 / 32767.0;
static const float tfx_color_unpack_scale = 128.0 / 32767.0;
//tfxReverseAnimationFlag in the emitter properties flags
static const uint tfx_reverse_animation_flag = 8192;

//The r16g16_sscaled vertex format: two sign extended 16 bit halves
float2 tfx_unpack_sscaled16x2(uint packed) {
	return float2(float(int(packed << 16) >> 16), float(int(packed) >> 16));
}

//Two unsigned 16 bit halves
float2 tfx_unpack_uscaled16x2(uint packed) {
	return float2(float(packed & 0xFFFFu), float(packed >> 16));
}

//GLSL's unpackSnorm4x8, xyz of it
float3 tfx_unpack_snorm8x3(uint packed) {
	const int3 value = int3(int(packed << 24) >> 24, int(packed << 16) >> 24, int(packed << 8) >> 24);
	return max(float3(value) / 127.0, -1.0);
}

//GLSL's unpackUnorm4x8
float4 tfx_unpack_unorm8x4(uint packed) {
	return float4(float(packed & 0xFFu), float((packed >> 8) & 0xFFu), float((packed >> 16) & 0xFFu), float(packed >> 24)) / 255.0;
}

//GLSL's packUnorm4x8
uint tfx_pack_unorm8x4(float4 value) {
	const uint4 scaled = uint4(round(saturate(value) * 255.0));
	return scaled.x | (scaled.y << 8) | (scaled.z << 16) | (scaled.w << 24);
}

//GLSL's packSnorm4x8
uint tfx_pack_snorm8x4(float4 value) {
	const int4 scaled = int4(round(clamp(value, -1.0, 1.0) * 127.0));
	return (uint(scaled.x) & 0xFFu) | ((uint(scaled.y) & 0xFFu) << 8) | ((uint(scaled.z) & 0xFFu) << 16) | ((uint(scaled.w) & 0xFFu) << 24);
}

//GLSL's packSnorm2x16
uint tfx_pack_snorm16x2(float2 value) {
	const int2 scaled = int2(round(clamp(value, -1.0, 1.0) * 32767.0));
	return (uint(scaled.x) & 0xFFFFu) | (uint(scaled.y) << 16);
}

//GLSL's mat3(a, b, c) takes columns, Slang's float3x3(a, b, c) takes rows
float3x3 tfx_matrix_from_columns(float3 first, float3 second, float3 third) {
	return transpose(float3x3(first, second, third));
}

float3x3 tfx_view_rotation(float4x4 view) {
	return float3x3(view[0].xyz, view[1].xyz, view[2].xyz);
}

float3x3 tfx_quaternion_to_matrix(float4 q) {
	const float xx = q.x * q.x; const float xy = q.x * q.y; const float xz = q.x * q.z; const float xw = q.x * q.w;
	const float yy = q.y * q.y; const float yz = q.y * q.z; const float yw = q.y * q.w;
	const float zz = q.z * q.z; const float zw = q.z * q.w;
	return float3x3(
		1.0 - 2.0 * (yy + zz), 2.0 * (xy - zw), 2.0 * (xz + yw),
		2.0 * (xy + zw), 1.0 - 2.0 * (xx + zz), 2.0 * (yz - xw),
		2.0 * (xz - yw), 2.0 * (yz + xw), 1.0 - 2.0 * (xx + yy));
}

//The world frame of a quad aligned to the vector: tangent, bitangent and the alignment as its columns
float3x3 tfx_alignment_frame(float3 alignment) {
	const float3 normal = normalize(alignment);
	const float3 tangent = normalize(cross(tfx_billboard_up, normal));
	const float3 bitangent = normalize(cross(normal, tangent));
	return tfx_matrix_from_columns(tangent, bitangent, normal);
}

//The quad's rotation for its alignment mode, from the lowest 2 bits of the emitter properties flags:
//00 faces the camera, 01 doesn't, 10 faces the camera and rolls to the alignment, 11 aligns to the vector with no camera facing
float3x3 tfx_billboard_rotation(float4 quaternion, float3 alignment, float3 camera_relative_alignment, uint align_flags) {
	const float3x3 spin = tfx_quaternion_to_matrix(quaternion);
	const float3x3 align_matrix = tfx_alignment_frame(alignment);
	//The roll that turns the view's up, (0, 1), to the alignment projected onto the view plane
	const float2 align_xy = normalize(camera_relative_alignment.xy + 0.00001);
	const float roll_cos = align_xy.y;
	const float roll_sin = -align_xy.x;
	const float3x3 roll_matrix = tfx_matrix_from_columns(float3(roll_cos, roll_sin, 0.0), float3(-roll_sin, roll_cos, 0.0), float3(0.0, 0.0, 1.0));
	if (align_flags == 0x3u) {
		return mul(align_matrix, spin);
	} else if (align_flags == 0x2u) {
		return mul(roll_matrix, spin);
	}
	return spin;
}

//A quad corner before rotation, corners 0 to 3 being top left, top right, bottom left, bottom right
float3 tfx_billboard_corner(uint corner, float2 size, float2 handle) {
	const float2 corner_position = float2(float(corner & 1u), float(corner >> 1));
	return float3(size.x * (corner_position.x - handle.x), -size.y * (corner_position.y - handle.y), 0);
}

//--- Billboard_vertex

//tfx_instance_t through the pipeline's vertex formats, one per instance
struct tfx_billboard_input {
	[[vk::location(0)]] float4 position : POSITION;							//Stretch in w
	[[vk::location(1)]] float4 quaternion : QUATERNION;						//r16g16b16a16_snorm
	[[vk::location(2)]] float2 size : SIZE;									//r16g16_sscaled
	[[vk::location(3)]] float3 alignment : ALIGNMENT;						//r8g8b8_snorm
	[[vk::location(4)]] float2 intensity_gradient_map : INTENSITY;			//r16g16_sscaled
	[[vk::location(5)]] float3 curved_alpha_life : CURVED_ALPHA;			//r8g8b8_unorm
	[[vk::location(6)]] uint indexes : INDEXES;								//[gpu property index << 16 | capture flag << 15 | image frame]
	[[vk::location(7)]] uint captured_index : CAPTURED_INDEX;
};

//Explicit locations, so the editor's GLSL fragment shaders can take the same varyings. Only the uv is interpolated, the rest is per instance
struct tfx_billboard_output {
	float4 position : SV_Position;
	[[vk::location(0)]] float2 uv : TEXCOORD0;
	[[vk::location(1)]] nointerpolation uint2 indexes : TEXCOORD1;						//.x = [bindless image << 16 | image layer], .y = [color ramp layer << 8 | color ramp row]
	[[vk::location(2)]] nointerpolation float4 intensity_curved_alpha_map : TEXCOORD2;	//Intensity, curved alpha, sharpness, gradient map
};

//Everything but the position, which the live and baked billboards share
tfx_billboard_output tfx_billboard_varyings(float4 clip_position, uint corner, tfx_image_data image, tfx_emitter_properties properties, float2 intensity_gradient_map, float3 curved_alpha_life) {
	tfx_billboard_output output;
	output.position = clip_position;
	output.uv = float2((corner & 1u) != 0 ? image.uv.z : image.uv.x, (corner >> 1) != 0 ? image.uv.w : image.uv.y);
	output.indexes = uint2(image.texture_array_index, properties.color_ramp_indexes & 0xFFFFu);
	output.intensity_curved_alpha_map = float4(intensity_gradient_map.x * tfx_color_unpack_scale, curved_alpha_life.x, curved_alpha_life.y, intensity_gradient_map.y * tfx_color_unpack_scale);
	return output;
}

//The live billboard draw. Each instance is interpolated from its captured self in the previous buffer by the uniform buffer's timer_lerp
tfx_billboard_output tfx_live_billboard(tfx_billboard_push push, tfx_billboard_input input, uint corner) {
	const tfx_view_uniforms camera = tfx_load_view(push.uniform_index);
	const tfx_emitter_properties properties = tfx_load_emitter_properties(push.properties_index, (input.indexes >> 16) & 0xFFFFu);

	const float2 size = input.size * tfx_size_unpack_scale;
	const float4 quaternion = normalize(input.quaternion);

	//Not interpolated on a particle's first frame, when it loops back to the start of a path, or when too long has passed, and hidden that frame
	const float interpolate_is_active = float((input.indexes & tfx_capture_bit) >> 15);
	//The captured index means nothing without the capture bit, so row 0 keeps that read in range
	const uint previous_index = interpolate_is_active > 0.0 ? (input.captured_index & 0x0FFFFFFFu) + push.index_offset : 0;
	const tfx_instance previous = tfx_load_instance_row(push.prev_billboards_index, previous_index);

	const float2 lerped_size = lerp(tfx_unpack_sscaled16x2(previous.size) * tfx_size_unpack_scale, size, camera.timer_lerp) * interpolate_is_active;
	const float3 lerped_position = lerp(previous.position.xyz, input.position.xyz, camera.timer_lerp);
	float4 previous_quaternion = tfx_unpack16bit_quaternion(previous.quaternion);
	//q and -q are the same rotation, and lerping from one to the other passes through zero
	previous_quaternion = dot(previous_quaternion, quaternion) < 0.0 ? -previous_quaternion : previous_quaternion;
	const float4 lerped_quaternion = normalize(lerp(previous_quaternion, quaternion, camera.timer_lerp));
	float3 motion = input.position.xyz - previous.position.xyz;
	motion.z += 0.000001;

	const bool has_alignment = dot(input.alignment, input.alignment) > 0;
	const float travel_distance = length(motion);
	const float stretch_factor = input.position.w * interpolate_is_active;
	motion = normalize(motion);
	const float3 final_alignment = has_alignment ? input.alignment : motion;

	const uint align_flags = properties.flags & 0x3u;
	const bool billboarding = (align_flags & 0x1u) != 0;

	//GLSL's v * inverse(mat3(view)), which is the view rotation times v because a look at view is orthonormal
	const float3x3 view_rotation = tfx_view_rotation(camera.view);
	const float3x3 rotation = tfx_billboard_rotation(lerped_quaternion, final_alignment, mul(view_rotation, final_alignment), align_flags);

	float3 position = mul(rotation, tfx_billboard_corner(corner, lerped_size, properties.image_handle));

	//Stretch follows the direction of travel, and a stretch of 1 is the distance travelled in a 60th of a second (position.w has the update length taken out).
	//Each corner moves along it by its own projection over the quad's projected half extent, which keeps the silhouette continuous as
	//the direction sweeps past a corner where sign(dot) would snap, and keeps the total elongation at travel_distance
	const float3 stretch_direction = billboarding ? motion : normalize(mul(view_rotation, motion));
	float projected_half_extent = 0.0;
	for (uint quad_corner = 0; quad_corner < 4; ++quad_corner) {
		projected_half_extent = max(projected_half_extent, abs(dot(mul(rotation, tfx_billboard_corner(quad_corner, lerped_size, properties.image_handle)), stretch_direction)));
	}
	const float3 stretch_offset = stretch_direction * (dot(position, stretch_direction) / max(projected_half_extent, 0.00001)) * 0.5 * travel_distance * stretch_factor;
	position += (travel_distance > 0.00001 && stretch_factor > 0.0) ? stretch_offset : float3(0.0, 0.0, 0.0);

	//A camera facing quad's corners are offsets in view space, the rest are rotated by the view with the centre
	const float4 view_position = mul(camera.view, float4(lerped_position, 1.0));
	const float3 view_offset = billboarding ? mul(view_rotation, position) : position;
	const float4 clip_position = mul(camera.projection, float4(view_position.xyz + view_offset, view_position.w));

	//The sprite holds the animation frame only, the emitter's properties hold where its image starts
	const tfx_image_data image = tfx_load_image_data(push.image_data_index, properties.start_frame_index + (input.indexes & 0x00001FFFu));
	return tfx_billboard_varyings(clip_position, corner, image, properties, input.intensity_gradient_map, input.curved_alpha_life);
}

#if defined(TFX_FIXED_BINDINGS)
//What the vertex formats do to tfx_instance_t, for targets without sscaled formats that read the instance from a buffer
tfx_billboard_input tfx_unpack_billboard_input(tfx_instance instance) {
	tfx_billboard_input input;
	input.position = instance.position;
	input.quaternion = max(tfx_unpack16bit_quaternion(instance.quaternion), -1.0);
	input.size = tfx_unpack_sscaled16x2(instance.size);
	input.alignment = tfx_unpack_snorm8x3(instance.alignment);
	input.intensity_gradient_map = tfx_unpack_sscaled16x2(instance.intensity_gradient_map);
	input.curved_alpha_life = tfx_unpack_unorm8x4(instance.curved_alpha_life).xyz;
	input.indexes = instance.indexes;
	input.captured_index = instance.captured_index;
	return input;
}

[shader("vertex")]
tfx_billboard_output tfx_billboard_vertex(uint vertex_id : SV_VertexID, uint instance_id : SV_InstanceID) {
	return tfx_live_billboard(tfx_billboard_entry_push(), tfx_unpack_billboard_input(tfx_load_draw_instance(instance_id)), tfx_quad_indexes[vertex_id]);
}
#else
[shader("vertex")]
tfx_billboard_output tfx_billboard_vertex(tfx_billboard_input input, uint vertex_id : SV_VertexID) {
	return tfx_live_billboard(tfx_billboard_entry_push(), input, tfx_quad_indexes[vertex_id]);
}
#endif

//--- Billboard_baked_vertex

//The baked billboard draw. tfx_sprite_playback has already interpolated, resolved the image frame with its start frame, and baked
//the stretch into position.w, so nothing here reads a previous buffer
tfx_billboard_output tfx_baked_billboard(tfx_baked_billboard_push push, tfx_billboard_input input, uint corner) {
	const tfx_view_uniforms camera = tfx_load_view(push.uniform_index);
	const tfx_emitter_properties properties = tfx_load_emitter_properties(push.properties_index, (input.indexes >> 16) & 0xFFFFu);

	const uint align_flags = properties.flags & 0x3u;
	const bool billboarding = (align_flags & 0x1u) != 0;
	const float3x3 view_rotation = tfx_view_rotation(camera.view);
	const float3 camera_relative_alignment = mul(view_rotation, input.alignment);
	const float3x3 rotation = tfx_billboard_rotation(normalize(input.quaternion), input.alignment, camera_relative_alignment, align_flags);

	float3 position = mul(rotation, tfx_billboard_corner(corner, input.size * tfx_size_unpack_scale, properties.image_handle));

	//position.w is the baked stretch, half of it forwards along the alignment and half backwards
	const float3 stretch_direction = billboarding ? input.alignment : normalize(camera_relative_alignment);
	position += input.position.w > 0.00001 ? stretch_direction * sign(dot(position, stretch_direction)) * 0.5 * input.position.w : float3(0.0, 0.0, 0.0);

	const float4 view_position = mul(camera.view, float4(input.position.xyz, 1.0));
	const float3 view_offset = billboarding ? mul(view_rotation, position) : position;
	float4 clip_position = mul(camera.projection, float4(view_position.xyz + view_offset, view_position.w));
	clip_position.xy += float2(push.ndc_offset_x, push.ndc_offset_y) * clip_position.w;

	const tfx_image_data image = tfx_load_image_data(push.image_data_index, input.indexes & 0x00001FFFu);
	return tfx_billboard_varyings(clip_position, corner, image, properties, input.intensity_gradient_map, input.curved_alpha_life);
}

#if defined(TFX_FIXED_BINDINGS)
[shader("vertex")]
tfx_billboard_output tfx_billboard_baked_vertex(uint vertex_id : SV_VertexID, uint instance_id : SV_InstanceID) {
	return tfx_baked_billboard(tfx_baked_billboard_entry_push(), tfx_unpack_billboard_input(tfx_load_draw_instance(instance_id)), tfx_quad_indexes[vertex_id]);
}
#else
[shader("vertex")]
tfx_billboard_output tfx_billboard_baked_vertex(tfx_billboard_input input, uint vertex_id : SV_VertexID) {
	return tfx_baked_billboard(tfx_baked_billboard_entry_push(), input, tfx_quad_indexes[vertex_id]);
}
#endif

//--- Particle_shading

//A billboard or ribbon pixel: the image's red picks the color ramp column, and curved alpha fades the image's alpha. Writes premultiplied
//alpha, so blend with ONE, ONE_MINUS_SRC_ALPHA
float4 tfx_shade_particle(float3 tex_coord, int3 texture_indexes, float4 intensity_curved_alpha_map, uint image_index, uint color_ramp_texture_index, uint sampler_index) {
	const float4 texel = tfx_sample_particle_image(image_index, sampler_index, tex_coord);
	const int ramp_x = int(saturate(texel.r * intensity_curved_alpha_map.w) * 255);
	float4 ramp_texel = tfx_load_color_ramp(color_ramp_texture_index, int3(ramp_x, texture_indexes.x, texture_indexes.y)) * intensity_curved_alpha_map.x;
	ramp_texel.a = min(1.0, ramp_texel.a);
	const float curved_alpha = 1.0 - smoothstep(texel.a * intensity_curved_alpha_map.z, texel.a, 1.0 - intensity_curved_alpha_map.y);
	return float4(texel.rgb * ramp_texel.rgb * curved_alpha * texel.a, texel.a * ramp_texel.a * curved_alpha);
}

[shader("fragment")]
float4 tfx_billboard_fragment(tfx_billboard_output input) : SV_Target {
	const tfx_billboard_push push = tfx_billboard_entry_push();
	const float3 tex_coord = float3(input.uv, float(input.indexes.x & 0xFFFFu));
	const int3 texture_indexes = int3(int(input.indexes.y & 0xFFu), int(input.indexes.y >> 8), 0);
	return tfx_shade_particle(tex_coord, texture_indexes, input.intensity_curved_alpha_map, input.indexes.x >> 16, push.color_ramp_texture_index, push.sampler_index);
}

//Below this many pixels across, a ribbon starts dropping fragments because the rasteriser only covers pixels whose centre falls
//inside it, and a continuous ribbon breaks into dashes. Fading it out as it nears that width trades a tail that vanishes early for one that disintegrates
static const float tfx_min_ribbon_pixels = 1.5;

[shader("fragment")]
float4 tfx_ribbon_fragment(tfx_ribbon_output input) : SV_Target {
	const tfx_ribbon_bucket_globals push = tfx_ribbon_push_constants;
	const float4 color = tfx_shade_particle(input.tex_coord, input.texture_indexes, input.intensity_curved_alpha_map, input.image_index, push.color_ramp_texture_index, push.sampler_index);
	//fwidth is how far across the ribbon one pixel steps, so its reciprocal is the ribbon's width in pixels
	const float across_per_pixel = fwidth(input.ribbon_across);
	const float ribbon_pixels = across_per_pixel > 0.0 ? 1.0 / across_per_pixel : tfx_min_ribbon_pixels;
	return color * smoothstep(0.0, tfx_min_ribbon_pixels, ribbon_pixels);
}

//--- Ribbon_shape

//The ribbon graphs, noise and path sampling that the ribbon build and the ribbon bounds share

//Graph slots in an emitter's block of the ribbon graph lookups
static const uint tfx_ribbon_graph_intensity = 0;
static const uint tfx_ribbon_graph_alpha_sharpness = 1;
static const uint tfx_ribbon_graph_curved_alpha = 2;
static const uint tfx_ribbon_graph_gradient_map = 3;
static const uint tfx_ribbon_graph_width = 4;
static const uint tfx_ribbon_graph_scale = 5;
static const uint tfx_ribbon_graph_uv_offset_y = 6;
static const uint tfx_ribbon_graph_uv_scale_y = 7;
static const uint tfx_ribbon_graph_clip_offset = 8;
static const uint tfx_ribbon_graph_clip_size = 9;
static const uint tfx_ribbon_graph_morph_amount = 10;
static const uint tfx_ribbon_graph_noise_amount = 11;
static const uint tfx_ribbon_graph_noise_scroll = 12;
static const uint tfx_ribbon_graph_lag_amount = 13;
static const uint tfx_ribbon_graph_intensity_overlength = 14;
static const uint tfx_ribbon_graph_alpha_sharpness_overlength = 15;
static const uint tfx_ribbon_graph_curved_alpha_overlength = 16;
static const uint tfx_ribbon_graph_gradient_map_overlength = 17;
static const uint tfx_ribbon_graph_width_overlength = 18;
static const uint tfx_ribbon_graph_fixed_angle = 19;
static const uint tfx_ribbon_graph_morph_bias = 20;
static const uint tfx_ribbon_graph_noise_envelope = 21;
static const uint tfx_ribbon_graph_lag_profile = 22;

//tfx_ribbon_flag_bits
static const uint tfx_ribbon_active = 1 << 0;
static const uint tfx_ribbon_relative = 1 << 1;
static const uint tfx_ribbon_flags_shift = 24;
static const uint tfx_ribbon_start_index_mask = 0x00FFFFFFu;

static const uint tfx_angle_type_face_camera = 0;
static const uint tfx_angle_type_fixed_angle = 1;
static const uint tfx_angle_type_frenet_serret = 2;

//Hash axes of the ribbon noise phase offset, one per field axis. They must differ from each other and from every other axis, or the offsets correlate
static const uint tfx_ribbon_noise_phase_x_axis = 0x27d4eb2du;
static const uint tfx_ribbon_noise_phase_y_axis = 0x165667b1u;
static const uint tfx_ribbon_noise_phase_z_axis = 0x9e3779b1u;

//The noise algorithm in the low 4 bits of a ribbon's packed noise, the octaves in the next 4
static const uint tfx_ribbon_noise_value = 4u;

float tfx_ribbon_emitter_graph(uint graphs_index, uint lookup_offset, uint slot, float t) {
	return tfx_sample_graph(tfx_load_graph_row(graphs_index, lookup_offset + slot), t);
}

//Perlin's gradient select: one of the cube's edge vectors, picked by the low bits of the hash
float tfx_gradient_dot(uint hash, float3 distance_to_corner) {
	const uint h = hash & 15u;
	const float u = h < 8u ? distance_to_corner.x : distance_to_corner.y;
	const float v = h < 4u ? distance_to_corner.y : ((h == 12u || h == 14u) ? distance_to_corner.x : distance_to_corner.z);
	return ((h & 1u) == 0u ? u : -u) + ((h & 2u) == 0u ? v : -v);
}

static const int3 tfx_noise_primes = int3(1619, 31337, 6971);

//One scalar field of gradient noise. The lattice holds zero and the gradient is what's random, so the field crosses a cell boundary
//with a non zero slope. Value noise stalls at every boundary to stay continuous, which reads as morphing between fixed poses
float tfx_gradient_noise_field(float3 position, uint seed) {
	const float3 floor_position = floor(position);
	const float3 relative = position - floor_position;
	//Quintic fade: zero first and second derivative at the boundary, so the gradient stays continuous
	const float3 fade = relative * relative * relative * (relative * (relative * 6.0 - 15.0) + 10.0);

	const int3 lattice0 = int3(floor_position) * tfx_noise_primes;
	const int3 lattice1 = lattice0 + tfx_noise_primes;

	const float near_lower = lerp(tfx_gradient_dot(tfx_seedgen(uint(lattice0.x + lattice0.y) ^ uint(lattice0.z) ^ seed), relative),
	                              tfx_gradient_dot(tfx_seedgen(uint(lattice1.x + lattice0.y) ^ uint(lattice0.z) ^ seed), relative - float3(1.0, 0.0, 0.0)), fade.x);
	const float far_lower = lerp(tfx_gradient_dot(tfx_seedgen(uint(lattice0.x + lattice1.y) ^ uint(lattice0.z) ^ seed), relative - float3(0.0, 1.0, 0.0)),
	                             tfx_gradient_dot(tfx_seedgen(uint(lattice1.x + lattice1.y) ^ uint(lattice0.z) ^ seed), relative - float3(1.0, 1.0, 0.0)), fade.x);
	const float lower_slice = lerp(near_lower, far_lower, fade.y);

	const float near_upper = lerp(tfx_gradient_dot(tfx_seedgen(uint(lattice0.x + lattice0.y) ^ uint(lattice1.z) ^ seed), relative - float3(0.0, 0.0, 1.0)),
	                              tfx_gradient_dot(tfx_seedgen(uint(lattice1.x + lattice0.y) ^ uint(lattice1.z) ^ seed), relative - float3(1.0, 0.0, 1.0)), fade.x);
	const float far_upper = lerp(tfx_gradient_dot(tfx_seedgen(uint(lattice0.x + lattice1.y) ^ uint(lattice1.z) ^ seed), relative - float3(0.0, 1.0, 1.0)),
	                             tfx_gradient_dot(tfx_seedgen(uint(lattice1.x + lattice1.y) ^ uint(lattice1.z) ^ seed), relative - float3(1.0, 1.0, 1.0)), fade.x);
	const float upper_slice = lerp(near_upper, far_upper, fade.y);

	return lerp(lower_slice, upper_slice, fade.z);
}

//Measured so gradient noise lands at the value noise's rms, otherwise switching algorithm would read as an amplitude change
static const float tfx_gradient_noise_scale = 1.397;

//Three decorrelated fields, three times the work of value noise, which gets all three axes out of one hash
float3 tfx_gradient_noise_3d(float3 position) {
	return float3(tfx_gradient_noise_field(position, 0u),
	              tfx_gradient_noise_field(position, 0x9e3779b9u),
	              tfx_gradient_noise_field(position, 0x85ebca6bu)) * tfx_gradient_noise_scale;
}

//Octaves at doubling frequency and halving amplitude, normalised by the root of the summed squares: independent octaves add in
//variance, so that keeps the apparent amplitude constant as the octave count changes
float3 tfx_ribbon_noise_octaves(float3 position, uint algorithm, uint octaves) {
	float3 sum = float3(0.0, 0.0, 0.0);
	float amplitude = 1.0;
	float frequency = 1.0;
	float normalisation = 0.0;
	for (uint octave = 0u; octave < octaves; ++octave) {
		sum += amplitude * (algorithm == tfx_ribbon_noise_value ? tfx_value_noise_3d(position * frequency) : tfx_gradient_noise_3d(position * frequency));
		normalisation += amplitude * amplitude;
		amplitude *= 0.5;
		frequency *= 2.0;
	}
	return sum * rsqrt(max(normalisation, 1e-8));
}

//Where along the noise field a ribbon samples: its own phase offset, scrolled by the scroll graph and moved through time
float3 tfx_ribbon_noise_sample_offset(float age, uint phase_seed, float phase_range, float noise_speed, float noise_lock_rate, float noise_scroll) {
	float noise_time = age * noise_speed;
	if (noise_lock_rate > 0.0) {
		noise_time = floor(noise_time * noise_lock_rate) / noise_lock_rate;
	}
	const float3 phase_offset = float3(tfx_white_unit(phase_seed, tfx_ribbon_noise_phase_x_axis),
	                                   tfx_white_unit(phase_seed, tfx_ribbon_noise_phase_y_axis),
	                                   tfx_white_unit(phase_seed, tfx_ribbon_noise_phase_z_axis)) * phase_range;
	return phase_offset + float3(noise_scroll, noise_time, 0.0);
}

//With a clip size only the clip window from the clip offset is sampled, and all of the ribbon's segments are squeezed into it
float3 tfx_sample_ribbon_path(uint segments_index, int segment, uint segment_count, uint base_index, uint sample_count, float clip_offset, float clip_size) {
	const float segment_lerp = float(segment) / float(segment_count - 1);
	const float stored = clamp(clip_offset + segment_lerp * clip_size, 0.0, 1.0) * float(sample_count - 1);
	const uint first_sample = min(uint(stored), sample_count - 2);
	const float fraction = clamp(stored - float(first_sample), 0.0, 1.0);
	return lerp(tfx_load_ribbon_segment(segments_index, base_index + first_sample), tfx_load_ribbon_segment(segments_index, base_index + first_sample + 1), fraction);
}

//Path, then morph, then noise. Both paths are stored at the same density, so one sample pair addresses either. The bias is looked up
//per segment rather than reused from the centre one, otherwise a steep bias flattens the forward vector across the transition band
float3 tfx_sample_ribbon_shape(uint segments_index, int segment, uint segment_count, uint base_index, uint morph_base_index, uint sample_count,
                               float clip_offset, float clip_size, float morph_amount, tfx_gpu_graph_data morph_bias_graph, tfx_gpu_graph_data noise_envelope_graph,
                               float noise_amplitude, float noise_frequency, float3 noise_sample_offset, uint noise_packed) {
	float3 position = tfx_sample_ribbon_path(segments_index, segment, segment_count, base_index, sample_count, clip_offset, clip_size);
	const float segment_lerp = float(segment) / float(segment_count - 1);
	if (morph_base_index != tfx_invalid_index) {
		const float morph = clamp(morph_amount + tfx_sample_graph(morph_bias_graph, clamp(segment_lerp, 0.0, 1.0)), 0.0, 1.0);
		position = lerp(position, tfx_sample_ribbon_path(segments_index, segment, segment_count, morph_base_index, sample_count, clip_offset, clip_size), morph);
	}
	if (noise_amplitude == 0.0 || noise_packed == 0u) {
		return position;
	}
	const float envelope = tfx_sample_graph(noise_envelope_graph, clamp(segment_lerp, 0.0, 1.0));
	const float3 sample_position = noise_sample_offset + float3(segment_lerp * noise_frequency, 0.0, 0.0);
	return position + noise_amplitude * envelope * tfx_ribbon_noise_octaves(sample_position, noise_packed & 0xFu, (noise_packed >> 4) & 0xFu);
}

//--- Ribbon_build

//tfx_ribbon_compute builds the vertex and index buffers of a bucket's ribbons, one thread per segment. Every ribbon in a bucket has the
//same segment count, which is always a multiple of 32, so a warp never spans two ribbons and the angle type and easing branches stay uniform

//Measured 2.5x faster than 1024 groups
static const uint tfx_ribbon_group_size = 256;

struct tfx_lag_transform {
	float3 position;
	float4 quaternion;
};

//The emitter's transform depth along the ribbon back in time, from the lag spine of recent emitter transforms
tfx_lag_transform tfx_emitter_transform_at(uint emitters_index, uint emitter_index, float frame_lerp, float lag_time, float lag_span, float depth) {
	const float frame_length = lag_span - lag_time;
	const float distance_back = (1.0 - frame_lerp) * frame_length + depth * lag_time;
	const float spine_position = clamp(distance_back / lag_span, 0.0, 1.0) * float(tfx_lag_spine_samples - 1);
	const uint first_sample = min(uint(spine_position), uint(tfx_lag_spine_samples - 2));
	const float fraction = spine_position - float(first_sample);

	const float4 first_quaternion = tfx_load_lag_spine_quaternion(emitters_index, emitter_index, first_sample);
	float4 second_quaternion = tfx_load_lag_spine_quaternion(emitters_index, emitter_index, first_sample + 1);
	second_quaternion *= dot(first_quaternion, second_quaternion) < 0.0 ? -1.0 : 1.0;

	tfx_lag_transform result;
	result.position = lerp(tfx_load_lag_spine_position(emitters_index, emitter_index, first_sample), tfx_load_lag_spine_position(emitters_index, emitter_index, first_sample + 1), fraction);
	result.quaternion = normalize(lerp(first_quaternion, second_quaternion, fraction));
	return result;
}

float3 tfx_ribbon_world_position(float3 segment_position, float scale, float4 ribbon_quaternion, float3 ribbon_position, float4 emitter_quaternion, float3 emitter_position) {
	return tfx_rotate_vector(tfx_rotate_vector(segment_position * scale, ribbon_quaternion) + ribbon_position, emitter_quaternion) + emitter_position;
}

float3 tfx_safe_normalize(float3 v) {
	return v * rsqrt(dot(v, v) + 1e-16);
}

uint tfx_pack_segment_data(uint segment_index, uint ribbon_lerp, uint tessellation, uint tessellation_index, uint side) {
	return (tessellation << 29) | (tessellation_index << 26) | (side << 25) | (ribbon_lerp << 12) | segment_index;
}

void tfx_build_ribbon_segment(tfx_ribbon_bucket_globals push, uint global_index) {
	//The library counts no indexes for an untessellated ribbon, so there's nothing to build
	if (global_index >= push.segment_count * push.ribbon_count || push.tessellation == 0) return;

	const uint local_ribbon_index = global_index / push.segment_count;
	const uint global_ribbon_index = local_ribbon_index + push.ribbon_offset;
	const uint local_segment_index = global_index % push.segment_count;
	const float ribbon_lerp = float(local_segment_index) / float(push.segment_count - 1);

	const tfx_ribbon ribbon = tfx_load_ribbon(push.ribbons_index, global_ribbon_index);
	const float ribbon_age = float(ribbon.age_width & 0xFFFFu) * (1.0 / 65535.0);
	const float ribbon_width = tfx_half_to_float(ribbon.age_width >> 16);
	const float effect_scale = tfx_half_to_float(ribbon.scale & 0xFFFFu);
	const uint ribbon_flags = ribbon.start_index_flags >> tfx_ribbon_flags_shift;
	const uint emitter_index = (ribbon.emitter_index_phase_seed & 0xFFFFu) + push.emitter_offset;
	const tfx_ribbon_emitter_frame emitter = tfx_load_ribbon_emitter_frame(push.emitters_index, emitter_index);
	const uint lookup_offset = emitter.lookup_offset;

	//All of a bucket's ribbon segments are in one array, so the segment offset is where this bucket's start
	const uint base_segment_index = push.segment_offset + (ribbon.start_index_flags & tfx_ribbon_start_index_mask);
	//A unique id per segment for the packed vertex attribute
	const uint global_segment_index = base_segment_index + local_segment_index;

	const float clip_offset = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_clip_offset, ribbon_age);
	const float clip_size = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_clip_size, ribbon_age);
	const float morph_amount = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_morph_amount, ribbon_age);
	const tfx_gpu_graph_data morph_bias_graph = tfx_load_graph_row(push.graphs_index, lookup_offset + tfx_ribbon_graph_morph_bias);
	//The morph block lives in the same bucket as the static one, so it rebases with the same offset
	const uint morph_base_segment_index = emitter.morph_segment_start_index == tfx_invalid_index ? tfx_invalid_index : push.segment_offset + emitter.morph_segment_start_index;

	const tfx_gpu_graph_data noise_envelope_graph = tfx_load_graph_row(push.graphs_index, lookup_offset + tfx_ribbon_graph_noise_envelope);
	const float noise_amplitude = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_noise_amount, ribbon_age);
	const float noise_scroll = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_noise_scroll, ribbon_age);
	const float3 noise_sample_offset = tfx_ribbon_noise_sample_offset(ribbon_age, ribbon.emitter_index_phase_seed >> 16, emitter.noise_phase_range, emitter.noise_speed, emitter.noise_lock_rate, noise_scroll);

	const bool relative = (ribbon_flags & tfx_ribbon_relative) != 0;
	const float4 ribbon_quaternion = tfx_unpack16bit_quaternion(ribbon.quaternion);
	float4 emitter_quaternion = relative ? tfx_unpack16bit_quaternion(emitter.quaternion) : float4(0, 0, 0, 1.0);
	float3 emitter_position = relative ? lerp(emitter.captured_position, emitter.position, push.lerp) : float3(0.0, 0.0, 0.0);
	const bool has_next = local_segment_index < push.segment_count - 1;
	const int neighbour_segment_index = has_next ? int(local_segment_index) + 1 : int(local_segment_index) - 1;
	float4 neighbour_emitter_quaternion = emitter_quaternion;
	float3 neighbour_emitter_position = emitter_position;
	bool lagging = false;
	if (relative && emitter.lag_time > 0.0) {
		const float lag_amount = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_lag_amount, ribbon_age);
		if (lag_amount > 0.0) {
			const tfx_gpu_graph_data lag_profile_graph = tfx_load_graph_row(push.graphs_index, lookup_offset + tfx_ribbon_graph_lag_profile);
			const float neighbour_lerp = clamp(float(neighbour_segment_index) / float(push.segment_count - 1), 0.0, 1.0);
			const tfx_lag_transform centre_transform = tfx_emitter_transform_at(push.emitters_index, emitter_index, push.lerp, emitter.lag_time, emitter.lag_span, lag_amount * tfx_sample_graph(lag_profile_graph, ribbon_lerp));
			const tfx_lag_transform neighbour_transform = tfx_emitter_transform_at(push.emitters_index, emitter_index, push.lerp, emitter.lag_time, emitter.lag_span, lag_amount * tfx_sample_graph(lag_profile_graph, neighbour_lerp));
			emitter_quaternion = centre_transform.quaternion;
			emitter_position = centre_transform.position;
			neighbour_emitter_quaternion = neighbour_transform.quaternion;
			neighbour_emitter_position = neighbour_transform.position;
			lagging = true;
		}
	}

	const uint sample_count = emitter.sample_count_noise & 0xFFFFFFu;
	const uint noise_packed = emitter.sample_count_noise >> 24;
	const float3 current_segment_position = tfx_sample_ribbon_shape(push.ribbon_segments_index, int(local_segment_index), push.segment_count, base_segment_index, morph_base_segment_index, sample_count,
		clip_offset, clip_size, morph_amount, morph_bias_graph, noise_envelope_graph, noise_amplitude, emitter.noise_frequency, noise_sample_offset, noise_packed);
	//Towards the next segment, or from the previous one at the end
	const float3 neighbour_segment_position = tfx_sample_ribbon_shape(push.ribbon_segments_index, neighbour_segment_index, push.segment_count, base_segment_index, morph_base_segment_index, sample_count,
		clip_offset, clip_size, morph_amount, morph_bias_graph, noise_envelope_graph, noise_amplitude, emitter.noise_frequency, noise_sample_offset, noise_packed);

	const float scale = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_scale, ribbon_age) * effect_scale;
	const float3 current_position = tfx_ribbon_world_position(current_segment_position, scale, ribbon_quaternion, ribbon.position, emitter_quaternion, emitter_position);
	const float3 neighbour_position = tfx_ribbon_world_position(neighbour_segment_position, scale, ribbon_quaternion, ribbon.position, neighbour_emitter_quaternion, neighbour_emitter_position);
	const float3 world_forward = tfx_safe_normalize(has_next ? neighbour_position - current_position : current_position - neighbour_position);
	float3 right;

	if (emitter.angle_type == tfx_angle_type_face_camera) {
		const float3 to_camera = normalize(push.camera_position.xyz - current_position);
		right = tfx_safe_normalize(cross(world_forward, normalize(to_camera + world_forward)));
	} else if (emitter.angle_type == tfx_angle_type_fixed_angle) {
		const float fixed_angle = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_fixed_angle, ribbon_lerp);
		const float c = cos(fixed_angle);
		const float s = sin(fixed_angle);
		if (lagging) {
			const float3 up = tfx_rotate_vector(emitter.fixed_angle_normal, tfx_multiply_quaternions(ribbon_quaternion, emitter_quaternion));
			right = tfx_safe_normalize(cross(world_forward, up));
			right = tfx_safe_normalize(right * c + cross(world_forward, right) * s + world_forward * dot(world_forward, right) * (1.0 - c));
		} else {
			//Built in the ribbon's own space and rotated into the world after
			const float3 forward = tfx_safe_normalize(has_next ? neighbour_segment_position - current_segment_position : current_segment_position - neighbour_segment_position);
			right = tfx_safe_normalize(cross(forward, emitter.fixed_angle_normal));
			right = right * c + cross(forward, right) * s + forward * dot(forward, right) * (1.0 - c);
			right = tfx_safe_normalize(tfx_rotate_vector(right, tfx_multiply_quaternions(ribbon_quaternion, emitter_quaternion)));
		}
	} else {
		float3 up = float3(0.0, 1.0, 0.0);
		if (abs(dot(world_forward, up)) > 0.999) {
			up = float3(1.0, 0.0, 0.0);
		}
		right = tfx_safe_normalize(cross(tfx_safe_normalize(cross(world_forward, up)), world_forward));
	}

	const uint ribbon_lerp_packed = uint(ribbon_lerp * 2047);
	const float ribbon_is_active = float(ribbon_flags & tfx_ribbon_active);
	//The ribbon's effect scale, which is what the positions above are already scaled by
	const float width = ribbon_width * tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_width_overlength, ribbon_lerp) * effect_scale
		* tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_width, ribbon_age) * ribbon_is_active;

	tfx_ribbon_vertex_row vertex;
	vertex.uv_offset_scale = float2(tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_uv_offset_y, ribbon_age),
	                                tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_uv_scale_y, ribbon_age));
	vertex.ribbon_index = global_ribbon_index;
	vertex.clipped = 0;

	const uint vertices_per_segment = (push.tessellation + 1) * 2;
	const uint vertex_base_index = local_ribbon_index * push.segment_count * vertices_per_segment + local_segment_index * vertices_per_segment + push.vertex_offset;
	//One side's vertexes from the edge in to the centre line, then the other's from the centre line out
	for (uint tessellation_step = 0; tessellation_step <= push.tessellation; tessellation_step++) {
		const float t = float(tessellation_step) / float(push.tessellation);
		vertex.position = current_position + right * lerp(0.5, 0.0, t) * width;
		vertex.segment_index = tfx_pack_segment_data(global_segment_index, ribbon_lerp_packed, push.tessellation, tessellation_step, 0);
		tfx_store_ribbon_vertex(push.vertexes_index, vertex_base_index + tessellation_step, vertex);

		vertex.position = current_position + right * lerp(0.0, -0.5, t) * width;
		vertex.segment_index = tfx_pack_segment_data(global_segment_index, ribbon_lerp_packed, push.tessellation, tessellation_step, 1);
		tfx_store_ribbon_vertex(push.vertexes_index, vertex_base_index + push.tessellation + 1 + tessellation_step, vertex);
	}

	if (has_next) {
		const uint indexes_per_segment = push.tessellation * 12;
		const uint indexes_per_ribbon = push.segment_count * indexes_per_segment - indexes_per_segment;
		const uint index_base = push.index_offset + local_ribbon_index * indexes_per_ribbon + local_segment_index * indexes_per_segment;
		const uint next_base_vertex = vertex_base_index + vertices_per_segment;
		for (uint quad = 0; quad < push.tessellation; quad++) {
			const uint quad_index = index_base + quad * 12;
			const uint4 first_side = uint4(vertex_base_index + quad, vertex_base_index + quad + 1, next_base_vertex + quad, next_base_vertex + quad + 1);
			const uint4 second_side = first_side + push.tessellation + 1;
			tfx_store_ribbon_triangle(push.indexes_index, quad_index, first_side.xyz);
			tfx_store_ribbon_triangle(push.indexes_index, quad_index + 3, first_side.zyw);
			tfx_store_ribbon_triangle(push.indexes_index, quad_index + 6, second_side.xyz);
			tfx_store_ribbon_triangle(push.indexes_index, quad_index + 9, second_side.zyw);
		}
	}
}

[shader("compute")]
[numthreads(tfx_ribbon_group_size, 1, 1)]
void tfx_ribbon_compute(uint3 thread_id : SV_DispatchThreadID) {
	tfx_build_ribbon_segment(tfx_ribbon_push_constants, thread_id.x);
}

//--- Ribbon_draw

//tfx_ribbon_vertex_t through the pipeline's vertex formats
struct tfx_ribbon_input {
	[[vk::location(0)]] float3 position : POSITION;
	[[vk::location(1)]] uint segment_index : SEGMENT_INDEX;
	[[vk::location(2)]] float2 uv_offset_scale : UV_OFFSET_SCALE;
	[[vk::location(3)]] uint ribbon_index : RIBBON_INDEX;
};

struct tfx_ribbon_output {
	float4 position : SV_Position;
	[[vk::location(0)]] float3 tex_coord : TEXCOORD0;
	[[vk::location(1)]] nointerpolation int3 texture_indexes : TEXCOORD1;
	[[vk::location(2)]] float4 intensity_curved_alpha_map : TEXCOORD2;
	[[vk::location(3)]] nointerpolation uint image_index : TEXCOORD3;
	//0 to 1 across the ribbon's width, so its screen space derivative gives the ribbon's width in pixels
	[[vk::location(4)]] float ribbon_across : TEXCOORD4;
};

tfx_ribbon_output tfx_ribbon_draw(tfx_ribbon_bucket_globals push, tfx_ribbon_input input) {
	const tfx_view_uniforms camera = tfx_load_view(push.uniform_index);
	tfx_ribbon_output output;
	output.position = mul(camera.projection, mul(camera.view, float4(input.position, 1.0)));
	output.position.xy += float2(push.ndc_offset_x, push.ndc_offset_y) * output.position.w;

	const float ribbon_position = float((input.segment_index & 0x007FF000u) >> 12) / 2047.0;
	const tfx_ribbon ribbon = tfx_load_ribbon(push.ribbons_index, input.ribbon_index);
	const uint2 emitter = tfx_load_ribbon_emitter_draw(push.emitters_index, (ribbon.emitter_index_phase_seed & 0xFFFFu) + push.emitter_offset);
	const uint lookup_offset = emitter.x;
	//The ribbon holds the animation frame, its emitter holds where the image starts
	const tfx_image_data image = tfx_load_image_data(push.image_data_index, emitter.y + (ribbon.texture_indexes & 0x00001FFFu));

	const float tessellation = float((input.segment_index & 0xE0000000u) >> 29);
	const float tessellation_index = float((input.segment_index & 0x1C000000u) >> 26);
	const float side = float((input.segment_index & 0x02000000u) >> 25);

	const float sampled_gradient_map = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_gradient_map_overlength, ribbon_position);
	const float sampled_intensity = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_intensity_overlength, ribbon_position);
	const float sampled_curved_alpha = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_curved_alpha_overlength, ribbon_position);
	const float sampled_alpha_sharpness = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_alpha_sharpness_overlength, ribbon_position);

	//The images are sprite sheets in a texture array, so a texture repeating along the ribbon ping pongs: every other wrap runs backwards,
	//which leaves no flipped one segment seam where it wraps. The uv rect is the whole of the shape's own image, so a repeat sampler could replace this
	const float scaled_position = ribbon_position * input.uv_offset_scale.y + input.uv_offset_scale.x;
	const float wrap_index = floor(scaled_position);
	float wrap_fraction = frac(scaled_position);
	//GLSL's mod, which unlike fmod is never negative
	if (wrap_index - 2.0 * floor(wrap_index * 0.5) == 1.0) {
		wrap_fraction = 1.0 - wrap_fraction;
	}

	const float ribbon_across = (tessellation_index / tessellation) * 0.5 + side * 0.5;
	output.ribbon_across = ribbon_across;
	const float2 ribbon_intensity_gradient_map = tfx_unpack_sscaled16x2(ribbon.intensity_gradient_map) * tfx_color_unpack_scale;
	const float2 ribbon_curved_alpha = tfx_unpack_sscaled16x2(ribbon.curved_alpha) * tfx_color_unpack_scale;

	output.image_index = image.texture_array_index >> 16;
	output.tex_coord = float3(wrap_fraction * (image.uv.z - image.uv.x) + image.uv.x, ribbon_across * (image.uv.w - image.uv.y) + image.uv.y, float(image.texture_array_index & 0xFFFFu));
	output.texture_indexes = int3(int((ribbon.texture_indexes & 0xFF000000u) >> 24), int((ribbon.texture_indexes & 0x00FF0000u) >> 16), int(ribbon_position * 255));
	output.intensity_curved_alpha_map = float4(sampled_intensity * ribbon_intensity_gradient_map.x, sampled_curved_alpha * ribbon_curved_alpha.x,
		sampled_alpha_sharpness * ribbon_curved_alpha.y, sampled_gradient_map * ribbon_intensity_gradient_map.y);
	return output;
}

#if defined(TFX_FIXED_BINDINGS)
//The vertexes tfx_ribbon_compute built, read from the buffer it wrote
[shader("vertex")]
tfx_ribbon_output tfx_ribbon_vertex(uint vertex_id : SV_VertexID) {
	const tfx_ribbon_vertex_row vertex = tfx_load_draw_ribbon_vertex(vertex_id);
	tfx_ribbon_input input;
	input.position = vertex.position;
	input.segment_index = vertex.segment_index;
	input.uv_offset_scale = vertex.uv_offset_scale;
	input.ribbon_index = vertex.ribbon_index;
	return tfx_ribbon_draw(tfx_ribbon_push_constants, input);
}
#else
[shader("vertex")]
tfx_ribbon_output tfx_ribbon_vertex(tfx_ribbon_input input) {
	return tfx_ribbon_draw(tfx_ribbon_push_constants, input);
}
#endif

//--- Sprite_playback

//The playback passes turn baked frames into live instances, one thread per sprite or ribbon in the frame, interpolating each from its
//captured self in the frame before

static const uint tfx_playback_group_size = 128;

//The offsets are a running total per animation instance, so the instance a thread is in is the count of totals at or below it
uint tfx_playback_instance_index(uint offsets_index, uint animation_instances_total, uint index) {
	uint instance_index = 0;
	for (uint instance = 0; instance < animation_instances_total; ++instance) {
		instance_index += uint(tfx_load_offset(offsets_index, instance_index) <= index);
	}
	return instance_index;
}

//How far between this baked frame and the next the animation is, 0 to 1
float tfx_playback_frame_lerp(tfx_animation_instance animation) {
	const float frame_time = (animation.current_time / animation.animation_length_in_time) * float(animation.frame_count);
	return frame_time - floor(frame_time);
}

[shader("compute")]
[numthreads(tfx_playback_group_size, 1, 1)]
void tfx_sprite_playback(uint3 thread_id : SV_DispatchThreadID) {
	const tfx_sprite_data_push push = tfx_sprite_data_push_constants;
	const uint index = thread_id.x;
	if (index >= push.billboards_total) return;

	const uint instance_index = tfx_playback_instance_index(push.offsets_index, push.animation_instances_total, index);
	const uint sprite_index_offset = instance_index == 0 ? 0 : tfx_load_offset(push.offsets_index, instance_index - 1);
	const tfx_animation_instance animation = tfx_load_animation_instance(push.animation_instances_index, instance_index);
	float lerp_amount = tfx_playback_frame_lerp(animation);

	const uint sprite_data_index = (index - sprite_index_offset) + animation.offset_into_sprite_data;
	const tfx_instance current = tfx_load_instance_row(push.sprite_data_index, sprite_data_index);
	//A sprite in its first frame interpolates with itself
	const uint captured_index = current.captured_index == tfx_invalid_index ? sprite_data_index : current.captured_index;
	const tfx_instance captured = tfx_load_instance_row(push.sprite_data_index, captured_index);

	const uint property_index = (current.indexes & 0xFFFF0000u) >> 16;
	const float4 curved_alpha_life = tfx_unpack_unorm8x4(current.curved_alpha_life);
	//lerp_offset is how far into the frame the sprite spawns, held on the spawn record, which is the captured sprite of the frame after.
	//Until then it has no size
	const float lerp_offset = tfx_unpack_unorm8x4(captured.curved_alpha_life).w;
	const float spawn_scale = (lerp_offset >= 1.0 - lerp_amount && captured_index != sprite_data_index) ? 1.0 : 0.0;
	//The captured state is the sprite at its spawn time, so interpolate from the spawn point rather than the start of the frame
	lerp_amount = saturate((lerp_amount - (1.0 - lerp_offset)) / max(lerp_offset, 1.0 / 255.0));

	const float3 lerped_position = lerp(captured.position.xyz, current.position.xyz, lerp_amount) * animation.scale + animation.position;
	const float4 lerped_quaternion = lerp(tfx_unpack16bit_quaternion(captured.quaternion), tfx_unpack16bit_quaternion(current.quaternion), lerp_amount);
	const float2 lerped_scale = lerp(tfx_unpack_uscaled16x2(captured.size), tfx_unpack_uscaled16x2(current.size), lerp_amount) * tfx_size_unpack_scale * spawn_scale * animation.scale;
	const float3 lerped_alignment = lerp(tfx_unpack_snorm8x3(captured.alignment), tfx_unpack_snorm8x3(current.alignment), lerp_amount);
	const float lerped_stretch = lerp(captured.position.w, current.position.w, lerp_amount) * spawn_scale;
	const float2 lerped_intensity_gradient_map = lerp(tfx_unpack_uscaled16x2(captured.intensity_gradient_map), tfx_unpack_uscaled16x2(current.intensity_gradient_map), lerp_amount) * tfx_color_unpack_scale;
	const float3 lerped_curved_alpha_life = lerp(tfx_unpack_unorm8x4(captured.curved_alpha_life).xyz, curved_alpha_life.xyz, lerp_amount);

	const tfx_emitter_properties properties = tfx_load_emitter_properties(push.emitter_properties_index, property_index);
	uint lerped_frame = properties.start_frame_index;
	uint indexes;
	if (push.animated_shapes != 0) {
		//The same for every thread, so it never diverges. Interpolates the shape's frame, allowing for it wrapping round at the end
		const float reverse_direction = float((properties.flags & tfx_reverse_animation_flag) > 1);
		const float current_frame = float(current.indexes & 0x00001FFFu);
		const float previous_frame = float(captured.indexes & 0x00001FFFu);
		const float frames = properties.animation_frames;
		const float wrap = (1.0 - reverse_direction) * float(current_frame < previous_frame);
		const float reverse_wrap = reverse_direction * float(previous_frame < current_frame);
		const float frame = lerp(previous_frame + frames * reverse_wrap, current_frame + frames * wrap, lerp_amount);
		//GLSL's mod
		lerped_frame += uint(frame - frames * floor(frame / frames));
		indexes = (current.indexes & 0xFFFFE000u) | lerped_frame;
	} else {
		indexes = current.indexes + lerped_frame;
	}

	tfx_instance sprite;
	sprite.position = float4(lerped_position, lerped_stretch);
	sprite.quaternion = uint2(tfx_pack_snorm16x2(lerped_quaternion.xy), tfx_pack_snorm16x2(lerped_quaternion.zw));
	sprite.size = (uint(lerped_scale.y * tfx_size_packed_scale) << 16) + uint(lerped_scale.x * tfx_size_packed_scale);
	sprite.alignment = tfx_pack_snorm8x4(float4(lerped_alignment, 0));
	sprite.intensity_gradient_map = (uint(lerped_intensity_gradient_map.y * tfx_color_packed_scale) << 16) + uint(lerped_intensity_gradient_map.x * tfx_color_packed_scale);
	sprite.curved_alpha_life = tfx_pack_unorm8x4(float4(lerped_curved_alpha_life, 0));
	sprite.indexes = indexes;
	sprite.captured_index = index;
	tfx_store_instance_row(push.billboards_index, index, sprite);
}

//--- Ribbon_playback

[shader("compute")]
[numthreads(tfx_playback_group_size, 1, 1)]
void tfx_ribbon_playback(uint3 thread_id : SV_DispatchThreadID) {
	const tfx_ribbon_data_push push = tfx_ribbon_data_push_constants;
	const uint index = thread_id.x;
	if (index >= push.ribbons_total) return;

	const uint instance_index = tfx_playback_instance_index(push.offsets_index, push.animation_instances_total, index);
	const uint ribbon_index_offset = instance_index == 0 ? 0 : tfx_load_offset(push.offsets_index, instance_index - 1);
	const tfx_animation_instance animation = tfx_load_animation_instance(push.animation_instances_index, instance_index);
	float lerp_amount = tfx_playback_frame_lerp(animation);

	const uint ribbon_record_index = (index - ribbon_index_offset) + animation.offset_into_ribbon_data;
	const tfx_ribbon_instance_data current = tfx_load_baked_ribbon(push.ribbon_data_index, ribbon_record_index);
	//A ribbon in its first frame interpolates with itself
	const uint captured_index = current.captured_index == tfx_invalid_index ? ribbon_record_index : current.captured_index;
	const tfx_ribbon_instance_data captured = tfx_load_baked_ribbon(push.ribbon_data_index, captured_index);

	//The spawn record holds the lerp_offset, and it's the captured ribbon of the frame after, the first one that interpolates from it
	const float lerp_offset = float(captured.additional & 0xFFFFu) * (1.0 / 65535.0);
	const float spawn_scale = (lerp_offset >= 1.0 - lerp_amount && captured_index != ribbon_record_index) ? 1.0 : 0.0;
	//The captured state is the ribbon at its spawn time, so interpolate from the spawn point rather than the start of the frame
	lerp_amount = saturate((lerp_amount - (1.0 - lerp_offset)) / max(lerp_offset, 1.0 / 65535.0));

	const float lerped_age = lerp(captured.position.w, current.position.w, lerp_amount);
	//No width until it's spawned
	const float lerped_width = lerp(captured.width, current.width, lerp_amount) * spawn_scale;
	const float4 lerped_quaternion = lerp(tfx_unpack16bit_quaternion(captured.quaternion), tfx_unpack16bit_quaternion(current.quaternion), lerp_amount);
	const float2 lerped_intensity_gradient_map = lerp(tfx_unpack_uscaled16x2(captured.intensity_gradient_map), tfx_unpack_uscaled16x2(current.intensity_gradient_map), lerp_amount) * tfx_color_unpack_scale;
	//Packed as live ribbons pack it, sscaled, and lerped in whole steps so a ribbon that doesn't change comes out bit for bit
	const float2 lerped_curved_alpha = round(lerp(tfx_unpack_sscaled16x2(captured.curved_alpha), tfx_unpack_sscaled16x2(current.curved_alpha), lerp_amount));
	const uint property_index = (current.additional & 0xFFFF0000u) >> 16;

	tfx_ribbon ribbon;
	ribbon.position = lerp(captured.position.xyz, current.position.xyz, lerp_amount) * animation.scale + animation.position;
	ribbon.age_width = uint(saturate(lerped_age) * 65535.0 + 0.5) | (tfx_float_to_half(lerped_width) << 16);
	ribbon.quaternion = uint2(tfx_pack_snorm16x2(lerped_quaternion.xy), tfx_pack_snorm16x2(lerped_quaternion.zw));
	ribbon.start_index_flags = (current.start_index & tfx_ribbon_start_index_mask) | (current.flags << tfx_ribbon_flags_shift);
	//The baked emitter buffer has one row per ribbon property, so the property index is the emitter index. The seed is a per ribbon constant
	ribbon.emitter_index_phase_seed = (property_index & 0xFFFFu) | (current.phase_seed << 16);
	ribbon.texture_indexes = (current.texture_indexes & 0xFFFFE000u) | tfx_load_ribbon_properties(push.ribbon_properties_index, property_index).start_frame_index;
	ribbon.intensity_gradient_map = (uint(lerped_intensity_gradient_map.y * tfx_color_packed_scale) << 16) + uint(lerped_intensity_gradient_map.x * tfx_color_packed_scale);
	ribbon.curved_alpha = (uint(int(lerped_curved_alpha.x)) & 0xFFFFu) | (uint(int(lerped_curved_alpha.y)) << 16);
	ribbon.scale = tfx_float_to_half(lerp(captured.scale, current.scale, lerp_amount));
	tfx_store_ribbon(push.ribbon_output_index, index, ribbon);
}

//--- Bounds_reduction

//The bounding box of a baked frame, for export and the auto fit camera. Each thread bounds two sprites or ribbons, half a group apart,
//so a dispatch is half as many groups as there are group sized blocks, then the group reduces to one box and the CPU reduces the groups' boxes

static const uint tfx_bounds_group_size = 128;

groupshared float3 tfx_bounds_min[tfx_bounds_group_size];
groupshared float3 tfx_bounds_max[tfx_bounds_group_size];

//A barrier at every step: warps don't run in lockstep on NVIDIA since Volta, and without one a compiler is free to keep the shared values in registers
void tfx_reduce_bounds(uint local_index, uint group_index, float3 min_vertex, float3 max_vertex, uint bounding_boxes_index) {
	tfx_bounds_min[local_index] = min_vertex;
	tfx_bounds_max[local_index] = max_vertex;
	GroupMemoryBarrierWithGroupSync();
	[unroll]
	for (uint stride = tfx_bounds_group_size / 2; stride > 0; stride >>= 1) {
		if (local_index < stride) {
			tfx_bounds_min[local_index] = min(tfx_bounds_min[local_index], tfx_bounds_min[local_index + stride]);
			tfx_bounds_max[local_index] = max(tfx_bounds_max[local_index], tfx_bounds_max[local_index + stride]);
		}
		GroupMemoryBarrierWithGroupSync();
	}
	if (local_index == 0) {
		tfx_bounding_box box;
		box.min_corner = tfx_bounds_min[0];
		box.max_stretch = 0;
		box.max_corner = tfx_bounds_max[0];
		box.padding = 0;
		tfx_store_bounding_box(bounding_boxes_index, group_index, box);
	}
}

//Where the baked billboard draw puts a sprite's corners. A camera facing sprite (alignment 00 or 10) depends on a camera the bounds don't
//have, so it gets a conservative sphere instead
void tfx_add_sprite_bounds(tfx_sprite_data_push push, uint sprite_data_index, inout float3 min_vertex, inout float3 max_vertex) {
	const tfx_instance sprite = tfx_load_instance_row(push.sprite_data_index, sprite_data_index);
	const tfx_emitter_properties properties = tfx_load_emitter_properties(push.emitter_properties_index, (sprite.indexes & 0xFFFF0000u) >> 16);
	const uint align_flags = properties.flags & 0x3u;
	const float2 size = tfx_unpack_uscaled16x2(sprite.size) * tfx_size_unpack_scale;
	const float2 handle = properties.image_handle;
	const float stretch = sprite.position.w;
	const float3 alignment = tfx_unpack_snorm8x3(sprite.alignment);

	if ((align_flags & 0x1u) == 0 || align_flags == 0x2u) {
		//From the quad's centre to its furthest corner, then the handle's offset from the centre, then the stretch on every axis
		const float radius = length(tfx_billboard_corner(0, size, handle) - tfx_billboard_corner(3, size, handle)) * 0.5 + length((float2(0.5, 0.5) - handle) * size) + abs(stretch) * 0.5;
		min_vertex = min(min_vertex, sprite.position.xyz - radius);
		max_vertex = max(max_vertex, sprite.position.xyz + radius);
	} else {
		const float3x3 spin = tfx_quaternion_to_matrix(normalize(tfx_unpack16bit_quaternion(sprite.quaternion)));
		const float3x3 rotation = align_flags == 0x3u ? mul(tfx_alignment_frame(alignment), spin) : spin;
		for (uint corner = 0; corner < 4; ++corner) {
			float3 vertex = mul(rotation, tfx_billboard_corner(corner, size, handle));
			vertex += alignment * sign(dot(vertex, alignment)) * 0.5 * abs(stretch) + sprite.position.xyz;
			min_vertex = min(min_vertex, vertex);
			max_vertex = max(max_vertex, vertex);
		}
	}
}

[shader("compute")]
[numthreads(tfx_bounds_group_size, 1, 1)]
void tfx_sprite_bounds(uint3 group_thread_id : SV_GroupThreadID, uint3 group_id : SV_GroupID) {
	const tfx_sprite_data_push push = tfx_sprite_data_push_constants;
	const uint local_index = group_thread_id.x;
	const uint frame_offset = tfx_load_animation_instance(push.animation_instances_index, 0).offset_into_sprite_data;
	const uint last_sprite_in_frame = push.billboards_total > 0 ? push.billboards_total - 1 : 0;
	const uint sprite_in_frame = group_id.x * tfx_bounds_group_size * 2 + local_index;

	//Clamped rather than returning early, which would leave threads short of the barriers
	float3 min_vertex = float3(99999999.0, 99999999.0, 99999999.0);
	float3 max_vertex = float3(-99999999.0, -99999999.0, -99999999.0);
	tfx_add_sprite_bounds(push, frame_offset + min(sprite_in_frame, last_sprite_in_frame), min_vertex, max_vertex);
	tfx_add_sprite_bounds(push, frame_offset + min(sprite_in_frame + tfx_bounds_group_size, last_sprite_in_frame), min_vertex, max_vertex);
	tfx_reduce_bounds(local_index, group_id.x, min_vertex, max_vertex, push.bounding_boxes_index);
}

//Every drawn segment position, built as tfx_ribbon_compute builds it, expanded by the half width. Only the right vector needs the camera
//or the angle mode, so the half width is the one term spread over every axis. The emitter transform is already folded into the record
void tfx_add_ribbon_bounds(tfx_ribbon_bounds_push push, uint ribbon_index, inout float3 min_vertex, inout float3 max_vertex) {
	const tfx_ribbon_instance_data ribbon = tfx_load_baked_ribbon(push.ribbon_data_index, ribbon_index);
	const tfx_animation_ribbon_properties properties = tfx_load_ribbon_properties(push.ribbon_properties_index, (ribbon.additional >> 16) & 0xFFFFu);
	if (properties.segment_count < 2 || properties.sample_count < 2 || ribbon.start_index == tfx_invalid_index) {
		return;
	}
	const uint lookup_offset = properties.graph_lookup_offset;
	const float age = ribbon.position.w;
	const float clip_offset = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_clip_offset, age);
	const float clip_size = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_clip_size, age);
	const float scale = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_scale, age) * ribbon.scale;
	const float morph_amount = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_morph_amount, age);
	const float width_overtime = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_width, age);
	const float ribbon_is_active = float(ribbon.flags & tfx_ribbon_active);
	const float4 quaternion = tfx_unpack16bit_quaternion(ribbon.quaternion);
	const tfx_gpu_graph_data morph_bias_graph = tfx_load_graph_row(push.graphs_index, lookup_offset + tfx_ribbon_graph_morph_bias);
	const tfx_gpu_graph_data noise_envelope_graph = tfx_load_graph_row(push.graphs_index, lookup_offset + tfx_ribbon_graph_noise_envelope);
	const float noise_amplitude = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_noise_amount, age);
	const float noise_scroll = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_noise_scroll, age);
	const float3 noise_sample_offset = tfx_ribbon_noise_sample_offset(age, ribbon.phase_seed, properties.noise_phase_range, properties.noise_speed, properties.noise_lock_rate, noise_scroll);

	for (uint segment = 0; segment < properties.segment_count; ++segment) {
		const float3 segment_position = tfx_sample_ribbon_shape(push.ribbon_segments_index, int(segment), properties.segment_count, ribbon.start_index, properties.morph_segment_start_index, properties.sample_count,
			clip_offset, clip_size, morph_amount, morph_bias_graph, noise_envelope_graph, noise_amplitude, properties.noise_frequency, noise_sample_offset, properties.noise_packed);
		const float3 world_position = tfx_rotate_vector(segment_position * scale, quaternion) + ribbon.position.xyz;
		const float width_overlength = tfx_ribbon_emitter_graph(push.graphs_index, lookup_offset, tfx_ribbon_graph_width_overlength, float(segment) / float(properties.segment_count - 1));
		//A width graph is free to go negative, the quad just flips, so the half extent is the magnitude
		const float half_width = abs(ribbon.width * width_overlength * width_overtime) * ribbon_is_active * 0.5;
		min_vertex = min(min_vertex, world_position - half_width);
		max_vertex = max(max_vertex, world_position + half_width);
	}
}

[shader("compute")]
[numthreads(tfx_bounds_group_size, 1, 1)]
void tfx_ribbon_bounds(uint3 group_thread_id : SV_GroupThreadID, uint3 group_id : SV_GroupID) {
	const tfx_ribbon_bounds_push push = tfx_ribbon_bounds_push_constants;
	const uint local_index = group_thread_id.x;
	const uint frame_offset = tfx_load_animation_instance(push.animation_instances_index, 0).offset_into_ribbon_data;
	const uint last_ribbon_in_frame = push.ribbons_total > 0 ? push.ribbons_total - 1 : 0;
	const uint ribbon_in_frame = group_id.x * tfx_bounds_group_size * 2 + local_index;

	float3 min_vertex = float3(99999999.0, 99999999.0, 99999999.0);
	float3 max_vertex = float3(-99999999.0, -99999999.0, -99999999.0);
	tfx_add_ribbon_bounds(push, frame_offset + min(ribbon_in_frame, last_ribbon_in_frame), min_vertex, max_vertex);
	tfx_add_ribbon_bounds(push, frame_offset + min(ribbon_in_frame + tfx_bounds_group_size, last_ribbon_in_frame), min_vertex, max_vertex);
	tfx_reduce_bounds(local_index, group_id.x, min_vertex, max_vertex, push.bounding_boxes_index);
}
