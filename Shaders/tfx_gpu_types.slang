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
