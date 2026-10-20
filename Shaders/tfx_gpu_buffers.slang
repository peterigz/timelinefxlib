//Every buffer access the GPU particle kernel makes, so a port to another binding model replaces only this module
//Rows are relative to their own array and each accessor adds that array's base from the push constants, as tfx_GetGPUParticleFrame lays them out

__exported import tfx_gpu_types;
import tfx_gpu_spawn;
import tfx_graphs;

//-D TFX_FIXED_BINDINGS gives every target the four fixed raw buffers WGSL always uses, for D3D11, OpenGL and anything else without descriptor arrays
#if defined(TFX_FIXED_BINDINGS)
#define TFX_FIXED_CASES case wgsl: case hlsl: case glsl: case spirv:
#define TFX_FIXED_REGISTER(slot) : register(slot)
#define TFX_BINDLESS_ACCESS
#else
#define TFX_FIXED_CASES case wgsl:
#define TFX_FIXED_REGISTER(slot)
//GLSL needs the extension to index runtime descriptor arrays, and Slang only emits it from a function, never from the declarations
#define TFX_BINDLESS_ACCESS __glsl_extension(GL_EXT_nonuniform_qualifier)
#endif

//The pool's hot, cold and noise sections can be views of one buffer, as can the frame's arrays. Sharing a binding is the point, so its overlap warning is off
#pragma warning(disable: 39001)
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

//Metal has no descriptor arrays: there each push constant buffer index is a slot in this argument buffer of gpuAddress values
#ifndef TFX_METAL_BUFFER_COUNT
#define TFX_METAL_BUFFER_COUNT 16
#endif
struct tfx_metal_buffer_table {
	RWByteAddressBuffer buffers[TFX_METAL_BUFFER_COUNT];
};
ParameterBlock<tfx_metal_buffer_table> tfx_metal_buffers;

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
