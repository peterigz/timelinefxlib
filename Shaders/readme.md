# Shader Examples

TimelineFX is render agnostic, so you will have to decide how to integrate it into your own renderer. (I'll be looking to add integrations into popular renderers over time).

The shaders in this folder are the exact shaders that the TimelineFX editor uses to draw its preview, copied here as a working example of what you need in order to render TimelineFX particles and ribbons yourself. They're not compiled or referenced by the library itself - it never talks to a GPU - they're a reference implementation you can port to your own renderer and pipeline setup. They're all GLSL (Vulkan dialect, `#version 450`), and they use bindless resources (`GL_EXT_nonuniform_qualifier`), so every texture and buffer is indexed out of a global array by an index passed through the push constants.

## Billboard particles

| Shader | What it does |
| --- | --- |
| `timelinefx3d.vert` | Builds the billboard quads for particle sprites |
| `timelinefx.frag` | Shades them |

The vertex shader is instanced: there are no vertex or index buffers, the quad corners are generated in the shader from `gl_VertexIndex` and the only thing you upload each frame is the instance buffer of `tfx_instance_t` sprites that the library produced. It's an all-in-one shader that handles every billboard alignment mode (camera facing, free align, vector align) plus the interpolation between the previous and current frame, so an effect whose emitters mix alignment modes still draws in a single draw call. If your effects only ever use one alignment mode you can strip the branches out for a little more speed.

The fragment shader samples the particle image and then uses the red channel (scaled by the sprite's gradient map value) as a lookup into the effect's color ramp texture, which is where the color over lifetime, the intensity and the heat response come from. The curved alpha / sharpness terms do the dissolve fade.

## Ribbons

| Shader | What it does |
| --- | --- |
| `ribbons.comp` | Builds the ribbon vertex and index buffers on the GPU |
| `ribbon_3d.vert` | Transforms those vertices and does the per-vertex shading lookups |
| `ribbon.frag` | Shades it, same color ramp scheme as the particle fragment shader |

Ribbons take an extra step compared to billboards. The library writes out ribbon and segment data, and `ribbons.comp` turns that into actual geometry: one invocation per segment, sampling the width/scale/clip/morph/noise/lag graphs, building the ribbon frame against the camera position and writing out the tessellated triangle strip. The segment index is packed into each vertex so the later stages can recover which ribbon, which segment and which side of the strip a vertex came from.

The vertex shader then transforms those vertices, unpacks that segment index, samples the over-length shading graphs (intensity, gradient map, curved alpha, sharpness) and ping-pongs the U coordinate so a texture can wrap repeatedly along the ribbon without a flipped seam at each wrap point.

## GPU simulated particles

| Shader | What it does |
| --- | --- |
| `tfx_gpu_particles.slang` | Simulates `run_on_gpu` emitters and writes their instances: the whole kernel in one file, with its own graph sampling, noise fields and spawn generation |

It's Slang rather than GLSL, so it can carry user hooks and compile for every API from one source. Compile `simulate_particles` from `tfx_gpu_particles.slang`, a single module like the library's single header. A contents list at its top names each section with a keyword that appears just once more, at the section itself, so searching for it (in Vim, `*` on the keyword) jumps there.

The kernel reaches its buffers only through the `tfx_load_*` / `tfx_store_*` functions in its buffer access section, which take rows relative to their own array and add the push constants' bases. By default they index bindless arrays, which suits Vulkan and D3D12. Every other binding model is handled in the same file, so nothing needs replacing.

### Fixed bindings: WGSL, D3D11 and OpenGL

WGSL always uses fixed bindings, because WebGPU has no descriptor arrays or push constants and allows 8 storage buffers per stage by default. Compile with `-D TFX_FIXED_BINDINGS` to use the same layout on every other target, which is what D3D11 and OpenGL need. The kernel binds the one-buffer setup and ignores the push constants' buffer indexes (the bases still apply). Rows are read as raw bytes in the C layout.

| Buffer | WGSL / GLSL / SPIR-V (group/set 0) | HLSL |
| --- | --- | --- |
| The pool, read write: the hot section at its start, cold and noise at `cold_base` / `noise_base` | binding 0 | `u0` |
| The instance buffer, read write | binding 1 | `u1` |
| The graph lookups, read only | binding 2 | `t0` |
| The frame's upload from `tfx_CopyGPUParticleUpload`, read only, dispatches at its start | binding 3 | `t1` |
| The push constants | binding 4, a uniform buffer | `b0`, a constant buffer |

- **WGSL:** each dispatch needs its own `dispatch_index`, so write one 64 byte block per dispatch at 256 byte steps and bind it with a dynamic offset.
- **D3D11:** compile with `-target dxbc -profile cs_5_0 -line-directive-mode none`. fxc fails with "internal error: flattened side effect" on the `#line` directives Slang otherwise emits. Create the pool, instance, graph and upload buffers with `D3D11_RESOURCE_MISC_BUFFER_ALLOW_RAW_VIEWS` and give them raw (`DXGI_FORMAT_R32_TYPELESS`) views. Update the constant buffer before each dispatch. fxc warns that 256 threads times the kernel's temp registers is over its recommended total, so measure there.
- **OpenGL:** Slang emits `#version 450`, so OpenGL 4.5. Bind the four buffers as SSBOs at bindings 0 to 3 (`glBindBufferBase`) and the push constants as a UBO at binding 4.
- **A host entry point with its own push constants** does what `simulate_particles` does: declare them twice, `[[vk::push_constant]]` and `[[vk::binding(4, 0)]] ... : register(b15)` (the register keeps it out of Metal's buffer slots), and read them through a function that picks with `__target_switch`: `case wgsl:` (plus `case glsl:` under `TFX_FIXED_BINDINGS`) returns the uniform one, `default:` the push constant one.

### Metal

Metal has its own path, chosen when Slang targets Metal. Each buffer index in the push constants is a slot in `tfx_metal_buffers`, an argument buffer of `TFX_METAL_BUFFER_COUNT` (default 16, override with `-D`) 64 bit `gpuAddress` values. The accessors read and write every row through it as raw bytes packed exactly like the C structs, because a Metal `float3` is 16 bytes rather than 12. The host:
- fills an `MTLBuffer` with the addresses and points each push constant index at its slot. A section can be its own slot with its address offset and a base of 0, or share one address with the bases.
- binds the push constants and the table at the `[[buffer(n)]]` slots Slang gives them. `simulate_particles` uses 0 and 1. A host entry point with its own push constants gets its own slots, so read them from Slang's reflection.
- calls `useResource` on every buffer the table points to: read for the frame's upload and the graph lookups, read and write for the pool and the instance buffer.

On a stage with a `gpu_particle_capacity`, `tfx_GetGPUParticleFrame` hands you dispatch records, tick records, the rows spawned since the last fetch and a table of spawn locations that relative particles follow. Dispatch the kernel once per record with `tfx_gpu_particle_push_t` as push constants and `(slot_count + 255) / 256` groups. Each thread owns one ring slot: it takes in any spawn that lands in its slot, runs every tick of the fetch, keeps the particle state in a pool you size once from `tfx_GetGPUParticleBufferSizes`, and writes its instance into the same instance buffer as the CPU particles, in the range after them. The vertex and fragment shaders above draw them unchanged.

The pool has hot, cold and noise sections, and each section starts on a 256 byte boundary. They can share one buffer or have one each, and the same goes for the frame's arrays: `upload_bytes` and the frame's upload offsets lay them all out in one buffer, which `tfx_CopyGPUParticleUpload` fills. Every buffer in the push constants has its own index plus, for those that can share, a `*_base` row the section starts at. On Vulkan and Metal, where storage buffers are untyped, one buffer and one index with the bases is enough. On D3D12 a structured view carries its row size, so give each section a view of its own, either over the section with a base of 0, or over the whole buffer with the base. The fixed binding targets above always use one pool buffer and one upload buffer with the bases.

## Compute shaders for pre-baked effects

For maximum speed you can pre-bake an effect to sprite data and update the instance buffer on the GPU instead of simulating on the CPU:

| Shader | What it does |
| --- | --- |
| `sprite_data_playback.comp` | Interpolates baked billboard sprite data between the previous and current frame and writes the instance buffer |
| `ribbon_data_playback.comp` | The same thing for baked ribbon data |

The instance data these write is consumed by the exact same vertex shaders above, so the draw side of your renderer doesn't need to know whether an effect is being simulated live or played back from baked data.

The bounding box shaders are optional extras rather than part of the render path:

| Shader | What it does |
| --- | --- |
| `sprite_data_bounding_box.comp` | Parallel reduction over baked sprites to find the min/max bounds of an effect |
| `ribbon_data_bounding_box.comp` | The same for baked ribbon segments |

The editor uses these to auto-fit the camera and to size sprite sheet exports; they're useful in a game for culling a baked effect.

[See here for an example of how to use the library to draw pre-baked effect data](https://github.com/peterigz/zest/tree/main/examples/SDL2/zest-timelinefx-prerecorded-effects)
