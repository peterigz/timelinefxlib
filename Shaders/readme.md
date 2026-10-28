# Shader Examples

TimelineFX is render agnostic, so you will have to decide how to integrate it into your own renderer. The [implementations](../implementations) folder has complete examples for specific graphics APIs, starting with Direct3D 11.

The shaders in this folder are the exact shaders that the TimelineFX editor uses to draw its preview. The editor compiles them straight from this folder, so they're always up to date with the library, and they're a working example of what you need in order to render TimelineFX particles and ribbons yourself. The library itself never compiles or references them - it never talks to a GPU - so they're a reference implementation you can port to your own renderer and pipeline setup.

The render and baked playback shaders are GLSL (Vulkan dialect, `#version 450`) and they use bindless resources (`GL_EXT_nonuniform_qualifier`), so every texture and buffer is indexed out of a global array by an index passed through the push constants. The GPU simulation kernel is Slang so that it can compile for every graphics API from one file.

Every shape in a library has its own image (a `texture2DArray` holding its animation frames), so the image a particle samples can change from one particle to the next within a draw. The image data buffer you upload from `tfx_GetLibraryGPUShapes` holds each frame's uv rect and the bindless index of its image, and the shaders read the image with `nonuniformEXT`.

## Billboard particles

| Shader | What it does |
| --- | --- |
| `timelinefx3d.vert` | Builds the billboard quads for particle sprites |
| `timelinefx.frag` | Shades them |

The vertex shader is instanced. There's no vertex or index buffer for the quad, the corners are generated in the shader from `gl_VertexIndex`, so you draw 6 vertices per instance with the instance buffer of `tfx_instance_t` sprites from the library bound as an instance rate vertex buffer. It also reads:
* the previous frame's instance buffer, as a storage buffer, so that it can interpolate the position, size and rotation of each sprite using its `captured_index`. A sprite without the capture flag set (the first frame of a particle, or one that has just looped back to the start of a path) isn't drawn that frame.
* the emitter properties buffer from `tfx_GetParticlePropertiesBuffer` for the image handle, the billboard alignment, the color ramp row and the start frame of the emitter's image.
* the image data buffer for the uvs and the image index.
* a uniform buffer with the view and projection matrices and `timer_lerp`, how far the current frame is between the previous update and the current one.

It's an all-in-one shader that handles every billboard alignment mode (camera facing, align to a vector, and camera facing that rolls to a vector) plus the stretch along the direction of travel, so an effect whose emitters mix alignment modes still draws in a single draw call. If your effects only ever use one alignment mode you can strip the branches out for a little more speed. Define `NO_LERP` if you update once per frame and don't need the interpolation.

The fragment shader samples the particle image and then uses the red channel (scaled by the sprite's gradient map value) as a lookup into the emitter's row of the color ramp texture. The color ramps are plotted by the library from the emitter's color graphs, with the heat response already baked in, in the format you choose with `tfx_SetColorRampFormat`, and you upload them from `tfx_GetColorRampBitmap`. The ramp color is multiplied by the sprite's intensity, and the curved alpha / sharpness terms do the dissolve fade.

## Ribbons

| Shader | What it does |
| --- | --- |
| `ribbons.comp` | Builds the ribbon vertex and index buffers on the GPU |
| `ribbon_3d.vert` | Transforms those vertices and does the per-vertex shading lookups |
| `ribbon.frag` | Shades it, same color ramp scheme as the particle fragment shader |

Ribbons take an extra step compared to billboards. The library writes out ribbon, segment and ribbon emitter data, and `ribbons.comp` turns that into actual geometry: one invocation per segment, sampling the width/scale/clip/morph/noise/lag graphs, building the ribbon frame (facing the camera, at a fixed angle, or following the curve of the path) and writing out the tessellated triangle strip. The segment index is packed into each vertex so the later stages can recover which ribbon, which segment and which side of the strip a vertex came from.

Ribbons are grouped into buckets, one for each segment count and effect sort key, and each bucket is one compute dispatch and one draw. After `tfx_UpdateStage`, copy the ribbon data into your staging buffers with `tfx_CopyRibbonDataToStagingBuffers` and loop over the buckets with `tfx_CreateRibbonDispatch` and `tfx_NextRibbonDispatch`, which give you the offsets, counts and push constants for each dispatch and draw. To draw the ribbons of several stages with one upload, one compute pass and one draw pass, use a ribbon batch instead: `tfx_BeginRibbonBatch`, `tfx_AddStageToRibbonBatch` for each stage, `tfx_GetRibbonBatchSizes` and `tfx_CreateRibbonBatchDispatch`. Both flows are described in full in timelinefx.h.

The vertex shader then transforms those vertices, unpacks that segment index, samples the over-length shading graphs (intensity, gradient map, curved alpha, sharpness) and ping-pongs the U coordinate so a texture can wrap repeatedly along the ribbon without a flipped seam at each wrap point.

The fragment shader also fades a ribbon out as it gets thinner than about 1.5 pixels. Below that the rasteriser starts dropping fragments and a thin ribbon breaks up into dashes, so fading it gives a tail that thins out smoothly instead.

## GPU simulated particles

| Shader | What it does |
| --- | --- |
| `tfx_gpu_particles.slang` | Simulates `run_on_gpu` emitters and writes their instances: the whole kernel in one file, with its own graph sampling, noise fields, emitter forces and spawn generation |

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
- **D3D11:** compile with `-target dxbc -profile cs_5_0 -line-directive-mode none` and no `-stage` (the entry point declares it, and slangc warns if you pass it as well). fxc fails with "internal error: flattened side effect" on the `#line` directives Slang otherwise emits. [implementations](../implementations) has a working D3D11 host with this build in its CMake. Notes from it:
  - Create the pool, instance, graph and upload buffers with `D3D11_RESOURCE_MISC_BUFFER_ALLOW_RAW_VIEWS` and give them raw (`DXGI_FORMAT_R32_TYPELESS`) views. The upload can be a `DYNAMIC` buffer that `tfx_CopyGPUParticleUpload` fills straight through `Map(WRITE_DISCARD)`.
  - Create the pool zeroed. The library never clears it and it never needs clearing after that.
  - Update the constant buffer before each dispatch, changing only `dispatch_index`.
  - The instance buffer the kernel writes is the one your vertex shader reads, so unbind it from the vertex shader after drawing. A buffer bound for reading and writing at once gets one of its binds cleared, with a debug layer warning.
  - If you keep two instance buffers for interpolation, copy both across when they grow, because the GPU instances in the previous one exist nowhere else.
  - fxc warns that 256 threads times the kernel's temp registers is over its recommended total (X4714), and that a few functions may use an uninitialized variable (X4000). Neither shows up at run time: about 106,000 particles take about 110 µs per dispatching fetch, on par with Vulkan.
- **OpenGL:** Slang emits `#version 450`, so OpenGL 4.5. Bind the four buffers as SSBOs at bindings 0 to 3 (`glBindBufferBase`) and the push constants as a UBO at binding 4.
- **A host entry point with its own push constants** does what `simulate_particles` does: declare them twice, `[[vk::push_constant]]` and `[[vk::binding(4, 0)]] ... : register(b15)` (the register keeps it out of Metal's buffer slots), and read them through a function that picks with `__target_switch`: `case wgsl:` (plus `case glsl:` under `TFX_FIXED_BINDINGS`) returns the uniform one, `default:` the push constant one.

### Metal

Metal has its own path, chosen when Slang targets Metal. Each buffer index in the push constants is a slot in `tfx_metal_buffers`, an argument buffer of `TFX_METAL_BUFFER_COUNT` (default 16, override with `-D`) 64 bit `gpuAddress` values. The accessors read and write every row through it as raw bytes packed exactly like the C structs, because a Metal `float3` is 16 bytes rather than 12. The host:
- fills an `MTLBuffer` with the addresses and points each push constant index at its slot. A section can be its own slot with its address offset and a base of 0, or share one address with the bases.
- binds the push constants and the table at the `[[buffer(n)]]` slots Slang gives them. `simulate_particles` uses 0 and 1. A host entry point with its own push constants gets its own slots, so read them from Slang's reflection.
- calls `useResource` on every buffer the table points to: read for the frame's upload and the graph lookups, read and write for the pool and the instance buffer.

### Dispatching

On a stage with a `gpu_particle_capacity`, `tfx_GetGPUParticleFrame` hands you dispatch records, tick records, the rows spawned since the last fetch and a table of spawn locations that relative particles follow. Dispatch the kernel once per record with `tfx_gpu_particle_push_t` as push constants and `(slot_count + 255) / 256` groups. Each thread owns one ring slot: it takes in any spawn that lands in its slot, runs every tick of the fetch, keeps the particle state in a pool you size once from `tfx_GetGPUParticleBufferSizes`, and writes its instance into the same instance buffer as the CPU particles, in the range after them. The vertex and fragment shaders above draw them unchanged.

Point, User, Area, Ellipse, Disc, Cylinder and Line emitters have their spawns generated in the kernel from a few rows per tick and location, so a burst of any size costs the CPU next to nothing. Icosphere, Path, Other Emitter and Spawn on Ribbon emitters spawn on the CPU and upload a row per particle. Every emitter force runs in the kernel from force rows the library resolves each tick, so they need nothing more from you.

The pool has hot, cold and noise sections, and each section starts on a 256 byte boundary. They can share one buffer or have one each, and the same goes for the frame's arrays: `upload_bytes` and the frame's upload offsets lay them all out in one buffer, which `tfx_CopyGPUParticleUpload` fills. Every buffer in the push constants has its own index plus, for those that can share, a `*_base` row the section starts at. On Vulkan and Metal, where storage buffers are untyped, one buffer and one index with the bases is enough. On D3D12 a structured view carries its row size, so give each section a view of its own, either over the section with a base of 0, or over the whole buffer with the base. The fixed binding targets above always use one pool buffer and one upload buffer with the bases.

### Adding your own behaviour

You can add custom forces, collision, spawn placement or changes to how particles are drawn without changing `tfx_gpu_particles.slang`. Import it into your own Slang module and implement any of the three hook interfaces:

| Interface | What it lets you do |
| --- | --- |
| `ITfxParticleUpdate` | `acceleration` adds to the particle's acceleration every tick like the built-in forces, and `after_move` runs once the particle has moved, where you can also kill it |
| `ITfxParticleSpawn` | `after_spawn` places the particles of User emitters that the kernel generates, the GPU side of `tfx_SetStageSpawnCallback` |
| `ITfxParticleWrite` | `before_write` changes a live particle's position, rotation, size, color ramp position, alpha or image frame just before its instance is written |

Then write your own `[numthreads(tfx_group_size, 1, 1)]` entry point with `tfx_gpu_particle_push` as the first member of your push constants and call `tfx_simulate` with your hooks, using `tfx_no_particle_update`, `tfx_no_particle_spawn` and `tfx_no_particle_write` for the ones you don't need. The GPU particle simulation section of timelinefx.h walks through it step by step. Your hooks only run on the GPU, so emitters that fall back to the CPU (during warmup, for example) simulate without them.

## Compute shaders for pre-baked effects

For maximum speed you can pre-bake an effect to sprite data and update the instance buffer on the GPU instead of simulating on the CPU:

| Shader | What it does |
| --- | --- |
| `sprite_data_playback.comp` | Interpolates baked billboard sprite data between the previous and current frame and writes the instance buffer |
| `static_billboard.vert` | Builds the billboard quads for the instances `sprite_data_playback.comp` writes |
| `static_sprite.frag` | Shades them, the same as `timelinefx.frag` |
| `ribbon_data_playback.comp` | The same thing for baked ribbon data, writing the ribbons that `ribbons.comp` builds into geometry |

Baked billboards need their own vertex shader. `sprite_data_playback.comp` has already done the interpolation between baked frames, resolved the image frame (the emitter's start frame included) and baked the stretch into `position.w`, so `static_billboard.vert` reads each instance as it is, with no previous buffer and no timer lerp. `timelinefx3d.vert` would add the start frame a second time and see no motion to stretch along. It also takes a `vec4 offset` and `ndc_offset_x/y` in its push constants, which the editor uses to draw seamless sprite sheets; leave them at zero otherwise.

Baked ribbons go through the same `ribbons.comp` and ribbon shaders as live ones. Every frame, after `tfx_UpdateAnimationManager`, dispatch `ribbon_data_playback.comp` with `tfx_ribbon_data_push_t`, then walk the animation manager's buckets with `tfx_CreateAnimationRibbonDispatch` and `tfx_NextAnimationRibbonDispatch` exactly as you walk a stage's: one `ribbons.comp` dispatch and one indexed draw per bucket. `tfx_CopyAnimationRibbonEmitters` gives you the emitter buffer `ribbons.comp` reads for baked ribbons, and `tfx_GetAnimationRibbonBufferSizes` the buffer sizes for the frame. Baked ribbons carry their own graph rows, so bind `tfx_GetAnimationRibbonGraphLookupsBufferPointer` in place of `tfx_GetGPUGraphLookupsBuffer` for their compute and draw. That way sprite data plays back without a library loaded.

The bounding box shaders are optional extras rather than part of the render path:

| Shader | What it does |
| --- | --- |
| `sprite_data_bounding_box.comp` | Parallel reduction over baked sprites to find the min/max bounds of an effect |
| `ribbon_data_bounding_box.comp` | The same for baked ribbon segments |

The editor uses these to auto-fit the camera and to size sprite sheet exports; they're useful in a game for culling a baked effect.

[See here for an example of how to use the library to draw pre-baked effect data](https://github.com/peterigz/zest/tree/main/examples/SDL2/zest-timelinefx-prerecorded-effects)
