# Shaders

TimelineFX is render agnostic, so you will have to decide how to integrate it into your own renderer. The [implementations](../implementations) folder has complete examples for specific graphics APIs, starting with Direct3D 11.

Every shader a TimelineFX renderer needs is in one Slang file, `timelinefx.slang`: the GPU particle kernel, the billboard and ribbon draws, baked playback and the bounding boxes. It's the exact source the TimelineFX editor draws its preview with, compiled straight from this folder, so it's always up to date with the library. The library itself never compiles or references it, because it never talks to a GPU.

Slang compiles one source for every API, so you compile the file for your target rather than porting it. A contents list at its top names each section with a keyword that appears just once more, at the section itself, so searching for it (in Vim, `*` on the keyword) jumps there.

`archive/` holds the GLSL shaders this file replaced, for reference only.

## Compiling

Compile each entry point on its own with slangc, which comes with the Vulkan SDK and with the Slang releases on GitHub, or through the Slang API.

| Target | Command |
| --- | --- |
| Vulkan | `slangc timelinefx.slang -entry <entry point> -target spirv -profile spirv_1_5 -O3 -o <entry point>.spv` |
| D3D11 | `slangc timelinefx.slang -entry <entry point> -target dxbc -profile <vs_5_0, ps_5_0 or cs_5_0> -D TFX_FIXED_BINDINGS -line-directive-mode none -o <entry point>.cso` |

- **Vulkan needs `shaderDrawParameters`**, a core Vulkan 1.1 feature. Slang reads `SV_VertexID` as `gl_VertexIndex - gl_BaseVertex`, which needs the `DrawParameters` capability.
- **D3D11:** fxc fails with "internal error: flattened side effect" on the `#line` directives Slang otherwise emits, hence `-line-directive-mode none`. Leave out `-stage`: the entry point declares it, and slangc warns if you pass it as well.
- **Matrices are column major** whatever your Slang session defaults to (the Slang API defaults to row major, slangc to column major), because every matrix in the file is marked `column_major`. Upload them column major.

## Entry points

| Entry point | Stage | Group size | What it does | Push constants |
| --- | --- | --- | --- | --- |
| `simulate_particles` | compute | 256 | Simulates `run_on_gpu` emitters and writes their instances | `tfx_gpu_particle_push_t` |
| `tfx_billboard_vertex` | vertex | | The live billboard quads, interpolated from the previous update | `tfx_billboard_push` |
| `tfx_billboard_fragment` | fragment | | Shades live and baked billboards | `tfx_billboard_push` |
| `tfx_ribbon_compute` | compute | 256 | Builds a ribbon bucket's vertex and index buffers, one thread per segment | `tfx_ribbon_bucket_globals_t` |
| `tfx_ribbon_vertex` | vertex | | Transforms the ribbon vertexes and samples the over length shading graphs | `tfx_ribbon_bucket_globals_t` |
| `tfx_ribbon_fragment` | fragment | | Shades ribbons, fading out the ones thinner than a pixel and a half | `tfx_ribbon_bucket_globals_t` |
| `tfx_sprite_playback` | compute | 128 | Interpolates baked sprite data into billboard instances | `tfx_sprite_data_push_t` |
| `tfx_billboard_baked_vertex` | vertex | | The billboard quads for the instances `tfx_sprite_playback` writes | `tfx_baked_billboard_push` |
| `tfx_ribbon_playback` | compute | 128 | Interpolates baked ribbon data into the ribbons `tfx_ribbon_compute` builds | `tfx_ribbon_data_push_t` |
| `tfx_sprite_bounds` | compute | 128 | Optional: the bounding box of a baked frame's sprites | `tfx_sprite_data_push_t` |
| `tfx_ribbon_bounds` | compute | 128 | Optional: the bounding box of a baked frame's ribbons | `tfx_ribbon_bounds_push` |

Dispatch the compute entry points with enough groups to cover their count, `(count + group size - 1) / group size`. A thread past the end returns. The bounds passes are the exception: each thread bounds two items half a group apart, so they take half as many groups.

## Push constants

The library declares `tfx_gpu_particle_push_t`, `tfx_ribbon_bucket_globals_t`, `tfx_sprite_data_push_t` and `tfx_ribbon_data_push_t` in timelinefx.h, and fills the ribbon globals for you per bucket. The other push constants are yours to declare. Their members are `uint`s unless the table says otherwise:

| Struct | Members, in order |
| --- | --- |
| `tfx_billboard_push` | `color_ramp_texture_index`, `image_data_index`, `properties_index`, `sampler_index`, `prev_billboards_index`, `index_offset`, `uniform_index` |
| `tfx_baked_billboard_push` | The same seven, then `ndc_offset_x` and `ndc_offset_y` as floats |
| `tfx_ribbon_bounds_push` | `ribbons_total`, `animation_instances_index`, `ribbon_data_index`, `ribbon_segments_index`, `ribbon_properties_index`, `bounding_boxes_index`, `graphs_index` |

`index_offset` is added to every captured index, for a draw whose previous instances don't start at 0. The ndc offsets move a baked draw in clip space, which the editor uses for the neighbouring tiles of a seamless sprite sheet; leave them at zero otherwise. `tfx_baked_billboard_push` starts with the live push's members so `tfx_billboard_fragment` can shade both.

## Binding models

Every buffer, image and uniform read and write goes through a `tfx_load_*` / `tfx_store_*` accessor, in the file's Buffer_access section for the kernel and its Render_access section for the rest. Each accessor holds every binding model, so nothing needs replacing.

### Bindless: Vulkan and D3D12

The default. Every resource is indexed out of a descriptor array in set 0 by an index in the push constants:

| Binding | Holds |
| --- | --- |
| 0 | Samplers |
| 3 | Images: each shape's `Texture2DArray` of animation frames, and the color ramps |
| 5 | Storage buffers |
| 7 | Uniform buffers: the view, with `view`, `projection` (both 4x4), `screen_size`, `milliseconds`, `timer_lerp` and `update_time` |

The image data from `tfx_GetLibraryGPUShapes` holds each frame's uv rect and `texture_array_index`, which is the shape's bindless image index in the upper 16 bits and the frame's layer in the lower 16. The image a particle samples can change from one particle to the next within a draw, so it's read with `NonUniformResourceIndex`.

The two vertex inputs are vertex buffers:

| Draw | Buffer | Locations and formats |
| --- | --- | --- |
| Billboards | `tfx_instance_t`, instance rate | 0 `position` r32g32b32a32_sfloat, 1 `quaternion` r16g16b16a16_snorm, 2 `size` r16g16_sscaled, 3 `alignment` r8g8b8_snorm, 4 `intensity_gradient_map` r16g16_sscaled, 5 `curved_alpha_life` r8g8b8_unorm, 6 `indexes` r32_uint, 7 `captured_index` r32_uint |
| Ribbons | `tfx_ribbon_vertex_t`, vertex rate | 0 `position` r32g32b32_sfloat, 1 `segment_index` r32_uint, 2 `uv_offset_scale` r32g32_sfloat, 3 `ribbon_index` r32_uint |

D3D12 hasn't been tried yet.

### Fixed bindings: D3D11, OpenGL and WGSL

Compile with `-D TFX_FIXED_BINDINGS` for an API without descriptor arrays. The accessors ignore the push constants' buffer indexes and read one resource per register. The kernel's bases still apply. Every shape comes from one atlas texture array, so `texture_array_index` is only the atlas layer and the uv rects are the shapes' places in it. There are no vertex inputs: the billboard draws read their instance from a buffer by `SV_InstanceID` and unpack it the way the vertex formats above would, and the ribbon draw reads its vertex by `SV_VertexID`.

The kernel uses the one-buffer setup. WGSL always does, define or not, because WebGPU has no descriptor arrays or push constants and allows 8 storage buffers per stage by default.

| Buffer | WGSL / GLSL / SPIR-V (group/set 0) | HLSL |
| --- | --- | --- |
| The pool, read write: the hot section at its start, cold and noise at `cold_base` / `noise_base` | binding 0 | `u0` |
| The instance buffer, read write | binding 1 | `u1` |
| The graph lookups, read only | binding 2 | `t0` |
| The frame's upload from `tfx_CopyGPUParticleUpload`, read only, dispatches at its start | binding 3 | `t1` |
| The push constants | binding 4, a uniform buffer | `b0`, a constant buffer |

The render passes use these registers. D3D11 binds per stage, so a register only has to mean the same thing in the passes that share it. Raw buffers are rows in the C layout; structured buffers have the C struct's size as their stride.

| Register | Holds | Read by |
| --- | --- | --- |
| `b0` | The view, the uniform buffer above | The vertex shaders |
| `b1` | `instance_offset`, `index_offset` (uints), `ndc_offset_x`, `ndc_offset_y` (floats) | The billboard vertex shaders. `SV_InstanceID` leaves out the draw's first instance on D3D11, so pass it in `instance_offset` |
| `b2` | `tfx_ribbon_bucket_globals_t` | The ribbon passes |
| `b3` | `tfx_sprite_data_push_t` | `tfx_sprite_playback` |
| `b4` | `tfx_ribbon_data_push_t` | `tfx_ribbon_playback` |
| `t0` | The instances being drawn, raw | The billboard vertex shaders |
| `t1` | The previous instances, or the baked sprites, raw | `tfx_billboard_vertex`, `tfx_sprite_playback` |
| `t2` | The image data, structured | The vertex shaders |
| `t3` | The emitter properties, structured | The billboard vertex shaders, `tfx_sprite_playback` |
| `t4`, `s0` | The shape atlas and its sampler | The fragment shaders |
| `t5` | The color ramps, a texture array read with an integer load | The fragment shaders |
| `t6` | The ribbon emitters, structured | `tfx_ribbon_compute`, `tfx_ribbon_vertex` |
| `t7` | The graph lookups, raw | `tfx_ribbon_compute`, `tfx_ribbon_vertex` |
| `t8` | The ribbons, structured | `tfx_ribbon_compute`, `tfx_ribbon_vertex` |
| `t9` | The ribbon segments, structured | `tfx_ribbon_compute` |
| `t10` | The ribbon vertexes `tfx_ribbon_compute` wrote, structured | `tfx_ribbon_vertex` |
| `t11`, `t12` | The animation offsets and the animation instances, raw | Both playbacks |
| `t13`, `t14` | The baked ribbons and the baked ribbon properties, structured | `tfx_ribbon_playback` |
| `u0`, `u1` | The ribbon vertexes (structured) and indexes (raw) | `tfx_ribbon_compute` writes them |
| `u2` | The billboard instances, raw | `tfx_sprite_playback` writes them |
| `u3` | The ribbons, structured | `tfx_ribbon_playback` writes them |

The bounds passes have no fixed bindings.

- **D3D11:** [implementations](../implementations) has a working host with the whole build in its CMake. Notes from it:
  - Create raw buffers with `D3D11_RESOURCE_MISC_BUFFER_ALLOW_RAW_VIEWS` and give them raw (`DXGI_FORMAT_R32_TYPELESS`) views. The kernel's upload can be a `DYNAMIC` buffer that `tfx_CopyGPUParticleUpload` fills straight through `Map(WRITE_DISCARD)`.
  - Create the pool zeroed. The library never clears it and it never needs clearing after that.
  - Update the kernel's constant buffer before each dispatch, changing only `dispatch_index`.
  - Unbind what the next pass writes. The instance buffer the kernel writes is the one the vertex shader reads, and a buffer bound for reading and writing at once gets one of its binds cleared, with a debug layer warning.
  - If you keep two instance buffers for interpolation, copy both across when they grow, because the GPU instances in the previous one exist nowhere else.
  - fxc warns that 256 threads times the kernel's temp registers is over its recommended total (X4714), and that `pow` won't work for a negative base (X3571). The base is a unit hash, never negative, and neither shows up at run time.
- **OpenGL:** Slang emits `#version 450`, so OpenGL 4.5. For the kernel, bind the four buffers as SSBOs at bindings 0 to 3 (`glBindBufferBase`) and the push constants as a UBO at binding 4. The render passes haven't been tried on OpenGL.
- **WGSL:** each kernel dispatch needs its own `dispatch_index`, so write one 64 byte block per dispatch at 256 byte steps and bind it with a dynamic offset. The render passes haven't been tried on WebGPU.

### Metal

The kernel has its own path, chosen when Slang targets Metal. Each buffer index in the push constants is a slot in `tfx_metal_buffers`, an argument buffer of `TFX_METAL_BUFFER_COUNT` (default 16, override with `-D`) 64 bit `gpuAddress` values. The accessors read and write every row through it as raw bytes packed exactly like the C structs, because a Metal `float3` is 16 bytes rather than 12. The host:
- fills an `MTLBuffer` with the addresses and points each push constant index at its slot. A section can be its own slot with its address offset and a base of 0, or share one address with the bases.
- binds the push constants and the table at the `[[buffer(n)]]` slots Slang gives them. `simulate_particles` uses 0 and 1. A host entry point with its own push constants gets its own slots, so read them from Slang's reflection.
- calls `useResource` on every buffer the table points to: read for the frame's upload and the graph lookups, read and write for the pool and the instance buffer.

The render passes have no Metal path yet.

## Billboard particles

Draw 6 vertices per instance with no index buffer: `tfx_billboard_vertex` makes the quad's corners from `SV_VertexID`. Each layer's instances from `tfx_GetInstanceBufferByLayer` (or the key ranges from `tfx_GetNextInstanceKeyRange`) are one instanced draw. Besides the instance, the vertex shader reads:
* the previous update's instance buffer, so that it can interpolate the position, size and rotation of each sprite using its `captured_index`. A sprite without the capture flag set (the first frame of a particle, or one that has just looped back to the start of a path) isn't drawn that frame.
* the emitter properties buffer from `tfx_GetParticlePropertiesBuffer` for the image handle, the billboard alignment, the color ramp row and the start frame of the emitter's image.
* the image data buffer for the uvs and the image index.
* the view, with `timer_lerp`, how far the current frame is between the previous update and the current one.

It's an all-in-one shader that handles every billboard alignment mode (camera facing, align to a vector, and camera facing that rolls to a vector) plus the stretch along the direction of travel, so an effect whose emitters mix alignment modes still draws in a single draw call.

`tfx_billboard_fragment` samples the particle image and then uses the red channel (scaled by the sprite's gradient map value) as a lookup into the emitter's row of the color ramp texture. The color ramps are plotted by the library from the emitter's color graphs, with the heat response already baked in, in the format you choose with `tfx_SetColorRampFormat`, and you upload them from `tfx_GetColorRampBitmap`. The ramp color is multiplied by the sprite's intensity, and the curved alpha / sharpness terms do the dissolve fade. It writes premultiplied alpha, so blend with ONE, ONE_MINUS_SRC_ALPHA.

## Ribbons

Ribbons take an extra step compared to billboards. The library writes out ribbon, segment and ribbon emitter data, and `tfx_ribbon_compute` turns that into actual geometry: one thread per segment, sampling the width/scale/clip/morph/noise/lag graphs, building the ribbon frame (facing the camera, at a fixed angle, or following the curve of the path) and writing out the tessellated triangle strip. The segment index is packed into each vertex so the later stages can recover which ribbon, which segment and which side of the strip a vertex came from.

Ribbons are grouped into buckets, one for each segment count and effect sort key, and each bucket is one compute dispatch and one indexed draw. After `tfx_UpdateStage`, copy the ribbon data into your staging buffers with `tfx_CopyRibbonDataToStagingBuffers` and loop over the buckets with `tfx_CreateRibbonDispatch` and `tfx_NextRibbonDispatch`, which give you the offsets, counts and push constants for each dispatch and draw. To draw the ribbons of several stages with one upload, one compute pass and one draw pass, use a ribbon batch instead: `tfx_BeginRibbonBatch`, `tfx_AddStageToRibbonBatch` for each stage, `tfx_GetRibbonBatchSizes` and `tfx_CreateRibbonBatchDispatch`. Both flows are described in full in timelinefx.h.

`tfx_ribbon_vertex` then transforms those vertices, unpacks that segment index, samples the over-length shading graphs (intensity, gradient map, curved alpha, sharpness) and ping-pongs the U coordinate so a texture can wrap repeatedly along the ribbon without a flipped seam at each wrap point.

`tfx_ribbon_fragment` shades a ribbon the same way as a billboard, and fades it out as it gets thinner than about 1.5 pixels. Below that the rasteriser starts dropping fragments and a thin ribbon breaks up into dashes, so fading it gives a tail that thins out smoothly instead.

## GPU simulated particles

`simulate_particles` simulates `run_on_gpu` emitters and writes their instances, with its own graph sampling, noise fields, emitter forces and spawn generation.

### Dispatching

On a stage with a `gpu_particle_capacity`, `tfx_GetGPUParticleFrame` hands you dispatch records, tick records, the rows spawned since the last fetch and a table of spawn locations that relative particles follow. Dispatch the kernel once per record with `tfx_gpu_particle_push_t` as push constants and `(slot_count + 255) / 256` groups. Each thread owns one ring slot: it takes in any spawn that lands in its slot, runs every tick of the fetch, keeps the particle state in a pool you size once from `tfx_GetGPUParticleBufferSizes`, and writes its instance into the same instance buffer as the CPU particles, in the range after them. The billboard shaders draw them unchanged.

Point, User, Area, Ellipse, Disc, Cylinder and Line emitters have their spawns generated in the kernel from a few rows per tick and location, so a burst of any size costs the CPU next to nothing. Icosphere, Path, Other Emitter and Spawn on Ribbon emitters spawn on the CPU and upload a row per particle. Every emitter force runs in the kernel from force rows the library resolves each tick, so they need nothing more from you.

The pool has hot, cold and noise sections, and each section starts on a 256 byte boundary. They can share one buffer or have one each, and the same goes for the frame's arrays: `upload_bytes` and the frame's upload offsets lay them all out in one buffer, which `tfx_CopyGPUParticleUpload` fills. Every buffer in the push constants has its own index plus, for those that can share, a `*_base` row the section starts at. On Vulkan and Metal, where storage buffers are untyped, one buffer and one index with the bases is enough. On D3D12 a structured view carries its row size, so give each section a view of its own, either over the section with a base of 0, or over the whole buffer with the base. The fixed binding targets always use one pool buffer and one upload buffer with the bases.

### Adding your own behaviour

You can add custom forces, collision, spawn placement or changes to how particles are drawn without changing `timelinefx.slang`. Import it into your own Slang module and implement any of the three hook interfaces:

| Interface | What it lets you do |
| --- | --- |
| `ITfxParticleUpdate` | `acceleration` adds to the particle's acceleration every tick like the built-in forces, and `after_move` runs once the particle has moved, where you can also kill it |
| `ITfxParticleSpawn` | `after_spawn` places the particles of User emitters that the kernel generates, the GPU side of `tfx_SetStageSpawnCallback` |
| `ITfxParticleWrite` | `before_write` changes a live particle's position, rotation, size, color ramp position, alpha or image frame just before its instance is written |

Then write your own `[numthreads(tfx_group_size, 1, 1)]` entry point with `tfx_gpu_particle_push` as the first member of your push constants and call `tfx_simulate` with your hooks, using `tfx_no_particle_update`, `tfx_no_particle_spawn` and `tfx_no_particle_write` for the ones you don't need. The GPU particle simulation section of timelinefx.h walks through it step by step. Your hooks only run on the GPU, so emitters that fall back to the CPU (during warmup, for example) simulate without them.

A host entry point with its own push constants does what `simulate_particles` does: declare them twice, `[[vk::push_constant]]` and `[[vk::binding(4, 0)]] ... : register(b15)` (the register keeps it out of Metal's buffer slots), and read them through a function that picks with `__target_switch`: `case wgsl:` (plus `case glsl:` under `TFX_FIXED_BINDINGS`) returns the uniform one, `default:` the push constant one.

## Pre-baked effects

For maximum speed you can pre-bake an effect to sprite data and update the instance buffer on the GPU instead of simulating on the CPU.

Every frame, after `tfx_UpdateAnimationManager`, dispatch `tfx_sprite_playback` with `tfx_sprite_data_push_t`. It interpolates each baked sprite between its frame and the one before, resolves the image frame (the emitter's start frame included) and bakes the stretch into `position.w`, then writes a billboard instance. Draw those with `tfx_billboard_baked_vertex` and `tfx_billboard_fragment`. The baked vertex shader reads each instance as it is, with no previous buffer and no timer lerp: `tfx_billboard_vertex` would add the start frame a second time and see no motion to stretch along.

Baked ribbons go through the same `tfx_ribbon_compute` and ribbon draw as live ones. Dispatch `tfx_ribbon_playback` with `tfx_ribbon_data_push_t`, then walk the animation manager's buckets with `tfx_CreateAnimationRibbonDispatch` and `tfx_NextAnimationRibbonDispatch` exactly as you walk a stage's: one `tfx_ribbon_compute` dispatch and one indexed draw per bucket. `tfx_CopyAnimationRibbonEmitters` gives you the emitter buffer `tfx_ribbon_compute` reads for baked ribbons, and `tfx_GetAnimationRibbonBufferSizes` the buffer sizes for the frame. Baked ribbons carry their own graph rows, so bind `tfx_GetAnimationRibbonGraphLookupsBufferPointer` in place of `tfx_GetGPUGraphLookupsBuffer` for their compute and draw. That way sprite data plays back without a library loaded.

`tfx_sprite_bounds` and `tfx_ribbon_bounds` are optional extras rather than part of the render path. Each group reduces its sprites, or its ribbons' segments, to one box, and you reduce the groups' boxes on the CPU. The editor uses them to auto-fit the camera and to size sprite sheet exports; in a game they're useful for culling a baked effect.

[See here for an example of how to use the library to draw pre-baked effect data](https://github.com/peterigz/zest/tree/main/examples/SDL2/zest-timelinefx-prerecorded-effects)
