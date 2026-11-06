# Shaders

TimelineFX is render agnostic, so you will have to decide how to integrate it into your own renderer. The [implementations](../implementations) folder has complete examples for specific graphics APIs, starting with Direct3D 11.

Every shader a TimelineFX renderer needs is in one Slang file, `timelinefx.slang`: the GPU particle kernel, the billboard and ribbon draws, baked playback and the bounding boxes. It's the exact source the TimelineFX editor draws its preview with, compiled straight from this folder, so it's always up to date with the library. The library itself never compiles or references it, because it never talks to a GPU.

Slang compiles one source for every API, so you compile the file for your target rather than porting it. A contents list at its top names each section with a keyword that appears just once more, at the section itself, so searching for it (in Vim, `*` on the keyword) jumps there.

`archive/` holds the GLSL shaders this file replaced, for reference only.

## Compiling

Compile each entry point on its own with slangc, which comes with the Vulkan SDK and with the Slang releases on GitHub, or through the Slang API.

| Target | Command |
| --- | --- |
| Vulkan, bindless | `slangc timelinefx.slang -entry <entry point> -target spirv -profile spirv_1_5 -O3 -fvk-use-entrypoint-name -o <entry point>.spv` |
| D3D12, bindless | `slangc timelinefx.slang -entry <entry point> -target dxil -profile <vs_6_0, ps_6_0 or cs_6_0> -o <entry point>.dxil` |
| Fixed bindings, SPIR-V | `slangc timelinefx.slang -entry <entry point> -target spirv -profile spirv_1_5 -O3 -fvk-use-entrypoint-name -D TFX_FIXED_BINDINGS -D TFX_ENTRY_<entry point> -o <entry point>.spv` |
| Fixed bindings, DXIL | `slangc timelinefx.slang -entry <entry point> -target dxil -profile <vs_6_0, ps_6_0 or cs_6_0> -D TFX_FIXED_BINDINGS -D TFX_ENTRY_<entry point> -o <entry point>.dxil` |
| D3D11 | `slangc timelinefx.slang -entry <entry point> -target dxbc -profile <vs_5_0, ps_5_0 or cs_5_0> -D TFX_FIXED_BINDINGS -D TFX_ENTRY_<entry point> -line-directive-mode none -o <entry point>.cso` |
| WebGPU | `slangc timelinefx.slang -entry <entry point> -target wgsl -D TFX_FIXED_BINDINGS -D TFX_ENTRY_<entry point> -o <entry point>.wgsl` |
| OpenGL | `slangc timelinefx.slang -entry <entry point> -target spirv -profile spirv_1_5 -O3 -D TFX_FIXED_BINDINGS -D TFX_COMBINED_SAMPLERS -D TFX_ENTRY_<entry point> -o <entry point>.spv`, then `spirv-cross <entry point>.spv --version 430 --no-es --no-support-nonzero-baseinstance --output <entry point>.glsl` |
| Metal | `slangc timelinefx.slang -entry <entry point> -target metal -D TFX_FIXED_BINDINGS -D TFX_METAL_SLOTS -D TFX_ENTRY_<entry point> -o <entry point>.metal`, and `simulate_particles` with no defines |

- **`TFX_ENTRY_<entry point>`** is the entry point's name as it is, for example `-D TFX_ENTRY_tfx_billboard_vertex`. Under the fixed bindings it picks that entry point's resources and slots, so the compile has to have it. The kernel's slots don't depend on it.
- **`-fvk-use-entrypoint-name`** keeps the entry point's own name in the SPIR-V. Without it slangc names every one `main`, which is the name to give your API instead. DXIL, DXBC and WGSL always keep the name.
- **`spirv_1_5` needs a Vulkan 1.2 device.** The bindless model needs 1.2 anyway, for descriptor indexing. Under the fixed bindings `-profile spirv_1_3` also validates for Vulkan 1.1 (checked with spirv-val, not yet on a device). SDL_gpu creates a Vulkan 1.0 instance unless you pass `SDL_GPUVulkanOptions` with `vulkan_api_version` through `SDL_PROP_GPU_DEVICE_CREATE_VULKAN_OPTIONS_POINTER`.
- **Vulkan needs `shaderDrawParameters`**, a core Vulkan 1.1 feature. Slang reads `SV_VertexID` as `gl_VertexIndex - gl_BaseVertex`, which needs the `DrawParameters` capability. SDL_gpu only turns it on when you ask for it in `feature_list`, with Vulkan 1.2 for `spirv_1_5`:

  ```c
  VkPhysicalDeviceVulkan11Features vulkan_11_features = { VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_1_FEATURES };
  vulkan_11_features.shaderDrawParameters = VK_TRUE;
  SDL_GPUVulkanOptions vulkan_options = { 0 };
  vulkan_options.vulkan_api_version = VK_MAKE_API_VERSION(0, 1, 2, 0);
  vulkan_options.feature_list = &vulkan_11_features;	//Only read at 1.2 and up
  SDL_SetPointerProperty(properties, SDL_PROP_GPU_DEVICE_CREATE_VULKAN_OPTIONS_POINTER, &vulkan_options);
  ```

  Both have to stay alive until `SDL_CreateGPUDeviceWithProperties` returns. The Vulkan SDK's `Include` folder has its own, older SDL3 headers, so adding it to your include path for `vulkan_core.h` can hide `SDL_GPUVulkanOptions`. Put SDL's include folder first, or take the Vulkan headers from somewhere else. Without these options the shaders still draw on some drivers, so only the validation layer tells you.
- **OpenGL goes through SPIR-V and SPIRV-Cross**, as sokol-shdc does. Slang's own `-target glsl` is Vulkan GLSL: the kernel's happens to compile for OpenGL, but the vertex shaders use `gl_VertexIndex` and `gl_InstanceIndex`, which OpenGL doesn't have. spirv-cross ships with the Vulkan SDK. See the OpenGL notes under the fixed bindings.
- **Metal takes the text** `-target metal` writes, compiled at run time with `newLibraryWithSource`, or ahead of time with `xcrun -sdk macosx metal`, which needs the Metal toolchain that recent Xcode versions download separately (`xcodebuild -downloadComponent MetalToolchain`). Slang keeps its `#line` directives, so Metal's errors point at lines of `timelinefx.slang`. See the Metal section for the defines.
- **D3D11:** fxc fails with "internal error: flattened side effect" on the `#line` directives Slang otherwise emits, hence `-line-directive-mode none`. Leave out `-stage`: the entry point declares it, and slangc warns if you pass it as well. fxc has no register spaces, so leave every `*_SPACE` macro at 0.
- **Matrices are column major** whatever your Slang session defaults to (the Slang API defaults to row major, slangc to column major), because every matrix in the file is marked `column_major`. Upload them column major.

## Entry points

| Entry point | Stage | Group size | What it does | Push constants |
| --- | --- | --- | --- | --- |
| `simulate_particles` | compute | 256 | Simulates `run_on_gpu` emitters and writes their instances | `tfx_gpu_particle_push_t` |
| `tfx_billboard_vertex` | vertex | | The live billboard quads, interpolated from the previous update | `tfx_billboard_push` |
| `tfx_billboard_fragment` | fragment | | Shades live and baked billboards | `tfx_billboard_push` |
| `tfx_billboard_cutout_fragment` | fragment | | Shades live and baked billboards on a layer that writes depth, cut at the emitter's alpha threshold | `tfx_billboard_push` |
| `tfx_billboard_lit_vertex` | vertex | | `tfx_billboard_vertex` for a lit layer, adding the quad's world axes and the emitter's lighting | `tfx_billboard_push` |
| `tfx_billboard_lit_fragment` | fragment | | Lights live and baked billboards on a lit layer | `tfx_billboard_push` |
| `tfx_billboard_lit_cutout_fragment` | fragment | | Lights live and baked billboards on a lit layer that writes depth, cut at the alpha threshold | `tfx_billboard_push` |
| `tfx_ribbon_compute` | compute | 256 | Builds a ribbon bucket's vertex and index buffers, one thread per segment | `tfx_ribbon_bucket_globals_t` |
| `tfx_ribbon_vertex` | vertex | | Transforms the ribbon vertexes and samples the over length shading graphs | `tfx_ribbon_bucket_globals_t` |
| `tfx_ribbon_fragment` | fragment | | Shades ribbons, fading out the ones thinner than a pixel and a half | `tfx_ribbon_bucket_globals_t` |
| `tfx_ribbon_cutout_fragment` | fragment | | Shades live and baked ribbons on a layer that writes depth, cut at the emitter's alpha threshold | `tfx_ribbon_bucket_globals_t` |
| `tfx_ribbon_lit_vertex` | vertex | | `tfx_ribbon_vertex` for a lit layer, adding the normal `tfx_ribbon_compute` wrote and the emitter's lighting | `tfx_ribbon_bucket_globals_t` |
| `tfx_ribbon_lit_fragment` | fragment | | Lights live and baked ribbons on a lit layer | `tfx_ribbon_bucket_globals_t` |
| `tfx_ribbon_lit_cutout_fragment` | fragment | | Lights live and baked ribbons on a lit layer that writes depth, cut at the alpha threshold | `tfx_ribbon_bucket_globals_t` |
| `tfx_sprite_playback` | compute | 128 | Interpolates baked sprite data into billboard instances | `tfx_sprite_data_push_t` |
| `tfx_billboard_baked_vertex` | vertex | | The billboard quads for the instances `tfx_sprite_playback` writes | `tfx_baked_billboard_push` |
| `tfx_billboard_baked_lit_vertex` | vertex | | `tfx_billboard_baked_vertex` for a lit layer | `tfx_baked_billboard_push` |
| `tfx_ribbon_playback` | compute | 128 | Interpolates baked ribbon data into the ribbons `tfx_ribbon_compute` builds | `tfx_ribbon_data_push_t` |
| `tfx_sprite_bounds` | compute | 128 | Optional: the bounding box of a baked frame's sprites | `tfx_sprite_data_push_t` |
| `tfx_ribbon_bounds` | compute | 128 | Optional: the bounding box of a baked frame's ribbons | `tfx_ribbon_bounds_push` |

Dispatch the compute entry points with enough groups to cover their count, `(count + group size - 1) / group size`. A thread past the end returns. The bounds passes are the exception: each thread bounds two items half a group apart, so they take half as many groups.

## Push constants

The library declares `tfx_gpu_particle_push_t`, `tfx_ribbon_bucket_globals_t`, `tfx_sprite_data_push_t` and `tfx_ribbon_data_push_t` in timelinefx.h, and fills the ribbon globals for you per bucket. Under the fixed bindings each of them is a uniform buffer instead, and the billboard draws take the draw constants in the slots table in place of theirs. The other push constants are yours to declare. Their members are `uint`s unless the table says otherwise:

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
| 7 | Uniform buffers: the view, with `view`, `projection` (both 4x4), `screen_size`, `milliseconds`, `timer_lerp` and `update_time`, which the lit fragments read with the light after it (see [Lit layers](#lit-layers)) |

The image data from `tfx_GetLibraryGPUShapes` holds each frame's uv rect and `texture_array_index`, which is the shape's bindless image index in the upper 16 bits and the frame's layer in the lower 16. The image a particle samples can change from one particle to the next within a draw, so it's read with `NonUniformResourceIndex`.

The two vertex inputs are vertex buffers:

| Draw | Buffer | Locations and formats |
| --- | --- | --- |
| Billboards | `tfx_instance_t`, instance rate | 0 `position` r32g32b32a32_sfloat, 1 `quaternion` r16g16b16a16_snorm, 2 `size` r16g16_sint, 3 `alignment` r8g8b8a8_snorm, 4 `intensity_gradient_map` r16g16_sint, 5 `curved_alpha_life` r8g8b8a8_unorm, 6 `indexes` r32_uint, 7 `captured_index` r32_uint |
| Ribbons | `tfx_ribbon_vertex_t`, vertex rate | 0 `position` r32g32b32_sfloat, 1 `segment_index` r32_uint, 2 `uv_offset_scale` r32g32_sfloat, 3 `ribbon_index` r32_uint, and for `tfx_ribbon_lit_vertex` only, 4 `normal` r32_uint |

These formats exist on every API. `size` and `intensity_gradient_map` are signed 16-bit integers the shader converts to float, because D3D, Metal and WebGPU have no scaled formats, and `alignment` and `curved_alpha_life` read four bytes because D3D and WebGPU have no three byte formats. The fourth byte is spare.

#### Bindless on Vulkan

Compile to SPIR-V with no defines, with the Vulkan row of the compile table. [implementations](../implementations) has a working host, `tfx_impl_vulkan`, and the editor draws with the same shaders.

- **One descriptor set.** Every storage buffer array, whatever its row type, is at binding 5, so one array of storage buffer descriptors serves them all, and a buffer needs one descriptor however many row types it's read as. If you add images or buffers while frames are in flight, give bindings 3 and 5 `UPDATE_AFTER_BIND`, `UPDATE_UNUSED_WHILE_PENDING` and `PARTIALLY_BOUND`.
- **Push constants** all start at offset 0. The largest block is `tfx_ribbon_bucket_globals_t` at 112 bytes, inside the 128 every device has, so one range for the vertex, fragment and compute stages serves every pipeline. The baked billboard draw pushes `tfx_baked_billboard_push`, and the fragment shader it shares with the live draw reads the first seven values as `tfx_billboard_push`.
- **Features:** Vulkan 1.2's `runtimeDescriptorArray` and `shaderSampledImageArrayNonUniformIndexing`, 1.1's `shaderDrawParameters` (see the notes under the compile table), and the update after bind features if you use them.
- **`tfx_billboard_baked_vertex` has no location 7**, because baked sprites don't read their previous instance. Leave `captured_index` out of its pipeline's vertex input, or the validation layer warns that it isn't consumed.

#### Bindless on D3D12

Compile to DXIL with no defines (`-target dxil -profile <vs_6_0, ps_6_0 or cs_6_0>`). [implementations](../implementations) has a working host, `tfx_impl_dx12`. HLSL can't put arrays of different types at one binding, so each descriptor array has register 0 of a space of its own:

| Register | Spaces | Holds |
| --- | --- | --- |
| `b0` | 1 | The view uniforms, a constant buffer view of 256 bytes |
| `b0` | 2 | The same view uniforms as the lit fragments read them, with the light, for a lit layer |
| `t0` | 1 | Images, `Texture2DArray` |
| `t0` | 2 to 25 | Read only structured buffers, one space per row type |
| `u0` | 1 to 8 | Read write structured buffers, one space per row type |
| `s0` | 1 | Samplers |

- **One descriptor table covers them all.** Make a range for each space, every one unbounded and starting at offset 0, so all of them see the whole shader visible heap and a descriptor's heap slot is the index the shaders take. Samplers get a table of their own. Unbounded constant buffer and read write ranges need resource binding tier 3.
- **The push constants are root constants** in space 0: `tfx_gpu_particle_push_t` at `b0`, `tfx_billboard_push` at `b1`, `tfx_baked_billboard_push` at `b2`, `tfx_ribbon_bucket_globals_t` at `b3`, `tfx_sprite_data_push_t` at `b4`, `tfx_ribbon_data_push_t` at `b5` and `tfx_ribbon_bounds_push` at `b6`. Together they're over the 64 DWORD limit, so make a root signature per block. The baked billboard draw takes two: its vertex shader reads `b2` and the fragment shader it shares with the live draw reads `b1`, so set both with the same first seven values. Set the kernel's block once per fetch and only `dispatch_index`, its first DWORD, per dispatch.
- **A structured view bakes in its stride**, so a buffer read as two row types needs a view of each. The kernel reads its pool through three (hot 32, cold 64 and noise 16 bytes) and the frame's upload through four (dispatches 128, the 64 byte tick rows that the tick, force, spawn block and run arrays all share, spawns 112 and locations 32). Every section starts on a multiple of its own stride, so a view over the whole buffer from element 0 works with the push constants' bases. The ribbon emitters are 416 bytes, the baked ribbons 64, their properties 76 and the animation instances 64.
- **Images** are read by the image index in the upper 16 bits of `texture_array_index`. With one atlas from `tfx_impl_GetUV` that index is 0, so put the atlas's view in heap slot 0.
- **The billboard draws** take the instances as a vertex buffer, so start each range with the draw's start instance. The input assembler adds it and the shader never reads `SV_InstanceID`.

The fixed bindings below also run on D3D12, under SDL_gpu.

### Fixed bindings: D3D11, SDL_gpu, OpenGL, WGSL and Metal

Compile with `-D TFX_FIXED_BINDINGS` for an API without descriptor arrays, or one that binds a fixed set of resources per pipeline stage. Every billboard, ribbon and playback entry point and the kernel compile this way; the bounds passes don't.

**It needs storage buffers in the vertex stage.** The billboard vertex shader reads the previous update's instance at `captured_index`, so the instances are a buffer it can index, not vertex attributes. That rules out WebGL2, most OpenGL ES 3.x devices (vertex stage storage buffers are optional there and often 0) and OpenGL below 4.3, which includes macOS.

What changes under the define:
- **One resource per slot.** The accessors ignore the push constants' buffer indexes. The kernel's bases still apply.
- **Every buffer is raw**, read with `Load<T>` at the C size of its rows, because some APIs only give raw views of a storage buffer (SDL_gpu on D3D12). On D3D11 create them with `D3D11_RESOURCE_MISC_BUFFER_ALLOW_RAW_VIEWS` and raw views.
- **No push constants on any target.** Each block is a uniform buffer.
- **One atlas.** Every shape comes from one atlas texture array, so `texture_array_index` is only the atlas layer and the uv rects are the shapes' places in it.
- **No vertex inputs.** The billboard draws read their instance from a buffer by `SV_InstanceID` and unpack it the way the vertex formats above would, and the ribbon draw reads its vertex by `SV_VertexID`. `SV_InstanceID` leaves out the draw's first instance on D3D11 and D3D12, and Slang subtracts `gl_BaseInstance` on Vulkan, so start every billboard draw at instance 0 and pass its first instance in `instance_offset`.
- **Samplers are separate** (a texture and a sampler at the same slot number) unless you define `TFX_COMBINED_SAMPLERS`, which gives SPIR-V combined image samplers. OpenGL is always combined. Under combined samplers the color ramps get a sampler too, though they're only loaded: bind a nearest one. Some layers want every texture paired with a sampler even when they're separate (sokol_gfx): give the color ramps a nearest sampler at slot 1, which the shader never reads.

#### Slots

Each entry point's resources in each class count up from 0. Slot n of a class is SPIR-V binding `TFX_<class>_BINDING + n` in set `TFX_<class>_SPACE`, and HLSL register `<letter>n` in space `TFX_<class>_SPACE`. The classes and their defaults:

| Class | Letter | Set and space | SPIR-V binding base |
| --- | --- | --- | --- |
| `VERTEX_READ`: vertex shader buffers | `t` | 0 | 0 |
| `VERTEX_UNIFORM`: vertex shader uniforms | `b` | 0 | 5 |
| `FRAGMENT_TEXTURE`: fragment shader textures | `t` | 0 | 7 |
| `FRAGMENT_SAMPLER`: fragment shader samplers, when they're separate | `s` | 0 | 9 |
| `FRAGMENT_UNIFORM`: fragment shader uniforms, only the lit fragments' view | `b` | 0 | 10 |
| `COMPUTE_WRITE`: compute read write buffers | `u` | 0 | 0 |
| `COMPUTE_READ`: compute read only buffers | `t` | 0 | 2 |
| `COMPUTE_UNIFORM`: compute uniforms | `b` | 0 | 6 |

The defaults put everything in set and space 0, with no binding shared by two classes that a vertex and fragment pair, or one compute pass, uses. That suits D3D11, OpenGL, WGSL and a Vulkan host with one descriptor set per pipeline. Override any of them with `-D`. For SDL_gpu, which numbers each stage's slots from 0 in its own set and space:

```
-D TFX_VERTEX_UNIFORM_SPACE=1 -D TFX_VERTEX_UNIFORM_BINDING=0
-D TFX_FRAGMENT_TEXTURE_SPACE=2 -D TFX_FRAGMENT_TEXTURE_BINDING=0 -D TFX_FRAGMENT_SAMPLER_SPACE=2 -D TFX_FRAGMENT_SAMPLER_BINDING=0
-D TFX_FRAGMENT_UNIFORM_SPACE=3 -D TFX_FRAGMENT_UNIFORM_BINDING=0
-D TFX_COMPUTE_READ_BINDING=0 -D TFX_COMPUTE_WRITE_SPACE=1 -D TFX_COMPUTE_UNIFORM_SPACE=2 -D TFX_COMPUTE_UNIFORM_BINDING=0
-D TFX_COMBINED_SAMPLERS
```

Each entry point's slots, in the order SDL_gpu wants (textures, then read only buffers, then read write buffers, then uniforms). Raw buffers are rows in the C layout:

| Entry point | Slot | Holds |
| --- | --- | --- |
| `simulate_particles` | `COMPUTE_READ` 0 | The graph lookups |
| | `COMPUTE_READ` 1 | The frame's upload from `tfx_CopyGPUParticleUpload`, dispatches at its start |
| | `COMPUTE_WRITE` 0 | The pool: the hot section at its start, cold and noise at `cold_base` / `noise_base` |
| | `COMPUTE_WRITE` 1 | The instance buffer |
| | `COMPUTE_UNIFORM` 0 | `tfx_gpu_particle_push_t` |
| `tfx_billboard_vertex` | `VERTEX_READ` 0 | The instances being drawn |
| | `VERTEX_READ` 1 | The previous update's instances |
| | `VERTEX_READ` 2 | The image data |
| | `VERTEX_READ` 3 | The emitter properties |
| | `VERTEX_UNIFORM` 0 | The view: `view`, `projection` (both 4x4), `screen_size`, `milliseconds`, `timer_lerp`, `update_time` |
| | `VERTEX_UNIFORM` 1 | The draw: `instance_offset`, `index_offset` (uints), `ndc_offset_x`, `ndc_offset_y` (floats) |
| `tfx_billboard_baked_vertex` | `VERTEX_READ` 0 | The instances `tfx_sprite_playback` wrote |
| | `VERTEX_READ` 1 | The image data |
| | `VERTEX_READ` 2 | The emitter properties |
| | `VERTEX_UNIFORM` 0, 1 | The view and the draw, as above |
| `tfx_billboard_fragment`, `tfx_billboard_cutout_fragment`, `tfx_ribbon_fragment`, `tfx_ribbon_cutout_fragment` | `FRAGMENT_TEXTURE` 0 | The shape atlas, sampled |
| | `FRAGMENT_TEXTURE` 1 | The color ramps, a texture array read with an integer load |
| | `FRAGMENT_SAMPLER` 0 | The atlas's sampler, when samplers are separate |
| `tfx_billboard_lit_fragment`, `tfx_billboard_lit_cutout_fragment`, `tfx_ribbon_lit_fragment`, `tfx_ribbon_lit_cutout_fragment` | `FRAGMENT_TEXTURE` 0, 1, `FRAGMENT_SAMPLER` 0 | As the unlit fragments |
| | `FRAGMENT_UNIFORM` 0 | The view with the light, `tfxVIEW_UNIFORMS_SIZE` (224) bytes |
| `tfx_billboard_lit_vertex`, `tfx_billboard_baked_lit_vertex`, `tfx_ribbon_lit_vertex` | | The slots of `tfx_billboard_vertex`, `tfx_billboard_baked_vertex` and `tfx_ribbon_vertex` |
| `tfx_ribbon_compute` | `COMPUTE_READ` 0 | The ribbon emitters |
| | `COMPUTE_READ` 1 | The graph lookups |
| | `COMPUTE_READ` 2 | The ribbons |
| | `COMPUTE_READ` 3 | The ribbon segments |
| | `COMPUTE_WRITE` 0 | The ribbon vertexes it writes |
| | `COMPUTE_WRITE` 1 | The ribbon indexes it writes |
| | `COMPUTE_UNIFORM` 0 | `tfx_ribbon_bucket_globals_t` |
| `tfx_ribbon_vertex` | `VERTEX_READ` 0 | The ribbon vertexes `tfx_ribbon_compute` wrote |
| | `VERTEX_READ` 1 | The ribbons |
| | `VERTEX_READ` 2 | The ribbon emitters |
| | `VERTEX_READ` 3 | The graph lookups |
| | `VERTEX_READ` 4 | The image data |
| | `VERTEX_UNIFORM` 0 | The view |
| | `VERTEX_UNIFORM` 1 | `tfx_ribbon_bucket_globals_t` |
| `tfx_sprite_playback` | `COMPUTE_READ` 0 | The animation offsets |
| | `COMPUTE_READ` 1 | The animation instances |
| | `COMPUTE_READ` 2 | The baked sprites |
| | `COMPUTE_READ` 3 | The baked emitter properties |
| | `COMPUTE_WRITE` 0 | The billboard instances it writes |
| | `COMPUTE_UNIFORM` 0 | `tfx_sprite_data_push_t` |
| `tfx_ribbon_playback` | `COMPUTE_READ` 0 | The ribbon offsets |
| | `COMPUTE_READ` 1 | The animation instances |
| | `COMPUTE_READ` 2 | The baked ribbons |
| | `COMPUTE_READ` 3 | The baked ribbon properties |
| | `COMPUTE_WRITE` 0 | The ribbons it writes |
| | `COMPUTE_UNIFORM` 0 | `tfx_ribbon_data_push_t` |

With the defaults, `tfx_billboard_vertex` is t0 to t3 and b0, b1 in HLSL, and SPIR-V bindings 0 to 3, 5 and 6 in set 0. With the SDL_gpu line its uniforms are b0, b1 in space1, and bindings 0, 1 in set 1.

- **D3D11:** [implementations](../implementations) has a working host with the whole build in its CMake. Notes from it:
  - The kernel's upload can be a `DYNAMIC` raw buffer that `tfx_CopyGPUParticleUpload` fills straight through `Map(WRITE_DISCARD)`.
  - Create the pool zeroed. The library never clears it and it never needs clearing after that.
  - Update the kernel's constant buffer before each dispatch, changing only `dispatch_index`.
  - Unbind what the next pass writes. The instance buffer the kernel writes is the one the vertex shader reads, and a buffer bound for reading and writing at once gets one of its binds cleared, with a debug layer warning.
  - If you keep two instance buffers for interpolation, copy both across when they grow, because the GPU instances in the previous one exist nowhere else.
  - fxc warns that 256 threads times the kernel's temp registers is over its recommended total (X4714), and that `pow` won't work for a negative base (X3571). The base is a unit hash, never negative, and neither shows up at run time.
  - fxc also warns of a potentially uninitialized variable (X4000) in `tfx_generate_shape_offset`, `tfx_generate_stepped_direction`, `tfx_grid_march_length`, `tfx_grid_march_point` and `tfx_particle_quaternion`. Every path through them returns a value: Slang turns their early returns into a result variable that fxc can't follow.
  - `simulate_particles` takes about 10 seconds through fxc, so a clean build of the D3D11 shaders takes a while.
- **SDL_gpu:** compile with the line above, SPIR-V for Vulkan and DXIL for D3D12. Push each uniform block with `SDL_PushGPUVertexUniformData`, `SDL_PushGPUComputeUniformData` and so on at its slot, the kernel's once per dispatch. SDL doesn't zero buffers, so create the pool by uploading zeros. [implementations](../implementations) has a working host, `tfx_impl_sdlgpu`, with the whole build in its CMake and more notes in its README.
- **OpenGL:** compile through SPIR-V with the line in the compile table. Every entry point comes out as `#version 430`, so it needs OpenGL 4.3. [implementations](../implementations) has a working host, `tfx_impl_gl`, with the whole build in its CMake, and it has also run under sokol_gfx. With the defaults, buffers are SSBOs and uniform blocks std140 UBOs at the bindings in the slots table (`glBindBufferBase`), and the two textures are `sampler2DArray`s at bindings 7 and 8.
  - The vertex shaders use `gl_BaseVertexARB` and `gl_BaseInstanceARB` when the driver has `GL_ARB_shader_draw_parameters`, and otherwise a plain `uniform int SPIRV_Cross_BaseVertex` and `SPIRV_Cross_BaseInstance`, which are 0 unless you set them. Either way, draw at base vertex 0 and first instance 0, and pass the first instance in `instance_offset`. `--no-support-nonzero-baseinstance` stops SPIRV-Cross from adding the base instance a second time. glslangValidator rejects `gl_BaseVertexARB`, but drivers accept it.
  - Each uniform block is rows of 16 bytes holding 32-bit words, uints and floats mixed. A layer that sets uniforms member by member with `glUniform*` and has no uint type (sokol_gfx's GL backend) can't take them as SPIRV-Cross writes them, and `--flatten-ubo` refuses blocks that mix types. Add `--glsl-emit-ubo-as-plain-uniforms`, rewrite each `uniform Block name;` as `uniform ivec4 name_raw[N];`, and rebuild the struct at the top of `main()` with `intBitsToFloat` and `uint()`. Then feed it as N `ivec4`s.
- **WGSL:** compile with the WebGPU line in the compile table. The kernel always uses these bindings, define or not, because WebGPU has no descriptor arrays or push constants and allows 8 storage buffers per stage by default. With the defaults every entry point's resources are in group 0 at the bindings in the slots table, with the atlas's sampler separate. [implementations](../implementations) has a working host, `tfx_impl_webgpu`, run on wgpu-native's Vulkan and D3D12 backends.
  - Each uniform block is a uniform buffer bound at a dynamic offset. The kernel needs its own `dispatch_index` per dispatch, so write one 64 byte block per dispatch, and one per draw range and ribbon bucket, at 256 byte steps.
  - The vertex shaders only read storage buffers, which is all WebGPU allows in the vertex stage. The ribbon compute uses the most storage buffers, 6.
  - Bind the color ramps as `unfilterable-float`. They're only loaded, so a 32 bit float ramp works without `float32-filterable`.
  - Make every texture view a 2D array, even of a one layer texture, whose default view is a plain 2D one.
  - In your own fragment shaders, pass a varyings struct only to `[ForceInline]` functions. Slang's WGSL turns any other into a copy of the struct that keeps `@interpolate(flat)` without a `@location`, and naga rejects the shader with "input/output binding is not consistent". `tfx_shade_billboard_parts` and `tfx_shade_ribbon_parts` are inlined.
- **Metal:** compile with the Metal line in the compile table. Slang makes HLSL registers `tN`, `uN` and `bN` all `[[buffer(N)]]` and ignores spaces, so the classes a stage uses would collide. `-D TFX_METAL_SLOTS` gives each buffer class its own tens on Metal and changes nothing on any other target. [implementations](../implementations) has a working host, `tfx_impl_metal`.

  | Class | Metal slots |
  | --- | --- |
  | `VERTEX_READ`, `COMPUTE_WRITE` | `[[buffer(n)]]` |
  | `VERTEX_UNIFORM`, `COMPUTE_READ` | `[[buffer(10 + n)]]` |
  | `COMPUTE_UNIFORM` | `[[buffer(20 + n)]]` |
  | `FRAGMENT_TEXTURE`, `FRAGMENT_SAMPLER` | `[[texture(n)]]`, `[[sampler(n)]]` |
  | `FRAGMENT_UNIFORM` | the fragment stage's `[[buffer(n)]]` |

  - Each uniform block is small enough for `setVertexBytes` and `setBytes`, which take up to 4 KB, so no uniform buffers are needed.
  - The buffers come out as `device` pointers, though the vertex shaders only read them, so Metal warns of writable resources in a vertex function when it compiles them. The warning is harmless.
  - The color ramps are read with `read()`, so they take no sampler.

### Metal

The kernel has its own path, chosen when Slang targets Metal. Each buffer index in the push constants is a slot in `tfx_metal_buffers`, an argument buffer of `TFX_METAL_BUFFER_COUNT` (default 16, override with `-D`) 64 bit `gpuAddress` values. The accessors read and write every row through it as raw bytes packed exactly like the C structs, because a Metal `float3` is 16 bytes rather than 12. The host:
- fills the table with the addresses and points each push constant index at its slot. A section can be its own slot with its address offset and a base of 0, or share one address with the bases. The table is 128 bytes, so it can go with `setBytes` rather than in an `MTLBuffer`.
- binds the push constants and the table at the `[[buffer(n)]]` slots Slang gives them. `simulate_particles` compiled without defines uses 0 and 1. Under `TFX_FIXED_BINDINGS` the unused fixed buffers take slots first and move them, so compile the kernel without the defines. A host entry point's push constants take the same slot when declared with `TFX_BINDLESS_REGISTER(b0, space0)`, see [Adding your own behaviour](#adding-your-own-behaviour).
- calls `useResource` on every buffer the table points to: read for the frame's upload and the graph lookups, read and write for the pool and the instance buffer. Metal's shader validation checks these too, so a buffer left out is reported.

`gpuAddress` needs Metal 3 (macOS 13).

The render passes use the fixed bindings with `-D TFX_METAL_SLOTS`, the Metal line in the compile table, at the slots in the Metal notes under the fixed bindings.

## Billboard particles

Draw 6 vertices per instance with no index buffer: `tfx_billboard_vertex` makes the quad's corners from `SV_VertexID`. Each layer's instances from `tfx_GetInstanceBufferByLayer` (or the key ranges from `tfx_GetNextInstanceKeyRange`) are one instanced draw. Besides the instance, the vertex shader reads:
* the previous update's instance buffer, so that it can interpolate the position, size and rotation of each sprite using its `captured_index`. A sprite without the capture flag set (the first frame of a particle, or one that has just looped back to the start of a path) isn't drawn that frame.
* the emitter properties buffer from `tfx_GetParticlePropertiesBuffer` for the image handle, the billboard alignment, the color ramp row, the start frame of the emitter's image and its alpha threshold.
* the image data buffer for the uvs and the image index.
* the view, with `timer_lerp`, how far the current frame is between the previous update and the current one. Get it from `tfx_GetStageInterpolation`, not straight from your timer.

It's an all-in-one shader that handles every billboard alignment mode (camera facing, align to a vector, and camera facing that rolls to a vector) plus the stretch along the direction of travel, so an effect whose emitters mix alignment modes still draws in a single draw call.

`tfx_billboard_fragment` samples the particle image and then uses the red channel (scaled by the sprite's gradient map value) as a lookup into the emitter's row of the color ramp texture. The color ramps are plotted by the library from the emitter's color graphs, with the heat response already baked in, as RGBA16F, and you upload them from `tfx_GetColorRampBitmap`. The ramp color is multiplied by the sprite's intensity, and the curved alpha / sharpness terms do the dissolve fade. It writes premultiplied alpha, so blend with ONE, ONE_MINUS_SRC_ALPHA.

`tfx_billboard_cutout_fragment` shades the same way, but treats the shape as a mask: a pixel whose final alpha is below the emitter's `alpha_threshold` is discarded, and what's left is written opaque. Use it for the layers whose `tfx_GetLayerFlags` has `tfxLayerFlags_writes_depth`, which always includes layer 1, Depth. Give its pipeline blending off and depth test and write on, and draw those layers first, in index order, then the other layers with `tfx_billboard_fragment`, so the blended particles behind them are hidden. Draw both against your scene's depth buffer. The flags belong to the library, and a refresh can change them, so read them when you draw. It takes the same vertex shaders, push constants and slots as `tfx_billboard_fragment`. Both read their varyings through `tfx_shade_billboard_parts`, which a fragment shader of your own can call too, see [Your own fragment shaders](#your-own-fragment-shaders).

The alpha threshold reaches the fragment in bits 16 to 23 of the `indexes.y` varying, as an 8 bit unorm, so read the color ramp layer as `(indexes.y >> 8) & 0xFF` in a fragment shader of your own.

## Ribbons

Ribbons take an extra step compared to billboards. The library writes out ribbon, segment and ribbon emitter data, and `tfx_ribbon_compute` turns that into actual geometry: one thread per segment, sampling the width/scale/clip/morph/noise/lag graphs, building the ribbon frame (facing the camera, at a fixed angle, or following the curve of the path) and writing out the tessellated triangle strip. The segment index is packed into each vertex so the later stages can recover which ribbon, which segment and which side of the strip a vertex came from.

Ribbons are grouped into buckets, one for each segment count, effect sort key and layer, and each bucket is one compute dispatch and one indexed draw. After `tfx_UpdateStage`, copy the ribbon data into your staging buffers with `tfx_CopyRibbonDataToStagingBuffers` and loop over the buckets with `tfx_CreateRibbonDispatch` and `tfx_NextRibbonDispatch`, which give you the offsets, counts and push constants for each dispatch and draw. To draw the ribbons of several stages with one upload, one compute pass and one draw pass, use a ribbon batch instead: `tfx_BeginRibbonBatch`, `tfx_AddStageToRibbonBatch` for each stage, `tfx_GetRibbonBatchSizes` and `tfx_CreateRibbonBatchDispatch`. Both flows are described in full in timelinefx.h.

`tfx_ribbon_vertex` then transforms those vertices, unpacks that segment index, samples the over-length shading graphs (intensity, gradient map, curved alpha, sharpness) and ping-pongs the U coordinate so a texture can wrap repeatedly along the ribbon without a flipped seam at each wrap point.

`tfx_ribbon_fragment` shades a ribbon the same way as a billboard, and fades it out as it gets thinner than about 1.5 pixels. Below that the rasteriser starts dropping fragments and a thin ribbon breaks up into dashes, so fading it gives a tail that thins out smoothly instead.

`tfx_ribbon_cutout_fragment` is the ribbon version of `tfx_billboard_cutout_fragment`, for the buckets whose `dispatch.layer` writes depth. It discards a pixel whose alpha, times the thin ribbon fade, is below the emitter's alpha threshold, and writes the rest opaque, so a tapering ribbon ends cleanly rather than breaking into dashes. Give its pipeline the ribbon pipeline's vertex input with blending off and depth test and write on, and draw those buckets after the depth layers' billboards and before every blended layer. Live and baked ribbons share it, as they share `tfx_ribbon_fragment`. Both read their varyings through `tfx_shade_ribbon_parts`, which a fragment shader of your own can call too, see [Your own fragment shaders](#your-own-fragment-shaders). The threshold reaches the fragment in the `texture_indexes.z` varying, 0 to 255.

## Lit layers

A layer whose `tfx_GetLayerFlags` has `tfxLayerFlags_lit` draws with the lit entry points: `tfx_billboard_lit_vertex` (live) or `tfx_billboard_baked_lit_vertex` (baked) with `tfx_billboard_lit_fragment`, and `tfx_ribbon_lit_vertex` with `tfx_ribbon_lit_fragment`, or the lit cutouts when the layer writes depth too. Pipeline state, draw order, push constants and the vertex shaders' inputs are the same as the unlit ones', except that the lit ribbon vertex shader also reads `normal` at location 4. The unlit entry points don't change, so a host that never draws a lit layer, or draws one with the unlit pipelines, needs nothing new.

The color the unlit fragments output becomes the albedo, which the fragment lights with wrapped Lambert from one directional light plus a hemisphere ambient: `sky_color` on normals facing +Y and `ground_color` on normals facing down. Color above 1, from heat response or intensity, is added unlit, so a hot particle still glows on its shaded side. Alpha is unchanged.

The light is in the view uniforms, after the fields every vertex shader reads, as `tfx_view_light_t` at `tfxVIEW_LIGHT_OFFSET` (160):

| Offset | Field |
| --- | --- |
| 148 | 12 bytes of padding |
| 160 | `direction`: unit length, towards the light. w unused |
| 176 | `color`: the light's color times its intensity. w unused |
| 192 | `sky_color` |
| 208 | `ground_color` |

The lit fragments read the uniforms as `tfx_lit_view_uniforms`, `tfxVIEW_UNIFORMS_SIZE` (224) bytes, so bind at least that much to them. Fill the light every frame: a zeroed light draws lit particles black. The pads are three floats, not a float3, because std140 and WGSL would put a float3 at 160 where D3D puts it at 148. The unlit shaders still declare the 148 byte `tfx_view_uniforms`, so a 148 byte buffer stays valid for them.

Each emitter's lighting word is `[wrap unorm8 << 16 | roundness unorm8 << 8 | normal mode]`. Billboards keep it in the emitter properties row's `lighting`, and ribbons above the angle type in the ribbon emitter row's `angle_type` (and `shader_type` for baked ribbons), so the angle type is `angle_type & 0xFF`. In a ribbon's word the library also sets bit 7 when the emitter's layer isn't lit, and `tfx_ribbon_compute` then skips the normals. Live ribbons read the layer's flags every update, and baked ribbons keep the flags their layer had when they were recorded.

- **Billboards:** the vertex shader passes the quad's world space axes, with the normal turned towards the camera. Roundness bends that normal into a sphere across the quad, so at roundness 0 a camera facing billboard lights evenly.
- **Ribbons:** `tfx_ribbon_compute` writes an octahedral unit normal into each vertex's `normal` on a lit layer, and 0 on the others. Radial (mode 0) points out from the ribbon's origin, so a ribbon that wraps around it is lit on the side facing the light. Surface (mode 1) is the ribbon's face, turned towards `camera_position`. Roundness bends both towards the edges, so the ribbon shades like a tube.

## Your own fragment shaders

A layer can mean whatever your game wants. Name it in the editor, match it in your renderer by index or `tfx_GetLayerName`, and draw its range with a fragment shader of your own. The vertex shader, push constants and slots stay the library's, and your fragment starts from the same shading the built-in ones use. The editor previews the layer by its flags, so set Writes Depth and Lit to the nearest look. This one draws a billboard as `tfx_billboard_fragment` does, and also writes the color that heat response and intensity push past 1 to a second target, for a bloom pass that only particles feed:

```slang
#include "timelinefx.slang"

struct my_glow_targets {
	float4 color : SV_Target0;
	float4 glow : SV_Target1;
};

[shader("fragment")]
my_glow_targets my_glow_billboard_fragment(tfx_billboard_output input) {
	const tfx_particle_shade shade = tfx_shade_billboard_parts(input);
	my_glow_targets targets;
	targets.color = tfx_premultiply(shade);
	//Alpha 0, so the color target's ONE, ONE_MINUS_SRC_ALPHA blend adds it
	targets.glow = float4(max(shade.color - 1.0, 0.0) * shade.coverage, 0.0);
	return targets;
}
```

Compile it with your target's line from the compile table, adding `-I` with this folder, with `-entry my_glow_billboard_fragment` and `-D TFX_ENTRY_tfx_billboard_fragment`: under the fixed bindings, the entry point named in `TFX_ENTRY_` gives yours its slots. Its pipeline is the `tfx_billboard_fragment` one with a second color target, blended the same way.

- **`#include`, not `import`.** Slang checks every entry point of an imported module against the one profile, which fails on DXIL and DXBC, and the library's macros, `TFX_FIXED_DECLARE` among them, don't reach a file that imports it.
- **Each draw has its own varyings and helpers:**

  | Your fragment takes | From | Shade with | `TFX_ENTRY_` |
  | --- | --- | --- | --- |
  | `tfx_billboard_output` | `tfx_billboard_vertex`, `tfx_billboard_baked_vertex` | `tfx_shade_billboard_parts(input)` | `tfx_billboard_fragment` |
  | `tfx_lit_billboard_output` | `tfx_billboard_lit_vertex`, `tfx_billboard_baked_lit_vertex` | `tfx_shade_billboard_parts(tfx_unlit_billboard(input))`, then `tfx_light_particle(tfx_lit_billboard_view(), shade.color, tfx_lit_billboard_normal(input), input.lighting)` | `tfx_billboard_lit_fragment` |
  | `tfx_ribbon_output` | `tfx_ribbon_vertex` | `tfx_shade_ribbon_parts(input)`, times `tfx_thin_ribbon_fade(input.ribbon_across)` | `tfx_ribbon_fragment` |
  | `tfx_lit_ribbon_output` | `tfx_ribbon_lit_vertex` | `tfx_shade_ribbon_parts(tfx_unlit_ribbon(input))`, then `tfx_light_particle(tfx_lit_ribbon_view(), shade.color, tfx_safe_normalize(input.normal), input.lighting)`, times the fade | `tfx_ribbon_lit_fragment` |

  The lit billboard varyings also carry the quad's world space `tangent`, `bitangent` and `normal`, the frame a normal map needs.
- **`tfx_particle_shade`** has `color`, the image times the ramp and intensity, which can pass 1, `coverage`, the image's alpha faded by curved alpha, and `alpha`, the coverage times the ramp's alpha. `tfx_premultiply` gives what the built-in fragments write. On a layer that writes depth, discard below the alpha threshold as the cutouts do: `float(input.indexes.y >> 16) / 255.0` for billboards and `float(input.texture_indexes.z) / 255.0` for ribbons.
- **WGSL:** pass a varyings struct only to `[ForceInline]` functions, see the WGSL notes under the fixed bindings. The helpers above that take varyings are inlined.
- **Your own resources** need bindings the library's don't use:
  - Fixed bindings on D3D11, Metal and SDL_gpu: `TFX_FIXED_DECLARE(Texture2D, my_mask, FRAGMENT_TEXTURE, t, 2);` takes the fragment texture slot after the library's two, and follows the class's `-D` overrides.
  - Fixed bindings with the default SPIR-V and WGSL bases: that slot's binding, 9, is the sampler's, so give yours one above 10, `[[vk::binding(11, 0)]] Texture2D my_mask : register(t2);`. On OpenGL make it a `Sampler2D`, because spirv-cross drops the binding of a texture it has to pair with a dummy sampler.
  - Bindless: set 1 on Vulkan and a space above 25 on D3D12, `[[vk::binding(0, 1)]] Texture2D my_mask : register(t0, space26);`.

## GPU simulated particles

`simulate_particles` simulates `run_on_gpu` emitters and writes their instances, with its own graph sampling, noise fields, emitter forces and spawn generation.

### Dispatching

On a stage with a `gpu_particle_capacity`, `tfx_GetGPUParticleFrame` hands you dispatch records, tick records, the rows spawned since the last fetch and a table of spawn locations that relative particles follow. Dispatch the kernel once per record with `tfx_gpu_particle_push_t` as push constants and `(slot_count + 255) / 256` groups. Each thread owns one ring slot: it takes in any spawn that lands in its slot, runs every tick of the fetch, keeps the particle state in a pool you size once from `tfx_GetGPUParticleBufferSizes`, and writes its instance into the same instance buffer as the CPU particles, in the range after them. The billboard shaders draw them unchanged.

Point, User, Area, Ellipse, Disc, Cylinder and Line emitters have their spawns generated in the kernel from a few rows per tick and location, so a burst of any size costs the CPU next to nothing. Icosphere, Path, Other Emitter and Spawn on Ribbon emitters spawn on the CPU and upload a row per particle. Every emitter force runs in the kernel from force rows the library resolves each tick, so they need nothing more from you.

The pool has hot, cold and noise sections, and each section starts on a 256 byte boundary. They can share one buffer or have one each, and the same goes for the frame's arrays: `upload_bytes` and the frame's upload offsets lay them all out in one buffer, which `tfx_CopyGPUParticleUpload` fills. Every buffer in the push constants has its own index plus, for those that can share, a `*_base` row the section starts at. On Vulkan and Metal, where storage buffers are untyped, one buffer and one index with the bases is enough. On D3D12 a structured view carries its row size, so give each section a view of its own, either over the section with a base of 0, or over the whole buffer with the base. The fixed binding targets always use one pool buffer and one upload buffer with the bases.

### Adding your own behaviour

You can add custom forces, collision, spawn placement or changes to how particles are drawn without changing `timelinefx.slang`. `#include` it in a Slang file of your own, not `import`, see [Your own fragment shaders](#your-own-fragment-shaders) for why, and implement any of the three hook interfaces:

| Interface | What it lets you do |
| --- | --- |
| `ITfxParticleUpdate` | `acceleration` adds to the particle's acceleration every tick like the built-in forces, and `after_move` runs once the particle has moved, where you can also kill it |
| `ITfxParticleSpawn` | `after_spawn` places the particles of User emitters that the kernel generates, the GPU side of `tfx_SetStageSpawnCallback` |
| `ITfxParticleWrite` | `before_write` changes a live particle's position, rotation, size, color ramp position, alpha or image frame just before its instance is written |

Then write your own `[numthreads(tfx_group_size, 1, 1)]` entry point that calls `tfx_simulate` with your hooks. `tfx_simulate(thread_index, push, update)` takes an update hook, and a fourth and fifth argument add a spawn hook and a write hook, with `tfx_no_particle_update`, `tfx_no_particle_spawn` and `tfx_no_particle_write` for the ones you don't need. The GPU particle simulation section of timelinefx.h walks through each hook's arguments. Your hooks only run on the GPU, so emitters that fall back to the CPU (during warmup, for example) simulate without them. This one bounces particles off a floor whose height and bounce come from its own push constants:

```slang
#include "timelinefx.slang"

struct my_push {
	tfx_gpu_particle_push tfx;
	float ground_height;
	float bounce;
};

[[vk::push_constant]] ConstantBuffer<my_push> my_push_constants TFX_BINDLESS_REGISTER(b0, space0);
TFX_UNIFORM_PUSH(my_push, my_uniform_push);

my_push my_entry_push() {
	__target_switch {
	TFX_UNIFORM_PUSH_CASES return my_uniform_push;
	default: return my_push_constants;
	}
}

//In world space, so the particles of relative rings land on the floor too
struct my_floor_update : ITfxParticleUpdate {
	float ground_height;
	float bounce;
	float3 acceleration(tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, tfx_particle_context context) {
		return float3(0.0, 0.0, 0.0);
	}
	void after_move(inout tfx_gpu_particle_hot hot, tfx_gpu_particle_cold cold, tfx_particle_context context) {
		const tfx_ring_transform transform = context.ring_transform(cold);
		float3 world = transform.to_world(hot.position);
		if (world.y >= ground_height) return;
		world.y = ground_height;
		float3 velocity = transform.direction_to_world(hot.velocity);
		velocity.y = -velocity.y * bounce;
		hot.position = transform.to_ring(world);
		hot.velocity = transform.direction_to_ring(velocity);
	}
};

[shader("compute")]
[numthreads(tfx_group_size, 1, 1)]
void my_simulate_particles(uint3 thread_id : SV_DispatchThreadID) {
	const my_push push = my_entry_push();
	my_floor_update update;
	update.ground_height = push.ground_height;
	update.bounce = push.bounce;
	tfx_simulate(thread_id.x, push.tfx, update);
}
```

- **Compile it** with the kernel's line from the compile table for your target, adding `-I` with this folder and `-entry my_simulate_particles`. It needs no `TFX_ENTRY_`, and on Metal no defines, as `simulate_particles` doesn't.
- **The push constants are declared twice**, as `simulate_particles`' are: the push constant block, and `TFX_UNIFORM_PUSH`'s uniform twin at the kernel's uniform slot for WGSL and the fixed bindings, with `__target_switch` picking between them. `tfx_gpu_particle_push` comes first, so your own members start at offset 64 on every target.
- **`TFX_BINDLESS_REGISTER(b0, space0)`** gives the block the register `simulate_particles`' has on D3D12, `b0`, and its slot on Metal, buffer 0, with the buffer table still at 1. Without it Slang gives your block the next free ones, `b7` and buffer 12.
- **Dispatch it as you dispatch `simulate_particles`**, with the same buffers and slots. Only the block's size changes: on D3D12 its root constants, and where it's a uniform buffer, each dispatch's block.

## Pre-baked effects

For maximum speed you can pre-bake an effect to sprite data and update the instance buffer on the GPU instead of simulating on the CPU.

Every frame, after `tfx_UpdateAnimationManager`, dispatch `tfx_sprite_playback` with `tfx_sprite_data_push_t`. It interpolates each baked sprite between its frame and the one before, resolves the image frame (the emitter's start frame included) and bakes the stretch into `position.w`, then writes a billboard instance. Draw those with `tfx_billboard_baked_vertex` and `tfx_billboard_fragment`. The output isn't grouped by layer across instances, so baked layers that write depth draw like Default. The baked vertex shader reads each instance as it is, with no previous buffer and no timer lerp: `tfx_billboard_vertex` would add the start frame a second time and see no motion to stretch along.

Baked ribbons go through the same `tfx_ribbon_compute` and ribbon draw as live ones. Dispatch `tfx_ribbon_playback` with `tfx_ribbon_data_push_t`, then walk the animation manager's buckets with `tfx_CreateAnimationRibbonDispatch` and `tfx_NextAnimationRibbonDispatch` exactly as you walk a stage's: one `tfx_ribbon_compute` dispatch and one indexed draw per bucket. `tfx_CopyAnimationRibbonEmitters` gives you the emitter buffer `tfx_ribbon_compute` reads for baked ribbons, and `tfx_GetAnimationRibbonBufferSizes` the buffer sizes for the frame. Baked ribbons carry their own graph rows, so bind `tfx_GetAnimationRibbonGraphLookupsBufferPointer` in place of `tfx_GetGPUGraphLookupsBuffer` for their compute and draw. That way sprite data plays back without a library loaded.

`tfx_sprite_bounds` and `tfx_ribbon_bounds` are optional extras rather than part of the render path. Each group reduces its sprites, or its ribbons' segments, to one box, and you reduce the groups' boxes on the CPU. The editor uses them to auto-fit the camera and to size sprite sheet exports; in a game they're useful for culling a baked effect.

[See here for an example of how to use the library to draw pre-baked effect data](https://github.com/peterigz/zest/tree/main/examples/SDL2/zest-timelinefx-particle-techniques)
