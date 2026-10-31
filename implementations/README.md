# TimelineFX implementation examples

Minimal examples of drawing TimelineFX effects with different graphics API backends. Each API has one backend `.h`/`.cpp` pair, and every backend shares one renderer agnostic pair:

| File | What it holds |
| --- | --- |
| `tfx_impl_common.h/.cpp` | Everything that would be written the same way whatever API draws the particles: shape loading, the shape atlas, the fixed step timer, the camera matrices and gathering each frame's instances from a stage. Nothing in it talks to a GPU. |
| `tfx_impl_dx11.h/.cpp` | The Direct3D 11 backend. |
| `../Shaders/timelinefx.slang` | Every shader: the GPU particle kernel and the render passes. CMake has slangc compile each entry point the example uses to a `.cso`. |
| `example_scene.h/.cpp` | The demo scene every example plays. It's example code showing how a game uses the library, separate from the integration code above. |
| `example.tfx` | The library the scene plays. |
| `example.tfxsd` | Burst, recorded from `example.tfx` in the editor, for baked playback. |
| `example_win32_dx11/main.cpp` | A Win32 window, a D3D11 device and swap chain, and the main loop that calls the scene and the backend. |
| `third_party/` | stb_image (PNG shapes), stb_rect_pack (the atlas), bcdec (BC4/5/7 KTX2 shapes) and zstd's single file decompressor (supercompressed KTX2 shapes). |

So there are three layers. `tfx_impl_common` is integration code you can copy into a game as it is, `example_scene` is
what this demo does with the library, and `main.cpp` is only the window, the device and the calls into the other two.
A new graphics API only needs its own `tfx_impl_<api>` pair and a `main.cpp`. The test for where code goes is whether
every backend would write it the same way: if it would, it belongs in common.

The code is C-style C++: plain structs and free functions, no STL. The headers are wrapped in `extern "C"`, so a C
project can call a backend that was built as C++. The examples only include `timelinefx.h`, never
`timelinefx_internal.h`, which is exactly what a game sees.

Status: D3D11 draws CPU simulated billboards and ribbons, GPU simulated particles, and baked billboard and ribbon
playback. Every other API is still to come, see [Porting to a new API](#porting-to-a-new-api).

## Building and running

You can use cmake to build the examples:

```
cmake -S TimelineFXLib/implementations -B build-examples
cmake --build build-examples --config Debug
```

Every shader is compiled from `../Shaders/timelinefx.slang` with slangc, which comes with the Vulkan SDK and with
the Slang releases on GitHub. CMake looks on the `PATH` and in the Vulkan SDK, and stops if it can't find it, so put
its bin folder on the `PATH` or set `TFX_SLANGC` to it.

```
example_win32_dx11 [updates per second]
```

Space fires a burst of sparks, O pauses the camera orbit, the mouse wheel zooms and Escape quits. The window title shows
the frame rate and the particle counts. In debug builds the D3D11 debug layer's messages are repeated on the console.
Passing an update rate (20 or 400, say) is a quick way to check the interpolation: the effects should look the same at
any rate.

## The scene

| Effect | What it shows |
| --- | --- |
| Spitting Torch | An ordinary effect moved with `tfx_SetEffectPosition` every frame. Animated flames, smoke and embers trail behind it |
| Burst Sparks | Fire and forget: every burst (every so often, and on Space) adds a new effect at a random spot, and each one frees itself once it's finished. Stretched sparks |
| Leaf Explosion | One effect that spawns at user locations: every burst adds a location at a random spot, and each location removes itself once the effect has played out there. Cheaper than an effect per burst, and what Run on GPU effects need. Freely aligned, tumbling leaves |
| Space Weed | Ribbons |
| GPU Simulated | 100,000 particles simulated on the GPU, spawning at one location the scene moves around |
| Burst | Baked sprite data playback, billboards and a ribbon: an animation instance every 250 ms at a random spot at the back, each one interpolated by a compute shader with no simulation at all |

A backend tells the scene what it can draw with `example_scene_features`, and the scene leaves out the effects it
can't draw yet. With no features at all it plays the first three, which only need billboards.

## Setting up a library

The library and the renderer have to be set up in this order (`example_scene_Load`, then the backend's uploads):

1. `tfx_LoadEffectLibrary(path, tfx_impl_ShapeLoader, tfx_impl_GetUV, &shapes)`. The shape loader decodes every shape
   to RGBA8 frames and hands the library a pointer to each one.
2. `tfx_CreateEffectTemplate` for each effect you'll play, so that everything the templates add is in the upload.
3. For baked playback, `tfx_CreateAnimationManager` and `tfx_LoadSpriteData(path, manager, tfx_impl_ShapeLoader,
   &shapes)`. The sprite data carries its own shapes, and loading them into the same shape list puts them in the same
   atlas.
4. `tfx_impl_BuildShapeAtlas` packs every frame into the layers of one texture array.
5. `tfx_UpdateLibraryGPUImageData`. The library calls the uv lookup again now that the atlas exists. The first call,
   during the load, happens before the atlas exists. For baked playback, `tfx_BuildAnimationManagerGPUShapeData` with
   `tfx_impl_GetUV` does the same for the sprite data's own `tfx_gpu_shapes` list.
6. The backend uploads the atlas, the color ramps (`tfx_GetColorRampBitmap`), the shape frame data
   (`tfx_GetLibraryGPUShapes`), the emitter properties (`tfx_GetParticlePropertiesBuffer`) and the graph lookups
   (`tfx_GetGPUGraphLookupsBuffer`, which ribbons and the GPU kernel read). The animation manager has its own color
   ramps, shape frames, emitter properties, baked sprites and ribbons, uploaded once
   (`tfx_impl_dx11_UploadSpriteData`). Only the atlas is shared.

Free the stage before the library, and the shapes after the library, because the library holds pointers to them
(`example_scene_Shutdown` does all three). Call `tfx_EndTimelineFX` after that.

### Why an atlas

The reference shaders give every shape its own texture array and pick the one a particle samples per pixel through a
bindless index. D3D11, OpenGL and WebGPU can't index textures per pixel, so the common layer packs every frame of
every shape into one RGBA8 texture array instead. The uv lookup writes each frame's rect and puts its layer in
`texture_array_index`, leaving the upper 16 bits (the bindless index in the reference shaders) at 0. Each frame gets 2
pixels of padding filled with its own edge pixels, so bilinear filtering stays clean. The atlas makes its own mips, so a
little bleed between neighbours at a distance is expected.

## The frame loop

The stage updates at a fixed rate, 60 times a second in the examples. The vertex shader interpolates every sprite
between the previous update and the latest one, so the frame rate and the update rate are independent of each other.
This is `example_scene_Update` and the backend calls after it, cut down:

```c
//Mutators aren't safe while an update is running, so add effects and set the camera and positions first
tfx_SetStageCamera(stage, camera_front, camera_position);

tfxU32 ticks = tfx_impl_AdvanceTimer(&timer, elapsed_milliseconds);
if (ticks > 0) {
	//One update however many ticks are due, so the previous instance buffer is always the one drawn last
	tfx_UpdateStage(stage, timer.update_length * ticks);
}
//Waits for the update and does nothing when ticks is 0. Only fetch GPU particles if the backend dispatches them
tfx_impl_GatherFrame(stage, ticks, &ribbon_staging, fetch_gpu_particles, &frame);
//Baked playback runs on the real frame time rather than the stage's ticks
tfx_UpdateAnimationManager(animation_manager, (float)elapsed_milliseconds);
tfx_impl_GatherBakedFrame(animation_manager, &frame.baked);
tfx_impl_SetUniforms(&uniforms, &camera, screen_width, screen_height, &timer);

tfx_impl_dx11_UploadFrame(&renderer, &frame);	//Baked data every frame, instances and ribbons only when ticks > 0
tfx_impl_dx11_Compute(&renderer, &uniforms);	//GPU particles on update frames, baked playback and ribbon geometry every frame
tfx_impl_dx11_Render(&renderer, &uniforms);		//Every frame, interpolating with uniforms.timer_lerp
```

The backend keeps two instance buffers and swaps them on every update, so the other buffer always holds what was drawn
last, which is what `captured_index` points into. That's true for CPU and GPU instances alike, as long as there's one
stage update per frame. `tfx_UpdateStageSubSteps` would put CPU instances on a previous layout of their own, which can
overlap the GPU instances the kernel interpolates from. On the frames between updates nothing is uploaded or
dispatched for the stage, and the GPU keeps drawing the same two buffers with a growing lerp.

Rules the loop depends on:
- **Threading.** `tfx_UpdateStage` runs on its own thread. Add effects, move them and set the camera before it, and
  after it only read through the getters that wait for the update (the ones `tfx_impl_GatherFrame` uses), or call
  `tfx_CompleteStageWork` first.
- **Every GPU particle fetch is dispatched.** Each fetch carries on from the particle state the last dispatch left in
  the pool, so never call `tfx_GetGPUParticleFrame` on a frame you won't dispatch.
- **Effect ids are generational.** Keep the id `tfx_AddEffectTemplateToStage` returns and check it with
  `tfx_EffectIDIsValid`. A fire and forget effect's id goes stale once it frees itself.

## Direct3D 11 notes

- **No vertex buffers or input layout.** `tfx_instance_t` stores its size as signed 16-bit integers that the shader
  scales, which Vulkan reads with an SSCALED format. D3D has no scaled formats, so the vertex shader reads the current
  and previous instance from raw buffers by `SV_InstanceID` and unpacks them itself. The quad corners come from
  `SV_VertexID`, with 6 vertices per instance.
- **`SV_InstanceID` doesn't include `StartInstanceLocation`**, so a draw's first instance goes in a constant buffer.
- **The push constants become constant buffers.** The frame constants (`tfx_impl_uniforms_t`) are at b0 and the
  per-draw instance offsets at b1. Every other resource has a fixed register, so the bindless indexes are gone. The
  shaders are compiled with `-D TFX_FIXED_BINDINGS`, and the `[Fixed_render_bindings]` section of `timelinefx.slang`
  lists every register.
- **Struct layouts are checked.** `timelinefx.slang` mirrors the library's GPU structs by hand, so `tfx_impl_dx11.cpp`
  static_asserts their sizes and names the struct to fix when one changes.
- **Instance buffers** are `DEFAULT` raw buffers, with UAVs when the GPU kernel writes them. The CPU instances go up in
  one `UpdateSubresource`, because the library lays the buffer out layer by layer. Each layer is then one
  `DrawInstanced(6, count)`. When the buffers grow, both copy their old contents across, since the previous buffer is
  still what `captured_index` points into.
- **Blending:** the pixel shader writes premultiplied alpha, so it blends with ONE / INV_SRC_ALPHA, with the depth test
  on and depth writes off.
- **Views are unbound after each pass.** D3D11 won't bind a buffer for writing while it's bound for reading, so every
  compute pass and draw unbinds what the next one writes. That keeps the debug layer quiet, which is the test every
  part of the example has passed.
- **Ribbons** are uploaded on update frames and rebuilt by `tfx_ribbon_compute` every frame, because relative ribbons follow
  their emitter with the frame's lerp. `tfx_impl_dx11_Compute` runs one dispatch per bucket with the bucket's
  `tfx_ribbon_bucket_globals_t` in a constant buffer. `tfx_impl_dx11_Render` then walks the same buckets again for the
  indexed draws. The compute shader writes the vertices to a structured buffer that the vertex shader pulls from, and
  writes the indices to a raw buffer that's also bound as the index buffer. It uses 256 threads per group rather than
  the GLSL's 1024, which fxc warns is over its register budget.
- **GPU particles:** the kernel is compiled with `-D TFX_FIXED_BINDINGS`, which binds the pool at u0, the current
  instance buffer at u1, the graph lookups at t0, the frame's upload at t1 and `tfx_gpu_particle_push_t` at b0. The
  buffer indexes in the push constants are unused, and the `*_base` rows say where each section starts in the pool
  and the upload.
  - The pool is created once from `tfx_GetGPUParticleBufferSizes`, zeroed.
  - The upload is a dynamic buffer filled by `tfx_CopyGPUParticleUpload` straight into mapped memory.
  - Each `tfx_gpu_particle_dispatch_t` is one dispatch of `(slot_count + 255) / 256` groups, with `dispatch_index`
    changed in the constant buffer each time.
  - The kernel writes after the CPU instances, so the instance buffers are sized for `gpu_instance_start +
    gpu_instance_count` and each layer is drawn as two ranges, its CPU instances then its GPU instances.
  - Every Run on GPU emitter gets a ring sized for one spawn location, so an effect bursting at several locations at
    once needs `tfx_SetTemplateGPUCapacityForLocations`, and the stage's `gpu_particle_capacity` has to hold every ring.
    That function returns the pool slots the template's rings take, rounded the way the pool hands them out, so adding
    up its returns gives a pool that fits exactly. `example_scene.cpp` does both. The capacities assume about 60
    updates a second; at 20 the rings drop spawns.
- **Baked playback:**
  - Each frame the offsets and animation instances go up.
  - `tfx_sprite_playback` interpolates every visible baked sprite into a billboard instance.
  - `tfx_billboard_baked_vertex` draws it. That's a vertex shader with no interpolation of its own, because the compute
    has already interpolated between baked frames, resolved the image frame and baked the stretch.
  - `tfx_ribbon_playback` interpolates the baked ribbons into the same `tfx_ribbon_t` rows a stage uploads, and they go
    through `tfx_ribbon_compute` and the ribbon draw like live ones. The buckets come from
    `tfx_NextAnimationRibbonDispatch` instead of `tfx_NextRibbonDispatch`, and the graph lookups from the animation
    manager (`tfx_GetAnimationRibbonGraphLookupsBufferPointer`), because a `.tfxsd` plays without a library.
- **Shaders** are compiled ahead of time by slangc, one `<entry point>.cso` each, and `tfx_impl_dx11_Init` loads them
  from the folder the CMake project builds them into.

## Clip space

`tfx_impl_Perspective` builds a right handed projection. Every backend passes the flags for its own clip space:

| API | Flags |
| --- | --- |
| D3D11, D3D12, Metal, WebGPU | `tfx_impl_clip_space_depth_zero_to_one` |
| Vulkan | `tfx_impl_clip_space_depth_zero_to_one \| tfx_impl_clip_space_flip_y` |
| OpenGL | none |

## Porting to a new API

Use the D3D11 backend as the template: keep its function names with your API's prefix, and its order of work.

**Setting up**
- [ ] Add `example_<platform>_<api>/main.cpp`: a window, a device, a swap chain and the loop. It calls
      `example_scene_Load`, your upload functions, `example_scene_Start`, then `example_scene_Update` and your upload,
      compute and render calls each frame, with the clip space flags from the table above.
- [ ] Add `tfx_impl_<api>.h/.cpp` with the same functions as `tfx_impl_dx11.h`: `Init`, `UploadLibrary`,
      `UploadSpriteData`, `UploadFrame`, `Compute`, `Render` and `Shutdown`. The backend draws into whatever target
      `main.cpp` has bound, so it never owns the swap chain.
- [ ] Add a target to `CMakeLists.txt`, with a slangc build of each entry point your backend uses.
- [ ] Turn on the API's validation layer in debug builds and print what it says. Every step below is done when it's
      silent.
- [ ] Static_assert the sizes of every library struct your shaders mirror, as `tfx_impl_dx11.cpp` does.

**Shaders.** Never ported by hand: compile `../Shaders/timelinefx.slang` for your target, as `../Shaders/readme.md`
describes. Vulkan and D3D12 use its default bindless model. D3D11 compiles it with `-D TFX_FIXED_BINDINGS`, which reads
every resource from a fixed register and every shape from one atlas, and reads the instances and ribbon vertexes from
buffers in place of vertex attributes. OpenGL and WebGPU are meant to use the same define, but neither has been tried.

**1. CPU billboards** (`example_scene_feature_none`)
- [ ] The atlas as an RGBA8 texture array with mips, and a clamped linear sampler with mip filtering.
- [ ] The color ramps as a texture array, one layer per `tfx_GetColorRampBitmap`, in the format of its bytes per
      pixel (RGBA16F unless the game called `tfx_SetColorRampFormat`). It's read with an integer load, never sampled.
- [ ] The shape frames and emitter properties as read only buffers.
- [ ] Two instance buffers that swap on every update, sized for the instance count and copied across when they grow.
- [ ] One instanced draw of 6 vertices per layer. If the API's instance id doesn't include the draw's first instance
      (D3D11, D3D12 and OpenGL's `gl_InstanceID`), pass the first instance in yourself. Vulkan's `gl_InstanceIndex`
      and WGSL's `instance_index` include it.
- [ ] Read the size as signed 16-bit values, through an SSCALED vertex format where the API has one and by unpacking
      in the shader where it doesn't.
- [ ] Premultiplied blending (ONE / ONE_MINUS_SRC_ALPHA), depth test on, depth write off.
- [ ] Check: the torch circles with its smoke trailing behind, the sparks stretch along their motion, the leaves tumble,
      and running at 20 and 400 updates a second looks the same as 60.

**2. Ribbons** (`example_scene_feature_ribbons`)
- [ ] Read only segment, ribbon and emitter buffers, uploaded on update frames.
- [ ] A vertex buffer and an index buffer the ribbon compute writes, sized from `tfx_ribbon_batch_sizes_t`, then read
      by the draw. Put a barrier between them on APIs that need one.
- [ ] Every frame, walk the buckets twice: once to dispatch the compute, once for the indexed draws. Set the bucket
      globals' `lerp` to the frame's lerp before each dispatch.
- [ ] A thread group size that's a multiple of 32.
- [ ] Check: Space Weed's ribbons wave smoothly, taper and fade.

**3. GPU particles** (`example_scene_feature_gpu_simulation`)
- [ ] The kernel compiled for your binding model, with slangc as part of the build.
- [ ] The pool created zeroed, once, from `tfx_GetGPUParticleBufferSizes`. Nothing else ever clears it.
- [ ] The frame's upload filled by `tfx_CopyGPUParticleUpload`, and the push constants' bases set from the pool sizes
      and the frame.
- [ ] Instance buffers that the kernel can write, sized for `gpu_instance_start + gpu_instance_count`.
- [ ] One dispatch per record, every fetch, before the draw that reads the current instance buffer.
- [ ] Each layer drawn as its CPU range then its GPU range.
- [ ] Check: `tfx_GetGPUParticleStats` reports no dropped spawns at 60 updates a second, and the leaves and the GPU
      Simulated fountain both draw.
- [ ] If they don't, print `tfx_DescribeStage` and read the `TimelineFX tip:` lines in the log. Then turn the GPU
      simulation off with `tfx_SetStageGPUSimulation(stage, false)`, clear the stage and add the effects again. Still
      missing means the effects aren't spawning; appearing means the bug is in your upload, dispatch or draw ranges.

**4. Baked playback** (`example_scene_feature_baked`)
- [ ] The animation manager's own color ramps, shape frames, emitter properties and sprites uploaded once.
- [ ] The offsets and animation instances uploaded every frame, then the sprite playback compute into a buffer of
      billboard instances, drawn with the static billboard vertex shader and the normal pixel shader.
- [ ] The ribbon playback compute into a second set of ribbon buffers, then the same ribbon compute and draw as live
      ribbons, with the animation manager's graph lookups and buckets.
- [ ] Check: the Bursts at the back of the scene play their billboards and their electricity ribbon.

**Done when** the new example looks the same as the D3D11 one side by side, with a silent validation layer.

## Engines with their own instance buffers

The checklist assumes the backend owns its two instance buffers. An engine that already double buffers its instances,
or records its passes in a cached render graph, has to keep these rules itself:

- **One buffer holds both.** The GPU instances go straight after the CPU ones. Upload only `[0, gpu_instance_start)`
  from the CPU, and make room for `gpu_instance_start + gpu_instance_count` instances before the compute writes.
- **Swap only on update frames.** `captured_index` points into the buffer drawn last, for CPU and GPU instances alike,
  so the buffers swap once per stage update and never on the frames in between.
- **Carry the GPU instances if the buffer turns over without a dispatch.** If the engine resets or swaps the instance
  buffer on a frame with no GPU dispatch, for example to redraw while paused, copy
  `[gpu_instance_start, + gpu_instance_count)` from the buffer the kernel last wrote into the new one. Otherwise the
  GPU particles vanish for that frame.
- **Dispatch each fetch exactly once.** A host that replays dispatches from a cached graph must set its pending
  dispatch count to zero once it has dispatched. Otherwise a frame with no update runs the same ticks again, and the
  particles jump ahead.
- **Idle frames are idle.** On frames without a stage update the upload and simulate passes have nothing to do. If the
  engine has a GPU profiler, mark them idle rather than timing them as empty, or the averages come out too low.
