# TimelineFX implementation examples

Minimal examples of drawing TimelineFX effects with different graphics APIs, laid out the way Dear ImGui lays out its
backends. Each API has one backend `.h`/`.cpp` pair, and every backend shares one renderer agnostic pair:

| File | What it holds |
| --- | --- |
| `tfx_impl_common.h/.cpp` | Everything that would be written the same way whatever API draws the particles: shape loading, the shape atlas, the fixed step timer, the camera matrices and gathering each frame's instances from a stage. Nothing in it talks to a GPU. |
| `tfx_impl_dx11.h/.cpp` | The Direct3D 11 backend. |
| `shaders/tfx_dx11.hlsl` | Every D3D11 shader, ported from the GLSL reference shaders in `../Shaders`. |
| `example_scene.h/.cpp` | The demo scene every example plays. It's example code showing how a game uses the library, separate from the integration code above. |
| `example.tfx` | The library the scene plays. |
| `example.tfxsd` | Burst, recorded from `example.tfx` in the editor, for baked playback. |
| `example_win32_dx11/main.cpp` | A Win32 window, a D3D11 device and swap chain, and the main loop that calls the scene and the backend. |
| `third_party/` | stb_image (PNG shapes), stb_rect_pack (the atlas), bcdec (BC4/5/7 KTX2 shapes) and zstd's single file decompressor (supercompressed KTX2 shapes). |

The code is C-style C++: plain structs and free functions, no STL. The headers are wrapped in `extern "C"`, so a C
project can call a backend that was built as C++. The examples only include `timelinefx.h`, never
`timelinefx_internal.h`, which is exactly what a game sees.

Status: D3D11 draws CPU simulated billboard particles, ribbons, GPU simulated particles and baked billboard playback.
Baked ribbon playback is waiting on library API for it.

## Building and running

The examples build as part of the TFXEditor CMake project (`TFX_BUILD_EXAMPLES`, on by default), or on their own:

```
cmake -S TimelineFXLib/implementations -B build-examples
cmake --build build-examples --config Debug
```

The GPU particle kernel is compiled from `../Shaders/tfx_gpu_particles.slang` with slangc, which comes with the Vulkan
SDK. Without slangc the examples still build, just without GPU simulation.

```
example_win32_dx11 [updates per second]
```

Space fires a burst of sparks, O pauses the camera orbit, the mouse wheel zooms and Escape quits. In debug builds the
D3D11 debug layer's messages are repeated on the console.

## The scene

| Effect | What it shows |
| --- | --- |
| Spitting Torch | An ordinary effect moved with `tfx_SetEffectPosition` every frame. Animated flames, smoke and embers trail behind it |
| Burst Sparks | Fire and forget: every burst (every so often, and on Space) adds a new effect at a random spot, and each one frees itself once it's finished. Stretched sparks |
| Leaf Explosion | One effect that spawns at user locations: every burst adds a location at a random spot, and each location removes itself once the effect has played out there. Cheaper than an effect per burst, and what Run on GPU effects need. Freely aligned, tumbling leaves |
| Space Weed | Ribbons |
| GPU Simulated | 100,000 particles simulated on the GPU, spawning at one location the scene moves around |
| Burst | Baked sprite data playback: an animation instance every 250 ms at a random spot at the back, each one interpolated by a compute shader with no simulation at all. Its ribbon isn't drawn yet |

A backend tells the scene what it can draw with `example_scene_features`, and the scene leaves out the effects it
can't draw yet.

## Setting up a library

The library and the renderer have to be set up in this order:

1. `tfx_LoadEffectLibrary(path, tfx_impl_ShapeLoader, tfx_impl_GetUV, &shapes)`. The shape loader decodes every shape
   to RGBA8 frames and hands the library a pointer to each one.
2. `tfx_CreateEffectTemplate` for each effect you'll play.
3. `tfx_impl_BuildShapeAtlas` packs every frame into the layers of one texture array.
4. `tfx_UpdateLibraryGPUImageData`. The library calls the uv lookup again now that the atlas exists. The first call,
   during the load, happens before the atlas exists.
5. The backend uploads the atlas, the color ramps (`tfx_GetColorRampBitmap`), the shape frame data
   (`tfx_GetLibraryGPUShapes`) and the emitter properties (`tfx_GetParticlePropertiesBuffer`).

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
tfx_impl_dx11_UploadFrame(&renderer, &frame);
tfx_impl_dx11_Compute(&renderer, &uniforms);	//GPU particles on update frames, ribbon geometry every frame
//Draw every frame, passing timer.lerp to the vertex shader as the interpolation value
tfx_impl_dx11_Render(&renderer, &uniforms);
```

The backend keeps two instance buffers and swaps them on every update, so the other buffer always holds what was drawn
last, which is what `captured_index` points into. That's true for CPU and GPU instances alike, as long as there's one
stage update per frame. `tfx_UpdateStageSubSteps` would put CPU instances on a previous layout of their own, which can
overlap the GPU instances the kernel interpolates from. On the frames between updates nothing is uploaded or
dispatched, and the GPU keeps drawing the same two buffers with a growing lerp.

## Direct3D 11 notes

- **No vertex buffers or input layout.** `tfx_instance_t` stores its size as signed 16-bit integers that the shader
  scales, which Vulkan reads with an SSCALED format. D3D has no scaled formats, so the vertex shader reads the current
  and previous instance from raw buffers by `SV_InstanceID` and unpacks them itself. The quad corners come from
  `SV_VertexID`, with 6 vertices per instance.
- **`SV_InstanceID` doesn't include `StartInstanceLocation`**, so a draw's first instance goes in a constant buffer.
- **The push constants become constant buffers.** The frame constants (`tfx_impl_uniforms_t`) are at b0 and the
  per-draw instance offsets at b1. Every other resource has a fixed register, so the bindless indexes are gone.
- **Blending:** the pixel shader writes premultiplied alpha, so it blends with ONE / INV_SRC_ALPHA, with the depth test
  on and depth writes off.
- **Ribbons** are uploaded on update frames and rebuilt by `ribbon_compute` every frame, because relative ribbons follow
  their emitter with the frame's lerp. `tfx_impl_dx11_Compute` runs one dispatch per bucket with the bucket's
  `tfx_ribbon_bucket_globals_t` in a constant buffer. `tfx_impl_dx11_Render` then walks the same buckets again for the
  indexed draws. The compute shader writes the vertices to a structured buffer that the vertex shader pulls from, and
  writes the indices to a raw buffer that's also bound as the index buffer.
- **GPU particles:** the kernel is compiled with `-D TFX_FIXED_BINDINGS`, which binds the pool at u0, the current
  instance buffer at u1, the graph lookups at t0, the frame's upload at t1 and `tfx_gpu_particle_push_t` at b0. The
  buffer indexes in the push constants are unused, and the `*_base` rows say where each section starts in the pool
  and the upload.
  - The pool is created once from `tfx_GetGPUParticleBufferSizes`, zeroed.
  - The upload is a dynamic buffer filled by `tfx_CopyGPUParticleUpload` straight into mapped memory.
  - Each `tfx_gpu_particle_dispatch_t` is one dispatch of `(slot_count + 255) / 256` groups, with `dispatch_index`
    changed in the constant buffer each time.
  - Every fetch must be dispatched, which is why `tfx_impl_GatherFrame` only fetches when asked.
  - The kernel writes after the CPU instances, so each layer is drawn as two ranges, its CPU instances then its GPU
    instances.
  - Every Run on GPU emitter gets a ring sized for one spawn location, so an effect bursting at several locations at
    once needs `tfx_SetTemplateEmitterGPUCapacity`, and the stage's `gpu_particle_capacity` has to hold every ring.
    `example_scene.cpp` works both out.
- **Baked playback:**
  - `tfx_impl_dx11_UploadSpriteData` uploads the animation manager's own color ramps, shape frames, emitter
    properties and baked sprites. Only the atlas is shared with the library, because the sprite data's shapes go
    through the same shape loader before the atlas is built.
  - Each frame the offsets and animation instances go up.
  - `sprite_playback` interpolates every visible baked sprite into a billboard instance.
  - `billboard_static_vertex` draws it. That's a vertex shader with no interpolation of its own, because the compute
    has already interpolated between baked frames, resolved the image frame and baked the stretch.
- **Shaders** compile at startup with `D3DCompileFromFile`, which comes with Windows. The GPU kernel is compiled ahead
  of time by slangc and loaded.

## Clip space

`tfx_impl_Perspective` builds a right handed projection. Every backend passes the flags for its own clip space:

| API | Flags |
| --- | --- |
| D3D11, D3D12, Metal, WebGPU | `tfx_impl_clip_space_depth_zero_to_one` |
| Vulkan | `tfx_impl_clip_space_depth_zero_to_one \| tfx_impl_clip_space_flip_y` |
| OpenGL | none |
