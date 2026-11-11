# TimelineFX implementation examples

Minimal examples of drawing TimelineFX effects with different graphics API backends. Each API has one backend `.h`/`.cpp` pair, and every backend shares one renderer agnostic pair:

| File | What it holds |
| --- | --- |
| `tfx_impl_common.h/.cpp` | Everything that would be written the same way whatever API draws the particles: shape loading, the shape atlas, the fixed step timer, the camera matrices and gathering each frame's instances from a stage. Nothing in it talks to a GPU. |
| `tfx_impl_dx11.h/.cpp` | The Direct3D 11 backend. |
| `tfx_impl_dx12.h/.cpp` | The Direct3D 12 backend, with the bindless shaders: one descriptor heap and root constants. |
| `tfx_impl_vulkan.h/.cpp` | The Vulkan 1.3 backend, with the bindless shaders: one descriptor set, push constants and dynamic rendering. |
| `tfx_impl_sdlgpu.h/.cpp` | The SDL3 GPU (SDL_gpu) backend, for its Vulkan and Direct3D 12 drivers. |
| `tfx_impl_webgpu.h/.cpp` | The WebGPU backend, written against the standard `webgpu.h` and run on wgpu-native. Uploads go through the queue, uniforms through dynamic offsets. |
| `tfx_impl_gl.h/.cpp` | The OpenGL 4.3 core backend. It includes no GL header and loads its functions through a loader you pass it, such as `glfwGetProcAddress`. |
| `tfx_impl_metal.h/.mm` | The Metal backend, Objective-C++ without ARC. Fixed bindings for the render passes, the kernel's `gpuAddress` table, and frames in flight with upload rings. |
| `../Shaders/timelinefx.slang` | Every shader: the GPU particle kernel and the render passes. CMake has slangc compile each entry point an example uses, to `.cso` for D3D11, bindless `.dxil` for D3D12, bindless `.spv` for Vulkan, `.wgsl` for WebGPU, `.spv` and `.dxil` for SDL_gpu, `.spv` then `.glsl` through spirv-cross for OpenGL, and `.metal` source for Metal. |
| `example_scene.h/.cpp` | The demo scene every example plays. It's example code showing how a game uses the library, separate from the integration code above. |
| `example.tfx` | The library the scene plays. |
| `example_second.tfx` | A second library, so the scene draws effects from two libraries with one set of buffers. The scene plays without it. |
| `example.tfxsd` | Burst, recorded from `example.tfx` in the editor, for baked playback. |
| `example_win32_dx11/main.cpp` | A Win32 window, a D3D11 device and swap chain, and the main loop that calls the scene and the backend. |
| `example_win32_dx12/main.cpp` | A Win32 window, a D3D12 device, queue, swap chain and frames in flight, and the same main loop. |
| `example_sdl3_gpu/main.cpp` | An SDL3 window and GPU device, Vulkan by default or Direct3D 12 with `--d3d12`, and the same main loop. |
| `example_glfw_opengl/main.cpp` | A GLFW window with a 4.3 core debug context, and the same main loop. |
| `example_glfw_vulkan/main.cpp` | A GLFW window, a Vulkan 1.3 instance, device and swap chain with frames in flight, and the same main loop. |
| `example_glfw_webgpu/main.cpp` | A GLFW window, a wgpu-native instance, adapter, device and surface, Vulkan or D3D12 underneath, and the same main loop. |
| `example_glfw_metal/main.mm` | A GLFW window with a `CAMetalLayer`, a Metal device and queue, frames in flight, and the same main loop. |
| `third_party/` | stb_image (PNG shapes), stb_rect_pack (the atlas), bcdec (BC4/5/7 KTX2 shapes) and zstd's single file decompressor (supercompressed KTX2 shapes). |

So there are three layers. `tfx_impl_common` is integration code you can copy into a game as it is, `example_scene` is
what this demo does with the library, and `main.cpp` is only the window, the device and the calls into the other two.
A new graphics API only needs its own `tfx_impl_<api>` pair and a `main.cpp`. The test for where code goes is whether
every backend would write it the same way: if it would, it belongs in common.

The code is C-style C++: plain structs and free functions, no STL. The headers are wrapped in `extern "C"`, so a C
project can call a backend that was built as C++. The examples only include `timelinefx.h`, never
`timelinefx_internal.h`, which is exactly what a game sees. `tfx_impl_common` and the backends report problems with
`tfx_Print`, so they reach the same `tfx_SetPrintCallback` log as the library's own messages. Do the same in yours.

Status: D3D11, D3D12 (bindless), Vulkan (bindless), WebGPU (wgpu-native), SDL_gpu (Vulkan and D3D12), OpenGL 4.3 and Metal draw CPU simulated billboards and ribbons, GPU simulated particles, and
baked billboard and ribbon playback, and layers that write depth and lit layers. Every other API is still to come, see [Porting to a new API](#porting-to-a-new-api).

## Building and running

You can use cmake to build the examples:

```
cmake -S TimelineFXLib/implementations -B build-examples
cmake --build build-examples --config Debug
```

Every shader is compiled from `../Shaders/timelinefx.slang` with slangc, which comes with the Vulkan SDK and with
the Slang releases on GitHub. CMake looks on the `PATH` and in the Vulkan SDK, and stops if it can't find it, so put
its bin folder on the `PATH` or set `TFX_SLANGC` to it.

The SDL3 example is off by default, because it downloads and builds SDL3 with FetchContent. Turn it on with
`-D TFX_EXAMPLE_SDL3=ON` when you configure. The OpenGL example is off by default too, because it needs GLFW, which it
downloads unless the parent project already has a `glfw` target, and spirv-cross, which comes with the Vulkan SDK. Turn
it on with `-D TFX_EXAMPLE_GL=ON`. The D3D12 example is off by default to keep the build short: turn it on with
`-D TFX_EXAMPLE_DX12=ON`. The Vulkan example needs the Vulkan SDK (CMake's `find_package(Vulkan)`) and GLFW, fetched
the same way: turn it on with `-D TFX_EXAMPLE_VULKAN=ON`. The WebGPU example downloads wgpu-native's release for Windows x64
or Linux x64 (or set `TFX_WGPU_NATIVE_DIR` to one you've unzipped) and uses GLFW: turn it on with `-D TFX_EXAMPLE_WEBGPU=ON`.
The Metal example builds on Apple platforms only and uses GLFW. It needs no Metal toolchain, because it compiles its shaders at
start up: turn it on with `-D TFX_EXAMPLE_METAL=ON`.

```
example_win32_dx11 [updates per second]
example_win32_dx12 [updates per second] [--frames N] [--test-debug-output] [--gpu-validation]
example_sdl3_gpu [--d3d12] [updates per second]
example_glfw_opengl [updates per second] [--frames N] [--test-debug-output]
example_glfw_vulkan [updates per second] [--frames N] [--test-debug-output] [--gpu-validation]
example_glfw_webgpu [updates per second] [--frames N] [--test-debug-output] [--gpu-validation] [--vulkan | --d3d12]
example_glfw_metal [updates per second] [--frames N] [--test-debug-output] [--gpu-validation]
```

Space fires a burst of sparks, O pauses the camera orbit, L frees and reloads the second library, the mouse wheel
zooms and Escape quits. The window title shows the frame rate and the particle counts. In debug builds the D3D11
debug layer's messages are repeated on the console, and so are the Vulkan validation layer's and the D3D12 debug
layer's in the SDL3 example. The OpenGL example asks for a debug context and the backend prints the driver's messages,
the D3D12 and Vulkan examples print their debug and validation layers' messages, the WebGPU example prints WebGPU's
validation errors and wgpu-native's log, and the Metal example turns on Metal's API validation and prints what shader
validation finds. In the examples that take them,
`--test-debug-output` makes one deliberate error to prove the messages arrive, `--gpu-validation` turns on the GPU
side of the layer, which checks what the shaders do with each descriptor, and `--frames N` quits after N frames.
Passing an update rate (20 or 400, say) is a quick way to check the interpolation: the effects should look the same at
any rate. The examples size the GPU rings for it with `tfx_SetGPUSizingUpdateRate` before the scene starts.

## The scene

| Effect | What it shows |
| --- | --- |
| Spitting Torch | Animated flames and smoke on the CPU, and embers on the GPU, trailing behind it. Because it has a Run on GPU emitter it spawns at a user location, which the scene moves around a circle with `tfx_UpdateSpawnLocation` every frame |
| Burst Sparks | Fire and forget: every burst (every so often, and on Space) adds a new effect at a random spot, and each one frees itself once it's finished. Stretched sparks |
| Leaf Explosion | One effect that spawns at user locations: every burst adds a location at a random spot, and each location removes itself once the effect has played out there. Cheaper than an effect per burst, and what Run on GPU effects need. Freely aligned, tumbling leaves, on layer 2, which writes depth and is lit |
| Space Weed | Ribbons, on layer 2 with the leaves |
| GPU Simulated | 100,000 particles simulated on the GPU, spawning at one location the scene moves around |
| Burst | Baked sprite data playback, billboards and a ribbon: an animation instance every 250 ms at a random spot at the back, each one interpolated by a compute shader with no simulation at all |
| Fireflies | From `example_second.tfx`. It shares the stage, the draw ranges and the buffers with the effects from `example.tfx`, see [Setting up a library](#setting-up-a-library) |

A backend tells the scene what it can draw with `example_scene_features`, and the scene leaves out the effects it
can't draw yet. With no features at all it plays the first three, which only need billboards.

## Setting up a library

The libraries and the renderer have to be set up in this order (`example_scene_Load`, then the backend's uploads):

1. `tfx_LoadEffectLibrary(path, tfx_impl_ShapeLoader, tfx_impl_GetUV, &shapes)` for each library. The shape loader
   decodes every shape to RGBA8 frames and hands the library a pointer to each one. Give each library its own
   `tfx_impl_shapes_t`, so its shapes get their own layers of the atlas and it can be freed and reloaded on its own.
2. `tfx_CreateEffectTemplate` for each effect you'll play. Templates can be made at any time: they share the library's
   GPU data until a template customises a graph the GPU reads.
3. For baked playback, `tfx_CreateAnimationManager` and `tfx_LoadSpriteData(path, manager, tfx_impl_ShapeLoader,
   &sprite_shapes)`. The sprite data carries its own shapes, which get a shape list of their own too.
4. `tfx_impl_BuildShapeAtlas(&shapes, &atlas, first_layer, layer_size)` for each shape list packs its frames into
   layers from `first_layer` on, and the layers the uv lookup writes include it. Every list goes into one texture
   array, so they share a layer size: the largest `tfx_impl_GetAtlasLayerSize` of them. The scene puts the lists one
   after another and keeps spare layers after the second library, so a reload can grow.
5. `tfx_UpdateLibraryGPUImageData` for each library. The library calls the uv lookup again now that the atlas exists.
   The first call, during the load, happens before the atlas exists. For baked playback,
   `tfx_BuildAnimationManagerGPUShapeData` with `tfx_impl_GetUV` does the same for the sprite data's own
   `tfx_gpu_shapes` list.
6. `tfx_UpdateGlobalGPUData` lays every library out in one set of global buffers, and the stage is created with
   `global_library_indexes` on, so the indexes it writes point into them. It returns false if the libraries need more
   room than the global buffers have. Effects can only be added from libraries it has placed.
7. The backend creates one texture array with every list's layers (`tfx_impl_dx11_CreateAtlas`) and uploads each atlas
   into its own layers (`tfx_impl_dx11_UploadAtlas`), making mips for those layers only. Everything else goes up on update frames, whenever
   its version has changed since the last upload: the emitter properties (`tfx_GetGlobalParticlePropertiesBuffer`), the
   shape frame data (`tfx_GetGlobalGPUShapes`), the color ramps one layer at a time (`tfx_GetGlobalColorRampLayer`) and
   the graph lookups (`tfx_GetGPUGraphLookupsBuffer`, which ribbons and the GPU kernel read). So loading a library or
   customising a template later needs no extra code. The animation manager has its own color ramps, shape frames,
   emitter properties, baked sprites and ribbons, uploaded once (`tfx_impl_dx11_UploadSpriteData`). Only the atlas is
   shared.

The section "Drawing effects from more than one library" in `timelinefx.h` has the whole contract.

**One library only.** A game that only ever has one library loaded can leave `global_library_indexes` off, skip
`tfx_UpdateGlobalGPUData` and upload that library's own data instead: `tfx_GetColorRampBitmap`,
`tfx_GetLibraryGPUShapes` and `tfx_GetParticlePropertiesBuffer`, again when `tfx_GetColorRampsVersion` or
`tfx_GetParticlePropertiesVersion` change. With more libraries that way, each draw range can only hold effects from one
library, so use a sort key per library (`tfx_SetEffectSortKey`) and bind that library's buffers for its ranges.

Free the stage before the library, and the shapes after the library, because the library holds pointers to them
(`example_scene_Shutdown` does all three). Call `tfx_EndTimelineFX` after that.

**Freeing and reloading a library.** Press L in an example and the scene frees the second library and loads it again
while the first keeps playing (`example_scene__reload_second_library`). It waits for an update frame and does it before
the stage update, so the instances drawn next have the new layout. To free on a frame without an update (a paused game),
call `tfx_UpdateStage(stage, 0)` after `tfx_UpdateGlobalGPUData` and pass `stage_updated` true to `tfx_impl_GatherFrame`,
or the backend keeps drawing instances written with the old layout:

1. `tfx_FreeLibrary`, which removes its effects from every stage and frees its templates, then `tfx_impl_FreeShapes`.
2. Load it into the same shape list, then create its templates.
3. `tfx_impl_BuildShapeAtlas` at the same `first_layer`. If it needs more layers than were kept for it, the texture
   array has to be made again with every library's atlas rebuilt, and their pixels are gone by then, so keep spares.
4. `tfx_UpdateLibraryGPUImageData` and `tfx_UpdateGlobalGPUData`, then add its effects back to the stage.
5. The backend uploads the new atlas into its layers before the frame is drawn (`example_scene_GetPendingAtlas`). The
   global buffers go up with the frame as usual, when their versions changed. Reloading the same file into the same
   layers leaves them unchanged.

Upload between frames: an API that doesn't track hazards for you needs the GPU finished with the layers first.

### Why an atlas

The reference shaders give every shape its own texture array and pick the one a particle samples per pixel through a
bindless index. D3D11, OpenGL and WebGPU can't index textures per pixel, so the common layer packs every frame of
every shape into one RGBA8 texture array instead. The uv lookup writes each frame's rect and puts its layer in
`texture_array_index`, leaving the upper 16 bits (the bindless index in the reference shaders) at 0. Each frame gets 2
pixels of padding filled with its own edge pixels, so bilinear filtering stays clean. The atlas makes its own mips, so a
little bleed between neighbours at a distance is expected. They're made per layer after each upload, with D3D11's
`GenerateMips` on a view of just that atlas's layers. Their quality decides how particles a pixel or two across look,
so if small particles shimmer or vanish, build the mips on the CPU with a better filter than the API's box filter.

## The frame loop

The stage updates at a fixed rate, 60 times a second in the examples. The vertex shader interpolates every sprite
between the previous update and the latest one, so the frame rate and the update rate are independent of each other.
This is `example_scene_Update` and the backend calls after it, cut down:

```c
//Mutators aren't safe while an update is running, so add effects and set the camera and positions first
tfx_SetStageCamera(stage, camera_front, camera_position);

tfxU32 ticks = tfx_impl_AdvanceTimer(&timer, elapsed_milliseconds);
if (ticks > 0) {
	//Lays out any library or template change in the global buffers, free when nothing changed
	tfx_UpdateGlobalGPUData();
	//Every due tick at its own length, writing the instances once, so the previous instance buffer is always the one drawn last
	tfx_UpdateStageSubSteps(stage, timer.update_length, ticks);
}
//Waits for the update and does nothing without one. Pass true after tfx_UpdateStage(stage, 0) too, see "Freeing and reloading a library"
//Only fetch GPU particles if the backend dispatches them
//The layer flags come from library, see "Layers that write depth and lit layers" below
tfx_impl_GatherFrame(stage, library, ticks > 0, &ribbon_staging, fetch_gpu_particles, &frame);
//Baked playback runs on the real frame time rather than the stage's ticks
tfx_UpdateAnimationManager(animation_manager, (float)elapsed_milliseconds);
tfx_impl_GatherBakedFrame(animation_manager, &frame.baked);
//timer_lerp comes from tfx_GetStageInterpolation, every frame
tfx_impl_SetUniforms(&uniforms, &camera, screen_width, screen_height, &timer, stage);

tfx_impl_dx11_UploadFrame(&renderer, &frame);	//Baked data every frame, instances, ribbons and changed library data only after an update
tfx_impl_dx11_Compute(&renderer, &uniforms);	//GPU particles on update frames, baked playback and ribbon geometry every frame
tfx_impl_dx11_Render(&renderer, &uniforms);		//Every frame, interpolating with uniforms.timer_lerp
```

The backend keeps two instance buffers and swaps them on every update, so the other buffer always holds what was drawn
last, which is what `captured_index` points into. That's true for CPU and GPU instances alike, as long as each frame
makes one `tfx_UpdateStage` or `tfx_UpdateStageSubSteps` call, which both write instances once. On the frames between
updates nothing is uploaded or dispatched for the stage, and the GPU keeps drawing the same two buffers with a growing
lerp.

After an update that covers several ticks the previous buffer is all of them behind, so the timer's own `lerp` (how far
its accumulator is into the next tick) isn't the factor to draw with. `tfx_GetStageInterpolation(stage, timer.lerp,
timer.update_length)` turns it into the one that is, from the time the stage's last write actually covered, which
also allows for `max_frame_length` clamping a long frame. `tfx_impl_SetUniforms` does this. Without it, particles
stutter at high update rates and jump after a catch-up frame.

Rules the loop depends on:
- **Threading.** `tfx_UpdateStage` runs on its own thread. Add effects, move them and set the camera before it, and
  after it only read through the getters that wait for the update (the ones `tfx_impl_GatherFrame` uses), or call
  `tfx_CompleteStageWork` first.
- **Every GPU particle fetch is dispatched.** Each fetch carries on from the particle state the last dispatch left in
  the pool, so never call `tfx_GetGPUParticleFrame` on a frame you won't dispatch.
- **Effect ids are generational.** Keep the id `tfx_AddEffectTemplateToStage` returns and check it with
  `tfx_EffectIDIsValid`. A fire and forget effect's id goes stale once it frees itself.
- **Sync before updating stages.** Call `tfx_UpdateGlobalGPUData` before `tfx_UpdateStage`, on the thread that loads
  and frees libraries, and never while another thread is loading one. Upload the global data on the frames you upload
  instances, because the instances drawn between updates were written with the layout you uploaded last.

**Layers that write depth and lit layers.** Every example draws a layer whose `tfx_GetLayerFlags` has
`tfxLayerFlags_writes_depth` (layer 1 included) or `tfxLayerFlags_lit` as the editor does. A backend that draws every
layer with one blended, unlit pipeline still shows them, like Default. To add them to a backend, follow
`tfx_impl_dx11_Render`:
- Pick each layer's shaders from `frame.layer_flags`, which `tfx_impl_GatherFrame` reads from the library you pass it.
  The stage has one range per layer whichever library an effect came from, so keep every library's custom layers the
  same. The vertex shaders go by lit and the fragments by writes depth and lit: `tfx_impl_LayerWritesDepth` and
  `tfx_impl_LayerIsLit` give the indexes, and `tfx_impl_billboard_vertex_entry_points`,
  `tfx_impl_billboard_fragment_entry_points` and the ribbon pair name the entry points at each. On an API with pipeline
  objects that's four pipelines per draw type, plus the baked billboards' one.
- Draw in two groups against the depth buffer the scene drew into: the layers that write depth first, with blending off
  and depth writes on, then the rest blended over them. Walk the ribbon dispatches once for each group with a fresh
  dispatch and skip the buckets whose `dispatch.layer` isn't in it.
- Baked billboards aren't grouped by layer across instances, so they all draw blended and unlit. Baked ribbon buckets
  are, so they take the cutout and lit shaders.
- No example draws the blend modes in `tfxLayerFlags_blend_mask` yet, so those layers draw with the usual blend. They
  need only a different color blend in the blended group's pipelines, see Blend modes in `../Shaders/readme.md`.
- The lit fragments read `light` in `tfx_impl_uniforms_t`, which is at `tfxVIEW_LIGHT_OFFSET` because `camera_position`
  fills the 12 bytes before it. `tfx_impl_SetUniforms` sets the editor's default preview light, so set your own after
  it. See Lit layers in `../Shaders/readme.md`. Where each example gives the lit fragments the view:

  | Backend | The lit fragments' view |
  | --- | --- |
  | D3D11 | The frame constants at `b0` in the pixel shader as well |
  | D3D12 (bindless) | A second constant buffer range over the heap in the root signature, at `b0` space 2 |
  | Vulkan (bindless) | Nothing extra: the uniform buffer array at binding 7 is visible to the fragment stage |
  | WebGPU | Binding 10 in every render layout, fragment visible, the same block as binding 5 |
  | SDL_gpu | Fragment uniform slot 0 (space 3), pushed once a frame, and one uniform buffer on the lit fragment shaders |
  | OpenGL | The frame constants at uniform binding 10 as well |
  | Metal | `setFragmentBytes` at fragment buffer 0. Depth writes are a second `MTLDepthStencilState` |

  The bindless lit ribbon pipelines (Vulkan, D3D12) also read the `normal` attribute at location 4.

**Drawing by sort key.** `tfx_impl_GatherFrame` gathers each layer as one range, which is all you need while every
effect draws the same way. To draw some effects with their own pipeline or resources, give them a key with
`tfx_SetEffectSortKey` and walk each layer's keys yourself after the gather:

```c
for (tfxU32 layer = 0; layer != tfx_GetLayerCount(library); ++layer) {
	tfx_instance_key_range_t range;
	while (tfx_GetNextInstanceKeyRange(stage, layer, &range)) {
		//Bind what range.sort_key draws with, then draw the CPU instances and after them the GPU ones
		//draw(range.instance_start, range.instance_count);
		//draw(range.gpu_instance_start, range.gpu_instance_count);
	}
	tfx_ResetInstanceBufferLoopIndex(stage);
}
```

A library defines up to `tfxLAYERS` (64) layers, so loop them up to `tfx_GetLayerCount`, as `tfx_impl_GatherFrame` does
with the library you pass it. A stage matches layers by index whichever library an effect came from, so give libraries
that share a stage the same layer table. The editor's Copy Layers, in the import window, does that.

Ribbons aren't in those ranges. Each ribbon bucket holds one key's ribbons and reports it in `dispatch.sort_key`, and
the buckets come in the order they were made, not by key. Drawing every billboard and then every ribbon, as the D3D11
example does, puts all ribbons over all billboards whatever their keys. That's right until a key has to draw over
another key's ribbons. Then draw key by key, each key's billboards in every layer and then its ribbon buckets. Walk every
bucket and draw only the matching ones, because the dispatch offsets count up through the whole walk:

```c
//keys: the sort keys you gave out, lowest first
for (tfxU32 key_index = 0; key_index != key_count; ++key_index) {
	tfxU32 key = keys[key_index];
	for (tfxU32 layer = 0; layer != tfx_GetLayerCount(library); ++layer) {
		tfx_instance_key_range_t range;
		while (tfx_GetNextInstanceKeyRange(stage, layer, &range)) {
			if (range.sort_key == key) {
				//Draw the CPU then the GPU instances as above
			}
		}
		tfx_ResetInstanceBufferLoopIndex(stage);
	}
	tfx_ribbon_dispatch_t dispatch = tfx_CreateRibbonBatchDispatch(&frame.ribbons.stage_offsets);
	while (tfx_NextRibbonDispatch(stage, &dispatch)) {
		if (dispatch.sort_key == key) {
			//Push tfx_GetRibbonDispatchGlobals(&dispatch), then draw dispatch.index_count indexes from dispatch.index_offset
		}
	}
}
```

Every key switches from the billboard pipeline to the ribbon one and back, so keep this to the keys that need it.

## Direct3D 11 notes

- **No vertex buffers or input layout.** Under the fixed bindings the vertex shader reads the current and previous
  instance from raw buffers by `SV_InstanceID` and unpacks them itself. The quad corners come from
  `SV_VertexID`, with 6 vertices per instance.
- **`SV_InstanceID` doesn't include `StartInstanceLocation`**, so a draw's first instance goes in a constant buffer.
- **The push constants become constant buffers.** Every vertex shader takes the frame constants
  (`tfx_impl_uniforms_t`) at b0 and its draw's own constants at b1: the instance offset for billboards, the bucket
  globals for ribbons. Every other resource has a fixed register, so the bindless indexes are gone. The shaders are
  compiled with `-D TFX_FIXED_BINDINGS -D TFX_ENTRY_<entry point>`, which gives each entry point its own registers
  counting up from 0. Shaders/readme.md lists them.
- **Every buffer is raw.** They're created with `D3D11_RESOURCE_MISC_BUFFER_ALLOW_RAW_VIEWS` and given
  `DXGI_FORMAT_R32_TYPELESS` views, and the shaders read rows at their C size.
- **Struct layouts are checked.** `timelinefx.slang` mirrors the library's GPU structs by hand, so `tfx_impl_dx11.cpp`
  static_asserts their sizes and names the struct to fix when one changes.
- **Library data** goes up in `tfx_impl_dx11_UploadFrame` on update frames, each part only when its global version has
  changed. The properties, shape frames and graph lookups are `DEFAULT` buffers that grow by half again and are
  written with `UpdateSubresource`. The color ramps are a `DEFAULT` texture array with room for half as many layers
  again, and only the layers whose `tfx_GetGlobalColorRampLayerVersion` changed are written. It's recreated with every
  layer when the count outgrows it.
- **Instance buffers** are `DEFAULT` raw buffers, with UAVs when the GPU kernel writes them. The CPU instances go up in
  one `UpdateSubresource`, because the library lays the buffer out layer by layer. Each layer is then one
  `DrawInstanced(6, count)`. When the buffers grow, both copy their old contents across, since the previous buffer is
  still what `captured_index` points into.
- **Blending:** the pixel shaders write premultiplied alpha, so the blended layers use ONE / INV_SRC_ALPHA, with the
  depth test on and depth writes off. The layers that write depth draw before them with the default blend state, which
  has blending off, and a depth stencil state with writes on.
- **Shaders per layer:** `billboard_vertex_shaders[lit]` and `billboard_pixel_shaders[writes depth][lit]`, and the same
  for ribbons, set for each layer that has instances and each ribbon bucket. The frame constants are bound to the pixel
  shader's b0 too, where the lit fragments read the light.
- **Views are unbound after each pass.** D3D11 won't bind a buffer for writing while it's bound for reading, so every
  compute pass and draw unbinds what the next one writes. That keeps the debug layer quiet, which is the test every
  part of the example has passed.
- **Ribbons** are uploaded on update frames and rebuilt by `tfx_ribbon_compute` every frame, because relative ribbons follow
  their emitter with the frame's lerp. `tfx_impl_dx11_Compute` runs one dispatch per bucket with the bucket's
  `tfx_ribbon_bucket_globals_t` in a constant buffer. `tfx_impl_dx11_Render` then walks the same buckets again for the
  indexed draws. The compute shader writes the vertices to a raw buffer that the vertex shader pulls from, and
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
    up its returns gives a pool that fits exactly. `example_scene.cpp` does both. The capacities assume 60
    updates a second. At 20, or when one `tfx_UpdateStage` call catches up several ticks, rings drop spawns unless
    `tfx_SetGPUSizingUpdateRate` was set to the slowest rate first.
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

## Direct3D 12 notes

`tfx_impl_dx12` uses the bindless model rather than the fixed bindings, so it follows the D3D11 backend's functions and order
of work but binds nothing per draw: every buffer and texture is a view in one shader visible heap, and the shaders take each one's
heap slot from their push constants. The Shaders readme's "Bindless on D3D12" has the layout.

- **Frames in flight.** Call `tfx_impl_dx12_BeginFrame` with the frame slot once the GPU has finished the last frame that used it.
  Each slot has its own upload ring, and anything replaced while a frame may still read it (a grown buffer, its views, a big
  upload) is released when its slot comes round again. Upload, compute and render all record into your command list.
- **Every buffer a shader reads is in the default heap**, filled with `CopyBufferRegion` from the frame's upload ring, including
  the view uniforms and the kernel's upload. Buffers move between the copy, unordered access and one combined read state, and a
  buffer written again by the next frame gets a UAV barrier.
- **Root signatures.** One descriptor table of unbounded ranges, one per space, all at offset 0, a sampler table, and the push
  constant block as root constants. One root signature per block, because together they're over the 64 DWORD limit. Built with
  root signature version 1.0, where descriptors and data are volatile, because the heap gains views while frames are in flight.
- **The atlas is heap slot 0**, because `tfx_impl_GetUV` leaves the image index at 0. Mips are built on the CPU with the same 2 x 2
  box filter as the SDL_gpu backend, because D3D12 has no `GenerateMips`.
- **The billboard draws** use the instance buffer as a per instance vertex buffer and start each range with `DrawInstanced`'s
  start instance. The ribbon draws use the vertex and index buffers `tfx_ribbon_compute` wrote, which are also unordered access
  views.
- **Validation.** The example turns the debug layer on in debug builds and prints its messages every frame.
  `--test-debug-output` makes one deliberate error to prove they arrive, and `--gpu-validation` turns on GPU based validation,
  which checks what the shaders do with each descriptor. It reports live objects at exit, where only the device should be left.

## Vulkan notes

`tfx_impl_vulkan` uses the bindless model, the same SPIR-V the editor draws with. It follows the D3D12 backend's functions
and order of work, and binds one descriptor set for everything. The Shaders readme's "Bindless on Vulkan" has the layout.

- **Vulkan 1.3**, for dynamic rendering and synchronization2, plus the descriptor indexing features and
  `shaderDrawParameters` that `tfx_impl_vulkan_Init` lists in `tfx_impl_vulkan.h`. `example_glfw_vulkan` checks them
  when it picks a device and turns them on. The pipelines are made for the color and depth formats you pass, with no
  render pass. With a render pass of your own, give `VkGraphicsPipelineCreateInfo` its `renderPass` in place of the
  `VkPipelineRenderingCreateInfo`.
- **Frames in flight** work as on D3D12: call `tfx_impl_vulkan_BeginFrame` with the frame slot once its fence has
  signalled. Each slot has its own upload ring, and anything replaced while a frame may still read it (a grown buffer
  and its descriptor slot, a big upload) is destroyed when its slot comes round again.
- **Descriptors.** One set, with the sampler at binding 0, images at 3, storage buffers at 5 and uniform buffers at 7.
  Images and buffers are written into free slots while earlier frames are still in flight, so those two bindings are
  `UPDATE_AFTER_BIND`, `UPDATE_UNUSED_WHILE_PENDING` and `PARTIALLY_BOUND`. Storage buffers are untyped, so each buffer
  has one slot whatever row types it's read as: the pool's hot, cold and noise sections all take the pool's index with
  their bases, and the frame's upload arrays all take the upload's.
- **One pipeline layout** for every pipeline, with one push constant range of 112 bytes (the ribbon globals, the largest
  block) for the vertex, fragment and compute stages. The baked billboard draw pushes `tfx_baked_billboard_push`, whose
  first seven values are what the shared fragment shader reads. The kernel's block goes once per fetch, then only
  `dispatch_index`, its first 4 bytes, per dispatch.
- **Barriers are global memory barriers**, three of them each frame: one before the uploads' copies, one at the start
  of `tfx_impl_vulkan_Compute` and one at its end for the draws. So call `tfx_impl_vulkan_Compute` every frame, outside
  `vkCmdBeginRendering`, even when there's nothing to dispatch. An engine that tracks buffer state in a render graph can
  replace them with its own.
- **The view uniforms** are a small host visible uniform buffer per frame slot, written by `tfx_impl_vulkan_Render`,
  because copies aren't allowed inside rendering.
- **The atlas is image slot 0.** Its mips are linear blits of the level above, a 2 x 2 box filter. Every layer is moved
  to the shader read layout on the first upload, because the view covers layers no atlas has filled yet.
- **The billboard draws** take the instance buffer as a per instance vertex buffer, and `vkCmdDraw`'s first instance
  starts each range. The baked billboard pipeline leaves out location 7, `captured_index`, which its vertex shader
  doesn't read.
- **Memory.** Every buffer and image gets its own `vkAllocateMemory`, which is plenty for the few dozen the backend
  makes. An engine would sub-allocate from its own allocator or VMA.
- **Validation.** The example turns the validation layer on in debug builds when the Vulkan SDK has installed it, and
  prints its messages through a debug messenger. `--test-debug-output` creates a zero size buffer to prove they arrive,
  and `--gpu-validation` turns on GPU assisted validation. At exit the layer reports any object the backend leaked
  when the device is destroyed.

## WebGPU notes

`tfx_impl_webgpu` follows the SDL_gpu backend, with the same fixed bindings compiled to WGSL, and uses nothing but the
standard `webgpu.h`. It has run on wgpu-native v29 over Vulkan and D3D12. Dawn and the browser (Emscripten) haven't been
tried, and the example's window and device code is wgpu-native's.

- **Bind groups.** Everything is in group 0 at the bindings in the Shaders readme's slots table, with one explicit bind
  group layout per pipeline, since the uniform blocks take dynamic offsets. A bind group can't be changed once it's made,
  so each one is kept and made again only after a buffer or texture has been replaced, which a counter on the renderer
  tracks. The billboard draw and the kernel keep one for each instance buffer being the current one.
- **Uniforms.** WebGPU has no push constants, so every block lives in one uniform buffer at 256 byte steps and is bound
  at a dynamic offset: the view first, then each draw range's draw constants, then the compute blocks, one per kernel
  dispatch for its `dispatch_index`, the two playback blocks and each ribbon bucket's globals. The ribbon draw reads
  the globals again from the block its compute dispatch used. `tfx_impl_webgpu_Compute` and `tfx_impl_webgpu_Render` each
  write their blocks with one `wgpuQueueWriteBuffer`. A bind group holds the buffer itself, so the compute counts the
  ribbon buckets first and grows the buffer before taking any bind group.
- **Uploads go through the queue.** `wgpuQueueWriteBuffer` and `wgpuQueueWriteTexture` copy the data before returning
  and land after everything already submitted, so the backend has no upload rings, frame slots or deferred releases. A
  queue write lands before every command buffer of the next submit, though, so a copy recorded in the frame's encoder
  would run after the frame's uploads. Grown instance buffers are copied across in a command buffer of their own,
  submitted at once. `tfx_CopyGPUParticleUpload` fills a CPU copy that the queue copies again, and the pool needs no
  clearing because every new WebGPU buffer starts zeroed.
- **No barriers.** WebGPU puts them in between dispatches and passes. The compute is one compute pass, and the render
  pass binds nothing for writing.
- **Limits.** The defaults are enough: 6 storage buffers in a stage at most (the ribbon compute), read only in the
  vertex stage as WebGPU requires, and 256 invocations a workgroup. The example still asks for the adapter's
  `maxStorageBufferBindingSize` and `maxBufferSize`, because the default 128 MB binding can be too small for a big pool.
- **Textures.** Every view is made as a 2D array, even of a one layer texture, whose default view is plain 2D and doesn't
  fit the shaders' `texture_2d_array`. The color ramps are bound as unfilterable float, since they're only loaded, so a 32
  bit float ramp needs no `float32-filterable`. The atlas mips are built on the CPU, as for SDL_gpu, because WebGPU
  can't generate them.
- **Errors.** WebGPU validates every call and reports what no error scope catches to the device's uncaptured error
  callback, which the example prints. `tfx_impl_webgpu_Init` puts an error scope around the pipelines, so a shader or
  layout mistake fails Init. wgpu-native calls the scope back before `wgpuDevicePopErrorScope` returns. Where a WebGPU calls
  back later, the error is still printed but Init can't fail on it. In debug builds the example also turns on wgpu's
  backend validation (the Vulkan validation layer or the D3D12 debug layer) and prints wgpu-native's log.
  `--test-debug-output` creates a buffer that can't be both mapped and a storage buffer, and `--gpu-validation` turns on
  GPU based validation underneath.
- **D3D12 under wgpu** compiles naga's HLSL at pipeline creation, so the example asks for DXC by its file name,
  `dxcompiler.dll`, which it finds on the PATH (the Vulkan SDK has one). An empty path asks for a DXC built into
  wgpu-native, which the release builds don't have, and without the DLL wgpu leaves out D3D12 altogether. FXC works but
  takes much longer. The debug example starts in about 6 seconds on D3D12 and 3 on Vulkan.

## SDL_gpu notes

`tfx_impl_sdlgpu` follows the D3D11 backend call for call, with the same fixed bindings, so these are only the
differences.

- **The device.** SDL only turns on what you ask for, and the shaders need Vulkan 1.2 (`spirv_1_5`) and
  `shaderDrawParameters` (`SV_VertexID`). `CreateDevice` in `example_sdl3_gpu/main.cpp` passes both through
  `SDL_GPUVulkanOptions`. It takes `vulkan_core.h` from SDL's own source, because the Vulkan SDK's include folder has
  older SDL3 headers that hide `SDL_GPUVulkanOptions`. Without the options some drivers draw anyway, so only the
  validation layer tells you.
- **Validation output.** In debug mode SDL turns the validation layer on but installs no messenger, so nothing reaches
  you. The example sets `VK_KHRONOS_VALIDATION_DEBUG_ACTION=VK_DBG_LAYER_ACTION_LOG_MSG`, which prints the layer's
  messages to stdout (`VK_KHRONOS_VALIDATION_LOG_FILENAME` sends them to a file). D3D12's debug layer messages come
  through `SDL_Log`, which the example prints.
- **Shaders.** Each entry point is compiled twice with the Shaders readme's SDL `-D` line: SPIR-V with
  `-fvk-use-entrypoint-name`, so SDL finds it by name, and DXIL at shader model 6.0. SDL wants each shader's counts when
  it's created:

  | Entry point | Counts |
  | --- | --- |
  | `tfx_billboard_vertex` | 4 storage buffers, 2 uniform buffers |
  | `tfx_billboard_baked_vertex` | 3 storage buffers, 2 uniform buffers |
  | `tfx_ribbon_vertex` | 5 storage buffers, 2 uniform buffers |
  | `tfx_billboard_fragment`, `tfx_ribbon_fragment` | 2 samplers: the atlas and the color ramps |
  | `simulate_particles` | 2 read only, 2 read write, 1 uniform buffer, 256 threads |
  | `tfx_ribbon_compute` | 4 read only, 2 read write, 1 uniform buffer, 256 threads |
  | `tfx_sprite_playback`, `tfx_ribbon_playback` | 4 read only, 1 read write, 1 uniform buffer, 128 threads |

- **Uniforms** go through `SDL_PushGPUVertexUniformData` and `SDL_PushGPUComputeUniformData`: the frame's
  `tfx_impl_uniforms_t` at vertex slot 0 once per frame, each draw's first instance or bucket globals at slot 1, and
  each compute dispatch's block at compute slot 0, the kernel's before every dispatch.
- **Everything records into one command buffer** in call order: `tfx_impl_sdlgpu_UploadFrame` (one copy pass),
  `tfx_impl_sdlgpu_Compute` (compute passes), then `tfx_impl_sdlgpu_Render` inside your render pass. Upload and compute
  run whether or not you got a swapchain texture, because every GPU particle fetch has to be dispatched.
- **Uploads.** A frame's uploads are queued, copied into one transfer buffer mapped with cycling, and recorded in one
  copy pass, with grown instance buffers copied across first. Uploads into buffers never cycle, because the previous
  instance buffer and the particle pool have to keep what's in them. SDL doesn't zero buffers, so the pool's first
  upload is zeros.
- **Compute writes.** The kernel binds the pool and the instance buffer for writing without cycling, since the pool
  carries the particles over and the instance buffer already holds the frame's CPU instances. The playback and ribbon
  geometry passes rewrite everything they write, so they cycle and never wait on a frame still drawing.
- **Mips are built on the CPU.** `SDL_GenerateMipmapsForGPUTexture` does every layer, and its Vulkan and D3D12 mips
  differ from about the eighth level, which shows on particles a pixel or two across. A 2 x 2 box filter per layer
  matches D3D11's and is the same on both drivers.
- **D3D12 depth clears.** The depth texture is created with `SDL_PROP_GPU_TEXTURE_CREATE_D3D12_CLEAR_DEPTH_FLOAT` at 1,
  or the debug layer warns on every clear.

## OpenGL notes

`tfx_impl_gl` follows the D3D11 backend call for call, with the same fixed bindings, so these are only the differences.

- **The context.** 4.3 core, for storage buffers in the vertex stage and compute shaders, so not macOS. `tfx_impl_gl_Init`
  checks the version and loads every function it calls through your loader, which has to return the GL 1.1 functions
  too (`glfwGetProcAddress` and `SDL_GL_GetProcAddress` do, `wglGetProcAddress` alone doesn't).
- **Debug output.** With a debug context, Init installs a `glDebugMessageCallback` that prints through `tfx_Print`,
  leaving out notifications. Nothing else installs one, so ask for the debug context and prove the messages arrive.
- **Shaders.** Each entry point goes through SPIR-V and spirv-cross to `#version 430` GLSL, as the Shaders readme
  describes, with `TFX_COMBINED_SAMPLERS`. The backend compiles the text at start up and links the three render
  programs and four compute programs. Storage buffers, uniform blocks and textures keep the bindings from the slots table:
  GL numbers each of the three separately, so the defaults never collide.
- **Uniforms** are one small uniform buffer per kind of block, rewritten with `glBufferSubData` before each draw or
  dispatch: the frame's `tfx_impl_uniforms_t` at binding 5, a draw's first instance or a ribbon bucket's globals at 6, and
  each compute dispatch's block at 6.
- **Draws start at instance 0 and base vertex 0**, with the first instance in the draw constants. The vertex shaders
  subtract `gl_BaseInstanceARB` and `gl_BaseVertexARB` where the driver has `GL_ARB_shader_draw_parameters`, and a
  uniform left at 0 where it doesn't.
- **One empty vertex array.** Core profile draws need a vertex array bound, though nothing reads vertex attributes. The
  ribbon index buffer, which the compute writes as a storage buffer, is bound to it as the element buffer.
- **Barriers.** Storage buffer writes aren't visible to later draws without `glMemoryBarrier`. `tfx_impl_gl_Compute`
  issues one storage and element array barrier after all its dispatches, plus one between the ribbon playback and the
  ribbon compute that reads its output. The kernel's dispatches write separate rings and instance ranges, so they need
  none between them.
- **Textures.** The atlas is immutable storage (`glTexStorage3D`), so a `glTextureView` of an atlas's layers lets
  `glGenerateMipmap` rebuild only those, like the D3D11 backend's range view. The color ramps are read with
  `texelFetch`, but a texture still has to be complete for that, so they're created with one level and nearest filters.
- **Uploads** are `glBufferSubData` into buffers that grow by half again when they're too small, and grown instance
  buffers are copied across with `glCopyBufferSubData`. The kernel's upload is mapped with
  `GL_MAP_INVALIDATE_BUFFER_BIT` and filled by `tfx_CopyGPUParticleUpload`. The pool is created from zeroed memory, and
  checked against `GL_MAX_SHADER_STORAGE_BLOCK_SIZE`, which 4.3 only promises to 16 MB.
- **Depth writes.** `tfx_impl_gl_Render` turns them off for the particles and back on at the end, because `glClear`
  doesn't clear depth while they're off.

## Metal notes

`tfx_impl_metal` follows the WebGPU backend, with the same fixed bindings compiled to Metal source, and the Vulkan backend's frames
in flight. It's Objective-C++ built without ARC, so the renderer is a plain struct that releases what it owns, and its header needs
`<Metal/Metal.h>`. It has run on an Apple M2 under macOS 26. Intel and AMD Macs haven't been tried.

- **Shaders.** slangc writes each entry point as `.metal` text and `tfx_impl_metal_Init` compiles it with `newLibraryWithSource`,
  so the build needs no Metal toolchain, which recent Xcode versions download separately. All nine compile in about half a second,
  and an error points at its line in `timelinefx.slang`, because Slang keeps its `#line` directives.
- **Slots.** Slang makes HLSL registers `tN`, `uN` and `bN` all `[[buffer(N)]]` and ignores spaces, so under the fixed bindings
  a stage's classes would collide. `-D TFX_METAL_SLOTS` gives them their own tens: a vertex shader's buffers from 0, its view at
  10 and its draw block at 11; a compute pass's written buffers from 0, its read buffers from 10 and its block at 20. Textures
  and samplers have namespaces of their own and keep 0 and 1.
- **The kernel** takes its own Metal path whatever the defines, so it's compiled without them: its push constants at
  `[[buffer(0)]]` and a table of 16 `gpuAddress` values at `[[buffer(1)]]`, both sent with `setBytes`. The table holds the pool,
  the current instance buffer, the graph lookups and the frame's upload, and the push constants' bases do the rest. Metal only
  sees those buffers through `useResource`, read for the upload and the graph lookups, read and write for the pool and the
  instance buffer, which keeps them resident and ordered against the passes around them. `gpuAddress` needs Metal 3 (macOS 13),
  which Init checks.
- **Uniforms.** `setVertexBytes` and `setBytes` replace every uniform block: the view once per frame, each draw's first instance
  or each bucket's globals per draw, and each compute dispatch's block, the kernel's before every dispatch.
- **Frames in flight.** Call `tfx_impl_metal_BeginFrame` with the frame slot once the GPU has finished the last frame that used
  it. The example waits on a semaphore that its command buffers' completion handlers signal. Each slot has a shared upload ring.
  What the GPU keeps (library data, instances, ribbons, the baked data) is in private buffers filled by blits from the ring,
  because on Apple silicon a shared buffer is memory the GPU reads directly, and writing one that a frame in flight reads is a
  race. What a frame reads once (the kernel's upload, the baked playback's offsets and animation instances) is read in place
  from the ring. Command buffers keep what they use alive, so a grown private buffer is released at once, but a big upload's
  buffer of its own and an outgrown ring wait for their slot to come round, because later passes of the same frame still bind them.
- **Passes.** A frame's uploads are one blit pass, after a blit pass of their own that copies grown instance buffers across, so
  the new instances land on top. Compute is one concurrent compute pass: the kernel's dispatches, the playback passes and the
  ribbon buckets all write their own ranges, so its only barrier is `memoryBarrierWithScope` between the ribbon playback and
  the baked ribbon geometry that reads it. Metal tracks hazards between passes and command buffers, so there are no others.
  Wrap each frame in an `@autoreleasepool`, because encoders and drawables are autoreleased.
- **Draws** start at instance 0 with the first instance in the draw constants, as on the other fixed binding backends. Ribbons
  are drawn indexed from the index buffer `tfx_ribbon_compute` wrote.
- **Textures.** The atlas and the color ramps are private texture arrays. The atlas mips are built on the CPU with the same
  2 x 2 box filter as the WebGPU, SDL_gpu and D3D12 backends, straight into the upload memory after the layer they come from.
  The color ramps are only read, so they need no sampler.
- **The window.** GLFW with `GLFW_NO_API` and a `CAMetalLayer` on the Cocoa window's content view, BGRA8Unorm like the other
  examples, with its `drawableSize` set from the framebuffer on every resize. The depth texture is memoryless, since it's
  cleared and never stored.
- **Validation.** In debug builds the example sets `MTL_DEBUG_LAYER=1` and `MTL_DEBUG_LAYER_ERROR_MODE=nslog` before it creates
  the device, which works as well as setting them in the shell, and `--gpu-validation` sets `MTL_SHADER_VALIDATION=1`. Metal
  prints "Metal API Validation Enabled" and "Metal GPU Validation Enabled" when they take. API errors go to NSLog on stderr. A
  failed command buffer's error and shader validation's reports reach its completion handler, which prints them: shader
  validation leaves the command buffer completed and reports in its `logs`, with the line in `timelinefx.slang`. It also checks
  what the kernel reaches through its table, so a buffer left out of `useResource` is reported. `--test-debug-output` creates a
  buffer of zero length. With shader validation on the scene runs at about a quarter of its frame rate.
- **CMake** turns on OBJC as well as OBJCXX before anything adds GLFW, or a downloaded GLFW's `.m` files build as Objective-C++.

## sokol_gfx notes

There's no sokol backend here. An outside integration ported the D3D11 backend function for function onto sokol_gfx,
with one source building for `SOKOL_D3D11` and `SOKOL_GLCORE`. These are the differences it found.

- **Shaders.** D3D11 takes the DXBC line from the Shaders readme, GL takes the SPIR-V and SPIRV-Cross line, both with the
  fixed bindings. Each `sg_shader_desc` comes straight from the slots table: `VERTEX_READ` n is `views[n].storage_buffer`
  (read only, `hlsl_register_t_n` n, `glsl_binding_n` n), `COMPUTE_WRITE` n is `u`n and binding n, `COMPUTE_READ` n is `t`n
  and binding 2 + n, and the two fragment textures are t0 and t1, with GL names `tfx_fixed_combined_particle_images` and
  `tfx_fixed_combined_color_ramps`.
- **Uniforms.** Each block is a sokol uniform block at its b slot, applied with `sg_apply_uniforms` per draw and per
  dispatch. On GL, sokol sets each member with `glUniform*`, has no uint type and allows 16 members per block, so the GLSL
  is post-processed into `ivec4` arrays, as the Shaders readme's OpenGL notes describe, and fed as `SG_UNIFORMTYPE_INT4`.
  The size given to `sg_apply_uniforms` must match the declared one exactly.
- **The color ramps need a sampler.** sokol wants every texture in a texture/sampler pair on every backend, so declare
  a nearest sampler at slot 1 for them, which the D3D11 shader never reads.
- **Buffers.** Persistent buffers are `storage|copy_dst`, and the instance buffers `copy_src` too, so they can be copied
  across when they grow. Each frame's uploads go into one `write_transient` staging buffer and are copied out with
  `sg_copy_buffer_to_buffer` and `sg_copy_buffer_to_image` outside the passes. The kernel's upload can be a
  `write_transient` storage buffer itself, used in the same frame. Ribbon indexes are `storage|index_buffer`. The pool
  is created with zeroed data.
- **One compute pass a frame** holds the kernel's dispatches on update frames and the ribbon geometry every frame. sokol
  unbinds between passes, so none of the D3D11 backend's unbinding is needed.
- **Mips are built on the CPU** with a 2 x 2 box filter, as in the SDL_gpu backend: sokol has no `GenerateMips`.
- **Debug output.** sokol turns on the D3D11 debug layer and a GL debug context with `SOKOL_DEBUG`, but installs no GL
  callback. Load `glDebugMessageCallback` yourself and prove it reaches your log with a deliberate error.
- **No readback.** sokol_gfx can't read a buffer back. For checks, `sg_d3d11_query_buffer_info` and
  `sg_gl_query_buffer_info` give the native buffer to map or `glGetBufferSubData`.

## Clip space

`tfx_impl_Perspective` builds a right handed projection. Every backend passes the flags for its own clip space:

| API | Flags |
| --- | --- |
| D3D11, D3D12, Metal, WebGPU | `tfx_impl_clip_space_depth_zero_to_one` |
| Vulkan | `tfx_impl_clip_space_depth_zero_to_one \| tfx_impl_clip_space_flip_y` |
| OpenGL | none |

Pass the flags for what you draw through, not the API under it: a layer that gives every API one clip space flips
Vulkan for you (SDL_gpu does), so it takes the D3D row.

## Porting to a new API

Use the D3D11 backend as the template: keep its function names with your API's prefix, and its order of work.

The checks below are by eye, against `example_scene`. A port that plays its own scene, or can't see its window, has a
numeric version of each in `../agents/integration-checks.md`, and each check names its section.

**Setting up**
- [ ] Check the API can read storage buffers from the vertex stage. Both binding models need it, because the billboard
      vertex shader reads the previous update's instance at `captured_index`. WebGL2, most OpenGL ES 3.x devices and
      OpenGL below 4.3 (macOS included) can't.
- [ ] Add `example_<platform>_<api>/main.cpp`: a window, a device, a swap chain and the loop. It calls
      `example_scene_Load`, your upload functions, `example_scene_Start`, then `example_scene_Update` and your upload,
      compute and render calls each frame, with the clip space flags from the table above.
- [ ] Add `tfx_impl_<api>.h/.cpp` with the same functions as `tfx_impl_dx11.h`: `Init`, `CreateAtlas`, `UploadAtlas`,
      `UploadSpriteData`, `UploadFrame`, `Compute`, `Render` and `Shutdown`. The backend draws into whatever target
      `main.cpp` has bound, so it never owns the swap chain.
- [ ] Add a target to `CMakeLists.txt`, with a slangc build of each entry point your backend uses.
- [ ] Turn on the API's validation layer in debug builds and print what it says. Every step below is done when it's
      silent, so first prove its output reaches you with a deliberate error. Some layers (SDL_gpu, and sokol_gfx on
      OpenGL) turn validation on without installing a messenger or callback, and then silence means nothing is listening.
- [ ] Static_assert the sizes of every library struct your shaders mirror, as `tfx_impl_dx11.cpp` does.

**Shaders.** Never ported by hand: compile `../Shaders/timelinefx.slang` for your target, as `../Shaders/readme.md`
describes. Where the API has descriptor arrays (Vulkan, D3D12), use its default bindless model. Everywhere else, and
for any layer that binds a fixed set of resources per stage (SDL_gpu), compile with `-D TFX_FIXED_BINDINGS
-D TFX_ENTRY_<entry point>`: each entry point reads its own resources from slots counting up from 0, every shape
from one atlas, every buffer as raw rows, and the instances and ribbon vertexes from buffers in place of vertex
attributes. Macros set the set, space and binding base of each kind of slot, so the slots fit your API's or engine's
convention. D3D11 uses it with the defaults. OpenGL uses it through SPIR-V and SPIRV-Cross, as the Shaders readme
describes and the OpenGL backend does. WebGPU uses it as WGSL, as the WebGPU backend does. Metal uses it as Metal source
with `-D TFX_METAL_SLOTS`, as the Metal backend does.

**1. CPU billboards** (`example_scene_feature_none`)
- [ ] The atlas as an RGBA8 texture array with mips and room for every shape list's layers, and a clamped linear
      sampler with mip filtering. Upload each atlas into its own layers and make mips for just those layers, or for
      every layer if your API can't take a range (SDL_gpu can't): the other layers are rebuilt from their own pixels,
      which costs time but changes nothing. A layer over several APIs can make different mips on each (SDL_gpu's
      Vulkan and D3D12 do), which shows on particles a pixel or two across.
- [ ] The color ramps as a texture array, one layer per `tfx_GetGlobalColorRampLayer`, in the format of its bytes per
      pixel (always RGBA16F). It's read with an integer load, never sampled, but where textures and samplers are
      combined (OpenGL, SDL_gpu on Vulkan), or where your layer wants every texture paired (sokol_gfx), bind a nearest
      sampler with it. Upload only the layers whose version
      changed, and recreate it with headroom when `tfx_GetGlobalColorRampLayerCount` outgrows it.
- [ ] The shape frames (`tfx_GetGlobalGPUShapes`) and emitter properties (`tfx_GetGlobalParticlePropertiesBuffer`) as
      read only buffers, uploaded on update frames when their versions change.
- [ ] Two instance buffers that swap on every update, sized for the instance count and copied across when they grow.
- [ ] One instanced draw of 6 vertices per layer. If the API's instance id doesn't include the draw's first instance
      (D3D11, D3D12 and OpenGL's `gl_InstanceID`), pass the first instance in yourself. Vulkan's `gl_InstanceIndex`
      and WGSL's `instance_index` include it.
- [ ] Under the bindless model, use the vertex formats in the Shaders readme. `size` and `intensity_gradient_map` are
      `r16g16_sint`, which the shader converts to float.
- [ ] If you write your own billboard shader, decode `indexes` and `captured_index` as their comments in `tfx_instance_t`
      say: read the previous instance only when the capture flag (bit 15) is set, and draw the instance at zero size
      when it isn't: a particle's first frame, a loop back to the start of its path, or after too long a gap.
- [ ] If you write your own billboard or ribbon fragment shader, rotate the colour ramp texel's hue with `tfx_apply_hue`
      and the `hue_coefficients` varying. The billboard hue is `alignment.w * pi` radians, as `tfx_billboard_varyings` does.
- [ ] Premultiplied blending (ONE / ONE_MINUS_SRC_ALPHA), depth test on, depth write off.
- [ ] Check: the torch circles with its smoke trailing behind, the sparks stretch along their motion, the leaves tumble,
      and running at 20 and 400 updates a second looks the same as 60 (integration checks 1 and 5).
- [ ] Check: the Fireflies from the second library have their own shapes and colors. Swap the order the two libraries
      load in and they should look the same, because that moves every library's place in the global buffers
      (integration check 6).
- [ ] Check: press L to free and reload the second library while the first plays. The Fireflies come back with their
      own shapes and colors, nothing else flickers, and the validation layer stays silent. Check 6 of
      `../agents/integration-checks.md` has the procedure.

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
- [ ] Check: `tfx_GetGPUParticleStats` reports no dropped spawns at 60 updates a second (integration checks 2 to 4), and the leaves and the GPU
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

**5. Layers that write depth and lit layers** (any time after step 1, see "Layers that write depth and lit layers" in
[The frame loop](#the-frame-loop))
- [ ] The cutout and lit entry points compiled, and each layer's shaders or pipelines picked from `frame.layer_flags`.
- [ ] Two groups a frame: the layers that write depth first, blending off and depth writes on, then the rest.
- [ ] The frame uniforms, all `tfxVIEW_UNIFORMS_SIZE` bytes, bound to the lit fragments as well as the vertex shaders.
- [ ] Check: Space Weed's ribbons shade like tubes, bright where they face the light at the upper right and dark on the
      far side, and the leaves are opaque and lit. Flat ribbons mean the layer is drawn unlit, and black ones a zeroed
      light.

**Done when** the new example looks the same as the D3D11 one side by side, with a silent validation layer. Without
step 5, Space Weed and the leaves draw flat and blended.

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

## Level of Detail

Effects can turn on Level of Detail in the editor, which scales how many particles their emitters spawn, and how big they spawn, by distance from the camera. The library measures that distance from the camera position you pass to `tfx_SetStageCamera`, so set the camera every frame on any stage that plays effects with LoD. An effect at user spawn locations measures each location separately.

`tfx_SetStageLODDistanceScale(stage, scale)` multiplies every effect's LoD distances on the stage: 2 keeps full detail twice as far out, which suits a high quality setting or a zoomed in camera. 0 turns LoD off for the stage, and new stages start at 1. Call it between updates, as with `tfx_SetStageCamera`.

LoD only changes how many new particles spawn and their size at spawn, so nothing pops when the camera moves, and an effect spawning nothing in the distance stays alive and starts again when the camera comes back. Run on GPU emitters do less work when LoD thins them out, but their rings keep the same size, so LoD saves GPU time, not GPU memory.
