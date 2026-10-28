# Timeline FX C++ library

This is an alpha version (currently Alpha 0.38) of the TimelineFX library that you can use to incorporate particle effects directly into your games and applications. Effects are fully 3D, and for 2D games you just render them flat. You can use the [TimelineFX Particle Editor](https://www.rigzsoft.co.uk/timelinefx-alpha-version/) to design your effects.

**This is still very much a work in progress but feel free to use it and give some feedback on how to improve the interface and make it easier to use and implement.**

## Usage
The library is 1 cpp file and 2 header files, with no other dependencies:
* `timelinefx.h` - the public API. This is the only header your game needs to include.
* `timelinefx_internal.h` - everything internal to the library (allocator, containers, SIMD, threading, file IO). It's included by timelinefx.cpp, you don't need to include it yourself.
* `timelinefx.cpp` - the implementation.

Options for including the library in your project would be one of the following:
* Copy the three files into your project and `#include "timelinefx.h"`.
* Add the timelinefx folder to your list of Include Directories, add timelinefx.cpp to your build and `#include "timelinefx.h"`.
* Build timelinefx.cpp into your own static or dynamic library, link it and `#include "timelinefx.h"`.

The API is C compatible, so it can be used from C as well as C++.

The library is render agnostic, so you will have to provide your own means to render the particles and to load the images that particles use to render. The [Shaders](Shaders) folder has the exact shaders the editor uses to draw its preview, along with a [readme](Shaders/readme.md) explaining what each one does and how to port them to your own renderer.

The basic flow is:
1. `tfx_BeginTimelineFX` to start up the library and its worker threads.
2. `tfx_LoadEffectLibrary` to load a .tfx package (or a library saved as a folder), passing a callback to load the shape images into your renderer.
3. `tfx_CreateStage` to create a stage where your effects will be updated.
4. `tfx_CreateEffectTemplate` and `tfx_AddEffectTemplateToStage` to start an effect playing.
5. `tfx_UpdateStage` every frame, then upload the instance buffer (`tfx_GetInstanceBuffer`/`tfx_GetInstanceBufferByLayer`) and ribbon buffers and draw them.
6. `tfx_EndTimelineFX` when you're done.

Every function is documented in timelinefx.h, and there's a table of contents at the top so you can jump to each section.

## Features
* Simple C compatible API which is easy to use.
* Fast. Utilises multithreading and SIMD intrinsics (SSE4.1 minimum, optional AVX and FMA, NEON on ARM) to simulate particles as fast as possible. Work is always on going to make things faster.
* Particles are fully 3d, and for 2d games you can simply render everything flat. Along with particles there are ribbons for trails, beams, lightning and similar effects. Ribbons are built into geometry on the GPU with a compute shader and support noise, lag, morphing paths and clipping.
* Plenty of emission types: point, area, line, ellipse, disc, cylinder, icosphere, path (catmull-rom splines), spawning from another emitter's particles, spawning along a ribbon, and a user type where you place the particles yourself with a spawn callback. Emission can be distributed over a grid or stepped through evenly spaced angles (fans, rings and caps) rather than randomly.
* Forces on emitters: wind, attract/repel, vortex (with optional fluctuation), shockwave and noise (simplex, curl and value), each with a falloff profile graph and set in world, effect or emitter space.
* Simulate particles on the GPU. Point, user, area, ellipse, disc, cylinder and line emitters set to Run on GPU have their spawns generated in the kernel as well, the other emission types spawn on the CPU and upload their spawns. A compute shader simulates them, every emitter force included, and writes their instances after the CPU ones so the same shaders draw both. The kernel is a single Slang file, Shaders/tfx_gpu_particles.slang, which compiles for Vulkan, D3D12, D3D11, OpenGL, WebGPU and Metal, and you can add your own forces, collision or spawn placement to it from your own Slang module without changing it. See the GPU particle simulation section in timelinefx.h.
* User spawn locations. One running effect can spawn at a whole list of world positions that you add, move, rotate, pause and remove by id, so the engine trails of a whole fleet of ships can share a single effect instead of one effect each. Each location can have its own multipliers for life, size, velocity, spin, intensity and more, and forces work per location too.
* Effect templates let you tweak an effect from the library at runtime (scale graphs, change spawn amounts, disable emitters, set callbacks) without changing the base effect.
* Can pre-bake effects so that they can be animated entirely in a compute shader for ultimate speed. The baked data can also be compressed to save up to 90% memory depending on the effect. The compute shader does the job of interpolating the sprite data to keep things smooth between frames.
* Hot reload. `tfx_RefreshLibrary` checks the library on disk and merges any changes into a running game, restarting the effects that changed, and `tfx_RefreshShape` reloads a single image. The library can also live in a hot reloaded DLL: get and set the context with `tfx_GetContext`/`tfx_SetContext` and use `tfx_SuspendTimelineFX`/`tfx_ResumeTimelineFX` around the reload.
* Bookmarks. Place time markers on an effect in the editor and use `tfx_IsBookmarkCrossed` to time your own game events (sounds, screen shake, damage) off them.
* Warm up effects so that they start in a steady state rather than empty on the first frame.
* Fixed time step sub steps with `tfx_UpdateStageSubSteps` and `tfx_UpdateStageSubStep`.
* Instance buffers are grouped by layer, and effects can be given sort keys so that you can draw them in your own order or batch them by pipeline.
* Plug in your own memory allocation callbacks so that you own all of the memory the library uses.
* Color ramps can be stored as RGBA8, RGBA16F or RGBA32F, with the heat response baked into the ramp.
* Shape images can be PNG (rgba, luminance/alpha, luminance or alpha only) or block compressed KTX2 (BC4, BC5 and BC7) with mipmaps.
* Seeded random number generation so that effects can play out exactly the same every time.
* Optional profiling with [Tracy](https://github.com/wolfpld/tracy), define `tfxTRACY` to turn it on. Without it the profile markers compile away to nothing.
* Render agnostic, just a few simple integrations to get it drawing in your own renderer. (See the Shaders folder for example shaders)
* Runs on Intel/AMD and ARM CPUs, and is built and tested on Windows, Linux and MacOS.
* Many different types of effects are possible with more coming as features are added.

## Examples
The [implementations](implementations) folder has minimal examples of drawing TimelineFX with a graphics API directly, starting with Direct3D 11. Each API has one small backend, and they all share a renderer agnostic layer that loads the shapes, builds a texture atlas and gathers each frame's instances. They play the same scene, covering CPU particles, ribbons, GPU simulated particles and baked playback, and the [readme](implementations/README.md) walks through the setup order, the frame loop and a checklist for porting to another API.

There are more examples in [my own renderer, Zest](https://github.com/peterigz/zest) in the examples/SDL2 folder:
* zest-timelinefx - updating and drawing effects simulated in real time.
* zest-timelinefx-c - the same thing using the C API.
* zest-timelinefx-prerecorded-effects - drawing pre-baked effects with the compute shaders.
* zest-timelinefx-reload - running TimelineFX inside a DLL that is unloaded and reloaded while the simulation keeps running.

More graphics APIs (D3D12, Vulkan, OpenGL, WebGPU and Metal) will follow in the implementations folder, so watch this space.
