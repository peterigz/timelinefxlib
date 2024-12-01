# Timeline FX C++ library

This is an early alpha version of the TimelineFX library that you can use to incorporate 2D and 3D particle effects directly into your games and applications. You can use the [TimelineFX Particle Editor](https://www.rigzsoft.co.uk/timelinefx-alpha-version/) to design your effects.

**This is still very much a work in progress but feel free to use it and give some feedback on how to improve the interface and make it easier to use and implement.**

## Usage
The main library consists of 1 header file and 1 cpp file, there are no other dependencies. Options for including the library in your project would be one of the following three:
* Copy the timelinefx.h and .cpp files into your project and `#include "timelinefx.h"`.
* Add the timelinefx folder to your list of Includes Directories and `#include "timelinefx.h"`.
* Add TimelineFX.lib as a dependency and `#include "timelinefx.h"`

The library is render agnostic, so you will have to provide your own means to render the particles and to load the images that particles uses to render.

So far some of the features of the Library:
* Simple C compatible API which is easy to use.
* Fast. Utilises multithreading and SIMD intrinsics to simulate particles as fast as possible. Work is always on going to make things faster.
* Can pre-bake effects so that they can be animated entirely in a compute shader for ultimate speed. The baked data can also be compressed to save up to 90% memory depending on the effect. The compute shader does the job of interpolating the sprite data to keep things smooth between frames.
* Render both 2d or 3d particles.
* Render agnostic, just a few simple integrations to get it drawing in your own renderer. (See shaders folder for example shaders)
* Runs on Intel/ARM CPUs, currently tested on Windows and Mac. Should build on Linux but I don't actively test for that yet.
* Many different types of effects are possible with more coming as features are added.

I will build out proper examples of how to integrate this library with any renderer very soon so watch this space. In the meantime you will find some examples in by [own renderer here](https://github.com/peterigz/zest).
