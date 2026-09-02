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

[See here for an example of how to use the library to draw pre-baked effect data](https://github.com/peterigz/zest/blob/main/examples/GLFW/zest-timelinefx-prerecorded-effects/zest-timelinefx-prerecorded-effects.cpp)
