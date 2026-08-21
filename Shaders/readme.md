# Shader Examples

Timelinefx is render agnostic, so you will have to decide how to integrate it into your own renderer. (I'll be looking to add integrations into popular renders over time).

In this folder are the exact shaders used to render the particles in the TimelineFX editor, so should be a good starting point. Shaders are currently all GLSL format.

## 2d and 3d vertex shaders
Both of the vertex shaders for 2d and 3d sprites use the same image.frag fragment shader.

Both of the vertex shaders are for rendering instanced sprites which is highly recommended for max speed. This means that you only need to pass the position, rotation and other sprite properties to the shader. The vertex data is then constructed inside the shader, no need for vertex/index buffers.

## Compute shaders for pre-baked effects

For maximum speed you can use a compute shader update the sprite instance buffer from pre-baked effect data. These compute shaders simply interpolate the sprite data from the current frame to the previous frame. This instance data is then passed to the exact same vertex shaders above.

[See here for an example of how to use the library to draw pre-baked effect data](https://github.com/peterigz/zest/blob/main/examples/GLFW/zest-timelinefx-prerecorded-effects/zest-timelinefx-prerecorded-effects.cpp)
