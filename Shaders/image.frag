#version 450
#extension GL_ARB_separate_shader_objects : enable

layout(location = 0) in vec4 in_frag_color;
layout(location = 1) in vec3 in_tex_coord;

layout(location = 0) out vec4 outColor;

layout(binding = 1) uniform sampler2DArray texSampler;

layout(push_constant) uniform quad_index
{
    mat4 model;
    vec4 parameters1;
	vec4 parameters2;
	vec4 parameters3;
	vec4 camera;
} pc;

void main() {
	vec4 texel = texture(texSampler, in_tex_coord);

	//Pre multiply alpha
	outColor.rgb = texel.rgb * in_frag_color.rgb * texel.a;
	//If in_frag_color.a is 0 then color will be additive. The higher the value of a the more alpha blended the color will be.
	outColor.a = texel.a * in_frag_color.a;

}
