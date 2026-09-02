#version 450
#extension GL_EXT_nonuniform_qualifier : require

layout(location = 0) in vec3 in_tex_coord;
layout(location = 1) in flat ivec3 in_color_ramp_coords;
layout(location = 2) in vec4 in_intensity_curved_alpha_map;

layout(location = 0) out vec4 out_color;

layout(set = 0, binding = 3) uniform texture2DArray images[];
layout(set = 0, binding = 0) uniform sampler samplers[];

layout(push_constant) uniform quad_index
{
    uint particle_texture_index;
    uint color_ramp_texture_index;
    uint image_data_index;
    uint properties_index;
	uint sampler_index;
    uint prev_billboards_index;
    uint index_offset;
	uint uniform_index;
	uint disable_interpolation;
} pc;

void main() {
	vec4 texel = texture(sampler2DArray(images[nonuniformEXT(pc.particle_texture_index)], samplers[nonuniformEXT(pc.sampler_index)]), in_tex_coord);
	float lookup = clamp(texel.r * in_intensity_curved_alpha_map.w, 0.0, 1.0);
	int ramp_x = int(lookup * 255);
	ivec3 ramp = ivec3(ramp_x, in_color_ramp_coords.x, in_color_ramp_coords.y);
	vec4 ramp_texel = texelFetch(sampler2DArray(images[nonuniformEXT(pc.color_ramp_texture_index)], samplers[nonuniformEXT(pc.sampler_index)]), ramp, 0);
	ramp_texel *= in_intensity_curved_alpha_map.x;
	ramp_texel.a = min(1, ramp_texel.a);
	float curved_alpha = 1 - smoothstep(texel.a * in_intensity_curved_alpha_map.z, texel.a, 1 - in_intensity_curved_alpha_map.y);
	out_color.rgb = texel.rgb * ramp_texel.rgb * curved_alpha * texel.a;
	out_color.a = texel.a * ramp_texel.a * curved_alpha;
	//out_color.a = pow( texel.a * ramp_texel.a * curved_alpha, 2.0 );
	//out_color.rgb = ramp_texel.rgb * texel.a * curved_alpha;
	//out_color.a   = ramp_texel.a * texel.a * curved_alpha;
}