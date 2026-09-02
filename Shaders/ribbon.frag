#version 450
#extension GL_EXT_nonuniform_qualifier : require

layout(location = 0) in vec3 in_tex_coord;
layout(location = 1) in flat ivec3 in_color_ramp_coords;
layout(location = 2) in vec4 in_intensity_curved_alpha_map;

layout(location = 0) out vec4 out_color;

layout(binding = 3) uniform texture2DArray images[];
layout(binding = 0) uniform sampler samplers[];

layout(push_constant) uniform push_constants
{
    vec4 camera_position;
    uint segment_count;
    uint tessellation;
    uint index_offset;
    uint vertex_offset;
    uint ribbon_count;
    uint ribbon_offset;
    uint segment_offset;
	uint uniform_index;
	uint emitters_index;
	uint graphs_index;
	uint ribbons_index;
	uint ribbon_segments_index;
	uint vertexes_index;
	uint indexes_index;
	uint image_data_index;
	uint sampler_index;
	uint particle_texture_index;
	uint color_ramp_texture_index;
    float lerp;
    float time;
    float ndc_offset_x;
    float ndc_offset_y;
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
}