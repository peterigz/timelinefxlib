#version 450
#extension GL_ARB_separate_shader_objects : enable
#extension GL_EXT_nonuniform_qualifier : require

const float packed_scale_max = 128.0 / 32767.0;

const int tfx_intensity_offset                  = 0;
const int tfx_alpha_sharpness_offset            = 1;
const int tfx_curved_alpha_offset               = 2;
const int tfx_gradient_map_offset               = 3;
const int tfx_width_offset                      = 4;
const int tfx_scale_offset                      = 5;
const int tfx_uv_offset_y_offset                = 6;
const int tfx_uv_scale_y_offset                 = 7;
const int tfx_clip_offset_offset                = 8;
const int tfx_clip_size_offset                  = 9;
const int tfx_morph_amount_offset               = 10;
const int tfx_intensity_overlength_offset       = 14;
const int tfx_alpha_sharpness_overlength_offset = 15;
const int tfx_curved_alpha_overlength_offset    = 16;
const int tfx_gradient_map_overlength_offset    = 17;
const int tfx_width_overlength_offset           = 18;
const int tfx_fixed_angle_offset                = 19;
const int tfx_morph_bias_offset                 = 20;

layout(binding = 7) uniform UboView
{
    mat4 view;
    mat4 proj;
    vec2 screen_size;
	uint millisecs;
    float timer_lerp;
    float update_time;
} ub[];

struct ImageData {
	vec4 uv;
    vec2 padding;
	vec2 image_size;
	uint texture_array_index;
	float animation_frames;
    float max_radius;
};

struct tfx_ribbon {
    vec4 position;                      //normalised age is in w
    float width;
	uint segment_start_index;
	uint flags;
	uint padding_pre_quat;
	uvec2 quaternion;                   //16-bit snorm quaternion: .x = X|Y, .y = Z|W
    uint emitter_index;
	uint texture_indexes;
	uint intensity_gradient_map;		//Multiplier for the color of the ribbon
	uint curved_alpha;	    			//Sharpness and dissolve amount value for fading the image
	uint phase_seed;					//Ribbon uid, hash it to draw per ribbon phase offsets
	float scale;                        //Per frame overall_scale of the parent effect
};

const int tfx_lag_spine_samples = 16;

struct tfx_emitter {
	uvec2 quaternion;
	uint lookup_offset;
	uint angle_type;
	vec3 position;
	float age;
	vec3 captured_position;
	uint morph_segment_start_index;
	float noise_frequency;
	float noise_speed;
	float noise_phase_range;
	float noise_lock_rate;
	vec3 fixed_angle_normal;
	uint sample_count;
	uint noise_packed;
	float lag_time;
	float lag_span;
	uint padding;
	float lag_spine_position_x[tfx_lag_spine_samples];
	float lag_spine_position_y[tfx_lag_spine_samples];
	float lag_spine_position_z[tfx_lag_spine_samples];
	uvec2 lag_spine_quaternion[tfx_lag_spine_samples];
};

struct tfx_gpu_graph_data {
	vec4 node_data;
	vec4 oscillator;
	int easing_type;
	int flags;
	int padding[2];
};

layout(binding = 5) readonly buffer InImageData {
	ImageData data[];
} images[];

layout(binding = 5) readonly buffer InRibbonInstances {
	tfx_ribbon data[];
} ribbons[];

layout(binding = 5) readonly buffer InEmitters {
	tfx_emitter data[];
} emitters[];

layout (std430, binding = 5) readonly buffer InLookups {
    tfx_gpu_graph_data data[];
} graphs[];

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

layout(location = 0) in vec3 vertex_position;
layout(location = 1) in uint segment_index;
layout(location = 2) in vec2 uv_offset_scale;
layout(location = 3) in uint ribbon_index;

layout(location = 0) out vec3 out_tex_coord;
layout(location = 1) out ivec3 out_texture_indexes;
layout(location = 2) out vec4 out_intensity_curved_alpha_map;

vec2 unpack16bit_sscaled(uint packed) {
    int x_scaled = (int(packed) << 16) >> 16;
    int y_scaled = int(packed) >> 16;
    return vec2(float(x_scaled) * packed_scale_max, float(y_scaled) * packed_scale_max);
}

float tfx_easing_constant       (float t) { return 0;}
float tfx_easing_linear         (float t) { return t;}
float tfx_easing_in_cubic       (float t) { return t * t * t; }
float tfx_easing_out_cubic      (float t) { t = 1 - t; return 1 - (t * t * t); }
float tfx_easing_in_out_cubic   (float t) { float t3 = -2 * t + 2; t3 = t3 * t3 * t3; return t < 0.5 ? 4 * t * t * t : 1 - t3 * 0.5; }
float tfx_easing_smoothstep     (float t) { return t * t * (3.0 - 2.0 * t); }
float tfx_easing_out_in         (float t) { return t < 0.5 ? t * (2.0 - 2.0 * t) : 2.0 * (t - 0.5) * (t - 0.5) + 0.5; }

float tfx_apply_easing(int easing_type, float t) {
    switch(easing_type) {
        case 0:  return 0.0;                         // constant
        case 4:  return tfx_easing_in_cubic(t);      // in
        case 5:  return tfx_easing_out_cubic(t);     // out
        case 6:  return tfx_easing_in_out_cubic(t);  // in_out
        case 16: return t;                           // linear
        case 17: return tfx_easing_smoothstep(t);    // smoothstep
        case 18: return tfx_easing_out_in(t);        // out_in
        default: return t;                           // fallback to linear
    }
}

float tfx_do_bezier_lerp(float t, float first_node, float curve1, float curve2, float second_node) {
	float u = 1 - t;
	float w1 = u * u * u;
	float w2 = 3 * u * u * t;
	float w3 = 3 * u * t * t;
	float w4 = t * t * t;
	return w1 * first_node + w2 * curve1 + w3 * curve2 + w4 * second_node;
}

float tfx_do_linear_lerp(float t, float first_node, float second_node) {
	return first_node + (second_node - first_node) * t;
}

const int tfx_use_bezier_sampling_flag = 1;
const int tfx_enable_oscillator_flag   = 2;
const float tfx_two_pi = 6.283185307;

float tfx_graph_lookup(tfx_gpu_graph_data graph, float t) {
    float eased_t = tfx_apply_easing(graph.easing_type, t);
    float result;
    if ((graph.flags & tfx_use_bezier_sampling_flag) != 0) {
        result = tfx_do_bezier_lerp(eased_t, graph.node_data.x, graph.node_data.y, graph.node_data.z, graph.node_data.w);
    } else {
        result = tfx_do_linear_lerp(eased_t, graph.node_data.x, graph.node_data.w);
    }
    if ((graph.flags & tfx_enable_oscillator_flag) != 0) {
        //oscillator: .x = frequency, .y = amplitude, .z = offset_x, .w = offset_y
        //offset_x shifts the phase and offset_y lifts the multiplier, matching tfxOSCILLATOR_WIDE_APPLY.
        float oscillator_sin = 0.5 + sin((eased_t + graph.oscillator.z) * graph.oscillator.x * tfx_two_pi) * graph.oscillator.y;
        result = result * (oscillator_sin + graph.oscillator.w);
    }
    return result;
}

void main() {
	gl_Position = (ub[pc.uniform_index].proj * ub[pc.uniform_index].view * vec4(vertex_position, 1.0));
	gl_Position.x += pc.ndc_offset_x * gl_Position.w;
	gl_Position.y += pc.ndc_offset_y * gl_Position.w;
	float ribbon_position = float((segment_index & 0x007FF000) >> 12) / 2047.0;
	tfx_ribbon ribbon = ribbons[pc.ribbons_index].data[ribbon_index];
	tfx_emitter emitter = emitters[pc.emitters_index].data[ribbon.emitter_index];
	uint image_index = ribbon.texture_indexes & 0x00001FFF;

	//Calculate the uv coords across the width of the ribbon
	float tessellation = float((segment_index & 0xE0000000) >> 29);
	float tessellation_index = float((segment_index & 0x1C000000) >> 26);
	float side = float((segment_index & 0x02000000) >> 25);
	float t = tessellation_index / tessellation;  
	vec4 uv = images[pc.image_data_index].data[image_index].uv;
	float texture_uv_width = uv.z - uv.x;

	//Sample the shading graphs
	tfx_gpu_graph_data gradient_graph = graphs[pc.graphs_index].data[tfx_gradient_map_overlength_offset + emitter.lookup_offset];
	tfx_gpu_graph_data intensity_graph = graphs[pc.graphs_index].data[tfx_intensity_overlength_offset + emitter.lookup_offset];
	tfx_gpu_graph_data curved_alpha_graph = graphs[pc.graphs_index].data[tfx_curved_alpha_overlength_offset + emitter.lookup_offset];
	tfx_gpu_graph_data alpha_sharpness_graph = graphs[pc.graphs_index].data[tfx_alpha_sharpness_overlength_offset + emitter.lookup_offset];
	float sampled_gradient_map = tfx_graph_lookup(gradient_graph, ribbon_position);
	float sampled_intensity = tfx_graph_lookup(intensity_graph, ribbon_position);
	float sampled_curved_alpha = tfx_graph_lookup(curved_alpha_graph, ribbon_position);
	float sampled_alpha_sharpness = tfx_graph_lookup(alpha_sharpness_graph, ribbon_position);

	// In order to wrap the texture using images that are stored on in a texture array of sprite sheets we have to
	// ping pong the texture coordinates. This is the only solution I found to avoid the points along the ribbon
	// where the texture wraps and you end up with a small 1 ribbon segment sized texture that is flipped due to 
	// the wrapping of coordinates. Other alternatives mean you either have to add extra vertices when building the
	// ribbon or just using a separate texture array, neither of which I liked. Maybe there is another solution I'm 
	// overlooking though of course!
	float scaled_position = ribbon_position * uv_offset_scale.y + uv_offset_scale.x;
	float wrap_index = floor(scaled_position);		// Which wrap iteration we're on
	float wrap_fraction = fract(scaled_position);

	// If wrap_index is odd, reverse the direction (ping-pong)
	if(mod(wrap_index, 2.0) == 1.0) {
		wrap_fraction = 1.0 - wrap_fraction;  // Reverse direction for odd wraps
	}

	/*
	// Using multiply and abs to avoid branching(?) for future profiling experiments
	float ping_pong = abs(mod(wrap_index, 2.0) - wrap_fraction);
	// or alternatively:
	// float ping_pong = abs(fract(scaled_position * 0.5) * 2.0 - 1.0);
	*/

	float uv_x = wrap_fraction * texture_uv_width + uv.x;
	float uv_y = (t * 0.5 + (side * 0.5)) * (uv.w - uv.y) + uv.y;
	vec2 ribbon_intensity_gradient_map = unpack16bit_sscaled(ribbon.intensity_gradient_map);
	vec2 ribbon_curved_alpha = unpack16bit_sscaled(ribbon.curved_alpha);
	vec2 intensity_gradient_map = vec2(sampled_intensity, sampled_gradient_map);
	vec2 curved_alpha = vec2(sampled_curved_alpha, sampled_alpha_sharpness);

	int life = int(ribbon_position * 255);
	out_tex_coord = vec3(vec2(uv_x, uv_y), images[pc.image_data_index].data[image_index].texture_array_index);
	out_texture_indexes = ivec3((ribbon.texture_indexes & 0xFF000000) >> 24, (ribbon.texture_indexes & 0x00FF0000) >> 16, life);
	out_intensity_curved_alpha_map = vec4(intensity_gradient_map.x * ribbon_intensity_gradient_map.x, curved_alpha.x * ribbon_curved_alpha.x, curved_alpha.y * ribbon_curved_alpha.y, intensity_gradient_map.y * ribbon_intensity_gradient_map.y);
}