//Graph lookups for the GPU kernels. Mirrors tfx__sample_wide_graph: one segment, easing, then bezier or linear, then the oscillator

static const int tfx_graph_use_bezier_sampling = 1;
static const int tfx_graph_enable_oscillator = 2;

static const int tfx_easing_constant = 0;
static const int tfx_easing_in = 4;
static const int tfx_easing_out = 5;
static const int tfx_easing_in_out = 6;
static const int tfx_easing_smoothstep = 17;
static const int tfx_easing_out_in = 18;

static const float tfx_two_pi = 6.283185307;

//Mirrors tfx_gpu_graph_data_t, one per graph from tfx_GetGPUGraphLookupsBuffer
struct tfx_gpu_graph_data {
	float4 node_data;					//from, curve1, curve2, to
	float4 oscillator;					//frequency, amplitude, offset_x, offset_y
	int easing_type;
	int flags;
	int padding[2];
};

float tfx_apply_easing(int easing_type, float t) {
	switch (easing_type) {
		case tfx_easing_constant: return 0.0;
		case tfx_easing_in: return t * t * t;
		case tfx_easing_out: {
			float inverse = 1.0 - t;
			return 1.0 - inverse * inverse * inverse;
		}
		case tfx_easing_in_out: {
			float t3 = -2.0 * t + 2.0;
			t3 = t3 * t3 * t3;
			return t <= 0.5 ? 4.0 * t * t * t : 1.0 - t3 * 0.5;
		}
		case tfx_easing_smoothstep: return t * t * (3.0 - 2.0 * t);
		case tfx_easing_out_in: return t <= 0.5 ? t * (2.0 - 2.0 * t) : 2.0 * (t - 0.5) * (t - 0.5) + 0.5;
		default: return t;
	}
}

//Constant easing samples t = 0, which is why the CPU skips the bezier and the oscillator for it, and so does this
float tfx_sample_graph(tfx_gpu_graph_data graph, float t) {
	const float time = tfx_apply_easing(graph.easing_type, t);
	const bool constant = graph.easing_type == tfx_easing_constant;
	float value;
	if ((graph.flags & tfx_graph_use_bezier_sampling) != 0 && !constant) {
		const float inverse = 1.0 - time;
		const float weight1 = inverse * inverse * inverse;
		const float weight2 = 3.0 * inverse * inverse * time;
		const float weight3 = 3.0 * inverse * time * time;
		const float weight4 = time * time * time;
		value = weight1 * graph.node_data.x + weight2 * graph.node_data.y + weight3 * graph.node_data.z + weight4 * graph.node_data.w;
	} else {
		value = graph.node_data.x + (graph.node_data.w - graph.node_data.x) * time;
	}
	if ((graph.flags & tfx_graph_enable_oscillator) != 0 && !constant && graph.oscillator.x != 0.0 && graph.oscillator.y != 0.0) {
		const float oscillator = 0.5 + sin((time + graph.oscillator.z) * graph.oscillator.x * tfx_two_pi) * graph.oscillator.y;
		value *= oscillator + graph.oscillator.w;
	}
	return value;
}
