//Noise fields for the GPU particle kernel, scalar ports of the CPU's wide noise in timelinefx_internal.h. Every hash, prime and
//operation order follows the CPU so the two only differ by rounding.

import tfx_gpu_common;

static const float tfx_noise_epsilon = 1e-12;
static const float tfx_noise_offset = 100.0;			//tfxWIDENOISEOFFSET
static const float tfx_simplex_kernel_radius_squared = 0.5;	//tfxSIMPLEX_KERNEL_RADIUS_SQUARED
static const float tfx_curl_scale = 30.857;

//---- Value noise, tfx__wide_value_noise_3d

float3 tfx_value_noise_corner(uint lattice_x, uint lattice_y, uint lattice_z) {
	const uint h = tfx_seedgen((lattice_x + lattice_y) ^ lattice_z);
	return float3(float((h >> 20) & 0x3FFu), float((h >> 10) & 0x3FFu), float(h & 0x3FFu)) * (1.0 / 511.0) - 1.0;
}

float3 tfx_value_noise_lerp(float3 a, float3 b, float t) {
	return t * (b - a) + a;
}

float3 tfx_value_noise_3d(float3 position) {
	const float3 floor_position = floor(position);
	float3 fade = position - floor_position;
	fade = fade * fade * (3.0 - 2.0 * fade);

	const uint3 primes = uint3(1619u, 31337u, 6971u);
	const uint3 lattice0 = uint3(int3(floor_position)) * primes;
	const uint3 lattice1 = lattice0 + primes;

	float3 edge = tfx_value_noise_lerp(tfx_value_noise_corner(lattice0.x, lattice0.y, lattice0.z), tfx_value_noise_corner(lattice1.x, lattice0.y, lattice0.z), fade.x);
	const float3 slice = tfx_value_noise_lerp(edge, tfx_value_noise_lerp(tfx_value_noise_corner(lattice0.x, lattice1.y, lattice0.z), tfx_value_noise_corner(lattice1.x, lattice1.y, lattice0.z), fade.x), fade.y);
	edge = tfx_value_noise_lerp(tfx_value_noise_corner(lattice0.x, lattice0.y, lattice1.z), tfx_value_noise_corner(lattice1.x, lattice0.y, lattice1.z), fade.x);
	return tfx_value_noise_lerp(slice, tfx_value_noise_lerp(edge, tfx_value_noise_lerp(tfx_value_noise_corner(lattice0.x, lattice1.y, lattice1.z), tfx_value_noise_corner(lattice1.x, lattice1.y, lattice1.z), fade.x), fade.y), fade.z);
}

//---- Simplex noise, tfx__wide_simplex_cell and friends

//Named corners rather than arrays, so nothing is left for a driver to index dynamically or spill to scratch
struct tfx_simplex_cell {
	float3 offset0, offset1, offset2, offset3;
	uint hash0, hash1, hash2, hash3;
};

tfx_simplex_cell tfx_simplex_cell_at(float3 position) {
	const float skew = (position.x + (position.y + position.z)) * (1.0 / 3.0);
	const float3 cell = floor(position + skew);
	const float unskew = ((cell.x + cell.y) + cell.z) * (1.0 / 6.0);
	const float3 offset0 = position - (cell - unskew);

	const bool x_ge_y = offset0.x >= offset0.y;
	const bool y_ge_z = offset0.y >= offset0.z;
	const bool x_ge_z = offset0.x >= offset0.z;
	const uint3 first_step = uint3(x_ge_y && x_ge_z ? 1u : 0u, !x_ge_y && y_ge_z ? 1u : 0u, !(x_ge_z || y_ge_z) ? 1u : 0u);
	const uint3 second_step = uint3(x_ge_y || x_ge_z ? 1u : 0u, !x_ge_y || y_ge_z ? 1u : 0u, !(y_ge_z && x_ge_z) ? 1u : 0u);

	tfx_simplex_cell result;
	result.offset0 = offset0;
	result.offset1 = (offset0 - float3(first_step)) + (1.0 / 6.0);
	result.offset2 = (offset0 - float3(second_step)) + (1.0 / 3.0);
	result.offset3 = offset0 - 0.5;

	//floor() is already integral, so the truncating convert is exact for negative cells too
	const uint3 primes = uint3(1619u, 31337u, 6971u);
	const uint3 lattice = uint3(int3(cell)) * primes;
	const uint3 first_corner = lattice + first_step * primes;
	const uint3 second_corner = lattice + second_step * primes;
	const uint3 last_corner = lattice + primes;
	result.hash0 = tfx_seedgen((lattice.x + lattice.y) ^ lattice.z);
	result.hash1 = tfx_seedgen((first_corner.x + first_corner.y) ^ first_corner.z);
	result.hash2 = tfx_seedgen((second_corner.x + second_corner.y) ^ second_corner.z);
	result.hash3 = tfx_seedgen((last_corner.x + last_corner.y) ^ last_corner.z);
	return result;
}

float tfx_simplex_falloff(float3 offset) {
	return max(((tfx_simplex_kernel_radius_squared - offset.x * offset.x) - offset.y * offset.y) - offset.z * offset.z, 0.0);
}

//Perlin's 12 edge gradient, picked by the low 4 hash bits as tfx__wide_simplex_corner_edge does
float3 tfx_simplex_corner_edge(uint hash) {
	const bool u_is_x = (hash & 8u) == 0u;
	const bool v_is_y = (hash & 12u) == 0u;
	const bool v_is_x = (hash & 13u) == 12u;
	const bool v_is_z = !(v_is_y || v_is_x);
	const float u_signed_one = (hash & 1u) != 0u ? -1.0 : 1.0;
	const float v_signed_one = (hash & 2u) != 0u ? -1.0 : 1.0;
	return float3(u_is_x ? u_signed_one : (v_is_x ? v_signed_one : 0.0),
				  !u_is_x ? u_signed_one : (v_is_y ? v_signed_one : 0.0),
				  v_is_z ? v_signed_one : 0.0);
}

float tfx_simplex_corner(uint hash, float3 offset) {
	const float3 edge = tfx_simplex_corner_edge(hash);
	float falloff = tfx_simplex_falloff(offset);
	falloff *= falloff;
	falloff *= falloff;
	return falloff * (edge.x * offset.x + (edge.y * offset.y + edge.z * offset.z));
}

float tfx_simplex_noise_3d(float3 position) {
	const tfx_simplex_cell cell = tfx_simplex_cell_at(position);
	float sum = tfx_simplex_corner(cell.hash0, cell.offset0);
	sum += tfx_simplex_corner(cell.hash1, cell.offset1);
	sum += tfx_simplex_corner(cell.hash2, cell.offset2);
	sum += tfx_simplex_corner(cell.hash3, cell.offset3);
	return sum * 32.0;
}

//tfx__wide_simplex_noise_3d_gradient: without the x32, which tfx_curl_scale carries
//d/dp of falloff^4 * dot(edge, offset), which is falloff^4 * edge - 8 * falloff^3 * dot * offset
float3 tfx_simplex_corner_gradient(uint hash, float3 offset) {
	const float3 edge = tfx_simplex_corner_edge(hash);
	const float falloff = tfx_simplex_falloff(offset);
	const float falloff_squared = falloff * falloff;
	const float falloff_fourth = falloff_squared * falloff_squared;
	const float dot_edge = edge.x * offset.x + (edge.y * offset.y + edge.z * offset.z);
	const float offset_scale = (-8.0 * (falloff_squared * falloff)) * dot_edge;
	return falloff_fourth * edge + offset_scale * offset;
}

float3 tfx_simplex_noise_3d_gradient(float3 position) {
	const tfx_simplex_cell cell = tfx_simplex_cell_at(position);
	float3 gradient = tfx_simplex_corner_gradient(cell.hash0, cell.offset0);
	gradient += tfx_simplex_corner_gradient(cell.hash1, cell.offset1);
	gradient += tfx_simplex_corner_gradient(cell.hash2, cell.offset2);
	gradient += tfx_simplex_corner_gradient(cell.hash3, cell.offset3);
	return gradient;
}

//---- The four algorithms as a unit direction (curl is left unnormalised, as on the CPU)

float3 tfx_noise_normalise(float3 noise) {
	return noise * rsqrt(max(dot(noise, noise), tfx_noise_epsilon));
}

//tfx__wide_apply_noise_simplex_force: three decorrelated samples of one field
float3 tfx_simplex_noise_direction(float3 sample_position, float life) {
	const float shift = tfx_noise_offset + life;
	return tfx_noise_normalise(float3(tfx_simplex_noise_3d(sample_position),
									  tfx_simplex_noise_3d(float3(sample_position.x, sample_position.y + shift, sample_position.z)),
									  tfx_simplex_noise_3d(float3(sample_position.x, sample_position.y, sample_position.z + shift))));
}

//tfx__wide_apply_noise_curl_force: Bridson curl from the analytic gradients of three potentials
float3 tfx_curl_noise_direction(float3 sample_position, float life) {
	const float psi2_shift = tfx_noise_offset + life;
	const float psi3_shift = (tfx_noise_offset + tfx_noise_offset) + life;
	const float3 psi1 = tfx_simplex_noise_3d_gradient(sample_position);
	const float3 psi2 = tfx_simplex_noise_3d_gradient(sample_position + psi2_shift);
	const float3 psi3 = tfx_simplex_noise_3d_gradient(sample_position + psi3_shift);
	return float3(psi3.y - psi2.z, psi1.z - psi3.x, psi2.x - psi1.y) * tfx_curl_scale;
}

float3 tfx_value_noise_direction(float3 sample_position) {
	return tfx_noise_normalise(tfx_value_noise_3d(sample_position));
}

//tfx__wide_apply_white_noise_force's direction, already scaled by the inverse length the CPU folds into its amplitude
float3 tfx_white_noise_direction(uint uid, float age, float noise_resolution) {
	const float inv_period_scale = 1.0 / (250.0 * 300.0);
	const uint axis_x = 0x1b56c4f9u;
	const uint axis_y = 0x68e31da4u;
	const uint axis_z = 0xb5297a4du;
	const uint axis_phase = 0x27d4eb2fu;

	const float phase_offset = (tfx_white_unit(uid, axis_phase) * 0.5 + 0.5) * 256.0;
	const float phase = phase_offset + age * (noise_resolution * inv_period_scale);
	const uint bucket = uint(int(phase));
	float fraction = phase - float(int(bucket));
	fraction = fraction * fraction * (3.0 - 2.0 * fraction);

	const uint direction_seed = tfx_seedgen(bucket ^ uid);
	const uint next_direction_seed = tfx_seedgen((bucket + 1u) ^ uid);
	const float one_minus_fraction = 1.0 - fraction;
	const float3 noise = float3(tfx_white_unit(direction_seed, axis_x) * one_minus_fraction + tfx_white_unit(next_direction_seed, axis_x) * fraction,
								tfx_white_unit(direction_seed, axis_y) * one_minus_fraction + tfx_white_unit(next_direction_seed, axis_y) * fraction,
								tfx_white_unit(direction_seed, axis_z) * one_minus_fraction + tfx_white_unit(next_direction_seed, axis_z) * fraction);
	return noise * rsqrt(max(dot(noise, noise), tfx_noise_epsilon));
}

//tfx__wide_split_medium_velocity: steer and surge against the particle's heading by the two noise gains
float3 tfx_split_medium_velocity(float3 travel, float3 medium, float swerve_gain_minus_one, float speed_gain_minus_swerve_gain) {
	const float length_squared = dot(travel, travel);
	const float inverse_length = rsqrt(max(length_squared, 1e-8));
	travel *= inverse_length;
	const float movement_fade = length_squared * (inverse_length * inverse_length);
	const float swerve_gain = movement_fade * swerve_gain_minus_one + 1.0;
	const float along_scale = dot(medium, travel) * (movement_fade * speed_gain_minus_swerve_gain);
	return along_scale * travel + swerve_gain * medium;
}
