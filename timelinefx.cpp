#include "timelinefx.h"

/**
 * Computes the largest integer value not greater than the float one
 *
 * This method is faster than using (int32_t)std::floor(fp).
 *
 * I measured it to be approximately twice as fast:
 *  float:  ~18.4ns instead of ~39.6ns on an AMD APU),
 *  double: ~20.6ns instead of ~36.6ns on an AMD APU),
 * Reference: http://www.codeproject.com/Tips/700780/Fast-floor-ceiling-functions
 *
 * @param[in] fp    float input value
 *
 * @return largest integer value not greater than fp
 */
static inline int32_t fastfloor(float fp) {
	int32_t i = static_cast<int32_t>(fp);
	return (fp < i) ? (i - 1) : (i);
}

/**
 * Permutation table. This is just a random jumble of all numbers 0-255.
 *
 * This produce a repeatable pattern of 256, but Ken Perlin stated
 * that it is not a problem for graphic texture as the noise features disappear
 * at a distance far enough to be able to see a repeatable pattern of 256.
 *
 * This needs to be exactly the same for all instances on all platforms,
 * so it's easiest to just keep it as static explicit data.
 * This also removes the need for any initialisation of this class.
 *
 * Note that making this an uint32_t[] instead of a uint8_t[] might make the
 * code run faster on platforms with a high penalty for unaligned single
 * byte addressing. Intel x86 is generally single-byte-friendly, but
 * some other CPUs are faster with 4-aligned reads.
 * However, a char[] is smaller, which avoids cache trashing, and that
 * is probably the most important aspect on most architectures.
 * This array is accessed a *lot* by the noise functions.
 * A vector-valued noise over 3D accesses it 96 times, and a
 * float-valued 4D noise 64 times. We want this to fit in the cache!
 */
static const uint8_t perm[256] = {
	151, 160, 137, 91, 90, 15,
	131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23,
	190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
	88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166,
	77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244,
	102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196,
	135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123,
	5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
	223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
	129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228,
	251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107,
	49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
	138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
};

/**
 * Helper function to hash an integer using the above permutation table
 *
 *  This inline function costs around 1ns, and is called N+1 times for a noise of N dimension.
 *
 *  Using a real hash function would be better to improve the "repeatability of 256" of the above permutation table,
 * but fast integer Hash functions uses more time and have bad random properties.
 *
 * @param[in] i Integer value to hash
 *
 * @return 8-bits hashed value
 */
static inline uint8_t hash(int32_t i) {
	return perm[static_cast<uint8_t>(i)];
}

/* NOTE Gradient table to test if lookup-table are more efficient than calculs
static const float gradients1D[16] = {
		-8.f, -7.f, -6.f, -5.f, -4.f, -3.f, -2.f, -1.f,
		 1.f,  2.f,  3.f,  4.f,  5.f,  6.f,  7.f,  8.f
};
*/

/**
 * Helper function to compute gradients-dot-residual vectors (1D)
 *
 * @note that these generate gradients of more than unit length. To make
 * a close match with the value range of classic Perlin noise, the final
 * noise values need to be rescaled to fit nicely within [-1,1].
 * (The simplex noise functions as such also have different scaling.)
 * Note also that these noise functions are the most practical and useful
 * signed version of Perlin noise.
 *
 * @param[in] hash  hash value
 * @param[in] x     distance to the corner
 *
 * @return gradient value
 */
static float grad(int32_t hash, float x) {
	const int32_t h = hash & 0x0F;  // Convert low 4 bits of hash code
	float grad = 1.0f + (h & 7);    // Gradient value 1.0, 2.0, ..., 8.0
	if ((h & 8) != 0) grad = -grad; // Set a random sign for the gradient
//  float grad = gradients1D[h];    // NOTE : Test of Gradient look-up table instead of the above
	return (grad * x);              // Multiply the gradient with the distance
}

/**
 * Helper functions to compute gradients-dot-residual vectors (2D)
 *
 * @param[in] hash  hash value
 * @param[in] x     x coord of the distance to the corner
 * @param[in] y     y coord of the distance to the corner
 *
 * @return gradient value
 */
static float grad(int32_t hash, float x, float y) {
	const int32_t h = hash & 0x3F;  // Convert low 3 bits of hash code
	const float u = h < 4 ? x : y;  // into 8 simple gradient directions,
	const float v = h < 4 ? y : x;
	return ((h & 1) ? -u : u) + ((h & 2) ? -2.0f * v : 2.0f * v); // and compute the dot product with (x,y).
}

/**
 * Helper functions to compute gradients-dot-residual vectors (3D)
 *
 * @param[in] hash  hash value
 * @param[in] x     x coord of the distance to the corner
 * @param[in] y     y coord of the distance to the corner
 * @param[in] z     z coord of the distance to the corner
 *
 * @return gradient value
 */
static float grad(int32_t hash, float x, float y, float z) {
	int h = hash & 15;     // Convert low 4 bits of hash code into 12 simple
	float u = h < 8 ? x : y; // gradient directions, and compute dot product.
	float v = h < 4 ? y : h == 12 || h == 14 ? x : z; // Fix repeats at h = 12 to 15
	return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
}

/**
 * 1D Perlin simplex noise
 *
 *  Takes around 74ns on an AMD APU.
 *
 * @param[in] x float coordinate
 *
 * @return Noise value in the range[-1; 1], value of 0 on all integer coordinates.
 */
float SimplexNoise::noise(float x) {
	float n0, n1;   // Noise contributions from the two "corners"

	// No need to skew the input space in 1D

	// Corners coordinates (nearest integer values):
	int32_t i0 = fastfloor(x);
	int32_t i1 = i0 + 1;
	// Distances to corners (between 0 and 1):
	float x0 = x - i0;
	float x1 = x0 - 1.0f;

	// Calculate the contribution from the first corner
	float t0 = 1.0f - x0 * x0;
	//  if(t0 < 0.0f) t0 = 0.0f; // not possible
	t0 *= t0;
	n0 = t0 * t0 * grad(hash(i0), x0);

	// Calculate the contribution from the second corner
	float t1 = 1.0f - x1 * x1;
	//  if(t1 < 0.0f) t1 = 0.0f; // not possible
	t1 *= t1;
	n1 = t1 * t1 * grad(hash(i1), x1);

	// The maximum value of this noise is 8*(3/4)^4 = 2.53125
	// A factor of 0.395 scales to fit exactly within [-1,1]
	return 0.395f * (n0 + n1);
}

/**
 * 2D Perlin simplex noise
 *
 *  Takes around 150ns on an AMD APU.
 *
 * @param[in] x float coordinate
 * @param[in] y float coordinate
 *
 * @return Noise value in the range[-1; 1], value of 0 on all integer coordinates.
 */
float SimplexNoise::noise(float x, float y) {
	float n0, n1, n2;   // Noise contributions from the three corners

	// Skewing/Unskewing factors for 2D
	static const float F2 = 0.366025403f;  // F2 = (sqrt(3) - 1) / 2
	static const float G2 = 0.211324865f;  // G2 = (3 - sqrt(3)) / 6   = F2 / (1 + 2 * K)

	// Skew the input space to determine which simplex cell we're in
	const float s = (x + y) * F2;  // Hairy factor for 2D
	const float xs = x + s;
	const float ys = y + s;
	const int32_t i = fastfloor(xs);
	const int32_t j = fastfloor(ys);

	// Unskew the cell origin back to (x,y) space
	const float t = static_cast<float>(i + j) * G2;
	const float X0 = i - t;
	const float Y0 = j - t;
	const float x0 = x - X0;  // The x,y distances from the cell origin
	const float y0 = y - Y0;

	// For the 2D case, the simplex shape is an equilateral triangle.
	// Determine which simplex we are in.
	int32_t i1, j1;  // Offsets for second (middle) corner of simplex in (i,j) coords
	if (x0 > y0) {   // lower triangle, XY order: (0,0)->(1,0)->(1,1)
		i1 = 1;
		j1 = 0;
	}
	else {   // upper triangle, YX order: (0,0)->(0,1)->(1,1)
		i1 = 0;
		j1 = 1;
	}

	// A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
	// a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
	// c = (3-sqrt(3))/6

	const float x1 = x0 - i1 + G2;            // Offsets for middle corner in (x,y) unskewed coords
	const float y1 = y0 - j1 + G2;
	const float x2 = x0 - 1.0f + 2.0f * G2;   // Offsets for last corner in (x,y) unskewed coords
	const float y2 = y0 - 1.0f + 2.0f * G2;

	// Work out the hashed gradient indices of the three simplex corners
	const int gi0 = hash(i + hash(j));
	const int gi1 = hash(i + i1 + hash(j + j1));
	const int gi2 = hash(i + 1 + hash(j + 1));

	// Calculate the contribution from the first corner
	float t0 = 0.5f - x0 * x0 - y0 * y0;
	if (t0 < 0.0f) {
		n0 = 0.0f;
	}
	else {
		t0 *= t0;
		n0 = t0 * t0 * grad(gi0, x0, y0);
	}

	// Calculate the contribution from the second corner
	float t1 = 0.5f - x1 * x1 - y1 * y1;
	if (t1 < 0.0f) {
		n1 = 0.0f;
	}
	else {
		t1 *= t1;
		n1 = t1 * t1 * grad(gi1, x1, y1);
	}

	// Calculate the contribution from the third corner
	float t2 = 0.5f - x2 * x2 - y2 * y2;
	if (t2 < 0.0f) {
		n2 = 0.0f;
	}
	else {
		t2 *= t2;
		n2 = t2 * t2 * grad(gi2, x2, y2);
	}

	// Add contributions from each corner to get the final noise value.
	// The result is scaled to return values in the interval [-1,1].
	return 45.23065f * (n0 + n1 + n2);
}


/**
 * 3D Perlin simplex noise
 *
 * @param[in] x float coordinate
 * @param[in] y float coordinate
 * @param[in] z float coordinate
 *
 * @return Noise value in the range[-1; 1], value of 0 on all integer coordinates.
 */
float SimplexNoise::noise(float x, float y, float z) {
	float n0, n1, n2, n3; // Noise contributions from the four corners

	// Skewing/Unskewing factors for 3D
	static const float F3 = 1.0f / 3.0f;
	static const float G3 = 1.0f / 6.0f;

	// Skew the input space to determine which simplex cell we're in
	float s = (x + y + z) * F3; // Very nice and simple skew factor for 3D
	int i = fastfloor(x + s);
	int j = fastfloor(y + s);
	int k = fastfloor(z + s);
	float t = (i + j + k) * G3;
	float X0 = i - t; // Unskew the cell origin back to (x,y,z) space
	float Y0 = j - t;
	float Z0 = k - t;
	float x0 = x - X0; // The x,y,z distances from the cell origin
	float y0 = y - Y0;
	float z0 = z - Z0;

	// For the 3D case, the simplex shape is a slightly irregular tetrahedron.
	// Determine which simplex we are in.
	int i1, j1, k1; // Offsets for second corner of simplex in (i,j,k) coords
	int i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords
	if (x0 >= y0) {
		if (y0 >= z0) {
			i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 1; k2 = 0; // X Y Z order
		}
		else if (x0 >= z0) {
			i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 0; k2 = 1; // X Z Y order
		}
		else {
			i1 = 0; j1 = 0; k1 = 1; i2 = 1; j2 = 0; k2 = 1; // Z X Y order
		}
	}
	else { // x0<y0
		if (y0 < z0) {
			i1 = 0; j1 = 0; k1 = 1; i2 = 0; j2 = 1; k2 = 1; // Z Y X order
		}
		else if (x0 < z0) {
			i1 = 0; j1 = 1; k1 = 0; i2 = 0; j2 = 1; k2 = 1; // Y Z X order
		}
		else {
			i1 = 0; j1 = 1; k1 = 0; i2 = 1; j2 = 1; k2 = 0; // Y X Z order
		}
	}

	// A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
	// a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
	// a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
	// c = 1/6.
	float x1 = x0 - i1 + G3; // Offsets for second corner in (x,y,z) coords
	float y1 = y0 - j1 + G3;
	float z1 = z0 - k1 + G3;
	float x2 = x0 - i2 + 2.0f * G3; // Offsets for third corner in (x,y,z) coords
	float y2 = y0 - j2 + 2.0f * G3;
	float z2 = z0 - k2 + 2.0f * G3;
	float x3 = x0 - 1.0f + 3.0f * G3; // Offsets for last corner in (x,y,z) coords
	float y3 = y0 - 1.0f + 3.0f * G3;
	float z3 = z0 - 1.0f + 3.0f * G3;

	// Work out the hashed gradient indices of the four simplex corners
	int gi0 = hash(i + hash(j + hash(k)));
	int gi1 = hash(i + i1 + hash(j + j1 + hash(k + k1)));
	int gi2 = hash(i + i2 + hash(j + j2 + hash(k + k2)));
	int gi3 = hash(i + 1 + hash(j + 1 + hash(k + 1)));

	// Calculate the contribution from the four corners
	float t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0;
	if (t0 < 0) {
		n0 = 0.0;
	}
	else {
		t0 *= t0;
		n0 = t0 * t0 * grad(gi0, x0, y0, z0);
	}
	float t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1;
	if (t1 < 0) {
		n1 = 0.0;
	}
	else {
		t1 *= t1;
		n1 = t1 * t1 * grad(gi1, x1, y1, z1);
	}
	float t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2;
	if (t2 < 0) {
		n2 = 0.0;
	}
	else {
		t2 *= t2;
		n2 = t2 * t2 * grad(gi2, x2, y2, z2);
	}
	float t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3;
	if (t3 < 0) {
		n3 = 0.0;
	}
	else {
		t3 *= t3;
		n3 = t3 * t3 * grad(gi3, x3, y3, z3);
	}
	// Add contributions from each corner to get the final noise value.
	// The result is scaled to stay just inside [-1,1]
	return 32.0f*(n0 + n1 + n2 + n3);
}


/**
 * Fractal/Fractional Brownian Motion (fBm) summation of 1D Perlin Simplex noise
 *
 * @param[in] octaves   number of fraction of noise to sum
 * @param[in] x         float coordinate
 *
 * @return Noise value in the range[-1; 1], value of 0 on all integer coordinates.
 */
float SimplexNoise::fractal(size_t octaves, float x) const {
	float output = 0.f;
	float denom = 0.f;
	float frequency = mFrequency;
	float amplitude = mAmplitude;

	for (size_t i = 0; i < octaves; i++) {
		output += (amplitude * noise(x * frequency));
		denom += amplitude;

		frequency *= mLacunarity;
		amplitude *= mPersistence;
	}

	return (output / denom);
}

/**
 * Fractal/Fractional Brownian Motion (fBm) summation of 2D Perlin Simplex noise
 *
 * @param[in] octaves   number of fraction of noise to sum
 * @param[in] x         x float coordinate
 * @param[in] y         y float coordinate
 *
 * @return Noise value in the range[-1; 1], value of 0 on all integer coordinates.
 */
float SimplexNoise::fractal(size_t octaves, float x, float y) const {
	float output = 0.f;
	float denom = 0.f;
	float frequency = mFrequency;
	float amplitude = mAmplitude;

	for (size_t i = 0; i < octaves; i++) {
		output += (amplitude * noise(x * frequency, y * frequency));
		denom += amplitude;

		frequency *= mLacunarity;
		amplitude *= mPersistence;
	}

	return (output / denom);
}

/**
 * Fractal/Fractional Brownian Motion (fBm) summation of 3D Perlin Simplex noise
 *
 * @param[in] octaves   number of fraction of noise to sum
 * @param[in] x         x float coordinate
 * @param[in] y         y float coordinate
 * @param[in] z         z float coordinate
 *
 * @return Noise value in the range[-1; 1], value of 0 on all integer coordinates.
 */
float SimplexNoise::fractal(size_t octaves, float x, float y, float z) const {
	float output = 0.f;
	float denom = 0.f;
	float frequency = mFrequency;
	float amplitude = mAmplitude;

	for (size_t i = 0; i < octaves; i++) {
		output += (amplitude * noise(x * frequency, y * frequency, z * frequency));
		denom += amplitude;

		frequency *= mLacunarity;
		amplitude *= mPersistence;
	}

	return (output / denom);
}

namespace tfx {

	//these Variables determine the timing resolution that particles are updated at. So an Update frequency of 60 would mean that the particles are updated at 60 frames per second.
	float UPDATE_FREQUENCY = 60.f;
	float UPDATE_TIME = 1.f / UPDATE_FREQUENCY;
	float FRAME_LENGTH = 1000.f / UPDATE_FREQUENCY;

	//Set the udpate frequency for all particle effects - There may be options in the future for individual effects to be updated at their own specific frequency.
	void SetUpdateFrequency(float fps) {
		UPDATE_FREQUENCY = fps;
		UPDATE_TIME = 1.f / UPDATE_FREQUENCY;
		FRAME_LENGTH = 1000.f / UPDATE_FREQUENCY;
	}

	int FormatString(char* buf, size_t buf_size, const char* fmt, va_list args) {
		int w = vsnprintf(buf, buf_size, fmt, args);
		if (buf == NULL)
			return w;
		if (w == -1 || w >= (int)buf_size)
			w = (int)buf_size - 1;
		buf[w] = 0;
		return w;
	}

	void tfxText::Appendv(const char *format, va_list args) {
		va_list args_copy;
		va_copy(args_copy, args);

		int len = FormatString(NULL, 0, format, args);         // FIXME-OPT: could do a first pass write attempt, likely successful on first pass.
		if (len <= 0)
		{
			va_end(args_copy);
			return;
		}

		const int write_off = (string.size() != 0) ? string.size() : 1;
		const int needed_sz = write_off + len;
		if (write_off + (unsigned int)len >= string.capacity)
		{
			int new_capacity = string.capacity * 2;
			string.reserve(needed_sz > new_capacity ? needed_sz : new_capacity);
		}

		string.resize(needed_sz);
		FormatString(&string[write_off - 1], (size_t)len + 1, format, args);
		va_end(args_copy);

	}

	int tfxText::Find(const char *needle) {
		tfxText compare = needle;
		tfxText lower = Lower();
		compare = compare.Lower();
		if (compare.Length() > Length()) return -1;
		u32 pos = 0;
		int found = 0;
		while (compare.Length() + pos <= Length()) {
			if (strncmp(lower.string.data + pos, compare.string.data, compare.Length()) == 0) {
				return pos;
			}
			++pos;
		}
		return -1;
	}

	tfxText tfxText::Lower() {
		tfxText convert = *this;
		for (auto &c : convert.string) {
			c = tolower(c);
		}
		return convert;
	}

	void tfxText::AddLine(const char *format, ...) {
		va_list args;
		va_start(args, format);
		Appendv(format, args);
		va_end(args);
		Append('\n');
	}

	bool tfxText::SaveToFile(const char *file_name) {
		FILE *file = fopen(file_name, "wb");
		if (!file)
			return false;

		u32 l = Length();
		if (fwrite(string.data, 1, Length(), file) != Length())
			return false;

		fclose(file);
		return true;
	}

	const bool tfxText::IsInt() const {
		if (!Length()) return false;
		for (auto c : string) {
			if (c >= '0' && c <= '9');
			else {
				return false;
			}
		}
		return true;
	}

	const bool tfxText::IsFloat() const {
		if (!Length()) return false;
		int dot_count = 0;
		for (auto c : string) {
			if (c >= '0' && c <= '9');
			else if (c == '.' && dot_count == 0) {
				dot_count++;
			}
			else {
				return false;
			}
		}
		return dot_count == 1;
	}

	void tfxText::Appendf(const char *format, ...) {
		va_list args;
		va_start(args, format);

		va_list args_copy;
		va_copy(args_copy, args);

		int len = FormatString(NULL, 0, format, args);         // FIXME-OPT: could do a first pass write attempt, likely successful on first pass.
		if (len <= 0)
		{
			va_end(args_copy);
			return;
		}

		const int write_off = (string.size() != 0) ? string.size() : 1;
		const int needed_sz = write_off + len;
		if (write_off + (unsigned int)len >= string.capacity)
		{
			int new_capacity = string.capacity * 2;
			string.reserve(needed_sz > new_capacity ? needed_sz : new_capacity);
		}

		string.resize(needed_sz);
		FormatString(&string[write_off - 1], (size_t)len + 1, format, args_copy);
		va_end(args_copy);

		va_end(args);
	}

	tfxText tfxstream::ReadLine() {
		tfxText line;
		if (EoF()) return line;

		while (!EoF()) {
			char current = *(data + position);
			if (current == '\n') {
				position++;
				break;
			}
			line.Append(current);
			position++;
		}
		return line;
	}

	tfxPackage CreatePackage(const char *file_path) {
		tfxPackage package;
		package.header.magic_number = tfxMAGIC_NUMBER;
		package.header.flags = 0;
		package.header.offset_to_inventory = sizeof(tfxHeader);
		package.header.file_version = tfxFILE_VERSION;

		package.inventory.magic_number = tfxMAGIC_NUMBER_INVENTORY;
		package.inventory.entry_count = 0;

		package.file_path = file_path;
		return package;
	}

	bool ValidatePackage(tfxPackage &package) {
		if (package.header.magic_number != tfxMAGIC_NUMBER) return false;			//Package hasn't been initialised

		FILE *file = fopen(package.file_path.c_str(), "rb");
		if (!file) {
			return false;
		}

		if (file == NULL || _fseeki64(file, 0, SEEK_END)) {
			return false;
		}

		u64 length = _ftelli64(file);
		rewind(file);
		if (length == -1 || length >= SIZE_MAX) {
			fclose(file);
			return false;
		}

		if (length != package.file_size) return false;							//The file on disk is no longer the same size as the package file size since it was loaded

		//Everything seems ok
		return true;
	}

	void tfxEntryInfo::FreeData() {
		data.FreeAll();
	}

	tfxPackage::~tfxPackage() {
		Free();
	}

	tfxEntryInfo *tfxPackage::GetFile(const char *name) {
		if (!inventory.entries.ValidName(name))
			return nullptr;									//File not found in inventory
		assert(ValidatePackage(*this));						//The file on disk has changed since the package was loaded! Maybe this should return null instead?
		tfxEntryInfo *entry = &inventory.entries.At(name);
		if (entry->data.Size() != entry->file_size) {
			FILE *file = fopen(file_path.c_str(), "rb");
			_fseeki64(file, entry->offset_from_start_of_file, SEEK_SET);
			entry->data.Resize(entry->file_size);
			fread(entry->data.data, 1, entry->file_size, file);
			fclose(file);
		}
		return entry;
	}

	void tfxPackage::AddFile(tfxEntryInfo file) {
		inventory.entries.Insert(file.file_name, file);
		inventory.entry_count++;
	}

	void tfxPackage::AddFile(const char *file_name, tfxstream &data) {
		tfxEntryInfo entry;
		entry.file_name = file_name;
		entry.data = data;
		entry.file_size = data.size;

		inventory.entries.Insert(entry.file_name, entry);
		inventory.entry_count++;
	}

	void tfxPackage::Free() {
		for (auto &entry : inventory.entries.data) {
			entry.data.FreeAll();
		}
		inventory.entries.data.free_all();
		inventory.entries.map.free_all();
		file_data.FreeAll();
	}

	// Reads the whole file on disk into memory and returns the pointer
	tfxstream ReadEntireFile(const char *file_name, bool terminate) {
		tfxstream buffer;
		FILE *file = fopen(file_name, "rb");
		if (!file) {
			return buffer;
		}

		// file invalid? fseek() fail?
		if (file == NULL || fseek(file, 0, SEEK_END)) {
			return buffer;
		}

		long length = ftell(file);
		rewind(file);
		// Did ftell() fail?  Is the length too long?
		if (length == -1 || (unsigned long)length >= SIZE_MAX) {
			fclose(file);
			return buffer;
		}

		if (terminate)
			buffer.Resize(length + 1);
		else
			buffer.Resize(length);
		if (buffer.data == NULL || fread(buffer.data, 1, length, file) != length) {
			buffer.FreeAll();
			fclose(file);
			return buffer;
		}
		if (terminate)
			buffer.NullTerminate();

		fclose(file);
		return buffer;
	}

	bool SavePackageDisk(tfxPackage &package) {
		if (!package.file_path.Length()) return false;											//Package must have a file path
		if (package.header.magic_number != tfxMAGIC_NUMBER) return false;						//Header of package must contain correct magic number. Use CreatePackage to correctly initialise a package.
		if (package.inventory.magic_number != tfxMAGIC_NUMBER_INVENTORY) return false;			//Inventory of package must contain correct magic number

		FILE * file = fopen(package.file_path.c_str(), "wb");
		if (!file)
			return false;

		//Calculate the offset to the inventory which is stored at the end of the file after the contents
		u64 inventory_offset = sizeof(tfxHeader);
		for (auto &entry : package.inventory.entries.data) {
			entry.offset_from_start_of_file = inventory_offset;
			entry.file_size = entry.data.Size();
			inventory_offset += entry.data.Size();
		}

		//Sanity check, make sure that the entry count is the correct value
		package.inventory.entry_count = package.inventory.entries.Size();

		//Write the header, updating the inventory offset before hand
		package.header.offset_to_inventory = inventory_offset;
		fwrite((char*)&package.header, 1, sizeof(tfxHeader), file);

		//Write the file contents
		for (auto &entry : package.inventory.entries.data) {
			fwrite(entry.data.data, 1, entry.data.Size(), file);
		}

		//Write the inventory
		fwrite((char*)&package.inventory.magic_number, 1, sizeof(u32), file);
		fwrite((char*)&package.inventory.entry_count, 1, sizeof(u32), file);
		for (auto &entry : package.inventory.entries.data) {
			fwrite((char*)&entry.file_name.string.current_size, 1, sizeof(u32), file);
			fwrite(entry.file_name.c_str(), 1, entry.file_name.string.current_size, file);
			fwrite((char*)&entry.file_size, 1, sizeof(u64), file);
			fwrite((char*)&entry.offset_from_start_of_file, 1, sizeof(u64), file);
		}

		fclose(file);
		return true;
	}

	tfxstream SavePackageMemory(tfxPackage &package) {
		if (!package.file_path.Length()) return false;											//Package must have a file path
		if (package.header.magic_number != tfxMAGIC_NUMBER) return false;						//Header of package must contain correct magic number. CreatePackage to correctly initialise a package.
		if (package.inventory.magic_number != tfxMAGIC_NUMBER_INVENTORY) return false;			//Inventory of package must contain correct magic number

		//char *file = (char*)malloc(GetPackageSize(package));
		tfxstream file(GetPackageSize(package));
		if (!file.Size())
			return file;

		//Calculate the offset to the inventory which is stored at the end of the file after the contents
		u64 inventory_offset = sizeof(tfxHeader);
		for (auto &entry : package.inventory.entries.data) {
			entry.offset_from_start_of_file = inventory_offset;
			entry.file_size = entry.data.Size();
			inventory_offset += entry.data.Size();
		}

		//Sanity check, make sure that the entry count is the correct value
		package.inventory.entry_count = package.inventory.entries.Size();

		//Write the header, updating the inventory offset before hand
		package.header.offset_to_inventory = inventory_offset;
		file.Write(&package.header, sizeof(tfxHeader));

		//Write the file contents
		for (auto &entry : package.inventory.entries.data) {
			//fwrite(entry.data.data, 1, entry.data.Size(), file);
			file.Write(entry.data.data, entry.data.Size());
		}

		//Write the inventory
		file.Write(&package.inventory.magic_number, sizeof(u32));
		file.Write(&package.inventory.entry_count, sizeof(u32));
		for (auto &entry : package.inventory.entries.data) {
			file.Write(&entry.file_name.string.current_size, sizeof(u32));
			file.Write(entry.file_name.string.data, entry.file_name.string.current_size);
			file.Write(&entry.file_size, sizeof(u64));
			file.Write(&entry.offset_from_start_of_file, sizeof(u64));
		}

		return file;
	}

	u64 GetPackageSize(tfxPackage &package) {
		u64 space = 0;
		space += sizeof(tfxHeader);

		//Write the file contents
		for (auto &entry : package.inventory.entries.data) {
			space += entry.data.Size();
		}

		//Write the inventory
		space += sizeof(u32);
		space += sizeof(u32);
		for (auto &entry : package.inventory.entries.data) {
			space += sizeof(u32);
			space += entry.file_name.string.current_size;
			space += sizeof(u64);
			space += sizeof(u64);
		}

		return space;
	}

	int LoadPackage(const char *file_name, tfxPackage &package) {

		package.file_data = ReadEntireFile(file_name);
		if (package.file_data.Size() == 0)
			return tfxPackageErrorCode_unable_to_read_file;			//the file size is smaller then the expected header size

		package.file_size = package.file_data.Size();

		if (package.file_size < sizeof(tfxHeader))
			return tfxPackageErrorCode_wrong_file_size;				//the file size is smaller then the expected header size

		package.file_data.Read((char*)&package.header, sizeof(tfxHeader));

		if (package.header.magic_number != tfxMAGIC_NUMBER)
			return tfxPackageErrorCode_invalid_format;				//The header doesn't not contain the expected magic number "TFX!", incorrect file format;

		if (package.header.offset_to_inventory > package.file_size)
			return tfxPackageErrorCode_no_inventory;				//The offset to the inventory is beyond the size of the file

		package.file_data.Seek(package.header.offset_to_inventory);
		package.file_data.Read((char*)&package.inventory.magic_number, sizeof(u32));

		if (package.inventory.magic_number != tfxMAGIC_NUMBER_INVENTORY)
			return tfxPackageErrorCode_invalid_inventory;			//The value at the inventory offset does not equal the expected magic number "INV!"

		package.file_data.Read((char*)&package.inventory.entry_count, sizeof(u32));
		for (int i = 0; i != package.inventory.entry_count; ++i) {
			tfxEntryInfo entry;
			u32 file_name_size;
			package.file_data.Read((char*)&file_name_size, sizeof(u32));
			entry.file_name.string.resize(file_name_size);
			package.file_data.Read(entry.file_name.string.data, file_name_size);
			package.file_data.Read((char*)&entry.file_size, sizeof(u64));
			package.file_data.Read((char*)&entry.offset_from_start_of_file, sizeof(u64));
			package.inventory.entries.Insert(entry.file_name, entry);
		}

		package.file_path = file_name;

		return 0;
	}

	int LoadPackage(tfxstream &stream, tfxPackage &package) {
		//Note: tfxstream does not copy the memory, only the pointer, so if you FreeAll on the stream you pass in it will also free the file_data here as well
		package.file_data = stream;
		if (package.file_data.Size() == 0)
			return tfxPackageErrorCode_unable_to_read_file;			//the file size is smaller then the expected header size

		package.file_size = package.file_data.Size();

		if (package.file_size < sizeof(tfxHeader))
			return tfxPackageErrorCode_wrong_file_size;				//the file size is smaller then the expected header size

		package.file_data.Read((char*)&package.header, sizeof(tfxHeader));

		if (package.header.magic_number != tfxMAGIC_NUMBER)
			return tfxPackageErrorCode_invalid_format;				//The header doesn't not contain the expected magic number "TFX!", incorrect file format;

		if (package.header.offset_to_inventory > package.file_size)
			return tfxPackageErrorCode_no_inventory;				//The offset to the inventory is beyond the size of the file

		package.file_data.Seek(package.header.offset_to_inventory);
		package.file_data.Read((char*)&package.inventory.magic_number, sizeof(u32));

		if (package.inventory.magic_number != tfxMAGIC_NUMBER_INVENTORY)
			return tfxPackageErrorCode_invalid_inventory;			//The value at the inventory offset does not equal the expected magic number "INV!"

		package.file_data.Read((char*)&package.inventory.entry_count, sizeof(u32));
		for (int i = 0; i != package.inventory.entry_count; ++i) {
			tfxEntryInfo entry;
			u32 file_name_size;
			package.file_data.Read((char*)&file_name_size, sizeof(u32));
			entry.file_name.string.resize(file_name_size);
			package.file_data.Read(entry.file_name.string.data, file_name_size);
			package.file_data.Read((char*)&entry.file_size, sizeof(u64));
			package.file_data.Read((char*)&entry.offset_from_start_of_file, sizeof(u64));
			package.inventory.entries.Insert(entry.file_name, entry);
		}

		return 0;
	}

	EffectEmitter::~EffectEmitter() {
		sub_effectors.free_all();
	}

	void EffectEmitter::SoftExpire() {
		flags |= tfxEmitterStateFlags_stop_spawning;
	}

	void EffectEmitter::UpdateMaxLife() {
		max_life = GetMaxLife(*this);
		GetGraphByType(tfxOvertime_red)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_green)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_blue)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_blendfactor)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_intensity)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_velocity)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_width)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_height)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_weight)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_spin)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_stretch)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_spin)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_velocity_turbulance)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_direction_turbulance)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_velocity_adjuster)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_direction)->lookup.life = max_life;
	}

	void EffectEmitter::ResetAllBufferSizes() {
		tfxvec<EffectEmitter*> stack;
		stack.push_back(this);
		while (!stack.empty()) {
			EffectEmitter &current = *stack.pop_back();
			current.max_sub_emitters = 0;
			for (EachLayer) {
				current.max_particles[layer] = 0;
			}
			for (auto &sub : current.sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	void EffectEmitter::UpdateAllBufferSizes() {
		tfxvec<EffectEmitter*> stack;
		stack.push_back(this);
		for (int l = 0; l != tfxLAYERS; ++l) {
			max_particles[l] = 0;
		}
		max_sub_emitters = 0;
		while (!stack.empty()) {
			EffectEmitter &current = *stack.pop_back();
			if (current.type == tfxEmitterType) {
				unsigned int particle_count = 0;
				if (current.parent->parent && current.parent->parent->type != tfxFolder) {
					particle_count = current.GetHighestQty(current.parent->parent->max_life);
					max_sub_emitters += (current.parent->parent->max_particles[current.parent->parent->properties.layer] * (current.sub_effectors.size() + 1) + (current.sub_effectors.size() + 1)) * unsigned int(current.max_sub_emitter_life / 1000.f);
					if (max_sub_emitters > 10000) {
						int debug = 0;
					}
				}
				else {
					particle_count = current.GetHighestQty(0.f);
				}
				current.max_particles[current.properties.layer] = particle_count;
			}
			for (auto &sub : current.sub_effectors) {
				stack.push_back(&sub);
			}
		}
		//UpdateAllSpriteAmounts();
	}

	void EffectEmitter::UpdateAllSpriteAmounts() {
		for (int layer = 0; layer != tfxLAYERS; ++layer) {
			for (auto &emitter : sub_effectors) {
				max_particles[layer] += emitter.max_particles[layer];
				if (emitter.sub_effectors.size()) {
					max_particles[layer] += emitter.GetSubEffectSpriteCounts(layer, emitter.max_particles[emitter.properties.layer]);
				}
			}
		}
	}

	unsigned int EffectEmitter::GetSubEffectSpriteCounts(unsigned int layer, unsigned int multiplier) {
		float total = 0;
		for (auto &effect : sub_effectors) {
			for (auto &emitter : effect.sub_effectors) {
				if (properties.flags & tfxEmitterPropertyFlags_single && !(properties.flags & tfxEmitterPropertyFlags_one_shot))
					total += emitter.max_particles[layer] * multiplier;
				else
					total += emitter.max_particles[layer] * multiplier * (max_sub_emitter_life / 1000.f);
				if (emitter.sub_effectors.size()) {
					if (properties.flags & tfxEmitterPropertyFlags_single && !(properties.flags & tfxEmitterPropertyFlags_one_shot))
						total += emitter.GetSubEffectSpriteCounts(layer, emitter.max_particles[emitter.properties.layer] * multiplier);
					else
						total += emitter.GetSubEffectSpriteCounts(layer, emitter.max_particles[emitter.properties.layer] * multiplier * int(max_sub_emitter_life / 1000.f));
				}
			}
		}
		return (unsigned int)total;
	}

	float EffectEmitter::GetSubEffectLength() {
		if (!parent) return 0.f;
		float max_life = 0.f;
		for (auto &e : sub_effectors) {
			max_life = std::fmaxf(GetMaxLife(e), max_life);
		}
		max_life += GetMaxLife(*parent);
		return max_life;
	}

	unsigned int EffectEmitter::GetHighestQty(float parent_age) {
		max_sub_emitter_life = 0.f;
		if (!(properties.flags & tfxEmitterPropertyFlags_single) && !(properties.flags & tfxEmitterPropertyFlags_one_shot)) {
			float max_qty = 0.f;
			float amount_remainder = 0.f;
			float highest_age = 0.f;
			if (library->base_graphs[base].amount.nodes.size() > 1 || library->variation_graphs[variation].amount.nodes.size() > 1 || library->global_graphs[parent->global].amount.nodes.size() > 1
				|| library->base_graphs[base].life.nodes.size() > 1 || library->variation_graphs[variation].life.nodes.size() > 1 || library->global_graphs[parent->global].life.nodes.size() > 1) {
				tfxvec<tfxVec2> particles;
				float max_frames = std::fmaxf(library->base_graphs[base].amount.GetLastFrame(), library->variation_graphs[base].amount.GetLastFrame());
				max_frames = std::fmaxf(max_frames, library->base_graphs[base].life.GetLastFrame());
				max_frames = std::fmaxf(max_frames, library->variation_graphs[base].life.GetLastFrame());
				max_frames = std::fmaxf(max_frames, library->global_graphs[parent->global].amount.GetLastFrame());
				max_frames = std::fmaxf(max_frames, library->global_graphs[parent->global].life.GetLastFrame());
				max_frames = std::fmaxf(max_frames, properties.loop_length);
				float max_last_life = library->base_graphs[base].life.GetLastValue() + library->variation_graphs[variation].life.GetLastValue();
				float current_frame = 0.f;
				unsigned start_index = 0;
				if (parent_age > 0.f)
					max_frames = std::fminf(max_frames, parent_age);
				for (float frame = 0; frame <= max_frames + max_last_life; frame += FRAME_LENGTH) {
					float qty;
					qty = LookupFast(library->base_graphs[base].amount, current_frame);
					qty += LookupFast(library->variation_graphs[variation].amount, current_frame);
					float life = LookupFast(library->base_graphs[base].life, current_frame) + lookup_callback(library->variation_graphs[variation].life, current_frame);

					if (properties.flags & tfxEmitterPropertyFlags_use_spawn_ratio && (properties.emission_type == tfxArea || properties.emission_type == tfxEllipse)) {
						float area = LookupFast(library->property_graphs[property].emitter_width, current_frame) * LookupFast(library->property_graphs[property].emitter_height, current_frame);
						qty = (qty / 10000.f) * area;
					}
					else if (properties.flags & tfxEmitterPropertyFlags_use_spawn_ratio && properties.emission_type == tfxLine) {
						qty = (qty / 100.f) * LookupFast(library->property_graphs[property].emitter_height, current_frame);
					}

					qty *= LookupFast(parent->library->global_graphs[parent->global].amount, current_frame);
					qty *= UPDATE_TIME;
					qty += amount_remainder;

					float total_qty = 0;
					for (auto &p : particles) {
						if (parent_age && current_frame <= parent_age)
							highest_age = std::fmaxf(highest_age, p.y);
						p.y -= FRAME_LENGTH;
						if (p.y < 0) p.x = 0;
						total_qty += p.x;
					}
					max_qty = std::fmaxf(max_qty, total_qty);

					if (qty >= 1) {
						amount_remainder = qty - (int)qty;
						particles.push_back(tfxVec2(float((int)qty), life));
					}
					else {
						amount_remainder = qty;
					}
					current_frame += FRAME_LENGTH / tfxLOOKUP_FREQUENCY;
					if (properties.loop_length > 0 && current_frame > properties.loop_length)
						current_frame = 0;
				}

				float total_qty = 0;
				for (auto &p : particles) {
					p.y -= FRAME_LENGTH;
					if (parent_age && current_frame <= parent_age)
						highest_age = std::fmaxf(highest_age, p.y);
					if (p.y < 0) p.x = 0;
					total_qty += p.x;
				}
				max_qty = std::fmaxf(max_qty, total_qty);
				max_sub_emitter_life = highest_age + current_frame;

				return (unsigned int)max_qty;
			}
			else {
				float max_life = GetMaxLife(*this);
				max_sub_emitter_life = max_life + parent_age;
				float qty = (library->base_graphs[base].amount.GetFirstValue() + library->variation_graphs[variation].amount.GetFirstValue()) * (max_life / 1000.f) + 1.f;
				if (properties.flags & tfxEmitterPropertyFlags_use_spawn_ratio && (properties.emission_type == tfxArea || properties.emission_type == tfxEllipse)) {
					float area = library->property_graphs[property].emitter_width.GetFirstValue() * library->property_graphs[property].emitter_height.GetFirstValue();
					qty = (qty / 10000.f) * area;
				}
				else if (properties.flags & tfxEmitterPropertyFlags_use_spawn_ratio && properties.emission_type == tfxLine) {
					qty = (qty / 100.f) * library->property_graphs[property].emitter_height.GetFirstValue();
				}
				return unsigned int(qty);
			}
		}
		else {
			max_sub_emitter_life = max_life + parent_age;
			return (unsigned int)properties.spawn_amount;
		}
	}

	void tfxParticleMemoryTools::GetEffectMaxFrames(EffectEmitter &effect) {
		if (effect.type == tfxEffectType && effect.IsRootEffect()) {
			max_frames = std::fmaxf(max_frames, effect.library->global_graphs[effect.global].amount.GetLastFrame());
			max_frames = std::fmaxf(max_frames, effect.library->global_graphs[effect.global].life.GetLastFrame());
			max_frames = std::fmaxf(max_frames, effect.properties.loop_length);
		}
		else if (effect.type == tfxEmitterType) {
			max_frames = std::fmaxf(max_frames, effect.library->base_graphs[effect.base].amount.GetLastFrame());
			max_frames = std::fmaxf(max_frames, effect.library->variation_graphs[effect.base].amount.GetLastFrame());
			max_frames = std::fmaxf(max_frames, effect.library->base_graphs[effect.base].life.GetLastFrame());
			max_frames = std::fmaxf(max_frames, effect.library->variation_graphs[effect.base].life.GetLastFrame());
			max_frames = std::fmaxf(max_frames, effect.properties.loop_length);
			max_last_life = std::fmaxf(max_last_life, effect.library->base_graphs[effect.base].life.GetLastValue() + effect.library->variation_graphs[effect.variation].life.GetLastValue());
		}
		for (auto &sub : effect.sub_effectors) {
			GetEffectMaxFrames(sub);
		}
	}

	void tfxParticleMemoryTools::ProcessEffect(EffectEmitter &effect) {
		effects[0].reserve(10000);
		effects[1].reserve(10000);
		effects[0].clear();
		effects[1].clear();
		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			particles[layer][0].clear();
			particles[layer][1].clear();
			sprite_count[layer] = 0;
		}
		sub_effect_count = 0;
		emitters_removed = 0;
		GetEffectMaxFrames(effect);
		AddEffect(effect);
		initial_effect_size = effect.sub_effectors.size() + 1;
		Process();
		for (EachLayer) {
			effect.max_particles[layer] = sprite_count[layer];
			effect.max_sub_emitters = sub_effect_count;
		}
	}

	void tfxParticleMemoryTools::Process() {
		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			sprite_count[layer] = 0;
		}
		for (float frame = 0; frame <= max_frames + max_last_life; frame += FRAME_LENGTH) {
			current_effect.highest_particle_age += FRAME_LENGTH;
			for (int i = effects[current_buffer].current_size - 1; i >= 0; --i) {
				auto &emitter = effects[current_buffer][i];
				MockUpdateEmitter(emitter);
				if (emitter.library_link->type == tfxEmitterType && emitter.timeout_counter < emitter.timeout) {
					effects[!current_buffer].push_back(emitter);
				}
				else if (emitter.library_link->type == tfxEffectType && emitter.emitter_count > 0) {
					effects[!current_buffer].push_back(emitter);
				}
				else if (emitter.library_link->type == tfxEmitterType) {
					emitters_removed++;
				}
			}
			MockUpdateParticles();
			for (EachLayer) {
				sprite_count[layer] = (unsigned int)std::fmaxf((float)sprite_count[layer], (float)particles[layer][!current_buffer].current_size);
			}
			sub_effect_count = (unsigned int)std::fmaxf((float)effects[!current_buffer].current_size, (float)sub_effect_count);

			effects[current_buffer].clear();
			current_buffer = !current_buffer;
		}
		sub_effect_count -= initial_effect_size;
	}

	void tfxParticleMemoryTools::MockUpdateEmitter(tfxMockEffect &emitter) {
		if (emitter.library_link->type == tfxEmitterType) {

			for (auto p : emitter.particles[current_buffer]) {
				p -= FRAME_LENGTH;
				if (p > 0) emitter.particles[!current_buffer].push_back(p);
			}
			emitter.library_link->max_particles[emitter.library_link->properties.layer] = (unsigned int)std::fmaxf((float)emitter.library_link->max_particles[emitter.library_link->properties.layer], (float)emitter.particles[current_buffer].current_size);
			emitter.particles[current_buffer].clear();

			float life = LookupFast(emitter.library->base_graphs[emitter.library_link->base].life, emitter.frame) + lookup_callback(emitter.library->variation_graphs[emitter.library_link->variation].life, emitter.frame) * 0.75f;
			if (!(emitter.library_link->properties.flags & tfxEmitterPropertyFlags_single) && !(emitter.library_link->properties.flags & tfxEmitterPropertyFlags_one_shot)) {
				emitter.qty = LookupFast(emitter.library->base_graphs[emitter.library_link->base].amount, emitter.frame);
				emitter.qty += LookupFast(emitter.library->variation_graphs[emitter.library_link->variation].amount, emitter.frame);

				if (emitter.library_link->properties.flags & tfxEmitterPropertyFlags_use_spawn_ratio && (emitter.library_link->properties.emission_type == tfxArea || emitter.library_link->properties.emission_type == tfxEllipse)) {
					float area = LookupFast(emitter.library->property_graphs[emitter.library_link->property].emitter_width, emitter.frame) * LookupFast(emitter.library->property_graphs[emitter.library_link->property].emitter_height, emitter.frame);
					emitter.qty = (emitter.qty / 10000.f) * area;
				}
				else if (emitter.library_link->properties.flags & tfxEmitterPropertyFlags_use_spawn_ratio && emitter.library_link->properties.emission_type == tfxLine) {
					emitter.qty = (emitter.qty / 100.f) * LookupFast(emitter.library->property_graphs[emitter.library_link->property].emitter_height, emitter.frame);
				}

				emitter.qty *= LookupFast(emitter.library->global_graphs[emitter.library_link->parent->global].amount, emitter.frame);
				emitter.qty *= UPDATE_TIME;
				emitter.qty += emitter.amount_remainder;

				if (!emitter.started_spawning && emitter.qty < 1)
					emitter.qty = 1.f;

				if (emitter.qty >= 1) {
					while (emitter.qty > 1) {
						emitter.qty -= 1.f;
						particles[emitter.library_link->properties.layer][!current_buffer].push_back(life + FRAME_LENGTH);
						emitter.particles[!current_buffer].push_back(life + FRAME_LENGTH);
						for (auto &sub : emitter.library_link->sub_effectors) {
							AddEffect(sub);
						}
					}
					emitter.started_spawning = true;
					emitter.highest_particle_age = std::fmaxf(emitter.highest_particle_age, life);
				}
				emitter.amount_remainder = emitter.qty;
				if (emitter.library_link->properties.loop_length > 0 && emitter.frame > emitter.library_link->properties.loop_length)
					emitter.frame = 0;
			}
			else if(!emitter.single_shot_done) {
				emitter.started_spawning = true;
				emitter.library_link->max_particles[emitter.library_link->properties.layer] = emitter.library_link->properties.spawn_amount;
				for (int q = 0; q != emitter.library_link->properties.spawn_amount; ++q) {
					if (emitter.library_link->properties.flags & tfxEmitterPropertyFlags_single && !(emitter.library_link->properties.flags & tfxEmitterPropertyFlags_one_shot)) {
						particles[emitter.library_link->properties.layer][!current_buffer].push_back(max_frames + max_last_life);
						emitter.particles[!current_buffer].push_back(max_frames + max_last_life);
					}
					else {
						particles[emitter.library_link->properties.layer][!current_buffer].push_back(life + FRAME_LENGTH);
						emitter.particles[!current_buffer].push_back(life + FRAME_LENGTH);
					}
					emitter.highest_particle_age = std::fmaxf(emitter.highest_particle_age, life);
					for (auto &sub : emitter.library_link->sub_effectors) {
						AddEffect(sub);
					}
				}
				emitter.single_shot_done = true;
			}
		}
		else {
			if (emitters_removed > 0 && emitter.emitter_count > 0) {
				emitter.emitter_count--;
				emitters_removed--;
				if (emitter.emitter_count == 0)
					emitter.timeout_counter = emitter.timeout;
			}
			return;
		}

		emitter.age += FRAME_LENGTH;
		emitter.frame += FRAME_LENGTH / tfxLOOKUP_FREQUENCY;
		emitter.highest_particle_age -= FRAME_LENGTH;

		if (emitter.highest_particle_age <= 0 && emitter.started_spawning && emitter.qty == 0.f) {
			emitter.timeout_counter++;
		}
		else {
			emitter.timeout_counter = 0;
		}
	}

	void tfxParticleMemoryTools::MockUpdateParticles() {
		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			for (auto p : particles[layer][current_buffer]) {
				p -= FRAME_LENGTH;
				if (p > 0) particles[layer][!current_buffer].push_back(p);
			}
			particles[layer][current_buffer].clear();
		}
	}

	void tfxParticleMemoryTools::AddEffect(EffectEmitter &effect) {
		tfxMockEffect new_effect;
		new_effect.library_link = effect.library->GetEffect(effect.path_hash);
		new_effect.library = effect.library;
		new_effect.emitter_count = effect.sub_effectors.size();
		effects[!current_buffer].push_back(new_effect);
		for (auto &emitter : effect.sub_effectors) {
			tfxMockEffect new_emitter;
			new_emitter.library_link = effect.library->GetEffect(emitter.path_hash);
			new_emitter.library = effect.library;
			effects[!current_buffer].push_back(new_emitter);
		}
	}

	EffectEmitter& EffectEmitter::AddEmitter(EffectEmitter &e) {
		assert(e.name.Length());				//Emitter must have a name so that a hash can be generated
		e.type = EffectEmitterType::tfxEmitterType;
		e.library = library;
		e.uid = ++library->uid;
		sub_effectors.push_back(e);
		library->UpdateEffectPaths();
		ReIndex();
		return sub_effectors.back();
	}

	EffectEmitter& EffectEmitter::AddEffect(EffectEmitter &e) {
		assert(e.name.Length());				//Effect must have a name so that a hash can be generated
		e.type = EffectEmitterType::tfxEffectType;
		e.library = library;
		e.parent = this;
		e.uid = ++library->uid;
		sub_effectors.push_back(e);
		library->UpdateEffectPaths();
		ReIndex();
		return sub_effectors.back();
	}

	EffectEmitter& EffectEmitter::AddEffect() {
		EffectEmitter e;
		e.library = library;
		e.uid = ++library->uid;
		e.type = EffectEmitterType::tfxEffectType;
		e.name = "New Effect";
		sub_effectors.push_back(e);
		library->UpdateEffectPaths();
		ReIndex();
		return sub_effectors.back();
	}

	EffectEmitter &EffectEmitter::AddEffector(EffectEmitterType type) {
		EffectEmitter e;
		//e.parent_effect = this;
		e.type = type;
		e.library = library;
		e.uid = ++library->uid;
		if(e.type == tfxEffectType)
			e.name = "New Effect";
		else
			e.name = "New Emitter";
		sub_effectors.push_back(e);
		library->UpdateEffectPaths();
		ReIndex();
		return sub_effectors.back();
	}

	tfxEffect &GetEffect(tfxEffectPool &effect_pool, tfxEffectID &effect_id) {
		return effect_pool.GetEffect(effect_id);
	}

	bool ValidEffect(tfxEffectPool &effect_pool, tfxEffectID effect_id) {
		if (effect_pool.effect_memory.ranges.current_size <= effect_id)
			return false;
		tfxEffect &effect = effect_pool.GetEffect(effect_id);
		if (effect.id == effect_id)
			return true;
		return false;
	}

	void HardStopEffect(tfxEffectPool &effect_pool, tfxEffectID effect_id) {
		if (!ValidEffect(effect_pool, effect_id))
			return;
		tfxEffect &effect = GetEffect(effect_pool, effect_id);
		for (auto &emitter : effect.sub_emitters) {
			emitter.common.state_flags |= tfxEmitterStateFlags_stop_spawning;
			emitter.particles.clear();
		}
		for (auto &emitter : effect.sub_effects) {
			emitter.common.state_flags |= tfxEmitterStateFlags_stop_spawning;
			emitter.particles.clear();
		}
		effect.sub_effects.clear();
		for (EachLayer) {
			effect.sprites[layer].clear();
		}
	}

	void SoftStopEffect(tfxEffectPool &effect_pool, tfxEffectID effect_id) {
		if (!ValidEffect(effect_pool, effect_id))
			return;
		tfxEffect &effect = GetEffect(effect_pool, effect_id);
		for (auto &emitter : effect.sub_emitters) {
			emitter.common.state_flags |= tfxEmitterStateFlags_stop_spawning;
			emitter.common.state_flags |= tfxEmitterStateFlags_single_shot_done;
		}
		for (auto &emitter : effect.sub_effects) {
			emitter.common.state_flags |= tfxEmitterStateFlags_stop_spawning;
			emitter.common.state_flags |= tfxEmitterStateFlags_single_shot_done;
		}
	}

	void StartEffect(tfxEffectPool &effect_pool, tfxEffectID effect_id) {
		if (!ValidEffect(effect_pool, effect_id))
			return;
		HardStopEffect(effect_pool, effect_id);
		tfxEffect &effect = GetEffect(effect_pool, effect_id);
		for (auto &emitter : effect.sub_emitters) {
			emitter.common.state_flags &= ~tfxEmitterStateFlags_stop_spawning;
			emitter.common.state_flags &= ~tfxEmitterStateFlags_single_shot_done;
			emitter.common.age = 0.f;
			emitter.common.frame = 0.f;
			emitter.common.timeout_counter = 0.f;
		}
		effect.common.age = 0.f;
		effect.common.frame = 0.f;
		effect.common.timeout_counter = 0.f;
	}

	void ClearEffectPool(tfxEffectPool &effect_pool) {
		effect_pool.effect_memory.clear();
		effect_pool.emitter_memory.clear();
		effect_pool.sprite_memory.clear();
		effect_pool.particle_memory.clear();
		memset(effect_pool.effect_memory.data, tfxINVALID, effect_pool.effect_memory.capacity);
	}

	void FreeEffectPool(tfxEffectPool &effect_pool) {
		ClearEffectPool(effect_pool);
		effect_pool.effect_memory.free_all();
		effect_pool.emitter_memory.free_all();
		effect_pool.sprite_memory.free_all();
		effect_pool.particle_memory.free_all();
	}

	void UpdateEffect(tfxEffectPool &effect_pool, tfxEffectID effect_id) {
		if (effect_id > effect_pool.effect_memory.ranges.current_size)
			return;
		tfxEffect &e = effect_pool.GetEffect(effect_id);
		if (e.sub_emitters.empty()) 
			return;

		e.transform.captured = e.transform.world;

		if (e.lookup_mode == tfxPrecise) {
			e.common.frame = e.common.age;
		}
		else {
			e.common.frame = e.common.age / tfxLOOKUP_FREQUENCY;
		}


		float (*effect_lookup_callback)(Graph &graph, float age) = e.lookup_mode == tfxPrecise ? LookupPrecise : LookupFast;

		//Update the effect state
		e.current.life = effect_lookup_callback(e.common.library->global_graphs[e.library_link->global].life, e.common.frame);
		e.current.amount = effect_lookup_callback(e.common.library->global_graphs[e.library_link->global].amount, e.common.frame);
		if (!(e.common.property_flags & tfxEmitterPropertyFlags_global_uniform_size)) {
			e.current.size.x = effect_lookup_callback(e.common.library->global_graphs[e.library_link->global].width, e.common.frame);
			e.current.size.y = effect_lookup_callback(e.common.library->global_graphs[e.library_link->global].height, e.common.frame);
		}
		else {
			e.current.size.x = effect_lookup_callback(e.common.library->global_graphs[e.library_link->global].width, e.common.frame);
			e.current.size.y = e.current.size.x;
		}
		e.current.velocity = effect_lookup_callback(e.common.library->global_graphs[e.library_link->global].velocity, e.common.frame);
		e.current.spin = effect_lookup_callback(e.common.library->global_graphs[e.library_link->global].spin, e.common.frame);
		e.current.intensity = effect_lookup_callback(e.common.library->global_graphs[e.library_link->global].intensity, e.common.frame);
		e.current.splatter = effect_lookup_callback(e.common.library->global_graphs[e.library_link->global].splatter, e.common.frame);
		//We don't want to scale twice when the sub effect is transformed, so the values here are set to 1. That means that the root effect will only control the global scale.
		e.current.overal_scale = effect_lookup_callback(e.common.library->global_graphs[e.library_link->global].overal_scale, e.common.frame);
		if (!e.parent_particle) {
			e.transform.world.scale.x = e.current.overal_scale;
			e.transform.world.scale.y = e.current.overal_scale;
			e.transform.local.rotation = effect_lookup_callback(e.common.library->global_graphs[e.library_link->global].effect_angle, e.common.frame);
		}
		else {
			e.transform.world.scale.x = e.current.overal_scale;
			e.transform.world.scale.y = e.current.overal_scale;
			e.transform.local.rotation = 0.f;
		}
		e.current.stretch = effect_lookup_callback(e.common.library->global_graphs[e.library_link->global].stretch, e.common.frame);
		e.current.weight = effect_lookup_callback(e.common.library->global_graphs[e.library_link->global].weight, e.common.frame);

		if (!(e.common.state_flags & tfxEmitterStateFlags_retain_matrix)) {
			e.transform.world.position = e.transform.local.position;
			e.transform.world.rotation = e.transform.local.rotation;
			e.transform.world.position += e.common.handle * e.current.overal_scale;
			float s = sin(e.transform.local.rotation);
			float c = cos(e.transform.local.rotation);
			e.transform.matrix.Set(c, s, -s, c);
		}

		if (e.common.state_flags & tfxEmitterStateFlags_no_tween_this_update) {
			e.transform.captured = e.transform.world;
		}

		e.common.age += FRAME_LENGTH;
		if (!(e.common.property_flags & tfxEmitterPropertyFlags_single) || e.common.property_flags & tfxEmitterPropertyFlags_one_shot)
			e.common.highest_particle_age -= FRAME_LENGTH;

		if (e.common.loop_length && e.common.age > e.common.loop_length)
			e.common.age = 0;

		if (e.common.highest_particle_age <= 0 ) {
			e.common.timeout_counter += FRAME_LENGTH;
		}
		else {
			e.common.timeout_counter = 0;
		}

		e.common.state_flags &= ~tfxEmitterStateFlags_no_tween_this_update;

		//Update sub effects
		unsigned int index_offset = 0;
		for(int i = 0 ; i != e.sub_effects.current_size ; ++i) {
			tfxEmitter &sub_effect = e.sub_effects[i];
			if(sub_effect.type == tfxEffectType)
				sub_effect.UpdateAsSubEffect();
			else
				sub_effect.UpdateEmitter();

			if (sub_effect.common.timeout_counter >= sub_effect.common.timeout) {
				index_offset++;
				//sub_effect.next_emitter = nullptr;
				if (sub_effect.type == tfxEmitterType) {
					sub_effect.parent->common.active_children--;
					sub_effect.particles.free_range(e.storage->particle_memory);
				}
				continue;
			}
			if (index_offset > 0) {
				unsigned int tmp_offset = e.sub_effects[i - index_offset].offset;
				tfxEmitter *tmp_emitter = e.sub_effects[i - index_offset].next_emitter;
				sub_effect.offset = index_offset;
				sub_effect.next_emitter = &e.sub_effects[i - index_offset];
				e.sub_effects[i - index_offset] = sub_effect;
				e.sub_effects[i - index_offset].offset = tmp_offset;
				e.sub_effects[i - index_offset].next_emitter = tmp_emitter;
			}
			else {
				sub_effect.next_emitter = &sub_effect;
			}
		}
		e.sub_effects.current_size -= index_offset;

		//Update sub emitters
		for (auto &emitter : e.sub_emitters) {
			emitter.UpdateEmitter();
		}

		e.CompressSprites();
	}

	void tfxEmitter::UpdateAsSubEffect() {
		if (lookup_mode == tfxPrecise) {
			common.frame = common.age;
		}
		else {
			common.frame = common.age / tfxLOOKUP_FREQUENCY;
		}
		common.age += FRAME_LENGTH;

		if (common.loop_length && common.age > common.loop_length)
			common.age = 0;

		if (common.active_children == 0) {
			common.timeout_counter = common.timeout;
		}

		if (parent_particle && parent_particle->next_ptr) {
			parent_particle = parent_particle->next_ptr;
			current.overal_scale = common.root_effect->current.overal_scale;
			transform.world.scale.x = transform.world.scale.y = current.overal_scale;
			common.state_flags |= parent_particle->flags & tfxParticleFlags_remove;
			Transform(*this, *parent_particle);
			transform.world.position += common.handle * current.overal_scale;

			if (common.state_flags & tfxEmitterStateFlags_no_tween_this_update) {
				transform.captured = transform.world;
			}

		}
		else {
			parent_particle = nullptr;
			common.state_flags |= tfxEmitterStateFlags_retain_matrix;
			transform.local.position = transform.world.position;
			transform.local.rotation = transform.world.rotation;
			common.state_flags |= tfxEmitterStateFlags_stop_spawning;
		}
	}

	void tfxEmitter::UpdateEmitter() {
		transform.captured = transform.world;

		if (lookup_mode == tfxPrecise) {
			common.frame = common.age;
		}
		else {
			common.frame = common.age / tfxLOOKUP_FREQUENCY;
		}
		common.age += FRAME_LENGTH;

		common.property_flags = library_link->properties.flags;
		if (!(common.property_flags & tfxEmitterPropertyFlags_single) || common.property_flags & tfxEmitterPropertyFlags_one_shot || common.state_flags & tfxEmitterStateFlags_stop_spawning)
			common.highest_particle_age -= FRAME_LENGTH;

		if (common.loop_length && common.age > common.loop_length)
			common.age = 0;

		if (common.highest_particle_age <= 0 && current.qty == 0.f) {
			common.timeout_counter += FRAME_LENGTH;
		}
		else {
			common.timeout_counter = 0.f;
		}

		float (*effect_lookup_callback)(Graph &graph, float age) = common.root_effect->lookup_mode == tfxPrecise ? LookupPrecise : LookupFast;

		common.state_flags |= (common.root_effect->common.state_flags & tfxEmitterStateFlags_remove);
		transform.local.rotation = effect_lookup_callback(common.library->property_graphs[library_link->property].emitter_angle, common.frame);
		current.velocity_adjuster = effect_lookup_callback(common.library->overtime_graphs[library_link->overtime].velocity_adjuster, common.frame);
		current.overal_scale = common.root_effect->current.overal_scale;
		current.stretch = common.root_effect->current.stretch;

		current.emitter_size.y = effect_lookup_callback(common.library->property_graphs[library_link->property].emitter_height, common.frame);
		if (library_link->properties.emission_type == EmissionType::tfxArea || library_link->properties.emission_type == EmissionType::tfxEllipse) {
			current.emitter_size.x = effect_lookup_callback(common.library->property_graphs[library_link->property].emitter_width, common.frame);
		}
		else
			current.emitter_size.x = 0.f;
		if (common.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center && library_link->properties.emission_type != EmissionType::tfxPoint) {
			common.handle = current.emitter_size * -0.5f;
		}
		else if (!(common.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center)) {
			common.handle = library_link->properties.emitter_handle;
		}
		if (common.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center && library_link->properties.emission_type == EmissionType::tfxLine) {
			common.handle = current.emitter_size * 0.5f;
		}

		if(library_link->emitter_update_callback)
			library_link->emitter_update_callback(*this);
	
		if (parent) {
			parent = parent->next_emitter;
			Transform(*this, *parent);
		} else
			Transform(*this, *common.root_effect);

		if (common.state_flags & tfxEmitterStateFlags_no_tween_this_update) {
			transform.captured = transform.world;
		}

		SpawnParticles();
		
		ControlParticles();

		common.state_flags &= ~tfxEmitterStateFlags_no_tween_this_update;

	}

	void tfxEmitter::RefreshFromLibrary() {
		tfxEmitterStateFlags single_done = common.state_flags & tfxEmitterStateFlags_single_shot_done;
		tfxEmitterStateFlags spawn_done = common.state_flags & tfxEmitterStateFlags_stop_spawning;
		common.property_flags = library_link->properties.flags;
		common.state_flags = library_link->flags;
		common.loop_length = library_link->properties.loop_length;
		common.handle = library_link->properties.emitter_handle;
		if (library_link->properties.flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
			image_handle = tfxVec2(0.5f, 0.5f);
		}
		else {
			image_handle = library_link->properties.image_handle;
		}
		common.state_flags &= ~tfxEmitterStateFlags_retain_matrix;

		common.state_flags |= library_link->properties.flags & tfxEmitterPropertyFlags_single && !(library_link->properties.flags & tfxEmitterPropertyFlags_one_shot) ? tfxEmitterStateFlags_is_single : 0;
		common.state_flags |= (library_link->properties.emission_type != tfxLine && !(library_link->properties.flags & tfxEmitterPropertyFlags_edge_traversal)) || library_link->properties.emission_type == tfxLine && !(library_link->properties.flags & tfxEmitterPropertyFlags_edge_traversal) ? tfxEmitterStateFlags_not_line : 0;
		common.state_flags |= library_link->properties.flags & tfxEmitterPropertyFlags_random_color;
		common.state_flags |= library_link->properties.flags & tfxEmitterPropertyFlags_lifetime_uniform_size;
		common.state_flags |= library_link->properties.angle_setting != AngleSetting::tfxAlign && !(library_link->properties.flags & tfxEmitterPropertyFlags_relative_angle) ? tfxEmitterStateFlags_can_spin : 0;
		common.state_flags |= library_link->properties.angle_setting == AngleSetting::tfxAlign ? tfxEmitterStateFlags_align_with_velocity : 0;
		common.state_flags |= library_link->properties.emission_type == tfxLine && library_link->properties.flags & tfxEmitterPropertyFlags_edge_traversal ? tfxEmitterStateFlags_is_line_traversal : 0;
		common.state_flags |= library_link->properties.flags & tfxEmitterPropertyFlags_play_once;
		common.state_flags |= library_link->properties.end_behaviour == tfxLoop ? tfxEmitterStateFlags_loop : 0;
		common.state_flags |= library_link->properties.end_behaviour == tfxKill ? tfxEmitterStateFlags_kill : 0;
		common.state_flags |= library_link->properties.emission_type == tfxArea || library_link->properties.emission_type == tfxEllipse ? tfxEmitterStateFlags_is_area : 0;
		common.state_flags |= library_link->properties.emission_type == tfxLine ? tfxEmitterStateFlags_is_line : 0;
		common.state_flags |= single_done;
		common.state_flags |= spawn_done;
	}

	void tfxEffect::RefreshFromLibrary() {
		common.property_flags = library_link->properties.flags;
		common.state_flags = library_link->flags;
		common.loop_length = library_link->properties.loop_length;
		common.handle = library_link->properties.emitter_handle;
		for (auto &e : sub_emitters) {
			e.RefreshFromLibrary();
		}
		for (auto &e : sub_effects) {
			if (e.type == tfxEmitterType) {
				e.RefreshFromLibrary();
			}
			else {
				e.common.property_flags = e.library_link->properties.flags;
				e.common.state_flags = e.library_link->flags;
				e.common.loop_length = e.library_link->properties.loop_length;
				e.common.handle = e.library_link->properties.emitter_handle;
			}
		}
	}

	void tfxEffect::UpdateSpritePointers() {
		for (auto &emitter : sub_emitters) {
			for (auto &p : emitter.particles) {
				auto &s = emitter.common.root_effect->sprites[emitter.library_link->properties.layer][p.sprite_index];
				s.particle = &p;
			}
		}
		for (auto &emitter : sub_effects) {
			if (emitter.type == tfxEmitterType) {
				for (auto &p : emitter.particles) {
					auto &s = emitter.common.root_effect->sprites[emitter.library_link->properties.layer][p.sprite_index];
					s.particle = &p;
				}
			}
		}
	}
		
	bool tfxEmitter::GrowParticles(unsigned int min_amount) {
		auto &sprites = common.root_effect->sprites[library_link->properties.layer];
		unsigned int req_particle_size = particles.capacity + particles.capacity / 2;
		unsigned int req_sprite_size = sprites.capacity + sprites.capacity / 2;
		req_particle_size = req_particle_size < min_amount ? min_amount : req_particle_size;
		req_sprite_size = req_sprite_size < min_amount ? min_amount : req_sprite_size;
		unsigned int req_particle_mem = req_particle_size * sizeof(Particle);
		unsigned int req_sprite_mem = req_sprite_size * sizeof(ParticleSprite);
		if (!common.root_effect->storage->particle_memory.has_free_range_available(req_particle_mem) && common.root_effect->storage->particle_memory.free_unused_space() < req_particle_mem)
			return false;
		if (!common.root_effect->storage->sprite_memory.has_free_range_available(req_sprite_mem) && common.root_effect->storage->sprite_memory.free_unused_space() < req_sprite_mem)
			return false;

		tfxfixedvec<Particle> new_particle_memory;
		tfxfixedvec<ParticleSprite> new_sprite_memory;
		new_particle_memory.assign_memory(common.root_effect->storage->particle_memory, sizeof(Particle), req_particle_size);
		new_sprite_memory.assign_memory(common.root_effect->storage->sprite_memory, sizeof(ParticleSprite), req_sprite_size);

		if (particles.data && particles.capacity) {
			particles.copyto(common.root_effect->storage->particle_memory, new_particle_memory);
			particles.free_range(common.root_effect->storage->particle_memory);
		}
		particles = new_particle_memory;
		if (sprites.data && common.root_effect->sprites[library_link->properties.layer].capacity) {
			sprites.copyto(common.root_effect->storage->sprite_memory, new_sprite_memory);
			sprites.free_range(common.root_effect->storage->sprite_memory);
		}
		sprites = new_sprite_memory;
		common.root_effect->UpdateSpritePointers();
		return true;
	}

	void tfxEmitter::SpawnParticles() {
		if (common.state_flags & tfxEmitterStateFlags_single_shot_done || common.root_effect->common.state_flags & tfxEmitterStateFlags_stop_spawning ||
			(parent && parent->common.state_flags & tfxEmitterStateFlags_stop_spawning) ) {
			current.qty = 0.f;
			return;
		}

		float (*effect_lookup_callback)(Graph &graph, float age) = common.root_effect->lookup_mode == tfxPrecise ? LookupPrecise : LookupFast;
		if (!(common.property_flags & tfxEmitterPropertyFlags_single) && !(common.property_flags & tfxEmitterPropertyFlags_one_shot)) {
			current.qty = effect_lookup_callback(common.library->base_graphs[library_link->base].amount, common.frame);
			current.qty += random_generation.Range(effect_lookup_callback(common.library->variation_graphs[library_link->variation].amount, common.frame));

			if (common.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && (common.state_flags & tfxEmitterStateFlags_is_area)) {
				float area = current.emitter_size.x * current.emitter_size.y;
				current.qty = (current.qty / 10000.f) * area;
			}
			else if (common.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && common.state_flags & tfxEmitterStateFlags_is_line) {
				current.qty = (current.qty / 100.f) * current.emitter_size.y;
			}

			current.qty *= common.root_effect->current.amount;
			current.qty *= UPDATE_TIME;
			current.qty += current.amount_remainder;
		}
		else {
			current.qty = (float)library_link->properties.spawn_amount;
		}

		float tween = 0.f;
		float interpolate = (float)(int)current.qty;
		float count = 0;

		tfxEmitterSpawnControls spawn_values;
		if (current.qty >= 1) {
			spawn_values.life = effect_lookup_callback(common.library->base_graphs[library_link->base].life, common.frame) * common.root_effect->current.life;
			spawn_values.life_variation = effect_lookup_callback(common.library->variation_graphs[library_link->variation].life, common.frame) * common.root_effect->current.life;

			spawn_values.arc_size = 0.f;
			spawn_values.arc_offset = 0.f;
			if (library_link->properties.emission_type == EmissionType::tfxEllipse) {
				spawn_values.arc_size = effect_lookup_callback(common.library->property_graphs[library_link->property].arc_size, common.frame);
				spawn_values.arc_offset = effect_lookup_callback(common.library->property_graphs[library_link->property].arc_offset, common.frame);
			}
			spawn_values.weight = effect_lookup_callback(common.library->base_graphs[library_link->base].weight, common.frame) * common.root_effect->current.weight;
			spawn_values.weight_variation = effect_lookup_callback(common.library->variation_graphs[library_link->variation].weight, common.frame) * common.root_effect->current.weight;
			spawn_values.velocity = effect_lookup_callback(common.library->base_graphs[library_link->base].velocity, common.frame) * common.root_effect->current.velocity;
			spawn_values.velocity_variation = effect_lookup_callback(common.library->variation_graphs[library_link->variation].velocity, common.frame) * common.root_effect->current.velocity;
			if (!(common.property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
				spawn_values.size.x = effect_lookup_callback(common.library->base_graphs[library_link->base].width, common.frame) * common.root_effect->current.size.x;
				spawn_values.size.y = effect_lookup_callback(common.library->base_graphs[library_link->base].height, common.frame) * common.root_effect->current.size.y;
			}
			else {
				spawn_values.size.x = effect_lookup_callback(common.library->base_graphs[library_link->base].width, common.frame);
				if (common.root_effect->common.property_flags & tfxEmitterPropertyFlags_global_uniform_size)
					spawn_values.size.y = spawn_values.size.x * common.root_effect->current.size.x;
				else
					spawn_values.size.y = spawn_values.size.x * common.root_effect->current.size.y;
				spawn_values.size.x *= common.root_effect->current.size.x;
			}
			spawn_values.size_variation.x = effect_lookup_callback(common.library->variation_graphs[library_link->variation].width, common.frame) * common.root_effect->current.size.x;
			spawn_values.size_variation.y = effect_lookup_callback(common.library->variation_graphs[library_link->variation].height, common.frame) * common.root_effect->current.size.y;
			spawn_values.spin = effect_lookup_callback(common.library->base_graphs[library_link->base].spin, common.frame) * common.root_effect->current.spin;
			spawn_values.spin_variation = effect_lookup_callback(common.library->variation_graphs[library_link->variation].spin, common.frame) * common.root_effect->current.spin;
			spawn_values.splatter = effect_lookup_callback(common.library->property_graphs[library_link->property].splatter, common.frame) * common.root_effect->current.splatter;
			spawn_values.noise_offset_variation = effect_lookup_callback(common.library->variation_graphs[library_link->variation].noise_offset, common.frame);
			spawn_values.noise_offset = effect_lookup_callback(common.library->base_graphs[library_link->variation].noise_offset, common.frame);
			spawn_values.noise_resolution = effect_lookup_callback(common.library->variation_graphs[library_link->variation].noise_resolution, common.frame);

			if (common.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {
				if (library_link->properties.emission_type == EmissionType::tfxArea) {
					if (library_link->properties.grid_points.x > 1)
						spawn_values.grid_segment_size.x = current.emitter_size.x / (library_link->properties.grid_points.x - 1);
					if (library_link->properties.grid_points.y > 1)
						spawn_values.grid_segment_size.y = current.emitter_size.y / (library_link->properties.grid_points.y - 1);
				}
				else if (library_link->properties.emission_type == EmissionType::tfxEllipse) {
					if (library_link->properties.grid_points.x > 0)
						spawn_values.grid_segment_size.x = spawn_values.arc_size / (library_link->properties.grid_points.x);
				}
				else if (library_link->properties.emission_type == EmissionType::tfxLine) {
					if (library_link->properties.grid_points.x > 1)
						spawn_values.grid_segment_size.y = current.emitter_size.y / (library_link->properties.grid_points.x - 1);
				}
			}

			if (library_link->spawn_update_callback)
				library_link->spawn_update_callback(spawn_values, *this);
		}

		while (current.qty >= 1.f) {
			if (!FreeCapacity()) {
				if (common.root_effect->common.property_flags & tfxEmitterPropertyFlags_can_grow_particle_memory) {
					if (!GrowParticles(particles.current_size + (unsigned int)current.qty)) {
						current.amount_remainder = 0;
						break;
					}
				}
				else {
					current.amount_remainder = 0;
					break;
				}
			}
			else if (common.root_effect->common.property_flags & tfxEmitterPropertyFlags_can_grow_particle_memory) {
			}
			tween = count / interpolate;
			count++;
			current.qty -= 1.f;
			common.timeout_counter = 0;

			bool is_single = common.property_flags & tfxEmitterPropertyFlags_single && !(common.property_flags & tfxEmitterPropertyFlags_one_shot);

			Particle &p = GrabParticle();
			InitCPUParticle(p, spawn_values, tween);

			ParticleSprite &s = common.root_effect->GrabSprite(library_link->properties.layer);
			s.parameters = (unsigned int)p.image_frame;
			s.color = p.color;
			s.world = p.world;
			s.captured = p.captured;
			s.ptr = library_link->properties.image->ptr;
			s.intensity = p.intensity;
			s.handle = image_handle;
			s.particle = &p;
			p.sprite_index = common.root_effect->sprites[library_link->properties.layer].last_index();

			common.root_effect->common.highest_particle_age = std::fmaxf(common.root_effect->common.highest_particle_age, p.max_age);
			common.highest_particle_age = std::fmaxf(common.highest_particle_age, p.max_age);
		}

		current.amount_remainder = current.qty;
	}

	void EffectEmitter::NoTweenNextUpdate() {
		tfxvec<EffectEmitter*> stack;
		stack.push_back(this);
		while (!stack.empty()) {
			auto &current = stack.pop_back();
			current->flags |= tfxEmitterStateFlags_no_tween_this_update;
			for (auto &sub : current->sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}


	float tfxEmitter::GetEmissionDirection(tfxVec2 &local_position, tfxVec2 &world_position, tfxVec2 &emitter_size) {
		float (*effect_lookup_callback)(Graph &graph, float age) = common.root_effect->lookup_mode == tfxPrecise ? LookupPrecise : LookupFast;
		float emission_angle = effect_lookup_callback(common.library->property_graphs[library_link->property].emission_angle, common.frame);
		float emission_angle_variation = effect_lookup_callback(common.library->property_graphs[library_link->property].emission_range, common.frame);
		//----Emission
		float range = emission_angle_variation *.5f;
		float direction = 0;

		if (library_link->properties.emission_type == EmissionType::tfxPoint)
			return direction + emission_angle + random_generation.Range(-range, range);

		tfxVec2 tmp_position;
		if (local_position.x == 0 && local_position.y == 0)
			tmp_position = emitter_size;
		else
			tmp_position = local_position;

		if (library_link->properties.emission_direction == EmissionDirection::tfxOutwards) {

			tfxVec2 to_handle;

			if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = (tmp_position);
			else
				to_handle = (world_position - transform.world.position);

			direction = GetVectorAngle(to_handle.x, to_handle.y);

		}
		else if (library_link->properties.emission_direction == EmissionDirection::tfxInwards) {

			tfxVec2 to_handle;

			if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = (-tmp_position);
			else
				to_handle = (transform.world.position - world_position);

			direction = GetVectorAngle(to_handle.x, to_handle.y);

		}
		else if (library_link->properties.emission_direction == EmissionDirection::tfxBothways) {

			if (current.emission_alternator) {

				tfxVec2 to_handle;

				if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = (tmp_position);
				else
					to_handle = (world_position - transform.world.position);

				direction = GetVectorAngle(to_handle.x, to_handle.y);

			}
			else {

				tfxVec2 to_handle;

				if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = (-tmp_position);
				else
					to_handle = (transform.world.position - world_position);

				direction = GetVectorAngle(to_handle.x, to_handle.y);

			}

			current.emission_alternator = !current.emission_alternator;
		}

		if (std::isnan(direction))
			direction = 0.f;
		return direction + emission_angle + random_generation.Range(-range, range);
	}

	Particle& tfxEmitter::GrabParticle() {
		//Must check for free capacity before calling this function. Internal use only
		assert(particles.current_size != particles.capacity);
		return particles[particles.current_size++];
	}

	ParticleSprite& tfxEffect::GrabSprite(unsigned int layer) {
		//Must check for free capacity before calling this function. Internal use only
		assert(sprites[layer].current_size != sprites[layer].capacity);
		return sprites[layer][sprites[layer].current_size++];
	}

	tfxEmitter& tfxEffect::GrabSubEffect() {
		assert(sub_effects.current_size != sub_effects.capacity);
		return sub_effects[sub_effects.current_size++];
	}

	tfxEmitter& tfxEffect::AddSubEffect(tfxEmitter &sub_effect) {
		tfxEmitter &e = GrabSubEffect();
		e = sub_effect;
		e.next_emitter = &e;
		return e;
	}

	bool tfxEmitter::FreeCapacity() {
		return !particles.full() && !common.root_effect->sprites[library_link->properties.layer].full();
	}

	void *tfxEmitter::UserData() { 
		return common.root_effect->user_data; 
	}

	void tfxEmitter::InitCPUParticle(Particle &p, tfxEmitterSpawnControls &spawn_values, float tween) {
		p.flags = tfxParticleFlags_fresh;
		p.next_ptr = &p;

		if (common.property_flags & (tfxEmitterPropertyFlags_single | tfxEmitterPropertyFlags_one_shot))
			common.state_flags |= tfxEmitterStateFlags_single_shot_done;

		//----Life
		p.max_age = spawn_values.life + random_generation.Range(spawn_values.life_variation);
		p.age = 0.f;

		//----Position
		if (library_link->properties.emission_type == EmissionType::tfxPoint) {
			if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
				p.local.position = -common.handle;
			else {
				if (common.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) {
					p.local.position = InterpolateVec2(tween, transform.captured.position, transform.world.position);
				}
				else {
					tfxVec2 rotvec = transform.matrix.TransformVector(-common.handle);
					tfxVec2 spawn_position = InterpolateVec2(tween, transform.captured.position, transform.world.position) * transform.world.scale;
					p.local.position = rotvec + spawn_position;
				}
			}
		}
		else if (library_link->properties.emission_type == EmissionType::tfxArea) {
			tfxVec2 position = tfxVec2(0.f, 0.f);

			if (common.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				if (common.property_flags & tfxEmitterPropertyFlags_fill_area) {
					if (!(common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
						current.grid_coords.x--;
						if (current.grid_coords.x < 0.f) {
							current.grid_coords.y--;
							current.grid_coords.x = library_link->properties.grid_points.x - 1;
							if (current.grid_coords.y < 0.f)
								current.grid_coords.y = library_link->properties.grid_points.y - 1;
						}
					}

					p.local.position = position + (current.grid_coords * spawn_values.grid_segment_size) + common.handle;

					if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
						current.grid_coords.x++;
						if (current.grid_coords.x == library_link->properties.grid_points.x) {
							current.grid_coords.y++;
							current.grid_coords.x = 0.f;
							if (current.grid_coords.y >= library_link->properties.grid_points.y)
								current.grid_coords.y = 0.f;
						}
					}
				}
				else {

					if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {

						current.grid_direction.x = 1;
						current.grid_direction.y = 0;
						if (current.grid_coords.x == library_link->properties.grid_points.x - 1 && current.grid_coords.y >= 0 && current.grid_coords.y < library_link->properties.grid_points.y - 1) {
							current.grid_direction.x = 0;
							current.grid_direction.y = 1;
						}
						else if (current.grid_coords.x > 0 && current.grid_coords.x < library_link->properties.grid_points.x && current.grid_coords.y == library_link->properties.grid_points.y - 1) {
							current.grid_direction.x = -1;
							current.grid_direction.y = 0;
						}
						else if (current.grid_coords.x == 0 && current.grid_coords.y > 0 && current.grid_coords.y < library_link->properties.grid_points.y) {
							current.grid_direction.x = 0;
							current.grid_direction.y = -1;
						}

					}
					else {

						current.grid_direction.x = -1;
						current.grid_direction.y = 0;
						if (current.grid_coords.x == library_link->properties.grid_points.x - 1 && current.grid_coords.y > 0 && current.grid_coords.y < library_link->properties.grid_points.y) {
							current.grid_direction.x = 0;
							current.grid_direction.y = -1;
						}
						else if (current.grid_coords.x >= 0 && current.grid_coords.x < library_link->properties.grid_points.x - 1 && current.grid_coords.y == library_link->properties.grid_points.y - 1) {
							current.grid_direction.x = 1;
							current.grid_direction.y = 0;
						}
						else if (current.grid_coords.x == 0 && current.grid_coords.y >= 0 && current.grid_coords.y < library_link->properties.grid_points.y - 1) {
							current.grid_direction.x = 0;
							current.grid_direction.y = 1;
						}

					}

					current.grid_coords += current.grid_direction;
					tfxBound(current.grid_coords, library_link->properties.grid_points);
					p.local.position = position + (current.grid_coords * spawn_values.grid_segment_size) + common.handle;
				}
			}
			else {
				if (common.property_flags & tfxEmitterPropertyFlags_fill_area) {
					position.x = random_generation.Range(current.emitter_size.x);
					position.y = random_generation.Range(current.emitter_size.y);
				}
				else {
					//Spawn on one of 4 edges of the area
					unsigned int side = random_generation.RangeUInt(4);
					if (side == 0) {
						//left side
						position.x = 0.f;
						position.y = random_generation.Range(current.emitter_size.y);
					}
					else if (side == 1) {
						//right side
						position.x = current.emitter_size.x;
						position.y = random_generation.Range(current.emitter_size.y);
					}
					else if (side == 2) {
						//top side
						position.x = random_generation.Range(current.emitter_size.x);
						position.y = 0.f;
					}
					else if (side == 3) {
						//bottom side
						position.x = random_generation.Range(current.emitter_size.x);
						position.y = current.emitter_size.y;
					}
				}

				p.local.position = position + common.handle;
			}

			//----TForm and Emission
			if (!(common.property_flags & tfxEmitterPropertyFlags_relative_position)) {
				p.local.position = transform.matrix.TransformVector(tfxVec2(p.local.position.x, p.local.position.y));
				p.local.position = transform.world.position + p.local.position * transform.world.scale;
			}

		}
		else if (library_link->properties.emission_type == EmissionType::tfxEllipse) {
			tfxVec2 emitter_size = (current.emitter_size * .5f);
			tfxVec2 position = tfxVec2(0.f, 0.f);

			if (common.property_flags & tfxEmitterPropertyFlags_spawn_on_grid && !(common.property_flags & tfxEmitterPropertyFlags_fill_area)) {

				current.grid_coords.y = 0.f;

				if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					current.grid_coords.x--;
					if (current.grid_coords.x < 0.f) {
						current.grid_coords.x = library_link->properties.grid_points.x - 1;
					}
				}

				float th = current.grid_coords.x * spawn_values.grid_segment_size.x + spawn_values.arc_offset;
				p.local.position = tfxVec2(std::cosf(th) * emitter_size.x + common.handle.x + emitter_size.x,
					-std::sinf(th) * emitter_size.y + common.handle.y + emitter_size.y);

				if (!(common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					current.grid_coords.x++;
					if (current.grid_coords.x >= library_link->properties.grid_points.x) {
						current.grid_coords.x = 0.f;
					}
				}

			}
			else if (!(common.property_flags & tfxEmitterPropertyFlags_fill_area)) {
				float th = random_generation.Range(spawn_values.arc_size) + spawn_values.arc_offset;

				p.local.position = tfxVec2(std::cosf(th) * emitter_size.x + common.handle.x + emitter_size.x,
					-std::sinf(th) * emitter_size.y + common.handle.y + emitter_size.y);

			}
			else {
				p.local.position.x = random_generation.Range(-emitter_size.x, emitter_size.x);
				p.local.position.y = random_generation.Range(-emitter_size.y, emitter_size.y);

				while ((std::pow(p.local.position.x, 2) / std::pow(emitter_size.x, 2)) + (std::pow(p.local.position.y, 2) / std::pow(emitter_size.y, 2)) > 1) {
					p.local.position.x = random_generation.Range(-emitter_size.x, emitter_size.x);
					p.local.position.y = random_generation.Range(-emitter_size.y, emitter_size.y);
				}
			}

			//----TForm and Emission
			if (!(common.property_flags & tfxEmitterPropertyFlags_relative_position)) {
				p.local.position = transform.matrix.TransformVector(tfxVec2(p.local.position.x, p.local.position.y));
				p.local.position = transform.world.position + p.local.position * transform.world.scale;
			}

		}
		else if (library_link->properties.emission_type == EmissionType::tfxLine) {
			if (common.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				current.grid_coords.x = 0.f;

				if (!(common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					current.grid_coords.y--;
					if (current.grid_coords.y < 0.f) {
						current.grid_coords.y = library_link->properties.grid_points.x - 1;
					}
				}

				p.local.position = tfxVec2(current.grid_coords * -spawn_values.grid_segment_size);
				p.local.position += common.handle;

				if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					current.grid_coords.y++;
					if (current.grid_coords.y >= library_link->properties.grid_points.x) {
						current.grid_coords.y = 0.f;
					}
				}

			}
			else {
				p.local.position.x = 0.f;
				p.local.position.y = random_generation.Range(-current.emitter_size.y, 0.f);

				p.local.position += common.handle;

			}

			//----TForm and Emission
			if (!(common.property_flags & tfxEmitterPropertyFlags_relative_position) && !(common.property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
				p.local.position = transform.matrix.TransformVector(tfxVec2(p.local.position.x, p.local.position.y));
				p.local.position = transform.world.position + p.local.position * transform.world.scale;
			}

		}

		//----Weight
		if (spawn_values.weight) {
			p.base.weight = spawn_values.weight * common.library->overtime_graphs[library_link->overtime].weight.GetFirstValue();
			if (spawn_values.weight_variation > 0) {
				p.base.weight += random_generation.Range(-spawn_values.weight_variation, spawn_values.weight_variation) * common.library->overtime_graphs[library_link->overtime].weight.GetFirstValue();
			}
		}
		else {
			p.base.weight = 0;
		}
		p.weight_acceleration = p.base.weight * common.library->overtime_graphs[library_link->overtime].weight.GetFirstValue() * UPDATE_TIME;

		//----Velocity
		float velocity_scale = common.library->overtime_graphs[library_link->overtime].velocity.GetFirstValue() * current.velocity_adjuster;
		p.base.velocity = spawn_values.velocity + random_generation.Range(-spawn_values.velocity_variation, spawn_values.velocity_variation);

		//----Size
		if (!(common.property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			p.base.random_size.x = random_generation.Range(spawn_values.size_variation.x);
			p.base.random_size.y = random_generation.Range(spawn_values.size_variation.y);
			p.base.size.y = p.base.random_size.y + spawn_values.size.y;
			p.base.size.x = (p.base.random_size.x + spawn_values.size.x) / library_link->properties.image->image_size.x;
			float height = p.base.size.y / library_link->properties.image->image_size.y;

			p.world.scale.x = p.base.size.x * common.library->overtime_graphs[library_link->overtime].width.GetFirstValue();

			if (common.library->overtime_graphs[library_link->overtime].stretch.GetFirstValue()) {
				float velocity = std::fabsf(velocity_scale * p.base.velocity) * UPDATE_TIME;
				velocity += p.weight_acceleration * UPDATE_TIME;
				p.world.scale.y = (common.library->overtime_graphs[library_link->overtime].height.GetFirstValue() * common.root_effect->current.size.y * (p.base.size.y + (velocity * common.library->overtime_graphs[library_link->overtime].stretch.GetFirstValue() * common.root_effect->current.stretch))) / library_link->properties.image->image_size.y;
			}
			else {
				if (common.property_flags & tfxEmitterPropertyFlags_lifetime_uniform_size) {
					p.world.scale.y = p.world.scale.x;
				}
				else {
					p.world.scale.y = height * common.library->overtime_graphs[library_link->overtime].height.GetFirstValue();
				}
			}
		}
		else {
			p.base.random_size.x = random_generation.Range(spawn_values.size_variation.x);
			p.base.random_size.y = p.base.random_size.x;
			p.base.size.y = p.base.random_size.y + spawn_values.size.y;
			p.base.size.x = (p.base.random_size.x + spawn_values.size.x) / library_link->properties.image->image_size.x;
			float height = p.base.size.y / library_link->properties.image->image_size.y;

			p.world.scale.x = p.base.size.x * common.library->overtime_graphs[library_link->overtime].width.GetFirstValue();

			if (common.library->overtime_graphs[library_link->overtime].stretch.GetFirstValue()) {
				float velocity = std::fabsf(velocity_scale * p.base.velocity) * UPDATE_TIME;
				velocity += p.weight_acceleration * UPDATE_TIME;
				p.world.scale.y = (common.library->overtime_graphs[library_link->overtime].width.GetFirstValue() * common.root_effect->current.size.y * (p.base.size.y + (velocity * common.library->overtime_graphs[library_link->overtime].stretch.GetFirstValue() * common.root_effect->current.stretch))) / library_link->properties.image->image_size.y;
			}
			else {
				p.world.scale.y = p.world.scale.x;
			}
		}

		//----Spin
		p.base.spin = random_generation.Range(-spawn_values.spin_variation, std::abs(spawn_values.spin_variation)) + spawn_values.spin;

		switch (library_link->properties.angle_setting) {
		case AngleSetting::tfxRandom:
			p.captured.rotation = p.local.rotation = random_generation.Range(library_link->properties.angle_offset);
			break;
		case AngleSetting::tfxSpecify:
			p.captured.rotation = p.local.rotation = library_link->properties.angle_offset;
			break;
		default:
			p.captured.rotation = p.local.rotation = 0;
			break;
		}

		//----Splatter
		if (spawn_values.splatter) {
			float splattertemp = spawn_values.splatter;
			float splatx = random_generation.Range(-spawn_values.splatter, spawn_values.splatter);
			float splaty = random_generation.Range(-spawn_values.splatter, spawn_values.splatter);

			while (GetDistance(0, 0, splatx, splaty) >= splattertemp && splattertemp > 0) {
				splatx = random_generation.Range(-spawn_values.splatter, spawn_values.splatter);
				splaty = random_generation.Range(-spawn_values.splatter, spawn_values.splatter);
			}

			if (!(common.property_flags & tfxEmitterPropertyFlags_relative_position)) {
				p.local.position.x += splatx * transform.world.scale.x;
				p.local.position.y += splaty * transform.world.scale.y;
			}
			else {
				p.local.position.x += splatx;
				p.local.position.y += splaty;
			}
		}

		float direction = 0;

		if (library_link->properties.angle_setting == AngleSetting::tfxAlign && common.property_flags & tfxEmitterPropertyFlags_edge_traversal)
			p.world.rotation = p.local.rotation = direction + library_link->properties.angle_offset;

		bool line = common.property_flags & tfxEmitterPropertyFlags_edge_traversal && library_link->properties.emission_type == EmissionType::tfxLine;

		TransformParticle(p, *this);
		p.captured = p.world;
		p.captured.scale = p.world.scale;

		//----Motion randomness
		p.noise_offset = random_generation.Range(spawn_values.noise_offset_variation) + spawn_values.noise_offset;
		p.noise_resolution = spawn_values.noise_resolution + 0.01f;

		if (!(common.property_flags & tfxEmitterPropertyFlags_edge_traversal) || library_link->properties.emission_type != EmissionType::tfxLine) {
			direction = p.emission_angle = GetEmissionDirection(p.local.position, p.world.position, current.emitter_size) + common.library->overtime_graphs[library_link->overtime].direction.GetFirstValue();
		}

		//----Normalize Velocity to direction
		tfxVec2 velocity_normal;
		velocity_normal.x = std::sinf(direction);
		velocity_normal.y = -std::cosf(direction);

		//p.velocity = p.velocity_normal * p.base.velocity * p.velocity_scale * UPDATE_TIME;

		if ((library_link->properties.angle_setting == AngleSetting::tfxAlign || library_link->properties.angle_setting == tfxAlignWithEmission) && !line) {
			p.world.rotation = p.local.rotation = GetVectorAngle(velocity_normal.x, velocity_normal.y) + library_link->properties.angle_offset;
			if (common.property_flags & tfxEmitterPropertyFlags_relative_angle)
				p.world.rotation += transform.world.rotation;
			p.captured.rotation = p.world.rotation;
			//Reset the matrix again so that any child particles spawn in the correct place
			if (library_link->sub_effectors.size()) {
				float s = sin(p.local.rotation);
				float c = cos(p.local.rotation);
				p.matrix.Set(c, s, -s, c);
			}
		}

		//----Handle
		/*if (common.property_flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
			p.handle = tfxVec2(0.5f, 0.5f);
		}
		else {
			p.handle = properties.image_handle;
		}*/

		//----Image
		//p.image = properties.image;
		if (common.property_flags & tfxEmitterPropertyFlags_random_start_frame && library_link->properties.image->animation_frames > 1) {
			p.image_frame = random_generation.Range(library_link->properties.image->animation_frames);
		}
		else {
			p.image_frame = library_link->properties.start_frame;
		}

		//----Color
		p.color.a = unsigned char(255.f * common.library->overtime_graphs[library_link->overtime].blendfactor.GetFirstValue());
		p.intensity = common.library->overtime_graphs[library_link->overtime].intensity.GetFirstValue() * common.root_effect->current.intensity;
		if (common.property_flags & tfxEmitterPropertyFlags_random_color) {
			float age = random_generation.Range(p.max_age);
			p.color.r = unsigned char(255.f * lookup_overtime_callback(common.library->overtime_graphs[library_link->overtime].red, age, p.max_age));
			p.color.g = unsigned char(255.f * lookup_overtime_callback(common.library->overtime_graphs[library_link->overtime].green, age, p.max_age));
			p.color.b = unsigned char(255.f * lookup_overtime_callback(common.library->overtime_graphs[library_link->overtime].blue, age, p.max_age));
		}
		else {
			p.color.r = unsigned char(255.f * common.library->overtime_graphs[library_link->overtime].red.GetFirstValue());
			p.color.g = unsigned char(255.f * common.library->overtime_graphs[library_link->overtime].green.GetFirstValue());
			p.color.b = unsigned char(255.f * common.library->overtime_graphs[library_link->overtime].blue.GetFirstValue());
		}

		for(auto &sub_effect : library_link->sub_effectors) {
			tfxEmitter &new_sub_effect = common.root_effect->GrabSubEffect();
			new_sub_effect.type = tfxEffectType;
			new_sub_effect.common.root_effect = common.root_effect;
			new_sub_effect.common.highest_particle_age = 0;
			new_sub_effect.common.timeout = 50.f;
			new_sub_effect.common.timeout_counter = 0.f;
			new_sub_effect.common.active_children = 0;
			new_sub_effect.common.handle = sub_effect.properties.emitter_handle;
			new_sub_effect.common.state_flags = sub_effect.flags;
			new_sub_effect.common.age = 0.f;
			new_sub_effect.transform.local.position = tfxVec2();
			new_sub_effect.transform.local.rotation = 0.f;
			new_sub_effect.next_emitter = &new_sub_effect;
			new_sub_effect.library_link = &sub_effect;
			new_sub_effect.offset = 0;
			for (auto &sub_emitter : sub_effect.sub_effectors) {
				tfxEmitter &new_emitter = common.root_effect->GrabSubEffect();
				new_emitter.type = tfxEmitterType;
				new_emitter.particles.reset_to_null();
				new_emitter.common.root_effect = common.root_effect;
				new_emitter.next_emitter = &new_emitter;
				new_emitter.common.highest_particle_age = 0;
				new_emitter.common.timeout = 50.f;
				new_emitter.common.timeout_counter = 0.f;
				new_emitter.transform.local.position = tfxVec2();
				new_emitter.transform.local.rotation = 0.f;
				new_emitter.current.amount_remainder = 0.f;
				new_emitter.current.emission_alternator = 0;
				new_emitter.parent = &new_sub_effect;
				new_emitter.offset = 0;
				new_sub_effect.common.active_children++;
				sub_emitter.CopyToEmitter(new_emitter, *common.root_effect->storage, true);
			}	
			new_sub_effect.parent_particle = &p;
		}

		if (library_link->particle_onspawn_callback)
			library_link->particle_onspawn_callback(p);
	}

	void ReloadBaseValues(Particle &p, EffectEmitter &e) {
		//----Life
		//std::uniform_real_distribution<float> random_life(0, e.current.life_variation);
		//p.max_age = e.current.life + random_life(random_generation.engine);

		//----Velocity
		//std::uniform_real_distribution<float> random_velocity;
		//if (e.current.velocity_variation > 0)
			//random_velocity = std::uniform_real_distribution<float>(0, e.current.velocity_variation);
		//else
			//random_velocity = std::uniform_real_distribution<float>(e.current.velocity_variation, 0);
		//p.velocity_scale = e.library->overtime_graphs[e.overtime].velocity.GetFirstValue() * e.current.velocity_adjuster;
		//p.base.velocity = e.current.velocity + random_velocity(random_generation.engine);
		float velocity_scale = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].velocity, p.age, p.max_age) * e.current.velocity_adjuster;

		//----Size
		if (!(e.properties.flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			//std::uniform_real_distribution<float> random_width(0, e.current.size_variation.x);
			//std::uniform_real_distribution<float> random_height(0, e.current.size_variation.y);

			p.base.size.y = p.base.random_size.y + e.current.size.y;
			p.base.size.x = (p.base.random_size.x + e.current.size.x) / e.properties.image->image_size.x;
			//p.base.size.y = p.base.height / e.properties.image->image_size.y;

			//p.local.scale.x = p.base.size.x * e.library->overtime_graphs[e.overtime].width.GetFirstValue();

			if (e.library->overtime_graphs[e.overtime].stretch.GetFirstValue()) {
				float velocity = std::fabsf(velocity_scale * p.base.velocity);
				//p.local.scale.y = (e.library->overtime_graphs[e.overtime].height.GetFirstValue() * e.parent->current.size.y * (p.base.height + (velocity * e.library->overtime_graphs[e.overtime].stretch.GetFirstValue() * e.parent->current.stretch))) / e.properties.image.image_size.y;
			}
			else {
				if (e.properties.flags & tfxEmitterPropertyFlags_lifetime_uniform_size) {
					//p.local.scale.y = p.base.size.y * e.library->overtime_graphs[e.overtime].width.GetFirstValue();
				}
				else {
					//p.local.scale.y = p.base.size.y * e.library->overtime_graphs[e.overtime].height.GetFirstValue();
				}
			}
		}
		else {
			//std::uniform_real_distribution<float> random_width(0, e.current.size_variation.x);

			p.base.size.y = p.base.random_size.y + e.current.size.x;
			p.base.size.x = (p.base.random_size.x + e.current.size.x) / e.properties.image->image_size.x;
			//p.base.size.y = p.base.height / e.properties.image->image_size.y;

			//p.local.scale.x = p.base.size.x * e.library->overtime_graphs[e.overtime].width.GetFirstValue();

			if (e.library->overtime_graphs[e.overtime].stretch.GetFirstValue()) {
				float velocity = std::fabsf(velocity_scale * p.base.velocity);
				//p.local.scale.y = (e.library->overtime_graphs[e.overtime].height.GetFirstValue() * e.parent->current.size.y * (p.base.height + (velocity * e.library->overtime_graphs[e.overtime].stretch.GetFirstValue() * e.parent->current.stretch))) / e.properties.image.image_size.y;
			}
			else {
				//p.local.scale.y = p.local.scale.x;
			}
		}

		//----Spin
		p.base.spin = random_generation.Range(-e.current.spin_variation, std::abs(e.current.spin_variation)) + e.current.spin;

		//std::uniform_real_distribution<float> random_angle(0, e.properties.angle_offset);
		//switch (e.properties.angle_setting) {
		//case AngleSetting::kRandom:
			//p.local.rotation = random_angle(random_generation.engine);
			//break;
		//case AngleSetting::kSpecify:
			//p.local.rotation = e.properties.angle_offset;
			//break;
		//default:
			//break;
		//}

		//----Weight
		if (e.current.weight) {
			if (e.current.weight_variation > 0) {
				p.base.weight = (random_generation.Range(e.current.weight_variation) + e.current.weight) * e.library->overtime_graphs[e.overtime].weight.GetFirstValue();
			}
			else {
				p.base.weight = (random_generation.Range(e.current.weight_variation, 0) + e.current.weight) * e.library->overtime_graphs[e.overtime].weight.GetFirstValue();
			}
			//p.weight_acceleration = 0;
		}
		else {
			//p.weight_acceleration = 0;
			p.base.weight = 0;
		}

		//TransformParticle(p, e);

		//----Motion randomness
		//p.noise_offset = e.current.noise_offset;
		//float mr = tfxRadians(p.direction_turbulance * e.library->overtime_graphs[e.overtime].direction_turbulance.GetFirstValue());
		//std::uniform_real_distribution<float> random_motion(-mr, mr);
		//p.motion_randomness_direction = tfxRadians(22.5f * random_motion(random_generation.engine));
		//p.motion_randomness_speed = 30.f * random_motion(random_generation.engine);
		//p.motion_tracker = 0;

		//if (!e.properties.edge_traversal || e.properties.emission_type != EmissionType::kLine) {
			//p.direction = p.emission_angle = e.GetEmissionDirection(p) + p.motion_randomness_direction;
		//}

		//----Normalize Velocity to direction
		//p.velocity_normal.x = std::sinf(p.direction);
		//p.velocity_normal.y = -std::cosf(p.direction);

		//p.velocity = p.velocity_normal * p.base.velocity * p.velocity_scale * e.timer->UpdateTime();
		//bool line = e.properties.edge_traversal && e.properties.emission_type == EmissionType::kLine;

		//if (e.properties.angle_setting == AngleSetting::kAlign && !line) {
			//p.world.rotation = p.local.rotation = GetVectorAngle(p.velocity_normal.x, p.velocity_normal.y) + e.properties.angle_offset + e.world.rotation;
			//p.captured.rotation = p.world.rotation;
		//}

		//----Handle
		//if (e.properties.flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
			//p.handle = tfxVec2(0.5f, 0.5f);
		//}
		//else {
			//p.handle = e.properties.image_handle;
		//}

		//----Image
		//p.image = properties.image;
		//if (e.properties.image.random_start_frame && e.properties.image.animation_frames > 1) {
			//std::uniform_real_distribution<float> random_start_frame(0.f, e.properties.image.animation_frames - 1);
			//p.image_frame = random_start_frame(random_generation.engine);
		//}
		//else {
			//p.image_frame = e.properties.image.start_frame;
		//}

		//----Color
		//p.color.a = unsigned char(255.f * e.library->overtime_graphs[e.overtime].opacity.GetFirstValue());
		//p.intensity = e.library->overtime_graphs[e.overtime].opacity.GetFirstValue();
		//if (e.properties.random_color) {
			//std::uniform_real_distribution<float> random_color(0.f, p.max_age);
			//float age = random_color(random_generation.engine);
			//p.color.r = unsigned char(255.f * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].red, age, p.max_age));
			//p.color.g = unsigned char(255.f * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].green, age, p.max_age));
			//p.color.b = unsigned char(255.f * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].blue, age, p.max_age));
		//}
		//else {
			//p.color.r = unsigned char(255.f * e.library->overtime_graphs[e.overtime].red.GetFirstValue());
			//p.color.g = unsigned char(255.f * e.library->overtime_graphs[e.overtime].green.GetFirstValue());
			//p.color.b = unsigned char(255.f * e.library->overtime_graphs[e.overtime].blue.GetFirstValue());
		//}

	}

	void tfxEmitter::ControlParticles() {
		bool can_bump = true;

		unsigned int index_offset = 0;

		tfxEffect &root_effect = *common.root_effect;
		unsigned int flags = common.state_flags;
		float velocity_adjuster = current.velocity_adjuster;
		float global_intensity = root_effect.current.intensity;
		float image_size_y = library_link->properties.image->image_size.y;
		float image_frame_rate = library_link->properties.image->animation_frames > 1 && common.property_flags & tfxEmitterPropertyFlags_animate ? library_link->properties.frame_rate : 0.f;
		float stretch = current.stretch;
		float emitter_size_y = current.emitter_size.y;
		float emitter_handle_y = common.handle.y;
		float overal_scale = current.overal_scale;
		float angle_offset = library_link->properties.angle_offset;
		OvertimeAttributes &graphs = common.library->overtime_graphs[library_link->overtime];

		unsigned int offset_index = 0;

		for(unsigned int i = 0 ; i != particles.current_size; ++i) {
			auto &p = particles[i];

			if (!(p.flags & tfxParticleFlags_fresh)) {
				p.captured = p.world;
			}

			p.age += FRAME_LENGTH;

			if (common.state_flags & tfxEmitterStateFlags_remove) {
				index_offset++;
				ParticleSprite &s = root_effect.sprites[library_link->properties.layer][p.sprite_index];
				s.parameters = tfxINVALID;
				p.next_ptr = nullptr;
				continue;
			}
			if (!(common.state_flags & tfxEmitterStateFlags_stop_spawning) && p.age >= p.max_age && common.state_flags & tfxEmitterStateFlags_is_single) {
				p.age = 0;
			}
			else if (p.age >= p.max_age) {
				index_offset++;
				ParticleSprite &s = root_effect.sprites[library_link->properties.layer][p.sprite_index];
				s.parameters = tfxINVALID;
				p.next_ptr = nullptr;
				continue;
			}

			//-------------------------------------------------------
			//Control what the particle does over the course of
			//it's lifetime
			//-------------------------------------------------------

			u32 lookup_frame = static_cast<u32>((p.age / p.max_age * graphs.velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

			float lookup_velocity = graphs.velocity.lookup.values[std::min<u32>(lookup_frame, graphs.velocity.lookup.last_frame)] * velocity_adjuster;
			float lookup_velocity_turbulance = graphs.velocity_turbulance.lookup.values[std::min<u32>(lookup_frame, graphs.velocity_turbulance.lookup.last_frame)];
			float lookup_direction_turbulance = graphs.direction_turbulance.lookup.values[std::min<u32>(lookup_frame, graphs.direction_turbulance.lookup.last_frame)];
			float lookup_direction = graphs.direction.lookup.values[std::min<u32>(lookup_frame, graphs.direction.lookup.last_frame)] + p.emission_angle;
			float lookup_noise_resolution = graphs.noise_resolution.lookup.values[std::min<u32>(lookup_frame, graphs.noise_resolution.lookup.last_frame)] * p.noise_resolution;
			float lookup_weight = graphs.weight.lookup.values[std::min<u32>(lookup_frame, graphs.weight.lookup.last_frame)];

			float lookup_width = graphs.width.lookup.values[std::min<u32>(lookup_frame, graphs.width.lookup.last_frame)];
			float lookup_height = graphs.height.lookup.values[std::min<u32>(lookup_frame, graphs.height.lookup.last_frame)];
			float lookup_stretch = graphs.stretch.lookup.values[std::min<u32>(lookup_frame, graphs.stretch.lookup.last_frame)];

			float lookup_spin = graphs.spin.lookup.values[std::min<u32>(lookup_frame, graphs.spin.lookup.last_frame)] * p.base.spin;

			float lookup_red = graphs.red.lookup.values[std::min<u32>(lookup_frame, graphs.red.lookup.last_frame)];
			float lookup_green = graphs.green.lookup.values[std::min<u32>(lookup_frame, graphs.green.lookup.last_frame)];
			float lookup_blue = graphs.blue.lookup.values[std::min<u32>(lookup_frame, graphs.blue.lookup.last_frame)];
			float lookup_opacity = graphs.blendfactor.lookup.values[std::min<u32>(lookup_frame, graphs.blendfactor.lookup.last_frame)];
			float lookup_intensity = graphs.intensity.lookup.values[std::min<u32>(lookup_frame, graphs.intensity.lookup.last_frame)];

			float direction = 0;
			float mr_angle = 0;
			float mr_speed = 0;

			tfxVec2 mr_vec;
			if (flags & tfxEmitterStateFlags_not_line) {
				direction = lookup_direction;
			}

			if (lookup_velocity_turbulance + lookup_direction_turbulance) {
				float noise = SimplexNoise::noise(p.local.position.x / lookup_noise_resolution + p.noise_offset, p.local.position.y / lookup_noise_resolution + p.noise_offset);
				mr_speed = noise * lookup_velocity_turbulance;
				mr_angle = noise * 3.14159265f * 2 * lookup_direction_turbulance;
				mr_vec.x = std::sinf(mr_angle) * mr_speed;
				mr_vec.y = -std::cosf(mr_angle) * mr_speed;
			}

			//----Weight Changes
			p.weight_acceleration += p.base.weight * lookup_weight * UPDATE_TIME;

			//----Velocity Changes
			tfxVec2 velocity_normal;
			velocity_normal.x = std::sinf(direction);
			velocity_normal.y = -std::cosf(direction);

			tfxVec2 current_velocity = (p.base.velocity * lookup_velocity) * velocity_normal;
			current_velocity += mr_vec;
			current_velocity.y += p.weight_acceleration;
			current_velocity *= UPDATE_TIME;

			//----Color changes
			p.color.a = unsigned char(255.f * lookup_opacity * global_intensity);
			p.intensity = lookup_intensity;
			if (!(flags & tfxEmitterStateFlags_random_color)) {
				p.color.r = unsigned char(255.f * lookup_red);
				p.color.g = unsigned char(255.f * lookup_green);
				p.color.b = unsigned char(255.f * lookup_blue);
			}

			p.color = tfxRGBA8(p.color.r, p.color.g, p.color.b, p.color.a);

			//----Size Changes
			tfxVec2 scale;
			scale.x = p.base.size.x * lookup_width;
			if (scale.x < 0.f)
				scale.x = scale.x;

			//----Stretch Changes
			float velocity = std::fabsf(lookup_velocity * p.base.velocity + mr_speed + p.weight_acceleration);
			if (flags & tfxEmitterStateFlags_lifetime_uniform_size) {
				scale.y = (lookup_width * (p.base.size.y + (velocity * lookup_stretch * stretch))) / image_size_y;
				if (scale.y < scale.x)
					scale.y = scale.x;
			}
			else
				scale.y = (lookup_height * (p.base.size.y + (velocity * lookup_stretch * stretch))) / image_size_y;

			//----Spin and angle Changes
			float spin = 0;
			if (flags & tfxEmitterStateFlags_can_spin) {
				spin = lookup_spin;
			}

			//---------------
			//Now that the latest changes are applied, affect the particle state
			//---------------

			//----Rotation
			if (flags & tfxEmitterStateFlags_align_with_velocity) {
				tfxVec2 vd = current_velocity.IsNill() ? velocity_normal : current_velocity;
				p.local.rotation = GetVectorAngle(vd.x, vd.y) + angle_offset;
			}
			else {
				p.local.rotation += spin * UPDATE_TIME;
			}

			//----Position
			p.local.position += current_velocity * overal_scale;

			//----Scale
			p.world.scale = scale;

			//Lines - Reposition if the particle is travelling along a line
			tfxVec2 offset = velocity_normal * emitter_size_y;
			float length = std::fabsf(p.local.position.y - emitter_handle_y);
			float emitter_length = emitter_size_y;
			bool line_and_kill = (flags & tfxEmitterStateFlags_is_line_traversal) && (flags & tfxEmitterStateFlags_kill) && length > emitter_length;
			bool line_and_loop = (flags & tfxEmitterStateFlags_is_line_traversal) && (flags & tfxEmitterStateFlags_loop) && length > emitter_length;
			if (line_and_loop) {
				p.local.position.y -= offset.y;
				p.flags |= tfxParticleFlags_capture_after_transform;
			}
			else if (line_and_kill) {
				p.flags |= tfxParticleFlags_remove;
			}

			//----Image animation
			p.image_frame += image_frame_rate * UPDATE_TIME;
			p.image_frame = (flags & tfxEmitterStateFlags_play_once) && p.image_frame > library_link->properties.end_frame ? p.image_frame = library_link->properties.end_frame : p.image_frame;
			p.image_frame = (flags & tfxEmitterStateFlags_play_once) && p.image_frame < 0 ? p.image_frame = 0 : p.image_frame;
			p.image_frame = std::fmodf(p.image_frame, library_link->properties.end_frame + 1);

			if (!(p.flags & tfxParticleFlags_fresh)) {
				TransformParticle(p, *this);
				if (p.flags & tfxParticleFlags_capture_after_transform) {
					p.captured.position = p.world.position;
					p.flags &= ~tfxParticleFlags_capture_after_transform;
				}

				ParticleSprite &s = root_effect.sprites[library_link->properties.layer][p.sprite_index];
				s.intensity = p.intensity;
				s.color = p.color;
				s.world = p.world;
				s.captured = p.captured;
				s.ptr = library_link->properties.image->ptr;
				s.handle = image_handle;
				s.particle = &particles[i - index_offset];
				s.parameters = (unsigned int)p.image_frame;
			}
			else {
				ParticleSprite &s = root_effect.sprites[library_link->properties.layer][p.sprite_index];
				s.particle = &particles[i - index_offset];
				p.flags &= ~tfxParticleFlags_fresh;
			}

			if (index_offset) {
				unsigned int tmp_offset = particles[i - index_offset].offset;
				Particle *tmp_particle = particles[i - index_offset].next_ptr;
				p.offset = index_offset;
				p.next_ptr = &particles[i - index_offset];
				particles[i - index_offset] = p;
				particles[i - index_offset].offset = tmp_offset;
				particles[i - index_offset].next_ptr = tmp_particle;
			}
			else {
				p.next_ptr = &p;
			}
		}
		particles.current_size -= index_offset;
	}

	EffectEmitter CreateEffector(float x, float y) {
		EffectEmitter effector;
		effector.transform.local.position = tfxVec2(x, y);

		return effector;
	}

	void Transform(tfxEmitter &e, Particle &parent) {
		float s = sin(e.transform.local.rotation);
		float c = cos(e.transform.local.rotation);

		e.transform.matrix.Set(c, s, -s, c);

		e.transform.world.rotation = parent.world.rotation + e.transform.local.rotation;

		e.transform.matrix = e.transform.matrix.Transform(parent.matrix);
		tfxVec2 rotatevec = parent.matrix.TransformVector(tfxVec2(e.transform.local.position.x, e.transform.local.position.y));

		e.transform.world.position = parent.world.position + rotatevec * parent.world.scale;

	}

	void Transform(tfxEmitter &e, tfxEffect &parent) {
		float s = sin(e.transform.local.rotation);
		float c = cos(e.transform.local.rotation);

		e.transform.matrix.Set(c, s, -s, c);
		e.transform.world.scale = parent.transform.world.scale;

		e.transform.world.rotation = parent.transform.world.rotation + e.transform.local.rotation;

		e.transform.matrix = e.transform.matrix.Transform(parent.transform.matrix);
		tfxVec2 rotatevec = parent.transform.matrix.TransformVector(tfxVec2(e.transform.local.position.x, e.transform.local.position.y));

		e.transform.world.position = parent.transform.world.position + rotatevec * parent.transform.world.scale;

	}

	void Transform(tfxEmitter &e, tfxEmitter &parent) {
		float s = sin(e.transform.local.rotation);
		float c = cos(e.transform.local.rotation);

		e.transform.matrix.Set(c, s, -s, c);
		e.transform.world.scale = parent.transform.world.scale;

		e.transform.world.rotation = parent.transform.world.rotation + e.transform.local.rotation;

		e.transform.matrix = e.transform.matrix.Transform(parent.transform.matrix);
		tfxVec2 rotatevec = parent.transform.matrix.TransformVector(tfxVec2(e.transform.local.position.x, e.transform.local.position.y));

		e.transform.world.position = parent.transform.world.position + rotatevec * parent.transform.world.scale;

	}

	void TransformParticle(Particle &p, EffectEmitter &e) {
		//The Particle matrix is only needed for sub effect transformations
		float s = sin(p.local.rotation);
		float c = cos(p.local.rotation);
		p.matrix.Set(c, s, -s, c);
		bool line = (e.properties.flags & tfxEmitterPropertyFlags_edge_traversal && e.properties.emission_type == tfxLine);

		if (e.properties.flags & tfxEmitterPropertyFlags_relative_position || line) {
			p.world.scale = p.world.scale;

			if (e.properties.flags & tfxEmitterPropertyFlags_relative_angle || line)
				p.world.rotation = e.transform.world.rotation + p.local.rotation;
			else
				p.world.rotation = p.local.rotation;

			p.matrix = p.matrix.Transform(e.transform.matrix);
			tfxVec2 rotatevec = e.transform.matrix.TransformVector(tfxVec2(p.local.position.x, p.local.position.y));

			p.world.position = e.transform.world.position + rotatevec * e.transform.world.scale;

		}
		else {
			p.world.position = p.local.position;
			p.world.scale = p.world.scale;
			if (e.properties.flags & tfxEmitterPropertyFlags_relative_angle)
				p.world.rotation = e.transform.world.rotation + p.local.rotation;
			else
				p.world.rotation = p.local.rotation;
		}

	}

	void TransformParticle(Particle &p, tfxEmitter &e) {
		//The Particle matrix is only needed for sub effect transformations
		float s = sin(p.local.rotation);
		float c = cos(p.local.rotation);
		p.matrix.Set(c, s, -s, c);
		bool line = (e.common.property_flags & tfxEmitterPropertyFlags_edge_traversal && e.library_link->properties.emission_type == tfxLine);

		if (e.common.property_flags & tfxEmitterPropertyFlags_relative_position || line) {
			p.world.scale = p.world.scale;

			if (e.common.property_flags & tfxEmitterPropertyFlags_relative_angle || line)
				p.world.rotation = e.transform.world.rotation + p.local.rotation;
			else
				p.world.rotation = p.local.rotation;

			p.matrix = p.matrix.Transform(e.transform.matrix);
			tfxVec2 rotatevec = e.transform.matrix.TransformVector(tfxVec2(p.local.position.x, p.local.position.y));

			p.world.position = e.transform.world.position + rotatevec * e.transform.world.scale;

		}
		else {
			p.world.position = p.local.position;
			p.world.scale = p.world.scale;
			if (e.common.property_flags & tfxEmitterPropertyFlags_relative_angle)
				p.world.rotation = e.transform.world.rotation + p.local.rotation;
			else
				p.world.rotation = p.local.rotation;
		}

	}

	void Transform(FormState &local, FormState &world, EffectEmitter &e) {
		//The Particle matrix is only needed for sub effect transformations
		bool line = (e.properties.flags & tfxEmitterPropertyFlags_edge_traversal && e.properties.emission_type == tfxLine);

		if (e.properties.flags & tfxEmitterPropertyFlags_relative_position || line) {
			world.scale = local.scale;

			if (e.properties.flags & tfxEmitterPropertyFlags_relative_angle || line)
				world.rotation = e.transform.world.rotation + local.rotation;
			else
				world.rotation = local.rotation;

			tfxVec2 rotatevec = e.transform.matrix.TransformVector(tfxVec2(local.position.x, local.position.y));

			world.position = e.transform.world.position + rotatevec * e.transform.world.scale;

		}
		else {
			world.position = local.position;
			world.scale = local.scale;
			if (e.properties.flags & tfxEmitterPropertyFlags_relative_angle)
				world.rotation = e.transform.world.rotation + local.rotation;
			else
				world.rotation = local.rotation;
		}
	}

	void EffectEmitter::ResetGlobalGraphs(bool add_node) {
		library->global_graphs[global].life.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].life.type = tfxGlobal_life;
		library->global_graphs[global].amount.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].amount.type = tfxGlobal_amount;
		library->global_graphs[global].velocity.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].velocity.type = tfxGlobal_velocity;
		library->global_graphs[global].width.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].width.type = tfxGlobal_width;
		library->global_graphs[global].height.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].height.type = tfxGlobal_height;
		library->global_graphs[global].weight.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].weight.type = tfxGlobal_weight;
		library->global_graphs[global].spin.Reset(1.f, tfxGlobalPercentPresetSigned, add_node); library->global_graphs[global].spin.type = tfxGlobal_spin;
		library->global_graphs[global].stretch.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].stretch.type = tfxGlobal_stretch;
		library->global_graphs[global].overal_scale.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].overal_scale.type = tfxGlobal_overal_scale;
		library->global_graphs[global].intensity.Reset(1.f, tfxGlobalOpacityPreset, add_node); library->global_graphs[global].intensity.type = tfxGlobal_intensity;
		library->global_graphs[global].frame_rate.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].frame_rate.type = tfxGlobal_frame_rate;
		library->global_graphs[global].splatter.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].splatter.type = tfxGlobal_splatter;
		library->global_graphs[global].effect_angle.Reset(0.f, tfxAnglePreset, add_node); library->global_graphs[global].effect_angle.type = tfxGlobal_effect_angle;
		library->CompileGlobalGraph(global);
	}

	void EffectEmitter::ResetBaseGraphs(bool add_node) {
		library->base_graphs[base].life.Reset(1000.f, tfxLifePreset, add_node); library->base_graphs[base].life.type = tfxBase_life;
		library->base_graphs[base].amount.Reset(1.f, tfxAmountPreset, add_node); library->base_graphs[base].amount.type = tfxBase_amount;
		library->base_graphs[base].velocity.Reset(0.f, tfxVelocityPreset, add_node); library->base_graphs[base].velocity.type = tfxBase_velocity;
		library->base_graphs[base].width.Reset(128.f, tfxDimensionsPreset, add_node); library->base_graphs[base].width.type = tfxBase_width;
		library->base_graphs[base].height.Reset(128.f, tfxDimensionsPreset, add_node); library->base_graphs[base].height.type = tfxBase_height;
		library->base_graphs[base].weight.Reset(0.f, tfxWeightPreset, add_node); library->base_graphs[base].weight.type = tfxBase_weight;
		library->base_graphs[base].spin.Reset(0.f, tfxSpinPreset, add_node); library->base_graphs[base].spin.type = tfxBase_spin;
		library->base_graphs[base].noise_offset.Reset(0.f, tfxGlobalPercentPreset, add_node); library->base_graphs[base].noise_offset.type = tfxBase_noise_offset;
		library->CompileBaseGraph(base);
	}

	void EffectEmitter::ResetPropertyGraphs(bool add_node) {
		library->property_graphs[property].emission_angle.Reset(0.f, tfxAnglePreset, add_node); library->property_graphs[property].emission_angle.type = tfxProperty_emission_angle;
		library->property_graphs[property].emission_range.Reset(0.f, tfxEmissionRangePreset, add_node); library->property_graphs[property].emission_range.type = tfxProperty_emission_range;
		library->property_graphs[property].emitter_angle.Reset(0.f, tfxAnglePreset, add_node); library->property_graphs[property].emitter_angle.type = tfxProperty_emitter_angle;
		library->property_graphs[property].splatter.Reset(0.f, tfxDimensionsPreset, add_node); library->property_graphs[property].splatter.type = tfxProperty_splatter;
		library->property_graphs[property].emitter_width.Reset(0.f, tfxDimensionsPreset, add_node); library->property_graphs[property].emitter_width.type = tfxProperty_emitter_width;
		library->property_graphs[property].emitter_height.Reset(0.f, tfxDimensionsPreset, add_node); library->property_graphs[property].emitter_height.type = tfxProperty_emitter_height;
		library->property_graphs[property].arc_size.Reset(tfxRadians(360.f), tfxArcPreset, add_node); library->property_graphs[property].arc_size.type = tfxProperty_arc_size;
		library->property_graphs[property].arc_offset.Reset(0.f, tfxArcPreset, add_node); library->property_graphs[property].arc_offset.type = tfxProperty_arc_offset;
		library->CompilePropertyGraph(property);
	}

	void EffectEmitter::ResetVariationGraphs(bool add_node) {
		library->variation_graphs[variation].life.Reset(0.f, tfxLifePreset, add_node); library->variation_graphs[variation].life.type = tfxVariation_life;
		library->variation_graphs[variation].amount.Reset(0.f, tfxAmountPreset, add_node); library->variation_graphs[variation].amount.type = tfxVariation_amount;
		library->variation_graphs[variation].velocity.Reset(0.f, tfxVelocityPreset, add_node); library->variation_graphs[variation].velocity.type = tfxVariation_velocity;
		library->variation_graphs[variation].width.Reset(0.f, tfxDimensionsPreset, add_node); library->variation_graphs[variation].width.type = tfxVariation_width;
		library->variation_graphs[variation].height.Reset(0.f, tfxDimensionsPreset, add_node); library->variation_graphs[variation].height.type = tfxVariation_height;
		library->variation_graphs[variation].weight.Reset(0.f, tfxWeightVariationPreset, add_node); library->variation_graphs[variation].weight.type = tfxVariation_weight;
		library->variation_graphs[variation].spin.Reset(0.f, tfxSpinVariationPreset, add_node); library->variation_graphs[variation].spin.type = tfxVariation_spin;
		library->variation_graphs[variation].noise_offset.Reset(0.f, tfxNoiseOffsetVariationPreset, add_node); library->variation_graphs[variation].noise_offset.type = tfxVariation_noise_offset;
		library->variation_graphs[variation].noise_resolution.Reset(300.f, tfxNoiseResolutionPreset, add_node); library->variation_graphs[variation].noise_resolution.type = tfxVariation_noise_resolution;
		library->CompileVariationGraph(variation);
	}

	void EffectEmitter::ResetOvertimeGraphs(bool add_node) {
		library->overtime_graphs[overtime].velocity.Reset(1.f, tfxVelocityOvertimePreset, add_node); library->overtime_graphs[overtime].velocity.type = tfxOvertime_velocity;
		library->overtime_graphs[overtime].velocity_adjuster.Reset(1.f, tfxGlobalPercentPreset, add_node); library->overtime_graphs[overtime].velocity_adjuster.type = tfxOvertime_velocity_adjuster;
		library->overtime_graphs[overtime].width.Reset(1.f, tfxPercentOvertime, add_node); library->overtime_graphs[overtime].width.type = tfxOvertime_width;
		library->overtime_graphs[overtime].height.Reset(1.f, tfxPercentOvertime, add_node); library->overtime_graphs[overtime].height.type = tfxOvertime_height;
		library->overtime_graphs[overtime].weight.Reset(1.f, tfxPercentOvertime, add_node); library->overtime_graphs[overtime].weight.type = tfxOvertime_weight;
		library->overtime_graphs[overtime].spin.Reset(0.f, tfxSpinOvertimePreset, add_node); library->overtime_graphs[overtime].spin.type = tfxOvertime_spin;
		library->overtime_graphs[overtime].stretch.Reset(0.f, tfxPercentOvertime, add_node); library->overtime_graphs[overtime].stretch.type = tfxOvertime_stretch;
		library->overtime_graphs[overtime].red.Reset(1.f, tfxColorPreset, add_node); library->overtime_graphs[overtime].red.type = tfxOvertime_red;
		library->overtime_graphs[overtime].green.Reset(1.f, tfxColorPreset, add_node); library->overtime_graphs[overtime].green.type = tfxOvertime_green;
		library->overtime_graphs[overtime].blue.Reset(1.f, tfxColorPreset, add_node); library->overtime_graphs[overtime].blue.type = tfxOvertime_blue;
		library->overtime_graphs[overtime].blendfactor.Reset(1.f, tfxOpacityOvertimePreset, add_node); library->overtime_graphs[overtime].blendfactor.type = tfxOvertime_blendfactor;
		library->overtime_graphs[overtime].intensity.Reset(1.f, tfxIntensityOvertimePreset, add_node); library->overtime_graphs[overtime].intensity.type = tfxOvertime_intensity;
		library->overtime_graphs[overtime].velocity_turbulance.Reset(30.f, tfxFrameratePreset, add_node); library->overtime_graphs[overtime].velocity_turbulance.type = tfxOvertime_velocity_turbulance;
		library->overtime_graphs[overtime].stretch.Reset(0.f, tfxPercentOvertime, add_node); library->overtime_graphs[overtime].stretch.type = tfxOvertime_stretch;
		library->overtime_graphs[overtime].direction_turbulance.Reset(0.f, tfxPercentOvertime, add_node); library->overtime_graphs[overtime].direction_turbulance.type = tfxOvertime_direction_turbulance;
		library->overtime_graphs[overtime].direction.Reset(0.f, tfxDirectionOvertimePreset, add_node); library->overtime_graphs[overtime].direction.type = tfxOvertime_direction;
		library->overtime_graphs[overtime].noise_resolution.Reset(1.f, tfxPercentOvertime, add_node); library->overtime_graphs[overtime].noise_resolution.type = tfxOvertime_noise_resolution;
		library->CompileOvertimeGraph(overtime);
	}

	void EffectEmitter::ResetEffectGraphs(bool add_node) {
		ResetGlobalGraphs(add_node);
	}

	void EffectEmitter::ResetEmitterGraphs(bool add_node) {
		ResetBaseGraphs(add_node);
		ResetPropertyGraphs(add_node);
		ResetVariationGraphs(add_node);
		UpdateMaxLife();
		ResetOvertimeGraphs(add_node);
	}

	void EffectEmitter::InitialiseUninitialisedGraphs() {
		if (type == tfxEffectType) {
			if (library->global_graphs[global].life.nodes.size() == 0) library->global_graphs[global].life.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].amount.nodes.size() == 0) library->global_graphs[global].amount.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].velocity.nodes.size() == 0) library->global_graphs[global].velocity.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].width.nodes.size() == 0) library->global_graphs[global].width.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].height.nodes.size() == 0) library->global_graphs[global].height.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].weight.nodes.size() == 0) library->global_graphs[global].weight.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].spin.nodes.size() == 0) library->global_graphs[global].spin.Reset(1.f, tfxGlobalPercentPresetSigned);
			if (library->global_graphs[global].effect_angle.nodes.size() == 0) library->global_graphs[global].effect_angle.Reset(0.f, tfxAnglePreset);
			if (library->global_graphs[global].stretch.nodes.size() == 0) library->global_graphs[global].stretch.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].overal_scale.nodes.size() == 0) library->global_graphs[global].overal_scale.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].intensity.nodes.size() == 0) library->global_graphs[global].intensity.Reset(1.f, tfxOpacityOvertimePreset);
			if (library->global_graphs[global].frame_rate.nodes.size() == 0) library->global_graphs[global].frame_rate.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].splatter.nodes.size() == 0) library->global_graphs[global].splatter.Reset(1.f, tfxGlobalPercentPreset);
		}

		if (type == tfxEmitterType) {
			if (library->base_graphs[base].life.nodes.size() == 0) library->base_graphs[base].life.Reset(1000.f, tfxLifePreset);
			if (library->base_graphs[base].amount.nodes.size() == 0) library->base_graphs[base].amount.Reset(1.f, tfxAmountPreset);
			if (library->base_graphs[base].velocity.nodes.size() == 0) library->base_graphs[base].velocity.Reset(0.f, tfxVelocityPreset);
			if (library->base_graphs[base].width.nodes.size() == 0) library->base_graphs[base].width.Reset(128.f, tfxDimensionsPreset);
			if (library->base_graphs[base].height.nodes.size() == 0) library->base_graphs[base].height.Reset(128.f, tfxDimensionsPreset);
			if (library->base_graphs[base].weight.nodes.size() == 0) library->base_graphs[base].weight.Reset(0.f, tfxWeightPreset);
			if (library->base_graphs[base].spin.nodes.size() == 0) library->base_graphs[base].spin.Reset(0.f, tfxSpinPreset);
			if (library->base_graphs[base].noise_offset.nodes.size() == 0) library->base_graphs[base].noise_offset.Reset(0.f, tfxGlobalPercentPreset);

			if (library->property_graphs[property].emitter_angle.nodes.size() == 0) library->property_graphs[property].emitter_angle.Reset(0.f, tfxAnglePreset);
			if (library->property_graphs[property].emission_angle.nodes.size() == 0) library->property_graphs[property].emission_angle.Reset(0.f, tfxAnglePreset);
			if (library->property_graphs[property].emission_range.nodes.size() == 0) library->property_graphs[property].emission_range.Reset(0.f, tfxEmissionRangePreset);
			if (library->property_graphs[property].splatter.nodes.size() == 0) library->property_graphs[property].splatter.Reset(0.f, tfxDimensionsPreset);
			if (library->property_graphs[property].emitter_width.nodes.size() == 0) library->property_graphs[property].emitter_width.Reset(0.f, tfxDimensionsPreset);
			if (library->property_graphs[property].emitter_height.nodes.size() == 0) library->property_graphs[property].emitter_height.Reset(0.f, tfxDimensionsPreset);
			if (library->property_graphs[property].arc_size.nodes.size() == 0) library->property_graphs[property].arc_size.Reset(tfxRadians(360.f), tfxArcPreset);
			if (library->property_graphs[property].arc_offset.nodes.size() == 0) library->property_graphs[property].arc_offset.Reset(0.f, tfxArcPreset);

			if (library->variation_graphs[variation].life.nodes.size() == 0) library->variation_graphs[variation].life.Reset(0.f, tfxLifePreset);
			if (library->variation_graphs[variation].amount.nodes.size() == 0) library->variation_graphs[variation].amount.Reset(0.f, tfxAmountPreset);
			if (library->variation_graphs[variation].velocity.nodes.size() == 0) library->variation_graphs[variation].velocity.Reset(0.f, tfxVelocityPreset);
			if (library->variation_graphs[variation].weight.nodes.size() == 0) library->variation_graphs[variation].weight.Reset(0.f, tfxVelocityPreset);
			if (library->variation_graphs[variation].width.nodes.size() == 0) library->variation_graphs[variation].width.Reset(0.f, tfxDimensionsPreset);
			if (library->variation_graphs[variation].height.nodes.size() == 0) library->variation_graphs[variation].height.Reset(0.f, tfxDimensionsPreset);
			if (library->variation_graphs[variation].weight.nodes.size() == 0) library->variation_graphs[variation].weight.Reset(0.f, tfxWeightVariationPreset);
			if (library->variation_graphs[variation].spin.nodes.size() == 0) library->variation_graphs[variation].spin.Reset(0.f, tfxSpinVariationPreset);
			if (library->variation_graphs[variation].noise_offset.nodes.size() == 0) library->variation_graphs[variation].noise_offset.Reset(0.f, tfxNoiseOffsetVariationPreset);
			if (library->variation_graphs[variation].noise_resolution.nodes.size() == 0) library->variation_graphs[variation].noise_resolution.Reset(300.f, tfxNoiseResolutionPreset);

			if (library->overtime_graphs[overtime].velocity.nodes.size() == 0) library->overtime_graphs[overtime].velocity.Reset(1.f, tfxVelocityOvertimePreset);
			if (library->overtime_graphs[overtime].width.nodes.size() == 0) library->overtime_graphs[overtime].width.Reset(1.f, tfxPercentOvertime);
			if (library->overtime_graphs[overtime].height.nodes.size() == 0) library->overtime_graphs[overtime].height.Reset(1.f, tfxPercentOvertime);
			if (library->overtime_graphs[overtime].weight.nodes.size() == 0) library->overtime_graphs[overtime].weight.Reset(1.f, tfxPercentOvertime);
			if (library->overtime_graphs[overtime].spin.nodes.size() == 0) library->overtime_graphs[overtime].spin.Reset(1.f, tfxSpinOvertimePreset);
			if (library->overtime_graphs[overtime].stretch.nodes.size() == 0) library->overtime_graphs[overtime].stretch.Reset(0.f, tfxPercentOvertime);
			if (library->overtime_graphs[overtime].red.nodes.size() == 0) library->overtime_graphs[overtime].red.Reset(1.f, tfxColorPreset);
			if (library->overtime_graphs[overtime].green.nodes.size() == 0) library->overtime_graphs[overtime].green.Reset(1.f, tfxColorPreset);
			if (library->overtime_graphs[overtime].blue.nodes.size() == 0) library->overtime_graphs[overtime].blue.Reset(1.f, tfxColorPreset);
			if (library->overtime_graphs[overtime].blendfactor.nodes.size() == 0) library->overtime_graphs[overtime].blendfactor.Reset(1.f, tfxOpacityOvertimePreset);
			if (library->overtime_graphs[overtime].intensity.nodes.size() == 0) library->overtime_graphs[overtime].intensity.Reset(1.f, tfxIntensityOvertimePreset);
			if (library->overtime_graphs[overtime].velocity_turbulance.nodes.size() == 0) library->overtime_graphs[overtime].velocity_turbulance.Reset(0.f, tfxFrameratePreset);
			if (library->overtime_graphs[overtime].direction_turbulance.nodes.size() == 0) library->overtime_graphs[overtime].direction_turbulance.Reset(0.f, tfxPercentOvertime);
			if (library->overtime_graphs[overtime].velocity_adjuster.nodes.size() == 0) library->overtime_graphs[overtime].velocity_adjuster.Reset(1.f, tfxGlobalPercentPreset);
			if (library->overtime_graphs[overtime].direction.nodes.size() == 0) library->overtime_graphs[overtime].direction.Reset(0.f, tfxDirectionOvertimePreset);
			if (library->overtime_graphs[overtime].noise_resolution.nodes.size() == 0) library->overtime_graphs[overtime].noise_resolution.Reset(1.f, tfxPercentOvertime);
		}
	}

	void EffectEmitter::SetName(const char *n) {
		name = n;
	}

	FormState Tween(float tween, FormState &world, FormState &captured) {
		FormState tweened;
		tweened.position = world.position * tween + captured.position * (1.f - tween);
		tweened.scale = world.scale * tween + captured.scale * (1.f - tween);
		//Not tweening rotation for now, need to figure out when it tweens over 180 degrees.
		//tweened.rotation = world.rotation * tween + captured.rotation * (1.f - tween);
		tweened.rotation = world.rotation;

		return tweened;
	}

	void EffectEmitter::ClearColors() {
		library->overtime_graphs[overtime].red.Clear();
		library->overtime_graphs[overtime].green.Clear();
		library->overtime_graphs[overtime].blue.Clear();
	}

	void EffectEmitter::AddColorOvertime(float frame, tfxRGB color) {
		library->overtime_graphs[overtime].red.AddNode(frame, color.r);
		library->overtime_graphs[overtime].green.AddNode(frame, color.g);
		library->overtime_graphs[overtime].blue.AddNode(frame, color.b);
	}

	void EffectEmitter::ReSeed(uint64_t seed) {
		if (seed == 0) {
			seed = 0xFFFFFFFFFFF;
		}
		random_generation.ReSeed(seed, seed / 2);
	}

	 void EffectEmitter::SetUserData(void *data) {
		user_data = data;
	}

	void* EffectEmitter::GetUserData() {
		return user_data;
	}

	void EffectEmitter::SetTimeout(unsigned int frames) {
		timeout = frames;
		for (auto &sub : sub_effectors) {
			sub.SetTimeout(frames);
		}
	}

	bool EffectEmitter::HasSingle() {
		for (auto &e : sub_effectors) {
			if (e.properties.flags & (tfxEmitterPropertyFlags_single | tfxEmitterPropertyFlags_one_shot))
				return true;
		}
		return false;
	}

	bool EffectEmitter::RenameSubEffector(EffectEmitter &emitter, const char *new_name) {
		if (!NameExists(emitter, new_name) && strlen(new_name) > 0) {
			emitter.SetName(new_name);
			library->UpdateEffectPaths();
			return true;
		}

		return false;
	}

	bool EffectEmitter::NameExists(EffectEmitter &emitter, const char *name) {
		for (auto &e : sub_effectors) {
			if (&emitter != &e) {
				if (e.name == name) {
					return true;
				}
			}
		}

		return false;
	}

	void EffectEmitter::ReIndex() {
		unsigned int index = 0;
		for (auto &e : sub_effectors) {
			e.library_index = index++;
			e.parent = this;
			e.ReIndex();
		}
	}

	void EffectEmitter::CountChildren(int &emitters, int &effects) {
		tfxvec<EffectEmitter*> stack;
		stack.push_back(this);
		emitters = 0;
		effects = 0;
		while (!stack.empty()) {
			EffectEmitter *current = stack.pop_back();
			if (current->type == tfxEffectType)
				effects++;
			else if (current->type == tfxEmitterType)
				emitters++;
			for (auto &sub : current->sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	EffectEmitter* EffectEmitter::GetRootEffect() {
		if (!parent || parent->type == tfxFolder)
			return nullptr;
		EffectEmitter *p = parent;
		unsigned int timeout = 0;
		while (p || ++timeout < 100) {
			if (!p->parent)
				return p;
			p = p->parent;
		}
		return nullptr;
	}

	bool EffectEmitter::IsRootEffect() {
		if (type != tfxEffectType) return false;
		if (type == tfxEffectType && !parent) return true;
		if (parent && parent->type == tfxFolder) return true;
		return false;
	}

	void EffectEmitter::ResetParents() {
		parent = nullptr;
		for (auto &e : sub_effectors) {
			e.ResetParents();
		}
	}

	EffectEmitter* EffectEmitter::MoveUp(EffectEmitter &emitter) {
		if (emitter.library_index > 0) {
			unsigned int new_index = emitter.library_index - 1;
			std::swap(sub_effectors[emitter.library_index], sub_effectors[new_index]);
			ReIndex();
			emitter.library->UpdateEffectPaths();
			return &sub_effectors[new_index];
		}

		return nullptr;
	}

	EffectEmitter* EffectEmitter::MoveDown(EffectEmitter &emitter) {
		if (emitter.library_index < sub_effectors.size() - 1) {
			unsigned int new_index = emitter.library_index + 1;
			std::swap(sub_effectors[emitter.library_index], sub_effectors[new_index]);
			ReIndex();
			emitter.library->UpdateEffectPaths();
			return &sub_effectors[new_index];
		}
		return nullptr;
	}

	void EffectEmitter::DeleteEmitter(EffectEmitter *emitter) {
		EffectLibrary *library = emitter->library;
		tfxvec<EffectEmitter> stack;
		stack.push_back(*emitter);
		while (stack.size()) {
			EffectEmitter &current = stack.pop_back();
			if (current.type == tfxEffectType && !current.parent) {
				library->FreeGlobal(current.global);
			}
			else if(current.type == tfxEmitterType) {
				library->FreeProperty(current.property);
				library->FreeBase(current.base);
				library->FreeVariation(current.variation);
				library->FreeOvertime(current.overtime);
			}
			for (auto &sub : current.sub_effectors) {
				stack.push_back(sub);
			}
		}
		sub_effectors.erase(emitter);

		ReIndex();
		if(library)
			library->UpdateEffectPaths();
	}

	void EffectEmitter::CleanUp() {
		if (sub_effectors.size()) {
			tfxvec<EffectEmitter> stack;
			stack.push_back(*this);
			while (stack.size()) {
				EffectEmitter current = stack.pop_back();
				if (current.type == tfxEffectType && !current.parent) {
					library->FreeGlobal(current.global);
				}
				else if(current.type == tfxEmitterType) {
					library->FreeProperty(current.property);
					library->FreeBase(current.base);
					library->FreeVariation(current.variation);
					library->FreeOvertime(current.overtime);
				}
				for (auto &sub : current.sub_effectors) {
					stack.push_back(sub);
				}
				current.sub_effectors.clear();
			}
		}

		ReIndex();
	}

	void EffectEmitter::Clone(EffectEmitter &clone, EffectEmitter *root_parent, EffectLibrary *destination_library, bool keep_user_data, bool force_clone_global) {
		clone = *this;
		clone.sub_effectors.clear();
		clone.flags |= tfxEmitterStateFlags_enabled;
		if(!keep_user_data)
			clone.user_data = nullptr;
		clone.library = destination_library;

		if (type == tfxEffectType) {
			if (root_parent == &clone) {
				clone.global = library->CloneGlobal(global, destination_library);
				clone.library->CompileGlobalGraph(clone.global);
			}
			else {
				if (!force_clone_global) {
					clone.global = root_parent->global;
				}
				else {
					clone.global = library->CloneGlobal(root_parent->global, destination_library);
					clone.library->CompileGlobalGraph(clone.global);
				}
			}
		}
		else if(type == tfxEmitterType) {
			clone.property = library->CloneProperty(property, destination_library);
			clone.base = library->CloneBase(base, destination_library);
			clone.variation = library->CloneVariation(variation, destination_library);
			clone.overtime = library->CloneOvertime(overtime, destination_library);
			clone.UpdateMaxLife();
			clone.library->CompilePropertyGraph(clone.property);
			clone.library->CompileBaseGraph(clone.base);
			clone.library->CompileVariationGraph(clone.variation);
			clone.library->CompileOvertimeGraph(clone.overtime);
		}

		for (auto &e : sub_effectors) {
			if (e.type == tfxEmitterType) {
				EffectEmitter emitter_copy;
				e.Clone(emitter_copy, root_parent, destination_library);
				if(!keep_user_data)
					emitter_copy.user_data = nullptr;
				clone.AddEmitter(emitter_copy);
			}
			else if(e.type == tfxEffectType) {
				EffectEmitter effect_copy;
				if(clone.type == tfxFolder)
					e.Clone(effect_copy, &effect_copy, destination_library);
				else
					e.Clone(effect_copy, root_parent, destination_library);
				if(!keep_user_data)
					effect_copy.user_data = nullptr;
				clone.AddEffect(effect_copy);
			}
		}
	}

	void InitEffectPool(tfxEffectPool &effect_pool, unsigned int max_effects, unsigned int max_emitters, unsigned int max_particles) {
		effect_pool.Init(max_effects, max_emitters, max_particles);
	}

	bool PoolEffectFromLibrary(EffectLibrary &library, tfxEffectPool &storage, const char *name, tfxEffectID &out) {
		assert(library.effect_paths.ValidName(name));
		EffectEmitter *effect = library.GetEffect(name);
		if (!Copy(storage, *effect, out))
			return false;
		return true;
	}

	bool PoolEffectFromLibrary(EffectLibrary &library, tfxEffectPool &storage, tfxKey path_hash, tfxEffectID &out) {
		assert(library.effect_paths.ValidKey(path_hash));
		EffectEmitter *effect = library.GetEffect(path_hash);
		if (!Copy(storage, *effect, out))
			return false;
		return true;
	}

	bool PoolEffectFromLibrary(EffectLibrary &library, tfxEffectPool &storage, unsigned int index, tfxEffectID &out) {
		assert(library.effects.size() > index);
		//if overwriting an effect that is already in storage, we need to add it to the cache
		if (!Copy(storage, library.effects[index], out))
			return false;
		return true;
	}

	bool PoolEffect(tfxEffectPool &storage, EffectEmitter &effect, tfxEffectID &out) {
		if (!Copy(storage, effect, out))
			return false;
		return true;
	}

	bool PrepareEffectTemplate(EffectLibrary &library, const char *name, tfxEffectTemplate &effect_template) {
		if (library.effect_paths.ValidName(name)) {
			library.PrepareEffectTemplate(name, effect_template);
			return true;
		}
		return false;
	}

	bool PoolEffectFromTemplate(tfxEffectPool &effect_pool, tfxEffectTemplate &effect_template, tfxEffectID &effect_id) {
		if (!Copy(effect_pool, effect_template.effect_template, effect_id))
			return false;
		return true;
	}

	void tfxEffect::Reset() {
		common.age = 0.f;
		common.highest_particle_age = 0.f;
		common.timeout_counter = 0;
		ClearSprites();
		for (auto &e : sub_effects) {
			if (e.type == tfxEmitterType) {
				e.particles.clear();
				e.particles.free_range(storage->particle_memory);
			}
		}
		sub_effects.clear();
		for (auto &e : sub_emitters) {
			e.Reset();
		}
	}

	void tfxEffect::ClearSprites() {
		for (int i = 0; i != tfxLAYERS; ++i) {
			sprites[i].clear();
		}
	}

	void tfxEffect::CompressSprites() {
		for (EachLayer) {
			unsigned int offset = 0;
			for (int si = 0; si != sprites[layer].current_size;++si) {
				if (sprites[layer][si].parameters == tfxINVALID) {
					offset++;
				}
				else if(offset > 0) {
					sprites[layer][si].particle->sprite_index -= offset;
					sprites[layer][si - offset] = sprites[layer][si];
				}
			}
			sprites[layer].current_size -= offset;
		}
	}

	void tfxEffect::ReleaseMemory() {
		common.age = 0.f;
		common.highest_particle_age = 0.f;
		common.timeout_counter = 0;
		for (int i = 0; i != tfxLAYERS; ++i) {
			sprites[i].free_range(storage->sprite_memory);
		}
		for (auto &e : sub_effects) {
			if (e.type == tfxEmitterType) {
				e.particles.free_range(storage->particle_memory);
			}
		}
		sub_effects.free_range(storage->emitter_memory);
		for (auto &e : sub_emitters) {
			e.particles.free_range(storage->particle_memory);
			e.Reset();
		}
		sub_emitters.free_range(storage->emitter_memory);
		storage->effect_memory.free_range(id);
	}

	void tfxEmitter::Reset() {
		transform.world.scale.x = transform.world.scale.y = 1.f;
		common.age = 0.f;
		common.state_flags &= ~tfxEmitterStateFlags_single_shot_done;
		common.state_flags &= ~tfxEmitterStateFlags_stop_spawning;
		current.emission_alternator = 0.f;
		current.amount_remainder = 0.f;
		current.grid_coords = tfxVec2();
		current.grid_direction = tfxVec2();
		parent_particle = nullptr;
		particles.clear();
	}

	bool Copy(tfxEffectPool &storage, EffectEmitter &in, tfxEffectID &out) {
		//If it's an effect that's in use then we need to check if there's enough memory to store it, otherwise
		//we can just overwrite it.
		if (!ValidEffect(storage, out)) {
			int emitters = 0;
			int effects = 0;
			int particles = 0;
			for (int layer = 0; layer != tfxLAYERS; ++layer) {
				particles += in.max_particles[layer];
			}
			unsigned int emitter_mem_req = emitters * in.max_sub_emitters * sizeof(tfxEmitter);
			unsigned int effect_mem_req = sizeof(tfxEffect);
			unsigned int particle_mem_req = particles * sizeof(Particle);
			unsigned int sprite_mem_req = particles * sizeof(ParticleSprite);
			if (storage.emitter_memory.free_unused_space() < emitter_mem_req) 
				return false;
			if (storage.effect_memory.free_unused_space() < effect_mem_req)
				return false;
			if (storage.particle_memory.free_unused_space() < particle_mem_req)
				return false;
			if (storage.sprite_memory.free_unused_space() < sprite_mem_req)
				return false;
		}
		in.CopyToEffect(out, storage);
		return true;
	}

	void EffectEmitter::CopyToEffect(tfxEffectID &effect_id, tfxEffectPool &storage) {
		assert(type == tfxEffectType);		//Must be called on an effect type only
		if (ValidEffect(storage, effect_id)) {
			tfxEffect &tmp_effect = storage.GetEffect(effect_id);
			tmp_effect.ReleaseMemory();
		}
		else {
			effect_id = storage.InsertEffect();
		}
		tfxEffect &effect = storage.GetEffect(effect_id);
		effect.id = effect_id;
		effect.library_link = this;
		effect.common.property_flags = properties.flags;
		effect.common.state_flags = flags;
		effect.common.library = library;
		effect.path_hash = path_hash;
		effect.common.loop_length = properties.loop_length;
		effect.common.handle = properties.emitter_handle;
		effect.lookup_mode = tfxFast;
		effect.storage = &storage;
		effect.user_data = user_data;
		effect.update_callback = root_effect_update_callback;
		for (int i = 0; i != tfxLAYERS; ++i) {
			effect.max_particles[i] = max_particles[i];
			effect.sprites[i].assign_memory(storage.sprite_memory, sizeof(ParticleSprite), effect.max_particles[i]);
		}
		effect.sub_emitters.assign_memory(storage.emitter_memory, sizeof(tfxEmitter), sub_effectors.size());
		effect.sub_effects.assign_memory(storage.emitter_memory, sizeof(tfxEmitter), max_sub_emitters);
		for (auto &emitter : sub_effectors) {
			tfxEmitter new_emitter;
			new_emitter.common.root_effect = &effect;
			emitter.CopyToEmitter(new_emitter, storage, true);
			assert(effect.sub_emitters.push_back(new_emitter));
		}
		effect.Reset();
	}

	void EffectEmitter::CopyToEmitter(tfxEmitter &e, tfxEffectPool &storage, bool assign_memory) {
		//Internal use only, CopyToEffect must be called first
		e.library_link = this;
		e.type = tfxEmitterType;
		e.common.property_flags = properties.flags;
		e.common.state_flags = flags;
		e.common.loop_length = properties.loop_length;	
		e.common.library = library;
		e.common.handle = properties.emitter_handle;
		e.common.age = 0.f;
		if(assign_memory)
			e.particles.assign_memory(storage.particle_memory, sizeof(Particle), max_particles[properties.layer]);
		
		e.lookup_mode = tfxFast;
		if (properties.flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
			e.image_handle = tfxVec2(0.5f, 0.5f);
		}
		else {
			e.image_handle = properties.image_handle;
		}
		e.common.state_flags &= ~tfxEmitterStateFlags_retain_matrix;

		e.common.state_flags |= properties.flags & tfxEmitterPropertyFlags_single && !(properties.flags & tfxEmitterPropertyFlags_one_shot) ? tfxEmitterStateFlags_is_single : 0;
		e.common.state_flags |= (properties.emission_type != tfxLine && !(properties.flags & tfxEmitterPropertyFlags_edge_traversal)) || properties.emission_type == tfxLine && !(properties.flags & tfxEmitterPropertyFlags_edge_traversal) ? tfxEmitterStateFlags_not_line : 0;
		e.common.state_flags |= properties.flags & tfxEmitterPropertyFlags_random_color;
		e.common.state_flags |= properties.flags & tfxEmitterPropertyFlags_lifetime_uniform_size;
		e.common.state_flags |= properties.angle_setting != AngleSetting::tfxAlign && !(properties.flags & tfxEmitterPropertyFlags_relative_angle) ? tfxEmitterStateFlags_can_spin : 0;
		e.common.state_flags |= properties.angle_setting == AngleSetting::tfxAlign ? tfxEmitterStateFlags_align_with_velocity : 0;
		e.common.state_flags |= properties.emission_type == tfxLine && properties.flags & tfxEmitterPropertyFlags_edge_traversal ? tfxEmitterStateFlags_is_line_traversal : 0;
		e.common.state_flags |= properties.flags & tfxEmitterPropertyFlags_play_once;
		e.common.state_flags |= properties.end_behaviour == tfxLoop ? tfxEmitterStateFlags_loop : 0;
		e.common.state_flags |= properties.end_behaviour == tfxKill ? tfxEmitterStateFlags_kill : 0;
		e.common.state_flags |= properties.emission_type == tfxArea || properties.emission_type == tfxEllipse ? tfxEmitterStateFlags_is_area : 0;
		e.common.state_flags |= properties.emission_type == tfxLine ? tfxEmitterStateFlags_is_line : 0;
		e.current.grid_coords = tfxVec2();
		
	}

	void EffectEmitter::EnableAllEmitters() {
		for (auto &e : sub_effectors) {
			e.flags |= tfxEmitterStateFlags_enabled;
			e.EnableAllEmitters();
		}
	}

	void EffectEmitter::EnableEmitter() {
		flags |= tfxEmitterStateFlags_enabled;
	}

	void EffectEmitter::DisableAllEmitters() {
		for (auto &e : sub_effectors) {
			e.flags &= ~tfxEmitterStateFlags_enabled;
			e.DisableAllEmitters();
		}
	}

	void EffectEmitter::DisableAllEmittersExcept(EffectEmitter &emitter) {
		for (auto &e : sub_effectors) {
			if (e.library_index == emitter.library_index)
				e.flags |= tfxEmitterStateFlags_enabled;
			else
				e.flags &= ~tfxEmitterStateFlags_enabled;
		}
	}

	Graph* EffectEmitter::GetGraphByType(GraphType type) {

		if (type < tfxGlobalCount) {
			return &((Graph*)&library->global_graphs[global])[type];
		}
		else if (type >= tfxPropertyStart && type < tfxBaseStart) {
			int ref = type - tfxPropertyStart;
			return &((Graph*)&library->property_graphs[property])[ref];
		}
		else if (type >= tfxBaseStart && type < tfxVariationStart) {
			int ref = type - tfxBaseStart;
			return &((Graph*)&library->base_graphs[base])[ref];
		}
		else if (type >= tfxVariationStart && type < tfxOvertimeStart) {
			int ref = type - tfxVariationStart;
			return &((Graph*)&library->variation_graphs[variation])[ref];
		}
		else if (type >= tfxOvertimeStart) {
			int ref = type - tfxOvertimeStart;
			return &((Graph*)&library->overtime_graphs[overtime])[ref];
		}

		return nullptr;

	}

	unsigned int EffectEmitter::GetGraphIndexByType(GraphType type) {

		if (type < tfxGlobalCount) {
			return global;
		}
		else if (type >= tfxPropertyStart && type < tfxBaseStart) {
			return property;
		}
		else if (type >= tfxBaseStart && type < tfxVariationStart) {
			return base;
		}
		else if (type >= tfxVariationStart && type < tfxOvertimeStart) {
			return variation;
		}
		else if (type >= tfxOvertimeStart) {
			return overtime;
		}

		assert(0);	//Unable to find a graph of that type (shouldn't happen)

		return 0;

	}

	void EffectEmitter::FreeGraphs() {
		if (type == tfxEffectType) {
			library->global_graphs[global].life.Free();
			library->global_graphs[global].amount.Free();
			library->global_graphs[global].velocity.Free();
			library->global_graphs[global].width.Free();
			library->global_graphs[global].height.Free();
			library->global_graphs[global].weight.Free();
			library->global_graphs[global].spin.Free();
			library->global_graphs[global].effect_angle.Free();
			library->global_graphs[global].stretch.Free();
			library->global_graphs[global].overal_scale.Free();
			library->global_graphs[global].intensity.Free();
			library->global_graphs[global].frame_rate.Free();
			library->global_graphs[global].splatter.Free();
		}

		if (type == tfxEmitterType) {
			library->property_graphs[property].emission_angle.Free();
			library->property_graphs[property].emission_range.Free();
			library->property_graphs[property].emitter_angle.Free();
			library->property_graphs[property].splatter.Free();
			library->property_graphs[property].emitter_width.Free();
			library->property_graphs[property].emitter_height.Free();
			library->property_graphs[property].arc_size.Free();
			library->property_graphs[property].arc_offset.Free();

			library->base_graphs[base].life.Free();
			library->base_graphs[base].amount.Free();
			library->base_graphs[base].velocity.Free();
			library->base_graphs[base].width.Free();
			library->base_graphs[base].height.Free();
			library->base_graphs[base].weight.Free();
			library->base_graphs[base].spin.Free();
			library->base_graphs[base].noise_offset.Free();

			library->variation_graphs[variation].life.Free();
			library->variation_graphs[variation].amount.Free();
			library->variation_graphs[variation].velocity.Free();
			library->variation_graphs[variation].width.Free();
			library->variation_graphs[variation].height.Free();
			library->variation_graphs[variation].weight.Free();
			library->variation_graphs[variation].spin.Free();
			library->variation_graphs[variation].noise_offset.Free();
			library->variation_graphs[variation].noise_resolution.Free();

			library->overtime_graphs[overtime].velocity.Free();
			library->overtime_graphs[overtime].width.Free();
			library->overtime_graphs[overtime].height.Free();
			library->overtime_graphs[overtime].weight.Free();
			library->overtime_graphs[overtime].spin.Free();
			library->overtime_graphs[overtime].stretch.Free();
			library->overtime_graphs[overtime].red.Free();
			library->overtime_graphs[overtime].green.Free();
			library->overtime_graphs[overtime].blue.Free();
			library->overtime_graphs[overtime].blendfactor.Free();
			library->overtime_graphs[overtime].intensity.Free();
			library->overtime_graphs[overtime].velocity_turbulance.Free();
			library->overtime_graphs[overtime].direction_turbulance.Free();
			library->overtime_graphs[overtime].velocity_adjuster.Free();
			library->overtime_graphs[overtime].direction.Free();
			library->overtime_graphs[overtime].noise_resolution.Free();
		}
	}

	void EffectEmitter::CompileGraphs() {
		if (type == tfxEffectType) {
			for (unsigned int t = (unsigned int)tfxGlobal_life; t != (unsigned int)tfxProperty_emission_angle; ++t) {
				CompileGraph(*GetGraphByType(GraphType(t)));
			}

			for (auto &emitter : sub_effectors) {
				for (unsigned int t = (unsigned int)tfxProperty_emission_angle; t != (unsigned int)tfxOvertime_velocity; ++t) {
					CompileGraph(*GetGraphByType((GraphType)t));
				}

				for (unsigned int t = (unsigned int)tfxOvertime_velocity; t != (unsigned int)tfxGraphMaxIndex; ++t) {
					CompileGraphOvertime(*emitter.GetGraphByType((GraphType)t));
				}
			}
		}
	}

	EffectEmitter& EffectLibrary::operator[] (uint32_t index) {
		return effects[index];
	}

	bool EffectLibrary::RenameEffect(EffectEmitter &effect, const char *new_name) {
		if (!NameExists(effect, new_name) && strlen(new_name) > 0) {
			effect.SetName(new_name);
			UpdateEffectPaths();
			return true;
		}

		return false;
	}

	bool EffectLibrary::NameExists(EffectEmitter& effect, const char *name) {
		for (auto &e : effects) {
			if (effect.library_index != e.library_index) {
				if (e.name == name) {
					return true;
				}
			}
		}

		return false;
	}

	bool EffectLibrary::NameExists2(EffectEmitter& effect, const char *name) {
		for (auto &e : effects) {
			if (e.name == name) {
				return true;
			}
		}
		return false;
	}

	void EffectLibrary::UpdateEffectPaths() {
		effect_paths.Clear();
		for (auto &e : effects) {
			tfxText path = e.name;
			e.path_hash = XXHash64::hash(path.c_str(), path.Length(), 0);
			AddPath(e, path);
		}
	}

	void EffectLibrary::AddPath(EffectEmitter &effectemitter, tfxText path) {
		effect_paths.Insert(path, &effectemitter);
		for (auto &sub : effectemitter.sub_effectors) {
			tfxText sub_path = path;
			sub_path.Appendf("/%s", sub.name.c_str());
			sub.path_hash = XXHash64::hash(sub_path.c_str(), sub_path.Length(), 0);
			AddPath(sub, sub_path);
		}
	}

	EffectEmitter &EffectLibrary::AddEffect(EffectEmitter &effect) {
		effect.library_index = effects.current_size;
		effect.type = tfxEffectType;
		effect.uid = ++uid;
		effect.library = this;
		effects.push_back(effect);
		ReIndex();
		UpdateEffectPaths();
		return effects.back();
	}

	EffectEmitter &EffectLibrary::AddFolder(tfxText name) {
		EffectEmitter folder;
		folder.name = name;
		folder.type = tfxFolder;
		folder.library = this;
		folder.uid = ++uid;
		effects.push_back(folder);
		ReIndex();
		UpdateEffectPaths();
		return effects.back();
	}

	EffectEmitter &EffectLibrary::AddFolder(EffectEmitter &folder) {
		assert(folder.type == tfxFolder);			//Must be type tfxFolder if adding a folder
		folder.library = this;
		folder.uid = ++uid;
		effects.push_back(folder);
		ReIndex();
		UpdateEffectPaths();
		return effects.back();
	}

	EffectEmitter* EffectLibrary::GetEffect(tfxText &path) {
		assert(effect_paths.ValidName(path));		//Effect was not found by that name
		return effect_paths.At(path);
	}

	EffectEmitter* EffectLibrary::GetEffect(const char *path) {
		assert(effect_paths.ValidName(path));		//Effect was not found by that name
		return effect_paths.At(path);
	}

	EffectEmitter* EffectLibrary::GetEffect(tfxKey key) {
		assert(effect_paths.ValidKey(key));			//Effect was not found by that key
		return effect_paths.At(key);
	}

	void EffectLibrary::PrepareEffectTemplate(tfxText path, tfxEffectTemplate &effect_template) {
		EffectEmitter *effect = GetEffect(path);
		assert(effect);
		assert(effect->type == tfxEffectType);
		effect->Clone(effect_template.effect_template, &effect_template.effect_template, this);
		effect_template.AddPath(effect_template.effect_template, effect_template.effect_template.name);
	}

	void EffectLibrary::ReIndex() {
		unsigned int index = 0;
		for (auto &e : effects) {
			e.library_index = index++;
			e.parent = nullptr;
			e.ReIndex();
		}
	}

	void EffectLibrary::UpdateParticleShapeReferences(tfxvec<EffectEmitter> &effects, unsigned int default_index) {
		for (auto &effect : effects) {
			if (effect.type == tfxFolder) {
				UpdateParticleShapeReferences(effect.sub_effectors, default_index);
			}
			else {
				for (auto &emitter : effect.sub_effectors) {
					if (particle_shapes.ValidIntName(emitter.properties.shape_index))
						emitter.properties.image = &particle_shapes.AtInt(emitter.properties.shape_index);
					else
						emitter.properties.image = &particle_shapes.AtInt(default_index);
					UpdateParticleShapeReferences(emitter.sub_effectors, default_index);
				}
			}
		}
	}

	EffectEmitter* EffectLibrary::MoveUp(EffectEmitter &effect) {
		if (effect.library_index > 0) {
			unsigned int new_index = effect.library_index - 1;
			std::swap(effects[effect.library_index], effects[new_index]);
			UpdateEffectPaths();
			ReIndex();
			return &effects[new_index];
		}
		return nullptr;
	}

	EffectEmitter* EffectLibrary::MoveDown(EffectEmitter &effect) {
		if (effect.library_index < effects.size() - 1) {
			unsigned int new_index = effect.library_index + 1;
			std::swap(effects[effect.library_index], effects[new_index]);
			UpdateEffectPaths();
			ReIndex();
			return &effects[new_index];
		}
		return nullptr;
	}

	void EffectLibrary::BuildComputeShapeData(void* dst, tfxVec4(uv_lookup)(void *ptr, ComputeImageData &image_data, int offset)) {
		assert(dst);	//must be a valid pointer to a space in memory
		assert(particle_shapes.Size());		//There are no shapes to copy!
		unsigned int index = 0;
		for (auto &shape : particle_shapes.data) {
			if (shape.animation_frames == 1) {
				ComputeImageData cs;
				cs.animation_frames = shape.animation_frames;
				cs.image_size = shape.image_size;
				cs.uv = uv_lookup(shape.ptr, cs, 0);
				shape_data.push_back(cs);
				shape.compute_shape_index = index++;
			}
			else {
				shape.compute_shape_index = index;
				for (int f = 0; f != shape.animation_frames; ++f) {
					ComputeImageData cs;
					cs.animation_frames = shape.animation_frames;
					cs.image_size = shape.image_size;
					cs.uv = uv_lookup(shape.ptr, cs, f);
					shape_data.push_back(cs);
					index++;
				}
			}
		}
	}

	void EffectLibrary::CopyComputeShapeData(void* dst) {
		assert(shape_data.size());	//You must call BuildComputeShapeData first
		memcpy(dst, shape_data.data, shape_data.size() * sizeof(ComputeImageData));
	}

	void EffectLibrary::CopyLookupIndexesData(void* dst) {
		assert(dst);	//must be a valid pointer to a space in memory
		assert(compiled_lookup_indexes.size());		//There is no data to copy, make sure a library has been loaded properly and it contains effects with emitters
		GraphLookupIndex *test = static_cast<GraphLookupIndex*>(dst);
		memcpy(dst, compiled_lookup_indexes.data, GetLookupIndexesSizeInBytes());
	}

	void EffectLibrary::CopyLookupValuesData(void* dst) {
		assert(dst);	//must be a valid pointer to a space in memory
		assert(compiled_lookup_indexes.size());		//There is no data to copy, make sure a library has been loaded properly and it contains effects with emitters
		memcpy(dst, compiled_lookup_values.data, GetLookupValuesSizeInBytes());
	}

	u32 EffectLibrary::GetComputeShapeDataSizeInBytes() {
		unsigned int frame_count = 0;
		for (auto &shape : particle_shapes.data) {
			frame_count += (unsigned int)shape.animation_frames;
		}
		return frame_count * sizeof(ComputeImageData);
	}

	u32 EffectLibrary::GetComputeShapeCount() {
		unsigned int frame_count = 0;
		for (auto &shape : particle_shapes.data) {
			frame_count += (unsigned int)shape.animation_frames;
		}
		return frame_count;
	}

	u32 EffectLibrary::GetLookupIndexCount() {
		return compiled_lookup_indexes.size() * tfxOvertimeCount;
	}

	u32 EffectLibrary::GetLookupValueCount() {
		return compiled_lookup_values.size();
	}

	u32 EffectLibrary::GetLookupIndexesSizeInBytes() {
		return sizeof(GraphLookupIndex) * tfxOvertimeCount * compiled_lookup_indexes.size();
	}

	u32 EffectLibrary::GetLookupValuesSizeInBytes() {
		return sizeof(float) * compiled_lookup_values.size();
	}

	void EffectLibrary::RemoveShape(unsigned int shape_index) {
		particle_shapes.RemoveInt(shape_index);
		for (auto &m : particle_shapes.map) {
			particle_shapes[m.index].shape_index = (unsigned int)m.key;
		}
	}

	void EffectLibrary::DeleteEffect(EffectEmitter *effect) {
		effects[effect->library_index].CleanUp();
		effects.erase(&effects[effect->library_index]);

		UpdateEffectPaths();
		ReIndex();
	}

	unsigned int EffectLibrary::AddGlobal() {
		if (free_global_graphs.size())
			return free_global_graphs.pop_back();
		GlobalAttributes global;
		global_graphs.push_back(global);
		return global_graphs.size() - 1;
	}
	unsigned int EffectLibrary::AddProperty() {
		if (free_property_graphs.size())
			return free_property_graphs.pop_back();
		PropertyAttributes property;
		property_graphs.push_back(property);
		return property_graphs.size() - 1;
	}
	unsigned int EffectLibrary::AddBase() {
		if (free_base_graphs.size())
			return free_base_graphs.pop_back();
		BaseAttributes base;
		base_graphs.push_back(base);
		return base_graphs.size() - 1;
	}
	unsigned int EffectLibrary::AddVariation() {
		if (free_variation_graphs.size())
			return free_variation_graphs.pop_back();
		VariationAttributes variation;
		variation_graphs.push_back(variation);
		return variation_graphs.size() - 1;
	}
	unsigned int EffectLibrary::AddOvertime() {
		if (free_overtime_graphs.size())
			return free_overtime_graphs.pop_back();
		OvertimeAttributes overtime;
		overtime_graphs.push_back(overtime);
		return overtime_graphs.size() - 1;
	}

	void EffectLibrary::FreeGlobal(unsigned int index) {
		assert(index < global_graphs.size());
		free_global_graphs.push_back(index);
	}
	void EffectLibrary::FreeProperty(unsigned int index) {
		assert(index < property_graphs.size());
		free_property_graphs.push_back(index);
	}
	void EffectLibrary::FreeBase(unsigned int index) {
		assert(index < base_graphs.size());
		free_base_graphs.push_back(index);
	}
	void EffectLibrary::FreeVariation(unsigned int index) {
		assert(index < variation_graphs.size());
		free_variation_graphs.push_back(index);
	}
	void EffectLibrary::FreeOvertime(unsigned int index) {
		assert(index < overtime_graphs.size());
		free_overtime_graphs.push_back(index);
	}

	unsigned int EffectLibrary::CloneGlobal(unsigned int source_index, EffectLibrary *destination_library) {
		unsigned int index = destination_library->AddGlobal();
		destination_library->global_graphs[index] = global_graphs[source_index];
		return index;
	}

	unsigned int EffectLibrary::CloneProperty(unsigned int source_index, EffectLibrary *destination_library) {
		unsigned int index = destination_library->AddProperty();
		destination_library->property_graphs[index] = property_graphs[source_index];
		return index;
	}

	unsigned int EffectLibrary::CloneBase(unsigned int source_index, EffectLibrary *destination_library) {
		unsigned int index = destination_library->AddBase();
		destination_library->base_graphs[index] = base_graphs[source_index];
		return index;
	}

	unsigned int EffectLibrary::CloneVariation(unsigned int source_index, EffectLibrary *destination_library) {
		unsigned int index = destination_library->AddVariation();
		destination_library->variation_graphs[index] = variation_graphs[source_index];
		return index;
	}

	unsigned int EffectLibrary::CloneOvertime(unsigned int source_index, EffectLibrary *destination_library) {
		unsigned int index = destination_library->AddOvertime();
		destination_library->overtime_graphs[index] = overtime_graphs[source_index];
		return index;
	}

	void EffectLibrary::AddEmitterGraphs(EffectEmitter& emitter) {
		emitter.property = AddProperty();
		emitter.base = AddBase();
		emitter.variation = AddVariation();
		emitter.overtime = AddOvertime();
	}

	void EffectLibrary::AddEffectGraphs(EffectEmitter& effect) {
		EffectEmitter *root_effect = effect.GetRootEffect();
		if (!root_effect)
			effect.global = AddGlobal();
		else
			effect.global = root_effect->global;
	}

	unsigned int EffectLibrary::AddAnimationSettings(EffectEmitter& effect) {
		assert(effect.type == tfxEffectType);
		AnimationSettings a;
		a.frames = 32;
		a.current_frame = 1;
		a.frame_offset = 0;
		a.position = tfxVec2(0.f, 0.f);
		a.frame_size = tfxVec2(256.f, 256.f);
		a.loop = false;
		a.seamless = false;
		a.seed = 0;
		a.zoom = 1.f;
		a.scale = 1.f;
		a.needs_recording = true;
		a.needs_exporting = 0;
		a.color_option = ExportColorOptions::tfxFullColor;
		a.export_option = ExportOptions::tfxSpriteSheet;
		a.export_with_transparency = true;
		animation_settings.push_back(a);
		effect.animation_settings = animation_settings.size() - 1;
		return effect.animation_settings;
	}

	void EffectLibrary::Clear() {
		for (auto &e : effects) {
			e.FreeGraphs();
		}
		effects.free_all();
		effect_paths.FreeAll();
		particle_shapes.FreeAll();
		global_graphs.free_all();
		property_graphs.free_all();
		base_graphs.free_all();
		variation_graphs.free_all();
		overtime_graphs.free_all();

		free_global_graphs.free_all();
		free_property_graphs.free_all();
		free_base_graphs.free_all();
		free_variation_graphs.free_all();
		free_overtime_graphs.free_all();
		uid = 0;
	}

	void EffectLibrary::UpdateEffectParticleStorage() {
		tfxParticleMemoryTools tools;
		for (auto &e : effects) {
			if (e.type != tfxFolder) {
				e.ResetAllBufferSizes();
				tools.ProcessEffect(e);
			}
			else {
				for (auto &sub : e.sub_effectors) {
					sub.ResetAllBufferSizes();
					tools.ProcessEffect(e);
				}
			}
		}
	}

	void EffectLibrary::UpdateComputeNodes() {
		unsigned int running_node_index = 0;
		unsigned int running_value_index = 0;
		tfxvec<EffectEmitter*> stack;
		all_nodes.clear();
		node_lookup_indexes.clear();
		compiled_lookup_values.clear();
		compiled_lookup_indexes.clear();
		for (auto &effect : effects) {
			stack.push_back(&effect);
			while(!stack.empty()) {
				EffectEmitter *current = stack.pop_back();
				if (current->type == tfxFolder) {
					for (auto &sub : current->sub_effectors) {
						stack.push_back(&sub);
					}
					continue;
				}
				EffectLookUpData lookup_data;
				EffectLookUpData value_lookup_data;
				memset(&lookup_data, 0, sizeof(EffectLookUpData));
				memset(&value_lookup_data, 0, sizeof(EffectLookUpData));
				if (current->type == tfxEmitterType) {

					int offset = tfxGlobalCount + tfxPropertyCount + tfxBaseCount + tfxVariationCount;

					current->lookup_value_index = compiled_lookup_indexes.size();
					for (int i = 0; i != tfxOvertimeCount; ++i) {
						Graph &graph = ((Graph*)(&overtime_graphs[current->overtime]))[i];
						GraphLookupIndex &index = ((GraphLookupIndex*)&lookup_data)[i];
						index.start_index = running_node_index;
						index.length = graph.nodes.size();
						index.max_life = graph.lookup.life;
						for (auto &node : graph.nodes) {
							all_nodes.push_back(node);
							running_node_index++;
						}

						GraphLookupIndex value_index;
						value_index.start_index = running_value_index;
						value_index.length = graph.lookup.values.size();
						value_index.max_life = graph.lookup.life;
						compiled_lookup_indexes.push_back(value_index);
						for (auto value : graph.lookup.values) {
							compiled_lookup_values.push_back(value);
							running_value_index++;
						}
					}

					node_lookup_indexes.push_back(lookup_data);
					current->lookup_node_index = node_lookup_indexes.size() - 1;

				}

				for (auto &sub : current->sub_effectors) {
					stack.push_back(&sub);
				}
			}
		}
	}

	void EffectLibrary::CompileAllGraphs() {
		for (auto &g : global_graphs) {
			CompileGraph(g.amount);
			CompileGraph(g.effect_angle);
			CompileGraph(g.frame_rate);
			CompileGraph(g.height);
			CompileGraph(g.width);
			CompileGraph(g.life);
			CompileGraph(g.intensity);
			CompileGraph(g.overal_scale);
			CompileGraph(g.spin);
			CompileGraph(g.splatter);
			CompileGraph(g.stretch);
			CompileGraph(g.velocity);
			CompileGraph(g.weight);
		}
		for (auto &g : property_graphs) {
			CompileGraph(g.arc_offset);
			CompileGraph(g.arc_size);
			CompileGraph(g.emission_angle);
			CompileGraph(g.emission_range);
			CompileGraph(g.emitter_angle);
			CompileGraph(g.emitter_width);
			CompileGraph(g.emitter_height);
			CompileGraph(g.splatter);
		}
		for (auto &g : base_graphs) {
			CompileGraph(g.amount);
			CompileGraph(g.width);
			CompileGraph(g.height);
			CompileGraph(g.life);
			CompileGraph(g.spin);
			CompileGraph(g.noise_offset);
			CompileGraph(g.velocity);
			CompileGraph(g.weight);
		}
		for (auto &g : variation_graphs) {
			CompileGraph(g.amount);
			CompileGraph(g.width);
			CompileGraph(g.height);
			CompileGraph(g.life);
			CompileGraph(g.noise_offset);
			CompileGraph(g.noise_resolution);
			CompileGraph(g.spin);
			CompileGraph(g.velocity);
			CompileGraph(g.weight);
		}
		for (auto &g : overtime_graphs) {
			CompileGraphOvertime(g.red);
			CompileGraphOvertime(g.green);
			CompileGraphOvertime(g.blue);
			CompileGraphOvertime(g.blendfactor);
			CompileGraphOvertime(g.intensity);
			CompileGraphOvertime(g.velocity_turbulance);
			CompileGraphOvertime(g.width);
			CompileGraphOvertime(g.height);
			CompileGraphOvertime(g.direction_turbulance);
			CompileGraphOvertime(g.spin);
			CompileGraphOvertime(g.stretch);
			CompileGraphOvertime(g.velocity);
			CompileGraph(g.velocity_adjuster);
			CompileGraphOvertime(g.weight);
			CompileGraphOvertime(g.direction);
			CompileGraphOvertime(g.noise_resolution);
		}
	}

	void EffectLibrary::CompileGlobalGraph(unsigned int index) {
		GlobalAttributes &g = global_graphs[index];
		CompileGraph(g.amount);
		CompileGraph(g.effect_angle);
		CompileGraph(g.frame_rate);
		CompileGraph(g.height);
		CompileGraph(g.width);
		CompileGraph(g.life);
		CompileGraph(g.intensity);
		CompileGraph(g.overal_scale);
		CompileGraph(g.spin);
		CompileGraph(g.splatter);
		CompileGraph(g.stretch);
		CompileGraph(g.velocity);
		CompileGraph(g.weight);
	}
	void EffectLibrary::CompilePropertyGraph(unsigned int index) {
		PropertyAttributes &g = property_graphs[index];
		CompileGraph(g.arc_offset);
		CompileGraph(g.arc_size);
		CompileGraph(g.emission_angle);
		CompileGraph(g.emission_range);
		CompileGraph(g.emitter_angle);
		CompileGraph(g.emitter_width);
		CompileGraph(g.emitter_height);
		CompileGraph(g.splatter);
	}
	void EffectLibrary::CompileBaseGraph(unsigned int index) {
		BaseAttributes &g = base_graphs[index];
		CompileGraph(g.amount);
		CompileGraph(g.width);
		CompileGraph(g.height);
		CompileGraph(g.life);
		CompileGraph(g.spin);
		CompileGraph(g.noise_offset);
		CompileGraph(g.velocity);
		CompileGraph(g.weight);
	}
	void EffectLibrary::CompileVariationGraph(unsigned int index) {
		VariationAttributes &g = variation_graphs[index];
		CompileGraph(g.amount);
		CompileGraph(g.width);
		CompileGraph(g.height);
		CompileGraph(g.life);
		CompileGraph(g.noise_offset);
		CompileGraph(g.noise_resolution);
		CompileGraph(g.spin);
		CompileGraph(g.velocity);
		CompileGraph(g.weight);
	}
	void EffectLibrary::CompileOvertimeGraph(unsigned int index) {
		OvertimeAttributes &g = overtime_graphs[index];
		CompileGraphOvertime(g.red);
		CompileGraphOvertime(g.green);
		CompileGraphOvertime(g.blue);
		CompileGraphOvertime(g.blendfactor);
		CompileGraphOvertime(g.intensity);
		CompileGraphOvertime(g.velocity_turbulance);
		CompileGraphOvertime(g.width);
		CompileGraphOvertime(g.height);
		CompileGraphOvertime(g.direction_turbulance);
		CompileGraphOvertime(g.spin);
		CompileGraphOvertime(g.stretch);
		CompileGraphOvertime(g.velocity);
		CompileGraph(g.velocity_adjuster);
		CompileGraphOvertime(g.weight);
		CompileGraphOvertime(g.direction);
		CompileGraphOvertime(g.noise_resolution);
	}
	void EffectLibrary::CompileColorGraphs(unsigned int index) {
		OvertimeAttributes &g = overtime_graphs[index];
		CompileGraphOvertime(g.red);
		CompileGraphOvertime(g.green);
		CompileGraphOvertime(g.blue);
	}

	void EffectLibrary::SetMinMaxData() {
		graph_min_max.clear();
		graph_min_max.create_pool(tfxGraphMaxIndex);

		graph_min_max[tfxGlobal_life] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_amount] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_velocity] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_width] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_height] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_weight] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_spin] = GetMinMaxGraphValues(tfxGlobalPercentPresetSigned);
		graph_min_max[tfxGlobal_stretch] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_overal_scale] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_intensity] = GetMinMaxGraphValues(tfxOpacityOvertimePreset);
		graph_min_max[tfxGlobal_frame_rate] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_splatter] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_effect_angle] = GetMinMaxGraphValues(tfxAnglePreset);

		graph_min_max[tfxProperty_emitter_angle] = GetMinMaxGraphValues(tfxAnglePreset);
		graph_min_max[tfxProperty_emission_angle] = GetMinMaxGraphValues(tfxAnglePreset);
		graph_min_max[tfxProperty_emission_range] = GetMinMaxGraphValues(tfxEmissionRangePreset);
		graph_min_max[tfxProperty_splatter] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxProperty_emitter_width] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxProperty_emitter_height] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxProperty_arc_size] = GetMinMaxGraphValues(tfxArcPreset);
		graph_min_max[tfxProperty_arc_offset] = GetMinMaxGraphValues(tfxArcPreset);

		graph_min_max[tfxBase_life] = GetMinMaxGraphValues(tfxLifePreset);
		graph_min_max[tfxBase_amount] = GetMinMaxGraphValues(tfxAmountPreset);
		graph_min_max[tfxBase_velocity] = GetMinMaxGraphValues(tfxVelocityPreset);
		graph_min_max[tfxBase_width] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxBase_height] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxBase_weight] = GetMinMaxGraphValues(tfxWeightPreset);
		graph_min_max[tfxBase_spin] = GetMinMaxGraphValues(tfxSpinPreset);
		graph_min_max[tfxBase_noise_offset] = GetMinMaxGraphValues(tfxDimensionsPreset);

		graph_min_max[tfxVariation_life] = GetMinMaxGraphValues(tfxLifePreset);
		graph_min_max[tfxVariation_amount] = GetMinMaxGraphValues(tfxAmountPreset);
		graph_min_max[tfxVariation_velocity] = GetMinMaxGraphValues(tfxVelocityPreset);
		graph_min_max[tfxVariation_width] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxVariation_height] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxVariation_weight] = GetMinMaxGraphValues(tfxWeightVariationPreset);
		graph_min_max[tfxVariation_spin] = GetMinMaxGraphValues(tfxSpinVariationPreset);
		graph_min_max[tfxVariation_noise_offset] = GetMinMaxGraphValues(tfxGlobalPercentPreset);

		graph_min_max[tfxOvertime_velocity] = GetMinMaxGraphValues(tfxVelocityOvertimePreset);
		graph_min_max[tfxOvertime_width] = GetMinMaxGraphValues(tfxPercentOvertime);
		graph_min_max[tfxOvertime_height] = GetMinMaxGraphValues(tfxPercentOvertime);
		graph_min_max[tfxOvertime_weight] = GetMinMaxGraphValues(tfxPercentOvertime);
		graph_min_max[tfxOvertime_spin] = GetMinMaxGraphValues(tfxSpinOvertimePreset);
		graph_min_max[tfxOvertime_stretch] = GetMinMaxGraphValues(tfxPercentOvertime);
		graph_min_max[tfxOvertime_red] = GetMinMaxGraphValues(tfxColorPreset);
		graph_min_max[tfxOvertime_green] = GetMinMaxGraphValues(tfxColorPreset);
		graph_min_max[tfxOvertime_blue] = GetMinMaxGraphValues(tfxColorPreset);
		graph_min_max[tfxOvertime_blendfactor] = GetMinMaxGraphValues(tfxOpacityOvertimePreset);
		graph_min_max[tfxOvertime_intensity] = GetMinMaxGraphValues(tfxIntensityOvertimePreset);
		graph_min_max[tfxOvertime_velocity_turbulance] = GetMinMaxGraphValues(tfxFrameratePreset);
		graph_min_max[tfxOvertime_direction_turbulance] = GetMinMaxGraphValues(tfxPercentOvertime);
		graph_min_max[tfxOvertime_velocity_adjuster] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxOvertime_direction] = GetMinMaxGraphValues(tfxDirectionOvertimePreset);

	}

	float EffectLibrary::LookupPreciseOvertimeNodeList(GraphType graph_type, int lookup_node_index, float age, float life) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		AttributeNode *lastec = nullptr;
		GraphLookupIndex &lookup_data = ((GraphLookupIndex*)&node_lookup_indexes[lookup_node_index])[graph_type];
		float min_y = graph_min_max[graph_type].y;
		float max_y = graph_min_max[graph_type].w;
		for (int i = lookup_data.start_index; i != lookup_data.start_index + lookup_data.length; ++i) {
			AttributeNode &a = all_nodes[i];
			float frame = a.frame * life;
			if (age < frame) {
				p = (age - lastf) / (frame - lastf);
				float bezier_value = GetBezierValue(lastec, a, p, min_y, max_y);
				if (bezier_value) {
					return bezier_value;
				}
				else {
					return lastv - p * (lastv - a.value);
				}
			}
			lastv = a.value;
			lastf = frame - 1;
			lastec = &a;
		}
		return lastv;
	}

	float EffectLibrary::LookupPreciseNodeList(GraphType graph_type, int lookup_node_index, float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		AttributeNode *lastec = nullptr;
		GraphLookupIndex &lookup_data = ((GraphLookupIndex*)&node_lookup_indexes[lookup_node_index])[graph_type];
		float min_y = graph_min_max[graph_type].y;
		float max_y = graph_min_max[graph_type].w;
		for (int i = lookup_data.start_index; i != lookup_data.start_index + lookup_data.length; ++i) {
			AttributeNode &a = all_nodes[i];
			if (age < a.frame) {
				p = (age - lastf) / (a.frame - lastf);
				float bezier_value = GetBezierValue(lastec, a, p, min_y, max_y);
				if (bezier_value) {
					return bezier_value;
				}
				else {
					return lastv - p * (lastv - a.value);
				}
			}
			lastv = a.value;
			lastf = a.frame - 1;
			lastec = &a;
		}
		return lastv;
	}

	float EffectLibrary::LookupFastValueList(GraphType graph_type, int lookup_node_index, float frame) {
		GraphLookupIndex &lookup_data = ((GraphLookupIndex*)&compiled_lookup_indexes[lookup_node_index])[graph_type];
		frame += lookup_data.start_index;
		unsigned int end_frame = lookup_data.start_index + lookup_data.length - 1;
		frame = frame > end_frame ? end_frame : frame;
		return compiled_lookup_values[(unsigned int)frame];
	}

	float EffectLibrary::LookupFastOvertimeValueList(GraphType graph_type, int lookup_value_index, float age, float lifetime) {
		GraphLookupIndex &lookup_data = ((GraphLookupIndex*)&compiled_lookup_indexes[lookup_value_index])[graph_type - tfxOvertime_velocity];
		float frame = (float)lookup_data.start_index;
		if (lifetime)
			frame += (age / lifetime * lookup_data.max_life) / tfxLOOKUP_FREQUENCY_OVERTIME;
		if (frame < lookup_data.start_index + lookup_data.length - 1)
			return compiled_lookup_values[(unsigned int)frame];
		return compiled_lookup_values[lookup_data.start_index + lookup_data.length - 1];
	}

	unsigned int EffectLibrary::CountOfGraphsInUse() {
		return global_graphs.size() + property_graphs.size() + base_graphs.size() + variation_graphs.size() + overtime_graphs.size() - CountOfFreeGraphs();
	}

	unsigned int EffectLibrary::CountOfFreeGraphs() {
		return free_global_graphs.size() + free_property_graphs.size() + free_base_graphs.size() + free_variation_graphs.size() + free_overtime_graphs.size();
	}

	DataTypesDictionary::DataTypesDictionary() {
		eff.Insert("name", tfxString);
		eff.Insert("image_index", tfxUint);
		eff.Insert("image_handle_x", tfxFloat);
		eff.Insert("image_handle_y", tfxFloat);
		eff.Insert("spawn_amount", tfxUint);
		eff.Insert("blend_mode", tfxSInt);
		eff.Insert("image_start_frame", tfxFloat);
		eff.Insert("image_end_frame", tfxFloat);
		eff.Insert("image_frame_rate", tfxFloat);

		eff.Insert("emission_type", tfxSInt);
		eff.Insert("emission_direction", tfxSInt);
		eff.Insert("grid_rows", tfxFloat);
		eff.Insert("grid_columns", tfxFloat);
		eff.Insert("loop_length", tfxFloat);
		eff.Insert("emitter_handle_x", tfxFloat);
		eff.Insert("emitter_handle_y", tfxFloat);
		eff.Insert("end_behaviour", tfxSInt);
		eff.Insert("angle_setting", tfxSInt);
		eff.Insert("angle_offset", tfxFloat);
		eff.Insert("multiply_blend_factor", tfxFloat);

		eff.Insert("random_color", tfxBool);
		eff.Insert("relative_position", tfxBool);
		eff.Insert("relative_angle", tfxBool);
		eff.Insert("image_handle_auto_center", tfxBool);
		eff.Insert("single", tfxBool);
		eff.Insert("one_shot", tfxBool);
		eff.Insert("spawn_on_grid", tfxBool);
		eff.Insert("grid_spawn_clockwise", tfxBool);
		eff.Insert("fill_area", tfxBool);
		eff.Insert("emitter_handle_auto_center", tfxBool);
		eff.Insert("edge_traversal", tfxBool);
		eff.Insert("image_reverse_animation", tfxBool);
		eff.Insert("image_play_once", tfxBool);
		eff.Insert("image_animate", tfxBool);
		eff.Insert("image_random_start_frame", tfxBool);
		eff.Insert("global_uniform_size", tfxBool);
		eff.Insert("base_uniform_size", tfxBool);
		eff.Insert("lifetime_uniform_size", tfxBool);
		eff.Insert("use_spawn_ratio", tfxBool);

		eff.Insert("frames", tfxUint);
		eff.Insert("current_frame", tfxUint);
		eff.Insert("frame_offset", tfxUint);
		eff.Insert("layer", tfxUint);
		eff.Insert("position_x", tfxFloat);
		eff.Insert("position_y", tfxFloat);
		eff.Insert("frame_width", tfxFloat);
		eff.Insert("frame_height", tfxFloat);
		eff.Insert("loop", tfxBool);
		eff.Insert("seamless", tfxBool);
		eff.Insert("seed", tfxUint);
		eff.Insert("zoom", tfxFloat);
		eff.Insert("scale", tfxFloat);
		eff.Insert("color_option", tfxSInt);
		eff.Insert("export_option", tfxSInt);
		eff.Insert("export_with_transparency", tfxBool);

		//Editor config, move this to the editor
		eff.Insert("only_play_selected_emitter", tfxBool);
		eff.Insert("load_examples", tfxBool);
		eff.Insert("load_last_file", tfxBool);
		eff.Insert("load_last_file_path", tfxString);
		eff.Insert("recent1", tfxString);
		eff.Insert("recent2", tfxString);
		eff.Insert("recent3", tfxString);
		eff.Insert("recent4", tfxString);
		eff.Insert("background_color_red", tfxFloat);
		eff.Insert("background_color_green", tfxFloat);
		eff.Insert("background_color_blue", tfxFloat);
		eff.Insert("use_checker_background", tfxBool);
		eff.Insert("preview_zoom", tfxFloat);
		eff.Insert("updates_per_second", tfxFloat);
		eff.Insert("background_image", tfxString);
		eff.Insert("use_background_image", tfxBool);
		eff.Insert("background_image_scale_x", tfxFloat);
		eff.Insert("background_image_scale_y", tfxFloat);
		eff.Insert("background_image_offset_x", tfxFloat);
		eff.Insert("background_image_offset_y", tfxFloat);
		eff.Insert("autoplay_effect", tfxSInt);
		eff.Insert("sync_refresh_rate", tfxBool);
		eff.Insert("window_maximised", tfxBool);
		eff.Insert("window_width", tfxSInt);
		eff.Insert("window_height", tfxSInt);
		eff.Insert("show_emitter_positions", tfxBool);
		eff.Insert("dpi_factor", tfxFloat);
		eff.Insert("graph_lookup_mode", tfxSInt);
		eff.Insert("show_tool_tips", tfxBool);
		eff.Insert("preview_trail_mode", tfxBool);
		eff.Insert("try_autorecover", tfxBool);
		eff.Insert("autorecovery_file", tfxString);
	}

	int ValidateEffectPackage(const char *filename) {
		tfxPackage package;
		int status = LoadPackage(filename, package);
		if (status) return status;					//returns 1 to 4 if it's an invalid package format

		tfxEntryInfo *data_txt = package.GetFile("data.txt");
		if (!data_txt) return -5;					//Unable to load the the data.txt file in the package

		return 0;
	}

	void AssignGraphData(EffectEmitter &effect, tfxvec<tfxText> &values) {
		if (values.size() > 0) {
			if (values[0] == "global_amount") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].amount.AddNode(n); }
			if (values[0] == "global_effect_angle") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].effect_angle.AddNode(n); }
			if (values[0] == "global_frame_rate") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].frame_rate.AddNode(n); }
			if (values[0] == "global_height") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].height.AddNode(n); }
			if (values[0] == "global_width") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].width.AddNode(n); }
			if (values[0] == "global_life") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].life.AddNode(n); }
			if (values[0] == "global_opacity") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].intensity.AddNode(n); }
			if (values[0] == "global_spin") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].spin.AddNode(n); }
			if (values[0] == "global_splatter") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].splatter.AddNode(n); }
			if (values[0] == "global_stretch") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].stretch.AddNode(n); }
			if (values[0] == "global_overal_scale") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].overal_scale.AddNode(n); }
			if (values[0] == "global_weight") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].weight.AddNode(n); }
			if (values[0] == "global_velocity") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].velocity.AddNode(n); }

			if (values[0] == "base_arc_offset") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].arc_offset.AddNode(n); }
			if (values[0] == "base_arc_size") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].arc_size.AddNode(n); }
			if (values[0] == "base_emission_angle") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emission_angle.AddNode(n); }
			if (values[0] == "base_emission_range") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emission_range.AddNode(n); }
			if (values[0] == "base_emitter_height") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emitter_height.AddNode(n); }
			if (values[0] == "base_emitter_width") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emitter_width.AddNode(n); }
			if (values[0] == "base_splatter") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].splatter.AddNode(n); }

			if (values[0] == "property_arc_offset") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].arc_offset.AddNode(n); }
			if (values[0] == "property_arc_size") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].arc_size.AddNode(n); }
			if (values[0] == "property_emitter_angle") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emitter_angle.AddNode(n); }
			if (values[0] == "property_emission_angle") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emission_angle.AddNode(n); }
			if (values[0] == "property_emission_range") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emission_range.AddNode(n); }
			if (values[0] == "property_emitter_height") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emitter_height.AddNode(n); }
			if (values[0] == "property_emitter_width") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emitter_width.AddNode(n); }
			if (values[0] == "property_splatter") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].splatter.AddNode(n); }

			if (values[0] == "base_amount") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].amount.AddNode(n); }
			if (values[0] == "base_life") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].life.AddNode(n); }
			if (values[0] == "base_height") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].height.AddNode(n); }
			if (values[0] == "base_width") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].width.AddNode(n); }
			if (values[0] == "base_spin") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].spin.AddNode(n); }
			if (values[0] == "base_noise_offset") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].noise_offset.AddNode(n); }
			if (values[0] == "base_velocity") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].velocity.AddNode(n); }
			if (values[0] == "base_weight") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].weight.AddNode(n); }

			if (values[0] == "variation_amount") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].amount.AddNode(n); }
			if (values[0] == "variation_height") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].height.AddNode(n); }
			if (values[0] == "variation_width") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].width.AddNode(n); }
			if (values[0] == "variation_life") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].life.AddNode(n); }
			if (values[0] == "variation_velocity") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].velocity.AddNode(n); }
			if (values[0] == "variation_weight") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].weight.AddNode(n); }
			if (values[0] == "variation_spin") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].spin.AddNode(n); }
			if (values[0] == "variation_motion_randomness") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].noise_offset.AddNode(n); }
			if (values[0] == "variation_noise_offset") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].noise_offset.AddNode(n); }
			if (values[0] == "variation_noise_resolution") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].noise_resolution.AddNode(n); }

			if (values[0] == "overtime_red") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].red.AddNode(n); }
			if (values[0] == "overtime_green") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].green.AddNode(n); }
			if (values[0] == "overtime_blue") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].blue.AddNode(n); }
			if (values[0] == "overtime_opacity") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].blendfactor.AddNode(n); }
			if (values[0] == "overtime_intensity") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].intensity.AddNode(n); }
			if (values[0] == "overtime_velocity_turbulance") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].velocity_turbulance.AddNode(n); }
			if (values[0] == "overtime_spin") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].spin.AddNode(n); }
			if (values[0] == "overtime_stretch") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].stretch.AddNode(n); }
			if (values[0] == "overtime_velocity") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].velocity.AddNode(n); }
			if (values[0] == "overtime_weight") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].weight.AddNode(n); }
			if (values[0] == "overtime_width") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].width.AddNode(n); }
			if (values[0] == "overtime_height") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].height.AddNode(n); }
			if (values[0] == "overtime_direction_turbulance") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].direction_turbulance.AddNode(n); }
			if (values[0] == "overtime_direction") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].direction.AddNode(n); }
			if (values[0] == "overtime_velocity_adjuster") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].velocity_adjuster.AddNode(n); }
			if (values[0] == "overtime_noise_resolution") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].noise_resolution.AddNode(n); }
		}
	}

	void AssignNodeData(AttributeNode &n, tfxvec<tfxText> &values) {
		n.frame = (float)atof(values[1].c_str());
		n.value = (float)atof(values[2].c_str());
		n.flags = (bool)atoi(values[3].c_str()) ? tfxAttributeNodeFlags_is_curve : 0;
		n.left.x = (float)atof(values[4].c_str());
		n.left.y = (float)atof(values[5].c_str());
		n.right.x = (float)atof(values[6].c_str());
		n.right.y = (float)atof(values[7].c_str());
		if (n.flags & tfxAttributeNodeFlags_is_curve)
			n.flags |= tfxAttributeNodeFlags_curves_initialised;
	}

	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, uint32_t value) {
		if (field == "image_index")
			effect.properties.shape_index = value;
		if (field == "spawn_amount")
			effect.properties.spawn_amount = value;
		if (field == "frames")
			effect.library->animation_settings[effect.animation_settings].frames = value;
		if (field == "current_frame")
			effect.library->animation_settings[effect.animation_settings].current_frame = value;
		if (field == "seed")
			effect.library->animation_settings[effect.animation_settings].seed = value;
		if (field == "layer")
			effect.properties.layer = value;
		if (field == "frame_offset")
			effect.library->animation_settings[effect.animation_settings].frame_offset = value;
	}
	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, int value) {
		if (field == "emission_type")
			effect.properties.emission_type = (EmissionType)value;
		if (field == "emission_direction")
			effect.properties.emission_direction = (EmissionDirection)value;
		if (field == "angle_setting")
			effect.properties.angle_setting = (AngleSetting)value;
		if (field == "color_option")
			effect.library->animation_settings[effect.animation_settings].color_option = (ExportColorOptions)value;
		if (field == "export_option")
			effect.library->animation_settings[effect.animation_settings].export_option = (ExportOptions)value;
		if (field == "end_behaviour")
			effect.properties.end_behaviour = (LineTraversalEndBehaviour)value;
		if (field == "frame_offset")
			effect.library->animation_settings[effect.animation_settings].frame_offset = value;
	}
	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, tfxText &value) {
		if (field == "name")
			effect.name = value;
	}
	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, float value) {
		if (field == "position_x")
			effect.library->animation_settings[effect.animation_settings].position.x = value;
		if (field == "position_y")
			effect.library->animation_settings[effect.animation_settings].position.y = value;
		if (field == "frame_width")
			effect.library->animation_settings[effect.animation_settings].frame_size.x = value;
		if (field == "frame_height")
			effect.library->animation_settings[effect.animation_settings].frame_size.y = value;
		if (field == "zoom")
			effect.library->animation_settings[effect.animation_settings].zoom = value;
		if (field == "scale")
			effect.library->animation_settings[effect.animation_settings].scale = value;
		if (field == "image_handle_x")
			effect.properties.image_handle.x = value;
		if (field == "image_handle_y")
			effect.properties.image_handle.y = value;
		if (field == "grid_rows")
			effect.properties.grid_points.x = value;
		if (field == "grid_columns")
			effect.properties.grid_points.y = value;
		if (field == "loop_length")
			effect.properties.loop_length = value;
		if (field == "emitter_handle_x")
			effect.properties.emitter_handle.x = value;
		if (field == "emitter_handle_y")
			effect.properties.emitter_handle.y = value;
		if (field == "image_start_frame")
			effect.properties.start_frame = value;
		if (field == "image_end_frame")
			effect.properties.end_frame = value;
		if (field == "image_frame_rate")
			effect.properties.frame_rate = value;
		if (field == "angle_offset")
			effect.properties.angle_offset = value;
	}
	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, bool value) {
		if (field == "loop")
			effect.library->animation_settings[effect.animation_settings].loop = value;
		if (field == "seamless")
			effect.library->animation_settings[effect.animation_settings].seamless = value;
		if (field == "export_with_transparency")
			effect.library->animation_settings[effect.animation_settings].export_with_transparency = value;
		if (field == "random_color")
			if (value) effect.properties.flags |= tfxEmitterPropertyFlags_random_color; else effect.properties.flags &= ~tfxEmitterPropertyFlags_random_color;
		if (field == "relative_position")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_relative_position; else effect.properties.flags &= ~tfxEmitterPropertyFlags_relative_position;
		if (field == "relative_angle")
			if (value) effect.properties.flags |= tfxEmitterPropertyFlags_relative_angle; else effect.properties.flags &= ~tfxEmitterPropertyFlags_relative_angle;
		if (field == "image_handle_auto_center")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_image_handle_auto_center; else effect.properties.flags &= ~tfxEmitterPropertyFlags_image_handle_auto_center;
		if (field == "single")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_single; else effect.properties.flags &= ~tfxEmitterPropertyFlags_single;
		if (field == "one_shot")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_one_shot; else effect.properties.flags &= ~tfxEmitterPropertyFlags_one_shot;
		if (field == "spawn_on_grid")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_spawn_on_grid; else effect.properties.flags &= ~tfxEmitterPropertyFlags_spawn_on_grid;
		if (field == "grid_spawn_clockwise")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_grid_spawn_clockwise; else effect.properties.flags &= ~tfxEmitterPropertyFlags_grid_spawn_clockwise;
		if (field == "fill_area")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_fill_area; else effect.properties.flags &= ~tfxEmitterPropertyFlags_fill_area;
		if (field == "emitter_handle_auto_center")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_emitter_handle_auto_center; else effect.properties.flags &= ~tfxEmitterPropertyFlags_emitter_handle_auto_center;
		if (field == "edge_traversal")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_edge_traversal; else effect.properties.flags &= ~tfxEmitterPropertyFlags_edge_traversal;
		if (field == "image_reverse_animation")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_reverse_animation; else effect.properties.flags &= ~tfxEmitterPropertyFlags_reverse_animation;
		if (field == "image_play_once")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_play_once; else effect.properties.flags &= ~tfxEmitterPropertyFlags_play_once;
		if (field == "image_animate")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_animate; else effect.properties.flags &= ~tfxEmitterPropertyFlags_animate;
		if (field == "image_random_start_frame")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_random_start_frame; else effect.properties.flags &= ~tfxEmitterPropertyFlags_random_start_frame;
		if (field == "global_uniform_size")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_global_uniform_size; else effect.properties.flags &= ~tfxEmitterPropertyFlags_global_uniform_size;
		if (field == "base_uniform_size")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_base_uniform_size; else effect.properties.flags &= ~tfxEmitterPropertyFlags_base_uniform_size;
		if (field == "lifetime_uniform_size")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_lifetime_uniform_size; else effect.properties.flags &= ~tfxEmitterPropertyFlags_lifetime_uniform_size;
		if (field == "use_spawn_ratio")
			if (value) effect.properties.flags |= tfxEmitterPropertyFlags_use_spawn_ratio; else effect.properties.flags &= ~tfxEmitterPropertyFlags_use_spawn_ratio;
	}

	void StreamProperties(EmitterProperties &property, tfxText &file) {

		file.AddLine("image_index=%i", property.shape_index);
		file.AddLine("image_handle_x=%f", property.image_handle.x);
		file.AddLine("image_handle_y=%f", property.image_handle.y);
		file.AddLine("image_start_frame=%f", property.start_frame);
		file.AddLine("image_end_frame=%f", property.end_frame);
		file.AddLine("image_frame_rate=%f", property.frame_rate);
		file.AddLine("image_play_once=%i", (property.flags & tfxEmitterPropertyFlags_play_once));
		file.AddLine("image_reverse_animation=%i", (property.flags & tfxEmitterPropertyFlags_reverse_animation));
		file.AddLine("image_animate=%i", (property.flags & tfxEmitterPropertyFlags_animate));
		file.AddLine("image_random_start_frame=%i", (property.flags & tfxEmitterPropertyFlags_random_start_frame));
		file.AddLine("spawn_amount=%i", property.spawn_amount);
		file.AddLine("emission_type=%i", property.emission_type);
		file.AddLine("emission_direction=%i", property.emission_direction);
		file.AddLine("grid_rows=%f", property.grid_points.x);
		file.AddLine("grid_columns=%f", property.grid_points.y);
		file.AddLine("loop_length=%f", property.loop_length);
		file.AddLine("emitter_handle_x=%f", property.emitter_handle.x);
		file.AddLine("emitter_handle_y=%f", property.emitter_handle.y);
		file.AddLine("end_behaviour=%i", property.end_behaviour);
		file.AddLine("random_color=%i", (property.flags & tfxEmitterPropertyFlags_random_color));
		file.AddLine("relative_position=%i", (property.flags & tfxEmitterPropertyFlags_relative_position));
		file.AddLine("relative_angle=%i", (property.flags & tfxEmitterPropertyFlags_relative_angle));
		file.AddLine("image_handle_auto_center=%i", (property.flags & tfxEmitterPropertyFlags_image_handle_auto_center));
		file.AddLine("single=%i", (property.flags & tfxEmitterPropertyFlags_single));
		file.AddLine("one_shot=%i", (property.flags & tfxEmitterPropertyFlags_one_shot));
		file.AddLine("spawn_on_grid=%i", (property.flags & tfxEmitterPropertyFlags_spawn_on_grid));
		file.AddLine("grid_spawn_clockwise=%i", (property.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise));
		file.AddLine("fill_area=%i", (property.flags & tfxEmitterPropertyFlags_fill_area));
		file.AddLine("emitter_handle_auto_center=%i", (property.flags & tfxEmitterPropertyFlags_emitter_handle_auto_center));
		file.AddLine("edge_traversal=%i", (property.flags & tfxEmitterPropertyFlags_edge_traversal));
		file.AddLine("angle_setting=%i", property.angle_setting);
		file.AddLine("angle_offset=%f", property.angle_offset);
		file.AddLine("global_uniform_size=%i", (property.flags & tfxEmitterPropertyFlags_global_uniform_size));
		file.AddLine("base_uniform_size=%i", (property.flags & tfxEmitterPropertyFlags_base_uniform_size));
		file.AddLine("lifetime_uniform_size=%i", (property.flags & tfxEmitterPropertyFlags_lifetime_uniform_size));
		file.AddLine("use_spawn_ratio=%i", (property.flags & tfxEmitterPropertyFlags_use_spawn_ratio));
		file.AddLine("layer=%i", property.layer);

	}

	void StreamGraph(const char * name, Graph &graph, tfxText &file) {

		for (auto &n : graph.nodes) {
			file.AddLine("%s,%f,%f,%i,%f,%f,%f,%f", name, n.frame, n.value, (n.flags & tfxAttributeNodeFlags_is_curve), n.left.x, n.left.y, n.right.x, n.right.y);
		}

	}

	Random::Random() {
		ReSeed();
	}

	void Random::ReSeed() {
		seeds[0] = (uint64_t)Millisecs();
		seeds[1] = (uint64_t)Millisecs() * 2;
	}

	void Random::ReSeed(uint64_t seed1, uint64_t seed2) {
		seeds[0] = seed1;
		seeds[1] = seed2;
	}

	double Random::Millisecs() {
		auto now = Clock::now().time_since_epoch();
		auto m = std::chrono::duration_cast<std::chrono::milliseconds>(now).count();
		return double(m);
	}

	static bool CompareNodes(AttributeNode &left, AttributeNode &right) {
		return left.frame < right.frame;
	}

	bool Graph::IsOvertimeGraph() {
		return type >= tfxOvertime_velocity && type != tfxOvertime_velocity_adjuster;
	}

	bool Graph::IsGlobalGraph() {
		return type >= tfxGlobal_life && type <= tfxGlobal_splatter;
	}

	bool Graph::IsAngleGraph() {
		return (type == tfxGlobal_effect_angle || type == tfxProperty_emission_angle || type == tfxProperty_emission_range || type == tfxProperty_emitter_angle ||
			type == tfxProperty_arc_offset || type == tfxProperty_arc_size || type == tfxBase_spin || type == tfxVariation_spin || type == tfxOvertime_direction);
	}

	float AttributeNode::GetX() {
		return frame;
	}
	float AttributeNode::GetY() {
		return value;
	}

	bool SetNode(Graph &graph, AttributeNode &node, float _frame, float _value, tfxAttributeNodeFlags flags, float _c0x, float _c0y, float _c1x, float _c1y) {
		node.frame = _frame;
		node.value = _value;
		node.flags = flags;
		node.left.x = _c0x;
		node.left.y = _c0y;
		node.right.x = _c1x;
		node.right.y = _c1y;
		if (graph.nodes[0] == node) {
			node.frame = graph.min.x;
			ClampNode(graph, node);
		}
		else {
			ClampNode(graph, node);
		}

		if (graph.Sort()) {
			graph.ReIndex();
			return true;
		}

		return false;
	}

	bool SetNodeFrame(Graph &graph, AttributeNode &node, float &_frame) {
		node.frame = _frame;
		if (graph.nodes[0] == node) {
			node.frame = graph.min.x;
			ClampNode(graph, node);
		}
		else {
			ClampNode(graph, node);
		}
		_frame = node.frame;

		if (graph.Sort()) {
			graph.ReIndex();
			return true;
		}

		return false;

	}

	bool SetNodeValue(Graph &graph, AttributeNode &node, float &_value) {
		node.value = _value;
		ClampNode(graph, node);
		_value = node.value;

		return false;
	}

	bool SetNode(Graph &graph, AttributeNode &node, float &frame, float &value) {
		float old_frame = node.frame;
		float old_value = node.value;

		node.frame = frame;
		node.value = value;

		if (graph.nodes[0] == node) {
			node.frame = graph.min.x;
			ClampNode(graph, node);
		}
		else {
			ClampNode(graph, node);
		}

		if (node.flags & tfxAttributeNodeFlags_curves_initialised) {
			node.left.y += node.value - old_value;
			node.left.x += node.frame - old_frame;
			node.right.y += node.value - old_value;
			node.right.x += node.frame - old_frame;
		}

		frame = node.frame;
		value = node.value;

		if (graph.Sort()) {
			graph.ReIndex();
			return true;
		}

		return false;
	}

	void SetCurve(Graph &graph, AttributeNode &node, bool is_left_curve, float &frame, float &value) {
		if (is_left_curve) {
			node.left.x = frame;
			node.left.y = value;
			if (node.left.x > node.frame)
				node.left.x = node.frame;
			else
				ClampCurve(graph, node.left, node);
			frame = node.left.x;
			value = node.left.y;
		}
		else {
			node.right.x = frame;
			node.right.y = value;
			if (node.right.x < node.frame)
				node.right.x = node.frame;
			else
				ClampCurve(graph, node.right, node);
			frame = node.right.x;
			value = node.right.y;
		}
	}

	bool MoveNode(Graph &graph, AttributeNode &node, float frame, float value, bool sort) {
		float old_frame = node.frame;
		float old_value = node.value;

		node.frame += frame;
		node.value += value;

		if (graph.nodes[0] == node) {
			node.frame = graph.min.x;
			ClampNode(graph, node);
		}
		else {
			ClampNode(graph, node);
		}

		if (node.flags & tfxAttributeNodeFlags_curves_initialised) {
			node.left.y += node.value - old_value;
			node.left.x += node.frame - old_frame;
			node.right.y += node.value - old_value;
			node.right.x += node.frame - old_frame;
		}

		if (sort) {
			if (graph.Sort()) {
				graph.ReIndex();
				return true;
			}
		}

		return false;
	}

	void ClampNode(Graph &graph, AttributeNode &node) {
		if (node.value < graph.min.y) node.value = graph.min.y;
		if (node.frame < graph.min.x) node.frame = graph.min.x;
		if (node.value > graph.max.y) node.value = graph.max.y;
		if (node.frame > graph.max.x) node.frame = graph.max.x;
	}

	void ClampGraph(Graph &graph) {
		for (auto &node : graph.nodes) {
			ClampNode(graph, node);
			if (node.flags & tfxAttributeNodeFlags_is_curve) {
				ClampCurve(graph, node.left, node);
				ClampCurve(graph, node.right, node);
			}
		}
	}

	void ClampCurve(Graph &graph, Point &p, AttributeNode &node) {
		if (p.y < graph.min.y) p.y = graph.min.y;
		if (p.x < graph.min.x) p.x = graph.min.x;
		//if (p.y > graph.max.y) p.y = graph.max.y;
		if (p.x > graph.max.x) p.x = graph.max.x;

		AttributeNode *next = graph.GetNextNode(node);
		if (next) {
			if (p.x > next->frame) p.x = next->frame;
		}

		AttributeNode *prev = graph.GetPrevNode(node);
		if (prev) {
			if (p.x < prev->frame) p.x = prev->frame;
		}
	}

	Graph::Graph() {
		min.x = 0.f;
		min.y = 0.f;
		max.x = 1000.f;
		max.y = 1000.f;

		effector = nullptr;
	}

	Graph::~Graph() {
		Free();
	}

	const float BEZIER_ACCURACY = 0.01f;

	float GetDistance(float fromx, float fromy, float tox, float toy) {

		float w = tox - fromx;
		float h = toy - fromy;

		return std::sqrt(w * w + h * h);

	}

	float GetVectorAngle(float x, float y) {
		return std::atan2f(x, -y);
	}

	Point GetQuadBezier(Point p0, Point p1, Point p2, float t, float ymin, float ymax, bool clamp) {
		Point b;
		b.x = std::powf(1.f - t, 2.f) * p0.x + 2.f * t * (1.f - t) * p1.x + std::powf(t, 2.f) * p2.x;
		b.y = std::powf(1.f - t, 2.f) * p0.y + 2.f * t * (1.f - t) * p1.y + std::powf(t, 2.f) * p2.y;
		if (b.x < p0.x) b.x = p0.x;
		if (b.x > p2.x) b.x = p2.x;
		if (clamp) {
			if (b.y < ymin) b.y = ymin;
			if (b.y > ymax) b.y = ymax;
		}
		return b;
	}

	Point GetCubicBezier(Point p0, Point p1, Point p2, Point p3, float t, float ymin, float ymax, bool clamp) {
		Point b;
		b.x = std::powf(1.f - t, 3.f) * p0.x + 3.f * t * std::powf(1.f - t, 2.f) * p1.x + 3.f * std::powf(t, 2.f) * (1.f - t) * p2.x + std::powf(t, 3.f) * p3.x;
		b.y = std::powf(1.f - t, 3.f) * p0.y + 3.f * t * std::powf(1.f - t, 2.f) * p1.y + 3.f * std::powf(t, 2.f) * (1.f - t) * p2.y + std::powf(t, 3.f) * p3.y;
		if (b.x < p0.x) b.x = p0.x;
		if (b.x > p3.x) b.x = p3.x;
		if (clamp) {
			if (b.y < ymin) b.y = ymin;
			if (b.y > ymax) b.y = ymax;
		}
		return b;
	}

	float GetBezierValue(const AttributeNode *lastec, const AttributeNode &a, float t, float ymin, float ymax) {
		if (lastec) {
			if (a.flags & tfxAttributeNodeFlags_is_curve) {
				if (lastec->flags & tfxAttributeNodeFlags_is_curve) {
					Point p0(lastec->frame, lastec->value);
					Point p1(lastec->right.x, lastec->right.y);
					Point p2(a.left.x, a.left.y);
					Point p3(a.frame, a.value);
					Point value = GetCubicBezier(p0, p1, p2, p3, t, ymin, ymax);
					return value.y;
				}
				else {
					Point p0(lastec->frame, lastec->value);
					Point p1(a.left.x, a.left.y);
					Point p2(a.frame, a.value);
					Point value = GetQuadBezier(p0, p1, p2, t, ymin, ymax);
					return value.y;
				}
			}
			else if (lastec->flags & tfxAttributeNodeFlags_is_curve) {
				Point p0(lastec->frame, lastec->value);
				Point p1(lastec->right.x, lastec->right.y);
				Point p2(a.frame, a.value);
				Point value = GetQuadBezier(p0, p1, p2, t, ymin, ymax);
				return value.y;
			}
		}
		else {
			return 0;
		}

		return 0;
	}

	AttributeNode* Graph::AddNode(float _frame, float _value, tfxAttributeNodeFlags flags, float _c0x, float _c0y, float _c1x, float _c1y) {
		AttributeNode node;

		if (nodes.size())
			node.frame = _frame;
		else
			node.frame = 0.f;

		node.value = _value;
		node.flags = flags;
		node.left.x = _c0x;
		node.left.y = _c0y;
		node.right.x = _c1x;
		node.right.y = _c1y;
		ClampNode(*this, node);
		nodes.push_back(node);
		Sort();

		ReIndex();
		return &nodes.back();
	}

	void Graph::AddNode(AttributeNode &node) {
		for (auto &n : nodes) {
			if (n.frame == node.frame)
				return;
		}
		nodes.push_back(node);
		Sort();
		ReIndex();
	}

	AttributeNode* Graph::AddCoordNode(float _frame, float _value) {
		AttributeNode node;

		if (nodes.size())
			node.frame = _frame;
		else
			node.frame = 0.f;

		node.value = _value;
		node.flags = 0;
		node.left.x = 0.f;
		node.left.y = 0.f;
		node.right.x = 0.f;
		node.right.y = 0.f;
		ClampNode(*this, node);
		AttributeNode &n = nodes.push_back(node);
		if (Sort()) {
			ReIndex();
			return nodes.find(n);
		}

		ReIndex();
		return &n;
	}

	AttributeNode* Graph::InsertCoordNode(float _frame, float _value) {
		AttributeNode node;

		if (nodes.size())
			node.frame = _frame;
		else
			node.frame = 0.f;

		node.value = _value;
		node.flags = 0;
		node.left.x = 0.f;
		node.left.y = 0.f;
		node.right.x = 0.f;
		node.right.y = 0.f;
		ClampNode(*this, node);

		if (nodes.size() > 1) {
			AttributeNode *last_node = nullptr;
			for (auto *n = nodes.begin() + 1; n != nodes.end(); ++n) {
				if (node.frame < n->frame)
					last_node = n;
				else
					break;
			}

			if (last_node) {
				AttributeNode *r_value = nodes.insert(last_node, node);
				ReIndex();
				return r_value;
			}
		}

		AttributeNode *r_value = &nodes.push_back(node);
		ReIndex();
		return r_value;
	}

	AttributeNode* Graph::InsertNode(float _frame, float _value) {
		AttributeNode node;

		if (nodes.size())
			node.frame = _frame;
		else
			node.frame = 0.f;

		node.value = _value;
		node.flags = 0;
		node.left.x = 0.f;
		node.left.y = 0.f;
		node.right.x = 0.f;
		node.right.y = 0.f;
		ClampNode(*this, node);

		if (nodes.size() > 1) {
			AttributeNode *last_node = nullptr;
			for (auto *n = nodes.begin() + 1; n != nodes.end(); ++n) {
				if (node.frame < n->frame)
					last_node = n;
				else
					break;
			}

			if (last_node) {
				AttributeNode *r_value = nodes.insert(last_node, node);
				ReIndex();
				return r_value;
			}
		}

		AttributeNode *r_value = &nodes.push_back(node);
		ReIndex();
		return r_value;
	}

	void Graph::SetNode(uint32_t i, float _frame, float _value, tfxAttributeNodeFlags flags, float _c0x, float _c0y, float _c1x, float _c1y) {
		if (!nodes.empty() && i < nodes.size()) {
			nodes[i].frame = _frame;
			nodes[i].value = _value;
			nodes[i].flags = flags;
			nodes[i].left.x = _c0x;
			nodes[i].left.y = _c0y;
			nodes[i].right.x = _c1x;
			nodes[i].right.y = _c1y;
			if (Sort())
				ReIndex();
		}
	}

	tfxvec<AttributeNode>& Graph::Nodes() {
		return nodes;
	}

	float Graph::GetValue(float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		AttributeNode *lastec = nullptr;
		for (auto &a : nodes) {
			if (age < a.frame) {
				p = (age - lastf) / (a.frame - lastf);
				float bezier_value = GetBezierValue(lastec, a, p, min.y, max.y);
				if (bezier_value) {
					return bezier_value;
				}
				else {
					return lastv - p * (lastv - a.value);
				}
			}
			lastv = a.value;
			lastf = a.frame - 1;
			lastec = &a;
		}
		return lastv;

	}

	AttributeNode *Graph::GetNextNode(AttributeNode &node) {
		if (node.index < nodes.size() - 1) {
			return &nodes[node.index + 1];
		}

		return nullptr;
	}

	AttributeNode *Graph::GetPrevNode(AttributeNode &node) {
		if (node.index > 0) {
			return &nodes[node.index - 1];
		}

		return nullptr;
	}

	AttributeNode *Graph::GetLastNode() {
		return &nodes.back();
	}

	float Graph::GetRandomValue(float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		AttributeNode *lastec = nullptr;
		for (auto &a : nodes) {
			if (age < a.frame) {
				p = (age - lastf) / (a.frame - lastf);
				float bezier_value = GetBezierValue(lastec, a, p, min.y, max.y);
				if (bezier_value) {
					return random_generation.Range(bezier_value);
				}
				else {
					return random_generation.Range(lastv - p * (lastv - a.value));
				}
			}
			lastv = a.value;
			lastf = a.frame - 1;
			lastec = &a;
		}
		return random_generation.Range(lastv);

	}

	float Graph::GetValue(float age, float life) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		AttributeNode *lastec = nullptr;
		for (auto &a : nodes) {
			float frame = a.frame * life;
			if (age < frame) {
				p = (age - lastf) / (frame - lastf);
				float bezier_value = GetBezierValue(lastec, a, p, min.y, max.y);
				if (bezier_value) {
					return bezier_value;
				}
				else {
					return lastv - p * (lastv - a.value);
				}
			}
			lastv = a.value;
			lastf = frame - 1;
			lastec = &a;
		}
		return lastv;
	}

	float Graph::GetFirstValue() {
		if (nodes.size())
			return nodes.front().value;
		return 0.f;
	}

	float* Graph::LinkFirstValue() {
		if (nodes.size())
			return &nodes.front().value;
		return nullptr;
	}

	float Graph::GetLastValue() {
		if (nodes.size())
			return nodes.back().value;

		return 0.f;
	}

	float Graph::GetMaxValue() {
		if (nodes.size()) {
			float value = tfxMIN_FLOAT;
			for (auto &n : nodes) {
				if (value < n.value)
					value = n.value;
			}
			return value;
		}
		return 0.f;
	}

	float Graph::GetMinValue() {
		if (nodes.size()) {
			float value = tfxMAX_FLOAT;
			for (auto &n : nodes) {
				if (value > n.value)
					value = n.value;
			}
			return value;
		}
		return 0.f;
	}

	float Graph::GetLastFrame() {
		if (nodes.size())
			return nodes.back().frame;

		return 0.f;
	}

	AttributeNode* Graph::FindNode(const AttributeNode &n) {
		return nodes.find(n);
	}

	void Graph::ValidateCurves() {
		unsigned int index = 0;
		unsigned int last_index = nodes.size() - 1;
		for(auto &n : nodes) {
			if (n.flags & tfxAttributeNodeFlags_is_curve) {
				if (index < last_index) {
					if (nodes[index + 1].frame < n.right.x)
						n.right.x = nodes[index + 1].frame;
				}
				if (index > 0) {
					if (nodes[index - 1].frame > n.left.x)
						n.left.x = nodes[index - 1].frame;
				}
				if (n.left.x > n.frame)
					n.left.x = n.frame;
				if (n.right.x < n.frame)
					n.right.x = n.frame;
			}
			index++;
		}
	}

	void Graph::DeleteNode(const AttributeNode &n) {
		nodes.erase(&n);
	}

	void Graph::Reset(float v, GraphPreset preset, bool add_node) {
		nodes.clear();
		if (add_node)
			AddNode(0.f, v);
		switch (preset) {
		case GraphPreset::tfxGlobalPercentPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 20.f };
			break;
		case GraphPreset::tfxGlobalOpacityPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 1.f };
			break;
		case GraphPreset::tfxGlobalPercentPresetSigned:
			min = { 0.f, -20.f }; max = { tfxMAX_FRAME, 20.f };
			break;
		case GraphPreset::tfxAnglePreset:
			min = { 0.f, -1080.f }; max = { tfxMAX_FRAME, 1080.f };
			break;
		case GraphPreset::tfxArcPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 360.f };
			break;
		case GraphPreset::tfxEmissionRangePreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 360.f };
			break;
		case GraphPreset::tfxDimensionsPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 4000.f };
			break;
		case GraphPreset::tfxLifePreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 100000.f };
			break;
		case GraphPreset::tfxAmountPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 5000.f };
			break;
		case GraphPreset::tfxVelocityPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 10000.f };
			break;
		case GraphPreset::tfxWeightPreset:
			min = { 0.f, -2500.f }; max = { tfxMAX_FRAME, 2500.f };
			break;
		case GraphPreset::tfxWeightVariationPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 2500.f };
			break;
		case GraphPreset::tfxNoiseOffsetVariationPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 1000.f };
			break;
		case GraphPreset::tfxNoiseResolutionPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 10000.f };
			break;
		case GraphPreset::tfxSpinPreset:
			min = { 0.f, -2000.f }; max = { tfxMAX_FRAME, 2000.f };
			break;
		case GraphPreset::tfxSpinVariationPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 2000.f };
			break;
		case GraphPreset::tfxDirectionVariationPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 22.5f };
			break;
		case GraphPreset::tfxWeightOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 20.f };
			break;
		case GraphPreset::tfxDirectionOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 4320.f };
			break;
		case GraphPreset::tfxSpinOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 20.f };
			break;
		case GraphPreset::tfxVelocityOvertimePreset:
			min = { 0.f, -20.f }; max = { 1.f, 20.f };
			break;
		case GraphPreset::tfxPercentOvertime:
			min = { 0.f, 0.f }; max = { 1.f, 20.f };
			break;
		case GraphPreset::tfxFrameratePreset:
			min = { 0.f, 0.f }; max = { 1.f, 200.f };
			break;
		case GraphPreset::tfxOpacityOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 1.f };
			break;
		case GraphPreset::tfxColorPreset:
			min = { 0.f, 0.f }; max = { 1.f, 255.f };
			break;
		case GraphPreset::tfxIntensityOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 5.f };
			break;
		}

		graph_preset = preset;
	}

	tfxVec4 GetMinMaxGraphValues(GraphPreset preset) {
		tfxVec4 mm;
		switch (preset) {
		case GraphPreset::tfxGlobalPercentPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 20.f };
			break;
		case GraphPreset::tfxGlobalOpacityPreset:
			mm = { 0.f, 0.f , tfxMAX_FRAME, 1.f };
			break;
		case GraphPreset::tfxGlobalPercentPresetSigned:
			mm = { 0.f, -20.f, tfxMAX_FRAME, 20.f };
			break;
		case GraphPreset::tfxAnglePreset:
			mm = { 0.f, -1080.f, tfxMAX_FRAME, 1080.f };
			break;
		case GraphPreset::tfxArcPreset:
			mm = { 0.f, 0.f , tfxMAX_FRAME, 360.f };
			break;
		case GraphPreset::tfxEmissionRangePreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 360.f };
			break;
		case GraphPreset::tfxDimensionsPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 4000.f };
			break;
		case GraphPreset::tfxLifePreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 100000.f };
			break;
		case GraphPreset::tfxAmountPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 5000.f };
			break;
		case GraphPreset::tfxVelocityPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 10000.f };
			break;
		case GraphPreset::tfxWeightPreset:
			mm = { 0.f, -2500.f, tfxMAX_FRAME, 2500.f };
			break;
		case GraphPreset::tfxWeightVariationPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 2500.f };
			break;
		case GraphPreset::tfxSpinPreset:
			mm = { 0.f, -2000.f, tfxMAX_FRAME, 2000.f };
			break;
		case GraphPreset::tfxSpinVariationPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 2000.f };
			break;
		case GraphPreset::tfxDirectionVariationPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 22.5f };
			break;
		case GraphPreset::tfxWeightOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 20.f };
			break;
		case GraphPreset::tfxDirectionOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 4320.f };
			break;
		case GraphPreset::tfxSpinOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 20.f };
			break;
		case GraphPreset::tfxVelocityOvertimePreset:
			mm = { 0.f, -20.f, 1.f, 20.f };
			break;
		case GraphPreset::tfxPercentOvertime:
			mm = { 0.f, 0.f, 1.f, 20.f };
			break;
		case GraphPreset::tfxFrameratePreset:
			mm = { 0.f, 0.f, 1.f, 200.f };
			break;
		case GraphPreset::tfxOpacityOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 1.f };
			break;
		case GraphPreset::tfxColorPreset:
			mm = { 0.f, 0.f, 1.f, 255.f };
			break;
		case GraphPreset::tfxIntensityOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 5.f };
			break;
		}
		
		return mm;
	}

	void Graph::DragValues(GraphPreset preset, float &frame, float &value) {
		switch (preset) {
		case GraphPreset::tfxOpacityOvertimePreset:
		case GraphPreset::tfxGlobalPercentPreset:
		case GraphPreset::tfxIntensityOvertimePreset:
			frame = 0.001f;
			value = 0.001f;
			break;
		case GraphPreset::tfxDirectionOvertimePreset:
			frame = 0.001f;
			value = 0.1f;
			break;
		case GraphPreset::tfxLifePreset:
			frame = 5;
			value = 5;
			break;
		case GraphPreset::tfxAnglePreset:
		case GraphPreset::tfxArcPreset:
		case GraphPreset::tfxEmissionRangePreset:
			frame = 5;
			value = 0.1f;
			break;
		case GraphPreset::tfxDimensionsPreset:
		case GraphPreset::tfxAmountPreset:
		case GraphPreset::tfxVelocityPreset:
		case GraphPreset::tfxWeightPreset:
		case GraphPreset::tfxWeightVariationPreset:
		case GraphPreset::tfxSpinPreset:
		case GraphPreset::tfxSpinVariationPreset:
		case GraphPreset::tfxFrameratePreset:
			frame = 5.f;
			value = 1.f;
			break;
		case GraphPreset::tfxWeightOvertimePreset:
		case GraphPreset::tfxVelocityOvertimePreset:
		case GraphPreset::tfxSpinOvertimePreset:
		case GraphPreset::tfxDirectionVariationPreset:
			frame = 0.001f;
			value = 0.01f;
			break;
		case GraphPreset::tfxColorPreset:
			frame = 0.001f;
			value = 1.f;
			break;
		case GraphPreset::tfxPercentOvertime:
			frame = 0.05f;
			value = 0.05f;
			break;
		default:
			frame = 1;
			value = 0.1f;
			break;
		}
	}

	void Graph::Clear() {
		nodes.clear();
	}

	void Graph::Free() {
		nodes.free_all();
	}

	void Graph::Copy(Graph &to) {
		to.nodes.reserve(nodes.size());
		std::copy(nodes.begin(), nodes.end(), to.nodes.begin());
		to.nodes.current_size = nodes.current_size;
	}

	bool Graph::Sort() {
		if (!std::is_sorted(nodes.begin(), nodes.end(), CompareNodes)) {
			std::sort(nodes.begin(), nodes.end(), CompareNodes);
			return true;
		}
		return false;
	}

	void Graph::ReIndex() {
		unsigned int i = 0;
		for (auto &a : nodes) {
			a.index = i++;
		}
	}

	tfxVec2 Graph::GetInitialZoom() {
		switch (graph_preset) {
		case GraphPreset::tfxOpacityOvertimePreset:
			return tfxVec2(.0017f, 0.00275f);
		case GraphPreset::tfxGlobalPercentPreset:
			return tfxVec2(10.f, 0.005f);
			break;
		case GraphPreset::tfxGlobalPercentPresetSigned:
			return tfxVec2(10.f, 0.006f);
			break;
		case GraphPreset::tfxGlobalOpacityPreset:
			return tfxVec2(10.f, 0.003f);
			break;
		case GraphPreset::tfxLifePreset:
			return tfxVec2(10.f, 3.5f);
			break;
		case GraphPreset::tfxAnglePreset:
			return tfxVec2(10.f, 1.f);
			break;
		case GraphPreset::tfxArcPreset:
			return tfxVec2(10.f, 1.f);
			break;
		case GraphPreset::tfxEmissionRangePreset:
			return tfxVec2(10.f, .5f);
			break;
		case GraphPreset::tfxAmountPreset:
			return tfxVec2(10.f, 1.25f);
			break;
		case GraphPreset::tfxFrameratePreset:
			return tfxVec2(0.0017f, .5f);
			break;
		case GraphPreset::tfxDimensionsPreset:
		case GraphPreset::tfxVelocityPreset:
		case GraphPreset::tfxWeightPreset:
		case GraphPreset::tfxWeightVariationPreset:
		case GraphPreset::tfxSpinPreset:
		case GraphPreset::tfxSpinVariationPreset:
			return tfxVec2(10.f, 2.5f);
			break;
		case GraphPreset::tfxNoiseResolutionPreset:
			return tfxVec2(10.f, 1.f);
			break;
		case GraphPreset::tfxNoiseOffsetVariationPreset:
			return tfxVec2(10.f, .01f);
			break;
		case GraphPreset::tfxDirectionOvertimePreset:
			return tfxVec2(0.0017f, 1.f);
			break;
		case GraphPreset::tfxWeightOvertimePreset:
		case GraphPreset::tfxVelocityOvertimePreset:
		case GraphPreset::tfxSpinOvertimePreset:
		case GraphPreset::tfxDirectionVariationPreset:
		case GraphPreset::tfxPercentOvertime:
			return tfxVec2(0.0017f, 0.0035f);
			break;
		case GraphPreset::tfxIntensityOvertimePreset:
			return tfxVec2(0.0017f, 0.01115f);
			break;
		case GraphPreset::tfxColorPreset:
			break;
		default:
			return tfxVec2(0.1f, 0.1f);
			break;
		}

		return tfxVec2(0.1f, 0.1f);
	}

	void CompileGraph(Graph &graph) {
		float last_frame = graph.GetLastFrame();
		graph.lookup.last_frame = unsigned int(last_frame / tfxLOOKUP_FREQUENCY);
		graph.lookup.values.clear();
		if (graph.lookup.last_frame) {
			graph.lookup.values.resize(graph.lookup.last_frame + 1);
			for (unsigned int f = 0; f != graph.lookup.last_frame + 1; ++f) {
				graph.lookup.values[f] = graph.GetValue((float)f * tfxLOOKUP_FREQUENCY);
			}
			graph.lookup.values[graph.lookup.last_frame] = graph.GetLastValue();
		}
		else {
			graph.lookup.values.push_back(graph.GetFirstValue());
		}
	}

	void CompileGraphOvertime(Graph &graph) {
		graph.lookup.values.clear();
		if (graph.nodes.size() > 1) {
			graph.lookup.last_frame = unsigned int(graph.lookup.life / tfxLOOKUP_FREQUENCY_OVERTIME);
			graph.lookup.values.resize(graph.lookup.last_frame + 1);
			for (unsigned int f = 0; f != graph.lookup.last_frame + 1; ++f) {
				graph.lookup.values[f] = graph.GetValue((float)f * tfxLOOKUP_FREQUENCY_OVERTIME, graph.lookup.life);
			}
			graph.lookup.values[graph.lookup.last_frame] = graph.GetLastValue();
		}
		else {
			graph.lookup.last_frame = 0;
			graph.lookup.values.push_back(graph.GetFirstValue());
		}
	}

	float LookupFastOvertime(Graph &graph, float age, float lifetime) {
		u32 frame = static_cast<u32>((age / lifetime * graph.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);
		frame = std::min<u32>(frame, graph.lookup.last_frame);
		return graph.lookup.values[frame];
	}

	float LookupFast(Graph &graph, float frame) {
		if ((unsigned int)frame < graph.lookup.last_frame)
			return graph.lookup.values[(unsigned int)frame];
		return graph.lookup.values[graph.lookup.last_frame];
	}

	float LookupPreciseOvertime(Graph &graph, float age, float life) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		AttributeNode *lastec = nullptr;
		for (auto &a : graph.nodes) {
			float frame = a.frame * life;
			if (age < frame) {
				p = (age - lastf) / (frame - lastf);
				float bezier_value = GetBezierValue(lastec, a, p, graph.min.y, graph.max.y);
				if (bezier_value) {
					return bezier_value;
				}
				else {
					return lastv - p * (lastv - a.value);
				}
			}
			lastv = a.value;
			lastf = frame - 1;
			lastec = &a;
		}
		return lastv;
	}

	float LookupPrecise(Graph &graph, float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		AttributeNode *lastec = nullptr;
		for (auto &a : graph.nodes) {
			if (age < a.frame) {
				p = (age - lastf) / (a.frame - lastf);
				float bezier_value = GetBezierValue(lastec, a, p, graph.min.y, graph.max.y);
				if (bezier_value) {
					return bezier_value;
				}
				else {
					return lastv - p * (lastv - a.value);
				}
			}
			lastv = a.value;
			lastf = a.frame - 1;
			lastec = &a;
		}
		return lastv;
	}

	float GetRandomFast(Graph &graph, float frame) {
		float value = 0;
		if ((unsigned int)frame < graph.lookup.last_frame)
			value = graph.lookup.values[(unsigned int)frame];
		value = graph.lookup.values[graph.lookup.last_frame];
		return random_generation.Range(value);
	}

	float GetRandomPrecise(Graph &graph, float frame) {
		return graph.GetRandomValue(frame);
	}

	float GetMaxLife(EffectEmitter &e) {
		Graph &life = *e.GetGraphByType(tfxBase_life);
		Graph &life_variation = *e.GetGraphByType(tfxVariation_life);
		float templife = 0;
		float max_life = 0;
		float life_last_frame = life.GetLastFrame();
		float life_variation_last_frame = life_variation.GetLastFrame();
		float global_adjust = 1.f;
		if (life_last_frame + life_variation_last_frame > 0) {
			for (float f = 0; f < std::fmaxf(life_last_frame, life_variation_last_frame); ++f) {
				if (e.parent)
					global_adjust = e.parent->GetGraphByType(tfxGlobal_life)->GetValue(f);
				templife = life.GetValue(f) + life_variation.GetValue(f);
				templife *= global_adjust;
				if (max_life < templife)
					max_life = templife;
			}
		}
		else {
			max_life = life.GetFirstValue() + life_variation.GetFirstValue();
		}

		return max_life;
	}

	float GetMaxAmount(EffectEmitter &e) {
		if (e.properties.flags & tfxEmitterPropertyFlags_single || e.properties.flags & tfxEmitterPropertyFlags_one_shot)
			return (float)e.properties.spawn_amount;
		Graph &amount = *e.GetGraphByType(tfxBase_amount);
		Graph &amount_variation = *e.GetGraphByType(tfxVariation_amount);
		float tempamount = 0;
		float max_amount = 0;
		float amount_last_frame = amount.GetLastFrame();
		float amount_variation_last_frame = amount_variation.GetLastFrame();
		float global_adjust = 1.f;
		if (amount_last_frame + amount_variation_last_frame > 0) {
			for (float f = 0; f < std::fmaxf(amount_last_frame, amount_variation_last_frame); ++f) {
				if (e.parent)
					global_adjust = e.parent->GetGraphByType(tfxGlobal_amount)->GetMaxValue();
				tempamount = amount.GetValue(f) + amount_variation.GetValue(f);
				tempamount *= global_adjust;
				if (max_amount < tempamount)
					max_amount = tempamount;
			}
		}
		else {
			max_amount = amount.GetFirstValue() + amount_variation.GetFirstValue();
		}

		return max_amount;
	}

	bool IsOvertimeGraph(GraphType type) {
		return type >= tfxOvertime_velocity && type != tfxOvertime_velocity_adjuster;
	}

	bool IsOvertimePercentageGraph(GraphType type) {
		return type >= tfxOvertime_velocity && type != tfxOvertime_velocity_adjuster && type != tfxOvertime_direction;
	}

	bool IsGlobalGraph(GraphType type) {
		return type >= tfxGlobal_life && type <= tfxGlobal_effect_angle;
	}

	bool IsGlobalPercentageGraph(GraphType type) {
		return type >= tfxGlobal_life && type <= tfxGlobal_splatter;
	}

	bool IsAngleGraph(GraphType type) {
		return (type == tfxGlobal_effect_angle || type == tfxProperty_emission_angle || type == tfxProperty_emission_range || type == tfxProperty_emitter_angle ||
			type == tfxProperty_arc_offset || type == tfxProperty_arc_size || type == tfxBase_spin || type == tfxVariation_spin);
	}

	bool IsAngleOvertimeGraph(GraphType type) {
		return type == tfxOvertime_direction;
	}

	bool IsEverythingElseGraph(GraphType type) {
		return !IsOvertimeGraph(type) && !IsOvertimePercentageGraph(type) && !IsGlobalGraph(type) && !IsAngleGraph(type) && !IsOvertimeGraph(type);
	}

	void tfxEffectPool::Init(unsigned int max_effects, unsigned int max_emitters, unsigned int max_particles) {
		effect_memory.reserve(sizeof(tfxEffect) * max_effects);
		emitter_memory.reserve(sizeof(tfxEmitter) * max_emitters);
		particle_memory.reserve(sizeof(Particle) * max_particles);
		sprite_memory.reserve(sizeof(ParticleSprite) * max_particles);
	}

	tfxEffectID tfxEffectPool::InsertEffect() {
		tfxEffectID id = effect_memory.get_range(sizeof(tfxEffect));
		void *ptr = (char*)effect_memory.data + effect_memory.ranges[id].offset_into_memory;
		new((void*)(ptr)) tfxEffect();
		return id;
	}

	tfxEffect &tfxEffectPool::GetEffect(tfxEffectID id) {
		assert(effect_memory.ranges.current_size > id);		//Not a valid Effect ID
		void *ptr = (char*)effect_memory.data + effect_memory.ranges[id].offset_into_memory;
		return *static_cast<tfxEffect*>(ptr);
	}

	bool HasDataValue(tfxStorageMap<DataEntry> &config, tfxText key) {
		return config.ValidName(key);
	}

	void AddDataValue(tfxStorageMap<DataEntry> &map, tfxText key, const char *value) {
		DataEntry entry;
		entry.type = tfxString;
		entry.key = key;
		entry.str_value = value;
		map.Insert(key, entry);
	}

	void AddDataValue(tfxStorageMap<DataEntry> &map, tfxText key, int value) {
		DataEntry entry;
		entry.type = tfxSInt;
		entry.key = key;
		entry.int_value = value;
		entry.bool_value = (bool)value;
		map.Insert(key, entry);
	}

	void AddDataValue(tfxStorageMap<DataEntry> &map, tfxText key, bool value) {
		DataEntry entry;
		entry.type = tfxBool;
		entry.key = key;
		entry.bool_value = value;
		entry.int_value = (int)value;
		map.Insert(key, entry);
	}

	void AddDataValue(tfxStorageMap<DataEntry> &map, tfxText key, double value) {
		DataEntry entry;
		entry.type = tfxDouble;
		entry.key = key;
		entry.double_value = value;
		map.Insert(key, entry);
	}

	void AddDataValue(tfxStorageMap<DataEntry> &map, tfxText key, float value) {
		DataEntry entry;
		entry.type = tfxFloat;
		entry.key = key;
		entry.float_value = value;
		map.Insert(key, entry);
	}

	tfxText &GetDataStrValue(tfxStorageMap<DataEntry> &map, const char* key) {
		return map.At(key).str_value;
	}
	int& GetDataIntValue(tfxStorageMap<DataEntry> &map, const char* key) {
		return map.At(key).int_value;
	}
	float& GetDataFloatValue(tfxStorageMap<DataEntry> &map, const char* key) {
		return map.At(key).float_value;
	}

	bool SaveDataFile(tfxStorageMap<DataEntry> &map, const char* path) {
		FILE *file = fopen(path, "w");

		if (file == NULL)
			return false;

		if (map.Size()) {
			for (auto &entry : map.data) {
				tfxText ini_line = entry.key;
				ini_line.Appendf("=");
				switch (entry.type) {
				case tfxString:
					ini_line.Appendf(entry.str_value.c_str());
					break;
				case tfxSInt:
					ini_line.Appendf("%i", entry.int_value);
					break;
				case tfxFloat:
					ini_line.Appendf("%f", entry.float_value);
					break;
				case tfxBool:
					ini_line.Appendf("%i", (int)entry.bool_value);
					break;
				}
				ini_line.Appendf("\n");
				fwrite(ini_line.c_str(), 1, ini_line.Length(), file);
			}
		}

		fclose(file);
		return true;

	}

	bool LoadDataFile(tfxStorageMap<DataEntry> &map, const char* path) {
		FILE* fp;
		fp = fopen(path, "r");
		if (fp == NULL) {
			return false;
		}

		const size_t max_line_length = 256;
		char buffer[max_line_length];

		while (fgets(buffer, max_line_length, fp)) {
			buffer[strcspn(buffer, "\n")] = 0;
			tfxText str = buffer;
			tfxvec<tfxText> pair = SplitString(str, 61);
			if (pair.size() == 2) {
				tfxText key = pair[0];
				DataType t = data_types.eff.At(pair[0]);
				if (t == tfxBool) {
					AddDataValue(map, key, (bool)atoi(pair[1].c_str()));
				}
				if (t == tfxSInt) {
					AddDataValue(map, key, atoi(pair[1].c_str()));
				}
				else if(t == tfxFloat) {
					AddDataValue(map, key, (float)atof(pair[1].c_str()));
				}
				else if (t == tfxString) {
					AddDataValue(map, key, pair[1].c_str());
				}
			}
		}

		int close = fclose(fp);
		return true;

	}

	tfxvec<tfxText> SplitString(const tfx::tfxText &str, char delim) {
		tfxvec<tfxText> ret;

		tfxText line;
		for (char c : str.string) {
			if (c == delim && line.Length() && c != NULL) {
				ret.push_back(line);
				line.Clear();
			}
			else if(c != NULL) {
				line.Append(c);
			}
		}

		if (line.Length()) {
			ret.push_back(line);
		}

		return ret;
	}

	bool StringIsUInt(const tfxText &s) {

		for (auto c : s.string) {
			if (!std::isdigit(c) && c != 0)
				return false;
		}

		return true;
	}

	int GetDataType(const tfxText &s) {
		if (s.Length() == 0)
			return tfxString;

		if(s.IsInt())
			return tfxSInt;

		if (s.IsFloat())
			return tfxFloat;

		return tfxString;
	}

	//Get a graph by GraphID
	Graph &GetGraph(EffectLibrary &library, GraphID &graph_id) {
		GraphType type = graph_id.type;

		if (type < tfxGlobalCount) {
			return ((Graph*)&library.global_graphs[graph_id.graph_id])[type];
		}
		else if (type >= tfxPropertyStart && type < tfxBaseStart) {
			int ref = type - tfxPropertyStart;
			return ((Graph*)&library.property_graphs[graph_id.graph_id])[ref];
		}
		else if (type >= tfxBaseStart && type < tfxVariationStart) {
			int ref = type - tfxBaseStart;
			return ((Graph*)&library.base_graphs[graph_id.graph_id])[ref];
		}
		else if (type >= tfxVariationStart && type < tfxOvertimeStart) {
			int ref = type - tfxVariationStart;
			return ((Graph*)&library.variation_graphs[graph_id.graph_id])[ref];
		}
		else if (type >= tfxOvertimeStart) {
			int ref = type - tfxOvertimeStart;
			return ((Graph*)&library.overtime_graphs[graph_id.graph_id])[ref];
		}

		assert(0);	//This function must return a value, make sure the graph_id is valid

		return((Graph*)&library.overtime_graphs[graph_id.graph_id])[type];

	}

	//API Functions
	//Get the number of shapes that are stored in an effects library saved on disk. This can be useful if you need to reserve the space in 
	//a list to store them in your custom ShapeLoader function.
	int GetShapesInPackage(const char *filename) {
		int context = 0;
		int error = 0;

		tfxPackage package;
		error = LoadPackage(filename, package);

		tfxEntryInfo *data = package.GetFile("data.txt");

		if (!data)
			error = -5;


		if (error < 0) {
			package.Free();
			return error;
		}

		int shape_count = 0;

		while (!data->data.EoF()) {
			tfxText line = data->data.ReadLine();
			bool context_set = false;
			if (StringIsUInt(line.c_str())) {
				context = atoi(line.c_str());
				if (context == tfxEndShapes)
					break;
				context_set = true;
			}
			if (context_set == false) {
				tfxvec<tfxText> pair = SplitString(line.c_str());
				if (pair.size() != 2) {
					pair = SplitString(line.c_str(), 44);
					if (pair.size() < 2) {
						error = 1;
						break;
					}
				}
				if (context == tfxStartShapes) {
					if (pair.size() >= 5) {
						int frame_count = atoi(pair[2].c_str());
						shape_count += frame_count;
					}
				}
			}
		}

		return shape_count;
	}

	int LoadEffectLibraryPackage(const char *filename, EffectLibrary &lib, void(*shape_loader)(const char* filename, ImageData &image_data, void *raw_image_data, int image_size, void *user_data), void *user_data) {
		assert(shape_loader);
		lib.Clear();

		tfxvec<EffectEmitter> effect_stack;
		int context = 0;
		int error = 0;
		int uid = 0;
		unsigned int current_global_graph = 0;

		tfxPackage package;
		error = LoadPackage(filename, package);

		tfxEntryInfo *data = package.GetFile("data.txt");

		if (!data)
			error = -5;

		if (error < 0) {
			package.Free();
			return error;
		}

		int first_shape_index = -1;

		while (!data->data.EoF()) {
			tfxText line = data->data.ReadLine();
			bool context_set = false;

			if (StringIsUInt(line.c_str())) {
				context = atoi(line.c_str());
				if (context == tfxEndOfFile)
					break;

				context_set = true;
				if (context == tfxStartFolder) {
					EffectEmitter effect;
					effect.library = &lib;
					effect.uid = uid++;
					effect.properties = EmitterProperties();
					effect.type = EffectEmitterType::tfxFolder;
					effect_stack.push_back(effect);
				}
				if (context == tfxStartEffect) {
					EffectEmitter effect;
					effect.library = &lib;
					if (effect_stack.size() <= 1) { //Only root effects get the global graphs
						lib.AddEffectGraphs(effect);
						effect.ResetEffectGraphs(false);
						current_global_graph = effect.global;
					}
					effect.uid = uid++;
					effect.properties = EmitterProperties();
					effect.type = EffectEmitterType::tfxEffectType;
					effect_stack.push_back(effect);
					lib.AddAnimationSettings(effect_stack.back());
				}
				if (context == tfxStartEmitter) {
					EffectEmitter emitter;
					emitter.library = &lib;
					lib.AddEmitterGraphs(emitter);
					emitter.uid = uid++;
					emitter.type = EffectEmitterType::tfxEmitterType;
					emitter.ResetEmitterGraphs(false);
					effect_stack.push_back(emitter);
				}
			}

			if (context_set == false) {
				tfxvec<tfxText> pair = SplitString(line.c_str());
				if (pair.size() != 2) {
					pair = SplitString(line.c_str(), 44);
					if (pair.size() < 2) {
						error = 1;
						break;
					}
				}

				if (context == tfxStartAnimationSettings || context == tfxStartEmitter || context == tfxStartEffect || context == tfxStartFolder) {
					switch (data_types.eff.At(pair[0])) {
					case tfxUint:
						AssignEffectorProperty(effect_stack.back(), pair[0], (unsigned int)atoi(pair[1].c_str()));
						break;
					case tfxFloat:
						AssignEffectorProperty(effect_stack.back(), pair[0], (float)atof(pair[1].c_str()));
						break;
					case tfxSInt:
						AssignEffectorProperty(effect_stack.back(), pair[0], atoi(pair[1].c_str()));
						break;
					case tfxBool:
						AssignEffectorProperty(effect_stack.back(), pair[0], (bool)(atoi(pair[1].c_str())));
						break;
					case tfxString:
						AssignEffectorProperty(effect_stack.back(), pair[0], pair[1]);
						break;
					}
				}

				if (context == tfxStartGraphs && effect_stack.back().type == tfxEmitterType) {
					AssignGraphData(effect_stack.back(), pair);
				}
				else if (context == tfxStartGraphs && effect_stack.back().type == tfxEffectType) {
					if (effect_stack.size() <= 2)
						AssignGraphData(effect_stack.back(), pair);
				}

				if (context == tfxStartShapes) {
					if (pair.size() >= 5) {
						ShapeData s;
						strcpy_s(s.name, pair[0].c_str());
						s.shape_index = atoi(pair[1].c_str());
						s.frame_count = atoi(pair[2].c_str());
						s.width = atoi(pair[3].c_str());
						s.height = atoi(pair[4].c_str());
						if (pair.size() > 5)
							s.import_filter = atoi(pair[5].c_str());
						if (s.import_filter < 0 || s.import_filter>1)
							s.import_filter = 0;

						tfxEntryInfo *shape_entry = package.GetFile(s.name);
						if (shape_entry) {
							ImageData image_data;
							image_data.animation_frames = (float)s.frame_count;
							image_data.image_size = tfxVec2((float)s.width, (float)s.height);
							image_data.shape_index = s.shape_index;
							image_data.import_filter = s.import_filter;

							shape_loader(s.name, image_data, shape_entry->data.data, (u32)shape_entry->file_size, user_data);

							if (!image_data.ptr) {
								//uid = -6;
							}
							else {
								lib.particle_shapes.InsertByInt(s.shape_index, image_data);
								if (first_shape_index == -1)
									first_shape_index = s.shape_index;
							}
						}
						else {
							//Maybe don't actually need to break here, just means for some a reason a shaped couldn't be loaded, but no reason not to load the effects anyway
							//uid = -7;
							//break;
						}
					}
				}
			}

			if (context == tfxEndEmitter) {
				effect_stack.back().InitialiseUninitialisedGraphs();
				effect_stack.back().UpdateMaxLife();
				effect_stack.parent().sub_effectors.push_back(effect_stack.back());
				effect_stack.pop();
			}

			if (context == tfxEndEffect) {
				effect_stack.back().ReIndex();
				if (effect_stack.size() > 1)
					effect_stack.parent().sub_effectors.push_back(effect_stack.back());
				else {
					lib.effects.push_back(effect_stack.back());
					effect_stack.back().InitialiseUninitialisedGraphs();
				}
				effect_stack.pop();
			}

			if (context == tfxEndFolder) {
				assert(effect_stack.size() == 1);			//Folders should not be contained within anything
				lib.effects.push_back(effect_stack.back());
				effect_stack.pop();
			}

		}

		package.Free();
		//Returning anything over 0 means that effects were loaded ok
		//-1 to -4 = Package not in correct format
		//-5 = Data in package could not be loaded
		//-6 = ShapeLoader did not add a pointer into the image data
		//-7 = A shape image could not be loaded from the package

		if (uid >= 0) {
			lib.CompileAllGraphs();
			lib.ReIndex();
			if(first_shape_index != -1)
				lib.UpdateParticleShapeReferences(lib.effects, first_shape_index);
			lib.UpdateEffectPaths();
			lib.UpdateComputeNodes();
			lib.UpdateEffectParticleStorage();
			lib.SetMinMaxData();
		}

		return uid - 1;

	}

	void tfxEffectTemplate::SetUserDataAll(void *data) {
		tfxvec<EffectEmitter*> stack;
		stack.push_back(&effect_template);
		while (stack.size()) {
			EffectEmitter *current = stack.pop_back();
			current->user_data = data;
			for (auto &sub : current->sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	void tfxEffectTemplate::SetParticleOnSpawnCallback(tfxText path, void(*particle_onspawn_callback)(Particle &particle)) {
		assert(paths.ValidName(path));
		EffectEmitter &e = *paths.At(path);
		assert(e.type == tfxEmitterType);
		e.particle_onspawn_callback = particle_onspawn_callback;
	}

}