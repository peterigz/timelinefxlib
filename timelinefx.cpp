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
		particles.free_range(pm->particle_memory);
	}

	void EffectEmitter::SoftExpire() {
		flags |= tfxEmitterStateFlags_stop_spawning;
	}

	void EffectEmitter::Rotate(float r) {
		transform.local.rotation += r;
	}

	void EffectEmitter::SetAngle(float a) {
		transform.local.rotation = a;
	}

	void EffectEmitter::Scale(const tfxVec2& s) {
		transform.world.scale = s;
	}

	void EffectEmitter::Scale(float x, float y) {
		transform.world.scale.x = x;
		transform.world.scale.y = y;
	}

	void EffectEmitter::Move(const tfxVec2& m) {
		transform.local.position += m;
	}

	void EffectEmitter::Move(float x, float y) {
		transform.local.position.x += x;
		transform.local.position.y += y;
	}

	void EffectEmitter::Position(const tfxVec2& p) {
		transform.local.position = p;
	}

	void EffectEmitter::Position(float x, float y) {
		transform.local.position.x = x;
		transform.local.position.y = y;
	}

	void EffectEmitter::UpdateMaxLife() {
		max_life = GetMaxLife(*this);
		GetGraphByType(tfxOvertime_red)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_green)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_blue)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_opacity)->lookup.life = max_life;
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

	void EffectEmitter::UpdateAllBufferSizes(unsigned int totals[tfxLAYERS]) {
		if (type == tfxEmitterType) {
			UpdateBufferSize(totals);
		}
		else {
			for (auto &emitter : sub_effectors) {
				emitter.UpdateAllBufferSizes(totals);
			}
		}
		if (!parent) {
			for (int i = 0; i != tfxLAYERS; ++i) {
				max_particles[i] = totals[i];
			}
		}
	}

	void EffectEmitter::UpdateBufferSize(unsigned int totals[tfxLAYERS]) {
		if (properties.flags & tfxEmitterPropertyFlags_single || properties.flags & tfxEmitterPropertyFlags_one_shot)
			max_particles[properties.layer] += properties.spawn_amount;
		else
			max_particles[properties.layer] += unsigned int(GetMaxAmount(*this) * (GetMaxLife(*this) / 1000.f)) + 1;
		totals[properties.layer] += max_particles[properties.layer];
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

	void UpdateEffect(tfxEffect &e) {
		if (e.sub_emitters.empty()) return;
		e.transform.captured = e.transform.world;

		if (e.common.lookup_mode == tfxPrecise) {
			e.common.frame = e.common.age;
		}
		else {
			e.common.frame = e.common.age / tfxLOOKUP_FREQUENCY;
		}

		//Update the effect state
		e.current.life = lookup_callback(e.common.library->global_graphs[e.global].life, e.common.frame);
		e.current.amount = lookup_callback(e.common.library->global_graphs[e.global].amount, e.common.frame);
		if (!(e.common.property_flags & tfxEmitterPropertyFlags_global_uniform_size)) {
			e.current.size.x = lookup_callback(e.common.library->global_graphs[e.global].width, e.common.frame);
			e.current.size.y = lookup_callback(e.common.library->global_graphs[e.global].height, e.common.frame);
		}
		else {
			e.current.size.x = lookup_callback(e.common.library->global_graphs[e.global].width, e.common.frame);
			e.current.size.y = e.current.size.x;
		}
		e.current.velocity = lookup_callback(e.common.library->global_graphs[e.global].velocity, e.common.frame);
		e.current.spin = lookup_callback(e.common.library->global_graphs[e.global].spin, e.common.frame);
		e.current.opacity = lookup_callback(e.common.library->global_graphs[e.global].opacity, e.common.frame);
		e.current.splatter = lookup_callback(e.common.library->global_graphs[e.global].splatter, e.common.frame);
		//We don't want to scale twice when the sub effect is transformed, so the values here are set to 1. That means that the root effect will only control the global scale.
		e.current.overal_scale = lookup_callback(e.common.library->global_graphs[e.global].overal_scale, e.common.frame);
		if (!e.parent_particle) {
			e.transform.world.scale.x = e.current.overal_scale;
			e.transform.world.scale.y = e.current.overal_scale;
			e.transform.local.rotation = lookup_callback(e.common.library->global_graphs[e.global].effect_angle, e.common.frame);
		}
		else {
			e.transform.world.scale.x = e.current.overal_scale;
			e.transform.world.scale.y = e.current.overal_scale;
			e.transform.local.rotation = 0.f;
		}
		e.current.stretch = lookup_callback(e.common.library->global_graphs[e.global].stretch, e.common.frame);
		e.current.weight = lookup_callback(e.common.library->global_graphs[e.global].weight, e.common.frame);

		if (e.parent_particle) {
			e.common.state_flags |= e.parent_particle->flags & tfxParticleFlags_remove;
			if (e.parent_particle->next_ptr) {
				e.parent_particle = e.parent_particle->next_ptr;
				Transform(e, *e.parent_particle);
				e.transform.world.position += e.common.handle * e.current.overal_scale;

				if (e.common.state_flags & tfxEmitterStateFlags_no_tween_this_update) {
					e.transform.captured = e.transform.world;
				}

			}
			else {
				e.parent_particle = nullptr;
				e.common.state_flags |= tfxEmitterStateFlags_retain_matrix;
				e.transform.local.position = e.transform.world.position;
				e.transform.local.rotation = e.transform.world.rotation;
				e.common.state_flags |= tfxEmitterStateFlags_stop_spawning;
			}
		}
		else {
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
		}

		e.common.age += FRAME_LENGTH;
		if (!(e.common.property_flags & tfxEmitterPropertyFlags_single) || e.common.property_flags & tfxEmitterPropertyFlags_one_shot)
			e.common.highest_particle_age -= FRAME_LENGTH;

		if (e.common.loop_length && e.common.age > e.common.loop_length)
			e.common.age = 0;

		if (e.common.highest_particle_age <= 0) {
			e.common.timeout_counter++;
		}
		else {
			e.common.timeout_counter = 0;
		}

		e.common.state_flags &= ~tfxEmitterStateFlags_no_tween_this_update;

		//Now udpate the Children
		e.sub_emitters.reset();
		while (!e.sub_emitters.eob()) {
			tfxEmitter &emitter = e.sub_emitters.next();
			emitter.UpdateEmitter();
		}

		for (int layer = 0; layer != tfxLAYERS; ++layer) {
			if(e.particles[layer][e.current_buffer].current_size)
				e.ControlParticles(layer);
			e.particles[layer][e.current_buffer].clear();
			e.sprites[layer][e.current_buffer].clear();
		}
		e.current_buffer = !e.current_buffer;

	}

	void tfxEmitter::UpdateEmitter () {
		transform.captured = transform.world;

		if (common.lookup_mode == tfxPrecise) {
			common.frame = common.age;
		}
		else {
			common.frame = common.age / tfxLOOKUP_FREQUENCY;
		}

		common.state_flags |= (common.root_effect->common.state_flags & tfxEmitterStateFlags_remove);
		transform.local.rotation = lookup_callback(common.library->property_graphs[property].emitter_angle, common.frame);
		current.velocity_adjuster = lookup_callback(common.library->overtime_graphs[overtime].velocity_adjuster, common.frame);
		current.overal_scale = common.root_effect->current.overal_scale;
		current.stretch = common.root_effect->current.stretch;

		current.emitter_size.y = lookup_callback(common.library->property_graphs[property].emitter_height, common.frame);
		if (emission_type == EmissionType::tfxArea || emission_type == EmissionType::tfxEllipse) {
			current.emitter_size.x = lookup_callback(common.library->property_graphs[property].emitter_width, common.frame);
		}
		else
			current.emitter_size.x = 0.f;
		if (common.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center && emission_type != EmissionType::tfxPoint) {
			common.handle = current.emitter_size * -0.5f;
		}
		else if (!(common.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center)) {
			common.handle = emitter_handle;
		}
		if (common.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center && emission_type == EmissionType::tfxLine) {
			common.handle = current.emitter_size * 0.5f;
		}
		Transform(*this, *parent_effect);

		if (common.state_flags & tfxEmitterStateFlags_no_tween_this_update) {
			transform.captured = transform.world;
		}

		SpawnParticles();

		parent_effect->common.highest_particle_age = common.highest_particle_age;

		common.age += FRAME_LENGTH;
		if (!(common.property_flags & tfxEmitterPropertyFlags_single) || common.property_flags & tfxEmitterPropertyFlags_one_shot)
			common.highest_particle_age -= FRAME_LENGTH;

		if (common.loop_length && common.age > common.loop_length)
			common.age = 0;

		if (common.highest_particle_age <= 0) {
			common.timeout_counter++;
		}
		else {
			common.timeout_counter = 0;
		}

		common.state_flags &= ~tfxEmitterStateFlags_no_tween_this_update;

	}

	void tfxEmitter::SpawnParticles() {
		if (common.state_flags & tfxEmitterStateFlags_single_shot_done || parent_effect->common.state_flags & tfxEmitterStateFlags_stop_spawning)
			return;

		float qty;
		if (!(common.property_flags & tfxEmitterPropertyFlags_single) && !(common.property_flags & tfxEmitterPropertyFlags_one_shot)) {
			qty = lookup_callback(common.library->base_graphs[base].amount, common.frame);
			qty += random_generation.Range(lookup_callback(common.library->variation_graphs[variation].amount, common.frame));

			if (common.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && (common.state_flags & tfxEmitterStateFlags_is_area)) {
				float area = current.emitter_size.x * current.emitter_size.y;
				qty = (qty / 10000.f) * area;
			}
			else if (common.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && common.state_flags & tfxEmitterStateFlags_is_line) {
				qty = (qty / 100.f) * current.emitter_size.y;
			}

			qty *= common.root_effect->current.amount;
			qty *= UPDATE_TIME;
			qty += current.amount_remainder;
		}
		else {
			qty = (float)spawn_amount;
		}

		float tween = 0.f;
		float interpolate = (float)(int)qty;
		float count = 0;

		unsigned int current_buffer = common.root_effect->current_buffer;

		while (qty >= 1.f) {
			if (!common.root_effect->FreeCapacity(layer)) {
				current.amount_remainder = 0;
				break;
			}
			tween = count / interpolate;
			count++;
			qty -= 1.f;

			bool is_single = common.property_flags & tfxEmitterPropertyFlags_single && !(common.property_flags & tfxEmitterPropertyFlags_one_shot);

			Particle &p = common.root_effect->GrabParticle(layer);
			InitCPUParticle(p, tween);

			ParticleSprite &s = common.root_effect->sprites[layer][current_buffer].back();
			s.parameters = (blend_mode << 28) + (unsigned int)p.image_frame;
			s.color = p.color;
			s.world = p.world;
			s.captured = p.captured;
			s.ptr = image_ptr;
			s.intensity = p.intensity;
			s.handle = image_handle;
			p.sprite_index = common.root_effect->sprites[layer][current_buffer].last_index();

			common.highest_particle_age = std::fmaxf(common.highest_particle_age, p.max_age);
		}

		current.amount_remainder = qty;
	}

	void EffectEmitter::Update() {
		transform.captured = transform.world;

		if (pm->lookup_mode == tfxPrecise) {
			current.frame = current.age;
		}
		else {
			current.frame = current.age / tfxLOOKUP_FREQUENCY;
		}

		if (type == EffectEmitterType::tfxEffectType) {
			UpdateEffectState();
		}

		if (parent && parent->type != tfxFolder) {
			parent = parent->next_ptr;
			flags |= (parent->flags & tfxEmitterStateFlags_remove);
			UpdateEmitterState();
			TransformEffector(*parent, true, properties.flags & tfxEmitterPropertyFlags_relative_angle);
			if (pm->use_compute_shader && properties.flags & tfxEmitterPropertyFlags_is_bottom_emitter)
				UpdateComputeController();

			if (flags & tfxEmitterStateFlags_no_tween_this_update) {
				transform.captured = transform.world;
			}

			bool is_compute = properties.flags & tfxEmitterPropertyFlags_is_bottom_emitter && pm->use_compute_shader;

			if (!is_compute && pm->contain_particles_in_emitters)
				ControlParticles(*this);

			GetRootEffect()->BumpSprites();

			if (!pm->disable_spawing && pm->FreeCapacity(properties.layer, is_compute))
				SpawnParticles();

			parent->highest_particle_age = highest_particle_age;
		}
		else if (parent_particle) {
			flags |= parent_particle->flags & tfxParticleFlags_remove;
			if (parent_particle->next_ptr) {
				parent_particle = parent_particle->next_ptr;
				TransformEffector(*parent_particle, true, properties.flags & tfxEmitterPropertyFlags_relative_angle);
				transform.world.position += properties.emitter_handle * current.overal_scale;

				if (flags & tfxEmitterStateFlags_no_tween_this_update) {
					transform.captured = transform.world;
				}

			}
			else {
				parent_particle = nullptr;
				flags |= tfxEmitterStateFlags_retain_matrix;
				transform.local.position = transform.world.position;
				transform.local.rotation = transform.world.rotation;
				flags |= tfxEmitterStateFlags_stop_spawning;
			}
		}
		else {
			if (!(flags & tfxEmitterStateFlags_retain_matrix)) {
				transform.world.position = transform.local.position;
				transform.world.rotation = transform.local.rotation;
				transform.world.position += properties.emitter_handle * current.overal_scale;
				float s = sin(transform.local.rotation);
				float c = cos(transform.local.rotation);
				transform.matrix.Set(c, s, -s, c);
			}

			if (flags & tfxEmitterStateFlags_no_tween_this_update) {
				transform.captured = transform.world;
			}
		}

		current.age += FRAME_LENGTH;
		if(!(properties.flags & tfxEmitterPropertyFlags_single) || properties.flags & tfxEmitterPropertyFlags_one_shot)
			highest_particle_age -= FRAME_LENGTH;

		if (properties.loop_length && current.age > properties.loop_length)
			current.age = 0;

		if (highest_particle_age <= 0) {
			timeout_counter++;
			if (parent && parent->type != tfxFolder && timeout_counter >= timeout)
				parent->active_children--;
		}
		else {
			timeout_counter = 0;
		}

		flags &= ~tfxEmitterStateFlags_no_tween_this_update;
	}

	void EffectEmitter::UpdateParticles() {
		bool can_bump = true;
		particles.reset();
		while(!particles.eob()) {
			Particle &p = particles.next();
			if (!(p.flags & tfxParticleFlags_fresh)) {
				p.captured = p.world;

				if (ControlParticleFast(p, *this)) {
					TransformParticle(p, *this);
					if (p.flags & tfxParticleFlags_capture_after_transform) {
						p.captured.position = p.world.position;
						p.flags &= ~tfxParticleFlags_capture_after_transform;
					}
					can_bump = false;
				}
				else if (can_bump) {
					particles.bump();
				}
			}
			else {
				p.flags &= ~tfxParticleFlags_fresh;
			}
		}
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

	void EffectEmitter::SpawnParticles() {
		if (current.single_shot_done || parent->flags & tfxEmitterStateFlags_stop_spawning)
			return;

		float qty;
		if (!(properties.flags & tfxEmitterPropertyFlags_single) && !(properties.flags & tfxEmitterPropertyFlags_one_shot)) {
			qty = current.amount;
			qty += random_generation.Range(current.amount_variation);

			if (properties.flags & tfxEmitterPropertyFlags_use_spawn_ratio && (properties.emission_type == tfxArea || properties.emission_type == tfxEllipse)) {
				float area = current.emitter_size.x * current.emitter_size.y;
				qty = (qty / 10000.f) * area;
			}
			else if (properties.flags & tfxEmitterPropertyFlags_use_spawn_ratio && properties.emission_type == tfxLine) {
				qty = (qty / 100.f) * current.emitter_size.y;
			}

			qty *= lookup_callback(parent->library->global_graphs[parent->global].amount, current.frame);
			qty *= UPDATE_TIME;
			qty += current.amount_remainder;
		}
		else {
			qty = (float)properties.spawn_amount;
		}

		float tween = 0.f;
		float interpolate = (float)(int)qty;
		float count = 0;
		bool is_compute = properties.flags & tfxEmitterPropertyFlags_is_bottom_emitter && pm->use_compute_shader;
		EffectEmitter &root_effect = *GetRootEffect();

		while (qty >= 1.f) {
			if (!pm->contain_particles_in_emitters && !pm->FreeCapacity(properties.layer, is_compute)) {
				current.amount_remainder = 0;
				break;
			}
			else if (pm->contain_particles_in_emitters && !FreeCapacity()) {
				current.amount_remainder = 0;
				break;
			}
			tween = count / interpolate;
			count++;
			qty -= 1.f;

			bool is_single = properties.flags & tfxEmitterPropertyFlags_single && !(properties.flags & tfxEmitterPropertyFlags_one_shot);

			if (properties.flags & tfxEmitterPropertyFlags_is_bottom_emitter && pm->use_compute_shader && !is_single){
				ComputeParticle &p = pm->GrabComputeParticle(properties.layer);
				InitComputeParticle(p, tween);
				pm->new_particles_count++;
				highest_particle_age = std::fmaxf(highest_particle_age, p.max_age);
			}
			else {
				Particle &p = pm->contain_particles_in_emitters ? GrabParticle() : pm->GrabCPUParticle(properties.layer);
				InitCPUParticle(p, tween);

				if (pm->contain_particles_in_emitters) {
					if (root_effect.GrabSprite(properties.layer)) {
						ParticleSprite &s = root_effect.sprites[properties.layer].back();
						s.parameters = (properties.blend_mode << 28) + (unsigned int)p.image_frame;
						s.color = p.color;
						s.world = p.world;
						s.captured = p.captured;
						s.ptr = properties.image->ptr;
						s.intensity = p.intensity;
						s.handle = current.image_handle;
						p.sprite_index = root_effect.sprites[properties.layer].last_index();
					}
				}

				highest_particle_age = std::fmaxf(highest_particle_age, p.max_age);
			}
		}

		current.amount_remainder = qty;
	}

	float tfxEmitter::GetEmissionDirection(tfxVec2 &local_position, tfxVec2 &world_position, tfxVec2 &emitter_size) {
		float emission_angle = lookup_callback(common.library->property_graphs[property].emission_angle, common.frame);
		float emission_angle_variation = lookup_callback(common.library->property_graphs[property].emission_range, common.frame);
		//----Emission
		float range = emission_angle_variation *.5f;
		float direction = 0;

		if (emission_type == EmissionType::tfxPoint)
			return direction + emission_angle + random_generation.Range(-range, range);

		tfxVec2 tmp_position;
		if (local_position.x == 0 && local_position.y == 0)
			tmp_position = emitter_size;
		else
			tmp_position = local_position;

		if (emission_direction == EmissionDirection::tfxOutwards) {

			tfxVec2 to_handle;

			if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = (tmp_position);
			else
				to_handle = (world_position - transform.world.position);

			direction = GetVectorAngle(to_handle.x, to_handle.y);

		}
		else if (emission_direction == EmissionDirection::tfxInwards) {

			tfxVec2 to_handle;

			if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = (-tmp_position);
			else
				to_handle = (transform.world.position - world_position);

			direction = GetVectorAngle(to_handle.x, to_handle.y);

		}
		else if (emission_direction == EmissionDirection::tfxBothways) {

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

	Particle& tfxEffect::GrabParticle(unsigned int layer) {
		assert(particles[layer][current_buffer].current_size != particles[layer][current_buffer].capacity);
		unsigned int size = particles[layer][current_buffer].current_size++;
		sprites[layer][current_buffer].current_size++;
		return particles[layer][current_buffer][size];
	}

	bool tfxEffect::FreeCapacity(unsigned int layer) {
		return !particles[layer][current_buffer].full();
	}

	void tfxEmitter::InitCPUParticle(Particle &p, float tween) {
		p.flags = tfxParticleFlags_fresh;
		p.next_ptr = &p;
		p.emitter = this;

		if (common.property_flags & (tfxEmitterPropertyFlags_single | tfxEmitterPropertyFlags_one_shot))
			common.state_flags |= tfxEmitterStateFlags_single_shot_done;

		float life = lookup_callback(common.library->base_graphs[base].life, common.frame) * common.root_effect->current.life;
		float life_variation = lookup_callback(common.library->variation_graphs[variation].life, common.frame) * common.root_effect->current.life;

		float arc_size = 0.f;
		float arc_offset = 0.f;
		if (emission_type == EmissionType::tfxEllipse) {
			arc_size = lookup_callback(common.library->property_graphs[property].arc_size, common.frame);
			arc_offset = lookup_callback(common.library->property_graphs[property].arc_offset, common.frame);
		}
		float weight = lookup_callback(common.library->base_graphs[base].weight, common.frame) * common.root_effect->current.weight;
		float weight_variation = lookup_callback(common.library->variation_graphs[variation].weight, common.frame) * common.root_effect->current.weight;
		float velocity = lookup_callback(common.library->base_graphs[base].velocity, common.frame) * common.root_effect->current.velocity;
		float velocity_variation = lookup_callback(common.library->variation_graphs[variation].velocity, common.frame) * common.root_effect->current.velocity;
		tfxVec2 size;
		if (!(common.property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			size.x = lookup_callback(common.library->base_graphs[base].width, common.frame) * common.root_effect->current.size.x;
			size.y = lookup_callback(common.library->base_graphs[base].height, common.frame) * common.root_effect->current.size.y;
		}
		else {
			size.x = lookup_callback(common.library->base_graphs[base].width, common.frame);
			if (common.root_effect->common.property_flags & tfxEmitterPropertyFlags_global_uniform_size)
				size.y = size.x * common.root_effect->current.size.x;
			else
				size.y = size.x * common.root_effect->current.size.y;
			size.x *= common.root_effect->current.size.x;
		}
		tfxVec2 size_variation;
		size_variation.x = lookup_callback(common.library->variation_graphs[variation].width, common.frame) * common.root_effect->current.size.x;
		size_variation.y = lookup_callback(common.library->variation_graphs[variation].height, common.frame) * common.root_effect->current.size.y;
		float spin = lookup_callback(common.library->base_graphs[base].spin, common.frame) * common.root_effect->current.spin;
		float spin_variation = lookup_callback(common.library->variation_graphs[variation].spin, common.frame) * common.root_effect->current.spin;
		float splatter = lookup_callback(common.library->property_graphs[property].splatter, common.frame) * common.root_effect->current.splatter;
		float noise_offset_variation = lookup_callback(common.library->variation_graphs[variation].noise_offset, common.frame);
		float noise_offset = lookup_callback(common.library->base_graphs[variation].noise_offset, common.frame);
		float noise_resolution = lookup_callback(common.library->variation_graphs[variation].noise_resolution, common.frame);

		tfxVec2 grid_segment_size;
		if (common.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {
			if (emission_type == EmissionType::tfxArea) {
				if (grid_points.x > 1)
					grid_segment_size.x = current.emitter_size.x / (grid_points.x - 1);
				if (grid_points.y > 1)
					grid_segment_size.y = current.emitter_size.y / (grid_points.y - 1);
			}
			else if (emission_type == EmissionType::tfxEllipse) {
				if (grid_points.x > 0)
					grid_segment_size.x = arc_size / (grid_points.x);
			}
			else if (emission_type == EmissionType::tfxLine) {
				if (grid_points.x > 1)
					grid_segment_size.y = current.emitter_size.y / (grid_points.x - 1);
			}
		}

		//----Life
		p.max_age = life + random_generation.Range(life_variation);
		p.age = 0.f;

		//----Position
		if (emission_type == EmissionType::tfxPoint) {
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
		else if (emission_type == EmissionType::tfxArea) {
			tfxVec2 position = tfxVec2(0.f, 0.f);

			if (common.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				if (common.property_flags & tfxEmitterPropertyFlags_fill_area) {
					if (!(common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
						current.grid_coords.x--;
						if (current.grid_coords.x < 0.f) {
							current.grid_coords.y--;
							current.grid_coords.x = grid_points.x - 1;
							if (current.grid_coords.y < 0.f)
								current.grid_coords.y = grid_points.y - 1;
						}
					}

					p.local.position = position + (current.grid_coords * grid_segment_size) + common.handle;

					if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
						current.grid_coords.x++;
						if (current.grid_coords.x == grid_points.x) {
							current.grid_coords.y++;
							current.grid_coords.x = 0.f;
							if (current.grid_coords.y >= grid_points.y)
								current.grid_coords.y = 0.f;
						}
					}
				}
				else {

					if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {

						current.grid_direction.x = 1;
						current.grid_direction.y = 0;
						if (current.grid_coords.x == grid_points.x - 1 && current.grid_coords.y >= 0 && current.grid_coords.y < grid_points.y - 1) {
							current.grid_direction.x = 0;
							current.grid_direction.y = 1;
						}
						else if (current.grid_coords.x > 0 && current.grid_coords.x < grid_points.x && current.grid_coords.y == grid_points.y - 1) {
							current.grid_direction.x = -1;
							current.grid_direction.y = 0;
						}
						else if (current.grid_coords.x == 0 && current.grid_coords.y > 0 && current.grid_coords.y < grid_points.y) {
							current.grid_direction.x = 0;
							current.grid_direction.y = -1;
						}

					}
					else {

						current.grid_direction.x = -1;
						current.grid_direction.y = 0;
						if (current.grid_coords.x == grid_points.x - 1 && current.grid_coords.y > 0 && current.grid_coords.y < grid_points.y) {
							current.grid_direction.x = 0;
							current.grid_direction.y = -1;
						}
						else if (current.grid_coords.x >= 0 && current.grid_coords.x < grid_points.x - 1 && current.grid_coords.y == grid_points.y - 1) {
							current.grid_direction.x = 1;
							current.grid_direction.y = 0;
						}
						else if (current.grid_coords.x == 0 && current.grid_coords.y >= 0 && current.grid_coords.y < grid_points.y - 1) {
							current.grid_direction.x = 0;
							current.grid_direction.y = 1;
						}

					}

					current.grid_coords += current.grid_direction;
					tfxBound(current.grid_coords, grid_points);
					p.local.position = position + (current.grid_coords * grid_segment_size) + common.handle;
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
		else if (emission_type == EmissionType::tfxEllipse) {
			tfxVec2 emitter_size = (current.emitter_size * .5f);
			tfxVec2 position = tfxVec2(0.f, 0.f);

			if (common.property_flags & tfxEmitterPropertyFlags_spawn_on_grid && !(common.property_flags & tfxEmitterPropertyFlags_fill_area)) {

				current.grid_coords.y = 0.f;

				if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					current.grid_coords.x--;
					if (current.grid_coords.x < 0.f) {
						current.grid_coords.x = grid_points.x - 1;
					}
				}

				float th = current.grid_coords.x * grid_segment_size.x + arc_offset;
				p.local.position = tfxVec2(std::cosf(th) * emitter_size.x + common.handle.x + emitter_size.x,
					-std::sinf(th) * emitter_size.y + common.handle.y + emitter_size.y);

				if (!(common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					current.grid_coords.x++;
					if (current.grid_coords.x >= grid_points.x) {
						current.grid_coords.x = 0.f;
					}
				}

			}
			else if (!(common.property_flags & tfxEmitterPropertyFlags_fill_area)) {
				float th = random_generation.Range(arc_size) + arc_offset;

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
		else if (emission_type == EmissionType::tfxLine) {
			if (common.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				current.grid_coords.x = 0.f;

				if (!(common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					current.grid_coords.y--;
					if (current.grid_coords.y < 0.f) {
						current.grid_coords.y = grid_points.x - 1;
					}
				}

				p.local.position = tfxVec2(current.grid_coords * -grid_segment_size);
				p.local.position += common.handle;

				if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					current.grid_coords.y++;
					if (current.grid_coords.y >= grid_points.x) {
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
		if (weight) {
			p.base.weight = weight * common.library->overtime_graphs[overtime].weight.GetFirstValue();
			if (weight_variation > 0) {
				p.base.weight += random_generation.Range(-weight_variation, weight_variation) * common.library->overtime_graphs[overtime].weight.GetFirstValue();
			}
		}
		else {
			p.base.weight = 0;
		}
		p.weight_acceleration = p.base.weight * common.library->overtime_graphs[overtime].weight.GetFirstValue() * UPDATE_TIME;

		//----Velocity
		float velocity_scale = common.library->overtime_graphs[overtime].velocity.GetFirstValue() * current.velocity_adjuster;
		p.base.velocity = velocity + random_generation.Range(-velocity_variation, velocity_variation);

		//----Size
		if (!(common.property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			p.base.random_size.x = random_generation.Range(size_variation.x);
			p.base.random_size.y = random_generation.Range(size_variation.y);
			p.base.size.y = p.base.random_size.y + size.y;
			p.base.size.x = (p.base.random_size.x + size.x) / image_size.x;
			float height = p.base.size.y / image_size.y;

			p.world.scale.x = p.base.size.x * common.library->overtime_graphs[overtime].width.GetFirstValue();

			if (common.library->overtime_graphs[overtime].stretch.GetFirstValue()) {
				float velocity = std::fabsf(velocity_scale * p.base.velocity) * UPDATE_TIME;
				velocity += p.weight_acceleration * UPDATE_TIME;
				p.world.scale.y = (common.library->overtime_graphs[overtime].height.GetFirstValue() * common.root_effect->current.size.y * (p.base.size.y + (velocity * common.library->overtime_graphs[overtime].stretch.GetFirstValue() * common.root_effect->current.stretch))) / image_size.y;
			}
			else {
				if (common.property_flags & tfxEmitterPropertyFlags_lifetime_uniform_size) {
					p.world.scale.y = p.world.scale.x;
				}
				else {
					p.world.scale.y = height * common.library->overtime_graphs[overtime].height.GetFirstValue();
				}
			}
		}
		else {
			p.base.random_size.x = random_generation.Range(size_variation.x);
			p.base.random_size.y = p.base.random_size.x;
			p.base.size.y = p.base.random_size.y + size.y;
			p.base.size.x = (p.base.random_size.x + size.x) / image_size.x;
			float height = p.base.size.y / image_size.y;

			p.world.scale.x = p.base.size.x * common.library->overtime_graphs[overtime].width.GetFirstValue();

			if (common.library->overtime_graphs[overtime].stretch.GetFirstValue()) {
				float velocity = std::fabsf(velocity_scale * p.base.velocity) * UPDATE_TIME;
				velocity += p.weight_acceleration * UPDATE_TIME;
				p.world.scale.y = (common.library->overtime_graphs[overtime].width.GetFirstValue() * common.root_effect->current.size.y * (p.base.size.y + (velocity * common.library->overtime_graphs[overtime].stretch.GetFirstValue() * common.root_effect->current.stretch))) / image_size.y;
			}
			else {
				p.world.scale.y = p.world.scale.x;
			}
		}

		//----Spin
		p.base.spin = random_generation.Range(-spin_variation, std::abs(spin_variation)) + spin;

		switch (angle_setting) {
		case AngleSetting::tfxRandom:
			p.captured.rotation = p.local.rotation = random_generation.Range(angle_offset);
			break;
		case AngleSetting::tfxSpecify:
			p.captured.rotation = p.local.rotation = angle_offset;
			break;
		default:
			p.captured.rotation = p.local.rotation = 0;
			break;
		}

		//----Splatter
		if (splatter) {
			float splattertemp = splatter;
			float splatx = random_generation.Range(-splatter, splatter);
			float splaty = random_generation.Range(-splatter, splatter);

			while (GetDistance(0, 0, splatx, splaty) >= splattertemp && splattertemp > 0) {
				splatx = random_generation.Range(-splatter, splatter);
				splaty = random_generation.Range(-splatter, splatter);
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

		if (angle_setting == AngleSetting::tfxAlign && common.property_flags & tfxEmitterPropertyFlags_edge_traversal)
			p.world.rotation = p.local.rotation = direction + angle_offset;

		bool line = common.property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == EmissionType::tfxLine;

		TransformParticle(p, *this);
		p.captured = p.world;
		p.captured.scale = p.world.scale;

		//----Motion randomness
		p.noise_offset = random_generation.Range(noise_offset_variation) + noise_offset;
		p.noise_resolution = noise_resolution + 0.01f;

		if (!(common.property_flags & tfxEmitterPropertyFlags_edge_traversal) || emission_type != EmissionType::tfxLine) {
			direction = p.emission_angle = GetEmissionDirection(p.local.position, p.world.position, current.emitter_size) + common.library->overtime_graphs[overtime].direction.GetFirstValue();
		}

		//----Normalize Velocity to direction
		tfxVec2 velocity_normal;
		velocity_normal.x = std::sinf(direction);
		velocity_normal.y = -std::cosf(direction);

		//p.velocity = p.velocity_normal * p.base.velocity * p.velocity_scale * UPDATE_TIME;

		if ((angle_setting == AngleSetting::tfxAlign || angle_setting == tfxAlignWithEmission) && !line) {
			p.world.rotation = p.local.rotation = GetVectorAngle(velocity_normal.x, velocity_normal.y) + angle_offset;
			if (common.property_flags & tfxEmitterPropertyFlags_relative_angle)
				p.world.rotation += transform.world.rotation;
			p.captured.rotation = p.world.rotation;
			//Reset the matrix again so that any child particles spawn in the correct place
			if (sub_effects.size()) {
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
		if (common.property_flags & tfxEmitterPropertyFlags_random_start_frame && animation_frames > 1) {
			p.image_frame = random_generation.Range(animation_frames);
		}
		else {
			p.image_frame = start_frame;
		}

		//----Color
		p.color.a = unsigned char(255.f * common.library->overtime_graphs[overtime].opacity.GetFirstValue() * common.root_effect->current.opacity);
		p.intensity = common.library->overtime_graphs[overtime].opacity.GetFirstValue();
		if (common.property_flags & tfxEmitterPropertyFlags_random_color) {
			float age = random_generation.Range(p.max_age);
			p.color.r = unsigned char(255.f * lookup_overtime_callback(common.library->overtime_graphs[overtime].red, age, p.max_age));
			p.color.g = unsigned char(255.f * lookup_overtime_callback(common.library->overtime_graphs[overtime].green, age, p.max_age));
			p.color.b = unsigned char(255.f * lookup_overtime_callback(common.library->overtime_graphs[overtime].blue, age, p.max_age));
		}
		else {
			p.color.r = unsigned char(255.f * common.library->overtime_graphs[overtime].red.GetFirstValue());
			p.color.g = unsigned char(255.f * common.library->overtime_graphs[overtime].green.GetFirstValue());
			p.color.b = unsigned char(255.f * common.library->overtime_graphs[overtime].blue.GetFirstValue());
		}

		sub_effects.reset();
		while (!sub_effects.eob()) {
			tfxEffect &e = sub_effects.next();
			//Todo: add the sub effects into particles
			//e.parent = nullptr;
			//e.parent_particle = &p;
			//e.highest_particle_age = p.max_age;
			//e.current.overal_scale = current.overal_scale;
			//pm->AddEffect(e, !pm->current_ebuff);
		}

		//if (particle_onspawn_callback)
			//particle_onspawn_callback(p);
	}

	bool EffectEmitter::FreeCapacity() {
		return !particles.full() && !GetRootEffect()->sprites[properties.layer].full();
	}
	
	Particle &EffectEmitter::GrabParticle() {
		Particle p;
		particles.push_back(p);
		return particles.back();
	}
	
	bool EffectEmitter::GrabSprite(unsigned int layer) {
		ParticleSprite p;
		if (sprites[layer].push_back(p))
			return true;

		return false;
	}

	void EffectEmitter::BumpSprites() {
		for (int i = 0; i != tfxLAYERS; ++i) {
			if(!sprites[i].reset()) continue;
			int start_index = sprites[i].start_index;
			while (!sprites[i].eob()) {
				ParticleSprite &s = sprites[i].next(start_index);
				if (s.parameters == tfxINVALID)
					sprites[i].bump();
				else
					return;
			}
		}
	}

	void EffectEmitter::InitCPUParticle(Particle &p, float tween) {
		p.flags = tfxParticleFlags_fresh;
		p.parent = this;
		p.next_ptr = &p;

		if (properties.flags & (tfxEmitterPropertyFlags_single | tfxEmitterPropertyFlags_one_shot))
			current.single_shot_done = true;

		//----Properties

		//Set base values-------------------------------

		//----Life
		p.max_age = current.life + random_generation.Range(current.life_variation);
		p.age = 0.f;

		//----Position
		if (properties.emission_type == EmissionType::tfxPoint) {
			if (properties.flags & tfxEmitterPropertyFlags_relative_position)
				p.local.position = -current.emitter_handle;
			else {
				if (properties.flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) {
					p.local.position = InterpolateVec2(tween, transform.captured.position, transform.world.position);
				}
				else {
					tfxVec2 rotvec = transform.matrix.TransformVector(-current.emitter_handle);
					tfxVec2 spawn_position = InterpolateVec2(tween, transform.captured.position, transform.world.position) * transform.world.scale;
					p.local.position = rotvec + spawn_position;
				}
			}
		}
		else if (properties.emission_type == EmissionType::tfxArea) {
			tfxVec2 position = tfxVec2(0.f, 0.f);

			if (properties.flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				if (properties.flags & tfxEmitterPropertyFlags_fill_area) {
					if (!(properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
						current.grid_coords.x--;
						if (current.grid_coords.x < 0.f) {
							current.grid_coords.y--;
							current.grid_coords.x = properties.grid_points.x - 1;
							if (current.grid_coords.y < 0.f)
								current.grid_coords.y = properties.grid_points.y - 1;
						}
					}

					p.local.position = position + (current.grid_coords * current.grid_segment_size) + current.emitter_handle;

					if (properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
						current.grid_coords.x++;
						if (current.grid_coords.x == properties.grid_points.x) {
							current.grid_coords.y++;
							current.grid_coords.x = 0.f;
							if (current.grid_coords.y >= properties.grid_points.y)
								current.grid_coords.y = 0.f;
						}
					}
				}
				else {

					if (properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {

						current.grid_direction.x = 1;
						current.grid_direction.y = 0;
						if (current.grid_coords.x == properties.grid_points.x - 1 && current.grid_coords.y >= 0 && current.grid_coords.y < properties.grid_points.y - 1) {
							current.grid_direction.x = 0;
							current.grid_direction.y = 1;
						}
						else if (current.grid_coords.x > 0 && current.grid_coords.x < properties.grid_points.x && current.grid_coords.y == properties.grid_points.y - 1) {
							current.grid_direction.x = -1;
							current.grid_direction.y = 0;
						}
						else if (current.grid_coords.x == 0 && current.grid_coords.y > 0 && current.grid_coords.y < properties.grid_points.y) {
							current.grid_direction.x = 0;
							current.grid_direction.y = -1;
						}

					}
					else {

						current.grid_direction.x = -1;
						current.grid_direction.y = 0;
						if (current.grid_coords.x == properties.grid_points.x - 1 && current.grid_coords.y > 0 && current.grid_coords.y < properties.grid_points.y) {
							current.grid_direction.x = 0;
							current.grid_direction.y = -1;
						}
						else if (current.grid_coords.x >= 0 && current.grid_coords.x < properties.grid_points.x - 1 && current.grid_coords.y == properties.grid_points.y - 1) {
							current.grid_direction.x = 1;
							current.grid_direction.y = 0;
						}
						else if (current.grid_coords.x == 0 && current.grid_coords.y >= 0 && current.grid_coords.y < properties.grid_points.y - 1) {
							current.grid_direction.x = 0;
							current.grid_direction.y = 1;
						}

					}

					current.grid_coords += current.grid_direction;
					tfxBound(current.grid_coords, properties.grid_points);
					p.local.position = position + (current.grid_coords * current.grid_segment_size) + current.emitter_handle;
				}
			}
			else {
				if (properties.flags & tfxEmitterPropertyFlags_fill_area) {
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

				p.local.position = position + current.emitter_handle;
			}

			//----TForm and Emission
			if (!(properties.flags & tfxEmitterPropertyFlags_relative_position)) {
				p.local.position = transform.matrix.TransformVector(tfxVec2(p.local.position.x, p.local.position.y));
				p.local.position = transform.world.position + p.local.position * transform.world.scale;
			}

		}
		else if (properties.emission_type == EmissionType::tfxEllipse) {
			tfxVec2 emitter_size = (current.emitter_size * .5f);
			tfxVec2 position = tfxVec2(0.f, 0.f);

			if (properties.flags & tfxEmitterPropertyFlags_spawn_on_grid && !(properties.flags & tfxEmitterPropertyFlags_fill_area)) {

				current.grid_coords.y = 0.f;

				if (properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					current.grid_coords.x--;
					if (current.grid_coords.x < 0.f) {
						current.grid_coords.x = properties.grid_points.x - 1;
					}
				}

				float th = current.grid_coords.x * current.grid_segment_size.x + current.arc_offset;
				p.local.position = tfxVec2(std::cosf(th) * emitter_size.x + current.emitter_handle.x + emitter_size.x,
					-std::sinf(th) * emitter_size.y + current.emitter_handle.y + emitter_size.y);

				if (!(properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					current.grid_coords.x++;
					if (current.grid_coords.x >= properties.grid_points.x) {
						current.grid_coords.x = 0.f;
					}
				}

			}
			else if (!(properties.flags & tfxEmitterPropertyFlags_fill_area)) {
				float th = random_generation.Range(current.arc_size) + current.arc_offset;

				p.local.position = tfxVec2(std::cosf(th) * emitter_size.x + current.emitter_handle.x + emitter_size.x,
					-std::sinf(th) * emitter_size.y + current.emitter_handle.y + emitter_size.y);

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
			if (!(properties.flags & tfxEmitterPropertyFlags_relative_position)) {
				p.local.position = transform.matrix.TransformVector(tfxVec2(p.local.position.x, p.local.position.y));
				p.local.position = transform.world.position + p.local.position * transform.world.scale;
			}

		}
		else if (properties.emission_type == EmissionType::tfxLine) {
			if (properties.flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				current.grid_coords.x = 0.f;

				if (!(properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					current.grid_coords.y--;
					if (current.grid_coords.y < 0.f) {
						current.grid_coords.y = properties.grid_points.x - 1;
					}
				}

				p.local.position = tfxVec2(current.grid_coords * -current.grid_segment_size);
				p.local.position += current.emitter_handle;

				if (properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					current.grid_coords.y++;
					if (current.grid_coords.y >= properties.grid_points.x) {
						current.grid_coords.y = 0.f;
					}
				}

			}
			else {
				p.local.position.x = 0.f;
				p.local.position.y = random_generation.Range(-current.emitter_size.y, 0.f);

				p.local.position += current.emitter_handle;

			}

			//----TForm and Emission
			if (!(properties.flags & tfxEmitterPropertyFlags_relative_position) && !(properties.flags & tfxEmitterPropertyFlags_edge_traversal)) {
				p.local.position = transform.matrix.TransformVector(tfxVec2(p.local.position.x, p.local.position.y));
				p.local.position = transform.world.position + p.local.position * transform.world.scale;
			}

		}

		//----Weight
		if (current.weight) {
			p.base.weight = current.weight * library->overtime_graphs[overtime].weight.GetFirstValue();
			if (current.weight_variation > 0) {
				p.base.weight += random_generation.Range(-current.weight_variation, current.weight_variation) * library->overtime_graphs[overtime].weight.GetFirstValue();
			}
		}
		else {
			p.base.weight = 0;
		}
		p.weight_acceleration = p.base.weight * library->overtime_graphs[overtime].weight.GetFirstValue() * UPDATE_TIME;

		//----Velocity
		float velocity_scale = library->overtime_graphs[overtime].velocity.GetFirstValue() * current.velocity_adjuster;
		p.base.velocity = current.velocity + random_generation.Range(-current.velocity_variation, current.velocity_variation);

		//----Size
		if (!(properties.flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			p.base.random_size.x = random_generation.Range(current.size_variation.x);
			p.base.random_size.y = random_generation.Range(current.size_variation.y);
			p.base.size.y = p.base.random_size.y + current.size.y;
			p.base.size.x = (p.base.random_size.x + current.size.x) / properties.image->image_size.x;
			float height = p.base.size.y / properties.image->image_size.y;

			p.world.scale.x = p.base.size.x * library->overtime_graphs[overtime].width.GetFirstValue();

			if (library->overtime_graphs[overtime].stretch.GetFirstValue()) {
				float velocity = std::fabsf(velocity_scale * p.base.velocity) * UPDATE_TIME;
				velocity += p.weight_acceleration * UPDATE_TIME;
				p.world.scale.y = (library->overtime_graphs[overtime].height.GetFirstValue() * parent->current.size.y * (p.base.size.y + (velocity * library->overtime_graphs[overtime].stretch.GetFirstValue() * parent->current.stretch))) / properties.image->image_size.y;
			}
			else {
				if (properties.flags & tfxEmitterPropertyFlags_lifetime_uniform_size) {
					p.world.scale.y = p.world.scale.x;
				}
				else {
					p.world.scale.y = height * library->overtime_graphs[overtime].height.GetFirstValue();
				}
			}
		}
		else {
			p.base.random_size.x = random_generation.Range(current.size_variation.x);
			p.base.random_size.y = p.base.random_size.x;
			p.base.size.y = p.base.random_size.y + current.size.y;
			p.base.size.x = (p.base.random_size.x + current.size.x) / properties.image->image_size.x;
			float height = p.base.size.y / properties.image->image_size.y;

			p.world.scale.x = p.base.size.x * library->overtime_graphs[overtime].width.GetFirstValue();

			if (library->overtime_graphs[overtime].stretch.GetFirstValue()) {
				float velocity = std::fabsf(velocity_scale * p.base.velocity) * UPDATE_TIME;
				velocity += p.weight_acceleration * UPDATE_TIME;
				p.world.scale.y = (library->overtime_graphs[overtime].width.GetFirstValue() * parent->current.size.y * (p.base.size.y + (velocity * library->overtime_graphs[overtime].stretch.GetFirstValue() * parent->current.stretch))) / properties.image->image_size.y;
			}
			else {
				p.world.scale.y = p.world.scale.x;
			}
		}

		//----Spin
		p.base.spin = random_generation.Range(-current.spin_variation, std::abs(current.spin_variation)) + current.spin;

		switch (properties.angle_setting) {
		case AngleSetting::tfxRandom:
			p.captured.rotation = p.local.rotation = random_generation.Range(properties.angle_offset);
			break;
		case AngleSetting::tfxSpecify:
			p.captured.rotation = p.local.rotation = properties.angle_offset;
			break;
		default:
			p.captured.rotation = p.local.rotation = 0;
			break;
		}

		//----Splatter
		if (current.splatter) {
			float splattertemp = current.splatter;
			float splatx = random_generation.Range(-current.splatter, current.splatter);
			float splaty = random_generation.Range(-current.splatter, current.splatter);

			while (GetDistance(0, 0, splatx, splaty) >= splattertemp && splattertemp > 0) {
				splatx = random_generation.Range(-current.splatter, current.splatter);
				splaty = random_generation.Range(-current.splatter, current.splatter);
			}

			if (!(properties.flags & tfxEmitterPropertyFlags_relative_position)) {
				p.local.position.x += splatx * transform.world.scale.x;
				p.local.position.y += splaty * transform.world.scale.y;
			}
			else {
				p.local.position.x += splatx;
				p.local.position.y += splaty;
			}
		}

		float direction = 0;

		if (properties.angle_setting == AngleSetting::tfxAlign && properties.flags & tfxEmitterPropertyFlags_edge_traversal)
			p.world.rotation = p.local.rotation = direction + properties.angle_offset;

		bool line = properties.flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type == EmissionType::tfxLine;

		TransformParticle(p, *this);
		p.captured = p.world;
		p.captured.scale = p.world.scale;

		//----Motion randomness
		p.noise_offset = random_generation.Range(current.noise_offset_variation) + current.noise_offset;
		p.noise_resolution = current.noise_resolution + 0.01f;

		if (!(properties.flags & tfxEmitterPropertyFlags_edge_traversal) || properties.emission_type != EmissionType::tfxLine) {
			direction = p.emission_angle = GetEmissionDirection(p.local.position, p.world.position) + library->overtime_graphs[overtime].direction.GetFirstValue();
		}

		//----Normalize Velocity to direction
		tfxVec2 velocity_normal;
		velocity_normal.x = std::sinf(direction);
		velocity_normal.y = -std::cosf(direction);

		//p.velocity = p.velocity_normal * p.base.velocity * p.velocity_scale * UPDATE_TIME;

		if ((properties.angle_setting == AngleSetting::tfxAlign || properties.angle_setting == tfxAlignWithEmission) && !line) {
			p.world.rotation = p.local.rotation = GetVectorAngle(velocity_normal.x, velocity_normal.y) + properties.angle_offset;
			if (properties.flags & tfxEmitterPropertyFlags_relative_angle)
				p.world.rotation += transform.world.rotation;
			p.captured.rotation = p.world.rotation;
			//Reset the matrix again so that any child particles spawn in the correct place
			if (sub_effectors.size()) {
				float s = sin(p.local.rotation);
				float c = cos(p.local.rotation);
				p.matrix.Set(c, s, -s, c);
			}
		}

		//----Handle
		/*if (properties.flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
			p.handle = tfxVec2(0.5f, 0.5f);
		}
		else {
			p.handle = properties.image_handle;
		}*/

		//----Image
		//p.image = properties.image;
		if (properties.flags & tfxEmitterPropertyFlags_random_start_frame && properties.image->animation_frames > 1) {
			p.image_frame = random_generation.Range(properties.image->animation_frames);
		}
		else {
			p.image_frame = properties.start_frame;
		}

		//----Color
		p.color.a = unsigned char(255.f * library->overtime_graphs[overtime].opacity.GetFirstValue() * parent->current.color.a);
		p.intensity = library->overtime_graphs[overtime].opacity.GetFirstValue();
		if (properties.flags & tfxEmitterPropertyFlags_random_color) {
			float age = random_generation.Range(p.max_age);
			p.color.r = unsigned char(255.f * lookup_overtime_callback(library->overtime_graphs[overtime].red, age, p.max_age));
			p.color.g = unsigned char(255.f * lookup_overtime_callback(library->overtime_graphs[overtime].green, age, p.max_age));
			p.color.b = unsigned char(255.f * lookup_overtime_callback(library->overtime_graphs[overtime].blue, age, p.max_age));
		}
		else {
			p.color.r = unsigned char(255.f * library->overtime_graphs[overtime].red.GetFirstValue());
			p.color.g = unsigned char(255.f * library->overtime_graphs[overtime].green.GetFirstValue());
			p.color.b = unsigned char(255.f * library->overtime_graphs[overtime].blue.GetFirstValue());
		}

		if (sub_effectors.size()) {
			for (auto &e : sub_effectors) {
				if (!pm->FreeEffectCapacity())
					break;
				e.parent = nullptr;
				e.parent_particle = &p;
				e.highest_particle_age = p.max_age;
				e.current.overal_scale = current.overal_scale;
				pm->AddEffect(e, !pm->current_ebuff);
			}
		}

		if (particle_onspawn_callback)
			particle_onspawn_callback(p);
	}

	void EffectEmitter::InitComputeParticle(ComputeParticle &p, float tween) {
		if (properties.flags & (tfxEmitterPropertyFlags_single | tfxEmitterPropertyFlags_one_shot))
			current.single_shot_done = true;

		//----Properties

		//Set base values-------------------------------

		//----Life
		p.max_age = current.life + random_generation.Range(current.life_variation);
		p.age = 0.f;

		//----Position
		if (properties.emission_type == EmissionType::tfxPoint) {
			if (properties.flags & tfxEmitterPropertyFlags_relative_position)
				p.local_position = -current.emitter_handle;
			else {
				if (properties.flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) {
					p.local_position = InterpolateVec2(tween, transform.captured.position, transform.world.position);
				}
				else {
					tfxVec2 rotvec = transform.matrix.TransformVector(-current.emitter_handle);
					tfxVec2 spawn_position = InterpolateVec2(tween, transform.captured.position, transform.world.position) * transform.world.scale;
					p.local_position = rotvec + spawn_position;
				}
			}
		}
		else if (properties.emission_type == EmissionType::tfxArea) {
			tfxVec2 position = tfxVec2(0.f, 0.f);

			if (properties.flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				if (properties.flags & tfxEmitterPropertyFlags_fill_area) {
					if (!(properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
						current.grid_coords.x--;
						if (current.grid_coords.x < 0.f) {
							current.grid_coords.y--;
							current.grid_coords.x = properties.grid_points.x - 1;
							if (current.grid_coords.y < 0.f)
								current.grid_coords.y = properties.grid_points.y - 1;
						}
					}

					p.local_position = position + (current.grid_coords * current.grid_segment_size) + current.emitter_handle;

					if (properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
						current.grid_coords.x++;
						if (current.grid_coords.x == properties.grid_points.x) {
							current.grid_coords.y++;
							current.grid_coords.x = 0.f;
							if (current.grid_coords.y >= properties.grid_points.y)
								current.grid_coords.y = 0.f;
						}
					}
				}
				else {

					if (properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {

						current.grid_direction.x = 1;
						current.grid_direction.y = 0;
						if (current.grid_coords.x == properties.grid_points.x - 1 && current.grid_coords.y >= 0 && current.grid_coords.y < properties.grid_points.y - 1) {
							current.grid_direction.x = 0;
							current.grid_direction.y = 1;
						}
						else if (current.grid_coords.x > 0 && current.grid_coords.x < properties.grid_points.x && current.grid_coords.y == properties.grid_points.y - 1) {
							current.grid_direction.x = -1;
							current.grid_direction.y = 0;
						}
						else if (current.grid_coords.x == 0 && current.grid_coords.y > 0 && current.grid_coords.y < properties.grid_points.y) {
							current.grid_direction.x = 0;
							current.grid_direction.y = -1;
						}

					}
					else {

						current.grid_direction.x = -1;
						current.grid_direction.y = 0;
						if (current.grid_coords.x == properties.grid_points.x - 1 && current.grid_coords.y > 0 && current.grid_coords.y < properties.grid_points.y) {
							current.grid_direction.x = 0;
							current.grid_direction.y = -1;
						}
						else if (current.grid_coords.x >= 0 && current.grid_coords.x < properties.grid_points.x - 1 && current.grid_coords.y == properties.grid_points.y - 1) {
							current.grid_direction.x = 1;
							current.grid_direction.y = 0;
						}
						else if (current.grid_coords.x == 0 && current.grid_coords.y >= 0 && current.grid_coords.y < properties.grid_points.y - 1) {
							current.grid_direction.x = 0;
							current.grid_direction.y = 1;
						}

					}

					current.grid_coords += current.grid_direction;
					tfxBound(current.grid_coords, properties.grid_points);
					p.local_position = position + (current.grid_coords * current.grid_segment_size) + current.emitter_handle;
				}
			}
			else {
				if (properties.flags & tfxEmitterPropertyFlags_fill_area) {
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

				p.local_position = position + current.emitter_handle;
			}

			//----TForm and Emission
			if (!(properties.flags & tfxEmitterPropertyFlags_relative_position)) {
				p.local_position = transform.matrix.TransformVector(tfxVec2(p.local_position.x, p.local_position.y));
				p.local_position = transform.world.position + p.local_position * transform.world.scale;
			}

		}
		else if (properties.emission_type == EmissionType::tfxEllipse) {
			tfxVec2 emitter_size = (current.emitter_size * .5f);
			tfxVec2 position = tfxVec2(0.f, 0.f);

			if (properties.flags & tfxEmitterPropertyFlags_spawn_on_grid && !(properties.flags & tfxEmitterPropertyFlags_fill_area)) {

				current.grid_coords.y = 0.f;

				if (properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					current.grid_coords.x--;
					if (current.grid_coords.x < 0.f) {
						current.grid_coords.x = properties.grid_points.x - 1;
					}
				}

				float th = current.grid_coords.x * current.grid_segment_size.x + current.arc_offset;
				p.local_position = tfxVec2(std::cosf(th) * emitter_size.x + current.emitter_handle.x + emitter_size.x,
					-std::sinf(th) * emitter_size.y + current.emitter_handle.y + emitter_size.y);

				if (!(properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					current.grid_coords.x++;
					if (current.grid_coords.x >= properties.grid_points.x) {
						current.grid_coords.x = 0.f;
					}
				}

			}
			else if (!(properties.flags & tfxEmitterPropertyFlags_fill_area)) {
				float th = random_generation.Range(current.arc_size) + current.arc_offset;

				p.local_position = tfxVec2(std::cosf(th) * emitter_size.x + current.emitter_handle.x + emitter_size.x,
					-std::sinf(th) * emitter_size.y + current.emitter_handle.y + emitter_size.y);

			}
			else {
				p.local_position.x = random_generation.Range(-emitter_size.x, emitter_size.x);
				p.local_position.y = random_generation.Range(-emitter_size.y, emitter_size.y);

				while ((std::pow(p.local_position.x, 2) / std::pow(emitter_size.x, 2)) + (std::pow(p.local_position.y, 2) / std::pow(emitter_size.y, 2)) > 1) {
					p.local_position.x = random_generation.Range(-emitter_size.x, emitter_size.x);
					p.local_position.y = random_generation.Range(-emitter_size.y, emitter_size.y);
				}
			}

			//----TForm and Emission
			if (!(properties.flags & tfxEmitterPropertyFlags_relative_position)) {
				p.local_position = transform.matrix.TransformVector(tfxVec2(p.local_position.x, p.local_position.y));
				p.local_position = transform.world.position + p.local_position * transform.world.scale;
			}

		}
		else if (properties.emission_type == EmissionType::tfxLine) {
			if (properties.flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				current.grid_coords.x = 0.f;

				if (!(properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					current.grid_coords.y--;
					if (current.grid_coords.y < 0.f) {
						current.grid_coords.y = properties.grid_points.x - 1;
					}
				}

				p.local_position = tfxVec2(current.grid_coords * -current.grid_segment_size);
				p.local_position += current.emitter_handle;

				if (properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					current.grid_coords.y++;
					if (current.grid_coords.y >= properties.grid_points.x) {
						current.grid_coords.y = 0.f;
					}
				}

			}
			else {
				p.local_position.x = 0.f;
				p.local_position.y = random_generation.Range(-current.emitter_size.y, 0.f);

				p.local_position += current.emitter_handle;

			}

			//----TForm and Emission
			if (!(properties.flags & tfxEmitterPropertyFlags_relative_position) && !(properties.flags & tfxEmitterPropertyFlags_edge_traversal)) {
				p.local_position = transform.matrix.TransformVector(tfxVec2(p.local_position.x, p.local_position.y));
				p.local_position = transform.world.position + p.local_position * transform.world.scale;
			}

		}

		//----Weight
		if (current.weight) {
			p.base_weight = current.weight * library->overtime_graphs[overtime].weight.GetFirstValue();
			if (current.weight_variation > 0) {
				p.base_weight += random_generation.Range(-current.weight_variation, current.weight_variation) * library->overtime_graphs[overtime].weight.GetFirstValue();
			}
		}
		else {
			p.base_weight = 0;
		}
		p.weight_acceleration = p.base_weight * library->overtime_graphs[overtime].weight.GetFirstValue() * UPDATE_TIME;

		//----Velocity
		float velocity_scale = library->overtime_graphs[overtime].velocity.GetFirstValue() * current.velocity_adjuster;
		p.base_velocity = current.velocity + random_generation.Range(-current.velocity_variation, current.velocity_variation);

		//----Size
		if (!(properties.flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			float base_random_size_x = random_generation.Range(current.size_variation.x);
			float base_random_size_y = random_generation.Range(current.size_variation.y);
			p.base_size.y = base_random_size_y + current.size.y;
			p.base_size.x = (base_random_size_x + current.size.x) / properties.image->image_size.x;
			//p.base_size.y = p.base_height / properties.image->image_size.y;

			//p.scale_rotation.x = p.base_size.x * library->overtime_graphs[overtime].width.GetFirstValue();

			if (library->overtime_graphs[overtime].stretch.GetFirstValue()) {
				float velocity = std::fabsf(velocity_scale * p.base_velocity) * UPDATE_TIME;
				velocity += p.weight_acceleration * UPDATE_TIME;
				//p.scale_rotation.y = (library->overtime_graphs[overtime].height.GetFirstValue() * parent->current.size.y * (p.base_height + (velocity * library->overtime_graphs[overtime].stretch.GetFirstValue() * parent->current.stretch))) / properties.image->image_size.y;
			}
			else {
				if (properties.flags & tfxEmitterPropertyFlags_lifetime_uniform_size) {
					//p.scale_rotation.y = p.base_size.y * library->overtime_graphs[overtime].width.GetFirstValue();
				}
				else {
					//p.scale_rotation.y = p.base_size.y * library->overtime_graphs[overtime].height.GetFirstValue();
				}
			}
		}
		else {
			float base_random_size_x = random_generation.Range(current.size_variation.x);
			float base_random_size_y = base_random_size_x;
			p.base_size.y = base_random_size_y + current.size.y;
			p.base_size.x = (base_random_size_x + current.size.x) / properties.image->image_size.x;
			//p.base_size.y = p.base_height / properties.image->image_size.y;

			//p.scale_rotation.x = p.base_size.x * library->overtime_graphs[overtime].width.GetFirstValue();

			if (library->overtime_graphs[overtime].stretch.GetFirstValue()) {
				float velocity = std::fabsf(velocity_scale * p.base_velocity) * UPDATE_TIME;
				velocity += p.weight_acceleration * UPDATE_TIME;
				//p.scale_rotation.y = (library->overtime_graphs[overtime].width.GetFirstValue() * parent->current.size.y * (p.base_height + (velocity * library->overtime_graphs[overtime].stretch.GetFirstValue() * parent->current.stretch))) / properties.image->image_size.y;
			}
			else {
				//p.scale_rotation.y = p.scale_rotation.x;
			}
		}

		//----Spin
		p.base_spin = random_generation.Range(-current.spin_variation, std::abs(current.spin_variation)) + current.spin;

		switch (properties.angle_setting) {
		case AngleSetting::tfxRandom:
			p.local_rotation = random_generation.Range(properties.angle_offset);
			break;
		case AngleSetting::tfxSpecify:
			p.local_rotation = properties.angle_offset;
			break;
		default:
			p.local_rotation = 0;
			break;
		}

		//----Splatter
		if (current.splatter) {
			float splattertemp = current.splatter;
			float splatx = random_generation.Range(-current.splatter, current.splatter);
			float splaty = random_generation.Range(-current.splatter, current.splatter);

			while (GetDistance(0, 0, splatx, splaty) >= splattertemp && splattertemp > 0) {
				splatx = random_generation.Range(-current.splatter, current.splatter);
				splaty = random_generation.Range(-current.splatter, current.splatter);
			}

			if (!(properties.flags & tfxEmitterPropertyFlags_relative_position)) {
				p.local_position.x += splatx * transform.world.scale.x;
				p.local_position.y += splaty * transform.world.scale.y;
			}
			else {
				p.local_position.x += splatx;
				p.local_position.y += splaty;
			}
		}

		float direction = 0;

		//if (properties.angle_setting == AngleSetting::tfxAlign && properties.flags & tfxEmitterPropertyFlags_edge_traversal)
			//p.scale_rotation.w = p.scale_rotation.z = direction + properties.angle_offset;

		bool line = properties.flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type == EmissionType::tfxLine;

		FormState local;
		local.position = p.local_position;
		local.rotation = p.local_rotation;
		FormState world;
		Transform(local, world, *this);

		//----Motion randomness
		p.noise_offset = random_generation.Range(current.noise_offset_variation) + current.noise_offset;
		p.noise_resolution = current.noise_resolution + 0.01f;

		if (!(properties.flags & tfxEmitterPropertyFlags_edge_traversal) || properties.emission_type != EmissionType::tfxLine) {
			direction = p.emission_angle = GetEmissionDirection(p.local_position, world.position) + library->overtime_graphs[overtime].direction.GetFirstValue();
		}

		//----Normalize Velocity to direction
		tfxVec2 velocity_normal;
		velocity_normal.x = std::sinf(direction);
		velocity_normal.y = -std::cosf(direction);

		//p.velocity = p.velocity_normal * p.base_velocity * p.velocity_scale * UPDATE_TIME;

		if ((properties.angle_setting == AngleSetting::tfxAlign || properties.angle_setting == tfxAlignWithEmission) && !line) {
			p.local_rotation = GetVectorAngle(velocity_normal.x, velocity_normal.y) + properties.angle_offset;
			//if (properties.flags & tfxEmitterPropertyFlags_relative_angle)
				//p.scale_rotation.w += world.rotation;
		}

		//----Handle
		/*if (properties.flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
			p.handle = tfxVec2(0.5f, 0.5f);
		}
		else {
			p.handle = properties.image_handle;
		}*/

		//----Image
		//p.image = properties.image;
		if (properties.flags & tfxEmitterPropertyFlags_random_start_frame && properties.image->animation_frames > 1) {
			p.image_frame = random_generation.Range(properties.image->animation_frames);
		}
		else {
			p.image_frame = properties.start_frame;
		}

		//----Color
		/*p.color.a = unsigned char(255.f * library->overtime_graphs[overtime].opacity.GetFirstValue() * parent->current.color.a);
		p.intensity = library->overtime_graphs[overtime].opacity.GetFirstValue();
		if (properties.flags & tfxEmitterPropertyFlags_random_color) {
			float age = random_generation.Range(p.max_age);
			p.color.r = unsigned char(255.f * lookup_overtime_callback(library->overtime_graphs[overtime].red, age, p.max_age));
			p.color.g = unsigned char(255.f * lookup_overtime_callback(library->overtime_graphs[overtime].green, age, p.max_age));
			p.color.b = unsigned char(255.f * lookup_overtime_callback(library->overtime_graphs[overtime].blue, age, p.max_age));
		}
		else {
			p.color.r = unsigned char(255.f * library->overtime_graphs[overtime].red.GetFirstValue());
			p.color.g = unsigned char(255.f * library->overtime_graphs[overtime].green.GetFirstValue());
			p.color.b = unsigned char(255.f * library->overtime_graphs[overtime].blue.GetFirstValue());
		}*/

		p.control_slot_and_layer = (properties.layer << 24) + compute_slot_id;

	}

	float EffectEmitter::GetEmissionDirection(tfxVec2 &local_position, tfxVec2 &world_position) {

		//----Emission
		float range = current.emission_angle_variation *.5f;
		float direction = 0;

		if (properties.emission_type == EmissionType::tfxPoint)
			return direction + current.emission_angle + random_generation.Range(-range, range);

		tfxVec2 tmp_position;
		if (local_position.x == 0 && local_position.y == 0)
			tmp_position = current.emitter_size;
		else
			tmp_position = local_position;

		if (properties.emission_direction == EmissionDirection::tfxOutwards) {

			tfxVec2 to_handle;

			if (properties.flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = (tmp_position);
			else
				to_handle = (world_position - transform.world.position);

			direction = GetVectorAngle(to_handle.x, to_handle.y);

		}
		else if (properties.emission_direction == EmissionDirection::tfxInwards) {

			tfxVec2 to_handle;

			if (properties.flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = (-tmp_position);
			else
				to_handle = (transform.world.position - world_position);

			direction = GetVectorAngle(to_handle.x, to_handle.y);

		}
		else if (properties.emission_direction == EmissionDirection::tfxBothways) {

			if (current.emission_alternator) {

				tfxVec2 to_handle;

				if (properties.flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = (tmp_position);
				else
					to_handle = (world_position - transform.world.position);

				direction = GetVectorAngle(to_handle.x, to_handle.y);

			}
			else {

				tfxVec2 to_handle;

				if (properties.flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = (-tmp_position);
				else
					to_handle = (transform.world.position - world_position);

				direction = GetVectorAngle(to_handle.x, to_handle.y);

			}

			current.emission_alternator = !current.emission_alternator;
		}

		if (std::isnan(direction))
			direction = 0.f;
		return direction + current.emission_angle + random_generation.Range(-range, range);
	}

	void EffectEmitter::UpdateEmitterState() {
		EffectEmitter &e = *parent;
		current.amount = lookup_callback(library->base_graphs[base].amount, current.frame);
		current.amount_variation = lookup_callback(library->variation_graphs[base].amount, current.frame);
		current.life = lookup_callback(library->base_graphs[base].life, current.frame) * e.current.life;
		current.life_variation = lookup_callback(library->variation_graphs[variation].life, current.frame) * e.current.life;
		if (!(properties.flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			current.size.x = lookup_callback(library->base_graphs[base].width, current.frame) * e.current.size.x;
			current.size.y = lookup_callback(library->base_graphs[base].height, current.frame) * e.current.size.y;
		}
		else {
			current.size.x = lookup_callback(library->base_graphs[base].width, current.frame);
			if (e.properties.flags & tfxEmitterPropertyFlags_global_uniform_size)
				current.size.y = current.size.x * e.current.size.x;
			else
				current.size.y = current.size.x * e.current.size.y;
			current.size.x *= e.current.size.x;
		}
		current.size_variation.x = lookup_callback(library->variation_graphs[variation].width, current.frame) * e.current.size.x;
		current.size_variation.y = lookup_callback(library->variation_graphs[variation].height, current.frame) * e.current.size.y;
		current.velocity = lookup_callback(library->base_graphs[base].velocity, current.frame) * e.current.velocity;
		current.velocity_variation = lookup_callback(library->variation_graphs[variation].velocity, current.frame) * e.current.velocity;
		current.velocity_adjuster = lookup_callback(library->overtime_graphs[overtime].velocity_adjuster, current.frame);
		current.spin = lookup_callback(library->base_graphs[base].spin, current.frame) * e.current.spin;
		current.spin_variation = lookup_callback(library->variation_graphs[variation].spin, current.frame) * e.current.spin;
		transform.local.rotation = lookup_callback(library->property_graphs[property].emitter_angle, current.frame);
		current.emission_angle = lookup_callback(library->property_graphs[property].emission_angle, current.frame);
		current.emission_angle_variation = lookup_callback(library->property_graphs[property].emission_range, current.frame);
		current.color.r = library->overtime_graphs[overtime].red.GetFirstValue();
		current.color.g = library->overtime_graphs[overtime].green.GetFirstValue();
		current.color.b = library->overtime_graphs[overtime].blue.GetFirstValue();
		current.color.a = e.current.color.a;
		current.splatter = lookup_callback(library->property_graphs[property].splatter, current.frame) * e.current.splatter;
		current.emitter_size.y = lookup_callback(library->property_graphs[property].emitter_height, current.frame);
		current.weight = lookup_callback(library->base_graphs[base].weight, current.frame) * e.current.weight;
		current.weight_variation = lookup_callback(library->variation_graphs[variation].weight, current.frame) * e.current.weight;
		current.noise_offset_variation = lookup_callback(library->variation_graphs[variation].noise_offset, current.frame);
		current.noise_offset = lookup_callback(library->base_graphs[variation].noise_offset, current.frame);
		current.noise_resolution = lookup_callback(library->variation_graphs[variation].noise_resolution, current.frame);
		current.stretch = e.current.stretch;
		current.overal_scale = e.current.overal_scale;
		transform.world.scale = e.transform.world.scale;

		//----Handle
		if (properties.flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
			current.image_handle = tfxVec2(0.5f, 0.5f);
		}
		else {
			current.image_handle = properties.image_handle;
		}

		if (properties.emission_type == EmissionType::tfxArea || properties.emission_type == EmissionType::tfxEllipse) {
			current.emitter_size.x = lookup_callback(library->property_graphs[property].emitter_width, current.frame);
		}
		else
			current.emitter_size.x = 0.f;

		if (properties.emission_type == EmissionType::tfxEllipse) {
			current.arc_size = lookup_callback(library->property_graphs[property].arc_size, current.frame);
			current.arc_offset = lookup_callback(library->property_graphs[property].arc_offset, current.frame);
		}

		if (properties.flags & tfxEmitterPropertyFlags_emitter_handle_auto_center && properties.emission_type != EmissionType::tfxPoint) {
			current.emitter_handle = current.emitter_size * -0.5f;
		}
		else if (!(properties.flags & tfxEmitterPropertyFlags_emitter_handle_auto_center)) {
			current.emitter_handle = properties.emitter_handle;
		}

		else {
			current.emitter_handle.x = current.emitter_handle.y = 0;
		}

		if (properties.flags & tfxEmitterPropertyFlags_emitter_handle_auto_center && properties.emission_type == EmissionType::tfxLine) {
			current.emitter_handle = current.emitter_size * 0.5f;
		}

		if (properties.flags & tfxEmitterPropertyFlags_spawn_on_grid) {
			if (properties.emission_type == EmissionType::tfxArea) {
				if (properties.grid_points.x > 1)
					current.grid_segment_size.x = current.emitter_size.x / (properties.grid_points.x - 1);
				if (properties.grid_points.y > 1)
					current.grid_segment_size.y = current.emitter_size.y / (properties.grid_points.y - 1);
			}
			else if (properties.emission_type == EmissionType::tfxEllipse) {
				if (properties.grid_points.x > 0)
					current.grid_segment_size.x = current.arc_size / (properties.grid_points.x);
			}
			else if (properties.emission_type == EmissionType::tfxLine) {
				if (properties.grid_points.x > 1)
					current.grid_segment_size.y = current.emitter_size.y / (properties.grid_points.x - 1);
			}
		}

		if (update_callback)
			update_callback(*this);

	}

	void EffectEmitter::UpdateComputeController() {
		ComputeController &c = *(static_cast<ComputeController*>(pm->compute_controller_ptr) + compute_slot_id);
		c.position = transform.world.position;
		c.scale_rotation.x = transform.world.scale.x;
		c.scale_rotation.y = transform.world.scale.y;
		c.scale_rotation.z = transform.world.rotation;
		c.scale_rotation.w = current.velocity_adjuster;
		c.line_length = current.emitter_size.y;
		c.angle_offset = properties.angle_offset;
		c.end_frame = properties.end_frame;
		int age_rate = (properties.flags & tfxEmitterPropertyFlags_single && !(properties.flags & tfxEmitterPropertyFlags_one_shot)) ? 0 : 1;
		int line_negator = (properties.flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type == tfxLine) ? 0 : 1;
		int spin_negator = (properties.angle_setting == AngleSetting::tfxAlign || properties.flags & tfxEmitterPropertyFlags_relative_angle) ? 0 : 1;
		int position_negator = (properties.flags & tfxEmitterPropertyFlags_relative_position) ? 1 : 0;
		int additive = (properties.blend_mode == tfxAdditive ? 1 : 0);
		c.normalised_values = (age_rate << 31) + (line_negator << 30) + (spin_negator << 29) + (position_negator << 28) + (additive << 27) + (properties.layer << 24) + (int(current.color.a * 255) << 16) + (unsigned short)lookup_value_index;
		//std::cout << ((((c.normalised_values & (u32(0x08000000)) >> 27 ) + 1) << 24) + 2) << std::endl;
		//std::cout << (((c.normalised_values & (u32(0x08000000))) >> 27) + 1) << std::endl;
		c.flags = properties.compute_flags;
		c.image_handle = current.image_handle;
		c.emitter_handle = current.emitter_handle;
		c.stretch = current.stretch;
		c.noise_offset = current.noise_offset;
		c.noise_resolution = current.noise_resolution;
		c.image_data_index = properties.image->compute_shape_index;
		c.frame_rate = properties.image->animation_frames > 1 && properties.flags & tfxEmitterPropertyFlags_animate ? properties.frame_rate : 0.f;
	}

	void EffectEmitter::UpdateEffectState() {
		//If this effect is a sub effect then the graph index will reference the global graphs for the root parent effect
		current.life = lookup_callback(library->global_graphs[global].life, current.frame);
		current.amount = lookup_callback(library->global_graphs[global].amount, current.frame);
		if (!(properties.flags & tfxEmitterPropertyFlags_global_uniform_size)) {
			current.size.x = lookup_callback(library->global_graphs[global].width, current.frame);
			current.size.y = lookup_callback(library->global_graphs[global].height, current.frame);
		}
		else {
			current.size.x = lookup_callback(library->global_graphs[global].width, current.frame);
			current.size.y = current.size.x;
		}
		current.velocity = lookup_callback(library->global_graphs[global].velocity, current.frame);
		current.spin = lookup_callback(library->global_graphs[global].spin, current.frame);
		current.color.a = lookup_callback(library->global_graphs[global].opacity, current.frame);
		current.splatter = lookup_callback(library->global_graphs[global].splatter, current.frame);
		//We don't want to scale twice when the sub effect is transformed, so the values here are set to 1. That means that the root effect will only control the global scale.
		current.overal_scale = lookup_callback(library->global_graphs[global].overal_scale, current.frame);
		if (!parent_particle) {
			transform.world.scale.x = current.overal_scale;
			transform.world.scale.y = current.overal_scale;
			transform.local.rotation = lookup_callback(library->global_graphs[global].effect_angle, current.frame);
		}
		else {
			transform.world.scale.x = current.overal_scale;
			transform.world.scale.y = current.overal_scale;
			transform.local.rotation = 0.f;
		}
		current.stretch = lookup_callback(library->global_graphs[global].stretch, current.frame);
		current.weight = lookup_callback(library->global_graphs[global].weight, current.frame);

		if (update_callback)
			update_callback(*this);
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

	bool ControlParticle(Particle &p, EffectEmitter &e) {
		if (e.pm->update_base_values)
			ReloadBaseValues(p, e);

		p.age += FRAME_LENGTH;

		//-------------------------------------------------------
		//Controll what the particle does over the course of
		//it's lifetime
		//-------------------------------------------------------

		//Before we do anything, see if the particle should be removed/end of life
		p.flags |= e.flags & tfxParticleFlags_remove;
		if (p.flags & tfxParticleFlags_remove)
			return false;

		if (p.age >= p.max_age) {
			if (e.properties.flags & tfxEmitterPropertyFlags_single && !(e.properties.flags & tfxEmitterPropertyFlags_one_shot) && !e.pm->disable_spawing)
				p.age = 0.f;
			else {
				return false;
			}
		}

		float lookup_velocity = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].velocity, p.age, p.max_age);
		float lookup_width = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].width, p.age, p.max_age);
		float lookup_height = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].height, p.age, p.max_age);
		float lookup_weight = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].weight, p.age, p.max_age);
		float lookup_spin = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].spin, p.age, p.max_age);
		float lookup_stretch = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].stretch, p.age, p.max_age);
		float lookup_red = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].red, p.age, p.max_age);
		float lookup_green = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].green, p.age, p.max_age);
		float lookup_blue = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].blue, p.age, p.max_age);
		float lookup_opacity = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].opacity, p.age, p.max_age);
		float lookup_velocity_turbulance = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].velocity_turbulance, p.age, p.max_age);
		float lookup_direction_turbulance = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].direction_turbulance, p.age, p.max_age);
		float lookup_intensity = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].intensity, p.age, p.max_age);
		float lookup_direction = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].direction, p.age, p.max_age);
		float lookup_noise_resolution = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].noise_resolution, p.age, p.max_age);

		float direction = 0;
		float mr_angle = 0;
		float mr_speed = 0;

		tfxVec2 mr_vec;
		if ((e.properties.emission_type != tfxLine && !(e.properties.flags & tfxEmitterPropertyFlags_edge_traversal))
			|| e.properties.emission_type == tfxLine && !(e.properties.flags & tfxEmitterPropertyFlags_edge_traversal)) {
			//----Motion randomness
			if (lookup_velocity_turbulance + lookup_velocity_turbulance) {
				float noise_resolution = lookup_noise_resolution * p.noise_resolution;
				float noise = SimplexNoise::noise(p.local.position.x / noise_resolution + p.noise_offset, p.local.position.y / noise_resolution + p.noise_offset);
				mr_speed = noise * lookup_velocity_turbulance;
				mr_angle = noise * 3.14159265f * 2 * lookup_direction_turbulance;
				mr_vec.x = std::sinf(mr_angle) * mr_speed;
				mr_vec.y = -std::cosf(mr_angle) * mr_speed;
			}

			direction = p.emission_angle + lookup_direction;
		}

		//----Weight Changes
		p.weight_acceleration += p.base.weight * lookup_weight * UPDATE_TIME;

		//----Velocity Changes
		tfxVec2 velocity_normal;
		velocity_normal.x = std::sinf(direction);
		velocity_normal.y = -std::cosf(direction);
		float velocity_scale = lookup_velocity * e.current.velocity_adjuster;

		tfxVec2 current_velocity = ((p.base.velocity * velocity_scale)) * velocity_normal * UPDATE_TIME;
		current_velocity += mr_vec * UPDATE_TIME;
		current_velocity.y += p.weight_acceleration * UPDATE_TIME;

		//----Color changes
		p.color.a = unsigned char(255.f * lookup_opacity * e.parent->current.color.a);
		p.intensity = lookup_intensity;
		if (!(e.properties.flags & tfxEmitterPropertyFlags_random_color)) {
			p.color.r = unsigned char(255.f * lookup_red);
			p.color.g = unsigned char(255.f * lookup_green);
			p.color.b = unsigned char(255.f * lookup_blue);
		}

		p.color = tfxRGBA8(p.color.r, p.color.g, p.color.b, p.color.a);

		//----Size Changes
		tfxVec2 scale;
		float width_overtime = lookup_width;
		scale.x = p.base.size.x * width_overtime;
		//Just here to test:
		//float test1 = p.base.size.x * e.library->LookupPreciseOvertimeNodeList(tfxOvertime_width, e.lookup_node_index, p.age, p.max_age);
		//p.local.scale.x = p.base.size.x * e.library->LookupFastOvertimeValueList(tfxOvertime_width, e.lookup_value_index, p.age, p.max_age);
		if (scale.x < 0.f)
			scale.x = scale.x;

		//----Stretch Changes
		float stretch = lookup_stretch;
		float velocity = std::fabsf(velocity_scale * p.base.velocity + mr_speed + p.weight_acceleration);
		if (e.properties.flags & tfxEmitterPropertyFlags_lifetime_uniform_size) {
			scale.y = (width_overtime * (p.base.size.y + (velocity * stretch * e.current.stretch))) / e.properties.image->image_size.y;
			if (scale.y < scale.x)
				scale.y = scale.x;
		}
		else
			scale.y = (lookup_height * (p.base.size.y + (velocity * stretch * e.current.stretch))) / e.properties.image->image_size.y;

		//----Spin and angle Changes
		float spin = 0;
		if (e.properties.angle_setting != AngleSetting::tfxAlign && !(e.properties.flags & tfxEmitterPropertyFlags_relative_angle)) {
			spin = lookup_spin;
			spin *= p.base.spin;
		}

		//---------------
		//Now that the latest changes are applied, affect the particle state
		//---------------

		//Before applying the behaviour to the particle position, scale and rotation, you have the chance to override them here
		if (e.particle_update_callback)
			e.particle_update_callback(p);

		//----Rotation
		p.local.rotation += spin * UPDATE_TIME;
		if (e.properties.angle_setting == AngleSetting::tfxAlign) {
			tfxVec2 vd = current_velocity.IsNill() ? velocity_normal : current_velocity;
			p.local.rotation = GetVectorAngle(vd.x, vd.y) + e.properties.angle_offset;
		}

		//----Position
		p.local.position += current_velocity * e.current.overal_scale;

		//----Scale
		p.world.scale = scale;

		//Lines - reposition if the particle is travelling along a line
		tfxVec2 offset = velocity_normal * e.current.emitter_size.y;
		float length = std::fabsf(p.local.position.y - e.current.emitter_handle.y);
		float emitter_length = e.current.emitter_size.y;
		bool is_line = e.properties.emission_type == tfxLine && e.properties.flags & tfxEmitterPropertyFlags_edge_traversal;
		bool line_and_loop = is_line && e.properties.end_behaviour == tfxLoop && length > emitter_length;
		bool line_and_kill = is_line && e.properties.end_behaviour == tfxKill && length > emitter_length;
		if (line_and_loop) {
			p.local.position.y -= offset.y;
			p.flags |= tfxParticleFlags_capture_after_transform;
		}
		else if (line_and_kill) {
			return false;
		}

		//----Image animation
		if (e.properties.flags & tfxEmitterPropertyFlags_animate) {
			if (e.properties.flags & tfxEmitterPropertyFlags_reverse_animation)
				p.image_frame -= e.properties.frame_rate * UPDATE_TIME;
			else
				p.image_frame += e.properties.frame_rate * UPDATE_TIME;

			if (p.image_frame >= e.properties.end_frame + 1) {
				if (e.properties.flags & tfxEmitterPropertyFlags_play_once)
					p.image_frame = e.properties.end_frame;
				else
					p.image_frame -= e.properties.end_frame + 1;
			}
			else if (p.image_frame < 0) {
				if (e.properties.flags & tfxEmitterPropertyFlags_play_once)
					p.image_frame = 0;
				else
					p.image_frame += e.properties.end_frame;
			}
		}
		return true;
	}

	/*
	float direction = 0;
	float mr_angle = 0;
	float mr_speed = 0;
	tfxVec2 mr_vec;

		float o_mr = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].direction_turbulance, p.age, p.max_age);
		float v_mr = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].velocity_turbulance, p.age, p.max_age);

			float noise_resolution = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].noise_resolution, p.age, p.max_age) * p.noise_resolution + 0.01f;
			float noise = SimplexNoise::noise(p.local.position.x / noise_resolution + p.noise_offset, p.local.position.y / noise_resolution + p.noise_offset);

	tfxVec2 velocity_normal;
	float velocity_scale = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].velocity, p.age, p.max_age) * e.current.velocity_adjuster;
	tfxVec2 current_velocity = ((p.base.velocity * velocity_scale)) * velocity_normal * UPDATE_TIME;

	tfxVec2 scale;

	float stretch = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].stretch, p.age, p.max_age);
	float velocity = std::fabsf(velocity_scale * p.base.velocity + mr_speed + p.weight_acceleration);

	float spin = 0;

	tfxVec2 offset = velocity_normal * e.current.emitter_size.y;
	float length = std::fabsf(p.local.position.y - e.current.emitter_handle.y);
	float emitter_length = e.current.emitter_size.y;
	bool is_line = e.properties.emission_type == tfxLine && e.properties.flags & tfxEmitterPropertyFlags_edge_traversal;
	bool line_and_loop = is_line && e.properties.end_behaviour == tfxLoop && length > emitter_length;
	bool line_and_kill = is_line && e.properties.end_behaviour == tfxKill && length > emitter_length;

*/

	void tfxEffect::ControlParticles(unsigned int layer) {
		bool can_bump = true;

		unsigned int next_index = 0;
		unsigned int next_buffer = !current_buffer;

		for(unsigned int i = 0 ; i != particles[layer][current_buffer].current_size; ++i) {
			auto &p = particles[layer][current_buffer][i];
			tfxEmitter &e = *p.emitter;

			if (!(p.flags & tfxParticleFlags_fresh)) {
				p.captured = p.world;
			}

			p.age += FRAME_LENGTH;

			if (common.state_flags & tfxEmitterStateFlags_remove) {
				continue;
			}
			if (p.age >= p.max_age && common.state_flags & tfxEmitterStateFlags_is_single) {
				p.age = 0;
			}
			else if (p.age >= p.max_age) {
				continue;
			}

			tfxEffect &root_effect = *e.common.root_effect;
			unsigned int flags = e.common.state_flags;
			float velocity_adjuster = e.current.velocity_adjuster;
			float global_opacity = root_effect.current.opacity;
			float image_size_y = e.image_size.y;
			float frame_rate = e.animation_frames > 1 && e.common.property_flags & tfxEmitterPropertyFlags_animate ? e.frame_rate : 0.f;
			float stretch = e.current.stretch;
			float emitter_size_y = e.current.emitter_size.y;
			float emitter_handle_y = e.common.handle.y;
			float overal_scale = e.current.overal_scale;
			float angle_offset = e.angle_offset;
			OvertimeAttributes &graphs = e.common.library->overtime_graphs[e.overtime];

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
			float lookup_opacity = graphs.opacity.lookup.values[std::min<u32>(lookup_frame, graphs.opacity.lookup.last_frame)];
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
			p.color.a = unsigned char(255.f * lookup_opacity * global_opacity);
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
			p.image_frame += frame_rate * UPDATE_TIME;
			p.image_frame = (flags & tfxEmitterStateFlags_play_once) && p.image_frame > e.end_frame ? p.image_frame = e.end_frame : p.image_frame;
			p.image_frame = (flags & tfxEmitterStateFlags_play_once) && p.image_frame < 0 ? p.image_frame = 0 : p.image_frame;
			p.image_frame = std::fmodf(p.image_frame, e.end_frame + 1);

			if (!(p.flags & tfxParticleFlags_fresh)) {
				TransformParticle(p, e);
				if (p.flags & tfxParticleFlags_capture_after_transform) {
					p.captured.position = p.world.position;
					p.flags &= ~tfxParticleFlags_capture_after_transform;
				}

				root_effect.sprites[layer][next_buffer].current_size++;
				ParticleSprite &s = root_effect.sprites[layer][next_buffer][next_index];
				s.intensity = p.intensity;
				s.color = p.color;
				s.world = p.world;
				s.captured = p.captured;
				s.ptr = e.image_ptr;
				s.handle = e.image_handle;
				s.parameters = (e.blend_mode << 28) + (unsigned int)p.image_frame;
			}
			else {
				p.flags &= ~tfxParticleFlags_fresh;
			}

			root_effect.particles[layer][next_buffer].current_size++;
			root_effect.particles[layer][next_buffer][next_index++] = p;
		}
	}

	void ControlParticles(EffectEmitter &e) {
		bool can_bump = true;

		unsigned int flags = e.flags;
		unsigned int layer = e.properties.layer;
		unsigned int blend_mode = e.properties.blend_mode;
		float velocity_adjuster = e.current.velocity_adjuster;
		float global_opacity = e.parent->current.color.a;
		float image_size_y = e.properties.image->image_size.y;
		float frame_rate = e.properties.image->animation_frames > 1 && e.properties.flags & tfxEmitterPropertyFlags_animate ? e.properties.frame_rate : 0.f;
		float end_frame = e.properties.end_frame;
		float stretch = e.current.stretch;
		float emitter_size_y = e.current.emitter_size.y;
		float emitter_handle_y = e.current.emitter_handle.y;
		float angle_offset = e.properties.angle_offset;
		float overal_scale = e.current.overal_scale;
		void *image_ptr = e.properties.image->ptr;
		tfxVec2 image_handle = e.current.image_handle;
		OvertimeAttributes &graphs = e.library->overtime_graphs[e.overtime];
		EffectEmitter &root_effect = *e.GetRootEffect();

		e.particles.reset();
		unsigned int bump_amount = 0;
		while (!e.particles.eob()) {
			auto &p = e.particles.next();

			if (!(p.flags & tfxParticleFlags_fresh)) {
				p.captured = p.world;
			}

			p.age += FRAME_LENGTH;

			//-------------------------------------------------------
			//Controll what the particle does over the course of
			//it's lifetime
			//-------------------------------------------------------

			//Before we do anything, see if the particle should be removed/end of life
			if (flags & tfxEmitterStateFlags_remove) {
				if(can_bump) bump_amount++;

				ParticleSprite &s = root_effect.sprites[layer].AtAbs(p.sprite_index);
				s.world.scale = s.captured.scale = tfxVec2();
				s.parameters = tfxINVALID;

				continue;
			}

			if (p.age >= p.max_age && flags & tfxEmitterStateFlags_is_single)
				p.age = 0.f;
			else if (p.age >= p.max_age) {
				if (can_bump) {
					bump_amount++;

					ParticleSprite &s = root_effect.sprites[layer].AtAbs(p.sprite_index);
					s.world.scale = s.captured.scale = tfxVec2();
					s.parameters = tfxINVALID;
				}

				continue;
			}
			else {
				can_bump = false;
			}

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
			float lookup_opacity = graphs.opacity.lookup.values[std::min<u32>(lookup_frame, graphs.opacity.lookup.last_frame)];
			float lookup_intensity = graphs.intensity.lookup.values[std::min<u32>(lookup_frame, graphs.intensity.lookup.last_frame)];

			float direction = 0;
			float mr_angle = 0;
			float mr_speed = 0;

			tfxVec2 mr_vec;
			if (flags & tfxEmitterStateFlags_not_line) {
				direction = lookup_direction;
			}

			if (lookup_velocity_turbulance + lookup_velocity_turbulance) {
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
			p.color.a = unsigned char(255.f * lookup_opacity * global_opacity);
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
				if(can_bump) bump_amount++;

				ParticleSprite &s = root_effect.sprites[layer].AtAbs(p.sprite_index);
				s.world.scale = s.captured.scale = tfxVec2();
				s.parameters = tfxINVALID;
				
				continue;
			}

			//----Image animation
			p.image_frame += frame_rate * UPDATE_TIME;
			p.image_frame = (flags & tfxEmitterStateFlags_play_once) && p.image_frame > end_frame ? p.image_frame = end_frame : p.image_frame;
			p.image_frame = (flags & tfxEmitterStateFlags_play_once) && p.image_frame < 0 ? p.image_frame = 0 : p.image_frame;
			p.image_frame = std::fmodf(p.image_frame, end_frame + 1);

			if (!(p.flags & tfxParticleFlags_fresh)) {
				TransformParticle(p, e);
				if (p.flags & tfxParticleFlags_capture_after_transform) {
					p.captured.position = p.world.position;
					p.flags &= ~tfxParticleFlags_capture_after_transform;
				}

				unsigned int index = p.sprite_index;
				if (index > 0) index--;
				ParticleSprite &prev = root_effect.sprites[layer].AtAbs(index);
				if (prev.parameters == tfxINVALID) {
					root_effect.sprites[layer].AtAbs(p.sprite_index).parameters = tfxINVALID;
					p.sprite_index--;
				}
				ParticleSprite &s = root_effect.sprites[layer].AtAbs(p.sprite_index);
				s.intensity = p.intensity;
				s.color = p.color;
				s.world = p.world;
				s.captured = p.captured;
				s.ptr = image_ptr;
				s.handle = image_handle;
				s.parameters = (blend_mode << 28) + (unsigned int)p.image_frame;
			}
			else {
				p.flags &= ~tfxParticleFlags_fresh;
			}
		}
		if (bump_amount) e.particles.bump(bump_amount);
	}

	bool ControlParticleFast(Particle &p, EffectEmitter &e) {
		if (e.pm->update_base_values)
			ReloadBaseValues(p, e);

		p.age += FRAME_LENGTH;

		//-------------------------------------------------------
		//Controll what the particle does over the course of
		//it's lifetime
		//-------------------------------------------------------

		//Before we do anything, see if the particle should be removed/end of life
		p.flags |= e.flags & tfxParticleFlags_remove;
		if (p.flags & tfxParticleFlags_remove)
			return false;

		bool is_single = e.properties.flags & tfxEmitterPropertyFlags_single && !(e.properties.flags & tfxEmitterPropertyFlags_one_shot) && !e.pm->disable_spawing;

		if (p.age >= p.max_age && is_single)
				p.age = 0.f;
		else if(p.age >= p.max_age)
			return false;

		u32 lookup_frame = static_cast<u32>((p.age / p.max_age * e.library->overtime_graphs[e.overtime].velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);
		OvertimeAttributes &graphs = e.library->overtime_graphs[e.overtime];

		float lookup_velocity = graphs.velocity.lookup.values[std::min<u32>(lookup_frame, graphs.velocity.lookup.last_frame)] * e.current.velocity_adjuster;
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
		float lookup_opacity = graphs.opacity.lookup.values[std::min<u32>(lookup_frame, graphs.opacity.lookup.last_frame)];
		float lookup_intensity = graphs.intensity.lookup.values[std::min<u32>(lookup_frame, graphs.intensity.lookup.last_frame)];

		float direction = 0;
		float mr_angle = 0;
		float mr_speed = 0;

		tfxVec2 mr_vec;
		if ((e.properties.emission_type != tfxLine && !(e.properties.flags & tfxEmitterPropertyFlags_edge_traversal))
			|| e.properties.emission_type == tfxLine && !(e.properties.flags & tfxEmitterPropertyFlags_edge_traversal)) {

			direction = lookup_direction;
		}
			
		if (lookup_velocity_turbulance + lookup_velocity_turbulance) {
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
		p.color.a = unsigned char(255.f * lookup_opacity * e.parent->current.color.a);
		p.intensity = lookup_intensity;
		if (!(e.properties.flags & tfxEmitterPropertyFlags_random_color)) {
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
		if (e.properties.flags & tfxEmitterPropertyFlags_lifetime_uniform_size) {
			scale.y = (lookup_width * (p.base.size.y + (velocity * lookup_stretch * e.current.stretch))) / e.properties.image->image_size.y;
			if (scale.y < scale.x)
				scale.y = scale.x;
		}
		else
			scale.y = (lookup_height * (p.base.size.y + (velocity * lookup_stretch * e.current.stretch))) / e.properties.image->image_size.y;

		//----Spin and angle Changes
		float spin = 0;
		if (e.properties.angle_setting != AngleSetting::tfxAlign && !(e.properties.flags & tfxEmitterPropertyFlags_relative_angle)) {
			spin = lookup_spin;
		}

		//---------------
		//Now that the latest changes are applied, affect the particle state
		//---------------

		//Before applying the behaviour to the particle position, scale and rotation, you have the chance to override them here
		if (e.particle_update_callback)
			e.particle_update_callback(p);

		//----Rotation
		if (e.properties.angle_setting == AngleSetting::tfxAlign) {
			tfxVec2 vd = current_velocity.IsNill() ? velocity_normal : current_velocity;
			p.local.rotation = GetVectorAngle(vd.x, vd.y) + e.properties.angle_offset;
		}
		else {
			p.local.rotation += spin * UPDATE_TIME;
		}

		//----Position
		p.local.position += current_velocity * e.current.overal_scale;

		//----Scale
		p.world.scale = scale;

		//Lines - Reposition if the particle is travelling along a line
		tfxVec2 offset = velocity_normal * e.current.emitter_size.y;
		float length = std::fabsf(p.local.position.y - e.current.emitter_handle.y);
		float emitter_length = e.current.emitter_size.y;
		bool is_line = e.properties.emission_type == tfxLine && e.properties.flags & tfxEmitterPropertyFlags_edge_traversal;
		bool line_and_loop = is_line && e.properties.end_behaviour == tfxLoop && length > emitter_length;
		bool line_and_kill = is_line && e.properties.end_behaviour == tfxKill && length > emitter_length;
		if (line_and_loop) {
			p.local.position.y -= offset.y;
			p.flags |= tfxParticleFlags_capture_after_transform;
		}
		else if (line_and_kill) {
			return false;
		}

		//----Image animation
		float frame_rate = e.properties.image->animation_frames > 1 && e.properties.flags & tfxEmitterPropertyFlags_animate ? e.properties.frame_rate : 0.f;
		p.image_frame += frame_rate * UPDATE_TIME;
		p.image_frame = (e.properties.flags & tfxEmitterPropertyFlags_play_once) && p.image_frame > e.properties.end_frame ? p.image_frame = e.properties.end_frame : p.image_frame;
		p.image_frame = (e.properties.flags & tfxParticleControlFlags_play_once) && p.image_frame < 0 ? p.image_frame = 0 : p.image_frame;
		p.image_frame = std::fmodf(p.image_frame, e.properties.end_frame + 1);

		return true;
	}

	EffectEmitter CreateEffector(float x, float y) {
		EffectEmitter effector;
		effector.transform.local.position = tfxVec2(x, y);

		return effector;
	}

	void EffectEmitter::TransformEffector(EffectEmitter &parent, bool relative_position, bool relative_angle) {
		float s = sin(transform.local.rotation);
		float c = cos(transform.local.rotation);
		transform.matrix.Set(c, s, -s, c);

		if (relative_position) {
			transform.world.rotation = parent.transform.world.rotation + transform.local.rotation;

			transform.matrix = transform.matrix.Transform(parent.transform.matrix);
			tfxVec2 rotatevec = parent.transform.matrix.TransformVector(tfxVec2(transform.local.position.x, transform.local.position.y));

			transform.world.position = parent.transform.world.position + rotatevec * parent.transform.world.scale;

		}
		else {
			transform.world.position = transform.local.position;
			transform.world.rotation = transform.local.rotation;
		}
	}

	void EffectEmitter::TransformEffector(Particle &parent, bool relative_position, bool relative_angle) {
		float s = sin(transform.local.rotation);
		float c = cos(transform.local.rotation);

		transform.matrix.Set(c, s, -s, c);

		if (relative_position) {
			transform.world.rotation = parent.world.rotation + transform.local.rotation;

			transform.matrix = transform.matrix.Transform(parent.matrix);
			tfxVec2 rotatevec = parent.matrix.TransformVector(tfxVec2(transform.local.position.x, transform.local.position.y));

			transform.world.position = parent.world.position + rotatevec * parent.world.scale;

		}
		else {
			transform.world.position = transform.local.position;
			transform.world.rotation = transform.local.rotation;
		}

	}

	void Transform(tfxEffect &e, Particle &parent) {
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
		bool line = (e.common.property_flags & tfxEmitterPropertyFlags_edge_traversal && e.emission_type == tfxLine);

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
		library->global_graphs[global].opacity.Reset(1.f, tfxGlobalOpacityPreset, add_node); library->global_graphs[global].opacity.type = tfxGlobal_opacity;
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
		library->overtime_graphs[overtime].opacity.Reset(1.f, tfxOpacityOvertimePreset, add_node); library->overtime_graphs[overtime].opacity.type = tfxOvertime_opacity;
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
			if (library->global_graphs[global].opacity.nodes.size() == 0) library->global_graphs[global].opacity.Reset(1.f, tfxOpacityOvertimePreset);
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
			if (library->overtime_graphs[overtime].opacity.nodes.size() == 0) library->overtime_graphs[overtime].opacity.Reset(1.f, tfxOpacityOvertimePreset);
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

	 void EffectEmitter::SetUpdateCallback(void(*callback)(EffectEmitter &effectemitter)) {
		update_callback = callback;
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

	void EffectEmitter::Clone(EffectEmitter &clone, EffectEmitter *root_parent, EffectLibrary *destination_library, bool keep_user_data) {
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
				clone.global = root_parent->global;
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

	bool GetEffect(EffectLibrary &library, tfxEffectStorage &storage, const char *name, tfxEffect &out) {
		assert(library.effect_paths.ValidName(name));
		EffectEmitter *effect = library.GetEffect(name);
		if (!Copy(storage, *effect, out))
			return false;
		return true;
	}

	bool GetEffect(EffectLibrary &library, tfxEffectStorage &storage, unsigned int index, tfxEffect &out) {
		assert(library.effects.size() > index);
		//if overwriting an effect that is already in storage, we need to add it to the cache
		if (out.sub_emitters.size() && out.path_hash) {
			storage.PutEffectInCache(out);
		}
		if (!Copy(storage, library.effects[index], out))
			return false;
		return true;
	}

	void tfxEffect::Reset() {
		common.age = 0.f;
		common.highest_particle_age = 0.f;
		common.active_children = sub_emitters.size();
		common.timeout_counter = 0;
		ClearParticles();
		sub_emitters.reset();
		while (!sub_emitters.eob()) {
			tfxEmitter &e = sub_emitters.next();
			e.Reset();
		}
	}

	void tfxEffect::ClearParticles() {
		for (int i = 0; i != tfxLAYERS; ++i) {
			particles[i][0].clear();
			particles[i][1].clear();
			sprites[i][0].clear();
			sprites[i][1].clear();
		}
	}

	void tfxEmitter::Reset() {
		common.age = 0.f;
		common.highest_particle_age = 0.f;
		common.active_children = sub_effects.size();
		common.timeout_counter = 0;
		common.state_flags &= ~tfxEmitterStateFlags_single_shot_done;
		common.state_flags &= ~tfxEmitterStateFlags_stop_spawning;
		current.emission_alternator = 0.f;
		current.amount_remainder = 0.f;
		current.grid_coords = tfxVec2();
		current.grid_direction = tfxVec2();
		sub_effects.reset();
		while (!sub_effects.eob()) {
			tfxEffect &e = sub_effects.next();
			e.Reset();
		}
	}

	bool Copy(tfxEffectStorage &storage, EffectEmitter &in, tfxEffect &out) {
		int emitters = 0;
		int effects = 0;
		if (!storage.GetEffectFromCache(in.path_hash, out)) {
			in.CountChildren(emitters, effects);
			unsigned int emitter_mem_req = emitters * sizeof(tfxEmitter);
			unsigned int effect_mem_req = effects * sizeof(tfxEffect);
			if (storage.emitter_memory.free_unused_space() < emitter_mem_req)
				return false;
			if (storage.effect_memory.free_unused_space() < effect_mem_req) 
				return false;
			in.CopyToEffect(out, storage);
		}
		else {
			out.ClearParticles();
			out.Reset();
		}
		return true;
	}

	void EffectEmitter::CopyToEffect(tfxEffect &e, tfxEffectStorage &storage) {
		assert(type == tfxEffectType);		//Must be called on an effect type only
		tfxEffectID effect_id = storage.InsertEffect();
		e = storage.GetEffect(effect_id);
		e.common.property_flags = properties.flags;
		e.common.state_flags = flags;
		e.common.library = library;
		e.path_hash = path_hash;
		e.global = global;
		e.common.loop_length = properties.loop_length;
		e.common.handle = properties.emitter_handle;
		e.common.lookup_mode = tfxFast;
		for (int i = 0; i != tfxLAYERS; ++i) {
			e.common.max_particles[i] = max_particles[i];
			e.sprites[i][0].assign_memory(storage.sprite_memory, sizeof(ParticleSprite), e.common.max_particles[i]);
			e.sprites[i][1].assign_memory(storage.sprite_memory, sizeof(ParticleSprite), e.common.max_particles[i]);
			e.particles[i][0].assign_memory(storage.particle_memory, sizeof(Particle), e.common.max_particles[i]);
			e.particles[i][1].assign_memory(storage.particle_memory, sizeof(Particle), e.common.max_particles[i]);
		}
		e.sub_emitters.assign_memory(storage.effect_memory, sizeof(tfxEmitter), sub_effectors.size());
		for (auto &emitter : sub_effectors) {
			tfxEmitter new_emitter;
			new_emitter.common.root_effect = &e;
			new_emitter.parent_effect = &e;
			emitter.CopyToEmitter(new_emitter, storage);
			assert(e.sub_emitters.push_back(new_emitter));
		}
	}

	void EffectEmitter::CopyToEmitter(tfxEmitter &e, tfxEffectStorage &storage) {
		//Internal use only, CopyToEffect must be called first
		e.common.property_flags = properties.flags;
		e.common.state_flags = flags;
		e.common.loop_length = properties.loop_length;	
		e.common.library = library;
		e.base = base;
		e.property = property;
		e.variation = variation;
		e.overtime = overtime;
		e.layer = properties.layer;
		e.blend_mode = properties.blend_mode;
		e.angle_setting = properties.angle_setting;
		e.emission_direction = properties.emission_direction;
		e.emission_type = properties.emission_type;
		e.spawn_amount = properties.spawn_amount;
		e.image_size = properties.image->image_size;
		e.end_frame = properties.end_frame;
		e.animation_frames = properties.image->animation_frames;
		e.start_frame = properties.start_frame;
		e.frame_rate = properties.frame_rate;
		e.image_ptr = properties.image->ptr;
		e.angle_offset = properties.angle_offset;
		e.image_handle = current.image_handle;
		e.common.handle = properties.emitter_handle;
		
		e.common.lookup_mode = tfxFast;
		if (properties.flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
			e.image_handle = tfxVec2(0.5f, 0.5f);
		}
		else {
			e.image_handle = properties.image_handle;
		}
		e.grid_points = properties.grid_points;
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

		e.common.max_particles[e.layer] = max_particles[e.layer];
		e.sub_effects.assign_memory(storage.effect_memory, sizeof(tfxEffect), sub_effectors.size());
		for (auto &effect : sub_effectors) {
			tfxEffect new_effect;
			new_effect.common.root_effect = e.common.root_effect;
			effect.CopyToSubEffect(new_effect, storage);
			assert(e.sub_effects.push_back(new_effect));
		}
	}

	void EffectEmitter::CopyToSubEffect(tfxEffect &e, tfxEffectStorage &storage) {
		e.common.property_flags = properties.flags;
		e.common.state_flags = flags;
		e.common.library = library;
		e.global = global;
		e.common.loop_length = properties.loop_length;
		e.sub_emitters.assign_memory(storage.effect_memory, sizeof(tfxEmitter), sub_effectors.size());
		for (auto &emitter : sub_effectors) {
			tfxEmitter new_emitter;
			new_emitter.common.root_effect = e.common.root_effect;
			emitter.CopyToEmitter(new_emitter, storage);
			assert(e.sub_emitters.push_back(new_emitter));
		}
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
			library->global_graphs[global].opacity.Free();
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
			library->overtime_graphs[overtime].opacity.Free();
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
		assert(effect_paths.ValidName(path));
		return effect_paths.At(path);
	}

	EffectEmitter* EffectLibrary::GetEffect(const char *path) {
		assert(effect_paths.ValidName(path));
		return effect_paths.At(path);
	}

	EffectEmitter* EffectLibrary::GetEffect(tfxKey key) {
		assert(effect_paths.ValidKey(key));
		return effect_paths.At(key);
	}

	void EffectLibrary::PrepareEffectTemplate(tfxText path, EffectEmitterTemplate &effect_template) {
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
		for (auto &e : effects) {
			if(e.type != tfxFolder)
				e.UpdateAllBufferSizes(e.max_particles);
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
			CompileGraph(g.opacity);
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
			CompileGraphOvertime(g.opacity);
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
		CompileGraph(g.opacity);
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
		CompileGraphOvertime(g.opacity);
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
		graph_min_max[tfxGlobal_opacity] = GetMinMaxGraphValues(tfxOpacityOvertimePreset);
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
		graph_min_max[tfxOvertime_opacity] = GetMinMaxGraphValues(tfxOpacityOvertimePreset);
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
			if (values[0] == "global_opacity") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].opacity.AddNode(n); }
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
			if (values[0] == "overtime_opacity") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].opacity.AddNode(n); }
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
	}
	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, int value) {
		if (field == "emission_type")
			effect.properties.emission_type = (EmissionType)value;
		if (field == "emission_direction")
			effect.properties.emission_direction = (EmissionDirection)value;
		if (field == "blend_mode")
			effect.properties.blend_mode = (BlendMode)value;
		if (field == "angle_setting")
			effect.properties.angle_setting = (AngleSetting)value;
		if (field == "color_option")
			effect.library->animation_settings[effect.animation_settings].color_option = (ExportColorOptions)value;
		if (field == "export_option")
			effect.library->animation_settings[effect.animation_settings].export_option = (ExportOptions)value;
		if (field == "end_behaviour")
			effect.properties.end_behaviour = (LineTraversalEndBehaviour)value;
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
		file.AddLine("blend_mode=%i", property.blend_mode);
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

	ParticleManager::~ParticleManager() {
	}

	EffectEmitter &ParticleManager::operator[] (unsigned int index) {
		return effects[current_ebuff][index];
	}
	
	void ParticleManager::AddEffect(EffectEmitter &effect, unsigned int buffer) {
		if (effects[buffer].current_size == effects[buffer].capacity)
			return;
		if (use_compute_shader && highest_compute_controller_index >= max_compute_controllers && free_compute_controllers.empty())
			return;
		unsigned int parent_index = effects[buffer].current_size++;
		effects[buffer][parent_index] = effect;
		effects[buffer][parent_index].active_children = 0;
		effects[buffer][parent_index].flags &= ~tfxEmitterStateFlags_retain_matrix;
		effects[buffer][parent_index].pm = this;
		effects[buffer][parent_index].ResetParents();
		for (int i = 0; i != tfxLAYERS; ++i) {
			EffectEmitter &effect = effects[buffer][parent_index];
			effects[buffer][parent_index].sprites[i].assign_memory(sprite_memory, sizeof(ParticleSprite), effects[buffer][parent_index].max_particles[i]);
		}
		for (auto &e : effect.sub_effectors) {
			if (!FreeEffectCapacity())
				break;
			if (e.flags & tfxEmitterStateFlags_enabled) {
				unsigned int index = effects[buffer].current_size++;
				effects[buffer][index] = e;
				EffectEmitter &emitter = effects[buffer].back();
				emitter.parent = &effects[buffer][parent_index];
				emitter.next_ptr = emitter.parent;
				emitter.pm = this;
				emitter.active_children = 0;
				emitter.flags &= ~tfxEmitterStateFlags_retain_matrix;

				emitter.flags |= e.properties.flags & tfxEmitterPropertyFlags_single && !(e.properties.flags & tfxEmitterPropertyFlags_one_shot) && !disable_spawing ? tfxEmitterStateFlags_is_single : 0;
				emitter.flags |= (e.properties.emission_type != tfxLine && !(e.properties.flags & tfxEmitterPropertyFlags_edge_traversal)) || e.properties.emission_type == tfxLine && !(e.properties.flags & tfxEmitterPropertyFlags_edge_traversal) ? tfxEmitterStateFlags_not_line : 0;
				emitter.flags |= e.properties.flags & tfxEmitterPropertyFlags_random_color;
				emitter.flags |= e.properties.flags & tfxEmitterPropertyFlags_lifetime_uniform_size;
				emitter.flags |= e.properties.angle_setting != AngleSetting::tfxAlign && !(e.properties.flags & tfxEmitterPropertyFlags_relative_angle) ? tfxEmitterStateFlags_can_spin : 0;
				emitter.flags |= e.properties.angle_setting == AngleSetting::tfxAlign ? tfxEmitterStateFlags_align_with_velocity : 0;
				emitter.flags |= e.properties.emission_type == tfxLine && e.properties.flags & tfxEmitterPropertyFlags_edge_traversal ? tfxEmitterStateFlags_is_line_traversal : 0;
				emitter.flags |= e.properties.flags & tfxEmitterPropertyFlags_play_once;
				emitter.flags |= e.properties.end_behaviour == tfxLoop ? tfxEmitterStateFlags_loop : 0;
				emitter.flags |= e.properties.end_behaviour == tfxKill ? tfxEmitterStateFlags_kill : 0;

				effects[buffer][parent_index].active_children++;
				if (contain_particles_in_emitters) {
					emitter.particles.assign_memory(particle_memory, sizeof(Particle), emitter.max_particles[emitter.properties.layer]);
				}
				if (use_compute_shader && e.sub_effectors.empty()) {
					int free_slot = AddComputeController();
					if (free_slot != -1) {
						emitter.compute_slot_id = free_slot;
						emitter.properties.flags |= tfxEmitterPropertyFlags_is_bottom_emitter;
						ComputeController &c = *(static_cast<ComputeController*>(compute_controller_ptr) + free_slot);
						c.flags = 0;
						c.flags |= emitter.properties.flags & tfxParticleControlFlags_random_color;
						c.flags |= emitter.properties.flags & tfxParticleControlFlags_relative_position;
						c.flags |= emitter.properties.flags & tfxParticleControlFlags_relative_angle;
						c.flags |= emitter.properties.flags & tfxParticleControlFlags_edge_traversal;
						c.flags |= emitter.properties.flags & tfxParticleControlFlags_base_uniform_size;
						c.flags |= emitter.properties.flags & tfxParticleControlFlags_lifetime_uniform_size;
						c.flags |= emitter.properties.flags & tfxParticleControlFlags_reverse_animation;
						c.flags |= emitter.properties.flags & tfxParticleControlFlags_play_once;
						c.flags |= ((emitter.properties.emission_type == tfxPoint) << 3);
						c.flags |= ((emitter.properties.emission_type == tfxArea) << 4);
						c.flags |= ((emitter.properties.emission_type == tfxLine) << 5);
						c.flags |= ((emitter.properties.emission_type == tfxEllipse) << 6);
						c.flags |= ((emitter.properties.end_behaviour == tfxLoop) << 7);
						c.flags |= ((emitter.properties.end_behaviour == tfxKill) << 8);
						c.flags |= ((emitter.properties.end_behaviour == tfxLetFree) << 9);
						c.flags |= ((emitter.properties.angle_setting == tfxAlign) << 17);
						c.flags |= ((emitter.properties.angle_setting == tfxAlignWithEmission) << 18);
						c.flags |= ((emitter.properties.angle_setting == tfxRandom) << 19);
						c.flags |= ((emitter.properties.angle_setting == tfxSpecify) << 20);
						c.flags |= ((emitter.properties.blend_mode == tfxAlpha) << 21);
						c.flags |= ((emitter.properties.blend_mode == tfxAdditive) << 22);
						emitter.properties.compute_flags = c.flags;
						if (emitter.properties.flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
							c.image_handle = tfxVec2(0.5f, 0.5f);
						}
						else {
							c.image_handle = emitter.properties.image_handle;
						}
						c.image_handle = emitter.properties.image_handle;
					}
				}
			}
		}
		effects[buffer][parent_index].NoTweenNextUpdate();
	}

	void ParticleManager::AddEffect(EffectEmitterTemplate &effect, unsigned int buffer) {
		AddEffect(effect.effect_template, current_ebuff);
	}

	int ParticleManager::AddComputeController() {
		//Compute slots should only ever be added for the bottom emitter that has no sub effects
		unsigned int free_slot;
		if (!free_compute_controllers.empty()) {
			free_slot = free_compute_controllers.pop_back();
		}
		else {
			free_slot = highest_compute_controller_index++;
		}
		if (free_slot >= max_compute_controllers)
			return -1;
		return free_slot;
	}

	void ParticleManager::ResetParticlePtr(void *ptr) {
		new_compute_particle_ptr = ptr;
		new_compute_particle_index = 0;
	}

	void ParticleManager::ResetControllerPtr(void *ptr) {
		compute_controller_ptr = ptr;
	}

	void ParticleManager::UpdateCompute(void *sampled_particles, unsigned int sample_size) {
		for (int i = 0; i != sample_size; ++i) {
			if (compute_global_state.current_length == 0)
				break;
			ComputeParticle *sample = static_cast<ComputeParticle*>(sampled_particles) + i;
			if (sample->age > sample->max_age) {
				compute_global_state.start_index++;
				compute_global_state.start_index %= compute_global_state.end_index;
				compute_global_state.current_length--;
				sample->age = 0;
			}
			else {
				break;
			}
		}
	}

	uint32_t ParticleManager::AddParticle(unsigned int layer, Particle &p) {
		assert(particles[layer][current_pbuff].current_size != particles[layer][current_pbuff].capacity);
		particles[layer][current_pbuff][particles[layer][current_pbuff].current_size] = p;
		particles[layer][current_pbuff].current_size++;
		return (uint32_t)particles[layer][current_pbuff].current_size - 1;
	}

	Particle& ParticleManager::GrabCPUParticle(unsigned int layer) {
		assert(particles[layer][current_pbuff].current_size != particles[layer][current_pbuff].capacity);
		particles[layer][current_pbuff].current_size++;
		//particles[layer][current_pbuff][particles[layer][current_pbuff].current_size -1].particle_id = particle_id++;
		return particles[layer][current_pbuff][particles[layer][current_pbuff].current_size - 1];
	}

	ComputeParticle& ParticleManager::GrabComputeParticle(unsigned int layer) {
		assert(new_compute_particle_ptr);		//Use must assign the compute ptr to point to an area in memory where you can stage new particles for uploading to the GPU - See ResetComputePtr
		return *(static_cast<ComputeParticle*>(new_compute_particle_ptr) + new_compute_particle_index++);
	}

	void ParticleManager::Update() {
		new_compute_particle_index = 0;
		unsigned int index = 0;

		unsigned int next_buffer = !current_ebuff;
		effects[next_buffer].clear();

		for (auto &e : effects[current_ebuff]) {
			e.Update();
			if (e.type == tfxEffectType) {
				if (e.active_children > 0) {
					e.next_ptr = SetNextEffect(e, next_buffer);
				}
				else {
					e.next_ptr = nullptr;
					e.sub_effectors.free_all();
				}
			}
			else {
				if (e.timeout_counter < e.timeout) {
					e.next_ptr = SetNextEffect(e, next_buffer);
				}
				else {
					e.next_ptr = nullptr;
					e.sub_effectors.free_all();
					e.particles.free_range(particle_memory);
					if (use_compute_shader && e.properties.flags & tfxEmitterPropertyFlags_is_bottom_emitter)
						FreeComputeSlot(e.compute_slot_id);
				}
			}
			index++;
		}

		current_ebuff = next_buffer;

		next_buffer = !current_pbuff;

		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			particles[layer][next_buffer].clear();

			index = 0;
			for (auto &p : particles[layer][current_pbuff]) {
				p.parent = p.parent->next_ptr;

				if (!(p.flags & tfxParticleFlags_fresh)) {

					p.captured = p.world;

					if (ControlParticleFast(p, *p.parent)) {
						TransformParticle(p, *p.parent);
						if (p.flags & tfxParticleFlags_capture_after_transform) {
							p.captured.position = p.world.position;
							p.flags &= ~tfxParticleFlags_capture_after_transform;
						}

						p.next_ptr = SetNextParticle(layer, p, next_buffer);
					}
					else {
						p.next_ptr = nullptr;
					}
				}
				else {
					p.flags &= ~tfxParticleFlags_fresh;
					p.next_ptr = SetNextParticle(layer, p, next_buffer);
				}
				index++;
			}
		}

		current_pbuff = next_buffer;

		update_base_values = false;

	}

	inline Particle* ParticleManager::SetNextParticle(unsigned int layer, Particle &p, unsigned int buffer) {
		unsigned int index = particles[layer][buffer].current_size++;
		assert(index < particles[layer][buffer].capacity);
		particles[layer][buffer][index] = p;
		return &particles[layer][buffer][index];
	}

	inline EffectEmitter* ParticleManager::SetNextEffect(EffectEmitter &e, unsigned int buffer) {
		unsigned int index = effects[buffer].current_size++;
		assert(index < effects[buffer].capacity);
		effects[buffer][index] = e;
		return &effects[buffer][index];
	}

	void ParticleManager::Render(float tween, void *data) {
		if (!render_func)
			return;
		for (auto &p : particles[current_pbuff]) {
			render_func(tween, &p, data);
		}
	}
	
	tfxvec<Particle> *ParticleManager::GetParticleBuffer(unsigned int layer) {
		return &particles[layer][current_pbuff];
	}
	
	tfxvec<EffectEmitter> *ParticleManager::GetEffectBuffer() {
		return &effects[current_ebuff];
	}

	void ParticleManager::SetRenderCallback(void func(float, void*, void*)) {
		render_func = func;
	}
	
	void tfxEffectStorage::Init(unsigned int max_effects, unsigned int max_emitters, unsigned int max_particles) {
		effect_memory.reserve(sizeof(tfxEffect) * max_effects);
		emitter_memory.reserve(sizeof(tfxEmitter) * max_emitters);
		particle_memory.reserve(sizeof(Particle) * max_particles);
		sprite_memory.reserve(sizeof(ParticleSprite) * max_particles);
	}

	void tfxEffectStorage::FreeEffect(tfxEffect &effect) {
		if (effect_cache.ValidKey(effect.path_hash)) {
			tfxvec<tfxEffect> &bucket = effect_cache.At(effect.path_hash);
			bucket.push_back(effect);
		}
		else {
			tfxvec<tfxEffect> new_bucket;
			effect_cache.Insert(effect.path_hash, new_bucket);
			tfxvec<tfxEffect> &bucket = effect_cache.At(effect.path_hash);
			bucket.push_back(effect);
		}
	}

	bool tfxEffectStorage::EffectInCache(tfxKey path_hash) {
		if (effect_cache.ValidKey(path_hash)) {
			tfxvec<tfxEffect> &bucket = effect_cache.At(path_hash);
			return bucket.size() > 0;
		}
		return false;
	}

	bool tfxEffectStorage::GetEffectFromCache(tfxKey path_hash, tfxEffect &effect) {
		if (effect_cache.ValidKey(path_hash)) {
			tfxvec<tfxEffect> &bucket = effect_cache.At(path_hash);
			if (bucket.size() > 0) {
				effect = bucket.pop_back();
				return true;
			}
		}
		return false;
	}

	void tfxEffectStorage::PutEffectInCache(tfxEffect &effect) {
		if (effect_cache.ValidKey(effect.path_hash)) {
			tfxvec<tfxEffect> &bucket = effect_cache.At(effect.path_hash);
			bucket.push_back(effect);
		}
		else {
			tfxvec<tfxEffect> new_bucket;
			new_bucket.push_back(effect);
			effect_cache.Insert(effect.path_hash, new_bucket);
		}
	}

	tfxEffectID tfxEffectStorage::InsertEffect() {
		tfxEffectID id = effect_memory.get_range(sizeof(tfxEffect));
		void *ptr = (char*)effect_memory.data + effect_memory.ranges[id].offset_into_memory;
		new((void*)(ptr)) tfxEffect();
		return id;
	}

	tfxEffect &tfxEffectStorage::GetEffect(tfxEffectID id) {
		void *ptr = (char*)effect_memory.data + effect_memory.ranges[id].offset_into_memory;
		return *static_cast<tfxEffect*>(ptr);
	}

	void ParticleManager::Init(unsigned int effects_limit, unsigned int particle_limit_per_layer) {
		max_cpu_particles_per_layer = particle_limit_per_layer;
		max_effects = effects_limit;

		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			particles[layer][0].resize(max_cpu_particles_per_layer);
			particles[layer][1].resize(max_cpu_particles_per_layer);
			sprite_memory.reserve(max_cpu_particles_per_layer * sizeof(ComputeSprite));
			particles[layer][0].clear();
			particles[layer][1].clear();
		}
		particle_memory.reserve(max_cpu_particles_per_layer * sizeof(Particle));
		effects[0].create_pool(max_effects);
		effects[1].create_pool(max_effects);
		effects[0].clear();
		effects[1].clear();

	}

	uint32_t ParticleManager::ParticleCount() {
		unsigned int count = 0;
		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			count += particles[layer][current_pbuff].size();
		}
		return count;
	}

	void ParticleManager::ClearAll() {
		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			particles[layer][0].clear();
			particles[layer][1].clear();
		}
		for (unsigned int i = 0; i != 2; ++i) {
			for (auto &e : effects[i]) {
				e.sub_effectors.free_all();
			}
			effects[i].clear();
		}
		particle_id = 0;
	}
	void ParticleManager::SoftExpireAll() {
		for (auto &e : effects[current_ebuff]) {
			e.flags |= tfxEmitterStateFlags_stop_spawning;
		}
	}
	void ParticleManager::ClearDepths() {
		for (unsigned int i = 0; i != 2; ++i) {
			for (auto &e : effects[i]) {
				e.sub_effectors.free_all();
			}
			effects[i].clear();
		}
	}
	void ParticleManager::SetLookUpMode(LookupMode mode) {
		if (mode == tfxPrecise) {
			lookup_overtime_callback = LookupPreciseOvertime;
			lookup_callback = LookupPrecise;
		}
		else {
			lookup_overtime_callback = LookupFastOvertime;
			lookup_callback = LookupFast;
		}
		lookup_mode = mode;
	}

	void ParticleManager::UpdateBaseValues() {
		update_base_values = true;
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

	void StopSpawning(ParticleManager &pm) {
		pm.SoftExpireAll();
	}

	 void RemoveAllEffects(ParticleManager &pm) {
		pm.ClearAll();
	}

	 void InitParticleManager(ParticleManager &pm, unsigned int effects_limit, unsigned int particle_limit_per_layer) {
		pm.Init(effects_limit, particle_limit_per_layer);
	}

	void AddEffect(ParticleManager &pm, EffectEmitter &effect, float x, float y) {
		effect.Position(x, y);
		pm.AddEffect(effect, pm.current_ebuff);
	}

	void AddEffect(ParticleManager &pm, EffectEmitterTemplate &effect, float x, float y) {
		effect.effect_template.Position(x, y);
		pm.AddEffect(effect, pm.current_ebuff);
	}

	void EffectEmitterTemplate::SetUserDataAll(void *data) {
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

	 void EffectEmitterTemplate::SetUpdateCallbackAll(void(*update_callback)(EffectEmitter &effectemitter)) {
		tfxvec<EffectEmitter*> stack;
		stack.push_back(&effect_template);
		while (stack.size()) {
			EffectEmitter *current = stack.pop_back();
			current->update_callback = update_callback;
			for (auto &sub : current->sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	void EffectEmitterTemplate::SetParticleUpdateCallback(tfxText path, void(*particle_update_callback)(Particle &particle)) {
		assert(paths.ValidName(path));
		EffectEmitter &e = *paths.At(path);
		assert(e.type == tfxEmitterType);
		e.particle_update_callback = particle_update_callback;
	}

	void EffectEmitterTemplate::SetParticleOnSpawnCallback(tfxText path, void(*particle_onspawn_callback)(Particle &particle)) {
		assert(paths.ValidName(path));
		EffectEmitter &e = *paths.At(path);
		assert(e.type == tfxEmitterType);
		e.particle_onspawn_callback = particle_onspawn_callback;
	}

}