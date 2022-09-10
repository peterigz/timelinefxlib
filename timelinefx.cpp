#include "timelinefx.h"

namespace tfx {

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

	//simd floor function thanks to St�phanie Rancourt: http://dss.stephanierct.com/DevBlog/?p=8
	inline __m128 _mm_floor_ps2(const __m128& x) {
		__m128i v0 = _mm_setzero_si128();
		__m128i v1 = _mm_cmpeq_epi32(v0, v0);
		__m128i ji = _mm_srli_epi32(v1, 25);
		__m128 j = *(__m128*)&_mm_slli_epi32(ji, 23); //create vector 1.0f
		__m128i i = _mm_cvttps_epi32(x);
		__m128 fi = _mm_cvtepi32_ps(i);
		__m128 igx = _mm_cmpgt_ps(fi, x);
		j = _mm_and_ps(igx, j);
		return _mm_sub_ps(fi, j);
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
	const int32_t perm[] =
	{ 151,160,137,91,90,15,
	131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
	190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
	88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
	77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
	102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
	135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
	5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
	223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
	129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
	251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
	49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
	138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
	151,160,137,91,90,15,
	131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
	190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
	88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
	77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
	102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
	135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
	5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
	223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
	129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
	251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
	49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
	138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
	};

	static const int permMOD12[] =
	{
	7, 4, 5, 7, 6, 3, 11, 1, 9, 11, 0, 5, 2, 5, 7, 9, 8, 0, 7, 6, 9, 10, 8, 3,
	1, 0, 9, 10, 11, 10, 6, 4, 7, 0, 6, 3, 0, 2, 5, 2, 10, 0, 3, 11, 9, 11, 11,
	8, 9, 9, 9, 4, 9, 5, 8, 3, 6, 8, 5, 4, 3, 0, 8, 7, 2, 9, 11, 2, 7, 0, 3, 10,
	5, 2, 2, 3, 11, 3, 1, 2, 0, 7, 1, 2, 4, 9, 8, 5, 7, 10, 5, 4, 4, 6, 11, 6,
	5, 1, 3, 5, 1, 0, 8, 1, 5, 4, 0, 7, 4, 5, 6, 1, 8, 4, 3, 10, 8, 8, 3, 2, 8,
	4, 1, 6, 5, 6, 3, 4, 4, 1, 10, 10, 4, 3, 5, 10, 2, 3, 10, 6, 3, 10, 1, 8, 3,
	2, 11, 11, 11, 4, 10, 5, 2, 9, 4, 6, 7, 3, 2, 9, 11, 8, 8, 2, 8, 10, 7, 10, 5,
	9, 5, 11, 11, 7, 4, 9, 9, 10, 3, 1, 7, 2, 0, 2, 7, 5, 8, 4, 10, 5, 4, 8, 2, 6,
	1, 0, 11, 10, 2, 1, 10, 6, 0, 0, 11, 11, 6, 1, 9, 3, 1, 7, 9, 2, 11, 11, 1, 0,
	10, 7, 1, 7, 10, 1, 4, 0, 0, 8, 7, 1, 2, 9, 7, 4, 6, 2, 6, 8, 1, 9, 6, 6, 7, 5,
	0, 0, 3, 9, 8, 3, 6, 6, 11, 1, 0, 0,
	7, 4, 5, 7, 6, 3, 11, 1, 9, 11, 0, 5, 2, 5, 7, 9, 8, 0, 7, 6, 9, 10, 8, 3,
	1, 0, 9, 10, 11, 10, 6, 4, 7, 0, 6, 3, 0, 2, 5, 2, 10, 0, 3, 11, 9, 11, 11,
	8, 9, 9, 9, 4, 9, 5, 8, 3, 6, 8, 5, 4, 3, 0, 8, 7, 2, 9, 11, 2, 7, 0, 3, 10,
	5, 2, 2, 3, 11, 3, 1, 2, 0, 7, 1, 2, 4, 9, 8, 5, 7, 10, 5, 4, 4, 6, 11, 6,
	5, 1, 3, 5, 1, 0, 8, 1, 5, 4, 0, 7, 4, 5, 6, 1, 8, 4, 3, 10, 8, 8, 3, 2, 8,
	4, 1, 6, 5, 6, 3, 4, 4, 1, 10, 10, 4, 3, 5, 10, 2, 3, 10, 6, 3, 10, 1, 8, 3,
	2, 11, 11, 11, 4, 10, 5, 2, 9, 4, 6, 7, 3, 2, 9, 11, 8, 8, 2, 8, 10, 7, 10, 5,
	9, 5, 11, 11, 7, 4, 9, 9, 10, 3, 1, 7, 2, 0, 2, 7, 5, 8, 4, 10, 5, 4, 8, 2, 6,
	1, 0, 11, 10, 2, 1, 10, 6, 0, 0, 11, 11, 6, 1, 9, 3, 1, 7, 9, 2, 11, 11, 1, 0,
	10, 7, 1, 7, 10, 1, 4, 0, 0, 8, 7, 1, 2, 9, 7, 4, 6, 2, 6, 8, 1, 9, 6, 6, 7, 5,
	0, 0, 3, 9, 8, 3, 6, 6, 11, 1, 0, 0
	};

	/**
	 * Helper function to hash an integer using the above permutation table
	 *
	 *  This inline function costs around 1ns, and is called N+1 times for a noise of N dimension.
	 *
	 *  Using a real hash function would be better to improve the "repeatability of 256" of the above permutation table,
	 * but fast integer Hash functions uses more time and have bad random GetProperties().
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
		float test = ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
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
	float tfxSimplexNoise::noise(float x) {
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
	float tfxSimplexNoise::noise(float x, float y) {
		tfxPROFILE;
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
		tfxU32 ii = i & 255;
		tfxU32 jj = j & 255;
		const int gi0 = permMOD12[perm[ii + perm[jj]]];
		const int gi1 = permMOD12[perm[ii + i1 + perm[jj + j1]]];
		const int gi2 = permMOD12[perm[ii + 1 + perm[jj + 1]]];

		// Calculate the contribution from the first corner
		float t0 = 0.5f - x0 * x0 - y0 * y0;
		if (t0 < 0.0f) {
			n0 = 0.0f;
		}
		else {
			t0 *= t0;
			n0 = t0 * t0 * Dot(gradX[gi0], gradY[gi0], x0, y0);
		}

		// Calculate the contribution from the second corner
		float t1 = 0.5f - x1 * x1 - y1 * y1;
		if (t1 < 0.0f) {
			n1 = 0.0f;
		}
		else {
			t1 *= t1;
			n1 = t1 * t1 * Dot(gradX[gi1], gradY[gi1],  x1, y1);
		}

		// Calculate the contribution from the third corner
		float t2 = 0.5f - x2 * x2 - y2 * y2;
		if (t2 < 0.0f) {
			n2 = 0.0f;
		}
		else {
			t2 *= t2;
			n2 = t2 * t2 * Dot(gradX[gi2], gradY[gi2], x2, y2);
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
	float tfxSimplexNoise::noise(float x, float y, float z) {
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

	//A Simd (SSE3) version of simplex noise allowing you to do 4 samples with 1 call for a speed boost
	tfxVec4 tfxSimplexNoise::noise4(const __m128 &x4, const __m128 &y4, const __m128 &z4) {
		tfxPROFILE;
		// Skewing/Unskewing factors for 3D

		// Skew the input space to determine which simplex cell we're in
		//float s = (v1.x + v1.y + v1.z) * F3; // Very nice and simple skew factor for 3D
		__m128 s4 = _mm_mul_ps(_mm_add_ps(x4, _mm_add_ps(y4, z4)), tfxF3_4);
		__m128 x4_s4 = _mm_add_ps(x4, s4);
		__m128 y4_s4 = _mm_add_ps(y4, s4);
		__m128 z4_s4 = _mm_add_ps(z4, s4);
		__m128 i = _mm_floor_ps2(x4_s4);
		__m128 j = _mm_floor_ps2(y4_s4);
		__m128 k = _mm_floor_ps2(z4_s4);
		__m128 t = _mm_add_ps(i, j);
		t = _mm_add_ps(t, k);
		t = _mm_mul_ps(t, tfxG3_4);

		__m128 X0 = _mm_sub_ps(i, t); // Unskew the cell origin back to (v1.x,v1.y,v1.z) space
		__m128 Y0 = _mm_sub_ps(j, t);
		__m128 Z0 = _mm_sub_ps(k, t);
		__m128 x0 = _mm_sub_ps(x4, X0); // The v1.x,v1.y,v1.z distances from the cell origin
		__m128 y0 = _mm_sub_ps(y4, Y0);
		__m128 z0 = _mm_sub_ps(z4, Z0);

		// For the 3D case, the simplex shape is a slightly irregular tetrahedron.
		// Determine which simplex we are in.
		m128i i1, i2, j1, j2, k1, k2;

		i1.m = _mm_and_si128(tfxONE, _mm_and_si128(_mm_castps_si128(_mm_cmpge_ps(x0, y0)), _mm_castps_si128(_mm_cmpge_ps(x0, z0))));
		j1.m = _mm_and_si128(tfxONE, _mm_and_si128(_mm_castps_si128(_mm_cmpgt_ps(y0, x0)), _mm_castps_si128(_mm_cmpge_ps(y0, z0))));
		k1.m = _mm_and_si128(tfxONE, _mm_and_si128(_mm_castps_si128(_mm_cmpgt_ps(z0, x0)), _mm_castps_si128(_mm_cmpgt_ps(z0, y0))));

		//for i2
		__m128i yx_xz = _mm_and_si128(_mm_castps_si128(_mm_cmpge_ps(x0, y0)), _mm_castps_si128(_mm_cmplt_ps(x0, z0)));
		__m128i zx_xy = _mm_and_si128(_mm_castps_si128(_mm_cmpge_ps(x0, z0)), _mm_castps_si128(_mm_cmplt_ps(x0, y0)));

		//for j2
		__m128i xy_yz = _mm_and_si128(_mm_castps_si128(_mm_cmplt_ps(x0, y0)), _mm_castps_si128(_mm_cmplt_ps(y0, z0)));
		__m128i zy_yx = _mm_and_si128(_mm_castps_si128(_mm_cmpge_ps(y0, z0)), _mm_castps_si128(_mm_cmpge_ps(x0, y0)));

		//for k2
		__m128i yz_zx = _mm_and_si128(_mm_castps_si128(_mm_cmplt_ps(y0, z0)), _mm_castps_si128(_mm_cmpge_ps(x0, z0)));
		__m128i xz_zy = _mm_and_si128(_mm_castps_si128(_mm_cmplt_ps(x0, z0)), _mm_castps_si128(_mm_cmpge_ps(y0, z0)));

		i2.m = _mm_and_si128(tfxONE, _mm_or_si128(i1.m, _mm_or_si128(yx_xz, zx_xy)));
		j2.m = _mm_and_si128(tfxONE, _mm_or_si128(j1.m, _mm_or_si128(xy_yz, zy_yx)));
		k2.m = _mm_and_si128(tfxONE, _mm_or_si128(k1.m, _mm_or_si128(yz_zx, xz_zy)));

		__m128 x1 = _mm_add_ps(_mm_sub_ps(x0, _mm_cvtepi32_ps(i1.m)), tfxG3_4);
		__m128 y1 = _mm_add_ps(_mm_sub_ps(y0, _mm_cvtepi32_ps(j1.m)), tfxG3_4);
		__m128 z1 = _mm_add_ps(_mm_sub_ps(z0, _mm_cvtepi32_ps(k1.m)), tfxG3_4);
		__m128 x2 = _mm_add_ps(_mm_sub_ps(x0, _mm_cvtepi32_ps(i2.m)), tfxG32_4);
		__m128 y2 = _mm_add_ps(_mm_sub_ps(y0, _mm_cvtepi32_ps(j2.m)), tfxG32_4);
		__m128 z2 = _mm_add_ps(_mm_sub_ps(z0, _mm_cvtepi32_ps(k2.m)), tfxG32_4);
		__m128 x3 = _mm_add_ps(_mm_sub_ps(x0, tfxONEF), tfxG33_4);
		__m128 y3 = _mm_add_ps(_mm_sub_ps(y0, tfxONEF), tfxG33_4);
		__m128 z3 = _mm_add_ps(_mm_sub_ps(z0, tfxONEF), tfxG33_4);

		// Work out the hashed gradient indices of the four simplex corners
		m128i ii;
		ii.m = _mm_and_si128(_mm_cvttps_epi32(i), tfxFF);
		m128i jj;
		jj.m = _mm_and_si128(_mm_cvttps_epi32(j), tfxFF);
		m128i kk;
		kk.m = _mm_and_si128(_mm_cvttps_epi32(k), tfxFF);
		m128i gi0, gi1, gi2, gi3;

		for (int i = 0; i < 4; ++i)
		{
			gi0.a[i] = permMOD12[ii.a[i] + perm[jj.a[i] + perm[kk.a[i]]]];
			gi1.a[i] = permMOD12[ii.a[i] + i1.a[i] + perm[jj.a[i] + j1.a[i] + perm[kk.a[i] + k1.a[i]]]];
			gi2.a[i] = permMOD12[ii.a[i] + i2.a[i] + perm[jj.a[i] + j2.a[i] + perm[kk.a[i] + k2.a[i]]]];
			gi3.a[i] = permMOD12[ii.a[i] + 1 + perm[jj.a[i] + 1 + perm[kk.a[i] + 1]]];
		}

		__m128 t0 = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(tfxPSIX, _mm_mul_ps(x0, x0)), _mm_mul_ps(y0, y0)), _mm_mul_ps(z0, z0));
		__m128 t1 = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(tfxPSIX, _mm_mul_ps(x1, x1)), _mm_mul_ps(y1, y1)), _mm_mul_ps(z1, z1));
		__m128 t2 = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(tfxPSIX, _mm_mul_ps(x2, x2)), _mm_mul_ps(y2, y2)), _mm_mul_ps(z2, z2));
		__m128 t3 = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(tfxPSIX, _mm_mul_ps(x3, x3)), _mm_mul_ps(y3, y3)), _mm_mul_ps(z3, z3));

		__m128 t0q = _mm_mul_ps(t0, t0);
		t0q = _mm_mul_ps(t0q, t0q);
		__m128 t1q = _mm_mul_ps(t1, t1);
		t1q = _mm_mul_ps(t1q, t1q);
		__m128 t2q = _mm_mul_ps(t2, t2);
		t2q = _mm_mul_ps(t2q, t2q);
		__m128 t3q = _mm_mul_ps(t3, t3);
		t3q = _mm_mul_ps(t3q, t3q);

		m128 gi0x, gi0y, gi0z, gi1x, gi1y, gi1z, gi2x, gi2y, gi2z, gi3x, gi3y, gi3z; 

		for (int i = 0; i < 4; i++)
		{
			gi0x.a[i] = gradX[gi0.a[i]];
			gi0y.a[i] = gradY[gi0.a[i]];
			gi0z.a[i] = gradZ[gi0.a[i]];

			gi1x.a[i] = gradX[gi1.a[i]];
			gi1y.a[i] = gradY[gi1.a[i]];
			gi1z.a[i] = gradZ[gi1.a[i]];

			gi2x.a[i] = gradX[gi2.a[i]];
			gi2y.a[i] = gradY[gi2.a[i]];
			gi2z.a[i] = gradZ[gi2.a[i]];

			gi3x.a[i] = gradX[gi3.a[i]];
			gi3y.a[i] = gradY[gi3.a[i]];
			gi3z.a[i] = gradZ[gi3.a[i]];

		}

		__m128 n0 = _mm_mul_ps(t0q, DotProductSIMD(gi0x.m, gi0y.m, gi0z.m, x0, y0, z0));
		__m128 n1 = _mm_mul_ps(t1q, DotProductSIMD(gi1x.m, gi1y.m, gi1z.m, x1, y1, z1));
		__m128 n2 = _mm_mul_ps(t2q, DotProductSIMD(gi2x.m, gi2y.m, gi2z.m, x2, y2, z2));
		__m128 n3 = _mm_mul_ps(t3q, DotProductSIMD(gi3x.m, gi3y.m, gi3z.m, x3, y3, z3));

		__m128 cond;

		cond = _mm_cmplt_ps(t0, tfxZERO);
		n0 = _mm_or_ps(_mm_andnot_ps(cond, n0), _mm_and_ps(cond, tfxZERO));
		cond = _mm_cmplt_ps(t1, tfxZERO);
		n1 = _mm_or_ps(_mm_andnot_ps(cond, n1), _mm_and_ps(cond, tfxZERO));
		cond = _mm_cmplt_ps(t2, tfxZERO);
		n2 = _mm_or_ps(_mm_andnot_ps(cond, n2), _mm_and_ps(cond, tfxZERO));
		cond = _mm_cmplt_ps(t3, tfxZERO);
		n3 = _mm_or_ps(_mm_andnot_ps(cond, n3), _mm_and_ps(cond, tfxZERO));

		tfxVec4 result;
		_mm_store_ps(&result.x, _mm_mul_ps(tfxTHIRTYTWO, _mm_add_ps(n0, _mm_add_ps(n1, _mm_add_ps(n2, n3)))));
		return result;
	}

	inline float simplex3d(float x, float y, float z)
	{
		static const float f3 = 1.0f / 3.0f;
		static const float g3 = 1.0f / 6.0f;
		static const float g32 = g3 * 2.f;
		static const float g33 = g3 * 3.f;

		float n0, n1, n2, n3; // Noise contributions from the four corners
							   // Skew the input space to determine which simplex cell we're in
		float s = (x + y + z)*f3; // Very nice and simple skew factor for 3D
		int i = fastfloor((x + s));
		int j = fastfloor(y + s);
		int k = fastfloor(z + s);
		float t = (i + j + k)*g3;
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
				i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 1; k2 = 0;
			} // X Y Z order
			else if (x0 >= z0) {
				i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 0; k2 = 1;
			} // X Z Y order
			else {
				i1 = 0; j1 = 0; k1 = 1; i2 = 1; j2 = 0; k2 = 1;
			} // Z X Y order
		}
		else { // x0<y0
			if (y0 < z0) {
				i1 = 0; j1 = 0; k1 = 1; i2 = 0; j2 = 1; k2 = 1;
			} // Z Y X order
			else if (x0 < z0) {
				i1 = 0; j1 = 1; k1 = 0; i2 = 0; j2 = 1; k2 = 1;
			} // Y Z X order
			else {
				i1 = 0; j1 = 1; k1 = 0; i2 = 1; j2 = 1; k2 = 0;
			} // Y X Z order
		}
		// A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
		// a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
		// a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
		// c = 1/6.
		float x1 = x0 - i1 + g3; // Offsets for second corner in (x,y,z) coords
		float y1 = y0 - j1 + g3;
		float z1 = z0 - k1 + g3;
		float x2 = x0 - i2 + g32; // Offsets for third corner in (x,y,z) coords
		float y2 = y0 - j2 + g32;
		float z2 = z0 - k2 + g32;
		float x3 = x0 - 1.0f + g33; // Offsets for last corner in (x,y,z) coords
		float y3 = y0 - 1.0f + g33;
		float z3 = z0 - 1.0f + g33;
		// Work out the hashed gradient indices of the four simplex corners
		int ii = i & 255;
		int jj = j & 255;
		int kk = k & 255;
		int gi0 = permMOD12[ii + perm[jj + perm[kk]]];
		int gi1 = permMOD12[ii + i1 + perm[jj + j1 + perm[kk + k1]]];
		int gi2 = permMOD12[ii + i2 + perm[jj + j2 + perm[kk + k2]]];
		int gi3 = permMOD12[ii + 1 + perm[jj + 1 + perm[kk + 1]]];
		// Calculate the contribution from the four corners
		float t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0;
		if (t0 < 0) n0 = 0.0f;
		else {
			t0 *= t0;
			n0 = t0 * t0 * Dot(gradX[gi0], gradY[gi0], gradZ[gi0], x0, y0, z0);
		}
		float t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1;
		if (t1 < 0) n1 = 0.0f;
		else {
			t1 *= t1;
			n1 = t1 * t1 * Dot(gradX[gi1], gradY[gi1], gradZ[gi1], x1, y1, z1);
		}
		float t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2;
		if (t2 < 0) n2 = 0.0f;
		else {
			t2 *= t2;
			n2 = t2 * t2 * Dot(gradX[gi2], gradY[gi2], gradZ[gi2], x2, y2, z2);
		}
		float t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3;
		if (t3 < 0) n3 = 0.0f;
		else {
			t3 *= t3;
			n3 = t3 * t3 * Dot(gradX[gi3], gradY[gi3], gradZ[gi3], x3, y3, z3);
		}
		// Add contributions from each corner to get the final noise value.	
		return 32.f * (n0 + n1 + n2 + n3);
	}


	/**
	 * Fractal/Fractional Brownian Motion (fBm) summation of 1D Perlin Simplex noise
	 *
	 * @param[in] octaves   number of fraction of noise to sum
	 * @param[in] x         float coordinate
	 *
	 * @return Noise value in the range[-1; 1], value of 0 on all integer coordinates.
	 */
	float tfxSimplexNoise::fractal(size_t octaves, float x) const {
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
	 * @param[in] v1.y         y float coordinate
	 *
	 * @return Noise value in the range[-1; 1], value of 0 on all integer coordinates.
	 */
	float tfxSimplexNoise::fractal(size_t octaves, float x, float y) const {
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
	float tfxSimplexNoise::fractal(size_t octaves, float x, float y, float z) const {
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

	//these Variables determine the timing resolution that particles are updated at. So an Update frequency of 60 would mean that the particles are updated at 60 frames per second.
	float tfxUPDATE_FREQUENCY = 60.f;
	float tfxUPDATE_TIME = 1.f / tfxUPDATE_FREQUENCY;
	float tfxFRAME_LENGTH = 1000.f / tfxUPDATE_FREQUENCY;

	//Set the udpate frequency for all particle effects - There may be options in the future for individual effects to be updated at their own specific frequency.
	void SetUpdateFrequency(float fps) {
		tfxUPDATE_FREQUENCY = fps;
		tfxUPDATE_TIME = 1.f / tfxUPDATE_FREQUENCY;
		tfxFRAME_LENGTH = 1000.f / tfxUPDATE_FREQUENCY;
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

	void tfxStr::Appendv(const char *format, va_list args) {
		va_list args_copy;
		va_copy(args_copy, args);

		int len = FormatString(NULL, 0, format, args);         // FIXME-OPT: could do a first pass write attempt, likely successful on first pass.
		if (len <= 0)
		{
			va_end(args_copy);
			return;
		}

		const int write_off = (current_size != 0) ? current_size : 1;
		const int needed_sz = write_off + len;
		if (write_off + (tfxU32)len >= capacity)
		{
			int new_capacity = capacity * 2;
			reserve(needed_sz > new_capacity ? needed_sz : new_capacity);
		}

		resize(needed_sz);
		FormatString(&data[write_off - 1], (size_t)len + 1, format, args);
		va_end(args_copy);

	}

	int tfxStr::Find(const char *needle) {
		tfxStr compare = needle;
		tfxStr lower = Lower();
		compare = compare.Lower();
		if (compare.Length() > Length()) return -1;
		tfxU32 pos = 0;
		int found = 0;
		while (compare.Length() + pos <= Length()) {
			if (strncmp(lower.data + pos, compare.data, compare.Length()) == 0) {
				return pos;
			}
			++pos;
		}
		return -1;
	}

	tfxStr tfxStr::Lower() {
		tfxStr convert = *this;
		for (auto &c : convert) {
			c = tolower(c);
		}
		return convert;
	}

	void tfxStr::AddLine(const char *format, ...) {
		va_list args;
		va_start(args, format);
		Appendv(format, args);
		va_end(args);
		Append('\n');
	}

	bool tfxStr::SaveToFile(const char *file_name) {
		FILE *file = fopen(file_name, "wb");
		if (!file)
			return false;

		tfxU32 l = Length();
		if (fwrite(data, 1, Length(), file) != Length())
			return false;

		fclose(file);
		return true;
	}

	const bool tfxStr::IsInt() const {
		if (!Length()) return false;
		for (auto c : *this) {
			if (c >= '0' && c <= '9');
			else {
				return false;
			}
		}
		return true;
	}

	const bool tfxStr::IsFloat() const {
		if (!Length()) return false;
		int dot_count = 0;
		for (auto c : *this) {
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

	void tfxStr::Appendf(const char *format, ...) {
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

		const int write_off = (current_size != 0) ? current_size : 1;
		const int needed_sz = write_off + len;
		if (write_off + (tfxU32)len >= capacity)
		{
			int new_capacity = capacity * 2;
			reserve(needed_sz > new_capacity ? needed_sz : new_capacity);
		}

		resize(needed_sz);
		FormatString(&strbuffer()[write_off - 1], (size_t)len + 1, format, args_copy);
		va_end(args_copy);

		va_end(args);
	}

	tfxStr512 tfxstream::ReadLine() {
		tfxStr512 line;
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

		tfxU64 length = _ftelli64(file);
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
		tfxU64 inventory_offset = sizeof(tfxHeader);
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
		fwrite((char*)&package.inventory.magic_number, 1, sizeof(tfxU32), file);
		fwrite((char*)&package.inventory.entry_count, 1, sizeof(tfxU32), file);
		for (auto &entry : package.inventory.entries.data) {
			fwrite((char*)&entry.file_name.current_size, 1, sizeof(tfxU32), file);
			fwrite(entry.file_name.c_str(), 1, entry.file_name.current_size, file);
			fwrite((char*)&entry.file_size, 1, sizeof(tfxU64), file);
			fwrite((char*)&entry.offset_from_start_of_file, 1, sizeof(tfxU64), file);
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
		tfxU64 inventory_offset = sizeof(tfxHeader);
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
		file.Write(&package.inventory.magic_number, sizeof(tfxU32));
		file.Write(&package.inventory.entry_count, sizeof(tfxU32));
		for (auto &entry : package.inventory.entries.data) {
			file.Write(&entry.file_name.current_size, sizeof(tfxU32));
			file.Write(entry.file_name.data, entry.file_name.current_size);
			file.Write(&entry.file_size, sizeof(tfxU64));
			file.Write(&entry.offset_from_start_of_file, sizeof(tfxU64));
		}

		return file;
	}

	tfxU64 GetPackageSize(tfxPackage &package) {
		tfxU64 space = 0;
		space += sizeof(tfxHeader);

		//Write the file contents
		for (auto &entry : package.inventory.entries.data) {
			space += entry.data.Size();
		}

		//Write the inventory
		space += sizeof(tfxU32);
		space += sizeof(tfxU32);
		for (auto &entry : package.inventory.entries.data) {
			space += sizeof(tfxU32);
			space += entry.file_name.current_size;
			space += sizeof(tfxU64);
			space += sizeof(tfxU64);
		}

		return space;
	}

	tfxErrorFlags LoadPackage(const char *file_name, tfxPackage &package) {

		package.file_data = ReadEntireFile(file_name);
		if (package.file_data.Size() == 0)
			return tfxErrorCode_unable_to_read_file;			//the file size is smaller then the expected header size

		package.file_size = package.file_data.Size();

		if (package.file_size < sizeof(tfxHeader))
			return tfxErrorCode_wrong_file_size;				//the file size is smaller then the expected header size

		package.file_data.Read((char*)&package.header, sizeof(tfxHeader));

		if (package.header.magic_number != tfxMAGIC_NUMBER)
			return tfxErrorCode_invalid_format;				//The header doesn't not contain the expected magic number "TFX!", incorrect file format;

		if (package.header.offset_to_inventory > package.file_size)
			return tfxErrorCode_no_inventory;				//The offset to the inventory is beyond the size of the file

		package.file_data.Seek(package.header.offset_to_inventory);
		package.file_data.Read((char*)&package.inventory.magic_number, sizeof(tfxU32));

		if (package.inventory.magic_number != tfxMAGIC_NUMBER_INVENTORY)
			return tfxErrorCode_invalid_inventory;			//The value at the inventory offset does not equal the expected magic number "INV!"

		package.file_data.Read((char*)&package.inventory.entry_count, sizeof(tfxU32));
		for (int i = 0; i != package.inventory.entry_count; ++i) {
			tfxEntryInfo entry;
			tfxU32 file_name_size;
			package.file_data.Read((char*)&file_name_size, sizeof(tfxU32));
			entry.file_name.resize(file_name_size);
			package.file_data.Read(entry.file_name.data, file_name_size);
			package.file_data.Read((char*)&entry.file_size, sizeof(tfxU64));
			package.file_data.Read((char*)&entry.offset_from_start_of_file, sizeof(tfxU64));
			package.inventory.entries.Insert(entry.file_name, entry);
		}

		package.file_path = file_name;

		return 0;
	}

	tfxErrorFlags LoadPackage(tfxstream &stream, tfxPackage &package) {
		//Note: tfxstream does not copy the memory, only the pointer, so if you FreeAll on the stream you pass in it will also free the file_data here as well
		package.file_data = stream;
		if (package.file_data.Size() == 0)
			return tfxErrorCode_unable_to_read_file;			//the file size is smaller then the expected header size

		package.file_size = package.file_data.Size();

		if (package.file_size < sizeof(tfxHeader))
			return tfxErrorCode_wrong_file_size;				//the file size is smaller then the expected header size

		package.file_data.Read((char*)&package.header, sizeof(tfxHeader));

		if (package.header.magic_number != tfxMAGIC_NUMBER)
			return tfxErrorCode_invalid_format;				//The header doesn't not contain the expected magic number "TFX!", incorrect file format;

		if (package.header.offset_to_inventory > package.file_size)
			return tfxErrorCode_no_inventory;				//The offset to the inventory is beyond the size of the file

		package.file_data.Seek(package.header.offset_to_inventory);
		package.file_data.Read((char*)&package.inventory.magic_number, sizeof(tfxU32));

		if (package.inventory.magic_number != tfxMAGIC_NUMBER_INVENTORY)
			return tfxErrorCode_invalid_inventory;			//The value at the inventory offset does not equal the expected magic number "INV!"

		package.file_data.Read((char*)&package.inventory.entry_count, sizeof(tfxU32));
		for (int i = 0; i != package.inventory.entry_count; ++i) {
			tfxEntryInfo entry;
			tfxU32 file_name_size;
			package.file_data.Read((char*)&file_name_size, sizeof(tfxU32));
			entry.file_name.resize(file_name_size);
			package.file_data.Read(entry.file_name.data, file_name_size);
			package.file_data.Read((char*)&entry.file_size, sizeof(tfxU64));
			package.file_data.Read((char*)&entry.offset_from_start_of_file, sizeof(tfxU64));
			package.inventory.entries.Insert(entry.file_name, entry);
		}

		return 0;
	}

	tfxEffectEmitter::~tfxEffectEmitter() {
	}

	void tfxEffectEmitter::SoftExpire() {
		flags |= tfxEmitterStateFlags_stop_spawning;
	}

	void tfxEffectEmitter::UpdateMaxLife() {
		GetInfo().max_life = GetMaxLife(*this);
		GetGraphByType(tfxOvertime_red)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_green)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_blue)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_blendfactor)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_intensity)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_velocity)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_width)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_height)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_weight)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_spin)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_stretch)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_spin)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_velocity_turbulance)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_direction_turbulance)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_velocity_adjuster)->lookup.life = GetInfo().max_life;
		GetGraphByType(tfxOvertime_direction)->lookup.life = GetInfo().max_life;
	}

	void tfxEffectEmitter::ResetAllBufferSizes() {
		tmpStack(tfxEffectEmitter*, stack);
		stack.push_back(this);
		while (!stack.empty()) {
			tfxEffectEmitter &current = *stack.pop_back();
			current.GetInfo().max_sub_emitters = 0;
			for (tfxEachLayer) {
				current.GetInfo().max_particles[layer] = 0;
			}
			for (auto &sub : current.GetInfo().sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	bool tfxEffectEmitter::IsFiniteEffect() {
		for (auto &e : GetInfo().sub_effectors) {
			if (e.common.property_flags & tfxEmitterPropertyFlags_single)
				return true;
			float qty = e.common.library->emitter_attributes[e.emitter_attributes].base.amount.GetLastValue() + e.common.library->emitter_attributes[e.emitter_attributes].variation.amount.GetLastValue();
			if (!(e.common.property_flags & tfxEmitterPropertyFlags_single) && qty > 0)
				return false;
		}
		return true;
	}

	void tfxEffectEmitter::FlagAs3D(bool flag) {
		if(flag)
			common.property_flags |= tfxEmitterPropertyFlags_is_3d;
		else
			common.property_flags &= ~tfxEmitterPropertyFlags_is_3d;
		for (auto &sub : GetInfo().sub_effectors) {
			sub.FlagAs3D(flag);
		}
	}

	bool tfxEffectEmitter::Is3DEffect() {
		return common.property_flags & tfxEmitterPropertyFlags_is_3d;
	}

	tfxParticleManagerModes tfxEffectEmitter::GetRequiredParticleManagerMode() {
		if (effect_flags & tfxEffectPropertyFlags_guaranteed_order && effect_flags & tfxEffectPropertyFlags_depth_draw_order) {
			return tfxParticleManagerMode_ordered_by_depth_guaranteed;
		}
		else if (effect_flags & tfxEffectPropertyFlags_depth_draw_order) {
			return tfxParticleManagerMode_ordered_by_depth;
		}
		else if (effect_flags & tfxEffectPropertyFlags_age_order) {
			return tfxParticleManagerMode_ordered_by_age;
		}
		return tfxParticleManagerMode_unordered;
	}

	void tfxEffectEmitter::UpdateAllBufferSizes() {
		tmpStack(tfxEffectEmitter*, stack);
		stack.push_back(this);
		for (int l = 0; l != tfxLAYERS; ++l) {
			GetInfo().max_particles[l] = 0;
		}
		GetInfo().max_sub_emitters = 0;
		while (!stack.empty()) {
			tfxEffectEmitter &current = *stack.pop_back();
			if (current.type == tfxEmitterType) {
				tfxU32 particle_count = 0;
				if (current.parent->parent && current.parent->parent->type != tfxFolder) {
					particle_count = current.GetHighestQty(current.parent->parent->GetInfo().max_life);
					GetInfo().max_sub_emitters += (current.parent->parent->GetInfo().max_particles[current.parent->parent->GetProperties().layer] * (current.GetInfo().sub_effectors.size() + 1) + (current.GetInfo().sub_effectors.size() + 1)) * tfxU32(current.GetInfo().max_sub_emitter_life / 1000.f);
				}
				else {
					particle_count = current.GetHighestQty(0.f);
				}
				current.GetInfo().max_particles[current.GetProperties().layer] = particle_count;
			}
			for (auto &sub : current.GetInfo().sub_effectors) {
				stack.push_back(&sub);
			}
		}
		//UpdateAllSpriteAmounts();
	}

	void tfxEffectEmitter::UpdateAllSpriteAmounts() {
		for (int layer = 0; layer != tfxLAYERS; ++layer) {
			for (auto &emitter : GetInfo().sub_effectors) {
				GetInfo().max_particles[layer] += emitter.GetInfo().max_particles[layer];
				if (emitter.GetInfo().sub_effectors.size()) {
					GetInfo().max_particles[layer] += emitter.GetSubEffectSpriteCounts(layer, emitter.GetInfo().max_particles[emitter.GetProperties().layer]);
				}
			}
		}
	}

	tfxU32 tfxEffectEmitter::GetSubEffectSpriteCounts(tfxU32 layer, tfxU32 multiplier) {
		float total = 0;
		for (auto &effect :GetInfo().sub_effectors) {
			for (auto &emitter : effect.GetInfo().sub_effectors) {
				if (common.property_flags & tfxEmitterPropertyFlags_single)
					total += emitter.GetInfo().max_particles[layer] * multiplier;
				else
					total += emitter.GetInfo().max_particles[layer] * multiplier * (GetInfo().max_sub_emitter_life / 1000.f);
				if (emitter.GetInfo().sub_effectors.size()) {
					if (common.property_flags & tfxEmitterPropertyFlags_single)
						total += emitter.GetSubEffectSpriteCounts(layer, emitter.GetInfo().max_particles[emitter.GetProperties().layer] * multiplier);
					else
						total += emitter.GetSubEffectSpriteCounts(layer, emitter.GetInfo().max_particles[emitter.GetProperties().layer] * multiplier * int(GetInfo().max_sub_emitter_life / 1000.f));
				}
			}
		}
		return (tfxU32)total;
	}

	float tfxEffectEmitter::GetSubEffectLength() {
		if (!parent) return 0.f;
		float max_life = 0.f;
		for (auto &e : GetInfo().sub_effectors) {
			max_life = std::fmaxf(GetMaxLife(e), max_life);
		}
		max_life += GetMaxLife(*parent);
		return max_life;
	}

	tfxU32 tfxEffectEmitter::GetHighestQty(float parent_age) {
		GetInfo().max_sub_emitter_life = 0.f;
		if (!(common.property_flags & tfxEmitterPropertyFlags_single)) {
			float max_qty = 0.f;
			float amount_remainder = 0.f;
			float highest_age = 0.f;
			if (common.library->emitter_attributes[emitter_attributes].base.amount.nodes.size() > 1 || common.library->emitter_attributes[emitter_attributes].variation.amount.nodes.size() > 1 || common.library->global_graphs[parent->global].amount.nodes.size() > 1
				|| common.library->emitter_attributes[emitter_attributes].base.life.nodes.size() > 1 || common.library->emitter_attributes[emitter_attributes].variation.life.nodes.size() > 1 || common.library->global_graphs[parent->global].life.nodes.size() > 1) {
				tmpStack(tfxVec2, particles);
				float max_frames = std::fmaxf(common.library->emitter_attributes[emitter_attributes].base.amount.GetLastFrame(), common.library->emitter_attributes[emitter_attributes].variation.amount.GetLastFrame());
				max_frames = std::fmaxf(max_frames, common.library->emitter_attributes[emitter_attributes].base.life.GetLastFrame());
				max_frames = std::fmaxf(max_frames, common.library->emitter_attributes[emitter_attributes].variation.life.GetLastFrame());
				max_frames = std::fmaxf(max_frames, common.library->global_graphs[parent->global].amount.GetLastFrame());
				max_frames = std::fmaxf(max_frames, common.library->global_graphs[parent->global].life.GetLastFrame());
				max_frames = std::fmaxf(max_frames, GetProperties().loop_length);
				float max_last_life = common.library->emitter_attributes[emitter_attributes].base.life.GetLastValue() + common.library->emitter_attributes[emitter_attributes].variation.life.GetLastValue();
				float current_frame = 0.f;
				unsigned start_index = 0;
				if (parent_age > 0.f)
					max_frames = std::fminf(max_frames, parent_age);
				for (float frame = 0; frame <= max_frames + max_last_life; frame += tfxFRAME_LENGTH) {
					float qty;
					qty = LookupFast(common.library->emitter_attributes[emitter_attributes].base.amount, current_frame);
					qty += LookupFast(common.library->emitter_attributes[emitter_attributes].variation.amount, current_frame);
					float life = LookupFast(common.library->emitter_attributes[emitter_attributes].base.life, current_frame) + lookup_callback(common.library->emitter_attributes[emitter_attributes].variation.life, current_frame);

					if (common.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && (GetProperties().emission_type == tfxArea || GetProperties().emission_type == tfxEllipse)) {
						float area = LookupFast(common.library->emitter_attributes[emitter_attributes].properties.emitter_width, current_frame) * LookupFast(common.library->emitter_attributes[emitter_attributes].properties.emitter_height, current_frame);
						qty = (qty / 10000.f) * area;
					}
					else if (common.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && GetProperties().emission_type == tfxLine) {
						qty = (qty / 100.f) * LookupFast(common.library->emitter_attributes[emitter_attributes].properties.emitter_height, current_frame);
					}

					qty *= LookupFast(parent->common.library->global_graphs[parent->global].amount, current_frame);
					qty *= tfxUPDATE_TIME;
					qty += amount_remainder;

					float total_qty = 0;
					for (auto &p : particles) {
						if (parent_age && current_frame <= parent_age)
							highest_age = std::fmaxf(highest_age, p.y);
						p.y -= tfxFRAME_LENGTH;
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
					current_frame += tfxFRAME_LENGTH / tfxLOOKUP_FREQUENCY;
					if (GetProperties().loop_length > 0 && current_frame > GetProperties().loop_length)
						current_frame = 0;
				}

				float total_qty = 0;
				for (auto &p : particles) {
					p.y -= tfxFRAME_LENGTH;
					if (parent_age && current_frame <= parent_age)
						highest_age = std::fmaxf(highest_age, p.y);
					if (p.y < 0) p.x = 0;
					total_qty += p.x;
				}
				max_qty = std::fmaxf(max_qty, total_qty);
				GetInfo().max_sub_emitter_life = highest_age + current_frame;

				return (tfxU32)max_qty;
			}
			else {
				float max_life = GetMaxLife(*this);
				GetInfo().max_sub_emitter_life = max_life + parent_age;
				float qty = (common.library->emitter_attributes[emitter_attributes].base.amount.GetFirstValue() + common.library->emitter_attributes[emitter_attributes].variation.amount.GetFirstValue()) * (max_life / 1000.f) + 1.f;
				if (common.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && (GetProperties().emission_type == tfxArea || GetProperties().emission_type == tfxEllipse)) {
					float area = common.library->emitter_attributes[emitter_attributes].properties.emitter_width.GetFirstValue() * common.library->emitter_attributes[emitter_attributes].properties.emitter_height.GetFirstValue();
					qty = (qty / 10000.f) * area;
				}
				else if (common.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && GetProperties().emission_type == tfxLine) {
					qty = (qty / 100.f) * common.library->emitter_attributes[emitter_attributes].properties.emitter_height.GetFirstValue();
				}
				return tfxU32(qty);
			}
		}
		else {
			GetInfo().max_sub_emitter_life = GetInfo().max_life + parent_age;
			return (tfxU32)GetProperties().spawn_amount;
		}
	}

	void tfxParticleMemoryTools::GetEffectMaxFrames(tfxEffectEmitter &effect) {
		if (effect.type == tfxEffectType && effect.IsRootEffect()) {
			max_frames = std::fmaxf(max_frames, effect.common.library->global_graphs[effect.global].amount.GetLastFrame());
			max_frames = std::fmaxf(max_frames, effect.common.library->global_graphs[effect.global].life.GetLastFrame());
			max_frames = std::fmaxf(max_frames, effect.GetProperties().loop_length);
		}
		else if (effect.type == tfxEmitterType) {
			max_frames = std::fmaxf(max_frames, effect.common.library->emitter_attributes[effect.emitter_attributes].base.amount.GetLastFrame());
			max_frames = std::fmaxf(max_frames, effect.common.library->emitter_attributes[effect.emitter_attributes].variation.amount.GetLastFrame());
			max_frames = std::fmaxf(max_frames, effect.common.library->emitter_attributes[effect.emitter_attributes].base.life.GetLastFrame());
			max_frames = std::fmaxf(max_frames, effect.common.library->emitter_attributes[effect.emitter_attributes].variation.life.GetLastFrame());
			max_frames = std::fmaxf(max_frames, effect.GetProperties().loop_length);
			max_last_life = std::fmaxf(max_last_life, effect.common.library->emitter_attributes[effect.emitter_attributes].base.life.GetLastValue() + effect.common.library->emitter_attributes[effect.emitter_attributes].variation.life.GetLastValue());
		}
		for (auto &sub : effect.GetInfo().sub_effectors) {
			GetEffectMaxFrames(sub);
		}
	}

	void tfxParticleMemoryTools::ProcessEffect(tfxEffectEmitter &effect) {
		effects[0].reserve(10000);
		effects[1].reserve(10000);
		effects[0].clear();
		effects[1].clear();
		for (tfxU32 layer = 0; layer != tfxLAYERS; ++layer) {
			particles[layer][0].clear();
			particles[layer][1].clear();
			sprite_count[layer] = 0;
		}
		sub_effect_count = 0;
		emitters_removed = 0;
		GetEffectMaxFrames(effect);
		AddEffect(effect);
		initial_effect_size = effect.GetInfo().sub_effectors.size() + 1;
		Process();
		for (tfxEachLayer) {
			effect.GetInfo().max_particles[layer] = sprite_count[layer];
			effect.GetInfo().max_sub_emitters = sub_effect_count;
		}
	}

	void tfxParticleMemoryTools::Process() {
		for (tfxU32 layer = 0; layer != tfxLAYERS; ++layer) {
			sprite_count[layer] = 0;
		}
		for (float frame = 0; frame <= max_frames + max_last_life; frame += tfxFRAME_LENGTH) {
			current_effect.highest_particle_age += tfxFRAME_LENGTH;
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
			for (tfxEachLayer) {
				sprite_count[layer] = (tfxU32)std::fmaxf((float)sprite_count[layer], (float)particles[layer][!current_buffer].current_size);
			}
			sub_effect_count = (tfxU32)std::fmaxf((float)effects[!current_buffer].current_size, (float)sub_effect_count);

			effects[current_buffer].clear();
			current_buffer = !current_buffer;
		}
		sub_effect_count -= initial_effect_size;
	}

	void tfxParticleMemoryTools::MockUpdateEmitter(tfxMockEffect &emitter) {
		if (emitter.library_link->type == tfxEmitterType) {

			for (auto p : emitter.particles[current_buffer]) {
				p -= tfxFRAME_LENGTH;
				if (p > 0) emitter.particles[!current_buffer].push_back(p);
			}
			tfxEmitterProperties &properties = emitter.library_link->GetProperties();
			emitter.library_link->GetInfo().max_particles[properties.layer] = (tfxU32)std::fmaxf((float)emitter.library_link->GetInfo().max_particles[properties.layer], (float)emitter.particles[current_buffer].current_size);
			emitter.particles[current_buffer].clear();

			float life = LookupFast(emitter.library->emitter_attributes[emitter.library_link->emitter_attributes].base.life, emitter.frame) + lookup_callback(emitter.library->emitter_attributes[emitter.library_link->emitter_attributes].variation.life, emitter.frame) * 0.75f;
			if (!(emitter.library_link->common.property_flags & tfxEmitterPropertyFlags_single)) {
				emitter.qty = LookupFast(emitter.library->emitter_attributes[emitter.library_link->emitter_attributes].base.amount, emitter.frame);
				emitter.qty += LookupFast(emitter.library->emitter_attributes[emitter.library_link->emitter_attributes].variation.amount, emitter.frame);

				if (emitter.library_link->common.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && (properties.emission_type == tfxArea || properties.emission_type == tfxEllipse)) {
					float area = LookupFast(emitter.library->emitter_attributes[emitter.library_link->emitter_attributes].properties.emitter_width, emitter.frame) * LookupFast(emitter.library->emitter_attributes[emitter.library_link->emitter_attributes].properties.emitter_height, emitter.frame);
					emitter.qty = (emitter.qty / 10000.f) * area;
				}
				else if (emitter.library_link->common.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && properties.emission_type == tfxLine) {
					emitter.qty = (emitter.qty / 100.f) * LookupFast(emitter.library->emitter_attributes[emitter.library_link->emitter_attributes].properties.emitter_height, emitter.frame);
				}

				emitter.qty *= LookupFast(emitter.library->global_graphs[emitter.library_link->parent->global].amount, emitter.frame);
				emitter.qty *= tfxUPDATE_TIME;
				emitter.qty += emitter.amount_remainder;

				if (!emitter.started_spawning && emitter.qty < 1)
					emitter.qty = 1.f;

				if (emitter.qty >= 1) {
					while (emitter.qty > 1) {
						emitter.qty -= 1.f;
						particles[properties.layer][!current_buffer].push_back(life + tfxFRAME_LENGTH);
						emitter.particles[!current_buffer].push_back(life + tfxFRAME_LENGTH);
						for (auto &sub : emitter.library_link->GetInfo().sub_effectors) {
							AddEffect(sub);
						}
					}
					emitter.started_spawning = true;
					emitter.highest_particle_age = std::fmaxf(emitter.highest_particle_age, life);
				}
				emitter.amount_remainder = emitter.qty;
				if (properties.loop_length > 0 && emitter.frame > properties.loop_length)
					emitter.frame = 0;
			}
			else if(!emitter.single_shot_done) {
				emitter.started_spawning = true;
				emitter.library_link->GetInfo().max_particles[properties.layer] = properties.spawn_amount;
				for (int q = 0; q != properties.spawn_amount; ++q) {
					if (emitter.library_link->common.property_flags & tfxEmitterPropertyFlags_single) {
						particles[properties.layer][!current_buffer].push_back(max_frames + max_last_life);
						emitter.particles[!current_buffer].push_back(max_frames + max_last_life);
					}
					else {
						particles[properties.layer][!current_buffer].push_back(life + tfxFRAME_LENGTH);
						emitter.particles[!current_buffer].push_back(life + tfxFRAME_LENGTH);
					}
					emitter.highest_particle_age = std::fmaxf(emitter.highest_particle_age, life);
					for (auto &sub : emitter.library_link->GetInfo().sub_effectors) {
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

		emitter.age += tfxFRAME_LENGTH;
		emitter.frame += tfxFRAME_LENGTH / tfxLOOKUP_FREQUENCY;
		emitter.highest_particle_age -= tfxFRAME_LENGTH;

		if (emitter.highest_particle_age <= 0 && emitter.started_spawning && emitter.qty == 0.f) {
			emitter.timeout_counter++;
		}
		else {
			emitter.timeout_counter = 0;
		}
	}

	void tfxParticleMemoryTools::MockUpdateParticles() {
		for (tfxU32 layer = 0; layer != tfxLAYERS; ++layer) {
			for (auto p : particles[layer][current_buffer]) {
				p -= tfxFRAME_LENGTH;
				if (p > 0) particles[layer][!current_buffer].push_back(p);
			}
			particles[layer][current_buffer].clear();
		}
	}

	void tfxParticleMemoryTools::AddEffect(tfxEffectEmitter &effect) {
		tfxMockEffect new_effect;
		new_effect.library_link = effect.common.library->GetEffect(effect.path_hash);
		new_effect.library = effect.common.library;
		new_effect.emitter_count = effect.GetInfo().sub_effectors.size();
		effects[!current_buffer].push_back(new_effect);
		for (auto &emitter : effect.GetInfo().sub_effectors) {
			tfxMockEffect new_emitter;
			new_emitter.library_link = effect.common.library->GetEffect(emitter.path_hash);
			new_emitter.library = effect.common.library;
			effects[!current_buffer].push_back(new_emitter);
		}
	}

	tfxEffectEmitter& tfxEffectEmitter::AddEmitter(tfxEffectEmitter &e) {
		assert(e.GetInfo().name.Length());				//Emitter must have a name so that a hash can be generated
		e.type = tfxEffectEmitterType::tfxEmitterType;
		e.common.library = common.library;
		e.GetInfo().uid = ++common.library->uid;
		GetInfo().sub_effectors.push_back(e);
		common.library->UpdateEffectPaths();
		ReIndex();
		return GetInfo().sub_effectors.back();
	}

	tfxEffectEmitter& tfxEffectEmitter::AddEffect(tfxEffectEmitter &e) {
		assert(e.GetInfo().name.Length());				//Effect must have a name so that a hash can be generated
		e.type = tfxEffectEmitterType::tfxEffectType;
		e.common.library = common.library;
		e.parent = this;
		e.GetInfo().uid = ++common.library->uid;
		GetInfo().sub_effectors.push_back(e);
		common.library->UpdateEffectPaths();
		ReIndex();
		return GetInfo().sub_effectors.back();
	}

	tfxEffectEmitter& tfxEffectEmitter::AddEffect() {
		tfxEffectEmitter e;
		e.common.library = common.library;
		e.GetInfo().uid = ++common.library->uid;
		e.type = tfxEffectEmitterType::tfxEffectType;
		e.GetInfo().name = "New Effect";
		GetInfo().sub_effectors.push_back(e);
		common.library->UpdateEffectPaths();
		ReIndex();
		return GetInfo().sub_effectors.back();
	}

	tfxEffectEmitter &tfxEffectEmitter::AddEffector(tfxEffectEmitterType type) {
		tfxEffectEmitter e;
		//e.parent_effect = this;
		e.type = type;
		e.common.library = common.library;
		e.GetInfo().uid = ++common.library->uid;
		if(e.type == tfxEffectType)
			e.GetInfo().name = "New Effect";
		else
			e.GetInfo().name = "New Emitter";
		GetInfo().sub_effectors.push_back(e);
		common.library->UpdateEffectPaths();
		ReIndex();
		return GetInfo().sub_effectors.back();
	}

	void tfxEffectEmitter::NoTweenNextUpdate() {
		tmpStack(tfxEffectEmitter*, stack);
		stack.push_back(this);
		while (!stack.empty()) {
			auto &current = stack.pop_back();
			current->flags |= tfxEmitterStateFlags_no_tween_this_update;
			for (auto &sub : current->GetInfo().sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	float GetEmissionDirection2d(tfxCommon &common, tfxEmitterState &current, tfxEffectEmitter *library_link, tfxVec2 local_position, tfxVec2 world_position, tfxVec2 emitter_size) {
		//float (*effect_lookup_callback)(tfxGraph &graph, float age) = common.root_effect->lookup_mode == tfxPrecise ? LookupPrecise : LookupFast;
		float emission_angle = lookup_callback(common.library->emitter_attributes[library_link->emitter_attributes].properties.emission_pitch, common.frame);
		float emission_angle_variation = lookup_callback(common.library->emitter_attributes[library_link->emitter_attributes].properties.emission_range, common.frame);
		//----Emission
		float range = emission_angle_variation *.5f;
		float direction = 0;
		tfxEmitterProperties &properties = library_link->GetProperties();

		if (properties.emission_type == tfxEmissionType::tfxPoint)
			return direction + emission_angle + random_generation.Range(-range, range);

		tfxVec2 tmp_position;
		if (common.handle.x + local_position.x == 0 && common.handle.y + local_position.y == 0)
			tmp_position = emitter_size;
		else
			tmp_position = local_position + common.handle.xy();

		if (properties.emission_direction == tfxEmissionDirection::tfxOutwards) {

			tfxVec2 to_handle;

			if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = tmp_position;
			else
				to_handle = world_position - common.transform.world_position.xy();

			direction = GetVectorAngle(to_handle.x, to_handle.y);

		}
		else if (properties.emission_direction == tfxEmissionDirection::tfxInwards) {

			tfxVec2 to_handle;

			if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = (-tmp_position);

			else
				to_handle = (common.transform.world_position.xy() - world_position - common.handle.xy());

			direction = GetVectorAngle(to_handle.x, to_handle.y);

		}
		else if (properties.emission_direction == tfxEmissionDirection::tfxBothways) {

			if (current.emission_alternator) {

				tfxVec2 to_handle;

				if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = (tmp_position);
				else
					to_handle = (world_position - common.transform.world_position.xy() - common.handle.xy());

				direction = GetVectorAngle(to_handle.x, to_handle.y);

			}
			else {

				tfxVec2 to_handle;

				if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = (-tmp_position);
				else
					to_handle = (common.transform.world_position.xy() - world_position);

				direction = GetVectorAngle(to_handle.x, to_handle.y);

			}

			current.emission_alternator = !current.emission_alternator;
		}

		if (std::isnan(direction))
			direction = 0.f;
		return direction + emission_angle + random_generation.Range(-range, range);
	}

	tfxVec3 GetEmissionDirection3d(tfxCommon &common, tfxEmitterState &current, tfxEffectEmitter *library_link, float emission_pitch, float emission_yaw, tfxVec3 local_position, tfxVec3 world_position, tfxVec3 emitter_size) {
		//float (*effect_lookup_callback)(tfxGraph &graph, float age) = common.root_effect->lookup_mode == tfxPrecise ? LookupPrecise : LookupFast;
		float emission_angle_variation = lookup_callback(common.library->emitter_attributes[library_link->emitter_attributes].properties.emission_range, common.frame);
		//----Emission
		float range = emission_angle_variation * .5f;

		tfxVec3 result;
		tfxVec3 tmp_position;
		if (common.handle.x + local_position.x == 0 && common.handle.y + local_position.y == 0)
			tmp_position = emitter_size;
		else
			tmp_position = local_position + common.handle;

		tfxEmitterProperties &properties = library_link->GetProperties();

		tfxVec3 to_handle(0.f, 1.f, 0.f);
		float parent_pitch = 0.f;
		float parent_yaw = 0.f;
		if (properties.emission_type != tfxEmissionType::tfxPoint) {
			if (properties.emission_direction == tfxEmissionDirection::tfxOutwards) {

				if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = tmp_position;
				else
					to_handle = world_position - common.transform.world_position;

				to_handle = FastNormalizeVec(to_handle);

			}
			else if (properties.emission_direction == tfxEmissionDirection::tfxInwards) {

				if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = -tmp_position;
				else
					to_handle = common.transform.world_position - world_position;

				to_handle = FastNormalizeVec(to_handle);

			}
			else if (properties.emission_direction == tfxEmissionDirection::tfxBothways) {

				if (current.emission_alternator) {

					if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
						to_handle = tmp_position;
					else
						to_handle = world_position - common.transform.world_position;
				}
				else {

					if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
						to_handle = -tmp_position;
					else
						to_handle = common.transform.world_position - world_position;
				}

				current.emission_alternator = !current.emission_alternator;
				to_handle = FastNormalizeVec(to_handle);
			}
			else {
				parent_pitch = common.transform.world_rotations.pitch;
				parent_yaw = common.transform.world_rotations.yaw;
			}
		}
		else {
			parent_pitch = common.transform.world_rotations.pitch;
			parent_yaw = common.transform.world_rotations.yaw;
		}

		float pitch = asinf(-to_handle.y);
		float yaw = atan2(to_handle.x, to_handle.z);
		tfxVec3 direction;
		direction.z = cos(emission_yaw + yaw + parent_yaw) * cos(emission_pitch + pitch + parent_pitch);
		direction.y = -sin(emission_pitch + pitch + parent_pitch);
		direction.x = sin(emission_yaw + yaw + parent_yaw) * cos(emission_pitch + pitch + parent_pitch);
		tfxVec3 v = direction;
		if (range != 0) {
			result.y = random_generation.Range(1.f) * (1.f - cos(range)) + cos(range);
			float phi = random_generation.Range(1.f) * 2.f * tfxPI;
			float s = sqrt(1.f - (result.y * result.y));
			result.x = s * cos(phi);
			result.z = s * sin(phi);

			v = result;

			if (direction.y != 1.f) {
				tfxVec3 u = Cross(tfxVec3(0.f, 1.f, 0.f), direction);
				float rot = acosf(DotProduct(direction, tfxVec3(0.f, 1.f, 0.f)));
				tfxMatrix4 handle_mat = M4();
				handle_mat = mmRotate(handle_mat, rot, u);
				v = mmTransformVector(handle_mat, result).xyz();
				v.x = -v.x;
				v.z = -v.z;
			}
		}

		return v;
	}

	void InitialiseParticle2d(tfxParticleData &data, tfxSpriteTransform2d &sprite_transform, tfxEmitterState &current, tfxCommon &common, tfxEmitterSpawnControls &spawn_values, tfxEffectEmitter *library_link, float tween) {
		//----Position
		data.local_position = 0;
		sprite_transform.position = 0;
		data.captured_position = 0;
		tfxEmitterProperties &properties = library_link->GetProperties();
		tfxVec2 lerp_position = InterpolateVec2(tween, common.transform.captured_position.xy(), common.transform.world_position.xy());
		if (properties.emission_type == tfxEmissionType::tfxPoint) {
			if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
				data.local_position = 0;
			else {
				if (common.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) {
					data.local_position = lerp_position;
				}
				else {
					tfxVec2 rotvec = mmTransformVector(common.transform.matrix, -common.handle.xy());
					data.local_position = rotvec + lerp_position;
				}
			}
		}
		else if (properties.emission_type == tfxEmissionType::tfxArea) {
			tfxVec2 position = tfxVec2(0.f, 0.f);

			if (common.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				if (common.property_flags & tfxEmitterPropertyFlags_fill_area) {
					if (!(common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
						current.grid_coords.x--;
						if (current.grid_coords.x < 0.f) {
							current.grid_coords.y--;
							current.grid_coords.x = properties.grid_points.x - 1;
							if (current.grid_coords.y < 0.f)
								current.grid_coords.y = properties.grid_points.y - 1;
						}
					}

					data.local_position = position + (current.grid_coords.xy() * spawn_values.grid_segment_size.xy());

					if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
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

					if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {

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
					tfxBound(current.grid_coords.xy(), properties.grid_points.xy());
					data.local_position = position + (current.grid_coords.xy() * spawn_values.grid_segment_size.xy());
				}
			}
			else {
				if (common.property_flags & tfxEmitterPropertyFlags_fill_area) {
					position.x = random_generation.Range(current.emitter_size.x);
					position.y = random_generation.Range(current.emitter_size.y);
				}
				else {
					//Spawn on one of 4 edges of the area
					tfxU32 side = random_generation.RangeUInt(4);
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

				data.local_position = position;
			}

			//----TForm and Emission
			if (!(common.property_flags & tfxEmitterPropertyFlags_relative_position)) {
				data.local_position = mmTransformVector(common.transform.matrix, data.local_position.xy() + common.handle.xy());
				data.local_position = lerp_position + data.local_position.xy() * common.transform.scale.xy();
			}

		}
		else if (properties.emission_type == tfxEmissionType::tfxEllipse) {
			tfxVec2 emitter_size = (current.emitter_size.xy() * .5f);
			tfxVec2 position = tfxVec2(0.f, 0.f);

			if (common.property_flags & tfxEmitterPropertyFlags_spawn_on_grid && !(common.property_flags & tfxEmitterPropertyFlags_fill_area)) {

				current.grid_coords.y = 0.f;

				if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					current.grid_coords.x--;
					if (current.grid_coords.x < 0.f) {
						current.grid_coords.x = properties.grid_points.x - 1;
					}
				}

				float th = current.grid_coords.x * spawn_values.grid_segment_size.x + spawn_values.arc_offset;
				data.local_position = tfxVec2(std::cosf(th) * emitter_size.x + emitter_size.x, -std::sinf(th) * emitter_size.y + emitter_size.y);

				if (!(common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					current.grid_coords.x++;
					if (current.grid_coords.x >= properties.grid_points.x) {
						current.grid_coords.x = 0.f;
					}
				}

			}
			else if (!(common.property_flags & tfxEmitterPropertyFlags_fill_area)) {
				float th = random_generation.Range(spawn_values.arc_size) + spawn_values.arc_offset;

				data.local_position = tfxVec2(std::cosf(th) * emitter_size.x + emitter_size.x, -std::sinf(th) * emitter_size.y + emitter_size.y);

			}
			else {
				data.local_position.x = random_generation.Range(0, current.emitter_size.x);
				data.local_position.y = random_generation.Range(0, current.emitter_size.y);

				while ((std::pow(data.local_position.x - emitter_size.x, 2) / std::pow(emitter_size.x, 2)) + (std::pow(data.local_position.y - emitter_size.y, 2) / std::pow(emitter_size.y, 2)) > 1) {
					data.local_position.x = random_generation.Range(0, current.emitter_size.x);
					data.local_position.y = random_generation.Range(0, current.emitter_size.y);
				}
			}

			//----TForm and Emission
			if (!(common.property_flags & tfxEmitterPropertyFlags_relative_position)) {
				data.local_position = mmTransformVector(common.transform.matrix, data.local_position.xy() + common.handle.xy());
				data.local_position = lerp_position + data.local_position.xy() * common.transform.scale.xy();
			}

		}
		else if (properties.emission_type == tfxEmissionType::tfxLine) {
			if (common.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				current.grid_coords.x = 0.f;

				if (!(common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					current.grid_coords.y--;
					if (current.grid_coords.y < 0.f) {
						current.grid_coords.y = properties.grid_points.x - 1;
					}
				}

				data.local_position = tfxVec2(current.grid_coords.xy() * -spawn_values.grid_segment_size.xy());

				if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					current.grid_coords.y++;
					if (current.grid_coords.y >= properties.grid_points.x) {
						current.grid_coords.y = 0.f;
					}
				}

			}
			else {
				data.local_position.x = 0.f;
				data.local_position.y = random_generation.Range(-current.emitter_size.y, 0.f);

			}

			//----TForm and Emission
			if (!(common.property_flags & tfxEmitterPropertyFlags_relative_position) && !(common.property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
				data.local_position = mmTransformVector(common.transform.matrix, data.local_position.xy() + common.handle.xy());
				data.local_position = lerp_position + data.local_position.xy() * common.transform.scale.xy();
			}
		}

		//----Weight
		if (spawn_values.weight) {
			data.base.weight = spawn_values.weight;
			if (spawn_values.weight_variation > 0) {
				data.base.weight += random_generation.Range(-spawn_values.weight_variation, spawn_values.weight_variation);
			}
		}
		else {
			data.base.weight = 0;
		}
		data.weight_acceleration = data.base.weight * common.library->emitter_attributes[library_link->emitter_attributes].overtime.weight.GetFirstValue() * tfxUPDATE_TIME;

		//----Velocity
		float velocity_scale = common.library->emitter_attributes[library_link->emitter_attributes].overtime.velocity.GetFirstValue() * current.velocity_adjuster;
		data.base.velocity = spawn_values.velocity + random_generation.Range(-spawn_values.velocity_variation, spawn_values.velocity_variation);

		//----Size
		if (!(common.property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			float random_size_x = random_generation.Range(spawn_values.size_variation.x);
			float random_size_y = random_generation.Range(spawn_values.size_variation.y);
			data.base.size.y = random_size_y + spawn_values.size.y;
			data.base.size.x = (random_size_x + spawn_values.size.x) / properties.image->image_size.x;
			float height = data.base.size.y / properties.image->image_size.y;

			sprite_transform.scale.x = data.base.size.x * common.library->emitter_attributes[library_link->emitter_attributes].overtime.width.GetFirstValue();

			if (common.property_flags & tfxEmitterPropertyFlags_lifetime_uniform_size) {
				sprite_transform.scale.y = height * common.library->emitter_attributes[library_link->emitter_attributes].overtime.width.GetFirstValue();
			}
			else {
				sprite_transform.scale.y = height * common.library->emitter_attributes[library_link->emitter_attributes].overtime.height.GetFirstValue();
			}
		}
		else {
			float random_size_x = random_generation.Range(spawn_values.size_variation.x);
			float random_size_y = random_size_x;
			data.base.size.y = random_size_y + spawn_values.size.y;
			data.base.size.x = (random_size_x + spawn_values.size.x) / properties.image->image_size.x;
			float height = data.base.size.y / properties.image->image_size.y;

			sprite_transform.scale.x = data.base.size.x * common.library->emitter_attributes[library_link->emitter_attributes].overtime.width.GetFirstValue();
			sprite_transform.scale.y = sprite_transform.scale.x;
		}

		sprite_transform.scale *= current.overal_scale;

		//----Spin
		data.base.spin = random_generation.Range(-spawn_values.spin_variation, std::abs(spawn_values.spin_variation)) + spawn_values.spin;

		sprite_transform.rotation = 0;
		data.local_rotations = 0;
		if (properties.angle_settings & tfxAngleSettingFlags_random_roll) {
			sprite_transform.rotation = data.local_rotations.roll = random_generation.Range(properties.angle_offsets.roll);
		}
		else if (properties.angle_settings & tfxAngleSettingFlags_specify_roll) {
			sprite_transform.rotation = data.local_rotations.roll = properties.angle_offsets.roll;
		}else{
			sprite_transform.rotation = data.local_rotations.roll = 0;
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
				data.local_position.x += splatx * common.transform.scale.x;
				data.local_position.y += splaty * common.transform.scale.y;
			}
			else {
				data.local_position.x += splatx;
				data.local_position.y += splaty;
			}
		}

		float direction = 0;

		if (properties.angle_settings & tfxAngleSettingFlags_align_roll && common.property_flags & tfxEmitterPropertyFlags_edge_traversal)
			sprite_transform.rotation = data.local_rotations.roll = direction + properties.angle_offsets.roll;

		bool line = common.property_flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type == tfxEmissionType::tfxLine;

		if (!line && !(common.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			current.transform_particle_callback2d(data, sprite_transform.position, sprite_transform.rotation, common, common.transform.world_position);
			data.captured_position = sprite_transform.captured_position = sprite_transform.position;
		}

		if (!line) {
			direction = data.velocity_normal.x = GetEmissionDirection2d(common, current, library_link, data.local_position.xy(), sprite_transform.position, current.emitter_size.xy()) + common.library->emitter_attributes[library_link->emitter_attributes].overtime.direction.GetFirstValue();
		}
		//Do a micro update
		float micro_time = tfxUPDATE_TIME * (1.f - tween);
		data.weight_acceleration += data.base.weight * common.library->emitter_attributes[library_link->emitter_attributes].overtime.weight.GetFirstValue() * micro_time;
		//----Velocity Changes
		tfxVec2 velocity_normal;
		velocity_normal.x = std::sinf(direction);
		velocity_normal.y = -std::cosf(direction);
		tfxVec2 current_velocity = (data.base.velocity * common.library->emitter_attributes[library_link->emitter_attributes].overtime.velocity.GetFirstValue()) * velocity_normal;
		current_velocity.y += data.weight_acceleration;
		current_velocity *= micro_time;
		data.local_position += current_velocity;
		if (line || common.property_flags & tfxEmitterPropertyFlags_relative_position) {
			tfxVec2 rotatevec = mmTransformVector(common.transform.matrix, tfxVec2(data.local_position.x, data.local_position.y) + common.handle.xy());
			data.captured_position = sprite_transform.captured_position = common.transform.captured_position.xy() + rotatevec * common.transform.scale.xy();
			current.transform_particle_callback2d(data, sprite_transform.position, sprite_transform.rotation, common, tfxVec3(common.transform.world_position.x, common.transform.world_position.y, 0.f));
		}
		else {
			sprite_transform.position += current_velocity;
		}
		//end micro update

		//data.velocity = data.velocity_normal * data.base.velocity * data.velocity_scale * tfxUPDATE_TIME;

		if ((properties.angle_settings & tfxAngleSettingFlags_align_roll || properties.angle_settings & tfxAngleSettingFlags_align_with_emission) && !line) {
			//----Normalize Velocity to direction
			tfxVec2 velocity_normal;
			velocity_normal.x = std::sinf(direction);
			velocity_normal.y = -std::cosf(direction);
			sprite_transform.rotation = data.local_rotations.roll = GetVectorAngle(velocity_normal.x, velocity_normal.y) + properties.angle_offsets.roll;
			if (common.property_flags & tfxEmitterPropertyFlags_relative_angle)
				sprite_transform.rotation += common.transform.world_rotations.roll;
		}

		//----Motion randomness
		data.noise_offset = random_generation.Range(spawn_values.noise_offset_variation) + spawn_values.noise_offset;
		data.noise_resolution = spawn_values.noise_resolution + 0.01f;

		//----Handle
		/*if (common.property_flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
			data.handle = tfxVec2(0.5f, 0.5f);
		}
		else {
			data.handle = GetProperties().image_handle;
		}*/

		//----Image
		//data.image = GetProperties().image;
		if (common.property_flags & tfxEmitterPropertyFlags_random_start_frame && properties.image->animation_frames > 1) {
			data.image_frame = random_generation.Range(properties.image->animation_frames);
		}
		else {
			data.image_frame = properties.start_frame;
		}

		//----Color
		data.color.a = unsigned char(255.f * common.library->emitter_attributes[library_link->emitter_attributes].overtime.blendfactor.GetFirstValue());
		data.intensity = common.library->emitter_attributes[library_link->emitter_attributes].overtime.intensity.GetFirstValue() * current.intensity;
		if (common.property_flags & tfxEmitterPropertyFlags_random_color) {
			float age = random_generation.Range(data.max_age);
			data.color.r = unsigned char(255.f * lookup_overtime_callback(common.library->emitter_attributes[library_link->emitter_attributes].overtime.red, age, data.max_age));
			data.color.g = unsigned char(255.f * lookup_overtime_callback(common.library->emitter_attributes[library_link->emitter_attributes].overtime.green, age, data.max_age));
			data.color.b = unsigned char(255.f * lookup_overtime_callback(common.library->emitter_attributes[library_link->emitter_attributes].overtime.blue, age, data.max_age));
		}
		else {
			data.color.r = unsigned char(255.f * common.library->emitter_attributes[library_link->emitter_attributes].overtime.red.GetFirstValue());
			data.color.g = unsigned char(255.f * common.library->emitter_attributes[library_link->emitter_attributes].overtime.green.GetFirstValue());
			data.color.b = unsigned char(255.f * common.library->emitter_attributes[library_link->emitter_attributes].overtime.blue.GetFirstValue());
		}

	}

	tfxSpawnPosition InitialisePosition3d(tfxEmitterState &current, tfxCommon &common, tfxEmitterSpawnControls &spawn_values, tfxEffectEmitter *library_link, float tween) {
		tfxPROFILE;
		//----Position
		tfxSpawnPosition out;
		tfxVec3 lerp_position = InterpolateVec3(tween, common.transform.captured_position, common.transform.world_position);
		tfxEmitterProperties &properties = library_link->GetProperties();
		if (properties.emission_type == tfxEmissionType::tfxPoint) {
			if (common.property_flags & tfxEmitterPropertyFlags_relative_position)
				out.local_position = 0;
			else {
				if (common.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) {
					out.local_position = lerp_position;
				}
				else {
					tfxVec3 rotvec = mmTransformVector3(common.transform.matrix, -common.handle);
					out.local_position = rotvec + lerp_position;
				}
			}
		}
		else if (properties.emission_type == tfxEmissionType::tfxArea) {
			tfxVec3 position;

			if (common.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				if (common.property_flags & tfxEmitterPropertyFlags_fill_area) {
					if (!(common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
						current.grid_coords.x--;
						if (current.grid_coords.x < 0.f) {
							current.grid_coords.y--;
							current.grid_coords.x = properties.grid_points.x - 1;
							if (current.grid_coords.y < 0.f) {
								current.grid_coords.z--;
								current.grid_coords.y = properties.grid_points.y - 1;
								if (current.grid_coords.z < 0.f)
									current.grid_coords.z = properties.grid_points.z;
							}
						}
					}

					out.local_position = position + (current.grid_coords * spawn_values.grid_segment_size);

					if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
						current.grid_coords.x++;
						if (current.grid_coords.x == properties.grid_points.x) {
							current.grid_coords.y++;
							current.grid_coords.x = 0.f;
							if (current.grid_coords.y >= properties.grid_points.y) {
								current.grid_coords.z++;
								current.grid_coords.y = 0.f;
								if (current.grid_coords.z >= properties.grid_points.z)
									current.grid_coords.z = 0.f;
							}
						}
					}
				}
				else {
					if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
						if (current.grid_direction.z == 0) {
							//left side
							current.grid_coords.z--;
							current.grid_coords.x = 0.f;
							if (current.grid_coords.z < 0.f) {
								current.grid_coords.y--;
								current.grid_coords.z = properties.grid_points.z - 1;
								if (current.grid_coords.y < 0.f) {
									current.grid_coords.y = properties.grid_points.y - 1;
									current.grid_direction.z++;
								}
							}
						}
						else if (current.grid_direction.z == 1) {
							//right side
							current.grid_coords.z--;
							current.grid_coords.x = properties.grid_points.x - 1;
							if (current.grid_coords.z < 0.f) {
								current.grid_coords.y--;
								current.grid_coords.z = properties.grid_points.z - 1;
								if (current.grid_coords.y < 0.f) {
									current.grid_coords.y = properties.grid_points.y - 1;
									current.grid_direction.z++;
								}
							}
						}
						else if (current.grid_direction.z == 2) {
							//top side
							current.grid_coords.z--;
							current.grid_coords.y = 0.f;
							if (current.grid_coords.z < 0.f) {
								current.grid_coords.x--;
								current.grid_coords.z = properties.grid_points.z - 1;
								if (current.grid_coords.x < 0.f) {
									current.grid_coords.x = properties.grid_points.x - 1;
									current.grid_direction.z++;
								}
							}
						}
						else if (current.grid_direction.z == 3) {
							//bottom side
							current.grid_coords.z--;
							current.grid_coords.y = properties.grid_points.y - 1;
							if (current.grid_coords.z < 0.f) {
								current.grid_coords.x--;
								current.grid_coords.z = properties.grid_points.z - 1;
								if (current.grid_coords.x < 0.f) {
									current.grid_coords.x = properties.grid_points.x - 1;
									current.grid_direction.z++;
								}
							}
						}
						else if (current.grid_direction.z == 4) {
							//End far
							current.grid_coords.x--;
							current.grid_coords.z = 0.f;
							if (current.grid_coords.x < 0.f) {
								current.grid_coords.y--;
								current.grid_coords.x = properties.grid_points.x - 1;
								if (current.grid_coords.y < 0.f) {
									current.grid_coords.y = properties.grid_points.y - 1;
									current.grid_direction.z++;
								}
							}
						}
						else if (current.grid_direction.z == 5) {
							//End near
							current.grid_coords.x--;
							current.grid_coords.z = properties.grid_points.z - 1;
							if (current.grid_coords.x < 0.f) {
								current.grid_coords.y--;
								current.grid_coords.x = properties.grid_points.x - 1;
								if (current.grid_coords.y < 0.f) {
									current.grid_coords.y = properties.grid_points.y - 1;
									current.grid_direction.z = 0;
								}
							}
						}
					}
					else {
						if (current.grid_direction.z == 0) {
							//left side
							current.grid_coords.z++;
							current.grid_coords.x = 0.f;
							if (current.grid_coords.z >= properties.grid_points.z) {
								current.grid_coords.y++;
								current.grid_coords.z = 0.f;
								if (current.grid_coords.y >= properties.grid_points.y) {
									current.grid_coords.y = 0.f;
									current.grid_direction.z++;
								}
							}
						}
						else if (current.grid_direction.z == 1) {
							//right side
							current.grid_coords.z++;
							current.grid_coords.x = properties.grid_points.x - 1;
							if (current.grid_coords.z >= properties.grid_points.z) {
								current.grid_coords.y++;
								current.grid_coords.z = 0.f;
								if (current.grid_coords.y >= properties.grid_points.y) {
									current.grid_coords.y = 0.f;
									current.grid_coords.x = 0.f;
									current.grid_direction.z++;
								}
							}
						}
						else if (current.grid_direction.z == 2) {
							//top side
							current.grid_coords.z++;
							current.grid_coords.y = 0.f;
							if (current.grid_coords.z >= properties.grid_points.z) {
								current.grid_coords.x++;
								current.grid_coords.z = 0.f;
								if (current.grid_coords.x >= properties.grid_points.x) {
									current.grid_coords.x = 0.f;
									current.grid_direction.z++;
								}
							}
						}
						else if (current.grid_direction.z == 3) {
							//bottom side
							current.grid_coords.z++;
							current.grid_coords.y = properties.grid_points.y - 1;
							if (current.grid_coords.z >= properties.grid_points.z) {
								current.grid_coords.x++;
								current.grid_coords.z = 0.f;
								if (current.grid_coords.x >= properties.grid_points.x) {
									current.grid_coords.x = 0.f;
									current.grid_coords.y = 0.f;
									current.grid_direction.z++;
								}
							}
						}
						else if (current.grid_direction.z == 4) {
							//End far
							current.grid_coords.x++;
							current.grid_coords.z = 0.f;
							if (current.grid_coords.x >= properties.grid_points.x) {
								current.grid_coords.y++;
								current.grid_coords.x = 0.f;
								if (current.grid_coords.y >= properties.grid_points.y) {
									current.grid_coords.y = 0.f;
									current.grid_direction.z++;
								}
							}
						}
						else if (current.grid_direction.z == 5) {
							//End near
							current.grid_coords.x++;
							current.grid_coords.z = properties.grid_points.z - 1;
							if (current.grid_coords.x >= properties.grid_points.x) {
								current.grid_coords.y++;
								current.grid_coords.x = 0.f;
								if (current.grid_coords.y >= properties.grid_points.y) {
									current.grid_coords.y = 0.f;
									current.grid_coords.z = 0.f;
									current.grid_direction.z = 0.f;
								}
							}
						}
					}

					tfxBound3d(current.grid_coords, properties.grid_points);
					out.local_position = position + (current.grid_coords * spawn_values.grid_segment_size);
				}
			}
			else {
				if (common.property_flags & tfxEmitterPropertyFlags_fill_area) {
					position.x = random_generation.Range(current.emitter_size.x);
					position.y = random_generation.Range(current.emitter_size.y);
					position.z = random_generation.Range(current.emitter_size.z);
				}
				else {
					//Spawn on one of 6 edges of the cuboid
					tfxU32 side = random_generation.RangeUInt(6);
					if (side == 0) {
						//left side
						position.x = 0.f;
						position.y = random_generation.Range(current.emitter_size.y);
						position.z = random_generation.Range(current.emitter_size.z);
					}
					else if (side == 1) {
						//right side
						position.x = current.emitter_size.x;
						position.y = random_generation.Range(current.emitter_size.y);
						position.z = random_generation.Range(current.emitter_size.z);
					}
					else if (side == 2) {
						//top side
						position.x = random_generation.Range(current.emitter_size.x);
						position.y = 0.f;
						position.z = random_generation.Range(current.emitter_size.z);
					}
					else if (side == 3) {
						//bottom side
						position.x = random_generation.Range(current.emitter_size.x);
						position.y = current.emitter_size.y;
						position.z = random_generation.Range(current.emitter_size.z);
					}
					else if (side == 4) {
						//End far
						position.x = random_generation.Range(current.emitter_size.x);
						position.y = random_generation.Range(current.emitter_size.y);
						position.z = current.emitter_size.z;
					}
					else if (side == 5) {
						//End near
						position.x = random_generation.Range(current.emitter_size.x);
						position.y = random_generation.Range(current.emitter_size.y);
						position.z = 0.f;
					}
				}

				out.local_position = position;
			}

			//----TForm and Emission
			if (!(common.property_flags & tfxEmitterPropertyFlags_relative_position)) {
				out.local_position = mmTransformVector3(common.transform.matrix, out.local_position + common.handle);
				out.local_position = lerp_position + out.local_position * common.transform.scale;
			}

		}
		else if (properties.emission_type == tfxEmissionType::tfxEllipse) {
			tfxVec3 emitter_size = current.emitter_size * .5f;
			tfxVec3 position;

			if (!(common.property_flags & tfxEmitterPropertyFlags_fill_area)) {
				float u = random_generation.Range(1.f);
				float v = random_generation.Range(1.f);
				float theta = u * 2.f * tfxPI;
				float phi = std::acosf(2.f * v - 1.f);
				float sin_theta = std::sinf(theta);
				float cos_theta = std::cosf(theta);
				float sin_phi = std::sinf(phi);
				float cos_phi = std::cosf(phi);
				out.local_position.x = emitter_size.x * sin_phi * cos_theta;
				out.local_position.y = emitter_size.y * sin_phi * sin_theta;
				out.local_position.z = emitter_size.z * cos_phi;
			}
			else {
				position.x = random_generation.Range(-emitter_size.x, emitter_size.x);
				position.y = random_generation.Range(-emitter_size.y, emitter_size.y);
				position.z = random_generation.Range(-emitter_size.z, emitter_size.z);

				while (std::powf(position.x / emitter_size.x, 2.f) + std::powf(position.y / emitter_size.y, 2.f) + std::powf(position.z / emitter_size.z, 2.f) > 1.f) {
					position.x = random_generation.Range(-emitter_size.x, emitter_size.x);
					position.y = random_generation.Range(-emitter_size.y, emitter_size.y);
					position.z = random_generation.Range(-emitter_size.z, emitter_size.z);
				}

				out.local_position = position;
			}

			//----TForm and Emission
			if (!(common.property_flags & tfxEmitterPropertyFlags_relative_position)) {
				out.local_position = mmTransformVector3(common.transform.matrix, out.local_position + common.handle);
				out.local_position = lerp_position + out.local_position * common.transform.scale;
			}
		}
		else if (properties.emission_type == tfxEmissionType::tfxLine) {
			if (common.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				current.grid_coords.x = 0.f;
				current.grid_coords.z = 0.f;

				if (!(common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					current.grid_coords.y--;
					if (current.grid_coords.y < 0.f) {
						current.grid_coords.y = properties.grid_points.x - 1;
					}
				}

				out.local_position = current.grid_coords * spawn_values.grid_segment_size;

				if (common.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					current.grid_coords.y++;
					if (current.grid_coords.y >= properties.grid_points.x) {
						current.grid_coords.y = 0.f;
					}
				}

			}
			else {
				out.local_position.x = 0.f;
				out.local_position.y = random_generation.Range(0.f, current.emitter_size.y);
				out.local_position.z = 0.f;
			}

			//----TForm and Emission
			if (!(common.property_flags & tfxEmitterPropertyFlags_relative_position) && !(common.property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
				out.local_position = mmTransformVector3(common.transform.matrix, out.local_position + common.handle);
				out.local_position = lerp_position + out.local_position * common.transform.scale;
			}
		}

		//----Splatter
		if (spawn_values.splatter) {
			float splatx = random_generation.Range(-spawn_values.splatter, spawn_values.splatter);
			float splaty = random_generation.Range(-spawn_values.splatter, spawn_values.splatter);
			float splatz = random_generation.Range(-spawn_values.splatter, spawn_values.splatter);

			while (std::powf(splatx / spawn_values.splatter, 2.f) + std::powf(splaty / spawn_values.splatter, 2.f) + std::powf(splatz / spawn_values.splatter, 2.f) > 1.f) {
				splatx = random_generation.Range(-spawn_values.splatter, spawn_values.splatter);
				splaty = random_generation.Range(-spawn_values.splatter, spawn_values.splatter);
				splatz = random_generation.Range(-spawn_values.splatter, spawn_values.splatter);
			}

			if (!(common.property_flags & tfxEmitterPropertyFlags_relative_position)) {
				out.local_position.x += splatx * common.transform.scale.x;
				out.local_position.y += splaty * common.transform.scale.y;
				out.local_position.z += splatz * common.transform.scale.z;
			}
			else {
				out.local_position.x += splatx;
				out.local_position.y += splaty;
				out.local_position.z += splatz;
			}
		}

		bool line = common.property_flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type == tfxEmissionType::tfxLine;

		if (!line && !(common.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			if (!(common.property_flags & tfxEmitterPropertyFlags_relative_position) && !(common.property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
				out.world_position = out.local_position;
			}
			else {
				tfxVec4 rotatevec = mmTransformVector(common.transform.matrix, out.local_position + common.handle);
				out.world_position = common.transform.world_position + rotatevec.xyz() * common.transform.scale;
			}
			out.captured_position = out.world_position;
		}

		//----Weight
		if (spawn_values.weight) {
			out.base_weight = spawn_values.weight;
			if (spawn_values.weight_variation > 0) {
				out.base_weight += random_generation.Range(-spawn_values.weight_variation, spawn_values.weight_variation);
			}
		}
		else {
			out.base_weight = 0;
		}
		out.weight_acceleration = out.base_weight * common.library->emitter_attributes[library_link->emitter_attributes].overtime.weight.GetFirstValue() * tfxUPDATE_TIME;

		//----Velocity
		float emission_pitch = lookup_callback(common.library->emitter_attributes[library_link->emitter_attributes].properties.emission_pitch, common.frame);
		float emission_yaw = lookup_callback(common.library->emitter_attributes[library_link->emitter_attributes].properties.emission_yaw, common.frame);

		if (!(common.property_flags & tfxEmitterPropertyFlags_edge_traversal) || properties.emission_type != tfxEmissionType::tfxLine) {
			out.velocity_normal = GetEmissionDirection3d(common, current, library_link, emission_pitch, emission_yaw, out.local_position, out.world_position, current.emitter_size);
		}
		else if (common.property_flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type == tfxEmissionType::tfxLine) {
			out.velocity_normal = tfxVec3(0, 1.f, 0.f);
		}
		out.base_velocity = spawn_values.velocity + random_generation.Range(-spawn_values.velocity_variation, spawn_values.velocity_variation);
		float velocity_scale = common.library->emitter_attributes[library_link->emitter_attributes].overtime.velocity.GetFirstValue() * current.velocity_adjuster * out.base_velocity;

		//data.velocity = data.velocity_normal * data.base.velocity * data.velocity_scale * tfxUPDATE_TIME;

		//Do a micro update
		//A bit hacky but the epsilon after tween just ensures that theres a guaranteed small difference between captured/world positions so that
		//the alignment on the first frame can be calculated
		float micro_time = tfxUPDATE_TIME * (1 - tween + 0.001f);
		out.weight_acceleration += out.base_weight * common.library->emitter_attributes[library_link->emitter_attributes].overtime.weight.GetFirstValue() * micro_time;
		//----Velocity Changes
		tfxVec3 current_velocity = out.velocity_normal.xyz() * (out.base_velocity * common.library->emitter_attributes[library_link->emitter_attributes].overtime.velocity.GetFirstValue());
		current_velocity.y -= out.weight_acceleration;
		if (properties.vector_align_type == tfxVectorAlignType_motion) {
			float l = FastLength(current_velocity * tfxUPDATE_TIME);
			out.velocity_normal.w = common.library->emitter_attributes[library_link->emitter_attributes].overtime.stretch.GetFirstValue() * l * 10.f;
		}
		current_velocity *= micro_time;
		out.local_position += current_velocity;
		if (line || common.property_flags & tfxEmitterPropertyFlags_relative_position) {
			if (!(common.property_flags & tfxEmitterPropertyFlags_relative_position) && !(common.property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
				out.world_position = out.local_position;
			}
			else {
				tfxVec4 rotatevec = mmTransformVector(common.transform.matrix, out.local_position + common.handle);
				out.captured_position = common.transform.captured_position + rotatevec.xyz() * common.transform.scale;
				out.world_position = common.transform.world_position + rotatevec.xyz() * common.transform.scale;
			}
		}
		else {
			out.world_position += current_velocity;
		}

		return out;
	}

	void InitialiseParticle3d(tfxParticleData &data, tfxSpriteTransform3d &sprite_transform, tfxEmitterState &current, tfxCommon &common, tfxEmitterSpawnControls &spawn_values, tfxEffectEmitter *library_link, float tween) {
		tfxPROFILE;
		//----Spin
		data.base.spin = random_generation.Range(-spawn_values.spin_variation, std::abs(spawn_values.spin_variation)) + spawn_values.spin;

		tfxEmitterProperties &properties = library_link->GetProperties();

		sprite_transform.rotations.roll = data.local_rotations.pitch = 0;
		sprite_transform.rotations.roll = data.local_rotations.yaw = 0;
		sprite_transform.rotations.roll = data.local_rotations.roll = 0;
		if (properties.angle_settings & tfxAngleSettingFlags_specify_roll)
			sprite_transform.rotations.roll = data.local_rotations.roll = properties.angle_offsets.roll;
		if (properties.angle_settings & tfxAngleSettingFlags_specify_pitch)
			sprite_transform.rotations.pitch = data.local_rotations.pitch = properties.angle_offsets.pitch;
		if (properties.angle_settings & tfxAngleSettingFlags_specify_yaw)
			sprite_transform.rotations.yaw = data.local_rotations.yaw = properties.angle_offsets.yaw;
		if (properties.angle_settings & tfxAngleSettingFlags_random_pitch)
			sprite_transform.rotations.pitch = data.local_rotations.pitch = random_generation.Range(properties.angle_offsets.pitch);
		if (properties.angle_settings & tfxAngleSettingFlags_random_yaw)
			sprite_transform.rotations.yaw = data.local_rotations.yaw = random_generation.Range(properties.angle_offsets.yaw);
		if (properties.angle_settings & tfxAngleSettingFlags_random_roll)
			sprite_transform.rotations.roll = data.local_rotations.roll = random_generation.Range(properties.angle_offsets.roll);

		if (common.property_flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type == tfxLine) {
			sprite_transform.rotations = data.local_rotations;
			float s = sin(data.local_rotations.roll);
			float c = cos(data.local_rotations.roll);
			tfxMatrix2 pmat;
			pmat.Set(c, s, -s, c);
			pmat = pmat.Transform(common.transform.matrix);
		}
		else if (common.property_flags & tfxEmitterPropertyFlags_relative_position) {
			sprite_transform.rotations = data.local_rotations;
			float s = sin(data.local_rotations.roll);
			float c = cos(data.local_rotations.roll);
			tfxMatrix2 pmat;
			pmat.Set(c, s, -s, c);
			pmat = pmat.Transform(common.transform.matrix);
		}
		else if (common.property_flags & tfxEmitterPropertyFlags_relative_angle) {
			sprite_transform.rotations = common.transform.world_rotations + data.local_rotations;
		}

		//----Size
		if (!(common.property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			float random_size_x = random_generation.Range(spawn_values.size_variation.x);
			float random_size_y = random_generation.Range(spawn_values.size_variation.y);
			data.base.size.y = random_size_y + spawn_values.size.y;
			data.base.size.x = random_size_x + spawn_values.size.x;
			float height = data.base.size.y;

			sprite_transform.scale.x = data.base.size.x * common.library->emitter_attributes[library_link->emitter_attributes].overtime.width.GetFirstValue();

			if (common.property_flags & tfxEmitterPropertyFlags_lifetime_uniform_size) {
				sprite_transform.scale.y = height * common.library->emitter_attributes[library_link->emitter_attributes].overtime.width.GetFirstValue();
			}
			else {
				sprite_transform.scale.y = height * common.library->emitter_attributes[library_link->emitter_attributes].overtime.height.GetFirstValue();
			}
		}
		else {
			float random_size_x = random_generation.Range(spawn_values.size_variation.x);
			float random_size_y = random_size_x;
			data.base.size.y = random_size_y + spawn_values.size.y;
			data.base.size.x = random_size_x + spawn_values.size.x;
			float height = data.base.size.y;

			sprite_transform.scale.x = data.base.size.x * common.library->emitter_attributes[library_link->emitter_attributes].overtime.width.GetFirstValue();

			if (common.property_flags & tfxEmitterPropertyFlags_lifetime_uniform_size) {
				sprite_transform.scale.y = sprite_transform.scale.x;
			}
			else {
				sprite_transform.scale.y = height * common.library->emitter_attributes[library_link->emitter_attributes].overtime.height.GetFirstValue();
			}
		}

		sprite_transform.scale *= current.overal_scale;

		//----Motion randomness
		data.noise_offset = random_generation.Range(spawn_values.noise_offset_variation) + spawn_values.noise_offset;
		data.noise_resolution = spawn_values.noise_resolution + 0.01f;

		//end micro update

		//----Handle
		/*if (common.property_flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
			data.handle = tfxVec2(0.5f, 0.5f);
		}
		else {
			data.handle = GetProperties().image_handle;
		}*/

		//----Image
		//data.image = GetProperties().image;
		if (common.property_flags & tfxEmitterPropertyFlags_random_start_frame && properties.image->animation_frames > 1) {
			data.image_frame = random_generation.Range(properties.image->animation_frames);
		}
		else {
			data.image_frame = properties.start_frame;
		}

		//----Color
		data.color.a = unsigned char(255.f * common.library->emitter_attributes[library_link->emitter_attributes].overtime.blendfactor.GetFirstValue());
		data.intensity = common.library->emitter_attributes[library_link->emitter_attributes].overtime.intensity.GetFirstValue() * current.intensity;
		if (common.property_flags & tfxEmitterPropertyFlags_random_color) {
			float age = random_generation.Range(data.max_age);
			data.color.r = unsigned char(255.f * lookup_overtime_callback(common.library->emitter_attributes[library_link->emitter_attributes].overtime.red, age, data.max_age));
			data.color.g = unsigned char(255.f * lookup_overtime_callback(common.library->emitter_attributes[library_link->emitter_attributes].overtime.green, age, data.max_age));
			data.color.b = unsigned char(255.f * lookup_overtime_callback(common.library->emitter_attributes[library_link->emitter_attributes].overtime.blue, age, data.max_age));
		}
		else {
			data.color.r = unsigned char(255.f * common.library->emitter_attributes[library_link->emitter_attributes].overtime.red.GetFirstValue());
			data.color.g = unsigned char(255.f * common.library->emitter_attributes[library_link->emitter_attributes].overtime.green.GetFirstValue());
			data.color.b = unsigned char(255.f * common.library->emitter_attributes[library_link->emitter_attributes].overtime.blue.GetFirstValue());
		}
	}

	void ReloadBaseValues(tfxParticle &p, tfxEffectEmitter &e) {
		//----Life
		//std::uniform_real_distribution<float> random_life(0, e.current.life_variation);
		//p.data.max_age = e.current.life + random_life(random_generation.engine);

		//----Velocity
		//std::uniform_real_distribution<float> random_velocity;
		//if (e.current.velocity_variation > 0)
			//random_velocity = std::uniform_real_distribution<float>(0, e.current.velocity_variation);
		//else
			//random_velocity = std::uniform_real_distribution<float>(e.current.velocity_variation, 0);
		//p.data.velocity_scale = e.common.library->emitter_attributes[e.emitter_attributes].overtime.velocity.GetFirstValue() * e.current.velocity_adjuster;
		//p.data.base.velocity = e.current.velocity + random_velocity(random_generation.engine);
		float velocity_scale = lookup_overtime_callback(e.common.library->emitter_attributes[e.emitter_attributes].overtime.velocity, p.data.age, p.data.max_age) * e.current.velocity_adjuster;

		//----Size
		if (!(e.common.property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			//std::uniform_real_distribution<float> random_width(0, e.current.size_variation.x);
			//std::uniform_real_distribution<float> random_height(0, e.current.size_variation.y);

			//p.data.base.size.y = p.data.base.random_size.y + e.current..y;
			//p.data.base.size.x = (p.data.base.random_size.x + e.current.size.x) / e.GetProperties().image->image_size.x;
			//p.data.base.size.y = p.data.base.height / e.GetProperties().image->image_size.y;

			//p.data.local.scale.x = p.data.base.size.x * e.common.library->emitter_attributes[e.emitter_attributes].overtime.width.GetFirstValue();

			if (e.common.library->emitter_attributes[e.emitter_attributes].overtime.stretch.GetFirstValue()) {
				float velocity = std::fabsf(velocity_scale * p.data.base.velocity);
				//p.data.local.scale.y = (e.common.library->emitter_attributes[e.emitter_attributes].overtime.height.GetFirstValue() * e.parent->current.size.y * (p.data.base.height + (velocity * e.common.library->emitter_attributes[e.emitter_attributes].overtime.stretch.GetFirstValue() * e.parent->current.stretch))) / e.GetProperties().image.image_size.y;
			}
			else {
				if (e.common.property_flags & tfxEmitterPropertyFlags_lifetime_uniform_size) {
					//p.data.local.scale.y = p.data.base.size.y * e.common.library->emitter_attributes[e.emitter_attributes].overtime.width.GetFirstValue();
				}
				else {
					//p.data.local.scale.y = p.data.base.size.y * e.common.library->emitter_attributes[e.emitter_attributes].overtime.height.GetFirstValue();
				}
			}
		}
		else {
			//std::uniform_real_distribution<float> random_width(0, e.current.size_variation.x);

			//p.data.base.size.y = p.data.base.random_size.y + e.current.size.x;
			//p.data.base.size.x = (p.data.base.random_size.x + e.current.size.x) / e.GetProperties().image->image_size.x;
			//p.data.base.size.y = p.data.base.height / e.GetProperties().image->image_size.y;

			//p.data.local.scale.x = p.data.base.size.x * e.common.library->emitter_attributes[e.emitter_attributes].overtime.width.GetFirstValue();

			if (e.common.library->emitter_attributes[e.emitter_attributes].overtime.stretch.GetFirstValue()) {
				float velocity = std::fabsf(velocity_scale * p.data.base.velocity);
				//p.data.local.scale.y = (e.common.library->emitter_attributes[e.emitter_attributes].overtime.height.GetFirstValue() * e.parent->current.size.y * (p.data.base.height + (velocity * e.common.library->emitter_attributes[e.emitter_attributes].overtime.stretch.GetFirstValue() * e.parent->current.stretch))) / e.GetProperties().image.image_size.y;
			}
			else {
				//p.data.local.scale.y = p.data.local.scale.x;
			}
		}

		//----Spin
		//p.data.base.spin = random_generation.Range(-e.current.spin_variation, std::abs(e.current.spin_variation)) + e.current.spin;

		//std::uniform_real_distribution<float> random_angle(0, e.GetProperties().angle_offset);
		//switch (e.GetProperties().angle_setting) {
		//case AngleSetting::kRandom:
			//p.data.rotations.roll = random_angle(random_generation.engine);
			//break;
		//case AngleSetting::kSpecify:
			//p.data.rotations.roll = e.GetProperties().angle_offset;
			//break;
		//default:
			//break;
		//}

		//----Weight
		//if (e.current.weight) {
			//if (e.current.weight_variation > 0) {
				//p.data.base.weight = (random_generation.Range(e.current.weight_variation) + e.current.weight) * e.common.library->emitter_attributes[e.emitter_attributes].overtime.weight.GetFirstValue();
			//}
			//else {
				//p.data.base.weight = (random_generation.Range(e.current.weight_variation, 0) + e.current.weight) * e.common.library->emitter_attributes[e.emitter_attributes].overtime.weight.GetFirstValue();
			//}
			//p.data.weight_acceleration = 0;
		//}
		//else {
			//p.data.weight_acceleration = 0;
			//p.data.base.weight = 0;
		//}

		//TransformParticle(p, e);

		//----Motion randomness
		//p.data.noise_offset = e.current.noise_offset;
		//float mr = tfxRadians(p.data.direction_turbulance * e.common.library->emitter_attributes[e.emitter_attributes].overtime.direction_turbulance.GetFirstValue());
		//std::uniform_real_distribution<float> random_motion(-mr, mr);
		//p.data.motion_randomness_direction = tfxRadians(22.5f * random_motion(random_generation.engine));
		//p.data.motion_randomness_speed = 30.f * random_motion(random_generation.engine);
		//p.data.motion_tracker = 0;

		//if (!e.GetProperties().edge_traversal || e.GetProperties().emission_type != tfxEmissionType::kLine) {
			//p.data.direction = p.data.emission_angle = e.GetEmissionDirection2d(p) + p.data.motion_randomness_direction;
		//}

		//----Normalize Velocity to direction
		//p.data.velocity_normal.x = std::sinf(p.data.direction);
		//p.data.velocity_normal.y = -std::cosf(p.data.direction);

		//p.data.velocity = p.data.velocity_normal * p.data.base.velocity * p.data.velocity_scale * e.timer->UpdateTime();
		//bool line = e.GetProperties().edge_traversal && e.GetProperties().emission_type == tfxEmissionType::kLine;

		//if (e.GetProperties().angle_setting == AngleSetting::kAlign && !line) {
			//p.data.rotations.roll = p.data.rotations.roll = GetVectorAngle(p.data.velocity_normal.x, p.data.velocity_normal.y) + e.GetProperties().angle_offset + e.world.position.w;
			//p.data.captured_position.w = p.data.rotations.roll;
		//}

		//----Handle
		//if (e.common.property_flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
			//p.data.handle = tfxVec2(0.5f, 0.5f);
		//}
		//else {
			//p.data.handle = e.GetProperties().image_handle;
		//}

		//----Image
		//p.data.image = GetProperties().image;
		//if (e.GetProperties().image.random_start_frame && e.GetProperties().image.animation_frames > 1) {
			//std::uniform_real_distribution<float> random_start_frame(0.f, e.GetProperties().image.animation_frames - 1);
			//p.data.image_frame = random_start_frame(random_generation.engine);
		//}
		//else {
			//p.data.image_frame = e.GetProperties().image.start_frame;
		//}

		//----Color
		//p.data.color.a = unsigned char(255.f * e.common.library->emitter_attributes[e.emitter_attributes].overtime.opacity.GetFirstValue());
		//p.data.intensity = e.common.library->emitter_attributes[e.emitter_attributes].overtime.opacity.GetFirstValue();
		//if (e.GetProperties().random_color) {
			//std::uniform_real_distribution<float> random_color(0.f, p.data.max_age);
			//float age = random_color(random_generation.engine);
			//p.data.color.r = unsigned char(255.f * lookup_overtime_callback(e.common.library->emitter_attributes[e.emitter_attributes].overtime.red, age, p.data.max_age));
			//p.data.color.g = unsigned char(255.f * lookup_overtime_callback(e.common.library->emitter_attributes[e.emitter_attributes].overtime.green, age, p.data.max_age));
			//p.data.color.b = unsigned char(255.f * lookup_overtime_callback(e.common.library->emitter_attributes[e.emitter_attributes].overtime.blue, age, p.max_age));
		//}
		//else {
			//p.data.color.r = unsigned char(255.f * e.common.library->emitter_attributes[e.emitter_attributes].overtime.red.GetFirstValue());
			//p.data.color.g = unsigned char(255.f * e.common.library->emitter_attributes[e.emitter_attributes].overtime.green.GetFirstValue());
			//p.data.color.b = unsigned char(255.f * e.common.library->emitter_attributes[e.emitter_attributes].overtime.blue.GetFirstValue());
		//}

	}

	void UpdateParticle2d(tfxParticleData &data, tfxVec2 &sprite_scale, tfxControlData &c) {
		tfxPROFILE;
		tfxU32 lookup_frame = static_cast<tfxU32>((data.age / data.max_age * c.graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

		float lookup_velocity = c.graphs->velocity.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->velocity.lookup.last_frame)] * c.velocity_adjuster;
		float lookup_velocity_turbulance = c.graphs->velocity_turbulance.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->velocity_turbulance.lookup.last_frame)];
		float lookup_direction = c.graphs->direction.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->direction.lookup.last_frame)] + data.velocity_normal.x;
		float lookup_noise_resolution = c.graphs->noise_resolution.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->noise_resolution.lookup.last_frame)] * data.noise_resolution;
		float lookup_stretch = c.graphs->stretch.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->stretch.lookup.last_frame)];
		float lookup_weight = c.graphs->weight.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->weight.lookup.last_frame)];

		float lookup_width = c.graphs->width.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->width.lookup.last_frame)];
		float lookup_height = c.graphs->height.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->height.lookup.last_frame)];
		float lookup_spin = c.graphs->spin.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->spin.lookup.last_frame)] * data.base.spin;
		float lookup_red = c.graphs->red.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->red.lookup.last_frame)];
		float lookup_green = c.graphs->green.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->green.lookup.last_frame)];
		float lookup_blue = c.graphs->blue.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->blue.lookup.last_frame)];
		float lookup_opacity = c.graphs->blendfactor.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->blendfactor.lookup.last_frame)];
		float lookup_intensity = c.graphs->intensity.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->intensity.lookup.last_frame)];

		float direction = 0;

		tfxVec2 mr_vec;
		if (c.flags & tfxEmitterStateFlags_not_line) {
			direction = lookup_direction;
		}

		if (lookup_velocity_turbulance) {
			float eps = 0.0001f;
			float eps2 = 0.0001f * 2.f;

			float x = data.local_position.x / lookup_noise_resolution + data.noise_offset;
			float y = data.local_position.y / lookup_noise_resolution + data.noise_offset;

			//Find rate of change in YZ plane
			float n1 = tfxSimplexNoise::noise(x, y + eps);
			float n2 = tfxSimplexNoise::noise(x, y - eps);
			//Average to find approximate derivative
			float a = (n1 - n2) / eps2;
			n1 = tfxSimplexNoise::noise(x, y);
			n2 = tfxSimplexNoise::noise(x, y);
			//Average to find approximate derivative
			float b = (n1 - n2) / eps2;
			mr_vec.x = a - b;

			//Find rate of change in XZ plane
			n1 = tfxSimplexNoise::noise(x, y);
			n2 = tfxSimplexNoise::noise(x, y);
			a = (n1 - n2) / eps2;
			n1 = tfxSimplexNoise::noise(x + eps, y);
			n2 = tfxSimplexNoise::noise(x - eps, y);
			b = (n1 - n2) / eps2;
			mr_vec.y = a - b;
			mr_vec *= lookup_velocity_turbulance;
		}

		//----Weight Changes
		data.weight_acceleration += data.base.weight * lookup_weight * tfxUPDATE_TIME;

		//----Velocity Changes
		tfxVec2 velocity_normal;
		velocity_normal.x = std::sinf(direction);
		velocity_normal.y = -std::cosf(direction);

		tfxVec2 current_velocity = (data.base.velocity * lookup_velocity) * velocity_normal;
		current_velocity += mr_vec;
		current_velocity.y += data.weight_acceleration;
		current_velocity *= tfxUPDATE_TIME;

		//----Color changes
		data.color.a = unsigned char(255.f * lookup_opacity);
		data.intensity = lookup_intensity * c.global_intensity;
		if (!(c.flags & tfxEmitterStateFlags_random_color)) {
			data.color.r = unsigned char(255.f * lookup_red);
			data.color.g = unsigned char(255.f * lookup_green);
			data.color.b = unsigned char(255.f * lookup_blue);
		}

		data.color = tfxRGBA8(data.color.r, data.color.g, data.color.b, data.color.a);

		//----Size Changes
		tfxVec2 scale;
		scale.x = data.base.size.x * lookup_width;
		if (scale.x < 0.f)
			scale.x = scale.x;

		//----Stretch Changes
		float velocity = std::fabsf(lookup_velocity * data.base.velocity + data.weight_acceleration);
		if (c.flags & tfxEmitterStateFlags_lifetime_uniform_size) {
			scale.y = (lookup_width * (data.base.size.y + (velocity * lookup_stretch * c.stretch))) / c.image_size_y;
			if (c.flags & tfxEmitterPropertyFlags_base_uniform_size && scale.y < scale.x)
				scale.y = scale.x;
		}
		else
			scale.y = (lookup_height * (data.base.size.y + (velocity * lookup_stretch * c.stretch))) / c.image_size_y;

		//----Spin and angle Changes
		float spin = 0;
		if (c.flags & tfxEmitterStateFlags_can_spin) {
			spin = lookup_spin;
		}

		//---------------
		//Now that the latest changes are applied, affect the particle state
		//---------------

		//Before applying the behaviour to the particle position, scale and rotation, you have the chance to override them here
		if (c.particle_update_callback)
			c.particle_update_callback(data, c.user_data);

		//----Rotation
		if (c.flags & tfxEmitterStateFlags_align_with_velocity) {
			tfxVec2 vd = current_velocity.IsNill() ? velocity_normal : current_velocity;
			data.local_rotations.roll = GetVectorAngle(vd.x, vd.y) + c.angle_offset;
		}
		else {
			data.local_rotations.roll += spin * tfxUPDATE_TIME;
		}

		//----Position
		data.local_position += current_velocity * c.overal_scale;

		//----Scale
		sprite_scale = scale * c.overal_scale;

		//Lines - Reposition if the particle is travelling along a line
		tfxVec2 offset = velocity_normal * c.emitter_size_y;
		float length = std::fabsf(data.local_position.y);
		float emitter_length = c.emitter_size_y;
		bool line_and_kill = (c.flags & tfxEmitterStateFlags_is_line_traversal) && (c.flags & tfxEmitterStateFlags_kill) && length > emitter_length;
		bool line_and_loop = (c.flags & tfxEmitterStateFlags_is_line_traversal) && (c.flags & tfxEmitterStateFlags_loop) && length > emitter_length;
		if (line_and_loop) {
			data.local_position.y -= offset.y;
			data.flags |= tfxParticleFlags_capture_after_transform;
		}
		else if (line_and_kill) {
			data.flags |= tfxParticleFlags_remove;
		}

		//----Image animation
		data.image_frame += c.image_frame_rate * tfxUPDATE_TIME;
		data.image_frame = (c.flags & tfxEmitterStateFlags_play_once) && data.image_frame > c.end_frame ? data.image_frame = c.end_frame : data.image_frame;
		data.image_frame = (c.flags & tfxEmitterStateFlags_play_once) && data.image_frame < 0 ? data.image_frame = 0 : data.image_frame;
		data.image_frame = std::fmodf(data.image_frame, c.end_frame + 1);

	}

	void UpdateParticle3d(tfxParticleData &data, tfxVec2 &sprite_scale, tfxControlData &c) {
		tfxPROFILE;
		tfxU32 lookup_frame = static_cast<tfxU32>((data.age / data.max_age * c.graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

		float lookup_velocity = c.graphs->velocity.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->velocity.lookup.last_frame)] * c.velocity_adjuster;
		float lookup_velocity_turbulance = c.graphs->velocity_turbulance.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->velocity_turbulance.lookup.last_frame)];
		float lookup_noise_resolution = c.graphs->noise_resolution.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->noise_resolution.lookup.last_frame)] * data.noise_resolution * c.overal_scale;
		float lookup_stretch = c.graphs->stretch.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->stretch.lookup.last_frame)];
		float lookup_weight = c.graphs->weight.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->weight.lookup.last_frame)];

		float lookup_width = c.graphs->width.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->width.lookup.last_frame)];
		float lookup_height = c.graphs->height.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->height.lookup.last_frame)];
		float lookup_spin = c.graphs->spin.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->spin.lookup.last_frame)] * data.base.spin;
		float lookup_red = c.graphs->red.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->red.lookup.last_frame)];
		float lookup_green = c.graphs->green.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->green.lookup.last_frame)];
		float lookup_blue = c.graphs->blue.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->blue.lookup.last_frame)];
		float lookup_opacity = c.graphs->blendfactor.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->blendfactor.lookup.last_frame)];
		float lookup_intensity = c.graphs->intensity.lookup.values[std::min<tfxU32>(lookup_frame, c.graphs->intensity.lookup.last_frame)];

		float direction = 0.f;
		float mr_angle = 0.f;

		tfxVec3 mr_vec;

		if (lookup_velocity_turbulance) {
			float eps = 0.001f;
			float eps2 = 0.001f * 2.f;

			float noise_offset = data.noise_offset * c.overal_scale;

			float x = data.local_position.x / lookup_noise_resolution + noise_offset;
			float y = data.local_position.y / lookup_noise_resolution + noise_offset;
			float z = data.local_position.z / lookup_noise_resolution + noise_offset;

			__m128 x4 = _mm_set1_ps(x);
			__m128 y4 = _mm_set1_ps(y);
			__m128 z4 = _mm_set1_ps(z);

			__m128 xeps4 = _mm_set_ps(x - eps, x + eps, x, x);
			__m128 xeps4r = _mm_set_ps(x, x, x - eps, x + eps);
			__m128 yeps4 = _mm_set_ps(y, y, y - eps, y + eps);
			__m128 yeps4r = _mm_set_ps(y - eps, y + eps, y, y);
			__m128 zeps4 = _mm_set_ps(z - eps, z + eps, z, z);
			__m128 zeps4r = _mm_set_ps(z, z, z - eps, z + eps);

			//Find rate of change in YZ plane
			tfxVec4 yz = tfxSimplexNoise::noise4(x4, yeps4, zeps4);
			//Average to find approximate derivative
			float a = (yz.x - yz.y) / eps2;
			//Average to find approximate derivative
			float b = (yz.z - yz.w) / eps2;
			mr_vec.x = a - b;

			y += 100.f;
			yeps4 = _mm_set_ps(y, y, y - eps, y + eps);
			yeps4r = _mm_set_ps(y - eps, y + eps, y, y);
			//Find rate of change in XZ plane
			tfxVec4 xz = tfxSimplexNoise::noise4(xeps4, y4, zeps4r);
			a = (xz.x - xz.y) / eps2;
			b = (xz.z - xz.w) / eps2;
			mr_vec.y = a - b;

			z += 100.f;
			zeps4 = _mm_set_ps(z - eps, z + eps, z, z);
			zeps4r = _mm_set_ps(z, z, z - eps, z + eps);
			//Find rate of change in XY plane
			tfxVec4 xy = tfxSimplexNoise::noise4(xeps4r, yeps4r, z4);
			a = (xy.x - xy.y) / eps2;
			b = (xy.z - xy.w) / eps2;
			mr_vec.z = a - b;

			mr_vec *= lookup_velocity_turbulance;

		}

		//----Weight Changes
		data.weight_acceleration += data.base.weight * lookup_weight * tfxUPDATE_TIME;

		//----Velocity Changes
		float velocity_scalar = data.base.velocity * lookup_velocity;
		tfxVec3 current_velocity = data.velocity_normal.xyz() * velocity_scalar;
		current_velocity += mr_vec;
		current_velocity.y -= data.weight_acceleration;
		data.velocity_normal.w = lookup_stretch * c.stretch;
		current_velocity *= tfxUPDATE_TIME;

		//----Color changes
		data.color.a = unsigned char(255.f * lookup_opacity);
		data.intensity = lookup_intensity * c.global_intensity;
		if (!(c.flags & tfxEmitterStateFlags_random_color)) {
			data.color.r = unsigned char(255.f * lookup_red);
			data.color.g = unsigned char(255.f * lookup_green);
			data.color.b = unsigned char(255.f * lookup_blue);
		}

		data.color = tfxRGBA8(data.color.r, data.color.g, data.color.b, data.color.a);

		//----Size Changes
		tfxVec2 scale;
		scale.x = data.base.size.x * lookup_width;
		if (scale.x < 0.f)
			scale.x = scale.x;

		//----Stretch Changes
		if (c.flags & tfxEmitterStateFlags_lifetime_uniform_size) {
			scale.y = lookup_width * data.base.size.y;
			if (c.flags & tfxEmitterPropertyFlags_base_uniform_size && scale.y < scale.x)
				scale.y = scale.x;
		}
		else
			scale.y = lookup_height * data.base.size.y;

		//----Spin and angle Changes
		float spin = 0;
		if (c.flags & tfxEmitterStateFlags_can_spin) {
			spin = lookup_spin;
		}

		//---------------
		//Now that the latest changes are applied, affect the particle state
		//---------------

		//Before applying the behaviour to the particle position, scale and rotation, you have the chance to override them here
		if (c.particle_update_callback)
			c.particle_update_callback(data, c.user_data);

		//----Rotation
		if (c.flags & tfxEmitterStateFlags_align_with_velocity) {
			tfxVec3 vd = current_velocity.IsNill() ? data.velocity_normal.xyz() : current_velocity;
			data.local_rotations.roll = GetVectorAngle(vd.x, vd.y) + c.angle_offset;
		}
		else {
			data.local_rotations.roll += spin * tfxUPDATE_TIME;
		}

		//----Position
		data.local_position += current_velocity * c.overal_scale;

		//----Scale
		sprite_scale = scale * c.overal_scale;

		//Lines - Reposition if the particle is travelling along a line
		tfxVec3 offset = data.velocity_normal.xyz() * c.emitter_size_y;
		float length = std::fabsf(data.local_position.y);
		float emitter_length = c.emitter_size_y;
		bool line_and_kill = (c.flags & tfxEmitterStateFlags_is_line_traversal) && (c.flags & tfxEmitterStateFlags_kill) && length > emitter_length;
		bool line_and_loop = (c.flags & tfxEmitterStateFlags_is_line_traversal) && (c.flags & tfxEmitterStateFlags_loop) && length > emitter_length;
		if (line_and_loop) {
			data.local_position.y -= offset.y;
			data.flags |= tfxParticleFlags_capture_after_transform;
		}
		else if (line_and_kill) {
			data.flags |= tfxParticleFlags_remove;
		}

		//----Image animation
		data.image_frame += c.image_frame_rate * tfxUPDATE_TIME;
		data.image_frame = (c.flags & tfxEmitterStateFlags_play_once) && data.image_frame > c.end_frame ? data.image_frame = c.end_frame : data.image_frame;
		data.image_frame = (c.flags & tfxEmitterStateFlags_play_once) && data.image_frame < 0 ? data.image_frame = 0 : data.image_frame;
		data.image_frame = std::fmodf(data.image_frame, c.end_frame + 1);

		data.flags = current_velocity.IsNill() ? data.flags &= ~tfxParticleFlags_has_velocity : data.flags |= tfxParticleFlags_has_velocity;
	}

	void tfxEffectEmitter::ResetGlobalGraphs(bool add_node, bool compile) {
		common.library->global_graphs[global].life.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->global_graphs[global].life.type = tfxGlobal_life;
		common.library->global_graphs[global].amount.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->global_graphs[global].amount.type = tfxGlobal_amount;
		common.library->global_graphs[global].velocity.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->global_graphs[global].velocity.type = tfxGlobal_velocity;
		common.library->global_graphs[global].width.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->global_graphs[global].width.type = tfxGlobal_width;
		common.library->global_graphs[global].height.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->global_graphs[global].height.type = tfxGlobal_height;
		common.library->global_graphs[global].weight.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->global_graphs[global].weight.type = tfxGlobal_weight;
		common.library->global_graphs[global].spin.Reset(1.f, tfxGlobalPercentPresetSigned, add_node); common.library->global_graphs[global].spin.type = tfxGlobal_spin;
		common.library->global_graphs[global].stretch.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->global_graphs[global].stretch.type = tfxGlobal_stretch;
		common.library->global_graphs[global].overal_scale.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->global_graphs[global].overal_scale.type = tfxGlobal_overal_scale;
		common.library->global_graphs[global].intensity.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->global_graphs[global].intensity.type = tfxGlobal_intensity;
		common.library->global_graphs[global].frame_rate.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->global_graphs[global].frame_rate.type = tfxGlobal_frame_rate;
		common.library->global_graphs[global].splatter.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->global_graphs[global].splatter.type = tfxGlobal_splatter;
		common.library->global_graphs[global].roll.Reset(0.f, tfxAnglePreset, add_node); common.library->global_graphs[global].roll.type = tfxGlobal_effect_roll;
		common.library->global_graphs[global].pitch.Reset(0.f, tfxAnglePreset, add_node); common.library->global_graphs[global].pitch.type = tfxGlobal_effect_pitch;
		common.library->global_graphs[global].yaw.Reset(0.f, tfxAnglePreset, add_node); common.library->global_graphs[global].yaw.type = tfxGlobal_effect_yaw;
		common.library->global_graphs[global].emitter_width.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->global_graphs[global].emitter_width.type = tfxGlobal_emitter_width;
		common.library->global_graphs[global].emitter_height.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->global_graphs[global].emitter_height.type = tfxGlobal_emitter_height;
		common.library->global_graphs[global].emitter_depth.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->global_graphs[global].emitter_depth.type = tfxGlobal_emitter_depth;
		if(compile)
			common.library->CompileGlobalGraph(global);
	}

	void tfxEffectEmitter::ResetBaseGraphs(bool add_node, bool compile) {
		common.library->emitter_attributes[emitter_attributes].base.life.Reset(1000.f, tfxLifePreset, add_node); common.library->emitter_attributes[emitter_attributes].base.life.type = tfxBase_life;
		common.library->emitter_attributes[emitter_attributes].base.amount.Reset(1.f, tfxAmountPreset, add_node); common.library->emitter_attributes[emitter_attributes].base.amount.type = tfxBase_amount;
		common.library->emitter_attributes[emitter_attributes].base.velocity.Reset(0.f, tfxVelocityPreset, add_node); common.library->emitter_attributes[emitter_attributes].base.velocity.type = tfxBase_velocity;
		common.library->emitter_attributes[emitter_attributes].base.width.Reset(128.f, tfxDimensionsPreset, add_node); common.library->emitter_attributes[emitter_attributes].base.width.type = tfxBase_width;
		common.library->emitter_attributes[emitter_attributes].base.height.Reset(128.f, tfxDimensionsPreset, add_node); common.library->emitter_attributes[emitter_attributes].base.height.type = tfxBase_height;
		common.library->emitter_attributes[emitter_attributes].base.weight.Reset(0.f, tfxWeightPreset, add_node); common.library->emitter_attributes[emitter_attributes].base.weight.type = tfxBase_weight;
		common.library->emitter_attributes[emitter_attributes].base.spin.Reset(0.f, tfxSpinPreset, add_node); common.library->emitter_attributes[emitter_attributes].base.spin.type = tfxBase_spin;
		common.library->emitter_attributes[emitter_attributes].base.noise_offset.Reset(0.f, tfxGlobalPercentPreset, add_node); common.library->emitter_attributes[emitter_attributes].base.noise_offset.type = tfxBase_noise_offset;
		if(compile)
			common.library->CompileBaseGraph(emitter_attributes);
	}

	void tfxEffectEmitter::ResetPropertyGraphs(bool add_node, bool compile) {
		common.library->emitter_attributes[emitter_attributes].properties.emission_pitch.Reset(0.f, tfxAnglePreset, add_node); common.library->emitter_attributes[emitter_attributes].properties.emission_pitch.type = tfxProperty_emission_pitch;
		common.library->emitter_attributes[emitter_attributes].properties.emission_yaw.Reset(0.f, tfxAnglePreset, add_node); common.library->emitter_attributes[emitter_attributes].properties.emission_yaw.type = tfxProperty_emission_yaw;
		common.library->emitter_attributes[emitter_attributes].properties.emission_range.Reset(0.f, tfxEmissionRangePreset, add_node); common.library->emitter_attributes[emitter_attributes].properties.emission_range.type = tfxProperty_emission_range;
		common.library->emitter_attributes[emitter_attributes].properties.roll.Reset(0.f, tfxAnglePreset, add_node); common.library->emitter_attributes[emitter_attributes].properties.roll.type = tfxProperty_emitter_roll;
		common.library->emitter_attributes[emitter_attributes].properties.pitch.Reset(0.f, tfxAnglePreset, add_node); common.library->emitter_attributes[emitter_attributes].properties.pitch.type = tfxProperty_emitter_roll;
		common.library->emitter_attributes[emitter_attributes].properties.yaw.Reset(0.f, tfxAnglePreset, add_node); common.library->emitter_attributes[emitter_attributes].properties.yaw.type = tfxProperty_emitter_roll;
		common.library->emitter_attributes[emitter_attributes].properties.splatter.Reset(0.f, tfxDimensionsPreset, add_node); common.library->emitter_attributes[emitter_attributes].properties.splatter.type = tfxProperty_splatter;
		common.library->emitter_attributes[emitter_attributes].properties.emitter_width.Reset(0.f, tfxDimensionsPreset, add_node); common.library->emitter_attributes[emitter_attributes].properties.emitter_width.type = tfxProperty_emitter_width;
		common.library->emitter_attributes[emitter_attributes].properties.emitter_height.Reset(0.f, tfxDimensionsPreset, add_node); common.library->emitter_attributes[emitter_attributes].properties.emitter_height.type = tfxProperty_emitter_height;
		common.library->emitter_attributes[emitter_attributes].properties.emitter_depth.Reset(0.f, tfxDimensionsPreset, add_node); common.library->emitter_attributes[emitter_attributes].properties.emitter_depth.type = tfxProperty_emitter_depth;
		common.library->emitter_attributes[emitter_attributes].properties.arc_size.Reset(tfxRadians(360.f), tfxArcPreset, add_node); common.library->emitter_attributes[emitter_attributes].properties.arc_size.type = tfxProperty_arc_size;
		common.library->emitter_attributes[emitter_attributes].properties.arc_offset.Reset(0.f, tfxArcPreset, add_node); common.library->emitter_attributes[emitter_attributes].properties.arc_offset.type = tfxProperty_arc_offset;
		if(compile)
			common.library->CompilePropertyGraph(emitter_attributes);
	}

	void tfxEffectEmitter::ResetVariationGraphs(bool add_node, bool compile) {
		common.library->emitter_attributes[emitter_attributes].variation.life.Reset(0.f, tfxLifePreset, add_node); common.library->emitter_attributes[emitter_attributes].variation.life.type = tfxVariation_life;
		common.library->emitter_attributes[emitter_attributes].variation.amount.Reset(0.f, tfxAmountPreset, add_node); common.library->emitter_attributes[emitter_attributes].variation.amount.type = tfxVariation_amount;
		common.library->emitter_attributes[emitter_attributes].variation.velocity.Reset(0.f, tfxVelocityPreset, add_node); common.library->emitter_attributes[emitter_attributes].variation.velocity.type = tfxVariation_velocity;
		common.library->emitter_attributes[emitter_attributes].variation.width.Reset(0.f, tfxDimensionsPreset, add_node); common.library->emitter_attributes[emitter_attributes].variation.width.type = tfxVariation_width;
		common.library->emitter_attributes[emitter_attributes].variation.height.Reset(0.f, tfxDimensionsPreset, add_node); common.library->emitter_attributes[emitter_attributes].variation.height.type = tfxVariation_height;
		common.library->emitter_attributes[emitter_attributes].variation.weight.Reset(0.f, tfxWeightVariationPreset, add_node); common.library->emitter_attributes[emitter_attributes].variation.weight.type = tfxVariation_weight;
		common.library->emitter_attributes[emitter_attributes].variation.spin.Reset(0.f, tfxSpinVariationPreset, add_node); common.library->emitter_attributes[emitter_attributes].variation.spin.type = tfxVariation_spin;
		common.library->emitter_attributes[emitter_attributes].variation.noise_offset.Reset(0.f, tfxNoiseOffsetVariationPreset, add_node); common.library->emitter_attributes[emitter_attributes].variation.noise_offset.type = tfxVariation_noise_offset;
		common.library->emitter_attributes[emitter_attributes].variation.noise_resolution.Reset(300.f, tfxNoiseResolutionPreset, add_node); common.library->emitter_attributes[emitter_attributes].variation.noise_resolution.type = tfxVariation_noise_resolution;
		if(compile)
			common.library->CompileVariationGraph(emitter_attributes);
	}

	void tfxEffectEmitter::ResetOvertimeGraphs(bool add_node, bool compile) {
		common.library->emitter_attributes[emitter_attributes].overtime.velocity.Reset(1.f, tfxVelocityOvertimePreset, add_node); common.library->emitter_attributes[emitter_attributes].overtime.velocity.type = tfxOvertime_velocity;
		common.library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.Reset(1.f, tfxGlobalPercentPreset, add_node); common.library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.type = tfxOvertime_velocity_adjuster;
		common.library->emitter_attributes[emitter_attributes].overtime.width.Reset(1.f, tfxPercentOvertime, add_node); common.library->emitter_attributes[emitter_attributes].overtime.width.type = tfxOvertime_width;
		common.library->emitter_attributes[emitter_attributes].overtime.height.Reset(1.f, tfxPercentOvertime, add_node); common.library->emitter_attributes[emitter_attributes].overtime.height.type = tfxOvertime_height;
		common.library->emitter_attributes[emitter_attributes].overtime.weight.Reset(1.f, tfxPercentOvertime, add_node); common.library->emitter_attributes[emitter_attributes].overtime.weight.type = tfxOvertime_weight;
		common.library->emitter_attributes[emitter_attributes].overtime.spin.Reset(0.f, tfxSpinOvertimePreset, add_node); common.library->emitter_attributes[emitter_attributes].overtime.spin.type = tfxOvertime_spin;
		common.library->emitter_attributes[emitter_attributes].overtime.stretch.Reset(0.f, tfxPercentOvertime, add_node); common.library->emitter_attributes[emitter_attributes].overtime.stretch.type = tfxOvertime_stretch;
		common.library->emitter_attributes[emitter_attributes].overtime.red.Reset(1.f, tfxColorPreset, add_node); common.library->emitter_attributes[emitter_attributes].overtime.red.type = tfxOvertime_red;
		common.library->emitter_attributes[emitter_attributes].overtime.green.Reset(1.f, tfxColorPreset, add_node); common.library->emitter_attributes[emitter_attributes].overtime.green.type = tfxOvertime_green;
		common.library->emitter_attributes[emitter_attributes].overtime.blue.Reset(1.f, tfxColorPreset, add_node); common.library->emitter_attributes[emitter_attributes].overtime.blue.type = tfxOvertime_blue;
		common.library->emitter_attributes[emitter_attributes].overtime.blendfactor.Reset(1.f, tfxOpacityOvertimePreset, add_node); common.library->emitter_attributes[emitter_attributes].overtime.blendfactor.type = tfxOvertime_blendfactor;
		common.library->emitter_attributes[emitter_attributes].overtime.intensity.Reset(1.f, tfxIntensityOvertimePreset, add_node); common.library->emitter_attributes[emitter_attributes].overtime.intensity.type = tfxOvertime_intensity;
		common.library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance.Reset(30.f, tfxVelocityTurbulancePreset, add_node); common.library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance.type = tfxOvertime_velocity_turbulance;
		common.library->emitter_attributes[emitter_attributes].overtime.stretch.Reset(0.f, tfxPercentOvertime, add_node); common.library->emitter_attributes[emitter_attributes].overtime.stretch.type = tfxOvertime_stretch;
		common.library->emitter_attributes[emitter_attributes].overtime.direction_turbulance.Reset(0.f, tfxPercentOvertime, add_node); common.library->emitter_attributes[emitter_attributes].overtime.direction_turbulance.type = tfxOvertime_direction_turbulance;
		common.library->emitter_attributes[emitter_attributes].overtime.direction.Reset(0.f, tfxDirectionOvertimePreset, add_node); common.library->emitter_attributes[emitter_attributes].overtime.direction.type = tfxOvertime_direction;
		common.library->emitter_attributes[emitter_attributes].overtime.noise_resolution.Reset(1.f, tfxPercentOvertime, add_node); common.library->emitter_attributes[emitter_attributes].overtime.noise_resolution.type = tfxOvertime_noise_resolution;
		if(compile)
			common.library->CompileOvertimeGraph(emitter_attributes);
	}

	void tfxEffectEmitter::ResetEffectGraphs(bool add_node, bool compile) {
		ResetGlobalGraphs(add_node, compile);
	}

	void tfxEffectEmitter::ResetEmitterGraphs(bool add_node, bool compile) {
		ResetBaseGraphs(add_node, compile);
		ResetPropertyGraphs(add_node, compile);
		ResetVariationGraphs(add_node, compile);
		UpdateMaxLife();
		ResetOvertimeGraphs(add_node, compile);
	}

	void tfxEffectEmitter::InitialiseUninitialisedGraphs() {
		if (type == tfxEffectType) {
			if (common.library->global_graphs[global].life.nodes.size() == 0) common.library->global_graphs[global].life.Reset(1.f, tfxGlobalPercentPreset);
			if (common.library->global_graphs[global].amount.nodes.size() == 0) common.library->global_graphs[global].amount.Reset(1.f, tfxGlobalPercentPreset);
			if (common.library->global_graphs[global].velocity.nodes.size() == 0) common.library->global_graphs[global].velocity.Reset(1.f, tfxGlobalPercentPreset);
			if (common.library->global_graphs[global].width.nodes.size() == 0) common.library->global_graphs[global].width.Reset(1.f, tfxGlobalPercentPreset);
			if (common.library->global_graphs[global].height.nodes.size() == 0) common.library->global_graphs[global].height.Reset(1.f, tfxGlobalPercentPreset);
			if (common.library->global_graphs[global].weight.nodes.size() == 0) common.library->global_graphs[global].weight.Reset(1.f, tfxGlobalPercentPreset);
			if (common.library->global_graphs[global].spin.nodes.size() == 0) common.library->global_graphs[global].spin.Reset(1.f, tfxGlobalPercentPresetSigned);
			if (common.library->global_graphs[global].roll.nodes.size() == 0) common.library->global_graphs[global].roll.Reset(0.f, tfxAnglePreset);
			if (common.library->global_graphs[global].pitch.nodes.size() == 0) common.library->global_graphs[global].pitch.Reset(0.f, tfxAnglePreset);
			if (common.library->global_graphs[global].yaw.nodes.size() == 0) common.library->global_graphs[global].yaw.Reset(0.f, tfxAnglePreset);
			if (common.library->global_graphs[global].stretch.nodes.size() == 0) common.library->global_graphs[global].stretch.Reset(1.f, tfxGlobalPercentPreset);
			if (common.library->global_graphs[global].overal_scale.nodes.size() == 0) common.library->global_graphs[global].overal_scale.Reset(1.f, tfxGlobalPercentPreset);
			if (common.library->global_graphs[global].intensity.nodes.size() == 0) common.library->global_graphs[global].intensity.Reset(1.f, tfxGlobalPercentPreset);
			if (common.library->global_graphs[global].frame_rate.nodes.size() == 0) common.library->global_graphs[global].frame_rate.Reset(1.f, tfxGlobalPercentPreset);
			if (common.library->global_graphs[global].splatter.nodes.size() == 0) common.library->global_graphs[global].splatter.Reset(1.f, tfxGlobalPercentPreset);
			if (common.library->global_graphs[global].emitter_width.nodes.size() == 0) common.library->global_graphs[global].emitter_width.Reset(1.f, tfxGlobalPercentPreset);
			if (common.library->global_graphs[global].emitter_height.nodes.size() == 0) common.library->global_graphs[global].emitter_height.Reset(1.f, tfxGlobalPercentPreset);
			if (common.library->global_graphs[global].emitter_depth.nodes.size() == 0) common.library->global_graphs[global].emitter_depth.Reset(1.f, tfxGlobalPercentPreset);
		}

		if (type == tfxEmitterType) {
			if (common.library->emitter_attributes[emitter_attributes].base.life.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].base.life.Reset(1000.f, tfxLifePreset);
			if (common.library->emitter_attributes[emitter_attributes].base.amount.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].base.amount.Reset(1.f, tfxAmountPreset);
			if (common.library->emitter_attributes[emitter_attributes].base.velocity.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].base.velocity.Reset(0.f, tfxVelocityPreset);
			if (common.library->emitter_attributes[emitter_attributes].base.width.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].base.width.Reset(128.f, tfxDimensionsPreset);
			if (common.library->emitter_attributes[emitter_attributes].base.height.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].base.height.Reset(128.f, tfxDimensionsPreset);
			if (common.library->emitter_attributes[emitter_attributes].base.weight.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].base.weight.Reset(0.f, tfxWeightPreset);
			if (common.library->emitter_attributes[emitter_attributes].base.spin.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].base.spin.Reset(0.f, tfxSpinPreset);
			if (common.library->emitter_attributes[emitter_attributes].base.noise_offset.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].base.noise_offset.Reset(0.f, tfxGlobalPercentPreset);

			if (common.library->emitter_attributes[emitter_attributes].properties.roll.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].properties.roll.Reset(0.f, tfxAnglePreset);
			if (common.library->emitter_attributes[emitter_attributes].properties.pitch.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].properties.pitch.Reset(0.f, tfxAnglePreset);
			if (common.library->emitter_attributes[emitter_attributes].properties.yaw.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].properties.yaw.Reset(0.f, tfxAnglePreset);
			if (common.library->emitter_attributes[emitter_attributes].properties.emission_pitch.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].properties.emission_pitch.Reset(0.f, tfxAnglePreset);
			if (common.library->emitter_attributes[emitter_attributes].properties.emission_yaw.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].properties.emission_yaw.Reset(0.f, tfxAnglePreset);
			if (common.library->emitter_attributes[emitter_attributes].properties.emission_range.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].properties.emission_range.Reset(0.f, tfxEmissionRangePreset);
			if (common.library->emitter_attributes[emitter_attributes].properties.splatter.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].properties.splatter.Reset(0.f, tfxDimensionsPreset);
			if (common.library->emitter_attributes[emitter_attributes].properties.emitter_width.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].properties.emitter_width.Reset(0.f, tfxDimensionsPreset);
			if (common.library->emitter_attributes[emitter_attributes].properties.emitter_height.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].properties.emitter_height.Reset(0.f, tfxDimensionsPreset);
			if (common.library->emitter_attributes[emitter_attributes].properties.emitter_depth.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].properties.emitter_depth.Reset(0.f, tfxDimensionsPreset);
			if (common.library->emitter_attributes[emitter_attributes].properties.arc_size.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].properties.arc_size.Reset(tfxRadians(360.f), tfxArcPreset);
			if (common.library->emitter_attributes[emitter_attributes].properties.arc_offset.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].properties.arc_offset.Reset(0.f, tfxArcPreset);

			if (common.library->emitter_attributes[emitter_attributes].variation.life.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].variation.life.Reset(0.f, tfxLifePreset);
			if (common.library->emitter_attributes[emitter_attributes].variation.amount.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].variation.amount.Reset(0.f, tfxAmountPreset);
			if (common.library->emitter_attributes[emitter_attributes].variation.velocity.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].variation.velocity.Reset(0.f, tfxVelocityPreset);
			if (common.library->emitter_attributes[emitter_attributes].variation.weight.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].variation.weight.Reset(0.f, tfxVelocityPreset);
			if (common.library->emitter_attributes[emitter_attributes].variation.width.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].variation.width.Reset(0.f, tfxDimensionsPreset);
			if (common.library->emitter_attributes[emitter_attributes].variation.height.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].variation.height.Reset(0.f, tfxDimensionsPreset);
			if (common.library->emitter_attributes[emitter_attributes].variation.weight.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].variation.weight.Reset(0.f, tfxWeightVariationPreset);
			if (common.library->emitter_attributes[emitter_attributes].variation.spin.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].variation.spin.Reset(0.f, tfxSpinVariationPreset);
			if (common.library->emitter_attributes[emitter_attributes].variation.noise_offset.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].variation.noise_offset.Reset(0.f, tfxNoiseOffsetVariationPreset);
			if (common.library->emitter_attributes[emitter_attributes].variation.noise_resolution.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].variation.noise_resolution.Reset(300.f, tfxNoiseResolutionPreset);

			if (common.library->emitter_attributes[emitter_attributes].overtime.velocity.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.velocity.Reset(1.f, tfxVelocityOvertimePreset);
			if (common.library->emitter_attributes[emitter_attributes].overtime.width.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.width.Reset(1.f, tfxPercentOvertime);
			if (common.library->emitter_attributes[emitter_attributes].overtime.height.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.height.Reset(1.f, tfxPercentOvertime);
			if (common.library->emitter_attributes[emitter_attributes].overtime.weight.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.weight.Reset(1.f, tfxPercentOvertime);
			if (common.library->emitter_attributes[emitter_attributes].overtime.spin.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.spin.Reset(1.f, tfxSpinOvertimePreset);
			if (common.library->emitter_attributes[emitter_attributes].overtime.stretch.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.stretch.Reset(0.f, tfxPercentOvertime);
			if (common.library->emitter_attributes[emitter_attributes].overtime.red.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.red.Reset(1.f, tfxColorPreset);
			if (common.library->emitter_attributes[emitter_attributes].overtime.green.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.green.Reset(1.f, tfxColorPreset);
			if (common.library->emitter_attributes[emitter_attributes].overtime.blue.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.blue.Reset(1.f, tfxColorPreset);
			if (common.library->emitter_attributes[emitter_attributes].overtime.blendfactor.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.blendfactor.Reset(1.f, tfxOpacityOvertimePreset);
			if (common.library->emitter_attributes[emitter_attributes].overtime.intensity.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.intensity.Reset(1.f, tfxIntensityOvertimePreset);
			if (common.library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance.Reset(0.f, tfxFrameratePreset);
			if (common.library->emitter_attributes[emitter_attributes].overtime.direction_turbulance.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.direction_turbulance.Reset(0.f, tfxPercentOvertime);
			if (common.library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.Reset(1.f, tfxGlobalPercentPreset);
			if (common.library->emitter_attributes[emitter_attributes].overtime.direction.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.direction.Reset(0.f, tfxDirectionOvertimePreset);
			if (common.library->emitter_attributes[emitter_attributes].overtime.noise_resolution.nodes.size() == 0) common.library->emitter_attributes[emitter_attributes].overtime.noise_resolution.Reset(1.f, tfxPercentOvertime);
		}
	}

	void tfxEffectEmitter::SetName(const char *n) {
		GetInfo().name = n;
	}

	void tfxEffectEmitter::ClearColors() {
		common.library->emitter_attributes[emitter_attributes].overtime.red.Clear();
		common.library->emitter_attributes[emitter_attributes].overtime.green.Clear();
		common.library->emitter_attributes[emitter_attributes].overtime.blue.Clear();
	}

	void tfxEffectEmitter::AddColorOvertime(float frame, tfxRGB color) {
		common.library->emitter_attributes[emitter_attributes].overtime.red.AddNode(frame, color.r);
		common.library->emitter_attributes[emitter_attributes].overtime.green.AddNode(frame, color.g);
		common.library->emitter_attributes[emitter_attributes].overtime.blue.AddNode(frame, color.b);
	}

	void tfxEffectEmitter::ReSeed(uint64_t seed) {
		if (seed == 0) {
			seed = 0xFFFFFFFFFFF;
		}
		random_generation.ReSeed(seed, seed / 2);
	}

	 void tfxEffectEmitter::SetUserData(void *data) {
		user_data = data;
	}

	void* tfxEffectEmitter::GetUserData() {
		return user_data;
	}

	void tfxEffectEmitter::SetTimeout(float frames) {
		common.timeout = frames;
		for (auto &sub : GetInfo().sub_effectors) {
			sub.SetTimeout(frames);
		}
	}

	inline tfxEffectEmitterInfo &tfxEffectEmitter::GetInfo() {
		return common.library->GetInfo(*this);
	}

	inline tfxEmitterProperties &tfxEffectEmitter::GetProperties() {
		return common.library->GetProperties(property_index);
	}

	bool tfxEffectEmitter::HasSingle() {
		for (auto &e : GetInfo().sub_effectors) {
			if (e.common.property_flags & tfxEmitterPropertyFlags_single)
				return true;
		}
		return false;
	}

	bool tfxEffectEmitter::RenameSubEffector(tfxEffectEmitter &emitter, const char *new_name) {
		if (!NameExists(emitter, new_name) && strlen(new_name) > 0) {
			emitter.SetName(new_name);
			common.library->UpdateEffectPaths();
			return true;
		}

		return false;
	}

	bool tfxEffectEmitter::NameExists(tfxEffectEmitter &emitter, const char *name) {
		for (auto &e : GetInfo().sub_effectors) {
			if (&emitter != &e) {
				if (e.GetInfo().name == name) {
					return true;
				}
			}
		}

		return false;
	}

	void tfxEffectEmitter::ReIndex() {
		tfxU32 index = 0;
		for (auto &e : GetInfo().sub_effectors) {
			e.library_index = index++;
			e.parent = this;
			e.ReIndex();
		}
	}

	void tfxEffectEmitter::CountChildren(int &emitters, int &effects) {
		tmpStack(tfxEffectEmitter*, stack);
		stack.push_back(this);
		emitters = 0;
		effects = 0;
		while (!stack.empty()) {
			tfxEffectEmitter *current = stack.pop_back();
			if (current->type == tfxEffectType)
				effects++;
			else if (current->type == tfxEmitterType)
				emitters++;
			for (auto &sub : current->GetInfo().sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	tfxEffectEmitter* tfxEffectEmitter::GetRootEffect() {
		if (!parent || parent->type == tfxFolder)
			return this;
		tfxEffectEmitter *p = parent;
		tfxU32 timeout = 0;
		while (p || ++timeout < 100) {
			if (!p->parent)
				return p;
			p = p->parent;
		}
		return nullptr;
	}

	bool tfxEffectEmitter::IsRootEffect() {
		if (type != tfxEffectType) return false;
		if (type == tfxEffectType && !parent) return true;
		if (parent && parent->type == tfxFolder) return true;
		return false;
	}

	void tfxEffectEmitter::ResetParents() {
		parent = nullptr;
		for (auto &e : GetInfo().sub_effectors) {
			e.ResetParents();
		}
	}

	tfxEffectEmitter* tfxEffectEmitter::MoveUp(tfxEffectEmitter &emitter) {
		if (emitter.library_index > 0) {
			tfxU32 new_index = emitter.library_index - 1;
			std::swap(GetInfo().sub_effectors[emitter.library_index], GetInfo().sub_effectors[new_index]);
			ReIndex();
			emitter.common.library->UpdateEffectPaths();
			return &GetInfo().sub_effectors[new_index];
		}

		return nullptr;
	}

	tfxEffectEmitter* tfxEffectEmitter::MoveDown(tfxEffectEmitter &emitter) {
		if (emitter.library_index < GetInfo().sub_effectors.size() - 1) {
			tfxU32 new_index = emitter.library_index + 1;
			std::swap(GetInfo().sub_effectors[emitter.library_index], GetInfo().sub_effectors[new_index]);
			ReIndex();
			emitter.common.library->UpdateEffectPaths();
			return &GetInfo().sub_effectors[new_index];
		}
		return nullptr;
	}

	void tfxEffectEmitter::DeleteEmitter(tfxEffectEmitter *emitter) {
		tfxEffectLibrary *library = emitter->common.library;
		tmpStack(tfxEffectEmitter, stack);
		stack.push_back(*emitter);
		while (stack.size()) {
			tfxEffectEmitter &current = stack.pop_back();
			if (current.type == tfxEffectType && !current.parent) {
				library->FreeGlobal(current.global);
			}
			else if(current.type == tfxEmitterType) {
				library->FreeEmitterAttributes(current.emitter_attributes);
			}
			for (auto &sub : current.GetInfo().sub_effectors) {
				stack.push_back(sub);
			}
		}
		GetInfo().sub_effectors.erase(emitter);

		ReIndex();
		if(library)
			library->UpdateEffectPaths();
	}

	void tfxEffectEmitter::CleanUp() {
		if (GetInfo().sub_effectors.size()) {
			tmpStack(tfxEffectEmitter, stack);
			stack.push_back(*this);
			while (stack.size()) {
				tfxEffectEmitter current = stack.pop_back();
				if (current.type == tfxEffectType && !current.parent) {
					common.library->FreeGlobal(current.global);
				}
				else if(current.type == tfxEmitterType) {
					common.library->FreeEmitterAttributes(current.emitter_attributes);
				}
				for (auto &sub : current.GetInfo().sub_effectors) {
					stack.push_back(sub);
				}
				current.GetInfo().sub_effectors.clear();
				common.library->FreeProperties(current.property_index);
				common.library->FreeInfo(current.info_index);
			}
		}

		ReIndex();
	}

	void tfxEffectEmitter::Clone(tfxEffectEmitter &clone, tfxEffectEmitter *root_parent, tfxEffectLibrary *destination_library, tfxEffectCloningFlags flags) {
		//tfxU32 size = common.library->global_graphs[0].amount.lookup.values.capacity;
		clone = *this;
		clone.info_index = clone.common.library->CloneInfo(info_index, destination_library);
		clone.property_index = clone.common.library->CloneProperties(property_index, destination_library);
		clone.flags |= tfxEmitterStateFlags_enabled;
		if(!(flags & tfxEffectCloningFlags_keep_user_data))
			clone.user_data = nullptr;
		clone.common.library = destination_library;
		clone.GetInfo().sub_effectors.clear();

		if (type == tfxEffectType) {
			if (root_parent == &clone) {
				clone.global = flags & tfxEffectCloningFlags_clone_graphs ? common.library->CloneGlobal(global, destination_library) : clone.global = global;
				if(flags & tfxEffectCloningFlags_compile_graphs)
					clone.common.library->CompileGlobalGraph(clone.global);
			}
			else {
				if (!(flags & tfxEffectCloningFlags_force_clone_global)) {
					clone.global = root_parent->global;
				}
				else {
					clone.global = common.library->CloneGlobal(root_parent->global, destination_library);
					if(flags & tfxEffectCloningFlags_compile_graphs)
						clone.common.library->CompileGlobalGraph(clone.global);
				}
			}
		}
		else if(type == tfxEmitterType) {
			clone.emitter_attributes = flags & tfxEffectCloningFlags_clone_graphs ? common.library->CloneEmitterAttributes(emitter_attributes, destination_library) : emitter_attributes;
			clone.UpdateMaxLife();
			if (flags & tfxEffectCloningFlags_compile_graphs) {
				clone.common.library->CompilePropertyGraph(clone.emitter_attributes);
				clone.common.library->CompileBaseGraph(clone.emitter_attributes);
				clone.common.library->CompileVariationGraph(clone.emitter_attributes);
				clone.common.library->CompileOvertimeGraph(clone.emitter_attributes);
			}
		}

		for (auto &e : GetInfo().sub_effectors) {
			if (e.type == tfxEmitterType) {
				tfxEffectEmitter emitter_copy;
				e.Clone(emitter_copy, root_parent, destination_library, flags);
				if(!(flags & tfxEffectCloningFlags_keep_user_data))
					emitter_copy.user_data = nullptr;
				clone.AddEmitter(emitter_copy);
			}
			else if(e.type == tfxEffectType) {
				tfxEffectEmitter effect_copy;
				if(clone.type == tfxFolder)
					e.Clone(effect_copy, &effect_copy, destination_library, flags);
				else
					e.Clone(effect_copy, root_parent, destination_library, flags);
				if(!(flags & tfxEffectCloningFlags_keep_user_data))
					effect_copy.user_data = nullptr;
				clone.AddEffect(effect_copy);
			}
		}

	}

	bool PrepareEffectTemplate(tfxEffectLibrary &library, const char *name, tfxEffectTemplate &effect_template) {
		if (library.effect_paths.ValidName(name)) {
			library.PrepareEffectTemplate(name, effect_template);
			return true;
		}
		return false;
	}

	void tfxEffectEmitter::EnableAllEmitters() {
		for (auto &e : GetInfo().sub_effectors) {
			e.flags |= tfxEmitterStateFlags_enabled;
			e.EnableAllEmitters();
		}
	}

	void tfxEffectEmitter::EnableEmitter() {
		flags |= tfxEmitterStateFlags_enabled;
	}

	void tfxEffectEmitter::DisableAllEmitters() {
		for (auto &e : GetInfo().sub_effectors) {
			e.flags &= ~tfxEmitterStateFlags_enabled;
			e.DisableAllEmitters();
		}
	}

	void tfxEffectEmitter::DisableAllEmittersExcept(tfxEffectEmitter &emitter) {
		for (auto &e : GetInfo().sub_effectors) {
			if (e.library_index == emitter.library_index)
				e.flags |= tfxEmitterStateFlags_enabled;
			else
				e.flags &= ~tfxEmitterStateFlags_enabled;
		}
	}

	tfxGraph* tfxEffectEmitter::GetGraphByType(tfxGraphType type) {

		if (type < tfxGlobalCount) {
			return &((tfxGraph*)&common.library->global_graphs[global])[type];
		}
		else if (type >= tfxPropertyStart && type < tfxBaseStart) {
			int ref = type - tfxPropertyStart;
			return &((tfxGraph*)&common.library->emitter_attributes[emitter_attributes].properties)[ref];
		}
		else if (type >= tfxBaseStart && type < tfxVariationStart) {
			int ref = type - tfxBaseStart;
			return &((tfxGraph*)&common.library->emitter_attributes[emitter_attributes].base)[ref];
		}
		else if (type >= tfxVariationStart && type < tfxOvertimeStart) {
			int ref = type - tfxVariationStart;
			return &((tfxGraph*)&common.library->emitter_attributes[emitter_attributes].variation)[ref];
		}
		else if (type >= tfxOvertimeStart) {
			int ref = type - tfxOvertimeStart;
			return &((tfxGraph*)&common.library->emitter_attributes[emitter_attributes].overtime)[ref];
		}

		return nullptr;

	}

	tfxU32 tfxEffectEmitter::GetGraphIndexByType(tfxGraphType type) {

		if (type < tfxGlobalCount) {
			return global;
		}
		else {
			return emitter_attributes;
		}

	}

	void tfxEffectEmitter::FreeGraphs() {
		if (type == tfxEffectType) {
			common.library->global_graphs[global].life.Free();
			common.library->global_graphs[global].amount.Free();
			common.library->global_graphs[global].velocity.Free();
			common.library->global_graphs[global].width.Free();
			common.library->global_graphs[global].height.Free();
			common.library->global_graphs[global].weight.Free();
			common.library->global_graphs[global].spin.Free();
			common.library->global_graphs[global].roll.Free();
			common.library->global_graphs[global].pitch.Free();
			common.library->global_graphs[global].yaw.Free();
			common.library->global_graphs[global].stretch.Free();
			common.library->global_graphs[global].overal_scale.Free();
			common.library->global_graphs[global].intensity.Free();
			common.library->global_graphs[global].frame_rate.Free();
			common.library->global_graphs[global].splatter.Free();
			common.library->global_graphs[global].emitter_width.Free();
			common.library->global_graphs[global].emitter_height.Free();
			common.library->global_graphs[global].emitter_depth.Free();
		}

		if (type == tfxEmitterType) {
			common.library->emitter_attributes[emitter_attributes].properties.emission_pitch.Free();
			common.library->emitter_attributes[emitter_attributes].properties.emission_yaw.Free();
			common.library->emitter_attributes[emitter_attributes].properties.emission_range.Free();
			common.library->emitter_attributes[emitter_attributes].properties.roll.Free();
			common.library->emitter_attributes[emitter_attributes].properties.pitch.Free();
			common.library->emitter_attributes[emitter_attributes].properties.yaw.Free();
			common.library->emitter_attributes[emitter_attributes].properties.splatter.Free();
			common.library->emitter_attributes[emitter_attributes].properties.emitter_width.Free();
			common.library->emitter_attributes[emitter_attributes].properties.emitter_height.Free();
			common.library->emitter_attributes[emitter_attributes].properties.emitter_depth.Free();
			common.library->emitter_attributes[emitter_attributes].properties.arc_size.Free();
			common.library->emitter_attributes[emitter_attributes].properties.arc_offset.Free();

			common.library->emitter_attributes[emitter_attributes].base.life.Free();
			common.library->emitter_attributes[emitter_attributes].base.amount.Free();
			common.library->emitter_attributes[emitter_attributes].base.velocity.Free();
			common.library->emitter_attributes[emitter_attributes].base.width.Free();
			common.library->emitter_attributes[emitter_attributes].base.height.Free();
			common.library->emitter_attributes[emitter_attributes].base.weight.Free();
			common.library->emitter_attributes[emitter_attributes].base.spin.Free();
			common.library->emitter_attributes[emitter_attributes].base.noise_offset.Free();

			common.library->emitter_attributes[emitter_attributes].variation.life.Free();
			common.library->emitter_attributes[emitter_attributes].variation.amount.Free();
			common.library->emitter_attributes[emitter_attributes].variation.velocity.Free();
			common.library->emitter_attributes[emitter_attributes].variation.width.Free();
			common.library->emitter_attributes[emitter_attributes].variation.height.Free();
			common.library->emitter_attributes[emitter_attributes].variation.weight.Free();
			common.library->emitter_attributes[emitter_attributes].variation.spin.Free();
			common.library->emitter_attributes[emitter_attributes].variation.noise_offset.Free();
			common.library->emitter_attributes[emitter_attributes].variation.noise_resolution.Free();

			common.library->emitter_attributes[emitter_attributes].overtime.velocity.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.width.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.height.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.weight.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.spin.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.stretch.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.red.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.green.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.blue.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.blendfactor.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.intensity.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.direction_turbulance.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.direction.Free();
			common.library->emitter_attributes[emitter_attributes].overtime.noise_resolution.Free();
		}
	}

	tfxU32 tfxEffectEmitter::CountAllLookupValues() {
		tfxU32 count = 0;
		if (type == tfxEffectType) {
			count += common.library->CountGlobalLookUpValues(global);
			for (auto &emitter : GetInfo().sub_effectors) {
				count += common.library->CountEmitterLookUpValues(emitter.emitter_attributes);
			}
		}
		else if (type = tfxEmitterType) {
			count += common.library->CountEmitterLookUpValues(emitter_attributes);
		}
		return count;
	}

	void tfxEffectEmitter::CompileGraphs() {
		if (type == tfxEffectType) {
			for (tfxU32 t = (tfxU32)tfxGlobal_life; t != (tfxU32)tfxProperty_emission_pitch; ++t) {
				CompileGraph(*GetGraphByType(tfxGraphType(t)));
			}

			for (auto &emitter : GetInfo().sub_effectors) {
				for (tfxU32 t = (tfxU32)tfxProperty_emission_pitch; t != (tfxU32)tfxOvertime_velocity; ++t) {
					CompileGraph(*GetGraphByType((tfxGraphType)t));
				}

				for (tfxU32 t = (tfxU32)tfxOvertime_velocity; t != (tfxU32)tfxGraphMaxIndex; ++t) {
					CompileGraphOvertime(*emitter.GetGraphByType((tfxGraphType)t));
				}
			}
		}
	}

	tfxEffectEmitter& tfxEffectLibrary::operator[] (uint32_t index) {
		return effects[index];
	}

	bool tfxEffectLibrary::RenameEffect(tfxEffectEmitter &effect, const char *new_name) {
		if (!NameExists(effect, new_name) && strlen(new_name) > 0) {
			effect.SetName(new_name);
			UpdateEffectPaths();
			return true;
		}

		return false;
	}

	bool tfxEffectLibrary::NameExists(tfxEffectEmitter& effect, const char *name) {
		for (auto &e : effects) {
			if (effect.library_index != e.library_index) {
				if (e.GetInfo().name == name) {
					return true;
				}
			}
		}

		return false;
	}

	bool tfxEffectLibrary::NameExists2(tfxEffectEmitter& effect, const char *name) {
		for (auto &e : effects) {
			if (e.GetInfo().name == name) {
				return true;
			}
		}
		return false;
	}

	void tfxEffectLibrary::UpdateEffectPaths() {
		effect_paths.Clear();
		for (auto &e : effects) {
			tfxStr256 path = e.GetInfo().name;
			e.path_hash = tfxXXHash64::hash(path.c_str(), path.Length(), 0);
			AddPath(e, path);
		}
	}

	void tfxEffectLibrary::AddPath(tfxEffectEmitter &effect_emitter, tfxStr256 &path) {
		effect_paths.Insert(path, &effect_emitter);
		for (auto &sub : effect_emitter.GetInfo().sub_effectors) {
			tfxStr256 sub_path = path;
			sub_path.Appendf("/%s", sub.GetInfo().name.c_str());
			sub.path_hash = tfxXXHash64::hash(sub_path.c_str(), sub_path.Length(), 0);
			AddPath(sub, sub_path);
		}
	}

	tfxEffectEmitter &tfxEffectLibrary::AddEffect(tfxEffectEmitter &effect) {
		effect.library_index = effects.current_size;
		effect.type = tfxEffectType;
		effect.GetInfo().uid = ++uid;
		effect.common.library = this;
		effects.push_back(effect);
		ReIndex();
		UpdateEffectPaths();
		return effects.back();
	}

	tfxEffectEmitter &tfxEffectLibrary::AddFolder(tfxStr64 &name) {
		tfxEffectEmitter folder;
		folder.info_index = AddEffectEmitterInfo();
		folder.common.library = this;
		folder.GetInfo().name = name;
		folder.type = tfxFolder;
		folder.common.library = this;
		folder.GetInfo().uid = ++uid;
		effects.push_back(folder);
		ReIndex();
		UpdateEffectPaths();
		return effects.back();
	}

	tfxEffectEmitter &tfxEffectLibrary::AddFolder(tfxEffectEmitter &folder) {
		assert(folder.type == tfxFolder);			//Must be type tfxFolder if adding a folder
		folder.common.library = this;
		folder.GetInfo().uid = ++uid;
		effects.push_back(folder);
		ReIndex();
		UpdateEffectPaths();
		return effects.back();
	}

	tfxEffectEmitter* tfxEffectLibrary::GetEffect(tfxStr256 &path) {
		assert(effect_paths.ValidName(path));		//Effect was not found by that name
		return effect_paths.At(path);
	}

	tfxEffectEmitter* tfxEffectLibrary::GetEffect(const char *path) {
		assert(effect_paths.ValidName(path));		//Effect was not found by that name
		return effect_paths.At(path);
	}

	tfxEffectEmitter* tfxEffectLibrary::GetEffect(tfxKey key) {
		assert(effect_paths.ValidKey(key));			//Effect was not found by that key
		return effect_paths.At(key);
	}

	void tfxEffectLibrary::PrepareEffectTemplate(tfxStr256 path, tfxEffectTemplate &effect_template) {
		tfxEffectEmitter *effect = GetEffect(path);
		assert(effect);
		assert(effect->type == tfxEffectType);
		effect->Clone(effect_template.effect, &effect_template.effect, this);
		effect_template.AddPath(effect_template.effect, effect_template.effect.GetInfo().name.c_str());
	}

	void tfxEffectLibrary::ReIndex() {
		tfxU32 index = 0;
		for (auto &e : effects) {
			e.library_index = index++;
			e.parent = nullptr;
			e.ReIndex();
		}
	}

	void tfxEffectLibrary::UpdateParticleShapeReferences(tfxvec<tfxEffectEmitter> &effects, tfxU32 default_index) {
		for (auto &effect : effects) {
			if (effect.type == tfxFolder) {
				UpdateParticleShapeReferences(effect.GetInfo().sub_effectors, default_index);
			}
			else {
				for (auto &emitter : effect.GetInfo().sub_effectors) {
					if (particle_shapes.ValidIntName(emitter.GetProperties().shape_index))
						emitter.GetProperties().image = &particle_shapes.AtInt(emitter.GetProperties().shape_index);
					else
						emitter.GetProperties().image = &particle_shapes.AtInt(default_index);
					UpdateParticleShapeReferences(emitter.GetInfo().sub_effectors, default_index);
				}
			}
		}
	}

	tfxEffectEmitter* tfxEffectLibrary::MoveUp(tfxEffectEmitter &effect) {
		if (effect.library_index > 0) {
			tfxU32 new_index = effect.library_index - 1;
			std::swap(effects[effect.library_index], effects[new_index]);
			UpdateEffectPaths();
			ReIndex();
			return &effects[new_index];
		}
		return nullptr;
	}

	tfxEffectEmitter* tfxEffectLibrary::MoveDown(tfxEffectEmitter &effect) {
		if (effect.library_index < effects.size() - 1) {
			tfxU32 new_index = effect.library_index + 1;
			std::swap(effects[effect.library_index], effects[new_index]);
			UpdateEffectPaths();
			ReIndex();
			return &effects[new_index];
		}
		return nullptr;
	}

	void tfxEffectLibrary::BuildComputeShapeData(void* dst, tfxVec4(uv_lookup)(void *ptr, tfxComputeImageData &image_data, int offset)) {
		assert(dst);	//must be a valid pointer to a space in memory
		assert(particle_shapes.Size());		//There are no shapes to copy!
		tfxU32 index = 0;
		for (auto &shape : particle_shapes.data) {
			if (shape.animation_frames == 1) {
				tfxComputeImageData cs;
				cs.animation_frames = shape.animation_frames;
				cs.image_size = shape.image_size;
				cs.uv = uv_lookup(shape.ptr, cs, 0);
				shape_data.push_back(cs);
				shape.compute_shape_index = index++;
			}
			else {
				shape.compute_shape_index = index;
				for (int f = 0; f != shape.animation_frames; ++f) {
					tfxComputeImageData cs;
					cs.animation_frames = shape.animation_frames;
					cs.image_size = shape.image_size;
					cs.uv = uv_lookup(shape.ptr, cs, f);
					shape_data.push_back(cs);
					index++;
				}
			}
		}
	}

	void tfxEffectLibrary::CopyComputeShapeData(void* dst) {
		assert(shape_data.size());	//You must call BuildComputeShapeData first
		memcpy(dst, shape_data.data, shape_data.size() * sizeof(tfxComputeImageData));
	}

	void tfxEffectLibrary::CopyLookupIndexesData(void* dst) {
		assert(dst);	//must be a valid pointer to a space in memory
		assert(compiled_lookup_indexes.size());		//There is no data to copy, make sure a library has been loaded properly and it contains effects with emitters
		tfxGraphLookupIndex *test = static_cast<tfxGraphLookupIndex*>(dst);
		memcpy(dst, compiled_lookup_indexes.data, GetLookupIndexesSizeInBytes());
	}

	void tfxEffectLibrary::CopyLookupValuesData(void* dst) {
		assert(dst);	//must be a valid pointer to a space in memory
		assert(compiled_lookup_indexes.size());		//There is no data to copy, make sure a library has been loaded properly and it contains effects with emitters
		memcpy(dst, compiled_lookup_values.data, GetLookupValuesSizeInBytes());
	}

	tfxU32 tfxEffectLibrary::GetComputeShapeDataSizeInBytes() {
		tfxU32 frame_count = 0;
		for (auto &shape : particle_shapes.data) {
			frame_count += (tfxU32)shape.animation_frames;
		}
		return frame_count * sizeof(tfxComputeImageData);
	}

	tfxU32 tfxEffectLibrary::GetComputeShapeCount() {
		tfxU32 frame_count = 0;
		for (auto &shape : particle_shapes.data) {
			frame_count += (tfxU32)shape.animation_frames;
		}
		return frame_count;
	}

	tfxU32 tfxEffectLibrary::GetLookupIndexCount() {
		return compiled_lookup_indexes.size() * tfxOvertimeCount;
	}

	tfxU32 tfxEffectLibrary::GetLookupValueCount() {
		return compiled_lookup_values.size();
	}

	tfxU32 tfxEffectLibrary::GetLookupIndexesSizeInBytes() {
		return sizeof(tfxGraphLookupIndex) * tfxOvertimeCount * compiled_lookup_indexes.size();
	}

	tfxU32 tfxEffectLibrary::GetLookupValuesSizeInBytes() {
		return sizeof(float) * compiled_lookup_values.size();
	}

	void tfxEffectLibrary::RemoveShape(tfxU32 shape_index) {
		particle_shapes.RemoveInt(shape_index);
		for (auto &m : particle_shapes.map) {
			particle_shapes[m.index].shape_index = (tfxU32)m.key;
		}
	}

	void tfxEffectLibrary::DeleteEffect(tfxEffectEmitter *effect) {
		effects[effect->library_index].CleanUp();
		effects.erase(&effects[effect->library_index]);

		UpdateEffectPaths();
		ReIndex();
	}

	tfxU32 tfxEffectLibrary::AddGlobal() {
		if (free_global_graphs.size())
			return free_global_graphs.pop_back();
		tfxGlobalAttributes global;
		global.Initialise(&graph_node_allocator, &graph_lookup_allocator);
		global_graphs.push_back(global);
		return global_graphs.size() - 1;
	}

	tfxU32 tfxEffectLibrary::AddEmitterAttributes() {
		if (free_emitter_attributes.size())
			return free_emitter_attributes.pop_back();
		tfxEmitterAttributes attributes;
		attributes.Initialise(&graph_node_allocator, &graph_lookup_allocator);
		emitter_attributes.push_back(attributes);
		return emitter_attributes.size() - 1;
	}

	void tfxEffectLibrary::FreeGlobal(tfxU32 index) {
		assert(index < global_graphs.size());
		free_global_graphs.push_back(index);
		global_graphs[index].Free();
	}
	void tfxEffectLibrary::FreeEmitterAttributes(tfxU32 index) {
		assert(index < emitter_attributes.size());
		free_emitter_attributes.push_back(index);
		emitter_attributes[index].Free();
	}
	void tfxEffectLibrary::FreeProperties(tfxU32 index) {
		assert(index < emitter_properties.size());
		free_properties.push_back(index);
	}
	void tfxEffectLibrary::FreeInfo(tfxU32 index) {
		assert(index < effect_infos.size());
		free_infos.push_back(index);
	}

	tfxU32 tfxEffectLibrary::CountGlobalLookUpValues(tfxU32 index) {
		auto &global = global_graphs[index];
		tfxU32 count = 0;
		count += global.life.lookup.values.capacity;
		count += global.amount.lookup.values.capacity;
		count += global.velocity.lookup.values.capacity;
		count += global.width.lookup.values.capacity;
		count += global.height.lookup.values.capacity;
		count += global.weight.lookup.values.capacity;
		count += global.spin.lookup.values.capacity;
		count += global.stretch.lookup.values.capacity;
		count += global.overal_scale.lookup.values.capacity;
		count += global.intensity.lookup.values.capacity;
		count += global.frame_rate.lookup.values.capacity;
		count += global.splatter.lookup.values.capacity;
		count += global.roll.lookup.values.capacity;
		count += global.pitch.lookup.values.capacity;
		count += global.yaw.lookup.values.capacity;
		count += global.emitter_width.lookup.values.capacity;
		count += global.emitter_height.lookup.values.capacity;
		count += global.emitter_depth.lookup.values.capacity;
		return count;
	}

	tfxU32 tfxEffectLibrary::CountEmitterLookUpValues(tfxU32 index) {
		auto &attributes = emitter_attributes[index];
		tfxU32 count = 0;

		count += attributes.properties.emission_pitch.lookup.values.capacity;
		count += attributes.properties.emission_yaw.lookup.values.capacity;
		count += attributes.properties.emission_range.lookup.values.capacity;
		count += attributes.properties.roll.lookup.values.capacity;
		count += attributes.properties.pitch.lookup.values.capacity;
		count += attributes.properties.yaw.lookup.values.capacity;
		count += attributes.properties.splatter.lookup.values.capacity;
		count += attributes.properties.emitter_width.lookup.values.capacity;
		count += attributes.properties.emitter_height.lookup.values.capacity;
		count += attributes.properties.emitter_depth.lookup.values.capacity;
		count += attributes.properties.arc_size.lookup.values.capacity;
		count += attributes.properties.arc_offset.lookup.values.capacity;

		count += attributes.base.life.lookup.values.capacity;
		count += attributes.base.amount.lookup.values.capacity;
		count += attributes.base.velocity.lookup.values.capacity;
		count += attributes.base.width.lookup.values.capacity;
		count += attributes.base.height.lookup.values.capacity;
		count += attributes.base.weight.lookup.values.capacity;
		count += attributes.base.spin.lookup.values.capacity;
		count += attributes.base.noise_offset.lookup.values.capacity;

		count += attributes.variation.life.lookup.values.capacity;
		count += attributes.variation.amount.lookup.values.capacity;
		count += attributes.variation.velocity.lookup.values.capacity;
		count += attributes.variation.width.lookup.values.capacity;
		count += attributes.variation.height.lookup.values.capacity;
		count += attributes.variation.weight.lookup.values.capacity;
		count += attributes.variation.spin.lookup.values.capacity;
		count += attributes.variation.noise_offset.lookup.values.capacity;
		count += attributes.variation.noise_resolution.lookup.values.capacity;

		count += attributes.overtime.velocity.lookup.values.capacity;
		count += attributes.overtime.width.lookup.values.capacity;
		count += attributes.overtime.height.lookup.values.capacity;
		count += attributes.overtime.weight.lookup.values.capacity;
		count += attributes.overtime.spin.lookup.values.capacity;
		count += attributes.overtime.stretch.lookup.values.capacity;
		count += attributes.overtime.red.lookup.values.capacity;
		count += attributes.overtime.green.lookup.values.capacity;
		count += attributes.overtime.blue.lookup.values.capacity;
		count += attributes.overtime.blendfactor.lookup.values.capacity;
		count += attributes.overtime.velocity_turbulance.lookup.values.capacity;
		count += attributes.overtime.direction_turbulance.lookup.values.capacity;
		count += attributes.overtime.velocity_adjuster.lookup.values.capacity;
		count += attributes.overtime.intensity.lookup.values.capacity;
		count += attributes.overtime.direction.lookup.values.capacity;
		count += attributes.overtime.noise_resolution.lookup.values.capacity;

		return count;
	}

	tfxU32 tfxEffectLibrary::CloneGlobal(tfxU32 source_index, tfxEffectLibrary *destination_library) {
		tfxU32 index = destination_library->AddGlobal();
		global_graphs[source_index].CopyToNoLookups(&destination_library->global_graphs[index]);
		return index;
	}

	tfxU32 tfxEffectLibrary::CloneEmitterAttributes(tfxU32 source_index, tfxEffectLibrary *destination_library) {
		tfxU32 index = destination_library->AddEmitterAttributes();
		emitter_attributes[source_index].properties.CopyToNoLookups(&destination_library->emitter_attributes[index].properties);
		emitter_attributes[source_index].base.CopyToNoLookups(&destination_library->emitter_attributes[index].base);
		emitter_attributes[source_index].variation.CopyToNoLookups(&destination_library->emitter_attributes[index].variation);
		emitter_attributes[source_index].overtime.CopyToNoLookups(&destination_library->emitter_attributes[index].overtime);
		return index;
	}

	tfxU32 tfxEffectLibrary::CloneInfo(tfxU32 source_index, tfxEffectLibrary *destination_library) {
		tfxU32 index = destination_library->AddEffectEmitterInfo();
		destination_library->effect_infos[index] = effect_infos[source_index];
		return index;
	}

	tfxU32 tfxEffectLibrary::CloneProperties(tfxU32 source_index, tfxEffectLibrary *destination_library) {
		tfxU32 index = destination_library->AddEmitterProperties();
		destination_library->emitter_properties[index] = emitter_properties[source_index];
		return index;
	}

	void tfxEffectLibrary::AddEmitterGraphs(tfxEffectEmitter& emitter) {
		emitter.emitter_attributes = AddEmitterAttributes();
	}

	void tfxEffectLibrary::AddEffectGraphs(tfxEffectEmitter& effect) {
		tfxEffectEmitter *root_effect = effect.GetRootEffect();
		if (root_effect == &effect)
			effect.global = AddGlobal();
		else
			effect.global = root_effect->global;
	}

	tfxU32 tfxEffectLibrary::AddAnimationSettings(tfxEffectEmitter& effect) {
		assert(effect.type == tfxEffectType);
		tfxAnimationSettings a;
		a.frames = 32;
		a.current_frame = 1;
		a.frame_offset = 0;
		a.extra_frames_count = 0;
		a.position = tfxVec2(0.f, 0.f);
		a.frame_size = tfxVec2(256.f, 256.f);
		a.playback_speed = 1.f;
		a.loop = false;
		a.seamless = false;
		a.seed = 0;
		a.zoom = 1.f;
		a.scale = 1.f;
		a.needs_recording = true;
		a.needs_exporting = 0;
		a.color_option = tfxExportColorOptions::tfxFullColor;
		a.export_option = tfxExportOptions::tfxSpriteSheet;
		a.export_with_transparency = true;
		a.camera_settings.camera_floor_height = -10.f;
		a.camera_settings.camera_fov = tfxRadians(60);
		a.camera_settings.camera_pitch = tfxRadians(-30.f);
		a.camera_settings.camera_yaw = tfxRadians(-90.f);
		a.camera_settings.camera_position = tfxVec3(0.f, 3.5f, 7.5f);
		a.camera_settings.camera_isometric = false;
		a.camera_settings.camera_isometric_scale = 5.f;
		a.camera_settings.camera_hide_floor = false;
		animation_settings.push_back(a);
		effect.GetInfo().animation_settings = animation_settings.size() - 1;
		return effect.GetInfo().animation_settings;
	}

	tfxU32 tfxEffectLibrary::AddPreviewCameraSettings(tfxEffectEmitter& effect) {
		assert(effect.type == tfxEffectType);
		tfxPreviewCameraSettings a;
		a.camera_settings.camera_floor_height = -10.f;
		a.camera_settings.camera_fov = tfxRadians(60);
		a.camera_settings.camera_pitch = tfxRadians(-30.f);
		a.camera_settings.camera_yaw = tfxRadians(-90.f);
		a.camera_settings.camera_position = tfxVec3(0.f, 3.5f, 7.5f);
		a.camera_settings.camera_isometric = false;
		a.camera_settings.camera_isometric_scale = 5.f;
		a.camera_settings.camera_hide_floor = false;
		a.effect_z_offset = 5.f;
		a.camera_speed = 6.f;
		a.attach_effect_to_camera = false;
		preview_camera_settings.push_back(a);
		effect.GetInfo().preview_camera_settings = preview_camera_settings.size() - 1;
		return effect.GetInfo().preview_camera_settings;
	}

	tfxU32 tfxEffectLibrary::AddPreviewCameraSettings() {
		tfxPreviewCameraSettings a;
		a.camera_settings.camera_floor_height = -10.f;
		a.camera_settings.camera_fov = tfxRadians(60);
		a.camera_settings.camera_pitch = tfxRadians(-30.f);
		a.camera_settings.camera_yaw = tfxRadians(-90.f);
		a.camera_settings.camera_position = tfxVec3(0.f, 3.5f, 7.5f);
		a.camera_settings.camera_isometric = false;
		a.camera_settings.camera_isometric_scale = 5.f;
		a.camera_settings.camera_hide_floor = false;
		a.effect_z_offset = 5.f;
		a.camera_speed = 6.f;
		a.attach_effect_to_camera = false;
		preview_camera_settings.push_back(a);
		return preview_camera_settings.size() - 1;
	}

	tfxU32 tfxEffectLibrary::AddEffectEmitterInfo() {
		tfxEffectEmitterInfo info;
		if (free_infos.size()) {
			return free_infos.pop_back();
		}
		effect_infos.push_back(info);
		return effect_infos.size() - 1;
	}

	tfxU32 tfxEffectLibrary::AddEmitterProperties() {
		tfxEmitterProperties properties;
		if (free_properties.size()) {
			return free_properties.pop_back();
		}
		emitter_properties.push_back(properties);
		return emitter_properties.size() - 1;
	}

	void tfxEffectLibrary::Clear() {
		for (auto &e : effects) {
			e.FreeGraphs();
		}
		effects.free_all();
		effect_paths.FreeAll();
		particle_shapes.FreeAll();
		global_graphs.free_all();
		emitter_attributes.free_all();
		animation_settings.free_all();
		preview_camera_settings.free_all();
		emitter_properties.free_all();
		effect_infos.free_all();
		AddPreviewCameraSettings();

		graph_node_allocator.FreeAll();
		graph_lookup_allocator.FreeAll();

		free_global_graphs.free_all();
		free_emitter_attributes.free_all();
		free_animation_settings.free_all();
		free_preview_camera_settings.free_all();
		free_infos.free_all();
		free_properties.free_all();

		uid = 0;
	}

	void tfxEffectLibrary::UpdateEffectParticleStorage() {
		tfxParticleMemoryTools tools;
		for (auto &e : effects) {
			if (e.type != tfxFolder) {
				e.ResetAllBufferSizes();
				tools.ProcessEffect(e);
			}
			else {
				for (auto &sub : e.GetInfo().sub_effectors) {
					sub.ResetAllBufferSizes();
					tools.ProcessEffect(e);
				}
			}
		}
	}

	void tfxEffectLibrary::UpdateComputeNodes() {
		tfxU32 running_node_index = 0;
		tfxU32 running_value_index = 0;
		tmpStack(tfxEffectEmitter*, stack);
		all_nodes.clear();
		node_lookup_indexes.clear();
		compiled_lookup_values.clear();
		compiled_lookup_indexes.clear();
		for (auto &effect : effects) {
			stack.push_back(&effect);
			while(!stack.empty()) {
				tfxEffectEmitter *current = stack.pop_back();
				if (current->type == tfxFolder) {
					for (auto &sub : current->GetInfo().sub_effectors) {
						stack.push_back(&sub);
					}
					continue;
				}
				tfxEffectLookUpData lookup_data;
				tfxEffectLookUpData value_lookup_data;
				memset(&lookup_data, 0, sizeof(tfxEffectLookUpData));
				memset(&value_lookup_data, 0, sizeof(tfxEffectLookUpData));
				if (current->type == tfxEmitterType) {

					int offset = tfxGlobalCount + tfxPropertyCount + tfxBaseCount + tfxVariationCount;

					current->GetInfo().lookup_value_index = compiled_lookup_indexes.size();
					for (int i = 0; i != tfxOvertimeCount; ++i) {
						tfxGraph &graph = ((tfxGraph*)(&emitter_attributes[current->emitter_attributes].overtime))[i];
						tfxGraphLookupIndex &index = ((tfxGraphLookupIndex*)&lookup_data)[i];
						index.start_index = running_node_index;
						index.length = graph.nodes.size();
						index.max_life = graph.lookup.life;
						for (auto &node : graph.nodes) {
							all_nodes.push_back(node);
							running_node_index++;
						}

						tfxGraphLookupIndex value_index;
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
					current->GetInfo().lookup_node_index = node_lookup_indexes.size() - 1;

				}

				for (auto &sub : current->GetInfo().sub_effectors) {
					stack.push_back(&sub);
				}
			}
		}
	}

	void tfxEffectLibrary::CompileGraphsOfEffect(tfxEffectEmitter &effect, tfxU32 depth) {
		auto &info = effect.GetInfo();
		if (effect.type == tfxEffectType && depth == 0) {
			CompileGlobalGraph(effect.global);
		}
		else {
			CompileEmitterGraphs(effect.emitter_attributes);
		}
		for (auto &sub : info.sub_effectors) {
			CompileGraphsOfEffect(sub, ++depth);
		}
	}

	void tfxEffectLibrary::CompileAllGraphs() {
		for (auto &g : global_graphs) {
			CompileGraph(g.amount);
			CompileGraph(g.roll);
			CompileGraph(g.pitch);
			CompileGraph(g.yaw);
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
			CompileGraph(g.emitter_width);
			CompileGraph(g.emitter_height);
			CompileGraph(g.emitter_depth);
		}
		for (auto &g : emitter_attributes) {
			CompileGraph(g.properties.arc_offset);
			CompileGraph(g.properties.arc_size);
			CompileGraph(g.properties.emission_pitch);
			CompileGraph(g.properties.emission_yaw);
			CompileGraph(g.properties.emission_range);
			CompileGraph(g.properties.roll);
			CompileGraph(g.properties.pitch);
			CompileGraph(g.properties.yaw);
			CompileGraph(g.properties.emitter_width);
			CompileGraph(g.properties.emitter_height);
			CompileGraph(g.properties.emitter_depth);
			CompileGraph(g.properties.splatter);

			CompileGraph(g.base.amount);
			CompileGraph(g.base.width);
			CompileGraph(g.base.height);
			CompileGraph(g.base.life);
			CompileGraph(g.base.spin);
			CompileGraph(g.base.noise_offset);
			CompileGraph(g.base.velocity);
			CompileGraph(g.base.weight);

			CompileGraph(g.variation.amount);
			CompileGraph(g.variation.width);
			CompileGraph(g.variation.height);
			CompileGraph(g.variation.life);
			CompileGraph(g.variation.noise_offset);
			CompileGraph(g.variation.noise_resolution);
			CompileGraph(g.variation.spin);
			CompileGraph(g.variation.velocity);
			CompileGraph(g.variation.weight);

			CompileGraphOvertime(g.overtime.red);
			CompileGraphOvertime(g.overtime.green);
			CompileGraphOvertime(g.overtime.blue);
			CompileGraphOvertime(g.overtime.blendfactor);
			CompileGraphOvertime(g.overtime.intensity);
			CompileGraphOvertime(g.overtime.velocity_turbulance);
			CompileGraphOvertime(g.overtime.width);
			CompileGraphOvertime(g.overtime.height);
			CompileGraphOvertime(g.overtime.direction_turbulance);
			CompileGraphOvertime(g.overtime.spin);
			CompileGraphOvertime(g.overtime.stretch);
			CompileGraphOvertime(g.overtime.velocity);
			CompileGraph(g.overtime.velocity_adjuster);
			CompileGraphOvertime(g.overtime.weight);
			CompileGraphOvertime(g.overtime.direction);
			CompileGraphOvertime(g.overtime.noise_resolution);
		}
	}

	void tfxEffectLibrary::CompileGlobalGraph(tfxU32 index) {
		tfxGlobalAttributes &g = global_graphs[index];
		CompileGraph(g.amount);
		CompileGraph(g.roll);
		CompileGraph(g.pitch);
		CompileGraph(g.yaw);
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
		CompileGraph(g.emitter_width);
		CompileGraph(g.emitter_height);
		CompileGraph(g.emitter_depth);
	}

	void tfxEffectLibrary::CompileEmitterGraphs(tfxU32 index) {
		CompilePropertyGraph(index);
		CompileBaseGraph(index);
		CompileVariationGraph(index);
		CompileOvertimeGraph(index);
	}

	void tfxEffectLibrary::CompilePropertyGraph(tfxU32 index) {
		tfxPropertyAttributes &g = emitter_attributes[index].properties;
		CompileGraph(g.arc_offset);
		CompileGraph(g.arc_size);
		CompileGraph(g.emission_pitch);
		CompileGraph(g.emission_yaw);
		CompileGraph(g.emission_range);
		CompileGraph(g.roll);
		CompileGraph(g.pitch);
		CompileGraph(g.yaw);
		CompileGraph(g.emitter_width);
		CompileGraph(g.emitter_height);
		CompileGraph(g.splatter);
	}
	void tfxEffectLibrary::CompileBaseGraph(tfxU32 index) {
		tfxBaseAttributes &g = emitter_attributes[index].base;
		CompileGraph(g.amount);
		CompileGraph(g.width);
		CompileGraph(g.height);
		CompileGraph(g.life);
		CompileGraph(g.spin);
		CompileGraph(g.noise_offset);
		CompileGraph(g.velocity);
		CompileGraph(g.weight);
	}
	void tfxEffectLibrary::CompileVariationGraph(tfxU32 index) {
		tfxVariationAttributes &g = emitter_attributes[index].variation;
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
	void tfxEffectLibrary::CompileOvertimeGraph(tfxU32 index) {
		tfxOvertimeAttributes &g = emitter_attributes[index].overtime;
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
	void tfxEffectLibrary::CompileColorGraphs(tfxU32 index) {
		tfxOvertimeAttributes &g = emitter_attributes[index].overtime;
		CompileGraphOvertime(g.red);
		CompileGraphOvertime(g.green);
		CompileGraphOvertime(g.blue);
	}

	void tfxEffectLibrary::SetMinMaxData() {
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
		graph_min_max[tfxGlobal_effect_roll] = GetMinMaxGraphValues(tfxAnglePreset);
		graph_min_max[tfxGlobal_effect_pitch] = GetMinMaxGraphValues(tfxAnglePreset);
		graph_min_max[tfxGlobal_effect_yaw] = GetMinMaxGraphValues(tfxAnglePreset);
		graph_min_max[tfxGlobal_emitter_width] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_emitter_height] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_emitter_depth] = GetMinMaxGraphValues(tfxGlobalPercentPreset);

		graph_min_max[tfxProperty_emitter_roll] = GetMinMaxGraphValues(tfxAnglePreset);
		graph_min_max[tfxProperty_emission_pitch] = GetMinMaxGraphValues(tfxAnglePreset);
		graph_min_max[tfxProperty_emission_yaw] = GetMinMaxGraphValues(tfxAnglePreset);
		graph_min_max[tfxProperty_emission_range] = GetMinMaxGraphValues(tfxEmissionRangePreset);
		graph_min_max[tfxProperty_splatter] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxProperty_emitter_width] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxProperty_emitter_height] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxProperty_emitter_depth] = GetMinMaxGraphValues(tfxDimensionsPreset);
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

	float tfxEffectLibrary::LookupPreciseOvertimeNodeList(tfxGraphType graph_type, int lookup_node_index, float age, float life) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		tfxAttributeNode *lastec = nullptr;
		tfxGraphLookupIndex &lookup_data = ((tfxGraphLookupIndex*)&node_lookup_indexes[lookup_node_index])[graph_type];
		float min_y = graph_min_max[graph_type].y;
		float max_y = graph_min_max[graph_type].w;
		for (int i = lookup_data.start_index; i != lookup_data.start_index + lookup_data.length; ++i) {
			tfxAttributeNode &a = all_nodes[i];
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

	float tfxEffectLibrary::LookupPreciseNodeList(tfxGraphType graph_type, int lookup_node_index, float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		tfxAttributeNode *lastec = nullptr;
		tfxGraphLookupIndex &lookup_data = ((tfxGraphLookupIndex*)&node_lookup_indexes[lookup_node_index])[graph_type];
		float min_y = graph_min_max[graph_type].y;
		float max_y = graph_min_max[graph_type].w;
		for (int i = lookup_data.start_index; i != lookup_data.start_index + lookup_data.length; ++i) {
			tfxAttributeNode &a = all_nodes[i];
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

	float tfxEffectLibrary::LookupFastValueList(tfxGraphType graph_type, int lookup_node_index, float frame) {
		tfxGraphLookupIndex &lookup_data = ((tfxGraphLookupIndex*)&compiled_lookup_indexes[lookup_node_index])[graph_type];
		frame += lookup_data.start_index;
		tfxU32 end_frame = lookup_data.start_index + lookup_data.length - 1;
		frame = frame > end_frame ? end_frame : frame;
		return compiled_lookup_values[(tfxU32)frame];
	}

	float tfxEffectLibrary::LookupFastOvertimeValueList(tfxGraphType graph_type, int lookup_value_index, float age, float lifetime) {
		tfxGraphLookupIndex &lookup_data = ((tfxGraphLookupIndex*)&compiled_lookup_indexes[lookup_value_index])[graph_type - tfxOvertime_velocity];
		float frame = (float)lookup_data.start_index;
		if (lifetime)
			frame += (age / lifetime * lookup_data.max_life) / tfxLOOKUP_FREQUENCY_OVERTIME;
		if (frame < lookup_data.start_index + lookup_data.length - 1)
			return compiled_lookup_values[(tfxU32)frame];
		return compiled_lookup_values[lookup_data.start_index + lookup_data.length - 1];
	}

	tfxU32 tfxEffectLibrary::CountOfGraphsInUse() {
		return global_graphs.size() + emitter_attributes.size() - CountOfFreeGraphs();
	}

	tfxU32 tfxEffectLibrary::CountOfFreeGraphs() {
		return free_global_graphs.size() + free_emitter_attributes.size();
	}

	void tfxDataTypesDictionary::Init() {
		names_and_types.data.reserve(200);
		names_and_types.map.reserve(200);
		names_and_types.Insert("name", tfxString);
		names_and_types.Insert("image_index", tfxUint);
		names_and_types.Insert("image_handle_x", tfxFloat);
		names_and_types.Insert("image_handle_y", tfxFloat);
		names_and_types.Insert("spawn_amount", tfxUint);
		names_and_types.Insert("single_shot_limit", tfxUint);
		names_and_types.Insert("blend_mode", tfxSInt);
		names_and_types.Insert("image_start_frame", tfxFloat);
		names_and_types.Insert("image_end_frame", tfxFloat);
		names_and_types.Insert("image_frame_rate", tfxFloat);
		names_and_types.Insert("playback_speed", tfxFloat);

		names_and_types.Insert("emission_type", tfxSInt);
		names_and_types.Insert("emission_direction", tfxSInt);
		names_and_types.Insert("delay_spawning", tfxFloat);
		names_and_types.Insert("grid_rows", tfxFloat);
		names_and_types.Insert("grid_columns", tfxFloat);
		names_and_types.Insert("grid_depth", tfxFloat);
		names_and_types.Insert("loop_length", tfxFloat);
		names_and_types.Insert("emitter_handle_x", tfxFloat);
		names_and_types.Insert("emitter_handle_y", tfxFloat);
		names_and_types.Insert("emitter_handle_z", tfxFloat);
		names_and_types.Insert("end_behaviour", tfxSInt);
		names_and_types.Insert("angle_setting", tfxUint);
		names_and_types.Insert("angle_offset", tfxFloat);
		names_and_types.Insert("angle_offset_pitch", tfxFloat);
		names_and_types.Insert("angle_offset_yaw", tfxFloat);
		names_and_types.Insert("disable_billboard", tfxBool);
		names_and_types.Insert("billboard_option", tfxUint);
		names_and_types.Insert("vector_align_type", tfxUint);
		names_and_types.Insert("multiply_blend_factor", tfxFloat);
		names_and_types.Insert("sort_passes", tfxUint);

		names_and_types.Insert("random_color", tfxBool);
		names_and_types.Insert("relative_position", tfxBool);
		names_and_types.Insert("relative_angle", tfxBool);
		names_and_types.Insert("image_handle_auto_center", tfxBool);
		names_and_types.Insert("single", tfxBool);
		names_and_types.Insert("one_shot", tfxBool);
		names_and_types.Insert("spawn_on_grid", tfxBool);
		names_and_types.Insert("grid_spawn_clockwise", tfxBool);
		names_and_types.Insert("fill_area", tfxBool);
		names_and_types.Insert("emitter_handle_auto_center", tfxBool);
		names_and_types.Insert("edge_traversal", tfxBool);
		names_and_types.Insert("image_reverse_animation", tfxBool);
		names_and_types.Insert("image_play_once", tfxBool);
		names_and_types.Insert("image_animate", tfxBool);
		names_and_types.Insert("image_random_start_frame", tfxBool);
		names_and_types.Insert("global_uniform_size", tfxBool);
		names_and_types.Insert("base_uniform_size", tfxBool);
		names_and_types.Insert("lifetime_uniform_size", tfxBool);
		names_and_types.Insert("use_spawn_ratio", tfxBool);
		names_and_types.Insert("is_3d", tfxBool);
		names_and_types.Insert("draw_order_by_age", tfxBool);
		names_and_types.Insert("draw_order_by_depth", tfxBool);
		names_and_types.Insert("guaranteed_draw_order", tfxBool);

		//Animation settings
		names_and_types.Insert("animation_magenta_mask", tfxBool);
		names_and_types.Insert("frames", tfxUint);
		names_and_types.Insert("current_frame", tfxUint);
		names_and_types.Insert("frame_offset", tfxUint);
		names_and_types.Insert("extra_frames_count", tfxSInt);
		names_and_types.Insert("layer", tfxUint);
		names_and_types.Insert("position_x", tfxFloat);
		names_and_types.Insert("position_y", tfxFloat);
		names_and_types.Insert("frame_width", tfxFloat);
		names_and_types.Insert("frame_height", tfxFloat);
		names_and_types.Insert("loop", tfxBool);
		names_and_types.Insert("seamless", tfxBool);
		names_and_types.Insert("seed", tfxUint);
		names_and_types.Insert("zoom", tfxFloat);
		names_and_types.Insert("scale", tfxFloat);
		names_and_types.Insert("color_option", tfxSInt);
		names_and_types.Insert("export_option", tfxSInt);
		names_and_types.Insert("export_with_transparency", tfxBool);
		names_and_types.Insert("camera_position_x", tfxFloat);
		names_and_types.Insert("camera_position_y", tfxFloat);
		names_and_types.Insert("camera_position_z", tfxFloat);
		names_and_types.Insert("camera_pitch", tfxFloat);
		names_and_types.Insert("camera_yaw", tfxFloat);
		names_and_types.Insert("camera_fov", tfxFloat);
		names_and_types.Insert("camera_floor_height", tfxFloat);
		names_and_types.Insert("camera_isometric", tfxBool);
		names_and_types.Insert("camera_isometric_scale", tfxFloat);
		names_and_types.Insert("camera_hide_floor", tfxBool);
		names_and_types.Insert("camera_free_speed", tfxFloat);
		names_and_types.Insert("camera_ray_offset", tfxFloat);
		names_and_types.Insert("preview_camera_position_x", tfxFloat);
		names_and_types.Insert("preview_camera_position_y", tfxFloat);
		names_and_types.Insert("preview_camera_position_z", tfxFloat);
		names_and_types.Insert("preview_camera_pitch", tfxFloat);
		names_and_types.Insert("preview_camera_yaw", tfxFloat);
		names_and_types.Insert("preview_camera_fov", tfxFloat);
		names_and_types.Insert("preview_camera_floor_height", tfxFloat);
		names_and_types.Insert("preview_camera_isometric", tfxBool);
		names_and_types.Insert("preview_camera_isometric_scale", tfxFloat);
		names_and_types.Insert("preview_camera_speed", tfxFloat);
		names_and_types.Insert("preview_effect_z_offset", tfxFloat);
		names_and_types.Insert("preview_camera_hide_floor", tfxBool);
		names_and_types.Insert("preview_attach_effect_to_camera", tfxBool);

		//Editor config, move this to the editor
		names_and_types.Insert("only_play_selected_emitter", tfxBool);
		names_and_types.Insert("load_examples", tfxBool);
		names_and_types.Insert("load_last_file", tfxBool);
		names_and_types.Insert("load_last_file_path", tfxString);
		names_and_types.Insert("recent1", tfxString);
		names_and_types.Insert("recent2", tfxString);
		names_and_types.Insert("recent3", tfxString);
		names_and_types.Insert("recent4", tfxString);
		names_and_types.Insert("background_color_red", tfxFloat);
		names_and_types.Insert("background_color_green", tfxFloat);
		names_and_types.Insert("background_color_blue", tfxFloat);
		names_and_types.Insert("use_checker_background", tfxBool);
		names_and_types.Insert("preview_zoom", tfxFloat);
		names_and_types.Insert("updates_per_second", tfxFloat);
		names_and_types.Insert("background_image", tfxString);
		names_and_types.Insert("use_background_image", tfxBool);
		names_and_types.Insert("background_image_scale_x", tfxFloat);
		names_and_types.Insert("background_image_scale_y", tfxFloat);
		names_and_types.Insert("background_image_offset_x", tfxFloat);
		names_and_types.Insert("background_image_offset_y", tfxFloat);
		names_and_types.Insert("autoplay_effect", tfxSInt);
		names_and_types.Insert("sync_refresh_rate", tfxBool);
		names_and_types.Insert("window_maximised", tfxBool);
		names_and_types.Insert("window_width", tfxSInt);
		names_and_types.Insert("window_height", tfxSInt);
		names_and_types.Insert("window_x", tfxSInt);
		names_and_types.Insert("window_y", tfxSInt);
		names_and_types.Insert("show_emitter_positions", tfxBool);
		names_and_types.Insert("dpi_factor", tfxFloat);
		names_and_types.Insert("graph_lookup_mode", tfxSInt);
		names_and_types.Insert("show_tool_tips", tfxBool);
		names_and_types.Insert("preview_trail_mode", tfxBool);
		names_and_types.Insert("try_autorecover", tfxBool);
		names_and_types.Insert("autorecovery_file", tfxString);
		names_and_types.Insert("draw_outlines", tfxBool);
		initialised = true;
	}

	int ValidateEffectPackage(const char *filename) {
		tfxPackage package;
		tfxErrorFlags status = LoadPackage(filename, package);
		if (status) return status;					//returns 1 to 4 if it's an invalid package format

		tfxEntryInfo *data_txt = package.GetFile("data.txt");
		if (!data_txt) return tfxErrorCode_data_could_not_be_loaded;					//Unable to load the the data.txt file in the package

		return 0;
	}

	void AssignGraphData(tfxEffectEmitter &effect, tfxStack<tfxStr64> &values) {
		if (values.size() > 0) {
			if (values[0] == "global_amount") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].amount.AddNode(n); }
			if (values[0] == "global_effect_angle") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].roll.AddNode(n); }
			if (values[0] == "global_effect_roll") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].roll.AddNode(n); }
			if (values[0] == "global_effect_pitch") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].pitch.AddNode(n); }
			if (values[0] == "global_effect_yaw") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].yaw.AddNode(n); }
			if (values[0] == "global_frame_rate") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].frame_rate.AddNode(n); }
			if (values[0] == "global_height") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].height.AddNode(n); }
			if (values[0] == "global_width") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].width.AddNode(n); }
			if (values[0] == "global_life") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].life.AddNode(n); }
			if (values[0] == "global_opacity") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].intensity.AddNode(n); }
			if (values[0] == "global_spin") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].spin.AddNode(n); }
			if (values[0] == "global_splatter") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].splatter.AddNode(n); }
			if (values[0] == "global_stretch") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].stretch.AddNode(n); }
			if (values[0] == "global_overal_scale") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].overal_scale.AddNode(n); }
			if (values[0] == "global_weight") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].weight.AddNode(n); }
			if (values[0] == "global_velocity") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].velocity.AddNode(n); }
			if (values[0] == "global_emitter_width") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].emitter_width.AddNode(n); }
			if (values[0] == "global_emitter_height") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].emitter_height.AddNode(n); }
			if (values[0] == "global_emitter_depth") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->global_graphs[effect.global].emitter_depth.AddNode(n); }

			if (values[0] == "base_arc_offset") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.arc_offset.AddNode(n); }
			if (values[0] == "base_arc_size") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.arc_size.AddNode(n); }
			if (values[0] == "base_emission_angle") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.emission_pitch.AddNode(n); }
			if (values[0] == "base_emission_range") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.emission_range.AddNode(n); }
			if (values[0] == "base_emitter_height") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.emitter_height.AddNode(n); }
			if (values[0] == "base_emitter_width") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.emitter_width.AddNode(n); }
			if (values[0] == "base_splatter") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.splatter.AddNode(n); }

			if (values[0] == "property_arc_offset") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.arc_offset.AddNode(n); }
			if (values[0] == "property_arc_size") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.arc_size.AddNode(n); }
			if (values[0] == "property_emitter_angle") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.roll.AddNode(n); }
			if (values[0] == "property_emitter_roll") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.roll.AddNode(n); }
			if (values[0] == "property_emitter_pitch") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.pitch.AddNode(n); }
			if (values[0] == "property_emitter_yaw") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.yaw.AddNode(n); }
			if (values[0] == "property_emission_angle") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.emission_pitch.AddNode(n); }
			if (values[0] == "property_emission_pitch") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.emission_pitch.AddNode(n); }
			if (values[0] == "property_emission_yaw") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.emission_yaw.AddNode(n); }
			if (values[0] == "property_emission_range") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.emission_range.AddNode(n); }
			if (values[0] == "property_emitter_height") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.emitter_height.AddNode(n); }
			if (values[0] == "property_emitter_width") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.emitter_width.AddNode(n); }
			if (values[0] == "property_emitter_depth") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.emitter_depth.AddNode(n); }
			if (values[0] == "property_splatter") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].properties.splatter.AddNode(n); }

			if (values[0] == "base_amount") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].base.amount.AddNode(n); }
			if (values[0] == "base_life") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].base.life.AddNode(n); }
			if (values[0] == "base_height") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].base.height.AddNode(n); }
			if (values[0] == "base_width") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].base.width.AddNode(n); }
			if (values[0] == "base_spin") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].base.spin.AddNode(n); }
			if (values[0] == "base_noise_offset") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].base.noise_offset.AddNode(n); }
			if (values[0] == "base_velocity") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].base.velocity.AddNode(n); }
			if (values[0] == "base_weight") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].base.weight.AddNode(n); }

			if (values[0] == "variation_amount") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].variation.amount.AddNode(n); }
			if (values[0] == "variation_height") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].variation.height.AddNode(n); }
			if (values[0] == "variation_width") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].variation.width.AddNode(n); }
			if (values[0] == "variation_life") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].variation.life.AddNode(n); }
			if (values[0] == "variation_velocity") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].variation.velocity.AddNode(n); }
			if (values[0] == "variation_weight") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].variation.weight.AddNode(n); }
			if (values[0] == "variation_spin") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].variation.spin.AddNode(n); }
			if (values[0] == "variation_motion_randomness") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].variation.noise_offset.AddNode(n); }
			if (values[0] == "variation_noise_offset") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].variation.noise_offset.AddNode(n); }
			if (values[0] == "variation_noise_resolution") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].variation.noise_resolution.AddNode(n); }

			if (values[0] == "overtime_red") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.red.AddNode(n); }
			if (values[0] == "overtime_green") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.green.AddNode(n); }
			if (values[0] == "overtime_blue") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.blue.AddNode(n); }
			if (values[0] == "overtime_opacity") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.blendfactor.AddNode(n); }
			if (values[0] == "overtime_intensity") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.intensity.AddNode(n); }
			if (values[0] == "overtime_velocity_turbulance") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.velocity_turbulance.AddNode(n); }
			if (values[0] == "overtime_spin") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.spin.AddNode(n); }
			if (values[0] == "overtime_stretch") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.stretch.AddNode(n); }
			if (values[0] == "overtime_velocity") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.velocity.AddNode(n); }
			if (values[0] == "overtime_weight") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.weight.AddNode(n); }
			if (values[0] == "overtime_width") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.width.AddNode(n); }
			if (values[0] == "overtime_height") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.height.AddNode(n); }
			if (values[0] == "overtime_direction_turbulance") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.direction_turbulance.AddNode(n); }
			if (values[0] == "overtime_direction") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.direction.AddNode(n); }
			if (values[0] == "overtime_velocity_adjuster") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.velocity_adjuster.AddNode(n); }
			if (values[0] == "overtime_noise_resolution") { tfxAttributeNode n; AssignNodeData(n, values); effect.common.library->emitter_attributes[effect.emitter_attributes].overtime.noise_resolution.AddNode(n); }
		}
	}

	void AssignNodeData(tfxAttributeNode &n, tfxStack<tfxStr64> &values) {
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

	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, uint32_t value) {
		if (field == "image_index")
			effect.GetProperties().shape_index = value;
		if (field == "spawn_amount")
			effect.GetProperties().spawn_amount = value;
		if (field == "frames")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].frames = value;
		if (field == "current_frame")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].current_frame = value;
		if (field == "seed")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].seed = value;
		if (field == "layer")
			effect.GetProperties().layer = value;
		if (field == "frame_offset")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].frame_offset = value;
		if (field == "single_shot_limit")
			effect.GetProperties().single_shot_limit = value;
		if (field == "billboard_option")
			effect.GetProperties().billboard_option = (tfxBillboardingOptions)value;
		if (field == "vector_align_type")
			effect.GetProperties().vector_align_type = (tfxVectorAlignType)value;
		if (field == "angle_setting")
			effect.GetProperties().angle_settings = (tfxAngleSettingFlags)value;
		if (field == "sort_passes")
			effect.sort_passes = value;
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, int value) {
		if (field == "emission_type")
			effect.GetProperties().emission_type = (tfxEmissionType)value;
		if (field == "emission_direction")
			effect.GetProperties().emission_direction = (tfxEmissionDirection)value;
		if (field == "color_option")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].color_option = (tfxExportColorOptions)value;
		if (field == "export_option")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].export_option = (tfxExportOptions)value;
		if (field == "end_behaviour")
			effect.GetProperties().end_behaviour = (tfxLineTraversalEndBehaviour)value;
		if (field == "frame_offset")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].frame_offset = value;
		if (field == "extra_frames_count")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].extra_frames_count = value;
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, tfxStr &value) {
		if (field == "name") {
			effect.GetInfo().name = value;
		}
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, float value) {
		if (field == "position_x")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].position.x = value;
		if (field == "position_y")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].position.y = value;
		if (field == "frame_width")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].frame_size.x = value;
		if (field == "frame_height")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].frame_size.y = value;
		if (field == "zoom")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].zoom = value;
		if (field == "scale")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].scale = value;
		if (field == "playback_speed")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].playback_speed = value;
		if (field == "camera_position_x")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_position.x = value;
		if (field == "camera_position_y")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_position.y = value;
		if (field == "camera_position_z")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_position.z = value;
		if (field == "camera_pitch")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_pitch = value;
		if (field == "camera_yaw")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_yaw = value;
		if (field == "camera_fov")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_fov = value;
		if (field == "camera_floor_height")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_floor_height = value;
		if (field == "camera_isometric_scale")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_isometric_scale = value;
		if (field == "preview_camera_position_x")
			effect.common.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_position.x = value;
		if (field == "preview_camera_position_y")
			effect.common.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_position.y = value;
		if (field == "preview_camera_position_z")
			effect.common.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_position.z = value;
		if (field == "preview_camera_pitch")
			effect.common.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_pitch = value;
		if (field == "preview_camera_yaw")
			effect.common.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_yaw = value;
		if (field == "preview_camera_fov")
			effect.common.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_fov = value;
		if (field == "preview_camera_floor_height")
			effect.common.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_floor_height = value;
		if (field == "preview_camera_isometric_scale")
			effect.common.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_isometric_scale = value == 0 ? 5.f : value;
		if (field == "preview_effect_z_offset")
			effect.common.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].effect_z_offset = value;
		if (field == "preview_camera_speed")
			effect.common.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_speed = value;
		if (field == "image_handle_x")
			effect.GetProperties().image_handle.x = value;
		if (field == "image_handle_y")
			effect.GetProperties().image_handle.y = value;
		if (field == "delay_spawning")
			effect.GetProperties().delay_spawning = value;
		if (field == "grid_rows")
			effect.GetProperties().grid_points.x = value;
		if (field == "grid_columns")
			effect.GetProperties().grid_points.y = value;
		if (field == "grid_depth")
			effect.GetProperties().grid_points.z = value;
		if (field == "loop_length")
			effect.GetProperties().loop_length = value;
		if (field == "emitter_handle_x")
			effect.GetProperties().emitter_handle.x = value;
		if (field == "emitter_handle_y")
			effect.GetProperties().emitter_handle.y = value;
		if (field == "emitter_handle_z")
			effect.GetProperties().emitter_handle.z = value;
		if (field == "image_start_frame")
			effect.GetProperties().start_frame = value;
		if (field == "image_end_frame")
			effect.GetProperties().end_frame = value;
		if (field == "image_frame_rate")
			effect.GetProperties().frame_rate = value;
		if (field == "angle_offset")
			effect.GetProperties().angle_offsets.roll = value;
		if (field == "angle_offset_pitch")
			effect.GetProperties().angle_offsets.pitch = value;
		if (field == "angle_offset_yaw")
			effect.GetProperties().angle_offsets.yaw = value;
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, bool value) {
		if (field == "loop")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].loop = value;
		if (field == "seamless")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].seamless = value;
		if (field == "export_with_transparency")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].export_with_transparency = value;
		if (field == "camera_isometric")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_isometric = value;
		if (field == "camera_hide_floor")
			effect.common.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_hide_floor = value;
		if (field == "preview_attach_effect_to_camera")
			effect.common.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].attach_effect_to_camera = value;
		if (field == "preview_camera_hide_floor")
			effect.common.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_hide_floor = value;
		if (field == "preview_camera_isometric")
			effect.common.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_isometric = value;
		if (field == "random_color")
			if (value) effect.common.property_flags |= tfxEmitterPropertyFlags_random_color; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_random_color;
		if (field == "relative_position")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_relative_position; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_relative_position;
		if (field == "relative_angle")
			if (value) effect.common.property_flags |= tfxEmitterPropertyFlags_relative_angle; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_relative_angle;
		if (field == "image_handle_auto_center")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_image_handle_auto_center; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_image_handle_auto_center;
		if (field == "single")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_single; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_single;
		//if (field == "one_shot")
			//if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_one_shot; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_one_shot;
		if (field == "spawn_on_grid")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_spawn_on_grid; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_spawn_on_grid;
		if (field == "grid_spawn_clockwise")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_grid_spawn_clockwise; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_grid_spawn_clockwise;
		if (field == "fill_area")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_fill_area; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_fill_area;
		if (field == "emitter_handle_auto_center")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_emitter_handle_auto_center; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_emitter_handle_auto_center;
		if (field == "edge_traversal")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_edge_traversal; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_edge_traversal;
		if (field == "image_reverse_animation")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_reverse_animation; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_reverse_animation;
		if (field == "image_play_once")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_play_once; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_play_once;
		if (field == "image_animate")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_animate; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_animate;
		if (field == "image_random_start_frame")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_random_start_frame; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_random_start_frame;
		if (field == "global_uniform_size")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_global_uniform_size; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_global_uniform_size;
		if (field == "base_uniform_size")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_base_uniform_size; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_base_uniform_size;
		if (field == "lifetime_uniform_size")
			if(value) effect.common.property_flags |= tfxEmitterPropertyFlags_lifetime_uniform_size; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_lifetime_uniform_size;
		if (field == "use_spawn_ratio")
			if (value) effect.common.property_flags |= tfxEmitterPropertyFlags_use_spawn_ratio; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_use_spawn_ratio;
		if (field == "is_3d")
			if (value) effect.common.property_flags |= tfxEmitterPropertyFlags_is_3d; else effect.common.property_flags &= ~tfxEmitterPropertyFlags_is_3d;
		if (field == "draw_order_by_age")
			if (value) effect.effect_flags |= tfxEffectPropertyFlags_age_order; else effect.effect_flags &= ~tfxEffectPropertyFlags_age_order;
		if (field == "draw_order_by_depth")
			if (value) effect.effect_flags |= tfxEffectPropertyFlags_depth_draw_order; else effect.effect_flags &= ~tfxEffectPropertyFlags_depth_draw_order;
		if (field == "guaranteed_draw_order")
			if (value) effect.effect_flags |= tfxEffectPropertyFlags_guaranteed_order; else effect.effect_flags &= ~tfxEffectPropertyFlags_guaranteed_order;
	}

	void StreamProperties(tfxEmitterProperties &property, tfxEmitterPropertyFlags &flags, tfxStr &file) {

		file.AddLine("image_index=%i", property.shape_index);
		file.AddLine("image_handle_x=%f", property.image_handle.x);
		file.AddLine("image_handle_y=%f", property.image_handle.y);
		file.AddLine("image_start_frame=%f", property.start_frame);
		file.AddLine("image_end_frame=%f", property.end_frame);
		file.AddLine("image_frame_rate=%f", property.frame_rate);
		file.AddLine("image_play_once=%i", (flags & tfxEmitterPropertyFlags_play_once));
		file.AddLine("image_reverse_animation=%i", (flags & tfxEmitterPropertyFlags_reverse_animation));
		file.AddLine("image_animate=%i", (flags & tfxEmitterPropertyFlags_animate));
		file.AddLine("image_random_start_frame=%i", (flags & tfxEmitterPropertyFlags_random_start_frame));
		file.AddLine("spawn_amount=%i", property.spawn_amount);
		file.AddLine("emission_type=%i", property.emission_type);
		file.AddLine("emission_direction=%i", property.emission_direction);
		file.AddLine("grid_rows=%f", property.grid_points.x);
		file.AddLine("grid_columns=%f", property.grid_points.y);
		file.AddLine("grid_depth=%f", property.grid_points.z);
		file.AddLine("delay_spawning=%f", property.delay_spawning);
		file.AddLine("loop_length=%f", property.loop_length);
		file.AddLine("emitter_handle_x=%f", property.emitter_handle.x);
		file.AddLine("emitter_handle_y=%f", property.emitter_handle.y);
		file.AddLine("emitter_handle_z=%f", property.emitter_handle.z);
		file.AddLine("end_behaviour=%i", property.end_behaviour);
		file.AddLine("random_color=%i", (flags & tfxEmitterPropertyFlags_random_color));
		file.AddLine("relative_position=%i", (flags & tfxEmitterPropertyFlags_relative_position));
		file.AddLine("relative_angle=%i", (flags & tfxEmitterPropertyFlags_relative_angle));
		file.AddLine("image_handle_auto_center=%i", (flags & tfxEmitterPropertyFlags_image_handle_auto_center));
		file.AddLine("single=%i", (flags & tfxEmitterPropertyFlags_single));
		file.AddLine("single_shot_limit=%i", property.single_shot_limit);
		file.AddLine("spawn_on_grid=%i", (flags & tfxEmitterPropertyFlags_spawn_on_grid));
		file.AddLine("grid_spawn_clockwise=%i", (flags & tfxEmitterPropertyFlags_grid_spawn_clockwise));
		file.AddLine("fill_area=%i", (flags & tfxEmitterPropertyFlags_fill_area));
		file.AddLine("emitter_handle_auto_center=%i", (flags & tfxEmitterPropertyFlags_emitter_handle_auto_center));
		file.AddLine("edge_traversal=%i", (flags & tfxEmitterPropertyFlags_edge_traversal));
		file.AddLine("angle_setting=%i", property.angle_settings);
		file.AddLine("angle_offset=%f", property.angle_offsets.roll);
		file.AddLine("angle_offset_pitch=%f", property.angle_offsets.pitch);
		file.AddLine("angle_offset_yaw=%f", property.angle_offsets.yaw);
		file.AddLine("global_uniform_size=%i", (flags & tfxEmitterPropertyFlags_global_uniform_size));
		file.AddLine("base_uniform_size=%i", (flags & tfxEmitterPropertyFlags_base_uniform_size));
		file.AddLine("lifetime_uniform_size=%i", (flags & tfxEmitterPropertyFlags_lifetime_uniform_size));
		file.AddLine("use_spawn_ratio=%i", (flags & tfxEmitterPropertyFlags_use_spawn_ratio));
		file.AddLine("is_3d=%i", (flags & tfxEmitterPropertyFlags_is_3d));
		file.AddLine("billboard_option=%i", property.billboard_option);
		file.AddLine("vector_align_type=%i", property.vector_align_type);
		file.AddLine("layer=%i", property.layer);

	}

	void StreamProperties(tfxEffectEmitter &effect, tfxStr &file) {
		file.AddLine("draw_order_by_age=%i", effect.effect_flags & tfxEffectPropertyFlags_age_order);
		file.AddLine("draw_order_by_depth=%i", effect.effect_flags & tfxEffectPropertyFlags_depth_draw_order);
		file.AddLine("guaranteed_draw_order=%i", effect.effect_flags & tfxEffectPropertyFlags_guaranteed_order);
		file.AddLine("sort_passes=%i", effect.sort_passes);
	}

	void StreamGraph(const char * name, tfxGraph &graph, tfxStr &file) {

		for (auto &n : graph.nodes) {
			file.AddLine("%s,%f,%f,%i,%f,%f,%f,%f", name, n.frame, n.value, (n.flags & tfxAttributeNodeFlags_is_curve), n.left.x, n.left.y, n.right.x, n.right.y);
		}

	}

	tfxRandom::tfxRandom() {
		ReSeed();
	}

	void tfxRandom::ReSeed() {
		seeds[0] = Millisecs();
		seeds[1] = Millisecs() * 2;
	}

	void tfxRandom::ReSeed(uint64_t seed1, uint64_t seed2) {
		seeds[0] = seed1;
		seeds[1] = seed2;
	}

	static bool CompareNodes(tfxAttributeNode &left, tfxAttributeNode &right) {
		return left.frame < right.frame;
	}

	bool tfxGraph::IsOvertimeGraph() {
		return type >= tfxOvertime_velocity && type != tfxOvertime_velocity_adjuster;
	}

	bool tfxGraph::IsGlobalGraph() {
		return type >= tfxGlobal_life && type <= tfxGlobal_splatter;
	}

	bool tfxGraph::IsAngleGraph() {
		return (type == tfxGlobal_effect_roll || type == tfxGlobal_effect_pitch || type == tfxGlobal_effect_yaw || type == tfxProperty_emission_pitch || type == tfxProperty_emission_yaw || type == tfxProperty_emission_range ||
			type == tfxProperty_emitter_roll || type == tfxProperty_emitter_pitch || type == tfxProperty_emitter_yaw || type == tfxProperty_arc_offset || type == tfxProperty_arc_size || type == tfxBase_spin || type == tfxVariation_spin || type == tfxOvertime_direction);
	}

	void tfxGraph::MultiplyAllValues(float scalar) {
		for (auto &node : nodes) {
			node.value *= scalar;
			node.left.y *= scalar;
			node.right.y *= scalar;
		}
	}

	void tfxGraph::CopyToNoLookups(tfxGraph *graph) {
		graph->min = min;
		graph->max = max;
		graph->graph_preset = graph_preset;
		graph->type = type;
		graph->effector = effector;
		graph->nodes = nodes;
		graph->index = index;
		graph->lookup.life = lookup.life;
	}

	float tfxAttributeNode::GetX() {
		return frame;
	}
	float tfxAttributeNode::GetY() {
		return value;
	}

	bool SetNode(tfxGraph &graph, tfxAttributeNode &node, float _frame, float _value, tfxAttributeNodeFlags flags, float _c0x, float _c0y, float _c1x, float _c1y) {
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

	bool SetNodeFrame(tfxGraph &graph, tfxAttributeNode &node, float &_frame) {
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

	bool SetNodeValue(tfxGraph &graph, tfxAttributeNode &node, float &_value) {
		node.value = _value;
		ClampNode(graph, node);
		_value = node.value;

		return false;
	}

	bool SetNode(tfxGraph &graph, tfxAttributeNode &node, float &frame, float &value) {
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

	void SetCurve(tfxGraph &graph, tfxAttributeNode &node, bool is_left_curve, float &frame, float &value) {
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

	bool MoveNode(tfxGraph &graph, tfxAttributeNode &node, float frame, float value, bool sort) {
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

	void ClampNode(tfxGraph &graph, tfxAttributeNode &node) {
		if (node.value < graph.min.y) node.value = graph.min.y;
		if (node.frame < graph.min.x) node.frame = graph.min.x;
		if (node.value > graph.max.y) node.value = graph.max.y;
		if (node.frame > graph.max.x) node.frame = graph.max.x;
	}

	void ClampGraph(tfxGraph &graph) {
		for (auto &node : graph.nodes) {
			ClampNode(graph, node);
			if (node.flags & tfxAttributeNodeFlags_is_curve) {
				ClampCurve(graph, node.left, node);
				ClampCurve(graph, node.right, node);
			}
		}
	}

	void ClampCurve(tfxGraph &graph, tfxVec2 &p, tfxAttributeNode &node) {
		if (p.y < graph.min.y) p.y = graph.min.y;
		if (p.x < graph.min.x) p.x = graph.min.x;
		//if (p.y > graph.max.y) p.y = graph.max.y;
		if (p.x > graph.max.x) p.x = graph.max.x;

		tfxAttributeNode *next = graph.GetNextNode(node);
		if (next) {
			if (p.x > next->frame) p.x = next->frame;
		}

		tfxAttributeNode *prev = graph.GetPrevNode(node);
		if (prev) {
			if (p.x < prev->frame) p.x = prev->frame;
		}
	}

	tfxGraph::tfxGraph() {
		min.x = 0.f;
		min.y = 0.f;
		max.x = 1000.f;
		max.y = 1000.f;
		effector = nullptr;
	}

	tfxGraph::tfxGraph(tfxMemoryArenaManager *node_allocator, tfxU32 bucket_size) { 
		min.x = 0.f;
		min.y = 0.f;
		max.x = 1000.f;
		max.y = 1000.f;
		effector = nullptr;
		nodes.allocator = node_allocator; 
		nodes.size_of_each_bucket = bucket_size; 
	}

	tfxGraph::~tfxGraph() {
		//Free();
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

	tfxVec2 GetQuadBezier(tfxVec2 p0, tfxVec2 p1, tfxVec2 p2, float t, float ymin, float ymax, bool clamp) {
		tfxVec2 b;
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

	tfxVec2 GetCubicBezier(tfxVec2 p0, tfxVec2 p1, tfxVec2 p2, tfxVec2 p3, float t, float ymin, float ymax, bool clamp) {
		tfxVec2 b;
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

	float GetBezierValue(const tfxAttributeNode *lastec, const tfxAttributeNode &a, float t, float ymin, float ymax) {
		if (lastec) {
			if (a.flags & tfxAttributeNodeFlags_is_curve) {
				if (lastec->flags & tfxAttributeNodeFlags_is_curve) {
					tfxVec2 p0(lastec->frame, lastec->value);
					tfxVec2 p1(lastec->right.x, lastec->right.y);
					tfxVec2 p2(a.left.x, a.left.y);
					tfxVec2 p3(a.frame, a.value);
					tfxVec2 value = GetCubicBezier(p0, p1, p2, p3, t, ymin, ymax);
					return value.y;
				}
				else {
					tfxVec2 p0(lastec->frame, lastec->value);
					tfxVec2 p1(a.left.x, a.left.y);
					tfxVec2 p2(a.frame, a.value);
					tfxVec2 value = GetQuadBezier(p0, p1, p2, t, ymin, ymax);
					return value.y;
				}
			}
			else if (lastec->flags & tfxAttributeNodeFlags_is_curve) {
				tfxVec2 p0(lastec->frame, lastec->value);
				tfxVec2 p1(lastec->right.x, lastec->right.y);
				tfxVec2 p2(a.frame, a.value);
				tfxVec2 value = GetQuadBezier(p0, p1, p2, t, ymin, ymax);
				return value.y;
			}
		}
		else {
			return 0;
		}

		return 0;
	}

	tfxAttributeNode* tfxGraph::AddNode(float _frame, float _value, tfxAttributeNodeFlags flags, float _c0x, float _c0y, float _c1x, float _c1y) {
		tfxAttributeNode node;

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

	void tfxGraph::AddNode(tfxAttributeNode &node) {
		for (auto &n : nodes) {
			if (n.frame == node.frame)
				return;
		}
		nodes.push_back(node);
		Sort();
		ReIndex();
	}

	tfxAttributeNode* tfxGraph::AddCoordNode(float _frame, float _value) {
		tfxAttributeNode node;

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
		tfxAttributeNode &n = nodes.push_back(node);
		if (Sort()) {
			ReIndex();
			return nodes.find(n);
		}

		ReIndex();
		return &n;
	}

	tfxAttributeNode* tfxGraph::InsertCoordNode(float _frame, float _value) {
		tfxAttributeNode node;

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
			tfxAttributeNode *last_node = nullptr;
			for (auto *n = nodes.begin() + 1; n != nodes.end(); ++n) {
				if (node.frame < n->frame)
					last_node = n;
				else
					break;
			}

			if (last_node) {
				tfxAttributeNode *r_value = nodes.insert(last_node, node);
				ReIndex();
				return r_value;
			}
		}

		tfxAttributeNode *r_value = &nodes.push_back(node);
		ReIndex();
		return r_value;
	}

	tfxAttributeNode* tfxGraph::InsertNode(float _frame, float _value) {
		tfxAttributeNode node;

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
			tfxAttributeNode *last_node = nullptr;
			for (auto *n = nodes.begin() + 1; n != nodes.end(); ++n) {
				if (node.frame < n->frame)
					last_node = n;
				else
					break;
			}

			if (last_node) {
				tfxAttributeNode *r_value = nodes.insert(last_node, node);
				ReIndex();
				return r_value;
			}
		}

		tfxAttributeNode *r_value = &nodes.push_back(node);
		ReIndex();
		return r_value;
	}

	void tfxGraph::SetNode(uint32_t i, float _frame, float _value, tfxAttributeNodeFlags flags, float _c0x, float _c0y, float _c1x, float _c1y) {
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

	tfxBucketArray<tfxAttributeNode>& tfxGraph::Nodes() {
		return nodes;
	}

	float tfxGraph::GetValue(float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		tfxAttributeNode *lastec = nullptr;
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

	tfxAttributeNode *tfxGraph::GetNextNode(tfxAttributeNode &node) {
		if (node.index < nodes.size() - 1) {
			return &nodes[node.index + 1];
		}

		return nullptr;
	}

	tfxAttributeNode *tfxGraph::GetPrevNode(tfxAttributeNode &node) {
		if (node.index > 0) {
			return &nodes[node.index - 1];
		}

		return nullptr;
	}

	tfxAttributeNode *tfxGraph::GetLastNode() {
		return &nodes.back();
	}

	float tfxGraph::GetRandomValue(float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		tfxAttributeNode *lastec = nullptr;
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

	float tfxGraph::GetValue(float age, float life) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		tfxAttributeNode *lastec = nullptr;
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

	float tfxGraph::GetFirstValue() {
		if (nodes.size())
			return nodes.front().value;
		return 0.f;
	}

	float* tfxGraph::LinkFirstValue() {
		if (nodes.size())
			return &nodes.front().value;
		return nullptr;
	}

	float tfxGraph::GetLastValue() {
		if (nodes.size())
			return nodes.back().value;

		return 0.f;
	}

	float tfxGraph::GetMaxValue() {
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

	float tfxGraph::GetMinValue() {
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

	float tfxGraph::GetLastFrame() {
		if (nodes.size())
			return nodes.back().frame;

		return 0.f;
	}

	tfxAttributeNode* tfxGraph::FindNode(const tfxAttributeNode &n) {
		return nodes.find(n);
	}

	void tfxGraph::ValidateCurves() {
		tfxU32 index = 0;
		tfxU32 last_index = nodes.size() - 1;
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

	void tfxGraph::DeleteNode(const tfxAttributeNode &n) {
		nodes.erase(&n);
	}

	void tfxGraph::Reset(float v, tfxGraphPreset preset, bool add_node) {
		nodes.clear();
		if (add_node)
			AddNode(0.f, v);
		switch (preset) {
		case tfxGraphPreset::tfxGlobalPercentPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 20.f };
			break;
		case tfxGraphPreset::tfxGlobalOpacityPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 1.f };
			break;
		case tfxGraphPreset::tfxGlobalPercentPresetSigned:
			min = { 0.f, -20.f }; max = { tfxMAX_FRAME, 20.f };
			break;
		case tfxGraphPreset::tfxAnglePreset:
			min = { 0.f, -1080.f }; max = { tfxMAX_FRAME, 1080.f };
			break;
		case tfxGraphPreset::tfxArcPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 360.f };
			break;
		case tfxGraphPreset::tfxEmissionRangePreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 360.f };
			break;
		case tfxGraphPreset::tfxDimensionsPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 4000.f };
			break;
		case tfxGraphPreset::tfxLifePreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 100000.f };
			break;
		case tfxGraphPreset::tfxAmountPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 5000.f };
			break;
		case tfxGraphPreset::tfxVelocityPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 10000.f };
			break;
		case tfxGraphPreset::tfxWeightPreset:
			min = { 0.f, -2500.f }; max = { tfxMAX_FRAME, 2500.f };
			break;
		case tfxGraphPreset::tfxWeightVariationPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 2500.f };
			break;
		case tfxGraphPreset::tfxNoiseOffsetVariationPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 1000.f };
			break;
		case tfxGraphPreset::tfxNoiseResolutionPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 10000.f };
			break;
		case tfxGraphPreset::tfxSpinPreset:
			min = { 0.f, -2000.f }; max = { tfxMAX_FRAME, 2000.f };
			break;
		case tfxGraphPreset::tfxSpinVariationPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 2000.f };
			break;
		case tfxGraphPreset::tfxDirectionVariationPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 22.5f };
			break;
		case tfxGraphPreset::tfxWeightOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 20.f };
			break;
		case tfxGraphPreset::tfxDirectionOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 4320.f };
			break;
		case tfxGraphPreset::tfxSpinOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 20.f };
			break;
		case tfxGraphPreset::tfxVelocityOvertimePreset:
			min = { 0.f, -20.f }; max = { 1.f, 20.f };
			break;
		case tfxGraphPreset::tfxPercentOvertime:
			min = { 0.f, 0.f }; max = { 1.f, 20.f };
			break;
		case tfxGraphPreset::tfxFrameratePreset:
			min = { 0.f, 0.f }; max = { 1.f, 200.f };
			break;
		case tfxGraphPreset::tfxVelocityTurbulancePreset:
			min = { 0.f, 0.f }; max = { 1.f, 2000.f };
			break;
		case tfxGraphPreset::tfxOpacityOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 1.f };
			break;
		case tfxGraphPreset::tfxColorPreset:
			min = { 0.f, 0.f }; max = { 1.f, 255.f };
			break;
		case tfxGraphPreset::tfxIntensityOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 5.f };
			break;
		}

		graph_preset = preset;
	}

	tfxVec4 GetMinMaxGraphValues(tfxGraphPreset preset) {
		tfxVec4 mm;
		switch (preset) {
		case tfxGraphPreset::tfxGlobalPercentPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 20.f };
			break;
		case tfxGraphPreset::tfxGlobalOpacityPreset:
			mm = { 0.f, 0.f , tfxMAX_FRAME, 1.f };
			break;
		case tfxGraphPreset::tfxGlobalPercentPresetSigned:
			mm = { 0.f, -20.f, tfxMAX_FRAME, 20.f };
			break;
		case tfxGraphPreset::tfxAnglePreset:
			mm = { 0.f, -1080.f, tfxMAX_FRAME, 1080.f };
			break;
		case tfxGraphPreset::tfxArcPreset:
			mm = { 0.f, 0.f , tfxMAX_FRAME, 360.f };
			break;
		case tfxGraphPreset::tfxEmissionRangePreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 360.f };
			break;
		case tfxGraphPreset::tfxDimensionsPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 4000.f };
			break;
		case tfxGraphPreset::tfxLifePreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 100000.f };
			break;
		case tfxGraphPreset::tfxAmountPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 5000.f };
			break;
		case tfxGraphPreset::tfxVelocityPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 10000.f };
			break;
		case tfxGraphPreset::tfxWeightPreset:
			mm = { 0.f, -2500.f, tfxMAX_FRAME, 2500.f };
			break;
		case tfxGraphPreset::tfxWeightVariationPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 2500.f };
			break;
		case tfxGraphPreset::tfxSpinPreset:
			mm = { 0.f, -2000.f, tfxMAX_FRAME, 2000.f };
			break;
		case tfxGraphPreset::tfxSpinVariationPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 2000.f };
			break;
		case tfxGraphPreset::tfxDirectionVariationPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 22.5f };
			break;
		case tfxGraphPreset::tfxWeightOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 20.f };
			break;
		case tfxGraphPreset::tfxDirectionOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 4320.f };
			break;
		case tfxGraphPreset::tfxSpinOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 20.f };
			break;
		case tfxGraphPreset::tfxVelocityOvertimePreset:
			mm = { 0.f, -20.f, 1.f, 20.f };
			break;
		case tfxGraphPreset::tfxPercentOvertime:
			mm = { 0.f, 0.f, 1.f, 20.f };
			break;
		case tfxGraphPreset::tfxFrameratePreset:
			mm = { 0.f, 0.f, 1.f, 200.f };
			break;
		case tfxGraphPreset::tfxVelocityTurbulancePreset:
			mm = { 0.f, 0.f, 1.f, 2000.f };
			break;
		case tfxGraphPreset::tfxOpacityOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 1.f };
			break;
		case tfxGraphPreset::tfxColorPreset:
			mm = { 0.f, 0.f, 1.f, 255.f };
			break;
		case tfxGraphPreset::tfxIntensityOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 5.f };
			break;
		}
		
		return mm;
	}

	void tfxGraph::DragValues(tfxGraphPreset preset, float &frame, float &value) {
		switch (preset) {
		case tfxGraphPreset::tfxOpacityOvertimePreset:
		case tfxGraphPreset::tfxGlobalPercentPreset:
		case tfxGraphPreset::tfxIntensityOvertimePreset:
			frame = 0.001f;
			value = 0.001f;
			break;
		case tfxGraphPreset::tfxDirectionOvertimePreset:
			frame = 0.001f;
			value = 0.1f;
			break;
		case tfxGraphPreset::tfxLifePreset:
			frame = 5;
			value = 5;
			break;
		case tfxGraphPreset::tfxAnglePreset:
		case tfxGraphPreset::tfxArcPreset:
		case tfxGraphPreset::tfxEmissionRangePreset:
			frame = 5;
			value = 0.1f;
			break;
		case tfxGraphPreset::tfxDimensionsPreset:
		case tfxGraphPreset::tfxAmountPreset:
		case tfxGraphPreset::tfxVelocityPreset:
		case tfxGraphPreset::tfxWeightPreset:
		case tfxGraphPreset::tfxWeightVariationPreset:
		case tfxGraphPreset::tfxSpinPreset:
		case tfxGraphPreset::tfxSpinVariationPreset:
		case tfxGraphPreset::tfxFrameratePreset:
			frame = 5.f;
			value = 1.f;
			break;
		case tfxGraphPreset::tfxVelocityTurbulancePreset:
			frame = .001f;
			value = .01f;
			break;
		case tfxGraphPreset::tfxWeightOvertimePreset:
		case tfxGraphPreset::tfxVelocityOvertimePreset:
		case tfxGraphPreset::tfxSpinOvertimePreset:
		case tfxGraphPreset::tfxDirectionVariationPreset:
			frame = 0.001f;
			value = 0.01f;
			break;
		case tfxGraphPreset::tfxColorPreset:
			frame = 0.001f;
			value = 1.f;
			break;
		case tfxGraphPreset::tfxPercentOvertime:
			frame = 0.05f;
			value = 0.05f;
			break;
		default:
			frame = 1;
			value = 0.1f;
			break;
		}
	}

	void tfxGraph::Clear() {
		nodes.clear();
	}

	void tfxGraph::Free() {
		//Explicitly free the nodes
		nodes.free_all();
		lookup.values.free();
	}

	void tfxGraph::Copy(tfxGraph &to) {
		to.Clear();
		do {
			for (auto &n : nodes) {
				to.nodes.push_back(n);
			}
		} while (!nodes.EndOfBuckets());
		if(IsOvertimeGraph())
			CompileGraphOvertime(to);
		else
			CompileGraph(to);
	}

	bool tfxGraph::Sort() {
		bool needed_sorting = false;
		for (tfxU32 i = 1; i < nodes.current_size; ++i) {
			tfxAttributeNode key = nodes[i];
			int j = i - 1;
			while (j >= 0 && key.frame < nodes[j].frame) {
				nodes[j + 1] = nodes[j];
				--j;
				needed_sorting = true;
			}
			nodes[j + 1] = key;
		}
		return needed_sorting;
	}

	void tfxGraph::ReIndex() {
		tfxU32 i = 0;
		for (auto &a : nodes) {
			a.index = i++;
		}
	}

	tfxVec2 tfxGraph::GetInitialZoom() {
		switch (graph_preset) {
		case tfxGraphPreset::tfxOpacityOvertimePreset:
			return tfxVec2(.0017f, 0.00275f);
		case tfxGraphPreset::tfxGlobalPercentPreset:
			return tfxVec2(10.f, 0.005f);
			break;
		case tfxGraphPreset::tfxGlobalPercentPresetSigned:
			return tfxVec2(10.f, 0.006f);
			break;
		case tfxGraphPreset::tfxGlobalOpacityPreset:
			return tfxVec2(10.f, 0.003f);
			break;
		case tfxGraphPreset::tfxLifePreset:
			return tfxVec2(10.f, 3.5f);
			break;
		case tfxGraphPreset::tfxAnglePreset:
			return tfxVec2(10.f, 1.f);
			break;
		case tfxGraphPreset::tfxArcPreset:
			return tfxVec2(10.f, 1.f);
			break;
		case tfxGraphPreset::tfxEmissionRangePreset:
			return tfxVec2(10.f, .5f);
			break;
		case tfxGraphPreset::tfxAmountPreset:
			return tfxVec2(10.f, 1.25f);
			break;
		case tfxGraphPreset::tfxFrameratePreset:
			return tfxVec2(0.0017f, .5f);
			break;
		case tfxGraphPreset::tfxVelocityTurbulancePreset:
			return tfxVec2(0.0017f, .1f);
			break;
		case tfxGraphPreset::tfxDimensionsPreset:
		case tfxGraphPreset::tfxVelocityPreset:
		case tfxGraphPreset::tfxWeightPreset:
		case tfxGraphPreset::tfxWeightVariationPreset:
		case tfxGraphPreset::tfxSpinPreset:
		case tfxGraphPreset::tfxSpinVariationPreset:
			return tfxVec2(10.f, 2.5f);
			break;
		case tfxGraphPreset::tfxNoiseResolutionPreset:
			return tfxVec2(10.f, 1.f);
			break;
		case tfxGraphPreset::tfxNoiseOffsetVariationPreset:
			return tfxVec2(10.f, .01f);
			break;
		case tfxGraphPreset::tfxDirectionOvertimePreset:
			return tfxVec2(0.0017f, 1.f);
			break;
		case tfxGraphPreset::tfxWeightOvertimePreset:
		case tfxGraphPreset::tfxVelocityOvertimePreset:
		case tfxGraphPreset::tfxSpinOvertimePreset:
		case tfxGraphPreset::tfxDirectionVariationPreset:
		case tfxGraphPreset::tfxPercentOvertime:
			return tfxVec2(0.0017f, 0.0035f);
			break;
		case tfxGraphPreset::tfxIntensityOvertimePreset:
			return tfxVec2(0.0017f, 0.01115f);
			break;
		case tfxGraphPreset::tfxColorPreset:
			break;
		default:
			return tfxVec2(0.1f, 0.1f);
			break;
		}

		return tfxVec2(0.1f, 0.1f);
	}

	tfxVec2 tfxGraph::GetInitialZoom3d() {
		switch (graph_preset) {
		case tfxGraphPreset::tfxOpacityOvertimePreset:
			return tfxVec2(.0017f, 0.00275f);
		case tfxGraphPreset::tfxGlobalPercentPreset:
			return tfxVec2(10.f, 0.005f);
			break;
		case tfxGraphPreset::tfxGlobalPercentPresetSigned:
			return tfxVec2(10.f, 0.006f);
			break;
		case tfxGraphPreset::tfxGlobalOpacityPreset:
			return tfxVec2(10.f, 0.003f);
			break;
		case tfxGraphPreset::tfxLifePreset:
			return tfxVec2(10.f, 3.5f);
			break;
		case tfxGraphPreset::tfxAnglePreset:
			return tfxVec2(10.f, 1.f);
			break;
		case tfxGraphPreset::tfxArcPreset:
			return tfxVec2(10.f, 1.f);
			break;
		case tfxGraphPreset::tfxEmissionRangePreset:
			return tfxVec2(10.f, .5f);
			break;
		case tfxGraphPreset::tfxAmountPreset:
			return tfxVec2(10.f, 1.25f);
			break;
		case tfxGraphPreset::tfxFrameratePreset:
			return tfxVec2(0.0017f, .01f);
			break;
		case tfxGraphPreset::tfxDimensionsPreset:
		case tfxGraphPreset::tfxVelocityPreset:
		case tfxGraphPreset::tfxWeightPreset:
		case tfxGraphPreset::tfxWeightVariationPreset:
			return tfxVec2(10.f, 0.01f);
			break;
		case tfxGraphPreset::tfxVelocityTurbulancePreset:
			return tfxVec2(0.0017f, .01f);
			break;
		case tfxGraphPreset::tfxSpinPreset:
		case tfxGraphPreset::tfxSpinVariationPreset:
			return tfxVec2(10.f, 2.5f);
			break;
		case tfxGraphPreset::tfxNoiseResolutionPreset:
			return tfxVec2(10.f, .01f);
			break;
		case tfxGraphPreset::tfxNoiseOffsetVariationPreset:
			return tfxVec2(10.f, .01f);
			break;
		case tfxGraphPreset::tfxDirectionOvertimePreset:
			return tfxVec2(0.0017f, 1.f);
			break;
		case tfxGraphPreset::tfxWeightOvertimePreset:
		case tfxGraphPreset::tfxVelocityOvertimePreset:
		case tfxGraphPreset::tfxSpinOvertimePreset:
		case tfxGraphPreset::tfxDirectionVariationPreset:
		case tfxGraphPreset::tfxPercentOvertime:
			return tfxVec2(0.0017f, 0.0035f);
			break;
		case tfxGraphPreset::tfxIntensityOvertimePreset:
			return tfxVec2(0.0017f, 0.01115f);
			break;
		case tfxGraphPreset::tfxColorPreset:
			break;
		default:
			return tfxVec2(0.1f, 0.1f);
			break;
		}

		return tfxVec2(0.1f, 0.1f);
	}

	void CompileGraph(tfxGraph &graph) {
		float last_frame = graph.GetLastFrame();
		graph.lookup.last_frame = tfxU32(last_frame / tfxLOOKUP_FREQUENCY);
		if (graph.lookup.last_frame) {
			assert(graph.lookup.values.resize(graph.lookup.last_frame + 1));
			for (tfxU32 f = 0; f != graph.lookup.last_frame + 1; ++f) {
				graph.lookup.values[f] = graph.GetValue((float)f * tfxLOOKUP_FREQUENCY);
			}
			graph.lookup.values[graph.lookup.last_frame] = graph.GetLastValue();
		}
		else {
			graph.lookup.values.resize(1);
			graph.lookup.values[0] = graph.GetFirstValue();
		}
	}

	void CompileGraphOvertime(tfxGraph &graph) {
		if (graph.nodes.size() > 1) {
			graph.lookup.last_frame = tfxU32(graph.lookup.life / tfxLOOKUP_FREQUENCY_OVERTIME);
			assert(graph.lookup.values.resize(graph.lookup.last_frame + 1));
			for (tfxU32 f = 0; f != graph.lookup.last_frame + 1; ++f) {
				graph.lookup.values[f] = graph.GetValue((float)f * tfxLOOKUP_FREQUENCY_OVERTIME, graph.lookup.life);
			}
			graph.lookup.values[graph.lookup.last_frame] = graph.GetLastValue();
		}
		else {
			graph.lookup.last_frame = 0;
			graph.lookup.values.resize(1);
			graph.lookup.values[0] = graph.GetFirstValue();
		}
	}

	float LookupFastOvertime(tfxGraph &graph, float age, float lifetime) {
		tfxU32 frame = static_cast<tfxU32>((age / lifetime * graph.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);
		frame = std::min<tfxU32>(frame, graph.lookup.last_frame);
		return graph.lookup.values[frame];
	}

	float LookupFast(tfxGraph &graph, float frame) {
		if ((tfxU32)frame < graph.lookup.last_frame)
			return graph.lookup.values[(tfxU32)frame];
		return graph.lookup.values[graph.lookup.last_frame];
	}

	float LookupPreciseOvertime(tfxGraph &graph, float age, float life) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		tfxAttributeNode *lastec = nullptr;
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

	float LookupPrecise(tfxGraph &graph, float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		tfxAttributeNode *lastec = nullptr;
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

	float GetRandomFast(tfxGraph &graph, float frame) {
		float value = 0;
		if ((tfxU32)frame < graph.lookup.last_frame)
			value = graph.lookup.values[(tfxU32)frame];
		value = graph.lookup.values[graph.lookup.last_frame];
		return random_generation.Range(value);
	}

	float GetRandomPrecise(tfxGraph &graph, float frame) {
		return graph.GetRandomValue(frame);
	}

	float GetMaxLife(tfxEffectEmitter &e) {
		tfxGraph &life = *e.GetGraphByType(tfxBase_life);
		tfxGraph &life_variation = *e.GetGraphByType(tfxVariation_life);
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

	float GetMaxAmount(tfxEffectEmitter &e) {
		if (e.common.property_flags & tfxEmitterPropertyFlags_single)
			return (float)e.GetProperties().spawn_amount;
		tfxGraph &amount = *e.GetGraphByType(tfxBase_amount);
		tfxGraph &amount_variation = *e.GetGraphByType(tfxVariation_amount);
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

	bool IsOvertimeGraph(tfxGraphType type) {
		return type >= tfxOvertime_velocity && type != tfxOvertime_velocity_adjuster;
	}

	bool IsOvertimePercentageGraph(tfxGraphType type) {
		return type >= tfxOvertime_velocity && type != tfxOvertime_velocity_adjuster && type != tfxOvertime_direction;
	}

	bool IsGlobalGraph(tfxGraphType type) {
		return type >= tfxGlobal_life && type <= tfxGlobal_effect_yaw;
	}

	bool IsGlobalPercentageGraph(tfxGraphType type) {
		return type >= tfxGlobal_life && type <= tfxGlobal_splatter;
	}

	bool IsAngleGraph(tfxGraphType type) {
		return (type == tfxGlobal_effect_roll || type == tfxGlobal_effect_pitch || type == tfxGlobal_effect_yaw || type == tfxProperty_emission_pitch || type == tfxProperty_emission_yaw || type == tfxProperty_emission_range ||
			type == tfxProperty_emitter_roll || type == tfxProperty_emitter_pitch || type == tfxProperty_emitter_yaw || type == tfxProperty_arc_offset || type == tfxProperty_arc_size || type == tfxBase_spin || type == tfxVariation_spin);
	}

	bool IsAngleOvertimeGraph(tfxGraphType type) {
		return type == tfxOvertime_direction;
	}

	bool IsEverythingElseGraph(tfxGraphType type) {
		return !IsOvertimeGraph(type) && !IsOvertimePercentageGraph(type) && !IsGlobalGraph(type) && !IsAngleGraph(type) && !IsOvertimeGraph(type);
	}

	bool HasDataValue(tfxStorageMap<tfxDataEntry> &config, tfxStr32 key) {
		return config.ValidName(key);
	}

	void AddDataValue(tfxStorageMap<tfxDataEntry> &map, tfxStr32 key, const char *value) {
		tfxDataEntry entry;
		entry.type = tfxString;
		entry.key = key;
		entry.str_value = value;
		map.Insert(key, entry);
	}

	void AddDataValue(tfxStorageMap<tfxDataEntry> &map, tfxStr32 key, int value) {
		tfxDataEntry entry;
		entry.type = tfxSInt;
		entry.key = key;
		entry.int_value = value;
		entry.bool_value = (bool)value;
		map.Insert(key, entry);
	}

	void AddDataValue(tfxStorageMap<tfxDataEntry> &map, tfxStr32 key, bool value) {
		tfxDataEntry entry;
		entry.type = tfxBool;
		entry.key = key;
		entry.bool_value = value;
		entry.int_value = (int)value;
		map.Insert(key, entry);
	}

	void AddDataValue(tfxStorageMap<tfxDataEntry> &map, tfxStr32 key, double value) {
		tfxDataEntry entry;
		entry.type = tfxDouble;
		entry.key = key;
		entry.double_value = value;
		map.Insert(key, entry);
	}

	void AddDataValue(tfxStorageMap<tfxDataEntry> &map, tfxStr32 key, float value) {
		tfxDataEntry entry;
		entry.type = tfxFloat;
		entry.key = key;
		entry.float_value = value;
		map.Insert(key, entry);
	}

	tfxStr &GetDataStrValue(tfxStorageMap<tfxDataEntry> &map, const char* key) {
		return map.At(key).str_value;
	}
	int& GetDataIntValue(tfxStorageMap<tfxDataEntry> &map, const char* key) {
		return map.At(key).int_value;
	}
	float& GetDataFloatValue(tfxStorageMap<tfxDataEntry> &map, const char* key) {
		return map.At(key).float_value;
	}

	bool SaveDataFile(tfxStorageMap<tfxDataEntry> &map, const char* path) {
		FILE *file = fopen(path, "w");

		if (file == NULL)
			return false;

		if (map.Size()) {
			for (auto &entry : map.data) {
				tfxStr512 ini_line = entry.key;
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

	bool LoadDataFile(tfxStorageMap<tfxDataEntry> &map, const char* path) {
		if (!data_types.initialised) data_types.Init();
		FILE* fp;
		fp = fopen(path, "r");
		if (fp == NULL) {
			return false;
		}

		const size_t max_line_length = 512;
		char buffer[max_line_length];

		tmpStack(tfxStr64, pair);
		while (fgets(buffer, max_line_length, fp)) {
			buffer[strcspn(buffer, "\n")] = 0;
			tfxStr512 str = buffer;
			pair.clear();
			SplitStringStack(str, pair, 61);
			if (pair.size() == 2) {
				tfxStr64 key = pair[0];
				if (data_types.names_and_types.ValidName(pair[0])) {
					tfxDataType t = data_types.names_and_types.At(pair[0]);
					if (t == tfxBool) {
						AddDataValue(map, key, (bool)atoi(pair[1].c_str()));
					}
					if (t == tfxSInt) {
						AddDataValue(map, key, atoi(pair[1].c_str()));
					}
					else if (t == tfxFloat) {
						AddDataValue(map, key, (float)atof(pair[1].c_str()));
					}
					else if (t == tfxString) {
						AddDataValue(map, key, pair[1].c_str());
					}
				}
			}
		}

		int close = fclose(fp);
		return true;

	}

	void SplitStringStack(const tfxStr &str, tfxStack<tfxStr64> &pair, char delim) {
		tfxStr64 line;
		for (char c : str) {
			if (c == delim && line.Length() && c != NULL) {
				pair.push_back(line);
				line.Clear();
			}
			else if (c != NULL) {
				line.Append(c);
			}
		}

		if (line.Length()) {
			pair.push_back(line);
		}
	}

	bool StringIsUInt(const tfxStr &s) {

		for (auto c : s) {
			if (!std::isdigit(c) && c != 0)
				return false;
		}

		return true;
	}

	int GetDataType(const tfxStr &s) {
		if (s.Length() == 0)
			return tfxString;

		if(s.IsInt())
			return tfxSInt;

		if (s.IsFloat())
			return tfxFloat;

		return tfxString;
	}

	//Get a graph by tfxGraphID
	tfxGraph &GetGraph(tfxEffectLibrary &library, tfxGraphID &graph_id) {
		tfxGraphType type = graph_id.type;

		if (type < tfxGlobalCount) {
			return ((tfxGraph*)&library.global_graphs[graph_id.graph_id])[type];
		}
		else if (type >= tfxPropertyStart && type < tfxBaseStart) {
			int ref = type - tfxPropertyStart;
			return ((tfxGraph*)&library.emitter_attributes[graph_id.graph_id].properties)[ref];
		}
		else if (type >= tfxBaseStart && type < tfxVariationStart) {
			int ref = type - tfxBaseStart;
			return ((tfxGraph*)&library.emitter_attributes[graph_id.graph_id].base)[ref];
		}
		else if (type >= tfxVariationStart && type < tfxOvertimeStart) {
			int ref = type - tfxVariationStart;
			return ((tfxGraph*)&library.emitter_attributes[graph_id.graph_id].variation)[ref];
		}
		else if (type >= tfxOvertimeStart) {
			int ref = type - tfxOvertimeStart;
			return ((tfxGraph*)&library.emitter_attributes[graph_id.graph_id].overtime)[ref];
		}

		assert(0);	//This function must return a value, make sure the graph_id is valid

		return((tfxGraph*)&library.emitter_attributes[graph_id.graph_id].overtime)[type];

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
		tmpStack(tfxStr64, pair);

		while (!data->data.EoF()) {
			pair.clear();
			tfxStr128 line = data->data.ReadLine();
			bool context_set = false;
			if (StringIsUInt(line.c_str())) {
				context = atoi(line.c_str());
				if (context == tfxEndShapes)
					break;
				context_set = true;
			}
			if (context_set == false) {
				SplitStringStack(line.c_str(), pair);
				if (pair.size() != 2) {
					pair.clear();
					SplitStringStack(line.c_str(), pair, 44);
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

	int GetEffectLibraryStats(const char *filename, tfxEffectLibraryStats &stats) {
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

		memset(&stats, 0, sizeof(tfxEffectLibraryStats));
		bool inside_emitter = false;

		tmpStack(tfxStr64, pair);
		while (!data->data.EoF()) {
			pair.clear();
			tfxStr128 line = data->data.ReadLine();
			bool context_set = false;
			if (StringIsUInt(line.c_str())) {
				context_set = true;
				if (context == tfxEndEmitter) {
					inside_emitter = false;
				}
			}
			if (context_set == false) {
				SplitStringStack(line.c_str(), pair);
				if (pair.size() != 2) {
					pair.clear();
					SplitStringStack(line.c_str(), pair, 44);
					if (pair.size() < 2) {
						error = 1;
						break;
					}
				}
				if (context == tfxStartShapes) {
					if (pair.size() >= 5) {
						int frame_count = atoi(pair[2].c_str());
						stats.total_shapes += frame_count;
					}
				}
				else if (context == tfxStartEmitter) {
					inside_emitter = true;
					stats.total_emitters++;
				}
				else if (context == tfxStartEffect) {
					if (inside_emitter)
						stats.total_sub_effects++;
					else
						stats.total_effects++;
				}
			}
		}

		return error;
	}

	tfxEffectLibraryStats CreateLibraryStats(tfxEffectLibrary &lib) {
		tfxEffectLibraryStats stats;
		memset(&stats, 0, sizeof(stats));
		stats.total_effects = lib.effects.size();
		stats.total_node_lookup_indexes = lib.node_lookup_indexes.size();
		stats.total_attribute_nodes = lib.all_nodes.size();
		tmpStack(tfxEffectEmitter, stack);
		for (auto &effect : lib.effects) {
			stack.push_back(effect);
		}
		while (!stack.empty()) {
			tfxEffectEmitter &current = stack.pop_back();
			if (current.parent) {
				if (current.type == tfxEffectType) {
					stats.total_sub_effects++;
				}
				else if(current.type == tfxEmitterType) {
					stats.total_emitters++;
				}
			}
			for (auto &sub : current.GetInfo().sub_effectors) {
				stack.push_back(sub);
			}
		}
		stats.total_shapes = lib.particle_shapes.data.size();
		stats.required_graph_node_memory = lib.graph_node_allocator.TotalMemoryInUse();
		stats.required_graph_lookup_memory = lib.graph_lookup_allocator.TotalMemoryInUse();

		return stats;
	}

	tfxErrorFlags LoadEffectLibraryPackage(const char *filename, tfxEffectLibrary &lib, void(*shape_loader)(const char* filename, tfxImageData &image_data, void *raw_image_data, int image_size, void *user_data), void *user_data, bool read_only) {
		assert(shape_loader);
		if (!data_types.initialised) data_types.Init();
		lib.Clear();

		tmpStack(tfxEffectEmitter, effect_stack);
		int context = 0;
		tfxErrorFlags error = 0;
		int uid = 0;
		tfxU32 current_global_graph = 0;

		tfxPackage package;
		error = LoadPackage(filename, package);

		tfxEntryInfo *data = package.GetFile("data.txt");
		tfxEntryInfo *stats_struct = package.GetFile("stats.struct");

		if (!stats_struct) {
			lib.graph_node_allocator = CreateArenaManager(tfxMegabyte(2), 8);
			lib.graph_lookup_allocator = CreateArenaManager(tfxMegabyte(4), 256);
		}
		else {
			tfxEffectLibraryStats stats;
			memcpy(&stats, stats_struct->data.data, stats_struct->file_size);
			if (read_only) {
				lib.graph_node_allocator = CreateArenaManager(NearestMultiple((size_t)stats.required_graph_node_memory, tfxMegabyte(2)), 8);
				lib.graph_lookup_allocator = CreateArenaManager(NearestMultiple((size_t)stats.required_graph_lookup_memory, tfxMegabyte(4)), 256);
			}
			else {
				lib.graph_node_allocator = CreateArenaManager(tfxMegabyte(2), 8);
				lib.graph_lookup_allocator = CreateArenaManager(tfxMegabyte(4), 256);
			}
		}

		if (!data)
			error |= tfxErrorCode_data_could_not_be_loaded;

		if (error != 0) {
			package.Free();
			return error;
		}

		int first_shape_index = -1;
		tmpStack(tfxStr64, pair);

		while (!data->data.EoF()) {
			tfxStr512 line = data->data.ReadLine();
			bool context_set = false;

			if (StringIsUInt(line.c_str())) {
				context = atoi(line.c_str());
				if (context == tfxEndOfFile)
					break;

				context_set = true;
				if (context == tfxStartFolder) {
					tfxEffectEmitter effect;
					effect.common.library = &lib;
					effect.type = tfxEffectEmitterType::tfxFolder;
					effect.info_index = lib.AddEffectEmitterInfo();
					effect.property_index = lib.AddEmitterProperties();
					lib.GetInfo(effect).uid = uid++;
					effect_stack.push_back(effect);
				}
				if (context == tfxStartEffect) {
					tfxEffectEmitter effect;
					effect.common.library = &lib;
					effect.info_index = lib.AddEffectEmitterInfo();
					effect.property_index = lib.AddEmitterProperties();
					if (effect_stack.size() <= 1) { //Only root effects get the global graphs
						lib.AddEffectGraphs(effect);
						effect.ResetEffectGraphs(false, false);
						current_global_graph = effect.global;
					}
					effect.type = tfxEffectEmitterType::tfxEffectType;
					lib.AddAnimationSettings(effect);
					lib.AddPreviewCameraSettings(effect);
					lib.GetInfo(effect).uid = uid++;
					effect_stack.push_back(effect);

				}
				if (context == tfxStartEmitter) {
					tfxEffectEmitter emitter;
					emitter.common.library = &lib;
					emitter.info_index = lib.AddEffectEmitterInfo();
					emitter.property_index = lib.AddEmitterProperties();
					lib.AddEmitterGraphs(emitter);
					emitter.type = tfxEffectEmitterType::tfxEmitterType;
					emitter.ResetEmitterGraphs(false, false);
					lib.GetInfo(emitter).uid = uid++;
					effect_stack.push_back(emitter);
				}
			}

			if (context_set == false) {
				pair.clear();
				SplitStringStack(line.c_str(), pair);
				if (pair.size() != 2) {
					pair.clear();
					SplitStringStack(line.c_str(), pair, 44);
					if (pair.size() < 2) {
						error |= tfxErrorCode_some_data_not_loaded;
						break;
					}
				}

				if (context == tfxStartAnimationSettings || context == tfxStartEmitter || context == tfxStartEffect || context == tfxStartFolder || context == tfxStartPreviewCameraSettings) {
					if (data_types.names_and_types.ValidName(pair[0])) {
						switch (data_types.names_and_types.At(pair[0])) {
						case tfxUint:
							AssignEffectorProperty(effect_stack.back(), pair[0], (tfxU32)atoi(pair[1].c_str()));
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
					else {
						error |= tfxErrorCode_some_data_not_loaded;
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
						tfxShapeData s;
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
							tfxImageData image_data;
							image_data.animation_frames = (float)s.frame_count;
							image_data.image_size = tfxVec2((float)s.width, (float)s.height);
							image_data.shape_index = s.shape_index;
							image_data.import_filter = s.import_filter;

							shape_loader(s.name, image_data, shape_entry->data.data, (tfxU32)shape_entry->file_size, user_data);

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
#ifdef tfxTRACK_MEMORY
				tfxvec<tfxEffectEmitter> &sub_effectors = effect_stack.parent().GetInfo().sub_effectors;
				memcpy(sub_effectors.name, "emitter_sub_effects\0", 20);
#endif
				effect_stack.parent().GetInfo().sub_effectors.push_back(effect_stack.back());
				effect_stack.pop();
			}

			if (context == tfxEndEffect) {
				effect_stack.back().ReIndex();
				if (effect_stack.size() > 1) {
#ifdef tfxTRACK_MEMORY
					tfxvec<tfxEffectEmitter> &sub_effectors = effect_stack.parent().GetInfo().sub_effectors;
					memcpy(sub_effectors.name, "effect_sub_emitters\0", 20);
#endif
					effect_stack.parent().GetInfo().sub_effectors.push_back(effect_stack.back());
					effect_stack.back().InitialiseUninitialisedGraphs();
				}
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

		if (uid >= 0) {
			//Effects were loaded so let's compile them
			lib.CompileAllGraphs();
			lib.ReIndex();
			if(first_shape_index != -1)
				lib.UpdateParticleShapeReferences(lib.effects, first_shape_index);
			lib.UpdateEffectPaths();
			lib.UpdateComputeNodes();
			//lib.UpdateEffectParticleStorage();
			lib.SetMinMaxData();
		}
		lib.uid = uid;

		return error;

	}

	void SetTimeOut(tfxEffectTemplate &effect_template, float frames) {
		if (frames < 1.f) frames = 1.f;
		tmpStack(tfxEffectEmitter*, stack);
		stack.push_back(&effect_template.effect);
		while (stack.size()) {
			tfxEffectEmitter *current = stack.pop_back();
			current->common.timeout = frames;
			for (auto &sub : current->GetInfo().sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	void tfxEffectTemplate::SetUserDataAll(void *data) {
		tmpStack(tfxEffectEmitter*, stack);
		stack.push_back(&effect);
		while (stack.size()) {
			tfxEffectEmitter *current = stack.pop_back();
			current->user_data = data;
			for (auto &sub : current->GetInfo().sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	void tfxEffectTemplate::SetParticleUpdateCallback(tfxStr256 path, void(*particle_update_callback)(tfxParticleData &particle, void *user_data)) {
		assert(paths.ValidName(path));
		tfxEffectEmitter &e = *paths.At(path);
		assert(e.type == tfxEmitterType);
		e.particle_update_callback = particle_update_callback;
	}

	void tfxEffectTemplate::SetParticleOnSpawnCallback(tfxStr256 path, void(*particle_onspawn_callback)(tfxParticle &particle, void *user_data)) {
		assert(paths.ValidName(path));
		tfxEffectEmitter &e = *paths.At(path);
		assert(e.type == tfxEmitterType);
		e.particle_onspawn_callback = particle_onspawn_callback;
	}

	tfxParticleManager::~tfxParticleManager() {
	}

	tfxEffectEmitter &tfxParticleManager::operator[] (unsigned int index) {
		return effects[current_ebuff][index];
	}

	void tfxParticleManager::AddEffect(tfxEffectEmitter &effect, unsigned int buffer, bool is_sub_emitter) {
		tfxPROFILE;
		if (effects[buffer].current_size == effects[buffer].capacity)
			return;
		if (flags & tfxEffectManagerFlags_use_compute_shader && highest_compute_controller_index >= max_compute_controllers && free_compute_controllers.empty())
			return;
		unsigned int parent_index = effects[buffer].current_size++;
		effects[buffer][parent_index] = effect;
		if(!is_sub_emitter)
			effects[buffer][parent_index].parent_particle = nullptr;
		effects[buffer][parent_index].flags &= ~tfxEmitterStateFlags_retain_matrix;
		effects[buffer][parent_index].ResetParents();
		if (!effect.Is3DEffect()) {
			flags &= ~tfxEffectManagerFlags_3d_effects;
		}
		for (auto &e : effect.GetInfo().sub_effectors) {
			if (!FreeEffectCapacity())
				break;
			if (e.flags & tfxEmitterStateFlags_enabled) {
				unsigned int index = effects[buffer].current_size++;
				effects[buffer][index] = e;
				tfxEffectEmitter &emitter = effects[buffer].back();
				tfxEmitterProperties &properties = emitter.GetProperties();
				emitter.parent = &effects[buffer][parent_index];
				emitter.next_ptr = &emitter;
				if (emitter.particles_index == tfxINVALID && flags & tfxEffectManagerFlags_unorderd) {
					emitter.particles_index = GrabParticleBank(*this, emitter.path_hash, 100);
				}
				else {
					emitter.particles_index = properties.layer * 2 + current_pbuff;
				}
				emitter.flags &= ~tfxEmitterStateFlags_retain_matrix;
				emitter.flags |= emitter.parent->flags & tfxEmitterStateFlags_no_tween;
				emitter.flags |= e.common.property_flags & tfxEmitterPropertyFlags_single && !(flags & tfxEffectManagerFlags_disable_spawning) ? tfxEmitterStateFlags_is_single : 0;
				emitter.flags |= emitter.common.property_flags & tfxParticleControlFlags_base_uniform_size;
				emitter.flags |= (properties.emission_type != tfxLine && !(e.common.property_flags & tfxEmitterPropertyFlags_edge_traversal)) || properties.emission_type == tfxLine && !(e.common.property_flags & tfxEmitterPropertyFlags_edge_traversal) ? tfxEmitterStateFlags_not_line : 0;
				emitter.flags |= e.common.property_flags & tfxEmitterPropertyFlags_random_color;
				emitter.flags |= e.common.property_flags & tfxEmitterPropertyFlags_lifetime_uniform_size;
				emitter.flags |= properties.angle_settings != tfxAngleSettingFlags_align_roll && !(e.common.property_flags & tfxEmitterPropertyFlags_relative_angle) ? tfxEmitterStateFlags_can_spin : 0;
				emitter.flags |= properties.angle_settings == tfxAngleSettingFlags_align_roll ? tfxEmitterStateFlags_align_with_velocity : 0;
				emitter.flags |= properties.emission_type == tfxLine && e.common.property_flags & tfxEmitterPropertyFlags_edge_traversal ? tfxEmitterStateFlags_is_line_traversal : 0;
				emitter.flags |= e.common.property_flags & tfxEmitterPropertyFlags_play_once;
				emitter.flags |= properties.end_behaviour == tfxLoop ? tfxEmitterStateFlags_loop : 0;
				emitter.flags |= properties.end_behaviour == tfxKill ? tfxEmitterStateFlags_kill : 0;
				if (is_sub_emitter)
					emitter.flags |= tfxEmitterStateFlags_is_sub_emitter;

				if (effect.Is3DEffect()) {
					if (e.common.property_flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type == tfxLine) {
						emitter.current.transform_particle_callback3d = TransformParticle3dRelativeLine;
					}
					else if (e.common.property_flags & tfxEmitterPropertyFlags_relative_position) {
						emitter.current.transform_particle_callback3d = TransformParticle3dRelative;
					}
					else if (e.common.property_flags & tfxEmitterPropertyFlags_relative_angle) {
						emitter.current.transform_particle_callback3d = TransformParticle3dAngle;
					}
					else {
						emitter.current.transform_particle_callback3d = TransformParticle3d;
					}
				}
				else {
					if (e.common.property_flags & tfxEmitterPropertyFlags_relative_position && (e.common.property_flags & tfxEmitterPropertyFlags_relative_angle || (e.common.property_flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type == tfxLine))) {
						emitter.current.transform_particle_callback2d = TransformParticleRelativeLine;
					}
					else if (e.common.property_flags & tfxEmitterPropertyFlags_relative_position) {
						emitter.current.transform_particle_callback2d = TransformParticleRelative;
					}
					else if (e.common.property_flags & tfxEmitterPropertyFlags_relative_angle) {
						emitter.current.transform_particle_callback2d = TransformParticleAngle;
					}
					else {
						emitter.current.transform_particle_callback2d = TransformParticle;
					}
				}

				if (flags & tfxEffectManagerFlags_use_compute_shader && e.GetInfo().sub_effectors.empty()) {
					int free_slot = AddComputeController();
					if (free_slot != -1) {
						emitter.compute_slot_id = free_slot;
						emitter.common.property_flags |= tfxEmitterPropertyFlags_is_bottom_emitter;
						tfxComputeController &c = *(static_cast<tfxComputeController*>(compute_controller_ptr) + free_slot);
						c.flags = 0;
						c.flags |= emitter.common.property_flags & tfxParticleControlFlags_random_color;
						c.flags |= emitter.common.property_flags & tfxParticleControlFlags_relative_position;
						c.flags |= emitter.common.property_flags & tfxParticleControlFlags_relative_angle;
						c.flags |= emitter.common.property_flags & tfxParticleControlFlags_edge_traversal;
						c.flags |= emitter.common.property_flags & tfxParticleControlFlags_base_uniform_size;
						c.flags |= emitter.common.property_flags & tfxParticleControlFlags_lifetime_uniform_size;
						c.flags |= emitter.common.property_flags & tfxParticleControlFlags_reverse_animation;
						c.flags |= emitter.common.property_flags & tfxParticleControlFlags_play_once;
						c.flags |= ((properties.emission_type == tfxPoint) << 3);
						c.flags |= ((properties.emission_type == tfxArea) << 4);
						c.flags |= ((properties.emission_type == tfxLine) << 5);
						c.flags |= ((properties.emission_type == tfxEllipse) << 6);
						c.flags |= ((properties.end_behaviour == tfxLoop) << 7);
						c.flags |= ((properties.end_behaviour == tfxKill) << 8);
						c.flags |= ((properties.end_behaviour == tfxLetFree) << 9);
						properties.compute_flags = c.flags;
						if (emitter.common.property_flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
							c.image_handle = tfxVec2(0.5f, 0.5f);
						}
						else {
							c.image_handle = properties.image_handle;
						}
						c.image_handle = properties.image_handle;
					}
				}
			}
		}
		effects[buffer][parent_index].NoTweenNextUpdate(); 
	}

	void tfxParticleManager::AddEffect(tfxEffectTemplate &effect) {
		AddEffect(effect.effect, current_ebuff);
	}

	int tfxParticleManager::AddComputeController() {
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

	void tfxParticleManager::ResetParticlePtr(void *ptr) {
		new_compute_particle_ptr = ptr;
		new_compute_particle_index = 0;
	}

	void tfxParticleManager::ResetControllerPtr(void *ptr) {
		compute_controller_ptr = ptr;
	}

	void tfxParticleManager::UpdateCompute(void *sampled_particles, unsigned int sample_size) {
		for (int i = 0; i != sample_size; ++i) {
			if (compute_global_state.current_length == 0)
				break;
			tfxComputeParticle *sample = static_cast<tfxComputeParticle*>(sampled_particles) + i;
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

	tfxParticle& tfxParticleManager::GrabCPUParticle(unsigned int index) {
		return particle_banks[index].grab();
	}

	tfxComputeParticle& tfxParticleManager::GrabComputeParticle(unsigned int layer) {
		assert(new_compute_particle_ptr);		//Use must assign the compute ptr to point to an area in memory where you can stage new particles for uploading to the GPU - See ResetComputePtr
		return *(static_cast<tfxComputeParticle*>(new_compute_particle_ptr) + new_compute_particle_index++);
	}

	void tfxParticleManager::FreeParticleBank(tfxEffectEmitter &emitter) {
		if (free_particle_banks.ValidKey(emitter.path_hash)) {
			free_particle_banks.At(emitter.path_hash).push_back(emitter.particles_index);
		}
		else {
			tfxvec<tfxU32> new_indexes(tfxCONSTRUCTOR_VEC_INIT("FreeParticleBanks::new_indexes"));
			new_indexes.push_back(emitter.particles_index);
			free_particle_banks.Insert(emitter.path_hash, new_indexes);
		}
	}

	void tfxParticleManager::Update() {
		tfxPROFILE;
		new_compute_particle_index = 0;

		unsigned int next_buffer = !current_ebuff;
		effects[next_buffer].clear();

		memset(new_particles_index_start, tfxMAX_UINT, 4 * tfxLAYERS);
		memset(sprite_index_point, 0, 4 * tfxLAYERS);

		if (!(flags & tfxEffectManagerFlags_unorderd)) {
			for (tfxEachLayer) {
				int layer_offset = layer * 2;
				int current_buffer_index = current_pbuff + layer_offset;
				sprite_index_point[layer] = particle_banks[current_buffer_index].current_size;
			}
		}

		for (tfxEachLayer) {
			sprites2d[layer].clear();
			sprites3d[layer].clear();
		}

		tfxU32 start_size = effects[current_ebuff].current_size;
		for (auto &e : effects[current_ebuff]) {
			UpdatePMEmitter(*this, e);
			if (e.type == tfxEffectType) {
				if (e.common.timeout_counter <= e.common.timeout) {
					e.next_ptr = SetNextEffect(e, next_buffer);
				}
				else {
					e.next_ptr = nullptr;
				}
			}
			else {
				if (e.common.timeout_counter <= e.common.timeout) {
					e.next_ptr = SetNextEffect(e, next_buffer);
				}
				else {
					e.next_ptr = nullptr;
					if (flags & tfxEffectManagerFlags_use_compute_shader && e.common.property_flags & tfxEmitterPropertyFlags_is_bottom_emitter)
						FreeComputeSlot(e.compute_slot_id);
					if (flags & tfxEffectManagerFlags_unorderd) {
						FreeParticleBank(e);
					}
				}
			}
		}

		for (int i = start_size; i != effects[current_ebuff].current_size; ++i) {
			tfxEffectEmitter &e = effects[current_ebuff][i];
			UpdatePMEmitter(*this, e);
			e.next_ptr = SetNextEffect(e, next_buffer);
		}

		if (!(flags & tfxEffectManagerFlags_unorderd)) {
			for (tfxEachLayer) {
				int layer_offset = layer * 2;
				int a_buffer = !current_pbuff + layer_offset;
				int b_buffer = current_pbuff + layer_offset;
				int a_cap = particle_banks[a_buffer].capacity;
				int b_cap = particle_banks[b_buffer].capacity;
				if (particle_banks[a_buffer].capacity != particle_banks[b_buffer].capacity) {
					if(particle_banks[a_buffer].capacity > particle_banks[b_buffer].capacity) {
						particle_banks[b_buffer].reserve(particle_banks[a_buffer].capacity);
					}
					else {
						particle_banks[a_buffer].reserve(particle_banks[b_buffer].capacity);
					}
				}
			}
		}

		current_ebuff = next_buffer;

		if (!(flags & tfxEffectManagerFlags_unorderd)) {
			if(flags & tfxEffectManagerFlags_3d_effects && tfxEffectManagerFlags_order_by_depth) 
				ControlParticlesDepthOrdered3d(*this);
			else if(flags & tfxEffectManagerFlags_3d_effects)
				ControlParticlesOrdered3d(*this);
			else
				ControlParticlesOrdered2d(*this);
		}

		flags &= ~tfxEffectManagerFlags_update_base_values;

	}

	void ControlParticlesOrdered2d(tfxParticleManager &pm) {
		tfxU32 next_buffer = !pm.current_pbuff;

		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			int layer_offset = layer * 2;
			int next_buffer_index = next_buffer + layer_offset;
			int current_buffer_index = pm.current_pbuff + layer_offset;

			pm.particle_banks[next_buffer_index].clear();

			tfxU32 index = 0;
			for (tfxU32 i = 0; i != pm.particle_banks[current_buffer_index].current_size; ++i) {
				tfxParticle &p = pm.particle_banks[current_buffer_index][i];
				p.parent = p.parent->next_ptr;
				tfxEmitterProperties &properties = p.parent->GetProperties();

				//All of the new particles will be at the end of the buffer so might not need this condition
				//and instead have 2 separate for loops? It is highly predictable though.
				if (!(p.data.flags & tfxParticleFlags_fresh)) {

					tfxParticleSprite2d &s = pm.sprites2d[properties.layer][index];

					p.sprite_index = (properties.layer << 28) + index;
					if (p.parent && ControlParticle(pm, p, s.transform.scale, *p.parent)) {
						if (p.data.flags & tfxParticleFlags_capture_after_transform) {
							p.parent->current.transform_particle_callback2d(p.data, s.transform.position, s.transform.rotation, p.parent->common, p.parent->common.transform.captured_position);
							p.data.captured_position = s.transform.position;
							p.parent->current.transform_particle_callback2d(p.data, s.transform.position, s.transform.rotation, p.parent->common, p.parent->common.transform.world_position);
							p.data.flags &= ~tfxParticleFlags_capture_after_transform;
						}
						else {
							p.parent->current.transform_particle_callback2d(p.data, s.transform.position, s.transform.rotation, p.parent->common, p.parent->common.transform.world_position);
						}
						p.next_ptr = pm.SetNextParticle(next_buffer_index, p);

						s.color = p.data.color;
						s.image_ptr = properties.image->ptr;
						s.image_frame = (tfxU32)p.data.image_frame;
						s.transform.captured_position = p.data.captured_position.xy();
						s.intensity = p.data.intensity;
						s.handle = p.parent->current.image_handle;
						p.next_ptr->data.captured_position = s.transform.position;

					}
					else {
						s.transform.scale = 0;
						s.intensity = 0;
						s.image_frame = 0;
						s.image_ptr = properties.image->ptr;
						p.next_ptr = nullptr;
					}

				}
				else {
					p.sprite_index = (properties.layer << 28) + index;
					p.data.flags &= ~tfxParticleFlags_fresh;
					p.next_ptr = pm.SetNextParticle(next_buffer_index, p);
				}

				index++;
			}
		}

		pm.current_pbuff = next_buffer;

		pm.flags &= ~tfxEffectManagerFlags_update_base_values;
	}

	void ControlParticlesOrdered3d(tfxParticleManager &pm) {
		tfxU32 next_buffer = !pm.current_pbuff;

		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			int layer_offset = layer * 2;
			int next_buffer_index = next_buffer + layer_offset;
			int current_buffer_index = pm.current_pbuff + layer_offset;

			pm.particle_banks[next_buffer_index].clear();

			tfxU32 index = 0;
			for (tfxU32 i = 0; i != pm.particle_banks[current_buffer_index].current_size; ++i) {
				tfxParticle &p = pm.particle_banks[current_buffer_index][i];
				p.parent = p.parent->next_ptr;
				tfxEmitterProperties &properties = p.parent->GetProperties();

				//All of the new particles will be at the end of the buffer so might not need this condition
				//and instead have 2 separate for loops? It is highly predictable though.
				if (!(p.data.flags & tfxParticleFlags_fresh)) {

					tfxParticleSprite3d &s = pm.sprites3d[properties.layer][index];
					p.sprite_index = (properties.layer << 28) + index;

					if (p.parent && ControlParticle(pm, p, s.transform.scale, *p.parent)) {
						if (p.data.flags & tfxParticleFlags_capture_after_transform) {
							p.parent->current.transform_particle_callback3d(p.data, s.transform.position, s.transform.rotations, p.parent->common, p.parent->common.transform.captured_position);
							p.data.captured_position = s.transform.position;
							p.parent->current.transform_particle_callback3d(p.data, s.transform.position, s.transform.rotations, p.parent->common, p.parent->common.transform.world_position);
							p.data.flags &= ~tfxParticleFlags_capture_after_transform;
						}
						else {
							p.parent->current.transform_particle_callback3d(p.data, s.transform.position, s.transform.rotations, p.parent->common, p.parent->common.transform.world_position);
						}

						tfxVec3 alignment_vector = SetParticleAlignment(p, s.transform.position, properties);
						p.next_ptr = pm.SetNextParticle(next_buffer_index, p);

						s.color = p.data.color;
						s.image_frame_plus = (properties.billboard_option << 24) + (tfxU32)p.data.image_frame;
						s.transform.captured_position = p.data.captured_position;
						s.stretch = p.data.velocity_normal.w;
						s.intensity = p.data.intensity;
						s.handle = p.parent->current.image_handle;
						alignment_vector.y += 0.002f;	//We don't want a 0 alignment normal
						s.alignment = Pack10bit(alignment_vector, properties.billboard_option & 0x00000003);
						s.image_ptr = properties.image->ptr;

						p.next_ptr->data.captured_position = s.transform.position;

					}
					else {
						s.transform.scale = 0;
						s.intensity = 0;
						s.image_frame_plus = 0;
						s.image_ptr = properties.image->ptr;
						p.next_ptr = nullptr;
					}

				}
				else {
					p.sprite_index = (properties.layer << 28) + index;
					p.data.flags &= ~tfxParticleFlags_fresh;
					p.next_ptr = pm.SetNextParticle(next_buffer_index, p);
				}

				index++;
			}
		}

		pm.current_pbuff = next_buffer;

		pm.flags &= ~tfxEffectManagerFlags_update_base_values;
	}


	void ControlParticlesDepthOrdered3d(tfxParticleManager &pm) {
		tfxU32 next_buffer = !pm.current_pbuff;

		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			int layer_offset = layer * 2;
			int next_buffer_index = next_buffer + layer_offset;
			int current_buffer_index = pm.current_pbuff + layer_offset;
			pm.particle_banks[next_buffer_index].clear();
			tfxU32 new_index = pm.new_particles_index_start[layer];
			tfxParticle *new_particle = nullptr;
			if (new_index < pm.particle_banks[current_buffer_index].size())
				new_particle = &pm.particle_banks[current_buffer_index][new_index];

			tfxU32 index = 0;
			tfxU32 sprite_index = 0;
			tfxU32 next_index = 0;
			for (tfxU32 i = 0; i != pm.particle_banks[current_buffer_index].current_size; ++i) {
				tfxParticle &p = pm.particle_banks[current_buffer_index][i];
				p.parent = p.parent->next_ptr;
				p.prev_index = index;
				tfxEmitterProperties &properties = p.parent->GetProperties();

				if (!(p.data.flags & tfxParticleFlags_fresh)) {

					tfxParticleSprite3d s;

					if (p.parent && ControlParticle(pm, p, s.transform.scale, *p.parent)) {
						if (p.data.flags & tfxParticleFlags_capture_after_transform) {
							p.parent->current.transform_particle_callback3d(p.data, s.transform.position, s.transform.rotations, p.parent->common, p.parent->common.transform.captured_position);
							p.data.captured_position = s.transform.position;
							p.parent->current.transform_particle_callback3d(p.data, s.transform.position, s.transform.rotations, p.parent->common, p.parent->common.transform.world_position);
							p.data.flags &= ~tfxParticleFlags_capture_after_transform;
						}
						else {
							p.parent->current.transform_particle_callback3d(p.data, s.transform.position, s.transform.rotations, p.parent->common, p.parent->common.transform.world_position);
						}
						p.data.depth = LengthVec3NoSqR(s.transform.position - pm.camera_position);

						tfxVec3 alignment_vector = SetParticleAlignment(p, s.transform.position, properties);
						while (new_particle && new_particle->data.depth > p.data.depth) {
							tfxEmitterProperties &new_properties = new_particle->parent->GetProperties();
							tfxParticleSprite3d &ns = pm.sprites3d[new_properties.layer][new_particle->sprite_index & 0x0FFFFFFF];
							pm.sprites3d[new_properties.layer][sprite_index++] = ns;
							new_particle->sprite_index = (properties.layer << 28) + next_index;
							new_particle->prev_index = new_index;
							new_particle->data.flags &= ~tfxParticleFlags_fresh;
							new_particle->parent = new_particle->parent->next_ptr;
							new_particle->next_ptr = pm.SetNextParticle(next_buffer_index, *new_particle);
							next_index++;
							new_particle = nullptr;
							if (++new_index < pm.particle_banks[current_buffer_index].size())
								new_particle = &pm.particle_banks[current_buffer_index][new_index];
						}

						p.sprite_index = (properties.layer << 28) + sprite_index;
						s.color = p.data.color;
						s.image_frame_plus = (properties.billboard_option << 24) + (tfxU32)p.data.image_frame;
						s.transform.captured_position = p.data.captured_position;
						s.stretch = p.data.velocity_normal.w;
						s.intensity = p.data.intensity;
						s.handle = p.parent->current.image_handle;
						alignment_vector.y += 0.002f;	//We don't want a 0 alignment normal
						s.alignment = Pack10bit(alignment_vector, properties.billboard_option & 0x00000003);
						s.image_ptr = properties.image->ptr;
						pm.sprites3d[properties.layer][sprite_index++] = s;

						p.data.captured_position = s.transform.position;
						p.next_ptr = pm.SetNextParticle(next_buffer_index, p);
						next_index++;
					}
					else {
						s.transform.scale = 0;
						s.intensity = 0;
						s.image_frame_plus = 0;
						s.image_ptr = properties.image->ptr;
						pm.sprites3d[properties.layer][sprite_index++] = s;
						p.next_ptr = nullptr;
					}

				}
				else {
					p.data.flags &= ~tfxParticleFlags_fresh;
					p.next_ptr = pm.SetNextParticle(next_buffer_index, p);
					next_index++;
				}

				if (!(pm.flags & tfxEffectManagerFlags_guarantee_order) && next_index > 1 && pm.particle_banks[next_buffer_index][next_index - 2].data.depth < pm.particle_banks[next_buffer_index][next_index - 1].data.depth) {
					tfxParticle *p1 = &pm.particle_banks[next_buffer_index][next_index - 2];
					tfxParticle *p2 = &pm.particle_banks[next_buffer_index][next_index - 1];
					tfxParticle tmp = *p1;
					tfxParticleSprite3d *s1 = &pm.sprites3d[(p1->sprite_index & 0xF0000000) >> 28][p1->sprite_index & 0x0FFFFFFF];
					tfxParticleSprite3d *s2 = &pm.sprites3d[(p2->sprite_index & 0xF0000000) >> 28][p2->sprite_index & 0x0FFFFFFF];
					tfxParticleSprite3d sprite_tmp = *s1;
					*p1 = *p2;
					*p2 = tmp;
					tfxU32 tmp_index = p1->sprite_index;
					p1->sprite_index = p2->sprite_index;
					p2->sprite_index = tmp_index;
					*s1 = *s2;
					*s2 = sprite_tmp;
					pm.particle_banks[current_buffer_index][pm.particle_banks[next_buffer_index][next_index - 2].prev_index].next_ptr = p1;
					pm.particle_banks[current_buffer_index][pm.particle_banks[next_buffer_index][next_index - 1].prev_index].next_ptr = p2;

					//tfxParticle tmp = pm.particle_banks[next_buffer_index][next_index - 2];
					//pm.particle_banks[next_buffer_index][next_index - 2] = pm.particle_banks[next_buffer_index][next_index - 1];
					//pm.particle_banks[next_buffer_index][next_index - 1] = tmp;
					//particles[current_buffer_index][particles[next_buffer_index][next_index - 2].prev_index].next_ptr = next_index - 2;
					//particles[current_buffer_index][particles[next_buffer_index][next_index - 1].prev_index].next_ptr = next_index - 1;
				}

				index++;
				if (index == pm.new_particles_index_start[layer] && new_index != 0)
					break;
			}
			if (new_particle && new_index != 0) {
				while (new_index < pm.particle_banks[current_buffer_index].current_size) {
					new_particle = &pm.particle_banks[current_buffer_index][new_index];
					tfxEmitterProperties &new_properties = new_particle->parent->GetProperties();
					tfxParticleSprite3d &ns = pm.sprites3d[new_properties.layer][new_particle->sprite_index & 0x0FFFFFFF];
					pm.sprites3d[new_properties.layer][next_index] = ns;
					new_particle->prev_index = new_index;
					new_particle->data.flags &= ~tfxParticleFlags_fresh;
					new_particle->parent = new_particle->parent->next_ptr;
					new_particle->next_ptr = pm.SetNextParticle(next_buffer_index, *new_particle);
					new_index++;
				}
			}
			if (!(pm.flags & tfxEffectManagerFlags_guarantee_order) && pm.sort_passes > 0) {
				for (tfxU32 sorts = 0; sorts != pm.sort_passes; ++sorts) {
					for (tfxU32 i = 1; i < pm.particle_banks[next_buffer_index].size(); ++i) {
						tfxParticle *p1 = &pm.particle_banks[next_buffer_index][i - 1];
						tfxParticle *p2 = &pm.particle_banks[next_buffer_index][i];
						if (p1->data.depth < p2->data.depth) {
							tfxParticleSprite3d *s1 = &pm.sprites3d[(p1->sprite_index & 0xF0000000) >> 28][i - 1];
							tfxParticleSprite3d *s2 = &pm.sprites3d[(p2->sprite_index & 0xF0000000) >> 28][i];
							tfxParticleSprite3d sprite_tmp = *s1;
							tfxParticle tmp = *p1;
							*p1 = *p2;
							*s1 = *s2;

							*p2 = tmp;
							*s2 = sprite_tmp;

							tfxU32 tmp_index = p1->sprite_index;
							p1->sprite_index = p2->sprite_index;
							p2->sprite_index = tmp_index;

							pm.particle_banks[current_buffer_index][p1->prev_index].next_ptr = p1;
							pm.particle_banks[current_buffer_index][p2->prev_index].next_ptr = p2;
						}
					}
				}
			}
			if (pm.flags & tfxEffectManagerFlags_guarantee_order) {
				InsertionSortParticles(pm.particle_banks[next_buffer_index], pm.particle_banks[current_buffer_index]);
			}
		}

		pm.current_pbuff = next_buffer;

		pm.flags &= ~tfxEffectManagerFlags_update_base_values;
	}

	void tfxParticleManager::UpdateParticleOrderOnly() {
		/*if (!(flags & tfxEffectManagerFlags_order_by_depth))
			return;
		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			int current_buffer_index = layer * 2 + pm.current_pbuff;
			for (auto &p : particles[current_buffer_index]) {
				p.data.depth = LengthVec3NoSqR(p.data.world_position - camera_position);
			}
			InsertionSortParticles(particles[current_buffer_index], particles[layer * 2 + !pm.current_pbuff]);
		}*/
	}

	inline tfxEffectEmitter* tfxParticleManager::SetNextEffect(tfxEffectEmitter &e, unsigned int buffer) {
		unsigned int index = effects[buffer].current_size++;
		assert(index < effects[buffer].capacity);
		effects[buffer][index] = e;
		return &effects[buffer][index];
	}

	tfxvec<tfxEffectEmitter> *tfxParticleManager::GetEffectBuffer() {
		return &effects[current_ebuff];
	}

	void tfxParticleManager::InitFor2d(tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfxParticleManagerModes mode) {
		max_effects = effects_limit;

		flags = 0;

		if (mode == tfxParticleManagerMode_unordered)
			flags = tfxEffectManagerFlags_unorderd;
		else if (mode == tfxParticleManagerMode_ordered_by_depth)
			flags = tfxEffectManagerFlags_order_by_depth;
		else if (mode == tfxParticleManagerMode_ordered_by_depth_guaranteed)
			flags = tfxEffectManagerFlags_order_by_depth | tfxEffectManagerFlags_guarantee_order;

		for (tfxEachLayer) {
			max_cpu_particles_per_layer[layer] = layer_max_values[layer];

#ifdef tfxTRACK_MEMORY
			memcpy(sprites2d[layer].name, "ParticleManager::sprites2d\0", 27);
#endif

			sprites2d[layer].reserve(layer_max_values[layer]);
		}

		if (!(flags & tfxEffectManagerFlags_unorderd)) {
			for (tfxEachLayer) {
				tfxring<tfxParticle> particles(tfxCONSTRUCTOR_VEC_INIT("Ordered ring particle list"));
				particle_banks.push_back(particles);
				particle_banks.back().reserve(layer_max_values[layer]);
				particle_banks.push_back(particles);
				particle_banks.back().reserve(layer_max_values[layer]);
			}
		}

#ifdef tfxTRACK_MEMORY
		memcpy(effects[0].name, "ParticleManager::effects0\0", 26);
		memcpy(effects[1].name, "ParticleManager::effects1\0", 26);
#endif

		effects[0].create_pool(max_effects);
		effects[1].create_pool(max_effects);
		effects[0].clear();
		effects[1].clear();
	}

	void tfxParticleManager::InitFor3d(tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfxParticleManagerModes mode) {
		max_effects = effects_limit;

		flags = 0;

		if (mode == tfxParticleManagerMode_unordered)
			flags = tfxEffectManagerFlags_unorderd;
		else if (mode == tfxParticleManagerMode_ordered_by_depth)
			flags = tfxEffectManagerFlags_order_by_depth;
		else if (mode == tfxParticleManagerMode_ordered_by_depth_guaranteed)
			flags = tfxEffectManagerFlags_order_by_depth | tfxEffectManagerFlags_guarantee_order;

		flags |= tfxEffectManagerFlags_3d_effects;

		for (tfxEachLayer) {
			max_cpu_particles_per_layer[layer] = layer_max_values[layer];

#ifdef tfxTRACK_MEMORY
			memcpy(sprites3d[layer].name, "ParticleManager::sprites3d\0", 27);
#endif

			sprites3d[layer].reserve(layer_max_values[layer]);
		}

		if (!(flags & tfxEffectManagerFlags_unorderd)) {
			tfxring<tfxParticle> particles(tfxCONSTRUCTOR_VEC_INIT("Ordered ring particle list"));
			for (tfxEachLayer) {
				particle_banks.push_back(particles);
				particle_banks.push_back(particles);
			}
		}

#ifdef tfxTRACK_MEMORY
		memcpy(effects[0].name, "ParticleManager::effects0\0", 26);
		memcpy(effects[1].name, "ParticleManager::effects1\0", 26);
#endif

		effects[0].create_pool(max_effects);
		effects[1].create_pool(max_effects);
		effects[0].clear();
		effects[1].clear();
	}

	void tfxParticleManager::InitFor2d(unsigned int effects_limit, tfxParticleManagerModes mode) {
		tfxU32 layer_max_values[tfxLAYERS];
		memset(layer_max_values, 0, 16);
		InitFor2d(layer_max_values, effects_limit, mode);
		flags |= tfxEffectManagerFlags_dynamic_sprite_allocation;
	}

	void tfxParticleManager::InitFor3d(unsigned int effects_limit, tfxParticleManagerModes mode) {
		tfxU32 layer_max_values[tfxLAYERS];
		memset(layer_max_values, 0, 16);
		InitFor3d(layer_max_values, effects_limit, mode);
		flags |= tfxEffectManagerFlags_dynamic_sprite_allocation;
	}

	void tfxParticleManager::CreateParticleBanksForEachLayer() {
		FreeParticleBanks();
		for (tfxEachLayer) {
			tfxring<tfxParticle> particles(tfxCONSTRUCTOR_VEC_INIT("Ordered ring particle list"));
			particle_banks.push_back(particles);
			particle_banks.back().reserve(max_cpu_particles_per_layer[layer]);
			particle_banks.push_back(particles);
			particle_banks.back().reserve(max_cpu_particles_per_layer[layer]);
		}
	}

	void tfxParticleManager::Reconfigure(tfxParticleManagerModes mode, bool is_3d) {
		if (flags & tfxEffectManagerFlags_unorderd && mode != tfxParticleManagerMode_unordered) {
			FreeParticleBanks();
			CreateParticleBanksForEachLayer();
		}
		else if (!(flags & tfxEffectManagerFlags_unorderd) && mode == tfxParticleManagerMode_unordered) {
			FreeParticleBanks();
		}

		if (mode == tfxParticleManagerMode_unordered)
			flags = tfxEffectManagerFlags_unorderd;
		else if (mode == tfxParticleManagerMode_ordered_by_depth)
			flags = tfxEffectManagerFlags_order_by_depth;
		else if (mode == tfxParticleManagerMode_ordered_by_depth_guaranteed)
			flags = tfxEffectManagerFlags_order_by_depth | tfxEffectManagerFlags_guarantee_order;
		else if (mode == tfxParticleManagerMode_ordered_by_age)
			flags = tfxEffectManagerFlags_ordered_by_age;

		if (is_3d)
			flags |= tfxEffectManagerFlags_3d_effects;

	}

	void tfxParticleManager::InitForBoth(tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfxParticleManagerModes mode) {
		max_effects = effects_limit;

		flags = 0;

		if (mode == tfxParticleManagerMode_unordered)
			flags = tfxEffectManagerFlags_unorderd;
		else if (mode == tfxParticleManagerMode_ordered_by_depth)
			flags = tfxEffectManagerFlags_order_by_depth;
		else if (mode == tfxParticleManagerMode_ordered_by_depth_guaranteed)
			flags = tfxEffectManagerFlags_order_by_depth | tfxEffectManagerFlags_guarantee_order;
		else if (mode == tfxParticleManagerMode_ordered_by_age)
			flags = tfxEffectManagerFlags_ordered_by_age;

		for (tfxEachLayer) {
			max_cpu_particles_per_layer[layer] = layer_max_values[layer];

#ifdef tfxTRACK_MEMORY
			memcpy(sprites2d[layer].name, "ParticleManager::sprites2d\0", 27);
			memcpy(sprites3d[layer].name, "ParticleManager::sprites3d\0", 27);
#endif

			sprites2d[layer].reserve(layer_max_values[layer]);
			sprites3d[layer].reserve(layer_max_values[layer]);
		}

		if (!(flags & tfxEffectManagerFlags_unorderd)) {
			CreateParticleBanksForEachLayer();
		}

#ifdef tfxTRACK_MEMORY
		memcpy(effects[0].name, "ParticleManager::effects0\0", 26);
		memcpy(effects[1].name, "ParticleManager::effects1\0", 26);
#endif

		effects[0].create_pool(max_effects);
		effects[1].create_pool(max_effects);
		effects[0].clear();
		effects[1].clear();
		flags |= tfxEffectManagerFlags_dynamic_sprite_allocation;
	}

	void tfxParticleManager::ClearAll(bool free_memory) {
		for (tfxEachLayer) {
			sprites2d[layer].clear();
			sprites3d[layer].clear();
			if (!(flags & tfxEffectManagerFlags_unorderd)) {
				particle_banks[layer].clear();
				particle_banks[layer + 1].clear();
			}
		}
		if (flags & tfxEffectManagerFlags_unorderd) {
			if (free_memory) {
				FreeParticleBanks();
				for (auto &bank : free_particle_banks.data) {
					bank.free_all();
				}
				free_particle_banks.FreeAll();
			}
			else {
				for (auto &bank : particle_banks) {
					bank.clear();
				}
				//particle_banks.clear();
			}
			for (auto &e : effects[current_ebuff]) {
				if (e.type == tfxEmitterType) {
					FreeParticleBank(e);
				}
			}
		}
		effects[0].clear();
		effects[1].clear();
		particle_id = 0;
	}
	void tfxParticleManager::FreeParticleBanks() {
		for (auto &bank : particle_banks) {
			bank.free_all();
		}
		particle_banks.clear();
	}
	void tfxParticleManager::SoftExpireAll() {
		for (auto &e : effects[current_ebuff]) {
			e.flags |= tfxEmitterStateFlags_stop_spawning;
		}
	}
	void tfxParticleManager::SetLookUpMode(tfxLookupMode mode) {
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

	void tfxParticleManager::UpdateBaseValues() {
		flags |= tfxEffectManagerFlags_update_base_values;
	}

	tfxU32 GrabParticleBank(tfxParticleManager &pm, tfxKey emitter_hash, tfxU32 reserve_amount) {
		if (pm.free_particle_banks.ValidKey(emitter_hash)) {
			tfxvec<tfxU32> &free_banks = pm.free_particle_banks.At(emitter_hash);
			if (free_banks.current_size) {
				return free_banks.pop_back();
			}
		}
		tfxring<tfxParticle> particles(tfxCONSTRUCTOR_VEC_INIT("Unordered ring particle list"));
		particles.reserve(reserve_amount);
		pm.particle_banks.push_back(particles);
		return pm.particle_banks.current_size - 1;
	}

	void StopSpawning(tfxParticleManager &pm) {
		pm.SoftExpireAll();
	}

	void RemoveAllEffects(tfxParticleManager &pm) {
		pm.ClearAll();
	}

	void AddEffect(tfxParticleManager &pm, tfxEffectEmitter &effect, float x, float y) {
		effect.common.transform.local_position.x = x;
		effect.common.transform.local_position.y = y;
		pm.AddEffect(effect, pm.current_ebuff);
	}

	void AddEffect(tfxParticleManager &pm, tfxEffectTemplate &effect_template, float x, float y) {
		effect_template.effect.common.transform.local_position.x = x;
		effect_template.effect.common.transform.local_position.y = y;
		pm.AddEffect(effect_template);
	}

	void Rotate(tfxEffectEmitter &e, float r) {
		e.common.transform.local_rotations.roll += r;
	}

	void SetAngle(tfxEffectEmitter &e, float a) {
		e.common.transform.local_rotations.roll = a;
	}

	void Scale(tfxEffectEmitter &e, const tfxVec3& s) {
		e.common.transform.scale = s;
	}

	void Scale(tfxEffectEmitter &e, float x, float y, float z) {
		e.common.transform.scale = tfxVec3(x, y, z);
	}

	void Position(tfxEffectEmitter &e, const tfxVec2 &p) {
		e.common.transform.local_position = p;
	}

	void Position(tfxEffectEmitter &e, const tfxVec3 &p) {
		e.common.transform.local_position = p;
	}

	void UpdatePMEmitter(tfxParticleManager &pm, tfxEffectEmitter &e) {
		tfxPROFILE;
		e.common.transform.captured_position = e.common.transform.world_position;

		if (pm.lookup_mode == tfxPrecise) {
			e.common.frame = e.common.age;
		}
		else {
			e.common.frame = e.common.age / tfxLOOKUP_FREQUENCY;
		}

		if (e.type == tfxEffectEmitterType::tfxEffectType) {
			pm.parent_spawn_controls = UpdateEffectState(e);
		}

		tfxEmitterProperties &properties = e.GetProperties();

		if (e.parent && e.parent->type != tfxFolder && e.parent->next_ptr) {
			e.parent = e.parent->next_ptr;

			e.parent->common.timeout_counter = 0;
			if (e.parent->common.age < properties.delay_spawning) {
				return;
			}

			//e.flags |= e.parent->flags & tfxEmitterStateFlags_stop_spawning;
			e.flags |= e.parent->flags & tfxEmitterStateFlags_no_tween;
			tfxEmitterSpawnControls spawn_controls = UpdateEmitterState(e, pm.parent_spawn_controls);

			bool is_compute = e.common.property_flags & tfxEmitterPropertyFlags_is_bottom_emitter && pm.flags & tfxEffectManagerFlags_use_compute_shader;
			tfxU32 amount_spawned = 0;
			tfxU32 max_spawn_count = NewSpritesNeeded(pm, e);

			if (pm.flags & tfxEffectManagerFlags_unorderd) {

				if (e.common.property_flags & tfxEmitterPropertyFlags_is_3d) {
					tfxring<tfxParticleSprite3d> &sprite_buffer = pm.sprites3d[properties.layer];

					Transform3d(e.common.transform, e.parent->common.transform);

					if (e.flags & tfxEmitterStateFlags_no_tween_this_update || e.flags & tfxEmitterStateFlags_no_tween) {
						e.common.transform.captured_position = e.common.transform.world_position;
					}

					e.sprites_count = pm.particle_banks[e.particles_index].current_size;
					if (pm.flags & tfxEffectManagerFlags_dynamic_sprite_allocation && e.sprites_count + max_spawn_count > sprite_buffer.free_space()) {
						sprite_buffer.reserve(sprite_buffer._grow_capacity(sprite_buffer.capacity + (e.sprites_count + max_spawn_count - sprite_buffer.free_space()) + 1));
						sprite_buffer.current_size += e.sprites_count;
					}
					else if (!(pm.flags & tfxEffectManagerFlags_dynamic_sprite_allocation)) {
						e.sprites_count = e.sprites_count > sprite_buffer.free_space() ? sprite_buffer.free_space() : e.sprites_count;
						sprite_buffer.current_size += e.sprites_count;
						max_spawn_count = max_spawn_count > sprite_buffer.free_space() ? sprite_buffer.free_space() : max_spawn_count;
					}
					else {
						sprite_buffer.current_size += e.sprites_count;
					}
					e.sprites_count += max_spawn_count;
					sprite_buffer.current_size += max_spawn_count;
					e.sprites_index = pm.sprite_index_point[properties.layer];
					pm.sprite_index_point[properties.layer] += e.sprites_count;

					if (e.flags & tfxEmitterStateFlags_is_sub_emitter) {
						if (e.common.age > 0 && !(pm.flags & tfxEffectManagerFlags_disable_spawning))
							amount_spawned = SpawnParticles3d(pm, e, spawn_controls, max_spawn_count);
					}
					else {
						if (!(pm.flags & tfxEffectManagerFlags_disable_spawning))
							amount_spawned = SpawnParticles3d(pm, e, spawn_controls, max_spawn_count);
					}
					sprite_buffer.current_size -= (max_spawn_count - amount_spawned);
					pm.sprite_index_point[properties.layer] -= (max_spawn_count - amount_spawned);

					ControlParticles3d(pm, e, amount_spawned);
				}
				else {
					tfxring<tfxParticleSprite2d> &sprite_buffer = pm.sprites2d[properties.layer];

					Transform(e.common.transform, e.parent->common.transform);

					if (e.flags & tfxEmitterStateFlags_no_tween_this_update || e.flags & tfxEmitterStateFlags_no_tween) {
						e.common.transform.captured_position = e.common.transform.world_position;
					}

					e.sprites_count = pm.particle_banks[e.particles_index].current_size;
					if (pm.flags & tfxEffectManagerFlags_dynamic_sprite_allocation && e.sprites_count + max_spawn_count > sprite_buffer.free_space()) {
						sprite_buffer.reserve(sprite_buffer._grow_capacity(sprite_buffer.capacity + (e.sprites_count + max_spawn_count - sprite_buffer.free_space()) + 1));
						sprite_buffer.current_size += e.sprites_count;
					}
					else if (!(pm.flags & tfxEffectManagerFlags_dynamic_sprite_allocation)) {
						e.sprites_count = e.sprites_count > sprite_buffer.free_space() ? sprite_buffer.free_space() : e.sprites_count;
						sprite_buffer.current_size += e.sprites_count;
						max_spawn_count = max_spawn_count > sprite_buffer.free_space() ? sprite_buffer.free_space() : max_spawn_count;
					}
					else {
						sprite_buffer.current_size += e.sprites_count;
					}
					e.sprites_count += max_spawn_count;
					sprite_buffer.current_size += max_spawn_count;
					e.sprites_index = pm.sprite_index_point[properties.layer];
					pm.sprite_index_point[properties.layer] += e.sprites_count;

					if (e.flags & tfxEmitterStateFlags_is_sub_emitter) {
						if (e.common.age > 0 && !(pm.flags & tfxEffectManagerFlags_disable_spawning))
							amount_spawned = SpawnParticles2d(pm, e, spawn_controls, max_spawn_count);
					}
					else {
						if (!(pm.flags & tfxEffectManagerFlags_disable_spawning))
							amount_spawned = SpawnParticles2d(pm, e, spawn_controls, max_spawn_count);
					}
					sprite_buffer.current_size -= (max_spawn_count - amount_spawned);
					pm.sprite_index_point[properties.layer] -= (max_spawn_count - amount_spawned);

					ControlParticles2d(pm, e, amount_spawned);
				}
			}
			else {
				e.sprites_index = pm.sprite_index_point[properties.layer];
				pm.sprite_index_point[properties.layer] += max_spawn_count;
				e.particles_index = properties.layer * 2 + pm.current_pbuff;
				if (e.common.property_flags & tfxEmitterPropertyFlags_is_3d) {
					Transform3d(e.common.transform, e.parent->common.transform);

					if (e.flags & tfxEmitterStateFlags_no_tween_this_update || e.flags & tfxEmitterStateFlags_no_tween) {
						e.common.transform.captured_position = e.common.transform.world_position;
					}

					tfxring<tfxParticleSprite3d> &sprite_buffer = pm.sprites3d[properties.layer];
					if (pm.flags & tfxEffectManagerFlags_dynamic_sprite_allocation && pm.sprite_index_point[properties.layer] > sprite_buffer.free_space()) {
						sprite_buffer.reserve(sprite_buffer._grow_capacity(sprite_buffer.capacity + (pm.sprite_index_point[properties.layer] - sprite_buffer.free_space()) + 1));
					}
					sprite_buffer.current_size = pm.sprite_index_point[properties.layer];

					if (e.flags & tfxEmitterStateFlags_is_sub_emitter) {
						if (e.common.age > 0 && e.common.property_flags & tfxEmitterPropertyFlags_is_3d && !(pm.flags & tfxEffectManagerFlags_disable_spawning))
							amount_spawned = SpawnParticles3d(pm, e, spawn_controls, max_spawn_count);
					}
					else {
						if (e.common.property_flags & tfxEmitterPropertyFlags_is_3d && !(pm.flags & tfxEffectManagerFlags_disable_spawning))
							amount_spawned = SpawnParticles3d(pm, e, spawn_controls, max_spawn_count);
					}
					sprite_buffer.current_size -= (max_spawn_count - amount_spawned);
					pm.sprite_index_point[properties.layer] -= (max_spawn_count - amount_spawned);
				}
				else {
					Transform(e.common.transform, e.parent->common.transform);

					if (e.flags & tfxEmitterStateFlags_no_tween_this_update || e.flags & tfxEmitterStateFlags_no_tween) {
						e.common.transform.captured_position = e.common.transform.world_position;
					}

					tfxring<tfxParticleSprite2d> &sprite_buffer = pm.sprites2d[properties.layer];
					if (pm.flags & tfxEffectManagerFlags_dynamic_sprite_allocation && pm.sprite_index_point[properties.layer] > sprite_buffer.free_space()) {
						sprite_buffer.reserve(sprite_buffer._grow_capacity(sprite_buffer.capacity + (pm.sprite_index_point[properties.layer] - sprite_buffer.free_space()) + 1));
					}
					sprite_buffer.current_size = pm.sprite_index_point[properties.layer];

					if (e.flags & tfxEmitterStateFlags_is_sub_emitter) {
						if (e.common.age > 0 && !(pm.flags & tfxEffectManagerFlags_disable_spawning))
							amount_spawned = SpawnParticles2d(pm, e, spawn_controls, max_spawn_count);
					}
					else {
						if (!(pm.flags & tfxEffectManagerFlags_disable_spawning))
							amount_spawned = SpawnParticles2d(pm, e, spawn_controls, max_spawn_count);
					}
					sprite_buffer.current_size -= (max_spawn_count - amount_spawned);
					pm.sprite_index_point[properties.layer] -= (max_spawn_count - amount_spawned);
				}
			}

			e.parent->highest_particle_age = e.highest_particle_age;
		}
		else if (e.parent && e.parent->type != tfxFolder && !e.parent->next_ptr) {
			e.common.timeout_counter = e.common.timeout;
			return;
		}
		else if (e.parent_particle) {
			//e.flags |= e.parent_particle->data.flags & tfxParticleFlags_remove;
			if (e.parent_particle->next_ptr) {
				e.parent_particle = e.parent_particle->next_ptr;
				tfxU32 sprite_layer = (e.parent_particle->sprite_index & 0xF0000000) >> 28;
				tfxU32 sprite_index = e.parent_particle->sprite_index & 0x0FFFFFFF;
				if (e.common.property_flags & tfxEmitterPropertyFlags_is_3d)
					TransformEffector3d(e, pm.sprites3d[sprite_layer][sprite_index].transform, true, e.common.property_flags & tfxEmitterPropertyFlags_relative_angle);
				else
					TransformEffector(e, pm.sprites2d[properties.layer][sprite_index].transform, true, e.common.property_flags & tfxEmitterPropertyFlags_relative_angle);

				e.common.transform.world_position += properties.emitter_handle * e.current.overal_scale;
				if (e.flags & tfxEmitterStateFlags_no_tween_this_update || e.flags & tfxEmitterStateFlags_no_tween) {
					e.common.transform.captured_position = e.common.transform.world_position;
				}
			}
			else {
				e.parent_particle = nullptr;
				e.flags |= tfxEmitterStateFlags_retain_matrix;
				e.common.transform.local_position = e.common.transform.world_position;
				e.common.transform.local_rotations.roll = e.common.transform.world_rotations.roll;
				e.flags |= tfxEmitterStateFlags_stop_spawning;
			}
		}
		else {
			if (!(e.flags & tfxEmitterStateFlags_retain_matrix)) {
				e.common.transform.world_position = e.common.transform.local_position;
				e.common.transform.world_position += properties.emitter_handle * e.current.overal_scale;
				if (e.common.property_flags & tfxEmitterPropertyFlags_is_3d) {
					e.common.transform.world_rotations = e.common.transform.local_rotations;
					tfxMatrix4 roll = mmZRotate(e.common.transform.local_rotations.roll);
					tfxMatrix4 pitch = mmXRotate(e.common.transform.local_rotations.pitch);
					tfxMatrix4 yaw = mmYRotate(e.common.transform.local_rotations.yaw);
					e.common.transform.matrix = mmTransform(yaw, pitch);
					e.common.transform.matrix = mmTransform(e.common.transform.matrix, roll);
				}
				else {
					e.common.transform.world_rotations.roll = e.common.transform.local_rotations.roll;
					float s = sin(e.common.transform.local_rotations.roll);
					float c = cos(e.common.transform.local_rotations.roll);
					e.common.transform.matrix.Set2(c, s, -s, c);
				}
			}

			if (e.flags & tfxEmitterStateFlags_no_tween_this_update || e.flags & tfxEmitterStateFlags_no_tween) {
				e.common.transform.captured_position = e.common.transform.world_position;
			}

		}

		e.common.age += tfxFRAME_LENGTH;
		if (!(e.common.property_flags & tfxEmitterPropertyFlags_single) || properties.single_shot_limit > 0)
			e.highest_particle_age -= tfxFRAME_LENGTH;

		if (properties.loop_length && e.common.age > properties.loop_length)
			e.common.age = 0;

		if (e.highest_particle_age <= 0 && e.common.age > tfxFRAME_LENGTH * 5.f) {
			e.common.timeout_counter++;
		}
		else {
			e.common.timeout_counter = 0;
		}

		e.flags &= ~tfxEmitterStateFlags_no_tween_this_update;
	}

	tfxU32 NewSpritesNeeded(tfxParticleManager &pm, tfxEffectEmitter &e) {

		if (e.flags & tfxEmitterStateFlags_single_shot_done || e.parent->flags & tfxEmitterStateFlags_stop_spawning)
			return 0;
		if (e.current.qty == 0)
			return 0;

		tfxEmitterProperties &properties = e.GetProperties();

		float qty_step_size = 0;
		if (!(e.common.property_flags & tfxEmitterPropertyFlags_single)) {
			if (e.common.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && (properties.emission_type == tfxArea || properties.emission_type == tfxEllipse)) {
				if (e.Is3DEffect()) {
					float area = std::fmaxf(0.1f, e.current.emitter_size.x) * std::fmaxf(0.1f, e.current.emitter_size.y) * std::fmaxf(0.1f, e.current.emitter_size.z);
					e.current.qty = (e.current.qty / 50.f) * area;
				}
				else {
					float area = e.current.emitter_size.x * e.current.emitter_size.y;
					e.current.qty = (e.current.qty / 10000.f) * area;
				}
			}
			else if (e.common.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && properties.emission_type == tfxLine) {
				e.current.qty = (e.current.qty / 100.f) * e.current.emitter_size.y;
			}

			e.current.qty *= tfxUPDATE_TIME;
			qty_step_size = 1.f / e.current.qty;
		}
		else {
			qty_step_size = 1.f / e.current.qty;
		}

		float tween = e.current.amount_remainder;
		return tween >= 1.f ? 0 : tfxU32((1.f - e.current.amount_remainder) / qty_step_size) + 1;
	}

	tfxU32 SpawnParticles2d(tfxParticleManager &pm, tfxEffectEmitter &e, tfxEmitterSpawnControls &spawn_controls, tfxU32 max_spawn_amount) {
		tfxEmitterProperties &properties = e.GetProperties();

		if (e.flags & tfxEmitterStateFlags_single_shot_done || e.parent->flags & tfxEmitterStateFlags_stop_spawning)
			return 0;
		if (e.current.qty == 0)
			return 0;

		float qty_step_size = 1.f / e.current.qty;
		float tween = e.current.amount_remainder;
		bool is_compute = e.common.property_flags & tfxEmitterPropertyFlags_is_bottom_emitter && pm.flags & tfxEffectManagerFlags_use_compute_shader;
		tfxU32 amount_spawned = 0;

		while (tween < 1.f) {
			if (amount_spawned >= max_spawn_amount) {
				e.current.amount_remainder = 0;
				break;
			}

			bool is_single = e.common.property_flags & tfxEmitterPropertyFlags_single;

			tfxParticle *p = &pm.GrabCPUParticle(e.particles_index);
			assert(e.sprites_index < pm.sprites2d[properties.layer].capacity);
			p->sprite_index = (properties.layer << 28) + e.sprites_index;
			tfxParticleSprite2d &s = pm.sprites2d[properties.layer][e.sprites_index++];
			InitCPUParticle2d(pm, e, *p, s.transform, spawn_controls, tween);

			e.highest_particle_age = std::fmaxf(e.highest_particle_age, p->data.max_age);
			e.parent->highest_particle_age = e.highest_particle_age + tfxFRAME_LENGTH;

			tween += qty_step_size;

			s.color = p->data.color;
			s.image_frame = (tfxU32)p->data.image_frame;
			s.intensity = p->data.intensity;
			s.handle = e.current.image_handle;
			s.image_ptr = properties.image->ptr;

			amount_spawned++;
		}

		e.current.amount_remainder = tween - 1.f;
		return amount_spawned;
	}

	tfxU32 SpawnParticles3d(tfxParticleManager &pm, tfxEffectEmitter &e, tfxEmitterSpawnControls &spawn_controls, tfxU32 max_spawn_amount) {
		tfxPROFILE;

		if (e.flags & tfxEmitterStateFlags_single_shot_done || e.parent->flags & tfxEmitterStateFlags_stop_spawning)
			return 0;
		if (e.current.qty == 0)
			return 0;

		tfxEmitterProperties &properties = e.GetProperties();

		float qty_step_size = 1.f / e.current.qty;
		float tween = e.current.amount_remainder;
		bool is_compute = e.common.property_flags & tfxEmitterPropertyFlags_is_bottom_emitter && pm.flags & tfxEffectManagerFlags_use_compute_shader;
		float positions_qty = e.current.qty;
		tfxU32 amount_spawned = 0;

		while (tween < 1.f) {
			if (amount_spawned >= max_spawn_amount) {
				e.current.amount_remainder = 0;
				break;
			}
			tfxSpawnPosition new_position = InitialisePosition3d(e.current, e.common, spawn_controls, &e, tween);
			if (pm.flags & tfxEffectManagerFlags_order_by_depth) {
				new_position.distance_to_camera = LengthVec3NoSqR(new_position.world_position - pm.camera_position);
			}
			pm.new_positions.push_back(new_position);
			tween += qty_step_size;

			amount_spawned++;
		}

		if (pm.flags & tfxEffectManagerFlags_order_by_depth) {
			std::qsort(pm.new_positions.data, pm.new_positions.current_size, sizeof(tfxSpawnPosition), SortParticles);
		}

		e.current.amount_remainder = tween - 1.f;

		pm.new_particles_index_start[properties.layer] = std::min(pm.new_particles_index_start[properties.layer], pm.particle_banks[e.particles_index].current_size);

		for (auto &position : pm.new_positions) {
			bool is_single = e.common.property_flags & tfxEmitterPropertyFlags_single;

			tfxParticle *p = &pm.GrabCPUParticle(e.particles_index);

			assert(e.sprites_index < pm.sprites3d[properties.layer].capacity);
			p->sprite_index = (properties.layer << 28) + e.sprites_index;
			tfxParticleSprite3d &s = pm.sprites3d[properties.layer][e.sprites_index++];

			p->data.local_position = position.local_position;
			p->data.captured_position = s.transform.captured_position = position.captured_position;
			s.transform.position = position.world_position;
			p->data.weight_acceleration = position.weight_acceleration;
			p->data.base.weight = position.base_weight;
			p->data.base.velocity = position.base_velocity;
			p->data.velocity_normal = position.velocity_normal;
			p->data.depth = position.distance_to_camera;
			InitCPUParticle3d(pm, e, *p, s.transform, spawn_controls, tween);

			e.highest_particle_age = std::fmaxf(e.highest_particle_age, p->data.max_age);
			e.parent->highest_particle_age = e.highest_particle_age + tfxFRAME_LENGTH;

			bool line = e.common.property_flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type == tfxEmissionType::tfxLine;
			tfxVec3 alignment_vector(0.f, 0.002f, 0.f);
			if (p->data.flags & tfxParticleFlags_fresh && properties.vector_align_type == tfxVectorAlignType_motion) {
				alignment_vector = p->data.velocity_normal.xyz();
				alignment_vector = FastNormalizeVec(alignment_vector);
			}
			else if (properties.vector_align_type == tfxVectorAlignType_motion) {
				alignment_vector = s.transform.position - p->data.captured_position;
				float l = FastLength(alignment_vector);
				p->data.velocity_normal.w *= l * 10.f;
				alignment_vector = FastNormalizeVec(alignment_vector);
			}
			else if (properties.vector_align_type == tfxVectorAlignType_emission) {
				alignment_vector = p->data.velocity_normal.xyz();
			}
			else if (properties.vector_align_type == tfxVectorAlignType_emitter) {
				alignment_vector = mmTransformVector(e.common.transform.matrix, tfxVec4(0.f, 1.f, 0.f, 0.f)).xyz();
			}
			else {
				if (!(e.common.property_flags & tfxEmitterPropertyFlags_fill_area)) {
					if (properties.vector_align_type == tfxVectorAlignType_surface_normal) {
						tfxVec3 emitter_size = e.current.emitter_size * .5f;
						alignment_vector = tfxVec3(p->data.local_position.x / (emitter_size.x * emitter_size.x),
							p->data.local_position.y / (emitter_size.y * emitter_size.y),
							p->data.local_position.z / (emitter_size.z * emitter_size.z));
						alignment_vector *= 2.f;
					}
				}
			}
			//p->data.flags &= ~tfxParticleFlags_fresh;

			s.color = p->data.color;
			s.image_frame_plus = (properties.billboard_option << 24) + (tfxU32)p->data.image_frame;
			s.stretch = p->data.velocity_normal.w;
			s.intensity = p->data.intensity;
			s.handle = e.current.image_handle;
			s.alignment = Pack10bit(alignment_vector, properties.billboard_option & 0x00000003);
			s.image_ptr = properties.image->ptr;

		}

		pm.new_positions.clear();

		return amount_spawned;
	}

	void InitCPUParticle2d(tfxParticleManager &pm, tfxEffectEmitter &e, tfxParticle &p, tfxSpriteTransform2d &sprite_transform, tfxEmitterSpawnControls &spawn_controls, float tween) {
		tfxPROFILE;
		p.data.flags = tfxParticleFlags_fresh;
		p.next_ptr = &p;
		p.parent = &e;

		if (e.common.property_flags & tfxEmitterPropertyFlags_single)
			e.flags |= tfxEmitterStateFlags_single_shot_done;

		//----Properties

		//Set base values-------------------------------

		//----Life
		p.data.max_age = spawn_controls.life + random_generation.Range(spawn_controls.life_variation);
		p.data.age = 0.f;
		p.data.single_loop_count = 1;

		InitialiseParticle2d(p.data, sprite_transform, e.current, e.common, spawn_controls, &e, tween);

		if (e.GetInfo().sub_effectors.size()) {
			for (auto &sub : e.GetInfo().sub_effectors) {
				if (!pm.FreeEffectCapacity())
					break;
				sub.parent = nullptr;
				sub.parent_particle = &p;
				sub.highest_particle_age = p.data.max_age;
				sub.current.overal_scale = e.current.overal_scale;
				sub.flags |= e.flags & tfxEmitterStateFlags_no_tween;
				sub.SetTimeout(5);
				pm.AddEffect(sub, pm.current_ebuff, true);
			}
		}
	}

	void InitCPUParticle3d(tfxParticleManager &pm, tfxEffectEmitter &e, tfxParticle &p, tfxSpriteTransform3d &sprite_transform, tfxEmitterSpawnControls &spawn_controls, float tween) {
		tfxPROFILE;
		p.data.flags = tfxParticleFlags_fresh;
		p.next_ptr = &p;
		p.parent = &e;

		if (e.common.property_flags & tfxEmitterPropertyFlags_single)
			e.flags |= tfxEmitterStateFlags_single_shot_done;

		//----Properties

		//Set base values-------------------------------

		//----Life
		p.data.max_age = spawn_controls.life + random_generation.Range(spawn_controls.life_variation);
		p.data.age = 0.f;
		p.data.single_loop_count = 1;

		InitialiseParticle3d(p.data, sprite_transform, e.current, e.common, spawn_controls, &e, tween);

		if (e.GetInfo().sub_effectors.size()) {
			for (auto &sub : e.GetInfo().sub_effectors) {
				if (!pm.FreeEffectCapacity())
					break;
				sub.parent = nullptr;
				sub.parent_particle = &p;
				sub.highest_particle_age = p.data.max_age;
				sub.current.overal_scale = e.current.overal_scale;
				sub.flags |= e.flags & tfxEmitterStateFlags_no_tween;
				pm.AddEffect(sub, pm.current_ebuff, true);
			}
		}
	}

	tfxEmitterSpawnControls UpdateEmitterState(tfxEffectEmitter &e, tfxParentSpawnControls &parent_spawn_controls) {
		tfxPROFILE;
		tfxEffectEmitter &parent = *e.parent;
		tfxEmitterProperties &properties = e.GetProperties();
		tfxEmitterSpawnControls spawn_controls;

		e.current.overal_scale = parent.current.overal_scale;
		spawn_controls.life = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].base.life, e.common.frame) * parent_spawn_controls.life;
		spawn_controls.life_variation = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].variation.life, e.common.frame) * parent_spawn_controls.life;
		if (!(e.common.property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			spawn_controls.size.x = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].base.width, e.common.frame) * parent_spawn_controls.size_x;
			spawn_controls.size.y = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].base.height, e.common.frame) * parent_spawn_controls.size_y;
		}
		else {
			spawn_controls.size.x = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].base.width, e.common.frame);
			if (parent.common.property_flags & tfxEmitterPropertyFlags_global_uniform_size)
				spawn_controls.size.y = spawn_controls.size.x * parent_spawn_controls.size_x;
			else
				spawn_controls.size.y = spawn_controls.size.x * parent_spawn_controls.size_y;
			spawn_controls.size.x *= parent_spawn_controls.size_x;
		}
		spawn_controls.size_variation.x = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].variation.width, e.common.frame) * parent_spawn_controls.size_x;
		spawn_controls.size_variation.y = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].variation.height, e.common.frame) * parent_spawn_controls.size_y;
		spawn_controls.velocity = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].base.velocity, e.common.frame) * parent_spawn_controls.velocity;
		spawn_controls.velocity_variation = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].variation.velocity, e.common.frame) * parent_spawn_controls.velocity;
		e.current.velocity_adjuster = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].overtime.velocity_adjuster, e.common.frame);
		spawn_controls.spin = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].base.spin, e.common.frame) * parent_spawn_controls.spin;
		spawn_controls.spin_variation = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].variation.spin, e.common.frame) * parent_spawn_controls.spin;
		e.common.transform.local_rotations.roll = LookupPrecise(e.common.library->emitter_attributes[e.emitter_attributes].properties.roll, e.common.age);
		e.common.transform.local_rotations.pitch = LookupPrecise(e.common.library->emitter_attributes[e.emitter_attributes].properties.pitch, e.common.age);
		e.common.transform.local_rotations.yaw = LookupPrecise(e.common.library->emitter_attributes[e.emitter_attributes].properties.yaw, e.common.age);
		e.current.intensity = parent_spawn_controls.intensity;
		spawn_controls.splatter = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].properties.splatter, e.common.frame) * parent_spawn_controls.splatter;
		e.current.emitter_size.y = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].properties.emitter_height, e.common.frame);
		spawn_controls.weight = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].base.weight, e.common.frame) * parent_spawn_controls.weight;
		spawn_controls.weight_variation = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].variation.weight, e.common.frame) * parent_spawn_controls.weight;
		spawn_controls.noise_offset_variation = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].variation.noise_offset, e.common.frame);
		spawn_controls.noise_offset = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].base.noise_offset, e.common.frame);
		spawn_controls.noise_resolution = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].variation.noise_resolution, e.common.frame);
		e.current.stretch = parent.current.stretch;
		e.common.transform.scale = parent.common.transform.scale;

		if (!(e.common.property_flags & tfxEmitterPropertyFlags_single)) {
			e.current.qty = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].base.amount, e.common.frame);
			float amount_variation = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].variation.amount, e.common.frame);
			e.current.qty += amount_variation > 0.f ? random_generation.Range(1.f, amount_variation) : 0.f;

			e.current.qty *= lookup_callback(e.parent->common.library->global_graphs[e.parent->global].amount, e.common.frame);
		}
		else {
			e.current.qty = (float)properties.spawn_amount;
		}

		//----Handle
		if (e.common.property_flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
			e.current.image_handle = tfxVec2(0.5f, 0.5f);
		}
		else {
			e.current.image_handle = properties.image_handle;
		}

		bool is_area = properties.emission_type == tfxEmissionType::tfxArea || properties.emission_type == tfxEmissionType::tfxEllipse;

		if (is_area) {
			e.current.emitter_size.x = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].properties.emitter_width, e.common.frame);
		}
		else
			e.current.emitter_size.x = 0.f;

		if (e.common.property_flags & tfxEmitterPropertyFlags_is_3d && is_area) {
			e.current.emitter_size.z = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].properties.emitter_depth, e.common.frame);
		}

		e.current.emitter_size *= parent.current.emitter_size;

		if (properties.emission_type == tfxEmissionType::tfxEllipse) {
			spawn_controls.arc_size = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].properties.arc_size, e.common.frame);
			spawn_controls.arc_offset = lookup_callback(e.common.library->emitter_attributes[e.emitter_attributes].properties.arc_offset, e.common.frame);
		}

		if (e.common.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center && properties.emission_type != tfxEmissionType::tfxPoint) {
			if (properties.emission_type == tfxEmissionType::tfxEllipse && e.common.property_flags & tfxEmitterPropertyFlags_is_3d)
				e.common.handle = e.current.emitter_size * 0.f;
			else if (e.common.property_flags & tfxEmitterPropertyFlags_is_3d)
				e.common.handle = e.current.emitter_size * -0.5f;
			else if (properties.emission_type == tfxLine)
				e.common.handle = e.current.emitter_size * 0.5f;
			else
				e.common.handle = e.current.emitter_size * -0.5f;
		}
		else if (!(e.common.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center)) {
			e.common.handle = properties.emitter_handle;
		}
		else {
			e.common.handle = 0.f;
		}

		if (e.common.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {
			if (properties.emission_type == tfxEmissionType::tfxArea) {
				if (properties.grid_points.x > 1)
					spawn_controls.grid_segment_size.x = e.current.emitter_size.x / (properties.grid_points.x - 1);
				if (properties.grid_points.y > 1)
					spawn_controls.grid_segment_size.y = e.current.emitter_size.y / (properties.grid_points.y - 1);
				if (properties.grid_points.z > 1)
					spawn_controls.grid_segment_size.z = e.current.emitter_size.z / (properties.grid_points.z - 1);
			}
			else if (properties.emission_type == tfxEmissionType::tfxEllipse) {
				if (properties.grid_points.x > 0)
					spawn_controls.grid_segment_size.x = spawn_controls.arc_size / (properties.grid_points.x);
			}
			else if (properties.emission_type == tfxEmissionType::tfxLine) {
				if (properties.grid_points.x > 1)
					spawn_controls.grid_segment_size.y = e.current.emitter_size.y / (properties.grid_points.x - 1);
			}
		}

		if (e.update_emitter_callback)
			e.update_emitter_callback(e, spawn_controls);

		return spawn_controls;
	}

	tfxParentSpawnControls UpdateEffectState(tfxEffectEmitter &e) {
		tfxPROFILE;
		//If this effect is a sub effect then the graph index will reference the global graphs for the root parent effect
		tfxParentSpawnControls spawn_controls;
		spawn_controls.life = lookup_callback(e.common.library->global_graphs[e.global].life, e.common.frame);
		if (!(e.common.property_flags & tfxEmitterPropertyFlags_global_uniform_size)) {
			spawn_controls.size_x = lookup_callback(e.common.library->global_graphs[e.global].width, e.common.frame);
			spawn_controls.size_y = lookup_callback(e.common.library->global_graphs[e.global].height, e.common.frame);
		}
		else {
			spawn_controls.size_x = lookup_callback(e.common.library->global_graphs[e.global].width, e.common.frame);
			spawn_controls.size_y = spawn_controls.size_x;
		}
		spawn_controls.velocity = lookup_callback(e.common.library->global_graphs[e.global].velocity, e.common.frame);
		spawn_controls.spin = lookup_callback(e.common.library->global_graphs[e.global].spin, e.common.frame);
		spawn_controls.intensity = lookup_callback(e.common.library->global_graphs[e.global].intensity, e.common.frame);
		spawn_controls.splatter = lookup_callback(e.common.library->global_graphs[e.global].splatter, e.common.frame);
		e.current.emitter_size.x = lookup_callback(e.common.library->global_graphs[e.global].emitter_width, e.common.frame);
		e.current.emitter_size.y = lookup_callback(e.common.library->global_graphs[e.global].emitter_height, e.common.frame);
		e.current.emitter_size.z = lookup_callback(e.common.library->global_graphs[e.global].emitter_depth, e.common.frame);
		//We don't want to scale twice when the sub effect is transformed, so the values here are set to 1. That means that the root effect will only control the global scale.
		e.current.overal_scale = lookup_callback(e.common.library->global_graphs[e.global].overal_scale, e.common.frame);
		if (!e.parent_particle) {
			e.common.transform.scale.x = e.current.overal_scale;
			e.common.transform.scale.y = e.current.overal_scale;
			e.common.transform.scale.z = e.current.overal_scale;
			e.common.transform.local_rotations.roll = LookupPrecise(e.common.library->global_graphs[e.global].roll, e.common.age);
			e.common.transform.local_rotations.pitch = LookupPrecise(e.common.library->global_graphs[e.global].pitch, e.common.age);
			e.common.transform.local_rotations.yaw = LookupPrecise(e.common.library->global_graphs[e.global].yaw, e.common.age);
		}
		else {
			e.common.transform.scale.x = e.current.overal_scale;
			e.common.transform.scale.y = e.current.overal_scale;
			e.common.transform.scale.z = e.current.overal_scale;
			e.common.transform.local_rotations.roll = 0.f;
			e.common.transform.local_rotations.pitch = 0.f;
			e.common.transform.local_rotations.yaw = 0.f;
		}
		e.current.stretch = lookup_callback(e.common.library->global_graphs[e.global].stretch, e.common.frame);
		spawn_controls.weight = lookup_callback(e.common.library->global_graphs[e.global].weight, e.common.frame);

		if (e.update_effect_callback)
			e.update_effect_callback(e, spawn_controls);

		return spawn_controls;
	}

	void ControlParticles2d(tfxParticleManager &pm, tfxEffectEmitter &e, tfxU32 amount_spawned) {
		tfxPROFILE;

		tfxEmitterProperties &properties = e.GetProperties();

		//-------------------------------------------------------
		//Controll what the particle does over the course of
		//it's lifetime
		//-------------------------------------------------------

		tfxControlData c;
		c.flags = e.flags;
		c.velocity_adjuster = e.current.velocity_adjuster;
		c.global_intensity = e.current.intensity;
		c.image_size_y = properties.image->image_size.y;
		c.image_frame_rate = properties.image->animation_frames > 1 && e.common.property_flags & tfxEmitterPropertyFlags_animate ? properties.frame_rate : 0.f;
		c.end_frame = properties.end_frame;
		c.stretch = e.current.stretch;
		c.emitter_size_y = e.current.emitter_size.y;
		c.overal_scale = e.current.overal_scale;
		c.angle_offset = properties.angle_offsets.roll;
		c.graphs = &e.common.library->emitter_attributes[e.emitter_attributes].overtime;
		c.particle_update_callback = e.particle_update_callback;
		c.user_data = e.user_data;

		tfxring<tfxParticle> &bank = pm.particle_banks[e.particles_index];
		int offset = 0;
		int particles_index = e.particles_index;
		tfxU32 amount_to_update = bank.current_size;
		if (e.sprites_count < bank.current_size) {
			amount_to_update = e.sprites_count;
		}
		for (int i = amount_to_update - 1 - amount_spawned; i >= 0; --i) {
			tfxParticle *p = &bank[i];

			//Before we do anything, see if the particle should be removed/end of life
			p->data.age += tfxFRAME_LENGTH;
			p->data.flags |= e.flags & tfxParticleFlags_remove;

			if (p->data.flags & tfxParticleFlags_remove || p->data.age >= p->data.max_age) {
				if (e.common.property_flags & tfxEmitterPropertyFlags_single && !(pm.flags & tfxEffectManagerFlags_disable_spawning))
					if (p->data.single_loop_count++ != properties.single_shot_limit) {
						p->data.age = 0;
					}
					else {
						p->data.flags |= tfxParticleFlags_remove;
						p->next_ptr = nullptr;
					}
				else {
					p->data.flags |= tfxParticleFlags_remove;
					p->next_ptr = nullptr;
				}
			}

			if (p->data.flags & tfxParticleFlags_remove) {
				//Draw a blank sprite for one frame so there's no gaps with invalid data
				assert(e.sprites_index < pm.sprites2d[properties.layer].capacity);
				tfxParticleSprite2d &s = pm.sprites2d[properties.layer][e.sprites_index++];
				s.transform.scale = 0;
				s.color.a = 0;
				s.image_ptr = properties.image->ptr;
				s.image_frame = 0;
				offset++;
				continue;
			}
			else if (offset > 0) {
				tfxParticle *next_ptr = bank[i + offset].next_ptr;
				bank[i + offset] = bank[i];
				bank[i].next_ptr = &bank[i + offset];
				p = &bank[i + offset];
				p->next_ptr = next_ptr;
			}

			assert(e.sprites_index < pm.sprites2d[properties.layer].capacity);
			p->sprite_index = (properties.layer << 28) + e.sprites_index;
			tfxParticleSprite2d &s = pm.sprites2d[properties.layer][e.sprites_index++];

			UpdateParticle2d(p->data, s.transform.scale, c);

			if (p->data.flags & tfxParticleFlags_capture_after_transform) {
				e.current.transform_particle_callback2d(p->data, s.transform.position, s.transform.rotation, e.common, e.common.transform.captured_position);
				p->data.captured_position = s.transform.position;
				e.current.transform_particle_callback2d(p->data, s.transform.position, s.transform.rotation, e.common, e.common.transform.world_position);
				p->data.flags &= ~tfxParticleFlags_capture_after_transform;
			}
			else {
				e.current.transform_particle_callback2d(p->data, s.transform.position, s.transform.rotation, e.common, e.common.transform.world_position);
			}

			p->data.flags &= ~tfxParticleFlags_fresh;

			s.color = p->data.color;
			s.image_ptr = properties.image->ptr;
			s.image_frame = (tfxU32)p->data.image_frame;
			s.transform.captured_position = p->data.captured_position.xy();
			s.intensity = p->data.intensity;
			s.handle = e.current.image_handle;

			p->data.captured_position = s.transform.position;

		}
		if (offset) {
			bank.bump(offset);
			//pm.sprites2d->current_size -= offset; //I put this in for a reason to fix another bug, but it causes current size to wrap around with some effects.
		}
	}

	void ControlParticles3d(tfxParticleManager &pm, tfxEffectEmitter &e, tfxU32 amount_spawned) {
		tfxPROFILE;

		tfxEmitterProperties &properties = e.GetProperties();

		//-------------------------------------------------------
		//Controll what the particle does over the course of
		//it's lifetime
		//-------------------------------------------------------

		tfxControlData c;
		c.flags = e.flags;
		c.velocity_adjuster = e.current.velocity_adjuster;
		c.global_intensity = e.current.intensity;
		c.image_size_y = properties.image->image_size.y;
		c.image_frame_rate = properties.image->animation_frames > 1 && e.common.property_flags & tfxEmitterPropertyFlags_animate ? properties.frame_rate : 0.f;
		c.end_frame = properties.end_frame;
		c.stretch = e.current.stretch;
		c.emitter_size_y = e.current.emitter_size.y;
		c.overal_scale = e.current.overal_scale;
		c.angle_offset = properties.angle_offsets.roll;
		c.graphs = &e.common.library->emitter_attributes[e.emitter_attributes].overtime;
		c.particle_update_callback = e.particle_update_callback;
		c.user_data = e.user_data;

		tfxring<tfxParticle> &bank = pm.particle_banks[e.particles_index];
		int offset = 0;
		int particles_index = e.particles_index;
		tfxU32 amount_to_update = bank.current_size;
		if (e.sprites_count < bank.current_size) {
			amount_to_update = e.sprites_count;
		}
		for (int i = amount_to_update - 1 - amount_spawned; i >= 0; --i) {
			tfxParticle *p = &bank[i];
			p->next_ptr = p;

			//Before we do anything, see if the particle should be removed/end of life
			p->data.age += tfxFRAME_LENGTH;
			p->data.flags |= e.flags & tfxParticleFlags_remove;

			if (p->data.flags & tfxParticleFlags_remove || p->data.age >= p->data.max_age) {
				if (e.common.property_flags & tfxEmitterPropertyFlags_single && !(pm.flags & tfxEffectManagerFlags_disable_spawning))
					if (p->data.single_loop_count++ != properties.single_shot_limit) {
						p->data.age = 0;
					}
					else {
						p->data.flags |= tfxParticleFlags_remove;
						p->next_ptr = nullptr;
					}
				else {
					p->data.flags |= tfxParticleFlags_remove;
					p->next_ptr = nullptr;
				}
			}
			if (p->data.flags & tfxParticleFlags_remove) {
				//Draw a blank sprite for one frame so there's no gaps with invalid data
				tfxParticleSprite3d &s = pm.sprites3d[properties.layer][e.sprites_index++];
				s.transform.scale = 0;
				s.color.a = 0;
				s.image_ptr = properties.image->ptr;
				s.image_frame_plus = 0;
				offset++;
				continue;
			}
			else if (offset > 0) {
				tfxParticle *next_ptr = bank[i + offset].next_ptr;
				bank[i + offset] = bank[i];
				bank[i].next_ptr = &bank[i + offset];
				p = &bank[i + offset];
				p->next_ptr = next_ptr;
			}
			p->sprite_index = (properties.layer << 28) + e.sprites_index;
			tfxParticleSprite3d &s = pm.sprites3d[properties.layer][e.sprites_index++];

			UpdateParticle3d(p->data, s.transform.scale, c);

			if (p->data.flags & tfxParticleFlags_capture_after_transform) {
				e.current.transform_particle_callback3d(p->data, s.transform.position, s.transform.rotations, e.common, e.common.transform.captured_position);
				p->data.captured_position = s.transform.position;
				e.current.transform_particle_callback3d(p->data, s.transform.position, s.transform.rotations, e.common, e.common.transform.world_position);
				p->data.flags &= ~tfxParticleFlags_capture_after_transform;
			}
			else {
				e.current.transform_particle_callback3d(p->data, s.transform.position, s.transform.rotations, e.common, e.common.transform.world_position);
			}

			tfxVec3 alignment_vector;
			bool line = e.common.property_flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type == tfxEmissionType::tfxLine;
			if (properties.vector_align_type == tfxVectorAlignType_motion) {
				alignment_vector = s.transform.position - p->data.captured_position;
				float l = FastLength(alignment_vector);
				p->data.velocity_normal.w *= l * 10.f;
				alignment_vector = FastNormalizeVec(alignment_vector);
			}
			else if (properties.vector_align_type == tfxVectorAlignType_emission) {
				alignment_vector = p->data.velocity_normal.xyz();
			}
			else if (properties.vector_align_type == tfxVectorAlignType_emitter) {
				alignment_vector = mmTransformVector(e.common.transform.matrix, tfxVec4(0.f, 1.f, 0.f, 0.f)).xyz();
			}

			s.color = p->data.color;
			s.image_frame_plus = (properties.billboard_option << 24) + (tfxU32)p->data.image_frame;
			s.transform.captured_position = p->data.captured_position;
			s.stretch = p->data.velocity_normal.w;
			s.intensity = p->data.intensity;
			s.handle = e.current.image_handle;
			alignment_vector.y += 0.002f;	//We don't want a 0 alignment normal
			s.alignment = Pack10bit(alignment_vector, properties.billboard_option & 0x00000003);
			s.image_ptr = properties.image->ptr;

			p->data.captured_position = s.transform.position;
		}
		if (offset) {
			bank.bump(offset);
			//pm.sprites3d->current_size -= offset; //I put this in for a reason to fix another bug, but it causes current size to wrap around with some effects.
		}
	}

	bool ControlParticle(tfxParticleManager &pm, tfxParticle &p, tfxVec2 &sprite_scale, tfxEffectEmitter &e) {
		tfxPROFILE;
		if (pm.flags & tfxEffectManagerFlags_update_base_values)
			ReloadBaseValues(p, e);

		p.data.age += tfxFRAME_LENGTH;
		tfxEmitterProperties &properties = e.GetProperties();

		//-------------------------------------------------------
		//Controll what the particle does over the course of
		//it's lifetime
		//-------------------------------------------------------

		//Before we do anything, see if the particle should be removed/end of life
		p.data.flags |= e.flags & tfxParticleFlags_remove;
		if (p.data.flags & tfxParticleFlags_remove)
			return false;

		if (p.data.age >= p.data.max_age) {
			if (e.common.property_flags & tfxEmitterPropertyFlags_single && !(pm.flags & tfxEffectManagerFlags_disable_spawning))
				if (p.data.single_loop_count++ != properties.single_shot_limit) {
					p.data.age = 0;
				}
				else {
					return false;
				}
			else {
				return false;
			}
		}

		tfxControlData c;
		c.flags = e.flags;
		c.velocity_adjuster = e.current.velocity_adjuster;
		c.global_intensity = e.current.intensity;
		c.image_size_y = properties.image->image_size.y;
		c.image_frame_rate = properties.image->animation_frames > 1 && e.common.property_flags & tfxEmitterPropertyFlags_animate ? properties.frame_rate : 0.f;
		c.end_frame = properties.end_frame;
		c.stretch = e.current.stretch;
		c.emitter_size_y = e.current.emitter_size.y;
		c.overal_scale = e.current.overal_scale;
		c.angle_offset = properties.angle_offsets.roll;
		c.graphs = &e.common.library->emitter_attributes[e.emitter_attributes].overtime;
		c.particle_update_callback = e.particle_update_callback;
		c.user_data = e.user_data;

		if (e.common.property_flags & tfxEmitterPropertyFlags_is_3d)
			UpdateParticle3d(p.data, sprite_scale, c);
		else
			UpdateParticle2d(p.data, sprite_scale, c);

		return p.data.flags & tfxParticleFlags_remove ? false : true;
	}

	void TransformEffector(tfxEffectEmitter &e, tfxSpriteTransform2d &parent, bool relative_position, bool relative_angle) {

		if (relative_position) {
			e.common.transform.world_rotations.roll = parent.rotation + e.common.transform.local_rotations.roll;
			e.common.transform.world_position = parent.position;
		}
		else {
			e.common.transform.world_position = e.common.transform.local_position;
			e.common.transform.world_rotations.roll = e.common.transform.local_rotations.roll;
		}

		float s = sin(e.common.transform.world_rotations.roll);
		float c = cos(e.common.transform.world_rotations.roll);
		e.common.transform.matrix.Set2(c, s, -s, c);

	}

	void TransformEffector3d(tfxEffectEmitter &e, tfxSpriteTransform3d &parent, bool relative_position, bool relative_angle) {

		if (relative_position) {
			e.common.transform.world_rotations = parent.rotations + e.common.transform.local_rotations;
			e.common.transform.world_position = parent.position;
		}
		else {
			e.common.transform.world_position = e.common.transform.local_position;
			e.common.transform.world_rotations = e.common.transform.local_rotations;
		}

		tfxMatrix4 roll = mmZRotate(e.common.transform.world_rotations.roll);
		tfxMatrix4 pitch = mmXRotate(e.common.transform.world_rotations.pitch);
		tfxMatrix4 yaw = mmYRotate(e.common.transform.world_rotations.yaw);

		e.common.transform.matrix = mmTransform(yaw, pitch);
		e.common.transform.matrix = mmTransform(e.common.transform.matrix, roll);

	}

	tfxU32 tfxCURRENT_PROFILE_OFFSET = 0;
	const tfxU32 tfxPROFILE_COUNT = __COUNTER__;
	tfxProfile tfxPROFILE_ARRAY[tfxPROFILE_COUNT];
	tfxMemoryTrackerLog tfxMEMORY_TRACKER;
	char tfxMEMORY_CONTEXT[64];
	tfxDataTypesDictionary data_types;
	tfxMemoryArenaManager tfxSTACK_ALLOCATOR;

	void InitialiseTimelineFX() {
		tfxSTACK_ALLOCATOR = CreateArenaManager(tfxSTACK_SIZE, 8);
	}

	bool EndOfProfiles() {
		assert(tfxPROFILE_COUNT);	//there must be tfxPROFILE used in the code
		if (tfxCURRENT_PROFILE_OFFSET == tfxPROFILE_COUNT) {
			tfxCURRENT_PROFILE_OFFSET = 0;
			memset(tfxPROFILE_ARRAY, 0, sizeof(tfxPROFILE_ARRAY));
			return true;
		}
		return false;
	}
	tfxProfile* NextProfile() {
		assert(tfxPROFILE_COUNT && tfxCURRENT_PROFILE_OFFSET < tfxPROFILE_COUNT);	//there must be tfxPROFILE used in the code
		tfxProfile *profile = tfxPROFILE_ARRAY + tfxCURRENT_PROFILE_OFFSET++;
		return profile;
	}

}