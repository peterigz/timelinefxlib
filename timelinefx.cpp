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
		__m128i tmp = _mm_slli_epi32(ji, 23);
		__m128 j = *(__m128*)&tmp; //create vector 1.0f
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

	tfxvec<tfxVec3> tfxIcospherePoints[6];

	void MakeIcospheres() {
		const float x = .525731112119133606f;
		const float z = .850650808352039932f;
		const float n = 0.f;
		const int subdivisions = 6;

		tfxIcospherePoints[0].reserve(42);
		tfxIcospherePoints[1].reserve(162);
		tfxIcospherePoints[2].reserve(642);
		tfxIcospherePoints[3].reserve(2562);
		tfxIcospherePoints[4].reserve(10242);
		tfxIcospherePoints[5].reserve(40962);

		tfxvec<tfxVec3> vertices;
		vertices.push_back({-x, n, z});
		vertices.push_back({ x, n, z});
		vertices.push_back({-x, n,-z});
		vertices.push_back({ x, n,-z});
		vertices.push_back({ n, z, x});
		vertices.push_back({ n, z,-x});
		vertices.push_back({ n,-z, x});
		vertices.push_back({ n,-z,-x});
		vertices.push_back({ z, x, n});
		vertices.push_back({-z, x, n});
		vertices.push_back({ z,-x, n});
		vertices.push_back({-z,-x, n});

		tfxvec<tfxFace> triangles;
		triangles.push_back({0,  4, 1});
		triangles.push_back({0 , 9, 4});
		triangles.push_back({9 , 5, 4});
		triangles.push_back({4 , 5, 8});
		triangles.push_back({4 , 8, 1});
		triangles.push_back({8 ,10, 1});
		triangles.push_back({8 , 3,10});
		triangles.push_back({5 , 3, 8});
		triangles.push_back({5 , 2, 3});
		triangles.push_back({2 , 7, 3});
		triangles.push_back({7 ,10, 3});
		triangles.push_back({7 , 6,10});
		triangles.push_back({7 ,11, 6});
		triangles.push_back({11, 0, 6});
		triangles.push_back({0 , 1, 6});
		triangles.push_back({6 , 1,10});
		triangles.push_back({9 , 0,11});
		triangles.push_back({9 ,11, 2});
		triangles.push_back({9 , 2, 5});
		triangles.push_back({7 , 2,11});

		tfxStorageMap<int> point_cache;

		for (int i = 0; i < subdivisions; ++i)
		{
			triangles = SubDivideIcosphere(point_cache, vertices, triangles);
			assert(tfxIcospherePoints[i].capacity == vertices.current_size);	//Must be the same size
			memcpy(tfxIcospherePoints[i].data, vertices.data, vertices.current_size * sizeof(tfxVec3));
			tfxIcospherePoints[i].current_size = vertices.current_size;
			std::qsort(tfxIcospherePoints[i].data, tfxIcospherePoints[i].current_size, sizeof(tfxVec3), SortIcospherePoints);
		}

	}

	int VertexForEdge(tfxStorageMap<int> &point_cache, tfxvec<tfxVec3>& vertices, int first, int second)
	{
		tfxKey key = ((tfxKey)first << 32) + second;
		if (first > second)
			key = ((tfxKey)second << 32) + first;

		if (point_cache.ValidKey(key))
			return point_cache.At(key);

		point_cache.Insert(key, vertices.size());

		auto& edge0 = vertices[first];
		auto& edge1 = vertices[second];
		auto point = NormalizeVec(edge0 + edge1);
		vertices.push_back(point);

		return vertices.size() - 1;
	}

	tfxvec<tfxFace> SubDivideIcosphere(tfxStorageMap<int> &point_cache, tfxvec<tfxVec3> &vertices, tfxvec<tfxFace> &triangles)
	{
		tfxvec<tfxFace> result;

		for (auto&& each : triangles)
		{
			int mid[3];
			for (int edge = 0; edge < 3; ++edge)
			{
				mid[edge] = VertexForEdge(point_cache, vertices, each.v[edge], each.v[(edge + 1) % 3]);
			}

			result.push_back({ each.v[0], mid[0], mid[2] });
			result.push_back({ each.v[1], mid[1], mid[0] });
			result.push_back({ each.v[2], mid[2], mid[1] });
			result.push_back({ mid[0], mid[1], mid[2] });
		}

		return result;
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

	void SoftExpire(tfxParticleManager &pm, tfxU32 effect_id) {
		pm.emitters.state_flags[effect_id] |= tfxEmitterStateFlags_stop_spawning;
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
			if (e.property_flags & tfxEmitterPropertyFlags_single)
				return true;
			float qty = e.library->emitter_attributes[e.emitter_attributes].base.amount.GetLastValue() + e.library->emitter_attributes[e.emitter_attributes].variation.amount.GetLastValue();
			if (!(e.property_flags & tfxEmitterPropertyFlags_single) && qty > 0)
				return false;
		}
		return true;
	}

	void tfxEffectEmitter::FlagAs3D(bool flag) {
		if(flag)
			property_flags |= tfxEmitterPropertyFlags_is_3d;
		else
			property_flags &= ~tfxEmitterPropertyFlags_is_3d;
		for (auto &sub : GetInfo().sub_effectors) {
			sub.FlagAs3D(flag);
		}
	}

	bool tfxEffectEmitter::Is3DEffect() {
		return property_flags & tfxEmitterPropertyFlags_is_3d;
	}

	tfxParticleManagerModes tfxEffectEmitter::GetRequiredParticleManagerMode() {
		if (type == tfxEffectType) {
			if (effect_flags & tfxEffectPropertyFlags_guaranteed_order && effect_flags & tfxEffectPropertyFlags_depth_draw_order) {
				return tfxParticleManagerMode_ordered_by_depth_guaranteed;
			}
			else if (effect_flags & tfxEffectPropertyFlags_depth_draw_order) {
				return tfxParticleManagerMode_ordered_by_depth;
			}
			else if (effect_flags & tfxEffectPropertyFlags_age_order) {
				return tfxParticleManagerMode_ordered_by_age;
			}
		}
		else if (type == tfxStage) {
			tfxParticleManagerModes result = tfxParticleManagerMode_unordered;
			for (auto &effect : GetInfo().sub_effectors) {
				if (effect_flags & tfxEffectPropertyFlags_guaranteed_order && effect_flags & tfxEffectPropertyFlags_depth_draw_order) {
					return tfxParticleManagerMode_ordered_by_depth_guaranteed;
				}
				else if (effect_flags & tfxEffectPropertyFlags_depth_draw_order) {
					result = tfxParticleManagerMode_ordered_by_depth;
				}
				else if (result != tfxParticleManagerMode_ordered_by_depth && effect_flags & tfxEffectPropertyFlags_age_order) {
					result = tfxParticleManagerMode_ordered_by_age;
				}
			}
			return result;
		}
		return tfxParticleManagerMode_unordered;
	}

	tfxPreviewCameraSettings &tfxEffectEmitter::GetCameraSettings() {
		return library->preview_camera_settings[GetInfo().preview_camera_settings];
	}

	tfxEffectEmitter& tfxEffectEmitter::AddEmitter(tfxEffectEmitter &e) {
		assert(e.GetInfo().name.Length());				//Emitter must have a name so that a hash can be generated
		e.type = tfxEffectEmitterType::tfxEmitterType;
		e.library = library;
		e.GetInfo().uid = ++library->uid;
		GetInfo().sub_effectors.push_back(e);
		library->UpdateEffectPaths();
		ReIndex();
		return GetInfo().sub_effectors.back();
	}

	tfxEffectEmitter& tfxEffectEmitter::AddEffect(tfxEffectEmitter &e) {
		assert(e.GetInfo().name.Length());				//Effect must have a name so that a hash can be generated
		e.type = tfxEffectEmitterType::tfxEffectType;
		e.library = library;
		e.parent = this;
		e.GetInfo().uid = ++library->uid;
		GetInfo().sub_effectors.push_back(e);
		library->UpdateEffectPaths();
		ReIndex();
		return GetInfo().sub_effectors.back();
	}

	tfxEffectEmitter& tfxEffectEmitter::AddEffect() {
		tfxEffectEmitter e;
		e.library = library;
		e.GetInfo().uid = ++library->uid;
		e.type = tfxEffectEmitterType::tfxEffectType;
		e.GetInfo().name = "New Effect";
		GetInfo().sub_effectors.push_back(e);
		library->UpdateEffectPaths();
		ReIndex();
		return GetInfo().sub_effectors.back();
	}

	tfxEffectEmitter &tfxEffectEmitter::AddEffector(tfxEffectEmitterType type) {
		tfxEffectEmitter e;
		//e.parent_effect = this;
		e.type = type;
		e.library = library;
		e.GetInfo().uid = ++library->uid;
		if(e.type == tfxEffectType)
			e.GetInfo().name = "New Effect";
		else
			e.GetInfo().name = "New Emitter";
		GetInfo().sub_effectors.push_back(e);
		library->UpdateEffectPaths();
		ReIndex();
		return GetInfo().sub_effectors.back();
	}

	float GetEmissionDirection2d(tfxParticleManager &pm, tfxEffectLibrary *library, tfxU32 property_index, tfxU32 emitter_index, tfxVec2 local_position, tfxVec2 world_position, tfxVec2 emitter_size) {
		//float (*effect_lookup_callback)(tfxGraph &graph, float age) = common.root_effect->lookup_mode == tfxPrecise ? LookupPrecise : LookupFast;
		const float frame = pm.emitters.frame[emitter_index];
		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[emitter_index];
		float emission_angle = lookup_callback(library->emitter_attributes[emitter_attributes].properties.emission_pitch, frame);
		float emission_angle_variation = lookup_callback(library->emitter_attributes[emitter_attributes].properties.emission_range, frame);
		//----Emission
		float range = emission_angle_variation *.5f;
		float direction = 0;
		tfxEmissionType emission_type = library->emitter_properties.emission_type[property_index];
		tfxEmissionDirection emission_direction = library->emitter_properties.emission_direction[property_index];

		if (emission_type == tfxPoint)
			return direction + emission_angle + random_generation.Range(-range, range);

		const tfxVec3 &handle = pm.emitters.handle[emitter_index];
		const tfxEmitterPropertyFlags &property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 &emitter_world_position = pm.emitters.world_position[emitter_index];
		float &emission_alternator = pm.emitters.emission_alternator[emitter_index];

		tfxVec2 tmp_position;
		if (handle.x + local_position.x == 0 && handle.y + local_position.y == 0)
			tmp_position = emitter_size;
		else
			tmp_position = local_position + handle.xy();

		if (emission_direction == tfxEmissionDirection::tfxOutwards) {

			tfxVec2 to_handle;

			if (property_flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = tmp_position;
			else
				to_handle = world_position - emitter_world_position.xy();

			direction = GetVectorAngle(to_handle.x, to_handle.y);

		}
		else if (emission_direction == tfxEmissionDirection::tfxInwards) {

			tfxVec2 to_handle;

			if (property_flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = -tmp_position;

			else
				to_handle = emitter_world_position.xy() - world_position;

			direction = GetVectorAngle(to_handle.x, to_handle.y);

		}
		else if (emission_direction == tfxEmissionDirection::tfxBothways) {

			//todo: replace these if statements
			if (emission_alternator) {

				tfxVec2 to_handle;

				if (property_flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = (tmp_position);
				else
					to_handle = world_position - emitter_world_position.xy();

				direction = GetVectorAngle(to_handle.x, to_handle.y);

			}
			else {

				tfxVec2 to_handle;

				if (property_flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = -tmp_position;
				else
					to_handle = (emitter_world_position.xy() - world_position);

				direction = GetVectorAngle(to_handle.x, to_handle.y);

			}

			emission_alternator = !emission_alternator;
		}

		if (std::isnan(direction))
			direction = 0.f;
		return direction + emission_angle + random_generation.Range(-range, range);
	}

	tfxVec3 GetEmissionDirection3d( tfxParticleManager &pm,	tfxEffectLibrary *library, tfxU32 property_index, tfxU32 emitter_index, float emission_pitch, float emission_yaw, tfxVec3 local_position, tfxVec3 world_position, tfxVec3 emitter_size) {
		//float (*effect_lookup_callback)(tfxGraph &graph, float age) = common.root_effect->lookup_mode == tfxPrecise ? LookupPrecise : LookupFast;
		const float frame = pm.emitters.frame[emitter_index];
		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[emitter_index];
		float emission_angle_variation = lookup_callback(library->emitter_attributes[emitter_attributes].properties.emission_range, frame);
		//----Emission
		float range = emission_angle_variation * .5f;

		const tfxVec3 &handle = pm.emitters.handle[emitter_index];
		const tfxVec3 &emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxVec3 &emitter_world_rotations = pm.emitters.world_rotations[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		float &emission_alternator = pm.emitters.emission_alternator[emitter_index];

		tfxVec3 result;
		tfxVec3 tmp_position;
		if (handle.x + local_position.x == 0 && handle.y + local_position.y == 0 && handle.z + local_position.z == 0)
			tmp_position = emitter_size;
		else
			tmp_position = local_position + handle;

		tfxEmissionType emission_type = library->emitter_properties.emission_type[property_index];
		tfxEmissionDirection emission_direction = library->emitter_properties.emission_direction[property_index];

		tfxVec3 to_handle(0.f, 1.f, 0.f);
		float parent_pitch = 0.f;
		float parent_yaw = 0.f;
		if (emission_type != tfxPoint) {
			if (emission_direction == tfxEmissionDirection::tfxOutwards) {

				if (property_flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = tmp_position;
				else
					to_handle = world_position - emitter_world_position;

				to_handle = FastNormalizeVec(to_handle);

			}
			else if (emission_direction == tfxEmissionDirection::tfxInwards) {

				if (property_flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = -tmp_position;
				else
					to_handle = emitter_world_position - world_position;

				to_handle = FastNormalizeVec(to_handle);

			}
			else if (emission_direction == tfxEmissionDirection::tfxBothways) {

				if (emission_alternator) {

					if (property_flags & tfxEmitterPropertyFlags_relative_position)
						to_handle = tmp_position;
					else
						to_handle = world_position - emitter_world_position;
				}
				else {

					if (property_flags & tfxEmitterPropertyFlags_relative_position)
						to_handle = -tmp_position;
					else
						to_handle = emitter_world_position - world_position;
				}

				emission_alternator = !emission_alternator;
				to_handle = FastNormalizeVec(to_handle);
			}
			else {
				parent_pitch = emitter_world_rotations.pitch;
				parent_yaw = emitter_world_rotations.yaw;
			}
		}
		else {
			parent_pitch = emitter_world_rotations.pitch;
			parent_yaw = emitter_world_rotations.yaw;
		}

		float pitch = asinf(-to_handle.y);
		float yaw = atan2(to_handle.x, to_handle.z);
		tfxVec3 direction;
		direction.z = cos(emission_yaw + yaw) * cos(emission_pitch + pitch);
		direction.y = -sin(emission_pitch + pitch);
		direction.x = sin(emission_yaw + yaw) * cos(emission_pitch + pitch);
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

	void ReloadBaseValues(tfxParticle &p, tfxEffectEmitter &e) {
		//Todo: let's just rethink this whole function
	}

	void tfxEffectEmitter::ResetGlobalGraphs(bool add_node, bool compile) {
		library->global_graphs[global].life.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].life.type = tfxGlobal_life;
		library->global_graphs[global].amount.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].amount.type = tfxGlobal_amount;
		library->global_graphs[global].velocity.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].velocity.type = tfxGlobal_velocity;
		library->global_graphs[global].width.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].width.type = tfxGlobal_width;
		library->global_graphs[global].height.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].height.type = tfxGlobal_height;
		library->global_graphs[global].weight.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].weight.type = tfxGlobal_weight;
		library->global_graphs[global].spin.Reset(1.f, tfxGlobalPercentPresetSigned, add_node); library->global_graphs[global].spin.type = tfxGlobal_spin;
		library->global_graphs[global].stretch.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].stretch.type = tfxGlobal_stretch;
		library->global_graphs[global].overal_scale.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].overal_scale.type = tfxGlobal_overal_scale;
		library->global_graphs[global].intensity.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].intensity.type = tfxGlobal_intensity;
		library->global_graphs[global].frame_rate.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].frame_rate.type = tfxGlobal_frame_rate;
		library->global_graphs[global].splatter.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].splatter.type = tfxGlobal_splatter;
		library->global_graphs[global].emitter_width.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].emitter_width.type = tfxGlobal_emitter_width;
		library->global_graphs[global].emitter_height.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].emitter_height.type = tfxGlobal_emitter_height;
		library->global_graphs[global].emitter_depth.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].emitter_depth.type = tfxGlobal_emitter_depth;
		if(compile)
			library->CompileGlobalGraph(global);
	}

	void tfxEffectEmitter::ResetTransformGraphs(bool add_node, bool compile) {
		library->transform_attributes[transform_attributes].roll.Reset(0.f, tfxAnglePreset, add_node); library->transform_attributes[transform_attributes].roll.type = tfxTransform_roll;
		library->transform_attributes[transform_attributes].pitch.Reset(0.f, tfxAnglePreset, add_node); library->transform_attributes[transform_attributes].pitch.type = tfxTransform_pitch;
		library->transform_attributes[transform_attributes].yaw.Reset(0.f, tfxAnglePreset, add_node); library->transform_attributes[transform_attributes].yaw.type = tfxTransform_yaw;
		library->transform_attributes[transform_attributes].translation_x.Reset(0.f, tfxTranslationPreset, add_node); library->transform_attributes[transform_attributes].translation_x.type = tfxTransform_translate_x;
		library->transform_attributes[transform_attributes].translation_y.Reset(0.f, tfxTranslationPreset, add_node); library->transform_attributes[transform_attributes].translation_y.type = tfxTransform_translate_y;
		library->transform_attributes[transform_attributes].translation_z.Reset(0.f, tfxTranslationPreset, add_node); library->transform_attributes[transform_attributes].translation_z.type = tfxTransform_translate_z;
		if(compile)
			library->CompileKeyframeGraph(transform_attributes);
	}

	void tfxEffectEmitter::ResetBaseGraphs(bool add_node, bool compile) {
		library->emitter_attributes[emitter_attributes].base.life.Reset(1000.f, tfxLifePreset, add_node); library->emitter_attributes[emitter_attributes].base.life.type = tfxBase_life;
		library->emitter_attributes[emitter_attributes].base.amount.Reset(1.f, tfxAmountPreset, add_node); library->emitter_attributes[emitter_attributes].base.amount.type = tfxBase_amount;
		library->emitter_attributes[emitter_attributes].base.velocity.Reset(0.f, tfxVelocityPreset, add_node); library->emitter_attributes[emitter_attributes].base.velocity.type = tfxBase_velocity;
		library->emitter_attributes[emitter_attributes].base.width.Reset(128.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].base.width.type = tfxBase_width;
		library->emitter_attributes[emitter_attributes].base.height.Reset(128.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].base.height.type = tfxBase_height;
		library->emitter_attributes[emitter_attributes].base.weight.Reset(0.f, tfxWeightPreset, add_node); library->emitter_attributes[emitter_attributes].base.weight.type = tfxBase_weight;
		library->emitter_attributes[emitter_attributes].base.spin.Reset(0.f, tfxSpinPreset, add_node); library->emitter_attributes[emitter_attributes].base.spin.type = tfxBase_spin;
		library->emitter_attributes[emitter_attributes].base.noise_offset.Reset(0.f, tfxGlobalPercentPreset, add_node); library->emitter_attributes[emitter_attributes].base.noise_offset.type = tfxBase_noise_offset;
		if(compile)
			library->CompileBaseGraph(emitter_attributes);
	}

	void tfxEffectEmitter::ResetPropertyGraphs(bool add_node, bool compile) {
		library->emitter_attributes[emitter_attributes].properties.emission_pitch.Reset(0.f, tfxAnglePreset, add_node); library->emitter_attributes[emitter_attributes].properties.emission_pitch.type = tfxProperty_emission_pitch;
		library->emitter_attributes[emitter_attributes].properties.emission_yaw.Reset(0.f, tfxAnglePreset, add_node); library->emitter_attributes[emitter_attributes].properties.emission_yaw.type = tfxProperty_emission_yaw;
		library->emitter_attributes[emitter_attributes].properties.emission_range.Reset(0.f, tfxEmissionRangePreset, add_node); library->emitter_attributes[emitter_attributes].properties.emission_range.type = tfxProperty_emission_range;
		library->emitter_attributes[emitter_attributes].properties.splatter.Reset(0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].properties.splatter.type = tfxProperty_splatter;
		library->emitter_attributes[emitter_attributes].properties.emitter_width.Reset(0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].properties.emitter_width.type = tfxProperty_emitter_width;
		library->emitter_attributes[emitter_attributes].properties.emitter_height.Reset(0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].properties.emitter_height.type = tfxProperty_emitter_height;
		library->emitter_attributes[emitter_attributes].properties.emitter_depth.Reset(0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].properties.emitter_depth.type = tfxProperty_emitter_depth;
		library->emitter_attributes[emitter_attributes].properties.arc_size.Reset(tfxRadians(360.f), tfxArcPreset, add_node); library->emitter_attributes[emitter_attributes].properties.arc_size.type = tfxProperty_arc_size;
		library->emitter_attributes[emitter_attributes].properties.arc_offset.Reset(0.f, tfxArcPreset, add_node); library->emitter_attributes[emitter_attributes].properties.arc_offset.type = tfxProperty_arc_offset;
		if(compile)
			library->CompilePropertyGraph(emitter_attributes);
	}

	void tfxEffectEmitter::ResetVariationGraphs(bool add_node, bool compile) {
		library->emitter_attributes[emitter_attributes].variation.life.Reset(0.f, tfxLifePreset, add_node); library->emitter_attributes[emitter_attributes].variation.life.type = tfxVariation_life;
		library->emitter_attributes[emitter_attributes].variation.amount.Reset(0.f, tfxAmountPreset, add_node); library->emitter_attributes[emitter_attributes].variation.amount.type = tfxVariation_amount;
		library->emitter_attributes[emitter_attributes].variation.velocity.Reset(0.f, tfxVelocityPreset, add_node); library->emitter_attributes[emitter_attributes].variation.velocity.type = tfxVariation_velocity;
		library->emitter_attributes[emitter_attributes].variation.width.Reset(0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].variation.width.type = tfxVariation_width;
		library->emitter_attributes[emitter_attributes].variation.height.Reset(0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].variation.height.type = tfxVariation_height;
		library->emitter_attributes[emitter_attributes].variation.weight.Reset(0.f, tfxWeightVariationPreset, add_node); library->emitter_attributes[emitter_attributes].variation.weight.type = tfxVariation_weight;
		library->emitter_attributes[emitter_attributes].variation.spin.Reset(0.f, tfxSpinVariationPreset, add_node); library->emitter_attributes[emitter_attributes].variation.spin.type = tfxVariation_spin;
		library->emitter_attributes[emitter_attributes].variation.noise_offset.Reset(0.f, tfxNoiseOffsetVariationPreset, add_node); library->emitter_attributes[emitter_attributes].variation.noise_offset.type = tfxVariation_noise_offset;
		library->emitter_attributes[emitter_attributes].variation.noise_resolution.Reset(300.f, tfxNoiseResolutionPreset, add_node); library->emitter_attributes[emitter_attributes].variation.noise_resolution.type = tfxVariation_noise_resolution;
		if(compile)
			library->CompileVariationGraph(emitter_attributes);
	}

	void tfxEffectEmitter::ResetOvertimeGraphs(bool add_node, bool compile) {
		library->emitter_attributes[emitter_attributes].overtime.velocity.Reset(1.f, tfxVelocityOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.velocity.type = tfxOvertime_velocity;
		library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.Reset(1.f, tfxGlobalPercentPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.type = tfxOvertime_velocity_adjuster;
		library->emitter_attributes[emitter_attributes].overtime.width.Reset(1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.width.type = tfxOvertime_width;
		library->emitter_attributes[emitter_attributes].overtime.height.Reset(1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.height.type = tfxOvertime_height;
		library->emitter_attributes[emitter_attributes].overtime.weight.Reset(1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.weight.type = tfxOvertime_weight;
		library->emitter_attributes[emitter_attributes].overtime.spin.Reset(0.f, tfxSpinOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.spin.type = tfxOvertime_spin;
		library->emitter_attributes[emitter_attributes].overtime.stretch.Reset(0.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.stretch.type = tfxOvertime_stretch;
		library->emitter_attributes[emitter_attributes].overtime.red.Reset(1.f, tfxColorPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.red.type = tfxOvertime_red;
		library->emitter_attributes[emitter_attributes].overtime.green.Reset(1.f, tfxColorPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.green.type = tfxOvertime_green;
		library->emitter_attributes[emitter_attributes].overtime.blue.Reset(1.f, tfxColorPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.blue.type = tfxOvertime_blue;
		library->emitter_attributes[emitter_attributes].overtime.blendfactor.Reset(1.f, tfxOpacityOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.blendfactor.type = tfxOvertime_blendfactor;
		library->emitter_attributes[emitter_attributes].overtime.intensity.Reset(1.f, tfxIntensityOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.intensity.type = tfxOvertime_intensity;
		library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance.Reset(30.f, tfxVelocityTurbulancePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance.type = tfxOvertime_velocity_turbulance;
		library->emitter_attributes[emitter_attributes].overtime.stretch.Reset(0.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.stretch.type = tfxOvertime_stretch;
		library->emitter_attributes[emitter_attributes].overtime.direction_turbulance.Reset(0.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.direction_turbulance.type = tfxOvertime_direction_turbulance;
		library->emitter_attributes[emitter_attributes].overtime.direction.Reset(0.f, tfxDirectionOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.direction.type = tfxOvertime_direction;
		library->emitter_attributes[emitter_attributes].overtime.noise_resolution.Reset(1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.noise_resolution.type = tfxOvertime_noise_resolution;
		if(compile)
			library->CompileOvertimeGraph(emitter_attributes);
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
		if (library->transform_attributes[transform_attributes].translation_x.nodes.size() == 0) library->transform_attributes[transform_attributes].translation_x.Reset(0.f, tfxTranslationPreset);
		if (library->transform_attributes[transform_attributes].translation_y.nodes.size() == 0) library->transform_attributes[transform_attributes].translation_y.Reset(0.f, tfxTranslationPreset);
		if (library->transform_attributes[transform_attributes].translation_z.nodes.size() == 0) library->transform_attributes[transform_attributes].translation_z.Reset(0.f, tfxTranslationPreset);
		if (library->transform_attributes[transform_attributes].roll.nodes.size() == 0) library->transform_attributes[transform_attributes].roll.Reset(0.f, tfxAnglePreset);
		if (library->transform_attributes[transform_attributes].pitch.nodes.size() == 0) library->transform_attributes[transform_attributes].pitch.Reset(0.f, tfxAnglePreset);
		if (library->transform_attributes[transform_attributes].yaw.nodes.size() == 0) library->transform_attributes[transform_attributes].yaw.Reset(0.f, tfxAnglePreset);

		if (type == tfxEffectType) {
			if (library->global_graphs[global].life.nodes.size() == 0) library->global_graphs[global].life.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].amount.nodes.size() == 0) library->global_graphs[global].amount.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].velocity.nodes.size() == 0) library->global_graphs[global].velocity.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].width.nodes.size() == 0) library->global_graphs[global].width.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].height.nodes.size() == 0) library->global_graphs[global].height.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].weight.nodes.size() == 0) library->global_graphs[global].weight.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].spin.nodes.size() == 0) library->global_graphs[global].spin.Reset(1.f, tfxGlobalPercentPresetSigned);
			if (library->global_graphs[global].stretch.nodes.size() == 0) library->global_graphs[global].stretch.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].overal_scale.nodes.size() == 0) library->global_graphs[global].overal_scale.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].intensity.nodes.size() == 0) library->global_graphs[global].intensity.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].frame_rate.nodes.size() == 0) library->global_graphs[global].frame_rate.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].splatter.nodes.size() == 0) library->global_graphs[global].splatter.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].emitter_width.nodes.size() == 0) library->global_graphs[global].emitter_width.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].emitter_height.nodes.size() == 0) library->global_graphs[global].emitter_height.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].emitter_depth.nodes.size() == 0) library->global_graphs[global].emitter_depth.Reset(1.f, tfxGlobalPercentPreset);
		}

		if (type == tfxEmitterType) {
			if (library->emitter_attributes[emitter_attributes].base.life.nodes.size() == 0) library->emitter_attributes[emitter_attributes].base.life.Reset(1000.f, tfxLifePreset);
			if (library->emitter_attributes[emitter_attributes].base.amount.nodes.size() == 0) library->emitter_attributes[emitter_attributes].base.amount.Reset(1.f, tfxAmountPreset);
			if (library->emitter_attributes[emitter_attributes].base.velocity.nodes.size() == 0) library->emitter_attributes[emitter_attributes].base.velocity.Reset(0.f, tfxVelocityPreset);
			if (library->emitter_attributes[emitter_attributes].base.width.nodes.size() == 0) library->emitter_attributes[emitter_attributes].base.width.Reset(128.f, tfxDimensionsPreset);
			if (library->emitter_attributes[emitter_attributes].base.height.nodes.size() == 0) library->emitter_attributes[emitter_attributes].base.height.Reset(128.f, tfxDimensionsPreset);
			if (library->emitter_attributes[emitter_attributes].base.weight.nodes.size() == 0) library->emitter_attributes[emitter_attributes].base.weight.Reset(0.f, tfxWeightPreset);
			if (library->emitter_attributes[emitter_attributes].base.spin.nodes.size() == 0) library->emitter_attributes[emitter_attributes].base.spin.Reset(0.f, tfxSpinPreset);
			if (library->emitter_attributes[emitter_attributes].base.noise_offset.nodes.size() == 0) library->emitter_attributes[emitter_attributes].base.noise_offset.Reset(0.f, tfxGlobalPercentPreset);

			if (library->emitter_attributes[emitter_attributes].properties.emission_pitch.nodes.size() == 0) library->emitter_attributes[emitter_attributes].properties.emission_pitch.Reset(0.f, tfxAnglePreset);
			if (library->emitter_attributes[emitter_attributes].properties.emission_yaw.nodes.size() == 0) library->emitter_attributes[emitter_attributes].properties.emission_yaw.Reset(0.f, tfxAnglePreset);
			if (library->emitter_attributes[emitter_attributes].properties.emission_range.nodes.size() == 0) library->emitter_attributes[emitter_attributes].properties.emission_range.Reset(0.f, tfxEmissionRangePreset);
			if (library->emitter_attributes[emitter_attributes].properties.splatter.nodes.size() == 0) library->emitter_attributes[emitter_attributes].properties.splatter.Reset(0.f, tfxDimensionsPreset);
			if (library->emitter_attributes[emitter_attributes].properties.emitter_width.nodes.size() == 0) library->emitter_attributes[emitter_attributes].properties.emitter_width.Reset(0.f, tfxDimensionsPreset);
			if (library->emitter_attributes[emitter_attributes].properties.emitter_height.nodes.size() == 0) library->emitter_attributes[emitter_attributes].properties.emitter_height.Reset(0.f, tfxDimensionsPreset);
			if (library->emitter_attributes[emitter_attributes].properties.emitter_depth.nodes.size() == 0) library->emitter_attributes[emitter_attributes].properties.emitter_depth.Reset(0.f, tfxDimensionsPreset);
			if (library->emitter_attributes[emitter_attributes].properties.arc_size.nodes.size() == 0) library->emitter_attributes[emitter_attributes].properties.arc_size.Reset(tfxRadians(360.f), tfxArcPreset);
			if (library->emitter_attributes[emitter_attributes].properties.arc_offset.nodes.size() == 0) library->emitter_attributes[emitter_attributes].properties.arc_offset.Reset(0.f, tfxArcPreset);

			if (library->emitter_attributes[emitter_attributes].variation.life.nodes.size() == 0) library->emitter_attributes[emitter_attributes].variation.life.Reset(0.f, tfxLifePreset);
			if (library->emitter_attributes[emitter_attributes].variation.amount.nodes.size() == 0) library->emitter_attributes[emitter_attributes].variation.amount.Reset(0.f, tfxAmountPreset);
			if (library->emitter_attributes[emitter_attributes].variation.velocity.nodes.size() == 0) library->emitter_attributes[emitter_attributes].variation.velocity.Reset(0.f, tfxVelocityPreset);
			if (library->emitter_attributes[emitter_attributes].variation.weight.nodes.size() == 0) library->emitter_attributes[emitter_attributes].variation.weight.Reset(0.f, tfxVelocityPreset);
			if (library->emitter_attributes[emitter_attributes].variation.width.nodes.size() == 0) library->emitter_attributes[emitter_attributes].variation.width.Reset(0.f, tfxDimensionsPreset);
			if (library->emitter_attributes[emitter_attributes].variation.height.nodes.size() == 0) library->emitter_attributes[emitter_attributes].variation.height.Reset(0.f, tfxDimensionsPreset);
			if (library->emitter_attributes[emitter_attributes].variation.weight.nodes.size() == 0) library->emitter_attributes[emitter_attributes].variation.weight.Reset(0.f, tfxWeightVariationPreset);
			if (library->emitter_attributes[emitter_attributes].variation.spin.nodes.size() == 0) library->emitter_attributes[emitter_attributes].variation.spin.Reset(0.f, tfxSpinVariationPreset);
			if (library->emitter_attributes[emitter_attributes].variation.noise_offset.nodes.size() == 0) library->emitter_attributes[emitter_attributes].variation.noise_offset.Reset(0.f, tfxNoiseOffsetVariationPreset);
			if (library->emitter_attributes[emitter_attributes].variation.noise_resolution.nodes.size() == 0) library->emitter_attributes[emitter_attributes].variation.noise_resolution.Reset(300.f, tfxNoiseResolutionPreset);

			if (library->emitter_attributes[emitter_attributes].overtime.velocity.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.velocity.Reset(1.f, tfxVelocityOvertimePreset);
			if (library->emitter_attributes[emitter_attributes].overtime.width.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.width.Reset(1.f, tfxPercentOvertime);
			if (library->emitter_attributes[emitter_attributes].overtime.height.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.height.Reset(1.f, tfxPercentOvertime);
			if (library->emitter_attributes[emitter_attributes].overtime.weight.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.weight.Reset(1.f, tfxPercentOvertime);
			if (library->emitter_attributes[emitter_attributes].overtime.spin.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.spin.Reset(1.f, tfxSpinOvertimePreset);
			if (library->emitter_attributes[emitter_attributes].overtime.stretch.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.stretch.Reset(0.f, tfxPercentOvertime);
			if (library->emitter_attributes[emitter_attributes].overtime.red.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.red.Reset(1.f, tfxColorPreset);
			if (library->emitter_attributes[emitter_attributes].overtime.green.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.green.Reset(1.f, tfxColorPreset);
			if (library->emitter_attributes[emitter_attributes].overtime.blue.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.blue.Reset(1.f, tfxColorPreset);
			if (library->emitter_attributes[emitter_attributes].overtime.blendfactor.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.blendfactor.Reset(1.f, tfxOpacityOvertimePreset);
			if (library->emitter_attributes[emitter_attributes].overtime.intensity.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.intensity.Reset(1.f, tfxIntensityOvertimePreset);
			if (library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance.Reset(0.f, tfxFrameratePreset);
			if (library->emitter_attributes[emitter_attributes].overtime.direction_turbulance.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.direction_turbulance.Reset(0.f, tfxPercentOvertime);
			if (library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.Reset(1.f, tfxGlobalPercentPreset);
			if (library->emitter_attributes[emitter_attributes].overtime.direction.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.direction.Reset(0.f, tfxDirectionOvertimePreset);
			if (library->emitter_attributes[emitter_attributes].overtime.noise_resolution.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.noise_resolution.Reset(1.f, tfxPercentOvertime);
		}
	}

	void tfxEffectEmitter::SetName(const char *n) {
		GetInfo().name = n;
	}

	void tfxEffectEmitter::ClearColors() {
		library->emitter_attributes[emitter_attributes].overtime.red.Clear();
		library->emitter_attributes[emitter_attributes].overtime.green.Clear();
		library->emitter_attributes[emitter_attributes].overtime.blue.Clear();
	}

	void tfxEffectEmitter::AddColorOvertime(float frame, tfxRGB color) {
		library->emitter_attributes[emitter_attributes].overtime.red.AddNode(frame, color.r);
		library->emitter_attributes[emitter_attributes].overtime.green.AddNode(frame, color.g);
		library->emitter_attributes[emitter_attributes].overtime.blue.AddNode(frame, color.b);
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

	tfxEffectEmitterInfo &tfxEffectEmitter::GetInfo() {
		return library->GetInfo(*this);
	}

	tfxEmitterPropertiesSoA &tfxEffectEmitter::GetProperties() {
		return library->emitter_properties;
	}

	bool tfxEffectEmitter::HasSingle() {
		for (auto &e : GetInfo().sub_effectors) {
			if (e.property_flags & tfxEmitterPropertyFlags_single)
				return true;
		}
		return false;
	}

	bool tfxEffectEmitter::RenameSubEffector(tfxEffectEmitter &emitter, const char *new_name) {
		if (!NameExists(emitter, new_name) && strlen(new_name) > 0) {
			emitter.SetName(new_name);
			library->UpdateEffectPaths();
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
			emitter.library->UpdateEffectPaths();
			return &GetInfo().sub_effectors[new_index];
		}

		return nullptr;
	}

	tfxEffectEmitter* tfxEffectEmitter::MoveDown(tfxEffectEmitter &emitter) {
		if (emitter.library_index < GetInfo().sub_effectors.size() - 1) {
			tfxU32 new_index = emitter.library_index + 1;
			std::swap(GetInfo().sub_effectors[emitter.library_index], GetInfo().sub_effectors[new_index]);
			ReIndex();
			emitter.library->UpdateEffectPaths();
			return &GetInfo().sub_effectors[new_index];
		}
		return nullptr;
	}

	void tfxEffectEmitter::DeleteEmitter(tfxEffectEmitter *emitter) {
		tfxEffectLibrary *library = emitter->library;
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
					library->FreeGlobal(current.global);
				}
				else if(current.type == tfxEmitterType) {
					library->FreeEmitterAttributes(current.emitter_attributes);
				}
				for (auto &sub : current.GetInfo().sub_effectors) {
					stack.push_back(sub);
				}
				current.GetInfo().sub_effectors.clear();
				library->FreeProperties(current.property_index);
				library->FreeInfo(current.info_index);
			}
		}

		ReIndex();
	}

	void tfxEffectEmitter::Clone(tfxEffectEmitter &clone, tfxEffectEmitter *root_parent, tfxEffectLibrary *destination_library, tfxEffectCloningFlags flags) {
		//tfxU32 size = library->global_graphs[0].amount.lookup.values.capacity;
		clone = *this;
		clone.info_index = clone.library->CloneInfo(info_index, destination_library);
		if (clone.type != tfxFolder) {
			clone.property_index = clone.library->CloneProperties(property_index, destination_library);
		}
		clone.property_flags |= tfxEmitterPropertyFlags_enabled;
		if(!(flags & tfxEffectCloningFlags_keep_user_data))
			clone.user_data = nullptr;
		clone.library = destination_library;
		clone.GetInfo().sub_effectors.clear();

		if (type == tfxEffectType) {
			if (root_parent == &clone) {
				clone.global = flags & tfxEffectCloningFlags_clone_graphs ? library->CloneGlobal(global, destination_library) : clone.global = global;
				clone.transform_attributes = flags & tfxEffectCloningFlags_clone_graphs ? library->CloneKeyframes(transform_attributes, destination_library) : clone.transform_attributes = transform_attributes;
				if (flags & tfxEffectCloningFlags_compile_graphs) {
					clone.library->CompileGlobalGraph(clone.global);
					clone.library->CompileKeyframeGraph(clone.transform_attributes);
				}
			}
			else {
				clone.transform_attributes = flags & tfxEffectCloningFlags_clone_graphs ? library->CloneKeyframes(transform_attributes, destination_library) : clone.transform_attributes = transform_attributes;
				if (flags & tfxEffectCloningFlags_compile_graphs) {
					clone.library->CompileKeyframeGraph(clone.transform_attributes);
				}
				if (!(flags & tfxEffectCloningFlags_force_clone_global)) {
					clone.global = root_parent->global;
				}
				else {
					clone.global = library->CloneGlobal(root_parent->global, destination_library);
					if(flags & tfxEffectCloningFlags_compile_graphs)
						clone.library->CompileGlobalGraph(clone.global);
				}
			}
		}
		else if(type == tfxEmitterType) {
			clone.emitter_attributes = flags & tfxEffectCloningFlags_clone_graphs ? library->CloneEmitterAttributes(emitter_attributes, destination_library) : emitter_attributes;
			clone.transform_attributes = flags & tfxEffectCloningFlags_clone_graphs ? library->CloneKeyframes(transform_attributes, destination_library) : transform_attributes;
			clone.UpdateMaxLife();
			if (flags & tfxEffectCloningFlags_compile_graphs) {
				clone.library->CompileKeyframeGraph(clone.transform_attributes);
				clone.library->CompilePropertyGraph(clone.emitter_attributes);
				clone.library->CompileBaseGraph(clone.emitter_attributes);
				clone.library->CompileVariationGraph(clone.emitter_attributes);
				clone.library->CompileOvertimeGraph(clone.emitter_attributes);
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

	void SetEffectPosition(tfxParticleManager &pm, tfxU32 effect_index, float x, float y) {
		tfxVec2 position(x, y);
		pm.effects.local_position[effect_index] = position;
	}

	void SetEffectPosition(tfxParticleManager &pm, tfxU32 effect_index, float x, float y, float z) {
		tfxVec3 position(x, y, z);
		pm.effects.local_position[effect_index] = position;
	}

	void tfxEffectEmitter::EnableAllEmitters() {
		for (auto &e : GetInfo().sub_effectors) {
			e.property_flags |= tfxEmitterPropertyFlags_enabled;
			e.EnableAllEmitters();
		}
	}

	void tfxEffectEmitter::EnableEmitter() {
		property_flags |= tfxEmitterPropertyFlags_enabled;
	}

	void tfxEffectEmitter::DisableAllEmitters() {
		for (auto &e : GetInfo().sub_effectors) {
			e.property_flags &= ~tfxEmitterPropertyFlags_enabled;
			e.DisableAllEmitters();
		}
	}

	void tfxEffectEmitter::DisableAllEmittersExcept(tfxEffectEmitter &emitter) {
		for (auto &e : GetInfo().sub_effectors) {
			if (e.library_index == emitter.library_index)
				e.property_flags |= tfxEmitterPropertyFlags_enabled;
			else
				e.property_flags &= ~tfxEmitterPropertyFlags_enabled;
		}
	}

	tfxGraph* tfxEffectEmitter::GetGraphByType(tfxGraphType type) {

		if (type < tfxGlobalCount) {
			return &((tfxGraph*)&library->global_graphs[global])[type];
		}
		else if (type >= tfxPropertyStart && type < tfxBaseStart) {
			int ref = type - tfxPropertyStart;
			return &((tfxGraph*)&library->emitter_attributes[emitter_attributes].properties)[ref];
		}
		else if (type >= tfxBaseStart && type < tfxVariationStart) {
			int ref = type - tfxBaseStart;
			return &((tfxGraph*)&library->emitter_attributes[emitter_attributes].base)[ref];
		}
		else if (type >= tfxVariationStart && type < tfxOvertimeStart) {
			int ref = type - tfxVariationStart;
			return &((tfxGraph*)&library->emitter_attributes[emitter_attributes].variation)[ref];
		}
		else if (type >= tfxOvertimeStart && type < tfxTransformStart) {
			int ref = type - tfxOvertimeStart;
			return &((tfxGraph*)&library->emitter_attributes[emitter_attributes].overtime)[ref];
		}
		else if (type >= tfxTransformStart) {
			int ref = type - tfxTransformStart;
			return &((tfxGraph*)&library->transform_attributes[transform_attributes].roll)[ref];
		}

		return nullptr;

	}

	tfxU32 tfxEffectEmitter::GetGraphIndexByType(tfxGraphType type) {

		if (type < tfxGlobalCount) {
			return global;
		}
		else if (type < tfxTransformStart) {
			return emitter_attributes;
		}
		else {
			return transform_attributes;
		}

	}

	void tfxEffectEmitter::FreeGraphs() {

		if (type == tfxEffectType) {
			library->global_graphs[global].life.Free();
			library->global_graphs[global].amount.Free();
			library->global_graphs[global].velocity.Free();
			library->global_graphs[global].width.Free();
			library->global_graphs[global].height.Free();
			library->global_graphs[global].weight.Free();
			library->global_graphs[global].spin.Free();
			library->global_graphs[global].stretch.Free();
			library->global_graphs[global].overal_scale.Free();
			library->global_graphs[global].intensity.Free();
			library->global_graphs[global].frame_rate.Free();
			library->global_graphs[global].splatter.Free();
			library->global_graphs[global].emitter_width.Free();
			library->global_graphs[global].emitter_height.Free();
			library->global_graphs[global].emitter_depth.Free();

			library->transform_attributes[transform_attributes].roll.Free();
			library->transform_attributes[transform_attributes].pitch.Free();
			library->transform_attributes[transform_attributes].yaw.Free();
			library->transform_attributes[transform_attributes].translation_x.Free();
			library->transform_attributes[transform_attributes].translation_y.Free();
			library->transform_attributes[transform_attributes].translation_z.Free();
		}

		if (type == tfxEmitterType) {
			library->transform_attributes[transform_attributes].translation_x.Free();
			library->transform_attributes[transform_attributes].translation_y.Free();
			library->transform_attributes[transform_attributes].translation_z.Free();

			library->emitter_attributes[emitter_attributes].properties.emission_pitch.Free();
			library->emitter_attributes[emitter_attributes].properties.emission_yaw.Free();
			library->emitter_attributes[emitter_attributes].properties.emission_range.Free();
			library->emitter_attributes[emitter_attributes].properties.splatter.Free();
			library->emitter_attributes[emitter_attributes].properties.emitter_width.Free();
			library->emitter_attributes[emitter_attributes].properties.emitter_height.Free();
			library->emitter_attributes[emitter_attributes].properties.emitter_depth.Free();
			library->emitter_attributes[emitter_attributes].properties.arc_size.Free();
			library->emitter_attributes[emitter_attributes].properties.arc_offset.Free();

			library->emitter_attributes[emitter_attributes].base.life.Free();
			library->emitter_attributes[emitter_attributes].base.amount.Free();
			library->emitter_attributes[emitter_attributes].base.velocity.Free();
			library->emitter_attributes[emitter_attributes].base.width.Free();
			library->emitter_attributes[emitter_attributes].base.height.Free();
			library->emitter_attributes[emitter_attributes].base.weight.Free();
			library->emitter_attributes[emitter_attributes].base.spin.Free();
			library->emitter_attributes[emitter_attributes].base.noise_offset.Free();

			library->emitter_attributes[emitter_attributes].variation.life.Free();
			library->emitter_attributes[emitter_attributes].variation.amount.Free();
			library->emitter_attributes[emitter_attributes].variation.velocity.Free();
			library->emitter_attributes[emitter_attributes].variation.width.Free();
			library->emitter_attributes[emitter_attributes].variation.height.Free();
			library->emitter_attributes[emitter_attributes].variation.weight.Free();
			library->emitter_attributes[emitter_attributes].variation.spin.Free();
			library->emitter_attributes[emitter_attributes].variation.noise_offset.Free();
			library->emitter_attributes[emitter_attributes].variation.noise_resolution.Free();

			library->emitter_attributes[emitter_attributes].overtime.velocity.Free();
			library->emitter_attributes[emitter_attributes].overtime.width.Free();
			library->emitter_attributes[emitter_attributes].overtime.height.Free();
			library->emitter_attributes[emitter_attributes].overtime.weight.Free();
			library->emitter_attributes[emitter_attributes].overtime.spin.Free();
			library->emitter_attributes[emitter_attributes].overtime.stretch.Free();
			library->emitter_attributes[emitter_attributes].overtime.red.Free();
			library->emitter_attributes[emitter_attributes].overtime.green.Free();
			library->emitter_attributes[emitter_attributes].overtime.blue.Free();
			library->emitter_attributes[emitter_attributes].overtime.blendfactor.Free();
			library->emitter_attributes[emitter_attributes].overtime.intensity.Free();
			library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance.Free();
			library->emitter_attributes[emitter_attributes].overtime.direction_turbulance.Free();
			library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.Free();
			library->emitter_attributes[emitter_attributes].overtime.direction.Free();
			library->emitter_attributes[emitter_attributes].overtime.noise_resolution.Free();
		}
	}

	tfxU32 tfxEffectEmitter::CountAllLookupValues() {
		tfxU32 count = 0;
		if (type == tfxEffectType) {
			count += library->CountGlobalLookUpValues(global);
			for (auto &emitter : GetInfo().sub_effectors) {
				count += library->CountEmitterLookUpValues(emitter.emitter_attributes);
			}
		}
		else if (type = tfxEmitterType) {
			count += library->CountEmitterLookUpValues(emitter_attributes);
		}
		return count;
	}

	void tfxEffectEmitter::CompileGraphs() {
		for (tfxU32 t = (tfxU32)tfxTransform_translate_x; t != (tfxU32)tfxGraphMaxIndex; ++t) {
			CompileGraph(*GetGraphByType(tfxGraphType(t)));
		}
		if (type == tfxEffectType) {
			for (tfxU32 t = (tfxU32)tfxGlobal_life; t != (tfxU32)tfxProperty_emission_pitch; ++t) {
				CompileGraph(*GetGraphByType(tfxGraphType(t)));
			}
		}
		else if (type == tfxEmitterType) {
			for (tfxU32 t = (tfxU32)tfxProperty_emission_pitch; t != (tfxU32)tfxOvertime_velocity; ++t) {
				CompileGraph(*GetGraphByType((tfxGraphType)t));
			}
			for (tfxU32 t = (tfxU32)tfxOvertime_velocity; t != (tfxU32)tfxTransform_translate_x; ++t) {
				CompileGraphOvertime(*GetGraphByType((tfxGraphType)t));
			}
		}
		for (auto &sub : GetInfo().sub_effectors) {
			sub.CompileGraphs();
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
		effect.library = this;
		effects.push_back(effect);
		ReIndex();
		UpdateEffectPaths();
		return effects.back();
	}

	tfxEffectEmitter &tfxEffectLibrary::AddFolder(tfxStr64 &name) {
		tfxEffectEmitter folder;
		folder.info_index = AddEffectEmitterInfo();
		folder.library = this;
		folder.GetInfo().name = name;
		folder.type = tfxFolder;
		folder.library = this;
		folder.GetInfo().uid = ++uid;
		effects.push_back(folder);
		ReIndex();
		UpdateEffectPaths();
		return effects.back();
	}

	tfxEffectEmitter &tfxEffectLibrary::AddFolder(tfxEffectEmitter &folder) {
		assert(folder.type == tfxFolder);			//Must be type tfxFolder if adding a folder
		folder.library = this;
		folder.GetInfo().uid = ++uid;
		effects.push_back(folder);
		ReIndex();
		UpdateEffectPaths();
		return effects.back();
	}

	tfxEffectEmitter &tfxEffectLibrary::AddStage(tfxStr64 &name) {
		tfxEffectEmitter stage;
		stage.info_index = AddEffectEmitterInfo();
		stage.library = this;
		stage.GetInfo().name = name;
		stage.type = tfxStage;
		stage.library = this;
		stage.GetInfo().uid = ++uid;
		effects.push_back(stage);
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

	void tfxEffectLibrary::PrepareEffectTemplate(tfxEffectEmitter &effect, tfxEffectTemplate &effect_template) {
		assert(effect.type == tfxEffectType);
		effect.Clone(effect_template.effect, &effect_template.effect, this);
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
			if (effect.type == tfxFolder || effect.type == tfxStage) {
				UpdateParticleShapeReferences(effect.GetInfo().sub_effectors, default_index);
			}
			else {
				for (auto &emitter : effect.GetInfo().sub_effectors) {
					if (particle_shapes.ValidIntName(emitter_properties.shape_index[emitter.property_index])) {
						emitter_properties.image[emitter.property_index] = &particle_shapes.AtInt(emitter_properties.shape_index[emitter.property_index]);
						emitter_properties.end_frame[emitter.property_index] = particle_shapes.AtInt(emitter_properties.shape_index[emitter.property_index]).animation_frames - 1;
					}
					else {
						emitter_properties.image[emitter.property_index] = &particle_shapes.AtInt(default_index);
						emitter_properties.end_frame[emitter.property_index] = particle_shapes.AtInt(default_index).animation_frames - 1;
					}
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

	tfxU32 tfxEffectLibrary::AddKeyframes() {
		if (free_keyframes.size())
			return free_keyframes.pop_back();
		tfxTransformAttributes keyframes;
		keyframes.Initialise(&graph_node_allocator, &graph_lookup_allocator);
		transform_attributes.push_back(keyframes);
		return transform_attributes.size() - 1;
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

	void tfxEffectLibrary::FreeKeyframes(tfxU32 index) {
		assert(index < transform_attributes.size());
		free_keyframes.push_back(index);
		transform_attributes[index].Free();
	}

	void tfxEffectLibrary::FreeEmitterAttributes(tfxU32 index) {
		assert(index < emitter_attributes.size());
		free_emitter_attributes.push_back(index);
		emitter_attributes[index].Free();
	}

	void tfxEffectLibrary::FreeProperties(tfxU32 index) {
		assert(index < emitter_properties_buffer.current_size);
		free_properties.push_back(index);
	}

	void tfxEffectLibrary::FreeInfo(tfxU32 index) {
		assert(index < effect_infos.size());
		free_infos.push_back(index);
	}

	tfxU32 tfxEffectLibrary::CountKeyframeLookUpValues(tfxU32 index) {
		auto &transform = transform_attributes[index];
		tfxU32 count = 0;
		count += transform.roll.lookup.values.capacity;
		count += transform.pitch.lookup.values.capacity;
		count += transform.yaw.lookup.values.capacity;
		count += transform.translation_x.lookup.values.capacity;
		count += transform.translation_y.lookup.values.capacity;
		count += transform.translation_z.lookup.values.capacity;
		return count;
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

	tfxU32 tfxEffectLibrary::CloneKeyframes(tfxU32 source_index, tfxEffectLibrary *destination_library) {
		tfxU32 index = destination_library->AddKeyframes();
		transform_attributes[source_index].CopyToNoLookups(&destination_library->transform_attributes[index]);
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
		CopyEmitterProperites(emitter_properties, source_index, destination_library->emitter_properties, index);
		return index;
	}

	void tfxEffectLibrary::AddEmitterGraphs(tfxEffectEmitter& emitter) {
		emitter.emitter_attributes = AddEmitterAttributes();
	}

	void tfxEffectLibrary::AddTransformGraphs(tfxEffectEmitter& emitter) {
		emitter.transform_attributes = AddKeyframes();
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
		assert(effect.type == tfxEffectType || effect.type == tfxStage);
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
		if (free_properties.size()) {
			return free_properties.pop_back();
		}
		return AddRow(&emitter_properties_buffer, true);
	}

	void tfxEffectLibrary::Init() {
		graph_node_allocator = CreateArenaManager(tfxMegabyte(2), 8);
		graph_lookup_allocator = CreateArenaManager(tfxMegabyte(4), 8); 
		InitEmitterProperties();
	}

	void tfxEffectLibrary::InitEmitterProperties() {
		AddStructArray(&emitter_properties_buffer, sizeof(tfxVec3), offsetof(tfxEmitterPropertiesSoA, angle_offsets));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxVectorAlignType), offsetof(tfxEmitterPropertiesSoA, vector_align_type));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxEmissionType), offsetof(tfxEmitterPropertiesSoA, emission_type));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxU32), offsetof(tfxEmitterPropertiesSoA, single_shot_limit));
		AddStructArray(&emitter_properties_buffer, sizeof(float), offsetof(tfxEmitterPropertiesSoA, frame_rate));
		AddStructArray(&emitter_properties_buffer, sizeof(float), offsetof(tfxEmitterPropertiesSoA, end_frame));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxImageData*), offsetof(tfxEmitterPropertiesSoA, image));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxBillboardingOptions), offsetof(tfxEmitterPropertiesSoA, billboard_option));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxVec3), offsetof(tfxEmitterPropertiesSoA, grid_points));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxAngleSettingFlags), offsetof(tfxEmitterPropertiesSoA, angle_settings));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxU32), offsetof(tfxEmitterPropertiesSoA, layer));
		AddStructArray(&emitter_properties_buffer, sizeof(float), offsetof(tfxEmitterPropertiesSoA, delay_spawning));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxEmissionDirection), offsetof(tfxEmitterPropertiesSoA, emission_direction));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxLineTraversalEndBehaviour), offsetof(tfxEmitterPropertiesSoA, end_behaviour));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxParticleControlFlags), offsetof(tfxEmitterPropertiesSoA, compute_flags));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxVec2), offsetof(tfxEmitterPropertiesSoA, image_handle));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxVec3), offsetof(tfxEmitterPropertiesSoA, emitter_handle));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxU32), offsetof(tfxEmitterPropertiesSoA, spawn_amount));
		AddStructArray(&emitter_properties_buffer, sizeof(tfxU32), offsetof(tfxEmitterPropertiesSoA, shape_index));
		AddStructArray(&emitter_properties_buffer, sizeof(float), offsetof(tfxEmitterPropertiesSoA, loop_length));
		AddStructArray(&emitter_properties_buffer, sizeof(float), offsetof(tfxEmitterPropertiesSoA, start_frame));
		FinishSoABufferSetup(&emitter_properties_buffer, &emitter_properties, 100);
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
		transform_attributes.free_all();
		animation_settings.free_all();
		preview_camera_settings.free_all();
		FreeSoABuffer(&emitter_properties_buffer);
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
						graph.nodes.ResetIteratorIndex();
						do {
							for (auto &node : graph.nodes) {
								all_nodes.push_back(node);
								running_node_index++;
							}
						} while (!graph.nodes.EndOfBuckets());

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
			CompileKeyframeGraph(effect.transform_attributes);
			CompileGlobalGraph(effect.global);
		}
		else if (effect.type == tfxEffectType) {
			CompileKeyframeGraph(effect.transform_attributes);
		}
		else if(effect.type == tfxEmitterType) {
			CompileKeyframeGraph(effect.transform_attributes);
			CompileEmitterGraphs(effect.emitter_attributes);
			for (auto &sub : info.sub_effectors) {
				CompileGraphsOfEffect(sub, ++depth);
			}
		}
		else if (effect.type == tfxFolder) {
			for (auto &sub : info.sub_effectors) {
				CompileGraphsOfEffect(sub, 0);
			}
		}
		if (effect.type == tfxEffectType) {
			for (auto &sub : info.sub_effectors) {
				CompileGraphsOfEffect(sub, ++depth);
			}
		}
	}

	void tfxEffectLibrary::CompileAllGraphs() {
		for (auto &g : global_graphs) {
			CompileGraph(g.amount);
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
		for (auto &g : transform_attributes) {
			CompileGraph(g.roll);
			CompileGraph(g.pitch);
			CompileGraph(g.yaw);
			CompileGraph(g.translation_x);
			CompileGraph(g.translation_y);
			CompileGraph(g.translation_z);
		}
		for (auto &g : emitter_attributes) {
			CompileGraph(g.properties.arc_offset);
			CompileGraph(g.properties.arc_size);
			CompileGraph(g.properties.emission_pitch);
			CompileGraph(g.properties.emission_yaw);
			CompileGraph(g.properties.emission_range);
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

	void tfxEffectLibrary::CompileKeyframeGraph(tfxU32 index) {
		tfxTransformAttributes &g = transform_attributes[index];
		CompileGraph(g.roll);
		CompileGraph(g.pitch);
		CompileGraph(g.yaw);
		CompileGraph(g.translation_x);
		CompileGraph(g.translation_y);
		CompileGraph(g.translation_z);
	}

	void tfxEffectLibrary::CompileEmitterGraphs(tfxU32 index) {
		CompilePropertyGraph(index);
		CompileKeyframeGraph(index);
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
		CompileGraph(g.emitter_width);
		CompileGraph(g.emitter_height);
		CompileGraph(g.emitter_depth);
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
		graph_min_max[tfxGlobal_emitter_width] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_emitter_height] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_emitter_depth] = GetMinMaxGraphValues(tfxGlobalPercentPreset);

		graph_min_max[tfxTransform_roll] = GetMinMaxGraphValues(tfxAnglePreset);
		graph_min_max[tfxTransform_pitch] = GetMinMaxGraphValues(tfxAnglePreset);
		graph_min_max[tfxTransform_yaw] = GetMinMaxGraphValues(tfxAnglePreset);
		graph_min_max[tfxTransform_translate_x] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxTransform_translate_y] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxTransform_translate_z] = GetMinMaxGraphValues(tfxDimensionsPreset);

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
		names_and_types.Insert("grid_spawn_random", tfxBool);
		names_and_types.Insert("area_open_ends", tfxBool);
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
			if (values[0] == "global_amount") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].amount.AddNode(n); }
			if (values[0] == "global_frame_rate") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].frame_rate.AddNode(n); }
			if (values[0] == "global_height") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].height.AddNode(n); }
			if (values[0] == "global_width") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].width.AddNode(n); }
			if (values[0] == "global_life") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].life.AddNode(n); }
			if (values[0] == "global_opacity") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].intensity.AddNode(n); }
			if (values[0] == "global_spin") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].spin.AddNode(n); }
			if (values[0] == "global_splatter") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].splatter.AddNode(n); }
			if (values[0] == "global_stretch") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].stretch.AddNode(n); }
			if (values[0] == "global_overal_scale") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].overal_scale.AddNode(n); }
			if (values[0] == "global_weight") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].weight.AddNode(n); }
			if (values[0] == "global_velocity") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].velocity.AddNode(n); }
			if (values[0] == "global_emitter_width") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].emitter_width.AddNode(n); }
			if (values[0] == "global_emitter_height") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].emitter_height.AddNode(n); }
			if (values[0] == "global_emitter_depth") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].emitter_depth.AddNode(n); }

			if (values[0] == "global_effect_angle") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->transform_attributes[effect.transform_attributes].roll.AddNode(n); }
			if (values[0] == "global_effect_roll") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->transform_attributes[effect.transform_attributes].roll.AddNode(n); }
			if (values[0] == "global_effect_pitch") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->transform_attributes[effect.transform_attributes].pitch.AddNode(n); }
			if (values[0] == "global_effect_yaw") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->transform_attributes[effect.transform_attributes].yaw.AddNode(n); }
			if (values[0] == "keyframe_translate_x") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->transform_attributes[effect.transform_attributes].translation_x.AddNode(n); }
			if (values[0] == "keyframe_translate_y") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->transform_attributes[effect.transform_attributes].translation_y.AddNode(n); }
			if (values[0] == "keyframe_translate_z") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->transform_attributes[effect.transform_attributes].translation_z.AddNode(n); }
			if (values[0] == "property_emitter_angle") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->transform_attributes[effect.transform_attributes].roll.AddNode(n); }
			if (values[0] == "property_emitter_roll") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->transform_attributes[effect.transform_attributes].roll.AddNode(n); }
			if (values[0] == "property_emitter_pitch") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->transform_attributes[effect.transform_attributes].pitch.AddNode(n); }
			if (values[0] == "property_emitter_yaw") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->transform_attributes[effect.transform_attributes].yaw.AddNode(n); }

			if (values[0] == "base_arc_offset") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.arc_offset.AddNode(n); }
			if (values[0] == "base_arc_size") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.arc_size.AddNode(n); }
			if (values[0] == "base_emission_angle") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.emission_pitch.AddNode(n); }
			if (values[0] == "base_emission_range") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.emission_range.AddNode(n); }
			if (values[0] == "base_emitter_height") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.emitter_height.AddNode(n); }
			if (values[0] == "base_emitter_width") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.emitter_width.AddNode(n); }
			if (values[0] == "base_splatter") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.splatter.AddNode(n); }

			if (values[0] == "property_arc_offset") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.arc_offset.AddNode(n); }
			if (values[0] == "property_arc_size") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.arc_size.AddNode(n); }
			if (values[0] == "property_emission_angle") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.emission_pitch.AddNode(n); }
			if (values[0] == "property_emission_pitch") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.emission_pitch.AddNode(n); }
			if (values[0] == "property_emission_yaw") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.emission_yaw.AddNode(n); }
			if (values[0] == "property_emission_range") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.emission_range.AddNode(n); }
			if (values[0] == "property_emitter_height") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.emitter_height.AddNode(n); }
			if (values[0] == "property_emitter_width") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.emitter_width.AddNode(n); }
			if (values[0] == "property_emitter_depth") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.emitter_depth.AddNode(n); }
			if (values[0] == "property_splatter") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].properties.splatter.AddNode(n); }

			if (values[0] == "base_amount") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].base.amount.AddNode(n); }
			if (values[0] == "base_life") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].base.life.AddNode(n); }
			if (values[0] == "base_height") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].base.height.AddNode(n); }
			if (values[0] == "base_width") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].base.width.AddNode(n); }
			if (values[0] == "base_spin") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].base.spin.AddNode(n); }
			if (values[0] == "base_noise_offset") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].base.noise_offset.AddNode(n); }
			if (values[0] == "base_velocity") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].base.velocity.AddNode(n); }
			if (values[0] == "base_weight") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].base.weight.AddNode(n); }

			if (values[0] == "variation_amount") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].variation.amount.AddNode(n); }
			if (values[0] == "variation_height") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].variation.height.AddNode(n); }
			if (values[0] == "variation_width") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].variation.width.AddNode(n); }
			if (values[0] == "variation_life") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].variation.life.AddNode(n); }
			if (values[0] == "variation_velocity") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].variation.velocity.AddNode(n); }
			if (values[0] == "variation_weight") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].variation.weight.AddNode(n); }
			if (values[0] == "variation_spin") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].variation.spin.AddNode(n); }
			if (values[0] == "variation_motion_randomness") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].variation.noise_offset.AddNode(n); }
			if (values[0] == "variation_noise_offset") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].variation.noise_offset.AddNode(n); }
			if (values[0] == "variation_noise_resolution") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].variation.noise_resolution.AddNode(n); }

			if (values[0] == "overtime_red") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.red.AddNode(n); }
			if (values[0] == "overtime_green") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.green.AddNode(n); }
			if (values[0] == "overtime_blue") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.blue.AddNode(n); }
			if (values[0] == "overtime_opacity") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.blendfactor.AddNode(n); }
			if (values[0] == "overtime_intensity") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.intensity.AddNode(n); }
			if (values[0] == "overtime_velocity_turbulance") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.velocity_turbulance.AddNode(n); }
			if (values[0] == "overtime_spin") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.spin.AddNode(n); }
			if (values[0] == "overtime_stretch") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.stretch.AddNode(n); }
			if (values[0] == "overtime_velocity") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.velocity.AddNode(n); }
			if (values[0] == "overtime_weight") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.weight.AddNode(n); }
			if (values[0] == "overtime_width") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.width.AddNode(n); }
			if (values[0] == "overtime_height") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.height.AddNode(n); }
			if (values[0] == "overtime_direction_turbulance") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.direction_turbulance.AddNode(n); }
			if (values[0] == "overtime_direction") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.direction.AddNode(n); }
			if (values[0] == "overtime_velocity_adjuster") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.velocity_adjuster.AddNode(n); }
			if (values[0] == "overtime_noise_resolution") { tfxAttributeNode n; AssignNodeData(n, values); effect.library->emitter_attributes[effect.emitter_attributes].overtime.noise_resolution.AddNode(n); }
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

	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, uint32_t value) {
	}

	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, float value) {
	}

	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, bool value) {
	}

	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, int value) {
	}

	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, tfxStr &value) {
		if (field == "name") {
			effect.GetInfo().name = value;
		}
	}

	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, uint32_t value) {
		tfxEmitterPropertiesSoA &emitter_properties = effect.library->emitter_properties;
		if (field == "image_index")
			emitter_properties.shape_index[effect.property_index] = value;
		if (field == "spawn_amount")
			emitter_properties.spawn_amount[effect.property_index] = value;
		if (field == "frames")
			effect.library->animation_settings[effect.GetInfo().animation_settings].frames = value;
		if (field == "current_frame")
			effect.library->animation_settings[effect.GetInfo().animation_settings].current_frame = value;
		if (field == "seed")
			effect.library->animation_settings[effect.GetInfo().animation_settings].seed = value;
		if (field == "layer")
			emitter_properties.layer[effect.property_index] = value;
		if (field == "frame_offset")
			effect.library->animation_settings[effect.GetInfo().animation_settings].frame_offset = value;
		if (field == "single_shot_limit")
			emitter_properties.single_shot_limit[effect.property_index] = value;
		if (field == "billboard_option")
			emitter_properties.billboard_option[effect.property_index] = (tfxBillboardingOptions)value;
		if (field == "vector_align_type")
			emitter_properties.vector_align_type[effect.property_index] = (tfxVectorAlignType)value;
		if (field == "angle_setting")
			emitter_properties.angle_settings[effect.property_index] = (tfxAngleSettingFlags)value;
		if (field == "sort_passes")
			effect.sort_passes = value;
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, int value) {
		tfxEmitterPropertiesSoA &emitter_properties = effect.library->emitter_properties;
		if (field == "emission_type")
			emitter_properties.emission_type[effect.property_index] = (tfxEmissionType)value;
		if (field == "emission_direction")
			emitter_properties.emission_direction[effect.property_index] = (tfxEmissionDirection)value;
		if (field == "color_option")
			effect.library->animation_settings[effect.GetInfo().animation_settings].color_option = (tfxExportColorOptions)value;
		if (field == "export_option")
			effect.library->animation_settings[effect.GetInfo().animation_settings].export_option = (tfxExportOptions)value;
		if (field == "end_behaviour")
			emitter_properties.end_behaviour[effect.property_index] = (tfxLineTraversalEndBehaviour)value;
		if (field == "frame_offset")
			effect.library->animation_settings[effect.GetInfo().animation_settings].frame_offset = value;
		if (field == "extra_frames_count")
			effect.library->animation_settings[effect.GetInfo().animation_settings].extra_frames_count = value;
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, tfxStr &value) {
		if (field == "name") {
			effect.GetInfo().name = value;
		}
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, float value) {
		tfxEmitterPropertiesSoA &emitter_properties = effect.library->emitter_properties;
		if (field == "position_x")
			effect.library->animation_settings[effect.GetInfo().animation_settings].position.x = value;
		if (field == "position_y")
			effect.library->animation_settings[effect.GetInfo().animation_settings].position.y = value;
		if (field == "frame_width")
			effect.library->animation_settings[effect.GetInfo().animation_settings].frame_size.x = value;
		if (field == "frame_height")
			effect.library->animation_settings[effect.GetInfo().animation_settings].frame_size.y = value;
		if (field == "zoom")
			effect.library->animation_settings[effect.GetInfo().animation_settings].zoom = value;
		if (field == "scale")
			effect.library->animation_settings[effect.GetInfo().animation_settings].scale = value;
		if (field == "playback_speed")
			effect.library->animation_settings[effect.GetInfo().animation_settings].playback_speed = value;
		if (field == "camera_position_x")
			effect.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_position.x = value;
		if (field == "camera_position_y")
			effect.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_position.y = value;
		if (field == "camera_position_z")
			effect.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_position.z = value;
		if (field == "camera_pitch")
			effect.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_pitch = value;
		if (field == "camera_yaw")
			effect.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_yaw = value;
		if (field == "camera_fov")
			effect.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_fov = value;
		if (field == "camera_floor_height")
			effect.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_floor_height = value;
		if (field == "camera_isometric_scale")
			effect.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_isometric_scale = value;
		if (field == "preview_camera_position_x")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_position.x = value;
		if (field == "preview_camera_position_y")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_position.y = value;
		if (field == "preview_camera_position_z")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_position.z = value;
		if (field == "preview_camera_pitch")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_pitch = value;
		if (field == "preview_camera_yaw")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_yaw = value;
		if (field == "preview_camera_fov")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_fov = value;
		if (field == "preview_camera_floor_height")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_floor_height = value;
		if (field == "preview_camera_isometric_scale")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_isometric_scale = value == 0 ? 5.f : value;
		if (field == "preview_effect_z_offset")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].effect_z_offset = value;
		if (field == "preview_camera_speed")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_speed = value;
		if (field == "image_handle_x")
			emitter_properties.image_handle[effect.property_index].x = value;
		if (field == "image_handle_y")
			emitter_properties.image_handle[effect.property_index].y = value;
		if (field == "delay_spawning")
			emitter_properties.delay_spawning[effect.property_index] = value;
		if (field == "grid_rows")
			emitter_properties.grid_points[effect.property_index].x = value;
		if (field == "grid_columns")
			emitter_properties.grid_points[effect.property_index].y = value;
		if (field == "grid_depth")
			emitter_properties.grid_points[effect.property_index].z = value;
		if (field == "loop_length")
			emitter_properties.loop_length[effect.property_index] = value < 0 ? 0.f : value;
		if (field == "emitter_handle_x")
			emitter_properties.emitter_handle[effect.property_index].x = value;
		if (field == "emitter_handle_y")
			emitter_properties.emitter_handle[effect.property_index].y = value;
		if (field == "emitter_handle_z")
			emitter_properties.emitter_handle[effect.property_index].z = value;
		if (field == "image_start_frame")
			emitter_properties.start_frame[effect.property_index] = value;
		if (field == "image_end_frame")
			emitter_properties.end_frame[effect.property_index] = value;
		if (field == "image_frame_rate")
			emitter_properties.frame_rate[effect.property_index] = value;
		if (field == "angle_offset")
			emitter_properties.angle_offsets[effect.property_index].roll = value;
		if (field == "angle_offset_pitch")
			emitter_properties.angle_offsets[effect.property_index].pitch = value;
		if (field == "angle_offset_yaw")
			emitter_properties.angle_offsets[effect.property_index].yaw = value;
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, bool value) {
		if (field == "loop")
			effect.library->animation_settings[effect.GetInfo().animation_settings].loop = value;
		if (field == "seamless")
			effect.library->animation_settings[effect.GetInfo().animation_settings].seamless = value;
		if (field == "export_with_transparency")
			effect.library->animation_settings[effect.GetInfo().animation_settings].export_with_transparency = value;
		if (field == "camera_isometric")
			effect.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_isometric = value;
		if (field == "camera_hide_floor")
			effect.library->animation_settings[effect.GetInfo().animation_settings].camera_settings.camera_hide_floor = value;
		if (field == "preview_attach_effect_to_camera")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].attach_effect_to_camera = value;
		if (field == "preview_camera_hide_floor")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_hide_floor = value;
		if (field == "preview_camera_isometric")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_isometric = value;
		if (field == "random_color")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_random_color; else effect.property_flags &= ~tfxEmitterPropertyFlags_random_color;
		if (field == "relative_position")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_relative_position; else effect.property_flags &= ~tfxEmitterPropertyFlags_relative_position;
		if (field == "relative_angle")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_relative_angle; else effect.property_flags &= ~tfxEmitterPropertyFlags_relative_angle;
		if (field == "image_handle_auto_center")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_image_handle_auto_center; else effect.property_flags &= ~tfxEmitterPropertyFlags_image_handle_auto_center;
		if (field == "single")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_single; else effect.property_flags &= ~tfxEmitterPropertyFlags_single;
		//if (field == "one_shot")
			//if(value) effect.property_flags |= tfxEmitterPropertyFlags_one_shot; else effect.property_flags &= ~tfxEmitterPropertyFlags_one_shot;
		if (field == "spawn_on_grid")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_spawn_on_grid; else effect.property_flags &= ~tfxEmitterPropertyFlags_spawn_on_grid;
		if (field == "grid_spawn_clockwise")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_grid_spawn_clockwise; else effect.property_flags &= ~tfxEmitterPropertyFlags_grid_spawn_clockwise;
		if (field == "fill_area")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_fill_area; else effect.property_flags &= ~tfxEmitterPropertyFlags_fill_area;
		if (field == "grid_spawn_random")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_grid_spawn_random; else effect.property_flags &= ~tfxEmitterPropertyFlags_grid_spawn_random;
		if (field == "area_open_ends")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_area_open_ends; else effect.property_flags &= ~tfxEmitterPropertyFlags_area_open_ends;
		if (field == "emitter_handle_auto_center")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_emitter_handle_auto_center; else effect.property_flags &= ~tfxEmitterPropertyFlags_emitter_handle_auto_center;
		if (field == "edge_traversal")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_edge_traversal; else effect.property_flags &= ~tfxEmitterPropertyFlags_edge_traversal;
		if (field == "image_reverse_animation")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_reverse_animation; else effect.property_flags &= ~tfxEmitterPropertyFlags_reverse_animation;
		if (field == "image_play_once")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_play_once; else effect.property_flags &= ~tfxEmitterPropertyFlags_play_once;
		if (field == "image_animate")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_animate; else effect.property_flags &= ~tfxEmitterPropertyFlags_animate;
		if (field == "image_random_start_frame")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_random_start_frame; else effect.property_flags &= ~tfxEmitterPropertyFlags_random_start_frame;
		if (field == "global_uniform_size")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_global_uniform_size; else effect.property_flags &= ~tfxEmitterPropertyFlags_global_uniform_size;
		if (field == "base_uniform_size")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_base_uniform_size; else effect.property_flags &= ~tfxEmitterPropertyFlags_base_uniform_size;
		if (field == "lifetime_uniform_size")
			if(value) effect.property_flags |= tfxEmitterPropertyFlags_lifetime_uniform_size; else effect.property_flags &= ~tfxEmitterPropertyFlags_lifetime_uniform_size;
		if (field == "use_spawn_ratio")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_use_spawn_ratio; else effect.property_flags &= ~tfxEmitterPropertyFlags_use_spawn_ratio;
		if (field == "is_3d")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_is_3d; else effect.property_flags &= ~tfxEmitterPropertyFlags_is_3d;
		if (field == "draw_order_by_age")
			if (value) effect.effect_flags |= tfxEffectPropertyFlags_age_order; else effect.effect_flags &= ~tfxEffectPropertyFlags_age_order;
		if (field == "draw_order_by_depth")
			if (value) effect.effect_flags |= tfxEffectPropertyFlags_depth_draw_order; else effect.effect_flags &= ~tfxEffectPropertyFlags_depth_draw_order;
		if (field == "guaranteed_draw_order")
			if (value) effect.effect_flags |= tfxEffectPropertyFlags_guaranteed_order; else effect.effect_flags &= ~tfxEffectPropertyFlags_guaranteed_order;
	}

	void StreamProperties(tfxEmitterPropertiesSoA &property, tfxU32 index, tfxEmitterPropertyFlags &flags, tfxStr &file) {

		file.AddLine("image_index=%i", property.shape_index[index]);
		file.AddLine("image_handle_x=%f", property.image_handle[index].x);
		file.AddLine("image_handle_y=%f", property.image_handle[index].y);
		file.AddLine("image_start_frame=%f", property.start_frame[index]);
		file.AddLine("image_end_frame=%f", property.end_frame[index]);
		file.AddLine("image_frame_rate=%f", property.frame_rate[index]);
		file.AddLine("image_play_once=%i", (flags & tfxEmitterPropertyFlags_play_once));
		file.AddLine("image_reverse_animation=%i", (flags & tfxEmitterPropertyFlags_reverse_animation));
		file.AddLine("image_animate=%i", (flags & tfxEmitterPropertyFlags_animate));
		file.AddLine("image_random_start_frame=%i", (flags & tfxEmitterPropertyFlags_random_start_frame));
		file.AddLine("image_handle_auto_center=%i", (flags & tfxEmitterPropertyFlags_image_handle_auto_center));
		file.AddLine("spawn_amount=%i", property.spawn_amount[index]);
		file.AddLine("emission_type=%i", property.emission_type[index]);
		file.AddLine("emission_direction=%i", property.emission_direction[index]);
		file.AddLine("grid_rows=%f", property.grid_points[index].x);
		file.AddLine("grid_columns=%f", property.grid_points[index].y);
		file.AddLine("grid_depth=%f", property.grid_points[index].z);
		file.AddLine("delay_spawning=%f", property.delay_spawning[index]);
		file.AddLine("loop_length=%f", property.loop_length[index]);
		file.AddLine("emitter_handle_x=%f", property.emitter_handle[index].x);
		file.AddLine("emitter_handle_y=%f", property.emitter_handle[index].y);
		file.AddLine("emitter_handle_z=%f", property.emitter_handle[index].z);
		file.AddLine("end_behaviour=%i", property.end_behaviour[index]);
		file.AddLine("random_color=%i", (flags & tfxEmitterPropertyFlags_random_color));
		file.AddLine("relative_position=%i", (flags & tfxEmitterPropertyFlags_relative_position));
		file.AddLine("relative_angle=%i", (flags & tfxEmitterPropertyFlags_relative_angle));
		file.AddLine("single=%i", (flags & tfxEmitterPropertyFlags_single));
		file.AddLine("single_shot_limit=%i", property.single_shot_limit);
		file.AddLine("spawn_on_grid=%i", (flags & tfxEmitterPropertyFlags_spawn_on_grid));
		file.AddLine("grid_spawn_clockwise=%i", (flags & tfxEmitterPropertyFlags_grid_spawn_clockwise));
		file.AddLine("fill_area=%i", (flags & tfxEmitterPropertyFlags_fill_area));
		file.AddLine("grid_spawn_random=%i", (flags & tfxEmitterPropertyFlags_grid_spawn_random));
		file.AddLine("area_open_ends=%i", (flags & tfxEmitterPropertyFlags_area_open_ends));
		file.AddLine("emitter_handle_auto_center=%i", (flags & tfxEmitterPropertyFlags_emitter_handle_auto_center));
		file.AddLine("edge_traversal=%i", (flags & tfxEmitterPropertyFlags_edge_traversal));
		file.AddLine("angle_setting=%i", property.angle_settings[index]);
		file.AddLine("angle_offset=%f", property.angle_offsets[index].roll);
		file.AddLine("angle_offset_pitch=%f", property.angle_offsets[index].pitch);
		file.AddLine("angle_offset_yaw=%f", property.angle_offsets[index].yaw);
		file.AddLine("global_uniform_size=%i", (flags & tfxEmitterPropertyFlags_global_uniform_size));
		file.AddLine("base_uniform_size=%i", (flags & tfxEmitterPropertyFlags_base_uniform_size));
		file.AddLine("lifetime_uniform_size=%i", (flags & tfxEmitterPropertyFlags_lifetime_uniform_size));
		file.AddLine("use_spawn_ratio=%i", (flags & tfxEmitterPropertyFlags_use_spawn_ratio));
		file.AddLine("is_3d=%i", (flags & tfxEmitterPropertyFlags_is_3d));
		file.AddLine("billboard_option=%i", property.billboard_option[index]);
		file.AddLine("vector_align_type=%i", property.vector_align_type[index]);
		file.AddLine("layer=%i", property.layer[index]);

	}

	void StreamProperties(tfxEffectEmitter &effect, tfxStr &file) {
		file.AddLine("draw_order_by_age=%i", effect.effect_flags & tfxEffectPropertyFlags_age_order);
		file.AddLine("draw_order_by_depth=%i", effect.effect_flags & tfxEffectPropertyFlags_depth_draw_order);
		file.AddLine("guaranteed_draw_order=%i", effect.effect_flags & tfxEffectPropertyFlags_guaranteed_order);
		file.AddLine("sort_passes=%i", effect.sort_passes);
	}

	void StreamGraph(const char *name, tfxGraph &graph, tfxStr &file) {

		graph.nodes.ResetIteratorIndex();
		do {
			for (auto &n : graph.nodes) {
				file.AddLine("%s,%f,%f,%i,%f,%f,%f,%f", name, n.frame, n.value, (n.flags & tfxAttributeNodeFlags_is_curve), n.left.x, n.left.y, n.right.x, n.right.y);
			}
		} while (!graph.nodes.EndOfBuckets());

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
		return (type == tfxTransform_roll || type == tfxTransform_pitch || type == tfxTransform_yaw || type == tfxProperty_emission_pitch || type == tfxProperty_emission_yaw 
			|| type == tfxProperty_emission_range || type == tfxProperty_arc_offset || type == tfxProperty_arc_size || type == tfxBase_spin || type == tfxVariation_spin || type == tfxOvertime_direction);
	}

	bool tfxGraph::IsTranslationGraph() {
		return type == tfxTransform_translate_x || type == tfxTransform_translate_y || type == tfxTransform_translate_z;
	}

	void tfxGraph::MultiplyAllValues(float scalar) {
		nodes.ResetIteratorIndex();
		do {
			for (auto &node : nodes) {
				node.value *= scalar;
				node.left.y *= scalar;
				node.right.y *= scalar;
			}
		} while (!nodes.EndOfBuckets());
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
		graph.nodes.ResetIteratorIndex();
		do {
			for (auto &node : graph.nodes) {
				ClampNode(graph, node);
				if (node.flags & tfxAttributeNodeFlags_is_curve) {
					ClampCurve(graph, node.left, node);
					ClampCurve(graph, node.right, node);
				}
			}
		} while (!graph.nodes.EndOfBuckets());
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
		nodes.ResetIteratorIndex();
		do {
			for (auto &n : nodes) {
				if (n.frame == node.frame && n.value == node.value)
					return;
			}
		} while (!nodes.EndOfBuckets());
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
			nodes.ResetIteratorIndex();
			do {
				for (auto *n = nodes.begin() + 1; n != nodes.end(); ++n) {
					if (node.frame < n->frame)
						last_node = n;
					else
						break;
				}
			} while (!nodes.EndOfBuckets());

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
		nodes.ResetIteratorIndex();
		do {
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
		} while (!nodes.EndOfBuckets());
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
		nodes.ResetIteratorIndex();
		do {
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
		} while (!nodes.EndOfBuckets());
		return random_generation.Range(lastv);

	}

	float tfxGraph::GetValue(float age, float life) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		tfxAttributeNode *lastec = nullptr;
		nodes.ResetIteratorIndex();
		do {
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
		} while (!nodes.EndOfBuckets());
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
			nodes.ResetIteratorIndex();
			do {
				for (auto &n : nodes) {
					if (value < n.value)
						value = n.value;
				}
			} while (!nodes.EndOfBuckets());
			return value;
		}
		return 0.f;
	}

	float tfxGraph::GetMinValue() {
		if (nodes.size()) {
			float value = tfxMAX_FLOAT;
			nodes.ResetIteratorIndex();
			do {
				for (auto &n : nodes) {
					if (value > n.value)
						value = n.value;
				}
			} while (!nodes.EndOfBuckets());
			return value;
		}
		return 0.f;
	}

	float tfxGraph::GetLastFrame() {
		if (nodes.size()) {
			return nodes.size() > 1 && nodes.back().frame == 0 ? tfxFRAME_LENGTH : nodes.back().frame;
		}

		return 0.f;
	}

	tfxAttributeNode* tfxGraph::FindNode(const tfxAttributeNode &n) {
		return nodes.find(n);
	}

	void tfxGraph::ValidateCurves() {
		tfxU32 index = 0;
		tfxU32 last_index = nodes.size() - 1;
		nodes.ResetIteratorIndex();
		do {
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
		} while (!nodes.EndOfBuckets());
	}

	void tfxGraph::DeleteNode(const tfxAttributeNode &n) {
		nodes.erase(&n);
	}

	void tfxGraph::DeleteNodeAtFrame(float frame) {
		nodes.ResetIteratorIndex();
		do {
			for (auto &n : nodes) {
				if (n.frame == frame) {
					nodes.erase(&n);
					return;
				}
			}
		} while (!nodes.EndOfBuckets());
		
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
		case tfxGraphPreset::tfxTranslationPreset:
			min = { 0.f, -4000.f }; max = { tfxMAX_FRAME, 4000.f };
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
		case tfxGraphPreset::tfxTranslationPreset:
			mm = { 0.f, -4000.f, tfxMAX_FRAME, 4000.f };
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
		nodes.ResetIteratorIndex();
		do {
			for (auto &a : nodes) {
				a.index = i++;
			}
		} while (!nodes.EndOfBuckets());
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
		case tfxGraphPreset::tfxTranslationPreset:
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
		case tfxGraphPreset::tfxTranslationPreset:
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
		graph.nodes.ResetIteratorIndex();
		do {
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
		} while (!graph.nodes.EndOfBuckets());
		return lastv;
	}

	float LookupPrecise(tfxGraph &graph, float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		tfxAttributeNode *lastec = nullptr;
		graph.nodes.ResetIteratorIndex();
		do {
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
		} while (!graph.nodes.EndOfBuckets());
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

	bool IsOvertimeGraph(tfxGraphType type) {
		return type >= tfxOvertime_velocity && type != tfxOvertime_noise_resolution && type <= tfxOvertime_noise_resolution;
	}

	bool IsOvertimePercentageGraph(tfxGraphType type) {
		return type >= tfxOvertime_velocity && type != tfxOvertime_velocity_adjuster && type != tfxOvertime_direction && type <= tfxOvertime_noise_resolution;
	}

	bool IsGlobalGraph(tfxGraphType type) {
		return type >= tfxGlobal_life && type <= tfxGlobal_emitter_depth;
	}

	bool IsTransformGraph(tfxGraphType type) {
		return type >= tfxTransform_roll && type <= tfxTransform_translate_z;
	}

	bool IsGlobalPercentageGraph(tfxGraphType type) {
		return type >= tfxGlobal_life && type <= tfxGlobal_splatter;
	}

	bool IsAngleGraph(tfxGraphType type) {
		return (type == tfxTransform_roll || type == tfxTransform_pitch || type == tfxTransform_yaw || type == tfxProperty_emission_pitch || type == tfxProperty_emission_yaw || type == tfxProperty_emission_range ||
			type == tfxProperty_arc_offset || type == tfxProperty_arc_size || type == tfxBase_spin || type == tfxVariation_spin);
	}

	bool IsAngleOvertimeGraph(tfxGraphType type) {
		return type == tfxOvertime_direction;
	}

	bool IsEverythingElseGraph(tfxGraphType type) {
		return !IsOvertimeGraph(type) && !IsOvertimePercentageGraph(type) && !IsGlobalGraph(type) && !IsAngleGraph(type) && !IsOvertimeGraph(type);
	}

	bool HasKeyframes(tfxEffectEmitter &e) {
		assert(e.transform_attributes < e.library->transform_attributes.size());		//Must be a valid keyframes index into the library
		tfxTransformAttributes &keyframes = e.library->transform_attributes[e.transform_attributes];
		tfxU32 size =	keyframes.translation_x.nodes.size() +
						keyframes.translation_y.nodes.size() +
						keyframes.translation_z.nodes.size();
		return size > 0;
	}

	bool HasMoreThanOneKeyframe(tfxEffectEmitter &e) {
		assert(e.transform_attributes < e.library->transform_attributes.size());		//Must be a valid keyframes index into the library
		tfxTransformAttributes &keyframes = e.library->transform_attributes[e.transform_attributes];
		return	keyframes.translation_x.nodes.size() > 1 ||
				keyframes.translation_y.nodes.size() > 1 ||
				keyframes.translation_z.nodes.size() > 1;
	}

	void PushTranslationPoints(tfxEffectEmitter &e, tfxStack<tfxVec3> &points, float frame) {
		assert(e.transform_attributes < e.library->transform_attributes.size());		//Must be a valid keyframes index into the library
		tfxTransformAttributes &keyframes = e.library->transform_attributes[e.transform_attributes];
		tfxVec3 point(	lookup_callback(keyframes.translation_x, frame),
						lookup_callback(keyframes.translation_y, frame),
						lookup_callback(keyframes.translation_z, frame));
		points.push_back(point);
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
		else if (type >= tfxOvertimeStart && type < tfxTransformStart) {
			int ref = type - tfxOvertimeStart;
			return ((tfxGraph*)&library.emitter_attributes[graph_id.graph_id].overtime)[ref];
		}
		else if (type >= tfxTransformStart) {
			int ref = type - tfxTransformStart;
			return ((tfxGraph*)&library.transform_attributes[graph_id.graph_id])[ref];
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

	tfxErrorFlags LoadEffectLibraryPackage(tfxPackage &package, tfxEffectLibrary &lib, void(*shape_loader)(const char *filename, tfxImageData &image_data, void *raw_image_data, int image_size, void *user_data), void *user_data, bool read_only) {

		assert(shape_loader);
		if (!data_types.initialised) data_types.Init();
		lib.Clear();
		if (tfxIcospherePoints[0].current_size == 0) {
			MakeIcospheres();
		}

		tfxEntryInfo *data = package.GetFile("data.txt");
		tfxEntryInfo *stats_struct = package.GetFile("stats.struct");
		tmpStack(tfxEffectEmitter, effect_stack);
		tfxErrorFlags error = 0;

		int context = 0;
		int uid = 0;
		tfxU32 current_global_graph = 0;

		lib.property_array_allocator = CreateArenaManager(tfxMegabyte(4));
		lib.InitEmitterProperties();

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
					effect.library = &lib;
					effect.type = tfxEffectEmitterType::tfxFolder;
					effect.info_index = lib.AddEffectEmitterInfo();
					lib.GetInfo(effect).uid = uid++;
					effect_stack.push_back(effect);
				}
				else if (context == tfxStartStage) {
					tfxEffectEmitter effect;
					effect.library = &lib;
					effect.type = tfxEffectEmitterType::tfxStage;
					effect.info_index = lib.AddEffectEmitterInfo();
					lib.AddPreviewCameraSettings(effect);
					effect.transform_attributes = lib.AddKeyframes();
					lib.GetInfo(effect).uid = uid++;
					effect_stack.push_back(effect);
				}
				else if (context == tfxStartEffect) {
					tfxEffectEmitter effect;
					effect.library = &lib;
					effect.info_index = lib.AddEffectEmitterInfo();
					effect.property_index = lib.AddEmitterProperties();
					effect.transform_attributes = lib.AddKeyframes();
					if (effect_stack.size() <= 1) { //Only root effects get the global graphs
						lib.AddEffectGraphs(effect);
						effect.ResetEffectGraphs(false, false);
						current_global_graph = effect.global;
					}
					lib.AddTransformGraphs(effect);
					effect.ResetTransformGraphs(false, false);
					effect.type = tfxEffectEmitterType::tfxEffectType;
					lib.AddAnimationSettings(effect);
					lib.AddPreviewCameraSettings(effect);
					lib.GetInfo(effect).uid = uid++;
					effect_stack.push_back(effect);

				}
				else if (context == tfxStartEmitter) {
					tfxEffectEmitter emitter;
					emitter.library = &lib;
					emitter.info_index = lib.AddEffectEmitterInfo();
					emitter.property_index = lib.AddEmitterProperties();
					emitter.transform_attributes = lib.AddKeyframes();
					lib.AddEmitterGraphs(emitter);
					lib.AddTransformGraphs(emitter);
					emitter.type = tfxEffectEmitterType::tfxEmitterType;
					emitter.ResetEmitterGraphs(false, false);
					emitter.ResetTransformGraphs(false, false);
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
				else if (context == tfxStartGraphs && effect_stack.back().type == tfxEmitterType) {
					AssignGraphData(effect_stack.back(), pair);
				}
				else if (context == tfxStartGraphs && effect_stack.back().type == tfxEffectType) {
					if (effect_stack.size() <= 2)
						AssignGraphData(effect_stack.back(), pair);
				}
				else if (context == tfxStartStage) {
					if (data_types.names_and_types.ValidName(pair[0])) {
						switch (data_types.names_and_types.At(pair[0])) {
						case tfxUint:
							AssignStageProperty(effect_stack.back(), pair[0], (tfxU32)atoi(pair[1].c_str()));
							break;
						case tfxFloat:
							AssignStageProperty(effect_stack.back(), pair[0], (float)atof(pair[1].c_str()));
							break;
						case tfxSInt:
							AssignStageProperty(effect_stack.back(), pair[0], atoi(pair[1].c_str()));
							break;
						case tfxBool:
							AssignStageProperty(effect_stack.back(), pair[0], (bool)(atoi(pair[1].c_str())));
							break;
						case tfxString:
							AssignStageProperty(effect_stack.back(), pair[0], pair[1]);
							break;
						}
					}
					else {
						error |= tfxErrorCode_some_data_not_loaded;
					}
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
					if (effect_stack.parent().type == tfxStage && effect_stack.parent().GetInfo().sub_effectors.size() == 0) {
						tfxEffectPropertyFlags tmp = effect_stack.parent().property_flags;
						if (effect_stack.back().Is3DEffect()) {
							effect_stack.parent().property_flags |= tfxEmitterPropertyFlags_is_3d;
						}
						tmp = effect_stack.parent().property_flags;
					}
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

			if (context == tfxEndStage) {
				assert(effect_stack.size() == 1);			//Stages should not be contained within anything
				lib.effects.push_back(effect_stack.back());
				effect_stack.pop();
			}

		}

		package.Free();

		if (uid >= 0) {
			//Effects were loaded so let's compile them
			lib.CompileAllGraphs();
			lib.ReIndex();
			if (first_shape_index != -1)
				lib.UpdateParticleShapeReferences(lib.effects, first_shape_index);
			lib.UpdateEffectPaths();
			lib.UpdateComputeNodes();
			lib.SetMinMaxData();
		}
		lib.uid = uid;

		return error;
	}

	tfxErrorFlags LoadEffectLibraryPackage(const char *filename, tfxEffectLibrary &lib, void(*shape_loader)(const char* filename, tfxImageData &image_data, void *raw_image_data, int image_size, void *user_data), void *user_data, bool read_only) {

		tfxErrorFlags error = 0;

		tfxPackage package;
		error = LoadPackage(filename, package);
		if (error != 0)
			return error;
		error = LoadEffectLibraryPackage(package, lib, shape_loader, user_data, read_only);

		return error;
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

	tfxParticleManager::~tfxParticleManager() {
		FreeSoABuffer(&effect_buffers);
		FreeSoABuffer(&emitter_buffers);
		particle_array_allocator.FreeAll();
	}

	tfxU32 tfxParticleManager::AddEffect(tfxEffectEmitter &effect, int buffer, int hierarchy_depth, bool is_sub_emitter, float add_delayed_spawning) {
		tfxPROFILE;
		assert(effect.type == tfxEffectType);
		if (flags & tfxEffectManagerFlags_use_compute_shader && highest_compute_controller_index >= max_compute_controllers && free_compute_controllers.empty())
			return tfxINVALID;
		unsigned int parent_index = GetEffectSlot();
		if(parent_index == tfxINVALID)
			return tfxINVALID;
		if (!is_sub_emitter) {
			effects.highest_particle_age[parent_index] = tfxFRAME_LENGTH * 3.f;
		}
		effects.global_attributes[parent_index] = effect.global;
		effects.transform_attributes[parent_index] = effect.transform_attributes;
		effects.age[parent_index] = -add_delayed_spawning;
		effects.state_flags[parent_index] = 0;
		effects.frame[parent_index] = 0.f;
		effects.property_flags[parent_index] = effect.property_flags;
		effects.local_position[parent_index] = tfxVec3();
		effects.timeout[parent_index] = 100.f;
		effects.library[parent_index] = effect.library;
		effects.parent_particle_index[parent_index] = tfxINVALID;
		effects.info_index[parent_index] = effect.info_index;
		effects.properties_index[parent_index] = effect.property_index;
		effects.timeout_counter[parent_index] = 0;
		effects_in_use[hierarchy_depth][buffer].push_back(parent_index);
		sort_passes = tfxMax(effect.sort_passes, sort_passes);
		if (!effect.Is3DEffect()) {
			flags &= ~tfxEffectManagerFlags_3d_effects;
		}
		tfxEmitterPropertiesSoA &properties = effect.library->emitter_properties;
		for (auto &e : effect.GetInfo().sub_effectors) {
			if (!FreeEffectCapacity())
				break;
			if (e.property_flags & tfxEmitterPropertyFlags_enabled) {
				unsigned int index = GetEmitterSlot();
				if (index == tfxINVALID)
					break;
				emitters.particles_index[index] = tfxINVALID;
				emitters_in_use[hierarchy_depth][buffer].push_back(index);
				emitters.parent_index[index] = parent_index;
				if (emitters.particles_index[index] == tfxINVALID && flags & tfxEffectManagerFlags_unordered) {
					if(!is_sub_emitter)
						emitters.particles_index[index] = GrabParticleLists(*this, e.path_hash, 100);
				}
				else {
					emitters.particles_index[index] = properties.layer[e.property_index];
				}
				emitters.path_hash[index] = e.path_hash;
				emitters.info_index[index] = e.info_index;
				emitters.properties_index[index] = e.property_index;
				emitters.emitter_attributes[index] = e.emitter_attributes;
				emitters.transform_attributes[index] = e.transform_attributes;
				emitters.library[index] = effect.library;
				emitters.delay_spawning[index] = properties.delay_spawning[e.property_index];
				emitters.age[index] = 0.f;
				emitters.frame[index] = 0.f;
				emitters.local_position[index] = tfxVec3();
				emitters.grid_coords[index] = tfxVec3();
				emitters.grid_direction[index] = tfxVec3();
				emitters.property_flags[index] = e.property_flags;
				emitters.image_size[index] = properties.image[e.property_index]->image_size;
				emitters.image_frame_rate[index] = properties.image[e.property_index]->animation_frames > 1 && e.property_flags & tfxEmitterPropertyFlags_animate ? properties.frame_rate[e.property_index] : 0.f;
				emitters.end_frame[index] = properties.end_frame[e.property_index];
				emitters.angle_offsets[index] = properties.angle_offsets[e.property_index];
				emitters.timeout[index] = 100.f;
				emitters.amount_remainder[index] = 0.f;
				emitters.qty_step_size[index] = 0.f;
				emitters.timeout_counter[index] = 0;
				emitters.emitter_size[index] = 0.f;
				emitters.hierarchy_depth[index] = hierarchy_depth;
				//----Handle
				if (e.property_flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
					emitters.image_handle[index] = tfxVec2(0.5f, 0.5f);
				}
				else {
					emitters.image_handle[index] = properties.image_handle[e.property_index];
				}
				tfxEmitterStateFlags &state_flags = emitters.state_flags[index];
				const tfxEmitterStateFlags &parent_state_flags = emitters.state_flags[parent_index];
				
				state_flags = tfxEmitterStateFlags_no_tween_this_update;
				state_flags &= ~tfxEmitterStateFlags_retain_matrix;
				state_flags |= parent_state_flags & tfxEmitterStateFlags_no_tween;
				state_flags |= e.property_flags & tfxEmitterPropertyFlags_single && !(flags & tfxEffectManagerFlags_disable_spawning) ? tfxEmitterStateFlags_is_single : 0;
				state_flags |= e.property_flags & tfxEmitterPropertyFlags_base_uniform_size;
				state_flags |= (properties.emission_type[e.property_index] != tfxLine && !(e.property_flags & tfxEmitterPropertyFlags_edge_traversal)) || properties.emission_type[e.property_index] == tfxLine && !(e.property_flags & tfxEmitterPropertyFlags_edge_traversal) ? tfxEmitterStateFlags_not_line : 0;
				state_flags |= e.property_flags & tfxEmitterPropertyFlags_random_color;
				state_flags |= e.property_flags & tfxEmitterPropertyFlags_lifetime_uniform_size;
				state_flags |= properties.angle_settings[e.property_index] != tfxAngleSettingFlags_align_roll && !(e.property_flags & tfxEmitterPropertyFlags_relative_angle) ? tfxEmitterStateFlags_can_spin : 0;
				state_flags |= properties.angle_settings[e.property_index] == tfxAngleSettingFlags_align_roll ? tfxEmitterStateFlags_align_with_velocity : 0;
				state_flags |= properties.emission_type[e.property_index] == tfxLine && e.property_flags & tfxEmitterPropertyFlags_edge_traversal ? tfxEmitterStateFlags_is_line_traversal : 0;
				state_flags |= e.property_flags & tfxEmitterPropertyFlags_play_once;
				state_flags |= properties.end_behaviour[e.property_index] == tfxLoop ? tfxEmitterStateFlags_loop : 0;
				state_flags |= properties.end_behaviour[e.property_index] == tfxKill ? tfxEmitterStateFlags_kill : 0;

				if (is_sub_emitter) {
					state_flags |= tfxEmitterStateFlags_is_sub_emitter;
				}
				else {
					emitters.highest_particle_age[index] = tfxFRAME_LENGTH * 2.f;
				}

				if (effect.Is3DEffect()) {
					if (e.property_flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type[e.property_index] == tfxLine) {
						emitters.transform_particle_callback3d[index] = TransformParticlePositionRelativeLine3d;
					}
					else if (e.property_flags & tfxEmitterPropertyFlags_relative_position) {
						emitters.transform_particle_callback3d[index] = TransformParticlePositionRelative3d;
					}
					else if (e.property_flags & tfxEmitterPropertyFlags_relative_angle) {
						emitters.transform_particle_callback3d[index] = TransformParticlePositionAngle3d;
					}
					else {
						emitters.transform_particle_callback3d[index] = TransformParticlePosition3d;
					}
				}
				else {
					if (e.property_flags & tfxEmitterPropertyFlags_relative_position && (e.property_flags & tfxEmitterPropertyFlags_relative_angle || (e.property_flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type[e.property_index] == tfxLine))) {
						emitters.transform_particle_callback2d[index] = TransformParticlePositionRelativeLine;
					}
					else if (e.property_flags & tfxEmitterPropertyFlags_relative_position) {
						emitters.transform_particle_callback2d[index] = TransformParticlePositionRelative;
					}
					else if (e.property_flags & tfxEmitterPropertyFlags_relative_angle) {
						emitters.transform_particle_callback2d[index] = TransformParticlePositionAngle;
					}
					else {
						emitters.transform_particle_callback2d[index] = TransformParticlePosition;
					}
				}

				/*if (flags & tfxEffectManagerFlags_use_compute_shader && e.GetInfo().sub_effectors.empty()) {
					int free_slot = AddComputeController();
					if (free_slot != -1) {
						emitter.compute_slot_id = free_slot;
						emitter.property_flags |= tfxEmitterPropertyFlags_is_bottom_emitter;
						tfxComputeController &c = *(static_cast<tfxComputeController*>(compute_controller_ptr) + free_slot);
						c.flags = 0;
						c.flags |= emitter.property_flags & tfxParticleControlFlags_random_color;
						c.flags |= emitter.property_flags & tfxParticleControlFlags_relative_position;
						c.flags |= emitter.property_flags & tfxParticleControlFlags_relative_angle;
						c.flags |= emitter.property_flags & tfxParticleControlFlags_edge_traversal;
						c.flags |= emitter.property_flags & tfxParticleControlFlags_base_uniform_size;
						c.flags |= emitter.property_flags & tfxParticleControlFlags_lifetime_uniform_size;
						c.flags |= emitter.property_flags & tfxParticleControlFlags_reverse_animation;
						c.flags |= emitter.property_flags & tfxParticleControlFlags_play_once;
						c.flags |= ((properties.emission_type[emitter.property_index] == tfxPoint) << 3);
						c.flags |= ((properties.emission_type[emitter.property_index] == tfxArea) << 4);
						c.flags |= ((properties.emission_type[emitter.property_index] == tfxLine) << 5);
						c.flags |= ((properties.emission_type[emitter.property_index] == tfxEllipse) << 6);
						c.flags |= ((properties.end_behaviour[emitter.property_index] == tfxLoop) << 7);
						c.flags |= ((properties.end_behaviour[emitter.property_index] == tfxKill) << 8);
						c.flags |= ((properties.end_behaviour[emitter.property_index] == tfxLetFree) << 9);
						properties.compute_flags[emitter.property_index] = c.flags;
						if (emitter.property_flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
							c.image_handle = tfxVec2(0.5f, 0.5f);
						}
						else {
							c.image_handle = properties.image_handle[emitter.property_index];
						}
						c.image_handle = properties.image_handle[emitter.property_index];
					}
				}*/
			}
		}
		emitters.state_flags[parent_index] |= tfxEmitterStateFlags_no_tween_this_update; 
		return parent_index;
	}

	tfxU32 tfxParticleManager::AddEffect(tfxEffectTemplate &effect) {
		return AddEffect(effect.effect, current_ebuff, 0);
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

	tfxComputeParticle& tfxParticleManager::GrabComputeParticle(unsigned int layer) {
		assert(new_compute_particle_ptr);		//Use must assign the compute ptr to point to an area in memory where you can stage new particles for uploading to the GPU - See ResetComputePtr
		return *(static_cast<tfxComputeParticle*>(new_compute_particle_ptr) + new_compute_particle_index++);
	}

	void tfxParticleManager::FreeParticleList(tfxU32 index) {
		if (free_particle_lists.ValidKey(emitters.path_hash[index])) {
			free_particle_lists.At(emitters.path_hash[index]).push_back(emitters.particles_index[index]);
		}
		else {
#ifdef tfxTRACK_MEMORY
			tfxvec<tfxU32> new_indexes(tfxCONSTRUCTOR_VEC_INIT("FreeParticleBanks::new_indexes"));
#else
			tfxvec<tfxU32> new_indexes;
#endif
			new_indexes.push_back(emitters.particles_index[index]);
			free_particle_lists.Insert(emitters.path_hash[index], new_indexes);
		}
	}

	void tfxParticleManager::Update() {
		tfxPROFILE;
		tfxCompleteAllWork(&work_queue);
		new_compute_particle_index = 0;

		unsigned int next_buffer = !current_ebuff;

		memset(new_particles_index_start, tfxMAX_UINT, sizeof(tfxU32) * tfxLAYERS);
		memset(sprite_index_point, 0, sizeof(tfxU32) * tfxLAYERS);

		if (flags & tfxEffectManagerFlags_ordered_by_age) {
			for (tfxEachLayer) {
				sprite_index_point[layer] = particle_array_buffers[layer].current_size;
			}
		}
		else if (flags & tfxEffectManagerFlags_order_by_depth) {
			for (tfxEachLayer) {
				int layer_offset = layer * 2;
				int current_buffer_index = current_pbuff + layer_offset;
				sprite_index_point[layer] = particle_array_buffers[current_buffer_index].current_size;
			}
		}

		for (tfxEachLayer) {
			sprites2d[layer].clear();
			sprites3d[layer].clear();
		}
		
		tfxU32 effects_start_size[tfxMAXDEPTH];
		tfxU32 emitter_start_size[tfxMAXDEPTH];
		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			effects_start_size[depth] = effects_in_use[depth][current_ebuff].current_size;
			emitter_start_size[depth] = emitters_in_use[depth][current_ebuff].current_size;
		}

		tmpMTStack(tfxSpawnWorkEntry, spawn_work);
		//Loop over all the effects and emitters, depth by depth, and add spawn jobs to the worker queue
		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			effects_in_use[depth][next_buffer].clear();
			emitters_in_use[depth][next_buffer].clear();

			for (int i = 0; i != effects_start_size[depth]; ++i) {
				tfxU32 current_index = effects_in_use[depth][current_ebuff][i];
				float &timeout_counter = effects.timeout_counter[current_index];

				UpdatePMEffect(*this, current_index);
				if (timeout_counter <= effects.timeout[current_index]) {
					effects_in_use[depth][next_buffer].push_back(current_index);
				}
				else {
					free_effects.push_back(current_index);
				}
			}

			for (int i = 0; i != emitter_start_size[depth]; ++i) {
				tfxSpawnWorkEntry *spawn_work_entry = &spawn_work.next();
				tfxU32 current_index = emitters_in_use[depth][current_ebuff][i];
				spawn_work_entry->depth = depth;
				spawn_work_entry->emitter_index = current_index;
				spawn_work_entry->next_buffer = next_buffer;
				spawn_work_entry->properties = &emitters.library[current_index]->emitter_properties;
				spawn_work_entry->sub_effects = &emitters.library[current_index]->effect_infos[emitters.info_index[current_index]].sub_effectors;
				spawn_work_entry->amount_to_spawn = 0;
				spawn_work_entry->end_index = 0;

				float &timeout_counter = emitters.timeout_counter[current_index];

				UpdatePMEmitter(*this, spawn_work_entry);
				if (timeout_counter <= emitters.timeout[current_index]) {
					emitters_in_use[depth][next_buffer].push_back(current_index);
				}
				else {
					free_emitters.push_back(current_index);
					//if (flags & tfxEffectManagerFlags_use_compute_shader && emitters.property_flags[current_index] & tfxEmitterPropertyFlags_is_bottom_emitter)
						//FreeComputeSlot(emitters.compute_slot_id[current_index]);
					if (flags & tfxEffectManagerFlags_unordered) {
						FreeParticleList(current_index);
					}
				}
			}
		}

		tfxCompleteAllWork(&work_queue);

		for (auto &work_entry : spawn_work) {
			tfxU32 index = work_entry.emitter_index;
			emitters.highest_particle_age[index] = std::fmaxf(emitters.highest_particle_age[index], work_entry.highest_particle_age);
			effects.highest_particle_age[emitters.parent_index[index]] = emitters.highest_particle_age[index] + tfxFRAME_LENGTH;
		}
		spawn_work.free();

		if (flags & tfxEffectManagerFlags_unordered) {
			for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
				{
					tmpMTStack(tfxControlWorkEntry, work);
					for (int index : emitters_in_use[depth][next_buffer]) {
						tfxSoABuffer &bank = particle_array_buffers[emitters.particles_index[index]];
						int particles_to_update = bank.current_size;
						tfxU32 running_start_index = 0;
						while (particles_to_update > 0) {
							tfxControlWorkEntry &work_entry = work.next();
							work_entry.properties = &emitters.library[index]->emitter_properties;
							work_entry.emitter_index = index;
							work_entry.start_index = running_start_index;
							work_entry.end_index = particles_to_update > mt_batch_size ? running_start_index + mt_batch_size :  running_start_index + particles_to_update;
							particles_to_update -= mt_batch_size;
							running_start_index += mt_batch_size;
							if (emitters.property_flags[index] & tfxEmitterPropertyFlags_is_3d)
								ControlParticles3d(*this, index, work_entry);
							else
								ControlParticles2d(*this, index, work_entry);
						}
					}
					tfxCompleteAllWork(&work_queue);
					work.free();
				}
				{
					tmpMTStack(tfxParticleAgeWorkEntry, work);
					for (int index : emitters_in_use[depth][next_buffer]) {
						tfxSoABuffer &bank = particle_array_buffers[emitters.particles_index[index]];
						tfxParticleAgeWorkEntry &work_entry = work.next();
						work_entry.properties = &emitters.library[index]->emitter_properties;
						work_entry.start_index = bank.current_size - 1;
						work_entry.emitter_index = index;
						work_entry.pm = this;
#if tfxMULTITHREADED
						tfxAddWorkQueueEntry(&work_queue, &work_entry, ControlParticleAge);
#else
						ControlParticleAge(&work_queue, &work_entry);
#endif
					}
					tfxCompleteAllWork(&work_queue);
					work.free();
				}
			}
		}
		else if (flags & tfxEffectManagerFlags_ordered_by_age) {

			//if(flags & tfxEffectManagerFlags_3d_effects && flags & tfxEffectManagerFlags_order_by_depth) 
				//ControlParticlesDepthOrdered3d(*this);
			//else if(flags & tfxEffectManagerFlags_3d_effects)
				//ControlParticlesOrdered3d(*this);
			//else

			{
				tmpMTStack(tfxControlWorkEntryOrdered, work);
				for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
					int particles_to_update = particle_array_buffers[layer].current_size;
					tfxU32 running_start_index = 0;
					while(particles_to_update > 0) {
						tfxControlWorkEntryOrdered &work_entry = work.next();
						work_entry.pm = this;
						work_entry.sprite_layer = layer;
						work_entry.current_buffer_index = layer;
						work_entry.start_index = running_start_index;
						work_entry.end_index = particles_to_update > mt_batch_size ? running_start_index + mt_batch_size :  running_start_index + particles_to_update;
						work_entry.amount_to_update = work_entry.end_index - work_entry.start_index;
						particles_to_update -= mt_batch_size;
						running_start_index += mt_batch_size;
						if(flags & tfxEffectManagerFlags_3d_effects)
							ControlParticlesOrdered3d(*this, work_entry);
						else
							ControlParticlesOrdered2d(*this, work_entry);
					}
				}
		
				tfxCompleteAllWork(&work_queue);
				work.free();
			}
			for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
				tfxControlWorkEntryOrdered &work_entry = ordered_age_work_entry[layer];
				work_entry.amount_to_update = particle_array_buffers[layer].current_size;
				if (work_entry.amount_to_update == 0)
					continue;
				work_entry.pm = this;
				work_entry.sprite_layer = layer;
				work_entry.current_buffer_index = layer;
				work_entry.start_index = work_entry.amount_to_update - 1;
				work_entry.end_index = 0;
#if tfxMULTITHREADED
				tfxAddWorkQueueEntry(&work_queue, &work_entry, ControlParticleOrderedAge);
#else
				ControlParticleOrderedAge(&work_queue, &work_entry);
#endif
			}

		}
		else if (flags & tfxEffectManagerFlags_order_by_depth) {
			int next_particle_buffer = !current_pbuff;

			for (tfxEachLayer) {
				int layer_offset = layer * 2;
				int next_buffer_index = next_particle_buffer + layer_offset;
				int current_buffer_index = current_pbuff + layer_offset;
				ClearSoABuffer(&particle_array_buffers[next_buffer_index]);
				tfxControlWorkEntryOrdered &work_entry = ordered_age_work_entry[current_buffer_index];
				work_entry.amount_to_update = particle_array_buffers[current_buffer_index].current_size;
				if (work_entry.amount_to_update == 0)
					continue;
				//QuickSortSoAParticles(particle_arrays[current_buffer_index], new_particles_index_start[layer] + 1, particle_array_buffers[current_buffer_index].current_size - 1);
				work_entry.pm = this;
				work_entry.sprite_layer = layer;
				work_entry.current_buffer_index = current_buffer_index;
				work_entry.next_buffer_index = next_buffer_index;
				work_entry.start_index = work_entry.amount_to_update - 1;
				work_entry.end_index = 0;
				ControlParticleOrderedDepth(&work_queue, &work_entry);
			}
			current_pbuff = next_particle_buffer;

			{
				tmpMTStack(tfxControlWorkEntryOrdered, work);
				for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
					int layer_offset = layer * 2;
					int current_buffer_index = current_pbuff + layer_offset;
					int particles_to_update = particle_array_buffers[current_buffer_index].current_size;
					sprites3d[layer].current_size = particles_to_update;
					tfxU32 running_start_index = 0;
					while (particles_to_update > 0) {
						tfxControlWorkEntryOrdered &work_entry = work.next();
						work_entry.pm = this;
						work_entry.sprite_layer = layer;
						work_entry.current_buffer_index = current_buffer_index;
						work_entry.start_index = running_start_index;
						work_entry.end_index = particles_to_update > mt_batch_size ? running_start_index + mt_batch_size : running_start_index + particles_to_update;
						work_entry.amount_to_update = work_entry.end_index - work_entry.start_index;
						particles_to_update -= mt_batch_size;
						running_start_index += mt_batch_size;
						ControlParticlesOrdered3d(*this, work_entry);
					}
				}

				tfxCompleteAllWork(&work_queue);
				work.free();
				if (flags & tfxEffectManagerFlags_guarantee_order) {
					for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
						int layer_offset = layer * 2;
						int current_buffer_index = current_pbuff + layer_offset;
						InsertionSortSoAParticles(*this, particle_arrays[current_buffer_index], current_buffer_index, sprites3d[layer].current_size);
					}
				}
			}

		}

		//Add Subeffects to the next buffer
		for (int depth = 1; depth != tfxMAXDEPTH; ++depth) {
			for (int i = effects_start_size[depth]; i != effects_in_use[depth][current_ebuff].current_size; ++i) {
				tfxU32 current_index = effects_in_use[depth][current_ebuff][i];
				effects_in_use[depth][next_buffer].push_back(current_index);
			}
			for (int i = emitter_start_size[depth]; i != emitters_in_use[depth][current_ebuff].current_size; ++i) {
				tfxU32 current_index = emitters_in_use[depth][current_ebuff][i];
				if (flags & tfxEffectManagerFlags_unordered) {
					emitters.particles_index[current_index] = GrabParticleLists(*this, emitters.path_hash[current_index], 100);
				}
				emitters_in_use[depth][next_buffer].push_back(current_index);
			}
		}

		current_ebuff = next_buffer;

		flags &= ~tfxEffectManagerFlags_update_base_values;

	}

	void ControlParticlesOrdered2d(tfxParticleManager &pm, tfxControlWorkEntryOrdered &work_entry) {
		tfxPROFILE;

		//-------------------------------------------------------
		//Controll what the particle does over the course of
		//it's lifetime
		//-------------------------------------------------------

		work_entry.sprites2d = &pm.sprites2d[work_entry.sprite_layer];

#if tfxMULTITHREADED
		tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticlePositionOrdered2d);
		tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleSizeOrdered2d);
		tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleColorOrdered2d);
		tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleImageFrameOrdered2d);
#else
		ControlParticlePositionOrdered2d(&pm.work_queue, &work_entry);
		ControlParticleSizeOrdered2d(&pm.work_queue, &work_entry);
		ControlParticleColorOrdered2d(&pm.work_queue, &work_entry);
		ControlParticleImageFrameOrdered2d(&pm.work_queue, &work_entry);
#endif
	}

	void ControlParticlesOrdered3d(tfxParticleManager &pm, tfxControlWorkEntryOrdered &work_entry) {
		tfxPROFILE;

		//-------------------------------------------------------
		//Controll what the particle does over the course of
		//it's lifetime
		//-------------------------------------------------------

		work_entry.sprites3d = &pm.sprites3d[work_entry.sprite_layer];

#if tfxMULTITHREADED
		tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticlePositionOrdered3d);
		tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleSizeOrdered3d);
		tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleColorOrdered3d);
		tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleImageFrameOrdered3d);
#else
		ControlParticlePositionOrdered3d(&pm.work_queue, &work_entry);
		ControlParticleSizeOrdered3d(&pm.work_queue, &work_entry);
		ControlParticleColorOrdered3d(&pm.work_queue, &work_entry);
		ControlParticleImageFrameOrdered3d(&pm.work_queue, &work_entry);
#endif
	}

	void ControlParticleOrderedAge(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntryOrdered *work_entry = static_cast<tfxControlWorkEntryOrdered*>(data);
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[work_entry->current_buffer_index];

		tfxU32 offset = 0;
		for (int i = work_entry->start_index; i >= 0; --i) {
			const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[work_entry->current_buffer_index], i);
			const tfxU32 parent_index = bank.parent_index[index];
			tfxEffectLibrary *library = pm.emitters.library[parent_index];

			const tfxU32 property_index = pm.emitters.properties_index[parent_index];
			const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[parent_index];
			const tfxEmitterStateFlags state_flags = pm.emitters.state_flags[parent_index];
			const tfxU32 single_shot_limit = library->emitter_properties.single_shot_limit[property_index];

			float &age = bank.age[index];
			const float max_age = bank.max_age[index];
			tfxU32 &single_loop_count = bank.single_loop_count[index];
			tfxParticleFlags &flags = bank.flags[index];
			age += tfxFRAME_LENGTH;
			flags |= state_flags & tfxParticleFlags_remove;

			if (flags & tfxParticleFlags_remove || age >= max_age) {
				if (property_flags & tfxEmitterPropertyFlags_single && !(work_entry->pm->flags & tfxEffectManagerFlags_disable_spawning))
					if (++single_loop_count != single_shot_limit) {
						age = 0;
					}
					else {
						flags |= tfxParticleFlags_remove;
					}
				else {
					flags |= tfxParticleFlags_remove;
				}
			}

			if (flags & tfxParticleFlags_remove) {
				offset++;
				if (flags & tfxParticleFlags_has_sub_effects) {
					pm.FreeParticleIndex(bank.particle_index[index]);
				}
			}
			else if (offset > 0) {
				tfxU32 next_index = GetCircularIndex(&work_entry->pm->particle_array_buffers[work_entry->current_buffer_index], i + offset);
				if(flags & tfxParticleFlags_has_sub_effects)
					pm.particle_indexes[bank.particle_index[index]] = MakeParticleID(work_entry->current_buffer_index, next_index);

				bank.parent_index[next_index] = bank.parent_index[index];
				bank.sprite_index[next_index] = bank.sprite_index[index];
				bank.particle_index[next_index] = bank.particle_index[index];
				bank.flags[next_index] = bank.flags[index];
				bank.age[next_index] = bank.age[index];
				bank.max_age[next_index] = bank.max_age[index];
				bank.local_position[next_index] = bank.local_position[index];
				bank.captured_position[next_index] = bank.captured_position[index];
				bank.local_rotations[next_index] = bank.local_rotations[index];
				bank.velocity_normal[next_index] = bank.velocity_normal[index];
				bank.stretch[next_index] = bank.stretch[index];
				bank.weight_acceleration[next_index] = bank.weight_acceleration[index];
				bank.base_weight[next_index] = bank.base_weight[index];
				bank.base_velocity[next_index] = bank.base_velocity[index];
				bank.base_spin[next_index] = bank.base_spin[index];
				bank.noise_offset[next_index] = bank.noise_offset[index];
				bank.noise_resolution[next_index] = bank.noise_resolution[index];
				bank.color[next_index] = bank.color[index];
				bank.intensity[next_index] = bank.intensity[index];
				bank.image_frame[next_index] = bank.image_frame[index];
				bank.base_size[next_index] = bank.base_size[index];
				bank.single_loop_count[next_index] = bank.single_loop_count[index];

			}

		}

		if (offset) {
			Bump(&work_entry->pm->particle_array_buffers[work_entry->current_buffer_index], offset);
		}
	}

	void ControlParticleOrderedDepth(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntryOrdered *work_entry = static_cast<tfxControlWorkEntryOrdered*>(data);
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[work_entry->current_buffer_index];

		tfxU32 new_index = pm.new_particles_index_start[work_entry->sprite_layer];
		if (new_index >= pm.particle_array_buffers[work_entry->current_buffer_index].current_size) {
			new_index = tfxINVALID;
		}

		tfxU32 amount_to_update = tfxMin(pm.new_particles_index_start[work_entry->sprite_layer], work_entry->amount_to_update);

		for(int index = 0; index != amount_to_update; ++index) {
			const float depth = bank.depth[index];

			while (new_index != tfxINVALID && bank.depth[new_index] > depth) {
				tfxU32 next_index = pm.SetNextParticle(work_entry->next_buffer_index, work_entry->current_buffer_index, new_index);
				if (bank.flags[new_index] & tfxParticleFlags_has_sub_effects) {
					pm.particle_indexes[bank.particle_index[new_index]] = next_index;
				}
				new_index++;
				new_index = new_index >= work_entry->amount_to_update ? tfxINVALID : new_index;
			}
			
			const tfxU32 parent_index = bank.parent_index[index];
			tfxEffectLibrary *library = pm.emitters.library[parent_index];

			tfxParticleFlags &flags = bank.flags[index];

			const tfxU32 property_index = pm.emitters.properties_index[parent_index];
			const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[parent_index];
			const tfxEmitterStateFlags state_flags = pm.emitters.state_flags[parent_index];
			const tfxU32 single_shot_limit = library->emitter_properties.single_shot_limit[property_index];

			float &age = bank.age[index];
			const float &max_age = bank.max_age[index];
			tfxU32 &single_loop_count = bank.single_loop_count[index];
			age += tfxFRAME_LENGTH;
			flags |= state_flags & tfxParticleFlags_remove;

			if (flags & tfxParticleFlags_remove || age >= max_age) {
				if (property_flags & tfxEmitterPropertyFlags_single && !(work_entry->pm->flags & tfxEffectManagerFlags_disable_spawning))
					if (++single_loop_count != single_shot_limit) {
						age = 0;
					}
					else {
						flags |= tfxParticleFlags_remove;
					}
				else {
					flags |= tfxParticleFlags_remove;
				}
			}

			if (!(flags & tfxParticleFlags_remove)) {
				tfxU32 next_index = pm.SetNextParticle(work_entry->next_buffer_index, work_entry->current_buffer_index, index);
				if (bank.flags[index] & tfxParticleFlags_has_sub_effects) {
					pm.particle_indexes[bank.particle_index[index]] = next_index;
				}
			}
			else if (flags & tfxParticleFlags_has_sub_effects) {
				pm.FreeParticleIndex(bank.particle_index[index]);
			}

		}

		if (new_index != tfxINVALID) {
			for (tfxU32 index = new_index; index < work_entry->amount_to_update; ++index) {
				tfxU32 next_index = pm.SetNextParticle(work_entry->next_buffer_index, work_entry->current_buffer_index, index);
				if (bank.flags[index] & tfxParticleFlags_has_sub_effects) {
					pm.particle_indexes[bank.particle_index[index]] = next_index;
				}
			}
		}

		if (!(pm.flags & tfxEffectManagerFlags_guarantee_order) && pm.sort_passes > 0) {
			for (tfxU32 sorts = 0; sorts != pm.sort_passes; ++sorts) {
				tfxParticleSoA &next_bank = pm.particle_arrays[work_entry->next_buffer_index];
				for (tfxU32 i = 1; i < work_entry->amount_to_update; ++i) {
					float depth1 = next_bank.depth[i - 1];
					float depth2 = next_bank.depth[i];
					if (depth1 < depth2) {
						SwapSoAParticle(next_bank, i - 1, i);
						if (next_bank.flags[i - 1] & tfxParticleFlags_has_sub_effects) pm.particle_indexes[next_bank.particle_index[i - 1]] = MakeParticleID(work_entry->next_buffer_index, i - 1);
						if (next_bank.flags[i] & tfxParticleFlags_has_sub_effects) pm.particle_indexes[next_bank.particle_index[i]] = MakeParticleID(work_entry->next_buffer_index, i);
					}
				}
			}
		}

	}

	void ControlParticlePositionOrdered2d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntryOrdered *work_entry = static_cast<tfxControlWorkEntryOrdered*>(data);
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[work_entry->current_buffer_index];

		for (tfxU32 i = work_entry->start_index; i != work_entry->end_index; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[work_entry->current_buffer_index], i);
			const tfxU32 parent_index = bank.parent_index[index];
			tfxEffectLibrary *library = pm.emitters.library[parent_index];

			const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[parent_index];
			const tfxVec3 emitter_size = pm.emitters.emitter_size[parent_index];
			const float overal_scale = pm.emitters.overal_scale[parent_index];
			const float angle_offset = pm.emitters.angle_offsets[parent_index].roll;
			const float velocity_adjuster = pm.emitters.velocity_adjuster[parent_index];
			const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[parent_index];

			const float life = bank.age[index] / bank.max_age[index];
			tfxOvertimeAttributes *graphs = &library->emitter_attributes[emitter_attributes].overtime;
			const tfxU32 lookup_frame = static_cast<tfxU32>((life * graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

			const float base_weight = bank.base_weight[index];
			const float base_velocity = bank.base_velocity[index];
			const float base_spin = bank.base_spin[index];
			const float noise_offset = bank.noise_offset[index];
			const float noise_resolution = bank.noise_resolution[index];
			const float angle = bank.velocity_normal[index].x;
			float &weight_acceleration = bank.weight_acceleration[index];
			tfxVec3 &local_position = bank.local_position[index];
			tfxVec3 &captured_position = bank.captured_position[index];
			tfxVec3 &local_rotations = bank.local_rotations[index];
			tfxParticleFlags &flags = bank.flags[index];

			const float lookup_velocity = graphs->velocity.lookup.values[std::min<tfxU32>(lookup_frame, graphs->velocity.lookup.last_frame)] * velocity_adjuster;
			const float lookup_velocity_turbulance = graphs->velocity_turbulance.lookup.values[std::min<tfxU32>(lookup_frame, graphs->velocity_turbulance.lookup.last_frame)];
			const float lookup_direction = graphs->direction.lookup.values[std::min<tfxU32>(lookup_frame, graphs->direction.lookup.last_frame)] + angle;
			const float lookup_noise_resolution = graphs->noise_resolution.lookup.values[std::min<tfxU32>(lookup_frame, graphs->noise_resolution.lookup.last_frame)] * noise_resolution;
			const float lookup_weight = graphs->weight.lookup.values[std::min<tfxU32>(lookup_frame, graphs->weight.lookup.last_frame)];
			const float lookup_spin = graphs->spin.lookup.values[std::min<tfxU32>(lookup_frame, graphs->spin.lookup.last_frame)] * base_spin;

			float direction = 0;

			tfxVec2 mr_vec;
			if (emitter_flags & tfxEmitterStateFlags_not_line) {
				direction = lookup_direction;
			}

			if (lookup_velocity_turbulance) {
				float eps = 0.0001f;
				float eps2 = 0.0001f * 2.f;

				float x = local_position.x / lookup_noise_resolution + noise_offset;
				float y = local_position.y / lookup_noise_resolution + noise_offset;

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
			weight_acceleration += base_weight * lookup_weight * tfxUPDATE_TIME;

			//----Velocity Changes
			tfxVec2 velocity_normal;
			velocity_normal.x = std::sinf(direction);
			velocity_normal.y = -std::cosf(direction);

			tfxVec2 current_velocity = (base_velocity * lookup_velocity) * velocity_normal;
			current_velocity += mr_vec;
			current_velocity.y += weight_acceleration;
			current_velocity *= tfxUPDATE_TIME;

			//----Spin and angle Changes
			float spin = 0;
			if (emitter_flags & tfxEmitterStateFlags_can_spin) {
				spin = lookup_spin;
			}

			//----Rotation
			if (emitter_flags & tfxEmitterStateFlags_align_with_velocity) {
				tfxVec2 vd = current_velocity.IsNill() ? velocity_normal : current_velocity;
				local_rotations.roll = GetVectorAngle(vd.x, vd.y) + angle_offset;
			}
			else {
				local_rotations.roll += spin * tfxUPDATE_TIME;
			}

			local_position += current_velocity * overal_scale;

			//Lines - Reposition if the particle is travelling along a line
			tfxVec2 offset = velocity_normal * emitter_size.y;
			float length = std::fabsf(local_position.y);
			float emitter_length = emitter_size.y;
			bool line_and_kill = (emitter_flags & tfxEmitterStateFlags_is_line_traversal) && (emitter_flags & tfxEmitterStateFlags_kill) && length > emitter_length;
			bool line_and_loop = (emitter_flags & tfxEmitterStateFlags_is_line_traversal) && (emitter_flags & tfxEmitterStateFlags_loop) && length > emitter_length;
			if (line_and_loop) {
				local_position.y -= offset.y;
				flags |= tfxParticleFlags_capture_after_transform;
			}
			else if (line_and_kill) {
				flags |= tfxParticleFlags_remove;
			}

		}

		tfxU32 running_sprite_index = work_entry->start_index;

		for (tfxU32 i = work_entry->start_index; i != work_entry->end_index; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[work_entry->current_buffer_index], i);
			const tfxU32 parent_index = bank.parent_index[index];

			tfxVec3 &e_captured_position = pm.emitters.captured_position[parent_index];
			tfxVec3 &e_world_position = pm.emitters.world_position[parent_index];
			tfxVec3 &e_world_rotations = pm.emitters.world_rotations[parent_index];
			tfxVec3 &e_handle = pm.emitters.handle[parent_index];
			tfxMatrix4 &e_matrix = pm.emitters.matrix[parent_index];
			tfxVec3 &e_scale = pm.emitters.scale[parent_index];

			const tfxVec3 &local_position = bank.local_position[index];
			const tfxVec3 &local_rotations = bank.local_rotations[index];
			tfxVec3 &captured_position = bank.captured_position[index];
			tfxParticleFlags &flags = bank.flags[index];

			tfxParticleSprite2d &s = (*work_entry->sprites2d)[running_sprite_index++];
			auto transform_particle_callback2d = pm.emitters.transform_particle_callback2d[parent_index];
			if (flags & tfxParticleFlags_capture_after_transform) {
				transform_particle_callback2d(local_position.xy(), local_rotations.roll, s.transform.position, s.transform.rotation, e_world_rotations, e_matrix, e_handle, e_scale, tfxVec3(e_captured_position.x, e_captured_position.y, 0.f));
				captured_position = s.transform.position;
				transform_particle_callback2d(local_position.xy(), local_rotations.roll, s.transform.position, s.transform.rotation, e_world_rotations, e_matrix, e_handle, e_scale, tfxVec3(e_world_position.x, e_world_position.y, 0.f));
				flags &= ~tfxParticleFlags_capture_after_transform;
			}
			else {
				transform_particle_callback2d(local_position.xy(), local_rotations.roll, s.transform.position, s.transform.rotation, e_world_rotations, e_matrix, e_handle, e_scale, tfxVec3(e_world_position.x, e_world_position.y, 0.f));
			}
			s.transform.captured_position = captured_position.xy();
			captured_position = s.transform.position;

		}
	}


	void ControlParticleSizeOrdered2d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntryOrdered *work_entry = static_cast<tfxControlWorkEntryOrdered*>(data);
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[work_entry->current_buffer_index];

		tfxU32 running_sprite_index = work_entry->start_index;

		for (tfxU32 i = work_entry->start_index; i != work_entry->end_index; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[work_entry->current_buffer_index], i);
			const tfxU32 parent_index = bank.parent_index[index];
			tfxEffectLibrary *library = pm.emitters.library[parent_index];

			const float overal_scale = pm.emitters.overal_scale[parent_index];
			const float velocity_adjuster = pm.emitters.velocity_adjuster[parent_index];
			const float stretch = pm.emitters.stretch[parent_index];
			const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[parent_index];
			const tfxVec2 image_size = pm.emitters.image_size[parent_index];
			const tfxVec2 image_handle = pm.emitters.image_handle[parent_index];
			const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[parent_index];

			const float life = bank.age[index] / bank.max_age[index];
			tfxOvertimeAttributes *graphs = &library->emitter_attributes[emitter_attributes].overtime;
			const tfxU32 lookup_frame = static_cast<tfxU32>((life * graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

			const tfxVec2 base_size = bank.base_size[index];
			const float base_velocity = bank.base_velocity[index];
			const float weight_acceleration = bank.weight_acceleration[index];

			float lookup_velocity = graphs->velocity.lookup.values[std::min<tfxU32>(lookup_frame, graphs->velocity.lookup.last_frame)] * velocity_adjuster;
			float lookup_stretch = graphs->stretch.lookup.values[std::min<tfxU32>(lookup_frame, graphs->stretch.lookup.last_frame)];
			float lookup_width = graphs->width.lookup.values[std::min<tfxU32>(lookup_frame, graphs->width.lookup.last_frame)];
			float lookup_height = graphs->height.lookup.values[std::min<tfxU32>(lookup_frame, graphs->height.lookup.last_frame)];

			//----Size Changes
			tfxVec2 scale;
			scale.x = base_size.x * lookup_width;
			if (scale.x < 0.f)
				scale.x = scale.x;

			//----Stretch Changes
			float velocity = std::fabsf(lookup_velocity * base_velocity + weight_acceleration);
			if (emitter_flags & tfxEmitterStateFlags_lifetime_uniform_size) {
				scale.y = (lookup_width * (base_size.y + (velocity * lookup_stretch * stretch))) / image_size.y;
				if (emitter_flags & tfxEmitterPropertyFlags_base_uniform_size && scale.y < scale.x)
					scale.y = scale.x;
			}
			else
				scale.y = (lookup_height * (base_size.y + (velocity * lookup_stretch * stretch))) / image_size.y;

			tfxParticleSprite2d &s = (*work_entry->sprites2d)[running_sprite_index++];
			s.transform.scale = scale * overal_scale;
			s.handle = image_handle;
		}
	}

	void ControlParticleColorOrdered2d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntryOrdered *work_entry = static_cast<tfxControlWorkEntryOrdered*>(data);
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[work_entry->current_buffer_index];

		tfxU32 running_sprite_index = work_entry->start_index;

		for (tfxU32 i = work_entry->start_index; i != work_entry->end_index; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[work_entry->current_buffer_index], i);
			const tfxU32 parent_index = bank.parent_index[index];
			tfxEffectLibrary *library = pm.emitters.library[parent_index];

			const float global_intensity = pm.emitters.intensity[parent_index];
			const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[parent_index];
			const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[parent_index];

			const float life = bank.age[index] / bank.max_age[index];
			tfxOvertimeAttributes *graphs = &library->emitter_attributes[emitter_attributes].overtime;
			const tfxU32 lookup_frame = static_cast<tfxU32>((life * graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

			const float lookup_red = graphs->red.lookup.values[std::min<tfxU32>(lookup_frame, graphs->red.lookup.last_frame)];
			const float lookup_green = graphs->green.lookup.values[std::min<tfxU32>(lookup_frame, graphs->green.lookup.last_frame)];
			const float lookup_blue = graphs->blue.lookup.values[std::min<tfxU32>(lookup_frame, graphs->blue.lookup.last_frame)];
			const float lookup_opacity = graphs->blendfactor.lookup.values[std::min<tfxU32>(lookup_frame, graphs->blendfactor.lookup.last_frame)];
			const float lookup_intensity = graphs->intensity.lookup.values[std::min<tfxU32>(lookup_frame, graphs->intensity.lookup.last_frame)];

			tfxRGBA8 &color = bank.color[index];
			float &intensity = bank.intensity[index];

			//----Color changes
			color.a = unsigned char(255.f * lookup_opacity);
			intensity = lookup_intensity * global_intensity;
			if (!(emitter_flags & tfxEmitterStateFlags_random_color)) {
				color.r = unsigned char(255.f * lookup_red);
				color.g = unsigned char(255.f * lookup_green);
				color.b = unsigned char(255.f * lookup_blue);
			}

			color = tfxRGBA8(color.r, color.g, color.b, color.a);

			tfxParticleSprite2d &s = (*work_entry->sprites2d)[running_sprite_index++];
			s.color = color;
			s.intensity = intensity;
		}

	}

	void ControlParticleImageFrameOrdered2d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntryOrdered *work_entry = static_cast<tfxControlWorkEntryOrdered*>(data);
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[work_entry->current_buffer_index];

		tfxU32 running_sprite_index = work_entry->start_index;

		for (tfxU32 i = work_entry->start_index; i != work_entry->end_index; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[work_entry->current_buffer_index], i);
			const tfxU32 parent_index = bank.parent_index[index];
			tfxEffectLibrary *library = pm.emitters.library[parent_index];

			const tfxU32 property_index = pm.emitters.properties_index[parent_index];
			const float image_frame_rate = pm.emitters.image_frame_rate[parent_index];
			const float end_frame = pm.emitters.end_frame[parent_index];
			const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[parent_index];
			const tfxImageData *image = library->emitter_properties.image[property_index];

			float &image_frame = bank.image_frame[index];
			tfxU32 &sprites_index = bank.sprite_index[index];
			sprites_index = (work_entry->sprite_layer << 28) + running_sprite_index;
			tfxParticleSprite2d &s = (*work_entry->sprites2d)[running_sprite_index++];

			//----Image animation
			image_frame += image_frame_rate * tfxUPDATE_TIME;
			image_frame = (emitter_flags & tfxEmitterStateFlags_play_once) && image_frame > end_frame ? image_frame = end_frame : image_frame;
			image_frame = (emitter_flags & tfxEmitterStateFlags_play_once) && image_frame < 0 ? image_frame = 0 : image_frame;
			image_frame = std::fmodf(image_frame, end_frame + 1);

			s.image_frame = (tfxU32)image_frame;
			s.image_ptr = image->ptr;
		}

	}

	void ControlParticlePositionOrdered3d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntryOrdered *work_entry = static_cast<tfxControlWorkEntryOrdered*>(data);
		tfxParticleManager &pm = *work_entry->pm;
		const tfxU32 layer = work_entry->current_buffer_index;
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[work_entry->current_buffer_index];

		for (tfxU32 i = work_entry->start_index; i != work_entry->end_index; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[work_entry->current_buffer_index], i);
			const tfxU32 parent_index = bank.parent_index[index];

			tfxEffectLibrary *library = pm.emitters.library[parent_index];
			const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[parent_index];
			const tfxVec3 emitter_size = pm.emitters.emitter_size[parent_index];
			const float overal_scale = pm.emitters.overal_scale[parent_index];
			const tfxVec3 angle_offsets = pm.emitters.angle_offsets[parent_index];
			const float velocity_adjuster = pm.emitters.velocity_adjuster[parent_index];
			const float stretch = pm.emitters.stretch[parent_index];
			const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[parent_index];

			const float life = bank.age[index] / bank.max_age[index];
			tfxOvertimeAttributes *graphs = &library->emitter_attributes[emitter_attributes].overtime;
			const tfxU32 lookup_frame = static_cast<tfxU32>((life * graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

			const float base_weight = bank.base_weight[index];
			const float base_velocity = bank.base_velocity[index];
			const float base_spin = bank.base_spin[index];
			const float base_noise_offset = bank.noise_offset[index];
			const float noise_resolution = bank.noise_resolution[index];
			const float angle = bank.velocity_normal[index].x;
			float &weight_acceleration = bank.weight_acceleration[index];
			tfxVec3 &local_position = bank.local_position[index];
			tfxVec3 &captured_position = bank.captured_position[index];
			tfxVec3 &local_rotations = bank.local_rotations[index];
			tfxVec4 &velocity_normal = bank.velocity_normal[index];
			tfxParticleFlags &flags = bank.flags[index];

			const float lookup_velocity = graphs->velocity.lookup.values[std::min<tfxU32>(lookup_frame, graphs->velocity.lookup.last_frame)] * velocity_adjuster;
			const float lookup_velocity_turbulance = graphs->velocity_turbulance.lookup.values[std::min<tfxU32>(lookup_frame, graphs->velocity_turbulance.lookup.last_frame)];
			const float lookup_direction = graphs->direction.lookup.values[std::min<tfxU32>(lookup_frame, graphs->direction.lookup.last_frame)] + angle;
			const float lookup_noise_resolution = graphs->noise_resolution.lookup.values[std::min<tfxU32>(lookup_frame, graphs->noise_resolution.lookup.last_frame)] * noise_resolution;
			const float lookup_weight = graphs->weight.lookup.values[std::min<tfxU32>(lookup_frame, graphs->weight.lookup.last_frame)];
			const float lookup_spin = graphs->spin.lookup.values[std::min<tfxU32>(lookup_frame, graphs->spin.lookup.last_frame)] * base_spin;
			const float lookup_stretch = graphs->stretch.lookup.values[std::min<tfxU32>(lookup_frame, graphs->stretch.lookup.last_frame)];

			float direction = 0.f;
			float mr_angle = 0.f;

			tfxVec3 mr_vec;

			if (lookup_velocity_turbulance) {
				float eps = 0.001f;
				float eps2 = 0.001f * 2.f;

				float noise_offset = base_noise_offset * overal_scale;

				float x = local_position.x / lookup_noise_resolution + noise_offset;
				float y = local_position.y / lookup_noise_resolution + noise_offset;
				float z = local_position.z / lookup_noise_resolution + noise_offset;

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
			weight_acceleration += base_weight * lookup_weight * tfxUPDATE_TIME;

			//----Velocity Changes
			float velocity_scalar = base_velocity * lookup_velocity;
			tfxVec3 current_velocity = velocity_normal.xyz() * velocity_scalar;
			current_velocity += mr_vec;
			current_velocity.y -= weight_acceleration;
			velocity_normal.w = lookup_stretch * stretch;
			current_velocity *= tfxUPDATE_TIME;

			//----Spin and angle Changes
			float spin = 0;
			if (emitter_flags & tfxEmitterStateFlags_can_spin) {
				spin = lookup_spin;
			}

			//----Rotation
			if (flags & tfxEmitterStateFlags_align_with_velocity) {
				tfxVec3 vd = current_velocity.IsNill() ? velocity_normal.xyz() : current_velocity;
				local_rotations.roll = GetVectorAngle(vd.x, vd.y) + angle_offsets.roll;
			}
			else {
				local_rotations.roll += spin * tfxUPDATE_TIME;
			}

			//----Position
			local_position += current_velocity * overal_scale;

			//Lines - Reposition if the particle is travelling along a line
			tfxVec3 offset = velocity_normal.xyz() * emitter_size.y;
			float length = std::fabsf(local_position.y);
			float emitter_length = emitter_size.y;
			bool line_and_kill = (flags & tfxEmitterStateFlags_is_line_traversal) && (flags & tfxEmitterStateFlags_kill) && length > emitter_length;
			bool line_and_loop = (flags & tfxEmitterStateFlags_is_line_traversal) && (flags & tfxEmitterStateFlags_loop) && length > emitter_length;
			if (line_and_loop) {
				local_position.y -= offset.y;
				flags |= tfxParticleFlags_capture_after_transform;
			}
			else if (line_and_kill) {
				flags |= tfxParticleFlags_remove;
			}

		}

		tfxU32 running_sprite_index = work_entry->start_index;

		for (tfxU32 i = work_entry->start_index; i != work_entry->end_index; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[work_entry->current_buffer_index], i);
			const tfxU32 parent_index = bank.parent_index[index];

			tfxVec3 &e_captured_position = pm.emitters.captured_position[parent_index];
			tfxVec3 &e_world_position = pm.emitters.world_position[parent_index];
			tfxVec3 &e_world_rotations = pm.emitters.world_rotations[parent_index];
			tfxVec3 &e_handle = pm.emitters.handle[parent_index];
			tfxMatrix4 &e_matrix = pm.emitters.matrix[parent_index];
			tfxVec3 &e_scale = pm.emitters.scale[parent_index];
			const tfxU32 property_index = pm.emitters.properties_index[parent_index];
			tfxEffectLibrary *library = pm.emitters.library[parent_index];
			const tfxEmissionType emission_type = library->emitter_properties.emission_type[property_index];
			const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[parent_index];
			const tfxVectorAlignType vector_align_type = library->emitter_properties.vector_align_type[property_index];
			const tfxBillboardingOptions billboard_option = library->emitter_properties.billboard_option[property_index];

			const tfxVec3 local_position = bank.local_position[index];
			const tfxVec3 local_rotations = bank.local_rotations[index];
			tfxVec4 &velocity_normal = bank.velocity_normal[index];
			tfxVec3 &captured_position = bank.captured_position[index];
			tfxParticleFlags &flags = bank.flags[index];

			tfxParticleSprite3d &s = (*work_entry->sprites3d)[running_sprite_index++];
			auto transform_particle_callback3d = pm.emitters.transform_particle_callback3d[parent_index];
			if (flags & tfxParticleFlags_capture_after_transform) {
				transform_particle_callback3d(local_position, local_rotations, s.transform.position, s.transform.rotations, e_world_rotations, e_matrix, e_handle, e_scale, e_captured_position);
				captured_position = s.transform.position;
				transform_particle_callback3d(local_position, local_rotations, s.transform.position, s.transform.rotations, e_world_rotations, e_matrix, e_handle, e_scale, e_world_position);
				flags &= ~tfxParticleFlags_capture_after_transform;
			}
			else {
				transform_particle_callback3d(local_position, local_rotations, s.transform.position, s.transform.rotations, e_world_rotations, e_matrix, e_handle, e_scale, e_world_position);
			}

			tfxVec3 alignment_vector;
			bool line = property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == tfxLine;
			if (vector_align_type == tfxVectorAlignType_motion) {
				alignment_vector = s.transform.position - captured_position;
				float l = FastLength(alignment_vector);
				velocity_normal.w *= l * 10.f;
				alignment_vector = FastNormalizeVec(alignment_vector);
			}
			else if (vector_align_type == tfxVectorAlignType_emission && property_flags & tfxEmitterPropertyFlags_relative_position) {
				alignment_vector = mmTransformVector3(e_matrix, velocity_normal.xyz());
			}
			else if (vector_align_type == tfxVectorAlignType_emission) {
				alignment_vector = velocity_normal.xyz();
			}
			else if (vector_align_type == tfxVectorAlignType_emitter) {
				alignment_vector = mmTransformVector(e_matrix, tfxVec4(0.f, 1.f, 0.f, 0.f)).xyz();
			}

			s.transform.captured_position = captured_position;
			alignment_vector.y += 0.002f;	//We don't want a 0 alignment normal
			s.alignment = Pack10bit(alignment_vector, billboard_option & 0x00000003);
			s.stretch = velocity_normal.w;

			captured_position = s.transform.position;
			bank.depth[index] = LengthVec3NoSqR(s.transform.position - pm.camera_position);
		}
	}


	void ControlParticleSizeOrdered3d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntryOrdered *work_entry = static_cast<tfxControlWorkEntryOrdered*>(data);
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = pm.particle_arrays[work_entry->current_buffer_index];

		tfxU32 running_sprite_index = work_entry->start_index;

		for (tfxU32 i = work_entry->start_index; i != work_entry->end_index; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[work_entry->current_buffer_index], i);
			const tfxU32 parent_index = bank.parent_index[index];

			tfxEffectLibrary *library = pm.emitters.library[parent_index];
			const float overal_scale = pm.emitters.overal_scale[parent_index];
			const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[parent_index];
			const tfxVec2 image_handle = pm.emitters.image_handle[parent_index];
			const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[parent_index];
			tfxOvertimeAttributes *graphs = &library->emitter_attributes[emitter_attributes].overtime;
			const tfxU32 width_last_frame = graphs->width.lookup.last_frame;
			const tfxU32 height_last_frame = graphs->height.lookup.last_frame;

			const float life = bank.age[index] / bank.max_age[index];
			const tfxU32 lookup_frame = static_cast<tfxU32>((life * graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

			const tfxVec2 base_size = bank.base_size[index];

			float lookup_width = graphs->width.lookup.values[std::min<tfxU32>(lookup_frame, width_last_frame)];
			float lookup_height = graphs->height.lookup.values[std::min<tfxU32>(lookup_frame, height_last_frame)];

			//----Size Changes
			tfxVec2 scale;
			scale.x = base_size.x * lookup_width;
			if (scale.x < 0.f)
				scale.x = scale.x;

			if (emitter_flags & tfxEmitterStateFlags_lifetime_uniform_size) {
				scale.y = lookup_width * base_size.y;
				if (emitter_flags & tfxEmitterPropertyFlags_base_uniform_size && scale.y < scale.x)
					scale.y = scale.x;
			}
			else
				scale.y = lookup_height * base_size.y;

			tfxParticleSprite3d &s = (*work_entry->sprites3d)[running_sprite_index++];
			s.transform.scale = scale * overal_scale;
			s.handle = image_handle;
		}

	}

	void ControlParticleColorOrdered3d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntryOrdered *work_entry = static_cast<tfxControlWorkEntryOrdered*>(data);
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[work_entry->current_buffer_index];

		tfxU32 running_sprite_index = work_entry->start_index;

		for (tfxU32 i = work_entry->start_index; i != work_entry->end_index; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[work_entry->current_buffer_index], i);
			const tfxU32 parent_index = bank.parent_index[index];
			tfxEffectLibrary *library = pm.emitters.library[parent_index];

			const float global_intensity = pm.emitters.intensity[parent_index];
			const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[parent_index];
			const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[parent_index];

			const float life = bank.age[index] / bank.max_age[index];
			tfxOvertimeAttributes *graphs = &library->emitter_attributes[emitter_attributes].overtime;
			const tfxU32 lookup_frame = static_cast<tfxU32>((life * graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

			const float lookup_opacity = graphs->blendfactor.lookup.values[std::min<tfxU32>(lookup_frame, graphs->blendfactor.lookup.last_frame)];
			const float lookup_intensity = graphs->intensity.lookup.values[std::min<tfxU32>(lookup_frame, graphs->intensity.lookup.last_frame)];

			tfxRGBA8 &color = bank.color[index];
			float &intensity = bank.intensity[index];

			//----Color changes
			color.a = unsigned char(255.f * lookup_opacity);
			intensity = lookup_intensity * global_intensity;
			if (!(emitter_flags & tfxEmitterStateFlags_random_color)) {
				const float lookup_red = graphs->red.lookup.values[std::min<tfxU32>(lookup_frame, graphs->red.lookup.last_frame)];
				const float lookup_green = graphs->green.lookup.values[std::min<tfxU32>(lookup_frame, graphs->green.lookup.last_frame)];
				const float lookup_blue = graphs->blue.lookup.values[std::min<tfxU32>(lookup_frame, graphs->blue.lookup.last_frame)];
				color.r = unsigned char(255.f * lookup_red);
				color.g = unsigned char(255.f * lookup_green);
				color.b = unsigned char(255.f * lookup_blue);
			}

			color = tfxRGBA8(color.r, color.g, color.b, color.a);

			tfxParticleSprite3d &s = (*work_entry->sprites3d)[running_sprite_index++];
			s.color = color;
			s.intensity = intensity;
		}
	}

	void ControlParticleImageFrameOrdered3d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntryOrdered *work_entry = static_cast<tfxControlWorkEntryOrdered*>(data);
		
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = pm.particle_arrays[work_entry->current_buffer_index];

		tfxU32 running_sprite_index = work_entry->start_index;

		for (tfxU32 i = work_entry->start_index; i != work_entry->end_index; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[work_entry->current_buffer_index], i);
			const tfxU32 parent_index = bank.parent_index[index];
			tfxEffectLibrary *library = pm.emitters.library[parent_index];

			const tfxU32 property_index = work_entry->pm->emitters.properties_index[parent_index];
			tfxImageData *image = library->emitter_properties.image[property_index];
			const tfxBillboardingOptions billboard_option = library->emitter_properties.billboard_option[property_index];
			float image_frame_rate = pm.emitters.image_frame_rate[parent_index];
			float end_frame = pm.emitters.end_frame[parent_index];
			tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[parent_index];
			const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[parent_index];

			float &image_frame = bank.image_frame[index];
			tfxU32 &sprites_index = bank.sprite_index[index];
			sprites_index = (work_entry->sprite_layer << 28) + running_sprite_index;
			tfxParticleSprite3d &s = (*work_entry->sprites3d)[running_sprite_index++];

			//----Image animation
			image_frame += image_frame_rate * tfxUPDATE_TIME;
			image_frame = (emitter_flags & tfxEmitterStateFlags_play_once) && image_frame > end_frame ? image_frame = end_frame : image_frame;
			image_frame = (emitter_flags & tfxEmitterStateFlags_play_once) && image_frame < 0 ? image_frame = 0 : image_frame;
			image_frame = std::fmodf(image_frame, end_frame + 1);

			s.image_frame_plus = (billboard_option << 24) + (tfxU32)image_frame;
			s.image_ptr = image->ptr;
		}

	}

	void tfxParticleManager::UpdateParticleOrderOnly() {
		//todo:
		return;
	}

	tfxvec<tfxU32> *tfxParticleManager::GetEffectBuffer() {
		return &effects_in_use[0][current_ebuff];
	}

	void tfxParticleManager::InitFor2d(tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfxParticleManagerModes mode, tfxU32 multi_threaded_batch_size) {
		max_effects = effects_limit;
		mt_batch_size = multi_threaded_batch_size;
		tfxInitialiseWorkQueue(&work_queue);

		if (particle_array_allocator.arenas.current_size == 0) {
			//todo need to be able to adjust the arena size
			particle_array_allocator = CreateArenaManager(tfxMegabyte(2), 8);
			particle_arrays = tfxBucketArray<tfxParticleSoA>(&particle_array_allocator, 32);
		}

		flags = 0;

		if (mode == tfxParticleManagerMode_unordered)
			flags = tfxEffectManagerFlags_unordered;
		else if (mode == tfxParticleManagerMode_ordered_by_age)
			flags = tfxEffectManagerFlags_ordered_by_age;
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

		if (!(flags & tfxEffectManagerFlags_unordered)) {
			for (tfxEachLayer) {
				tfxParticleSoA lists;
				tfxU32 index = particle_arrays.locked_push_back(lists);
				tfxSoABuffer buffer;
				particle_array_buffers.push_back(buffer);
				InitParticleSoA(&particle_array_buffers[index], &particle_arrays.back(), layer_max_values[layer]);
				assert(index == particle_array_buffers.current_size - 1);
			}
		}

#ifdef tfxTRACK_MEMORY
		memcpy(effects[0].name, "ParticleManager::effects0\0", 26);
		memcpy(effects[1].name, "ParticleManager::effects1\0", 26);
#endif

		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			effects_in_use[depth][0].reserve(max_effects);
			effects_in_use[depth][1].reserve(max_effects);

			emitters_in_use[depth][0].reserve(max_effects);
			emitters_in_use[depth][1].reserve(max_effects);
		}

		free_effects.reserve(max_effects);
		InitEffectSoA(&effect_buffers, &effects, max_effects);
		InitEmitterSoA(&emitter_buffers, &emitters, max_effects);
	}

	void tfxParticleManager::InitFor3d(tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfxParticleManagerModes mode, tfxU32 multi_threaded_batch_size) {
		max_effects = effects_limit;
		mt_batch_size = multi_threaded_batch_size;
		tfxInitialiseWorkQueue(&work_queue);

		if (particle_array_allocator.arenas.current_size == 0) {
			//todo need to be able to adjust the arena size
			particle_array_allocator = CreateArenaManager(tfxMegabyte(2), 8);
			particle_arrays = tfxBucketArray<tfxParticleSoA>(&particle_array_allocator, 32);
		}

		flags = 0;

		if (mode == tfxParticleManagerMode_unordered)
			flags = tfxEffectManagerFlags_unordered;
		else if (mode == tfxParticleManagerMode_ordered_by_age)
			flags = tfxEffectManagerFlags_ordered_by_age;
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

		if (flags & tfxEffectManagerFlags_ordered_by_age) {
			for (tfxEachLayer) {
				tfxParticleSoA lists;
				tfxU32 index = particle_arrays.locked_push_back(lists);
				tfxSoABuffer buffer;
				particle_array_buffers.push_back(buffer);
				assert(index == particle_array_buffers.current_size - 1);
				InitParticleSoA(&particle_array_buffers[index], &particle_arrays.back(), max_cpu_particles_per_layer[layer]);
			}
		}
		else if (flags & tfxEffectManagerFlags_order_by_depth) {
			for (tfxEachLayerDB) {
				tfxParticleSoA lists;
				tfxU32 index = particle_arrays.locked_push_back(lists);
				tfxSoABuffer buffer;
				particle_array_buffers.push_back(buffer);
				assert(index == particle_array_buffers.current_size - 1);
				InitParticleSoA(&particle_array_buffers[index], &particle_arrays.back(), max_cpu_particles_per_layer[layer / 2]);
			}
		}

#ifdef tfxTRACK_MEMORY
		memcpy(effects[0].name, "ParticleManager::effects0\0", 26);
		memcpy(effects[1].name, "ParticleManager::effects1\0", 26);
#endif

		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			effects_in_use[depth][0].reserve(max_effects);
			effects_in_use[depth][1].reserve(max_effects);

			emitters_in_use[depth][0].reserve(max_effects);
			emitters_in_use[depth][1].reserve(max_effects);
		}
		free_effects.reserve(max_effects);
		InitEffectSoA(&effect_buffers, &effects, max_effects);
		InitEmitterSoA(&emitter_buffers, &emitters, max_effects);
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
#ifdef tfxTRACK_MEMORY
			tfxring<tfxParticle> particles(tfxCONSTRUCTOR_VEC_INIT("Ordered ring particle list"));
#else
			tfxring<tfxParticle> particles;
#endif

			tfxParticleSoA lists;
			tfxU32 index = particle_arrays.locked_push_back(lists);
			tfxSoABuffer buffer;
			particle_array_buffers.push_back(buffer);
			assert(index == particle_array_buffers.current_size - 1);
			InitParticleSoA(&particle_array_buffers[index], &particle_arrays.back(), max_cpu_particles_per_layer[layer]);
		}
	}

	void tfxParticleManager::Reconfigure(tfxParticleManagerModes mode, tfxU32 req_sort_passes, bool is_3d) {
		if (flags & tfxEffectManagerFlags_unordered && mode != tfxParticleManagerMode_unordered) {
			FreeParticleBanks();
			CreateParticleBanksForEachLayer();
		}
		else if (!(flags & tfxEffectManagerFlags_unordered) && mode == tfxParticleManagerMode_unordered) {
			FreeParticleBanks();
			for (auto &bank : free_particle_lists.data) {
				bank.free_all();
			}
			free_particle_lists.FreeAll();
		}

		tfxEffectManagerFlags dynamic = flags & tfxEffectManagerFlags_dynamic_sprite_allocation;

		if (mode == tfxParticleManagerMode_unordered)
			flags = tfxEffectManagerFlags_unordered;
		else if (mode == tfxParticleManagerMode_ordered_by_depth)
			flags = tfxEffectManagerFlags_order_by_depth;
		else if (mode == tfxParticleManagerMode_ordered_by_depth_guaranteed)
			flags = tfxEffectManagerFlags_order_by_depth | tfxEffectManagerFlags_guarantee_order;
		else if (mode == tfxParticleManagerMode_ordered_by_age)
			flags = tfxEffectManagerFlags_ordered_by_age;

		if (is_3d)
			flags |= tfxEffectManagerFlags_3d_effects;

		flags |= dynamic;

		for (tfxEachLayer) {
			sprites2d[layer].clear();
			sprites3d[layer].clear();
		}

		memset(new_particles_index_start, tfxMAX_UINT, 4 * tfxLAYERS);
		memset(sprite_index_point, 0, 4 * tfxLAYERS);

		ClearSoABuffer(&emitter_buffers);
		ClearSoABuffer(&effect_buffers);
		
		sort_passes = req_sort_passes;
	}

	void tfxParticleManager::InitForBoth(tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfxParticleManagerModes mode, bool dynamic_sprite_allocation) {
		max_effects = effects_limit;

		if (particle_array_allocator.arenas.current_size == 0) {
			//todo need to be able to adjust the arena size
			particle_array_allocator = CreateArenaManager(tfxMegabyte(2), 8);
		}

		flags = 0;

		if (mode == tfxParticleManagerMode_unordered)
			flags = tfxEffectManagerFlags_unordered;
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

		if (!(flags & tfxEffectManagerFlags_unordered)) {
			CreateParticleBanksForEachLayer();
		}

#ifdef tfxTRACK_MEMORY
		memcpy(effects[0].name, "ParticleManager::effects0\0", 26);
		memcpy(effects[1].name, "ParticleManager::effects1\0", 26);
#endif

		flags |= dynamic_sprite_allocation ? tfxEffectManagerFlags_dynamic_sprite_allocation : 0;
	}

	void tfxParticleManager::ClearAll(bool free_memory) {
		tfxCompleteAllWork(&work_queue);
		for (tfxEachLayer) {
			sprites2d[layer].clear();
			sprites3d[layer].clear();
			if (!(flags & tfxEffectManagerFlags_unordered)) {
				ClearSoABuffer(&particle_array_buffers[layer * 2]);
				ClearSoABuffer(&particle_array_buffers[layer * 2 + 1]);
			}
		}
		if (flags & tfxEffectManagerFlags_unordered) {
			if (free_memory) {
				FreeParticleBanks();
			}
			else {
				for (auto &bank : particle_array_buffers) {
					ClearSoABuffer(&bank);
				}
			}
			for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
				for (auto index : emitters_in_use[depth][current_ebuff]) {
					FreeParticleList(index);
				}
			}
		}
		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			effects_in_use[depth][0].clear();
			effects_in_use[depth][1].clear();

			emitters_in_use[depth][0].clear();
			emitters_in_use[depth][1].clear();
		}
		free_effects.clear();
		free_emitters.clear();
		particle_indexes.clear();
		free_particle_indexes.clear();
		ClearSoABuffer(&emitter_buffers);
		ClearSoABuffer(&effect_buffers);
		particle_id = 0;
	}

	void tfxParticleManager::FreeParticleBanks() {
		for (auto &bank : particle_array_buffers) {
			FreeSoABuffer(&bank);
		}
		particle_array_buffers.clear();
		particle_arrays.clear();
	}

	void tfxParticleManager::SoftExpireAll() {
		for (auto index : effects_in_use[0][current_ebuff]) {
			emitters.state_flags[index] |= tfxEmitterStateFlags_stop_spawning;
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

	tfxU32 GrabParticleLists(tfxParticleManager &pm, tfxKey emitter_hash, tfxU32 reserve_amount) {
		if (pm.free_particle_lists.ValidKey(emitter_hash)) {
			tfxvec<tfxU32> &free_banks = pm.free_particle_lists.At(emitter_hash);
			if (free_banks.current_size) {
				return free_banks.pop_back();
			}
		}

		tfxParticleSoA lists;
		tfxU32 index = pm.particle_arrays.locked_push_back(lists);
		tfxSoABuffer buffer;
		pm.particle_array_buffers.push_back(buffer);
		assert(index == pm.particle_array_buffers.current_size - 1);
		InitParticleSoA(&pm.particle_array_buffers[index], &pm.particle_arrays.back(), 100);
		return index;
	}

	void StopSpawning(tfxParticleManager &pm) {
		pm.SoftExpireAll();
	}

	void RemoveAllEffects(tfxParticleManager &pm) {
		pm.ClearAll();
	}

	void AddEffect(tfxParticleManager &pm, tfxEffectEmitter &effect, tfxVec3 position) {
		tfxU32 buffer_index = pm.AddEffect(effect, pm.current_ebuff, 0);
		pm.emitters.local_position[buffer_index] = position;
	}

	void UpdatePMEffect(tfxParticleManager &pm, tfxU32 index, tfxU32 parent_index) {
		tfxPROFILE;

		const tfxU32 property_index = pm.effects.properties_index[index];
		tfxU32 &parent_particle_index = pm.effects.parent_particle_index[index];
		tfxVec3 &captured_position = pm.effects.captured_position[index];
		tfxVec3 &world_position = pm.effects.world_position[index];
		tfxVec3 &local_position = pm.effects.local_position[index];
		tfxVec3 &world_rotations = pm.effects.world_rotations[index];
		tfxVec3 &local_rotations = pm.effects.local_rotations[index];
		tfxVec3 &scale = pm.effects.scale[index];
		tfxVec3 &translation = pm.effects.translation[index];
		tfxMatrix4 &matrix = pm.effects.matrix[index];
		float &frame = pm.effects.frame[index];
		float &age = pm.effects.age[index];
		float &timeout_counter = pm.effects.timeout_counter[index];
		float &highest_particle_age = pm.effects.highest_particle_age[index];
		tfxEmitterPropertyFlags &property_flags = pm.effects.property_flags[index];
		tfxEmitterStateFlags &state_flags = pm.effects.state_flags[index];
		tfxEffectLibrary *library = pm.effects.library[index];

		captured_position = world_position;

		if (pm.lookup_mode == tfxPrecise) {
			frame = age;
		}
		else {
			frame = age / tfxLOOKUP_FREQUENCY;
		}

		UpdateEffectState(pm, index);

		tfxEmitterPropertiesSoA &properties = library->emitter_properties;

		if (parent_particle_index != tfxINVALID) {
			tfxParticleID parent_particle_id = pm.particle_indexes[parent_particle_index];
			if (parent_particle_id != tfxINVALID) {
				tfxU32 bank_index = ParticleBank(parent_particle_id);
				const float overal_scale = pm.effects.overal_scale[index];
				tfxU32 sprite_id = pm.GetParticleSpriteIndex(parent_particle_id);
				tfxU32 sprite_layer = (sprite_id & 0xF0000000) >> 28;
				tfxU32 sprite_index = sprite_id & 0x0FFFFFFF;
				if (property_flags & tfxEmitterPropertyFlags_is_3d)
					TransformEffector3d(world_rotations, local_rotations, world_position, local_position, matrix, pm.sprites3d[sprite_layer][sprite_index].transform, true, property_flags & tfxEmitterPropertyFlags_relative_angle);
				else
					TransformEffector2d(world_rotations, local_rotations, world_position, local_position, matrix, pm.sprites2d[sprite_layer][sprite_index].transform, true, property_flags & tfxEmitterPropertyFlags_relative_angle);

				world_position += properties.emitter_handle[property_index] * overal_scale;
				if (state_flags & tfxEmitterStateFlags_no_tween_this_update || state_flags & tfxEmitterStateFlags_no_tween) {
					captured_position = world_position;
				}
			}
			else {
				parent_particle_index = tfxINVALID;
				state_flags |= tfxEmitterStateFlags_retain_matrix;
				local_position = world_position;
				local_rotations.roll = world_rotations.roll;
				state_flags |= tfxEmitterStateFlags_stop_spawning;
			}
		}
		else {
			if (!(state_flags & tfxEmitterStateFlags_retain_matrix)) {
				const float overal_scale = pm.effects.overal_scale[index];
				world_position = local_position + translation;
				world_position += properties.emitter_handle[property_index] * overal_scale;
				if (property_flags & tfxEmitterPropertyFlags_is_3d) {
					world_rotations = local_rotations;
					tfxMatrix4 roll = mmZRotate(local_rotations.roll);
					tfxMatrix4 pitch = mmXRotate(local_rotations.pitch);
					tfxMatrix4 yaw = mmYRotate(local_rotations.yaw);
					matrix = mmTransform(yaw, pitch);
					matrix = mmTransform(matrix, roll);
				}
				else {
					world_rotations.roll = local_rotations.roll;
					float s = sin(local_rotations.roll);
					float c = cos(local_rotations.roll);
					matrix.Set2(c, s, -s, c);
				}
			}

			if (state_flags & tfxEmitterStateFlags_no_tween_this_update || state_flags & tfxEmitterStateFlags_no_tween) {
				captured_position = world_position;
			}
		}

		age += tfxFRAME_LENGTH;
		if (!(property_flags & tfxEmitterPropertyFlags_single) || properties.single_shot_limit[property_index] > 0) {
			highest_particle_age -= tfxFRAME_LENGTH;
		}

		if (properties.loop_length && age > properties.loop_length[property_index])
			age -= properties.loop_length[property_index];

		if (highest_particle_age <= 0 && age > tfxFRAME_LENGTH * 5.f) {
			timeout_counter++;
		}
		else {
			timeout_counter = 0;
		}

		state_flags &= ~tfxEmitterStateFlags_no_tween_this_update;
	}

	void UpdatePMEmitter(tfxParticleManager &pm, tfxSpawnWorkEntry *spawn_work_entry) {
		tfxPROFILE;
		tfxU32 index = spawn_work_entry->emitter_index;

		const tfxU32 parent_index = pm.emitters.parent_index[index];
		const tfxU32 property_index = pm.emitters.properties_index[index];
		tfxU32 &sprites_count = pm.emitters.sprites_count[index];
		tfxU32 &sprites_index = pm.emitters.sprites_index[index];
		tfxVec3 &captured_position = pm.emitters.captured_position[index];
		tfxVec3 &world_position = pm.emitters.world_position[index];
		tfxVec3 &local_position = pm.emitters.local_position[index];
		tfxVec3 &world_rotations = pm.emitters.world_rotations[index];
		tfxVec3 &local_rotations = pm.emitters.local_rotations[index];
		tfxVec3 &scale = pm.emitters.scale[index];
		tfxVec3 &translation = pm.emitters.translation[index];
		tfxMatrix4 &matrix = pm.emitters.matrix[index];
		float &frame = pm.emitters.frame[index];
		float &age = pm.emitters.age[index];
		float &timeout_counter = pm.emitters.timeout_counter[index];
		float &delay_spawning = pm.emitters.delay_spawning[index];
		float &highest_particle_age = pm.emitters.highest_particle_age[index];
		tfxEmitterPropertyFlags &property_flags = pm.emitters.property_flags[index];
		tfxEmitterStateFlags &state_flags = pm.emitters.state_flags[index];
		tfxU32 &particles_index = pm.emitters.particles_index[index];
		tfxEffectLibrary *library = pm.emitters.library[index];

		captured_position = world_position;

		if (pm.lookup_mode == tfxPrecise) {
			frame = age;
		}
		else {
			frame = age / tfxLOOKUP_FREQUENCY;
		}

		spawn_work_entry->highest_particle_age = pm.emitters.highest_particle_age[index];

		tfxEmitterPropertiesSoA &properties = library->emitter_properties;

		if (parent_index != tfxINVALID) {
			tfxU32 layer = properties.layer[property_index];

			float &parent_timeout_counter = pm.effects.timeout_counter[parent_index];
			const float parent_age = pm.effects.age[parent_index];
			const tfxEmitterStateFlags parent_state_flags = pm.effects.state_flags[parent_index];

			parent_timeout_counter = 0;
			if (parent_age < delay_spawning) {
				return;
			}
			delay_spawning = -tfxFRAME_LENGTH;

			//e.state_flags |= e.parent->state_flags & tfxEmitterStateFlags_stop_spawning;
			state_flags |= parent_state_flags & tfxEmitterStateFlags_no_tween;
			UpdateEmitterState(pm, index, parent_index, pm.effects.spawn_controls[parent_index]);

			bool is_compute = property_flags & tfxEmitterPropertyFlags_is_bottom_emitter && pm.flags & tfxEffectManagerFlags_use_compute_shader;
			tfxU32 amount_spawned = 0;
			tfxU32 max_spawn_count = NewSpritesNeeded(pm, index, parent_index, properties);

			tfxVec3 &parent_captured_position = pm.effects.captured_position[parent_index];
			tfxVec3 &parent_world_position = pm.effects.world_position[parent_index];
			tfxVec3 &parent_world_rotations = pm.effects.world_rotations[parent_index];
			tfxVec3 &parent_scale = pm.effects.scale[parent_index];
			tfxMatrix4 &parent_matrix = pm.effects.matrix[parent_index];

			if (pm.flags & tfxEffectManagerFlags_unordered) {

				if (property_flags & tfxEmitterPropertyFlags_is_3d) {
					tfxring<tfxParticleSprite3d> &sprite_buffer = pm.sprites3d[layer];

					Transform3d(world_rotations, local_rotations, scale, world_position, local_position, translation, matrix, parent_world_rotations, parent_scale, parent_world_position, parent_matrix);

					if (state_flags & tfxEmitterStateFlags_no_tween_this_update || state_flags & tfxEmitterStateFlags_no_tween) {
						captured_position = world_position;
					}

					sprites_count = pm.particle_array_buffers[particles_index].current_size;
					if (pm.flags & tfxEffectManagerFlags_dynamic_sprite_allocation && sprites_count + max_spawn_count > sprite_buffer.free_space()) {
						sprite_buffer.reserve(sprite_buffer._grow_capacity(sprite_buffer.capacity + (sprites_count + max_spawn_count - sprite_buffer.free_space()) + 1));
						sprite_buffer.current_size += sprites_count;
					}
					else if (!(pm.flags & tfxEffectManagerFlags_dynamic_sprite_allocation)) {
						sprites_count = sprites_count > sprite_buffer.free_space() ? sprite_buffer.free_space() : sprites_count;
						sprite_buffer.current_size += sprites_count;
						max_spawn_count = max_spawn_count > sprite_buffer.free_space() ? sprite_buffer.free_space() : max_spawn_count;
					}
					else {
						sprite_buffer.current_size += sprites_count;
					}
					sprites_count += max_spawn_count;
					sprites_index = pm.sprite_index_point[layer];
					pm.sprite_index_point[layer] += sprites_count;
					sprite_buffer.current_size += max_spawn_count;

					if (state_flags & tfxEmitterStateFlags_is_sub_emitter) {
						if (age > 0 && !(pm.flags & tfxEffectManagerFlags_disable_spawning))
							amount_spawned = SpawnWideParticles3d(pm, *spawn_work_entry, max_spawn_count);
					}
					else {
						if (!(pm.flags & tfxEffectManagerFlags_disable_spawning))
							amount_spawned = SpawnWideParticles3d(pm, *spawn_work_entry, max_spawn_count);
					}
					sprite_buffer.current_size -= (max_spawn_count - amount_spawned);
					pm.sprite_index_point[layer] -= (max_spawn_count - amount_spawned);
				}
				else {
					tfxring<tfxParticleSprite2d> &sprite_buffer = pm.sprites2d[layer];

					Transform2d(world_rotations, local_rotations, scale, world_position, local_position, translation, matrix, parent_world_rotations, parent_scale, parent_world_position, parent_matrix);

					if (state_flags & tfxEmitterStateFlags_no_tween_this_update || state_flags & tfxEmitterStateFlags_no_tween) {
						captured_position = world_position;
					}

					sprites_count = pm.particle_array_buffers[particles_index].current_size;
					if (pm.flags & tfxEffectManagerFlags_dynamic_sprite_allocation && sprites_count + max_spawn_count > sprite_buffer.free_space()) {
						sprite_buffer.reserve(sprite_buffer._grow_capacity(sprite_buffer.capacity + (sprites_count + max_spawn_count - sprite_buffer.free_space()) + 1));
						sprite_buffer.current_size += sprites_count;
					}
					else if (!(pm.flags & tfxEffectManagerFlags_dynamic_sprite_allocation)) {
						sprites_count = sprites_count > sprite_buffer.free_space() ? sprite_buffer.free_space() : sprites_count;
						sprite_buffer.current_size += sprites_count;
						max_spawn_count = max_spawn_count > sprite_buffer.free_space() ? sprite_buffer.free_space() : max_spawn_count;
					}
					else {
						sprite_buffer.current_size += sprites_count;
					}

					sprites_count += max_spawn_count;
					sprites_index = pm.sprite_index_point[layer];
					pm.sprite_index_point[layer] += sprites_count;
					sprite_buffer.current_size += max_spawn_count;

					if (state_flags & tfxEmitterStateFlags_is_sub_emitter) {
						if (age > 0 && !(pm.flags & tfxEffectManagerFlags_disable_spawning)) {
							amount_spawned = SpawnWideParticles2d(pm, *spawn_work_entry, max_spawn_count);
						}
					}
					else {
						if (!(pm.flags & tfxEffectManagerFlags_disable_spawning)) {
							amount_spawned = SpawnWideParticles2d(pm, *spawn_work_entry, max_spawn_count);
						}
					}
					sprite_buffer.current_size -= (max_spawn_count - amount_spawned);
					pm.sprite_index_point[layer] -= (max_spawn_count - amount_spawned);
				}
			}
			else {
				sprites_index = pm.sprite_index_point[layer];
				pm.sprite_index_point[layer] += max_spawn_count;
				particles_index = pm.flags & tfxEffectManagerFlags_order_by_depth ? layer * 2 + pm.current_pbuff : layer;
				if (property_flags & tfxEmitterPropertyFlags_is_3d) {
					Transform3d(world_rotations, local_rotations, scale, world_position, local_position, translation, matrix, parent_world_rotations, parent_scale, parent_world_position, parent_matrix);

					if (state_flags & tfxEmitterStateFlags_no_tween_this_update || state_flags & tfxEmitterStateFlags_no_tween) {
						captured_position = world_position;
					}

					tfxring<tfxParticleSprite3d> &sprite_buffer = pm.sprites3d[layer];
					if (pm.flags & tfxEffectManagerFlags_dynamic_sprite_allocation && pm.sprite_index_point[layer] > sprite_buffer.free_space()) {
						sprite_buffer.reserve(sprite_buffer._grow_capacity(sprite_buffer.capacity + (pm.sprite_index_point[layer] - sprite_buffer.free_space()) + 1));
					}
					else if(sprites_index + max_spawn_count > sprite_buffer.capacity) {
						pm.sprite_index_point[layer] = sprite_buffer.capacity;
						max_spawn_count = sprite_buffer.capacity - sprites_index;
					}
					sprite_buffer.current_size = pm.sprite_index_point[layer];

					if (state_flags & tfxEmitterStateFlags_is_sub_emitter) {
						if (age > 0 && property_flags & tfxEmitterPropertyFlags_is_3d && !(pm.flags & tfxEffectManagerFlags_disable_spawning))
							amount_spawned = SpawnWideParticles3d(pm, *spawn_work_entry, max_spawn_count);
					}
					else {
						if (property_flags & tfxEmitterPropertyFlags_is_3d && !(pm.flags & tfxEffectManagerFlags_disable_spawning))
							amount_spawned = SpawnWideParticles3d(pm, *spawn_work_entry, max_spawn_count);
					}
					sprite_buffer.current_size -= (max_spawn_count - amount_spawned);
					pm.sprite_index_point[layer] -= (max_spawn_count - amount_spawned);
				}
				else {
					Transform2d(world_rotations, local_rotations, scale, world_position, local_position, translation, matrix, parent_world_rotations, parent_scale, parent_world_position, parent_matrix);

					if (state_flags & tfxEmitterStateFlags_no_tween_this_update || state_flags & tfxEmitterStateFlags_no_tween) {
						captured_position = world_position;
					}

					tfxring<tfxParticleSprite2d> &sprite_buffer = pm.sprites2d[layer];
					if (pm.flags & tfxEffectManagerFlags_dynamic_sprite_allocation && pm.sprite_index_point[layer] > sprite_buffer.free_space()) {
						sprite_buffer.reserve(sprite_buffer._grow_capacity(sprite_buffer.capacity + (pm.sprite_index_point[layer] - sprite_buffer.free_space()) + 1));
					}
					else if(sprites_index + max_spawn_count > sprite_buffer.capacity) {
						pm.sprite_index_point[layer] = sprite_buffer.capacity;
						max_spawn_count = sprite_buffer.capacity - sprites_index;
					}
					sprite_buffer.current_size = pm.sprite_index_point[layer];

					if (state_flags & tfxEmitterStateFlags_is_sub_emitter) {
						if (age > 0 && !(pm.flags & tfxEffectManagerFlags_disable_spawning))
							amount_spawned = SpawnWideParticles2d(pm, *spawn_work_entry, max_spawn_count);
					}
					else {
						if (!(pm.flags & tfxEffectManagerFlags_disable_spawning))
							amount_spawned = SpawnWideParticles2d(pm, *spawn_work_entry, max_spawn_count);
					}
					sprite_buffer.current_size -= (max_spawn_count - amount_spawned);
					pm.sprite_index_point[layer] -= (max_spawn_count - amount_spawned);
				}
			}

		}
		else {
			assert(false);	//Emitter must have a valid parent index
		}

		age += tfxFRAME_LENGTH;
		if (!(property_flags & tfxEmitterPropertyFlags_single) || properties.single_shot_limit[property_index] > 0) {
			highest_particle_age -= tfxFRAME_LENGTH;
			spawn_work_entry->highest_particle_age -= tfxFRAME_LENGTH;
		}

		if (properties.loop_length && age > properties.loop_length[property_index])
			age -= properties.loop_length[property_index];

		if (highest_particle_age <= 0 && age > tfxFRAME_LENGTH * 5.f) {
			timeout_counter++;
		}
		else {
			timeout_counter = 0;
		}

		state_flags &= ~tfxEmitterStateFlags_no_tween_this_update;
	}

	tfxU32 NewSpritesNeeded(tfxParticleManager &pm, tfxU32 index, tfxU32 parent_index, tfxEmitterPropertiesSoA &properties) {

		const tfxEmitterStateFlags state_flags = pm.emitters.state_flags[index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[index];
		const tfxEmitterStateFlags parent_state_flags = pm.effects.state_flags[parent_index];
		const tfxU32 property_index = pm.effects.properties_index[parent_index];
		const float amount_remainder = pm.emitters.amount_remainder[index];
		float &spawn_quantity = pm.emitters.spawn_quantity[index];

		if (state_flags & tfxEmitterStateFlags_single_shot_done || parent_state_flags & tfxEmitterStateFlags_stop_spawning)
			return 0;
		if (spawn_quantity == 0)
			return 0;

		const tfxVec3 emitter_size = pm.emitters.emitter_size[index];

		float step_size = 0;
		if (!(property_flags & tfxEmitterPropertyFlags_single)) {
			if (property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && (properties.emission_type[property_index] == tfxArea || properties.emission_type[property_index] == tfxEllipse)) {
				if (property_flags & tfxEmitterPropertyFlags_is_3d) {
					float area = std::fmaxf(0.1f, emitter_size.x) * std::fmaxf(0.1f, emitter_size.y) * std::fmaxf(0.1f, emitter_size.z);
					spawn_quantity = (spawn_quantity / 50.f) * area;
				}
				else {
					float area = emitter_size.x * emitter_size.y;
					spawn_quantity = (spawn_quantity / 10000.f) * area;
				}
			}
			else if (property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && properties.emission_type[property_index] == tfxLine) {
				spawn_quantity = (spawn_quantity / 100.f) * emitter_size.y;
			}

			spawn_quantity *= tfxUPDATE_TIME;
			step_size = 1.f / spawn_quantity;
		}
		else {
			step_size = 1.f / spawn_quantity;
		}

		float tween = amount_remainder;
		return tween >= 1.f ? 0 : tfxU32((1.f - amount_remainder) / step_size) + 1;
	}

	tfxU32 SpawnWideParticles2d(tfxParticleManager &pm, tfxSpawnWorkEntry &work_entry, tfxU32 max_spawn_count) {
		const tfxEmitterPropertiesSoA &properties = *work_entry.properties;
		const tfxEmitterStateFlags state_flags = pm.emitters.state_flags[work_entry.emitter_index];
		const tfxEmitterStateFlags parent_state_flags = pm.emitters.state_flags[work_entry.emitter_index];
		const float spawn_quantity = pm.emitters.spawn_quantity[work_entry.emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[work_entry.emitter_index];
		float &qty_step_size = pm.emitters.qty_step_size[work_entry.emitter_index];
		float &amount_remainder = pm.emitters.amount_remainder[work_entry.emitter_index];

		if (state_flags & tfxEmitterStateFlags_single_shot_done || parent_state_flags & tfxEmitterStateFlags_stop_spawning)
			return 0;
		if (spawn_quantity == 0)
			return 0;

		float step_size = 1.f / spawn_quantity;
		float tween = 0;
		if (step_size == qty_step_size || property_flags & tfxEmitterPropertyFlags_single)
			tween = amount_remainder;
		else
			tween = amount_remainder - (qty_step_size - step_size);
		qty_step_size = step_size;
		//bool is_compute = work_entry.e->property_flags & tfxEmitterPropertyFlags_is_bottom_emitter && pm.state_flags & tfxEffectManagerFlags_use_compute_shader;

		if (tween >= 1) {
			tween -= spawn_quantity;
		}

		work_entry.pm = &pm;
		work_entry.tween = tween;
		work_entry.max_spawn_count = max_spawn_count;
		work_entry.qty_step_size = step_size;
		work_entry.amount_to_spawn = 0;
		work_entry.particle_data = &pm.particle_arrays[pm.emitters.particles_index[work_entry.emitter_index]];

		if (tween >= 1) {
			amount_remainder = tween - 1.f;
		}
		else {
			float amount_that_will_spawn = (1.f - tween) / step_size;
			work_entry.amount_to_spawn = (tfxU32)std::ceilf(amount_that_will_spawn);
			if (work_entry.amount_to_spawn > max_spawn_count) {
				work_entry.amount_to_spawn = max_spawn_count;
				amount_remainder = 0.f;
			}
			else {
				amount_remainder = amount_that_will_spawn - (tfxU32)amount_that_will_spawn;
				amount_remainder = (1.f - amount_remainder) * step_size;
			}
		}

		work_entry.spawn_start_index = AddRows(&pm.particle_array_buffers[pm.emitters.particles_index[work_entry.emitter_index]], work_entry.amount_to_spawn, true);
		tfxEmissionType &emission_type = properties.emission_type[pm.emitters.properties_index[work_entry.emitter_index]];

#if tfxMULTITHREADED
		if (work_entry.amount_to_spawn > 0) {
			work_entry.end_index = work_entry.amount_to_spawn;
			//tfxBumpCompletionCount(&pm.work_queue);
			if (emission_type == tfxPoint) {
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticlePoint2d);
			}
			else if (emission_type == tfxArea) {
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleArea2d);
			}
			else if (emission_type == tfxEllipse) {
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleEllipse2d);
			}
			else if (emission_type == tfxLine) {
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleLine2d);
			}
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleWeight);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleVelocity);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleRoll);
			//Can maybe revisit this. We have to complete the above work before doing the micro update. I would like to add the micro update from one of the above threads
			//when all 4 have finished but synchronisation is hard to get right. Would have to rethink for a multi producer work queue. For now though this is working
			//fine and is stable
			tfxCompleteAllWork(&pm.work_queue);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleMicroUpdate2d);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleAge);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleNoise);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleImageFrame);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleSize2d);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleSpin2d);
		}
#else
		SpawnParticleAge(&pm.work_queue, &work_entry);

		if (emission_type == tfxPoint) {
			SpawnParticlePoint2d(&pm.work_queue, &work_entry);
		}
		else if (emission_type == tfxArea) {
			SpawnParticleArea2d(&pm.work_queue, &work_entry);
		}
		else if (emission_type == tfxEllipse) {
			SpawnParticleEllipse2d(&pm.work_queue, &work_entry);
		}
		else if (emission_type == tfxLine) {
			SpawnParticleLine2d(&pm.work_queue, &work_entry);
		}
		SpawnParticleWeight(&pm.work_queue, &work_entry);
		SpawnParticleVelocity(&pm.work_queue, &work_entry);
		SpawnParticleRoll(&pm.work_queue, &work_entry);
		SpawnParticleMicroUpdate2d(&pm.work_queue, &work_entry);
		SpawnParticleNoise(&pm.work_queue, &work_entry);
		SpawnParticleImageFrame(&pm.work_queue, &work_entry);
		SpawnParticleSize2d(&pm.work_queue, &work_entry);
		SpawnParticleSpin2d(&pm.work_queue, &work_entry);
#endif

		if (work_entry.amount_to_spawn > 0 && property_flags & tfxEmitterPropertyFlags_single)
			pm.emitters.state_flags[work_entry.emitter_index] |= tfxEmitterStateFlags_single_shot_done;

		return work_entry.amount_to_spawn;
	}

	tfxU32 SpawnWideParticles3d(tfxParticleManager &pm, tfxSpawnWorkEntry &work_entry, tfxU32 max_spawn_count) {
		const tfxEmitterPropertiesSoA &properties = *work_entry.properties;
		const tfxEmitterStateFlags state_flags = pm.emitters.state_flags[work_entry.emitter_index];
		const tfxEmitterStateFlags parent_state_flags = pm.emitters.state_flags[work_entry.emitter_index];
		const float spawn_quantity = pm.emitters.spawn_quantity[work_entry.emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[work_entry.emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[work_entry.emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[work_entry.emitter_index];
		const tfxU32 layer = properties.layer[property_index];
		float &qty_step_size = pm.emitters.qty_step_size[work_entry.emitter_index];
		float &amount_remainder = pm.emitters.amount_remainder[work_entry.emitter_index];

		if (state_flags & tfxEmitterStateFlags_single_shot_done || parent_state_flags & tfxEmitterStateFlags_stop_spawning)
			return 0;
		if (spawn_quantity == 0)
			return 0;

		float step_size = 1.f / spawn_quantity;
		float tween = 0;
		if (step_size == qty_step_size || property_flags & tfxEmitterPropertyFlags_single)
			tween = amount_remainder;
		else
			tween = amount_remainder - (qty_step_size - step_size);
		qty_step_size = step_size;
		//bool is_compute = work_entry.e->property_flags & tfxEmitterPropertyFlags_is_bottom_emitter && pm.state_flags & tfxEffectManagerFlags_use_compute_shader;

		if (tween >= 1) {
			tween -= spawn_quantity;
		}

		work_entry.pm = &pm;
		work_entry.tween = tween;
		work_entry.max_spawn_count = max_spawn_count;
		work_entry.qty_step_size = step_size;
		work_entry.amount_to_spawn = 0;
		work_entry.particle_data = &pm.particle_arrays[particles_index];
		if (pm.flags & tfxEffectManagerFlags_order_by_depth)
			pm.new_particles_index_start[layer] = tfxMin(pm.new_particles_index_start[layer], pm.particle_array_buffers[particles_index].current_size);

		if (tween >= 1) {
			amount_remainder = tween - 1.f;
		}
		else {
			float amount_that_will_spawn = (1.f - tween) / step_size;
			work_entry.amount_to_spawn = (tfxU32)std::ceilf(amount_that_will_spawn);
			if (work_entry.amount_to_spawn > max_spawn_count) {
				work_entry.amount_to_spawn = max_spawn_count;
				amount_remainder = 0.f;
			}
			else {
				amount_remainder = amount_that_will_spawn - (tfxU32)amount_that_will_spawn;
				amount_remainder = (1.f - amount_remainder) * step_size;
			}
		}

		work_entry.spawn_start_index = AddRows(&pm.particle_array_buffers[pm.emitters.particles_index[work_entry.emitter_index]], work_entry.amount_to_spawn, true);
		tfxEmissionType &emission_type = properties.emission_type[property_index];

#if tfxMULTITHREADED
		if (work_entry.amount_to_spawn > 0) {
			if (emission_type == tfxPoint) {
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticlePoint3d);
			}
			else if (emission_type == tfxArea) {
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleArea3d);
			}
			else if (emission_type == tfxEllipse) {
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleEllipse3d);
			}
			else if (emission_type == tfxLine) {
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleLine3d);
			}
			else if (emission_type == tfxCylinder) {
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleCylinder3d);
			}
			else if (emission_type == tfxIcosphere && property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleIcosphereRandom3d);
			}
			else if (emission_type == tfxIcosphere && !(property_flags & tfxEmitterPropertyFlags_grid_spawn_random)) {
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleIcosphere3d);
			}
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleWeight);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleVelocity);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleRoll);
			//Can maybe revisit this. We have to complete the above work before doing the micro update. I would like to add the micro update from one of the above threads
			//when all 4 have finished but synchronisation is hard to get right. Would have to rethink for a multi producer work queue. For now though this is working
			//fine and is stable
			tfxCompleteAllWork(&pm.work_queue);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleMicroUpdate3d);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleAge);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleNoise);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleImageFrame);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleSize3d);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleSpin3d);
		}
#else
		SpawnParticleAge(&pm.work_queue, &work_entry);
		if (emission_type == tfxPoint) {
			SpawnParticlePoint3d(&pm.work_queue, &work_entry);
		}
		else if (emission_type == tfxArea) {
			SpawnParticleArea3d(&pm.work_queue, &work_entry);
		}
		else if (emission_type == tfxEllipse) {
			SpawnParticleEllipse3d(&pm.work_queue, &work_entry);
		}
		else if (emission_type == tfxLine) {
			SpawnParticleLine3d(&pm.work_queue, &work_entry);
		}
		else if (emission_type == tfxIcosphere && property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
			SpawnParticleIcosphereRandom3d(&pm.work_queue, &work_entry);
		}
		else if (emission_type == tfxIcosphere && !(property_flags & tfxEmitterPropertyFlags_grid_spawn_random)) {
			SpawnParticleIcosphere3d(&pm.work_queue, &work_entry);
		}
		else if (emission_type == tfxCylinder) {
			SpawnParticleCylinder3d(&pm.work_queue, &work_entry);
		}
		SpawnParticleWeight(&pm.work_queue, &work_entry);
		SpawnParticleVelocity(&pm.work_queue, &work_entry);
		SpawnParticleRoll(&pm.work_queue, &work_entry);
		SpawnParticleMicroUpdate3d(&pm.work_queue, &work_entry);
		SpawnParticleNoise(&pm.work_queue, &work_entry);
		SpawnParticleImageFrame(&pm.work_queue, &work_entry);
		SpawnParticleSize3d(&pm.work_queue, &work_entry);
		SpawnParticleSpin3d(&pm.work_queue, &work_entry);
#endif

		if (work_entry.amount_to_spawn > 0 && property_flags & tfxEmitterPropertyFlags_single)
			pm.emitters.state_flags[work_entry.emitter_index] |= tfxEmitterStateFlags_single_shot_done;

		return work_entry.amount_to_spawn;
	}

	void SpawnParticleAge(tfxWorkQueue *queue, void *data) {

		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;

		const float life = pm.emitters.life[emitter_index];
		const float life_variation = pm.emitters.life_variation[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxU32 sprites_index = pm.emitters.sprites_index[emitter_index];
		const float highest_particle_age = pm.emitters.highest_particle_age[emitter_index];
		const tfxU32 loop_count = entry->properties->single_shot_limit[property_index] + 1;
		tfxEffectLibrary *library = pm.emitters.library[emitter_index];
		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[emitter_index];
		const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];
		const float emitter_intensity = pm.emitters.intensity[emitter_index];
		const float first_red_value = library->emitter_attributes[emitter_attributes].overtime.red.GetFirstValue();
		const float first_green_value = library->emitter_attributes[emitter_attributes].overtime.green.GetFirstValue();
		const float first_blue_value = library->emitter_attributes[emitter_attributes].overtime.blue.GetFirstValue();
		const float first_intensity_value = library->emitter_attributes[emitter_attributes].overtime.intensity.GetFirstValue();

		for(int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &age = entry->particle_data->age[index];
			float &max_age = entry->particle_data->max_age[index];
			tfxRGBA8 &color = entry->particle_data->color[index];
			float &intensity = entry->particle_data->intensity[index];
			tfxU32 &single_loop_count = entry->particle_data->single_loop_count[index];

			tfxParticleFlags &flags = entry->particle_data->flags[index];
			tfxU32 &parent = entry->particle_data->parent_index[index];
			tfxU32 &particle_index = entry->particle_data->particle_index[index];
			parent = emitter_index;
			particle_index = tfxINVALID;

			flags = 0;

			//Max age
			//Todo: should age be set to the tween value?
			age = 0.f;
			max_age = life + random_generation.Range(life_variation);
			single_loop_count = 0;

			color.a = unsigned char(255.f * library->emitter_attributes[emitter_attributes].overtime.blendfactor.GetFirstValue());
			intensity = first_intensity_value * emitter_intensity;
			//intensity = 0.f;
			if (emitter_flags & tfxEmitterStateFlags_random_color) {
				float age = random_generation.Range(max_age);
				color.r = unsigned char(255.f * lookup_overtime_callback(library->emitter_attributes[emitter_attributes].overtime.red, age, max_age));
				color.g = unsigned char(255.f * lookup_overtime_callback(library->emitter_attributes[emitter_attributes].overtime.green, age, max_age));
				color.b = unsigned char(255.f * lookup_overtime_callback(library->emitter_attributes[emitter_attributes].overtime.blue, age, max_age));
			}
			else {
				color.r = unsigned char(255.f * first_red_value);
				color.g = unsigned char(255.f * first_green_value);
				color.b = unsigned char(255.f * first_blue_value);
			}

			entry->highest_particle_age = std::fmaxf(highest_particle_age, (max_age * loop_count) + tfxFRAME_LENGTH + 1);

			if (entry->sub_effects->current_size > 0) {
				particle_index = pm.GetParticleIndexSlot(MakeParticleID(particles_index, index));
				flags |= tfxParticleFlags_has_sub_effects;
				for (auto &sub : *entry->sub_effects) {
					if (!pm.FreeEffectCapacity())
						break;
					assert(entry->depth < tfxMAXDEPTH - 1);
					tfxU32 added_index = pm.AddEffect(sub, pm.current_ebuff, entry->depth + 1, true);
					pm.effects.overal_scale[added_index] = pm.emitters.overal_scale[index];
					pm.effects.parent_particle_index[added_index] = particle_index;
				}
			}
		}

	}

	void SpawnParticleImageFrame(tfxWorkQueue *queue, void *data) {

		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		const tfxU32 index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const tfxU32 particles_index = pm.emitters.particles_index[index];
		const tfxU32 property_index = pm.emitters.properties_index[index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[index];

		for(int i = 0; i != entry->amount_to_spawn; ++i) {

			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &image_frame = entry->particle_data->image_frame[index];

			//----Image
			//data.image = GetProperties().image;
			if (property_flags & tfxEmitterPropertyFlags_random_start_frame && properties.image[property_index]->animation_frames > 1) {
				image_frame = random_generation.Range(properties.image[property_index]->animation_frames);
			}
			else {
				image_frame = properties.start_frame[property_index];
			}

		}

	}

	void SpawnParticleSize2d(tfxWorkQueue *queue, void *data) {

		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		const tfxEmitterPropertiesSoA &properties = *entry->properties;

		const tfxVec2 size = pm.emitters.size[index];
		const tfxVec2 size_variation = pm.emitters.size_variation[index];
		const tfxU32 particles_index = pm.emitters.particles_index[index];
		const tfxU32 property_index = pm.emitters.properties_index[index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[index];
		tfxVec2 &image_size = properties.image[property_index]->image_size;

		for(int i = 0; i != entry->amount_to_spawn; ++i) {

			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			tfxVec2 &base_size = entry->particle_data->base_size[index];

			//----Size
			if (!(property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
				float random_size_x = random_generation.Range(size_variation.x);
				float random_size_y = random_generation.Range(size_variation.y);
				base_size.y = random_size_y + size.y;
				base_size.x = (random_size_x + size.x) / image_size.x;
			}
			else {
				float random_size_x = random_generation.Range(size_variation.x);
				float random_size_y = random_size_x;
				base_size.y = random_size_y + size.y;
				base_size.x = (random_size_x + size.x) / image_size.x;
			}
		}

	}

	void SpawnParticleSize3d(tfxWorkQueue *queue, void *data) {

		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		const tfxEmitterPropertiesSoA &properties = *entry->properties;

		const tfxVec2 size = pm.emitters.size[index];
		const tfxVec2 size_variation = pm.emitters.size_variation[index];
		const tfxU32 particles_index = pm.emitters.particles_index[index];
		const tfxU32 property_index = pm.emitters.properties_index[index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[index];
		tfxVec2 &image_size = properties.image[property_index]->image_size;

		for (int i = 0; i != entry->amount_to_spawn; ++i) {

			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			tfxVec2 &base_size = entry->particle_data->base_size[index];

			//----Size
			if (!(property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
				float random_size_x = random_generation.Range(size_variation.x);
				float random_size_y = random_generation.Range(size_variation.y);
				base_size.y = random_size_y + size.y;
				base_size.x = random_size_x + size.x;
			}
			else {
				float random_size_x = random_generation.Range(size_variation.x);
				float random_size_y = random_size_x;
				base_size.y = random_size_y + size.y;
				base_size.x = random_size_x + size.x;
			}

		}

	}

	void SpawnParticleNoise(tfxWorkQueue *queue, void *data) {

		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		const float emitter_noise_offset_variation = pm.emitters.noise_offset_variation[index];
		const float emitter_noise_offset = pm.emitters.noise_offset[index];
		const float emitter_noise_resolution = pm.emitters.noise_resolution[index];
		const tfxU32 particles_index = pm.emitters.particles_index[index];

		for(int i = 0; i != entry->amount_to_spawn; ++i) {

			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &noise_offset = entry->particle_data->noise_offset[index];
			float &noise_resolution = entry->particle_data->noise_resolution[index];

			//----Motion randomness
			noise_offset = random_generation.Range(emitter_noise_offset_variation) + emitter_noise_offset;
			noise_resolution = emitter_noise_resolution + 0.01f;

		}
	}

	void SpawnParticleSpin2d(tfxWorkQueue *queue, void *data) {

		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;

		const float spin_variation = pm.emitters.spin_variation[index];
		const float spin = pm.emitters.spin[index];
		const tfxU32 particles_index = pm.emitters.particles_index[index];

		for(int i = 0; i != entry->amount_to_spawn; ++i) {

			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &base_spin = entry->particle_data->base_spin[index];

			//----Spin
			base_spin = random_generation.Range(-spin_variation, spin_variation) + spin;

		}

	}

	void SpawnParticleSpin3d(tfxWorkQueue *queue, void *data) {

		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;

		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const float spin_variation = pm.emitters.spin_variation[emitter_index];
		const float spin = pm.emitters.spin[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxVec3 angle_offsets = pm.emitters.angle_offsets[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxAngleSettingFlags angle_settings = entry->properties->angle_settings[property_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {

			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &base_spin = entry->particle_data->base_spin[index];
			tfxVec3 &local_rotations = entry->particle_data->local_rotations[index];

			//----Spin
			base_spin = random_generation.Range(-spin_variation, spin_variation) + spin;
			if (angle_settings & tfxAngleSettingFlags_specify_roll)
				local_rotations.roll = angle_offsets.roll;
			if (angle_settings & tfxAngleSettingFlags_specify_pitch)
				local_rotations.pitch = angle_offsets.pitch;
			if (angle_settings & tfxAngleSettingFlags_specify_yaw)
				local_rotations.yaw = angle_offsets.yaw;
			if (angle_settings & tfxAngleSettingFlags_random_pitch)
				local_rotations.pitch = random_generation.Range(angle_offsets.pitch);
			if (angle_settings & tfxAngleSettingFlags_random_yaw)
				local_rotations.yaw = random_generation.Range(angle_offsets.yaw);
			if (angle_settings & tfxAngleSettingFlags_random_roll)
				local_rotations.roll = random_generation.Range(angle_offsets.roll);

		}

	}

	void SpawnParticlePoint2d(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		
		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			tfxVec3 &local_position = entry->particle_data->local_position[index];

			local_position = 0;
			tfxVec2 lerp_position = InterpolateVec2(tween, emitter_captured_position.xy(), emitter_world_position.xy());
			if (property_flags & tfxEmitterPropertyFlags_relative_position)
				local_position = 0;
			else {
				if (property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) {
					local_position = lerp_position;
				}
				else {
					tfxVec2 rotvec = mmTransformVector(matrix, -handle.xy());
					local_position = rotvec + lerp_position;
				}
			}
			tween += entry->qty_step_size;
		}

	}

	void SpawnParticlePoint3d(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			tfxVec3 &local_position = entry->particle_data->local_position[index];

			local_position = 0;
			tfxVec3 lerp_position = InterpolateVec3(tween, emitter_captured_position, emitter_world_position);
			if (property_flags & tfxEmitterPropertyFlags_relative_position)
				local_position = 0;
			else {
				if (property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) {
					local_position = lerp_position;
				}
				else {
					tfxVec3 rotvec = mmTransformVector3(matrix, -handle);
					local_position = rotvec + lerp_position;
				}
			}
			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleLine2d(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxVec3 &grid_points = properties.grid_points[property_index];
		const tfxVec3 &grid_segment_size = pm.emitters.grid_segment_size[emitter_index];
		const tfxVec2 emitter_size = pm.emitters.emitter_size[emitter_index].xy();
		const tfxVec3 scale = pm.emitters.scale[emitter_index];
		tfxVec3 &grid_coords = pm.emitters.grid_coords[emitter_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			tfxVec3 &local_position = entry->particle_data->local_position[index];

			local_position = 0;
			tfxVec2 lerp_position = InterpolateVec2(tween, emitter_captured_position.xy(), emitter_world_position.xy());

			if (property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				grid_coords.x = 0.f;

				if (!(property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					grid_coords.y--;
					if (grid_coords.y < 0.f) {
						grid_coords.y = grid_points.x - 1;
					}
				}

				local_position = tfxVec2(grid_coords.xy() * -grid_segment_size.xy());

				if (property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					grid_coords.y++;
					if (grid_coords.y >= grid_points.x) {
						grid_coords.y = 0.f;
					}
				}

			}
			else {
				local_position.x = 0.f;
				local_position.y = random_generation.Range(-emitter_size.y, 0.f);

			}

			//----TForm and Emission
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position) && !(property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
				local_position = mmTransformVector(matrix, local_position.xy() + handle.xy());
				local_position = lerp_position + local_position.xy() * scale.xy();
			}

			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleLine3d(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxVec3 &grid_points = properties.grid_points[property_index];
		const tfxVec3 &grid_segment_size = pm.emitters.grid_segment_size[emitter_index];
		const tfxVec2 emitter_size = pm.emitters.emitter_size[emitter_index].xy();
		const tfxVec3 scale = pm.emitters.scale[emitter_index];
		tfxVec3 &grid_coords = pm.emitters.grid_coords[emitter_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			tfxVec3 &local_position = entry->particle_data->local_position[index];

			local_position = 0;
			tfxVec3 lerp_position = InterpolateVec3(tween, emitter_captured_position, emitter_world_position);

			if (property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				grid_coords.x = 0.f;
				grid_coords.z = 0.f;

				if (property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
					grid_coords.y = (float)random_generation.RangeUInt((tfxU32)grid_points.x);
					local_position = grid_coords * grid_segment_size;
				}
				else {

					if (!(property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
						grid_coords.y--;
						if (grid_coords.y < 0.f) {
							grid_coords.y = grid_points.x - 1;
						}
					}

					local_position = grid_coords * grid_segment_size;

					if (property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
						grid_coords.y++;
						if (grid_coords.y >= grid_points.x) {
							grid_coords.y = 0.f;
						}
					}
				}

			}
			else {
				local_position.x = 0.f;
				local_position.y = random_generation.Range(0.f, emitter_size.y);
				local_position.z = 0.f;
			}

			//----TForm and Emission
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position) && !(property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
				local_position = mmTransformVector3(matrix, local_position + handle);
				local_position = lerp_position + local_position * scale;
			}
	
			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleArea2d(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxVec3 &grid_points = properties.grid_points[property_index];
		const tfxVec3 &grid_segment_size = pm.emitters.grid_segment_size[emitter_index];
		const tfxVec2 emitter_size = pm.emitters.emitter_size[emitter_index].xy();
		const tfxVec3 scale = pm.emitters.scale[emitter_index];
		tfxVec3 &grid_direction = pm.emitters.grid_direction[emitter_index];
		tfxVec3 &grid_coords = pm.emitters.grid_coords[emitter_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			tfxVec3 &local_position = entry->particle_data->local_position[index];

			local_position = 0;
			tfxVec2 lerp_position = InterpolateVec2(tween, emitter_captured_position.xy(), emitter_world_position.xy());

			tfxVec2 position = tfxVec2(0.f, 0.f);

			if (property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				if (property_flags & tfxEmitterPropertyFlags_fill_area) {
					if (property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
						grid_coords.x = (float)random_generation.RangeUInt((tfxU32)grid_points.x);
						grid_coords.y = (float)random_generation.RangeUInt((tfxU32)grid_points.y);
						local_position = grid_coords.xy() * grid_segment_size.xy();
					}
					else {
						if (!(property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
							grid_coords.x--;
							if (grid_coords.x < 0.f) {
								grid_coords.y--;
								grid_coords.x = grid_points.x - 1;
								if (grid_coords.y < 0.f)
									grid_coords.y = grid_points.y - 1;
							}
						}

						local_position = grid_coords.xy() * grid_segment_size.xy();

						if (property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
							grid_coords.x++;
							if (grid_coords.x == grid_points.x) {
								grid_coords.y++;
								grid_coords.x = 0.f;
								if (grid_coords.y >= grid_points.y)
									grid_coords.y = 0.f;
							}
						}
					}
				}
				else {
					if (property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
						tfxU32 side = random_generation.RangeUInt(4);
						if (side == 0) {
							//left side
							grid_coords.x = 0.f;
							grid_coords.y = (float)random_generation.RangeUInt((tfxU32)grid_points.y);
						}
						else if (side == 1) {
							//right side
							grid_coords.x = grid_points.x - 1;
							grid_coords.y = (float)random_generation.RangeUInt((tfxU32)grid_points.y);
						}
						else if (side == 2) {
							//top side
							grid_coords.x = (float)random_generation.RangeUInt((tfxU32)grid_points.x);
							grid_coords.y = 0.f;
						}
						else if (side == 3) {
							//bottom side
							grid_coords.x = (float)random_generation.RangeUInt((tfxU32)grid_points.x);
							grid_coords.y = grid_points.y - 1;
						}
						local_position = grid_coords.xy() * grid_segment_size.xy();
					}
					else {
						if (property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {

							grid_direction.x = 1;
							grid_direction.y = 0;
							if (grid_coords.x == grid_points.x - 1 && grid_coords.y >= 0 && grid_coords.y < grid_points.y - 1) {
								grid_direction.x = 0;
								grid_direction.y = 1;
							}
							else if (grid_coords.x > 0 && grid_coords.x < grid_points.x && grid_coords.y == grid_points.y - 1) {
								grid_direction.x = -1;
								grid_direction.y = 0;
							}
							else if (grid_coords.x == 0 && grid_coords.y > 0 && grid_coords.y < grid_points.y) {
								grid_direction.x = 0;
								grid_direction.y = -1;
							}

						}
						else {

							grid_direction.x = -1;
							grid_direction.y = 0;
							if (grid_coords.x == grid_points.x - 1 && grid_coords.y > 0 && grid_coords.y < grid_points.y) {
								grid_direction.x = 0;
								grid_direction.y = -1;
							}
							else if (grid_coords.x >= 0 && grid_coords.x < grid_points.x - 1 && grid_coords.y == grid_points.y - 1) {
								grid_direction.x = 1;
								grid_direction.y = 0;
							}
							else if (grid_coords.x == 0 && grid_coords.y >= 0 && grid_coords.y < grid_points.y - 1) {
								grid_direction.x = 0;
								grid_direction.y = 1;
							}

						}

						grid_coords += grid_direction;
						tfxBound(grid_coords.xy(), grid_points.xy());
						local_position = position + (grid_coords.xy() * grid_segment_size.xy());
					}
				}
			}
			else {
				if (property_flags & tfxEmitterPropertyFlags_fill_area) {
					position.x = random_generation.Range(emitter_size.x);
					position.y = random_generation.Range(emitter_size.y);
				}
				else {
					//Spawn on one of 4 edges of the area
					tfxU32 side = random_generation.RangeUInt(4);
					if (side == 0) {
						//left side
						position.x = 0.f;
						position.y = random_generation.Range(emitter_size.y);
					}
					else if (side == 1) {
						//right side
						position.x = emitter_size.x;
						position.y = random_generation.Range(emitter_size.y);
					}
					else if (side == 2) {
						//top side
						position.x = random_generation.Range(emitter_size.x);
						position.y = 0.f;
					}
					else if (side == 3) {
						//bottom side
						position.x = random_generation.Range(emitter_size.x);
						position.y = emitter_size.y;
					}
				}

				local_position = position;
			}

			//----TForm and Emission
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				local_position = mmTransformVector(matrix, local_position.xy() + handle.xy());
				local_position = lerp_position + local_position.xy() * scale.xy();
			}

			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleArea3d(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxVec3 &grid_points = properties.grid_points[property_index];
		const tfxVec3 &grid_segment_size = pm.emitters.grid_segment_size[emitter_index];
		const tfxVec3 emitter_size = pm.emitters.emitter_size[emitter_index];
		const tfxVec3 scale = pm.emitters.scale[emitter_index];
		tfxVec3 &grid_direction = pm.emitters.grid_direction[emitter_index];
		tfxVec3 &grid_coords = pm.emitters.grid_coords[emitter_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			tfxVec3 &local_position = entry->particle_data->local_position[index];

			local_position = 0;
			tfxVec3 lerp_position = InterpolateVec3(tween, emitter_captured_position, emitter_world_position);

			tfxVec3 position;

			if (property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				if (property_flags & tfxEmitterPropertyFlags_fill_area) {

					if (property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
						grid_coords.x = (float)random_generation.RangeUInt((tfxU32)grid_points.x);
						grid_coords.y = (float)random_generation.RangeUInt((tfxU32)grid_points.y);
						grid_coords.z = (float)random_generation.RangeUInt((tfxU32)grid_points.z);

						local_position = grid_coords * grid_segment_size;
					}
					else {
						if (!(property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
							grid_coords.x--;
							if (grid_coords.x < 0.f) {
								grid_coords.y--;
								grid_coords.x = grid_points.x - 1;
								if (grid_coords.y < 0.f) {
									grid_coords.z--;
									grid_coords.y = grid_points.y - 1;
									if (grid_coords.z < 0.f)
										grid_coords.z = grid_points.z;
								}
							}
						}

						local_position = position + (grid_coords * grid_segment_size);

						if (property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
							grid_coords.x++;
							if (grid_coords.x == grid_points.x) {
								grid_coords.y++;
								grid_coords.x = 0.f;
								if (grid_coords.y >= grid_points.y) {
									grid_coords.z++;
									grid_coords.y = 0.f;
									if (grid_coords.z >= grid_points.z)
										grid_coords.z = 0.f;
								}
							}
						}
					}

				}
				else {
					if (property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
						//Spawn on one of 6 edges of the cuboid
						tfxU32 side = random_generation.RangeUInt((property_flags & tfxEmitterPropertyFlags_area_open_ends) ? 4 : 6);
						if (side == 0) {
							//left side
							grid_coords.x = 0.f;
							grid_coords.y = (float)random_generation.RangeUInt((tfxU32)grid_points.y);
							grid_coords.z = (float)random_generation.RangeUInt((tfxU32)grid_points.z);
						}
						else if (side == 1) {
							//right side
							grid_coords.x = grid_points.x - 1;
							grid_coords.y = (float)random_generation.RangeUInt((tfxU32)grid_points.y);
							grid_coords.z = (float)random_generation.RangeUInt((tfxU32)grid_points.z);
						}
						else if (side == 2) {
							//top side
							grid_coords.x = (float)random_generation.RangeUInt((tfxU32)grid_points.x);
							grid_coords.y = 0.f;
							grid_coords.z = (float)random_generation.RangeUInt((tfxU32)grid_points.z);
						}
						else if (side == 3) {
							//bottom side
							grid_coords.x = (float)random_generation.RangeUInt((tfxU32)grid_points.x);
							grid_coords.y = grid_points.y - 1;
							grid_coords.z = (float)random_generation.RangeUInt((tfxU32)grid_points.z);
						}
						else if (side == 4) {
							//End far
							grid_coords.x = (float)random_generation.RangeUInt((tfxU32)grid_points.x);
							grid_coords.y = (float)random_generation.RangeUInt((tfxU32)grid_points.y);
							grid_coords.z = grid_points.z - 1;
						}
						else if (side == 5) {
							//End near
							grid_coords.x = (float)random_generation.RangeUInt((tfxU32)grid_points.x);
							grid_coords.y = (float)random_generation.RangeUInt((tfxU32)grid_points.y);
							grid_coords.z = 0.f;
						}
						local_position = grid_coords * grid_segment_size;
					}
					else {
						if (property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
							if (grid_direction.z == 0) {
								//right side
								grid_coords.z--;
								grid_coords.x = 0.f;
								if (grid_coords.z < 0.f) {
									grid_coords.y++;
									grid_coords.z = grid_points.z - 1;
									if (grid_coords.y >= grid_points.y - 1) {
										grid_coords.y = grid_points.y - 1;
										grid_direction.z = 2;
									}
								}
							}
							else if (grid_direction.z == 1) {
								//left side
								grid_coords.z--;
								grid_coords.x = grid_points.x - 1;
								if (grid_coords.z < 0.f) {
									grid_coords.y--;
									grid_coords.z = grid_points.z - 1;
									if (grid_coords.y < 0) {
										grid_coords.y = 0;
										grid_coords.x--;
										grid_direction.z = 3;
									}
								}
							}
							else if (grid_direction.z == 2) {
								//top side
								grid_coords.z--;
								grid_coords.y = grid_points.y - 1;
								if (grid_coords.z < 0.f) {
									grid_coords.x++;
									grid_coords.z = grid_points.z - 1;
									if (grid_coords.x >= grid_points.x - 1) {
										grid_coords.x = grid_points.x - 1;
										grid_direction.z = 1;
									}
								}
							}
							else if (grid_direction.z == 3) {
								//bottom side
								grid_coords.z--;
								grid_coords.y = 0.f;
								if (grid_coords.z < 0.f) {
									grid_coords.x--;
									grid_coords.z = grid_points.z - 1;
									if (grid_coords.x < 0) {
										grid_coords.x = 0.f;
										grid_coords.y = 1.f;
										grid_direction.z = (property_flags & tfxEmitterPropertyFlags_area_open_ends) ? 0.f : 4.f;
									}
								}
							}
							else if (grid_direction.z == 4) {
								//End far
								grid_coords.x++;
								grid_coords.z = 0.f;
								if (grid_coords.x >= grid_points.x) {
									grid_coords.y++;
									grid_coords.x = 0.f;
									if (grid_coords.y >= grid_points.y - 1) {
										grid_coords.y = grid_points.y - 1;
										grid_coords.x = grid_points.x - 1;
										grid_direction.z = 5;
									}
								}
							}
							else if (grid_direction.z == 5) {
								//End near
								grid_coords.x--;
								grid_coords.z = grid_points.z - 1;
								if (grid_coords.x < 0.f) {
									grid_coords.y--;
									grid_coords.x = grid_points.x - 1;
									if (grid_coords.y < 0.f) {
										grid_coords.y = 0.f;
										grid_direction.z = 0;
									}
								}
							}
						}
						else {
							if (grid_direction.z == 0) {
								//right side
								grid_coords.z--;
								grid_coords.x = 0.f;
								if (grid_coords.z < 0.f) {
									grid_coords.y--;
									grid_coords.z = grid_points.z - 1;
									if (grid_coords.y < 1) {
										grid_coords.y = 0;
										grid_direction.z = 3;
									}
								}
							}
							else if (grid_direction.z == 1) {
								//left side
								grid_coords.z--;
								grid_coords.x = grid_points.x - 1;
								if (grid_coords.z < 0.f) {
									grid_coords.y++;
									grid_coords.z = grid_points.z - 1;
									if (grid_coords.y >= grid_points.y) {
										grid_coords.y = grid_points.y - 1;
										grid_coords.x--;
										grid_direction.z = 2;
									}
								}
							}
							else if (grid_direction.z == 2) {
								//top side
								grid_coords.z--;
								grid_coords.y = grid_points.y - 1;
								if (grid_coords.z < 0.f) {
									grid_coords.x--;
									grid_coords.z = grid_points.z - 1;
									if (grid_coords.x < 1) {
										grid_coords.x = 0.f;
										grid_direction.z = (property_flags & tfxEmitterPropertyFlags_area_open_ends) ? 0.f : 4.f;
									}
								}
							}
							else if (grid_direction.z == 3) {
								//bottom side
								grid_coords.z--;
								grid_coords.y = 0.f;
								if (grid_coords.z < 0.f) {
									grid_coords.x++;
									grid_coords.z = grid_points.z - 1;
									if (grid_coords.x >= grid_points.x - 1) {
										grid_coords.x = grid_points.x - 1;
										grid_coords.y = 0.f;
										grid_direction.z = 1;
									}
								}
							}
							else if (grid_direction.z == 4) {
								//End far
								grid_coords.x++;
								grid_coords.z = 0.f;
								if (grid_coords.x >= grid_points.x) {
									grid_coords.y--;
									grid_coords.x = 0.f;
									if (grid_coords.y < 0) {
										grid_coords.y = 0.f;
										grid_coords.x = grid_points.x - 1;
										grid_direction.z = 5;
									}
								}
							}
							else if (grid_direction.z == 5) {
								//End near
								grid_coords.x--;
								grid_coords.z = grid_points.z - 1;
								if (grid_coords.x < 0.f) {
									grid_coords.y++;
									grid_coords.x = grid_points.x - 1;
									if (grid_coords.y >= grid_points.y - 1) {
										grid_coords.y = grid_points.y - 1;
										grid_direction.z = 0.f;
									}
								}
							}
						}
						tfxBound3d(grid_coords, grid_points);
						local_position = position + (grid_coords * grid_segment_size);
					}

				}
			}
			else {
				if (property_flags & tfxEmitterPropertyFlags_fill_area) {
					position.x = random_generation.Range(emitter_size.x);
					position.y = random_generation.Range(emitter_size.y);
					position.z = random_generation.Range(emitter_size.z);
				}
				else {
					//Spawn on one of 6 edges of the cuboid
					tfxU32 side = random_generation.RangeUInt((property_flags & tfxEmitterPropertyFlags_area_open_ends) ? 4 : 6);
					if (side == 0) {
						//left side
						position.x = 0.f;
						position.y = random_generation.Range(emitter_size.y);
						position.z = random_generation.Range(emitter_size.z);
					}
					else if (side == 1) {
						//right side
						position.x = emitter_size.x;
						position.y = random_generation.Range(emitter_size.y);
						position.z = random_generation.Range(emitter_size.z);
					}
					else if (side == 2) {
						//top side
						position.x = random_generation.Range(emitter_size.x);
						position.y = 0.f;
						position.z = random_generation.Range(emitter_size.z);
					}
					else if (side == 3) {
						//bottom side
						position.x = random_generation.Range(emitter_size.x);
						position.y = emitter_size.y;
						position.z = random_generation.Range(emitter_size.z);
					}
					else if (side == 4) {
						//End far
						position.x = random_generation.Range(emitter_size.x);
						position.y = random_generation.Range(emitter_size.y);
						position.z = emitter_size.z;
					}
					else if (side == 5) {
						//End near
						position.x = random_generation.Range(emitter_size.x);
						position.y = random_generation.Range(emitter_size.y);
						position.z = 0.f;
					}
				}

				local_position = position;
			}

			//----TForm and Emission
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				local_position = mmTransformVector3(matrix, local_position + handle);
				local_position = lerp_position + local_position * scale;
			}
	
			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleEllipse2d(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxVec3 &grid_points = properties.grid_points[property_index];
		const tfxVec3 &grid_segment_size = pm.emitters.grid_segment_size[emitter_index];
		const tfxVec2 emitter_size = pm.emitters.emitter_size[emitter_index].xy();
		const tfxVec3 scale = pm.emitters.scale[emitter_index];
		const float arc_offset = pm.emitters.arc_offset[emitter_index];
		const float arc_size = pm.emitters.arc_size[emitter_index];
		tfxVec3 &grid_coords = pm.emitters.grid_coords[emitter_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			tfxVec3 &local_position = entry->particle_data->local_position[index];

			local_position = 0;
			tfxVec2 lerp_position = InterpolateVec2(tween, emitter_captured_position.xy(), emitter_world_position.xy());

			tfxVec2 half_emitter_size = (emitter_size * .5f);
			tfxVec2 position = tfxVec2(0.f, 0.f);

			if (property_flags & tfxEmitterPropertyFlags_spawn_on_grid && !(property_flags & tfxEmitterPropertyFlags_fill_area)) {

				grid_coords.y = 0.f;

				if (property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					grid_coords.x--;
					if (grid_coords.x < 0.f) {
						grid_coords.x = grid_points.x - 1;
					}
				}

				float th = grid_coords.x * grid_segment_size.x + arc_offset;
				local_position = tfxVec2(std::cosf(th) * half_emitter_size.x + half_emitter_size.x, -std::sinf(th) * half_emitter_size.y + half_emitter_size.y);

				if (!(property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					grid_coords.x++;
					if (grid_coords.x >= grid_points.x) {
						grid_coords.x = 0.f;
					}
				}

			}
			else if (!(property_flags & tfxEmitterPropertyFlags_fill_area)) {
				float th = random_generation.Range(arc_size) + arc_offset;

				local_position = tfxVec2(std::cosf(th) * half_emitter_size.x + half_emitter_size.x, -std::sinf(th) * half_emitter_size.y + half_emitter_size.y);

			}
			else {
				local_position.x = random_generation.Range(0, emitter_size.x);
				local_position.y = random_generation.Range(0, emitter_size.y);

				while ((std::pow(local_position.x - half_emitter_size.x, 2) / std::pow(half_emitter_size.x, 2)) + (std::pow(local_position.y - half_emitter_size.y, 2) / std::pow(half_emitter_size.y, 2)) > 1) {
					local_position.x = random_generation.Range(0, emitter_size.x);
					local_position.y = random_generation.Range(0, emitter_size.y);
				}
			}

			//----TForm and Emission
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				local_position = mmTransformVector(matrix, local_position.xy() + handle.xy());
				local_position = lerp_position + local_position.xy() * scale.xy();
			}


			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleEllipse3d(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxVec3 emitter_size = pm.emitters.emitter_size[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxVec3 &grid_points = properties.grid_points[property_index];
		const tfxVec3 &grid_segment_size = pm.emitters.grid_segment_size[emitter_index];
		const tfxVec3 scale = pm.emitters.scale[emitter_index];
		const float arc_offset = pm.emitters.arc_offset[emitter_index];
		const float arc_size = pm.emitters.arc_size[emitter_index];
		tfxVec3 &grid_coords = pm.emitters.grid_coords[emitter_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			tfxVec3 &local_position = entry->particle_data->local_position[index];

			local_position = 0;
			tfxVec3 lerp_position = InterpolateVec3(tween, emitter_captured_position, emitter_world_position);

			tfxVec3 half_emitter_size = emitter_size * .5f;
			tfxVec3 position;

			if (!(property_flags & tfxEmitterPropertyFlags_fill_area)) {
				float u = random_generation.Range(1.f);
				float v = random_generation.Range(1.f);
				float theta = u * 2.f * tfxPI;
				float phi = std::acosf(2.f * v - 1.f);
				float sin_theta = std::sinf(theta);
				float cos_theta = std::cosf(theta);
				float sin_phi = std::sinf(phi);
				float cos_phi = std::cosf(phi);
				local_position.x = half_emitter_size.x * sin_phi * cos_theta;
				local_position.y = half_emitter_size.y * sin_phi * sin_theta;
				local_position.z = half_emitter_size.z * cos_phi;
			}
			else {
				position.x = random_generation.Range(-half_emitter_size.x, half_emitter_size.x);
				position.y = random_generation.Range(-half_emitter_size.y, half_emitter_size.y);
				position.z = random_generation.Range(-half_emitter_size.z, half_emitter_size.z);

				while (std::powf(position.x / half_emitter_size.x, 2.f) + std::powf(position.y / half_emitter_size.y, 2.f) + std::powf(position.z / half_emitter_size.z, 2.f) > 1.f) {
					position.x = random_generation.Range(-half_emitter_size.x, half_emitter_size.x);
					position.y = random_generation.Range(-half_emitter_size.y, half_emitter_size.y);
					position.z = random_generation.Range(-half_emitter_size.z, half_emitter_size.z);
				}

				local_position = position;
			}

			//----TForm and Emission
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				local_position = mmTransformVector3(matrix, local_position + handle);
				local_position = lerp_position + local_position * scale;
			}

			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleIcosphere3d(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxVec3 &grid_points = properties.grid_points[property_index];
		const tfxVec3 emitter_size = pm.emitters.emitter_size[emitter_index];
		const tfxVec3 scale = pm.emitters.scale[emitter_index];
		tfxVec3 &grid_coords = pm.emitters.grid_coords[emitter_index];
		tfxVec3 half_emitter_size = emitter_size * .5f;

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			tfxVec3 &local_position = entry->particle_data->local_position[index];

			local_position = 0;
			tfxVec3 lerp_position = InterpolateVec3(tween, emitter_captured_position, emitter_world_position);

			tfxU32 sub_division = (tfxU32)grid_points.x;
			local_position = tfxIcospherePoints[sub_division][(tfxU32)grid_coords.x] * half_emitter_size;
			if (++grid_coords.x >= tfxIcospherePoints[sub_division].current_size) {
				grid_coords.x = 0;
			}
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				local_position = mmTransformVector3(matrix, local_position + handle);
				local_position = lerp_position + local_position * scale;
			}

			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleIcosphereRandom3d(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxVec3 &grid_points = properties.grid_points[property_index];
		const tfxVec3 emitter_size = pm.emitters.emitter_size[emitter_index];
		const tfxVec3 scale = pm.emitters.scale[emitter_index];
		tfxVec3 &grid_coords = pm.emitters.grid_coords[emitter_index];
		tfxVec3 half_emitter_size = emitter_size * .5f;

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			tfxVec3 &local_position = entry->particle_data->local_position[index];

			local_position = 0;
			tfxVec3 lerp_position = InterpolateVec3(tween, emitter_captured_position, emitter_world_position);

			tfxVec3 emitter_size = emitter_size * .5f;
			tfxU32 sub_division = (tfxU32)grid_points.x;
			int ico_point = random_generation.RangeUInt(tfxIcospherePoints[sub_division].current_size);
			local_position = tfxIcospherePoints[sub_division][ico_point] * emitter_size;
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				local_position = mmTransformVector3(matrix, local_position + handle);
				local_position = lerp_position + local_position * scale;
			}

			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleCylinder3d(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxVec3 &grid_points = properties.grid_points[property_index];
		const tfxVec3 &grid_segment_size = pm.emitters.grid_segment_size[emitter_index];
		const tfxVec3 emitter_size = pm.emitters.emitter_size[emitter_index];
		const tfxVec3 scale = pm.emitters.scale[emitter_index];
		const float arc_offset = pm.emitters.arc_offset[emitter_index];
		const float arc_size = pm.emitters.arc_size[emitter_index];
		tfxVec3 &grid_coords = pm.emitters.grid_coords[emitter_index];
		tfxVec3 half_emitter_size = emitter_size * .5f;

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			tfxVec3 &local_position = entry->particle_data->local_position[index];

			local_position = 0;
			tfxVec3 lerp_position = InterpolateVec3(tween, emitter_captured_position, emitter_world_position);

			if (property_flags & tfxEmitterPropertyFlags_spawn_on_grid && !(property_flags & tfxEmitterPropertyFlags_fill_area)) {

				grid_coords.z = 0.f;

				if (property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
					grid_coords.x = (float)random_generation.RangeUInt((tfxU32)grid_points.x);
					grid_coords.y = (float)random_generation.RangeUInt((tfxU32)grid_points.y);

					float th = grid_coords.x * grid_segment_size.x + arc_offset;
					local_position = tfxVec3(std::cosf(th) * half_emitter_size.x + half_emitter_size.x, grid_coords.y * grid_segment_size.y, -std::sinf(th) * half_emitter_size.z + half_emitter_size.z);
				}
				else {
					if (property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
						grid_coords.x--;
						if (grid_coords.x < 0.f) {
							grid_coords.x = grid_points.x - 1;
							grid_coords.y--;
							if (grid_coords.y < 0.f) {
								grid_coords.y = grid_points.y - 1;
							}
						}
					}

					float th = grid_coords.x * grid_segment_size.x + arc_offset;
					local_position = tfxVec3(std::cosf(th) * half_emitter_size.x + half_emitter_size.x, grid_coords.y * grid_segment_size.y, -std::sinf(th) * half_emitter_size.z + half_emitter_size.z);

					if (!(property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
						grid_coords.x++;
						if (grid_coords.x >= grid_points.x) {
							grid_coords.x = 0.f;
							grid_coords.y++;
							if (grid_coords.y >= grid_points.y) {
								grid_coords.y = 0.f;
							}
						}
					}
				}

			}
			else if (!(property_flags & tfxEmitterPropertyFlags_fill_area)) {
				float th = random_generation.Range(arc_size) + arc_offset;

				local_position = tfxVec3(std::cosf(th) * half_emitter_size.x + half_emitter_size.x, random_generation.Range(half_emitter_size.y), -std::sinf(th) * half_emitter_size.z + half_emitter_size.z);
			}
			else {
				local_position.x = random_generation.Range(0, half_emitter_size.x);
				local_position.y = random_generation.Range(0, half_emitter_size.y);
				local_position.z = random_generation.Range(0, half_emitter_size.z);

				while ((std::pow(local_position.x - half_emitter_size.x, 2) / std::pow(half_emitter_size.x, 2)) + (std::pow(local_position.z - half_emitter_size.z, 2) / std::pow(half_emitter_size.z, 2)) > 1) {
					local_position.x = random_generation.Range(0, half_emitter_size.x);
					local_position.z = random_generation.Range(0, half_emitter_size.z);
				}
			}

			//----TForm and Emission
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				local_position = mmTransformVector3(matrix, local_position + handle);
				local_position = lerp_position + local_position * scale;
			}

			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleWeight(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const float weight = pm.emitters.weight[emitter_index];
		const float weight_variation = pm.emitters.weight_variation[emitter_index];
		tfxEffectLibrary *library = pm.emitters.library[emitter_index];
		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[emitter_index];
		const float first_weight_value = library->emitter_attributes[emitter_attributes].overtime.weight.GetFirstValue() * tfxUPDATE_TIME;

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &weight_acceleration = entry->particle_data->weight_acceleration[index];
			float &base_weight = entry->particle_data->base_weight[index];

			//----Weight
			if (weight) {
				base_weight = weight;
				if (weight_variation > 0) {
					base_weight += random_generation.Range(-weight_variation, weight_variation);
				}
			}
			else {
				base_weight = 0;
			}
			weight_acceleration = base_weight * first_weight_value;
		}

	}

	void SpawnParticleVelocity(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		tfxEffectLibrary *library = pm.emitters.library[emitter_index];
		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[emitter_index];
		const float velocity = pm.emitters.velocity[emitter_index];
		const float velocity_variation = pm.emitters.velocity_variation[emitter_index];
		const float velocity_adjuster = pm.emitters.velocity_adjuster[emitter_index];
		const float first_velocity_value = library->emitter_attributes[emitter_attributes].overtime.velocity.GetFirstValue() * tfxUPDATE_TIME;

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &base_velocity = entry->particle_data->base_velocity[index];

			//----Velocity
			float velocity_scale = first_velocity_value * velocity_adjuster;
			base_velocity = velocity + random_generation.Range(-velocity_variation, velocity_variation);
		}

	}

	void SpawnParticleRoll(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const tfxAngleSettingFlags angle_settings = properties.angle_settings[property_index];
		const float angle_roll_offset = properties.angle_offsets[property_index].roll;

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &roll = entry->particle_data->local_rotations[index].roll;

			roll = 0;
			if (angle_settings & tfxAngleSettingFlags_random_roll) {
				roll = random_generation.Range(angle_roll_offset);
			}
			else if (angle_settings & tfxAngleSettingFlags_specify_roll) {
				roll = angle_roll_offset;
			}
			else {
				roll = 0;
			}
		}

	}

	void SpawnParticleMicroUpdate2d(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const float splatter = pm.emitters.splatter[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 scale = pm.emitters.scale[emitter_index];

		if (splatter) {
			for (int i = 0; i != entry->amount_to_spawn; ++i) {
				tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
				tfxVec3 &local_position = entry->particle_data->local_position[index];

				//----Splatter
				float splattertemp = splatter;
				float splatx = random_generation.Range(-splatter, splatter);
				float splaty = random_generation.Range(-splatter, splatter);

				while (GetDistance(0, 0, splatx, splaty) >= splattertemp && splattertemp > 0) {
					splatx = random_generation.Range(-splatter, splatter);
					splaty = random_generation.Range(-splatter, splatter);
				}

				if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
					local_position.x += splatx * scale.x;
					local_position.y += splaty * scale.y;
				}
				else {
					local_position.x += splatx;
					local_position.y += splaty;
				}
			}
		}

		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[emitter_index];
		tfxEffectLibrary *library = pm.emitters.library[emitter_index];
		const float first_velocity_value = library->emitter_attributes[emitter_attributes].overtime.velocity.GetFirstValue() * tfxUPDATE_TIME;
		const float first_weight_value = library->emitter_attributes[emitter_attributes].overtime.weight.GetFirstValue() * tfxUPDATE_TIME;
		const tfxAngleSettingFlags angle_settings = properties.angle_settings[property_index];
		const float angle_roll_offset = properties.angle_offsets[property_index].roll;
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxVec3 emitter_world_rotations = pm.emitters.world_rotations[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxEmissionType emission_type = properties.emission_type[property_index];
		const tfxVec2 emitter_size = pm.emitters.emitter_size[emitter_index].xy();
		auto transform_particle_callback2d = pm.emitters.transform_particle_callback2d[emitter_index];

		//Micro Update
		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			const float base_weight = entry->particle_data->base_weight[index];
			float &weight_acceleration = entry->particle_data->weight_acceleration[index];
			float &roll = entry->particle_data->local_rotations[index].roll;
			tfxVec3 &local_position = entry->particle_data->local_position[index];
			tfxVec3 &captured_position = entry->particle_data->captured_position[index];
			float &velocity_normal = entry->particle_data->velocity_normal[index].x;
			float &base_velocity = entry->particle_data->base_velocity[index];

			float micro_time = tfxUPDATE_TIME * (1.f - tween);

			tfxSpriteTransform2d sprite_transform;
			float direction = 0;

			if (angle_settings & tfxAngleSettingFlags_align_roll && property_flags & tfxEmitterPropertyFlags_edge_traversal)
				sprite_transform.rotation = roll = angle_roll_offset;

			bool line = property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == tfxLine;

			if (!line && !(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				TransformParticlePosition(local_position.xy(), roll, sprite_transform.position, sprite_transform.rotation, emitter_world_rotations, matrix, handle, scale, emitter_world_position);
				captured_position = sprite_transform.position;
			}

			if (!line) {
				direction = velocity_normal = GetEmissionDirection2d(pm, library, property_index, emitter_index, local_position.xy(), sprite_transform.position, emitter_size) + library->emitter_attributes[emitter_attributes].overtime.direction.GetFirstValue();
			}

			weight_acceleration += base_weight * first_weight_value * micro_time;
			//----Velocity Changes
			tfxVec2 velocity_normal_tmp;
			velocity_normal_tmp.x = std::sinf(direction);
			velocity_normal_tmp.y = -std::cosf(direction);
			tfxVec2 current_velocity = base_velocity * first_velocity_value * velocity_normal_tmp;
			current_velocity.y += weight_acceleration;
			local_position += current_velocity * micro_time;
			if (line || property_flags & tfxEmitterPropertyFlags_relative_position) {
				tfxVec2 rotatevec = mmTransformVector(matrix, tfxVec2(local_position.x, local_position.y) + handle.xy());
				captured_position = sprite_transform.captured_position = emitter_captured_position.xy() + rotatevec * scale.xy();
				transform_particle_callback2d(local_position.xy(), roll, sprite_transform.position, sprite_transform.rotation, emitter_world_rotations, matrix, handle, scale, tfxVec3(emitter_world_position.x, emitter_world_position.y, 0.f));
			}
			//end micro update

			if ((angle_settings & tfxAngleSettingFlags_align_roll || angle_settings & tfxAngleSettingFlags_align_with_emission) && !line) {
				roll = GetVectorAngle(velocity_normal_tmp.x, velocity_normal_tmp.y) + angle_roll_offset;
			}

			tween += entry->qty_step_size;
		}
	}

	void SpawnParticleMicroUpdate3d(tfxWorkQueue *queue, void *data) {
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const float splatter = pm.emitters.splatter[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 scale = pm.emitters.scale[emitter_index];

		if (splatter) {
			for (int i = 0; i != entry->amount_to_spawn; ++i) {
				tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
				tfxVec3 &local_position = entry->particle_data->local_position[index];

				//----Splatter
				if (splatter) {
					float splatx = random_generation.Range(-splatter, splatter);
					float splaty = random_generation.Range(-splatter, splatter);
					float splatz = random_generation.Range(-splatter, splatter);

					while (std::powf(splatx / splatter, 2.f) + std::powf(splaty / splatter, 2.f) + std::powf(splatz / splatter, 2.f) > 1.f) {
						splatx = random_generation.Range(-splatter, splatter);
						splaty = random_generation.Range(-splatter, splatter);
						splatz = random_generation.Range(-splatter, splatter);
					}

					if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
						local_position.x += splatx * scale.x;
						local_position.y += splaty * scale.y;
						local_position.z += splatz * scale.z;
					}
					else {
						local_position.x += splatx;
						local_position.y += splaty;
						local_position.z += splatz;
					}
				}

			}
		}

		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[emitter_index];
		tfxEffectLibrary *library = pm.emitters.library[emitter_index];
		const float first_velocity_value = library->emitter_attributes[emitter_attributes].overtime.velocity.GetFirstValue() * tfxUPDATE_TIME;
		const float first_weight_value = library->emitter_attributes[emitter_attributes].overtime.weight.GetFirstValue() * tfxUPDATE_TIME;
		const float first_stretch_value = library->emitter_attributes[emitter_attributes].overtime.stretch.GetFirstValue() * tfxUPDATE_TIME;
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxVec3 emitter_size = pm.emitters.emitter_size[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxEmissionType emission_type = properties.emission_type[property_index];
		const tfxVectorAlignType vector_align_type = properties.vector_align_type[property_index];
		const bool line = property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == tfxLine;
		const float velocity_adjuster = pm.emitters.velocity_adjuster[emitter_index];
		const float frame = pm.emitters.frame[emitter_index];

		//Micro Update
		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			const float base_weight = entry->particle_data->base_weight[index];
			float &weight_acceleration = entry->particle_data->weight_acceleration[index];
			float &roll = entry->particle_data->local_rotations[index].roll;
			tfxVec3 &local_position = entry->particle_data->local_position[index];
			tfxVec3 &captured_position = entry->particle_data->captured_position[index];
			tfxVec4 &velocity_normal = entry->particle_data->velocity_normal[index];
			const float base_velocity = entry->particle_data->base_velocity[index];

			tfxVec3 world_position;
			if (!line && !(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				if (!(property_flags & tfxEmitterPropertyFlags_relative_position) && !(property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
					world_position = local_position;
				}
				else {
					tfxVec4 rotatevec = mmTransformVector(matrix, local_position + handle);
					world_position = emitter_world_position + rotatevec.xyz() * scale;
				}
				captured_position = world_position;
			}

			//----Velocity
			float emission_pitch = lookup_callback(library->emitter_attributes[emitter_attributes].properties.emission_pitch, frame);
			float emission_yaw = lookup_callback(library->emitter_attributes[emitter_attributes].properties.emission_yaw, frame);

			if (!(property_flags & tfxEmitterPropertyFlags_edge_traversal) || emission_type != tfxLine) {
				velocity_normal = GetEmissionDirection3d(pm, library, property_index, emitter_index, emission_pitch, emission_yaw, local_position, world_position, emitter_size);
			}
			else if (property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == tfxLine) {
				velocity_normal = tfxVec3(0, 1.f, 0.f);
			}
			float velocity_scale = first_velocity_value * velocity_adjuster * base_velocity;

			//Do a micro update
			//A bit hacky but the epsilon after tween just ensures that theres a guaranteed small difference between captured/world positions so that
			//the alignment on the first frame can be calculated
			float micro_time = tfxUPDATE_TIME * (1.f - tween + 0.001f);
			weight_acceleration += base_weight * first_weight_value * micro_time;
			//----Velocity Changes
			tfxVec3 current_velocity = velocity_normal.xyz() * base_velocity * first_velocity_value;
			current_velocity.y -= weight_acceleration;
			if (vector_align_type == tfxVectorAlignType_motion) {
				float l = FastLength(current_velocity * tfxUPDATE_TIME);
				velocity_normal.w = first_stretch_value * l * 10.f;
			}
			current_velocity *= micro_time;
			local_position += current_velocity;
			if (line || property_flags & tfxEmitterPropertyFlags_relative_position) {
				if (!(property_flags & tfxEmitterPropertyFlags_relative_position) && !(property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
					world_position = local_position;
				}
				else {
					tfxVec4 rotatevec = mmTransformVector(matrix, local_position + handle);
					world_position = emitter_world_position + rotatevec.xyz() * scale;
					captured_position = world_position;
				}
			}
			else {
				world_position += current_velocity;
				captured_position = world_position;
			}
			if (pm.flags & tfxEffectManagerFlags_order_by_depth) {
				entry->particle_data->depth[index] = LengthVec3NoSqR(world_position - pm.camera_position);
			}
			tween += entry->qty_step_size;
		}
	}

	void UpdateEmitterState(tfxParticleManager &pm, tfxU32 index, tfxU32 parent_index, const tfxParentSpawnControls &parent_spawn_controls) {
		tfxPROFILE;
		tfxEffectLibrary *library = pm.emitters.library[index];
		tfxEmitterPropertiesSoA &properties = pm.emitters.library[index]->emitter_properties;

		const float frame = pm.emitters.frame[index];
		const float age = pm.emitters.age[index];
		const tfxU32 global_attributes = pm.effects.global_attributes[parent_index];
		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[index];
		const tfxU32 transform_attributes = pm.emitters.transform_attributes[index];
		const tfxU32 property_index = pm.emitters.properties_index[index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[index];
		const tfxEmitterPropertyFlags parent_property_flags = pm.effects.property_flags[parent_index];
		tfxVec3 &translation = pm.emitters.translation[index];
		tfxVec3 &local_rotations = pm.emitters.local_rotations[index];
		tfxVec3 &scale = pm.emitters.scale[index];
		tfxVec3 &handle = pm.emitters.handle[index];
		tfxVec3 &emitter_size = pm.emitters.emitter_size[index];
		float &overal_scale = pm.emitters.overal_scale[index];
		float &velocity_adjuster = pm.emitters.velocity_adjuster[index];
		float &stretch = pm.emitters.stretch[index];
		float &intensity = pm.emitters.intensity[index];
		float &spawn_quantity = pm.emitters.spawn_quantity[index];

		tfxVec3 &parent_scale = pm.effects.scale[parent_index];

		pm.emitters.life[index] = lookup_callback(library->emitter_attributes[emitter_attributes].base.life, frame) * parent_spawn_controls.life;
		pm.emitters.life_variation[index] = lookup_callback(library->emitter_attributes[emitter_attributes].variation.life, frame) * parent_spawn_controls.life;
		if (!(property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			pm.emitters.size[index].x = lookup_callback(library->emitter_attributes[emitter_attributes].base.width, frame) * parent_spawn_controls.size_x;
			pm.emitters.size[index].y = lookup_callback(library->emitter_attributes[emitter_attributes].base.height, frame) * parent_spawn_controls.size_y;
		}
		else {
			pm.emitters.size[index].x = lookup_callback(library->emitter_attributes[emitter_attributes].base.width, frame);
			if (parent_property_flags & tfxEmitterPropertyFlags_global_uniform_size)
				pm.emitters.size[index].y = pm.emitters.size[index].x * parent_spawn_controls.size_x;
			else
				pm.emitters.size[index].y = pm.emitters.size[index].x * parent_spawn_controls.size_y;
			pm.emitters.size[index].x *= parent_spawn_controls.size_x;
		}
		pm.emitters.size_variation[index].x = lookup_callback(library->emitter_attributes[emitter_attributes].variation.width, frame) * parent_spawn_controls.size_x;
		pm.emitters.size_variation[index].y = lookup_callback(library->emitter_attributes[emitter_attributes].variation.height, frame) * parent_spawn_controls.size_y;
		pm.emitters.velocity[index] = lookup_callback(library->emitter_attributes[emitter_attributes].base.velocity, frame) * parent_spawn_controls.velocity;
		pm.emitters.velocity_variation[index] = lookup_callback(library->emitter_attributes[emitter_attributes].variation.velocity, frame) * parent_spawn_controls.velocity;
		pm.emitters.spin[index] = lookup_callback(library->emitter_attributes[emitter_attributes].base.spin, frame) * parent_spawn_controls.spin;
		pm.emitters.spin_variation[index] = lookup_callback(library->emitter_attributes[emitter_attributes].variation.spin, frame) * parent_spawn_controls.spin;
		pm.emitters.splatter[index] = lookup_callback(library->emitter_attributes[emitter_attributes].properties.splatter, frame) * parent_spawn_controls.splatter;
		pm.emitters.weight[index] = lookup_callback(library->emitter_attributes[emitter_attributes].base.weight, frame) * parent_spawn_controls.weight;
		pm.emitters.weight_variation[index] = lookup_callback(library->emitter_attributes[emitter_attributes].variation.weight, frame) * parent_spawn_controls.weight;
		pm.emitters.noise_offset_variation[index] = lookup_callback(library->emitter_attributes[emitter_attributes].variation.noise_offset, frame);
		pm.emitters.noise_offset[index] = lookup_callback(library->emitter_attributes[emitter_attributes].base.noise_offset, frame);
		pm.emitters.noise_resolution[index] = lookup_callback(library->emitter_attributes[emitter_attributes].variation.noise_resolution, frame);

		translation.x = lookup_callback(library->transform_attributes[transform_attributes].translation_x, frame);
		translation.y = lookup_callback(library->transform_attributes[transform_attributes].translation_y, frame);
		translation.z = lookup_callback(library->transform_attributes[transform_attributes].translation_z, frame);
		overal_scale = pm.effects.overal_scale[parent_index];
		local_rotations.roll = LookupPrecise(library->transform_attributes[transform_attributes].roll, age);
		local_rotations.pitch = LookupPrecise(library->transform_attributes[transform_attributes].pitch, age);
		local_rotations.yaw = LookupPrecise(library->transform_attributes[transform_attributes].yaw, age);
		velocity_adjuster = lookup_callback(library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster, frame);
		intensity = parent_spawn_controls.intensity;
		stretch = pm.effects.stretch[parent_index];
		scale = parent_scale;

		if (!(property_flags & tfxEmitterPropertyFlags_single)) {
			spawn_quantity = lookup_callback(library->emitter_attributes[emitter_attributes].base.amount, frame);
			float amount_variation = lookup_callback(library->emitter_attributes[emitter_attributes].variation.amount, frame);
			spawn_quantity += amount_variation > 0.f ? random_generation.Range(1.f, amount_variation) : 0.f;

			spawn_quantity *= lookup_callback(library->global_graphs[global_attributes].amount, frame);
		}
		else {
			spawn_quantity = (float)properties.spawn_amount[property_index];
		}

		bool is_area = properties.emission_type[property_index] == tfxArea || properties.emission_type[property_index] == tfxEllipse || properties.emission_type[property_index] == tfxCylinder || properties.emission_type[property_index] == tfxIcosphere;

		emitter_size.y = lookup_callback(library->emitter_attributes[emitter_attributes].properties.emitter_height, frame);
		if (is_area) {
			emitter_size.x = lookup_callback(library->emitter_attributes[emitter_attributes].properties.emitter_width, frame);
		}
		else
			emitter_size.x = 0.f;

		if (property_flags & tfxEmitterPropertyFlags_is_3d && is_area) {
			emitter_size.z = lookup_callback(library->emitter_attributes[emitter_attributes].properties.emitter_depth, frame);
		}

		emitter_size *= pm.effects.emitter_size[parent_index];

		if (properties.emission_type[property_index] == tfxEllipse || properties.emission_type[property_index] == tfxCylinder) {
			pm.emitters.arc_size[index] = lookup_callback(library->emitter_attributes[emitter_attributes].properties.arc_size, frame);
			pm.emitters.arc_offset[index] = lookup_callback(library->emitter_attributes[emitter_attributes].properties.arc_offset, frame);
		}

		if (property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center && properties.emission_type[property_index] != tfxPoint) {
			if ((properties.emission_type[property_index] == tfxEllipse || properties.emission_type[property_index] == tfxIcosphere) && property_flags & tfxEmitterPropertyFlags_is_3d)
				handle = emitter_size * 0.f;
			else if (property_flags & tfxEmitterPropertyFlags_is_3d)
				handle = emitter_size * -0.5f;
			else if (properties.emission_type[property_index] == tfxLine)
				handle = emitter_size * 0.5f;
			else
				handle = emitter_size * -0.5f;
		}
		else if (!(property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center)) {
			handle = properties.emitter_handle[property_index];
		}
		else {
			handle = 0.f;
		}

		if (property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {
			if (properties.emission_type[property_index] == tfxArea) {
				if (properties.grid_points[property_index].x > 1)
					pm.emitters.grid_segment_size[index].x = emitter_size.x / (properties.grid_points[property_index].x - 1);
				if (properties.grid_points[property_index].y > 1)
					pm.emitters.grid_segment_size[index].y = emitter_size.y / (properties.grid_points[property_index].y - 1);
				if (properties.grid_points[property_index].z > 1)
					pm.emitters.grid_segment_size[index].z = emitter_size.z / (properties.grid_points[property_index].z - 1);
			}
			else if (properties.emission_type[property_index] == tfxEllipse) {
				if (properties.grid_points[property_index].x > 0)
					pm.emitters.grid_segment_size[index].x = pm.emitters.arc_size[index] / (properties.grid_points[property_index].x);
			}
			else if (properties.emission_type[property_index] == tfxCylinder) {
				if (properties.grid_points[property_index].x > 0)
					pm.emitters.grid_segment_size[index].x = pm.emitters.arc_size[index] / (properties.grid_points[property_index].x);
				if (properties.grid_points[property_index].y > 1)
					pm.emitters.grid_segment_size[index].y = emitter_size.y / (properties.grid_points[property_index].y - 1);
			}
			else if (properties.emission_type[property_index] == tfxLine) {
				if (properties.grid_points[property_index].x > 1)
					pm.emitters.grid_segment_size[index].y = emitter_size.y / (properties.grid_points[property_index].x - 1);
			}
		}

		//if (e.update_emitter_callback)
			//e.update_emitter_callback(pm.emitters, index);

	}

	void UpdateEffectState(tfxParticleManager &pm, tfxU32 index) {
		tfxPROFILE;

		const float frame = pm.effects.frame[index];
		const float age = pm.effects.age[index];
		const tfxU32 global_attributes = pm.effects.global_attributes[index];
		const tfxU32 transform_attributes = pm.effects.transform_attributes[index];
		tfxEffectLibrary *library = pm.effects.library[index];
		tfxVec3 &translation = pm.effects.translation[index];
		tfxVec3 &local_rotations = pm.effects.local_rotations[index];
		tfxVec3 &scale = pm.effects.scale[index];
		tfxVec3 &emitter_size = pm.effects.emitter_size[index];
		float &overal_scale = pm.effects.overal_scale[index];
		float &stretch = pm.effects.stretch[index];

		//If this effect is a sub effect then the graph index will reference the global graphs for the root parent effect
		tfxParentSpawnControls &spawn_controls = pm.effects.spawn_controls[index];
		spawn_controls.life = lookup_callback(library->global_graphs[global_attributes].life, frame);
		if (!(pm.effects.property_flags[index] & tfxEmitterPropertyFlags_global_uniform_size)) {
			spawn_controls.size_x = lookup_callback(library->global_graphs[global_attributes].width, frame);
			spawn_controls.size_y = lookup_callback(library->global_graphs[global_attributes].height, frame);
		}
		else {
			spawn_controls.size_x = lookup_callback(library->global_graphs[global_attributes].width, frame);
			spawn_controls.size_y = spawn_controls.size_x;
		}
		spawn_controls.velocity = lookup_callback(library->global_graphs[global_attributes].velocity, frame);
		spawn_controls.spin = lookup_callback(library->global_graphs[global_attributes].spin, frame);
		spawn_controls.intensity = lookup_callback(library->global_graphs[global_attributes].intensity, frame);
		spawn_controls.splatter = lookup_callback(library->global_graphs[global_attributes].splatter, frame);
		spawn_controls.weight = lookup_callback(library->global_graphs[global_attributes].weight, frame);
		emitter_size.x = lookup_callback(library->global_graphs[global_attributes].emitter_width, frame);
		emitter_size.y = lookup_callback(library->global_graphs[global_attributes].emitter_height, frame);
		emitter_size.z = lookup_callback(library->global_graphs[global_attributes].emitter_depth, frame);
		//We don't want to scale twice when the sub effect is transformed, so the values here are set to 1. That means that the root effect will only control thglobal scale.
		overal_scale = lookup_callback(library->global_graphs[global_attributes].overal_scale, frame);
		if (pm.effects.parent_particle_index[index] == tfxINVALID) {
			scale.x = overal_scale;
			scale.y = overal_scale;
			scale.z = overal_scale;
			local_rotations.roll = LookupPrecise(library->transform_attributes[transform_attributes].roll, age);
			local_rotations.pitch = LookupPrecise(library->transform_attributes[transform_attributes].pitch, age);
			local_rotations.yaw = LookupPrecise(library->transform_attributes[transform_attributes].yaw, age);
		}
		else {
			scale.x = overal_scale;
			scale.y = overal_scale;
			scale.z = overal_scale;
			local_rotations.roll = 0.f;
			local_rotations.pitch = 0.f;
			local_rotations.yaw = 0.f;
		}
		stretch = lookup_callback(library->global_graphs[global_attributes].stretch, frame);
		translation.x = lookup_callback(library->transform_attributes[transform_attributes].translation_x, frame);
		translation.y = lookup_callback(library->transform_attributes[transform_attributes].translation_y, frame);
		translation.z = lookup_callback(library->transform_attributes[transform_attributes].translation_z, frame);

		//if (e.update_effect_callback)
			//e.update_effect_callback(pm, e, spawn_controls);

	}

	void ControlParticleAge(tfxWorkQueue *queue, void *data) {
		tfxParticleAgeWorkEntry *work_entry = static_cast<tfxParticleAgeWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		tfxParticleManager &pm = *work_entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		tfxParticleSoA &bank = pm.particle_arrays[particles_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		tfxU32 single_shot_limit = work_entry->properties->single_shot_limit[property_index];

		const tfxEmitterStateFlags state_flags = pm.emitters.state_flags[emitter_index];

		tfxU32 offset = 0;
		for (int i = work_entry->start_index; i >= 0; --i) {
			const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i);
			float &age = bank.age[index];
			const float &max_age = bank.max_age[index];
			tfxU32 &single_loop_count = bank.single_loop_count[index];
			tfxParticleFlags &flags = bank.flags[index];
			age += tfxFRAME_LENGTH;
			flags |= state_flags & tfxParticleFlags_remove;

			if (flags & tfxParticleFlags_remove || age >= max_age) {
				if (property_flags & tfxEmitterPropertyFlags_single && !(work_entry->pm->flags & tfxEffectManagerFlags_disable_spawning))
					if (++single_loop_count != single_shot_limit) {
						age = 0;
					}
					else {
						flags |= tfxParticleFlags_remove;
					}
				else {
					flags |= tfxParticleFlags_remove;
				}
			}

			if (flags & tfxParticleFlags_remove) {
				offset++;
				if (flags & tfxParticleFlags_has_sub_effects) {
					pm.FreeParticleIndex(bank.particle_index[index]);
				}
			}
			else if (offset > 0) {
				tfxU32 next_index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i + offset);
				if(flags & tfxParticleFlags_has_sub_effects)
					pm.particle_indexes[bank.particle_index[index]] = MakeParticleID(particles_index, next_index);

				bank.parent_index[next_index] = bank.parent_index[index];
				bank.sprite_index[next_index] = bank.sprite_index[index];
				bank.particle_index[next_index] = bank.particle_index[index];
				bank.flags[next_index] = bank.flags[index];
				bank.age[next_index] = bank.age[index];	
				bank.max_age[next_index] = bank.max_age[index];
				bank.local_position[next_index] = bank.local_position[index];
				bank.captured_position[next_index] = bank.captured_position[index];
				bank.local_rotations[next_index] = bank.local_rotations[index];
				bank.velocity_normal[next_index] = bank.velocity_normal[index];
				bank.stretch[next_index] = bank.stretch[index];
				bank.weight_acceleration[next_index] = bank.weight_acceleration[index];
				bank.base_weight[next_index] = bank.base_weight[index];
				bank.base_velocity[next_index] = bank.base_velocity[index];
				bank.base_spin[next_index] = bank.base_spin[index];
				bank.noise_offset[next_index] = bank.noise_offset[index];		
				bank.noise_resolution[next_index] = bank.noise_resolution[index];
				bank.color[next_index] = bank.color[index];			
				bank.intensity[next_index] = bank.intensity[index];
				bank.image_frame[next_index] = bank.image_frame[index];
				bank.base_size[next_index] = bank.base_size[index];
				bank.single_loop_count[next_index] = bank.single_loop_count[index];	
			}

		}

		if (offset) {
			Bump(&work_entry->pm->particle_array_buffers[particles_index], offset);
		}

	}

	void ControlParticlePosition2d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		tfxParticleManager &pm = *work_entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[particles_index];

		const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];
		const tfxVec3 emitter_size = pm.emitters.emitter_size[emitter_index];
		const float overal_scale = pm.emitters.overal_scale[emitter_index];
		const float angle_offset = pm.emitters.angle_offsets[emitter_index].roll;
		const float velocity_adjuster = pm.emitters.velocity_adjuster[emitter_index];

		for (int i = work_entry->start_index; i != work_entry->end_index; ++i) {
			const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i);

			const float age = bank.age[index];
			const float max_age = bank.max_age[index];
			const tfxU32 lookup_frame = static_cast<tfxU32>((age / max_age * work_entry->graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

			const float base_weight = bank.base_weight[index];
			const float base_velocity = bank.base_velocity[index];
			const float base_spin = bank.base_spin[index];
			const float noise_offset = bank.noise_offset[index];
			const float noise_resolution = bank.noise_resolution[index];
			const float angle = bank.velocity_normal[index].x;
			float &weight_acceleration = bank.weight_acceleration[index];
			tfxVec3 &local_position = bank.local_position[index];
			tfxVec3 &captured_position = bank.captured_position[index];
			tfxVec3 &local_rotations = bank.local_rotations[index];
			tfxParticleFlags &flags = bank.flags[index];

			const float lookup_velocity = work_entry->graphs->velocity.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->velocity.lookup.last_frame)] * velocity_adjuster;
			const float lookup_velocity_turbulance = work_entry->graphs->velocity_turbulance.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->velocity_turbulance.lookup.last_frame)];
			const float lookup_direction = work_entry->graphs->direction.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->direction.lookup.last_frame)] + angle;
			const float lookup_noise_resolution = work_entry->graphs->noise_resolution.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->noise_resolution.lookup.last_frame)] * noise_resolution;
			const float lookup_weight = work_entry->graphs->weight.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->weight.lookup.last_frame)];
			const float lookup_spin = work_entry->graphs->spin.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->spin.lookup.last_frame)] * base_spin;

			float direction = 0;

			tfxVec2 mr_vec;
			if (emitter_flags & tfxEmitterStateFlags_not_line) {
				direction = lookup_direction;
			}

			if (lookup_velocity_turbulance) {
				float eps = 0.0001f;
				float eps2 = 0.0001f * 2.f;

				float x = local_position.x / lookup_noise_resolution + noise_offset;
				float y = local_position.y / lookup_noise_resolution + noise_offset;

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
			weight_acceleration += base_weight * lookup_weight * tfxUPDATE_TIME;

			//----Velocity Changes
			tfxVec2 velocity_normal;
			velocity_normal.x = std::sinf(direction);
			velocity_normal.y = -std::cosf(direction);

			tfxVec2 current_velocity = (base_velocity * lookup_velocity) * velocity_normal;
			current_velocity += mr_vec;
			current_velocity.y += weight_acceleration;
			current_velocity *= tfxUPDATE_TIME;

			//----Spin and angle Changes
			float spin = 0;
			if (emitter_flags & tfxEmitterStateFlags_can_spin) {
				spin = lookup_spin;
			}

			//----Rotation
			if (emitter_flags & tfxEmitterStateFlags_align_with_velocity) {
				tfxVec2 vd = current_velocity.IsNill() ? velocity_normal : current_velocity;
				local_rotations.roll = GetVectorAngle(vd.x, vd.y) + angle_offset;
			}
			else {
				local_rotations.roll += spin * tfxUPDATE_TIME;
			}

			local_position += current_velocity * overal_scale;

			//Lines - Reposition if the particle is travelling along a line
			tfxVec2 offset = velocity_normal * emitter_size.y;
			float length = std::fabsf(local_position.y);
			float emitter_length = emitter_size.y;
			bool line_and_kill = (emitter_flags & tfxEmitterStateFlags_is_line_traversal) && (emitter_flags & tfxEmitterStateFlags_kill) && length > emitter_length;
			bool line_and_loop = (emitter_flags & tfxEmitterStateFlags_is_line_traversal) && (emitter_flags & tfxEmitterStateFlags_loop) && length > emitter_length;
			if (line_and_loop) {
				local_position.y -= offset.y;
				flags |= tfxParticleFlags_capture_after_transform;
			}
			else if (line_and_kill) {
				flags |= tfxParticleFlags_remove;
			}

		}

		tfxU32 running_sprite_index = work_entry->sprites_index;

		tfxVec3 &e_captured_position = pm.emitters.captured_position[emitter_index];
		tfxVec3 &e_world_position = pm.emitters.world_position[emitter_index];
		tfxVec3 &e_world_rotations = pm.emitters.world_rotations[emitter_index];
		tfxVec3 &e_handle = pm.emitters.handle[emitter_index];
		tfxMatrix4 &e_matrix = pm.emitters.matrix[emitter_index];
		tfxVec3 &e_scale = pm.emitters.scale[emitter_index];

		for (int i = work_entry->start_index; i != work_entry->end_index; ++i) {
			const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i);
			const tfxVec3 &local_position = bank.local_position[index];
			const tfxVec3 &local_rotations = bank.local_rotations[index];
			tfxVec3 &captured_position = bank.captured_position[index];
			tfxParticleFlags &flags = bank.flags[index];

			tfxParticleSprite2d &s = (*work_entry->sprites2d)[running_sprite_index++];
			auto transform_particle_callback2d = pm.emitters.transform_particle_callback2d[emitter_index];
			if (flags & tfxParticleFlags_capture_after_transform) {
				transform_particle_callback2d(local_position.xy(), local_rotations.roll, s.transform.position, s.transform.rotation, e_world_rotations, e_matrix, e_handle, e_scale, tfxVec3(e_captured_position.x, e_captured_position.y, 0.f));
				captured_position = s.transform.position;
				transform_particle_callback2d(local_position.xy(), local_rotations.roll, s.transform.position, s.transform.rotation, e_world_rotations, e_matrix, e_handle, e_scale, tfxVec3(e_world_position.x, e_world_position.y, 0.f));
				flags &= ~tfxParticleFlags_capture_after_transform;
			}
			else {
				transform_particle_callback2d(local_position.xy(), local_rotations.roll, s.transform.position, s.transform.rotation, e_world_rotations, e_matrix, e_handle, e_scale, tfxVec3(e_world_position.x, e_world_position.y, 0.f));
			}
			s.transform.captured_position = captured_position.xy();
			captured_position = s.transform.position;

		}
	}

	void ControlParticlePosition3d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		tfxParticleManager &pm = *work_entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[particles_index];

		const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];
		const tfxVec3 emitter_size = pm.emitters.emitter_size[emitter_index];
		const float overal_scale = pm.emitters.overal_scale[emitter_index];
		const tfxVec3 angle_offsets = pm.emitters.angle_offsets[emitter_index];
		const float velocity_adjuster = pm.emitters.velocity_adjuster[emitter_index];
		const float stretch = pm.emitters.stretch[emitter_index];

		tfxU32 running_sprite_index = work_entry->sprites_index;

		for (int i = work_entry->start_index; i != work_entry->end_index; ++i) {
			const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i);

			const float age = bank.age[index];
			const float max_age = bank.max_age[index];
			const tfxU32 lookup_frame = static_cast<tfxU32>((age / max_age * work_entry->graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

			const float base_weight = bank.base_weight[index];
			const float base_velocity = bank.base_velocity[index];
			const float base_spin = bank.base_spin[index];
			const float base_noise_offset = bank.noise_offset[index];
			const float noise_resolution = bank.noise_resolution[index];
			const float angle = bank.velocity_normal[index].x;
			float &weight_acceleration = bank.weight_acceleration[index];
			tfxVec3 &local_position = bank.local_position[index];
			tfxVec3 &captured_position = bank.captured_position[index];
			tfxVec3 &local_rotations = bank.local_rotations[index];
			tfxVec4 &velocity_normal = bank.velocity_normal[index];
			tfxParticleFlags &flags = bank.flags[index];

			const float lookup_velocity = work_entry->graphs->velocity.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->velocity.lookup.last_frame)] * velocity_adjuster;
			const float lookup_velocity_turbulance = work_entry->graphs->velocity_turbulance.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->velocity_turbulance.lookup.last_frame)];
			const float lookup_direction = work_entry->graphs->direction.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->direction.lookup.last_frame)] + angle;
			const float lookup_noise_resolution = work_entry->graphs->noise_resolution.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->noise_resolution.lookup.last_frame)] * noise_resolution;
			const float lookup_weight = work_entry->graphs->weight.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->weight.lookup.last_frame)];
			const float lookup_spin = work_entry->graphs->spin.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->spin.lookup.last_frame)] * base_spin;
			const float lookup_stretch = work_entry->graphs->stretch.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->stretch.lookup.last_frame)];

			float direction = 0.f;
			float mr_angle = 0.f;

			tfxVec3 mr_vec;

			if (lookup_velocity_turbulance) {
				float eps = 0.001f;
				float eps2 = 0.001f * 2.f;

				float noise_offset = base_noise_offset * overal_scale;

				float x = local_position.x / lookup_noise_resolution + noise_offset;
				float y = local_position.y / lookup_noise_resolution + noise_offset;
				float z = local_position.z / lookup_noise_resolution + noise_offset;

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
			weight_acceleration += base_weight * lookup_weight * tfxUPDATE_TIME;

			//----Velocity Changes
			float velocity_scalar = base_velocity * lookup_velocity;
			tfxVec3 current_velocity = velocity_normal.xyz() * velocity_scalar;
			current_velocity += mr_vec;
			current_velocity.y -= weight_acceleration;
			velocity_normal.w = lookup_stretch * stretch;
			current_velocity *= tfxUPDATE_TIME;

			//----Spin and angle Changes
			float spin = 0;
			if (emitter_flags & tfxEmitterStateFlags_can_spin) {
				spin = lookup_spin;
			}

			//----Rotation
			if (flags & tfxEmitterStateFlags_align_with_velocity) {
				tfxVec3 vd = current_velocity.IsNill() ? velocity_normal.xyz() : current_velocity;
				local_rotations.roll = GetVectorAngle(vd.x, vd.y) + angle_offsets.roll;
			}
			else {
				local_rotations.roll += spin * tfxUPDATE_TIME;
			}

			//----Position
			local_position += current_velocity * overal_scale;

			//Lines - Reposition if the particle is travelling along a line
			tfxVec3 offset = velocity_normal.xyz() * emitter_size.y;
			float length = std::fabsf(local_position.y);
			float emitter_length = emitter_size.y;
			bool line_and_kill = (flags & tfxEmitterStateFlags_is_line_traversal) && (flags & tfxEmitterStateFlags_kill) && length > emitter_length;
			bool line_and_loop = (flags & tfxEmitterStateFlags_is_line_traversal) && (flags & tfxEmitterStateFlags_loop) && length > emitter_length;
			if (line_and_loop) {
				local_position.y -= offset.y;
				flags |= tfxParticleFlags_capture_after_transform;
			}
			else if (line_and_kill) {
				flags |= tfxParticleFlags_remove;
			}

		}

		running_sprite_index = work_entry->sprites_index;

		tfxVec3 &e_captured_position = pm.emitters.captured_position[emitter_index];
		tfxVec3 &e_world_position = pm.emitters.world_position[emitter_index];
		tfxVec3 &e_world_rotations = pm.emitters.world_rotations[emitter_index];
		tfxVec3 &e_handle = pm.emitters.handle[emitter_index];
		tfxMatrix4 &e_matrix = pm.emitters.matrix[emitter_index];
		tfxVec3 &e_scale = pm.emitters.scale[emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmissionType emission_type = work_entry->properties->emission_type[property_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVectorAlignType vector_align_type = work_entry->properties->vector_align_type[property_index];
		const tfxBillboardingOptions billboard_option = work_entry->properties->billboard_option[property_index];

		for (int i = work_entry->start_index; i != work_entry->end_index; ++i) {
			const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i);
			const tfxVec3 local_position = bank.local_position[index];
			const tfxVec3 local_rotations = bank.local_rotations[index];
			tfxVec4 &velocity_normal = bank.velocity_normal[index];
			tfxVec3 &captured_position = bank.captured_position[index];
			tfxParticleFlags &flags = bank.flags[index];

			tfxParticleSprite3d &s = (*work_entry->sprites3d)[running_sprite_index++];
			auto transform_particle_callback3d = pm.emitters.transform_particle_callback3d[emitter_index];
			if (flags & tfxParticleFlags_capture_after_transform) {
				transform_particle_callback3d(local_position, local_rotations, s.transform.position, s.transform.rotations, e_world_rotations, e_matrix, e_handle, e_scale, e_captured_position);
				captured_position = s.transform.position;
				transform_particle_callback3d(local_position, local_rotations, s.transform.position, s.transform.rotations, e_world_rotations, e_matrix, e_handle, e_scale, e_world_position);
				flags &= ~tfxParticleFlags_capture_after_transform;
			}
			else {
				transform_particle_callback3d(local_position, local_rotations, s.transform.position, s.transform.rotations, e_world_rotations, e_matrix, e_handle, e_scale, e_world_position);
			}

			tfxVec3 alignment_vector;
			bool line = property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == tfxLine;
			if (vector_align_type == tfxVectorAlignType_motion) {
				alignment_vector = s.transform.position - captured_position;
				float l = FastLength(alignment_vector);
				velocity_normal.w *= l * 10.f;
				alignment_vector = FastNormalizeVec(alignment_vector);
			}
			else if (vector_align_type == tfxVectorAlignType_emission && property_flags & tfxEmitterPropertyFlags_relative_position) {
				alignment_vector = mmTransformVector3(e_matrix, velocity_normal.xyz());
			}
			else if (vector_align_type == tfxVectorAlignType_emission) {
				alignment_vector = velocity_normal.xyz();
			}
			else if (vector_align_type == tfxVectorAlignType_emitter) {
				alignment_vector = mmTransformVector(e_matrix, tfxVec4(0.f, 1.f, 0.f, 0.f)).xyz();
			}

			s.transform.captured_position = captured_position;
			alignment_vector.y += 0.002f;	//We don't want a 0 alignment normal
			s.alignment = Pack10bit(alignment_vector, billboard_option & 0x00000003);
			s.stretch = velocity_normal.w;

			captured_position = s.transform.position;

		}
	}

	void ControlParticleSize2d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		const tfxU32 particles_index = work_entry->pm->emitters.particles_index[emitter_index];
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = pm.particle_arrays[particles_index];

		const float overal_scale = pm.emitters.overal_scale[emitter_index];
		const float velocity_adjuster = pm.emitters.velocity_adjuster[emitter_index];
		const float stretch = pm.emitters.stretch[emitter_index];
		const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];
		const tfxVec2 image_size = pm.emitters.image_size[emitter_index];
		const tfxVec2 image_handle = pm.emitters.image_handle[emitter_index];

		tfxU32 running_sprite_index = work_entry->sprites_index;

		for (int i = work_entry->start_index; i != work_entry->end_index; ++i) {
			const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i);
			const float age = bank.age[index];
			const float max_age = bank.max_age[index];
			const tfxU32 lookup_frame = static_cast<tfxU32>((age / max_age * work_entry->graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

			const tfxVec2 base_size = bank.base_size[index];
			const float base_velocity = bank.base_velocity[index];
			const float weight_acceleration = bank.weight_acceleration[index];

			float lookup_velocity = work_entry->graphs->velocity.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->velocity.lookup.last_frame)] * velocity_adjuster;
			float lookup_stretch = work_entry->graphs->stretch.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->stretch.lookup.last_frame)];
			float lookup_width = work_entry->graphs->width.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->width.lookup.last_frame)];
			float lookup_height = work_entry->graphs->height.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->height.lookup.last_frame)];

			//----Size Changes
			tfxVec2 scale;
			scale.x = base_size.x * lookup_width;
			if (scale.x < 0.f)
				scale.x = scale.x;

			//----Stretch Changes
			float velocity = std::fabsf(lookup_velocity * base_velocity + weight_acceleration);
			if (emitter_flags & tfxEmitterStateFlags_lifetime_uniform_size) {
				scale.y = (lookup_width * (base_size.y + (velocity * lookup_stretch * stretch))) / image_size.y;
				if (emitter_flags & tfxEmitterPropertyFlags_base_uniform_size && scale.y < scale.x)
					scale.y = scale.x;
			}
			else
				scale.y = (lookup_height * (base_size.y + (velocity * lookup_stretch * stretch))) / image_size.y;

			tfxParticleSprite2d &s = (*work_entry->sprites2d)[running_sprite_index++];
			s.transform.scale = scale * overal_scale;
			s.handle = image_handle;
		}

	}

	void ControlParticleColor2d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		const tfxU32 particles_index = work_entry->pm->emitters.particles_index[emitter_index];
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[particles_index];

		const float global_intensity = pm.emitters.intensity[emitter_index];
		const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];

		tfxU32 running_sprite_index = work_entry->sprites_index;

		for (int i = work_entry->start_index; i != work_entry->end_index; ++i) {
			const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i);
			const float age = bank.age[index];
			const float max_age = bank.max_age[index];
			const tfxU32 lookup_frame = static_cast<tfxU32>((age / max_age * work_entry->graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

			const float lookup_opacity = work_entry->graphs->blendfactor.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->blendfactor.lookup.last_frame)];
			const float lookup_intensity = work_entry->graphs->intensity.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->intensity.lookup.last_frame)];

			tfxRGBA8 &color = bank.color[index];
			float &intensity = bank.intensity[index];

			//----Color changes
			color.a = unsigned char(255.f * lookup_opacity);
			intensity = lookup_intensity * global_intensity;
			if (!(emitter_flags & tfxEmitterStateFlags_random_color)) {
				const float lookup_red = work_entry->graphs->red.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->red.lookup.last_frame)];
				const float lookup_green = work_entry->graphs->green.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->green.lookup.last_frame)];
				const float lookup_blue = work_entry->graphs->blue.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->blue.lookup.last_frame)];
				color.r = unsigned char(255.f * lookup_red);
				color.g = unsigned char(255.f * lookup_green);
				color.b = unsigned char(255.f * lookup_blue);
			}

			color = tfxRGBA8(color.r, color.g, color.b, color.a);

			tfxParticleSprite2d &s = (*work_entry->sprites2d)[running_sprite_index++];
			s.color = color;
			s.intensity = intensity;
		}

	}

	void ControlParticleImageFrame2d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		const tfxU32 particles_index = work_entry->pm->emitters.particles_index[emitter_index];
		const tfxU32 property_index = work_entry->pm->emitters.properties_index[emitter_index];
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = pm.particle_arrays[particles_index];
		tfxImageData *image = work_entry->properties->image[property_index];

		float image_frame_rate = pm.emitters.image_frame_rate[emitter_index];
		float end_frame = pm.emitters.end_frame[emitter_index];
		tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];

		tfxU32 running_sprite_index = work_entry->sprites_index;

		for (int i = work_entry->start_index; i != work_entry->end_index; ++i) {
			const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i);
			const float age = bank.age[index];
			const float max_age = bank.max_age[index];
			const tfxU32 lookup_frame = static_cast<tfxU32>((age / max_age * work_entry->graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

			float &image_frame = bank.image_frame[index];
			tfxU32 &sprites_index = bank.sprite_index[index];
			sprites_index = (work_entry->layer << 28) + running_sprite_index;
			tfxParticleSprite2d &s = (*work_entry->sprites2d)[running_sprite_index++];

			//----Image animation
			image_frame += image_frame_rate * tfxUPDATE_TIME;
			image_frame = (emitter_flags & tfxEmitterStateFlags_play_once) && image_frame > end_frame ? image_frame = end_frame : image_frame;
			image_frame = (emitter_flags & tfxEmitterStateFlags_play_once) && image_frame < 0 ? image_frame = 0 : image_frame;
			image_frame = std::fmodf(image_frame, end_frame + 1);

			s.image_frame = (tfxU32)image_frame;
			s.image_ptr = image->ptr;
		}

	}

	void ControlParticleSize3d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		const tfxU32 particles_index = work_entry->pm->emitters.particles_index[emitter_index];
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = pm.particle_arrays[particles_index];

		const float overal_scale = pm.emitters.overal_scale[emitter_index];
		const float velocity_adjuster = pm.emitters.velocity_adjuster[emitter_index];
		const float stretch = pm.emitters.stretch[emitter_index];
		const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];
		const tfxVec2 image_size = pm.emitters.image_size[emitter_index];
		const tfxVec2 image_handle = pm.emitters.image_handle[emitter_index];
		const tfxU32 width_last_frame = work_entry->graphs->width.lookup.last_frame;
		const tfxU32 height_last_frame = work_entry->graphs->height.lookup.last_frame;

		tfxU32 running_sprite_index = work_entry->sprites_index;

		for (int i = work_entry->start_index; i != work_entry->end_index; ++i) {
			const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i);
			const float age = bank.age[index];
			const float max_age = bank.max_age[index];
			const tfxU32 lookup_frame = static_cast<tfxU32>((age / max_age * work_entry->graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

			const tfxVec2 base_size = bank.base_size[index];

			float lookup_width = work_entry->graphs->width.lookup.values[std::min<tfxU32>(lookup_frame, width_last_frame)];
			float lookup_height = work_entry->graphs->height.lookup.values[std::min<tfxU32>(lookup_frame, height_last_frame)];

			//----Size Changes
			tfxVec2 scale;
			scale.x = base_size.x * lookup_width;
			if (scale.x < 0.f)
				scale.x = scale.x;

			if (emitter_flags & tfxEmitterStateFlags_lifetime_uniform_size) {
				scale.y = lookup_width * base_size.y;
				if (emitter_flags & tfxEmitterPropertyFlags_base_uniform_size && scale.y < scale.x)
					scale.y = scale.x;
			}
			else
				scale.y = lookup_height * base_size.y;

			tfxParticleSprite3d &s = (*work_entry->sprites3d)[running_sprite_index++];
			s.transform.scale = scale * overal_scale;
			s.handle = image_handle;
		}

	}

	void ControlParticleColor3d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		const tfxU32 particles_index = work_entry->pm->emitters.particles_index[emitter_index];
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[particles_index];

		const float global_intensity = pm.emitters.intensity[emitter_index];
		const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];

		tfxU32 running_sprite_index = work_entry->sprites_index;

		for (int i = work_entry->start_index; i != work_entry->end_index; ++i) {
			const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i);
			const float age = bank.age[index];
			const float max_age = bank.max_age[index];
			const tfxU32 lookup_frame = static_cast<tfxU32>((age / max_age * work_entry->graphs->velocity.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);

			const float lookup_red = work_entry->graphs->red.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->red.lookup.last_frame)];
			const float lookup_green = work_entry->graphs->green.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->green.lookup.last_frame)];
			const float lookup_blue = work_entry->graphs->blue.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->blue.lookup.last_frame)];
			const float lookup_opacity = work_entry->graphs->blendfactor.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->blendfactor.lookup.last_frame)];
			const float lookup_intensity = work_entry->graphs->intensity.lookup.values[std::min<tfxU32>(lookup_frame, work_entry->graphs->intensity.lookup.last_frame)];

			tfxRGBA8 &color = bank.color[index];
			float &intensity = bank.intensity[index];

			//----Color changes
			color.a = unsigned char(255.f * lookup_opacity);
			intensity = lookup_intensity * global_intensity;
			if (!(emitter_flags & tfxEmitterStateFlags_random_color)) {
				color.r = unsigned char(255.f * lookup_red);
				color.g = unsigned char(255.f * lookup_green);
				color.b = unsigned char(255.f * lookup_blue);
			}

			color = tfxRGBA8(color.r, color.g, color.b, color.a);

			tfxParticleSprite3d &s = (*work_entry->sprites3d)[running_sprite_index++];
			s.color = color;
			s.intensity = intensity;
		}

	}

	void ControlParticleImageFrame3d(tfxWorkQueue *queue, void *data) {
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		const tfxU32 particles_index = work_entry->pm->emitters.particles_index[emitter_index];
		const tfxU32 property_index = work_entry->pm->emitters.properties_index[emitter_index];
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = pm.particle_arrays[particles_index];
		tfxImageData *image = work_entry->properties->image[property_index];
		const tfxBillboardingOptions billboard_option = work_entry->properties->billboard_option[property_index];

		float image_frame_rate = pm.emitters.image_frame_rate[emitter_index];
		float end_frame = pm.emitters.end_frame[emitter_index];
		tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];

		tfxU32 running_sprite_index = work_entry->sprites_index;

		for (int i = work_entry->start_index; i != work_entry->end_index; ++i) {
			const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i);

			float &image_frame = bank.image_frame[index];
			tfxU32 &sprites_index = bank.sprite_index[index];
			sprites_index = (work_entry->layer << 28) + running_sprite_index;
			tfxParticleSprite3d &s = (*work_entry->sprites3d)[running_sprite_index++];

			//----Image animation
			image_frame += image_frame_rate * tfxUPDATE_TIME;
			image_frame = (emitter_flags & tfxEmitterStateFlags_play_once) && image_frame > end_frame ? image_frame = end_frame : image_frame;
			image_frame = (emitter_flags & tfxEmitterStateFlags_play_once) && image_frame < 0 ? image_frame = 0 : image_frame;
			image_frame = std::fmodf(image_frame, end_frame + 1);

			s.image_frame_plus = (billboard_option << 24) + (tfxU32)image_frame;
			s.image_ptr = image->ptr;
		}

	}

	void ControlParticles2d(tfxParticleManager &pm, tfxU32 emitter_index, tfxControlWorkEntry &work_entry) {
		tfxPROFILE;

		tfxEffectLibrary *library = pm.emitters.library[emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[emitter_index];
		const tfxU32 sprites_count = pm.emitters.sprites_count[emitter_index];
		const tfxU32 sprites_index = pm.emitters.sprites_index[emitter_index];
		tfxEmitterPropertiesSoA &properties = library->emitter_properties;

		//-------------------------------------------------------
		//Controll what the particle does over the course of
		//it's lifetime
		//-------------------------------------------------------

		work_entry.graphs = &library->emitter_attributes[emitter_attributes].overtime;
		//work_entry.c.particle_update_callback = e.particle_update_callback;
		//work_entry.c.user_data = e.user_data;

		tfxSoABuffer &buffer = pm.particle_array_buffers[pm.emitters.particles_index[emitter_index]];
		int offset = 0;
		tfxU32 amount_to_update = buffer.current_size;
		if (sprites_count < buffer.current_size) {
			amount_to_update = sprites_count;
		}
		
		work_entry.pm = &pm;
		work_entry.sprites_index = sprites_index + work_entry.start_index;
		work_entry.layer = properties.layer[property_index];
		work_entry.sprites2d = &pm.sprites2d[work_entry.layer];

		if (amount_to_update > 0) {
#if tfxMULTITHREADED
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticlePosition2d);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleSize2d);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleColor2d);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleImageFrame2d);
#else
			ControlParticlePosition2d(&pm.work_queue, &work_entry);
			ControlParticleSize2d(&pm.work_queue, &work_entry);
			ControlParticleColor2d(&pm.work_queue, &work_entry);
			ControlParticleImageFrame2d(&pm.work_queue, &work_entry);
#endif
		}

	}

	void ControlParticles3d(tfxParticleManager &pm, tfxU32 emitter_index, tfxControlWorkEntry &work_entry) {
		tfxPROFILE;

		tfxEffectLibrary *library = pm.emitters.library[emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[emitter_index];
		const tfxU32 sprites_count = pm.emitters.sprites_count[emitter_index];
		const tfxU32 sprites_index = pm.emitters.sprites_index[emitter_index];
		tfxEmitterPropertiesSoA &properties = library->emitter_properties;

		//-------------------------------------------------------
		//Controll what the particle does over the course of
		//it's lifetime
		//-------------------------------------------------------

		work_entry.graphs = &library->emitter_attributes[emitter_attributes].overtime;
		//work_entry.c.particle_update_callback = e.particle_update_callback;
		//work_entry.c.user_data = e.user_data;

		tfxSoABuffer &buffer = pm.particle_array_buffers[pm.emitters.particles_index[emitter_index]];
		int offset = 0;
		tfxU32 amount_to_update = buffer.current_size;
		if (sprites_count < buffer.current_size) {
			amount_to_update = sprites_count;
		}

		work_entry.pm = &pm;
		work_entry.sprites_index = sprites_index + work_entry.start_index;
		work_entry.layer = properties.layer[property_index];
		work_entry.sprites3d = &pm.sprites3d[work_entry.layer];

		if (amount_to_update > 0) {
#if tfxMULTITHREADED
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticlePosition3d);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleSize3d);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleColor3d);
			tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleImageFrame3d);
#else
			ControlParticlePosition3d(&pm.work_queue, &work_entry);
			ControlParticleSize3d(&pm.work_queue, &work_entry);
			ControlParticleColor3d(&pm.work_queue, &work_entry);
			ControlParticleImageFrame3d(&pm.work_queue, &work_entry);
#endif
		}
	}

	void TransformEffector2d(tfxVec3 &world_rotations, tfxVec3 &local_rotations, tfxVec3 &world_position, tfxVec3 &local_position, tfxMatrix4 &matrix, tfxSpriteTransform2d &parent, bool relative_position, bool relative_angle) {

		if (relative_position) {
			world_rotations.roll = parent.rotation + local_rotations.roll;
			world_position = parent.position;
		}
		else {
			world_position = local_position;
			world_rotations.roll = local_rotations.roll;
		}

		float s = sin(world_rotations.roll);
		float c = cos(world_rotations.roll);
		matrix.Set2(c, s, -s, c);

	}

	void TransformEffector3d(tfxVec3 &world_rotations, tfxVec3 &local_rotations, tfxVec3 &world_position, tfxVec3 &local_position, tfxMatrix4 &matrix, tfxSpriteTransform3d &parent, bool relative_position, bool relative_angle) {

		if (relative_position) {
			world_rotations = parent.rotations + local_rotations;
			world_position = parent.position;
		}
		else {
			world_position = local_position;
			world_rotations = local_rotations;
		}

		tfxMatrix4 roll = mmZRotate(world_rotations.roll);
		tfxMatrix4 pitch = mmXRotate(world_rotations.pitch);
		tfxMatrix4 yaw = mmYRotate(world_rotations.yaw);

		matrix = mmTransform(yaw, pitch);
		matrix = mmTransform(matrix, roll);

	}

	tfxU32 tfxCURRENT_PROFILE_OFFSET = 0;
	const tfxU32 tfxPROFILE_COUNT = __COUNTER__;
	tfxProfile tfxPROFILE_ARRAY[tfxPROFILE_COUNT];
	tfxMemoryTrackerLog tfxMEMORY_TRACKER;
	char tfxMEMORY_CONTEXT[64];
	tfxDataTypesDictionary data_types;
	void *tfxDeferred_data_for_freeing[256];
	tfxU32 tfxDeferred_index = 0;
	HANDLE tfxThreadSemaphore;
	tfxU32 tfxNumberOfThreadsInAdditionToMain;
	tfxQueueProcessor tfxThreadQueues;
	tfxMemoryArenaManager tfxSTACK_ALLOCATOR;
	tfxMemoryArenaManager tfxMT_STACK_ALLOCATOR;

	void InitialiseTimelineFX(float percent_of_available_threads_to_use) {
		tfxSTACK_ALLOCATOR = CreateArenaManager(tfxSTACK_SIZE, 8);
		tfxMT_STACK_ALLOCATOR = CreateArenaManager(tfxMT_STACK_SIZE, 8);
		tfxU32 max_threads = tfxU32((float)std::thread::hardware_concurrency() * percent_of_available_threads_to_use);
		tfxNumberOfThreadsInAdditionToMain = tfxMin(max_threads - 1 < 0 ? 0 : max_threads - 1, std::thread::hardware_concurrency() - 1);
		tfxInitialiseThreads(&tfxThreadQueues);
		lookup_callback = LookupPrecise;
		lookup_overtime_callback = LookupPreciseOvertime;
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