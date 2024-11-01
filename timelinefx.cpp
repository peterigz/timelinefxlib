#define TFX_ALLOCATOR_IMPLEMENTATION
#include "timelinefx.h"

namespace tfx {

#ifdef _WIN32
#if defined (_MSC_VER) && (_MSC_VER >= 1400) && (defined (_M_IX86) || defined (_M_X64))
	FILE *tfx__open_file(const char *file_name, const char *mode) {
		FILE *file = NULL;
		errno_t err = fopen_s(&file, file_name, mode);
		if (err != 0 || file == NULL) {
			char errMessage[100];
			if (strerror_s(errMessage, sizeof(errMessage), err) == 0) {
				printf("strerror_s says open failed: %s\n", errMessage);
			}
			else {
				printf("Error retrieving error message\n");
			}
			return NULL;
		}
		return file;
	}
#elif defined(__GNUC__) && ((__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 8)) && \
      (defined(__i386__) || defined(__x86_64__)) || defined(__clang__)
	FILE *tfx__open_file(const char *file_name, const char *mode) {
		return fopen(file_name, mode);
	}
#endif
#else
	FILE *tfx__open_file(const char *file_name, const char *mode) {
		return fopen(file_name, mode);
	}
#endif
}

size_t tfxGetNextPower(size_t n) {
	return 1ULL << (tfx__scan_reverse(n) + 1);
}

void tfxAddHostMemoryPool(size_t size) {
	TFX_ASSERT(tfx::tfxStore->memory_pool_count < 32);    //Reached the max number of memory pools
	size_t pool_size = tfx::tfxStore->default_memory_pool_size;
	if (pool_size <= size) {
		pool_size = tfxGetNextPower(size);
	}
	TFX_PRINT_NOTICE(TFX_NOTICE_COLOR"%s: Ran out of memory, creating a new pool of size %zu. \n", TFX_NOTICE_NAME, pool_size);
	tfx::tfxStore->memory_pools[tfx::tfxStore->memory_pool_count] = (tfx_pool *)tfxALLOCATE_POOL(pool_size);
	TFX_ASSERT(tfx::tfxStore->memory_pools[tfx::tfxStore->memory_pool_count]);    //Unable to allocate more memory. Out of memory?
	tfx_AddPool(tfx::tfxMemoryAllocator, (tfx_pool *)tfx::tfxStore->memory_pools[tfx::tfxStore->memory_pool_count], pool_size);
	tfx::tfxStore->memory_pool_sizes[tfx::tfxStore->memory_pool_count] = pool_size;
	tfx::tfxStore->memory_pool_count++;
}

void *tfxAllocate(size_t size) {
	void *allocation = tfx_Allocate(tfx::tfxMemoryAllocator, size);
	if (!allocation) {
		tfxAddHostMemoryPool(size);
		allocation = tfx_Allocate(tfx::tfxMemoryAllocator, size);
		TFX_ASSERT(allocation);    //Unable to allocate even after adding a pool
	}
	return allocation;
}

void *tfxReallocate(void *memory, size_t size) {
	void *allocation = tfx_Reallocate(tfx::tfxMemoryAllocator, memory, size);
	if (!allocation) {
		tfxAddHostMemoryPool(size);
		allocation = tfx_Reallocate(tfx::tfxMemoryAllocator, memory, size);
		TFX_ASSERT(allocation);    //Unable to allocate even after adding a pool
	}
	return allocation;
}

void *tfxAllocateAligned(size_t size, size_t alignment) {
	void *allocation = tfx_AllocateAligned(tfx::tfxMemoryAllocator, size, alignment);
	if (!allocation) {
		tfxAddHostMemoryPool(size);
		allocation = tfx_AllocateAligned(tfx::tfxMemoryAllocator, size, alignment);
		TFX_ASSERT(allocation);    //Unable to allocate even after adding a pool
	}
	return allocation;
}

tfx_allocator *tfxGetAllocator() {
	return tfx::tfxMemoryAllocator;
}

tfx_bool tfx_SafeCopy(void *dst, void *src, tfx_size size) {
	tfx_header *block = tfx__block_from_allocation(dst);
	if (size > block->size) {
		return 0;
	}
	tfx_header *next_physical_block = tfx__next_physical_block(block);
	ptrdiff_t diff_check = (ptrdiff_t)((char *)dst + size) - (ptrdiff_t)next_physical_block;
	if (diff_check > 0) {
		return 0;
	}
	memcpy(dst, src, size);
	return 1;
}

tfx_bool tfx_SafeCopyBlock(void *dst_block_start, void *dst, void *src, tfx_size size) {
	tfx_header *block = tfx__block_from_allocation(dst_block_start);
	tfx_header *next_physical_block = tfx__next_physical_block(block);
	ptrdiff_t diff_check = (ptrdiff_t)((char *)dst + size) - (ptrdiff_t)next_physical_block;
	if (diff_check > 0) {
		return 0;
	}
	memcpy(dst, src, size);
	return 1;
}

tfx_bool tfx_SafeMemset(void *allocation, void *dst, int value, tfx_size size) {
	tfx_header *block = tfx__block_from_allocation(allocation);
	tfx_header *next_physical_block = tfx__next_physical_block(block);
	ptrdiff_t diff_check = (ptrdiff_t)((char *)dst + size) - (ptrdiff_t)next_physical_block;
	if (diff_check > 0) {
		return 0;
	}
	memset(dst, value, size);
	next_physical_block = tfx__next_physical_block(block);
	return 1;
}

namespace tfx {

	tfx_storage_t *GetGlobals() {
		return tfxStore;
	}

#define tfxNoise2dPermMOD12LoopUnroll(i)    \
    gi0[i] = tfx_perm_mod12[tfx_permutation_table[ii.a[i] + tfx_permutation_table[jj.a[i]]]];    \
    gi1[i] = tfx_perm_mod12[tfx_permutation_table[ii.a[i] + i1.a[i] + tfx_permutation_table[jj.a[i] + j1.a[i]]]];    \
    gi2[i] = tfx_perm_mod12[tfx_permutation_table[ii.a[i] + 1 + tfx_permutation_table[jj.a[i] + 1]]];    

#define tfxNoise3dGradientLoopUnroll(i) \
    gi0x.a[i] = gradX[gi0.a[i]];    \
    gi0y.a[i] = gradY[gi0.a[i]];    \
    gi0z.a[i] = gradZ[gi0.a[i]];    \
    gi1x.a[i] = gradX[gi1.a[i]];    \
    gi1y.a[i] = gradY[gi1.a[i]];    \
    gi1z.a[i] = gradZ[gi1.a[i]];    \
    gi2x.a[i] = gradX[gi2.a[i]];    \
    gi2y.a[i] = gradY[gi2.a[i]];    \
    gi2z.a[i] = gradZ[gi2.a[i]];    \
    gi3x.a[i] = gradX[gi3.a[i]];    \
    gi3y.a[i] = gradY[gi3.a[i]];    \
    gi3z.a[i] = gradZ[gi3.a[i]];

#define tfxNoise3dPermModLoopUnroll(i) \
    gi0.a[i] = tfx_perm_mod12[ii.a[i] + tfx_permutation_table[jj.a[i] + tfx_permutation_table[kk.a[i]]]];    \
    gi1.a[i] = tfx_perm_mod12[ii.a[i] + i1.a[i] + tfx_permutation_table[jj.a[i] + j1.a[i] + tfx_permutation_table[kk.a[i] + k1.a[i]]]];    \
    gi2.a[i] = tfx_perm_mod12[ii.a[i] + i2.a[i] + tfx_permutation_table[jj.a[i] + j2.a[i] + tfx_permutation_table[kk.a[i] + k2.a[i]]]];    \
    gi3.a[i] = tfx_perm_mod12[ii.a[i] + 1 + tfx_permutation_table[jj.a[i] + 1 + tfx_permutation_table[kk.a[i] + 1]]];    \

#ifdef tfxINTEL
	//A 2d Simd (SSE3) version of simplex noise allowing you to do 4 samples with 1 call for a speed boost
	tfx128Array tfxNoise4_2d(const tfx128 &x4, const tfx128 &y4) {
		tfxPROFILE;

		tfx128 s4 = _mm_mul_ps(_mm_add_ps(x4, y4), tfxF2_4);
		tfx128 x4_s4 = _mm_add_ps(x4, s4);
		tfx128 y4_s4 = _mm_add_ps(y4, s4);
		tfx128 i = tfxFloor128(x4_s4);
		tfx128 j = tfxFloor128(y4_s4);
		tfx128 t = _mm_add_ps(i, j);
		t = _mm_mul_ps(t, tfxG2_4);

		tfx128 X0 = _mm_sub_ps(i, t);
		tfx128 Y0 = _mm_sub_ps(j, t);
		tfx128 x0 = _mm_sub_ps(x4, X0);
		tfx128 y0 = _mm_sub_ps(y4, Y0);

		tfx128iArray i1, j1;

		i1.m = _mm_and_si128(tfxONE, _mm_castps_si128(_mm_cmpgt_ps(x0, y0)));
		j1.m = _mm_and_si128(tfxONE, _mm_castps_si128(_mm_cmpge_ps(y0, x0)));

		const tfx128 x1 = _mm_add_ps(_mm_sub_ps(x0, _mm_cvtepi32_ps(i1.m)), tfxG2_4);
		const tfx128 y1 = _mm_add_ps(_mm_sub_ps(y0, _mm_cvtepi32_ps(j1.m)), tfxG2_4);
		const tfx128 x2 = _mm_add_ps(_mm_sub_ps(x0, _mm_set1_ps(1.f)), tfxG2_4x2);
		const tfx128 y2 = _mm_add_ps(_mm_sub_ps(y0, _mm_set1_ps(1.f)), tfxG2_4x2);

		tfx128iArray ii, jj;
		ii.m = _mm_and_si128(_mm_cvttps_epi32(i), tfxFF);
		jj.m = _mm_and_si128(_mm_cvttps_epi32(j), tfxFF);

		int gi0[4], gi1[4], gi2[4];

		tfxNoise2dPermMOD12LoopUnroll(0);
		tfxNoise2dPermMOD12LoopUnroll(1);
		tfxNoise2dPermMOD12LoopUnroll(2);
		tfxNoise2dPermMOD12LoopUnroll(3);

		tfx128 n0, n1, n2;
		tfx128 gx0, gy0, gx1, gy1, gx2, gy2;
		gx0 = _mm_set_ps(gradX[gi0[3]], gradX[gi0[2]], gradX[gi0[1]], gradX[gi0[0]]);
		gy0 = _mm_set_ps(gradY[gi0[3]], gradY[gi0[2]], gradY[gi0[1]], gradY[gi0[0]]);
		gx1 = _mm_set_ps(gradX[gi1[3]], gradX[gi1[2]], gradX[gi1[1]], gradX[gi1[0]]);
		gy1 = _mm_set_ps(gradY[gi1[3]], gradY[gi1[2]], gradY[gi1[1]], gradY[gi1[0]]);
		gx2 = _mm_set_ps(gradX[gi2[3]], gradX[gi2[2]], gradX[gi2[1]], gradX[gi2[0]]);
		gy2 = _mm_set_ps(gradY[gi2[3]], gradY[gi2[2]], gradY[gi2[1]], gradY[gi2[0]]);

		tfx128 t0 = _mm_sub_ps(_mm_sub_ps(_mm_set1_ps(0.5f), _mm_mul_ps(x0, x0)), _mm_mul_ps(y0, y0));
		tfx128 t02 = _mm_mul_ps(t0, t0);
		n0 = _mm_and_ps(_mm_mul_ps(_mm_mul_ps(t02, t02), tfx__dot128_xy(&gx0, &gy0, &x0, &y0)), _mm_cmpge_ps(t0, _mm_setzero_ps()));

		tfx128 t1 = _mm_sub_ps(_mm_sub_ps(_mm_set1_ps(0.5f), _mm_mul_ps(x1, x1)), _mm_mul_ps(y1, y1));
		tfx128 t12 = _mm_mul_ps(t1, t1);
		n1 = _mm_and_ps(_mm_mul_ps(_mm_mul_ps(t12, t12), tfx__dot128_xy(&gx1, &gy1, &x1, &y1)), _mm_cmpge_ps(t1, _mm_setzero_ps()));

		tfx128 t2 = _mm_sub_ps(_mm_sub_ps(_mm_set1_ps(0.5f), _mm_mul_ps(x2, x2)), _mm_mul_ps(y2, y2));
		tfx128 t22 = _mm_mul_ps(t2, t2);
		n2 = _mm_and_ps(_mm_mul_ps(_mm_mul_ps(t22, t22), tfx__dot128_xy(&gx2, &gy2, &x2, &y2)), _mm_cmpge_ps(t2, _mm_setzero_ps()));

		tfx128Array result;
		result.m = _mm_mul_ps(_mm_set1_ps(45.23065f), _mm_add_ps(n0, _mm_add_ps(n1, n2)));
		return result;
	}

	//A 3d Simd (SSE3) version of simplex noise allowing you to do 4 samples with 1 call for a speed boost
	tfx128Array tfxNoise4_3d(const tfx128 &x4, const tfx128 &y4, const tfx128 &z4) {
		tfxPROFILE;
		// Skewing/Unskewing factors for 3D

		// Skew the input space to determine which simplex cell we're in
		//float s = (v1.x + v1.y + v1.z) * F3; // Very nice and simple skew factor for 3D
		tfx128 s4 = _mm_mul_ps(_mm_add_ps(x4, _mm_add_ps(y4, z4)), tfxF3_4);
		tfx128 x4_s4 = _mm_add_ps(x4, s4);
		tfx128 y4_s4 = _mm_add_ps(y4, s4);
		tfx128 z4_s4 = _mm_add_ps(z4, s4);
		tfx128 i = tfxFloor128(x4_s4);
		tfx128 j = tfxFloor128(y4_s4);
		tfx128 k = tfxFloor128(z4_s4);
		tfx128 t = _mm_add_ps(i, j);
		t = _mm_add_ps(t, k);
		t = _mm_mul_ps(t, tfxG3_4);

		tfx128 X0 = _mm_sub_ps(i, t); // Unskew the cell origin back to (v1.x,v1.y,v1.z) space
		tfx128 Y0 = _mm_sub_ps(j, t);
		tfx128 Z0 = _mm_sub_ps(k, t);
		tfx128 x0 = _mm_sub_ps(x4, X0); // The v1.x,v1.y,v1.z distances from the cell origin
		tfx128 y0 = _mm_sub_ps(y4, Y0);
		tfx128 z0 = _mm_sub_ps(z4, Z0);

		// For the 3D case, the simplex shape is a slightly irregular tetrahedron.
		// Determine which simplex we are in.
		tfx128iArray i1, i2, j1, j2, k1, k2;

		i1.m = _mm_and_si128(tfxONE, _mm_and_si128(_mm_castps_si128(_mm_cmpge_ps(x0, y0)), _mm_castps_si128(_mm_cmpge_ps(x0, z0))));
		j1.m = _mm_and_si128(tfxONE, _mm_and_si128(_mm_castps_si128(_mm_cmpgt_ps(y0, x0)), _mm_castps_si128(_mm_cmpge_ps(y0, z0))));
		k1.m = _mm_and_si128(tfxONE, _mm_and_si128(_mm_castps_si128(_mm_cmpgt_ps(z0, x0)), _mm_castps_si128(_mm_cmpgt_ps(z0, y0))));

		//for i2
		tfx128i yx_xz = _mm_and_si128(_mm_castps_si128(_mm_cmpge_ps(x0, y0)), _mm_castps_si128(_mm_cmplt_ps(x0, z0)));
		tfx128i zx_xy = _mm_and_si128(_mm_castps_si128(_mm_cmpge_ps(x0, z0)), _mm_castps_si128(_mm_cmplt_ps(x0, y0)));

		//for j2
		tfx128i xy_yz = _mm_and_si128(_mm_castps_si128(_mm_cmplt_ps(x0, y0)), _mm_castps_si128(_mm_cmplt_ps(y0, z0)));
		tfx128i zy_yx = _mm_and_si128(_mm_castps_si128(_mm_cmpge_ps(y0, z0)), _mm_castps_si128(_mm_cmpge_ps(x0, y0)));

		//for k2
		tfx128i yz_zx = _mm_and_si128(_mm_castps_si128(_mm_cmplt_ps(y0, z0)), _mm_castps_si128(_mm_cmpge_ps(x0, z0)));
		tfx128i xz_zy = _mm_and_si128(_mm_castps_si128(_mm_cmplt_ps(x0, z0)), _mm_castps_si128(_mm_cmpge_ps(y0, z0)));

		i2.m = _mm_and_si128(tfxONE, _mm_or_si128(i1.m, _mm_or_si128(yx_xz, zx_xy)));
		j2.m = _mm_and_si128(tfxONE, _mm_or_si128(j1.m, _mm_or_si128(xy_yz, zy_yx)));
		k2.m = _mm_and_si128(tfxONE, _mm_or_si128(k1.m, _mm_or_si128(yz_zx, xz_zy)));

		tfx128 x1 = _mm_add_ps(_mm_sub_ps(x0, _mm_cvtepi32_ps(i1.m)), tfxG3_4);
		tfx128 y1 = _mm_add_ps(_mm_sub_ps(y0, _mm_cvtepi32_ps(j1.m)), tfxG3_4);
		tfx128 z1 = _mm_add_ps(_mm_sub_ps(z0, _mm_cvtepi32_ps(k1.m)), tfxG3_4);
		tfx128 x2 = _mm_add_ps(_mm_sub_ps(x0, _mm_cvtepi32_ps(i2.m)), tfxG32_4);
		tfx128 y2 = _mm_add_ps(_mm_sub_ps(y0, _mm_cvtepi32_ps(j2.m)), tfxG32_4);
		tfx128 z2 = _mm_add_ps(_mm_sub_ps(z0, _mm_cvtepi32_ps(k2.m)), tfxG32_4);
		tfx128 x3 = _mm_add_ps(_mm_sub_ps(x0, tfxONEF), tfxG33_4);
		tfx128 y3 = _mm_add_ps(_mm_sub_ps(y0, tfxONEF), tfxG33_4);
		tfx128 z3 = _mm_add_ps(_mm_sub_ps(z0, tfxONEF), tfxG33_4);

		// Work out the hashed gradient indices of the four simplex corners
		tfx128iArray ii;
		ii.m = _mm_and_si128(_mm_cvttps_epi32(i), tfxFF);
		tfx128iArray jj;
		jj.m = _mm_and_si128(_mm_cvttps_epi32(j), tfxFF);
		tfx128iArray kk;
		kk.m = _mm_and_si128(_mm_cvttps_epi32(k), tfxFF);
		tfx128iArray gi0, gi1, gi2, gi3;

		tfxNoise3dPermModLoopUnroll(0);
		tfxNoise3dPermModLoopUnroll(1);
		tfxNoise3dPermModLoopUnroll(2);
		tfxNoise3dPermModLoopUnroll(3);

		tfx128 t0 = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(tfxPSIX, _mm_mul_ps(x0, x0)), _mm_mul_ps(y0, y0)), _mm_mul_ps(z0, z0));
		tfx128 t1 = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(tfxPSIX, _mm_mul_ps(x1, x1)), _mm_mul_ps(y1, y1)), _mm_mul_ps(z1, z1));
		tfx128 t2 = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(tfxPSIX, _mm_mul_ps(x2, x2)), _mm_mul_ps(y2, y2)), _mm_mul_ps(z2, z2));
		tfx128 t3 = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(tfxPSIX, _mm_mul_ps(x3, x3)), _mm_mul_ps(y3, y3)), _mm_mul_ps(z3, z3));

		tfx128 t0q = _mm_mul_ps(t0, t0);
		t0q = _mm_mul_ps(t0q, t0q);
		tfx128 t1q = _mm_mul_ps(t1, t1);
		t1q = _mm_mul_ps(t1q, t1q);
		tfx128 t2q = _mm_mul_ps(t2, t2);
		t2q = _mm_mul_ps(t2q, t2q);
		tfx128 t3q = _mm_mul_ps(t3, t3);
		t3q = _mm_mul_ps(t3q, t3q);

		tfx128Array gi0x, gi0y, gi0z, gi1x, gi1y, gi1z, gi2x, gi2y, gi2z, gi3x, gi3y, gi3z;

		tfxNoise3dGradientLoopUnroll(0)
		tfxNoise3dGradientLoopUnroll(1)
		tfxNoise3dGradientLoopUnroll(2)
		tfxNoise3dGradientLoopUnroll(3)

		tfx128 n0 = _mm_mul_ps(t0q, tfx__dot128_xyz(&gi0x.m, &gi0y.m, &gi0z.m, &x0, &y0, &z0));
		tfx128 n1 = _mm_mul_ps(t1q, tfx__dot128_xyz(&gi1x.m, &gi1y.m, &gi1z.m, &x1, &y1, &z1));
		tfx128 n2 = _mm_mul_ps(t2q, tfx__dot128_xyz(&gi2x.m, &gi2y.m, &gi2z.m, &x2, &y2, &z2));
		tfx128 n3 = _mm_mul_ps(t3q, tfx__dot128_xyz(&gi3x.m, &gi3y.m, &gi3z.m, &x3, &y3, &z3));

		tfx128 cond;

		cond = _mm_cmplt_ps(t0, tfxZERO);
		n0 = _mm_or_ps(_mm_andnot_ps(cond, n0), _mm_and_ps(cond, tfxZERO));
		cond = _mm_cmplt_ps(t1, tfxZERO);
		n1 = _mm_or_ps(_mm_andnot_ps(cond, n1), _mm_and_ps(cond, tfxZERO));
		cond = _mm_cmplt_ps(t2, tfxZERO);
		n2 = _mm_or_ps(_mm_andnot_ps(cond, n2), _mm_and_ps(cond, tfxZERO));
		cond = _mm_cmplt_ps(t3, tfxZERO);
		n3 = _mm_or_ps(_mm_andnot_ps(cond, n3), _mm_and_ps(cond, tfxZERO));

		tfx128Array result;
		result.m = _mm_mul_ps(tfxTHIRTYTWO, _mm_add_ps(n0, _mm_add_ps(n1, _mm_add_ps(n2, n3))));
		return result;
	}

	tfx128 tfx__dot128_xyz(const tfx128 *x1, const tfx128 *y1, const tfx128 *z1, const tfx128 *x2, const tfx128 *y2, const tfx128 *z2)
	{
		tfx128 xx = _mm_mul_ps(*x1, *x2);
		tfx128 yy = _mm_mul_ps(*y1, *y2);
		tfx128 zz = _mm_mul_ps(*z1, *z2);
		return _mm_add_ps(xx, _mm_add_ps(yy, zz));
	}

	tfx128 tfx__dot128_xy(const tfx128 *x1, const tfx128 *y1, const tfx128 *x2, const tfx128 *y2)
	{
		tfx128 xx = _mm_mul_ps(*x1, *x2);
		tfx128 yy = _mm_mul_ps(*y1, *y2);
		return _mm_add_ps(xx, yy);
	}

	tfx_mat4_t tfx__transform_matrix4(const tfx_mat4_t *in, const tfx_mat4_t *m) {
		tfx_mat4_t res = tfx__create_matrix4(0.f);

		tfx128 in_row[4];
		in_row[0] = _mm_load_ps(&in->v[0].x);
		in_row[1] = _mm_load_ps(&in->v[1].x);
		in_row[2] = _mm_load_ps(&in->v[2].x);
		in_row[3] = _mm_load_ps(&in->v[3].x);

		tfx128 m_row1 = _mm_set_ps(m->v[3].x, m->v[2].x, m->v[1].x, m->v[0].x);
		tfx128 m_row2 = _mm_set_ps(m->v[3].y, m->v[2].y, m->v[1].y, m->v[0].y);
		tfx128 m_row3 = _mm_set_ps(m->v[3].z, m->v[2].z, m->v[1].z, m->v[0].z);
		tfx128 m_row4 = _mm_set_ps(m->v[3].w, m->v[2].w, m->v[1].w, m->v[0].w);

		for (int r = 0; r <= 3; ++r)
		{

			tfx128 row1result = _mm_mul_ps(in_row[r], m_row1);
			tfx128 row2result = _mm_mul_ps(in_row[r], m_row2);
			tfx128 row3result = _mm_mul_ps(in_row[r], m_row3);
			tfx128 row4result = _mm_mul_ps(in_row[r], m_row4);

			float tmp[4];
			_mm_store_ps(tmp, row1result);
			res.v[r].x = tmp[0] + tmp[1] + tmp[2] + tmp[3];
			_mm_store_ps(tmp, row2result);
			res.v[r].y = tmp[0] + tmp[1] + tmp[2] + tmp[3];
			_mm_store_ps(tmp, row3result);
			res.v[r].z = tmp[0] + tmp[1] + tmp[2] + tmp[3];
			_mm_store_ps(tmp, row4result);
			res.v[r].w = tmp[0] + tmp[1] + tmp[2] + tmp[3];

		}
		return res;
	}

	tfx_vec4_t tfx__transform_matrix4_vec4(const tfx_mat4_t *mat, const tfx_vec4_t vec) {
		tfx_vec4_t v;

		tfx128 v4 = _mm_set_ps(vec.w, vec.z, vec.y, vec.x);

		tfx__readbarrier;

		tfx128 mrow1 = _mm_load_ps(&mat->v[0].c0);
		tfx128 mrow2 = _mm_load_ps(&mat->v[1].c0);
		tfx128 mrow3 = _mm_load_ps(&mat->v[2].c0);
		tfx128 mrow4 = _mm_load_ps(&mat->v[3].c0);

		tfx__readbarrier;

		tfx128 row1result = _mm_mul_ps(v4, mrow1);
		tfx128 row2result = _mm_mul_ps(v4, mrow2);
		tfx128 row3result = _mm_mul_ps(v4, mrow3);
		tfx128 row4result = _mm_mul_ps(v4, mrow4);

		float tmp[4];
		_mm_store_ps(tmp, row1result);
		v.x = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		_mm_store_ps(tmp, row2result);
		v.y = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		_mm_store_ps(tmp, row3result);
		v.z = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		_mm_store_ps(tmp, row4result);
		v.w = tmp[0] + tmp[1] + tmp[2] + tmp[3];

		return v;
	}

	tfx_vec4_t tfx_InterpolateVec4(float tween, tfx_vec4_t *from, tfx_vec4_t *to) {
		tfx128 l4 = _mm_set_ps1(tween);
		tfx128 l4minus1 = _mm_set_ps1(1.f - tween);
		tfx128 f4 = _mm_set_ps(from->x, from->y, from->z, from->w);
		tfx128 t4 = _mm_set_ps(to->x, to->y, to->z, to->w);
		tfx128 from_lerp = _mm_mul_ps(f4, l4);
		tfx128 to_lerp = _mm_mul_ps(f4, l4minus1);
		tfx128 result = _mm_add_ps(from_lerp, to_lerp);
		tfx_vec4_t vec;
		_mm_store_ps(&vec.x, result);
		return vec;
	}

#elif defined(tfxARM)

	tfx128Array tfxNoise4_2d(const tfx128 &x4, const tfx128 &y4) {
		tfxPROFILE;

		tfx128 s4 = vmulq_f32(vaddq_f32(x4, y4), tfxF2_4);
		tfx128 x4_s4 = vaddq_f32(x4, s4);
		tfx128 y4_s4 = vaddq_f32(y4, s4);
		tfx128 i = tfxFloor128(x4_s4);
		tfx128 j = tfxFloor128(y4_s4);
		tfx128 t = vmulq_f32(vaddq_f32(i, j), tfxG2_4);

		tfx128 X0 = vsubq_f32(i, t);
		tfx128 Y0 = vsubq_f32(j, t);
		tfx128 x0 = vsubq_f32(x4, X0);
		tfx128 y0 = vsubq_f32(y4, Y0);

		tfx128iArray i1, j1;

		i1.m = vandq_s32(tfxONE, vreinterpretq_s32_f32(vcgtq_f32(x0, y0)));
		j1.m = vandq_s32(tfxONE, vreinterpretq_s32_f32(vcgeq_f32(y0, x0)));

		const tfx128 x1 = vaddq_f32(vsubq_f32(x0, vcvtq_f32_s32(i1.m)), tfxG2_4);
		const tfx128 y1 = vaddq_f32(vsubq_f32(y0, vcvtq_f32_s32(j1.m)), tfxG2_4);
		const tfx128 x2 = vaddq_f32(vsubq_f32(x0, vdupq_n_f32(1.f)), tfxG2_4x2);
		const tfx128 y2 = vaddq_f32(vsubq_f32(y0, vdupq_n_f32(1.f)), tfxG2_4x2);

		tfx128iArray ii, jj;
		ii.m = vandq_s32(vcvtq_s32_f32(i), tfxFF);
		jj.m = vandq_s32(vcvtq_s32_f32(j), tfxFF);

		int gi0[4], gi1[4], gi2[4];

		tfxNoise2dPermMOD12LoopUnroll(0);
		tfxNoise2dPermMOD12LoopUnroll(1);
		tfxNoise2dPermMOD12LoopUnroll(2);
		tfxNoise2dPermMOD12LoopUnroll(3);

		tfx128 n0, n1, n2;
		tfx128 gx0, gy0, gx1, gy1, gx2, gy2;
		gx0 = tfx128Set(gradX[gi0[3]], gradX[gi0[2]], gradX[gi0[1]], gradX[gi0[0]]);
		gy0 = tfx128Set(gradY[gi0[3]], gradY[gi0[2]], gradY[gi0[1]], gradY[gi0[0]]);
		gx1 = tfx128Set(gradX[gi1[3]], gradX[gi1[2]], gradX[gi1[1]], gradX[gi1[0]]);
		gy1 = tfx128Set(gradY[gi1[3]], gradY[gi1[2]], gradY[gi1[1]], gradY[gi1[0]]);
		gx2 = tfx128Set(gradX[gi2[3]], gradX[gi2[2]], gradX[gi2[1]], gradX[gi2[0]]);
		gy2 = tfx128Set(gradY[gi2[3]], gradY[gi2[2]], gradY[gi2[1]], gradY[gi2[0]]);

		tfx128 t0 = vsubq_f32(vsubq_f32(vdupq_n_f32(0.5f), vmulq_f32(x0, x0)), vmulq_f32(y0, y0));
		tfx128 t02 = vmulq_f32(t0, t0);
		n0 = tfxSIMD_AND(vmulq_f32(vmulq_f32(t02, t02), tfx__dot128_xy(&gx0, &gy0, &x0, &y0)), vcgeq_f32(t0, vdupq_n_f32(0.f)));

		tfx128 t1 = vsubq_f32(vsubq_f32(vdupq_n_f32(0.5f), vmulq_f32(x1, x1)), vmulq_f32(y1, y1));
		tfx128 t12 = vmulq_f32(t1, t1);
		n1 = tfxSIMD_AND(vmulq_f32(vmulq_f32(t12, t12), tfx__dot128_xy(&gx1, &gy1, &x1, &y1)), vcgeq_f32(t1, vdupq_n_f32(0.f)));

		tfx128 t2 = vsubq_f32(vsubq_f32(vdupq_n_f32(0.5f), vmulq_f32(x2, x2)), vmulq_f32(y2, y2));
		tfx128 t22 = vmulq_f32(t2, t2);
		n2 = tfxSIMD_AND(vmulq_f32(vmulq_f32(t22, t22), tfx__dot128_xy(&gx2, &gy2, &x2, &y2)), vcgeq_f32(t2, vdupq_n_f32(0.f)));

		tfx128Array result;
		result.m = vmulq_f32(vdupq_n_f32(45.23065f), vaddq_f32(n0, vaddq_f32(n1, n2)));
		return result;
	}

	tfx128Array tfxNoise4_3d(const tfx128 &x4, const tfx128 &y4, const tfx128 &z4) {
		tfxPROFILE;

		// Skewing/Unskewing factors for 3D
		tfx128 s4 = vmulq_f32(vaddq_f32(x4, vaddq_f32(y4, z4)), tfxF3_4);
		tfx128 x4_s4 = vaddq_f32(x4, s4);
		tfx128 y4_s4 = vaddq_f32(y4, s4);
		tfx128 z4_s4 = vaddq_f32(z4, s4);
		tfx128 i = tfxFloor128(x4_s4);
		tfx128 j = tfxFloor128(y4_s4);
		tfx128 k = tfxFloor128(z4_s4);
		tfx128 t = vmulq_f32(vaddq_f32(vaddq_f32(i, j), k), tfxG3_4);

		tfx128 X0 = vsubq_f32(i, t); // Unskew the cell origin back to (v1.x,v1.y,v1.z) space
		tfx128 Y0 = vsubq_f32(j, t);
		tfx128 Z0 = vsubq_f32(k, t);
		tfx128 x0 = vsubq_f32(x4, X0); // The v1.x,v1.y,v1.z distances from the cell origin
		tfx128 y0 = vsubq_f32(y4, Y0);
		tfx128 z0 = vsubq_f32(z4, Z0);

		// For the 3D case, the simplex shape is a slightly irregular tetrahedron.
		// Determine which simplex we are in.
		tfx128iArray i1, i2, j1, j2, k1, k2;

		i1.m = vandq_s32(tfxONE, vandq_s32(vreinterpretq_s32_f32(vcgeq_f32(x0, y0)), vreinterpretq_s32_f32(vcgeq_f32(x0, z0))));
		j1.m = vandq_s32(tfxONE, vandq_s32(vreinterpretq_s32_f32(vcgtq_f32(y0, x0)), vreinterpretq_s32_f32(vcgeq_f32(y0, z0))));
		k1.m = vandq_s32(tfxONE, vandq_s32(vreinterpretq_s32_f32(vcgtq_f32(z0, x0)), vreinterpretq_s32_f32(vcgtq_f32(z0, y0))));

		//for i2
		tfx128i yx_xz = vandq_s32(vreinterpretq_s32_f32(vcgeq_f32(x0, y0)), vreinterpretq_s32_f32(vcltq_f32(x0, z0)));
		tfx128i zx_xy = vandq_s32(vreinterpretq_s32_f32(vcgeq_f32(x0, z0)), vreinterpretq_s32_f32(vcltq_f32(x0, y0)));

		//for j2
		tfx128i xy_yz = vandq_s32(vreinterpretq_s32_f32(vcltq_f32(x0, y0)), vreinterpretq_s32_f32(vcltq_f32(y0, z0)));
		tfx128i zy_yx = vandq_s32(vreinterpretq_s32_f32(vcgeq_f32(y0, z0)), vreinterpretq_s32_f32(vcgeq_f32(x0, y0)));

		//for k2
		tfx128i yz_zx = vandq_s32(vreinterpretq_s32_f32(vcltq_f32(y0, z0)), vreinterpretq_s32_f32(vcgeq_f32(x0, z0)));
		tfx128i xz_zy = vandq_s32(vreinterpretq_s32_f32(vcltq_f32(x0, z0)), vreinterpretq_s32_f32(vcgeq_f32(y0, z0)));

		i2.m = vandq_s32(tfxONE, vorrq_s32(i1.m, vorrq_s32(yx_xz, zx_xy)));
		j2.m = vandq_s32(tfxONE, vorrq_s32(j1.m, vorrq_s32(xy_yz, zy_yx)));
		k2.m = vandq_s32(tfxONE, vorrq_s32(k1.m, vorrq_s32(yz_zx, xz_zy)));

		tfx128 x1 = vaddq_f32(vsubq_f32(x0, vcvtq_f32_s32(i1.m)), tfxG3_4);
		tfx128 y1 = vaddq_f32(vsubq_f32(y0, vcvtq_f32_s32(j1.m)), tfxG3_4);
		tfx128 z1 = vaddq_f32(vsubq_f32(z0, vcvtq_f32_s32(k1.m)), tfxG3_4);
		tfx128 x2 = vaddq_f32(vsubq_f32(x0, vcvtq_f32_s32(i2.m)), tfxG32_4);
		tfx128 y2 = vaddq_f32(vsubq_f32(y0, vcvtq_f32_s32(j2.m)), tfxG32_4);
		tfx128 z2 = vaddq_f32(vsubq_f32(z0, vcvtq_f32_s32(k2.m)), tfxG32_4);
		tfx128 x3 = vaddq_f32(vsubq_f32(x0, tfxONEF), tfxG33_4);
		tfx128 y3 = vaddq_f32(vsubq_f32(y0, tfxONEF), tfxG33_4);
		tfx128 z3 = vaddq_f32(vsubq_f32(z0, tfxONEF), tfxG33_4);

		// Work out the hashed gradient indices of the four simplex corners
		tfx128iArray ii, jj, kk;
		ii.m = vandq_s32(vcvtq_s32_f32(i), tfxFF);
		jj.m = vandq_s32(vcvtq_s32_f32(j), tfxFF);
		kk.m = vandq_s32(vcvtq_s32_f32(k), tfxFF);
		tfx128iArray gi0, gi1, gi2, gi3;

		tfxNoise3dPermModLoopUnroll(0);
		tfxNoise3dPermModLoopUnroll(1);
		tfxNoise3dPermModLoopUnroll(2);
		tfxNoise3dPermModLoopUnroll(3);

		tfx128 t0 = vsubq_f32(vsubq_f32(vsubq_f32(tfxPSIX, vmulq_f32(x0, x0)), vmulq_f32(y0, y0)), vmulq_f32(z0, z0));
		tfx128 t1 = vsubq_f32(vsubq_f32(vsubq_f32(tfxPSIX, vmulq_f32(x1, x1)), vmulq_f32(y1, y1)), vmulq_f32(z1, z1));
		tfx128 t2 = vsubq_f32(vsubq_f32(vsubq_f32(tfxPSIX, vmulq_f32(x2, x2)), vmulq_f32(y2, y2)), vmulq_f32(z2, z2));
		tfx128 t3 = vsubq_f32(vsubq_f32(vsubq_f32(tfxPSIX, vmulq_f32(x3, x3)), vmulq_f32(y3, y3)), vmulq_f32(z3, z3));

		tfx128 t0q = vmulq_f32(t0, t0);
		t0q = vmulq_f32(t0q, t0q);
		tfx128 t1q = vmulq_f32(t1, t1);
		t1q = vmulq_f32(t1q, t1q);
		tfx128 t2q = vmulq_f32(t2, t2);
		t2q = vmulq_f32(t2q, t2q);
		tfx128 t3q = vmulq_f32(t3, t3);
		t3q = vmulq_f32(t3q, t3q);

		tfx128Array gi0x, gi0y, gi0z, gi1x, gi1y, gi1z, gi2x, gi2y, gi2z, gi3x, gi3y, gi3z;

		tfxNoise3dGradientLoopUnroll(0)
			tfxNoise3dGradientLoopUnroll(1)
			tfxNoise3dGradientLoopUnroll(2)
			tfxNoise3dGradientLoopUnroll(3)

			tfx128 n0 = vmulq_f32(t0q, tfx__dot128_xyz(&gi0x.m, &gi0y.m, &gi0z.m, &x0, &y0, &z0));
		tfx128 n1 = vmulq_f32(t1q, tfx__dot128_xyz(&gi1x.m, &gi1y.m, &gi1z.m, &x1, &y1, &z1));
		tfx128 n2 = vmulq_f32(t2q, tfx__dot128_xyz(&gi2x.m, &gi2y.m, &gi2z.m, &x2, &y2, &z2));
		tfx128 n3 = vmulq_f32(t3q, tfx__dot128_xyz(&gi3x.m, &gi3y.m, &gi3z.m, &x3, &y3, &z3));

		tfx128 cond;

		cond = vcltq_f32(t0, tfxZERO);
		n0 = vbslq_f32(cond, tfxZERO, n0);
		cond = vcltq_f32(t1, tfxZERO);
		n1 = vbslq_f32(cond, tfxZERO, n1);
		cond = vcltq_f32(t2, tfxZERO);
		n2 = vbslq_f32(cond, tfxZERO, n2);
		cond = vcltq_f32(t3, tfxZERO);
		n3 = vbslq_f32(cond, tfxZERO, n3);

		tfx128Array result;
		result.m = vmulq_f32(tfxTHIRTYTWO, vaddq_f32(n0, vaddq_f32(n1, vaddq_f32(n2, n3))));
		return result;
	}

	tfx128 tfx__dot128_xyz(const tfx128 *x1, const tfx128 *y1, const tfx128 *z1, const tfx128 *x2, const tfx128 *y2, const tfx128 *z2)
	{
		tfx128 xx = vmulq_f32(*x1, *x2);
		tfx128 yy = vmulq_f32(*y1, *y2);
		tfx128 zz = vmulq_f32(*z1, *z2);
		return vaddq_f32(xx, vaddq_f32(yy, zz));
	}

	tfx128 tfx__dot128_xy(const tfx128 *x1, const tfx128 *y1, const tfx128 *x2, const tfx128 *y2)
	{
		tfx128 xx = vmulq_f32(*x1, *x2);
		tfx128 yy = vmulq_f32(*y1, *y2);
		return vaddq_f32(xx, yy);
	}

	tfx_mat4_t tfx__transform_matrix4(const tfx_mat4_t *in, const tfx_mat4_t *m) {
		tfx_mat4_t res = tfx__create_matrix4(0.f);

		tfx128 in_row[4];
		in_row[0] = vld1q_f32(&in->v[0].x);
		in_row[1] = vld1q_f32(&in->v[1].x);
		in_row[2] = vld1q_f32(&in->v[2].x);
		in_row[3] = vld1q_f32(&in->v[3].x);

		tfx128 m_row1 = tfx128Set(m->v[3].x, m->v[2].x, m->v[1].x, m->v[0].x);
		tfx128 m_row2 = tfx128Set(m->v[3].y, m->v[2].y, m->v[1].y, m->v[0].y);
		tfx128 m_row3 = tfx128Set(m->v[3].z, m->v[2].z, m->v[1].z, m->v[0].z);
		tfx128 m_row4 = tfx128Set(m->v[3].w, m->v[2].w, m->v[1].w, m->v[0].w);

		for (int r = 0; r <= 3; ++r)
		{
			tfx128 row1result = vmulq_f32(in_row[r], m_row1);
			tfx128 row2result = vmulq_f32(in_row[r], m_row2);
			tfx128 row3result = vmulq_f32(in_row[r], m_row3);
			tfx128 row4result = vmulq_f32(in_row[r], m_row4);

			float tmp[4];
			vst1q_f32(tmp, row1result);
			res.v[r].x = tmp[0] + tmp[1] + tmp[2] + tmp[3];
			vst1q_f32(tmp, row2result);
			res.v[r].y = tmp[0] + tmp[1] + tmp[2] + tmp[3];
			vst1q_f32(tmp, row3result);
			res.v[r].z = tmp[0] + tmp[1] + tmp[2] + tmp[3];
			vst1q_f32(tmp, row4result);
			res.v[r].w = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		}
		return res;
	}

	tfx_vec4_t tfx__transform_matrix4_vec4(const tfx_mat4_t *mat, const tfx_vec4_t vec) {
		tfx_vec4_t v;

		tfx128 v4 = vld1q_f32(&vec.x);

		tfx128 mrow1 = vld1q_f32(&mat->v[0].x);
		tfx128 mrow2 = vld1q_f32(&mat->v[1].x);
		tfx128 mrow3 = vld1q_f32(&mat->v[2].x);
		tfx128 mrow4 = vld1q_f32(&mat->v[3].x);

		tfx128 row1result = vmulq_f32(v4, mrow1);
		tfx128 row2result = vmulq_f32(v4, mrow2);
		tfx128 row3result = vmulq_f32(v4, mrow3);
		tfx128 row4result = vmulq_f32(v4, mrow4);

		float tmp[4];
		vst1q_f32(tmp, row1result);
		v.x = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		vst1q_f32(tmp, row2result);
		v.y = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		vst1q_f32(tmp, row3result);
		v.z = tmp[0] + tmp[1] + tmp[2] + tmp[3];
		vst1q_f32(tmp, row4result);
		v.w = tmp[0] + tmp[1] + tmp[2] + tmp[3];

		return v;
	}

	tfx_vec4_t tfx_InterpolateVec4(float tween, tfx_vec4_t *from, tfx_vec4_t *to) {
		tfx128 l4 = vdupq_n_f32(tween);
		tfx128 l4minus1 = vdupq_n_f32(1.f - tween);
		tfx128 f4 = tfx128Set(from->x, from->y, from->z, from->w);
		tfx128 t4 = tfx128Set(to->x, to->y, to->z, to->w);
		tfx128 from_lerp = vmulq_f32(f4, l4);
		tfx128 to_lerp = vmulq_f32(f4, l4minus1);
		tfx128 result = vaddq_f32(from_lerp, to_lerp);
		tfx_vec4_t vec;
		vst1q_f32(&vec.x, result);
		return vec;
	}

#endif

tfx_random_t NewRandom(tfxU32 seed) {
	tfx_random_t random;
	memset(random.seeds, 0, sizeof(tfxU64) * 2);
	RandomReSeed(&random, seed);
	return random;
}

void AdvanceRandom(tfx_random_t *random) {
	tfxU64 s1 = random->seeds[0];
	tfxU64 s0 = random->seeds[1];
	random->seeds[0] = s0;
	s1 ^= s1 << 23; // a
	random->seeds[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5); // b, c
}

void RandomReSeed(tfx_random_t *random) {
	random->seeds[0] = tfx_Millisecs(); random->seeds[1] = (tfxU64)tfx_Millisecs() * 2;
	AdvanceRandom(random);
}

void RandomReSeed(tfx_random_t *random, tfxU64 seed1, tfxU64 seed2) {
	random->seeds[0] = seed1;
	random->seeds[1] = seed2;
	AdvanceRandom(random);
}

void RandomReSeed(tfx_random_t *random, tfxU64 seed) {
	random->seeds[0] = seed;
	random->seeds[1] = seed * 2;
	AdvanceRandom(random);
}

float GenerateRandom(tfx_random_t *random) {
	tfxU64 s1 = random->seeds[0];
	tfxU64 s0 = random->seeds[1];
	tfxU64 result = s0 + s1;
	random->seeds[0] = s0;
	s1 ^= s1 << 23; // a
	random->seeds[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5); // b, c
	return float((double)result / tfxTWO64f);
}

float RandomRange(tfx_random_t *random, float max) {
	return GenerateRandom(random) * max;
};

float RandomRange(tfx_random_t *random, float from, float to) {
	float a = GenerateRandom(random);
	float range = to - from;
	return to - range * a;
};

int RandomRange(tfx_random_t *random, int from, int to) {
	float a = (to - from + 1) * GenerateRandom(random) + from;
	return (int)a;
};

tfxU32 RandomRange(tfx_random_t *random, tfxU32 max) {
	float g = GenerateRandom(random);
	float a = g * (float)max;
	return tfxU32(a);
};

int RandomRange(tfx_random_t *random, int max) {
	float g = GenerateRandom(random);
	float a = g * (float)max;
	return tfxU32(a);
};

void AlterRandomSeed(tfx_random_t *random, tfxU64 amount) {
	random->seeds[0] *= amount;
	random->seeds[1] += amount;
};

void AlterRandomSeed(tfx_random_t *random, tfxU32 amount) {
	random->seeds[0] *= amount;
	random->seeds[1] += amount;
};

tfx_rgba8_t tfx__convert_float_color(float color_array[4]) {
	tfx_rgba8_t color;
	color.r = (char)tfx__Min(255, (int)(color_array[0] * 255.f));
	color.g = (char)tfx__Min(255, (int)(color_array[1] * 255.f));
	color.b = (char)tfx__Min(255, (int)(color_array[2] * 255.f));
	color.a = (char)tfx__Min(255, (int)(color_array[3] * 255.f));
	return color;
}

tfx_vector_t<tfx_vec3_t> tfxIcospherePoints[6];

void tfx__make_icospheres() {
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

	tfx_vector_t<tfx_vec3_t> vertices;
	vertices.push_back({ -x, n, z });
	vertices.push_back({ x, n, z });
	vertices.push_back({ -x, n,-z });
	vertices.push_back({ x, n,-z });
	vertices.push_back({ n, z, x });
	vertices.push_back({ n, z,-x });
	vertices.push_back({ n,-z, x });
	vertices.push_back({ n,-z,-x });
	vertices.push_back({ z, x, n });
	vertices.push_back({ -z, x, n });
	vertices.push_back({ z,-x, n });
	vertices.push_back({ -z,-x, n });

	tfx_vector_t<tfx_face_t> triangles;
	triangles.push_back({ 0,  4, 1 });
	triangles.push_back({ 0 , 9, 4 });
	triangles.push_back({ 9 , 5, 4 });
	triangles.push_back({ 4 , 5, 8 });
	triangles.push_back({ 4 , 8, 1 });
	triangles.push_back({ 8 ,10, 1 });
	triangles.push_back({ 8 , 3,10 });
	triangles.push_back({ 5 , 3, 8 });
	triangles.push_back({ 5 , 2, 3 });
	triangles.push_back({ 2 , 7, 3 });
	triangles.push_back({ 7 ,10, 3 });
	triangles.push_back({ 7 , 6,10 });
	triangles.push_back({ 7 ,11, 6 });
	triangles.push_back({ 11, 0, 6 });
	triangles.push_back({ 0 , 1, 6 });
	triangles.push_back({ 6 , 1,10 });
	triangles.push_back({ 9 , 0,11 });
	triangles.push_back({ 9 ,11, 2 });
	triangles.push_back({ 9 , 2, 5 });
	triangles.push_back({ 7 , 2,11 });

	tfx_storage_map_t<int> point_cache;

	for (int i = 0; i < subdivisions; ++i)
	{
		triangles = tfx__sub_divide_icosphere(&point_cache, &vertices, &triangles);
		TFX_ASSERT(tfxIcospherePoints[i].capacity == vertices.current_size);    //Must be the same size
		memcpy(tfxIcospherePoints[i].data, vertices.data, vertices.current_size * sizeof(tfx_vec3_t));
		tfxIcospherePoints[i].current_size = vertices.current_size;
		std::qsort(tfxIcospherePoints[i].data, tfxIcospherePoints[i].current_size, sizeof(tfx_vec3_t), tfx__sort_icosphere_points);
	}

}

int tfx__vertex_for_edge(tfx_storage_map_t<int> *point_cache, tfx_vector_t<tfx_vec3_t> *vertices, int first, int second)
{
	tfxKey key = ((tfxKey)first << 32) + second;
	if (first > second)
		key = ((tfxKey)second << 32) + first;

	if (point_cache->ValidKey(key))
		return point_cache->At(key);

	point_cache->Insert(key, vertices->size());

	tfx_vec3_t edge_sum = (*vertices)[first] + (*vertices)[second];
	tfx_vec3_t point = tfx__normalize_vec3(&edge_sum);
	vertices->push_back(point);

	return vertices->size() - 1;
}

tfx_vector_t<tfx_face_t> tfx__sub_divide_icosphere(tfx_storage_map_t<int> *point_cache, tfx_vector_t<tfx_vec3_t> *vertices, tfx_vector_t<tfx_face_t> *triangles)
{
	tfx_vector_t<tfx_face_t> result;

	for (tfx_face_t face : *triangles)
	{
		int mid[3];
		for (int edge = 0; edge < 3; ++edge)
		{
			mid[edge] = tfx__vertex_for_edge(point_cache, vertices, face.v[edge], face.v[(edge + 1) % 3]);
		}

		result.push_back({ face.v[0], mid[0], mid[2] });
		result.push_back({ face.v[1], mid[1], mid[0] });
		result.push_back({ face.v[2], mid[2], mid[1] });
		result.push_back({ mid[0], mid[1], mid[2] });
	}

	return result;
}

int tfx__sort_icosphere_points(void const *left, void const *right) {
	float d1 = static_cast<const tfx_vec3_t *>(left)->y;
	float d2 = static_cast<const tfx_vec3_t *>(right)->y;
	return (d2 > d1) - (d2 < d1);
}

int tfx__sort_depth(void const *left, void const *right) {
	float d1 = static_cast<const tfx_depth_index_t *>(left)->depth;
	float d2 = static_cast<const tfx_depth_index_t *>(right)->depth;
	return (d2 > d1) - (d2 < d1);
}

void tfx__insertion_sort_depth(tfx_work_queue_t *queue, void *work_entry) {
	tfx_bucket_array_t<tfx_particle_soa_t> &bank = *static_cast<tfx_sort_work_entry_t *>(work_entry)->bank;
	tfx_vector_t<tfx_depth_index_t> &depth_indexes = *static_cast<tfx_sort_work_entry_t *>(work_entry)->depth_indexes;
	for (tfxU32 i = 1; i < depth_indexes.current_size; ++i) {
		tfx_depth_index_t key = depth_indexes[i];
		int j = i - 1;
		while (j >= 0 && key.depth > depth_indexes[j].depth) {
			depth_indexes[j + 1] = depth_indexes[j];
			bank[tfx__particle_bank(depth_indexes[j + 1].particle_id)].depth_index[tfx__particle_index(depth_indexes[j + 1].particle_id)] = j + 1;
			--j;
		}
		depth_indexes[j + 1] = key;
		bank[tfx__particle_bank(depth_indexes[j + 1].particle_id)].depth_index[tfx__particle_index(depth_indexes[j + 1].particle_id)] = j + 1;
	}
}

tfx_hsv_t tfx__rgb_to_hsv(tfx_rgb_t in)
{
	tfx_hsv_t      out;
	float      min, max, delta;

	min = in.r < in.g ? in.r : in.g;
	min = min < in.b ? min : in.b;

	max = in.r > in.g ? in.r : in.g;
	max = max > in.b ? max : in.b;

	out.v = max;                                // v
	delta = max - min;
	if (delta < 0.00001f)
	{
		out.s = 0;
		out.h = 0; // undefined, maybe nan?
		return out;
	}
	if (max > 0.0f) { // NOTE: if Max is == 0, this divide would cause a crash
		out.s = (delta / max);                  // s
	}
	else {
		// if max is 0, then r = g = b = 0              
		// s = 0, h is undefined
		out.s = 0.0f;
		out.h = NAN;                            // its now undefined
		return out;
	}
	if (in.r >= max)                           // > is bogus, just keeps compilor happy
		out.h = (in.g - in.b) / delta;        // between yellow & magenta
	else
		if (in.g >= max)
			out.h = 2.0f + (in.b - in.r) / delta;  // between cyan & yellow
		else
			out.h = 4.0f + (in.r - in.g) / delta;  // between magenta & cyan

	out.h *= 60.0f;                              // degrees

	if (out.h < 0.0f)
		out.h += 360.0f;

	return out;
}


tfx_rgb_t tfx__hsv_to_rgb(tfx_hsv_t in)
{
	float      hh, p, q, t, ff;
	long        i;
	tfx_rgb_t      out;

	if (in.s <= 0.0f) {       // < is bogus, just shuts up warnings
		out.r = in.v;
		out.g = in.v;
		out.b = in.v;
		return out;
	}
	hh = in.h;
	if (hh >= 360.0f) hh = 0.0f;
	hh /= 60.0f;
	i = (long)hh;
	ff = hh - i;
	p = in.v * (1.0f - in.s);
	q = in.v * (1.0f - (in.s * ff));
	t = in.v * (1.0f - (in.s * (1.0f - ff)));

	switch (i) {
	case 0:
		out.r = in.v;
		out.g = t;
		out.b = p;
		break;
	case 1:
		out.r = q;
		out.g = in.v;
		out.b = p;
		break;
	case 2:
		out.r = p;
		out.g = in.v;
		out.b = t;
		break;

	case 3:
		out.r = p;
		out.g = q;
		out.b = in.v;
		break;
	case 4:
		out.r = t;
		out.g = p;
		out.b = in.v;
		break;
	case 5:
	default:
		out.r = in.v;
		out.g = p;
		out.b = q;
		break;
	}
	return out;
}

float tfx__degrees_to_radians(float degrees) { return degrees * 0.01745329251994329576923690768489f; }
float tfx__radians_to_degrees(float radians) { return radians * 57.295779513082320876798154814105f; }

float tfx__length_vec3_nosqr(tfx_vec3_t const *v) {
	return v->x * v->x + v->y * v->y + v->z * v->z;
}

float tfx__length_vec4_nosqr(tfx_vec4_t const *v) {
	return v->x * v->x + v->y * v->y + v->z * v->z + v->w * v->w;
}

float tfx__length_vec3(tfx_vec3_t const *v) {
	return sqrtf(tfx__length_vec3_nosqr(v));
}

float tfx__length_vec4(tfx_vec4_t const *v) {
	return sqrtf(tfx__length_vec4_nosqr(v));
}

float tfx__has_length_vec3(tfx_vec3_t const *v) {
	return (v->x == 0 && v->y == 0 && v->z == 0) ? 0.f : 1.f;
}

tfx_vec3_t tfx__normalize_vec3(tfx_vec3_t const *v) {
	float length = tfx__length_vec3(v);
	return length > 0.f ? tfx_vec3_t(v->x / length, v->y / length, v->z / length) : *v;
}

tfx_vec4_t tfx__normalize_vec4(tfx_vec4_t const *v) {
	if (v->x == 0 && v->y == 0 && v->z == 0 && v->w == 0) return tfx_vec4_t(1.f, 0.f, 0.f, 0.f);
	float length = tfx__length_vec4(v);
	return tfx_vec4_t(v->x / length, v->y / length, v->z / length, v->w / length);
}

tfx_vec3_t tfx__cross_product_vec3(tfx_vec3_t *a, tfx_vec3_t *b) {
	tfx_vec3_t result;
	result.x = a->y * b->z - a->z * b->y;
	result.y = a->z * b->x - a->x * b->z;
	result.z = a->x * b->y - a->y * b->x;
	return(result);
}

void tfx__wide_cross_product(tfxWideFloat ax, tfxWideFloat ay, tfxWideFloat az, tfxWideFloat *bx, tfxWideFloat *by, tfxWideFloat *bz, tfxWideFloat *rx, tfxWideFloat *ry, tfxWideFloat *rz) {
	*rx = tfxWideSub(tfxWideMul(ay, *bz), tfxWideMul(az, *by));
	*ry = tfxWideSub(tfxWideMul(az, *bx), tfxWideMul(ax, *bz));
	*rz = tfxWideSub(tfxWideMul(ax, *by), tfxWideMul(ay, *bx));
}

float tfx__dot_product_vec4(const tfx_vec4_t *a, const tfx_vec4_t *b)
{
	return (a->x * b->x + a->y * b->y + a->z * b->z + a->w * b->w);
}

float tfx__dot_product_vec3(const tfx_vec3_t *a, const tfx_vec3_t *b)
{
	return (a->x * b->x + a->y * b->y + a->z * b->z);
}

float tfx__dot_product_vec2(const tfx_vec2_t *a, const tfx_vec2_t *b)
{
	return (a->x * b->x + a->y * b->y);
}

tfx_quaternion_t QuaternionFromAxisAngle(float x, float y, float z, float angle) {
	float half_angle = angle * .5f;
	float sin_half_angle = sinf(half_angle);

	return tfx_quaternion_t(cos(half_angle), x * sin_half_angle, y * sin_half_angle, z * sin_half_angle);
}

tfx_quaternion_t QuaternionFromDirection(tfx_vec3_t *normalised_dir) {
	// Initial direction (default y-axis) because this is how paths are generated
	tfx_vec3_t initial_dir = { 0.0f, 1.0f, 0.0f };

	// Calculate rotation axis
	tfx_vec3_t rotation_axis = tfx__cross_product_vec3(&initial_dir, normalised_dir);
	rotation_axis = tfx__normalize_vec3_fast(&rotation_axis);

	// Calculate dot product
	float dot_product = tfx__dot_product_vec3(&initial_dir, normalised_dir);

	// Calculate rotation angle
	float angle = acosf(dot_product);  // Angle between initial_dir and normalised_dir
	float half_angle = angle / 2.0f;
	float sin_half_angle = sinf(half_angle);

	// Construct quaternion
	tfx_quaternion_t result(
		cosf(half_angle),                    // w
		rotation_axis.x * sin_half_angle,    // x
		rotation_axis.y * sin_half_angle,    // y
		rotation_axis.z * sin_half_angle     // z
	);

	return NormalizeQuaternion(&result);
}

tfx_vec3_t tfx__cylinder_surface_normal(float x, float z, float width, float depth) {
	// Calculate the gradient of the ellipse equation
	float dx = 2.f * x / (width * width);
	float dz = 2.f * z / (depth * depth);

	// Normalize the gradient vector to obtain the surface normal
	float length = 1.f / tfx__quake_sqrt(dx * dx + dz * dz);
	tfx_vec3_t normal;
	normal.x = dx / length;
	normal.y = 0.f;
	normal.z = dz / length;

	return normal;
}

tfx_vec3_t tfx__ellipse_surface_normal(float x, float y, float z, float width, float height, float depth) {
	float dx = 2.f * x / (width * width);
	float dy = 2.f * y / (height * height);
	float dz = 2.f * z / (depth * depth);

	// Normalize the gradient vector to obtain the surface normal
	float length = 1.f / tfx__quake_sqrt(dx * dx + dy * dy + dz * dz);
	tfx_vec3_t normal;
	normal.x = dx / length;
	normal.y = dy / length;
	normal.z = dz / length;

	return normal;
}

void tfx__wide_ellipse_surface_normal(const tfxWideFloat *x, const tfxWideFloat *y, const tfxWideFloat *z, const tfxWideFloat *width, const tfxWideFloat *height, const tfxWideFloat *depth, tfxWideFloat *normal_x, tfxWideFloat *normal_y, tfxWideFloat *normal_z) {
	// Calculate the gradient of the ellipse equation
	tfxWideFloat w2 = tfxWideMul(*width, *width);
	tfxWideFloat h2 = tfxWideMul(*height, *height);
	tfxWideFloat d2 = tfxWideMul(*depth, *depth);

	tfxWideFloat dx = tfxWideDiv(*x, w2);
	tfxWideFloat dy = tfxWideDiv(*y, h2);
	tfxWideFloat dz = tfxWideDiv(*z, d2);

	// Normalize the gradient vector to obtain the surface normal
	tfxWideFloat length = tfxWideMul(dx, dx);
	length = tfxWideAdd(length, tfxWideMul(dy, dy));
	length = tfxWideAdd(length, tfxWideMul(dz, dz));
#ifdef tfxARM
	length = tfxWideMul(tfxWideRSqrt(length), length);
#else
	length = tfxWideSqrt(length);
#endif
	tfxWideFloat scale = tfxWideDiv(tfxWIDEONE, length);

	*normal_x = tfxWideDiv(tfxWideMul(scale, *x), w2);
	*normal_y = tfxWideDiv(tfxWideMul(scale, *y), h2);
	*normal_z = tfxWideDiv(tfxWideMul(scale, *z), d2);
}

tfx_vec2_t tfx__catmull_rom_spline_2d(const tfx_vec4_t *p0, const tfx_vec4_t *p1, const tfx_vec4_t *p2, const tfx_vec4_t *p3, float t) {
	float t2 = t * t;
	float t3 = t2 * t;

	float b0 = -t3 + 2.0f * t2 - t;
	float b1 = 3.0f * t3 - 5.0f * t2 + 2.0f;
	float b2 = -3.0f * t3 + 4.0f * t2 + t;
	float b3 = t3 - t2;

	float x = p0->x * b0 + p1->x * b1 + p2->x * b2 + p3->x * b3;
	float y = p0->y * b0 + p1->y * b1 + p2->y * b2 + p3->y * b3;

	return { x * 0.5f, y * 0.5f };
}

tfx_vec2_t tfx__catmull_rom_spline_2d_soa(const float *p_x, const float *p_y, int p0, float t) {
	float t2 = t * t;
	float t3 = t2 * t;

	int p1 = p0 + 1;
	int p2 = p0 + 2;
	int p3 = p0 + 3;

	float b0 = -t3 + 2.0f * t2 - t;
	float b1 = 3.0f * t3 - 5.0f * t2 + 2.0f;
	float b2 = -3.0f * t3 + 4.0f * t2 + t;
	float b3 = t3 - t2;

	float x = p_x[p0] * b0 + p_x[p1] * b1 + p_x[p2] * b2 + p_x[p3] * b3;
	float y = p_y[p0] * b0 + p_y[p1] * b1 + p_y[p2] * b2 + p_y[p3] * b3;

	return { x * .5f, y * .5f };
}

tfx_vec2_t tfx__catmull_rom_spline_gradient_2d_soa(const float *px, const float *py, float t) {
	float t2 = t * t;

	float b0 = -3.f * t2 + 4.f * t - 1.f;
	float b1 = 9.f * t2 - 10.f * t;
	float b2 = -9.f * t2 + 8.f * t + 1.f;
	float b3 = 3.f * t2 - 2.f * t;

	float x = px[0] * b0 + px[1] * b1 + px[2] * b2 + px[3] * b3;
	float y = py[0] * b0 + py[1] * b1 + py[2] * b2 + py[3] * b3;

	return { x * 0.5f, y * 0.5f };
}

tfx_vec3_t tfx__catmull_rom_spline_3d_soa(const float *p_x, const float *p_y, const float *p_z, int p0, float t) {
	float t2 = t * t;
	float t3 = t2 * t;

	int p1 = p0 + 1;
	int p2 = p0 + 2;
	int p3 = p0 + 3;

	float b0 = -t3 + 2.0f * t2 - t;
	float b1 = 3.0f * t3 - 5.0f * t2 + 2.0f;
	float b2 = -3.0f * t3 + 4.0f * t2 + t;
	float b3 = t3 - t2;

	float x = p_x[p0] * b0 + p_x[p1] * b1 + p_x[p2] * b2 + p_x[p3] * b3;
	float y = p_y[p0] * b0 + p_y[p1] * b1 + p_y[p2] * b2 + p_y[p3] * b3;
	float z = p_z[p0] * b0 + p_z[p1] * b1 + p_z[p2] * b2 + p_z[p3] * b3;

	return { x * .5f, y * .5f, z * .5f };
}

tfx_vec3_t tfx__catmull_rom_spline_3d(const tfx_vec4_t *p0, const tfx_vec4_t *p1, const tfx_vec4_t *p2, const tfx_vec4_t *p3, float t) {
	float t2 = t * t;
	float t3 = t2 * t;

	float b0 = -t3 + 2.0f * t2 - t;
	float b1 = 3.0f * t3 - 5.0f * t2 + 2.0f;
	float b2 = -3.0f * t3 + 4.0f * t2 + t;
	float b3 = t3 - t2;

	float x = p0->x * b0 + p1->x * b1 + p2->x * b2 + p3->x * b3;
	float y = p0->y * b0 + p1->y * b1 + p2->y * b2 + p3->y * b3;
	float z = p0->z * b0 + p1->z * b1 + p2->z * b2 + p3->z * b3;

	return { x * 0.5f, y * 0.5f, z * 0.5f };
}

tfx_vec3_t tfx__catmull_rom_spline_gradient_3d(const tfx_vec4_t *p0, const tfx_vec4_t *p1, const tfx_vec4_t *p2, const tfx_vec4_t *p3, float t) {
	float t2 = t * t;

	float b0 = -3.f * t2 + 4.f * t - 1.f;
	float b1 = 9.f * t2 - 10.f * t;
	float b2 = -9.f * t2 + 8.f * t + 1.f;
	float b3 = 3.f * t2 - 2.f * t;

	float x = p0->x * b0 + p1->x * b1 + p2->x * b2 + p3->x * b3;
	float y = p0->y * b0 + p1->y * b1 + p2->y * b2 + p3->y * b3;
	float z = p0->z * b0 + p1->z * b1 + p2->z * b2 + p3->z * b3;

	return { x * 0.5f, y * 0.5f, z * 0.5f };
}

tfx_vec3_t tfx__catmull_rom_spline_gradient_3d_soa(const float *px, const float *py, const float *pz, float t) {
	float t2 = t * t;

	float b0 = -3.f * t2 + 4.f * t - 1.f;
	float b1 = 9.f * t2 - 10.f * t;
	float b2 = -9.f * t2 + 8.f * t + 1.f;
	float b3 = 3.f * t2 - 2.f * t;

	float x = px[0] * b0 + px[1] * b1 + px[2] * b2 + px[3] * b3;
	float y = py[0] * b0 + py[1] * b1 + py[2] * b2 + py[3] * b3;
	float z = pz[0] * b0 + pz[1] * b1 + pz[2] * b2 + pz[3] * b3;

	return { x * 0.5f, y * 0.5f, z * 0.5f };
}

void tfx__wide_catmull_rom_spline_2d(tfxWideArrayi *pi, tfxWideFloat t, float *x, float *y, tfxWideFloat *vx, tfxWideFloat *vy) {
	//This calculates the position on a catmull rom spline for 4 (sse) or 8 (avx) particles at a time.
	//pi contains the first index in the path node list, t is the % of the segment on the path to calcuate for. 
	tfxWideFloat t2 = tfxWideMul(t, t);
	tfxWideFloat t3 = tfxWideMul(t2, t);

	//Calculate the weigths
	tfxWideFloat b0 = tfxWideSub(tfxWideAdd(tfxWideMul(t3, tfxWideSetSingle(-1.f)), tfxWideMul(tfxWideSetSingle(2.f), t2)), t);
	tfxWideFloat b1 = tfxWideAdd(tfxWideSub(tfxWideMul(tfxWideSetSingle(3.f), t3), tfxWideMul(tfxWideSetSingle(5.f), t2)), tfxWideSetSingle(2.f));
	tfxWideFloat b2 = tfxWideAdd(tfxWideAdd(tfxWideMul(tfxWideSetSingle(-3.f), t3), tfxWideMul(tfxWideSetSingle(4.f), t2)), t);
	tfxWideFloat b3 = tfxWideSub(t3, t2);

	//Load in the node data for 4 nodes to to calculate the path
	tfxWideFloat px0 = tfxWideLookupSet(x, (*pi));
	tfxWideFloat py0 = tfxWideLookupSet(y, (*pi));
	tfxWideFloat px1 = tfxWideLookupSetOffset(x, (*pi), 1);
	tfxWideFloat py1 = tfxWideLookupSetOffset(y, (*pi), 1);
	tfxWideFloat px2 = tfxWideLookupSetOffset(x, (*pi), 2);
	tfxWideFloat py2 = tfxWideLookupSetOffset(y, (*pi), 2);
	tfxWideFloat px3 = tfxWideLookupSetOffset(x, (*pi), 3);
	tfxWideFloat py3 = tfxWideLookupSetOffset(y, (*pi), 3);

	*vx = tfxWideMul(tfxWideAdd(tfxWideAdd(tfxWideAdd(tfxWideMul(px0, b0), tfxWideMul(px1, b1)), tfxWideMul(px2, b2)), tfxWideMul(px3, b3)), tfxWideSetSingle(.5f));
	*vy = tfxWideMul(tfxWideAdd(tfxWideAdd(tfxWideAdd(tfxWideMul(py0, b0), tfxWideMul(py1, b1)), tfxWideMul(py2, b2)), tfxWideMul(py3, b3)), tfxWideSetSingle(.5f));
}

void tfx__wide_catmull_tom_spline_3d(tfxWideArrayi *pi, tfxWideFloat t, float *x, float *y, float *z, tfxWideFloat *vx, tfxWideFloat *vy, tfxWideFloat *vz) {
	//This calculates the position on a catmull rom spline for 4 (sse) or 8 (avx) particles at a time.
	//pi contains the first index in the path node list, t is the % of the segment on the path to calcuate for. 
	tfxWideFloat t2 = tfxWideMul(t, t);
	tfxWideFloat t3 = tfxWideMul(t2, t);

	//Calculate the weigths
	tfxWideFloat b0 = tfxWideSub(tfxWideAdd(tfxWideMul(t3, tfxWideSetSingle(-1.f)), tfxWideMul(tfxWideSetSingle(2.f), t2)), t);
	tfxWideFloat b1 = tfxWideAdd(tfxWideSub(tfxWideMul(tfxWideSetSingle(3.f), t3), tfxWideMul(tfxWideSetSingle(5.f), t2)), tfxWideSetSingle(2.f));
	tfxWideFloat b2 = tfxWideAdd(tfxWideAdd(tfxWideMul(tfxWideSetSingle(-3.f), t3), tfxWideMul(tfxWideSetSingle(4.f), t2)), t);
	tfxWideFloat b3 = tfxWideSub(t3, t2);

	//Load in the node data for 4 nodes to to calculate the path
	tfxWideFloat px0 = tfxWideLookupSet(x, (*pi));
	tfxWideFloat py0 = tfxWideLookupSet(y, (*pi));
	tfxWideFloat pz0 = tfxWideLookupSet(z, (*pi));
	tfxWideFloat px1 = tfxWideLookupSetOffset(x, (*pi), 1);
	tfxWideFloat py1 = tfxWideLookupSetOffset(y, (*pi), 1);
	tfxWideFloat pz1 = tfxWideLookupSetOffset(z, (*pi), 1);
	tfxWideFloat px2 = tfxWideLookupSetOffset(x, (*pi), 2);
	tfxWideFloat py2 = tfxWideLookupSetOffset(y, (*pi), 2);
	tfxWideFloat pz2 = tfxWideLookupSetOffset(z, (*pi), 2);
	tfxWideFloat px3 = tfxWideLookupSetOffset(x, (*pi), 3);
	tfxWideFloat py3 = tfxWideLookupSetOffset(y, (*pi), 3);
	tfxWideFloat pz3 = tfxWideLookupSetOffset(z, (*pi), 3);

	*vx = tfxWideMul(tfxWideAdd(tfxWideAdd(tfxWideAdd(tfxWideMul(px0, b0), tfxWideMul(px1, b1)), tfxWideMul(px2, b2)), tfxWideMul(px3, b3)), tfxWideSetSingle(.5f));
	*vy = tfxWideMul(tfxWideAdd(tfxWideAdd(tfxWideAdd(tfxWideMul(py0, b0), tfxWideMul(py1, b1)), tfxWideMul(py2, b2)), tfxWideMul(py3, b3)), tfxWideSetSingle(.5f));
	*vz = tfxWideMul(tfxWideAdd(tfxWideAdd(tfxWideAdd(tfxWideMul(pz0, b0), tfxWideMul(pz1, b1)), tfxWideMul(pz2, b2)), tfxWideMul(pz3, b3)), tfxWideSetSingle(.5f));
}

void tfx__wide_catmull_rom_spline_gradient_2d(tfxWideArrayi *pi, tfxWideFloat t, float *x, float *y, tfxWideFloat *vx, tfxWideFloat *vy) {
	tfxWideFloat t2 = tfxWideMul(t, t);

	//Calculate the weigths
	tfxWideFloat b0 = tfxWideSub(tfxWideAdd(tfxWideMul(tfxWideSetSingle(-3.f), t2), tfxWideMul(tfxWideSetSingle(4.f), t)), tfxWideSetSingle(1.f));
	tfxWideFloat b1 = tfxWideSub(tfxWideMul(tfxWideSetSingle(9.f), t2), tfxWideMul(tfxWideSetSingle(10.f), t));
	tfxWideFloat b2 = tfxWideAdd(tfxWideAdd(tfxWideMul(tfxWideSetSingle(-9.f), t2), tfxWideMul(tfxWideSetSingle(8.f), t)), tfxWideSetSingle(1.f));
	tfxWideFloat b3 = tfxWideSub(tfxWideMul(tfxWideSetSingle(3.f), t2), tfxWideMul(tfxWideSetSingle(2.f), t));

	//Load in the node data for 4 nodes to to calculate the path
	tfxWideFloat px0 = tfxWideLookupSet(x, (*pi));
	tfxWideFloat py0 = tfxWideLookupSet(y, (*pi));
	tfxWideFloat px1 = tfxWideLookupSetOffset(x, (*pi), 1);
	tfxWideFloat py1 = tfxWideLookupSetOffset(y, (*pi), 1);
	tfxWideFloat px2 = tfxWideLookupSetOffset(x, (*pi), 2);
	tfxWideFloat py2 = tfxWideLookupSetOffset(y, (*pi), 2);
	tfxWideFloat px3 = tfxWideLookupSetOffset(x, (*pi), 3);
	tfxWideFloat py3 = tfxWideLookupSetOffset(y, (*pi), 3);

	tfxWideFloat half = tfxWideSetSingle(0.5f);
	*vx = tfxWideMul(tfxWideAdd(tfxWideAdd(tfxWideAdd(tfxWideMul(px0, b0), tfxWideMul(px1, b1)), tfxWideMul(px2, b2)), tfxWideMul(px3, b3)), half);
	*vy = tfxWideMul(tfxWideAdd(tfxWideAdd(tfxWideAdd(tfxWideMul(py0, b0), tfxWideMul(py1, b1)), tfxWideMul(py2, b2)), tfxWideMul(py3, b3)), half);
}

void tfx__wide_catmull_rom_spline_gradient_3d(tfxWideArrayi *pi, tfxWideFloat t, float *x, float *y, float *z, tfxWideFloat *vx, tfxWideFloat *vy, tfxWideFloat *vz) {
	tfxWideFloat t2 = tfxWideMul(t, t);

	//Calculate the weigths
	tfxWideFloat b0 = tfxWideSub(tfxWideAdd(tfxWideMul(tfxWideSetSingle(-3.f), t2), tfxWideMul(tfxWideSetSingle(4.f), t)), tfxWideSetSingle(1.f));
	tfxWideFloat b1 = tfxWideSub(tfxWideMul(tfxWideSetSingle(9.f), t2), tfxWideMul(tfxWideSetSingle(10.f), t));
	tfxWideFloat b2 = tfxWideAdd(tfxWideAdd(tfxWideMul(tfxWideSetSingle(-9.f), t2), tfxWideMul(tfxWideSetSingle(8.f), t)), tfxWideSetSingle(1.f));
	tfxWideFloat b3 = tfxWideSub(tfxWideMul(tfxWideSetSingle(3.f), t2), tfxWideMul(tfxWideSetSingle(2.f), t));

	//Load in the node data for 4 nodes to to calculate the path
	tfxWideFloat px0 = tfxWideLookupSet(x, (*pi));
	tfxWideFloat py0 = tfxWideLookupSet(y, (*pi));
	tfxWideFloat pz0 = tfxWideLookupSet(z, (*pi));
	tfxWideFloat px1 = tfxWideLookupSetOffset(x, (*pi), 1);
	tfxWideFloat py1 = tfxWideLookupSetOffset(y, (*pi), 1);
	tfxWideFloat pz1 = tfxWideLookupSetOffset(z, (*pi), 1);
	tfxWideFloat px2 = tfxWideLookupSetOffset(x, (*pi), 2);
	tfxWideFloat py2 = tfxWideLookupSetOffset(y, (*pi), 2);
	tfxWideFloat pz2 = tfxWideLookupSetOffset(z, (*pi), 2);
	tfxWideFloat px3 = tfxWideLookupSetOffset(x, (*pi), 3);
	tfxWideFloat py3 = tfxWideLookupSetOffset(y, (*pi), 3);
	tfxWideFloat pz3 = tfxWideLookupSetOffset(z, (*pi), 3);

	*vx = tfxWideMul(tfxWideAdd(tfxWideAdd(tfxWideAdd(tfxWideMul(px0, b0), tfxWideMul(px1, b1)), tfxWideMul(px2, b2)), tfxWideMul(px3, b3)), tfxWideSetSingle(.5f));
	*vy = tfxWideMul(tfxWideAdd(tfxWideAdd(tfxWideAdd(tfxWideMul(py0, b0), tfxWideMul(py1, b1)), tfxWideMul(py2, b2)), tfxWideMul(py3, b3)), tfxWideSetSingle(.5f));
	*vz = tfxWideMul(tfxWideAdd(tfxWideAdd(tfxWideAdd(tfxWideMul(pz0, b0), tfxWideMul(pz1, b1)), tfxWideMul(pz2, b2)), tfxWideMul(pz3, b3)), tfxWideSetSingle(.5f));

}

//Quake 3 inverse square root
float tfx__quake_sqrt(float number)
{
	union {
		float f;
		uint32_t i;
	} conv;

	float x2;
	const float threehalfs = 1.5F;

	x2 = number * 0.5F;
	conv.f = number;
	conv.i = 0x5f3759df - (conv.i >> 1);
	conv.f = conv.f * (threehalfs - (x2 * conv.f * conv.f));
	return conv.f;
}

float tfx__vec2_length_fast(tfx_vec2_t const *v) {
	return 1.f / tfx__quake_sqrt(tfx__dot_product_vec2(v, v));
}

float tfx__vec3_length_fast(tfx_vec3_t const *v) {
	return 1.f / tfx__quake_sqrt(tfx__dot_product_vec3(v, v));
}

tfx_vec3_t tfx__normalize_vec3_fast(tfx_vec3_t const *v) {
	return *v * tfx__quake_sqrt(tfx__dot_product_vec3(v, v));
}

tfx_vec2_t tfx__normalize_vec2(tfx_vec2_t const *v) {
	float length = tfx__vec2_length_fast(v);
	return tfx_vec2_t(v->x / length, v->y / length);
}

tfx_mat3_t tfx__create_matrix3(float v) {
	tfx_mat3_t R =
	{ {
		{v, 0, 0},
		{0, v, 0},
		{0, 0, v}},
	};
	return(R);
}

tfx_mat3_t tfx__translate_matrix3_vec3(tfx_mat3_t const *m, tfx_vec3_t const *v) {
	tfx_mat3_t result;
	result.v[2] = m->v[0] * v->x + m->v[1] * v->y + m->v[2];
	return result;
}

tfx_mat3_t tfx__rotate_matrix3(tfx_mat3_t const *m, float r) {
	float const a = r;
	float const c = cosf(a);
	float const s = sinf(a);

	tfx_mat3_t result;
	result.v[0] = m->v[0] * c + m->v[1] * s;
	result.v[1] = m->v[0] * -s + m->v[1] * c;
	result.v[2] = m->v[2];
	return result;
}

tfx_mat3_t tfx__scale_matrix3_vec2(tfx_mat3_t const *m, tfx_vec2_t const &v) {
	tfx_mat3_t result;
	result.v[0] = m->v[0] * v.x;
	result.v[1] = m->v[1] * v.y;
	result.v[2] = m->v[2];
	return result;
}

tfx_vec2_t tfx__transform_vec2_matrix4(const tfx_mat4_t *mat, const tfx_vec2_t v) {
	tfx_vec2_t tv = tfx_vec2_t(0.f, 0.f);
	tv.x = v.x * mat->v[0].x + v.y * mat->v[1].x;
	tv.y = v.x * mat->v[0].y + v.y * mat->v[1].y;
	return tv;
}

tfx_mat4_t tfx__create_matrix4(float v) {
	tfx_mat4_t R =
	{ {
		{v, 0, 0, 0},
		{0, v, 0, 0},
		{0, 0, v, 0},
		{0, 0, 0, v}},
	};
	return(R);
}

tfx_mat4_t tfx__matrix4_from_vecs(tfx_vec4_t a, tfx_vec4_t b, tfx_vec4_t c, tfx_vec4_t d) {
	tfx_mat4_t R =
	{ {
		{a.x, a.y, a.z, a.w},
		{b.x, b.y, b.z, b.w},
		{c.x, c.y, c.z, c.w},
		{d.x, d.y, d.z, d.w}},
	};
	return(R);
}

tfx_mat4_t tfx__matrix4_rotate_x(float angle) {
	float c = cosf(angle);
	float s = sinf(angle);
	tfx_mat4_t r =
	{ {
		{1, 0, 0, 0},
		{0, c,-s, 0},
		{0, s, c, 0},
		{0, 0, 0, 1}},
	};
	return r;
}

tfx_mat4_t tfx__matrix4_rotate_y(float angle) {
	float c = cosf(angle);
	float s = sinf(angle);
	tfx_mat4_t r =
	{ {
		{ c, 0, s, 0},
		{ 0, 1, 0, 0},
		{-s, 0, c, 0},
		{ 0, 0, 0, 1}},
	};
	return r;
}

tfx_mat4_t tfx__matrix4_rotate_z(float angle) {
	float c = cosf(angle);
	float s = sinf(angle);
	tfx_mat4_t r =
	{ {
		{c, -s, 0, 0},
		{s,  c, 0, 0},
		{0,  0, 1, 0},
		{0,  0, 0, 1}},
	};
	return r;
}

tfx_mat4_t tfx__transpose_matrix4(tfx_mat4_t *mat) {
	return tfx__matrix4_from_vecs(
		tfx_vec4_t(mat->v[0].x, mat->v[1].x, mat->v[2].x, mat->v[3].x),
		tfx_vec4_t(mat->v[0].y, mat->v[1].y, mat->v[2].y, mat->v[3].y),
		tfx_vec4_t(mat->v[0].z, mat->v[1].z, mat->v[2].z, mat->v[3].z),
		tfx_vec4_t(mat->v[0].w, mat->v[1].w, mat->v[2].w, mat->v[3].w)
	);
}

tfx_mat4_t tfx__transform_matrix4_matrix4_as_2(const tfx_mat4_t *in, const tfx_mat4_t *m) {
	tfx_mat4_t r;
	r.v[0].x = in->v[0].x * m->v[0].x + in->v[0].y * m->v[1].x; r.v[0].y = in->v[0].x * m->v[0].y + in->v[0].y * m->v[1].y;
	r.v[1].x = in->v[1].x * m->v[0].x + in->v[1].y * m->v[1].x; r.v[1].y = in->v[1].x * m->v[0].y + in->v[1].y * m->v[1].y;
	return r;
}

tfx_mat4_t tfx__transform_matrix4_matrix2(const tfx_mat4_t *in, const tfx_mat2_t *m) {
	tfx_mat4_t r;
	r.v[0].x = in->v[0].x * m->aa + in->v[0].y * m->ba; r.v[0].y = in->v[0].x * m->ab + in->v[0].y * m->bb;
	r.v[1].x = in->v[1].x * m->aa + in->v[1].y * m->ba; r.v[1].y = in->v[1].x * m->ab + in->v[1].y * m->bb;
	return r;
}

void tfx__wide_transform_quaternion_vec3(const tfx_quaternion_t *q, tfxWideFloat *x, tfxWideFloat *y, tfxWideFloat *z) {
	tfxWideFloat qv_x = *x;
	tfxWideFloat qv_y = *y;
	tfxWideFloat qv_z = *z;
	tfxWideFloat qv_w = tfxWideSetZero;

	tfxWideFloat q_x = tfxWideSetSingle(q->x);
	tfxWideFloat q_y = tfxWideSetSingle(q->y);
	tfxWideFloat q_z = tfxWideSetSingle(q->z);
	tfxWideFloat q_w = tfxWideSetSingle(q->w);

	tfxWideFloat two = tfxWideSetSingle(2.0f);

	tfxWideFloat t_x = tfxWideSub(tfxWideMul(q_y, qv_z), tfxWideMul(q_z, qv_y));
	t_x = tfxWideMul(t_x, two);
	tfxWideFloat t_y = tfxWideSub(tfxWideMul(q_z, qv_x), tfxWideMul(q_x, qv_z));
	t_y = tfxWideMul(t_y, two);
	tfxWideFloat t_z = tfxWideSub(tfxWideMul(q_x, qv_y), tfxWideMul(q_y, qv_x));
	t_z = tfxWideMul(t_z, two);

	tfxWideFloat qw_t_x = tfxWideMul(q_w, t_x);
	tfxWideFloat qw_t_y = tfxWideMul(q_w, t_y);
	tfxWideFloat qw_t_z = tfxWideMul(q_w, t_z);

	tfxWideFloat t_cross_x = tfxWideSub(tfxWideMul(q_y, t_z), tfxWideMul(q_z, t_y));
	tfxWideFloat t_cross_y = tfxWideSub(tfxWideMul(q_z, t_x), tfxWideMul(q_x, t_z));
	tfxWideFloat t_cross_z = tfxWideSub(tfxWideMul(q_x, t_y), tfxWideMul(q_y, t_x));

	*x = tfxWideAdd(tfxWideAdd(qv_x, qw_t_x), t_cross_x);
	*y = tfxWideAdd(tfxWideAdd(qv_y, qw_t_y), t_cross_y);
	*z = tfxWideAdd(tfxWideAdd(qv_z, qw_t_z), t_cross_z);
}

void tfx__wide_transform_packed_quaternion_vec3(tfxWideInt *quaternion, tfxWideFloat *x, tfxWideFloat *y, tfxWideFloat *z) {
	tfxWideFloat q_x, q_y, q_z, q_w;
	tfx__wide_unpack8bit(*quaternion, q_x, q_y, q_z, q_w);

	tfxWideFloat qv_x = *x;
	tfxWideFloat qv_y = *y;
	tfxWideFloat qv_z = *z;
	tfxWideFloat qv_w = tfxWideSetZero;

	tfxWideFloat two = tfxWideSetSingle(2.0f);

	tfxWideFloat t_x = tfxWideSub(tfxWideMul(q_y, qv_z), tfxWideMul(q_z, qv_y));
	t_x = tfxWideMul(t_x, two);
	tfxWideFloat t_y = tfxWideSub(tfxWideMul(q_z, qv_x), tfxWideMul(q_x, qv_z));
	t_y = tfxWideMul(t_y, two);
	tfxWideFloat t_z = tfxWideSub(tfxWideMul(q_x, qv_y), tfxWideMul(q_y, qv_x));
	t_z = tfxWideMul(t_z, two);

	tfxWideFloat qw_t_x = tfxWideMul(q_w, t_x);
	tfxWideFloat qw_t_y = tfxWideMul(q_w, t_y);
	tfxWideFloat qw_t_z = tfxWideMul(q_w, t_z);

	tfxWideFloat t_cross_x = tfxWideSub(tfxWideMul(q_y, t_z), tfxWideMul(q_z, t_y));
	tfxWideFloat t_cross_y = tfxWideSub(tfxWideMul(q_z, t_x), tfxWideMul(q_x, t_z));
	tfxWideFloat t_cross_z = tfxWideSub(tfxWideMul(q_x, t_y), tfxWideMul(q_y, t_x));

	*x = tfxWideAdd(tfxWideAdd(qv_x, qw_t_x), t_cross_x);
	*y = tfxWideAdd(tfxWideAdd(qv_y, qw_t_y), t_cross_y);
	*z = tfxWideAdd(tfxWideAdd(qv_z, qw_t_z), t_cross_z);
}

void tfx__wide_transform_packed_quaternion_vec2(tfxWideInt *quaternion, tfxWideFloat *x, tfxWideFloat *y) {
	tfxWideFloat q_x, q_y, q_z, q_w;
	tfx__wide_unpack8bit(*quaternion, q_x, q_y, q_z, q_w);

	tfxWideFloat s2 = tfxWideMul(q_z, q_z);
	tfxWideFloat c2 = tfxWideMul(q_w, q_w);
	tfxWideFloat sc = tfxWideMul(tfxWideSetSingle(2.f), tfxWideMul(q_z, q_w));

	tfxWideFloat rx = tfxWideSub(tfxWideMul(c2, *x), tfxWideAdd(tfxWideMul(sc, *y), tfxWideMul(s2, *x)));
	tfxWideFloat ry = tfxWideAdd(tfxWideMul(sc, *x), tfxWideSub(tfxWideMul(c2, *y), tfxWideMul(s2, *y)));

	*x = rx;
	*y = ry;
}

void tfx__wide_transform_quaternion_vec2(const tfx_quaternion_t *q, tfxWideFloat *x, tfxWideFloat *y) {
	tfxWideFloat c = tfxWideSetSingle(q->w);
	tfxWideFloat s = tfxWideSetSingle(q->z);

	tfxWideFloat s2 = tfxWideMul(s, s);
	tfxWideFloat c2 = tfxWideMul(c, c);
	tfxWideFloat sc = tfxWideMul(tfxWideSetSingle(2.f), tfxWideMul(s, c));

	tfxWideFloat rx = tfxWideSub(tfxWideMul(c2, *x), tfxWideAdd(tfxWideMul(sc, *y), tfxWideMul(s2, *x)));
	tfxWideFloat ry = tfxWideAdd(tfxWideMul(sc, *x), tfxWideSub(tfxWideMul(c2, *y), tfxWideMul(s2, *y)));

	*x = rx;
	*y = ry;
}

void tfx__transform_matrix4_vec3(const tfx_mat4_t *mat, tfxWideFloat *x, tfxWideFloat *y, tfxWideFloat *z) {
	tfxWideFloat xr = tfxWideMul(*x, tfxWideSetSingle(mat->v[0].c0));
	xr = tfxWideAdd(tfxWideMul(*y, tfxWideSetSingle(mat->v[0].c1)), xr);
	xr = tfxWideAdd(tfxWideMul(*z, tfxWideSetSingle(mat->v[0].c2)), xr);
	tfxWideFloat yr = tfxWideMul(*x, tfxWideSetSingle(mat->v[1].c0));
	yr = tfxWideAdd(tfxWideMul(*y, tfxWideSetSingle(mat->v[1].c1)), yr);
	yr = tfxWideAdd(tfxWideMul(*z, tfxWideSetSingle(mat->v[1].c2)), yr);
	tfxWideFloat zr = tfxWideMul(*x, tfxWideSetSingle(mat->v[2].c0));
	zr = tfxWideAdd(tfxWideMul(*y, tfxWideSetSingle(mat->v[2].c1)), zr);
	zr = tfxWideAdd(tfxWideMul(*z, tfxWideSetSingle(mat->v[2].c2)), zr);
	*x = xr;
	*y = yr;
	*z = zr;
}

void tfx__transform_matrix4_vec2(const tfx_mat4_t *mat, tfxWideFloat *x, tfxWideFloat *y) {
	tfxWideFloat xr = tfxWideMul(*x, tfxWideSetSingle(mat->v[0].c0));
	xr = tfxWideAdd(tfxWideMul(*y, tfxWideSetSingle(mat->v[1].c0)), xr);
	tfxWideFloat yr = tfxWideMul(*x, tfxWideSetSingle(mat->v[0].c1));
	yr = tfxWideAdd(tfxWideMul(*y, tfxWideSetSingle(mat->v[1].c1)), yr);
	*x = xr;
	*y = yr;
}

tfxU32 tfx__pack10bit_unsigned(tfx_vec3_t const *v) {
	tfx_vec3_t converted = *v * 511.f + 511.f;
	tfxUInt10bit result;
	result.pack = 0;
	result.data.x = (tfxU32)converted.z;
	result.data.y = (tfxU32)converted.y;
	result.data.z = (tfxU32)converted.x;
	result.data.w = 0;
	return result.pack;
}

tfxU32 tfx__pack16bit_sscaled(float x, float y, float max_value) {
	int16_t x_scaled = (int16_t)(x * 32767.0f / max_value);
	int16_t y_scaled = (int16_t)(y * 32767.0f / max_value);
	return ((tfxU64)x_scaled) | ((tfxU64)y_scaled << 16);
}

tfx_vec2_t tfx__unpack16bit_2sscaled(tfxU32 packed, float max_value) {
	int16_t x_scaled = (int16_t)(packed & 0xFFFF);
	int16_t y_scaled = (int16_t)((packed >> 16) & 0xFFFF);
	tfx_vec2_t unpacked;
	unpacked.x = (float)x_scaled * max_value / 32767.0f;
	unpacked.y = (float)y_scaled * max_value / 32767.0f;
	return unpacked;
}

tfxU32 tfx__pack8bit(tfx_vec3_t v) {
	union
	{
		signed char in[3];
		tfxU32 out;
	} u;

	v.x = tfx__Clamp(-1.f, 1.f, v.x) * 127.f;
	v.y = tfx__Clamp(-1.f, 1.f, v.y) * 127.f;
	v.z = tfx__Clamp(-1.f, 1.f, v.z) * 127.f;

	u.in[0] = (signed char)v.x;
	u.in[1] = (signed char)v.y;
	u.in[2] = (signed char)v.z;

	return u.out;
}

tfxU32 tfx__pack8bit_quaternion(tfx_quaternion_t q) {
	uint8_t x = static_cast<uint8_t>((q.x * 0.5f + 0.5f) * 255.0f);
	uint8_t y = static_cast<uint8_t>((q.y * 0.5f + 0.5f) * 255.0f);
	uint8_t z = static_cast<uint8_t>((q.z * 0.5f + 0.5f) * 255.0f);
	uint8_t w = static_cast<uint8_t>((q.w * 0.5f + 0.5f) * 255.0f);

	// Pack into a single 32-bit unsigned integer
	tfxU32 result = (w << 24) | (z << 16) | (y << 8) | x;

	return result;
}

tfx_vec2_t tfx__unpack16bit(tfxU32 in) {
	float one_div_32k = 1.f / 32767.f;

	tfx_vec2_t result;
	result.x = ((signed short)(in & 0x0000FFFF)) * one_div_32k;
	result.y = ((signed short)((in & 0xFFFF0000) >> 16)) * one_div_32k;

	return result;
}

tfxWideInt tfx__wide_pack16bit(tfxWideFloat v_x, tfxWideFloat v_y) {
	tfxWideFloat w32k = tfxWideSetSingle(32767.f);
	tfxWideInt bits16 = tfxWideSetSinglei(0xFFFF);
	tfxWideInt converted_y = tfxWideConverti(tfxWideMul(v_y, w32k));
	converted_y = tfxWideAndi(converted_y, bits16);
	converted_y = tfxWideShiftLeft(converted_y, 16);
	tfxWideInt converted_x = tfxWideConverti(tfxWideMul(v_x, w32k));
	converted_x = tfxWideAndi(converted_x, bits16);
	return tfxWideOri(converted_x, converted_y);
}

tfxWideInt tfx__wide_pack16bit_2sscaled(tfxWideFloat v_x, tfxWideFloat v_y, float max_value) {
	tfxWideFloat w32k = tfxWideSetSingle(32767.f / max_value);
	tfxWideInt bits16 = tfxWideSetSinglei(0xFFFF);
	tfxWideInt converted_y = tfxWideConverti(tfxWideMul(v_y, w32k));
	converted_y = tfxWideAndi(converted_y, bits16);
	converted_y = tfxWideShiftLeft(converted_y, 16);
	tfxWideInt converted_x = tfxWideConverti(tfxWideMul(v_x, w32k));
	converted_x = tfxWideAndi(converted_x, bits16);
	return tfxWideOri(converted_x, converted_y);
}

tfxWideInt tfx__wide_pack8bit_xyz(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z) {
	tfxWideFloat w127 = tfxWideSetSingle(127.f);
	tfxWideInt bits8 = tfxWideSetSinglei(0xFF);
	tfxWideInt converted_x = tfxWideConverti(tfxWideMul(v_x, w127));
	converted_x = tfxWideAndi(converted_x, bits8);
	//converted_x = tfxWideShiftLeft(converted_x, 0);
	tfxWideInt converted_y = tfxWideConverti(tfxWideMul(v_y, w127));
	converted_y = tfxWideAndi(converted_y, bits8);
	converted_y = tfxWideShiftLeft(converted_y, 8);
	tfxWideInt converted_z = tfxWideConverti(tfxWideMul(v_z, w127));
	converted_z = tfxWideAndi(converted_z, bits8);
	converted_z = tfxWideShiftLeft(converted_z, 16);
	return tfxWideOri(tfxWideOri(converted_x, converted_y), converted_z);
}

tfxWideInt tfx__wide_pack10bit_unsigned(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z) {
	tfxWideFloat w511 = tfxWideSetSingle(511.f);
	tfxWideInt bits10 = tfxWideSetSinglei(0x3FF);
	tfxWideInt converted_x = tfxWideConverti(tfxWideAdd(tfxWideMul(v_x, w511), w511));
	converted_x = tfxWideAndi(converted_x, bits10);
	converted_x = tfxWideShiftLeft(converted_x, 20);
	tfxWideInt converted_y = tfxWideConverti(tfxWideAdd(tfxWideMul(v_y, w511), w511));
	converted_y = tfxWideAndi(converted_y, bits10);
	converted_y = tfxWideShiftLeft(converted_y, 10);
	tfxWideInt converted_z = tfxWideConverti(tfxWideAdd(tfxWideMul(v_z, w511), w511));
	converted_z = tfxWideAndi(converted_z, bits10);
	return tfxWideOri(tfxWideOri(converted_x, converted_y), converted_z);
}

void tfx__wide_unpack10bit(tfxWideInt in, tfxWideFloat &x, tfxWideFloat &y, tfxWideFloat &z) {
	tfxWideInt w511 = tfxWideSetSinglei(511);
	x = tfxWideConvert(tfxWideSubi(tfxWideShiftRight(tfxWideAndi(in, tfxWideSetSinglei(0x3FF00000)), 20), w511));
	y = tfxWideConvert(tfxWideSubi(tfxWideShiftRight(tfxWideAndi(in, tfxWideSetSinglei(0x000FFC00)), 10), w511));
	z = tfxWideConvert(tfxWideSubi(tfxWideAndi(in, tfxWideSetSinglei(0x000003FF)), w511));
	x = tfxWideMul(x, one_div_511_wide);
	y = tfxWideMul(y, one_div_511_wide);
	z = tfxWideMul(z, one_div_511_wide);
}

void tfx__wide_unpack8bit(tfxWideInt in, tfxWideFloat &x, tfxWideFloat &y, tfxWideFloat &z, tfxWideFloat &w) {
	tfxWideInt mask_w = tfxWideSetSinglei(0xFF000000);
	tfxWideInt mask_z = tfxWideSetSinglei(0x00FF0000);
	tfxWideInt mask_y = tfxWideSetSinglei(0x0000FF00);
	tfxWideInt mask_x = tfxWideSetSinglei(0x000000FF);
	tfxWideInt w127 = tfxWideSetSinglei(127);

	w = tfxWideConvert(tfxWideSubi(tfxWideShiftRight(tfxWideAndi(in, mask_w), 24), w127));
	z = tfxWideConvert(tfxWideSubi(tfxWideShiftRight(tfxWideAndi(in, mask_z), 16), w127));
	y = tfxWideConvert(tfxWideSubi(tfxWideShiftRight(tfxWideAndi(in, mask_y), 8), w127));
	x = tfxWideConvert(tfxWideSubi(tfxWideAndi(in, mask_x), w127));

	x = tfxWideMul(x, one_div_127_wide);
	y = tfxWideMul(y, one_div_127_wide);
	z = tfxWideMul(z, one_div_127_wide);
	w = tfxWideMul(w, one_div_127_wide);
}

tfxWideFloat tfx__wide_unpack10bit_y(tfxWideInt in) {
	return tfxWideMul(tfxWideConvert(tfxWideSubi(tfxWideShiftRight(tfxWideAndi(in, tfxWideSetSinglei(0x000FFC00)), 10), tfxWideSetSinglei(511))), one_div_511_wide);
}

tfx_quaternion_t tfx__unpack8bit_quaternion(tfxU32 packed) {
	// Extract each component by shifting and masking
	uint8_t x = packed & 0xFF;
	uint8_t y = (packed >> 8) & 0xFF;
	uint8_t z = (packed >> 16) & 0xFF;
	uint8_t w = (packed >> 24) & 0xFF;

	// Convert from 0-255 back to -1 to 1 range
	tfx_quaternion_t q;
	q.x = (x / 255.0f) * 2.0f - 1.0f;
	q.y = (y / 255.0f) * 2.0f - 1.0f;
	q.z = (z / 255.0f) * 2.0f - 1.0f;
	q.w = (w / 255.0f) * 2.0f - 1.0f;

	// Return the unpacked quaternion
	return q;
}

float tfx__distance_2d(float fromx, float fromy, float tox, float toy) {
	float w = tox - fromx;
	float h = toy - fromy;
	return sqrtf(w * w + h * h);
}

tfx_vec2_t tfx_InterpolateVec2(float tween, tfx_vec2_t from, tfx_vec2_t to) {
	return to * tween + from * (1.f - tween);
}

tfx_vec3_t tfx_InterpolateVec3(float tween, tfx_vec3_t from, tfx_vec3_t to) {
	return to * tween + from * (1.f - tween);
}

tfx_rgba8_t tfx_InterpolateRGBA(float tween, tfx_rgba8_t from, tfx_rgba8_t to) {
	tfx_rgba8_t out;
	out.r = char((float)to.r * tween + (float)from.r * (1 - tween));
	out.g = char((float)to.g * tween + (float)from.g * (1 - tween));
	out.b = char((float)to.b * tween + (float)from.b * (1 - tween));
	out.a = char((float)to.a * tween + (float)from.a * (1 - tween));
	return out;
}

float tfx__gamma_correct(float color, float gamma) {
	return powf(color, gamma);
}

tfxWideFloat tfx_WideInterpolate(tfxWideFloat tween, tfxWideFloat *from, tfxWideFloat *to) {
	tfxWideFloat one_minus_tween = tfxWideSub(tfxWIDEONE, tween);
	tfxWideFloat to_lerp = tfxWideMul(*to, tween);
	tfxWideFloat from_lerp = tfxWideMul(*from, one_minus_tween);
	tfxWideFloat result = tfxWideAdd(from_lerp, to_lerp);
	return result;
}

float tfx_InterpolateFloat(float tween, float from, float to) {
	return to * tween + from * (1.f - tween);
}

void tfx__transform_2d(tfx_vec3_t *out_rotations, tfx_vec3_t *out_local_rotations, float *out_scale, tfx_vec3_t *out_position, tfx_vec3_t *out_local_position, tfx_vec3_t *out_translation, tfx_quaternion_t *out_q, tfx_effect_state_t *parent) {
	ToQuaternion2d(out_q, out_local_rotations->roll);
	*out_scale = parent->overal_scale;

	out_rotations->roll = parent->world_rotations.roll + out_local_rotations->roll;

	*out_q = *out_q * parent->rotation;
	tfx_vec2_t rotatevec = RotateVectorQuaternion2d(&parent->rotation, tfx_vec2_t(out_local_position->x + out_translation->x, out_local_position->y + out_translation->y));

	*out_position = parent->world_position.xy() + rotatevec * parent->overal_scale;
}

void tfx__transform_3d(tfx_vec3_t *out_rotations, tfx_vec3_t *out_local_rotations, float *out_scale, tfx_vec3_t *out_position, tfx_vec3_t *out_local_position, tfx_vec3_t *out_translation, tfx_quaternion_t *out_q, const tfx_effect_state_t *parent) {
	*out_q = EulerToQuaternion(out_local_rotations->pitch, out_local_rotations->yaw, out_local_rotations->roll);
	*out_scale = parent->overal_scale;

	*out_rotations = parent->world_rotations + *out_local_rotations;

	*out_q = *out_q * parent->rotation;
	tfx_vec3_t translated_vec = *out_local_position + *out_translation;
	tfx_vec3_t rotatevec = RotateVectorQuaternion(&parent->rotation, translated_vec);

	*out_position = parent->world_position + rotatevec;
}
//-------------------------------------------------
//--New transform_3d particle functions for SoA data--
//--------------------------2d---------------------
void tfx__transform_particle_position(const float local_position_x, const float local_position_y, const float roll, tfx_vec2_t *world_position, float *world_rotations) {
	world_position->x = local_position_x;
	world_position->y = local_position_y;
	*world_rotations = roll;
}

int tfx_FormatString(char *buf, size_t buf_size, const char *fmt, va_list args) {
	int w = vsnprintf(buf, buf_size, fmt, args);
	if (buf == NULL)
		return w;
	if (w == -1 || w >= (int)buf_size)
		w = (int)buf_size - 1;
	buf[w] = 0;
	return w;
}

void tfx_str_t::Appendv(const char *format, va_list args) {
	va_list args_copy;
	va_copy(args_copy, args);

	int len = tfx_FormatString(NULL, 0, format, args);         // FIXME-OPT: could do a first pass write attempt, likely successful on first pass.
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
	tfx_FormatString(&data[write_off - 1], (size_t)len + 1, format, args_copy);
	va_end(args_copy);

}

int tfx_str_t::Find(const char *needle) {
	tfx_str_t compare = needle;
	tfx_str_t lower = Lower();
	compare = compare.Lower();
	if (compare.Length() > Length()) return -1;
	tfxU32 pos = 0;
	while (compare.Length() + pos <= Length()) {
		if (strncmp(lower.data + pos, compare.data, compare.Length()) == 0) {
			return pos;
		}
		++pos;
	}
	return -1;
}

tfx_str_t tfx_str_t::Lower() {
	tfx_str_t convert = *this;
	for (auto &c : convert) {
		c = tolower(c);
	}
	return convert;
}

void tfx_str_t::AddLine(const char *format, ...) {
	va_list args;
	va_start(args, format);
	Appendv(format, args);
	va_end(args);
	Append('\n');
}

bool tfx_str_t::SaveToFile(const char *file_name) {
	FILE *file = tfx__open_file(file_name, "wb");
	if (!file)
		return false;

	tfxU32 l = Length();
	if (fwrite(data, sizeof(char), Length(), file) != Length()) {
		fclose(file);
		return false;
	}

	fclose(file);
	return true;
}

const bool tfx_str_t::IsInt() const {
	if (!Length()) return false;
	for (auto c : *this) {
		if (c >= '0' && c <= '9');
		else {
			return false;
		}
	}
	return true;
}

const bool tfx_str_t::IsFloat() const {
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

void tfx_str_t::Setf(const char *format, ...) {
	Clear();
	va_list args;
	va_start(args, format);

	va_list args_copy;
	va_copy(args_copy, args);

	int len = tfx_FormatString(NULL, 0, format, args);         // FIXME-OPT: could do a first pass write attempt, likely successful on first pass.
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
	tfx_FormatString(&strbuffer()[write_off - 1], (size_t)len + 1, format, args_copy);
	va_end(args_copy);

	va_end(args);
}

void tfx_str_t::Appendf(const char *format, ...) {
	va_list args;
	va_start(args, format);

	va_list args_copy;
	va_copy(args_copy, args);

	int len = tfx_FormatString(NULL, 0, format, args);         // FIXME-OPT: could do a first pass write attempt, likely successful on first pass.
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
	tfx_FormatString(&strbuffer()[write_off - 1], (size_t)len + 1, format, args_copy);
	va_end(args_copy);

	va_end(args);
}

tfx_str512_t tfx_stream_t::ReadLine() {
	tfx_str512_t line;
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

tfx_package_t tfx__create_package(const char *file_path) {
	tfx_package_t package;
	package.header.magic_number = tfxMAGIC_NUMBER;
	package.header.flags = 0;
	package.header.offset_to_inventory = sizeof(tfx_package_header_t);
	package.header.file_version = tfxFILE_VERSION;

	package.inventory.magic_number = tfxMAGIC_NUMBER_INVENTORY;
	package.inventory.entry_count = 0;

	package.file_path = file_path;
	return package;
}

bool tfx__validate_package(tfx_package_t *package) {
	if (package->header.magic_number != tfxMAGIC_NUMBER) return false;            //Package hasn't been initialised

	if (package->flags & tfxPackageFlags_loaded_from_memory) {
		return true;
	}

	FILE *file = tfx__open_file(package->file_path.c_str(), "rb");
	if (!file) {
		return false;
	}

	if (file == NULL || tfx__fseek(file, 0, SEEK_END)) {
		return false;
	}

	tfxU64 length = tfx__ftell(file);
	rewind(file);
	if (length == -1 || length >= SIZE_MAX) {
		fclose(file);
		return false;
	}

	if (length != package->file_size) return false;                            //The file on disk is no longer the same size as the package file size since it was loaded

	//Everything seems ok
	fclose(file);
	return true;
}

void tfx_package_entry_info_t::FreeData() {
	data.FreeAll();
}

tfx_package_t::~tfx_package_t() {
	tfx__free_package(this);
}

tfx_package_entry_info_t *tfx__get_package_file(tfx_package_t *package, const char *name) {
	if (!package->inventory.entries.ValidName(name)) {
		return nullptr;                                        //File not found in inventory
	}
	TFX_ASSERT(tfx__validate_package(package));                    //The file on disk has changed since the package was loaded! Maybe this should return null instead?
	//Also: function call in assert, sort this out!
	tfx_package_entry_info_t *entry = &package->inventory.entries.At(name);
	if (entry->data.Size() != entry->file_size) {
		if (!(package->flags & tfxPackageFlags_loaded_from_memory)) {
			FILE *file = tfx__open_file(package->file_path.c_str(), "rb");
			TFX_ASSERT(file);        //couldn't open the file!
			tfx__fseek(file, entry->offset_from_start_of_file, SEEK_SET);
			entry->data.Resize(entry->file_size);
			fread(entry->data.data, 1, entry->file_size, file);
			fclose(file);
		}
		else {
			entry->data.Resize(entry->file_size);
			char *point_in_file = package->file_data.data + entry->offset_from_start_of_file;
			TFX_ASSERT(entry->offset_from_start_of_file + entry->file_size < package->file_size);    //Invalid entry/package data, make sure the package file is not corrupt
			memcpy(entry->data.data, point_in_file, entry->file_size);
		}
	}
	return entry;
}

bool tfx__file_exists_in_package(tfx_package_t *package, const char *file_name) {
	if (package->inventory.entries.ValidName(file_name)) {
		return true;
	}
	return false;
}

void tfx__add_entry_to_package(tfx_package_t *package, tfx_package_entry_info_t file) {
	package->inventory.entries.Insert(file.file_name, file);
	package->inventory.entry_count++;
}

void tfx__add_file_to_package(tfx_package_t *package, const char *file_name, tfx_stream_t *data) {
	tfx_package_entry_info_t entry;
	entry.file_name = file_name;
	entry.data = *data;
	entry.file_size = data->size;

	package->inventory.entries.Insert(entry.file_name, entry);
	package->inventory.entry_count++;
}

void tfx__free_package(tfx_package_t *package) {
	for (auto &entry : package->inventory.entries.data) {
		entry.data.FreeAll();
	}
	package->inventory.entries.data.free_all();
	package->inventory.entries.map.free_all();
	package->file_data.FreeAll();
}

void tfx__copy_stream(tfx_stream_t *dst, tfx_stream_t *src) {
	dst->FreeAll();
	dst->Resize(src->size);
	memcpy(dst->data, src->data, src->size);
}

void tfx__copy_stream_to_string(tfx_str_t *dst, tfx_stream_t *src) {
	dst->free_all();
	dst->resize((tfxU32)src->size);
	memcpy(dst->data, src->data, src->size);
}

void tfx__copy_string_to_stream(tfx_stream_t *dst, tfx_str_t *src) {
	dst->FreeAll();
	dst->Resize(src->current_size);
	memcpy(dst->data, src->data, src->current_size);
}

void tfx__copy_data_to_stream(tfx_stream_t *dst, const void *src, tfxU64 size) {
	dst->FreeAll();
	dst->Resize(size);
	memcpy(dst->data, src, size);
}

// Reads the whole file on disk into memory and returns the pointer
tfx_stream_t tfx__read_entire_file(const char *file_name, bool terminate) {
	tfx_stream_t buffer;
	FILE *file = tfx__open_file(file_name, "rb");
	if (!file) {
		return buffer;
	}

	// file invalid? fseek() fail?
	if (file == NULL || fseek(file, 0, SEEK_END)) {
		fclose(file);
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
		buffer.Resize((tfxU64)length + 1);
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

bool tfx__save_package_disk(tfx_package_t *package) {
	if (!package->file_path.Length()) return false;                                            //Package must have a file path
	if (package->header.magic_number != tfxMAGIC_NUMBER) return false;                        //Header of package must contain correct magic number. Use tfx__create_package to correctly initialise a package.
	if (package->inventory.magic_number != tfxMAGIC_NUMBER_INVENTORY) return false;            //Inventory of package must contain correct magic number

	FILE *file = tfx__open_file(package->file_path.c_str(), "wb");
	if (!file)
		return false;

	//Calculate the offset to the inventory which is stored at the end of the file after the contents
	tfxU64 inventory_offset = sizeof(tfx_package_header_t);
	for (auto &entry : package->inventory.entries.data) {
		entry.offset_from_start_of_file = inventory_offset;
		entry.file_size = entry.data.Size();
		inventory_offset += entry.data.Size();
	}

	//Sanity check, make sure that the entry count is the correct value
	package->inventory.entry_count = package->inventory.entries.Size();

	//Write the header, updating the inventory offset before hand
	package->header.offset_to_inventory = inventory_offset;
	fwrite((char *)&package->header, sizeof(char), sizeof(tfx_package_header_t), file);

	//Write the file contents
	for (auto &entry : package->inventory.entries.data) {
		fwrite(entry.data.data, sizeof(char), entry.data.Size(), file);
	}

	//Write the inventory
	fwrite((char *)&package->inventory.magic_number, sizeof(char), sizeof(tfxU32), file);
	fwrite((char *)&package->inventory.entry_count, sizeof(char), sizeof(tfxU32), file);
	for (auto &entry : package->inventory.entries.data) {
		fwrite((char *)&entry.file_name.current_size, sizeof(char), sizeof(tfxU32), file);
		fwrite(entry.file_name.c_str(), 1, entry.file_name.current_size, file);
		fwrite((char *)&entry.file_size, 1, sizeof(tfxU64), file);
		fwrite((char *)&entry.offset_from_start_of_file, sizeof(char), sizeof(tfxU64), file);
	}

	fclose(file);
	return true;
}

tfx_stream_t tfx__save_package_memory(tfx_package_t *package) {
	if (!package->file_path.Length()) return false;                                            //Package must have a file path
	if (package->header.magic_number != tfxMAGIC_NUMBER) return false;                        //Header of package must contain correct magic number. tfx__create_package to correctly initialise a package.
	if (package->inventory.magic_number != tfxMAGIC_NUMBER_INVENTORY) return false;            //Inventory of package must contain correct magic number

	//char *file = (char*)malloc(tfx__get_package_size(package));
	tfx_stream_t file(tfx__get_package_size(package));
	if (!file.Size())
		return file;

	//Calculate the offset to the inventory which is stored at the end of the file after the contents
	tfxU64 inventory_offset = sizeof(tfx_package_header_t);
	for (auto &entry : package->inventory.entries.data) {
		entry.offset_from_start_of_file = inventory_offset;
		entry.file_size = entry.data.Size();
		inventory_offset += entry.data.Size();
	}

	//Sanity check, make sure that the entry count is the correct value
	package->inventory.entry_count = package->inventory.entries.Size();

	//Write the header, updating the inventory offset before hand
	package->header.offset_to_inventory = inventory_offset;
	file.Write(&package->header, sizeof(tfx_package_header_t));

	//Write the file contents
	for (auto &entry : package->inventory.entries.data) {
		//fwrite(entry.data.data, sizeof(char), entry.data.Size(), file);
		file.Write(entry.data.data, entry.data.Size());
	}

	//Write the inventory
	file.Write(&package->inventory.magic_number, sizeof(tfxU32));
	file.Write(&package->inventory.entry_count, sizeof(tfxU32));
	for (auto &entry : package->inventory.entries.data) {
		file.Write(&entry.file_name.current_size, sizeof(tfxU32));
		file.Write(entry.file_name.data, entry.file_name.current_size);
		file.Write(&entry.file_size, sizeof(tfxU64));
		file.Write(&entry.offset_from_start_of_file, sizeof(tfxU64));
	}

	return file;
}

tfxU64 tfx__get_package_size(tfx_package_t *package) {
	tfxU64 space = 0;
	space += sizeof(tfx_package_header_t);

	//Write the file contents
	for (auto &entry : package->inventory.entries.data) {
		space += entry.data.Size();
	}

	//Write the inventory
	space += sizeof(tfxU32);
	space += sizeof(tfxU32);
	for (auto &entry : package->inventory.entries.data) {
		space += sizeof(tfxU32);
		space += entry.file_name.current_size;
		space += sizeof(tfxU64);
		space += sizeof(tfxU64);
	}

	return space;
}

tfxErrorFlags tfx__load_package_file(const char *file_name, tfx_package_t *package) {

	package->file_data = tfx__read_entire_file(file_name);
	if (package->file_data.Size() == 0)
		return tfxErrorCode_unable_to_read_file;            //the file size is smaller then the expected header size

	package->file_size = package->file_data.Size();

	if (package->file_size < sizeof(tfx_package_header_t))
		return tfxErrorCode_wrong_file_size;                //the file size is smaller then the expected header size

	package->file_data.Read((char *)&package->header, sizeof(tfx_package_header_t));

	if (package->header.magic_number != tfxMAGIC_NUMBER)
		return tfxErrorCode_invalid_format;                //The header doesn't not contain the expected magic number "TFX!", incorrect file format;

	if (package->header.offset_to_inventory > package->file_size)
		return tfxErrorCode_no_inventory;                //The offset to the inventory is beyond the size of the file

	package->file_data.Seek(package->header.offset_to_inventory);
	package->file_data.Read((char *)&package->inventory.magic_number, sizeof(tfxU32));

	if (package->inventory.magic_number != tfxMAGIC_NUMBER_INVENTORY)
		return tfxErrorCode_invalid_inventory;            //The value at the inventory offset does not equal the expected magic number "INV!"

	package->file_data.Read((char *)&package->inventory.entry_count, sizeof(tfxU32));
	for (int i = 0; i != package->inventory.entry_count; ++i) {
		tfx_package_entry_info_t entry;
		tfxU32 file_name_size;
		package->file_data.Read((char *)&file_name_size, sizeof(tfxU32));
		entry.file_name.resize(file_name_size);
		package->file_data.Read(entry.file_name.data, file_name_size);
		package->file_data.Read((char *)&entry.file_size, sizeof(tfxU64));
		package->file_data.Read((char *)&entry.offset_from_start_of_file, sizeof(tfxU64));
		package->inventory.entries.Insert(entry.file_name, entry);
	}

	package->file_path = file_name;

	return 0;
}

tfxErrorFlags tfx__load_package_stream(tfx_stream_t *stream, tfx_package_t *package) {
	//Note: tfx_stream_t does not copy the memory, only the pointer, so if you FreeAll on the stream you pass in it will also free the file_data here as well
	package->flags |= tfxPackageFlags_loaded_from_memory;
	package->file_data = *stream;
	package->file_data.Seek(0);
	if (package->file_data.Size() == 0)
		return tfxErrorCode_unable_to_read_file;            //the file size is smaller then the expected header size

	package->file_size = package->file_data.Size();

	if (package->file_size < sizeof(tfx_package_header_t))
		return tfxErrorCode_wrong_file_size;                //the file size is smaller then the expected header size

	package->file_data.Read((char *)&package->header, sizeof(tfx_package_header_t));

	if (package->header.magic_number != tfxMAGIC_NUMBER)
		return tfxErrorCode_invalid_format;                //The header doesn't not contain the expected magic number "TFX!", incorrect file format;

	if (package->header.offset_to_inventory > package->file_size)
		return tfxErrorCode_no_inventory;                //The offset to the inventory is beyond the size of the file

	package->file_data.Seek(package->header.offset_to_inventory);
	package->file_data.Read((char *)&package->inventory.magic_number, sizeof(tfxU32));

	if (package->inventory.magic_number != tfxMAGIC_NUMBER_INVENTORY)
		return tfxErrorCode_invalid_inventory;            //The value at the inventory offset does not equal the expected magic number "INV!"

	package->file_data.Read((char *)&package->inventory.entry_count, sizeof(tfxU32));
	for (int i = 0; i != package->inventory.entry_count; ++i) {
		tfx_package_entry_info_t entry;
		tfxU32 file_name_size;
		package->file_data.Read((char *)&file_name_size, sizeof(tfxU32));
		entry.file_name.resize(file_name_size);
		package->file_data.Read(entry.file_name.data, file_name_size);
		package->file_data.Read((char *)&entry.file_size, sizeof(tfxU64));
		package->file_data.Read((char *)&entry.offset_from_start_of_file, sizeof(tfxU64));
		package->inventory.entries.Insert(entry.file_name, entry);
	}

	return 0;
}

tfx_effect_emitter_t::~tfx_effect_emitter_t() {
}

void UpdateEffectMaxLife(tfx_effect_emitter_t *effect) {
	tfx_effect_emitter_info_t *info = GetEffectInfo(effect);
	info->max_life = GetMaxLife(effect);
	GetEffectGraphByType(effect, tfxOvertime_red)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_green)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_blue)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_blendfactor)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_intensity)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_alpha_sharpness)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_curved_alpha)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_velocity)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_width)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_height)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_weight)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_roll_spin)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_stretch)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_roll_spin)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_pitch_spin)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_yaw_spin)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_velocity_turbulance)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_direction_turbulance)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_velocity_adjuster)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_direction)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxOvertime_motion_randomness)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxFactor_life)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxFactor_velocity)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxFactor_size)->lookup.life = info->max_life;
	GetEffectGraphByType(effect, tfxFactor_intensity)->lookup.life = info->max_life;
}

bool IsFiniteEmitter(tfx_effect_emitter_t *emitter) {
	if (emitter->property_flags & tfxEmitterPropertyFlags_single && GetEffectProperties(emitter)->single_shot_limit == 0) {
		return false;
	}
	float qty = GetGraphLastValue(&emitter->library->emitter_attributes[emitter->emitter_attributes].base.amount) + GetGraphLastValue(&emitter->library->emitter_attributes[emitter->emitter_attributes].variation.amount);
	if (!(emitter->property_flags & tfxEmitterPropertyFlags_single) && qty > 0) {
		return false;
	}
	return true;
}

bool IsFiniteEffect(tfx_effect_emitter_t *effect) {
	TFX_ASSERT(effect->type == tfxEffectType);
	for (auto &e : GetEffectInfo(effect)->sub_effectors) {
		float qty = GetGraphLastValue(&e.library->emitter_attributes[e.emitter_attributes].base.amount) + GetGraphLastValue(&e.library->emitter_attributes[e.emitter_attributes].variation.amount);
		if (e.path_attributes != tfxINVALID && GetEffectProperties(&e)->emission_type == tfxPath) {
			tfx_emitter_path_t *path = GetEmitterPath(&e);
			if (path->rotation_range > 0 && path->maximum_paths > 0) {
				continue;
			}
		}
		if (!(e.property_flags & tfxEmitterPropertyFlags_single) && qty > 0) {
			return false;
		}
		else if (e.property_flags & tfxEmitterPropertyFlags_single && GetEffectProperties(&e)->single_shot_limit == 0) {
			return false;
		}
	}
	return true;
}

void FlagEffectAs3D(tfx_effect_emitter_t *effect, bool flag) {
	if (flag) {
		effect->property_flags |= tfxEmitterPropertyFlags_effect_is_3d;
	}
	else {
		effect->property_flags &= ~tfxEmitterPropertyFlags_effect_is_3d;
	}
	for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
		FlagEffectAs3D(&sub, flag);
	}
}

void FlagEffectsAs3D(tfx_library_t *library) {
	for (tfx_effect_emitter_t &effect : library->effects) {
		if (effect.type == tfxEffectType) {
			FlagEffectAs3D(&effect, Is3DEffect(&effect));
		}
		else {
			for (tfx_effect_emitter_t &sub : GetEffectInfo(&effect)->sub_effectors) {
				FlagEffectAs3D(&effect, Is3DEffect(&sub));
			}
		}
	}
}

bool Is3DEffect(tfx_effect_emitter_t *effect) {
	return effect->property_flags & tfxEmitterPropertyFlags_effect_is_3d;
}

bool IsOrderedEffect(tfx_effect_emitter_t *effect) {
	tfxEffectPropertyFlags ordered_flags = tfxEffectPropertyFlags_age_order | tfxEffectPropertyFlags_depth_draw_order;
	return (effect->effect_flags & ordered_flags) > 0;
}

bool IsOrderedEffectState(tfx_effect_state_t *effect) {
	tfxEffectPropertyFlags ordered_flags = tfxEffectPropertyFlags_age_order | tfxEffectPropertyFlags_depth_draw_order;
	return (effect->effect_flags & ordered_flags) > 0;
}

tfx_particle_manager_mode GetRequiredParticleManagerMode(tfx_effect_emitter_t *effect) {
	if (effect->type == tfxEffectType) {
		if (effect->effect_flags & tfxEffectPropertyFlags_guaranteed_order && effect->effect_flags & tfxEffectPropertyFlags_depth_draw_order) {
			return tfxParticleManagerMode_ordered_by_depth_guaranteed;
		}
		else if (effect->effect_flags & tfxEffectPropertyFlags_depth_draw_order) {
			return tfxParticleManagerMode_ordered_by_depth;
		}
		else if (effect->effect_flags & tfxEffectPropertyFlags_age_order) {
			return tfxParticleManagerMode_ordered_by_age;
		}
	}
	else if (effect->type == tfxStage) {
		tfx_particle_manager_mode result = tfxParticleManagerMode_unordered;
		for (auto &sub_effect : GetEffectInfo(effect)->sub_effectors) {
			if (sub_effect.effect_flags & tfxEffectPropertyFlags_guaranteed_order && sub_effect.effect_flags & tfxEffectPropertyFlags_depth_draw_order) {
				return tfxParticleManagerMode_ordered_by_depth_guaranteed;
			}
			else if (sub_effect.effect_flags & tfxEffectPropertyFlags_depth_draw_order) {
				result = tfxParticleManagerMode_ordered_by_depth;
			}
			else if (result != tfxParticleManagerMode_ordered_by_depth && sub_effect.effect_flags & tfxEffectPropertyFlags_age_order) {
				result = tfxParticleManagerMode_ordered_by_age;
			}
		}
		return result;
	}
	return tfxParticleManagerMode_unordered;
}

tfx_preview_camera_settings_t *GetEffectCameraSettings(tfx_effect_emitter_t *effect) {
	return &effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings];
}

float GetEffectLoopLength(tfx_effect_emitter_t *effect) {
	return GetEffectProperties(effect)->loop_length;
}

float GetEffectHighestLoopLength(tfx_effect_emitter_t *effect) {
	float loop_length = GetEffectLoopLength(effect);
	for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
		loop_length = tfxMax(GetEffectHighestLoopLength(&sub), loop_length);
	}
	return loop_length;
}

tfx_effect_emitter_t *AddEmitterToEffect(tfx_effect_emitter_t *effect, tfx_effect_emitter_t *emitter) {
	TFX_ASSERT(GetEffectInfo(emitter)->name.Length());                //Emitter must have a name so that a hash can be generated
	emitter->type = tfx_effect_emitter_type::tfxEmitterType;
	emitter->library = effect->library;
	GetEffectInfo(emitter)->uid = ++effect->library->uid;
	GetEffectInfo(effect)->sub_effectors.push_back(*emitter);
	UpdateLibraryEffectPaths(effect->library);
	ReIndexEffect(effect);
	return &GetEffectInfo(effect)->sub_effectors.back();
}

tfx_effect_emitter_t *AddEffectToEmitter(tfx_effect_emitter_t *emitter, tfx_effect_emitter_t *effect) {
	TFX_ASSERT(GetEffectInfo(effect)->name.Length());                //Effect must have a name so that a hash can be generated
	effect->type = tfx_effect_emitter_type::tfxEffectType;
	effect->library = emitter->library;
	effect->parent = emitter;
	GetEffectInfo(effect)->uid = ++emitter->library->uid;
	GetEffectInfo(emitter)->sub_effectors.push_back(*effect);
	UpdateLibraryEffectPaths(emitter->library);
	ReIndexEffect(emitter);
	return &GetEffectInfo(emitter)->sub_effectors.back();
}

tfx_effect_emitter_t *AddEffect(tfx_effect_emitter_t *e) {
	tfx_effect_emitter_t new_effect;
	new_effect.library = e->library;
	GetEffectInfo(&new_effect)->uid = ++e->library->uid;
	new_effect.type = tfx_effect_emitter_type::tfxEffectType;
	GetEffectInfo(&new_effect)->name = "New Effect";
	GetEffectInfo(e)->sub_effectors.push_back(new_effect);
	UpdateLibraryEffectPaths(e->library);
	ReIndexEffect(e);
	return &GetEffectInfo(e)->sub_effectors.back();
}

tfxU32 CountAllEffects(tfx_effect_emitter_t *effect, tfxU32 amount) {
	for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
		amount = CountAllEffects(&sub, amount);
	}
	return ++amount;
}

int GetEffectDepth(tfx_effect_emitter_t *e) {
	tfx_effect_emitter_t *current_parent = e->parent;
	int depth = 0;
	while (current_parent) {
		if (current_parent->type == tfxEmitterType) {
			depth++;
		}
		current_parent = current_parent->parent;
	}
	return depth;
}

template<typename T>
void ClearWrapBit(T* instance, tfx_sprite_data_t *sprite_data) {
    for (int i = 0; i != sprite_data->normal.frame_count; ++i) {
        for (tfxEachLayer) {
            for (int j = SpriteDataIndexOffset(sprite_data, i, layer); j != SpriteDataEndIndex(sprite_data, i, layer); ++j) {
                if (instance[j].captured_index == tfxINVALID) continue;
                tfxU32 wrap_bit = 0x80000000;
                instance[j].captured_index &= ~wrap_bit;
            }
        }
    }
}

template<typename T>
void SpriteDataOffsetCapturedIndexes(T* instance, tfx_sprite_data_t *sprite_data, tfxU32 previous_frame, tfxU32 current_frame) {
    tfxU32 layer_offset = 0;
    for (tfxEachLayer) {
        for (int j = SpriteDataIndexOffset(sprite_data, current_frame, layer); j != SpriteDataEndIndex(sprite_data, current_frame, layer); ++j) {
            if (instance[j].captured_index == tfxINVALID) {
                //If the captured index of the sprite is invalid then this is the first frame of the sprite and therefore nothing to interpolate with in the previous
                //frame
                continue;
            } else if (instance[j].captured_index != tfxINVALID && sprite_data->real_time_sprites.uid[j].age == 0) {
                //This deals with wrapped single instance_data where captured index is not invalid yet it's age is 0. Usually only captured indexes that have an age of 0
                //are invalid because there is no previous frame of the sprite to interpolate with but we do want that to happen with a wrapped sprite. So it means
                //that we have to manually find the wrapped sprite in the previous frame.
                for (int k = SpriteDataIndexOffset(sprite_data, previous_frame, layer); k != SpriteDataEndIndex(sprite_data, previous_frame, layer); ++k) {
                    if (sprite_data->real_time_sprites.uid[k].uid == sprite_data->real_time_sprites.uid[j].uid) {
                        tfxU32 wrap_bit = instance[j].captured_index & 0x80000000;
                        instance[j].captured_index = k;
                        TFX_ASSERT(sprite_data->real_time_sprites.uid[instance[j].captured_index].uid == sprite_data->real_time_sprites.uid[j].uid);
                        instance[j].captured_index |= wrap_bit;
                        break;
                    }
                }
            } else {
                //Add the index offset of the frame to realign the captured index of the sprite.
                tfxU32 wrap_bit = instance[j].captured_index & 0x80000000;
                instance[j].captured_index = (instance[j].captured_index & 0x0FFFFFFF) + sprite_data->normal.frame_meta[previous_frame].index_offset[layer];
                TFX_ASSERT(sprite_data->real_time_sprites.uid[instance[j].captured_index].uid == sprite_data->real_time_sprites.uid[j].uid);
                instance[j].captured_index |= wrap_bit;
            }
        }
        if (layer > 0) {
            layer_offset += sprite_data->normal.frame_meta[previous_frame].cumulative_offset[layer];
        }
    }
}

template<typename T>
void LinkSpriteDataCapturedIndexes(T* instance, int entry_frame, tfx_sprite_data_t *sprite_data) {
    tfx_sprite_data_soa_t &c_sprites = sprite_data->compressed_sprites;
    int frame = entry_frame - 1;
    int frame_pair[2];
    frame_pair[0] = entry_frame;
    frame_pair[1] = frame < 0 ? sprite_data->compressed.frame_count - 1 : frame;
    for (tfxEachLayer) {
        for (int i = sprite_data->compressed.frame_meta[frame_pair[0]].index_offset[layer]; i != sprite_data->compressed.frame_meta[frame_pair[0]].index_offset[layer] + sprite_data->compressed.frame_meta[frame_pair[0]].sprite_count[layer]; ++i) {
            if (instance[i].captured_index != tfxINVALID) {
                int age_diff = 0x7FFFFFFF;
                bool captured_found = false;
                for (int j = sprite_data->compressed.frame_meta[frame_pair[1]].index_offset[layer]; j != sprite_data->compressed.frame_meta[frame_pair[1]].index_offset[layer] + sprite_data->compressed.frame_meta[frame_pair[1]].sprite_count[layer]; ++j) {
                    //Find particles of the same uid and only if they haven't already been linked up to another particle already (the first bit in age will be flagged in this case)
                    if (c_sprites.uid[j].uid == c_sprites.uid[i].uid && !(c_sprites.uid[j].age & 0x80000000)) {
                        int diff = (int)(0x7FFFFFFF & c_sprites.uid[i].age) - (int)(0x7FFFFFFF & c_sprites.uid[j].age);
                        if (diff < 0) continue;
                        age_diff = diff < age_diff ? diff : age_diff;
                        //We link up the captured index to the particle in the previous frame with the smallest difference in age
                        //Make sure that the particle we're linking to is younger than the current particle
                        //instance[i].captured_index = diff < 0 ? tfxINVALID : instance[i].captured_index;
                        captured_found = false;
                        instance[i].captured_index = age_diff == diff ? j : instance[i].captured_index;
						//tfxPrint("%i (%i)) UID: %u: CI: %i, SI: %i - CIAge: %i, SAge: %i - Age Diff: %u, Diff: %u, Cap.index: %u, CPosy: %.2f SPosy: %.2f, H: %u",
						//	entry_frame, sprite_data->compressed.frame_meta[frame_pair[0]].sprite_count[layer], c_sprites.uid[i].uid, i, j, c_sprites.uid[i].age, c_sprites.uid[j].age, age_diff, diff, instance[i].captured_index, instance[i].position.y, instance[j].position.y, tfxU32(instance[j].size_handle.packed >> 32));
                        if (age_diff < 2) {    //We can just break if the age is less than 2 but not if we found and older particle. It's possible to find an older particle if the animation has been looped
                            //If the compression is high this won't be hit because the distance in time between the compressed frames will be high
                            //tfxPrint("\t Linked %i to %i: uid, %u, captured index: %u", i, j, c_sprites.uid[j].uid, instance[i].captured_index);
                            break;
                        }
                    }
                }
                if (captured_found) {
                    c_sprites.uid[instance[i].captured_index].age |= 0x80000000;
                }
            } else {
				//tfxPrint("%i (%i)) UID: %u: CI: %i, Cap: %u, H: %u", entry_frame, sprite_data->compressed.frame_meta[frame_pair[0]].sprite_count[layer], c_sprites.uid[i].uid, i, instance[i].captured_index, tfxU32(instance[i].size_handle.packed >> 32));
            }
        }
    }
}


float GetEmissionDirection2d(tfx_particle_manager_t *pm, tfx_library_t *library, tfx_random_t *random, tfx_emitter_state_t &emitter, tfx_vec2_t local_position, tfx_vec2_t world_position) {
	//float (*effect_lookup_callback)(tfx_graph_t &graph, float age) = common.root_effect->lookup_mode == tfxPrecise ? LookupPrecise : LookupFast;
	float emission_angle = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].properties.emission_pitch, emitter.frame);
	float emission_angle_variation = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].properties.emission_range, emitter.frame);
	//----Emission
	float range = emission_angle_variation * .5f;
	float direction = 0;
	tfx_emission_type emission_type = library->emitter_properties[emitter.properties_index].emission_type;
	tfx_emission_direction emission_direction = library->emitter_properties[emitter.properties_index].emission_direction;

	if (emission_type == tfxPoint)
		return direction + emission_angle + RandomRange(random, -range, range);

	tfx_vec2_t tmp_position;
	if (emitter.handle.x + local_position.x == 0 && emitter.handle.y + local_position.y == 0)
		tmp_position = emitter.emitter_size.xy();
	else
		tmp_position = local_position + emitter.handle.xy();

	if (emission_direction == tfx_emission_direction::tfxOutwards) {

		tfx_vec2_t to_handle;

		if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position)
			to_handle = tmp_position;
		else
			to_handle = world_position - emitter.world_position.xy();

		direction = GetVectorAngle(to_handle.x, to_handle.y);

	}
	else if (emission_direction == tfx_emission_direction::tfxInwards) {

		tfx_vec2_t to_handle;

		if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position)
			to_handle = -tmp_position;

		else
			to_handle = emitter.world_position.xy() - world_position;

		direction = GetVectorAngle(to_handle.x, to_handle.y);

	}
	else if (emission_direction == tfx_emission_direction::tfxBothways) {

		//todo: replace these if statements
		if (emitter.emission_alternator) {

			tfx_vec2_t to_handle;

			if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = (tmp_position);
			else
				to_handle = world_position - emitter.world_position.xy();

			direction = GetVectorAngle(to_handle.x, to_handle.y);

		}
		else {

			tfx_vec2_t to_handle;

			if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = -tmp_position;
			else
				to_handle = (emitter.world_position.xy() - world_position);

			direction = GetVectorAngle(to_handle.x, to_handle.y);

		}

		emitter.emission_alternator = !emitter.emission_alternator;
	}

	if (isnan(direction)) {
		direction = 0.f;
	}
	return direction + emission_angle + RandomRange(random, -range, range);
}

tfx_vec3_t RandomVectorInCone(tfx_random_t *random, tfx_vec3_t cone_direction, float cone_angle) {
	// Convert cone angle to radians

	// Calculate the minimum z value for the cone
	float min_z = cosf(cone_angle);

	// Randomly sample z in [min_z, 1]
	float z = RandomRange(random, min_z, 1.0f);

	// Randomly sample ϕ in [0, 2π)
	float phi = RandomRange(random, 0.0f, 2.0f * tfxPI);

	// Calculate the corresponding x and y for the random point on the unit sphere
	float sqrt_one_minus_z_squared = sqrtf(1.0f - z * z);
	float x = sqrt_one_minus_z_squared * cosf(phi);
	float y = sqrt_one_minus_z_squared * sinf(phi);

	// Create the random vector in the cone's local space
	tfx_vec3_t random_vector = { x, y, z };

	// If the cone is centered around the north pole (0, 0, 1), return the vector as is
	if (cone_direction.x == 0 && cone_direction.y == 0) {
		return random_vector;
	}

	// Calculate the rotation axis (cross product of (0, 0, 1) and cone_direction)
	tfx_vec3_t north_pole = { 0, 0, 1 };
	tfx_vec3_t rotation_axis(-cone_direction.y, cone_direction.x, 0.f);
	rotation_axis = tfx__normalize_vec3_fast(&rotation_axis);

	// Calculate the rotation angle (acos of dot product of (0, 0, 1) and cone_direction)
	float rotation_angle = acosf(cone_direction.z);
	float cos = cosf(rotation_angle);
	float sin = sinf(rotation_angle);
	//return (a->x * b->x + a->y * b->y + a->z * b->z);
	float dot = (rotation_axis.x * x + rotation_axis.y * y);

	// Rotate the random vector to align with the cone direction
	// Use Rodrigues' rotation formula
	tfx_vec3_t rotated_vector = random_vector * cos + tfx__cross_product_vec3(&rotation_axis, &random_vector) * sin + rotation_axis * dot * (1.f - cos);

	return rotated_vector;
}

void RandomVectorInConeWide(tfxWideInt seed, tfxWideFloat velocity_normal_x, tfxWideFloat velocity_normal_y, tfxWideFloat velocity_normal_z, tfxWideFloat cone_angle, tfxWideFloat *random_x, tfxWideFloat *random_y, tfxWideFloat *random_z) {
	// Convert cone angle to radians

	cone_angle = tfxWideMin(cone_angle, tfx180RadiansWide);
	// Calculate the minimum z value for the cone
	tfxWideFloat min_z = tfxWideCos52s(cone_angle);

	tfxWideFloat max_uint = tfxWideSetSingle((float)UINT32_MAX);
	// Randomly sample z in [min_z, 1]
	tfxWideFloat z = tfxWideAdd(tfxWideDiv(SeedGenWide(seed), max_uint), tfxWideSetSingle(0.5f));
	z = tfxWideSub(tfxWIDEONE, tfxWideMul(tfxWideSub(tfxWIDEONE, min_z), z));

	// Randomly sample ϕ in [0, 2π)
	tfxWideFloat phi = tfxWideMul(tfxWideAdd(tfxWideDiv(SeedGenWide(seed), max_uint), tfxWideSetSingle(0.5f)), tfxWideSetSingle(2.f * tfxPI));

	// Calculate the corresponding x and y for the random point on the unit sphere
	tfxWideFloat sqrt_one_minus_z_squared = tfxWideSub(tfxWIDEONE, tfxWideMul(z, z));
	sqrt_one_minus_z_squared = tfxWideMul(tfxWideRSqrt(sqrt_one_minus_z_squared), sqrt_one_minus_z_squared);
	tfxWideFloat sin;
	tfxWideFloat cos;
	tfxWideSinCos(phi, &sin, &cos);
	tfxWideFloat x = tfxWideMul(sqrt_one_minus_z_squared, cos);
	tfxWideFloat y = tfxWideMul(sqrt_one_minus_z_squared, sin);

	// Calculate the rotation axis (cross product of (0, 0, 1) and cone_direction)
	tfx_vec3_t north_pole = { 0, 0, 1.f };
	tfxWideFloat rotation_axis_x, rotation_axis_y, rotation_axis_z;
	rotation_axis_x = tfxWideSub(tfxWideSetZero, velocity_normal_y);
	rotation_axis_y = velocity_normal_x;
	rotation_axis_z = tfxWideSetZero;
	tfxWideFloat length = tfxWideMul(rotation_axis_x, rotation_axis_x);
	length = tfxWideAdd(length, tfxWideMul(rotation_axis_y, rotation_axis_y));
	length = tfxWideMul(tfxWideRSqrt(length), length);
	rotation_axis_x = tfxWideDiv(rotation_axis_x, length);
	rotation_axis_y = tfxWideDiv(rotation_axis_y, length);

	tfxWideArray rotation_angle;
	tfxWideArray dir_z;
	dir_z.m = velocity_normal_z;
	rotation_angle.a[0] = acosf(dir_z.a[0]);
	rotation_angle.a[1] = acosf(dir_z.a[1]);
	rotation_angle.a[2] = acosf(dir_z.a[2]);
	rotation_angle.a[3] = acosf(dir_z.a[3]);
#if defined(tfxUSEAVX)
	rotation_angle.a[4] = acosf(dir_z.a[4]);
	rotation_angle.a[5] = acosf(dir_z.a[5]);
	rotation_angle.a[6] = acosf(dir_z.a[6]);
	rotation_angle.a[7] = acosf(dir_z.a[7]);
#endif

	// Rotate the random vector to align with the cone direction
	// Use Rodrigues' rotation formula
	tfxWideSinCos(rotation_angle.m, &sin, &cos);
	tfxWideFloat dot = tfxWideAdd(tfxWideMul(rotation_axis_x, x), tfxWideMul(rotation_axis_y, y));
	tfxWideFloat cx, cy, cz;
	tfx__wide_cross_product(rotation_axis_x, rotation_axis_y, rotation_axis_z, &x, &y, &z, &cx, &cy, &cz);
	//tfx_vec3_t rotated_vector = random_vector * cos + tfx__cross_product_vec3(&rotation_axis, &random_vector) * sin + rotation_axis * tfx__dot_product_vec3(&rotation_axis, &random_vector) * (1.f - cos);
	*random_x = tfxWideAdd(tfxWideAdd(tfxWideMul(x, cos), tfxWideMul(cx, sin)), tfxWideMul(tfxWideMul(rotation_axis_x, dot), tfxWideSub(tfxWIDEONE, cos)));
	*random_y = tfxWideAdd(tfxWideAdd(tfxWideMul(y, cos), tfxWideMul(cy, sin)), tfxWideMul(tfxWideMul(rotation_axis_y, dot), tfxWideSub(tfxWIDEONE, cos)));
	*random_z = tfxWideAdd(tfxWideAdd(tfxWideMul(z, cos), tfxWideMul(cz, sin)), tfxWideMul(tfxWideMul(rotation_axis_z, dot), tfxWideSub(tfxWIDEONE, cos)));
}

tfx_vec3_t GetEmissionDirection3d(tfx_particle_manager_t *pm, tfx_library_t *library, tfx_random_t *random, tfx_emitter_state_t &emitter, float emission_pitch, float emission_yaw, tfx_vec3_t local_position, tfx_vec3_t world_position) {
	float emission_angle_variation = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].properties.emission_range, emitter.frame);
	//----Emission
	float range = emission_angle_variation * .5f;

	tfx_vec3_t result;
	tfx_vec3_t tmp_position;
	if (emitter.handle.x + local_position.x == 0 && emitter.handle.y + local_position.y == 0 && emitter.handle.z + local_position.z == 0)
		tmp_position = emitter.emitter_size;
	else
		tmp_position = local_position + emitter.handle;

	tfx_emission_type emission_type = library->emitter_properties[emitter.properties_index].emission_type;
	tfx_emission_direction emission_direction = library->emitter_properties[emitter.properties_index].emission_direction;

	tfx_vec3_t to_handle(0.f, 1.f, 0.f);
	float parent_pitch = 0.f;
	float parent_yaw = 0.f;
	if (emission_type != tfxPoint) {
		if (emission_direction == tfxOutwards) {

			if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = tmp_position;
			else
				to_handle = world_position - emitter.world_position;

			to_handle = tfx__normalize_vec3(&to_handle);

		}
		else if (emission_direction == tfxInwards) {

			if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = -tmp_position;
			else
				to_handle = emitter.world_position - world_position;

			to_handle = tfx__normalize_vec3(&to_handle);

		}
		else if (emission_direction == tfxBothways) {

			if (emitter.emission_alternator) {

				if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = tmp_position;
				else
					to_handle = world_position - emitter.world_position;
			}
			else {

				if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = -tmp_position;
				else
					to_handle = emitter.world_position - world_position;
			}

			emitter.emission_alternator = !emitter.emission_alternator;
			to_handle = tfx__normalize_vec3(&to_handle);
		}
		else if (emission_direction == tfxSurface && emitter.property_flags & tfxEmitterPropertyFlags_relative_position) {
			if (emission_type == tfxEllipse || emission_type == tfxIcosphere) {
				to_handle = tfx__ellipse_surface_normal(local_position.x, local_position.y, local_position.z, emitter.emitter_size.x * .5f, emitter.emitter_size.y * .5f, emitter.emitter_size.z * .5f);
			}
			else if (emission_type == tfxCylinder) {
				float radius_x = emitter.emitter_size.x * .5f;
				float radius_z = emitter.emitter_size.z * .5f;
				to_handle = tfx__cylinder_surface_normal(local_position.x - emitter.handle.x, local_position.z - emitter.handle.z, radius_x, radius_z);
			}
		}
		else if (emission_direction == tfxSurface) {
			if (emission_type == tfxEllipse || emission_type == tfxIcosphere) {
				to_handle = tfx__ellipse_surface_normal(local_position.x - emitter.world_position.x - emitter.handle.x, local_position.y - emitter.world_position.y - emitter.handle.y, local_position.z - emitter.world_position.z - emitter.handle.z, emitter.emitter_size.x * .5f, emitter.emitter_size.y * .5f, emitter.emitter_size.z * .5f);
			}
			else if (emission_type == tfxCylinder) {
				float radius_x = emitter.emitter_size.x * .5f;
				float radius_z = emitter.emitter_size.z * .5f;
				to_handle = tfx__cylinder_surface_normal(local_position.x - emitter.world_position.x - emitter.handle.x, local_position.z - emitter.world_position.z - emitter.handle.z, radius_x, radius_z);
			}
		}
		else {
			parent_pitch = emitter.world_rotations.pitch;
			parent_yaw = emitter.world_rotations.yaw;
		}
	}
	else {
		parent_pitch = emitter.world_rotations.pitch;
		parent_yaw = emitter.world_rotations.yaw;
	}

	float pitch = asinf(-to_handle.y);
	float yaw = atan2(to_handle.x, to_handle.z);
	tfx_vec3_t direction;
	direction.z = cosf(emission_yaw + yaw) * cosf(emission_pitch + pitch);
	direction.y = -sinf(emission_pitch + pitch);
	direction.x = sinf(emission_yaw + yaw) * cosf(emission_pitch + pitch);
	tfx_vec3_t v = direction;
	if (range != 0) {
		v = RandomVectorInCone(random, v, range);
	}

	return tfx__normalize_vec3_fast(&v);
}

tfx_quaternion_t GetPathRotation2d(tfx_random_t *random, float range, float angle) {
	range *= .5f;
	float v = RandomRange(random, -range, range) + angle;
	tfx_quaternion_t quaternion;
	ToQuaternion2d(&quaternion, v);
	return quaternion;
}

tfx_quaternion_t GetPathRotation(tfx_random_t *random, float range, float pitch, float yaw, bool y_axis_only) {
	tfx_vec3_t direction;
	if (y_axis_only) {
		range *= 0.5f;
		yaw += RandomRange(random, -range, range);
		tfx_quaternion_t yaw_quaternion = QuaternionFromAxisAngle(0.0f, 1.0f, 0.0f, yaw);
		tfx_quaternion_t pitch_quaternion = QuaternionFromAxisAngle(1.0f, 0.0f, 0.0f, pitch);
		tfx_quaternion_t combined_quaternion = yaw_quaternion * pitch_quaternion;
		return combined_quaternion;
	}
	else {
		pitch -= tfx90Radians;
		direction.z = cosf(yaw) * cosf(pitch);
		direction.y = -sinf(pitch);
		direction.x = sinf(yaw) * cosf(pitch);
		tfx_vec3_t v = direction;
		if (range != 0) {
			v = RandomVectorInCone(random, v, range * .5f);
		}
		return QuaternionFromDirection(&v);
	}
}

void ResetEffectGraphs(tfx_effect_emitter_t *effect, bool add_node, bool compile) {
	tfx_library_t *library = effect->library;
	tfxU32 global = effect->global;
	ResetGraph(&library->global_graphs[global].life, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].life.type = tfxGlobal_life;
	ResetGraph(&library->global_graphs[global].amount, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].amount.type = tfxGlobal_amount;
	ResetGraph(&library->global_graphs[global].velocity, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].velocity.type = tfxGlobal_velocity;
	ResetGraph(&library->global_graphs[global].noise, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].noise.type = tfxGlobal_noise;
	ResetGraph(&library->global_graphs[global].width, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].width.type = tfxGlobal_width;
	ResetGraph(&library->global_graphs[global].height, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].height.type = tfxGlobal_height;
	ResetGraph(&library->global_graphs[global].weight, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].weight.type = tfxGlobal_weight;
	ResetGraph(&library->global_graphs[global].spin, 1.f, tfxGlobalPercentPresetSigned, add_node); library->global_graphs[global].spin.type = tfxGlobal_roll_spin;
	ResetGraph(&library->global_graphs[global].pitch_spin, 1.f, tfxGlobalPercentPresetSigned, add_node); library->global_graphs[global].pitch_spin.type = tfxGlobal_pitch_spin;
	ResetGraph(&library->global_graphs[global].yaw_spin, 1.f, tfxGlobalPercentPresetSigned, add_node); library->global_graphs[global].yaw_spin.type = tfxGlobal_yaw_spin;
	ResetGraph(&library->global_graphs[global].stretch, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].stretch.type = tfxGlobal_stretch;
	ResetGraph(&library->global_graphs[global].overal_scale, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].overal_scale.type = tfxGlobal_overal_scale;
	ResetGraph(&library->global_graphs[global].intensity, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].intensity.type = tfxGlobal_intensity;
	ResetGraph(&library->global_graphs[global].splatter, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].splatter.type = tfxGlobal_splatter;
	ResetGraph(&library->global_graphs[global].emitter_width, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].emitter_width.type = tfxGlobal_emitter_width;
	ResetGraph(&library->global_graphs[global].emitter_height, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].emitter_height.type = tfxGlobal_emitter_height;
	ResetGraph(&library->global_graphs[global].emitter_depth, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].emitter_depth.type = tfxGlobal_emitter_depth;
	if (compile) {
		CompileLibraryGlobalGraph(library, global);
	}
}

void ResetTransformGraphs(tfx_effect_emitter_t *effect, bool add_node, bool compile) {
	tfx_library_t *library = effect->library;
	tfxU32 transform_attributes = effect->transform_attributes;
	ResetGraph(&library->transform_attributes[transform_attributes].roll, 0.f, tfxAnglePreset, add_node); library->transform_attributes[transform_attributes].roll.type = tfxTransform_roll;
	ResetGraph(&library->transform_attributes[transform_attributes].pitch, 0.f, tfxAnglePreset, add_node); library->transform_attributes[transform_attributes].pitch.type = tfxTransform_pitch;
	ResetGraph(&library->transform_attributes[transform_attributes].yaw, 0.f, tfxAnglePreset, add_node); library->transform_attributes[transform_attributes].yaw.type = tfxTransform_yaw;
	ResetGraph(&library->transform_attributes[transform_attributes].translation_x, 0.f, tfxTranslationPreset, add_node); library->transform_attributes[transform_attributes].translation_x.type = tfxTransform_translate_x;
	ResetGraph(&library->transform_attributes[transform_attributes].translation_y, 0.f, tfxTranslationPreset, add_node); library->transform_attributes[transform_attributes].translation_y.type = tfxTransform_translate_y;
	ResetGraph(&library->transform_attributes[transform_attributes].translation_z, 0.f, tfxTranslationPreset, add_node); library->transform_attributes[transform_attributes].translation_z.type = tfxTransform_translate_z;
	if (compile) {
		CompileLibraryKeyframeGraph(library, transform_attributes);
	}
}

void ResetEmitterBaseGraphs(tfx_effect_emitter_t *effect, bool add_node, bool compile) {
	tfx_library_t *library = effect->library;
	tfxU32 emitter_attributes = effect->emitter_attributes;
	ResetGraph(&library->emitter_attributes[emitter_attributes].base.life, 1000.f, tfxLifePreset, add_node); library->emitter_attributes[emitter_attributes].base.life.type = tfxBase_life;
	ResetGraph(&library->emitter_attributes[emitter_attributes].base.amount, 1.f, tfxAmountPreset, add_node); library->emitter_attributes[emitter_attributes].base.amount.type = tfxBase_amount;
	ResetGraph(&library->emitter_attributes[emitter_attributes].base.velocity, 0.f, tfxVelocityPreset, add_node); library->emitter_attributes[emitter_attributes].base.velocity.type = tfxBase_velocity;
	ResetGraph(&library->emitter_attributes[emitter_attributes].base.width, 128.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].base.width.type = tfxBase_width;
	ResetGraph(&library->emitter_attributes[emitter_attributes].base.height, 128.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].base.height.type = tfxBase_height;
	ResetGraph(&library->emitter_attributes[emitter_attributes].base.weight, 0.f, tfxWeightPreset, add_node); library->emitter_attributes[emitter_attributes].base.weight.type = tfxBase_weight;
	ResetGraph(&library->emitter_attributes[emitter_attributes].base.spin, 0.f, tfxSpinPreset, add_node); library->emitter_attributes[emitter_attributes].base.spin.type = tfxBase_roll_spin;
	ResetGraph(&library->emitter_attributes[emitter_attributes].base.pitch_spin, 0.f, tfxSpinPreset, add_node); library->emitter_attributes[emitter_attributes].base.pitch_spin.type = tfxBase_pitch_spin;
	ResetGraph(&library->emitter_attributes[emitter_attributes].base.yaw_spin, 0.f, tfxSpinPreset, add_node); library->emitter_attributes[emitter_attributes].base.yaw_spin.type = tfxBase_yaw_spin;
	ResetGraph(&library->emitter_attributes[emitter_attributes].base.noise_offset, 0.f, tfxGlobalPercentPreset, add_node); library->emitter_attributes[emitter_attributes].base.noise_offset.type = tfxBase_noise_offset;
	if (compile) {
		CompileLibraryBaseGraph(library, emitter_attributes);
	}
}

void ResetEmitterPropertyGraphs(tfx_effect_emitter_t *effect, bool add_node, bool compile) {
	tfx_library_t *library = effect->library;
	tfxU32 emitter_attributes = effect->emitter_attributes;
	ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emission_pitch, 0.f, tfxAnglePreset, add_node); library->emitter_attributes[emitter_attributes].properties.emission_pitch.type = tfxProperty_emission_pitch;
	ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emission_yaw, 0.f, tfxAnglePreset, add_node); library->emitter_attributes[emitter_attributes].properties.emission_yaw.type = tfxProperty_emission_yaw;
	ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emission_range, 0.f, tfxEmissionRangePreset, add_node); library->emitter_attributes[emitter_attributes].properties.emission_range.type = tfxProperty_emission_range;
	ResetGraph(&library->emitter_attributes[emitter_attributes].properties.splatter, 0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].properties.splatter.type = tfxProperty_splatter;
	ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_width, 0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].properties.emitter_width.type = tfxProperty_emitter_width;
	ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_height, 0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].properties.emitter_height.type = tfxProperty_emitter_height;
	ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_depth, 0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].properties.emitter_depth.type = tfxProperty_emitter_depth;
	ResetGraph(&library->emitter_attributes[emitter_attributes].properties.extrusion, 0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].properties.extrusion.type = tfxProperty_extrusion;
	ResetGraph(&library->emitter_attributes[emitter_attributes].properties.arc_size, tfx__degrees_to_radians(360.f), tfxArcPreset, add_node); library->emitter_attributes[emitter_attributes].properties.arc_size.type = tfxProperty_arc_size;
	ResetGraph(&library->emitter_attributes[emitter_attributes].properties.arc_offset, 0.f, tfxArcPreset, add_node); library->emitter_attributes[emitter_attributes].properties.arc_offset.type = tfxProperty_arc_offset;
	if (compile) {
		CompileLibraryPropertyGraph(library, emitter_attributes);
	}
}

void ResetEmitterVariationGraphs(tfx_effect_emitter_t *effect, bool add_node, bool compile) {
	tfx_library_t *library = effect->library;
	tfxU32 emitter_attributes = effect->emitter_attributes;
	ResetGraph(&library->emitter_attributes[emitter_attributes].variation.life, 0.f, tfxLifePreset, add_node); library->emitter_attributes[emitter_attributes].variation.life.type = tfxVariation_life;
	ResetGraph(&library->emitter_attributes[emitter_attributes].variation.amount, 0.f, tfxAmountPreset, add_node); library->emitter_attributes[emitter_attributes].variation.amount.type = tfxVariation_amount;
	ResetGraph(&library->emitter_attributes[emitter_attributes].variation.velocity, 0.f, tfxVelocityPreset, add_node); library->emitter_attributes[emitter_attributes].variation.velocity.type = tfxVariation_velocity;
	ResetGraph(&library->emitter_attributes[emitter_attributes].variation.width, 0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].variation.width.type = tfxVariation_width;
	ResetGraph(&library->emitter_attributes[emitter_attributes].variation.height, 0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].variation.height.type = tfxVariation_height;
	ResetGraph(&library->emitter_attributes[emitter_attributes].variation.weight, 0.f, tfxWeightVariationPreset, add_node); library->emitter_attributes[emitter_attributes].variation.weight.type = tfxVariation_weight;
	ResetGraph(&library->emitter_attributes[emitter_attributes].variation.spin, 0.f, tfxSpinVariationPreset, add_node); library->emitter_attributes[emitter_attributes].variation.spin.type = tfxVariation_roll_spin;
	ResetGraph(&library->emitter_attributes[emitter_attributes].variation.pitch_spin, 0.f, tfxSpinVariationPreset, add_node); library->emitter_attributes[emitter_attributes].variation.pitch_spin.type = tfxVariation_pitch_spin;
	ResetGraph(&library->emitter_attributes[emitter_attributes].variation.yaw_spin, 0.f, tfxSpinVariationPreset, add_node); library->emitter_attributes[emitter_attributes].variation.yaw_spin.type = tfxVariation_yaw_spin;
	ResetGraph(&library->emitter_attributes[emitter_attributes].variation.noise_offset, 0.f, tfxNoiseOffsetVariationPreset, add_node); library->emitter_attributes[emitter_attributes].variation.noise_offset.type = tfxVariation_noise_offset;
	ResetGraph(&library->emitter_attributes[emitter_attributes].variation.noise_resolution, 300.f, tfxNoiseResolutionPreset, add_node); library->emitter_attributes[emitter_attributes].variation.noise_resolution.type = tfxVariation_noise_resolution;
	ResetGraph(&library->emitter_attributes[emitter_attributes].variation.motion_randomness, 1.f, tfxNoiseResolutionPreset, add_node); library->emitter_attributes[emitter_attributes].variation.motion_randomness.type = tfxVariation_motion_randomness;
	if (compile) {
		CompileLibraryVariationGraph(library, emitter_attributes);
	}
}

void ResetEmitterOvertimeGraphs(tfx_effect_emitter_t *effect, bool add_node, bool compile) {
	tfx_library_t *library = effect->library;
	tfxU32 emitter_attributes = effect->emitter_attributes;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity, 1.f, tfxVelocityOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.velocity.type = tfxOvertime_velocity;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster, 1.f, tfxGlobalPercentPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.type = tfxOvertime_velocity_adjuster;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.width, 1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.width.type = tfxOvertime_width;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.height, 1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.height.type = tfxOvertime_height;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.weight, 1.f, tfxWeightOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.weight.type = tfxOvertime_weight;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.spin, 0.f, tfxSpinOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.spin.type = tfxOvertime_roll_spin;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.pitch_spin, 0.f, tfxSpinOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.pitch_spin.type = tfxOvertime_pitch_spin;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.yaw_spin, 0.f, tfxSpinOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.yaw_spin.type = tfxOvertime_yaw_spin;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.stretch, 0.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.stretch.type = tfxOvertime_stretch;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.red, 1.f, tfxColorPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.red.type = tfxOvertime_red;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.green, 1.f, tfxColorPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.green.type = tfxOvertime_green;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.blue, 1.f, tfxColorPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.blue.type = tfxOvertime_blue;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.blendfactor, 1.f, tfxOpacityOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.blendfactor.type = tfxOvertime_blendfactor;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.intensity, 1.f, tfxIntensityOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.intensity.type = tfxOvertime_intensity;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.red_hint, 1.f, tfxColorPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.red_hint.type = tfxOvertime_red_hint;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.green_hint, 1.f, tfxColorPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.green_hint.type = tfxOvertime_green_hint;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.blue_hint, 1.f, tfxColorPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.blue_hint.type = tfxOvertime_blue_hint;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.blendfactor_hint, 1.f, tfxOpacityOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.blendfactor_hint.type = tfxOvertime_blendfactor_hint;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.alpha_sharpness, 1.f, tfxOpacityOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.alpha_sharpness.type = tfxOvertime_alpha_sharpness;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.curved_alpha, 1.f, tfxOpacityOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.curved_alpha.type = tfxOvertime_curved_alpha;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance, 30.f, tfxVelocityTurbulancePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance.type = tfxOvertime_velocity_turbulance;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.stretch, 0.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.stretch.type = tfxOvertime_stretch;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.direction_turbulance, 0.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.direction_turbulance.type = tfxOvertime_direction_turbulance;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.direction, 0.f, tfxDirectionOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.direction.type = tfxOvertime_direction;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.noise_resolution, 1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.noise_resolution.type = tfxOvertime_noise_resolution;
	ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.motion_randomness, 0.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.motion_randomness.type = tfxOvertime_motion_randomness;
	if (compile) {
		CompileLibraryOvertimeGraph(library, emitter_attributes);
	}
}

void ResetEmitterFactorGraphs(tfx_effect_emitter_t *effect, bool add_node, bool compile) {
	tfx_library_t *library = effect->library;
	tfxU32 emitter_attributes = effect->emitter_attributes;
	ResetGraph(&library->emitter_attributes[emitter_attributes].factor.life, 1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].factor.life.type = tfxFactor_life;
	ResetGraph(&library->emitter_attributes[emitter_attributes].factor.velocity, 1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].factor.velocity.type = tfxFactor_velocity;
	ResetGraph(&library->emitter_attributes[emitter_attributes].factor.size, 1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].factor.size.type = tfxFactor_size;
	ResetGraph(&library->emitter_attributes[emitter_attributes].factor.intensity, 1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].factor.intensity.type = tfxFactor_intensity;
	if (compile) {
		CompileLibraryFactorGraph(library, emitter_attributes);
	}
}

void ResetEmitterGraphs(tfx_effect_emitter_t *effect, bool add_node, bool compile) {
	ResetEmitterBaseGraphs(effect, add_node, compile);
	ResetEmitterPropertyGraphs(effect, add_node, compile);
	ResetEmitterVariationGraphs(effect, add_node, compile);
	UpdateEffectMaxLife(effect);
	ResetEmitterOvertimeGraphs(effect, add_node, compile);
	ResetEmitterFactorGraphs(effect, add_node, compile);
}

void InitialiseUninitialisedGraphs(tfx_effect_emitter_t *effect) {
	tfx_library_t *library = effect->library;
	tfxU32 transform_attributes = effect->transform_attributes;
	if (library->transform_attributes[transform_attributes].translation_x.nodes.size() == 0) ResetGraph(&library->transform_attributes[transform_attributes].translation_x, 0.f, tfxTranslationPreset);
	if (library->transform_attributes[transform_attributes].translation_y.nodes.size() == 0) ResetGraph(&library->transform_attributes[transform_attributes].translation_y, 0.f, tfxTranslationPreset);
	if (library->transform_attributes[transform_attributes].translation_z.nodes.size() == 0) ResetGraph(&library->transform_attributes[transform_attributes].translation_z, 0.f, tfxTranslationPreset);
	if (library->transform_attributes[transform_attributes].roll.nodes.size() == 0) ResetGraph(&library->transform_attributes[transform_attributes].roll, 0.f, tfxAnglePreset);
	if (library->transform_attributes[transform_attributes].pitch.nodes.size() == 0) ResetGraph(&library->transform_attributes[transform_attributes].pitch, 0.f, tfxAnglePreset);
	if (library->transform_attributes[transform_attributes].yaw.nodes.size() == 0) ResetGraph(&library->transform_attributes[transform_attributes].yaw, 0.f, tfxAnglePreset);

	if (effect->type == tfxEffectType && effect->global != tfxINVALID) {
		tfxU32 global = effect->global;
		if (library->global_graphs[global].life.nodes.size() == 0) ResetGraph(&library->global_graphs[global].life, 1.f, tfxGlobalPercentPreset);
		if (library->global_graphs[global].amount.nodes.size() == 0) ResetGraph(&library->global_graphs[global].amount, 1.f, tfxGlobalPercentPreset);
		if (library->global_graphs[global].velocity.nodes.size() == 0) ResetGraph(&library->global_graphs[global].velocity, 1.f, tfxGlobalPercentPreset);
		if (library->global_graphs[global].noise.nodes.size() == 0) ResetGraph(&library->global_graphs[global].noise, 1.f, tfxGlobalPercentPreset);
		if (library->global_graphs[global].width.nodes.size() == 0) ResetGraph(&library->global_graphs[global].width, 1.f, tfxGlobalPercentPreset);
		if (library->global_graphs[global].height.nodes.size() == 0) ResetGraph(&library->global_graphs[global].height, 1.f, tfxGlobalPercentPreset);
		if (library->global_graphs[global].weight.nodes.size() == 0) ResetGraph(&library->global_graphs[global].weight, 1.f, tfxGlobalPercentPreset);
		if (library->global_graphs[global].spin.nodes.size() == 0) ResetGraph(&library->global_graphs[global].spin, 1.f, tfxGlobalPercentPresetSigned);
		if (library->global_graphs[global].pitch_spin.nodes.size() == 0) ResetGraph(&library->global_graphs[global].pitch_spin, 1.f, tfxGlobalPercentPresetSigned);
		if (library->global_graphs[global].yaw_spin.nodes.size() == 0) ResetGraph(&library->global_graphs[global].yaw_spin, 1.f, tfxGlobalPercentPresetSigned);
		if (library->global_graphs[global].stretch.nodes.size() == 0) ResetGraph(&library->global_graphs[global].stretch, 1.f, tfxGlobalPercentPreset);
		if (library->global_graphs[global].overal_scale.nodes.size() == 0) ResetGraph(&library->global_graphs[global].overal_scale, 1.f, tfxGlobalPercentPreset);
		if (library->global_graphs[global].intensity.nodes.size() == 0) ResetGraph(&library->global_graphs[global].intensity, 1.f, tfxGlobalPercentPreset);
		if (library->global_graphs[global].splatter.nodes.size() == 0) ResetGraph(&library->global_graphs[global].splatter, 1.f, tfxGlobalPercentPreset);
		if (library->global_graphs[global].emitter_width.nodes.size() == 0) ResetGraph(&library->global_graphs[global].emitter_width, 1.f, tfxGlobalPercentPreset);
		if (library->global_graphs[global].emitter_height.nodes.size() == 0) ResetGraph(&library->global_graphs[global].emitter_height, 1.f, tfxGlobalPercentPreset);
		if (library->global_graphs[global].emitter_depth.nodes.size() == 0) ResetGraph(&library->global_graphs[global].emitter_depth, 1.f, tfxGlobalPercentPreset);
	}

	if (effect->type == tfxEmitterType && effect->emitter_attributes != tfxINVALID) {
		tfxU32 emitter_attributes = effect->emitter_attributes;
		if (library->emitter_attributes[emitter_attributes].base.life.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].base.life, 1000.f, tfxLifePreset);
		if (library->emitter_attributes[emitter_attributes].base.amount.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].base.amount, 1.f, tfxAmountPreset);
		if (library->emitter_attributes[emitter_attributes].base.velocity.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].base.velocity, 0.f, tfxVelocityPreset);
		if (library->emitter_attributes[emitter_attributes].base.width.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].base.width, 128.f, tfxDimensionsPreset);
		if (library->emitter_attributes[emitter_attributes].base.height.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].base.height, 128.f, tfxDimensionsPreset);
		if (library->emitter_attributes[emitter_attributes].base.weight.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].base.weight, 0.f, tfxWeightPreset);
		if (library->emitter_attributes[emitter_attributes].base.spin.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].base.spin, 0.f, tfxSpinPreset);
		if (library->emitter_attributes[emitter_attributes].base.pitch_spin.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].base.pitch_spin, 0.f, tfxSpinPreset);
		if (library->emitter_attributes[emitter_attributes].base.yaw_spin.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].base.yaw_spin, 0.f, tfxSpinPreset);
		if (library->emitter_attributes[emitter_attributes].base.noise_offset.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].base.noise_offset, 0.f, tfxGlobalPercentPreset);

		if (library->emitter_attributes[emitter_attributes].properties.emission_pitch.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emission_pitch, 0.f, tfxAnglePreset);
		if (library->emitter_attributes[emitter_attributes].properties.emission_yaw.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emission_yaw, 0.f, tfxAnglePreset);
		if (library->emitter_attributes[emitter_attributes].properties.emission_range.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emission_range, 0.f, tfxEmissionRangePreset);
		if (library->emitter_attributes[emitter_attributes].properties.splatter.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.splatter, 0.f, tfxDimensionsPreset);
		if (library->emitter_attributes[emitter_attributes].properties.emitter_width.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_width, 0.f, tfxDimensionsPreset);
		if (library->emitter_attributes[emitter_attributes].properties.emitter_height.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_height, 0.f, tfxDimensionsPreset);
		if (library->emitter_attributes[emitter_attributes].properties.emitter_depth.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_depth, 0.f, tfxDimensionsPreset);
		if (library->emitter_attributes[emitter_attributes].properties.extrusion.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.extrusion, 0.f, tfxDimensionsPreset);
		if (library->emitter_attributes[emitter_attributes].properties.arc_size.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.arc_size, tfx__degrees_to_radians(360.f), tfxArcPreset);
		if (library->emitter_attributes[emitter_attributes].properties.arc_offset.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.arc_offset, 0.f, tfxArcPreset);

		if (library->emitter_attributes[emitter_attributes].variation.life.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.life, 0.f, tfxLifePreset);
		if (library->emitter_attributes[emitter_attributes].variation.amount.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.amount, 0.f, tfxAmountPreset);
		if (library->emitter_attributes[emitter_attributes].variation.velocity.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.velocity, 0.f, tfxVelocityPreset);
		if (library->emitter_attributes[emitter_attributes].variation.width.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.width, 0.f, tfxDimensionsPreset);
		if (library->emitter_attributes[emitter_attributes].variation.height.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.height, 0.f, tfxDimensionsPreset);
		if (library->emitter_attributes[emitter_attributes].variation.weight.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.weight, 0.f, tfxWeightVariationPreset);
		if (library->emitter_attributes[emitter_attributes].variation.spin.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.spin, 0.f, tfxSpinVariationPreset);
		if (library->emitter_attributes[emitter_attributes].variation.pitch_spin.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.pitch_spin, 0.f, tfxSpinVariationPreset);
		if (library->emitter_attributes[emitter_attributes].variation.yaw_spin.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.yaw_spin, 0.f, tfxSpinVariationPreset);
		if (library->emitter_attributes[emitter_attributes].variation.noise_offset.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.noise_offset, 0.f, tfxNoiseOffsetVariationPreset);
		if (library->emitter_attributes[emitter_attributes].variation.noise_resolution.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.noise_resolution, 300.f, tfxNoiseResolutionPreset);
		if (library->emitter_attributes[emitter_attributes].variation.motion_randomness.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.motion_randomness, 1.f, tfxNoiseResolutionPreset);

		if (library->emitter_attributes[emitter_attributes].overtime.velocity.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity, 1.f, tfxVelocityOvertimePreset);
		if (library->emitter_attributes[emitter_attributes].overtime.width.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.width, 1.f, tfxPercentOvertime);
		if (library->emitter_attributes[emitter_attributes].overtime.height.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.height, 1.f, tfxPercentOvertime);
		if (library->emitter_attributes[emitter_attributes].overtime.weight.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.weight, 1.f, tfxWeightOvertimePreset);
		if (library->emitter_attributes[emitter_attributes].overtime.spin.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.spin, 1.f, tfxSpinOvertimePreset);
		if (library->emitter_attributes[emitter_attributes].overtime.pitch_spin.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.pitch_spin, 1.f, tfxSpinOvertimePreset);
		if (library->emitter_attributes[emitter_attributes].overtime.yaw_spin.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.yaw_spin, 1.f, tfxSpinOvertimePreset);
		if (library->emitter_attributes[emitter_attributes].overtime.stretch.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.stretch, 0.f, tfxPercentOvertime);
		if (library->emitter_attributes[emitter_attributes].overtime.red.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.red, 1.f, tfxColorPreset);
		if (library->emitter_attributes[emitter_attributes].overtime.green.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.green, 1.f, tfxColorPreset);
		if (library->emitter_attributes[emitter_attributes].overtime.blue.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.blue, 1.f, tfxColorPreset);
		if (library->emitter_attributes[emitter_attributes].overtime.blendfactor.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.blendfactor, 1.f, tfxOpacityOvertimePreset);
		if (library->emitter_attributes[emitter_attributes].overtime.intensity.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.intensity, 1.f, tfxIntensityOvertimePreset);
		if (library->emitter_attributes[emitter_attributes].overtime.red_hint.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.red_hint, 1.f, tfxColorPreset);
		if (library->emitter_attributes[emitter_attributes].overtime.green_hint.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.green_hint, 1.f, tfxColorPreset);
		if (library->emitter_attributes[emitter_attributes].overtime.blue_hint.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.blue_hint, 1.f, tfxColorPreset);
		if (library->emitter_attributes[emitter_attributes].overtime.blendfactor_hint.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.blendfactor_hint, 1.f, tfxOpacityOvertimePreset);
		if (library->emitter_attributes[emitter_attributes].overtime.alpha_sharpness.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.alpha_sharpness, 0.f, tfxOpacityOvertimePreset);
		if (library->emitter_attributes[emitter_attributes].overtime.curved_alpha.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.curved_alpha, 1.f, tfxOpacityOvertimePreset);
		if (library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance, 0.f, tfxFrameratePreset);
		if (library->emitter_attributes[emitter_attributes].overtime.direction_turbulance.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.direction_turbulance, 0.f, tfxPercentOvertime);
		if (library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster, 1.f, tfxGlobalPercentPreset);
		if (library->emitter_attributes[emitter_attributes].overtime.direction.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.direction, 0.f, tfxDirectionOvertimePreset);
		if (library->emitter_attributes[emitter_attributes].overtime.noise_resolution.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.noise_resolution, 1.f, tfxPercentOvertime);
		if (library->emitter_attributes[emitter_attributes].overtime.motion_randomness.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.motion_randomness, 0.f, tfxPercentOvertime);

		if (library->emitter_attributes[emitter_attributes].factor.life.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].factor.life, 1.f, tfxPercentOvertime);
		if (library->emitter_attributes[emitter_attributes].factor.velocity.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].factor.velocity, 1.f, tfxPercentOvertime);
		if (library->emitter_attributes[emitter_attributes].factor.size.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].factor.size, 1.f, tfxPercentOvertime);
		if (library->emitter_attributes[emitter_attributes].factor.intensity.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].factor.intensity, 1.f, tfxPercentOvertime);
	}
}

void SetEffectName(tfx_effect_emitter_t *effect, const char *n) {
	GetEffectInfo(effect)->name = n;
}

void AddEmitterColorOvertime(tfx_effect_emitter_t *effect, float frame, tfx_rgb_t color) {
	AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.red, frame, color.r);
	AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.green, frame, color.g);
	AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.blue, frame, color.b);
}

void AddEmitterColorHintOvertime(tfx_effect_emitter_t *effect, float frame, tfx_rgb_t color) {
	AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.red_hint, frame, color.r);
	AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.green_hint, frame, color.g);
	AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.blue_hint, frame, color.b);
}

void SetEffectUserData(tfx_effect_emitter_t *effect, void *data) {
	effect->user_data = data;
}

void *GetEffectUserData(tfx_effect_emitter_t *effect) {
	return effect->user_data;
}

tfx_emitter_properties_t *GetEffectProperties(tfx_effect_emitter_t *effect) {
	return &effect->library->emitter_properties[effect->property_index];
}

bool RenameSubEffector(tfx_effect_emitter_t *emitter, const char *new_name) {
	TFX_ASSERT(emitter->parent);    //Must be an emitter or sub effect with a parent
	if (!EffectNameExists(emitter->parent, emitter, new_name) && strlen(new_name) > 0) {
		SetEffectName(emitter, new_name);
		UpdateLibraryEffectPaths(emitter->library);
		return true;
	}

	return false;
}

bool EffectNameExists(tfx_effect_emitter_t *in_effect, tfx_effect_emitter_t *excluding_effect, const char *name) {
	for (auto &e : GetEffectInfo(in_effect)->sub_effectors) {
		if (excluding_effect != &e) {
			if (GetEffectInfo(&e)->name == name) {
				return true;
			}
		}
	}

	return false;
}

void ReIndexEffect(tfx_effect_emitter_t *effect) {
	tfxU32 index = 0;
	for (auto &e : GetEffectInfo(effect)->sub_effectors) {
		e.library_index = index++;
		e.parent = effect;
		ReIndexEffect(&e);
	}
}

void CountEffectChildren(tfx_effect_emitter_t *effect, int *emitters, int *effects) {
	tmpStack(tfx_effect_emitter_t *, stack);
	stack.push_back(effect);
	*emitters = 0;
	*effects = 0;
	while (!stack.empty()) {
		tfx_effect_emitter_t *current = stack.pop_back();
		if (current->type == tfxEffectType)
			(*effects)++;
		else if (current->type == tfxEmitterType)
			(*emitters)++;
		for (auto &sub : GetEffectInfo(current)->sub_effectors) {
			stack.push_back(&sub);
		}
	}
}

tfx_effect_emitter_t *tfx_GetRootEffect(tfx_effect_emitter_t *effect) {
	if (!effect->parent || effect->parent->type == tfxFolder) {
		return effect;
	}
	tfx_effect_emitter_t *p = effect->parent;
	tfxU32 timeout = 0;
	while (p || ++timeout < 100) {
		if (!p->parent)
			return p;
		p = p->parent;
	}
	return nullptr;
}

bool IsRootEffect(tfx_effect_emitter_t *effect) {
	if (effect->type != tfxEffectType) return false;
	if (effect->type == tfxEffectType && !effect->parent) return true;
	if (effect->parent && effect->parent->type == tfxFolder) return true;
	return false;
}

void ResetEffectParents(tfx_effect_emitter_t *effect) {
	effect->parent = nullptr;
	for (auto &e : GetEffectInfo(effect)->sub_effectors) {
		ResetEffectParents(effect);
	}
}

tfx_effect_emitter_t *MoveEffectUp(tfx_effect_emitter_t *emitter) {
	tfx_effect_emitter_t *parent = emitter->parent;
	if (emitter->library_index > 0) {
		tfxU32 new_index = emitter->library_index - 1;
		std::swap(GetEffectInfo(parent)->sub_effectors[emitter->library_index], GetEffectInfo(parent)->sub_effectors[new_index]);
		ReIndexEffect(parent);
		UpdateLibraryEffectPaths(parent->library);
		return &GetEffectInfo(parent)->sub_effectors[new_index];
	}

	return nullptr;
}

tfx_effect_emitter_t *MoveEffectDown(tfx_effect_emitter_t *emitter) {
	tfx_effect_emitter_t *parent = emitter->parent;
	if (emitter->library_index < GetEffectInfo(parent)->sub_effectors.size() - 1) {
		tfxU32 new_index = emitter->library_index + 1;
		std::swap(GetEffectInfo(parent)->sub_effectors[emitter->library_index], GetEffectInfo(parent)->sub_effectors[new_index]);
		ReIndexEffect(parent);
		UpdateLibraryEffectPaths(parent->library);
		return &GetEffectInfo(parent)->sub_effectors[new_index];
	}
	return nullptr;
}

void DeleteEmitterFromEffect(tfx_effect_emitter_t *emitter) {
	tfx_effect_emitter_t *parent = emitter->parent;
	tfx_library_t *library = emitter->library;
	tmpStack(tfx_effect_emitter_t, stack);
	stack.push_back(*emitter);
	while (stack.size()) {
		tfx_effect_emitter_t &current = stack.pop_back();
		if (current.type == tfxEffectType && !current.parent) {
			FreeLibraryGlobal(library, current.global);
			FreeLibraryKeyframes(library, current.transform_attributes);
		} else if (current.type == tfxEmitterType) {
			FreeLibraryEmitterAttributes(library, current.emitter_attributes);
			FreeLibraryKeyframes(library, current.transform_attributes);
		}
		for (auto &sub : GetEffectInfo(&current)->sub_effectors) {
			stack.push_back(sub);
		}
	}
	GetEffectInfo(parent)->sub_effectors.erase(emitter);

	ReIndexEffect(parent);
	if (library) {
		UpdateLibraryEffectPaths(library);
	}
}

void CleanUpEffect(tfx_effect_emitter_t *effect) {
	if (GetEffectInfo(effect)->sub_effectors.size()) {
		tmpStack(tfx_effect_emitter_t, stack);
		stack.push_back(*effect);
		while (stack.size()) {
			tfx_effect_emitter_t current = stack.pop_back();
			if (current.type == tfxEffectType && !current.parent) {
				FreeLibraryGlobal(effect->library, current.global);
				FreeLibraryKeyframes(effect->library, current.transform_attributes);
			} else if (current.type == tfxEmitterType) {
				FreeLibraryEmitterAttributes(effect->library, current.emitter_attributes);
				FreeLibraryKeyframes(effect->library, current.transform_attributes);
			}
			for (auto &sub : GetEffectInfo(&current)->sub_effectors) {
				stack.push_back(sub);
			}
			GetEffectInfo(&current)->sub_effectors.free_all();
			GetEffectInfo(&current)->path.free_all();
			FreeLibraryProperties(effect->library, current.property_index);
			FreeLibraryInfo(effect->library, current.info_index);
		}
	}

	ReIndexEffect(effect);
}

void CloneEffect(tfx_effect_emitter_t *effect_to_clone, tfx_effect_emitter_t *clone, tfx_effect_emitter_t *root_parent, tfx_library_t *destination_library, tfxEffectCloningFlags flags) {
	//tfxU32 size = library->global_graphs[0].amount.lookup.values.capacity;
	*clone = *effect_to_clone;
	clone->info_index = CloneLibraryInfo(clone->library, effect_to_clone->info_index, destination_library);
	if (clone->type != tfxFolder) {
		clone->property_index = CloneLibraryProperties(clone->library, GetEffectProperties(effect_to_clone), destination_library);
	}
	clone->property_flags |= tfxEmitterPropertyFlags_enabled;
	if (!(flags & tfxEffectCloningFlags_keep_user_data))
		clone->user_data = nullptr;
	clone->library = destination_library;
	GetEffectInfo(clone)->sub_effectors.clear();

	tfx_library_t *library = effect_to_clone->library;

	if (effect_to_clone->type == tfxEffectType) {
		if (root_parent == clone) {
			clone->global = flags & tfxEffectCloningFlags_clone_graphs ? CloneLibraryGlobal(library, effect_to_clone->global, destination_library) : clone->global = effect_to_clone->global;
			clone->transform_attributes = flags & tfxEffectCloningFlags_clone_graphs ? CloneLibraryKeyframes(library, effect_to_clone->transform_attributes, destination_library) : clone->transform_attributes = effect_to_clone->transform_attributes;
			if (flags & tfxEffectCloningFlags_compile_graphs) {
				CompileLibraryGlobalGraph(clone->library, clone->global);
				CompileLibraryKeyframeGraph(clone->library, clone->transform_attributes);
			}
		} else {
			clone->transform_attributes = flags & tfxEffectCloningFlags_clone_graphs ? CloneLibraryKeyframes(library, effect_to_clone->transform_attributes, destination_library) : clone->transform_attributes = effect_to_clone->transform_attributes;
			if (flags & tfxEffectCloningFlags_compile_graphs) {
				CompileLibraryKeyframeGraph(clone->library, clone->transform_attributes);
			}
			if (!(flags & tfxEffectCloningFlags_force_clone_global)) {
				clone->global = root_parent->global;
			} else {
				clone->global = CloneLibraryGlobal(library, root_parent->global, destination_library);
				if (flags & tfxEffectCloningFlags_compile_graphs)
					CompileLibraryGlobalGraph(clone->library, clone->global);
			}
		}
	} else if (effect_to_clone->type == tfxEmitterType) {
		clone->emitter_attributes = flags & tfxEffectCloningFlags_clone_graphs ? CloneLibraryEmitterAttributes(library, effect_to_clone->emitter_attributes, destination_library) : effect_to_clone->emitter_attributes;
		clone->transform_attributes = flags & tfxEffectCloningFlags_clone_graphs ? CloneLibraryKeyframes(library, effect_to_clone->transform_attributes, destination_library) : effect_to_clone->transform_attributes;
		UpdateEffectMaxLife(clone);
		if (flags & tfxEffectCloningFlags_compile_graphs) {
			CompileLibraryKeyframeGraph(clone->library, clone->transform_attributes);
			CompileLibraryPropertyGraph(clone->library, clone->emitter_attributes);
			CompileLibraryBaseGraph(clone->library, clone->emitter_attributes);
			CompileLibraryVariationGraph(clone->library, clone->emitter_attributes);
			CompileLibraryOvertimeGraph(clone->library, clone->emitter_attributes, false);
			CompileLibraryFactorGraph(clone->library, clone->emitter_attributes);
			if (destination_library != effect_to_clone->library) {
				MaybeInsertColorRampBitmap(destination_library, &destination_library->emitter_attributes[clone->emitter_attributes].overtime, 0);
				MaybeInsertColorRampBitmap(destination_library, &destination_library->emitter_attributes[clone->emitter_attributes].overtime, 1);
			}
		}
		if (clone->path_attributes != tfxINVALID) {
			tfx_emitter_path_t path_copy = CopyPath(&library->paths[clone->path_attributes], "");
			clone->path_attributes = destination_library->paths.size();
			destination_library->paths.push_back(path_copy);
			if (destination_library->paths.back().flags & tfxPathFlags_2d) {
				BuildPathNodes2d(&destination_library->paths.back());
			}
			else {
				BuildPathNodes3d(&destination_library->paths.back());
			}
		}
	}

	for (auto &e : GetEffectInfo(effect_to_clone)->sub_effectors) {
		if (e.type == tfxEmitterType) {
			tfx_effect_emitter_t emitter_copy;
			CloneEffect(&e, &emitter_copy, root_parent, destination_library, flags);
			if (!(flags & tfxEffectCloningFlags_keep_user_data))
				emitter_copy.user_data = nullptr;
			AddEmitterToEffect(clone, &emitter_copy);
		}
		else if (e.type == tfxEffectType) {
			tfx_effect_emitter_t effect_copy;
			if (clone->type == tfxFolder)
				CloneEffect(&e, &effect_copy, &effect_copy, destination_library, flags);
			else
				CloneEffect(&e, &effect_copy, root_parent, destination_library, flags);
			if (!(flags & tfxEffectCloningFlags_keep_user_data)) {
				effect_copy.user_data = nullptr;
			}
			AddEffectToEmitter(clone, &effect_copy);
		}
	}
}

void AddTemplatePath(tfx_effect_template_t *effect_template, tfx_effect_emitter_t *effect_emitter, tfx_str256_t path) {
	effect_template->paths.Insert(path, effect_emitter);
	for (auto &sub : GetEffectInfo(effect_emitter)->sub_effectors) {
		tfx_str256_t sub_path = path;
		sub_path.Appendf("/%s", GetEffectInfo(&sub)->name.c_str());
		AddTemplatePath(effect_template, &sub, sub_path);
	}
}

bool PrepareEffectTemplate(tfx_library_t *library, const char *name, tfx_effect_template_t *effect_template) {
	ResetTemplate(effect_template);
	if (library->effect_paths.ValidName(name)) {
		PrepareLibraryEffectTemplate(library, name, effect_template);
		return true;
	}
	else {
		TFX_ASSERT(0);    //Not a valid effect name, make sure the effect exists in the library, name is case sensitive.
	}
	return false;
}

void EnableAllEmitters(tfx_effect_emitter_t *effect) {
	for (auto &e : GetEffectInfo(effect)->sub_effectors) {
		e.property_flags |= tfxEmitterPropertyFlags_enabled;
		EnableAllEmitters(&e);
	}
}

void EnableEmitter(tfx_effect_emitter_t *effect) {
	effect->property_flags |= tfxEmitterPropertyFlags_enabled;
}

void DisableAllEmitters(tfx_effect_emitter_t *effect) {
	for (auto &e : GetEffectInfo(effect)->sub_effectors) {
		e.property_flags &= ~tfxEmitterPropertyFlags_enabled;
		DisableAllEmitters(&e);
	}
}

void DisableAllEmittersExcept(tfx_effect_emitter_t *effect, tfx_effect_emitter_t *emitter) {
	for (auto &e : GetEffectInfo(effect)->sub_effectors) {
		if (e.library_index == emitter->library_index)
			e.property_flags |= tfxEmitterPropertyFlags_enabled;
		else
			e.property_flags &= ~tfxEmitterPropertyFlags_enabled;
	}
}

tfx_emitter_attributes_t *GetEmitterAttributes(tfx_effect_emitter_t *emitter) {
	TFX_ASSERT(emitter->type == tfxEmitterType);            //Must be an emitter
	TFX_ASSERT(emitter->emitter_attributes != tfxINVALID);    //Must be a valid emitter_attributes index into the library;
	return &emitter->library->emitter_attributes[emitter->emitter_attributes];
}

tfx_graph_t *GetEffectGraphByType(tfx_effect_emitter_t *effect, tfx_graph_type type) {
	tfx_library_t *library = effect->library;

	if (type < TFX_GLOBAL_COUNT) {
		return &((tfx_graph_t *)&library->global_graphs[effect->global])[type];
	}
	else if (type >= TFX_PROPERTY_START && type < TFX_BASE_START) {
		int ref = type - TFX_PROPERTY_START;
		return &((tfx_graph_t *)&library->emitter_attributes[effect->emitter_attributes].properties)[ref];
	}
	else if (type >= TFX_BASE_START && type < TFX_VARIATION_START) {
		int ref = type - TFX_BASE_START;
		return &((tfx_graph_t *)&library->emitter_attributes[effect->emitter_attributes].base)[ref];
	}
	else if (type >= TFX_VARIATION_START && type < TFX_OVERTIME_START) {
		int ref = type - TFX_VARIATION_START;
		return &((tfx_graph_t *)&library->emitter_attributes[effect->emitter_attributes].variation)[ref];
	}
	else if (type >= TFX_OVERTIME_START && type < TFX_FACTOR_START) {
		int ref = type - TFX_OVERTIME_START;
		return &((tfx_graph_t *)&library->emitter_attributes[effect->emitter_attributes].overtime)[ref];
	}
	else if (type >= TFX_FACTOR_START && type < TFX_TRANSFORM_START) {
		int ref = type - TFX_FACTOR_START;
		return &((tfx_graph_t *)&library->emitter_attributes[effect->emitter_attributes].factor)[ref];
	}
	else if (type >= TFX_TRANSFORM_START) {
		int ref = type - TFX_TRANSFORM_START;
		return &((tfx_graph_t *)&library->transform_attributes[effect->transform_attributes].roll)[ref];
	}

	return nullptr;

}

tfxU32 GetEffectGraphIndexByType(tfx_effect_emitter_t *effect, tfx_graph_type type) {

	if (type < TFX_GLOBAL_COUNT) {
		return effect->global;
	}
	else if (type < TFX_TRANSFORM_START) {
		return effect->emitter_attributes;
	}
	else {
		return effect->transform_attributes;
	}

}

void FreeEffectGraphs(tfx_effect_emitter_t *effect) {

	tfx_library_t *library = effect->library;

	if (effect->type == tfxEffectType) {
		tfxU32 global = effect->global;
		FreeGraph(&library->global_graphs[global].life);
		FreeGraph(&library->global_graphs[global].amount);
		FreeGraph(&library->global_graphs[global].velocity);
		FreeGraph(&library->global_graphs[global].width);
		FreeGraph(&library->global_graphs[global].height);
		FreeGraph(&library->global_graphs[global].weight);
		FreeGraph(&library->global_graphs[global].spin);
		FreeGraph(&library->global_graphs[global].stretch);
		FreeGraph(&library->global_graphs[global].overal_scale);
		FreeGraph(&library->global_graphs[global].intensity);
		FreeGraph(&library->global_graphs[global].splatter);
		FreeGraph(&library->global_graphs[global].emitter_width);
		FreeGraph(&library->global_graphs[global].emitter_height);
		FreeGraph(&library->global_graphs[global].emitter_depth);

		tfxU32 transform_attributes = effect->transform_attributes;
		FreeGraph(&library->transform_attributes[transform_attributes].roll);
		FreeGraph(&library->transform_attributes[transform_attributes].pitch);
		FreeGraph(&library->transform_attributes[transform_attributes].yaw);
		FreeGraph(&library->transform_attributes[transform_attributes].translation_x);
		FreeGraph(&library->transform_attributes[transform_attributes].translation_y);
		FreeGraph(&library->transform_attributes[transform_attributes].translation_z);
	}

	if (effect->type == tfxEmitterType) {
		tfxU32 transform_attributes = effect->transform_attributes;
		tfxU32 emitter_attributes = effect->emitter_attributes;

		FreeGraph(&library->transform_attributes[transform_attributes].roll);
		FreeGraph(&library->transform_attributes[transform_attributes].pitch);
		FreeGraph(&library->transform_attributes[transform_attributes].yaw);
		FreeGraph(&library->transform_attributes[transform_attributes].translation_x);
		FreeGraph(&library->transform_attributes[transform_attributes].translation_y);
		FreeGraph(&library->transform_attributes[transform_attributes].translation_z);

		FreeGraph(&library->emitter_attributes[emitter_attributes].properties.emission_pitch);
		FreeGraph(&library->emitter_attributes[emitter_attributes].properties.emission_yaw);
		FreeGraph(&library->emitter_attributes[emitter_attributes].properties.emission_range);
		FreeGraph(&library->emitter_attributes[emitter_attributes].properties.splatter);
		FreeGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_width);
		FreeGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_height);
		FreeGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_depth);
		FreeGraph(&library->emitter_attributes[emitter_attributes].properties.extrusion);
		FreeGraph(&library->emitter_attributes[emitter_attributes].properties.arc_size);
		FreeGraph(&library->emitter_attributes[emitter_attributes].properties.arc_offset);

		FreeGraph(&library->emitter_attributes[emitter_attributes].base.life);
		FreeGraph(&library->emitter_attributes[emitter_attributes].base.amount);
		FreeGraph(&library->emitter_attributes[emitter_attributes].base.velocity);
		FreeGraph(&library->emitter_attributes[emitter_attributes].base.width);
		FreeGraph(&library->emitter_attributes[emitter_attributes].base.height);
		FreeGraph(&library->emitter_attributes[emitter_attributes].base.weight);
		FreeGraph(&library->emitter_attributes[emitter_attributes].base.spin);
		FreeGraph(&library->emitter_attributes[emitter_attributes].base.noise_offset);

		FreeGraph(&library->emitter_attributes[emitter_attributes].variation.life);
		FreeGraph(&library->emitter_attributes[emitter_attributes].variation.amount);
		FreeGraph(&library->emitter_attributes[emitter_attributes].variation.velocity);
		FreeGraph(&library->emitter_attributes[emitter_attributes].variation.width);
		FreeGraph(&library->emitter_attributes[emitter_attributes].variation.height);
		FreeGraph(&library->emitter_attributes[emitter_attributes].variation.weight);
		FreeGraph(&library->emitter_attributes[emitter_attributes].variation.spin);
		FreeGraph(&library->emitter_attributes[emitter_attributes].variation.noise_offset);
		FreeGraph(&library->emitter_attributes[emitter_attributes].variation.noise_resolution);
		FreeGraph(&library->emitter_attributes[emitter_attributes].variation.motion_randomness);

		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.width);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.height);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.weight);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.spin);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.stretch);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.red);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.green);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.blue);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.blendfactor);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.intensity);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.direction_turbulance);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.direction);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.noise_resolution);
		FreeGraph(&library->emitter_attributes[emitter_attributes].overtime.motion_randomness);
	}
}

tfxU32 CountAllEffectLookupValues(tfx_effect_emitter_t *effect) {
	tfxU32 count = 0;
	if (effect->type == tfxEffectType) {
		count += CountLibraryGlobalLookUpValues(effect->library, effect->global);
		for (auto &emitter : GetEffectInfo(effect)->sub_effectors) {
			count += CountLibraryEmitterLookUpValues(effect->library, emitter.emitter_attributes);
		}
	}
	else if (effect->type == tfxEmitterType) {
		count += CountLibraryEmitterLookUpValues(effect->library, effect->emitter_attributes);
	}
	return count;
}

void InitialisePathGraphs(tfx_emitter_path_t *path, tfxU32 bucket_size) {
	path->angle_x.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	path->angle_x.type = tfxPath_angle_x;
	path->angle_x.graph_preset = tfxPathDirectionOvertimePreset;
	ResetGraph(&path->angle_x, 0.f, path->angle_x.graph_preset, true, 1.f);
	path->angle_y.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	path->angle_y.type = tfxPath_angle_y;
	path->angle_y.graph_preset = tfxPathDirectionOvertimePreset;
	ResetGraph(&path->angle_y, 0.f, path->angle_y.graph_preset, true, 1.f);
	path->angle_z.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	path->angle_z.type = tfxPath_angle_z;
	path->angle_z.graph_preset = tfxPathDirectionOvertimePreset;
	ResetGraph(&path->angle_z, 0.f, path->angle_z.graph_preset, true, 1.f);
	path->offset_x.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	path->offset_x.type = tfxPath_offset_z;
	path->offset_x.graph_preset = tfxPathTranslationOvertimePreset;
	ResetGraph(&path->offset_x, 0.f, path->offset_x.graph_preset, true, 1.f);
	path->offset_y.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	path->offset_y.type = tfxPath_offset_z;
	path->offset_y.graph_preset = tfxPathTranslationOvertimePreset;
	ResetGraph(&path->offset_y, 0.f, path->offset_y.graph_preset, true, 1.f);
	path->offset_z.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	path->offset_z.type = tfxPath_offset_z;
	path->offset_z.graph_preset = tfxPathTranslationOvertimePreset;
	ResetGraph(&path->offset_z, 0.f, path->offset_z.graph_preset, true, 1.f);
	path->distance.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	path->distance.type = tfxPath_distance;
	path->distance.graph_preset = tfxPathTranslationOvertimePreset;
	ResetGraph(&path->distance, 0.f, path->distance.graph_preset, true, 1.f);
}

void ResetPathGraphs(tfx_emitter_path_t *path, tfx_path_generator_type generator) {
	ResetGraph(&path->angle_x, 0.f, path->angle_x.graph_preset, true, 1.f);
	ResetGraph(&path->angle_y, 0.f, path->angle_y.graph_preset, true, 1.f);
	ResetGraph(&path->angle_z, 0.f, path->angle_z.graph_preset, true, 1.f);
	ResetGraph(&path->offset_x, 0.f, path->offset_x.graph_preset, true, 1.f);
	ResetGraph(&path->offset_y, 0.f, path->offset_y.graph_preset, true, 1.f);
	ResetGraph(&path->offset_z, 0.f, path->offset_z.graph_preset, true, 1.f);
	ResetGraph(&path->distance, 0.f, path->distance.graph_preset, true, 1.f);
	switch (generator) {
	case tfxPathGenerator_spiral:
		if (path->flags & tfxPathFlags_2d) {
			AddGraphNode(&path->angle_x, 1.f, tfxPI2);
			ResetGraph(&path->offset_x, 0.f, path->offset_x.graph_preset, true, 1.f);
			AddGraphNode(&path->offset_x, 1.f, 200.f);
		}
		else {
			AddGraphNode(&path->angle_y, 1.f, tfxPI2);
			ResetGraph(&path->offset_x, 2.f, path->offset_x.graph_preset, true, 1.f);
			AddGraphNode(&path->offset_x, 1.f, 2.f);
			AddGraphNode(&path->offset_y, 1.f, 5.f);
		}
		break;
	case tfxPathGenerator_arc:
		if (path->flags & tfxPathFlags_2d) {
			ResetGraph(&path->distance, 25.f, path->distance.graph_preset, true, 1.f);
			AddGraphNode(&path->angle_x, 1.f, tfx180Radians);
		}
		else {
			ResetGraph(&path->distance, .4f, path->distance.graph_preset, true, 1.f);
			AddGraphNode(&path->angle_x, 1.f, tfx180Radians);
		}
		break;
	case tfxPathGenerator_loop:
		if (path->flags & tfxPathFlags_2d) {
			ResetGraph(&path->distance, 25.f, path->distance.graph_preset, true, 1.f);
			ResetGraph(&path->angle_x, tfx90Radians, path->angle_x.graph_preset, true, 1.f);
			AddGraphNode(&path->angle_x, .25f, tfx90Radians);
			AddGraphNode(&path->angle_x, .75f, -tfx270Radians);
		}
		else {
			ResetGraph(&path->distance, .4f, path->distance.graph_preset, true, 1.f);
			ResetGraph(&path->angle_x, -tfx90Radians, path->angle_x.graph_preset, true, 1.f);
			AddGraphNode(&path->angle_x, .25f, -tfx90Radians);
			AddGraphNode(&path->angle_x, .75f, tfx270Radians);
		}
		break;
	case tfxPathGenerator_s_curve:
		if (path->flags & tfxPathFlags_2d) {
			ResetGraph(&path->distance, 25.f, path->distance.graph_preset, true, 1.f);
			ResetGraph(&path->angle_x, tfx__degrees_to_radians(90.f), path->angle_x.graph_preset, true, 1.f);
			AddGraphNode(&path->angle_x, .25f, tfx__degrees_to_radians(110.f));
			AddGraphNode(&path->angle_x, .75f, tfx__degrees_to_radians(70.f));
			AddGraphNode(&path->angle_x, 1.f, tfx__degrees_to_radians(90.f));
		}
		else {
			ResetGraph(&path->distance, .4f, path->distance.graph_preset, true, 1.f);
			ResetGraph(&path->angle_x, -0.f, path->angle_x.graph_preset, true, 1.f);
			AddGraphNode(&path->angle_x, .25f, tfx__degrees_to_radians(20.f));
			AddGraphNode(&path->angle_x, .75f, tfx__degrees_to_radians(-20.f));
			AddGraphNode(&path->angle_x, 1.f, tfx__degrees_to_radians(0.f));
		}
		break;
	case tfxPathGenerator_bend:
		if (path->flags & tfxPathFlags_2d) {
			ResetGraph(&path->distance, 25.f, path->distance.graph_preset, true, 1.f);
			AddGraphNode(&path->angle_x, .5f, 0.f);
			AddGraphNode(&path->angle_x, .6f, tfx90Radians);
			path->builder_parameters.x = .5f;
			path->builder_parameters.y = .1f;
		}
		else {
			ResetGraph(&path->distance, 0.4f, path->distance.graph_preset, true, 1.f);
			AddGraphNode(&path->angle_x, .5f, 0.f);
			AddGraphNode(&path->angle_x, .6f, tfx90Radians);
			path->builder_parameters.x = .5f;
			path->builder_parameters.y = .1f;
		}
		break;
	default:
		break;
	}
}

void FreePathGraphs(tfx_emitter_path_t *path) {
	FreeGraph(&path->angle_x);
	FreeGraph(&path->angle_y);
	FreeGraph(&path->angle_z);
	FreeGraph(&path->offset_x);
	FreeGraph(&path->offset_y);
	FreeGraph(&path->offset_z);
	FreeGraph(&path->distance);
}

void CopyPathGraphs(tfx_emitter_path_t *src, tfx_emitter_path_t *dst) {
	if (src == dst) return;
	CopyGraphNoLookups(&src->angle_x, &dst->angle_x);
	CopyGraphNoLookups(&src->angle_y, &dst->angle_y);
	CopyGraphNoLookups(&src->angle_z, &dst->angle_z);
	CopyGraphNoLookups(&src->offset_x, &dst->offset_x);
	CopyGraphNoLookups(&src->offset_y, &dst->offset_y);
	CopyGraphNoLookups(&src->offset_z, &dst->offset_z);
	CopyGraphNoLookups(&src->distance, &dst->distance);
}

tfx_emitter_path_t CopyPath(tfx_emitter_path_t *src, const char *name) {
	tfx_emitter_path_t path = {};
	path.flags = src->flags;
	path.name = name;
	path.node_count = src->node_count;
	path.generator_type = src->generator_type;
	path.extrusion_type = src->extrusion_type;
	path.offset = src->offset;
	path.maximum_active_paths = src->maximum_active_paths;
	path.maximum_paths = src->maximum_paths;
	path.rotation_cycle_length = src->rotation_cycle_length;
	path.rotation_stagger = src->rotation_stagger;
	path.rotation_range = src->rotation_range;
	path.rotation_pitch = src->rotation_pitch;
	path.rotation_yaw = src->rotation_yaw;
	InitialisePathGraphs(&path);
	CopyPathGraphs(src, &path);
	return path;
}

tfxU32 CreateEmitterPathAttributes(tfx_effect_emitter_t *emitter, bool add_node) {
	if (emitter->path_attributes == tfxINVALID) {
		tfx_emitter_path_t path = {};
		path.flags = 0;
		path.name = "";
		path.node_count = 32;
		path.extrusion_type = tfxExtrusionArc;
		path.generator_type = tfxPathGenerator_spiral;
		path.maximum_active_paths = 1;
		path.maximum_paths = 1;
		path.offset = { 0 };
		path.rotation_cycle_length = 0.f;
		path.rotation_range = 0.f;
		path.rotation_pitch = 0.f;
		path.rotation_yaw = 0.f;
		path.rotation_stagger = 0.f;
		InitialisePathGraphs(&path);
		ResetGraph(&path.angle_x, 0.f, path.angle_x.graph_preset, add_node, 1.f);
		ResetGraph(&path.angle_y, 0.f, path.angle_y.graph_preset, add_node, 1.f);
		ResetGraph(&path.angle_z, 0.f, path.angle_z.graph_preset, add_node, 1.f);
		ResetGraph(&path.offset_x, 1.f, path.offset_x.graph_preset, add_node, 1.f);
		ResetGraph(&path.offset_y, 0.f, path.offset_y.graph_preset, add_node, 1.f);
		ResetGraph(&path.offset_z, 0.f, path.offset_z.graph_preset, add_node, 1.f);
		ResetGraph(&path.distance, 0.f, path.offset_z.graph_preset, add_node, 1.f);
		emitter->path_attributes = emitter->library->paths.size();
		emitter->library->paths.push_back(path);
	}
	return emitter->path_attributes;
}

tfxU32 AddEmitterPathAttributes(tfx_library_t *library) {
	tfx_emitter_path_t path = {};
	path.flags = 0;
	path.name = "";
	path.node_count = 32;
	path.extrusion_type = tfxExtrusionArc;
	path.generator_type = tfxPathGenerator_spiral;
	path.maximum_active_paths = 1;
	path.maximum_paths = 1;
	path.offset = { 0 };
	path.rotation_cycle_length = 0.f;
	path.rotation_range = 0.f;
	path.rotation_pitch = 0.f;
	path.rotation_yaw = 0.f;
	path.rotation_stagger = 0.f;
	InitialisePathGraphs(&path);
	library->paths.push_back(path);
	return library->paths.size() - 1;
}

float tfx__catmull_rom_segment(tfx_vector_t<tfx_vec4_t> *nodes, float length) {
	tfxU32 i = 0;
	while (length > (*nodes)[i].w && i < nodes->current_size - 3) {
		length -= (*nodes)[i].w;
		i++;
	}
	return (float)i + ((*nodes)[i].w > 0 ? (length / (*nodes)[i].w) : 0.f);
}

void BuildPathNodesComplex(tfx_emitter_path_t *path) {
	//This is currently unused and can probably be removed at some point
	if (!path->node_buffer.capacity) {
		InitPathsSoA3d(&path->node_buffer, &path->node_soa, path->node_count);
	}
	if (path->nodes.current_size != path->node_count) {
		path->nodes.resize(path->node_count);
		if ((tfxU32)path->node_count > path->node_buffer.capacity) {
			GrowArrays(&path->node_buffer, path->node_buffer.capacity, path->node_count, false);
		}
		path->node_buffer.current_size = path->node_count;
	}
	tfx_vector_t<tfx_vec4_t> path_nodes;
	path_nodes.resize(path->node_count);
	float pitch, yaw, roll;
	tfx_mat4_t pitch_mat, yaw_mat, roll_mat, matrix;
	float node_count = (float)path->node_count - 3.f;
	if (path->flags & tfxPathFlags_mode_origin) {
		tfx_vec4_t offset, position;
		float age_inc = 1.f / node_count; float age = 0.f; int i = 1;
		while (i < path->node_count) {
			pitch_mat = tfx__matrix4_rotate_x(GetGraphValue(&path->angle_x, age));
			yaw_mat = tfx__matrix4_rotate_y(GetGraphValue(&path->angle_y, age));
			roll_mat = tfx__matrix4_rotate_z(GetGraphValue(&path->angle_z, age));
			matrix = tfx__transform_matrix4(&yaw_mat, &pitch_mat);
			matrix = tfx__transform_matrix4(&matrix, &roll_mat);
			offset = { GetGraphValue(&path->offset_x, age), GetGraphValue(&path->offset_y, age), GetGraphValue(&path->offset_z, age), 0.f };
			position = tfx__transform_matrix4_vec4(&matrix, offset);
			position += path->offset;
			age += age_inc;
			path_nodes[i++] = position;
		}
		path_nodes[0] = path_nodes[1];
		path_nodes[path->node_count - 1] = path_nodes[path->node_count - 2];
	}
	else if (path->flags & tfxPathFlags_mode_node) {
		tfx_vec4_t distance = { 0.f, GetGraphValue(&path->distance, 0.f), 0.f, 0.f };
		tfx_vec4_t position = tfx__transform_matrix4_vec4(&matrix, distance);
		float age_inc = 1.f / node_count; float age = 0.f; int i = 1;
		while (i < path->node_count) {
			pitch = GetGraphValue(&path->angle_x, age);
			yaw = GetGraphValue(&path->angle_y, age);
			roll = GetGraphValue(&path->angle_z, age);
			pitch_mat = tfx__matrix4_rotate_x(pitch);
			yaw_mat = tfx__matrix4_rotate_y(yaw);
			roll_mat = tfx__matrix4_rotate_z(roll);
			matrix = tfx__transform_matrix4(&yaw_mat, &pitch_mat);
			matrix = tfx__transform_matrix4(&matrix, &roll_mat);
			distance = { 0.f, GetGraphValue(&path->distance, age), 0.f, 0.f };
			position += tfx__transform_matrix4_vec4(&matrix, distance);
			age += age_inc;
			path_nodes[i++] = position + path->offset;
		}
		path_nodes[0] = path_nodes[1];
		path_nodes[path->node_count - 1] = path_nodes[path->node_count - 2];
	}
	if (path->flags & tfxPathFlags_space_nodes_evenly) {
		float length = 0.f;
		for (int i = 0; i != path->node_count - 3; ++i) {
			float step = 0.05f;
			path_nodes[i].w = 0.f;
			for (float t = 0.0f; t <= 1.0f - step; t += step) {
				tfx_vec3_t p1 = tfx__catmull_rom_spline_3d(&path_nodes[i], &path_nodes[i + 1], &path_nodes[i + 2], &path_nodes[i + 3], t);
				tfx_vec3_t p2 = tfx__catmull_rom_spline_3d(&path_nodes[i], &path_nodes[i + 1], &path_nodes[i + 2], &path_nodes[i + 3], t + step);
				tfx_vec3_t segment = p2 - p1;
				path_nodes[i].w += tfx__length_vec3(&segment);
			}
			length += path_nodes[i].w;
		}
		float segment_length = length / (node_count);
		float segment = 0.f;
		int i = 1;
		int c = path->node_count - 1;
		path->nodes[1] = tfx__catmull_rom_spline_3d(&path_nodes[0], &path_nodes[1], &path_nodes[2], &path_nodes[3], 0.f);
		float ni = tfx__catmull_rom_segment(&path_nodes, length);
		ni = ni == (int)ni && ni > 0 ? ni - 0.0001f : ni;
		path->nodes[c - 1] = tfx__catmull_rom_spline_3d(&path_nodes[(int)ni], &path_nodes[(int)ni + 1], &path_nodes[(int)ni + 2], &path_nodes[(int)ni + 3], ni - int(ni));
		path->node_soa.x[1] = path->nodes[1].x;
		path->node_soa.y[1] = path->nodes[1].y;
		path->node_soa.z[1] = path->nodes[1].z;
		path->node_soa.x[c - 1] = path->nodes[c - 1].x;
		path->node_soa.y[c - 1] = path->nodes[c - 1].y;
		path->node_soa.z[c - 1] = path->nodes[c - 1].z;
		while (i != c - 1) {
			ni = tfx__catmull_rom_segment(&path_nodes, segment);
			tfx_vec3_t position = tfx__catmull_rom_spline_3d(&path_nodes[(int)ni], &path_nodes[(int)ni + 1], &path_nodes[(int)ni + 2], &path_nodes[(int)ni + 3], ni - int(ni));
			path->node_soa.x[i] = position.x;
			path->node_soa.y[i] = position.y;
			path->node_soa.z[i] = position.z;
			path->nodes[i++] = position;
			segment += segment_length;
		}
		path->nodes[0] = path->nodes[1];
		path->nodes[c] = path->nodes[c - 1];
		path->node_soa.x[0] = path->node_soa.x[1];
		path->node_soa.y[0] = path->node_soa.y[1];
		path->node_soa.z[0] = path->node_soa.z[1];
		path->node_soa.length[0] = path->node_soa.length[1];
		path->node_soa.x[c] = path->node_soa.x[c - 1];
		path->node_soa.y[c] = path->node_soa.y[c - 1];
		path->node_soa.z[c] = path->node_soa.z[c - 1];
		path->node_soa.length[c] = path->node_soa.length[c - 1];
	}
	else {
		for (int i = 0; i != path->node_count; ++i) {
			path->nodes[i] = path_nodes[i];
			path->node_soa.x[i] = path_nodes[i].x;
			path->node_soa.y[i] = path_nodes[i].y;
			path->node_soa.z[i] = path_nodes[i].z;
			path->node_soa.length[i] = path_nodes[i].w;
		}
	}
	path_nodes.free_all();
}

void BuildPathNodes3d(tfx_emitter_path_t *path) {
	if (!path->node_buffer.capacity) {
		InitPathsSoA3d(&path->node_buffer, &path->node_soa, path->node_count);
	}
	if (path->nodes.current_size != path->node_count) {
		path->nodes.resize(path->node_count);
		if ((tfxU32)path->node_count > path->node_buffer.capacity) {
			GrowArrays(&path->node_buffer, path->node_buffer.capacity, path->node_count, false);
		}
		path->node_buffer.current_size = path->node_count;
	}
	tfx_vector_t<tfx_vec4_t> path_nodes;
	path_nodes.resize(path->node_count);
	tfx_mat4_t matrix;
	float node_count = (float)path->node_count - 1.f;
	if (path->generator_type == tfxPathGenerator_spiral) {
		tfx_vec4_t offset, position;
		float age_inc = 1.f / node_count; float age = 0.f; int i = 0;
		while (i < path->node_count) {
			matrix = tfx__matrix4_rotate_y(GetGraphValue(&path->angle_y, age));
			offset = { GetGraphValue(&path->offset_x, age), GetGraphValue(&path->offset_y, age), 0.f, 0.f };
			position = tfx__transform_matrix4_vec4(&matrix, offset);
			position += path->offset;
			age += age_inc;
			path_nodes[i++] = position;
		}
	}
	else if (path->generator_type == tfxPathGenerator_arc) {
		tfx_vec4_t offset, position;
		float age_inc = 1.f / node_count; float age = 0.f; int i = 0;
		tfx_vec4_t distance = {};
		while (i < path->node_count) {
			matrix = tfx__matrix4_rotate_x(GetGraphValue(&path->angle_x, age));
			distance = { 0.f, GetGraphValue(&path->distance, age), 0.f, 0.f };
			position += tfx__transform_matrix4_vec4(&matrix, distance);
			age += age_inc;
			path_nodes[i++] = position + path->offset;
		}
	}
	else if (path->generator_type == tfxPathGenerator_loop) {
		tfx_vec4_t offset, position;
		float age_inc = 1.f / node_count; float age = 0.f; int i = 0;
		tfx_vec4_t distance = {};
		tfx_mat4_t z_mat;
		while (i < path->node_count) {
			matrix = tfx__matrix4_rotate_x(GetGraphValue(&path->angle_x, age));
			z_mat = tfx__matrix4_rotate_z(GetGraphValue(&path->angle_z, age));
			matrix = tfx__transform_matrix4(&matrix, &z_mat);
			distance = { 0.f, GetGraphValue(&path->distance, age), 0.f, 0.f };
			position += tfx__transform_matrix4_vec4(&matrix, distance);
			age += age_inc;
			path_nodes[i++] = position + path->offset;
		}
	}
	else if (path->generator_type == tfxPathGenerator_s_curve) {
		tfx_vec4_t offset, position;
		float age_inc = 1.f / node_count; float age = 0.f; int i = 0;
		tfx_vec4_t distance = {};
		tfx_mat4_t z_mat;
		while (i < path->node_count) {
			matrix = tfx__matrix4_rotate_x(GetGraphValue(&path->angle_x, age));
			z_mat = tfx__matrix4_rotate_z(GetGraphValue(&path->angle_z, age));
			matrix = tfx__transform_matrix4(&matrix, &z_mat);
			distance = { 0.f, GetGraphValue(&path->distance, age), 0.f, 0.f };
			position += tfx__transform_matrix4_vec4(&matrix, distance);
			age += age_inc;
			path_nodes[i++] = position + path->offset;
		}
	}
	else if (path->generator_type == tfxPathGenerator_bend) {
		tfx_vec4_t offset, position;
		float age_inc = 1.f / node_count; float age = 0.f; int i = 0;
		tfx_vec4_t distance = {};
		tfx_mat4_t z_mat;
		if (path->builder_parameters.x + path->builder_parameters.y == 0) {
			for (tfxBucketLoop(path->angle_x.nodes)) {
				if (i == 1) {
					path->builder_parameters.x = path->angle_x.nodes[i].frame;
				}
				else if (i == 2) {
					path->builder_parameters.y = path->angle_x.nodes[i].frame - path->angle_x.nodes[1].frame;
				}
			}
		}
		while (i < path->node_count) {
			matrix = tfx__matrix4_rotate_x(GetGraphValue(&path->angle_x, age));
			distance = { 0.f, GetGraphValue(&path->distance, age), 0.f, 0.f };
			position += tfx__transform_matrix4_vec4(&matrix, distance);
			age += age_inc;
			path_nodes[i++] = position + path->offset;
		}
	}
	else if (path->generator_type == tfxPathGenerator_free_mode_origin) {
		tfx_vec4_t offset, position;
		tfx_mat4_t pitch_mat, yaw_mat, roll_mat, matrix;
		float age_inc = 1.f / node_count; float age = 0.f; int i = 0;
		while (i < path->node_count) {
			pitch_mat = tfx__matrix4_rotate_x(GetGraphValue(&path->angle_x, age));
			yaw_mat = tfx__matrix4_rotate_y(GetGraphValue(&path->angle_y, age));
			roll_mat = tfx__matrix4_rotate_z(GetGraphValue(&path->angle_z, age));
			matrix = tfx__transform_matrix4(&yaw_mat, &pitch_mat);
			matrix = tfx__transform_matrix4(&matrix, &roll_mat);
			offset = { GetGraphValue(&path->offset_x, age), GetGraphValue(&path->offset_y, age), GetGraphValue(&path->offset_z, age), 0.f };
			position = tfx__transform_matrix4_vec4(&matrix, offset);
			position += path->offset;
			age += age_inc;
			path_nodes[i++] = position;
		}
	}
	else if (path->generator_type == tfxPathGenerator_free_mode_distance) {
		tfx_mat4_t pitch_mat, yaw_mat, roll_mat, matrix;
		float pitch, yaw, roll;
		tfx_vec4_t distance = { 0.f, GetGraphValue(&path->distance, 0.f), 0.f, 0.f };
		tfx_vec4_t position = tfx__transform_matrix4_vec4(&matrix, distance);
		float age_inc = 1.f / node_count; float age = 0.f; int i = 0;
		while (i < path->node_count) {
			pitch = GetGraphValue(&path->angle_x, age);
			yaw = GetGraphValue(&path->angle_y, age);
			roll = GetGraphValue(&path->angle_z, age);
			pitch_mat = tfx__matrix4_rotate_x(pitch);
			yaw_mat = tfx__matrix4_rotate_y(yaw);
			roll_mat = tfx__matrix4_rotate_z(roll);
			matrix = tfx__transform_matrix4(&yaw_mat, &pitch_mat);
			matrix = tfx__transform_matrix4(&matrix, &roll_mat);
			distance = { 0.f, GetGraphValue(&path->distance, age), 0.f, 0.f };
			position += tfx__transform_matrix4_vec4(&matrix, distance);
			age += age_inc;
			path_nodes[i++] = position + path->offset;
		}
	}
	if (path->flags & tfxPathFlags_space_nodes_evenly) {
		float length = 0.f;
		for (int i = 0; i != path->node_count - 3; ++i) {
			float step = 0.05f;
			path_nodes[i].w = 0.f;
			for (float t = 0.0f; t <= 1.0f - step; t += step) {
				tfx_vec3_t p1 = tfx__catmull_rom_spline_3d(&path_nodes[i], &path_nodes[i + 1], &path_nodes[i + 2], &path_nodes[i + 3], t);
				tfx_vec3_t p2 = tfx__catmull_rom_spline_3d(&path_nodes[i], &path_nodes[i + 1], &path_nodes[i + 2], &path_nodes[i + 3], t + step);
				tfx_vec3_t segment = p2 - p1;
				path_nodes[i].w += tfx__length_vec3(&segment);
			}
			length += path_nodes[i].w;
		}
		float segment_length = length / node_count;
		float segment = 0.f;
		int i = 0;
		int c = path->node_count;
		path->nodes[1] = tfx__catmull_rom_spline_3d(&path_nodes[0], &path_nodes[1], &path_nodes[2], &path_nodes[3], 0.f);
		while (i != c) {
			float ni = tfx__catmull_rom_segment(&path_nodes, segment);
			if (ni >= path->node_count - 3) {
				ni = (float)path->node_count - 3.f - 0.0001f;
			}
			tfx_vec3_t position = tfx__catmull_rom_spline_3d(&path_nodes[(int)ni], &path_nodes[(int)ni + 1], &path_nodes[(int)ni + 2], &path_nodes[(int)ni + 3], ni - int(ni));
			if (path->flags & tfxPathFlags_reverse_direction) {
				int node_count = path->node_count - 1;
				path->node_soa.x[node_count - i] = position.x;
				path->node_soa.y[node_count - i] = position.y;
				path->node_soa.z[node_count - i] = position.z;
				path->nodes[node_count - i] = position;
			}
			else {
				path->node_soa.x[i] = position.x;
				path->node_soa.y[i] = position.y;
				path->node_soa.z[i] = position.z;
				path->nodes[i] = position;
			}
			segment += segment_length;
			i++;
		}
	}
	else {
		if (path->flags & tfxPathFlags_reverse_direction) {
			int node_count = path->node_count - 1;
			for (int i = 0; i != path->node_count; ++i) {
				path->nodes[node_count - i] = path_nodes[i];
				path->node_soa.x[node_count - i] = path_nodes[i].x;
				path->node_soa.y[node_count - i] = path_nodes[i].y;
				path->node_soa.z[node_count - i] = path_nodes[i].z;
				path->node_soa.length[node_count - i] = path_nodes[i].w;
			}
		}
		else {
			for (int i = 0; i != path->node_count; ++i) {
				path->nodes[i] = path_nodes[i];
				path->node_soa.x[i] = path_nodes[i].x;
				path->node_soa.y[i] = path_nodes[i].y;
				path->node_soa.z[i] = path_nodes[i].z;
				path->node_soa.length[i] = path_nodes[i].w;
			}
		}
	}
	path_nodes.free_all();
}

void BuildPathNodes2d(tfx_emitter_path_t *path) {
	if (!path->node_buffer.capacity) {
		InitPathsSoA2d(&path->node_buffer, &path->node_soa, path->node_count);
	}
	if (path->nodes.current_size != path->node_count) {
		path->nodes.resize(path->node_count);
		if ((tfxU32)path->node_count > path->node_buffer.capacity) {
			GrowArrays(&path->node_buffer, path->node_buffer.capacity, path->node_count, false);
		}
		path->node_buffer.current_size = path->node_count;
	}
	tfx_vector_t<tfx_vec4_t> path_nodes;
	path_nodes.resize(path->node_count);
	tfx_mat4_t matrix;
	float node_count = (float)path->node_count - 1.f;
	if (path->generator_type == tfxPathGenerator_spiral) {
		tfx_vec2_t position;
		float age_inc = 1.f / node_count; float age = 0.f; int i = 0;
		while (i < path->node_count) {
			float angle = GetGraphValue(&path->angle_x, age);
			float radius = GetGraphValue(&path->offset_x, age);
			position = { sinf(angle) * radius, -cosf(angle) * radius };
			age += age_inc;
			path_nodes[i++] = position + path->offset.xy();
		}
	}
	else if (path->generator_type == tfxPathGenerator_arc) {
		tfx_vec2_t offset, position;
		float age_inc = 1.f / node_count; float age = 0.f; int i = 0;
		float distance;
		while (i < path->node_count) {
			float angle = GetGraphValue(&path->angle_x, age);
			distance = GetGraphValue(&path->distance, age);
			position += {sinf(angle) *distance, -cosf(angle) * distance};
			age += age_inc;
			path_nodes[i++] = position + path->offset.xy();
		}
	}
	else if (path->generator_type == tfxPathGenerator_loop) {
		tfx_vec2_t offset, position;
		float age_inc = 1.f / node_count; float age = 0.f; int i = 0;
		float distance;
		while (i < path->node_count) {
			float angle = GetGraphValue(&path->angle_x, age);
			distance = GetGraphValue(&path->distance, age);
			position += {sinf(angle) *distance, -cosf(angle) * distance};
			age += age_inc;
			path_nodes[i++] = position + path->offset.xy();
		}
	}
	else if (path->generator_type == tfxPathGenerator_s_curve) {
		tfx_vec2_t offset, position;
		float age_inc = 1.f / node_count; float age = 0.f; int i = 0;
		float distance;
		while (i < path->node_count) {
			float angle = GetGraphValue(&path->angle_x, age);
			distance = GetGraphValue(&path->distance, age);
			position += {sinf(angle) *distance, -cosf(angle) * distance};
			age += age_inc;
			path_nodes[i++] = position + path->offset.xy();
		}
	}
	else if (path->generator_type == tfxPathGenerator_bend) {
		tfx_vec2_t offset, position;
		float age_inc = 1.f / node_count; float age = 0.f; int i = 0;
		float distance;
		if (path->builder_parameters.x + path->builder_parameters.y == 0) {
			for (tfxBucketLoop(path->angle_x.nodes)) {
				if (i == 1) {
					path->builder_parameters.x = path->angle_x.nodes[i].frame;
				}
				else if (i == 2) {
					path->builder_parameters.y = path->angle_x.nodes[i].frame - path->angle_x.nodes[1].frame;
				}
			}
		}
		while (i < path->node_count) {
			float angle = GetGraphValue(&path->angle_x, age);
			distance = GetGraphValue(&path->distance, age);
			position += {sinf(angle) *distance, -cosf(angle) * distance};
			age += age_inc;
			path_nodes[i++] = position + path->offset.xy();
		}
	}
	else if (path->generator_type == tfxPathGenerator_free_mode_origin) {
		tfx_vec2_t offset, position;
		float age_inc = 1.f / node_count; float age = 0.f; int i = 0;
		while (i < path->node_count) {
			float angle = GetGraphValue(&path->angle_x, age);
			offset = { GetGraphValue(&path->offset_x, age), GetGraphValue(&path->offset_y, age) };
			position = { sinf(angle) * offset.x, -cosf(angle) * offset.y };
			age += age_inc;
			path_nodes[i++] = position;
		}
	}
	else if (path->generator_type == tfxPathGenerator_free_mode_distance) {
		float distance = GetGraphValue(&path->distance, 0.f);
		tfx_vec2_t position = { distance, 0.f };
		float age_inc = 1.f / node_count; float age = 0.f; int i = 0;
		while (i < path->node_count) {
			float angle = GetGraphValue(&path->angle_x, age);
			distance = GetGraphValue(&path->distance, age);
			position += {sinf(angle) *distance, -cosf(angle) * distance};
			age += age_inc;
			path_nodes[i++] = position + path->offset.xy();
		}
	}
	if (path->flags & tfxPathFlags_space_nodes_evenly) {
		float length = 0.f;
		for (int i = 0; i != path->node_count - 3; ++i) {
			float step = 0.05f;
			path_nodes[i].w = 0.f;
			for (float t = 0.0f; t <= 1.0f - step; t += step) {
				tfx_vec2_t p1 = tfx__catmull_rom_spline_2d(&path_nodes[i], &path_nodes[i + 1], &path_nodes[i + 2], &path_nodes[i + 3], t);
				tfx_vec2_t p2 = tfx__catmull_rom_spline_2d(&path_nodes[i], &path_nodes[i + 1], &path_nodes[i + 2], &path_nodes[i + 3], t + step);
				tfx_vec2_t segment = p2 - p1;
				path_nodes[i].w += tfx__vec2_length_fast(&segment);
			}
			length += path_nodes[i].w;
		}
		float segment_length = length / node_count;
		float segment = 0.f;
		int i = 0;
		int c = path->node_count;
		path->nodes[1] = tfx__catmull_rom_spline_3d(&path_nodes[0], &path_nodes[1], &path_nodes[2], &path_nodes[3], 0.f);
		while (i != c) {
			float ni = tfx__catmull_rom_segment(&path_nodes, segment);
			if (ni >= path->node_count - 3) {
				ni = (float)path->node_count - 3.f - 0.0001f;
			}
			tfx_vec2_t position = tfx__catmull_rom_spline_2d(&path_nodes[(int)ni], &path_nodes[(int)ni + 1], &path_nodes[(int)ni + 2], &path_nodes[(int)ni + 3], ni - int(ni));
			if (path->flags & tfxPathFlags_reverse_direction) {
				int node_count = path->node_count - 1;
				path->node_soa.x[node_count - i] = position.x;
				path->node_soa.y[node_count - i] = position.y;
				path->nodes[node_count - i] = position;
			}
			else {
				path->node_soa.x[i] = position.x;
				path->node_soa.y[i] = position.y;
				path->nodes[i] = position;
			}
			segment += segment_length;
			i++;
		}
	}
	else {
		if (path->flags & tfxPathFlags_reverse_direction) {
			int node_count = path->node_count - 1;
			for (int i = 0; i != path->node_count; ++i) {
				path->nodes[node_count - i] = path_nodes[i];
				path->node_soa.x[node_count - i] = path_nodes[i].x;
				path->node_soa.y[node_count - i] = path_nodes[i].y;
				path->node_soa.length[node_count - i] = path_nodes[i].w;
			}
		}
		else {
			for (int i = 0; i != path->node_count; ++i) {
				path->nodes[i] = path_nodes[i];
				path->node_soa.x[i] = path_nodes[i].x;
				path->node_soa.y[i] = path_nodes[i].y;
				path->node_soa.length[i] = path_nodes[i].w;
			}
		}
	}
	path_nodes.free_all();
}

void InitialiseGlobalAttributes(tfx_global_attributes_t *attributes, tfxU32 bucket_size) {
	attributes->life.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->amount.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->velocity.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->noise.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->width.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->height.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->weight.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->spin.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->pitch_spin.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->yaw_spin.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->stretch.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->overal_scale.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->intensity.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->splatter.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->emitter_width.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->emitter_height.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->emitter_depth.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
}

void FreeGlobalAttributes(tfx_global_attributes_t *attributes) {
	FreeGraph(&attributes->life);
	FreeGraph(&attributes->amount);
	FreeGraph(&attributes->velocity);
	FreeGraph(&attributes->noise);
	FreeGraph(&attributes->width);
	FreeGraph(&attributes->height);
	FreeGraph(&attributes->weight);
	FreeGraph(&attributes->spin);
	FreeGraph(&attributes->pitch_spin);
	FreeGraph(&attributes->yaw_spin);
	FreeGraph(&attributes->stretch);
	FreeGraph(&attributes->overal_scale);
	FreeGraph(&attributes->intensity);
	FreeGraph(&attributes->splatter);
	FreeGraph(&attributes->emitter_width);
	FreeGraph(&attributes->emitter_height);
	FreeGraph(&attributes->emitter_depth);
}

void CopyGlobalAttributesNoLookups(tfx_global_attributes_t *src, tfx_global_attributes_t *dst) {
	if (src == dst) return;
	CopyGraphNoLookups(&src->life, &dst->life);
	CopyGraphNoLookups(&src->amount, &dst->amount);
	CopyGraphNoLookups(&src->velocity, &dst->velocity);
	CopyGraphNoLookups(&src->noise, &dst->noise);
	CopyGraphNoLookups(&src->width, &dst->width);
	CopyGraphNoLookups(&src->height, &dst->height);
	CopyGraphNoLookups(&src->weight, &dst->weight);
	CopyGraphNoLookups(&src->spin, &dst->spin);
	CopyGraphNoLookups(&src->pitch_spin, &dst->pitch_spin);
	CopyGraphNoLookups(&src->yaw_spin, &dst->yaw_spin);
	CopyGraphNoLookups(&src->stretch, &dst->stretch);
	CopyGraphNoLookups(&src->overal_scale, &dst->overal_scale);
	CopyGraphNoLookups(&src->intensity, &dst->intensity);
	CopyGraphNoLookups(&src->splatter, &dst->splatter);
	CopyGraphNoLookups(&src->emitter_width, &dst->emitter_width);
	CopyGraphNoLookups(&src->emitter_height, &dst->emitter_height);
	CopyGraphNoLookups(&src->emitter_depth, &dst->emitter_depth);
}

void CopyGlobalAttributes(tfx_global_attributes_t *src, tfx_global_attributes_t *dst) {
	if (src == dst) return;
	CopyGraph(&src->life, &dst->life);
	CopyGraph(&src->amount, &dst->amount);
	CopyGraph(&src->velocity, &dst->velocity);
	CopyGraph(&src->noise, &dst->noise);
	CopyGraph(&src->width, &dst->width);
	CopyGraph(&src->height, &dst->height);
	CopyGraph(&src->weight, &dst->weight);
	CopyGraph(&src->spin, &dst->spin);
	CopyGraph(&src->pitch_spin, &dst->pitch_spin);
	CopyGraph(&src->yaw_spin, &dst->yaw_spin);
	CopyGraph(&src->stretch, &dst->stretch);
	CopyGraph(&src->overal_scale, &dst->overal_scale);
	CopyGraph(&src->intensity, &dst->intensity);
	CopyGraph(&src->splatter, &dst->splatter);
	CopyGraph(&src->emitter_width, &dst->emitter_width);
	CopyGraph(&src->emitter_height, &dst->emitter_height);
	CopyGraph(&src->emitter_depth, &dst->emitter_depth);
}

void InitialiseTransformAttributes(tfx_transform_attributes_t *attributes, tfxU32 bucket_size) {
	attributes->roll.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->pitch.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->yaw.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->translation_x.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->translation_y.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->translation_z.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
}

void FreeTransformAttributes(tfx_transform_attributes_t *attributes) {
	FreeGraph(&attributes->roll);
	FreeGraph(&attributes->pitch);
	FreeGraph(&attributes->yaw);
	FreeGraph(&attributes->translation_x);
	FreeGraph(&attributes->translation_y);
	FreeGraph(&attributes->translation_z);
}

void CopyTransformAttributesNoLookups(tfx_transform_attributes_t *src, tfx_transform_attributes_t *dst) {
	if (src == dst) return;
	CopyGraphNoLookups(&src->roll, &dst->roll);
	CopyGraphNoLookups(&src->pitch, &dst->pitch);
	CopyGraphNoLookups(&src->yaw, &dst->yaw);
	CopyGraphNoLookups(&src->translation_x, &dst->translation_x);
	CopyGraphNoLookups(&src->translation_y, &dst->translation_y);
	CopyGraphNoLookups(&src->translation_z, &dst->translation_z);
}

void CopyTransformAttributes(tfx_transform_attributes_t *src, tfx_transform_attributes_t *dst) {
	if (src == dst) return;
	CopyGraph(&src->roll, &dst->roll);
	CopyGraph(&src->pitch, &dst->pitch);
	CopyGraph(&src->yaw, &dst->yaw);
	CopyGraph(&src->translation_x, &dst->translation_x);
	CopyGraph(&src->translation_y, &dst->translation_y);
	CopyGraph(&src->translation_z, &dst->translation_z);
}

bool HasTranslationKeyframes(tfx_transform_attributes_t *graphs) {
	return graphs->translation_x.nodes.size() || graphs->translation_y.nodes.size() || graphs->translation_z.nodes.size();
}

void AddTranslationNodes(tfx_transform_attributes_t *keyframes, float frame) {
	if (keyframes->translation_x.nodes.size()) {
		if (!HasNodeAtFrame(&keyframes->translation_x, frame))
			AddGraphCoordNode(&keyframes->translation_x, frame, 0.f);
		if (!HasNodeAtFrame(&keyframes->translation_y, frame))
			AddGraphCoordNode(&keyframes->translation_y, frame, 0.f);
		if (!HasNodeAtFrame(&keyframes->translation_z, frame))
			AddGraphCoordNode(&keyframes->translation_z, frame, 0.f);
	}
	else {
		AddGraphCoordNode(&keyframes->translation_x, 0.f, 0.f);
		AddGraphCoordNode(&keyframes->translation_y, 0.f, 0.f);
		AddGraphCoordNode(&keyframes->translation_z, 0.f, 0.f);
		if (frame != 0) {
			AddGraphCoordNode(&keyframes->translation_x, frame, 0.f);
			AddGraphCoordNode(&keyframes->translation_y, frame, 0.f);
			AddGraphCoordNode(&keyframes->translation_z, frame, 0.f);
		}
	}
}

void InitialisePropertyAttributes(tfx_property_attributes_t *attributes, tfxU32 bucket_size) {
	attributes->emission_pitch.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->emission_yaw.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->emission_range.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->splatter.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->emitter_width.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->emitter_height.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->emitter_depth.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->extrusion.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->arc_size.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->arc_offset.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
}

void FreePropertyAttributes(tfx_property_attributes_t *attributes) {
	FreeGraph(&attributes->emission_pitch);
	FreeGraph(&attributes->emission_yaw);
	FreeGraph(&attributes->emission_range);
	FreeGraph(&attributes->splatter);
	FreeGraph(&attributes->emitter_width);
	FreeGraph(&attributes->emitter_height);
	FreeGraph(&attributes->emitter_depth);
	FreeGraph(&attributes->extrusion);
	FreeGraph(&attributes->arc_size);
	FreeGraph(&attributes->arc_offset);
}

void CopyPropertyAttributesNoLookups(tfx_property_attributes_t *src, tfx_property_attributes_t *dst) {
	if (src == dst) return;
	CopyGraphNoLookups(&src->emission_pitch, &dst->emission_pitch);
	CopyGraphNoLookups(&src->emission_yaw, &dst->emission_yaw);
	CopyGraphNoLookups(&src->emission_range, &dst->emission_range);
	CopyGraphNoLookups(&src->splatter, &dst->splatter);
	CopyGraphNoLookups(&src->emitter_width, &dst->emitter_width);
	CopyGraphNoLookups(&src->emitter_height, &dst->emitter_height);
	CopyGraphNoLookups(&src->emitter_depth, &dst->emitter_depth);
	CopyGraphNoLookups(&src->extrusion, &dst->extrusion);
	CopyGraphNoLookups(&src->arc_size, &dst->arc_size);
	CopyGraphNoLookups(&src->arc_offset, &dst->arc_offset);
}

void CopyPropertyAttributes(tfx_property_attributes_t *src, tfx_property_attributes_t *dst) {
	if (src == dst) return;
	CopyGraph(&src->emission_pitch, &dst->emission_pitch);
	CopyGraph(&src->emission_yaw, &dst->emission_yaw);
	CopyGraph(&src->emission_range, &dst->emission_range);
	CopyGraph(&src->splatter, &dst->splatter);
	CopyGraph(&src->emitter_width, &dst->emitter_width);
	CopyGraph(&src->emitter_height, &dst->emitter_height);
	CopyGraph(&src->emitter_depth, &dst->emitter_depth);
	CopyGraph(&src->extrusion, &dst->extrusion);
	CopyGraph(&src->arc_size, &dst->arc_size);
	CopyGraph(&src->arc_offset, &dst->arc_offset);
}

void InitialiseBaseAttributes(tfx_base_attributes_t *attributes, tfxU32 bucket_size) {
	attributes->life.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->amount.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->velocity.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->width.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->height.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->weight.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->spin.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->pitch_spin.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->yaw_spin.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->noise_offset.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
}

void InitialiseVariationAttributes(tfx_variation_attributes_t *attributes, tfxU32 bucket_size) {
	attributes->life.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->amount.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->velocity.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->width.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->height.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->weight.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->spin.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->pitch_spin.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->yaw_spin.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->noise_offset.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->noise_resolution.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->motion_randomness.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
}

void InitialiseOvertimeAttributes(tfx_overtime_attributes_t *attributes, tfxU32 bucket_size) {
	attributes->velocity.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->width.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->height.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->weight.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->spin.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->pitch_spin.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->yaw_spin.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->stretch.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->red.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->blue.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->green.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->red_hint.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->blue_hint.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->green_hint.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->blendfactor.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->blendfactor_hint.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->velocity_turbulance.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->direction_turbulance.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->velocity_adjuster.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->intensity.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->alpha_sharpness.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->curved_alpha.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->direction.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->noise_resolution.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->motion_randomness.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
}

void InitialiseFactorAttributes(tfx_factor_attributes_t *attributes, tfxU32 bucket_size) {
	attributes->life.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->velocity.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->size.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
	attributes->intensity.nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
}

void FreeOvertimeAttributes(tfx_overtime_attributes_t *attributes) {
	FreeGraph(&attributes->velocity);
	FreeGraph(&attributes->width);
	FreeGraph(&attributes->height);
	FreeGraph(&attributes->weight);
	FreeGraph(&attributes->spin);
	FreeGraph(&attributes->pitch_spin);
	FreeGraph(&attributes->yaw_spin);
	FreeGraph(&attributes->stretch);
	FreeGraph(&attributes->red);
	FreeGraph(&attributes->green);
	FreeGraph(&attributes->blue);
	FreeGraph(&attributes->red_hint);
	FreeGraph(&attributes->green_hint);
	FreeGraph(&attributes->blue_hint);
	FreeGraph(&attributes->blendfactor);
	FreeGraph(&attributes->blendfactor_hint);
	FreeGraph(&attributes->velocity);
	FreeGraph(&attributes->direction_turbulance);
	FreeGraph(&attributes->velocity_adjuster);
	FreeGraph(&attributes->intensity);
	FreeGraph(&attributes->alpha_sharpness);
	FreeGraph(&attributes->curved_alpha);
	FreeGraph(&attributes->direction);
	FreeGraph(&attributes->noise_resolution);
	FreeGraph(&attributes->motion_randomness);
	attributes->color_ramp_bitmap_indexes[0] = 0;
	attributes->color_ramp_bitmap_indexes[1] = 0;
}

void CopyOvertimeAttributesNoLookups(tfx_overtime_attributes_t *src, tfx_overtime_attributes_t *dst) {
	if (src == dst) return;
	CopyGraphNoLookups(&src->velocity, &dst->velocity);
	CopyGraphNoLookups(&src->width, &dst->width);
	CopyGraphNoLookups(&src->height, &dst->height);
	CopyGraphNoLookups(&src->weight, &dst->weight);
	CopyGraphNoLookups(&src->spin, &dst->spin);
	CopyGraphNoLookups(&src->pitch_spin, &dst->pitch_spin);
	CopyGraphNoLookups(&src->yaw_spin, &dst->yaw_spin);
	CopyGraphNoLookups(&src->stretch, &dst->stretch);
	CopyGraphNoLookups(&src->red, &dst->red);
	CopyGraphNoLookups(&src->green, &dst->green);
	CopyGraphNoLookups(&src->blue, &dst->blue);
	CopyGraphNoLookups(&src->red_hint, &dst->red_hint);
	CopyGraphNoLookups(&src->green_hint, &dst->green_hint);
	CopyGraphNoLookups(&src->blue_hint, &dst->blue_hint);
	CopyGraphNoLookups(&src->blendfactor, &dst->blendfactor);
	CopyGraphNoLookups(&src->blendfactor_hint, &dst->blendfactor_hint);
	CopyGraphNoLookups(&src->velocity_turbulance, &dst->velocity_turbulance);
	CopyGraphNoLookups(&src->direction_turbulance, &dst->direction_turbulance);
	CopyGraphNoLookups(&src->velocity_adjuster, &dst->velocity_adjuster);
	CopyGraphNoLookups(&src->intensity, &dst->intensity);
	CopyGraphNoLookups(&src->alpha_sharpness, &dst->alpha_sharpness);
	CopyGraphNoLookups(&src->curved_alpha, &dst->curved_alpha);
	CopyGraphNoLookups(&src->direction, &dst->direction);
	CopyGraphNoLookups(&src->noise_resolution, &dst->noise_resolution);
	CopyGraphNoLookups(&src->motion_randomness, &dst->motion_randomness);
	dst->color_ramps[0] = src->color_ramps[0];
	dst->color_ramps[1] = src->color_ramps[1];
	dst->color_ramp_bitmap_indexes[0] = src->color_ramp_bitmap_indexes[0];
	dst->color_ramp_bitmap_indexes[1] = src->color_ramp_bitmap_indexes[1];
	tfxUnFlagColorRampIDAsEdited(dst->color_ramp_bitmap_indexes[0]);
	tfxUnFlagColorRampIDAsEdited(dst->color_ramp_bitmap_indexes[1]);
}

void CopyOvertimeAttributes(tfx_overtime_attributes_t *src, tfx_overtime_attributes_t *dst) {
	if (src == dst) return;
	CopyGraph(&src->velocity, &dst->velocity);
	CopyGraph(&src->width, &dst->width);
	CopyGraph(&src->height, &dst->height);
	CopyGraph(&src->weight, &dst->weight);
	CopyGraph(&src->spin, &dst->spin);
	CopyGraph(&src->pitch_spin, &dst->pitch_spin);
	CopyGraph(&src->yaw_spin, &dst->yaw_spin);
	CopyGraph(&src->stretch, &dst->stretch);
	CopyGraphColor(src, dst);
	CopyGraphColorHint(src, dst);
	CopyGraph(&src->velocity_turbulance, &dst->velocity_turbulance);
	CopyGraph(&src->direction_turbulance, &dst->direction_turbulance);
	CopyGraph(&src->velocity_adjuster, &dst->velocity_adjuster);
	CopyGraph(&src->intensity, &dst->intensity);
	CopyGraph(&src->alpha_sharpness, &dst->alpha_sharpness);
	CopyGraph(&src->curved_alpha, &dst->curved_alpha);
	CopyGraph(&src->direction, &dst->direction);
	CopyGraph(&src->noise_resolution, &dst->noise_resolution);
	CopyGraph(&src->motion_randomness, &dst->motion_randomness);
	dst->color_ramp_bitmap_indexes[0] = src->color_ramp_bitmap_indexes[0];
	dst->color_ramp_bitmap_indexes[1] = src->color_ramp_bitmap_indexes[1];
	tfxUnFlagColorRampIDAsEdited(dst->color_ramp_bitmap_indexes[0]);
	tfxUnFlagColorRampIDAsEdited(dst->color_ramp_bitmap_indexes[1]);
}

void FreeFactorAttributes(tfx_factor_attributes_t *attributes) {
	FreeGraph(&attributes->life);
	FreeGraph(&attributes->velocity);
	FreeGraph(&attributes->size);
	FreeGraph(&attributes->intensity);
}

void CopyFactorAttributesNoLookups(tfx_factor_attributes_t *src, tfx_factor_attributes_t *dst) {
	if (src == dst) return;
	CopyGraphNoLookups(&src->life, &dst->life);
	CopyGraphNoLookups(&src->velocity, &dst->velocity);
	CopyGraphNoLookups(&src->size, &dst->size);
	CopyGraphNoLookups(&src->intensity, &dst->intensity);
}

void CopyFactorAttributes(tfx_factor_attributes_t *src, tfx_factor_attributes_t *dst) {
	if (src == dst) return;
	CopyGraph(&src->life, &dst->life);
	CopyGraph(&src->velocity, &dst->velocity);
	CopyGraph(&src->size, &dst->size);
	CopyGraph(&src->intensity, &dst->intensity);
}

void FreeVariationAttributes(tfx_variation_attributes_t *attributes) {
	FreeGraph(&attributes->life);
	FreeGraph(&attributes->amount);
	FreeGraph(&attributes->velocity);
	FreeGraph(&attributes->width);
	FreeGraph(&attributes->height);
	FreeGraph(&attributes->weight);
	FreeGraph(&attributes->spin);
	FreeGraph(&attributes->pitch_spin);
	FreeGraph(&attributes->yaw_spin);
	FreeGraph(&attributes->noise_offset);
	FreeGraph(&attributes->noise_resolution);
	FreeGraph(&attributes->motion_randomness);
}

void CopyVariationAttributesNoLookups(tfx_variation_attributes_t *src, tfx_variation_attributes_t *dst) {
	if (src == dst) return;
	CopyGraphNoLookups(&src->life, &dst->life);
	CopyGraphNoLookups(&src->amount, &dst->amount);
	CopyGraphNoLookups(&src->velocity, &dst->velocity);
	CopyGraphNoLookups(&src->width, &dst->width);
	CopyGraphNoLookups(&src->height, &dst->height);
	CopyGraphNoLookups(&src->weight, &dst->weight);
	CopyGraphNoLookups(&src->spin, &dst->spin);
	CopyGraphNoLookups(&src->pitch_spin, &dst->pitch_spin);
	CopyGraphNoLookups(&src->yaw_spin, &dst->yaw_spin);
	CopyGraphNoLookups(&src->noise_offset, &dst->noise_offset);
	CopyGraphNoLookups(&src->noise_resolution, &dst->noise_resolution);
	CopyGraphNoLookups(&src->motion_randomness, &dst->motion_randomness);
}

void CopyVariationAttributes(tfx_variation_attributes_t *src, tfx_variation_attributes_t *dst) {
	CopyGraph(&src->life, &dst->life);
	CopyGraph(&src->amount, &dst->amount);
	CopyGraph(&src->velocity, &dst->velocity);
	CopyGraph(&src->width, &dst->width);
	CopyGraph(&src->height, &dst->height);
	CopyGraph(&src->weight, &dst->weight);
	CopyGraph(&src->spin, &dst->spin);
	CopyGraph(&src->pitch_spin, &dst->pitch_spin);
	CopyGraph(&src->yaw_spin, &dst->yaw_spin);
	CopyGraph(&src->noise_offset, &dst->noise_offset);
	CopyGraph(&src->noise_resolution, &dst->noise_resolution);
	CopyGraph(&src->motion_randomness, &dst->motion_randomness);
}

void FreeBaseAttributes(tfx_base_attributes_t *attributes) {
	FreeGraph(&attributes->life);
	FreeGraph(&attributes->amount);
	FreeGraph(&attributes->velocity);
	FreeGraph(&attributes->width);
	FreeGraph(&attributes->height);
	FreeGraph(&attributes->weight);
	FreeGraph(&attributes->spin);
	FreeGraph(&attributes->pitch_spin);
	FreeGraph(&attributes->yaw_spin);
	FreeGraph(&attributes->noise_offset);
}

void CopyBaseAttributesNoLookups(tfx_base_attributes_t *src, tfx_base_attributes_t *dst) {
	if (src == dst) return;
	CopyGraphNoLookups(&src->life, &dst->life);
	CopyGraphNoLookups(&src->amount, &dst->amount);
	CopyGraphNoLookups(&src->velocity, &dst->velocity);
	CopyGraphNoLookups(&src->width, &dst->width);
	CopyGraphNoLookups(&src->height, &dst->height);
	CopyGraphNoLookups(&src->weight, &dst->weight);
	CopyGraphNoLookups(&src->spin, &dst->spin);
	CopyGraphNoLookups(&src->pitch_spin, &dst->pitch_spin);
	CopyGraphNoLookups(&src->yaw_spin, &dst->yaw_spin);
	CopyGraphNoLookups(&src->noise_offset, &dst->noise_offset);
}

void CopyBaseAttributes(tfx_base_attributes_t *src, tfx_base_attributes_t *dst) {
	if (src == dst) return;
	CopyGraph(&src->life, &dst->life);
	CopyGraph(&src->amount, &dst->amount);
	CopyGraph(&src->velocity, &dst->velocity);
	CopyGraph(&src->width, &dst->width);
	CopyGraph(&src->height, &dst->height);
	CopyGraph(&src->weight, &dst->weight);
	CopyGraph(&src->spin, &dst->spin);
	CopyGraph(&src->pitch_spin, &dst->pitch_spin);
	CopyGraph(&src->yaw_spin, &dst->yaw_spin);
	CopyGraph(&src->noise_offset, &dst->noise_offset);
}

void InitialiseEmitterAttributes(tfx_emitter_attributes_t *attributes, tfxU32 bucket_size) {
	InitialisePropertyAttributes(&attributes->properties, bucket_size);
	InitialiseBaseAttributes(&attributes->base, bucket_size);
	InitialiseVariationAttributes(&attributes->variation, bucket_size);
	InitialiseOvertimeAttributes(&attributes->overtime, bucket_size);
	InitialiseFactorAttributes(&attributes->factor, bucket_size);
}

void FreeEmitterAttributes(tfx_emitter_attributes_t *attributes) {
	FreePropertyAttributes(&attributes->properties);
	FreeBaseAttributes(&attributes->base);
	FreeVariationAttributes(&attributes->variation);
	FreeOvertimeAttributes(&attributes->overtime);
	FreeFactorAttributes(&attributes->factor);
}

tfx_effect_emitter_t &tfx_library_t::operator[] (tfxU32 index) {
	return effects[index];
}

void MaybeGrowLibraryProperties(tfx_library_t *library, tfxU32 size_offset) {
	if (library->emitter_properties.current_size >= library->emitter_properties.capacity - size_offset) {
		library->emitter_properties.reserve(library->emitter_properties.capacity + library->emitter_properties.capacity / 2);
	}
}

void MaybeGrowLibraryInfos(tfx_library_t *library) {
	if (library->effect_infos.current_size >= library->effect_infos.capacity - 4) {
		library->effect_infos.reserve(library->effect_infos._grow_capacity(library->effect_infos.current_size + 1));
	}
}

tfx_effect_emitter_info_t *GetEffectInfo(tfx_effect_emitter_t *e) {
	TFX_ASSERT(e->library->effect_infos.size() > e->info_index);
	return &e->library->effect_infos[e->info_index];
}

bool RenameLibraryEffect(tfx_library_t *library, tfx_effect_emitter_t *effect, const char *new_name) {
	if (!LibraryNameExists(library, effect, new_name) && strlen(new_name) > 0) {
		SetEffectName(effect, new_name);
		UpdateLibraryEffectPaths(library);
		return true;
	}

	return false;
}

bool LibraryNameExists(tfx_library_t *library, tfx_effect_emitter_t *effect, const char *name) {
	for (auto &e : library->effects) {
		if (effect->library_index != e.library_index) {
			if (GetEffectInfo(&e)->name == name) {
				return true;
			}
		}
	}
	return false;
}

void UpdateLibraryEffectPaths(tfx_library_t *library) {
	library->effect_paths.Clear();
	for (auto &e : library->effects) {
		tfx_str256_t path = GetEffectInfo(&e)->name;
		GetEffectInfo(&e)->path = path;
		e.path_hash = tfxXXHash64::hash(path.c_str(), path.Length(), 0);
		AddLibraryPath(library, &e, &path);
	}
}

void BuildAllLibraryPaths(tfx_library_t *library) {
	for (tfxBucketLoop(library->paths)) {
		if (library->paths[i].flags & tfxPathFlags_2d) {
			BuildPathNodes2d(&library->paths[i]);
		}
		else {
			BuildPathNodes3d(&library->paths[i]);
		}
	}
}

void UpdateLibraryGPUImageData(tfx_library_t *library) {
	library->gpu_shapes.list.free();
    if(library->particle_shapes.Size() > 0) {
        BuildGPUShapeData(&library->particle_shapes.data, &library->gpu_shapes, library->uv_lookup);
    }
}

void AddLibraryPath(tfx_library_t *library, tfx_effect_emitter_t *effect_emitter, tfx_str256_t *path) {
	if (library->effect_paths.ValidName(path->c_str())) {
		*path = FindNewPathName(library, *path);
		GetEffectInfo(effect_emitter)->path = *path;
		GetEffectInfo(effect_emitter)->name = GetNameFromPath(path);
		effect_emitter->path_hash = tfxXXHash64::hash(path->c_str(), path->Length(), 0);
	}
	library->effect_paths.Insert(*path, effect_emitter);
	for (auto &sub : GetEffectInfo(effect_emitter)->sub_effectors) {
		tfx_str256_t sub_path = *path;
		sub_path.Appendf("/%s", GetEffectInfo(&sub)->name.c_str());
		GetEffectInfo(&sub)->path = sub_path;
		sub.path_hash = tfxXXHash64::hash(sub_path.c_str(), sub_path.Length(), 0);
		AddLibraryPath(library, &sub, &sub_path);
	}
}

tfx_effect_emitter_t *InsertLibraryEffect(tfx_library_t *library, tfx_effect_emitter_t *effect, tfx_effect_emitter_t *position) {
	effect->library_index = library->effects.current_size;
	effect->type = tfxEffectType;
	GetEffectInfo(effect)->uid = ++library->uid;
	effect->library = library;
	tfx_effect_emitter_t *inserted_effect = library->effects.insert_after(position, *effect);
	ReIndexLibrary(library);
	UpdateLibraryEffectPaths(library);
	return inserted_effect;
}

tfx_effect_emitter_t *AddLibraryEffect(tfx_library_t *library, tfx_effect_emitter_t *effect) {
	effect->library_index = library->effects.current_size;
	effect->type = tfxEffectType;
	GetEffectInfo(effect)->uid = ++library->uid;
	effect->library = library;
	library->effects.push_back(*effect);
	ReIndexLibrary(library);
	UpdateLibraryEffectPaths(library);
	return &library->effects.back();
}

tfx_effect_emitter_t *AddLibraryFolder(tfx_library_t *library, tfx_str64_t *name) {
	tfx_effect_emitter_t folder;
	folder.info_index = AddLibraryEffectEmitterInfo(library);
	folder.library = library;
	GetEffectInfo(&folder)->name = *name;
	folder.type = tfxFolder;
	folder.library = library;
	GetEffectInfo(&folder)->uid = ++library->uid;
	library->effects.push_back(folder);
	ReIndexLibrary(library);
	UpdateLibraryEffectPaths(library);
	return &library->effects.back();
}

tfx_effect_emitter_t *AddLibraryFolder(tfx_library_t *library, tfx_effect_emitter_t *folder) {
	TFX_ASSERT(folder->type == tfxFolder);            //Must be type tfxFolder if adding a folder
	folder->library = library;
	GetEffectInfo(folder)->uid = ++library->uid;
	library->effects.push_back(*folder);
	ReIndexLibrary(library);
	UpdateLibraryEffectPaths(library);
	return &library->effects.back();
}

tfx_effect_emitter_t *AddLibraryStage(tfx_library_t *library, tfx_str64_t *name) {
	tfx_effect_emitter_t stage;
	stage.info_index = AddLibraryEffectEmitterInfo(library);
	stage.library = library;
	GetEffectInfo(&stage)->name = *name;
	stage.type = tfxStage;
	GetEffectInfo(&stage)->uid = ++library->uid;
	library->effects.push_back(stage);
	ReIndexLibrary(library);
	UpdateLibraryEffectPaths(library);
	return &library->effects.back();
}

tfx_effect_emitter_t *GetLibraryEffect(tfx_library_t *library, const char *path) {
	TFX_ASSERT(library->effect_paths.ValidName(path));        //Effect was not found by that name
	return library->effect_paths.At(path);
}

bool IsValidEffectPath(tfx_library_t *library, const char *path) {
	return library->effect_paths.ValidName(path);
}

bool IsValidEffectKey(tfx_library_t *library, tfxKey key) {
	return library->effect_paths.ValidKey(key);
}

tfx_effect_emitter_t *GetLibraryEffect(tfx_library_t *library, tfxKey key) {
	TFX_ASSERT(library->effect_paths.ValidKey(key));            //Effect was not found by that key
	return library->effect_paths.At(key);
}

void PrepareLibraryEffectTemplate(tfx_library_t *library, tfx_str256_t path, tfx_effect_template_t *effect_template) {
	tfx_effect_emitter_t *effect = GetLibraryEffect(library, path.c_str());
	TFX_ASSERT(effect);                                //Effect was not found, make sure the path exists
	TFX_ASSERT(effect->type == tfxEffectType);        //The effect must be an effect type, not an emitter
	effect_template->original_effect_hash = effect->path_hash;
	CloneEffect(effect, &effect_template->effect, &effect_template->effect, library, tfxEffectCloningFlags_clone_graphs | tfxEffectCloningFlags_compile_graphs);
	AddTemplatePath(effect_template, &effect_template->effect, GetEffectInfo(&effect_template->effect)->name.c_str());
}

void PrepareLibraryEffectTemplate(tfx_library_t *library, tfx_effect_emitter_t *effect, tfx_effect_template_t *effect_template) {
	TFX_ASSERT(effect->type == tfxEffectType);
	CloneEffect(effect, &effect_template->effect, &effect_template->effect, library);
	AddTemplatePath(effect_template, &effect_template->effect, GetEffectInfo(&effect_template->effect)->name.c_str());
}

void ReIndexLibrary(tfx_library_t *library) {
	tfxU32 index = 0;
	for (auto &e : library->effects) {
		e.library_index = index++;
		e.parent = nullptr;
		ReIndexEffect(&e);
	}
}

void UpdateLibraryParticleShapeReferences(tfx_library_t *library, tfxKey default_hash) {
	tfx_vector_t<tfx_effect_emitter_t *> stack;
	for (auto &effect : library->effects) {
		stack.push_back(&effect);
	}
	while (stack.size()) {
		tfx_effect_emitter_t &current = *stack.pop_back();
		if (current.type == tfxEmitterType) {
			bool shape_found = false;
			tfxKey hash = library->emitter_properties[current.property_index].image_hash;
			if (hash == 0) {
				//Try to match index instead might be a converted eff file
				tfxU32 image_index = library->emitter_properties[current.property_index].image_index;
				for (tfx_image_data_t &image_data : library->particle_shapes.data) {
					if (image_data.shape_index == image_index) {
						library->emitter_properties[current.property_index].image = &image_data;
						library->emitter_properties[current.property_index].end_frame = image_data.animation_frames - 1;
						library->emitter_properties[current.property_index].image_hash = image_data.image_hash;
					}
				}
			}
			if (library->particle_shapes.ValidKey(library->emitter_properties[current.property_index].image_hash)) {
				library->emitter_properties[current.property_index].image = &library->particle_shapes.At(library->emitter_properties[current.property_index].image_hash);
				library->emitter_properties[current.property_index].end_frame = library->particle_shapes.At(library->emitter_properties[current.property_index].image_hash).animation_frames - 1;
				shape_found = true;
			}
			else {
				for (auto &shape : library->particle_shapes.data) {
					if (shape.image_hash == library->emitter_properties[current.property_index].image_hash) {
						library->emitter_properties[current.property_index].image_hash = shape.image_hash;
						library->emitter_properties[current.property_index].image = &library->particle_shapes.At(library->emitter_properties[current.property_index].image_hash);
						library->emitter_properties[current.property_index].end_frame = library->particle_shapes.At(library->emitter_properties[current.property_index].image_hash).animation_frames - 1;
						shape_found = true;
						break;
					}
				}
			}
			if (!shape_found) {
				library->emitter_properties[current.property_index].image = &library->particle_shapes.At(default_hash);
				library->emitter_properties[current.property_index].end_frame = library->particle_shapes.At(default_hash).animation_frames - 1;
			}
		}
		for (auto &sub : GetEffectInfo(&current)->sub_effectors) {
			stack.push_back(&sub);
		}
	}
}

tfx_effect_emitter_t *LibraryMoveUp(tfx_library_t *library, tfx_effect_emitter_t *effect) {
	if (effect->library_index > 0) {
		tfxU32 new_index = effect->library_index - 1;
		std::swap(library->effects[effect->library_index], library->effects[new_index]);
		UpdateLibraryEffectPaths(library);
		ReIndexLibrary(library);
		return &library->effects[new_index];
	}
	return nullptr;
}

tfx_effect_emitter_t *LibraryMoveDown(tfx_library_t *library, tfx_effect_emitter_t *effect) {
	if (effect->library_index < library->effects.size() - 1) {
		tfxU32 new_index = effect->library_index + 1;
		std::swap(library->effects[effect->library_index], library->effects[new_index]);
		UpdateLibraryEffectPaths(library);
		ReIndexLibrary(library);
		return &library->effects[new_index];
	}
	return nullptr;
}

void BuildGPUShapeData(tfx_vector_t<tfx_image_data_t> *particle_shapes, tfx_gpu_shapes_t *shape_data, void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset)) {
	TFX_ASSERT(particle_shapes->size());        //There are no shapes to copy!
    TFX_ASSERT(uv_lookup);  //You must set a function that applies the uv coordinates for each image you load
	tfxU32 index = 0;
	shape_data->list.clear();
	for (tfx_image_data_t &shape : *particle_shapes) {
		if (shape.animation_frames == 1) {
			tfx_gpu_image_data_t cs;
			cs.animation_frames = shape.animation_frames;
			cs.image_size = shape.image_size;
			uv_lookup(shape.ptr, &cs, 0);
			shape_data->list.push_back_copy(cs);
			shape.compute_shape_index = index++;
		}
		else {
			shape.compute_shape_index = index;
			for (int f = 0; f != shape.animation_frames; ++f) {
				tfx_gpu_image_data_t cs;
				cs.animation_frames = shape.animation_frames;
				cs.image_size = shape.image_size;
				uv_lookup(shape.ptr, &cs, f);
				shape_data->list.push_back_copy(cs);
				index++;
			}
		}
	}
}

void CopyLibraryLookupIndexesData(tfx_library_t *library, void *dst) {
	TFX_ASSERT(dst);    //must be a valid pointer to a space in memory
	TFX_ASSERT(library->compiled_lookup_indexes.size());        //There is no data to copy, make sure a library has been loaded properly and it contains effects with emitters
	memcpy(dst, library->compiled_lookup_indexes.data, GetLibraryLookupIndexesSizeInBytes(library));
}

void CopyLibraryLookupValuesData(tfx_library_t *library, void *dst) {
	TFX_ASSERT(dst);    //must be a valid pointer to a space in memory
	TFX_ASSERT(library->compiled_lookup_indexes.size());        //There is no data to copy, make sure a library has been loaded properly and it contains effects with emitters
	memcpy(dst, library->compiled_lookup_values.data, GetLibraryLookupValuesSizeInBytes(library));
}

tfxU32 GetLibraryComputeShapeDataSizeInBytes(tfx_library_t *library) {
	tfxU32 frame_count = 0;
	for (auto &shape : library->particle_shapes.data) {
		frame_count += (tfxU32)shape.animation_frames;
	}
	return frame_count * sizeof(tfx_gpu_image_data_t);
}

tfxU32 GetLibraryComputeShapeCount(tfx_library_t *library) {
	tfxU32 frame_count = 0;
	for (auto &shape : library->particle_shapes.data) {
		frame_count += (tfxU32)shape.animation_frames;
	}
	return frame_count;
}

tfxU32 GetLibraryLookupIndexCount(tfx_library_t *library) {
	return library->compiled_lookup_indexes.size() * TFX_OVERTIME_COUNT;
}

tfxU32 GetLibraryLookupValueCount(tfx_library_t *library) {
	return library->compiled_lookup_values.size();
}

tfxU32 GetLibraryLookupIndexesSizeInBytes(tfx_library_t *library) {
	return sizeof(tfx_graph_lookup_index_t) * TFX_OVERTIME_COUNT * library->compiled_lookup_indexes.size();
}

tfxU32 GetLibraryLookupValuesSizeInBytes(tfx_library_t *library) {
	return sizeof(float) * library->compiled_lookup_values.size();
}

bool IsLibraryShapeUsed(tfx_library_t *library, tfxKey image_hash) {
	tmpStack(tfx_effect_emitter_t *, effect_stack);
	for (auto &effect : library->effects) {
		effect_stack.push_back(&effect);
	}
	while (effect_stack.size()) {
		auto current = effect_stack.pop_back();
		if (current->type == tfxEmitterType) {
			if (GetEffectProperties(current)->image->image_hash == image_hash) {
				return true;
			}
		}
		for (auto &sub : GetEffectInfo(current)->sub_effectors) {
			effect_stack.push_back(&sub);
		}
	}
	return false;
}

bool LibraryShapeExists(tfx_library_t *library, tfxKey image_hash) {
	return library->particle_shapes.ValidKey(image_hash);
}

bool RemoveLibraryShape(tfx_library_t *library, tfxKey image_hash) {
	if (!library->particle_shapes.ValidKey(image_hash)) {
		return false;
	}
	library->particle_shapes.Remove(image_hash);
	for (auto &m : library->particle_shapes.map) {
		library->particle_shapes[m.index].image_hash = m.key;
	}
	return true;
}

void DeleteLibraryEffect(tfx_library_t *library, tfx_effect_emitter_t *effect) {
	CleanUpEffect(&library->effects[effect->library_index]);
	library->effects.erase(&library->effects[effect->library_index]);

	UpdateLibraryEffectPaths(library);
	ReIndexLibrary(library);
}

tfxU32 AddLibraryGlobal(tfx_library_t *library) {
	if (library->free_global_graphs.size()) {
		return library->free_global_graphs.pop_back();
	}
	tfx_global_attributes_t global;
	InitialiseGlobalAttributes(&global);
	library->global_graphs.push_back(global);
	return library->global_graphs.size() - 1;
}

tfxU32 AddLibraryKeyframes(tfx_library_t *library) {
	if (library->free_keyframes.size())
		return library->free_keyframes.pop_back();
	tfx_transform_attributes_t keyframes;
	InitialiseTransformAttributes(&keyframes);
	library->transform_attributes.push_back(keyframes);
	return library->transform_attributes.size() - 1;
}

tfxU32 AddLibraryEmitterAttributes(tfx_library_t *library) {
	if (library->free_emitter_attributes.size()) {
		return library->free_emitter_attributes.pop_back();
	}
	tfx_emitter_attributes_t attributes;
	InitialiseEmitterAttributes(&attributes);
	library->emitter_attributes.push_back(attributes);
	return library->emitter_attributes.size() - 1;
}

void FreeLibraryGlobal(tfx_library_t *library, tfxU32 index) {
	TFX_ASSERT(index < library->global_graphs.size());
	library->free_global_graphs.push_back(index);
	FreeGlobalAttributes(&library->global_graphs[index]);
}

void FreeLibraryKeyframes(tfx_library_t *library, tfxU32 index) {
	TFX_ASSERT(index < library->transform_attributes.size());
	library->free_keyframes.push_back(index);
	FreeTransformAttributes(&library->transform_attributes[index]);
}

void FreeLibraryEmitterAttributes(tfx_library_t *library, tfxU32 index) {
	TFX_ASSERT(index < library->emitter_attributes.size());
	library->free_emitter_attributes.push_back(index);
	FreeEmitterAttributes(&library->emitter_attributes[index]);
}

void FreeLibraryProperties(tfx_library_t *library, tfxU32 index) {
	TFX_ASSERT(index < library->emitter_properties.current_size);
	library->free_properties.push_back(index);
}

void FreeLibraryInfo(tfx_library_t *library, tfxU32 index) {
	TFX_ASSERT(index < library->effect_infos.size());
	library->free_infos.push_back(index);
}

tfxU32 CountLibraryKeyframeLookUpValues(tfx_library_t *library, tfxU32 index) {
	auto &transform = library->transform_attributes[index];
	tfxU32 count = 0;
	count += transform.roll.lookup.values.capacity;
	count += transform.pitch.lookup.values.capacity;
	count += transform.yaw.lookup.values.capacity;
	count += transform.translation_x.lookup.values.capacity;
	count += transform.translation_y.lookup.values.capacity;
	count += transform.translation_z.lookup.values.capacity;
	return count;
}

tfxU32 CountLibraryGlobalLookUpValues(tfx_library_t *library, tfxU32 index) {
	auto &global = library->global_graphs[index];
	tfxU32 count = 0;
	count += global.life.lookup.values.capacity;
	count += global.amount.lookup.values.capacity;
	count += global.velocity.lookup.values.capacity;
	count += global.width.lookup.values.capacity;
	count += global.height.lookup.values.capacity;
	count += global.weight.lookup.values.capacity;
	count += global.spin.lookup.values.capacity;
	count += global.pitch_spin.lookup.values.capacity;
	count += global.yaw_spin.lookup.values.capacity;
	count += global.stretch.lookup.values.capacity;
	count += global.overal_scale.lookup.values.capacity;
	count += global.intensity.lookup.values.capacity;
	count += global.splatter.lookup.values.capacity;
	count += global.emitter_width.lookup.values.capacity;
	count += global.emitter_height.lookup.values.capacity;
	count += global.emitter_depth.lookup.values.capacity;
	return count;
}

tfxU32 CountLibraryEmitterLookUpValues(tfx_library_t *library, tfxU32 index) {
	auto &attributes = library->emitter_attributes[index];
	tfxU32 count = 0;

	count += attributes.properties.emission_pitch.lookup.values.capacity;
	count += attributes.properties.emission_yaw.lookup.values.capacity;
	count += attributes.properties.emission_range.lookup.values.capacity;
	count += attributes.properties.splatter.lookup.values.capacity;
	count += attributes.properties.emitter_width.lookup.values.capacity;
	count += attributes.properties.emitter_height.lookup.values.capacity;
	count += attributes.properties.emitter_depth.lookup.values.capacity;
	count += attributes.properties.extrusion.lookup.values.capacity;
	count += attributes.properties.arc_size.lookup.values.capacity;
	count += attributes.properties.arc_offset.lookup.values.capacity;

	count += attributes.base.life.lookup.values.capacity;
	count += attributes.base.amount.lookup.values.capacity;
	count += attributes.base.velocity.lookup.values.capacity;
	count += attributes.base.width.lookup.values.capacity;
	count += attributes.base.height.lookup.values.capacity;
	count += attributes.base.weight.lookup.values.capacity;
	count += attributes.base.spin.lookup.values.capacity;
	count += attributes.base.pitch_spin.lookup.values.capacity;
	count += attributes.base.yaw_spin.lookup.values.capacity;
	count += attributes.base.noise_offset.lookup.values.capacity;

	count += attributes.variation.life.lookup.values.capacity;
	count += attributes.variation.amount.lookup.values.capacity;
	count += attributes.variation.velocity.lookup.values.capacity;
	count += attributes.variation.width.lookup.values.capacity;
	count += attributes.variation.height.lookup.values.capacity;
	count += attributes.variation.weight.lookup.values.capacity;
	count += attributes.variation.spin.lookup.values.capacity;
	count += attributes.variation.pitch_spin.lookup.values.capacity;
	count += attributes.variation.yaw_spin.lookup.values.capacity;
	count += attributes.variation.noise_offset.lookup.values.capacity;
	count += attributes.variation.noise_resolution.lookup.values.capacity;
	count += attributes.variation.motion_randomness.lookup.values.capacity;

	count += attributes.overtime.velocity.lookup.values.capacity;
	count += attributes.overtime.width.lookup.values.capacity;
	count += attributes.overtime.height.lookup.values.capacity;
	count += attributes.overtime.weight.lookup.values.capacity;
	count += attributes.overtime.spin.lookup.values.capacity;
	count += attributes.overtime.pitch_spin.lookup.values.capacity;
	count += attributes.overtime.yaw_spin.lookup.values.capacity;
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
	count += attributes.overtime.motion_randomness.lookup.values.capacity;

	count += attributes.factor.life.lookup.values.capacity;
	count += attributes.factor.velocity.lookup.values.capacity;
	count += attributes.factor.size.lookup.values.capacity;
	count += attributes.factor.intensity.lookup.values.capacity;

	return count;
}

tfxU32 CloneLibraryGlobal(tfx_library_t *library, tfxU32 source_index, tfx_library_t *destination_library) {
	tfxU32 index = AddLibraryGlobal(destination_library);
	CopyGlobalAttributesNoLookups(&library->global_graphs[source_index], &destination_library->global_graphs[index]);
	return index;
}

tfxU32 CloneLibraryKeyframes(tfx_library_t *library, tfxU32 source_index, tfx_library_t *destination_library) {
	tfxU32 index = AddLibraryKeyframes(destination_library);
	CopyTransformAttributesNoLookups(&library->transform_attributes[source_index], &destination_library->transform_attributes[index]);
	return index;
}

tfxU32 CloneLibraryEmitterAttributes(tfx_library_t *library, tfxU32 source_index, tfx_library_t *destination_library) {
	tfxU32 index = AddLibraryEmitterAttributes(destination_library);
	CopyPropertyAttributesNoLookups(&library->emitter_attributes[source_index].properties, &destination_library->emitter_attributes[index].properties);
	CopyBaseAttributesNoLookups(&library->emitter_attributes[source_index].base, &destination_library->emitter_attributes[index].base);
	CopyVariationAttributesNoLookups(&library->emitter_attributes[source_index].variation, &destination_library->emitter_attributes[index].variation);
	CopyOvertimeAttributesNoLookups(&library->emitter_attributes[source_index].overtime, &destination_library->emitter_attributes[index].overtime);
	CopyFactorAttributesNoLookups(&library->emitter_attributes[source_index].factor, &destination_library->emitter_attributes[index].factor);
	return index;
}

tfxU32 CloneLibraryInfo(tfx_library_t *library, tfxU32 source_index, tfx_library_t *destination_library) {
	tfxU32 index = AddLibraryEffectEmitterInfo(destination_library);
	destination_library->effect_infos[index] = library->effect_infos[source_index];
	return index;
}

tfxU32 CloneLibraryProperties(tfx_library_t *library, tfx_emitter_properties_t *source, tfx_library_t *destination_library) {
	tfxU32 dst_index = AddLibraryEmitterProperties(destination_library);
	CopyEmitterProperites(source, &destination_library->emitter_properties[dst_index]);
	return dst_index;
}

void AddLibraryEmitterGraphs(tfx_library_t *library, tfx_effect_emitter_t *emitter) {
	emitter->emitter_attributes = AddLibraryEmitterAttributes(library);
}

void AddLibraryTransformGraphs(tfx_library_t *library, tfx_effect_emitter_t *emitter) {
	emitter->transform_attributes = AddLibraryKeyframes(library);
}

void AddLibraryEffectGraphs(tfx_library_t *library, tfx_effect_emitter_t *effect) {
	tfx_effect_emitter_t *root_effect = tfx_GetRootEffect(effect);
	if (root_effect == effect)
		effect->global = AddLibraryGlobal(library);
	else
		effect->global = root_effect->global;
}

tfxU32 AddLibrarySpriteSheetSettings(tfx_library_t *library, tfx_effect_emitter_t *effect) {
	TFX_ASSERT(effect->type == tfxEffectType);
	tfx_sprite_sheet_settings_t a{};
	a.frames = 32;
	a.current_frame = 1;
	a.frame_offset = 0;
	a.extra_frames_count = 0;
	a.position = tfx_vec2_t(0.f, 0.f);
	a.frame_size = tfx_vec2_t(256.f, 256.f);
	a.playback_speed = 1.f;
	a.animation_flags = tfxAnimationFlags_needs_recording | tfxAnimationFlags_export_with_transparency;
	a.seed = 0;
	a.zoom = 1.f;
	a.scale = 1.f;
	a.needs_exporting = 0;
	a.color_option = tfx_export_color_options::tfxFullColor;
	a.export_option = tfx_export_options::tfxSpriteSheet;
	a.camera_settings.camera_floor_height = -10.f;
	a.camera_settings.camera_fov = tfx__degrees_to_radians(60);
	a.camera_settings.camera_pitch = tfx__degrees_to_radians(-30.f);
	a.camera_settings.camera_yaw = tfx__degrees_to_radians(-90.f);
	a.camera_settings.camera_position = tfx_vec3_t(0.f, 3.5f, 7.5f);
	a.camera_settings.camera_isometric = false;
	a.camera_settings.camera_isometric_scale = 5.f;
	a.camera_settings.camera_hide_floor = false;
	a.camera_settings_orthographic.camera_floor_height = -10.f;
	a.camera_settings_orthographic.camera_fov = tfx__degrees_to_radians(60);
	a.camera_settings_orthographic.camera_pitch = tfx__degrees_to_radians(-30.f);
	a.camera_settings_orthographic.camera_yaw = tfx__degrees_to_radians(-90.f);
	a.camera_settings_orthographic.camera_position = tfx_vec3_t(0.f, 3.5f, 7.5f);
	a.camera_settings_orthographic.camera_isometric = true;
	a.camera_settings_orthographic.camera_isometric_scale = 5.f;
	a.camera_settings_orthographic.camera_hide_floor = false;
	library->sprite_sheet_settings.push_back(a);
	GetEffectInfo(effect)->sprite_sheet_settings_index = library->sprite_sheet_settings.size() - 1;
	return GetEffectInfo(effect)->sprite_sheet_settings_index;
}

void AddLibrarySpriteSheetSettingsSub(tfx_library_t *library, tfx_effect_emitter_t *effect) {
	if (effect->type == tfxEffectType) {
		tfx_sprite_sheet_settings_t a{};
		a.frames = 32;
		a.current_frame = 1;
		a.frame_offset = 0;
		a.extra_frames_count = 0;
		a.position = tfx_vec2_t(0.f, 0.f);
		a.frame_size = tfx_vec2_t(256.f, 256.f);
		a.playback_speed = 1.f;
		a.animation_flags = tfxAnimationFlags_needs_recording | tfxAnimationFlags_export_with_transparency;
		a.seed = 0;
		a.zoom = 1.f;
		a.scale = 1.f;
		a.needs_exporting = 0;
		a.color_option = tfx_export_color_options::tfxFullColor;
		a.export_option = tfx_export_options::tfxSpriteSheet;
		a.camera_settings.camera_floor_height = -10.f;
		a.camera_settings.camera_fov = tfx__degrees_to_radians(60);
		a.camera_settings.camera_pitch = tfx__degrees_to_radians(-30.f);
		a.camera_settings.camera_yaw = tfx__degrees_to_radians(-90.f);
		a.camera_settings.camera_position = tfx_vec3_t(0.f, 3.5f, 7.5f);
		a.camera_settings.camera_isometric = false;
		a.camera_settings.camera_isometric_scale = 5.f;
		a.camera_settings.camera_hide_floor = false;
		a.camera_settings_orthographic.camera_floor_height = -10.f;
		a.camera_settings_orthographic.camera_fov = tfx__degrees_to_radians(60);
		a.camera_settings_orthographic.camera_pitch = tfx__degrees_to_radians(-30.f);
		a.camera_settings_orthographic.camera_yaw = tfx__degrees_to_radians(-90.f);
		a.camera_settings_orthographic.camera_position = tfx_vec3_t(0.f, 3.5f, 7.5f);
		a.camera_settings_orthographic.camera_isometric = true;
		a.camera_settings_orthographic.camera_isometric_scale = 5.f;
		a.camera_settings_orthographic.camera_hide_floor = false;
		library->sprite_sheet_settings.push_back(a);
		GetEffectInfo(effect)->sprite_sheet_settings_index = library->sprite_sheet_settings.size() - 1;
		for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
			AddLibrarySpriteSheetSettingsSub(effect->library, &sub);
		}
	}
	else {
		for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
			AddLibrarySpriteSheetSettingsSub(effect->library, &sub);
		}
	}
}

tfxU32 AddLibrarySpriteDataSettings(tfx_library_t *library, tfx_effect_emitter_t *effect) {
	TFX_ASSERT(effect->type == tfxEffectType);
	tfx_sprite_data_settings_t a{};
	a.real_frames = 32;
	a.frames_after_compression = 32;
	a.current_frame = 1;
	a.frame_offset = 0;
	a.extra_frames_count = 0;
	a.playback_speed = 1.f;
	a.animation_flags = tfxAnimationFlags_needs_recording | tfxAnimationFlags_export_with_transparency;
	a.seed = 0;
	a.needs_exporting = 0;
	a.recording_frame_rate = 60.f;
	//a.camera_settings.camera_floor_height = -10.f;
	//a.camera_settings.camera_fov = tfx__degrees_to_radians(60);
	//a.camera_settings.camera_pitch = tfx__degrees_to_radians(-30.f);
	//a.camera_settings.camera_yaw = tfx__degrees_to_radians(-90.f);
	//a.camera_settings.camera_position = tfx_vec3_t(0.f, 3.5f, 7.5f);
	//a.camera_settings.camera_isometric = false;
	//a.camera_settings.camera_isometric_scale = 5.f;
	//a.camera_settings.camera_hide_floor = false;
	library->sprite_data_settings.push_back(a);
	GetEffectInfo(effect)->sprite_data_settings_index = library->sprite_data_settings.size() - 1;
	return GetEffectInfo(effect)->sprite_data_settings_index;
}

void AddLibrarySpriteDataSettingsSub(tfx_library_t *library, tfx_effect_emitter_t *effect) {
	if (effect->type == tfxEffectType) {
		tfx_sprite_data_settings_t a{};
		a.real_frames = 32;
		a.frames_after_compression = 32;
		a.current_frame = 1;
		a.frame_offset = 0;
		a.extra_frames_count = 0;
		a.playback_speed = 1.f;
		a.recording_frame_rate = 60.f;
		a.animation_flags = tfxAnimationFlags_needs_recording | tfxAnimationFlags_export_with_transparency;
		a.seed = 0;
		a.needs_exporting = 0;
		library->sprite_data_settings.push_back(a);
		GetEffectInfo(effect)->sprite_data_settings_index = library->sprite_data_settings.size() - 1;
		for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
			AddLibrarySpriteDataSettingsSub(effect->library, &sub);
		}
	}
	else {
		for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
			AddLibrarySpriteDataSettingsSub(effect->library, &sub);
		}
	}
}

tfxU32 AddLibraryPreviewCameraSettings(tfx_library_t *library, tfx_effect_emitter_t *effect) {
	TFX_ASSERT(effect->type == tfxEffectType || effect->type == tfxStage);
	tfx_preview_camera_settings_t a{};
	a.camera_settings.camera_floor_height = -10.f;
	a.camera_settings.camera_fov = tfx__degrees_to_radians(60);
	a.camera_settings.camera_pitch = tfx__degrees_to_radians(-30.f);
	a.camera_settings.camera_yaw = tfx__degrees_to_radians(-90.f);
	a.camera_settings.camera_position = tfx_vec3_t(0.f, 3.5f, 7.5f);
	a.camera_settings.camera_isometric = false;
	a.camera_settings.camera_isometric_scale = 5.f;
	a.camera_settings.camera_hide_floor = false;
	a.effect_z_offset = 5.f;
	a.camera_speed = 6.f;
	a.attach_effect_to_camera = false;
	library->preview_camera_settings.push_back(a);
	GetEffectInfo(effect)->preview_camera_settings = library->preview_camera_settings.size() - 1;
	return GetEffectInfo(effect)->preview_camera_settings;
}

void AddLibraryPreviewCameraSettingsSub(tfx_library_t *library, tfx_effect_emitter_t *effect) {
	if (effect->type == tfxEffectType) {
		tfx_preview_camera_settings_t a{};
		a.camera_settings.camera_floor_height = -10.f;
		a.camera_settings.camera_fov = tfx__degrees_to_radians(60);
		a.camera_settings.camera_pitch = tfx__degrees_to_radians(-30.f);
		a.camera_settings.camera_yaw = tfx__degrees_to_radians(-90.f);
		a.camera_settings.camera_position = tfx_vec3_t(0.f, 3.5f, 7.5f);
		a.camera_settings.camera_isometric = false;
		a.camera_settings.camera_isometric_scale = 5.f;
		a.camera_settings.camera_hide_floor = false;
		a.effect_z_offset = 5.f;
		a.camera_speed = 6.f;
		a.attach_effect_to_camera = false;
		library->preview_camera_settings.push_back(a);
		GetEffectInfo(effect)->preview_camera_settings = library->preview_camera_settings.size() - 1;
		for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
			AddLibraryPreviewCameraSettingsSub(effect->library, &sub);
		}
	}
	else {
		for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
			AddLibraryPreviewCameraSettingsSub(effect->library, &sub);
		}
	}
}

tfxU32 AddLibraryPreviewCameraSettings(tfx_library_t *library) {
	tfx_preview_camera_settings_t a{};
	a.camera_settings.camera_floor_height = -10.f;
	a.camera_settings.camera_fov = tfx__degrees_to_radians(60);
	a.camera_settings.camera_pitch = tfx__degrees_to_radians(-30.f);
	a.camera_settings.camera_yaw = tfx__degrees_to_radians(-90.f);
	a.camera_settings.camera_position = tfx_vec3_t(0.f, 3.5f, 7.5f);
	a.camera_settings.camera_isometric = false;
	a.camera_settings.camera_isometric_scale = 5.f;
	a.camera_settings.camera_hide_floor = false;
	a.effect_z_offset = 5.f;
	a.camera_speed = 6.f;
	a.attach_effect_to_camera = false;
	library->preview_camera_settings.push_back(a);
	return library->preview_camera_settings.size() - 1;
}

tfxU32 AddLibraryEffectEmitterInfo(tfx_library_t *library) {
	tfx_effect_emitter_info_t info{};
	if (library->free_infos.size()) {
		return library->free_infos.pop_back();
	}
	library->effect_infos.push_back(info);
	return library->effect_infos.size() - 1;
}

tfxU32 AddLibraryEmitterProperties(tfx_library_t *library) {
	if (library->free_properties.size()) {
		return library->free_properties.pop_back();
	}
	tfx_emitter_properties_t properties{};
	library->emitter_properties.push_back(properties);
	return library->emitter_properties.current_size - 1;
}

void InitLibrary(tfx_library_t *library) {
    tfx_color_ramp_t ramp{};
    for(int x = 0; x != tfxCOLOR_RAMP_WIDTH; ++x) {
        ramp.colors[x].color = 0xFFFFFFFF;
    }
    AddColorRampToBitmaps(&library->color_ramps, &ramp);
    library->gpu_shapes.list.set_alignment(16);
}

tfx_str64_t GetNameFromPath(tfx_str256_t *path) {
	tfx_str16_t extension;
	tmpStack(tfx_str256_t, file_split);
	tfx__split_string_stack(*path, &file_split, '/');
	if (file_split.size()) {
		extension = file_split.back();
	}
	return extension;
}

tfx_str256_t FindNewPathName(tfx_library_t *library, const tfx_str256_t &path) {
	tmpStack(tfx_str256_t, name);
	tfx__split_string_stack(path, &name, 46);
	tfx_str256_t new_path;
	if (name.size() > 2) {
		for (int i = 0; i != name.size() - 2; ++i) {
			new_path.Appendf(name[i].c_str());
		}
	}
	else if (name.size() == 2) {
		new_path.Appendf(name[0].c_str());
	}
	else {
		new_path = path;
	}
	tfx_str256_t find_name = new_path;
	int index = 1;
	while (library->effect_paths.ValidName(find_name.c_str())) {
		find_name = new_path;
		find_name.Appendf(".%i", index++);
	}
	return find_name;
}

void ClearLibrary(tfx_library_t *library) {
	for (auto &e : library->effects) {
		FreeEffectGraphs(&e);
	}
	library->effects.free_all();
	library->effect_paths.FreeAll();
	library->particle_shapes.FreeAll();
	for (tfx_global_attributes_t &global : library->global_graphs) {
		FreeGlobalAttributes(&global);
	}
	for (tfx_emitter_attributes_t &emitter_attributes : library->emitter_attributes) {
		FreeEmitterAttributes(&emitter_attributes);
	}
	for (tfx_transform_attributes_t &transform_attributes : library->transform_attributes) {
		FreeTransformAttributes(&transform_attributes);
	}
	for (tfx_bitmap_t &bitmap : library->color_ramps.color_ramp_bitmaps) {
		tfxFreeBitmap(&bitmap);
	}
	library->color_ramps.color_ramp_ids.FreeAll();
	library->color_ramps.color_ramp_bitmaps.free();
	library->global_graphs.free();
	library->emitter_attributes.free();
	library->transform_attributes.free();
	library->paths.free_all();
	library->sprite_sheet_settings.free_all();
	library->preview_camera_settings.free_all();
	library->all_nodes.free_all();
	library->node_lookup_indexes.free_all();
	library->compiled_lookup_indexes.free_all();
	library->compiled_lookup_values.free_all();
	library->emitter_properties.free();
	library->graph_min_max.free_all();
	for (auto &info : library->effect_infos) {
		info.sub_effectors.free_all();
	}
	library->effect_infos.free_all();
	AddLibraryPreviewCameraSettings(library);
	library->pre_recorded_effects.FreeAll();

	library->free_global_graphs.free_all();
	library->free_keyframe_graphs.free_all();
	library->free_emitter_attributes.free_all();
	library->free_animation_settings.free_all();
	library->free_preview_camera_settings.free_all();
	library->free_infos.free_all();
	library->free_properties.free_all();
	library->free_keyframes.free_all();

	library->uid = 0;
	library->color_ramps.color_ramp_count = 0;
}

void UpdateLibraryComputeNodes(tfx_library_t *library) {
	tfxU32 running_node_index = 0;
	tfxU32 running_value_index = 0;
	tmpStack(tfx_effect_emitter_t *, stack);
	library->all_nodes.clear();
	library->node_lookup_indexes.clear();
	library->compiled_lookup_values.clear();
	library->compiled_lookup_indexes.clear();
	for (auto &effect : library->effects) {
		stack.push_back(&effect);
		while (!stack.empty()) {
			tfx_effect_emitter_t *current = stack.pop_back();
			if (current->type == tfxFolder) {
				for (auto &sub : GetEffectInfo(current)->sub_effectors) {
					stack.push_back(&sub);
				}
				continue;
			}
			tfx_effect_lookup_data_t lookup_data;
			tfx_effect_lookup_data_t value_lookup_data;
			memset(&lookup_data, 0, sizeof(tfx_effect_lookup_data_t));
			memset(&value_lookup_data, 0, sizeof(tfx_effect_lookup_data_t));
			if (current->type == tfxEmitterType) {

				int offset = TFX_GLOBAL_COUNT + TFX_PROPERTY_COUNT + TFX_BASE_COUNT + TFX_VARIATION_COUNT;

				GetEffectInfo(current)->lookup_value_index = library->compiled_lookup_indexes.size();
				for (int i = 0; i != TFX_OVERTIME_COUNT; ++i) {
					tfx_graph_t &graph = ((tfx_graph_t *)(&library->emitter_attributes[current->emitter_attributes].overtime))[i];
					tfx_graph_lookup_index_t &index = ((tfx_graph_lookup_index_t *)&lookup_data)[i];
					index.start_index = running_node_index;
					index.length = graph.nodes.size();
					index.max_life = graph.lookup.life;
					for (tfxBucketLoop(graph.nodes)) {
						library->all_nodes.push_back(graph.nodes[i]);
						running_node_index++;
					}

					tfx_graph_lookup_index_t value_index;
					value_index.start_index = running_value_index;
					value_index.length = graph.lookup.values.size();
					value_index.max_life = graph.lookup.life;
					library->compiled_lookup_indexes.push_back(value_index);
					for (auto value : graph.lookup.values) {
						library->compiled_lookup_values.push_back(value);
						running_value_index++;
					}
				}

				library->node_lookup_indexes.push_back(lookup_data);
				GetEffectInfo(current)->lookup_node_index = library->node_lookup_indexes.size() - 1;

			}

			for (auto &sub : GetEffectInfo(current)->sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}
}

void CompileLibraryGraphsOfEffect(tfx_library_t *library, tfx_effect_emitter_t *effect, tfxU32 depth) {
	tfx_effect_emitter_info_t *info = GetEffectInfo(effect);
	if (effect->type == tfxEffectType && depth == 0) {
		CompileLibraryKeyframeGraph(library, effect->transform_attributes);
		CompileLibraryGlobalGraph(library, effect->global);
	}
	else if (effect->type == tfxEffectType) {
		CompileLibraryKeyframeGraph(library, effect->transform_attributes);
	}
	else if (effect->type == tfxEmitterType) {
		CompileLibraryKeyframeGraph(library, effect->transform_attributes);
		CompileLibraryEmitterGraphs(library, effect->emitter_attributes);
		for (auto &sub : info->sub_effectors) {
			CompileLibraryGraphsOfEffect(library, &sub, ++depth);
		}
	}
	else if (effect->type == tfxFolder) {
		for (auto &sub : info->sub_effectors) {
			CompileLibraryGraphsOfEffect(library, &sub, 0);
		}
	}
	if (effect->type == tfxEffectType) {
		for (auto &sub : info->sub_effectors) {
			CompileLibraryGraphsOfEffect(library, &sub, ++depth);
		}
	}
}

void CompileAllLibraryGraphs(tfx_library_t *library) {
	for (auto &g : library->global_graphs) {
		CompileGraph(&g.amount);
		CompileGraph(&g.height);
		CompileGraph(&g.width);
		CompileGraph(&g.life);
		CompileGraph(&g.intensity);
		CompileGraph(&g.overal_scale);
		CompileGraph(&g.spin);
		CompileGraph(&g.pitch_spin);
		CompileGraph(&g.yaw_spin);
		CompileGraph(&g.splatter);
		CompileGraph(&g.stretch);
		CompileGraph(&g.velocity);
		CompileGraph(&g.noise);
		CompileGraph(&g.weight);
		CompileGraph(&g.emitter_width);
		CompileGraph(&g.emitter_height);
		CompileGraph(&g.emitter_depth);
	}
	for (auto &g : library->transform_attributes) {
		CompileGraph(&g.roll);
		CompileGraph(&g.pitch);
		CompileGraph(&g.yaw);
		CompileGraph(&g.translation_x);
		CompileGraph(&g.translation_y);
		CompileGraph(&g.translation_z);
	}
	for (auto &g : library->emitter_attributes) {
		CompileGraph(&g.properties.arc_offset);
		CompileGraph(&g.properties.arc_size);
		CompileGraph(&g.properties.extrusion);
		CompileGraph(&g.properties.emission_pitch);
		CompileGraph(&g.properties.emission_yaw);
		CompileGraph(&g.properties.emission_range);
		CompileGraph(&g.properties.emitter_width);
		CompileGraph(&g.properties.emitter_height);
		CompileGraph(&g.properties.emitter_depth);
		CompileGraph(&g.properties.splatter);

		CompileGraph(&g.base.amount);
		CompileGraph(&g.base.width);
		CompileGraph(&g.base.height);
		CompileGraph(&g.base.life);
		CompileGraph(&g.base.spin);
		CompileGraph(&g.base.pitch_spin);
		CompileGraph(&g.base.yaw_spin);
		CompileGraph(&g.base.noise_offset);
		CompileGraph(&g.base.velocity);
		CompileGraph(&g.base.weight);

		CompileGraph(&g.variation.amount);
		CompileGraph(&g.variation.width);
		CompileGraph(&g.variation.height);
		CompileGraph(&g.variation.life);
		CompileGraph(&g.variation.noise_offset);
		CompileGraph(&g.variation.noise_resolution);
		CompileGraph(&g.variation.motion_randomness);
		CompileGraph(&g.variation.spin);
		CompileGraph(&g.variation.pitch_spin);
		CompileGraph(&g.variation.yaw_spin);
		CompileGraph(&g.variation.velocity);
		CompileGraph(&g.variation.weight);

		CompileGraphOvertime(&g.overtime.blendfactor);
		CompileGraphOvertime(&g.overtime.intensity);
		CompileGraphOvertime(&g.overtime.alpha_sharpness);
		CompileGraphOvertime(&g.overtime.curved_alpha);
		CompileGraphOvertime(&g.overtime.velocity_turbulance);
		CompileGraphOvertime(&g.overtime.width);
		CompileGraphOvertime(&g.overtime.height);
		CompileGraphOvertime(&g.overtime.direction_turbulance);
		CompileGraphOvertime(&g.overtime.spin);
		CompileGraphOvertime(&g.overtime.pitch_spin);
		CompileGraphOvertime(&g.overtime.yaw_spin);
		CompileGraphOvertime(&g.overtime.stretch);
		CompileGraphOvertime(&g.overtime.velocity);
		CompileGraph(&g.overtime.velocity_adjuster);
		CompileGraphOvertime(&g.overtime.weight);
		CompileGraphOvertime(&g.overtime.direction);
		CompileGraphOvertime(&g.overtime.noise_resolution);
		CompileGraphOvertime(&g.overtime.motion_randomness);

		CompileGraphOvertime(&g.factor.life);
		CompileGraphOvertime(&g.factor.velocity);
		CompileGraphOvertime(&g.factor.size);
		CompileGraphOvertime(&g.factor.intensity);
		g.overtime.color_ramps[0] = CompileColorRamp(&g.overtime);
		g.overtime.color_ramps[1] = CompileColorRampHint(&g.overtime);
	}
	CreateColorRampBitmaps(library);
}

void CompileLibraryGlobalGraph(tfx_library_t *library, tfxU32 index) {
	tfx_global_attributes_t &g = library->global_graphs[index];
	CompileGraph(&g.amount);
	CompileGraph(&g.height);
	CompileGraph(&g.width);
	CompileGraph(&g.life);
	CompileGraph(&g.intensity);
	CompileGraph(&g.overal_scale);
	CompileGraph(&g.spin);
	CompileGraph(&g.pitch_spin);
	CompileGraph(&g.yaw_spin);
	CompileGraph(&g.splatter);
	CompileGraph(&g.stretch);
	CompileGraph(&g.velocity);
	CompileGraph(&g.noise);
	CompileGraph(&g.weight);
	CompileGraph(&g.emitter_width);
	CompileGraph(&g.emitter_height);
	CompileGraph(&g.emitter_depth);
}

void CompileLibraryKeyframeGraph(tfx_library_t *library, tfxU32 index) {
	tfx_transform_attributes_t &g = library->transform_attributes[index];
	CompileGraph(&g.roll);
	CompileGraph(&g.pitch);
	CompileGraph(&g.yaw);
	CompileGraph(&g.translation_x);
	CompileGraph(&g.translation_y);
	CompileGraph(&g.translation_z);
}

void CompileLibraryEmitterGraphs(tfx_library_t *library, tfxU32 index) {
	CompileLibraryPropertyGraph(library, index);
	CompileLibraryKeyframeGraph(library, index);
	CompileLibraryBaseGraph(library, index);
	CompileLibraryVariationGraph(library, index);
	CompileLibraryOvertimeGraph(library, index);
	CompileLibraryFactorGraph(library, index);
}

void CompileLibraryPropertyGraph(tfx_library_t *library, tfxU32 index) {
	tfx_property_attributes_t &g = library->emitter_attributes[index].properties;
	CompileGraph(&g.arc_offset);
	CompileGraph(&g.arc_size);
	CompileGraph(&g.extrusion);
	CompileGraph(&g.emission_pitch);
	CompileGraph(&g.emission_yaw);
	CompileGraph(&g.emission_range);
	CompileGraph(&g.emitter_width);
	CompileGraph(&g.emitter_height);
	CompileGraph(&g.emitter_depth);
	CompileGraph(&g.splatter);
}

void CompileLibraryBaseGraph(tfx_library_t *library, tfxU32 index) {
	tfx_base_attributes_t &g = library->emitter_attributes[index].base;
	CompileGraph(&g.amount);
	CompileGraph(&g.width);
	CompileGraph(&g.height);
	CompileGraph(&g.life);
	CompileGraph(&g.spin);
	CompileGraph(&g.pitch_spin);
	CompileGraph(&g.yaw_spin);
	CompileGraph(&g.noise_offset);
	CompileGraph(&g.velocity);
	CompileGraph(&g.weight);
}

void CompileLibraryVariationGraph(tfx_library_t *library, tfxU32 index) {
	tfx_variation_attributes_t &g = library->emitter_attributes[index].variation;
	CompileGraph(&g.amount);
	CompileGraph(&g.width);
	CompileGraph(&g.height);
	CompileGraph(&g.life);
	CompileGraph(&g.noise_offset);
	CompileGraph(&g.noise_resolution);
	CompileGraph(&g.motion_randomness);
	CompileGraph(&g.spin);
	CompileGraph(&g.pitch_spin);
	CompileGraph(&g.yaw_spin);
	CompileGraph(&g.velocity);
	CompileGraph(&g.weight);
}

void CompileLibraryOvertimeGraph(tfx_library_t *library, tfxU32 index, bool including_color_ramps) {
	tfx_overtime_attributes_t &g = library->emitter_attributes[index].overtime;
	if (including_color_ramps) {
		g.color_ramps[0] = CompileColorRamp(&g);
		g.color_ramps[1] = CompileColorRampHint(&g);
		EditColorRampBitmap(library, &g, 0);
		EditColorRampBitmap(library, &g, 1);
	}
	CompileGraphOvertime(&g.intensity);
	CompileGraphOvertime(&g.alpha_sharpness);
	CompileGraphOvertime(&g.curved_alpha);
	CompileGraphOvertime(&g.velocity_turbulance);
	CompileGraphOvertime(&g.width);
	CompileGraphOvertime(&g.height);
	CompileGraphOvertime(&g.direction_turbulance);
	CompileGraphOvertime(&g.spin);
	CompileGraphOvertime(&g.pitch_spin);
	CompileGraphOvertime(&g.yaw_spin);
	CompileGraphOvertime(&g.stretch);
	CompileGraphOvertime(&g.velocity);
	CompileGraph(&g.velocity_adjuster);
	CompileGraphOvertime(&g.weight);
	CompileGraphOvertime(&g.direction);
	CompileGraphOvertime(&g.noise_resolution);
	CompileGraphOvertime(&g.motion_randomness);
}

void CompileLibraryFactorGraph(tfx_library_t *library, tfxU32 index) {
	tfx_factor_attributes_t &g = library->emitter_attributes[index].factor;
	CompileGraphOvertime(&g.life);
	CompileGraphOvertime(&g.velocity);
	CompileGraphOvertime(&g.size);
	CompileGraphOvertime(&g.intensity);
}

void CompileLibraryColorGraphs(tfx_library_t *library, tfxU32 index) {
	tfx_overtime_attributes_t &g = library->emitter_attributes[index].overtime;
	g.color_ramps[0] = CompileColorRamp(&g);
	g.color_ramps[1] = CompileColorRampHint(&g);
	EditColorRampBitmap(library, &g, 0);
	EditColorRampBitmap(library, &g, 1);
}

void SetLibraryMinMaxData(tfx_library_t *library) {
	library->graph_min_max.clear();
	library->graph_min_max.create_pool(tfxEmitterGraphMaxIndex);

	library->graph_min_max[tfxGlobal_life] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
	library->graph_min_max[tfxGlobal_amount] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
	library->graph_min_max[tfxGlobal_velocity] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
	library->graph_min_max[tfxGlobal_noise] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
	library->graph_min_max[tfxGlobal_width] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
	library->graph_min_max[tfxGlobal_height] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
	library->graph_min_max[tfxGlobal_weight] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
	library->graph_min_max[tfxGlobal_roll_spin] = GetMinMaxGraphValues(tfxGlobalPercentPresetSigned);
	library->graph_min_max[tfxGlobal_pitch_spin] = GetMinMaxGraphValues(tfxGlobalPercentPresetSigned);
	library->graph_min_max[tfxGlobal_yaw_spin] = GetMinMaxGraphValues(tfxGlobalPercentPresetSigned);
	library->graph_min_max[tfxGlobal_stretch] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
	library->graph_min_max[tfxGlobal_overal_scale] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
	library->graph_min_max[tfxGlobal_intensity] = GetMinMaxGraphValues(tfxOpacityOvertimePreset);
	library->graph_min_max[tfxGlobal_splatter] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
	library->graph_min_max[tfxGlobal_emitter_width] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
	library->graph_min_max[tfxGlobal_emitter_height] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
	library->graph_min_max[tfxGlobal_emitter_depth] = GetMinMaxGraphValues(tfxGlobalPercentPreset);

	library->graph_min_max[tfxTransform_roll] = GetMinMaxGraphValues(tfxAnglePreset);
	library->graph_min_max[tfxTransform_pitch] = GetMinMaxGraphValues(tfxAnglePreset);
	library->graph_min_max[tfxTransform_yaw] = GetMinMaxGraphValues(tfxAnglePreset);
	library->graph_min_max[tfxTransform_translate_x] = GetMinMaxGraphValues(tfxDimensionsPreset);
	library->graph_min_max[tfxTransform_translate_y] = GetMinMaxGraphValues(tfxDimensionsPreset);
	library->graph_min_max[tfxTransform_translate_z] = GetMinMaxGraphValues(tfxDimensionsPreset);

	library->graph_min_max[tfxProperty_emission_pitch] = GetMinMaxGraphValues(tfxAnglePreset);
	library->graph_min_max[tfxProperty_emission_yaw] = GetMinMaxGraphValues(tfxAnglePreset);
	library->graph_min_max[tfxProperty_emission_range] = GetMinMaxGraphValues(tfxEmissionRangePreset);
	library->graph_min_max[tfxProperty_splatter] = GetMinMaxGraphValues(tfxDimensionsPreset);
	library->graph_min_max[tfxProperty_emitter_width] = GetMinMaxGraphValues(tfxDimensionsPreset);
	library->graph_min_max[tfxProperty_emitter_height] = GetMinMaxGraphValues(tfxDimensionsPreset);
	library->graph_min_max[tfxProperty_emitter_depth] = GetMinMaxGraphValues(tfxDimensionsPreset);
	library->graph_min_max[tfxProperty_extrusion] = GetMinMaxGraphValues(tfxDimensionsPreset);
	library->graph_min_max[tfxProperty_arc_size] = GetMinMaxGraphValues(tfxArcPreset);
	library->graph_min_max[tfxProperty_arc_offset] = GetMinMaxGraphValues(tfxArcPreset);

	library->graph_min_max[tfxBase_life] = GetMinMaxGraphValues(tfxLifePreset);
	library->graph_min_max[tfxBase_amount] = GetMinMaxGraphValues(tfxAmountPreset);
	library->graph_min_max[tfxBase_velocity] = GetMinMaxGraphValues(tfxVelocityPreset);
	library->graph_min_max[tfxBase_width] = GetMinMaxGraphValues(tfxDimensionsPreset);
	library->graph_min_max[tfxBase_height] = GetMinMaxGraphValues(tfxDimensionsPreset);
	library->graph_min_max[tfxBase_weight] = GetMinMaxGraphValues(tfxWeightPreset);
	library->graph_min_max[tfxBase_roll_spin] = GetMinMaxGraphValues(tfxSpinPreset);
	library->graph_min_max[tfxBase_pitch_spin] = GetMinMaxGraphValues(tfxSpinPreset);
	library->graph_min_max[tfxBase_yaw_spin] = GetMinMaxGraphValues(tfxSpinPreset);
	library->graph_min_max[tfxBase_noise_offset] = GetMinMaxGraphValues(tfxDimensionsPreset);

	library->graph_min_max[tfxVariation_life] = GetMinMaxGraphValues(tfxLifePreset);
	library->graph_min_max[tfxVariation_amount] = GetMinMaxGraphValues(tfxAmountPreset);
	library->graph_min_max[tfxVariation_velocity] = GetMinMaxGraphValues(tfxVelocityPreset);
	library->graph_min_max[tfxVariation_width] = GetMinMaxGraphValues(tfxDimensionsPreset);
	library->graph_min_max[tfxVariation_height] = GetMinMaxGraphValues(tfxDimensionsPreset);
	library->graph_min_max[tfxVariation_weight] = GetMinMaxGraphValues(tfxWeightVariationPreset);
	library->graph_min_max[tfxVariation_roll_spin] = GetMinMaxGraphValues(tfxSpinVariationPreset);
	library->graph_min_max[tfxVariation_pitch_spin] = GetMinMaxGraphValues(tfxSpinVariationPreset);
	library->graph_min_max[tfxVariation_yaw_spin] = GetMinMaxGraphValues(tfxSpinVariationPreset);
	library->graph_min_max[tfxVariation_noise_offset] = GetMinMaxGraphValues(tfxGlobalPercentPreset);

	library->graph_min_max[tfxOvertime_velocity] = GetMinMaxGraphValues(tfxVelocityOvertimePreset);
	library->graph_min_max[tfxOvertime_width] = GetMinMaxGraphValues(tfxPercentOvertime);
	library->graph_min_max[tfxOvertime_height] = GetMinMaxGraphValues(tfxPercentOvertime);
	library->graph_min_max[tfxOvertime_weight] = GetMinMaxGraphValues(tfxPercentOvertime);
	library->graph_min_max[tfxOvertime_roll_spin] = GetMinMaxGraphValues(tfxSpinOvertimePreset);
	library->graph_min_max[tfxOvertime_pitch_spin] = GetMinMaxGraphValues(tfxSpinOvertimePreset);
	library->graph_min_max[tfxOvertime_yaw_spin] = GetMinMaxGraphValues(tfxSpinOvertimePreset);
	library->graph_min_max[tfxOvertime_stretch] = GetMinMaxGraphValues(tfxPercentOvertime);
	library->graph_min_max[tfxOvertime_red] = GetMinMaxGraphValues(tfxColorPreset);
	library->graph_min_max[tfxOvertime_green] = GetMinMaxGraphValues(tfxColorPreset);
	library->graph_min_max[tfxOvertime_blue] = GetMinMaxGraphValues(tfxColorPreset);
	library->graph_min_max[tfxOvertime_blendfactor] = GetMinMaxGraphValues(tfxOpacityOvertimePreset);
	library->graph_min_max[tfxOvertime_intensity] = GetMinMaxGraphValues(tfxIntensityOvertimePreset);
	library->graph_min_max[tfxOvertime_red_hint] = GetMinMaxGraphValues(tfxColorPreset);
	library->graph_min_max[tfxOvertime_green_hint] = GetMinMaxGraphValues(tfxColorPreset);
	library->graph_min_max[tfxOvertime_blue_hint] = GetMinMaxGraphValues(tfxColorPreset);
	library->graph_min_max[tfxOvertime_blendfactor_hint] = GetMinMaxGraphValues(tfxOpacityOvertimePreset);
	library->graph_min_max[tfxOvertime_alpha_sharpness] = GetMinMaxGraphValues(tfxIntensityOvertimePreset);
	library->graph_min_max[tfxOvertime_curved_alpha] = GetMinMaxGraphValues(tfxOpacityOvertimePreset);
	library->graph_min_max[tfxOvertime_velocity_turbulance] = GetMinMaxGraphValues(tfxFrameratePreset);
	library->graph_min_max[tfxOvertime_direction_turbulance] = GetMinMaxGraphValues(tfxPercentOvertime);
	library->graph_min_max[tfxOvertime_velocity_adjuster] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
	library->graph_min_max[tfxOvertime_direction] = GetMinMaxGraphValues(tfxDirectionOvertimePreset);
}

float LookupLibraryPreciseOvertimeNodeList(tfx_library_t *library, tfx_graph_type graph_type, int lookup_node_index, float age, float life) {
	float lastv = 0;
	float lastf = 0;
	float p = 0;
	tfx_attribute_node_t *lastec = nullptr;
	tfx_graph_lookup_index_t &lookup_data = ((tfx_graph_lookup_index_t *)&library->node_lookup_indexes[lookup_node_index])[graph_type];
	float min_y = library->graph_min_max[graph_type].y;
	float max_y = library->graph_min_max[graph_type].w;
	for (int i = lookup_data.start_index; i != lookup_data.start_index + lookup_data.length; ++i) {
		tfx_attribute_node_t &a = library->all_nodes[i];
		float frame = a.frame * life;
		if (age < frame) {
			p = (age - lastf) / (frame - lastf);
			float bezier_value = GetBezierValue(lastec, &a, p, min_y, max_y);
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

float LookupLibraryPreciseNodeList(tfx_library_t *library, tfx_graph_type graph_type, int lookup_node_index, float age) {
	float lastv = 0;
	float lastf = 0;
	float p = 0;
	tfx_attribute_node_t *lastec = nullptr;
	tfx_graph_lookup_index_t &lookup_data = ((tfx_graph_lookup_index_t *)&library->node_lookup_indexes[lookup_node_index])[graph_type];
	float min_y = library->graph_min_max[graph_type].y;
	float max_y = library->graph_min_max[graph_type].w;
	for (int i = lookup_data.start_index; i != lookup_data.start_index + lookup_data.length; ++i) {
		tfx_attribute_node_t &a = library->all_nodes[i];
		if (age < a.frame) {
			p = (age - lastf) / (a.frame - lastf);
			float bezier_value = GetBezierValue(lastec, &a, p, min_y, max_y);
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

float LookupLibraryFastValueList(tfx_library_t *library, tfx_graph_type graph_type, int lookup_node_index, float frame) {
	tfx_graph_lookup_index_t &lookup_data = ((tfx_graph_lookup_index_t *)&library->compiled_lookup_indexes[lookup_node_index])[graph_type];
	frame += lookup_data.start_index;
	tfxU32 end_frame = lookup_data.start_index + lookup_data.length - 1;
	frame = frame > end_frame ? end_frame : frame;
	return library->compiled_lookup_values[(tfxU32)frame];
}

float LookupLibraryFastOvertimeValueList(tfx_library_t *library, tfx_graph_type graph_type, int lookup_value_index, float age, float lifetime) {
	tfx_graph_lookup_index_t &lookup_data = ((tfx_graph_lookup_index_t *)&library->compiled_lookup_indexes[lookup_value_index])[graph_type - tfxOvertime_velocity];
	float frame = (float)lookup_data.start_index;
	if (lifetime)
		frame += (age / lifetime * lookup_data.max_life) / tfxLOOKUP_FREQUENCY_OVERTIME;
	if (frame < lookup_data.start_index + lookup_data.length - 1)
		return library->compiled_lookup_values[(tfxU32)frame];
	return library->compiled_lookup_values[lookup_data.start_index + lookup_data.length - 1];
}

tfxU32 CountOfGraphsInUse(tfx_library_t *library) {
	return library->global_graphs.size() + library->emitter_attributes.size() - CountOfFreeGraphs(library);
}

tfxU32 CountOfFreeGraphs(tfx_library_t *library) {
	return library->free_global_graphs.size() + library->free_emitter_attributes.size();
}

void tfx_data_types_dictionary_t::Init() {
	names_and_types.data.reserve(300);
	names_and_types.map.reserve(300);
	names_and_types.Insert("name", tfxString);
	names_and_types.Insert("image_index", tfxUint);
	names_and_types.Insert("image_hash", tfxUInt64);
	names_and_types.Insert("image_handle_x", tfxFloat);
	names_and_types.Insert("image_handle_y", tfxFloat);
	names_and_types.Insert("spawn_amount", tfxUint);
	names_and_types.Insert("spawn_amount_variation", tfxUint);
	names_and_types.Insert("single_shot_limit", tfxUint);
	names_and_types.Insert("blend_mode", tfxSInt);
	names_and_types.Insert("image_start_frame", tfxFloat);
	names_and_types.Insert("image_end_frame", tfxFloat);
	names_and_types.Insert("image_frame_rate", tfxFloat);

	names_and_types.Insert("emission_type", tfxSInt);
	names_and_types.Insert("emission_direction", tfxSInt);
	names_and_types.Insert("delay_spawning", tfxFloat);
	names_and_types.Insert("grid_rows", tfxFloat);
	names_and_types.Insert("grid_columns", tfxFloat);
	names_and_types.Insert("grid_depth", tfxFloat);
	names_and_types.Insert("loop_length", tfxFloat);
	names_and_types.Insert("noise_base_offset_range", tfxFloat);
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
	names_and_types.Insert("paired_emitter_hash", tfxUInt64);

	names_and_types.Insert("random_color", tfxBool);
	names_and_types.Insert("exclude_from_global_hue", tfxBool);
	names_and_types.Insert("relative_position", tfxBool);
	names_and_types.Insert("relative_angle", tfxBool);
	names_and_types.Insert("image_handle_auto_center", tfxBool);
	names_and_types.Insert("single", tfxBool);
	names_and_types.Insert("wrap_single_sprite", tfxBool);
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
	names_and_types.Insert("include_in_sprite_data_export", tfxBool);
	names_and_types.Insert("use_path_for_direction", tfxBool);
	names_and_types.Insert("alt_velocity_lifetime_sampling", tfxBool);
	names_and_types.Insert("alt_color_lifetime_sampling", tfxBool);
	names_and_types.Insert("alt_size_lifetime_sampling", tfxBool);
	names_and_types.Insert("use_simple_motion_randomness", tfxBool);
	names_and_types.Insert("spawn_location_source", tfxBool);
	names_and_types.Insert("use_color_hint", tfxBool);
	//names_and_types.Insert("simple_motion_smoothstep", tfxBool);

	//Graphs
	names_and_types.Insert("global_life", tfxFloat);
	names_and_types.Insert("global_amount", tfxFloat);
	names_and_types.Insert("global_velocity", tfxFloat);
	names_and_types.Insert("global_noise", tfxFloat);
	names_and_types.Insert("global_width", tfxFloat);
	names_and_types.Insert("global_height", tfxFloat);
	names_and_types.Insert("global_weight", tfxFloat);
	names_and_types.Insert("global_spin", tfxFloat);
	names_and_types.Insert("global_roll_spin", tfxFloat);
	names_and_types.Insert("global_pitch_spin", tfxFloat);
	names_and_types.Insert("global_yaw_spin", tfxFloat);
	names_and_types.Insert("global_stretch", tfxFloat);
	names_and_types.Insert("global_overal_scale", tfxFloat);
	names_and_types.Insert("global_intensity", tfxFloat);
	names_and_types.Insert("global_splatter", tfxFloat);
	names_and_types.Insert("global_emitter_width", tfxFloat);
	names_and_types.Insert("global_emitter_height", tfxFloat);
	names_and_types.Insert("global_emitter_depth", tfxFloat);

	names_and_types.Insert("property_emission_pitch", tfxFloat);
	names_and_types.Insert("property_emission_yaw", tfxFloat);
	names_and_types.Insert("property_emission_range", tfxFloat);
	names_and_types.Insert("property_splatter", tfxFloat);
	names_and_types.Insert("property_emitter_width", tfxFloat);
	names_and_types.Insert("property_emitter_height", tfxFloat);
	names_and_types.Insert("property_emitter_depth", tfxFloat);
	names_and_types.Insert("property_extrusion", tfxFloat);
	names_and_types.Insert("property_arc_size", tfxFloat);
	names_and_types.Insert("property_arc_offset", tfxFloat);

	names_and_types.Insert("base_life", tfxFloat);
	names_and_types.Insert("base_amount", tfxFloat);
	names_and_types.Insert("base_velocity", tfxFloat);
	names_and_types.Insert("base_width", tfxFloat);
	names_and_types.Insert("base_height", tfxFloat);
	names_and_types.Insert("base_weight", tfxFloat);
	names_and_types.Insert("base_spin", tfxFloat);
	names_and_types.Insert("base_roll_spin", tfxFloat);
	names_and_types.Insert("base_pitch_spin", tfxFloat);
	names_and_types.Insert("base_yaw_spin", tfxFloat);
	names_and_types.Insert("base_noise_offset", tfxFloat);

	names_and_types.Insert("variation_life", tfxFloat);
	names_and_types.Insert("variation_amount", tfxFloat);
	names_and_types.Insert("variation_velocity", tfxFloat);
	names_and_types.Insert("variation_width", tfxFloat);
	names_and_types.Insert("variation_height", tfxFloat);
	names_and_types.Insert("variation_weight", tfxFloat);
	names_and_types.Insert("variation_spin", tfxFloat);
	names_and_types.Insert("variation_roll_spin", tfxFloat);
	names_and_types.Insert("variation_pitch_spin", tfxFloat);
	names_and_types.Insert("variation_yaw_spin", tfxFloat);
	names_and_types.Insert("variation_noise_offset", tfxFloat);
	names_and_types.Insert("variation_noise_resolution", tfxFloat);
	names_and_types.Insert("variation_motion_randomness", tfxFloat);

	names_and_types.Insert("overtime_velocity", tfxFloat);
	names_and_types.Insert("overtime_width", tfxFloat);
	names_and_types.Insert("overtime_height", tfxFloat);
	names_and_types.Insert("overtime_weight", tfxFloat);
	names_and_types.Insert("overtime_spin", tfxFloat);
	names_and_types.Insert("overtime_roll_spin", tfxFloat);
	names_and_types.Insert("overtime_pitch_spin", tfxFloat);
	names_and_types.Insert("overtime_yaw_spin", tfxFloat);
	names_and_types.Insert("overtime_stretch", tfxFloat);
	names_and_types.Insert("overtime_red", tfxFloat);
	names_and_types.Insert("overtime_green", tfxFloat);
	names_and_types.Insert("overtime_blue", tfxFloat);
	names_and_types.Insert("overtime_opacity", tfxFloat);    //Legacy
	names_and_types.Insert("overtime_blendfactor", tfxFloat);
	names_and_types.Insert("overtime_intensity", tfxFloat);
	names_and_types.Insert("overtime_red_hint", tfxFloat);
	names_and_types.Insert("overtime_green_hint", tfxFloat);
	names_and_types.Insert("overtime_blue_hint", tfxFloat);
	names_and_types.Insert("overtime_blendfactor_hint", tfxFloat);
	names_and_types.Insert("overtime_alpha_sharpness", tfxFloat);
	names_and_types.Insert("overtime_curved_alpha", tfxFloat);
	names_and_types.Insert("overtime_velocity_turbulance", tfxFloat);
	names_and_types.Insert("overtime_direction_turbulance", tfxFloat);
	names_and_types.Insert("overtime_velocity_adjuster", tfxFloat);
	names_and_types.Insert("overtime_direction", tfxFloat);
	names_and_types.Insert("overtime_noise_resolution", tfxFloat);
	names_and_types.Insert("overtime_motion_randomness", tfxFloat);

	names_and_types.Insert("factor_life", tfxFloat);
	names_and_types.Insert("factor_velocity", tfxFloat);
	names_and_types.Insert("factor_size", tfxFloat);
	names_and_types.Insert("factor_intensity", tfxFloat);

	names_and_types.Insert("transform_roll", tfxFloat);
	names_and_types.Insert("transform_pitch", tfxFloat);
	names_and_types.Insert("transform_yaw", tfxFloat);
	names_and_types.Insert("transform_translate_x", tfxFloat);
	names_and_types.Insert("transform_translate_y", tfxFloat);
	names_and_types.Insert("transform_translate_z", tfxFloat);

	names_and_types.Insert("path_pitch", tfxFloat);
	names_and_types.Insert("path_yaw", tfxFloat);
	names_and_types.Insert("path_roll", tfxFloat);
	names_and_types.Insert("path_offset_x", tfxFloat);
	names_and_types.Insert("path_offset_y", tfxFloat);
	names_and_types.Insert("path_offset_z", tfxFloat);
	names_and_types.Insert("distance", tfxFloat);
	names_and_types.Insert("path_mode_origin", tfxBool);
	names_and_types.Insert("path_mode_node", tfxBool);
	names_and_types.Insert("path_node_count", tfxUint);
	names_and_types.Insert("path_is_2d", tfxBool);
	names_and_types.Insert("path_is_3d", tfxBool);    //Not used
	names_and_types.Insert("path_space_nodes_evenly", tfxBool);
	names_and_types.Insert("path_reverse_direction", tfxBool);
	names_and_types.Insert("path_extrusion_type", tfxSInt);
	names_and_types.Insert("path_generator_type", tfxSInt);
	names_and_types.Insert("path_rotation_range", tfxFloat);
	names_and_types.Insert("path_rotation_pitch", tfxFloat);
	names_and_types.Insert("path_rotation_yaw", tfxFloat);
	names_and_types.Insert("maximum_active_paths", tfxUint);
	names_and_types.Insert("maximum_path_cycles", tfxUint);
	names_and_types.Insert("path_rotation_stagger", tfxFloat);
	names_and_types.Insert("path_rotation_range_yaw_only", tfxBool);
	names_and_types.Insert("path_handle_x", tfxFloat);
	names_and_types.Insert("path_handle_y", tfxFloat);
	names_and_types.Insert("path_handle_z", tfxFloat);

	//Sprite data settings
	names_and_types.Insert("start_offset", tfxUint);
	names_and_types.Insert("frames_after_compression", tfxUint);
	names_and_types.Insert("real_frames", tfxUint);
	names_and_types.Insert("frame_count", tfxUint);
	names_and_types.Insert("total_sprites", tfxUint);
	names_and_types.Insert("total_memory_for_sprites", tfxUint);
	names_and_types.Insert("flags", tfxUint);
	names_and_types.Insert("animation_flags", tfxUint);
	names_and_types.Insert("animation_time", tfxFloat);
	names_and_types.Insert("animation_length_in_time", tfxFloat);
	names_and_types.Insert("name", tfxString);
	names_and_types.Insert("path_hash", tfxUInt64);

	//Frame meta
	names_and_types.Insert("total_sprites", tfxUint);
	names_and_types.Insert("min_corner", tfxFloat3);
	names_and_types.Insert("max_corner", tfxFloat3);

	//Sprite data emitter properties
	names_and_types.Insert("animation_frames", tfxFloat);
	names_and_types.Insert("handle", tfxFloat2);
	names_and_types.Insert("start_frame_index", tfxUint);

	//Animation settings
	names_and_types.Insert("playback_speed", tfxFloat);
	names_and_types.Insert("animation_magenta_mask", tfxBool);
	names_and_types.Insert("animation_magenta_mask_always", tfxBool);
	names_and_types.Insert("frames", tfxUint);
	names_and_types.Insert("current_frame", tfxUint);
	names_and_types.Insert("frame_offset", tfxUint);
	names_and_types.Insert("extra_frames_count", tfxSInt);
	names_and_types.Insert("layer", tfxUint);
	names_and_types.Insert("position_x", tfxFloat);
	names_and_types.Insert("position_y", tfxFloat);
	names_and_types.Insert("position_z", tfxFloat);
	names_and_types.Insert("frame_width", tfxFloat);
	names_and_types.Insert("frame_height", tfxFloat);
	names_and_types.Insert("animation_flags", tfxUint);
	names_and_types.Insert("loop", tfxBool);
	names_and_types.Insert("seamless", tfxBool);
	names_and_types.Insert("seed", tfxUint);
	names_and_types.Insert("zoom", tfxFloat);
	names_and_types.Insert("scale", tfxFloat);
	names_and_types.Insert("sprite_data_seed", tfxUint);
	names_and_types.Insert("sprite_data_flags", tfxUint);
	names_and_types.Insert("sprite_data_frames", tfxUint);
	names_and_types.Insert("sprite_data_frame_offset", tfxUint);
	names_and_types.Insert("sprite_data_extra_frames_count", tfxUint);
	names_and_types.Insert("sprite_data_playback_speed", tfxFloat);
	names_and_types.Insert("sprite_data_recording_frame_rate", tfxFloat);
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
	names_and_types.Insert("orthographic_camera_position_x", tfxFloat);
	names_and_types.Insert("orthographic_camera_position_y", tfxFloat);
	names_and_types.Insert("orthographic_camera_position_z", tfxFloat);
	names_and_types.Insert("orthographic_camera_pitch", tfxFloat);
	names_and_types.Insert("orthographic_camera_yaw", tfxFloat);
	names_and_types.Insert("orthographic_camera_fov", tfxFloat);
	names_and_types.Insert("orthographic_camera_floor_height", tfxFloat);
	names_and_types.Insert("orthographic_camera_isometric", tfxBool);
	names_and_types.Insert("orthographic_camera_isometric_scale", tfxFloat);
	names_and_types.Insert("orthographic_camera_hide_floor", tfxBool);
	names_and_types.Insert("orthographic_camera_free_speed", tfxFloat);
	names_and_types.Insert("orthographic_camera_ray_offset", tfxFloat);
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

	initialised = true;
}

int ValidateEffectPackage(const char *filename) {
	tfx_package_t package;
	tfxErrorFlags status = tfx__load_package_file(filename, &package);
	if (status) {
		tfx__free_package(&package);
		return status;                    //returns 1 to 4 if it's an invalid package format
	}

	tfx_package_entry_info_t *data_txt = tfx__get_package_file(&package, "data.txt");
	tfx__free_package(&package);
	if (!data_txt) return tfxErrorCode_data_could_not_be_loaded;                    //Unable to load the the data.txt file in the package

	return 0;
}

void tfx__assign_graph_data(tfx_effect_emitter_t *effect, tfx_vector_t<tfx_str256_t> *values) {
	if (values->size() > 0) {
		if ((*values)[0] == "global_amount") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].amount, &n); }
		if ((*values)[0] == "global_height") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].height, &n); }
		if ((*values)[0] == "global_width") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].width, &n); }
		if ((*values)[0] == "global_life") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].life, &n); }
		if ((*values)[0] == "global_opacity") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].intensity, &n); }
		if ((*values)[0] == "global_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].spin, &n); }
		if ((*values)[0] == "global_roll_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].spin, &n); }
		if ((*values)[0] == "global_pitch_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].pitch_spin, &n); }
		if ((*values)[0] == "global_yaw_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].yaw_spin, &n); }
		if ((*values)[0] == "global_splatter") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].splatter, &n); }
		if ((*values)[0] == "global_stretch") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].stretch, &n); }
		if ((*values)[0] == "global_overal_scale") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].overal_scale, &n); }
		if ((*values)[0] == "global_weight") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].weight, &n); }
		if ((*values)[0] == "global_velocity") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].velocity, &n); }
		if ((*values)[0] == "global_noise") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].noise, &n); }
		if ((*values)[0] == "global_emitter_width") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].emitter_width, &n); }
		if ((*values)[0] == "global_emitter_height") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].emitter_height, &n); }
		if ((*values)[0] == "global_emitter_depth") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->global_graphs[effect->global].emitter_depth, &n); }

		if ((*values)[0] == "global_effect_angle") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->transform_attributes[effect->transform_attributes].roll, &n); }
		if ((*values)[0] == "global_effect_roll") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->transform_attributes[effect->transform_attributes].roll, &n); }
		if ((*values)[0] == "global_effect_pitch") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->transform_attributes[effect->transform_attributes].pitch, &n); }
		if ((*values)[0] == "global_effect_yaw") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->transform_attributes[effect->transform_attributes].yaw, &n); }
		if ((*values)[0] == "keyframe_translate_x") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->transform_attributes[effect->transform_attributes].translation_x, &n); }
		if ((*values)[0] == "keyframe_translate_y") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->transform_attributes[effect->transform_attributes].translation_y, &n); }
		if ((*values)[0] == "keyframe_translate_z") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->transform_attributes[effect->transform_attributes].translation_z, &n); }
		if ((*values)[0] == "property_emitter_angle") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->transform_attributes[effect->transform_attributes].roll, &n); }
		if ((*values)[0] == "property_emitter_roll") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->transform_attributes[effect->transform_attributes].roll, &n); }
		if ((*values)[0] == "property_emitter_pitch") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->transform_attributes[effect->transform_attributes].pitch, &n); }
		if ((*values)[0] == "property_emitter_yaw") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->transform_attributes[effect->transform_attributes].yaw, &n); }

		if ((*values)[0] == "base_arc_offset") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.arc_offset, &n); }
		if ((*values)[0] == "base_arc_size") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.arc_size, &n); }
		if ((*values)[0] == "base_emission_angle") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.emission_pitch, &n); }
		if ((*values)[0] == "base_emission_range") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.emission_range, &n); }
		if ((*values)[0] == "base_emitter_height") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.emitter_height, &n); }
		if ((*values)[0] == "base_emitter_width") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.emitter_width, &n); }
		if ((*values)[0] == "base_splatter") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.splatter, &n); }

		if ((*values)[0] == "property_arc_offset") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.arc_offset, &n); }
		if ((*values)[0] == "property_arc_size") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.arc_size, &n); }
		if ((*values)[0] == "property_extrusion") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.extrusion, &n); }
		if ((*values)[0] == "property_emission_angle") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.emission_pitch, &n); }
		if ((*values)[0] == "property_emission_pitch") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.emission_pitch, &n); }
		if ((*values)[0] == "property_emission_yaw") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.emission_yaw, &n); }
		if ((*values)[0] == "property_emission_range") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.emission_range, &n); }
		if ((*values)[0] == "property_emitter_height") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.emitter_height, &n); }
		if ((*values)[0] == "property_emitter_width") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.emitter_width, &n); }
		if ((*values)[0] == "property_emitter_depth") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.emitter_depth, &n); }
		if ((*values)[0] == "property_splatter") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].properties.splatter, &n); }

		if ((*values)[0] == "base_amount") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].base.amount, &n); }
		if ((*values)[0] == "base_life") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].base.life, &n); }
		if ((*values)[0] == "base_height") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].base.height, &n); }
		if ((*values)[0] == "base_width") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].base.width, &n); }
		if ((*values)[0] == "base_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].base.spin, &n); }
		if ((*values)[0] == "base_roll_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].base.spin, &n); }
		if ((*values)[0] == "base_pitch_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].base.pitch_spin, &n); }
		if ((*values)[0] == "base_yaw_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].base.yaw_spin, &n); }
		if ((*values)[0] == "base_noise_offset") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].base.noise_offset, &n); }
		if ((*values)[0] == "base_velocity") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].base.velocity, &n); }
		if ((*values)[0] == "base_weight") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].base.weight, &n); }

		if ((*values)[0] == "variation_amount") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].variation.amount, &n); }
		if ((*values)[0] == "variation_height") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].variation.height, &n); }
		if ((*values)[0] == "variation_width") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].variation.width, &n); }
		if ((*values)[0] == "variation_life") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].variation.life, &n); }
		if ((*values)[0] == "variation_velocity") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].variation.velocity, &n); }
		if ((*values)[0] == "variation_weight") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].variation.weight, &n); }
		if ((*values)[0] == "variation_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].variation.spin, &n); }
		if ((*values)[0] == "variation_roll_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].variation.spin, &n); }
		if ((*values)[0] == "variation_pitch_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].variation.pitch_spin, &n); }
		if ((*values)[0] == "variation_yaw_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].variation.yaw_spin, &n); }
		if ((*values)[0] == "variation_noise_offset") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].variation.noise_offset, &n); }
		if ((*values)[0] == "variation_noise_resolution") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].variation.noise_resolution, &n); }
		if ((*values)[0] == "variation_motion_randomness") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].variation.motion_randomness, &n); }

		if ((*values)[0] == "overtime_red") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.red, &n); }
		if ((*values)[0] == "overtime_green") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.green, &n); }
		if ((*values)[0] == "overtime_blue") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.blue, &n); }
		if ((*values)[0] == "overtime_blendfactor") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.blendfactor, &n); }
		if ((*values)[0] == "overtime_opacity") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.blendfactor, &n); }    //Legacy
		if ((*values)[0] == "overtime_intensity") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.intensity, &n); }
		if ((*values)[0] == "overtime_red_hint") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.red_hint, &n); }
		if ((*values)[0] == "overtime_green_hint") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.green_hint, &n); }
		if ((*values)[0] == "overtime_blue_hint") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.blue_hint, &n); }
		if ((*values)[0] == "overtime_blendfactor_hint") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.blendfactor_hint, &n); }
		if ((*values)[0] == "overtime_alpha_sharpness") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.alpha_sharpness, &n); }
		if ((*values)[0] == "overtime_curved_alpha") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.curved_alpha, &n); }
		if ((*values)[0] == "overtime_intensity_hint") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.alpha_sharpness, &n); }
		if ((*values)[0] == "overtime_color_mix_balance") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.curved_alpha, &n); }
		if ((*values)[0] == "overtime_velocity_turbulance") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.velocity_turbulance, &n); }
		if ((*values)[0] == "overtime_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.spin, &n); }
		if ((*values)[0] == "overtime_roll_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.spin, &n); }
		if ((*values)[0] == "overtime_pitch_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.pitch_spin, &n); }
		if ((*values)[0] == "overtime_yaw_spin") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.yaw_spin, &n); }
		if ((*values)[0] == "overtime_stretch") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.stretch, &n); }
		if ((*values)[0] == "overtime_velocity") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.velocity, &n); }
		if ((*values)[0] == "overtime_weight") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.weight, &n); }
		if ((*values)[0] == "overtime_width") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.width, &n); }
		if ((*values)[0] == "overtime_height") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.height, &n); }
		if ((*values)[0] == "overtime_direction_turbulance") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.direction_turbulance, &n); }
		if ((*values)[0] == "overtime_direction") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.direction, &n); }
		if ((*values)[0] == "overtime_velocity_adjuster") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.velocity_adjuster, &n); }
		if ((*values)[0] == "overtime_noise_resolution") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.noise_resolution, &n); }
		if ((*values)[0] == "overtime_motion_randomness") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.motion_randomness, &n); }

		if ((*values)[0] == "factor_life") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].factor.life, &n); }
		if ((*values)[0] == "factor_velocity") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].factor.velocity, &n); }
		if ((*values)[0] == "factor_size") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].factor.size, &n); }
		if ((*values)[0] == "factor_intensity") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].factor.intensity, &n); }

		if ((*values)[0] == "path_pitch") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->paths[CreateEmitterPathAttributes(effect, false)].angle_x, &n); }
		if ((*values)[0] == "path_yaw") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->paths[CreateEmitterPathAttributes(effect, false)].angle_y, &n); }
		if ((*values)[0] == "path_roll") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->paths[CreateEmitterPathAttributes(effect, false)].angle_z, &n); }
		if ((*values)[0] == "path_offset_x") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->paths[CreateEmitterPathAttributes(effect, false)].offset_x, &n); }
		if ((*values)[0] == "path_offset_y") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->paths[CreateEmitterPathAttributes(effect, false)].offset_y, &n); }
		if ((*values)[0] == "path_offset_z") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->paths[CreateEmitterPathAttributes(effect, false)].offset_z, &n); }
		if ((*values)[0] == "path_distance") { tfx_attribute_node_t n; tfx__assign_node_data(&n, values); AddGraphNode(&effect->library->paths[CreateEmitterPathAttributes(effect, false)].distance, &n); }
	}
}

void tfx__assign_node_data(tfx_attribute_node_t *n, tfx_vector_t<tfx_str256_t> *values) {
	n->frame = (float)atof((*values)[1].c_str());
	n->value = (float)atof((*values)[2].c_str());
	n->flags = (bool)atoi((*values)[3].c_str()) ? tfxAttributeNodeFlags_is_curve : 0;
	n->left.x = (float)atof((*values)[4].c_str());
	n->left.y = (float)atof((*values)[5].c_str());
	n->right.x = (float)atof((*values)[6].c_str());
	n->right.y = (float)atof((*values)[7].c_str());
	if (n->flags & tfxAttributeNodeFlags_is_curve)
		n->flags |= tfxAttributeNodeFlags_curves_initialised;
}

void tfx__assign_stage_property_u32(tfx_effect_emitter_t *effect, tfx_str_t *field, tfxU32 value) {
}

void tfx__assign_stage_property_float(tfx_effect_emitter_t *effect, tfx_str_t *field, float value) {
}

void tfx__assign_stage_property_bool(tfx_effect_emitter_t *effect, tfx_str_t *field, bool value) {
}

void tfx__assign_stage_property_int(tfx_effect_emitter_t *effect, tfx_str_t *field, int value) {
}

void tfx__assign_stage_property_str(tfx_effect_emitter_t *effect, tfx_str_t *field, tfx_str_t *value) {
	if (*field == "name") {
		GetEffectInfo(effect)->name = *value;
	}
}

void tfx__assign_frame_meta_property_u32(tfx_frame_meta_t *metrics, tfx_str_t *field, tfxU32 value, tfxU32 file_version) {
	if (*field == "total_sprites")
		metrics->total_sprites = value;
}

tfx_vec3_t tfx__str_to_vec3(tfx_vector_t<tfx_str256_t> *str) {
	TFX_ASSERT(str->size() == 3);    //array must be size 3
	return tfx_vec3_t((float)atof((*str)[0].c_str()), (float)atof((*str)[1].c_str()), (float)atof((*str)[2].c_str()));
}

tfx_vec2_t tfx__str_to_vec2(tfx_vector_t<tfx_str256_t> *str) {
	TFX_ASSERT(str->size() == 2);    //array must be size 2
	return tfx_vec2_t((float)atof((*str)[0].c_str()), (float)atof((*str)[1].c_str()));
}

void tfx__assign_frame_meta_property_vec3(tfx_frame_meta_t *metrics, tfx_str_t *field, tfx_vec3_t value, tfxU32 file_version) {
	if (*field == "min_corner") {
		metrics->min_corner = value;
	}
	if (*field == "max_corner") {
		metrics->max_corner = value;
		//Max corner should always be read after min_corner so can put this here.
		tfx_vec3_t half_extents = (metrics->max_corner - metrics->min_corner) * 0.5f;
		metrics->radius = tfx__length_vec3(&half_extents);
		metrics->bb_center_point = (metrics->max_corner + metrics->min_corner) * 0.5f;
	}
}

void tfx__assign_animation_emitter_property_vec2(tfx_animation_emitter_properties_t *properties, tfx_str_t *field, tfx_vec2_t value, tfxU32 file_version) {
	if (*field == "handle") {
		properties->handle = value;
	}
}

void tfx__assign_animation_emitter_property_float(tfx_animation_emitter_properties_t *properties, tfx_str_t *field, float value, tfxU32 file_version) {
	if (*field == "animation_frames")
		properties->animation_frames = value;
}

void tfx__assign_animation_emitter_property_u32(tfx_animation_emitter_properties_t *properties, tfx_str_t *field, tfxU32 value, tfxU32 file_version) {
	if (*field == "flags") properties->flags = value;
	if (*field == "start_frame_index") properties->start_frame_index = value;
}

void tfx__assign_sprite_data_metrics_property_u32(tfx_sprite_data_metrics_t *metrics, tfx_str_t *field, tfxU32 value, tfxU32 file_version) {
	if (*field == "start_offset")
		metrics->start_offset = value;
	if (*field == "frames_after_compression")
		metrics->frames_after_compression = value;
	if (*field == "real_frames")
		metrics->real_frames = value;
	if (*field == "frame_count")
		metrics->frame_count = value;
	if (*field == "total_sprites")
		metrics->total_sprites = value;
	if (*field == "total_memory_for_sprites")
		metrics->total_memory_for_sprites = value;
	if (*field == "flags")
		metrics->flags = value;
	if (*field == "animation_flags")
		metrics->animation_flags = value;
}

void tfx__assign_property_line(tfx_effect_emitter_t *effect, tfx_vector_t<tfx_str256_t> *pair, tfxU32 file_version) {
	switch (tfxStore->data_types.names_and_types.At((*pair)[0])) {
	case tfxUInt64:
		tfx__assign_effector_property_u64(effect, &(*pair)[0], (tfxU64)strtoull((*pair)[1].c_str(), NULL, 10), file_version);
		break;
	case tfxUint:
		tfx__assign_effector_property_u32(effect, &(*pair)[0], (tfxU32)atoi((*pair)[1].c_str()), file_version);
		break;
	case tfxFloat:
		tfx__assign_effector_property(effect, &(*pair)[0], (float)atof((*pair)[1].c_str()));
		break;
	case tfxSInt:
		tfx__assign_effector_property_int(effect, &(*pair)[0], atoi((*pair)[1].c_str()));
		break;
	case tfxBool:
		tfx__assign_effector_property_bool(effect, &(*pair)[0], (bool)(atoi((*pair)[1].c_str())));
		break;
	case tfxString:
		tfx__assign_effector_property_str(effect, &(*pair)[0], (*pair)[1]);
		break;
	default:
		break;
	}
}

void tfx__assign_sprite_data_metrics_property_u64(tfx_sprite_data_metrics_t *metrics, tfx_str_t *field, tfxU64 value, tfxU32 file_version) {
	if (*field == "path_hash") metrics->path_hash = value;
}

void tfx__assign_sprite_data_metrics_property_float(tfx_sprite_data_metrics_t *metrics, tfx_str_t *field, float value, tfxU32 file_version) {
	if (*field == "animation_length_in_time") metrics->animation_length_in_time = value;
}

void tfx__assign_sprite_data_metrics_property_str(tfx_sprite_data_metrics_t *metrics, tfx_str_t *field, tfx_str_t value, tfxU32 file_version) {
	if (*field == "name") metrics->name = value;
}

void tfx__assign_effector_property_u64(tfx_effect_emitter_t *effect, tfx_str_t *field, tfxU64 value, tfxU32 file_version) {
	if (*field == "image_hash") GetEffectProperties(effect)->image_hash = value;
	if (*field == "paired_emitter_hash") GetEffectProperties(effect)->paired_emitter_hash = value;
}

void tfx__assign_effector_property_u32(tfx_effect_emitter_t *effect, tfx_str_t *field, tfxU32 value, tfxU32 file_version) {
	tfx_emitter_properties_t *emitter_properties = GetEffectProperties(effect);
	if (*field == "image_index") emitter_properties->image_index = value;
	if (*field == "spawn_amount") emitter_properties->spawn_amount = value;
	if (*field == "spawn_amount_variation") emitter_properties->spawn_amount_variation = value;
	if (*field == "frames") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].frames = value;
	if (*field == "current_frame") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].current_frame = value;
	if (*field == "seed") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].seed = value;
	if (*field == "layer") emitter_properties->layer = value >= tfxLAYERS ? value = tfxLAYERS - 1 : value;
	if (*field == "frame_offset") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].frame_offset = value;
	if (*field == "single_shot_limit") emitter_properties->single_shot_limit = value;
	if (*field == "billboard_option") {
		//billboard options were changed so I added this to at least update the align to camera and vector values.
		//0 and 1 should still be ok, 4 now maps to 2, and 2 should now be 3 but I'll just manually update the effect
		//libs for that. This can be removed at some point.
		if (file_version == 1) {
			if (value == 4)
				value = 2;
			else if (value == 2) {
				value = 3;
			}
		}
		emitter_properties->billboard_option = (tfx_billboarding_option)value;
	}
	if (*field == "vector_align_type") emitter_properties->vector_align_type = value >= 0 && value < tfxVectorAlignType_max ? (tfx_vector_align_type)value : (tfx_vector_align_type)0;
	if (*field == "angle_setting") emitter_properties->angle_settings = (tfxAngleSettingFlags)value;
	if (*field == "sort_passes") effect->sort_passes = tfxMin(5, value);
	if (*field == "animation_flags") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].animation_flags = value | tfxAnimationFlags_needs_recording;
	if (*field == "sprite_data_flags") effect->library->sprite_data_settings[GetEffectInfo(effect)->sprite_data_settings_index].animation_flags = value | tfxAnimationFlags_needs_recording;
	if (*field == "sprite_data_seed") effect->library->sprite_data_settings[GetEffectInfo(effect)->sprite_data_settings_index].seed = value;
	if (*field == "sprite_data_frame_offset") effect->library->sprite_data_settings[GetEffectInfo(effect)->sprite_data_settings_index].frame_offset = value;
	if (*field == "sprite_data_frames") effect->library->sprite_data_settings[GetEffectInfo(effect)->sprite_data_settings_index].real_frames = value;
	if (*field == "sprite_data_extra_frames_count") effect->library->sprite_data_settings[GetEffectInfo(effect)->sprite_data_settings_index].extra_frames_count = value;
	if (*field == "maximum_active_paths") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; path->maximum_active_paths = value;
	}
	if (*field == "maximum_path_cycles") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; path->maximum_paths = value;
	}
	if (*field == "path_node_count") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; path->node_count = value;
	}
}
void tfx__assign_effector_property_int(tfx_effect_emitter_t *effect, tfx_str_t *field, int value) {
	tfx_emitter_properties_t *emitter_properties = GetEffectProperties(effect);
	if (*field == "emission_type") emitter_properties->emission_type = (tfx_emission_type)value;
	if (*field == "emission_direction") emitter_properties->emission_direction = (tfx_emission_direction)value;
	if (*field == "color_option") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].color_option = value > 3 ? tfxFullColor : (tfx_export_color_options)value;
	if (*field == "export_option") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].export_option = (tfx_export_options)value;
	if (*field == "end_behaviour") emitter_properties->end_behaviour = (tfx_line_traversal_end_behaviour)value;
	if (*field == "frame_offset") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].frame_offset = value;
	if (*field == "extra_frames_count") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].extra_frames_count = value;
	if (*field == "path_extrusion_type") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)];  path->extrusion_type = (tfx_path_extrusion_type)value;
	}
	if (*field == "path_generator_type") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)];  path->generator_type = (tfx_path_generator_type)value;
	}
}
void tfx__assign_effector_property_str(tfx_effect_emitter_t *effect, tfx_str_t *field, tfx_str_t &value) {
	if (*field == "name") {
		GetEffectInfo(effect)->name = value;
	}
}
void tfx__assign_effector_property(tfx_effect_emitter_t *effect, tfx_str_t *field, float value) {
	tfx_emitter_properties_t *emitter_properties = GetEffectProperties(effect);
	if (*field == "position_x") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].position.x = value;
	if (*field == "position_y") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].position.y = value;
	if (*field == "position_z") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].position.z = value;
	if (*field == "frame_width") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].frame_size.x = value;
	if (*field == "frame_height") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].frame_size.y = value;
	if (*field == "zoom") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].zoom = value;
	if (*field == "scale") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].scale = value;
	if (*field == "playback_speed") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].playback_speed = value;
	if (*field == "camera_position_x") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings.camera_position.x = value;
	if (*field == "camera_position_y") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings.camera_position.y = value;
	if (*field == "camera_position_z") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings.camera_position.z = value;
	if (*field == "camera_pitch") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings.camera_pitch = value;
	if (*field == "camera_yaw") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings.camera_yaw = value;
	if (*field == "camera_fov") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings.camera_fov = value;
	if (*field == "camera_floor_height") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings.camera_floor_height = value;
	if (*field == "camera_isometric_scale") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings.camera_isometric_scale = value;
	if (*field == "orthographic_camera_position_x") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_position.x = value;
	if (*field == "orthographic_camera_position_y") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_position.y = value;
	if (*field == "orthographic_camera_position_z") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_position.z = value;
	if (*field == "orthographic_camera_pitch") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_pitch = value;
	if (*field == "orthographic_camera_yaw") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_yaw = value;
	if (*field == "orthographic_camera_fov") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_fov = value;
	if (*field == "orthographic_camera_floor_height") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_floor_height = value;
	if (*field == "orthographic_camera_isometric_scale") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_isometric_scale = value;
	if (*field == "preview_camera_position_x") effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings].camera_settings.camera_position.x = value;
	if (*field == "preview_camera_position_y") effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings].camera_settings.camera_position.y = value;
	if (*field == "preview_camera_position_z") effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings].camera_settings.camera_position.z = value;
	if (*field == "preview_camera_pitch") effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings].camera_settings.camera_pitch = value;
	if (*field == "preview_camera_yaw") effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings].camera_settings.camera_yaw = value;
	if (*field == "preview_camera_fov") effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings].camera_settings.camera_fov = value;
	if (*field == "preview_camera_floor_height") effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings].camera_settings.camera_floor_height = value;
	if (*field == "preview_camera_isometric_scale") effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings].camera_settings.camera_isometric_scale = value == 0 ? 5.f : value;
	if (*field == "preview_effect_z_offset") effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings].effect_z_offset = value;
	if (*field == "preview_camera_speed") effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings].camera_speed = value;
	if (*field == "image_handle_x") emitter_properties->image_handle.x = value;
	if (*field == "image_handle_y") emitter_properties->image_handle.y = value;
	if (*field == "delay_spawning") emitter_properties->delay_spawning = value;
	if (*field == "grid_rows") emitter_properties->grid_points.x = value;
	if (*field == "grid_columns") emitter_properties->grid_points.y = value;
	if (*field == "grid_depth") emitter_properties->grid_points.z = value;
	if (*field == "loop_length") emitter_properties->loop_length = value < 0 ? 0.f : value;
	if (*field == "noise_base_offset_range") emitter_properties->noise_base_offset_range = value < 0 ? 0.f : value;
	if (*field == "emitter_handle_x") emitter_properties->emitter_handle.x = value;
	if (*field == "emitter_handle_y") emitter_properties->emitter_handle.y = value;
	if (*field == "emitter_handle_z") emitter_properties->emitter_handle.z = value;
	if (*field == "image_start_frame") emitter_properties->start_frame = value;
	if (*field == "image_end_frame") emitter_properties->end_frame = value;
	if (*field == "image_frame_rate") emitter_properties->frame_rate = value;
	if (*field == "angle_offset") emitter_properties->angle_offsets.roll = value;
	if (*field == "angle_offset_pitch") emitter_properties->angle_offsets.pitch = value;
	if (*field == "angle_offset_yaw") emitter_properties->angle_offsets.yaw = value;
	if (*field == "sprite_data_playback_speed") effect->library->sprite_data_settings[GetEffectInfo(effect)->sprite_data_settings_index].playback_speed = value;
	if (*field == "sprite_data_recording_frame_rate") effect->library->sprite_data_settings[GetEffectInfo(effect)->sprite_data_settings_index].recording_frame_rate = value;
	if (*field == "path_rotation_range") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; path->rotation_range = value;
	}
	if (*field == "path_rotation_pitch") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; path->rotation_pitch = value;
	}
	if (*field == "path_rotation_yaw") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; path->rotation_yaw = value;
	}
	if (*field == "path_rotation_stagger") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; path->rotation_stagger = value;
	}
	if (*field == "path_handle_x") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)];  path->offset.x = value;
	}
	if (*field == "path_handle_y") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; path->offset.y = value;
	}
	if (*field == "path_handle_z") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; path->offset.z = value;
	}
}
void tfx__assign_effector_property_bool(tfx_effect_emitter_t *effect, tfx_str_t *field, bool value) {
	if (*field == "loop") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].animation_flags |= value ? tfxAnimationFlags_loop : 0;
	if (*field == "seamless") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].animation_flags |= value ? tfxAnimationFlags_seamless : 0;
	if (*field == "export_with_transparency") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].animation_flags |= value ? tfxAnimationFlags_export_with_transparency : 0;
	if (*field == "camera_isometric") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings.camera_isometric = false;
	if (*field == "camera_hide_floor") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings.camera_hide_floor = value;
	if (*field == "orthographic_camera_isometric") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_isometric = true;
	if (*field == "orthographic_camera_hide_floor") effect->library->sprite_sheet_settings[GetEffectInfo(effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_hide_floor = value;
	if (*field == "preview_attach_effect_to_camera") effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings].attach_effect_to_camera = value;
	if (*field == "preview_camera_hide_floor") effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings].camera_settings.camera_hide_floor = value;
	if (*field == "preview_camera_isometric") effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings].camera_settings.camera_isometric = value;
	if (*field == "random_color") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_random_color; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_random_color; }
	}
	if (*field == "relative_position") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_relative_position; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_relative_position; }
	}
	if (*field == "relative_angle") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_relative_angle; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_relative_angle; }
	}
	if (*field == "image_handle_auto_center") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_image_handle_auto_center; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_image_handle_auto_center; }
	}
	if (*field == "single") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_single; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_single; }
	}
	if (*field == "wrap_single_sprite") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_wrap_single_sprite; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_wrap_single_sprite; }
	}
	if (*field == "spawn_on_grid") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_spawn_on_grid; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_spawn_on_grid; }
	}
	if (*field == "grid_spawn_clockwise") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_grid_spawn_clockwise; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_grid_spawn_clockwise; }
	}
	if (*field == "fill_area") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_fill_area; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_fill_area; }
	}
	if (*field == "grid_spawn_random") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_grid_spawn_random; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_grid_spawn_random; }
	}
	if (*field == "area_open_ends") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_area_open_ends; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_area_open_ends; }
	}
	if (*field == "emitter_handle_auto_center") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_emitter_handle_auto_center; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_emitter_handle_auto_center; }
	}
	if (*field == "edge_traversal") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_edge_traversal; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_edge_traversal; }
	}
	if (*field == "image_reverse_animation") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_reverse_animation; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_reverse_animation; }
	}
	if (*field == "image_play_once") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_play_once; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_play_once; }
	}
	if (*field == "image_animate") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_animate; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_animate; }
	}
	if (*field == "image_random_start_frame") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_random_start_frame; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_random_start_frame; }
	}
	if (*field == "global_uniform_size") {
		if (value) { effect->effect_flags |= tfxEffectPropertyFlags_global_uniform_size; }
		else { effect->property_flags &= ~tfxEffectPropertyFlags_global_uniform_size; }
	}
	if (*field == "base_uniform_size") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_base_uniform_size; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_base_uniform_size; }
	}
	if (*field == "lifetime_uniform_size") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_lifetime_uniform_size; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_lifetime_uniform_size; }
	}
	if (*field == "use_spawn_ratio") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_use_spawn_ratio; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_use_spawn_ratio; }
	}
	if (*field == "is_3d") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_effect_is_3d; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_effect_is_3d; }
	}
	if (*field == "draw_order_by_age") {
		if (value) { effect->effect_flags |= tfxEffectPropertyFlags_age_order; }
		else { effect->effect_flags &= ~tfxEffectPropertyFlags_age_order; }
	}
	if (*field == "draw_order_by_depth") {
		if (value) { effect->effect_flags |= tfxEffectPropertyFlags_depth_draw_order; }
		else { effect->effect_flags &= ~tfxEffectPropertyFlags_depth_draw_order; }
	}
	if (*field == "guaranteed_draw_order") {
		if (value) { effect->effect_flags |= tfxEffectPropertyFlags_guaranteed_order; }
		else { effect->effect_flags &= ~tfxEffectPropertyFlags_guaranteed_order; }
	}
	if (*field == "include_in_sprite_data_export") {
		if (value) { effect->effect_flags |= tfxEffectPropertyFlags_include_in_sprite_data_export; }
		else { effect->effect_flags &= ~tfxEffectPropertyFlags_include_in_sprite_data_export; }
	}
	if (*field == "use_path_for_direction") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_use_path_for_direction; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_use_path_for_direction; }
	}
	if (*field == "alt_velocity_lifetime_sampling") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_alt_velocity_lifetime_sampling; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_alt_velocity_lifetime_sampling; }
	}
	if (*field == "alt_color_lifetime_sampling") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_alt_color_lifetime_sampling; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_alt_color_lifetime_sampling; }
	}
	if (*field == "alt_size_lifetime_sampling") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_alt_size_lifetime_sampling; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_alt_size_lifetime_sampling; }
	}
	if (*field == "use_simple_motion_randomness") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_use_simple_motion_randomness; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_use_simple_motion_randomness; }
	}
	if (*field == "spawn_location_source") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_spawn_location_source; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_spawn_location_source; }
	}
	if (*field == "use_color_hint") {
		if (value) { effect->property_flags |= tfxEmitterPropertyFlags_use_color_hint; }
		else { effect->property_flags &= ~tfxEmitterPropertyFlags_use_color_hint; }
	}
	//if (*field == "simple_motion_smoothstep") {
		//if (value) { effect->property_flags |= tfxEmitterPropertyFlags_simple_motion_smoothstep; } else { effect->property_flags &= ~tfxEmitterPropertyFlags_simple_motion_smoothstep; }
	//}
	if (*field == "path_is_2d") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; if (value) { path->flags |= tfxPathFlags_2d; }
	}
	if (*field == "path_mode_origin") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; if (value) { path->flags |= tfxPathFlags_mode_origin; }
	}
	if (*field == "path_mode_node") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; if (value) { path->flags |= tfxPathFlags_mode_node; }
	}
	if (*field == "path_space_nodes_evenly") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; if (value) { path->flags |= tfxPathFlags_space_nodes_evenly; }
	}
	if (*field == "path_rotation_range_yaw_only") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; if (value) { path->flags |= tfxPathFlags_rotation_range_yaw_only; }
	}
	if (*field == "path_reverse_direction") {
		tfx_emitter_path_t *path = &effect->library->paths[CreateEmitterPathAttributes(effect, false)]; if (value) { path->flags |= tfxPathFlags_reverse_direction; }
	}

}

void tfx__stream_emitter_properties(tfx_emitter_properties_t *property, tfxEmitterPropertyFlags flags, tfx_str_t *file) {
	file->AddLine("image_hash=%llu", property->image_hash);
	file->AddLine("image_handle_x=%f", property->image_handle.x);
	file->AddLine("image_handle_y=%f", property->image_handle.y);
	file->AddLine("image_start_frame=%f", property->start_frame);
	file->AddLine("image_end_frame=%f", property->end_frame);
	file->AddLine("image_frame_rate=%f", property->frame_rate);
	file->AddLine("image_play_once=%i", (flags & tfxEmitterPropertyFlags_play_once));
	file->AddLine("image_reverse_animation=%i", (flags & tfxEmitterPropertyFlags_reverse_animation));
	file->AddLine("image_animate=%i", (flags & tfxEmitterPropertyFlags_animate));
	file->AddLine("image_random_start_frame=%i", (flags & tfxEmitterPropertyFlags_random_start_frame));
	file->AddLine("image_handle_auto_center=%i", (flags & tfxEmitterPropertyFlags_image_handle_auto_center));
	file->AddLine("paired_emitter_hash=%llu", property->paired_emitter_hash);
	file->AddLine("spawn_amount=%i", property->spawn_amount);
	file->AddLine("spawn_amount_variation=%i", property->spawn_amount_variation);
	file->AddLine("emission_type=%i", property->emission_type);
	file->AddLine("emission_direction=%i", property->emission_direction);
	file->AddLine("grid_rows=%f", property->grid_points.x);
	file->AddLine("grid_columns=%f", property->grid_points.y);
	file->AddLine("grid_depth=%f", property->grid_points.z);
	file->AddLine("delay_spawning=%f", property->delay_spawning);
	file->AddLine("loop_length=%f", property->loop_length);
	file->AddLine("emitter_handle_x=%f", property->emitter_handle.x);
	file->AddLine("emitter_handle_y=%f", property->emitter_handle.y);
	file->AddLine("emitter_handle_z=%f", property->emitter_handle.z);
	file->AddLine("end_behaviour=%i", property->end_behaviour);
	file->AddLine("random_color=%i", (flags & tfxEmitterPropertyFlags_random_color));
	file->AddLine("relative_position=%i", (flags & tfxEmitterPropertyFlags_relative_position));
	file->AddLine("relative_angle=%i", (flags & tfxEmitterPropertyFlags_relative_angle));
	file->AddLine("single=%i", (flags & tfxEmitterPropertyFlags_single));
	file->AddLine("wrap_single_sprite=%i", (flags & tfxEmitterPropertyFlags_wrap_single_sprite));
	file->AddLine("single_shot_limit=%i", property->single_shot_limit);
	file->AddLine("spawn_on_grid=%i", (flags & tfxEmitterPropertyFlags_spawn_on_grid));
	file->AddLine("grid_spawn_clockwise=%i", (flags & tfxEmitterPropertyFlags_grid_spawn_clockwise));
	file->AddLine("fill_area=%i", (flags & tfxEmitterPropertyFlags_fill_area));
	file->AddLine("grid_spawn_random=%i", (flags & tfxEmitterPropertyFlags_grid_spawn_random));
	file->AddLine("area_open_ends=%i", (flags & tfxEmitterPropertyFlags_area_open_ends));
	file->AddLine("emitter_handle_auto_center=%i", (flags & tfxEmitterPropertyFlags_emitter_handle_auto_center));
	file->AddLine("edge_traversal=%i", (flags & tfxEmitterPropertyFlags_edge_traversal));
	file->AddLine("angle_setting=%i", property->angle_settings);
	file->AddLine("angle_offset=%f", property->angle_offsets.roll);
	file->AddLine("angle_offset_pitch=%f", property->angle_offsets.pitch);
	file->AddLine("angle_offset_yaw=%f", property->angle_offsets.yaw);
	file->AddLine("base_uniform_size=%i", (flags & tfxEmitterPropertyFlags_base_uniform_size));
	file->AddLine("lifetime_uniform_size=%i", (flags & tfxEmitterPropertyFlags_lifetime_uniform_size));
	file->AddLine("use_spawn_ratio=%i", (flags & tfxEmitterPropertyFlags_use_spawn_ratio));
	file->AddLine("billboard_option=%i", property->billboard_option);
	file->AddLine("vector_align_type=%i", property->vector_align_type);
	file->AddLine("layer=%i", property->layer);
	file->AddLine("use_path_for_direction=%i", (flags & tfxEmitterPropertyFlags_use_path_for_direction));
	file->AddLine("alt_velocity_lifetime_sampling=%i", (flags & tfxEmitterPropertyFlags_alt_velocity_lifetime_sampling));
	file->AddLine("alt_color_lifetime_sampling=%i", (flags & tfxEmitterPropertyFlags_alt_color_lifetime_sampling));
	file->AddLine("alt_size_lifetime_sampling=%i", (flags & tfxEmitterPropertyFlags_alt_size_lifetime_sampling));
	file->AddLine("use_simple_motion_randomness=%i", (flags & tfxEmitterPropertyFlags_use_simple_motion_randomness));
	file->AddLine("spawn_location_source=%i", (flags & tfxEmitterPropertyFlags_spawn_location_source));
	file->AddLine("use_color_hint=%i", (flags & tfxEmitterPropertyFlags_use_color_hint));
	//file->AddLine("simple_motion_smoothstep=%i", (flags & tfxEmitterPropertyFlags_simple_motion_smoothstep));
}

void tfx__stream_effect_properties(tfx_effect_emitter_t *effect, tfx_str_t *file) {
	tfx_emitter_properties_t *properties = GetEffectProperties(effect);
	file->AddLine("is_3d=%i", (effect->property_flags & tfxEmitterPropertyFlags_effect_is_3d));
	file->AddLine("draw_order_by_age=%i", effect->effect_flags & tfxEffectPropertyFlags_age_order);
	file->AddLine("draw_order_by_depth=%i", effect->effect_flags & tfxEffectPropertyFlags_depth_draw_order);
	file->AddLine("guaranteed_draw_order=%i", effect->effect_flags & tfxEffectPropertyFlags_guaranteed_order);
	file->AddLine("include_in_sprite_data_export=%i", effect->effect_flags & tfxEffectPropertyFlags_include_in_sprite_data_export);
	file->AddLine("sort_passes=%i", effect->sort_passes);
	file->AddLine("noise_base_offset_range=%f", properties->noise_base_offset_range);
	file->AddLine("loop_length=%f", properties->loop_length);
	file->AddLine("emitter_handle_x=%f", properties->emitter_handle.x);
	file->AddLine("emitter_handle_y=%f", properties->emitter_handle.y);
	file->AddLine("emitter_handle_z=%f", properties->emitter_handle.z);
	file->AddLine("global_uniform_size=%i", (effect->effect_flags & tfxEffectPropertyFlags_global_uniform_size));
}

void tfx__stream_path_properties(tfx_effect_emitter_t *effect, tfx_str_t *file) {
	if (effect->path_attributes != tfxINVALID) {
		tfx_emitter_path_t *path = &effect->library->paths[effect->path_attributes];
		file->AddLine("path_is_2d=%i", (path->flags & tfxPathFlags_2d));
		file->AddLine("path_mode_origin=%i", (path->flags & tfxPathFlags_mode_origin));
		file->AddLine("path_mode_node=%i", (path->flags & tfxPathFlags_mode_node));
		file->AddLine("path_space_nodes_evenly=%i", (path->flags & tfxPathFlags_space_nodes_evenly));
		file->AddLine("path_generator_type=%i", (path->generator_type));
		file->AddLine("path_extrusion_type=%i", (path->extrusion_type));
		file->AddLine("path_rotation_range=%f", (path->rotation_range));
		file->AddLine("path_rotation_pitch=%f", (path->rotation_pitch));
		file->AddLine("path_rotation_yaw=%f", (path->rotation_yaw));
		file->AddLine("maximum_active_paths=%i", (path->maximum_active_paths));
		file->AddLine("maximum_path_cycles=%i", (path->maximum_paths));
		file->AddLine("path_node_count=%i", (path->node_count));
		file->AddLine("path_rotation_stagger=%f", (path->rotation_stagger));
		file->AddLine("path_rotation_range_yaw_only=%i", (path->flags & tfxPathFlags_rotation_range_yaw_only));
		file->AddLine("path_reverse_direction=%i", (path->flags & tfxPathFlags_reverse_direction));
		file->AddLine("path_handle_x=%f", (path->offset.x));
		file->AddLine("path_handle_y=%f", (path->offset.y));
		file->AddLine("path_handle_z=%f", (path->offset.z));
	}
}

void tfx__stream_graph(const char *name, tfx_graph_t *graph, tfx_str_t *file) {

	for (tfxBucketLoop(graph->nodes)) {
		file->AddLine("%s,%f,%f,%i,%f,%f,%f,%f", name, graph->nodes[i].frame, graph->nodes[i].value, (graph->nodes[i].flags & tfxAttributeNodeFlags_is_curve), graph->nodes[i].left.x, graph->nodes[i].left.y, graph->nodes[i].right.x, graph->nodes[i].right.y);
	}

}

bool CompareNodes(tfx_attribute_node_t *left, tfx_attribute_node_t *right) {
	return left->frame < right->frame;
}

bool IsOvertimeGraph(tfx_graph_t *graph) {
	return (graph->type >= tfxOvertime_velocity && graph->type <= tfxOvertime_motion_randomness && graph->type != tfxOvertime_velocity_adjuster) || graph->type > tfxEmitterGraphMaxIndex;
}

bool IsFactorGraph(tfx_graph_t *graph) {
	return graph->type >= tfxFactor_life && graph->type <= tfxFactor_intensity;
}

bool IsColorGraph(tfx_graph_t *graph) {
	return graph->type >= tfxOvertime_red && graph->type <= tfxOvertime_blendfactor_hint;
}

bool IsBlendfactorGraph(tfx_graph_t *graph) {
	return graph->type == tfxOvertime_blendfactor || graph->type == tfxOvertime_blendfactor_hint;
}

bool IsGlobalGraph(tfx_graph_t *graph) {
	return graph->type >= tfxGlobal_life && graph->type <= tfxGlobal_splatter;
}

bool IsAngleGraph(tfx_graph_t *graph) {
	return (graph->type == tfxTransform_roll || graph->type == tfxTransform_pitch || graph->type == tfxTransform_yaw || graph->type == tfxProperty_emission_pitch || graph->type == tfxProperty_emission_yaw
		|| graph->type == tfxProperty_emission_range || graph->type == tfxProperty_arc_offset || graph->type == tfxProperty_arc_size || graph->type == tfxBase_roll_spin || graph->type == tfxBase_pitch_spin || graph->type == tfxBase_yaw_spin
		|| graph->type == tfxVariation_roll_spin || graph->type == tfxVariation_pitch_spin || graph->type == tfxVariation_yaw_spin || graph->type == tfxOvertime_direction);
}

bool IsTranslationGraph(tfx_graph_t *graph) {
	return graph->type == tfxTransform_translate_x || graph->type == tfxTransform_translate_y || graph->type == tfxTransform_translate_z;
}

void MultiplyAllGraphValues(tfx_graph_t *graph, float scalar) {
	for (tfxBucketLoop(graph->nodes)) {
		graph->nodes[i].value *= scalar;
		graph->nodes[i].left.y *= scalar;
		graph->nodes[i].right.y *= scalar;
	}
}

void CopyGraphNoLookups(tfx_graph_t *src_graph, tfx_graph_t *dst_graph) {
	dst_graph->min = src_graph->min;
	dst_graph->max = src_graph->max;
	dst_graph->graph_preset = src_graph->graph_preset;
	dst_graph->type = src_graph->type;
	dst_graph->effector = src_graph->effector;
	tfxCopyBucketArray<tfx_attribute_node_t>(&dst_graph->nodes, &src_graph->nodes);
	dst_graph->index = src_graph->index;
	dst_graph->lookup.life = src_graph->lookup.life;
}

bool SetNode(tfx_graph_t *graph, tfx_attribute_node_t *node, float _frame, float _value, tfxAttributeNodeFlags flags, float _c0x, float _c0y, float _c1x, float _c1y) {
	node->frame = _frame;
	node->value = _value;
	node->flags = flags;
	node->left.x = _c0x;
	node->left.y = _c0y;
	node->right.x = _c1x;
	node->right.y = _c1y;
	if (&graph->nodes[0] == node) {
		node->frame = graph->min.x;
		ClampNode(graph, node);
	}
	else {
		ClampNode(graph, node);
	}

	if (SortGraph(graph)) {
		ReIndexGraph(graph);
		return true;
	}

	return false;
}

bool SetNodeFrame(tfx_graph_t *graph, tfx_attribute_node_t *node, float *frame) {
	node->frame = *frame;
	if (&graph->nodes[0] == node) {
		node->frame = graph->min.x;
		ClampNode(graph, node);
	}
	else {
		ClampNode(graph, node);
	}
	*frame = node->frame;

	if (SortGraph(graph)) {
		ReIndexGraph(graph);
		return true;
	}

	return false;

}

bool IsNodeCurve(tfx_attribute_node_t *node) {
	return node->flags & tfxAttributeNodeFlags_is_curve;
}

bool NodeCurvesAreInitialised(tfx_attribute_node_t *node) {
	return node->flags & tfxAttributeNodeFlags_curves_initialised;
}

bool SetNodeCurveInitialised(tfx_attribute_node_t *node) {
	return node->flags |= tfxAttributeNodeFlags_curves_initialised;
}

bool SetNodeValue(tfx_graph_t *graph, tfx_attribute_node_t *node, float *value) {
	node->value = *value;
	ClampNode(graph, node);
	*value = node->value;

	return false;
}

bool SetNode(tfx_graph_t *graph, tfx_attribute_node_t *node, float *frame, float *value) {
	float old_frame = node->frame;
	float old_value = node->value;

	node->frame = *frame;
	node->value = *value;

	if (&graph->nodes[0] == node) {
		node->frame = graph->min.x;
		ClampNode(graph, node);
	}
	else {
		ClampNode(graph, node);
	}

	if (node->flags & tfxAttributeNodeFlags_curves_initialised) {
		node->left.y += node->value - old_value;
		node->left.x += node->frame - old_frame;
		node->right.y += node->value - old_value;
		node->right.x += node->frame - old_frame;
		ClampCurve(graph, &node->right, node);
		ClampCurve(graph, &node->left, node);
	}

	*frame = node->frame;
	*value = node->value;

	if (SortGraph(graph)) {
		ReIndexGraph(graph);
		return true;
	}

	return false;
}

void SetCurve(tfx_graph_t *graph, tfx_attribute_node_t *node, bool is_left_curve, float *frame, float *value) {
	if (is_left_curve) {
		node->left.x = *frame;
		node->left.y = *value;
		if (node->left.x > node->frame)
			node->left.x = node->frame;
		ClampCurve(graph, &node->left, node);
		*frame = node->left.x;
		*value = node->left.y;
	}
	else {
		node->right.x = *frame;
		node->right.y = *value;
		if (node->right.x < node->frame)
			node->right.x = node->frame;
		ClampCurve(graph, &node->right, node);
		*frame = node->right.x;
		*value = node->right.y;
	}
}

bool MoveNode(tfx_graph_t *graph, tfx_attribute_node_t *node, float frame, float value, bool sort) {
	float old_frame = node->frame;
	float old_value = node->value;

	node->frame += frame;
	node->value += value;

	if (&graph->nodes[0] == node) {
		node->frame = graph->min.x;
		ClampNode(graph, node);
	}
	else {
		ClampNode(graph, node);
	}

	if (node->flags & tfxAttributeNodeFlags_curves_initialised) {
		node->left.y += node->value - old_value;
		node->left.x += node->frame - old_frame;
		node->right.y += node->value - old_value;
		node->right.x += node->frame - old_frame;
		ClampCurve(graph, &node->right, node);
		ClampCurve(graph, &node->left, node);
	}

	if (sort) {
		if (SortGraph(graph)) {
			ReIndexGraph(graph);
			return true;
		}
	}

	return false;
}

void ClampNode(tfx_graph_t *graph, tfx_attribute_node_t *node) {
	if (node->value < graph->min.y) node->value = graph->min.y;
	if (node->frame < graph->min.x) node->frame = graph->min.x;
	if (node->value > graph->max.y) node->value = graph->max.y;
	if (node->frame > graph->max.x) node->frame = graph->max.x;
}

void ClampGraph(tfx_graph_t *graph) {
	for (tfxBucketLoop(graph->nodes)) {
		ClampNode(graph, &graph->nodes[i]);
		if (graph->nodes[i].flags & tfxAttributeNodeFlags_is_curve) {
			ClampCurve(graph, &graph->nodes[i].left, &graph->nodes[i]);
			ClampCurve(graph, &graph->nodes[i].right, &graph->nodes[i]);
		}
	}
}

void ClampCurve(tfx_graph_t *graph, tfx_vec2_t *p, tfx_attribute_node_t *node) {
	if (p->y < graph->min.y) p->y = graph->min.y;
	if (p->x < graph->min.x) p->x = graph->min.x;
	if (p->y > graph->max.y) p->y = graph->max.y;
	if (p->x > graph->max.x) p->x = graph->max.x;

	tfx_attribute_node_t *next = GetGraphNextNode(graph, node);
	if (next) {
		if (p->x > next->frame) p->x = next->frame;
	}

	tfx_attribute_node_t *prev = GetGraphPrevNode(graph, node);
	if (prev) {
		if (p->x < prev->frame) p->x = prev->frame;
	}
}

tfx_graph_t::tfx_graph_t() {
	graph_preset = tfxGlobalPercentPreset;
	index = 0;
	gamma = 0.f;
	type = tfxEmitterGraphMaxIndex;
	min.x = 0.f;
	min.y = 0.f;
	max.x = 1000.f;
	max.y = 1000.f;
	effector = nullptr;
	nodes = tfxCreateBucketArray<tfx_attribute_node_t>(8);
}

tfx_graph_t::tfx_graph_t(tfxU32 bucket_size) {
	type = tfxEmitterGraphMaxIndex;
	graph_preset = tfxGlobalPercentPreset;
	index = 0;
	gamma = 0.f;
	min.x = 0.f;
	min.y = 0.f;
	max.x = 1000.f;
	max.y = 1000.f;
	effector = nullptr;
	nodes = tfxCreateBucketArray<tfx_attribute_node_t>(bucket_size);
}

tfx_graph_t::~tfx_graph_t() {
	//nodes.free_all();
	//Free();
}

const float BEZIER_ACCURACY = 0.01f;

float GetDistance(float fromx, float fromy, float tox, float toy) {

	float w = tox - fromx;
	float h = toy - fromy;

	return sqrtf(w * w + h * h);

}

tfx_vec2_t GetQuadBezierClamp(tfx_vec2_t p0, tfx_vec2_t p1, tfx_vec2_t p2, float t, float ymin, float ymax) {
	tfx_vec2_t b;
	float ti = 1.f - t;
	float ti2 = ti * ti;
	float t2 = t * t;
	b.x = ti2 * p0.x + 2.f * t * ti * p1.x + t2 * p2.x;
	b.y = ti2 * p0.y + 2.f * t * ti * p1.y + t2 * p2.y;
	b.x = tfx__Clamp(p0.x, p2.x, b.x);
	b.y = tfx__Clamp(ymin, ymax, b.y);
	return b;
}

tfx_vec2_t GetQuadBezier(tfx_vec2_t p0, tfx_vec2_t p1, tfx_vec2_t p2, float t) {
	tfx_vec2_t b;
	float ti = 1.f - t;
	float ti2 = ti * ti;
	float t2 = t * t;
	b.x = ti2 * p0.x + 2.f * t * ti * p1.x + t2 * p2.x;
	b.y = ti2 * p0.y + 2.f * t * ti * p1.y + t2 * p2.y;
	return b;
}

tfx_vec2_t GetCubicBezierClamp(tfx_vec2_t p0, tfx_vec2_t p1, tfx_vec2_t p2, tfx_vec2_t p3, float t, float ymin, float ymax) {
	tfx_vec2_t b;
	float ti = 1.f - t;
	float ti3 = ti * ti * ti;
	float ti2 = ti * ti;
	float t3 = t * t * t;
	float t2 = t * t;
	b.x = ti3 * p0.x + 3.f * t * ti2 * p1.x + 3.f * t2 * ti * p2.x + t3 * p3.x;
	b.y = ti3 * p0.y + 3.f * t * ti2 * p1.y + 3.f * t2 * ti * p2.y + t3 * p3.y;
	b.x = tfx__Clamp(p0.x, p2.x, b.x);
	b.y = tfx__Clamp(ymin, ymax, b.y);
	return b;
}

tfx_vec2_t GetCubicBezier(tfx_vec2_t p0, tfx_vec2_t p1, tfx_vec2_t p2, tfx_vec2_t p3, float t) {
	tfx_vec2_t b;
	float ti = 1.f - t;
	float ti3 = ti * ti * ti;
	float ti2 = ti * ti;
	float t3 = t * t * t;
	float t2 = t * t;
	b.x = ti3 * p0.x + 3.f * t * ti2 * p1.x + 3.f * t2 * ti * p2.x + t3 * p3.x;
	b.y = ti3 * p0.y + 3.f * t * ti2 * p1.y + 3.f * t2 * ti * p2.y + t3 * p3.y;
	return b;
}

tfx_vec3_t GetCubicBezier3d(tfx_vec4_t *p0, tfx_vec4_t *p1, tfx_vec4_t *p2, tfx_vec4_t *p3, float t) {
	tfx_vec3_t b;
	float ti = 1.f - t;
	float ti3 = ti * ti * ti;
	float ti2 = ti * ti;
	float t3 = t * t * t;
	float t2 = t * t;
	b.x = ti3 * p0->x + 3.f * t * ti2 * p1->x + 3.f * t2 * ti * p2->x + t3 * p3->x;
	b.y = ti3 * p0->y + 3.f * t * ti2 * p1->y + 3.f * t2 * ti * p2->y + t3 * p3->y;
	b.z = ti3 * p0->z + 3.f * t * ti2 * p1->z + 3.f * t2 * ti * p2->z + t3 * p3->z;
	return b;
}

float GetBezierValue(const tfx_attribute_node_t *lastec, const tfx_attribute_node_t *a, float t, float ymin, float ymax) {
	if (lastec) {
		if (a->flags & tfxAttributeNodeFlags_is_curve) {
			if (lastec->flags & tfxAttributeNodeFlags_is_curve) {
				tfx_vec2_t p0(lastec->frame, lastec->value);
				tfx_vec2_t p1(lastec->right.x, lastec->right.y);
				tfx_vec2_t p2(a->left.x, a->left.y);
				tfx_vec2_t p3(a->frame, a->value);
				tfx_vec2_t value = GetCubicBezierClamp(p0, p1, p2, p3, t, ymin, ymax);
				return value.y;
			}
			else {
				tfx_vec2_t p0(lastec->frame, lastec->value);
				tfx_vec2_t p1(a->left.x, a->left.y);
				tfx_vec2_t p2(a->frame, a->value);
				tfx_vec2_t value = GetQuadBezierClamp(p0, p1, p2, t, ymin, ymax);
				return value.y;
			}
		}
		else if (lastec->flags & tfxAttributeNodeFlags_is_curve) {
			tfx_vec2_t p0(lastec->frame, lastec->value);
			tfx_vec2_t p1(lastec->right.x, lastec->right.y);
			tfx_vec2_t p2(a->frame, a->value);
			tfx_vec2_t value = GetQuadBezierClamp(p0, p1, p2, t, ymin, ymax);
			return value.y;
		}
	}
	else {
		return 0;
	}

	return 0;
}

tfx_attribute_node_t *AddGraphNode(tfx_graph_t *graph, float _frame, float _value, tfxAttributeNodeFlags flags, float _c0x, float _c0y, float _c1x, float _c1y) {
	tfx_attribute_node_t node;

	if (graph->nodes.size())
		node.frame = _frame;
	else
		node.frame = 0.f;

	node.value = _value;
	node.flags = flags;
	node.left.x = _c0x;
	node.left.y = _c0y;
	node.right.x = _c1x;
	node.right.y = _c1y;
	ClampNode(graph, &node);
	graph->nodes.push_back(node);
	SortGraph(graph);

	ReIndexGraph(graph);
	return &graph->nodes.back();
}

void AddGraphNode(tfx_graph_t *graph, tfx_attribute_node_t *node) {
	for (tfxBucketLoop(graph->nodes)) {
		if (graph->nodes[i].frame == node->frame && graph->nodes[i].value == node->value)
			return;
	}
	graph->nodes.push_back(*node);
	SortGraph(graph);
	ReIndexGraph(graph);
}

tfx_attribute_node_t *AddGraphCoordNode(tfx_graph_t *graph, float _frame, float _value) {
	tfx_attribute_node_t node;

	if (graph->nodes.size())
		node.frame = _frame;
	else
		node.frame = 0.f;

	node.value = _value;
	node.flags = 0;
	node.left.x = 0.f;
	node.left.y = 0.f;
	node.right.x = 0.f;
	node.right.y = 0.f;
	ClampNode(graph, &node);
	tfx_attribute_node_t &n = graph->nodes.push_back(node);
	if (SortGraph(graph)) {
		ReIndexGraph(graph);
		return graph->nodes.find(&n);
	}

	ReIndexGraph(graph);
	return &n;
}

tfx_attribute_node_t *InsertGraphCoordNode(tfx_graph_t *graph, float _frame, float _value) {
	tfx_attribute_node_t node;

	if (graph->nodes.size())
		node.frame = _frame;
	else
		node.frame = 0.f;

	node.value = _value;
	node.flags = 0;
	node.left.x = 0.f;
	node.left.y = 0.f;
	node.right.x = 0.f;
	node.right.y = 0.f;
	ClampNode(graph, &node);

	if (graph->nodes.size() > 1) {
		tfx_attribute_node_t *last_node = nullptr;
		for (auto *n = graph->nodes.begin() + 1; n != graph->nodes.end(); ++n) {
			if (node.frame < n->frame)
				last_node = n;
			else
				break;
		}

		if (last_node) {
			tfx_attribute_node_t *r_value = graph->nodes.insert(last_node, node);
			ReIndexGraph(graph);
			return r_value;
		}
	}

	tfx_attribute_node_t *r_value = &graph->nodes.push_back(node);
	ReIndexGraph(graph);
	return r_value;
}

tfx_attribute_node_t *InsertGraphNode(tfx_graph_t *graph, float _frame, float _value) {
	tfx_attribute_node_t node;

	if (graph->nodes.size())
		node.frame = _frame;
	else
		node.frame = 0.f;

	node.value = _value;
	node.flags = 0;
	node.left.x = 0.f;
	node.left.y = 0.f;
	node.right.x = 0.f;
	node.right.y = 0.f;
	ClampNode(graph, &node);

	if (graph->nodes.size() > 1) {
		tfx_attribute_node_t *last_node = nullptr;
		for (tfxBucketLoop(graph->nodes)) {
			if (node.frame < graph->nodes[i].frame) {
				last_node = &graph->nodes[i];
			}
			else {
				break;
			}
		}

		if (last_node) {
			tfx_attribute_node_t *r_value = graph->nodes.insert(last_node, node);
			ReIndexGraph(graph);
			return r_value;
		}
	}

	tfx_attribute_node_t *r_value = &graph->nodes.push_back(node);
	ReIndexGraph(graph);
	return r_value;
}

void SetGraphNode(tfx_graph_t *graph, tfxU32 i, float _frame, float _value, tfxAttributeNodeFlags flags, float _c0x, float _c0y, float _c1x, float _c1y) {
	if (!graph->nodes.empty() && i < graph->nodes.size()) {
		graph->nodes[i].frame = _frame;
		graph->nodes[i].value = _value;
		graph->nodes[i].flags = flags;
		graph->nodes[i].left.x = _c0x;
		graph->nodes[i].left.y = _c0y;
		graph->nodes[i].right.x = _c1x;
		graph->nodes[i].right.y = _c1y;
		if (SortGraph(graph))
			ReIndexGraph(graph);
	}
}

tfx_attribute_node_t *GraphNodeByIndex(tfx_graph_t *graph, tfxU32 index) {
	TFX_ASSERT(graph->nodes.current_size > index);    //Index is out of bounds
	return &graph->nodes[index];
}

float GraphValueByIndex(tfx_graph_t *graph, tfxU32 index) {
	TFX_ASSERT(graph->nodes.current_size > index);    //Index is out of bounds
	return graph->nodes[index].value;
}

float GraphFrameByIndex(tfx_graph_t *graph, tfxU32 index) {
	TFX_ASSERT(graph->nodes.current_size > index);    //Index is out of bounds
	return graph->nodes[index].frame;
}

float GetGraphValue(tfx_graph_t *graph, float age) {
	float lastv = 0;
	float lastf = 0;
	float p = 0;
	tfx_attribute_node_t *lastec = nullptr;
	for (tfxBucketLoop(graph->nodes)) {
		if (age < graph->nodes[i].frame) {
			p = (age - lastf) / (graph->nodes[i].frame - lastf);
			float bezier_value = GetBezierValue(lastec, &graph->nodes[i], p, graph->min.y, graph->max.y);
			if (bezier_value) {
				return bezier_value;
			}
			else {
				return lastv - p * (lastv - graph->nodes[i].value);
			}
		}
		lastv = graph->nodes[i].value;
		lastf = graph->nodes[i].frame;
		lastec = &graph->nodes[i];
	}
	return lastv;

}

tfx_attribute_node_t *GetGraphNextNode(tfx_graph_t *graph, tfx_attribute_node_t *node) {
	if (node->index < graph->nodes.size() - 1) {
		return &graph->nodes[node->index + 1];
	}

	return nullptr;
}

tfx_attribute_node_t *GetGraphPrevNode(tfx_graph_t *graph, tfx_attribute_node_t *node) {
	if (node->index > 0) {
		return &graph->nodes[node->index - 1];
	}

	return nullptr;
}

tfx_attribute_node_t *GetGraphLastNode(tfx_graph_t *graph) {
	return &graph->nodes.back();
}

float GetGraphRandomValue(tfx_graph_t *graph, float age, tfx_random_t *random) {
	float lastv = 0;
	float lastf = 0;
	float p = 0;
	tfx_attribute_node_t *lastec = nullptr;
	for (tfxBucketLoop(graph->nodes)) {
		if (age < graph->nodes[i].frame) {
			p = (age - lastf) / (graph->nodes[i].frame - lastf);
			float bezier_value = GetBezierValue(lastec, &graph->nodes[i], p, graph->min.y, graph->max.y);
			if (bezier_value) {
				return RandomRange(random, bezier_value);
			}
			else {
				return RandomRange(random, lastv - p * (lastv - graph->nodes[i].value));
			}
		}
		lastv = graph->nodes[i].value;
		lastf = graph->nodes[i].frame - 1;
		lastec = &graph->nodes[i];
	}
	return RandomRange(random, lastv);

}

float GetGraphValue(tfx_graph_t *graph, float age, float life) {
	float lastv = 0;
	float lastf = 0;
	float p = 0;
	tfx_attribute_node_t *lastec = nullptr;
	for (tfxBucketLoop(graph->nodes)) {
		float frame = graph->nodes[i].frame * life;
		if (age < frame) {
			p = (age - lastf) / (frame - lastf);
			float bezier_value = GetBezierValue(lastec, &graph->nodes[i], p, graph->min.y, graph->max.y);
			if (bezier_value) {
				return bezier_value;
			}
			else {
				return lastv - p * (lastv - graph->nodes[i].value);
			}
		}
		lastv = graph->nodes[i].value;
		lastf = frame - 1;
		lastec = &graph->nodes[i];
	}
	return lastv;
}

float GetGraphFirstValue(tfx_graph_t *graph) {
	if (graph->nodes.size())
		return graph->nodes.front().value;
	return 0.f;
}

float *LinkGraphFirstValue(tfx_graph_t *graph) {
	if (graph->nodes.size())
		return &graph->nodes.front().value;
	return nullptr;
}

float *LinkGraphLastValue(tfx_graph_t *graph) {
	if (graph->nodes.size())
		return &graph->nodes.back().value;
	return nullptr;
}

float GetGraphLastValue(tfx_graph_t *graph) {
	if (graph->nodes.size())
		return graph->nodes.back().value;

	return 0.f;
}

float GetGraphMaxValue(tfx_graph_t *graph) {
	if (graph->nodes.size()) {
		float value = tfxMIN_FLOAT;
		for (tfxBucketLoop(graph->nodes)) {
			if (value < graph->nodes[i].value)
				value = graph->nodes[i].value;
		}
		return value;
	}
	return 0.f;
}

float GetGraphMinValue(tfx_graph_t *graph) {
	if (graph->nodes.size()) {
		float value = tfxMAX_FLOAT;
		for (tfxBucketLoop(graph->nodes)) {
			if (value > graph->nodes[i].value)
				value = graph->nodes[i].value;
		}
		return value;
	}
	return 0.f;
}

float GetGraphLastFrame(tfx_graph_t *graph, float update_frequency) {
	float frame_length = 1000.f / update_frequency;
	if (graph->nodes.size()) {
		return graph->nodes.size() > 1 && graph->nodes.back().frame == 0 ? frame_length : graph->nodes.back().frame;
	}

	return 0.f;
}

tfx_attribute_node_t *FindGraphNode(tfx_graph_t *graph, tfx_attribute_node_t *n) {
	return graph->nodes.find(n);
}

void ValidateGraphCurves(tfx_graph_t *graph) {
	tfxU32 index = 0;
	tfxU32 last_index = graph->nodes.size() - 1;
	for (tfxBucketLoop(graph->nodes)) {
		if (graph->nodes[i].flags & tfxAttributeNodeFlags_is_curve) {
			if (index < last_index) {
				if (graph->nodes[index + 1].frame < graph->nodes[i].right.x)
					graph->nodes[i].right.x = graph->nodes[index + 1].frame;
			}
			if (index > 0) {
				if (graph->nodes[index - 1].frame > graph->nodes[i].left.x)
					graph->nodes[i].left.x = graph->nodes[index - 1].frame;
			}
			if (graph->nodes[i].left.x > graph->nodes[i].frame)
				graph->nodes[i].left.x = graph->nodes[i].frame;
			if (graph->nodes[i].right.x < graph->nodes[i].frame)
				graph->nodes[i].right.x = graph->nodes[i].frame;
		}
		index++;
	}
}

void DeleteGraphNode(tfx_graph_t *graph, tfx_attribute_node_t *n) {
	graph->nodes.erase(n);
	ReIndexGraph(graph);
}

void DeleteGraphNodeAtFrame(tfx_graph_t *graph, float frame) {
	for (tfxBucketLoop(graph->nodes)) {
		if (graph->nodes[i].frame == frame) {
			graph->nodes.erase(&graph->nodes[i]);
			ReIndexGraph(graph);
			return;
		}
	}
}

void ResetGraphNodes(tfx_graph_t *graph, float v, tfx_graph_preset preset, bool add_node) {
	graph->nodes.clear();
	graph->nodes.trim_buckets();
	if (add_node && preset == tfxWeightOvertimePreset) {
		AddGraphNode(graph, 0.f, 0.f, 0);
		tfx_attribute_node_t *node = AddGraphNode(graph, 1.f, 1.f, tfxAttributeNodeFlags_is_curve, 0.f, 1.f, 1.f, 1.f);
		SetNodeCurveInitialised(node);
	}
	else if (add_node) {
		if (preset == tfxWeightOvertimePreset) {
			AddGraphNode(graph, 0.f, 0.f, 0);
			tfx_attribute_node_t *node = AddGraphNode(graph, 1.f, 1.f, tfxAttributeNodeFlags_is_curve, 0.f, 1.f, 1.f, 1.f);
			SetNodeCurveInitialised(node);
		}
		else {
			AddGraphNode(graph, 0.f, v);
		}
	}
}

void ResetGraph(tfx_graph_t *graph, float v, tfx_graph_preset preset, bool add_node, float max_frames) {
	graph->nodes.clear();
	graph->nodes.trim_buckets();
	if (add_node && preset == tfxWeightOvertimePreset) {
		AddGraphNode(graph, 0.f, 0.f, 0);
		tfx_attribute_node_t *node = AddGraphNode(graph, 1.f, 1.f, tfxAttributeNodeFlags_is_curve, 0.f, 1.f, 1.f, 1.f);
		SetNodeCurveInitialised(node);
	}
	else if (add_node) {
		if (preset == tfxWeightOvertimePreset) {
			AddGraphNode(graph, 0.f, 0.f, 0);
			tfx_attribute_node_t *node = AddGraphNode(graph, 1.f, 1.f, tfxAttributeNodeFlags_is_curve, 0.f, 1.f, 1.f, 1.f);
			SetNodeCurveInitialised(node);
		}
		else {
			AddGraphNode(graph, 0.f, v);
		}
	}
	if (!max_frames) {
		max_frames = tfxMAX_FRAME;
	}
	switch (preset) {
	case tfx_graph_preset::tfxGlobalPercentPreset:
		//We have a epsilon to prevent divide by 0 here
		graph->min = { 0.f, 0.0001f }; graph->max = { max_frames, 20.f };
		break;
	case tfx_graph_preset::tfxGlobalOpacityPreset:
		graph->min = { 0.f, 0.f }; graph->max = { max_frames, 1.f };
		break;
	case tfx_graph_preset::tfxGlobalPercentPresetSigned:
		graph->min = { 0.f, -20.f }; graph->max = { max_frames, 20.f };
		break;
	case tfx_graph_preset::tfxAnglePreset:
		graph->min = { 0.f, -1080.f }; graph->max = { max_frames, 1080.f };
		break;
	case tfx_graph_preset::tfxArcPreset:
		graph->min = { 0.f, 0.f }; graph->max = { max_frames, 360.f };
		break;
	case tfx_graph_preset::tfxEmissionRangePreset:
		graph->min = { 0.f, 0.f }; graph->max = { max_frames, 360.f };
		break;
	case tfx_graph_preset::tfxDimensionsPreset:
		graph->min = { 0.f, 0.f }; graph->max = { max_frames, 4000.f };
		break;
	case tfx_graph_preset::tfxTranslationPreset:
		graph->min = { 0.f, -4000.f }; graph->max = { max_frames, 4000.f };
		break;
	case tfx_graph_preset::tfxLifePreset:
		//We have a epsilon to prevent divide by 0 here. The divide by zero occurrs in control functions (ControlParticleImageFrame3d etc.) when the current % life of the particle is calculated
		graph->min = { 0.f, 0.0001f }; graph->max = { max_frames, 100000.f };
		break;
	case tfx_graph_preset::tfxAmountPreset:
		graph->min = { 0.f, 0.f }; graph->max = { max_frames, 5000.f };
		break;
	case tfx_graph_preset::tfxVelocityPreset:
		graph->min = { 0.f, 0.f }; graph->max = { max_frames, 10000.f };
		break;
	case tfx_graph_preset::tfxWeightPreset:
		graph->min = { 0.f, -10000.f }; graph->max = { max_frames, 10000.f };
		break;
	case tfx_graph_preset::tfxWeightVariationPreset:
		graph->min = { 0.f, 0.f }; graph->max = { max_frames, 20000.f };
		break;
	case tfx_graph_preset::tfxNoiseOffsetVariationPreset:
		graph->min = { 0.f, 0.f }; graph->max = { max_frames, 1000.f };
		break;
	case tfx_graph_preset::tfxNoiseResolutionPreset:
		graph->min = { 0.f, 0.f }; graph->max = { max_frames, 10000.f };
		break;
	case tfx_graph_preset::tfxSpinPreset:
		graph->min = { 0.f, -2000.f }; graph->max = { max_frames, 2000.f };
		break;
	case tfx_graph_preset::tfxSpinVariationPreset:
		graph->min = { 0.f, 0.f }; graph->max = { max_frames, 2000.f };
		break;
	case tfx_graph_preset::tfxDirectionVariationPreset:
		graph->min = { 0.f, 0.f }; graph->max = { max_frames, 22.5f };
		break;
	case tfx_graph_preset::tfxWeightOvertimePreset:
		graph->min = { 0.f, -20.f }; graph->max = { 1.f, 20.f };
		break;
	case tfx_graph_preset::tfxDirectionOvertimePreset:
		graph->min = { 0.f, 0.f }; graph->max = { 1.f, 4320.f };
		break;
	case tfx_graph_preset::tfxSpinOvertimePreset:
		graph->min = { 0.f, 0.f }; graph->max = { 1.f, 20.f };
		break;
	case tfx_graph_preset::tfxVelocityOvertimePreset:
		graph->min = { 0.f, -20.f }; graph->max = { 1.f, 20.f };
		break;
	case tfx_graph_preset::tfxPercentOvertime:
		graph->min = { 0.f, 0.f }; graph->max = { 1.f, 20.f };
		break;
	case tfx_graph_preset::tfxFrameratePreset:
		graph->min = { 0.f, 0.f }; graph->max = { 1.f, 200.f };
		break;
	case tfx_graph_preset::tfxVelocityTurbulancePreset:
		graph->min = { 0.f, 0.f }; graph->max = { 1.f, 2000.f };
		break;
	case tfx_graph_preset::tfxOpacityOvertimePreset:
		graph->min = { 0.f, 0.f }; graph->max = { 1.f, 1.f };
		break;
	case tfx_graph_preset::tfxColorPreset:
		graph->min = { 0.f, 0.f }; graph->max = { 1.f, 255.f };
		break;
	case tfx_graph_preset::tfxIntensityOvertimePreset:
		graph->min = { 0.f, 0.f }; graph->max = { 1.f, 10.f };
		break;
	case tfx_graph_preset::tfxPathDirectionOvertimePreset:
		graph->min = { 0.f, -4320.f }; graph->max = { 1.f, 4320.f };
		break;
	case tfx_graph_preset::tfxPathTranslationOvertimePreset:
		graph->min = { 0.f, -1000.f }; graph->max = { 1.f, 1000.f };
		break;
	}

	graph->graph_preset = preset;
}

tfx_vec4_t GetMinMaxGraphValues(tfx_graph_preset preset) {
	tfx_vec4_t mm;
	switch (preset) {
	case tfx_graph_preset::tfxGlobalPercentPreset:
		mm = { 0.f, 0.f, tfxMAX_FRAME, 20.f };
		break;
	case tfx_graph_preset::tfxGlobalOpacityPreset:
		mm = { 0.f, 0.f , tfxMAX_FRAME, 1.f };
		break;
	case tfx_graph_preset::tfxGlobalPercentPresetSigned:
		mm = { 0.f, -20.f, tfxMAX_FRAME, 20.f };
		break;
	case tfx_graph_preset::tfxAnglePreset:
		mm = { 0.f, -1080.f, tfxMAX_FRAME, 1080.f };
		break;
	case tfx_graph_preset::tfxArcPreset:
		mm = { 0.f, 0.f , tfxMAX_FRAME, 360.f };
		break;
	case tfx_graph_preset::tfxEmissionRangePreset:
		mm = { 0.f, 0.f, tfxMAX_FRAME, 360.f };
		break;
	case tfx_graph_preset::tfxDimensionsPreset:
		mm = { 0.f, 0.f, tfxMAX_FRAME, 4000.f };
		break;
	case tfx_graph_preset::tfxTranslationPreset:
		mm = { 0.f, -4000.f, tfxMAX_FRAME, 4000.f };
		break;
	case tfx_graph_preset::tfxLifePreset:
		mm = { 0.f, 0.f, tfxMAX_FRAME, 100000.f };
		break;
	case tfx_graph_preset::tfxAmountPreset:
		mm = { 0.f, 0.f, tfxMAX_FRAME, 5000.f };
		break;
	case tfx_graph_preset::tfxVelocityPreset:
		mm = { 0.f, 0.f, tfxMAX_FRAME, 10000.f };
		break;
	case tfx_graph_preset::tfxWeightPreset:
		mm = { 0.f, -10000.f, tfxMAX_FRAME, 10000.f };
		break;
	case tfx_graph_preset::tfxWeightVariationPreset:
		mm = { 0.f, 0.f, tfxMAX_FRAME, 20000.f };
		break;
	case tfx_graph_preset::tfxSpinPreset:
		mm = { 0.f, -2000.f, tfxMAX_FRAME, 2000.f };
		break;
	case tfx_graph_preset::tfxSpinVariationPreset:
		mm = { 0.f, 0.f, tfxMAX_FRAME, 2000.f };
		break;
	case tfx_graph_preset::tfxDirectionVariationPreset:
		mm = { 0.f, 0.f, tfxMAX_FRAME, 22.5f };
		break;
	case tfx_graph_preset::tfxWeightOvertimePreset:
		mm = { 0.f, -20.f, 1.f, 20.f };
		break;
	case tfx_graph_preset::tfxDirectionOvertimePreset:
		mm = { 0.f, 0.f, 1.f, 4320.f };
		break;
	case tfx_graph_preset::tfxSpinOvertimePreset:
		mm = { 0.f, 0.f, 1.f, 20.f };
		break;
	case tfx_graph_preset::tfxVelocityOvertimePreset:
		mm = { 0.f, -20.f, 1.f, 20.f };
		break;
	case tfx_graph_preset::tfxPercentOvertime:
		mm = { 0.f, 0.f, 1.f, 20.f };
		break;
	case tfx_graph_preset::tfxFrameratePreset:
		mm = { 0.f, 0.f, 1.f, 200.f };
		break;
	case tfx_graph_preset::tfxVelocityTurbulancePreset:
		mm = { 0.f, 0.f, 1.f, 2000.f };
		break;
	case tfx_graph_preset::tfxOpacityOvertimePreset:
		mm = { 0.f, 0.f, 1.f, 1.f };
		break;
	case tfx_graph_preset::tfxColorPreset:
		mm = { 0.f, 0.f, 1.f, 255.f };
		break;
	case tfx_graph_preset::tfxIntensityOvertimePreset:
		mm = { 0.f, 0.f, 1.f, 10.f };
		break;
	default:
		mm = { 0.f, 0.f, tfxMAX_FRAME, 20.f };
		break;
	}

	return mm;
}

void DragGraphValues(tfx_graph_preset preset, float *frame, float *value) {
	switch (preset) {
	case tfx_graph_preset::tfxOpacityOvertimePreset:
	case tfx_graph_preset::tfxGlobalPercentPreset:
	case tfx_graph_preset::tfxIntensityOvertimePreset:
		*frame = 0.001f;
		*value = 0.001f;
		break;
	case tfx_graph_preset::tfxDirectionOvertimePreset:
		*frame = 0.001f;
		*value = 0.1f;
		break;
	case tfx_graph_preset::tfxLifePreset:
		*frame = 5;
		*value = 5;
		break;
	case tfx_graph_preset::tfxAnglePreset:
	case tfx_graph_preset::tfxArcPreset:
	case tfx_graph_preset::tfxEmissionRangePreset:
		*frame = 5;
		*value = 0.1f;
		break;
	case tfx_graph_preset::tfxDimensionsPreset:
	case tfx_graph_preset::tfxAmountPreset:
	case tfx_graph_preset::tfxVelocityPreset:
	case tfx_graph_preset::tfxWeightPreset:
	case tfx_graph_preset::tfxWeightVariationPreset:
	case tfx_graph_preset::tfxSpinPreset:
	case tfx_graph_preset::tfxSpinVariationPreset:
	case tfx_graph_preset::tfxFrameratePreset:
		*frame = 5.f;
		*value = 1.f;
		break;
	case tfx_graph_preset::tfxVelocityTurbulancePreset:
		*frame = .001f;
		*value = .01f;
		break;
	case tfx_graph_preset::tfxWeightOvertimePreset:
	case tfx_graph_preset::tfxVelocityOvertimePreset:
	case tfx_graph_preset::tfxSpinOvertimePreset:
	case tfx_graph_preset::tfxDirectionVariationPreset:
		*frame = 0.001f;
		*value = 0.01f;
		break;
	case tfx_graph_preset::tfxColorPreset:
		*frame = 0.001f;
		*value = 1.f;
		break;
	case tfx_graph_preset::tfxPercentOvertime:
		*frame = 0.05f;
		*value = 0.05f;
		break;
	default:
		*frame = 1;
		*value = 0.1f;
		break;
	}
}

void ClearGraphToOne(tfx_graph_t *graph, float value) {
	graph->nodes.clear();
	AddGraphNode(graph, 0.f, value);
}

void ClearGraph(tfx_graph_t *graph) {
	graph->nodes.clear();
}

void FreeGraph(tfx_graph_t *graph) {
	//Explicitly free the nodes
	graph->nodes.free_all();
	graph->lookup.values.free();
}

void CopyGraph(tfx_graph_t *from, tfx_graph_t *to, bool compile) {
	ClearGraph(to);
	for (tfxBucketLoop(from->nodes)) {
		to->nodes.push_back(from->nodes[i]);
	}
	if (compile) {
		if (IsColorGraph(from)) {
			CompileColorOvertime(to);
		}
		else if (IsOvertimeGraph(from)) {
			CompileGraphOvertime(to);
		}
		else {
			CompileGraph(to);
		}
	}
}

void CopyGraphColor(tfx_overtime_attributes_t *from, tfx_overtime_attributes_t *to) {
	ClearGraph(&to->red);
	ClearGraph(&to->green);
	ClearGraph(&to->blue);
	for (tfxBucketLoop(from->red.nodes)) {
		to->red.nodes.push_back(from->red.nodes[i]);
	}
	for (tfxBucketLoop(from->green.nodes)) {
		to->green.nodes.push_back(from->green.nodes[i]);
	}
	for (tfxBucketLoop(from->blue.nodes)) {
		to->blue.nodes.push_back(from->blue.nodes[i]);
	}
}

void CopyGraphColorHint(tfx_overtime_attributes_t *from, tfx_overtime_attributes_t *to) {
	ClearGraph(&to->red_hint);
	ClearGraph(&to->green_hint);
	ClearGraph(&to->blue_hint);
	for (tfxBucketLoop(from->red_hint.nodes)) {
		to->red.nodes.push_back(from->red_hint.nodes[i]);
	}
	for (tfxBucketLoop(from->green_hint.nodes)) {
		to->green_hint.nodes.push_back(from->green_hint.nodes[i]);
	}
	for (tfxBucketLoop(from->blue.nodes)) {
		to->blue_hint.nodes.push_back(from->blue_hint.nodes[i]);
	}
}
void CopyGraphColors(tfx_graph_t *from_red, tfx_graph_t *from_green, tfx_graph_t *from_blue, tfx_graph_t *to_red, tfx_graph_t *to_green, tfx_graph_t *to_blue) {
	ClearGraph(to_red);
	ClearGraph(to_green);
	ClearGraph(to_blue);
	for (tfxBucketLoop(from_red->nodes)) {
		to_red->nodes.push_back(from_red->nodes[i]);
	}
	for (tfxBucketLoop(from_green->nodes)) {
		to_green->nodes.push_back(from_green->nodes[i]);
	}
	for (tfxBucketLoop(from_blue->nodes)) {
		to_blue->nodes.push_back(from_blue->nodes[i]);
	}
}

bool SortGraph(tfx_graph_t *graph) {
	bool needed_sorting = false;
	for (tfxU32 i = 1; i < graph->nodes.current_size; ++i) {
		tfx_attribute_node_t key = graph->nodes[i];
		int j = i - 1;
		while (j >= 0 && key.frame < graph->nodes[j].frame) {
			graph->nodes[j + 1] = graph->nodes[j];
			--j;
			needed_sorting = true;
		}
		graph->nodes[j + 1] = key;
	}
	return needed_sorting;
}

void FlipGraph(tfx_graph_t *graph) {
	if (graph->nodes.size() == 1) {
		return;
	}
	int left = 0;
	int right = graph->nodes.size() - 1;
	while (left < right) {
		float left_value = graph->nodes[left].value;
		graph->nodes[left++].value = graph->nodes[right].value;
		graph->nodes[right--].value = left_value;
	}
}

void ReIndexGraph(tfx_graph_t *graph) {
	tfxU32 index = 0;
	for (tfxBucketLoop(graph->nodes)) {
		graph->nodes[i].index = index++;
	}
}

tfx_vec2_t GetGraphInitialZoom(tfx_graph_t *graph) {
	switch (graph->graph_preset) {
	case tfx_graph_preset::tfxOpacityOvertimePreset:
		return tfx_vec2_t(.0017f, 0.00275f);
	case tfx_graph_preset::tfxGlobalPercentPreset:
		return tfx_vec2_t(10.f, 0.005f);
		break;
	case tfx_graph_preset::tfxGlobalPercentPresetSigned:
		return tfx_vec2_t(10.f, 0.006f);
		break;
	case tfx_graph_preset::tfxGlobalOpacityPreset:
		return tfx_vec2_t(10.f, 0.003f);
		break;
	case tfx_graph_preset::tfxLifePreset:
		return tfx_vec2_t(10.f, 3.5f);
		break;
	case tfx_graph_preset::tfxAnglePreset:
		return tfx_vec2_t(10.f, 1.f);
		break;
	case tfx_graph_preset::tfxArcPreset:
		return tfx_vec2_t(10.f, 1.f);
		break;
	case tfx_graph_preset::tfxEmissionRangePreset:
		return tfx_vec2_t(10.f, .5f);
		break;
	case tfx_graph_preset::tfxAmountPreset:
		return tfx_vec2_t(10.f, 1.25f);
		break;
	case tfx_graph_preset::tfxFrameratePreset:
		return tfx_vec2_t(0.0017f, .5f);
		break;
	case tfx_graph_preset::tfxVelocityTurbulancePreset:
		return tfx_vec2_t(0.0017f, .1f);
		break;
	case tfx_graph_preset::tfxDimensionsPreset:
	case tfx_graph_preset::tfxTranslationPreset:
	case tfx_graph_preset::tfxVelocityPreset:
	case tfx_graph_preset::tfxWeightPreset:
	case tfx_graph_preset::tfxWeightVariationPreset:
	case tfx_graph_preset::tfxSpinPreset:
	case tfx_graph_preset::tfxSpinVariationPreset:
		return tfx_vec2_t(10.f, 2.5f);
		break;
	case tfx_graph_preset::tfxNoiseResolutionPreset:
		return tfx_vec2_t(10.f, 1.f);
		break;
	case tfx_graph_preset::tfxNoiseOffsetVariationPreset:
		return tfx_vec2_t(10.f, .01f);
		break;
	case tfx_graph_preset::tfxDirectionOvertimePreset:
		return tfx_vec2_t(0.0017f, 1.f);
		break;
	case tfx_graph_preset::tfxWeightOvertimePreset:
	case tfx_graph_preset::tfxVelocityOvertimePreset:
	case tfx_graph_preset::tfxSpinOvertimePreset:
	case tfx_graph_preset::tfxDirectionVariationPreset:
	case tfx_graph_preset::tfxPercentOvertime:
		return tfx_vec2_t(0.0017f, 0.0035f);
		break;
	case tfx_graph_preset::tfxIntensityOvertimePreset:
		return tfx_vec2_t(0.0017f, 0.01115f);
		break;
	case tfx_graph_preset::tfxColorPreset:
		break;
	case tfx_graph_preset::tfxPathDirectionOvertimePreset:
		return tfx_vec2_t(0.0017f, 5.f);
		break;
	case tfx_graph_preset::tfxPathTranslationOvertimePreset:
		return tfx_vec2_t(0.0017f, 2.5f);
		break;
	default:
		return tfx_vec2_t(0.1f, 0.1f);
		break;
	}

	return tfx_vec2_t(0.1f, 0.1f);
}

tfx_vec2_t GetGraphInitialZoom3d(tfx_graph_t *graph) {
	switch (graph->graph_preset) {
	case tfx_graph_preset::tfxOpacityOvertimePreset:
		return tfx_vec2_t(.0017f, 0.00275f);
	case tfx_graph_preset::tfxGlobalPercentPreset:
		return tfx_vec2_t(10.f, 0.005f);
		break;
	case tfx_graph_preset::tfxGlobalPercentPresetSigned:
		return tfx_vec2_t(10.f, 0.006f);
		break;
	case tfx_graph_preset::tfxGlobalOpacityPreset:
		return tfx_vec2_t(10.f, 0.003f);
		break;
	case tfx_graph_preset::tfxLifePreset:
		return tfx_vec2_t(10.f, 3.5f);
		break;
	case tfx_graph_preset::tfxAnglePreset:
		return tfx_vec2_t(10.f, 1.f);
		break;
	case tfx_graph_preset::tfxArcPreset:
		return tfx_vec2_t(10.f, 1.f);
		break;
	case tfx_graph_preset::tfxEmissionRangePreset:
		return tfx_vec2_t(10.f, .5f);
		break;
	case tfx_graph_preset::tfxAmountPreset:
		return tfx_vec2_t(10.f, 1.25f);
		break;
	case tfx_graph_preset::tfxFrameratePreset:
		return tfx_vec2_t(0.0017f, .01f);
		break;
	case tfx_graph_preset::tfxDimensionsPreset:
	case tfx_graph_preset::tfxTranslationPreset:
	case tfx_graph_preset::tfxVelocityPreset:
		return tfx_vec2_t(10.f, 0.01f);
		break;
	case tfx_graph_preset::tfxWeightPreset:
	case tfx_graph_preset::tfxWeightVariationPreset:
		return tfx_vec2_t(10.f, 1.f);
		break;
	case tfx_graph_preset::tfxVelocityTurbulancePreset:
		return tfx_vec2_t(0.0017f, .01f);
		break;
	case tfx_graph_preset::tfxSpinPreset:
	case tfx_graph_preset::tfxSpinVariationPreset:
		return tfx_vec2_t(10.f, 2.5f);
		break;
	case tfx_graph_preset::tfxNoiseResolutionPreset:
		return tfx_vec2_t(10.f, .01f);
		break;
	case tfx_graph_preset::tfxNoiseOffsetVariationPreset:
		return tfx_vec2_t(10.f, .01f);
		break;
	case tfx_graph_preset::tfxDirectionOvertimePreset:
		return tfx_vec2_t(0.0017f, 1.f);
		break;
	case tfx_graph_preset::tfxWeightOvertimePreset:
	case tfx_graph_preset::tfxVelocityOvertimePreset:
	case tfx_graph_preset::tfxSpinOvertimePreset:
	case tfx_graph_preset::tfxDirectionVariationPreset:
	case tfx_graph_preset::tfxPercentOvertime:
		return tfx_vec2_t(0.0017f, 0.0035f);
		break;
	case tfx_graph_preset::tfxIntensityOvertimePreset:
		return tfx_vec2_t(0.0017f, 0.01115f);
		break;
	case tfx_graph_preset::tfxColorPreset:
		break;
	case tfx_graph_preset::tfxPathDirectionOvertimePreset:
		return tfx_vec2_t(0.0017f, 5.f);
		break;
	case tfx_graph_preset::tfxPathTranslationOvertimePreset:
		return tfx_vec2_t(0.0017f, .05f);
		break;
	default:
		return tfx_vec2_t(0.1f, 0.1f);
		break;
	}

	return tfx_vec2_t(0.1f, 0.1f);
}

void CompileGraph(tfx_graph_t *graph) {
	float last_frame = GetGraphLastFrame(graph, 60.f);
	graph->lookup.last_frame = tfxU32(last_frame / tfxLOOKUP_FREQUENCY);
	if (graph->lookup.last_frame) {
		graph->lookup.values.resize(graph->lookup.last_frame + 1);
		for (tfxU32 f = 0; f != graph->lookup.last_frame + 1; ++f) {
			graph->lookup.values[f] = GetGraphValue(graph, (float)f * tfxLOOKUP_FREQUENCY);
		}
		graph->lookup.values[graph->lookup.last_frame] = GetGraphLastValue(graph);
	}
	else {
		graph->lookup.values.resize(1);
		graph->lookup.values[0] = GetGraphFirstValue(graph);
	}
}

void CompileGraphOvertime(tfx_graph_t *graph) {
	if (graph->type == tfxOvertime_intensity || graph->type == tfxOvertime_curved_alpha || graph->type == tfxOvertime_alpha_sharpness) {
		CompileGraphRampOvertime(graph);
		return;
	}
	if (graph->nodes.size() > 1) {
		graph->lookup.last_frame = tfxU32(graph->lookup.life / tfxLOOKUP_FREQUENCY_OVERTIME);
		graph->lookup.values.resize(graph->lookup.last_frame + 1);
		for (tfxU32 f = 0; f != graph->lookup.last_frame + 1; ++f) {
			graph->lookup.values[f] = GetGraphValue(graph, (float)f * tfxLOOKUP_FREQUENCY_OVERTIME, graph->lookup.life);
		}
		graph->lookup.values[graph->lookup.last_frame] = GetGraphLastValue(graph);
	}
	else {
		graph->lookup.last_frame = 0;
		graph->lookup.values.resize(1);
		graph->lookup.values[0] = GetGraphFirstValue(graph);
	}
}

void CompileGraphRampOvertime(tfx_graph_t *graph) {
	//if (graph->nodes.size() > 1) {
	graph->lookup.last_frame = tfxCOLOR_RAMP_WIDTH - 1;
	graph->lookup.values.resize(tfxCOLOR_RAMP_WIDTH);
	for (tfxU32 f = 0; f != tfxCOLOR_RAMP_WIDTH; ++f) {
		float age = ((float)f / tfxCOLOR_RAMP_WIDTH) * graph->lookup.life;
		graph->lookup.values[f] = GetGraphValue(graph, age, graph->lookup.life);
	}
	graph->lookup.values[graph->lookup.last_frame] = GetGraphLastValue(graph);
	//}
	//else {
		//graph->lookup.last_frame = 0;
		//graph->lookup.values.resize(1);
		//graph->lookup.values[0] = GetGraphFirstValue(graph);
	//}
}

void CompileColorOvertime(tfx_graph_t *graph, float gamma) {
	if (graph->nodes.size() > 1) {
		graph->lookup.last_frame = tfxU32(graph->lookup.life / tfxLOOKUP_FREQUENCY_OVERTIME);
		graph->lookup.values.resize(graph->lookup.last_frame + 1);
		for (tfxU32 f = 0; f != graph->lookup.last_frame + 1; ++f) {
			graph->lookup.values[f] = tfx__gamma_correct(GetGraphValue(graph, (float)f * tfxLOOKUP_FREQUENCY_OVERTIME, graph->lookup.life), gamma);
		}
		graph->lookup.values[graph->lookup.last_frame] = tfx__gamma_correct(GetGraphLastValue(graph), gamma);
	}
	else {
		graph->lookup.last_frame = 0;
		graph->lookup.values.resize(1);
		graph->lookup.values[0] = tfx__gamma_correct(GetGraphFirstValue(graph), gamma);
	}
}

tfxKey HashColorRamp(tfx_color_ramp_t *ramp) {
	tfxKey hash = tfxXXHash64::hash(ramp->colors, tfxCOLOR_RAMP_WIDTH * sizeof(tfx_rgba8_t), 0);
	return hash;
}

tfx_bitmap_t tfxCreateBitmap(int width, int height, int channels) {
	TFX_ASSERT(width > 0 && height > 0 && channels > 0);
	tfx_bitmap_t bitmap = {};
    bitmap.size = width * height * channels;
    if (bitmap.size > 0) {
        bitmap.data = (tfx_byte*)tfxALLOCATE_ALIGNED(bitmap.size, 16);
        bitmap.width = width;
        bitmap.height = height;
        bitmap.channels = channels;
        bitmap.stride = width * channels;
    }
	return bitmap;
}
void tfxPlotBitmap(tfx_bitmap_t *image, int x, int y, tfx_rgba8_t color) {

    tfx_size pos = y * image->stride + (x * image->channels);

    if (pos >= image->size) {
        return;
    }

    if (image->channels == 4) {
		*(image->data + pos) = color.r;
		*(image->data + pos + 1) = color.g;
		*(image->data + pos + 2) = color.b;
		*(image->data + pos + 3) = color.a;
    }
    else if (image->channels == 3) {
		*(image->data + pos) = color.r;
		*(image->data + pos + 1) = color.g;
		*(image->data + pos + 2) = color.b;
    }
    else if (image->channels == 2) {
		*(image->data + pos) = color.r;
		*(image->data + pos + 3) = color.a;
    }
    else if (image->channels == 1) {
		*(image->data + pos) = color.r;
    }

}

void tfxFreeBitmap(tfx_bitmap_t *bitmap) {
	if (bitmap->data) {
		tfxFREE(bitmap->data);
	}
}

tfx_color_ramp_t CompileColorRamp(tfx_overtime_attributes_t *attributes, float gamma) {
	float r, g, b, a;
	tfx_color_ramp_t color_ramp = attributes->color_ramps[0];
	if (color_ramp.flags & tfxColorRampFlags_use_sinusoidal_ramp_generation) {
		for (tfxU32 f = 0; f != tfxCOLOR_RAMP_WIDTH; ++f) {
			float x = (float)f / tfxCOLOR_RAMP_WIDTH;
			float age = x * attributes->red.lookup.life;

			r = color_ramp.brightness.x + color_ramp.contrast.x * cosf(tfxPI2 * (color_ramp.frequency.x * x + color_ramp.offsets.x));
			g = color_ramp.brightness.y + color_ramp.contrast.y * cosf(tfxPI2 * (color_ramp.frequency.y * x + color_ramp.offsets.y));
			b = color_ramp.brightness.z + color_ramp.contrast.z * cosf(tfxPI2 * (color_ramp.frequency.z * x + color_ramp.offsets.z));
			a = tfx__gamma_correct(GetGraphValue(&attributes->blendfactor, age, attributes->blendfactor.lookup.life), gamma);
			color_ramp.colors[f].r = tfxU32(r * 255.f);
			color_ramp.colors[f].g = tfxU32(g * 255.f);
			color_ramp.colors[f].b = tfxU32(b * 255.f);
			color_ramp.colors[f].a = tfxU32(a * 255.f);
		}
	} else {
		for (tfxU32 f = 0; f != tfxCOLOR_RAMP_WIDTH; ++f) {
			float age = ((float)f / tfxCOLOR_RAMP_WIDTH) * attributes->red.lookup.life;
			r = tfx__gamma_correct(GetGraphValue(&attributes->red, age, attributes->red.lookup.life), gamma);
			g = tfx__gamma_correct(GetGraphValue(&attributes->green, age, attributes->green.lookup.life), gamma);
			b = tfx__gamma_correct(GetGraphValue(&attributes->blue, age, attributes->blue.lookup.life), gamma);
			a = tfx__gamma_correct(GetGraphValue(&attributes->blendfactor, age, attributes->blendfactor.lookup.life), gamma);
			color_ramp.colors[f].r = tfxU32(r * 255.f);
			color_ramp.colors[f].g = tfxU32(g * 255.f);
			color_ramp.colors[f].b = tfxU32(b * 255.f);
			color_ramp.colors[f].a = tfxU32(a * 255.f);
		}
		color_ramp.colors[tfxCOLOR_RAMP_WIDTH - 1].r = tfxU32(GetGraphLastValue(&attributes->red) * 255.f);
		color_ramp.colors[tfxCOLOR_RAMP_WIDTH - 1].g = tfxU32(GetGraphLastValue(&attributes->green) * 255.f);
		color_ramp.colors[tfxCOLOR_RAMP_WIDTH - 1].b = tfxU32(GetGraphLastValue(&attributes->blue) * 255.f);
		color_ramp.colors[tfxCOLOR_RAMP_WIDTH - 1].a = tfxU32(GetGraphLastValue(&attributes->blendfactor) * 255.f);
	}
	return color_ramp;
}

tfx_color_ramp_t CompileColorRampHint(tfx_overtime_attributes_t *attributes, float gamma) {
	float r, g, b, a;
	tfx_color_ramp_t color_ramp_hint;
	for (tfxU32 f = 0; f != tfxCOLOR_RAMP_WIDTH; ++f) {
		float age = ((float)f / tfxCOLOR_RAMP_WIDTH) * attributes->red.lookup.life;
		r = tfx__gamma_correct(GetGraphValue(&attributes->red_hint, age, attributes->red.lookup.life), gamma);
		g = tfx__gamma_correct(GetGraphValue(&attributes->green_hint, age, attributes->green.lookup.life), gamma);
		b = tfx__gamma_correct(GetGraphValue(&attributes->blue_hint, age, attributes->blue.lookup.life), gamma);
		a = tfx__gamma_correct(GetGraphValue(&attributes->blendfactor_hint, age, attributes->blendfactor_hint.lookup.life), gamma);
		color_ramp_hint.colors[f].r = tfxU32(r * 255.f);
		color_ramp_hint.colors[f].g = tfxU32(g * 255.f);
		color_ramp_hint.colors[f].b = tfxU32(b * 255.f);
		color_ramp_hint.colors[f].a = tfxU32(a * 255.f);
	}
	color_ramp_hint.colors[tfxCOLOR_RAMP_WIDTH - 1].r = tfxU32(GetGraphLastValue(&attributes->red_hint) * 255.f);
	color_ramp_hint.colors[tfxCOLOR_RAMP_WIDTH - 1].g = tfxU32(GetGraphLastValue(&attributes->green_hint) * 255.f);
	color_ramp_hint.colors[tfxCOLOR_RAMP_WIDTH - 1].b = tfxU32(GetGraphLastValue(&attributes->blue_hint) * 255.f);
	color_ramp_hint.colors[tfxCOLOR_RAMP_WIDTH - 1].a = tfxU32(GetGraphLastValue(&attributes->blendfactor_hint) * 255.f);
	return color_ramp_hint;
}

void PlotColorRamp(tfx_bitmap_t *bitmap, tfx_color_ramp_t *ramp, tfxU32 y) {
	for (int x = 0; x != tfxCOLOR_RAMP_WIDTH; ++x) {
		tfx_rgba8_t color = ramp->colors[x];
		if (tfxStore->color_ramp_format == tfx_color_format_bgra) {
			color.color = tfx_SWIZZLE_RGBA_TO_BGRA(color.color);
		}
		tfxPlotBitmap(bitmap, x, y, color);
	}
}

tfxU32 AddColorRampToBitmaps(tfx_color_ramp_bitmap_data_t *ramp_data, tfx_color_ramp_t *ramp) {
	tfxU32 layer = ramp_data->color_ramp_count / 256;
	tfxU32 y = ramp_data->color_ramp_count % 256;
	if (ramp_data->color_ramp_bitmaps.size() <= layer) {
		tfx_bitmap_t bitmap = tfxCreateBitmap(tfxCOLOR_RAMP_WIDTH, 256, 4);
		ramp_data->color_ramp_bitmaps.push_back(bitmap);
	}
	tfx_bitmap_t &current_bitmap = ramp_data->color_ramp_bitmaps[layer];
	PlotColorRamp(&current_bitmap, ramp, y);
	tfxU32 ramp_id = tfxMakeColorRampIndex(layer, y);
	ramp_data->color_ramp_count++;
	return ramp_id;
}

void CreateColorRampBitmaps(tfx_library_t *library) {
	tfxU32 y = 0;
	library->color_ramps.color_ramp_bitmaps.clear();
	library->color_ramps.color_ramp_ids.Clear();
	library->color_ramps.color_ramp_count = 0;
	for (tfx_emitter_attributes_t &a : library->emitter_attributes) {
		for (int i = 0; i != 2; ++i) {
			tfx_color_ramp_t &ramp = a.overtime.color_ramps[i];
			tfxKey hash = tfxXXHash64::hash(ramp.colors, sizeof(tfx_rgba8_t) * tfxCOLOR_RAMP_WIDTH, 0);
			if (library->color_ramps.color_ramp_ids.ValidKey(hash)) {	//0 = main color, 1 = color hint
				a.overtime.color_ramp_bitmap_indexes[i] = library->color_ramps.color_ramp_ids.At(hash);
			}
			else {
				a.overtime.color_ramp_bitmap_indexes[i] = AddColorRampToBitmaps(&library->color_ramps, &ramp);
				library->color_ramps.color_ramp_ids.Insert(hash, a.overtime.color_ramp_bitmap_indexes[i]);
			}
		}
	}
}

void EditColorRampBitmap(tfx_library_t *library, tfx_overtime_attributes_t *a, tfxU32 ramp_id) {
	tfx_color_ramp_t &ramp = a->color_ramps[ramp_id];
	if (!tfxColorRampIsEdited(a->color_ramp_bitmap_indexes[ramp_id])) {
		a->color_ramp_bitmap_indexes[ramp_id] = AddColorRampToBitmaps(&library->color_ramps, &ramp);
		tfxFlagColorRampIDAsEdited(a->color_ramp_bitmap_indexes[ramp_id]);
	}
	else {
		tfxU32 layer = tfxColorRampLayer(a->color_ramp_bitmap_indexes[ramp_id]);
		tfxU32 index = tfxColorRampIndex(a->color_ramp_bitmap_indexes[ramp_id]);
		tfx_bitmap_t &bitmap = library->color_ramps.color_ramp_bitmaps[layer];
		PlotColorRamp(&bitmap, &ramp, index);
	}
}

void MaybeInsertColorRampBitmap(tfx_library_t *library, tfx_overtime_attributes_t *a, tfxU32 ramp_id) {
	tfxKey hash = tfxXXHash64::hash(a->color_ramps, sizeof(tfx_rgba8_t) * tfxCOLOR_RAMP_WIDTH, 0);
	if (library->color_ramps.color_ramp_ids.ValidKey(hash)) {
		a->color_ramp_bitmap_indexes[ramp_id] = library->color_ramps.color_ramp_ids.At(hash);
	} else {
		a->color_ramp_bitmap_indexes[ramp_id] = AddColorRampToBitmaps(&library->color_ramps, &a->color_ramps[ramp_id]);
		library->color_ramps.color_ramp_ids.Insert(hash, a->color_ramp_bitmap_indexes[ramp_id]);
	}
}

void CopyColorRampToAnimationManager(tfx_animation_manager_t *animation_manager, tfxU32 properties_index, tfx_color_ramp_t *ramp) {
	tfxKey hash = tfxXXHash64::hash(ramp->colors, sizeof(tfx_rgba8_t) * tfxCOLOR_RAMP_WIDTH, 0);
	if (animation_manager->color_ramps.color_ramp_ids.ValidKey(hash)) {
		animation_manager->emitter_properties[properties_index].color_ramp_index = animation_manager->color_ramps.color_ramp_ids.At(hash);
	} else {
		tfxU32 ramp_id = AddColorRampToBitmaps(&animation_manager->color_ramps, ramp);
		animation_manager->emitter_properties[properties_index].color_ramp_index = ramp_id;
		animation_manager->color_ramps.color_ramp_ids.Insert(hash, ramp_id);
	}
}

float LookupFastOvertime(tfx_graph_t *graph, float age, float lifetime) {
	tfxU32 frame = static_cast<tfxU32>((age / lifetime * graph->lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME);
	frame = std::min<tfxU32>(frame, graph->lookup.last_frame);
	return graph->lookup.values[frame];
}

float LookupFast(tfx_graph_t *graph, float frame) {
	if ((tfxU32)frame < graph->lookup.last_frame)
		return graph->lookup.values[(tfxU32)frame];
	return graph->lookup.values[graph->lookup.last_frame];
}

float LookupPreciseOvertime(tfx_graph_t *graph, float age, float life) {
	float lastv = 0;
	float lastf = 0;
	float p = 0;
	tfx_attribute_node_t *lastec = nullptr;
	for (tfxBucketLoop(graph->nodes)) {
		float frame = graph->nodes[i].frame * life;
		if (age < frame) {
			p = (age - lastf) / (frame - lastf);
			float bezier_value = GetBezierValue(lastec, &graph->nodes[i], p, graph->min.y, graph->max.y);
			if (bezier_value) {
				return bezier_value;
			}
			else {
				return lastv - p * (lastv - graph->nodes[i].value);
			}
		}
		lastv = graph->nodes[i].value;
		lastf = frame - 1;
		lastec = &graph->nodes[i];
	}
	return lastv;
}

float LookupPrecise(tfx_graph_t *graph, float age) {
	float lastv = 0;
	float lastf = 0;
	float p = 0;
	tfx_attribute_node_t *lastec = nullptr;
	for (tfxBucketLoop(graph->nodes)) {
		if (age < graph->nodes[i].frame) {
			p = (age - lastf) / (graph->nodes[i].frame - lastf);
			float bezier_value = GetBezierValue(lastec, &graph->nodes[i], p, graph->min.y, graph->max.y);
			if (bezier_value) {
				return bezier_value;
			}
			else {
				return lastv - p * (lastv - graph->nodes[i].value);
			}
		}
		lastv = graph->nodes[i].value;
		lastf = graph->nodes[i].frame - 1;
		lastec = &graph->nodes[i];
	}
	return lastv;
}

float GetRandomFast(tfx_graph_t *graph, float frame, tfx_random_t *random) {
	float value = 0;
	if ((tfxU32)frame < graph->lookup.last_frame)
		value = graph->lookup.values[(tfxU32)frame];
	value = graph->lookup.values[graph->lookup.last_frame];
	return RandomRange(random, value);
}

float GetRandomPrecise(tfx_graph_t *graph, float frame, tfx_random_t *random) {
	return GetGraphRandomValue(graph, frame, random);
}

float GetMaxLife(tfx_effect_emitter_t *e) {
	tfx_graph_t &life = *GetEffectGraphByType(e, tfxBase_life);
	tfx_graph_t &life_variation = *GetEffectGraphByType(e, tfxVariation_life);
	float templife = 0;
	float max_life = 0;
	float life_last_frame = GetGraphLastFrame(&life, 60.f);
	float life_variation_last_frame = GetGraphLastFrame(&life_variation, 60.f);
	float global_adjust = 1.f;
	if (life_last_frame + life_variation_last_frame > 0) {
		for (float f = 0; f < fmaxf(life_last_frame, life_variation_last_frame); ++f) {
			if (e->parent)
				global_adjust = GetGraphValue(GetEffectGraphByType(e->parent, tfxGlobal_life), f);
			templife = GetGraphValue(&life, f) + GetGraphValue(&life_variation, f);
			templife *= global_adjust;
			if (max_life < templife)
				max_life = templife;
		}
	}
	else {
		max_life = GetGraphFirstValue(&life) + GetGraphFirstValue(&life_variation);
	}

	return max_life;
}

bool IsOvertimeGraph(tfx_graph_type type) {
	return type >= tfxOvertime_velocity && type != tfxOvertime_noise_resolution && type <= tfxOvertime_noise_resolution;
}

bool IsColorGraph(tfx_graph_type type) {
	return type >= tfxOvertime_red && type <= tfxOvertime_blue;
}

bool IsOvertimePercentageGraph(tfx_graph_type type) {
	return type >= tfxOvertime_velocity && type != tfxOvertime_velocity_adjuster && type != tfxOvertime_direction && type <= tfxOvertime_noise_resolution;
}

bool IsGlobalGraph(tfx_graph_type type) {
	return type >= tfxGlobal_life && type <= tfxGlobal_emitter_depth;
}

bool IsEmitterGraph(tfx_graph_type type) {
	return type >= TFX_PROPERTY_START && type < TFX_TRANSFORM_START;
}

bool IsTransformGraph(tfx_graph_type type) {
	return type >= tfxTransform_roll && type <= tfxTransform_translate_z;
}

bool IsGlobalPercentageGraph(tfx_graph_type type) {
	return type >= tfxGlobal_life && type <= tfxGlobal_splatter;
}

bool IsEmitterSizeGraph(tfx_graph_type type) {
	return type >= tfxProperty_emitter_width && type <= tfxProperty_emitter_depth;
}

bool IsAngleGraph(tfx_graph_type type) {
	return (type == tfxTransform_roll || type == tfxTransform_pitch || type == tfxTransform_yaw || type == tfxProperty_emission_pitch || type == tfxProperty_emission_yaw || type == tfxProperty_emission_range ||
		type == tfxProperty_arc_offset || type == tfxProperty_arc_size || type == tfxBase_roll_spin || type == tfxVariation_roll_spin || type == tfxBase_pitch_spin || type == tfxVariation_pitch_spin ||
		type == tfxBase_yaw_spin || type == tfxVariation_yaw_spin);
}

bool IsAngleOvertimeGraph(tfx_graph_type type) {
	return type == tfxOvertime_direction;
}

bool IsEverythingElseGraph(tfx_graph_type type) {
	return !IsOvertimeGraph(type) && !IsOvertimePercentageGraph(type) && !IsGlobalGraph(type) && !IsAngleGraph(type) && !IsOvertimeGraph(type);
}

bool HasNodeAtFrame(tfx_graph_t *graph, float frame) {
	for (tfxBucketLoop(graph->nodes)) {
		if (graph->nodes[i].frame == frame) return true;
	}
	return false;
}

bool HasKeyframes(tfx_effect_emitter_t *e) {
	TFX_ASSERT(e->transform_attributes < e->library->transform_attributes.size());        //Must be a valid keyframes index into the library
	tfx_transform_attributes_t &keyframes = e->library->transform_attributes[e->transform_attributes];
	tfxU32 size = keyframes.translation_x.nodes.size() +
		keyframes.translation_y.nodes.size() +
		keyframes.translation_z.nodes.size();
	return size > 0;
}

bool HasMoreThanOneKeyframe(tfx_effect_emitter_t *e) {
	TFX_ASSERT(e->transform_attributes < e->library->transform_attributes.size());        //Must be a valid keyframes index into the library
	tfx_transform_attributes_t &keyframes = e->library->transform_attributes[e->transform_attributes];
	return    keyframes.translation_x.nodes.size() > 1 ||
		keyframes.translation_y.nodes.size() > 1 ||
		keyframes.translation_z.nodes.size() > 1;
}

void PushTranslationPoints(tfx_effect_emitter_t *e, tfx_vector_t<tfx_vec3_t> *points, float frame) {
	TFX_ASSERT(e->transform_attributes < e->library->transform_attributes.size());        //Must be a valid keyframes index into the library
	tfx_transform_attributes_t *keyframes = &e->library->transform_attributes[e->transform_attributes];
	tfx_vec3_t point(lookup_callback(&keyframes->translation_x, frame),
		lookup_callback(&keyframes->translation_y, frame),
		lookup_callback(&keyframes->translation_z, frame));
	points->push_back(point);
}

bool tfx__has_data_value(tfx_storage_map_t<tfx_data_entry_t> *config, tfx_str32_t key) {
	return config->ValidName(key);
}

void tfx__add_data_value_str(tfx_storage_map_t<tfx_data_entry_t> *config, tfx_str32_t key, const char *value) {
	tfx_data_entry_t entry;
	entry.type = tfxString;
	entry.key = key;
	entry.str_value = value;
	entry.str_value.SanitizeLineFeeds();
	config->Insert(key, entry);
}

void tfx__add_data_value_int(tfx_storage_map_t<tfx_data_entry_t> *config, tfx_str32_t key, int value) {
	tfx_data_entry_t entry;
	entry.type = tfxSInt;
	entry.key = key;
	entry.int_value = value;
	entry.bool_value = (bool)value;
	config->Insert(key, entry);
}

void tfx__add_color_value(tfx_storage_map_t<tfx_data_entry_t> *config, tfx_str32_t key, tfx_rgba8_t color) {
	tfx_data_entry_t entry;
	entry.type = tfxColor;
	entry.key = key;
	entry.int_value = (int)color.color;
	entry.color_value = color;
	entry.bool_value = (bool)color.color;
	config->Insert(key, entry);
}

void tfx__add_color_value_from_int(tfx_storage_map_t<tfx_data_entry_t> *config, tfx_str32_t key, tfxU32 color) {
	tfx_data_entry_t entry;
	entry.type = tfxColor;
	entry.key = key;
	entry.int_value = color;
	entry.color_value.color = color;
	entry.bool_value = (bool)color;
	config->Insert(key, entry);
}

void tfx__add_data_value_bool(tfx_storage_map_t<tfx_data_entry_t> *config, tfx_str32_t key, bool value) {
	tfx_data_entry_t entry;
	entry.type = tfxBool;
	entry.key = key;
	entry.bool_value = value;
	entry.int_value = (int)value;
	config->Insert(key, entry);
}

void tfx__add_data_value_double(tfx_storage_map_t<tfx_data_entry_t> *config, tfx_str32_t key, double value) {
	tfx_data_entry_t entry;
	entry.type = tfxDouble;
	entry.key = key;
	entry.double_value = value;
	config->Insert(key, entry);
}

void tfx__add_data_value_float(tfx_storage_map_t<tfx_data_entry_t> *config, tfx_str32_t key, float value) {
	tfx_data_entry_t entry;
	entry.type = tfxFloat;
	entry.key = key;
	entry.float_value = value;
	config->Insert(key, entry);
}

tfx_str_t tfx__get_data_str_value(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key) {
	return config->At(key).str_value;
}
int tfx__get_data_int_value(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key) {
	return config->At(key).int_value;
}
tfx_rgba8_t tfx__get_data_color_value(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key) {
	return config->At(key).color_value;
}
float tfx__get_data_float_value(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key) {
	return config->At(key).float_value;
}

bool tfx__save_data_file(tfx_storage_map_t<tfx_data_entry_t> *config, const char *path) {
	FILE *file = tfx__open_file(path, "wb");

	if (file == NULL)
		return false;

	if (config->Size()) {
		for (auto &entry : config->data) {
			tfx_str512_t ini_line = entry.key;
			ini_line.Appendf("=");
			switch (entry.type) {
			case tfxString:
				ini_line.Appendf(entry.str_value.c_str());
				break;
			case tfxSInt:
				ini_line.Appendf("%i", entry.int_value);
				break;
			case tfxColor:
				ini_line.Appendf("%ull", entry.color_value.color);
				break;
			case tfxFloat:
				ini_line.Appendf("%f", entry.float_value);
				break;
			case tfxBool:
				ini_line.Appendf("%i", (int)entry.bool_value);
				break;
			default:
				break;
			}
			ini_line.Appendf("\n");
			fwrite(ini_line.c_str(), sizeof(char), ini_line.Length(), file);
		}
	}

	fclose(file);
	return true;

}

bool tfx__load_data_file(tfx_data_types_dictionary_t *data_types, tfx_storage_map_t<tfx_data_entry_t> *config, const char *path) {
	FILE *fp = tfx__open_file(path, "rb");
	if (fp == NULL) {
		return false;
	}

	const size_t max_line_length = 512;
	char buffer[max_line_length];

	tfx_vector_t<tfx_str256_t> pair;
	while (fgets(buffer, max_line_length, fp)) {
		buffer[strcspn(buffer, "\n")] = 0;
		tfx_str512_t str = buffer;
		pair.clear();
		tfx__split_string_vec(str, &pair, 61);
		if (pair.size() == 2) {
			tfx_str256_t key = pair[0];
			if (data_types->names_and_types.ValidName(pair[0])) {
				tfx_data_type t = data_types->names_and_types.At(pair[0]);
				if (t == tfxBool) {
					tfx__add_data_value_bool(config, key, (bool)atoi(pair[1].c_str()));
				}
				if (t == tfxSInt) {
					tfx__add_data_value_int(config, key, atoi(pair[1].c_str()));
				}
				else if (t == tfxFloat) {
					tfx__add_data_value_float(config, key, (float)atof(pair[1].c_str()));
				}
				else if (t == tfxString) {
					tfx__add_data_value_str(config, key, pair[1].c_str());
				}
				else if (t == tfxColor) {
					char *endptr;
					tfx__add_color_value_from_int(config, key, (tfxU32)strtoul(pair[1].c_str(), &endptr, 10));
				}
			}
		}
	}

	fclose(fp);
	return true;

}

void tfx__split_string_stack(const tfx_str_t str, tfx_vector_t<tfx_str256_t> *pair, char delim) {
	tfx_str256_t line;
	for (char c : str) {
		if (c == delim && line.Length() && c != 0) {
			pair->push_back(line);
			line.Clear();
		}
		else if (c != 0) {
			line.Append(c);
		}
	}

	if (line.Length()) {
		pair->push_back(line);
	}
}

void tfx__split_string_vec(const tfx_str_t str, tfx_vector_t<tfx_str256_t> *pair, char delim) {
	tfx_str256_t line;
	for (char c : str) {
		if (c == delim && line.Length() && c != 0) {
			pair->push_back(line);
			line.Clear();
		}
		else if (c != 0) {
			line.Append(c);
		}
	}

	if (line.Length()) {
		pair->push_back(line);
	}
}

bool tfx__string_is_uint(const tfx_str_t s) {

	for (auto c : s) {
		if (!std::isdigit(c) && c != 0)
			return false;
	}

	return true;
}

int tfx__get_data_type(const tfx_str_t &s) {
	if (s.Length() == 0)
		return tfxString;

	if (s.IsInt())
		return tfxSInt;

	if (s.IsFloat())
		return tfxFloat;

	return tfxString;
}

//Get a graph by tfx_graph_id_t
tfx_graph_t *tfxGetGraph(tfx_library_t *library, tfx_graph_id_t graph_id) {
	tfx_graph_type type = graph_id.type;

	if (type < TFX_GLOBAL_COUNT) {
		return &((tfx_graph_t *)&library->global_graphs[graph_id.graph_id])[type];
	}
	else if (type >= TFX_PROPERTY_START && type < TFX_BASE_START) {
		int ref = type - TFX_PROPERTY_START;
		return &((tfx_graph_t *)&library->emitter_attributes[graph_id.graph_id].properties)[ref];
	}
	else if (type >= TFX_BASE_START && type < TFX_VARIATION_START) {
		int ref = type - TFX_BASE_START;
		return &((tfx_graph_t *)&library->emitter_attributes[graph_id.graph_id].base)[ref];
	}
	else if (type >= TFX_VARIATION_START && type < TFX_OVERTIME_START) {
		int ref = type - TFX_VARIATION_START;
		return &((tfx_graph_t *)&library->emitter_attributes[graph_id.graph_id].variation)[ref];
	}
	else if (type >= TFX_OVERTIME_START && type < TFX_FACTOR_START) {
		int ref = type - TFX_OVERTIME_START;
		return &((tfx_graph_t *)&library->emitter_attributes[graph_id.graph_id].overtime)[ref];
	}
	else if (type >= TFX_FACTOR_START && type < TFX_TRANSFORM_START) {
		int ref = type - TFX_FACTOR_START;
		return &((tfx_graph_t *)&library->emitter_attributes[graph_id.graph_id].factor)[ref];
	}
	else if (type >= TFX_TRANSFORM_START) {
		int ref = type - TFX_TRANSFORM_START;
		return &((tfx_graph_t *)&library->transform_attributes[graph_id.graph_id])[ref];
	}

	TFX_ASSERT(0);    //This function must return a value, make sure the graph_id is valid

	return nullptr;

}

//API Functions
//Get the number of shapes that are stored in an effects library saved on disk. This can be useful if you need to reserve the space in 
//a list to store them in your custom ShapeLoader function.
int GetShapeCountInLibrary(const char *filename) {
	int context = 0;
	int error = 0;

	tfx_package_t package;
	error = tfx__load_package_file(filename, &package);

	tfx_package_entry_info_t *data = tfx__get_package_file(&package, "data.txt");

	if (!data)
		error = -5;


	if (error < 0) {
		tfx__free_package(&package);
		return error;
	}

	int shape_count = 0;
	tmpStack(tfx_str256_t, pair);

	while (!data->data.EoF()) {
		pair.clear();
		tfx_str128_t line = data->data.ReadLine();
		bool context_set = false;
		if (tfx__string_is_uint(line.c_str())) {
			context = atoi(line.c_str());
			if (context == tfxEndShapes)
				break;
			context_set = true;
		}
		if (context_set == false) {
			tfx__split_string_stack(line.c_str(), &pair);
			if (pair.size() != 2) {
				pair.clear();
				tfx__split_string_stack(line.c_str(), &pair, 44);
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

	tfx__free_package(&package);
	return shape_count;
}

int GetEffectLibraryStats(const char *filename, tfx_effect_library_stats_t *stats) {
	int context = 0;
	int error = 0;

	tfx_package_t package;
	error = tfx__load_package_file(filename, &package);

	tfx_package_entry_info_t *data = tfx__get_package_file(&package, "data.txt");

	if (!data)
		error = -5;


	if (error < 0) {
		tfx__free_package(&package);
		return error;
	}

	memset(stats, 0, sizeof(tfx_effect_library_stats_t));
	bool inside_emitter = false;

	tmpStack(tfx_str256_t, pair);
	while (!data->data.EoF()) {
		pair.clear();
		tfx_str128_t line = data->data.ReadLine();
		bool context_set = false;
		if (tfx__string_is_uint(line.c_str())) {
			context_set = true;
			if (context == tfxEndEmitter) {
				inside_emitter = false;
			}
		}
		if (context_set == false) {
			tfx__split_string_stack(line.c_str(), &pair);
			if (pair.size() != 2) {
				pair.clear();
				tfx__split_string_stack(line.c_str(), &pair, 44);
				if (pair.size() < 2) {
					error = 1;
					break;
				}
			}
			if (context == tfxStartShapes) {
				if (pair.size() >= 5) {
					int frame_count = atoi(pair[2].c_str());
					stats->total_shapes += frame_count;
				}
			}
			else if (context == tfxStartEmitter) {
				inside_emitter = true;
				stats->total_emitters++;
			}
			else if (context == tfxStartEffect) {
				if (inside_emitter)
					stats->total_sub_effects++;
				else
					stats->total_effects++;
			}
		}
	}

	return error;
}

tfx_effect_library_stats_t CreateLibraryStats(tfx_library_t *lib) {
	tfx_effect_library_stats_t stats;
	memset(&stats, 0, sizeof(stats));
	stats.total_effects = lib->effects.size();
	stats.total_node_lookup_indexes = lib->node_lookup_indexes.size();
	stats.total_attribute_nodes = lib->all_nodes.size();
	tmpStack(tfx_effect_emitter_t, stack);
	for (auto &effect : lib->effects) {
		stack.push_back(effect);
	}
	while (!stack.empty()) {
		tfx_effect_emitter_t &current = stack.pop_back();
		if (current.parent) {
			if (current.type == tfxEffectType) {
				stats.total_sub_effects++;
			}
			else if (current.type == tfxEmitterType) {
				stats.total_emitters++;
			}
		}
		for (auto &sub : GetEffectInfo(&current)->sub_effectors) {
			stack.push_back(sub);
		}
	}
	stats.total_shapes = lib->particle_shapes.data.size();

	return stats;
}

tfxAPI tfxErrorFlags LoadSpriteData(const char *filename, tfx_animation_manager_t *animation_manager, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void *user_data) {
	//TFX_ASSERT(shape_loader);            //Must have a shape_loader function to load your shapes with. This will be a custom user function suited for whichever renderer you're using
	if (!tfxStore->data_types.initialised)
		tfxStore->data_types.Init();

	tfx_package_t package;
	tfxErrorFlags error = tfx__load_package_file(filename, &package);
	if (error != 0) {
		tfx__free_package(&package);
		return error;
	}

	tfx_package_entry_info_t *data = tfx__get_package_file(&package, "data.txt");

	if (!data) {
		error |= tfxErrorCode_data_could_not_be_loaded;
	}

	tfx_package_entry_info_t *sprite_data = tfx__get_package_file(&package, "sprite_data");

	if (!sprite_data) {
		error |= tfxErrorCode_data_could_not_be_loaded;
	}

	if (error != 0) {
		tfx__free_package(&package);
		return error;
	}

	if (package.header.user_data2 == 1 && !(animation_manager->flags & tfxAnimationManagerFlags_is_3d)) {
		return tfxErrorCode_sprite_data_is_3d_but_animation_manager_is_2d;
	}

	if (package.header.user_data2 == 0 && animation_manager->flags & tfxAnimationManagerFlags_is_3d) {
		return tfxErrorCode_sprite_data_is_2d_but_animation_manager_is_3d;
	}

	if (animation_manager->flags & tfxAnimationManagerFlags_is_3d) {
		animation_manager->sprite_data_3d.resize((tfxU32)(sprite_data->file_size / package.header.user_data1));
		memcpy(animation_manager->sprite_data_3d.data, sprite_data->data.data, sprite_data->file_size);
	}
	else {
		animation_manager->sprite_data_2d.resize((tfxU32)(sprite_data->file_size / package.header.user_data1));
		memcpy(animation_manager->sprite_data_2d.data, sprite_data->data.data, sprite_data->file_size);
	}

	tmpStack(tfx_sprite_data_metrics_t, metrics_stack);
	tmpStack(tfx_frame_meta_t, frame_meta_stack);
	tmpStack(tfx_animation_emitter_properties_t, emitter_properties_stack);
	tmpStack(tfx_str256_t, pair);
	tmpStack(tfx_str256_t, multi);

	tfxKey first_shape_hash = 0;
	int context = 0;

	while (!data->data.EoF()) {
		tfx_str512_t line = data->data.ReadLine();
		bool context_set = false;

		if (tfx__string_is_uint(line.c_str())) {
			context = atoi(line.c_str());
			if (context == tfxEndOfFile) {
				break;
			}

			context_set = true;
			if (context == tfxStartEffectAnimationInfo) {
				tfx_sprite_data_metrics_t metrics;
				metrics_stack.push_back_copy(metrics);
			}
			else if (context == tfxStartFrameMeta) {
				tfx_frame_meta_t frame_meta;
				frame_meta_stack.push_back_copy(frame_meta);
			}
			else if (context == tfxStartEmitter) {
				tfx_animation_emitter_properties_t emitter_properties;
				emitter_properties_stack.push_back_copy(emitter_properties);
			}
		}

		if (context_set == false) {
			pair.clear();
			tfx__split_string_stack(line.c_str(), &pair);
			if (pair.size() != 2) {
				pair.clear();
				tfx__split_string_stack(line.c_str(), &pair, ',');
				if (pair.size() < 2) {
					error |= tfxErrorCode_some_data_not_loaded;
					continue;
				}
			}

			if (context == tfxStartEffectAnimationInfo) {
				if (tfxStore->data_types.names_and_types.ValidName(pair[0])) {
					switch (tfxStore->data_types.names_and_types.At(pair[0])) {
					case tfxUint:
						tfx__assign_sprite_data_metrics_property_u32(&metrics_stack.back(), &pair[0], (tfxU32)atoi(pair[1].c_str()), package.header.file_version);
						break;
					case tfxFloat:
						tfx__assign_sprite_data_metrics_property_float(&metrics_stack.back(), &pair[0], (float)atof(pair[1].c_str()), package.header.file_version);
						break;
					case tfxString:
						tfx__assign_sprite_data_metrics_property_str(&metrics_stack.back(), &pair[0], pair[1], package.header.file_version);
						break;
					default:
						break;
					}
				}
				else {
					error |= tfxErrorCode_some_data_not_loaded;
				}
			}
			else if (context == tfxStartFrameMeta) {
				if (tfxStore->data_types.names_and_types.ValidName(pair[0])) {
					switch (tfxStore->data_types.names_and_types.At(pair[0])) {
					case tfxUint:
						tfx__assign_frame_meta_property_u32(&frame_meta_stack.back(), &pair[0], (tfxU32)atoi(pair[1].c_str()), package.header.file_version);
						break;
					case tfxFloat3:
						multi.clear();
						tfx__split_string_stack(pair[1], &multi, ',');
						tfx__assign_frame_meta_property_vec3(&frame_meta_stack.back(), &pair[0], tfx__str_to_vec3(&multi), package.header.file_version);
						break;
					default:
						break;
					}
				}
			}
			else if (context == tfxStartEmitter) {
				if (tfxStore->data_types.names_and_types.ValidName(pair[0])) {
					switch (tfxStore->data_types.names_and_types.At(pair[0])) {
					case tfxUint:
						tfx__assign_animation_emitter_property_u32(&emitter_properties_stack.back(), &pair[0], (tfxU32)atoi(pair[1].c_str()), package.header.file_version);
						break;
					case tfxFloat:
						tfx__assign_animation_emitter_property_float(&emitter_properties_stack.back(), &pair[0], (float)atof(pair[1].c_str()), package.header.file_version);
						break;
					case tfxFloat2:
						multi.clear();
						tfx__split_string_stack(pair[1], &multi, ',');
						tfx__assign_animation_emitter_property_vec2(&emitter_properties_stack.back(), &pair[0], tfx__str_to_vec2(&multi), package.header.file_version);
						break;
					default:
						break;
					}
				}
			}
			else if (context == tfxStartFrameOffsets) {
				if (pair.size() == tfxLAYERS + 1) {
					if (pair[0] == "index_offset") {
						for (int i = 1; i != tfxLAYERS + 1; ++i) {
							frame_meta_stack.back().index_offset[i - 1] = atoi(pair[i].c_str());
						}
					}
					else if (pair[0] == "sprite_count") {
						for (int i = 1; i != tfxLAYERS + 1; ++i) {
							frame_meta_stack.back().sprite_count[i - 1] = atoi(pair[i].c_str());
						}
					}
				}
				else {
					//Not enougn layers, set tfxLAYERS to the required amount. The default is 4.
				}
			}

			if (context == tfxStartShapes && shape_loader != nullptr) {
				if (pair.size() >= 5) {
					tfx_shape_data_t s;
					tfx__strcpy(s.name, pair[0].c_str());
					s.shape_index = atoi(pair[1].c_str());
					s.frame_count = atoi(pair[2].c_str());
					s.width = atoi(pair[3].c_str());
					s.height = atoi(pair[4].c_str());
					s.import_filter = atoi(pair[5].c_str());
					if (pair.current_size > 6) {
						s.image_hash = strtoull(pair[6].c_str(), NULL, 10);
					}
					if (s.import_filter < 0 || s.import_filter>1) {
						s.import_filter = 0;
					}

					tfx_package_entry_info_t *shape_entry = tfx__get_package_file(&package, s.name);
					if (shape_entry) {
						tfx_image_data_t image_data;
						image_data.shape_index = s.shape_index;
						image_data.animation_frames = (float)s.frame_count;
						image_data.image_size = tfx_vec2_t((float)s.width, (float)s.height);
						image_data.name = s.name;
						image_data.import_filter = s.import_filter;
						image_data.image_hash = tfxXXHash64::hash(shape_entry->data.data, shape_entry->file_size, 0);
						if (s.image_hash == 0) {
							s.image_hash = image_data.image_hash;
						}
						TFX_ASSERT(s.image_hash == image_data.image_hash);

						shape_loader(s.name, &image_data, shape_entry->data.data, (tfxU32)shape_entry->file_size, user_data);

						if (!image_data.ptr) {
							//uid = -6;
						}
						else {
							animation_manager->particle_shapes.Insert(image_data.image_hash, image_data);
							if (first_shape_hash == 0)
								first_shape_hash = s.image_hash;
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

		if (context == tfxEndFrameOffsets) {
			context = tfxStartFrameMeta;
		}
		else if (context == tfxEndFrameMeta) {
			TFX_ASSERT(metrics_stack.current_size);
			TFX_ASSERT(frame_meta_stack.current_size);
			metrics_stack.back().frame_meta.push_back_copy(frame_meta_stack.pop_back());
		}
		else if (context == tfxEndEffectAnimationInfo) {
			TFX_ASSERT(metrics_stack.current_size);
			animation_manager->effect_animation_info.Insert(metrics_stack.back().name, metrics_stack.back());
			metrics_stack.pop();
		}
		else if (context == tfxEndEmitter) {
			TFX_ASSERT(emitter_properties_stack.current_size);
			emitter_properties_stack.back().handle_packed = tfx__pack16bit_sscaled(emitter_properties_stack.back().handle.x, emitter_properties_stack.back().handle.y, 128.f);
			animation_manager->emitter_properties.push_back_copy(emitter_properties_stack.pop_back());
		}

	}

	tfxU32 bitmap_index = 0;
	tfx_str_t bitmap_file_name = "color_ramp_0";

	while (tfx__file_exists_in_package(&package, bitmap_file_name.c_str())) {
		tfx_package_entry_info_t *bitmap_file = tfx__get_package_file(&package, bitmap_file_name.c_str());
		tfx_bitmap_t bitmap = tfxCreateBitmap(tfxCOLOR_RAMP_WIDTH, tfxCOLOR_RAMP_WIDTH, 4);
		memcpy(bitmap.data, bitmap_file->data.data, bitmap.size);
		bitmap_file_name.Clear();
		bitmap_index++;
		bitmap_file_name.Appendf("color_ramp_%u", bitmap_index);
		animation_manager->color_ramps.color_ramp_bitmaps.push_back(bitmap);
		animation_manager->color_ramps.color_ramp_count++;
	}

	tfx__free_package(&package);

	if (!shape_loader) {
		error |= tfxErrorCode_library_loaded_without_shape_loader;
	}

	return error;
}

tfxErrorFlags LoadEffectLibraryPackage(tfx_package_t *package, tfx_library_t *lib, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset), void *user_data) {

	TFX_ASSERT(shape_loader);            //Must have a shape_loader function to load your shapes with. This will be a custom user function suited for whichever renderer you're using
	if (!tfxStore->data_types.initialised)
		tfxStore->data_types.Init();

	ClearLibrary(lib);
	if (tfxIcospherePoints[0].current_size == 0) {
		tfx__make_icospheres();
	}

	tfx_package_entry_info_t *data = tfx__get_package_file(package, "data.txt");
	tfx_package_entry_info_t *stats_struct = tfx__get_package_file(package, "stats.struct");
	tfxErrorFlags error = 0;

	int context = 0;
	int uid = 0;
	tfxU32 current_global_graph = 0;

	if (!data)
		error |= tfxErrorCode_data_could_not_be_loaded;

	if (error != 0) {
		tfx__free_package(package);
		return error;
	}

	tfxKey first_shape_hash = 0;

	//You must call InitialiseTimelineFX() before doing anything!    
	tmpStack(tfx_effect_emitter_t, effect_stack);
	tmpStack(tfx_str256_t, pair);

	while (!data->data.EoF()) {
		tfx_str512_t line = data->data.ReadLine();
		bool context_set = false;

		if (tfx__string_is_uint(line.c_str())) {
			context = atoi(line.c_str());
			if (context == tfxEndOfFile)
				break;

			context_set = true;
			if (context == tfxStartFolder) {
				tfx_effect_emitter_t effect{};
				effect.library = lib;
				effect.type = tfx_effect_emitter_type::tfxFolder;
				effect.info_index = AddLibraryEffectEmitterInfo(lib);
				GetEffectInfo(&effect)->uid = uid++;
				effect_stack.push_back(effect);
			}
			else if (context == tfxStartStage) {
				tfx_effect_emitter_t effect{};
				effect.library = lib;
				effect.type = tfx_effect_emitter_type::tfxStage;
				effect.info_index = AddLibraryEffectEmitterInfo(lib);
				AddLibraryPreviewCameraSettings(lib, &effect);
				effect.transform_attributes = AddLibraryKeyframes(lib);
				GetEffectInfo(&effect)->uid = uid++;
				effect_stack.push_back(effect);
			}
			else if (context == tfxStartEffect) {
				tfx_effect_emitter_t effect{};
				effect.library = lib;
				effect.info_index = AddLibraryEffectEmitterInfo(lib);
				effect.property_index = AddLibraryEmitterProperties(lib);
				effect.transform_attributes = AddLibraryKeyframes(lib);
				if (effect_stack.size() <= 1) { //Only root effects get the global graphs
					AddLibraryEffectGraphs(lib, &effect);
					ResetEffectGraphs(&effect, false, false);
					current_global_graph = effect.global;
				}
				AddLibraryTransformGraphs(lib, &effect);
				ResetTransformGraphs(&effect, false, false);
				effect.type = tfx_effect_emitter_type::tfxEffectType;
				AddLibrarySpriteSheetSettings(lib, &effect);
				AddLibrarySpriteDataSettings(lib, &effect);
				AddLibraryPreviewCameraSettings(lib, &effect);
				GetEffectInfo(&effect)->uid = uid++;
				effect_stack.push_back(effect);

			}
			else if (context == tfxStartEmitter) {
				tfx_effect_emitter_t emitter = {};
				emitter.effect_flags = 0;
				emitter.property_flags = 0;
				emitter.library = lib;
				emitter.info_index = AddLibraryEffectEmitterInfo(lib);
				emitter.property_index = AddLibraryEmitterProperties(lib);
				emitter.transform_attributes = AddLibraryKeyframes(lib);
				AddLibraryEmitterGraphs(lib, &emitter);
				AddLibraryTransformGraphs(lib, &emitter);
				emitter.type = tfx_effect_emitter_type::tfxEmitterType;
				ResetEmitterGraphs(&emitter, false, false);
				ResetTransformGraphs(&emitter, false, false);
				GetEffectInfo(&emitter)->uid = uid++;
				effect_stack.push_back(emitter);
			}
		}

		if (context_set == false) {
			pair.clear();
			tfx__split_string_stack(line.c_str(), &pair);
			if (pair.size() != 2) {
				pair.clear();
				tfx__split_string_stack(line.c_str(), &pair, 44);
				if (pair.size() < 2) {
					error |= tfxErrorCode_some_data_not_loaded;
					continue;
				}
			}

			if (context == tfxStartAnimationSettings || context == tfxStartEmitter || context == tfxStartEffect || context == tfxStartFolder || context == tfxStartPreviewCameraSettings) {
				if (tfxStore->data_types.names_and_types.ValidName(pair[0])) {
					tfx__assign_property_line(&effect_stack.back(), &pair, package->header.file_version);
				}
				else {
					error |= tfxErrorCode_some_data_not_loaded;
				}
			}
			else if (context == tfxStartGraphs && effect_stack.back().type == tfxEmitterType) {
				tfx__assign_graph_data(&effect_stack.back(), &pair);
			}
			else if (context == tfxStartGraphs && effect_stack.back().type == tfxEffectType) {
				if (effect_stack.size() <= 2)
					tfx__assign_graph_data(&effect_stack.back(), &pair);
			}
			else if (context == tfxStartStage) {
				if (tfxStore->data_types.names_and_types.ValidName(pair[0])) {
					switch (tfxStore->data_types.names_and_types.At(pair[0])) {
					case tfxUint:
						tfx__assign_stage_property_u32(&effect_stack.back(), &pair[0], (tfxU32)atoi(pair[1].c_str()));
						break;
					case tfxFloat:
						tfx__assign_stage_property_float(&effect_stack.back(), &pair[0], (float)atof(pair[1].c_str()));
						break;
					case tfxSInt:
						tfx__assign_stage_property_int(&effect_stack.back(), &pair[0], atoi(pair[1].c_str()));
						break;
					case tfxBool:
						tfx__assign_stage_property_bool(&effect_stack.back(), &pair[0], (bool)(atoi(pair[1].c_str())));
						break;
					case tfxString:
						tfx__assign_stage_property_str(&effect_stack.back(), &pair[0], &pair[1]);
						break;
					default:
						break;
					}
				}
				else {
					error |= tfxErrorCode_some_data_not_loaded;
				}
			}

			if (context == tfxStartShapes) {
				if (pair.size() >= 5) {
					tfx_shape_data_t s;
					tfx__strcpy(s.name, pair[0].c_str());
					s.shape_index = atoi(pair[1].c_str());
					s.frame_count = atoi(pair[2].c_str());
					s.width = atoi(pair[3].c_str());
					s.height = atoi(pair[4].c_str());
					s.import_filter = atoi(pair[5].c_str());
					if (pair.current_size > 6) {
						s.image_hash = strtoull(pair[6].c_str(), NULL, 10);
					}
					if (s.import_filter < 0 || s.import_filter>1) {
						s.import_filter = 0;
					}

					tfx_package_entry_info_t *shape_entry = tfx__get_package_file(package, s.name);
					if (shape_entry) {
						tfx_image_data_t image_data;
						image_data.shape_index = s.shape_index;
						image_data.animation_frames = (float)s.frame_count;
						image_data.image_size = tfx_vec2_t((float)s.width, (float)s.height);
						image_data.name = s.name;
						image_data.import_filter = s.import_filter;
						image_data.image_hash = tfxXXHash64::hash(shape_entry->data.data, shape_entry->file_size, 0);
						if (s.image_hash == 0) {
							s.image_hash = image_data.image_hash;
						}
						TFX_ASSERT(s.image_hash == image_data.image_hash);

						shape_loader(s.name, &image_data, shape_entry->data.data, (tfxU32)shape_entry->file_size, user_data);

						if (!image_data.ptr) {
							//uid = -6;
						}
						else {
							lib->particle_shapes.Insert(image_data.image_hash, image_data);
							if (first_shape_hash == 0)
								first_shape_hash = s.image_hash;
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
			InitialiseUninitialisedGraphs(&effect_stack.back());
			UpdateEffectMaxLife(&effect_stack.back());
			UpdateEmitterControlProfile(&effect_stack.back());
			if (effect_stack.back().property_flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
				lib->emitter_properties[effect_stack.back().property_index].image_handle = { .5f, .5f };
			}
			tfx_vec2_t handle = lib->emitter_properties[effect_stack.back().property_index].image_handle;
			lib->emitter_properties[effect_stack.back().property_index].image_handle_packed = tfx__pack16bit_sscaled(handle.x, handle.y, 128.f);
			effect_stack.back().property_flags |= tfxEmitterPropertyFlags_enabled;
			if (effect_stack.back().path_attributes != tfxINVALID) {
				if (Is3DEffect(&effect_stack.parent())) {
					GetEmitterPath(&effect_stack.back())->flags &= ~tfxPathFlags_2d;
				}
			}
			GetEffectInfo(&effect_stack.parent())->sub_effectors.push_back(effect_stack.back());
			effect_stack.pop();
		}

		if (context == tfxEndEffect) {
			ReIndexEffect(&effect_stack.back());
			if (effect_stack.size() > 1) {
				if (effect_stack.parent().type == tfxStage && GetEffectInfo(&effect_stack.parent())->sub_effectors.size() == 0) {
					tfxEffectPropertyFlags tmp = effect_stack.parent().property_flags;
					if (Is3DEffect(&effect_stack.back())) {
						effect_stack.parent().property_flags |= tfxEmitterPropertyFlags_effect_is_3d;
					}
					tmp = effect_stack.parent().property_flags;
				}
				else if (effect_stack.parent().type == tfxEmitterType) {
					effect_stack.back().global = current_global_graph;
				}
				if (effect_stack.parent().type == tfxFolder) {
					InitialiseUninitialisedGraphs(&effect_stack.back());
				}
				GetEffectInfo(&effect_stack.parent())->sub_effectors.push_back(effect_stack.back());
			}
			else {
				InitialiseUninitialisedGraphs(&effect_stack.back());
				lib->effects.push_back(effect_stack.back());
			}
			effect_stack.pop();
		}

		if (context == tfxEndFolder) {
			TFX_ASSERT(effect_stack.size() == 1);            //Folders should not be contained within anything
			lib->effects.push_back(effect_stack.back());
			effect_stack.pop();
		}

		if (context == tfxEndStage) {
			TFX_ASSERT(effect_stack.size() == 1);            //Stages should not be contained within anything
			lib->effects.push_back(effect_stack.back());
			effect_stack.pop();
		}

	}

	tfx__free_package(package);

	if (uid >= 0) {
		//Effects were loaded so let's compile them
		CompileAllLibraryGraphs(lib);
		ReIndexLibrary(lib);
		if (first_shape_hash != 0) {
			UpdateLibraryParticleShapeReferences(lib, first_shape_hash);
		}
		UpdateLibraryEffectPaths(lib);
		UpdateLibraryComputeNodes(lib);
		BuildAllLibraryPaths(lib);
		SetLibraryMinMaxData(lib);
		lib->uv_lookup = uv_lookup;
		UpdateLibraryGPUImageData(lib);
	}
	lib->uid = uid;

	return error;
}

tfxErrorFlags LoadEffectLibrary(const char *filename, tfx_library_t *lib, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset), void *user_data) {

	tfxErrorFlags error = 0;

	tfx_package_t package;
	error = tfx__load_package_file(filename, &package);
	if (error != 0) {
		tfx__free_package(&package);
		return error;
	}
	error = LoadEffectLibraryPackage(&package, lib, shape_loader, uv_lookup, user_data);

	tfx__free_package(&package);
	return error;
}

tfxErrorFlags LoadEffectLibrary(const void *data, tfxU32 size, tfx_library_t *lib, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset), void *user_data) {
	tfxErrorFlags error = 0;

	tfx_stream_t data_stream;
	tfx__copy_data_to_stream(&data_stream, data, size);
	tfx_package_t package;
	error = tfx__load_package_stream(&data_stream, &package);
	if (error != 0) {
		tfx__free_package(&package);
		return error;
	}
	error = LoadEffectLibraryPackage(&package, lib, shape_loader, uv_lookup, user_data);

	tfx__free_package(&package);
	return error;
}

void SetTemplateUserDataAll(tfx_effect_template_t *t, void *data) {
	tmpStack(tfx_effect_emitter_t *, stack);
	stack.push_back(&t->effect);
	while (stack.size()) {
		tfx_effect_emitter_t *current = stack.pop_back();
		current->user_data = data;
		for (auto &sub : GetEffectInfo(current)->sub_effectors) {
			stack.push_back(&sub);
		}
	}
}

void ResetSpriteDataLerpOffset(tfx_sprite_data_t *sprite_data) {
	tfx_sprite_data_soa_t &sprites = sprite_data->real_time_sprites;
	for (int i = 0; i != sprite_data->real_time_sprites_buffer.current_size; ++i) {
		sprites.lerp_offset[i] = 1.f;
	}
}

void RecordSpriteData(tfx_particle_manager_t *pm, tfx_effect_emitter_t *effect, float update_frequency, float camera_position[3], std::atomic_int *progress) {
	TFX_ASSERT(update_frequency > 0); //Update frequency must be greater then 0. 60 is recommended for best results
	tfx_sprite_data_settings_t &anim = effect->library->sprite_data_settings[GetEffectInfo(effect)->sprite_data_settings_index];
	float frame_length = 1000.f / update_frequency;
	tfxU32 frames = anim.real_frames;
	tfxU32 start_frame = anim.frame_offset;
	int extra_frames = anim.extra_frames_count;
	tfxU32 frame = 0;
	int extra_frame_count = 0;
	tfxU32 offset = 0;
	bool particles_started = false;
	bool start_counting_extra_frames = false;
	pm->animation_length_in_time = (float)frames * frame_length - 1.f;
	bool is_3d = Is3DEffect(effect);
	progress->store(tfxCalculateFrames);

	anim.recording_frame_rate = tfxMin(tfxMax(30.f, anim.recording_frame_rate), 240.f);

	bool auto_set_length = false;
	if (anim.animation_flags & tfxAnimationFlags_auto_set_length && !(anim.animation_flags & tfxAnimationFlags_loop) && IsFiniteEffect(effect)) {
		frames = 99999;
		auto_set_length = true;
	}

	//First pass to count the number of instance_data in each frame
	tfxU32 total_sprites = 0;
	tfx_vector_t<tfx_frame_meta_t> tmp_frame_meta;
	tfxU32 sprites_in_layers = 0;
	tfx_vec3_t pm_camera_position = pm->camera_position;
	tfxEffectID preview_effect_index;

	//For now, record sprite data with instance_data stored in the particle manager rather then each effect.
	tfxParticleManagerFlags pm_effect_sprites_flag = pm->flags & tfxParticleManagerFlags_use_effect_sprite_buffers;
	pm->flags &= ~tfxParticleManagerFlags_use_effect_sprite_buffers;
	ReconfigureParticleManager(pm, GetRequiredParticleManagerMode(effect), effect->sort_passes, Is3DEffect(effect));
	tfxU32 struct_size = sizeof(tfx_2d_instance_t);
	if (is_3d) struct_size = sizeof(tfx_3d_instance_t);
	for (tfxEachLayer) {
		pm->instance_buffer_for_recording[0][layer].free();
		pm->instance_buffer_for_recording[0][layer] = tfxCreateBuffer(struct_size, 16);
		if (pm->flags & tfxParticleManagerFlags_double_buffer_sprites) {
			pm->instance_buffer_for_recording[1][layer].free();
			pm->instance_buffer_for_recording[1][layer] = tfxCreateBuffer(struct_size, 16);
			
		}
	}
	pm->flags |= tfxParticleManagerFlags_recording_sprites;
	if (!(pm->flags & tfxParticleManagerFlags_using_uids)) {
		ToggleSpritesWithUID(pm, true);
	}
	pm->unique_particle_id = 0;
	SetSeed(pm, anim.seed);
	preview_effect_index = AddEffectToParticleManager(pm, effect, pm->current_ebuff, 0, false, 0, 0.f);
	pm->camera_position = tfx_vec3_t(camera_position[0], camera_position[1], camera_position[2]);
	SetEffectPosition(pm, preview_effect_index, tfx_vec3_t(0.f, 0.f, 0.f));
	if (is_3d) {
		tfx__transform_3d(&pm->effects[preview_effect_index].world_rotations,
			&pm->effects[preview_effect_index].local_rotations,
			&pm->effects[preview_effect_index].overal_scale,
			&pm->effects[preview_effect_index].world_position,
			&pm->effects[preview_effect_index].local_position,
			&pm->effects[preview_effect_index].translation,
			&pm->effects[preview_effect_index].rotation,
			&pm->effects[preview_effect_index]
		);
	}
	else {
		tfx__transform_2d(&pm->effects[preview_effect_index].world_rotations,
			&pm->effects[preview_effect_index].local_rotations,
			&pm->effects[preview_effect_index].overal_scale,
			&pm->effects[preview_effect_index].world_position,
			&pm->effects[preview_effect_index].local_position,
			&pm->effects[preview_effect_index].translation,
			&pm->effects[preview_effect_index].rotation,
			&pm->effects[preview_effect_index]
		);
	}

	frame = 0;
	total_sprites = 0;
	sprites_in_layers = 0;
	tmp_frame_meta.clear();
	while (frame < frames &&offset < 99999) {
		tfxU32 count_this_frame = 0;
		UpdateParticleManager(pm, frame_length);
		bool particles_processed_last_frame = false;

		if (offset >= start_frame) {
			sprites_in_layers = 0;
			for (tfxEachLayer) {
				if (frame >= tmp_frame_meta.size()) {
					tfx_frame_meta_t meta;
					memset(&meta, 0, sizeof(tfx_frame_meta_t));
					tmp_frame_meta.push_back(meta);
				}
				tmp_frame_meta[frame].sprite_count[layer] += pm->layer_sizes[layer];
				tmp_frame_meta[frame].total_sprites += pm->layer_sizes[layer];
				total_sprites += pm->layer_sizes[layer];
				sprites_in_layers += pm->layer_sizes[layer];
				particles_started = total_sprites > 0;
				particles_processed_last_frame |= pm->layer_sizes[layer] > 0;
			}
		}

		if (auto_set_length && !(anim.animation_flags & tfxAnimationFlags_loop) && particles_started && sprites_in_layers == 0) {
			break;
		}

		offset++;

		if (particles_started)
			frame++;

		if (anim.animation_flags & tfxAnimationFlags_loop && particles_started) {
			if (frame == frames) {
				frame = 0;
				start_counting_extra_frames = true;
			}
			if (!particles_processed_last_frame)
				break;
		}
		if (start_counting_extra_frames && extra_frame_count++ >= extra_frames) {
			SoftExpireEffect(pm, preview_effect_index);
		}

	}

	progress->store(tfxBakeSpriteData);

	frames = tmp_frame_meta.size();

	tfx_sprite_data_t *sprite_data = nullptr;
	if (effect->library->pre_recorded_effects.ValidKey(effect->path_hash)) {
		sprite_data = &effect->library->pre_recorded_effects.At(effect->path_hash);
		FreeSpriteData(sprite_data);
		sprite_data->normal.frame_count = frames;
		sprite_data->normal.animation_length_in_time = frames * frame_length;
		sprite_data->normal.frame_meta.resize(frames);
		sprite_data->normal.frame_meta.zero();
	}
	else {
		tfx_sprite_data_t data;
		data.normal.frame_count = frames;
		data.normal.animation_length_in_time = frames * frame_length;
		data.normal.frame_meta.resize(frames);
		data.normal.frame_meta.zero();
		effect->library->pre_recorded_effects.Insert(effect->path_hash, data);
		sprite_data = &effect->library->pre_recorded_effects.At(effect->path_hash);
	}

	anim.real_frames = frames;
	anim.animation_length_in_time = sprite_data->normal.animation_length_in_time;
	sprite_data->frame_compression = anim.playback_speed;

	tfx_vector_t<tfx_frame_meta_t> &frame_meta = sprite_data->normal.frame_meta;
	memcpy(frame_meta.data, tmp_frame_meta.data, tmp_frame_meta.size_in_bytes());
	tmp_frame_meta.free_all();

	tfxU32 last_count = 0;
	for (tfx_frame_meta_t &meta : frame_meta) {
		meta.captured_offset = last_count;
		tfxU32 layer_counts = 0;
		for (tfxEachLayer) {
			meta.index_offset[layer] = last_count;
			meta.cumulative_offset[layer] = 0;
			last_count += meta.sprite_count[layer];
			layer_counts = last_count;
		}
	}

	ReconfigureParticleManager(pm, GetRequiredParticleManagerMode(effect), effect->sort_passes, Is3DEffect(effect));
	pm->flags |= tfxParticleManagerFlags_recording_sprites;
	if (!(pm->flags & tfxParticleManagerFlags_using_uids)) {
		ToggleSpritesWithUID(pm, true);
	}
	SetSeed(pm, anim.seed);
	preview_effect_index = AddEffectToParticleManager(pm, effect, pm->current_ebuff, 0, false, 0, 0.f);
	SetEffectPosition(pm, preview_effect_index, tfx_vec3_t(0.f, 0.f, 0.f));
	if (is_3d) {
		tfx__transform_3d(&pm->effects[preview_effect_index].world_rotations,
			&pm->effects[preview_effect_index].local_rotations,
			&pm->effects[preview_effect_index].overal_scale,
			&pm->effects[preview_effect_index].world_position,
			&pm->effects[preview_effect_index].local_position,
			&pm->effects[preview_effect_index].translation,
			&pm->effects[preview_effect_index].rotation,
			&pm->effects[preview_effect_index]
		);
	}
	else {
		tfx__transform_2d(&pm->effects[preview_effect_index].world_rotations,
			&pm->effects[preview_effect_index].local_rotations,
			&pm->effects[preview_effect_index].overal_scale,
			&pm->effects[preview_effect_index].world_position,
			&pm->effects[preview_effect_index].local_position,
			&pm->effects[preview_effect_index].translation,
			&pm->effects[preview_effect_index].rotation,
			&pm->effects[preview_effect_index]
		);
	}

	if (total_sprites == 0) {
		return;
	}

	sprite_data->normal.total_sprites = total_sprites;
	tfx_soa_buffer_t temp_sprites_buffer[tfxLAYERS];
	tfx_sprite_data_soa_t temp_sprites[tfxLAYERS];
	if (is_3d) {
		InitSpriteData3dSoA(&sprite_data->real_time_sprites_buffer, &sprite_data->real_time_sprites, total_sprites);
		for (tfxEachLayer) {
			InitSpriteData3dSoA(&temp_sprites_buffer[layer], &temp_sprites[layer], 100);
		}
	}
	else {
		InitSpriteData2dSoA(&sprite_data->real_time_sprites_buffer, &sprite_data->real_time_sprites, total_sprites);
		for (tfxEachLayer) {
			InitSpriteData2dSoA(&temp_sprites_buffer[layer], &temp_sprites[layer], 100);
		}
	}
	sprite_data->normal.total_memory_for_sprites = total_sprites * (tfxU32)sprite_data->real_time_sprites_buffer.struct_size;

	tfx_vector_t<tfxU32> running_count[tfxLAYERS];

	for (tfxEachLayer) {
		running_count[layer].resize(frames);
		running_count[layer].zero();
	}

	frame = 0;
	offset = 0;
	extra_frame_count = 0;
	particles_started = false;
	start_counting_extra_frames = false;
	DisablePMSpawning(pm, false);
	total_sprites = 0;
	tfxU32 captured_offset[tfxLAYERS] = { 0, 0, 0, 0 };
	pm->unique_particle_id = 0;

	while (frame < frames && offset < 99999) {
		tfxU32 count_this_frame = 0;
		UpdateParticleManager(pm, frame_length);
		for (tfxEachLayer) {
			if (is_3d) {
				InvalidateNewSpriteCapturedIndex(tfxCastBufferRef(tfx_3d_instance_t, pm->instance_buffer_for_recording[pm->current_sprite_buffer][layer]), pm->unique_sprite_ids[pm->current_sprite_buffer][layer], pm, layer);
			} else {
				InvalidateNewSpriteCapturedIndex(tfxCastBufferRef(tfx_2d_instance_t, pm->instance_buffer_for_recording[pm->current_sprite_buffer][layer]), pm->unique_sprite_ids[pm->current_sprite_buffer][layer], pm, layer);
			}
		}
		bool particles_processed_last_frame = false;

		if (offset >= start_frame) {
			if (offset == start_frame && start_frame > 0) {
				for (tfxEachLayer) {
					if (is_3d) {
						InvalidateOffsettedSpriteCapturedIndex(tfxCastBufferRef(tfx_3d_instance_t, pm->instance_buffer_for_recording[pm->current_sprite_buffer][layer]), pm->unique_sprite_ids[pm->current_sprite_buffer][layer], pm, layer);
					} else {
						InvalidateOffsettedSpriteCapturedIndex(tfxCastBufferRef(tfx_2d_instance_t, pm->instance_buffer_for_recording[pm->current_sprite_buffer][layer]), pm->unique_sprite_ids[pm->current_sprite_buffer][layer], pm, layer);
					}
				}
			}
			for (tfxEachLayer) {
				if (running_count[layer][frame] > 0 && pm->layer_sizes[layer] > 0) {
					//Copy instance_data that have looped round (for looped effects) into a temporary buffer, to be copied back after the fresh instance_data have been copied
					Resize(&temp_sprites_buffer[layer], running_count[layer][frame]);
					memcpy(temp_sprites[layer].uid, sprite_data->real_time_sprites.uid + frame_meta[frame].index_offset[layer], sizeof(tfx_unique_sprite_id_t) * running_count[layer][frame]);
					if (is_3d) {
						memcpy(temp_sprites[layer].billboard_instance, sprite_data->real_time_sprites.billboard_instance + frame_meta[frame].index_offset[layer], sizeof(tfx_3d_instance_t) * running_count[layer][frame]);
						if (captured_offset[layer] > 0) {
							for (int temp_i = 0; temp_i != temp_sprites_buffer[layer].current_size; ++temp_i) {
								if (temp_sprites[layer].billboard_instance[temp_i].captured_index != tfxINVALID) {
									temp_sprites[layer].billboard_instance[temp_i].captured_index += captured_offset[layer];
								}
							}
						}
					} else {
						memcpy(temp_sprites[layer].sprite_instance, sprite_data->real_time_sprites.sprite_instance + frame_meta[frame].index_offset[layer], sizeof(tfx_2d_instance_t) * running_count[layer][frame]);
						if (captured_offset[layer] > 0) {
							for (int temp_i = 0; temp_i != temp_sprites_buffer[layer].current_size; ++temp_i) {
								if (temp_sprites[layer].sprite_instance[temp_i].captured_index != tfxINVALID) {
									tfxU32 ci = temp_sprites[layer].sprite_instance[temp_i].captured_index & 0xFFFFFFF;
									ci += captured_offset[layer];
									temp_sprites[layer].sprite_instance[temp_i].captured_index += captured_offset[layer];
								}
							}
						}
					}
				} else if (captured_offset[layer] > 0 && pm->layer_sizes[layer] == 0) {
					if (is_3d) {
						for (int index = SpriteDataIndexOffset(sprite_data, frame, layer); index != SpriteDataEndIndex(sprite_data, frame, layer); ++index) {
							if (sprite_data->real_time_sprites.billboard_instance[index].captured_index != tfxINVALID) {
								sprite_data->real_time_sprites.billboard_instance[index].captured_index += captured_offset[layer];
							}
						}
					} else {
						for (int index = SpriteDataIndexOffset(sprite_data, frame, layer); index != SpriteDataEndIndex(sprite_data, frame, layer); ++index) {
							if (sprite_data->real_time_sprites.sprite_instance[index].captured_index != tfxINVALID) {
								tfxU32 ci = sprite_data->real_time_sprites.sprite_instance[index].captured_index & 0xFFFFFFF;
								ci += captured_offset[layer];
								sprite_data->real_time_sprites.sprite_instance[index].captured_index += captured_offset[layer];
							}
						}
					}
				}
			}
			for(tfxEachLayer) {
				tfx_2d_instance_t *sprites = tfxCastBufferRef(tfx_2d_instance_t, pm->instance_buffer_for_recording[pm->current_sprite_buffer][layer]);
				tfx_3d_instance_t *billboards = tfxCastBufferRef(tfx_3d_instance_t, pm->instance_buffer_for_recording[pm->current_sprite_buffer][layer]);
				//Copy fresh instance_data this frame
				memcpy(sprite_data->real_time_sprites.uid + frame_meta[frame].index_offset[layer], pm->unique_sprite_ids[pm->current_sprite_buffer][layer].data, sizeof(tfx_unique_sprite_id_t) * pm->layer_sizes[layer]);
				int index_offset = frame_meta[frame].index_offset[layer];
				int current_size = pm->layer_sizes[layer];
				if (is_3d) {
					memcpy(sprite_data->real_time_sprites.billboard_instance + frame_meta[frame].index_offset[layer], billboards, sizeof(tfx_3d_instance_t) * pm->layer_sizes[layer]);
				}
				else {
					memcpy(sprite_data->real_time_sprites.sprite_instance + frame_meta[frame].index_offset[layer], sprites, sizeof(tfx_2d_instance_t) * pm->layer_sizes[layer]);
				}

				if (running_count[layer][frame] > 0 && pm->layer_sizes[layer] > 0) {
					//Copy instance_data that have looped round (for looped effects)
					memcpy(sprite_data->real_time_sprites.uid + frame_meta[frame].index_offset[layer] + pm->layer_sizes[layer], temp_sprites[layer].uid, sizeof(tfx_unique_sprite_id_t) *temp_sprites_buffer[layer].current_size);
					if (is_3d) {
						memcpy(sprite_data->real_time_sprites.billboard_instance + frame_meta[frame].index_offset[layer] + pm->layer_sizes[layer], temp_sprites[layer].billboard_instance, sizeof(tfx_3d_instance_t) * temp_sprites_buffer[layer].current_size);
					}
					else {
						memcpy(sprite_data->real_time_sprites.sprite_instance + frame_meta[frame].index_offset[layer] + pm->layer_sizes[layer], temp_sprites[layer].sprite_instance, sizeof(tfx_2d_instance_t) * temp_sprites_buffer[layer].current_size);
					}
					tfxU32 layer_offset = 0;
					for (int l = 0; l <= layer; ++l) {
						layer_offset += pm->layer_sizes[l];
					}
					captured_offset[layer] = pm->layer_sizes[layer];
					if (layer + 1 < tfxLAYERS) {
						sprite_data->normal.frame_meta[frame].cumulative_offset[layer + 1] += pm->layer_sizes[layer];
					}
				}
				else if (pm->layer_sizes[layer] == 0) {
					captured_offset[layer] = 0;
				}
				running_count[layer][frame] += pm->layer_sizes[layer];
				total_sprites += pm->layer_sizes[layer];
				particles_started = total_sprites > 0;
				particles_processed_last_frame |= pm->layer_sizes[layer] > 0;
			}
		}

		offset++;

		if (particles_started)
			frame++;

		if (anim.animation_flags & tfxAnimationFlags_loop && particles_started) {
			if (frame == frames) {
				frame = 0;
				start_counting_extra_frames = true;
			}
			if (!particles_processed_last_frame)
				break;
		}
		if (start_counting_extra_frames && extra_frame_count++ >= extra_frames) {
			SoftExpireEffect(pm, preview_effect_index);
		}

	}

	sprite_data->real_time_sprites_buffer.current_size = total_sprites;
	//total instance_data should not exceed the capacity of the sprite buffer
	TFX_ASSERT(sprite_data->real_time_sprites_buffer.current_size <= sprite_data->real_time_sprites_buffer.capacity);
	ResetSpriteDataLerpOffset(sprite_data);

	for (int i = 0; i != anim.real_frames; ++i) {
		int previous_frame = i - 1;
		previous_frame = previous_frame < 0 ? frames - 1 : previous_frame;
		if (is_3d) {
			SpriteDataOffsetCapturedIndexes(sprite_data->real_time_sprites.billboard_instance, sprite_data, previous_frame, i);
		} else {
			SpriteDataOffsetCapturedIndexes(sprite_data->real_time_sprites.sprite_instance, sprite_data, previous_frame, i);
		}
	}

	if (is_3d) {
		WrapSingleParticleSprites(sprite_data->real_time_sprites.billboard_instance, sprite_data);
		ClearWrapBit(sprite_data->real_time_sprites.billboard_instance, sprite_data);
	}
	else {
		WrapSingleParticleSprites(sprite_data->real_time_sprites.sprite_instance, sprite_data);
		ClearWrapBit(sprite_data->real_time_sprites.sprite_instance, sprite_data);
	}

	for (tfxEachLayer) {
		FreeSoABuffer(&temp_sprites_buffer[layer]);
	}
	DisablePMSpawning(pm, false);
	ClearParticleManager(pm, false, false);
	pm->flags &= ~tfxParticleManagerFlags_recording_sprites;
	pm->flags |= pm_effect_sprites_flag;

	if (anim.playback_speed < 1.f) {
		CompressSpriteData(pm, effect, is_3d, frame_length, progress);
	}
	else {
		sprite_data->compressed_sprites = sprite_data->real_time_sprites;
		sprite_data->compressed_sprites_buffer = sprite_data->real_time_sprites_buffer;
		sprite_data->compressed = sprite_data->normal;
		anim.frames_after_compression = sprite_data->normal.frame_count;
	}

	pm->camera_position = pm_camera_position;
	ToggleSpritesWithUID(pm, false);
	for (tfxEachLayer) {
		pm->instance_buffer_for_recording[0][layer].free();
		pm->instance_buffer_for_recording[1][layer].free();
	}
}

void CompressSpriteData(tfx_particle_manager_t *pm, tfx_effect_emitter_t *effect, bool is_3d, float frame_length, std::atomic_int *progress) {
	progress->store(tfxLinkUpSprites);
	tfx_sprite_data_settings_t &anim = effect->library->sprite_data_settings[GetEffectInfo(effect)->sprite_data_settings_index];
	tfx_sprite_data_t *sprite_data = &effect->library->pre_recorded_effects.At(effect->path_hash);
	if (is_3d) {
		InitSpriteData3dSoACompression(&sprite_data->compressed_sprites_buffer, &sprite_data->compressed_sprites, tfxU32((float)sprite_data->real_time_sprites_buffer.current_size * sprite_data->frame_compression));
	}
	else {
		InitSpriteData2dSoACompression(&sprite_data->compressed_sprites_buffer, &sprite_data->compressed_sprites, tfxU32((float)sprite_data->real_time_sprites_buffer.current_size * sprite_data->frame_compression));
	}

	sprite_data->compressed.frame_meta.resize(tfxU32((float)anim.real_frames * anim.playback_speed) + 1);
	sprite_data->compressed.frame_meta.zero();

	float frequency = frame_length * (1.f / anim.playback_speed);
	float real_time = 0.f;
	float compressed_time = 0.f;
	int compressed_frame = 0;

	//First pass, add the instance_data from the real time sprite data to the compressed data
	int f = 0;
	tfx_sprite_data_soa_t &sprites = sprite_data->real_time_sprites;
	tfx_sprite_data_soa_t &c_sprites = sprite_data->compressed_sprites;
	tfxU32 ci = 0;
	bool frame_done = false;
	tfxU32 running_offset = 0;
	tfxU32 layer = 0;
	int start_frame = 0;
	bool finished = false;
	do {
		real_time = f * frame_length;
		tfxU32 next_compressed_frame = compressed_frame + 1;
		float next_compressed_time = next_compressed_frame * frequency;
		float lerp_offset = (next_compressed_time - real_time) / (next_compressed_time - compressed_time);
		if (real_time >= compressed_time && frame_done == false) {
			for (tfxU32 i = SpriteDataIndexOffset(sprite_data, f, layer); i != SpriteDataEndIndex(sprite_data, f, layer); ++i) {
				//Add to compress instance_data but make the invalid captured indexed create the offset
				sprite_data->compressed.frame_meta[compressed_frame].sprite_count[layer]++;
				sprite_data->compressed.frame_meta[compressed_frame].total_sprites++;
				c_sprites.uid[ci] = sprites.uid[i];
				if (is_3d) {
					c_sprites.lerp_offset[ci] = sprites.billboard_instance[i].captured_index == tfxINVALID ? lerp_offset : 1.f;
					c_sprites.billboard_instance[ci] = sprites.billboard_instance[i];
				}
				else {
					c_sprites.lerp_offset[ci] = sprites.sprite_instance[i].captured_index == tfxINVALID ? lerp_offset : 1.f;
					c_sprites.sprite_instance[ci] = sprites.sprite_instance[i];
				}
				ci++;
				if (ci >= sprite_data->compressed_sprites_buffer.capacity) {
					bool result = GrowArrays(&sprite_data->compressed_sprites_buffer, ci, sprite_data->compressed_sprites_buffer.capacity + 1, true);
					TFX_ASSERT(result);        //Failed to grow sprite compression array
				}
			}
			frame_done = true;
			f++;
		}
		else if (real_time > compressed_time && real_time < next_compressed_time && f < anim.real_frames) {
			for (tfxU32 i = SpriteDataIndexOffset(sprite_data, f, layer); i != SpriteDataEndIndex(sprite_data, f, layer); ++i) {
				bool captured = false;
				if (is_3d) {
					captured = sprites.billboard_instance[i].captured_index == tfxINVALID;
				}
				else {
					captured = sprites.sprite_instance[i].captured_index == tfxINVALID;
				}
				if (captured) {
					//Add to compressed instance_data frame but add the lerp offset
					sprite_data->compressed.frame_meta[compressed_frame].sprite_count[layer]++;
					sprite_data->compressed.frame_meta[compressed_frame].total_sprites++;
					c_sprites.uid[ci] = sprites.uid[i];
					c_sprites.lerp_offset[ci] = lerp_offset;
					if (is_3d) {
						c_sprites.billboard_instance[ci].captured_index = tfxINVALID;
						c_sprites.billboard_instance[ci] = sprites.billboard_instance[i];
					}
					else {
						c_sprites.sprite_instance[ci].captured_index = tfxINVALID;
						c_sprites.sprite_instance[ci] = sprites.sprite_instance[i];
					}
					ci++;
					if (ci >= sprite_data->compressed_sprites_buffer.capacity) {
						bool result = GrowArrays(&sprite_data->compressed_sprites_buffer, ci, sprite_data->compressed_sprites_buffer.capacity + 1, true);
						TFX_ASSERT(result);        //Failed to grow sprite compression array
					}
				}

			}
			f++;
		}
		else {
			sprite_data->compressed.frame_meta[compressed_frame].index_offset[layer] = running_offset;
			running_offset += sprite_data->compressed.frame_meta[compressed_frame].sprite_count[layer];
			layer++;
			if (layer == tfxLAYERS) {
				layer = 0;
				start_frame = f;
				if (start_frame >= anim.real_frames) {
					finished = true;
				}
				else {
					compressed_frame++;
					compressed_time = frequency * compressed_frame;
				}
			}
			else {
				f = start_frame;
			}
			frame_done = false;
		}
	} while (!finished);

	sprite_data->compressed.total_sprites = ci;
	sprite_data->compressed_sprites_buffer.current_size = ci;
	sprite_data->compressed.total_memory_for_sprites = ci * (tfxU32)sprite_data->compressed_sprites_buffer.struct_size;

	f = 0;
	//Second pass, link up the captured indexes using the UIDs
	sprite_data->compressed.frame_count = compressed_frame + 1;
	anim.animation_length_in_time = sprite_data->compressed.animation_length_in_time = sprite_data->compressed.frame_count * frequency;
	anim.frames_after_compression = sprite_data->compressed.frame_count;
	tfx_vector_t<tfx_compress_work_entry_t> compress_work;
	compress_work.reserve((int)sprite_data->compressed.frame_count);
	progress->store(tfxCompressFrames);
	while (f < (int)sprite_data->compressed.frame_count) {
		tfx_compress_work_entry_t *entry = &compress_work.next();
		entry->is_3d = is_3d;
		entry->sprite_data = sprite_data;
		entry->frame = f;
		if (tfxNumberOfThreadsInAdditionToMain > 0) {
			tfxAddWorkQueueEntry(&pm->work_queue, entry, LinkUpSpriteCapturedIndexes);
		}
		else {
			LinkUpSpriteCapturedIndexes(&pm->work_queue, entry);
		}
		f++;
	}
	tfxCompleteAllWork(&pm->work_queue);
	compress_work.free();
	progress->store(tfxBakingDone);

	TrimSoABuffer(&sprite_data->compressed_sprites_buffer);
	sprite_data->compressed.total_memory_for_sprites = (tfxU32)sprite_data->compressed_sprites_buffer.current_arena_size;

}

void LinkUpSpriteCapturedIndexes(tfx_work_queue_t *queue, void *work_entry) {
	tfx_compress_work_entry_t *entry = static_cast<tfx_compress_work_entry_t *>(work_entry);
	tfx_sprite_data_t *sprite_data = static_cast<tfx_sprite_data_t *>(entry->sprite_data);
	tfx_sprite_data_soa_t &c_sprites = sprite_data->compressed_sprites;

	if (entry->is_3d) {
		LinkSpriteDataCapturedIndexes(sprite_data->compressed_sprites.billboard_instance, entry->frame, sprite_data);
	}
	else {
		LinkSpriteDataCapturedIndexes(sprite_data->compressed_sprites.sprite_instance, entry->frame, sprite_data);
	}
}

void InitialiseAnimationManager(tfx_animation_manager_t *animation_manager, tfxU32 max_instances) {
	animation_manager->instances.reserve(max_instances);
	animation_manager->free_instances.reserve(max_instances);
	animation_manager->render_queue.reserve(max_instances);
	animation_manager->offsets.reserve(max_instances);
	animation_manager->instances_in_use[0].reserve(max_instances);
	animation_manager->instances_in_use[1].reserve(max_instances);
	animation_manager->current_in_use_buffer = 0;
	animation_manager->buffer_metrics.instances_size = 0;
	animation_manager->buffer_metrics.instances_size_in_bytes = 0;
	animation_manager->buffer_metrics.offsets_size = 0;
	animation_manager->buffer_metrics.offsets_size_in_bytes = 0;
	animation_manager->buffer_metrics.sprite_data_size = 0;
	animation_manager->buffer_metrics.total_sprites_to_draw = 0;
	animation_manager->update_frequency = 60.f;
	animation_manager->user_data = nullptr;
	animation_manager->maybe_render_instance_callback = nullptr;
    animation_manager->sprite_data_2d.set_alignment(16);
    animation_manager->sprite_data_3d.set_alignment(16);
	animation_manager->color_ramps.color_ramp_count = 0;
}

void InitialiseAnimationManagerFor3d(tfx_animation_manager_t *animation_manager, tfxU32 max_instances, tfxU32 initial_sprite_data_capacity) {
	InitialiseAnimationManager(animation_manager, max_instances);
	animation_manager->sprite_data_3d.reserve(initial_sprite_data_capacity);
	animation_manager->flags = tfxAnimationManagerFlags_initialised | tfxAnimationManagerFlags_is_3d;
}

void InitialiseAnimationManagerFor2d(tfx_animation_manager_t *animation_manager, tfxU32 max_instances, tfxU32 initial_sprite_data_capacity) {
	InitialiseAnimationManager(animation_manager, max_instances);
	animation_manager->sprite_data_2d.reserve(initial_sprite_data_capacity);
	animation_manager->flags = tfxAnimationManagerFlags_initialised;
}

tfxAnimationID AddAnimationInstance(tfx_animation_manager_t *animation_manager) {
	if (animation_manager->free_instances.current_size > 0) {
		tfxU32 index = animation_manager->free_instances.pop_back();
		animation_manager->instances_in_use[animation_manager->current_in_use_buffer].push_back(index);
		return index;
	}
	tfx_animation_instance_t instance;
	tfxU32 index = animation_manager->instances.current_size;
	TFX_ASSERT(animation_manager->instances.capacity != animation_manager->instances.current_size);        //At capacity! not enough room to add another instance.
	animation_manager->instances.push_back(instance);
	animation_manager->instances_in_use[animation_manager->current_in_use_buffer].push_back(index);
	return index;
}

void FreeAnimationInstance(tfx_animation_manager_t *animation_manager, tfxU32 index) {
	animation_manager->free_instances.push_back(index);
}

void AddEffectEmitterProperties(tfx_animation_manager_t *animation_manager, tfx_effect_emitter_t *effect, bool *has_animated_shape) {
	if (effect->type != tfxEmitterType) {
		for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
			AddEffectEmitterProperties(animation_manager, &sub, has_animated_shape);
		}
	}
	else {
		if (effect->library->emitter_properties[effect->property_index].animation_property_index != tfxINVALID) {
			tfx_animation_emitter_properties_t properties;
			properties.handle = effect->library->emitter_properties[effect->property_index].image_handle;
			properties.handle_packed = effect->library->emitter_properties[effect->property_index].image_handle_packed;
			properties.flags = effect->property_flags;
			tfx_image_data_t &image = *effect->library->emitter_properties[effect->property_index].image;
			properties.animation_frames = image.animation_frames;
			if (properties.animation_frames > 1 && effect->property_flags & tfxEmitterPropertyFlags_animate) {
				*has_animated_shape = true;
			}
			if (animation_manager->particle_shapes.ValidKey(image.image_hash)) {
				properties.start_frame_index = animation_manager->particle_shapes.At(image.image_hash).compute_shape_index;
			}
			else {
				properties.start_frame_index = image.compute_shape_index;
			}
			effect->library->emitter_properties[effect->property_index].animation_property_index = animation_manager->emitter_properties.current_size;
			tfx__readbarrier;
			tfxU32 index = animation_manager->emitter_properties.current_size;
			tfx_overtime_attributes_t &a = effect->library->emitter_attributes[effect->emitter_attributes].overtime;
			animation_manager->emitter_properties.push_back_copy(properties);
			CopyColorRampToAnimationManager(animation_manager, index, &a.color_ramps[0]);
			for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
				AddEffectEmitterProperties(animation_manager, &sub, has_animated_shape);
			}
		}
	}
}

void AddEffectShapes(tfx_animation_manager_t *animation_manager, tfx_effect_emitter_t *effect) {
	if (effect->type == tfxEmitterType) {
		tfx_image_data_t *image_data = GetEffectProperties(effect)->image;
		if (!animation_manager->particle_shapes.ValidKey(image_data->image_hash)) {
			animation_manager->particle_shapes.Insert(image_data->image_hash, *image_data);
		}
		for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
			AddEffectShapes(animation_manager, &sub);
		}
	}
	else {
		for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
			AddEffectShapes(animation_manager, &sub);
		}
	}
}

void AddSpriteData(tfx_animation_manager_t *animation_manager, tfx_effect_emitter_t *effect, tfx_particle_manager_t *pm, tfx_vec3_t camera_position) {
	if (Is3DEffect(effect)) {
		//If you're adding 3d effect sprite data then the animation manager must have been initialised with InitialiseAnimationManagerFor3d
		TFX_ASSERT(animation_manager->flags & tfxAnimationManagerFlags_is_3d);
	}
	else {
		//If you're adding 2d effect sprite data then the animation manager must have been initialised with InitialiseAnimationManagerFor2d
		TFX_ASSERT(!(animation_manager->flags & tfxAnimationManagerFlags_is_3d));
	}
	tfx_sprite_data_settings_t &anim = effect->library->sprite_data_settings[GetEffectInfo(effect)->sprite_data_settings_index];
	if (!effect->library->pre_recorded_effects.ValidKey(effect->path_hash)) {
		TFX_ASSERT(pm);        //You must pass an appropriate particle manager if the animation needs recording
		std::atomic_int progress;
		RecordSpriteData(pm, effect, animation_manager->update_frequency, &camera_position.x, &progress);
	}

	bool has_animated_shape = false;
	AddEffectEmitterProperties(animation_manager, effect, &has_animated_shape);

	bool is_3d = Is3DEffect(effect);
	tfx_sprite_data_t &sprite_data = effect->library->pre_recorded_effects.At(effect->path_hash);
	animation_manager->effect_animation_info.Insert(effect->path_hash, sprite_data.compressed);
	tfx_sprite_data_metrics_t &metrics = animation_manager->effect_animation_info.At(effect->path_hash);
	metrics.name = GetEffectInfo(effect)->name;
	metrics.frames_after_compression = anim.frames_after_compression;
	metrics.real_frames = anim.real_frames;
	metrics.animation_length_in_time = anim.animation_length_in_time;
	metrics.animation_flags = anim.animation_flags;
	metrics.flags = has_animated_shape ? tfxAnimationManagerFlags_has_animated_shapes : 0;
	tfx_sprite_data_soa_t &sprites = sprite_data.compressed_sprites;
	if (is_3d) {
		metrics.start_offset = animation_manager->sprite_data_3d.current_size;
		for (int i = 0; i != metrics.total_sprites; ++i) {
			tfx_sprite_data3d_t sprite;
			memcpy(&sprite, &sprites.billboard_instance[i], sizeof(tfx_3d_instance_t));
			sprite.captured_index += sprite.captured_index == tfxINVALID ? 0 : metrics.start_offset;
			sprite.additional = tfxU32(sprites.lerp_offset[i] * 65535.f);
			tfxU32 animation_property_index = effect->library->emitter_properties[sprites.uid[i].property_index].animation_property_index;
			sprite.additional |= (animation_property_index << 16);
			sprite.indexes &= 0x0000FFFF;
			sprite.indexes |= (tfxColorRampIndex(animation_manager->emitter_properties[animation_property_index].color_ramp_index) << 24);
			sprite.indexes |= (tfxColorRampLayer(animation_manager->emitter_properties[animation_property_index].color_ramp_index) << 16);
			animation_manager->sprite_data_3d.push_back_copy(sprite);
		}
		metrics.total_memory_for_sprites = sizeof(tfx_sprite_data3d_t) * metrics.total_sprites;
	}
	else {
		metrics.start_offset = animation_manager->sprite_data_2d.current_size;
		for (int i = 0; i != metrics.total_sprites; ++i) {
			tfx_sprite_data2d_t sprite;
			memcpy(&sprite, &sprites.sprite_instance[i], sizeof(tfx_2d_instance_t));
			sprite.captured_index += sprite.captured_index == tfxINVALID ? 0 : metrics.start_offset;
			sprite.additional = tfxU32(sprites.lerp_offset[i] * 65535.f);
			tfxU32 animation_property_index = effect->library->emitter_properties[sprites.uid[i].property_index].animation_property_index;
			sprite.additional |= (animation_property_index << 16);
			sprite.indexes &= 0x0000FFFF;
			sprite.indexes |= (tfxColorRampIndex(animation_manager->emitter_properties[animation_property_index].color_ramp_index) << 24);
			sprite.indexes |= (tfxColorRampLayer(animation_manager->emitter_properties[animation_property_index].color_ramp_index) << 16);
			animation_manager->sprite_data_2d.push_back_copy(sprite);
		}
		metrics.total_memory_for_sprites = sizeof(tfx_sprite_data2d_t) * metrics.total_sprites;
	}

	//Update the index offset frame meta in the metrics depending on where this sprite data is being inserted into
	//the list (the list may contain many different effect animations)
	for (tfxEachLayer) {
		for (auto &meta : metrics.frame_meta) {
			meta.index_offset[layer] += metrics.start_offset;
		}
	}

	animation_manager->buffer_metrics.sprite_data_size += metrics.total_memory_for_sprites;
}

void SetAnimationManagerInstanceCallback(tfx_animation_manager_t *animation_manager, bool((*maybe_render_instance_callback)(tfx_animation_manager_t *animation_manager, tfx_animation_instance_t *instance, tfx_frame_meta_t *meta, void *user_data))) {
	animation_manager->maybe_render_instance_callback = maybe_render_instance_callback;
}

void SetAnimationManagerUserData(tfx_animation_manager_t *animation_manager, void *user_data) {
	animation_manager->user_data = user_data;
}

tfxAnimationID AddAnimationInstance(tfx_animation_manager_t *animation_manager, tfxKey path, tfxU32 start_frame) {
	TFX_ASSERT(animation_manager->effect_animation_info.ValidKey(path));                //You must have added the effect sprite data to the animation manager
	//Call AddSpriteData to do so
	if (animation_manager->instances_in_use->current_size >= animation_manager->instances_in_use->capacity) {
		return tfxINVALID;
	}
	tfxU32 info_index = animation_manager->effect_animation_info.GetIndex(path);
	tfx_sprite_data_metrics_t &metrics = animation_manager->effect_animation_info.data[info_index];
	TFX_ASSERT(start_frame < metrics.frames_after_compression);
	tfxAnimationID index = AddAnimationInstance(animation_manager);
	tfx_animation_instance_t &instance = animation_manager->instances[index];
	instance.scale = 1.f;
	float frame_length = 1000.f / animation_manager->update_frequency;
	float frame_length_total = float(metrics.real_frames) / float(metrics.frames_after_compression) * frame_length;
	instance.current_time = start_frame * frame_length_total;
	instance.animation_length_in_time = metrics.animation_length_in_time;
	instance.tween = 0.f;
	instance.flags = metrics.animation_flags;
	instance.info_index = info_index;
	instance.offset_into_sprite_data = metrics.start_offset;
	instance.sprite_count = metrics.frame_meta[start_frame].total_sprites;
	instance.frame_count = metrics.frame_count;
	return index;
}

tfxAnimationID AddAnimationInstance(tfx_animation_manager_t *animation_manager, const char *path, tfxU32 start_frame) {
	tfxKey path_hash = tfxXXHash64::hash(path, strlen(path), 0);
	return AddAnimationInstance(animation_manager, path_hash, start_frame);
}

void UpdateAnimationManager(tfx_animation_manager_t *animation_manager, float elapsed) {
	TFX_ASSERT(animation_manager->instances_in_use[animation_manager->current_in_use_buffer].capacity > 0);    //You must call InitialiseAnimationManager before trying to update one
	tfxU32 next_buffer = animation_manager->current_in_use_buffer ^ 1;
	animation_manager->instances_in_use[next_buffer].clear();
	animation_manager->render_queue.clear();
	animation_manager->offsets.clear();
	tfxU32 running_sprite_count = 0;
	animation_manager->flags &= ~tfxAnimationManagerFlags_has_animated_shapes;
	for (auto i : animation_manager->instances_in_use[animation_manager->current_in_use_buffer]) {
		tfx_animation_instance_t &instance = animation_manager->instances[i];
		tfx_sprite_data_metrics_t &metrics = animation_manager->effect_animation_info.data[instance.info_index];
		instance.current_time += elapsed;
		float frame_time = (instance.current_time / instance.animation_length_in_time) * (float)instance.frame_count;
		tfxU32 frame = tfxU32(frame_time);
		frame++;
		frame = frame >= metrics.frame_count ? 0 : frame;
		if (instance.current_time >= instance.animation_length_in_time) {
			if (instance.flags & tfxAnimationInstanceFlags_loop) {
				instance.sprite_count = metrics.frame_meta[0].total_sprites;
				instance.offset_into_sprite_data = metrics.frame_meta[0].index_offset[0];
				instance.current_time -= instance.animation_length_in_time;
				animation_manager->instances_in_use[next_buffer].push_back(i);
				if (animation_manager->maybe_render_instance_callback && !animation_manager->maybe_render_instance_callback(animation_manager, &instance, &metrics.frame_meta[0], animation_manager->user_data)) {
					continue;
				}
				running_sprite_count += instance.sprite_count;
				animation_manager->offsets.push_back(running_sprite_count);
				animation_manager->flags |= metrics.flags & tfxAnimationManagerFlags_has_animated_shapes;
				animation_manager->render_queue.push_back(instance);
				UpdateAnimationManagerBufferMetrics(animation_manager);
			}
			else {
				FreeAnimationInstance(animation_manager, i);
			}
		}
		else {
			instance.sprite_count = metrics.frame_meta[frame].total_sprites;
			instance.offset_into_sprite_data = metrics.frame_meta[frame].index_offset[0];
			animation_manager->instances_in_use[next_buffer].push_back(i);
			if (animation_manager->maybe_render_instance_callback && !animation_manager->maybe_render_instance_callback(animation_manager, &instance, &metrics.frame_meta[frame], animation_manager->user_data)) {
				continue;
			}
			animation_manager->render_queue.push_back(instance);
			running_sprite_count += instance.sprite_count;
			animation_manager->flags |= metrics.flags & tfxAnimationManagerFlags_has_animated_shapes;
			animation_manager->offsets.push_back(running_sprite_count);
			UpdateAnimationManagerBufferMetrics(animation_manager);
		}
	}
	animation_manager->buffer_metrics.total_sprites_to_draw = running_sprite_count;
	animation_manager->current_in_use_buffer = animation_manager->current_in_use_buffer ^ 1;
}

void CycleAnimationManager(tfx_animation_manager_t *animation_manager) {
	TFX_ASSERT(animation_manager->instances_in_use[animation_manager->current_in_use_buffer].capacity > 0);    //You must call InitialiseAnimationManager before trying to update one
	tfxU32 next_buffer = animation_manager->current_in_use_buffer ^ 1;
	animation_manager->instances_in_use[next_buffer].clear();
	animation_manager->render_queue.clear();
	animation_manager->offsets.clear();
	tfxU32 running_sprite_count = 0;
	animation_manager->flags &= ~tfxAnimationManagerFlags_has_animated_shapes;
	for (auto i : animation_manager->instances_in_use[animation_manager->current_in_use_buffer]) {
		tfx_animation_instance_t &instance = animation_manager->instances[i];
		tfx_sprite_data_metrics_t &metrics = animation_manager->effect_animation_info.data[instance.info_index];
		float frame_time = (instance.current_time / instance.animation_length_in_time) * (float)instance.frame_count;
		tfxU32 frame = tfxU32(frame_time) + 1;
		frame = frame >= metrics.frame_count ? 0 : frame;
		instance.sprite_count = metrics.frame_meta[frame].total_sprites;
		instance.offset_into_sprite_data = metrics.frame_meta[frame].index_offset[0];
		animation_manager->instances_in_use[next_buffer].push_back(i);
		animation_manager->render_queue.push_back(instance);
		running_sprite_count += instance.sprite_count;
		animation_manager->flags |= metrics.flags & tfxAnimationManagerFlags_has_animated_shapes;
		animation_manager->offsets.push_back(running_sprite_count);
		UpdateAnimationManagerBufferMetrics(animation_manager);
	}
	animation_manager->buffer_metrics.total_sprites_to_draw = running_sprite_count;
	animation_manager->current_in_use_buffer = animation_manager->current_in_use_buffer ^ 1;
}

void ClearAllAnimationInstances(tfx_animation_manager_t *animation_manager) {
	animation_manager->free_instances.clear();
	animation_manager->instances_in_use[0].clear();
	animation_manager->instances_in_use[1].clear();
	animation_manager->render_queue.clear();
	animation_manager->instances.clear();
	animation_manager->offsets.clear();
}

void ResetAnimationManager(tfx_animation_manager_t *animation_manager) {
	animation_manager->free_instances.clear();
	animation_manager->instances_in_use[0].clear();
	animation_manager->instances_in_use[1].clear();
	animation_manager->render_queue.clear();
	animation_manager->instances.clear();
	animation_manager->offsets.clear();
	animation_manager->sprite_data_2d.clear();
	animation_manager->sprite_data_3d.clear();
	animation_manager->emitter_properties.clear();
	animation_manager->effect_animation_info.Clear();
	animation_manager->particle_shapes.Clear();
	animation_manager->color_ramps.color_ramp_bitmaps.clear();
	animation_manager->color_ramps.color_ramp_ids.Clear();
	animation_manager->color_ramps.color_ramp_count = 0;

	animation_manager->buffer_metrics.instances_size = 0;
	animation_manager->buffer_metrics.instances_size_in_bytes = 0;
	animation_manager->buffer_metrics.offsets_size = 0;
	animation_manager->buffer_metrics.offsets_size_in_bytes = 0;
	animation_manager->buffer_metrics.sprite_data_size = 0;
	animation_manager->buffer_metrics.total_sprites_to_draw = 0;
}

void FreeAnimationManager(tfx_animation_manager_t *animation_manager) {
	animation_manager->free_instances.free_all();
	animation_manager->instances_in_use[0].free_all();
	animation_manager->instances_in_use[1].free_all();
	animation_manager->render_queue.free_all();
	animation_manager->instances.free_all();
	animation_manager->offsets.free_all();
	animation_manager->sprite_data_2d.free_all();
	animation_manager->sprite_data_3d.free_all();
	animation_manager->emitter_properties.free_all();
	animation_manager->effect_animation_info.FreeAll();
	animation_manager->flags = 0;
}

void UpdateAnimationManagerBufferMetrics(tfx_animation_manager_t *animation_manager) {
	animation_manager->buffer_metrics.instances_size = animation_manager->render_queue.current_size;
	animation_manager->buffer_metrics.offsets_size = animation_manager->offsets.current_size;
	animation_manager->buffer_metrics.instances_size_in_bytes = animation_manager->buffer_metrics.instances_size * sizeof(tfx_animation_instance_t) * animation_manager->buffer_metrics.instances_size;
	animation_manager->buffer_metrics.offsets_size_in_bytes = animation_manager->buffer_metrics.offsets_size * sizeof(tfxU32) * animation_manager->buffer_metrics.instances_size;
}

void RecordTemplateEffect(tfx_effect_template_t *t, tfx_particle_manager_t *pm, float update_frequency, float camera_position[3]) {
	std::atomic_int progress;
	RecordSpriteData(pm, &t->effect, update_frequency, camera_position, &progress);
}

void DisableTemplateEmitter(tfx_effect_template_t *t, const char *path) {
	TFX_ASSERT(t->paths.ValidName(path));            //Must be a valid path to the emitter
	tfx_effect_emitter_t *emitter = t->paths.At(path);
	TFX_ASSERT(emitter->type == tfxEmitterType);    //Must be an emitter that you're trying to remove. Use RemoveSubEffect if you're trying to remove one of those. 
	emitter->property_flags &= ~tfxEmitterPropertyFlags_enabled;
}

void EnableTemplateEmitter(tfx_effect_template_t *t, const char *path) {
	TFX_ASSERT(t->paths.ValidName(path));            //Must be a valid path to the emitter
	tfx_effect_emitter_t *emitter = t->paths.At(path);
	TFX_ASSERT(emitter->type == tfxEmitterType);    //Must be an emitter that you're trying to remove. Use RemoveSubEffect if you're trying to remove one of those
	emitter->property_flags |= tfxEmitterPropertyFlags_enabled;
}

void ScaleTemplateGlobalMultiplier(tfx_effect_template_t *t, tfx_graph_type global_type, float amount) {
	TFX_ASSERT(IsGlobalGraph(global_type));
	tfx_graph_t *graph = GetEffectGraphByType(&t->effect, global_type);
	tfx_effect_emitter_t *original_effect = GetLibraryEffect(t->effect.library, t->original_effect_hash);
	tfx_graph_t *original_graph = GetEffectGraphByType(original_effect, global_type);
	CopyGraph(original_graph, graph, false);
	MultiplyAllGraphValues(graph, amount);
	CompileGraph(graph);
}

void ScaleTemplateEmitterGraph(tfx_effect_template_t *t, const char *emitter_path, tfx_graph_type graph_type, float amount) {
	TFX_ASSERT(IsEmitterGraph(graph_type));        //Must be an emitter graph type. This is any property, base, variaion or overtime graph
	TFX_ASSERT(t->paths.ValidName(emitter_path));            //Must be a valid path to the emitter
	tfx_effect_emitter_t *emitter = t->paths.At(emitter_path);
	tfx_graph_t *graph = GetEffectGraphByType(emitter, graph_type);
	tfx_effect_emitter_t *original_emitter = GetLibraryEffect(t->effect.library, emitter_path);
	tfx_graph_t *original_graph = GetEffectGraphByType(original_emitter, graph_type);
	CopyGraph(original_graph, graph, false);
	MultiplyAllGraphValues(graph, amount);
	CompileGraph(graph);
}

void SetTemplateSingleSpawnAmount(tfx_effect_template_t *t, const char *emitter_path, tfxU32 amount) {
	TFX_ASSERT(amount >= 0);                            //Amount must not be less than 0
	TFX_ASSERT(t->paths.ValidName(emitter_path));            //Must be a valid path to the emitter
	tfx_effect_emitter_t *emitter = t->paths.At(emitter_path);
	GetEffectProperties(emitter)->spawn_amount = amount;
}

void *GetAnimationEmitterPropertiesBufferPointer(tfx_animation_manager_t *animation_manager) {
	return animation_manager->emitter_properties.data;
}

void ResetTemplate(tfx_effect_template_t *t) {
	if (t->paths.Size()) {
		t->paths.Clear();
		CleanUpEffect(&t->effect);
	}
}

tfx_effect_emitter_t *GetEffectFromTemplate(tfx_effect_template_t *t) {
	return &t->effect;
}

tfx_effect_emitter_t *GetEmitterFromTemplate(tfx_effect_template_t *t, tfx_str256_t *path) {
	if (t->paths.ValidName(*path)) return t->paths.At(*path); return nullptr;
}

tfx_emitter_path_t *GetEmitterPath(tfx_effect_emitter_t *e) {
	if (e->path_attributes != tfxINVALID) {
		TFX_ASSERT(e->library->paths.size() > e->path_attributes); //The emitter path attributes is out of bounds. This really shouldn't happen, either a bug in the library or the path attributes was set manually and incorrectly.
		return &e->library->paths[e->path_attributes];
	}
	return nullptr;
}

void SetTemplateUserData(tfx_effect_template_t *t, tfx_str256_t *path, void *data) {
	if (t->paths.ValidName(*path)) t->paths.At(*path)->user_data = data;
}

void SetTemplateEffectUserData(tfx_effect_template_t *t, void *data) {
	t->effect.user_data = data;
}

void SetTemplateUserDataAll(tfx_effect_template_t *t, void *data);

void SetTemplateEffectUpdateCallback(tfx_effect_template_t *t, void(*update_callback)(tfx_particle_manager_t *pm, tfxEffectID effect_index)) {
	t->effect.update_callback = update_callback;
}

tfx_particle_manager_t::~tfx_particle_manager_t() {
}

bool AddEffectToParticleManager(tfx_particle_manager_t *pm, tfx_effect_template_t *effect_template, tfxEffectID *effect_id) {
	tfxEffectID id;
	id = AddEffectToParticleManager(pm, &effect_template->effect, pm->current_ebuff, 0, false, 0, 0.f);
	if (effect_id) {
		*effect_id = id;
	}
	return id != tfxINVALID;
}

bool AddEffectToParticleManager(tfx_particle_manager_t *pm, tfx_effect_emitter_t *effect, tfxEffectID *effect_id) {
	tfxEffectID id;
	id = AddEffectToParticleManager(pm, effect, pm->current_ebuff, 0, false, 0, 0.f);
	if (effect_id) {
		*effect_id = id;
	}
	return id != tfxINVALID;
}

tfxEffectID AddEffectToParticleManager(tfx_particle_manager_t *pm, tfx_effect_emitter_t *effect, int buffer, int hierarchy_depth, bool is_sub_emitter, tfxU32 root_effect_index, float add_delayed_spawning) {
	tfxPROFILE;
	std::lock_guard<std::mutex> lock(pm->add_effect_mutex);
	TFX_ASSERT(effect->type == tfxEffectType);
	TFX_ASSERT(effect->library == pm->library);    //The effect must belong to the same library that is assigned to the particle manager
	if (pm->flags & tfxParticleManagerFlags_use_compute_shader && pm->highest_compute_controller_index >= pm->max_compute_controllers && pm->free_compute_controllers.empty()) {
		return tfxINVALID;
	}
	tfx_effect_index_t parent_index = GetPMEffectSlot(pm);
	if (parent_index.index == tfxINVALID) {
		return tfxINVALID;
	}
	if (!is_sub_emitter) {
		pm->effects[parent_index.index].highest_particle_age = pm->frame_length * 3.f;
	}
	tfx_emitter_properties_t *properties = GetEffectProperties(effect);
	tfx_effect_state_t &new_effect = pm->effects[parent_index.index];
	new_effect.path_hash = effect->path_hash;
	new_effect.global_attributes = effect->global;
	new_effect.transform_attributes = effect->transform_attributes;
	new_effect.age = -add_delayed_spawning;
	new_effect.state_flags = 0;
	new_effect.frame = 0.f;
	new_effect.property_flags = effect->property_flags;
	new_effect.effect_flags = effect->effect_flags;
	new_effect.local_position = tfx_vec3_t();
	new_effect.timeout = 1000.f;
	new_effect.library = effect->library;
	new_effect.parent_particle_index = tfxINVALID;
	new_effect.info_index = effect->info_index;
	new_effect.properties_index = effect->property_index;
	new_effect.timeout_counter = 0;
	new_effect.user_data = effect->user_data;
	new_effect.update_callback = effect->update_callback;
	float range = properties->noise_base_offset_range;
	new_effect.noise_base_offset = RandomRange(&pm->random, range);
	pm->effects_in_use[hierarchy_depth][buffer].push_back(parent_index);
	pm->sort_passes = tfxMax(effect->sort_passes, pm->sort_passes);
	pm->sort_passes = tfxMin(5, pm->sort_passes);
	new_effect.sort_passes = pm->sort_passes;
	new_effect.instance_data.instance_start_index = tfxINVALID;
	new_effect.emitter_indexes[0].clear();
	new_effect.emitter_indexes[1].clear();
	new_effect.emitter_start_size = 0;

	tfxU32 seed_index = 0;
	struct hash_index_pair_t {
		tfxKey hash;
		tfxU32 index;
	};
	tmpStack(hash_index_pair_t, source_emitters);
	tmpStack(hash_index_pair_t, target_emitters);
	for (auto &e : GetEffectInfo(effect)->sub_effectors) {
		if (e.property_flags & tfxEmitterPropertyFlags_enabled) {
			unsigned int index = GetPMEmitterSlot(pm);
			if (index == tfxINVALID) {
				break;
			}
			new_effect.emitter_indexes[pm->current_ebuff].locked_push_back(index);
			tfx_emitter_state_t &emitter = pm->emitters[index];
			emitter.particles_index = tfxINVALID;
			pm->emitters[index].parent_index = parent_index.index;
			tfx_emitter_properties_t *emitter_properties = GetEffectProperties(&e);
			emitter.path_quaternions = nullptr;
			emitter.path_hash = e.path_hash;
			emitter.info_index = e.info_index;
			emitter.properties_index = e.property_index;
			emitter.emitter_attributes = e.emitter_attributes;
			emitter.transform_attributes = e.transform_attributes;
			emitter.path_attributes = e.path_attributes;
			emitter.delay_spawning = emitter_properties->delay_spawning;
			emitter.age = 0.f;
			emitter.frame = 0.f;
			emitter.local_position = tfx_vec3_t();
			emitter.grid_coords = tfx_vec3_t();
			emitter.grid_direction = tfx_vec3_t();
			emitter.property_flags = e.property_flags;
			emitter.image_size = emitter_properties->image->image_size;
			emitter.image_frame_rate = emitter_properties->image->animation_frames > 1 && e.property_flags & tfxEmitterPropertyFlags_animate ? emitter_properties->frame_rate : 0.f;
			//emitter.image_frame_rate = e.property_flags & tfxEmitterPropertyFlags_reverse_animation ? emitter.image_frame_rate * -1.f : emitter.image_frame_rate;
			emitter.end_frame = emitter_properties->end_frame;
			emitter.angle_offsets = emitter_properties->angle_offsets;
			emitter.timeout = 1000.f;
			emitter.amount_remainder = 0.f;
			emitter.qty_step_size = 0.f;
			emitter.timeout_counter = 0;
			emitter.emitter_size = 0.f;
			emitter.hierarchy_depth = hierarchy_depth;
			emitter.world_rotations = 0.f;
			emitter.seed_index = seed_index++;
			emitter.control_profile = e.control_profile;
			emitter.spawn_locations_index = tfxINVALID;
			//----Handle
			if (e.property_flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
				emitter.image_handle_packed = (tfxU64)tfx__pack16bit_sscaled(0.5f, 0.5f, 128.f) << 32;
			}
			else {
				emitter.image_handle_packed = (tfxU64)emitter_properties->image_handle_packed << 32;
			}
			tfxEmitterStateFlags &state_flags = emitter.state_flags;
			const tfxEmitterStateFlags &parent_state_flags = new_effect.state_flags;

			state_flags = tfxEmitterStateFlags_no_tween_this_update;
			state_flags |= parent_state_flags & tfxEffectStateFlags_no_tween;
			state_flags |= (e.property_flags & tfxEmitterPropertyFlags_wrap_single_sprite) && emitter_properties->single_shot_limit == 0 ? tfxEmitterStateFlags_wrap_single_sprite : 0;
			state_flags |= e.property_flags & tfxEmitterPropertyFlags_play_once;
			state_flags |= e.property_flags & tfxEmitterPropertyFlags_single && !(pm->flags & tfxParticleManagerFlags_disable_spawning) ? tfxEmitterStateFlags_is_single : 0;
			state_flags |= e.property_flags & tfxEmitterPropertyFlags_base_uniform_size;
			state_flags |= (emitter_properties->emission_type != tfxLine && !(e.property_flags & tfxEmitterPropertyFlags_edge_traversal)) || (emitter_properties->emission_type == tfxLine && !(e.property_flags & tfxEmitterPropertyFlags_edge_traversal)) ? tfxEmitterStateFlags_not_line : 0;
			state_flags |= e.property_flags & tfxEmitterPropertyFlags_random_color;
			state_flags |= e.property_flags & tfxEmitterPropertyFlags_lifetime_uniform_size;
			state_flags |= emitter_properties->angle_settings != tfxAngleSettingFlags_align_roll && !(e.property_flags & tfxEmitterPropertyFlags_relative_angle) ? tfxEmitterStateFlags_can_spin : 0;
			state_flags |= (emitter_properties->angle_settings & tfxAngleSettingFlags_align_roll) ? tfxEmitterStateFlags_align_with_velocity : 0;
			state_flags |= emitter_properties->emission_type == tfxLine && e.property_flags & tfxEmitterPropertyFlags_edge_traversal ? tfxEmitterStateFlags_is_edge_traversal : 0;
			state_flags |= emitter_properties->emission_type == tfxPath && e.property_flags & tfxEmitterPropertyFlags_edge_traversal ? tfxEmitterStateFlags_is_edge_traversal : 0;
			state_flags |= emitter_properties->end_behaviour == tfxLoop ? tfxEmitterStateFlags_loop : 0;
			state_flags |= emitter_properties->end_behaviour == tfxKill ? tfxEmitterStateFlags_kill : 0;
			state_flags |= emitter_properties->emission_type == tfxLine && e.property_flags & tfxEmitterPropertyFlags_edge_traversal && (state_flags & tfxEmitterStateFlags_loop || state_flags & tfxEmitterStateFlags_kill) ? tfxEmitterStateFlags_is_line_loop_or_kill : 0;
			state_flags |= (effect->property_flags & tfxEmitterPropertyFlags_effect_is_3d) && (emitter_properties->billboard_option == tfxBillboarding_free_align || emitter_properties->billboard_option == tfxBillboarding_align_to_vector) ? tfxEmitterStateFlags_can_spin_pitch_and_yaw : 0;
			state_flags |= emitter_properties->emission_type == tfxPath ? tfxEmitterStateFlags_has_path : 0;
			if (emitter_properties->emission_type == tfxPath) {
				tfx_emitter_path_t *path = &pm->library->paths[emitter.path_attributes];
				state_flags |= (path->rotation_range > 0) ? tfxEmitterStateFlags_has_rotated_path : 0;
				emitter.last_path_index = 0;
				emitter.active_paths = (emitter.state_flags & tfxEmitterStateFlags_has_rotated_path) && path->rotation_stagger == 0 ? path->maximum_active_paths : 1;
				emitter.path_stagger_counter = 0.f;
				emitter.path_quaternion_index = AllocatePathQuaterion(pm, path->maximum_active_paths);
				emitter.path_quaternions = pm->path_quaternions[emitter.path_quaternion_index];
				emitter.path_quaternions[0].grid_coord = (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) ? 0.f : (float)path->node_count - 4;
				emitter.path_quaternions[0].cycles = 0;
				emitter.path_cycle_count = path->maximum_paths;
				emitter.path_start_index = 0;
				if (emitter.state_flags & tfxEmitterStateFlags_has_rotated_path) {
					for (int qi = 0; qi != path->maximum_active_paths; ++qi) {
						emitter.path_quaternions[qi].cycles = tfxINVALID;
					}
					for (int qi = 0; qi != emitter.active_paths; ++qi) {
						if (path->flags & tfxPathFlags_2d) {
							tfx_quaternion_t q = GetPathRotation2d(&pm->random, path->rotation_range, path->rotation_pitch);
							emitter.path_quaternions[qi].quaternion = tfx__pack8bit_quaternion(q);
						}
						else {
							tfx_quaternion_t q = GetPathRotation(&pm->random, path->rotation_range, path->rotation_pitch, path->rotation_yaw, ((path->flags & tfxPathFlags_rotation_range_yaw_only) > 0));
							emitter.path_quaternions[qi].quaternion = tfx__pack8bit_quaternion(q);
						}
						emitter.path_quaternions[qi].grid_coord = (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) ? 0.f : (float)path->node_count - 4;
						emitter.path_quaternions[qi].age = 0.f;
						emitter.path_quaternions[qi].cycles = 0;
						if (emitter.path_cycle_count > 0) {
							emitter.path_cycle_count--;
						}
						else if (path->maximum_paths > 0) {
							emitter.path_quaternions[qi].cycles = tfxINVALID;
						}
					}
				}
			}

			if (emitter.particles_index == tfxINVALID) {
				if (!is_sub_emitter) {
					emitter.particles_index = tfx__grab_particle_lists(pm, e.path_hash, (effect->property_flags & tfxEmitterPropertyFlags_effect_is_3d), 100, e.control_profile);
				}
			}

			emitter.property_flags |= (effect->property_flags & tfxEmitterPropertyFlags_effect_is_3d);
			if (state_flags & tfxEmitterStateFlags_is_edge_traversal) {
				emitter.property_flags |= tfxEmitterPropertyFlags_relative_position;
			}

			if (is_sub_emitter) {
				state_flags |= tfxEmitterStateFlags_is_sub_emitter;
				emitter.root_index = root_effect_index;
				state_flags |= IsOrderedEffectState(&pm->effects[root_effect_index]) ? tfxEmitterStateFlags_is_in_ordered_effect : 0;
			}
			else {
				emitter.root_index = parent_index.index;
				emitter.highest_particle_age = pm->frame_length * 2.f;
				state_flags |= IsOrderedEffect(effect) ? tfxEmitterStateFlags_is_in_ordered_effect : 0;
			}

			if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_location_source) {
				source_emitters.push_back({ emitter.path_hash, index });
				emitter.spawn_locations_index = tfx__grab_particle_location_lists(pm, e.path_hash, (effect->property_flags & tfxEmitterPropertyFlags_effect_is_3d), 100);
			}
			else if (emitter_properties->emission_type == tfxOtherEmitter) {
				target_emitters.push_back({ emitter_properties->paired_emitter_hash, index });
			}

			/*if (pm->flags & tfxParticleManagerFlags_use_compute_shader && GetEffectInfo(e)->sub_effectors.empty()) {
				int free_slot = AddComputeController();
				if (free_slot != -1) {
					emitter.compute_slot_id = free_slot;
					emitter.property_flags |= tfxEmitterPropertyFlags_is_bottom_emitter;
					tfx_compute_controller_t &c = *(static_cast<tfx_compute_controller_t*>(compute_controller_ptr) + free_slot);
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
						c.image_handle = tfx_vec2_t(0.5f, 0.5f);
					}
					else {
						c.image_handle = properties.image_handle[emitter.property_index];
					}
					c.image_handle = properties.image_handle[emitter.property_index];
				}
			}*/

		}
	}

	if (target_emitters.current_size && source_emitters.current_size) {
		for (hash_index_pair_t target_pair : target_emitters) {
			for (hash_index_pair_t source_pair : source_emitters) {
				if (target_pair.hash == source_pair.hash) {
					pm->emitters[target_pair.index].spawn_locations_index = pm->emitters[source_pair.index].spawn_locations_index;
					break;
				}
			}
		}
	}

	new_effect.state_flags |= tfxEmitterStateFlags_no_tween_this_update;
	return parent_index.index;
}

void UpdateEmitterControlProfile(tfx_effect_emitter_t *emitter) {
	tfx_emitter_properties_t *props = GetEffectProperties(emitter);
	emitter->control_profile = 0;
	if (emitter->property_flags & tfxEmitterPropertyFlags_use_simple_motion_randomness && GetGraphMaxValue(&emitter->library->emitter_attributes[emitter->emitter_attributes].overtime.motion_randomness) > 0.f) {
		emitter->control_profile |= tfxEmitterControlProfile_motion_randomness;
	}
	if (props->emission_direction == tfxOrbital && props->emission_type != tfxPoint) {
		emitter->control_profile |= tfxEmitterControlProfile_orbital;
	}
	if ((!(emitter->property_flags & tfxEmitterPropertyFlags_use_simple_motion_randomness) && GetGraphMaxValue(&emitter->library->emitter_attributes[emitter->emitter_attributes].overtime.velocity_turbulance) && GetGraphMaxValue(&emitter->library->emitter_attributes[emitter->emitter_attributes].overtime.noise_resolution))) {
		emitter->control_profile |= tfxEmitterControlProfile_noise;
	}
	if (props->emission_type == tfxPath) {
		emitter->control_profile |= tfxEmitterControlProfile_path;
	}
	if (emitter->property_flags & tfxEmitterPropertyFlags_edge_traversal && (props->emission_type == tfxPath || props->emission_type == tfxLine)) {
		emitter->control_profile |= tfxEmitterControlProfile_edge_traversal;
		if (props->end_behaviour == tfxLoop) {
			emitter->control_profile |= tfxEmitterControlProfile_edge_loop;
		}
		else if (props->end_behaviour == tfxKill) {
			emitter->control_profile |= tfxEmitterControlProfile_edge_kill;
		}
	}
	if (GetGraphMaxValue(&emitter->library->emitter_attributes[emitter->emitter_attributes].overtime.stretch)) {
		emitter->control_profile |= tfxEmitterControlProfile_stretch;
	}
}

int AddComputeController(tfx_particle_manager_t *pm) {
	//Compute slots should only ever be added for the bottom emitter that has no sub effects
	unsigned int free_slot;
	if (!pm->free_compute_controllers.empty()) {
		free_slot = pm->free_compute_controllers.pop_back();
	}
	else {
		free_slot = pm->highest_compute_controller_index++;
	}
	if (free_slot >= pm->max_compute_controllers)
		return -1;
	return free_slot;
}

void ResetParticlePtr(tfx_particle_manager_t *pm, void *ptr) {
	pm->new_compute_particle_ptr = ptr;
	pm->new_compute_particle_index = 0;
}

void ResetControllerPtr(tfx_particle_manager_t *pm, void *ptr) {
	pm->compute_controller_ptr = ptr;
}

void UpdateCompute(tfx_particle_manager_t *pm, void *sampled_particles, unsigned int sample_size) {
	for (int i = 0; i != sample_size; ++i) {
		if (pm->compute_global_state.current_length == 0)
			break;
		tfx_compute_particle_t *sample = static_cast<tfx_compute_particle_t *>(sampled_particles) + i;
		if (sample->age > sample->max_age) {
			pm->compute_global_state.start_index++;
			pm->compute_global_state.start_index %= pm->compute_global_state.end_index;
			pm->compute_global_state.current_length--;
			sample->age = 0;
		}
		else {
			break;
		}
	}
}

tfx_compute_particle_t *GrabComputeParticle(tfx_particle_manager_t *pm, unsigned int layer) {
	TFX_ASSERT(pm->new_compute_particle_ptr);        //Use must assign the compute ptr to point to an area in memory where you can stage new particles for uploading to the GPU - See ResetComputePtr
	return (static_cast<tfx_compute_particle_t *>(pm->new_compute_particle_ptr) + pm->new_compute_particle_index++);
}

void FreeParticleList(tfx_particle_manager_t *pm, tfxU32 index) {
	if (pm->free_particle_lists.ValidKey(pm->emitters[index].path_hash) && pm->emitters[index].particles_index != tfxINVALID) {
		pm->free_particle_lists.At(pm->emitters[index].path_hash).push_back(pm->emitters[index].particles_index);
	}
	else if (pm->emitters[index].particles_index != tfxINVALID) {
		tfx_vector_t<tfxU32> new_indexes;
		new_indexes.push_back(pm->emitters[index].particles_index);
		pm->free_particle_lists.Insert(pm->emitters[index].path_hash, new_indexes);
	}
}

void FreeSpawnLocationList(tfx_particle_manager_t *pm, tfxU32 index) {
	if (pm->free_particle_location_lists.ValidKey(pm->emitters[index].path_hash) && pm->emitters[index].spawn_locations_index != tfxINVALID) {
		ClearSoABuffer(&pm->particle_location_buffers[pm->emitters[index].spawn_locations_index]);
		pm->free_particle_location_lists.At(pm->emitters[index].path_hash).push_back(pm->emitters[index].spawn_locations_index);
	}
	else if (pm->emitters[index].spawn_locations_index != tfxINVALID) {
		tfx_vector_t<tfxU32> new_indexes;
		ClearSoABuffer(&pm->particle_location_buffers[pm->emitters[index].spawn_locations_index]);
		new_indexes.push_back(pm->emitters[index].spawn_locations_index);
		pm->free_particle_location_lists.Insert(pm->emitters[index].path_hash, new_indexes);
	}
}

void OrderEffectSprites(tfx_effect_instance_data_t *sprites, tfxU32 layer, tfx_particle_manager_t *pm) {
	tfxU32 depth_starting_index = sprites->depth_starting_index[layer];
	tfxU32 current_depth_buffer = sprites->current_depth_buffer_index[layer];
	tfx_vector_t<tfx_depth_index_t> &current_depth_indexes = sprites->depth_indexes[layer][current_depth_buffer];
	if (depth_starting_index < current_depth_indexes.current_size) {
		tfxU32 next_depth_buffer = current_depth_buffer ^ 1;
		tfx_vector_t<tfx_depth_index_t> &next_depth_indexes = sprites->depth_indexes[layer][next_depth_buffer];
		if (next_depth_indexes.capacity < current_depth_indexes.capacity) {
			next_depth_indexes.reserve(current_depth_indexes.capacity);
		}
		std::qsort(&current_depth_indexes[depth_starting_index], current_depth_indexes.current_size - depth_starting_index, sizeof(tfx_depth_index_t), tfx__sort_depth);
		tfxU32 current_depth_index = 0;
		tfxU32 second_index = depth_starting_index;
		for (auto &depth_index : current_depth_indexes) {
			if (depth_starting_index != 0) {
				while (second_index < current_depth_indexes.current_size && depth_index.depth < current_depth_indexes[second_index].depth) {
					pm->particle_arrays[tfx__particle_bank(current_depth_indexes[second_index].particle_id)].depth_index[tfx__particle_index(current_depth_indexes[second_index].particle_id)] = next_depth_indexes.current_size;
					next_depth_indexes.push_back(current_depth_indexes[second_index++]);
				}
			}
			tfxU32 bank = tfx__particle_bank(depth_index.particle_id);
			tfxU32 index = tfx__particle_index(depth_index.particle_id);
			TFX_ASSERT(index < pm->particle_array_buffers[bank].capacity);
			pm->particle_arrays[bank].depth_index[index] = next_depth_indexes.current_size;
			next_depth_indexes.push_back(depth_index);
			if (++current_depth_index == depth_starting_index)
				break;
		}
		if (depth_starting_index != 0 && second_index < current_depth_indexes.current_size) {
			while (second_index < current_depth_indexes.current_size) {
				tfxU32 bank = tfx__particle_bank(current_depth_indexes[second_index].particle_id);
				tfxU32 index = tfx__particle_index(current_depth_indexes[second_index].particle_id);
				TFX_ASSERT(index < pm->particle_array_buffers[bank].capacity);
				pm->particle_arrays[bank].depth_index[index] = next_depth_indexes.current_size;
				next_depth_indexes.push_back(current_depth_indexes[second_index++]);
			}
		}
		TFX_ASSERT(next_depth_indexes.current_size == current_depth_indexes.current_size);
		current_depth_indexes.clear();
		sprites->current_depth_buffer_index[layer] = next_depth_buffer;
	}
}

void UpdateParticleManager(tfx_particle_manager_t *pm, float elapsed_time) {
	tfxPROFILE;

	if (pm->flags & tfxParticleManagerFlags_direct_to_staging_buffer) {
		TFX_ASSERT(pm->instance_buffer.data);	//You must call SetStagingBuffer if flagging the particle manager to write direct to staging buffer
	}
	if(elapsed_time <= 0) return;

	tfxCompleteAllWork(&pm->work_queue);

	if (pm->flags & tfxParticleManagerFlags_use_effect_sprite_buffers && pm->flags & tfxParticleManagerFlags_auto_order_effects) {
		tfx_vector_t<tfx_effect_index_t> &effects_in_use = pm->effects_in_use[0][pm->current_ebuff];
		for (tfxU32 i = 1; i < effects_in_use.current_size; ++i) {
			tfx_effect_index_t key = effects_in_use[i];
			int j = i - 1;
			while (j >= 0 && key.depth > effects_in_use[j].depth) {
				effects_in_use[j + 1] = effects_in_use[j];
				--j;
			}
			effects_in_use[j + 1] = key;
		}
	}

	pm->frame_length = tfx__Min(elapsed_time, pm->max_frame_length);
	pm->frame_length_wide = tfxWideSetSingle(pm->frame_length);
	pm->update_frequency = 1000.f / elapsed_time;
	pm->update_time = 1.f / pm->update_frequency;
	pm->update_time_wide = tfxWideSetSingle(pm->update_time);
	pm->new_compute_particle_index = 0;
	tfxU32 next_buffer = pm->current_ebuff ^ 1;

	//tfxU32 depth_starting_index[tfxLAYERS];
	//for (tfxEachLayer) {
		//depth_starting_index[layer] = pm->depth_indexes[layer][pm->current_depth_buffer_index[layer]].current_size;
	//}
	pm->current_sprite_buffer = pm->flags & tfxParticleManagerFlags_double_buffer_sprites ? pm->current_sprite_buffer ^ 1 : 0;

	memset(pm->sprite_index_point, 0, sizeof(tfxU32) * tfxLAYERS);
	memset(pm->layer_sizes, 0, sizeof(tfxU32) * tfxLAYERS);
	memset(pm->cumulative_index_point, 0, sizeof(tfxU32) * tfxLAYERS);

	pm->instance_buffer.clear();
	pm->highest_depth_index = 0;

	tfxU32 effects_start_size[tfxMAXDEPTH];
	for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
		effects_start_size[depth] = pm->effects_in_use[depth][pm->current_ebuff].current_size;
	}

	pm->control_emitter_queue.clear();

	tfxU32 last_instance_count = 0;

	//Loop over all the effects and emitters, depth by depth, and add spawn jobs to the worker queue
	for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {

		pm->effects_in_use[depth][next_buffer].clear();

		for (int i = 0; i != effects_start_size[depth]; ++i) {
			tfx_effect_index_t &effect_index = pm->effects_in_use[depth][pm->current_ebuff][i];
			tfx_effect_state_t &effect = pm->effects[effect_index.index];
			float &timeout_counter = effect.timeout_counter;
			effect.emitter_indexes[next_buffer].clear();
			effect.emitter_start_size = effect.emitter_indexes[pm->current_ebuff].current_size;

			if (depth == 0) {
				tfx_effect_instance_data_t &instance_data = pm->effects[effect_index.index].instance_data;
				instance_data.instance_start_index = tfxINVALID;
				memset(instance_data.sprite_index_point, 0, sizeof(tfxU32) * tfxLAYERS);
				instance_data.instance_count = 0;
			}

			UpdatePMEffect(pm, effect_index.index);
			if (pm->flags & tfxParticleManagerFlags_auto_order_effects) {
				tfx_vec3_t effect_to_camera = pm->effects[effect_index.index].world_position - pm->camera_position;
				effect_index.depth = (pm->flags & tfxParticleManagerFlags_3d_effects) ? tfx__vec3_length_fast(&effect_to_camera) : effect_index.depth = pm->effects[effect_index.index].world_position.y;
			}
			if (timeout_counter <= pm->effects[effect_index.index].timeout) {
				pm->effects_in_use[depth][next_buffer].push_back(effect_index);
			}
			else {
				pm->free_effects.push_back(effect_index);
			}

			for (int emitter_index : pm->effects[effect_index.index].emitter_indexes[pm->current_ebuff]) {
				//If you hit this assert it means there are more then the default amount of work entries being created for updating particles. You can increase the amount
				//by calling SetPMWorkQueueSizes. It could also hit the limit if you have a small multithreaded_batch_size (set when you created the particle manager) which
				//would cause more work entries to be created.
				TFX_ASSERT(pm->spawn_work.current_size != pm->spawn_work.capacity);
				tfx_spawn_work_entry_t *spawn_work_entry = &pm->spawn_work.next();
				spawn_work_entry->random = pm->threaded_random;
				spawn_work_entry->depth = depth;
				spawn_work_entry->emitter_index = emitter_index;
				spawn_work_entry->next_buffer = next_buffer;
				spawn_work_entry->properties = &pm->library->emitter_properties[pm->emitters[emitter_index].properties_index];
				spawn_work_entry->sub_effects = &pm->library->effect_infos[pm->emitters[emitter_index].info_index].sub_effectors;
				spawn_work_entry->amount_to_spawn = 0;
				spawn_work_entry->highest_particle_age = pm->emitters[emitter_index].highest_particle_age;
				spawn_work_entry->pm = pm;
				spawn_work_entry->depth_indexes = nullptr;

				float &timeout_counter = pm->emitters[emitter_index].timeout_counter;

				UpdatePMEmitter(&pm->work_queue, spawn_work_entry);
				if (timeout_counter <= pm->emitters[emitter_index].timeout) {
					effect.emitter_indexes[next_buffer].push_back(emitter_index);
					pm->control_emitter_queue.push_back(emitter_index);
				} else {
					if (pm->emitters[emitter_index].path_quaternions) {
						FreePathQuaternion(pm, pm->emitters[emitter_index].path_quaternion_index);
					}
					//if (pm->flags & tfxParticleManagerFlags_use_compute_shader && pm->emitters[emitter_index].property_flags & tfxEmitterPropertyFlags_is_bottom_emitter)
						//FreeComputeSlot(pm->emitters[emitter_index].compute_slot_id);
					if (pm->flags & tfxParticleManagerFlags_unordered) {
						FreeParticleList(pm, emitter_index);
					}
					if (pm->emitters[emitter_index].spawn_locations_index != tfxINVALID) {
						FreeSpawnLocationList(pm, emitter_index);
					}
					pm->free_emitters.push_back(emitter_index);
				}
			}

			if (depth == 0) {
				tfx_effect_instance_data_t &instance_data = pm->effects[effect_index.index].instance_data;
				if (!(pm->flags & tfxParticleManagerFlags_recording_sprites)) {
					instance_data.cumulative_index_point[0] = 0;
					instance_data.cumulative_index_point[1] = instance_data.sprite_index_point[0];
					instance_data.cumulative_index_point[2] = instance_data.cumulative_index_point[1] + instance_data.sprite_index_point[1];
					instance_data.cumulative_index_point[3] = instance_data.cumulative_index_point[2] + instance_data.sprite_index_point[2];
				} else {
					memset(instance_data.cumulative_index_point, 0, sizeof(tfxU32) * tfxLAYERS);
				}
				instance_data.instance_start_index = last_instance_count;
				last_instance_count += instance_data.instance_count;
			}
		}
	}

	for (tfx_spawn_work_entry_t *spawn_work : pm->deffered_spawn_work) {
		//We deffer any spawn work to here for any emitters that have ordered effects so that the required buffer space can be calculated
		//before doing any spawning.
		tfxAddWorkQueueEntry(&pm->work_queue, spawn_work, pm->flags & tfxParticleManagerFlags_3d_effects ? DoSpawnWork3d : DoSpawnWork2d);
	}
	pm->deffered_spawn_work.clear();

	tfxCompleteAllWork(&pm->work_queue);
	AdvanceRandom(&pm->threaded_random);

	for (auto &work_entry : pm->spawn_work) {
		tfxU32 index = work_entry.emitter_index;
		pm->emitters[index].highest_particle_age = fmaxf(pm->emitters[index].highest_particle_age, work_entry.highest_particle_age);
		pm->effects[pm->emitters[index].parent_index].highest_particle_age = pm->emitters[index].highest_particle_age + pm->frame_length;
	}
	pm->spawn_work.clear();

	if (!(pm->flags & tfxParticleManagerFlags_recording_sprites)) {
		pm->cumulative_index_point[0] = 0;
		pm->cumulative_index_point[1] = pm->sprite_index_point[0];
		pm->cumulative_index_point[2] = pm->cumulative_index_point[1] + pm->sprite_index_point[1];
		pm->cumulative_index_point[3] = pm->cumulative_index_point[2] + pm->sprite_index_point[2];
	} 

	for (tfx_effect_index_t effect_index : pm->effects_in_use[0][next_buffer]) {
		tfx_effect_instance_data_t &sprites = pm->effects[effect_index.index].instance_data;
		if (IsOrderedEffectState(&pm->effects[effect_index.index])) {
			for (tfxEachLayer) {
				OrderEffectSprites(&sprites, layer, pm);
			}
		}
	}

	bool is_recording = (pm->flags & tfxParticleManagerFlags_recording_sprites) > 0 && (pm->flags & tfxParticleManagerFlags_using_uids) > 0;
	for (int index : pm->control_emitter_queue) {
		tfx_soa_buffer_t &bank = pm->particle_array_buffers[pm->emitters[index].particles_index];
		int particles_to_update = bank.current_size;
		tfxU32 running_start_index = 0;
		if (pm->emitters[index].property_flags & tfxEmitterPropertyFlags_spawn_location_source && pm->emitters[index].spawn_locations_index != tfxINVALID) {
			ClearSoABuffer(&pm->particle_location_buffers[pm->emitters[index].spawn_locations_index]);
		}
		while (particles_to_update > 0) {
			//If you hit this assert it means there are more then the default amount of work entries being created for updating particles. You can increase the amount
			//by calling SetPMWorkQueueSizes. It could also hit the limit if you have a small multithreaded_batch_size (set when you created the particle manager) which
			//would cause more work entries to be created.
			TFX_ASSERT(pm->control_work.current_size != pm->control_work.capacity);
			tfx_control_work_entry_t &work_entry = pm->control_work.next();
			work_entry.properties = &pm->library->emitter_properties[pm->emitters[index].properties_index];
			work_entry.pm = pm;
			work_entry.emitter_index = index;
			work_entry.start_index = running_start_index;
			work_entry.end_index = particles_to_update > pm->mt_batch_size ? running_start_index + pm->mt_batch_size : running_start_index + particles_to_update;
			work_entry.sprite_instances = !is_recording ? &pm->instance_buffer : &pm->instance_buffer_for_recording[pm->current_sprite_buffer][work_entry.properties->layer];
			tfxU32 circular_start = GetCircularIndex(&pm->particle_array_buffers[pm->emitters[index].particles_index], work_entry.start_index);
			tfxU32 block_start_index = (circular_start / tfxDataWidth) * tfxDataWidth;
			work_entry.wide_end_index = (tfxU32)(ceilf((float)work_entry.end_index / tfxDataWidth)) * tfxDataWidth;
			work_entry.start_diff = circular_start - block_start_index;
			work_entry.wide_end_index = work_entry.wide_end_index - work_entry.start_diff < work_entry.end_index ? work_entry.wide_end_index + tfxDataWidth : work_entry.wide_end_index;
			tfx_effect_state_t &parent_effect = pm->effects[pm->emitters[index].parent_index];
			work_entry.global_stretch = parent_effect.stretch;
			work_entry.global_noise = parent_effect.noise;
			work_entry.global_intensity = parent_effect.spawn_controls.intensity;
			particles_to_update -= pm->mt_batch_size;
			running_start_index += pm->mt_batch_size;
			tfxAddWorkQueueEntry(&pm->work_queue, &work_entry, ControlParticles);
		}
	}

	tfxCompleteAllWork(&pm->work_queue);
	pm->control_work.clear();

	{
		for (int index : pm->control_emitter_queue) {
			tfx_soa_buffer_t &bank = pm->particle_array_buffers[pm->emitters[index].particles_index];
			//If you hit this assert it means there are more then the default amount of work entries being created for updating particles. You can increase the amount
			//by calling SetPMWorkQueueSizes. It could also hit the limit if you have a small multithreaded_batch_size (set when you created the particle manager) which
			//would cause more work entries to be created.
			TFX_ASSERT(pm->age_work.current_size != pm->age_work.capacity);
			tfx_particle_age_work_entry_t &work_entry = pm->age_work.next();
			work_entry.properties = &pm->library->emitter_properties[pm->emitters[index].properties_index];
			work_entry.start_index = bank.current_size - 1;
			work_entry.emitter_index = index;
			tfxU32 circular_start = GetCircularIndex(&pm->particle_array_buffers[pm->emitters[index].particles_index], 0);
			tfxU32 block_start_index = (circular_start / tfxDataWidth) * tfxDataWidth;
			work_entry.wide_end_index = (tfxU32)(ceilf((float)bank.current_size / tfxDataWidth)) * tfxDataWidth;
			work_entry.start_diff = circular_start - block_start_index;
			work_entry.wide_end_index += work_entry.wide_end_index - work_entry.start_diff < bank.current_size ? tfxDataWidth : 0;
			work_entry.pm = pm;
			if (!(pm->flags & tfxParticleManagerFlags_single_threaded) && tfxNumberOfThreadsInAdditionToMain) {
				tfxAddWorkQueueEntry(&pm->work_queue, &work_entry, ControlParticleAge);
			}
			else {
				ControlParticleAge(&pm->work_queue, &work_entry);
			}
		}
	}
	tfxCompleteAllWork(&pm->work_queue);
	pm->age_work.clear();

	//Todo work queue this for each layer
	/*
	if (!(pm->flags & tfxParticleManagerFlags_use_effect_sprite_buffers) && !(pm->flags & tfxParticleManagerFlags_unordered)) {
		for (tfxEachLayer) {
			for (auto &depth_index : pm->depth_indexes[layer][pm->current_depth_buffer_index[layer]]) {
				if (depth_index.particle_id != tfxINVALID) {
					pm->particle_arrays[tfx__particle_bank(depth_index.particle_id)].depth_index[tfx__particle_index(depth_index.particle_id)] = pm->depth_indexes[layer][!pm->current_depth_buffer_index[layer]].current_size;
					pm->depth_indexes[layer][!pm->current_depth_buffer_index[layer]].push_back(depth_index);
				}
			}
			pm->depth_indexes[layer][pm->current_depth_buffer_index[layer]].clear();
			pm->current_depth_buffer_index[layer] = pm->current_depth_buffer_index[layer] ^ 1;
		}
	}
	else if (pm->flags & tfxParticleManagerFlags_use_effect_sprite_buffers) {
	*/
		for (tfx_effect_index_t effect_index : pm->effects_in_use[0][next_buffer]) {
			tfx_effect_instance_data_t &sprites = pm->effects[effect_index.index].instance_data;
			for (tfxEachLayer) {
				for (auto &depth : sprites.depth_indexes[layer][sprites.current_depth_buffer_index[layer]]) {
					if (depth.particle_id != tfxINVALID) {
						tfxU32 bank = tfx__particle_bank(depth.particle_id);
						tfxU32 index = tfx__particle_index(depth.particle_id);
						pm->particle_arrays[bank].depth_index[index] = sprites.depth_indexes[layer][sprites.current_depth_buffer_index[layer] ^ 1].current_size;
						sprites.depth_indexes[layer][sprites.current_depth_buffer_index[layer] ^ 1].push_back(depth);
					}
				}
				sprites.depth_indexes[layer][sprites.current_depth_buffer_index[layer]].clear();
				sprites.current_depth_buffer_index[layer] = sprites.current_depth_buffer_index[layer] ^ 1;
			}
		}
	//}

	/*
	if (!(pm->flags & tfxParticleManagerFlags_use_effect_sprite_buffers) && pm->flags & tfxParticleManagerFlags_order_by_depth && pm->flags & tfxParticleManagerFlags_3d_effects) {
		pm->sorting_work_entry.clear();
		if (pm->flags & tfxParticleManagerFlags_guarantee_order) {
			for (tfxEachLayer) {
				tfx_sort_work_entry_t &work_entry = pm->sorting_work_entry.next();
				work_entry.bank = &pm->particle_arrays;
				work_entry.depth_indexes = &pm->depth_indexes[layer][pm->current_depth_buffer_index[layer]];
				if (!(pm->flags & tfxParticleManagerFlags_single_threaded) && tfxNumberOfThreadsInAdditionToMain > 0) {
					tfxAddWorkQueueEntry(&pm->work_queue, &work_entry, tfx__insertion_sort_depth);
				}
				else {
					tfx__insertion_sort_depth(&pm->work_queue, &work_entry);
				}
			}
		}
		else if (pm->sort_passes > 0) {
			for (tfxEachLayer) {
				tfx_vector_t<tfx_depth_index_t> &depth_index = pm->depth_indexes[layer][pm->current_depth_buffer_index[layer]];
				//Add this to a work queue
				for (tfxU32 sorts = 0; sorts != pm->sort_passes; ++sorts) {
					for (tfxU32 i = 1; i < depth_index.current_size; ++i) {
						float depth1 = depth_index[i - 1].depth;
						float depth2 = depth_index[i].depth;
						if (depth1 < depth2) {
							pm->particle_arrays[tfx__particle_bank(depth_index[i].particle_id)].depth_index[tfx__particle_index(depth_index[i].particle_id)] = i - 1;
							pm->particle_arrays[tfx__particle_bank(depth_index[i - 1].particle_id)].depth_index[tfx__particle_index(depth_index[i - 1].particle_id)] = i;
							std::swap(depth_index[i], depth_index[i - 1]);
						}
					}
				}
			}
		}
	}
	else if (pm->flags & tfxParticleManagerFlags_use_effect_sprite_buffers) {
	*/
		pm->sorting_work_entry.clear();
		for (tfx_effect_index_t effect_index : pm->effects_in_use[0][next_buffer]) {
			tfx_effect_state_t &effect = pm->effects[effect_index.index];
			if (effect.effect_flags & tfxEffectPropertyFlags_depth_draw_order) {
				if (effect.effect_flags & tfxEffectPropertyFlags_guaranteed_order) {
					for (tfxEachLayer) {
						tfx_sort_work_entry_t &work_entry = pm->sorting_work_entry.next();
						work_entry.bank = &pm->particle_arrays;
						work_entry.depth_indexes = &effect.instance_data.depth_indexes[layer][effect.instance_data.current_depth_buffer_index[layer]];
						if (!(pm->flags & tfxParticleManagerFlags_single_threaded) && tfxNumberOfThreadsInAdditionToMain > 0) {
							tfxAddWorkQueueEntry(&pm->work_queue, &work_entry, tfx__insertion_sort_depth);
						}
						else {
							tfx__insertion_sort_depth(&pm->work_queue, &work_entry);
						}
					}
				}
				else if (effect.sort_passes > 0) {
					for (tfxEachLayer) {
						tfx_vector_t<tfx_depth_index_t> &depth_index = effect.instance_data.depth_indexes[layer][effect.instance_data.current_depth_buffer_index[layer]];
						//Add this to a work queue
						for (tfxU32 sorts = 0; sorts != pm->sort_passes; ++sorts) {
							for (tfxU32 i = 1; i < depth_index.current_size; ++i) {
								float depth1 = depth_index[i - 1].depth;
								float depth2 = depth_index[i].depth;
								if (depth1 < depth2) {
									pm->particle_arrays[tfx__particle_bank(depth_index[i].particle_id)].depth_index[tfx__particle_index(depth_index[i].particle_id)] = i - 1;
									pm->particle_arrays[tfx__particle_bank(depth_index[i - 1].particle_id)].depth_index[tfx__particle_index(depth_index[i - 1].particle_id)] = i;
									std::swap(depth_index[i], depth_index[i - 1]);
								}
							}
						}
					}
				}
			}
		}
	//}

	//Add Subeffects to the next buffer
	for (int depth = 1; depth != tfxMAXDEPTH; ++depth) {
		for (int i = effects_start_size[depth]; i != pm->effects_in_use[depth][pm->current_ebuff].current_size; ++i) {
			tfx_effect_index_t current_index = pm->effects_in_use[depth][pm->current_ebuff][i];
			pm->effects_in_use[depth][next_buffer].push_back(current_index);
			tfx_effect_state_t &effect = pm->effects[current_index.index];
			for (int i = effect.emitter_start_size; i != effect.emitter_indexes[pm->current_ebuff].current_size; ++i) {
				tfxU32 emitter_index = effect.emitter_indexes[pm->current_ebuff][i];
				//Make sure to grab a particle list for the sub effect emitters as this doesn't happen when calling AddEffectToParticleManager
				pm->emitters[emitter_index].particles_index = tfx__grab_particle_lists(pm, pm->emitters[emitter_index].path_hash, pm->flags & tfxParticleManagerFlags_3d_effects, 100, pm->emitters[emitter_index].control_profile);
				effect.emitter_indexes[next_buffer].push_back(emitter_index);
			}
		}
	}

	pm->current_ebuff = next_buffer;

	if (pm->flags & tfxParticleManagerFlags_update_bounding_boxes) {
		for (int i = 0; i != pm->effects_in_use[0][pm->current_ebuff].size(); ++i) {
			tfx_effect_index_t effect_index = pm->effects_in_use[0][pm->current_ebuff][i];
			tfx_effect_state_t &effect = pm->effects[effect_index.index];
			tfx_bounding_box_t &effect_bb = effect.bounding_box;
			effect_bb.max_corner.x = -FLT_MAX;
			effect_bb.max_corner.y = -FLT_MAX;
			effect_bb.max_corner.z = -FLT_MAX;
			effect_bb.min_corner.x = FLT_MAX;
			effect_bb.min_corner.y = FLT_MAX;
			effect_bb.min_corner.z = FLT_MAX;
		}
		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			for (int i = 0; i != pm->effects_in_use[depth][pm->current_ebuff].size(); ++i) {
				tfx_effect_index_t effect_index = pm->effects_in_use[0][pm->current_ebuff][i];
				tfx_effect_state_t &effect = pm->effects[effect_index.index];
				for (tfxU32 emitter_index : effect.emitter_indexes[pm->current_ebuff]) {
					tfx_emitter_state_t &emitter = pm->emitters[emitter_index];
					tfx_bounding_box_t &effect_bb = pm->effects[emitter.root_index].bounding_box;
					effect_bb.max_corner.x = tfx__Max(effect_bb.max_corner.x, emitter.bounding_box.max_corner.x);
					effect_bb.max_corner.y = tfx__Max(effect_bb.max_corner.y, emitter.bounding_box.max_corner.y);
					effect_bb.max_corner.z = tfx__Max(effect_bb.max_corner.z, emitter.bounding_box.max_corner.z);
					effect_bb.min_corner.x = tfx__Min(effect_bb.min_corner.x, emitter.bounding_box.min_corner.x);
					effect_bb.min_corner.y = tfx__Min(effect_bb.min_corner.y, emitter.bounding_box.min_corner.y);
					effect_bb.min_corner.z = tfx__Min(effect_bb.min_corner.z, emitter.bounding_box.min_corner.z);
				}
			}
		}
	}

	pm->flags &= ~tfxParticleManagerFlags_update_base_values;
}

#define tfxParticleNoise2dLoopUnroll(n)        \
x4 = tfx128SetSingle(x.a[n]);    \
y4 = tfx128SetSingle(y.a[n]);    \
xeps4 = tfx128Set(x.a[n] - eps, x.a[n] + eps, x.a[n], x.a[n]);    \
sample = tfxNoise4_2d(xeps4, y4);    \
a = (sample.a[0] - sample.a[1]) / eps2;    \
b = (sample.a[2] - sample.a[3]) / eps2;    \
noise_x.a[n] = a - b;    \
y.a[n] += 100.f;    \
yeps4r = tfx128Set(y.a[n] - eps, y.a[n] + eps, y.a[n], y.a[n]);    \
sample = tfxNoise4_2d(x4, yeps4r);    \
a = (sample.a[0] - sample.a[1]) / eps2;    \
b = (sample.a[2] - sample.a[3]) / eps2;    \
noise_y.a[n] = a - b;    \

void ControlParticlePositionPath2d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfxU32 emitter_index = work_entry->emitter_index;
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[emitter_index];
	tfx_particle_soa_t &bank = work_entry->pm->particle_arrays[emitter.particles_index];

	//There must be a path setup for the emitter.
	TFX_ASSERT(emitter.path_attributes != tfxINVALID);
	tfx_emitter_path_t *path = &pm.library->paths[emitter.path_attributes];
	tfxWideFloat node_count = tfxWideSetSingle((float)path->node_count - 3.f);
	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
	tfxWideFloat life;
	const tfxWideInt velocity_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity.lookup.last_frame);
	const tfxWideFloat emitter_width = tfxWideSetSingle(emitter.emitter_size.x);
	const tfxWideFloat emitter_height = tfxWideSetSingle(emitter.emitter_size.y);
	const tfxWideFloat e_handle_x = tfxWideSetSingle(emitter.handle.x);
	const tfxWideFloat e_handle_y = tfxWideSetSingle(emitter.handle.y);
	const tfxWideFloat velocity_adjuster = tfxWideSetSingle(lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].overtime.velocity_adjuster, emitter.frame));
	const tfxWideInt velocity_turbulance_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity_turbulance.lookup.last_frame);
	const tfxWideInt noise_resolution_last_frame = tfxWideSetSinglei(work_entry->graphs->noise_resolution.lookup.last_frame);
	const tfxWideFloat overal_scale_wide = tfxWideSetSingle(work_entry->overal_scale);

	const tfxWideInt capture_after_transform_flag = tfxWideSetSinglei(tfxParticleFlags_capture_after_transform);
	const tfxWideInt stretch_last_frame = tfxWideSetSinglei(work_entry->graphs->stretch.lookup.last_frame);
	const tfxWideFloat stretch = tfxWideSetSingle(work_entry->global_stretch);
	tfxWideArray p_stretch;
	tfxU32 start_diff = work_entry->start_diff;
	const float eps = 0.0001f;
	const float eps2 = 0.0002f;
	tfxU32 running_sprite_index = work_entry->sprites_index;
	tfx_2d_instance_t *sprites = tfxCastBuffer(tfx_2d_instance_t, work_entry->sprite_instances);

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;

		tfxWideFloat local_position_x = tfxWideLoad(&bank.position_x[index]);
		tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
		tfxWideFloat path_position = tfxWideLoad(&bank.path_position[index]);
		tfxWideFloat path_offset = tfxWideLoad(&bank.path_offset[index]);
		const tfxWideFloat base_velocity = tfxWideLoad(&bank.base_velocity[index]);
		tfxWideArrayi lookup_frame;

		//----Velocity Changes

		const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
		const tfxWideFloat age = tfxWideLoad(&bank.age[index]);

		tfx__readbarrier;

		if (emitter.property_flags & tfxEmitterPropertyFlags_alt_velocity_lifetime_sampling) {
			life = tfxWideDiv(path_position, node_count);
			life = tfxWideMul(life, max_life);
		}
		else {
			life = tfxWideDiv(age, max_age);
			life = tfxWideMul(life, max_life);
		}

		lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_last_frame);
		const tfxWideFloat lookup_velocity = tfxWideLookupSet(work_entry->graphs->velocity.lookup.values, lookup_frame);
		tfxWideFloat velocity_scalar = tfxWideMul(tfxWideMul(base_velocity, lookup_velocity), velocity_adjuster);
		path_position = tfxWideAdd(path_position, tfxWideMul(velocity_scalar, pm.update_time_wide));

		tfxWideInt flags = tfxWideLoadi((tfxWideIntLoader *)&bank.flags[index]);
		if (emitter.state_flags & tfxEmitterStateFlags_kill) {
			//Kill if the particle has reached the end of the path
			tfxWideInt remove_flag = tfxWideSetSinglei(tfxParticleFlags_remove);
			tfxWideInt remove_flags = tfxWideAndi(remove_flag, tfxWideOri(tfxWideCasti(tfxWideLess(path_position, tfxWideSetZero)), tfxWideCasti(tfxWideGreaterEqual(path_position, node_count))));
			path_position = tfxWideMax(path_position, tfxWideSetZero);
			flags = tfxWideOri(flags, remove_flags);
			tfxWideStorei((tfxWideIntLoader *)&bank.flags[index], flags);
		}
		else {
			//Reposition if the particle is travelling along the path
			tfxWideFloat at_end = tfxWideGreaterEqual(path_position, node_count);
			path_position = tfxWideSub(path_position, tfxWideAnd(at_end, node_count));
			flags = tfxWideOri(flags, tfxWideAndi(capture_after_transform_flag, tfxWideCasti(at_end)));
			at_end = tfxWideLess(path_position, tfxWideSetZero);
			path_position = tfxWideAdd(path_position, tfxWideAnd(at_end, node_count));
			flags = tfxWideOri(flags, tfxWideAndi(capture_after_transform_flag, tfxWideCasti(at_end)));
			tfxWideStorei((tfxWideIntLoader *)&bank.flags[index], flags);
		}

		tfxWideArrayi node_index;
		node_index.m = tfxWideConverti(path_position);
		tfxWideArray t;
		t.m = tfxWideSub(path_position, tfxWideConvert(node_index.m));
		tfxWideArray point_x;
		tfxWideArray point_y;
		tfx__wide_catmull_rom_spline_2d(&node_index, t.m, path->node_soa.x, path->node_soa.y, &point_x.m, &point_y.m);
		tfxWideFloat last_position_x = local_position_x;
		tfxWideFloat last_position_y = local_position_y;
		if (path->extrusion_type == tfxExtrusionArc) {
			tfxWideFloat radius = tfxWideAdd(tfxWideMul(point_x.m, point_x.m), tfxWideMul(point_y.m, point_y.m));
			tfxWideFloat length_mask = tfxWideGreater(radius, tfxWideSetZero);
			radius = tfxWideMul(tfxWideRSqrt(radius), radius);
			tfxWideArray angle;
			tfxWideArray rx;
			tfxWideArray ry;
			angle.m = tfxWideAtan2(point_y.m, point_x.m);
			angle.m = tfxWideAnd(length_mask, angle.m);
			angle.m = tfxWideAdd(angle.m, path_offset);
			tfxWideSinCos(angle.m, &ry.m, &rx.m);
			local_position_x = tfxWideMul(rx.m, radius);
			local_position_y = tfxWideMul(ry.m, radius);
		}
		else {
			local_position_x = tfxWideAdd(point_x.m, path_offset);
			local_position_y = point_y.m;
		}

		local_position_x = tfxWideAdd(local_position_x, e_handle_x);
		local_position_y = tfxWideAdd(local_position_y, e_handle_y);

		//Emission rotation
		if (emitter.state_flags & tfxEmitterStateFlags_has_rotated_path) {
			tfxWideInt quaternion = tfxWideLoadi((tfxWideIntLoader *)&bank.quaternion[index]);
			tfx__wide_transform_packed_quaternion_vec2(&quaternion, &local_position_x, &local_position_y);
		}
		//---

		if (emitter.control_profile & tfxEmitterControlProfile_noise) {
			tfxWideArray noise_x, noise_y;

			const tfxWideFloat noise_resolution = tfxWideLoad(&bank.noise_resolution[index]);
			const tfxWideFloat noise_offset = tfxWideLoad(&bank.noise_offset[index]);

			tfx__readbarrier;

			lookup_frame.m = tfxWideMini(tfxWideConverti(life), noise_resolution_last_frame);
			const tfxWideFloat lookup_noise_resolution = tfxWideMul(tfxWideLookupSet(work_entry->graphs->noise_resolution.lookup.values, lookup_frame), noise_resolution);
			lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_turbulance_last_frame);
			const tfxWideFloat lookup_velocity_turbulance = tfxWideLookupSet(work_entry->graphs->velocity_turbulance.lookup.values, lookup_frame);

			tfxWideArray x, y;
			x.m = tfxWideAdd(tfxWideDiv(local_position_x, lookup_noise_resolution), noise_offset);
			y.m = tfxWideAdd(tfxWideDiv(local_position_y, lookup_noise_resolution), noise_offset);

			tfx128 x4, y4, xeps4, yeps4r;
			float a, b;
			tfx128Array sample;

			tfxParticleNoise2dLoopUnroll(0);
			tfxParticleNoise2dLoopUnroll(1);
			tfxParticleNoise2dLoopUnroll(2);
			tfxParticleNoise2dLoopUnroll(3);

#if defined(tfxUSEAVX)
			tfxParticleNoise2dLoopUnroll(4);
			tfxParticleNoise2dLoopUnroll(5);
			tfxParticleNoise2dLoopUnroll(6);
			tfxParticleNoise2dLoopUnroll(7);
#endif

			noise_x.m = tfxWideMul(lookup_velocity_turbulance, noise_x.m);
			noise_y.m = tfxWideMul(lookup_velocity_turbulance, noise_y.m);

			local_position_x = tfxWideAdd(local_position_x, noise_x.m);
			local_position_y = tfxWideAdd(local_position_y, noise_y.m);
		}

		local_position_x = tfxWideMul(local_position_x, emitter_width);
		local_position_y = tfxWideMul(local_position_y, emitter_height);
		tfxWideStore(&bank.position_x[index], local_position_x);
		tfxWideStore(&bank.position_y[index], local_position_y);

		tfxWideStore(&bank.path_position[index], path_position);

		if (emitter.control_profile & tfxEmitterControlProfile_stretch) {
			//Stretch
			lookup_frame.m = tfxWideMini(tfxWideConverti(life), stretch_last_frame);
			const tfxWideFloat lookup_stretch = tfxWideLookupSet(work_entry->graphs->stretch.lookup.values, lookup_frame);
			p_stretch.m = tfxWideMul(lookup_stretch, stretch);

			tfxWideFloat stretch_velocity_x;
			tfxWideFloat stretch_velocity_y;

			stretch_velocity_x = tfxWideSub(local_position_x, last_position_x);
			stretch_velocity_y = tfxWideSub(local_position_y, last_position_y);

			stretch_velocity_y = tfxWideAdd(stretch_velocity_y, tfxWideSetSingle(0.000001f));
			tfxWideFloat l = tfxWideMul(stretch_velocity_x, stretch_velocity_x);
			l = tfxWideAdd(l, tfxWideMul(stretch_velocity_y, stretch_velocity_y));
			l = tfxWideMul(tfxWideRSqrt(l), l);

			p_stretch.m = tfxWideMul(p_stretch.m, tfxWideMul(velocity_scalar, tfxWideSetSingle(0.015f)));
			stretch_velocity_x = tfxWideDiv(stretch_velocity_x, l);
			stretch_velocity_y = tfxWideDiv(stretch_velocity_y, l);

			if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position) {
				tfx__wide_transform_quaternion_vec2(&emitter.rotation, &stretch_velocity_x, &stretch_velocity_y);
			}

			tfxWideArrayi packed;
			packed.m = tfx__wide_pack16bit(tfxWideMax(tfxWIDEMINUSONE, tfxWideMin(tfxWIDEONE, stretch_velocity_x)), tfxWideMax(tfxWIDEMINUSONE, tfxWideMin(tfxWIDEONE, stretch_velocity_y)));

			tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
			if (!(pm.flags & tfxParticleManagerFlags_unordered)) {    //Predictable
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 sprite_depth_index = bank.depth_index[index + j] + work_entry->cumulative_index_point + work_entry->effect_instance_offset;
					sprites[sprite_depth_index].position.z = p_stretch.a[j];
					sprites[sprite_depth_index].alignment.packed = packed.a[j];
					running_sprite_index++;
				}
			}
			else {
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					sprites[running_sprite_index].position.z = p_stretch.a[j];
					sprites[running_sprite_index].alignment.packed = packed.a[j];
					running_sprite_index++;
				}
			}
			start_diff = 0;
		}
	}
}

#define tfxParticleNoise3dLoopUnroll(n)        \
x4 = tfx128SetSingle(x.a[n]);    \
y4 = tfx128SetSingle(y.a[n]);    \
z4 = tfx128SetSingle(z.a[n]);    \
xeps4 = tfx128Set(x.a[n] - eps, x.a[n] + eps, x.a[n], x.a[n]);    \
xeps4r = tfx128Set(x.a[n], x.a[n], x.a[n] - eps, x.a[n] + eps);    \
yeps4 = tfx128Set(y.a[n], y.a[n], y.a[n] - eps, y.a[n] + eps);    \
zeps4 = tfx128Set(z.a[n] - eps, z.a[n] + eps, z.a[n], z.a[n]);    \
zeps4r = tfx128Set(z.a[n], z.a[n], z.a[n] - eps, z.a[n] + eps);    \
sample = tfxNoise4_3d(x4, yeps4, zeps4);    \
a = (sample.a[0] - sample.a[1]) / eps2;    \
b = (sample.a[2] - sample.a[3]) / eps2;    \
noise_x.a[n] = a - b;    \
y.a[n] += 100.f;    \
yeps4r = tfx128Set(y.a[n] - eps, y.a[n] + eps, y.a[n], y.a[n]);    \
sample = tfxNoise4_3d(xeps4, y4, zeps4r);    \
a = (sample.a[0] - sample.a[1]) / eps2;    \
b = (sample.a[2] - sample.a[3]) / eps2;    \
noise_y.a[n] = a - b;    \
sample = tfxNoise4_3d(xeps4r, yeps4r, z4);    \
a = (sample.a[0] - sample.a[1]) / eps2;    \
b = (sample.a[2] - sample.a[3]) / eps2;    \
noise_z.a[n] = a - b;

void ControlParticlePositionPath3d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfxU32 emitter_index = work_entry->emitter_index;
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[emitter_index];
	tfx_particle_soa_t &bank = work_entry->pm->particle_arrays[emitter.particles_index];

	//There must be a path setup for the emitter.
	TFX_ASSERT(emitter.path_attributes != tfxINVALID);
	tfx_emitter_path_t *path = &pm.library->paths[emitter.path_attributes];
	tfxWideFloat node_count = tfxWideSetSingle((float)path->node_count - 3.f);
	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
	tfxWideFloat life;
	const tfxWideInt velocity_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity.lookup.last_frame);
	const tfxWideFloat emitter_width = tfxWideSetSingle(emitter.emitter_size.x);
	const tfxWideFloat emitter_height = tfxWideSetSingle(emitter.emitter_size.y);
	const tfxWideFloat emitter_depth = tfxWideSetSingle(emitter.emitter_size.z);
	const tfxWideFloat e_handle_x = tfxWideSetSingle(emitter.handle.x);
	const tfxWideFloat e_handle_y = tfxWideSetSingle(emitter.handle.y);
	const tfxWideFloat e_handle_z = tfxWideSetSingle(emitter.handle.z);
	const tfxWideFloat velocity_adjuster = tfxWideSetSingle(lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].overtime.velocity_adjuster, emitter.frame));
	const tfxWideInt velocity_turbulance_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity_turbulance.lookup.last_frame);
	const tfxWideInt noise_resolution_last_frame = tfxWideSetSinglei(work_entry->graphs->noise_resolution.lookup.last_frame);
	const tfxWideFloat overal_scale_wide = tfxWideSetSingle(work_entry->overal_scale);
	const tfxWideInt capture_after_transform_flag = tfxWideSetSinglei(tfxParticleFlags_capture_after_transform);

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;

		tfxWideFloat local_position_x = tfxWideLoad(&bank.position_x[index]);
		tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
		tfxWideFloat local_position_z = tfxWideLoad(&bank.position_z[index]);
		tfxWideFloat path_position = tfxWideLoad(&bank.path_position[index]);
		tfxWideFloat path_offset = tfxWideLoad(&bank.path_offset[index]);
		const tfxWideFloat base_velocity = tfxWideLoad(&bank.base_velocity[index]);
		tfxWideArrayi lookup_frame;

		//----Velocity Changes

		const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
		const tfxWideFloat age = tfxWideLoad(&bank.age[index]);

		tfx__readbarrier;

		if (emitter.property_flags & tfxEmitterPropertyFlags_alt_velocity_lifetime_sampling) {
			life = tfxWideDiv(path_position, node_count);
			life = tfxWideMul(life, max_life);
		}
		else {
			life = tfxWideDiv(age, max_age);
			life = tfxWideMul(life, max_life);
		}

		lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_last_frame);
		const tfxWideFloat lookup_velocity = tfxWideLookupSet(work_entry->graphs->velocity.lookup.values, lookup_frame);
		tfxWideFloat velocity_scalar = tfxWideMul(tfxWideMul(base_velocity, lookup_velocity), velocity_adjuster);
		path_position = tfxWideAdd(path_position, tfxWideMul(velocity_scalar, pm.update_time_wide));

		tfxWideInt flags = tfxWideLoadi((tfxWideIntLoader *)&bank.flags[index]);
		if (emitter.state_flags & tfxEmitterStateFlags_kill) {
			//Kill if the particle has reached the end of the path
			tfxWideInt remove_flag = tfxWideSetSinglei(tfxParticleFlags_remove);
			tfxWideInt remove_flags = tfxWideAndi(remove_flag, tfxWideOri(tfxWideCasti(tfxWideLess(path_position, tfxWideSetZero)), tfxWideCasti(tfxWideGreaterEqual(path_position, node_count))));
			path_position = tfxWideMax(path_position, tfxWideSetZero);
			flags = tfxWideOri(flags, remove_flags);
			tfxWideStorei((tfxWideIntLoader *)&bank.flags[index], flags);
		}
		else {
			//Reposition if the particle is travelling along the path
			tfxWideFloat at_end = tfxWideGreaterEqual(path_position, node_count);
			path_position = tfxWideSub(path_position, tfxWideAnd(at_end, node_count));
			flags = tfxWideOri(flags, tfxWideAndi(capture_after_transform_flag, tfxWideCasti(at_end)));
			at_end = tfxWideLess(path_position, tfxWideSetZero);
			path_position = tfxWideAdd(path_position, tfxWideAnd(at_end, node_count));
			flags = tfxWideOri(flags, tfxWideAndi(capture_after_transform_flag, tfxWideCasti(at_end)));
			tfxWideStorei((tfxWideIntLoader *)&bank.flags[index], flags);
		}

		tfxWideArrayi node_index;
		node_index.m = tfxWideConverti(path_position);
		tfxWideArray t;
		t.m = tfxWideSub(path_position, tfxWideConvert(node_index.m));
		tfxWideArray point_x;
		tfxWideArray point_z;
		tfx__wide_catmull_tom_spline_3d(&node_index, t.m, path->node_soa.x, path->node_soa.y, path->node_soa.z, &point_x.m, &local_position_y, &point_z.m);
		if (path->extrusion_type == tfxExtrusionArc) {
			tfxWideFloat radius = tfxWideAdd(tfxWideMul(point_x.m, point_x.m), tfxWideMul(point_z.m, point_z.m));
			tfxWideFloat length_mask = tfxWideGreater(radius, tfxWideSetZero);
			radius = tfxWideMul(tfxWideRSqrt(radius), radius);
			tfxWideArray angle;
			tfxWideArray rx;
			tfxWideArray rz;
			angle.m = tfxWideAtan2(point_z.m, point_x.m);
			angle.m = tfxWideAnd(length_mask, angle.m);
			angle.m = tfxWideAdd(angle.m, path_offset);
			tfxWideSinCos(angle.m, &rz.m, &rx.m);
			local_position_x = tfxWideMul(rx.m, radius);
			local_position_z = tfxWideMul(rz.m, radius);
		}
		else {
			local_position_x = tfxWideAdd(point_x.m, path_offset);
			local_position_z = point_z.m;
		}

		local_position_x = tfxWideAdd(local_position_x, e_handle_x);
		local_position_y = tfxWideAdd(local_position_y, e_handle_y);
		local_position_z = tfxWideAdd(local_position_z, e_handle_z);

		//Emission rotation
		if (emitter.state_flags & tfxEmitterStateFlags_has_rotated_path) {
			tfxWideInt quaternion = tfxWideLoadi((tfxWideIntLoader *)&bank.quaternion[index]);
			tfx__wide_transform_packed_quaternion_vec3(&quaternion, &local_position_x, &local_position_y, &local_position_z);
		}
		//---

		if (emitter.control_profile & tfxEmitterControlProfile_noise) {
			tfxWideArray noise_x;
			tfxWideArray noise_y;
			tfxWideArray noise_z;

			float eps = 0.001f;
			float eps2 = 0.001f * 2.f;
			const tfxWideFloat noise_resolution = tfxWideLoad(&bank.noise_resolution[index]);
			const tfxWideFloat base_noise_offset = tfxWideLoad(&bank.noise_offset[index]);
			tfxWideFloat noise_offset = tfxWideMul(base_noise_offset, overal_scale_wide);

			tfx__readbarrier;

			tfxWideArrayi lookup_frame;
			lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_turbulance_last_frame);
			const tfxWideFloat lookup_velocity_turbulance = tfxWideLookupSet(work_entry->graphs->velocity_turbulance.lookup.values, lookup_frame);

			lookup_frame.m = tfxWideMini(tfxWideConverti(life), noise_resolution_last_frame);
			const tfxWideFloat lookup_noise_resolution = tfxWideMul(tfxWideLookupSet(work_entry->graphs->noise_resolution.lookup.values, lookup_frame), noise_resolution);

			tfxWideArray x, y, z;
			x.m = tfxWideAdd(tfxWideDiv(local_position_x, lookup_noise_resolution), noise_offset);
			y.m = tfxWideAdd(tfxWideDiv(local_position_y, lookup_noise_resolution), noise_offset);
			z.m = tfxWideAdd(tfxWideDiv(local_position_z, lookup_noise_resolution), noise_offset);

			tfx128 x4, y4, z4, xeps4, xeps4r, yeps4, zeps4, zeps4r, yeps4r;
			tfx128Array sample;
			float a, b;

			tfxParticleNoise3dLoopUnroll(0)
				tfxParticleNoise3dLoopUnroll(1)
				tfxParticleNoise3dLoopUnroll(2)
				tfxParticleNoise3dLoopUnroll(3)

#if defined(tfxUSEAVX)
				tfxParticleNoise3dLoopUnroll(4)
				tfxParticleNoise3dLoopUnroll(5)
				tfxParticleNoise3dLoopUnroll(6)
				tfxParticleNoise3dLoopUnroll(7)
#endif

				noise_x.m = tfxWideMul(lookup_velocity_turbulance, noise_x.m);
			noise_y.m = tfxWideMul(lookup_velocity_turbulance, noise_y.m);
			noise_z.m = tfxWideMul(lookup_velocity_turbulance, noise_z.m);

			local_position_x = tfxWideAdd(local_position_x, noise_x.m);
			local_position_y = tfxWideAdd(local_position_y, noise_y.m);
			local_position_z = tfxWideAdd(local_position_z, noise_z.m);
		}

		tfxWideStore(&bank.position_x[index], tfxWideMul(local_position_x, emitter_width));
		tfxWideStore(&bank.position_y[index], tfxWideMul(local_position_y, emitter_height));
		tfxWideStore(&bank.position_z[index], tfxWideMul(local_position_z, emitter_depth));

		tfxWideStore(&bank.path_position[index], path_position);
	}
}

//I'm not a fan of this at all, but currently it's the best I can think off to create multiple versions of the control particle function
//without having to repeat a whole bunch of code and make it harder to maintain whilst avoiding funtion callbacks and such.
//This define sets up variables needed in all the control particle position (3d) functions
#define tfxControlParticleLoopSetup    \
tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;    \
tfxWideFloat life;    \
const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);    \
const tfxWideFloat age = tfxWideLoad(&bank.age[index]);    \
tfxWideFloat velocity_normal_x;    \
tfxWideFloat velocity_normal_y;    \
tfxWideFloat velocity_normal_z;    \
tfx__readbarrier;

//Apply the velocity changes and update the particle position
#define tfxControlParticleUpdatePosition \
tfxWideArrayi lookup_frame_weight;    \
lookup_frame_weight.m = tfxWideMini(tfxWideConverti(life), weight_last_frame);    \
const tfxWideFloat lookup_weight = tfxWideLookupSet(work_entry->graphs->weight.lookup.values, lookup_frame_weight);    \
tfxWideFloat age_fraction = tfxWideMin(tfxWideDiv(age, pm.frame_length_wide), tfxWIDEONE);    \
current_velocity_y = tfxWideSub(current_velocity_y, tfxWideMul(base_weight, lookup_weight));    \
current_velocity_x = tfxWideMul(tfxWideMul(tfxWideMul(current_velocity_x, pm.update_time_wide), velocity_adjuster), age_fraction);    \
current_velocity_y = tfxWideMul(tfxWideMul(tfxWideMul(current_velocity_y, pm.update_time_wide), velocity_adjuster), age_fraction);    \
current_velocity_z = tfxWideMul(tfxWideMul(tfxWideMul(current_velocity_z, pm.update_time_wide), velocity_adjuster), age_fraction);    \
local_position_x = tfxWideAdd(local_position_x, tfxWideMul(current_velocity_x, overal_scale_wide));    \
local_position_y = tfxWideAdd(local_position_y, tfxWideMul(current_velocity_y, overal_scale_wide));    \
local_position_z = tfxWideAdd(local_position_z, tfxWideMul(current_velocity_z, overal_scale_wide));    \
tfxWideStore(&bank.position_x[index], local_position_x);    \
tfxWideStore(&bank.position_y[index], local_position_y);    \
tfxWideStore(&bank.position_z[index], local_position_z);    

//Simplex noise
#define tfxControlParticleUpdateNoise \
tfxWideFloat current_velocity_x = tfxWideMul(velocity_normal_x, velocity_scalar);    \
tfxWideFloat current_velocity_y = tfxWideMul(velocity_normal_y, velocity_scalar);    \
tfxWideFloat current_velocity_z = tfxWideMul(velocity_normal_z, velocity_scalar);    \
tfxWideArray noise_x;    \
tfxWideArray noise_y;    \
tfxWideArray noise_z;    \
float eps = 0.001f;    \
float eps2 = 0.001f * 2.f;    \
const tfxWideFloat noise_resolution = tfxWideLoad(&bank.noise_resolution[index]);    \
const tfxWideFloat base_noise_offset = tfxWideLoad(&bank.noise_offset[index]);    \
tfxWideFloat noise_offset = tfxWideMul(base_noise_offset, overal_scale_wide);    \
tfx__readbarrier;    \
lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_turbulance_last_frame);    \
const tfxWideFloat lookup_velocity_turbulance = tfxWideLookupSet(work_entry->graphs->velocity_turbulance.lookup.values, lookup_frame);    \
lookup_frame.m = tfxWideMini(tfxWideConverti(life), noise_resolution_last_frame);    \
const tfxWideFloat lookup_noise_resolution = tfxWideMul(tfxWideLookupSet(work_entry->graphs->noise_resolution.lookup.values, lookup_frame), noise_resolution);    \
tfxWideArray x, y, z;    \
x.m = tfxWideAdd(tfxWideDiv(local_position_x, lookup_noise_resolution), noise_offset);    \
y.m = tfxWideAdd(tfxWideDiv(local_position_y, lookup_noise_resolution), noise_offset);    \
z.m = tfxWideAdd(tfxWideDiv(local_position_z, lookup_noise_resolution), noise_offset);    \
tfx128 x4, y4, z4, xeps4, xeps4r, yeps4, zeps4, zeps4r, yeps4r;    \
tfx128Array sample;    \
float a, b;    

//Sample over life of path
#define tfxControlParticleSampleOverPathLife \
if (work_entry->sample_path_life) {    \
tfxWideFloat path_position = tfxWideLoad(&bank.path_position[index]);    \
life = tfxWideDiv(path_position, node_count);    \
}    \
else {    \
life = tfxWideDiv(age, max_age);    \
}    

#define tfxControlParticleOrbital \
velocity_normal_z = tfxWideSub(local_position_x, emitter_x);    \
velocity_normal_y = tfxWideSetZero;    \
velocity_normal_x = tfxWideMul(tfxWideSub(local_position_z, emitter_z), tfxWideSetSingle(-1.f));    \
tfxWideFloat l = tfxWideMul(velocity_normal_x, velocity_normal_x);    \
l = tfxWideAdd(l, tfxWideMul(velocity_normal_z, velocity_normal_z));    \
l = tfxWideMul(tfxWideRSqrt(l), l);    \
velocity_normal_x = tfxWideDiv(velocity_normal_x, l);    \
velocity_normal_z = tfxWideDiv(velocity_normal_z, l);    

//Simple motion randomness
#define tfxControlParticleMotionRandomness \
tfxWideInt uid = tfxWideLoadi((tfxWideIntLoader*)&bank.uid[index]);    \
tfxWideInt seed = SeedGenWide(time_step, uid);    \
tfxWideFloat speed = tfxWideLoad(&bank.noise_offset[index]);    \
tfxWideArrayi lookup_motion_randomness = { tfxWideMini(tfxWideConverti(life), motion_randomness_last_frame) };    \
const tfxWideFloat influence = tfxWideMul(tfxWideMul(motion_randomness_base, global_noise), tfxWideLookupSet(work_entry->graphs->motion_randomness.lookup.values, lookup_motion_randomness));    \
tfxWideFloat point_one_influence = tfxWideMul(tfxWideSetSingle(0.1f), influence);    \
tfxWideFloat random_speed = tfxWideMul(tfxWideDiv(SeedGenWide(seed), tfxMAXUINTf), tfxWideMul(tfxWideSetSingle(0.01f), influence));    \
tfxWideFloat random_x, random_y, random_z;    \
RandomVectorInConeWide(seed, velocity_normal_x, velocity_normal_y, velocity_normal_z, tfxWideMul(tfxDEGREERANGEMR, influence), &random_x, &random_y, &random_z);    \
speed = tfxWideAdd(speed, random_speed);    \
tfxWideFloat length = tfxWideMul(random_x, random_x);    \
length = tfxWideAdd(length, tfxWideMul(random_y, random_y));    \
length = tfxWideAdd(length, tfxWideMul(random_z, random_z));    \
length = tfxWideMul(tfxWideRSqrt(length), length);    \
tfxWideFloat length_one = tfxWideDiv(tfxWIDEONE, length);    \
random_x = tfxWideMul(random_x, length_one);    \
random_y = tfxWideMul(random_y, length_one);    \
random_z = tfxWideMul(random_z, length_one);    \
velocity_scalar = tfxWideAdd(velocity_scalar, tfxWideMul(speed, global_noise));    

void ControlParticlePosition3dBasic(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfxU32 emitter_index = work_entry->emitter_index;
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[emitter_index];
	tfx_particle_soa_t &bank = work_entry->pm->particle_arrays[emitter.particles_index];
	const tfxWideFloat overal_scale_wide = tfxWideSetSingle(work_entry->overal_scale);
	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
	const tfxWideInt velocity_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity.lookup.last_frame);
	const tfxWideFloat velocity_adjuster = tfxWideSetSingle(lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].overtime.velocity_adjuster, emitter.frame));
	const tfxWideInt weight_last_frame = tfxWideSetSinglei(work_entry->graphs->weight.lookup.last_frame);
	const tfxWideFloat node_count = tfxWideSetSingle(work_entry->node_count);

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxControlParticleLoopSetup;
		tfxControlParticleSampleOverPathLife;

		life = tfxWideMul(life, max_life);
		life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

		tfxWideInt velocity_normal = tfxWideLoadi((tfxWideIntLoader *)&bank.velocity_normal[index]);
		tfx__wide_unpack10bit(velocity_normal, velocity_normal_x, velocity_normal_y, velocity_normal_z);

		const tfxWideFloat base_velocity = tfxWideLoad(&bank.base_velocity[index]);
		const tfxWideFloat base_weight = tfxWideLoad(&bank.base_weight[index]);

		tfxWideArrayi lookup_frame;
		lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_last_frame);
		const tfxWideFloat lookup_velocity = tfxWideLookupSet(work_entry->graphs->velocity.lookup.values, lookup_frame);

		//----Velocity Changes
		tfxWideFloat velocity_scalar = tfxWideMul(base_velocity, lookup_velocity);
		tfxWideFloat current_velocity_x = tfxWideMul(velocity_normal_x, velocity_scalar);
		tfxWideFloat current_velocity_y = tfxWideMul(velocity_normal_y, velocity_scalar);
		tfxWideFloat current_velocity_z = tfxWideMul(velocity_normal_z, velocity_scalar);

		tfxWideFloat local_position_x = tfxWideLoad(&bank.position_x[index]);
		tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
		tfxWideFloat local_position_z = tfxWideLoad(&bank.position_z[index]);

		tfxControlParticleUpdatePosition;
	}
}

void ControlParticlePosition3dOrbital(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfxU32 emitter_index = work_entry->emitter_index;
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[emitter_index];
	tfx_particle_soa_t &bank = work_entry->pm->particle_arrays[emitter.particles_index];
	const tfxWideFloat overal_scale_wide = tfxWideSetSingle(work_entry->overal_scale);
	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
	const tfxWideInt velocity_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity.lookup.last_frame);
	const tfxWideFloat velocity_adjuster = tfxWideSetSingle(lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].overtime.velocity_adjuster, emitter.frame));
	const tfxWideInt weight_last_frame = tfxWideSetSinglei(work_entry->graphs->weight.lookup.last_frame);
	const tfxWideFloat node_count = tfxWideSetSingle(work_entry->node_count);

	tfxWideFloat emitter_x = {};
	tfxWideFloat emitter_z = {};
	if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
		emitter_x = tfxWideSetSingle(emitter.world_position.x);
		emitter_z = tfxWideSetSingle(emitter.world_position.z);
	}
	else if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position) {
		emitter_x = tfxWideSetSingle(emitter.handle.x);
		emitter_z = tfxWideSetSingle(emitter.handle.z);
	}

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxControlParticleLoopSetup;
		tfxControlParticleSampleOverPathLife;

		tfxWideFloat local_position_x = tfxWideLoad(&bank.position_x[index]);
		tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
		tfxWideFloat local_position_z = tfxWideLoad(&bank.position_z[index]);

		//Calculate orbital trajectory that is relative to the emitter
		tfxControlParticleOrbital;
		//---

		life = tfxWideMul(life, max_life);
		life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

		const tfxWideFloat base_velocity = tfxWideLoad(&bank.base_velocity[index]);
		const tfxWideFloat base_weight = tfxWideLoad(&bank.base_weight[index]);

		tfxWideArrayi lookup_frame;
		lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_last_frame);
		const tfxWideFloat lookup_velocity = tfxWideLookupSet(work_entry->graphs->velocity.lookup.values, lookup_frame);

		//----Velocity Changes
		tfxWideFloat velocity_scalar = tfxWideMul(base_velocity, lookup_velocity);
		tfxWideFloat current_velocity_x = tfxWideMul(velocity_normal_x, velocity_scalar);
		tfxWideFloat current_velocity_y = tfxWideMul(velocity_normal_y, velocity_scalar);
		tfxWideFloat current_velocity_z = tfxWideMul(velocity_normal_z, velocity_scalar);

		tfxControlParticleUpdatePosition;
	}

}

//Used for emitters that have simplex noise only.
void ControlParticlePosition3dNoise(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfxU32 emitter_index = work_entry->emitter_index;
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[emitter_index];
	tfx_particle_soa_t &bank = work_entry->pm->particle_arrays[emitter.particles_index];
	const tfxWideFloat overal_scale_wide = tfxWideSetSingle(work_entry->overal_scale);
	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
	const tfxWideInt velocity_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity.lookup.last_frame);
	const tfxWideFloat velocity_adjuster = tfxWideSetSingle(lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].overtime.velocity_adjuster, emitter.frame));
	const tfxWideFloat global_noise = tfxWideSetSingle(work_entry->global_noise);
	const tfxWideInt weight_last_frame = tfxWideSetSinglei(work_entry->graphs->weight.lookup.last_frame);
	const tfxWideInt velocity_turbulance_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity_turbulance.lookup.last_frame);
	const tfxWideInt noise_resolution_last_frame = tfxWideSetSinglei(work_entry->graphs->noise_resolution.lookup.last_frame);
	const tfxWideFloat node_count = tfxWideSetSingle(work_entry->node_count);

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxControlParticleLoopSetup;
		tfxControlParticleSampleOverPathLife;

		tfxWideInt velocity_normal = tfxWideLoadi((tfxWideIntLoader *)&bank.velocity_normal[index]);
		tfx__wide_unpack10bit(velocity_normal, velocity_normal_x, velocity_normal_y, velocity_normal_z);

		life = tfxWideMul(life, max_life);
		life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

		const tfxWideFloat base_velocity = tfxWideLoad(&bank.base_velocity[index]);
		const tfxWideFloat base_weight = tfxWideLoad(&bank.base_weight[index]);

		tfxWideArrayi lookup_frame;
		lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_last_frame);
		const tfxWideFloat lookup_velocity = tfxWideLookupSet(work_entry->graphs->velocity.lookup.values, lookup_frame);

		//----Velocity Changes
		tfxWideFloat velocity_scalar = tfxWideMul(base_velocity, lookup_velocity);

		tfxWideFloat local_position_x = tfxWideLoad(&bank.position_x[index]);
		tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
		tfxWideFloat local_position_z = tfxWideLoad(&bank.position_z[index]);

		tfxControlParticleUpdateNoise;

		tfxParticleNoise3dLoopUnroll(0)
			tfxParticleNoise3dLoopUnroll(1)
			tfxParticleNoise3dLoopUnroll(2)
			tfxParticleNoise3dLoopUnroll(3)
#if defined(tfxUSEAVX)    
			tfxParticleNoise3dLoopUnroll(4)
			tfxParticleNoise3dLoopUnroll(5)
			tfxParticleNoise3dLoopUnroll(6)
			tfxParticleNoise3dLoopUnroll(7)
#endif    

			noise_x.m = tfxWideMul(global_noise, tfxWideMul(lookup_velocity_turbulance, noise_x.m));
		noise_y.m = tfxWideMul(global_noise, tfxWideMul(lookup_velocity_turbulance, noise_y.m));
		noise_z.m = tfxWideMul(global_noise, tfxWideMul(lookup_velocity_turbulance, noise_z.m));

		current_velocity_x = tfxWideAdd(current_velocity_x, noise_x.m);
		current_velocity_y = tfxWideAdd(current_velocity_y, noise_y.m);
		current_velocity_z = tfxWideAdd(current_velocity_z, noise_z.m);

		tfxControlParticleUpdatePosition;
	}

}

//Used for emitters that have simplex noise only and orbital emission type non relative positioning
void ControlParticlePosition3dNoiseOrbital(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfxU32 emitter_index = work_entry->emitter_index;
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[emitter_index];
	tfx_particle_soa_t &bank = work_entry->pm->particle_arrays[emitter.particles_index];
	const tfxWideFloat overal_scale_wide = tfxWideSetSingle(work_entry->overal_scale);
	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
	const tfxWideInt velocity_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity.lookup.last_frame);
	const tfxWideFloat velocity_adjuster = tfxWideSetSingle(lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].overtime.velocity_adjuster, emitter.frame));
	const tfxWideInt weight_last_frame = tfxWideSetSinglei(work_entry->graphs->weight.lookup.last_frame);
	const tfxWideInt velocity_turbulance_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity_turbulance.lookup.last_frame);
	const tfxWideFloat global_noise = tfxWideSetSingle(work_entry->global_noise);
	const tfxWideInt noise_resolution_last_frame = tfxWideSetSinglei(work_entry->graphs->noise_resolution.lookup.last_frame);
	const tfxWideFloat node_count = tfxWideSetSingle(work_entry->node_count);
	tfxWideFloat emitter_x = {};
	tfxWideFloat emitter_z = {};
	if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
		emitter_x = tfxWideSetSingle(emitter.world_position.x);
		emitter_z = tfxWideSetSingle(emitter.world_position.z);
	}
	else if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position) {
		emitter_x = tfxWideSetSingle(emitter.handle.x);
		emitter_z = tfxWideSetSingle(emitter.handle.z);
	}

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxControlParticleLoopSetup;
		tfxControlParticleSampleOverPathLife;

		tfxWideFloat local_position_x = tfxWideLoad(&bank.position_x[index]);
		tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
		tfxWideFloat local_position_z = tfxWideLoad(&bank.position_z[index]);

		//Calculate orbital trajectory that is relative to the emitter
		tfxControlParticleOrbital;
		//---

		life = tfxWideMul(life, max_life);
		life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

		const tfxWideFloat base_velocity = tfxWideLoad(&bank.base_velocity[index]);
		const tfxWideFloat base_weight = tfxWideLoad(&bank.base_weight[index]);

		tfxWideArrayi lookup_frame;
		lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_last_frame);
		const tfxWideFloat lookup_velocity = tfxWideLookupSet(work_entry->graphs->velocity.lookup.values, lookup_frame);

		//----Velocity Changes
		tfxWideFloat velocity_scalar = tfxWideMul(base_velocity, lookup_velocity);

		tfxControlParticleUpdateNoise;

		tfxParticleNoise3dLoopUnroll(0)
			tfxParticleNoise3dLoopUnroll(1)
			tfxParticleNoise3dLoopUnroll(2)
			tfxParticleNoise3dLoopUnroll(3)
#if defined(tfxUSEAVX)    
			tfxParticleNoise3dLoopUnroll(4)
			tfxParticleNoise3dLoopUnroll(5)
			tfxParticleNoise3dLoopUnroll(6)
			tfxParticleNoise3dLoopUnroll(7)
#endif    

			noise_x.m = tfxWideMul(global_noise, tfxWideMul(lookup_velocity_turbulance, noise_x.m));
		noise_y.m = tfxWideMul(global_noise, tfxWideMul(lookup_velocity_turbulance, noise_y.m));
		noise_z.m = tfxWideMul(global_noise, tfxWideMul(lookup_velocity_turbulance, noise_z.m));

		current_velocity_x = tfxWideAdd(current_velocity_x, noise_x.m);
		current_velocity_y = tfxWideAdd(current_velocity_y, noise_y.m);
		current_velocity_z = tfxWideAdd(current_velocity_z, noise_z.m);

		tfxControlParticleUpdatePosition;
	}

}

//Used for emitters that have simple motion randomness only
void ControlParticlePosition3dMotionRandomness(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfxU32 emitter_index = work_entry->emitter_index;
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[emitter_index];
	tfx_particle_soa_t &bank = work_entry->pm->particle_arrays[emitter.particles_index];
	const tfxWideFloat overal_scale_wide = tfxWideSetSingle(work_entry->overal_scale);
	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
	const tfxWideInt velocity_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity.lookup.last_frame);
	const tfxWideFloat velocity_adjuster = tfxWideSetSingle(lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].overtime.velocity_adjuster, emitter.frame));
	const tfxWideFloat global_noise = tfxWideSetSingle(work_entry->global_noise);
	const tfxWideInt weight_last_frame = tfxWideSetSinglei(work_entry->graphs->weight.lookup.last_frame);
	const tfxWideInt motion_randomness_last_frame = tfxWideSetSinglei(work_entry->graphs->motion_randomness.lookup.last_frame);
	const tfxWideFloat motion_randomness_base = tfxWideSetSingle(lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].variation.motion_randomness, emitter.frame));
	const tfxWideFloat node_count = tfxWideSetSingle(work_entry->node_count);

	tfxWideInt time_step = tfxWideConverti(tfxWideSetSingle(emitter.age / 250.f));
	tfxWideInt next_time_step = tfxWideConverti(tfxWideSetSingle((emitter.age + pm.frame_length) / 250.f));
	tfxWideInt time_changed_mask = tfxWideLessi(time_step, next_time_step);
	tfxWideFloat time_step_fraction = tfxWideSub(tfxWideSetSingle(emitter.age / 250.f), tfxWideConvert(time_step));
	time_step_fraction = tfxWideMul(tfxWideMul(time_step_fraction, time_step_fraction), tfxWideSub(tfxWideSetSingle(3.f), tfxWideMul(tfxWideSetSingle(2.f), time_step_fraction)));

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxControlParticleLoopSetup;
		tfxControlParticleSampleOverPathLife;

		tfxWideInt velocity_normal = tfxWideLoadi((tfxWideIntLoader *)&bank.velocity_normal[index]);
		tfx__wide_unpack10bit(velocity_normal, velocity_normal_x, velocity_normal_y, velocity_normal_z);

		life = tfxWideMul(life, max_life);
		life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

		const tfxWideFloat base_velocity = tfxWideLoad(&bank.base_velocity[index]);
		const tfxWideFloat base_weight = tfxWideLoad(&bank.base_weight[index]);

		tfxWideArrayi lookup_frame;
		lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_last_frame);
		const tfxWideFloat lookup_velocity = tfxWideLookupSet(work_entry->graphs->velocity.lookup.values, lookup_frame);

		//----Velocity Changes
		tfxWideFloat velocity_scalar = tfxWideMul(base_velocity, lookup_velocity);

		tfxControlParticleMotionRandomness;

		//Non Orbit emission direction
		//Add the random direction to the current velocity
		//return current * tween + captured * (1.f - tween);
		velocity_normal_x = tfxWideAdd(tfxWideMul(random_x, time_step_fraction), tfxWideMul(velocity_normal_x, tfxWideSub(tfxWIDEONE, time_step_fraction)));
		velocity_normal_y = tfxWideAdd(tfxWideMul(random_y, time_step_fraction), tfxWideMul(velocity_normal_y, tfxWideSub(tfxWIDEONE, time_step_fraction)));
		velocity_normal_z = tfxWideAdd(tfxWideMul(random_z, time_step_fraction), tfxWideMul(velocity_normal_z, tfxWideSub(tfxWIDEONE, time_step_fraction)));
		length = tfxWideMul(velocity_normal_x, velocity_normal_x);
		length = tfxWideAdd(length, tfxWideMul(velocity_normal_y, velocity_normal_y));
		length = tfxWideAdd(length, tfxWideMul(velocity_normal_z, velocity_normal_z));
		length = tfxWideMul(tfxWideRSqrt(length), length);
		velocity_normal_x = tfxWideDiv(velocity_normal_x, length);
		velocity_normal_y = tfxWideDiv(velocity_normal_y, length);
		velocity_normal_z = tfxWideDiv(velocity_normal_z, length);
		//--

		tfxWideFloat current_velocity_x = tfxWideMul(velocity_normal_x, velocity_scalar);
		tfxWideFloat current_velocity_y = tfxWideMul(velocity_normal_y, velocity_scalar);
		tfxWideFloat current_velocity_z = tfxWideMul(velocity_normal_z, velocity_scalar);

		tfxWideInt packed_normal = tfx__wide_pack10bit_unsigned(velocity_normal_x, velocity_normal_y, velocity_normal_z);
		tfxWideInt normal_to_store = tfxWideOri(tfxWideAndi(packed_normal, time_changed_mask), tfxWideAndi(velocity_normal, tfxWideXOri(time_changed_mask, tfxWIDEMINUSONEi)));
		tfxWideStorei((tfxWideIntLoader *)&bank.velocity_normal[index], normal_to_store);
		tfxWideStore(&bank.noise_offset[index], speed);

		tfxWideFloat local_position_x = tfxWideLoad(&bank.position_x[index]);
		tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
		tfxWideFloat local_position_z = tfxWideLoad(&bank.position_z[index]);

		tfxControlParticleUpdatePosition;
	}
}

//Used for emitters that have simple motion randomness only with orbital emission type
void ControlParticlePosition3dMotionRandomnessOrbital(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfxU32 emitter_index = work_entry->emitter_index;
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[emitter_index];
	tfx_particle_soa_t &bank = work_entry->pm->particle_arrays[emitter.particles_index];
	const tfxWideFloat overal_scale_wide = tfxWideSetSingle(work_entry->overal_scale);
	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
	const tfxWideInt velocity_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity.lookup.last_frame);
	const tfxWideFloat velocity_adjuster = tfxWideSetSingle(lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].overtime.velocity_adjuster, emitter.frame));
	const tfxWideInt weight_last_frame = tfxWideSetSinglei(work_entry->graphs->weight.lookup.last_frame);
	const tfxWideInt motion_randomness_last_frame = tfxWideSetSinglei(work_entry->graphs->motion_randomness.lookup.last_frame);
	const tfxWideFloat motion_randomness_base = tfxWideSetSingle(lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].variation.motion_randomness, emitter.frame));
	const tfxWideFloat global_noise = tfxWideSetSingle(work_entry->global_noise);
	const tfxWideFloat node_count = tfxWideSetSingle(work_entry->node_count);

	tfxWideInt time_step = tfxWideConverti(tfxWideSetSingle(emitter.age / 250.f));
	tfxWideInt next_time_step = tfxWideConverti(tfxWideSetSingle((emitter.age + pm.frame_length) / 250.f));
	tfxWideInt time_changed_mask = tfxWideLessi(time_step, next_time_step);
	tfxWideFloat time_step_fraction = tfxWideSub(tfxWideSetSingle(emitter.age / 250.f), tfxWideConvert(time_step));
	time_step_fraction = tfxWideMul(tfxWideMul(time_step_fraction, time_step_fraction), tfxWideSub(tfxWideSetSingle(3.f), tfxWideMul(tfxWideSetSingle(2.f), time_step_fraction)));

	tfxWideFloat emitter_x = {};
	tfxWideFloat emitter_z = {};
	if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
		emitter_x = tfxWideSetSingle(emitter.world_position.x);
		emitter_z = tfxWideSetSingle(emitter.world_position.z);
	}
	else if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position) {
		emitter_x = tfxWideSetSingle(emitter.handle.x);
		emitter_z = tfxWideSetSingle(emitter.handle.z);
	}

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxControlParticleLoopSetup;
		tfxControlParticleSampleOverPathLife;

		tfxWideFloat local_position_x = tfxWideLoad(&bank.position_x[index]);
		tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
		tfxWideFloat local_position_z = tfxWideLoad(&bank.position_z[index]);

		//Calculate orbital trajectory that is relative to the emitter
		tfxControlParticleOrbital;
		//---

		life = tfxWideMul(life, max_life);
		life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

		const tfxWideFloat base_velocity = tfxWideLoad(&bank.base_velocity[index]);
		const tfxWideFloat base_weight = tfxWideLoad(&bank.base_weight[index]);

		tfxWideArrayi lookup_frame;
		lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_last_frame);
		const tfxWideFloat lookup_velocity = tfxWideLookupSet(work_entry->graphs->velocity.lookup.values, lookup_frame);

		//----Velocity Changes
		tfxWideFloat velocity_scalar = tfxWideMul(base_velocity, lookup_velocity);

		tfxControlParticleMotionRandomness;

		//Orbit emission direction for both relative and non relative particles
		tfxWideFloat vx, vy, vz;
		tfxWideInt velocity_normal = tfxWideLoadi((tfxWideIntLoader *)&bank.velocity_normal[index]);
		tfx__wide_unpack10bit(velocity_normal, vx, vy, vz);
		vx = tfxWideAdd(tfxWideMul(random_x, time_step_fraction), tfxWideMul(vx, tfxWideSub(tfxWIDEONE, time_step_fraction)));
		vy = tfxWideAdd(tfxWideMul(random_y, time_step_fraction), tfxWideMul(vy, tfxWideSub(tfxWIDEONE, time_step_fraction)));
		vz = tfxWideAdd(tfxWideMul(random_z, time_step_fraction), tfxWideMul(vz, tfxWideSub(tfxWIDEONE, time_step_fraction)));
		velocity_normal_x = tfxWideAdd(velocity_normal_x, vx);
		velocity_normal_y = tfxWideAdd(velocity_normal_y, vy);
		velocity_normal_z = tfxWideAdd(velocity_normal_z, vz);
		length = tfxWideMul(velocity_normal_x, velocity_normal_x);
		length = tfxWideAdd(length, tfxWideMul(velocity_normal_y, velocity_normal_y));
		length = tfxWideAdd(length, tfxWideMul(velocity_normal_z, velocity_normal_z));
		length = tfxWideMul(tfxWideRSqrt(length), length);
		velocity_normal_x = tfxWideDiv(velocity_normal_x, length);
		velocity_normal_y = tfxWideDiv(velocity_normal_y, length);
		velocity_normal_z = tfxWideDiv(velocity_normal_z, length);
		length = tfxWideMul(vx, vx);
		length = tfxWideAdd(length, tfxWideMul(vy, vy));
		length = tfxWideAdd(length, tfxWideMul(vz, vz));
		length = tfxWideMul(tfxWideRSqrt(length), length);
		vx = tfxWideDiv(vx, length);
		vy = tfxWideDiv(vy, length);
		vz = tfxWideDiv(vz, length);
		tfxWideInt packed_normal = tfx__wide_pack10bit_unsigned(vx, vy, vz);
		//--

		tfxWideFloat current_velocity_x = tfxWideMul(velocity_normal_x, velocity_scalar);
		tfxWideFloat current_velocity_y = tfxWideMul(velocity_normal_y, velocity_scalar);
		tfxWideFloat current_velocity_z = tfxWideMul(velocity_normal_z, velocity_scalar);

		tfxWideInt normal_to_store = tfxWideOri(tfxWideAndi(packed_normal, time_changed_mask), tfxWideAndi(velocity_normal, tfxWideXOri(time_changed_mask, tfxWIDEMINUSONEi)));
		tfxWideStorei((tfxWideIntLoader *)&bank.velocity_normal[index], normal_to_store);
		tfxWideStore(&bank.noise_offset[index], speed);

		tfxControlParticleUpdatePosition;
	}

}

void ControlParticleLineBehaviourKill(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfxU32 emitter_index = work_entry->emitter_index;
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[emitter_index];
	tfx_particle_soa_t &bank = work_entry->pm->particle_arrays[emitter.particles_index];
	const tfxWideFloat emitter_size_y = tfxWideSetSingle(emitter.emitter_size.y);

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;
		tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
		tfxWideInt flags = tfxWideLoadi((tfxWideIntLoader *)&bank.flags[index]);

		tfx__readbarrier;

		//Lines - Kill if the particle has reached the end of the line
		tfxWideInt remove_flags = tfxWideAndi(tfxWideSetSinglei(tfxParticleFlags_remove), tfxWideCasti(tfxWideGreater(local_position_y, emitter_size_y)));
		flags = tfxWideOri(flags, remove_flags);
		remove_flags = tfxWideAndi(tfxWideSetSinglei(tfxParticleFlags_remove), tfxWideCasti(tfxWideLess(local_position_y, tfxWideSetZero)));
		flags = tfxWideOri(flags, remove_flags);
		tfxWideStorei((tfxWideIntLoader *)&bank.flags[index], flags);
	}
}

void ControlParticleLineBehaviourLoop(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfxU32 emitter_index = work_entry->emitter_index;
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[emitter_index];
	tfx_particle_soa_t &bank = work_entry->pm->particle_arrays[emitter.particles_index];
	const tfxWideFloat emitter_size_y = tfxWideSetSingle(emitter.emitter_size.y);

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;
		tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
		tfxWideInt flags = tfxWideLoadi((tfxWideIntLoader *)&bank.flags[index]);

		tfx__readbarrier;

		//Lines - Reposition if the particle is travelling along a line
		tfxWideFloat at_end = tfxWideGreater(local_position_y, emitter_size_y);
		local_position_y = tfxWideSub(local_position_y, tfxWideAnd(at_end, emitter_size_y));
		flags = tfxWideOri(flags, tfxWideAndi(tfxWideSetSinglei(tfxParticleFlags_capture_after_transform), tfxWideCasti(at_end)));
		at_end = tfxWideLess(local_position_y, tfxWideSetZero);
		local_position_y = tfxWideAdd(local_position_y, tfxWideAnd(at_end, emitter_size_y));
		flags = tfxWideOri(flags, tfxWideAndi(tfxWideSetSinglei(tfxParticleFlags_capture_after_transform), tfxWideCasti(at_end)));
		tfxWideStorei((tfxWideIntLoader *)&bank.flags[index], flags);
		tfxWideStore(&bank.position_y[index], local_position_y);
	}
}

void ControlParticleTransform3d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[work_entry->emitter_index];
	tfx_bounding_box_t &bounding_box = emitter.bounding_box;
	tfx_particle_soa_t &bank = work_entry->pm->particle_arrays[emitter.particles_index];

	tfxU32 running_sprite_index = work_entry->sprites_index;

	const tfxWideFloat e_world_position_x = tfxWideSetSingle(emitter.world_position.x);
	const tfxWideFloat e_world_position_y = tfxWideSetSingle(emitter.world_position.y);
	const tfxWideFloat e_world_position_z = tfxWideSetSingle(emitter.world_position.z);
	const tfxWideFloat e_handle_x = tfxWideSetSingle(emitter.handle.x);
	const tfxWideFloat e_handle_y = tfxWideSetSingle(emitter.handle.y);
	const tfxWideFloat e_handle_z = tfxWideSetSingle(emitter.handle.z);
	const tfxWideFloat e_scale = tfxWideSetSingle(work_entry->overal_scale);
	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
	const tfxWideInt stretch_last_frame = tfxWideSetSinglei(work_entry->graphs->stretch.lookup.last_frame);
	const tfxWideFloat stretch = tfxWideSetSingle(work_entry->global_stretch);
	tfxWideArray p_stretch;
	tfxU32 start_diff = work_entry->start_diff;

	const tfxWideInt capture_after_transform = tfxWideSetSinglei(tfxParticleFlags_capture_after_transform);
	const tfxEmitterPropertyFlags property_flags = emitter.property_flags;
	const tfx_vector_align_type vector_align_type = work_entry->properties->vector_align_type;
	const tfx_emission_type emission_type = work_entry->properties->emission_type;
	const tfxU32 sprite_layer = work_entry->properties->layer;
	tfx_3d_instance_t *sprites = tfxCastBuffer(tfx_3d_instance_t, work_entry->sprite_instances);

	bool is_ordered = (!(pm.flags & tfxParticleManagerFlags_unordered) || (IsOrderedEffectState(&pm.effects[emitter.root_index])));

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;
		const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
		const tfxWideFloat age = tfxWideLoad(&bank.age[index]);
		tfx__readbarrier;
		tfxWideFloat life = tfxWideDiv(age, max_age);
		life = tfxWideMul(life, max_life);
		life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

		tfxWideArrayi lookup_frame;
		lookup_frame.m = tfxWideMini(tfxWideConverti(life), stretch_last_frame);
		const tfxWideFloat lookup_stretch = tfxWideLookupSet(work_entry->graphs->stretch.lookup.values, lookup_frame);
		p_stretch.m = tfxWideMul(lookup_stretch, stretch);

		tfxWideArray position_x;
		tfxWideArray position_y;
		tfxWideArray position_z;
		position_x.m = tfxWideLoad(&bank.position_x[index]);
		position_y.m = tfxWideLoad(&bank.position_y[index]);
		position_z.m = tfxWideLoad(&bank.position_z[index]);
		tfxWideArray captured_position_x;
		tfxWideArray captured_position_y;
		tfxWideArray captured_position_z;
		captured_position_x.m = tfxWideLoad(&bank.captured_position_x[index]);
		captured_position_y.m = tfxWideLoad(&bank.captured_position_y[index]);
		captured_position_z.m = tfxWideLoad(&bank.captured_position_z[index]);
		tfxWideInt flags = tfxWideLoadi((tfxWideIntLoader *)&bank.flags[index]);
		tfxWideArray capture_flag;
		tfx__readbarrier;
		capture_flag.m = tfxWideCast(tfxWideGreateri(tfxWideAndi(flags, capture_after_transform), tfxWideSetZeroi));
		tfxWideFloat xor_capture_flag = tfxWideEquals(capture_flag.m, tfxWideSetZero);
		tfx__readbarrier;

		tfxWideFloat alignment_vector_x;
		tfxWideFloat alignment_vector_y;
		tfxWideFloat alignment_vector_z;

		if (property_flags & tfxEmitterPropertyFlags_relative_position && emission_type != tfxPath) {
			position_x.m = tfxWideAdd(position_x.m, e_handle_x);
			position_y.m = tfxWideAdd(position_y.m, e_handle_y);
			position_z.m = tfxWideAdd(position_z.m, e_handle_z);
			tfx__wide_transform_quaternion_vec3(&emitter.rotation, &position_x.m, &position_y.m, &position_z.m);
			position_x.m = tfxWideAdd(tfxWideMul(position_x.m, e_scale), e_world_position_x);
			position_y.m = tfxWideAdd(tfxWideMul(position_y.m, e_scale), e_world_position_y);
			position_z.m = tfxWideAdd(tfxWideMul(position_z.m, e_scale), e_world_position_z);
		}
		else if (property_flags & tfxEmitterPropertyFlags_relative_position && emission_type == tfxPath) {
			tfx__wide_transform_quaternion_vec3(&emitter.rotation, &position_x.m, &position_y.m, &position_z.m);
			position_x.m = tfxWideAdd(tfxWideMul(position_x.m, e_scale), e_world_position_x);
			position_y.m = tfxWideAdd(tfxWideMul(position_y.m, e_scale), e_world_position_y);
			position_z.m = tfxWideAdd(tfxWideMul(position_z.m, e_scale), e_world_position_z);
		}

		captured_position_x.m = tfxWideAdd(tfxWideAnd(position_x.m, capture_flag.m), tfxWideAnd(captured_position_x.m, xor_capture_flag));
		captured_position_y.m = tfxWideAdd(tfxWideAnd(position_y.m, capture_flag.m), tfxWideAnd(captured_position_y.m, xor_capture_flag));
		captured_position_z.m = tfxWideAdd(tfxWideAnd(position_z.m, capture_flag.m), tfxWideAnd(captured_position_z.m, xor_capture_flag));

		if (vector_align_type == tfxVectorAlignType_motion) {
			alignment_vector_x = tfxWideSub(position_x.m, captured_position_x.m);
			alignment_vector_y = tfxWideSub(position_y.m, captured_position_y.m);
			alignment_vector_z = tfxWideAdd(tfxWideSub(position_z.m, captured_position_z.m), tfxWideSetSingle(0.000001f)); //epsilon to prevent divide by 0
			tfxWideFloat l = tfxWideMul(alignment_vector_x, alignment_vector_x);
			l = tfxWideAdd(l, tfxWideMul(alignment_vector_y, alignment_vector_y));
			l = tfxWideAdd(l, tfxWideMul(alignment_vector_z, alignment_vector_z));
			l = tfxWideMul(tfxWideRSqrt(l), l);

			//We divide the length by the overal scale because we don't want the stretch to be effected by it.
			p_stretch.m = tfxWideMul(p_stretch.m, tfxWideDiv(tfxWideDiv(l, e_scale), pm.update_time_wide));    //This is too arbitrary, think up a better solution!
			alignment_vector_x = tfxWideDiv(alignment_vector_x, l);
			alignment_vector_y = tfxWideDiv(alignment_vector_y, l);
			alignment_vector_z = tfxWideDiv(alignment_vector_z, l);
		}
		else if (vector_align_type == tfxVectorAlignType_emission && property_flags & tfxEmitterPropertyFlags_relative_position) {
			const tfxWideInt velocity_normal = tfxWideLoadi((tfxWideIntLoader *)&bank.velocity_normal[index]);
			tfxWideFloat velocity_normal_x;
			tfxWideFloat velocity_normal_y;
			tfxWideFloat velocity_normal_z;
			tfx__wide_unpack10bit(velocity_normal, velocity_normal_x, velocity_normal_y, velocity_normal_z);
			alignment_vector_x = velocity_normal_x;
			alignment_vector_y = velocity_normal_y;
			alignment_vector_z = velocity_normal_z;
			tfx__wide_transform_quaternion_vec3(&emitter.rotation, &alignment_vector_x, &alignment_vector_y, &alignment_vector_z);
		}
		else if (vector_align_type == tfxVectorAlignType_emission) {
			const tfxWideInt velocity_normal = tfxWideLoadi((tfxWideIntLoader *)&bank.velocity_normal[index]);
			tfxWideFloat velocity_normal_x;
			tfxWideFloat velocity_normal_y;
			tfxWideFloat velocity_normal_z;
			tfx__wide_unpack10bit(velocity_normal, velocity_normal_x, velocity_normal_y, velocity_normal_z);
			alignment_vector_x = velocity_normal_x;
			alignment_vector_y = velocity_normal_y;
			alignment_vector_z = velocity_normal_z;
		}
		else if (vector_align_type == tfxVectorAlignType_emitter) {
			alignment_vector_x = tfxWideSetZero;
			alignment_vector_y = tfxWideSetSingle(1.f);
			alignment_vector_z = tfxWideSetZero;
			tfx__wide_transform_quaternion_vec3(&emitter.rotation, &alignment_vector_x, &alignment_vector_y, &alignment_vector_z);
		}

		//instance_data.transform_3d.captured_position = captured_position;
		//alignment_vector_y.m = tfxWideAdd(alignment_vector_y.m, tfxWideSetSingle(0.002f));    //We don't want a 0 alignment normal
		tfxWideArrayi alignment_packed;
		alignment_packed.m = tfx__wide_pack8bit_xyz(alignment_vector_x, alignment_vector_y, alignment_vector_z);

		if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_location_source && emitter.spawn_locations_index != tfxINVALID) {
			tfx_spawn_points_soa_t &locations = work_entry->pm->particle_location_arrays[emitter.spawn_locations_index];
			tfxWideStore(&locations.position_x[index], position_x.m);
			tfxWideStore(&locations.position_y[index], position_y.m);
			tfxWideStore(&locations.position_z[index], position_z.m);
			tfxWideStore(&locations.captured_position_x[index], captured_position_x.m);
			tfxWideStore(&locations.captured_position_y[index], captured_position_y.m);
			tfxWideStore(&locations.captured_position_z[index], captured_position_z.m);
			tfxWideStore(&locations.age[index], tfxWideDiv(age, max_age));
		}

		tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
		if (is_ordered) {
			for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
				int index_j = index + j;
				tfxU32 sprite_depth_index = bank.depth_index[index_j] + work_entry->cumulative_index_point + work_entry->effect_instance_offset;
				sprites[sprite_depth_index].alignment.packed = alignment_packed.a[j];
				sprites[sprite_depth_index].position.w = p_stretch.a[j];
				sprites[sprite_depth_index].position.x = position_x.a[j];
				sprites[sprite_depth_index].position.y = position_y.a[j];
				sprites[sprite_depth_index].position.z = position_z.a[j];
				bank.captured_position_x[index_j] = sprites[sprite_depth_index].position.x;
				bank.captured_position_y[index_j] = sprites[sprite_depth_index].position.y;
				bank.captured_position_z[index_j] = sprites[sprite_depth_index].position.z;
				tfx_vec3_t sprite_plus_camera_position = sprites[sprite_depth_index].position.xyz() - pm.camera_position;
				(*work_entry->depth_indexes)[sprite_depth_index - work_entry->cumulative_index_point - work_entry->effect_instance_offset].depth = tfx__length_vec3_nosqr(&sprite_plus_camera_position);
				if (pm.flags & tfxParticleManagerFlags_update_bounding_boxes) {
					bounding_box.min_corner.x = tfx__Min(position_x.a[j], bounding_box.min_corner.x);
					bounding_box.min_corner.y = tfx__Min(position_y.a[j], bounding_box.min_corner.y);
					bounding_box.min_corner.z = tfx__Min(position_z.a[j], bounding_box.min_corner.z);
					bounding_box.max_corner.x = tfx__Max(position_x.a[j], bounding_box.max_corner.x);
					bounding_box.max_corner.y = tfx__Max(position_y.a[j], bounding_box.max_corner.y);
					bounding_box.max_corner.z = tfx__Max(position_z.a[j], bounding_box.max_corner.z);
				}
				running_sprite_index++;
			}
		}
		else {
			for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
				int index_j = index + j;
				sprites[running_sprite_index].position.w = p_stretch.a[j];
				sprites[running_sprite_index].alignment.packed = alignment_packed.a[j];
				sprites[running_sprite_index].position.x = position_x.a[j];
				sprites[running_sprite_index].position.y = position_y.a[j];
				sprites[running_sprite_index].position.z = position_z.a[j];
				bank.captured_position_x[index_j] = sprites[running_sprite_index].position.x;
				bank.captured_position_y[index_j] = sprites[running_sprite_index].position.y;
				bank.captured_position_z[index_j] = sprites[running_sprite_index].position.z;
				if (pm.flags & tfxParticleManagerFlags_update_bounding_boxes) {
					bounding_box.min_corner.x = tfx__Min(position_x.a[j], bounding_box.min_corner.x);
					bounding_box.min_corner.y = tfx__Min(position_y.a[j], bounding_box.min_corner.y);
					bounding_box.min_corner.z = tfx__Min(position_z.a[j], bounding_box.min_corner.z);
					bounding_box.max_corner.x = tfx__Max(position_x.a[j], bounding_box.max_corner.x);
					bounding_box.max_corner.y = tfx__Max(position_y.a[j], bounding_box.max_corner.y);
					bounding_box.max_corner.z = tfx__Max(position_z.a[j], bounding_box.max_corner.z);
				}
				running_sprite_index++;
			}
		}

		start_diff = 0;
	}
}

void ControlParticlePosition2d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[work_entry->emitter_index];
	tfx_particle_soa_t &bank = work_entry->pm->particle_arrays[emitter.particles_index];

	const tfxWideFloat overal_scale_wide = tfxWideSetSingle(work_entry->overal_scale);

	tfxU32 running_sprite_index = work_entry->sprites_index;
	tfx_2d_instance_t *sprites = tfxCastBuffer(tfx_2d_instance_t, work_entry->sprite_instances);

	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
	const tfxWideInt velocity_turbulance_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity_turbulance.lookup.last_frame);
	const tfxWideInt noise_resolution_last_frame = tfxWideSetSinglei(work_entry->graphs->noise_resolution.lookup.last_frame);
	tfxWideInt time_step = tfxWideConverti(tfxWideSetSingle(emitter.age / 100.f));

	//Noise
	const tfxWideInt velocity_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity.lookup.last_frame);
	const tfxWideInt spin_last_frame = tfxWideSetSinglei(work_entry->graphs->spin.lookup.last_frame);
	const tfxWideFloat velocity_adjuster = tfxWideSetSingle(lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].overtime.velocity_adjuster, emitter.frame));
	const tfxWideInt weight_last_frame = tfxWideSetSinglei(work_entry->graphs->weight.lookup.last_frame);
	const tfxWideInt direction_last_frame = tfxWideSetSinglei(work_entry->graphs->direction.lookup.last_frame);
	const tfxWideInt stretch_last_frame = tfxWideSetSinglei(work_entry->graphs->stretch.lookup.last_frame);
	const tfxWideFloat stretch = tfxWideSetSingle(work_entry->global_stretch);
	const tfxWideInt motion_randomness_last_frame = tfxWideSetSinglei(work_entry->graphs->motion_randomness.lookup.last_frame);
	const tfxWideFloat motion_randomness_base = tfxWideSetSingle(lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].variation.motion_randomness, emitter.frame));
	tfxWideFloat emitter_x;
	tfxWideFloat emitter_y;
	tfxWideArray p_stretch;

	tfxU32 start_diff = work_entry->start_diff;

	const float eps = 0.0001f;
	const float eps2 = 0.0002f;

	bool orbit_non_relative = false;
	bool orbit_relative = false;
	if (work_entry->properties->emission_direction == tfxOrbital && !(emitter.property_flags & tfxEmitterPropertyFlags_relative_position) && work_entry->properties->emission_type != tfxPoint) {
		emitter_x = tfxWideSetSingle(emitter.world_position.x);
		emitter_y = tfxWideSetSingle(emitter.world_position.y);
		orbit_non_relative = true;
	}
	else if (work_entry->properties->emission_direction == tfxOrbital && (emitter.property_flags & tfxEmitterPropertyFlags_relative_position) && work_entry->properties->emission_type != tfxPoint) {
		emitter_x = tfxWideSetSingle(emitter.handle.x);
		emitter_y = tfxWideSetSingle(emitter.handle.y);
		orbit_relative = true;
	}

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;

		const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
		const tfxWideFloat age = tfxWideLoad(&bank.age[index]);

		tfx__readbarrier;

		tfxWideFloat life = tfxWideDiv(age, max_age);
		life = tfxWideMul(life, max_life);
		life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

		tfxWideFloat velocity_normal_x = tfxWideSetZero;
		tfxWideFloat velocity_normal_y = tfxWideSetZero;
		tfxWideFloat local_position_x = tfxWideLoad(&bank.position_x[index]);
		tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
		tfxWideArrayi lookup_frame;

		if (orbit_non_relative) {
			velocity_normal_y = tfxWideSub(local_position_x, emitter_x);
			velocity_normal_x = tfxWideMul(tfxWideSub(local_position_y, emitter_y), tfxWideSetSingle(-1.f));
			tfxWideFloat l = tfxWideMul(velocity_normal_x, velocity_normal_x);
			l = tfxWideAdd(l, tfxWideMul(velocity_normal_y, velocity_normal_y));
#ifdef tfxARM
			l = tfxWideMul(tfxWideRSqrt(l), l);
#else
			l = tfxWideSqrt(l);
#endif
			velocity_normal_x = tfxWideDiv(velocity_normal_x, l);
			velocity_normal_y = tfxWideDiv(velocity_normal_y, l);
		}
		else if (orbit_relative) {
			velocity_normal_y = tfxWideAdd(emitter_x, local_position_x);
			velocity_normal_x = tfxWideMul(tfxWideAdd(emitter_y, local_position_y), tfxWideSetSingle(-1.f));
			tfxWideFloat l = tfxWideMul(velocity_normal_x, velocity_normal_x);
			l = tfxWideAdd(l, tfxWideMul(velocity_normal_y, velocity_normal_y));
#ifdef tfxARM
			l = tfxWideMul(tfxWideRSqrt(l), l);
#else
			l = tfxWideSqrt(l);
#endif
			velocity_normal_x = tfxWideDiv(velocity_normal_x, l);
			velocity_normal_y = tfxWideDiv(velocity_normal_y, l);
		}
		else if (!(emitter.control_profile & tfxEmitterControlProfile_motion_randomness)) {
			tfxWideFloat angle = tfxWideLoad(&bank.local_rotations_x[index]);
			tfxWideArray lookup_direction;
			lookup_frame.m = tfxWideMini(tfxWideConverti(life), direction_last_frame);
			lookup_direction.m = tfxWideLookupSet(work_entry->graphs->direction.lookup.values, lookup_frame);
			lookup_direction.m = tfxWideAdd(lookup_direction.m, angle);
			tfxWideSinCos(lookup_direction.m, &velocity_normal_x, &velocity_normal_y);
			velocity_normal_y = tfxWideMul(velocity_normal_y, tfxWideSetSingle(-1.f));
		}
		else {
			velocity_normal_x = tfxWideSetZero;
			velocity_normal_y = tfxWideSetZero;
		}

		const tfxWideFloat base_velocity = tfxWideLoad(&bank.base_velocity[index]);
		const tfxWideFloat base_spin = tfxWideLoad(&bank.base_spin[index]);
		const tfxWideFloat base_weight = tfxWideLoad(&bank.base_weight[index]);
		const tfxWideFloat base_size_y = tfxWideLoad(&bank.base_size_y[index]);

		//----Velocity Changes
		lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_last_frame);
		const tfxWideFloat lookup_velocity = tfxWideLookupSet(work_entry->graphs->velocity.lookup.values, lookup_frame);
		tfxWideFloat velocity_scalar = tfxWideMul(base_velocity, lookup_velocity);

		if (emitter.control_profile & tfxEmitterControlProfile_noise) {
			velocity_normal_x = tfxWideMul(velocity_normal_x, velocity_scalar);
			velocity_normal_y = tfxWideMul(velocity_normal_y, velocity_scalar);

			tfxWideArray noise_x;
			tfxWideArray noise_y;

			const tfxWideFloat noise_resolution = tfxWideLoad(&bank.noise_resolution[index]);
			const tfxWideFloat noise_offset = tfxWideLoad(&bank.noise_offset[index]);

			tfx__readbarrier;

			lookup_frame.m = tfxWideMini(tfxWideConverti(life), noise_resolution_last_frame);
			const tfxWideFloat lookup_noise_resolution = tfxWideMul(tfxWideLookupSet(work_entry->graphs->noise_resolution.lookup.values, lookup_frame), noise_resolution);
			lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_turbulance_last_frame);
			const tfxWideFloat lookup_velocity_turbulance = tfxWideLookupSet(work_entry->graphs->velocity_turbulance.lookup.values, lookup_frame);

			tfxWideArray x, y;
			x.m = tfxWideAdd(tfxWideDiv(local_position_x, lookup_noise_resolution), noise_offset);
			y.m = tfxWideAdd(tfxWideDiv(local_position_y, lookup_noise_resolution), noise_offset);

			tfx128 x4, y4, xeps4, yeps4r;
			float a, b;
			tfx128Array sample;

			tfxParticleNoise2dLoopUnroll(0);
			tfxParticleNoise2dLoopUnroll(1);
			tfxParticleNoise2dLoopUnroll(2);
			tfxParticleNoise2dLoopUnroll(3);

#if defined(tfxUSEAVX)
			tfxParticleNoise2dLoopUnroll(4);
			tfxParticleNoise2dLoopUnroll(5);
			tfxParticleNoise2dLoopUnroll(6);
			tfxParticleNoise2dLoopUnroll(7);
#endif

			noise_x.m = tfxWideMul(lookup_velocity_turbulance, noise_x.m);
			noise_y.m = tfxWideMul(lookup_velocity_turbulance, noise_y.m);

			velocity_normal_x = tfxWideAdd(velocity_normal_x, noise_x.m);
			velocity_normal_y = tfxWideAdd(velocity_normal_y, noise_y.m);
		}
		else if (emitter.control_profile & tfxEmitterControlProfile_motion_randomness) {
			tfxWideInt uid = tfxWideLoadi((tfxWideIntLoader *)&bank.uid[index]);
			tfxWideInt seed = SeedGenWide(time_step, uid);
			tfxWideFloat speed = tfxWideLoad(&bank.noise_offset[index]);
			tfxWideFloat direction = tfxWideLoad(&bank.noise_resolution[index]);
			tfxWideFloat two = tfxWideSetSingle(2.f);
			tfxWideFloat max_uint = tfxWideSetSingle((float)UINT32_MAX);

			tfxWideArrayi lookup_motion_randomness;
			lookup_motion_randomness.m = tfxWideMini(tfxWideConverti(life), motion_randomness_last_frame);
			const tfxWideFloat influence = tfxWideMul(motion_randomness_base, tfxWideLookupSet(work_entry->graphs->motion_randomness.lookup.values, lookup_motion_randomness));

			tfxWideFloat speed_influence = tfxWideMul(tfxWideSetSingle(10.f), influence);
			tfxWideFloat max_degrees = tfxWideMul(tfxWideSetSingle(0.3926991f), influence); //22.5 degrees

			tfxWideFloat random_direction = tfxWideMul(tfxWideDiv(SeedGenWide(seed), max_uint), max_degrees);
			tfxWideFloat random_speed = tfxWideMul(tfxWideDiv(SeedGenWide(tfxWideShiftLeft(seed, 1)), max_uint), speed_influence);

			tfxWideFloat angle = tfxWideLoad(&bank.local_rotations_x[index]);
			tfxWideArray lookup_direction;
			lookup_frame.m = tfxWideMini(tfxWideConverti(life), direction_last_frame);
			lookup_direction.m = tfxWideLookupSet(work_entry->graphs->direction.lookup.values, lookup_frame);
			lookup_direction.m = tfxWideAdd(lookup_direction.m, angle);
			direction = tfxWideAdd(direction, random_direction);
			if (orbit_relative + orbit_non_relative) {
				tfxWideFloat multiplier = tfxWideSetSingle(2.f);
				velocity_normal_x = tfxWideMul(velocity_normal_x, multiplier);
				velocity_normal_y = tfxWideMul(velocity_normal_y, multiplier);
				tfxWideFloat mx, my;
				tfxWideSinCos(tfxWideAdd(direction, lookup_direction.m), &mx, &my);
				velocity_normal_x = tfxWideAdd(velocity_normal_x, tfxWideMul(mx, influence));
				velocity_normal_y = tfxWideAdd(velocity_normal_y, tfxWideMul(my, influence));
				tfxWideFloat l = tfxWideMul(velocity_normal_x, velocity_normal_x);
				l = tfxWideAdd(l, tfxWideMul(velocity_normal_y, velocity_normal_y));
#ifdef tfxARM
				l = tfxWideMul(tfxWideRSqrt(l), l);
#else
				l = tfxWideSqrt(l);
#endif
				velocity_normal_x = tfxWideDiv(velocity_normal_x, l);
				velocity_normal_y = tfxWideDiv(velocity_normal_y, l);
			}
			else {
				tfxWideSinCos(tfxWideAdd(direction, lookup_direction.m), &velocity_normal_x, &velocity_normal_y);
				velocity_normal_y = tfxWideMul(velocity_normal_y, tfxWideSetSingle(-1.f));
			}

			speed = tfxWideAdd(random_speed, speed);
			velocity_normal_x = tfxWideMul(velocity_normal_x, tfxWideAdd(speed, velocity_scalar));
			velocity_normal_y = tfxWideMul(velocity_normal_y, tfxWideAdd(speed, velocity_scalar));
			tfxWideStore(&bank.noise_offset[index], speed);
			tfxWideStore(&bank.noise_resolution[index], direction);
		}
		else {
			velocity_normal_x = tfxWideMul(velocity_normal_x, velocity_scalar);
			velocity_normal_y = tfxWideMul(velocity_normal_y, velocity_scalar);
		}

		lookup_frame.m = tfxWideMini(tfxWideConverti(life), spin_last_frame);
		const tfxWideFloat lookup_spin = tfxWideMul(tfxWideLookupSet(work_entry->graphs->spin.lookup.values, lookup_frame), base_spin);
		lookup_frame.m = tfxWideMini(tfxWideConverti(life), weight_last_frame);
		const tfxWideFloat lookup_weight = tfxWideLookupSet(work_entry->graphs->weight.lookup.values, lookup_frame);
		lookup_frame.m = tfxWideMini(tfxWideConverti(life), stretch_last_frame);
		const tfxWideFloat lookup_stretch = tfxWideLookupSet(work_entry->graphs->stretch.lookup.values, lookup_frame);
		p_stretch.m = tfxWideMul(lookup_stretch, stretch);

		tfxWideFloat stretch_velocity_x;
		tfxWideFloat stretch_velocity_y;

		velocity_normal_y = tfxWideAdd(velocity_normal_y, tfxWideMul(lookup_weight, base_weight));
		stretch_velocity_x = velocity_normal_x;
		stretch_velocity_y = velocity_normal_y;
		tfxWideFloat age_fraction = tfxWideMin(tfxWideDiv(age, pm.frame_length_wide), tfxWIDEONE);
		velocity_normal_x = tfxWideMul(tfxWideMul(tfxWideMul(velocity_normal_x, pm.update_time_wide), velocity_adjuster), age_fraction);
		velocity_normal_y = tfxWideMul(tfxWideMul(tfxWideMul(velocity_normal_y, pm.update_time_wide), velocity_adjuster), age_fraction);

		//----Position
		local_position_x = tfxWideAdd(local_position_x, tfxWideMul(velocity_normal_x, overal_scale_wide));
		local_position_y = tfxWideAdd(local_position_y, tfxWideMul(velocity_normal_y, overal_scale_wide));

		tfxWideStore(&bank.position_x[index], local_position_x);
		tfxWideStore(&bank.position_y[index], local_position_y);

		stretch_velocity_y = tfxWideAdd(stretch_velocity_y, tfxWideSetSingle(0.000001f));
		tfxWideFloat l = tfxWideMul(stretch_velocity_x, stretch_velocity_x);
		l = tfxWideAdd(l, tfxWideMul(stretch_velocity_y, stretch_velocity_y));
		l = tfxWideMul(tfxWideRSqrt(l), l);

		//We divide the length by the overal scale because we don't want the stretch to be effected by it.
		p_stretch.m = tfxWideMul(p_stretch.m, tfxWideMul(tfxWideDiv(l, overal_scale_wide), tfxWideSetSingle(0.001f)));
		stretch_velocity_x = tfxWideDiv(stretch_velocity_x, l);
		stretch_velocity_y = tfxWideDiv(stretch_velocity_y, l);

		if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position) {
			tfx__wide_transform_quaternion_vec2(&emitter.rotation, &stretch_velocity_x, &stretch_velocity_y);
		}

		tfxWideArrayi packed;
		packed.m = tfx__wide_pack16bit(tfxWideMax(tfxWIDEMINUSONE, tfxWideMin(tfxWIDEONE, stretch_velocity_x)), tfxWideMax(tfxWIDEMINUSONE, tfxWideMin(tfxWIDEONE, stretch_velocity_y)));

		tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
		if (!(pm.flags & tfxParticleManagerFlags_unordered)) {    //Predictable
			for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
				tfxU32 sprite_depth_index = bank.depth_index[index + j] + work_entry->cumulative_index_point + work_entry->effect_instance_offset;
				TFX_ASSERT(sprite_depth_index < work_entry->sprite_instances->current_size);
				sprites[sprite_depth_index].position.z = p_stretch.a[j];
				sprites[sprite_depth_index].alignment.packed = packed.a[j];
				running_sprite_index++;
			}
		}
		else {
			for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
				TFX_ASSERT(running_sprite_index < work_entry->sprite_instances->current_size);
				sprites[running_sprite_index].position.z = p_stretch.a[j];
				sprites[running_sprite_index].alignment.packed = packed.a[j];
				running_sprite_index++;
			}
		}
		start_diff = 0;
	}

	const tfxWideFloat emitter_size_y = tfxWideSetSingle(emitter.emitter_size.y);
	const tfxWideInt emitter_flags_wide = tfxWideSetSinglei(emitter.state_flags);

	if (emitter.state_flags & tfxEmitterStateFlags_is_line_loop_or_kill) {
		if (emitter.state_flags & tfxEmitterStateFlags_kill) {
			for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
				tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;
				const tfxWideFloat offset_y = tfxWideMul(tfx__wide_unpack10bit_y(tfxWideLoadi((tfxWideIntLoader *)&bank.velocity_normal[index])), emitter_size_y);
				tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
				tfxWideInt flags = tfxWideLoadi((tfxWideIntLoader *)&bank.flags[index]);

				tfx__readbarrier;

				//Lines - Reposition if the particle is travelling along a line
				tfxWideFloat length = tfxWideAbs(local_position_y);
				tfxWideInt remove_flags = tfxWideAndi(tfxWideSetSinglei(tfxParticleFlags_remove), tfxWideCasti(tfxWideGreater(length, emitter_size_y)));
				flags = tfxWideOri(flags, remove_flags);
				tfxWideStorei((tfxWideIntLoader *)&bank.flags[index], flags);
			}
		}
		else {
			for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
				tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;
				const tfxWideFloat offset_y = tfxWideMul(tfx__wide_unpack10bit_y(tfxWideLoadi((tfxWideIntLoader *)&bank.velocity_normal[index])), emitter_size_y);
				tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
				tfxWideInt flags = tfxWideLoadi((tfxWideIntLoader *)&bank.flags[index]);

				//Lines - Reposition if the particle is travelling along a line
				tfxWideFloat length = tfxWideAbs(local_position_y);
				tfxWideFloat at_end = tfxWideGreater(length, emitter_size_y);

				tfx__readbarrier;

				local_position_y = tfxWideSub(local_position_y, tfxWideAnd(at_end, offset_y));
				flags = tfxWideOri(flags, tfxWideAndi(tfxWideSetSinglei(tfxParticleFlags_capture_after_transform), tfxWideCasti(at_end)));
				tfxWideStorei((tfxWideIntLoader *)&bank.flags[index], flags);
				tfxWideStore(&bank.position_y[index], local_position_y);
			}
		}
	}

	ControlParticleTransform2d(&pm.work_queue, data);
}

void ControlParticleTransform2d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[work_entry->emitter_index];
	tfx_bounding_box_t &bounding_box = emitter.bounding_box;
	tfx_particle_soa_t &bank = work_entry->pm->particle_arrays[emitter.particles_index];

	tfx_2d_instance_t *sprites = tfxCastBuffer(tfx_2d_instance_t, work_entry->sprite_instances);
	tfxU32 running_sprite_index = work_entry->sprites_index;

	const tfxWideFloat e_world_position_x = tfxWideSetSingle(emitter.world_position.x);
	const tfxWideFloat e_world_position_y = tfxWideSetSingle(emitter.world_position.y);
	const tfxWideFloat e_handle_x = tfxWideSetSingle(emitter.handle.x);
	const tfxWideFloat e_handle_y = tfxWideSetSingle(emitter.handle.y);
	const tfxWideFloat e_scale = tfxWideSetSingle(work_entry->overal_scale);
	const tfxWideInt capture_after_transform = tfxWideSetSinglei(tfxParticleFlags_capture_after_transform);
	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);

	tfxU32 start_diff = work_entry->start_diff;

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;

		tfxWideArray position_x;
		tfxWideArray position_y;
		tfxWideArray captured_position_x;
		tfxWideArray captured_position_y;
		position_x.m = tfxWideLoad(&bank.position_x[index]);
		position_y.m = tfxWideLoad(&bank.position_y[index]);
		captured_position_x.m = tfxWideLoad(&bank.captured_position_x[index]);
		captured_position_y.m = tfxWideLoad(&bank.captured_position_y[index]);
		tfxWideInt flags = tfxWideLoadi((tfxWideIntLoader *)&bank.flags[index]);
		tfxWideFloat capture_flag = tfxWideCast(tfxWideGreateri(tfxWideAndi(flags, capture_after_transform), tfxWideSetZeroi));
		tfxWideFloat xor_capture_flag = tfxWideEquals(capture_flag, tfxWideSetZero);

		const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
		const tfxWideFloat age = tfxWideLoad(&bank.age[index]);
		tfx__readbarrier;
		tfxWideFloat life = tfxWideDiv(age, max_age);
		life = tfxWideMul(life, max_life);
		life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

		if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position) {
			position_x.m = tfxWideAdd(position_x.m, e_handle_x);
			position_y.m = tfxWideAdd(position_y.m, e_handle_y);
			tfx__wide_transform_quaternion_vec2(&emitter.rotation, &position_x.m, &position_y.m);
			position_x.m = tfxWideAdd(tfxWideMul(position_x.m, e_scale), e_world_position_x);
			position_y.m = tfxWideAdd(tfxWideMul(position_y.m, e_scale), e_world_position_y);
		}

		captured_position_x.m = tfxWideAdd(tfxWideAnd(position_x.m, capture_flag), tfxWideAnd(captured_position_x.m, xor_capture_flag));
		captured_position_y.m = tfxWideAdd(tfxWideAnd(position_y.m, capture_flag), tfxWideAnd(captured_position_y.m, xor_capture_flag));

		if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_location_source && emitter.spawn_locations_index != tfxINVALID) {
			tfx_spawn_points_soa_t &locations = work_entry->pm->particle_location_arrays[emitter.spawn_locations_index];
			tfxWideStore(&locations.position_x[index], position_x.m);
			tfxWideStore(&locations.position_y[index], position_y.m);
			tfxWideStore(&locations.captured_position_x[index], captured_position_x.m);
			tfxWideStore(&locations.captured_position_y[index], captured_position_y.m);
			tfxWideStore(&locations.age[index], tfxWideDiv(age, max_age));
		}

		tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
		if (!(pm.flags & tfxParticleManagerFlags_unordered)) {    //Predictable
			for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
				int index_j = index + j;
				tfxU32 sprite_depth_index = bank.depth_index[index_j] + work_entry->cumulative_index_point + work_entry->effect_instance_offset;
				TFX_ASSERT(sprite_depth_index < work_entry->sprite_instances->current_size);
				sprites[sprite_depth_index].position.x = position_x.a[j];
				sprites[sprite_depth_index].position.y = position_y.a[j];
				bank.captured_position_x[index_j] = sprites[sprite_depth_index].position.x;
				bank.captured_position_y[index_j] = sprites[sprite_depth_index].position.y;
				/*
				//Not sure what I'm doing for this yet
				bounding_box.min_corner.x = tfx__Min(position_x.a[j], bounding_box.min_corner.x);
				bounding_box.min_corner.y = tfx__Min(position_y.a[j], bounding_box.min_corner.y);
				bounding_box.max_corner.x = tfx__Max(position_x.a[j], bounding_box.max_corner.x);
				bounding_box.max_corner.y = tfx__Max(position_y.a[j], bounding_box.max_corner.y);
				*/
				running_sprite_index++;
			}
		}
		else {
			for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
				int index_j = index + j;
				TFX_ASSERT(running_sprite_index < work_entry->sprite_instances->current_size);
				sprites[running_sprite_index].position.x = position_x.a[j];
				sprites[running_sprite_index].position.y = position_y.a[j];
				bank.captured_position_x[index_j] = sprites[running_sprite_index].position.x;
				bank.captured_position_y[index_j] = sprites[running_sprite_index].position.y;
				/*
				//Not sure what I'm doing for this yet
				bounding_box.min_corner.x = tfx__Min(position_x.a[j], bounding_box.min_corner.x);
				bounding_box.min_corner.y = tfx__Min(position_y.a[j], bounding_box.min_corner.y);
				bounding_box.max_corner.x = tfx__Max(position_x.a[j], bounding_box.max_corner.x);
				bounding_box.max_corner.y = tfx__Max(position_y.a[j], bounding_box.max_corner.y);
				*/
				running_sprite_index++;
			}
		}
		start_diff = 0;
	}
}

void ControlParticleBoundingBox(tfx_work_queue_t *queue, void *data) {

}

void ControlParticleSpin(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfxU32 emitter_index = work_entry->emitter_index;
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[work_entry->emitter_index];
	const tfx_emission_type emission_type = work_entry->properties->emission_type;
	tfx_particle_soa_t &bank = pm.particle_arrays[emitter.particles_index];

	const tfxWideInt spin_last_frame = tfxWideSetSinglei(work_entry->graphs->spin.lookup.last_frame);

	tfxU32 running_sprite_index = work_entry->sprites_index;

	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);

	tfxU32 start_diff = work_entry->start_diff;

	tfxWideArrayi lookup_frame;
	tfx_2d_instance_t *sprites_2d = tfxCastBuffer(tfx_2d_instance_t, work_entry->sprite_instances);
	tfx_3d_instance_t *sprites_3d = tfxCastBuffer(tfx_3d_instance_t, work_entry->sprite_instances);

	const tfxWideFloat e_world_rotations_z = tfxWideSetSingle(emitter.world_rotations.z);
	bool relative_position = emitter.property_flags & tfxEmitterPropertyFlags_relative_position || (emitter.property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == tfxLine);
	bool is_ordered = (!(pm.flags & tfxParticleManagerFlags_unordered) || (IsOrderedEffectState(&pm.effects[emitter.root_index])));

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;

		const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
		const tfxWideFloat age = tfxWideLoad(&bank.age[index]);
		tfx__readbarrier;
		tfxWideFloat life = tfxWideDiv(age, max_age);
		life = tfxWideMul(life, max_life);
		life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

		const tfxWideFloat base_spin = tfxWideLoad(&bank.base_spin[index]);

		tfxWideArray rotations_z;
		rotations_z.m = tfxWideLoad(&bank.local_rotations_z[index]);

		lookup_frame.m = tfxWideMini(tfxWideConverti(life), spin_last_frame);
		const tfxWideFloat lookup_roll_spin = tfxWideMul(tfxWideLookupSet(work_entry->graphs->spin.lookup.values, lookup_frame), base_spin);

		//----Spin and angle Changes
		rotations_z.m = tfxWideAdd(rotations_z.m, tfxWideMul(lookup_roll_spin, pm.update_time_wide));
		tfxWideStore(&bank.local_rotations_z[index], rotations_z.m);

		if (emitter.property_flags & tfxEmitterPropertyFlags_relative_angle) {
			rotations_z.m = tfxWideAdd(rotations_z.m, e_world_rotations_z);
		}

		tfx__readbarrier;

		tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
		if (pm.flags & tfxParticleManagerFlags_3d_effects) { //Predictable
			if (is_ordered) {    //Predictable
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 sprite_depth_index = bank.depth_index[index + j] + work_entry->cumulative_index_point + work_entry->effect_instance_offset;
					TFX_ASSERT(sprite_depth_index < work_entry->sprite_instances->current_size);
					sprites_3d[sprite_depth_index].rotations.x = 0.f;
					sprites_3d[sprite_depth_index].rotations.y = 0.f;
					sprites_3d[sprite_depth_index].rotations.z = rotations_z.a[j];
					running_sprite_index++;
				}
			}
			else {
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					TFX_ASSERT(running_sprite_index < work_entry->sprite_instances->current_size);
					sprites_3d[running_sprite_index].rotations.x = 0.f;
					sprites_3d[running_sprite_index].rotations.y = 0.f;
					sprites_3d[running_sprite_index++].rotations.z = rotations_z.a[j];
				}
			}
		}
		else {
			if (is_ordered) {        //Predictable
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 sprite_depth_index = bank.depth_index[index + j] + work_entry->cumulative_index_point + work_entry->effect_instance_offset;
					TFX_ASSERT(sprite_depth_index < work_entry->sprite_instances->current_size);
					sprites_2d[sprite_depth_index].position.w = rotations_z.a[j];
					running_sprite_index++;
				}
			}
			else {
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					TFX_ASSERT(running_sprite_index < work_entry->sprite_instances->current_size);
					sprites_2d[running_sprite_index++].position.w = rotations_z.a[j];
				}
			}
		}
		start_diff = 0;
	}
}

void ControlParticleSpin3d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfxU32 emitter_index = work_entry->emitter_index;
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[work_entry->emitter_index];
	tfx_particle_soa_t &bank = pm.particle_arrays[emitter.particles_index];
	const tfx_emission_type emission_type = work_entry->properties->emission_type;

	const tfxWideInt spin_last_frame = tfxWideSetSinglei(work_entry->graphs->spin.lookup.last_frame);
	const tfxWideInt spin_pitch_last_frame = tfxWideSetSinglei(work_entry->graphs->pitch_spin.lookup.last_frame);
	const tfxWideInt spin_yaw_last_frame = tfxWideSetSinglei(work_entry->graphs->yaw_spin.lookup.last_frame);

	tfxU32 running_sprite_index = work_entry->sprites_index;

	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);

	tfxU32 start_diff = work_entry->start_diff;

	tfxWideArrayi lookup_frame;
	tfx_3d_instance_t *sprites = tfxCastBuffer(tfx_3d_instance_t, work_entry->sprite_instances);

	const tfxWideFloat e_world_rotations_x = tfxWideSetSingle(emitter.world_rotations.x);
	const tfxWideFloat e_world_rotations_y = tfxWideSetSingle(emitter.world_rotations.y);
	const tfxWideFloat e_world_rotations_z = tfxWideSetSingle(emitter.world_rotations.z);
	bool relative_position = emitter.property_flags & tfxEmitterPropertyFlags_relative_position || (emitter.property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == tfxLine);

	bool is_ordered = (!(pm.flags & tfxParticleManagerFlags_unordered) || (IsOrderedEffectState(&pm.effects[emitter.root_index])));

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;

		const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
		const tfxWideFloat age = tfxWideLoad(&bank.age[index]);
		tfx__readbarrier;
		tfxWideFloat life = tfxWideDiv(age, max_age);
		life = tfxWideMul(life, max_life);
		life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

		const tfxWideFloat base_spin = tfxWideLoad(&bank.base_spin[index]);
		const tfxWideFloat base_pitch_spin = tfxWideLoad(&bank.base_pitch_spin[index]);
		const tfxWideFloat base_yaw_spin = tfxWideLoad(&bank.base_yaw_spin[index]);

		tfxWideArray rotations_x;
		tfxWideArray rotations_y;
		tfxWideArray rotations_z;
		rotations_x.m = tfxWideLoad(&bank.local_rotations_x[index]);
		rotations_y.m = tfxWideLoad(&bank.local_rotations_y[index]);
		rotations_z.m = tfxWideLoad(&bank.local_rotations_z[index]);

		lookup_frame.m = tfxWideMini(tfxWideConverti(life), spin_last_frame);
		const tfxWideFloat lookup_roll_spin = tfxWideMul(tfxWideLookupSet(work_entry->graphs->spin.lookup.values, lookup_frame), base_spin);
		lookup_frame.m = tfxWideMini(tfxWideConverti(life), spin_pitch_last_frame);
		const tfxWideFloat lookup_pitch_spin = tfxWideMul(tfxWideLookupSet(work_entry->graphs->pitch_spin.lookup.values, lookup_frame), base_pitch_spin);
		lookup_frame.m = tfxWideMini(tfxWideConverti(life), spin_yaw_last_frame);
		const tfxWideFloat lookup_yaw_spin = tfxWideMul(tfxWideLookupSet(work_entry->graphs->yaw_spin.lookup.values, lookup_frame), base_yaw_spin);

		//----Spin and angle Changes
		rotations_x.m = tfxWideAdd(rotations_x.m, tfxWideMul(lookup_pitch_spin, pm.update_time_wide));
		rotations_y.m = tfxWideAdd(rotations_y.m, tfxWideMul(lookup_yaw_spin, pm.update_time_wide));
		rotations_z.m = tfxWideAdd(rotations_z.m, tfxWideMul(lookup_roll_spin, pm.update_time_wide));

		tfxWideStore(&bank.local_rotations_x[index], rotations_x.m);
		tfxWideStore(&bank.local_rotations_y[index], rotations_y.m);
		tfxWideStore(&bank.local_rotations_z[index], rotations_z.m);

		if (!relative_position && emitter.property_flags & tfxEmitterPropertyFlags_relative_angle) {
			rotations_x.m = tfxWideAdd(rotations_x.m, e_world_rotations_x);
			rotations_y.m = tfxWideAdd(rotations_y.m, e_world_rotations_y);
			rotations_z.m = tfxWideAdd(rotations_z.m, e_world_rotations_z);
		}

		tfx__readbarrier;

		tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
		if (is_ordered) {    //Predictable
			for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
				tfxU32 sprite_depth_index = bank.depth_index[index + j] + work_entry->cumulative_index_point + work_entry->effect_instance_offset;
				TFX_ASSERT(sprite_depth_index < work_entry->sprite_instances->current_size);
				sprites[sprite_depth_index].rotations.x = rotations_x.a[j];
				sprites[sprite_depth_index].rotations.y = rotations_y.a[j];
				sprites[sprite_depth_index].rotations.z = rotations_z.a[j];
				running_sprite_index++;
			}
		}
		else {
			for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
				TFX_ASSERT(running_sprite_index < work_entry->sprite_instances->current_size);
				sprites[running_sprite_index].rotations.x = rotations_x.a[j];
				sprites[running_sprite_index].rotations.y = rotations_y.a[j];
				sprites[running_sprite_index].rotations.z = rotations_z.a[j];
				running_sprite_index++;
			}
		}
		start_diff = 0;
	}
}

void ControlParticleSize(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfxU32 emitter_index = work_entry->emitter_index;
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[work_entry->emitter_index];
	tfx_particle_soa_t &bank = pm.particle_arrays[emitter.particles_index];

	const tfxWideInt width_last_frame = tfxWideSetSinglei(work_entry->graphs->width.lookup.last_frame);
	const tfxWideInt height_last_frame = tfxWideSetSinglei(work_entry->graphs->height.lookup.last_frame);
	const tfxWideFloat overal_scale = tfxWideSetSingle(work_entry->overal_scale);

	tfxU32 running_sprite_index = work_entry->sprites_index;

	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);

	tfxU32 start_diff = work_entry->start_diff;

	tfxWideArrayi lookup_frame;
	tfxWideFloat scale_x;
	tfxWideFloat scale_y;

	tfx_emitter_path_t *path;
	tfxWideFloat life;
	float packed_scale_amount = pm.flags & tfxParticleManagerFlags_3d_effects ? 256.f : 8192.f;

	bool sample_based_on_path_position = emitter.property_flags & tfxEmitterPropertyFlags_alt_size_lifetime_sampling && work_entry->properties->emission_type == tfxPath;

	if (sample_based_on_path_position) {
		path = &pm.library->paths[emitter.path_attributes];
	}
	bool is_ordered = (!(pm.flags & tfxParticleManagerFlags_unordered) || (IsOrderedEffectState(&pm.effects[emitter.root_index])));

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;

		const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
		const tfxWideFloat age = tfxWideLoad(&bank.age[index]);
		tfx__readbarrier;

		if (sample_based_on_path_position) {
			const tfxWideFloat path_position = tfxWideLoad(&bank.path_position[index]);
			life = tfxWideDiv(path_position, tfxWideSetSingle(path->node_count - 3.f));
			life = tfxWideMul(life, max_life);
		}
		else {
			life = tfxWideDiv(age, max_age);
			life = tfxWideMul(life, max_life);
		}

		lookup_frame.m = tfxWideMini(tfxWideConverti(life), width_last_frame);
		const tfxWideFloat lookup_width = tfxWideLookupSet(work_entry->graphs->width.lookup.values, lookup_frame);

		lookup_frame.m = tfxWideMini(tfxWideConverti(life), height_last_frame);
		const tfxWideFloat lookup_height = tfxWideLookupSet(work_entry->graphs->height.lookup.values, lookup_frame);

		const tfxWideFloat base_size_x = tfxWideLoad(&bank.base_size_x[index]);
		const tfxWideFloat base_size_y = tfxWideLoad(&bank.base_size_y[index]);

		tfx__readbarrier;

		//----Size Changes
		scale_x = tfxWideMul(base_size_x, lookup_width);

		if (emitter.state_flags & tfxEmitterStateFlags_lifetime_uniform_size) {
			scale_y = tfxWideMul(lookup_width, base_size_y);
			if (emitter.state_flags & tfxEmitterPropertyFlags_base_uniform_size)
				scale_y = tfxWideMin(scale_x, scale_y);
		}
		else
			scale_y = tfxWideMul(lookup_height, base_size_y);

		scale_x = tfxWideMul(scale_x, overal_scale);
		scale_y = tfxWideMul(scale_y, overal_scale);

		tfxWideArrayi packed_scale = { tfx__wide_pack16bit_2sscaled(scale_x, scale_y, packed_scale_amount) };

		tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
		if (pm.flags & tfxParticleManagerFlags_3d_effects) { //Predictable
			tfx_3d_instance_t *sprites = tfxCastBuffer(tfx_3d_instance_t, work_entry->sprite_instances);
			if (is_ordered) {    //Predictable
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 sprite_depth_index = bank.depth_index[index + j] + work_entry->cumulative_index_point + work_entry->effect_instance_offset;
					TFX_ASSERT(sprite_depth_index < work_entry->sprite_instances->current_size);
					sprites[sprite_depth_index].size_handle.packed = packed_scale.a[j] | emitter.image_handle_packed;
					running_sprite_index++;
				}
			}
			else {
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					TFX_ASSERT(running_sprite_index < work_entry->sprite_instances->current_size);
					tfx_float16x4_t scale_handle;
					scale_handle.packed = packed_scale.a[j] | emitter.image_handle_packed;
					sprites[running_sprite_index].size_handle.packed = scale_handle.packed;
					running_sprite_index++;
				}
			}
		}
		else {
			tfx_2d_instance_t *sprites = tfxCastBuffer(tfx_2d_instance_t, work_entry->sprite_instances);
			if (is_ordered) {    //Predictable
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 sprite_depth_index = bank.depth_index[index + j] + work_entry->cumulative_index_point + work_entry->effect_instance_offset;
					TFX_ASSERT(sprite_depth_index < work_entry->sprite_instances->current_size);
					sprites[sprite_depth_index].size_handle.packed = packed_scale.a[j] | emitter.image_handle_packed;
					running_sprite_index++;
				}
			}
			else {
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					TFX_ASSERT(running_sprite_index < work_entry->sprite_instances->current_size);
					sprites[running_sprite_index++].size_handle.packed = packed_scale.a[j] | emitter.image_handle_packed;
				}
			}
		}
		start_diff = 0;
	}
}

void ControlParticleColor(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[work_entry->emitter_index];
	tfx_particle_soa_t &bank = work_entry->pm->particle_arrays[emitter.particles_index];

	const tfxWideFloat global_intensity = tfxWideSetSingle(work_entry->global_intensity);

	tfxU32 running_sprite_index = work_entry->sprites_index;

	tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
	tfxU32 start_diff = work_entry->start_diff;

	tfx_color_ramp_t &ramp = work_entry->graphs->color_ramps[0];
	tfxWideFloat life;
	tfxWideArrayi ramp_index;
	tfxWideFloat color_ramp_size = tfxWideSetSingle(tfxCOLOR_RAMP_WIDTH);
	tfxWideInt color_ramp_sizei = tfxWideSetSinglei(tfxCOLOR_RAMP_WIDTH - 1);

	bool sample_based_on_path_position = emitter.property_flags & tfxEmitterPropertyFlags_alt_color_lifetime_sampling && work_entry->properties->emission_type == tfxPath;

	tfx_emitter_path_t *path;
	if (sample_based_on_path_position) {
		path = &pm.library->paths[emitter.path_attributes];
	}
	bool is_ordered = (!(pm.flags & tfxParticleManagerFlags_unordered) || (IsOrderedEffectState(&pm.effects[emitter.root_index])));
	bool is_mixed_color = emitter.property_flags & tfxEmitterPropertyFlags_use_color_hint;
	tfxWideArrayi curved_alpha;

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;

		const tfxWideFloat intensity_factor = tfxWideLoad(&bank.intensity_factor[index]);
		const tfxWideFloat age = tfxWideLoad(&bank.age[index]);
		const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
		tfx__readbarrier;

		if (sample_based_on_path_position) {
			const tfxWideFloat path_position = tfxWideLoad(&bank.path_position[index]);
			life = tfxWideDiv(path_position, tfxWideSetSingle(path->node_count - 3.f));
		}
		else {
			life = tfxWideDiv(age, max_age);
		}
		ramp_index.m = tfxWideMini(tfxWideConverti(tfxWideMul(life, color_ramp_size)), color_ramp_sizei);
		tfxWideFloat lookup_intensity = tfxWideLookupSet(work_entry->graphs->intensity.lookup.values, ramp_index);
		tfxWideFloat dissolve_lerp = tfxWideLookupSet(work_entry->graphs->curved_alpha.lookup.values, ramp_index);
		tfxWideFloat sharpness = tfxWideLookupSet(work_entry->graphs->alpha_sharpness.lookup.values, ramp_index);

		//----Color changes
		lookup_intensity = tfxWideMul(tfxWideMul(global_intensity, lookup_intensity), intensity_factor);

		tfxWideArrayi packed_intensity_life = { tfx__wide_pack16bit_2sscaled(lookup_intensity, life, 128.f) };
        curved_alpha.m = tfx__wide_pack16bit(dissolve_lerp, sharpness);

		tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
		if (pm.flags & tfxParticleManagerFlags_3d_effects) { //Predictable
			tfx_3d_instance_t *sprites = tfxCastBuffer(tfx_3d_instance_t, work_entry->sprite_instances);
			if (is_ordered) {
				WriteParticleColorSpriteDataOrdered(sprites, pm, work_entry->layer, start_diff, limit_index, bank.depth_index, index, packed_intensity_life, curved_alpha, running_sprite_index, work_entry->effect_instance_offset + work_entry->cumulative_index_point);
			}
			else {
				WriteParticleColorSpriteData(sprites, start_diff, limit_index,bank.depth_index, index, packed_intensity_life, curved_alpha, running_sprite_index);
			}
		}
		else {
			tfx_2d_instance_t *sprites = tfxCastBuffer(tfx_2d_instance_t, work_entry->sprite_instances);
			if (is_ordered) {
				WriteParticleColorSpriteDataOrdered(sprites, pm, work_entry->layer, start_diff, limit_index, bank.depth_index, index, packed_intensity_life, curved_alpha, running_sprite_index, work_entry->effect_instance_offset + work_entry->cumulative_index_point);
			}
			else {
				WriteParticleColorSpriteData(sprites, start_diff, limit_index, bank.depth_index, index, packed_intensity_life, curved_alpha, running_sprite_index);
			}
		}
		start_diff = 0;
	}

}

void ControlParticleImageFrame(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[work_entry->emitter_index];
	tfx_particle_soa_t &bank = pm.particle_arrays[emitter.particles_index];
	tfx_image_data_t *image = work_entry->properties->image;
	const tfx_billboarding_option billboard_option = work_entry->properties->billboard_option;
	tfxU32 layer = work_entry->layer << 28;

	tfxU32 start_diff = work_entry->start_diff;

	tfxWideFloat image_frame_rate = tfxWideSetSingle(emitter.image_frame_rate);
	image_frame_rate = tfxWideMul(image_frame_rate, pm.update_time_wide);
	tfxWideFloat end_frame = tfxWideSetSingle(emitter.end_frame);
	tfxWideFloat frames = tfxWideSetSingle(emitter.end_frame + 1);
	tfxEmitterStateFlags emitter_flags = emitter.state_flags;
	tfxEmitterStateFlags property_flags = emitter.property_flags;
	const tfxWideInt xor_capture_after_transform_flag = tfxWideXOri(tfxWideSetSinglei(tfxParticleFlags_capture_after_transform), tfxWideSetSinglei(-1));
	const tfxWideInt capture_after_transform_flag = tfxWideSetSinglei(tfxParticleFlags_capture_after_transform);

	tfxWideInt color_ramp_indexes;
	color_ramp_indexes = tfxWideSetSinglei(tfxColorRampIndex(work_entry->graphs->color_ramp_bitmap_indexes[0]) << 24);
	color_ramp_indexes = tfxWideOri(color_ramp_indexes, tfxWideSetSinglei(tfxColorRampLayer(work_entry->graphs->color_ramp_bitmap_indexes[0]) << 16));
	tfxWideInt image_start_index = tfxWideSetSinglei((pm.flags & tfxParticleManagerFlags_recording_sprites) && (pm.flags & tfxParticleManagerFlags_using_uids) ? 0 : image->compute_shape_index);

	tfxU32 running_sprite_index = work_entry->sprites_index;
	bool is_ordered = (!(pm.flags & tfxParticleManagerFlags_unordered) || (IsOrderedEffectState(&pm.effects[emitter.root_index])));

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;

		tfxWideArray image_frame;
		image_frame.m = tfxWideLoad(&bank.image_frame[index]);
		tfxWideArrayi flags;
		tfxWideArrayi single_loop_count;
		single_loop_count.m = tfxWideLoadi((tfxWideIntLoader *)&bank.single_loop_count[index]);
		//We only want to not capture if single loop count is 0.
		flags.m = tfxWideLoadi((tfxWideIntLoader *)&bank.flags[index]);
		flags.m = tfxWideXOri(tfxWideAndi(flags.m, capture_after_transform_flag), capture_after_transform_flag);
		//flags.m = tfxWideOri(flags.m, tfxWideAndi(capture_after_transform_flag, tfxWideGreateri(single_loop_count.m, tfxWideSetZeroi)));

		tfx__readbarrier;

		//----Image animation
		image_frame.m = tfxWideAdd(image_frame.m, image_frame_rate);
		tfxWideStore(&bank.image_frame[index], image_frame.m);
		if (property_flags & tfxEmitterPropertyFlags_reverse_animation && emitter.state_flags & tfxEmitterStateFlags_play_once) {
			image_frame.m = tfxWideSub(end_frame, image_frame.m);
			image_frame.m = tfxWideMax(image_frame.m, tfxWideSetZero);
		}
		else if (emitter.state_flags & tfxEmitterStateFlags_play_once) {
			image_frame.m = tfxWideMin(image_frame.m, end_frame);
		}
		else if (property_flags & tfxEmitterPropertyFlags_reverse_animation) {
			image_frame.m = tfxWideMod(image_frame.m, frames);
			image_frame.m = tfxWideSub(end_frame, image_frame.m);
		}
		else {
			image_frame.m = tfxWideMod(image_frame.m, frames);
		}

		tfxWideArrayi image_indexes = { tfxWideOri(color_ramp_indexes, tfxWideAddi(tfxWideConverti(image_frame.m), image_start_index)) };

		tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
		if (pm.flags & tfxParticleManagerFlags_3d_effects) { //Predictable
			tfx_3d_instance_t *sprites = tfxCastBuffer(tfx_3d_instance_t, work_entry->sprite_instances);
			if (is_ordered) {
				WriteParticleImageSpriteDataOrdered(sprites, pm, layer, start_diff, limit_index, bank, flags, image_indexes, emitter_flags, billboard_option, index, running_sprite_index, work_entry->effect_instance_offset + work_entry->cumulative_index_point);
			} else {
				WriteParticleImageSpriteData(sprites, pm, layer, start_diff, limit_index, bank, flags, image_indexes, emitter_flags, billboard_option, index, running_sprite_index);
			}
		} else {
			tfx_2d_instance_t *sprites = tfxCastBuffer(tfx_2d_instance_t, work_entry->sprite_instances);
			if (is_ordered) {
				WriteParticleImageSpriteDataOrdered(sprites, pm, layer, start_diff, limit_index, bank, flags, image_indexes, emitter_flags, billboard_option, index, running_sprite_index, work_entry->effect_instance_offset + work_entry->cumulative_index_point);
			} else {
				WriteParticleImageSpriteData(sprites, pm, layer, start_diff, limit_index, bank, flags, image_indexes, emitter_flags, billboard_option, index, running_sprite_index);
			}
		}

		//We can't alter the flags here, there was a bug where storing a block of 4 here would unflag a particle when it shouldn't but I'm not entirely sure why.
		//Something to do with limit_index being less then 4 and so a particle in the bank (which in theory shouldn't exist yet) gets unflagged. limit_index is
		//only ever less than 4 on the final loop so maybe something carries over to the next frame. 
		//It generally only happened when the spawn rate is increased in the editor which may have been throwing things off. For now the flagging is done inside the sprite writing
		//so we know we're only marking the flag for each particle that is definitely in play. The alternative that also works is a separate loop below, need to 
		//test which is faster but it would be nice to know exactly why this happens in the first place.
		//flags.m = tfxWideAndi(flags.m, xor_capture_after_transform_flag);
		//tfxWideStorei((tfxWideIntLoader*)&bank.flags[index], flags.m);

		start_diff = 0;
	}

	/*
	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;
		tfxWideInt flags = tfxWideLoadi((tfxWideInt*)&bank.flags[index]);
		flags = tfxWideAndi(flags, xor_capture_after_transform_flag);
		tfxWideStorei((tfxWideIntLoader*)&bank.flags[index], flags);
	}
	*/

}

void ControlParticleUID(tfx_work_queue_t *queue, void *data) {
	//This function is only run when recording sprite data
	tfxPROFILE;
	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[work_entry->emitter_index];
	tfx_particle_soa_t &bank = pm.particle_arrays[emitter.particles_index];

	tfxU32 start_diff = work_entry->start_diff;

	tfxU32 running_sprite_index = work_entry->sprites_index;
	tfx_vector_t<tfx_unique_sprite_id_t> &sprite_uids = pm.unique_sprite_ids[pm.current_sprite_buffer][work_entry->layer];
	tfx_2d_instance_t *instance = tfxCastBuffer(tfx_2d_instance_t, work_entry->sprite_instances);
	bool is_ordered = (!(pm.flags & tfxParticleManagerFlags_unordered) || (IsOrderedEffectState(&pm.effects[emitter.root_index])));
	bool is_wrapped = emitter.state_flags & tfxEmitterStateFlags_wrap_single_sprite;

	for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;

		tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
		if (is_ordered) {                //Predictable
			for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
				int index_j = index + j;
				tfxU32 sprite_depth_index = bank.depth_index[index_j] + work_entry->cumulative_index_point + work_entry->effect_instance_offset;
				bool new_id = bank.age[index_j] == 0 && bank.single_loop_count[index_j] > 0 && !is_wrapped ? true : false;
				sprite_uids[sprite_depth_index].uid = new_id ? (tfxU32)tfx__rdtsc() : bank.uid[index_j];
				bank.uid[index_j] = sprite_uids[sprite_depth_index].uid;
				sprite_uids[sprite_depth_index].age = tfxU32((bank.age[index_j] + 0.1f) / pm.frame_length);
				sprite_uids[sprite_depth_index].property_index = emitter.properties_index;
				running_sprite_index++;
			}
		}
		else {
			for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
				int index_j = index + j;
				bool new_id = bank.age[index_j] == 0 && bank.single_loop_count[index_j] > 0 && !is_wrapped ? true : false;
				sprite_uids[running_sprite_index].uid = new_id ? (tfxU32)tfx__rdtsc() : bank.uid[index_j];
				bank.uid[index_j] = sprite_uids[running_sprite_index].uid;
				sprite_uids[running_sprite_index].age = tfxU32((bank.age[index_j] + 0.1f) / pm.frame_length);
				sprite_uids[running_sprite_index].property_index = emitter.properties_index;
				running_sprite_index++;
			}
		}
		start_diff = 0;
	}
}

tfx_vector_t<tfx_effect_index_t> *GetPMEffectBuffer(tfx_particle_manager_t *pm, tfxU32 depth) {
	return &pm->effects_in_use[depth][pm->current_ebuff];
}

tfx_vector_t<tfxU32> *GetPMEmitterBuffer(tfx_particle_manager_t *pm, tfxU32 depth) {
	return &pm->control_emitter_queue;
}

void ToggleSpritesWithUID(tfx_particle_manager_t *pm, bool switch_on) {
	if (switch_on) {
		for (tfxEachLayer) {
			pm->unique_sprite_ids[0][layer].reserve(pm->max_cpu_particles_per_layer[layer]);
			if (!(pm->flags & tfxParticleManagerFlags_use_effect_sprite_buffers)) {
				pm->unique_sprite_ids[1][layer].reserve(pm->max_cpu_particles_per_layer[layer]);
			}
		}
		pm->flags |= tfxParticleManagerFlags_using_uids;
	}
	else {
		pm->flags &= ~tfxParticleManagerFlags_using_uids;
	}
}

void ReconfigureParticleManager(tfx_particle_manager_t *pm, tfx_particle_manager_mode mode, tfxU32 req_sort_passes, bool is_3d) {
	ClearParticleManager(pm, true, true);
	for (auto &bank : pm->free_particle_lists.data) {
		bank.free_all();
	}
	pm->free_particle_lists.FreeAll();

	tfxParticleManagerFlags current_flags = (pm->flags & tfxParticleManagerFlags_dynamic_sprite_allocation) | 
											(pm->flags & tfxParticleManagerFlags_double_buffer_sprites) | 
											(pm->flags & tfxParticleManagerFlags_2d_and_3d) | 
											(pm->flags & tfxParticleManagerFlags_use_effect_sprite_buffers) |
											(pm->flags & tfxParticleManagerFlags_direct_to_staging_buffer);

	if (mode == tfxParticleManagerMode_unordered)
		pm->flags = tfxParticleManagerFlags_unordered;
	else if (mode == tfxParticleManagerMode_ordered_by_depth)
		pm->flags = tfxParticleManagerFlags_order_by_depth;
	else if (mode == tfxParticleManagerMode_ordered_by_depth_guaranteed)
		pm->flags = tfxParticleManagerFlags_order_by_depth | tfxParticleManagerFlags_guarantee_order;
	else if (mode == tfxParticleManagerMode_ordered_by_age)
		pm->flags = tfxParticleManagerFlags_ordered_by_age;

	tfxU32 size_in_bytes = pm->instance_buffer.capacity * pm->instance_buffer.struct_size;
	if (is_3d) {
		pm->flags |= tfxParticleManagerFlags_3d_effects;
		tfxReconfigureBuffer(&pm->instance_buffer, sizeof(tfx_3d_instance_t));
	} else {
		pm->flags &= ~tfxParticleManagerFlags_3d_effects;
		tfxReconfigureBuffer(&pm->instance_buffer, sizeof(tfx_2d_instance_t));
	}

	pm->flags |= current_flags;

	pm->instance_buffer.clear();
	for (tfxEachLayer) {
		pm->unique_sprite_ids[0][layer].clear();
		if (pm->flags & tfxParticleManagerFlags_double_buffer_sprites) {
			pm->unique_sprite_ids[0][layer].clear();
		}
	}

	memset(pm->sprite_index_point, 0, 4 * tfxLAYERS);

	pm->emitters.clear();
	pm->effects.clear();

	pm->sort_passes = req_sort_passes;
}

void SetStagingBuffer(tfx_particle_manager_t *pm, void *staging_buffer, tfxU32 size_in_bytes) {
	TFX_ASSERT(pm->flags & tfxParticleManagerFlags_direct_to_staging_buffer);		//Particle manager must be flagged to write direct to staging before on creation
	TFX_ASSERT(staging_buffer);		//Staging buffer is null!
	tfxU32 unit_size = size_in_bytes / pm->instance_buffer.struct_size;
	pm->instance_buffer.data = staging_buffer;
	pm->instance_buffer.capacity = unit_size;
}

void SetPMWorkQueueSizes(tfx_particle_manager_t *pm, tfxU32 spawn_work_max, tfxU32 control_work_max, tfxU32 age_work_max) {
	pm->spawn_work.reserve(spawn_work_max);
	pm->control_work.reserve(control_work_max);
	pm->age_work.reserve(age_work_max);
	pm->sorting_work_entry.reserve(age_work_max);
}

void ClearParticleManager(tfx_particle_manager_t *pm, bool free_particle_banks, bool free_sprite_buffers) {
	tfxCompleteAllWork(&pm->work_queue);
	if (free_particle_banks) {
		FreeAllParticleLists(pm);
		FreeAllSpawnLocationLists(pm);
		for (auto &list : pm->free_particle_lists.data) {
			list.free_all();
		}
		pm->free_particle_lists.FreeAll();
		for (auto &list : pm->free_particle_location_lists.data) {
			list.free_all();
		}
		pm->free_particle_location_lists.FreeAll();
	}
	else {
		for (tfx_soa_buffer_t &buffer : pm->particle_array_buffers) {
			ClearSoABuffer(&buffer);
		}
		for (tfx_soa_buffer_t &buffer : pm->particle_location_buffers) {
			ClearSoABuffer(&buffer);
		}
		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			for (tfx_effect_index_t index : pm->effects_in_use[depth][pm->current_ebuff]) {
				tfx_effect_state_t &effect = pm->effects[index.index];
				for (tfxU32 emitter_index : effect.emitter_indexes[pm->current_ebuff]) {
					FreeParticleList(pm, emitter_index);
					if (pm->emitters[emitter_index].spawn_locations_index != tfxINVALID) {
						FreeSpawnLocationList(pm, emitter_index);
					}
				}
				effect.emitter_indexes[0].clear();
				effect.emitter_indexes[1].clear();
			}
		}
	}
	for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
		pm->effects_in_use[depth][0].clear();
		pm->effects_in_use[depth][1].clear();
	}
	pm->control_emitter_queue.clear();
	for (tfx_effect_state_t &effect : pm->effects) {
		for (tfxEachLayer) {
			effect.instance_data.depth_indexes[layer][0].free();
			effect.instance_data.depth_indexes[layer][1].free();
		}
		effect.emitter_indexes[0].clear();
		effect.emitter_indexes[1].clear();
	}
	pm->free_effects.clear();
	pm->free_emitters.clear();
	pm->particle_indexes.clear();
	pm->free_particle_indexes.clear();
	pm->emitters.clear();
	pm->effects.clear();
	pm->particle_id = 0;
	pm->spawn_work.clear();
	pm->control_work.clear();
	pm->age_work.clear();
	for (int i = 0; i != pm->path_quaternions.current_size; ++i) {
		if (pm->path_quaternions[i]) {
			tfxFREE(pm->path_quaternions[i]);
		}
	}
	pm->path_quaternions.clear();
	pm->free_path_quaternions.clear();
	pm->instance_buffer.clear();
	for (tfxEachLayer) {
		pm->instance_buffer_for_recording[0][layer].clear();
		if (pm->flags & tfxParticleManagerFlags_double_buffer_sprites) {
			pm->instance_buffer_for_recording[1][layer].clear();
		}
	}
}

void FreeParticleManager(tfx_particle_manager_t *pm) {
	FreeAllParticleLists(pm);
	for (auto &list : pm->free_particle_lists.data) {
		list.free_all();
	}
	pm->free_particle_lists.FreeAll();
	for (auto &list : pm->free_particle_location_lists.data) {
		list.free_all();
	}
	pm->free_particle_location_lists.FreeAll();

	for (tfxEachLayer) {
		pm->unique_sprite_ids[0][layer].free();
		if (pm->flags & tfxParticleManagerFlags_double_buffer_sprites) {
			pm->unique_sprite_ids[1][layer].free();
		}
	}
	pm->instance_buffer.free();
	for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
		pm->effects_in_use[depth][0].free();
		pm->effects_in_use[depth][1].free();
	}
	pm->control_emitter_queue.free();
	for (tfx_effect_state_t &effect : pm->effects) {
		for (tfxEachLayer) {
			effect.instance_data.depth_indexes[layer][0].free();
			effect.instance_data.depth_indexes[layer][1].free();
		}
		effect.emitter_indexes[0].free();
		effect.emitter_indexes[1].free();
	}
	pm->emitters_check_capture.free();
	pm->free_effects.free();
	pm->free_emitters.free();
	pm->particle_indexes.free();
	pm->free_particle_indexes.free();
	pm->effects.free();
	pm->emitters.free();
	pm->free_compute_controllers.free();
	pm->particle_id = 0;
	pm->spawn_work.free();
	pm->control_work.free();
	pm->age_work.free();
	for (int i = 0; i != pm->path_quaternions.current_size; ++i) {
		if (pm->path_quaternions[i]) {
			tfxFREE(pm->path_quaternions[i]);
		}
	}
	pm->path_quaternions.free();
	pm->free_path_quaternions.free();
}

void FreeAllParticleLists(tfx_particle_manager_t *pm) {
	for (auto &buffer : pm->particle_array_buffers) {
		FreeSoABuffer(&buffer);
	}
	pm->particle_array_buffers.clear();
	pm->particle_arrays.clear();
}

void FreeAllSpawnLocationLists(tfx_particle_manager_t *pm) {
	for (auto &buffer : pm->particle_location_buffers) {
		FreeSoABuffer(&buffer);
	}
	pm->particle_location_buffers.clear();
	pm->particle_location_arrays.clear();
}

void SoftExpireAll(tfx_particle_manager_t *pm) {
	for (tfx_effect_index_t effect_index : pm->effects_in_use[0][pm->current_ebuff]) {
		pm->emitters[effect_index.index].state_flags |= tfxEmitterStateFlags_stop_spawning;
	}
}

void SetPMLookUpMode(tfx_particle_manager_t *pm, tfx_lookup_mode mode) {
	if (mode == tfxPrecise) {
		lookup_overtime_callback = LookupPreciseOvertime;
		lookup_callback = LookupPrecise;
	}
	else {
		lookup_overtime_callback = LookupFastOvertime;
		lookup_callback = LookupFast;
	}
	pm->lookup_mode = mode;
}

void KeepBoundingBoxesUpdated(tfx_particle_manager_t *pm, bool yesno) {
	if (yesno) {
		pm->flags |= tfxParticleManagerFlags_update_bounding_boxes;
	}
	else {
		pm->flags &= ~tfxParticleManagerFlags_update_bounding_boxes;
	}
}

void UpdatePMBaseValues(tfx_particle_manager_t *pm) {
	pm->flags |= tfxParticleManagerFlags_update_base_values;
}

bool FreePMEffectCapacity(tfx_particle_manager_t *pm) {
	return pm->effects.current_size < pm->max_effects;
}

tfx_effect_index_t GetPMEffectSlot(tfx_particle_manager_t *pm) {
	if (!pm->free_effects.empty()) {
		return pm->free_effects.pop_back();
	}
	if (pm->effects.current_size == pm->effects.capacity) {
		return { tfxINVALID, 0.f };
	}
	pm->effects.current_size++;
	tfx_effect_instance_data_t *instance_data = &pm->effects[pm->effects.current_size - 1].instance_data;
	instance_data->instance_start_index = tfxINVALID;
	memset(instance_data->depth_starting_index, 0, sizeof(tfxU32) * tfxLAYERS);
	memset(instance_data->current_depth_buffer_index, 0, sizeof(tfxU32) * tfxLAYERS);
	instance_data->instance_count = 0;
	for (tfxEachLayer) {
		memset(instance_data->depth_indexes, 0, sizeof(tfx_vector_t<tfx_depth_index_t>) * 2 * tfxLAYERS);
	}
	tfx_effect_state_t &effect = pm->effects.back();
	effect.emitter_indexes[0].init();
	effect.emitter_indexes[1].init();
	return { pm->effects.current_size - 1, 0.f };
}

tfxU32 GetPMEmitterSlot(tfx_particle_manager_t *pm) {
	if (!pm->free_emitters.empty()) {
		return pm->free_emitters.pop_back();
	}
	if (pm->emitters.current_size == pm->emitters.capacity) {
		return tfxINVALID;
	}
	pm->emitters.current_size++;
	return pm->emitters.current_size - 1;
}

tfxU32 AllocatePathQuaterion(tfx_particle_manager_t *pm, tfxU32 amount) {
	tfx_path_quaternion_t *q = (tfx_path_quaternion_t *)tfxALLOCATE(sizeof(tfx_path_quaternion_t) * amount);
	if (!pm->free_path_quaternions.empty()) {
		tfxU32 free_index = pm->free_path_quaternions.pop_back();
		TFX_ASSERT(pm->path_quaternions[free_index] == nullptr);        //Free path quaternion should be null! For some reason the path was not freed before being added to the the list of free path quaternions
		//or the path was allocated outside of this function.
		pm->path_quaternions[free_index] = q;
		return free_index;
	}
	else {
		pm->path_quaternions.locked_push_back(q);
	}
	return pm->path_quaternions.current_size - 1;
}

void FreePathQuaternion(tfx_particle_manager_t *pm, tfxU32 index) {
	if (pm->path_quaternions[index] != nullptr) {
		tfxFREE(pm->path_quaternions[index]);
		pm->path_quaternions[index] = nullptr;
		pm->free_path_quaternions.push_back(index);
	}
}

tfxU32 GetPMParticleIndexSlot(tfx_particle_manager_t *pm, tfxParticleID particle_id) {
	//Todo: ideally we want a better thread safe container for this
	std::lock_guard<std::mutex> lock(pm->particle_index_mutex);
	if (!pm->free_particle_indexes.empty()) {
		pm->particle_indexes[pm->free_particle_indexes.back()] = particle_id;
		return pm->free_particle_indexes.pop_back();
	}
	pm->particle_indexes.push_back(particle_id);
	return pm->particle_indexes.current_size - 1;
}

void FreePMParticleIndex(tfx_particle_manager_t *pm, tfxU32 *index) {
	//Todo: ideally we want a better thread safe container for this
	//Note: since changing the threading to thread per emitter this is probably not needed now?
	std::lock_guard<std::mutex> lock(pm->particle_index_mutex);
	pm->particle_indexes[*index] = tfxINVALID;
	pm->free_particle_indexes.push_back(*index);
	*index = tfxINVALID;
}

tfxU32 PushPMDepthIndex(tfx_vector_t<tfx_depth_index_t> *depth_indexes, tfx_depth_index_t depth_index) {
	(*depth_indexes).push_back(depth_index);
	return (*depth_indexes).current_size - 1;
}

void SetPMLibrary(tfx_particle_manager_t *pm, tfx_library_t *lib) {
	pm->library = lib;
}

void SetPMCamera(tfx_particle_manager_t *pm, float front[3], float position[3]) {
	pm->camera_front.x = front[0];
	pm->camera_front.y = front[1];
	pm->camera_front.z = front[2];
	pm->camera_position.x = position[0];
	pm->camera_position.y = position[1];
	pm->camera_position.z = position[2];
}

void ResetPMFlags(tfx_particle_manager_t *pm) {
	pm->flags = 0;
}

tfxU32 GetParticleSpriteIndex(tfx_particle_manager_t *pm, tfxParticleID id) {
	return pm->particle_arrays[tfx__particle_bank(id)].sprite_index[tfx__particle_index(id)];
}

unsigned int GetControllerMemoryUsage(tfx_particle_manager_t *pm) {
	return pm->highest_compute_controller_index * sizeof(tfx_compute_controller_t);
}

unsigned int GetParticleMemoryUsage(tfx_particle_manager_t *pm) {
	return pm->new_compute_particle_index * sizeof(tfx_compute_particle_t);
}

void FreeComputeSlot(tfx_particle_manager_t *pm, unsigned int slot_id) {
	pm->free_compute_controllers.push_back(slot_id);
}

tfxU32 ParticleCount(tfx_particle_manager_t *pm) {
	tfxU32 count = 0;
	return pm->instance_buffer.current_size;
}

tfxU32 EffectCount(tfx_particle_manager_t *pm) {
	tfxU32 count = 0;
	for (int d = 0; d != tfxMAXDEPTH; ++d) {
		count += pm->effects_in_use[d][pm->current_ebuff].current_size;
	}
	return count;
}

tfxU32 EmitterCount(tfx_particle_manager_t *pm) {
	tfxU32 count = 0;
	for (int d = 0; d != tfxMAXDEPTH; ++d) {
		count += pm->control_emitter_queue.current_size;
	}
	return count;
}

void ResizeParticleSoACallback(tfx_soa_buffer_t *buffer, tfxU32 index) {
	tfx_particle_soa_t *particles = static_cast<tfx_particle_soa_t *>(buffer->user_data);
	for (int i = index; i != buffer->capacity; ++i) {
		particles->max_age[i] = 1.f;
		particles->age[i] = 1.f;
		particles->flags[i] = 0;
	}
}

tfxU32 tfx__grab_particle_lists(tfx_particle_manager_t *pm, tfxKey emitter_hash, bool is_3d, tfxU32 reserve_amount, tfxEmitterControlProfileFlags flags) {
	if (pm->free_particle_lists.ValidKey(emitter_hash)) {
		tfx_vector_t<tfxU32> &free_banks = pm->free_particle_lists.At(emitter_hash);
		if (free_banks.current_size) {
			pm->particle_array_buffers[free_banks.back()].current_size = 0;
			return free_banks.pop_back();
		}
	}
	tfx_particle_soa_t lists;
	pm->particle_arrays.push_back(lists);
	tfxU32 index = pm->particle_arrays.current_size - 1;
	tfx_soa_buffer_t buffer;
	buffer.resize_callback = ResizeParticleSoACallback;
	buffer.user_data = &pm->particle_arrays.back();
	pm->particle_array_buffers.push_back(buffer);
	TFX_ASSERT(index == pm->particle_array_buffers.current_size - 1);
	if (is_3d) {
		InitParticleSoA3d(&pm->particle_array_buffers[index], &pm->particle_arrays.back(), reserve_amount, flags);
	}
	else {
		InitParticleSoA2d(&pm->particle_array_buffers[index], &pm->particle_arrays.back(), reserve_amount, flags);
	}
	return index;
}

void FreeParticleListsMemory(tfx_particle_manager_t *pm, tfx_effect_emitter_t *emitter) {
	if (pm->free_particle_lists.ValidKey(emitter->path_hash)) {
		tfx_vector_t<tfxU32> &free_banks = pm->free_particle_lists.At(emitter->path_hash);
		for (tfxU32 i : free_banks) {
			FreeSoABuffer(&pm->particle_array_buffers[i]);
		}
		free_banks.free();
	}
	if (pm->free_particle_location_lists.ValidKey(emitter->path_hash)) {
		tfx_vector_t<tfxU32> &free_banks = pm->free_particle_location_lists.At(emitter->path_hash);
		for (tfxU32 i : free_banks) {
			FreeSoABuffer(&pm->particle_location_buffers[i]);
		}
		free_banks.free();
	}
	for (tfx_effect_emitter_t &effect : GetEffectInfo(emitter)->sub_effectors) {
		FreeEffectListsMemory(pm, &effect);
	}
}

tfxAPI void FreeEffectListsMemory(tfx_particle_manager_t *pm, tfx_effect_emitter_t *effect) {
	if (effect->type == tfxFolder) {
		for (tfx_effect_emitter_t &sub_effect : GetEffectInfo(effect)->sub_effectors) {
			FreeEffectListsMemory(pm, &sub_effect);
		}
	}
	for (tfx_effect_emitter_t &emitter : GetEffectInfo(effect)->sub_effectors) {
		FreeParticleListsMemory(pm, &emitter);
	}
}

tfxINTERNAL tfxU32 tfx__grab_particle_location_lists(tfx_particle_manager_t *pm, tfxKey emitter_hash, bool is_3d, tfxU32 reserve_amount) {
	if (pm->free_particle_location_lists.ValidKey(emitter_hash)) {
		tfx_vector_t<tfxU32> &free_banks = pm->free_particle_location_lists.At(emitter_hash);
		if (free_banks.current_size) {
			pm->particle_location_buffers[free_banks.back()].current_size = 0;
			return free_banks.pop_back();
		}
	}

	tfx_spawn_points_soa_t lists;
	tfxU32 index = pm->particle_location_arrays.locked_push_back(lists);
	tfx_soa_buffer_t buffer;
	buffer.resize_callback = nullptr;
	buffer.user_data = nullptr;
	pm->particle_location_buffers.push_back(buffer);
	TFX_ASSERT(index == pm->particle_location_buffers.current_size - 1);
	if (is_3d) {
		InitParticleLocationSoA3d(&pm->particle_location_buffers[index], &pm->particle_location_arrays.back(), reserve_amount);
	}
	else {
		InitParticleLocationSoA2d(&pm->particle_location_buffers[index], &pm->particle_location_arrays.back(), reserve_amount);
	}
	return index;
}

void tfx__gather_stats(tfx_profile_t *profile, tfx_profile_stats_t *stat) {
	stat->cycle_high = 0;
	stat->cycle_low = tfxMAX_64u;
	stat->time_high = 0;
	stat->time_low = tfxMAX_64u;
	stat->hit_count = 0;
	stat->cycle_average = 0;
	stat->time_average = 0;
	for (int i = 0; i != tfxPROFILER_SAMPLES; ++i) {
		tfx_profile_snapshot_t *snap = profile->snapshots + i;
		stat->cycle_high = tfxMax(snap->cycle_count, stat->cycle_high);
		stat->cycle_low = tfxMin(snap->cycle_count, stat->cycle_low);
		stat->cycle_average += snap->cycle_count;
		stat->time_high = tfxMax(snap->run_time, stat->time_high);
		stat->time_low = tfxMin(snap->run_time, stat->time_low);
		stat->time_average += snap->run_time;
		stat->hit_count += snap->hit_count;
	}
	stat->cycle_average /= tfxPROFILER_SAMPLES;
	stat->time_average /= tfxPROFILER_SAMPLES;
	stat->hit_count /= tfxPROFILER_SAMPLES;
}

void tfx__reset_snap_shot(tfx_profile_snapshot_t *snapshot) {
	snapshot->cycle_count = 0;
	snapshot->hit_count = 0;
	snapshot->run_time = 0;
}

void tfx__reset_snap_shots() {
	for (int i = 0; i != tfxPROFILE_COUNT; ++i) {
		tfx_profile_t *profile = tfxProfileArray + i;
		memset(profile->snapshots, 0, tfxPROFILER_SAMPLES * sizeof(tfx_profile_snapshot_t));
	}
}

void tfx__dump_snapshots(tfx_storage_map_t<tfx_vector_t<tfx_profile_snapshot_t>> *profile_snapshots, tfxU32 amount) {
	for (int i = 0; i != tfxPROFILE_COUNT; ++i) {
		tfx_profile_t *profile = tfxProfileArray + i;
		if (!profile->name) {
			tfx__reset_snap_shot(profile->snapshots + i);
			continue;
		}
		if (!profile_snapshots->ValidName(profile->name)) {
			tfx_vector_t<tfx_profile_snapshot_t> snapshots;
			profile_snapshots->Insert(profile->name, snapshots);
		}
		tfx_vector_t<tfx_profile_snapshot_t> &snapshots = profile_snapshots->At(profile->name);
		tfxU32 offset = snapshots.current_size;
		snapshots.resize(snapshots.current_size + tfxPROFILER_SAMPLES);
		memcpy(snapshots.data + offset, profile->snapshots, amount * sizeof(tfx_profile_snapshot_t));
		memset(profile->snapshots, 0, sizeof(tfx_profile_snapshot_t) * tfxPROFILER_SAMPLES);
	}
}

void SetEffectUserData(tfx_particle_manager_t &pm, tfxU32 effect_index, void *data) {
	TFX_ASSERT(effect_index < pm.effects.current_size);    //effect index is out of bounds of the array
	pm.effects[effect_index].user_data = data;
}

void UpdatePMEffect(tfx_particle_manager_t *pm, tfxU32 index, tfxU32 parent_index) {
	tfxPROFILE;

	tfx_effect_state_t &effect = pm->effects[index];

	effect.captured_position = effect.world_position;

	if (pm->lookup_mode == tfxPrecise) {
		effect.frame = effect.age;
	}
	else {
		effect.frame = effect.age / tfxLOOKUP_FREQUENCY;
	}

	UpdateEffectState(pm, index);

	if (IsOrderedEffectState(&effect)) {
		for (tfxEachLayer) {
			tfxU32 &starting_index = effect.instance_data.depth_starting_index[layer];
			tfxU32 current_depth_buffer_index = effect.instance_data.current_depth_buffer_index[layer];
			starting_index = effect.instance_data.depth_indexes[layer][current_depth_buffer_index].current_size;
		}
	}

	tfx_emitter_properties_t &properties = effect.library->emitter_properties[effect.properties_index];

	if (effect.parent_particle_index != tfxINVALID) {
		tfxParticleID parent_particle_id = pm->particle_indexes[effect.parent_particle_index];
		if (parent_particle_id != tfxINVALID) {
			tfxU32 sprite_id = GetParticleSpriteIndex(pm, parent_particle_id);
			tfxU32 sprite_index = sprite_id & 0x0FFFFFFF;
			tfxU32 sprite_layer = (sprite_id & 0xF0000000) >> 28;
			tfx_buffer_t *instance_buffer = pm->flags & tfxParticleManagerFlags_recording_sprites ? &pm->instance_buffer_for_recording[pm->current_sprite_buffer ^ 1][sprite_layer] : &pm->instance_buffer;
			if (sprite_id != tfxINVALID) {
				if (effect.property_flags & tfxEmitterPropertyFlags_effect_is_3d) {
					tfx_3d_instance_t *sprites = tfxCastBuffer(tfx_3d_instance_t, instance_buffer);
					tfx_sprite_transform3d_t transform = { sprites[sprite_index].position.xyz(), sprites[sprite_index].rotations};
					TransformEffector3d(&effect.world_rotations, &effect.local_rotations, &effect.world_position, &effect.local_position, &effect.rotation, &transform, true, effect.property_flags &tfxEmitterPropertyFlags_relative_angle);
				}
				else {
					tfx_2d_instance_t *sprites = tfxCastBuffer(tfx_2d_instance_t, instance_buffer);
					tfx_sprite_transform2d_t transform = { sprites[sprite_index].position.xy(), {0.f, 0.f}, sprites[sprite_index].position.w };
					TransformEffector2d(&effect.world_rotations, &effect.local_rotations, &effect.world_position, &effect.local_position, &effect.rotation, &transform, true, effect.property_flags & tfxEmitterPropertyFlags_relative_angle);
				}

				effect.world_position += properties.emitter_handle * effect.overal_scale;
				if (effect.state_flags & tfxEffectStateFlags_no_tween_this_update || effect.state_flags & tfxEffectStateFlags_no_tween) {
					effect.captured_position = effect.world_position;
				}
			}
		}
		else {
			effect.parent_particle_index = tfxINVALID;
			effect.state_flags |= tfxEffectStateFlags_retain_matrix;
			effect.local_position = effect.world_position;
			effect.local_rotations.roll = effect.world_rotations.roll;
			effect.state_flags |= tfxEffectStateFlags_stop_spawning;
		}
	}
	else {
		if (!(effect.state_flags & tfxEffectStateFlags_retain_matrix)) {
			effect.world_position = effect.local_position + effect.translation;
			effect.world_position += properties.emitter_handle * effect.overal_scale;
			if (effect.property_flags & tfxEmitterPropertyFlags_effect_is_3d) {
				effect.world_rotations = effect.local_rotations;
				effect.rotation = EulerToQuaternion(effect.local_rotations.pitch, effect.local_rotations.yaw, effect.local_rotations.roll);
			}
			else {
				effect.world_rotations.roll = effect.local_rotations.roll;
				ToQuaternion2d(&effect.rotation, effect.local_rotations.roll);
			}
		}

		if (effect.state_flags & tfxEffectStateFlags_no_tween_this_update || effect.state_flags & tfxEffectStateFlags_no_tween) {
			effect.captured_position = effect.world_position;
		}
	}

	effect.age += pm->frame_length;
	effect.highest_particle_age -= pm->frame_length;

	if (properties.loop_length && effect.age > properties.loop_length)
		effect.age -= properties.loop_length;

	if (effect.highest_particle_age <= 0 && effect.age > pm->frame_length * 5.f) {
		effect.timeout_counter += pm->frame_length;
	}
	else {
		effect.timeout_counter = 0;
	}

	if (effect.state_flags & tfxEffectStateFlags_remove) {
		effect.timeout = 1;
		effect.highest_particle_age = 0;
	}

	effect.state_flags &= ~tfxEffectStateFlags_no_tween_this_update;
}

void UpdatePMEmitter(tfx_work_queue_t *work_queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *spawn_work_entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfxU32 emitter_index = spawn_work_entry->emitter_index;

	tfx_particle_manager_t *pm = spawn_work_entry->pm;

	tfx_emitter_state_t &emitter = pm->emitters[emitter_index];

	tfx_effect_state_t &parent_effect = pm->effects[emitter.parent_index];
	emitter.state_flags |= parent_effect.state_flags & tfxEmitterStateFlags_remove;

	tfx_vec3_t local_rotations;
	tfx_library_t *library = pm->library;
	tfx_vec3_t translation;
	translation.x = LookupPrecise(&pm->library->transform_attributes[emitter.transform_attributes].translation_x, emitter.age);
	translation.y = LookupPrecise(&pm->library->transform_attributes[emitter.transform_attributes].translation_y, emitter.age);
	translation.z = LookupPrecise(&pm->library->transform_attributes[emitter.transform_attributes].translation_z, emitter.age);

	emitter.captured_position = emitter.world_position;

	if (pm->lookup_mode == tfxPrecise) {
		emitter.frame = emitter.age;
	}
	else {
		emitter.frame = emitter.age / tfxLOOKUP_FREQUENCY;
	}

	local_rotations.roll = LookupPrecise(&pm->library->transform_attributes[emitter.transform_attributes].roll, emitter.age);
	local_rotations.pitch = LookupPrecise(&pm->library->transform_attributes[emitter.transform_attributes].pitch, emitter.age);
	local_rotations.yaw = LookupPrecise(&pm->library->transform_attributes[emitter.transform_attributes].yaw, emitter.age);

	spawn_work_entry->highest_particle_age = emitter.highest_particle_age;

	tfx_emitter_properties_t &properties = *spawn_work_entry->properties;

	TFX_ASSERT(emitter.parent_index != tfxINVALID);    //Emitter must have a valid parent (an effect)

	tfxU32 layer = properties.layer;

	parent_effect.timeout_counter = 0;
	if (parent_effect.age < emitter.delay_spawning) {
		return;
	}
	spawn_work_entry->root_effect_flags = pm->effects[emitter.root_index].effect_flags;
	bool ordered_effect = (spawn_work_entry->root_effect_flags & tfxEffectPropertyFlags_age_order) || (spawn_work_entry->root_effect_flags & tfxEffectPropertyFlags_depth_draw_order) > 0;
	emitter.delay_spawning = -pm->frame_length;

	emitter.state_flags |= parent_effect.state_flags & tfxEmitterStateFlags_no_tween;
	if (properties.emission_type != tfxOtherEmitter) {
		emitter.state_flags |= parent_effect.state_flags & tfxEmitterStateFlags_stop_spawning;
	}
	spawn_work_entry->parent_spawn_controls = &parent_effect.spawn_controls;
	spawn_work_entry->parent_property_flags = parent_effect.property_flags;
	spawn_work_entry->parent_index = emitter.parent_index;
	spawn_work_entry->overal_scale = parent_effect.overal_scale;
	UpdateEmitterState(pm, emitter, emitter.parent_index, spawn_work_entry->parent_spawn_controls, spawn_work_entry);

	if ((emitter.state_flags & tfxEmitterStateFlags_is_line_loop_or_kill && emitter.state_flags & tfxEmitterStateFlags_loop) || (emitter.state_flags & tfxEmitterStateFlags_is_single && properties.single_shot_limit > 1)) {
		pm->emitters_check_capture.push_back(emitter_index);
	}

	//bool is_compute = emitter.property_flags & tfxEmitterPropertyFlags_is_bottom_emitter && pm->flags & tfxParticleManagerFlags_use_compute_shader;
	tfxU32 amount_spawned = 0;
	tfxU32 max_spawn_count = NewSpritesNeeded(pm, &spawn_work_entry->random, emitter_index, &parent_effect, &properties);
	tfx_effect_instance_data_t &instance_data = pm->effects[emitter.root_index].instance_data;

	tfx_vector_t<tfx_unique_sprite_id_t> &uid_buffer = pm->unique_sprite_ids[pm->current_sprite_buffer][layer];
	bool is_recording = (pm->flags & tfxParticleManagerFlags_recording_sprites) > 0 && (pm->flags & tfxParticleManagerFlags_using_uids) > 0;
	tfx_buffer_t &instance_buffer = !is_recording ? pm->instance_buffer : pm->instance_buffer_for_recording[pm->current_sprite_buffer][layer];
	tfxU32 &effect_instance_index_point = instance_data.sprite_index_point[layer];
	tfxU32 &pm_instance_index_point =  pm->sprite_index_point[layer];
	tfxU32 &instance_index_point = (pm->flags & tfxParticleManagerFlags_use_effect_sprite_buffers) ? instance_data.sprite_index_point[layer] : pm->sprite_index_point[layer];
	if (ordered_effect) {
		spawn_work_entry->depth_indexes = &instance_data.depth_indexes[layer][instance_data.current_depth_buffer_index[layer]];
	}

	if (emitter.property_flags & tfxEmitterPropertyFlags_effect_is_3d) {
		tfx__transform_3d(&emitter.world_rotations, &local_rotations, &spawn_work_entry->overal_scale, &emitter.world_position, &emitter.local_position, &translation, &emitter.rotation, &parent_effect);
	}
	else {
		tfx__transform_2d(&emitter.world_rotations, &local_rotations, &spawn_work_entry->overal_scale, &emitter.world_position, &emitter.local_position, &translation, &emitter.rotation, &parent_effect);
	}

	if (emitter.state_flags & tfxEmitterStateFlags_no_tween_this_update || emitter.state_flags & tfxEmitterStateFlags_no_tween) {
		emitter.captured_position = emitter.world_position;
	}

	emitter.sprites_count = pm->particle_array_buffers[emitter.particles_index].current_size;
	if (pm->flags & tfxParticleManagerFlags_dynamic_sprite_allocation || pm->flags & tfxParticleManagerFlags_use_effect_sprite_buffers) {
		if (emitter.sprites_count + instance_buffer.current_size + max_spawn_count >= instance_buffer.capacity) {
			tfxU32 new_size = instance_buffer.capacity + (emitter.sprites_count + max_spawn_count) + 1;
			if (pm->flags & tfxParticleManagerFlags_recording_sprites && pm->flags & tfxParticleManagerFlags_using_uids) {
				uid_buffer.reserve(new_size);
			}
			if (pm->flags & tfxParticleManagerFlags_direct_to_staging_buffer && pm->info.grow_staging_buffer_callback) {
				if (!pm->info.grow_staging_buffer_callback(new_size, pm, pm->info.user_data)) {
					max_spawn_count = 0;
				}
			} else if(pm->flags & tfxParticleManagerFlags_direct_to_staging_buffer) {
				max_spawn_count = 0;
			} else {
				instance_buffer.reserve(new_size);
			}
		}
	}
	else {
		tfxU32 required_space = emitter.sprites_count + instance_buffer.current_size + max_spawn_count;
		if (required_space >= instance_buffer.capacity) {
			tfxU32 free_space = instance_buffer.capacity - instance_buffer.current_size;
			if (free_space > emitter.sprites_count) {
				max_spawn_count = free_space - emitter.sprites_count;
			}
			else {
				max_spawn_count = tfxMin(free_space, max_spawn_count);
			}
			TFX_ASSERT(free_space >= max_spawn_count);    //Trying to spawn particles when no space left in sprite buffer. If this is hit then there's a bug in TimelineFX!
		}
	}

	instance_buffer.current_size += max_spawn_count + emitter.sprites_count;
	emitter.sprites_count += max_spawn_count;
	emitter.sprites_index = instance_index_point;
	effect_instance_index_point += emitter.sprites_count;
	pm_instance_index_point += emitter.sprites_count;

	spawn_work_entry->max_spawn_count = max_spawn_count;

	if (emitter.state_flags & tfxEmitterStateFlags_is_sub_emitter) {
		if (emitter.age > 0 && !(pm->flags & tfxParticleManagerFlags_disable_spawning)) {
			amount_spawned = SpawnParticles(pm, spawn_work_entry);
		}
	}
	else {
		if (!(pm->flags & tfxParticleManagerFlags_disable_spawning)) {
			amount_spawned = SpawnParticles(pm, spawn_work_entry);
		}
	}

	TFX_ASSERT(amount_spawned <= max_spawn_count);
	tfxU32 spawn_difference = max_spawn_count - amount_spawned;
	instance_buffer.current_size -= spawn_difference;
	TFX_ASSERT(instance_buffer.current_size < instance_buffer.capacity);
	pm_instance_index_point -= spawn_difference;
	effect_instance_index_point -= spawn_difference;
	pm->layer_sizes[layer] += emitter.sprites_count - spawn_difference;
	instance_data.instance_count += emitter.sprites_count - spawn_difference;

	if (pm->flags & tfxParticleManagerFlags_recording_sprites && pm->flags & tfxParticleManagerFlags_using_uids) {
		uid_buffer.current_size = instance_buffer.current_size;
	}

	emitter.age += pm->frame_length;
	if (!(emitter.property_flags & tfxEmitterPropertyFlags_single) || (emitter.property_flags & tfxEmitterPropertyFlags_single && properties.single_shot_limit > 0) || emitter.state_flags & tfxEmitterStateFlags_stop_spawning) {
		emitter.highest_particle_age -= pm->frame_length;
		spawn_work_entry->highest_particle_age -= pm->frame_length;
	}

	if (properties.loop_length && emitter.age > properties.loop_length)
		emitter.age -= properties.loop_length;

	if (emitter.highest_particle_age <= 0 && emitter.age > pm->frame_length * 5.f) {
		emitter.timeout_counter += pm->frame_length;
	}
	else {
		emitter.timeout_counter = 0;
	}

	if (emitter.state_flags & tfxEmitterStateFlags_remove) {
		emitter.timeout = 0;
		emitter.highest_particle_age = 0;
		return;
	}

	emitter.state_flags &= ~tfxEmitterStateFlags_no_tween_this_update;
}

tfxU32 NewSpritesNeeded(tfx_particle_manager_t *pm, tfx_random_t *random, tfxU32 index, tfx_effect_state_t *parent, tfx_emitter_properties_t *properties) {
	tfx_emitter_state_t &emitter = pm->emitters[index];
	AlterRandomSeed(random, 25 + emitter.seed_index);
	if (!(emitter.property_flags & tfxEmitterPropertyFlags_single)) {
		emitter.spawn_quantity = lookup_callback(&pm->library->emitter_attributes[emitter.emitter_attributes].base.amount, emitter.frame);
		float amount_variation = lookup_callback(&pm->library->emitter_attributes[emitter.emitter_attributes].variation.amount, emitter.frame);
		emitter.spawn_quantity += amount_variation > 0.f ? RandomRange(random, 1.f, amount_variation) : 0.f;
		emitter.spawn_quantity *= lookup_callback(&pm->library->global_graphs[parent->global_attributes].amount, emitter.frame);
	}
	else {
		emitter.spawn_quantity = (float)properties->spawn_amount + RandomRange(random, (float)properties->spawn_amount_variation);
		emitter.spawn_quantity *= lookup_callback(&pm->library->global_graphs[parent->global_attributes].amount, emitter.frame);
	}

	if (emitter.state_flags & tfxEmitterStateFlags_single_shot_done || emitter.state_flags & tfxEmitterStateFlags_stop_spawning) {
		return 0;
	}

	if (emitter.spawn_quantity == 0) {
		return 0;
	}

	if (properties->emission_type == tfxPath && emitter.state_flags & tfxEmitterStateFlags_has_rotated_path) {
		tfx_emitter_path_t *path = &pm->library->paths[emitter.path_attributes];
		emitter.spawn_quantity *= (float)emitter.active_paths / (float)path->maximum_active_paths;
		emitter.path_stagger_counter += pm->frame_length;
	}

	if (emitter.spawn_locations_index != tfxINVALID && properties->emission_type == tfxOtherEmitter && emitter.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio) {
		tfx_soa_buffer_t &spawn_point_buffer = pm->particle_array_buffers[emitter.spawn_locations_index];
		emitter.spawn_quantity *= spawn_point_buffer.current_size;
	}

	float step_size = 0;
	if (!(emitter.property_flags & tfxEmitterPropertyFlags_single)) {
		if (emitter.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && (properties->emission_type == tfxArea || properties->emission_type == tfxEllipse)) {
			if (emitter.property_flags & tfxEmitterPropertyFlags_effect_is_3d) {
				float area = tfxMax(0.1f, emitter.emitter_size.x) * tfxMax(0.1f, emitter.emitter_size.y) * tfxMax(0.1f, emitter.emitter_size.z);
				emitter.spawn_quantity = (emitter.spawn_quantity / 50.f) * area;
			}
			else {
				float area = emitter.emitter_size.x * emitter.emitter_size.y;
				emitter.spawn_quantity = (emitter.spawn_quantity / 10000.f) * area;
			}
		}
		else if (emitter.property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && properties->emission_type == tfxLine) {
			emitter.spawn_quantity = (emitter.spawn_quantity / 100.f) * emitter.emitter_size.y;
		}

		emitter.spawn_quantity *= pm->update_time;
		step_size = 1.f / emitter.spawn_quantity;
	}
	else if (emitter.property_flags & tfxEmitterPropertyFlags_match_amount_to_grid_points) {
		if (emitter.property_flags & tfxEmitterPropertyFlags_match_amount_to_grid_points && emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {
			float x = tfxMax(properties->grid_points.x, 1.f);
			float y = tfxMax(properties->grid_points.y, 1.f);
			float z = tfxMax(properties->grid_points.z, 1.f);
			if (properties->emission_type == tfxArea) {
			}
			else if (properties->emission_type == tfxCylinder) {
				emitter.spawn_quantity = x * y;
			}
			switch (properties->emission_type) {
			case tfx_emission_type::tfxArea:
				if (emitter.property_flags & tfxEmitterPropertyFlags_effect_is_3d) {
					if (emitter.property_flags & tfxEmitterPropertyFlags_fill_area) {
						emitter.spawn_quantity = x * y * z;
					}
					else if (emitter.property_flags & tfxEmitterPropertyFlags_area_open_ends) {
						emitter.spawn_quantity = x * z * 2 + y * z * 2 - 4 * z;
					}
					else {
						emitter.spawn_quantity = x * z * 2 + (y - 2) * (x * 2 + z * 2 - 4);
					}
				}
				else {
					if (emitter.property_flags & tfxEmitterPropertyFlags_fill_area) {
						emitter.spawn_quantity = x * y;
					}
					else {
						emitter.spawn_quantity = x * 2 + (y - 2) * 2;
					}
				}
				break;
			case tfx_emission_type::tfxCylinder:
				emitter.spawn_quantity = x * y;
				break;
			case tfx_emission_type::tfxEllipse:
				if (!(emitter.property_flags & tfxEmitterPropertyFlags_effect_is_3d)) {
					emitter.spawn_quantity = x;
				}
				break;
			case tfx_emission_type::tfxLine:
				emitter.spawn_quantity = x;
				break;
			case tfx_emission_type::tfxIcosphere:
				emitter.spawn_quantity = (float)tfxIcospherePoints[tfxMin((tfxU32)x, 5)].current_size;
				break;
			default:
				break;
			}
		}
		step_size = 1.f / emitter.spawn_quantity;
	}
	else {
		step_size = 1.f / emitter.spawn_quantity;
	}

	float tween = emitter.amount_remainder;
	return tween >= 1.f ? 0 : tfxU32((1.f - emitter.amount_remainder) / step_size) + 1;
}

void CompletePMWork(tfx_particle_manager_t *pm) {
	tfxCompleteAllWork(&pm->work_queue);
}

tfxU32 SpawnParticles(tfx_particle_manager_t *pm, tfx_spawn_work_entry_t *work_entry) {
	tfx_emitter_state_t &emitter = pm->emitters[work_entry->emitter_index];
	const tfx_emitter_properties_t &properties = *work_entry->properties;
	const tfxU32 layer = properties.layer;
	// Note that adding this in will mess up other emitter emission types and looped sprite recording || pm->effects[work_entry->parent_index].state_flags & tfxEffectStateFlags_stop_spawning
	if (emitter.state_flags & tfxEmitterStateFlags_single_shot_done || emitter.state_flags & tfxEmitterStateFlags_stop_spawning)
		return 0;
	if (emitter.spawn_quantity == 0)
		return 0;

	float step_size = 1.f / emitter.spawn_quantity;
	float tween = 0;
	if (step_size == emitter.qty_step_size || emitter.property_flags & tfxEmitterPropertyFlags_single) {
		tween = emitter.amount_remainder;
	}
	else {
		tween = emitter.amount_remainder - (emitter.qty_step_size - step_size);
	}
	emitter.qty_step_size = step_size;
	//bool is_compute = work_entry->e->property_flags & tfxEmitterPropertyFlags_is_bottom_emitter && pm->emitter.state_flags & tfxParticleManagerFlags_use_compute_shader;

	if (tween >= 1) {
		tween -= emitter.spawn_quantity;
	}

	work_entry->tween = tween;
	work_entry->qty_step_size = step_size;
	work_entry->amount_to_spawn = 0;
	work_entry->particle_data = &pm->particle_arrays[emitter.particles_index];
	work_entry->property_flags = emitter.property_flags;

	if (tween >= 1) {
		emitter.amount_remainder = tween - 1.f;
	}
	else if (!(emitter.property_flags & tfxEmitterPropertyFlags_single)) {
		float amount_that_will_spawn = (1.f - tween) / step_size;
		work_entry->amount_to_spawn = (tfxU32)ceilf(amount_that_will_spawn);
		if (work_entry->amount_to_spawn > work_entry->max_spawn_count) {
			work_entry->amount_to_spawn = work_entry->max_spawn_count;
			emitter.amount_remainder = 0.f;
		}
		else {
			emitter.amount_remainder = amount_that_will_spawn - (tfxU32)amount_that_will_spawn;
			emitter.amount_remainder = (1.f - emitter.amount_remainder) * step_size;
		}
	}
	else {
		work_entry->amount_to_spawn = (tfxU32)emitter.spawn_quantity;
	}

	work_entry->emission_type = properties.emission_type;
	if (work_entry->emission_type == tfxOtherEmitter) {
		if (emitter.spawn_locations_index == tfxINVALID) {
			work_entry->amount_to_spawn = 0;
		}
		else {
			tfx_soa_buffer_t &spawn_point_buffer = pm->particle_location_buffers[emitter.spawn_locations_index];
			if (spawn_point_buffer.current_size == 0) {
				work_entry->amount_to_spawn = 0;
			}
		}
	}

	tfx_soa_buffer_t &buffer = pm->particle_array_buffers[emitter.particles_index];

	bool grew = false;
	tfxU32 start_index = buffer.start_index;
	work_entry->spawn_start_index = AddRows(&pm->particle_array_buffers[emitter.particles_index], work_entry->amount_to_spawn, true, grew);
	if (grew && work_entry->depth_indexes && start_index > 0) {
		//Todo: This should be avoided by allocating the correct amount for the particle buffer ahead of time
		//If the particle buffer is allocated a larger memory size then the ring buffer index has to be reset in the depth buffer list
		//to align them to the correct particle ids again.
		tfx_particle_soa_t &bank = pm->particle_arrays[emitter.particles_index];
		//TFX_ASSERT(pm->particle_array_buffers[particles_index].start_index == 0); //If start_index isn't 0 after the arrays grew then something went wrong with the allocation
		for (int i = 0; i != pm->particle_array_buffers[emitter.particles_index].current_size - work_entry->amount_to_spawn; ++i) {
			(*work_entry->depth_indexes)[bank.depth_index[i]].particle_id = tfx__make_particle_id(emitter.particles_index, i);
		}
	}

	work_entry->depth_index_start = work_entry->depth_indexes ? work_entry->depth_indexes->current_size : 0;

	if (work_entry->amount_to_spawn > 0) {
		if (emitter.state_flags & tfxEmitterStateFlags_is_in_ordered_effect) {
			TFX_ASSERT(work_entry->depth_indexes);	//This must have been set if the effect is ordered!
			if (work_entry->depth_indexes && work_entry->depth_indexes->capacity <= work_entry->depth_indexes->current_size + work_entry->amount_to_spawn) {
				tfxU32 new_capacity = work_entry->depth_indexes->current_size + work_entry->amount_to_spawn;
				work_entry->depth_indexes->reserve(tfx__Max(new_capacity + new_capacity / 2, 8));
			}
			work_entry->depth_indexes->current_size += work_entry->amount_to_spawn;
			TFX_ASSERT(work_entry->depth_indexes->current_size < work_entry->depth_indexes->capacity);
			pm->deffered_spawn_work.push_back(work_entry);
		}
		else {
			tfxAddWorkQueueEntry(&pm->work_queue, work_entry, pm->flags & tfxParticleManagerFlags_3d_effects ? DoSpawnWork3d : DoSpawnWork2d);
		}
	}

	if (work_entry->amount_to_spawn > 0 && emitter.property_flags & tfxEmitterPropertyFlags_single) {
		emitter.state_flags |= tfxEmitterStateFlags_single_shot_done;
	}

	return work_entry->amount_to_spawn;
}

void DoSpawnWork3d(tfx_work_queue_t *queue, void *data) {
	tfx_spawn_work_entry_t *work_entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_particle_manager_t *pm = work_entry->pm;
	tfx_emitter_state_t &emitter = pm->emitters[work_entry->emitter_index];
	SpawnParticleAge(&pm->work_queue, work_entry);
	if (work_entry->emission_type == tfxOtherEmitter) {
		SpawnParticleOtherEmitter3d(&pm->work_queue, work_entry);
	}
	else if (work_entry->emission_type == tfxPoint) {
		SpawnParticlePoint3d(&pm->work_queue, work_entry);
	}
	else if (work_entry->emission_type == tfxArea) {
		SpawnParticleArea3d(&pm->work_queue, work_entry);
	}
	else if (work_entry->emission_type == tfxEllipse) {
		SpawnParticleEllipsoid(&pm->work_queue, work_entry);
	}
	else if (work_entry->emission_type == tfxLine) {
		SpawnParticleLine3d(&pm->work_queue, work_entry);
	}
	else if (work_entry->emission_type == tfxIcosphere && work_entry->property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
		SpawnParticleIcosphereRandom3d(&pm->work_queue, work_entry);
	}
	else if (work_entry->emission_type == tfxIcosphere && !(work_entry->property_flags & tfxEmitterPropertyFlags_grid_spawn_random)) {
		SpawnParticleIcosphere3d(&pm->work_queue, work_entry);
	}
	else if (work_entry->emission_type == tfxCylinder) {
		SpawnParticleCylinder3d(&pm->work_queue, work_entry);
	}
	else if (work_entry->emission_type == tfxPath) {
		SpawnParticlePath3d(&pm->work_queue, work_entry);
	}

	SpawnParticleWeight(&pm->work_queue, work_entry);
	SpawnParticleVelocity(&pm->work_queue, work_entry);
	SpawnParticleRoll(&pm->work_queue, work_entry);
	SpawnParticleMicroUpdate3d(&pm->work_queue, work_entry);
	if (emitter.control_profile & tfxEmitterControlProfile_noise) {
		SpawnParticleNoise(&pm->work_queue, work_entry);
	}
	else if (emitter.control_profile & tfxEmitterControlProfile_motion_randomness) {
		SpawnParticleMotionRandomness(&pm->work_queue, work_entry);
	}
	SpawnParticleImageFrame(&pm->work_queue, work_entry);
	SpawnParticleSize3d(&pm->work_queue, work_entry);
	SpawnParticleSpin3d(&pm->work_queue, work_entry);
}

void DoSpawnWork2d(tfx_work_queue_t *queue, void *data) {
	tfx_spawn_work_entry_t *work_entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_particle_manager_t *pm = work_entry->pm;
	tfx_emitter_state_t &emitter = pm->emitters[work_entry->emitter_index];
	SpawnParticleAge(&pm->work_queue, work_entry);
	if (work_entry->emission_type == tfxOtherEmitter) {
		SpawnParticleOtherEmitter2d(&pm->work_queue, work_entry);
	}
	else if (work_entry->emission_type == tfxPoint) {
		SpawnParticlePoint2d(&pm->work_queue, work_entry);
	}
	else if (work_entry->emission_type == tfxArea) {
		SpawnParticleArea2d(&pm->work_queue, work_entry);
	}
	else if (work_entry->emission_type == tfxEllipse) {
		SpawnParticleEllipse2d(&pm->work_queue, work_entry);
	}
	else if (work_entry->emission_type == tfxLine) {
		SpawnParticleLine2d(&pm->work_queue, work_entry);
	}
	else if (work_entry->emission_type == tfxPath) {
		SpawnParticlePath2d(&pm->work_queue, work_entry);
	}
	SpawnParticleWeight(&pm->work_queue, work_entry);
	SpawnParticleVelocity(&pm->work_queue, work_entry);
	SpawnParticleRoll(&pm->work_queue, work_entry);
	SpawnParticleMicroUpdate2d(&pm->work_queue, work_entry);
	if (emitter.control_profile & tfxEmitterControlProfile_noise) {
		SpawnParticleNoise(&pm->work_queue, work_entry);
	}
	else if (emitter.control_profile & tfxEmitterControlProfile_motion_randomness) {
		SpawnParticleMotionRandomness(&pm->work_queue, work_entry);
	}
	SpawnParticleImageFrame(&pm->work_queue, work_entry);
	SpawnParticleSize2d(&pm->work_queue, work_entry);
	SpawnParticleSpin2d(&pm->work_queue, work_entry);
}

void SpawnParticleAge(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	const tfx_emitter_properties_t &properties = *entry->properties;
	float tween = entry->tween;
	tfxU32 emitter_index = entry->emitter_index;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 2 + emitter.seed_index);

	//const float life = pm.emitters.life[emitter_index];
	//const float life_variation = pm.emitters.life_variation[emitter_index];
	tfx_library_t *library = pm.library;
	const float life = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].base.life, emitter.frame) * entry->parent_spawn_controls->life;
	const float life_variation = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].variation.life, emitter.frame) * entry->parent_spawn_controls->life;
	const tfxU32 loop_count = entry->properties->single_shot_limit + 1;
	const float emitter_intensity = entry->parent_spawn_controls->intensity;
	const float first_red_value = GetGraphFirstValue(&library->emitter_attributes[emitter.emitter_attributes].overtime.red);
	const float first_green_value = GetGraphFirstValue(&library->emitter_attributes[emitter.emitter_attributes].overtime.green);
	const float first_blue_value = GetGraphFirstValue(&library->emitter_attributes[emitter.emitter_attributes].overtime.blue);
	const float first_alpha_value = GetGraphFirstValue(&library->emitter_attributes[emitter.emitter_attributes].overtime.blendfactor);
	const tfx_color_ramp_t &color_ramp = library->emitter_attributes[emitter.emitter_attributes].overtime.color_ramps[0];

	TFX_ASSERT(random.seeds[0] > 0);

	float frame_fraction = entry->pm->frame_length / (float)entry->amount_to_spawn;
	float age_accumulator = 0.f;

	if (entry->emission_type == tfxOtherEmitter) {
		if ((tfxU32)emitter.grid_coords.x >= pm.particle_location_buffers[emitter.spawn_locations_index].current_size) {
			emitter.grid_coords.x = 0.f;
		}
		emitter.grid_coords.y = emitter.grid_coords.x;
	}

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &age = entry->particle_data->age[index];
		float &max_age = entry->particle_data->max_age[index];
		tfx_rgba8_t &color = entry->particle_data->color[index];
		tfxU32 &single_loop_count = entry->particle_data->single_loop_count[index];

		tfxParticleFlags &flags = entry->particle_data->flags[index];
		tfxU32 &particle_index = entry->particle_data->particle_index[index];
		particle_index = tfxINVALID;

		flags = tfxParticleFlags_capture_after_transform;

		//Set the age of the particle to an interpolated value between 0 and the length of the frame depending on how many particles are being spawned this frame
		age = age_accumulator;
		//I don't have any better idea than this for a unique id generation. Because of the way the simple random movement works it uses the uid of the particle
		//as a random seed to so that we don't have to store any extra values related to the random movement. This means that the uid needs to be staggered enough
		//so that particles don't share seeds. This might sound strange because they are UIDs but there's a few random numbers generated for each particle and we offset from the uid
		//to get different numbers. When the UID is incrementle then future particles will just get the same seed unless the UID can be staggered enough. Hence I 
		//just use clock cycles which serves the purpose well enough. It's kind of hard to explain but see more in ControlParticleSimpleRandomMovement
		entry->particle_data->uid[index] = (tfxU32)tfx__rdtsc();
		age_accumulator += frame_fraction;
		if (emitter.property_flags & tfxEmitterPropertyFlags_wrap_single_sprite && pm.flags & tfxParticleManagerFlags_recording_sprites) {
			max_age = tfx__Max(pm.animation_length_in_time, 1.f);
		}
		else {
			max_age = tfx__Max(life + RandomRange(&random, life_variation), 1.f);
		}
		single_loop_count = 0;

		float alpha = 255.f * first_alpha_value;

		if (entry->emission_type == tfxOtherEmitter) {
			int spawn_index = (int)emitter.grid_coords.x;
			spawn_index = GetCircularIndex(&pm.particle_location_buffers[emitter.spawn_locations_index], spawn_index);
			float age_lerp = pm.particle_location_arrays[emitter.spawn_locations_index].age[spawn_index];
			emitter.grid_coords.x++;
			if ((tfxU32)emitter.grid_coords.x >= pm.particle_location_buffers[emitter.spawn_locations_index].current_size) {
				emitter.grid_coords.x = 0.f;
			}
			float factor = lookup_overtime_callback(&library->emitter_attributes[emitter.emitter_attributes].factor.life, age_lerp * max_age, max_age);
			max_age *= factor;
			max_age = tfx__Max(max_age, pm.frame_length);
			factor = lookup_overtime_callback(&library->emitter_attributes[emitter.emitter_attributes].factor.intensity, age_lerp * max_age, max_age);
			entry->particle_data->intensity_factor[index] = factor;
		}
		else {
			entry->particle_data->intensity_factor[index] = 1.f;
		}

		if (emitter.state_flags & tfxEmitterStateFlags_random_color) {
			float age = RandomRange(&random, max_age) / max_age * (tfxCOLOR_RAMP_WIDTH - 1.f);
			color = color_ramp.colors[(int)age];
			color.a = (tfxU32)alpha;
		}
		else {
			color = color_ramp.colors[0];
			color.a = (tfxU32)alpha;
		}

		entry->highest_particle_age = fmaxf(entry->highest_particle_age, (max_age * loop_count) + pm.frame_length + 1);

		if (entry->sub_effects->current_size > 0) {
			particle_index = GetPMParticleIndexSlot(&pm, tfx__make_particle_id(emitter.particles_index, index));
			flags |= tfxParticleFlags_has_sub_effects;
			for (auto &sub : *entry->sub_effects) {
				if (!FreePMEffectCapacity(&pm))
					break;
				TFX_ASSERT(entry->depth < tfxMAXDEPTH - 1);
				tfxU32 added_index = AddEffectToParticleManager(&pm, &sub, pm.current_ebuff, entry->depth + 1, true, emitter.root_index, 0.f);
				pm.effects[added_index].overal_scale = entry->overal_scale;
				pm.effects[added_index].parent_particle_index = particle_index;
			}
		}
	}
}

void SpawnParticleImageFrame(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 3 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;

	for (int i = 0; i != entry->amount_to_spawn; ++i) {

		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &image_frame = entry->particle_data->image_frame[index];
		tfxU32 &sprites_index = entry->particle_data->sprite_index[index];
		sprites_index = tfxINVALID;

		//----Image
		//data.image = GetEffectProperties(this)->image;
		if (emitter.property_flags & tfxEmitterPropertyFlags_random_start_frame && properties.image->animation_frames > 1) {
			image_frame = RandomRange(&random, properties.image->animation_frames);
		}
		else {
			image_frame = properties.start_frame;
		}

	}

}

void SpawnParticleSize2d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_library_t *library = pm.library;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 4 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;

	tfx_vec2_t size;
	if (!(emitter.property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
		size.x = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].base.width, emitter.frame) * entry->parent_spawn_controls->size_x;
		size.y = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].base.height, emitter.frame) * entry->parent_spawn_controls->size_y;
	}
	else {
		size.x = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].base.width, emitter.frame);
		if (entry->parent_property_flags & tfxEffectPropertyFlags_global_uniform_size)
			size.y = size.x * entry->parent_spawn_controls->size_x;
		else
			size.y = size.x * entry->parent_spawn_controls->size_y;
		size.x *= entry->parent_spawn_controls->size_x;
	}
	tfx_vec2_t size_variation;
	size_variation.x = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].variation.width, emitter.frame) * entry->parent_spawn_controls->size_x;
	size_variation.y = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].variation.height, emitter.frame) * entry->parent_spawn_controls->size_y;

	if (entry->emission_type == tfxOtherEmitter) {
		emitter.grid_coords.x = emitter.grid_coords.y;
	}

	for (int i = 0; i != entry->amount_to_spawn; ++i) {

		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &base_size_x = entry->particle_data->base_size_x[index];
		float &base_size_y = entry->particle_data->base_size_y[index];

		//----Size
		if (!(emitter.property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			float random_size_x = RandomRange(&random, size_variation.x);
			float random_size_y = RandomRange(&random, size_variation.y);
			base_size_y = (random_size_y + size.y);
			base_size_x = (random_size_x + size.x);
		}
		else {
			float random_size_x = RandomRange(&random, size_variation.x);
			float random_size_y = random_size_x;
			base_size_y = (random_size_y + size.y);
			base_size_x = (random_size_x + size.x);
		}

		if (entry->emission_type == tfxOtherEmitter) {
			int spawn_index = (int)emitter.grid_coords.x;
			spawn_index = GetCircularIndex(&pm.particle_location_buffers[emitter.spawn_locations_index], spawn_index);
			float age_lerp = pm.particle_location_arrays[emitter.spawn_locations_index].age[spawn_index];
			emitter.grid_coords.x++;
			if ((tfxU32)emitter.grid_coords.x >= pm.particle_location_buffers[emitter.spawn_locations_index].current_size) {
				emitter.grid_coords.x = 0.f;
			}
			float max_age = entry->particle_data->max_age[index];
			float factor = lookup_overtime_callback(&library->emitter_attributes[emitter.emitter_attributes].factor.size, age_lerp * max_age, max_age);
			base_size_x *= factor;
			base_size_y *= factor;
		}
	}

}

void SpawnParticleSize3d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 5 + emitter.seed_index);
	tfx_library_t *library = pm.library;
	const tfx_emitter_properties_t &properties = *entry->properties;

	tfx_vec2_t size;
	if (!(emitter.property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
		size.x = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].base.width, emitter.frame) * entry->parent_spawn_controls->size_x;
		size.y = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].base.height, emitter.frame) * entry->parent_spawn_controls->size_y;
	}
	else {
		size.x = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].base.width, emitter.frame);
		if (entry->parent_property_flags & tfxEffectPropertyFlags_global_uniform_size)
			size.y = size.x * entry->parent_spawn_controls->size_x;
		else
			size.y = size.x * entry->parent_spawn_controls->size_y;
		size.x *= entry->parent_spawn_controls->size_x;
	}
	tfx_vec2_t size_variation;
	size_variation.x = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].variation.width, emitter.frame) * entry->parent_spawn_controls->size_x;
	size_variation.y = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].variation.height, emitter.frame) * entry->parent_spawn_controls->size_y;
	tfx_vec2_t &image_size = properties.image->image_size;

	if (entry->emission_type == tfxOtherEmitter) {
		emitter.grid_coords.x = emitter.grid_coords.y;
	}

	for (int i = 0; i != entry->amount_to_spawn; ++i) {

		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &base_size_x = entry->particle_data->base_size_x[index];
		float &base_size_y = entry->particle_data->base_size_y[index];

		//----Size
		if (!(emitter.property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			float random_size_x = RandomRange(&random, size_variation.x);
			float random_size_y = RandomRange(&random, size_variation.y);
			base_size_y = random_size_y + size.y;
			base_size_x = random_size_x + size.x;
		}
		else {
			float random_size_x = RandomRange(&random, size_variation.x);
			float random_size_y = random_size_x;
			base_size_y = random_size_y + size.y;
			base_size_x = random_size_x + size.x;
		}

		if (entry->emission_type == tfxOtherEmitter) {
			int spawn_index = (int)emitter.grid_coords.x;
			spawn_index = GetCircularIndex(&pm.particle_location_buffers[emitter.spawn_locations_index], spawn_index);
			float age_lerp = pm.particle_location_arrays[emitter.spawn_locations_index].age[spawn_index];
			emitter.grid_coords.x++;
			if ((tfxU32)emitter.grid_coords.x >= pm.particle_location_buffers[emitter.spawn_locations_index].current_size) {
				emitter.grid_coords.x = 0.f;
			}
			float max_age = entry->particle_data->max_age[index];
			float factor = lookup_overtime_callback(&library->emitter_attributes[emitter.emitter_attributes].factor.size, age_lerp * max_age, max_age);
			base_size_x *= factor;
			base_size_y *= factor;
		}

	}

}

void SpawnParticleNoise(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfxU32 emitter_index = entry->emitter_index;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_library_t *library = pm.library;
	tfx_emitter_state_t &emitter = pm.emitters[emitter_index];
	AlterRandomSeed(&random, 6 + emitter.seed_index);
	float emitter_noise_offset_variation = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].variation.noise_offset, emitter.frame);
	float emitter_noise_offset = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].base.noise_offset, emitter.frame);
	float emitter_noise_resolution = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].variation.noise_resolution, emitter.frame);
	const float parent_noise_base_offset = pm.effects[emitter.parent_index].noise_base_offset;

	for (int i = 0; i != entry->amount_to_spawn; ++i) {

		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &noise_offset = entry->particle_data->noise_offset[index];
		float &noise_resolution = entry->particle_data->noise_resolution[index];

		//----Motion randomness
		noise_offset = RandomRange(&random, emitter_noise_offset_variation) + emitter_noise_offset + parent_noise_base_offset;
		noise_resolution = emitter_noise_resolution + 0.01f;

	}
}

void SpawnParticleMotionRandomness(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfxU32 emitter_index = entry->emitter_index;
	tfx_emitter_state_t &emitter = entry->pm->emitters[emitter_index];

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&entry->pm->particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		entry->particle_data->noise_offset[index] = 0.f;
		entry->particle_data->noise_resolution[index] = 0.f;
	}
}

void SpawnParticleSpin2d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 7 + emitter.seed_index);

	const float spin = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].base.spin, emitter.frame) * entry->parent_spawn_controls->spin;
	const float spin_variation = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].variation.spin, emitter.frame) * entry->parent_spawn_controls->spin;

	for (int i = 0; i != entry->amount_to_spawn; ++i) {

		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &base_spin = entry->particle_data->base_spin[index];

		//----Spin
		base_spin = RandomRange(&random, -spin_variation, spin_variation) + spin;

	}

}

void SpawnParticleSpin3d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 8 + emitter.seed_index);

	const float roll_spin = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].base.spin, emitter.frame) * entry->parent_spawn_controls->spin;
	const float pitch_spin = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].base.pitch_spin, emitter.frame) * entry->parent_spawn_controls->pitch_spin;
	const float yaw_spin = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].base.yaw_spin, emitter.frame) * entry->parent_spawn_controls->yaw_spin;
	const float spin_variation = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].variation.spin, emitter.frame) * entry->parent_spawn_controls->spin;
	const float spin_pitch_variation = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].variation.pitch_spin, emitter.frame) * entry->parent_spawn_controls->pitch_spin;
	const float spin_yaw_variation = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].variation.yaw_spin, emitter.frame) * entry->parent_spawn_controls->yaw_spin;
	const tfxAngleSettingFlags angle_settings = entry->properties->angle_settings;

	for (int i = 0; i != entry->amount_to_spawn; ++i) {

		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &base_spin = entry->particle_data->base_spin[index];
		float &local_rotations_x = entry->particle_data->local_rotations_x[index];
		float &local_rotations_y = entry->particle_data->local_rotations_y[index];
		float &local_rotations_z = entry->particle_data->local_rotations_z[index];

		//----Spin
		base_spin = RandomRange(&random, -spin_variation, spin_variation) + roll_spin;
		if (emitter.state_flags & tfxEmitterStateFlags_can_spin_pitch_and_yaw) {
			float &base_pitch_spin = entry->particle_data->base_pitch_spin[index];
			float &base_yaw_spin = entry->particle_data->base_yaw_spin[index];
			base_pitch_spin = RandomRange(&random, -spin_pitch_variation, spin_pitch_variation) + pitch_spin;
			base_yaw_spin = RandomRange(&random, -spin_yaw_variation, spin_yaw_variation) + yaw_spin;
		}
		if (angle_settings & tfxAngleSettingFlags_specify_roll)
			local_rotations_z = emitter.angle_offsets.roll;
		if (angle_settings & tfxAngleSettingFlags_specify_pitch)
			local_rotations_x = emitter.angle_offsets.pitch;
		if (angle_settings & tfxAngleSettingFlags_specify_yaw)
			local_rotations_y = emitter.angle_offsets.yaw;
		if (angle_settings & tfxAngleSettingFlags_random_pitch)
			local_rotations_x = RandomRange(&random, emitter.angle_offsets.pitch);
		if (angle_settings & tfxAngleSettingFlags_random_yaw)
			local_rotations_y = RandomRange(&random, emitter.angle_offsets.yaw);
		if (angle_settings & tfxAngleSettingFlags_random_roll)
			local_rotations_z = RandomRange(&random, emitter.angle_offsets.roll);

	}

}

void SpawnParticlePoint2d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 9 + emitter.seed_index);

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];

		tween = 1.f - entry->particle_data->age[index] / pm.frame_length;
		local_position_x = local_position_y = 0;
		if (emitter.property_flags & tfxEmitterPropertyFlags_relative_position)
			local_position_x = local_position_y = 0;
		else {
			tfx_vec2_t lerp_position = tfx_InterpolateVec2(tween, emitter.captured_position.xy(), emitter.world_position.xy());
			if (emitter.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) {
				local_position_x = lerp_position.x;
				local_position_y = lerp_position.y;
			}
			else {
				tfx_vec2_t rotvec = RotateVectorQuaternion2d(&emitter.rotation, -emitter.handle.xy());
				local_position_x = rotvec.x + lerp_position.x;
				local_position_y = rotvec.y + lerp_position.y;
			}
		}
	}

}

void SpawnParticlePoint3d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = 0.f;
	tfxU32 emitter_index = entry->emitter_index;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 10 + emitter.seed_index);
	tfx_effect_state_t *parent = &pm.effects[entry->parent_index];

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];
		float &local_position_z = entry->particle_data->position_z[index];

		tween = 1.f - entry->particle_data->age[index] / pm.frame_length;

		local_position_x = local_position_y = local_position_z = 0;
		if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			tfx_vec3_t lerp_position = tfx_InterpolateVec3(tween, emitter.captured_position, emitter.world_position);
			if (emitter.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) {
				local_position_x = lerp_position.x;
				local_position_y = lerp_position.y;
				local_position_z = lerp_position.z;
			}
			else {
				tfx_vec3_t rotvec = RotateVectorQuaternion(&emitter.rotation, emitter.handle);
				local_position_x = rotvec.x + lerp_position.x;
				local_position_y = rotvec.y + lerp_position.y;
				local_position_z = rotvec.z + lerp_position.z;
			}
		}
	}

}

void SpawnParticleOtherEmitter3d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = 0.f;
	tfxU32 emitter_index = entry->emitter_index;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 10 + emitter.seed_index);
	tfx_effect_state_t *parent = &pm.effects[entry->parent_index];

	if (emitter.spawn_locations_index == tfxINVALID) {
		entry->amount_to_spawn = 0;
		return;
	}

	tfx_soa_buffer_t &spawn_point_buffer = pm.particle_location_buffers[emitter.spawn_locations_index];
	tfx_spawn_points_soa_t &spawn_points = pm.particle_location_arrays[emitter.spawn_locations_index];
	if (spawn_point_buffer.current_size == 0) {
		entry->amount_to_spawn = 0;
		return;
	}

	if (entry->emission_type == tfxOtherEmitter) {
		emitter.grid_coords.x = emitter.grid_coords.y;
	}

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];
		float &local_position_z = entry->particle_data->position_z[index];

		int spawn_index = (int)emitter.grid_coords.x;
		spawn_index = GetCircularIndex(&spawn_point_buffer, spawn_index);
		local_position_x = spawn_points.position_x[spawn_index];
		local_position_y = spawn_points.position_y[spawn_index];
		local_position_z = spawn_points.position_z[spawn_index];

		emitter.grid_coords.x++;
		if ((tfxU32)emitter.grid_coords.x >= spawn_point_buffer.current_size) {
			emitter.grid_coords.x = 0.f;
		}

		if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			tween = 1.f - entry->particle_data->age[index] / pm.frame_length;
			float lerp_position_x = tfx_InterpolateFloat(tween, spawn_points.captured_position_x[spawn_index], spawn_points.position_x[spawn_index]);
			float lerp_position_y = tfx_InterpolateFloat(tween, spawn_points.captured_position_y[spawn_index], spawn_points.position_y[spawn_index]);
			float lerp_position_z = tfx_InterpolateFloat(tween, spawn_points.captured_position_z[spawn_index], spawn_points.position_z[spawn_index]);
			if (emitter.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) {
				local_position_x = lerp_position_x;
				local_position_y = lerp_position_y;
				local_position_z = lerp_position_z;
			}
			else {
				tfx_vec3_t rotvec = RotateVectorQuaternion(&emitter.rotation, emitter.handle);
				local_position_x = rotvec.x + lerp_position_x;
				local_position_y = rotvec.y + lerp_position_y;
				local_position_z = rotvec.z + lerp_position_z;
			}
		}

	}
}

void SpawnParticleOtherEmitter2d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = 0.f;
	tfxU32 emitter_index = entry->emitter_index;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 10 + emitter.seed_index);
	tfx_effect_state_t *parent = &pm.effects[entry->parent_index];

	if (emitter.spawn_locations_index == tfxINVALID) {
		entry->amount_to_spawn = 0;
		return;
	}

	tfx_soa_buffer_t &spawn_point_buffer = pm.particle_location_buffers[emitter.spawn_locations_index];
	tfx_spawn_points_soa_t &spawn_points = pm.particle_location_arrays[emitter.spawn_locations_index];
	if (spawn_point_buffer.current_size == 0) {
		entry->amount_to_spawn = 0;
		return;
	}

	if (entry->emission_type == tfxOtherEmitter) {
		emitter.grid_coords.x = emitter.grid_coords.y;
	}

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];

		int spawn_index = (int)emitter.grid_coords.x;
		spawn_index = GetCircularIndex(&spawn_point_buffer, spawn_index);
		local_position_x = spawn_points.position_x[spawn_index];
		local_position_y = spawn_points.position_y[spawn_index];

		emitter.grid_coords.x++;
		if ((tfxU32)emitter.grid_coords.x >= spawn_point_buffer.current_size) {
			emitter.grid_coords.x = 0.f;
		}

		if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			tween = 1.f - entry->particle_data->age[index] / pm.frame_length;
			float lerp_position_x = tfx_InterpolateFloat(tween, spawn_points.captured_position_x[spawn_index], spawn_points.position_x[spawn_index]);
			float lerp_position_y = tfx_InterpolateFloat(tween, spawn_points.captured_position_y[spawn_index], spawn_points.position_y[spawn_index]);
			if (emitter.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) {
				local_position_x = lerp_position_x;
				local_position_y = lerp_position_y;
			}
			else {
				tfx_vec2_t rotvec = RotateVectorQuaternion2d(&emitter.rotation, { emitter.handle.x, emitter.handle.y });
				local_position_x = rotvec.x + lerp_position_x;
				local_position_y = rotvec.y + lerp_position_y;
			}
		}
	}
}


void SpawnParticleLine2d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 11 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;
	const tfx_vec3_t &grid_points = properties.grid_points;
	const float emitter_size = emitter.emitter_size.y;
	const float grid_segment_size_y = -emitter_size / tfxMax(grid_points.y - 1.f, 1.f);;

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];

		tween = 1.f - entry->particle_data->age[index] / pm.frame_length;
		local_position_x = local_position_y = 0;

		if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

			if (!(emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
				emitter.grid_coords.y--;
				if (emitter.grid_coords.y < 0.f) {
					emitter.grid_coords.y = grid_points.y - 1;
				}
			}

			local_position_y = emitter.grid_coords.y * grid_segment_size_y;

			if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
				emitter.grid_coords.y++;
				if (emitter.grid_coords.y >= grid_points.y) {
					emitter.grid_coords.y = 0.f;
				}
			}

		}
		else {
			local_position_y = RandomRange(&random, -emitter_size);
		}

		//----TForm and Emission
		if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position) && !(emitter.property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
			tfx_vec2_t lerp_position = tfx_InterpolateVec2(tween, emitter.captured_position.xy(), emitter.world_position.xy());
			tfx_vec2_t pos = RotateVectorQuaternion2d(&emitter.rotation, tfx_vec2_t(local_position_x, local_position_y) + emitter.handle.xy());
			local_position_x = lerp_position.x + pos.x * entry->overal_scale;
			local_position_y = lerp_position.y + pos.y * entry->overal_scale;
		}

	}

}

void SpawnParticleLine3d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 12 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;
	const tfx_vec3_t &grid_points = properties.grid_points;
	const float emitter_size = emitter.emitter_size.y;
	const float grid_segment_size_y = emitter_size / tfxMax(grid_points.y - 1.f, 1.f);

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];
		float &local_position_z = entry->particle_data->position_z[index];

		local_position_x = local_position_y = local_position_z = 0;

		if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

			if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
				emitter.grid_coords.y = (float)RandomRange(&random, (tfxU32)grid_points.y);
				local_position_y = emitter.grid_coords.y * grid_segment_size_y;
			}
			else {

				if (!(emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					emitter.grid_coords.y--;
					if (emitter.grid_coords.y < 0.f) {
						emitter.grid_coords.y = grid_points.y - 1;
					}
				}

				local_position_y = emitter.grid_coords.y * grid_segment_size_y;

				if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					emitter.grid_coords.y++;
					if (emitter.grid_coords.y >= grid_points.y) {
						emitter.grid_coords.y = 0.f;
					}
				}
			}

		}
		else {
			local_position_y = RandomRange(&random, 0.f, emitter_size);
		}

		//----TForm and Emission
		if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position) && !(emitter.property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
			tfx_vec3_t lerp_position = tfx_InterpolateVec3(tween, emitter.captured_position, emitter.world_position);
			tfx_vec3_t position_plus_handle = tfx_vec3_t(local_position_x, local_position_y, local_position_z) + emitter.handle;
			tfx_vec3_t pos = RotateVectorQuaternion(&emitter.rotation, position_plus_handle);
			local_position_x = lerp_position.x + pos.x * entry->overal_scale;
			local_position_y = lerp_position.y + pos.y * entry->overal_scale;
			local_position_z = lerp_position.z + pos.z * entry->overal_scale;
		}

		tween += entry->qty_step_size;
	}

}

void SpawnParticleArea2d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 13 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;
	const tfx_vec3_t &grid_points = properties.grid_points;
	const float grid_segment_size_x = emitter.emitter_size.x / tfxMax(grid_points.x - 1.f, 1.f);;
	const float grid_segment_size_y = emitter.emitter_size.y / tfxMax(grid_points.y - 1.f, 1.f);;
	tfx_vec3_t half_emitter_size = emitter.emitter_size * .5f;

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];

		local_position_x = local_position_y = 0;

		tfx_vec2_t position = tfx_vec2_t(0.f, 0.f);

		if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

			if (emitter.property_flags & tfxEmitterPropertyFlags_fill_area) {
				if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
					emitter.grid_coords.x = (float)RandomRange(&random, (tfxU32)grid_points.x);
					emitter.grid_coords.y = (float)RandomRange(&random, (tfxU32)grid_points.y);
					local_position_x = emitter.grid_coords.x * grid_segment_size_x;
					local_position_y = emitter.grid_coords.y * grid_segment_size_y;
				}
				else {
					if (!(emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
						emitter.grid_coords.x--;
						if (emitter.grid_coords.x < 0.f) {
							emitter.grid_coords.y--;
							emitter.grid_coords.x = grid_points.x - 1;
							if (emitter.grid_coords.y < 0.f)
								emitter.grid_coords.y = grid_points.y - 1;
						}
					}

					local_position_x = emitter.grid_coords.x * grid_segment_size_x;
					local_position_y = emitter.grid_coords.y * grid_segment_size_y;

					if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
						emitter.grid_coords.x++;
						if (emitter.grid_coords.x == grid_points.x) {
							emitter.grid_coords.y++;
							emitter.grid_coords.x = 0.f;
							if (emitter.grid_coords.y >= grid_points.y)
								emitter.grid_coords.y = 0.f;
						}
					}
				}
				local_position_x -= half_emitter_size.x;
				local_position_y -= half_emitter_size.y;
			}
			else {
				if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
					tfxU32 side = RandomRange(&random, (tfxU32)4);
					if (side == 0) {
						//left side
						emitter.grid_coords.x = 0.f;
						emitter.grid_coords.y = (float)RandomRange(&random, (tfxU32)grid_points.y);
					}
					else if (side == 1) {
						//right side
						emitter.grid_coords.x = grid_points.x - 1;
						emitter.grid_coords.y = (float)RandomRange(&random, (tfxU32)grid_points.y);
					}
					else if (side == 2) {
						//top side
						emitter.grid_coords.x = (float)RandomRange(&random, (tfxU32)grid_points.x);
						emitter.grid_coords.y = 0.f;
					}
					else if (side == 3) {
						//bottom side
						emitter.grid_coords.x = (float)RandomRange(&random, (tfxU32)grid_points.x);
						emitter.grid_coords.y = grid_points.y - 1;
					}
					local_position_x = emitter.grid_coords.x * grid_segment_size_x;
					local_position_y = emitter.grid_coords.y * grid_segment_size_y;
				}
				else {
					if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {

						emitter.grid_direction.x = 1;
						emitter.grid_direction.y = 0;
						if (emitter.grid_coords.x == grid_points.x - 1 && emitter.grid_coords.y >= 0 && emitter.grid_coords.y < grid_points.y - 1) {
							emitter.grid_direction.x = 0;
							emitter.grid_direction.y = 1;
						}
						else if (emitter.grid_coords.x > 0 && emitter.grid_coords.x < grid_points.x && emitter.grid_coords.y == grid_points.y - 1) {
							emitter.grid_direction.x = -1;
							emitter.grid_direction.y = 0;
						}
						else if (emitter.grid_coords.x == 0 && emitter.grid_coords.y > 0 && emitter.grid_coords.y < grid_points.y) {
							emitter.grid_direction.x = 0;
							emitter.grid_direction.y = -1;
						}

					}
					else {

						emitter.grid_direction.x = -1;
						emitter.grid_direction.y = 0;
						if (emitter.grid_coords.x == grid_points.x - 1 && emitter.grid_coords.y > 0 && emitter.grid_coords.y < grid_points.y) {
							emitter.grid_direction.x = 0;
							emitter.grid_direction.y = -1;
						}
						else if (emitter.grid_coords.x >= 0 && emitter.grid_coords.x < grid_points.x - 1 && emitter.grid_coords.y == grid_points.y - 1) {
							emitter.grid_direction.x = 1;
							emitter.grid_direction.y = 0;
						}
						else if (emitter.grid_coords.x == 0 && emitter.grid_coords.y >= 0 && emitter.grid_coords.y < grid_points.y - 1) {
							emitter.grid_direction.x = 0;
							emitter.grid_direction.y = 1;
						}

					}

					emitter.grid_coords += emitter.grid_direction;
					local_position_x = position.x + (emitter.grid_coords.x * grid_segment_size_x);
					local_position_y = position.y + (emitter.grid_coords.y * grid_segment_size_y);
				}
				local_position_x -= half_emitter_size.x;
				local_position_y -= half_emitter_size.y;
			}
		}
		else {
			if (emitter.property_flags & tfxEmitterPropertyFlags_fill_area) {
				position.x = RandomRange(&random, emitter.emitter_size.x);
				position.y = RandomRange(&random, emitter.emitter_size.y);
			}
			else {
				//Spawn on one of 4 edges of the area
				tfxU32 side = RandomRange(&random, (tfxU32)4);
				if (side == 0) {
					//left side
					position.x = 0.f;
					position.y = RandomRange(&random, emitter.emitter_size.y);
				}
				else if (side == 1) {
					//right side
					position.x = emitter.emitter_size.x;
					position.y = RandomRange(&random, emitter.emitter_size.y);
				}
				else if (side == 2) {
					//top side
					position.x = RandomRange(&random, emitter.emitter_size.x);
					position.y = 0.f;
				}
				else if (side == 3) {
					//bottom side
					position.x = RandomRange(&random, emitter.emitter_size.x);
					position.y = emitter.emitter_size.y;
				}
			}

			local_position_x = position.x - half_emitter_size.x;
			local_position_y = position.y - half_emitter_size.y;
		}

		//----TForm and Emission
		if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			tfx_vec2_t lerp_position = tfx_InterpolateVec2(tween, emitter.captured_position.xy(), emitter.world_position.xy());
			tfx_vec2_t pos = RotateVectorQuaternion2d(&emitter.rotation, tfx_vec2_t(local_position_x, local_position_y) + emitter.handle.xy());
			local_position_x = lerp_position.x + pos.x * entry->overal_scale;
			local_position_y = lerp_position.y + pos.y * entry->overal_scale;
		}

		tween += entry->qty_step_size;
	}

}

void SpawnParticleArea3d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 14 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;
	const tfx_vec3_t &grid_points = properties.grid_points;
	const float grid_segment_size_x = emitter.emitter_size.x / tfxMax(grid_points.x - 1.f, 1.f);
	const float grid_segment_size_y = emitter.emitter_size.y / tfxMax(grid_points.y - 1.f, 1.f);
	const float grid_segment_size_z = emitter.emitter_size.z / tfxMax(grid_points.z - 1.f, 1.f);

	tfx_vec3_t half_emitter_size = emitter.emitter_size * .5f;

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];
		float &local_position_z = entry->particle_data->position_z[index];
		tfxU32 &velocity_normal_packed = entry->particle_data->velocity_normal[index];

		local_position_x = local_position_y = local_position_z = 0;

		tfx_vec3_t position;

		if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

			if (emitter.property_flags & tfxEmitterPropertyFlags_fill_area) {

				if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
					emitter.grid_coords.x = (float)RandomRange(&random, (tfxU32)grid_points.x);
					emitter.grid_coords.y = (float)RandomRange(&random, (tfxU32)grid_points.y);
					emitter.grid_coords.z = (float)RandomRange(&random, (tfxU32)grid_points.z);

					local_position_x = emitter.grid_coords.x * grid_segment_size_x - half_emitter_size.x;
					local_position_y = emitter.grid_coords.y * grid_segment_size_y - half_emitter_size.y;
					local_position_z = emitter.grid_coords.z * grid_segment_size_z - half_emitter_size.z;
				}
				else {
					if (!(emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
						emitter.grid_coords.x--;
						if (emitter.grid_coords.x < 0.f) {
							emitter.grid_coords.y--;
							emitter.grid_coords.x = grid_points.x - 1;
							if (emitter.grid_coords.y < 0.f) {
								emitter.grid_coords.z--;
								emitter.grid_coords.y = grid_points.y - 1;
								if (emitter.grid_coords.z < 0.f)
									emitter.grid_coords.z = grid_points.z - 1;
							}
						}
					}

					local_position_x = position.x + (emitter.grid_coords.x * grid_segment_size_x) - half_emitter_size.x;
					local_position_y = position.y + (emitter.grid_coords.y * grid_segment_size_y) - half_emitter_size.y;
					local_position_z = position.z + (emitter.grid_coords.z * grid_segment_size_z) - half_emitter_size.z;

					if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
						emitter.grid_coords.x++;
						if (emitter.grid_coords.x == grid_points.x) {
							emitter.grid_coords.y++;
							emitter.grid_coords.x = 0.f;
							if (emitter.grid_coords.y >= grid_points.y) {
								emitter.grid_coords.z++;
								emitter.grid_coords.y = 0.f;
								if (emitter.grid_coords.z >= grid_points.z)
									emitter.grid_coords.z = 0.f;
							}
						}
					}
				}

			}
			else {
				if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
					//Spawn on one of 6 edges of the cuboid
					tfxU32 side = RandomRange(&random, (emitter.property_flags & tfxEmitterPropertyFlags_area_open_ends) ? (tfxU32)4 : (tfxU32)6);
					tfx_vec3_t velocity_normal;
					if (side == 0) {
						//left side
						emitter.grid_coords.x = 0.f;
						emitter.grid_coords.y = (float)RandomRange(&random, (tfxU32)grid_points.y);
						emitter.grid_coords.z = (float)RandomRange(&random, (tfxU32)grid_points.z);
						velocity_normal.x = -1.f;
					}
					else if (side == 1) {
						//right side
						emitter.grid_coords.x = grid_points.x - 1;
						emitter.grid_coords.y = (float)RandomRange(&random, (tfxU32)grid_points.y);
						emitter.grid_coords.z = (float)RandomRange(&random, (tfxU32)grid_points.z);
						velocity_normal.x = 1.f;
					}
					else if (side == 2) {
						//top side
						emitter.grid_coords.x = (float)RandomRange(&random, (tfxU32)grid_points.x);
						emitter.grid_coords.y = 0.f;
						emitter.grid_coords.z = (float)RandomRange(&random, (tfxU32)grid_points.z);
						velocity_normal.y = -1.f;
					}
					else if (side == 3) {
						//bottom side
						emitter.grid_coords.x = (float)RandomRange(&random, (tfxU32)grid_points.x);
						emitter.grid_coords.y = grid_points.y - 1;
						emitter.grid_coords.z = (float)RandomRange(&random, (tfxU32)grid_points.z);
						velocity_normal.y = 1.f;
					}
					else if (side == 4) {
						//End far
						emitter.grid_coords.x = (float)RandomRange(&random, (tfxU32)grid_points.x);
						emitter.grid_coords.y = (float)RandomRange(&random, (tfxU32)grid_points.y);
						emitter.grid_coords.z = grid_points.z - 1;
						velocity_normal.z = 1.f;
					}
					else if (side == 5) {
						//End near
						emitter.grid_coords.x = (float)RandomRange(&random, (tfxU32)grid_points.x);
						emitter.grid_coords.y = (float)RandomRange(&random, (tfxU32)grid_points.y);
						emitter.grid_coords.z = 0.f;
						velocity_normal.z = -1.f;
					}
					velocity_normal_packed = tfx__pack10bit_unsigned(&velocity_normal);
					local_position_x = emitter.grid_coords.x * grid_segment_size_x - half_emitter_size.x;
					local_position_y = emitter.grid_coords.y * grid_segment_size_y - half_emitter_size.y;
					local_position_z = emitter.grid_coords.z * grid_segment_size_z - half_emitter_size.z;
				}
				else {
					tfx_vec3_t velocity_normal;
					if (!(emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
						if ((emitter.grid_coords.z > 0 && emitter.grid_coords.z < grid_points.z - 1) || emitter.property_flags & tfxEmitterPropertyFlags_area_open_ends) {
							emitter.grid_coords.x -= emitter.grid_coords.y == 0 || emitter.grid_coords.y == grid_points.y - 1 ? 1.f : grid_points.x - 1;
							if (emitter.grid_coords.x < 0.f) {
								emitter.grid_coords.y--;
								emitter.grid_coords.x = grid_points.x - 1;
								if (emitter.grid_coords.y < 0.f) {
									emitter.grid_coords.z--;
									emitter.grid_coords.y = grid_points.y - 1;
									if (emitter.grid_coords.z < 0.f)
										emitter.grid_coords.z = grid_points.z - 1;
								}
							}
						}
						else {
							emitter.grid_coords.x--;
							if (emitter.grid_coords.x < 0.f) {
								emitter.grid_coords.y--;
								emitter.grid_coords.x = grid_points.x - 1;
								if (emitter.grid_coords.y < 0.f) {
									emitter.grid_coords.z--;
									emitter.grid_coords.y = grid_points.y - 1;
									if (emitter.grid_coords.z < 0.f)
										emitter.grid_coords.z = grid_points.z - 1;
								}
							}
						}
					}

					local_position_x = position.x + (emitter.grid_coords.x * grid_segment_size_x) - half_emitter_size.x;
					local_position_y = position.y + (emitter.grid_coords.y * grid_segment_size_y) - half_emitter_size.y;
					local_position_z = position.z + (emitter.grid_coords.z * grid_segment_size_z) - half_emitter_size.z;

					velocity_normal.x = emitter.grid_coords.x == 0 ? -1.f : (emitter.grid_coords.x == grid_points.x - 1 ? 1.f : 0.f);
					velocity_normal.y = emitter.grid_coords.y == 0 ? -1.f : (emitter.grid_coords.y == grid_points.y - 1 ? 1.f : 0.f);
					velocity_normal.z = emitter.grid_coords.z == 0 ? -1.f : (emitter.grid_coords.z == grid_points.z - 1 ? 1.f : 0.f);

					velocity_normal_packed = tfx__pack10bit_unsigned(&velocity_normal);

					if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
						if ((emitter.grid_coords.z > 0 && emitter.grid_coords.z < grid_points.z - 1) || emitter.property_flags & tfxEmitterPropertyFlags_area_open_ends) {
							emitter.grid_coords.x += emitter.grid_coords.y == 0 || emitter.grid_coords.y == grid_points.y - 1 ? 1.f : grid_points.x - 1;
							if (emitter.grid_coords.x >= grid_points.x) {
								emitter.grid_coords.y++;
								emitter.grid_coords.x = 0.f;
								if (emitter.grid_coords.y >= grid_points.y) {
									emitter.grid_coords.z++;
									emitter.grid_coords.y = 0.f;
									if (emitter.grid_coords.z >= grid_points.z) {
										emitter.grid_coords.z = 0.f;
									}
								}
							}
						}
						else {
							emitter.grid_coords.x++;
							if (emitter.grid_coords.x == grid_points.x) {
								emitter.grid_coords.y++;
								emitter.grid_coords.x = 0.f;
								if (emitter.grid_coords.y >= grid_points.y) {
									emitter.grid_coords.z++;
									emitter.grid_coords.y = 0.f;
									if (emitter.grid_coords.z >= grid_points.z) {
										emitter.grid_coords.z = 0.f;
									}
								}
							}
						}
					}
				}

			}
		}
		else {
			if (emitter.property_flags & tfxEmitterPropertyFlags_fill_area) {
				position.x = RandomRange(&random, -half_emitter_size.x, half_emitter_size.x);
				position.y = RandomRange(&random, -half_emitter_size.y, half_emitter_size.y);
				position.z = RandomRange(&random, -half_emitter_size.z, half_emitter_size.z);
			}
			else {
				tfx_vec3_t velocity_normal;
				//Spawn on one of 6 edges of the cuboid
				tfxU32 side = RandomRange(&random, (emitter.property_flags & tfxEmitterPropertyFlags_area_open_ends) ? (tfxU32)4 : (tfxU32)6);
				if (side == 0) {
					//left side
					position.x = -half_emitter_size.x;
					position.y = RandomRange(&random, -half_emitter_size.y, half_emitter_size.y);
					position.z = RandomRange(&random, -half_emitter_size.z, half_emitter_size.z);
					velocity_normal.x = -1.f;
				}
				else if (side == 1) {
					//right side
					position.x = half_emitter_size.x;
					position.y = RandomRange(&random, -half_emitter_size.y, half_emitter_size.y);
					position.z = RandomRange(&random, -half_emitter_size.z, half_emitter_size.z);
					velocity_normal.x = 1.f;
				}
				else if (side == 2) {
					//top side
					position.x = RandomRange(&random, -half_emitter_size.y, half_emitter_size.y);
					position.y = -half_emitter_size.y;
					position.z = RandomRange(&random, -half_emitter_size.z, half_emitter_size.z);
					velocity_normal.y = -1.f;
				}
				else if (side == 3) {
					//bottom side
					position.x = RandomRange(&random, -half_emitter_size.y, half_emitter_size.y);
					position.y = half_emitter_size.y;
					position.z = RandomRange(&random, -half_emitter_size.z, half_emitter_size.z);
					velocity_normal.y = 1.f;
				}
				else if (side == 4) {
					//End far
					position.x = RandomRange(&random, -half_emitter_size.y, half_emitter_size.y);
					position.y = RandomRange(&random, -half_emitter_size.y, half_emitter_size.y);
					position.z = half_emitter_size.z;
					velocity_normal.z = 1.f;
				}
				else if (side == 5) {
					//End near
					position.x = RandomRange(&random, -half_emitter_size.y, half_emitter_size.y);
					position.y = RandomRange(&random, -half_emitter_size.y, half_emitter_size.y);
					position.z = -half_emitter_size.z;
					velocity_normal.z = -1.f;
				}
				velocity_normal_packed = tfx__pack10bit_unsigned(&velocity_normal);
			}

			local_position_x = position.x;
			local_position_y = position.y;
			local_position_z = position.z;
		}

		//----TForm and Emission
		if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			tfx_vec3_t lerp_position = tfx_InterpolateVec3(tween, emitter.captured_position, emitter.world_position);
			tfx_vec3_t position_plus_handle = tfx_vec3_t(local_position_x, local_position_y, local_position_z) + emitter.handle;
			tfx_vec3_t pos = RotateVectorQuaternion(&emitter.rotation, position_plus_handle);
			local_position_x = lerp_position.x + pos.x * entry->overal_scale;
			local_position_y = lerp_position.y + pos.y * entry->overal_scale;
			local_position_z = lerp_position.z + pos.z * entry->overal_scale;
		}

		tween += entry->qty_step_size;
	}

}

void SpawnParticleEllipse2d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 15 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;
	const tfx_vec3_t &grid_points = properties.grid_points;
	const tfx_vec2_t &emitter_size = emitter.emitter_size.xy();
	float arc_size = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.arc_size, emitter.frame);
	float arc_offset = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.arc_offset, emitter.frame);
	const float grid_segment_size_x = arc_size / tfxMax(grid_points.x, 1.f);

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];

		local_position_x = local_position_y = 0;

		tfx_vec2_t half_emitter_size = (emitter_size * .5f);
		tfx_vec2_t position = tfx_vec2_t(0.f, 0.f);

		if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid && !(emitter.property_flags & tfxEmitterPropertyFlags_fill_area) && !(emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_random)) {

			emitter.grid_coords.y = 0.f;

			if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
				emitter.grid_coords.x--;
				if (emitter.grid_coords.x < 0.f) {
					emitter.grid_coords.x = grid_points.x - 1;
				}
			}

			float th = emitter.grid_coords.x * grid_segment_size_x + arc_offset;
			local_position_x = cosf(th) * half_emitter_size.x;
			local_position_y = -sinf(th) * half_emitter_size.y;

			if (!(emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
				emitter.grid_coords.x++;
				if (emitter.grid_coords.x >= grid_points.x) {
					emitter.grid_coords.x = 0.f;
				}
			}

		}
		else if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid && !(emitter.property_flags & tfxEmitterPropertyFlags_fill_area) && emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
			float th = (float)RandomRange(&random, (tfxU32)grid_points.x) * grid_segment_size_x + arc_offset;
			local_position_x = cosf(th) * half_emitter_size.x;
			local_position_y = -sinf(th) * half_emitter_size.y;
		}
		else if (!(emitter.property_flags & tfxEmitterPropertyFlags_fill_area)) {
			float th = RandomRange(&random, arc_size) + arc_offset;

			local_position_x = cosf(th) * half_emitter_size.x;
			local_position_y = -sinf(th) * half_emitter_size.y;

		}
		else {
			local_position_x = RandomRange(&random, 0.f, emitter_size.x);
			local_position_y = RandomRange(&random, 0.f, emitter_size.y);

			while ((powf(local_position_x - half_emitter_size.x, 2) / powf(half_emitter_size.x, 2)) + (powf(local_position_y - half_emitter_size.y, 2) / powf(half_emitter_size.y, 2)) > 1) {
				local_position_x = RandomRange(&random, 0.f, emitter_size.x);
				local_position_y = RandomRange(&random, 0.f, emitter_size.y);
			}
			local_position_x -= half_emitter_size.x;
			local_position_y -= half_emitter_size.y;
		}

		//----TForm and Emission
		if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			tfx_vec2_t lerp_position = tfx_InterpolateVec2(tween, emitter.captured_position.xy(), emitter.world_position.xy());
			tfx_vec2_t pos = RotateVectorQuaternion2d(&emitter.rotation, tfx_vec2_t(local_position_x, local_position_y) + emitter.handle.xy());
			local_position_x = lerp_position.x + pos.x * entry->overal_scale;
			local_position_y = lerp_position.y + pos.y * entry->overal_scale;
		}


		tween += entry->qty_step_size;
	}

}

void SpawnParticlePath2d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 26 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;
	tfx_vec3_t half_emitter_size = emitter.emitter_size * .5f;
	const tfx_vec3_t &grid_points = properties.grid_points;
	tfx_emitter_path_t *path = &pm.library->paths[emitter.path_attributes];
	float total_grid_points = (float)path->node_count - 3.f;
	float increment = 1.f / grid_points.x;
	float arc_size = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.arc_size, emitter.frame);
	float arc_offset = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.arc_offset, emitter.frame);
	float extrusion = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.extrusion, emitter.frame);
	tfx_vec2_t point;
	bool has_rotation = path->rotation_range > 0 || path->rotation_pitch != 0 || path->rotation_yaw != 0;

	float emission_direction;
	float emission_angle_variation;
	float range;
	tfx_vec2_t velocity_direction;

	if (properties.emission_direction == tfxPathGradient) {
		emission_direction = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.emission_pitch, emitter.frame);
		emission_angle_variation = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.emission_range, emitter.frame);
		range = emission_angle_variation * .5f;
	}

	if (path->rotation_cycle_length > 0) {
		if ((emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid && emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_random) ||
			!(emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid)
			) {
			for (int qi = 0; qi != emitter.active_paths; ++qi) {
				int index = (qi + emitter.path_start_index) % path->maximum_active_paths;
				emitter.path_quaternions[qi].age += pm.frame_length;
				if (emitter.path_quaternions[qi].age >= path->rotation_cycle_length) {
					tfx_quaternion_t q = GetPathRotation2d(&random, path->rotation_range, path->rotation_pitch);
					emitter.path_quaternions[qi].quaternion = tfx__pack8bit_quaternion(q);
					emitter.path_quaternions[qi].age = 0.f;
					emitter.path_cycle_count--;
				}
			}
		}
	}

	tfxU32 qi = (emitter.path_start_index + emitter.last_path_index) % path->maximum_active_paths;
	TFX_ASSERT(qi < path->maximum_active_paths);

	if (path->rotation_stagger > 0 && emitter.path_stagger_counter >= path->rotation_stagger) {
		if (emitter.active_paths < path->maximum_active_paths && (emitter.path_cycle_count > 0 || path->maximum_paths == 0)) {
			emitter.last_path_index = emitter.active_paths;
			qi = (emitter.path_start_index + emitter.active_paths++) % path->maximum_active_paths;
			TFX_ASSERT(qi < path->maximum_active_paths);
			tfx_quaternion_t q = GetPathRotation2d(&random, path->rotation_range, path->rotation_pitch);
			emitter.path_quaternions[qi].quaternion = tfx__pack8bit_quaternion(q);
			emitter.path_quaternions[qi].cycles = 0;
			if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
				emitter.path_quaternions[qi].grid_coord = 0.f;
			}
			else if (!(emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
				emitter.path_quaternions[qi].grid_coord = total_grid_points - increment;
			}
			emitter.path_cycle_count--;
			emitter.path_stagger_counter = 0.f;
		}
	}

	int dead_paths = 0;
	int node = 0;
	float t = 0.f;

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];
		float &path_position = entry->particle_data->path_position[index];
		float &path_offset = entry->particle_data->path_offset[index];

		if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid && emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
			node = RandomRange(&random, path->node_count - 3);
			t = (float)RandomRange(&random, (int)grid_points.x) * increment;
			path_position = (float)node + t;
			point = tfx__catmull_rom_spline_2d_soa(path->node_soa.x, path->node_soa.y, node, t);
		}
		else if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {
			float &grid_coord = emitter.path_quaternions[qi].grid_coord;
			bool new_path = false;
			if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
				grid_coord += increment;
				if (grid_coord >= total_grid_points) {
					grid_coord = 0.f;
					new_path = true;
				}
			}
			else if (!(emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
				grid_coord -= increment;
				if (grid_coord < 0) {
					grid_coord = total_grid_points - increment;
					new_path = true;
				}
			}
			if (new_path) {
				if (emitter.state_flags & tfxEmitterStateFlags_has_rotated_path && path->rotation_stagger == 0) {
					if (path->maximum_paths == 0 || emitter.path_cycle_count > 0) {
						tfx_quaternion_t q = GetPathRotation2d(&random, path->rotation_range, path->rotation_pitch);
						emitter.path_quaternions[qi].quaternion = tfx__pack8bit_quaternion(q);
						emitter.path_cycle_count--;
					}
					else {
						dead_paths++;
						emitter.path_quaternions[qi].cycles = tfxINVALID;
						entry->particle_data->flags[index] |= tfxParticleFlags_remove;
					}
				}
				else {
					emitter.path_quaternions[qi].cycles = tfxINVALID;
					dead_paths++;
				}
			}
			node = (int)grid_coord;
			t = grid_coord - (int)grid_coord;
			path_position = (float)node + t;
			point = tfx__catmull_rom_spline_2d_soa(path->node_soa.x, path->node_soa.y, node, t);
		}
		else {
			node = RandomRange(&random, path->node_count - 3);
			t = GenerateRandom(&random);
			path_position = (float)node + t;

			point = tfx__catmull_rom_spline_2d_soa(path->node_soa.x, path->node_soa.y, node, t);
		}

		if (path->extrusion_type == tfxExtrusionLinear) {
			if (properties.emission_direction == tfxPathGradient) {
				velocity_direction = tfx__catmull_rom_spline_gradient_2d_soa(&path->node_soa.x[node], &path->node_soa.y[node], t);
				velocity_direction = tfx__normalize_vec2(&velocity_direction);
			}
			float radius = extrusion * .5f;
			path_offset = RandomRange(&random, -radius, radius);
			local_position_x = point.x + path_offset;
			local_position_x *= emitter.emitter_size.x;
			local_position_y = point.y * emitter.emitter_size.y;
		}
		else {
			path_offset = RandomRange(&random, arc_size) + arc_offset;
			float length_squared = point.x * point.x + point.y * point.y;
			float radius = length_squared == 0.f ? 0.f : 1.f / tfx__quake_sqrt(length_squared);
			float angle = atan2f(point.y, point.x) + path_offset;
			float rx = cosf(angle);
			float ry = sinf(angle);
			local_position_x = rx * radius * emitter.emitter_size.x;
			local_position_y = ry * radius * emitter.emitter_size.y;
			if (properties.emission_direction == tfxPathGradient) {
				velocity_direction = tfx__catmull_rom_spline_gradient_2d_soa(&path->node_soa.x[node], &path->node_soa.y[node], t);
				velocity_direction = tfx__normalize_vec2(&velocity_direction);
				rx = cosf(path_offset);
				ry = sinf(path_offset);
				tfx_vec2_t v = velocity_direction;
				velocity_direction.x = v.x * rx - v.y * ry;
				velocity_direction.y = v.x * ry + v.y * rx;
			}
		}

		/*
		if (path->extrusion_type == tfxExtrusionArc) {
			tfxWideFloat radius = tfxWideAdd(tfxWideMul(point_x.m, point_x.m), tfxWideMul(point_y.m, point_y.m));
			tfxWideFloat length_mask = tfxWideGreater(radius, tfxWideSetZero);
			radius = tfxWideMul(tfxWideRSqrt(radius), radius);
			tfxWideArray angle;
			tfxWideArray rx;
			tfxWideArray ry;
			angle.m = tfxWideAtan2(point_y.m, point_x.m);
			angle.m = tfxWideAnd(length_mask, angle.m);
			angle.m = tfxWideAdd(angle.m, path_offset);
			tfxWideSinCos(angle.m, &ry.m, &rx.m);
			local_position_x = tfxWideMul(rx.m, radius);
			local_position_y = tfxWideMul(ry.m, radius);
		} else {
			local_position_x = tfxWideAdd(point_x.m, path_offset);
			local_position_y = point_y.m;
		}
		*/

		if (emitter.state_flags & tfxEmitterStateFlags_has_rotated_path && emitter.active_paths > 0) {
			if (emitter.path_quaternions[qi].cycles == tfxINVALID) {
				entry->particle_data->flags[index] |= tfxParticleFlags_remove;
				emitter.last_path_index++;
				emitter.last_path_index %= emitter.active_paths;
				qi = (emitter.path_start_index + emitter.last_path_index) % path->maximum_active_paths;
				TFX_ASSERT(qi < path->maximum_active_paths);
				continue;
			}
			tfx_quaternion_t q = tfx__unpack8bit_quaternion(emitter.path_quaternions[qi].quaternion);
			entry->particle_data->quaternion[index] = emitter.path_quaternions[qi].quaternion;
			tfx_vec2_t rp = { local_position_x, local_position_y };
			rp = RotateVectorQuaternion2d(&q, rp);
			local_position_x = rp.x;
			local_position_y = rp.y;
			emitter.last_path_index++;
			emitter.last_path_index %= emitter.active_paths;
			qi = (emitter.path_start_index + emitter.last_path_index) % path->maximum_active_paths;
			TFX_ASSERT(qi < path->maximum_active_paths);
		}

		if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			tfx_vec2_t lerp_position = tfx_InterpolateVec2(tween, emitter.captured_position.xy(), emitter.world_position.xy());
			tfx_vec2_t position_plus_handle = tfx_vec2_t(local_position_x, local_position_y) + emitter.handle.xy();
			tfx_vec2_t pos = RotateVectorQuaternion2d(&emitter.rotation, position_plus_handle);
			local_position_x = lerp_position.x + pos.x * entry->overal_scale;
			local_position_y = lerp_position.y + pos.y * entry->overal_scale;
			if (properties.emission_direction == tfxPathGradient) {
				tfx_quaternion_t offset_quaternion;
				ToQuaternion2d(&offset_quaternion, emission_direction);
				tfx_vec2_t rotated_normal = RotateVectorQuaternion2d(&offset_quaternion, velocity_direction);
				float direction = atan2f(rotated_normal.y, rotated_normal.x);
				if (range != 0.f) {
					direction = RandomRange(&random, -range, range) + direction;
				}
				entry->particle_data->local_rotations_x[index] = direction;
			}
		}
		else if (properties.emission_direction == tfxPathGradient) {
			tfx_quaternion_t offset_quaternion;
			ToQuaternion2d(&offset_quaternion, emission_direction);
			tfx_vec2_t rotated_normal = RotateVectorQuaternion2d(&offset_quaternion, velocity_direction);
			float direction = atan2f(rotated_normal.y, rotated_normal.x);
			if (range != 0.f) {
				direction = RandomRange(&random, -range, range) + direction;
			}
			entry->particle_data->local_rotations_x[index] = direction;
		}

		tween += entry->qty_step_size;
	}

	if (dead_paths > 0 && emitter.active_paths > 0) {
		tfxU32 offset = 0;
		for (int qi = emitter.active_paths - 1; qi >= 0; --qi) {
			int index = (emitter.path_start_index + qi) % path->maximum_active_paths;
			if (emitter.path_quaternions[index].cycles == tfxINVALID) {
				offset++;
			}
			else if (offset > 0) {
				tfxU32 next_index = (qi + offset + emitter.path_start_index) % path->maximum_active_paths;
				TFX_ASSERT(next_index < path->maximum_active_paths);
				emitter.path_quaternions[next_index] = emitter.path_quaternions[index];
			}
		}
		emitter.path_start_index = (emitter.path_start_index + offset) % path->maximum_active_paths;
		emitter.active_paths -= tfx__Min(emitter.active_paths, offset);
		emitter.last_path_index = 0;
	}

}

void SpawnParticleEllipsoid(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 16 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;
	float arc_size = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.arc_size, emitter.frame);
	float arc_offset = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.arc_offset, emitter.frame);

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];
		float &local_position_z = entry->particle_data->position_z[index];

		local_position_x = local_position_y = local_position_z = 0;
		tfx_vec3_t lerp_position = tfx_InterpolateVec3(tween, emitter.captured_position, emitter.world_position);

		tfx_vec3_t half_emitter_size = emitter.emitter_size * .5f;
		tfx_vec3_t position;

		if (!(emitter.property_flags & tfxEmitterPropertyFlags_fill_area)) {
			float theta = RandomRange(&random, tfxPI2);
			float v = RandomRange(&random, 1.f);
			float phi = acosf(2.f * v - 1.f);
			float sin_theta = sinf(theta);
			float cos_theta = cosf(theta);
			float sin_phi = sinf(phi);
			float cos_phi = cosf(phi);
			local_position_x = half_emitter_size.x * sin_phi * cos_theta;
			local_position_y = half_emitter_size.y * sin_phi * sin_theta;
			local_position_z = half_emitter_size.z * cos_phi;
		}
		else {
			position.x = RandomRange(&random, -half_emitter_size.x, half_emitter_size.x);
			position.y = RandomRange(&random, -half_emitter_size.y, half_emitter_size.y);
			position.z = RandomRange(&random, -half_emitter_size.z, half_emitter_size.z);

			while (powf(position.x / half_emitter_size.x, 2.f) + powf(position.y / half_emitter_size.y, 2.f) + powf(position.z / half_emitter_size.z, 2.f) > 1.f) {
				position.x = RandomRange(&random, -half_emitter_size.x, half_emitter_size.x);
				position.y = RandomRange(&random, -half_emitter_size.y, half_emitter_size.y);
				position.z = RandomRange(&random, -half_emitter_size.z, half_emitter_size.z);
			}

			local_position_x = position.x;
			local_position_y = position.y;
			local_position_z = position.z;
		}

		//----TForm and Emission
		if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			tfx_vec3_t lerp_position = tfx_InterpolateVec3(tween, emitter.captured_position, emitter.world_position);
			tfx_vec3_t position_plus_handle = tfx_vec3_t(local_position_x, local_position_y, local_position_z) + emitter.handle;
			tfx_vec3_t pos = RotateVectorQuaternion(&emitter.rotation, position_plus_handle);
			local_position_x = lerp_position.x + pos.x * entry->overal_scale;
			local_position_y = lerp_position.y + pos.y * entry->overal_scale;
			local_position_z = lerp_position.z + pos.z * entry->overal_scale;
		}

		tween += entry->qty_step_size;
	}

}

void SpawnParticleIcosphere3d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 17 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;
	tfx_vec3_t half_emitter_size = emitter.emitter_size * .5f;
	tfxU32 sub_division = tfx__Min((tfxU32)properties.grid_points.x, 5);

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];
		float &local_position_z = entry->particle_data->position_z[index];

		local_position_x = local_position_y = local_position_z = 0;

		if (emitter.grid_coords.x >= tfxIcospherePoints[sub_division].current_size) {
			emitter.grid_coords.x = 0;
		}
		local_position_x = tfxIcospherePoints[sub_division][(tfxU32)emitter.grid_coords.x].x * half_emitter_size.x;
		local_position_y = tfxIcospherePoints[sub_division][(tfxU32)emitter.grid_coords.x].y * half_emitter_size.y;
		local_position_z = tfxIcospherePoints[sub_division][(tfxU32)emitter.grid_coords.x].z * half_emitter_size.z;
		emitter.grid_coords.x++;

		if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			tfx_vec3_t lerp_position = tfx_InterpolateVec3(tween, emitter.captured_position, emitter.world_position);
			tfx_vec3_t position_plus_handle = tfx_vec3_t(local_position_x, local_position_y, local_position_z) + emitter.handle;
			tfx_vec3_t pos = RotateVectorQuaternion(&emitter.rotation, position_plus_handle);
			local_position_x = lerp_position.x + pos.x * entry->overal_scale;
			local_position_y = lerp_position.y + pos.y * entry->overal_scale;
			local_position_z = lerp_position.z + pos.z * entry->overal_scale;
		}

		tween += entry->qty_step_size;
	}

}

void SpawnParticlePath3d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 26 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;
	tfx_vec3_t half_emitter_size = emitter.emitter_size * .5f;
	const tfx_vec3_t &grid_points = properties.grid_points;
	tfx_emitter_path_t *path = &pm.library->paths[emitter.path_attributes];
	float total_grid_points = (float)path->node_count - 3.f;
	float increment = 1.f / grid_points.x;
	float arc_size = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.arc_size, emitter.frame);
	float arc_offset = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.arc_offset, emitter.frame);
	float extrusion = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.extrusion, emitter.frame);
	tfx_vec3_t point;
	bool has_rotation = path->rotation_range > 0 || path->rotation_pitch != 0 || path->rotation_yaw != 0;

	float emission_pitch;
	float emission_yaw;
	float emission_angle_variation;
	float range;
	tfx_vec3_t velocity_direction;

	if (properties.emission_direction == tfxPathGradient) {
		emission_pitch = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.emission_pitch, emitter.frame);
		emission_yaw = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.emission_yaw, emitter.frame);
		emission_angle_variation = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.emission_range, emitter.frame);
		range = emission_angle_variation * .5f;
	}

	if (path->rotation_cycle_length > 0) {
		if ((emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid && emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_random) ||
			!(emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid)
			) {
			for (int qi = 0; qi != emitter.active_paths; ++qi) {
				int index = (qi + emitter.path_start_index) % path->maximum_active_paths;
				emitter.path_quaternions[qi].age += pm.frame_length;
				if (emitter.path_quaternions[qi].age >= path->rotation_cycle_length) {
					tfx_quaternion_t q = GetPathRotation(&random, path->rotation_range, path->rotation_pitch, path->rotation_yaw, ((path->flags & tfxPathFlags_rotation_range_yaw_only) > 0));
					emitter.path_quaternions[qi].quaternion = tfx__pack8bit_quaternion(q);
					emitter.path_quaternions[qi].age = 0.f;
					emitter.path_cycle_count--;
				}
			}
		}
	}

	tfxU32 qi = (emitter.path_start_index + emitter.last_path_index) % path->maximum_active_paths;
	TFX_ASSERT(qi < path->maximum_active_paths);

	if (path->rotation_stagger > 0 && emitter.path_stagger_counter >= path->rotation_stagger) {
		if (emitter.active_paths < path->maximum_active_paths && (emitter.path_cycle_count > 0 || path->maximum_paths == 0)) {
			emitter.last_path_index = emitter.active_paths;
			qi = (emitter.path_start_index + emitter.active_paths++) % path->maximum_active_paths;
			TFX_ASSERT(qi < path->maximum_active_paths);
			tfx_quaternion_t q = GetPathRotation(&random, path->rotation_range, path->rotation_pitch, path->rotation_yaw, ((path->flags & tfxPathFlags_rotation_range_yaw_only) > 0));
			emitter.path_quaternions[qi].quaternion = tfx__pack8bit_quaternion(q);
			emitter.path_quaternions[qi].cycles = 0;
			if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
				emitter.path_quaternions[qi].grid_coord = 0.f;
			}
			else if (!(emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
				emitter.path_quaternions[qi].grid_coord = total_grid_points - increment;
			}
			emitter.path_cycle_count--;
			emitter.path_stagger_counter = 0.f;
		}
	}

	int dead_paths = 0;
	int node = 0;
	float t = 0.f;

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];
		float &local_position_z = entry->particle_data->position_z[index];
		float &path_position = entry->particle_data->path_position[index];
		float &path_offset = entry->particle_data->path_offset[index];

		if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid && emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
			node = RandomRange(&random, path->node_count - 3);
			t = (float)RandomRange(&random, (int)grid_points.x) * increment;
			path_position = (float)node + t;
			point = tfx__catmull_rom_spline_3d_soa(path->node_soa.x, path->node_soa.y, path->node_soa.z, node, t);
		}
		else if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {
			float &grid_coord = emitter.path_quaternions[qi].grid_coord;
			bool new_path = false;
			if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
				grid_coord += increment;
				if (grid_coord >= total_grid_points) {
					grid_coord = 0.f;
					new_path = true;
				}
			}
			else if (!(emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
				grid_coord -= increment;
				if (grid_coord < 0) {
					grid_coord = total_grid_points - increment;
					new_path = true;
				}
			}
			if (new_path) {
				if (emitter.state_flags & tfxEmitterStateFlags_has_rotated_path && path->rotation_stagger == 0) {
					if (path->maximum_paths == 0 || emitter.path_cycle_count > 0) {
						tfx_quaternion_t q = GetPathRotation(&random, path->rotation_range, path->rotation_pitch, path->rotation_yaw, ((path->flags & tfxPathFlags_rotation_range_yaw_only) > 0));
						emitter.path_quaternions[qi].quaternion = tfx__pack8bit_quaternion(q);
						emitter.path_cycle_count--;
					}
					else {
						dead_paths++;
						emitter.path_quaternions[qi].cycles = tfxINVALID;
						entry->particle_data->flags[index] |= tfxParticleFlags_remove;
					}
				}
				else {
					emitter.path_quaternions[qi].cycles = tfxINVALID;
					dead_paths++;
				}
			}
			node = (int)grid_coord;
			t = grid_coord - (int)grid_coord;
			path_position = (float)node + t;
			point = tfx__catmull_rom_spline_3d_soa(path->node_soa.x, path->node_soa.y, path->node_soa.z, node, t);
		}
		else {
			node = RandomRange(&random, path->node_count - 3);
			t = GenerateRandom(&random);
			path_position = (float)node + t;

			point = tfx__catmull_rom_spline_3d_soa(path->node_soa.x, path->node_soa.y, path->node_soa.z, node, t);
		}

		if (path->extrusion_type == tfxExtrusionLinear) {
			if (properties.emission_direction == tfxPathGradient) {
				velocity_direction = tfx__catmull_rom_spline_gradient_3d_soa(&path->node_soa.x[node], &path->node_soa.y[node], &path->node_soa.z[node], t);
				velocity_direction = tfx__normalize_vec3_fast(&velocity_direction);
			}
			float radius = extrusion * .5f;
			path_offset = RandomRange(&random, -radius, radius);
			local_position_x = point.x + path_offset;
			local_position_x *= emitter.emitter_size.x;
			local_position_y = point.y * emitter.emitter_size.y;
			local_position_z = point.z * emitter.emitter_size.z;
		}
		else {
			path_offset = RandomRange(&random, arc_size) + arc_offset;
			float length_squared = point.x * point.x + point.z * point.z;
			float radius = length_squared == 0.f ? 0.f : 1.f / tfx__quake_sqrt(length_squared);
			float angle = atan2f(point.z, point.x) + path_offset;
			float rx = cosf(angle);
			float rz = sinf(angle);
			local_position_x = rx * radius * emitter.emitter_size.x;
			local_position_z = rz * radius * emitter.emitter_size.z;
			local_position_y = point.y * emitter.emitter_size.y;
			if (properties.emission_direction == tfxPathGradient) {
				velocity_direction = tfx__catmull_rom_spline_gradient_3d_soa(&path->node_soa.x[node], &path->node_soa.y[node], &path->node_soa.z[node], t);
				velocity_direction = tfx__normalize_vec3_fast(&velocity_direction);
				rx = cosf(path_offset);
				rz = sinf(path_offset);
				tfx_vec3_t v = velocity_direction;
				velocity_direction.x = v.x * rx - v.z * rz;
				velocity_direction.z = v.x * rz + v.z * rx;
			}
		}

		if (emitter.state_flags & tfxEmitterStateFlags_has_rotated_path && emitter.active_paths > 0) {
			if (emitter.path_quaternions[qi].cycles == tfxINVALID) {
				entry->particle_data->flags[index] |= tfxParticleFlags_remove;
				emitter.last_path_index++;
				emitter.last_path_index %= emitter.active_paths;
				qi = (emitter.path_start_index + emitter.last_path_index) % path->maximum_active_paths;
				TFX_ASSERT(qi < path->maximum_active_paths);
				continue;
			}
			tfx_quaternion_t q = tfx__unpack8bit_quaternion(emitter.path_quaternions[qi].quaternion);
			entry->particle_data->quaternion[index] = emitter.path_quaternions[qi].quaternion;
			tfx_vec3_t rp = { local_position_x, local_position_y, local_position_z };
			rp = RotateVectorQuaternion(&q, rp);
			local_position_x = rp.x;
			local_position_y = rp.y;
			local_position_z = rp.z;
			emitter.last_path_index++;
			emitter.last_path_index %= emitter.active_paths;
			qi = (emitter.path_start_index + emitter.last_path_index) % path->maximum_active_paths;
			TFX_ASSERT(qi < path->maximum_active_paths);
		}

		if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			tfx_vec3_t lerp_position = tfx_InterpolateVec3(tween, emitter.captured_position, emitter.world_position);
			tfx_vec3_t position_plus_handle = tfx_vec3_t(local_position_x, local_position_y, local_position_z) + emitter.handle;
			tfx_vec3_t pos = RotateVectorQuaternion(&emitter.rotation, position_plus_handle);
			local_position_x = lerp_position.x + pos.x * entry->overal_scale;
			local_position_y = lerp_position.y + pos.y * entry->overal_scale;
			local_position_z = lerp_position.z + pos.z * entry->overal_scale;
			if (properties.emission_direction == tfxPathGradient) {
				tfx_quaternion_t offset_quaternion = EulerToQuaternion(emission_yaw, emission_pitch, 0.f);
				tfx_vec3_t rotated_normal = RotateVectorQuaternion(&offset_quaternion, velocity_direction);
				rotated_normal = RotateVectorQuaternion(&emitter.rotation, rotated_normal);
				if (range != 0.f) {
					rotated_normal = RandomVectorInCone(&random, rotated_normal, range);
				}
				entry->particle_data->velocity_normal[index] = tfx__pack10bit_unsigned(&rotated_normal);
			}
		}
		else if (properties.emission_direction == tfxPathGradient) {
			tfx_quaternion_t offset_quaternion = EulerToQuaternion(emission_yaw, emission_pitch, 0.f);
			tfx_vec3_t rotated_normal = RotateVectorQuaternion(&offset_quaternion, velocity_direction);
			if (range != 0.f) {
				rotated_normal = RandomVectorInCone(&random, rotated_normal, range);
			}
			entry->particle_data->velocity_normal[index] = tfx__pack10bit_unsigned(&rotated_normal);
		}

		tween += entry->qty_step_size;
	}

	if (dead_paths > 0 && emitter.active_paths > 0) {
		tfxU32 offset = 0;
		for (int qi = emitter.active_paths - 1; qi >= 0; --qi) {
			int index = (emitter.path_start_index + qi) % path->maximum_active_paths;
			if (emitter.path_quaternions[index].cycles == tfxINVALID) {
				offset++;
			}
			else if (offset > 0) {
				tfxU32 next_index = (qi + offset + emitter.path_start_index) % path->maximum_active_paths;
				TFX_ASSERT(next_index < path->maximum_active_paths);
				emitter.path_quaternions[next_index] = emitter.path_quaternions[index];
			}
		}
		emitter.path_start_index = (emitter.path_start_index + offset) % path->maximum_active_paths;
		emitter.active_paths -= tfx__Min(emitter.active_paths, offset);
		emitter.last_path_index = 0;
	}

}

void SpawnParticleIcosphereRandom3d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 18 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;
	tfxU32 sub_division = tfx__Min((tfxU32)properties.grid_points.x, 5);
	tfx_vec3_t half_emitter_size = emitter.emitter_size * .5f;

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];
		float &local_position_z = entry->particle_data->position_z[index];

		local_position_x = local_position_y = local_position_z = 0;

		tfx_vec3_t half_emitter_size = emitter.emitter_size * .5f;
		int ico_point = RandomRange(&random, tfxIcospherePoints[sub_division].current_size);
		local_position_x = tfxIcospherePoints[sub_division][ico_point].x * half_emitter_size.x;
		local_position_y = tfxIcospherePoints[sub_division][ico_point].y * half_emitter_size.y;
		local_position_z = tfxIcospherePoints[sub_division][ico_point].z * half_emitter_size.z;
		if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			tfx_vec3_t lerp_position = tfx_InterpolateVec3(tween, emitter.captured_position, emitter.world_position);
			tfx_vec3_t position_plus_handle = tfx_vec3_t(local_position_x, local_position_y, local_position_z) + emitter.handle;
			tfx_vec3_t pos = RotateVectorQuaternion(&emitter.rotation, position_plus_handle);
			local_position_x = lerp_position.x + pos.x * entry->overal_scale;
			local_position_y = lerp_position.y + pos.y * entry->overal_scale;
			local_position_z = lerp_position.z + pos.z * entry->overal_scale;
		}

		tween += entry->qty_step_size;
	}

}

void SpawnParticleCylinder3d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	float tween = entry->tween;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 19 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;
	const tfx_vec3_t &grid_points = properties.grid_points;
	float arc_size = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.arc_size, emitter.frame);
	float arc_offset = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.arc_offset, emitter.frame);
	const float grid_segment_size_x = arc_size / tfxMax(grid_points.x, 1.f);
	const float grid_segment_size_y = emitter.emitter_size.y / tfxMax(grid_points.y - 1.f, 1.f);
	tfx_vec3_t half_emitter_size = emitter.emitter_size * .5f;

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];
		float &local_position_z = entry->particle_data->position_z[index];

		local_position_x = local_position_y = local_position_z = 0;

		if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_on_grid && !(emitter.property_flags & tfxEmitterPropertyFlags_fill_area)) {
			if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
				emitter.grid_coords.x = (float)RandomRange(&random, (tfxU32)grid_points.x);
				emitter.grid_coords.y = (float)RandomRange(&random, (tfxU32)grid_points.y);

				float th = emitter.grid_coords.x * grid_segment_size_x + arc_offset;
				local_position_x = cosf(th) * half_emitter_size.x;
				local_position_y = emitter.grid_coords.y * grid_segment_size_y;
				local_position_z = -sinf(th) * half_emitter_size.z;
			}
			else {
				if (emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					emitter.grid_coords.x--;
					if (emitter.grid_coords.x < 0.f) {
						emitter.grid_coords.x = grid_points.x - 1;
						emitter.grid_coords.y--;
						if (emitter.grid_coords.y < 0.f) {
							emitter.grid_coords.y = grid_points.y - 1;
						}
					}
				}

				float th = emitter.grid_coords.x * grid_segment_size_x + arc_offset;
				local_position_x = cosf(th) * half_emitter_size.x;
				local_position_y = emitter.grid_coords.y * grid_segment_size_y;
				local_position_z = -sinf(th) * half_emitter_size.z;

				if (!(emitter.property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					emitter.grid_coords.x++;
					if (emitter.grid_coords.x >= grid_points.x) {
						emitter.grid_coords.x = 0.f;
						emitter.grid_coords.y++;
						if (emitter.grid_coords.y >= grid_points.y) {
							emitter.grid_coords.y = 0.f;
						}
					}
				}
			}

		}
		else if (!(emitter.property_flags & tfxEmitterPropertyFlags_fill_area)) {
			float th = RandomRange(&random, arc_size) + arc_offset;

			local_position_x = cosf(th) * half_emitter_size.x;
			local_position_y = RandomRange(&random, emitter.emitter_size.y);
			local_position_z = -sinf(th) * half_emitter_size.z;
		}
		else {
			local_position_x = RandomRange(&random, 0.f, emitter.emitter_size.x);
			local_position_y = RandomRange(&random, 0.f, emitter.emitter_size.y);
			local_position_z = RandomRange(&random, 0.f, emitter.emitter_size.z);

			while ((powf(local_position_x - half_emitter_size.x, 2) / powf(half_emitter_size.x, 2)) + (powf(local_position_z - half_emitter_size.z, 2) / powf(half_emitter_size.z, 2)) > 1) {
				local_position_x = RandomRange(&random, 0.f, emitter.emitter_size.x);
				local_position_z = RandomRange(&random, 0.f, emitter.emitter_size.z);
			}
			local_position_x -= half_emitter_size.x;
			local_position_z -= half_emitter_size.z;
		}

		//----TForm and Emission
		if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			tfx_vec3_t lerp_position = tfx_InterpolateVec3(tween, emitter.captured_position, emitter.world_position);
			tfx_vec3_t position_plus_handle = tfx_vec3_t(local_position_x, local_position_y, local_position_z) + emitter.handle;
			tfx_vec3_t pos = RotateVectorQuaternion(&emitter.rotation, position_plus_handle);
			local_position_x = lerp_position.x + pos.x * entry->overal_scale;
			local_position_y = lerp_position.y + pos.y * entry->overal_scale;
			local_position_z = lerp_position.z + pos.z * entry->overal_scale;
		}

		tween += entry->qty_step_size;
	}

}

void SpawnParticleWeight(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_library_t *library = pm.library;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 20 + emitter.seed_index);
	float weight = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].base.weight, emitter.frame) * entry->parent_spawn_controls->weight;
	float weight_variation = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].variation.weight, emitter.frame) * entry->parent_spawn_controls->weight;

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &base_weight = entry->particle_data->base_weight[index];

		//----Weight
		if (weight) {
			base_weight = weight;
			if (weight_variation > 0) {
				base_weight += RandomRange(&random, -weight_variation, weight_variation);
			}
		}
		else {
			base_weight = 0;
		}
	}

}

void SpawnParticleVelocity(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 21 + emitter.seed_index);
	tfx_library_t *library = pm.library;

	float velocity = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].base.velocity, emitter.frame) * entry->parent_spawn_controls->velocity;
	float velocity_variation = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].variation.velocity, emitter.frame) * entry->parent_spawn_controls->velocity;

	if (entry->emission_type == tfxOtherEmitter) {
		emitter.grid_coords.x = emitter.grid_coords.y;
	}

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &base_velocity = entry->particle_data->base_velocity[index];

		//----Velocity
		base_velocity = velocity + RandomRange(&random, -velocity_variation, velocity_variation);

		if (entry->emission_type == tfxOtherEmitter) {
			int spawn_index = (int)emitter.grid_coords.x;
			spawn_index = GetCircularIndex(&pm.particle_location_buffers[emitter.spawn_locations_index], spawn_index);
			float age_lerp = pm.particle_location_arrays[emitter.spawn_locations_index].age[spawn_index];

			emitter.grid_coords.x++;
			if ((tfxU32)emitter.grid_coords.x >= pm.particle_location_buffers[emitter.spawn_locations_index].current_size) {
				emitter.grid_coords.x = 0.f;
			}
			float max_age = entry->particle_data->max_age[index];
			float factor = lookup_overtime_callback(&library->emitter_attributes[emitter.emitter_attributes].factor.velocity, age_lerp * max_age, max_age);
			base_velocity *= factor;
		}

	}

}

void SpawnParticleRoll(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 22 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;
	const tfxAngleSettingFlags angle_settings = properties.angle_settings;
	const float angle_roll_offset = properties.angle_offsets.roll;

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		float &roll = entry->particle_data->local_rotations_z[index];

		roll = 0;
		if (angle_settings & tfxAngleSettingFlags_random_roll) {
			roll = RandomRange(&random, angle_roll_offset);
		}
		else if (angle_settings & tfxAngleSettingFlags_specify_roll) {
			roll = angle_roll_offset;
		}
		else {
			roll = 0;
		}
	}

}

void SpawnParticleMicroUpdate2d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 23 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;
	const float splatter = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.splatter, emitter.frame) * entry->parent_spawn_controls->splatter;

	if (splatter) {
		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
			float &local_position_x = entry->particle_data->position_x[index];
			float &local_position_y = entry->particle_data->position_y[index];

			//----Splatter
			float splattertemp = splatter;
			float splatx = RandomRange(&random, -splatter, splatter);
			float splaty = RandomRange(&random, -splatter, splatter);

			while (GetDistance(0, 0, splatx, splaty) >= splattertemp && splattertemp > 0) {
				splatx = RandomRange(&random, -splatter, splatter);
				splaty = RandomRange(&random, -splatter, splatter);
			}

			if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
				local_position_x += splatx * entry->overal_scale;
				local_position_y += splaty * entry->overal_scale;
			}
			else {
				local_position_x += splatx;
				local_position_y += splaty;
			}
		}
	}

	tfx_library_t *library = pm.library;
	const float first_velocity_value = GetGraphFirstValue(&library->emitter_attributes[emitter.emitter_attributes].overtime.velocity);
	const float first_weight_value = GetGraphFirstValue(&library->emitter_attributes[emitter.emitter_attributes].overtime.weight);
	const tfxAngleSettingFlags angle_settings = properties.angle_settings;
	const float angle_roll_offset = properties.angle_offsets.roll;
	const tfx_emission_type emission_type = properties.emission_type;
	const tfxU32 layer = properties.layer;
	bool is_recording = (pm.flags & tfxParticleManagerFlags_recording_sprites) > 0 && (pm.flags & tfxParticleManagerFlags_using_uids) > 0;
	tfx_buffer_t &instance_buffer = !is_recording ? pm.instance_buffer : pm.instance_buffer_for_recording[pm.current_sprite_buffer][layer];

	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		const float base_weight = entry->particle_data->base_weight[index];
		float &roll = entry->particle_data->local_rotations_z[index];
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];
		float &direction = entry->particle_data->local_rotations_x[index];

		tfx_vec2_t sprite_transform_position;
		float sprite_transform_rotation;

		bool line = emitter.property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == tfxLine;

		if (!line && !(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			tfx__transform_particle_position(local_position_x, local_position_y, roll, &sprite_transform_position, &sprite_transform_rotation);
		}

		if (!line) {
			if (properties.emission_direction != tfxPathGradient) {
				direction = GetEmissionDirection2d(&pm, library, &random, emitter, tfx_vec2_t(local_position_x, local_position_y), sprite_transform_position) + GetGraphFirstValue(&library->emitter_attributes[emitter.emitter_attributes].overtime.direction);
			}
		}
		else {
			direction = 0;
		}

		if (line || emitter.property_flags & tfxEmitterPropertyFlags_relative_position) {
			tfx_vec2_t rotatevec = RotateVectorQuaternion2d(&emitter.rotation, tfx_vec2_t(local_position_x, local_position_y) + emitter.handle.xy());
		}

		if ((angle_settings & tfxAngleSettingFlags_align_roll || angle_settings & tfxAngleSettingFlags_align_with_emission) && !line) {
			tfx_vec2_t velocity_normal;
			//Why are we doing this?
			velocity_normal.x = sinf(direction);
			velocity_normal.y = -cosf(direction);
			roll = GetVectorAngle(velocity_normal.x, velocity_normal.y) + angle_roll_offset;
		}

		if (pm.flags & tfxParticleManagerFlags_ordered_by_age) {
			tfx_depth_index_t depth_id;
			depth_id.particle_id = tfx__make_particle_id(emitter.particles_index, index);
			depth_id.depth = entry->particle_data->age[index];
			entry->particle_data->depth_index[index] = entry->depth_index_start;
			(*entry->depth_indexes)[entry->depth_index_start] = depth_id;
			TFX_ASSERT(entry->depth_index_start < instance_buffer.current_size);
			entry->depth_index_start++;
		}
	}
}

void SpawnParticleMicroUpdate3d(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_spawn_work_entry_t *entry = static_cast<tfx_spawn_work_entry_t *>(data);
	tfx_random_t random = entry->random;
	tfx_particle_manager_t &pm = *entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[entry->emitter_index];
	AlterRandomSeed(&random, 24 + emitter.seed_index);
	const tfx_emitter_properties_t &properties = *entry->properties;
	const float splatter = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].properties.splatter, emitter.frame) * entry->parent_spawn_controls->splatter;

	if (splatter) {
		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
			float &local_position_x = entry->particle_data->position_x[index];
			float &local_position_y = entry->particle_data->position_y[index];
			float &local_position_z = entry->particle_data->position_z[index];

			//----Splatter
			if (splatter) {
				float splatx = RandomRange(&random, -splatter, splatter);
				float splaty = RandomRange(&random, -splatter, splatter);
				float splatz = RandomRange(&random, -splatter, splatter);

				while (powf(splatx / splatter, 2.f) + powf(splaty / splatter, 2.f) + powf(splatz / splatter, 2.f) > 1.f) {
					splatx = RandomRange(&random, -splatter, splatter);
					splaty = RandomRange(&random, -splatter, splatter);
					splatz = RandomRange(&random, -splatter, splatter);
				}

				if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
					local_position_x += splatx * entry->overal_scale;
					local_position_y += splaty * entry->overal_scale;
					local_position_z += splatz * entry->overal_scale;
				}
				else {
					local_position_x += splatx;
					local_position_y += splaty;
					local_position_z += splatz;
				}
			}

		}
	}


	tfx_library_t *library = pm.library;
	float emission_pitch = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].properties.emission_pitch, emitter.frame);
	float emission_yaw = lookup_callback(&library->emitter_attributes[emitter.emitter_attributes].properties.emission_yaw, emitter.frame);
	const float first_velocity_value = GetGraphFirstValue(&library->emitter_attributes[emitter.emitter_attributes].overtime.velocity);
	const float first_weight_value = GetGraphFirstValue(&library->emitter_attributes[emitter.emitter_attributes].overtime.weight);
	const tfx_emission_type emission_type = properties.emission_type;
	const bool line = emitter.property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == tfxLine;
	const float velocity_adjuster = lookup_callback(&pm.library->emitter_attributes[emitter.emitter_attributes].overtime.velocity_adjuster, emitter.frame);
	const tfxU32 layer = properties.layer;
	tfx_emission_direction emission_direction = library->emitter_properties[emitter.properties_index].emission_direction;
	bool ordered_effect = (entry->parent_property_flags & tfxEffectPropertyFlags_age_order) || (entry->parent_property_flags & tfxEffectPropertyFlags_depth_draw_order) > 0;
	bool is_recording = (pm.flags & tfxParticleManagerFlags_recording_sprites) > 0 && (pm.flags & tfxParticleManagerFlags_using_uids) > 0;
	tfx_buffer_t &instance_buffer = !is_recording ? pm.instance_buffer : pm.instance_buffer_for_recording[pm.current_sprite_buffer][layer];

	//Micro Update
	for (int i = 0; i != entry->amount_to_spawn; ++i) {
		tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[emitter.particles_index], entry->spawn_start_index + i);
		const float age = entry->particle_data->age[index];
		const float base_weight = entry->particle_data->base_weight[index];
		float &roll = entry->particle_data->local_rotations_z[index];
		float &local_position_x = entry->particle_data->position_x[index];
		float &local_position_y = entry->particle_data->position_y[index];
		float &local_position_z = entry->particle_data->position_z[index];
		float &captured_position_x = entry->particle_data->captured_position_x[index];
		float &captured_position_y = entry->particle_data->captured_position_y[index];
		float &captured_position_z = entry->particle_data->captured_position_z[index];
		tfxU32 &velocity_normal_packed = entry->particle_data->velocity_normal[index];
		const float base_velocity = entry->particle_data->base_velocity[index];

		tfx_vec3_t world_position;
		if (!line && !(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
			if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position) && !(emitter.property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
				world_position.x = local_position_x;
				world_position.y = local_position_y;
				world_position.z = local_position_z;
			}
			else {
				tfx_vec3_t position_plus_handle = tfx_vec3_t(local_position_x, local_position_y, local_position_z) + emitter.handle;
				tfx_vec3_t rotatevec = RotateVectorQuaternion(&emitter.rotation, tfx_vec3_t(local_position_x, local_position_y, local_position_z) + emitter.handle);
				world_position = emitter.world_position + rotatevec * entry->overal_scale;
			}
			captured_position_x = world_position.x;
			captured_position_y = world_position.y;
			captured_position_z = world_position.z;
		}

		tfx_vec3_t velocity_normal;
		if (emission_type == tfxPoint) {
			velocity_normal = GetEmissionDirection3d(&pm, library, &random, emitter, emission_pitch, emission_yaw, tfx_vec3_t(local_position_x, local_position_y, local_position_z), world_position);
			velocity_normal_packed = tfx__pack10bit_unsigned(&velocity_normal);
		}
		else if (emission_direction != tfxPathGradient && emission_direction != tfxOrbital) {
			if (emitter.property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == tfxLine) {
				velocity_normal_packed = tfxPACKED_Y_NORMAL_3D;
			}
			else if ((emission_type != tfxArea || (emission_type == tfxArea && emission_direction != tfxSurface))) {
				//----Velocity
				velocity_normal = GetEmissionDirection3d(&pm, library, &random, emitter, emission_pitch, emission_yaw, tfx_vec3_t(local_position_x, local_position_y, local_position_z), world_position);
				velocity_normal_packed = tfx__pack10bit_unsigned(&velocity_normal);
			}
		}
		else if (emission_direction == tfxOrbital) {
			if (!(emitter.property_flags & tfxEmitterPropertyFlags_relative_position)) {
				velocity_normal.z = local_position_x - emitter.world_position.x;
				velocity_normal.y = 0.f;
				velocity_normal.x = local_position_z - emitter.world_position.z * -1.f;
				velocity_normal = tfx__normalize_vec3_fast(&velocity_normal);
			}
			else {
				velocity_normal.z = emitter.handle.x + local_position_x;
				velocity_normal.y = 0.f;
				velocity_normal.x = (emitter.handle.z + local_position_z) * -1.f;
				velocity_normal = tfx__normalize_vec3_fast(&velocity_normal);
			}
			velocity_normal_packed = tfx__pack10bit_unsigned(&velocity_normal);
		}
		if (emitter.control_profile & tfxEmitterControlProfile_motion_randomness) {
			entry->particle_data->noise_offset[index] = 0;
		}
		if (pm.flags & tfxParticleManagerFlags_order_by_depth || (pm.flags & tfxParticleManagerFlags_use_effect_sprite_buffers && entry->root_effect_flags & tfxEffectPropertyFlags_depth_draw_order)) {
			tfx_depth_index_t depth_index;
			depth_index.particle_id = tfx__make_particle_id(emitter.particles_index, index);
			tfx_vec3_t world_minus_camera = tfx_vec3_t(local_position_x, local_position_y, local_position_z) - pm.camera_position;
			depth_index.depth = tfx__length_vec3_nosqr(&world_minus_camera);
			entry->particle_data->depth_index[index] = entry->depth_index_start;
			(*entry->depth_indexes)[entry->depth_index_start] = depth_index;
			entry->depth_index_start++;
		}
		else if (pm.flags & tfxParticleManagerFlags_ordered_by_age || (pm.flags & tfxParticleManagerFlags_use_effect_sprite_buffers && entry->root_effect_flags & tfxEffectPropertyFlags_age_order)) {
			tfx_depth_index_t depth_index;
			depth_index.particle_id = tfx__make_particle_id(emitter.particles_index, index);
			depth_index.depth = entry->particle_data->age[index];
			entry->particle_data->depth_index[index] = entry->depth_index_start;
			(*entry->depth_indexes)[entry->depth_index_start] = depth_index;
			entry->depth_index_start++;
		}
	}
}

void UpdateEmitterState(tfx_particle_manager_t *pm, tfx_emitter_state_t &emitter, tfxU32 parent_index, const tfx_parent_spawn_controls_t *parent_spawn_controls, tfx_spawn_work_entry_t *entry) {
	tfxPROFILE;

	tfx_library_t *library = pm->library;
	tfx_emitter_properties_t &properties = *entry->properties;
	emitter.bounding_box.min_corner.x = FLT_MAX;
	emitter.bounding_box.min_corner.y = FLT_MAX;
	emitter.bounding_box.min_corner.z = FLT_MAX;
	emitter.bounding_box.max_corner.x = -FLT_MAX;
	emitter.bounding_box.max_corner.y = -FLT_MAX;
	emitter.bounding_box.max_corner.z = -FLT_MAX;

	bool is_area = properties.emission_type != tfxPoint && properties.emission_type != tfxLine;

	emitter.emitter_size = { 0 };
	if (is_area) {
		emitter.emitter_size.y = LookupPrecise(&library->emitter_attributes[emitter.emitter_attributes].properties.emitter_height, emitter.age);
		emitter.emitter_size.x = LookupPrecise(&library->emitter_attributes[emitter.emitter_attributes].properties.emitter_width, emitter.age);
	}
	else if (properties.emission_type == tfxLine) {
		emitter.emitter_size.y = LookupPrecise(&library->emitter_attributes[emitter.emitter_attributes].properties.emitter_height, emitter.age);
	}

	if (emitter.property_flags & tfxEmitterPropertyFlags_effect_is_3d) {
		emitter.emitter_size.z = LookupPrecise(&library->emitter_attributes[emitter.emitter_attributes].properties.emitter_depth, emitter.age);
	}

	emitter.emitter_size *= pm->effects[parent_index].emitter_size;

	emitter.handle = { 0 };
	if (!(emitter.property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center)) {
		emitter.handle = properties.emitter_handle;
	}

}

void UpdateEffectState(tfx_particle_manager_t *pm, tfxU32 index) {
	tfxPROFILE;

	const float frame = pm->effects[index].frame;
	const float age = pm->effects[index].age;
	const tfxU32 global_attributes = pm->effects[index].global_attributes;
	const tfxU32 transform_attributes = pm->effects[index].transform_attributes;
	tfx_library_t *library = pm->effects[index].library;
	tfx_vec3_t &translation = pm->effects[index].translation;
	tfx_vec3_t &local_rotations = pm->effects[index].local_rotations;
	tfx_vec3_t &emitter_size = pm->effects[index].emitter_size;
	float &overal_scale = pm->effects[index].overal_scale;
	float &noise = pm->effects[index].noise;
	tfxEffectStateFlags &state_flags = pm->effects[index].state_flags;
	//float &stretch = pm->effects[index].stretch;

	//If this effect is a sub effect then the graph index will reference the global graphs for the root parent effect
	tfx_parent_spawn_controls_t &spawn_controls = pm->effects[index].spawn_controls;
	spawn_controls.life = lookup_callback(&library->global_graphs[global_attributes].life, frame);
	if (!(pm->effects[index].property_flags & tfxEffectPropertyFlags_global_uniform_size)) {
		spawn_controls.size_x = lookup_callback(&library->global_graphs[global_attributes].width, frame);
		spawn_controls.size_y = lookup_callback(&library->global_graphs[global_attributes].height, frame);
	}
	else {
		spawn_controls.size_x = lookup_callback(&library->global_graphs[global_attributes].width, frame);
		spawn_controls.size_y = spawn_controls.size_x;
	}
	spawn_controls.velocity = lookup_callback(&library->global_graphs[global_attributes].velocity, frame);
	noise = lookup_callback(&library->global_graphs[global_attributes].noise, frame);
	spawn_controls.spin = lookup_callback(&library->global_graphs[global_attributes].spin, frame);
	spawn_controls.pitch_spin = lookup_callback(&library->global_graphs[global_attributes].pitch_spin, frame);
	spawn_controls.yaw_spin = lookup_callback(&library->global_graphs[global_attributes].yaw_spin, frame);
	spawn_controls.intensity = lookup_callback(&library->global_graphs[global_attributes].intensity, frame);
	spawn_controls.splatter = lookup_callback(&library->global_graphs[global_attributes].splatter, frame);
	spawn_controls.weight = lookup_callback(&library->global_graphs[global_attributes].weight, frame);
	if (!(state_flags & tfxEffectStateFlags_override_size_multiplier)) {
		emitter_size.x = lookup_callback(&library->global_graphs[global_attributes].emitter_width, frame);
		emitter_size.y = lookup_callback(&library->global_graphs[global_attributes].emitter_height, frame);
		emitter_size.z = lookup_callback(&library->global_graphs[global_attributes].emitter_depth, frame);
	}
	//We don't want to scale twice when the sub effect is transformed, so the values here are set to 1. That means that the root effect will only control the global scale.
	overal_scale = state_flags & tfxEffectStateFlags_override_overal_scale ? overal_scale : lookup_callback(&library->global_graphs[global_attributes].overal_scale, frame);
	if (pm->effects[index].parent_particle_index == tfxINVALID) {
		if (!(state_flags & tfxEffectStateFlags_override_orientiation)) {
			local_rotations.roll = LookupPrecise(&library->transform_attributes[transform_attributes].roll, age);
			local_rotations.pitch = LookupPrecise(&library->transform_attributes[transform_attributes].pitch, age);
			local_rotations.yaw = LookupPrecise(&library->transform_attributes[transform_attributes].yaw, age);
		}
	}
	else {
		local_rotations.roll = 0.f;
		local_rotations.pitch = 0.f;
		local_rotations.yaw = 0.f;
	}
	pm->effects[index].stretch = lookup_callback(&library->global_graphs[global_attributes].stretch, frame);
	translation.x = LookupPrecise(&library->transform_attributes[transform_attributes].translation_x, age);
	translation.y = LookupPrecise(&library->transform_attributes[transform_attributes].translation_y, age);
	translation.z = LookupPrecise(&library->transform_attributes[transform_attributes].translation_z, age);

	if (pm->effects[index].update_callback) {
		pm->effects[index].update_callback(pm, index);
	}

}

void ControlParticleAge(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;
	tfx_particle_age_work_entry_t *work_entry = static_cast<tfx_particle_age_work_entry_t *>(data);
	tfx_particle_manager_t &pm = *work_entry->pm;
	tfx_emitter_state_t &emitter = pm.emitters[work_entry->emitter_index];
	tfx_effect_state_t &effect = pm.effects[emitter.root_index];
	tfx_vector_t<tfx_depth_index_t> *depth_indexes;
	const tfxWideInt single_shot_limit = tfxWideSetSinglei(work_entry->properties->single_shot_limit);
	const tfxU32 layer = work_entry->properties->layer;

	const tfxWideInt remove_flag = tfxWideSetSinglei(tfxParticleFlags_remove);
	const tfxWideInt capture_after_transform = tfxWideSetSinglei(tfxParticleFlags_capture_after_transform);
	const tfxWideInt remove = tfxWideSetSinglei(emitter.state_flags & tfxEmitterStateFlags_remove);
	const tfxWideInt single = tfxWideGreateri(tfxWideSetSinglei(emitter.property_flags & tfxEmitterPropertyFlags_single), tfxWideSetZeroi);
	const tfxWideInt not_single = tfxWideXOri(single, tfxWideSetSinglei(-1));
	const tfxWideInt wrap = tfxWideEqualsi(tfxWideSetSinglei(emitter.property_flags & tfxEmitterPropertyFlags_wrap_single_sprite), tfxWideSetZeroi);
	tfxWideInt state_flags_no_spawning = tfxWideGreateri(tfxWideOri(tfxWideSetSinglei(emitter.state_flags & tfxEmitterStateFlags_stop_spawning), tfxWideSetSinglei(work_entry->pm->flags & tfxParticleManagerFlags_disable_spawning)), tfxWideSetZeroi);
	if (emitter.property_flags & tfxEmitterPropertyFlags_wrap_single_sprite && pm.flags & tfxParticleManagerFlags_recording_sprites) {
		state_flags_no_spawning = tfxWideGreateri(tfxWideSetSinglei(emitter.property_flags & tfxEmitterPropertyFlags_wrap_single_sprite), tfxWideSetZeroi);
	}
	const tfxWideInt xor_state_flags_no_spawning = tfxWideXOri(state_flags_no_spawning, tfxWideSetSinglei(-1));

	tfx_particle_soa_t &bank = pm.particle_arrays[emitter.particles_index];

	for (int i = 0; i != work_entry->wide_end_index; i += tfxDataWidth) {
		tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i) / tfxDataWidth * tfxDataWidth;

		const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
		tfxWideFloat age = tfxWideLoad(&bank.age[index]);
		tfxWideInt single_loop_count = tfxWideLoadi((tfxWideIntLoader *)&bank.single_loop_count[index]);
		tfxWideInt flags = tfxWideLoadi((tfxWideIntLoader *)&bank.flags[index]);
		age = tfxWideAdd(age, pm.frame_length_wide);

		tfx__readbarrier;

		tfxWideInt expired = tfxWideCasti(tfxWideGreaterEqual(age, max_age));
		single_loop_count = tfxWideAddi(single_loop_count, tfxWideAndi(tfxWIDEONEi, expired));
		tfxWideInt loop_limit = tfxWideEqualsi(single_loop_count, single_shot_limit);
		tfxWideInt loop_age = tfxWideXOri(tfxWideAndi(tfxWideAndi(single, expired), xor_state_flags_no_spawning), tfxWideSetSinglei(-1));
		age = tfxWideAnd(age, tfxWideCast(loop_age));
		flags = tfxWideOri(flags, tfxWideAndi(remove_flag, tfxWideGreateri(remove, tfxWideSetZeroi)));
		flags = tfxWideOri(flags, tfxWideAndi(remove_flag, tfxWideAndi(not_single, expired)));
		flags = tfxWideOri(flags, tfxWideAndi(remove_flag, tfxWideAndi(tfxWideOri(tfxWideAndi(single, loop_limit), state_flags_no_spawning), expired)));
		flags = tfxWideOri(flags, tfxWideAndi(capture_after_transform, tfxWideAndi(expired, wrap)));

		tfxWideStore(&bank.age[index], age);
		tfxWideStorei((tfxWideIntLoader *)&bank.flags[index], flags);
		tfxWideStorei((tfxWideIntLoader *)&bank.single_loop_count[index], single_loop_count);
	}

	if (IsOrderedEffectState(&effect)) { //&& pm.flags & tfxParticleManagerFlags_use_effect_sprite_buffers) {
		depth_indexes = &effect.instance_data.depth_indexes[layer][effect.instance_data.current_depth_buffer_index[layer]];
	}
	//else {
		//depth_indexes = &pm.depth_indexes[layer][pm.current_depth_buffer_index[layer]];
	//}

	tfxU32 offset = 0;
	bool has_random_movement = (emitter.control_profile & tfxEmitterControlProfile_noise) + (emitter.control_profile & tfxEmitterControlProfile_motion_randomness) > 0;
	bool is_ordered = !(pm.flags & tfxParticleManagerFlags_unordered) || (pm.flags & tfxParticleManagerFlags_use_effect_sprite_buffers && (effect.effect_flags & tfxEffectPropertyFlags_depth_draw_order || effect.effect_flags & tfxEffectPropertyFlags_age_order));
	tfxU32 max_index = 0;
	for (int i = work_entry->start_index; i >= 0; --i) {
		const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i);
		tfxParticleFlags &flags = bank.flags[index];
		if (flags & tfxParticleFlags_remove) {
			offset++;
			if (flags & tfxParticleFlags_has_sub_effects) {
				FreePMParticleIndex(&pm, &bank.particle_index[index]);
			}
			if (is_ordered) {
				(*depth_indexes)[bank.depth_index[index]].particle_id = tfxINVALID;
			}
		}
		else if (offset > 0) {
			//Can we eliminate this? The only reason that we do the below is to fill gaps in the buffer when life variation is used. Life variation means that 
			//particles expire at different rates so we have to manually fill the wholes in the ring buffer. What if we "pretend" that each particle spawns with the 
			//maximum life possible but sample the graphs based on the actual life that it has. Then we expire the particle based on the maximum life value so 
			//we don't have to do the below. This means that we'd still be processing particles that have already expired (but their scaling could be set to 0 so they're
			//not drawn) but only needing to bump the start index of the ring buffer rather then the memory moving that's done below and is pretty slow.
			//This would be an interesting experiment to run and profile because it would be nice to not have to do the below. My intuition tells me that below is
			//slower then just processing extra particles each frame because memory work tends to be.
			tfxU32 next_index = GetCircularIndex(&work_entry->pm->particle_array_buffers[emitter.particles_index], i + offset);
			max_index = tfx__Max(max_index, next_index);
			if (flags & tfxParticleFlags_has_sub_effects) {
				pm.particle_indexes[bank.particle_index[index]] = tfx__make_particle_id(emitter.particles_index, next_index);
			}

			if (pm.flags & tfxParticleManagerFlags_order_by_depth || effect.effect_flags & tfxEffectPropertyFlags_depth_draw_order) {
				(*depth_indexes)[bank.depth_index[index]].particle_id = tfx__make_particle_id(emitter.particles_index, next_index);
			}
			else if (pm.flags & tfxParticleManagerFlags_ordered_by_age || effect.effect_flags & tfxEffectPropertyFlags_age_order) {
				(*depth_indexes)[bank.depth_index[index]].particle_id = tfx__make_particle_id(emitter.particles_index, next_index);
				(*depth_indexes)[bank.depth_index[index]].depth = bank.age[index];
			}

			bank.sprite_index[next_index] = bank.sprite_index[index];
			bank.depth_index[next_index] = bank.depth_index[index];
			bank.particle_index[next_index] = bank.particle_index[index];
			bank.uid[next_index] = bank.uid[index];
			bank.uid[index] = 0;
			bank.flags[next_index] = bank.flags[index];
			bank.age[next_index] = bank.age[index];
			bank.max_age[next_index] = bank.max_age[index];
			bank.position_x[next_index] = bank.position_x[index];
			bank.position_y[next_index] = bank.position_y[index];
			bank.captured_position_x[next_index] = bank.captured_position_x[index];
			bank.captured_position_y[next_index] = bank.captured_position_y[index];
			bank.local_rotations_z[next_index] = bank.local_rotations_z[index];
			bank.local_rotations_x[next_index] = bank.local_rotations_x[index];
			bank.velocity_normal[next_index] = bank.velocity_normal[index];
			bank.base_weight[next_index] = bank.base_weight[index];
			bank.base_velocity[next_index] = bank.base_velocity[index];
			bank.base_spin[next_index] = bank.base_spin[index];
			bank.intensity_factor[next_index] = bank.intensity_factor[index];
			if (emitter.property_flags & tfxEmitterPropertyFlags_effect_is_3d) {
				bank.local_rotations_y[next_index] = bank.local_rotations_y[index];
				bank.position_z[next_index] = bank.position_z[index];
				bank.captured_position_z[next_index] = bank.captured_position_z[index];
			}
			if (emitter.property_flags & tfxEmitterPropertyFlags_effect_is_3d && emitter.state_flags & tfxEmitterStateFlags_can_spin_pitch_and_yaw) {
				bank.base_pitch_spin[next_index] = bank.base_pitch_spin[index];
				bank.base_yaw_spin[next_index] = bank.base_yaw_spin[index];
			}
			if (has_random_movement) {
				bank.noise_offset[next_index] = bank.noise_offset[index];
				bank.noise_resolution[next_index] = bank.noise_resolution[index];
			}
			if (emitter.state_flags & tfxEmitterStateFlags_has_path) {
				bank.path_position[next_index] = bank.path_position[index];
				bank.path_offset[next_index] = bank.path_offset[index];
			}
			if (emitter.state_flags & tfxEmitterStateFlags_has_rotated_path) {
				bank.quaternion[next_index] = bank.quaternion[index];
			}
			bank.color[next_index] = bank.color[index];
			bank.image_frame[next_index] = bank.image_frame[index];
			bank.base_size_x[next_index] = bank.base_size_x[index];
			bank.base_size_y[next_index] = bank.base_size_y[index];
			bank.single_loop_count[next_index] = bank.single_loop_count[index];
		}
	}

	if (offset) {
		Bump(&work_entry->pm->particle_array_buffers[emitter.particles_index], offset);
	}
}

void ControlParticles(tfx_work_queue_t *queue, void *data) {
	tfxPROFILE;

	tfx_control_work_entry_t *work_entry = static_cast<tfx_control_work_entry_t *>(data);

	tfx_particle_manager_t *pm = work_entry->pm;
	tfx_library_t *library = pm->library;
	tfx_emitter_state_t &emitter = pm->emitters[work_entry->emitter_index];
	tfx_emitter_properties_t &properties = *work_entry->properties;

	//-------------------------------------------------------
	//Controll what the particle does over the course of
	//it's lifetime
	//-------------------------------------------------------

	work_entry->graphs = &library->emitter_attributes[emitter.emitter_attributes].overtime;
	//work_entry->c.particle_update_callback = e.particle_update_callback;
	//work_entry->c.user_data = e.user_data;

	tfx_soa_buffer_t &buffer = pm->particle_array_buffers[emitter.particles_index];
	int offset = 0;
	tfxU32 amount_to_update = buffer.current_size;
	if (emitter.sprites_count < buffer.current_size) {
		amount_to_update = emitter.sprites_count;
	}

	work_entry->layer = properties.layer;
	tfx_effect_instance_data_t &instance_data = pm->effects[emitter.root_index].instance_data;
	work_entry->cumulative_index_point = (pm->flags & tfxParticleManagerFlags_use_effect_sprite_buffers) ? instance_data.cumulative_index_point[work_entry->layer] : pm->cumulative_index_point[work_entry->layer];
	work_entry->effect_instance_offset = (pm->flags & tfxParticleManagerFlags_use_effect_sprite_buffers) ? instance_data.instance_start_index : 0;
	work_entry->sprites_index = emitter.sprites_index + work_entry->start_index + work_entry->cumulative_index_point + work_entry->effect_instance_offset;
	work_entry->sprite_buffer_end_index = work_entry->sprites_index + (work_entry->end_index - work_entry->start_index);
	tfx_effect_instance_data_t &sprites = pm->effects[emitter.root_index].instance_data;
	work_entry->depth_indexes = &sprites.depth_indexes[work_entry->layer][sprites.current_depth_buffer_index[work_entry->layer]];
	work_entry->overal_scale = pm->effects[emitter.parent_index].overal_scale;
	work_entry->path = emitter.path_attributes != tfxINVALID ? &pm->library->paths[emitter.path_attributes] : nullptr;
	work_entry->node_count = work_entry->path ? work_entry->path->node_count : 0.f;
	work_entry->sample_path_life = work_entry->path ? work_entry->properties->emission_type == tfxPath && (emitter.property_flags & tfxEmitterPropertyFlags_alt_velocity_lifetime_sampling) > 0 : false;

	if (emitter.property_flags & tfxEmitterPropertyFlags_spawn_location_source && emitter.spawn_locations_index != tfxINVALID) {
		bool grew;
		tfx_soa_buffer_t &particle_buffer = pm->particle_array_buffers[emitter.particles_index];
		tfx_soa_buffer_t &spawn_point_buffer = pm->particle_location_buffers[emitter.spawn_locations_index];
		if (buffer.current_size > spawn_point_buffer.current_size) {
			AddRows(&pm->particle_location_buffers[emitter.spawn_locations_index], buffer.current_size - spawn_point_buffer.current_size, true, grew);
		}
		spawn_point_buffer.current_size = particle_buffer.current_size;
		TFX_ASSERT(spawn_point_buffer.current_size < spawn_point_buffer.capacity);
		spawn_point_buffer.start_index = particle_buffer.start_index;
	}

	if (amount_to_update > 0) {
		if (pm->flags & tfxParticleManagerFlags_recording_sprites && pm->flags & tfxParticleManagerFlags_using_uids) {
			ControlParticleUID(&pm->work_queue, work_entry);
		}
		if (pm->flags & tfxParticleManagerFlags_3d_effects && emitter.control_profile & tfxEmitterControlProfile_path && emitter.control_profile & tfxEmitterControlProfile_edge_traversal) {
			ControlParticlePositionPath3d(&pm->work_queue, work_entry);
			ControlParticleTransform3d(&pm->work_queue, work_entry);
		}
		else if (emitter.control_profile & tfxEmitterControlProfile_path && emitter.control_profile & tfxEmitterControlProfile_edge_traversal) {
			ControlParticlePositionPath2d(&pm->work_queue, work_entry);
			ControlParticleTransform2d(&pm->work_queue, work_entry);
		}
		else if (pm->flags & tfxParticleManagerFlags_3d_effects) {
			if (emitter.control_profile & tfxEmitterControlProfile_noise && emitter.control_profile & tfxEmitterControlProfile_orbital) {
				ControlParticlePosition3dNoiseOrbital(&pm->work_queue, work_entry);
			}
			else if (emitter.control_profile & tfxEmitterControlProfile_noise) {
				ControlParticlePosition3dNoise(&pm->work_queue, work_entry);
			}
			else if (emitter.control_profile & tfxEmitterControlProfile_motion_randomness && emitter.control_profile & tfxEmitterControlProfile_orbital) {
				ControlParticlePosition3dMotionRandomnessOrbital(&pm->work_queue, work_entry);
			}
			else if (emitter.control_profile & tfxEmitterControlProfile_motion_randomness) {
				ControlParticlePosition3dMotionRandomness(&pm->work_queue, work_entry);
			}
			else if (emitter.control_profile & tfxEmitterControlProfile_orbital) {
				ControlParticlePosition3dOrbital(&pm->work_queue, work_entry);
			}
			else {
				ControlParticlePosition3dBasic(&pm->work_queue, work_entry);
			}
			if (emitter.control_profile & tfxEmitterControlProfile_edge_kill && emitter.control_profile & tfxEmitterControlProfile_edge_traversal) {
				ControlParticleLineBehaviourKill(&pm->work_queue, work_entry);
			}
			else if (emitter.control_profile & tfxEmitterControlProfile_edge_loop && emitter.control_profile & tfxEmitterControlProfile_edge_traversal) {
				ControlParticleLineBehaviourLoop(&pm->work_queue, work_entry);
			}
			ControlParticleTransform3d(&pm->work_queue, work_entry);
		}
		else {
			ControlParticlePosition2d(&pm->work_queue, work_entry);
		}
		if (pm->flags & tfxParticleManagerFlags_3d_effects && emitter.state_flags & tfxEmitterStateFlags_can_spin_pitch_and_yaw) {
			ControlParticleSpin3d(&pm->work_queue, work_entry);
		}
		else {
			ControlParticleSpin(&pm->work_queue, work_entry);
		}
		ControlParticleSize(&pm->work_queue, work_entry);
		ControlParticleColor(&pm->work_queue, work_entry);
		ControlParticleImageFrame(&pm->work_queue, work_entry);
	}
}

void TransformEffector2d(tfx_vec3_t *world_rotations, tfx_vec3_t *local_rotations, tfx_vec3_t *world_position, tfx_vec3_t *local_position, tfx_quaternion_t *q, tfx_sprite_transform2d_t *parent, bool relative_position, bool relative_angle) {

	if (relative_position) {
		world_rotations->roll = parent->rotation + local_rotations->roll;
		*world_position = parent->position;
	}
	else {
		*world_position = *local_position;
		world_rotations->roll = local_rotations->roll;
	}

	ToQuaternion2d(q, world_rotations->roll);
}

void TransformEffector3d(tfx_vec3_t *world_rotations, tfx_vec3_t *local_rotations, tfx_vec3_t *world_position, tfx_vec3_t *local_position, tfx_quaternion_t *q, tfx_sprite_transform3d_t *parent, bool relative_position, bool relative_angle) {

	if (relative_position) {
		*world_rotations = parent->rotations + *local_rotations;
		*world_position = parent->position;
	}
	else {
		*world_position = *local_position;
		*world_rotations = *local_rotations;
	}

	*q = EulerToQuaternion(world_rotations->pitch, world_rotations->yaw, world_rotations->roll);
}

const tfxU32 tfxPROFILE_COUNT = __COUNTER__;
tfxU32 tfxCurrentSnapshot = 0;
tfx_profile_t tfxProfileArray[tfxPROFILE_COUNT];
int tfxNumberOfThreadsInAdditionToMain = 0;
bool tfxThreadUsage[32];

tfx_storage_t *tfxStore = 0;
tfx_allocator *tfxMemoryAllocator = 0;

void InitialiseTimelineFXMemory(size_t memory_pool_size) {
	if (tfxMemoryAllocator) return;
	void *memory_pool = tfxALLOCATE_POOL(memory_pool_size);
	TFX_ASSERT(memory_pool);    //unable to allocate initial memory pool
	tfxMemoryAllocator = tfx_InitialiseAllocatorWithPool(memory_pool, memory_pool_size, &tfxMemoryAllocator);
}

bool tfxInitialiseThreads(tfx_queue_processor_t *thread_queues) {
	//todo: create a function to close all the threads 
	for (int thread_index = 0; thread_index < tfxNumberOfThreadsInAdditionToMain; ++thread_index) {
		tfxStore->threads[thread_index] = std::thread(tfxThreadWorker, &tfxStore->tfxThreadQueues);
		tfxStore->thread_count++;
	}
	return true;
}

void tfxEndThread(tfx_work_queue_t *queue, void *data) {
	return;
}

//Passing a max_threads value of 0 or 1 will make timeline fx run in single threaded mode. 2 or more will be multithreaded.
//max_threads includes the main thread so for example if you set it to 4 then there will be the main thread plus an additional 3 threads.
void InitialiseTimelineFX(int max_threads, size_t memory_pool_size) {
	if (!tfxMemoryAllocator) {
		void *memory_pool = tfxALLOCATE_POOL(memory_pool_size);
		TFX_ASSERT(memory_pool);    //unable to allocate initial memory pool
		tfxMemoryAllocator = tfx_InitialiseAllocatorWithPool(memory_pool, memory_pool_size, &tfxMemoryAllocator);
	}
    tfx_storage_t store{};
	tfxStore = (tfx_storage_t *)tfx_AllocateAligned(tfxMemoryAllocator, sizeof(tfx_storage_t), 16);
    memcpy(tfxStore, &store, sizeof(tfx_storage_t));
	tfxStore->default_memory_pool_size = memory_pool_size;
	tfxStore->memory_pools[0] = (tfx_pool *)((char *)tfx__allocator_first_block(tfxMemoryAllocator) + tfx__POINTER_SIZE);
	tfxStore->memory_pool_count = 1;

	tfxNumberOfThreadsInAdditionToMain = max_threads = tfxMin(max_threads - 1 < 0 ? 0 : max_threads - 1, (int)std::thread::hardware_concurrency() - 1);
	lookup_callback = LookupFast;
	lookup_overtime_callback = LookupFastOvertime;
    
    std::unique_lock<std::mutex> lock(tfxStore->tfxThreadQueues.mutex);
	tfxInitialiseThreads(&tfxStore->tfxThreadQueues);
}

void EndTimelineFX() {
	tfxStore->tfxThreadQueues.end_all_threads = true;
	tfx__writebarrier;
	tfx_work_queue_t end_queue{};
	tfxU32 thread_count = tfxStore->thread_count;
	while (thread_count > 0) {
		tfxAddWorkQueueEntry(&end_queue, nullptr, tfxEndThread);
		tfxCompleteAllWork(&end_queue);
		thread_count--;
	}
	for (int i = 0; i != tfxStore->thread_count; ++i) {
		tfxStore->threads[i].join();
	}
	int pool_count = tfxStore->memory_pool_count;
	tfx__unlock_thread_access(tfxMemoryAllocator);
	for (int i = 0; i != 6; ++i) {
		tfxIcospherePoints[i].free();
	}
	for (int i = pool_count - 1; i != 0; --i) {
		free(tfxStore->memory_pools[i]);
	}
	free(tfxMemoryAllocator);
}

void SetColorRampFormat(tfx_color_ramp_format color_format) {
	tfxStore->color_ramp_format = color_format;
}

tfx_pool_stats_t CreateMemorySnapshot(tfx_header *first_block) {
	tfx_pool_stats_t stats = { 0 };
	tfx_header *current_block = first_block;
	while (!tfx__is_last_block_in_pool(current_block)) {
		if (tfx__is_free_block(current_block)) {
			stats.free_blocks++;
			stats.free_size += tfx__block_size(current_block);
		}
		else {
			stats.used_blocks++;
			stats.used_size += tfx__block_size(current_block);
		}
		current_block = tfx__next_physical_block(current_block);
	}
	if (tfx__is_free_block(current_block)) {
		stats.free_blocks++;
		stats.free_size += tfx__block_size(current_block);
	}
	else {
		stats.used_blocks++;
		stats.used_size += tfx__block_size(current_block);
	}
	return stats;
}

tfx_profile_tag_t::tfx_profile_tag_t(tfxU32 id, const char *name) {
	profile = tfxProfileArray + id;
	profile->name = name;
	snapshot = profile->snapshots + tfxCurrentSnapshot;
	start_time = tfx_Microsecs();
	start_cycles = tfx__rdtsc();
	tfx_AtomicAdd32(&snapshot->hit_count, 1);
}

void InitSpriteData3dSoACompression(tfx_soa_buffer_t *buffer, tfx_sprite_data_soa_t *soa, tfxU32 reserve_amount) {
	AddStructArray(buffer, sizeof(tfx_unique_sprite_id_t), offsetof(tfx_sprite_data_soa_t, uid));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_sprite_data_soa_t, lerp_offset));
	AddStructArray(buffer, sizeof(tfx_3d_instance_t), offsetof(tfx_sprite_data_soa_t, billboard_instance));
	FinishSoABufferSetup(buffer, soa, reserve_amount, 16);
}

void InitSpriteData3dSoA(tfx_soa_buffer_t *buffer, tfx_sprite_data_soa_t *soa, tfxU32 reserve_amount) {
	AddStructArray(buffer, sizeof(tfx_unique_sprite_id_t), offsetof(tfx_sprite_data_soa_t, uid));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_sprite_data_soa_t, lerp_offset));
	AddStructArray(buffer, sizeof(tfx_3d_instance_t), offsetof(tfx_sprite_data_soa_t, billboard_instance));
	FinishSoABufferSetup(buffer, soa, reserve_amount, 16);
}

void InitSpriteData2dSoACompression(tfx_soa_buffer_t *buffer, tfx_sprite_data_soa_t *soa, tfxU32 reserve_amount) {
	AddStructArray(buffer, sizeof(tfx_unique_sprite_id_t), offsetof(tfx_sprite_data_soa_t, uid));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_sprite_data_soa_t, lerp_offset));
	AddStructArray(buffer, sizeof(tfx_2d_instance_t), offsetof(tfx_sprite_data_soa_t, sprite_instance));
	FinishSoABufferSetup(buffer, soa, reserve_amount, 16);
}

void InitSpriteData2dSoA(tfx_soa_buffer_t *buffer, tfx_sprite_data_soa_t *soa, tfxU32 reserve_amount) {
	AddStructArray(buffer, sizeof(tfx_unique_sprite_id_t), offsetof(tfx_sprite_data_soa_t, uid));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_sprite_data_soa_t, lerp_offset));
	AddStructArray(buffer, sizeof(tfx_2d_instance_t), offsetof(tfx_sprite_data_soa_t, sprite_instance));
	FinishSoABufferSetup(buffer, soa, reserve_amount, 16);
}

void InitParticleSoA2d(tfx_soa_buffer_t *buffer, tfx_particle_soa_t *soa, tfxU32 reserve_amount, tfxEmitterControlProfileFlags control_profile) {
	AddStructArray(buffer, sizeof(tfxU32), offsetof(tfx_particle_soa_t, uid));
	AddStructArray(buffer, sizeof(tfxU32), offsetof(tfx_particle_soa_t, sprite_index));
	AddStructArray(buffer, sizeof(tfxParticleID), offsetof(tfx_particle_soa_t, particle_index));
	AddStructArray(buffer, sizeof(tfxParticleFlags), offsetof(tfx_particle_soa_t, flags));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, age));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, max_age));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, position_x));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, position_y));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, captured_position_x));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, captured_position_y));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, local_rotations_x));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, local_rotations_z));
	AddStructArray(buffer, sizeof(tfxU32), offsetof(tfx_particle_soa_t, velocity_normal));
	AddStructArray(buffer, sizeof(tfxU32), offsetof(tfx_particle_soa_t, depth_index));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, base_weight));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, base_velocity));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, base_spin));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, intensity_factor));
	if ((control_profile & tfxEmitterControlProfile_noise) + (control_profile & tfxEmitterControlProfile_motion_randomness) > 0) {
		AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, noise_offset));
		AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, noise_resolution));
	}
	if (control_profile & tfxEmitterControlProfile_path) {
		AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, path_position));
		AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, path_offset));
	}
	if (tfxEmitterStateFlags_has_rotated_path) {
		AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, quaternion));
	}
	AddStructArray(buffer, sizeof(tfx_rgba8_t), offsetof(tfx_particle_soa_t, color));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, image_frame));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, base_size_x));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, base_size_y));
	AddStructArray(buffer, sizeof(tfxU32), offsetof(tfx_particle_soa_t, single_loop_count));
	FinishSoABufferSetup(buffer, soa, reserve_amount, 16);
}

void InitParticleSoA3d(tfx_soa_buffer_t *buffer, tfx_particle_soa_t *soa, tfxU32 reserve_amount, tfxEmitterControlProfileFlags control_profile) {
	AddStructArray(buffer, sizeof(tfxU32), offsetof(tfx_particle_soa_t, uid));
	AddStructArray(buffer, sizeof(tfxU32), offsetof(tfx_particle_soa_t, sprite_index));
	AddStructArray(buffer, sizeof(tfxParticleID), offsetof(tfx_particle_soa_t, particle_index));
	AddStructArray(buffer, sizeof(tfxParticleFlags), offsetof(tfx_particle_soa_t, flags));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, age));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, max_age));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, position_x));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, position_y));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, position_z));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, captured_position_x));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, captured_position_y));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, captured_position_z));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, local_rotations_x));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, local_rotations_y));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, local_rotations_z));
	AddStructArray(buffer, sizeof(tfxU32), offsetof(tfx_particle_soa_t, velocity_normal));
	AddStructArray(buffer, sizeof(tfxU32), offsetof(tfx_particle_soa_t, depth_index));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, base_weight));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, base_velocity));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, base_spin));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, base_pitch_spin));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, base_yaw_spin));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, intensity_factor));
	if ((control_profile & tfxEmitterControlProfile_noise) + (control_profile & tfxEmitterControlProfile_motion_randomness) > 0) {
		AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, noise_offset));
		AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, noise_resolution));
	}
	if (control_profile & tfxEmitterControlProfile_path) {
		AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, path_position));
		AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, path_offset));
	}
	if (tfxEmitterStateFlags_has_rotated_path) {
		AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, quaternion));
	}
	AddStructArray(buffer, sizeof(tfx_rgba8_t), offsetof(tfx_particle_soa_t, color));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, image_frame));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, base_size_x));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_particle_soa_t, base_size_y));
	AddStructArray(buffer, sizeof(tfxU32), offsetof(tfx_particle_soa_t, single_loop_count));
	FinishSoABufferSetup(buffer, soa, reserve_amount, 16);
}

void InitPathsSoA2d(tfx_soa_buffer_t *buffer, tfx_path_nodes_soa_t *soa, tfxU32 reserve_amount) {
	AddStructArray(buffer, sizeof(float), offsetof(tfx_path_nodes_soa_t, x));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_path_nodes_soa_t, y));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_path_nodes_soa_t, length));
	FinishSoABufferSetup(buffer, soa, reserve_amount, 16);
}

void InitPathsSoA3d(tfx_soa_buffer_t *buffer, tfx_path_nodes_soa_t *soa, tfxU32 reserve_amount) {
	AddStructArray(buffer, sizeof(float), offsetof(tfx_path_nodes_soa_t, x));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_path_nodes_soa_t, y));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_path_nodes_soa_t, z));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_path_nodes_soa_t, length));
	FinishSoABufferSetup(buffer, soa, reserve_amount, 16);
}

void InitParticleLocationSoA3d(tfx_soa_buffer_t *buffer, tfx_spawn_points_soa_t *soa, tfxU32 reserve_amount) {
	AddStructArray(buffer, sizeof(float), offsetof(tfx_spawn_points_soa_t, position_x));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_spawn_points_soa_t, position_y));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_spawn_points_soa_t, position_z));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_spawn_points_soa_t, captured_position_x));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_spawn_points_soa_t, captured_position_y));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_spawn_points_soa_t, captured_position_z));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_spawn_points_soa_t, age));
	FinishSoABufferSetup(buffer, soa, reserve_amount, 16);
}

void InitParticleLocationSoA2d(tfx_soa_buffer_t *buffer, tfx_spawn_points_soa_t *soa, tfxU32 reserve_amount) {
	AddStructArray(buffer, sizeof(float), offsetof(tfx_spawn_points_soa_t, position_x));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_spawn_points_soa_t, position_y));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_spawn_points_soa_t, captured_position_x));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_spawn_points_soa_t, captured_position_y));
	AddStructArray(buffer, sizeof(float), offsetof(tfx_spawn_points_soa_t, age));
	FinishSoABufferSetup(buffer, soa, reserve_amount, 16);
}

void InitEmitterProperites(tfx_emitter_properties_t *properties) {
	properties->angle_offsets = { 0.f, 0.f, tfx360Radians };
	properties->image = nullptr;
	properties->image_handle = tfx_vec2_t();
	properties->spawn_amount = 1;
	properties->single_shot_limit = 0;
	properties->emission_type = tfx_emission_type::tfxPoint;
	properties->billboard_option = tfxBillboarding_align_to_camera;
	properties->vector_align_type = tfxVectorAlignType_motion;
	properties->emission_direction = tfx_emission_direction::tfxOutwards;
	properties->grid_points = { 10.f, 10.f, 10.f };
	properties->emitter_handle = { 0.f, 0.f, 0.f };
	properties->end_behaviour = tfx_line_traversal_end_behaviour::tfxLoop;
	properties->loop_length = 0.f;
	properties->layer = 0;
	properties->image_hash = 1;
	properties->start_frame = 0;
	properties->end_frame = 0;
	properties->frame_rate = 30.f;
	properties->angle_settings = tfxAngleSettingFlags_random_roll | tfxAngleSettingFlags_specify_pitch | tfxAngleSettingFlags_specify_yaw;
	properties->delay_spawning = 0.f;
	properties->noise_base_offset_range = 1000.f;
	properties->animation_property_index = tfxINVALID;
	properties->paired_emitter_hash = 0;
}

//Use with care, no checks for out of bounds
void CopyEmitterProperites(tfx_emitter_properties_t *from_properties, tfx_emitter_properties_t *to_properties) {
	*to_properties = *from_properties;
}

void FreeSpriteData(tfx_sprite_data_t *sprite_data) {
	if (sprite_data->compressed_sprites_buffer.data == sprite_data->real_time_sprites_buffer.data) {
		FreeSoABuffer(&sprite_data->real_time_sprites_buffer);
		sprite_data->normal.frame_meta.free_all();
		sprite_data->compressed_sprites_buffer = tfx_soa_buffer_t();
	}
	else {
		FreeSoABuffer(&sprite_data->compressed_sprites_buffer);
		FreeSoABuffer(&sprite_data->real_time_sprites_buffer);
		sprite_data->normal.frame_meta.free_all();
		sprite_data->compressed.frame_meta.free_all();
	}
}

bool ValidEffectID(tfx_particle_manager_t *pm, tfxEffectID id) {
	return id != tfxINVALID && pm->effects.capacity > id;
}

tfx_particle_manager_info_t CreateParticleManagerInfo(tfx_particle_manager_setup setup) {
	tfx_particle_manager_info_t info;
	info.max_particles = 5000;
	switch (setup) {
	case tfxParticleManagerSetup_2d_unordered:
		break;
	case tfxParticleManagerSetup_2d_ordered_by_age:
		info.order_mode = tfxParticleManagerMode_ordered_by_age;
		break;
	case tfxParticleManagerSetup_3d_ordered_by_age:
		info.order_mode = tfxParticleManagerMode_ordered_by_age;
		info.is_3d = true;
		break;
	case tfxParticleManagerSetup_2d_group_sprites_by_effect:
		info.group_sprites_by_effect = true;
		info.auto_order_effects = true;
		break;
	case tfxParticleManagerSetup_3d_group_sprites_by_effect:
		info.group_sprites_by_effect = true;
		info.auto_order_effects = true;
		info.is_3d = true;
		break;
	case tfxParticleManagerSetup_3d_unordered:
		info.is_3d = true;
		break;
	case tfxParticleManagerSetup_3d_ordered_by_depth:
		info.order_mode = tfxParticleManagerMode_ordered_by_depth;
		info.is_3d = true;
		break;
	case tfxParticleManagerSetup_3d_ordered_by_depth_guaranteed:
		info.order_mode = tfxParticleManagerMode_ordered_by_depth_guaranteed;
		info.is_3d = true;
		break;
	}

	return info;
}

void InitCommonParticleManager(tfx_particle_manager_t *pm, tfx_library_t *library, tfxU32 max_particles, unsigned int effects_limit, tfx_particle_manager_mode mode, bool double_buffered_sprites, bool dynamic_sprite_allocation, bool group_sprites_by_effect, tfxU32 mt_batch_size) {
	pm->random = NewRandom(tfx_Millisecs());
	pm->threaded_random = NewRandom(tfx_Millisecs());
	pm->max_effects = effects_limit;
	pm->mt_batch_size = mt_batch_size;
	tfxInitialiseWorkQueue(&pm->work_queue);
	pm->library = library;

	if (pm->particle_arrays.bucket_list.current_size == 0) {
		//todo need to be able to adjust the arena size
		pm->particle_arrays = tfxCreateBucketArray<tfx_particle_soa_t>(32);
	}

	pm->flags = 0;
	if (mode == tfxParticleManagerMode_unordered)
		pm->flags = tfxParticleManagerFlags_unordered;
	else if (mode == tfxParticleManagerMode_ordered_by_age)
		pm->flags = tfxParticleManagerFlags_ordered_by_age;
	else if (mode == tfxParticleManagerMode_ordered_by_depth)
		pm->flags = tfxParticleManagerFlags_order_by_depth;
	else if (mode == tfxParticleManagerMode_ordered_by_depth_guaranteed)
		pm->flags = tfxParticleManagerFlags_order_by_depth | tfxParticleManagerFlags_guarantee_order;

	if (double_buffered_sprites)
		pm->flags |= tfxParticleManagerFlags_double_buffer_sprites;

	pm->flags |= dynamic_sprite_allocation ? tfxParticleManagerFlags_dynamic_sprite_allocation : 0;
	pm->flags |= group_sprites_by_effect ? tfxParticleManagerFlags_use_effect_sprite_buffers : 0;

	for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
		pm->effects_in_use[depth][0].reserve(pm->max_effects);
		pm->effects_in_use[depth][1].reserve(pm->max_effects);
	}
	pm->control_emitter_queue.reserve(pm->max_effects);

	pm->free_effects.reserve(pm->max_effects);
	//InitEffectSoA(&pm->effect_buffers, &pm->effects, pm->max_effects);
	pm->emitters.set_alignment(16);
	pm->effects.set_alignment(16);
	pm->emitters.reserve(pm->max_effects);
	pm->effects.reserve(pm->max_effects);
	pm->particle_indexes.reserve(effects_limit);    //todo: Handle this better.
	pm->spawn_work.reserve(effects_limit);
	pm->control_work.reserve(effects_limit);
	pm->age_work.reserve(effects_limit);
	//for (tfxEachLayer) {
		//pm->current_depth_buffer_index[layer] = 0;
	//}
}

void InitializeParticleManager(tfx_particle_manager_t *pm, tfx_library_t *library, tfx_particle_manager_info_t info) {
	pm->info = info;
	InitCommonParticleManager(pm, library, info.max_particles, info.max_effects, info.order_mode, info.double_buffer_sprites, info.dynamic_sprite_allocation, info.group_sprites_by_effect, info.multi_threaded_batch_size);

	pm->flags |= info.is_3d ? tfxParticleManagerFlags_3d_effects : 0;
	pm->flags |= info.auto_order_effects ? tfxParticleManagerFlags_auto_order_effects : 0;

	for (tfxEachLayer) {
		pm->max_cpu_particles_per_layer[layer] = info.max_particles;
	}

	if (pm->flags & tfxParticleManagerFlags_3d_effects) {
		pm->instance_buffer = tfxCreateBuffer(sizeof(tfx_3d_instance_t), 16);
	}
	else {
		pm->instance_buffer = tfxCreateBuffer(sizeof(tfx_2d_instance_t), 16);
	}
	if (info.write_direct_to_staging_buffer) {
		pm->flags |= tfxParticleManagerFlags_direct_to_staging_buffer;
		pm->instance_buffer.data = nullptr;
	} else {
		pm->instance_buffer.reserve(tfxMax((info.max_particles / tfxDataWidth + 1) * tfxDataWidth, 8));
	}

	if (!(pm->flags & tfxParticleManagerFlags_use_effect_sprite_buffers) && (pm->flags & tfxParticleManagerFlags_ordered_by_age || pm->flags & tfxParticleManagerFlags_order_by_depth)) {
		FreeAllParticleLists(pm);
	}
}

void SetEffectPosition(tfx_particle_manager_t *pm, tfxEffectID effect_index, float x, float y) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	tfx_vec2_t position(x, y);
	pm->effects[effect_index].local_position = position;
}

void SetEffectPosition(tfx_particle_manager_t *pm, tfxEffectID effect_index, tfx_vec2_t position) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].local_position = position;
}

void SetEffectPosition(tfx_particle_manager_t *pm, tfxEffectID effect_index, float x, float y, float z) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	tfx_vec3_t position(x, y, z);
	pm->effects[effect_index].local_position = position;
}

void SetEffectPosition(tfx_particle_manager_t *pm, tfxEffectID effect_index, tfx_vec3_t position) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].local_position = position;
}

void SetAnimationPosition(tfx_animation_manager_t *animation_manager, tfxAnimationID effect_index, float position[3]) {
	animation_manager->instances[effect_index].position.x = position[0];
	animation_manager->instances[effect_index].position.y = position[1];
	animation_manager->instances[effect_index].position.z = position[2];
}

tfx_animation_instance_t *GetAnimationInstance(tfx_animation_manager_t *animation_manager, tfxAnimationID animation_id) {
	return &animation_manager->instances[animation_id];
}

void SetAnimationPosition(tfx_animation_manager_t *animation_manager, tfxAnimationID effect_index, float x, float y) {
	animation_manager->instances[effect_index].position.x = x;
	animation_manager->instances[effect_index].position.y = y;
}

void SetAnimationScale(tfx_animation_manager_t *animation_manager, tfxAnimationID effect_index, float scale) {
	animation_manager->instances[effect_index].scale = scale;
}

void MoveEffect(tfx_particle_manager_t *pm, tfxEffectID effect_index, tfx_vec3_t amount) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].local_position += amount;
}

void MoveEffect(tfx_particle_manager_t *pm, tfxEffectID effect_index, float x, float y, float z) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].local_position += {x, y, z};
}

tfxAPI tfx_vec3_t GetEffectPosition(tfx_particle_manager_t *pm, tfxEffectID effect_index) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	return pm->effects[effect_index].local_position;
}

tfxAPI tfx_2d_instance_t *GetEffectSpriteBuffer(tfx_particle_manager_t *pm, tfxEffectID effect_index, tfxU32 *sprite_count) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	tfx_effect_instance_data_t &instance_data = pm->effects[effect_index].instance_data;
	*sprite_count = instance_data.instance_count;
	return &tfxCastBufferRef(tfx_2d_instance_t, pm->instance_buffer)[instance_data.instance_start_index];
}

tfxAPI tfx_3d_instance_t *GetEffectBillboardBuffer(tfx_particle_manager_t *pm, tfxEffectID effect_index, tfxU32 *sprite_count) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	tfx_effect_instance_data_t &instance_data = pm->effects[effect_index].instance_data;
	*sprite_count = instance_data.instance_count;
	return &tfxCastBufferRef(tfx_3d_instance_t, pm->instance_buffer)[instance_data.instance_start_index];
}

bool GetNextSpriteBuffer(tfx_particle_manager_t *pm, tfx_2d_instance_t **instances, tfx_effect_instance_data_t **instance_data, tfxU32 *instance_count) {
	if (pm->effect_index_position >= pm->effects_in_use[0][pm->current_ebuff].current_size) {
		*instances = nullptr;
		*instance_data = nullptr;
		*instance_count = 0;
		return false;
	}
	tfx_effect_index_t effect_index = pm->effects_in_use[0][pm->current_ebuff][pm->effect_index_position];
	pm->effect_index_position++;
	tfx_effect_instance_data_t &data = pm->effects[effect_index.index].instance_data;
	*instance_count = data.instance_count;
	*instances = &tfxCastBufferRef(tfx_2d_instance_t, pm->instance_buffer)[data.instance_start_index];
	*instance_data = &data;
	return true;
}

bool GetNextBillboardBuffer(tfx_particle_manager_t *pm, tfx_3d_instance_t **instances, tfx_effect_instance_data_t **instance_data, tfxU32 *instance_count) {
	if (pm->effect_index_position >= pm->effects_in_use[0][pm->current_ebuff].current_size) {
		*instances = nullptr;
		*instance_data = nullptr;
		*instance_count = 0;
		return false;
	}
	tfx_effect_index_t effect_index = pm->effects_in_use[0][pm->current_ebuff][pm->effect_index_position];
	pm->effect_index_position++;
	tfx_effect_instance_data_t &data = pm->effects[effect_index.index].instance_data;
	*instance_count = data.instance_count;
	*instances = &tfxCastBufferRef(tfx_3d_instance_t, pm->instance_buffer)[data.instance_start_index];
	*instance_data = &data;
	return true;
}

void ResetInstanceBufferLoopIndex(tfx_particle_manager_t *pm) {
	pm->effect_index_position = 0;
}

void SetEffectRotation(tfx_particle_manager_t *pm, tfxEffectID effect_index, float rotation) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].local_rotations.roll = rotation;
	pm->effects[effect_index].state_flags |= tfxEffectStateFlags_override_orientiation;
}

void SetEffectRoll(tfx_particle_manager_t *pm, tfxEffectID effect_index, float roll) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].local_rotations.roll = roll;
	pm->effects[effect_index].state_flags |= tfxEffectStateFlags_override_orientiation;
}

void SetEffectPitch(tfx_particle_manager_t *pm, tfxEffectID effect_index, float pitch) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].local_rotations.pitch = pitch;
	pm->effects[effect_index].state_flags |= tfxEffectStateFlags_override_orientiation;
}

void SetEffectYaw(tfx_particle_manager_t *pm, tfxEffectID effect_index, float pitch) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].local_rotations.pitch = pitch;
	pm->effects[effect_index].state_flags |= tfxEffectStateFlags_override_orientiation;
}

void SetEffectWidthMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float width) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].emitter_size.x = width;
	pm->effects[effect_index].state_flags |= tfxEffectStateFlags_override_size_multiplier;
}

void SetEffectHeightMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float height) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].emitter_size.y = height;
	pm->effects[effect_index].state_flags |= tfxEffectStateFlags_override_size_multiplier;
}

void SetEffectDepthMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float depth) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].emitter_size.z = depth;
	pm->effects[effect_index].state_flags |= tfxEffectStateFlags_override_size_multiplier;
}

void SetEffectLifeMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float life) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].spawn_controls.life = life;
}

void SetEffectParticleWidthMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float width) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].spawn_controls.size_x = width;
}

void SetEffectParticleHeightMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float height) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].spawn_controls.size_y = height;
}

void SetEffectVelocityMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float velocity) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].spawn_controls.velocity = velocity;
}

void SetEffectSpinMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float spin) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].spawn_controls.spin = spin;
}

void SetEffectIntensityMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float intensity) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].spawn_controls.intensity = intensity;
}

void SetEffectSplatterMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float splatter) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].spawn_controls.splatter = splatter;
}

void SetEffectWeightMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float weight) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].spawn_controls.weight = weight;
}

void SetEffectOveralScale(tfx_particle_manager_t *pm, tfxEffectID effect_index, float overal_scale) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].overal_scale = overal_scale;
	pm->effects[effect_index].state_flags |= tfxEffectStateFlags_override_overal_scale;
}

void SetEffectBaseNoiseOffset(tfx_particle_manager_t *pm, tfxEffectID effect_index, float noise_offset) {
	TFX_ASSERT(ValidEffectID(pm, effect_index));    //Not a valid effect id. Make sure that when you call AddEffectToParticleManager you check that it returns true.
	pm->effects[effect_index].noise_base_offset = noise_offset;
}

}        //Namespace
