#include "timelinefx.h"

namespace tfx {

	//A 2d Simd (SSE3) version of simplex noise allowing you to do 4 samples with 1 call for a speed boost
	tfx128Array tfxNoise4(const tfx128 x4, const tfx128 y4) {
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

		for (int i = 0; i < 4; ++i)
		{
			gi0[i] = permMOD12[perm[ii.a[i] + perm[jj.a[i]]]];
			gi1[i] = permMOD12[perm[ii.a[i] + i1.a[i] + perm[jj.a[i] + j1.a[i]]]];
			gi2[i] = permMOD12[perm[ii.a[i] + 1 + perm[jj.a[i] + 1]]];
		}

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
		n0 = _mm_and_ps(_mm_mul_ps(_mm_mul_ps(t02, t02), Dot128XY(gx0, gy0, x0, y0)), _mm_cmpge_ps(t0, _mm_setzero_ps()));

		tfx128 t1 = _mm_sub_ps(_mm_sub_ps(_mm_set1_ps(0.5f), _mm_mul_ps(x1, x1)), _mm_mul_ps(y1, y1));
		tfx128 t12 = _mm_mul_ps(t1, t1);
		n1 = _mm_and_ps(_mm_mul_ps(_mm_mul_ps(t12, t12), Dot128XY(gx1, gy1, x1, y1)), _mm_cmpge_ps(t1, _mm_setzero_ps()));

		tfx128 t2 = _mm_sub_ps(_mm_sub_ps(_mm_set1_ps(0.5f), _mm_mul_ps(x2, x2)), _mm_mul_ps(y2, y2));
		tfx128 t22 = _mm_mul_ps(t2, t2);
		n2 = _mm_and_ps(_mm_mul_ps(_mm_mul_ps(t22, t22), Dot128XY(gx2, gy2, x2, y2)), _mm_cmpge_ps(t2, _mm_setzero_ps()));

		tfx128Array result;
		result.m = _mm_mul_ps(_mm_set1_ps(45.23065f), _mm_add_ps(n0, _mm_add_ps(n1, n2)));
		return result;
	}

	//A 3d Simd (SSE3) version of simplex noise allowing you to do 4 samples with 1 call for a speed boost
	tfx128Array tfxNoise4(const tfx128 &x4, const tfx128 &y4, const tfx128 &z4) {
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

		for (int i = 0; i < 4; ++i)
		{
			gi0.a[i] = permMOD12[ii.a[i] + perm[jj.a[i] + perm[kk.a[i]]]];
			gi1.a[i] = permMOD12[ii.a[i] + i1.a[i] + perm[jj.a[i] + j1.a[i] + perm[kk.a[i] + k1.a[i]]]];
			gi2.a[i] = permMOD12[ii.a[i] + i2.a[i] + perm[jj.a[i] + j2.a[i] + perm[kk.a[i] + k2.a[i]]]];
			gi3.a[i] = permMOD12[ii.a[i] + 1 + perm[jj.a[i] + 1 + perm[kk.a[i] + 1]]];
		}

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

		tfx128 n0 = _mm_mul_ps(t0q, Dot128XYZ(gi0x.m, gi0y.m, gi0z.m, x0, y0, z0));
		tfx128 n1 = _mm_mul_ps(t1q, Dot128XYZ(gi1x.m, gi1y.m, gi1z.m, x1, y1, z1));
		tfx128 n2 = _mm_mul_ps(t2q, Dot128XYZ(gi2x.m, gi2y.m, gi2z.m, x2, y2, z2));
		tfx128 n3 = _mm_mul_ps(t3q, Dot128XYZ(gi3x.m, gi3y.m, gi3z.m, x3, y3, z3));

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

		tfxvec<tfxFace> triangles;
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
	tfxWideFloat tfxUPDATE_TIME_WIDE = tfxWideSetSingle(tfxUPDATE_TIME);
	float tfxFRAME_LENGTH = 1000.f / tfxUPDATE_FREQUENCY;
	tfxWideFloat tfxFRAME_LENGTH_WIDE = tfxWideSetSingle(1000.f / tfxUPDATE_FREQUENCY);

	//Set the udpate frequency for all particle effects - There may be options in the future for individual effects to be updated at their own specific frequency.
	void SetUpdateFrequency(float fps) {
		tfxUPDATE_FREQUENCY = fps;
		tfxUPDATE_TIME = 1.f / tfxUPDATE_FREQUENCY;
		tfxUPDATE_TIME_WIDE = tfxWideSetSingle(tfxUPDATE_TIME);
		tfxFRAME_LENGTH_WIDE = tfxWideSetSingle(1000.f / tfxUPDATE_FREQUENCY);
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
		if (fwrite(data, 1, Length(), file) != Length()) {
			fclose(file);
			return false;
		}

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

	tfxStr512 tfxStream::ReadLine() {
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
		fclose(file);
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
			//FILE *file = fopen(file_path.c_str(), "rb");
			FILE *file;
			errno_t error = fopen_s(&file, file_path.c_str(), "rb");
			if (error != 0) {
				printf("strerror says open failed: %s\n", strerror(error));
			}
			assert(file);		//couldn't open the file!
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

	void tfxPackage::AddFile(const char *file_name, tfxStream &data) {
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
	tfxStream ReadEntireFile(const char *file_name, bool terminate) {
		tfxStream buffer;
		FILE *file = fopen(file_name, "rb");
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

	tfxStream SavePackageMemory(tfxPackage &package) {
		if (!package.file_path.Length()) return false;											//Package must have a file path
		if (package.header.magic_number != tfxMAGIC_NUMBER) return false;						//Header of package must contain correct magic number. CreatePackage to correctly initialise a package.
		if (package.inventory.magic_number != tfxMAGIC_NUMBER_INVENTORY) return false;			//Inventory of package must contain correct magic number

		//char *file = (char*)malloc(GetPackageSize(package));
		tfxStream file(GetPackageSize(package));
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

	tfxErrorFlags LoadPackage(tfxStream &stream, tfxPackage &package) {
		//Note: tfxStream does not copy the memory, only the pointer, so if you FreeAll on the stream you pass in it will also free the file_data here as well
		package.file_data = stream;
		package.file_data.Seek(0);
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

	bool tfxEffectEmitter::IsFinite() {
		for (auto &e : GetInfo().sub_effectors) {
			float qty = e.library->emitter_attributes[e.emitter_attributes].base.amount.GetLastValue() + e.library->emitter_attributes[e.emitter_attributes].variation.amount.GetLastValue();
			if (!(e.property_flags & tfxEmitterPropertyFlags_single) && qty > 0)
				return false;
			else if (e.property_flags & tfxEmitterPropertyFlags_single && e.GetProperties().single_shot_limit[e.property_index] == 0)
				return false;

		}
		return true;
	}

	void tfxEffectEmitter::FlagAs3D(bool flag) {
		if (flag)
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

	tfxU32 tfxEffectEmitter::GetCameraSettingsIndex() {
		return GetInfo().preview_camera_settings;
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
		if (e.type == tfxEffectType)
			e.GetInfo().name = "New Effect";
		else
			e.GetInfo().name = "New Emitter";
		GetInfo().sub_effectors.push_back(e);
		library->UpdateEffectPaths();
		ReIndex();
		return GetInfo().sub_effectors.back();
	}

	float GetEmissionDirection2d(tfxParticleManager &pm, tfxLibrary *library, tfxRandom &random, tfxU32 property_index, tfxU32 emitter_index, tfxVec2 local_position, tfxVec2 world_position, tfxVec2 emitter_size) {
		//float (*effect_lookup_callback)(tfxGraph &graph, float age) = common.root_effect->lookup_mode == tfxPrecise ? LookupPrecise : LookupFast;
		const float frame = pm.emitters.frame[emitter_index];
		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[emitter_index];
		float emission_angle = lookup_callback(library->emitter_attributes[emitter_attributes].properties.emission_pitch, frame);
		float emission_angle_variation = lookup_callback(library->emitter_attributes[emitter_attributes].properties.emission_range, frame);
		//----Emission
		float range = emission_angle_variation * .5f;
		float direction = 0;
		tfxEmissionType emission_type = library->emitter_properties.emission_type[property_index];
		tfxEmissionDirection emission_direction = library->emitter_properties.emission_direction[property_index];

		if (emission_type == tfxPoint)
			return direction + emission_angle + random.Range(-range, range);

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
		return direction + emission_angle + random.Range(-range, range);
	}

	tfxVec3 GetEmissionDirection3d(tfxParticleManager &pm, tfxLibrary *library, tfxRandom &random, tfxU32 property_index, tfxU32 emitter_index, float emission_pitch, float emission_yaw, tfxVec3 local_position, tfxVec3 world_position, tfxVec3 emitter_size) {
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
			result.y = random.Generate() * (1.f - cos(range)) + cos(range);
			float phi = random.Generate() * 2.f * tfxPI;
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
		if (compile)
			library->CompileGlobalGraph(global);
	}

	void tfxEffectEmitter::ResetTransformGraphs(bool add_node, bool compile) {
		library->transform_attributes[transform_attributes].roll.Reset(0.f, tfxAnglePreset, add_node); library->transform_attributes[transform_attributes].roll.type = tfxTransform_roll;
		library->transform_attributes[transform_attributes].pitch.Reset(0.f, tfxAnglePreset, add_node); library->transform_attributes[transform_attributes].pitch.type = tfxTransform_pitch;
		library->transform_attributes[transform_attributes].yaw.Reset(0.f, tfxAnglePreset, add_node); library->transform_attributes[transform_attributes].yaw.type = tfxTransform_yaw;
		library->transform_attributes[transform_attributes].translation_x.Reset(0.f, tfxTranslationPreset, add_node); library->transform_attributes[transform_attributes].translation_x.type = tfxTransform_translate_x;
		library->transform_attributes[transform_attributes].translation_y.Reset(0.f, tfxTranslationPreset, add_node); library->transform_attributes[transform_attributes].translation_y.type = tfxTransform_translate_y;
		library->transform_attributes[transform_attributes].translation_z.Reset(0.f, tfxTranslationPreset, add_node); library->transform_attributes[transform_attributes].translation_z.type = tfxTransform_translate_z;
		if (compile)
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
		if (compile)
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
		if (compile)
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
		if (compile)
			library->CompileVariationGraph(emitter_attributes);
	}

	void tfxEffectEmitter::ResetOvertimeGraphs(bool add_node, bool compile) {
		library->emitter_attributes[emitter_attributes].overtime.velocity.Reset(1.f, tfxVelocityOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.velocity.type = tfxOvertime_velocity;
		library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.Reset(1.f, tfxGlobalPercentPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.type = tfxOvertime_velocity_adjuster;
		library->emitter_attributes[emitter_attributes].overtime.width.Reset(1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.width.type = tfxOvertime_width;
		library->emitter_attributes[emitter_attributes].overtime.height.Reset(1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.height.type = tfxOvertime_height;
		library->emitter_attributes[emitter_attributes].overtime.weight.Reset(1.f, tfxWeightOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.weight.type = tfxOvertime_weight;
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
		if (compile)
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

		if (type == tfxEffectType && global != tfxINVALID) {
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

		if (type == tfxEmitterType && emitter_attributes != tfxINVALID) {
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
			if (library->emitter_attributes[emitter_attributes].overtime.weight.nodes.size() == 0) library->emitter_attributes[emitter_attributes].overtime.weight.Reset(1.f, tfxWeightOvertimePreset);
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
		tfxLibrary *library = emitter->library;
		tmpStack(tfxEffectEmitter, stack);
		stack.push_back(*emitter);
		while (stack.size()) {
			tfxEffectEmitter &current = stack.pop_back();
			if (current.type == tfxEffectType && !current.parent) {
				library->FreeGlobal(current.global);
			}
			else if (current.type == tfxEmitterType) {
				library->FreeEmitterAttributes(current.emitter_attributes);
			}
			for (auto &sub : current.GetInfo().sub_effectors) {
				stack.push_back(sub);
			}
		}
		GetInfo().sub_effectors.erase(emitter);

		ReIndex();
		if (library)
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
				else if (current.type == tfxEmitterType) {
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

	void tfxEffectEmitter::Clone(tfxEffectEmitter &clone, tfxEffectEmitter *root_parent, tfxLibrary *destination_library, tfxEffectCloningFlags flags) {
		//tfxU32 size = library->global_graphs[0].amount.lookup.values.capacity;
		clone = *this;
		clone.info_index = clone.library->CloneInfo(info_index, destination_library);
		if (clone.type != tfxFolder) {
			clone.property_index = clone.library->CloneProperties(property_index, destination_library);
		}
		clone.property_flags |= tfxEmitterPropertyFlags_enabled;
		if (!(flags & tfxEffectCloningFlags_keep_user_data))
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
					if (flags & tfxEffectCloningFlags_compile_graphs)
						clone.library->CompileGlobalGraph(clone.global);
				}
			}
		}
		else if (type == tfxEmitterType) {
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
				if (!(flags & tfxEffectCloningFlags_keep_user_data))
					emitter_copy.user_data = nullptr;
				clone.AddEmitter(emitter_copy);
			}
			else if (e.type == tfxEffectType) {
				tfxEffectEmitter effect_copy;
				if (clone.type == tfxFolder)
					e.Clone(effect_copy, &effect_copy, destination_library, flags);
				else
					e.Clone(effect_copy, root_parent, destination_library, flags);
				if (!(flags & tfxEffectCloningFlags_keep_user_data))
					effect_copy.user_data = nullptr;
				clone.AddEffect(effect_copy);
			}
		}
	}

	bool PrepareEffectTemplate(tfxLibrary &library, const char *name, tfxEffectTemplate &effect_template) {
		effect_template.Reset();
		if (library.effect_paths.ValidName(name)) {
			library.PrepareEffectTemplate(name, effect_template);
			return true;
		}
		return false;
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
			library->transform_attributes[transform_attributes].roll.Free();
			library->transform_attributes[transform_attributes].pitch.Free();
			library->transform_attributes[transform_attributes].yaw.Free();
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
		else if (type == tfxEmitterType) {
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

	tfxEffectEmitter& tfxLibrary::operator[] (uint32_t index) {
		return effects[index];
	}

	bool tfxLibrary::RenameEffect(tfxEffectEmitter &effect, const char *new_name) {
		if (!NameExists(effect, new_name) && strlen(new_name) > 0) {
			effect.SetName(new_name);
			UpdateEffectPaths();
			return true;
		}

		return false;
	}

	bool tfxLibrary::NameExists(tfxEffectEmitter& effect, const char *name) {
		for (auto &e : effects) {
			if (effect.library_index != e.library_index) {
				if (e.GetInfo().name == name) {
					return true;
				}
			}
		}

		return false;
	}

	bool tfxLibrary::NameExists2(tfxEffectEmitter& effect, const char *name) {
		for (auto &e : effects) {
			if (e.GetInfo().name == name) {
				return true;
			}
		}
		return false;
	}

	void tfxLibrary::UpdateEffectPaths() {
		effect_paths.Clear();
		for (auto &e : effects) {
			tfxStr256 path = e.GetInfo().name;
			e.path_hash = tfxXXHash64::hash(path.c_str(), path.Length(), 0);
			AddPath(e, path);
		}
	}

	void tfxLibrary::AddPath(tfxEffectEmitter &effect_emitter, tfxStr256 &path) {
		effect_paths.Insert(path, &effect_emitter);
		for (auto &sub : effect_emitter.GetInfo().sub_effectors) {
			tfxStr256 sub_path = path;
			sub_path.Appendf("/%s", sub.GetInfo().name.c_str());
			sub.path_hash = tfxXXHash64::hash(sub_path.c_str(), sub_path.Length(), 0);
			AddPath(sub, sub_path);
		}
	}

	tfxEffectEmitter &tfxLibrary::InsertEffect(tfxEffectEmitter &effect, tfxEffectEmitter *position) {
		effect.library_index = effects.current_size;
		effect.type = tfxEffectType;
		effect.GetInfo().uid = ++uid;
		effect.library = this;
		tfxEffectEmitter *inserted_effect = effects.insert_after(position, effect);
		ReIndex();
		UpdateEffectPaths();
		return *inserted_effect;
	}

	tfxEffectEmitter &tfxLibrary::AddEffect(tfxEffectEmitter &effect) {
		effect.library_index = effects.current_size;
		effect.type = tfxEffectType;
		effect.GetInfo().uid = ++uid;
		effect.library = this;
		effects.push_back(effect);
		ReIndex();
		UpdateEffectPaths();
		return effects.back();
	}

	tfxEffectEmitter &tfxLibrary::AddFolder(tfxStr64 &name) {
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

	tfxEffectEmitter &tfxLibrary::AddFolder(tfxEffectEmitter &folder) {
		assert(folder.type == tfxFolder);			//Must be type tfxFolder if adding a folder
		folder.library = this;
		folder.GetInfo().uid = ++uid;
		effects.push_back(folder);
		ReIndex();
		UpdateEffectPaths();
		return effects.back();
	}

	tfxEffectEmitter &tfxLibrary::AddStage(tfxStr64 &name) {
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

	tfxEffectEmitter* tfxLibrary::GetEffect(tfxStr256 &path) {
		assert(effect_paths.ValidName(path));		//Effect was not found by that name
		return effect_paths.At(path);
	}

	tfxEffectEmitter* tfxLibrary::GetEffect(const char *path) {
		assert(effect_paths.ValidName(path));		//Effect was not found by that name
		return effect_paths.At(path);
	}

	tfxEffectEmitter* tfxLibrary::GetEffect(tfxKey key) {
		assert(effect_paths.ValidKey(key));			//Effect was not found by that key
		return effect_paths.At(key);
	}

	void tfxLibrary::PrepareEffectTemplate(tfxStr256 path, tfxEffectTemplate &effect_template) {
		tfxEffectEmitter *effect = GetEffect(path);
		assert(effect);								//Effect was not found, make sure the path exists
		assert(effect->type == tfxEffectType);		//The effect must be an effect type, not an emitter
		effect_template.original_effect_hash = effect->path_hash;
		effect->Clone(effect_template.effect, &effect_template.effect, this, tfxEffectCloningFlags_clone_graphs | tfxEffectCloningFlags_compile_graphs);
		effect_template.AddPath(effect_template.effect, effect_template.effect.GetInfo().name.c_str());
	}

	void tfxLibrary::PrepareEffectTemplate(tfxEffectEmitter &effect, tfxEffectTemplate &effect_template) {
		assert(effect.type == tfxEffectType);
		effect.Clone(effect_template.effect, &effect_template.effect, this);
		effect_template.AddPath(effect_template.effect, effect_template.effect.GetInfo().name.c_str());
	}

	void tfxLibrary::ReIndex() {
		tfxU32 index = 0;
		for (auto &e : effects) {
			e.library_index = index++;
			e.parent = nullptr;
			e.ReIndex();
		}
	}

	void tfxLibrary::UpdateParticleShapeReferences(tfxvec<tfxEffectEmitter> &effects, tfxU32 default_index) {
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

	tfxEffectEmitter* tfxLibrary::MoveUp(tfxEffectEmitter &effect) {
		if (effect.library_index > 0) {
			tfxU32 new_index = effect.library_index - 1;
			std::swap(effects[effect.library_index], effects[new_index]);
			UpdateEffectPaths();
			ReIndex();
			return &effects[new_index];
		}
		return nullptr;
	}

	tfxEffectEmitter* tfxLibrary::MoveDown(tfxEffectEmitter &effect) {
		if (effect.library_index < effects.size() - 1) {
			tfxU32 new_index = effect.library_index + 1;
			std::swap(effects[effect.library_index], effects[new_index]);
			UpdateEffectPaths();
			ReIndex();
			return &effects[new_index];
		}
		return nullptr;
	}

	void tfxLibrary::BuildComputeShapeData(void* dst, tfxVec4(uv_lookup)(void *ptr, tfxComputeImageData &image_data, int offset)) {
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

	void tfxLibrary::CopyComputeShapeData(void* dst) {
		assert(shape_data.size());	//You must call BuildComputeShapeData first
		memcpy(dst, shape_data.data, shape_data.size() * sizeof(tfxComputeImageData));
	}

	void tfxLibrary::CopyLookupIndexesData(void* dst) {
		assert(dst);	//must be a valid pointer to a space in memory
		assert(compiled_lookup_indexes.size());		//There is no data to copy, make sure a library has been loaded properly and it contains effects with emitters
		tfxGraphLookupIndex *test = static_cast<tfxGraphLookupIndex*>(dst);
		memcpy(dst, compiled_lookup_indexes.data, GetLookupIndexesSizeInBytes());
	}

	void tfxLibrary::CopyLookupValuesData(void* dst) {
		assert(dst);	//must be a valid pointer to a space in memory
		assert(compiled_lookup_indexes.size());		//There is no data to copy, make sure a library has been loaded properly and it contains effects with emitters
		memcpy(dst, compiled_lookup_values.data, GetLookupValuesSizeInBytes());
	}

	tfxU32 tfxLibrary::GetComputeShapeDataSizeInBytes() {
		tfxU32 frame_count = 0;
		for (auto &shape : particle_shapes.data) {
			frame_count += (tfxU32)shape.animation_frames;
		}
		return frame_count * sizeof(tfxComputeImageData);
	}

	tfxU32 tfxLibrary::GetComputeShapeCount() {
		tfxU32 frame_count = 0;
		for (auto &shape : particle_shapes.data) {
			frame_count += (tfxU32)shape.animation_frames;
		}
		return frame_count;
	}

	tfxU32 tfxLibrary::GetLookupIndexCount() {
		return compiled_lookup_indexes.size() * tfxOvertimeCount;
	}

	tfxU32 tfxLibrary::GetLookupValueCount() {
		return compiled_lookup_values.size();
	}

	tfxU32 tfxLibrary::GetLookupIndexesSizeInBytes() {
		return sizeof(tfxGraphLookupIndex) * tfxOvertimeCount * compiled_lookup_indexes.size();
	}

	tfxU32 tfxLibrary::GetLookupValuesSizeInBytes() {
		return sizeof(float) * compiled_lookup_values.size();
	}

	void tfxLibrary::RemoveShape(tfxU32 shape_index) {
		particle_shapes.RemoveInt(shape_index);
		for (auto &m : particle_shapes.map) {
			particle_shapes[m.index].shape_index = (tfxU32)m.key;
		}
	}

	void tfxLibrary::DeleteEffect(tfxEffectEmitter *effect) {
		effects[effect->library_index].CleanUp();
		effects.erase(&effects[effect->library_index]);

		UpdateEffectPaths();
		ReIndex();
	}

	tfxU32 tfxLibrary::AddGlobal() {
		if (free_global_graphs.size())
			return free_global_graphs.pop_back();
		tfxGlobalAttributes global;
		global.Initialise(&graph_node_allocator, &graph_lookup_allocator);
		global_graphs.push_back(global);
		return global_graphs.size() - 1;
	}

	tfxU32 tfxLibrary::AddKeyframes() {
		if (free_keyframes.size())
			return free_keyframes.pop_back();
		tfxTransformAttributes keyframes;
		keyframes.Initialise(&graph_node_allocator, &graph_lookup_allocator);
		transform_attributes.push_back(keyframes);
		return transform_attributes.size() - 1;
	}

	tfxU32 tfxLibrary::AddEmitterAttributes() {
		if (free_emitter_attributes.size())
			return free_emitter_attributes.pop_back();
		tfxEmitterAttributes attributes;
		attributes.Initialise(&graph_node_allocator, &graph_lookup_allocator);
		emitter_attributes.push_back(attributes);
		return emitter_attributes.size() - 1;
	}

	void tfxLibrary::FreeGlobal(tfxU32 index) {
		assert(index < global_graphs.size());
		free_global_graphs.push_back(index);
		global_graphs[index].Free();
	}

	void tfxLibrary::FreeKeyframes(tfxU32 index) {
		assert(index < transform_attributes.size());
		free_keyframes.push_back(index);
		transform_attributes[index].Free();
	}

	void tfxLibrary::FreeEmitterAttributes(tfxU32 index) {
		assert(index < emitter_attributes.size());
		free_emitter_attributes.push_back(index);
		emitter_attributes[index].Free();
	}

	void tfxLibrary::FreeProperties(tfxU32 index) {
		assert(index < emitter_properties_buffer.current_size);
		free_properties.push_back(index);
	}

	void tfxLibrary::FreeInfo(tfxU32 index) {
		assert(index < effect_infos.size());
		free_infos.push_back(index);
	}

	tfxU32 tfxLibrary::CountKeyframeLookUpValues(tfxU32 index) {
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

	tfxU32 tfxLibrary::CountGlobalLookUpValues(tfxU32 index) {
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

	tfxU32 tfxLibrary::CountEmitterLookUpValues(tfxU32 index) {
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

	tfxU32 tfxLibrary::CloneGlobal(tfxU32 source_index, tfxLibrary *destination_library) {
		tfxU32 index = destination_library->AddGlobal();
		global_graphs[source_index].CopyToNoLookups(&destination_library->global_graphs[index]);
		return index;
	}

	tfxU32 tfxLibrary::CloneKeyframes(tfxU32 source_index, tfxLibrary *destination_library) {
		tfxU32 index = destination_library->AddKeyframes();
		transform_attributes[source_index].CopyToNoLookups(&destination_library->transform_attributes[index]);
		return index;
	}

	tfxU32 tfxLibrary::CloneEmitterAttributes(tfxU32 source_index, tfxLibrary *destination_library) {
		tfxU32 index = destination_library->AddEmitterAttributes();
		emitter_attributes[source_index].properties.CopyToNoLookups(&destination_library->emitter_attributes[index].properties);
		emitter_attributes[source_index].base.CopyToNoLookups(&destination_library->emitter_attributes[index].base);
		emitter_attributes[source_index].variation.CopyToNoLookups(&destination_library->emitter_attributes[index].variation);
		emitter_attributes[source_index].overtime.CopyToNoLookups(&destination_library->emitter_attributes[index].overtime);
		return index;
	}

	tfxU32 tfxLibrary::CloneInfo(tfxU32 source_index, tfxLibrary *destination_library) {
		tfxU32 index = destination_library->AddEffectEmitterInfo();
		destination_library->effect_infos[index] = effect_infos[source_index];
		return index;
	}

	tfxU32 tfxLibrary::CloneProperties(tfxU32 source_index, tfxLibrary *destination_library) {
		tfxU32 index = destination_library->AddEmitterProperties();
		CopyEmitterProperites(emitter_properties, source_index, destination_library->emitter_properties, index);
		return index;
	}

	void tfxLibrary::AddEmitterGraphs(tfxEffectEmitter& emitter) {
		emitter.emitter_attributes = AddEmitterAttributes();
	}

	void tfxLibrary::AddTransformGraphs(tfxEffectEmitter& emitter) {
		emitter.transform_attributes = AddKeyframes();
	}

	void tfxLibrary::AddEffectGraphs(tfxEffectEmitter& effect) {
		tfxEffectEmitter *root_effect = effect.GetRootEffect();
		if (root_effect == &effect)
			effect.global = AddGlobal();
		else
			effect.global = root_effect->global;
	}

	tfxU32 tfxLibrary::AddSpriteSheetSettings(tfxEffectEmitter& effect) {
		assert(effect.type == tfxEffectType);
		tfxSpriteSheetSettings a;
		a.frames = 32;
		a.current_frame = 1;
		a.frame_offset = 0;
		a.extra_frames_count = 0;
		a.position = tfxVec2(0.f, 0.f);
		a.frame_size = tfxVec2(256.f, 256.f);
		a.playback_speed = 1.f;
		a.animation_flags = tfxAnimationFlags_needs_recording | tfxAnimationFlags_export_with_transparency;
		a.seed = 0;
		a.zoom = 1.f;
		a.scale = 1.f;
		a.needs_exporting = 0;
		a.color_option = tfxExportColorOptions::tfxFullColor;
		a.export_option = tfxExportOptions::tfxSpriteSheet;
		a.camera_settings.camera_floor_height = -10.f;
		a.camera_settings.camera_fov = tfxRadians(60);
		a.camera_settings.camera_pitch = tfxRadians(-30.f);
		a.camera_settings.camera_yaw = tfxRadians(-90.f);
		a.camera_settings.camera_position = tfxVec3(0.f, 3.5f, 7.5f);
		a.camera_settings.camera_isometric = false;
		a.camera_settings.camera_isometric_scale = 5.f;
		a.camera_settings.camera_hide_floor = false;
		a.camera_settings_orthographic.camera_floor_height = -10.f;
		a.camera_settings_orthographic.camera_fov = tfxRadians(60);
		a.camera_settings_orthographic.camera_pitch = tfxRadians(-30.f);
		a.camera_settings_orthographic.camera_yaw = tfxRadians(-90.f);
		a.camera_settings_orthographic.camera_position = tfxVec3(0.f, 3.5f, 7.5f);
		a.camera_settings_orthographic.camera_isometric = true;
		a.camera_settings_orthographic.camera_isometric_scale = 5.f;
		a.camera_settings_orthographic.camera_hide_floor = false;
		sprite_sheet_settings.push_back(a);
		effect.GetInfo().sprite_sheet_settings_index = sprite_sheet_settings.size() - 1;
		return effect.GetInfo().sprite_sheet_settings_index;
	}

	tfxU32 tfxLibrary::AddSpriteDataSettings(tfxEffectEmitter& effect) {
		assert(effect.type == tfxEffectType);
		tfxSpriteDataSettings a;
		a.frames = 32;
		a.current_frame = 1;
		a.frame_offset = 0;
		a.extra_frames_count = 0;
		a.playback_speed = 1.f;
		a.animation_flags = tfxAnimationFlags_needs_recording | tfxAnimationFlags_export_with_transparency;
		a.seed = 0;
		a.needs_exporting = 0;
		//a.camera_settings.camera_floor_height = -10.f;
		//a.camera_settings.camera_fov = tfxRadians(60);
		//a.camera_settings.camera_pitch = tfxRadians(-30.f);
		//a.camera_settings.camera_yaw = tfxRadians(-90.f);
		//a.camera_settings.camera_position = tfxVec3(0.f, 3.5f, 7.5f);
		//a.camera_settings.camera_isometric = false;
		//a.camera_settings.camera_isometric_scale = 5.f;
		//a.camera_settings.camera_hide_floor = false;
		sprite_data_settings.push_back(a);
		effect.GetInfo().sprite_data_settings_index = sprite_data_settings.size() - 1;
		return effect.GetInfo().sprite_data_settings_index;
	}

	tfxU32 tfxLibrary::AddPreviewCameraSettings(tfxEffectEmitter& effect) {
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

	tfxU32 tfxLibrary::AddPreviewCameraSettings() {
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

	tfxU32 tfxLibrary::AddEffectEmitterInfo() {
		tfxEffectEmitterInfo info;
		if (free_infos.size()) {
			return free_infos.pop_back();
		}
		effect_infos.push_back(info);
		return effect_infos.size() - 1;
	}

	tfxU32 tfxLibrary::AddEmitterProperties() {
		if (free_properties.size()) {
			return free_properties.pop_back();
		}
		return AddRow(&emitter_properties_buffer, true);
	}

	void tfxLibrary::Init() {
		graph_node_allocator = CreateArenaManager(tfxMegabyte(2), 8);
		graph_lookup_allocator = CreateArenaManager(tfxMegabyte(4), 8);
		InitEmitterProperties();
	}

	void tfxLibrary::InitEmitterProperties() {
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
		AddStructArray(&emitter_properties_buffer, sizeof(float), offsetof(tfxEmitterPropertiesSoA, noise_base_offset_range));
		FinishSoABufferSetup(&emitter_properties_buffer, &emitter_properties, 100);
	}

	void tfxLibrary::Clear() {
		for (auto &e : effects) {
			e.FreeGraphs();
		}
		effects.free_all();
		effect_paths.FreeAll();
		particle_shapes.FreeAll();
		global_graphs.free_all();
		emitter_attributes.free_all();
		transform_attributes.free_all();
		sprite_sheet_settings.free_all();
		preview_camera_settings.free_all();
		all_nodes.free_all();
		node_lookup_indexes.free_all();
		compiled_lookup_indexes.free_all();
		compiled_lookup_values.free_all();
		FreeSoABuffer(&emitter_properties_buffer);
		shape_data.free_all();
		graph_min_max.free_all();
		for (auto &info : effect_infos) {
			info.sub_effectors.free_all();
		}
		effect_infos.free_all();
		AddPreviewCameraSettings();
		pre_recorded_effects.FreeAll();

		graph_node_allocator.FreeAll();
		graph_lookup_allocator.FreeAll();
		sprite_data_allocator.FreeAll();

		free_global_graphs.free_all();
		free_keyframe_graphs.free_all();
		free_emitter_attributes.free_all();
		free_animation_settings.free_all();
		free_preview_camera_settings.free_all();
		free_infos.free_all();
		free_properties.free_all();
		free_keyframes.free_all();

		uid = 0;
	}

	void tfxLibrary::UpdateComputeNodes() {
		tfxU32 running_node_index = 0;
		tfxU32 running_value_index = 0;
		tmpStack(tfxEffectEmitter*, stack);
		all_nodes.clear();
		node_lookup_indexes.clear();
		compiled_lookup_values.clear();
		compiled_lookup_indexes.clear();
		for (auto &effect : effects) {
			stack.push_back(&effect);
			while (!stack.empty()) {
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

	void tfxLibrary::CompileGraphsOfEffect(tfxEffectEmitter &effect, tfxU32 depth) {
		auto &info = effect.GetInfo();
		if (effect.type == tfxEffectType && depth == 0) {
			CompileKeyframeGraph(effect.transform_attributes);
			CompileGlobalGraph(effect.global);
		}
		else if (effect.type == tfxEffectType) {
			CompileKeyframeGraph(effect.transform_attributes);
		}
		else if (effect.type == tfxEmitterType) {
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

	void tfxLibrary::CompileAllGraphs() {
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

	void tfxLibrary::CompileGlobalGraph(tfxU32 index) {
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

	void tfxLibrary::CompileKeyframeGraph(tfxU32 index) {
		tfxTransformAttributes &g = transform_attributes[index];
		CompileGraph(g.roll);
		CompileGraph(g.pitch);
		CompileGraph(g.yaw);
		CompileGraph(g.translation_x);
		CompileGraph(g.translation_y);
		CompileGraph(g.translation_z);
	}

	void tfxLibrary::CompileEmitterGraphs(tfxU32 index) {
		CompilePropertyGraph(index);
		CompileKeyframeGraph(index);
		CompileBaseGraph(index);
		CompileVariationGraph(index);
		CompileOvertimeGraph(index);
	}

	void tfxLibrary::CompilePropertyGraph(tfxU32 index) {
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
	void tfxLibrary::CompileBaseGraph(tfxU32 index) {
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
	void tfxLibrary::CompileVariationGraph(tfxU32 index) {
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
	void tfxLibrary::CompileOvertimeGraph(tfxU32 index) {
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
	void tfxLibrary::CompileColorGraphs(tfxU32 index) {
		tfxOvertimeAttributes &g = emitter_attributes[index].overtime;
		CompileGraphOvertime(g.red);
		CompileGraphOvertime(g.green);
		CompileGraphOvertime(g.blue);
	}

	void tfxLibrary::SetMinMaxData() {
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

	float tfxLibrary::LookupPreciseOvertimeNodeList(tfxGraphType graph_type, int lookup_node_index, float age, float life) {
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

	float tfxLibrary::LookupPreciseNodeList(tfxGraphType graph_type, int lookup_node_index, float age) {
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

	float tfxLibrary::LookupFastValueList(tfxGraphType graph_type, int lookup_node_index, float frame) {
		tfxGraphLookupIndex &lookup_data = ((tfxGraphLookupIndex*)&compiled_lookup_indexes[lookup_node_index])[graph_type];
		frame += lookup_data.start_index;
		tfxU32 end_frame = lookup_data.start_index + lookup_data.length - 1;
		frame = frame > end_frame ? end_frame : frame;
		return compiled_lookup_values[(tfxU32)frame];
	}

	float tfxLibrary::LookupFastOvertimeValueList(tfxGraphType graph_type, int lookup_value_index, float age, float lifetime) {
		tfxGraphLookupIndex &lookup_data = ((tfxGraphLookupIndex*)&compiled_lookup_indexes[lookup_value_index])[graph_type - tfxOvertime_velocity];
		float frame = (float)lookup_data.start_index;
		if (lifetime)
			frame += (age / lifetime * lookup_data.max_life) / tfxLOOKUP_FREQUENCY_OVERTIME;
		if (frame < lookup_data.start_index + lookup_data.length - 1)
			return compiled_lookup_values[(tfxU32)frame];
		return compiled_lookup_values[lookup_data.start_index + lookup_data.length - 1];
	}

	tfxU32 tfxLibrary::CountOfGraphsInUse() {
		return global_graphs.size() + emitter_attributes.size() - CountOfFreeGraphs();
	}

	tfxU32 tfxLibrary::CountOfFreeGraphs() {
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
		names_and_types.Insert("max_threads", tfxSInt);
		names_and_types.Insert("use_texture_filtering", tfxBool);
		initialised = true;
	}

	int ValidateEffectPackage(const char *filename) {
		tfxPackage package;
		tfxErrorFlags status = LoadPackage(filename, package);
		if (status) {
			package.Free();
			return status;					//returns 1 to 4 if it's an invalid package format
		}

		tfxEntryInfo *data_txt = package.GetFile("data.txt");
		package.Free();
		if (!data_txt) return tfxErrorCode_data_could_not_be_loaded;					//Unable to load the the data.txt file in the package

		return 0;
	}

	void AssignGraphData(tfxEffectEmitter &effect, tfxStack<tfxStr256> &values) {
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

	void AssignNodeData(tfxAttributeNode &n, tfxStack<tfxStr256> &values) {
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
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].frames = value;
		if (field == "current_frame")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].current_frame = value;
		if (field == "seed")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].seed = value;
		if (field == "layer")
			emitter_properties.layer[effect.property_index] = value >= tfxLAYERS ? value = tfxLAYERS - 1 : value;
		if (field == "frame_offset")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].frame_offset = value;
		if (field == "single_shot_limit")
			emitter_properties.single_shot_limit[effect.property_index] = value;
		if (field == "billboard_option")
			emitter_properties.billboard_option[effect.property_index] = (tfxBillboardingOptions)value;
		if (field == "vector_align_type")
			emitter_properties.vector_align_type[effect.property_index] = value >= 0 && value < tfxVectorAlignType_max ? (tfxVectorAlignType)value : (tfxVectorAlignType)0;
		if (field == "angle_setting")
			emitter_properties.angle_settings[effect.property_index] = (tfxAngleSettingFlags)value;
		if (field == "sort_passes")
			effect.sort_passes = tfxMin(5, value);
		if (field == "animation_flags")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].animation_flags = value;
		if (field == "sprite_data_flags")
			effect.library->sprite_data_settings[effect.GetInfo().sprite_data_settings_index].animation_flags = value;
		if (field == "sprite_data_seed")
			effect.library->sprite_data_settings[effect.GetInfo().sprite_data_settings_index].seed = value;
		if (field == "sprite_data_frame_offset")
			effect.library->sprite_data_settings[effect.GetInfo().sprite_data_settings_index].frame_offset = value;
		if (field == "sprite_data_frames")
			effect.library->sprite_data_settings[effect.GetInfo().sprite_data_settings_index].frames = value;
		if (field == "sprite_data_extra_frames_count")
			effect.library->sprite_data_settings[effect.GetInfo().sprite_data_settings_index].extra_frames_count = value;
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, int value) {
		tfxEmitterPropertiesSoA &emitter_properties = effect.library->emitter_properties;
		if (field == "emission_type")
			emitter_properties.emission_type[effect.property_index] = (tfxEmissionType)value;
		if (field == "emission_direction")
			emitter_properties.emission_direction[effect.property_index] = (tfxEmissionDirection)value;
		if (field == "color_option")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].color_option = value > 3 ? tfxFullColor : (tfxExportColorOptions)value;
		if (field == "export_option")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].export_option = (tfxExportOptions)value;
		if (field == "end_behaviour")
			emitter_properties.end_behaviour[effect.property_index] = (tfxLineTraversalEndBehaviour)value;
		if (field == "frame_offset")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].frame_offset = value;
		if (field == "extra_frames_count")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].extra_frames_count = value;
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, tfxStr &value) {
		if (field == "name") {
			effect.GetInfo().name = value;
		}
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, float value) {
		tfxEmitterPropertiesSoA &emitter_properties = effect.library->emitter_properties;
		if (field == "position_x")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].position.x = value;
		if (field == "position_y")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].position.y = value;
		if (field == "position_z")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].position.z = value;
		if (field == "frame_width")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].frame_size.x = value;
		if (field == "frame_height")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].frame_size.y = value;
		if (field == "zoom")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].zoom = value;
		if (field == "scale")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].scale = value;
		if (field == "playback_speed")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].playback_speed = value;
		if (field == "camera_position_x")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings.camera_position.x = value;
		if (field == "camera_position_y")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings.camera_position.y = value;
		if (field == "camera_position_z")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings.camera_position.z = value;
		if (field == "camera_pitch")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings.camera_pitch = value;
		if (field == "camera_yaw")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings.camera_yaw = value;
		if (field == "camera_fov")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings.camera_fov = value;
		if (field == "camera_floor_height")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings.camera_floor_height = value;
		if (field == "camera_isometric_scale")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings.camera_isometric_scale = value;
		if (field == "orthographic_camera_position_x")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings_orthographic.camera_position.x = value;
		if (field == "orthographic_camera_position_y")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings_orthographic.camera_position.y = value;
		if (field == "orthographic_camera_position_z")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings_orthographic.camera_position.z = value;
		if (field == "orthographic_camera_pitch")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings_orthographic.camera_pitch = value;
		if (field == "orthographic_camera_yaw")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings_orthographic.camera_yaw = value;
		if (field == "orthographic_camera_fov")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings_orthographic.camera_fov = value;
		if (field == "orthographic_camera_floor_height")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings_orthographic.camera_floor_height = value;
		if (field == "orthographic_camera_isometric_scale")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings_orthographic.camera_isometric_scale = value;
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
		if (field == "noise_base_offset_range")
			emitter_properties.noise_base_offset_range[effect.property_index] = value < 0 ? 0.f : value;
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
		if (field == "sprite_data_playback_speed")
			effect.library->sprite_data_settings[effect.GetInfo().sprite_data_settings_index].playback_speed = value;
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, bool value) {
		if (field == "loop")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].animation_flags |= value ? tfxAnimationFlags_loop : 0;
		if (field == "seamless")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].animation_flags |= value ? tfxAnimationFlags_seamless : 0;
		if (field == "export_with_transparency")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].animation_flags |= value ? tfxAnimationFlags_export_with_transparency : 0;
		if (field == "camera_isometric")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings.camera_isometric = false;
		if (field == "camera_hide_floor")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings.camera_hide_floor = value;
		if (field == "orthographic_camera_isometric")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings_orthographic.camera_isometric = true;
		if (field == "orthographic_camera_hide_floor")
			effect.library->sprite_sheet_settings[effect.GetInfo().sprite_sheet_settings_index].camera_settings_orthographic.camera_hide_floor = value;
		if (field == "preview_attach_effect_to_camera")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].attach_effect_to_camera = value;
		if (field == "preview_camera_hide_floor")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_hide_floor = value;
		if (field == "preview_camera_isometric")
			effect.library->preview_camera_settings[effect.GetInfo().preview_camera_settings].camera_settings.camera_isometric = value;
		if (field == "random_color")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_random_color; else effect.property_flags &= ~tfxEmitterPropertyFlags_random_color;
		if (field == "relative_position")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_relative_position; else effect.property_flags &= ~tfxEmitterPropertyFlags_relative_position;
		if (field == "relative_angle")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_relative_angle; else effect.property_flags &= ~tfxEmitterPropertyFlags_relative_angle;
		if (field == "image_handle_auto_center")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_image_handle_auto_center; else effect.property_flags &= ~tfxEmitterPropertyFlags_image_handle_auto_center;
		if (field == "single")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_single; else effect.property_flags &= ~tfxEmitterPropertyFlags_single;
		//if (field == "one_shot")
			//if(value) effect.property_flags |= tfxEmitterPropertyFlags_one_shot; else effect.property_flags &= ~tfxEmitterPropertyFlags_one_shot;
		if (field == "spawn_on_grid")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_spawn_on_grid; else effect.property_flags &= ~tfxEmitterPropertyFlags_spawn_on_grid;
		if (field == "grid_spawn_clockwise")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_grid_spawn_clockwise; else effect.property_flags &= ~tfxEmitterPropertyFlags_grid_spawn_clockwise;
		if (field == "fill_area")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_fill_area; else effect.property_flags &= ~tfxEmitterPropertyFlags_fill_area;
		if (field == "grid_spawn_random")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_grid_spawn_random; else effect.property_flags &= ~tfxEmitterPropertyFlags_grid_spawn_random;
		if (field == "area_open_ends")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_area_open_ends; else effect.property_flags &= ~tfxEmitterPropertyFlags_area_open_ends;
		if (field == "emitter_handle_auto_center")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_emitter_handle_auto_center; else effect.property_flags &= ~tfxEmitterPropertyFlags_emitter_handle_auto_center;
		if (field == "edge_traversal")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_edge_traversal; else effect.property_flags &= ~tfxEmitterPropertyFlags_edge_traversal;
		if (field == "image_reverse_animation")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_reverse_animation; else effect.property_flags &= ~tfxEmitterPropertyFlags_reverse_animation;
		if (field == "image_play_once")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_play_once; else effect.property_flags &= ~tfxEmitterPropertyFlags_play_once;
		if (field == "image_animate")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_animate; else effect.property_flags &= ~tfxEmitterPropertyFlags_animate;
		if (field == "image_random_start_frame")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_random_start_frame; else effect.property_flags &= ~tfxEmitterPropertyFlags_random_start_frame;
		if (field == "global_uniform_size")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_global_uniform_size; else effect.property_flags &= ~tfxEmitterPropertyFlags_global_uniform_size;
		if (field == "base_uniform_size")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_base_uniform_size; else effect.property_flags &= ~tfxEmitterPropertyFlags_base_uniform_size;
		if (field == "lifetime_uniform_size")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_lifetime_uniform_size; else effect.property_flags &= ~tfxEmitterPropertyFlags_lifetime_uniform_size;
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
		file.AddLine("noise_base_offset_range=%f", property.noise_base_offset_range[index]);
		file.AddLine("emitter_handle_x=%f", property.emitter_handle[index].x);
		file.AddLine("emitter_handle_y=%f", property.emitter_handle[index].y);
		file.AddLine("emitter_handle_z=%f", property.emitter_handle[index].z);
		file.AddLine("end_behaviour=%i", property.end_behaviour[index]);
		file.AddLine("random_color=%i", (flags & tfxEmitterPropertyFlags_random_color));
		file.AddLine("relative_position=%i", (flags & tfxEmitterPropertyFlags_relative_position));
		file.AddLine("relative_angle=%i", (flags & tfxEmitterPropertyFlags_relative_angle));
		file.AddLine("single=%i", (flags & tfxEmitterPropertyFlags_single));
		file.AddLine("single_shot_limit=%i", property.single_shot_limit[index]);
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

	static bool CompareNodes(tfxAttributeNode &left, tfxAttributeNode &right) {
		return left.frame < right.frame;
	}

	bool tfxGraph::IsOvertimeGraph() {
		return type >= tfxOvertime_velocity && type <= tfxOvertime_noise_resolution && type != tfxOvertime_velocity_adjuster;
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
			ClampCurve(graph, node.right, node);
			ClampCurve(graph, node.left, node);
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
			ClampCurve(graph, node.right, node);
			ClampCurve(graph, node.left, node);
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
				for (auto &n : nodes) {
					if (node.frame < n.frame)
						last_node = &n;
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

	float tfxGraph::GetRandomValue(float age, tfxRandom &random) {
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
						return random.Range(bezier_value);
					}
					else {
						return random.Range(lastv - p * (lastv - a.value));
					}
				}
				lastv = a.value;
				lastf = a.frame - 1;
				lastec = &a;
			}
		} while (!nodes.EndOfBuckets());
		return random.Range(lastv);

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
			for (auto &n : nodes) {
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
		nodes.TrimBuckets();
		if (add_node && preset == tfxWeightOvertimePreset) {
			AddNode(0.f, 0.f, 0);
			tfxAttributeNode *node = AddNode(1.f, 1.f, tfxAttributeNodeFlags_is_curve, 0.f, 1.f, 1.f, 1.f);
			node->SetCurveInitialised();
		}
		else if (add_node) {
			AddNode(0.f, v);
		}
		switch (preset) {
		case tfxGraphPreset::tfxGlobalPercentPreset:
			//We have a epsilon to prevent divide by 0 here
			min = { 0.f, 0.0001f }; max = { tfxMAX_FRAME, 20.f };
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
			//We have a epsilon to prevent divide by 0 here. The divide by zero occurrs in control functions (ControlParticleImageFrame3d etc.) when the current % life of the particle is calculated
			min = { 0.f, 0.0001f }; max = { tfxMAX_FRAME, 100000.f };
			break;
		case tfxGraphPreset::tfxAmountPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 5000.f };
			break;
		case tfxGraphPreset::tfxVelocityPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 10000.f };
			break;
		case tfxGraphPreset::tfxWeightPreset:
			min = { 0.f, -10000.f }; max = { tfxMAX_FRAME, 10000.f };
			break;
		case tfxGraphPreset::tfxWeightVariationPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 20000.f };
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
			min = { 0.f, -20.f }; max = { 1.f, 20.f };
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
			mm = { 0.f, -10000.f, tfxMAX_FRAME, 10000.f };
			break;
		case tfxGraphPreset::tfxWeightVariationPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 20000.f };
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
			mm = { 0.f, -20.f, 1.f, 20.f };
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

	void tfxGraph::ClearToOne(float value) {
		nodes.clear();
		AddNode(0.f, value);
	}

	void tfxGraph::Clear() {
		nodes.clear();
	}

	void tfxGraph::Free() {
		//Explicitly free the nodes
		nodes.free_all();
		lookup.values.free();
	}

	void tfxGraph::Copy(tfxGraph &to, bool compile) {
		to.Clear();
		do {
			for (auto &n : nodes) {
				to.nodes.push_back(n);
			}
		} while (!nodes.EndOfBuckets());
		if (compile) {
			if (IsOvertimeGraph())
				CompileGraphOvertime(to);
			else
				CompileGraph(to);
		}
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
			return tfxVec2(10.f, 0.01f);
			break;
		case tfxGraphPreset::tfxWeightPreset:
		case tfxGraphPreset::tfxWeightVariationPreset:
			return tfxVec2(10.f, 1.f);
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

	float GetRandomFast(tfxGraph &graph, float frame, tfxRandom &random) {
		float value = 0;
		if ((tfxU32)frame < graph.lookup.last_frame)
			value = graph.lookup.values[(tfxU32)frame];
		value = graph.lookup.values[graph.lookup.last_frame];
		return random.Range(value);
	}

	float GetRandomPrecise(tfxGraph &graph, float frame, tfxRandom &random) {
		return graph.GetRandomValue(frame, random);
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

	bool IsEmitterGraph(tfxGraphType type) {
		return type >= tfxPropertyStart && type < tfxTransformStart;
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
		tfxU32 size = keyframes.translation_x.nodes.size() +
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
		tfxVec3 point(lookup_callback(keyframes.translation_x, frame),
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
		if (!tfxDataTypes.initialised) tfxDataTypes.Init();
		FILE* fp;
		errno_t error = fopen_s(&fp, path, "r");
		if (error != 0) {
			printf("strerror says open failed: %s\n", strerror(error));
		}
		if (fp == NULL) {
			return false;
		}

		const size_t max_line_length = 512;
		char buffer[max_line_length];

		tfxvec<tfxStr256> pair;
		while (fgets(buffer, max_line_length, fp)) {
			buffer[strcspn(buffer, "\n")] = 0;
			tfxStr512 str = buffer;
			pair.clear();
			SplitStringVec(str, pair, 61);
			if (pair.size() == 2) {
				tfxStr256 key = pair[0];
				if (tfxDataTypes.names_and_types.ValidName(pair[0])) {
					tfxDataType t = tfxDataTypes.names_and_types.At(pair[0]);
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

		fclose(fp);
		return true;

	}

	void SplitStringStack(const tfxStr &str, tfxStack<tfxStr256> &pair, char delim) {
		tfxStr256 line;
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

	void SplitStringVec(const tfxStr &str, tfxvec<tfxStr256> &pair, char delim) {
		tfxStr256 line;
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

		if (s.IsInt())
			return tfxSInt;

		if (s.IsFloat())
			return tfxFloat;

		return tfxString;
	}

	//Get a graph by tfxGraphID
	tfxGraph &GetGraph(tfxLibrary &library, tfxGraphID &graph_id) {
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
		tmpStack(tfxStr256, pair);

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

		package.Free();
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

		tmpStack(tfxStr256, pair);
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

	tfxEffectLibraryStats CreateLibraryStats(tfxLibrary &lib) {
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
				else if (current.type == tfxEmitterType) {
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

	tfxErrorFlags LoadEffectLibraryPackage(tfxPackage &package, tfxLibrary &lib, void(*shape_loader)(const char *filename, tfxImageData &image_data, void *raw_image_data, int image_size, void *user_data), void *user_data, bool read_only) {

		assert(shape_loader);			//Must have a shape_loader function to load your shapes with. This will be a custom user function suited for whichever renderer you're using
		if (!tfxDataTypes.initialised)
			tfxDataTypes.Init();
		lib.Clear();
		if (tfxIcospherePoints[0].current_size == 0) {
			MakeIcospheres();
		}

		tfxEntryInfo *data = package.GetFile("data.txt");
		tfxEntryInfo *stats_struct = package.GetFile("stats.struct");
		tfxErrorFlags error = 0;

		int context = 0;
		int uid = 0;
		tfxU32 current_global_graph = 0;

		lib.InitEmitterProperties();
		lib.sprite_data_allocator = CreateArenaManager(tfxMegabyte(4));

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
		tmpStack(tfxEffectEmitter, effect_stack);
		tmpStack(tfxStr256, pair);

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
					lib.AddSpriteSheetSettings(effect);
					lib.AddSpriteDataSettings(effect);
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
						continue;
					}
				}

				if (context == tfxStartAnimationSettings || context == tfxStartEmitter || context == tfxStartEffect || context == tfxStartFolder || context == tfxStartPreviewCameraSettings) {
					if (tfxDataTypes.names_and_types.ValidName(pair[0])) {
						switch (tfxDataTypes.names_and_types.At(pair[0])) {
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
					if (tfxDataTypes.names_and_types.ValidName(pair[0])) {
						switch (tfxDataTypes.names_and_types.At(pair[0])) {
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
				if (effect_stack.back().property_flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
					lib.emitter_properties.image_handle[effect_stack.back().property_index] = { .5f, .5f };
				}
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
					else if (effect_stack.parent().type == tfxEmitterType) {
						effect_stack.back().global = current_global_graph;
					}
					effect_stack.parent().GetInfo().sub_effectors.push_back(effect_stack.back());
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

	tfxErrorFlags LoadEffectLibraryPackage(const char *filename, tfxLibrary &lib, void(*shape_loader)(const char* filename, tfxImageData &image_data, void *raw_image_data, int image_size, void *user_data), void *user_data, bool read_only) {

		tfxErrorFlags error = 0;

		tfxPackage package;
		error = LoadPackage(filename, package);
		if (error != 0) {
			package.Free();
			return error;
		}
		error = LoadEffectLibraryPackage(package, lib, shape_loader, user_data, read_only);

		package.Free();
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

	void InvalidateNewSpriteCapturedIndex(tfxParticleManager &pm) {
		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			tfxSpriteSoA &sprites = pm.sprites[pm.current_sprite_buffer][layer];
			for (int i = 0; i != pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size; ++i) {
				if ((sprites.captured_index[i] & 0xF0000000) >> 28 == pm.current_sprite_buffer) {
					sprites.captured_index[i] = tfxINVALID;
				}
			}
		}
	}

	void RecordSpriteData2d(tfxParticleManager &pm, tfxEffectEmitter &effect) {

	}

	void RecordSpriteData3d(tfxParticleManager &pm, tfxEffectEmitter &effect) {
		tfxSpriteDataSettings &anim = effect.library->sprite_data_settings[effect.GetInfo().sprite_data_settings_index];
		tfxU32 frames = anim.frames;
		tfxU32 start_frame = anim.frame_offset;
		int extra_frames = anim.extra_frames_count;
		tfxU32 frame = 0;
		int extra_frame_count = 0;
		tfxU32 offset = 0;
		bool particles_started = false;
		bool start_counting_extra_frames = false;

		float update_freq = tfxUPDATE_FREQUENCY;
		//SetUpdateFrequency(60.f * (anim.playback_speed ? anim.playback_speed : 1.f));

		bool auto_set_length = false;
		if (anim.animation_flags & tfxAnimationFlags_auto_set_length && !(anim.animation_flags & tfxAnimationFlags_loop) && effect.IsFinite()) {
			frames = 99999;
			auto_set_length = true;
		}

		//First pass to count the number of sprites in each frame
		//pm.UpdateAgeOnly(false);
		//pm.ForceSingleThreaded(true);

		pm.ClearAll();
		SetSeed(&pm, anim.seed);
		tfxU32 preview_effect_index = pm.AddEffect(effect, pm.current_ebuff);
		SetEffectPosition(&pm, preview_effect_index, tfxVec3(0.f, 0.f, 0.f));
		Transform3d(pm.effects.world_rotations[preview_effect_index],
			pm.effects.local_rotations[preview_effect_index],
			pm.effects.scale[preview_effect_index],
			pm.effects.world_position[preview_effect_index],
			pm.effects.local_position[preview_effect_index],
			pm.effects.translation[preview_effect_index],
			pm.effects.matrix[preview_effect_index],
			pm.effects.world_rotations[preview_effect_index],
			pm.effects.scale[preview_effect_index],
			pm.effects.world_position[preview_effect_index],
			pm.effects.matrix[preview_effect_index]
		);

		tfxU32 total_sprites = 0;
		tfxvec<tfxFrameMeta> tmp_frame_meta;
		tfxU32 sprites_in_layers = 0;
		while (frame < frames && offset < 99999) {
			tfxU32 count_this_frame = 0;
			pm.Update();
			bool particles_processed_last_frame = false;

			if (offset >= start_frame) {
				sprites_in_layers = 0;
				for (tfxEachLayer) {
					if (frame >= tmp_frame_meta.size()) {
						tfxFrameMeta meta;
						memset(&meta, 0, sizeof(tfxFrameMeta));
						tmp_frame_meta.push_back(meta);
					}
					tmp_frame_meta[frame].sprite_count[layer] += pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size;
					total_sprites += pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size;
					sprites_in_layers += pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size;
					particles_started = total_sprites > 0;
					particles_processed_last_frame |= pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size > 0;
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
			if (start_counting_extra_frames && extra_frame_count++ >= extra_frames)
				pm.DisableSpawning(true);
		}

		frames = tmp_frame_meta.size();

		tfxSpriteData *sprite_data = nullptr;
		if (effect.library->pre_recorded_effects.ValidKey(effect.path_hash)) {
			sprite_data = &effect.library->pre_recorded_effects.At(effect.path_hash);
			FreeSpriteData(*sprite_data);
			sprite_data->frame_count = frames;
			sprite_data->animation_length_in_time = frames * tfxFRAME_LENGTH;
			sprite_data->frame_meta = tfxArray<tfxFrameMeta>(&effect.library->sprite_data_allocator, frames);
			sprite_data->frame_meta.zero();
		}
		else {
			tfxSpriteData data;
			data.frame_count = frames;
			data.animation_length_in_time = frames * tfxFRAME_LENGTH;
			data.frame_meta = tfxArray<tfxFrameMeta>(&effect.library->sprite_data_allocator, frames);
			data.frame_meta.zero();
			effect.library->pre_recorded_effects.Insert(effect.path_hash, data);
			sprite_data = &effect.library->pre_recorded_effects.At(effect.path_hash);
		}

		anim.frames = frames;
		anim.animation_time = sprite_data->animation_length_in_time;
		sprite_data->frame_compression = anim.playback_speed;

		tfxArray<tfxFrameMeta> &frame_meta = sprite_data->frame_meta;
		memcpy(frame_meta.block, tmp_frame_meta.data, tmp_frame_meta.size_in_bytes());
		tmp_frame_meta.free_all();

		tfxU32 last_count = 0;
		for (auto &meta : frame_meta) {
			for (tfxEachLayer) {
				meta.index_offset[layer] = last_count;
				last_count += meta.sprite_count[layer];
			}
		}

		//std::cout << "Total Sprites: " << total_sprites << std::endl;

		//for (auto &meta : frame_meta) {
			//std::cout << meta.index_offset[0] << ", " << meta.sprite_count[0] << std::endl;
		//}

		//pm.UpdateAgeOnly(false);

		pm.ClearAll();
		SetSeed(&pm, anim.seed);
		preview_effect_index = pm.AddEffect(effect, pm.current_ebuff);
		SetEffectPosition(&pm, preview_effect_index, tfxVec3(0.f, 0.f, 0.f));
		Transform3d(pm.effects.world_rotations[preview_effect_index],
			pm.effects.local_rotations[preview_effect_index],
			pm.effects.scale[preview_effect_index],
			pm.effects.world_position[preview_effect_index],
			pm.effects.local_position[preview_effect_index],
			pm.effects.translation[preview_effect_index],
			pm.effects.matrix[preview_effect_index],
			pm.effects.world_rotations[preview_effect_index],
			pm.effects.scale[preview_effect_index],
			pm.effects.world_position[preview_effect_index],
			pm.effects.matrix[preview_effect_index]
		);

		sprite_data->total_sprites = total_sprites;
		sprite_data->total_memory_for_sprites = total_sprites * sizeof(tfxSprite3d); 
		InitSpriteData3dSoACompression(&sprite_data->real_time_sprites_buffer, &sprite_data->real_time_sprites, total_sprites);
		InitSpriteData3dSoA(&sprite_data->compressed_sprites_buffer, &sprite_data->compressed_sprites, total_sprites);

		tfxSoABuffer temp_sprites_buffer;
		tfxSpriteDataSoA temp_sprites;
		InitSpriteData3dSoA(&temp_sprites_buffer, &temp_sprites, 100);
		tfxvec<tfxU32> running_count[tfxLAYERS];

		for (tfxEachLayer) {
			running_count[layer].resize(frames);
			running_count[layer].zero();
		}

		frame = 0;
		offset = 0;
		extra_frame_count = 0;
		particles_started = false;
		start_counting_extra_frames = false;
		pm.DisableSpawning(false);
		total_sprites = 0;
		tfxU32 captured_offset[tfxLAYERS] = { 0, 0, 0, 0 };

		while (frame < frames && offset < 99999) {
			tfxU32 count_this_frame = 0;
			pm.Update();
			InvalidateNewSpriteCapturedIndex(pm);
			bool particles_processed_last_frame = false;

			if (offset >= start_frame) {
				for (tfxEachLayer) {
					tfxU32 meta_count = frame_meta[frame].sprite_count[layer];
					tfxU32 pm_count = pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size;
					//assert(frame_meta[frame].sprite_count[layer] == pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size);
					if (running_count[layer][frame] > 0 && pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size > 0) {
						Resize(&temp_sprites_buffer, running_count[layer][frame]);
						memcpy(temp_sprites.alignment, sprite_data->real_time_sprites.alignment + frame_meta[frame].index_offset[layer], sizeof(tfxU32) * running_count[layer][frame]);
						memcpy(temp_sprites.captured_index, sprite_data->real_time_sprites.captured_index + frame_meta[frame].index_offset[layer], sizeof(tfxU32) * running_count[layer][frame]);
						memcpy(temp_sprites.color, sprite_data->real_time_sprites.color + frame_meta[frame].index_offset[layer], sizeof(tfxU32) * running_count[layer][frame]);
						memcpy(temp_sprites.image_frame_plus, sprite_data->real_time_sprites.image_frame_plus + frame_meta[frame].index_offset[layer], sizeof(tfxU32) * running_count[layer][frame]);
						memcpy(temp_sprites.intensity, sprite_data->real_time_sprites.intensity + frame_meta[frame].index_offset[layer], sizeof(float) * running_count[layer][frame]);
						memcpy(temp_sprites.stretch, sprite_data->real_time_sprites.stretch + frame_meta[frame].index_offset[layer], sizeof(float) * running_count[layer][frame]);
						memcpy(temp_sprites.transform_3d, sprite_data->real_time_sprites.transform_3d + frame_meta[frame].index_offset[layer], sizeof(tfxSpriteTransform3d) * running_count[layer][frame]);
						if (captured_offset[layer] > 0) {
							for (int temp_i = 0; temp_i != temp_sprites_buffer.current_size; ++temp_i) {
								if(temp_sprites.captured_index[temp_i] != tfxINVALID)
									temp_sprites.captured_index[temp_i] += captured_offset[layer];
							}
						}
					}
					else if (captured_offset[layer] > 0 && pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size == 0) {
						for (int index = SpriteDataIndexOffset(sprite_data, frame, layer); index != SpriteDataEndIndex(sprite_data, frame, layer); ++index) {
							if(sprite_data->real_time_sprites.captured_index[index] != tfxINVALID)
								sprite_data->real_time_sprites.captured_index[index] += captured_offset[layer];
						}
					}
					memcpy(sprite_data->real_time_sprites.alignment + frame_meta[frame].index_offset[layer], pm.sprites[pm.current_sprite_buffer][layer].alignment, sizeof(tfxU32) * pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size);
					memcpy(sprite_data->real_time_sprites.captured_index + frame_meta[frame].index_offset[layer], pm.sprites[pm.current_sprite_buffer][layer].captured_index, sizeof(tfxU32) * pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size);
					memcpy(sprite_data->real_time_sprites.color + frame_meta[frame].index_offset[layer], pm.sprites[pm.current_sprite_buffer][layer].color, sizeof(tfxU32) * pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size);
					memcpy(sprite_data->real_time_sprites.image_frame_plus + frame_meta[frame].index_offset[layer], pm.sprites[pm.current_sprite_buffer][layer].image_frame_plus, sizeof(tfxU32) * pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size);
					memcpy(sprite_data->real_time_sprites.intensity + frame_meta[frame].index_offset[layer], pm.sprites[pm.current_sprite_buffer][layer].intensity, sizeof(float) * pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size);
					memcpy(sprite_data->real_time_sprites.stretch + frame_meta[frame].index_offset[layer], pm.sprites[pm.current_sprite_buffer][layer].stretch, sizeof(float) * pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size);
					memcpy(sprite_data->real_time_sprites.transform_3d + frame_meta[frame].index_offset[layer], pm.sprites[pm.current_sprite_buffer][layer].transform_3d, sizeof(tfxSpriteTransform3d) * pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size);
					if (running_count[layer][frame] > 0 && pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size > 0) {
						memcpy(sprite_data->real_time_sprites.alignment + frame_meta[frame].index_offset[layer] + pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size, temp_sprites.alignment, sizeof(tfxU32) * temp_sprites_buffer.current_size);
						memcpy(sprite_data->real_time_sprites.captured_index + frame_meta[frame].index_offset[layer] + pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size, temp_sprites.captured_index, sizeof(tfxU32) * temp_sprites_buffer.current_size);
						memcpy(sprite_data->real_time_sprites.color + frame_meta[frame].index_offset[layer] + pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size, temp_sprites.color, sizeof(tfxU32) * temp_sprites_buffer.current_size);
						memcpy(sprite_data->real_time_sprites.image_frame_plus + frame_meta[frame].index_offset[layer] + pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size, temp_sprites.image_frame_plus, sizeof(tfxU32) * temp_sprites_buffer.current_size);
						memcpy(sprite_data->real_time_sprites.intensity + frame_meta[frame].index_offset[layer] + pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size, temp_sprites.intensity, sizeof(float) * temp_sprites_buffer.current_size);
						memcpy(sprite_data->real_time_sprites.stretch + frame_meta[frame].index_offset[layer] + pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size, temp_sprites.stretch, sizeof(float) * temp_sprites_buffer.current_size);
						memcpy(sprite_data->real_time_sprites.transform_3d + frame_meta[frame].index_offset[layer] + pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size, temp_sprites.transform_3d, sizeof(tfxSpriteTransform3d) * temp_sprites_buffer.current_size);
						captured_offset[layer] = pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size;
					}
					else if (pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size == 0) {
						captured_offset[layer] = 0;
					}
					running_count[layer][frame] += pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size;
					total_sprites += pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size;
					particles_started = total_sprites > 0;
					particles_processed_last_frame |= pm.sprite_buffer[pm.current_sprite_buffer][layer].current_size > 0;
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
			if (start_counting_extra_frames && extra_frame_count++ >= extra_frames)
				pm.DisableSpawning(true);
		}

		FreeSoABuffer(&temp_sprites_buffer);
		SetUpdateFrequency(update_freq);
		pm.DisableSpawning(false);
		//pm.ForceSingleThreaded(false);

	}

	void CompressSpriteData3d(tfxParticleManager &pm, tfxEffectEmitter &effect) {
		tfxSpriteDataSettings &anim = effect.library->sprite_data_settings[effect.GetInfo().sprite_data_settings_index];
		tfxSpriteData *sprite_data = &effect.library->pre_recorded_effects.At(effect.path_hash);

		sprite_data->compressed_frame_meta.free();
		sprite_data->compressed_frame_meta = tfxArray<tfxFrameMeta>(&effect.library->sprite_data_allocator, anim.frames * anim.playback_speed);
		sprite_data->compressed_frame_meta.zero();

		float frequency = tfxFRAME_LENGTH * (anim.playback_speed ? anim.playback_speed : 1.f);
		float real_time = 0.f;
		float compressed_time = 0.f;
		int compressed_frame = 0;

		int f = 0;
		tfxSpriteDataSoA &sprites = sprite_data->real_time_sprites;
		tfxSpriteDataSoA &c_sprites = sprite_data->compressed_sprites;
		int ci = 0;
		while (f < anim.frames) {
			real_time = f * tfxFRAME_LENGTH;
			int next_frame = f + 1;
			int next_compressed_frame = compressed_frame++;
			float next_compressed_time = next_compressed_frame * frequency;
			tfxU32 frame_pair[2];
			frame_pair[0] = next_frame >= sprite_data->frame_count ? 0 : next_frame;
			frame_pair[1] = next_frame >= sprite_data->frame_count ? sprite_data->frame_count - 1 : next_frame - 1;
			for (tfxEachLayer) {
				if (real_time >= compressed_time && sprite_data->compressed_frame_meta[compressed_frame].sprite_count[layer] == 0) {
					for (int i = SpriteDataIndexOffset(sprite_data, f, layer); i != SpriteDataEndIndex(sprite_data, f, layer); ++i) {
						//Add to compress sprites but make invalid captured indexed create the offset
						sprite_data->compressed_frame_meta[compressed_frame].sprite_count[layer]++;
						c_sprites.alignment[ci] = sprites.alignment[i];
						c_sprites.captured_index[ci] = sprites.captured_index[i];
						c_sprites.color[ci] = sprites.color[i];
						c_sprites.image_frame_plus[ci] = sprites.image_frame_plus[i];
						c_sprites.intensity[ci] = sprites.intensity[i];
						c_sprites.stretch[ci] = sprites.stretch[i];
						c_sprites.transform_3d[ci] = sprites.transform_3d[i];
						ci++;
					}
					f++;
				}
				else if (real_time > compressed_time && real_time < next_compressed_time) {
					for (int i = SpriteDataIndexOffset(sprite_data, f, layer); i != SpriteDataEndIndex(sprite_data, f, layer); ++i) {
						if (sprites.captured_index[i] == tfxINVALID) {
							//Add to compressed sprites frame but add the lerp offset
							sprite_data->compressed_frame_meta[compressed_frame].sprite_count[layer]++;
							c_sprites.alignment[ci] = sprites.alignment[i];
							c_sprites.captured_index[ci] = tfxINVALID;
							c_sprites.color[ci] = sprites.color[i];
							c_sprites.image_frame_plus[ci] = sprites.image_frame_plus[i];
							c_sprites.intensity[ci] = sprites.intensity[i];
							c_sprites.stretch[ci] = sprites.stretch[i];
							c_sprites.transform_3d[ci] = sprites.transform_3d[i];
							ci++;
						}
					}
					f++;
				}
				else {
					compressed_time += frequency;
					compressed_frame++;
				}
			}
		}

	}

	tfxAPI void tfxEffectTemplate::RecordSpriteData(tfxParticleManager &pm) {
		if (effect.Is3DEffect()) {
			RecordSpriteData3d(pm, effect);
		}
		else {
			RecordSpriteData2d(pm, effect);
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
		assert(effect.library == library);	//The effect must belong to the same library that is assigned to the particle manager
		if (flags & tfxEffectManagerFlags_use_compute_shader && highest_compute_controller_index >= max_compute_controllers && free_compute_controllers.empty())
			return tfxINVALID;
		unsigned int parent_index = GetEffectSlot();
		if (parent_index == tfxINVALID)
			return tfxINVALID;
		if (!is_sub_emitter) {
			effects.highest_particle_age[parent_index] = tfxFRAME_LENGTH * 3.f;
		}
		tfxEmitterPropertiesSoA &properties = effect.library->emitter_properties;
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
		effects.user_data[parent_index] = effect.user_data;
		effects.update_callback[parent_index] = effect.update_callback;
		float range = properties.noise_base_offset_range[effect.property_index];
		effects.noise_base_offset[parent_index] = random.Range(range);
		effects_in_use[hierarchy_depth][buffer].push_back(parent_index);
		effect.pm_index = parent_index;
		sort_passes = tfxMax(effect.sort_passes, sort_passes);
		sort_passes = tfxMin(5, sort_passes);

		tfxU32 seed_index = 0;
		for (auto &e : effect.GetInfo().sub_effectors) {
			if (e.property_flags & tfxEmitterPropertyFlags_enabled) {
				unsigned int index = GetEmitterSlot();
				if (index == tfxINVALID)
					break;
				emitters.particles_index[index] = tfxINVALID;
				emitters_in_use[hierarchy_depth][buffer].push_back(index);
				emitters.parent_index[index] = parent_index;
				if (emitters.particles_index[index] == tfxINVALID) {
					if (!is_sub_emitter)
						emitters.particles_index[index] = GrabParticleLists(*this, e.path_hash, 100);
				}
				emitters.path_hash[index] = e.path_hash;
				emitters.info_index[index] = e.info_index;
				emitters.properties_index[index] = e.property_index;
				emitters.emitter_attributes[index] = e.emitter_attributes;
				emitters.transform_attributes[index] = e.transform_attributes;
				emitters.delay_spawning[index] = properties.delay_spawning[e.property_index];
				emitters.age[index] = 0.f;
				emitters.frame[index] = 0.f;
				emitters.local_position[index] = tfxVec3();
				emitters.grid_coords[index] = tfxVec3();
				emitters.grid_direction[index] = tfxVec3();
				emitters.property_flags[index] = e.property_flags;
				emitters.image_size[index] = properties.image[e.property_index]->image_size;
				emitters.image_frame_rate[index] = properties.image[e.property_index]->animation_frames > 1 && e.property_flags & tfxEmitterPropertyFlags_animate ? properties.frame_rate[e.property_index] : 0.f;
				emitters.image_frame_rate[index] = e.property_flags & tfxEmitterPropertyFlags_reverse_animation ? -emitters.image_frame_rate[index] : emitters.image_frame_rate[index];
				emitters.end_frame[index] = properties.end_frame[e.property_index];
				emitters.angle_offsets[index] = properties.angle_offsets[e.property_index];
				emitters.timeout[index] = 100.f;
				emitters.amount_remainder[index] = 0.f;
				emitters.qty_step_size[index] = 0.f;
				emitters.timeout_counter[index] = 0;
				emitters.emitter_size[index] = 0.f;
				emitters.hierarchy_depth[index] = hierarchy_depth;
				emitters.world_rotations[index] = 0.f;
				emitters.seed_index[index] = seed_index++;
				e.pm_index = index;		//Doesn't have much use beyond the editor?
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
				state_flags |= parent_state_flags & tfxEffectStateFlags_no_tween;
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
				state_flags |= properties.emission_type[e.property_index] == tfxLine && e.property_flags & tfxEmitterPropertyFlags_edge_traversal && (state_flags & tfxEmitterStateFlags_loop || state_flags & tfxEmitterStateFlags_kill) ? tfxEmitterStateFlags_is_line_loop_or_kill : 0;
				state_flags |= e.library->emitter_attributes[e.emitter_attributes].overtime.velocity_turbulance.GetMaxValue() > 0 ? tfxEmitterStateFlags_has_noise : 0;

				if (state_flags & tfxEmitterStateFlags_is_line_traversal) {
					emitters.property_flags[index] |= tfxEmitterPropertyFlags_relative_position;
				}

				if (is_sub_emitter) {
					state_flags |= tfxEmitterStateFlags_is_sub_emitter;
				}
				else {
					emitters.highest_particle_age[index] = tfxFRAME_LENGTH * 2.f;
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

		tfxU32 next_buffer = !current_ebuff;
		tfxU32 depth_starting_index[tfxLAYERS];
		current_sprite_buffer = flags & tfxEffectManagerFlags_double_buffer_sprites ? !current_sprite_buffer : 0;

		memset(sprite_index_point, 0, sizeof(tfxU32) * tfxLAYERS);

		for (tfxEachLayer) {
			depth_starting_index[layer] = depth_indexes[layer][current_depth_index_buffer].current_size;
		}

		for (tfxEachLayer) {
			ClearSoABuffer(&sprite_buffer[current_sprite_buffer][layer]);
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
				spawn_work_entry->properties = &library->emitter_properties;
				spawn_work_entry->sub_effects = &library->effect_infos[emitters.info_index[current_index]].sub_effectors;
				spawn_work_entry->amount_to_spawn = 0;
				spawn_work_entry->end_index = 0;
				spawn_work_entry->highest_particle_age = 0;

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
		random.Advance();

		for (auto &work_entry : spawn_work) {
			tfxU32 index = work_entry.emitter_index;
			emitters.highest_particle_age[index] = std::fmaxf(emitters.highest_particle_age[index], work_entry.highest_particle_age);
			effects.highest_particle_age[emitters.parent_index[index]] = emitters.highest_particle_age[index] + tfxFRAME_LENGTH;
		}
		spawn_work.free();

		if (!(flags & tfxEffectManagerFlags_unordered)) {
			for (tfxEachLayer) {
				if (depth_starting_index[layer] < depth_indexes[layer][current_depth_index_buffer].current_size) {
					tfxU32 next_depth_buffer = !current_depth_index_buffer;
					if (depth_indexes[layer][next_depth_buffer].capacity < depth_indexes[layer][current_depth_index_buffer].capacity) {
						depth_indexes[layer][next_depth_buffer].reserve(depth_indexes[layer][current_depth_index_buffer].capacity);
					}
					if (flags & tfxEffectManagerFlags_order_by_depth) {
						//No need to qsort ordered by age as the depth with all be 0 (depth is particle age)
						std::qsort(&depth_indexes[layer][current_depth_index_buffer][depth_starting_index[layer]], depth_indexes[layer][current_depth_index_buffer].current_size - depth_starting_index[layer], sizeof(tfxDepthIndex), SortDepth);
					}
					tfxU32 current_depth_index = 0;
					tfxU32 second_index = depth_starting_index[layer];
					for (auto &depth_index : depth_indexes[layer][current_depth_index_buffer]) {
						if (depth_starting_index[layer] != 0) {
							while (second_index < depth_indexes[layer][current_depth_index_buffer].current_size && depth_index.depth < depth_indexes[layer][current_depth_index_buffer][second_index].depth) {
								particle_arrays[ParticleBank(depth_indexes[layer][current_depth_index_buffer][second_index].particle_id)].depth_index[ParticleIndex(depth_indexes[layer][current_depth_index_buffer][second_index].particle_id)] = depth_indexes[layer][next_depth_buffer].current_size;
								depth_indexes[layer][next_depth_buffer].push_back(depth_indexes[layer][current_depth_index_buffer][second_index++]);
							}
						}
						particle_arrays[ParticleBank(depth_index.particle_id)].depth_index[ParticleIndex(depth_index.particle_id)] = depth_indexes[layer][next_depth_buffer].current_size;
						depth_indexes[layer][next_depth_buffer].push_back(depth_index);
						if (++current_depth_index == depth_starting_index[layer])
							break;
					}
					if (depth_starting_index[layer] != 0 && second_index < depth_indexes[layer][current_depth_index_buffer].current_size) {
						while (second_index < depth_indexes[layer][current_depth_index_buffer].current_size) {
							tfxU32 bank = ParticleBank(depth_indexes[layer][current_depth_index_buffer][second_index].particle_id);
							tfxU32 index = ParticleIndex(depth_indexes[layer][current_depth_index_buffer][second_index].particle_id);
							particle_arrays[ParticleBank(depth_indexes[layer][current_depth_index_buffer][second_index].particle_id)].depth_index[ParticleIndex(depth_indexes[layer][current_depth_index_buffer][second_index].particle_id)] = depth_indexes[layer][next_depth_buffer].current_size;
							depth_indexes[layer][next_depth_buffer].push_back(depth_indexes[layer][current_depth_index_buffer][second_index++]);
						}
					}
					assert(depth_indexes[layer][next_depth_buffer].current_size == depth_indexes[layer][current_depth_index_buffer].current_size);
					depth_indexes[layer][current_depth_index_buffer].clear();
					current_depth_index_buffer = next_depth_buffer;
				}
			}
		}

		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			tmpMTStack(tfxControlWorkEntry, work);
			for (int index : emitters_in_use[depth][next_buffer]) {
				tfxSoABuffer &bank = particle_array_buffers[emitters.particles_index[index]];
				int particles_to_update = bank.current_size;
				tfxU32 running_start_index = 0;
				while (particles_to_update > 0) {
					tfxControlWorkEntry &work_entry = work.next();
					work_entry.properties = &library->emitter_properties;
					work_entry.emitter_index = index;
					work_entry.start_index = running_start_index;
					work_entry.end_index = particles_to_update > mt_batch_size ? running_start_index + mt_batch_size : running_start_index + particles_to_update;
					tfxU32 circular_start = GetCircularIndex(&particle_array_buffers[emitters.particles_index[index]], work_entry.start_index);
					tfxU32 block_start_index = (circular_start / tfxDataWidth) * tfxDataWidth;
					work_entry.wide_end_index = (tfxU32)(ceilf((float)work_entry.end_index / tfxDataWidth)) * tfxDataWidth;
					work_entry.start_diff = circular_start - block_start_index;
					work_entry.wide_end_index = work_entry.wide_end_index - work_entry.start_diff < work_entry.end_index ? work_entry.wide_end_index + tfxDataWidth : work_entry.wide_end_index;
					particles_to_update -= mt_batch_size;
					running_start_index += mt_batch_size;
					if (flags & tfxEffectManagerFlags_3d_effects) {
						ControlParticles3d(*this, index, work_entry);
					}
					else {
						ControlParticles2d(*this, index, work_entry);
					}
				}
			}

			tfxCompleteAllWork(&work_queue);
			work.free();
			{
				tmpMTStack(tfxParticleAgeWorkEntry, work);
				for (int index : emitters_in_use[depth][next_buffer]) {
					tfxSoABuffer &bank = particle_array_buffers[emitters.particles_index[index]];
					tfxParticleAgeWorkEntry &work_entry = work.next();
					work_entry.properties = &library->emitter_properties;
					work_entry.start_index = bank.current_size - 1;
					work_entry.emitter_index = index;
					tfxU32 circular_start = GetCircularIndex(&particle_array_buffers[emitters.particles_index[index]], 0);
					tfxU32 block_start_index = (circular_start / tfxDataWidth) * tfxDataWidth;
					work_entry.wide_end_index = (tfxU32)(ceilf((float)bank.current_size / tfxDataWidth)) * tfxDataWidth;
					work_entry.start_diff = circular_start - block_start_index;
					work_entry.wide_end_index += work_entry.wide_end_index - work_entry.start_diff < bank.current_size ? tfxDataWidth : 0;
					work_entry.pm = this;
					if (!(flags & tfxEffectManagerFlags_single_threaded) && tfxNumberOfThreadsInAdditionToMain) {
						tfxAddWorkQueueEntry(&work_queue, &work_entry, ControlParticleAge);
					}
					else {
						ControlParticleAge(&work_queue, &work_entry);
					}
				}

				tfxCompleteAllWork(&work_queue);
				work.free();
			}
		}

		//Todo work queue this for each layer
		if (!(flags & tfxEffectManagerFlags_unordered)) {
			for (tfxEachLayer) {
				for (auto &depth_index : depth_indexes[layer][current_depth_index_buffer]) {
					if (depth_index.particle_id != tfxINVALID) {
						particle_arrays[ParticleBank(depth_index.particle_id)].depth_index[ParticleIndex(depth_index.particle_id)] = depth_indexes[layer][!current_depth_index_buffer].current_size;
						depth_indexes[layer][!current_depth_index_buffer].push_back(depth_index);
					}
				}
				//if(depth_indexes[layer][!current_depth_index_buffer].current_size)
					//assert(sprite_buffer[current_sprite_buffer][layer].current_size == depth_indexes[layer][!current_depth_index_buffer].current_size);
			}
		}

		tfxCompleteAllWork(&work_queue);

		for (tfxEachLayer) {
			depth_indexes[layer][current_depth_index_buffer].clear();
		}
		current_depth_index_buffer = !current_depth_index_buffer;

		if (flags & tfxEffectManagerFlags_order_by_depth && flags & tfxEffectManagerFlags_3d_effects) {
			if (flags & tfxEffectManagerFlags_guarantee_order) {
				for (tfxEachLayer) {
					tfxSortWorkEntry &work_entry = sorting_work_entry[layer];
					work_entry.bank = &particle_arrays;
					work_entry.depth_indexes = &depth_indexes[layer][current_depth_index_buffer];
					if (!(flags & tfxEffectManagerFlags_single_threaded) && tfxNumberOfThreadsInAdditionToMain > 0) {
						tfxAddWorkQueueEntry(&work_queue, &work_entry, InsertionSortDepth);
					}
					else {
						InsertionSortDepth(&work_queue, &work_entry);
					}
				}
			}
			else if (sort_passes > 0) {
				for (tfxEachLayer) {
					tfxvec<tfxDepthIndex> &depth_index = depth_indexes[layer][current_depth_index_buffer];
					//Add this to a work queue
					for (tfxU32 sorts = 0; sorts != sort_passes; ++sorts) {
						for (tfxU32 i = 1; i < depth_index.current_size; ++i) {
							float depth1 = depth_index[i - 1].depth;
							float depth2 = depth_index[i].depth;
							if (depth1 < depth2) {
								particle_arrays[ParticleBank(depth_index[i].particle_id)].depth_index[ParticleIndex(depth_index[i].particle_id)] = i - 1;
								particle_arrays[ParticleBank(depth_index[i - 1].particle_id)].depth_index[ParticleIndex(depth_index[i - 1].particle_id)] = i;
								std::swap(depth_index[i], depth_index[i - 1]);
							}
						}
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
				emitters.particles_index[current_index] = GrabParticleLists(*this, emitters.path_hash[current_index], 100);
				emitters_in_use[depth][next_buffer].push_back(current_index);
			}
		}

		current_ebuff = next_buffer;

		flags &= ~tfxEffectManagerFlags_update_base_values;

	}

	void ControlParticlePosition3d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		tfxParticleManager &pm = *work_entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[particles_index];

		const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];
		const tfxWideFloat overal_scale_wide = tfxWideSetSingle(pm.emitters.overal_scale[emitter_index]);

		tfxU32 running_sprite_index = work_entry->sprites_index;

		tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
		const tfxWideInt velocity_turbulance_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity_turbulance.lookup.last_frame);
		const tfxWideInt noise_resolution_last_frame = tfxWideSetSinglei(work_entry->graphs->noise_resolution.lookup.last_frame);

		//Noise
		const tfxWideInt velocity_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity.lookup.last_frame);
		const tfxWideInt spin_last_frame = tfxWideSetSinglei(work_entry->graphs->spin.lookup.last_frame);
		const tfxWideFloat velocity_adjuster = tfxWideSetSingle(pm.emitters.velocity_adjuster[emitter_index]);
		const tfxWideFloat angle_offsets_z = tfxWideSetSingle(pm.emitters.angle_offsets[emitter_index].roll);
		const tfxWideInt weight_last_frame = tfxWideSetSinglei(work_entry->graphs->weight.lookup.last_frame);

		for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
			tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;

			const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
			const tfxWideFloat age = tfxWideLoad(&bank.age[index]);
			_ReadBarrier();
			tfxWideFloat life = tfxWideDiv(age, max_age);
			life = tfxWideMul(life, max_life);
			life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

			const tfxWideFloat base_velocity = tfxWideLoad(&bank.base_velocity[index]);
			const tfxWideFloat base_spin = tfxWideLoad(&bank.base_spin[index]);
			const tfxWideFloat base_weight = tfxWideLoad(&bank.base_weight[index]);
			tfxWideFloat local_position_x = tfxWideLoad(&bank.position_x[index]);
			tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
			tfxWideFloat local_position_z = tfxWideLoad(&bank.position_z[index]);
			tfxWideArray roll;
			roll.m = tfxWideLoad(&bank.local_rotations_z[index]);
			tfxWideInt velocity_normal = tfxWideLoadi((tfxWideInt*)&bank.velocity_normal[index]);
			tfxWideFloat velocity_normal_x;
			tfxWideFloat velocity_normal_y;
			tfxWideFloat velocity_normal_z;
			UnPackWide10bit(velocity_normal, velocity_normal_x, velocity_normal_y, velocity_normal_z);

			tfxWideArray noise_x;
			tfxWideArray noise_y;
			tfxWideArray noise_z;

			if (emitter_flags & tfxEmitterStateFlags_has_noise) {

				float eps = 0.001f;
				float eps2 = 0.001f * 2.f;
				const tfxWideFloat noise_resolution = tfxWideLoad(&bank.noise_resolution[index]);

				tfxWideArrayi lookup_frame;
				lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_turbulance_last_frame);
				const tfxWideFloat lookup_velocity_turbulance = tfxWideLookupSet(work_entry->graphs->velocity_turbulance.lookup.values, lookup_frame);

				lookup_frame.m = tfxWideMini(tfxWideConverti(life), noise_resolution_last_frame);
				const tfxWideFloat lookup_noise_resolution = tfxWideMul(tfxWideLookupSet(work_entry->graphs->noise_resolution.lookup.values, lookup_frame), noise_resolution);

				const tfxWideFloat base_noise_offset = tfxWideLoad(&bank.noise_offset[index]);
				tfxWideFloat noise_offset = tfxWideMul(base_noise_offset, overal_scale_wide);

				tfxWideArray x, y, z;
				x.m = tfxWideAdd(tfxWideDiv(local_position_x, lookup_noise_resolution), noise_offset);
				y.m = tfxWideAdd(tfxWideDiv(local_position_y, lookup_noise_resolution), noise_offset);
				z.m = tfxWideAdd(tfxWideDiv(local_position_z, lookup_noise_resolution), noise_offset);

				for (int n = 0; n != tfxDataWidth; ++n) {
					tfx128 x4 = _mm_set1_ps(x.a[n]);
					tfx128 y4 = _mm_set1_ps(y.a[n]);
					tfx128 z4 = _mm_set1_ps(z.a[n]);

					tfx128 xeps4 = _mm_set_ps(x.a[n] - eps, x.a[n] + eps, x.a[n], x.a[n]);
					tfx128 xeps4r = _mm_set_ps(x.a[n], x.a[n], x.a[n] - eps, x.a[n] + eps);
					tfx128 yeps4 = _mm_set_ps(y.a[n], y.a[n], y.a[n] - eps, y.a[n] + eps);
					tfx128 zeps4 = _mm_set_ps(z.a[n] - eps, z.a[n] + eps, z.a[n], z.a[n]);
					tfx128 zeps4r = _mm_set_ps(z.a[n], z.a[n], z.a[n] - eps, z.a[n] + eps);

					//Find rate of change in YZ plane
					tfx128Array sample = tfxNoise4(x4, yeps4, zeps4);
					float a = (sample.a[0] - sample.a[1]) / eps2;
					//Average to find approximate derivative
					float b = (sample.a[2] - sample.a[3]) / eps2;
					noise_x.a[n] = a - b;

					y.a[n] += 100.f;
					tfx128 yeps4r = _mm_set_ps(y.a[n] - eps, y.a[n] + eps, y.a[n], y.a[n]);
					//Find rate of change in XZ plane
					sample = tfxNoise4(xeps4, y4, zeps4r);
					a = (sample.a[0] - sample.a[1]) / eps2;
					b = (sample.a[2] - sample.a[3]) / eps2;
					noise_y.a[n] = a - b;

					//Find rate of change in XY plane
					sample = tfxNoise4(xeps4r, yeps4r, z4);
					a = (sample.a[0] - sample.a[1]) / eps2;
					b = (sample.a[2] - sample.a[3]) / eps2;
					noise_z.a[n] = a - b;
				}

				noise_x.m = tfxWideMul(lookup_velocity_turbulance, noise_x.m);
				noise_y.m = tfxWideMul(lookup_velocity_turbulance, noise_y.m);
				noise_z.m = tfxWideMul(lookup_velocity_turbulance, noise_z.m);
			}

			tfxWideArrayi lookup_frame;
			lookup_frame.m = tfxWideMini(tfxWideConverti(life), spin_last_frame);
			const tfxWideFloat lookup_spin = tfxWideMul(tfxWideLookupSet(work_entry->graphs->spin.lookup.values, lookup_frame), base_spin);
			lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_last_frame);
			const tfxWideFloat lookup_velocity = tfxWideLookupSet(work_entry->graphs->velocity.lookup.values, lookup_frame);
			tfxWideArrayi lookup_frame_weight;
			lookup_frame_weight.m = tfxWideMini(tfxWideConverti(life), weight_last_frame);
			const tfxWideFloat lookup_weight = tfxWideLookupSet(work_entry->graphs->weight.lookup.values, lookup_frame_weight);

			//----Velocity Changes
			tfxWideFloat velocity_scalar = tfxWideMul(base_velocity, lookup_velocity);
			tfxWideFloat current_velocity_x = tfxWideMul(velocity_normal_x, velocity_scalar);
			tfxWideFloat current_velocity_y = tfxWideMul(velocity_normal_y, velocity_scalar);
			tfxWideFloat current_velocity_z = tfxWideMul(velocity_normal_z, velocity_scalar);
			if (emitter_flags & tfxEmitterStateFlags_has_noise) {
				current_velocity_x = tfxWideAdd(current_velocity_x, noise_x.m);
				current_velocity_y = tfxWideAdd(current_velocity_y, noise_y.m);
				current_velocity_z = tfxWideAdd(current_velocity_z, noise_z.m);
			}
			current_velocity_y = tfxWideSub(current_velocity_y, tfxWideMul(base_weight, tfxWideMul(lookup_weight, tfxUPDATE_TIME_WIDE)));
			current_velocity_x = tfxWideMul(tfxWideMul(current_velocity_x, tfxUPDATE_TIME_WIDE), velocity_adjuster);
			current_velocity_y = tfxWideMul(tfxWideMul(current_velocity_y, tfxUPDATE_TIME_WIDE), velocity_adjuster);
			current_velocity_z = tfxWideMul(tfxWideMul(current_velocity_z, tfxUPDATE_TIME_WIDE), velocity_adjuster);

			//----Spin and angle Changes
			if (emitter_flags & tfxEmitterStateFlags_can_spin) {
				roll.m = tfxWideAdd(roll.m, tfxWideMul(lookup_spin, tfxUPDATE_TIME_WIDE));
			}

			//----Position
			local_position_x = tfxWideAdd(local_position_x, tfxWideMul(current_velocity_x, overal_scale_wide));
			local_position_y = tfxWideAdd(local_position_y, tfxWideMul(current_velocity_y, overal_scale_wide));
			local_position_z = tfxWideAdd(local_position_z, tfxWideMul(current_velocity_z, overal_scale_wide));

			tfxWideStore(&bank.position_x[index], local_position_x);
			tfxWideStore(&bank.position_y[index], local_position_y);
			tfxWideStore(&bank.position_z[index], local_position_z);
			tfxWideStore(&bank.local_rotations_z[index], roll.m);
		}

		const tfxWideFloat emitter_size_y = tfxWideSetSingle(pm.emitters.emitter_size[emitter_index].y);
		const tfxWideInt emitter_flags_wide = tfxWideSetSinglei(emitter_flags);

		if (emitter_flags & tfxEmitterStateFlags_is_line_loop_or_kill) {
			//Todo: this should also update the captured position as well solving the issue of interpolating from the end back to the beginning
			if (emitter_flags & tfxEmitterStateFlags_kill) {
				for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
					tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;
					const tfxWideFloat offset_y = tfxWideMul(UnPackWide10bitY(tfxWideLoadi((tfxWideInt*)&bank.velocity_normal[index])), emitter_size_y);
					tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
					tfxWideInt flags = tfxWideLoadi((tfxWideInt*)&bank.flags[index]);

					//Lines - Reposition if the particle is travelling along a line
					tfxWideFloat length = tfxWideAbs(local_position_y);
					tfxWideInt remove_flags = tfxWideAndi(tfxWideSetSinglei(tfxParticleFlags_remove), tfxWideCasti(tfxWideGreater(length, emitter_size_y)));
					flags = tfxWideOri(flags, remove_flags);
					tfxWideStorei((tfxWideInt*)&bank.flags[index], flags);
				}
			}
			else {
				for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
					tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;
					const tfxWideFloat offset_y = tfxWideMul(UnPackWide10bitY(tfxWideLoadi((tfxWideInt*)&bank.velocity_normal[index])), emitter_size_y);
					tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
					tfxWideInt flags = tfxWideLoadi((tfxWideInt*)&bank.flags[index]);

					//Lines - Reposition if the particle is travelling along a line
					tfxWideFloat length = tfxWideAbs(local_position_y);
					tfxWideFloat at_end = tfxWideGreater(length, emitter_size_y);
					local_position_y = tfxWideSub(local_position_y, tfxWideAnd(at_end, offset_y));
					flags = tfxWideOri(flags, tfxWideAndi(tfxWideSetSinglei(tfxParticleFlags_capture_after_transform), tfxWideCasti(at_end)));
					tfxWideStorei((tfxWideInt*)&bank.flags[index], flags);
					tfxWideStore(&bank.position_y[index], local_position_y);
				}
			}
		}

		ControlParticleTransform3d(&pm.work_queue, data);
	}

	void ControlParticleTransform3d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		tfxParticleManager &pm = *work_entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[particles_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];

		tfxU32 running_sprite_index = work_entry->sprites_index;

		const tfxWideFloat e_world_position_x = tfxWideSetSingle(pm.emitters.world_position[emitter_index].x);
		const tfxWideFloat e_world_position_y = tfxWideSetSingle(pm.emitters.world_position[emitter_index].y);
		const tfxWideFloat e_world_position_z = tfxWideSetSingle(pm.emitters.world_position[emitter_index].z);
		const tfxWideFloat e_world_rotations_x = tfxWideSetSingle(pm.emitters.world_rotations[emitter_index].x);
		const tfxWideFloat e_world_rotations_y = tfxWideSetSingle(pm.emitters.world_rotations[emitter_index].y);
		const tfxWideFloat e_world_rotations_z = tfxWideSetSingle(pm.emitters.world_rotations[emitter_index].z);
		const tfxWideFloat e_handle_x = tfxWideSetSingle(pm.emitters.handle[emitter_index].x);
		const tfxWideFloat e_handle_y = tfxWideSetSingle(pm.emitters.handle[emitter_index].y);
		const tfxWideFloat e_handle_z = tfxWideSetSingle(pm.emitters.handle[emitter_index].z);
		const tfxWideFloat e_scale_x = tfxWideSetSingle(pm.emitters.scale[emitter_index].x);
		const tfxWideFloat e_scale_y = tfxWideSetSingle(pm.emitters.scale[emitter_index].y);
		const tfxWideFloat e_scale_z = tfxWideSetSingle(pm.emitters.scale[emitter_index].z);
		tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
		const tfxWideInt stretch_last_frame = tfxWideSetSinglei(work_entry->graphs->stretch.lookup.last_frame);
		const tfxWideFloat stretch = tfxWideSetSingle(pm.emitters.stretch[emitter_index]);

		const tfxWideInt capture_after_transform = tfxWideSetSinglei(tfxParticleFlags_capture_after_transform);
		const tfxWideInt xor_capture_after_transform_flag = tfxWideXOri(tfxWideSetSinglei(tfxParticleFlags_capture_after_transform), tfxWideSetSinglei(-1));
		tfxMatrix4 &e_matrix = pm.emitters.matrix[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVectorAlignType vector_align_type = work_entry->properties->vector_align_type[property_index];
		const tfxBillboardingOptions billboard_option = work_entry->properties->billboard_option[property_index];
		const tfxEmissionType emission_type = work_entry->properties->emission_type[property_index];
		const tfxU32 sprite_layer = work_entry->properties->layer[property_index];
		tfxSpriteSoA &sprites = *work_entry->sprites;
		tfxU32 start_diff = work_entry->start_diff;
		tfxWideArray p_stretch;

		for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
			tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;
			const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
			const tfxWideFloat age = tfxWideLoad(&bank.age[index]);
			_ReadBarrier();
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
			tfxWideArray rotations_x;
			tfxWideArray rotations_y;
			tfxWideArray rotations_z;
			rotations_x.m = tfxWideLoad(&bank.local_rotations_x[index]);
			rotations_y.m = tfxWideLoad(&bank.local_rotations_y[index]);
			rotations_z.m = tfxWideLoad(&bank.local_rotations_z[index]);
			const tfxWideInt velocity_normal = tfxWideLoadi((tfxWideInt*)&bank.velocity_normal[index]);
			tfxWideFloat velocity_normal_x;
			tfxWideFloat velocity_normal_y;
			tfxWideFloat velocity_normal_z;
			UnPackWide10bit(velocity_normal, velocity_normal_x, velocity_normal_y, velocity_normal_z);
			tfxWideArray captured_position_x;
			tfxWideArray captured_position_y;
			tfxWideArray captured_position_z;
			captured_position_x.m = tfxWideLoad(&bank.captured_position_x[index]);
			captured_position_y.m = tfxWideLoad(&bank.captured_position_y[index]);
			captured_position_z.m = tfxWideLoad(&bank.captured_position_z[index]);
			tfxWideInt flags = tfxWideLoadi((tfxWideInt*)&bank.flags[index]);
			tfxWideFloat capture_flag = tfxWideCast(tfxWideGreateri(tfxWideAndi(flags, capture_after_transform), tfxWideSetZeroi()));
			tfxWideFloat xor_capture_flag = tfxWideEquals(capture_flag, tfxWideSetZero());
			_ReadBarrier();

			if (property_flags & tfxEmitterPropertyFlags_relative_position || (property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == tfxLine)) {
				position_x.m = tfxWideAdd(position_x.m, e_handle_x);
				position_y.m = tfxWideAdd(position_y.m, e_handle_y);
				position_z.m = tfxWideAdd(position_z.m, e_handle_z);
				mmWideTransformVector(e_matrix, position_x.m, position_y.m, position_z.m);
				position_x.m = tfxWideAdd(tfxWideMul(position_x.m, e_scale_x), e_world_position_x);
				position_y.m = tfxWideAdd(tfxWideMul(position_y.m, e_scale_y), e_world_position_y);
				position_z.m = tfxWideAdd(tfxWideMul(position_z.m, e_scale_z), e_world_position_z);
			}
			else if (property_flags & tfxEmitterPropertyFlags_relative_angle) {
				rotations_x.m = tfxWideAdd(rotations_x.m, e_world_rotations_x);
				rotations_y.m = tfxWideAdd(rotations_y.m, e_world_rotations_y);
				rotations_z.m = tfxWideAdd(rotations_z.m, e_world_rotations_z);
			}

			captured_position_x.m = tfxWideAdd(tfxWideAnd(position_x.m, capture_flag), tfxWideAnd(captured_position_x.m, xor_capture_flag));
			captured_position_y.m = tfxWideAdd(tfxWideAnd(position_y.m, capture_flag), tfxWideAnd(captured_position_y.m, xor_capture_flag));
			captured_position_z.m = tfxWideAdd(tfxWideAnd(position_z.m, capture_flag), tfxWideAnd(captured_position_z.m, xor_capture_flag));

			flags = tfxWideAndi(flags, xor_capture_after_transform_flag);

			tfxWideFloat alignment_vector_x;
			tfxWideFloat alignment_vector_y;
			tfxWideFloat alignment_vector_z;
			if (vector_align_type == tfxVectorAlignType_motion) {
				alignment_vector_x = tfxWideSub(position_x.m, captured_position_x.m);
				alignment_vector_y = tfxWideAdd(tfxWideSub(position_y.m, captured_position_y.m), tfxWideSetSingle(0.000001f)); //epsilon to prevent divide by 0
				alignment_vector_z = tfxWideSub(position_z.m, captured_position_z.m);
				tfxWideFloat l = tfxWideMul(alignment_vector_x, alignment_vector_x);
				l = tfxWideAdd(l, tfxWideMul(alignment_vector_y, alignment_vector_y));
				l = tfxWideAdd(l, tfxWideMul(alignment_vector_z, alignment_vector_z));
				l = tfxWideSqrt(l);
				p_stretch.m = tfxWideMul(p_stretch.m, tfxWideMul(l, tfxWideSetSingle(10.f)));	//This is too arbitrary, think up a better solution!
				alignment_vector_x = tfxWideDiv(alignment_vector_x, l);
				alignment_vector_y = tfxWideDiv(alignment_vector_y, l);
				alignment_vector_z = tfxWideDiv(alignment_vector_z, l);
			}
			else if (vector_align_type == tfxVectorAlignType_emission && property_flags & tfxEmitterPropertyFlags_relative_position) {
				alignment_vector_x = velocity_normal_x;
				alignment_vector_y = velocity_normal_y;
				alignment_vector_z = velocity_normal_z;
				mmWideTransformVector(e_matrix, alignment_vector_x, alignment_vector_y, alignment_vector_z);
			}
			else if (vector_align_type == tfxVectorAlignType_emission) {
				alignment_vector_x = velocity_normal_x;
				alignment_vector_y = velocity_normal_y;
				alignment_vector_z = velocity_normal_z;
			}
			else if (vector_align_type == tfxVectorAlignType_emitter) {
				alignment_vector_x = tfxWideSetSingle(0.f);
				alignment_vector_y = tfxWideSetSingle(1.f);
				alignment_vector_z = tfxWideSetSingle(0.f);
				mmWideTransformVector(e_matrix, alignment_vector_x, alignment_vector_y, alignment_vector_z);
			}

			//sprites.transform_3d.captured_position = captured_position;
			//alignment_vector_y.m = tfxWideAdd(alignment_vector_y.m, tfxWideSetSingle(0.002f));	//We don't want a 0 alignment normal
			tfxWideArrayi packed;
			packed.m = PackWide10bit(alignment_vector_x, alignment_vector_y, alignment_vector_z, billboard_option & 0x00000003);

			tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
			if (!(pm.flags & tfxEffectManagerFlags_unordered)) {	//Predictable
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 sprite_depth_index = bank.depth_index[index + j];
					sprites.stretch[sprite_depth_index] = p_stretch.a[j];
					sprites.transform_3d[sprite_depth_index].rotations.x = rotations_x.a[j];
					sprites.transform_3d[sprite_depth_index].rotations.y = rotations_y.a[j];
					sprites.transform_3d[sprite_depth_index].rotations.z = rotations_z.a[j];
					sprites.transform_3d[sprite_depth_index].position.x = position_x.a[j];
					sprites.transform_3d[sprite_depth_index].position.y = position_y.a[j];
					sprites.transform_3d[sprite_depth_index].position.z = position_z.a[j];
					sprites.alignment[sprite_depth_index] = packed.a[j];
					bank.captured_position_x[index + j] = sprites.transform_3d[sprite_depth_index].position.x;
					bank.captured_position_y[index + j] = sprites.transform_3d[sprite_depth_index].position.y;
					bank.captured_position_z[index + j] = sprites.transform_3d[sprite_depth_index].position.z;
					pm.depth_indexes[sprite_layer][pm.current_depth_index_buffer][sprite_depth_index].depth = LengthVec3NoSqR(sprites.transform_3d[sprite_depth_index].position - pm.camera_position);
					running_sprite_index++;
				}
			}
			else {
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					sprites.stretch[running_sprite_index] = p_stretch.a[j];
					sprites.transform_3d[running_sprite_index].rotations.x = rotations_x.a[j];
					sprites.transform_3d[running_sprite_index].rotations.y = rotations_y.a[j];
					sprites.transform_3d[running_sprite_index].rotations.z = rotations_z.a[j];
					sprites.transform_3d[running_sprite_index].position.x = position_x.a[j];
					sprites.transform_3d[running_sprite_index].position.y = position_y.a[j];
					sprites.transform_3d[running_sprite_index].position.z = position_z.a[j];
					sprites.alignment[running_sprite_index] = packed.a[j];
					bank.captured_position_x[index + j] = sprites.transform_3d[running_sprite_index].position.x;
					bank.captured_position_y[index + j] = sprites.transform_3d[running_sprite_index].position.y;
					bank.captured_position_z[index + j] = sprites.transform_3d[running_sprite_index].position.z;
					running_sprite_index++;
				}
			}
			tfxWideStorei((tfxWideInt*)&bank.flags[index], flags);
			start_diff = 0;
		}
	}

	void ControlParticleTransform2d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		tfxParticleManager &pm = *work_entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[particles_index];

		tfxSpriteSoA &sprites = *work_entry->sprites;
		tfxU32 running_sprite_index = work_entry->sprites_index;

		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxWideFloat e_world_position_x = tfxWideSetSingle(pm.emitters.world_position[emitter_index].x);
		const tfxWideFloat e_world_position_y = tfxWideSetSingle(pm.emitters.world_position[emitter_index].y);
		const tfxWideFloat e_world_rotations_roll = tfxWideSetSingle(pm.emitters.world_rotations[emitter_index].roll);
		const tfxWideFloat e_handle_x = tfxWideSetSingle(pm.emitters.handle[emitter_index].x);
		const tfxWideFloat e_handle_y = tfxWideSetSingle(pm.emitters.handle[emitter_index].y);
		const tfxWideFloat e_scale_x = tfxWideSetSingle(pm.emitters.scale[emitter_index].x);
		const tfxWideFloat e_scale_y = tfxWideSetSingle(pm.emitters.scale[emitter_index].y);
		const tfxWideInt capture_after_transform = tfxWideSetSinglei(tfxParticleFlags_capture_after_transform);
		const tfxWideInt xor_capture_after_transform_flag = tfxWideXOri(tfxWideSetSinglei(tfxParticleFlags_capture_after_transform), tfxWideSetSinglei(-1));
		tfxMatrix4 &e_matrix = pm.emitters.matrix[emitter_index];
		tfxVec3 &e_scale = pm.emitters.scale[emitter_index];
		tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);

		tfxU32 start_diff = work_entry->start_diff;

		for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
			tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;

			tfxWideArray position_x;
			tfxWideArray position_y;
			tfxWideArray captured_position_x;
			tfxWideArray captured_position_y;
			tfxWideArray roll;
			position_x.m = tfxWideLoad(&bank.position_x[index]);
			position_y.m = tfxWideLoad(&bank.position_y[index]);
			roll.m = tfxWideLoad(&bank.local_rotations_z[index]);
			captured_position_x.m = tfxWideLoad(&bank.captured_position_x[index]);
			captured_position_y.m = tfxWideLoad(&bank.captured_position_y[index]);
			tfxWideInt flags = tfxWideLoadi((tfxWideInt*)&bank.flags[index]);
			tfxWideFloat capture_flag = tfxWideCast(tfxWideGreateri(tfxWideAndi(flags, capture_after_transform), tfxWideSetZeroi()));
			tfxWideFloat xor_capture_flag = tfxWideEquals(capture_flag, tfxWideSetZero());

			const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
			const tfxWideFloat age = tfxWideLoad(&bank.age[index]);
			_ReadBarrier();
			tfxWideFloat life = tfxWideDiv(age, max_age);
			life = tfxWideMul(life, max_life);
			life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

			if (property_flags & tfxEmitterPropertyFlags_relative_position) {
				position_x.m = tfxWideAdd(position_x.m, e_handle_x);
				position_y.m = tfxWideAdd(position_y.m, e_handle_y);
				mmWideTransformVector(e_matrix, position_x.m, position_y.m);
				position_x.m = tfxWideAdd(tfxWideMul(position_x.m, e_scale_x), e_world_position_x);
				position_y.m = tfxWideAdd(tfxWideMul(position_y.m, e_scale_y), e_world_position_y);
			}
			if (property_flags & tfxEmitterPropertyFlags_relative_angle) {
				roll.m = tfxWideAdd(roll.m, e_world_rotations_roll);
			}

			flags = tfxWideAndi(flags, xor_capture_after_transform_flag);

			captured_position_x.m = tfxWideAdd(tfxWideAnd(position_x.m, capture_flag), tfxWideAnd(captured_position_x.m, xor_capture_flag));
			captured_position_y.m = tfxWideAdd(tfxWideAnd(position_y.m, capture_flag), tfxWideAnd(captured_position_y.m, xor_capture_flag));

			tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
			if (!(pm.flags & tfxEffectManagerFlags_unordered)) {	//Predictable
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 sprite_depth_index = bank.depth_index[index + j];
					sprites.transform_2d[sprite_depth_index].rotation = roll.a[j];
					sprites.transform_2d[sprite_depth_index].position.x = position_x.a[j];
					sprites.transform_2d[sprite_depth_index].position.y = position_y.a[j];
					bank.captured_position_x[index + j] = sprites.transform_2d[sprite_depth_index].position.x;
					bank.captured_position_y[index + j] = sprites.transform_2d[sprite_depth_index].position.y;
					running_sprite_index++;
				}
			}
			else {
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					sprites.transform_2d[running_sprite_index].rotation = roll.a[j];
					sprites.transform_2d[running_sprite_index].position.x = position_x.a[j];
					sprites.transform_2d[running_sprite_index].position.y = position_y.a[j];
					bank.captured_position_x[index + j] = sprites.transform_2d[running_sprite_index].position.x;
					bank.captured_position_y[index + j] = sprites.transform_2d[running_sprite_index].position.y;
					running_sprite_index++;
				}
			}
			tfxWideStorei((tfxWideInt*)&bank.flags[index], flags);
			start_diff = 0;
		}
	}

	void ControlParticleSize(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		const tfxU32 particles_index = work_entry->pm->emitters.particles_index[emitter_index];
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = pm.particle_arrays[particles_index];

		const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];
		const tfxWideInt width_last_frame = tfxWideSetSinglei(work_entry->graphs->width.lookup.last_frame);
		const tfxWideInt height_last_frame = tfxWideSetSinglei(work_entry->graphs->height.lookup.last_frame);
		const tfxWideFloat overal_scale = tfxWideSetSingle(pm.emitters.overal_scale[emitter_index]);

		tfxU32 running_sprite_index = work_entry->sprites_index;

		tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);

		tfxU32 start_diff = work_entry->start_diff;

		tfxWideArrayi lookup_frame;
		tfxSpriteSoA &sprites = *work_entry->sprites;
		tfxWideArray scale_x;
		tfxWideArray scale_y;

		for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
			tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;

			const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
			const tfxWideFloat age = tfxWideLoad(&bank.age[index]);
			_ReadBarrier();
			tfxWideFloat life = tfxWideDiv(age, max_age);
			life = tfxWideMul(life, max_life);
			life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

			lookup_frame.m = tfxWideMini(tfxWideConverti(life), width_last_frame);
			const tfxWideFloat lookup_width = tfxWideLookupSet(work_entry->graphs->width.lookup.values, lookup_frame);

			lookup_frame.m = tfxWideMini(tfxWideConverti(life), height_last_frame);
			const tfxWideFloat lookup_height = tfxWideLookupSet(work_entry->graphs->height.lookup.values, lookup_frame);

			const tfxWideFloat base_size_x = tfxWideLoad(&bank.base_size_x[index]);
			const tfxWideFloat base_size_y = tfxWideLoad(&bank.base_size_y[index]);

			//----Size Changes
			scale_x.m = tfxWideMul(base_size_x, lookup_width);

			if (emitter_flags & tfxEmitterStateFlags_lifetime_uniform_size) {
				scale_y.m = tfxWideMul(lookup_width, base_size_y);
				if (emitter_flags & tfxEmitterPropertyFlags_base_uniform_size)
					scale_y.m = tfxWideMin(scale_x.m, scale_y.m);
			}
			else
				scale_y.m = tfxWideMul(lookup_height, base_size_y);

			scale_x.m = tfxWideMul(scale_x.m, overal_scale);
			scale_y.m = tfxWideMul(scale_y.m, overal_scale);

			tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
			if (pm.flags & tfxEffectManagerFlags_3d_effects) { //Predictable
				if (!(pm.flags & tfxEffectManagerFlags_unordered)) {	//Predictable
					for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
						tfxU32 sprite_depth_index = bank.depth_index[index + j];
						sprites.transform_3d[sprite_depth_index].scale.x = scale_x.a[j];
						sprites.transform_3d[sprite_depth_index].scale.y = scale_y.a[j];
						running_sprite_index++;
					}
				}
				else {
					for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
						sprites.transform_3d[running_sprite_index].scale.x = scale_x.a[j];
						sprites.transform_3d[running_sprite_index++].scale.y = scale_y.a[j];
					}
				}
			}
			else {
				if (!(pm.flags & tfxEffectManagerFlags_unordered)) {	//Predictable
					for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
						tfxU32 sprite_depth_index = bank.depth_index[index + j];
						sprites.transform_2d[sprite_depth_index].scale.x = scale_x.a[j];
						sprites.transform_2d[sprite_depth_index].scale.y = scale_y.a[j];
						running_sprite_index++;
					}
				}
				else {
					for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
						sprites.transform_2d[running_sprite_index].scale.x = scale_x.a[j];
						sprites.transform_2d[running_sprite_index++].scale.y = scale_y.a[j];
					}
				}
			}
			start_diff = 0;
		}
	}

	void ControlParticleColor(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		const tfxU32 particles_index = work_entry->pm->emitters.particles_index[emitter_index];
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[particles_index];

		const tfxWideFloat global_intensity = tfxWideSetSingle(pm.emitters.intensity[emitter_index]);
		const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];

		tfxU32 running_sprite_index = work_entry->sprites_index;

		tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
		tfxU32 start_diff = work_entry->start_diff;

		const tfxWideInt last_frame_color = tfxWideSetSinglei(work_entry->graphs->red.lookup.last_frame);
		const tfxWideInt last_frame_intensity = tfxWideSetSinglei(work_entry->graphs->intensity.lookup.last_frame);
		const tfxWideInt last_frame_opacity = tfxWideSetSinglei(work_entry->graphs->blendfactor.lookup.last_frame);
		tfxWideArray wide_alpha;
		tfxWideArray wide_intensity;
		tfxWideArrayi lookup_frame;
		tfxWideArrayi packed_color;
		tfxSpriteSoA &sprites = *work_entry->sprites;

		for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
			tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;

			const tfxWideFloat age = tfxWideLoad(&bank.age[index]);
			const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
			_ReadBarrier();
			tfxWideFloat life = tfxWideDiv(age, max_age);
			life = tfxWideMul(life, max_life);
			life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

			lookup_frame.m = tfxWideMini(tfxWideConverti(life), last_frame_color);
			const tfxWideFloat lookup_red = tfxWideLookupSet(work_entry->graphs->red.lookup.values, lookup_frame);

			lookup_frame.m = tfxWideMini(tfxWideConverti(life), last_frame_color);
			const tfxWideFloat lookup_green = tfxWideLookupSet(work_entry->graphs->green.lookup.values, lookup_frame);

			lookup_frame.m = tfxWideMini(tfxWideConverti(life), last_frame_color);
			const tfxWideFloat lookup_blue = tfxWideLookupSet(work_entry->graphs->blue.lookup.values, lookup_frame);

			lookup_frame.m = tfxWideMini(tfxWideConverti(life), last_frame_intensity);
			const tfxWideFloat lookup_intensity = tfxWideLookupSet(work_entry->graphs->intensity.lookup.values, lookup_frame);

			lookup_frame.m = tfxWideMini(tfxWideConverti(life), last_frame_opacity);
			const tfxWideFloat lookup_opacity = tfxWideLookupSet(work_entry->graphs->blendfactor.lookup.values, lookup_frame);

			//----Color changes
			wide_alpha.m = tfxWideMul(tfxWIDE255, lookup_opacity);
			wide_intensity.m = tfxWideMul(global_intensity, lookup_intensity);

			if (!(emitter_flags & tfxEmitterStateFlags_random_color)) {
				packed_color.m = PackWideColor( tfxWideMul(tfxWIDE255, lookup_red), tfxWideMul(tfxWIDE255, lookup_green), tfxWideMul(tfxWIDE255, lookup_blue), wide_alpha.m);
			}
			else {
				packed_color.m = tfxWideLoadi((tfxWideInt*)&bank.color[index]);
			}

			tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
			if (!(pm.flags & tfxEffectManagerFlags_unordered)) {	//Predictable
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 sprite_depth_index = bank.depth_index[index + j];
					sprites.color[sprite_depth_index].color = packed_color.a[j];
					sprites.intensity[sprite_depth_index] = wide_intensity.a[j];
					running_sprite_index++;
				}
			}
			else {
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					sprites.color[running_sprite_index].color = packed_color.a[j];
					sprites.intensity[running_sprite_index++] = wide_intensity.a[j];
				}
			}
			start_diff = 0;
		}

	}

	void ControlParticleImageFrame(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		const tfxU32 particles_index = work_entry->pm->emitters.particles_index[emitter_index];
		const tfxU32 property_index = work_entry->pm->emitters.properties_index[emitter_index];
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = pm.particle_arrays[particles_index];
		tfxImageData *image = work_entry->properties->image[property_index];
		const tfxBillboardingOptions billboard_option = work_entry->properties->billboard_option[property_index];

		tfxU32 start_diff = work_entry->start_diff;
		float image_frames[tfxDataWidth];

		tfxWideFloat image_frame_rate = tfxWideSetSingle(pm.emitters.image_frame_rate[emitter_index]);
		image_frame_rate = tfxWideMul(image_frame_rate, tfxUPDATE_TIME_WIDE);
		tfxWideFloat end_frame = tfxWideSetSingle(pm.emitters.end_frame[emitter_index]);
		tfxWideFloat frames = tfxWideSetSingle(pm.emitters.end_frame[emitter_index] + 1);
		tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];
		tfxEmitterStateFlags property_flags = pm.emitters.property_flags[emitter_index];

		tfxU32 running_sprite_index = work_entry->sprites_index;
		tfxSpriteSoA &sprites = *work_entry->sprites;

		for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
			tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;

			tfxWideFloat image_frame = tfxWideLoad(&bank.image_frame[index]);

			//----Image animation
			image_frame = tfxWideAdd(image_frame, image_frame_rate);
			if (emitter_flags & tfxEmitterStateFlags_play_once) {
				image_frame = tfxWideMin(image_frame, end_frame);
				image_frame = tfxWideMax(image_frame, tfxWideSetZero());
			}
			else if(property_flags & tfxEmitterPropertyFlags_reverse_animation) {
				tfxWideFloat mask = tfxWideLess(image_frame, tfxWideSetZero());
				image_frame = tfxWideAdd(image_frame, tfxWideAnd(mask, frames));
			}
			else {
				tfxWideFloat mask = tfxWideGreaterEqual(image_frame, frames);
				image_frame = tfxWideSub(image_frame, tfxWideAnd(mask, frames));
			}

			tfxWideStore(&bank.image_frame[index], image_frame);
			tfxWideStore(image_frames, image_frame);
			tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
			if (!(pm.flags & tfxEffectManagerFlags_unordered)) {	//Predictable
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 sprite_depth_index = bank.depth_index[index + j];
					tfxU32 &sprites_index = bank.sprite_index[index + j];
					float &age = bank.age[index + j];
					sprites.captured_index[sprite_depth_index] = age == 0.f ? (pm.current_sprite_buffer << 28) + sprite_depth_index : (!pm.current_sprite_buffer << 28) + (sprites_index & 0x0FFFFFFF);
					sprites_index = (work_entry->layer << 28) + sprite_depth_index;
					sprites.image_frame_plus[sprite_depth_index] = (billboard_option << 24) + ((tfxU32)image_frames[j] << 16) + (property_index);
					tfxU32 ifp = sprites.image_frame_plus[sprite_depth_index];
					running_sprite_index++;
				}
			}
			else {
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 &sprites_index = bank.sprite_index[index + j];
					float &age = bank.age[index + j];
					sprites.captured_index[running_sprite_index] = age == 0.f ? (pm.current_sprite_buffer << 28) + running_sprite_index : (!pm.current_sprite_buffer << 28) + (sprites_index & 0x0FFFFFFF);
					sprites_index = (work_entry->layer << 28) + running_sprite_index;
					sprites.image_frame_plus[running_sprite_index++] = (billboard_option << 24) + ((tfxU32)image_frames[j] << 16) + (property_index);
				}
			}
			start_diff = 0;
		}

	}

	void tfxParticleManager::UpdateParticleOrderOnly() {
		//todo:
		return;
	}

	tfxvec<tfxU32> *tfxParticleManager::GetEffectBuffer(tfxU32 depth) {
		return &effects_in_use[depth][current_ebuff];
	}

	tfxvec<tfxU32> *tfxParticleManager::GetEmitterBuffer(tfxU32 depth) {
		return &emitters_in_use[depth][current_ebuff];
	}

	void tfxParticleManager::InitFor2d(tfxLibrary *lib, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfxParticleManagerModes mode, bool double_buffer_sprites, bool dynamic_sprite_allocation, tfxU32 multi_threaded_batch_size) {
		assert(mode == tfxParticleManagerMode_unordered || mode == tfxParticleManagerMode_ordered_by_age);	//Only these 2 modes are available for 2d effects
		max_effects = effects_limit;
		mt_batch_size = multi_threaded_batch_size;
		tfxInitialiseWorkQueue(&work_queue);
		library = lib;

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

			InitSprite2dSoA(&sprite_buffer[0][layer], &sprites[0][layer], tfxMax((layer_max_values[layer] / tfxDataWidth + 1) * tfxDataWidth, 8));
			if (flags & tfxEffectManagerFlags_double_buffer_sprites) {
				InitSprite2dSoA(&sprite_buffer[1][layer], &sprites[1][layer], tfxMax((layer_max_values[layer] / tfxDataWidth + 1) * tfxDataWidth, 8));
			}
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

		flags |= dynamic_sprite_allocation ? tfxEffectManagerFlags_dynamic_sprite_allocation : 0;

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

	void tfxParticleManager::InitFor3d(tfxLibrary *lib, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfxParticleManagerModes mode, bool double_buffer_sprites, bool dynamic_sprite_allocation, tfxU32 multi_threaded_batch_size) {
		max_effects = effects_limit;
		mt_batch_size = multi_threaded_batch_size;
		tfxInitialiseWorkQueue(&work_queue);
		library = lib;

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
		if (double_buffer_sprites)
			flags |= tfxEffectManagerFlags_double_buffer_sprites;

		for (tfxEachLayer) {
			max_cpu_particles_per_layer[layer] = layer_max_values[layer];

#ifdef tfxTRACK_MEMORY
			memcpy(sprites3d[layer].name, "ParticleManager::sprites3d\0", 27);
#endif

			InitSprite3dSoA(&sprite_buffer[0][layer], &sprites[0][layer], tfxMax((layer_max_values[layer] / tfxDataWidth + 1) * tfxDataWidth, 8));
			if (flags & tfxEffectManagerFlags_double_buffer_sprites) {
				InitSprite3dSoA(&sprite_buffer[1][layer], &sprites[1][layer], tfxMax((layer_max_values[layer] / tfxDataWidth + 1) * tfxDataWidth, 8));
			}

		}

		if (flags & tfxEffectManagerFlags_ordered_by_age || flags & tfxEffectManagerFlags_order_by_depth) {
			FreeParticleBanks();
			for (tfxEachLayer) {
				tfxParticleSoA lists;
				tfxU32 index = particle_arrays.locked_push_back(lists);
				tfxSoABuffer buffer;
				particle_array_buffers.push_back(buffer);
				assert(index == particle_array_buffers.current_size - 1);
				InitParticleSoA(&particle_array_buffers[index], &particle_arrays.back(), tfxMax(max_cpu_particles_per_layer[layer], 8));
				particle_array_buffers[index].user_data = &particle_arrays.back();
				tfxResizeParticleSoACallback(&particle_array_buffers[index], 0);
			}
		}
		/*
		else if (flags & tfxEffectManagerFlags_order_by_depth) {
			FreeParticleBanks();
			for (tfxEachLayerDB) {
				tfxParticleSoA lists;
				tfxU32 index = particle_arrays.locked_push_back(lists);
				tfxSoABuffer buffer;
				particle_array_buffers.push_back(buffer);
				assert(index == particle_array_buffers.current_size - 1);
				InitParticleSoA(&particle_array_buffers[index], &particle_arrays.back(), tfxMax(max_cpu_particles_per_layer[layer / 2], 8));
				particle_array_buffers[index].user_data = &particle_arrays.back();
				tfxResizeParticleSoACallback(&particle_array_buffers[index], 0);
			}
		}
		*/

		flags |= dynamic_sprite_allocation ? tfxEffectManagerFlags_dynamic_sprite_allocation : 0;

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

	void tfxParticleManager::InitFor2d(tfxLibrary *lib, unsigned int effects_limit, tfxParticleManagerModes mode) {
		tfxU32 layer_max_values[tfxLAYERS];
		memset(layer_max_values, 0, 16);
		InitFor2d(lib, layer_max_values, effects_limit, mode);
		flags |= tfxEffectManagerFlags_dynamic_sprite_allocation;
	}

	void tfxParticleManager::InitFor3d(tfxLibrary *lib, unsigned int effects_limit, tfxParticleManagerModes mode) {
		tfxU32 layer_max_values[tfxLAYERS];
		memset(layer_max_values, 0, 16);
		InitFor3d(lib, layer_max_values, effects_limit, mode);
		flags |= tfxEffectManagerFlags_dynamic_sprite_allocation;
	}

	void tfxParticleManager::CreateParticleBanksForEachLayer() {
		FreeParticleBanks();
		for (tfxEachLayer) {
			tfxParticleSoA lists;
			tfxU32 index = particle_arrays.locked_push_back(lists);
			tfxSoABuffer buffer;
			buffer.resize_callback = tfxResizeParticleSoACallback;
			buffer.user_data = &particle_arrays.back();
			particle_array_buffers.push_back(buffer);
			assert(index == particle_array_buffers.current_size - 1);
			InitParticleSoA(&particle_array_buffers[index], &particle_arrays.back(), tfxMax(max_cpu_particles_per_layer[layer], 8));
			particle_array_buffers[index].user_data = &particle_arrays.back();
			tfxResizeParticleSoACallback(&particle_array_buffers[index], 0);
		}
	}

	void tfxParticleManager::Reconfigure(tfxParticleManagerModes mode, tfxU32 req_sort_passes, bool is_3d) {
		FreeParticleBanks();
		for (auto &bank : free_particle_lists.data) {
			bank.free_all();
		}
		free_particle_lists.FreeAll();

		tfxParticleManagerFlags current_flags = flags & tfxEffectManagerFlags_dynamic_sprite_allocation | flags & tfxEffectManagerFlags_double_buffer_sprites;

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
		else
			flags &= ~tfxEffectManagerFlags_3d_effects;


		flags |= current_flags;

		for (tfxEachLayer) {
			ClearSoABuffer(&sprite_buffer[0][layer]);
			if (flags & tfxEffectManagerFlags_double_buffer_sprites) {
				ClearSoABuffer(&sprite_buffer[1][layer]);
			}
			depth_indexes[layer][0].clear();
			depth_indexes[layer][1].clear();
		}

		memset(sprite_index_point, 0, 4 * tfxLAYERS);

		ClearSoABuffer(&emitter_buffers);
		ClearSoABuffer(&effect_buffers);

		sort_passes = req_sort_passes;
	}

	void tfxParticleManager::InitForBoth(tfxLibrary *lib, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfxParticleManagerModes mode, bool double_buffer_sprites, bool dynamic_sprite_allocation, tfxU32 multi_threaded_batch_size) {
		max_effects = effects_limit;
		mt_batch_size = multi_threaded_batch_size;
		library = lib;

		if (particle_array_allocator.arenas.current_size == 0) {
			//todo need to be able to adjust the arena size
			particle_array_allocator = CreateArenaManager(tfxMegabyte(2), 8);
			particle_arrays = tfxBucketArray<tfxParticleSoA>(&particle_array_allocator, 32);
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

		if (double_buffer_sprites)
			flags |= tfxEffectManagerFlags_double_buffer_sprites;

		for (tfxEachLayer) {
			max_cpu_particles_per_layer[layer] = layer_max_values[layer];

#ifdef tfxTRACK_MEMORY
			memcpy(sprites2d[layer].name, "ParticleManager::sprites2d\0", 27);
			memcpy(sprites3d[layer].name, "ParticleManager::sprites3d\0", 27);
#endif

			InitSpriteBothSoA(&sprite_buffer[0][layer], &sprites[0][layer], tfxMax((layer_max_values[layer] / tfxDataWidth + 1) * tfxDataWidth, tfxDataWidth * 2));
			if (flags & tfxEffectManagerFlags_double_buffer_sprites) {
				InitSpriteBothSoA(&sprite_buffer[1][layer], &sprites[1][layer], tfxMax((layer_max_values[layer] / tfxDataWidth + 1) * tfxDataWidth, tfxDataWidth * 2));
			}

			depth_indexes[layer][0].reserve(layer_max_values[layer]);
			depth_indexes[layer][1].reserve(layer_max_values[layer]);
		}

		flags |= dynamic_sprite_allocation ? tfxEffectManagerFlags_dynamic_sprite_allocation : 0;

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

	void tfxParticleManager::ClearAll(bool free_memory) {
		tfxCompleteAllWork(&work_queue);
		for (tfxEachLayer) {
			ClearSoABuffer(&sprite_buffer[0][layer]);
			if (flags & tfxEffectManagerFlags_double_buffer_sprites) {
				ClearSoABuffer(&sprite_buffer[1][layer]);
			}
			depth_indexes[layer][0].clear();
			depth_indexes[layer][1].clear();
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
				pm.particle_array_buffers[free_banks.back()].current_size = 0;
				return free_banks.pop_back();
			}
		}

		tfxParticleSoA lists;
		tfxU32 index = pm.particle_arrays.locked_push_back(lists);
		tfxSoABuffer buffer;
		buffer.resize_callback = tfxResizeParticleSoACallback;
		buffer.user_data = &pm.particle_arrays.back();
		pm.particle_array_buffers.push_back(buffer);
		assert(index == pm.particle_array_buffers.current_size - 1);
		InitParticleSoA(&pm.particle_array_buffers[index], &pm.particle_arrays.back(), reserve_amount);
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

	void SetEffectUserData(tfxParticleManager &pm, tfxU32 effect_index, void *data) {
		assert(effect_index < pm.effect_buffers.current_size);	//effect index is out of bounds of the array
		pm.effects.user_data[effect_index] = data;
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
		tfxLibrary *library = pm.effects.library[index];

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
				const float overal_scale = pm.effects.overal_scale[index];
				tfxU32 sprite_id = pm.GetParticleSpriteIndex(parent_particle_id);
				tfxU32 sprite_layer = (sprite_id & 0xF0000000) >> 28;
				tfxU32 sprite_index = sprite_id & 0x0FFFFFFF;
				if (sprite_id != tfxINVALID) {
					if (property_flags & tfxEmitterPropertyFlags_is_3d)
						TransformEffector3d(world_rotations, local_rotations, world_position, local_position, matrix, pm.sprites[!pm.current_sprite_buffer][sprite_layer].transform_3d[sprite_index], true, property_flags & tfxEmitterPropertyFlags_relative_angle);
					else
						TransformEffector2d(world_rotations, local_rotations, world_position, local_position, matrix, pm.sprites[!pm.current_sprite_buffer][sprite_layer].transform_2d[sprite_index], true, property_flags & tfxEmitterPropertyFlags_relative_angle);

					world_position += properties.emitter_handle[property_index] * overal_scale;
					if (state_flags & tfxEffectStateFlags_no_tween_this_update || state_flags & tfxEffectStateFlags_no_tween) {
						captured_position = world_position;
					}
				}
			}
			else {
				parent_particle_index = tfxINVALID;
				state_flags |= tfxEffectStateFlags_retain_matrix;
				local_position = world_position;
				local_rotations.roll = world_rotations.roll;
				state_flags |= tfxEffectStateFlags_stop_spawning;
			}
		}
		else {
			if (!(state_flags & tfxEffectStateFlags_retain_matrix)) {
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

			if (state_flags & tfxEffectStateFlags_no_tween_this_update || state_flags & tfxEffectStateFlags_no_tween) {
				captured_position = world_position;
			}
		}

		age += tfxFRAME_LENGTH;
		highest_particle_age -= tfxFRAME_LENGTH;

		if (properties.loop_length && age > properties.loop_length[property_index])
			age -= properties.loop_length[property_index];

		if (highest_particle_age <= 0 && age > tfxFRAME_LENGTH * 5.f) {
			timeout_counter++;
		}
		else {
			timeout_counter = 0;
		}

		if (state_flags & tfxEffectStateFlags_remove) {
			pm.emitters.timeout[index] = 1;
			highest_particle_age = 0;
		}

		state_flags &= ~tfxEffectStateFlags_no_tween_this_update;
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
		const tfxU32 particles_index = pm.emitters.particles_index[index];
		tfxLibrary *library = pm.library;

		captured_position = world_position;

		if (pm.lookup_mode == tfxPrecise) {
			frame = age;
		}
		else {
			frame = age / tfxLOOKUP_FREQUENCY;
		}

		spawn_work_entry->highest_particle_age = pm.emitters.highest_particle_age[index];

		tfxEmitterPropertiesSoA &properties = library->emitter_properties;

		assert(parent_index != tfxINVALID);	//Emitter must have a valid parent (an effect)

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
		state_flags |= parent_state_flags & tfxEffectStateFlags_no_tween;
		state_flags |= parent_state_flags & tfxEffectStateFlags_stop_spawning;
		state_flags |= parent_state_flags & tfxEffectStateFlags_remove;
		UpdateEmitterState(pm, index, parent_index, pm.effects.spawn_controls[parent_index], spawn_work_entry);

		bool is_compute = property_flags & tfxEmitterPropertyFlags_is_bottom_emitter && pm.flags & tfxEffectManagerFlags_use_compute_shader;
		tfxU32 amount_spawned = 0;
		tfxU32 max_spawn_count = NewSpritesNeeded(pm, index, parent_index, properties);

		tfxVec3 &parent_captured_position = pm.effects.captured_position[parent_index];
		tfxVec3 &parent_world_position = pm.effects.world_position[parent_index];
		tfxVec3 &parent_world_rotations = pm.effects.world_rotations[parent_index];
		tfxVec3 &parent_scale = pm.effects.scale[parent_index];
		tfxMatrix4 &parent_matrix = pm.effects.matrix[parent_index];

		if (property_flags & tfxEmitterPropertyFlags_is_3d) {
			tfxSoABuffer &sprite_buffer = pm.sprite_buffer[pm.current_sprite_buffer][layer];

			Transform3d(world_rotations, local_rotations, scale, world_position, local_position, translation, matrix, parent_world_rotations, parent_scale, parent_world_position, parent_matrix);

			if (state_flags & tfxEmitterStateFlags_no_tween_this_update || state_flags & tfxEmitterStateFlags_no_tween) {
				captured_position = world_position;
			}

			sprites_count = pm.particle_array_buffers[particles_index].current_size;
			if (pm.flags & tfxEffectManagerFlags_dynamic_sprite_allocation) {
				if (sprites_count + max_spawn_count > FreeSpace(&sprite_buffer)) {
					AddRows(&sprite_buffer, sprite_buffer.capacity + (sprites_count + max_spawn_count - FreeSpace(&sprite_buffer)) + 1, true);
					if (!(pm.flags & tfxEffectManagerFlags_unordered)) {
						pm.depth_indexes[layer][pm.current_depth_index_buffer].reserve(sprite_buffer.capacity);
					}
					sprite_buffer.current_size -= sprite_buffer.current_size - (sprites_count + max_spawn_count);
				}
				else {
					sprite_buffer.current_size += max_spawn_count + sprites_count;
				}
			}
			else {
				sprites_count = sprites_count > FreeSpace(&sprite_buffer) ? FreeSpace(&sprite_buffer) : sprites_count;
				sprite_buffer.current_size += max_spawn_count + sprites_count;
				max_spawn_count = max_spawn_count > FreeSpace(&sprite_buffer) ? FreeSpace(&sprite_buffer) : max_spawn_count;
			}

			sprites_count += max_spawn_count;
			sprites_index = pm.sprite_index_point[layer];
			pm.sprite_index_point[layer] += sprites_count;

			if (state_flags & tfxEmitterStateFlags_is_sub_emitter) {
				if (age > 0 && !(pm.flags & tfxEffectManagerFlags_disable_spawning))
					amount_spawned = SpawnParticles3d(pm, *spawn_work_entry, max_spawn_count);
			}
			else {
				if (!(pm.flags & tfxEffectManagerFlags_disable_spawning))
					amount_spawned = SpawnParticles3d(pm, *spawn_work_entry, max_spawn_count);
			}
			sprite_buffer.current_size -= (max_spawn_count - amount_spawned);
			pm.sprite_index_point[layer] -= (max_spawn_count - amount_spawned);
		}
		else {
			tfxSoABuffer &sprite_buffer = pm.sprite_buffer[pm.current_sprite_buffer][layer];

			Transform2d(world_rotations, local_rotations, scale, world_position, local_position, translation, matrix, parent_world_rotations, parent_scale, parent_world_position, parent_matrix);

			if (state_flags & tfxEmitterStateFlags_no_tween_this_update || state_flags & tfxEmitterStateFlags_no_tween) {
				captured_position = world_position;
			}

			sprites_count = pm.particle_array_buffers[particles_index].current_size;
			if (pm.flags & tfxEffectManagerFlags_dynamic_sprite_allocation) {
				if (sprites_count + max_spawn_count > FreeSpace(&sprite_buffer)) {
					AddRows(&sprite_buffer, sprite_buffer.capacity + (sprites_count + max_spawn_count - FreeSpace(&sprite_buffer)) + 1, true);
					if (!(pm.flags & tfxEffectManagerFlags_unordered)) {
						pm.depth_indexes[layer][pm.current_depth_index_buffer].reserve(sprite_buffer.capacity);
					}
					sprite_buffer.current_size -= sprite_buffer.current_size - (sprites_count + max_spawn_count);
				}
				else {
					sprite_buffer.current_size += max_spawn_count + sprites_count;
				}
			}
			else {
				sprites_count = sprites_count > FreeSpace(&sprite_buffer) ? FreeSpace(&sprite_buffer) : sprites_count;
				sprite_buffer.current_size += max_spawn_count + sprites_count;
				max_spawn_count = max_spawn_count > FreeSpace(&sprite_buffer) ? FreeSpace(&sprite_buffer) : max_spawn_count;
			}

			sprites_count += max_spawn_count;
			sprites_index = pm.sprite_index_point[layer];
			pm.sprite_index_point[layer] += sprites_count;

			if (state_flags & tfxEmitterStateFlags_is_sub_emitter) {
				if (age > 0 && !(pm.flags & tfxEffectManagerFlags_disable_spawning)) {
					amount_spawned = SpawnParticles2d(pm, *spawn_work_entry, max_spawn_count);
				}
			}
			else {
				if (!(pm.flags & tfxEffectManagerFlags_disable_spawning)) {
					amount_spawned = SpawnParticles2d(pm, *spawn_work_entry, max_spawn_count);
				}
			}
			sprite_buffer.current_size -= (max_spawn_count - amount_spawned);
			pm.sprite_index_point[layer] -= (max_spawn_count - amount_spawned);
		}

		age += tfxFRAME_LENGTH;
		if (!(property_flags & tfxEmitterPropertyFlags_single) || (property_flags & tfxEmitterPropertyFlags_single && properties.single_shot_limit[property_index] > 0) || state_flags & tfxEmitterStateFlags_stop_spawning) {
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

		if (state_flags & tfxEmitterStateFlags_remove) {
			pm.emitters.timeout[index] = 1;
			highest_particle_age = 0;
		}

		state_flags &= ~tfxEmitterStateFlags_no_tween_this_update;
	}

	tfxU32 NewSpritesNeeded(tfxParticleManager &pm, tfxU32 index, tfxU32 parent_index, tfxEmitterPropertiesSoA &properties) {

		const tfxEmitterStateFlags state_flags = pm.emitters.state_flags[index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[index];
		const tfxEmitterStateFlags parent_state_flags = pm.effects.state_flags[parent_index];
		const tfxU32 property_index = pm.emitters.properties_index[index];
		const float amount_remainder = pm.emitters.amount_remainder[index];
		float &spawn_quantity = pm.emitters.spawn_quantity[index];

		if (state_flags & tfxEmitterStateFlags_single_shot_done || parent_state_flags & tfxEffectStateFlags_stop_spawning)
			return 0;
		if (spawn_quantity == 0)
			return 0;

		const tfxVec3 emitter_size = pm.emitters.emitter_size[index];

		float step_size = 0;
		if (!(property_flags & tfxEmitterPropertyFlags_single)) {
			if (property_flags & tfxEmitterPropertyFlags_use_spawn_ratio && (properties.emission_type[property_index] == tfxArea || properties.emission_type[property_index] == tfxEllipse)) {
				if (property_flags & tfxEmitterPropertyFlags_is_3d) {
					float area = tfxMax(0.1f, emitter_size.x) * tfxMax(0.1f, emitter_size.y) * tfxMax(0.1f, emitter_size.z);
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

	void CompletePMWork(tfxParticleManager &pm) {
		tfxCompleteAllWork(&pm.work_queue);
	}

	tfxU32 SpawnParticles2d(tfxParticleManager &pm, tfxSpawnWorkEntry &work_entry, tfxU32 max_spawn_count) {
		const tfxEmitterPropertiesSoA &properties = *work_entry.properties;
		const tfxEmitterStateFlags state_flags = pm.emitters.state_flags[work_entry.emitter_index];
		const tfxEmitterStateFlags parent_state_flags = pm.emitters.state_flags[work_entry.emitter_index];
		const float spawn_quantity = pm.emitters.spawn_quantity[work_entry.emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[work_entry.emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[work_entry.emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[work_entry.emitter_index];
		const tfxU32 layer = properties.layer[property_index];
		float &qty_step_size = pm.emitters.qty_step_size[work_entry.emitter_index];
		float &amount_remainder = pm.emitters.amount_remainder[work_entry.emitter_index];

		if (state_flags & tfxEmitterStateFlags_single_shot_done || parent_state_flags & tfxEffectStateFlags_stop_spawning)
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

		bool grew = false;
		work_entry.spawn_start_index = AddRows(&pm.particle_array_buffers[particles_index], work_entry.amount_to_spawn, true, grew);
		tfxEmissionType &emission_type = properties.emission_type[pm.emitters.properties_index[work_entry.emitter_index]];
		if (grew && !(pm.flags & tfxEffectManagerFlags_unordered)) {
			//Todo: This should be avoided by allocating the correct amount for the particle buffer ahead of time
			//If the particle buffer is allocated a larger memory size then the ring buffer index has to be reset in the depth buffer list
			//to align them to the correct particle ids again.
			tfxParticleSoA &bank = pm.particle_arrays[particles_index];
			assert(pm.particle_array_buffers[particles_index].start_index == 0); //If start_index isn't 0 after the arrays grew then something went wrong with the allocation
			for (int i = 0; i != work_entry.spawn_start_index; ++i) {
				pm.depth_indexes[layer][pm.current_depth_index_buffer][bank.depth_index[i]].particle_id = MakeParticleID(particles_index, i);
			}
		}

		if (!(pm.flags & tfxEffectManagerFlags_update_age_only) && !(pm.flags & tfxEffectManagerFlags_single_threaded) && tfxNumberOfThreadsInAdditionToMain) {
			if (work_entry.amount_to_spawn > 0) {
				work_entry.end_index = work_entry.amount_to_spawn;
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
				//Update: I did try this and had a work queue that you could have sempahores and get jobs to run when specific jobs had finished, but it was actually slower.
				//This is probably because my implementation was bad and the fact that I had to use mutexes into order to manage adding jobs working in a thread friendly manner
				//(multi producer/worker). Back to the drawing board, will give it another go at some point.
				tfxCompleteAllWork(&pm.work_queue);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleMicroUpdate2d);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleAge);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleNoise);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleImageFrame);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleSize2d);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleSpin2d);
			}
		}
		else if (!(pm.flags & tfxEffectManagerFlags_update_age_only)) {
			if (work_entry.amount_to_spawn > 0) {
				work_entry.end_index = work_entry.amount_to_spawn;
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
				SpawnParticleAge(&pm.work_queue, &work_entry);
				SpawnParticleNoise(&pm.work_queue, &work_entry);
				SpawnParticleImageFrame(&pm.work_queue, &work_entry);
				SpawnParticleSize2d(&pm.work_queue, &work_entry);
				SpawnParticleSpin2d(&pm.work_queue, &work_entry);
			}
		}
		else {
			SpawnParticleAge(&pm.work_queue, &work_entry);
		}

		if (work_entry.amount_to_spawn > 0 && property_flags & tfxEmitterPropertyFlags_single)
			pm.emitters.state_flags[work_entry.emitter_index] |= tfxEmitterStateFlags_single_shot_done;

		return work_entry.amount_to_spawn;
	}

	tfxU32 SpawnParticles3d(tfxParticleManager &pm, tfxSpawnWorkEntry &work_entry, tfxU32 max_spawn_count) {
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

		if (state_flags & tfxEmitterStateFlags_single_shot_done || parent_state_flags & tfxEffectStateFlags_stop_spawning)
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

		if (pm.flags & tfxEffectManagerFlags_order_by_depth) {
			//We must complete all work first before potentially growing the particle_array_buffers as some threads may still be working in the buffer
			tfxCompleteAllWork(&pm.work_queue);
		}
		bool grew = false;
		tfxU32 start_index = pm.particle_array_buffers[particles_index].start_index;
		work_entry.spawn_start_index = AddRows(&pm.particle_array_buffers[particles_index], work_entry.amount_to_spawn, true, grew);
		if (grew && !(pm.flags & tfxEffectManagerFlags_unordered) && start_index > 0) {
			//Todo: This should be avoided by allocating the correct amount for the particle buffer ahead of time
			//If the particle buffer is allocated a larger memory size then the ring buffer index has to be reset in the depth buffer list
			//to align them to the correct particle ids again.
			tfxParticleSoA &bank = pm.particle_arrays[particles_index];
			//assert(pm.particle_array_buffers[particles_index].start_index == 0); //If start_index isn't 0 after the arrays grew then something went wrong with the allocation
			for (int i = 0; i != pm.particle_array_buffers[particles_index].current_size - work_entry.amount_to_spawn; ++i) {
				tfxU32 depth_index = bank.depth_index[i];
				tfxU32 depth_id = pm.depth_indexes[layer][pm.current_depth_index_buffer][bank.depth_index[i]].particle_id;
				pm.depth_indexes[layer][pm.current_depth_index_buffer][bank.depth_index[i]].particle_id = MakeParticleID(particles_index, i);
			}
		}
		tfxEmissionType &emission_type = properties.emission_type[property_index];

		if (!(pm.flags & tfxEffectManagerFlags_update_age_only) && !(pm.flags & tfxEffectManagerFlags_single_threaded) && tfxNumberOfThreadsInAdditionToMain) {
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
				//Can maybe revisit this. We have to complete the above work before doing the micro update. I would like to add the micro update from one of the above threadsrandom.Advance();
				//when all 4 have finished but synchronisation is hard to get right. Would have to rethink for a multi producer work queue. For now though this is working
				//fine and is stable
				//Another idea: We could split this into 2 functions where where we push all the above spawn functions to the queue for all emitters, then in the second function
				//we do all of the micro updates for all particles
				tfxCompleteAllWork(&pm.work_queue);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleMicroUpdate3d);  
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleAge);  
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleNoise);  
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleImageFrame);  
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleSize3d);  
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, SpawnParticleSpin3d);  
			}
		}
		else if (!(pm.flags & tfxEffectManagerFlags_update_age_only)) {
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
		}
		else {
			SpawnParticleAge(&pm.work_queue, &work_entry);
		}

		if (work_entry.amount_to_spawn > 0 && property_flags & tfxEmitterPropertyFlags_single)
			pm.emitters.state_flags[work_entry.emitter_index] |= tfxEmitterStateFlags_single_shot_done;

		return work_entry.amount_to_spawn;
	}

	void SpawnParticleAge(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;

		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(2 + pm.emitters.seed_index[emitter_index]);

		const float life = pm.emitters.life[emitter_index];
		const float life_variation = pm.emitters.life_variation[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxU32 sprites_index = pm.emitters.sprites_index[emitter_index];
		const float highest_particle_age = pm.emitters.highest_particle_age[emitter_index];
		const tfxU32 loop_count = entry->properties->single_shot_limit[property_index] + 1;
		tfxLibrary *library = pm.library;
		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[emitter_index];
		const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];
		const float emitter_intensity = pm.emitters.intensity[emitter_index];
		const float first_red_value = library->emitter_attributes[emitter_attributes].overtime.red.GetFirstValue();
		const float first_green_value = library->emitter_attributes[emitter_attributes].overtime.green.GetFirstValue();
		const float first_blue_value = library->emitter_attributes[emitter_attributes].overtime.blue.GetFirstValue();
		const float first_alpha_value = library->emitter_attributes[emitter_attributes].overtime.blendfactor.GetFirstValue();
		const float first_intensity_value = library->emitter_attributes[emitter_attributes].overtime.intensity.GetFirstValue();

		assert(random.seeds[0] > 0);

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &age = entry->particle_data->age[index];
			float &max_age = entry->particle_data->max_age[index];
			tfxRGBA8 &color = entry->particle_data->color[index];
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
			max_age = life + random.Range(life_variation);
			single_loop_count = 0;

			float alpha = 255.f * first_alpha_value;
			float intensity = first_intensity_value * emitter_intensity;
			//intensity = 0.f;
			if (emitter_flags & tfxEmitterStateFlags_random_color) {
				float age = random.Range(max_age);
				color = tfxRGBA8(	255.f * lookup_overtime_callback(library->emitter_attributes[emitter_attributes].overtime.red, age, max_age),
									255.f * lookup_overtime_callback(library->emitter_attributes[emitter_attributes].overtime.green, age, max_age),
									255.f * lookup_overtime_callback(library->emitter_attributes[emitter_attributes].overtime.blue, age, max_age), alpha);
			}
			else {
				color = tfxRGBA8(255.f * first_red_value, 255.f * first_green_value, 255.f * first_blue_value, alpha);
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
		tfxPROFILE;

		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		const tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(3 + pm.emitters.seed_index[emitter_index]);
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {

			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &image_frame = entry->particle_data->image_frame[index];
			tfxU32 &sprites_index = entry->particle_data->sprite_index[index];
			sprites_index = tfxINVALID;

			//----Image
			//data.image = GetProperties().image;
			if (property_flags & tfxEmitterPropertyFlags_random_start_frame && properties.image[property_index]->animation_frames > 1) {
				image_frame = random.Range(properties.image[property_index]->animation_frames);
			}
			else {
				image_frame = properties.start_frame[property_index];
			}

		}

	}

	void SpawnParticleSize2d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;

		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(4 + pm.emitters.seed_index[emitter_index]);
		const tfxEmitterPropertiesSoA &properties = *entry->properties;

		const tfxVec2 size = pm.emitters.size[emitter_index];
		const tfxVec2 size_variation = pm.emitters.size_variation[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		tfxVec2 &image_size = properties.image[property_index]->image_size;

		for (int i = 0; i != entry->amount_to_spawn; ++i) {

			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &base_size_x = entry->particle_data->base_size_x[index];
			float &base_size_y = entry->particle_data->base_size_y[index];

			//----Size
			if (!(property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
				float random_size_x = random.Range(size_variation.x);
				float random_size_y = random.Range(size_variation.y);
				base_size_y = (random_size_y + size.y) / image_size.y;
				base_size_x = (random_size_x + size.x) / image_size.x;
			}
			else {
				float random_size_x = random.Range(size_variation.x);
				float random_size_y = random_size_x;
				base_size_y = (random_size_y + size.y) / image_size.y;
				base_size_x = (random_size_x + size.x) / image_size.x;
			}
		}

	}

	void SpawnParticleSize3d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;

		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(5 + pm.emitters.seed_index[emitter_index]);
		const tfxEmitterPropertiesSoA &properties = *entry->properties;

		const tfxVec2 size = pm.emitters.size[emitter_index];
		const tfxVec2 size_variation = pm.emitters.size_variation[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		tfxVec2 &image_size = properties.image[property_index]->image_size;

		for (int i = 0; i != entry->amount_to_spawn; ++i) {

			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &base_size_x = entry->particle_data->base_size_x[index];
			float &base_size_y = entry->particle_data->base_size_y[index];

			//----Size
			if (!(property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
				float random_size_x = random.Range(size_variation.x);
				float random_size_y = random.Range(size_variation.y);
				base_size_y = random_size_y + size.y;
				base_size_x = random_size_x + size.x;
			}
			else {
				float random_size_x = random.Range(size_variation.x);
				float random_size_y = random_size_x;
				base_size_y = random_size_y + size.y;
				base_size_x = random_size_x + size.x;
			}

		}

	}

	void SpawnParticleNoise(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;

		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(6 + pm.emitters.seed_index[emitter_index]);
		const float emitter_noise_offset_variation = pm.emitters.noise_offset_variation[emitter_index];
		const float emitter_noise_offset = pm.emitters.noise_offset[emitter_index];
		const float emitter_noise_resolution = pm.emitters.noise_resolution[emitter_index];
		const float parent_noise_base_offset = pm.effects.noise_base_offset[pm.emitters.parent_index[emitter_index]];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {

			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &noise_offset = entry->particle_data->noise_offset[index];
			float &noise_resolution = entry->particle_data->noise_resolution[index];

			//----Motion randomness
			noise_offset = random.Range(emitter_noise_offset_variation) + emitter_noise_offset + parent_noise_base_offset;
			noise_resolution = emitter_noise_resolution + 0.01f;

		}
	}

	void SpawnParticleSpin2d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;

		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(7 + pm.emitters.seed_index[emitter_index]);

		const float spin_variation = pm.emitters.spin_variation[emitter_index];
		const float spin = pm.emitters.spin[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {

			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &base_spin = entry->particle_data->base_spin[index];

			//----Spin
			base_spin = random.Range(-spin_variation, spin_variation) + spin;

		}

	}

	void SpawnParticleSpin3d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;

		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(8 + pm.emitters.seed_index[emitter_index]);

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
			float &local_rotations_x = entry->particle_data->local_rotations_x[index];
			float &local_rotations_y = entry->particle_data->local_rotations_y[index];
			float &local_rotations_z = entry->particle_data->local_rotations_z[index];

			//----Spin
			base_spin = random.Range(-spin_variation, spin_variation) + spin;
			if (angle_settings & tfxAngleSettingFlags_specify_roll)
				local_rotations_z = angle_offsets.roll;
			if (angle_settings & tfxAngleSettingFlags_specify_pitch)
				local_rotations_x = angle_offsets.pitch;
			if (angle_settings & tfxAngleSettingFlags_specify_yaw)
				local_rotations_y = angle_offsets.yaw;
			if (angle_settings & tfxAngleSettingFlags_random_pitch)
				local_rotations_x = random.Range(angle_offsets.pitch);
			if (angle_settings & tfxAngleSettingFlags_random_yaw)
				local_rotations_y = random.Range(angle_offsets.yaw);
			if (angle_settings & tfxAngleSettingFlags_random_roll)
				local_rotations_z = random.Range(angle_offsets.roll);

		}

	}

	void SpawnParticlePoint2d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(9 + pm.emitters.seed_index[emitter_index]);
		tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &local_position_x = entry->particle_data->position_x[index];
			float &local_position_y = entry->particle_data->position_y[index];

			local_position_x = local_position_y = 0;
			tfxVec2 lerp_position = InterpolateVec2(tween, emitter_captured_position.xy(), emitter_world_position.xy());
			if (property_flags & tfxEmitterPropertyFlags_relative_position)
				local_position_x = local_position_y = 0;
			else {
				if (property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) {
					local_position_x = lerp_position.x;
					local_position_y = lerp_position.y;
				}
				else {
					tfxVec2 rotvec = mmTransformVector(matrix, -handle.xy());
					local_position_x = rotvec.x + lerp_position.x;
					local_position_y = rotvec.y + lerp_position.y;
				}
			}
			tween += entry->qty_step_size;
		}

	}

	void SpawnParticlePoint3d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(10 + pm.emitters.seed_index[emitter_index]);
		tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &local_position_x = entry->particle_data->position_x[index];
			float &local_position_y = entry->particle_data->position_y[index];
			float &local_position_z = entry->particle_data->position_z[index];

			local_position_x = local_position_y = local_position_z = 0;
			tfxVec3 lerp_position = InterpolateVec3(tween, emitter_captured_position, emitter_world_position);
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				if (property_flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) {
					local_position_x = lerp_position.x;
					local_position_y = lerp_position.y;
					local_position_z = lerp_position.z;
				}
				else {
					tfxVec3 rotvec = mmTransformVector3(matrix, -handle);
					local_position_x = rotvec.x + lerp_position.x;
					local_position_y = rotvec.y + lerp_position.y;
					local_position_z = rotvec.z + lerp_position.z;
				}
			}
			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleLine2d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(11 + pm.emitters.seed_index[emitter_index]);
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
			float &local_position_x = entry->particle_data->position_x[index];
			float &local_position_y = entry->particle_data->position_y[index];

			local_position_x = local_position_y = 0;
			tfxVec2 lerp_position = InterpolateVec2(tween, emitter_captured_position.xy(), emitter_world_position.xy());

			if (property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				grid_coords.x = 0.f;

				if (!(property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					grid_coords.y--;
					if (grid_coords.y < 0.f) {
						grid_coords.y = grid_points.x - 1;
					}
				}

				local_position_x = grid_coords.x * -grid_segment_size.x;
				local_position_y = grid_coords.y * -grid_segment_size.y;

				if (property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					grid_coords.y++;
					if (grid_coords.y >= grid_points.x) {
						grid_coords.y = 0.f;
					}
				}

			}
			else {
				local_position_x = 0.f;
				local_position_y = random.Range(-emitter_size.y, 0.f);

			}

			//----TForm and Emission
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position) && !(property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
				tfxVec2 pos = mmTransformVector(matrix, tfxVec2(local_position_x, local_position_y) + handle.xy());
				local_position_x = lerp_position.x + pos.x * scale.x;
				local_position_y = lerp_position.y + pos.y * scale.y;
			}

			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleLine3d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(12 + pm.emitters.seed_index[emitter_index]);
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
			float &local_position_x = entry->particle_data->position_x[index];
			float &local_position_y = entry->particle_data->position_y[index];
			float &local_position_z = entry->particle_data->position_z[index];

			local_position_x = local_position_y = local_position_z = 0;
			tfxVec3 lerp_position = InterpolateVec3(tween, emitter_captured_position, emitter_world_position);

			if (property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				grid_coords.x = 0.f;
				grid_coords.z = 0.f;

				if (property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
					grid_coords.y = (float)random.Range((tfxU32)grid_points.x);
					local_position_x = grid_coords.x * grid_segment_size.x;
					local_position_y = grid_coords.y * grid_segment_size.y;
					local_position_z = grid_coords.z * grid_segment_size.z;
				}
				else {

					if (!(property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
						grid_coords.y--;
						if (grid_coords.y < 0.f) {
							grid_coords.y = grid_points.x - 1;
						}
					}

					local_position_x = grid_coords.x * grid_segment_size.x;
					local_position_y = grid_coords.y * grid_segment_size.y;
					local_position_z = grid_coords.z * grid_segment_size.z;

					if (property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
						grid_coords.y++;
						if (grid_coords.y >= grid_points.x) {
							grid_coords.y = 0.f;
						}
					}
				}

			}
			else {
				local_position_x = 0.f;
				local_position_y = random.Range(0.f, emitter_size.y);
				local_position_z = 0.f;
			}

			//----TForm and Emission
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position) && !(property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
				tfxVec3 pos = mmTransformVector3(matrix, tfxVec3(local_position_x, local_position_y, local_position_z) + handle);
				local_position_x = lerp_position.x + pos.x * scale.x;
				local_position_y = lerp_position.y + pos.y * scale.y;
				local_position_z = lerp_position.z + pos.z * scale.z;
			}

			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleArea2d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(13 + pm.emitters.seed_index[emitter_index]);
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
			float &local_position_x = entry->particle_data->position_x[index];
			float &local_position_y = entry->particle_data->position_y[index];

			local_position_x = local_position_y = 0;
			tfxVec2 lerp_position = InterpolateVec2(tween, emitter_captured_position.xy(), emitter_world_position.xy());

			tfxVec2 position = tfxVec2(0.f, 0.f);

			if (property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				if (property_flags & tfxEmitterPropertyFlags_fill_area) {
					if (property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
						grid_coords.x = (float)random.Range((tfxU32)grid_points.x);
						grid_coords.y = (float)random.Range((tfxU32)grid_points.y);
						local_position_x = grid_coords.x * grid_segment_size.x;
						local_position_y = grid_coords.y * grid_segment_size.y;
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

						local_position_x = grid_coords.x * grid_segment_size.x;
						local_position_y = grid_coords.y * grid_segment_size.y;

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
						tfxU32 side = random.Range((tfxU32)4);
						if (side == 0) {
							//left side
							grid_coords.x = 0.f;
							grid_coords.y = (float)random.Range((tfxU32)grid_points.y);
						}
						else if (side == 1) {
							//right side
							grid_coords.x = grid_points.x - 1;
							grid_coords.y = (float)random.Range((tfxU32)grid_points.y);
						}
						else if (side == 2) {
							//top side
							grid_coords.x = (float)random.Range((tfxU32)grid_points.x);
							grid_coords.y = 0.f;
						}
						else if (side == 3) {
							//bottom side
							grid_coords.x = (float)random.Range((tfxU32)grid_points.x);
							grid_coords.y = grid_points.y - 1;
						}
						local_position_x = grid_coords.x * grid_segment_size.x;
						local_position_y = grid_coords.y * grid_segment_size.y;
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
						local_position_x = position.x + (grid_coords.x * grid_segment_size.x);
						local_position_y = position.y + (grid_coords.y * grid_segment_size.y);
					}
				}
			}
			else {
				if (property_flags & tfxEmitterPropertyFlags_fill_area) {
					position.x = random.Range(emitter_size.x);
					position.y = random.Range(emitter_size.y);
				}
				else {
					//Spawn on one of 4 edges of the area
					tfxU32 side = random.Range((tfxU32)4);
					if (side == 0) {
						//left side
						position.x = 0.f;
						position.y = random.Range(emitter_size.y);
					}
					else if (side == 1) {
						//right side
						position.x = emitter_size.x;
						position.y = random.Range(emitter_size.y);
					}
					else if (side == 2) {
						//top side
						position.x = random.Range(emitter_size.x);
						position.y = 0.f;
					}
					else if (side == 3) {
						//bottom side
						position.x = random.Range(emitter_size.x);
						position.y = emitter_size.y;
					}
				}

				local_position_x = position.x;
				local_position_y = position.y;
			}

			//----TForm and Emission
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				tfxVec2 pos = mmTransformVector(matrix, tfxVec2(local_position_x, local_position_y) + handle.xy());
				local_position_x = lerp_position.x + pos.x * scale.x;
				local_position_y = lerp_position.y + pos.y * scale.y;
			}

			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleArea3d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(14 + pm.emitters.seed_index[emitter_index]);
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
			float &local_position_x = entry->particle_data->position_x[index];
			float &local_position_y = entry->particle_data->position_y[index];
			float &local_position_z = entry->particle_data->position_z[index];

			local_position_x = local_position_y = local_position_z = 0;
			tfxVec3 lerp_position = InterpolateVec3(tween, emitter_captured_position, emitter_world_position);

			tfxVec3 position;

			if (property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {

				if (property_flags & tfxEmitterPropertyFlags_fill_area) {

					if (property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
						grid_coords.x = (float)random.Range((tfxU32)grid_points.x);
						grid_coords.y = (float)random.Range((tfxU32)grid_points.y);
						grid_coords.z = (float)random.Range((tfxU32)grid_points.z);

						local_position_x = grid_coords.x * grid_segment_size.x;
						local_position_y = grid_coords.y * grid_segment_size.y;
						local_position_z = grid_coords.z * grid_segment_size.z;
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

						local_position_x = position.x + (grid_coords.x * grid_segment_size.x);
						local_position_y = position.y + (grid_coords.y * grid_segment_size.y);
						local_position_z = position.z + (grid_coords.z * grid_segment_size.z);

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
						tfxU32 side = random.Range((property_flags & tfxEmitterPropertyFlags_area_open_ends) ? (tfxU32)4 : (tfxU32)6);
						if (side == 0) {
							//left side
							grid_coords.x = 0.f;
							grid_coords.y = (float)random.Range((tfxU32)grid_points.y);
							grid_coords.z = (float)random.Range((tfxU32)grid_points.z);
						}
						else if (side == 1) {
							//right side
							grid_coords.x = grid_points.x - 1;
							grid_coords.y = (float)random.Range((tfxU32)grid_points.y);
							grid_coords.z = (float)random.Range((tfxU32)grid_points.z);
						}
						else if (side == 2) {
							//top side
							grid_coords.x = (float)random.Range((tfxU32)grid_points.x);
							grid_coords.y = 0.f;
							grid_coords.z = (float)random.Range((tfxU32)grid_points.z);
						}
						else if (side == 3) {
							//bottom side
							grid_coords.x = (float)random.Range((tfxU32)grid_points.x);
							grid_coords.y = grid_points.y - 1;
							grid_coords.z = (float)random.Range((tfxU32)grid_points.z);
						}
						else if (side == 4) {
							//End far
							grid_coords.x = (float)random.Range((tfxU32)grid_points.x);
							grid_coords.y = (float)random.Range((tfxU32)grid_points.y);
							grid_coords.z = grid_points.z - 1;
						}
						else if (side == 5) {
							//End near
							grid_coords.x = (float)random.Range((tfxU32)grid_points.x);
							grid_coords.y = (float)random.Range((tfxU32)grid_points.y);
							grid_coords.z = 0.f;
						}
						local_position_x = grid_coords.x * grid_segment_size.x;
						local_position_y = grid_coords.y * grid_segment_size.y;
						local_position_z = grid_coords.z * grid_segment_size.z;
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
						local_position_x = position.x + (grid_coords.x * grid_segment_size.x);
						local_position_y = position.y + (grid_coords.y * grid_segment_size.y);
						local_position_z = position.z + (grid_coords.z * grid_segment_size.z);
					}

				}
			}
			else {
				if (property_flags & tfxEmitterPropertyFlags_fill_area) {
					position.x = random.Range(emitter_size.x);
					position.y = random.Range(emitter_size.y);
					position.z = random.Range(emitter_size.z);
				}
				else {
					//Spawn on one of 6 edges of the cuboid
					tfxU32 side = random.Range((property_flags & tfxEmitterPropertyFlags_area_open_ends) ? (tfxU32)4 : (tfxU32)6);
					if (side == 0) {
						//left side
						position.x = 0.f;
						position.y = random.Range(emitter_size.y);
						position.z = random.Range(emitter_size.z);
					}
					else if (side == 1) {
						//right side
						position.x = emitter_size.x;
						position.y = random.Range(emitter_size.y);
						position.z = random.Range(emitter_size.z);
					}
					else if (side == 2) {
						//top side
						position.x = random.Range(emitter_size.x);
						position.y = 0.f;
						position.z = random.Range(emitter_size.z);
					}
					else if (side == 3) {
						//bottom side
						position.x = random.Range(emitter_size.x);
						position.y = emitter_size.y;
						position.z = random.Range(emitter_size.z);
					}
					else if (side == 4) {
						//End far
						position.x = random.Range(emitter_size.x);
						position.y = random.Range(emitter_size.y);
						position.z = emitter_size.z;
					}
					else if (side == 5) {
						//End near
						position.x = random.Range(emitter_size.x);
						position.y = random.Range(emitter_size.y);
						position.z = 0.f;
					}
				}

				local_position_x = position.x;
				local_position_y = position.y;
				local_position_z = position.z;
			}

			//----TForm and Emission
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				tfxVec3 pos = mmTransformVector3(matrix, tfxVec3(local_position_x, local_position_y, local_position_z) + handle);
				local_position_x = lerp_position.x + pos.x * scale.x;
				local_position_y = lerp_position.y + pos.y * scale.y;
				local_position_z = lerp_position.z + pos.z * scale.z;
			}

			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleEllipse2d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(15 + pm.emitters.seed_index[emitter_index]);
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
			float &local_position_x = entry->particle_data->position_x[index];
			float &local_position_y = entry->particle_data->position_y[index];

			local_position_x = local_position_y = 0;
			tfxVec2 lerp_position = InterpolateVec2(tween, emitter_captured_position.xy(), emitter_world_position.xy());

			tfxVec2 half_emitter_size = (emitter_size * .5f);
			tfxVec2 position = tfxVec2(0.f, 0.f);

			if (property_flags & tfxEmitterPropertyFlags_spawn_on_grid && !(property_flags & tfxEmitterPropertyFlags_fill_area) && !(property_flags & tfxEmitterPropertyFlags_grid_spawn_random)) {

				grid_coords.y = 0.f;

				if (property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
					grid_coords.x--;
					if (grid_coords.x < 0.f) {
						grid_coords.x = grid_points.x - 1;
					}
				}

				float th = grid_coords.x * grid_segment_size.x + arc_offset;
				local_position_x = std::cosf(th) * half_emitter_size.x + half_emitter_size.x;
				local_position_y = -std::sinf(th) * half_emitter_size.y + half_emitter_size.y;

				if (!(property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
					grid_coords.x++;
					if (grid_coords.x >= grid_points.x) {
						grid_coords.x = 0.f;
					}
				}

			}
			else if (property_flags & tfxEmitterPropertyFlags_spawn_on_grid && !(property_flags & tfxEmitterPropertyFlags_fill_area) && property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
				float th = (float)random.Range((tfxU32)grid_points.x) * grid_segment_size.x + arc_offset;
				local_position_x = std::cosf(th) * half_emitter_size.x + half_emitter_size.x;
				local_position_x = -std::sinf(th) * half_emitter_size.y + half_emitter_size.y;
			}
			else if (!(property_flags & tfxEmitterPropertyFlags_fill_area)) {
				float th = random.Range(arc_size) + arc_offset;

				local_position_x = std::cosf(th) * half_emitter_size.x + half_emitter_size.x;
				local_position_y = -std::sinf(th) * half_emitter_size.y + half_emitter_size.y;

			}
			else {
				local_position_x = random.Range(0.f, emitter_size.x);
				local_position_y = random.Range(0.f, emitter_size.y);

				while ((std::pow(local_position_x - half_emitter_size.x, 2) / std::pow(half_emitter_size.x, 2)) + (std::pow(local_position_y - half_emitter_size.y, 2) / std::pow(half_emitter_size.y, 2)) > 1) {
					local_position_x = random.Range(0.f, emitter_size.x);
					local_position_y = random.Range(0.f, emitter_size.y);
				}
			}

			//----TForm and Emission
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				tfxVec2 pos = mmTransformVector(matrix, tfxVec2(local_position_x, local_position_y) + handle.xy());
				local_position_x = lerp_position.x + pos.x * scale.x;
				local_position_y = lerp_position.y + pos.y * scale.y;
			}


			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleEllipse3d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		random.AlterSeed(16 + pm.emitters.seed_index[emitter_index]);
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
			float &local_position_x = entry->particle_data->position_x[index];
			float &local_position_y = entry->particle_data->position_y[index];
			float &local_position_z = entry->particle_data->position_z[index];

			local_position_x = local_position_y = local_position_z = 0;
			tfxVec3 lerp_position = InterpolateVec3(tween, emitter_captured_position, emitter_world_position);

			tfxVec3 half_emitter_size = emitter_size * .5f;
			tfxVec3 position;

			if (!(property_flags & tfxEmitterPropertyFlags_fill_area)) {
				float u = random.Range(1.f);
				float v = random.Range(1.f);
				float theta = u * 2.f * tfxPI;
				float phi = std::acosf(2.f * v - 1.f);
				float sin_theta = std::sinf(theta);
				float cos_theta = std::cosf(theta);
				float sin_phi = std::sinf(phi);
				float cos_phi = std::cosf(phi);
				local_position_x = half_emitter_size.x * sin_phi * cos_theta;
				local_position_y = half_emitter_size.y * sin_phi * sin_theta;
				local_position_z = half_emitter_size.z * cos_phi;
			}
			else {
				position.x = random.Range(-half_emitter_size.x, half_emitter_size.x);
				position.y = random.Range(-half_emitter_size.y, half_emitter_size.y);
				position.z = random.Range(-half_emitter_size.z, half_emitter_size.z);

				while (std::powf(position.x / half_emitter_size.x, 2.f) + std::powf(position.y / half_emitter_size.y, 2.f) + std::powf(position.z / half_emitter_size.z, 2.f) > 1.f) {
					position.x = random.Range(-half_emitter_size.x, half_emitter_size.x);
					position.y = random.Range(-half_emitter_size.y, half_emitter_size.y);
					position.z = random.Range(-half_emitter_size.z, half_emitter_size.z);
				}

				local_position_x = position.x;
				local_position_y = position.y;
				local_position_z = position.z;
			}

			//----TForm and Emission
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				tfxVec3 pos = mmTransformVector3(matrix, tfxVec3(local_position_x, local_position_y, local_position_z) + handle);
				local_position_x = lerp_position.x + pos.x * scale.x;
				local_position_y = lerp_position.y + pos.y * scale.y;
				local_position_z = lerp_position.z + pos.z * scale.z;
			}

			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleIcosphere3d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(17 + pm.emitters.seed_index[emitter_index]);
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
			float &local_position_x = entry->particle_data->position_x[index];
			float &local_position_y = entry->particle_data->position_y[index];
			float &local_position_z = entry->particle_data->position_z[index];

			local_position_x = local_position_y = local_position_z = 0;
			tfxVec3 lerp_position = InterpolateVec3(tween, emitter_captured_position, emitter_world_position);

			tfxU32 sub_division = (tfxU32)grid_points.x;
			assert(sub_division < 6);	//Make sure that grid_points.x is set to 0-5 as that is used for the sub divisions array index

			if (grid_coords.x >= tfxIcospherePoints[sub_division].current_size) {
				grid_coords.x = 0;
			}
			local_position_x = tfxIcospherePoints[sub_division][(tfxU32)grid_coords.x].x * half_emitter_size.x;
			local_position_y = tfxIcospherePoints[sub_division][(tfxU32)grid_coords.x].y * half_emitter_size.y;
			local_position_z = tfxIcospherePoints[sub_division][(tfxU32)grid_coords.x].z * half_emitter_size.z;
			grid_coords.x++;

			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				tfxVec3 pos = mmTransformVector3(matrix, tfxVec3(local_position_x, local_position_y, local_position_z) + handle);
				local_position_x = lerp_position.x + local_position_x * scale.x;
				local_position_y = lerp_position.y + local_position_y * scale.y;
				local_position_z = lerp_position.z + local_position_z * scale.z;
			}

			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleIcosphereRandom3d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(18 + pm.emitters.seed_index[emitter_index]);
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
			float &local_position_x = entry->particle_data->position_x[index];
			float &local_position_y = entry->particle_data->position_y[index];
			float &local_position_z = entry->particle_data->position_z[index];

			local_position_x = local_position_y = local_position_z = 0;
			tfxVec3 lerp_position = InterpolateVec3(tween, emitter_captured_position, emitter_world_position);

			tfxVec3 half_emitter_size = emitter_size * .5f;
			tfxU32 sub_division = (tfxU32)grid_points.x;
			assert(sub_division < 6);	//Make sure that grid_points.x is set to 0-5 as that is used for the sub divisions array index
			int ico_point = random.Range(tfxIcospherePoints[sub_division].current_size);
			local_position_x = tfxIcospherePoints[sub_division][ico_point].x * half_emitter_size.x;
			local_position_y = tfxIcospherePoints[sub_division][ico_point].y * half_emitter_size.y;
			local_position_z = tfxIcospherePoints[sub_division][ico_point].z * half_emitter_size.z;
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				tfxVec3 pos = mmTransformVector3(matrix, tfxVec3(local_position_x, local_position_y, local_position_z) + handle);
				local_position_x = lerp_position.x + local_position_x * scale.x;
				local_position_y = lerp_position.y + local_position_y * scale.y;
				local_position_z = lerp_position.z + local_position_z * scale.z;
			}

			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleCylinder3d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(19 + pm.emitters.seed_index[emitter_index]);
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
			float &local_position_x = entry->particle_data->position_x[index];
			float &local_position_y = entry->particle_data->position_y[index];
			float &local_position_z = entry->particle_data->position_z[index];

			local_position_x = local_position_y = local_position_z = 0;
			tfxVec3 lerp_position = InterpolateVec3(tween, emitter_captured_position, emitter_world_position);

			if (property_flags & tfxEmitterPropertyFlags_spawn_on_grid && !(property_flags & tfxEmitterPropertyFlags_fill_area)) {

				grid_coords.z = 0.f;

				if (property_flags & tfxEmitterPropertyFlags_grid_spawn_random) {
					grid_coords.x = (float)random.Range((tfxU32)grid_points.x);
					grid_coords.y = (float)random.Range((tfxU32)grid_points.y);

					float th = grid_coords.x * grid_segment_size.x + arc_offset;
					local_position_x = std::cosf(th) * half_emitter_size.x + half_emitter_size.x;
					local_position_y = grid_coords.y * grid_segment_size.y;
					local_position_z = -std::sinf(th) * half_emitter_size.z + half_emitter_size.z;
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
					local_position_x = std::cosf(th) * half_emitter_size.x + half_emitter_size.x;
					local_position_y = grid_coords.y * grid_segment_size.y;
					local_position_z = -std::sinf(th) * half_emitter_size.z + half_emitter_size.z;

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
				float th = random.Range(arc_size) + arc_offset;

				local_position_x = std::cosf(th) * half_emitter_size.x + half_emitter_size.x;
				local_position_y = random.Range(half_emitter_size.y);
				local_position_z = -std::sinf(th) * half_emitter_size.z + half_emitter_size.z;
			}
			else {
				local_position_x = random.Range(0.f, half_emitter_size.x);
				local_position_y = random.Range(0.f, half_emitter_size.y);
				local_position_z = random.Range(0.f, half_emitter_size.z);

				while ((std::pow(local_position_x - half_emitter_size.x, 2) / std::pow(half_emitter_size.x, 2)) + (std::pow(local_position_z - half_emitter_size.z, 2) / std::pow(half_emitter_size.z, 2)) > 1) {
					local_position_x = random.Range(0.f, half_emitter_size.x);
					local_position_z = random.Range(0.f, half_emitter_size.z);
				}
			}

			//----TForm and Emission
			if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				tfxVec3 pos = mmTransformVector3(matrix, tfxVec3(local_position_x, local_position_y, local_position_z) + handle);
				local_position_x = lerp_position.x + local_position_x * scale.x;
				local_position_y = lerp_position.y + local_position_y * scale.y;
				local_position_z = lerp_position.z + local_position_z * scale.z;
			}

			tween += entry->qty_step_size;
		}

	}

	void SpawnParticleWeight(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		random.AlterSeed(20 + pm.emitters.seed_index[emitter_index]);
		const float weight = pm.emitters.weight[emitter_index];
		const float weight_variation = pm.emitters.weight_variation[emitter_index];
		tfxLibrary *library = pm.library;

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &base_weight = entry->particle_data->base_weight[index];

			//----Weight
			if (weight) {
				base_weight = weight;
				if (weight_variation > 0) {
					base_weight += random.Range(-weight_variation, weight_variation);
				}
			}
			else {
				base_weight = 0;
			}
		}

	}

	void SpawnParticleVelocity(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		random.AlterSeed(21 + pm.emitters.seed_index[emitter_index]);
		tfxLibrary *library = pm.library;
		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[emitter_index];
		const float velocity = pm.emitters.velocity[emitter_index];
		const float velocity_variation = pm.emitters.velocity_variation[emitter_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &base_velocity = entry->particle_data->base_velocity[index];

			//----Velocity
			base_velocity = velocity + random.Range(-velocity_variation, velocity_variation);
		}

	}

	void SpawnParticleRoll(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		random.AlterSeed(22 + pm.emitters.seed_index[emitter_index]);
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const tfxAngleSettingFlags angle_settings = properties.angle_settings[property_index];
		const float angle_roll_offset = properties.angle_offsets[property_index].roll;

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &roll = entry->particle_data->local_rotations_z[index];

			roll = 0;
			if (angle_settings & tfxAngleSettingFlags_random_roll) {
				roll = random.Range(angle_roll_offset);
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
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(23 + pm.emitters.seed_index[emitter_index]);
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const float splatter = pm.emitters.splatter[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 scale = pm.emitters.scale[emitter_index];

		if (splatter) {
			for (int i = 0; i != entry->amount_to_spawn; ++i) {
				tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
				float &local_position_x = entry->particle_data->position_x[index];
				float &local_position_y = entry->particle_data->position_y[index];
				float &local_position_z = entry->particle_data->position_z[index];

				//----Splatter
				float splattertemp = splatter;
				float splatx = random.Range(-splatter, splatter);
				float splaty = random.Range(-splatter, splatter);

				while (GetDistance(0, 0, splatx, splaty) >= splattertemp && splattertemp > 0) {
					splatx = random.Range(-splatter, splatter);
					splaty = random.Range(-splatter, splatter);
				}

				if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
					local_position_x += splatx * scale.x;
					local_position_y += splaty * scale.y;
				}
				else {
					local_position_x += splatx;
					local_position_y += splaty;
				}
			}
		}

		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[emitter_index];
		tfxLibrary *library = pm.library;
		const float first_velocity_value = library->emitter_attributes[emitter_attributes].overtime.velocity.GetFirstValue();
		const float first_weight_value = library->emitter_attributes[emitter_attributes].overtime.weight.GetFirstValue();
		const tfxAngleSettingFlags angle_settings = properties.angle_settings[property_index];
		const float angle_roll_offset = properties.angle_offsets[property_index].roll;
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxVec3 emitter_world_rotations = pm.emitters.world_rotations[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxEmissionType emission_type = properties.emission_type[property_index];
		const tfxVec2 emitter_size = pm.emitters.emitter_size[emitter_index].xy();
		const tfxU32 layer = properties.layer[property_index];

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			const float base_weight = entry->particle_data->base_weight[index];
			float &roll = entry->particle_data->local_rotations_z[index];
			float &local_position_x = entry->particle_data->position_x[index];
			float &local_position_y = entry->particle_data->position_y[index];
			float &captured_position_x = entry->particle_data->captured_position_x[index];
			float &captured_position_y = entry->particle_data->captured_position_y[index];
			float &direction = entry->particle_data->local_rotations_x[index];
			float &base_velocity = entry->particle_data->base_velocity[index];

			float micro_time = tfxUPDATE_TIME * (1.f - tween);

			tfxVec2 sprite_transform_position;
			float sprite_transform_rotation;

			bool line = property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == tfxLine;

			if (!line && !(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				TransformParticlePosition(local_position_x, local_position_y, roll, sprite_transform_position, sprite_transform_rotation, emitter_world_rotations, matrix, handle, scale, emitter_world_position);
				captured_position_x = sprite_transform_position.x;
				captured_position_y = sprite_transform_position.y;
			}

			direction = 0;
			if (!line) {
				direction = GetEmissionDirection2d(pm, library, random, property_index, emitter_index, tfxVec2(local_position_x, local_position_y), sprite_transform_position, emitter_size) + library->emitter_attributes[emitter_attributes].overtime.direction.GetFirstValue();
			}

			tfxVec2 velocity_normal;
			velocity_normal.x = sinf(direction);
			velocity_normal.y = -cosf(direction);
			float weight_acceleration = base_weight * first_weight_value;
			//----Velocity Changes
			tfxVec2 current_velocity = velocity_normal * (base_velocity * first_velocity_value);
			current_velocity.y += weight_acceleration;
			local_position_x += current_velocity.x * micro_time;
			local_position_y += current_velocity.y * micro_time;
			if (line || property_flags & tfxEmitterPropertyFlags_relative_position) {
				tfxVec2 rotatevec = mmTransformVector(matrix, tfxVec2(local_position_x, local_position_y) + handle.xy());
				captured_position_x = emitter_captured_position.x + rotatevec.x * scale.x;
				captured_position_y = emitter_captured_position.y + rotatevec.y * scale.y;
			}

			if ((angle_settings & tfxAngleSettingFlags_align_roll || angle_settings & tfxAngleSettingFlags_align_with_emission) && !line) {
				roll = GetVectorAngle(velocity_normal.x, velocity_normal.y) + angle_roll_offset;
			}

			if (pm.flags & tfxEffectManagerFlags_ordered_by_age) {
				tfxDepthIndex depth_index;
				depth_index.particle_id = MakeParticleID(particles_index, index);
				depth_index.depth = 0.f;
				entry->particle_data->depth_index[index] = pm.PushDepthIndex(layer, depth_index);
			}

			tween += entry->qty_step_size;
		}
	}

	void SpawnParticleMicroUpdate3d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxSpawnWorkEntry *entry = static_cast<tfxSpawnWorkEntry*>(data);
		tfxRandom random = entry->pm->random;
		float tween = entry->tween;
		tfxU32 emitter_index = entry->emitter_index;
		tfxParticleManager &pm = *entry->pm;
		random.AlterSeed(24 + pm.emitters.seed_index[emitter_index]);
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxEmitterPropertiesSoA &properties = *entry->properties;
		const float splatter = pm.emitters.splatter[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxVec3 scale = pm.emitters.scale[emitter_index];

		if (splatter) {
			for (int i = 0; i != entry->amount_to_spawn; ++i) {
				tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
				float &local_position_x = entry->particle_data->position_x[index];
				float &local_position_y = entry->particle_data->position_y[index];
				float &local_position_z = entry->particle_data->position_z[index];

				//----Splatter
				if (splatter) {
					float splatx = random.Range(-splatter, splatter);
					float splaty = random.Range(-splatter, splatter);
					float splatz = random.Range(-splatter, splatter);

					while (std::powf(splatx / splatter, 2.f) + std::powf(splaty / splatter, 2.f) + std::powf(splatz / splatter, 2.f) > 1.f) {
						splatx = random.Range(-splatter, splatter);
						splaty = random.Range(-splatter, splatter);
						splatz = random.Range(-splatter, splatter);
					}

					if (!(property_flags & tfxEmitterPropertyFlags_relative_position)) {
						local_position_x += splatx * scale.x;
						local_position_y += splaty * scale.y;
						local_position_z += splatz * scale.z;
					}
					else {
						local_position_x += splatx;
						local_position_y += splaty;
						local_position_z += splatz;
					}
				}

			}
		}

		const tfxU32 emitter_attributes = pm.emitters.emitter_attributes[emitter_index];
		tfxLibrary *library = pm.library;
		const float first_velocity_value = library->emitter_attributes[emitter_attributes].overtime.velocity.GetFirstValue();
		const float first_weight_value = library->emitter_attributes[emitter_attributes].overtime.weight.GetFirstValue();
		const tfxVec3 emitter_world_position = pm.emitters.world_position[emitter_index];
		const tfxVec3 emitter_captured_position = pm.emitters.captured_position[emitter_index];
		const tfxVec3 emitter_size = pm.emitters.emitter_size[emitter_index];
		const tfxVec3 handle = pm.emitters.handle[emitter_index];
		const tfxMatrix4 matrix = pm.emitters.matrix[emitter_index];
		const tfxEmissionType emission_type = properties.emission_type[property_index];
		const tfxVectorAlignType vector_align_type = properties.vector_align_type[property_index];
		const bool line = property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == tfxLine;
		const float velocity_adjuster = pm.emitters.velocity_adjuster[emitter_index];
		const float frame = pm.emitters.frame[emitter_index];
		const tfxU32 layer = properties.layer[property_index];

		//Micro Update
		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
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

			tfxVec3 world_position;
			if (!line && !(property_flags & tfxEmitterPropertyFlags_relative_position)) {
				if (!(property_flags & tfxEmitterPropertyFlags_relative_position) && !(property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
					world_position.x = local_position_x;
					world_position.y = local_position_y;
					world_position.z = local_position_z;
				}
				else {
					tfxVec4 rotatevec = mmTransformVector(matrix, tfxVec3(local_position_x, local_position_y, local_position_z) + handle);
					world_position = emitter_world_position + rotatevec.xyz() * scale;
				}
				captured_position_x = world_position.x;
				captured_position_y = world_position.y;
				captured_position_z = world_position.z;
			}

			//----Velocity
			float emission_pitch = lookup_callback(library->emitter_attributes[emitter_attributes].properties.emission_pitch, frame);
			float emission_yaw = lookup_callback(library->emitter_attributes[emitter_attributes].properties.emission_yaw, frame);

			tfxVec3 velocity_normal;
			if (!(property_flags & tfxEmitterPropertyFlags_edge_traversal) || emission_type != tfxLine) {
				velocity_normal = GetEmissionDirection3d(pm, library, random, property_index, emitter_index, emission_pitch, emission_yaw, tfxVec3(local_position_x, local_position_y, local_position_z), world_position, emitter_size);
				velocity_normal_packed = Pack10bitUnsigned(velocity_normal);
			}
			else if (property_flags & tfxEmitterPropertyFlags_edge_traversal && emission_type == tfxLine) {
				velocity_normal_packed = tfxPACKED_Y_NORMAL_3D;
			}
			float velocity_scale = first_velocity_value * velocity_adjuster * base_velocity;

			//Do a micro update
			//A bit hacky but the epsilon after tween just ensures that theres a guaranteed small difference between captured/world positions so that
			//the alignment on the first frame can be calculated
			float micro_time = tfxUPDATE_TIME * (1.f - tween + 0.001f);
			float weight_acceleration = base_weight * first_weight_value;
			//----Velocity Changes
			tfxVec3 current_velocity = tfxVec3(velocity_normal.x, velocity_normal.y, velocity_normal.z) * base_velocity * first_velocity_value;
			current_velocity.y -= weight_acceleration;
			current_velocity *= micro_time;
			local_position_x += current_velocity.x;
			local_position_y += current_velocity.y;
			local_position_z += current_velocity.z;
			if (line || property_flags & tfxEmitterPropertyFlags_relative_position) {

				if (!(property_flags & tfxEmitterPropertyFlags_relative_position) && !(property_flags & tfxEmitterPropertyFlags_edge_traversal)) {
					world_position.x = local_position_x;
					world_position.y = local_position_y;
					world_position.z = local_position_z;
				}
				else {
					tfxVec4 rotatevec = mmTransformVector(matrix, tfxVec3(local_position_x, local_position_y, local_position_z) + handle);
					captured_position_x = world_position.x = emitter_captured_position.x + rotatevec.x * scale.x;
					captured_position_y = world_position.y = emitter_captured_position.y + rotatevec.y * scale.y;
					captured_position_z = world_position.z = emitter_captured_position.z + rotatevec.z * scale.z;
				}
			}
			else {
				world_position += current_velocity;
				captured_position_x = world_position.x;
				captured_position_y = world_position.y;
				captured_position_z = world_position.z;
			}
			if (pm.flags & tfxEffectManagerFlags_order_by_depth) {
				tfxDepthIndex depth_index;
				depth_index.particle_id = MakeParticleID(particles_index, index);
				depth_index.depth = LengthVec3NoSqR(world_position - pm.camera_position);
				entry->particle_data->depth_index[index] = pm.PushDepthIndex(layer, depth_index);
			}
			else if (pm.flags & tfxEffectManagerFlags_ordered_by_age) {
				tfxDepthIndex depth_index;
				depth_index.particle_id = MakeParticleID(particles_index, index);
				depth_index.depth = 0.f;
				entry->particle_data->depth_index[index] = pm.PushDepthIndex(layer, depth_index);
			}
			tween += entry->qty_step_size;
		}
	}

	void UpdateEmitterState(tfxParticleManager &pm, tfxU32 index, tfxU32 parent_index, const tfxParentSpawnControls &parent_spawn_controls, tfxSpawnWorkEntry *entry) {
		tfxPROFILE;
		tfxLibrary *library = pm.library;
		tfxEmitterPropertiesSoA &properties = library->emitter_properties;

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
			spawn_quantity += amount_variation > 0.f ? pm.random.Range(1.f, amount_variation) : 0.f;

			spawn_quantity *= lookup_callback(library->global_graphs[global_attributes].amount, frame);
		}
		else {
			spawn_quantity = (float)properties.spawn_amount[property_index];
			spawn_quantity *= lookup_callback(library->global_graphs[global_attributes].amount, frame);
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
		tfxLibrary *library = pm.effects.library[index];
		tfxVec3 &translation = pm.effects.translation[index];
		tfxVec3 &local_rotations = pm.effects.local_rotations[index];
		tfxVec3 &scale = pm.effects.scale[index];
		tfxVec3 &emitter_size = pm.effects.emitter_size[index];
		float &overal_scale = pm.effects.overal_scale[index];
		tfxEffectStateFlags &state_flags = pm.effects.state_flags[index];
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
		if (!(state_flags & tfxEffectStateFlags_override_size_multiplier)) {
			emitter_size.x = lookup_callback(library->global_graphs[global_attributes].emitter_width, frame);
			emitter_size.y = lookup_callback(library->global_graphs[global_attributes].emitter_height, frame);
			emitter_size.z = lookup_callback(library->global_graphs[global_attributes].emitter_depth, frame);
		}
		//We don't want to scale twice when the sub effect is transformed, so the values here are set to 1. That means that the root effect will only control the global scale.
		overal_scale = state_flags & tfxEffectStateFlags_override_overal_scale ? overal_scale : lookup_callback(library->global_graphs[global_attributes].overal_scale, frame);
		if (pm.effects.parent_particle_index[index] == tfxINVALID) {
			scale.x = overal_scale;
			scale.y = overal_scale;
			scale.z = overal_scale;
			if (!(state_flags & tfxEffectStateFlags_override_orientiation)) {
				local_rotations.roll = LookupPrecise(library->transform_attributes[transform_attributes].roll, age);
				local_rotations.pitch = LookupPrecise(library->transform_attributes[transform_attributes].pitch, age);
				local_rotations.yaw = LookupPrecise(library->transform_attributes[transform_attributes].yaw, age);
			}
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
		translation.x = LookupPrecise(library->transform_attributes[transform_attributes].translation_x, age);
		translation.y = LookupPrecise(library->transform_attributes[transform_attributes].translation_y, age);
		translation.z = LookupPrecise(library->transform_attributes[transform_attributes].translation_z, age);

		if (pm.effects.update_callback[index])
			pm.effects.update_callback[index](&pm, index);

	}

	void ControlParticleAge(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxParticleAgeWorkEntry *work_entry = static_cast<tfxParticleAgeWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		tfxParticleManager &pm = *work_entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		tfxParticleSoA &bank = pm.particle_arrays[particles_index];
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxU32 property_flags = pm.emitters.property_flags[emitter_index];
		const tfxWideInt single_shot_limit = tfxWideSetSinglei(work_entry->properties->single_shot_limit[property_index]);
		const tfxU32 layer = work_entry->properties->layer[property_index];

		const tfxWideInt remove_flag = tfxWideSetSinglei(tfxParticleFlags_remove);
		const tfxWideInt remove = tfxWideSetSinglei(pm.emitters.state_flags[emitter_index] & tfxParticleFlags_remove);
		const tfxWideInt single = tfxWideGreateri(tfxWideSetSinglei(property_flags & tfxEmitterPropertyFlags_single), tfxWideSetZeroi());
		const tfxWideInt not_single = tfxWideXOri(single, tfxWideSetSinglei(-1));
		tfxWideInt state_flags_no_spawning = tfxWideGreateri(tfxWideOri(tfxWideSetSinglei(pm.emitters.state_flags[emitter_index] & tfxEmitterStateFlags_stop_spawning), tfxWideSetSinglei(work_entry->pm->flags & tfxEffectManagerFlags_disable_spawning)), tfxWideSetZeroi());
		const tfxWideInt xor_state_flags_no_spawning = tfxWideXOri(state_flags_no_spawning, tfxWideSetSinglei(-1));

		for (int i = 0; i != work_entry->wide_end_index; i += tfxDataWidth) {
			tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;

			const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
			tfxWideFloat age = tfxWideLoad(&bank.age[index]);
			tfxWideInt single_loop_count = tfxWideLoadi((tfxWideInt*)&bank.single_loop_count[index]);
			tfxWideInt flags = tfxWideLoadi((tfxWideInt*)&bank.flags[index]);
			age = tfxWideAdd(age, tfxFRAME_LENGTH_WIDE);

			tfxWideInt expired = tfxWideCasti(tfxWideGreaterEqual(age, max_age));
			single_loop_count = tfxWideAddi(single_loop_count, tfxWideAndi(tfxWideSetSinglei(1), expired));
			tfxWideInt loop_limit = tfxWideEqualsi(single_loop_count, single_shot_limit);
			tfxWideInt loop_age = tfxWideXOri(tfxWideAndi(tfxWideAndi(single, expired), xor_state_flags_no_spawning), tfxWideSetSinglei(-1));
			age = tfxWideAnd(age, tfxWideCast(loop_age));
			flags = tfxWideOri(flags, tfxWideAndi(remove_flag, tfxWideGreateri(remove, tfxWideSetSinglei(0))));
			flags = tfxWideOri(flags, tfxWideAndi(remove_flag, tfxWideAndi(not_single, expired)));
			flags = tfxWideOri(flags, tfxWideAndi(remove_flag, tfxWideAndi(tfxWideOri(tfxWideAndi(single, loop_limit), state_flags_no_spawning), expired)));

			tfxWideStore(&bank.age[index], age);
			tfxWideStorei((tfxWideInt*)&bank.flags[index], flags);
			tfxWideStorei((tfxWideInt*)&bank.single_loop_count[index], single_loop_count);
		}


		tfxU32 offset = 0;
		for (int i = work_entry->start_index; i >= 0; --i) {
			const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i);
			tfxParticleFlags &flags = bank.flags[index];
			if (index == 16)
				int d = 0;
			if (flags & tfxParticleFlags_remove) {
				offset++;
				if (flags & tfxParticleFlags_has_sub_effects) {
					pm.FreeParticleIndex(bank.particle_index[index]);
				}
				if (!(pm.flags & tfxEffectManagerFlags_unordered)) {
					tfxU32 temp_depth = bank.depth_index[index];
					pm.depth_indexes[layer][pm.current_depth_index_buffer][bank.depth_index[index]].particle_id = tfxINVALID;
				}
			}
			else if (offset > 0) {
				tfxU32 next_index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i + offset);
				if (flags & tfxParticleFlags_has_sub_effects) {
					pm.particle_indexes[bank.particle_index[index]] = MakeParticleID(particles_index, next_index);
				}

				if (pm.flags & tfxEffectManagerFlags_order_by_depth) {
					pm.depth_indexes[layer][pm.current_depth_index_buffer][bank.depth_index[index]].particle_id = MakeParticleID(particles_index, next_index);
				}
				else if (pm.flags & tfxEffectManagerFlags_ordered_by_age) {
					pm.depth_indexes[layer][pm.current_depth_index_buffer][bank.depth_index[index]].particle_id = MakeParticleID(particles_index, next_index);
					pm.depth_indexes[layer][pm.current_depth_index_buffer][bank.depth_index[index]].depth = bank.age[index];
				}

				bank.parent_index[next_index] = bank.parent_index[index];
				bank.sprite_index[next_index] = bank.sprite_index[index];
				bank.depth_index[next_index] = bank.depth_index[index];
				bank.particle_index[next_index] = bank.particle_index[index];
				bank.flags[next_index] = bank.flags[index];
				bank.age[next_index] = bank.age[index];
				bank.max_age[next_index] = bank.max_age[index];
				bank.position_x[next_index] = bank.position_x[index];
				bank.position_y[next_index] = bank.position_y[index];
				bank.position_z[next_index] = bank.position_z[index];
				bank.captured_position_x[next_index] = bank.captured_position_x[index];
				bank.captured_position_y[next_index] = bank.captured_position_y[index];
				bank.captured_position_z[next_index] = bank.captured_position_z[index];
				bank.local_rotations_x[next_index] = bank.local_rotations_x[index];
				bank.local_rotations_y[next_index] = bank.local_rotations_y[index];
				bank.local_rotations_z[next_index] = bank.local_rotations_z[index];
				bank.velocity_normal[next_index] = bank.velocity_normal[index];
				bank.base_weight[next_index] = bank.base_weight[index];
				bank.base_velocity[next_index] = bank.base_velocity[index];
				bank.base_spin[next_index] = bank.base_spin[index];
				bank.noise_offset[next_index] = bank.noise_offset[index];
				bank.noise_resolution[next_index] = bank.noise_resolution[index];
				bank.color[next_index] = bank.color[index];
				bank.image_frame[next_index] = bank.image_frame[index];
				bank.base_size_x[next_index] = bank.base_size_x[index];
				bank.base_size_y[next_index] = bank.base_size_y[index];
				bank.single_loop_count[next_index] = bank.single_loop_count[index];
			}
		}

		if (offset) {
			Bump(&work_entry->pm->particle_array_buffers[particles_index], offset);
		}

	}

	void ControlParticlePosition2d(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		tfxParticleManager &pm = *work_entry->pm;
		const tfxU32 particles_index = pm.emitters.particles_index[emitter_index];
		tfxParticleSoA &bank = work_entry->pm->particle_arrays[particles_index];

		const tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const tfxMatrix4 &matrix = pm.emitters.matrix[emitter_index];
		const tfxWideFloat overal_scale_wide = tfxWideSetSingle(pm.emitters.overal_scale[emitter_index]);

		tfxU32 running_sprite_index = work_entry->sprites_index;
		tfxSpriteSoA &sprites = *work_entry->sprites;

		tfxWideFloat max_life = tfxWideSetSingle(work_entry->graphs->velocity.lookup.life);
		const tfxWideInt velocity_turbulance_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity_turbulance.lookup.last_frame);
		const tfxWideInt noise_resolution_last_frame = tfxWideSetSinglei(work_entry->graphs->noise_resolution.lookup.last_frame);

		//Noise
		const tfxWideInt velocity_last_frame = tfxWideSetSinglei(work_entry->graphs->velocity.lookup.last_frame);
		const tfxWideInt spin_last_frame = tfxWideSetSinglei(work_entry->graphs->spin.lookup.last_frame);
		const tfxWideFloat velocity_adjuster = tfxWideSetSingle(pm.emitters.velocity_adjuster[emitter_index]);
		const tfxWideInt weight_last_frame = tfxWideSetSinglei(work_entry->graphs->weight.lookup.last_frame);
		const tfxWideInt direction_last_frame = tfxWideSetSinglei(work_entry->graphs->direction.lookup.last_frame);
		const tfxWideInt stretch_last_frame = tfxWideSetSinglei(work_entry->graphs->stretch.lookup.last_frame);
		const tfxWideFloat stretch = tfxWideSetSingle(pm.emitters.stretch[emitter_index]);
		tfxWideArray p_stretch;

		tfxU32 start_diff = work_entry->start_diff;

		const float eps = 0.0001f;
		const float eps2 = 0.0002f;

		for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
			tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;

			const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
			const tfxWideFloat age = tfxWideLoad(&bank.age[index]);
			_ReadBarrier();
			tfxWideFloat life = tfxWideDiv(age, max_age);
			life = tfxWideMul(life, max_life);
			life = tfxWideDiv(life, tfxLOOKUP_FREQUENCY_OVERTIME_WIDE);

			const tfxWideFloat base_velocity = tfxWideLoad(&bank.base_velocity[index]);
			const tfxWideFloat base_spin = tfxWideLoad(&bank.base_spin[index]);
			const tfxWideFloat base_weight = tfxWideLoad(&bank.base_weight[index]);
			tfxWideFloat local_position_x = tfxWideLoad(&bank.position_x[index]);
			tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
			tfxWideFloat angle = tfxWideLoad(&bank.local_rotations_x[index]);
			tfxWideArray roll;
			roll.m = tfxWideLoad(&bank.local_rotations_z[index]);

			tfxWideArrayi lookup_frame;
			tfxWideArray noise_x;
			tfxWideArray noise_y;

			if (emitter_flags & tfxEmitterStateFlags_has_noise) {
				const tfxWideFloat noise_resolution = tfxWideLoad(&bank.noise_resolution[index]);
				const tfxWideFloat noise_offset = tfxWideLoad(&bank.noise_offset[index]);

				lookup_frame.m = tfxWideMini(tfxWideConverti(life), noise_resolution_last_frame);
				const tfxWideFloat lookup_noise_resolution = tfxWideMul(tfxWideLookupSet(work_entry->graphs->noise_resolution.lookup.values, lookup_frame), noise_resolution);
				lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_turbulance_last_frame);
				const tfxWideFloat lookup_velocity_turbulance = tfxWideLookupSet(work_entry->graphs->velocity_turbulance.lookup.values, lookup_frame);

				tfxWideArray x, y;
				x.m = tfxWideAdd(tfxWideDiv(local_position_x, lookup_noise_resolution), noise_offset);
				y.m = tfxWideAdd(tfxWideDiv(local_position_y, lookup_noise_resolution), noise_offset);

				for (int n = 0; n != tfxDataWidth; ++n) {
					tfx128 x4 = _mm_set1_ps(x.a[n]);
					tfx128 y4 = _mm_set1_ps(y.a[n]);

					tfx128 xeps4 = _mm_set_ps(x.a[n] - eps, x.a[n] + eps, x.a[n], x.a[n]);
					tfx128 yeps4 = _mm_set_ps(y.a[n], y.a[n], y.a[n] - eps, y.a[n] + eps);

					tfx128Array sample = tfxNoise4(x4, yeps4);
					float a = (sample.a[0] - sample.a[1]) / eps2;
					float b = (sample.a[2] - sample.a[3]) / eps2;
					noise_x.a[n] = a - b;

					y.a[n] += 100.f;
					tfx128 yeps4r = _mm_set_ps(y.a[n] - eps, y.a[n] + eps, y.a[n], y.a[n]);
					sample = tfxNoise4(xeps4, y4);
					a = (sample.a[0] - sample.a[1]) / eps2;
					b = (sample.a[2] - sample.a[3]) / eps2;
					noise_y.a[n] = a - b;
				}

				noise_x.m = tfxWideMul(lookup_velocity_turbulance, noise_x.m);
				noise_y.m = tfxWideMul(lookup_velocity_turbulance, noise_y.m);
			}

			lookup_frame.m = tfxWideMini(tfxWideConverti(life), spin_last_frame);
			const tfxWideFloat lookup_spin = tfxWideMul(tfxWideLookupSet(work_entry->graphs->spin.lookup.values, lookup_frame), base_spin);
			lookup_frame.m = tfxWideMini(tfxWideConverti(life), velocity_last_frame);
			const tfxWideFloat lookup_velocity = tfxWideLookupSet(work_entry->graphs->velocity.lookup.values, lookup_frame);
			lookup_frame.m = tfxWideMini(tfxWideConverti(life), weight_last_frame);
			const tfxWideFloat lookup_weight = tfxWideLookupSet(work_entry->graphs->weight.lookup.values, lookup_frame);
			lookup_frame.m = tfxWideMini(tfxWideConverti(life), direction_last_frame);
			tfxWideArray lookup_direction;
			lookup_direction.m = tfxWideLookupSet(work_entry->graphs->direction.lookup.values, lookup_frame);
			lookup_direction.m = tfxWideAdd(lookup_direction.m, angle);
			lookup_frame.m = tfxWideMini(tfxWideConverti(life), stretch_last_frame);
			const tfxWideFloat lookup_stretch = tfxWideLookupSet(work_entry->graphs->stretch.lookup.values, lookup_frame);
			p_stretch.m = tfxWideMul(lookup_stretch, stretch);

			//----Velocity Changes
			tfxWideFloat velocity_scalar = tfxWideMul(base_velocity, lookup_velocity);
			tfxWideArray current_velocity_x;
			tfxWideArray current_velocity_y;
			for (int j = 0; j != tfxDataWidth; ++j) {
				current_velocity_x.a[j] = sinf(lookup_direction.a[j]);
				current_velocity_y.a[j] = -cosf(lookup_direction.a[j]);
			}
			current_velocity_x.m = tfxWideMul(current_velocity_x.m, velocity_scalar);
			current_velocity_y.m = tfxWideMul(current_velocity_y.m, velocity_scalar);
			if (emitter_flags & tfxEmitterStateFlags_has_noise) {
				current_velocity_x.m = tfxWideAdd(current_velocity_x.m, noise_x.m);
				current_velocity_y.m = tfxWideAdd(current_velocity_y.m, noise_y.m);
			}
			current_velocity_y.m = tfxWideAdd(current_velocity_y.m, tfxWideMul(lookup_weight, base_weight));
			current_velocity_x.m = tfxWideMul(tfxWideMul(current_velocity_x.m, tfxUPDATE_TIME_WIDE), velocity_adjuster);
			current_velocity_y.m = tfxWideMul(tfxWideMul(current_velocity_y.m, tfxUPDATE_TIME_WIDE), velocity_adjuster);

			//----Spin and angle Changes
			if (emitter_flags & tfxEmitterStateFlags_can_spin) {
				roll.m = tfxWideAdd(roll.m, tfxWideMul(lookup_spin, tfxUPDATE_TIME_WIDE));
			}

			//----Position
			local_position_x = tfxWideAdd(local_position_x, tfxWideMul(current_velocity_x.m, overal_scale_wide));
			local_position_y = tfxWideAdd(local_position_y, tfxWideMul(current_velocity_y.m, overal_scale_wide));

			tfxWideStore(&bank.position_x[index], local_position_x);
			tfxWideStore(&bank.position_y[index], local_position_y);
			tfxWideStore(&bank.local_rotations_z[index], roll.m);

			current_velocity_y.m = tfxWideAdd(current_velocity_y.m, tfxWideSetSingle(0.000001f));
			tfxWideFloat l = tfxWideMul(current_velocity_x.m, current_velocity_x.m);
			l = tfxWideAdd(l, tfxWideMul(current_velocity_y.m, current_velocity_y.m));
			l = tfxWideSqrt(l);
			p_stretch.m = tfxWideMul(p_stretch.m, l);
			current_velocity_x.m = tfxWideDiv(current_velocity_x.m, l);
			current_velocity_y.m = tfxWideDiv(current_velocity_y.m, l);

			if (property_flags & tfxEmitterPropertyFlags_relative_position) {
				mmWideTransformVector(matrix, current_velocity_x.m, current_velocity_y.m);
			}

			tfxWideArrayi packed;
			packed.m = PackWide16bit(current_velocity_x.m, current_velocity_y.m);

			tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
			if (!(pm.flags & tfxEffectManagerFlags_unordered)) {	//Predictable
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 sprite_depth_index = bank.depth_index[index + j];
					sprites.stretch[sprite_depth_index] = p_stretch.a[j];
					sprites.alignment[sprite_depth_index] = packed.a[j];
					running_sprite_index++;
				}
			}else {
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					sprites.stretch[running_sprite_index] = p_stretch.a[j];
					sprites.alignment[running_sprite_index++] = packed.a[j];
				}
			}
			start_diff = 0;
		}

		const tfxWideFloat emitter_size_y = tfxWideSetSingle(pm.emitters.emitter_size[emitter_index].y);
		const tfxWideInt emitter_flags_wide = tfxWideSetSinglei(emitter_flags);

		if (emitter_flags & tfxEmitterStateFlags_is_line_loop_or_kill) {
			//Todo: this should also update the captured position as well solving the issue of interpolating from the end back to the beginning
			if (emitter_flags & tfxEmitterStateFlags_kill) {
				for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
					tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;
					const tfxWideFloat offset_y = tfxWideMul(UnPackWide10bitY(tfxWideLoadi((tfxWideInt*)&bank.velocity_normal[index])), emitter_size_y);
					tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
					tfxWideInt flags = tfxWideLoadi((tfxWideInt*)&bank.flags[index]);

					//Lines - Reposition if the particle is travelling along a line
					tfxWideFloat length = tfxWideAbs(local_position_y);
					tfxWideInt remove_flags = tfxWideAndi(tfxWideSetSinglei(tfxParticleFlags_remove), tfxWideCasti(tfxWideGreater(length, emitter_size_y)));
					flags = tfxWideOri(flags, remove_flags);
					tfxWideStorei((tfxWideInt*)&bank.flags[index], flags);
				}
			}
			else {
				for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
					tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;
					const tfxWideFloat offset_y = tfxWideMul(UnPackWide10bitY(tfxWideLoadi((tfxWideInt*)&bank.velocity_normal[index])), emitter_size_y);
					tfxWideFloat local_position_y = tfxWideLoad(&bank.position_y[index]);
					tfxWideInt flags = tfxWideLoadi((tfxWideInt*)&bank.flags[index]);

					//Lines - Reposition if the particle is travelling along a line
					tfxWideFloat length = tfxWideAbs(local_position_y);
					tfxWideFloat at_end = tfxWideGreater(length, emitter_size_y);
					local_position_y = tfxWideSub(local_position_y, tfxWideAnd(at_end, offset_y));
					flags = tfxWideOri(flags, tfxWideAndi(tfxWideSetSinglei(tfxParticleFlags_capture_after_transform), tfxWideCasti(at_end)));
					tfxWideStorei((tfxWideInt*)&bank.flags[index], flags);
					tfxWideStore(&bank.position_y[index], local_position_y);
				}
			}
		}

		ControlParticleTransform2d(&pm.work_queue, data);
	}

	void ControlParticles2d(tfxParticleManager &pm, tfxU32 emitter_index, tfxControlWorkEntry &work_entry) {
		tfxPROFILE;

		tfxLibrary *library = pm.library;
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
		work_entry.sprite_buffer_end_index = work_entry.sprites_index + work_entry.end_index;
		work_entry.layer = properties.layer[property_index];
		work_entry.sprites = &pm.sprites[pm.current_sprite_buffer][work_entry.layer];

		if (amount_to_update > 0) {
			if (!(pm.flags & tfxEffectManagerFlags_single_threaded) && tfxNumberOfThreadsInAdditionToMain) {
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticlePosition2d);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleSize);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleColor);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleImageFrame);
			}
			else {
				ControlParticlePosition2d(&pm.work_queue, &work_entry);
				ControlParticleSize(&pm.work_queue, &work_entry);
				ControlParticleColor(&pm.work_queue, &work_entry);
				ControlParticleImageFrame(&pm.work_queue, &work_entry);
			}
		}

	}

	void ControlParticles3d(tfxParticleManager &pm, tfxU32 emitter_index, tfxControlWorkEntry &work_entry) {
		tfxPROFILE;

		tfxLibrary *library = pm.library;
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
		work_entry.sprite_buffer_end_index = work_entry.sprites_index + (work_entry.end_index - work_entry.start_index);
		work_entry.layer = properties.layer[property_index];
		work_entry.sprites = &pm.sprites[pm.current_sprite_buffer][work_entry.layer];

		if (amount_to_update > 0) {
			if (!(pm.flags & tfxEffectManagerFlags_single_threaded) && tfxNumberOfThreadsInAdditionToMain) {
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticlePosition3d);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleSize);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleColor);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleImageFrame);
			}
			else {
				ControlParticlePosition3d(&pm.work_queue, &work_entry);
				ControlParticleSize(&pm.work_queue, &work_entry);
				ControlParticleColor(&pm.work_queue, &work_entry);
				ControlParticleImageFrame(&pm.work_queue, &work_entry);
			}
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

	const tfxU32 tfxPROFILE_COUNT = __COUNTER__;
	tfxU32 tfxCurrentSnapshot = 0;
	tfxProfile tfxProfileArray[tfxPROFILE_COUNT];
	tfxMemoryTrackerLog tfxMEMORY_TRACKER;
	char tfxMEMORY_CONTEXT[64];
	tfxDataTypesDictionary tfxDataTypes;
	void *tfxDeferred_data_for_freeing[256];
	tfxU32 tfxDeferred_index = 0;
	HANDLE tfxThreadSemaphore;
	int tfxNumberOfThreadsInAdditionToMain;
	tfxQueueProcessor tfxThreadQueues;
	tfxMemoryArenaManager tfxSTACK_ALLOCATOR;
	tfxMemoryArenaManager tfxMT_STACK_ALLOCATOR;

	//Passing a max_threads value of 0 or 1 will make timeline fx run in single threaded mode. 2 or more will be multithreaded.
	//max_threads includes the main thread so for example if you set it to 4 then there will be the main thread plus an additional 3 threads.
	void InitialiseTimelineFX(int max_threads) {
		tfxSTACK_ALLOCATOR = CreateArenaManager(tfxSTACK_SIZE, 8);
		tfxMT_STACK_ALLOCATOR = CreateArenaManager(tfxMT_STACK_SIZE, 8);
		tfxNumberOfThreadsInAdditionToMain = tfxMin(max_threads - 1 < 0 ? 0 : max_threads - 1, (int)std::thread::hardware_concurrency() - 1);
		lookup_callback = LookupFast;
		lookup_overtime_callback = LookupFastOvertime;
		memset(tfxProfileArray, 0, tfxPROFILE_COUNT * sizeof(tfxProfile));
		tfxInitialiseThreads(&tfxThreadQueues);
	}

	tfxProfileTag::tfxProfileTag(tfxU32 id, const char *name) {
		profile = tfxProfileArray + id;
		profile->name = name;
		snapshot = profile->snapshots + tfxCurrentSnapshot;
		start_time = Microsecs();
		start_cycles = __rdtsc();
		AtomicAdd32(&snapshot->hit_count, 1);
	}

	void InitParticleManagerFor3d(tfxParticleManager *pm, tfxLibrary *library, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfxParticleManagerModes mode, bool double_buffered_sprites, bool dynamic_allocation, tfxU32 mt_batch_size) {
		assert(pm->flags == 0);		//You must use a particle manager that has not been initialised already. You can call reconfigure if you want to re-initialise a particle manager
		pm->SetLibrary(library);
		pm->InitFor3d(library, layer_max_values, effects_limit, mode, double_buffered_sprites, dynamic_allocation, mt_batch_size);
	}

	void InitParticleManagerFor2d(tfxParticleManager *pm, tfxLibrary *library, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfxParticleManagerModes mode, bool double_buffered_sprites, bool dynamic_allocation, tfxU32 mt_batch_size) {
		assert(pm->flags == 0);		//You must use a particle manager that has not been initialised already. You can call reconfigure if you want to re-initialise a particle manager
		pm->SetLibrary(library);
		pm->InitFor2d(library, layer_max_values, effects_limit, mode, double_buffered_sprites, dynamic_allocation, mt_batch_size);
	}


	tfxU32 AddEffectToParticleManager(tfxParticleManager *pm, tfxEffectTemplate &effect) {
		return pm->AddEffect(effect.effect, pm->current_ebuff, 0, false, 0.f);
	}

	void SetEffectPosition(tfxParticleManager *pm, tfxEffectID effect_index, float x, float y) {
		tfxVec2 position(x, y);
		pm->effects.local_position[effect_index] = position;
	}

	void SetEffectPosition(tfxParticleManager *pm, tfxEffectID effect_index, tfxVec2 position) {
		pm->effects.local_position[effect_index] = position;
	}

	void SetEffectPosition(tfxParticleManager *pm, tfxEffectID effect_index, float x, float y, float z) {
		tfxVec3 position(x, y, z);
		pm->effects.local_position[effect_index] = position;
	}

	void SetEffectPosition(tfxParticleManager *pm, tfxEffectID effect_index, tfxVec3 position) {
		pm->effects.local_position[effect_index] = position;
	}

	void MoveEffect(tfxParticleManager *pm, tfxEffectID effect_index, tfxVec3 amount) {
		pm->effects.local_position[effect_index] += amount;
	}

	void MoveEffect(tfxParticleManager *pm, tfxEffectID effect_index, float x, float y, float z) {
		pm->effects.local_position[effect_index] += {x, y, z};
	}

	tfxAPI tfxVec3 GetEffectPosition(tfxParticleManager *pm, tfxEffectID effect_index) {
		return pm->effects.local_position[effect_index];
	}

	void SetEffectRotation(tfxParticleManager *pm, tfxEffectID effect_index, float rotation) {
		pm->effects.local_rotations[effect_index].roll = rotation;
		pm->effects.state_flags[effect_index] |= tfxEffectStateFlags_override_orientiation;
	}

	void SetEffectRoll(tfxParticleManager *pm, tfxEffectID effect_index, float roll) {
		pm->effects.local_rotations[effect_index].roll = roll;
		pm->effects.state_flags[effect_index] |= tfxEffectStateFlags_override_orientiation;
	}

	void SetEffectPitch(tfxParticleManager *pm, tfxEffectID effect_index, float pitch) {
		pm->effects.local_rotations[effect_index].pitch = pitch;
		pm->effects.state_flags[effect_index] |= tfxEffectStateFlags_override_orientiation;
	}

	void SetEffectYaw(tfxParticleManager *pm, tfxEffectID effect_index, float pitch) {
		pm->effects.local_rotations[effect_index].pitch = pitch;
		pm->effects.state_flags[effect_index] |= tfxEffectStateFlags_override_orientiation;
	}

	void SetEffectWidthMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float width) {
		pm->effects.emitter_size[effect_index].x = width;
		pm->effects.state_flags[effect_index] |= tfxEffectStateFlags_override_size_multiplier;
	}

	void SetEffectHeightMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float height) {
		pm->effects.emitter_size[effect_index].y = height;
		pm->effects.state_flags[effect_index] |= tfxEffectStateFlags_override_size_multiplier;
	}

	void SetEffectDepthMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float depth) {
		pm->effects.emitter_size[effect_index].z = depth;
		pm->effects.state_flags[effect_index] |= tfxEffectStateFlags_override_size_multiplier;
	}

	void SetEffectLifeMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float life) {
		pm->effects.spawn_controls[effect_index].life = life;
	}

	void SetEffectParticleWidthMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float width) {
		pm->effects.spawn_controls[effect_index].size_x = width;
	}

	void SetEffectParticleHeightMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float height) {
		pm->effects.spawn_controls[effect_index].size_y = height;
	}

	void SetEffectVelocityMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float velocity) {
		pm->effects.spawn_controls[effect_index].velocity = velocity;
	}

	void SetEffectSpinMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float spin) {
		pm->effects.spawn_controls[effect_index].spin = spin;
	}

	void SetEffectIntensityMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float intensity) {
		pm->effects.spawn_controls[effect_index].intensity = intensity;
	}

	void SetEffectSplatterMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float splatter) {
		pm->effects.spawn_controls[effect_index].splatter = splatter;
	}

	void SetEffectWeightMultiplier(tfxParticleManager *pm, tfxEffectID effect_index, float weight) {
		pm->effects.spawn_controls[effect_index].weight = weight;
	}

	void SetEffectOveralScale(tfxParticleManager *pm, tfxEffectID effect_index, float overal_scale) {
		pm->effects.overal_scale[effect_index] = overal_scale;
		pm->effects.state_flags[effect_index] |= tfxEffectStateFlags_override_overal_scale;
	}

	void SetEffectBaseNoiseOffset(tfxParticleManager *pm, tfxEffectID effect_index, float noise_offset) {
		pm->effects.noise_base_offset[effect_index] = noise_offset;
	}

}