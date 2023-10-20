#include "timelinefx.h"

namespace tfx {

#ifdef _WIN32
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
#else
	FILE *tfx__open_file(const char *file_name, const char *mode) {
		return fopen(file_name, mode);
	}
#endif

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
		FILE *file = tfx__open_file(file_name, "wb");
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

	void tfxStr::Setf(const char *format, ...) {
		Clear();
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

	bool ValidatePackage(tfxPackage *package) {
		if (package->header.magic_number != tfxMAGIC_NUMBER) return false;			//Package hasn't been initialised

		FILE *file = tfx__open_file(package->file_path.c_str(), "rb");
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

		if (length != package->file_size) return false;							//The file on disk is no longer the same size as the package file size since it was loaded

		//Everything seems ok
		fclose(file);
		return true;
	}

	void tfxEntryInfo::FreeData() {
		data.FreeAll();
	}

	tfxPackage::~tfxPackage() {
		FreePackage(this);
	}

	tfxEntryInfo *GetPackageFile(tfxPackage *package, const char *name) {
		if (!package->inventory.entries.ValidName(name)) {
			return nullptr;									//File not found in inventory
		}
		assert(ValidatePackage(package));						//The file on disk has changed since the package was loaded! Maybe this should return null instead?
		tfxEntryInfo *entry = &package->inventory.entries.At(name);
		if (entry->data.Size() != entry->file_size) {
			//FILE *file = tfx__open_file(file_path.c_str(), "rb");
			FILE *file = tfx__open_file(package->file_path.c_str(), "rb");
			assert(file);		//couldn't open the file!
			_fseeki64(file, entry->offset_from_start_of_file, SEEK_SET);
			entry->data.Resize(entry->file_size);
			fread(entry->data.data, 1, entry->file_size, file);
			fclose(file);
		}
		return entry;
	}

	bool FileExists(tfxPackage *package, const char *file_name) {
		if (package->inventory.entries.ValidName(file_name)) {
			return true;
		}
		return false;
	}

	void AddEntryToPackage(tfxPackage *package, tfxEntryInfo file) {
		package->inventory.entries.Insert(file.file_name, file);
		package->inventory.entry_count++;
	}

	void AddFileToPackage(tfxPackage *package, const char *file_name, tfxStream *data) {
		tfxEntryInfo entry;
		entry.file_name = file_name;
		entry.data = *data;
		entry.file_size = data->size;

		package->inventory.entries.Insert(entry.file_name, entry);
		package->inventory.entry_count++;
	}

	void FreePackage(tfxPackage *package) {
		for (auto &entry : package->inventory.entries.data) {
			entry.data.FreeAll();
		}
		package->inventory.entries.data.free_all();
		package->inventory.entries.map.free_all();
		package->file_data.FreeAll();
	}

	// Reads the whole file on disk into memory and returns the pointer
	tfxStream ReadEntireFile(const char *file_name, bool terminate) {
		tfxStream buffer;
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

	bool SavePackageDisk(tfxPackage &package) {
		if (!package.file_path.Length()) return false;											//Package must have a file path
		if (package.header.magic_number != tfxMAGIC_NUMBER) return false;						//Header of package must contain correct magic number. Use CreatePackage to correctly initialise a package.
		if (package.inventory.magic_number != tfxMAGIC_NUMBER_INVENTORY) return false;			//Inventory of package must contain correct magic number

		FILE * file = tfx__open_file(package.file_path.c_str(), "wb");
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

	void UpdateEffectMaxLife(tfxEffectEmitter *effect) {
		tfxEffectEmitterInfo *info = GetEffectInfo(effect);
		info->max_life = GetMaxLife(*effect);
		GetEffectGraphByType(effect, tfxOvertime_red)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_green)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_blue)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_blendfactor)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_intensity)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_velocity)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_width)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_height)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_weight)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_spin)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_stretch)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_spin)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_velocity_turbulance)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_direction_turbulance)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_velocity_adjuster)->lookup.life = info->max_life;
		GetEffectGraphByType(effect, tfxOvertime_direction)->lookup.life = info->max_life;
	}

	bool IsFiniteEffect(tfxEffectEmitter *effect) {
		for (auto &e : GetEffectInfo(effect)->sub_effectors) {
			float qty = GetGraphLastValue(&e.library->emitter_attributes[e.emitter_attributes].base.amount) + GetGraphLastValue(&e.library->emitter_attributes[e.emitter_attributes].variation.amount);
			if (!(e.property_flags & tfxEmitterPropertyFlags_single) && qty > 0)
				return false;
			else if (e.property_flags & tfxEmitterPropertyFlags_single && GetEffectProperties(&e)->single_shot_limit[e.property_index] == 0)
				return false;

		}
		return true;
	}

	void FlagEffectAs3D(tfxEffectEmitter *effect, bool flag) {
		if (flag) {
			effect->property_flags |= tfxEmitterPropertyFlags_is_3d;
		} else {
			effect->property_flags &= ~tfxEmitterPropertyFlags_is_3d;
		}
		for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
			FlagEffectAs3D(&sub, flag);
		}
	}

	bool Is3DEffect(tfxEffectEmitter *effect) {
		return effect->property_flags & tfxEmitterPropertyFlags_is_3d;
	}

	tfxParticleManagerModes GetRequiredParticleManagerMode(tfxEffectEmitter *effect) {
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
			tfxParticleManagerModes result = tfxParticleManagerMode_unordered;
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

	tfxPreviewCameraSettings *GetEffectCameraSettings(tfxEffectEmitter *effect) {
		return &effect->library->preview_camera_settings[GetEffectInfo(effect)->preview_camera_settings];
	}

	float GetEffectLoopLength(tfxEffectEmitter *effect) {
		return GetEffectProperties(effect)->loop_length[effect->property_index];
	}

	float GetEffectHighestLoopLength(tfxEffectEmitter *effect) {
		float loop_length = GetEffectLoopLength(effect);
		for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
			loop_length = tfxMax(GetEffectHighestLoopLength(&sub), loop_length);
		}
		return loop_length;
	}

	tfxEffectEmitter* AddEmitterToEffect(tfxEffectEmitter *effect, tfxEffectEmitter *emitter) {
		assert(GetEffectInfo(emitter)->name.Length());				//Emitter must have a name so that a hash can be generated
		emitter->type = tfxEffectEmitterType::tfxEmitterType;
		emitter->library = effect->library;
		GetEffectInfo(emitter)->uid = ++effect->library->uid;
		GetEffectInfo(effect)->sub_effectors.push_back(*emitter);
		UpdateLibraryEffectPaths(effect->library);
		ReIndexEffect(effect);
		return &GetEffectInfo(effect)->sub_effectors.back();
	}

	tfxEffectEmitter* AddEffectToEmitter(tfxEffectEmitter *emitter, tfxEffectEmitter *effect) {
		assert(GetEffectInfo(effect)->name.Length());				//Effect must have a name so that a hash can be generated
		effect->type = tfxEffectEmitterType::tfxEffectType;
		effect->library = emitter->library;
		effect->parent = emitter;
		GetEffectInfo(effect)->uid = ++emitter->library->uid;
		GetEffectInfo(emitter)->sub_effectors.push_back(*effect);
		UpdateLibraryEffectPaths(emitter->library);
		ReIndexEffect(emitter);
		return &GetEffectInfo(emitter)->sub_effectors.back();
	}

	tfxEffectEmitter* AddEffect(tfxEffectEmitter *e) {
		tfxEffectEmitter new_effect;
		new_effect.library = e->library;
		GetEffectInfo(&new_effect)->uid = ++e->library->uid;
		new_effect.type = tfxEffectEmitterType::tfxEffectType;
		GetEffectInfo(&new_effect)->name = "New Effect";
		GetEffectInfo(e)->sub_effectors.push_back(new_effect);
		UpdateLibraryEffectPaths(e->library);
		ReIndexEffect(e);
		return &GetEffectInfo(e)->sub_effectors.back();
	}

	tfxU32 CountAllEffects(tfxEffectEmitter *effect, tfxU32 amount) {
		for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
			amount = CountAllEffects(&sub, amount);
		}
		return ++amount;
	}

	int GetEffectDepth(tfxEffectEmitter *e) {
		tfxEffectEmitter *current_parent = e->parent;
		int depth = 0;
		while (current_parent) {
			if (current_parent->type == tfxEmitterType) {
				depth++;
			}
			current_parent = current_parent->parent;
		}
		return depth;
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

	void ResetEffectGraphs(tfxEffectEmitter *effect, bool add_node, bool compile) {
		tfxLibrary *library = effect->library;
		tfxU32 global = effect->global;
		ResetGraph(&library->global_graphs[global].life, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].life.type = tfxGlobal_life;
		ResetGraph(&library->global_graphs[global].amount, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].amount.type = tfxGlobal_amount;
		ResetGraph(&library->global_graphs[global].velocity, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].velocity.type = tfxGlobal_velocity;
		ResetGraph(&library->global_graphs[global].width, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].width.type = tfxGlobal_width;
		ResetGraph(&library->global_graphs[global].height, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].height.type = tfxGlobal_height;
		ResetGraph(&library->global_graphs[global].weight, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].weight.type = tfxGlobal_weight;
		ResetGraph(&library->global_graphs[global].spin, 1.f, tfxGlobalPercentPresetSigned, add_node); library->global_graphs[global].spin.type = tfxGlobal_spin;
		ResetGraph(&library->global_graphs[global].stretch, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].stretch.type = tfxGlobal_stretch;
		ResetGraph(&library->global_graphs[global].overal_scale, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].overal_scale.type = tfxGlobal_overal_scale;
		ResetGraph(&library->global_graphs[global].intensity, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].intensity.type = tfxGlobal_intensity;
		ResetGraph(&library->global_graphs[global].frame_rate, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].frame_rate.type = tfxGlobal_frame_rate;
		ResetGraph(&library->global_graphs[global].splatter, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].splatter.type = tfxGlobal_splatter;
		ResetGraph(&library->global_graphs[global].emitter_width, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].emitter_width.type = tfxGlobal_emitter_width;
		ResetGraph(&library->global_graphs[global].emitter_height, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].emitter_height.type = tfxGlobal_emitter_height;
		ResetGraph(&library->global_graphs[global].emitter_depth, 1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].emitter_depth.type = tfxGlobal_emitter_depth;
		if (compile) {
			CompileLibraryGlobalGraph(library, global);
		}
	}

	void ResetTransformGraphs(tfxEffectEmitter *effect, bool add_node, bool compile) {
		tfxLibrary *library = effect->library;
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

	void ResetEmitterBaseGraphs(tfxEffectEmitter *effect, bool add_node, bool compile) {
		tfxLibrary *library = effect->library;
		tfxU32 emitter_attributes = effect->emitter_attributes;
		ResetGraph(&library->emitter_attributes[emitter_attributes].base.life, 1000.f, tfxLifePreset, add_node); library->emitter_attributes[emitter_attributes].base.life.type = tfxBase_life;
		ResetGraph(&library->emitter_attributes[emitter_attributes].base.amount, 1.f, tfxAmountPreset, add_node); library->emitter_attributes[emitter_attributes].base.amount.type = tfxBase_amount;
		ResetGraph(&library->emitter_attributes[emitter_attributes].base.velocity, 0.f, tfxVelocityPreset, add_node); library->emitter_attributes[emitter_attributes].base.velocity.type = tfxBase_velocity;
		ResetGraph(&library->emitter_attributes[emitter_attributes].base.width, 128.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].base.width.type = tfxBase_width;
		ResetGraph(&library->emitter_attributes[emitter_attributes].base.height, 128.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].base.height.type = tfxBase_height;
		ResetGraph(&library->emitter_attributes[emitter_attributes].base.weight, 0.f, tfxWeightPreset, add_node); library->emitter_attributes[emitter_attributes].base.weight.type = tfxBase_weight;
		ResetGraph(&library->emitter_attributes[emitter_attributes].base.spin, 0.f, tfxSpinPreset, add_node); library->emitter_attributes[emitter_attributes].base.spin.type = tfxBase_spin;
		ResetGraph(&library->emitter_attributes[emitter_attributes].base.noise_offset, 0.f, tfxGlobalPercentPreset, add_node); library->emitter_attributes[emitter_attributes].base.noise_offset.type = tfxBase_noise_offset;
		if (compile) {
			CompileLibraryBaseGraph(library, emitter_attributes);
		}
	}

	void ResetEmitterPropertyGraphs(tfxEffectEmitter *effect, bool add_node, bool compile) {
		tfxLibrary *library = effect->library;
		tfxU32 emitter_attributes = effect->emitter_attributes;
		ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emission_pitch, 0.f, tfxAnglePreset, add_node); library->emitter_attributes[emitter_attributes].properties.emission_pitch.type = tfxProperty_emission_pitch;
		ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emission_yaw, 0.f, tfxAnglePreset, add_node); library->emitter_attributes[emitter_attributes].properties.emission_yaw.type = tfxProperty_emission_yaw;
		ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emission_range, 0.f, tfxEmissionRangePreset, add_node); library->emitter_attributes[emitter_attributes].properties.emission_range.type = tfxProperty_emission_range;
		ResetGraph(&library->emitter_attributes[emitter_attributes].properties.splatter, 0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].properties.splatter.type = tfxProperty_splatter;
		ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_width, 0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].properties.emitter_width.type = tfxProperty_emitter_width;
		ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_height, 0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].properties.emitter_height.type = tfxProperty_emitter_height;
		ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_depth, 0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].properties.emitter_depth.type = tfxProperty_emitter_depth;
		ResetGraph(&library->emitter_attributes[emitter_attributes].properties.arc_size, tfxRadians(360.f), tfxArcPreset, add_node); library->emitter_attributes[emitter_attributes].properties.arc_size.type = tfxProperty_arc_size;
		ResetGraph(&library->emitter_attributes[emitter_attributes].properties.arc_offset, 0.f, tfxArcPreset, add_node); library->emitter_attributes[emitter_attributes].properties.arc_offset.type = tfxProperty_arc_offset;
		if (compile) {
			CompileLibraryPropertyGraph(library, emitter_attributes);
		}
	}

	void ResetEmitterVariationGraphs(tfxEffectEmitter *effect, bool add_node, bool compile) {
		tfxLibrary *library = effect->library;
		tfxU32 emitter_attributes = effect->emitter_attributes;
		ResetGraph(&library->emitter_attributes[emitter_attributes].variation.life, 0.f, tfxLifePreset, add_node); library->emitter_attributes[emitter_attributes].variation.life.type = tfxVariation_life;
		ResetGraph(&library->emitter_attributes[emitter_attributes].variation.amount, 0.f, tfxAmountPreset, add_node); library->emitter_attributes[emitter_attributes].variation.amount.type = tfxVariation_amount;
		ResetGraph(&library->emitter_attributes[emitter_attributes].variation.velocity, 0.f, tfxVelocityPreset, add_node); library->emitter_attributes[emitter_attributes].variation.velocity.type = tfxVariation_velocity;
		ResetGraph(&library->emitter_attributes[emitter_attributes].variation.width, 0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].variation.width.type = tfxVariation_width;
		ResetGraph(&library->emitter_attributes[emitter_attributes].variation.height, 0.f, tfxDimensionsPreset, add_node); library->emitter_attributes[emitter_attributes].variation.height.type = tfxVariation_height;
		ResetGraph(&library->emitter_attributes[emitter_attributes].variation.weight, 0.f, tfxWeightVariationPreset, add_node); library->emitter_attributes[emitter_attributes].variation.weight.type = tfxVariation_weight;
		ResetGraph(&library->emitter_attributes[emitter_attributes].variation.spin, 0.f, tfxSpinVariationPreset, add_node); library->emitter_attributes[emitter_attributes].variation.spin.type = tfxVariation_spin;
		ResetGraph(&library->emitter_attributes[emitter_attributes].variation.noise_offset, 0.f, tfxNoiseOffsetVariationPreset, add_node); library->emitter_attributes[emitter_attributes].variation.noise_offset.type = tfxVariation_noise_offset;
		ResetGraph(&library->emitter_attributes[emitter_attributes].variation.noise_resolution, 300.f, tfxNoiseResolutionPreset, add_node); library->emitter_attributes[emitter_attributes].variation.noise_resolution.type = tfxVariation_noise_resolution;
		if (compile) {
			CompileLibraryVariationGraph(library, emitter_attributes);
		}
	}

	void ResetEmitterOvertimeGraphs(tfxEffectEmitter *effect, bool add_node, bool compile) {
		tfxLibrary *library = effect->library;
		tfxU32 emitter_attributes = effect->emitter_attributes;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity, 1.f, tfxVelocityOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.velocity.type = tfxOvertime_velocity;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster, 1.f, tfxGlobalPercentPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.type = tfxOvertime_velocity_adjuster;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.width, 1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.width.type = tfxOvertime_width;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.height, 1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.height.type = tfxOvertime_height;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.weight, 1.f, tfxWeightOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.weight.type = tfxOvertime_weight;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.spin, 0.f, tfxSpinOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.spin.type = tfxOvertime_spin;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.stretch, 0.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.stretch.type = tfxOvertime_stretch;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.red, 1.f, tfxColorPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.red.type = tfxOvertime_red;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.green, 1.f, tfxColorPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.green.type = tfxOvertime_green;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.blue, 1.f, tfxColorPreset, add_node); library->emitter_attributes[emitter_attributes].overtime.blue.type = tfxOvertime_blue;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.blendfactor, 1.f, tfxOpacityOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.blendfactor.type = tfxOvertime_blendfactor;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.intensity, 1.f, tfxIntensityOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.intensity.type = tfxOvertime_intensity;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance, 30.f, tfxVelocityTurbulancePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance.type = tfxOvertime_velocity_turbulance;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.stretch, 0.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.stretch.type = tfxOvertime_stretch;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.direction_turbulance, 0.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.direction_turbulance.type = tfxOvertime_direction_turbulance;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.direction, 0.f, tfxDirectionOvertimePreset, add_node); library->emitter_attributes[emitter_attributes].overtime.direction.type = tfxOvertime_direction;
		ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.noise_resolution, 1.f, tfxPercentOvertime, add_node); library->emitter_attributes[emitter_attributes].overtime.noise_resolution.type = tfxOvertime_noise_resolution;
		if (compile) {
			CompileLibraryOvertimeGraph(library, emitter_attributes);
		}
	}

	void ResetEmitterGraphs(tfxEffectEmitter *effect, bool add_node, bool compile) {
		ResetEmitterBaseGraphs(effect, add_node, compile);
		ResetEmitterPropertyGraphs(effect, add_node, compile);
		ResetEmitterVariationGraphs(effect, add_node, compile);
		UpdateEffectMaxLife(effect);
		ResetEmitterOvertimeGraphs(effect, add_node, compile);
	}

	void InitialiseUninitialisedGraphs(tfxEffectEmitter *effect) {
		tfxLibrary *library = effect->library;
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
			if (library->global_graphs[global].width.nodes.size() == 0) ResetGraph(&library->global_graphs[global].width, 1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].height.nodes.size() == 0) ResetGraph(&library->global_graphs[global].height, 1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].weight.nodes.size() == 0) ResetGraph(&library->global_graphs[global].weight, 1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].spin.nodes.size() == 0) ResetGraph(&library->global_graphs[global].spin, 1.f, tfxGlobalPercentPresetSigned);
			if (library->global_graphs[global].stretch.nodes.size() == 0) ResetGraph(&library->global_graphs[global].stretch, 1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].overal_scale.nodes.size() == 0) ResetGraph(&library->global_graphs[global].overal_scale, 1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].intensity.nodes.size() == 0) ResetGraph(&library->global_graphs[global].intensity, 1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].frame_rate.nodes.size() == 0) ResetGraph(&library->global_graphs[global].frame_rate, 1.f, tfxGlobalPercentPreset);
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
			if (library->emitter_attributes[emitter_attributes].base.noise_offset.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].base.noise_offset, 0.f, tfxGlobalPercentPreset);

			if (library->emitter_attributes[emitter_attributes].properties.emission_pitch.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emission_pitch, 0.f, tfxAnglePreset);
			if (library->emitter_attributes[emitter_attributes].properties.emission_yaw.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emission_yaw, 0.f, tfxAnglePreset);
			if (library->emitter_attributes[emitter_attributes].properties.emission_range.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emission_range, 0.f, tfxEmissionRangePreset);
			if (library->emitter_attributes[emitter_attributes].properties.splatter.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.splatter, 0.f, tfxDimensionsPreset);
			if (library->emitter_attributes[emitter_attributes].properties.emitter_width.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_width, 0.f, tfxDimensionsPreset);
			if (library->emitter_attributes[emitter_attributes].properties.emitter_height.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_height, 0.f, tfxDimensionsPreset);
			if (library->emitter_attributes[emitter_attributes].properties.emitter_depth.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.emitter_depth, 0.f, tfxDimensionsPreset);
			if (library->emitter_attributes[emitter_attributes].properties.arc_size.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.arc_size, tfxRadians(360.f), tfxArcPreset);
			if (library->emitter_attributes[emitter_attributes].properties.arc_offset.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].properties.arc_offset, 0.f, tfxArcPreset);

			if (library->emitter_attributes[emitter_attributes].variation.life.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.life, 0.f, tfxLifePreset);
			if (library->emitter_attributes[emitter_attributes].variation.amount.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.amount, 0.f, tfxAmountPreset);
			if (library->emitter_attributes[emitter_attributes].variation.velocity.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.velocity, 0.f, tfxVelocityPreset);
			if (library->emitter_attributes[emitter_attributes].variation.width.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.width, 0.f, tfxDimensionsPreset);
			if (library->emitter_attributes[emitter_attributes].variation.height.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.height, 0.f, tfxDimensionsPreset);
			if (library->emitter_attributes[emitter_attributes].variation.weight.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.weight, 0.f, tfxWeightVariationPreset);
			if (library->emitter_attributes[emitter_attributes].variation.spin.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.spin, 0.f, tfxSpinVariationPreset);
			if (library->emitter_attributes[emitter_attributes].variation.noise_offset.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.noise_offset, 0.f, tfxNoiseOffsetVariationPreset);
			if (library->emitter_attributes[emitter_attributes].variation.noise_resolution.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].variation.noise_resolution, 300.f, tfxNoiseResolutionPreset);

			if (library->emitter_attributes[emitter_attributes].overtime.velocity.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity, 1.f, tfxVelocityOvertimePreset);
			if (library->emitter_attributes[emitter_attributes].overtime.width.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.width, 1.f, tfxPercentOvertime);
			if (library->emitter_attributes[emitter_attributes].overtime.height.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.height, 1.f, tfxPercentOvertime);
			if (library->emitter_attributes[emitter_attributes].overtime.weight.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.weight, 1.f, tfxWeightOvertimePreset);
			if (library->emitter_attributes[emitter_attributes].overtime.spin.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.spin, 1.f, tfxSpinOvertimePreset);
			if (library->emitter_attributes[emitter_attributes].overtime.stretch.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.stretch, 0.f, tfxPercentOvertime);
			if (library->emitter_attributes[emitter_attributes].overtime.red.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.red, 1.f, tfxColorPreset);
			if (library->emitter_attributes[emitter_attributes].overtime.green.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.green, 1.f, tfxColorPreset);
			if (library->emitter_attributes[emitter_attributes].overtime.blue.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.blue, 1.f, tfxColorPreset);
			if (library->emitter_attributes[emitter_attributes].overtime.blendfactor.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.blendfactor, 1.f, tfxOpacityOvertimePreset);
			if (library->emitter_attributes[emitter_attributes].overtime.intensity.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.intensity, 1.f, tfxIntensityOvertimePreset);
			if (library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity_turbulance, 0.f, tfxFrameratePreset);
			if (library->emitter_attributes[emitter_attributes].overtime.direction_turbulance.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.direction_turbulance, 0.f, tfxPercentOvertime);
			if (library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.velocity_adjuster, 1.f, tfxGlobalPercentPreset);
			if (library->emitter_attributes[emitter_attributes].overtime.direction.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.direction, 0.f, tfxDirectionOvertimePreset);
			if (library->emitter_attributes[emitter_attributes].overtime.noise_resolution.nodes.size() == 0) ResetGraph(&library->emitter_attributes[emitter_attributes].overtime.noise_resolution, 1.f, tfxPercentOvertime);
		}
	}

	void SetEffectName(tfxEffectEmitter *effect, const char *n) {
		GetEffectInfo(effect)->name = n;
	}

	void AddEmitterColorOvertime(tfxEffectEmitter *effect, float frame, tfxRGB color) {
		AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.red, frame, color.r);
		AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.green, frame, color.g);
		AddGraphNode(&effect->library->emitter_attributes[effect->emitter_attributes].overtime.blue, frame, color.b);
	}

	void SetEffectUserData(tfxEffectEmitter *effect, void *data) {
		effect->user_data = data;
	}

	void* GetEffectUserData(tfxEffectEmitter *effect) {
		return effect->user_data;
	}

	tfxEmitterPropertiesSoA *GetEffectProperties(tfxEffectEmitter *effect) {
		return &effect->library->emitter_properties;
	}

	bool RenameSubEffector(tfxEffectEmitter *emitter, const char *new_name) {
		assert(emitter->parent);	//Must be an emitter or sub effect with a parent
		if (!EffectNameExists(emitter->parent, emitter, new_name) && strlen(new_name) > 0) {
			SetEffectName(emitter, new_name);
			UpdateLibraryEffectPaths(emitter->library);
			return true;
		}

		return false;
	}

	bool EffectNameExists(tfxEffectEmitter *in_effect, tfxEffectEmitter *excluding_effect, const char *name) {
		for (auto &e : GetEffectInfo(in_effect)->sub_effectors) {
			if (excluding_effect != &e) {
				if (GetEffectInfo(&e)->name == name) {
					return true;
				}
			}
		}

		return false;
	}

	void ReIndexEffect(tfxEffectEmitter *effect) {
		tfxU32 index = 0;
		for (auto &e : GetEffectInfo(effect)->sub_effectors) {
			e.library_index = index++;
			e.parent = effect;
			ReIndexEffect(&e);
		}
	}

	void CountEffectChildren(tfxEffectEmitter *effect, int *emitters, int *effects) {
		tmpStack(tfxEffectEmitter*, stack);
		stack.push_back(effect);
		*emitters = 0;
		*effects = 0;
		while (!stack.empty()) {
			tfxEffectEmitter *current = stack.pop_back();
			if (current->type == tfxEffectType)
				*effects++;
			else if (current->type == tfxEmitterType)
				*emitters++;
			for (auto &sub : GetEffectInfo(current)->sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	tfxEffectEmitter* GetRootEffect(tfxEffectEmitter *effect) {
		if (!effect->parent || effect->parent->type == tfxFolder) {
			return effect;
		}
		tfxEffectEmitter *p = effect->parent;
		tfxU32 timeout = 0;
		while (p || ++timeout < 100) {
			if (!p->parent)
				return p;
			p = p->parent;
		}
		return nullptr;
	}

	bool IsRootEffect(tfxEffectEmitter *effect) {
		if (effect->type != tfxEffectType) return false;
		if (effect->type == tfxEffectType && !effect->parent) return true;
		if (effect->parent && effect->parent->type == tfxFolder) return true;
		return false;
	}

	void ResetEffectParents(tfxEffectEmitter *effect) {
		effect->parent = nullptr;
		for (auto &e : GetEffectInfo(effect)->sub_effectors) {
			ResetEffectParents(effect);
		}
	}

	tfxEffectEmitter* MoveEffectUp(tfxEffectEmitter *emitter) {
		tfxEffectEmitter *parent = emitter->parent;
		if (emitter->library_index > 0) {
			tfxU32 new_index = emitter->library_index - 1;
			std::swap(GetEffectInfo(parent)->sub_effectors[emitter->library_index], GetEffectInfo(parent)->sub_effectors[new_index]);
			ReIndexEffect(parent);
			UpdateLibraryEffectPaths(parent->library);
			return &GetEffectInfo(parent)->sub_effectors[new_index];
		}

		return nullptr;
	}

	tfxEffectEmitter* MoveEffectDown(tfxEffectEmitter *emitter) {
		tfxEffectEmitter *parent = emitter->parent;
		if (emitter->library_index < GetEffectInfo(parent)->sub_effectors.size() - 1) {
			tfxU32 new_index = emitter->library_index + 1;
			std::swap(GetEffectInfo(parent)->sub_effectors[emitter->library_index], GetEffectInfo(parent)->sub_effectors[new_index]);
			ReIndexEffect(parent);
			UpdateLibraryEffectPaths(parent->library);
			return &GetEffectInfo(parent)->sub_effectors[new_index];
		}
		return nullptr;
	}

	void DeleteEmitterFromEffect(tfxEffectEmitter *emitter) {
		tfxEffectEmitter *parent = emitter->parent;
		tfxLibrary *library = emitter->library;
		tmpStack(tfxEffectEmitter, stack);
		stack.push_back(*emitter);
		while (stack.size()) {
			tfxEffectEmitter &current = stack.pop_back();
			if (current.type == tfxEffectType && !current.parent) {
				FreeLibraryGlobal(library, current.global);
			}
			else if (current.type == tfxEmitterType) {
				FreeLibraryEmitterAttributes(library, current.emitter_attributes);
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

	void CleanUpEffect(tfxEffectEmitter *effect) {
		if (GetEffectInfo(effect)->sub_effectors.size()) {
			tmpStack(tfxEffectEmitter, stack);
			stack.push_back(*effect);
			while (stack.size()) {
				tfxEffectEmitter current = stack.pop_back();
				if (current.type == tfxEffectType && !current.parent) {
					FreeLibraryGlobal(effect->library, current.global);
				}
				else if (current.type == tfxEmitterType) {
					FreeLibraryEmitterAttributes(effect->library, current.emitter_attributes);
				}
				for (auto &sub : GetEffectInfo(&current)->sub_effectors) {
					stack.push_back(sub);
				}
				GetEffectInfo(&current)->sub_effectors.clear();
				FreeLibraryProperties(effect->library, current.property_index);
				FreeLibraryInfo(effect->library, current.info_index);
			}
		}

		ReIndexEffect(effect);
	}

	void CloneEffect(tfxEffectEmitter *effect_to_clone, tfxEffectEmitter *clone, tfxEffectEmitter *root_parent, tfxLibrary *destination_library, tfxEffectCloningFlags flags) {
		//tfxU32 size = library->global_graphs[0].amount.lookup.values.capacity;
		*clone = *effect_to_clone;
		clone->info_index = CloneLibraryInfo(clone->library, effect_to_clone->info_index, destination_library);
		if (clone->type != tfxFolder) {
			clone->property_index = CloneLibraryProperties(clone->library, effect_to_clone->property_index, destination_library);
		}
		clone->property_flags |= tfxEmitterPropertyFlags_enabled;
		if (!(flags & tfxEffectCloningFlags_keep_user_data))
			clone->user_data = nullptr;
		clone->library = destination_library;
		GetEffectInfo(clone)->sub_effectors.clear();

		tfxLibrary *library = effect_to_clone->library;

		if (effect_to_clone->type == tfxEffectType) {
			if (root_parent == clone) {
				clone->global = flags & tfxEffectCloningFlags_clone_graphs ? CloneLibraryGlobal(library, effect_to_clone->global, destination_library) : clone->global = effect_to_clone->global;
				clone->transform_attributes = flags & tfxEffectCloningFlags_clone_graphs ? CloneLibraryKeyframes(library, effect_to_clone->transform_attributes, destination_library) : clone->transform_attributes = effect_to_clone->transform_attributes;
				if (flags & tfxEffectCloningFlags_compile_graphs) {
					CompileLibraryGlobalGraph(clone->library, clone->global);
					CompileLibraryKeyframeGraph(clone->library, clone->transform_attributes);
				}
			}
			else {
				clone->transform_attributes = flags & tfxEffectCloningFlags_clone_graphs ? CloneLibraryKeyframes(library, effect_to_clone->transform_attributes, destination_library) : clone->transform_attributes = effect_to_clone->transform_attributes;
				if (flags & tfxEffectCloningFlags_compile_graphs) {
					CompileLibraryKeyframeGraph(clone->library, clone->transform_attributes);
				}
				if (!(flags & tfxEffectCloningFlags_force_clone_global)) {
					clone->global = root_parent->global;
				}
				else {
					clone->global = CloneLibraryGlobal(library, root_parent->global, destination_library);
					if (flags & tfxEffectCloningFlags_compile_graphs)
						CompileLibraryGlobalGraph(clone->library, clone->global);
				}
			}
		}
		else if (effect_to_clone->type == tfxEmitterType) {
			clone->emitter_attributes = flags & tfxEffectCloningFlags_clone_graphs ? CloneLibraryEmitterAttributes(library, effect_to_clone->emitter_attributes, destination_library) : effect_to_clone->emitter_attributes;
			clone->transform_attributes = flags & tfxEffectCloningFlags_clone_graphs ? CloneLibraryKeyframes(library, effect_to_clone->transform_attributes, destination_library) : effect_to_clone->transform_attributes;
			UpdateEffectMaxLife(clone);
			if (flags & tfxEffectCloningFlags_compile_graphs) {
				CompileLibraryKeyframeGraph(clone->library, clone->transform_attributes);
				CompileLibraryPropertyGraph(clone->library, clone->emitter_attributes);
				CompileLibraryBaseGraph(clone->library, clone->emitter_attributes);
				CompileLibraryVariationGraph(clone->library, clone->emitter_attributes);
				CompileLibraryOvertimeGraph(clone->library, clone->emitter_attributes);
			}
		}

		for (auto &e : GetEffectInfo(effect_to_clone)->sub_effectors) {
			if (e.type == tfxEmitterType) {
				tfxEffectEmitter emitter_copy;
				CloneEffect(&e, &emitter_copy, root_parent, destination_library, flags);
				if (!(flags & tfxEffectCloningFlags_keep_user_data))
					emitter_copy.user_data = nullptr;
				AddEmitterToEffect(clone, &emitter_copy);
			}
			else if (e.type == tfxEffectType) {
				tfxEffectEmitter effect_copy;
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

	void AddTemplatePath(tfxEffectTemplate *effect_template, tfxEffectEmitter *effect_emitter, tfxStr256 path) {
		effect_template->paths.Insert(path, effect_emitter);
		for (auto &sub : GetEffectInfo(effect_emitter)->sub_effectors) {
			tfxStr256 sub_path = path;
			sub_path.Appendf("/%s", GetEffectInfo(&sub)->name.c_str());
			AddTemplatePath(effect_template, &sub, sub_path);
		}
	}

	bool PrepareEffectTemplate(tfxLibrary *library, const char *name, tfxEffectTemplate *effect_template) {
		ResetTemplate(effect_template);
		if (library->effect_paths.ValidName(name)) {
			PrepareLibraryEffectTemplate(library, name, *effect_template);
			return true;
		}
		return false;
	}

	void EnableAllEmitters(tfxEffectEmitter *effect) {
		for (auto &e : GetEffectInfo(effect)->sub_effectors) {
			e.property_flags |= tfxEmitterPropertyFlags_enabled;
			EnableAllEmitters(&e);
		}
	}

	void EnableEmitter(tfxEffectEmitter *effect) {
		effect->property_flags |= tfxEmitterPropertyFlags_enabled;
	}

	void DisableAllEmitters(tfxEffectEmitter *effect) {
		for (auto &e : GetEffectInfo(effect)->sub_effectors) {
			e.property_flags &= ~tfxEmitterPropertyFlags_enabled;
			DisableAllEmitters(&e);
		}
	}

	void DisableAllEmittersExcept(tfxEffectEmitter *effect, tfxEffectEmitter *emitter) {
		for (auto &e : GetEffectInfo(effect)->sub_effectors) {
			if (e.library_index == emitter->library_index)
				e.property_flags |= tfxEmitterPropertyFlags_enabled;
			else
				e.property_flags &= ~tfxEmitterPropertyFlags_enabled;
		}
	}

	tfxGraph* GetEffectGraphByType(tfxEffectEmitter *effect, tfxGraphType type) {
		tfxLibrary *library = effect->library;

		if (type < tfxGlobalCount) {
			return &((tfxGraph*)&library->global_graphs[effect->global])[type];
		}
		else if (type >= tfxPropertyStart && type < tfxBaseStart) {
			int ref = type - tfxPropertyStart;
			return &((tfxGraph*)&library->emitter_attributes[effect->emitter_attributes].properties)[ref];
		}
		else if (type >= tfxBaseStart && type < tfxVariationStart) {
			int ref = type - tfxBaseStart;
			return &((tfxGraph*)&library->emitter_attributes[effect->emitter_attributes].base)[ref];
		}
		else if (type >= tfxVariationStart && type < tfxOvertimeStart) {
			int ref = type - tfxVariationStart;
			return &((tfxGraph*)&library->emitter_attributes[effect->emitter_attributes].variation)[ref];
		}
		else if (type >= tfxOvertimeStart && type < tfxTransformStart) {
			int ref = type - tfxOvertimeStart;
			return &((tfxGraph*)&library->emitter_attributes[effect->emitter_attributes].overtime)[ref];
		}
		else if (type >= tfxTransformStart) {
			int ref = type - tfxTransformStart;
			return &((tfxGraph*)&library->transform_attributes[effect->transform_attributes].roll)[ref];
		}

		return nullptr;

	}

	tfxU32 GetEffectGraphIndexByType(tfxEffectEmitter *effect, tfxGraphType type) {

		if (type < tfxGlobalCount) {
			return effect->global;
		}
		else if (type < tfxTransformStart) {
			return effect->emitter_attributes;
		}
		else {
			return effect->transform_attributes;
		}

	}

	void FreeEffectGraphs(tfxEffectEmitter *effect) {

		tfxLibrary *library = effect->library;

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
			FreeGraph(&library->global_graphs[global].frame_rate);
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
		}
	}

	tfxU32 CountAllEffectLookupValues(tfxEffectEmitter *effect) {
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

	void CompileEffectGraphs(tfxEffectEmitter *effect) {
		for (tfxU32 t = (tfxU32)tfxTransform_translate_x; t != (tfxU32)tfxGraphMaxIndex; ++t) {
			CompileGraph(*GetEffectGraphByType(effect, tfxGraphType(t)));
		}
		if (effect->type == tfxEffectType) {
			for (tfxU32 t = (tfxU32)tfxGlobal_life; t != (tfxU32)tfxProperty_emission_pitch; ++t) {
				CompileGraph(*GetEffectGraphByType(effect, tfxGraphType(t)));
			}
		}
		else if (effect->type == tfxEmitterType) {
			for (tfxU32 t = (tfxU32)tfxProperty_emission_pitch; t != (tfxU32)tfxOvertime_velocity; ++t) {
				CompileGraph(*GetEffectGraphByType(effect, (tfxGraphType)t));
			}
			for (tfxU32 t = (tfxU32)tfxOvertime_velocity; t != (tfxU32)tfxTransform_translate_x; ++t) {
				if (IsColorGraph((tfxGraphType)t)) {
					CompileColorOvertime(*GetEffectGraphByType(effect, (tfxGraphType)t));
				}
				else {
					CompileGraphOvertime(*GetEffectGraphByType(effect, (tfxGraphType)t));
				}
			}
		}
		for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
			CompileEffectGraphs(&sub);
		}
	}

	void InitialiseGlobalAttributes(tfxGlobalAttributes *attributes, tfxMemoryArenaManager *allocator, tfxMemoryArenaManager *value_allocator, tfxU32 bucket_size) {
		attributes->life.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->amount.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->velocity.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->width.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->height.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->weight.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->spin.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->stretch.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->overal_scale.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->intensity.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->frame_rate.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->splatter.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->emitter_width.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->emitter_height.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->emitter_depth.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);

		attributes->life.lookup.values.allocator = value_allocator;
		attributes->amount.lookup.values.allocator = value_allocator;
		attributes->velocity.lookup.values.allocator = value_allocator;
		attributes->width.lookup.values.allocator = value_allocator;
		attributes->height.lookup.values.allocator = value_allocator;
		attributes->weight.lookup.values.allocator = value_allocator;
		attributes->spin.lookup.values.allocator = value_allocator;
		attributes->stretch.lookup.values.allocator = value_allocator;
		attributes->overal_scale.lookup.values.allocator = value_allocator;
		attributes->intensity.lookup.values.allocator = value_allocator;
		attributes->frame_rate.lookup.values.allocator = value_allocator;
		attributes->splatter.lookup.values.allocator = value_allocator;
		attributes->emitter_width.lookup.values.allocator = value_allocator;
		attributes->emitter_height.lookup.values.allocator = value_allocator;
		attributes->emitter_depth.lookup.values.allocator = value_allocator;
	}

	void FreeGlobalAttributes(tfxGlobalAttributes *attributes) {
		FreeGraph(&attributes->life);
		FreeGraph(&attributes->amount);
		FreeGraph(&attributes->velocity);
		FreeGraph(&attributes->width);
		FreeGraph(&attributes->height);
		FreeGraph(&attributes->weight);
		FreeGraph(&attributes->spin);
		FreeGraph(&attributes->stretch);
		FreeGraph(&attributes->overal_scale);
		FreeGraph(&attributes->intensity);
		FreeGraph(&attributes->frame_rate);
		FreeGraph(&attributes->splatter);
		FreeGraph(&attributes->emitter_width);
		FreeGraph(&attributes->emitter_height);
		FreeGraph(&attributes->emitter_depth);
	}

	void CopyGlobalAttributesNoLookups(tfxGlobalAttributes *src, tfxGlobalAttributes *dst) {
		CopyGraphNoLookups(&src->life, &dst->life);
		CopyGraphNoLookups(&src->amount, &dst->amount);
		CopyGraphNoLookups(&src->velocity, &dst->velocity);
		CopyGraphNoLookups(&src->width, &dst->width);
		CopyGraphNoLookups(&src->height, &dst->height);
		CopyGraphNoLookups(&src->weight, &dst->weight);
		CopyGraphNoLookups(&src->spin, &dst->spin);
		CopyGraphNoLookups(&src->stretch, &dst->stretch);
		CopyGraphNoLookups(&src->overal_scale, &dst->overal_scale);
		CopyGraphNoLookups(&src->intensity, &dst->intensity);
		CopyGraphNoLookups(&src->frame_rate, &dst->frame_rate);
		CopyGraphNoLookups(&src->splatter, &dst->splatter);
		CopyGraphNoLookups(&src->emitter_width, &dst->emitter_width);
		CopyGraphNoLookups(&src->emitter_height, &dst->emitter_height);
		CopyGraphNoLookups(&src->emitter_depth, &dst->emitter_depth);
	}

	void InitialiseTransformAttributes(tfxTransformAttributes *attributes, tfxMemoryArenaManager *allocator, tfxMemoryArenaManager *value_allocator, tfxU32 bucket_size) {
		attributes->roll.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->pitch.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->yaw.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->translation_x.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->translation_y.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->translation_z.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);

		attributes->roll.lookup.values.allocator = value_allocator;
		attributes->pitch.lookup.values.allocator = value_allocator;
		attributes->yaw.lookup.values.allocator = value_allocator;
		attributes->translation_x.lookup.values.allocator = value_allocator;
		attributes->translation_y.lookup.values.allocator = value_allocator;
		attributes->translation_z.lookup.values.allocator = value_allocator;
	}

	void FreeTransformAttributes(tfxTransformAttributes *attributes) {
		FreeGraph(&attributes->roll);
		FreeGraph(&attributes->pitch);
		FreeGraph(&attributes->yaw);
		FreeGraph(&attributes->translation_x);
		FreeGraph(&attributes->translation_y);
		FreeGraph(&attributes->translation_z);
	}

	void CopyTransformAttributesNoLookups(tfxTransformAttributes *src, tfxTransformAttributes *dst) {
		CopyGraphNoLookups(&src->roll, &dst->roll);
		CopyGraphNoLookups(&src->pitch, &dst->pitch);
		CopyGraphNoLookups(&src->yaw, &dst->yaw);
		CopyGraphNoLookups(&src->translation_x, &dst->translation_x);
		CopyGraphNoLookups(&src->translation_y, &dst->translation_y);
		CopyGraphNoLookups(&src->translation_z, &dst->translation_z);
	}

	bool HasTranslationKeyframes(tfxTransformAttributes &graphs) {
		return graphs.translation_x.nodes.size() || graphs.translation_y.nodes.size() || graphs.translation_z.nodes.size();
	}

	void AddTranslationNodes(tfxTransformAttributes &keyframes, float frame) {
		if (keyframes.translation_x.nodes.size()) {
			if (!HasNodeAtFrame(keyframes.translation_x, frame))
				AddGraphCoordNode(&keyframes.translation_x, frame, 0.f);
			if (!HasNodeAtFrame(keyframes.translation_y, frame))
				AddGraphCoordNode(&keyframes.translation_y, frame, 0.f);
			if (!HasNodeAtFrame(keyframes.translation_z, frame))
				AddGraphCoordNode(&keyframes.translation_z, frame, 0.f);
		}
		else {
			AddGraphCoordNode(&keyframes.translation_x, 0.f, 0.f);
			AddGraphCoordNode(&keyframes.translation_y, 0.f, 0.f);
			AddGraphCoordNode(&keyframes.translation_z, 0.f, 0.f);
			if (frame != 0) {
				AddGraphCoordNode(&keyframes.translation_x, frame, 0.f);
				AddGraphCoordNode(&keyframes.translation_y, frame, 0.f);
				AddGraphCoordNode(&keyframes.translation_z, frame, 0.f);
			}
		}
	}

	void InitialisePropertyAttributes(tfxPropertyAttributes *attributes, tfxMemoryArenaManager *allocator, tfxMemoryArenaManager *value_allocator, tfxU32 bucket_size) {
		attributes->emission_pitch.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->emission_yaw.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->emission_range.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->splatter.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->emitter_width.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->emitter_height.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->emitter_depth.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->arc_size.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->arc_offset.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);

		attributes->emission_pitch.lookup.values.allocator = value_allocator;
		attributes->emission_yaw.lookup.values.allocator = value_allocator;
		attributes->emission_range.lookup.values.allocator = value_allocator;
		attributes->splatter.lookup.values.allocator = value_allocator;
		attributes->emitter_width.lookup.values.allocator = value_allocator;
		attributes->emitter_height.lookup.values.allocator = value_allocator;
		attributes->emitter_depth.lookup.values.allocator = value_allocator;
		attributes->arc_size.lookup.values.allocator = value_allocator;
		attributes->arc_offset.lookup.values.allocator = value_allocator;
	}

	void FreePropertyAttributes(tfxPropertyAttributes *attributes) {
		FreeGraph(&attributes->emission_pitch);
		FreeGraph(&attributes->emission_yaw);
		FreeGraph(&attributes->emission_range);
		FreeGraph(&attributes->splatter);
		FreeGraph(&attributes->emitter_width);
		FreeGraph(&attributes->emitter_height);
		FreeGraph(&attributes->emitter_depth);
		FreeGraph(&attributes->arc_size);
		FreeGraph(&attributes->arc_offset);
	}

	void CopyPropertyAttributesNoLookups(tfxPropertyAttributes *src, tfxPropertyAttributes *dst) {
		CopyGraphNoLookups(&src->emission_pitch, &dst->emission_pitch);
		CopyGraphNoLookups(&src->emission_yaw, &dst->emission_yaw);
		CopyGraphNoLookups(&src->emission_range, &dst->emission_range);
		CopyGraphNoLookups(&src->splatter, &dst->splatter);
		CopyGraphNoLookups(&src->emitter_width, &dst->emitter_width);
		CopyGraphNoLookups(&src->emitter_height, &dst->emitter_height);
		CopyGraphNoLookups(&src->emitter_depth, &dst->emitter_depth);
		CopyGraphNoLookups(&src->arc_size, &dst->arc_size);
		CopyGraphNoLookups(&src->arc_offset, &dst->arc_offset);
	}

	void InitialiseBaseAttributes(tfxBaseAttributes *attributes, tfxMemoryArenaManager *allocator, tfxMemoryArenaManager *value_allocator, tfxU32 bucket_size) {
		attributes->life.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->amount.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->velocity.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->width.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->height.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->weight.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->spin.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->noise_offset.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);

		attributes->life.lookup.values.allocator = value_allocator;
		attributes->amount.lookup.values.allocator = value_allocator;
		attributes->velocity.lookup.values.allocator = value_allocator;
		attributes->width.lookup.values.allocator = value_allocator;
		attributes->height.lookup.values.allocator = value_allocator;
		attributes->weight.lookup.values.allocator = value_allocator;
		attributes->spin.lookup.values.allocator = value_allocator;
		attributes->noise_offset.lookup.values.allocator = value_allocator;
	}

	void InitialiseVariationAttributes(tfxVariationAttributes *attributes, tfxMemoryArenaManager *allocator, tfxMemoryArenaManager *value_allocator, tfxU32 bucket_size) {
		attributes->life.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->amount.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->velocity.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->width.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->height.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->weight.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->spin.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->noise_offset.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->noise_resolution.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);

		attributes->life.lookup.values.allocator = value_allocator;
		attributes->amount.lookup.values.allocator = value_allocator;
		attributes->velocity.lookup.values.allocator = value_allocator;
		attributes->width.lookup.values.allocator = value_allocator;
		attributes->height.lookup.values.allocator = value_allocator;
		attributes->weight.lookup.values.allocator = value_allocator;
		attributes->spin.lookup.values.allocator = value_allocator;
		attributes->noise_offset.lookup.values.allocator = value_allocator;
		attributes->noise_resolution.lookup.values.allocator = value_allocator;
	}

	void InitialiseOvertimeAttributes(tfxOvertimeAttributes *attributes, tfxMemoryArenaManager *allocator, tfxMemoryArenaManager *value_allocator, tfxU32 bucket_size) {
		attributes->velocity.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->width.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->height.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->weight.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->spin.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->stretch.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->red.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->blue.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->green.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->blendfactor.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->velocity_turbulance.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->direction_turbulance.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->velocity_adjuster.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->intensity.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->direction.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);
		attributes->noise_resolution.nodes = tfxBucketArray<tfxAttributeNode>(allocator, bucket_size);

		attributes->velocity.lookup.values.allocator = value_allocator;
		attributes->width.lookup.values.allocator = value_allocator;
		attributes->height.lookup.values.allocator = value_allocator;
		attributes->weight.lookup.values.allocator = value_allocator;
		attributes->spin.lookup.values.allocator = value_allocator;
		attributes->stretch.lookup.values.allocator = value_allocator;
		attributes->red.lookup.values.allocator = value_allocator;
		attributes->blue.lookup.values.allocator = value_allocator;
		attributes->green.lookup.values.allocator = value_allocator;
		attributes->blendfactor.lookup.values.allocator = value_allocator;
		attributes->velocity_turbulance.lookup.values.allocator = value_allocator;
		attributes->direction_turbulance.lookup.values.allocator = value_allocator;
		attributes->velocity_adjuster.lookup.values.allocator = value_allocator;
		attributes->intensity.lookup.values.allocator = value_allocator;
		attributes->direction.lookup.values.allocator = value_allocator;
		attributes->noise_resolution.lookup.values.allocator = value_allocator;
	}

	void FreeOvertimeAttributes(tfxOvertimeAttributes *attributes) {
		FreeGraph(&attributes->velocity);
		FreeGraph(&attributes->width);
		FreeGraph(&attributes->height);
		FreeGraph(&attributes->weight);
		FreeGraph(&attributes->spin);
		FreeGraph(&attributes->stretch);
		FreeGraph(&attributes->red);
		FreeGraph(&attributes->green);
		FreeGraph(&attributes->blue);
		FreeGraph(&attributes->blendfactor);
		FreeGraph(&attributes->velocity);
		FreeGraph(&attributes->direction_turbulance);
		FreeGraph(&attributes->velocity_adjuster);
		FreeGraph(&attributes->intensity);
		FreeGraph(&attributes->direction);
		FreeGraph(&attributes->noise_resolution);
	}

	void CopyOvertimeAttributesNoLookups(tfxOvertimeAttributes *src, tfxOvertimeAttributes *dst) {
		CopyGraphNoLookups(&src->velocity, &dst->velocity);
		CopyGraphNoLookups(&src->width, &dst->width);
		CopyGraphNoLookups(&src->height, &dst->height);
		CopyGraphNoLookups(&src->weight, &dst->weight);
		CopyGraphNoLookups(&src->spin, &dst->spin);
		CopyGraphNoLookups(&src->stretch, &dst->stretch);
		CopyGraphNoLookups(&src->red, &dst->red);
		CopyGraphNoLookups(&src->green, &dst->green);
		CopyGraphNoLookups(&src->blue, &dst->blue);
		CopyGraphNoLookups(&src->blendfactor, &dst->blendfactor);
		CopyGraphNoLookups(&src->velocity_turbulance, &dst->velocity_turbulance);
		CopyGraphNoLookups(&src->direction_turbulance, &dst->direction_turbulance);
		CopyGraphNoLookups(&src->velocity_adjuster, &dst->velocity_adjuster);
		CopyGraphNoLookups(&src->intensity, &dst->intensity);
		CopyGraphNoLookups(&src->direction, &dst->direction);
		CopyGraphNoLookups(&src->noise_resolution, &dst->noise_resolution);
	}

	void FreeVariationAttributes(tfxVariationAttributes *attributes) {
		FreeGraph(&attributes->life);
		FreeGraph(&attributes->amount);
		FreeGraph(&attributes->velocity);
		FreeGraph(&attributes->width);
		FreeGraph(&attributes->height);
		FreeGraph(&attributes->weight);
		FreeGraph(&attributes->spin);
		FreeGraph(&attributes->noise_offset);
		FreeGraph(&attributes->noise_resolution);
	}

	void CopyVariationAttributesNoLookups(tfxVariationAttributes *src, tfxVariationAttributes *dst) {
		CopyGraphNoLookups(&src->life, &dst->life);
		CopyGraphNoLookups(&src->amount, &dst->amount);
		CopyGraphNoLookups(&src->velocity, &dst->velocity);
		CopyGraphNoLookups(&src->width, &dst->width);
		CopyGraphNoLookups(&src->height, &dst->height);
		CopyGraphNoLookups(&src->weight, &dst->weight);
		CopyGraphNoLookups(&src->spin, &dst->spin);
		CopyGraphNoLookups(&src->noise_offset, &dst->noise_offset);
		CopyGraphNoLookups(&src->noise_resolution, &dst->noise_resolution);
	}

	void FreeBaseAttributes(tfxBaseAttributes *attributes) {
		FreeGraph(&attributes->life);
		FreeGraph(&attributes->amount);
		FreeGraph(&attributes->velocity);
		FreeGraph(&attributes->width);
		FreeGraph(&attributes->height);
		FreeGraph(&attributes->weight);
		FreeGraph(&attributes->spin);
		FreeGraph(&attributes->noise_offset);
	}

	void CopyBaseAttributesNoLookups(tfxBaseAttributes *src, tfxBaseAttributes *dst) {
		CopyGraphNoLookups(&src->life, &dst->life);
		CopyGraphNoLookups(&src->amount, &dst->amount);
		CopyGraphNoLookups(&src->velocity, &dst->velocity);
		CopyGraphNoLookups(&src->width, &dst->width);
		CopyGraphNoLookups(&src->height, &dst->height);
		CopyGraphNoLookups(&src->weight, &dst->weight);
		CopyGraphNoLookups(&src->spin, &dst->spin);
		CopyGraphNoLookups(&src->noise_offset, &dst->noise_offset);
	}

	void InitialiseEmitterAttributes(tfxEmitterAttributes *attributes, tfxMemoryArenaManager *allocator, tfxMemoryArenaManager *value_allocator, tfxU32 bucket_size) {
		InitialisePropertyAttributes(&attributes->properties, allocator, value_allocator, bucket_size);
		InitialiseBaseAttributes(&attributes->base, allocator, value_allocator, bucket_size);
		InitialiseVariationAttributes(&attributes->variation, allocator, value_allocator, bucket_size);
		InitialiseOvertimeAttributes(&attributes->overtime, allocator, value_allocator, bucket_size);
	}

	void FreeEmitterAttributes(tfxEmitterAttributes *attributes) {
		FreePropertyAttributes(&attributes->properties);
		FreeBaseAttributes(&attributes->base);
		FreeVariationAttributes(&attributes->variation);
		FreeOvertimeAttributes(&attributes->overtime);
	}

	tfxEffectEmitter& tfxLibrary::operator[] (tfxU32 index) {
		return effects[index];
	}

	void MaybeGrowLibraryProperties(tfxLibrary *library, tfxU32 size_offset) {
		if (library->emitter_properties_buffer.current_size >= library->emitter_properties_buffer.capacity - size_offset) {
			GrowArrays(&library->emitter_properties_buffer, library->emitter_properties_buffer.capacity, library->emitter_properties_buffer.capacity + 1);
		}
	}

	void MaybeGrowLibraryInfos(tfxLibrary *library) {
		if (library->effect_infos.current_size >= library->effect_infos.capacity - 4) {
			library->effect_infos.reserve(library->effect_infos._grow_capacity(library->effect_infos.current_size + 1));
		}
	}

	tfxEffectEmitterInfo *GetEffectInfo(tfxEffectEmitter *e) {
		assert(e->library->effect_infos.size() > e->info_index);
		return &e->library->effect_infos[e->info_index];
	}

	bool RenameLibraryEffect(tfxLibrary *library, tfxEffectEmitter &effect, const char *new_name) {
		if (!LibraryNameExists(library, effect, new_name) && strlen(new_name) > 0) {
			SetEffectName(&effect, new_name);
			UpdateLibraryEffectPaths(library);
			return true;
		}

		return false;
	}

	bool LibraryNameExists(tfxLibrary *library, tfxEffectEmitter& effect, const char *name) {
		for (auto &e : library->effects) {
			if (effect.library_index != e.library_index) {
				if (GetEffectInfo(&e)->name == name) {
					return true;
				}
			}
		}

		return false;
	}

	bool LibraryNameExists2(tfxLibrary *library, tfxEffectEmitter &effect, const char *name) {
		for (auto &e : library->effects) {
			if (GetEffectInfo(&e)->name == name) {
				return true;
			}
		}
		return false;
	}

	void UpdateLibraryEffectPaths(tfxLibrary *library) {
		library->effect_paths.Clear();
		for (auto &e : library->effects) {
			tfxStr256 path = GetEffectInfo(&e)->name;
			e.path_hash = tfxXXHash64::hash(path.c_str(), path.Length(), 0);
			AddLibraryPath(library, e, path);
		}
	}

	void AddLibraryPath(tfxLibrary *library, tfxEffectEmitter &effect_emitter, tfxStr256 &path) {
		library->effect_paths.Insert(path, &effect_emitter);
		for (auto &sub : GetEffectInfo(&effect_emitter)->sub_effectors) {
			tfxStr256 sub_path = path;
			sub_path.Appendf("/%s", GetEffectInfo(&sub)->name.c_str());
			sub.path_hash = tfxXXHash64::hash(sub_path.c_str(), sub_path.Length(), 0);
			AddLibraryPath(library, sub, sub_path);
		}
	}

	tfxEffectEmitter &InsertLibraryEffect(tfxLibrary *library, tfxEffectEmitter &effect, tfxEffectEmitter *position) {
		effect.library_index = library->effects.current_size;
		effect.type = tfxEffectType;
		GetEffectInfo(&effect)->uid = ++library->uid;
		effect.library = library;
		tfxEffectEmitter *inserted_effect = library->effects.insert_after(position, effect);
		ReIndexLibrary(library);
		UpdateLibraryEffectPaths(library);
		return *inserted_effect;
	}

	tfxEffectEmitter &AddLibraryEffect(tfxLibrary *library, tfxEffectEmitter &effect) {
		effect.library_index = library->effects.current_size;
		effect.type = tfxEffectType;
		GetEffectInfo(&effect)->uid = ++library->uid;
		effect.library = library;
		library->effects.push_back(effect);
		ReIndexLibrary(library);
		UpdateLibraryEffectPaths(library);
		return library->effects.back();
	}

	tfxEffectEmitter &AddLibraryFolder(tfxLibrary *library, tfxStr64 &name) {
		tfxEffectEmitter folder;
		folder.info_index = AddLibraryEffectEmitterInfo(library);
		folder.library = library;
		GetEffectInfo(&folder)->name = name;
		folder.type = tfxFolder;
		folder.library = library;
		GetEffectInfo(&folder)->uid = ++library->uid;
		library->effects.push_back(folder);
		ReIndexLibrary(library);
		UpdateLibraryEffectPaths(library);
		return library->effects.back();
	}

	tfxEffectEmitter &AddLibraryFolder(tfxLibrary *library, tfxEffectEmitter &folder) {
		assert(folder.type == tfxFolder);			//Must be type tfxFolder if adding a folder
		folder.library = library;
		GetEffectInfo(&folder)->uid = ++library->uid;
		library->effects.push_back(folder);
		ReIndexLibrary(library);
		UpdateLibraryEffectPaths(library);
		return library->effects.back();
	}

	tfxEffectEmitter &AddLibraryStage(tfxLibrary *library, tfxStr64 &name) {
		tfxEffectEmitter stage;
		stage.info_index = AddLibraryEffectEmitterInfo(library);
		stage.library = library;
		GetEffectInfo(&stage)->name = name;
		stage.type = tfxStage;
		GetEffectInfo(&stage)->uid = ++library->uid;
		library->effects.push_back(stage);
		ReIndexLibrary(library);
		UpdateLibraryEffectPaths(library);
		return library->effects.back();
	}

	tfxEffectEmitter* GetLibraryEffect(tfxLibrary *library, tfxStr256 &path) {
		assert(library->effect_paths.ValidName(path));		//Effect was not found by that name
		return library->effect_paths.At(path);
	}

	tfxEffectEmitter* GetLibraryEffect(tfxLibrary *library, const char *path) {
		assert(library->effect_paths.ValidName(path));		//Effect was not found by that name
		return library->effect_paths.At(path);
	}

	tfxEffectEmitter* GetLibraryEffect(tfxLibrary *library, tfxKey key) {
		assert(library->effect_paths.ValidKey(key));			//Effect was not found by that key
		return library->effect_paths.At(key);
	}

	void PrepareLibraryEffectTemplate(tfxLibrary *library, tfxStr256 path, tfxEffectTemplate &effect_template) {
		tfxEffectEmitter *effect = GetLibraryEffect(library, path);
		assert(effect);								//Effect was not found, make sure the path exists
		assert(effect->type == tfxEffectType);		//The effect must be an effect type, not an emitter
		effect_template.original_effect_hash = effect->path_hash;
		CloneEffect(effect, &effect_template.effect, &effect_template.effect, library, tfxEffectCloningFlags_clone_graphs | tfxEffectCloningFlags_compile_graphs);
		AddTemplatePath(&effect_template, &effect_template.effect, GetEffectInfo(&effect_template.effect)->name.c_str());
	}

	void PrepareLibraryEffectTemplate(tfxLibrary *library, tfxEffectEmitter &effect, tfxEffectTemplate &effect_template) {
		assert(effect.type == tfxEffectType);
		CloneEffect(&effect, &effect_template.effect, &effect_template.effect, library);
		AddTemplatePath(&effect_template, &effect_template.effect, GetEffectInfo(&effect_template.effect)->name.c_str());
	}

	void ReIndexLibrary(tfxLibrary *library) {
		tfxU32 index = 0;
		for (auto &e : library->effects) {
			e.library_index = index++;
			e.parent = nullptr;
			ReIndexEffect(&e);
		}
	}

	void UpdateLibraryParticleShapeReferences(tfxLibrary *library, tfxvec<tfxEffectEmitter> &effects, tfxKey default_hash) {
		tfxvec<tfxEffectEmitter*> stack;
		for (auto &effect : library->effects) {
			stack.push_back(&effect);
		}
		while (stack.size()) {
			tfxEffectEmitter &current = *stack.pop_back();
			if (current.type == tfxEmitterType) {
				bool shape_found = false;
				tfxKey hash = library->emitter_properties.image_hash[current.property_index];
				if (library->particle_shapes.ValidKey(library->emitter_properties.image_hash[current.property_index])) {
					library->emitter_properties.image[current.property_index] = &library->particle_shapes.At(library->emitter_properties.image_hash[current.property_index]);
					library->emitter_properties.end_frame[current.property_index] = library->particle_shapes.At(library->emitter_properties.image_hash[current.property_index]).animation_frames - 1;
					shape_found = true;
				}
				else {
					for (auto &shape : library->particle_shapes.data) {
						if (shape.shape_index == library->emitter_properties.image_index[current.property_index]) {
							library->emitter_properties.image_hash[current.property_index] = shape.image_hash;
							library->emitter_properties.image[current.property_index] = &library->particle_shapes.At(library->emitter_properties.image_hash[current.property_index]);
							library->emitter_properties.end_frame[current.property_index] = library->particle_shapes.At(library->emitter_properties.image_hash[current.property_index]).animation_frames - 1;
							shape_found = true;
							break;
						}
					}
				}
				if (!shape_found) {
					library->emitter_properties.image[current.property_index] = &library->particle_shapes.At(default_hash);
					library->emitter_properties.end_frame[current.property_index] = library->particle_shapes.At(default_hash).animation_frames - 1;
				}
			}
			for (auto &sub : GetEffectInfo(&current)->sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	tfxEffectEmitter* LibraryMoveUp(tfxLibrary *library, tfxEffectEmitter &effect) {
		if (effect.library_index > 0) {
			tfxU32 new_index = effect.library_index - 1;
			std::swap(library->effects[effect.library_index], library->effects[new_index]);
			UpdateLibraryEffectPaths(library);
			ReIndexLibrary(library);
			return &library->effects[new_index];
		}
		return nullptr;
	}

	tfxEffectEmitter* LibraryMoveDown(tfxLibrary *library, tfxEffectEmitter &effect) {
		if (effect.library_index < library->effects.size() - 1) {
			tfxU32 new_index = effect.library_index + 1;
			std::swap(library->effects[effect.library_index], library->effects[new_index]);
			UpdateLibraryEffectPaths(library);
			ReIndexLibrary(library);
			return &library->effects[new_index];
		}
		return nullptr;
	}

	tfxGPUShapes BuildGPUShapeData(tfxvec<tfxImageData> *particle_shapes, tfxVec4(uv_lookup)(void *ptr, tfxGPUImageData *image_data, int offset)) {
		assert(particle_shapes->size());		//There are no shapes to copy!
		tfxU32 index = 0;
		tfxGPUShapes shape_data;
		for (auto &shape : *particle_shapes) {
			if (shape.animation_frames == 1) {
				tfxGPUImageData cs;
				cs.animation_frames = shape.animation_frames;
				cs.image_size = shape.image_size;
				cs.uv = uv_lookup(shape.ptr, &cs, 0);
				shape_data.list.push_back(cs);
				shape.compute_shape_index = index++;
			}
			else {
				shape.compute_shape_index = index;
				for (int f = 0; f != shape.animation_frames; ++f) {
					tfxGPUImageData cs;
					cs.animation_frames = shape.animation_frames;
					cs.image_size = shape.image_size;
					cs.uv = uv_lookup(shape.ptr, &cs, f);
					shape_data.list.push_back(cs);
					index++;
				}
			}
		}
		return shape_data;
	}

	void CopyLibraryLookupIndexesData(tfxLibrary *library, void* dst) {
		assert(dst);	//must be a valid pointer to a space in memory
		assert(library->compiled_lookup_indexes.size());		//There is no data to copy, make sure a library has been loaded properly and it contains effects with emitters
		tfxGraphLookupIndex *test = static_cast<tfxGraphLookupIndex*>(dst);
		memcpy(dst, library->compiled_lookup_indexes.data, GetLibraryLookupIndexesSizeInBytes(library));
	}

	void CopyLibraryLookupValuesData(tfxLibrary *library, void* dst) {
		assert(dst);	//must be a valid pointer to a space in memory
		assert(library->compiled_lookup_indexes.size());		//There is no data to copy, make sure a library has been loaded properly and it contains effects with emitters
		memcpy(dst, library->compiled_lookup_values.data, GetLibraryLookupValuesSizeInBytes(library));
	}

	tfxU32 GetLibraryComputeShapeDataSizeInBytes(tfxLibrary *library) {
		tfxU32 frame_count = 0;
		for (auto &shape : library->particle_shapes.data) {
			frame_count += (tfxU32)shape.animation_frames;
		}
		return frame_count * sizeof(tfxGPUImageData);
	}

	tfxU32 GetLibraryComputeShapeCount(tfxLibrary *library) {
		tfxU32 frame_count = 0;
		for (auto &shape : library->particle_shapes.data) {
			frame_count += (tfxU32)shape.animation_frames;
		}
		return frame_count;
	}

	tfxU32 GetLibraryLookupIndexCount(tfxLibrary *library) {
		return library->compiled_lookup_indexes.size() * tfxOvertimeCount;
	}

	tfxU32 GetLibraryLookupValueCount(tfxLibrary *library) {
		return library->compiled_lookup_values.size();
	}

	tfxU32 GetLibraryLookupIndexesSizeInBytes(tfxLibrary *library) {
		return sizeof(tfxGraphLookupIndex) * tfxOvertimeCount * library->compiled_lookup_indexes.size();
	}

	tfxU32 GetLibraryLookupValuesSizeInBytes(tfxLibrary *library) {
		return sizeof(float) * library->compiled_lookup_values.size();
	}

	bool IsLibraryShapeUsed(tfxLibrary *library, tfxKey image_hash) {
		tmpStack(tfxEffectEmitter*, effect_stack);
		for (auto &effect : library->effects) {
			effect_stack.push_back(&effect);
		}
		while (effect_stack.size()) {
			auto current = effect_stack.pop_back();
			if (current->type == tfxEmitterType) {
				if (GetEffectProperties(current)->image[current->property_index]->image_hash == image_hash) {
					return true;
				}
			}
			for (auto &sub : GetEffectInfo(current)->sub_effectors) {
				effect_stack.push_back(&sub);
			}
		}
		return false;
	}

	bool LibraryShapeExists(tfxLibrary *library, tfxKey image_hash) {
		return library->particle_shapes.ValidKey(image_hash);
	}

	bool RemoveLibraryShape(tfxLibrary *library, tfxKey image_hash) {
		if (!library->particle_shapes.ValidKey(image_hash)) {
			return false;
		}
		library->particle_shapes.Remove(image_hash);
		for (auto &m : library->particle_shapes.map) {
			library->particle_shapes[m.index].image_hash = m.key;
		}
		return true;
	}

	void DeleteLibraryEffect(tfxLibrary *library, tfxEffectEmitter *effect) {
		CleanUpEffect(&library->effects[effect->library_index]);
		library->effects.erase(&library->effects[effect->library_index]);

		UpdateLibraryEffectPaths(library);
		ReIndexLibrary(library);
	}

	tfxU32 AddLibraryGlobal(tfxLibrary *library) {
		if (library->free_global_graphs.size())
			return library->free_global_graphs.pop_back();
		tfxGlobalAttributes global;
		InitialiseGlobalAttributes(&global, &library->graph_node_allocator, &library->graph_lookup_allocator);
		library->global_graphs.push_back(global);
		return library->global_graphs.size() - 1;
	}

	tfxU32 AddLibraryKeyframes(tfxLibrary *library) {
		if (library->free_keyframes.size())
			return library->free_keyframes.pop_back();
		tfxTransformAttributes keyframes;
		InitialiseTransformAttributes(&keyframes, &library->graph_node_allocator, &library->graph_lookup_allocator);
		library->transform_attributes.push_back(keyframes);
		return library->transform_attributes.size() - 1;
	}

	tfxU32 AddLibraryEmitterAttributes(tfxLibrary *library) {
		if (library->free_emitter_attributes.size())
			return library->free_emitter_attributes.pop_back();
		tfxEmitterAttributes attributes;
		InitialiseEmitterAttributes(&attributes, &library->graph_node_allocator, &library->graph_lookup_allocator);
		library->emitter_attributes.push_back(attributes);
		return library->emitter_attributes.size() - 1;
	}

	void FreeLibraryGlobal(tfxLibrary *library, tfxU32 index) {
		assert(index < library->global_graphs.size());
		library->free_global_graphs.push_back(index);
		FreeGlobalAttributes(&library->global_graphs[index]);
	}

	void FreeLibraryKeyframes(tfxLibrary *library, tfxU32 index) {
		assert(index < library->transform_attributes.size());
		library->free_keyframes.push_back(index);
		FreeTransformAttributes(&library->transform_attributes[index]);
	}

	void FreeLibraryEmitterAttributes(tfxLibrary *library, tfxU32 index) {
		assert(index < library->emitter_attributes.size());
		library->free_emitter_attributes.push_back(index);
		FreeEmitterAttributes(&library->emitter_attributes[index]);
	}

	void FreeLibraryProperties(tfxLibrary *library, tfxU32 index) {
		assert(index < library->emitter_properties_buffer.current_size);
		library->free_properties.push_back(index);
	}

	void FreeLibraryInfo(tfxLibrary *library, tfxU32 index) {
		assert(index < library->effect_infos.size());
		library->free_infos.push_back(index);
	}

	tfxU32 CountLibraryKeyframeLookUpValues(tfxLibrary *library, tfxU32 index) {
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

	tfxU32 CountLibraryGlobalLookUpValues(tfxLibrary *library, tfxU32 index) {
		auto &global = library->global_graphs[index];
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

	tfxU32 CountLibraryEmitterLookUpValues(tfxLibrary *library, tfxU32 index) {
		auto &attributes = library->emitter_attributes[index];
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

	tfxU32 CloneLibraryGlobal(tfxLibrary *library, tfxU32 source_index, tfxLibrary *destination_library) {
		tfxU32 index = AddLibraryGlobal(destination_library);
		CopyGlobalAttributesNoLookups(&library->global_graphs[source_index], &destination_library->global_graphs[index]);
		return index;
	}

	tfxU32 CloneLibraryKeyframes(tfxLibrary *library, tfxU32 source_index, tfxLibrary *destination_library) {
		tfxU32 index = AddLibraryKeyframes(destination_library);
		CopyTransformAttributesNoLookups(&library->transform_attributes[source_index], &destination_library->transform_attributes[index]);
		return index;
	}

	tfxU32 CloneLibraryEmitterAttributes(tfxLibrary *library, tfxU32 source_index, tfxLibrary *destination_library) {
		tfxU32 index = AddLibraryEmitterAttributes(destination_library);
		CopyPropertyAttributesNoLookups(&library->emitter_attributes[source_index].properties, &destination_library->emitter_attributes[index].properties);
		CopyBaseAttributesNoLookups(&library->emitter_attributes[source_index].base, &destination_library->emitter_attributes[index].base);
		CopyVariationAttributesNoLookups(&library->emitter_attributes[source_index].variation, &destination_library->emitter_attributes[index].variation);
		CopyOvertimeAttributesNoLookups(&library->emitter_attributes[source_index].overtime, &destination_library->emitter_attributes[index].overtime);
		return index;
	}

	tfxU32 CloneLibraryInfo(tfxLibrary *library, tfxU32 source_index, tfxLibrary *destination_library) {
		tfxU32 index = AddLibraryEffectEmitterInfo(destination_library);
		destination_library->effect_infos[index] = library->effect_infos[source_index];
		return index;
	}

	tfxU32 CloneLibraryProperties(tfxLibrary *library, tfxU32 source_index, tfxLibrary *destination_library) {
		tfxU32 index = AddLibraryEmitterProperties(destination_library);
		CopyEmitterProperites(library->emitter_properties, source_index, destination_library->emitter_properties, index);
		return index;
	}

	void AddLibraryEmitterGraphs(tfxLibrary *library, tfxEffectEmitter& emitter) {
		emitter.emitter_attributes = AddLibraryEmitterAttributes(library);
	}

	void AddLibraryTransformGraphs(tfxLibrary *library, tfxEffectEmitter& emitter) {
		emitter.transform_attributes = AddLibraryKeyframes(library);
	}

	void AddLibraryEffectGraphs(tfxLibrary *library, tfxEffectEmitter& effect) {
		tfxEffectEmitter *root_effect = GetRootEffect(&effect);
		if (root_effect == &effect)
			effect.global = AddLibraryGlobal(library);
		else
			effect.global = root_effect->global;
	}

	tfxU32 AddLibrarySpriteSheetSettings(tfxLibrary *library, tfxEffectEmitter& effect) {
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
		library->sprite_sheet_settings.push_back(a);
		GetEffectInfo(&effect)->sprite_sheet_settings_index = library->sprite_sheet_settings.size() - 1;
		return GetEffectInfo(&effect)->sprite_sheet_settings_index;
	}

	void AddLibrarySpriteSheetSettingsSub(tfxLibrary *library, tfxEffectEmitter& effect) {
		if (effect.type == tfxEffectType) {
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
			library->sprite_sheet_settings.push_back(a);
			GetEffectInfo(&effect)->sprite_sheet_settings_index = library->sprite_sheet_settings.size() - 1;
		}
		else {
			for (auto &sub : GetEffectInfo(&effect)->sub_effectors) {
				AddLibrarySpriteSheetSettingsSub(effect.library, sub);
			}
		}
	}

	tfxU32 AddLibrarySpriteDataSettings(tfxLibrary *library, tfxEffectEmitter& effect) {
		assert(effect.type == tfxEffectType);
		tfxSpriteDataSettings a;
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
		//a.camera_settings.camera_fov = tfxRadians(60);
		//a.camera_settings.camera_pitch = tfxRadians(-30.f);
		//a.camera_settings.camera_yaw = tfxRadians(-90.f);
		//a.camera_settings.camera_position = tfxVec3(0.f, 3.5f, 7.5f);
		//a.camera_settings.camera_isometric = false;
		//a.camera_settings.camera_isometric_scale = 5.f;
		//a.camera_settings.camera_hide_floor = false;
		library->sprite_data_settings.push_back(a);
		GetEffectInfo(&effect)->sprite_data_settings_index = library->sprite_data_settings.size() - 1;
		return GetEffectInfo(&effect)->sprite_data_settings_index;
	}

	void AddLibrarySpriteDataSettingsSub(tfxLibrary *library, tfxEffectEmitter& effect) {
		if (effect.type == tfxEffectType) {
			tfxSpriteDataSettings a;
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
			GetEffectInfo(&effect)->sprite_data_settings_index = library->sprite_data_settings.size() - 1;
		}
		else {
			for (auto &sub : GetEffectInfo(&effect)->sub_effectors) {
				AddLibrarySpriteDataSettingsSub(effect.library, sub);
			}
		}
	}

	tfxU32 AddLibraryPreviewCameraSettings(tfxLibrary *library, tfxEffectEmitter& effect) {
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
		library->preview_camera_settings.push_back(a);
		GetEffectInfo(&effect)->preview_camera_settings = library->preview_camera_settings.size() - 1;
		return GetEffectInfo(&effect)->preview_camera_settings;
	}

	tfxU32 AddLibraryPreviewCameraSettings(tfxLibrary *library) {
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
		library->preview_camera_settings.push_back(a);
		return library->preview_camera_settings.size() - 1;
	}

	tfxU32 AddLibraryEffectEmitterInfo(tfxLibrary *library) {
		tfxEffectEmitterInfo info;
		if (library->free_infos.size()) {
			return library->free_infos.pop_back();
		}
		library->effect_infos.push_back(info);
		return library->effect_infos.size() - 1;
	}

	tfxU32 AddLibraryEmitterProperties(tfxLibrary *library) {
		if (library->free_properties.size()) {
			return library->free_properties.pop_back();
		}
		return AddRow(&library->emitter_properties_buffer, true);
	}

	void InitLibrary(tfxLibrary *library) {
		library->graph_node_allocator = CreateArenaManager(tfxMegabyte(2), 8);
		library->graph_lookup_allocator = CreateArenaManager(tfxMegabyte(4), 8);
		InitLibraryEmitterProperties(library);
	}

	void InitLibraryEmitterProperties(tfxLibrary *library) {
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxVec3), offsetof(tfxEmitterPropertiesSoA, angle_offsets));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxVectorAlignType), offsetof(tfxEmitterPropertiesSoA, vector_align_type));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxEmissionType), offsetof(tfxEmitterPropertiesSoA, emission_type));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxU32), offsetof(tfxEmitterPropertiesSoA, single_shot_limit));
		AddStructArray(&library->emitter_properties_buffer, sizeof(float), offsetof(tfxEmitterPropertiesSoA, frame_rate));
		AddStructArray(&library->emitter_properties_buffer, sizeof(float), offsetof(tfxEmitterPropertiesSoA, end_frame));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxImageData*), offsetof(tfxEmitterPropertiesSoA, image));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxBillboardingOptions), offsetof(tfxEmitterPropertiesSoA, billboard_option));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxVec3), offsetof(tfxEmitterPropertiesSoA, grid_points));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxAngleSettingFlags), offsetof(tfxEmitterPropertiesSoA, angle_settings));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxU32), offsetof(tfxEmitterPropertiesSoA, layer));
		AddStructArray(&library->emitter_properties_buffer, sizeof(float), offsetof(tfxEmitterPropertiesSoA, delay_spawning));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxEmissionDirection), offsetof(tfxEmitterPropertiesSoA, emission_direction));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxLineTraversalEndBehaviour), offsetof(tfxEmitterPropertiesSoA, end_behaviour));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxParticleControlFlags), offsetof(tfxEmitterPropertiesSoA, compute_flags));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxVec2), offsetof(tfxEmitterPropertiesSoA, image_handle));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxVec3), offsetof(tfxEmitterPropertiesSoA, emitter_handle));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxU32), offsetof(tfxEmitterPropertiesSoA, spawn_amount));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxU32), offsetof(tfxEmitterPropertiesSoA, image_index));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxKey), offsetof(tfxEmitterPropertiesSoA, image_hash));
		AddStructArray(&library->emitter_properties_buffer, sizeof(float), offsetof(tfxEmitterPropertiesSoA, loop_length));
		AddStructArray(&library->emitter_properties_buffer, sizeof(float), offsetof(tfxEmitterPropertiesSoA, start_frame));
		AddStructArray(&library->emitter_properties_buffer, sizeof(float), offsetof(tfxEmitterPropertiesSoA, noise_base_offset_range));
		AddStructArray(&library->emitter_properties_buffer, sizeof(tfxU32), offsetof(tfxEmitterPropertiesSoA, animation_property_index));
		FinishSoABufferSetup(&library->emitter_properties_buffer, &library->emitter_properties, 100);
	}

	void ClearLibrary(tfxLibrary *library) {
		for (auto &e : library->effects) {
			FreeEffectGraphs(&e);
		}
		library->effects.free_all();
		library->effect_paths.FreeAll();
		library->particle_shapes.FreeAll();
		library->global_graphs.free_all();
		library->emitter_attributes.free_all();
		library->transform_attributes.free_all();
		library->sprite_sheet_settings.free_all();
		library->preview_camera_settings.free_all();
		library->all_nodes.free_all();
		library->node_lookup_indexes.free_all();
		library->compiled_lookup_indexes.free_all();
		library->compiled_lookup_values.free_all();
		FreeSoABuffer(&library->emitter_properties_buffer);
		library->graph_min_max.free_all();
		for (auto &info : library->effect_infos) {
			info.sub_effectors.free_all();
		}
		library->effect_infos.free_all();
		AddLibraryPreviewCameraSettings(library);
		library->pre_recorded_effects.FreeAll();

		library->graph_node_allocator.FreeAll();
		library->graph_lookup_allocator.FreeAll();

		library->free_global_graphs.free_all();
		library->free_keyframe_graphs.free_all();
		library->free_emitter_attributes.free_all();
		library->free_animation_settings.free_all();
		library->free_preview_camera_settings.free_all();
		library->free_infos.free_all();
		library->free_properties.free_all();
		library->free_keyframes.free_all();

		library->uid = 0;
	}

	void UpdateLibraryComputeNodes(tfxLibrary *library) {
		tfxU32 running_node_index = 0;
		tfxU32 running_value_index = 0;
		tmpStack(tfxEffectEmitter*, stack);
		library->all_nodes.clear();
		library->node_lookup_indexes.clear();
		library->compiled_lookup_values.clear();
		library->compiled_lookup_indexes.clear();
		for (auto &effect : library->effects) {
			stack.push_back(&effect);
			while (!stack.empty()) {
				tfxEffectEmitter *current = stack.pop_back();
				if (current->type == tfxFolder) {
					for (auto &sub : GetEffectInfo(current)->sub_effectors) {
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

					GetEffectInfo(current)->lookup_value_index = library->compiled_lookup_indexes.size();
					for (int i = 0; i != tfxOvertimeCount; ++i) {
						tfxGraph &graph = ((tfxGraph*)(&library->emitter_attributes[current->emitter_attributes].overtime))[i];
						tfxGraphLookupIndex &index = ((tfxGraphLookupIndex*)&lookup_data)[i];
						index.start_index = running_node_index;
						index.length = graph.nodes.size();
						index.max_life = graph.lookup.life;
						graph.nodes.ResetIteratorIndex();
						do {
							for (auto &node : graph.nodes) {
								library->all_nodes.push_back(node);
								running_node_index++;
							}
						} while (!graph.nodes.EndOfBuckets());

						tfxGraphLookupIndex value_index;
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

	void CompileLibraryGraphsOfEffect(tfxLibrary *library, tfxEffectEmitter &effect, tfxU32 depth) {
		tfxEffectEmitterInfo *info = GetEffectInfo(&effect);
		if (effect.type == tfxEffectType && depth == 0) {
			CompileLibraryKeyframeGraph(library, effect.transform_attributes);
			CompileLibraryGlobalGraph(library, effect.global);
		}
		else if (effect.type == tfxEffectType) {
			CompileLibraryKeyframeGraph(library, effect.transform_attributes);
		}
		else if (effect.type == tfxEmitterType) {
			CompileLibraryKeyframeGraph(library, effect.transform_attributes);
			CompileLibraryEmitterGraphs(library, effect.emitter_attributes);
			for (auto &sub : info->sub_effectors) {
				CompileLibraryGraphsOfEffect(library, sub, ++depth);
			}
		}
		else if (effect.type == tfxFolder) {
			for (auto &sub : info->sub_effectors) {
				CompileLibraryGraphsOfEffect(library, sub, 0);
			}
		}
		if (effect.type == tfxEffectType) {
			for (auto &sub : info->sub_effectors) {
				CompileLibraryGraphsOfEffect(library, sub, ++depth);
			}
		}
	}

	void CompileAllLibraryGraphs(tfxLibrary *library) {
		for (auto &g : library->global_graphs) {
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
		for (auto &g : library->transform_attributes) {
			CompileGraph(g.roll);
			CompileGraph(g.pitch);
			CompileGraph(g.yaw);
			CompileGraph(g.translation_x);
			CompileGraph(g.translation_y);
			CompileGraph(g.translation_z);
		}
		for (auto &g : library->emitter_attributes) {
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

			CompileColorOvertime(g.overtime.red);
			CompileColorOvertime(g.overtime.green);
			CompileColorOvertime(g.overtime.blue);
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

	void CompileLibraryGlobalGraph(tfxLibrary *library, tfxU32 index) {
		tfxGlobalAttributes &g = library->global_graphs[index];
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

	void CompileLibraryKeyframeGraph(tfxLibrary *library, tfxU32 index) {
		tfxTransformAttributes &g = library->transform_attributes[index];
		CompileGraph(g.roll);
		CompileGraph(g.pitch);
		CompileGraph(g.yaw);
		CompileGraph(g.translation_x);
		CompileGraph(g.translation_y);
		CompileGraph(g.translation_z);
	}

	void CompileLibraryEmitterGraphs(tfxLibrary *library, tfxU32 index) {
		CompileLibraryPropertyGraph(library, index);
		CompileLibraryKeyframeGraph(library, index);
		CompileLibraryBaseGraph(library, index);
		CompileLibraryVariationGraph(library, index);
		CompileLibraryOvertimeGraph(library, index);
	}

	void CompileLibraryPropertyGraph(tfxLibrary *library, tfxU32 index) {
		tfxPropertyAttributes &g = library->emitter_attributes[index].properties;
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
	void CompileLibraryBaseGraph(tfxLibrary *library, tfxU32 index) {
		tfxBaseAttributes &g = library->emitter_attributes[index].base;
		CompileGraph(g.amount);
		CompileGraph(g.width);
		CompileGraph(g.height);
		CompileGraph(g.life);
		CompileGraph(g.spin);
		CompileGraph(g.noise_offset);
		CompileGraph(g.velocity);
		CompileGraph(g.weight);
	}
	void CompileLibraryVariationGraph(tfxLibrary *library, tfxU32 index) {
		tfxVariationAttributes &g = library->emitter_attributes[index].variation;
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
	void CompileLibraryOvertimeGraph(tfxLibrary *library, tfxU32 index) {
		tfxOvertimeAttributes &g = library->emitter_attributes[index].overtime;
		CompileColorOvertime(g.red);
		CompileColorOvertime(g.green);
		CompileColorOvertime(g.blue);
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
	void CompileLibraryColorGraphs(tfxLibrary *library, tfxU32 index) {
		tfxOvertimeAttributes &g = library->emitter_attributes[index].overtime;
		CompileColorOvertime(g.red);
		CompileColorOvertime(g.green);
		CompileColorOvertime(g.blue);
	}

	void SetLibraryMinMaxData(tfxLibrary *library) {
		library->graph_min_max.clear();
		library->graph_min_max.create_pool(tfxGraphMaxIndex);

		library->graph_min_max[tfxGlobal_life] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		library->graph_min_max[tfxGlobal_amount] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		library->graph_min_max[tfxGlobal_velocity] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		library->graph_min_max[tfxGlobal_width] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		library->graph_min_max[tfxGlobal_height] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		library->graph_min_max[tfxGlobal_weight] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		library->graph_min_max[tfxGlobal_spin] = GetMinMaxGraphValues(tfxGlobalPercentPresetSigned);
		library->graph_min_max[tfxGlobal_stretch] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		library->graph_min_max[tfxGlobal_overal_scale] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		library->graph_min_max[tfxGlobal_intensity] = GetMinMaxGraphValues(tfxOpacityOvertimePreset);
		library->graph_min_max[tfxGlobal_frame_rate] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
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
		library->graph_min_max[tfxProperty_arc_size] = GetMinMaxGraphValues(tfxArcPreset);
		library->graph_min_max[tfxProperty_arc_offset] = GetMinMaxGraphValues(tfxArcPreset);

		library->graph_min_max[tfxBase_life] = GetMinMaxGraphValues(tfxLifePreset);
		library->graph_min_max[tfxBase_amount] = GetMinMaxGraphValues(tfxAmountPreset);
		library->graph_min_max[tfxBase_velocity] = GetMinMaxGraphValues(tfxVelocityPreset);
		library->graph_min_max[tfxBase_width] = GetMinMaxGraphValues(tfxDimensionsPreset);
		library->graph_min_max[tfxBase_height] = GetMinMaxGraphValues(tfxDimensionsPreset);
		library->graph_min_max[tfxBase_weight] = GetMinMaxGraphValues(tfxWeightPreset);
		library->graph_min_max[tfxBase_spin] = GetMinMaxGraphValues(tfxSpinPreset);
		library->graph_min_max[tfxBase_noise_offset] = GetMinMaxGraphValues(tfxDimensionsPreset);

		library->graph_min_max[tfxVariation_life] = GetMinMaxGraphValues(tfxLifePreset);
		library->graph_min_max[tfxVariation_amount] = GetMinMaxGraphValues(tfxAmountPreset);
		library->graph_min_max[tfxVariation_velocity] = GetMinMaxGraphValues(tfxVelocityPreset);
		library->graph_min_max[tfxVariation_width] = GetMinMaxGraphValues(tfxDimensionsPreset);
		library->graph_min_max[tfxVariation_height] = GetMinMaxGraphValues(tfxDimensionsPreset);
		library->graph_min_max[tfxVariation_weight] = GetMinMaxGraphValues(tfxWeightVariationPreset);
		library->graph_min_max[tfxVariation_spin] = GetMinMaxGraphValues(tfxSpinVariationPreset);
		library->graph_min_max[tfxVariation_noise_offset] = GetMinMaxGraphValues(tfxGlobalPercentPreset);

		library->graph_min_max[tfxOvertime_velocity] = GetMinMaxGraphValues(tfxVelocityOvertimePreset);
		library->graph_min_max[tfxOvertime_width] = GetMinMaxGraphValues(tfxPercentOvertime);
		library->graph_min_max[tfxOvertime_height] = GetMinMaxGraphValues(tfxPercentOvertime);
		library->graph_min_max[tfxOvertime_weight] = GetMinMaxGraphValues(tfxPercentOvertime);
		library->graph_min_max[tfxOvertime_spin] = GetMinMaxGraphValues(tfxSpinOvertimePreset);
		library->graph_min_max[tfxOvertime_stretch] = GetMinMaxGraphValues(tfxPercentOvertime);
		library->graph_min_max[tfxOvertime_red] = GetMinMaxGraphValues(tfxColorPreset);
		library->graph_min_max[tfxOvertime_green] = GetMinMaxGraphValues(tfxColorPreset);
		library->graph_min_max[tfxOvertime_blue] = GetMinMaxGraphValues(tfxColorPreset);
		library->graph_min_max[tfxOvertime_blendfactor] = GetMinMaxGraphValues(tfxOpacityOvertimePreset);
		library->graph_min_max[tfxOvertime_intensity] = GetMinMaxGraphValues(tfxIntensityOvertimePreset);
		library->graph_min_max[tfxOvertime_velocity_turbulance] = GetMinMaxGraphValues(tfxFrameratePreset);
		library->graph_min_max[tfxOvertime_direction_turbulance] = GetMinMaxGraphValues(tfxPercentOvertime);
		library->graph_min_max[tfxOvertime_velocity_adjuster] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		library->graph_min_max[tfxOvertime_direction] = GetMinMaxGraphValues(tfxDirectionOvertimePreset);
	}

	float LookupLibraryPreciseOvertimeNodeList(tfxLibrary *library, tfxGraphType graph_type, int lookup_node_index, float age, float life) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		tfxAttributeNode *lastec = nullptr;
		tfxGraphLookupIndex &lookup_data = ((tfxGraphLookupIndex*)&library->node_lookup_indexes[lookup_node_index])[graph_type];
		float min_y = library->graph_min_max[graph_type].y;
		float max_y = library->graph_min_max[graph_type].w;
		for (int i = lookup_data.start_index; i != lookup_data.start_index + lookup_data.length; ++i) {
			tfxAttributeNode &a = library->all_nodes[i];
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

	float LookupLibraryPreciseNodeList(tfxLibrary *library, tfxGraphType graph_type, int lookup_node_index, float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		tfxAttributeNode *lastec = nullptr;
		tfxGraphLookupIndex &lookup_data = ((tfxGraphLookupIndex*)&library->node_lookup_indexes[lookup_node_index])[graph_type];
		float min_y = library->graph_min_max[graph_type].y;
		float max_y = library->graph_min_max[graph_type].w;
		for (int i = lookup_data.start_index; i != lookup_data.start_index + lookup_data.length; ++i) {
			tfxAttributeNode &a = library->all_nodes[i];
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

	float LookupLibraryFastValueList(tfxLibrary *library, tfxGraphType graph_type, int lookup_node_index, float frame) {
		tfxGraphLookupIndex &lookup_data = ((tfxGraphLookupIndex*)&library->compiled_lookup_indexes[lookup_node_index])[graph_type];
		frame += lookup_data.start_index;
		tfxU32 end_frame = lookup_data.start_index + lookup_data.length - 1;
		frame = frame > end_frame ? end_frame : frame;
		return library->compiled_lookup_values[(tfxU32)frame];
	}

	float LookupLibraryFastOvertimeValueList(tfxLibrary *library, tfxGraphType graph_type, int lookup_value_index, float age, float lifetime) {
		tfxGraphLookupIndex &lookup_data = ((tfxGraphLookupIndex*)&library->compiled_lookup_indexes[lookup_value_index])[graph_type - tfxOvertime_velocity];
		float frame = (float)lookup_data.start_index;
		if (lifetime)
			frame += (age / lifetime * lookup_data.max_life) / tfxLOOKUP_FREQUENCY_OVERTIME;
		if (frame < lookup_data.start_index + lookup_data.length - 1)
			return library->compiled_lookup_values[(tfxU32)frame];
		return library->compiled_lookup_values[lookup_data.start_index + lookup_data.length - 1];
	}

	tfxU32 CountOfGraphsInUse(tfxLibrary *library) {
		return library->global_graphs.size() + library->emitter_attributes.size() - CountOfFreeGraphs(library);
	}

	tfxU32 CountOfFreeGraphs(tfxLibrary *library) {
		return library->free_global_graphs.size() + library->free_emitter_attributes.size();
	}

	void tfxDataTypesDictionary::Init() {
		names_and_types.data.reserve(200);
		names_and_types.map.reserve(200);
		names_and_types.Insert("name", tfxString);
		names_and_types.Insert("image_index", tfxUint);
		names_and_types.Insert("image_hash", tfxUInt64);
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
		names_and_types.Insert("corner1_3d", tfxFloat3);
		names_and_types.Insert("corner2_3d", tfxFloat3);

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
		names_and_types.Insert("background_image_tint_red", tfxFloat);
		names_and_types.Insert("background_image_tint_green", tfxFloat);
		names_and_types.Insert("background_image_tint_blue", tfxFloat);
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
		names_and_types.Insert("sprite_data_tab_use_compute", tfxSInt);
		names_and_types.Insert("update_mode", tfxSInt);
		names_and_types.Insert("camera_mouse_sensitivity", tfxFloat);
		initialised = true;
	}

	int ValidateEffectPackage(const char *filename) {
		tfxPackage package;
		tfxErrorFlags status = LoadPackage(filename, package);
		if (status) {
			FreePackage(&package);
			return status;					//returns 1 to 4 if it's an invalid package format
		}

		tfxEntryInfo *data_txt = GetPackageFile(&package, "data.txt");
		FreePackage(&package);
		if (!data_txt) return tfxErrorCode_data_could_not_be_loaded;					//Unable to load the the data.txt file in the package

		return 0;
	}

	void AssignGraphData(tfxEffectEmitter &effect, tfxStack<tfxStr256> &values) {
		if (values.size() > 0) {
			if (values[0] == "global_amount") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].amount, n); }
			if (values[0] == "global_frame_rate") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].frame_rate, n); }
			if (values[0] == "global_height") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].height, n); }
			if (values[0] == "global_width") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].width, n); }
			if (values[0] == "global_life") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].life, n); }
			if (values[0] == "global_opacity") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].intensity, n); }
			if (values[0] == "global_spin") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].spin, n); }
			if (values[0] == "global_splatter") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].splatter, n); }
			if (values[0] == "global_stretch") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].stretch, n); }
			if (values[0] == "global_overal_scale") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].overal_scale, n); }
			if (values[0] == "global_weight") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].weight, n); }
			if (values[0] == "global_velocity") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].velocity, n); }
			if (values[0] == "global_emitter_width") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].emitter_width, n); }
			if (values[0] == "global_emitter_height") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].emitter_height, n); }
			if (values[0] == "global_emitter_depth") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->global_graphs[effect.global].emitter_depth, n); }

			if (values[0] == "global_effect_angle") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->transform_attributes[effect.transform_attributes].roll, n); }
			if (values[0] == "global_effect_roll") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->transform_attributes[effect.transform_attributes].roll, n); }
			if (values[0] == "global_effect_pitch") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->transform_attributes[effect.transform_attributes].pitch, n); }
			if (values[0] == "global_effect_yaw") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->transform_attributes[effect.transform_attributes].yaw, n); }
			if (values[0] == "keyframe_translate_x") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->transform_attributes[effect.transform_attributes].translation_x, n); }
			if (values[0] == "keyframe_translate_y") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->transform_attributes[effect.transform_attributes].translation_y, n); }
			if (values[0] == "keyframe_translate_z") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->transform_attributes[effect.transform_attributes].translation_z, n); }
			if (values[0] == "property_emitter_angle") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->transform_attributes[effect.transform_attributes].roll, n); }
			if (values[0] == "property_emitter_roll") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->transform_attributes[effect.transform_attributes].roll, n); }
			if (values[0] == "property_emitter_pitch") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->transform_attributes[effect.transform_attributes].pitch, n); }
			if (values[0] == "property_emitter_yaw") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->transform_attributes[effect.transform_attributes].yaw, n); }

			if (values[0] == "base_arc_offset") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.arc_offset, n); }
			if (values[0] == "base_arc_size") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.arc_size, n); }
			if (values[0] == "base_emission_angle") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.emission_pitch, n); }
			if (values[0] == "base_emission_range") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.emission_range, n); }
			if (values[0] == "base_emitter_height") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.emitter_height, n); }
			if (values[0] == "base_emitter_width") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.emitter_width, n); }
			if (values[0] == "base_splatter") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.splatter, n); }

			if (values[0] == "property_arc_offset") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.arc_offset, n); }
			if (values[0] == "property_arc_size") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.arc_size, n); }
			if (values[0] == "property_emission_angle") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.emission_pitch, n); }
			if (values[0] == "property_emission_pitch") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.emission_pitch, n); }
			if (values[0] == "property_emission_yaw") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.emission_yaw, n); }
			if (values[0] == "property_emission_range") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.emission_range, n); }
			if (values[0] == "property_emitter_height") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.emitter_height, n); }
			if (values[0] == "property_emitter_width") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.emitter_width, n); }
			if (values[0] == "property_emitter_depth") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.emitter_depth, n); }
			if (values[0] == "property_splatter") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].properties.splatter, n); }

			if (values[0] == "base_amount") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].base.amount, n); }
			if (values[0] == "base_life") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].base.life, n); }
			if (values[0] == "base_height") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].base.height, n); }
			if (values[0] == "base_width") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].base.width, n); }
			if (values[0] == "base_spin") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].base.spin, n); }
			if (values[0] == "base_noise_offset") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].base.noise_offset, n); }
			if (values[0] == "base_velocity") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].base.velocity, n); }
			if (values[0] == "base_weight") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].base.weight, n); }

			if (values[0] == "variation_amount") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].variation.amount, n); }
			if (values[0] == "variation_height") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].variation.height, n); }
			if (values[0] == "variation_width") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].variation.width, n); }
			if (values[0] == "variation_life") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].variation.life, n); }
			if (values[0] == "variation_velocity") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].variation.velocity, n); }
			if (values[0] == "variation_weight") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].variation.weight, n); }
			if (values[0] == "variation_spin") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].variation.spin, n); }
			if (values[0] == "variation_motion_randomness") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].variation.noise_offset, n); }
			if (values[0] == "variation_noise_offset") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].variation.noise_offset, n); }
			if (values[0] == "variation_noise_resolution") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].variation.noise_resolution, n); }

			if (values[0] == "overtime_red") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.red, n); }
			if (values[0] == "overtime_green") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.green, n); }
			if (values[0] == "overtime_blue") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.blue, n); }
			if (values[0] == "overtime_opacity") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.blendfactor, n); }
			if (values[0] == "overtime_intensity") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.intensity, n); }
			if (values[0] == "overtime_velocity_turbulance") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.velocity_turbulance, n); }
			if (values[0] == "overtime_spin") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.spin, n); }
			if (values[0] == "overtime_stretch") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.stretch, n); }
			if (values[0] == "overtime_velocity") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.velocity, n); }
			if (values[0] == "overtime_weight") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.weight, n); }
			if (values[0] == "overtime_width") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.width, n); }
			if (values[0] == "overtime_height") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.height, n); }
			if (values[0] == "overtime_direction_turbulance") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.direction_turbulance, n); }
			if (values[0] == "overtime_direction") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.direction, n); }
			if (values[0] == "overtime_velocity_adjuster") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.velocity_adjuster, n); }
			if (values[0] == "overtime_noise_resolution") { tfxAttributeNode n; AssignNodeData(n, values); AddGraphNode(&effect.library->emitter_attributes[effect.emitter_attributes].overtime.noise_resolution, n); }
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

	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, tfxU32 value) {
	}

	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, float value) {
	}

	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, bool value) {
	}

	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, int value) {
	}

	void AssignStageProperty(tfxEffectEmitter &effect, tfxStr &field, tfxStr &value) {
		if (field == "name") {
			GetEffectInfo(&effect)->name = value;
		}
	}

	void AssignFrameMetaProperty(tfxFrameMeta &metrics, tfxStr &field, tfxU32 value, tfxU32 file_version) {
		if (field == "total_sprites")
			metrics.total_sprites = value;
	}

	tfxVec3 StrToVec3(tfxStack<tfxStr256> &str) {
		assert(str.size() == 3);	//array must be size 3
		return tfxVec3((float)atof(str[0].c_str()), (float)atof(str[1].c_str()), (float)atof(str[2].c_str()));
	}

	tfxVec2 StrToVec2(tfxStack<tfxStr256> &str) {
		assert(str.size() == 2);	//array must be size 2
		return tfxVec2((float)atof(str[0].c_str()), (float)atof(str[1].c_str()));
	}

	void AssignFrameMetaProperty(tfxFrameMeta &metrics, tfxStr &field, tfxVec3 value, tfxU32 file_version) {
		if (field == "corner1")
			metrics.corner1 = value;
		if (field == "corner2")
			metrics.corner2 = value;
	}

	void AssignAnimationEmitterProperty(tfxAnimationEmitterProperties &properties, tfxStr &field, tfxVec2 value, tfxU32 file_version) {
		if (field == "handle")
			properties.handle = value;
	}

	void AssignAnimationEmitterProperty(tfxAnimationEmitterProperties &properties, tfxStr &field, float value, tfxU32 file_version) {
		if (field == "animation_frames")
			properties.animation_frames = value;
	}

	void AssignAnimationEmitterProperty(tfxAnimationEmitterProperties &properties, tfxStr &field, tfxU32 value, tfxU32 file_version) {
		if (field == "flags")
			properties.flags = value;
		if (field == "start_frame_index")
			properties.start_frame_index = value;
	}

	void AssignSpriteDataMetricsProperty(tfxSpriteDataMetrics &metrics, tfxStr &field, tfxU32 value, tfxU32 file_version) {
		if (field == "start_offset")
			metrics.start_offset = value;
		if (field == "frames_after_compression")
			metrics.frames_after_compression = value;
		if (field == "real_frames")
			metrics.real_frames = value;
		if (field == "frame_count")
			metrics.frame_count = value;
		if (field == "total_sprites")
			metrics.total_sprites = value;
		if (field == "total_memory_for_sprites")
			metrics.total_memory_for_sprites = value;
		if (field == "flags")
			metrics.flags = value;
		if (field == "animation_flags")
			metrics.animation_flags = value;
	}

	void AssignSpriteDataMetricsProperty(tfxSpriteDataMetrics &metrics, tfxStr &field, tfxU64 value, tfxU32 file_version) {
		if (field == "path_hash")
			metrics.path_hash = value;
	}

	void AssignSpriteDataMetricsProperty(tfxSpriteDataMetrics &metrics, tfxStr &field, float value, tfxU32 file_version) {
		if (field == "animation_length_in_time")
			metrics.animation_length_in_time = value;
	}

	void AssignSpriteDataMetricsProperty(tfxSpriteDataMetrics &metrics, tfxStr &field, tfxStr value, tfxU32 file_version) {
		if (field == "name")
			metrics.name = value;
	}

	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, tfxU64 value, tfxU32 file_version) {
		tfxEmitterPropertiesSoA &emitter_properties = effect.library->emitter_properties;
		if (field == "image_hash")
			emitter_properties.image_hash[effect.property_index] = value;
	}

	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, tfxU32 value, tfxU32 file_version) {
		tfxEmitterPropertiesSoA &emitter_properties = effect.library->emitter_properties;
		if (field == "image_index")
			emitter_properties.image_index[effect.property_index] = value;
		if (field == "spawn_amount")
			emitter_properties.spawn_amount[effect.property_index] = value;
		if (field == "frames")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].frames = value;
		if (field == "current_frame")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].current_frame = value;
		if (field == "seed")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].seed = value;
		if (field == "layer")
			emitter_properties.layer[effect.property_index] = value >= tfxLAYERS ? value = tfxLAYERS - 1 : value;
		if (field == "frame_offset")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].frame_offset = value;
		if (field == "single_shot_limit")
			emitter_properties.single_shot_limit[effect.property_index] = value;
		if (field == "billboard_option") {
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
			emitter_properties.billboard_option[effect.property_index] = (tfxBillboardingOptions)value;
		}
		if (field == "vector_align_type")
			emitter_properties.vector_align_type[effect.property_index] = value >= 0 && value < tfxVectorAlignType_max ? (tfxVectorAlignType)value : (tfxVectorAlignType)0;
		if (field == "angle_setting")
			emitter_properties.angle_settings[effect.property_index] = (tfxAngleSettingFlags)value;
		if (field == "sort_passes")
			effect.sort_passes = tfxMin(5, value);
		if (field == "animation_flags")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].animation_flags = value | tfxAnimationFlags_needs_recording;
		if (field == "sprite_data_flags")
			effect.library->sprite_data_settings[GetEffectInfo(&effect)->sprite_data_settings_index].animation_flags = value | tfxAnimationFlags_needs_recording;
		if (field == "sprite_data_seed")
			effect.library->sprite_data_settings[GetEffectInfo(&effect)->sprite_data_settings_index].seed = value;
		if (field == "sprite_data_frame_offset")
			effect.library->sprite_data_settings[GetEffectInfo(&effect)->sprite_data_settings_index].frame_offset = value;
		if (field == "sprite_data_frames")
			effect.library->sprite_data_settings[GetEffectInfo(&effect)->sprite_data_settings_index].real_frames = value;
		if (field == "sprite_data_extra_frames_count")
			effect.library->sprite_data_settings[GetEffectInfo(&effect)->sprite_data_settings_index].extra_frames_count = value;
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, int value) {
		tfxEmitterPropertiesSoA &emitter_properties = effect.library->emitter_properties;
		if (field == "emission_type")
			emitter_properties.emission_type[effect.property_index] = (tfxEmissionType)value;
		if (field == "emission_direction")
			emitter_properties.emission_direction[effect.property_index] = (tfxEmissionDirection)value;
		if (field == "color_option")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].color_option = value > 3 ? tfxFullColor : (tfxExportColorOptions)value;
		if (field == "export_option")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].export_option = (tfxExportOptions)value;
		if (field == "end_behaviour")
			emitter_properties.end_behaviour[effect.property_index] = (tfxLineTraversalEndBehaviour)value;
		if (field == "frame_offset")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].frame_offset = value;
		if (field == "extra_frames_count")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].extra_frames_count = value;
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, tfxStr &value) {
		if (field == "name") {
			GetEffectInfo(&effect)->name = value;
		}
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, float value) {
		tfxEmitterPropertiesSoA &emitter_properties = effect.library->emitter_properties;
		if (field == "position_x")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].position.x = value;
		if (field == "position_y")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].position.y = value;
		if (field == "position_z")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].position.z = value;
		if (field == "frame_width")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].frame_size.x = value;
		if (field == "frame_height")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].frame_size.y = value;
		if (field == "zoom")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].zoom = value;
		if (field == "scale")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].scale = value;
		if (field == "playback_speed")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].playback_speed = value;
		if (field == "camera_position_x")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings.camera_position.x = value;
		if (field == "camera_position_y")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings.camera_position.y = value;
		if (field == "camera_position_z")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings.camera_position.z = value;
		if (field == "camera_pitch")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings.camera_pitch = value;
		if (field == "camera_yaw")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings.camera_yaw = value;
		if (field == "camera_fov")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings.camera_fov = value;
		if (field == "camera_floor_height")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings.camera_floor_height = value;
		if (field == "camera_isometric_scale")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings.camera_isometric_scale = value;
		if (field == "orthographic_camera_position_x")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_position.x = value;
		if (field == "orthographic_camera_position_y")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_position.y = value;
		if (field == "orthographic_camera_position_z")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_position.z = value;
		if (field == "orthographic_camera_pitch")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_pitch = value;
		if (field == "orthographic_camera_yaw")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_yaw = value;
		if (field == "orthographic_camera_fov")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_fov = value;
		if (field == "orthographic_camera_floor_height")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_floor_height = value;
		if (field == "orthographic_camera_isometric_scale")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_isometric_scale = value;
		if (field == "preview_camera_position_x")
			effect.library->preview_camera_settings[GetEffectInfo(&effect)->preview_camera_settings].camera_settings.camera_position.x = value;
		if (field == "preview_camera_position_y")
			effect.library->preview_camera_settings[GetEffectInfo(&effect)->preview_camera_settings].camera_settings.camera_position.y = value;
		if (field == "preview_camera_position_z")
			effect.library->preview_camera_settings[GetEffectInfo(&effect)->preview_camera_settings].camera_settings.camera_position.z = value;
		if (field == "preview_camera_pitch")
			effect.library->preview_camera_settings[GetEffectInfo(&effect)->preview_camera_settings].camera_settings.camera_pitch = value;
		if (field == "preview_camera_yaw")
			effect.library->preview_camera_settings[GetEffectInfo(&effect)->preview_camera_settings].camera_settings.camera_yaw = value;
		if (field == "preview_camera_fov")
			effect.library->preview_camera_settings[GetEffectInfo(&effect)->preview_camera_settings].camera_settings.camera_fov = value;
		if (field == "preview_camera_floor_height")
			effect.library->preview_camera_settings[GetEffectInfo(&effect)->preview_camera_settings].camera_settings.camera_floor_height = value;
		if (field == "preview_camera_isometric_scale")
			effect.library->preview_camera_settings[GetEffectInfo(&effect)->preview_camera_settings].camera_settings.camera_isometric_scale = value == 0 ? 5.f : value;
		if (field == "preview_effect_z_offset")
			effect.library->preview_camera_settings[GetEffectInfo(&effect)->preview_camera_settings].effect_z_offset = value;
		if (field == "preview_camera_speed")
			effect.library->preview_camera_settings[GetEffectInfo(&effect)->preview_camera_settings].camera_speed = value;
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
			effect.library->sprite_data_settings[GetEffectInfo(&effect)->sprite_data_settings_index].playback_speed = value;
		if (field == "sprite_data_recording_frame_rate")
			effect.library->sprite_data_settings[GetEffectInfo(&effect)->sprite_data_settings_index].recording_frame_rate = value;
	}
	void AssignEffectorProperty(tfxEffectEmitter &effect, tfxStr &field, bool value) {
		if (field == "loop")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].animation_flags |= value ? tfxAnimationFlags_loop : 0;
		if (field == "seamless")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].animation_flags |= value ? tfxAnimationFlags_seamless : 0;
		if (field == "export_with_transparency")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].animation_flags |= value ? tfxAnimationFlags_export_with_transparency : 0;
		if (field == "camera_isometric")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings.camera_isometric = false;
		if (field == "camera_hide_floor")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings.camera_hide_floor = value;
		if (field == "orthographic_camera_isometric")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_isometric = true;
		if (field == "orthographic_camera_hide_floor")
			effect.library->sprite_sheet_settings[GetEffectInfo(&effect)->sprite_sheet_settings_index].camera_settings_orthographic.camera_hide_floor = value;
		if (field == "preview_attach_effect_to_camera")
			effect.library->preview_camera_settings[GetEffectInfo(&effect)->preview_camera_settings].attach_effect_to_camera = value;
		if (field == "preview_camera_hide_floor")
			effect.library->preview_camera_settings[GetEffectInfo(&effect)->preview_camera_settings].camera_settings.camera_hide_floor = value;
		if (field == "preview_camera_isometric")
			effect.library->preview_camera_settings[GetEffectInfo(&effect)->preview_camera_settings].camera_settings.camera_isometric = value;
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
		if (field == "wrap_single_sprite")
			if (value) effect.property_flags |= tfxEmitterPropertyFlags_wrap_single_sprite; else effect.property_flags &= ~tfxEmitterPropertyFlags_wrap_single_sprite;
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
		if (field == "include_in_sprite_data_export")
			if (value) effect.effect_flags |= tfxEffectPropertyFlags_include_in_sprite_data_export; else effect.effect_flags &= ~tfxEffectPropertyFlags_include_in_sprite_data_export;
	}

	void StreamProperties(tfxEmitterPropertiesSoA &property, tfxU32 index, tfxEmitterPropertyFlags &flags, tfxStr &file) {

		file.AddLine("image_hash=%llu", property.image_hash[index]);
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
		file.AddLine("wrap_single_sprite=%i", (flags & tfxEmitterPropertyFlags_wrap_single_sprite));
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
		file.AddLine("include_in_sprite_data_export=%i", effect.effect_flags & tfxEffectPropertyFlags_include_in_sprite_data_export);
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

	bool IsOvertimeGraph(tfxGraph *graph) {
		return graph->type >= tfxOvertime_velocity && graph->type <= tfxOvertime_noise_resolution && graph->type != tfxOvertime_velocity_adjuster;
	}

	bool IsColorGraph(tfxGraph *graph) {
		return graph->type >= tfxOvertime_red && graph->type <= tfxOvertime_blue;
	}

	bool IsGlobalGraph(tfxGraph *graph) {
		return graph->type >= tfxGlobal_life && graph->type <= tfxGlobal_splatter;
	}

	bool IsAngleGraph(tfxGraph *graph) {
		return (graph->type == tfxTransform_roll || graph->type == tfxTransform_pitch || graph->type == tfxTransform_yaw || graph->type == tfxProperty_emission_pitch || graph->type == tfxProperty_emission_yaw
			|| graph->type == tfxProperty_emission_range || graph->type == tfxProperty_arc_offset || graph->type == tfxProperty_arc_size || graph->type == tfxBase_spin || graph->type == tfxVariation_spin || graph->type == tfxOvertime_direction);
	}

	bool IsTranslationGraph(tfxGraph *graph) {
		return graph->type == tfxTransform_translate_x || graph->type == tfxTransform_translate_y || graph->type == tfxTransform_translate_z;
	}

	void MultiplyAllGraphValues(tfxGraph *graph, float scalar) {
		graph->nodes.ResetIteratorIndex();
		do {
			for (auto &node : graph->nodes) {
				node.value *= scalar;
				node.left.y *= scalar;
				node.right.y *= scalar;
			}
		} while (!graph->nodes.EndOfBuckets());
	}

	void CopyGraphNoLookups(tfxGraph *src_graph, tfxGraph *dst_graph) {
		dst_graph->min = src_graph->min;
		dst_graph->max = src_graph->max;
		dst_graph->graph_preset = src_graph->graph_preset;
		dst_graph->type = src_graph->type;
		dst_graph->effector = src_graph->effector;
		dst_graph->nodes = src_graph->nodes;
		dst_graph->index = src_graph->index;
		dst_graph->lookup.life = src_graph->lookup.life;
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

		if (SortGraph(&graph)) {
			ReIndexGraph(&graph);
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

		if (SortGraph(&graph)) {
			ReIndexGraph(&graph);
			return true;
		}

		return false;

	}

	bool IsNodeCurve(tfxAttributeNode *node) {
		return node->flags & tfxAttributeNodeFlags_is_curve;
	}

	bool NodeCurvesAreInitialised(tfxAttributeNode *node) {
		return node->flags & tfxAttributeNodeFlags_curves_initialised;
	}

	bool SetNodeCurveInitialised(tfxAttributeNode *node) {
		return node->flags |= tfxAttributeNodeFlags_curves_initialised;
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

		if (SortGraph(&graph)) {
			ReIndexGraph(&graph);
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
			if (SortGraph(&graph)) {
				ReIndexGraph(&graph);
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

		tfxAttributeNode *next = GetGraphNextNode(&graph, node);
		if (next) {
			if (p.x > next->frame) p.x = next->frame;
		}

		tfxAttributeNode *prev = GetGraphPrevNode(&graph, node);
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

	tfxAttributeNode* AddGraphNode(tfxGraph *graph, float _frame, float _value, tfxAttributeNodeFlags flags, float _c0x, float _c0y, float _c1x, float _c1y) {
		tfxAttributeNode node;

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
		ClampNode(*graph, node);
		graph->nodes.push_back(node);
		SortGraph(graph);

		ReIndexGraph(graph);
		return &graph->nodes.back();
	}

	void AddGraphNode(tfxGraph *graph, tfxAttributeNode &node) {
		graph->nodes.ResetIteratorIndex();
		do {
			for (auto &n : graph->nodes) {
				if (n.frame == node.frame && n.value == node.value)
					return;
			}
		} while (!graph->nodes.EndOfBuckets());
		graph->nodes.push_back(node);
		SortGraph(graph);
		ReIndexGraph(graph);
	}

	tfxAttributeNode* AddGraphCoordNode(tfxGraph *graph, float _frame, float _value) {
		tfxAttributeNode node;

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
		ClampNode(*graph, node);
		tfxAttributeNode &n = graph->nodes.push_back(node);
		if (SortGraph(graph)) {
			ReIndexGraph(graph);
			return graph->nodes.find(n);
		}

		ReIndexGraph(graph);
		return &n;
	}

	tfxAttributeNode* InsertGraphCoordNode(tfxGraph *graph, float _frame, float _value) {
		tfxAttributeNode node;

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
		ClampNode(*graph, node);

		if (graph->nodes.size() > 1) {
			tfxAttributeNode *last_node = nullptr;
			for (auto *n = graph->nodes.begin() + 1; n != graph->nodes.end(); ++n) {
				if (node.frame < n->frame)
					last_node = n;
				else
					break;
			}

			if (last_node) {
				tfxAttributeNode *r_value = graph->nodes.insert(last_node, node);
				ReIndexGraph(graph);
				return r_value;
			}
		}

		tfxAttributeNode *r_value = &graph->nodes.push_back(node);
		ReIndexGraph(graph);
		return r_value;
	}

	tfxAttributeNode* InsertGraphNode(tfxGraph *graph, float _frame, float _value) {
		tfxAttributeNode node;

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
		ClampNode(*graph, node);

		if (graph->nodes.size() > 1) {
			tfxAttributeNode *last_node = nullptr;
			graph->nodes.ResetIteratorIndex();
			do {
				for (auto &n : graph->nodes) {
					if (node.frame < n.frame)
						last_node = &n;
					else
						break;
				}
			} while (!graph->nodes.EndOfBuckets());

			if (last_node) {
				tfxAttributeNode *r_value = graph->nodes.insert(last_node, node);
				ReIndexGraph(graph);
				return r_value;
			}
		}

		tfxAttributeNode *r_value = &graph->nodes.push_back(node);
		ReIndexGraph(graph);
		return r_value;
	}

	void SetGraphNode(tfxGraph *graph, tfxU32 i, float _frame, float _value, tfxAttributeNodeFlags flags, float _c0x, float _c0y, float _c1x, float _c1y) {
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

	tfxBucketArray<tfxAttributeNode>& GraphNodes(tfxGraph *graph) {
		return graph->nodes;
	}

	float GetGraphValue(tfxGraph *graph, float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		tfxAttributeNode *lastec = nullptr;
		graph->nodes.ResetIteratorIndex();
		do {
			for (auto &a : graph->nodes) {
				if (age < a.frame) {
					p = (age - lastf) / (a.frame - lastf);
					float bezier_value = GetBezierValue(lastec, a, p, graph->min.y, graph->max.y);
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
		} while (!graph->nodes.EndOfBuckets());
		return lastv;

	}

	tfxAttributeNode *GetGraphNextNode(tfxGraph *graph, tfxAttributeNode &node) {
		if (node.index < graph->nodes.size() - 1) {
			return &graph->nodes[node.index + 1];
		}

		return nullptr;
	}

	tfxAttributeNode *GetGraphPrevNode(tfxGraph *graph, tfxAttributeNode &node) {
		if (node.index > 0) {
			return &graph->nodes[node.index - 1];
		}

		return nullptr;
	}

	tfxAttributeNode *GetGraphLastNode(tfxGraph *graph) {
		return &graph->nodes.back();
	}

	float GetGraphRandomValue(tfxGraph *graph, float age, tfxRandom &random) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		tfxAttributeNode *lastec = nullptr;
		graph->nodes.ResetIteratorIndex();
		do {
			for (auto &a : graph->nodes) {
				if (age < a.frame) {
					p = (age - lastf) / (a.frame - lastf);
					float bezier_value = GetBezierValue(lastec, a, p, graph->min.y, graph->max.y);
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
		} while (!graph->nodes.EndOfBuckets());
		return random.Range(lastv);

	}

	float GetGraphValue(tfxGraph *graph, float age, float life) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		tfxAttributeNode *lastec = nullptr;
		graph->nodes.ResetIteratorIndex();
		do {
			for (auto &a : graph->nodes) {
				float frame = a.frame * life;
				if (age < frame) {
					p = (age - lastf) / (frame - lastf);
					float bezier_value = GetBezierValue(lastec, a, p, graph->min.y, graph->max.y);
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
		} while (!graph->nodes.EndOfBuckets());
		return lastv;
	}

	float GetGraphFirstValue(tfxGraph *graph) {
		if (graph->nodes.size())
			return graph->nodes.front().value;
		return 0.f;
	}

	float* LinkGraphFirstValue(tfxGraph *graph) {
		if (graph->nodes.size())
			return &graph->nodes.front().value;
		return nullptr;
	}

	float GetGraphLastValue(tfxGraph *graph) {
		if (graph->nodes.size())
			return graph->nodes.back().value;

		return 0.f;
	}

	float GetGraphMaxValue(tfxGraph *graph) {
		if (graph->nodes.size()) {
			float value = tfxMIN_FLOAT;
			graph->nodes.ResetIteratorIndex();
			do {
				for (auto &n : graph->nodes) {
					if (value < n.value)
						value = n.value;
				}
			} while (!graph->nodes.EndOfBuckets());
			return value;
		}
		return 0.f;
	}

	float GetGraphMinValue(tfxGraph *graph) {
		if (graph->nodes.size()) {
			float value = tfxMAX_FLOAT;
			graph->nodes.ResetIteratorIndex();
			do {
				for (auto &n : graph->nodes) {
					if (value > n.value)
						value = n.value;
				}
			} while (!graph->nodes.EndOfBuckets());
			return value;
		}
		return 0.f;
	}

	float GetGraphLastFrame(tfxGraph *graph, float update_frequency) {
		float frame_length = 1000.f / update_frequency;
		if (graph->nodes.size()) {
			return graph->nodes.size() > 1 && graph->nodes.back().frame == 0 ? frame_length : graph->nodes.back().frame;
		}

		return 0.f;
	}

	tfxAttributeNode* FindGraphNode(tfxGraph *graph, const tfxAttributeNode &n) {
		return graph->nodes.find(n);
	}

	void ValidateGraphCurves(tfxGraph *graph) {
		tfxU32 index = 0;
		tfxU32 last_index = graph->nodes.size() - 1;
		graph->nodes.ResetIteratorIndex();
		do {
			for (auto &n : graph->nodes) {
				if (n.flags & tfxAttributeNodeFlags_is_curve) {
					if (index < last_index) {
						if (graph->nodes[index + 1].frame < n.right.x)
							n.right.x = graph->nodes[index + 1].frame;
					}
					if (index > 0) {
						if (graph->nodes[index - 1].frame > n.left.x)
							n.left.x = graph->nodes[index - 1].frame;
					}
					if (n.left.x > n.frame)
						n.left.x = n.frame;
					if (n.right.x < n.frame)
						n.right.x = n.frame;
				}
				index++;
			}
		} while (!graph->nodes.EndOfBuckets());
	}

	void DeleteGraphNode(tfxGraph *graph, const tfxAttributeNode &n) {
		graph->nodes.erase(&n);
	}

	void DeleteGraphNodeAtFrame(tfxGraph *graph, float frame) {
		graph->nodes.ResetIteratorIndex();
		do {
			for (auto &n : graph->nodes) {
				if (n.frame == frame) {
					graph->nodes.erase(&n);
					return;
				}
			}
		} while (!graph->nodes.EndOfBuckets());

	}

	void ResetGraph(tfxGraph *graph, float v, tfxGraphPreset preset, bool add_node) {
		graph->nodes.clear();
		graph->nodes.TrimBuckets();
		if (add_node && preset == tfxWeightOvertimePreset) {
			AddGraphNode(graph, 0.f, 0.f, 0);
			tfxAttributeNode *node = AddGraphNode(graph, 1.f, 1.f, tfxAttributeNodeFlags_is_curve, 0.f, 1.f, 1.f, 1.f);
			SetNodeCurveInitialised(node);
		}
		else if (add_node) {
			AddGraphNode(graph, 0.f, v);
		}
		switch (preset) {
		case tfxGraphPreset::tfxGlobalPercentPreset:
			//We have a epsilon to prevent divide by 0 here
			graph->min = { 0.f, 0.0001f }; graph->max = { tfxMAX_FRAME, 20.f };
			break;
		case tfxGraphPreset::tfxGlobalOpacityPreset:
			graph->min = { 0.f, 0.f }; graph->max = { tfxMAX_FRAME, 1.f };
			break;
		case tfxGraphPreset::tfxGlobalPercentPresetSigned:
			graph->min = { 0.f, -20.f }; graph->max = { tfxMAX_FRAME, 20.f };
			break;
		case tfxGraphPreset::tfxAnglePreset:
			graph->min = { 0.f, -1080.f }; graph->max = { tfxMAX_FRAME, 1080.f };
			break;
		case tfxGraphPreset::tfxArcPreset:
			graph->min = { 0.f, 0.f }; graph->max = { tfxMAX_FRAME, 360.f };
			break;
		case tfxGraphPreset::tfxEmissionRangePreset:
			graph->min = { 0.f, 0.f }; graph->max = { tfxMAX_FRAME, 360.f };
			break;
		case tfxGraphPreset::tfxDimensionsPreset:
			graph->min = { 0.f, 0.f }; graph->max = { tfxMAX_FRAME, 4000.f };
			break;
		case tfxGraphPreset::tfxTranslationPreset:
			graph->min = { 0.f, -4000.f }; graph->max = { tfxMAX_FRAME, 4000.f };
			break;
		case tfxGraphPreset::tfxLifePreset:
			//We have a epsilon to prevent divide by 0 here. The divide by zero occurrs in control functions (ControlParticleImageFrame3d etc.) when the current % life of the particle is calculated
			graph->min = { 0.f, 0.0001f }; graph->max = { tfxMAX_FRAME, 100000.f };
			break;
		case tfxGraphPreset::tfxAmountPreset:
			graph->min = { 0.f, 0.f }; graph->max = { tfxMAX_FRAME, 5000.f };
			break;
		case tfxGraphPreset::tfxVelocityPreset:
			graph->min = { 0.f, 0.f }; graph->max = { tfxMAX_FRAME, 10000.f };
			break;
		case tfxGraphPreset::tfxWeightPreset:
			graph->min = { 0.f, -10000.f }; graph->max = { tfxMAX_FRAME, 10000.f };
			break;
		case tfxGraphPreset::tfxWeightVariationPreset:
			graph->min = { 0.f, 0.f }; graph->max = { tfxMAX_FRAME, 20000.f };
			break;
		case tfxGraphPreset::tfxNoiseOffsetVariationPreset:
			graph->min = { 0.f, 0.f }; graph->max = { tfxMAX_FRAME, 1000.f };
			break;
		case tfxGraphPreset::tfxNoiseResolutionPreset:
			graph->min = { 0.f, 0.f }; graph->max = { tfxMAX_FRAME, 10000.f };
			break;
		case tfxGraphPreset::tfxSpinPreset:
			graph->min = { 0.f, -2000.f }; graph->max = { tfxMAX_FRAME, 2000.f };
			break;
		case tfxGraphPreset::tfxSpinVariationPreset:
			graph->min = { 0.f, 0.f }; graph->max = { tfxMAX_FRAME, 2000.f };
			break;
		case tfxGraphPreset::tfxDirectionVariationPreset:
			graph->min = { 0.f, 0.f }; graph->max = { tfxMAX_FRAME, 22.5f };
			break;
		case tfxGraphPreset::tfxWeightOvertimePreset:
			graph->min = { 0.f, -20.f }; graph->max = { 1.f, 20.f };
			break;
		case tfxGraphPreset::tfxDirectionOvertimePreset:
			graph->min = { 0.f, 0.f }; graph->max = { 1.f, 4320.f };
			break;
		case tfxGraphPreset::tfxSpinOvertimePreset:
			graph->min = { 0.f, 0.f }; graph->max = { 1.f, 20.f };
			break;
		case tfxGraphPreset::tfxVelocityOvertimePreset:
			graph->min = { 0.f, -20.f }; graph->max = { 1.f, 20.f };
			break;
		case tfxGraphPreset::tfxPercentOvertime:
			graph->min = { 0.f, 0.f }; graph->max = { 1.f, 20.f };
			break;
		case tfxGraphPreset::tfxFrameratePreset:
			graph->min = { 0.f, 0.f }; graph->max = { 1.f, 200.f };
			break;
		case tfxGraphPreset::tfxVelocityTurbulancePreset:
			graph->min = { 0.f, 0.f }; graph->max = { 1.f, 2000.f };
			break;
		case tfxGraphPreset::tfxOpacityOvertimePreset:
			graph->min = { 0.f, 0.f }; graph->max = { 1.f, 1.f };
			break;
		case tfxGraphPreset::tfxColorPreset:
			graph->min = { 0.f, 0.f }; graph->max = { 1.f, 255.f };
			break;
		case tfxGraphPreset::tfxIntensityOvertimePreset:
			graph->min = { 0.f, 0.f }; graph->max = { 1.f, 5.f };
			break;
		}

		graph->graph_preset = preset;
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

	void DragGraphValues(tfxGraphPreset preset, float &frame, float &value) {
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

	void ClearGraphToOne(tfxGraph *graph, float value) {
		graph->nodes.clear();
		AddGraphNode(graph, 0.f, value);
	}

	void ClearGraph(tfxGraph *graph) {
		graph->nodes.clear();
	}

	void FreeGraph(tfxGraph *graph) {
		//Explicitly free the nodes
		graph->nodes.free_all();
		graph->lookup.values.free();
	}

	void CopyGraph(tfxGraph *from, tfxGraph *to, bool compile) {
		ClearGraph(to);
		do {
			for (auto &n : from->nodes) {
				to->nodes.push_back(n);
			}
		} while (!from->nodes.EndOfBuckets());
		if (compile) {
			if (IsColorGraph(from))
				CompileColorOvertime(*to);
			else if (IsOvertimeGraph(from))
				CompileGraphOvertime(*to);
			else
				CompileGraph(*to);
		}
	}

	bool SortGraph(tfxGraph *graph) {
		bool needed_sorting = false;
		for (tfxU32 i = 1; i < graph->nodes.current_size; ++i) {
			tfxAttributeNode key = graph->nodes[i];
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

	void ReIndexGraph(tfxGraph *graph) {
		tfxU32 i = 0;
		graph->nodes.ResetIteratorIndex();
		do {
			for (auto &a : graph->nodes) {
				a.index = i++;
			}
		} while (!graph->nodes.EndOfBuckets());
	}

	tfxVec2 GetGraphInitialZoom(tfxGraph *graph) {
		switch (graph->graph_preset) {
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

	tfxVec2 GetGraphInitialZoom3d(tfxGraph *graph) {
		switch (graph->graph_preset) {
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
		float last_frame = GetGraphLastFrame(&graph, 60.f);
		graph.lookup.last_frame = tfxU32(last_frame / tfxLOOKUP_FREQUENCY);
		if (graph.lookup.last_frame) {
			graph.lookup.values.resize(graph.lookup.last_frame + 1);
			for (tfxU32 f = 0; f != graph.lookup.last_frame + 1; ++f) {
				graph.lookup.values[f] = GetGraphValue(&graph, (float)f * tfxLOOKUP_FREQUENCY);
			}
			graph.lookup.values[graph.lookup.last_frame] = GetGraphLastValue(&graph);
		}
		else {
			graph.lookup.values.resize(1);
			graph.lookup.values[0] = GetGraphFirstValue(&graph);
		}
	}

	void CompileGraphOvertime(tfxGraph &graph) {
		if (graph.nodes.size() > 1) {
			graph.lookup.last_frame = tfxU32(graph.lookup.life / tfxLOOKUP_FREQUENCY_OVERTIME);
			graph.lookup.values.resize(graph.lookup.last_frame + 1);
			for (tfxU32 f = 0; f != graph.lookup.last_frame + 1; ++f) {
				graph.lookup.values[f] = GetGraphValue(&graph, (float)f * tfxLOOKUP_FREQUENCY_OVERTIME, graph.lookup.life);
			}
			graph.lookup.values[graph.lookup.last_frame] = GetGraphLastValue(&graph);
		}
		else {
			graph.lookup.last_frame = 0;
			graph.lookup.values.resize(1);
			graph.lookup.values[0] = GetGraphFirstValue(&graph);
		}
	}

	void CompileColorOvertime(tfxGraph &graph, float gamma) {
		if (graph.nodes.size() > 1) {
			graph.lookup.last_frame = tfxU32(graph.lookup.life / tfxLOOKUP_FREQUENCY_OVERTIME);
			graph.lookup.values.resize(graph.lookup.last_frame + 1);
			for (tfxU32 f = 0; f != graph.lookup.last_frame + 1; ++f) {
				graph.lookup.values[f] = GammaCorrect(GetGraphValue(&graph, (float)f * tfxLOOKUP_FREQUENCY_OVERTIME, graph.lookup.life), gamma);
			}
			graph.lookup.values[graph.lookup.last_frame] = GammaCorrect(GetGraphLastValue(&graph), gamma);
		}
		else {
			graph.lookup.last_frame = 0;
			graph.lookup.values.resize(1);
			graph.lookup.values[0] = GammaCorrect(GetGraphFirstValue(&graph), gamma);
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
		return GetGraphRandomValue(&graph, frame, random);
	}

	float GetMaxLife(tfxEffectEmitter &e) {
		tfxGraph &life = *GetEffectGraphByType(&e, tfxBase_life);
		tfxGraph &life_variation = *GetEffectGraphByType(&e, tfxVariation_life);
		float templife = 0;
		float max_life = 0;
		float life_last_frame = GetGraphLastFrame(&life, 60.f);
		float life_variation_last_frame = GetGraphLastFrame(&life_variation, 60.f);
		float global_adjust = 1.f;
		if (life_last_frame + life_variation_last_frame > 0) {
			for (float f = 0; f < std::fmaxf(life_last_frame, life_variation_last_frame); ++f) {
				if (e.parent)
					global_adjust = GetGraphValue(GetEffectGraphByType(e.parent, tfxGlobal_life), f);
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

	bool IsOvertimeGraph(tfxGraphType type) {
		return type >= tfxOvertime_velocity && type != tfxOvertime_noise_resolution && type <= tfxOvertime_noise_resolution;
	}

	bool IsColorGraph(tfxGraphType type) {
		return type >= tfxOvertime_red && type <= tfxOvertime_blue;
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

	bool HasNodeAtFrame(tfxGraph &graph, float frame) {
		graph.nodes.ResetIteratorIndex();
		do {
			for (auto &node : graph.nodes) {
				if (node.frame == frame) return true;
			}
		} while (!graph.nodes.EndOfBuckets());
		return false;
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
		FILE *file = tfx__open_file(path, "w");

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
		FILE* fp = tfx__open_file(path, "r");
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
	int GetShapeCountInLibrary(const char *filename) {
		int context = 0;
		int error = 0;

		tfxPackage package;
		error = LoadPackage(filename, package);

		tfxEntryInfo *data = GetPackageFile(&package, "data.txt");

		if (!data)
			error = -5;


		if (error < 0) {
			FreePackage(&package);
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

		FreePackage(&package);
		return shape_count;
	}

	int GetEffectLibraryStats(const char *filename, tfxEffectLibraryStats &stats) {
		int context = 0;
		int error = 0;

		tfxPackage package;
		error = LoadPackage(filename, package);

		tfxEntryInfo *data = GetPackageFile(&package, "data.txt");

		if (!data)
			error = -5;


		if (error < 0) {
			FreePackage(&package);
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
			for (auto &sub : GetEffectInfo(&current)->sub_effectors) {
				stack.push_back(sub);
			}
		}
		stats.total_shapes = lib.particle_shapes.data.size();
		stats.required_graph_node_memory = lib.graph_node_allocator.TotalMemoryInUse();
		stats.required_graph_lookup_memory = lib.graph_lookup_allocator.TotalMemoryInUse();

		return stats;
	}

	tfxAPI tfxErrorFlags LoadSpriteData(const char *filename, tfxAnimationManager &animation_manager, void(*shape_loader)(const char *filename, tfxImageData &image_data, void *raw_image_data, int image_size, void *user_data), void *user_data) {
		//assert(shape_loader);			//Must have a shape_loader function to load your shapes with. This will be a custom user function suited for whichever renderer you're using
		if (!tfxDataTypes.initialised)
			tfxDataTypes.Init();

		tfxPackage package;
		tfxErrorFlags error = LoadPackage(filename, package);
		if (error != 0) {
			FreePackage(&package);
			return error;
		}

		tfxEntryInfo *data = GetPackageFile(&package, "data.txt");

		if (!data) {
			error |= tfxErrorCode_data_could_not_be_loaded;
		}

		tfxEntryInfo *sprite_data = GetPackageFile(&package, "sprite_data");

		if (!sprite_data) {
			error |= tfxErrorCode_data_could_not_be_loaded;
		}

		if (error != 0) {
			FreePackage(&package);
			return error;
		}

		if (package.header.user_data2 == 1 && !(animation_manager.flags & tfxAnimationManagerFlags_is_3d)) {
			return tfxErrorCode_sprite_data_is_3d_but_animation_manager_is_2d;
		}

		if (package.header.user_data2 == 0 && animation_manager.flags & tfxAnimationManagerFlags_is_3d) {
			return tfxErrorCode_sprite_data_is_2d_but_animation_manager_is_3d;
		}

		if (animation_manager.flags & tfxAnimationManagerFlags_is_3d) {
			animation_manager.sprite_data_3d.resize((tfxU32)(sprite_data->file_size / package.header.user_data1));
			memcpy(animation_manager.sprite_data_3d.data, sprite_data->data.data, sprite_data->file_size);
		}
		else {
			animation_manager.sprite_data_2d.resize((tfxU32)(sprite_data->file_size / package.header.user_data1));
			memcpy(animation_manager.sprite_data_2d.data, sprite_data->data.data, sprite_data->file_size);
		}

		tmpStack(tfxSpriteDataMetrics, metrics_stack);
		tmpStack(tfxFrameMeta, frame_meta_stack);
		tmpStack(tfxAnimationEmitterProperties, emitter_properties_stack);
		tmpStack(tfxStr256, pair);
		tmpStack(tfxStr256, multi);

		tfxKey first_shape_hash = 0;
		int context = 0;

		while (!data->data.EoF()) {
			tfxStr512 line = data->data.ReadLine();
			bool context_set = false;

			if (StringIsUInt(line.c_str())) {
				context = atoi(line.c_str());
				if (context == tfxEndOfFile) {
					break;
				}

				context_set = true;
				if (context == tfxStartEffectAnimationInfo) {
					tfxSpriteDataMetrics metrics;
					metrics_stack.push_back(metrics);
				}
				else if (context == tfxStartFrameMeta) {
					tfxFrameMeta frame_meta;
					frame_meta_stack.push_back(frame_meta);
				}
				else if (context == tfxStartEmitter) {
					tfxAnimationEmitterProperties emitter_properties;
					emitter_properties_stack.push_back(emitter_properties);
				}
			}

			if (context_set == false) {
				pair.clear();
				SplitStringStack(line.c_str(), pair);
				if (pair.size() != 2) {
					pair.clear();
					SplitStringStack(line.c_str(), pair, ',');
					if (pair.size() < 2) {
						error |= tfxErrorCode_some_data_not_loaded;
						continue;
					}
				}

				if (context == tfxStartEffectAnimationInfo) {
					if (tfxDataTypes.names_and_types.ValidName(pair[0])) {
						switch (tfxDataTypes.names_and_types.At(pair[0])) {
						case tfxUint:
							AssignSpriteDataMetricsProperty(metrics_stack.back(), pair[0], (tfxU32)atoi(pair[1].c_str()), package.header.file_version);
							break;
						case tfxFloat:
							AssignSpriteDataMetricsProperty(metrics_stack.back(), pair[0], (float)atof(pair[1].c_str()), package.header.file_version);
							break;
						case tfxString:
							AssignSpriteDataMetricsProperty(metrics_stack.back(), pair[0], pair[1], package.header.file_version);
							break;
						}
					}
					else {
						error |= tfxErrorCode_some_data_not_loaded;
					}
				}
				else if (context == tfxStartFrameMeta) {
					if (tfxDataTypes.names_and_types.ValidName(pair[0])) {
						switch (tfxDataTypes.names_and_types.At(pair[0])) {
						case tfxUint:
							AssignFrameMetaProperty(frame_meta_stack.back(), pair[0], (tfxU32)atoi(pair[1].c_str()), package.header.file_version);
							break;
						case tfxFloat3:
							multi.clear();
							SplitStringStack(pair[1], multi, ',');
							AssignFrameMetaProperty(frame_meta_stack.back(), pair[0], StrToVec3(multi), package.header.file_version);
							break;
						}
					}
				}
				else if (context == tfxStartEmitter) {
					if (tfxDataTypes.names_and_types.ValidName(pair[0])) {
						switch (tfxDataTypes.names_and_types.At(pair[0])) {
						case tfxUint:
							AssignAnimationEmitterProperty(emitter_properties_stack.back(), pair[0], (tfxU32)atoi(pair[1].c_str()), package.header.file_version);
							break;
						case tfxFloat:
							AssignAnimationEmitterProperty(emitter_properties_stack.back(), pair[0], (float)atof(pair[1].c_str()), package.header.file_version);
							break;
						case tfxFloat2:
							multi.clear();
							SplitStringStack(pair[1], multi, ',');
							AssignAnimationEmitterProperty(emitter_properties_stack.back(), pair[0], StrToVec2(multi), package.header.file_version);
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
						tfxShapeData s;
						strcpy_s(s.name, pair[0].c_str());
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

						tfxEntryInfo *shape_entry = GetPackageFile(&package, s.name);
						if (shape_entry) {
							tfxImageData image_data;
							image_data.shape_index = s.shape_index;
							image_data.animation_frames = (float)s.frame_count;
							image_data.image_size = tfxVec2((float)s.width, (float)s.height);
							image_data.name = s.name;
							image_data.import_filter = s.import_filter;
							image_data.image_hash = tfxXXHash64::hash(shape_entry->data.data, shape_entry->file_size, 0);
							if (s.image_hash == 0) {
								s.image_hash = image_data.image_hash;
							}
							assert(s.image_hash == image_data.image_hash);

							shape_loader(s.name, image_data, shape_entry->data.data, (tfxU32)shape_entry->file_size, user_data);

							if (!image_data.ptr) {
								//uid = -6;
							}
							else {
								animation_manager.particle_shapes.Insert(image_data.image_hash, image_data);
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
				assert(metrics_stack.current_size);	
				assert(frame_meta_stack.current_size);
				metrics_stack.back().frame_meta.push_back(frame_meta_stack.pop_back());
			}
			else if (context == tfxEndEffectAnimationInfo) {
				assert(metrics_stack.current_size);
				animation_manager.effect_animation_info.Insert(metrics_stack.back().name, metrics_stack.back());
				metrics_stack.pop();
			}
			else if (context == tfxEndEmitter) {
				assert(emitter_properties_stack.current_size);
				animation_manager.emitter_properties.push_back(emitter_properties_stack.pop_back());
			}

		}

		FreePackage(&package);

		if (!shape_loader) {
			error |= tfxErrorCode_library_loaded_without_shape_loader;
		}

		return error;
	}

	tfxErrorFlags LoadEffectLibraryPackage(tfxPackage &package, tfxLibrary &lib, void(*shape_loader)(const char *filename, tfxImageData &image_data, void *raw_image_data, int image_size, void *user_data), void *user_data, bool read_only) {

		assert(shape_loader);			//Must have a shape_loader function to load your shapes with. This will be a custom user function suited for whichever renderer you're using
		if (!tfxDataTypes.initialised)
			tfxDataTypes.Init();

		ClearLibrary(&lib);
		if (tfxIcospherePoints[0].current_size == 0) {
			MakeIcospheres();
		}

		tfxEntryInfo *data = GetPackageFile(&package, "data.txt");
		tfxEntryInfo *stats_struct = GetPackageFile(&package, "stats.struct");
		tfxErrorFlags error = 0;

		int context = 0;
		int uid = 0;
		tfxU32 current_global_graph = 0;

		InitLibraryEmitterProperties(&lib);

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
			FreePackage(&package);
			return error;
		}

		tfxKey first_shape_hash = 0;

		//You must call InitialiseTimelineFX() before doing anything!	
		assert(tfxSTACK_ALLOCATOR.arena_size > 0);
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
					effect.info_index = AddLibraryEffectEmitterInfo(&lib);
					GetEffectInfo(&effect)->uid = uid++;
					effect_stack.push_back(effect);
				}
				else if (context == tfxStartStage) {
					tfxEffectEmitter effect;
					effect.library = &lib;
					effect.type = tfxEffectEmitterType::tfxStage;
					effect.info_index = AddLibraryEffectEmitterInfo(&lib);
					AddLibraryPreviewCameraSettings(&lib, effect);
					effect.transform_attributes = AddLibraryKeyframes(&lib);
					GetEffectInfo(&effect)->uid = uid++;
					effect_stack.push_back(effect);
				}
				else if (context == tfxStartEffect) {
					tfxEffectEmitter effect;
					effect.library = &lib;
					effect.info_index = AddLibraryEffectEmitterInfo(&lib);
					effect.property_index = AddLibraryEmitterProperties(&lib);
					effect.transform_attributes = AddLibraryKeyframes(&lib);
					if (effect_stack.size() <= 1) { //Only root effects get the global graphs
						AddLibraryEffectGraphs(&lib, effect);
						ResetEffectGraphs(&effect, false, false);
						current_global_graph = effect.global;
					}
					AddLibraryTransformGraphs(&lib, effect);
					ResetTransformGraphs(&effect, false, false);
					effect.type = tfxEffectEmitterType::tfxEffectType;
					AddLibrarySpriteSheetSettings(&lib, effect);
					AddLibrarySpriteDataSettings(&lib, effect);
					AddLibraryPreviewCameraSettings(&lib, effect);
					GetEffectInfo(&effect)->uid = uid++;
					effect_stack.push_back(effect);

				}
				else if (context == tfxStartEmitter) {
					tfxEffectEmitter emitter;
					emitter.library = &lib;
					emitter.info_index = AddLibraryEffectEmitterInfo(&lib);
					emitter.property_index = AddLibraryEmitterProperties(&lib);
					emitter.transform_attributes = AddLibraryKeyframes(&lib);
					AddLibraryEmitterGraphs(&lib, emitter);
					AddLibraryTransformGraphs(&lib, emitter);
					emitter.type = tfxEffectEmitterType::tfxEmitterType;
					ResetEmitterGraphs(&emitter, false, false);
					ResetTransformGraphs(&emitter, false, false);
					GetEffectInfo(&emitter)->uid = uid++;
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
						case tfxUInt64:
							AssignEffectorProperty(effect_stack.back(), pair[0], (tfxU64)strtoull(pair[1].c_str(), NULL, 10), package.header.file_version);
							break;
						case tfxUint:
							AssignEffectorProperty(effect_stack.back(), pair[0], (tfxU32)atoi(pair[1].c_str()), package.header.file_version);
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
						s.import_filter = atoi(pair[5].c_str());
						if (pair.current_size > 6) {
							s.image_hash = strtoull(pair[6].c_str(), NULL, 10);
						}
						if (s.import_filter < 0 || s.import_filter>1) {
							s.import_filter = 0;
						}

						tfxEntryInfo *shape_entry = GetPackageFile(&package, s.name);
						if (shape_entry) {
							tfxImageData image_data;
							image_data.shape_index = s.shape_index;
							image_data.animation_frames = (float)s.frame_count;
							image_data.image_size = tfxVec2((float)s.width, (float)s.height);
							image_data.name = s.name;
							image_data.import_filter = s.import_filter;
							image_data.image_hash = tfxXXHash64::hash(shape_entry->data.data, shape_entry->file_size, 0);
							if (s.image_hash == 0) {
								s.image_hash = image_data.image_hash;
							}
							assert(s.image_hash == image_data.image_hash);

							shape_loader(s.name, image_data, shape_entry->data.data, (tfxU32)shape_entry->file_size, user_data);

							if (!image_data.ptr) {
								//uid = -6;
							}
							else {
								lib.particle_shapes.Insert(image_data.image_hash, image_data);
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
#ifdef tfxTRACK_MEMORY
				tfxvec<tfxEffectEmitter> &sub_effectors = effect_stack.parent().GetEffectInfo(this)->sub_effectors;
				memcpy(sub_effectors.name, "emitter_sub_effects\0", 20);
#endif
				if (effect_stack.back().property_flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
					lib.emitter_properties.image_handle[effect_stack.back().property_index] = { .5f, .5f };
				}
				effect_stack.back().property_flags |= tfxEmitterPropertyFlags_enabled;
				GetEffectInfo(&effect_stack.parent())->sub_effectors.push_back(effect_stack.back());
				effect_stack.pop();
			}

			if (context == tfxEndEffect) {
				ReIndexEffect(&effect_stack.back());
				if (effect_stack.size() > 1) {
#ifdef tfxTRACK_MEMORY
					tfxvec<tfxEffectEmitter> &sub_effectors = effect_stack.parent().GetEffectInfo(this)->sub_effectors;
					memcpy(sub_effectors.name, "effect_sub_emitters\0", 20);
#endif
					if (effect_stack.parent().type == tfxStage && GetEffectInfo(&effect_stack.parent())->sub_effectors.size() == 0) {
						tfxEffectPropertyFlags tmp = effect_stack.parent().property_flags;
						if (Is3DEffect(&effect_stack.back())) {
							effect_stack.parent().property_flags |= tfxEmitterPropertyFlags_is_3d;
						}
						tmp = effect_stack.parent().property_flags;
					}
					else if (effect_stack.parent().type == tfxEmitterType) {
						effect_stack.back().global = current_global_graph;
					}
					GetEffectInfo(&effect_stack.parent())->sub_effectors.push_back(effect_stack.back());
				}
				else {
					lib.effects.push_back(effect_stack.back());
					InitialiseUninitialisedGraphs(&effect_stack.back());
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

		FreePackage(&package);

		if (uid >= 0) {
			//Effects were loaded so let's compile them
			CompileAllLibraryGraphs(&lib);
			ReIndexLibrary(&lib);
			if (first_shape_hash != 0) {
				UpdateLibraryParticleShapeReferences(&lib, lib.effects, first_shape_hash);
			}
			UpdateLibraryEffectPaths(&lib);
			UpdateLibraryComputeNodes(&lib);
			SetLibraryMinMaxData(&lib);
		}
		lib.uid = uid;

		return error;
	}

	tfxErrorFlags LoadEffectLibraryPackage(const char *filename, tfxLibrary &lib, void(*shape_loader)(const char* filename, tfxImageData &image_data, void *raw_image_data, int image_size, void *user_data), void *user_data, bool read_only) {

		tfxErrorFlags error = 0;

		tfxPackage package;
		error = LoadPackage(filename, package);
		if (error != 0) {
			FreePackage(&package);
			return error;
		}
		error = LoadEffectLibraryPackage(package, lib, shape_loader, user_data, read_only);

		FreePackage(&package);
		return error;
	}

	void SetTemplateUserDataAll(tfxEffectTemplate *t, void *data) {
		tmpStack(tfxEffectEmitter*, stack);
		stack.push_back(&t->effect);
		while (stack.size()) {
			tfxEffectEmitter *current = stack.pop_back();
			current->user_data = data;
			for (auto &sub : GetEffectInfo(current)->sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	void InvalidateNewSpriteCapturedIndex(tfxParticleManager *pm) {
		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			tfxSpriteSoA &sprites = pm->sprites[pm->current_sprite_buffer][layer];
			for (int i = 0; i != pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size; ++i) {
				if ((sprites.captured_index[i] & 0xC0000000) >> 30 == pm->current_sprite_buffer && !(sprites.captured_index[i] & 0x10000000)) {
					sprites.captured_index[i] = tfxINVALID;
				} 
			}
		}
	}

	void ResetSpriteDataLerpOffset(tfxSpriteData &sprite_data) {
		tfxSpriteDataSoA &sprites = sprite_data.real_time_sprites;
		for (int i = 0; i != sprite_data.real_time_sprites_buffer.current_size; ++i) {
			sprites.lerp_offset[i] = 1.f;
		}
	}

	void RecordSpriteData(tfxParticleManager *pm, tfxEffectEmitter *effect, float update_frequency, float camera_position[3]) {
		assert(update_frequency > 0); //Update frequency must be greater then 0. 60 is recommended for best results
		ReconfigureParticleManager(pm, GetRequiredParticleManagerMode(effect), effect->sort_passes, Is3DEffect(effect));
		tfxSpriteDataSettings &anim = effect->library->sprite_data_settings[GetEffectInfo(effect)->sprite_data_settings_index];
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

		anim.recording_frame_rate = tfxMin(tfxMax(30.f, anim.recording_frame_rate), 240.f);

		bool auto_set_length = false;
		if (anim.animation_flags & tfxAnimationFlags_auto_set_length && !(anim.animation_flags & tfxAnimationFlags_loop) && IsFiniteEffect(effect)) {
			frames = 99999;
			auto_set_length = true;
		}

		//First pass to count the number of sprites in each frame
		pm->flags |= tfxEffectManagerFlags_recording_sprites;
		ClearParticleManager(pm, false);
		if (!(pm->flags & tfxEffectManagerFlags_using_uids)) {
			ToggleSpritesWithUID(pm, true);
		}
		SetSeed(pm, anim.seed);
		tfxU32 preview_effect_index = AddEffectToParticleManager(pm, effect, pm->current_ebuff, 0, false, 0.f);
		tfxVec3 pm_camera_position = pm->camera_position;
		pm->camera_position = tfxVec3(camera_position[0], camera_position[1], camera_position[2]);
		SetEffectPosition(pm, preview_effect_index, tfxVec3(0.f, 0.f, 0.f));
		if (is_3d) {
			Transform3d(pm->effects.world_rotations[preview_effect_index],
				pm->effects.local_rotations[preview_effect_index],
				pm->effects.scale[preview_effect_index],
				pm->effects.world_position[preview_effect_index],
				pm->effects.local_position[preview_effect_index],
				pm->effects.translation[preview_effect_index],
				pm->effects.matrix[preview_effect_index],
				pm->effects.world_rotations[preview_effect_index],
				pm->effects.scale[preview_effect_index],
				pm->effects.world_position[preview_effect_index],
				pm->effects.matrix[preview_effect_index]
			);
		}
		else {
			Transform2d(pm->effects.world_rotations[preview_effect_index],
				pm->effects.local_rotations[preview_effect_index],
				pm->effects.scale[preview_effect_index],
				pm->effects.world_position[preview_effect_index],
				pm->effects.local_position[preview_effect_index],
				pm->effects.translation[preview_effect_index],
				pm->effects.matrix[preview_effect_index],
				pm->effects.world_rotations[preview_effect_index],
				pm->effects.scale[preview_effect_index],
				pm->effects.world_position[preview_effect_index],
				pm->effects.matrix[preview_effect_index]);
		}

		tfxU32 total_sprites = 0;
		tfxvec<tfxFrameMeta> tmp_frame_meta;
		tfxU32 sprites_in_layers = 0;
		while (frame < frames && offset < 99999) {
			tfxU32 count_this_frame = 0;
			UpdateParticleManager(pm, frame_length);
			bool particles_processed_last_frame = false;

			if (offset >= start_frame) {
				sprites_in_layers = 0;
				for (tfxEachLayer) {
					if (frame >= tmp_frame_meta.size()) {
						tfxFrameMeta meta;
						memset(&meta, 0, sizeof(tfxFrameMeta));
						tmp_frame_meta.push_back(meta);
					}
					tmp_frame_meta[frame].sprite_count[layer] += pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size;
					tmp_frame_meta[frame].total_sprites += pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size;
					total_sprites += pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size;
					sprites_in_layers += pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size;
					particles_started = total_sprites > 0;
					particles_processed_last_frame |= pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size > 0;
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
				DisablePMSpawning(pm, true);
		}

		frames = tmp_frame_meta.size();

		tfxSpriteData *sprite_data = nullptr;
		if (effect->library->pre_recorded_effects.ValidKey(effect->path_hash)) {
			sprite_data = &effect->library->pre_recorded_effects.At(effect->path_hash);
			FreeSpriteData(*sprite_data);
			sprite_data->normal.frame_count = frames;
			sprite_data->normal.animation_length_in_time = frames * frame_length;
			sprite_data->normal.frame_meta.resize( frames);
			sprite_data->normal.frame_meta.zero();
		}
		else {
			tfxSpriteData data;
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

		tfxvec<tfxFrameMeta> &frame_meta = sprite_data->normal.frame_meta;
		memcpy(frame_meta.data, tmp_frame_meta.data, tmp_frame_meta.size_in_bytes());
		tmp_frame_meta.free_all();

		tfxU32 last_count = 0;
		for (auto &meta : frame_meta) {
			for (tfxEachLayer) {
				meta.index_offset[layer] = last_count;
				last_count += meta.sprite_count[layer];
			}
		}

		ClearParticleManager(pm, false);
		SetSeed(pm, anim.seed);
		preview_effect_index = AddEffectToParticleManager(pm, effect, pm->current_ebuff, 0, false, 0.f);
		SetEffectPosition(pm, preview_effect_index, tfxVec3(0.f, 0.f, 0.f));
		if (is_3d) {
			Transform3d(pm->effects.world_rotations[preview_effect_index],
				pm->effects.local_rotations[preview_effect_index],
				pm->effects.scale[preview_effect_index],
				pm->effects.world_position[preview_effect_index],
				pm->effects.local_position[preview_effect_index],
				pm->effects.translation[preview_effect_index],
				pm->effects.matrix[preview_effect_index],
				pm->effects.world_rotations[preview_effect_index],
				pm->effects.scale[preview_effect_index],
				pm->effects.world_position[preview_effect_index],
				pm->effects.matrix[preview_effect_index]
			);
		}
		else {
			Transform2d(pm->effects.world_rotations[preview_effect_index],
				pm->effects.local_rotations[preview_effect_index],
				pm->effects.scale[preview_effect_index],
				pm->effects.world_position[preview_effect_index],
				pm->effects.local_position[preview_effect_index],
				pm->effects.translation[preview_effect_index],
				pm->effects.matrix[preview_effect_index],
				pm->effects.world_rotations[preview_effect_index],
				pm->effects.scale[preview_effect_index],
				pm->effects.world_position[preview_effect_index],
				pm->effects.matrix[preview_effect_index]);
		}

		if (total_sprites == 0) {
			return;
		}

		sprite_data->normal.total_sprites = total_sprites;
		tfxSoABuffer temp_sprites_buffer;
		tfxSpriteDataSoA temp_sprites;
		if (is_3d) {
			InitSpriteData3dSoA(&sprite_data->real_time_sprites_buffer, &sprite_data->real_time_sprites, total_sprites);
			InitSpriteData3dSoA(&temp_sprites_buffer, &temp_sprites, 100);
		}
		else {
			InitSpriteData2dSoA(&sprite_data->real_time_sprites_buffer, &sprite_data->real_time_sprites, total_sprites);
			InitSpriteData2dSoA(&temp_sprites_buffer, &temp_sprites, 100);
		}
		sprite_data->normal.total_memory_for_sprites = total_sprites * (tfxU32)sprite_data->real_time_sprites_buffer.struct_size; 

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
		DisablePMSpawning(pm, false);
		total_sprites = 0;
		tfxU32 captured_offset[tfxLAYERS] = { 0, 0, 0, 0 };

		while (frame < frames && offset < 99999) {
			tfxU32 count_this_frame = 0;
			UpdateParticleManager(pm, frame_length);
			InvalidateNewSpriteCapturedIndex(pm);
			bool particles_processed_last_frame = false;

			if (offset >= start_frame) {
				for (tfxEachLayer) {
					tfxU32 meta_count = frame_meta[frame].sprite_count[layer];
					tfxU32 pm_count = pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size;
					if (running_count[layer][frame] > 0 && pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size > 0) {
						Resize(&temp_sprites_buffer, running_count[layer][frame]);
						memcpy(temp_sprites.alignment, sprite_data->real_time_sprites.alignment + frame_meta[frame].index_offset[layer], sizeof(tfxU32) * running_count[layer][frame]);
						memcpy(temp_sprites.captured_index, sprite_data->real_time_sprites.captured_index + frame_meta[frame].index_offset[layer], sizeof(tfxU32) * running_count[layer][frame]);
						memcpy(temp_sprites.uid, sprite_data->real_time_sprites.uid + frame_meta[frame].index_offset[layer], sizeof(tfxUniqueSpriteID) * running_count[layer][frame]);
						memcpy(temp_sprites.color, sprite_data->real_time_sprites.color + frame_meta[frame].index_offset[layer], sizeof(tfxU32) * running_count[layer][frame]);
						memcpy(temp_sprites.property_indexes, sprite_data->real_time_sprites.property_indexes + frame_meta[frame].index_offset[layer], sizeof(tfxU32) * running_count[layer][frame]);
						memcpy(temp_sprites.intensity, sprite_data->real_time_sprites.intensity + frame_meta[frame].index_offset[layer], sizeof(float) * running_count[layer][frame]);
						memcpy(temp_sprites.stretch, sprite_data->real_time_sprites.stretch + frame_meta[frame].index_offset[layer], sizeof(float) * running_count[layer][frame]);
						if (is_3d) {
							memcpy(temp_sprites.transform_3d, sprite_data->real_time_sprites.transform_3d + frame_meta[frame].index_offset[layer], sizeof(tfxSpriteTransform3d) * running_count[layer][frame]);
						}
						else {
							memcpy(temp_sprites.transform_2d, sprite_data->real_time_sprites.transform_2d + frame_meta[frame].index_offset[layer], sizeof(tfxSpriteTransform2d) * running_count[layer][frame]);
						}
						if (captured_offset[layer] > 0) {
							for (int temp_i = 0; temp_i != temp_sprites_buffer.current_size; ++temp_i) {
								if(temp_sprites.captured_index[temp_i] != tfxINVALID)
									temp_sprites.captured_index[temp_i] += captured_offset[layer];
							}
						}
					}
					else if (captured_offset[layer] > 0 && pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size == 0) {
						for (int index = SpriteDataIndexOffset(sprite_data, frame, layer); index != SpriteDataEndIndex(sprite_data, frame, layer); ++index) {
							if(sprite_data->real_time_sprites.captured_index[index] != tfxINVALID)
								sprite_data->real_time_sprites.captured_index[index] += captured_offset[layer];
						}
					}

					memcpy(sprite_data->real_time_sprites.alignment + frame_meta[frame].index_offset[layer], pm->sprites[pm->current_sprite_buffer][layer].alignment, sizeof(tfxU32) * pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size);
					memcpy(sprite_data->real_time_sprites.captured_index + frame_meta[frame].index_offset[layer], pm->sprites[pm->current_sprite_buffer][layer].captured_index, sizeof(tfxU32) * pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size);
					memcpy(sprite_data->real_time_sprites.uid + frame_meta[frame].index_offset[layer], pm->sprites[pm->current_sprite_buffer][layer].uid, sizeof(tfxUniqueSpriteID) * pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size);
					memcpy(sprite_data->real_time_sprites.color + frame_meta[frame].index_offset[layer], pm->sprites[pm->current_sprite_buffer][layer].color, sizeof(tfxU32) * pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size);
					memcpy(sprite_data->real_time_sprites.property_indexes + frame_meta[frame].index_offset[layer], pm->sprites[pm->current_sprite_buffer][layer].property_indexes, sizeof(tfxU32) * pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size);
					memcpy(sprite_data->real_time_sprites.intensity + frame_meta[frame].index_offset[layer], pm->sprites[pm->current_sprite_buffer][layer].intensity, sizeof(float) * pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size);
					memcpy(sprite_data->real_time_sprites.stretch + frame_meta[frame].index_offset[layer], pm->sprites[pm->current_sprite_buffer][layer].stretch, sizeof(float) * pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size);
					if (is_3d) {
						memcpy(sprite_data->real_time_sprites.transform_3d + frame_meta[frame].index_offset[layer], pm->sprites[pm->current_sprite_buffer][layer].transform_3d, sizeof(tfxSpriteTransform3d) * pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size);
					}
					else {
						memcpy(sprite_data->real_time_sprites.transform_2d + frame_meta[frame].index_offset[layer], pm->sprites[pm->current_sprite_buffer][layer].transform_2d, sizeof(tfxSpriteTransform2d) * pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size);
					}

					if (running_count[layer][frame] > 0 && pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size > 0) {
						memcpy(sprite_data->real_time_sprites.alignment + frame_meta[frame].index_offset[layer] + pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size, temp_sprites.alignment, sizeof(tfxU32) * temp_sprites_buffer.current_size);
						memcpy(sprite_data->real_time_sprites.captured_index + frame_meta[frame].index_offset[layer] + pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size, temp_sprites.captured_index, sizeof(tfxU32) * temp_sprites_buffer.current_size);
						memcpy(sprite_data->real_time_sprites.uid + frame_meta[frame].index_offset[layer] + pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size, temp_sprites.uid, sizeof(tfxUniqueSpriteID) * temp_sprites_buffer.current_size);
						memcpy(sprite_data->real_time_sprites.color + frame_meta[frame].index_offset[layer] + pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size, temp_sprites.color, sizeof(tfxU32) * temp_sprites_buffer.current_size);
						memcpy(sprite_data->real_time_sprites.property_indexes + frame_meta[frame].index_offset[layer] + pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size, temp_sprites.property_indexes, sizeof(tfxU32) * temp_sprites_buffer.current_size);
						memcpy(sprite_data->real_time_sprites.intensity + frame_meta[frame].index_offset[layer] + pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size, temp_sprites.intensity, sizeof(float) * temp_sprites_buffer.current_size);
						memcpy(sprite_data->real_time_sprites.stretch + frame_meta[frame].index_offset[layer] + pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size, temp_sprites.stretch, sizeof(float) * temp_sprites_buffer.current_size);
						if (is_3d) {
							memcpy(sprite_data->real_time_sprites.transform_3d + frame_meta[frame].index_offset[layer] + pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size, temp_sprites.transform_3d, sizeof(tfxSpriteTransform3d) * temp_sprites_buffer.current_size);
						}
						else {
							memcpy(sprite_data->real_time_sprites.transform_2d + frame_meta[frame].index_offset[layer] + pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size, temp_sprites.transform_2d, sizeof(tfxSpriteTransform2d) * temp_sprites_buffer.current_size);
						}
						captured_offset[layer] = pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size;
					}
					else if (pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size == 0) {
						captured_offset[layer] = 0;
					}
					running_count[layer][frame] += pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size;
					total_sprites += pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size;
					particles_started = total_sprites > 0;
					particles_processed_last_frame |= pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size > 0;
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
				DisablePMSpawning(pm, true);
		}

		/*
		for (int i = 0; i != total_sprites; ++i) {
			tfxU32 property_index = sprite_data->real_time_sprites.property_indexes[i] & 0x0000FFFF;
			if (property_index == 0) {
				int d = 0;
			}
		}
		*/

		sprite_data->real_time_sprites_buffer.current_size = total_sprites;
		//total sprites should not exceed the capacity of the sprite buffer
		assert(sprite_data->real_time_sprites_buffer.current_size <= sprite_data->real_time_sprites_buffer.capacity);
		ResetSpriteDataLerpOffset(*sprite_data);
		tfxSpriteDataSoA &sprites = sprite_data->real_time_sprites;

		for (int i = 0; i != anim.real_frames; ++i) {
			for (tfxEachLayer) {
				for (int j = SpriteDataIndexOffset(sprite_data, i, layer); j != SpriteDataEndIndex(sprite_data, i, layer); ++j) {
					if (sprites.captured_index[j] == tfxINVALID) continue;
					int frame = i - 1;
					frame = frame < 0 ? anim.real_frames - 1 : frame;
					tfxU32 wrap_bit = sprites.captured_index[j] & 0x10000000;
					sprites.captured_index[j] = (sprites.captured_index[j] & 0x0FFFFFFF) + sprite_data->normal.frame_meta[frame].index_offset[layer];
					sprites.captured_index[j] |= wrap_bit;
				}
			}
		}

		WrapSingleParticleSprites(sprite_data);
		ClearWrapBit(sprite_data);

		FreeSoABuffer(&temp_sprites_buffer);
		DisablePMSpawning(pm, false);
		ClearParticleManager(pm, false);
		pm->flags &= ~tfxEffectManagerFlags_recording_sprites;

		if (anim.playback_speed < 1.f) {
			CompressSpriteData(pm, effect, is_3d, frame_length);
		}
		else {
			sprite_data->compressed_sprites = sprite_data->real_time_sprites;
			sprite_data->compressed_sprites_buffer = sprite_data->real_time_sprites_buffer;
			sprite_data->compressed = sprite_data->normal;
			anim.frames_after_compression = sprite_data->normal.frame_count;
		}

		pm->camera_position = pm_camera_position;
		ToggleSpritesWithUID(pm, false);
	}

	void CompressSpriteData(tfxParticleManager *pm, tfxEffectEmitter *effect, bool is_3d, float frame_length) {
		tfxSpriteDataSettings &anim = effect->library->sprite_data_settings[GetEffectInfo(effect)->sprite_data_settings_index];
		tfxSpriteData *sprite_data = &effect->library->pre_recorded_effects.At(effect->path_hash);
		if (is_3d) {
			InitSpriteData3dSoACompression(&sprite_data->compressed_sprites_buffer, &sprite_data->compressed_sprites, sprite_data->real_time_sprites_buffer.current_size);
		}
		else {
			InitSpriteData2dSoACompression(&sprite_data->compressed_sprites_buffer, &sprite_data->compressed_sprites, sprite_data->real_time_sprites_buffer.current_size);
		}

		sprite_data->compressed.frame_meta.resize(tfxU32((float)anim.real_frames * anim.playback_speed) + 1);
		sprite_data->compressed.frame_meta.zero();

		float frequency = frame_length * (1.f / anim.playback_speed);
		float real_time = 0.f;
		float compressed_time = 0.f;
		int compressed_frame = 0;

		//First pass, add the sprites from the real time sprite data to the compressed data
		int f = 0;
		tfxSpriteDataSoA &sprites = sprite_data->real_time_sprites;
		tfxSpriteDataSoA &c_sprites = sprite_data->compressed_sprites;
		int ci = 0;
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
				for (int i = SpriteDataIndexOffset(sprite_data, f, layer); i != SpriteDataEndIndex(sprite_data, f, layer); ++i) {
					//Add to compress sprites but make invalid captured indexed create the offset
					sprite_data->compressed.frame_meta[compressed_frame].sprite_count[layer]++;
					sprite_data->compressed.frame_meta[compressed_frame].total_sprites++;
					c_sprites.alignment[ci] = sprites.alignment[i];
					c_sprites.captured_index[ci] = sprites.captured_index[i];
					c_sprites.uid[ci] = sprites.uid[i];
					c_sprites.lerp_offset[ci] = sprites.captured_index[i] == tfxINVALID ? lerp_offset : 1.f;
					c_sprites.color[ci] = sprites.color[i];
					c_sprites.property_indexes[ci] = sprites.property_indexes[i];
					c_sprites.intensity[ci] = sprites.intensity[i];
					c_sprites.stretch[ci] = sprites.stretch[i];
					if (is_3d) {
						c_sprites.transform_3d[ci] = sprites.transform_3d[i];
					}
					else {
						c_sprites.transform_2d[ci] = sprites.transform_2d[i];
					}
					ci++;
				}
				frame_done = true;
				f++;
			}
			else if (real_time > compressed_time && real_time < next_compressed_time && f < anim.real_frames) {
				for (int i = SpriteDataIndexOffset(sprite_data, f, layer); i != SpriteDataEndIndex(sprite_data, f, layer); ++i) {
					if (sprites.captured_index[i] == tfxINVALID) {
						//Add to compressed sprites frame but add the lerp offset
						sprite_data->compressed.frame_meta[compressed_frame].sprite_count[layer]++;
						sprite_data->compressed.frame_meta[compressed_frame].total_sprites++;
						c_sprites.alignment[ci] = sprites.alignment[i];
						c_sprites.captured_index[ci] = tfxINVALID;
						c_sprites.uid[ci] = sprites.uid[i];
						c_sprites.lerp_offset[ci] = lerp_offset;
						c_sprites.color[ci] = sprites.color[i];
						c_sprites.property_indexes[ci] = sprites.property_indexes[i];
						c_sprites.intensity[ci] = sprites.intensity[i];
						c_sprites.stretch[ci] = sprites.stretch[i];
						if (is_3d) {
							c_sprites.transform_3d[ci] = sprites.transform_3d[i];
						}
						else {
							c_sprites.transform_2d[ci] = sprites.transform_2d[i];
						}
						ci++;
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
		} while(!finished);

		sprite_data->compressed.total_sprites = ci;
		sprite_data->compressed_sprites_buffer.current_size = ci;
		sprite_data->compressed.total_memory_for_sprites = ci * (tfxU32)sprite_data->compressed_sprites_buffer.struct_size;

		f = 0;
		//Second pass, link up the captured indexes using the UIDs
		sprite_data->compressed.frame_count = compressed_frame + 1;
		anim.animation_length_in_time = sprite_data->compressed.animation_length_in_time = sprite_data->compressed.frame_count * frequency;
		anim.frames_after_compression = sprite_data->compressed.frame_count;
		tmpMTStack(tfxCompressWorkEntry, compress_entry);
		while (f < (int)sprite_data->compressed.frame_count) {
			tfxCompressWorkEntry *entry = &compress_entry.next();
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
		compress_entry.free();

		TrimSoABuffer(&sprite_data->compressed_sprites_buffer);
		sprite_data->compressed.total_memory_for_sprites = (tfxU32)sprite_data->compressed_sprites_buffer.current_arena_size;

	}

	void LinkUpSpriteCapturedIndexes(tfxWorkQueue *queue, void *work_entry) {
		tfxCompressWorkEntry *entry = static_cast<tfxCompressWorkEntry*>(work_entry);
		tfxSpriteData *sprite_data = static_cast<tfxSpriteData*>(entry->sprite_data);
		tfxSpriteDataSoA &c_sprites = sprite_data->compressed_sprites;

		int frame = entry->frame - 1;
		int frame_pair[2];
		frame_pair[0] = entry->frame;
		frame_pair[1] = frame < 0 ? sprite_data->compressed.frame_count - 1 : frame;
		for (tfxEachLayer) {
			for (int i = sprite_data->compressed.frame_meta[frame_pair[0]].index_offset[layer]; i != sprite_data->compressed.frame_meta[frame_pair[0]].index_offset[layer] + sprite_data->compressed.frame_meta[frame_pair[0]].sprite_count[layer]; ++i) {
				if (c_sprites.captured_index[i] != tfxINVALID) {
					tfxU32 age_diff = 0xFFFFFFFF;
					for (int j = sprite_data->compressed.frame_meta[frame_pair[1]].index_offset[layer]; j != sprite_data->compressed.frame_meta[frame_pair[1]].index_offset[layer] + sprite_data->compressed.frame_meta[frame_pair[1]].sprite_count[layer]; ++j) {
						if (c_sprites.uid[j].uid == c_sprites.uid[i].uid) {
							tfxU32 diff = c_sprites.uid[i].age - c_sprites.uid[j].age;
							age_diff = diff < age_diff ? diff : age_diff;
							c_sprites.captured_index[i] = age_diff == diff ? j : c_sprites.captured_index[i];
							if (age_diff < 2) break;
						}
					}
				}
			}
		}
	}

	void WrapSingleParticleSprites(tfxSpriteData *sprite_data) {
		tfxSpriteDataSoA &sprites = sprite_data->real_time_sprites;
		for (tfxEachLayer) {
			for (int i = sprite_data->normal.frame_meta[0].index_offset[layer]; i != sprite_data->normal.frame_meta[0].index_offset[layer] + sprite_data->normal.frame_meta[0].sprite_count[layer]; ++i) {
				if (sprites.captured_index[i] != tfxINVALID && sprites.captured_index[i] & 0x10000000) {
					for (int j = sprite_data->normal.frame_meta[sprite_data->normal.frame_count - 1].index_offset[layer]; j != sprite_data->normal.frame_meta[sprite_data->normal.frame_count - 1].index_offset[layer] + sprite_data->normal.frame_meta[sprite_data->normal.frame_count - 1].sprite_count[layer]; ++j) {
						if (sprites.uid[j].uid == sprites.uid[i].uid) {
							sprites.captured_index[i] = j;
						}
					}
				}
			}
		}
	}

	void ClearWrapBit(tfxSpriteData *sprite_data) {
		tfxSpriteDataSoA &sprites = sprite_data->real_time_sprites;
		for (int i = 0; i != sprite_data->normal.frame_count; ++i) {
			for (tfxEachLayer) {
				for (int j = SpriteDataIndexOffset(sprite_data, i, layer); j != SpriteDataEndIndex(sprite_data, i, layer); ++j) {
					if (sprites.captured_index[j] == tfxINVALID) continue;
					tfxU32 wrap_bit = 0x10000000;
					sprites.captured_index[j] &= ~wrap_bit;
				}
			}
		}
	}

	void InitialiseAnimationManagerFor3d(tfxAnimationManager *animation_manager, tfxU32 max_instances, tfxU32 initial_sprite_data_capacity) {
		animation_manager->instances.reserve(max_instances);
		animation_manager->free_instances.reserve(max_instances);
		animation_manager->render_queue.reserve(max_instances);
		animation_manager->offsets.reserve(max_instances);
		animation_manager->instances_in_use[0].reserve(max_instances);
		animation_manager->instances_in_use[1].reserve(max_instances);
		animation_manager->sprite_data_3d.reserve(initial_sprite_data_capacity);
		animation_manager->current_in_use_buffer = 0;
		animation_manager->buffer_metrics.instances_size = 0;
		animation_manager->buffer_metrics.instances_size_in_bytes = 0;
		animation_manager->buffer_metrics.offsets_size = 0;
		animation_manager->buffer_metrics.offsets_size_in_bytes = 0;
		animation_manager->buffer_metrics.sprite_data_size = 0;
		animation_manager->buffer_metrics.total_sprites_to_draw = 0;
		animation_manager->update_frequency = 60.f;
		animation_manager->flags = tfxAnimationManagerFlags_initialised | tfxAnimationManagerFlags_is_3d;
	}

	void InitialiseAnimationManagerFor2d(tfxAnimationManager *animation_manager, tfxU32 max_instances, tfxU32 initial_sprite_data_capacity) {
		animation_manager->instances.reserve(max_instances);
		animation_manager->free_instances.reserve(max_instances);
		animation_manager->render_queue.reserve(max_instances);
		animation_manager->offsets.reserve(max_instances);
		animation_manager->instances_in_use[0].reserve(max_instances);
		animation_manager->instances_in_use[1].reserve(max_instances);
		animation_manager->sprite_data_2d.reserve(initial_sprite_data_capacity);
		animation_manager->current_in_use_buffer = 0;
		animation_manager->buffer_metrics.instances_size = 0;
		animation_manager->buffer_metrics.instances_size_in_bytes = 0;
		animation_manager->buffer_metrics.offsets_size = 0;
		animation_manager->buffer_metrics.offsets_size_in_bytes = 0;
		animation_manager->buffer_metrics.sprite_data_size = 0;
		animation_manager->buffer_metrics.total_sprites_to_draw = 0;
		animation_manager->update_frequency = 60.f;
		animation_manager->flags = tfxAnimationManagerFlags_initialised;
	}

	tfxAnimationID AddAnimationInstance(tfxAnimationManager *animation_manager) {
		if (animation_manager->free_instances.current_size > 0) {
			tfxU32 index = animation_manager->free_instances.pop_back();
			animation_manager->instances_in_use[animation_manager->current_in_use_buffer].push_back(index);
			return index;
		}
		tfxAnimationInstance instance;
		tfxU32 index = animation_manager->instances.current_size;
		assert(animation_manager->instances.capacity != animation_manager->instances.current_size);		//At capacity! not enough room to add another instance.
		animation_manager->instances.push_back(instance);
		animation_manager->instances_in_use[animation_manager->current_in_use_buffer].push_back(index);
		return index;
	}

	void FreeAnimationInstance(tfxAnimationManager *animation_manager, tfxU32 index) {
		animation_manager->free_instances.push_back(index);
	}

	void AddEffectEmitterProperties(tfxAnimationManager *animation_manager, tfxEffectEmitter *effect, bool *has_animated_shape) {
		if (effect->type != tfxEmitterType) {
			for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
				AddEffectEmitterProperties(animation_manager, &sub, has_animated_shape);
			}
		}
		else {
			if (effect->library->emitter_properties.animation_property_index[effect->property_index] != tfxINVALID) {
				tfxAnimationEmitterProperties properties;
				properties.handle = effect->library->emitter_properties.image_handle[effect->property_index];
				properties.flags = effect->property_flags;
				tfxImageData &image = *effect->library->emitter_properties.image[effect->property_index];
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
				effect->library->emitter_properties.animation_property_index[effect->property_index] = animation_manager->emitter_properties.current_size;
				properties.image_ptr = image.ptr;
				animation_manager->emitter_properties.push_back(properties);
				for (auto &sub : GetEffectInfo(effect)->sub_effectors) {
					AddEffectEmitterProperties(animation_manager, &sub, has_animated_shape);
				}
			}
		}
	}

	void AddEffectShapes(tfxAnimationManager *animation_manager, tfxEffectEmitter *effect) {
		if (effect->type == tfxEmitterType) {
			tfxImageData *image_data = GetEffectProperties(effect)->image[effect->property_index];
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

	void AddSpriteData(tfxAnimationManager *animation_manager, tfxEffectEmitter *effect, tfxParticleManager *pm, tfxVec3 camera_position) {
		if (Is3DEffect(effect)) {
			//If you're adding 3d effect sprite data then the animation manager must have been initialised with InitialiseAnimationManagerFor3d
			assert(animation_manager->flags & tfxAnimationManagerFlags_is_3d);
		}
		else {
			//If you're adding 2d effect sprite data then the animation manager must have been initialised with InitialiseAnimationManagerFor2d
			assert(!(animation_manager->flags & tfxAnimationManagerFlags_is_3d));
		}
		tfxSpriteDataSettings &anim = effect->library->sprite_data_settings[GetEffectInfo(effect)->sprite_data_settings_index];
		if (!effect->library->pre_recorded_effects.ValidKey(effect->path_hash)) {
			assert(pm);		//You must pass an appropriate particle manager if the animation needs recording
			RecordSpriteData(pm, effect, animation_manager->update_frequency, &camera_position.x);
		}

		bool has_animated_shape = false;
		AddEffectEmitterProperties(animation_manager, effect, &has_animated_shape);

		bool is_3d = Is3DEffect(effect);
		tfxSpriteData &sprite_data = effect->library->pre_recorded_effects.At(effect->path_hash);
		animation_manager->effect_animation_info.Insert(effect->path_hash, sprite_data.compressed);
		tfxSpriteDataMetrics &metrics = animation_manager->effect_animation_info.At(effect->path_hash);
		metrics.name = GetEffectInfo(effect)->name;
		metrics.frames_after_compression = anim.frames_after_compression;
		metrics.real_frames = anim.real_frames;
		metrics.animation_length_in_time = anim.animation_length_in_time;
		metrics.animation_flags = anim.animation_flags;
		metrics.flags = has_animated_shape ? tfxAnimationManagerFlags_has_animated_shapes : 0;
		tfxSpriteDataSoA &sprites = sprite_data.compressed_sprites;
		if (is_3d) {
			metrics.start_offset = animation_manager->sprite_data_3d.current_size;
			for (int i = 0; i != metrics.total_sprites; ++i) {
				tfxSpriteData3d sprite;
				sprite.alignment = sprites.alignment[i];
				sprite.captured_index = sprites.captured_index[i];
				sprite.captured_index += sprite.captured_index == tfxINVALID ? 0 : metrics.start_offset;
				sprite.color = sprites.color[i];
				sprite.property_indexes = sprites.property_indexes[i];
				tfxU32 property_index = sprite.property_indexes & 0x0000FFFF;
				//tfxImageData &image = *effect->library->emitter_properties.image[property_index];
				//Temporary while debugging:
				//sprite.lookup_indexes = image.compute_shape_index + ((sprites.property_indexes[i] & 0x00FF0000) >> 16);
				//sprite.lookup_indexes += (sprite.property_indexes & 0x0000FFFF) << 16;
				//-------------------------
				sprite.property_indexes &= ~0x0000FFFF;
				sprite.property_indexes += effect->library->emitter_properties.animation_property_index[property_index];
				tfxU32 alignment = (sprite.property_indexes & 0xFF000000) >> 24;
				sprite.intensity = sprites.intensity[i];
				sprite.lerp_offset = sprites.lerp_offset[i];
				sprite.stretch = sprites.stretch[i];
				sprite.position = sprites.transform_3d[i].position;
				sprite.rotations = sprites.transform_3d[i].rotations;
				sprite.scale = sprites.transform_3d[i].scale;
				animation_manager->sprite_data_3d.push_back(sprite);
			}
			metrics.total_memory_for_sprites = sizeof(tfxSpriteData3d) * metrics.total_sprites;
		}
		else {
			metrics.start_offset = animation_manager->sprite_data_2d.current_size;
			for (int i = 0; i != metrics.total_sprites; ++i) {
				tfxSpriteData2d sprite;
				sprite.alignment = sprites.alignment[i];
				sprite.captured_index = sprites.captured_index[i];
				sprite.captured_index += sprite.captured_index == tfxINVALID ? 0 : metrics.start_offset;
				sprite.color = sprites.color[i];
				sprite.property_indexes = sprites.property_indexes[i];
				tfxU32 property_index = sprite.property_indexes & 0x0000FFFF;
				//tfxImageData &image = *effect->library->emitter_properties.image[property_index];
				//Temporary while debugging:
				//sprite.lookup_indexes = image.compute_shape_index + ((sprites.property_indexes[i] & 0x00FF0000) >> 16);
				//sprite.lookup_indexes += (sprite.property_indexes & 0x0000FFFF) << 16;
				//-------------------------
				sprite.property_indexes &= ~0x0000FFFF;
				sprite.property_indexes += effect->library->emitter_properties.animation_property_index[property_index];
				tfxU32 alignment = (sprite.property_indexes & 0xFF000000) >> 24;
				sprite.intensity = sprites.intensity[i];
				sprite.lerp_offset = sprites.lerp_offset[i];
				sprite.stretch = sprites.stretch[i];
				sprite.position = sprites.transform_2d[i].position;
				sprite.rotation = sprites.transform_2d[i].rotation;
				sprite.scale = sprites.transform_2d[i].scale;
				animation_manager->sprite_data_2d.push_back(sprite);
			}
			metrics.total_memory_for_sprites = sizeof(tfxSpriteData2d) * metrics.total_sprites;
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

	tfxAnimationID AddAnimationInstance(tfxAnimationManager *animation_manager, tfxKey path, tfxU32 start_frame) {
		assert(animation_manager->effect_animation_info.ValidKey(path));				//You must have added the effect sprite data to the animation manager
																						//Call AddSpriteData to do so
		if(animation_manager->instances_in_use->current_size >= animation_manager->instances_in_use->capacity) {
			return tfxINVALID;
		}
		tfxU32 info_index = animation_manager->effect_animation_info.GetIndex(path);
		tfxSpriteDataMetrics &metrics = animation_manager->effect_animation_info.data[info_index];
		assert(start_frame < metrics.frames_after_compression);
		tfxAnimationID index = AddAnimationInstance(animation_manager);
		tfxAnimationInstance &instance = animation_manager->instances[index];
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

	tfxAnimationID AddAnimationInstance(tfxAnimationManager *animation_manager, const char *path, tfxU32 start_frame) {
		tfxKey path_hash = tfxXXHash64::hash(path, strlen(path), 0);
		return AddAnimationInstance(animation_manager, path_hash, start_frame);
	}

	void UpdateAnimationManager(tfxAnimationManager *animation_manager, float elapsed) {
		assert(animation_manager->instances_in_use[animation_manager->current_in_use_buffer].capacity > 0);	//You must call InitialiseAnimationManager before trying to update one
		tfxU32 next_buffer = !animation_manager->current_in_use_buffer;
		animation_manager->instances_in_use[next_buffer].clear();
		animation_manager->render_queue.clear();
		animation_manager->offsets.clear();
		tfxU32 running_sprite_count = 0;
		animation_manager->flags &= ~tfxAnimationManagerFlags_has_animated_shapes;
		for (auto i : animation_manager->instances_in_use[animation_manager->current_in_use_buffer]) {
			auto &instance = animation_manager->instances[i];
			tfxSpriteDataMetrics &metrics = animation_manager->effect_animation_info.data[instance.info_index];
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
				animation_manager->render_queue.push_back(instance);
				running_sprite_count += instance.sprite_count;
				animation_manager->flags |= metrics.flags & tfxAnimationManagerFlags_has_animated_shapes;
				animation_manager->offsets.push_back(running_sprite_count);
				UpdateAnimationManagerBufferMetrics(animation_manager);
			}
		}
		animation_manager->buffer_metrics.total_sprites_to_draw = running_sprite_count;
		animation_manager->current_in_use_buffer = !animation_manager->current_in_use_buffer;
	}

	void CycleAnimationManager(tfxAnimationManager *animation_manager) {
		assert(animation_manager->instances_in_use[animation_manager->current_in_use_buffer].capacity > 0);	//You must call InitialiseAnimationManager before trying to update one
		tfxU32 next_buffer = !animation_manager->current_in_use_buffer;
		animation_manager->instances_in_use[next_buffer].clear();
		animation_manager->render_queue.clear();
		animation_manager->offsets.clear();
		tfxU32 running_sprite_count = 0;
		animation_manager->flags &= ~tfxAnimationManagerFlags_has_animated_shapes;
		for (auto i : animation_manager->instances_in_use[animation_manager->current_in_use_buffer]) {
			auto &instance = animation_manager->instances[i];
			tfxSpriteDataMetrics &metrics = animation_manager->effect_animation_info.data[instance.info_index];
			float frame_time = (instance.current_time / instance.animation_length_in_time) * (float)instance.frame_count;
			tfxU32 frame = tfxU32(frame_time);
			frame++;
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
		animation_manager->current_in_use_buffer = !animation_manager->current_in_use_buffer;
	}

	void ClearAllAnimationInstances(tfxAnimationManager *animation_manager) {
		animation_manager->free_instances.clear();
		animation_manager->instances_in_use[0].clear();
		animation_manager->instances_in_use[1].clear();
		animation_manager->render_queue.clear();
		animation_manager->instances.clear();
		animation_manager->offsets.clear();
	}

	void ResetAnimationManager(tfxAnimationManager *animation_manager) {
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

		animation_manager->buffer_metrics.instances_size = 0;
		animation_manager->buffer_metrics.instances_size_in_bytes = 0;
		animation_manager->buffer_metrics.offsets_size = 0;
		animation_manager->buffer_metrics.offsets_size_in_bytes = 0;
		animation_manager->buffer_metrics.sprite_data_size = 0;
		animation_manager->buffer_metrics.total_sprites_to_draw = 0;
	}

	void FreeAnimationManager(tfxAnimationManager *animation_manager) {
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

	void UpdateAnimationManagerBufferMetrics(tfxAnimationManager *animation_manager) {
		animation_manager->buffer_metrics.instances_size = animation_manager->render_queue.current_size;
		animation_manager->buffer_metrics.offsets_size = animation_manager->offsets.current_size;
		animation_manager->buffer_metrics.instances_size_in_bytes = animation_manager->buffer_metrics.instances_size * sizeof(tfxAnimationInstance) * animation_manager->buffer_metrics.instances_size;
		animation_manager->buffer_metrics.offsets_size_in_bytes = animation_manager->buffer_metrics.offsets_size * sizeof(tfxU32) * animation_manager->buffer_metrics.instances_size;
	}

	void RecordTemplateEffect(tfxEffectTemplate *t, tfxParticleManager *pm, float update_frequency, float camera_position[3]) {
		RecordSpriteData(pm, &t->effect, update_frequency, camera_position);
	}

	void DisableTemplateEmitter(tfxEffectTemplate *t, const char *path) {
		assert(t->paths.ValidName(path));			//Must be a valid path to the emitter
		tfxEffectEmitter *emitter = t->paths.At(path);
		assert(emitter->type == tfxEmitterType);	//Must be an emitter that you're trying to remove. Use RemoveSubEffect if you're trying to remove one of those. 
		emitter->property_flags &= ~tfxEmitterPropertyFlags_enabled;
	}

	void EnableTemplateEmitter(tfxEffectTemplate *t, const char *path) {
		assert(t->paths.ValidName(path));			//Must be a valid path to the emitter
		tfxEffectEmitter *emitter = t->paths.At(path);
		assert(emitter->type == tfxEmitterType);	//Must be an emitter that you're trying to remove. Use RemoveSubEffect if you're trying to remove one of those
		emitter->property_flags |= tfxEmitterPropertyFlags_enabled;
	}

	void ScaleTemplateGlobalMultiplier(tfxEffectTemplate *t, tfxGraphType global_type, float amount) {
		assert(IsGlobalGraph(global_type));
		tfxGraph *graph = GetEffectGraphByType(&t->effect, global_type);
		tfxEffectEmitter *original_effect = GetLibraryEffect(t->effect.library, t->original_effect_hash);
		tfxGraph *original_graph = GetEffectGraphByType(original_effect, global_type);
		CopyGraph(original_graph, graph, false);
		MultiplyAllGraphValues(graph, amount);
		CompileGraph(*graph);
	}

	void ScaleTemplateEmitterGraph(tfxEffectTemplate *t, const char *emitter_path, tfxGraphType graph_type, float amount) {
		assert(IsEmitterGraph(graph_type));		//Must be an emitter graph type. This is any property, base, variaion or overtime graph
		assert(t->paths.ValidName(emitter_path));			//Must be a valid path to the emitter
		tfxEffectEmitter *emitter = t->paths.At(emitter_path);
		tfxGraph *graph = GetEffectGraphByType(emitter, graph_type);
		tfxEffectEmitter *original_emitter = GetLibraryEffect(t->effect.library, emitter_path);
		tfxGraph *original_graph = GetEffectGraphByType(original_emitter, graph_type);
		CopyGraph(original_graph, graph, false);
		MultiplyAllGraphValues(graph, amount);
		CompileGraph(*graph);
	}

	void SetTemplateSingleSpawnAmount(tfxEffectTemplate *t, const char *emitter_path, tfxU32 amount) {
		assert(amount >= 0);							//Amount must not be less than 0
		assert(t->paths.ValidName(emitter_path));			//Must be a valid path to the emitter
		tfxEffectEmitter *emitter = t->paths.At(emitter_path);
		GetEffectProperties(emitter)->spawn_amount[emitter->property_index] = amount;
	}

	void* GetAnimationEmitterPropertiesBufferPointer(tfxAnimationManager *animation_manager) {
		return animation_manager->emitter_properties.data;
	}

	void ResetTemplate(tfxEffectTemplate *t) {
		if (t->paths.Size()) {
			t->paths.Clear();
			CleanUpEffect(&t->effect);
		}
	}

	tfxEffectEmitter *GetEffectFromTemplate(tfxEffectTemplate *t) {
		return &t->effect;
	}

	tfxEffectEmitter *GetEmitterFromTemplate(tfxEffectTemplate *t, tfxStr256 &path) {
		if (t->paths.ValidName(path)) return t->paths.At(path); return nullptr;
	}

	void SetTemplateUserData(tfxEffectTemplate *t, tfxStr256 &path, void *data) {
		if (t->paths.ValidName(path)) t->paths.At(path)->user_data = data;
	}

	void SetTemplateEffectUserData(tfxEffectTemplate *t, void *data) {
		t->effect.user_data = data;
	}

	void SetTemplateUserDataAll(tfxEffectTemplate *t, void *data);

	void SetTemplateEffectUpdateCallback(tfxEffectTemplate *t, void(*update_callback)(tfxParticleManager *pm, tfxEffectID effect_index)) {
		t->effect.update_callback = update_callback;
	}

	tfxParticleManager::~tfxParticleManager() {
		FreeSoABuffer(&effect_buffers);
		FreeSoABuffer(&emitter_buffers);
		particle_array_allocator.FreeAll();
	}

	tfxEffectID AddEffectToParticleManager(tfxParticleManager *pm, tfxEffectTemplate *effect_template) {
		return AddEffectToParticleManager(pm, &effect_template->effect, pm->current_ebuff, 0, false, 0.f);
	}

	tfxEffectID AddEffectToParticleManager(tfxParticleManager *pm, tfxEffectEmitter *effect) {
		return AddEffectToParticleManager(pm, effect, pm->current_ebuff, 0, false, 0.f);
	}

	tfxEffectID AddEffectToParticleManager(tfxParticleManager *pm, tfxEffectEmitter *effect, int buffer, int hierarchy_depth, bool is_sub_emitter, float add_delayed_spawning) {
		tfxPROFILE;
		assert(effect->type == tfxEffectType);
		assert(effect->library == pm->library);	//The effect must belong to the same library that is assigned to the particle manager
		if (pm->flags & tfxEffectManagerFlags_use_compute_shader && pm->highest_compute_controller_index >= pm->max_compute_controllers && pm->free_compute_controllers.empty())
			return tfxINVALID;
		unsigned int parent_index = GetPMEffectSlot(pm);
		if (parent_index == tfxINVALID)
			return tfxINVALID;
		if (!is_sub_emitter) {
			pm->effects.highest_particle_age[parent_index] = pm->frame_length * 3.f;
		}
		tfxEmitterPropertiesSoA &properties = effect->library->emitter_properties;
		pm->effects.global_attributes[parent_index] = effect->global;
		pm->effects.transform_attributes[parent_index] = effect->transform_attributes;
		pm->effects.age[parent_index] = -add_delayed_spawning;
		pm->effects.state_flags[parent_index] = 0;
		pm->effects.frame[parent_index] = 0.f;
		pm->effects.property_flags[parent_index] = effect->property_flags;
		pm->effects.local_position[parent_index] = tfxVec3();
		pm->effects.timeout[parent_index] = 100.f;
		pm->effects.library[parent_index] = effect->library;
		pm->effects.parent_particle_index[parent_index] = tfxINVALID;
		pm->effects.info_index[parent_index] = effect->info_index;
		pm->effects.properties_index[parent_index] = effect->property_index;
		pm->effects.timeout_counter[parent_index] = 0;
		pm->effects.user_data[parent_index] = effect->user_data;
		pm->effects.update_callback[parent_index] = effect->update_callback;
		float range = properties.noise_base_offset_range[effect->property_index];
		pm->effects.noise_base_offset[parent_index] = pm->random.Range(range);
		pm->effects_in_use[hierarchy_depth][buffer].push_back(parent_index);
		effect->pm_index = parent_index;
		pm->sort_passes = tfxMax(effect->sort_passes, pm->sort_passes);
		pm->sort_passes = tfxMin(5, pm->sort_passes);

		tfxU32 seed_index = 0;
		for (auto &e : GetEffectInfo(effect)->sub_effectors) {
			if (e.property_flags & tfxEmitterPropertyFlags_enabled) {
				unsigned int index = GetPMEmitterSlot(pm);
				if (index == tfxINVALID)
					break;
				pm->emitters.particles_index[index] = tfxINVALID;
				pm->emitters_in_use[hierarchy_depth][buffer].push_back(index);
				pm->emitters.parent_index[index] = parent_index;
				if (pm->emitters.particles_index[index] == tfxINVALID) {
					if (!is_sub_emitter)
						pm->emitters.particles_index[index] = GrabParticleLists(*pm, e.path_hash, 100);
				}
				pm->emitters.path_hash[index] = e.path_hash;
				pm->emitters.info_index[index] = e.info_index;
				pm->emitters.properties_index[index] = e.property_index;
				pm->emitters.emitter_attributes[index] = e.emitter_attributes;
				pm->emitters.transform_attributes[index] = e.transform_attributes;
				pm->emitters.delay_spawning[index] = properties.delay_spawning[e.property_index];
				pm->emitters.age[index] = 0.f;
				pm->emitters.frame[index] = 0.f;
				pm->emitters.local_position[index] = tfxVec3();
				pm->emitters.grid_coords[index] = tfxVec3();
				pm->emitters.grid_direction[index] = tfxVec3();
				pm->emitters.property_flags[index] = e.property_flags;
				pm->emitters.image_size[index] = properties.image[e.property_index]->image_size;
				pm->emitters.image_frame_rate[index] = properties.image[e.property_index]->animation_frames > 1 && e.property_flags & tfxEmitterPropertyFlags_animate ? properties.frame_rate[e.property_index] : 0.f;
				//pm->emitters.image_frame_rate[index] = e.property_flags & tfxEmitterPropertyFlags_reverse_animation ? -pm->emitters.image_frame_rate[index] : pm->emitters.image_frame_rate[index];
				pm->emitters.end_frame[index] = properties.end_frame[e.property_index];
				pm->emitters.angle_offsets[index] = properties.angle_offsets[e.property_index];
				pm->emitters.timeout[index] = 2000.f;
				pm->emitters.amount_remainder[index] = 0.f;
				pm->emitters.qty_step_size[index] = 0.f;
				pm->emitters.timeout_counter[index] = 0;
				pm->emitters.emitter_size[index] = 0.f;
				pm->emitters.hierarchy_depth[index] = hierarchy_depth;
				pm->emitters.world_rotations[index] = 0.f;
				pm->emitters.seed_index[index] = seed_index++;
				e.pm_index = index;		//Doesn't have much use beyond the editor?
				//----Handle
				if (e.property_flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
					pm->emitters.image_handle[index] = tfxVec2(0.5f, 0.5f);
				}
				else {
					pm->emitters.image_handle[index] = properties.image_handle[e.property_index];
				}
				tfxEmitterStateFlags &state_flags = pm->emitters.state_flags[index];
				const tfxEmitterStateFlags &parent_state_flags = pm->emitters.state_flags[parent_index];

				state_flags = tfxEmitterStateFlags_no_tween_this_update;
				state_flags &= ~tfxEmitterStateFlags_retain_matrix;
				state_flags |= parent_state_flags & tfxEffectStateFlags_no_tween;
				state_flags |= e.property_flags & tfxEmitterPropertyFlags_single && !(pm->flags & tfxEffectManagerFlags_disable_spawning) ? tfxEmitterStateFlags_is_single : 0;
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
				state_flags |= GetGraphMaxValue(&e.library->emitter_attributes[e.emitter_attributes].overtime.velocity_turbulance) > 0 ? tfxEmitterStateFlags_has_noise : 0;

				if (state_flags & tfxEmitterStateFlags_is_line_traversal) {
					pm->emitters.property_flags[index] |= tfxEmitterPropertyFlags_relative_position;
				}

				if (is_sub_emitter) {
					state_flags |= tfxEmitterStateFlags_is_sub_emitter;
				}
				else {
					pm->emitters.highest_particle_age[index] = pm->frame_length * 2.f;
				}

				/*if (pm->flags & tfxEffectManagerFlags_use_compute_shader && GetEffectInfo(e)->sub_effectors.empty()) {
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
		pm->emitters.state_flags[parent_index] |= tfxEmitterStateFlags_no_tween_this_update;
		return parent_index;
	}

	int AddComputeController(tfxParticleManager *pm) {
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

	void ResetParticlePtr(tfxParticleManager *pm, void *ptr) {
		pm->new_compute_particle_ptr = ptr;
		pm->new_compute_particle_index = 0;
	}

	void ResetControllerPtr(tfxParticleManager *pm, void *ptr) {
		pm->compute_controller_ptr = ptr;
	}

	void UpdateCompute(tfxParticleManager *pm, void *sampled_particles, unsigned int sample_size) {
		for (int i = 0; i != sample_size; ++i) {
			if (pm->compute_global_state.current_length == 0)
				break;
			tfxComputeParticle *sample = static_cast<tfxComputeParticle*>(sampled_particles) + i;
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

	tfxComputeParticle *GrabComputeParticle(tfxParticleManager *pm, unsigned int layer) {
		assert(pm->new_compute_particle_ptr);		//Use must assign the compute ptr to point to an area in memory where you can stage new particles for uploading to the GPU - See ResetComputePtr
		return (static_cast<tfxComputeParticle*>(pm->new_compute_particle_ptr) + pm->new_compute_particle_index++);
	}

	void FreeParticleList(tfxParticleManager *pm, tfxU32 index) {
		if (pm->free_particle_lists.ValidKey(pm->emitters.path_hash[index])) {
			pm->free_particle_lists.At(pm->emitters.path_hash[index]).push_back(pm->emitters.particles_index[index]);
		}
		else {
#ifdef tfxTRACK_MEMORY
			tfxvec<tfxU32> new_indexes(tfxCONSTRUCTOR_VEC_INIT("FreeParticleBanks::new_indexes"));
#else
			tfxvec<tfxU32> new_indexes;
#endif
			new_indexes.push_back(pm->emitters.particles_index[index]);
			pm->free_particle_lists.Insert(pm->emitters.path_hash[index], new_indexes);
		}
	}

	void UpdateParticleManager(tfxParticleManager *pm, float elapsed_time) {
		tfxPROFILE;

		tfxCompleteAllWork(&pm->work_queue);

		pm->frame_length = elapsed_time;
		pm->frame_length_wide = tfxWideSetSingle(pm->frame_length);
		pm->update_frequency = 1000.f / elapsed_time;
		pm->update_time = 1.f / pm->update_frequency;
		pm->update_time_wide = tfxWideSetSingle(pm->update_time);
		pm->new_compute_particle_index = 0;

		tfxU32 next_buffer = !pm->current_ebuff;
		tfxU32 depth_starting_index[tfxLAYERS];
		pm->current_sprite_buffer = pm->flags & tfxEffectManagerFlags_double_buffer_sprites ? !pm->current_sprite_buffer : 0;

		memset(pm->sprite_index_point, 0, sizeof(tfxU32) * tfxLAYERS);

		for (tfxEachLayer) {
			depth_starting_index[layer] = pm->depth_indexes[layer][pm->current_depth_index_buffer].current_size;
		}

		for (tfxEachLayer) {
			ClearSoABuffer(&pm->sprite_buffer[pm->current_sprite_buffer][layer]);
		}

		tfxU32 effects_start_size[tfxMAXDEPTH];
		tfxU32 emitter_start_size[tfxMAXDEPTH];
		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			effects_start_size[depth] = pm->effects_in_use[depth][pm->current_ebuff].current_size;
			emitter_start_size[depth] = pm->emitters_in_use[depth][pm->current_ebuff].current_size;
		}

		tmpMTStack(tfxSpawnWorkEntry, spawn_work);
		//Loop over all the effects and emitters, depth by depth, and add spawn jobs to the worker queue
		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			pm->effects_in_use[depth][next_buffer].clear();
			pm->emitters_in_use[depth][next_buffer].clear();

			for (int i = 0; i != effects_start_size[depth]; ++i) {
				tfxU32 current_index = pm->effects_in_use[depth][pm->current_ebuff][i];
				float &timeout_counter = pm->effects.timeout_counter[current_index];

				UpdatePMEffect(*pm, current_index);
				if (timeout_counter <= pm->effects.timeout[current_index]) {
					pm->effects_in_use[depth][next_buffer].push_back(current_index);
				}
				else {
					pm->free_effects.push_back(current_index);
				}
			}

			for (int i = 0; i != emitter_start_size[depth]; ++i) {
				tfxSpawnWorkEntry *spawn_work_entry = &spawn_work.next();
				tfxU32 current_index = pm->emitters_in_use[depth][pm->current_ebuff][i];
				spawn_work_entry->depth = depth;
				spawn_work_entry->emitter_index = current_index;
				spawn_work_entry->next_buffer = next_buffer;
				spawn_work_entry->properties = &pm->library->emitter_properties;
				spawn_work_entry->sub_effects = &pm->library->effect_infos[pm->emitters.info_index[current_index]].sub_effectors;
				spawn_work_entry->amount_to_spawn = 0;
				spawn_work_entry->end_index = 0;
				spawn_work_entry->highest_particle_age = 0;

				float &timeout_counter = pm->emitters.timeout_counter[current_index];

				UpdatePMEmitter(*pm, spawn_work_entry);
				if (timeout_counter <= pm->emitters.timeout[current_index]) {
					pm->emitters_in_use[depth][next_buffer].push_back(current_index);
				}
				else {
					pm->free_emitters.push_back(current_index);
					//if (pm->flags & tfxEffectManagerFlags_use_compute_shader && pm->emitters.property_flags[current_index] & tfxEmitterPropertyFlags_is_bottom_emitter)
						//FreeComputeSlot(pm->emitters.compute_slot_id[current_index]);
					if (pm->flags & tfxEffectManagerFlags_unordered) {
						FreeParticleList(pm, current_index);
					}
				}
			}
		}

		tfxCompleteAllWork(&pm->work_queue);
		pm->random.Advance();

		for (auto &work_entry : spawn_work) {
			tfxU32 index = work_entry.emitter_index;
			pm->emitters.highest_particle_age[index] = std::fmaxf(pm->emitters.highest_particle_age[index], work_entry.highest_particle_age);
			pm->effects.highest_particle_age[pm->emitters.parent_index[index]] = pm->emitters.highest_particle_age[index] + pm->frame_length;
		}
		spawn_work.free();

		if (!(pm->flags & tfxEffectManagerFlags_unordered)) {
			for (tfxEachLayer) {
				if (depth_starting_index[layer] < pm->depth_indexes[layer][pm->current_depth_index_buffer].current_size) {
					tfxU32 next_depth_buffer = !pm->current_depth_index_buffer;
					if (pm->depth_indexes[layer][next_depth_buffer].capacity < pm->depth_indexes[layer][pm->current_depth_index_buffer].capacity) {
						pm->depth_indexes[layer][next_depth_buffer].reserve(pm->depth_indexes[layer][pm->current_depth_index_buffer].capacity);
					}
					if (pm->flags & tfxEffectManagerFlags_order_by_depth) {
						//No need to qsort ordered by age as the depth with all be 0 (depth is particle age)
						std::qsort(&pm->depth_indexes[layer][pm->current_depth_index_buffer][depth_starting_index[layer]], pm->depth_indexes[layer][pm->current_depth_index_buffer].current_size - depth_starting_index[layer], sizeof(tfxDepthIndex), SortDepth);
					}
					tfxU32 current_depth_index = 0;
					tfxU32 second_index = depth_starting_index[layer];
					for (auto &depth_index : pm->depth_indexes[layer][pm->current_depth_index_buffer]) {
						if (depth_starting_index[layer] != 0) {
							while (second_index < pm->depth_indexes[layer][pm->current_depth_index_buffer].current_size && depth_index.depth < pm->depth_indexes[layer][pm->current_depth_index_buffer][second_index].depth) {
								pm->particle_arrays[ParticleBank(pm->depth_indexes[layer][pm->current_depth_index_buffer][second_index].particle_id)].depth_index[ParticleIndex(pm->depth_indexes[layer][pm->current_depth_index_buffer][second_index].particle_id)] = pm->depth_indexes[layer][next_depth_buffer].current_size;
								pm->depth_indexes[layer][next_depth_buffer].push_back(pm->depth_indexes[layer][pm->current_depth_index_buffer][second_index++]);
							}
						}
						pm->particle_arrays[ParticleBank(depth_index.particle_id)].depth_index[ParticleIndex(depth_index.particle_id)] = pm->depth_indexes[layer][next_depth_buffer].current_size;
						pm->depth_indexes[layer][next_depth_buffer].push_back(depth_index);
						if (++current_depth_index == depth_starting_index[layer])
							break;
					}
					if (depth_starting_index[layer] != 0 && second_index < pm->depth_indexes[layer][pm->current_depth_index_buffer].current_size) {
						while (second_index < pm->depth_indexes[layer][pm->current_depth_index_buffer].current_size) {
							tfxU32 bank = ParticleBank(pm->depth_indexes[layer][pm->current_depth_index_buffer][second_index].particle_id);
							tfxU32 index = ParticleIndex(pm->depth_indexes[layer][pm->current_depth_index_buffer][second_index].particle_id);
							pm->particle_arrays[ParticleBank(pm->depth_indexes[layer][pm->current_depth_index_buffer][second_index].particle_id)].depth_index[ParticleIndex(pm->depth_indexes[layer][pm->current_depth_index_buffer][second_index].particle_id)] = pm->depth_indexes[layer][next_depth_buffer].current_size;
							pm->depth_indexes[layer][next_depth_buffer].push_back(pm->depth_indexes[layer][pm->current_depth_index_buffer][second_index++]);
						}
					}
					assert(pm->depth_indexes[layer][next_depth_buffer].current_size == pm->depth_indexes[layer][pm->current_depth_index_buffer].current_size);
					pm->depth_indexes[layer][pm->current_depth_index_buffer].clear();
					pm->current_depth_index_buffer = next_depth_buffer;
				}
			}
		}

		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			tmpMTStack(tfxControlWorkEntry, work);
			for (int index : pm->emitters_in_use[depth][next_buffer]) {
				tfxSoABuffer &bank = pm->particle_array_buffers[pm->emitters.particles_index[index]];
				int particles_to_update = bank.current_size;
				tfxU32 running_start_index = 0;
				while (particles_to_update > 0) {
					tfxControlWorkEntry &work_entry = work.next();
					work_entry.properties = &pm->library->emitter_properties;
					work_entry.emitter_index = index;
					work_entry.start_index = running_start_index;
					work_entry.end_index = particles_to_update > pm->mt_batch_size ? running_start_index + pm->mt_batch_size : running_start_index + particles_to_update;
					tfxU32 circular_start = GetCircularIndex(&pm->particle_array_buffers[pm->emitters.particles_index[index]], work_entry.start_index);
					tfxU32 block_start_index = (circular_start / tfxDataWidth) * tfxDataWidth;
					work_entry.wide_end_index = (tfxU32)(ceilf((float)work_entry.end_index / tfxDataWidth)) * tfxDataWidth;
					work_entry.start_diff = circular_start - block_start_index;
					work_entry.wide_end_index = work_entry.wide_end_index - work_entry.start_diff < work_entry.end_index ? work_entry.wide_end_index + tfxDataWidth : work_entry.wide_end_index;
					particles_to_update -= pm->mt_batch_size;
					running_start_index += pm->mt_batch_size;
					ControlParticles(*pm, index, work_entry);
				}
			}

			tfxCompleteAllWork(&pm->work_queue);
			work.free();
			if(pm->emitters_check_capture.current_size > 0) {
				tmpMTStack(tfxControlWorkEntry, work);
				for (int index : pm->emitters_check_capture) {
					//Really don't like this but is fine for now. For any line emitters where the particles loop back round to the beginning we need to set the captured index of the sprites
					//so that they don't interpolate the frame that they loop
					tfxSoABuffer &bank = pm->particle_array_buffers[pm->emitters.particles_index[index]];
					int particles_to_update = bank.current_size;
					tfxU32 running_start_index = 0;
					while (particles_to_update > 0) {
						tfxControlWorkEntry &work_entry = work.next();
						work_entry.properties = &pm->library->emitter_properties;
						work_entry.emitter_index = index;
						work_entry.start_index = running_start_index;
						work_entry.end_index = particles_to_update > pm->mt_batch_size ? running_start_index + pm->mt_batch_size : running_start_index + particles_to_update;
						tfxU32 circular_start = GetCircularIndex(&pm->particle_array_buffers[pm->emitters.particles_index[index]], work_entry.start_index);
						tfxU32 block_start_index = (circular_start / tfxDataWidth) * tfxDataWidth;
						work_entry.wide_end_index = (tfxU32)(ceilf((float)work_entry.end_index / tfxDataWidth)) * tfxDataWidth;
						work_entry.start_diff = circular_start - block_start_index;
						work_entry.wide_end_index = work_entry.wide_end_index - work_entry.start_diff < work_entry.end_index ? work_entry.wide_end_index + tfxDataWidth : work_entry.wide_end_index;
						particles_to_update -= pm->mt_batch_size;
						running_start_index += pm->mt_batch_size;
						work_entry.pm = pm;
						const tfxU32 property_index = pm->emitters.properties_index[index];
						const tfxU32 sprites_index = pm->emitters.sprites_index[index];
						tfxEmitterPropertiesSoA &properties = pm->library->emitter_properties;
						work_entry.sprites_index = sprites_index + work_entry.start_index;
						work_entry.sprite_buffer_end_index = work_entry.sprites_index + (work_entry.end_index - work_entry.start_index);
						work_entry.layer = properties.layer[property_index];
						work_entry.sprites = &pm->sprites[pm->current_sprite_buffer][work_entry.layer];
						if (!(pm->flags & tfxEffectManagerFlags_single_threaded) && tfxNumberOfThreadsInAdditionToMain) {
							tfxAddWorkQueueEntry(&pm->work_queue, &work_entry, ControlParticleCaptureFlag);
						}
						else {
							ControlParticleCaptureFlag(&pm->work_queue, &work_entry);
						}
					}
				}
				tfxCompleteAllWork(&pm->work_queue);
				work.free();
				pm->emitters_check_capture.clear();
			}

			{
				tmpMTStack(tfxParticleAgeWorkEntry, work);
				for (int index : pm->emitters_in_use[depth][next_buffer]) {
					tfxSoABuffer &bank = pm->particle_array_buffers[pm->emitters.particles_index[index]];
					tfxParticleAgeWorkEntry &work_entry = work.next();
					work_entry.properties = &pm->library->emitter_properties;
					work_entry.start_index = bank.current_size - 1;
					work_entry.emitter_index = index;
					tfxU32 circular_start = GetCircularIndex(&pm->particle_array_buffers[pm->emitters.particles_index[index]], 0);
					tfxU32 block_start_index = (circular_start / tfxDataWidth) * tfxDataWidth;
					work_entry.wide_end_index = (tfxU32)(ceilf((float)bank.current_size / tfxDataWidth)) * tfxDataWidth;
					work_entry.start_diff = circular_start - block_start_index;
					work_entry.wide_end_index += work_entry.wide_end_index - work_entry.start_diff < bank.current_size ? tfxDataWidth : 0;
					work_entry.pm = pm;
					if (!(pm->flags & tfxEffectManagerFlags_single_threaded) && tfxNumberOfThreadsInAdditionToMain) {
						tfxAddWorkQueueEntry(&pm->work_queue, &work_entry, ControlParticleAge);
					}
					else {
						ControlParticleAge(&pm->work_queue, &work_entry);
					}
				}

				tfxCompleteAllWork(&pm->work_queue);
				work.free();
			}
		}

		//Todo work queue this for each layer
		if (!(pm->flags & tfxEffectManagerFlags_unordered)) {
			for (tfxEachLayer) {
				for (auto &depth_index : pm->depth_indexes[layer][pm->current_depth_index_buffer]) {
					if (depth_index.particle_id != tfxINVALID) {
						pm->particle_arrays[ParticleBank(depth_index.particle_id)].depth_index[ParticleIndex(depth_index.particle_id)] = pm->depth_indexes[layer][!pm->current_depth_index_buffer].current_size;
						pm->depth_indexes[layer][!pm->current_depth_index_buffer].push_back(depth_index);
					}
				}
				//if(pm->depth_indexes[layer][!pm->current_depth_index_buffer].current_size)
					//assert(pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size == pm->depth_indexes[layer][!pm->current_depth_index_buffer].current_size);
			}
		}

		tfxCompleteAllWork(&pm->work_queue);

		for (tfxEachLayer) {
			pm->depth_indexes[layer][pm->current_depth_index_buffer].clear();
		}
		pm->current_depth_index_buffer = !pm->current_depth_index_buffer;

		if (pm->flags & tfxEffectManagerFlags_order_by_depth && pm->flags & tfxEffectManagerFlags_3d_effects) {
			if (pm->flags & tfxEffectManagerFlags_guarantee_order) {
				for (tfxEachLayer) {
					tfxSortWorkEntry &work_entry = pm->sorting_work_entry[layer];
					work_entry.bank = &pm->particle_arrays;
					work_entry.depth_indexes = &pm->depth_indexes[layer][pm->current_depth_index_buffer];
					if (!(pm->flags & tfxEffectManagerFlags_single_threaded) && tfxNumberOfThreadsInAdditionToMain > 0) {
						tfxAddWorkQueueEntry(&pm->work_queue, &work_entry, InsertionSortDepth);
					}
					else {
						InsertionSortDepth(&pm->work_queue, &work_entry);
					}
				}
			}
			else if (pm->sort_passes > 0) {
				for (tfxEachLayer) {
					tfxvec<tfxDepthIndex> &depth_index = pm->depth_indexes[layer][pm->current_depth_index_buffer];
					//Add this to a work queue
					for (tfxU32 sorts = 0; sorts != pm->sort_passes; ++sorts) {
						for (tfxU32 i = 1; i < depth_index.current_size; ++i) {
							float depth1 = depth_index[i - 1].depth;
							float depth2 = depth_index[i].depth;
							if (depth1 < depth2) {
								pm->particle_arrays[ParticleBank(depth_index[i].particle_id)].depth_index[ParticleIndex(depth_index[i].particle_id)] = i - 1;
								pm->particle_arrays[ParticleBank(depth_index[i - 1].particle_id)].depth_index[ParticleIndex(depth_index[i - 1].particle_id)] = i;
								std::swap(depth_index[i], depth_index[i - 1]);
							}
						}
					}
				}
			}
		}

		//Add Subeffects to the next buffer
		for (int depth = 1; depth != tfxMAXDEPTH; ++depth) {
			for (int i = effects_start_size[depth]; i != pm->effects_in_use[depth][pm->current_ebuff].current_size; ++i) {
				tfxU32 current_index = pm->effects_in_use[depth][pm->current_ebuff][i];
				pm->effects_in_use[depth][next_buffer].push_back(current_index);
			}
			for (int i = emitter_start_size[depth]; i != pm->emitters_in_use[depth][pm->current_ebuff].current_size; ++i) {
				tfxU32 current_index = pm->emitters_in_use[depth][pm->current_ebuff][i];
				pm->emitters.particles_index[current_index] = GrabParticleLists(*pm, pm->emitters.path_hash[current_index], 100);
				pm->emitters_in_use[depth][next_buffer].push_back(current_index);
			}
		}

		pm->current_ebuff = next_buffer;

		pm->flags &= ~tfxEffectManagerFlags_update_base_values;

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
			current_velocity_y = tfxWideSub(current_velocity_y, tfxWideMul(base_weight, lookup_weight));
			current_velocity_x = tfxWideMul(tfxWideMul(current_velocity_x, pm.update_time_wide), velocity_adjuster);
			current_velocity_y = tfxWideMul(tfxWideMul(current_velocity_y, pm.update_time_wide), velocity_adjuster);
			current_velocity_z = tfxWideMul(tfxWideMul(current_velocity_z, pm.update_time_wide), velocity_adjuster);

			//----Spin and angle Changes
			if (emitter_flags & tfxEmitterStateFlags_can_spin) {
				roll.m = tfxWideAdd(roll.m, tfxWideMul(lookup_spin, pm.update_time_wide));
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
			tfxWideArray capture_flag;
			capture_flag.m = tfxWideCast(tfxWideGreateri(tfxWideAndi(flags, capture_after_transform), tfxWideSetZeroi()));
			tfxWideFloat xor_capture_flag = tfxWideEquals(capture_flag.m, tfxWideSetZero());
			_ReadBarrier();

			if ((tfxU32)capture_flag.a[0] & tfxParticleFlags_capture_after_transform) {
				int d = 0;
			}

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

			captured_position_x.m = tfxWideAdd(tfxWideAnd(position_x.m, capture_flag.m), tfxWideAnd(captured_position_x.m, xor_capture_flag));
			captured_position_y.m = tfxWideAdd(tfxWideAnd(position_y.m, capture_flag.m), tfxWideAnd(captured_position_y.m, xor_capture_flag));
			captured_position_z.m = tfxWideAdd(tfxWideAnd(position_z.m, capture_flag.m), tfxWideAnd(captured_position_z.m, xor_capture_flag));

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
				p_stretch.m = tfxWideMul(p_stretch.m, tfxWideDiv(l, pm.update_time_wide));	//This is too arbitrary, think up a better solution!
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
			packed.m = PackWide10bit(alignment_vector_x, alignment_vector_y, alignment_vector_z);

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

			//flags = tfxWideAndi(flags, xor_capture_after_transform_flag);
			//tfxWideStorei((tfxWideInt*)&bank.flags[index], flags);

			start_diff = 0;
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
			const tfxWideFloat base_size_y = tfxWideLoad(&bank.base_size_y[index]);
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
			tfxWideFloat stretch_velocity_x;
			tfxWideFloat stretch_velocity_y;
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
			stretch_velocity_x = current_velocity_x.m;
			stretch_velocity_y = current_velocity_y.m;
			current_velocity_x.m = tfxWideMul(tfxWideMul(current_velocity_x.m, pm.update_time_wide), velocity_adjuster);
			current_velocity_y.m = tfxWideMul(tfxWideMul(current_velocity_y.m, pm.update_time_wide), velocity_adjuster);

			//----Spin and angle Changes
			if (emitter_flags & tfxEmitterStateFlags_can_spin) {
				roll.m = tfxWideAdd(roll.m, tfxWideMul(lookup_spin, pm.update_time_wide));
			}

			//----Position
			local_position_x = tfxWideAdd(local_position_x, tfxWideMul(current_velocity_x.m, overal_scale_wide));
			local_position_y = tfxWideAdd(local_position_y, tfxWideMul(current_velocity_y.m, overal_scale_wide));

			tfxWideStore(&bank.position_x[index], local_position_x);
			tfxWideStore(&bank.position_y[index], local_position_y);
			tfxWideStore(&bank.local_rotations_z[index], roll.m);

			stretch_velocity_y = tfxWideAdd(stretch_velocity_y, tfxWideSetSingle(0.000001f));
			tfxWideFloat l = tfxWideMul(stretch_velocity_x, stretch_velocity_x);
			l = tfxWideAdd(l, tfxWideMul(stretch_velocity_y, stretch_velocity_y));
			l = tfxWideSqrt(l);
			p_stretch.m = tfxWideMul(p_stretch.m, tfxWideMul(l, tfxWideSetSingle(0.02f)));
			stretch_velocity_x = tfxWideDiv(stretch_velocity_x, l);
			stretch_velocity_y = tfxWideDiv(stretch_velocity_y, l);

			if (property_flags & tfxEmitterPropertyFlags_relative_position) {
				mmWideTransformVector(matrix, stretch_velocity_x, stretch_velocity_y);
			}

			tfxWideArrayi packed;
			packed.m = PackWide16bit(stretch_velocity_x, stretch_velocity_y);

			tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
			if (!(pm.flags & tfxEffectManagerFlags_unordered)) {	//Predictable
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 sprite_depth_index = bank.depth_index[index + j];
					sprites.stretch[sprite_depth_index] = p_stretch.a[j];
					sprites.alignment[sprite_depth_index] = packed.a[j];
					running_sprite_index++;
				}
			}
			else {
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

		tfxWideFloat image_frame_rate = tfxWideSetSingle(pm.emitters.image_frame_rate[emitter_index]);
		image_frame_rate = tfxWideMul(image_frame_rate, pm.update_time_wide);
		tfxWideFloat end_frame = tfxWideSetSingle(pm.emitters.end_frame[emitter_index]);
		tfxWideFloat frames = tfxWideSetSingle(pm.emitters.end_frame[emitter_index] + 1);
		tfxEmitterStateFlags emitter_flags = pm.emitters.state_flags[emitter_index];
		tfxEmitterStateFlags property_flags = pm.emitters.property_flags[emitter_index];

		tfxU32 running_sprite_index = work_entry->sprites_index;
		tfxSpriteSoA &sprites = *work_entry->sprites;

		for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
			tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;

			tfxWideArray image_frame;
			image_frame.m = tfxWideLoad(&bank.image_frame[index]);

			//----Image animation
			image_frame.m = tfxWideAdd(image_frame.m, image_frame_rate);
			tfxWideStore(&bank.image_frame[index], image_frame.m);
			if (emitter_flags & tfxEmitterStateFlags_play_once) {
				image_frame.m = tfxWideMin(image_frame.m, end_frame);
				image_frame.m = tfxWideMax(image_frame.m, tfxWideSetZero());
			}
			else if(property_flags & tfxEmitterPropertyFlags_reverse_animation) {
				image_frame.m = tfxWideMod(image_frame.m, frames);
				image_frame.m = tfxWideSub(end_frame, image_frame.m);
			}
			else {
				image_frame.m = tfxWideMod(image_frame.m, frames);
			}

			tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
			if (!(pm.flags & tfxEffectManagerFlags_unordered)) {				//Predictable
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 sprite_depth_index = bank.depth_index[index + j];
					tfxU32 &sprites_index = bank.sprite_index[index + j];
					float &age = bank.age[index + j];
					sprites.captured_index[sprite_depth_index] = age == 0.f && bank.single_loop_count[index + j] == 0 ? (pm.current_sprite_buffer << 30) + sprite_depth_index : (!pm.current_sprite_buffer << 30) + (sprites_index & 0x0FFFFFFF);
					sprites.captured_index[sprite_depth_index] |= property_flags & tfxEmitterPropertyFlags_wrap_single_sprite ? 0x10000000 : 0;
					sprites_index = (work_entry->layer << 28) + sprite_depth_index;
					sprites.property_indexes[sprite_depth_index] = (billboard_option << 24) + ((tfxU32)image_frame.a[j] << 16) + (property_index);
					running_sprite_index++;
				}
			}
			else {
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 &sprites_index = bank.sprite_index[index + j];
					float &age = bank.age[index + j];
					sprites.captured_index[running_sprite_index] = age == 0.f && bank.single_loop_count[index + j] == 0 ? (pm.current_sprite_buffer << 30) + running_sprite_index : (!pm.current_sprite_buffer << 30) + (sprites_index & 0x0FFFFFFF);
					sprites.captured_index[running_sprite_index] |= property_flags & tfxEmitterPropertyFlags_wrap_single_sprite ? 0x10000000 : 0;
					sprites_index = (work_entry->layer << 28) + running_sprite_index;
					sprites.property_indexes[running_sprite_index++] = (billboard_option << 24) + ((tfxU32)image_frame.a[j] << 16) + (property_index);
				}
			}

			start_diff = 0;
		}

	}

	void ControlParticleCaptureFlag(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		const tfxU32 particles_index = work_entry->pm->emitters.particles_index[emitter_index];
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = pm.particle_arrays[particles_index];

		tfxU32 start_diff = work_entry->start_diff;

		const tfxWideInt xor_capture_after_transform_flag = tfxWideXOri(tfxWideSetSinglei(tfxParticleFlags_capture_after_transform), tfxWideSetSinglei(-1));
		tfxU32 running_sprite_index = work_entry->sprites_index;
		tfxSpriteSoA &sprites = *work_entry->sprites;

		for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
			tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;
			tfxWideArrayi flags;
			flags.m = tfxWideLoadi((tfxWideInt*)&bank.flags[index]);

			tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
			if (!(pm.flags & tfxEffectManagerFlags_unordered)) {				//Predictable
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 sprite_depth_index = bank.depth_index[index + j];
					sprites.captured_index[sprite_depth_index] = (flags.a[j] & tfxParticleFlags_capture_after_transform) ? (pm.current_sprite_buffer << 30) + sprite_depth_index : sprites.captured_index[sprite_depth_index];
				}
			}
			else {
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					sprites.captured_index[running_sprite_index] = (flags.a[j] & tfxParticleFlags_capture_after_transform) ? (pm.current_sprite_buffer << 30) + running_sprite_index : sprites.captured_index[running_sprite_index];
					running_sprite_index++;
				}
			}
			flags.m = tfxWideAndi(flags.m, xor_capture_after_transform_flag);
			tfxWideStorei((tfxWideInt*)&bank.flags[index], flags.m);
			start_diff = 0;
		}
	}

	void ControlParticleUID(tfxWorkQueue *queue, void *data) {
		tfxPROFILE;
		tfxControlWorkEntry *work_entry = static_cast<tfxControlWorkEntry*>(data);
		tfxU32 emitter_index = work_entry->emitter_index;
		const tfxU32 particles_index = work_entry->pm->emitters.particles_index[emitter_index];
		tfxEmitterPropertyFlags emitter_flags = work_entry->pm->emitters.property_flags[emitter_index];
		tfxParticleManager &pm = *work_entry->pm;
		tfxParticleSoA &bank = pm.particle_arrays[particles_index];

		tfxU32 start_diff = work_entry->start_diff;

		tfxU32 running_sprite_index = work_entry->sprites_index;
		tfxSpriteSoA &sprites = *work_entry->sprites;

		for (tfxU32 i = work_entry->start_index; i != work_entry->wide_end_index; i += tfxDataWidth) {
			tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;

			tfxU32 limit_index = running_sprite_index + tfxDataWidth > work_entry->sprite_buffer_end_index ? work_entry->sprite_buffer_end_index - running_sprite_index : tfxDataWidth;
			if (!(pm.flags & tfxEffectManagerFlags_unordered)) {				//Predictable
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					tfxU32 sprite_depth_index = bank.depth_index[index + j];
					sprites.uid[sprite_depth_index].uid = bank.uid[index + j];
					sprites.uid[sprite_depth_index].age = tfxU32((bank.age[index + j] + 0.1f) / pm.frame_length);
					running_sprite_index++;
				}
			}
			else {
				for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
					sprites.uid[running_sprite_index].uid = bank.uid[index + j];
					sprites.uid[running_sprite_index++].age = tfxU32((bank.age[index + j] + 0.1f) / pm.frame_length);
				}
			}
			start_diff = 0;
		}
	}

	tfxvec<tfxU32> *GetPMEffectBuffer(tfxParticleManager *pm, tfxU32 depth) {
		return &pm->effects_in_use[depth][pm->current_ebuff];
	}

	tfxvec<tfxU32> *GetPMEmitterBuffer(tfxParticleManager *pm, tfxU32 depth) {
		return &pm->emitters_in_use[depth][pm->current_ebuff];
	}

	void ToggleSpritesWithUID(tfxParticleManager *pm, bool switch_on) {
		if (switch_on) {
			for (tfxEachLayer) {
				FreeSoABuffer(&pm->sprite_buffer[0][layer]);
				if (pm->flags & tfxEffectManagerFlags_double_buffer_sprites) {
					FreeSoABuffer(&pm->sprite_buffer[1][layer]);
				}

#ifdef tfxTRACK_MEMORY
				memcpy(sprites3d[layer].name, "ParticleManager::sprites3d\0", 27);
#endif

				if (pm->flags & tfxEffectManagerFlags_2d_and_3d) {
					InitSpriteBufferSoA(&pm->sprite_buffer[0][layer], &pm->sprites[0][layer], tfxMax((pm->max_cpu_particles_per_layer[layer] / tfxDataWidth + 1) * tfxDataWidth, 8), tfxSpriteBufferMode_both, true);
					if (pm->flags & tfxEffectManagerFlags_double_buffer_sprites) {
						InitSpriteBufferSoA(&pm->sprite_buffer[1][layer], &pm->sprites[1][layer], tfxMax((pm->max_cpu_particles_per_layer[layer] / tfxDataWidth + 1) * tfxDataWidth, 8), tfxSpriteBufferMode_both, true);
					}
				}
				else if (pm->flags & tfxEffectManagerFlags_3d_effects) {
					InitSpriteBufferSoA(&pm->sprite_buffer[0][layer], &pm->sprites[0][layer], tfxMax((pm->max_cpu_particles_per_layer[layer] / tfxDataWidth + 1) * tfxDataWidth, 8), tfxSpriteBufferMode_3d, true);
					if (pm->flags & tfxEffectManagerFlags_double_buffer_sprites) {
						InitSpriteBufferSoA(&pm->sprite_buffer[1][layer], &pm->sprites[1][layer], tfxMax((pm->max_cpu_particles_per_layer[layer] / tfxDataWidth + 1) * tfxDataWidth, 8), tfxSpriteBufferMode_3d, true);
					}
				}
				else {
					InitSpriteBufferSoA(&pm->sprite_buffer[0][layer], &pm->sprites[0][layer], tfxMax((pm->max_cpu_particles_per_layer[layer] / tfxDataWidth + 1) * tfxDataWidth, 8), tfxSpriteBufferMode_2d, true);
					if (pm->flags & tfxEffectManagerFlags_double_buffer_sprites) {
						InitSpriteBufferSoA(&pm->sprite_buffer[1][layer], &pm->sprites[1][layer], tfxMax((pm->max_cpu_particles_per_layer[layer] / tfxDataWidth + 1) * tfxDataWidth, 8), tfxSpriteBufferMode_2d, true);
					}
				}
				pm->flags |= tfxEffectManagerFlags_using_uids;
			}
		}
		else {
			pm->flags &= ~tfxEffectManagerFlags_using_uids;
		}
	}

	void ReconfigureParticleManager(tfxParticleManager *pm, tfxParticleManagerModes mode, tfxU32 req_sort_passes, bool is_3d) {
		ClearParticleManager(pm, false);
		FreeParticleBanks(pm);
		for (auto &bank : pm->free_particle_lists.data) {
			bank.free_all();
		}
		pm->free_particle_lists.FreeAll();

		tfxParticleManagerFlags current_flags = pm->flags & tfxEffectManagerFlags_dynamic_sprite_allocation | pm->flags & tfxEffectManagerFlags_double_buffer_sprites | pm->flags & tfxEffectManagerFlags_2d_and_3d;

		if (mode == tfxParticleManagerMode_unordered)
			pm->flags = tfxEffectManagerFlags_unordered;
		else if (mode == tfxParticleManagerMode_ordered_by_depth)
			pm->flags = tfxEffectManagerFlags_order_by_depth;
		else if (mode == tfxParticleManagerMode_ordered_by_depth_guaranteed)
			pm->flags = tfxEffectManagerFlags_order_by_depth | tfxEffectManagerFlags_guarantee_order;
		else if (mode == tfxParticleManagerMode_ordered_by_age)
			pm->flags = tfxEffectManagerFlags_ordered_by_age;

		if (is_3d)
			pm->flags |= tfxEffectManagerFlags_3d_effects;
		else
			pm->flags &= ~tfxEffectManagerFlags_3d_effects;


		pm->flags |= current_flags;

		for (tfxEachLayer) {
			ClearSoABuffer(&pm->sprite_buffer[0][layer]);
			if (pm->flags & tfxEffectManagerFlags_double_buffer_sprites) {
				ClearSoABuffer(&pm->sprite_buffer[1][layer]);
			}
			if (pm->flags & tfxEffectManagerFlags_ordered_by_age || pm->flags & tfxEffectManagerFlags_order_by_depth) {
				if (pm->depth_indexes[layer][0].capacity < pm->max_cpu_particles_per_layer[layer]) {
					pm->depth_indexes[layer][0].reserve(pm->max_cpu_particles_per_layer[layer]);
				}
				if (pm->depth_indexes[layer][1].capacity < pm->max_cpu_particles_per_layer[layer]) {
					pm->depth_indexes[layer][1].reserve(pm->max_cpu_particles_per_layer[layer]);
				}
				pm->depth_indexes[layer][0].clear();
				pm->depth_indexes[layer][1].clear();
			}
		}

		memset(pm->sprite_index_point, 0, 4 * tfxLAYERS);

		ClearSoABuffer(&pm->emitter_buffers);
		ClearSoABuffer(&pm->effect_buffers);

		pm->sort_passes = req_sort_passes;
	}

	void InitParticleManagerForBoth(tfxParticleManager *pm, tfxLibrary *lib, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfxParticleManagerModes mode, bool double_buffer_sprites, bool dynamic_sprite_allocation, tfxU32 multi_threaded_batch_size) {
		pm->max_effects = effects_limit;
		pm->mt_batch_size = multi_threaded_batch_size;
		pm->library = lib;

		if (pm->particle_array_allocator.arenas.current_size == 0) {
			//todo need to be able to adjust the arena size
			pm->particle_array_allocator = CreateArenaManager(tfxMegabyte(2), 8);
			pm->particle_arrays = tfxBucketArray<tfxParticleSoA>(&pm->particle_array_allocator, 32);
		}

		pm->flags = 0;

		if (mode == tfxParticleManagerMode_unordered)
			pm->flags = tfxEffectManagerFlags_unordered;
		else if (mode == tfxParticleManagerMode_ordered_by_depth)
			pm->flags = tfxEffectManagerFlags_order_by_depth;
		else if (mode == tfxParticleManagerMode_ordered_by_depth_guaranteed)
			pm->flags = tfxEffectManagerFlags_order_by_depth | tfxEffectManagerFlags_guarantee_order;
		else if (mode == tfxParticleManagerMode_ordered_by_age)
			pm->flags = tfxEffectManagerFlags_ordered_by_age;

		if (double_buffer_sprites)
			pm->flags |= tfxEffectManagerFlags_double_buffer_sprites;

		pm->flags |= tfxEffectManagerFlags_2d_and_3d;

		for (tfxEachLayer) {
			pm->max_cpu_particles_per_layer[layer] = layer_max_values[layer];

#ifdef tfxTRACK_MEMORY
			memcpy(sprites2d[layer].name, "ParticleManager::sprites2d\0", 27);
			memcpy(sprites3d[layer].name, "ParticleManager::sprites3d\0", 27);
#endif

			InitSpriteBufferSoA(&pm->sprite_buffer[0][layer], &pm->sprites[0][layer], tfxMax((layer_max_values[layer] / tfxDataWidth + 1) * tfxDataWidth, tfxDataWidth * 2), tfxSpriteBufferMode_both);
			if (pm->flags & tfxEffectManagerFlags_double_buffer_sprites) {
				InitSpriteBufferSoA(&pm->sprite_buffer[1][layer], &pm->sprites[1][layer], tfxMax((layer_max_values[layer] / tfxDataWidth + 1) * tfxDataWidth, tfxDataWidth * 2), tfxSpriteBufferMode_both);
			}

			pm->depth_indexes[layer][0].reserve(layer_max_values[layer]);
			pm->depth_indexes[layer][1].reserve(layer_max_values[layer]);
		}

		pm->flags |= dynamic_sprite_allocation ? tfxEffectManagerFlags_dynamic_sprite_allocation : 0;

		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			pm->effects_in_use[depth][0].reserve(pm->max_effects);
			pm->effects_in_use[depth][1].reserve(pm->max_effects);

			pm->emitters_in_use[depth][0].reserve(pm->max_effects);
			pm->emitters_in_use[depth][1].reserve(pm->max_effects);
		}

		pm->free_effects.reserve(pm->max_effects);
		pm->particle_indexes.reserve(1000);	//todo: Handle this better.
		InitEffectSoA(&pm->effect_buffers, &pm->effects, pm->max_effects);
		InitEmitterSoA(&pm->emitter_buffers, &pm->emitters, pm->max_effects);
	}

	void ClearParticleManager(tfxParticleManager *pm, bool free_particle_banks) {
		tfxCompleteAllWork(&pm->work_queue);
		for (tfxEachLayer) {
			ClearSoABuffer(&pm->sprite_buffer[0][layer]);
			if (pm->flags & tfxEffectManagerFlags_double_buffer_sprites) {
				ClearSoABuffer(&pm->sprite_buffer[1][layer]);
			}
			pm->depth_indexes[layer][0].clear();
			pm->depth_indexes[layer][1].clear();
		}
		if (pm->flags & tfxEffectManagerFlags_unordered) {
			if (free_particle_banks) {
				FreeParticleBanks(pm);
			}
			else {
				for (auto &bank : pm->particle_array_buffers) {
					ClearSoABuffer(&bank);
				}
			}
			for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
				for (auto index : pm->emitters_in_use[depth][pm->current_ebuff]) {
					FreeParticleList(pm, index);
				}
			}
		}
		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			pm->effects_in_use[depth][0].clear();
			pm->effects_in_use[depth][1].clear();

			pm->emitters_in_use[depth][0].clear();
			pm->emitters_in_use[depth][1].clear();
		}
		pm->free_effects.clear();
		pm->free_emitters.clear();
		pm->particle_indexes.clear();
		pm->free_particle_indexes.clear();
		ClearSoABuffer(&pm->emitter_buffers);
		ClearSoABuffer(&pm->effect_buffers);
		pm->particle_id = 0;
	}

	void FreeParticleBanks(tfxParticleManager *pm) {
		for (auto &bank : pm->particle_array_buffers) {
			FreeSoABuffer(&bank);
		}
		pm->particle_array_buffers.clear();
		pm->particle_arrays.clear();
	}

	void SoftExpireAll(tfxParticleManager *pm) {
		for (auto index : pm->effects_in_use[0][pm->current_ebuff]) {
			pm->emitters.state_flags[index] |= tfxEmitterStateFlags_stop_spawning;
		}
	}

	void SetPMLookUpMode(tfxParticleManager *pm, tfxLookupMode mode) {
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

	void UpdatePMBaseValues(tfxParticleManager *pm) {
		pm->flags |= tfxEffectManagerFlags_update_base_values;
	}

	bool FreePMEffectCapacity(tfxParticleManager *pm) {
		return pm->emitter_buffers.current_size < pm->max_effects;
	}

	tfxU32 GetPMEffectSlot(tfxParticleManager *pm) {
		if (!pm->free_effects.empty()) {
			return pm->free_effects.pop_back();
		}
		if (pm->effect_buffers.current_size == pm->effect_buffers.capacity)
			return tfxINVALID;
		AddRow(&pm->effect_buffers);
		return pm->effect_buffers.current_size - 1;
	}

	tfxU32 GetPMEmitterSlot(tfxParticleManager *pm) {
		if (!pm->free_emitters.empty()) {
			return pm->free_emitters.pop_back();
		}
		if (pm->emitter_buffers.current_size == pm->emitter_buffers.capacity) {
			return tfxINVALID;
		}
		AddRow(&pm->emitter_buffers);
		return pm->emitter_buffers.current_size - 1;
	}

	tfxU32 GetPMParticleIndexSlot(tfxParticleManager *pm, tfxParticleID particle_id) {
		//Todo: ideally we want a better thread safe container for this
		std::lock_guard<std::mutex> lock(pm->particle_index_mutex);
		if (!pm->free_particle_indexes.empty()) {
			pm->particle_indexes[pm->free_particle_indexes.back()] = particle_id;
			return pm->free_particle_indexes.pop_back();
		}
		pm->particle_indexes.push_back(particle_id);
		return pm->particle_indexes.current_size - 1;
	}

	void FreePMParticleIndex(tfxParticleManager *pm, tfxU32 &index) {
		//Todo: ideally we want a better thread safe container for this
		std::lock_guard<std::mutex> lock(pm->particle_index_mutex);
		pm->particle_indexes[index] = tfxINVALID;
		pm->free_particle_indexes.push_back(index);
		index = tfxINVALID;
	}

	tfxU32 PushPMDepthIndex(tfxParticleManager *pm, tfxU32 layer, tfxDepthIndex depth_index) {
		pm->depth_indexes[layer][pm->current_depth_index_buffer].push_back(depth_index);
		return pm->depth_indexes[layer][pm->current_depth_index_buffer].current_size - 1;
	}

	void SetPMLibrary(tfxParticleManager *pm, tfxLibrary *lib) {
		pm->library = lib;
	}

	void SetPMCamera(tfxParticleManager *pm, float front[3], float position[3]) {
		pm->camera_front.x = front[0];
		pm->camera_front.y = front[1];
		pm->camera_front.z = front[2];
		pm->camera_position.x = position[0];
		pm->camera_position.y = position[1];
		pm->camera_position.z = position[2];
	}

	void ResetPMFlags(tfxParticleManager *pm) { 
		pm->flags = 0; 
	}

	tfxU32 &GetParticleSpriteIndex(tfxParticleManager *pm, tfxParticleID id) { 
		return pm->particle_arrays[ParticleBank(id)].sprite_index[ParticleIndex(id)]; 
	}

	unsigned int GetControllerMemoryUsage(tfxParticleManager *pm) { 
		return pm->highest_compute_controller_index * sizeof(tfxComputeController); 
	}
	
	unsigned int GetParticleMemoryUsage(tfxParticleManager *pm) { 
		return pm->new_compute_particle_index * sizeof(tfxComputeParticle); 
	}

	void FreeComputeSlot(tfxParticleManager *pm, unsigned int slot_id) { 
		pm->free_compute_controllers.push_back(slot_id); 
	}

	tfxU32 ParticleCount(tfxParticleManager *pm) {
		tfxU32 count = 0;
		for (tfxEachLayer) {
			count += pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size;
		}
		return count;
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
				tfxU32 sprite_id = GetParticleSpriteIndex(&pm, parent_particle_id);
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

		age += pm.frame_length;
		highest_particle_age -= pm.frame_length;

		if (properties.loop_length && age > properties.loop_length[property_index])
			age -= properties.loop_length[property_index];

		if (highest_particle_age <= 0 && age > pm.frame_length * 5.f) {
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
		delay_spawning = -pm.frame_length;

		//e.state_flags |= e.parent->state_flags & tfxEmitterStateFlags_stop_spawning;
		state_flags |= parent_state_flags & tfxEffectStateFlags_no_tween;
		state_flags |= parent_state_flags & tfxEffectStateFlags_stop_spawning;
		state_flags |= parent_state_flags & tfxEffectStateFlags_remove;
		UpdateEmitterState(pm, index, parent_index, pm.effects.spawn_controls[parent_index], spawn_work_entry);

		if (state_flags & tfxEmitterStateFlags_is_line_loop_or_kill && state_flags & tfxEmitterStateFlags_loop) {
			pm.emitters_check_capture.push_back(index);
		}

		if (state_flags & tfxEmitterStateFlags_is_single && properties.single_shot_limit[property_index] > 1) {
			pm.emitters_check_capture.push_back(index);
		}

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
					GrowArrays(&sprite_buffer, sprite_buffer.capacity, sprite_buffer.capacity + (sprites_count + max_spawn_count - FreeSpace(&sprite_buffer)) + 1);
					if (!(pm.flags & tfxEffectManagerFlags_unordered)) {
						pm.depth_indexes[layer][pm.current_depth_index_buffer].reserve(sprite_buffer.capacity);
					}
				}
			}
			else {
				sprites_count = sprites_count > FreeSpace(&sprite_buffer) ? FreeSpace(&sprite_buffer) : sprites_count;
				max_spawn_count = max_spawn_count > FreeSpace(&sprite_buffer) ? FreeSpace(&sprite_buffer) : max_spawn_count;
			}
			sprite_buffer.current_size += max_spawn_count + sprites_count;

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
					GrowArrays(&sprite_buffer, sprite_buffer.capacity, sprite_buffer.capacity + (sprites_count + max_spawn_count - FreeSpace(&sprite_buffer)) + 1);
					if (!(pm.flags & tfxEffectManagerFlags_unordered)) {
						pm.depth_indexes[layer][pm.current_depth_index_buffer].reserve(sprite_buffer.capacity);
					}
				}
			}
			else {
				sprites_count = sprites_count > FreeSpace(&sprite_buffer) ? FreeSpace(&sprite_buffer) : sprites_count;
				max_spawn_count = max_spawn_count > FreeSpace(&sprite_buffer) ? FreeSpace(&sprite_buffer) : max_spawn_count;
			}
			sprite_buffer.current_size += max_spawn_count + sprites_count;

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

		age += pm.frame_length;
		if (!(property_flags & tfxEmitterPropertyFlags_single) || (property_flags & tfxEmitterPropertyFlags_single && properties.single_shot_limit[property_index] > 0) || state_flags & tfxEmitterStateFlags_stop_spawning) {
			highest_particle_age -= pm.frame_length;
			spawn_work_entry->highest_particle_age -= pm.frame_length;
		}

		if (properties.loop_length && age > properties.loop_length[property_index])
			age -= properties.loop_length[property_index];

		if (highest_particle_age <= 0 && age > pm.frame_length * 5.f) {
			timeout_counter += pm.frame_length;
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

			spawn_quantity *= pm.update_time;
			step_size = 1.f / spawn_quantity;
		}
		else if (property_flags & tfxEmitterPropertyFlags_match_amount_to_grid_points) {
			if (property_flags & tfxEmitterPropertyFlags_match_amount_to_grid_points && property_flags & tfxEmitterPropertyFlags_spawn_on_grid) {
				float x = tfxMax(properties.grid_points[property_index].x, 1.f);
				float y = tfxMax(properties.grid_points[property_index].y, 1.f);
				float z = tfxMax(properties.grid_points[property_index].z, 1.f);
				if (properties.emission_type[property_index] == tfxArea) {
				}
				else if (properties.emission_type[property_index] == tfxCylinder) {
					spawn_quantity = x * y;
				}
				switch (properties.emission_type[property_index]) {
				case tfxEmissionType::tfxArea:
					if (property_flags & tfxEmitterPropertyFlags_is_3d) {
						if (property_flags & tfxEmitterPropertyFlags_fill_area) {
							spawn_quantity = x * y * z;
						}
						else if (property_flags & tfxEmitterPropertyFlags_area_open_ends) {
							spawn_quantity = x * z * 2 + y * z * 2 - 4 * z;
						}
						else {
							spawn_quantity = x * z * 2 + (y - 2) * (x * 2 + z * 2 - 4);
						}
					}
					else {
						if (property_flags & tfxEmitterPropertyFlags_fill_area) {
							spawn_quantity = x * y;
						}
						else {
							spawn_quantity = x * 2 + (y - 2) * 2;
						}
					}
					break;
				case tfxEmissionType::tfxCylinder:
					spawn_quantity = x * y;
					break;
				case tfxEmissionType::tfxEllipse:
					if (!(property_flags & tfxEmitterPropertyFlags_is_3d)) {
						spawn_quantity = x;
					}
					break;
				case tfxEmissionType::tfxLine:
					spawn_quantity = x;
					break;
				case tfxEmissionType::tfxIcosphere:
					spawn_quantity = (float)tfxIcospherePoints[tfxMin((tfxU32)x, 5)].current_size;
					break;
				}
			}
			step_size = 1.f / spawn_quantity;
		} else {
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
				//amount_remainder += step_size;
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
		else if(!(property_flags & tfxEmitterPropertyFlags_single)) {
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
		else {
			work_entry.amount_to_spawn = (tfxU32)pm.emitters.spawn_quantity[work_entry.emitter_index];
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
		const tfxEmitterStateFlags state_flags = pm.emitters.state_flags[emitter_index];
		const tfxEmitterPropertyFlags property_flags = pm.emitters.property_flags[emitter_index];
		const float emitter_intensity = pm.emitters.intensity[emitter_index];
		const float first_red_value = GetGraphFirstValue(&library->emitter_attributes[emitter_attributes].overtime.red);
		const float first_green_value = GetGraphFirstValue(&library->emitter_attributes[emitter_attributes].overtime.green);
		const float first_blue_value = GetGraphFirstValue(&library->emitter_attributes[emitter_attributes].overtime.blue);
		const float first_alpha_value = GetGraphFirstValue(&library->emitter_attributes[emitter_attributes].overtime.blendfactor);
		const float first_intensity_value = GetGraphFirstValue(&library->emitter_attributes[emitter_attributes].overtime.intensity);

		assert(random.seeds[0] > 0);

		for (int i = 0; i != entry->amount_to_spawn; ++i) {
			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			entry->particle_data->uid[index] = pm.unique_particle_id++;
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
			if (property_flags & tfxEmitterPropertyFlags_wrap_single_sprite && pm.flags & tfxEffectManagerFlags_recording_sprites) {
				max_age = pm.animation_length_in_time;
			}
			//Remove this unless I change my mind about it:
			//Only makes sense for single particles
			//else if(property_flags & tfxEmitterPropertyFlags_life_proportional_to_animation && pm.flags & tfxEffectManagerFlags_recording_sprites) {
			//	float life_proportion = life / (life + life_variation);
			//	float life_variation_proportion = life_variation / (life + life_variation);
			//	max_age = life_proportion * pm.animation_length_in_time + random.Range(life_variation_proportion * pm.animation_length_in_time);
			//}
			else {
				max_age = life + random.Range(life_variation);
			}
			single_loop_count = 0;

			float alpha = 255.f * first_alpha_value;
			float intensity = first_intensity_value * emitter_intensity;
			//intensity = 0.f;
			if (state_flags & tfxEmitterStateFlags_random_color) {
				float age = random.Range(max_age);
				color = tfxRGBA8(	255.f * lookup_overtime_callback(library->emitter_attributes[emitter_attributes].overtime.red, age, max_age),
									255.f * lookup_overtime_callback(library->emitter_attributes[emitter_attributes].overtime.green, age, max_age),
									255.f * lookup_overtime_callback(library->emitter_attributes[emitter_attributes].overtime.blue, age, max_age), alpha);
			}
			else {
				color = tfxRGBA8(255.f * first_red_value, 255.f * first_green_value, 255.f * first_blue_value, alpha);
			}

			entry->highest_particle_age = std::fmaxf(highest_particle_age, (max_age * loop_count) + pm.frame_length + 1);

			if (entry->sub_effects->current_size > 0) {
				particle_index = GetPMParticleIndexSlot(&pm, MakeParticleID(particles_index, index));
				flags |= tfxParticleFlags_has_sub_effects;
				for (auto &sub : *entry->sub_effects) {
					if (!FreePMEffectCapacity(&pm))
						break;
					assert(entry->depth < tfxMAXDEPTH - 1);
					tfxU32 added_index = AddEffectToParticleManager(&pm, &sub, pm.current_ebuff, entry->depth + 1, true, 0.f);
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
			//data.image = GetEffectProperties(this)->image;
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

		for (int i = 0; i != entry->amount_to_spawn; ++i) {

			tfxU32 index = GetCircularIndex(&pm.particle_array_buffers[particles_index], entry->spawn_start_index + i);
			float &base_size_x = entry->particle_data->base_size_x[index];
			float &base_size_y = entry->particle_data->base_size_y[index];

			//----Size
			if (!(property_flags & tfxEmitterPropertyFlags_base_uniform_size)) {
				float random_size_x = random.Range(size_variation.x);
				float random_size_y = random.Range(size_variation.y);
				base_size_y = (random_size_y + size.y);
				base_size_x = (random_size_x + size.x);
			}
			else {
				float random_size_x = random.Range(size_variation.x);
				float random_size_y = random_size_x;
				base_size_y = (random_size_y + size.y);
				base_size_x = (random_size_x + size.x);
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
										grid_coords.z = grid_points.z - 1;
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
						if (!(property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
							if ((grid_coords.z > 0 && grid_coords.z < grid_points.z - 1) || property_flags & tfxEmitterPropertyFlags_area_open_ends) {
								grid_coords.x -= grid_coords.y == 0 || grid_coords.y == grid_points.y - 1 ? 1.f : grid_points.x - 1;
								if (grid_coords.x < 0.f) {
									grid_coords.y--;
									grid_coords.x = grid_points.x - 1;
									if (grid_coords.y < 0.f) {
										grid_coords.z--;
										grid_coords.y = grid_points.y - 1;
										if (grid_coords.z < 0.f)
											grid_coords.z = grid_points.z - 1;
									}
								}
							}
							else {
								grid_coords.x--;
								if (grid_coords.x < 0.f) {
									grid_coords.y--;
									grid_coords.x = grid_points.x - 1;
									if (grid_coords.y < 0.f) {
										grid_coords.z--;
										grid_coords.y = grid_points.y - 1;
										if (grid_coords.z < 0.f)
											grid_coords.z = grid_points.z - 1;
									}
								}
							}
						}

						local_position_x = position.x + (grid_coords.x * grid_segment_size.x);
						local_position_y = position.y + (grid_coords.y * grid_segment_size.y);
						local_position_z = position.z + (grid_coords.z * grid_segment_size.z);

						if (property_flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
							if ((grid_coords.z > 0 && grid_coords.z < grid_points.z - 1) || property_flags & tfxEmitterPropertyFlags_area_open_ends) {
								grid_coords.x += grid_coords.y == 0 || grid_coords.y == grid_points.y - 1 ? 1.f : grid_points.x - 1;
								if (grid_coords.x >= grid_points.x) {
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
							else {
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
				local_position_y = random.Range(emitter_size.y);
				local_position_z = -std::sinf(th) * half_emitter_size.z + half_emitter_size.z;
			}
			else {
				local_position_x = random.Range(0.f, emitter_size.x);
				local_position_y = random.Range(0.f, emitter_size.y);
				local_position_z = random.Range(0.f, emitter_size.z);

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
		const float first_velocity_value = GetGraphFirstValue(&library->emitter_attributes[emitter_attributes].overtime.velocity);
		const float first_weight_value = GetGraphFirstValue(&library->emitter_attributes[emitter_attributes].overtime.weight);
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

			float micro_time = pm.update_time * (1.f - tween);

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
				direction = GetEmissionDirection2d(pm, library, random, property_index, emitter_index, tfxVec2(local_position_x, local_position_y), sprite_transform_position, emitter_size) + GetGraphFirstValue(&library->emitter_attributes[emitter_attributes].overtime.direction);
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
				entry->particle_data->depth_index[index] = PushPMDepthIndex(&pm, layer, depth_index);
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
		const float first_velocity_value = GetGraphFirstValue(&library->emitter_attributes[emitter_attributes].overtime.velocity);
		const float first_weight_value = GetGraphFirstValue(&library->emitter_attributes[emitter_attributes].overtime.weight);
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
			float micro_time = pm.update_time * (1.f - tween + 0.001f);
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
				entry->particle_data->depth_index[index] = PushPMDepthIndex(&pm, layer, depth_index);
			}
			else if (pm.flags & tfxEffectManagerFlags_ordered_by_age) {
				tfxDepthIndex depth_index;
				depth_index.particle_id = MakeParticleID(particles_index, index);
				depth_index.depth = 0.f;
				entry->particle_data->depth_index[index] = PushPMDepthIndex(&pm, layer, depth_index);
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
		const tfxU32 property_index = pm.emitters.properties_index[emitter_index];
		const tfxU32 property_flags = pm.emitters.property_flags[emitter_index];
		const tfxWideInt single_shot_limit = tfxWideSetSinglei(work_entry->properties->single_shot_limit[property_index]);
		const tfxU32 layer = work_entry->properties->layer[property_index];

		const tfxWideInt remove_flag = tfxWideSetSinglei(tfxParticleFlags_remove);
		const tfxWideInt capture_after_transform = tfxWideSetSinglei(tfxParticleFlags_capture_after_transform);
		const tfxWideInt remove = tfxWideSetSinglei(pm.emitters.state_flags[emitter_index] & tfxParticleFlags_remove);
		const tfxWideInt single = tfxWideGreateri(tfxWideSetSinglei(property_flags & tfxEmitterPropertyFlags_single), tfxWideSetZeroi());
		const tfxWideInt not_single = tfxWideXOri(single, tfxWideSetSinglei(-1));
		tfxWideInt state_flags_no_spawning = tfxWideGreateri(tfxWideOri(tfxWideSetSinglei(pm.emitters.state_flags[emitter_index] & tfxEmitterStateFlags_stop_spawning), tfxWideSetSinglei(work_entry->pm->flags & tfxEffectManagerFlags_disable_spawning)), tfxWideSetZeroi());
		if (property_flags & tfxEmitterPropertyFlags_wrap_single_sprite && pm.flags & tfxEffectManagerFlags_recording_sprites) {
			state_flags_no_spawning = tfxWideGreateri(tfxWideSetSinglei(property_flags & tfxEmitterPropertyFlags_wrap_single_sprite), tfxWideSetZeroi());
		}
		const tfxWideInt xor_state_flags_no_spawning = tfxWideXOri(state_flags_no_spawning, tfxWideSetSinglei(-1));

		tfxParticleSoA &bank = pm.particle_arrays[particles_index];

		for (int i = 0; i != work_entry->wide_end_index; i += tfxDataWidth) {
			tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i) / tfxDataWidth * tfxDataWidth;

			const tfxWideFloat max_age = tfxWideLoad(&bank.max_age[index]);
			tfxWideFloat age = tfxWideLoad(&bank.age[index]);
			tfxWideInt single_loop_count = tfxWideLoadi((tfxWideInt*)&bank.single_loop_count[index]);
			tfxWideInt flags = tfxWideLoadi((tfxWideInt*)&bank.flags[index]);
			age = tfxWideAdd(age, pm.frame_length_wide);

			tfxWideInt expired = tfxWideCasti(tfxWideGreaterEqual(age, max_age));
			single_loop_count = tfxWideAddi(single_loop_count, tfxWideAndi(tfxWideSetSinglei(1), expired));
			tfxWideInt loop_limit = tfxWideEqualsi(single_loop_count, single_shot_limit);
			tfxWideInt loop_age = tfxWideXOri(tfxWideAndi(tfxWideAndi(single, expired), xor_state_flags_no_spawning), tfxWideSetSinglei(-1));
			age = tfxWideAnd(age, tfxWideCast(loop_age));
			flags = tfxWideOri(flags, tfxWideAndi(remove_flag, tfxWideGreateri(remove, tfxWideSetSinglei(0))));
			flags = tfxWideOri(flags, tfxWideAndi(remove_flag, tfxWideAndi(not_single, expired)));
			flags = tfxWideOri(flags, tfxWideAndi(remove_flag, tfxWideAndi(tfxWideOri(tfxWideAndi(single, loop_limit), state_flags_no_spawning), expired)));
			flags = tfxWideOri(flags, tfxWideAndi(capture_after_transform, expired));

			tfxWideStore(&bank.age[index], age);
			tfxWideStorei((tfxWideInt*)&bank.flags[index], flags);
			tfxWideStorei((tfxWideInt*)&bank.single_loop_count[index], single_loop_count);
		}

		tfxU32 offset = 0;
		for (int i = work_entry->start_index; i >= 0; --i) {
			const tfxU32 index = GetCircularIndex(&work_entry->pm->particle_array_buffers[particles_index], i);
			tfxParticleFlags &flags = bank.flags[index];
			if (flags & tfxParticleFlags_remove) {
				offset++;
				if (flags & tfxParticleFlags_has_sub_effects) {
					FreePMParticleIndex(&pm, bank.particle_index[index]);
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
				if (pm.flags & tfxEffectManagerFlags_recording_sprites) {
					bank.uid[next_index] = bank.uid[index];
				}
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

	void ControlParticles(tfxParticleManager &pm, tfxU32 emitter_index, tfxControlWorkEntry &work_entry) {
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
				if (pm.flags & tfxEffectManagerFlags_3d_effects) {
					tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticlePosition3d);
				}
				else {
					tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticlePosition2d);
				}
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleSize);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleColor);
				tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleImageFrame);
				if (pm.flags & tfxEffectManagerFlags_recording_sprites && pm.flags & tfxEffectManagerFlags_using_uids) {
					tfxAddWorkQueueEntry(&pm.work_queue, &work_entry, ControlParticleUID);
				}
			}
			else {
				if (pm.flags & tfxEffectManagerFlags_3d_effects) {
					ControlParticlePosition3d(&pm.work_queue, &work_entry);
				}
				else {
					ControlParticlePosition2d(&pm.work_queue, &work_entry);
				}
				ControlParticleSize(&pm.work_queue, &work_entry);
				ControlParticleColor(&pm.work_queue, &work_entry);
				ControlParticleImageFrame(&pm.work_queue, &work_entry);
				if (pm.flags & tfxEffectManagerFlags_recording_sprites && pm.flags & tfxEffectManagerFlags_using_uids) {
					ControlParticleUID(&pm.work_queue, &work_entry);
				}
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
	tfxDataTypesDictionary tfxDataTypes;
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

	void InitSpriteData3dSoACompression(tfxSoABuffer *buffer, tfxSpriteDataSoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteDataSoA, property_indexes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteDataSoA, captured_index));
		AddStructArray(buffer, sizeof(tfxUniqueSpriteID), offsetof(tfxSpriteDataSoA, uid));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteDataSoA, lerp_offset));
		AddStructArray(buffer, sizeof(tfxSpriteTransform3d), offsetof(tfxSpriteDataSoA, transform_3d));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteDataSoA, alignment));
		AddStructArray(buffer, sizeof(tfxRGBA8), offsetof(tfxSpriteDataSoA, color));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteDataSoA, stretch));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteDataSoA, intensity));
		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	void InitSpriteData3dSoA(tfxSoABuffer *buffer, tfxSpriteDataSoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteDataSoA, property_indexes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteDataSoA, captured_index));
		AddStructArray(buffer, sizeof(tfxUniqueSpriteID), offsetof(tfxSpriteDataSoA, uid));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteDataSoA, lerp_offset));
		AddStructArray(buffer, sizeof(tfxSpriteTransform3d), offsetof(tfxSpriteDataSoA, transform_3d));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteDataSoA, alignment));
		AddStructArray(buffer, sizeof(tfxRGBA8), offsetof(tfxSpriteDataSoA, color));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteDataSoA, stretch));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteDataSoA, intensity));
		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	void InitSpriteData2dSoACompression(tfxSoABuffer *buffer, tfxSpriteDataSoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteDataSoA, property_indexes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteDataSoA, captured_index));
		AddStructArray(buffer, sizeof(tfxUniqueSpriteID), offsetof(tfxSpriteDataSoA, uid));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteDataSoA, lerp_offset));
		AddStructArray(buffer, sizeof(tfxSpriteTransform2d), offsetof(tfxSpriteDataSoA, transform_2d));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteDataSoA, alignment));
		AddStructArray(buffer, sizeof(tfxRGBA8), offsetof(tfxSpriteDataSoA, color));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteDataSoA, stretch));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteDataSoA, intensity));
		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	void InitSpriteBufferSoA(tfxSoABuffer *buffer, tfxSpriteSoA *soa, tfxU32 reserve_amount, tfxSpriteBufferMode mode, bool use_uid) {
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteSoA, property_indexes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteSoA, captured_index));
		if (use_uid)
			AddStructArray(buffer, sizeof(tfxUniqueSpriteID), offsetof(tfxSpriteSoA, uid));
		if (mode == tfxSpriteBufferMode_2d) {
			AddStructArray(buffer, sizeof(tfxSpriteTransform2d), offsetof(tfxSpriteSoA, transform_2d));
		}
		else if (mode == tfxSpriteBufferMode_3d) {
			AddStructArray(buffer, sizeof(tfxSpriteTransform3d), offsetof(tfxSpriteSoA, transform_3d));
		}
		else {
			AddStructArray(buffer, sizeof(tfxSpriteTransform2d), offsetof(tfxSpriteSoA, transform_2d));
			AddStructArray(buffer, sizeof(tfxSpriteTransform3d), offsetof(tfxSpriteSoA, transform_3d));
		}
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteSoA, alignment));
		AddStructArray(buffer, sizeof(tfxRGBA8), offsetof(tfxSpriteSoA, color));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteSoA, stretch));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteSoA, intensity));
		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	void InitSpriteData2dSoA(tfxSoABuffer *buffer, tfxSpriteDataSoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteDataSoA, property_indexes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteDataSoA, captured_index));
		AddStructArray(buffer, sizeof(tfxUniqueSpriteID), offsetof(tfxSpriteDataSoA, uid));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteDataSoA, lerp_offset));
		AddStructArray(buffer, sizeof(tfxSpriteTransform2d), offsetof(tfxSpriteDataSoA, transform_2d));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxSpriteDataSoA, alignment));
		AddStructArray(buffer, sizeof(tfxRGBA8), offsetof(tfxSpriteDataSoA, color));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteDataSoA, stretch));
		AddStructArray(buffer, sizeof(float), offsetof(tfxSpriteDataSoA, intensity));
		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	void InitParticleSoA(tfxSoABuffer *buffer, tfxParticleSoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxParticleSoA, uid));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxParticleSoA, parent_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxParticleSoA, sprite_index));
		AddStructArray(buffer, sizeof(tfxParticleID), offsetof(tfxParticleSoA, particle_index));
		AddStructArray(buffer, sizeof(tfxParticleFlags), offsetof(tfxParticleSoA, flags));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, age));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, max_age));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, position_x));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, position_y));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, position_z));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, captured_position_x));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, captured_position_y));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, captured_position_z));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, local_rotations_x));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, local_rotations_y));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, local_rotations_z));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxParticleSoA, velocity_normal));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxParticleSoA, depth_index));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, base_weight));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, base_velocity));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, base_spin));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, noise_offset));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, noise_resolution));
		AddStructArray(buffer, sizeof(tfxRGBA8), offsetof(tfxParticleSoA, color));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, image_frame));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, base_size_x));
		AddStructArray(buffer, sizeof(float), offsetof(tfxParticleSoA, base_size_y));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxParticleSoA, single_loop_count));
		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	void InitEffectSoA(tfxSoABuffer *buffer, tfxEffectSoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, frame));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, age));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, highest_particle_age));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, timeout_counter));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, timeout));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, handle));
		AddStructArray(buffer, sizeof(tfxEmitterPropertyFlags), offsetof(tfxEffectSoA, property_flags));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, loop_length));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, translation));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, local_position));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, world_position));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, captured_position));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, local_rotations));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, world_rotations));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, scale));
		//Todo: save space and use a quaternion here?
		AddStructArray(buffer, sizeof(tfxMatrix4), offsetof(tfxEffectSoA, matrix));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEffectSoA, global_attributes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEffectSoA, transform_attributes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEffectSoA, parent_particle_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEffectSoA, properties_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEffectSoA, info_index));
		AddStructArray(buffer, sizeof(void*), offsetof(tfxEffectSoA, library));
		AddStructArray(buffer, sizeof(tfxParentSpawnControls), offsetof(tfxEffectSoA, spawn_controls));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEffectSoA, emitter_size));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, stretch));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, overal_scale));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEffectSoA, noise_base_offset));
		AddStructArray(buffer, sizeof(tfxEffectStateFlags), offsetof(tfxEffectSoA, state_flags));
		AddStructArray(buffer, sizeof(void*), offsetof(tfxEffectSoA, user_data));
		AddStructArray(buffer, sizeof(void*), offsetof(tfxEffectSoA, update_callback));
		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	void InitEmitterSoA(tfxSoABuffer *buffer, tfxEmitterSoA *soa, tfxU32 reserve_amount) {
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, frame));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, age));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, highest_particle_age));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, delay_spawning));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, timeout_counter));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, timeout));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, handle));
		AddStructArray(buffer, sizeof(tfxEmitterPropertyFlags), offsetof(tfxEmitterSoA, property_flags));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, loop_length));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, translation));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, local_position));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, world_position));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, captured_position));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, local_rotations));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, world_rotations));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, scale));
		//Todo: save space and use a quaternion here?
		AddStructArray(buffer, sizeof(tfxMatrix4), offsetof(tfxEmitterSoA, matrix));
		AddStructArray(buffer, sizeof(tfxVec2), offsetof(tfxEmitterSoA, image_handle));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, amount_remainder));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, spawn_quantity));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, qty_step_size));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, emitter_attributes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, transform_attributes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, overtime_attributes));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, parent_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, sprites_count));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, sprites_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, seed_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, properties_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, info_index));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, hierarchy_depth));
		AddStructArray(buffer, sizeof(tfxKey), offsetof(tfxEmitterSoA, path_hash));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, life));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, life_variation));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, arc_size));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, arc_offset));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, weight));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, weight_variation));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, velocity));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, velocity_variation));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, spin));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, spin_variation));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, splatter));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, noise_offset_variation));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, noise_offset));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, noise_resolution));
		AddStructArray(buffer, sizeof(tfxVec2), offsetof(tfxEmitterSoA, size));
		AddStructArray(buffer, sizeof(tfxVec2), offsetof(tfxEmitterSoA, size_variation));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, grid_segment_size));
		AddStructArray(buffer, sizeof(tfxU32), offsetof(tfxEmitterSoA, particles_index));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, overal_scale));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, velocity_adjuster));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, intensity));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, image_frame_rate));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, stretch));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, end_frame));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, grid_coords));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, grid_direction));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, emitter_size));
		AddStructArray(buffer, sizeof(float), offsetof(tfxEmitterSoA, emission_alternator));
		AddStructArray(buffer, sizeof(tfxEmitterStateFlags), offsetof(tfxEmitterSoA, state_flags));
		AddStructArray(buffer, sizeof(tfxVec2), offsetof(tfxEmitterSoA, image_size));
		AddStructArray(buffer, sizeof(tfxVec3), offsetof(tfxEmitterSoA, angle_offsets));
		FinishSoABufferSetup(buffer, soa, reserve_amount);
	}

	void InitEmitterProperites(tfxEmitterPropertiesSoA &properties, tfxU32 i) {
		properties.angle_offsets[i] = { 0.f, 0.f, tfx360Radians };
		properties.image[i] = nullptr;
		properties.image_handle[i] = tfxVec2();
		properties.spawn_amount[i] = 1;
		properties.single_shot_limit[i] = 0;
		properties.emission_type[i] = tfxEmissionType::tfxPoint;
		properties.billboard_option[i] = tfxBillboarding_align_to_camera;
		properties.vector_align_type[i] = tfxVectorAlignType_motion;
		properties.emission_direction[i] = tfxEmissionDirection::tfxOutwards;
		properties.grid_points[i] = { 10.f, 10.f, 10.f };
		properties.emitter_handle[i] = { 0.f, 0.f, 0.f };
		properties.end_behaviour[i] = tfxLineTraversalEndBehaviour::tfxLoop;
		properties.loop_length[i] = 0.f;
		properties.layer[i] = 0;
		properties.image_hash[i] = 1;
		properties.start_frame[i] = 0;
		properties.end_frame[i] = 0;
		properties.frame_rate[i] = 30.f;
		properties.angle_settings[i] = tfxAngleSettingFlags_random_roll | tfxAngleSettingFlags_specify_pitch | tfxAngleSettingFlags_specify_yaw;
		properties.delay_spawning[i] = 0.f;
		properties.noise_base_offset_range[i] = 1000.f;
		properties.animation_property_index[i] = tfxINVALID;
	}

	//Use with care, no checks for out of bounds
	void CopyEmitterProperites(tfxEmitterPropertiesSoA &from_properties, tfxU32 from_i, tfxEmitterPropertiesSoA &to_properties, tfxU32 to_i) {
		to_properties.angle_offsets[to_i] = from_properties.angle_offsets[from_i];
		to_properties.image[to_i] = from_properties.image[from_i];
		to_properties.image_handle[to_i] = from_properties.image_handle[from_i];
		to_properties.spawn_amount[to_i] = from_properties.spawn_amount[from_i];
		to_properties.single_shot_limit[to_i] = from_properties.single_shot_limit[from_i];
		to_properties.emission_type[to_i] = from_properties.emission_type[from_i];
		to_properties.billboard_option[to_i] = from_properties.billboard_option[from_i];
		to_properties.vector_align_type[to_i] = from_properties.vector_align_type[from_i];
		to_properties.emission_direction[to_i] = from_properties.emission_direction[from_i];
		to_properties.grid_points[to_i] = from_properties.grid_points[from_i];
		to_properties.emitter_handle[to_i] = from_properties.emitter_handle[from_i];
		to_properties.end_behaviour[to_i] = from_properties.end_behaviour[from_i];
		to_properties.loop_length[to_i] = from_properties.loop_length[from_i];
		to_properties.layer[to_i] = from_properties.layer[from_i];
		to_properties.image_hash[to_i] = from_properties.image_hash[from_i];
		to_properties.start_frame[to_i] = from_properties.start_frame[from_i];
		to_properties.end_frame[to_i] = from_properties.end_frame[from_i];
		to_properties.frame_rate[to_i] = from_properties.frame_rate[from_i];
		to_properties.angle_settings[to_i] = from_properties.angle_settings[from_i];
		to_properties.delay_spawning[to_i] = from_properties.delay_spawning[from_i];
		to_properties.noise_base_offset_range[to_i] = from_properties.noise_base_offset_range[from_i];
		to_properties.animation_property_index[to_i] = from_properties.animation_property_index[from_i];
	}

	void FreeSpriteData(tfxSpriteData &sprite_data) {
		if (sprite_data.compressed_sprites_buffer.data == sprite_data.real_time_sprites_buffer.data) {
			FreeSoABuffer(&sprite_data.real_time_sprites_buffer);
			sprite_data.normal.frame_meta.free_all();
			sprite_data.compressed_sprites_buffer = tfxSoABuffer();
		}
		else {
			FreeSoABuffer(&sprite_data.compressed_sprites_buffer);
			FreeSoABuffer(&sprite_data.real_time_sprites_buffer);
			sprite_data.normal.frame_meta.free_all();
			sprite_data.compressed.frame_meta.free_all();
		}
	}

	void InitParticleManagerFor3d(tfxParticleManager *pm, tfxLibrary *library, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfxParticleManagerModes mode, bool double_buffered_sprites, bool dynamic_sprite_allocation, tfxU32 mt_batch_size) {
		assert(pm->flags == 0);		//You must use a particle manager that has not been initialised already. You can call reconfigure if you want to re-initialise a particle manager
		pm->max_effects = effects_limit;
		pm->mt_batch_size = mt_batch_size;
		tfxInitialiseWorkQueue(&pm->work_queue);
		pm->library = library;

		if (pm->particle_array_allocator.arenas.current_size == 0) {
			//todo need to be able to adjust the arena size
			pm->particle_array_allocator = CreateArenaManager(tfxMegabyte(2), 8);
			pm->particle_arrays = tfxBucketArray<tfxParticleSoA>(&pm->particle_array_allocator, 32);
		}

		pm->flags = 0;

		if (mode == tfxParticleManagerMode_unordered)
			pm->flags = tfxEffectManagerFlags_unordered;
		else if (mode == tfxParticleManagerMode_ordered_by_age)
			pm->flags = tfxEffectManagerFlags_ordered_by_age;
		else if (mode == tfxParticleManagerMode_ordered_by_depth)
			pm->flags = tfxEffectManagerFlags_order_by_depth;
		else if (mode == tfxParticleManagerMode_ordered_by_depth_guaranteed)
			pm->flags = tfxEffectManagerFlags_order_by_depth | tfxEffectManagerFlags_guarantee_order;

		pm->flags |= tfxEffectManagerFlags_3d_effects;
		if (double_buffered_sprites)
			pm->flags |= tfxEffectManagerFlags_double_buffer_sprites;

		for (tfxEachLayer) {
			pm->max_cpu_particles_per_layer[layer] = layer_max_values[layer];

#ifdef tfxTRACK_MEMORY
			memcpy(sprites3d[layer].name, "ParticleManager::sprites3d\0", 27);
#endif

			InitSpriteBufferSoA(&pm->sprite_buffer[0][layer], &pm->sprites[0][layer], tfxMax((layer_max_values[layer] / tfxDataWidth + 1) * tfxDataWidth, 8), tfxSpriteBufferMode_3d);
			if (pm->flags & tfxEffectManagerFlags_double_buffer_sprites) {
				InitSpriteBufferSoA(&pm->sprite_buffer[1][layer], &pm->sprites[1][layer], tfxMax((layer_max_values[layer] / tfxDataWidth + 1) * tfxDataWidth, 8), tfxSpriteBufferMode_3d);
			}

		}

		if (pm->flags & tfxEffectManagerFlags_ordered_by_age || pm->flags & tfxEffectManagerFlags_order_by_depth) {
			FreeParticleBanks(pm);
			for (tfxEachLayer) {
				tfxParticleSoA lists;
				tfxU32 index = pm->particle_arrays.locked_push_back(lists);
				tfxSoABuffer buffer;
				pm->particle_array_buffers.push_back(buffer);
				assert(index == pm->particle_array_buffers.current_size - 1);
				InitParticleSoA(&pm->particle_array_buffers[index], &pm->particle_arrays.back(), tfxMax(pm->max_cpu_particles_per_layer[layer], 8));
				pm->particle_array_buffers[index].user_data = &pm->particle_arrays.back();
				tfxResizeParticleSoACallback(&pm->particle_array_buffers[index], 0);
				pm->depth_indexes[layer][0].reserve(pm->max_cpu_particles_per_layer[layer]);
				pm->depth_indexes[layer][1].reserve(pm->max_cpu_particles_per_layer[layer]);
			}
		}

		pm->flags |= dynamic_sprite_allocation ? tfxEffectManagerFlags_dynamic_sprite_allocation : 0;

		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			pm->effects_in_use[depth][0].reserve(pm->max_effects);
			pm->effects_in_use[depth][1].reserve(pm->max_effects);

			pm->emitters_in_use[depth][0].reserve(pm->max_effects);
			pm->emitters_in_use[depth][1].reserve(pm->max_effects);
		}
		pm->particle_indexes.reserve(1000);	//todo: Handle this better.
		pm->free_effects.reserve(pm->max_effects);
		InitEffectSoA(&pm->effect_buffers, &pm->effects, pm->max_effects);
		InitEmitterSoA(&pm->emitter_buffers, &pm->emitters, pm->max_effects);
	}

	void InitParticleManagerFor2d(tfxParticleManager *pm, tfxLibrary *library, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfxParticleManagerModes mode, bool double_buffered_sprites, bool dynamic_sprite_allocation, tfxU32 mt_batch_size) {
		assert(pm->flags == 0);		//You must use a particle manager that has not been initialised already. You can call reconfigure if you want to re-initialise a particle manager
		assert(mode == tfxParticleManagerMode_unordered || mode == tfxParticleManagerMode_ordered_by_age);	//Only these 2 modes are available for 2d effects
		pm->max_effects = effects_limit;
		pm->mt_batch_size = mt_batch_size;
		tfxInitialiseWorkQueue(&pm->work_queue);
		pm->library = library;

		if (pm->particle_array_allocator.arenas.current_size == 0) {
			//todo need to be able to adjust the arena size
			pm->particle_array_allocator = CreateArenaManager(tfxMegabyte(2), 8);
			pm->particle_arrays = tfxBucketArray<tfxParticleSoA>(&pm->particle_array_allocator, 32);
		}

		pm->flags = 0;

		if (mode == tfxParticleManagerMode_unordered)
			pm->flags = tfxEffectManagerFlags_unordered;
		else if (mode == tfxParticleManagerMode_ordered_by_age)
			pm->flags = tfxEffectManagerFlags_ordered_by_age;
		else if (mode == tfxParticleManagerMode_ordered_by_depth)
			pm->flags = tfxEffectManagerFlags_order_by_depth;
		else if (mode == tfxParticleManagerMode_ordered_by_depth_guaranteed)
			pm->flags = tfxEffectManagerFlags_order_by_depth | tfxEffectManagerFlags_guarantee_order;

		for (tfxEachLayer) {
			pm->max_cpu_particles_per_layer[layer] = layer_max_values[layer];

#ifdef tfxTRACK_MEMORY
			memcpy(sprites2d[layer].name, "ParticleManager::sprites2d\0", 27);
#endif

			InitSpriteBufferSoA(&pm->sprite_buffer[0][layer], &pm->sprites[0][layer], tfxMax((layer_max_values[layer] / tfxDataWidth + 1) * tfxDataWidth, 8), tfxSpriteBufferMode_2d);
			if (pm->flags & tfxEffectManagerFlags_double_buffer_sprites) {
				InitSpriteBufferSoA(&pm->sprite_buffer[1][layer], &pm->sprites[1][layer], tfxMax((layer_max_values[layer] / tfxDataWidth + 1) * tfxDataWidth, 8), tfxSpriteBufferMode_2d);
			}
		}

		if (!(pm->flags & tfxEffectManagerFlags_unordered)) {
			for (tfxEachLayer) {
				tfxParticleSoA lists;
				tfxU32 index = pm->particle_arrays.locked_push_back(lists);
				tfxSoABuffer buffer;
				pm->particle_array_buffers.push_back(buffer);
				InitParticleSoA(&pm->particle_array_buffers[index], &pm->particle_arrays.back(), layer_max_values[layer]);
				assert(index == pm->particle_array_buffers.current_size - 1);
				pm->depth_indexes[layer][0].reserve(pm->max_cpu_particles_per_layer[layer]);
				pm->depth_indexes[layer][1].reserve(pm->max_cpu_particles_per_layer[layer]);
			}
		}

		pm->flags |= dynamic_sprite_allocation ? tfxEffectManagerFlags_dynamic_sprite_allocation : 0;

		for (int depth = 0; depth != tfxMAXDEPTH; ++depth) {
			pm->effects_in_use[depth][0].reserve(pm->max_effects);
			pm->effects_in_use[depth][1].reserve(pm->max_effects);

			pm->emitters_in_use[depth][0].reserve(pm->max_effects);
			pm->emitters_in_use[depth][1].reserve(pm->max_effects);
		}

		pm->free_effects.reserve(pm->max_effects);
		pm->particle_indexes.reserve(1000);	//todo: Handle this better.
		InitEffectSoA(&pm->effect_buffers, &pm->effects, pm->max_effects);
		InitEmitterSoA(&pm->emitter_buffers, &pm->emitters, pm->max_effects);

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

	void SetAnimationPosition(tfxAnimationManager *animation_manager, tfxAnimationID effect_index, float position[3]) {
		animation_manager->instances[effect_index].position.x = position[0];
		animation_manager->instances[effect_index].position.y = position[1];
		animation_manager->instances[effect_index].position.z = position[2];
	}

	void SetAnimationPosition(tfxAnimationManager *animation_manager, tfxAnimationID effect_index, float x, float y) {
		animation_manager->instances[effect_index].position.x = x;
		animation_manager->instances[effect_index].position.y = y;
	}

	void SetAnimationScale(tfxAnimationManager *animation_manager, tfxAnimationID effect_index, float scale) {
		animation_manager->instances[effect_index].scale = scale;
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