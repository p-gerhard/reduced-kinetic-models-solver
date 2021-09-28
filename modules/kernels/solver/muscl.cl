#ifndef MUSCL_CL
#define MUSCL_CL

static void muscl_reconstruct(const __global float *wn, const long id_cell,
							  const long id_neighbour, const float vn[DIM],
							  const float sL[DIM][M], const float sR[DIM][M],
							  float wL[M], float wR[M])
{
#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		long imemR = id_neighbour + iw * NGRID;
		long imemL = id_cell + iw * NGRID;

#ifdef IS_2D
		wL[iw] = wn[imemL];
		wL[iw] += vn[0] * 0.5f * sL[0][iw];
		wL[iw] += vn[1] * 0.5f * sL[1][iw];

		wR[iw] = wn[imemR];
		wR[iw] -= vn[0] * 0.5f * sR[0][iw];
		wR[iw] -= vn[1] * 0.5f * sR[1][iw];
#else
		wL[iw] = wn[imemL];
		wL[iw] += vn[0] * 0.5f * sL[0][iw];
		wL[iw] += vn[1] * 0.5f * sL[1][iw];
		wL[iw] += vn[2] * 0.5f * sL[2][iw];

		wR[iw] = wn[imemR];
		wR[iw] -= vn[0] * 0.5f * sR[0][iw];
		wR[iw] -= vn[1] * 0.5f * sR[1][iw];
		wR[iw] -= vn[2] * 0.5f * sR[2][iw];
#endif
	}
}

static void muscl_reconstruct_boundary(const __global float *wn,
									   const long id_cell, const float vn[DIM],
									   const float sL[DIM][M], float wL[M])

{
#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		long imemL = id_cell + iw * NGRID;

#ifdef IS_2D
		wL[iw] = wn[imemL];
		wL[iw] += vn[0] * 0.5f * sL[0][iw];
		wL[iw] += vn[1] * 0.5f * sL[1][iw];
#else
		wL[iw] = wn[imemL];
		wL[iw] += vn[0] * 0.5f * sL[0][iw];
		wL[iw] += vn[1] * 0.5f * sL[1][iw];
		wL[iw] += vn[2] * 0.5f * sL[2][iw];
#endif
	}
}

#ifdef IS_MOMENT_MODEL
static void muscl_m1_reconstruct(const __global float *wn, const long id_cell,
								 const long id_neighbour, const float vn[DIM],
								 const float sL[DIM][M], const float sR[DIM][M],
								 float wL[M], float wR[M])

{
	float r;
	long imemL = id_cell + 0 * NGRID;
	long imemR = id_neighbour + 0 * NGRID;
#ifdef IS_2D
	wL[0] = wn[imemL];
	wL[0] += vn[0] * 0.5f * sL[0][0];
	wL[0] += vn[1] * 0.5f * sL[1][0];

	wR[0] = wn[imemR];
	wR[0] -= vn[0] * 0.5f * sR[0][0];
	wR[0] -= vn[1] * 0.5f * sR[1][0];
#else
	wL[0] = wn[imemL];
	wL[0] += vn[0] * 0.5f * sL[0][0];
	wL[0] += vn[1] * 0.5f * sL[1][0];
	wL[0] += vn[2] * 0.5f * sL[2][0];

	wR[0] = wn[imemR];
	wR[0] -= vn[0] * 0.5f * sR[0][0];
	wR[0] -= vn[1] * 0.5f * sR[1][0];
	wR[0] -= vn[2] * 0.5f * sR[2][0];
#endif

#pragma unroll
	float rhoL = wn[id_cell + 0 * NGRID];
	float rhoR = wn[id_neighbour + 0 * NGRID];
	for (int k = 1; k < M; k++) {
		/* Reconstruction using current cell value */
		r = get_r(rhoL, wn[id_cell + k * NGRID]);

		wL[k] = wL[0] * (r + vn[0] * 0.5f * sL[0][k] + vn[1] * 0.5f * sL[1][k] +
						 vn[2] * 0.5f * sL[2][k]);
		/* Reconstruction using neighbour cell value */
		r = get_r(rhoR, wn[id_neighbour + k * NGRID]);

		wR[k] = wR[0] * (r - vn[0] * 0.5f * sR[0][k] - vn[1] * 0.5f * sR[1][k] -
						 vn[2] * 0.5f * sR[2][k]);
	}
}

static void muscl_m1_reconstruct_boundary(const __global float *wn,
										  const long id_cell,
										  const float vn[DIM],
										  const float sL[DIM][M], float wL[M])
{
	float r;

	long imemL = id_cell + 0 * NGRID;
#ifdef IS_2D
	wL[0] = wn[imemL];
	wL[0] += vn[0] * 0.5f * sL[0][0];
	wL[0] += vn[1] * 0.5f * sL[1][0];
#else
	wL[0] = wn[imemL];
	wL[0] += vn[0] * 0.5f * sL[0][0];
	wL[0] += vn[1] * 0.5f * sL[1][0];
	wL[0] += vn[2] * 0.5f * sL[2][0];
#endif

	const float rhoL = wn[id_cell + 0 * NGRID];
#pragma unroll
	for (int k = 1; k < M; k++) {
		r = get_r(rhoL, wn[id_cell + k * NGRID]);
		wL[k] = wL[0] * (r + vn[0] * 0.5f * sL[0][k] + vn[1] * 0.5f * sL[1][k] +
						 vn[2] * 0.5f * sL[2][k]);
	}
}
#endif

static float muscl_minmod_2(const float alpha, const float beta)
{
	if (alpha * beta < 0) {
		return 0.f;
	} else if (fabs(alpha) < fabs(beta)) {
		return alpha;
	} else {
		return beta;
	}
}

static float muscl_minmod(const float alpha, const float beta,
						  const float gamma)
{
	if (alpha > 0.f && beta > 0.f && gamma > 0.f) {
		return min(alpha, min(beta, gamma));
	}

	if (alpha < 0.f && beta < 0.f && gamma < 0.f) {
		return max(alpha, max(beta, gamma));
	}

	return 0.f;
}

static float muscl_minmod_3(const float alpha, const float beta,
							const float gamma)
{
	if (alpha > 0.f && beta > 0.f) {
		return min(alpha, min(beta, gamma));
	}

	if (alpha < 0.f && beta < 0.f) {
		return max(alpha, max(beta, gamma));
	}

	return 0.f;
}

static void muscl_get_slope(const long id_cell, const __global long *elem2elem,
							const __global float *wn, float s[DIM][M])
{
#ifdef IS_2D
	const int id_neighbours[FACE_PER_ELEM] = {
		elem2elem[FACE_PER_ELEM * id_cell + 0], // Right [ 1,  0]
		elem2elem[FACE_PER_ELEM * id_cell + 1], // Left  [-1,  0]
		elem2elem[FACE_PER_ELEM * id_cell + 2], // Front [ 0,  1]
		elem2elem[FACE_PER_ELEM * id_cell + 3], // Back  [ 0, -1]
	};

#else
	const int id_neighbours[FACE_PER_ELEM] = {
		elem2elem[FACE_PER_ELEM * id_cell + 0], // Right [ 1,  0,  0]
		elem2elem[FACE_PER_ELEM * id_cell + 1], // Left  [-1,  0,  0]
		elem2elem[FACE_PER_ELEM * id_cell + 2], // Front [ 0,  1,  0]
		elem2elem[FACE_PER_ELEM * id_cell + 3], // Back  [ 0, -1,  0]
		elem2elem[FACE_PER_ELEM * id_cell + 4], // North [ 0,  0,  1]
		elem2elem[FACE_PER_ELEM * id_cell + 5], // South [ 0,  0, -1]
	};
#endif

	float wM[M];
	float w1[M], w2[M];
	float a, b, g;
	long imem;
	long id_n1, id_n2;

/* Load current cell */
#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		imem = id_cell + iw * NGRID;
		wM[iw] = wn[imem];
	}

/* Load neighbours */
#pragma unroll
	for (int phy_dim = 0; phy_dim < DIM; phy_dim++) {
		/* Load two-by-two complementary neighbour value (R-L), (F-B), (N-S) */
		id_n1 = id_neighbours[2 * phy_dim + 0];
		id_n2 = id_neighbours[2 * phy_dim + 1];

		/* We use OR since if one is null s is null due to muscl_minmod */
#pragma unroll
		for (int iw = 0; iw < M; iw++) {
			/* w1 = R, F, N,  */
			if (id_n1 != -1) {
				imem = id_n1 + iw * NGRID;
				w1[iw] = wn[imem];
			} else {
				w1[iw] = wM[iw];
			}

			/* w2 = L, B, S */
			if (id_n2 != -1) {
				imem = id_n2 + iw * NGRID;
				w2[iw] = wn[imem];
			} else {
				w2[iw] = wM[iw];
			}

			a = (wM[iw] - w2[iw]);
			b = (w1[iw] - wM[iw]);
			s[phy_dim][iw] = muscl_minmod_2(a, b);

			// a = 2.f * (wM[iw] - w2[iw]);
			// b = 2.f * (w1[iw] - wM[iw]);
			// g = 0.5f * (w1[iw] - w2[iw]);
			// s[phy_dim][iw] = muscl_minmod_3(a, b, g);
		}

#ifdef IS_MOMENT_MODEL
		/* Slope rx, ry, rz */
		float RiM[3] = { 0, 0, 0 };
		float Ri1[3] = { 0, 0, 0 };
		float Ri2[3] = { 0, 0, 0 };
#pragma unroll
		for (int k = 0; k < 3; k++) {
			RiM[k] = get_r(wM[0], wM[k + 1]);
			Ri1[k] = get_r(w1[0], w1[k + 1]);
			Ri2[k] = get_r(w2[0], w2[k + 1]);

			a = (RiM[k] - Ri2[k]);
			b = (Ri1[k] - RiM[k]);
			s[phy_dim][k + 1] = muscl_minmod_2(a, b);

			// a = 2.f * (RiM[k] - Ri2[k]);
			// b = 2.f * (Ri1[k] - RiM[k]);
			// g = 0.5f * (Ri1[k] - Ri2[k]);
			// s[phy_dim][k + 1] = muscl_minmod_3(a, b, g);
		}
#endif
	}
}
#endif
