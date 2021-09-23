#ifndef SOLVER_3D_CL
#define SOLVER_3D_CL

// #define USE_MUSCL

#ifdef IS_2D
#define DIM 2
#define FACE_PER_ELEM 4
#else
#define DIM 3
#define FACE_PER_ELEM 6
#endif

#include "muscl.cl"

void vf_source(const float x[DIM], const float wn[M], const float t,
			   float s[M]);

void vf_init_cond(const float x[DIM], const float t, float s[M]);

void vf_num_flux(const float wL[M], const float wR[M], const float vn[DIM],
				 float flux[M]);

void vf_num_flux_boundary(const float wL[M], const float wR[M],
						  const float vn[DIM], float flux[M]);

__kernel void vf_init_sol(__global const float *x, __global float *wn);

__kernel void solver_3d_muscl_euler_time_step(const float tnow,
											  __global const float *x,
											  __global const long *elem2elem,
											  __global const float *wn,
											  __global float *wnp1);

__kernel void vf_init_sol(__global const float *x, __global float *wn)
{
	const long id = get_global_id(0);
	long imem;
	float wnow[M];

#ifdef IS_2D
	const float xyz[DIM] = { x[DIM * id], x[DIM * id + 1] };
#else
	const float xyz[DIM] = { x[DIM * id], x[DIM * id + 1], x[DIM * id + 2] };
#endif

	vf_init_cond(xyz, 0.f, wnow);

#pragma unroll
	for (unsigned int iw = 0; iw < M; iw++) {
		imem = id + iw * NGRID;
		wn[imem] = wnow[iw];
	}
}

__kernel void solver_muscl_euler_time_step(const float tnow,
										   __global const float *x,
										   __global const long *elem2elem,
										   __global const float *wn,
										   __global float *wnp1)
{
	const long id_cell = get_global_id(0);

#ifdef IS_2D
	const int dir[FACE_PER_ELEM][DIM] = {
		{ 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 }
	};
	const float ds[FACE_PER_ELEM] = { DY, DY, DX, DX };
	const float xyz[DIM] = { x[DIM * id_cell + 0], x[DIM * id_cell + 1] };
#else
	const int dir[FACE_PER_ELEM][DIM] = { { 1, 0, 0 }, { -1, 0, 0 },
										  { 0, 1, 0 }, { 0, -1, 0 },
										  { 0, 0, 1 }, { 0, 0, -1 } };

	const float ds[FACE_PER_ELEM] = { DY * DZ, DY * DZ, DX * DY,
									  DX * DY, DZ * DY, DZ * DY };

	const float xyz[DIM] = { x[DIM * id_cell + 0], x[DIM * id_cell + 1],
							 x[DIM * id_cell + 2] };
#endif

	float flux[M];
	float wL[M];
	float wR[M];
	float wnext[M];
	float src[M];

	/* MUSCL slopes */
#ifdef USE_MUSCL
	float sL[DIM][M];
	float sR[DIM][M];
#endif

#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		long imemL = id_cell + iw * NGRID;
		wnext[iw] = wn[imemL];
#ifdef USE_MUSCL
		sL[0][iw] = 0.f;
		sL[1][iw] = 0.f;

		sR[0][iw] = 0.f;
		sR[1][iw] = 0.f;

#ifndef IS_2D
		sL[2][iw] = 0.f;
		sR[2][iw] = 0.f;
#endif
#endif
	}

	/* MUSCL slope of the current element (L) */
#ifdef USE_MUSCL
	muscl_get_slope(id_cell, elem2elem, wn, sL);
#endif

// /* Stencil */
#pragma unroll
	for (int id_face = 0; id_face < FACE_PER_ELEM; id_face++) {
#ifdef IS_2D
		float vn[DIM] = { dir[id_face][0], dir[id_face][1] };
#else
		float vn[DIM] = { dir[id_face][0], dir[id_face][1], dir[id_face][2] };
#endif
		const long id_neighbour = elem2elem[FACE_PER_ELEM * id_cell + id_face];

		if (id_neighbour != -1) {
#ifdef USE_MUSCL
			/* MUSCL slope of the neighbour element (R) */
			muscl_get_slope(id_neighbour, elem2elem, wn, sR);
#ifdef IS_MOMENT_MODEL
			muscl_m1_reconstruct(wn, id_cell, id_neighbour, vn, sL, sR, wL, wR);
#else
			muscl_reconstruct(wn, id_cell, id_neighbour, vn, sL, sR, wL, wR);
#endif
#else
#pragma unroll
			for (int iw = 0; iw < M; iw++) {
				wL[iw] = wn[id_cell + iw * NGRID];
				wR[iw] = wn[id_neighbour + iw * NGRID];
			}
#endif
			vf_num_flux(wL, wR, vn, flux);
		} else {
#ifdef USE_MUSCL
#ifdef IS_MOMENT_MODEL
			muscl_m1_reconstruct_boundary(wn, id_cell, vn, sL, wL);
#else
			muscl_reconstruct_boundary(wn, id_cell, vn, sL, wL);
#endif
#else
#pragma unroll
			for (int iw = 0; iw < M; iw++) {
				wL[iw] = wn[id_cell + iw * NGRID];
			}
#endif
			vf_num_flux_boundary(wL, wL, vn, flux);
		}
		/* Time evolution */
		for (int iw = 0; iw < M; iw++) {
			wnext[iw] -= DT * (ds[id_face] / VOL) * flux[iw];
		}
	}

	/* Source term */
	vf_source(xyz, wL, tnow, src);
#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		wnext[iw] += src[iw] * DT;
	}

/* Update main global vector wnp1 */
#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		long imemL = id_cell + iw * NGRID;
		wnp1[imemL] = wnext[iw];
	}
}
#endif