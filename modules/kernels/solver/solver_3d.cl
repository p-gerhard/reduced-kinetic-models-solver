#ifndef SOLVER_3D_CL
#define SOLVER_3D_CL

#if USE_DOUBLE
#if defined(cl_khr_fp64)
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define DOUBLE_SUPPORT_AVAILABLE
#elif defined(cl_amd_fp64)
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#define DOUBLE_SUPPORT_AVAILABLE
#endif
#endif

#if defined(DOUBLE_SUPPORT_AVAILABLE)

/* Double */
typedef double real_t;
typedef double2 real2_t;
typedef double3 real3_t;
typedef double4 real4_t;
typedef double8 real8_t;
typedef double16 real16_t;
#define PI 3.14159265358979323846

#else

/* Float */
typedef float real_t;
typedef float2 real2_t;
typedef float3 real3_t;
typedef float4 real4_t;
typedef float8 real8_t;
typedef float16 real16_t;
#define PI 3.14159265359f

#endif

#ifdef IS_2D
#define DIM 2
#define FACE_PER_ELEM 4
#else
#define DIM 3
#define FACE_PER_ELEM 6
#endif

#define IS_M1_3D

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

/* Definition */
static float get_r_2(const float rho, const float In)
{
	float r;
	float tol = 1e-6;

	if (isless(rho, tol)) {
		r = 0;
	} else {
		r = In / rho;
		r = fmin(In / rho, 0.999f);
	}
	return r;
}

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
#else
	const int dir[FACE_PER_ELEM][DIM] = { { 1, 0, 0 }, { -1, 0, 0 },
										  { 0, 1, 0 }, { 0, -1, 0 },
										  { 0, 0, 1 }, { 0, 0, -1 } };

	const float ds[FACE_PER_ELEM] = { DY * DZ, DY * DZ, DX * DY,
									  DX * DY, DZ * DY, DZ * DY };
#endif

	float flux[M];
	float wL[M];
	float wR[M];
	float wnext[M];

	/* MUSCL slopes */
	float sL[DIM][M];
	float sR[DIM][M];

#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		long imemL = id_cell + iw * NGRID;
		wnext[iw] = wn[imemL];
		sL[0][iw] = 0.f;
		sL[1][iw] = 0.f;

		sR[0][iw] = 0.f;
		sR[1][iw] = 0.f;

#ifndef IS_2D
		sL[2][iw] = 0.f;
		sR[2][iw] = 0.f;
#endif
	}

	/* MUSCL slope of the current element (L) */
	muscl_get_slope(id_cell, elem2elem, wn, sL);

	/* Stencil */
#pragma unroll
	for (int id_face = 0; id_face < FACE_PER_ELEM; id_face++) {
#ifdef IS_2D
		float vn[DIM] = { dir[id_face][0], dir[id_face][1] };
#else
		float vn[DIM] = { dir[id_face][0], dir[id_face][1], dir[id_face][2] };
#endif

		const long id_neighbour = elem2elem[FACE_PER_ELEM * id_cell + id_face];

		if (id_neighbour != -1) {
			muscl_get_slope(id_neighbour, elem2elem, wn, sR);
#pragma unroll
			for (int iw = 0; iw < M; iw++) {
				long imemR = id_neighbour + iw * NGRID;
				long imemL = id_cell + iw * NGRID;

#ifdef IS_2D
				wR[iw] = wn[imemR] - vn[0] * 0.5f * sR[0][iw] -
						 vn[1] * 0.5f * sR[1][iw];

				wL[iw] = wn[imemL] + vn[0] * 0.5f * sL[0][iw] +
						 vn[1] * 0.5f * sL[1][iw];
#else
				wR[iw] = wn[imemR] - vn[0] * 0.5f * sR[0][iw] -
						 vn[1] * 0.5f * sR[1][iw] - vn[2] * 0.5f * sR[2][iw];

				wL[iw] = wn[imemL] + vn[0] * 0.5f * sL[0][iw] +
						 vn[1] * 0.5f * sL[1][iw] + vn[2] * 0.5f * sL[2][iw];
			}

#ifdef IS_M1_3D
			// 				/* If wL[0] < 0 -> wn[imemL] */
			// if(isless(wL[0], 0.f)) {
			// 	wL[0] = wn[id_cell + 0 * NGRID];
			// }

			// 				if(isless(wR[0], 0.f)) {
			// 					wR[0] = wn[id_neighbour + 0 * NGRID];
			// 				}

			/* Reconstruction on Ik */
			for (int k = 1; k < 4; k++) {
				float ri_locL =
					get_r_2(wn[id_cell + 0 * NGRID], wn[id_cell + k * NGRID]);
				// if(ri_locL > 0.5) {
				// 	printf("riL = %f\n", ri_locL);
				// }
				wL[k] =
					wL[0] * (ri_locL + vn[0] * 0.5f * sL[0][k] +
							 vn[1] * 0.5f * sL[1][k] + vn[2] * 0.5f * sL[2][k]);

				// if(isnan(wL[0])) {
				// 		printf("riL = %f\n", ri_locL);
				// }
				float ri_locR = get_r_2(wn[id_neighbour + 0 * NGRID],
										wn[id_neighbour + k * NGRID]);

				wR[k] =
					wR[0] * (ri_locR - vn[0] * 0.5f * sR[0][k] -
							 vn[1] * 0.5f * sR[1][k] - vn[2] * 0.5f * sR[2][k]);

				// if(isnan(wR[0])) {
				// 		printf("riR = %f\n", ri_locR);
				// }
			}
#endif
#endif
				vf_num_flux(wL, wR, vn, flux);
			}
			else
			{
#pragma unroll
				for (int iw = 0; iw < M; iw++) {
					long imemL = id_cell + iw * NGRID;

#ifdef IS_2D
					wL[iw] = wn[imemL] + vn[0] * 0.5f * sL[0][iw] +
							 vn[1] * 0.5f * sL[1][iw];
#else
				wL[iw] = wn[imemL] + vn[0] * 0.5f * sL[0][iw] +
						 vn[1] * 0.5f * sL[1][iw] + vn[2] * 0.5f * sL[2][iw];
#endif
				}
#ifdef IS_M1_3D
				/* Reconstruction on Ik */
				for (int k = 1; k < 4; k++) {
					float ri_locL = get_r_2(wn[id_cell + 0 * NGRID],
											wn[id_cell + k * NGRID]);
					// if(ri_locL > 0.5) {
					// 	printf("riL = %f\n", ri_locL);
					// }
					wL[k] = wL[0] *
							(ri_locL + vn[0] * 0.5f * sL[0][k] +
							 vn[1] * 0.5f * sL[1][k] + vn[2] * 0.5f * sL[2][k]);
				}
#endif
				vf_num_flux_boundary(wL, wL, vn, flux);
			}

			/* Time evolution */
#pragma unroll
			for (int iw = 0; iw < M; iw++) {
				wnext[iw] -= DT * (ds[id_face] / VOL) * flux[iw];
			}
		}

/* Source term */
#ifdef IS_2D
		const float xyz[DIM] = { x[DIM * id_cell + 0], x[DIM * id_cell + 1] };
#else
	const float xyz[DIM] = { x[DIM * id_cell + 0], x[DIM * id_cell + 1],
							 x[DIM * id_cell + 2] };
#endif
		float src[M];
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