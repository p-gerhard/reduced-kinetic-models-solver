#ifndef SN_CL
#define SN_CL

#define USE_S6_3D

#ifdef USE_S6_3D
#include "s6_3d.cl"
#endif

#ifdef USE_S14_3D
#include "s14_3d.cl"
#endif

#ifdef USE_S38_3D
#include "s38_3d.cl"
#endif

void sn_num_flux_upwind(const float wL[M], const float wR[M],
						const float vn[DIM], float flux[M])
{
#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		const float v_dot_n = vn[0] * sn_vel[DIM * iw] +
							  vn[1] * sn_vel[DIM * iw + 1] +
							  vn[2] * sn_vel[DIM * iw + 2];

		const float vp = v_dot_n > 0.f ? v_dot_n : 0.f;
		const float vm = v_dot_n - vp;
		flux[iw] = vp * wL[iw] + vm * wR[iw];
	}
}

void sn_num_flux_boundary(const float wL[M], const float wR[M],
						  const float vn[DIM], float flux[M])
{
	float flux_spec = 0.f;
	float flux_diff = 0.f;
	float flux_out = 0.f;
	float v_dot_n;
	float w_sum = 0.f;

#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		flux[iw] = 0.f;
	}

	/* Compute outgoing flux */
#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		v_dot_n = vn[0] * sn_vel[DIM * iw] + vn[1] * sn_vel[DIM * iw + 1] +
				  vn[2] * sn_vel[DIM * iw + 2];

		if (v_dot_n > 0.f) {
			flux_out += v_dot_n * wL[iw];
			flux[iw] += v_dot_n * wL[iw];
		} else {
			w_sum -= v_dot_n;
		}
	}

	/* Compute ingoing flux */
#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		v_dot_n = vn[0] * sn_vel[DIM * iw] + vn[1] * sn_vel[DIM * iw + 1] +
				  vn[2] * sn_vel[DIM * iw + 2];

		if (v_dot_n <= 0.f) {
			flux_diff = flux_out * v_dot_n / w_sum;

			/* Multiplication by the normal's component to avoid branching.
             * Assuming normal vector component is always 1 or 0.
             */
			flux_spec = v_dot_n * wL[sn_spec[0][iw]] * fabs(vn[0]) +
						v_dot_n * wL[sn_spec[1][iw]] * fabs(vn[1]) +
						v_dot_n * wL[sn_spec[2][iw]] * fabs(vn[2]);

			flux[iw] += ALPHA * (BETA * flux_diff + (1.f - BETA) * flux_spec);
		}
	}
}

__kernel void sn_compute_macro(__global const float *wn,
							   __global float *wn_macro)
{
	unsigned long id = get_global_id(0);
	float w = 0.f;
	float Ix = 0.f;
	float Iy = 0.f;
	float Iz = 0.f;

/* Integration over S1 or S2 */
#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		long imem = id + iw * NGRID;

		/* S2 Lebedev quadrature using sn_wgt */
		w += sn_wgt[iw] * wn[imem];
		Ix += w * sn_vel[DIM * iw];
		Iy += w * sn_vel[DIM * iw + 1];
		Iz += w * sn_vel[DIM * iw + 2];
	}

	wn_macro[id + 0 * NGRID] = w;
	wn_macro[id + 1 * NGRID] = Ix;
	wn_macro[id + 2 * NGRID] = Iy;
	wn_macro[id + 3 * NGRID] = Iz;
}
#endif