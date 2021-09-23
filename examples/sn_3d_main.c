#define M _m_
#define DX _dx_
#define DY _dy_
#define DZ _dz_
#define DT _dt_
#define NGRID _ngrid_
#define VOL (DX * DY * DZ)

/* User model parameters */
#define SIGMA _sigma_
#define ALPHA _alpha_
#define BETA _beta_
#define SRC_TOFF _src_toff_
#define SRC_X _src_x_
#define SRC_Y _src_y_
#define SRC_Z _src_z_
#define SRC_R _src_r_
#define SRC_V _src_v_

#include <solver/solver_3d.cl>
// #include <sn/sn.cl>

#define USE_S6_3D

#ifdef USE_S6_3D
#include <sn/s6_3d.cl>
#endif

#ifdef USE_S14_3D
#include "s14_3d.cl"
#endif

#ifdef USE_S38_3D
#include "s38_3d.cl"
#endif

void sn_num_flux_upwind(const float wL[M], const float wR[M], const float vn[3],
						float flux[M])
{
#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		const float v_dot_n = vn[0] * sn_vel[3 * iw + 0] +
							  vn[1] * sn_vel[3 * iw + 1] +
							  vn[2] * sn_vel[3 * iw + 2];

		const float vp = v_dot_n > 0.f ? v_dot_n : 0.f;
		const float vm = v_dot_n - vp;
		flux[iw] = vp * wL[iw] + vm * wR[iw];
	}
}

void sn_num_flux_boundary(const float wL[M], const float wR[M],
						  const float vn[3], float flux[M])
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
		v_dot_n = vn[0] * sn_vel[3 * iw] + vn[1] * sn_vel[3 * iw + 1] +
				  vn[2] * sn_vel[3 * iw + 2];

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
		v_dot_n = vn[0] * sn_vel[3 * iw] + vn[1] * sn_vel[3 * iw + 1] +
				  vn[2] * sn_vel[3 * iw + 2];

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
		Ix += w * sn_vel[3 * iw + 0];
		Iy += w * sn_vel[3 * iw + 1];
		Iz += w * sn_vel[3 * iw + 2];
	}

	wn_macro[id + 0 * NGRID] = w;
	wn_macro[id + 1 * NGRID] = Ix;
	wn_macro[id + 2 * NGRID] = Iy;
	wn_macro[id + 3 * NGRID] = Iz;
}

static void sn_src_circle_3d(const float x[3], const float t, float wn[M])
{
	const float d = (x[0] - SRC_X) * (x[0] - SRC_X) +
					(x[1] - SRC_Y) * (x[1] - SRC_Y) +
					(x[2] - SRC_Z) * (x[2] - SRC_Z);

	if (t < SRC_TOFF) {
		if (d > SRC_R * SRC_R) {
			for (int iw = 0; iw < M; iw++) {
				wn[iw] = 0.f;
			}
		} else {
			for (int iw = 0; iw < M; iw++) {
				wn[iw] = 1.f / SRC_V;
			}
		}
	} else {
		for (int iw = 0; iw < M; iw++) {
			wn[iw] = 0.f;
		}
	}
}

static void sn_source(const float x[3], const float t, float s[M])
{
	sn_src_circle_3d(x, t, s);
}

void vf_init_cond(const float x[3], const float t, float s[M])
{
	sn_source(x, t, s);
}

void vf_source(const float x[3], const float wn[M], const float t, float s[M])
{
#pragma unroll
	for (int iv = 0; iv < M; iv++) {
		s[iv] = 0;
	}
}

void vf_num_flux(const float wL[M], const float wR[M], const float vn[3],
				 float flux[M])
{
	sn_num_flux_upwind(wL, wR, vn, flux);
}

void vf_num_flux_boundary(const float wL[M], const float wR[M],
						  const float vn[3], float flux[M])
{
	// float wR_null[M];
	// for (int iw = 0; iw < M; iw++) {
	// 	wR_null[iw] = 0.f;
	// }
	// sn_num_flux_upwind(wL, wR_null, vn, flux);

	sn_num_flux_boundary(wL, wL, vn, flux);
}