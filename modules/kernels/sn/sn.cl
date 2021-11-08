#ifndef SN_CL
#define SN_CL
#define DIM _dim_
#define M _m_
#define DX _dx_
#define DY _dy_

#if DIM == 3
#define DZ _dz_
#endif

#define DT _dt_
#define NGRID _ngrid_

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

#include <solver/solver.cl>

#if defined(USE_QUAD_UNIFORM) && DIM == 2
#if M == 4
#define USE_QUAD_UNIFORM_NB_V_4
#elif M == 8
#define USE_QUAD_UNIFORM_NB_V_8
#elif M == 16
#define USE_QUAD_UNIFORM_NB_V_16
#elif M == 32
#define USE_QUAD_UNIFORM_NB_V_32
#elif M == 64
#define USE_QUAD_UNIFORM_NB_V_64
#elif M == 128
#define USE_QUAD_UNIFORM_NB_V_128
#elif M == 256
#define USE_QUAD_UNIFORM_NB_V_256
#elif M == 512
#define USE_QUAD_UNIFORM_NB_V_512
#endif
#include <maths/quad_uniform.cl>
#endif

#if defined(USE_QUAD_LEBEDEV) && DIM == 3
#if M == 6
#define USE_QUAD_LEBEDEV_NB_V_6
#elif M == 14
#define USE_QUAD_LEBEDEV_NB_V_14
#elif M == 26
#define USE_QUAD_LEBEDEV_NB_V_26
#elif M == 38
#define USE_QUAD_LEBEDEV_NB_V_38
#elif M == 50
#define USE_QUAD_LEBEDEV_NB_V_50
#elif M == 74
#define USE_QUAD_LEBEDEV_NB_V_74
#elif M == 86
#define USE_QUAD_LEBEDEV_NB_V_86
#elif M == 110
#define USE_QUAD_LEBEDEV_NB_V_110
#elif M == 146
#define USE_QUAD_LEBEDEV_NB_V_146
#elif M == 170
#define USE_QUAD_LEBEDEV_NB_V_170
#elif M == 194
#define USE_QUAD_LEBEDEV_NB_V_194
#elif M == 230
#define USE_QUAD_LEBEDEV_NB_V_230
#elif M == 266
#define USE_QUAD_LEBEDEV_NB_V_266
#elif M == 302
#define USE_QUAD_LEBEDEV_NB_V_302
#elif M == 350
#define USE_QUAD_LEBEDEV_NB_V_350
#elif M == 434
#define USE_QUAD_LEBEDEV_NB_V_434
#elif M == 590
#define USE_QUAD_LEBEDEV_NB_V_590
#elif M == 770
#define USE_QUAD_LEBEDEV_NB_V_770
#endif
#include <maths/quad_lebedev.cl>
#endif

void sn_num_flux_upwind(const float wL[M], const float wR[M],
						const float vn[DIM], float flux[M])
{
#pragma unroll
	for (int iw = 0; iw < M; iw++) {
#ifdef IS_2D
		const float v_dot_n =
			vn[0] * quad_vi[DIM * iw + 0] + vn[1] * quad_vi[DIM * iw + 1];
#else
		const float v_dot_n = vn[0] * quad_vi[DIM * iw + 0] +
							  vn[1] * quad_vi[DIM * iw + 1] +
							  vn[2] * quad_vi[DIM * iw + 2];
#endif

		const float vp = v_dot_n > 0.f ? v_dot_n : 0.f;
		const float vm = v_dot_n - vp;
		flux[iw] = vp * wL[iw] + vm * wR[iw];
	}
}

void sn_num_flux_boundary(const float wL[M], const float wR[M],
						  const float vn[DIM], float flux[M])
{
	float sum_v_dot_n = 0.f;
	float sum_flux_out = 0.f;

	
	/* Compute outgoing flux */
#pragma unroll
	for (int iw = 0; iw < M; iw++) {
#ifdef IS_2D
		const float v_dot_n =
			vn[0] * quad_vi[DIM * iw + 0] + vn[1] * quad_vi[DIM * iw + 1];
#else
		const float v_dot_n = vn[0] * quad_vi[DIM * iw + 0] +
							  vn[1] * quad_vi[DIM * iw + 1] +
							  vn[2] * quad_vi[DIM * iw + 2];
#endif
		if (v_dot_n > 0) {
			flux[iw] = v_dot_n * wL[iw];
			sum_flux_out += quad_wi[iw] * v_dot_n * wL[iw];
			sum_v_dot_n += v_dot_n * quad_wi[iw];
		} else {
			flux[iw] = 0.f;
		}
	}

	/* Compute ingoing flux */
#pragma unroll
	for (int iw = 0; iw < M; iw++) {
#ifdef IS_2D
		const float v_dot_n =
			vn[0] * quad_vi[DIM * iw + 0] + vn[1] * quad_vi[DIM * iw + 1];
#else
		const float v_dot_n = vn[0] * quad_vi[DIM * iw + 0] +
							  vn[1] * quad_vi[DIM * iw + 1] +
							  vn[2] * quad_vi[DIM * iw + 2];
#endif
		if (v_dot_n < 0) {
			/* Multiplication by the normal's component to avoid branching.
             * Assuming normal vector component is always 1 or 0.
             */
#ifdef IS_2D
			float flux_spec = v_dot_n * wL[quad_spec[0][iw]] * fabs(vn[0]) +
							  v_dot_n * wL[quad_spec[1][iw]] * fabs(vn[1]);
#else
			float flux_spec = v_dot_n * wL[quad_spec[0][iw]] * fabs(vn[0]) +
							  v_dot_n * wL[quad_spec[1][iw]] * fabs(vn[1]) +
							  v_dot_n * wL[quad_spec[2][iw]] * fabs(vn[2]);
#endif
			float flux_diff = v_dot_n * sum_flux_out / sum_v_dot_n;

			flux[iw] = ALPHA * (BETA * flux_diff + (1.f - BETA) * flux_spec);
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

#ifndef IS_2D
	float Iz = 0.f;
#endif

/* Integration over S1 or S2 */
#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		long imem = id + iw * NGRID;

#ifdef IS_2D
		w += QUAD_WI * wn[imem];
#else
		/* S2 Lebedev quadrature using quad_wi */
		w += quad_wi[iw] * wn[imem];
#endif
		Ix += w * quad_vi[DIM * iw + 0];
		Iy += w * quad_vi[DIM * iw + 1];

#ifndef IS_2D
		Iz += w * quad_vi[DIM * iw + 2];
#endif
	}

	wn_macro[id + 0 * NGRID] = w;
	wn_macro[id + 1 * NGRID] = Ix;
	wn_macro[id + 2 * NGRID] = Iy;

#ifndef IS_2D
	wn_macro[id + 3 * NGRID] = Iz;
#endif
}

static void sn_src_circle(const float x[DIM], const float t, float wn[M])
{
#ifdef IS_2D
	const float d =
		(x[0] - SRC_X) * (x[0] - SRC_X) + (x[1] - SRC_Y) * (x[1] - SRC_Y);
#else
	const float d = (x[0] - SRC_X) * (x[0] - SRC_X) +
					(x[1] - SRC_Y) * (x[1] - SRC_Y) +
					(x[2] - SRC_Z) * (x[2] - SRC_Z);
#endif

	if (t < SRC_TOFF) {
		if (d > SRC_R * SRC_R) {
			for (int iw = 0; iw < M; iw++) {
				wn[iw] = 0.f;
			}
		} else {
			for (int iw = 0; iw < M; iw++) {
					wn[iw] = quad_wi[iw];
			}
		}
	} else {
		for (int iw = 0; iw < M; iw++) {
			wn[iw] = 0.f;
		}
	}
}

static void sn_source(const float x[DIM], const float t, float s[M])
{
	sn_src_circle(x, t, s);
}

void vf_init_cond(const float x[DIM], const float t, float s[M])
{
	sn_source(x, t, s);
}

void vf_source(const float x[DIM], const float wn[M], const float t, float s[M])
{
#pragma unroll
	for (int iv = 0; iv < M; iv++) {
		s[iv] = 0.f;
	}
}

void vf_num_flux(const float wL[M], const float wR[M], const float vn[DIM],
				 float flux[M])
{
	sn_num_flux_upwind(wL, wR, vn, flux);
}

void vf_num_flux_boundary(const float wL[M], const float wR[M],
						  const float vn[DIM], float flux[M])
{
	sn_num_flux_boundary(wL, wL, vn, flux);
}
#endif