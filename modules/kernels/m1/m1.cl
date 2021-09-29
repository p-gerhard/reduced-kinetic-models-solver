#ifndef M1_CL
#define M1_CL

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

#ifdef USE_M1_S1_CLOSURE
#undef USE_M1_S2_CLOSURE
#include <m1/m1_s1.cl>
#endif

#ifdef USE_KINETIC_NUM_FLUX_SPLIT
#define USE_QUAD_TRIGAUSS
#undef USE_QUAD_LEBEDEV
#endif

#ifdef USE_M1_S2_CLOSURE
#undef USE_M1_S1_CLOSURE
#include <m1/m1_s2.cl>
#endif

#include <solver/solver.cl>

static void m1_src_gaussian_line_source(const float x[DIM], const float t,
										float w[M])
{
#ifdef IS_2D
	const float d =
		(x[0] - SRC_X) * (x[0] - SRC_X) + (x[1] - SRC_Y) * (x[1] - SRC_Y);
#else
	const float d = (x[0] - SRC_X) * (x[0] - SRC_X) +
					(x[1] - SRC_Y) * (x[1] - SRC_Y) +
					(x[2] - SRC_Z) * (x[2] - SRC_Z);
#endif

	const float vaccum = 1e-4f;
	const float omega = 0.03f;
	const float omega2 = omega * omega;
	const float gauss = exp(-10 * d / (omega2));

#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		w[iw] = 0.f;
	}

	w[0] = max(gauss, vaccum);
}

static void m1_src_circle(const float x[DIM], const float t, float wn[M])
{
#ifdef IS_2D
	const float d =
		(x[0] - SRC_X) * (x[0] - SRC_X) + (x[1] - SRC_Y) * (x[1] - SRC_Y);

#else
	const float d = (x[0] - SRC_X) * (x[0] - SRC_X) +
					(x[1] - SRC_Y) * (x[1] - SRC_Y) +
					(x[2] - SRC_Z) * (x[2] - SRC_Z);
#endif

	const float vaccum = 1e-4f;

	if (t < SRC_TOFF) {
		if (d > SRC_R * SRC_R) {
			wn[0] = vaccum;
			wn[1] = 0.f;
			wn[2] = 0.f;
#ifndef IS_2D
			wn[3] = 0.f;
#endif
		} else {
			wn[0] = 1.f / SRC_V;
			wn[1] = 0.f;
			wn[2] = 0.f;
#ifndef IS_2D
			wn[3] = 0.f;
#endif
		}
	} else {
		wn[0] = vaccum;
		wn[1] = 0.f;
		wn[2] = 0.f;
#ifndef IS_2D
		wn[3] = 0.f;
#endif
	}
}

static void m1_source(const float x[DIM], const float t, float s[M])
{
	m1_src_gaussian_line_source(x, t, s);
	// m1_src_circle_3d(x, t, s);
}

void vf_init_cond(const float x[DIM], const float t, float s[M])
{
	m1_source(x, t, s);
}

void vf_source(const float x[DIM], const float wn[M], const float t, float s[M])
{
// m1_source(x, t, s);
#pragma unroll
	for (int iv = 0; iv < M; iv++) {
		s[iv] = 0.f;
	}
}

void vf_num_flux(const float wL[M], const float wR[M], const float vn[DIM],
				 float flux[M])
{
#ifdef USE_KINETIC_NUM_FLUX
	m1_num_flux_kinetic(wL, wR, vn, flux);
#elif defined(USE_KINETIC_NUM_FLUX_SPLIT)
	if (vn[0] == 1) {
		m1_num_flux_RL(wL, wR, vn, flux);
	}

	if (vn[0] == -1) {
		m1_num_flux_RL(wR, wL, vn, flux);
	}

	if (vn[1] == 1) {
		m1_num_flux_FB(wL, wR, vn, flux);
	}

	if (vn[1] == -1) {
		m1_num_flux_FB(wR, wL, vn, flux);
	}

	if (vn[2] == 1) {
		m1_num_flux_NS(wL, wR, vn, flux);
	}

	if (vn[2] == -1) {
		m1_num_flux_NS(wR, wL, vn, flux);
	}
#else
	m1_num_flux_rusanov(wL, wR, vn, flux);
#endif
}

void vf_num_flux_boundary(const float wL[M], const float wR[M],
						  const float vn[DIM], float flux[M])
{
	vf_num_flux(wL, wL, vn, flux);
	m1_num_flux_boundary(wL, wL, vn, flux);
}

#endif