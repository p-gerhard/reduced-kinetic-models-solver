#ifndef PN_CL
#define PN_CL
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

#ifdef USE_SPHERICAL_HARMONICS_P1
#include <pn/p1.cl>
#endif

#ifdef USE_SPHERICAL_HARMONICS_P3
#include <pn/p3.cl>
#endif

#ifdef USE_SPHERICAL_HARMONICS_P5
#include <pn/p5.cl>
#endif

#ifdef USE_SPHERICAL_HARMONICS_P7
#include <pn/p7.cl>
#endif

#ifdef USE_SPHERICAL_HARMONICS_P9
#include <pn/p9.cl>
#endif

#ifdef USE_SPHERICAL_HARMONICS_P11
#include <pn/p11.cl>
#endif


static void pn_src_gaussian_line_source(const float x[DIM], const float t,
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

void pn_src_circle(const float x[DIM], const float t, float w[M])
{
#ifdef IS_2D
	const float d =
		(x[0] - SRC_X) * (x[0] - SRC_X) + (x[1] - SRC_Y) * (x[1] - SRC_Y);
#else
	const float d = (x[0] - SRC_X) * (x[0] - SRC_X) +
					(x[1] - SRC_Y) * (x[1] - SRC_Y) +
					(x[2] - SRC_Z) * (x[2] - SRC_Z);
#endif

#pragma unroll
	for (int iw = 0; iw < M; iw++) {
		w[iw] = 0.f;
	}

	if (t < SRC_TOFF) {
		if (d <= SRC_R * SRC_R) {
			w[0] = 1.f / SRC_V;
		}
	}
}

void pn_source(const float x[DIM], const float t, float s[M])
{
	pn_src_circle(x, t, s);
}

void vf_init_cond(const float x[DIM], const float t, float s[M])
{
	pn_source(x, t, s);
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
#ifdef USE_KINETIC_NUM_FLUX
	if (vn[0] == 1) {
		pn_num_flux_kinetic_R(wL, wR, flux);
	}

	if (vn[0] == -1) {
		pn_num_flux_kinetic_L(wL, wR, flux);
	}

	if (vn[1] == 1) {
		pn_num_flux_kinetic_F(wL, wR, flux);
	}

	if (vn[1] == -1) {
		pn_num_flux_kinetic_B(wL, wR, flux);
	}
#ifndef IS_2D
	if (vn[2] == 1) {
		pn_num_flux_kinetic_N(wL, wR, flux);
	}

	if (vn[2] == -1) {
		pn_num_flux_kinetic_S(wL, wR, flux);
	}
#endif
#else
	pn_num_flux_rusanov(wL, wR, vn, flux);
#endif
}

void vf_num_flux_boundary(const float wL[M], const float wR[M],
						  const float vn[DIM], float flux[M])
{
	if (vn[0] == 1) {
		pn_num_flux_boundary_kinetic_R(wL, flux);
	}

	if (vn[0] == -1) {
		pn_num_flux_boundary_kinetic_L(wL, flux);
	}
	if (vn[1] == 1) {
		pn_num_flux_boundary_kinetic_F(wL, flux);
	}

	if (vn[1] == -1) {
		pn_num_flux_boundary_kinetic_B(wL, flux);
	}
#ifdef USE_2D
	if (vn[2] == 1) {
		pn_num_flux_boundary_kinetic_N(wL, flux);
	}

	if (vn[2] == -1) {
		pn_num_flux_boundary_kinetic_S(wL, flux);
	}
#endif
}
#endif