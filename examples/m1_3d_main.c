#define M 4
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

#define IS_MOMENT_MODEL

#include <m1/m1_s2.cl>
#include <solver/solver_3d.cl>

static void src_gaussian_line_source_3d(const float x[DIM], const float t,
										float w[M])
{
	const float d = (x[0] - SRC_X) * (x[0] - SRC_X) +
					(x[1] - SRC_Y) * (x[1] - SRC_Y) +
					(x[2] - SRC_Z) * (x[2] - SRC_Z);

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

static void m1_src_circle_3d(const float x[DIM], const float t, float wn[M])
{
	const float d = (x[0] - SRC_X) * (x[0] - SRC_X) +
					(x[1] - SRC_Y) * (x[1] - SRC_Y) +
					(x[2] - SRC_Z) * (x[2] - SRC_Z);

	const float vaccum = 1e-4f;

	if (t < SRC_TOFF) {
		if (d > SRC_R * SRC_R) {
			wn[0] = vaccum;
			wn[1] = 0.f;
			wn[2] = 0.f;
			wn[3] = 0.f;
		} else {
			wn[0] = 1.f / SRC_V;
			wn[1] = 0.f;
			wn[2] = 0.f;
			wn[3] = 0.f;
		}
	} else {
		wn[0] = vaccum;
		wn[1] = 0.f;
		wn[2] = 0.f;
		wn[3] = 0.f;
	}
}

static void m1_source(const float x[DIM], const float t, float s[M])
{
	src_gaussian_line_source_3d(x, t, s);
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
		s[iv] = 0;
	}
}

void vf_num_flux(const float wL[M], const float wR[M], const float vn[DIM],
				 float flux[M])
{
	m1_num_flux_rusanov(wL, wR, vn, flux);
	// m1_num_flux_kinetic_lebedev(wL, wR, vn, flux);
}

void vf_num_flux_boundary(const float wL[M], const float wR[M],
						  const float vn[DIM], float flux[M])
{
	vf_num_flux(wL, wL, vn, flux);
}