#ifndef MODEL_M1_S1_SIMPLE_CL
#define MODEL_M1_S1_SIMPLE_CL

#pragma once

#include "bessel_i.cl"
#include "rational.cl"
#include "quad_circle_uniform.cl"
#include "quad_tanh_sinh.cl"
#define MAX_K_SIMPLE 60.F
#define TOL_I 1E-8F
#define TOL_R 1E-8F
#define VACCUM 1E-4F

__constant static const float chi_p[7] = { 0.500000947627, -1.944182261747,
																					 2.930698162883, -2.177908221701,
																					 0.921417673004, -0.311258505782,
																					 0.081997371844 };

__constant static const float chi_q[7] = { 1.000000000000,	-3.888203179227,
																					 5.358815439044,	-2.394912664257,
																					 -1.057275075050, 1.306311165212,
																					 -0.323970518048 };

__constant static const float log_bessel_i0_p[7] = {
	0.722493563857652, -1.369083176985558, 1.076125440349176, 0.336568032124439,
	0.018918132550194, 0.000278937058459,	 0.000000916466550
};

__constant static const float log_bessel_i0_q[7] = {
	1.000000000000000, 1.887041149652891, 0.390686560521634, 0.019873797772182,
	0.000282852601570, 0.000000917124043, -0.000000000000335
};

inline static float get_angle_normalized(const float u[2])
{
	/* Return the angle between a normalized vector u and e_x = [1, 0] */
	return atan2(u[1], u[0]) - atan2(0.f, 1.f);
}

inline static float get_angle_non_normalized(const float u[2])
{
	/* Return the angle between a vector u and e_x = [1, 0] */
	const float u_norm = sqrt(u[0] * u[0] + u[1] * u[1]);
	if (isless(u_norm, 1e-8f)) {
		return 0;
	}
	return atan2(u[1] / u_norm, u[0] / u_norm) - atan2(0.f, 1.f);
}

static void flux_int(const float x, const float kn, const float log_i0,
										 const float th_k, const float th_n, float res[3])
{
	const float t0 = kn * cos(x - th_k) - log_i0;
	const float t1 = cos(x - th_n);
	const float t2 = exp(t0) * t1;

	res[0] = t2;
	res[1] = t2 * cos(x);
	res[2] = t2 * sin(x);
}

static void flux_integrand(float theta, const float kn, const float log_i0,
													 const float th_k, const float th_n, float value[3])
{
	const float t0 = kn * cos(theta - th_k) - log_i0;
	const float t1 = cos(theta - th_n);
	const float t2 = exp(t0) * t1;

	value[0] = t2;
	value[1] = t2 * cos(theta);
	value[2] = t2 * sin(theta);
}

static float get_r(const float *w, const float tol, bool *is_max)
{
	float rho = fabs(w[0]);
	const float In = hypot(w[1], w[2]);

	if (isless(rho, tol)) {
		rho = tol;
	}

	if (isless(In, rho * (1.f - tol))) {
		return In / rho;
	} else {
		*is_max = true;
		return 1.f - tol;
	}
}

static float get_chi_interp(const float *wn)
{
	float chi;
	bool is_max = false;
	float r = get_r(wn, TOL_R, &is_max);

	if (is_max) {
		return 1.f - TOL_R;
	}

	chi = evaluate_rational_simple(chi_p, chi_q, 7, r);
	chi = fmin(chi, 1.f - TOL_R);
	return chi;
}

static void get_kappa(const float *w, float *k, float *k_norm)
{
	/* Get kappa norm using interpolation */

	bool is_max = false;
	const float r = get_r(w, TOL_R, &is_max);
	const float chi = get_chi_interp(w);

	*k_norm = fmin(r / (1.f - chi), MAX_K_SIMPLE);

	const float Ix = w[1];
	const float Iy = w[2];
	const float In = hypot(Ix, Iy);

	if (isgreater(In, TOL_I)) {
		k[0] = *k_norm * Ix / In;
		k[1] = *k_norm * Iy / In;
	} else {
		*k_norm = 0.f;
		k[0] = 0.f;
		k[1] = 0.f;
	}
}

inline static void get_kappa_specular(const float *kL, const float *vn,
																			float *kR)
{
	/* Construct kR according to the specular law */
	const float kL_dot_vn = kL[0] * vn[0] + kL[1] * vn[1];
	kR[0] = kL[0] - 2.f * kL_dot_vn * vn[0];
	kR[1] = kL[1] - 2.f * kL_dot_vn * vn[1];
}

inline static float get_log_bessel_i0(const float k)
{
	if (isless(k, 15.f)) {
		return log(bessel_i0_simple(k));
	} else {
		float val =
			evaluate_rational_simple(log_bessel_i0_p, log_bessel_i0_q, 7, k);
		return val;
	}
}

static void phy_flux(const float *w, const float *vn, float *flux)
{
	const float rho = w[0];
	const float Ix = w[1];
	const float Iy = w[2];

	const float I2 = Ix * Ix + Iy * Iy;
	float P11, P12, P21, P22;

	if (isless(I2, TOL_I)) {
		const float chi = get_chi_interp(w);
		// P11 = w * 0.5; 			*/ TOCHECK chi should be 0.5 always here */
		P11 = rho * (1.f - chi);
		flux[0] = 0.f;
		flux[1] = P11 * vn[0];
		flux[2] = P11 * vn[1];
	} else {
		const float chi = get_chi_interp(w);
		const float IxIx = (Ix * Ix) / I2;
		const float IxIy = (Ix * Iy) / I2;
		const float IyIy = (Iy * Iy) / I2;

		const float t1 = (1.f - chi);
		const float t2 = (2.f * chi - 1.f);

		P11 = rho * (t1 + t2 * IxIx);
		P22 = rho * (t1 + t2 * IyIy);
		P12 = rho * t2 * IxIy;
		P21 = P12;

		flux[0] = Ix * vn[0] + Iy * vn[1];
		flux[1] = P11 * vn[0] + P12 * vn[1];
		flux[2] = P21 * vn[0] + P22 * vn[1];
	}
}

/* Numerical fluxes */

void m1_s1_num_flux_rusanov(const float *wL, const float *wR, const float *vn,
														float *flux)
{
	float fL[3];
	float fR[3];
	const float lambda = 1.f;

	phy_flux(wL, vn, fL);
	phy_flux(wR, vn, fR);

	flux[0] = 0.5f * (fL[0] + fR[0]) - 0.5f * lambda * (wR[0] - wL[0]);
	flux[1] = 0.5f * (fL[1] + fR[1]) - 0.5f * lambda * (wR[1] - wL[1]);
	flux[2] = 0.5f * (fL[2] + fR[2]) - 0.5f * lambda * (wR[2] - wL[2]);
}

void m1_s1_num_flux_kin(const float *wL, const float *wR, const float *vn,
												float *flux)
{
	float vx, vy, t, v_dot_n;
	float rhoL, rhoR, kLn, kRn, wifLvn, wifRvn;

	float kL[2];
	float kR[2];

	flux[0] = 0.f;
	flux[1] = 0.f;
	flux[2] = 0.f;

	rhoL = wL[0];
	rhoR = wR[0];

	get_kappa(wL, kL, &kLn);
	get_kappa(wR, kR, &kRn);

	const float log_i0_L = get_log_bessel_i0(kLn);
	const float log_i0_R = get_log_bessel_i0(kRn);

	__constant const float *vi = quad_circle_uniform_64_vi;
	const float wi = quad_circle_uniform_64_wi;
	const int N = quad_circle_uniform_64_N;

	for (unsigned int k = 0; k < N; k++) {
		vx = vi[2 * k];
		vy = vi[2 * k + 1];
		v_dot_n = vx * vn[0] + vy * vn[1];

		if (isgreater(v_dot_n, 0.f)) {
			t = kL[0] * vx + kL[1] * vy - log_i0_L;
			wifLvn = wi * rhoL * exp(t) * v_dot_n;
			flux[0] += wifLvn;
			flux[1] += wifLvn * vx;
			flux[2] += wifLvn * vy;
		} else {
			t = kR[0] * vx + kR[1] * vy - log_i0_R;
			wifRvn = wi * rhoR * exp(t) * v_dot_n;
			flux[0] += wifRvn;
			flux[1] += wifRvn * vx;
			flux[2] += wifRvn * vy;
		}
	}
}

void m1_s1_num_flux_kin_simpson3(const float *wL, const float *wR,
																 const float *vn, float *flux)
{
	float kLn, kRn, th;
	float kL[2], kR[2];

	get_kappa(wL, kL, &kLn);
	get_kappa(wR, kR, &kRn);

	const float th_n = get_angle_normalized(vn);
	const float th_kL = get_angle_non_normalized(kL);
	const float th_kR = get_angle_non_normalized(kR);

	const float log_i0_L = get_log_bessel_i0(kLn);
	const float log_i0_R = get_log_bessel_i0(kRn);

	float val_a[3];
	float val_b[3];

	float fluxL[3];
	float fluxR[3];

	/* Simpson's 3 quadrature */
	// const int N = 512;
	// const float h = 0.006135923151543; /* M_PI / 512*/

	const int N = 64;
	const float h = 0.049087385212341; /* M_PI / 64*/

	float th_minL = -M_PI / 2 + th_n;
	float th_maxL = M_PI / 2 + th_n;

	float th_minR = M_PI / 2 + th_n;
	float th_maxR = 3 * M_PI / 2 + th_n;

	th = th_minL;
	flux_integrand(th, kLn, log_i0_L, th_kL, th_n, val_a);
	th = th_maxL;
	flux_integrand(th, kLn, log_i0_L, th_kL, th_n, val_b);

	fluxL[0] = val_a[0] - val_b[0];
	fluxL[1] = val_a[1] - val_b[1];
	fluxL[2] = val_a[2] - val_b[2];

	th = th_minR;
	flux_integrand(th, kRn, log_i0_R, th_kR, th_n, val_a);
	th = th_maxR;
	flux_integrand(th, kRn, log_i0_R, th_kR, th_n, val_b);

	fluxR[0] = val_a[0] - val_b[0];
	fluxR[1] = val_a[1] - val_b[1];
	fluxR[2] = val_a[2] - val_b[2];

	for (int i = 1; i < N - 1; i += 2) {
		th = th_minL + i * h;
		flux_integrand(th, kLn, log_i0_L, th_kL, th_n, val_a);
		th = th_minL + (i + 1) * h;
		flux_integrand(th, kLn, log_i0_L, th_kL, th_n, val_b);
		fluxL[0] += wL[0] * (4.f * val_a[0] + 2.f * val_b[0]);
		fluxL[1] += wL[0] * (4.f * val_a[1] + 2.f * val_b[1]);
		fluxL[2] += wL[0] * (4.f * val_a[2] + 2.f * val_b[2]);

		th = th_minR + i * h;
		flux_integrand(th, kRn, log_i0_R, th_kR, th_n, val_a);
		th = th_minR + (i + 1) * h;
		flux_integrand(th, kRn, log_i0_R, th_kR, th_n, val_b);
		fluxR[0] += wR[0] * (4.f * val_a[0] + 2.f * val_b[0]);
		fluxR[1] += wR[0] * (4.f * val_a[1] + 2.f * val_b[1]);
		fluxR[2] += wR[0] * (4.f * val_a[2] + 2.f * val_b[2]);
	}

	/* Normalize flux */
	flux[0] = h * (fluxL[0] + fluxR[0]) / (6 * M_PI);
	flux[1] = h * (fluxL[1] + fluxR[1]) / (6 * M_PI);
	flux[2] = h * (fluxL[2] + fluxR[2]) / (6 * M_PI);
}

void m1_s1_num_flux_boundary_kin_2(const float *wL, const float *wR,
																	 const float *vn, float *flux)

{
	float kn, th;
	const float rho = wL[0];

	/* Get kappa */
	float kL[2], kR[3];
	get_kappa(wL, kL, &kn);
	get_kappa_specular(kL, vn, kR);

	const float th_n = get_angle_normalized(vn);
	const float th_kL = get_angle_non_normalized(kL);
	const float th_kR = get_angle_non_normalized(kR);

	const float log_i0 = get_log_bessel_i0(kn);

	float val_a[3];
	float val_b[3];

	float fluxL[3];
	float fluxR[3];

	/* Simpson's 3 quadrature */
	// const int N = 512;
	// const float h = 0.006135923151543; /* M_PI / 512*/

	const int N = 64;
	const float h = 0.049087385212341; /* M_PI / 64*/

	float th_minL = -M_PI / 2 + th_n;
	float th_maxL = M_PI / 2 + th_n;

	float th_minR = M_PI / 2 + th_n;
	float th_maxR = -M_PI / 2 + th_n;

	th = th_minL;
	flux_integrand(th, kn, log_i0, th_kL, th_n, val_a);
	th = th_maxL;
	flux_integrand(th, kn, log_i0, th_kL, th_n, val_b);

	fluxL[0] = val_a[0] - val_b[0];
	fluxL[1] = val_a[1] - val_b[1];
	fluxL[2] = val_a[2] - val_b[2];

	th = th_minR;
	flux_integrand(th, kn, log_i0, th_kR, th_n, val_a);
	th = th_maxR;
	flux_integrand(th, kn, log_i0, th_kR, th_n, val_b);

	fluxR[0] = val_a[0] - val_b[0];
	fluxR[1] = val_a[1] - val_b[1];
	fluxR[2] = val_a[2] - val_b[2];

	for (int i = 1; i < N - 1; i += 2) {
		th = th_minL + i * h;
		flux_integrand(th, kn, log_i0, th_kL, th_n, val_a);
		th = th_minL + (i + 1) * h;
		flux_integrand(th, kn, log_i0, th_kL, th_n, val_b);
		fluxL[0] += (4.f * val_a[0] + 2.f * val_b[0]);
		fluxL[1] += (4.f * val_a[1] + 2.f * val_b[1]);
		fluxL[2] += (4.f * val_a[2] + 2.f * val_b[2]);

		th = th_minR + i * h;
		flux_integrand(th, kn, log_i0, th_kR, th_n, val_a);
		th = th_minR + (i + 1) * h;
		flux_integrand(th, kn, log_i0, th_kR, th_n, val_b);
		fluxR[0] += (4.f * val_a[0] + 2.f * val_b[0]);
		fluxR[1] += (4.f * val_a[1] + 2.f * val_b[1]);
		fluxR[2] += (4.f * val_a[2] + 2.f * val_b[2]);
	}

	/* Normlise flux */
	flux[0] = h * rho * (fluxL[0] + ALPHA * fluxR[0]) / (6 * M_PI);
	flux[1] = h * rho * (fluxL[1] + ALPHA * fluxR[1]) / (6 * M_PI);
	flux[2] = h * rho * (fluxL[2] + ALPHA * fluxR[2]) / (6 * M_PI);
}

void m1_s1_num_flux_boundary_kin(const float *wL, const float *wR,
																 const float *vn, float *flux)
{
	float kn, vx, vy, t, v_dot_n, wifLvn, wifRvn;
	float sum_vn, int_fL;

	flux[0] = 0.f;
	flux[1] = 0.f;
	flux[2] = 0.f;

	const float rho = wL[0];
	/* Get kappa left */
	float kL[2];
	get_kappa(wL, kL, &kn);

	/* Get kappa right */
	float kR[2];
	get_kappa_specular(kL, vn, kR);

	const float log_i0 = get_log_bessel_i0(kn);

	__constant const float *vi = quad_circle_uniform_64_vi;
	const float wi = quad_circle_uniform_64_wi;
	const int N = quad_circle_uniform_64_N;

	float flux_test[3] = { 0, 0, 0 };
	float flux_test_2[3] = { 0, 0, 0 };
	m1_s1_num_flux_boundary_kin_2(wL, wR, vn, flux_test);

	/* Specular part */
	sum_vn = 0.f;
	int_fL = 0.f;
	for (unsigned int k = 0; k < N; k++) {
		vx = vi[2 * k];
		vy = vi[2 * k + 1];
		v_dot_n = vx * vn[0] + vy * vn[1];

		if (isgreater(v_dot_n, 0.f)) {
			t = kL[0] * vx + kL[1] * vy - log_i0;
			wifLvn = wi * rho * exp(t) * v_dot_n;

			/* Used in diffusive BC */
			sum_vn += v_dot_n;
			int_fL += wifLvn;

			flux[0] += wifLvn;
			flux[1] += wifLvn * vx;
			flux[2] += wifLvn * vy;

			flux_test_2[0] += wifLvn;
			flux_test_2[1] += wifLvn * vx;
			flux_test_2[2] += wifLvn * vy;

		} else {
			t = kR[0] * vx + kR[1] * vy - log_i0;
			wifRvn = wi * rho * exp(t) * v_dot_n;

			flux[0] += (1.f - BETA) * ALPHA * wifRvn;
			flux[1] += (1.f - BETA) * ALPHA * wifRvn * vx;
			flux[2] += (1.f - BETA) * ALPHA * wifRvn * vy;

			// flux_test_2[0] += (1.f - BETA) * ALPHA * wifRvn;
			// flux_test_2[1] += (1.f - BETA) * ALPHA * wifRvn * vx;
			// flux_test_2[2] += (1.f - BETA) * ALPHA * wifRvn * vy;
		}
	}

	/* Diffusive part */
	// int_fL = int_fL / sum_vn;
	// for (unsigned int k = 0; k < N; k++) {
	// 	vx = vi[2 * k];
	// 	vy = vi[2 * k + 1];
	// 	v_dot_n = vx * vn[0] + vy * vn[1];

	// 	if (isless(v_dot_n, 0.f)) {
	// 		flux[0] += ALPHA * BETA * int_fL * v_dot_n;
	// 		flux[1] += ALPHA * BETA * int_fL * v_dot_n * vx;
	// 		flux[2] += ALPHA * BETA * int_fL * v_dot_n * vy;
	// 	}
	// }
	// printf("flux_0 = %f  \n", flux_test_2[0] - flux_test[0]);
	// printf("flux_1 = %f  \n", flux_test_2[1] - flux_test[1]);
	// printf("flux_2 = %f  \n", flux_test_2[2] - flux_test[2]);

	// flux[0] = flux_test[0];
	// flux[1] = flux_test[1];
	// flux[2] = flux_test[2];
}

void m1_s1_num_flux_kinetic_tanh_sinh(const float *wL, const float *wR,
																			const float *vn, float *flux)
{
	/* Get cell parameters */
	float kLn, kRn, th, kL[2], kR[2];

	get_kappa(wL, kL, &kLn);
	get_kappa(wR, kR, &kRn);

	const float th_n = get_angle_normalized(vn);
	const float th_kL = get_angle_non_normalized(kL);
	const float th_kR = get_angle_non_normalized(kR);

	const float log_i0_L = get_log_bessel_i0(kLn);
	const float log_i0_R = get_log_bessel_i0(kRn);

	float val_a[3], val_b[3], fluxL[3], fluxR[3], IL[3], IR[3];
	/* Integration limits */
	const float th_minL = -M_PI_2 + th_n;
	const float th_maxL = M_PI_2 + th_n;
	const float th_minR = M_PI_2 + th_n;
	const float th_maxR = 3 * M_PI_2 + th_n;

	const float cL = 0.5 * (th_maxL - th_minL);
	const float cR = 0.5 * (th_maxR - th_minR);

	/* Level 0*/
	/* Left cell */
	flux_int(th_minL + cL * tanh_sinh_node[0], kLn, log_i0_L, th_kL, th_n, val_a);
	fluxL[0] = tanh_sinh_weight[0] * val_a[0];
	fluxL[1] = tanh_sinh_weight[0] * val_a[1];
	fluxL[2] = tanh_sinh_weight[0] * val_a[2];

	/* Right cell */
	flux_int(th_minR + cR * tanh_sinh_node[0], kRn, log_i0_R, th_kR, th_n, val_a);
	fluxR[0] = tanh_sinh_weight[0] * val_a[0];
	fluxR[1] = tanh_sinh_weight[0] * val_a[1];
	fluxR[2] = tanh_sinh_weight[0] * val_a[2];

	/* Levels from 1 to 7 */
	for (int k = 1; k < tanh_sinh_max_level; k++) {
		const int id_min = tanh_sinh_offset[k];
		const int id_max = tanh_sinh_offset[k + 1];

		IL[0] = 0;
		IL[1] = 0;
		IL[2] = 0;

		IR[0] = 0;
		IR[1] = 0;
		IR[2] = 0;

		for (int i = id_min; i < id_max; i++) {
			/* Left cell */
			flux_int(th_minL + cL * tanh_sinh_node[i], kLn, log_i0_L, th_kL, th_n,
							 val_a);
			flux_int(th_maxL - cL * tanh_sinh_node[i], kLn, log_i0_L, th_kL, th_n,
							 val_b);

			IL[0] += tanh_sinh_weight[i] * (val_a[0] + val_b[0]);
			IL[1] += tanh_sinh_weight[i] * (val_a[1] + val_b[1]);
			IL[2] += tanh_sinh_weight[i] * (val_a[2] + val_b[2]);

			/* Right cell */
			flux_int(th_minR + cR * tanh_sinh_node[i], kRn, log_i0_R, th_kR, th_n,
							 val_a);
			flux_int(th_maxR - cR * tanh_sinh_node[i], kRn, log_i0_R, th_kR, th_n,
							 val_b);

			IR[0] += tanh_sinh_weight[i] * (val_a[0] + val_b[0]);
			IR[1] += tanh_sinh_weight[i] * (val_a[1] + val_b[1]);
			IR[2] += tanh_sinh_weight[i] * (val_a[2] + val_b[2]);
		}
		/* Left cell */
		fluxL[0] = 0.5 * fluxL[0] + IL[0];
		fluxL[1] = 0.5 * fluxL[1] + IL[1];
		fluxL[2] = 0.5 * fluxL[2] + IL[2];

		/* Right cell */
		fluxR[0] = 0.5 * fluxR[0] + IR[0];
		fluxR[1] = 0.5 * fluxR[1] + IR[1];
		fluxR[2] = 0.5 * fluxR[2] + IR[2];
	}

	flux[0] = (cL * wL[0] * fluxL[0] + cR * wR[0] * fluxR[0]) / (2 * M_PI);
	flux[1] = (cL * wL[0] * fluxL[1] + cR * wR[0] * fluxR[1]) / (2 * M_PI);
	flux[2] = (cL * wL[0] * fluxL[2] + cR * wR[0] * fluxR[2]) / (2 * M_PI);
}

void m1_s1_num_flux_boundary_kinetic_tanh_sinh(const float *wL, const float *wR,
																							 const float *vn, float *flux)
{
	/* Get cell parameters */
	float kn, kL[2], kR[2];

	get_kappa(wL, kL, &kn);
	get_kappa_specular(kL, vn, kR);

	const float th_n = get_angle_normalized(vn);
	const float th_kL = get_angle_non_normalized(kL);
	const float th_kR = get_angle_non_normalized(kR);

	const float log_i0 = get_log_bessel_i0(kn);

		/* TEST */
	// float kL_dot_v_1 = kL[0] * cos(tanh_sinh_node[1]) + kL[1] *  sin(tanh_sinh_node[1]);
	// float kL_dot_v_2 = kn * cos(tanh_sinh_node[1] - th_kL);
	
  // float n_dot_v_1 = vn[0] * cos(tanh_sinh_node[1]) + vn[1] *  sin(tanh_sinh_node[1]);
	// float n_dot_v_2 = cos(tanh_sinh_node[1] - th_n);

	
	
	// printf("k.v = %f, %f, %f\n",kL_dot_v_1, kL_dot_v_2, kL_dot_v_1 - kL_dot_v_2);
	// printf("n/v = %f, %f, %f\n",n_dot_v_1, n_dot_v_2, n_dot_v_1 - n_dot_v_2);






	float val_a[3], val_b[3], fluxL[3], fluxR[3], IL[3], IR[3];
	/* Integration limits */
	const float th_minL = -M_PI_2 + th_n;
	const float th_maxL = M_PI_2 + th_n;
	const float th_minR = M_PI_2 + th_n;
	const float th_maxR = 3 * M_PI_2 + th_n;

	const float cL = 0.5 * (th_maxL - th_minL);
	const float cR = 0.5 * (th_maxR - th_minR);

	/* Level 0*/
	/* Left cell */
	flux_int(th_minL + cL * tanh_sinh_node[0], kn, log_i0, th_kL, th_n, val_a);
	fluxL[0] = tanh_sinh_weight[0] * val_a[0];
	fluxL[1] = tanh_sinh_weight[0] * val_a[1];
	fluxL[2] = tanh_sinh_weight[0] * val_a[2];

	/* Right cell */
	flux_int(th_minR + cR * tanh_sinh_node[0], kn, log_i0, th_kR, th_n, val_a);
	fluxR[0] = tanh_sinh_weight[0] * val_a[0];
	fluxR[1] = tanh_sinh_weight[0] * val_a[1];
	fluxR[2] = tanh_sinh_weight[0] * val_a[2];

	/* Levels from 1 to 7 */
	for (int k = 1; k < tanh_sinh_max_level; k++) {
		const int id_min = tanh_sinh_offset[k];
		const int id_max = tanh_sinh_offset[k + 1];

		IL[0] = 0;
		IL[1] = 0;
		IL[2] = 0;

		IR[0] = 0;
		IR[1] = 0;
		IR[2] = 0;

		for (int i = id_min; i < id_max; i++) {
			/* Left cell */
			flux_int(th_minL + cL * tanh_sinh_node[i], kn, log_i0, th_kL, th_n,
							 val_a);
			flux_int(th_maxL - cL * tanh_sinh_node[i], kn, log_i0, th_kL, th_n,
							 val_b);

			IL[0] += tanh_sinh_weight[i] * (val_a[0] + val_b[0]);
			IL[1] += tanh_sinh_weight[i] * (val_a[1] + val_b[1]);
			IL[2] += tanh_sinh_weight[i] * (val_a[2] + val_b[2]);

			/* Right cell */
			flux_int(th_minR + cR * tanh_sinh_node[i], kn, log_i0, th_kR, th_n,
							 val_a);
			flux_int(th_maxR - cR * tanh_sinh_node[i], kn, log_i0, th_kR, th_n,
							 val_b);

			IR[0] += tanh_sinh_weight[i] * (val_a[0] + val_b[0]);
			IR[1] += tanh_sinh_weight[i] * (val_a[1] + val_b[1]);
			IR[2] += tanh_sinh_weight[i] * (val_a[2] + val_b[2]);
		}
		/* Left cell */
		fluxL[0] = 0.5 * fluxL[0] + IL[0];
		fluxL[1] = 0.5 * fluxL[1] + IL[1];
		fluxL[2] = 0.5 * fluxL[2] + IL[2];

		/* Right cell */
		fluxR[0] = 0.5 * fluxR[0] + IR[0];
		fluxR[1] = 0.5 * fluxR[1] + IR[1];
		fluxR[2] = 0.5 * fluxR[2] + IR[2];
	}

	flux[0] = (cL * wL[0] * fluxL[0] + ALPHA * cR * wL[0] * fluxR[0]) / (2 * M_PI);
	flux[1] = (cL * wL[0] * fluxL[1] + ALPHA * cR * wL[0] * fluxR[1]) / (2 * M_PI);
	flux[2] = (cL * wL[0] * fluxL[2] + ALPHA * cR * wL[0] * fluxR[2]) / (2 * M_PI);
}

void m1_s1_num_flux_fantom(float *wL, const float *vn, float *wR)
{
	float I[2] = { wL[1], wL[2] };
	float Is[2];
	get_kappa_specular(I, vn, Is);

	wR[0] = ALPHA * wL[0];
	wR[1] = ALPHA * (1.f - BETA) * Is[0];
	wR[2] = ALPHA * (1.f - BETA) * Is[1];

	return;
	float kn, vx, vy, tR, tL, wifR, wifL, coef_1, coef_2, coef_3;

	wR[0] = 0;
	wR[1] = 0;
	wR[2] = 0;

	const float rho = wL[0];
	float wL_new[3] = { 0, 0, 0 };
	/* Get kappa left */
	float kL[2];
	get_kappa(wL, kL, &kn);

	// wL[0] = 0;
	// wL[1] = 0;
	// wL[2] = 0;

	/* Get kappa right */
	float kR[2];
	get_kappa_specular(kL, vn, kR);

	const float log_i0 = get_log_bessel_i0(kn);

	__constant const float *vi = quad_circle_uniform_64_vi;
	const float wi = quad_circle_uniform_64_wi;
	const int N = quad_circle_uniform_64_N;

	/* Specular part */
	for (unsigned int k = 0; k < N; k++) {
		vx = vi[2 * k];
		vy = vi[2 * k + 1];

		tR = kR[0] * vx + kR[1] * vy - log_i0;
		tL = kL[0] * vx + kL[1] * vy - log_i0;
		wifR = wi * rho * exp(tR);
		wifL = wi * rho * exp(tL);

		wR[0] += BETA * wifR + (1.f - BETA) * wi * rho;
		wR[1] += BETA * wifR * vx + (1.f - BETA) * wi * rho * vx;
		wR[2] += BETA * wifR * vy + (1.f - BETA) * wi * rho * vy;

		// wL[0] += wifL;
		// wL[1] += wifL * vx;
		// wL[2] += wifL * vy;
	}

	// printf("diff = %.12f \n", fabs(wL[0] - wR[0]));
	// coef_1 = (wL[0] - wL_new[0]) / wL[0];
	// coef_2 = (wL[1] - wL_new[1]) / wL[1];
	// coef_3 = (wL[2] - wL_new[2]) / wL[2];

	// printf("foo %.12f\n", coef);

	/* Correct Alpha */
	// wR[0] *= (ALPHA);
	// wR[1] *= (ALPHA);
	// wR[2] *= (ALPHA);
}

#endif