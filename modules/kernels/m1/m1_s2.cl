#ifndef M1_S2_CL
#define M1_S2_CL

#pragma once

#define M1_MAX_KAPPA_SIMPLE 40.F
#define M1_TOL_R 1E-8F

#include <maths/rational.cl>

#define USE_QUAD_LEBEDEV_ORDER_7
#include <maths/quad_lebedev.cl>

__constant static const float chi_p[7] = { 0.333333397667f, -1.312060956673f,
										   2.271654209855f, -2.514343153802f,
										   2.116995591933f, -1.216531569506f,
										   0.325269105049f };

__constant static const float chi_q[7] = { 1.000000000000f,	 -3.936171252072f,
										   5.614892442165f,	 -2.820874826409f,
										   -0.780269390244f, 1.253825802778f,
										   -0.327086152123f };

static float get_r(const float rho, const float In)
{
	float r;

	if (In < rho) {
		r = In / rho;

	} else {
		printf("WARNING RATIO \n");
		r = 0.99f;
	}

	return r;
}

// static float get_r(const float rho, const float In)
// {
// 	float r;
// 	float tol = 1e-6;

// 	if (isless(rho, 0.f)) {
// 		printf("ERROR R, %2.12f, %2.12f %1.12f \n", rho, r);
// 	}

// 	if (isless(rho, tol)) {
// 		r = 0;
// 	} else {
// 		// r = fmax(0.f, fmin(In / rho, 1.f - tol));
// 		r = In / rho;
// 		// r = In / rho;
// 		r = fmin(In / rho, 0.999f);
// 	}

// 	if (isgreaterequal(r, 1.f)) {
// 		printf("ERROR I, %2.12f, %2.12f \n", rho, r);
// 	}
// 	return r;
// }

static float get_chi_interp(const float r)
{
	return evaluate_rational_simple(chi_p, chi_q, 7, r);
}

static void get_kappa(const float wn[4], float k[3], float *k_norm)
{
	float3 I = (float3)(wn[1], wn[2], wn[3]);
	/* Normalized intensity vector In */
	float3 In = normalize(I);
	const float intensity_norm = length(I);

	const float r = get_r(wn[0], intensity_norm);
	const float chi = get_chi_interp(r);

	*k_norm = fmin(2.f * r / (1.f - chi), M1_MAX_KAPPA_SIMPLE);

	k[0] = *k_norm * In.x;
	k[1] = *k_norm * In.y;
	k[2] = *k_norm * In.z;
}

static void get_kappa_specular(const float kL[3], const float vn[3],
							   float kR[3])
{
	/* Construct kR according to the specular law */
	const float kL_dot_vn = kL[0] * vn[0] + kL[1] * vn[1] + kL[2] * vn[2];
	kR[0] = kL[0] - 2.f * kL_dot_vn * vn[0];
	kR[1] = kL[1] - 2.f * kL_dot_vn * vn[1];
	kR[2] = kL[2] - 2.f * kL_dot_vn * vn[2];
}

static float get_a(const float rho, const float k_norm)
{
	float a;
	float tol = 1e-6;
	if (k_norm < tol) {
		float t1 = k_norm * k_norm;
		a = rho * (0.1e1 - (0.1e1 / 0.6e1) * t1 + (0.7e1 / 360e0) * t1 * t1);
	} else {
		a = rho * k_norm / (sinh(k_norm));
	}
	return a;
}

static void phy_flux(const float wn[4], const float vn[3], float flux[4])
{
	const float rho = wn[0];
	const float Ix = wn[1];
	const float Iy = wn[2];
	const float Iz = wn[3];

	float3 I = (float3)(Ix, Iy, Iz);

	/* Normalized intensity vector In */
	float3 In = normalize(I);
	const float intensity_norm = length(I);

	/* Get closure value */
	const float r = get_r(rho, intensity_norm);
	const float chi = get_chi_interp(r);

	const float t1 = 0.5f * (1.f - chi);
	const float t2 = 0.5f * (3.f * chi - 1.f);

	/* Diagonal terms */
	const float Pxx = rho * (t1 + t2 * In.x * In.x);
	const float Pyy = rho * (t1 + t2 * In.y * In.y);
	const float Pzz = rho * (t1 + t2 * In.z * In.z);

	/* Off diagonal terms : Pxy = Pyz, Pxz = Pzx, Pyz = Pzy */
	const float Pxy = rho * t2 * In.x * In.y;
	const float Pxz = rho * t2 * In.x * In.z;
	const float Pyz = rho * t2 * In.y * In.z;

	flux[0] = Ix * vn[0] + Iy * vn[1] + Iz * vn[2];
	flux[1] = Pxx * vn[0] + Pxy * vn[1] + Pxz * vn[2];
	flux[2] = Pxy * vn[0] + Pyy * vn[1] + Pyz * vn[2];
	flux[3] = Pxz * vn[0] + Pyz * vn[1] + Pzz * vn[2];
}

void m1_num_flux_rusanov(const float wL[4], const float wR[4],
						 const float vn[3], float flux[4])
{
	float fL[4];
	float fR[4];
	const float lambda = 1.f;

	phy_flux(wL, vn, fL);
	phy_flux(wR, vn, fR);

	flux[0] = 0.5f * ((fL[0] + fR[0]) - lambda * (wR[0] - wL[0]));
	flux[1] = 0.5f * ((fL[1] + fR[1]) - lambda * (wR[1] - wL[1]));
	flux[2] = 0.5f * ((fL[2] + fR[2]) - lambda * (wR[2] - wL[2]));
	flux[3] = 0.5f * ((fL[3] + fR[3]) - lambda * (wR[3] - wL[3]));
}

static void m1_num_flux_kinetic_lebedev(const float wL[4], const float wR[4],
										const float vn[3], float flux[4])
{
	float wi, vx, vy, vz, t, v_dot_n;
	float kLn, kRn, wifLvn, wifRvn;

	float kL[3];
	float kR[3];


	get_kappa(wL, kL, &kLn);
	get_kappa(wR, kR, &kRn);

	// const float log_sinh_L = get_log_sinh(kLn);
	// const float log_sinh_R = get_log_sinh(kRn);

	float aL = get_a(wL[0], kLn);
	float aR = get_a(wR[0], kRn);

	flux[0] = 0.f;
	flux[1] = 0.f;
	flux[2] = 0.f;
	flux[3] = 0.f;
	
	for (int iv = 0; iv < LEBEDEV_ORDER_N; iv++) {
		vx = quad_leb_vi[3 * iv + 0];
		vy = quad_leb_vi[3 * iv + 1];
		vz = quad_leb_vi[3 * iv + 2];
		wi = quad_leb_wi[iv];

		v_dot_n = vx * vn[0] + vy * vn[1] + vz * vn[2];

		if (v_dot_n > 0.f) {
			t = kL[0] * vx + kL[1] * vy + kL[2] * vz;
			wifLvn = wi * aL * exp(t) * v_dot_n;
			flux[0] += wifLvn;
			flux[1] += wifLvn * vx;
			flux[2] += wifLvn * vy;
			flux[3] += wifLvn * vz;
		} else {
			t = kR[0] * vx + kR[1] * vy + kR[2] * vz;
			wifRvn = wi * aR * exp(t) * v_dot_n;
			flux[0] += wifRvn;
			flux[1] += wifRvn * vx;
			flux[2] += wifRvn * vy;
			flux[3] += wifRvn * vz;
		}
	}
}

#endif