#ifndef M1_S2_CL
#define M1_S2_CL

#define M1_MAX_KAPPA_SIMPLE 65.F
#define M1_TOL_R 1E-8F

#ifdef USE_QUAD_LEBEDEV
#undef USE_QUAD_TRIGAUSS
#define USE_QUAD_LEBEDEV_NB_V_302
#include <maths/quad_lebedev.cl>
#endif

#ifdef USE_QUAD_TRIGAUSS
#undef USE_QUAD_LEBEDEV
#include <maths/quad_trigauss.cl>
#endif

#include <maths/rational.cl>

__constant static const float chi_p[7] = { 0.333333397667f, -1.312060956673f,
										   2.271654209855f, -2.514343153802f,
										   2.116995591933f, -1.216531569506f,
										   0.325269105049f };

__constant static const float chi_q[7] = { 1.000000000000f,	 -3.936171252072f,
										   5.614892442165f,	 -2.820874826409f,
										   -0.780269390244f, 1.253825802778f,
										   -0.327086152123f };

static void safe_normalise(float u[3], float u_normalised[3], float *u_norm)
{
	*u_norm = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);

	if ((*u_norm) > FLT_EPSILON) {
		u_normalised[0] = u[0] / (*u_norm);
		u_normalised[1] = u[1] / (*u_norm);
		u_normalised[2] = u[2] / (*u_norm);
	} else {
		u_normalised[0] = 0.f;
		u_normalised[1] = 0.f;
		u_normalised[2] = 0.f;
	}
}

static void safe_invert(float a, float *inv_a)
{
	if (fabs(a) > FLT_EPSILON) {
		*inv_a = 1.f / a;
	} else {
		*inv_a = 1.f / FLT_EPSILON;
	}
}

static float get_chi_interp(const float r)
{
	return evaluate_rational_simple(chi_p, chi_q, 7, r);
}

static float get_r(const float rho, const float In)
{
	float r;
	float inv_rho;

	safe_invert(rho, &inv_rho);

	r = fmin(In * inv_rho, 0.97f);

	return r;
}

static void get_kappa(const float wn[4], float k[3], float *k_norm)
{
	float I[3] = { wn[1], wn[2], wn[3] };
	float In[3];
	float intensity_norm;

	/* Normalized intensity vector In */
	safe_normalise(I, In, &intensity_norm);

	const float r = get_r(wn[0], intensity_norm);
	const float chi = get_chi_interp(r);

	*k_norm = fmin(2.f * r / (1.f - chi), M1_MAX_KAPPA_SIMPLE);

	k[0] = *k_norm * In[0];
	k[1] = *k_norm * In[1];
	k[2] = *k_norm * In[2];
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
	float I[3] = { wn[1], wn[2], wn[3] };
	float In[3];
	float intensity_norm;

	/* Normalized intensity vector In */
	safe_normalise(I, In, &intensity_norm);

	/* Get closure value */
	const float r = get_r(wn[0], intensity_norm);
	const float chi = get_chi_interp(r);

	const float t1 = 0.5f * (1.f - chi);
	const float t2 = 0.5f * (3.f * chi - 1.f);

	/* Diagonal terms */
	const float Pxx = wn[0] * (t1 + t2 * In[0] * In[0]);
	const float Pyy = wn[0] * (t1 + t2 * In[1] * In[1]);
	const float Pzz = wn[0] * (t1 + t2 * In[2] * In[2]);

	/* Off diagonal terms : Pxy = Pyz, Pxz = Pzx, Pyz = Pzy */
	const float Pxy = wn[0] * t2 * In[0] * In[1];
	const float Pxz = wn[0] * t2 * In[0] * In[2];
	const float Pyz = wn[0] * t2 * In[1] * In[2];

	flux[0] = wn[1] * vn[0] + wn[2] * vn[1] + wn[3] * vn[2];
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

static void m1_num_flux_boundary(const float wL[4], const float wR[4],
								 const float vn[3], float flux[4])
{
	float kLn, kRn;
	float kL[3];
	float kR[3];

	get_kappa(wL, kL, &kLn);
	get_kappa_specular(kL, vn, kR);

	// const float log_sinh_L = get_log_sinh(kLn);
	// const float log_sinh_R = get_log_sinh(kRn);

	float a = get_a(wL[0], kLn);

	float sum_vn = 0.f;
	float int_fL = 0.f;

	flux[0] = 0.f;
	flux[1] = 0.f;
	flux[2] = 0.f;
	flux[3] = 0.f;

	for (int iv = 0; iv < QUAD_NB_V; iv++) {
		float vx = quad_vi[3 * iv + 0];
		float vy = quad_vi[3 * iv + 1];
		float vz = quad_vi[3 * iv + 2];
#ifdef USE_QUAD_TRIGAUSS
		float wi = quad_wi[iv] / (4 * M_PI);
#else
		float wi = quad_wi[iv];
#endif
		float v_dot_n = vx * vn[0] + vy * vn[1] + vz * vn[2];

		if (v_dot_n > 0.f) {
			float tL = kL[0] * vx + kL[1] * vy + kL[2] * vz;
			float wifLvn = wi * a * exp(tL) * v_dot_n;
			sum_vn += v_dot_n * wi;
			int_fL += wifLvn;

			flux[0] += wifLvn;
			flux[1] += wifLvn * vx;
			flux[2] += wifLvn * vy;
			flux[3] += wifLvn * vz;
		} else {
			float tR = kR[0] * vx + kR[1] * vy + kR[2] * vz;
			float wifRvn = (1.f - BETA) * wi * a * exp(tR) * v_dot_n;
			flux[0] += wifRvn;
			flux[1] += wifRvn * vx;
			flux[2] += wifRvn * vy;
			flux[3] += wifRvn * vz;
		}
	}

	/* Diffusive part */
	int_fL = int_fL / sum_vn;
	for (int iv = 0; iv < QUAD_NB_V; iv++) {
		float vx = quad_vi[3 * iv + 0];
		float vy = quad_vi[3 * iv + 1];
		float vz = quad_vi[3 * iv + 2];
#ifdef USE_QUAD_TRIGAUSS
		float wi = quad_wi[iv] / (4 * M_PI);
#else
		float wi = quad_wi[iv];
#endif
		float v_dot_n = vx * vn[0] + vy * vn[1] + vz * vn[2];

		if (v_dot_n < 0.f) {
			float wifRvn = wi * ALPHA * BETA * int_fL * v_dot_n;
			flux[0] += wifRvn;
			flux[1] += wifRvn * vx;
			flux[2] += wifRvn * vy;
			flux[3] += wifRvn * vz;
		}
	}
}

#ifdef USE_KINETIC_NUM_FLUX
static void m1_num_flux_kinetic(const float wL[4], const float wR[4],
								const float vn[3], float flux[4])
{
	float kLn, kRn;
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

	for (int iv = 0; iv < QUAD_NB_V; iv++) {
		float vx = quad_vi[3 * iv + 0];
		float vy = quad_vi[3 * iv + 1];
		float vz = quad_vi[3 * iv + 2];
#ifdef USE_QUAD_TRIGAUSS
		float wi = quad_wi[iv] / (4 * M_PI);
#else
		float wi = quad_wi[iv];
#endif
		float v_dot_n = vx * vn[0] + vy * vn[1] + vz * vn[2];

		if (v_dot_n > 0.f) {
			float tL = kL[0] * vx + kL[1] * vy + kL[2] * vz;
			float wifLvn = wi * aL * exp(tL) * v_dot_n;

			flux[0] += wifLvn;
			flux[1] += wifLvn * vx;
			flux[2] += wifLvn * vy;
			flux[3] += wifLvn * vz;
		} else {
			float tR = kR[0] * vx + kR[1] * vy + kR[2] * vz;
			float wifRvn = wi * aR * exp(tR) * v_dot_n;

			flux[0] += wifRvn;
			flux[1] += wifRvn * vx;
			flux[2] += wifRvn * vy;
			flux[3] += wifRvn * vz;
		}
	}
}
#endif

#ifdef USE_KINETIC_NUM_FLUX_SPLIT
void static m1_num_flux_RL(const float wL[4], const float wR[4],
						   const float vn[3], float flux[4])
{
	float kLn, kRn;
	float kL[3];
	float kR[3];

	get_kappa(wL, kL, &kLn);
	get_kappa(wR, kR, &kRn);

	float aL = get_a(wL[0], kLn);
	float aR = get_a(wR[0], kRn);

	float flux_o[4] = { 0.f, 0.f, 0.f, 0.f };
	float flux_i[4] = { 0.f, 0.f, 0.f, 0.f };

	/* R->L interface */

	for (int iv = 0; iv < QUAD_NB_V; iv++) {
		/* Integration on outer half-sphere */
		float vx = quad_vi_R[3 * iv + 0];
		float vy = quad_vi_R[3 * iv + 1];
		float vz = quad_vi_R[3 * iv + 2];
		float wi = quad_wi_half[iv] / (4.f * M_PI);

		float v_dot_n = vn[0] * vx;
		float tL = kL[0] * vx + kL[1] * vy + kL[2] * vz;
		float wifLvn = aL * wi * exp(tL) * v_dot_n;

		flux_o[0] += wifLvn;
		flux_o[1] += wifLvn * vx;
		flux_o[2] += wifLvn * vy;
		flux_o[3] += wifLvn * vz;

		/* Integration on inner half-sphere : (flip x,z coordinates) */
		vx = -vx;
		vz = -vz;
		v_dot_n = -v_dot_n;

		float tR = kR[0] * vx + kR[1] * vy + kR[2] * vz;
		float wifRvn = aR * wi * exp(tR) * v_dot_n;

		flux_i[0] += wifRvn;
		flux_i[1] += wifRvn * vx;
		flux_i[2] += wifRvn * vy;
		flux_i[3] += wifRvn * vz;
	}

	flux[0] = flux_o[0] + flux_i[0];
	flux[1] = flux_o[1] + flux_i[1];
	flux[2] = flux_o[2] + flux_i[2];
	flux[3] = flux_o[3] + flux_i[3];
}

void static m1_num_flux_FB(const float wL[4], const float wR[4],
						   const float vn[3], float flux[4])
{
	float kLn, kRn;
	float kL[3];
	float kR[3];

	get_kappa(wL, kL, &kLn);
	get_kappa(wR, kR, &kRn);

	float aL = get_a(wL[0], kLn);
	float aR = get_a(wR[0], kRn);

	float flux_o[4] = { 0.f, 0.f, 0.f, 0.f };
	float flux_i[4] = { 0.f, 0.f, 0.f, 0.f };

	/* R->L interface */
	for (int iv = 0; iv < QUAD_NB_V; iv++) {
		/* Integration on outer half-sphere */
		float vx = quad_vi_F[3 * iv + 0];
		float vy = quad_vi_F[3 * iv + 1];
		float vz = quad_vi_F[3 * iv + 2];
		float wi = quad_wi_half[iv] / (4.f * M_PI);

		float v_dot_n = vn[1] * vy;
		float tL = kL[0] * vx + kL[1] * vy + kL[2] * vz;
		float wifLvn = aL * wi * exp(tL) * v_dot_n;

		flux_o[0] += wifLvn;
		flux_o[1] += wifLvn * vx;
		flux_o[2] += wifLvn * vy;
		flux_o[3] += wifLvn * vz;

		/* Integration on inner half-sphere : (flip y,z coordinates) */
		vy = -vy;
		vz = -vz;
		v_dot_n = -v_dot_n;

		float tR = kR[0] * vx + kR[1] * vy + kR[2] * vz;
		float wifRvn = aR * wi * exp(tR) * v_dot_n;

		flux_i[0] += wifRvn;
		flux_i[1] += wifRvn * vx;
		flux_i[2] += wifRvn * vy;
		flux_i[3] += wifRvn * vz;
	}

	flux[0] = flux_o[0] + flux_i[0];
	flux[1] = flux_o[1] + flux_i[1];
	flux[2] = flux_o[2] + flux_i[2];
	flux[3] = flux_o[3] + flux_i[3];
}

void static m1_num_flux_NS(const float wL[4], const float wR[4],
						   const float vn[3], float flux[4])
{
	float kLn, kRn;
	float kL[3];
	float kR[3];

	get_kappa(wL, kL, &kLn);
	get_kappa(wR, kR, &kRn);

	float aL = get_a(wL[0], kLn);
	float aR = get_a(wR[0], kRn);

	float flux_o[4] = { 0.f, 0.f, 0.f, 0.f };
	float flux_i[4] = { 0.f, 0.f, 0.f, 0.f };

	/* R->L interface */
	for (int iv = 0; iv < QUAD_NB_V; iv++) {
		/* Integration on outer half-sphere */
		float vx = quad_vi_N[3 * iv + 0];
		float vy = quad_vi_N[3 * iv + 1];
		float vz = quad_vi_N[3 * iv + 2];
		float wi = quad_wi_half[iv] / (4.f * M_PI);

		float v_dot_n = vn[2] * vz;
		float tL = kL[0] * vx + kL[1] * vy + kL[2] * vz;
		float wifLvn = aL * wi * exp(tL) * v_dot_n;

		flux_o[0] += wifLvn;
		flux_o[1] += wifLvn * vx;
		flux_o[2] += wifLvn * vy;
		flux_o[3] += wifLvn * vz;

		/* Integration on inner half-sphere : (flip x,y,z coordinates) */
		vx = -vx;
		vy = -vy;
		vz = -vz;
		v_dot_n = -v_dot_n;

		float tR = kR[0] * vx + kR[1] * vy + kR[2] * vz;
		float wifRvn = aR * wi * exp(tR) * v_dot_n;

		flux_i[0] += wifRvn;
		flux_i[1] += wifRvn * vx;
		flux_i[2] += wifRvn * vy;
		flux_i[3] += wifRvn * vz;
	}

	flux[0] = flux_o[0] + flux_i[0];
	flux[1] = flux_o[1] + flux_i[1];
	flux[2] = flux_o[2] + flux_i[2];
	flux[3] = flux_o[3] + flux_i[3];
}
#endif
#endif