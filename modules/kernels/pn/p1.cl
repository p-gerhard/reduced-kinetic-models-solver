#ifndef P1_CL
#define P1_CL

#ifdef IS_2D
/* P1 2D MODEL :
 *
 * Type 2 filter lookup datas :
 * - fpn_type_2_filter_den          = N / (N + 1),
 * - fpn_type_2_filter_num_tab[iw]  = l[iw] / (N + 1),
 * 
 * with l[iw] the associated Pn degree 'l' for a given conservative 
 * variable 'iw'.
 */
#ifdef USE_FPN_TYPE_2
__constant const float fpn_type_2_filter_den = 0.500000000f;
__constant const float fpn_type_2_filter_num_tab[3] = { 0.000000000f,
														0.500000000f,
														0.500000000f };
#endif

void pn_num_flux_rusanov(const float wL[3], const float wR[3],
						 const float vn[2], float flux[3])
{
	const float lambda = 1;
	const float x0 = 0.28867513F * vn[0];
	const float x1 = 0.28867513F * vn[1];
	const float x2 = 0.5F * lambda;
	const float x3 = wL[0] + wR[0];
	flux[0] =
		x0 * (wL[2] + wR[2]) + x1 * (wL[1] + wR[1]) - x2 * (-wL[0] + wR[0]);
	flux[1] = x1 * x3 - x2 * (-wL[1] + wR[1]);
	flux[2] = x0 * x3 - x2 * (-wL[2] + wR[2]);
}

void pn_num_flux_kinetic_R(const float wL[3], const float wR[3], float flux[3])
{
	flux[0] = 0.25000006F * wL[0] + 0.288675129F * wL[2] - 0.25000006F * wR[0] +
			  0.288675129F * wR[2];
	flux[1] = 0.187500015F * wL[1] - 0.187500015F * wR[1];
	flux[2] = 0.288675129F * wL[0] + 0.375F * wL[2] + 0.288675129F * wR[0] -
			  0.375F * wR[2];
}

void pn_num_flux_kinetic_L(const float wL[3], const float wR[3], float flux[3])
{
	flux[0] = 0.25000006F * wL[0] - 0.288675129F * wL[2] - 0.25000006F * wR[0] -
			  0.288675129F * wR[2];
	flux[1] = 0.187500015F * wL[1] - 0.187500015F * wR[1];
	flux[2] = -0.288675129F * wL[0] + 0.375F * wL[2] - 0.288675129F * wR[0] -
			  0.375F * wR[2];
}

void pn_num_flux_kinetic_N(const float wL[3], const float wR[3], float flux[3])
{
	flux[0] = 0.25000006F * wL[0] + 0.288675129F * wL[1] - 0.25000006F * wR[0] +
			  0.288675129F * wR[1];
	flux[1] = 0.288675129F * wL[0] + 0.375F * wL[1] + 0.288675129F * wR[0] -
			  0.375F * wR[1];
	flux[2] = 0.1875F * wL[2] - 0.187500015F * wR[2];
}

void pn_num_flux_kinetic_S(const float wL[3], const float wR[3], float flux[3])
{
	flux[0] = 0.25000006F * wL[0] - 0.288675129F * wL[1] - 0.25000006F * wR[0] -
			  0.288675129F * wR[1];
	flux[1] = -0.288675129F * wL[0] + 0.375F * wL[1] - 0.288675129F * wR[0] -
			  0.375F * wR[1];
	flux[2] = 0.187500015F * wL[2] - 0.1875F * wR[2];
}

void pn_num_flux_boundary_kinetic_R(const float wL[3], float flux[3])
{
	const float alpha = ALPHA;
	const float beta = BETA;
	const float x0 = 0.25000006F * wL[0];
	const float x1 = alpha * wL[2];
	const float x2 = 0.288675129F * wL[0];
	flux[0] = -alpha * x0 + 0.288675129F * wL[2] + x0 +
			  x1 * (-5.96046448e-8F * beta - 0.288675129F);
	flux[1] = alpha * wL[1] * (0.187500015F * beta - 0.187500015F) +
			  0.187500015F * wL[1];
	flux[2] =
		alpha * x2 + 0.375F * wL[2] + x1 * (0.375F - 0.0416666269F * beta) + x2;
}

void pn_num_flux_boundary_kinetic_L(const float wL[3], float flux[3])
{
	const float alpha = ALPHA;
	const float beta = BETA;
	const float x0 = 0.25000006F * wL[0];
	const float x1 = alpha * wL[2];
	const float x2 = 0.288675129F * wL[0];
	flux[0] = -alpha * x0 - 0.288675129F * wL[2] + x0 +
			  x1 * (2.98023224e-8F * beta + 0.288675159F);
	flux[1] = alpha * wL[1] * (0.187500015F * beta - 0.187500015F) +
			  0.187500015F * wL[1];
	flux[2] = -alpha * x2 + 0.375F * wL[2] +
			  x1 * (0.375F - 0.0416666269F * beta) - x2;
}

void pn_num_flux_boundary_kinetic_N(const float wL[3], float flux[3])
{
	const float alpha = ALPHA;
	const float beta = BETA;
	const float x0 = 0.25000006F * wL[0];
	const float x1 = alpha * wL[1];
	const float x2 = 0.288675129F * wL[0];
	flux[0] = -alpha * x0 + 0.288675129F * wL[1] + x0 +
			  x1 * (-5.96046448e-8F * beta - 0.288675129F);
	flux[1] = alpha * x2 + 0.375F * wL[1] +
			  x1 * (0.37500003F - 0.0416666567F * beta) + x2;
	flux[2] =
		alpha * wL[2] * (0.187500015F * beta - 0.187500015F) + 0.1875F * wL[2];
}

void pn_num_flux_boundary_kinetic_S(const float wL[3], float flux[3])
{
	const float alpha = ALPHA;
	const float beta = BETA;
	const float x0 = 0.25000006F * wL[0];
	const float x1 = alpha * wL[1];
	const float x2 = 0.288675129F * wL[0];
	flux[0] = -alpha * x0 - 0.288675129F * wL[1] + x0 +
			  x1 * (5.96046448e-8F * beta + 0.288675129F);
	flux[1] = -alpha * x2 + 0.375F * wL[1] +
			  x1 * (0.375F - 0.0416666269F * beta) - x2;
	flux[2] = alpha * wL[2] * (0.187500015F * beta - 0.187500015F) +
			  0.187500015F * wL[2];
}
#else
/* P1 3D MODEL :
 *
 * Type 2 filter lookup datas :
 * - fpn_type_2_filter_den          = N / (N + 1),
 * - fpn_type_2_filter_num_tab[iw]  = l[iw] / (N + 1),
 * 
 * with l[iw] the associated Pn degree 'l' for a given conservative 
 * variable 'iw'.
 */
#ifdef USE_FPN_TYPE_2
__constant const float fpn_type_2_filter_den = 0.500000000f;
__constant const float fpn_type_2_filter_num_tab[4] = {
	0.000000000f, 0.500000000f, 0.500000000f, 0.500000000f
};
#endif
void pn_num_flux_rusanov(const float wL[4], const float wR[4],
						 const float vn[3], float flux[4])
{
	const float lambda = 1;
	const float x0 = 0.28867513F * vn[0];
	const float x1 = 0.28867513F * vn[1];
	const float x2 = 0.288675135F * vn[2];
	const float x3 = 0.5F * lambda;
	const float x4 = wL[0] + wR[0];
	flux[0] = x0 * (wL[3] + wR[3]) + x1 * (wL[1] + wR[1]) +
			  x2 * (wL[2] + wR[2]) - x3 * (-wL[0] + wR[0]);
	flux[1] = x1 * x4 - x3 * (-wL[1] + wR[1]);
	flux[2] = x2 * x4 - x3 * (-wL[2] + wR[2]);
	flux[3] = x0 * x4 - x3 * (-wL[3] + wR[3]);
}

void pn_num_flux_kinetic_R(const float wL[4], const float wR[4], float flux[4])
{
	flux[0] = 0.25000006F * wL[0] + 0.288675129F * wL[3] - 0.25000006F * wR[0] +
			  0.288675129F * wR[3];
	flux[1] = 0.187500015F * wL[1] - 0.187500015F * wR[1];
	flux[2] = 0.18749997F * wL[2] - 0.1875F * wR[2];
	flux[3] = 0.288675129F * wL[0] + 0.375F * wL[3] + 0.288675129F * wR[0] -
			  0.375F * wR[3];
}

void pn_num_flux_kinetic_L(const float wL[4], const float wR[4], float flux[4])
{
	flux[0] = 0.25000006F * wL[0] - 0.288675129F * wL[3] - 0.25000006F * wR[0] -
			  0.288675129F * wR[3];
	flux[1] = 0.187500015F * wL[1] - 0.187500015F * wR[1];
	flux[2] = 0.1875F * wL[2] - 0.18749997F * wR[2];
	flux[3] = -0.288675129F * wL[0] + 0.375F * wL[3] - 0.288675129F * wR[0] -
			  0.375F * wR[3];
}

void pn_num_flux_kinetic_F(const float wL[4], const float wR[4], float flux[4])
{
	flux[0] = 0.25000006F * wL[0] + 0.288675129F * wL[1] - 0.25000006F * wR[0] +
			  0.288675129F * wR[1];
	flux[1] = 0.288675129F * wL[0] + 0.375F * wL[1] + 0.288675129F * wR[0] -
			  0.375F * wR[1];
	flux[2] = 0.1875F * wL[2] - 0.18749997F * wR[2];
	flux[3] = 0.1875F * wL[3] - 0.187500015F * wR[3];
}

void pn_num_flux_kinetic_B(const float wL[4], const float wR[4], float flux[4])
{
	flux[0] = 0.25000006F * wL[0] - 0.288675129F * wL[1] - 0.25000006F * wR[0] -
			  0.288675129F * wR[1];
	flux[1] = -0.288675129F * wL[0] + 0.375F * wL[1] - 0.288675129F * wR[0] -
			  0.375F * wR[1];
	flux[2] = 0.18749997F * wL[2] - 0.1875F * wR[2];
	flux[3] = 0.187500015F * wL[3] - 0.1875F * wR[3];
}

void pn_num_flux_kinetic_N(const float wL[4], const float wR[4], float flux[4])
{
	flux[0] = 0.25000006F * wL[0] + 0.288675129F * wL[2] - 0.25000006F * wR[0] +
			  0.288675129F * wR[2];
	flux[1] = 0.1875F * wL[1] - 0.1875F * wR[1];
	flux[2] = 0.288675129F * wL[0] + 0.37499997F * wL[2] +
			  0.288675129F * wR[0] - 0.375F * wR[2];
	flux[3] = 0.1875F * wL[3] - 0.1875F * wR[3];
}

void pn_num_flux_kinetic_S(const float wL[4], const float wR[4], float flux[4])
{
	flux[0] = 0.25000006F * wL[0] - 0.288675129F * wL[2] - 0.25000006F * wR[0] -
			  0.288675129F * wR[2];
	flux[1] = 0.1875F * wL[1] - 0.1875F * wR[1];
	flux[2] = -0.288675129F * wL[0] + 0.375F * wL[2] - 0.288675129F * wR[0] -
			  0.37499997F * wR[2];
	flux[3] = 0.1875F * wL[3] - 0.1875F * wR[3];
}

void pn_num_flux_boundary_kinetic_R(const float wL[4], float flux[4])
{
	const float alpha = ALPHA;
	const float beta = BETA;
	const float x0 = 0.25000006F * wL[0];
	const float x1 = alpha * wL[3];
	const float x2 = 0.288675129F * wL[0];
	flux[0] = -alpha * x0 + 0.288675129F * wL[3] + x0 +
			  x1 * (-5.96046448e-8F * beta - 0.288675129F);
	flux[1] = alpha * wL[1] * (0.187500015F * beta - 0.187500015F) +
			  0.187500015F * wL[1];
	flux[2] = alpha * wL[2] * (0.1875F * beta - 0.1875F) + 0.18749997F * wL[2];
	flux[3] =
		alpha * x2 + 0.375F * wL[3] + x1 * (0.375F - 0.0416666269F * beta) + x2;
}

void pn_num_flux_boundary_kinetic_L(const float wL[4], float flux[4])
{
	const float alpha = ALPHA;
	const float beta = BETA;
	const float x0 = 0.25000006F * wL[0];
	const float x1 = alpha * wL[3];
	const float x2 = 0.288675129F * wL[0];
	flux[0] = -alpha * x0 - 0.288675129F * wL[3] + x0 +
			  x1 * (2.98023224e-8F * beta + 0.288675159F);
	flux[1] = alpha * wL[1] * (0.187500015F * beta - 0.187500015F) +
			  0.187500015F * wL[1];
	flux[2] =
		alpha * wL[2] * (0.18749997F * beta - 0.18749997F) + 0.1875F * wL[2];
	flux[3] = -alpha * x2 + 0.375F * wL[3] +
			  x1 * (0.375F - 0.0416666269F * beta) - x2;
}

void pn_num_flux_boundary_kinetic_F(const float wL[4], float flux[4])
{
	const float alpha = ALPHA;
	const float beta = BETA;
	const float x0 = 0.25000006F * wL[0];
	const float x1 = alpha * wL[1];
	const float x2 = 0.288675129F * wL[0];
	flux[0] = -alpha * x0 + 0.288675129F * wL[1] + x0 +
			  x1 * (-5.96046448e-8F * beta - 0.288675129F);
	flux[1] = alpha * x2 + 0.375F * wL[1] +
			  x1 * (0.37500003F - 0.0416666567F * beta) + x2;
	flux[2] =
		alpha * wL[2] * (0.18749997F * beta - 0.18749997F) + 0.1875F * wL[2];
	flux[3] =
		alpha * wL[3] * (0.187500015F * beta - 0.187500015F) + 0.1875F * wL[3];
}

void pn_num_flux_boundary_kinetic_B(const float wL[4], float flux[4])
{
	const float alpha = ALPHA;
	const float beta = BETA;
	const float x0 = 0.25000006F * wL[0];
	const float x1 = alpha * wL[1];
	const float x2 = 0.288675129F * wL[0];
	flux[0] = -alpha * x0 - 0.288675129F * wL[1] + x0 +
			  x1 * (5.96046448e-8F * beta + 0.288675129F);
	flux[1] = -alpha * x2 + 0.375F * wL[1] +
			  x1 * (0.375F - 0.0416666269F * beta) - x2;
	flux[2] = alpha * wL[2] * (0.1875F * beta - 0.1875F) + 0.18749997F * wL[2];
	flux[3] = alpha * wL[3] * (0.187500015F * beta - 0.187500015F) +
			  0.187500015F * wL[3];
}

void pn_num_flux_boundary_kinetic_N(const float wL[4], float flux[4])
{
	const float alpha = ALPHA;
	const float beta = BETA;
	const float x0 = 0.25000006F * wL[0];
	const float x1 = alpha * wL[2];
	const float x2 = alpha * (0.1875F * beta - 0.1875F);
	const float x3 = 0.288675129F * wL[0];
	flux[0] = -alpha * x0 + 0.288675129F * wL[2] + x0 +
			  x1 * (-5.96046448e-8F * beta - 0.288675129F);
	flux[1] = wL[1] * x2 + 0.1875F * wL[1];
	flux[2] = alpha * x3 + 0.37499997F * wL[2] +
			  x1 * (0.375F - 0.0416666269F * beta) + x3;
	flux[3] = wL[3] * x2 + 0.1875F * wL[3];
}

void pn_num_flux_boundary_kinetic_S(const float wL[4], float flux[4])
{
	const float alpha = ALPHA;
	const float beta = BETA;
	const float x0 = 0.25000006F * wL[0];
	const float x1 = alpha * wL[2];
	const float x2 = alpha * (0.1875F * beta - 0.1875F);
	const float x3 = 0.288675129F * wL[0];
	flux[0] = -alpha * x0 - 0.288675129F * wL[2] + x0 +
			  x1 * (5.96046448e-8F * beta + 0.288675129F);
	flux[1] = wL[1] * x2 + 0.1875F * wL[1];
	flux[2] = -alpha * x3 + 0.375F * wL[2] +
			  x1 * (0.375F - 0.0416666269F * beta) - x3;
	flux[3] = wL[3] * x2 + 0.1875F * wL[3];
}
#endif
#endif