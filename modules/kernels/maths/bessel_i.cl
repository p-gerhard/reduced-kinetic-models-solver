#ifndef BESSEL_I_CL
#define BESSEL_I_CL

#pragma once

#include "rational.cl"

__constant static const float bessel_i0_simple_P1[9] = {
	1.00000003928615375e+00f, 2.49999576572179639e-01f, 2.77785268558399407e-02f,
	1.73560257755821695e-03f, 6.96166518788906424e-05f, 1.89645733877137904e-06f,
	4.29455004657565361e-08f, 3.90565476357034480e-10f, 1.48095934745267240e-11f
};

__constant static const float bessel_i0_simple_P2[5] = {
	3.98942651588301770e-01f, 4.98327234176892844e-02f, 2.91866904423115499e-02f,
	1.35614940793742178e-02f, 1.31409251787866793e-01f
};

__constant static const float bessel_i0_simple_P3[3] = {
	3.98942391532752700e-01f, 4.98455950638200020e-02f, 2.94835666900682535e-02f
};

__constant static const double bessel_i0_double_P1[15] = {
	1.00000000000000000e+00, 2.49999999999999909e-01, 2.77777777777782257e-02,
	1.73611111111023792e-03, 6.94444444453352521e-05, 1.92901234513219920e-06,
	3.93675991102510739e-08, 6.15118672704439289e-10, 7.59407002058973446e-12,
	7.59389793369836367e-14, 6.27767773636292611e-16, 4.34709704153272287e-18,
	2.63417742690109154e-20, 1.13943037744822825e-22, 9.07926920085624812e-25
};

__constant static const double bessel_i0_double_P2[22] = {
	3.98942280401425088e-01,	4.98677850604961985e-02,	2.80506233928312623e-02,
	2.92211225166047873e-02,	4.44207299493659561e-02,	1.30970574605856719e-01,
	-3.35052280231727022e+00, 2.33025711583514727e+02,	-1.13366350697172355e+04,
	4.24057674317867331e+05,	-1.23157028595698731e+07, 2.80231938155267516e+08,
	-5.01883999713777929e+09, 7.08029243015109113e+10,	-7.84261082124811106e+11,
	6.76825737854096565e+12,	-4.49034849696138065e+13, 2.24155239966958995e+14,
	-8.13426467865659318e+14, 2.02391097391687777e+15,	-3.08675715295370878e+15,
	2.17587543863819074e+15
};

__constant static const double bessel_i0_double_P3[5] = {
	3.98942280401432905e-01, 4.98677850491434560e-02, 2.80506308916506102e-02,
	2.92179096853915176e-02, 4.53371208762579442e-02
};

float bessel_i0_simple(const float x)
{
	if (x < 7.75) {
		/* Max error in interpolated form: 3.929e-08
		 * Max Error found at float precision = Poly: 1.991226e-07
		 */
		float a = x * x / 4;
		return a * evaluate_polynomial_simple(bessel_i0_simple_P1, 9, a) + 1;
	} else if (x < 50) {
		/* Max error in interpolated form: 5.195e-08
		 * Max Error found at float precision = Poly: 8.502534e-08 
		 */
		return exp(x) *
					 evaluate_polynomial_simple(bessel_i0_simple_P2, 5, (float)(1 / x)) /
					 sqrt(x);
	} else {
		/* Max error in interpolated form: 1.782e-09
		 * Max Error found at float precision = Poly: 6.473568e-08 
		 */
		float ex = exp(x / 2);
		float result =
			ex * evaluate_polynomial_simple(bessel_i0_simple_P3, 3, (float)(1 / x)) /
			sqrt(x);
		result *= ex;
		return result;
	}
}

double bessel_i0_double(const double x)
{
	/* Bessel I0 over[10 ^ -16, 7.75] */
	if (x < 7.75) {
		/* Max error in interpolated form : 3.042e-18 
		 * Max Error found at double precision = Poly : 5.106609e-16 
		 * Cheb : 5.239199e-16 
		 */

		double a = x * x / 4;
		return a * evaluate_polynomial_double(bessel_i0_double_P1, 15, a) + 1;
	} else if (x < 500) {
		/* Max error in interpolated form : 1.685e-16
		 * Max Error found at double precision = Poly : 2.575063e-16 
		 * Cheb : 2.247615e+00
		 */
		return exp(x) *
					 evaluate_polynomial_double(bessel_i0_double_P2, 22,
																			(double)(1 / x)) /
					 sqrt(x);
	} else {
		/* Max error in interpolated form : 2.437e-18
		 * Max Error found at double precision = Poly : 1.216719e-16 
		 */
		double ex = exp(x / 2);
		double result =
			ex * evaluate_polynomial_double(bessel_i0_double_P3, 5, (double)(1 / x)) /
			sqrt(x);
		result *= ex;
		return result;
	}
}

__constant static const float bessel_i1_simple_P1[8] = {
	8.333333221e-02f, 6.944453712e-03f, 3.472097211e-04f, 1.158047174e-05f,
	2.739745142e-07f, 5.135884609e-09f, 5.262251502e-11f, 1.331933703e-12f
};

__constant static const float bessel_i1_simple_P2[5] = {
	3.98942115977513013e-01f, -1.49581264836620262e-01f,
	-4.76475741878486795e-02f, -2.65157315524784407e-02f,
	-1.47148600683672014e-01f
};

__constant static const double bessel_i1_double_P1[13] = {
	8.333333333333333803e-02, 6.944444444444341983e-03, 3.472222222225921045e-04,
	1.157407407354987232e-05, 2.755731926254790268e-07, 4.920949692800671435e-09,
	6.834657311305621830e-11, 7.593969849687574339e-13, 6.904822652741917551e-15,
	5.220157095351373194e-17, 3.410720494727771276e-19, 1.625212890947171108e-21,
	1.332898928162290861e-23
};

__constant static const double bessel_i1_double_P2[22] = {
	3.989422804014406054e-01,	 -1.496033551613111533e-01,
	-4.675104253598537322e-02, -4.090895951581637791e-02,
	-5.719036414430205390e-02, -1.528189554374492735e-01,
	3.458284470977172076e+00,	 -2.426181371595021021e+02,
	1.178785865993440669e+04,	 -4.404655582443487334e+05,
	1.277677779341446497e+07,	 -2.903390398236656519e+08,
	5.192386898222206474e+09,	 -7.313784438967834057e+10,
	8.087824484994859552e+11,	 -6.967602516005787001e+12,
	4.614040809616582764e+13,	 -2.298849639457172489e+14,
	8.325554073334618015e+14,	 -2.067285045778906105e+15,
	3.146401654361325073e+15,	 -2.213318202179221945e+15
};

__constant static const double bessel_i1_double_P3[5] = {
	3.989422804014314820e-01, -1.496033551467584157e-01,
	-4.675105322571775911e-02, -4.090421597376992892e-02,
	-5.843630344778927582e-02
};

float bessel_i1_simple(const float x)
{
	if (x < 7.75) {
		/* Max error in interpolated form : 1.348e-08
     * Max Error found at float precision = Poly : 1.469121e-07
		 */
		float a = x * x / 4;
		float Q[3] = { 1, 0.5f,
									 evaluate_polynomial_simple(bessel_i1_simple_P1, 8, a) };

		return x * evaluate_polynomial_simple_no_constant(Q, 3, a) / 2;
	} else {
		/* Max error in interpolated form: 9.000e-08
		 * Max Error found at float precision = Poly: 1.044345e-07
		 */
		float ex = exp(x / 2);
		float result =
			ex * evaluate_polynomial_simple(bessel_i1_simple_P2, 5, (float)(1 / x)) /
			sqrt(x);
		result *= ex;
		return result;
	}
}

double bessel_i1_double(const double x)
{
	/* Bessel I1 over[10 ^ -16, 7.75] */
	if (x < 7.75) {
		/* Max error in interpolated form: 5.639e-17
		 * Max Error found at double precision = Poly: 1.795559e-16
		 */
		double a = x * x / 4;
		double Q[3] = { 1, 0.5,
										evaluate_polynomial_double(bessel_i1_double_P1, 13, a) };
		return x * evaluate_polynomial_double_no_constant(Q, 3, a) / 2;
	} else if (x < 500) {
		/* Max error in interpolated form: 1.796e-16
		 * Max Error found at double precision = Poly: 2.898731e-16
		 */
		return exp(x) *
					 evaluate_polynomial_double(bessel_i1_double_P2, 22,
																			(double)(1 / x)) /
					 sqrt(x);
	} else {
		/* Max error in interpolated form: 1.320e-19
		 * Max Error found at double precision = Poly: 7.065357e-17
		 */
		double ex = exp(x / 2);
		double result =
			ex * evaluate_polynomial_double(bessel_i1_double_P3, 5, (double)(1 / x)) /
			sqrt(x);
		result *= ex;
		return result;
	}
}

#endif