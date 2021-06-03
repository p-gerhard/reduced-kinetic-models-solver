#include "bessel.h"

static float polyeval(const float *p, const int n, const float x)
{
	float res = p[n - 1];

#pragma unroll
	for (int i = n - 2; i > -1; i--) {
		res = res * x + p[i];
	}

	return res;
}


/* 
 *Modified Bessel function of the first kind of order zero and one
 *we use the approximating forms derived in :
 *"Rational Approximations for the Modified Bessel Function of the First Kind 
 *- I0(x) for Computations with Double Precision" by Pavel Holoborodko. 
 *see http://www.advanpix.com/2015/11/11/rational-approximations-for-the- /
 *modified-bessel-function-of-the-first-kind-i0-computations-double-precision
 *This method is used in Boost Library 1.72
*/

float bessel_i0_pavel(const float x)
{
	if (x < 7.75f) {
		const float P[9] = { 1.00000003928615375e+00f, 2.49999576572179639e-01f,
												 2.77785268558399407e-02f, 1.73560257755821695e-03f,
												 6.96166518788906424e-05f, 1.89645733877137904e-06f,
												 4.29455004657565361e-08f, 3.90565476357034480e-10f,
												 1.48095934745267240e-11f };

		const float a = x * x / 4.f;
		return a * polyeval(P, 9, a) + 1.f;

	} else if (x < 50.f) {
		const float P[5] = { 3.98942651588301770e-01f, 4.98327234176892844e-02f,
												 2.91866904423115499e-02f, 1.35614940793742178e-02f,
												 1.31409251787866793e-01f };

		return (exp(x) * polyeval(P, 5, (1.f / x)) / sqrt(x));

	} else {
		const float P[3] = { 3.98942391532752700e-01f, 4.98455950638200020e-02f,
												 2.94835666900682535e-02f };

		const float ex = exp(x / 2.f);
		float result = (ex * polyeval(P, 3, (1.f / x)) / sqrt(x));
		result *= ex;
		return result;
	}
}

float bessel_i1_pavel(const float x)
{
	if (x < 7.75f) {
		const float P[8] = { 8.333333221e-02f, 6.944453712e-03f, 3.472097211e-04f,
												 1.158047174e-05f, 2.739745142e-07f, 5.135884609e-09f,
												 5.262251502e-11f, 1.331933703e-12f };

		const float a = x * x / 4.f;
		const float Q[3] = { 1, 0.5f, polyeval(P, 8, a) };
		return x * polyeval(Q, 3, a) / 2.f;
	} else {
		const float P[5] = { 3.98942115977513013e-01f, -1.49581264836620262e-01f,
												 -4.76475741878486795e-02f, -2.65157315524784407e-02f,
												 -1.47148600683672014e-01f };

		const float ex = exp(x / 2.f);
		float res = (ex * polyeval(P, 5, (1.f / x)) / sqrt(x));
		res *= ex;
		return res;
	}
}

/* 
 *Modified Bessel function of the first kind of order zero and one
 *we use the approximating forms derived in :
 *Abramowitz & Stegun, Handbook of mathematical functions
*/

float bessel_i0_abramowitz(const float x)
{
	float res;
	float y;
	const float ax = fabs(x);

	if (ax < 3.75f) {
		y = x / 3.75f;
		y *= y;
		res =
			1.0f +
			y * (3.5156229f +
					 y * (3.0899424f +
								y * (1.2067492f +
										 y * (0.2659732f + y * (0.360768e-1f + y * 0.45813e-2f)))));
	} else {
		y = 3.75f / ax;
		res =
			(exp(ax) / sqrt(ax)) *
			(0.39894228f +
			 y *
				 (0.1328592e-1f +
					y * (0.225319e-2f +
							 y * (-0.157565e-2f +
										y * (0.916281e-2f + y * (-0.2057706e-1f +
																						 y * (0.2635537e-1f +
																									y * (-0.1647633e-1f +
																											 y * 0.392377e-2f))))))));
	}
	return res;
}

float bessel_i1_abramowitz(const float x)
{
	float res;
	float y;
	const float ax = fabs(x);

	if (ax < 3.75f) {
		y = x / 3.75f;
		y *= y;
		res = ax *
					(0.5f + y * (0.87890594f +
											 y * (0.51498869f +
														y * (0.15084934f +
																 y * (0.2658733e-1f +
																			y * (0.301532e-2f + y * 0.32411e-3f))))));
	} else {
		y = 3.75f / ax;

		res = 0.2282967e-1f +
					y * (-0.2895312e-1f + y * (0.1787654e-1f - y * 0.420059e-2f));

		res = 0.39894228f +
					y * (-0.3988024e-1f +
							 y * (-0.362018e-2f +
										y * (0.163801e-2f + y * (-0.1031555e-1f + y * res))));

		res *= (exp(ax) / sqrt(ax));
	}
	return x < 0.f ? -res : res;
}