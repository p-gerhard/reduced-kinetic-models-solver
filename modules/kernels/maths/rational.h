#ifndef MATHS_RATIONAL_CL
#define MATHS_RATIONAL_CL


float evaluate_polynomial_simple(__constant const float *p, const int n,
										const float x);

float evaluate_rational_simple(__constant const float *num,
							   __constant const float *denom, const int count,
							   float z);

#endif





