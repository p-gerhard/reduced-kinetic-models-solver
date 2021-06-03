#ifndef MATHS_RATIONAL_CL
#define MATHS_RATIONAL_CL

float evaluate_polynomial_simple(__constant const float *p, const int n,
								 const float x)
{
	float res = p[n - 1];
#pragma unroll
	for (int i = n - 2; i >= 0; i--) {
		res = res * x + p[i];
	}
	return res;
}

float evaluate_rational_simple(__constant const float *num,
							   __constant const float *denom, const int count,
							   float z)
{
	float s1, s2;
	if (z <= 1) {
		s1 = num[count - 1];
		s2 = denom[count - 1];
#pragma unroll
		for (int i = count - 2; i >= 0; i--) {
			s1 *= z;
			s2 *= z;
			s1 += num[i];
			s2 += denom[i];
		}
	} else {
		z = 1 / z;
		s1 = num[0];
		s2 = denom[0];
#pragma unroll
		for (unsigned i = 1; i < count; i++) {
			s1 *= z;
			s2 *= z;
			s1 += num[i];
			s2 += denom[i];
		}
	}
	return s1 / s2;
}
#endif