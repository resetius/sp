#ifndef UTILS_H
#define UTILS_H

void transpose (double * out, const double * in, int n, int m);
void transpose1 (double * out, const double * in, double k, int n, int m);
void vec_sum1 (double * out, const double * a, const double * b, double k1, double k2, int size);

inline double
ipow (double x, int p)
{
	int i;
	double r = 1;
	for (i = 0; i < p; i++)
	{
		r *= x;
	}
	return r;
}

#ifdef WIN32
#include <float.h>
#define isnan _isnan
inline bool isinf(double x)
{
	int c = _fpclass(x);
	return (c == _FPCLASS_NINF || c == _FPCLASS_PINF);
}
#endif

#endif /* UTILS_H */
