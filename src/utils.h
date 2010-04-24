#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

void transpose (double * out, const double * in, int n, int m);
void transpose1 (double * out, const double * in, double k, int n, int m);
void vec_sum1 (double * out, const double * a, const double * b, double k1, double k2, int size);
void vec_sum (double * out, const double * a, const double * b, int size);
void vec_mult_scalar(double * out, const double * a, double k1, int size);
void fprintfwmatrix(FILE * f, const double * u, long nlat, long nlong, const char * format);
void fprintfwmatrix(const char * fname, const double * u, long nlat, long nlong, const char * format);

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
#ifndef ISINF
#define ISINF
#include <float.h>
#define isnan _isnan
#define snprintf _snprintf
inline bool isinf(double x)
{
	int c = _fpclass(x);
	return (c == _FPCLASS_NINF || c == _FPCLASS_PINF);
}
#endif
#endif

#endif /* UTILS_H */
