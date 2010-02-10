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

#endif /* UTILS_H */
