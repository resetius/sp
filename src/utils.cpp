#include <stdio.h>
#include <stdlib.h>

#include "utils.h"

void transpose(double * out, const double * in, int n, int m)
{
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			out[j * n + i] = in[i * m + j];
		}
	}
}

void transpose1(double * out, const double * in, double k, int n, int m)
{
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			out[j * n + i] = k * in[i * m + j];
		}
	}
}

void vec_sum1(double * out, const double * a, const double * b, double k1, double k2, int size)
{
	for (int i = 0; i < size; ++i) {
		out[i] = k1 * a[i] + k2 * b[i];
	}
}

void vec_mult_scalar(double * out, const double * a, double k1, int size)
{
	for (int i = 0; i < size; ++i) {
		out[i] = k1 * a[i];
	}
}
