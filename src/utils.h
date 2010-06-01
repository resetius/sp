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

#endif /* UTILS_H */
