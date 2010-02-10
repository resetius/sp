#ifndef UTILS_H
#define UTILS_H

void transpose(double * out, const double * in, int n, int m);
void transpose1(double * out, const double * in, double k, int n, int m);
void vec_sum1(double * out, const double * a, const double * b, double k1, double k2, int size);

#endif /* UTILS_H */
