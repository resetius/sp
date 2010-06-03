#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>

#include "jac.h"
#include "linal.h"

#ifdef max
#undef max
#endif

using namespace std;
using namespace linal;

static inline double max(double a, double b)
{
	return (a > b) ? a : b;
}

double f1(double x, double y) {
	double t = 0.0;
	return x*sin(y+t)*ipow(cos(x),4);
}

double f2(double x, double y) {
	return sin(y)*sin(2*x);
}

double ans(double x, double y) {
	double t = 0.0;
	return ipow(cos(x),2)*(2*x*cos(y+t)*cos(x)*sin(y)*cos(2*x)-cos(y)*sin(2*x)*sin(y+t)*cos(x)+4*cos(y)*sin(2*x)*x*sin(y+t)*sin(x));
}

void solve()
{
	long nlat = 5 * 19, nlon = 5 * 36;
	double dlat = M_PI / (nlat-1);
	double dlon = 2. * M_PI /nlon;
	int i, j;

	vector < double > u(nlat * nlon);
	vector < double > v(nlat * nlon);
	vector < double > a(nlat * nlon);
	vector < double > jac(nlat * nlon);

	SphereOperator op(nlat, nlon, 0);
	SphereJacobian jacobian(op);

	double nev1 = 0;

	for (i = 0; i < nlat; ++i) {
		for (j = 0; j < nlon; ++j) {
			double phi    = -0.5 * M_PI + i * dlat;
			double lambda = j * dlon;

			u[i * nlon + j] = f1(phi, lambda);
			v[i * nlon + j] = f2(phi, lambda);
			a[i * nlon + j] = ans(phi, lambda);

		//	fprintf(stderr, "%.16le %.16le\n", u[i * nlon + j], v[i * nlon + j]);
		}
	}

	jacobian.calc(&jac[0], &u[0], &v[0]);

	for (i = 0; i < nlat; ++i) {
		for (j = 0; j < nlon; ++j) {
			nev1 = max(nev1, fabs(jac[i * nlon + j] - a[i * nlon + j]));
			//fprintf(stderr, "%.16le %.16le %.16le\n", jac[i * nlon + j], a[i * nlon + j], nev1);
		}
	}

	fprintf(stderr, "nev1=%.16le \n", nev1);
}

int main(int argc, char * argv[])
{
	solve();
}

