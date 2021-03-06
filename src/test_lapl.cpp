#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>

#include "lapl.h"
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


static double rp(double x, double y) {
	return -6.0 * sin(y) * sin(2.0 * x);
}

static double ans(double x, double y) {
	return sin(y) * sin(2.0 * x);
}

static double ans2(double x, double y, double t)
{
	return x*sin(y+t)*ipow(cos(x),4);
}

static double rp2(double x, double y, double t, double mu, double sigma)
{
	return sin(y+t)*ipow(cos(x),2)*(9*mu*sin(x)*cos(x)+20*mu*x*ipow(cos(x),2)+sigma*x*ipow(cos(x),2)-15*mu*x);
}

static bool solve()
{
	long nlat = 3*19, nlon = 3*36;
	double dlat = M_PI / (nlat-1);
	double dlon = 2. * M_PI /nlon;
	int i, j;
	bool ret = true;

	double mu    = -0.5;
	double sigma = 1000;

	vector < double > u(nlat * nlon);
	vector < double > v(nlat * nlon);
	vector < double > r(nlat * nlon);

	SphereOperator op(nlat, nlon, 0);
	SphereLaplace lapl(op);

	double nev1 = 0;

	for (i = 0; i < nlat; ++i) {
		for (j = 0; j < nlon; ++j) {
			double phi    = -0.5 * M_PI + i * dlat;
			double lambda = j * dlon;

			r[i * nlon + j] = rp2(phi, lambda, 0, mu, sigma);
		}
	}

	lapl.solve(&u[0], &r[0], -mu, sigma);

	for (i = 0; i < nlat; ++i) {
		for (j = 0; j < nlon; ++j) {
			double phi    = -0.5 * M_PI + i * dlat;
			double lambda = j * dlon;

			nev1 = max(nev1, fabs(u[i * nlon + j] - ans2(phi, lambda, 0)));
		}
	}

	fprintf(stderr, "nev1=%.16le \n", nev1);
	ret = nev1 < 1e-5;

	lapl.calc(&v[0], &u[0]);
	nev1 = 0.0;

	for (i = 0; i < nlat; ++i) {
		for (j = 0; j < nlon; ++j) {
			double phi    = -0.5 * M_PI + i * dlat;
			double lambda = j * dlon;

			nev1 = max(nev1, fabs(5 * v[i * nlon + j] - r[i * nlon + j]));
		}
	}

	fprintf(stderr, "nev1=%.16le \n", nev1);
	return ret;
}

extern "C" int test_lapl(int argc, char * argv[])
{
	if (solve()) {
		return 0;
	}
	return -1;
}
