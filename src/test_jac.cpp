#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>

#include "jac.h"

#ifdef max
#undef max
#endif

using namespace std;

static inline double max(double a, double b)
{
	return (a > b) ? a : b;
}

double f1(double x, double y) {
	return cos (x) * sin (y);
}

double f2(double x, double y) {
	return cos (x) * cos (y);
}

double ans(double x, double y) {
	return (cos (x) *cos (y) *cos (x) *cos (y) + cos (x) *sin (y) *cos (x) *sin (y) );
}

void solve()
{
	long nlat=19, nlon=36;
	double dlat = M_PI / (nlat-1);
	double dlon = 2. * M_PI /nlon;
	int i, j;

	vector < double > u(nlat * nlon);
	vector < double > v(nlat * nlon);
	vector < double > a(nlat * nlon);
	vector < double > jac(nlat * nlon);

	SphereJacobian jacobian(nlat, nlon);

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
			fprintf(stderr, "%.16le %.16le %.16le\n", jac[i * nlon + j], a[i * nlon + j], nev1);
		}
	}

	fprintf(stderr, "nev1=%.16le \n", nev1);
}

int main(int argc, char * argv[])
{
	solve();
}

