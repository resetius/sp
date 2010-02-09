#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>

#include "lapl.h"

#ifdef max
#undef max
#endif

using namespace std;

static inline double max(double a, double b)
{
	return (a > b) ? a : b;
}


double rp(double x, double y) {
	return -6.0 * sin(y) * sin(2.0 * x);
}

double ans(double x, double y) {
	return sin(y) * sin(2.0 * x);
}

void solve()
{
	long nlat=19, nlon=36;
	double dlat = M_PI / (nlat-1);
	double dlon = 2. * M_PI /nlon;
	int i, j;

	vector < double > u(nlat * nlon);
	vector < double > v(nlat * nlon);
	vector < double > r(nlat * nlon);

	SphereLaplace lapl(nlat, nlon);

	double nev1 = 0;

	for (j = 0; j < nlon; ++j) {
		for (i = 0; i < nlat; ++i) {
			double theta = -0.5 * M_PI + i * dlat;
			double phi   = j * dlon;

			r[j * nlat + i] = ans(phi, theta);
		}
	}

	lapl.calc(&u[0], &r[0]);

	for (j = 0; j < nlon; ++j) {
		for (i = 0; i < nlat; ++i) {
			double theta = -0.5 * M_PI + i * dlat;
			double phi   = j * dlon;

			nev1 = max(nev1, fabs(u[j * nlat + i] - rp(phi, theta)));
		}
	}

	fprintf(stderr, "nev1=%.16le \n", nev1);
}

int main(int argc, char * argv[])
{
	solve();
}

