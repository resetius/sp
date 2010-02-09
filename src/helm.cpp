#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "lapl.h"

#ifdef max
#undef max
#endif

static inline double max(double a, double b)
{
	return (a > b) ? a : b;
}

void solve()
{
	long nlat=19, nlon=36;
	double dlat = M_PI / (nlat-1);
	double dlon = 2. * M_PI /nlon;
	int i, j;

	SphereLaplace lapl(nlat, nlon);

	double * u = (double*)calloc(nlat * nlon, sizeof(double));
	double * r = (double*)calloc(nlat * nlon, sizeof(double));
	double nev = 0;

	for (j = 0; j < nlon; ++j) {
		for (i = 0; i < nlat; ++i) {
			double theta = -0.5 * M_PI + i * dlat;
			double phi   = j * dlon;

			double x = cos(theta) * cos(phi);
			double y = cos(theta) * sin(phi);
			double z = sin(theta);

			r[j * nlat + i] = -(x * y * (z * z + 6. * (z + 1.)) + z * (z + 2.)) * exp(z);
		}
	}

	lapl.solve(u, r, 1.0);


	for (j = 0; j < nlon; ++j) {
		for (i = 0; i < nlat; ++i) {
			double theta = -0.5 * M_PI + i * dlat;
			double phi   = j * dlon;

			double x = cos(theta) * cos(phi);
			double y = cos(theta) * sin(phi);
			double z = sin(theta);

			double ue = (1.+x*y)*exp(z);

			nev = max(nev, fabs(u[j * nlat + i] - ue));
		}
	}


	fprintf(stderr, "nev=%.16le \n", nev);

	free(u); free(r);
}

int main(int argc, char * argv[])
{
	solve();
}

