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

	double * u  = (double*)calloc(nlat * nlon, sizeof(double));
	double * v  = (double*)calloc(nlat * nlon, sizeof(double));
	double * r  = (double*)calloc(nlat * nlon, sizeof(double));
	double nev1 = 0;

	for (i = 0; i < nlat; ++i) {
		for (j = 0; j < nlon; ++j) {
			double theta = -0.5 * M_PI + i * dlat;
			double phi   = j * dlon;

			double x = cos(theta) * cos(phi);
			double y = cos(theta) * sin(phi);
			double z = sin(theta);

			r[i * nlon + j] = -(x * y * (z * z + 6. * (z + 1.)) + z * (z + 2.)) * exp(z);
		}
	}

	lapl.solve(u, r, 1.0, -1.0);


	for (i = 0; i < nlat; ++i) {
		for (j = 0; j < nlon; ++j) {
			double theta = -0.5 * M_PI + i * dlat;
			double phi   = j * dlon;

			double x = cos(theta) * cos(phi);
			double y = cos(theta) * sin(phi);
			double z = sin(theta);

			double ue = (1.+x*y)*exp(z);

			nev1 = max(nev1, fabs(u[i * nlon + j] - ue));
		}
	}


	fprintf(stderr, "nev1=%.16le \n", nev1);

	nev1 = 0;
	lapl.calc(v, u);
	for (i = 0; i < nlat; ++i) {
		for (j = 0; j < nlon; ++j) {
			nev1 = max(nev1, fabs(v[i * nlon + j] - u[i * nlon + j] - r[i * nlon + j]));
		}
	}
	fprintf(stderr, "nev2=%.16le \n", nev1);

	free(u); free(v); free(r);
}

int main(int argc, char * argv[])
{
	solve();
}

