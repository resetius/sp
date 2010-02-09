#include <math.h>
#include <stdlib.h>
#include "spherepack.h"

void solve()
{
	long nlat=19, nlon=36;
	double dlat = M_PI / (nlat-1);
	double dlon = 2. * M_PI /nlon;
	int i, j;

	long llsave = nlat*(nlat+1)+3*((nlat-2)*(nlat-1)+nlon + 15);

	double * wshaec = malloc(llsave * sizeof(double));
	double *dwork   = malloc((nlat + 1) * sizeof(double));
	long lshaec = llsave;
	long ldwork = nlat + 1;
	long ierror = 0;

	double * r = malloc(nlat * nlon * sizeof(double));

	for (i = 0; i < nlat; ++i) {
		for (j = 0; j < nlon; ++j) {
			double theta = -0.5 * M_PI + (i - 1) * dlat;
			double phi   = (j - 1) * dlon;

			double x = cos(theta) * cos(phi);
			double y = cos(theta) * sin(phi);
			double z = sin(theta);

			r[i * nlon + j] = -(x * x * (z * z + 6. * (z + 1.)) + z * (z + 2.)) * exp(z);
		}
	}

	shaeci_(&nlat, &nlon, wshaec, &lshaec, dwork, &ldwork, &ierror);
}

int main()
{
}

