#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "spherepack.h"

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

	long llsave = nlat*(nlat+1)+3*((nlat-2)*(nlat-1)+nlon + 15);
	long lwork  = nlat*(2*nlon+3*(nlat+1)+2*nlat+1);

	double * wshaec = calloc(llsave, sizeof(double));
	double * wshsec = calloc(llsave, sizeof(double));
	double * work   = calloc(lwork, sizeof(double));
	double * dwork  = calloc((nlat + 1), sizeof(double));

	long lshaec = llsave;
	long lshsec = llsave;
	long ldwork = nlat + 1;
	long isym = 0;
	long nt = 1;
	long ierror = 0;

	double xlmbda = 1.0;
	double pertrb = 0.0;
	double * a = calloc(nlat * nlat, sizeof(double));
	double * b = calloc(nlat * nlat, sizeof(double));
	double * u = calloc(nlat * nlon, sizeof(double));
	double * r = calloc(nlat * nlon, sizeof(double));
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

	shaeci_(&nlat, &nlon, wshaec, &lshaec, dwork, &ldwork, &ierror);

	if (ierror != 0) {
		fprintf(stderr, "shaeci error \n");
	}

	shseci_(&nlat, &nlon, wshsec, &lshsec, dwork, &ldwork, &ierror);

	if (ierror != 0) {
		fprintf(stderr, "shseci error \n");
	}

	shaec_(&nlat, &nlon, &isym, &nt, r,
			&nlat, &nlon, a, b,
			&nlat, &nlat, 
			wshaec, &lshaec,
			work, &lwork, &ierror);

	if (ierror != 0) {
		fprintf(stderr, "shaec error \n");
	}

	islapec_(&nlat, &nlon, &isym, &nt, &xlmbda, 
			u, &nlat, &nlon,
			a, b, &nlat, &nlat,
			wshsec, &lshsec, work, &lwork, &pertrb, &ierror);

	if (ierror != 0) {
		fprintf(stderr, "islapec_ error \n");
	}

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


	fprintf(stderr, "nev=%.16le pertrb=%.16le xlmbda=%.16le\n", nev, pertrb, xlmbda);

	free(wshaec); free(wshsec); free(dwork); free(work);
	free(a); free(b); free(u); free(r);
}

int main(int argc, char * argv[])
{
	solve();
}

