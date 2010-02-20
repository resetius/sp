#include <stdio.h>
#include <stdlib.h>

#include "spherepack.h"
#include "grad.h"
#include "utils.h"

using namespace std;

SphereGrad::SphereGrad (long nlat, long nlon) :
		nlat (nlat), nlon (nlon), mdab ( (nlon + 2) / 2),
		slsave (5*nlat*nlat*nlon), swsave (slsave),
		sldwork (4*nlat*nlat), sdwork (sldwork),
		vlsave (5*nlat*nlat*nlon), vwsave (slsave),
		vldwork (4*nlat*nlat), vdwork (vldwork),
		lwork (5*nlat*nlat*nlon), work (lwork)
{
	long ierror = 0;
	shaeci_ (&nlat, &nlon, &swsave[0], &slsave, &sdwork[0], &sldwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "shaeci_ error\n");
		exit(1);
	}

	vhseci_ (&nlat, &nlon, &vwsave[0], &vlsave, &vdwork[0], &vldwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "vhseci_ error\n");
		exit(1);
	}
}

SphereGrad::~SphereGrad()
{
}

void SphereGrad::calc (double * duph, double * dulambda, const double * u1)
{
	long ierror = 0;
	long isym   = 0;
	long nt     = 1;
	long n      = nlat * nlon;

	array_t a (nlat * nlat);
	array_t b (nlat * nlat);

	array_t u (nlat * nlon);
	array_t dutheta (nlat * nlon);
	array_t duphi (nlat * nlon);

	transpose (&u[0], &u1[0], nlat, nlon);

	shaec_ (&nlat, &nlon, &isym, &nt, &u[0], &nlat, &nlon, 
		&a[0], &b[0], &mdab, &nlat, &swsave[0],
	        &slsave, &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "shaec_ error\n");
		exit(1);
	}

	gradec_ (&nlat, &nlon, &isym, &nt, &dutheta[0], &duphi[0], &nlat, &nlon,
	         &a[0], &b[0], &mdab, &nlat, &vwsave[0],
	         &vlsave, &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "gradec_ error\n");
		exit(1);
	}

	transpose (&duphi[0], &dutheta[0], nlon, nlat);
	transpose (&dulambda[0], &duphi[0], nlon, nlat);
}
