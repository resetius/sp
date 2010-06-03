#include <stdio.h>
#include <stdlib.h>

#include "spherepack.h"
#include "operator.h"

using namespace std;

SphereOperator::SphereOperator (long nlat, long nlon) :
		nlat (nlat), nlon (nlon), 
		mdab ( (nlon + 2) / 2), mdb((nlon+1)/2),
		slsave (5*nlat*nlat*nlon), swsave (slsave),
		sldwork (4*nlat*nlat), sdwork (sldwork),
		vlsave (5*nlat*nlat*nlon), vwsave (slsave),
		vldwork (4*nlat*nlat), vdwork (vldwork),

		islsave (5*nlat*nlat*nlon), iswsave (islsave),
		isldwork (4*nlat*nlat), isdwork (isldwork),
		ivlsave (5*nlat*nlat*nlon), ivwsave (islsave),
		ivldwork (4*nlat*nlat), ivdwork (ivldwork),

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

	shseci_ (&nlat, &nlon, &iswsave[0], &islsave, &isdwork[0], &isldwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "shseci_ error\n");
		exit(1);
	}

	vhaeci_ (&nlat, &nlon, &ivwsave[0], &ivlsave, &ivdwork[0], &ivldwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "vhseci_ error\n");
		exit(1);
	}
}

SphereOperator::~SphereOperator()
{
}

