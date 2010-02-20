#include <stdio.h>
#include <stdlib.h>

#include "spherepack.h"
#include "grad.h"
#include "utils.h"

using namespace std;

SphereGrad::SphereGrad (long nlat, long nlon) :
		nlat (nlat), nlon (nlon), mdab ( (nlon + 2) / 2), mdb((nlon+1)/2),
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

	transpose (&duph[0], &dutheta[0], nlon, nlat);
	transpose (&dulambda[0], &duphi[0], nlon, nlat);
}

void SphereGrad::solve (double * u1, const double * duph, const double * dulambda)
{
	long ierror = 0;
	long isym   = 0;
	long ityp   = 0;
	long nt     = 1;
	long n      = nlat * nlon;

	array_t br(mdb * nlat * nt), bi(mdb * nlat * nt);
	array_t cr(mdb * nlat * nt), ci(mdb * nlat * nt);

	array_t u (nlat * nlon);
	array_t dutheta (nlat * nlon);
	array_t duphi (nlat * nlon);

	transpose (&dutheta[0], &duph[0], nlat, nlon);
	transpose (&duphi[0], &dulambda[0], nlat, nlon);

	// vhaec dutheta, duphi -> a, b
	vhaec_(&nlat, &nlon, &ityp, &nt, &dutheta[0], &duphi[0], &nlat, &nlon,
		&br[0], &bi[0], &cr[0], &ci[0], &mdb, &nlat, 
		&ivwsave[0], &ivlsave, &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "vhaec_ error\n");
		exit(1);
	}

	igradec_ (&nlat, &nlon, &isym, &nt, &u[0], &nlat, &nlon,
		&br[0], &bi[0], &mdb, &nlat, 
		&iswsave[0], &islsave, &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "igradec_ error\n");
		exit(1);
	}

	transpose (&u1[0], &u[0], nlon, nlat);
}
