#include <stdio.h>
#include <stdlib.h>

#include "spherepack.h"
#include "grad.h"
#include "linal.h"

using namespace std;
using namespace linal;

SphereGrad::SphereGrad (const SphereOperator & op) : SphereOperator(op)
{
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
