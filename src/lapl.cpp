#include <stdio.h>
#include <stdlib.h>

#include "lapl.h"
#include "vorticity.h"
#include "utils.h"

SphereLaplace::SphereLaplace (long nlat, long nlon, long isym) :
		nlat (nlat), nlon (nlon), isym(isym),
		lshaec (nlat* (nlat + 1) + 3* ( (nlat - 2) * (nlat - 1) + nlon + 15) ),
		lshsec (lshaec), lwork(nlat*(2*nlon+3*(nlat+1)+2*nlat+1)), ldwork(nlat + 1),
		wshaec(lshaec), wshsec(lshsec), work(lwork), dwork(ldwork)
{
	init();
}

void SphereLaplace::init()
{
	long ierror = 0;
	shaeci_ (&nlat, &nlon, &wshaec[0], &lshaec, &dwork[0], &ldwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "shaeci_ error %ld\n", ierror);
		exit(1);
	}

	shseci_ (&nlat, &nlon, &wshsec[0], &lshsec, &dwork[0], &ldwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "shseci_ error %ld\n", ierror);
		exit(1);
	}
}

SphereLaplace::~SphereLaplace() {}

void SphereLaplace::solve (double * out, const double * in, double mult, double diag)
{
	long ierror = 0;
	long nt   = 1;

	double pertrb = 0;
	double koef   = -diag / mult;

	array_t a (nlat * nlat);
	array_t b (nlat * nlat);
	array_t t (nlat * nlon);

	transpose1(&t[0], in, 1.0 / mult, nlat, nlon);

	shaec_ (&nlat, &nlon, &isym, &nt, &t[0],
	        &nlat, &nlon, &a[0], &b[0],
	        &nlat, &nlat,
	        &wshaec[0], &lshaec,
	        &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "shaec_ error %ld\n", ierror);
		exit(1);
	}

	islapec_ (&nlat, &nlon, &isym, &nt, &koef,
	          &t[0], &nlat, &nlon,
	          &a[0], &b[0], &nlat, &nlat,
	          &wshsec[0], &lshsec, &work[0], &lwork, &pertrb, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "islapec_ error %ld\n", ierror);
		exit(1);
	}

	transpose(out, &t[0], nlon, nlat);
}

void SphereLaplace::calc(double * out, const double * in)
{
	long ierror = 0;
	long nt = 1;
	long isym = 0;

	array_t a (nlat * nlat);
	array_t b (nlat * nlat);
	array_t t (nlat * nlon);

	transpose(&t[0], in, nlat, nlon);

	shaec_ (&nlat, &nlon, &isym, &nt, &t[0],
	        &nlat, &nlon, &a[0], &b[0],
	        &nlat, &nlat,
	        &wshaec[0], &lshaec,
	        &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "shaec_ error %ld\n", ierror);
		exit(1);
	}
	slapec_ (&nlat, &nlon, &isym, &nt, 
	          &t[0], &nlat, &nlon,
	          &a[0], &b[0], &nlat, &nlat,
	          &wshsec[0], &lshsec, &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "slapec_ error %ld\n", ierror);
		exit(1);
	}

	transpose(out, &t[0], nlon, nlat);
}

void SphereLaplace::make_psi(double * psi, const double * u, const double * v)
{
	SphereVorticity vor (nlat, nlon);

	vor.calc (&psi[0], &u[0], &v[0]);
	vec_mult_scalar (&psi[0], &psi[0], -1.0, nlat * nlon);
}

