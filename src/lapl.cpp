#include <stdio.h>
#include <stdlib.h>

#include "lapl.h"
#include "vorticity.h"
#include "linal.h"

using namespace linal;

SphereLaplace::SphereLaplace (const SphereOperator & op) : SphereOperator(op)
{
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

	// находим разложение (a, b) по сферическим гармоникам
	shaec_ (&nlat, &nlon, &isym, &nt, &t[0], &nlat, &nlon, 
	        &a[0], &b[0], &nlat, &nlat, &swsave[0], &slsave,
	        &work[0], &lwork, &ierror);
	// чтобы по разложению (a, b) собрать назад функцию надо воспользоваться
	// функцией shsec_
	if (ierror != 0) {
		fprintf(stderr, "shaec_ error %ld\n", ierror);
		exit(1);
	}

	islapec_ (&nlat, &nlon, &isym, &nt, &koef,
	          &t[0], &nlat, &nlon,
	          &a[0], &b[0], &nlat, &nlat,
	          &iswsave[0], &islsave, &work[0], &lwork, &pertrb, &ierror);
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
	        &swsave[0], &slsave,
	        &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "shaec_ error %ld\n", ierror);
		exit(1);
	}
	slapec_ (&nlat, &nlon, &isym, &nt, 
	          &t[0], &nlat, &nlon,
	          &a[0], &b[0], &nlat, &nlat,
	          &iswsave[0], &islsave, &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "slapec_ error %ld\n", ierror);
		exit(1);
	}

	transpose(out, &t[0], nlon, nlat);
}

void SphereLaplace::make_psi(double * psi, const double * u, const double * v)
{
	SphereVorticity vor (*this);

	vor.calc (&psi[0], &u[0], &v[0]);
	vec_mult_scalar (&psi[0], &psi[0], -1.0, nlat * nlon);
}

