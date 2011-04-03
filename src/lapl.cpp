/* Copyright (c) 2010 Alexey Ozeritsky
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>

#include "lapl.h"
#include "vorticity.h"
#include "linal.h"

using namespace linal;
using namespace std;

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

	array_t a (mdab * nlat);
	array_t b (mdab * nlat);
	array_t t (nlat * nlon);

	mat_transpose1(&t[0], in, 1.0 / mult, nlat, nlon);

	// находим разложение (a, b) по сферическим гармоникам
	shaec_ (&nlat, &nlon, &isym, &nt, &t[0], &nlat, &nlon, 
	        &a[0], &b[0], &mdab, &nlat, &swsave[0], &slsave,
	        &work[0], &lwork, &ierror);
	// чтобы по разложению (a, b) собрать назад функцию надо воспользоваться
	// функцией shsec_
	if (ierror != 0) {
		char buf[1024];
		sprintf(buf, "shaec_ error %ld\n", ierror);
		throw runtime_error(buf);
	}

	islapec_ (&nlat, &nlon, &isym, &nt, &koef,
	          &t[0], &nlat, &nlon,
	          &a[0], &b[0], &mdab, &nlat,
	          &iswsave[0], &islsave, &work[0], &lwork, &pertrb, &ierror);
	if (ierror != 0) {
		char buf[1024];
		sprintf(buf, "islapec_ error %ld\n", ierror);
		throw runtime_error(buf);
	}

	mat_transpose(out, &t[0], nlon, nlat);
}

void SphereLaplace::calc(double * out, const double * in)
{
	long ierror = 0;
	long nt = 1;
	long isym = 0;

	array_t a (mdab * nlat);
	array_t b (mdab * nlat);
	array_t t (nlat * nlon);

	mat_transpose(&t[0], in, nlat, nlon);

	shaec_ (&nlat, &nlon, &isym, &nt, &t[0],
	        &nlat, &nlon, &a[0], &b[0],
	        &mdab, &nlat,
	        &swsave[0], &slsave,
	        &work[0], &lwork, &ierror);
	if (ierror != 0) {
		char buf[1024];
		sprintf(buf, "shaec_ error %ld\n", ierror);
		throw runtime_error(buf);
	}
	slapec_ (&nlat, &nlon, &isym, &nt, 
	          &t[0], &nlat, &nlon,
	          &a[0], &b[0], &mdab, &nlat,
	          &iswsave[0], &islsave, &work[0], &lwork, &ierror);
	if (ierror != 0) {
		char buf[1024];
		sprintf(buf, "slapec_ error %ld\n", ierror);
		throw runtime_error(buf);
	}

	mat_transpose(out, &t[0], nlon, nlat);
}

void SphereLaplace::make_psi(double * psi, const double * u, const double * v)
{
	SphereVorticity vor (*this);

	vor.calc (&psi[0], &u[0], &v[0]);
	vec_mult_scalar (&psi[0], &psi[0], -1.0, nlat * nlon);

	// TODO: lapl_1 ???
}

