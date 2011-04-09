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

void SphereGrad::calc (double * u, double * v, const double * u1)
{
	long ierror = 0;
	long isym   = 0;
	long nt     = 1;
	long n      = nlat * nlon;

	array_t a (mdab * nlat);
	array_t b (mdab * nlat);

	array_t fun (nlat * nlon);
	array_t vv (nlat * nlon);
	array_t ww (nlat * nlon);

	geo2math (&fun[0], &u1[0]);

	shaec_ (&nlat, &nlon, &isym, &nt, &fun[0], &nlat, &nlon, 
		&a[0], &b[0], &mdab, &nlat, &swsave[0],
	        &slsave, &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "shaec_ error\n");
		exit(1);
	}

	gradec_ (&nlat, &nlon, &isym, &nt, &vv[0], &ww[0], &nlat, &nlon,
	         &a[0], &b[0], &mdab, &nlat, &vwsave[0],
	         &vlsave, &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "gradec_ error\n");
		exit(1);
	}

	math2geov(&u[0], &v[0], &ww[0], &vv[0]);
}

void SphereGrad::solve (double * u1, const double * u, const double * v)
{
	long ierror = 0;
	long isym   = 0;
	long ityp   = 0;
	long nt     = 1;
	long n      = nlat * nlon;

	array_t br(mdb * nlat * nt), bi(mdb * nlat * nt);
	array_t cr(mdb * nlat * nt), ci(mdb * nlat * nt);

	array_t grad (nlat * nlon);
	array_t vv (nlat * nlon);
	array_t ww (nlat * nlon);

	geo2mathv(&ww[0], &vv[0], u, v);

	// vhaec dutheta, duphi -> a, b
	vhaec_(&nlat, &nlon, &ityp, &nt, &vv[0], &ww[0], &nlat, &nlon,
		&br[0], &bi[0], &cr[0], &ci[0], &mdb, &nlat, 
		&ivwsave[0], &ivlsave, &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "vhaec_ error\n");
		exit(1);
	}

	igradec_ (&nlat, &nlon, &isym, &nt, &grad[0], &nlat, &nlon,
		&br[0], &bi[0], &mdb, &nlat, 
		&iswsave[0], &islsave, &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "igradec_ error\n");
		exit(1);
	}

	math2geo (&u1[0], &grad[0]);
}
