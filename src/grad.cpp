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

void SphereGrad::calc (double * duph, double * dulambda, const double * u1)
{
	long ierror = 0;
	long isym   = 0;
	long nt     = 1;
	long n      = nlat * nlon;

	array_t a (mdab * nlat);
	array_t b (mdab * nlat);

	array_t u (nlat * nlon);
	array_t dutheta (nlat * nlon);
	array_t duphi (nlat * nlon);

	mat_transpose (&u[0], &u1[0], nlat, nlon);

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

	mat_transpose (&duph[0], &dutheta[0], nlon, nlat);
	mat_transpose (&dulambda[0], &duphi[0], nlon, nlat);
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

	mat_transpose (&dutheta[0], &duph[0], nlat, nlon);
	mat_transpose (&duphi[0], &dulambda[0], nlat, nlon);

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

	mat_transpose (&u1[0], &u[0], nlon, nlat);
}

