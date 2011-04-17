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

#include <stdlib.h>

#include "div.h"
#include "spherepack.h"
#include "linal.h"

using namespace std;
using namespace linal;

SphereDiv::SphereDiv(const SphereOperator & op): SphereOperator(op)
{
}

SphereDiv::~SphereDiv()
{
}

void SphereDiv::calc(double * div, const double * u, const double *v)
{
	// use divec_ here
	long ierror = 0;
	long isym   = 0;
	long ityp   = 0;
	long nt     = 1;
	long n      = nlat * nlon;

	array_t br(mdc * nlat * nt), bi(mdc * nlat * nt);
	array_t cr(mdc * nlat * nt), ci(mdc * nlat * nt);

	array_t vv(n);
	array_t ww(n);
	array_t vt(n);

	geo2mathv (&ww[0], &vv[0], &u[0], &v[0]);

	vhaec_(&nlat, &nlon, &ityp, &nt, &vv[0], &ww[0], &nlat, &nlon,
		&br[0], &bi[0], &cr[0], &ci[0], &mdc, &nlat,
		&ivwsave[0], &ivlsave, &work[0], &lwork, &ierror);

	if (ierror != 0) {
		fprintf(stderr, "vhaec_ error\n");
		exit(1);
	}

	divec_(&nlat, &nlon, &isym, &nt, &vt[0], &nlat, &nlon,
		&br[0], &bi[0], &mdc, &nlat,
		&iswsave[0], &islsave, &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "divec_ error\n");
		exit(1);
	}

	math2geo (&div[0], &vt[0]);
}
