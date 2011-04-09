/* Copyright (c) 2010-2011 Alexey Ozeritsky
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

#include <math.h>
#include <stdlib.h>

#include "vorticity.h"
#include "spherepack.h"
#include "linal.h"

using namespace std;
using namespace linal;

SphereVorticity::SphereVorticity(const SphereOperator & op): SphereOperator(op)
{
}

SphereVorticity::~SphereVorticity()
{
}

void SphereVorticity::calc(double * div, const double * u, const double *v)
{
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

	math2geov(&ww[0], &vv[0], &u[0], &v[0]);

	vhaec_(&nlat, &nlon, &ityp, &nt, &vv[0], &ww[0], &nlat, &nlon,
		&br[0], &bi[0], &cr[0], &ci[0], &mdc, &nlat,
		&ivwsave[0], &ivlsave, &work[0], &lwork, &ierror);

	if (ierror != 0) {
		fprintf(stderr, "vhaec_ error\n");
		exit(1);
	}

	vrtec_(&nlat, &nlon, &isym, &nt, &vt[0], &nlat, &nlon,
		&cr[0], &ci[0], &mdc, &nlat,
		&iswsave[0], &islsave, &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "vrtec_ error\n");
		exit(1);
	}

	geo2math (&div[0], &vt[0]);
}

bool SphereVorticity::test()
{
	long n      = nlat * nlon;
	int i, j;
	array_t u(n);
	array_t v(n);
	array_t psi(n);
	array_t psi1(n);
	double nr = 1e10;


	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;

	for (i = 0; i < nlat; ++i)
	{
		double phi    = -0.5 * M_PI + i * dlat;

		for (j = 0; j < nlon; ++j)
		{
			double lambda = j * dlon;
			u[i * nlon + j]   =  (M_PI / 2 + phi) * (M_PI / 2. - phi);
			psi[i * nlon + j] =
				- sin(phi) * (M_PI / 2 + phi) * (M_PI / 2. - phi)
				+ cos(phi) * (M_PI / 2 - phi) - cos(phi) * (M_PI / 2 + phi);
			psi[i * nlon + j] /= cos(phi);
		}
	}

	calc(&psi1[0], &u[0], &v[0]);

	nr = dist(&psi1[0], &psi[0]);
	fprintf(stderr, "vorticity test= %.16lf\n", nr);
	return nr < 1e-7;
}
