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
#include <math.h>

#include <list>
#include <algorithm>

#include "linal.h"
#include "spherepack.h"
#include "operator.h"

using namespace std;
using namespace linal;

#undef min

struct SphereOperator::Data
{
	int links;
	list < array_t > arrays;

	Data() : links (1)
	{
	}

	array_t & create (size_t size)
	{
		arrays.push_back (array_t (size) );
		return arrays.back();
	}
};

SphereOperator::SphereOperator (const SphereOperator & op) : SphereNorm (op.nlat, op.nlon),
		d (op.d),
		isym (op.isym), nlat (op.nlat), nlon (op.nlon),

/*     mdab   the first dimension of the arrays a and b as it appears */
/*            in the program that calls shsec. mdab must be at least */
/*            min0(nlat,(nlon+2)/2) if nlon is even or at least */ /* четное */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */ /* нечетное*/

		// scalar
		mdab (min(nlat, (nlon + 2) / 2)), 

		// vector
		mdb (min(nlat, (nlon + 2) / 2)), 
		mdc (min(nlat, (nlon + 2) / 2)),

		slsave (op.slsave), swsave (op.swsave),
		sldwork (op.sldwork), sdwork (op.sdwork),
		vlsave (op.vlsave), vwsave (op.vwsave),
		vldwork (op.vldwork), vdwork (op.vdwork),

		islsave (op.islsave), iswsave (op.iswsave),
		isldwork (op.isldwork), isdwork (op.isdwork),
		ivlsave (op.ivlsave), ivwsave (op.ivwsave),
		ivldwork (op.ivldwork), ivdwork (op.ivdwork),

		lwork (op.lwork), work (op.work)
{
	++d->links;
}

SphereOperator::SphereOperator (long nlat, long nlon, long isym) : SphereNorm (nlat, nlon),
		d (new Data),
		isym (isym), nlat (nlat), nlon (nlon),

		// scalar
		mdab (min(nlat, (nlon + 2) / 2)), 

		// vector
		mdb (min(nlat, (nlon + 2) / 2)), 
		mdc (min(nlat, (nlon + 2) / 2)),

		slsave (5*nlat*nlat*nlon), swsave (d->create (slsave) ),
		sldwork (4*nlat*nlat), sdwork (d->create (sldwork) ),
		vlsave (5*nlat*nlat*nlon), vwsave (d->create (slsave) ),
		vldwork (4*nlat*nlat), vdwork (d->create (vldwork) ),

		islsave (5*nlat*nlat*nlon), iswsave (d->create (islsave) ),
		isldwork (4*nlat*nlat), isdwork (d->create (isldwork) ),
		ivlsave (5*nlat*nlat*nlon), ivwsave (d->create (islsave) ),
		ivldwork (4*nlat*nlat), ivdwork (d->create (ivldwork) ),

		lwork (5*nlat*nlat*nlon), work (d->create (lwork) )
{
	long ierror = 0;
	shaeci_ (&nlat, &nlon, &swsave[0], &slsave, &sdwork[0], &sldwork, &ierror);
	if (ierror != 0)
	{
		fprintf (stderr, "shaeci_ error\n");
		exit (1);
	}

	vhseci_ (&nlat, &nlon, &vwsave[0], &vlsave, &vdwork[0], &vldwork, &ierror);
	if (ierror != 0)
	{
		fprintf (stderr, "vhseci_ error\n");
		exit (1);
	}

	shseci_ (&nlat, &nlon, &iswsave[0], &islsave, &isdwork[0], &isldwork, &ierror);
	if (ierror != 0)
	{
		fprintf (stderr, "shseci_ error\n");
		exit (1);
	}

	vhaeci_ (&nlat, &nlon, &ivwsave[0], &ivlsave, &ivdwork[0], &ivldwork, &ierror);
	if (ierror != 0)
	{
		fprintf (stderr, "vhseci_ error\n");
		exit (1);
	}
}

SphereOperator::~SphereOperator()
{
	--d->links;
	if (!d->links)
	{
		delete d;
	}
}

void SphereOperator::func2koef (double * k, const double * f)
{
	long ierror = 0;
	long nt   = 1;

	double * a = &k[0];
	double * b = &k[mdab * nlat];
	array_t t (nlat * nlon);
	//why ?
	//mat_transpose (&t[0], f, nlat, nlon);
	geo2math(&t[0], f);

	shaec_ (&nlat, &nlon, &isym, &nt, &t[0], &nlat, &nlon,
	        &a[0], &b[0], &mdab, &nlat, &swsave[0], &slsave,
	        &work[0], &lwork, &ierror);

	if (ierror != 0)
	{
		fprintf (stderr, "shaec_ error %ld\n", ierror);
		exit (1);
	}
}

void SphereOperator::koef2func (double * f, const double * k)
{
	long ierror = 0;
	long nt   = 1;

	const double * a = &k[0];
	const double * b = &k[mdab * nlat];
	array_t t (nlat * nlon);

	shsec_ (&nlat, &nlon, &isym, &nt, &t[0], &nlat, &nlon,
	        (double*)&a[0], (double*)&b[0], &mdab, &nlat, &iswsave[0], &islsave,
	        &work[0], &lwork, &ierror);

	if (ierror != 0)
	{
		fprintf (stderr, "shsec_ error %ld\n", ierror);
		exit (1);
	}

	math2geo(f, &t[0]);

	// why?
	// mat_transpose (f, &t[0], nlon, nlat);
}

void SphereOperator::filter(double * out, const double * in)
{
	array_t k(2 * mdab * nlat);
	func2koef(&k[0], in);
	k[0] = 0.0;
	koef2func(out, &k[0]);
}

void SphereOperator::geo2math(double * out, const double * in)
{
	int n = nlat;
	int m = nlon;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			out[j * n + /*i*/ (n - i - 1)] = in[i * m + j];
		}
	}
}

void SphereOperator::math2geo(double * out, const double * in)
{
	int n = nlon;
	int m = nlat;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			out[/*j */ (m - j - 1) * n + i] = in[i * m + j];
		}
	}
}

void SphereOperator::geo2mathv(
	double * dest_w, double * dest_v, 
	const double * source_u, const double * source_v)
{
	geo2math(dest_w, source_u);
	geo2math(dest_v, source_v);
	vec_mult_scalar(dest_v, dest_v, -1.0, nlat * nlon);
}

void SphereOperator::math2geov(
	double * dest_u, double * dest_v, 
	const double * source_w, const double * source_v)
{
	math2geo(dest_u, source_w);
	math2geo(dest_v, source_v);
	vec_mult_scalar(dest_v, dest_v, -1.0, nlat * nlon);
}

void sp_geo2math(double * dest, const double * source, long long * nlat, long long * nlon)
{
	SphereOperator op(*nlat, *nlon, 0);
	op.geo2math(dest, source);
}

void sp_math2geo(double * dest, const double * source, long long * nlat, long long * nlon)
{
	SphereOperator op(*nlat, *nlon, 0);
	op.math2geo(dest, source);
}

void sp_geo2mathv(
	double * dest_w, double * dest_v,
	const double * source_u, const double * source_v,
	long long * nlat, long long * nlon)
{
	SphereOperator op(*nlat, *nlon, 0);
	op.geo2mathv(dest_w, dest_v, source_u, source_v);
}

void sp_math2geov(
	double * dest_u, double * dest_v,
	const double * source_w, const double * source_v,
	long long * nlat, long long * nlon)
{
	SphereOperator op(*nlat, *nlon, 0);
	op.math2geov(dest_u, dest_v, source_w, source_v);
}

void sp_koef2func(double *f, const double *k, long long *nlat, long long *nlon)
{
	SphereOperator op(*nlat, *nlon, 0);
	op.koef2func(f, k);
}

void sp_func2koef(double *k, const double *f, long long *nlat, long long *nlon)
{
	SphereOperator op(*nlat, *nlon, 0);
	op.func2koef(k, f);
}
