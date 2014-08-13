/* Copyright (c) 2010-2014 Alexey Ozeritsky
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
#include "norm.h"
#include "linal.h"

using namespace linal;

SphereNorm::SphereNorm(long nlat, long nlon): 
	nlat(nlat), nlon(nlon), cosi(new double[nlat])
{
	double dlat = M_PI / (nlat - 1);
	for (int i = 0; i < nlat; ++i)
	{
		cosi[i] = cos(-0.5 * M_PI + i * dlat);
	}
}

SphereNorm::~SphereNorm()
{
	delete [] cosi;
}

double SphereNorm::scalar(const double * u, const double * v)
{
	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;

	double sum = 0.0;
	for (int i = 0; i < nlat; ++i) {
		for (int j = 0; j < nlon; ++j) {
			sum += cosi[i] * u[i * nlon + j] * v[i * nlon + j];
		}
	}
	return sum * dlat * dlon;
}

double SphereNorm::norm(const double * u)
{
	return sqrt(scalar(u, u));
}

double SphereNorm::dist(const double * u, const double * v)
{
	long n = nlat * nlon;
	std::vector < double > tmp(n);
	vec_sum1(&tmp[0], u, v, 1.0, -1.0, n);
	return norm(&tmp[0]);
}
