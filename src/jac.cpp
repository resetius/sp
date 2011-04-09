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
#include "jac.h"
#include "linal.h"

using namespace std;
using namespace linal;

SphereJacobian::SphereJacobian (const SphereOperator & op) : SphereOperator(op), grad(op)
{
}

SphereJacobian::~SphereJacobian()
{
}

void SphereJacobian::calc (double * out, const double * u1, const double * v1)
{
	long n      = nlat * nlon;

	array_t dulambda (n);
	array_t dvlambda (n);

	array_t duphi (n);
	array_t dvphi (n);

	grad.calc(&dulambda[0], &duphi[0], u1);
	grad.calc(&dvlambda[0], &dvphi[0], v1);
	
	for (int i = 0; i < n; ++i)
	{
		out[i] = dulambda[i] * dvphi[i] - duphi[i] * dvlambda[i];
	}
}

void SphereJacobian::calc_t (double * out, const double * u1, const double * v1)
{
	calc(out, u1, v1);
	vec_mult_scalar(out, out, -1, nlon * nlat);
}

