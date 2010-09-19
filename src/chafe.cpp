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

#include <math.h>

#include "chafe.h"
#include "linal.h"

using namespace linal;

SphereChafe::SphereChafe(const SphereChafeConf & conf)
	: conf(conf), op(conf.nlat, conf.nlon, 0), lapl(op)
{
}

void SphereChafe::calc(double * out, const double * u, double t)
{
	double tau   = conf.tau;
	double mu    = conf.mu;
	double sigma = conf.sigma;
	double theta = conf.theta;
	long nlat    = conf.nlat;
	long nlon    = conf.nlon;
	long n       = nlat * nlon;

	double dlat = M_PI / (nlat-1);
	double dlon = 2. * M_PI /nlon;

	array_t delta_u(n);
	array_t rp(n);

	lapl.calc(&delta_u[0], &u[0]);
	// u/dt + (1-\theta) mu \Delta u - (1-\theta) \sigma u
	vec_sum1(&rp[0], &u[0], &delta_u[0], 1.0 / tau - (1-theta) * sigma, mu * (1-theta), n);

	// u/dt + (1-\theta) mu \Delta u / 2 - (1-\theta) \sigma u / 2 + f
	for (int i = 0; i < nlat; i++) {
		double phi    = -0.5 * M_PI + i * dlat;
		for (int j = 0; j < nlon; ++j) {
			double lambda = j * dlon;
			if (conf.rp) {
				rp[i * nlon + j] += conf.rp(phi, lambda, t, mu, sigma);
			}
		}
	}

	if (theta == 0) {
		vec_mult_scalar(out, &rp[0], tau, n);
	} else {
		// u'/dt - \theta mu \Delta u' + \theta \sigma u'
		lapl.solve(out, &rp[0], -theta * mu, 1.0 / tau + theta * sigma);
	}
}

