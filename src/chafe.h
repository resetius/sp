#ifndef CHAFE_H
#define CHAFE_H
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

#include "lapl.h"

typedef double (*rp_t)(double phi, double lambda, double t, double mu, double sigma);

/**
 * Chafe-Infante equation on sphere.
 * Configuration.
 */
struct SphereChafeConf
{
	long nlat;    ///<latitude
	long nlon;    ///<longitude
	double tau;   ///<time step
	double mu;    ///<\f$\mu\f$
	double sigma; ///<\f$\sigma\f$
	double theta; ///<parameter 0-1. 0.5 for Crank-Nicolson, 1.0 for backward Euler
	rp_t rp;      ///<right part
};

/**
 * Chafe-Infante equation on sphere.
 * \f$\frac{du}{dt} = \mu \delta u - \sigma u + f (u)\f$
 */
class SPAPI SphereChafe {
	typedef std::vector < double > array_t;
	SphereChafeConf conf;
	SphereOperator op;
	SphereLaplace lapl;

public:
	SphereChafe(const SphereChafeConf & conf);

	void calc(double * out, const double * in, double t);
};

#endif /* CHAFE_H */

