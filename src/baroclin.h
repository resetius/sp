#ifndef BARVORTEX_H
#define BARVORTEX_H
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

#include <vector>

#include "lapl.h"
#include "jac.h"
#include "norm.h"
#include "solver.h"

class SphereBaroclin: public SphereNorm
{
public:
	typedef std::vector < double > array_t;

	struct Conf
	{
		typedef double (*rp_t ) (double phi, double lambda, 
			double t, const Conf * conf );

		long nlat;
		long nlon;

		double tau;
		double sigma;
		double sigma1;
		double mu;
		double mu1;
		double k1;
		double k2;
		double theta;
		double alpha;

		rp_t cor;
		rp_t rp1;
		rp_t rp2;

		double * cor2;
		double * rp12;
		double * rp22;

		Conf(): cor(0), rp1(0), rp2(0), cor2(0), rp12(0), rp22(0) {}
	};

private:
	Conf conf;
	SphereOperator op;
	SphereLaplace lapl;
	SphereJacobian jac;
	linal::Solver < double > A;

	array_t lh;

	void build_matrix();

public:

	SphereBaroclin(const Conf & conf);
	~SphereBaroclin();

	// u = &in[0]
	// v = &in[nlat * nlon]

	void S_step(double * out, const double * in, double t);
	void S_step(double * u1_o, double * u2_o, const double * u1, const double * u2, double t);

	void L_step(double * out, const double * in, const double * z);
	void LT_step(double * out, const double * in, const double * z);
	void L_1_step(double * out, const double * in, const double * z);

	void p2u(double * u, const double * p);
	void u2p(double * p, const double * u);
};

#endif /* BARVORTEX_H */

