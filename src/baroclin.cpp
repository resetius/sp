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
#include <string.h>
#include <assert.h>
#include <stdlib.h>

#include "baroclin.h"
#include "linal.h"

using namespace std;

SphereBaroclin::SphereBaroclin (const SphereBaroclin::Conf & conf) : SphereNorm (conf.nlat, conf.nlon),
		conf (conf), op (conf.nlat, conf.nlon, 0), lapl (op), jac (op),
		A (4 * 2 * conf.nlat * (std::min (conf.nlat, conf.nlon / 2 + 1) ) ),
		A1 (4 * 2 * conf.nlat * (std::min (conf.nlat, conf.nlon / 2 + 1) ) ),
		lh (conf.nlat * conf.nlon)
{
	long nlat = conf.nlat;
	long nlon = conf.nlon;
	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;


	for (int i = 0; i < nlat; ++i)
	{
		double phi    = -0.5 * M_PI + i * dlat;
		for (int j = 0; j < nlon; ++j)
		{
			double lambda = j * dlon;
			if (conf.cor)
			{
				lh[i * nlon + j] = conf.cor (phi, lambda, 0, &conf);
			}
			else if (conf.cor2)
			{
				lh[i * nlon + j] = conf.cor2[i * nlon + j];
			}
		}
	}

	build_matrix();
	build_inverted_matrix();
}

SphereBaroclin::~SphereBaroclin()
{
}

void SphereBaroclin::build_matrix()
{
	double tau   = conf.tau;
	double theta = conf.theta;
	double sigma = conf.sigma;
	double sigma1 = conf.sigma;
	double alpha = conf.alpha;
	double mu    = conf.mu;
	double mu1   = conf.mu1;

	long nlat = conf.nlat;
	long mmax = std::min (conf.nlat, conf.nlon / 2 + 1);
	long n1   = 2 * conf.nlat * (std::min (conf.nlat, conf.nlon / 2 + 1) );

	long i0, j0;

	i0 = 0;
	j0 = 0;

	for (int l = 0; l < 2; ++l)
	{
		for (long n = 0; n < nlat; ++n)
	//	for (long m = 0; m < mmax; ++m)
		{
			for (long m = 0; m < mmax; ++m)
	//		for (long n = 0; n < nlat; ++n)
			{
				int fnn = -n * (n + 1);

				// 1
				// w1
				A.add (i0, i0, 1. / tau + theta * sigma * 0.5 - theta * mu * fnn);
				// w2
				A.add (i0, i0 + n1, -theta * sigma * 0.5);

				// 2
				// w1
				A.add (i0 + n1, i0, theta * sigma * 0.5);
				// w2
				A.add (i0 + n1, i0 + n1,
				       1. / tau + theta * sigma * 0.5 - theta * mu * fnn +
				       + alpha * alpha * theta * mu1);
				// u2
				A.add (i0 + n1, i0 + 3 * n1, -alpha * alpha / tau - alpha * alpha * theta * sigma1);

				// 3
				// Delta u1 - w1 = 0
				// w1
				A.add (i0 + 2 * n1, i0, -1);
				// u1
				A.add (i0 + 2 * n1, i0 + 2 * n1, fnn ? fnn : -1);

				// 4
				// Delta u2 - w2 = 0
				// w2
				A.add (i0 + 3 * n1, i0 + n1, -1);
				// u2
				A.add (i0 + 3 * n1, i0 + 3 * n1, fnn ? fnn : -1);

				i0 += 1;
			}
		}
	}

	//A.print(stdout);
}

void SphereBaroclin::build_inverted_matrix()
{
	double tau   = conf.tau;
	double theta = conf.theta;
	double sigma = conf.sigma;
	double sigma1 = conf.sigma;
	double alpha = conf.alpha;
	double mu    = conf.mu;
	double mu1   = conf.mu1;

	long nlat = conf.nlat;
	long mmax = std::min (conf.nlat, conf.nlon / 2 + 1);
	long n1   = 2 * conf.nlat * (std::min (conf.nlat, conf.nlon / 2 + 1) );

	long i0, j0;

	i0 = 0;
	j0 = 0;

	for (int l = 0; l < 2; ++l)
	{
		for (long n = 0; n < nlat; ++n)
		{
			for (long m = 0; m < mmax; ++m)
			{
				int fnn = -n * (n + 1);

				// 1
				// w1
				A1.add (i0, i0, 1. / tau - (1-theta) * sigma * 0.5 + (1-theta) * mu * fnn);
				// w2
				A1.add (i0, i0 + n1, (1-theta) * sigma * 0.5);

				// 2
				// w1
				A1.add (i0 + n1, i0, -(1-theta) * sigma * 0.5);
				// w2
				A1.add (i0 + n1, i0 + n1,
				       1. / tau - (1-theta) * sigma * 0.5 + (1-theta) * mu * fnn
				       - alpha * alpha * (1-theta) * mu1);
				// u2
				A1.add (i0 + n1, i0 + 3 * n1, -alpha * alpha / tau + alpha * alpha * (1-theta) * sigma1);

				// 3
				// Delta u1 - w1 = 0
				// w1
				A1.add (i0 + 2 * n1, i0, -1);
				// u1
				A1.add (i0 + 2 * n1, i0 + 2 * n1, fnn ? fnn : -1);

				// 4
				// Delta u2 - w2 = 0
				// w2
				A1.add (i0 + 3 * n1, i0 + n1, -1);
				// u2
				A1.add (i0 + 3 * n1, i0 + 3 * n1, fnn ? fnn : -1);

				i0 += 1;
			}
		}
	}
}

void SphereBaroclin::S_step (double * out, const double * in, double t)
{
	long n       = conf.nlat * conf.nlon;
	const double * u1 = in;
	const double * u2 = &in[n];
	double * u11 = out;
	double * u21 = &out[n];

	S_step (u11, u21, u1, u2, t);
}

void SphereBaroclin::S_step (double * u11, double * u21, const double * u1, const double * u2, double t)
{
	using namespace linal;

	long nlat    = conf.nlat;
	long nlon    = conf.nlon;
	long n1      = 2 * conf.nlat * std::min (conf.nlat, conf.nlon / 2 + 1);
	long n       = conf.nlat * conf.nlon;
	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;
	double tau   = conf.tau;
	double theta = conf.theta;
	double mu    = conf.mu;
	double mu1   = conf.mu1;
	double sigma = conf.sigma;
	double sigma1 = conf.sigma1;
	double k1    = conf.k1;
	double k2    = conf.k2;
	double alpha = conf.alpha;

	array_t w1 (n);
	array_t w2 (n);
	array_t dw1 (n);
	array_t dw2 (n);

	array_t FC (n);
	array_t GC (n);

// next
	array_t u1_n (n);
	array_t u2_n (n);
	array_t w1_n (n);
	array_t w2_n (n);

	array_t u1_n1 (n);
	array_t u2_n1 (n);
	array_t w1_n1 (n);
	array_t w2_n1 (n);


// tmp
	array_t tmp1 (n);
	array_t tmp2 (n);

	// jac
	array_t jac1 (n);
	array_t jac2 (n);
	array_t jac3 (n);

	//
	array_t F (n);
	array_t G (n);
	array_t rp (4 * n1);
	array_t x (4 * n1);

	//

	lapl.calc (&w1[0], &u1[0]);
	lapl.calc (&w2[0], &u2[0]);

	lapl.calc (&dw1[0], &w1[0]);
	lapl.calc (&dw2[0], &w2[0]);

	vec_sum1 (&FC[0], &w1[0], &w2[0],
	          -0.5 * (1.0 - theta) * sigma,
	          0.5 * (1.0 - theta) * sigma, n);
	vec_sum1 (&FC[0], &FC[0], &dw1[0], 1.0, mu * (1.0 - theta), n);
	vec_sum1 (&FC[0], &FC[0], &w1[0], 1.0, 1.0 / tau, n);

	vec_sum1 (&GC[0], &w1[0], &w2[0],
	          -0.5 * (1.0 - theta) * sigma,
	          -0.5 * (1.0 - theta) * sigma, n);
	vec_sum1 (&GC[0], &GC[0], &dw2[0], 1.0, mu * (1.0 - theta), n);
	vec_sum1 (&GC[0], &GC[0], &w2[0], 1.0, 1.0 / tau, n);
	vec_sum1 (&GC[0], &GC[0], &u2[0], 1.0, -alpha * alpha / tau, n);
	vec_sum1 (&GC[0], &GC[0], &w2[0], 1.0, -alpha * alpha * mu1 * (1 - theta), n);
	vec_sum1 (&GC[0], &GC[0], &u2[0], 1.0, alpha * alpha * sigma1 * (1 - theta), n);

	memcpy (&u1_n[0], &u1[0], n * sizeof (double) );
	memcpy (&u2_n[0], &u2[0], n * sizeof (double) );
	memcpy (&w1_n[0], &w1[0], n * sizeof (double) );
	memcpy (&w2_n[0], &w2[0], n * sizeof (double) );

	for (int i = 0; i < nlat; ++i)
	{
		double phi    = -0.5 * M_PI + i * dlat;
		for (int j = 0; j < nlon; ++j)
		{
			double lambda = j * dlon;
			if (conf.rp1)
			{
				FC[i * nlon + j] += conf.rp1 (phi, lambda, t, &conf);
			}
			else if (conf.rp12)
			{
				FC[i * nlon + j] += conf.rp12[i * nlon + j];
			}

			if (conf.rp2)
			{
				if (alpha == 0)
				{
					GC[i * nlon + j] += conf.rp2 (phi, lambda, t, &conf);
				}
				else
				{
					GC[i * nlon + j] += alpha * alpha * conf.rp2 (phi, lambda, t, &conf);
				}
			}
			else if (conf.rp22)
			{
				GC[i * nlon + j] += alpha * alpha * conf.rp22[i * nlon + j];
			}
		}
	}

	for (int it = 0; it < 20; ++it)
	{
		// - J(0.5(u1+u1), 0.5(w1+w1)+l+h) - J(0.5(u2+u2),w2+w2)
		// J(0.5(u1+u1), 0.5(w1+w1)+l+h)
		vec_sum1 (&tmp1[0], &u1[0], &u1_n[0], 1.0 - theta, theta, n);
		vec_sum1 (&tmp2[0], &w1[0], &w1_n[0], 1.0 - theta, theta, n);
		vec_sum (&tmp2[0], &tmp2[0], &lh[0], n);
		jac.calc (&jac1[0], &tmp1[0], &tmp2[0]);
		// J(0.5(u2+u2),w2+w2)
		vec_sum1 (&tmp1[0], &u2[0], &u2_n[0], 1.0 - theta, theta, n);
		vec_sum1 (&tmp2[0], &w2[0], &w2_n[0], 1.0 - theta, theta, n);
		jac.calc (&jac2[0], &tmp1[0], &tmp2[0]);

		vec_sum1 (&F[0], &jac1[0], &jac2[0], -1.0, -1.0, n);

		// -J(0.5(u1+u1), 0.5(w2+w2)) - J(0.5(u2+u2), 0.5(w1+w1)+l+h) +
		// + alpha^2 J(0.5(u1+u1), 0.5(u2+u2))
		vec_sum1 (&tmp1[0], &u1[0], &u1_n[0], 1.0 - theta, theta, n);
		vec_sum1 (&tmp2[0], &w2[0], &w2_n[0], 1.0 - theta, theta, n);
		jac.calc (&jac1[0], &tmp1[0], &tmp2[0]);
		vec_sum1 (&tmp1[0], &u2[0], &u2_n[0], 1.0 - theta, theta, n);
		vec_sum1 (&tmp2[0], &w1[0], &w1_n[0], 1.0 - theta, theta, n);
		vec_sum (&tmp2[0], &tmp2[0], &lh[0], n);
		jac.calc (&jac2[0], &tmp1[0], &tmp2[0]);
		vec_sum1 (&tmp1[0], &u1[0], &u1_n[0], 1.0 - theta, theta, n);
		vec_sum1 (&tmp2[0], &u2[0], &u2_n[0], 1.0 - theta, theta, n);
		jac.calc (&jac3[0], &tmp1[0], &tmp2[0]);
		vec_sum1 (&G[0], &jac1[0], &jac2[0], -1.0, -1.0, n);
		vec_sum1 (&G[0], &G[0], &jac3[0], 1.0, alpha * alpha, n);

		vec_sum (&F[0], &F[0], &FC[0], n);
		vec_sum (&G[0], &G[0], &GC[0], n);

		// 1. build right part [F, G, 0, 0]
		// 2. build right part koefs
		memset (&rp[0], 0, 4 * n1);
		op.func2koef (&rp[0],    &F[0]);
		op.func2koef (&rp[n1],   &G[0]);
		//op.func2koef (&rp[2*n1], &w1_n[0]);
		//op.func2koef (&rp[3*n1], &w2_n[0]);
		// selftest
		if (0)
		{
			op.koef2func (&w1_n1[0], &rp[0]);
			op.koef2func (&w2_n1[0], &rp[n1]);
			op.koef2func (&u1_n1[0], &rp[2 * n1]);
			op.koef2func (&u2_n1[0], &rp[3 * n1]);
			fprintf (stderr, "1: %le\n", dist (&w1_n1[0], &F[0]) );
			fprintf (stderr, "2: %le\n", dist (&w2_n1[0], &G[0]) );
			fprintf (stderr, "3: %le\n", dist (&u1_n1[0], &w1_n[0]) );
			fprintf (stderr, "4: %le\n", dist (&u2_n1[0], &w2_n[0]) );
			exit (1);
		}

		// 3. solve equation and find koefs

		//mat_print("rp.txt", &rp[0], rp.size(), 1, "%8.3le ");
		//exit(1);
		A.solve (&x[0], &rp[0]);
		// 4. build functions from koefs
		op.koef2func (&w1_n[0],  &x[0]);
		op.koef2func (&w2_n[0],  &x[n1]);
		op.koef2func (&u1_n1[0], &x[2 * n1]);
		op.koef2func (&u2_n1[0], &x[3 * n1]);

		double nr1 = dist (&u1_n1[0], &u1_n[0]);
		double nr2 = dist (&u2_n1[0], &u2_n[0]);
		double nr  = std::max (nr1, nr2);
		u1_n1.swap (u1_n);
		u2_n1.swap (u2_n);

		if (nr < 1e-8)
		{
			break;
		}
	}

	memcpy (u11, &u1_n[0], n * sizeof (double) );
	memcpy (u21, &u2_n[0], n * sizeof (double) );
}

void SphereBaroclin::L_step (double *u11, double * u21,
                             const double * u1, const double * u2, 
                             const double * z1, const double * z2)
{
	using namespace linal;

	long nlat    = conf.nlat;
	long nlon    = conf.nlon;
	long n1      = 2 * conf.nlat * std::min (conf.nlat, conf.nlon / 2 + 1);
	long n       = conf.nlat * conf.nlon;
	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;
	double tau   = conf.tau;
	double theta = conf.theta;
	double mu    = conf.mu;
	double mu1   = conf.mu1;
	double sigma = conf.sigma;
	double sigma1 = conf.sigma1;
	double k1    = conf.k1;
	double k2    = conf.k2;
	double alpha = conf.alpha;

	// правая часть 1:
	// - J(z1, 0.5(w1+w1)) - J(u1, L(z1)+l+h) -
	// - J(z2, 0.5(w2+w2)) - J(0.5(u2+u2), L(z2)) +
	// + w1/tau - 0.5 (1-theta)sigma (w1-w2)+mu(1-theta)(L w1)
	// правая часть 2:
	// - J(z1, 0.5(w2+w2)) - J(0.5(u1+u1), L(z2)) -
	// - J(0.5(u2+u2), L(z1)+l+h) - J(z2, 0.5(w1+w1)) -
	// - 0.5 (1-theta)sigma (w1 + w2) + (1-theta) mu L w2
	// + w2/tau - alpha^2 u2/tau +
	// + alpha^2 J(z1, 0.5(u2+u2)) + alpha^2 J(0.5(u1+u1), z2) -
	// - alpha^2 (1-theta) mu1 w2 +
	// + alpha^2 sigma1 (1-theta) u2 

	array_t w1(n);
	array_t w2(n);
	array_t dw1(n);
	array_t dw2(n);
	array_t dz1(n);
	array_t dz2(n);
	array_t FC(n);
	array_t GC(n);

	// next
	array_t u1_n(n);
	array_t u2_n(n);
	array_t w1_n(n);
	array_t w2_n(n);

	array_t u1_n1(n);
	array_t u2_n1(n);

	// tmp
	array_t tmp1(n);
	array_t tmp2(n);

	// jac!
	array_t jac1(n);
	array_t jac2(n);
	array_t jac3(n);

	//
	array_t F(n);
	array_t G(n);

	array_t rp(4*n1);
	array_t x(4*n1);

	lapl.calc(&w1[0], &u1[0]);
	lapl.calc(&w2[0], &u2[0]);

	lapl.calc(&dw1[0], &w1[0]);
	lapl.calc(&dw2[0], &w2[0]);

	lapl.calc(&dz1[0], &z1[0]);
	lapl.calc(&dz2[0], &z2[0]);

	// w1/tau - 0.5 (1-theta)sigma(w1-w2) + mu(1-theta)(L w1)
	vec_sum1(&FC[0], &w1[0], &w2[0], 
		-0.5 * (1.0 - theta) * sigma, 0.5 * (1.0 - theta) * sigma, n);
	vec_sum1(&FC[0], &FC[0], &dw1[0], 1.0, mu * (1.0 - theta), n);
	vec_sum1(&FC[0], &FC[0], &w1[0], 1.0, 1.0 / tau, n);

	// w2/tau - 0.5 (1-theta)sigma (w1 + w2) + (1-theta) mu L w2 -
	// - alpha^2 u2/tau - alpha^2 (1-theta) mu1 w2 + alpha^2 sigma1 (1-theta) u2
	vec_sum1(&GC[0], &w1[0], &w2[0],
			-0.5 * (1.0 - theta) * sigma, -0.5 * (1.0 - theta) * sigma, n);
	vec_sum1(&GC[0], &GC[0], &dw2[0], 1.0, mu * (1.0 - theta), n);
	vec_sum1(&GC[0], &GC[0], &w2[0], 1.0, 1.0 / tau, n);
	vec_sum1(&GC[0], &GC[0], &u2[0], 1.0, -alpha * alpha / tau, n);
	vec_sum1(&GC[0], &GC[0], &w2[0], 1.0, -alpha * alpha * mu1 * (1-theta), n);
	vec_sum1(&GC[0], &GC[0], &u2[0], 1.0, alpha * alpha * sigma1 * (1-theta), n);

	memcpy(&u1_n[0], &u1[0], n * sizeof(double));
	memcpy(&u2_n[0], &u2[0], n * sizeof(double));
	memcpy(&w1_n[0], &w1[0], n * sizeof(double));
	memcpy(&w2_n[0], &w2[0], n * sizeof(double));

	for (int it = 0; it < 100; ++it) {
		// - J(0.5(u1+u1), L(z1)+l+h) - J(z1, 0.5(w1+w1)) -
		// - J(z2,0.5(w2+w2)) - J(0.5(u2+u2),L(z2))

		// J(0.5(u1+u1), L(z1)+l+h)
		vec_sum1(&tmp1[0], &u1[0], &u1_n[0], 1.0 - theta, theta, n);
		vec_sum(&tmp2[0], &dz1[0], &lh[0], n);
		jac.calc(&jac1[0], &tmp1[0], &tmp2[0]);

		// J(z1, 0.5(w1+w1))
		vec_sum1(&tmp1[0], &w1[0], &w1_n[0], 1.0 - theta, theta, n);
		jac.calc(&jac2[0], &z1[0], &tmp1[0]);
		vec_sum1(&F[0], &jac1[0], &jac2[0], -1.0, -1.0, n);

		// J(z2,0.5(w2+w2))
		vec_sum1(&tmp1[0], &w2[0], &w2_n[0], 1.0 - theta, theta, n);
		jac.calc(&jac1[0], &z2[0], &tmp1[0]);
		vec_sum1(&F[0], &F[0], &jac1[0], 1.0, -1.0, n);

		// J(0.5(u2+u2),L(z2))
		vec_sum1(&tmp1[0], &u2[0], &u2_n[0], 1.0 - theta, theta, n);
		jac.calc(&jac1[0], &tmp1[0], &dz2[0]);
		vec_sum1(&F[0], &F[0], &jac1[0], 1.0, -1.0, n);
	
		// - J(z1, 0.5(w2+w2)) - J(0.5(u1+u1), L(z2)) -
		// - J(0.5(u2+u2), L(z1)+l+h) - J(z2, 0.5(w1+w1)) +
		// + alpha^2 J(z1, 0.5(u2+u2)) + alpha^2 J(0.5(u1+u1), z2))

		// J(z1, 0.5(w2+w2))
		vec_sum1(&tmp1[0], &w2[0], &w2_n[0], 1.0 - theta, theta, n);
		jac.calc(&jac1[0], &dz1[0], &tmp1[0]);

		// J(0.5(u1+u1), L(z2))
		vec_sum1(&tmp1[0], &u1[0], &u1_n[0], 1.0 - theta, theta, n);
		jac.calc(&jac2[0], &dz2[0], &tmp1[0]);
		vec_sum1(&G[0], &jac1[0], &jac2[0], -1.0, -1.0, n);

		// J(0.5(u2+u2), L(z1)+l+h)
		vec_sum1(&tmp1[0], &u2[0], &u2_n[0], 1.0 - theta, theta, n);
		vec_sum(&tmp2[0], &dz1[0], &lh[0], n);
		jac.calc(&jac1[0], &tmp1[0], &tmp2[0]);
		vec_sum1(&G[0], &G[0], &jac1[0], 1.0, -1.0, n);

		// J(z2, 0.5(w1+w1))
		vec_sum1(&tmp1[0], &w1[0], &w1_n[0], 1.0 - theta, theta, n);
		jac.calc(&jac1[0], &z2[0], &tmp1[0]);
		vec_sum1(&G[0], &G[0], &jac1[0], 1.0, -1.0, n);

		// alpha^2 J(z1, 0.5(u2+u2))
		vec_sum1(&tmp1[0], &u2[0], &u2_n[0], 1.0 - theta, theta, n);
		jac.calc(&jac1[0], &z1[0], &tmp1[0]);
		vec_sum1(&G[0], &G[0], &jac1[0], 1.0, alpha * alpha, n);

		// alpha^2 J(0.5(u1+u1), z2))
		vec_sum1(&tmp1[0], &u1[0], &u1_n[0], 1.0 - theta, theta, n);
		jac.calc(&jac1[0], &tmp1[0], &z2[0]);
		vec_sum1(&G[0], &G[0], &jac1[0], 1.0, alpha * alpha, n);

		vec_sum (&F[0], &F[0], &FC[0], n);
		vec_sum (&G[0], &G[0], &GC[0], n);

		memset (&rp[0], 0, 4 * n1);
		op.func2koef (&rp[0],    &F[0]);
		op.func2koef (&rp[n1],   &G[0]);

		A.solve (&x[0], &rp[0]);

		op.koef2func (&w1_n[0],  &x[0]);
		op.koef2func (&w2_n[0],  &x[n1]);
		op.koef2func (&u1_n1[0], &x[2 * n1]);
		op.koef2func (&u2_n1[0], &x[3 * n1]);

		double nr1 = dist (&u1_n1[0], &u1_n[0]);
		double nr2 = dist (&u2_n1[0], &u2_n[0]);
		double nr  = std::max (nr1, nr2);
		u1_n1.swap (u1_n);
		u2_n1.swap (u2_n);

		if (nr < 1e-8) {
			break;
		}
	}

	memcpy(u11, &u1_n[0], n * sizeof(double));
	memcpy(u21, &u2_n[0], n * sizeof(double));
}

void SphereBaroclin::L_step (double *out, const double *in, const double * z)
{
	long n       = conf.nlat * conf.nlon;
	const double * u1 = in;
	const double * u2 = &in[n];
	const double * z1 = z;
	const double * z2 = &z[n];
	double * u11 = out;
	double * u21 = &out[n];

	L_step (u11, u21, u1, u2, z1, z2);
}

void SphereBaroclin::L_1_step (double *u11, double * u21,
                               const double * u1, const double * u2, 
                               const double * z1, const double * z2)
{
	using namespace linal;

	long nlat    = conf.nlat;
	long nlon    = conf.nlon;
	long n1      = 2 * conf.nlat * std::min (conf.nlat, conf.nlon / 2 + 1);
	long n       = conf.nlat * conf.nlon;
	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;
	double tau   = conf.tau;
	double theta = conf.theta;
	double mu    = conf.mu;
	double mu1   = conf.mu1;
	double sigma = conf.sigma;
	double sigma1 = conf.sigma1;
	double k1    = conf.k1;
	double k2    = conf.k2;
	double alpha = conf.alpha;

	// правая часть 1:
	//  J(z1, 0.5(w1+w1)) + J(u1, L(z1)+l+h) +
	// + J(z2, 0.5(w2+w2)) + J(0.5(u2+u2), L(z2)) +
	// + w1/tau + 0.5 theta sigma (w1-w2) - mu theta (L w1)
	// правая часть 2:
	//  J(z1, 0.5(w2+w2)) + J(0.5(u1+u1), L(z2)) +
	// + J(0.5(u2+u2), L(z1)+l+h) - J(z2, 0.5(w1+w1)) +
	// + 0.5 theta sigma (w1 + w2) - theta mu L w2 +
	// + w2/tau - alpha^2 u2/tau 
	// - alpha^2 J(z1, 0.5(u2+u2)) - alpha^2 J(0.5(u1+u1), z2) +
	// + alpha^2 theta mu1 w2 -
	// - alpha^2 sigma1 theta u2 

	array_t w1(n);
	array_t w2(n);
	array_t dw1(n);
	array_t dw2(n);
	array_t dz1(n);
	array_t dz2(n);
	array_t FC(n);
	array_t GC(n);

	// next
	array_t u1_n(n);
	array_t u2_n(n);
	array_t w1_n(n);
	array_t w2_n(n);

	array_t u1_n1(n);
	array_t u2_n1(n);

	// tmp
	array_t tmp1(n);
	array_t tmp2(n);

	// jac
	array_t jac1(n);
	array_t jac2(n);
	array_t jac3(n);

	//
	array_t F(n);
	array_t G(n);

	array_t rp(4*n1);
	array_t x(4*n1);

	lapl.calc(&w1[0], &u1[0]);
	lapl.calc(&w2[0], &u2[0]);

	lapl.calc(&dw1[0], &w1[0]);
	lapl.calc(&dw2[0], &w2[0]);

	lapl.calc(&dz1[0], &z1[0]);
	lapl.calc(&dz2[0], &z2[0]);

	// w1/tau + 0.5 theta sigma(w1-w2) - mu theta(L w1)
	vec_sum1(&FC[0], &w1[0], &w2[0], 
		0.5 * theta * sigma, -0.5 * theta * sigma, n);
	vec_sum1(&FC[0], &FC[0], &dw1[0], 1.0, -mu * theta, n);
	vec_sum1(&FC[0], &FC[0], &w1[0], 1.0, 1.0 / tau, n);

	// w2/tau + 0.5 theta sigma (w1 + w2) - theta mu L w2 -
	// - alpha^2 u2/tau + alpha^2 theta mu1 w2 - alpha^2 sigma1 theta u2
	vec_sum1(&GC[0], &w1[0], &w2[0],
			0.5 * theta * sigma, 0.5 * theta * sigma, n);
	vec_sum1(&GC[0], &GC[0], &dw2[0], 1.0, -mu * theta, n);
	vec_sum1(&GC[0], &GC[0], &w2[0], 1.0, 1.0 / tau, n);
	vec_sum1(&GC[0], &GC[0], &u2[0], 1.0, -alpha * alpha / tau, n);
	vec_sum1(&GC[0], &GC[0], &w2[0], 1.0, alpha * alpha * mu1 * theta, n);
	vec_sum1(&GC[0], &GC[0], &u2[0], 1.0, -alpha * alpha * sigma1 * theta, n);

	memcpy(&u1_n[0], &u1[0], n * sizeof(double));
	memcpy(&u2_n[0], &u2[0], n * sizeof(double));
	memcpy(&w1_n[0], &w1[0], n * sizeof(double));
	memcpy(&w2_n[0], &w2[0], n * sizeof(double));

	for (int it = 0; it < 1; ++it) {
		//  J(0.5(u1+u1), L(z1)+l+h) + J(z1, 0.5(w1+w1)) +
		// + J(z2,0.5(w2+w2)) + J(0.5(u2+u2),L(z2))

		// J(0.5(u1+u1), L(z1)+l+h)
		vec_sum1(&tmp1[0], &u1[0], &u1_n[0], theta, 1.0 - theta, n);
		vec_sum(&tmp2[0], &dz1[0], &lh[0], n);
		jac.calc(&jac1[0], &tmp1[0], &tmp2[0]);

		// J(z1, 0.5(w1+w1))
		vec_sum1(&tmp1[0], &w1[0], &w1_n[0], theta, 1.0 - theta, n);
		jac.calc(&jac2[0], &z1[0], &tmp1[0]);
		vec_sum(&F[0], &jac1[0], &jac2[0], n);

		// J(z2,0.5(w2+w2))
		vec_sum1(&tmp1[0], &w2[0], &w2_n[0], theta, 1.0 - theta, n);
		jac.calc(&jac1[0], &z2[0], &tmp1[0]);
		vec_sum(&F[0], &F[0], &jac1[0], n);

		// J(0.5(u2+u2),L(z2))
		vec_sum1(&tmp1[0], &u2[0], &u2_n[0], theta, 1.0 - theta, n);
		jac.calc(&jac1[0], &tmp1[0], &dz2[0]);
		vec_sum(&F[0], &F[0], &jac1[0], n);
	
		//  J(z1, 0.5(w2+w2)) + J(0.5(u1+u1), L(z2)) +
		// + J(0.5(u2+u2), L(z1)+l+h) + J(z2, 0.5(w1+w1)) -
		// - alpha^2 J(z1, 0.5(u2+u2)) - alpha^2 J(0.5(u1+u1), z2))

		// J(z1, 0.5(w2+w2))
		vec_sum1(&tmp1[0], &w2[0], &w2_n[0], theta, 1.0 - theta, n);
		jac.calc(&jac1[0], &dz1[0], &tmp1[0]);

		// J(0.5(u1+u1), L(z2))
		vec_sum1(&tmp1[0], &u1[0], &u1_n[0], theta, 1.0 - theta, n);
		jac.calc(&jac2[0], &dz2[0], &tmp1[0]);
		vec_sum(&G[0], &jac1[0], &jac2[0], n);

		// J(0.5(u2+u2), L(z1)+l+h)
		vec_sum1(&tmp1[0], &u2[0], &u2_n[0], theta, 1.0 - theta, n);
		vec_sum(&tmp2[0], &dz1[0], &lh[0], n);
		jac.calc(&jac1[0], &tmp1[0], &tmp2[0]);
		vec_sum(&G[0], &G[0], &jac1[0], n);

		// J(z2, 0.5(w1+w1))
		vec_sum1(&tmp1[0], &w1[0], &w1_n[0], theta, 1.0 - theta, n);
		jac.calc(&jac1[0], &z2[0], &tmp1[0]);
		vec_sum(&G[0], &G[0], &jac1[0], n);

		// alpha^2 J(z1, 0.5(u2+u2))
		vec_sum1(&tmp1[0], &u2[0], &u2_n[0], theta, 1.0 - theta, n);
		jac.calc(&jac1[0], &z1[0], &tmp1[0]);
		vec_sum1(&G[0], &G[0], &jac1[0], 1.0, -alpha * alpha, n);

		// alpha^2 J(0.5(u1+u1), z2))
		vec_sum1(&tmp1[0], &u1[0], &u1_n[0], theta, 1.0 - theta, n);
		jac.calc(&jac1[0], &tmp1[0], &z2[0]);
		vec_sum1(&G[0], &G[0], &jac1[0], 1.0, -alpha * alpha, n);

		vec_sum (&F[0], &F[0], &FC[0], n);
		vec_sum (&G[0], &G[0], &GC[0], n);

		memset (&rp[0], 0, 4 * n1);
		op.func2koef (&rp[0],    &F[0]);
		op.func2koef (&rp[n1],   &G[0]);

		A1.solve (&x[0], &rp[0]);

		op.koef2func (&w1_n[0],  &x[0]);
		op.koef2func (&w2_n[0],  &x[n1]);
		op.koef2func (&u1_n1[0], &x[2 * n1]);
		op.koef2func (&u2_n1[0], &x[3 * n1]);

		double nr1 = dist (&u1_n1[0], &u1_n[0]);
		double nr2 = dist (&u2_n1[0], &u2_n[0]);
		double nr  = std::max (nr1, nr2);
		u1_n1.swap (u1_n);
		u2_n1.swap (u2_n);

		if (nr < 1e-8) {
			break;
		}
	}

	memcpy(u11, &u1_n[0], n * sizeof(double));
	memcpy(u21, &u2_n[0], n * sizeof(double));
}

void SphereBaroclin::L_1_step (double *out, const double *in, const double * z)
{
	long n       = conf.nlat * conf.nlon;
	const double * u1 = in;
	const double * u2 = &in[n];
	const double * z1 = z;
	const double * z2 = &z[n];
	double * u11 = out;
	double * u21 = &out[n];

	L_1_step (u11, u21, u1, u2, z1, z2);
}

void SphereBaroclin::p2u (double * u, const double * p)
{
	memcpy (u, p, conf.nlat * conf.nlon * sizeof (double) );
}

void SphereBaroclin::u2p (double * p, const double * u)
{
	memcpy (p, u, conf.nlat * conf.nlon * sizeof (double) );
}

