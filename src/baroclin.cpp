#include <math.h>
#include <string.h>
#include <assert.h>

#include "baroclin.h"
#include "linal.h"

using namespace std;

SphereBaroclin::SphereBaroclin (const SphereBaroclin::Conf & conf) : SphereNorm (conf.nlat, conf.nlon),
		conf (conf), op (conf.nlat, conf.nlon, 0), lapl (op), jac (op),
		A (4 * conf.nlat * std::min(conf.nlat, conf.nlon/2+1)),
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
}

SphereBaroclin::~SphereBaroclin()
{
}

void SphereBaroclin::build_matrix()
{
	double tau   = conf.tau;
	double theta = conf.theta;
	double sigma = conf.sigma;
	double sigma1= conf.sigma;
	double alpha = conf.alpha;
	double mu    = conf.mu;
	double mu1   = conf.mu1;

	long nlat = conf.nlat;
	long mmax = std::min(conf.nlat, conf.nlon/2+1);
	long n1   = conf.nlat * std::min(conf.nlat, conf.nlon/2+1);

	long i0, j0;

	i0 = 0; j0 = 0;
	for (long m = 0; m < mmax; ++m)
	{
		for (long n = 0; n < nlat; ++n)
		{
			double fnn = -n * (n + 1);

			// 1
			// w1
			A.add(i0, i0, 1. / tau + theta * sigma * 0.5 - theta * mu * fnn);
			// w2
			A.add(i0, i0 + n1, -theta * sigma * 0.5);

			// 2
			// w1
			A.add(i0 + n1, i0, theta * sigma * 0.5);
			// w2
			A.add(i0 + n1, i0 + n1, 
			      1. / tau + theta * sigma * 0.5 - theta * mu * fnn + 
			      + alpha * alpha * theta * mu1);
			// u2
			A.add(i0 + n1, i0 + 3 * n1, -alpha * alpha / tau - alpha * alpha * sigma1);

			// 3
			// u1
			A.add(i0 + 2 * n1, i0 + 2 * n1, fnn);

			// u2
			A.add(i0 + 2 * n1, i0 + 3 * n1, fnn);

			i0 += 1;
		}
	}

	A.print();
}

void SphereBaroclin::S_step (double * out, const double * in, double t)
{
	long n       = conf.nlat * conf.nlon;
	const double * u1 = in;
	const double * u2 = &in[n];
	double * u11 = out;
	double * u21 = &out[n];

	S_step(u11, u21, u1, u2, t);
}

void SphereBaroclin::S_step (double * u11, double * u21, const double * u1, const double * u2, double t)
{
	using namespace linal;

	long nlat    = conf.nlat;
	long nlon    = conf.nlon;
	long n1      = 2 * conf.nlat * conf.nlat;
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
	array_t rp(4 * n1);
	array_t x(4 * n1);

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
		memset(&rp[0], 0, 4 * n1);
		op.func2koef(&rp[0],    &F[0]);
		op.func2koef(&rp[n1],   &G[0]);
		op.func2koef(&rp[2*n1], &w1_n[0]);
		op.func2koef(&rp[3*n1], &w2_n[0]);

		// 3. solve equation and find koefs
		A.solve(&x[0], &rp[0]);
		// 4. build functions from koefs
		op.koef2func(&w1_n[0], &x[0]);
		op.koef2func(&w2_n[0], &x[n1]);
		op.koef2func(&u1_n[0], &x[2 * n1]);
		op.koef2func(&u2_n[0], &x[3 & n1]);

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

void SphereBaroclin::L_step (double *u1, const double *u, const double * z)
{
	assert(0);
}

void SphereBaroclin::LT_step (double *v1, const double *v, const double * z)
{
	assert(0);
}

void SphereBaroclin::L_1_step (double *u1, const double *u, const double * z)
{
	assert(0);
}

void SphereBaroclin::p2u (double * u, const double * p)
{
	memcpy (u, p, conf.nlat * conf.nlon * sizeof (double) );
}

void SphereBaroclin::u2p (double * p, const double * u)
{
	memcpy (p, u, conf.nlat * conf.nlon * sizeof (double) );
}
