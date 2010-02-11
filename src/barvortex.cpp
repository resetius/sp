#include <math.h>
#include <string.h>

#include "barvortex.h"
#include "utils.h"

using namespace std;

SphereBarvortex::SphereBarvortex (const SphereBarvortexConf & conf) :
		conf (conf), lapl (conf.nlat, conf.nlon), jac (conf.nlat, conf.nlon),
		lh (conf.nlat * conf.nlon)
{
	if (conf.coriolis)
	{
		long nlat = conf.nlat;
		long nlon = conf.nlon;
		double dlat = M_PI / (nlat-1);
		double dlon = 2. * M_PI /nlon;

		for (int i = 0; i < nlat; ++i)
		{
			double phi    = -0.5 * M_PI + i * dlat;
			for (int j = 0; j < nlon; ++j)
			{
				double lambda = j * dlon;
				lh[i * nlon + j] = conf.coriolis(phi, lambda);
			}
		}
	}
}

SphereBarvortex::~SphereBarvortex()
{
}

void SphereBarvortex::calc (double * out, const double * u, double t)
{
	long nlat    = conf.nlat;
	long nlon    = conf.nlon;
	long n       = conf.nlat * conf.nlon;
	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;
	double tau   = conf.tau;
	double theta = conf.theta;
	double mu    = conf.mu;
	double sigma = conf.sigma;
	double k1    = conf.k1;
	double k2    = conf.k2;
	double nr, nr0;

	array_t w (n);  // w = L(u)
	array_t dw (n); // dw = L(w) = LL(u)

	// next
	array_t u_n (n);
	array_t w_n (n);
	array_t u_n1 (n);

	// jac
	array_t jac0(n);

	// tmp
	array_t tmp1(n);
	array_t tmp2(n);

	//
	array_t FC (n);
	array_t F (n);

	// генерируем правую часть
	// w/dt + mu (1-theta) L w - \sigma(1-theta) w -
	// - k1 J(0.5(u+u), 0.5(w+w)) - k2 J(0.5(u+u), l + h) + f(x, y)

	// w = L (u)
	lapl.solve (&w[0], &u[0]);
	// dw = L (w)
	lapl.solve (&dw[0], &w[0]);

	// w/dt + mu (1-theta) L w
	vec_sum1 (&FC[0], &w[0], &dw[0], 1.0 / tau,
	          mu * (1.0 - theta), n);
	// w/dt + mu (1-theta) L w - \sigma (1-theta) w
	vec_sum1 (&FC[0], &FC[0], &w[0], 1.0,
	          -sigma * (1.0 - theta), n);

	for (int i = 0; i < nlat; ++i)
	{
		double phi    = -0.5 * M_PI + i * dlat;
		for (int j = 0; j < nlon; ++j)
		{
			double lambda = j * dlon;
			if (conf.rp)
			{
				FC[i * nlon + j] += conf.rp (phi, lambda, t, &conf);
			}
		}
	}
	// теперь прибавим часть якобиана
	jac.calc(&jac0[0], &u[0], &lh[0]);
	vec_sum1(&FC[0], &FC[0], &jac0[0], 1.0, -(1.0 - theta), n);

	memcpy(&u_n[0], &u[0], n * sizeof(double));
	memcpy(&w_n[0], &w[0], n * sizeof(double));

	// в FC содержится правая часть, которая не меняется при итерациях!
	
	for (int it = 0; it < 1000; ++it) {
		// k1 J(0.5(u+u), 0.5(w+w)) + k2 J(0.5(u+u), l + h)   =
		// = J(0.5 (u+u), 0.5 k1 (w+w)) + J(0.5 (u+u), k2 (l + h)) =
		// = J(0.5 (u+u), 0.5 k1 (w+w) + k2 (l + h))
		vec_sum1(&tmp1[0], &w_n[0], &w[0], k1 * theta,
				k1 * (1.0 - theta), n);
		vec_sum1(&tmp1[0], &tmp1[0], &lh[0], 1.0, k2, n);
		// 0.5(u+u)
		vec_sum1(&tmp2[0], &u_n[0], &u[0], theta,
				1.0 - theta, n);
		// - k1 J(0.5(u+u), 0.5(w+w)) - k2 J(0.5(u+u), l + h)
		jac.calc(&jac0[0], &tmp2[0], &tmp1[0]);
		vec_sum1(&F[0], &FC[0], &jac0[0], 1.0, -1.0, n);
		lapl.solve(&w_n[0], &F[0], -0.5 * mu, 1.0 / tau + 0.5 * sigma);
		lapl.solve(&u_n1[0], &w_n[0]);

		nr = dist(&u_n1[0], &u_n[0]);
		u_n1.swap(u_n);
		if (nr / nr0 < 1e-14 || isnan(nr)) {
			break;
		}
	}

	memcpy(out, &u_n[0], n * sizeof(double));
}

