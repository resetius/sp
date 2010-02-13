#include <math.h>

#include "chafe.h"
#include "utils.h"

SphereChafe::SphereChafe(const SphereChafeConf & conf)
	: conf(conf), lapl(conf.nlat, conf.nlon)
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

