#include <math.h>

#include "chafe.h"
#include "utils.h"

SphereChafe::SphereChafe(const SphereChafeConf & conf)
	: conf(conf), lapl(conf.nlat, conf.nlon)
{
}

void SphereChafe::calc(double * out, const double * in)
{
	double tau   = conf.tau;
	double mu    = conf.mu;
	double sigma = conf.sigma;
	double nlat  = conf.nlat;
	double nlon  = conf.nlon;
	int n = nlat * nlon;

	double dlat = M_PI / (nlat-1);
        double dlon = 2. * M_PI /nlon;

	array_t u(n);
	array_t delta_u(n);
	array_t crp(n);

	lapl.calc(&delta_u[0], in);
	// u/dt + mu \Delta u / 2
	vec_sum1(&delta_u[0], &u[0], &delta_u[0], (1.0 / tau), (mu * 0.5), n);

	// u/dt + mu \Delta u / 2 - \sigma u / 2
	vec_sum1(&delta_u[0], &delta_u[0], &u[0], 1.0, (-sigma * 0.5), n);

	// u/dt + mu \Delta u / 2 - \sigma u / 2 + f
	for (int i = 0; i < nlat; i++) {
		double phi    = -0.5 * M_PI + i * dlat;
		for (int j = 0; j < nlon; ++j) {
			double lambda = j * dlon;
			if (conf.rp) {
				delta_u[i * nlon + j] += conf.rp(phi, lambda);
			}
		}
	}

	// u'/dt - mu \Delta u' / 2 + \sigma u' / 2
	lapl.solve(out, &delta_u[0], -0.5 * mu, 1.0 / tau + sigma);
}

