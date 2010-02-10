#ifndef CHAFE_H
#define CHAFE_H

#include "lapl.h"

typedef double (*rp_t)(double phi, double lambda, double t, double mu, double sigma);

struct SphereChafeConf
{
	long nlat;
	long nlon;
	double tau;
	double mu;
	double sigma;
	rp_t rp;
};

/**
 * \f$\frac{du}{dt} = \mu \delta u - \sigma u + f (u)\f$
 */
class SphereChafe {
	typedef std::vector < double > array_t;
	SphereChafeConf conf;
	SphereLaplace lapl;

public:
	SphereChafe(const SphereChafeConf & conf);

	void calc(double * out, const double * in, double t);
};

#endif /* CHAFE_H */

