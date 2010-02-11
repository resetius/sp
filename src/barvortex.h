#ifndef BARVORTEX_H
#define BARVORTEX_H

#include "lapl.h"
#include "jac.h"

struct SphereBarvortexConf;

typedef double (*coriolis_t) (double phi, double lambda);
typedef double (*rp_t ) (double phi, double lambda, double t, SphereBarvortexConf * conf );

struct SphereBarvortexConf
{
	long nlat;
	long nlon;
	double tau;
	double sigma;
	double mu;
	double k1;
	double k2;
	double theta;
};

class SphereBarvortex
{
	SphereBarvortexConf conf;
	SphereLaplace lapl;
	SphereJacobian jac;

public:

	SphereBarvortex(const SphereBarvortexConf & conf);
	~SphereBarvortex();

	void calc(double * out, const double * in, double t);
};

#endif /* BARVORTEX_H */

