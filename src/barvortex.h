#ifndef BARVORTEX_H
#define BARVORTEX_H

#include <vector>

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

	coriolis_t coriolis;
	rp_t rp;
};

class SphereBarvortex
{
	typedef std::vector < double > array_t;

	SphereBarvortexConf conf;
	SphereLaplace lapl;
	SphereJacobian jac;

	array_t lh;

public:

	SphereBarvortex(const SphereBarvortexConf & conf);
	~SphereBarvortex();

	void calc(double * out, const double * in, double t);
	double scalar(const double * u, const double * v);
	double dist(const double * u, const double * v);
	double norm(const double * u);
};

#endif /* BARVORTEX_H */
