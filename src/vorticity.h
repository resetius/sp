#ifndef VORTICITY_H
#define VORTICITY_H

#include <vector>
#include "operator.h"

class SphereVorticity: public SphereOperator
{
public:
	SphereVorticity(const SphereOperator & op);
	~SphereVorticity();

	void calc(double * div, const double * u, const double *v);
	void test();
};

#endif /* VORTICITY_H */
