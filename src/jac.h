#ifndef SPHERE_JAC_H
#define SPHERE_JAC_H

#include <vector>
#include "operator.h"

class SphereJacobian: public SphereOperator
{
public:
	SphereJacobian(const SphereOperator & op);
	~SphereJacobian();

	void calc(double * out, const double * u, const double * v);
	void calc_t(double * out, const double * u, const double * v);
};

#endif /* SPHERE_JAC_H */

