#ifndef GRAD_H
#define GRAD_H

#include <vector>
#include "operator.h"

class SphereGrad: public SphereOperator
{
public:
	SphereGrad(const SphereOperator & op);
	~SphereGrad();

	void calc (double * duphi, double * dulambda, const double * u);
	void solve (double * u, const double * duphi, const double * dulambda);
};

#endif /* GRAD_H */
