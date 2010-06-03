#ifndef DIV_H
#define DIV_H

#include <vector>
#include "operator.h"

class SphereDiv: public SphereOperator
{
public:
	SphereDiv(const SphereOperator & op);
	~SphereDiv();

	void calc(double * div, const double * u, const double *v);
};

#endif /* DIV_H */
