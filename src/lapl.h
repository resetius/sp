#ifndef LAPL_H
#define LAPL_H

#include <vector>
#include "spherepack.h"
#include "operator.h"

class SphereLaplace: public SphereOperator
{
public:
	SphereLaplace (const SphereOperator & op);
	~SphereLaplace();

	/**
	 * out[i][j]
	 * i - latitude  (from the sourth pole to the north pole)
	 * j - longitude
	 **/
	void solve (double * out, const double * in, double mult = 1.0, double diag = 0.0);
	void calc (double * out, const double * in);
	void make_psi (double * psi, const double * u, const double * v);
};

#endif /* LAPL_H */

