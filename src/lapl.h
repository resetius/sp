#ifndef LAPL_H
#define LAPL_H

#include <vector>
#include "spherepack.h"

class SphereLaplace {
	typedef std::vector < double > array_t;

	long nlat;
	long nlon;

	long lshaec;
	long lshsec;
	long lwork;
	long ldwork;

	array_t wshaec;
	array_t wshsec;
	array_t work;
	array_t dwork;

public:
	SphereLaplace(long nlat, long nlon);
	~SphereLaplace();

	void solve(double * out, const double * in, double diag);
};

#endif /* LAPL_H */
