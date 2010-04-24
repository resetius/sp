#ifndef LAPL_H
#define LAPL_H

#include <vector>
#include "spherepack.h"

class SphereLaplace {
	typedef std::vector < double > array_t;

	long nlat;
	long nlon;
	long isym;

	long lshaec;
	long lshsec;
	long lwork;
	long ldwork;

	array_t wshaec;
	array_t wshsec;
	array_t work;
	array_t dwork;

	void init();

public:
	SphereLaplace(long nlat, long nlon, long isym);
	~SphereLaplace();

	/**
	 * out[i][j]
	 * i - latitude  (from the sourth pole to the north pole)
	 * j - longitude
	 **/
	void solve(double * out, const double * in, double mult = 1.0, double diag = 0.0);
	void calc(double * out, const double * in);
	void make_psi(double * psi, const double * u, const double * v);
};

#endif /* LAPL_H */

