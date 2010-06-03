#ifndef LAPL_H
#define LAPL_H

#include <vector>
#include "spherepack.h"

class SphereLaplace
{
	typedef std::vector < double > array_t;

	long nlat;
	long nlon;
	long isym;

	long slsave;
	array_t swsave;

	long sldwork;
	array_t sdwork;

	long islsave;
	array_t iswsave;

	long isldwork;
	array_t isdwork;

	long lwork;
	array_t work;
	void init();

public:
	SphereLaplace (long nlat, long nlon, long isym);
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

