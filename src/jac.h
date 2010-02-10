#ifndef SPHERE_JAC_H
#define SPHERE_JAC_H

#include <vector>

class SphereJacobian
{
	typedef std::vector < double > array_t;
	long nlat;
	long nlon;
	long mdab;

	long slsave;
	array_t swsave;

	long sldwork;
	array_t sdwork;

	long vlsave;
	array_t vwsave;

	long vldwork;
	array_t vdwork;

	long lwork;
	array_t work;

public:
	SphereJacobian(long nlat, long nlon);
	~SphereJacobian();

	void calc(double * out, const double * u, const double * v);
};

#endif /* SPHERE_JAC_H */

