#ifndef GRAD_H
#define GRAD_H

#include <vector>

class SphereGrad
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
	SphereGrad(long nlat, long lon);
	~SphereGrad();

	void calc (double * duphi, double * dulambda, const double * u);
};

#endif /* GRAD_H */
