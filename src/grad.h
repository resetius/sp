#ifndef GRAD_H
#define GRAD_H

#include <vector>

class SphereGrad
{
	typedef std::vector < double > array_t;

	long nlat;
	long nlon;
	long mdab;
	long mdb;

	long slsave;
	array_t swsave;

	long sldwork;
	array_t sdwork;

	long vlsave;
	array_t vwsave;

	long vldwork;
	array_t vdwork;

	long islsave;
	array_t iswsave;

	long isldwork;
	array_t isdwork;

	long ivlsave;
	array_t ivwsave;

	long ivldwork;
	array_t ivdwork;

	long lwork;
	array_t work;

public:
	SphereGrad(long nlat, long lon);
	~SphereGrad();

	void calc (double * duphi, double * dulambda, const double * u);
	void solve (double * u, const double * duphi, const double * dulambda);
};

#endif /* GRAD_H */
