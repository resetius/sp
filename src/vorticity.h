#ifndef VORTICITY_H
#define VORTICITY_H

#include <vector>

class SphereVorticity
{
	typedef std::vector < double > array_t;

	long nlat;
	long nlon;
	long mdc;

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
	SphereVorticity(long nlat, long nlon);
	~SphereVorticity();

	void calc(double * div, const double * u, const double *v);
};

#endif /* VORTICITY_H */
