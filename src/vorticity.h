#ifndef VORTICITY_H
#define VORTICITY_H

#include <vector>

class SphereVorticity
{
	typedef std::vector < double > array_t;

	long nlat;
	long nlon;

public:
	SphereVorticity(long nlat, long nlon);
	~SphereVorticity();

	void calc(double * div, const double * u, const double *v);
};

#endif /* VORTICITY_H */
