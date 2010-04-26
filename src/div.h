#ifndef DIV_H
#define DIV_H

#include <vector>

class SphereDiv
{
	typedef std::vector < double > array_t;

	long nlat;
	long nlon;

public:
	SphereDiv(long nlat, long nlon);
	~SphereDiv();

	void calc(double * div, const double * u, const double *v);
};

#endif /* DIV_H */
