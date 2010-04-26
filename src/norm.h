#ifndef NORM_H
#define NORM_H

#include <vector>

class SphereNorm
{
	long nlat;
	long nlon;
	std::vector < double > cosi;

public:
	SphereNorm(long nlat, long nlon);

	double scalar(const double * u, const double * v);
	double dist(const double * u, const double * v);
	double norm(const double * u);
};

#endif /* NORM_H */
