
#include <math.h>
#include "norm.h"
#include "utils.h"

SphereNorm::SphereNorm(long nlat, long nlon): nlat(nlat), nlon(nlon), cosi(nlat)
{
	double dlat = M_PI / (nlat - 1);
	for (int i = 0; i < nlat; ++i)
	{
		cosi[i] = cos(-0.5 * M_PI + i * dlat);
	}
}

double SphereNorm::scalar(const double * u, const double * v)
{
	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;

	double sum = 0.0;
	for (int i = 0; i < nlat; ++i) {
		for (int j = 0; j < nlon; ++j) {
			sum += cosi[i] * u[i * nlon + j] * v[i * nlon + j];
		}
	}
	return sum * dlat * dlon;
}

double SphereNorm::norm(const double * u)
{
	return sqrt(scalar(u, u));
}

double SphereNorm::dist(const double * u, const double * v)
{
	long n = nlat * nlon;
	std::vector < double > tmp(n);
	vec_sum1(&tmp[0], u, v, 1.0, -1.0, n);
	return norm(&tmp[0]);
}
