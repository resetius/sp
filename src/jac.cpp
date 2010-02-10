#include "spherepack.h"
#include "jac.h"
#include "utils.h"

SphereJacobian::SphereJacobian(long nlat, long nlon): nlat(nlat), nlon(nlon)
{
}

SphereJacobian::~SphereJacobian()
{
}

void SphereJacobian::calc(double * out, const double * u, const double * v)
{
	// TODO: use gradec here
}

