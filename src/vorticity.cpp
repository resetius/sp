#include "vorticity.h"
#include "spherepack.h"

using namespace std;

SphereVorticity::SphereVorticity(long nlat, long nlon):
	nlat(nlat), nlon(nlon)
{
}

void SphereVorticity::calc(double * div, const double * u, const double *v)
{
	// use vrtec_ here
}
