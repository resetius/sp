#include "div.h"
#include "spherepack.h"

using namespace std;

SphereDiv::SphereDiv(long nlat, long nlon):
	nlat(nlat), nlon(nlon)
{
}

void SphereDiv::calc(double * div, const double * u, const double *v)
{
	// use divec_ here
}
