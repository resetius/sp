#include "div.h"
#include "spherepack.h"

using namespace std;

SphereDiv::SphereDiv(const SphereOperator & op): SphereOperator(op)
{
}

void SphereDiv::calc(double * div, const double * u, const double *v)
{
	// use divec_ here
}
