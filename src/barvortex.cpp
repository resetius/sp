#include "barvortex.h"

SphereBarvortex::SphereBarvortex(const SphereBarvortexConf & conf): 
	conf(conf), lapl(conf.nlat, conf.nlon), jac(conf.nlat, conf.nlon)
{
}

SphereBarvortex::~SphereBarvortex()
{
}

void SphereBarvortex::calc(double * out, const double * in, double t)
{
}

