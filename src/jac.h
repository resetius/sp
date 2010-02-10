#ifndef SPHERE_JAC_H
#define SPHERE_JAC_H

class SphereJacobian
{
	long nlat;
	long nlon;

public:
	SphereJacobian(long nlat, long nlon);
	~SphereJacobian();

	void calc(double * out, const double * u, const double * v);
};

#endif /* SPHERE_JAC_H */

