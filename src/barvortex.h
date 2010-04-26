#ifndef BARVORTEX_H
#define BARVORTEX_H

#include <vector>

#include "lapl.h"
#include "jac.h"
#include "norm.h"

class SphereBarvortex: public SphereNorm
{
public:
	typedef std::vector < double > array_t;

	struct Conf
	{
		typedef double (*rp_t ) (double phi, double lambda, 
			double t, const Conf * conf );

		long nlat;
		long nlon;

		double tau;
		double sigma;
		double mu;
		double k1;
		double k2;
		double theta;

		rp_t cor;
		rp_t rp;
		double * cor2;
		double * rp2;

		Conf(): rp(0), rp2(0), cor(0), cor2(0) {}
	};

private:
	Conf conf;
	SphereLaplace lapl;
	SphereJacobian jac;

	array_t lh;

public:

	SphereBarvortex(const Conf & conf);
	~SphereBarvortex();

	void S_step(double * out, const double * in, double t);

	void L_step(double * out, const double * in, const double * z);
	void LT_step(double * out, const double * in, const double * z);
	void L_1_step(double * out, const double * in, const double * z);
};

#endif /* BARVORTEX_H */

