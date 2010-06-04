#ifndef BARVORTEX_H
#define BARVORTEX_H

#include <vector>

#include "lapl.h"
#include "jac.h"
#include "norm.h"
#include "solver.h"

class SphereBaroclin: public SphereNorm
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
		double sigma1;
		double mu;
		double mu1;
		double k1;
		double k2;
		double theta;
		double alpha;

		rp_t cor;
		rp_t rp1;
		rp_t rp2;

		double * cor2;
		double * rp12;
		double * rp22;

		Conf(): cor(0), rp1(0), rp2(0), cor2(0), rp12(0), rp22(0) {}
	};

private:
	Conf conf;
	SphereOperator op;
	SphereLaplace lapl;
	SphereJacobian jac;
	linal::Solver < double > A;

	array_t lh;

	void build_matrix();

public:

	SphereBaroclin(const Conf & conf);
	~SphereBaroclin();

	// u = &in[0]
	// v = &in[nlat * nlon]

	void S_step(double * out, const double * in, double t);

	void L_step(double * out, const double * in, const double * z);
	void LT_step(double * out, const double * in, const double * z);
	void L_1_step(double * out, const double * in, const double * z);

	void p2u(double * u, const double * p);
	void u2p(double * p, const double * u);
};

#endif /* BARVORTEX_H */

