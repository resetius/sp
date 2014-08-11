#ifndef BARVORTEX_H
#define BARVORTEX_H
/* Copyright (c) 2010 Alexey Ozeritsky
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <vector>

#include "lapl.h"
#include "jac.h"
#include "norm.h"
#include "solver.h"

/**
 * Solve the two-dimensional baroclinic atmosphere equations.
 \f{eqnarray*}
  \frac{\partial \Delta u_1}{\partial t} + J(u_1, \Delta u_1 + l + h)
  + J(u_2, \Delta u_2) + \frac{\sigma}{2} \Delta (u_1 - u_2)
  - \mu \Delta^2 u_1 &=& f(\phi, \lambda)\\
  \frac{\partial \Delta u_2}{\partial t} + J(u_1, \Delta u_2)
  + J(u_2, \Delta u_1 + l + h) + \frac{\sigma}{2} \Delta (u_1 + u_2)
  - \mu \Delta^2 u_2
    &-& \\
	- \alpha^2 (\frac{\partial u_2}{\partial t} + J(u_1, u_2)
	- \mu_1 \Delta u_2
	+ \sigma_1 u_2 + g(\phi, \lambda)) &=& 0,
 \f}
 * where \f$J(\cdot,\cdot)\f$ is spherical jacobian and \f$\Delta\f$ is spherical 
 * Laplace operator.
 * @see SphereJacobian, SphereLaplace
 */
class SPAPI SphereBaroclin : public SphereNorm
{
public:
	typedef std::vector < double > array_t;

	struct Conf
	{
		typedef double (*rp_t ) (double phi, double lambda, 
			double t, const Conf * conf );

		long nlat;    ///< latitude
		long nlon;    ///< longitude

		double tau;   ///< time step
		double sigma; ///< \f$sigma\f$
		double sigma1;///< \f$sigma_1\f$
		double mu;    ///< \f$theta\f$
		double mu1;   ///< \f$mu_1\f$
		double k1;    ///< \f$k_1\f$
		double k2;    ///< \f$k_2\f$

		/**
		 * Time discretization scheme parameter \f$\theta\f$.
		 * Use 0.5 for Crank-Nicolson.
		 */
		double theta;
		double alpha; ///<\f$\alpha\f$

		/**
		 * Coriolis.
		 * l+h: coriolis function plus orographic function.
		 */
		rp_t cor;
		rp_t rp1; ///<\f$f\f$
		rp_t rp2; ///<\f$g\f$

		/**
		 * Coriolis.
		 * l+h: coriolis function plus orographic function.
		 */
		double * cor2;
		double * rp12; ///<\f$f\f$
		double * rp22; ///<\f$g\f$

		Conf(): cor(0), rp1(0), rp2(0), cor2(0), rp12(0), rp22(0) {}
	};

private:
	Conf conf;
	SphereOperator op;
	SphereLaplace lapl;
	SphereJacobian jac;
	linal::Solver < double > A;
	linal::Solver < double > A1; //inverted

	array_t lh;

	void build_matrix();
	void build_inverted_matrix();

public:

	SphereBaroclin(const Conf & conf);
	~SphereBaroclin();


	/**
	 * Solve the two-dimensional baroclinic atmosphere equations.
 \f{eqnarray*}
  \frac{\partial \Delta u_1}{\partial t} + J(u_1, \Delta u_1 + l + h)
  + J(u_2, \Delta u_2) + \frac{\sigma}{2} \Delta (u_1 - u_2)
  - \mu \Delta^2 u_1 &=& f(\phi, \lambda)\\
  \frac{\partial \Delta u_2}{\partial t} + J(u_1, \Delta u_2)
  + J(u_2, \Delta u_1 + l + h) + \frac{\sigma}{2} \Delta (u_1 + u_2)
  - \mu \Delta^2 u_2
    &-&\\
	- \alpha^2 (\frac{\partial u_2}{\partial t} + J(u_1, u_2)
	- \mu_1 \Delta u_2
	+ \sigma_1 u_2 + g(\phi, \lambda)) &=& 0,
 \f}
	 * @param u1_o  - output vector
	 * @param u2_o  - output vector
	 * @param u1 - input vector (previous time step)
	 * @param u2 - input vector (previous time step)
	 * @param t   - time
	 */
	void S_step(double * u1_o, double * u2_o, const double * u1, const double * u2, double t);

	/**
	 * Provided for convince.
	 * u1 = &in[0]
	 * u2 = &in[nlat * nlon]
	 * @param out - result
	 * @param in - input
	 * @param t - time
	 */
	void S_step(double * out, const double * in, double t);

	/**
	 * Solve the linearized two-dimensional baroclinic atmosphere equations in a neibourhood of point (z1, z2)
 \f{eqnarray*}
  \frac{\partial \Delta u_1}{\partial t} + J(u_1, \Delta z_1 + l + h) + J(z_1, \Delta u_1)
  + J(z_2, \Delta u_2) + J(u_2, \Delta z_2) +
  \frac{\sigma}{2} \Delta (u_1 - u_2)
  - \mu \Delta^2 u_1 &=& 0\\
  \frac{\partial \Delta u_2}{\partial t} + J(u_1, \Delta z_2) + J(z_1, \Delta u_2)
  + J(u_2, \Delta z_1 + l + h) + J(z_2, \Delta u_1)
  + \frac{\sigma}{2} \Delta (u_1 + u_2)
  - \mu \Delta^2 u_2
    &-&\\
	- \alpha^2 (\frac{\partial u_2}{\partial t} 
    + J(z_1, u_2) + J(u_1, z_2)
	- \mu_1 \Delta u_2
	+ \sigma_1 u_2) &=& 0,
 \f}
 
	 * @param u1_o  - output vector
	 * @param u2_o  - output vector
	 * @param u1 - input vector (previous time step)
	 * @param u2 - input vector (previous time step)
	 * @param z1  - vector z1
	 * @param z2  - vector z2
	 * @param t   - time
	 */
	void L_step(double * u1_o, double * u2_o, 
	            const double * u1, const double * u2, 
	            const double * z1, const double * z2);

	/**
	 * Solve the linearized two-dimensional baroclinic atmosphere equations in a neibourhood of point (z1, z2).
	 * Provided for convince.
	 * u1 = &in[0]
	 * u2 = &in[nlat * nlon]
	 * @param out - result
	 * @param in - input
	 * @param t - time
	 */
	void L_step(double * out, const double * in, const double * z);

	/**
	 * Solve the invert linearized two-dimensional baroclinic atmosphere equations in a neibourhood of point (z1, z2).
	 * @param u1_o  - output vector
	 * @param u2_o  - output vector
	 * @param u1 - input vector (previous time step)
	 * @param u2 - input vector (previous time step)
	 * @param z1  - vector z1
	 * @param z2  - vector z2
	 */
	void L_1_step(double * u1_o, double * u2_o, 
	              const double * u1, const double * u2, 
	              const double * z1, const double * z2);
	/**
	 * Solve the invert linearized two-dimensional baroclinic atmosphere equations in a neibourhood of point (z1, z2).
	 * Provided for convince.
	 * u1 = &in[0]
	 * u2 = &in[nlat * nlon]
	 * @param out - result
	 * @param in - input
	 * @param t - time
	 */
	void L_1_step(double * out, const double * in, const double * z);

	/**
	 * Add boundary conditions.
	 * Just copy p to u. Provided for convince.
	 */
	void p2u(double * u, const double * p);

	/**
	 * Remove boundary conditions.
	 * Just copy p to u. Provided for convince.
	 */
	void u2p(double * p, const double * u);
};

#endif /* BARVORTEX_H */

