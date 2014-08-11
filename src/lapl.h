#ifndef LAPL_H
#define LAPL_H
/* Copyright (c) 2010-2014 Alexey Ozeritsky
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
#include "spherepack.h"
#include "operator.h"

/**
 * Solve Laplace equation and calculate Laplace operator.
 */
class SPAPI SphereLaplace : public SphereOperator
{
public:
	SphereLaplace (const SphereOperator & op);
	~SphereLaplace();

	/**
	 * Solve the equation.
	 \f[
	   m \Delta u + d u = f
	 \f]
	 * out[i][j]
	 * i - latitude  (from the sourth pole to the north pole)
	 * j - longitude
	 * @param out  -- result function
	 * @param in   -- right part
	 * @param m    -- Laplace multiplier
	 * @param d    -- coefficient
	 **/
	void solve (double * out, const double * in, double m = 1.0, double d = 0.0);

	/**
	 * Solve the equation.
	 \f[
	 mu_2 \Delta^3 + mu \Delta u + d u = f
	 \f]
	 * out[i][j]
	 * i - latitude  (from the sourth pole to the north pole)
	 * j - longitude
	 * @param out  -- result function
	 * @param in   -- right part
	 * @param m    -- Laplace multiplier
	 * @param d    -- coefficient
	 **/
	void solve_l3 (double * out, const double * in, double mu2, double mu, double d);


	/**
	 * Calculate Laplace operator.
	 \f[
	  out = \Delta in
	 \f]
	 * @param out -- result function
	 * @param in  -- input function
	 */
	void calc (double * out, const double * in);
};

extern "C" {
	void SPAPI sp_lapl(double * out, const double * in, long long * nlat, long long * nlon);
}

#endif /* LAPL_H */

