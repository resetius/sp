#ifndef SPHERE_OP_H
#define SPHERE_OP_H

/**
 * @file 
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 *
 * @page License
 * @section LICENSE
 *
 * @verbatim
  Copyright (c) 2010 Alexey Ozeritsky
  All rights reserved.
 
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
  1. Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
  2. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  3. The name of the author may not be used to endorse or promote products
     derived from this software without specific prior written permission
 
  THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
  NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
  THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  @endverbatim
 *
 * @mainpage SP Documentation
 * @section into_sec Introduction
 * This library is to solve partial differential equations on sphere.
 * It extends SpherePack library and make it easy to use in C++.
 * It also implements some of the common partial differential equations
 * such as Laplace equation, Chafe-Infante equation, 
 * Barotropic vorticity equation and two-dimensional baroclinic atmosphere equations.
 *
 * @section ex Usage examples
    -# @ref test_lapl.cpp "Laplace equation on the sphere"
  \f{eqnarray*}
  \Delta \psi &=& f(\varphi, \lambda) \\
  \Delta \psi &=& \frac{1}{cos\varphi}\frac{\partial}{\partial\varphi}cos(\varphi)\frac{\partial}{\partial\varphi}\psi+
  \frac{1}{cos^2\varphi}\frac{\partial^2}{\partial\lambda^2}\psi\\
  \f}
    -# @ref test_chafe.cpp "Chafe-Infante equation on the sphere"
  \f{eqnarray*}
  \frac{du}{dt} &=& \mu \Delta u - \sigma u + f (u) \\
  u(x,y,t)|_{t=0} &=& u_0 \\
  \f}
    -# @ref test_barvortex.cpp "The Barotropic vorticity equation on the sphere"
  \f{eqnarray*}
  \frac{\partial \Delta \psi}{\partial t} + J(\psi, \Delta \psi) 
    + J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi &=& f(\varphi, \lambda) \\
	\psi|_{t=0}=\psi_0
  \f}
    -# @ref test_baroclin.cpp  "The two-dimensional baroclinic atmosphere equations on the sphere"
 \f{eqnarray*}
  \frac{\partial \Delta u_1}{\partial t} + J(u_1, \Delta u_1 + l + h)
  + J(u_2, \Delta u_2) + \frac{\sigma}{2} \Delta (u_1 - u_2)
  - \mu \Delta^2 u_1 &=& f(\phi, \lambda)\\
  \frac{\partial \Delta u_2}{\partial t} + J(u_1, \Delta u_2)
  + J(u_2, \Delta u_1 + l + h) + \frac{\sigma}{2} \Delta (u_1 + u_2)
  - \mu \Delta^2 u_2
    &-&\\
	-\alpha^2 (\frac{\partial u_2}{\partial t} + J(u_1, u_2)
	- \mu_1 \Delta u_2
	+ \sigma_1 u_2 + g(\phi, \lambda)) &=& 0,\\
	u_1|_{t=0}&=&u_{10}\\
	u_2|_{t=0}&=&u_{20}\\
 \f}
 * @page Build
 * @section build_sec Build
 * SP uses Linal library method to solve linear equations. 
 * And Spherepack library (my C port!) for Spherical Harmonics method.
 * So you must put those libraries into SP root directory in the 
 * following way:
 * @verbatim
  sp_directory/
  sp_directory/linal
  sp_directory/spherepack_c
 * 
 * @subsection Unix-like
 * @verbatim
  mkdir build-directory
  cd build-directory
  cmake -DCMAKE_BUILD_TYPE=Debug path-to-sources   # for Debug build
  cmake -DCMAKE_BUILD_TYPE=Release path-to-sources # for Release build
  make
  @endverbatim
 * @subsection Windows
 * @verbatim
  mkdir build-directory
  cd build-directory
  cmake -G "Visual Studio 2009" #place your version of Visual Studio here
  @endverbatim
 */

/**
 * 
  @example test_lapl.cpp
  Laplace equation on a sphere
  \f{eqnarray*}
  \Delta \psi &=& f(\varphi, \lambda) \\
  \Delta \psi &=& \frac{1}{cos\varphi}\frac{\partial}{\partial\varphi}cos(\varphi)\frac{\partial}{\partial\varphi}\psi+
  \frac{1}{cos^2\varphi}\frac{\partial^2}{\partial\lambda^2}\psi\\
  \f}
  @example test_chafe.cpp
  Chafe-Infante equation on a sphere
  \f{eqnarray*}
  \frac{du}{dt} &=& \mu \Delta u - \sigma u + f (u) \\
  u(x,y,t)|_{t=0} &=& u_0 \\
  \f}
  @example test_barvortex.cpp
  the Barotropic vorticity equation
  \f[
  \frac{\partial \Delta \psi}{\partial t} + J(\psi, \Delta \psi) 
    + J(\psi, l + h) + \sigma \Delta \psi - \mu \Delta^2 \psi = f(\varphi, \lambda)
  \f]
  @example test_baroclin.cpp
  The two-dimensional baroclinic atmosphere equations
 \f{eqnarray*}
  \frac{\partial \Delta u_1}{\partial t} + J(u_1, \Delta u_1 + l + h)
  + J(u_2, \Delta u_2) + \frac{\sigma}{2} \Delta (u_1 - u_2)
  - \mu \Delta^2 u_1 &=& f(\phi, \lambda)\\
  \frac{\partial \Delta u_2}{\partial t} + J(u_1, \Delta u_2)
  + J(u_2, \Delta u_1 + l + h) + \frac{\sigma}{2} \Delta (u_1 + u_2)
  - \mu \Delta^2 u_2
    &-&\\
	-\alpha^2 (\frac{\partial u_2}{\partial t} + J(u_1, u_2)
	- \mu_1 \Delta u_2
	+ \sigma_1 u_2 + g(\phi, \lambda)) &=& 0,
 \f}
  */

#include <vector>
#include "norm.h"

/**
 * The base class for all of spherical operators.
 * It allocates memory for Spherepack and implements some usefull methods.
 */
class SphereOperator: public SphereNorm
{
	struct Data;
	Data * d;

protected:
	typedef std::vector < double > array_t;

	long isym;
	long nlat;
	long nlon;
	long mdab;
	long mdb;
	long mdc;

	long slsave;
	array_t & swsave;

	long sldwork;
	array_t & sdwork;

	long vlsave;
	array_t & vwsave;

	long vldwork;
	array_t & vdwork;

	long islsave;
	array_t & iswsave;

	long isldwork;
	array_t & isdwork;

	long ivlsave;
	array_t & ivwsave;

	long ivldwork;
	array_t & ivdwork;

	long lwork;
	array_t & work;

public:
	/**
	 * Constructor.
	 * @param nlat - latitude
	 * @param nlon - longitude
	 * @param isym - unused yet
	 */
	SphereOperator(long nlat, long nlon, long isym);
	SphereOperator(const SphereOperator & op);
	~SphereOperator();

	/**
	 * Converts function to coefficients.
	 * @param k -- array 2*nlat*nlat
	 * @param f -- array nlat*nlon
	 */
	void func2koef(double * k, const double * f);

	/**
	 * Converts coefficients to function.
	 * @param k -- array 2*nlat*nlat
	 * @param f -- array nlat*nlon
	 */
	void koef2func(double * f, const double * k);

private:
	SphereOperator & operator = (const SphereOperator & op);
};

#endif /* SPHERE_OP_H */

