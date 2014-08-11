#ifndef NORM_H
#define NORM_H
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

#undef min

#ifdef WIN32
#if defined(SPLIBRARY_EXPORT)
#define SPAPI   __declspec( dllexport ) 
#else
#define SPAPI   __declspec( dllimport ) 
#endif
#endif

/**
 * Spherical norm implementation.
 */
class SPAPI SphereNorm
{
	long nlat;
	long nlon;
	std::vector < double > cosi;

public:
	/**
	 * Constructor.
	 * @param nlat - latitude
	 * @param nlon - longitude
	 */
	SphereNorm(long nlat, long nlon);

	/**
	 * Scalar product of two functions.
	 * @param u
	 * @param v
	 */
	double scalar(const double * u, const double * v);

	/**
	 * the Distance between two functions.
	 * @param u
	 * @param v
	 */
	double dist(const double * u, const double * v);

	/**
	 * the norm of functions.
	 * @param u
	 */
	double norm(const double * u);
};

#endif /* NORM_H */
