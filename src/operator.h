#ifndef SPHERE_OP_H
#define SPHERE_OP_H
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
#include "norm.h"

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
	SphereOperator(long nlat, long lon, long isym);
	SphereOperator(const SphereOperator & op);
	~SphereOperator();

	/**
	 * k -- array 2*nlat*nlat
	 * f -- array nlat*nlon
	 */
	void func2koef(double * k, const double * f);
	void koef2func(double * f, const double * k);

private:
	SphereOperator & operator = (const SphereOperator & op);
};

#endif /* SPHERE_OP_H */

