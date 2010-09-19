#ifndef STATISTICS_H
#define STATISTICS_H
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

template < typename T >
class ExpectedValue
{
	typedef std::vector < T > array_t;
	int n;
	int i;
	array_t cur;

public:
	ExpectedValue(int n): n(n), i(0), cur(n) {}

	// M = sum ai / n
	// M0 = ai
	// Mi = ((n-1) sum a(i-1) / (n+1) + ai)/n = ((n-1)M(i-1) + ai)/n
	void accumulate(const array_t & value)
	{
		i += 1;
		for (int i0 = 0; i0 < n; ++i0) {
			cur[i0] = ((i - 1) * cur[i0] + value[i0]) / (double) i;
		}
	}

	array_t current()
	{
		return cur;
	}
};

template < typename T >
class Variance
{
	typedef std::vector < T > array_t;
	int n;
	ExpectedValue < T > m2x;
	ExpectedValue < T > mx2;

public:
	Variance(int n): n(n), m2x(n), mx2(n) {}

	void accumulate(const array_t & value)
	{
		array_t tmp(n);
		for (int i0 = 0; i0 < n; ++i0) {
			tmp[i0] = value[i0] * value[i0];
		}
		m2x.accumulate(value);
		mx2.accumulate(tmp);
	}

	array_t current()
	{
		array_t m2x_cur = m2x.current();
		array_t mx2_cur = mx2.current();
		array_t cur(n);
		for (int i0 = 0; i0 < n; ++i0) {
			cur[i0] = mx2_cur[i0] - m2x_cur[i0] * m2x_cur[i0];
		}
		return cur;
	}

	array_t m_current()
	{
		return m2x.current();
	}
};

#endif /* STATISTICS_H */
