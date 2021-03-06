#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>

#include "chafe.h"
#include "linal.h"

#ifdef max
#undef max
#endif

using namespace std;
using namespace linal;

static inline double max (double a, double b)
{
	return (a > b) ? a : b;
}

static double ans (double x, double y, double t)
{
	return x*sin (y + t) *ipow (cos (x), 4);
}

static double f (double x, double y, double t,
          double mu, double sigma)
{
	return ipow (cos (x), 2) *
	       (x*cos (y + t) *ipow (cos (x), 2)
	        + sigma*sin (y + t) *ipow (cos (x), 2) *x + 9*mu*sin (y + t) *sin (x) *cos (x)
	        - 15*mu*sin (y + t) *x + 20*mu*sin (y + t) *x*ipow (cos (x), 2) );		
}

static bool solve()
{
	long nlat = 19;
	long nlon = 36;

	SphereChafeConf conf;
	conf.nlat  = nlat;
	conf.nlon  = nlon;
	conf.mu    = 1.0;
	conf.sigma = +70.0;
	conf.tau   = 0.01;
	conf.theta = 0.5;
	conf.rp    = f;

	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;
	double t = 0;

	int i, j, it = 0;

	vector < double > u (nlat * nlon);
	vector < double > v (nlat * nlon);
	vector < double > r (nlat * nlon);

	SphereChafe chafe (conf);

	double nev1 = 0;

	for (i = 0; i < nlat; ++i)
	{
		for (j = 0; j < nlon; ++j)
		{
			double phi    = -0.5 * M_PI + i * dlat;
			double lambda = j * dlon;

			r[i * nlon + j] = ans (phi, lambda, t);
		}
	}

	while (true)
	{
		chafe.calc (&u[0], &r[0], t);
		t += conf.tau;

		if (it % 100 == 0)
		{
			nev1 = 0.0;
			for (i = 0; i < nlat; ++i)
			{
				for (j = 0; j < nlon; ++j)
				{
					double phi    = -0.5 * M_PI + i * dlat;
					double lambda = j * dlon;

					v[i * nlon + j] = ans (phi, lambda, t);
					nev1 = max (nev1, fabs (u[i * nlon + j] - v[i * nlon + j] ) );

//					fprintf(stderr, "%.16lf %.16lf \n", u[i * nlon + j], v[i * nlon + j]);
				}
			}

			fprintf (stderr, "nev1=%.16le \n", nev1);
			break;
		}
		r.swap(u);

		it += 1;
	}

	return nev1 < 1e-3;
}

extern "C" int test_chafe (int argc, char * argv[])
{
	if (solve()) {
		return 0;
	}
	return -1;
}

