#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>

#include "barvortex.h"
#include "utils.h"
#include "statistics.h"

#ifdef max
#undef max
#endif

using namespace std;

static inline double max (double a, double b)
{
	return (a > b) ? a : b;
}

double ans (double x, double y, double t)
{
	return x*sin (y + t) *ipow (cos (x), 4);
}

double zero_coriolis (double phi, double lambda)
{
	return 0.0;
}

double f (double x, double y, double t, SphereBarvortexConf * conf)
{
	double mu    = conf->mu;
	double sigma = conf->sigma;
	return 390*mu*sin(y+t)*x*ipow(cos(x),2)
		+147*mu*sin(y+t)*sin(x)*cos(x)-
		400*mu*sin(y+t)*x*
		ipow(cos(x),4)-360*mu*sin(y+t)*
		sin(x)*ipow(cos(x),3)-20*sigma*sin(y+t)*
		ipow(cos(x),4)*x+15*sigma*sin(y+t)*
		ipow(cos(x),2)*x-9*sigma*sin(y+t)*
		ipow(cos(x),3)*sin(x)+30*cos(y+t)*
		ipow(cos(x),4)*sin(y+t)*x*x*sin(x)+9*cos(y+t)*
		ipow(cos(x),6)*sin(y+t)*sin(x)-45*mu*sin(y+t)*x-
		9*cos(y+t)*ipow(cos(x),5)*sin(y+t)*x-20*x*cos(y+t)*
		ipow(cos(x),4)-9*cos(y+t)*
		ipow(cos(x),3)*sin(x)+15*x*cos(y+t)*
		ipow(cos(x),2);
}

void solve()
{
	long nlat = 19;
	long nlon = 36;

	SphereBarvortexConf conf;
	conf.nlat     = nlat;
	conf.nlon     = nlon;
	conf.mu       = 8e-5;
	conf.sigma    = 1.6e-2;
	conf.tau      = 0.001;
	conf.theta    = 0.5;
	conf.k1       = 1.0;
	conf.k2       = 1.0;
	conf.rp       = f;
	conf.coriolis = zero_coriolis;

	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;
	double t = 0;

	int i, j, it = 0;

	vector < double > u (nlat * nlon);
	vector < double > v (nlat * nlon);
	vector < double > r (nlat * nlon);

	SphereBarvortex bv (conf);

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
		bv.S_step (&u[0], &r[0], t);
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
					//fprintf(stderr, "%.16lf %.16lf \n", u[i * nlon + j], v[i * nlon + j]);
				}
			}

			nev1 = bv.dist(&u[0], &v[0]);

			fprintf (stderr, "time=%lf nev1=%.16le \n", t, nev1);
		}
		r.swap(u);

		it += 1;
	}
}

double test_coriolis (double phi, double lambda)
{
	return 2 * sin(phi) + 0.5 * cos(2 * lambda) * ipow(sin(2 * phi), 2);
}

double kornev1_rp_(double phi, double lambda, double t, SphereBarvortexConf * conf)
{
	double x   = phi;

	double pt1 = -0.5 * (sin(x)*M_PI*x-2*sin(x)*x*x);
	if (fabs(pt1) > 1e-14) {
		pt1 /= cos(x);
	}

	double omg = 2.*M_PI/24./60./60.; // ?
	double TE  = 1./omg;
	double U0  = 6.371e+6/TE;

	double pt2 = -0.5*(-M_PI+4*x);
	double f = - 16.0 / M_PI / M_PI * 30.0 / U0 * (pt1 + pt2);
	return conf->sigma * f;
}

void run_test()
{
	long nlat = 19;
	long nlon = 36;

	SphereBarvortexConf conf;
	conf.nlat     = nlat;
	conf.nlon     = nlon;
	conf.mu       = 6.77e-5;
	conf.sigma    = 1.14e-2;
	conf.tau      = 2 * M_PI / 12.;
	conf.theta    = 0.5;
	conf.k1       = 1.0;
	conf.k2       = 1.0;
	conf.rp       = kornev1_rp_;
	conf.coriolis = test_coriolis;

	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;
	double t = 0;
	double T = 2 * M_PI * 1000.0;

	int i, j, it = 0;

	vector < double > u (nlat * nlon);
	vector < double > v (nlat * nlon);
	vector < double > r (nlat * nlon);

	SphereBarvortex bv (conf);

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

	Variance < double > var(u.size());

	while (t < T)
	{
		bv.S_step (&u[0], &r[0], t);
		t += conf.tau;

		var.accumulate(u);

		r.swap(u);

		it += 1;
	}

	{
		FILE * f1 = fopen("m.txt", "wb");
		FILE * f2 = fopen("d.txt", "wb");
		vector < double > m = var.m_current();
		vector < double > d = var.current();
		for (i = 0; i < nlat; ++i)
		{
			for (j = 0; j < nlon; ++j)
			{
				fprintf(f1, "%.16lf ", m[i * nlon + j]);
				fprintf(f2, "%.16lf ", d[i * nlon + j]);
			}
			fprintf(f1, "\n");
			fprintf(f2, "\n");
		}
		fclose(f1);fclose(f2);
	}
}


int main (int argc, char * argv[])
{
	//solve();
	run_test();
}

