#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>

#include "barvortex.h"
#include "grad.h"
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

void run_test()
{
	long nlat = 19;
	long nlon = 36;

	SphereBarvortexConf conf;
	conf.nlat     = nlat;
	conf.nlon     = nlon;
	conf.mu       = 6.77e-5;
	conf.sigma    = 1.14e-2;
	conf.tau      = M_PI / 12.;
	conf.theta    = 0.5;
	conf.k1       = 1.0;
	conf.k2       = 1.0;
	conf.rp       = 0;
	conf.coriolis = test_coriolis;

	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;
	double t = 0;
	double T = 2 * M_PI * 1000.0;

	int i, j, it = 0;

	vector < double > u (nlat * nlon);
	vector < double > U (nlat * nlon);
	vector < double > r (nlat * nlon);
	vector < double > f (nlat * nlon);

	vector < double > uu (nlat * nlon);
	vector < double > vv (nlat * nlon);

	double nr = 0;
	double omg = 2.*M_PI/24./60./60.; // ?
	double TE  = 1./omg;
	double RE  = 6.371e+6;
	double PSI0 = RE * RE / TE;
	double U0  = 6.371e+6/TE;

	for (i = 0; i < nlat; ++i)
	{
		for (j = 0; j < nlon; ++j)
		{
			double phi    = -0.5 * M_PI + i * dlat;
			double lambda = j * dlon;

			r[i * nlon + j] = -(M_PI / 4 * ipow(phi, 2) - ipow(phi, 3) / 3.0) * 16.0 / M_PI / M_PI * 3.0 / U0;
		}
	}

	SphereLaplace lapl(nlat, nlon);
	SphereGrad grad(nlat, nlon);

	lapl.calc(&f[0], &r[0]);
	vec_mult_scalar(&f[0], &f[0], conf.sigma, nlat * nlon);

	conf.rp2 = &f[0];

	SphereBarvortex bv (conf);

	Variance < double > var(u.size());

	while (t < T)
	{
		bv.S_step (&u[0], &r[0], t);
		t += conf.tau;

		nr = bv.norm(&u[0]);
		fprintf(stderr, "nr=%.16lf, t=%.16lf of %.16lf\n", nr, t, T);

		var.accumulate(u);

		if (i % 100 == 0) {
			grad.calc(&uu[0], &vv[0], &u[0]);
			vec_mult_scalar(&uu[0], &uu[0], -1.0, nlat * nlon);

			char ubuf[1024]; char vbuf[1024]; char psibuf[1024];
			char Ubuf[1024]; char Vbuf[1024]; char Psibuf[1024];

			snprintf(ubuf, 1024,   "out/norm_u_06d.txt", it);
			snprintf(vbuf, 1024,   "out/norm_v_06d.txt", it);
			snprintf(psibuf, 1024, "out/norm_psi_06d.txt", it);

			snprintf(Ubuf, 1024,   "out/orig_u_06d.txt", it);
			snprintf(Vbuf, 1024,   "out/orig_v_06d.txt", it);
			snprintf(Psibuf, 1024, "out/orig_psi_06d.txt", it);

			fprintfwmatrix(ubuf,   &uu[0], nlat, nlon, "%23.16lf ");
			fprintfwmatrix(vbuf,   &vv[0], nlat, nlon, "%23.16lf ");
			fprintfwmatrix(psibuf,  &u[0], nlat, nlon, "%23.16lf ");

			vec_mult_scalar(&uu[0], &uu[0], 1.0 / RE, nlon * nlat);
			vec_mult_scalar(&vv[0], &vv[0], 1.0 / RE, nlon * nlat);
			vec_mult_scalar(&U[0],  &u[0],  PSI0, nlon * nlat);

			fprintfwmatrix(Ubuf,   &uu[0], nlat, nlon, "%23.16lf ");
			fprintfwmatrix(Vbuf,   &vv[0], nlat, nlon, "%23.16lf ");
			fprintfwmatrix(Psibuf,  &U[0], nlat, nlon, "%23.16lf ");
		}

		r.swap(u);

		it += 1;
	}

	{
		vector < double > m = var.m_current();
		vector < double > d = var.current();

		fprintfwmatrix("m.txt", &m[0], nlat, nlon, "%23.16lf ");
		fprintfwmatrix("d.txt", &d[0], nlat, nlon, "%23.16lf ");
	}
}


int main (int argc, char * argv[])
{
	//solve();
	run_test();
}

