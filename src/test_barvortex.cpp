#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>

#include "barvortex.h"
#include "grad.h"
#include "vorticity.h"
#include "statistics.h"
#include "linal.h"

#ifdef max
#undef max
#endif

#ifdef WIN32
#define snprintf _snprintf
#endif

using namespace std;
using namespace linal;

static inline double max (double a, double b)
{
	return (a > b) ? a : b;
}

double ans (double x, double y, double t)
{
	return x*sin (y + t) *ipow (cos (x), 4);
}

double zero_coriolis (double phi, double lambda, double t, const SphereBarvortex::Conf * conf)
{
	return 0.0;
}

double f (double x, double y, double t, const SphereBarvortex::Conf * conf)
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

	SphereBarvortex::Conf conf;
	conf.nlat     = nlat;
	conf.nlon     = nlon;
	conf.mu       = 8e-5;
	conf.sigma    = 1.6e-2;
	conf.tau      = 0.001;
	conf.theta    = 0.5;
	conf.k1       = 1.0;
	conf.k2       = 1.0;
	conf.rp       = f;
	conf.cor      = zero_coriolis;

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

inline double sign(double k)
{
	if (k > 0) {
		return 1.0;
	} else {
		return -1.0;
	}
}

double test_coriolis (double phi, double lambda)
{
	return 2 * sin(phi) + 0.5 * cos(2 * lambda) * sign(2 * phi) * ipow(sin(2 * phi), 2);
}

void output_psi(const char * prefix, const char * suffix,
				const double * psi, long nlat, long nlon,
				double U0, double PSI0,
				SphereGrad & grad)
{
	vector < double > u (nlat * nlon);
	vector < double > v (nlat * nlon);
	vector < double > Psi (nlat * nlon);

	grad.calc(&u[0], &v[0], &psi[0]);

	vec_mult_scalar(&u[0], &u[0], -1.0, nlat * nlon);

	char ubuf[1024]; char vbuf[1024]; char psibuf[1024];
	char Ubuf[1024]; char Vbuf[1024]; char Psibuf[1024];

	snprintf(ubuf, 1024,   "out/%snorm_u%s.txt", prefix, suffix);
	snprintf(vbuf, 1024,   "out/%snorm_v%s.txt", prefix, suffix);
	snprintf(psibuf, 1024, "out/%snorm_psi%s.txt", prefix, suffix);

	snprintf(Ubuf, 1024,   "out/%sorig_u%s.txt", prefix, suffix);
	snprintf(Vbuf, 1024,   "out/%sorig_v%s.txt", prefix, suffix);
	snprintf(Psibuf, 1024, "out/%sorig_psi%s.txt", prefix, suffix);

	fprintfwmatrix(ubuf,   &u[0], nlat, nlon, "%23.16lf ");
	fprintfwmatrix(vbuf,   &v[0], nlat, nlon, "%23.16lf ");
	fprintfwmatrix(psibuf, &psi[0], nlat, nlon, "%23.16lf ");

	vec_mult_scalar(&u[0], &u[0], U0, nlon * nlat);
	vec_mult_scalar(&v[0], &v[0], U0, nlon * nlat);
	vec_mult_scalar(&Psi[0],  &psi[0],  PSI0, nlon * nlat);

	fprintfwmatrix(Ubuf,   &u[0], nlat, nlon, "%23.16le ");
	fprintfwmatrix(Vbuf,   &v[0], nlat, nlon, "%23.16le ");
	fprintfwmatrix(Psibuf, &Psi[0], nlat, nlon, "%23.16le ");
}

void run_test(const char * srtm)
{
	long nlat = 5 * 19;
	long nlon = 5 * 36;

	SphereBarvortex::Conf conf;
	conf.nlat     = nlat;
	conf.nlon     = nlon;
	conf.mu       = 1.5e-5;
	conf.sigma    = 1.14e-2;
	int part_of_the_day = 48;
	conf.tau      = 2 * M_PI / (double) part_of_the_day;
	conf.theta    = 0.5;
	conf.k1       = 1.0;
	conf.k2       = 1.0;
	conf.rp       = 0;
	conf.cor      = 0;//test_coriolis;

	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;
	double t = 0;
	double T = 2 * M_PI * 1000.0;

	int i, j, it = 0;

	vector < double > u (nlat * nlon);
	vector < double > v (nlat * nlon);

	vector < double > r (nlat * nlon);
	vector < double > f (nlat * nlon);
	vector < double > cor(nlat * nlon);
	vector < double > rel(nlat * nlon);

	double nr = 0;
	double omg = 2.*M_PI/24./60./60.; // ?
	double TE  = 1./omg;
	double RE  = 6.371e+6;
	double PSI0 = RE * RE / TE;
	double U0  = 6.371e+6/TE;
	const char * fn = srtm ? srtm : "";

	if (fn) {
		FILE * f = fopen(fn, "rb");
		if (f) {
			size_t size = nlat * nlon * sizeof(double);
			if (fread(&rel[0], 1, size, f) != size) {
				fprintf(stderr, "bad relief file format ! \n");
			}
		} else {
			fprintf(stderr, "relief file not found ! \n");
		}
		fclose(f);
	}

	double rel_max = 0.0;
	for (i = 0; i < nlat * nlon; ++i) {
		if (rel_max < rel[i]) rel_max = rel[i];
		//if (rel_max < fabs(rel[i])) rel_max = fabs(rel[i]);
	}

	for (i = 0; i < nlat; ++i)
	{
		for (j = 0; j < nlon; ++j)
		{
			double phi    = -0.5 * M_PI + i * dlat;
			double lambda = j * dlon;

			//double ff = -(M_PI / 4 * ipow(phi, 2) - fabs(ipow(phi, 3)) / 3.0) * 16.0 / M_PI / M_PI * 3.0 / U0;
			//r[i * nlon + j] = (phi > 0) ? ff : -ff;
			if (phi > 0) {
				u[i * nlon + j] = (phi * (M_PI / 2. - phi) * 16 / M_PI / M_PI * 30.0 / U0);
			} else {
				u[i * nlon + j] = (-phi * (M_PI / 2. + phi) * 16 / M_PI / M_PI * 30.0 / U0);
			}
			v[i * nlon + j] = 0;
			//cor[i * nlon + j] = conf.coriolis(phi, lambda);

			//cor[i * nlon + j] = 1000 * rel[i * nlon + j] / rel_max + 2 * sin(phi);
			//

		//	rel[i * nlon + j] = 0.5 * sign(phi) * cos(2 * lambda) * ipow(sin(2 * phi), 2);

			if (rel[i * nlon + j] > 0) {
				rel[i * nlon + j] = 1.0 * rel[i * nlon + j] / rel_max;
			} else {
				rel[i * nlon + j] = 0.0;
			}
			cor[i * nlon + j] = rel[i * nlon + j] + 2 * sin(phi);
		}
	}

	SphereOperator op(nlat, nlon, 0);
	SphereLaplace lapl(op);
	SphereGrad grad(op);
	SphereVorticity vor(op);

	vor.calc(&f[0], &u[0], &v[0]);
	vor.test();
	vec_mult_scalar(&f[0], &f[0], -1.0, nlat * nlon);
#if 0
	for (i = 0; i < nlat; ++i)
	{
		for (j = 0; j < nlon; ++j)
		{
			double phi    = -0.5 * M_PI + i * dlat;
			if (phi < 0) {
				f[i * nlon + j] *= -1.0;
			}
		}
	}
#endif
	lapl.solve(&r[0], &f[0]);
	vec_mult_scalar(&f[0], &f[0], conf.sigma, nlat * nlon);

	conf.rp2  = &f[0];
	conf.cor2 = &cor[0];

	SphereBarvortex bv (conf);

	Variance < double > var(u.size());

	fprintfwmatrix("out/cor.txt", &cor[0], nlat, nlon, "%23.16lf ");
	fprintfwmatrix("out/rel.txt", &rel[0], nlat, nlon, "%23.16lf ");
	fprintfwmatrix("out/rp.txt", &f[0], nlat, nlon, "%23.16lf ");
	fprintfwmatrix("out/u0.txt", &u[0], nlat, nlon, "%23.16lf ");
	fprintfwmatrix("out/v0.txt", &v[0], nlat, nlon, "%23.16lf ");

	//exit(1);

	while (t < T)
	{
		if (it % part_of_the_day == 0) {
			char buf[1024];
			nr = bv.norm(&r[0]);
			if (isnan(nr)) break;
			fprintf(stderr, "nr=%.16lf, t=%.16lf of %.16lf\n", nr, t, T);
			snprintf(buf, 1024, "_%06d", it);
			output_psi("", buf, &r[0], nlat, nlon, U0, PSI0, grad);

			vector < double > m = var.m_current();
			vector < double > d = var.current();

			output_psi("m_", "", &m[0], nlat, nlon, U0, PSI0, grad);
			output_psi("d_", "", &d[0], nlat, nlon, U0, PSI0, grad);
		}

		bv.S_step (&u[0], &r[0], t);
		t += conf.tau;

		var.accumulate(u);

		r.swap(u);

		it += 1;
	}

	if (!isnan(nr)) {
		vector < double > m = var.m_current();
		vector < double > d = var.current();

		output_psi("m_", "", &m[0], nlat, nlon, U0, PSI0, grad);
		output_psi("d_", "", &d[0], nlat, nlon, U0, PSI0, grad);
	}
}


int main (int argc, char * argv[])
{
	//solve();

	// exe [relief in binary format!]
	run_test(argv[1]);
}

