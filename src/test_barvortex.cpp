#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>

#include "barvortex.h"
#include "grad.h"
#include "vorticity.h"
#include "statistics.h"
#include "linal.h"
#include "config.h"

#ifdef max
#undef max
#endif

#ifdef WIN32
#define snprintf _snprintf
#endif

using namespace std;
using namespace linal;

void usage(const Config & config, const char * name);

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

	while (it < 1000)
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

	if (nev1 < 3e-4) {
		exit(0);
	} else {
		exit(-1);
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

	char ubuf[1024]; char vbuf[1024]; char psibuf[1024];
	char Ubuf[1024]; char Vbuf[1024]; char Psibuf[1024];

	snprintf(ubuf, 1024,   "out/%snorm_u%s.txt", prefix, suffix);
	snprintf(vbuf, 1024,   "out/%snorm_v%s.txt", prefix, suffix);
	snprintf(psibuf, 1024, "out/%snorm_psi%s.txt", prefix, suffix);

	snprintf(Ubuf, 1024,   "out/%sorig_u%s.txt", prefix, suffix);
	snprintf(Vbuf, 1024,   "out/%sorig_v%s.txt", prefix, suffix);
	snprintf(Psibuf, 1024, "out/%sorig_psi%s.txt", prefix, suffix);

	mat_print(ubuf,   &u[0], nlat, nlon, "%23.16lf ");
	mat_print(vbuf,   &v[0], nlat, nlon, "%23.16lf ");
	mat_print(psibuf, &psi[0], nlat, nlon, "%23.16lf ");

	vec_mult_scalar(&u[0], &u[0], U0, nlon * nlat);
	vec_mult_scalar(&v[0], &v[0], U0, nlon * nlat);
	vec_mult_scalar(&Psi[0],  &psi[0],  PSI0, nlon * nlat);

	mat_print(Ubuf,   &u[0], nlat, nlon, "%23.16le ");
	mat_print(Vbuf,   &v[0], nlat, nlon, "%23.16le ");
	mat_print(Psibuf, &Psi[0], nlat, nlon, "%23.16le ");
}

void run_test(Config & c, int argc, char * argv[])
{
	string config_name = "sp_test_barvortex.ini";

	ConfigSkeleton s;

	s.data["sp"]["nlat"]   = make_pair(ConfigSkeleton::OPTIONAL, "latitude");
	s.data["sp"]["nlon"]   = make_pair(ConfigSkeleton::OPTIONAL, "longitude");
	s.data["sp"]["relief"] = make_pair(ConfigSkeleton::REQUIRED, "relief (text format)");

	c.set_skeleton(s);

	for (int i = 0; i < argc; ++i) {
		if (!strcmp(argv[i], "-c")) {
			if (i == argc - 1) {
				usage(c, argv[0]);
			}
			config_name = argv[i + 1];
		}

		if (!strcmp(argv[i], "-h")) {
			usage(c, argv[0]);
		}
	}

	c.open(config_name);
	c.rewrite(argc, argv);
	if (!c.validate()) {
		usage(c, argv[0]);
	}

	string relief_fn = c.gets("sp", "relief");
	string rp_u, rp_v;
	int real_f = c.get("sp", "real_f", 0);
	if (real_f) {
		rp_u = c.gets("sp", "rp_u");
		rp_v = c.gets("sp", "rp_v");
	}
	long nlat = c.get("sp", "nlat", 19);
	long nlon = c.get("sp", "nlon", 36);

	SphereBarvortex::Conf conf;
	conf.nlat     = nlat;
	conf.nlon     = nlon;
	//conf.mu       = 1.5e-5;

	conf.mu       = c.get("sp", "mu", 7.5e-5);
	conf.sigma    = c.get("sp", "sigma", 1.14e-2);
	int part_of_the_day = c.get("sp", "part_of_the_day", 256);
	conf.tau      = 2 * M_PI / (double) part_of_the_day;
	conf.theta    = c.get("sp", "theta", 0.5);
	conf.k1       = c.get("sp", "k1", 1.0);
	conf.k2       = c.get("sp", "k2", 1.0);
	conf.rp       = 0;
	conf.cor      = 0;//test_coriolis;

	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;
	double t = 0;
	double T = 2 * M_PI * 365 * 10.0;

	int i, j, it = 0;

	vector < double > u (nlat * nlon);
	vector < double > v (nlat * nlon);

	vector < double > r (nlat * nlon);
	vector < double > f (nlat * nlon);
	vector < double > cor(nlat * nlon);
	vector < double > rel;

	double nr = 0;
	double omg = 2.*M_PI/24./60./60.; // ?
	double TE  = 1./omg;
	double RE  = 6.371e+6;
	double PSI0 = RE * RE / TE;
	double U0  = 6.371e+6/TE;

	int n1, n2;
	mat_load(relief_fn.c_str(), rel, &n1, &n2);
	if (n1 != nlat || n2 != nlon) {
		fprintf(stderr, "relief: bad file format! %d!=%d %d!=%d\n", n1, nlat, n2, nlon);
		exit(1);
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
	SphereJacobian jac(op);
	SphereVorticity vor(op);

	//vor.calc(&f[0], &u[0], &v[0]);
	//vor.test();
	//vec_mult_scalar(&f[0], &f[0], -1.0, nlat * nlon);
	//
	lapl.make_psi(&f[0], &u[0], &v[0]);
	lapl.solve(&r[0], &f[0]);


	if (real_f) {
		vector < double > u;
		int n1, n2;
		mat_load(rp_u.c_str(), u, &n1, &n2);
		if (n1 != nlat || n2 != nlon)
		{
			fprintf(stderr, "bad file format! %d!=%d %d!=%d\n", n1, nlat, n2, nlon);
			exit(1);
		}

		vector < double > v;
		mat_load(rp_v.c_str(), v, &n1, &n2);
		if (n1 != nlat || n2 != nlon)
		{
			fprintf(stderr, "bad file format! %d!=%d %d!=%d\n", n1, nlat, n2, nlon);
			exit(1);
		}

		vector < double > psi(nlat * nlon);
		vector < double > dpsi(nlat * nlon);
		vector < double > ddpsi(nlat * nlon);
		vector < double > jac1(nlat * nlon);
		vector < double > jac2(nlat * nlon);

		lapl.make_psi(&dpsi[0], &u[0], &v[0]);
		lapl.solve(&psi[0], &dpsi[0]);
		lapl.calc(&ddpsi[0], &dpsi[0]);

		jac.calc(&jac1[0], &psi[0], &dpsi[0]);
		jac.calc(&jac2[0], &psi[0], &cor[0]);

		memset(&f[0], 0, f.size() * sizeof(double));
		vec_sum1(&f[0], &jac1[0], &jac2[0], conf.k1, conf.k2, nlat * nlon);
		vec_sum2(&f[0], &f[0], &dpsi[0], conf.sigma, nlat * nlon);
		vec_sum2(&f[0], &f[0], &ddpsi[0], -conf.mu, nlat * nlon);

		mat_print("out/u0f.txt", &u[0], nlat, nlon, "%23.16lf ");
		mat_print("out/v0f.txt", &v[0], nlat, nlon, "%23.16lf ");
	} else {
		vec_mult_scalar(&f[0], &f[0], conf.sigma, nlat * nlon);
	}

	conf.rp2  = &f[0];
	conf.cor2 = &cor[0];

	SphereBarvortex bv (conf);

	Variance < double > var(u.size());

	mat_print("out/cor.txt", &cor[0], nlat, nlon, "%23.16lf ");
	mat_print("out/rel.txt", &rel[0], nlat, nlon, "%23.16lf ");
	mat_print("out/rp.txt", &f[0], nlat, nlon, "%23.16lf ");
	mat_print("out/u0.txt", &u[0], nlat, nlon, "%23.16lf ");
	mat_print("out/v0.txt", &v[0], nlat, nlon, "%23.16lf ");

	//exit(1);

	while (t < T)
	{
		if (it % part_of_the_day == 0) {
			char buf[1024];
			nr = bv.norm(&r[0]);
			if (isnan(nr)) {
				fprintf(stderr, "NAN !\n");
				break;
			}
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

void usage(const Config & config, const char * name)
{
	fprintf(stderr, "%s ...\n"
		"-c config_file\n", name);
	config.help();
	exit(1);
}

int main (int argc, char * argv[])
{
	if (argc == 2 && !strcmp(argv[1], "test")) {
		solve();
	}

	// exe [relief in binary format!]
	Config config;
	try {  
		run_test(config, argc, argv);
	} catch (std::exception & e) {
		fprintf(stderr, "exception: %s\n", e.what());
		usage(config, argv[0]);
	}
}

