#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <vector>
#include <string>

#include "linal.h"
#include "config.h"
#include "vorticity.h"
#include "div.h"
#include "grad.h"
#include "statistics.h"
#include "barvortex.h"

using namespace std;
using namespace linal;

#ifdef WIN32
#define snprintf _snprintf
#endif

void usage(const Config & config, const char * name)
{
	fprintf(stderr, "%s ...\n"
		"-c config_file\n", name);
	config.help();
	exit(1);
}

struct RealData
{
	int months;
	vector < vector < double > > monthly_u;
	vector < vector < double > > monthly_v;
	vector < vector < double > > monthly_omega;
	vector < vector < double > > monthly_psi;

	vector < double > avg_u;
	vector < double > avg_v;
	vector < double > avg_omega;
	vector < double > avg_psi;

	vector < double > dvomgcl;
};

RealData load_real_data(
	SphereVorticity & vor,
	SphereDiv & div,
	SphereLaplace & lapl,
	const string & prefix, 
	const string & u_pattern,
	const string & v_pattern,
	int start_month,
	int nlat, int nlon)
{
	char ubuf[1024];
	char vbuf[1024];
	int i = 0;

	double omg = 2.*M_PI/24./60./60.; // ?
	double TE  = 1./omg;
	double RE  = 6.371e+6;
	double PSI0 = RE * RE / TE;
	double U0  = 6.371e+6/TE;

	vector < double > u;
	vector < double > v;
	vector < double > omega(nlat * nlon);
	vector < double > psi(nlat * nlon);

	ExpectedValue < double > avg_u(nlat * nlon);
	ExpectedValue < double > avg_v(nlat * nlon);
	ExpectedValue < double > avg_omega(nlat * nlon);
	ExpectedValue < double > avg_psi(nlat * nlon);

	ExpectedValue < double > uomgcl(nlat * nlon);
	ExpectedValue < double > vomgcl(nlat * nlon);

	RealData data;

	i = start_month;
	while (1) {
		int un1, un2;
		int vn1, vn2;
		char pattern[1024];
		snprintf(pattern, sizeof(pattern), "%s/%s", prefix.c_str(), u_pattern.c_str());
		snprintf(ubuf, sizeof(ubuf), pattern, i);
		ubuf[sizeof(ubuf)-1] = 0;
		snprintf(pattern, sizeof(pattern), "%s/%s", prefix.c_str(), v_pattern.c_str());
		snprintf(vbuf, sizeof(vbuf), pattern, i);
		vbuf[sizeof(vbuf)-1] = 0;

		u.clear();
		v.clear();

		try {
			mat_load(ubuf, u, &un1, &un2);
			mat_load(vbuf, v, &vn1, &vn2);
		} catch (exception & ) {
			break;
		}

		if (un1 != nlat || un2 != nlon) {
			u.clear();
		}

		if (vn1 != nlat || vn2 != nlon) {
			v.clear();
		}

		if (u.empty() || v.empty()) {
			break;
		}

		/* normalize */
		vec_mult_scalar(&u[0], &u[0], 1.0/U0, nlat * nlon);
		vec_mult_scalar(&v[0], &u[0], 1.0/U0, nlat * nlon);

		avg_u.accumulate(u);
		avg_v.accumulate(v);

		vor.calc(&omega[0], &u[0], &v[0]);
		//vec_mult_scalar (&omega[0], &psi[0], -1.0, nlat * nlon);

		avg_omega.accumulate(omega);

		lapl.solve(&psi[0], &omega[0]);
		avg_psi.accumulate(psi);
		avg_omega.accumulate(omega);

		data.monthly_u.push_back(u);
		data.monthly_v.push_back(v);
		data.monthly_psi.push_back(psi);
		data.monthly_omega.push_back(omega);

		fprintf(stderr, "month %d loaded\n", i++);
	}

	assert(data.monthly_u.size() == data.monthly_v.size());
	assert(data.monthly_v.size() == data.monthly_psi.size());
	assert(data.monthly_psi.size() == data.monthly_omega.size());

	data.avg_u     = avg_u.current();
	data.avg_v     = avg_v.current();
	data.avg_psi   = avg_psi.current();
	data.avg_omega = avg_omega.current();

	data.months = (int)data.monthly_u.size();

	fprintf(stderr, "calculate average\n");
	for (int i = 0; i < data.months; ++i)
	{
		vector < double > u1(nlat * nlon);
		vector < double > v1(nlat * nlon);
		vector < double > omg1(nlat * nlon);

		vec_diff(&u1[0], &data.avg_u[0], &data.monthly_u[i][0], nlat * nlon);
		vec_diff(&v1[0], &data.avg_v[0], &data.monthly_v[i][0], nlat * nlon);
		vec_diff(&omg1[0], &data.avg_omega[0], &data.monthly_omega[i][0], nlat * nlon);

		vec_mult(&u1[0], &u1[0], &omg1[0], nlat * nlon);
		vec_mult(&v1[0], &v1[0], &omg1[0], nlat * nlon);

		uomgcl.accumulate(u1);
		vomgcl.accumulate(v1);
	}
	fprintf(stderr, "done\n");

	vector < double > tmp1 = uomgcl.current();
	vector < double > tmp2 = vomgcl.current();
	data.dvomgcl.resize(nlat * nlon);
	div.calc(&data.dvomgcl[0], &tmp1[0], &tmp2[0]);

	return data;
}

void load_relief(double * cor, double *rel, long nlat, long nlon, 
				 const string & fn)
{
	int full   = 1;
	int offset = 0;

	fprintf(stderr, "loading relief\n");

	FILE * f = fopen(fn.c_str(), "rb");
	if (f) {
		vector < double > vec;
		int n1, n2;
		mat_load(f, vec, &n1, &n2);
		if (n1 != nlat || n2 != nlon)
		{
			fprintf(stderr, "bad relief file format\n");
			exit(1);
		}
		memcpy(rel, &vec[0], nlat * nlon * sizeof(double));
		fclose(f);
	} else {
		fprintf(stderr, "file not found !\n");
		exit(1);
	}

	double rel_max = 0.0;
	for (int i = 0; i < nlat * nlon; ++i)
	{
		if (rel_max < rel[i]) rel_max = rel[i];
		//if (rel_max < fabs(rel[i])) rel_max = fabs(rel[i]);
	}

	double dlat, dlon, phi, lambda;

	if (full && !offset) {
		dlat = M_PI / (nlat - 1);
		dlon = 2. * M_PI / nlon;
	} else if (!full && !offset) {
		dlat = M_PI / (nlat - 1) / 2;
		dlon = 2. * M_PI / nlon;
	} else if (full && offset) {
		dlat = M_PI / (double) nlat;
		dlon = 2. * M_PI / nlon;
	} else if (!full && offset) {
		dlat = M_PI / (double)(2 * (nlat - 1) + 1);
		dlon = 2. * M_PI / nlon;
	}

	for (int i = 0; i < nlat; ++i)
	{
		for (int j = 0; j < nlon; ++j)
		{
			if (!offset) {
				phi    = -0.5 * M_PI + i * dlat;
				lambda = j * dlon;
			} else if (offset && !full) {
				phi    = i * dlat;
				lambda = j * dlon;
			} else if (offset && full) {
				phi    = (dlat - M_PI) * 0.5 + (double)i * dlat;
				lambda = j * dlon;
			}

			if (rel[i * nlon + j] > 0)
			{
				rel[i * nlon + j] = 1.0 * rel[i * nlon + j] / rel_max;
			}
			else
			{
				rel[i * nlon + j] = 0.0;
			}
			cor[i * nlon + j] = rel[i * nlon + j] + 2 * sin (phi);
		}
	}

	fprintf(stderr, "done\n");
}

void run_test(Config & c, int argc, char * argv[])
{
	string config_name = "sp_test_barvortex.ini";

	ConfigSkeleton s;

	s.data["sp"]["nlat"]   = make_pair(ConfigSkeleton::OPTIONAL, "latitude");
	s.data["sp"]["nlon"]   = make_pair(ConfigSkeleton::OPTIONAL, "longitude");
	s.data["sp"]["relief"] = make_pair(ConfigSkeleton::REQUIRED, "relief (text format)");
	s.data["sp"]["real_data"] = make_pair(ConfigSkeleton::REQUIRED, "directory with real data");
	s.data["sp"]["real_data_u_pattern"] 
		= make_pair(ConfigSkeleton::REQUIRED, "U pattern");
	s.data["sp"]["real_data_v_pattern"] 
		= make_pair(ConfigSkeleton::REQUIRED, "V pattern");
	s.data["sp"]["real_data_start_month"] 
		= make_pair(ConfigSkeleton::REQUIRED, "start month number");

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

	string real_data_dir = c.gets("sp", "real_data");
	string real_data_upattern = c.gets("sp", "real_data_u_pattern");
	string real_data_vpattern = c.gets("sp", "real_data_v_pattern");
	int real_data_start_month = c.geti("sp", "real_data_start_month");
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

	conf.mu       = c.get("sp", "mu", 1.5e-3);
	conf.sigma    = c.get("sp", "sigma", 1.14e-2);
	int part_of_the_day = c.get("sp", "part_of_the_day", 256);
	conf.tau      = 2 * M_PI / (double) part_of_the_day;
	conf.theta    = c.get("sp", "theta", 0.5);
	conf.k1       = c.get("sp", "k1", 1.0);
	conf.k2       = c.get("sp", "k2", 1.0);
	conf.rp       = 0;
	conf.cor      = 0;

	SphereOperator op(nlat, nlon, 0);
	SphereLaplace lapl(op);
	SphereVorticity vor(op);
	SphereDiv div(op);
	SphereGrad grad(op);

	long n = nlat * nlon;

	vector < double > rel(n);
	vector < double > cor(n);
	vector < double > force1(n);
	vector < double > forcen(n);
	vector < double > omg(n);
	vector < double > u(n);
	vector < double > v(n);

	vector < double > pt1(n);
	vector < double > pt2(n);

	RealData data = load_real_data(vor, div, lapl, 
			real_data_dir, 
			real_data_upattern, 
			real_data_vpattern, 
			real_data_start_month, 
			nlat, nlon);
	load_relief(&cor[0], &rel[0], nlat, nlon, relief_fn);

	conf.cor2 = &cor[0];

	// forcing1
	// pt1
	vec_sum(&omg[0], &data.avg_omega[0], &cor[0], n); 
	vec_mult(&u[0], &data.avg_u[0], &omg[0], n);
	vec_mult(&v[0], &data.avg_v[0], &omg[0], n);
	div.calc(&pt1[0], &u[0], &v[0]);
	// pt2
	lapl.calc(&pt2[0], &data.avg_omega[0]);
	vec_sum1(&pt2[0], &data.avg_omega[0], &pt2[0], conf.sigma, conf.mu, n);
	// pt3 = data.dvomgcl
	vec_sum(&force1[0], &force1[0], &pt1[0], n);
	vec_sum(&force1[0], &force1[0], &pt2[0], n);
	vec_sum(&force1[0], &force1[0], &data.dvomgcl[0], n);

	mat_print("forcing.txt", &force1[0], nlat, nlon, "%23.16lf ");

	for (int it = 0; it < 1 /*30*/; ++it) {

		double dist = 0;
		ExpectedValue < double > avg_total_psi(n);
		ExpectedValue < double > uomg(n);
		ExpectedValue < double > vomg(n);

		vector < double > tmp;

		fprintf(stderr, "iteration %d of 30\n", it);

		for (int i = 0; i < data.months; ++i) {
			// 30 days per month
			double T = 2 * M_PI * 30;
			double t = 0;

			// load initial here!
			vector < double > psi(n);
			ExpectedValue < double > avg_monthly_psi(n);

			conf.rp2 = &force1[0];
			SphereBarvortex bv(conf);

			avg_monthly_psi.accumulate(psi);
			avg_total_psi.accumulate(psi);

			while (t < T) {
				double nr;
				bv.S_step(&psi[0], &psi[0], t /* unused */);
				t += conf.tau;
				avg_monthly_psi.accumulate(psi);
				avg_total_psi.accumulate(psi);

				nr = bv.norm(&psi[0]);
				if (isnan(nr)) {
					fprintf(stderr, "NAN t=%lf of %lf\n", t, T);
					exit(1);
				}
			}

			vector < double > monthly_u(n);
			vector < double > monthly_v(n);
			vector < double > monthly_omega(n);
			vector < double > monthly_psi(n);

			monthly_psi = avg_monthly_psi.current();
			lapl.calc(&monthly_omega[0], &monthly_psi[0]);
			grad.calc(&monthly_v[0], &monthly_u[0], &monthly_psi[0]);
			vec_mult_scalar(&monthly_u[0], &monthly_u[0], -1.0, nlat * nlon);

			vector < double > u1(n);
			vector < double > v1(n);
			vector < double > omg1(n);

			vec_diff(&u1[0], &data.avg_u[0], &monthly_u[0], n);
			vec_diff(&v1[0], &data.avg_v[0], &monthly_v[0], n);
			vec_diff(&omg1[0], &data.avg_omega[0], &monthly_omega[0], n);

			vec_mult(&u1[0], &u1[0], &omg1[0], n);
			vec_mult(&v1[0], &v1[0], &omg1[0], n);

			uomg.accumulate(u1);
			vomg.accumulate(v1);
			fprintf(stderr ,"  month %d of %d done\n", i+1, data.months);
		}

		vector < double > avg_omega = avg_total_psi.current();
		lapl.calc(&avg_omega[0], &avg_omega[0]);
		vec_sum (&forcen[0], &force1[0], &data.avg_omega[0], n);
		vec_diff(&forcen[0], &forcen[0], &avg_omega[0], n);

		vector < double > tmp1 = uomg.current();
		vector < double > tmp2 = vomg.current();
		vector < double > dvomg(n);
		div.calc(&dvomg[0], &tmp1[0], &tmp2[0]);

		vec_sum (&forcen[0], &forcen[0], &data.dvomgcl[0], n);
		vec_diff(&forcen[0], &forcen[0], &dvomg[0], n);

		dist = op.dist(&forcen[0], &force1[0]);
		forcen.swap(force1);

		if (dist < 1e-7) {
			fprintf(stderr, "dist = %lf\n", dist);
		}
	}

	mat_print("forcing.txt", &force1[0], nlat, nlon, "%23.16lf ");
}

int main (int argc, char * argv[])
{
	//solve();

	// exe [relief in binary format!]
	Config config;
	try {  
		run_test(config, argc, argv);
	} catch (std::exception & e) {
		fprintf(stderr, "exception: %s\n", e.what());
		usage(config, argv[0]);
	}
}

