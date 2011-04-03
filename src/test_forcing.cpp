#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <string>

#include "linal.h"
#include "config.h"
#include "vorticity.h"
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
	vector < vector < double > > movthly_u;
	vector < vector < double > > monthly_v;

	vector < double > avg_u;
	vector < double > avg_v;
};

void load_real_data(const string & prefix, int nlat, int nlon)
{
	char ubuf[1024];
	char vbuf[1024];
	int i = 0;

	vector < double > tmp1;
	vector < double > tmp2;

	while (1) {
		int un1, un2;
		int vn1, vn2;

		snprintf(ubuf, sizeof(ubuf), "%s/u_%02d.txt", prefix.c_str(), i);
		ubuf[sizeof(ubuf)-1] = 0;
		snprintf(vbuf, sizeof(vbuf), "%s/v_%02d.txt", prefix.c_str(), i);
		vbuf[sizeof(vbuf)-1] = 0;

		tmp1.clear();
		tmp2.clear();

		mat_load(ubuf, tmp1, &un1, &un2);
		mat_load(vbuf, tmp2, &vn1, &vn2);

		if (un1 != nlat || un2 != nlon) {
			tmp1.clear();
		}

		if (vn1 != nlat || vn2 != nlon) {
			tmp2.clear();
		}

		if (tmp1.empty() || tmp2.empty()) {
			break;
		}
	}
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

