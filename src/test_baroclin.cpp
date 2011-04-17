/**/

#include <string>
#include <vector>
#include <math.h>

#include "baroclin.h"
#include "grad.h"
#include "vorticity.h"
#include "statistics.h"
#include "linal.h"

using namespace std;
using namespace linal;

#ifdef WIN32
#define snprintf _snprintf
#endif

double u1_t (double x, double y, double t)
{
	return x*sin(y+t)*ipow(cos(x),4);
}

double u2_t (double x, double y, double t)
{
	return x*cos(y+t)*ipow(cos(x),4);
}

double rp_f1(double x, double y, double t, const SphereBaroclin::Conf * conf)
{
	double sigma = conf->sigma;
	double mu    = conf->mu;
	double sigma1= conf->sigma;
	double mu1   = conf->mu;
	double alpha = conf->alpha;
	
	return -45*mu*sin(y+t)*x-(9./2.)*sigma*
		ipow(cos(x),3)*sin(y+t)*sin(x)+(9./2.)*sigma*
		ipow(cos(x),3)*sin(x)*cos(y+t)-10*sigma*
		ipow(cos(x),4)*x*sin(y+t)+10*sigma*
		ipow(cos(x),4)*x*cos(y+t)+(15./2.)*sigma*
		ipow(cos(x),2)*x*sin(y+t)-(15./2.)*sigma*
		ipow(cos(x),2)*x*cos(y+t)-360*mu*sin(y+t)*sin(x)*
		ipow(cos(x),3)+147*mu*sin(y+t)*sin(x)*cos(x)-400*mu*sin(y+t)*x*
		ipow(cos(x),4)+390*mu*sin(y+t)*x*
		ipow(cos(x),2)-20*x*cos(y+t)*
		ipow(cos(x),4)-9*cos(y+t)*
		ipow(cos(x),3)*sin(x)+15*x*cos(y+t)*
		ipow(cos(x),2);
}

double rp_g1(double x, double y, double t, const SphereBaroclin::Conf * conf)
{
	double sigma = conf->sigma;
	double mu    = conf->mu;
	double sigma1= conf->sigma;
	double mu1   = conf->mu;
	double alpha = conf->alpha;

	double alpha2 = alpha * alpha;
	double r = 20*x*sin(y+t)*
		ipow(cos(x),4)+9*sin(y+t)*
		ipow(cos(x),3)*sin(x)-15*x*sin(y+t)*
		ipow(cos(x),2)+390*mu*cos(y+t)*x*
		ipow(cos(x),2)-10*sigma*
		ipow(cos(x),4)*x*cos(y+t)-(9./2.)*sigma*
		ipow(cos(x),3)*sin(y+t)*sin(x)-alpha2*
		ipow(cos(x),7)*x-45*mu*cos(y+t)*x+18*
		ipow(cos(x),6)*
		ipow(cos(y+t),2)*sin(x)-9*
		ipow(cos(x),6)*sin(x)+9*x*
		ipow(cos(x),5)-(9./2.)*sigma*
		ipow(cos(x),3)*sin(x)*cos(y+t)-30*x*x*
		ipow(cos(x),4)*sin(x)-18*x*
		ipow(cos(x),5)*
		ipow(cos(y+t),2)+alpha2*
		ipow(cos(x),4)*x*sin(y+t)+(15./2.)*sigma*
		ipow(cos(x),2)*x*sin(y+t)+(15./2.)*sigma*
		ipow(cos(x),2)*x*cos(y+t)-400*mu*cos(y+t)*x*
		ipow(cos(x),4)+147*mu*cos(y+t)*sin(x)*cos(x)+60*x*x*
		ipow(cos(x),4)*
		ipow(cos(y+t),2)*sin(x)-360*mu*cos(y+t)*sin(x)*
		ipow(cos(x),3)-10*sigma*
		ipow(cos(x),4)*x*sin(y+t)+4*alpha2*
		ipow(cos(x),6)*x*x*sin(x)-9*alpha2*
		ipow(cos(x),3)*mu1*cos(y+t)*sin(x)-20*alpha2*
		ipow(cos(x),4)*mu1*cos(y+t)*x+15*alpha2*
		ipow(cos(x),2)*mu1*cos(y+t)*x-alpha2*
		ipow(cos(x),4)*sigma1*x*cos(y+t);
	r /= alpha2;
	return r;
}

double rp_f2(double x, double y, double t, const SphereBaroclin::Conf * conf)
{
	double sigma = conf->sigma;
	double mu    = conf->mu;
	double sigma1= conf->sigma;
	double mu1   = conf->mu;
	double alpha = conf->alpha;

	return -9*cos(y+t)*
		ipow(cos(x),3)*sin(x)-20*x*cos(y+t)*
		ipow(cos(x),4)+15*x*cos(y+t)*
		ipow(cos(x),2)-(9./2.)*sigma*
		ipow(cos(x),3)*sin(y+t)*sin(x)+(9./2.)*sigma*
		ipow(cos(x),3)*sin(x)*cos(y+t)-10*sigma*
		ipow(cos(x),4)*x*sin(y+t)+10*sigma*
		ipow(cos(x),4)*x*cos(y+t)+(15./2.)*sigma*
		ipow(cos(x),2)*x*sin(y+t)-(15./2.)*sigma*
		ipow(cos(x),2)*x*cos(y+t)-360*mu*sin(y+t)*sin(x)*
		ipow(cos(x),3)+147*mu*sin(y+t)*sin(x)*cos(x)-400*mu*sin(y+t)*x*
		ipow(cos(x),4)+390*mu*sin(y+t)*x*
		ipow(cos(x),2)-45*mu*sin(y+t)*x;
}

double rp_g2(double x, double y, double t, SphereBaroclin::Conf * conf)
{
	double sigma = conf->sigma;
	double mu    = conf->mu;
	double sigma1= conf->sigma;
	double mu1   = conf->mu;
	double alpha = conf->alpha;

	return -15*x*sin(y+t)*
		ipow(cos(x),2)+20*x*sin(y+t)*
		ipow(cos(x),4)+9*sin(y+t)*
		ipow(cos(x),3)*sin(x)-(9./2.)*sigma*
		ipow(cos(x),3)*sin(x)*cos(y+t)-10*sigma*
		ipow(cos(x),4)*x*sin(y+t)-10*sigma*
		ipow(cos(x),4)*x*cos(y+t)-(9./2.)*sigma*
		ipow(cos(x),3)*sin(y+t)*sin(x)+(15./2.)*sigma*
		ipow(cos(x),2)*x*sin(y+t)+(15./2.)*sigma*
		ipow(cos(x),2)*x*cos(y+t)-360*mu*cos(y+t)*sin(x)*
		ipow(cos(x),3)+147*mu*cos(y+t)*sin(x)*cos(x)-400*mu*cos(y+t)*x*
		ipow(cos(x),4)+390*mu*cos(y+t)*x*
		ipow(cos(x),2)-45*mu*cos(y+t)*x;
}

double rp1(double phi, double lambda, SphereBaroclin::Conf * conf)
{
	double omg = 2.*M_PI/24./60./60.; // ?
	double T0  = 1./omg;
	double R   = 6.371e+6;
	double c   = T0*T0/R/R;
	double x   = phi;

	double pt1 = -0.5 * (sin(x)*M_PI*x-2*sin(x)*x*x);
	if (fabs(pt1) > 1e-14) {
		pt1 /= cos(x);
	}

	double pt2 = -0.5*(-M_PI+4*x);
	return -T0/R * 16.0 / M_PI / M_PI * 30.0 * (pt1 + pt2);
}

double rp2(double phi, double lambda, SphereBaroclin::Conf * conf)
{
	double omg = 2.*M_PI/24./60./60.; // ?
	double T0  = 1./omg;
	double R   = 6.371e+6;
	double c   = T0*T0/R/R;
	double x   = phi;

	double pt1 = -0.5 * (sin(x)*M_PI*x-2*sin(x)*x*x);
	if (fabs(pt1) > 1e-14) {
		pt1 /= cos(x);
	}

	double pt2 = -0.5*(-M_PI+4*x);
	return -T0/R * 16.0 / M_PI / M_PI * 30.0 * (pt1 + pt2);
}

double cor(double phi, double lambda, double, const SphereBaroclin::Conf * conf)
{
	return 2.*sin(phi) +  // l
		0.5 * cos(2*lambda)*sin(2*phi)*sin(2*phi); //h
}

double zero_cor(double phi, double lambda, double, const SphereBaroclin::Conf * conf)
{
	return 0;
}

double u0(double phi, double lambda)
{
	double omg = 2.*M_PI/24./60./60.; // ?
	double T0  = 1./omg;
	double R   = 6.371e+6;

	return -T0/R * 16.0 / M_PI / M_PI * 30.0 * 
		(M_PI/4 * phi * phi - phi * phi * phi / 3);
}

#define  pOff(i, j) ( i ) * conf.nlon + ( j )

template < typename T >
void proj(double * u, SphereBaroclin & bv, SphereBaroclin::Conf & conf, T rp, double t)
{
	double dlat = M_PI / (conf.nlat - 1);
	double dlon = 2. * M_PI / conf.nlon;

	for (int i = 0; i < conf.nlat; ++i) {
		double phi    = -0.5 * M_PI + i * dlat;
		for (int j = 0; j < conf.nlon; ++j) {
			double lambda = j * dlon;
			u[pOff(i, j)] = rp(phi, lambda, t);
		}
	}
}

void output_psi(const char * prefix, const char * suffix,
				const double * psi, long nlat, long nlon,
				double U0, double PSI0,
				SphereGrad & grad)
{
	vector < double > u (nlat * nlon);
	vector < double > v (nlat * nlon);
	vector < double > Psi (nlat * nlon);

	grad.calc(&v[0], &u[0], &psi[0]);
	vec_mult_scalar(&u[0], &u[0], -1.0, nlat * nlon);

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

bool test_convergence_baroclin()
{
	SphereBaroclin::Conf conf;
	double R   = 6.371e+6;
	double H   = 5000;
	double omg = 2.*M_PI/24./60./60.; // ?
	double T0  = 1./omg;
	conf.k1    = 1.0;
	conf.k2    = 1.0;
	conf.tau   = 0.0001;
	conf.sigma = 1.14e-2;
	conf.mu    = 6.77e-5;

	conf.sigma1= conf.sigma;
	conf.mu1   = conf.mu;

	conf.nlat  = 19;
	conf.nlon  = 36;
	conf.theta = 0.5;

	conf.alpha  = 1.0;
	conf.rp1    = rp_f1;
	conf.rp2    = rp_g1;

//	conf.alpha  = 0.0;
//	conf.rp1    = rp_f2;
//	conf.rp2    = rp_g2;

	conf.cor    = zero_cor;

	int n = conf.nlat * conf.nlon;

	double t = 0;
	double T = 0.1; // 30 * 2.0 * M_PI;;
	double nr1;
	double nr2;
	int i = 0;

	SphereBaroclin bv(conf);
	vector < double > u1(n);
	vector < double > u2(n);
	vector < double > u11(n);
	vector < double > u21(n);
	vector < double > u1r(n);
	vector < double > u2r(n);

	fprintf(stderr, "#mesh_w:%ld\n", conf.nlon);
	fprintf(stderr, "#mesh_h:%ld\n", conf.nlat);
	fprintf(stderr, "#tau:%.16lf\n", conf.tau);
	fprintf(stderr, "#sigma:%.16lf\n", conf.sigma);
	fprintf(stderr, "#mu:%.16lf\n", conf.mu);
	fprintf(stderr, "#sigma1:%.16lf\n", conf.sigma1);
	fprintf(stderr, "#mu1:%.16lf\n", conf.mu1);
	fprintf(stderr, "#alpha:%.16lf\n", conf.alpha);
	fprintf(stderr, "#k1:%.16lf\n", conf.k1);
	fprintf(stderr, "#k2:%.16lf\n", conf.k2);
	fprintf(stderr, "#theta:%.16lf\n", conf.theta);
	fprintf(stderr, "#rp:kornev1\n");
	fprintf(stderr, "#coriolis:kornev1\n");
	fprintf(stderr, "#initial:kornev1\n");
	fprintf(stderr, "#build:$Id$\n");

	proj(&u1[0], bv, conf, u1_t, 0);
	proj(&u2[0], bv, conf, u2_t, 0);

	while (t < T) {
		bv.S_step(&u11[0], &u21[0], &u1[0], &u2[0], t);
		t += conf.tau;

		proj(&u1r[0], bv, conf, u1_t, t);
		proj(&u2r[0], bv, conf, u2_t, t);

		nr1 = bv.dist(&u11[0], &u1r[0]);
		nr2 = bv.dist(&u21[0], &u2r[0]);
		fprintf(stderr, "t=%le; nr=%le; nr=%le; min=%le; max=%le;\n", 
				t, nr1, nr2, 
				vec_find_min(&u11[0], n),
				vec_find_max(&u11[0], n));

		if (isnan(nr1) || isnan(nr2)) {
			return false;
		}

		u11.swap(u1);
		u21.swap(u2);
		i ++;
	}

	return nr1 < 1e-5;
}

void real_calc(const char * srtm)
{
	long nlat = 5 * 19;
	long nlon = 5 * 36;

	SphereBaroclin::Conf conf;
	conf.nlat     = nlat;
	conf.nlon     = nlon;
	conf.mu       = 1.5e-5;
	conf.mu1      = conf.mu;
	conf.sigma    = 1.14e-2;
	conf.sigma1   = conf.sigma;
	conf.alpha    = 1.0;
	int part_of_the_day = 48;
	conf.tau      = 2 * M_PI / (double) part_of_the_day;
	conf.theta    = 0.5;
	conf.k1       = 1.0;
	conf.k2       = 1.0;
	conf.rp1      = 0;
	conf.rp2      = 0;
	conf.cor      = 0;//test_coriolis;

	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;
	double t = 0;
	double T = 2 * M_PI * 1000.0;

	int i, j, it = 0;

	vector < double > u0 (nlat * nlon);
	vector < double > v0 (nlat * nlon);

	vector < double > u1 (nlat * nlon);
	vector < double > v1 (nlat * nlon);

	vector < double > psi0 (nlat * nlon);
	vector < double > psi1 (nlat * nlon);

	vector < double > psi0_n (nlat * nlon);
	vector < double > psi1_n (nlat * nlon);

	vector < double > f (nlat * nlon);
	vector < double > cor(nlat * nlon);
	vector < double > rel(nlat * nlon);

	double nr0 = 0, nr1 = 0;
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
				u1[i * nlon + j] = (phi * (M_PI / 2. - phi) * 16 / M_PI / M_PI * 30.0 / U0);
			} else {
				u1[i * nlon + j] = (-phi * (M_PI / 2. + phi) * 16 / M_PI / M_PI * 30.0 / U0);
			}
			v1[i * nlon + j] = 0;
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

	vor.calc(&f[0], &u1[0], &v1[0]);
	vor.test();
	vec_mult_scalar(&f[0], &f[0], -1.0, nlat * nlon);

	lapl.solve(&psi0[0], &f[0]);
	vec_mult_scalar(&f[0], &f[0], conf.sigma, nlat * nlon);

	conf.rp12 = 0;
	conf.rp22 = &f[0];
	conf.cor2 = &cor[0];

	SphereBaroclin bv (conf);

	Variance < double > var0(u0.size());
	Variance < double > var1(u1.size());

	mat_print("out/cor.txt", &cor[0], nlat, nlon, "%23.16lf ");
	mat_print("out/rel.txt", &rel[0], nlat, nlon, "%23.16lf ");
	mat_print("out/rp.txt", &f[0], nlat, nlon, "%23.16lf ");

	mat_print("out/u0.txt", &u0[0], nlat, nlon, "%23.16lf ");
	mat_print("out/v0.txt", &v0[0], nlat, nlon, "%23.16lf ");

	mat_print("out/u1.txt", &u1[0], nlat, nlon, "%23.16lf ");
	mat_print("out/v1.txt", &v1[0], nlat, nlon, "%23.16lf ");

	while (t < T)
	{
		if (it % part_of_the_day == 0) {
			char buf[1024];
			nr0 = bv.norm(&psi0[0]);
			nr1 = bv.norm(&psi1[0]);
			if (isnan(nr0)) break;
			if (isnan(nr1)) break;

			fprintf(stderr, "nr0=%.16lf, nr1=%.16lf, t=%.16lf of %.16lf\n", nr0, nr1, t, T);
			snprintf(buf, 1024, "_%06d", it);

			output_psi("0_", buf, &psi0[0], nlat, nlon, U0, PSI0, grad);
			output_psi("1_", buf, &psi1[0], nlat, nlon, U0, PSI0, grad);

			vector < double > m0 = var0.m_current();
			vector < double > d0 = var0.current();

			vector < double > m1 = var1.m_current();
			vector < double > d1 = var1.current();

			output_psi("0_m_", "", &m0[0], nlat, nlon, U0, PSI0, grad);
			output_psi("0_d_", "", &d0[0], nlat, nlon, U0, PSI0, grad);

			output_psi("1_m_", "", &m1[0], nlat, nlon, U0, PSI0, grad);
			output_psi("1_d_", "", &d1[0], nlat, nlon, U0, PSI0, grad);
		}

		bv.S_step (&psi0_n[0], &psi1_n[0], &psi0[0], &psi1[0], t);
		t += conf.tau;

		var0.accumulate(psi0_n);
		var1.accumulate(psi1_n);

		psi0.swap(psi0_n);
		psi1.swap(psi1_n);

		it += 1;
	}

	if (!isnan(nr0) && !isnan(nr1)) {
		vector < double > m0 = var0.m_current();
		vector < double > d0 = var0.current();

		vector < double > m1 = var1.m_current();
		vector < double > d1 = var1.current();

		output_psi("0_m_", "", &m0[0], nlat, nlon, U0, PSI0, grad);
		output_psi("0_d_", "", &d0[0], nlat, nlon, U0, PSI0, grad);

		output_psi("1_m_", "", &m1[0], nlat, nlon, U0, PSI0, grad);
		output_psi("1_d_", "", &d1[0], nlat, nlon, U0, PSI0, grad);
	}
}

int main(int argc, char ** argv)
{
	fprintf(stderr, "#cmd:");
	for (int i = 0; i < argc; ++i) {
		fprintf(stderr, "%s ", argv[i]);
	}
	fprintf(stderr, "\n");
	//set_fpe_except();
	if (argc > 1) {
		if (!strcmp(argv[1], "test")) {
			return test_convergence_baroclin() ? 0 : -1;
		} else {
			real_calc(argv[1]);
		}
	}
	return 0;
}

