#include <math.h>
#include <string.h>

#include "barvortex.h"
#include "utils.h"

using namespace std;

SphereBarvortex::SphereBarvortex (const SphereBarvortexConf & conf) :
		conf (conf), lapl (conf.nlat, conf.nlon, conf.isym), 
			jac (conf.nlat, conf.nlon, conf.isym),
		lh (conf.nlat * conf.nlon), cosi(conf.nlat)
{
	long nlat = conf.nlat;
	long nlon = conf.nlon;
	double dlat = M_PI / (nlat-1);
	double dlon = 2. * M_PI /nlon;


	if (conf.coriolis)
	{
		for (int i = 0; i < nlat; ++i)
		{
			double phi    = -0.5 * M_PI + i * dlat;
			for (int j = 0; j < nlon; ++j)
			{
				double lambda = j * dlon;
				if (conf.coriolis) {
					lh[i * nlon + j] = conf.coriolis(phi, lambda);
				} else if (conf.coriolis2) {
					lh[i * nlon + j] = conf.coriolis2[i * nlon + j];
				}
			}
		}
	}

	for (int i = 0; i < nlat; ++i)
	{
		cosi[i] = cos(-0.5 * M_PI + i * dlat);
	}
}

SphereBarvortex::~SphereBarvortex()
{
}

double SphereBarvortex::scalar(const double * u, const double * v)
{
	long nlat = conf.nlat;
	long nlon = conf.nlon;
	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;

	double sum = 0.0;
	for (int i = 0; i < nlat; ++i) {
		for (int j = 0; j < nlon; ++j) {
			sum += cosi[i] * u[i * nlon + j] * v[i * nlon + j];
		}
	}
	return sum * dlat * dlon;
}

double SphereBarvortex::norm(const double * u)
{
	return sqrt(scalar(u, u));
}

double SphereBarvortex::dist(const double * u, const double * v)
{
	long n = conf.nlat * conf.nlon;
	array_t tmp(n);
	vec_sum1(&tmp[0], u, v, 1.0, -1.0, n);
	return norm(&tmp[0]);
}

void SphereBarvortex::S_step (double * out, const double * u, double t)
{
	long nlat    = conf.nlat;
	long nlon    = conf.nlon;
	long n       = conf.nlat * conf.nlon;
	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;
	double tau   = conf.tau;
	double theta = conf.theta;
	double mu    = conf.mu;
	double sigma = conf.sigma;
	double k1    = conf.k1;
	double k2    = conf.k2;
	double nr, nr0 = norm(u);

	array_t w (n);  // w = L(u)
	array_t dw (n); // dw = L(w) = LL(u)

	// next
	array_t u_n (n);
	array_t w_n (n);
	array_t u_n1 (n);

	// jac
	array_t jac0(n);

	// tmp
	array_t tmp1(n);
	array_t tmp2(n);

	//
	array_t FC (n);
	array_t F (n);

	// генерируем правую часть
	// w/dt + mu (1-theta) L w - \sigma(1-theta) w -
	// - k1 J(0.5(u+u), 0.5(w+w)) - k2 J(0.5(u+u), l + h) + f(x, y)

	// w = L (u)
	lapl.calc (&w[0], &u[0]);
	// dw = L (w)
	lapl.calc (&dw[0], &w[0]);

	// w/dt + mu (1-theta) L w
	vec_sum1 (&FC[0], &w[0], &dw[0], 1.0 / tau,
	          mu * (1.0 - theta), n);
	// w/dt + mu (1-theta) L w - \sigma (1-theta) w
	vec_sum1 (&FC[0], &FC[0], &w[0], 1.0,
	          -sigma * (1.0 - theta), n);

	for (int i = 0; i < nlat; ++i)
	{
		double phi    = -0.5 * M_PI + i * dlat;
		for (int j = 0; j < nlon; ++j)
		{
			double lambda = j * dlon;
			if (conf.rp)
			{
				FC[i * nlon + j] += conf.rp (phi, lambda, t, &conf);
			} else if (conf.rp2) {
				FC[i * nlon + j] += conf.rp2[i * nlon + j];
			}
		}
	}

	memcpy(&u_n[0], &u[0], n * sizeof(double));
	memcpy(&w_n[0], &w[0], n * sizeof(double));

	// в FC содержится правая часть, которая не меняется при итерациях!
	
	for (int it = 0; it < 1000; ++it) {
		// k1 J(0.5(u+u), 0.5(w+w)) + k2 J(0.5(u+u), l + h)   =
		// = J(0.5 (u+u), 0.5 k1 (w+w)) + J(0.5 (u+u), k2 (l + h)) =
		// = J(0.5 (u+u), 0.5 k1 (w+w) + k2 (l + h))
		vec_sum1(&tmp1[0], &w_n[0], &w[0], k1 * theta,
				k1 * (1.0 - theta), n);
		vec_sum1(&tmp1[0], &tmp1[0], &lh[0], 1.0, k2, n);
		// 0.5(u+u)
		vec_sum1(&tmp2[0], &u_n[0], &u[0], theta,
				1.0 - theta, n);
		// - k1 J(0.5(u+u), 0.5(w+w)) - k2 J(0.5(u+u), l + h)
		jac.calc(&jac0[0], &tmp2[0], &tmp1[0]);
		vec_sum1(&F[0], &FC[0], &jac0[0], 1.0, -1.0, n);
		lapl.solve(&w_n[0], &F[0], -theta * mu, 1.0 / tau + theta * sigma);
		lapl.solve(&u_n1[0], &w_n[0]);

		nr = dist(&u_n1[0], &u_n[0]);
		u_n1.swap(u_n);
		if (nr / nr0 < 1e-14 || isnan(nr)) {
			break;
		}
	}

	memcpy(out, &u_n[0], n * sizeof(double));
}

void SphereBarvortex::L_step(double *u1, const double *u, const double * z)
{
	long nlat    = conf.nlat;
	long nlon    = conf.nlon;
	long n       = conf.nlat * conf.nlon;
	double theta = conf.theta;
	double mu    = conf.mu;
	double sigma = conf.sigma;
	double tau   = conf.tau;

	array_t z_lapl(n);

	array_t pt1 (n); //лаплас, умноженный на коэф
	array_t pt2 (n); //лаплас в квадрате, умноженный на коэф
	array_t pt3 (n); //якобиан, умноженный на коэф

	lapl.calc(&z_lapl[0], z);
	lapl.calc(&pt1[0], u); //первая часть - лаплас, умноженный на коэф, 

	//умножаем позже, так как лаплас пока нужен
	//		memset(pt3, 0, nn * sizeof(double));
	//funct1.calc(&pt3[0], &u[0], &cor[0], z, &pt1[0], &z_lapl[0]);
	//vector_mult_scalar(&pt3[0], &pt3[0], -conf->rho, nn);

	//		memset(pt2, 0, nn * sizeof(double));
	lapl.calc(&pt2[0], &pt1[0]);
	vec_mult_scalar(&pt2[0], &pt2[0], (1. - theta) * mu, n);

	vec_mult_scalar(&pt1[0], &pt1[0], 1.0 / tau - (1. - theta) * sigma, n);

	memset(u1, 0, n * sizeof(double));
	vec_sum(u1, u1, &pt1[0], n);
	vec_sum(u1, u1, &pt2[0], n);
	vec_sum(u1, u1, &pt3[0], n);

	lapl.solve(u1, u1, - theta * mu, 1.0 / tau + theta * sigma);
	lapl.solve(u1, u1);
}

void SphereBarvortex::LT_step(double *v1, const double *v, const double * z)
{
	long nlat    = conf.nlat;
	long nlon    = conf.nlon;
	long n       = conf.nlat * conf.nlon;
	double theta = conf.theta;
	double mu    = conf.mu;
	double sigma = conf.sigma;
	double tau   = conf.tau;

	array_t z_lapl(n);

	array_t pt1 (n); //лаплас, умноженный на коэф
	array_t pt2 (n); //лаплас в квадрате, умноженный на коэф
	array_t pt3 (n); //якобиан, умноженный на коэф

	lapl.calc(&z_lapl[0], z);

	lapl.solve(v1, v);
	lapl.solve(v1, v1, - theta * mu, 1.0 / tau + theta * sigma);

	lapl.calc(&pt1[0], v1);
	vec_mult_scalar(&pt1[0], &pt1[0], 1.0 / tau - (1. - theta) * sigma, n);

	lapl.calc(&pt2[0], v1);
	lapl.calc(&pt2[0], &pt2[0]);

	vec_mult_scalar(&pt2[0], &pt2[0], (1. - theta) * mu, n);

	// TODO:
	//funct2.calc(&pt3[0], v1, &cor[0], z, 0, &z_lapl[0]);
	//vector_mult_scalar(&pt3[0], &pt3[0], -conf->rho, nn);

	memset(v1, 0, n * sizeof(double));
	vec_sum(v1, v1, &pt1[0], n);
	vec_sum(v1, v1, &pt2[0], n);
	vec_sum(v1, v1, &pt3[0], n);
}

void SphereBarvortex::L_1_step(double *u1, const double *u, const double * z)
{
	long nlat    = conf.nlat;
	long nlon    = conf.nlon;
	long n       = conf.nlat * conf.nlon;
	double theta = conf.theta;
	double mu    = conf.mu;
	double sigma = conf.sigma;
	double tau   = conf.tau;

	array_t z_lapl(n);

	array_t pt1 (n); //лаплас, умноженный на коэф
	array_t pt2 (n); //лаплас в квадрате, умноженный на коэф
	array_t pt3 (n); //якобиан, умноженный на коэф

	lapl.calc(&z_lapl[0], z);
	lapl.calc(&pt1[0], u); //первая часть - лаплас, умноженный на коэф, 

	//умножаем позже, так как лаплас пока нужен
	// TODO: implement
	//jac.calc(&pt3[0], u, &cor[0], z, &pt1[0], &z_lapl[0]);
	//vec_mult_scalar(&pt3[0], &pt3[0], conf->rho, n);

	lapl.calc(&pt2[0], &pt1[0]);
	vec_mult_scalar(&pt2[0], &pt2[0], ( - theta) * mu, n);

	vec_mult_scalar(&pt1[0], &pt1[0], 1 / tau + theta * sigma, n);

	memset(u1, 0, n * sizeof(double));
	vec_sum(u1, u1, &pt1[0], n);
	vec_sum(u1, u1, &pt2[0], n);
	vec_sum(u1, u1, &pt3[0], n);

	lapl.solve(u1, u1, (1.0 - theta) * mu, 1 / tau - (1.0 - theta) * sigma);
	lapl.solve(u1, u1);
}
