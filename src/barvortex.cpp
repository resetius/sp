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

#include <math.h>
#include <string.h>

#include "barvortex.h"
#include "linal.h"

using namespace std;
using namespace linal;

SphereBarvortex::SphereBarvortex (const SphereBarvortex::Conf & conf) : SphereNorm(conf.nlat, conf.nlon),
		conf (conf), op (conf.nlat, conf.nlon, 0), lapl(op), jac (op),
		lh (conf.nlat * conf.nlon)
{
	long nlat = conf.nlat;
	long nlon = conf.nlon;
	double dlat = M_PI / (nlat-1);
	double dlon = 2. * M_PI /nlon;


	for (int i = 0; i < nlat; ++i)
	{
		double phi    = -0.5 * M_PI + i * dlat;
		for (int j = 0; j < nlon; ++j)
		{
			double lambda = j * dlon;
			if (conf.cor) {
				lh[i * nlon + j] = conf.cor(phi, lambda, 0, &conf);
			} else if (conf.cor2) {
				lh[i * nlon + j] = conf.cor2[i * nlon + j];
			}
		}
	}
}

SphereBarvortex::~SphereBarvortex()
{
}

void SphereBarvortex::S_step (double * out, const double * u, double t)
{
       if(1==0)
            S_step_AVO_KN (out,u,t);
       else 
            S_step_BVL_KN (out,u,t);
}





void SphereBarvortex::S_step_AVO_KN (double * out, const double * u, double t)
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
			if (conf.rp) {
				FC[i * nlon + j] += conf.rp (phi, lambda, t, &conf);
			} 
			if (conf.rp2) {
				FC[i * nlon + j] += conf.rp2[i * nlon + j];
			}
			if (conf.rp3) {
				FC[i * nlon + j] += conf.rp3[i * nlon + j];
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

		lapl.filter(&w_n[0], &w_n[0]);
		lapl.filter(&u_n1[0], &u_n1[0]);

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
	double k1    = conf.k1;
	double k2    = conf.k2;
	double nr, nr0 = norm(u);

	if (fabs(nr0) < 1e-10) nr0 = 1.0; // use absolute norm

	array_t tmp1(n);
	array_t tmp2(n);

	array_t dz(n);
	array_t dw(n);
	array_t w(n);

	// next
	array_t u_n (n);
	array_t w_n (n);
	array_t u_n1 (n);

	// jac
	array_t jac0(n);
	array_t jac1(n);
	array_t jac2(n);

	// rp
	array_t FC(n);
	array_t F(n);

	// dz = L(z)
	lapl.calc(&dz[0], z);
	// w = L(u)
	lapl.calc(&w[0], u);
	// dw = L(w)
	lapl.calc(&dw[0], &w[0]);

	// w/dt + mu (1-theta) L w
	vec_sum1 (&FC[0], &w[0], &dw[0], 1.0 / tau,
	          mu * (1.0 - theta), n);
	// w/dt + mu (1-theta) L w - \sigma (1-theta) w
	vec_sum1 (&FC[0], &FC[0], &w[0], 1.0,
	          -sigma * (1.0 - theta), n);

	memcpy(&u_n[0], &u[0], n * sizeof(double));
	memcpy(&w_n[0], &w[0], n * sizeof(double));

	// в FC содержится правая часть, которая не меняется при итерациях!

	for (int it = 0; it < 1000; ++it) {
		// k1 J(0.5(u+u), dz) + k1 J(z, 0.5(w+w)) + k2 J(0.5(u+u), l + h)

		vec_sum1(&tmp1[0], &w_n[0], &w[0], theta,
				(1.0 - theta), n);
		vec_sum1(&tmp2[0], &u_n[0], &u[0], theta,
				(1.0 - theta), n);

		// J(0.5(u+u), dz)
		jac.calc(&jac0[0], &tmp2[0],  &dz[0]);
		// J(z, 0.5(w+w))
		jac.calc(&jac1[0], &z[0],   &tmp1[0]);
		// J(0.5(u+u), l + h)
		jac.calc(&jac2[0], &tmp2[0], &lh[0]);

		vec_sum1(&jac0[0], &jac0[0], &jac1[0], k1,  k1, n);
		vec_sum1(&jac0[0], &jac0[0], &jac2[0], 1.0, k2, n);

		vec_sum1(&F[0], &FC[0], &jac0[0], 1.0, -1.0, n);
	
		lapl.solve(&w_n[0], &F[0], -theta * mu, 1.0 / tau + theta * sigma);
		lapl.solve(&u_n1[0], &w_n[0]);

		lapl.filter(&w_n[0], &w_n[0]);
		lapl.filter(&u_n1[0], &u_n1[0]);

		nr = dist(&u_n1[0], &u_n[0]);
		u_n1.swap(u_n);
		if (nr / nr0 < 1e-14 || isnan(nr)) {
			break;
		}
	}

	memcpy(u1, &u_n[0], n * sizeof(double));
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

	array_t tmp(n);
	array_t p_lapl(n);
	jac.calc_t(&tmp[0], v1, &lh[0]);
	jac.calc_t(&p_lapl[0], z, v1);
	lapl.calc(&p_lapl[0], &p_lapl[0]);
	jac.calc_t(&pt3[0], v1, &z_lapl[0]);
	vec_sum(&pt3[0], &pt3[0], &p_lapl[0], n);
	vec_sum(&pt3[0], &pt3[0], &tmp[0], n);
	
	//vector_mult_scalar(&pt3[0], &pt3[0], -conf->rho, nn); //TODO: k1, k2

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
	double k1    = conf.k1;
	double k2    = conf.k2;
	double nr, nr0 = norm(u);

	array_t tmp1(n);
	array_t tmp2(n);

	array_t dz(n);
	array_t dw(n);
	array_t w(n);

	// next
	array_t u_n (n);
	array_t w_n (n);
	array_t u_n1 (n);

	// jac
	array_t jac0(n);
	array_t jac1(n);
	array_t jac2(n);

	// rp
	array_t FC(n);
	array_t F(n);

	// dz = L(z)
	lapl.calc(&dz[0], z);
	// w = L(u)
	lapl.calc(&w[0], u);
	// dw = L(w)
	lapl.calc(&dw[0], &w[0]);

	// w/dt - mu theta L w
	vec_sum1 (&FC[0], &w[0], &dw[0], 1.0 / tau,
	          -mu * theta, n);
	// w/dt - mu theta L w + \sigma theta w
	vec_sum1 (&FC[0], &FC[0], &w[0], 1.0,
	          sigma * theta, n);

	memcpy(&u_n[0], &u[0], n * sizeof(double));
	memcpy(&w_n[0], &w[0], n * sizeof(double));

	// в FC содержится правая часть, которая не меняется при итерациях!

	for (int it = 0; it < 1000; ++it) {
		// k1 J(0.5(u+u), dz) + k1 J(z, 0.5(w+w)) + k2 J(0.5(u+u), l + h)

		vec_sum1(&tmp1[0], &w_n[0], &w[0], 1-theta,
				theta, n);
		vec_sum1(&tmp2[0], &u_n[0], &u[0], 1-theta,
				theta, n);

		// J(0.5(u+u), dz)
		jac.calc(&jac0[0], &tmp2[0],  &dz[0]);
		// J(z, 0.5(w+w))
		jac.calc(&jac1[0], &z[0],   &tmp1[0]);
		// J(0.5(u+u), l + h)
		jac.calc(&jac2[0], &tmp2[0], &lh[0]);

		vec_sum1(&jac0[0], &jac0[0], &jac1[0], k1,  k1, n);
		vec_sum1(&jac0[0], &jac0[0], &jac2[0], 1.0, k2, n);

		vec_sum1(&F[0], &FC[0], &jac0[0], 1.0, 1.0, n);
	
		lapl.solve(&w_n[0], &F[0], (1-theta) * mu, 1.0 / tau - (1-theta) * sigma);
		lapl.solve(&u_n1[0], &w_n[0]);

		lapl.filter(&w_n[0], &w_n[0]);
		lapl.filter(&u_n1[0], &u_n1[0]);

		nr = dist(&u_n1[0], &u_n[0]);
		u_n1.swap(u_n);
		if (nr / nr0 < 1e-14 || isnan(nr)) {
			break;
		}
	}

	memcpy(u1, &u_n[0], n * sizeof(double));
}

void SphereBarvortex::p2u(double * u, const double * p)
{
	memcpy(u, p, conf.nlat * conf.nlon * sizeof(double));
}

void SphereBarvortex::u2p(double * p, const double * u)
{
	memcpy(p, u, conf.nlat * conf.nlon * sizeof(double));
}



//------------------------------------------BVL
//
//

void SphereBarvortex::S_step_BVL_KN (double * psn, const double * psi, double t)
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
	double nr,nr1,nr0=norm(psi);
	int it;


	double mu2   = conf.mu*conf.mu/50.; // Evristic formula

	const int CHECK = 1;   // 0,1,2,3
	const int PRESTEP = 1; // 0,1
	const int MAX_IT = 1000; // 0,1
	const int PRN_IT = 300; // 0,1
	
	
	array_t omg (n);  // omg = Dleta(psi)
	array_t domg (n); // domg = Delta(omg)
	array_t Domg (n); // domg = Delta^3(omg)

	// next
	array_t omn (n);
	array_t domn (n);
	array_t Domn (n);

	// mean 
	array_t om2 (n); // om2 = th*omn+(1-th)omg
	array_t ps2 (n); // ps2 = Delta(om2)
	
	// jac
	array_t jac2(n);

	// tmp
	array_t tmp1(n);
	array_t tmp2(n);

	//
	array_t FC (n);
	array_t F (n);

	// Equation are
	// Delta psi = omg;
	// (omn - omg)/tau =       mu D (th psn + (1-th)psi) + 
	//			mu2 D^3 (th psn + (1-th)psi) 
	// 			- \sigma D (th psn - (1-th)psi)
	// - k1 J( th*psn + (1-th)*psi), th*omn+(1-th)omg)  
	// - k2 J( th*psn + (1-th)*psi), l + h) 
	// + f(x, y)


	lapl.calc (&omg[0], &psi[0]);

	lapl.calc (&domg[0], &omg[0]);
	lapl.calc (&tmp1[0], &domg[0]);
	lapl.calc (&Domg[0], &tmp1[0]);

/*	if(CHECK){
	lapl.solve_3(&tmp1[0], &Domg[0], 1.0, 0. );
	vec_sum1 (&tmp1[0], &tmp1[0], &omg[0], 1.0 , - 1.0, n);
	lapl.filter(&tmp1[0], &tmp1[0]);
	nr = norm (&tmp1[0]);
	fprintf(stderr, "err_solve_3=%g \n", nr);
	if(CHECK==2) exit(1);
	}
*/	
	if(CHECK){
	double mu2=10.0, mu = 20.,  sigma = 10.;
	vec_sum1 (&tmp1[0], &Domg[0], &domg[0], mu2 , mu, n);
	vec_sum1 (&tmp1[0], &tmp1[0], &omg[0], 1.0 , -sigma, n);

	lapl.solve_l3(&F[0], &tmp1[0], mu2, mu, sigma);
	vec_sum1 (&tmp1[0], &F[0], &omg[0], 1.0 , - 1.0, n);
	lapl.filter(&tmp1[0], &tmp1[0]);
	nr = norm (&tmp1[0]);
	if(CHECK>1)fprintf(stderr, "%s: err_solve_l3=%g \n", __FUNCTION__, nr);
	if(CHECK>2) exit(1);
	}
	
	
	
	
	
	
	
	vec_sum1 (&FC[0], &omg[0], &Domg[0], 1.0 / tau,
	          mu2 * (1.0 - theta), n);
	vec_sum1 (&FC[0], &FC[0], &domg[0], 1.0,
	          mu * (1.0 - theta), n);
	vec_sum1 (&FC[0], &FC[0], &omg[0], 1.0,
	          -sigma * (1.0 - theta), n);


	for (int i = 0; i < nlat; ++i)
	{
		double phi    = -0.5 * M_PI + i * dlat;
		for (int j = 0; j < nlon; ++j)
		{
			double lambda = j * dlon;
			if (conf.rp) {
				F[i * nlon + j] = conf.rp (phi, lambda, t, &conf);
			} 
			if (conf.rp2) {
				F[i * nlon + j] = conf.rp2[i * nlon + j];
			}
			if (conf.rp3) {
				F[i * nlon + j] = conf.rp3[i * nlon + j];
			}
		}
	}


	vec_sum1 (&FC[0], &FC[0], &F[0], 1.0, 1.0, n);




	if(PRESTEP){

		vec_sum1(&tmp1[0], &omg[0], &lh[0], k1, k2, n);
		jac.calc(&jac2[0], &psi[0], &tmp1[0]);
		vec_sum1(&tmp1[0], &domg[0], &omg[0], mu, -sigma, n);
		vec_sum1(&tmp1[0], &tmp1[0], &Domg[0], 1.0, mu2, n);
		vec_sum1(&tmp1[0], &tmp1[0], &jac2[0], 1.0, -1.0, n);
		vec_sum1(&tmp1[0], &tmp1[0], &F[0], 1.0, 1.0, n);
		vec_sum1(&omn[0] , &omg[0], &tmp1[0], 1.0, tau, n);

		lapl.solve(&psn[0], &omn[0]);
		lapl.calc(&omn[0], &psn[0]);
		lapl.calc (&domn[0], &omn[0]);
		lapl.calc (&tmp1[0], &domn[0]);
		lapl.calc (&Domn[0], &tmp1[0]);


		}else
		{
		memcpy(&psn[0], &psi[0], n * sizeof(double));
		memcpy(&omn[0], &omg[0], n * sizeof(double));
		memcpy(&domn[0], &domg[0], n * sizeof(double));
		memcpy(&Domn[0], &Domg[0], n * sizeof(double));
		}
	
	// Solve nonlinear equation by the simple iteration method

	for (it = 0; it < MAX_IT; ++it) {
		vec_sum1(&om2[0], &omn[0], &omg[0], theta,
				(1.0 - theta), n);
		vec_sum1(&tmp1[0], &om2[0], &lh[0], k1, k2, n);
		vec_sum1(&ps2[0], &psn[0], &psi[0], theta,
				1.0 - theta, n);
		jac.calc(&jac2[0], &ps2[0], &tmp1[0]);
		vec_sum1(&tmp2[0], &FC[0], &jac2[0], 1.0, -1.0, n);

		{
		vec_sum1(&tmp1[0], &tmp2[0], &omn[0], 1.0, - theta*sigma, n);
		vec_sum1(&tmp1[0], &tmp1[0], &domn[0], 1.0, theta*mu, n);
		vec_sum1(&tmp1[0], &tmp1[0], &Domn[0], 1.0, theta*mu2, n);
		vec_sum1(&tmp1[0], &tmp1[0], &omn[0], 1.0, -1.0/tau , n);

		lapl.filter(&tmp1[0], &tmp1[0]);

		nr = norm(&tmp1[0]);
		nr1 = norm(&omg[0]);

		if (nr/nr1 < 1e-10 || isnan(nr) || isnan(nr1)) {
				break;
				}

		if( (it+1)%PRN_IT==0) 
			fprintf(stderr, "%s :nr=%g nr1=%g it=%d \n", __FUNCTION__,nr,nr1,it);
		}

		lapl.solve_l3(&omn[0], &tmp2[0],-theta*mu2, 
				-theta * mu, - 1.0 / tau - theta * sigma);
		

		lapl.solve(&tmp1[0], &omn[0]);

		lapl.filter(&omn[0], &omn[0]);
		lapl.filter(&tmp1[0], &tmp1[0]);

		nr = dist(&tmp1[0], &psn[0]);
		memcpy(psn, &tmp1[0], n * sizeof(double));
		if (nr / nr0 < 1e-14 || isnan(nr)) {
			break;
		}


	lapl.calc (&domn[0], &omn[0]);
	lapl.calc (&tmp1[0], &domn[0]);
	lapl.calc (&Domn[0], &tmp1[0]);

	}

	if(CHECK){


		for (int i = 0; i < nlat; ++i)
		{
			double phi    = -0.5 * M_PI + i * dlat;
			for (int j = 0; j < nlon; ++j)
			{
				double lambda = j * dlon;
				if (conf.rp) {
					F[i * nlon + j] = conf.rp (phi, lambda, t, &conf);
				} 
				if (conf.rp2) {
					F[i * nlon + j] = conf.rp2[i * nlon + j];
				}
				if (conf.rp3) {
					F[i * nlon + j] = conf.rp3[i * nlon + j];
				}
			}
		}



		lapl.calc (&omn[0], &psn[0]);
		vec_sum1(&om2[0], &omn[0], &omg[0], theta,
				(1.0 - theta), n);
		vec_sum1(&tmp1[0], &om2[0], &lh[0], k1, k2, n);
		vec_sum1(&ps2[0], &psn[0], &psi[0], theta,
				1.0 - theta, n);
		jac.calc(&jac2[0], &ps2[0], &tmp1[0]);
		vec_sum1(&tmp2[0], &F[0], &jac2[0], 1.0, -1.0, n);


		lapl.calc (&domn[0], &om2[0]);
		lapl.calc (&tmp1[0], &domn[0]);
		lapl.calc (&Domn[0], &tmp1[0]);

		vec_sum1(&tmp1[0], &tmp2[0], &om2[0], 1.0, - sigma, n);
		vec_sum1(&tmp1[0], &tmp1[0], &domn[0], 1.0, mu, n);
		vec_sum1(&tmp1[0], &tmp1[0], &Domn[0], 1.0, mu2, n);
		vec_sum1(&tmp1[0], &tmp1[0], &omg[0], 1.0, 1.0/tau , n);
		vec_sum1(&tmp1[0], &tmp1[0], &omn[0], 1.0, -1.0/tau , n);

		lapl.filter(&tmp1[0], &tmp1[0]);

		nr = norm(&tmp1[0]);

		if( nr>1.e-5 ) fprintf(stderr, "%s :::nr=%g it=%d \n", __FUNCTION__,nr,it);



		}

		
}		



