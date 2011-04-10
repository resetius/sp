#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <vector>

#include "vorticity.h"
#include "lapl.h"
#include "linal.h"
#include "grad.h"

using namespace std;
using namespace linal;

typedef vector < double > array_t;

void random_array(array_t & u)
{
	for (size_t i = 0; i < u.size(); ++i) {
		u[i] = (double)rand() / (double)RAND_MAX;
	}
}

static bool run_test_grad(
	SphereLaplace & lapl,
	SphereVorticity & vor,
	SphereGrad & grad,
	long nlat, long nlon)
{
	long n = nlat * nlon;
	double nr;

	array_t u(n);
	array_t v(n);

	array_t u1(n);
	array_t v1(n);

	array_t omg(n);
	array_t psi(n);
	array_t psi1(n);

	random_array(u);
	random_array(v);

	vor.calc(&omg[0], &u[0], &v[0]);
	lapl.solve(&psi[0], &omg[0]);

	grad.calc(&u1[0], &v1[0], &psi[0]);
	grad.solve(&psi1[0], &u1[0], &v1[0]);

	nr = grad.dist(&psi1[0], &psi[0]);
	fprintf(stderr, "nr1 = %.16le\n", nr);

	return nr < 1e-12;
}

static bool run_test_vorgrad(
	SphereLaplace & lapl,
	SphereVorticity & vor,
	SphereGrad & grad,
	long nlat, long nlon)
{
	long n = nlat * nlon;
	int i, j;

	array_t u(n);
	array_t v(n);

	array_t u1(n);
	array_t v1(n);

	array_t omg(n);
	array_t psi(n);
	array_t psi1(n);

	double nr;

	//random_array(u);
	//random_array(v);

	double dlat = M_PI / (nlat - 1);
	double dlon = 2. * M_PI / nlon;

	for (i = 0; i < nlat; ++i)
	{
		double phi    = -0.5 * M_PI + i * dlat;

		for (j = 0; j < nlon; ++j)
		{
			double lambda = j * dlon;
			u[i * nlon + j]   =  (M_PI / 2 + phi) * (M_PI / 2. - phi);
		}
	}

	vor.calc(&omg[0], &u[0], &v[0]);
	lapl.solve(&psi[0], &omg[0]);

	grad.calc(&v1[0], &u1[0], &psi[0]);
	vec_mult_scalar(&u1[0], &u1[0], -1.0, n);

	nr = grad.dist(&u1[0], &u[0]);
	fprintf(stderr, "nr1 = %.16le\n", nr);

	return nr < 1e-12;
}

extern "C" int test_op(int argc, char ** argv) 
{
	long nlat = 19;
	long nlon = 36;
	SphereOperator op(nlat, nlon, 0);
	SphereVorticity vor(op);
	SphereLaplace lapl(op);
	SphereGrad grad(op);
	srand(3214562);

	bool test_vor     = false;
	bool test_vorgrad = false;
	bool test_grad    = false;
	bool res = true;

	if (argc == 0) {
		test_vor     = true;
		test_vorgrad = true;
		test_grad    = true;
	} else {
		for (int i = 0; i < argc; ++i) {
			if (!strcmp(argv[i], "vorticity")) {
				test_vor     = true;
			} else if (!strcmp(argv[i], "grad")) {
				test_grad    = true;
			} else if (!strcmp(argv[i], "vorgrad")) {
				test_vorgrad = true;
			}
		}
	}

	if (test_vor) {
		res &= vor.test();
	}
	if (test_vorgrad) {
		res &= run_test_vorgrad(lapl, vor, grad, nlat, nlon);
	}
	if (test_grad) {
		res &= run_test_grad(lapl, vor, grad, nlat, nlon);
	}

	if (res) {
		return 0;
	} else {
		return -1;
	}
}
