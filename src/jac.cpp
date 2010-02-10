#include "spherepack.h"
#include "jac.h"
#include "utils.h"

using namespace std;

SphereJacobian::SphereJacobian (long nlat, long nlon) : nlat (nlat), nlon (nlon)
{
	long ierror = 0;
	shaeci_ (&nlat, &nlon, &swsave[0], &slsave, &sdwork[0], &sldwork, &ierror);
	vhseci_ (&nlat, &nlon, &vwsave[0], &vlsave, &vdwork[0], &vldwork, &ierror);
}

SphereJacobian::~SphereJacobian()
{
}

void SphereJacobian::calc (double * out, const double * u1, const double * v1)
{
	// TODO: use gradec here
	long ierror = 0;
	long isym   = 0;
	long nt     = 0;
	long n      = nlat * nlon;

	array_t a(nlat * nlat);
	array_t b(nlat * nlat);

	array_t dutheta(nlat * nlon);
	array_t dvtheta(nlat * nlon);

	array_t duphi(nlat * nlon);
	array_t dvphi(nlat * nlon);

	array_t u(nlat * nlon);
	array_t v(nlat * nlon);
	array_t j(nlat * nlon);

	transpose(&u[0], &u1[0], nlat, nlon);
	transpose(&v[0], &v1[0], nlat, nlon);

	shaec_ (&nlat, &nlon, &isym, &nt, &u[0], &nlat, &nlon, &a[0], &b[0], &mdab, &nlat, &swsave[0],
	       &slsave, &work[0], &lwork, &ierror);

	gradec_ (&nlat, &nlon, &isym, &nt, &dutheta[0], &duphi[0], &nlat, &nlon,
	        &a[0], &b[0], &mdab, &nlat, &vwsave[0],
	        &vlsave, &work[0], &lwork, &ierror);

	shaec_ (&nlat, &nlon, &isym, &nt, &v[0], &nlat, &nlon, &a[0], &b[0], &mdab, &nlat, &swsave[0],
	       &slsave, &work[0], &lwork, &ierror);

	gradec_ (&nlat, &nlon, &isym, &nt, &dvtheta[0], &dvphi[0], &nlat, &nlon,
	        &a[0], &b[0], &mdab, &nlat, &vwsave[0],
	        &vlsave, &work[0], &lwork, &ierror);

	for (int i = 0; i < n; ++i) {
		j[i] = dutheta[i] * dvphi[i] - dvtheta[i] * duphi[i];
	}

	transpose(out, &j[0], nlon, nlat);
}

