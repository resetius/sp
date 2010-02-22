#include "vorticity.h"
#include "spherepack.h"
#include "utils.h"

using namespace std;

SphereVorticity::SphereVorticity(long nlat, long nlon):
	nlat(nlat), nlon(nlon),  mdc( (nlon+1)/2 ),

	islsave (5*nlat*nlat*nlon), iswsave (islsave),
	isldwork (4*nlat*nlat), isdwork (isldwork),
	ivlsave (5*nlat*nlat*nlon), ivwsave (islsave),
	ivldwork (4*nlat*nlat), ivdwork (ivldwork),

	lwork (5*nlat*nlat*nlon), work (lwork)
{
	long ierror = 0;

	shseci_ (&nlat, &nlon, &iswsave[0], &islsave, &isdwork[0], &isldwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "shseci_ error\n");
		exit(1);
	}

	vhaeci_ (&nlat, &nlon, &ivwsave[0], &ivlsave, &ivdwork[0], &ivldwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "vhseci_ error\n");
		exit(1);
	}
}

SphereVorticity::~SphereVorticity()
{
}

void SphereVorticity::calc(double * div, const double * u, const double *v)
{
	long ierror = 0;
	long isym   = 0;
	long ityp   = 0;
	long nt     = 1;
	long n      = nlat * nlon;

	array_t br(mdc * nlat * nt), bi(mdc * nlat * nt);
	array_t cr(mdc * nlat * nt), ci(mdc * nlat * nt);

	array_t vv(n);
	array_t ww(n);
	array_t vt(n);

	transpose (&vv[0], &v[0], nlat, nlon);
	transpose (&ww[0], &u[0], nlat, nlon);

	vhaec_(&nlat, &nlon, &ityp, &nt, &vv[0], &ww[0], &nlat, &nlon,
		&br[0], &bi[0], &cr[0], &ci[0], &mdc, &nlat, 
		&ivwsave[0], &ivlsave, &work[0], &lwork, &ierror);

	if (ierror != 0) {
		fprintf(stderr, "vhaec_ error\n");
		exit(1);
	}

	vrtec_(&nlat, &nlon, &isym, &nt, &vt[0], &nlat, &nlon,
		&cr[0], &ci[0], &mdc, &nlat,
		&iswsave[0], &islsave, &work[0], &lwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "vrtec_ error\n");
		exit(1);
	}

	transpose (&div[0], &vt[0], nlon, nlat);
}
