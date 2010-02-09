#include "lapl.h"

SphereLaplace::SphereLaplace (long nlat, long nlon) :
		nlat (nlat), nlon (nlon),
		lshaec (nlat* (nlat + 1) + 3* ( (nlat - 2) * (nlat - 1) + nlon + 15) ),
		lshsec (lshaec), lwork(nlat*(2*nlon+3*(nlat+1)+2*nlat+1)), ldwork(nlat + 1),
		wshaec(lshaec), wshsec(lshsec), work(lwork), dwork(ldwork)
{
	long ierror = 0;
	shaeci_ (&nlat, &nlon, &wshaec[0], &lshaec, &dwork[0], &ldwork, &ierror);
	shseci_ (&nlat, &nlon, &wshsec[0], &lshsec, &dwork[0], &ldwork, &ierror);
}

SphereLaplace::~SphereLaplace() {}

void SphereLaplace::solve (double * out, const double * in, double diag)
{
	long ierror = 0;
	long nt   = 1;
	long isym = 0;
	double pertrb = 0;

	array_t a (nlat * nlat);
	array_t b (nlat * nlat);

	shaec_ (&nlat, &nlon, &isym, &nt, (double*)in,
	        &nlat, &nlon, &a[0], &b[0],
	        &nlat, &nlat,
	        &wshaec[0], &lshaec,
	        &work[0], &lwork, &ierror);

	islapec_ (&nlat, &nlon, &isym, &nt, &diag,
	          out, &nlat, &nlon,
	          &a[0], &b[0], &nlat, &nlat,
	          &wshsec[0], &lshsec, &work[0], &lwork, &pertrb, &ierror);
}

