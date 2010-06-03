#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <list>

#include "linal.h"
#include "spherepack.h"
#include "operator.h"

using namespace std;
using namespace linal;

struct SphereOperator::Data
{
	int links;
	list < array_t > arrays;

	Data() : links(1)
	{
	}

	array_t & create(size_t size)
	{
		arrays.push_back(array_t(size));
		return arrays.back();
	}
};

SphereOperator::SphereOperator(const SphereOperator & op): SphereNorm(op.nlat, op.nlon),
		d(op.d),
		isym(op.isym), nlat (op.nlat), nlon (op.nlon), 
		mdab ( (nlon + 2) / 2), mdb((nlon+1)/2), mdc( (nlon+1)/2 ),
		slsave (op.slsave), swsave (op.swsave),
		sldwork (op.sldwork), sdwork (op.sdwork),
		vlsave (op.vlsave), vwsave (op.vwsave),
		vldwork (op.vldwork), vdwork (op.vdwork),

		islsave (op.islsave), iswsave (op.iswsave),
		isldwork (op.isldwork), isdwork (op.isdwork),
		ivlsave (op.ivlsave), ivwsave (op.ivwsave),
		ivldwork (op.ivldwork), ivdwork (op.ivdwork),

		lwork (op.lwork), work (op.work)
{
	++d->links;
}

SphereOperator::SphereOperator (long nlat, long nlon, long isym) : SphereNorm(nlat, nlon),
		d(new Data),
		isym(isym), nlat (nlat), nlon (nlon), 
		mdab ( (nlon + 2) / 2), mdb((nlon+1)/2), mdc( (nlon+1)/2 ),
		slsave (5*nlat*nlat*nlon), swsave (d->create(slsave)),
		sldwork (4*nlat*nlat), sdwork (d->create(sldwork)),
		vlsave (5*nlat*nlat*nlon), vwsave (d->create(slsave)),
		vldwork (4*nlat*nlat), vdwork (d->create(vldwork)),

		islsave (5*nlat*nlat*nlon), iswsave (d->create(islsave)),
		isldwork (4*nlat*nlat), isdwork (d->create(isldwork)),
		ivlsave (5*nlat*nlat*nlon), ivwsave (d->create(islsave)),
		ivldwork (4*nlat*nlat), ivdwork (d->create(ivldwork)),

		lwork (5*nlat*nlat*nlon), work (d->create(lwork))
{
	long ierror = 0;
	shaeci_ (&nlat, &nlon, &swsave[0], &slsave, &sdwork[0], &sldwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "shaeci_ error\n");
		exit(1);
	}

	vhseci_ (&nlat, &nlon, &vwsave[0], &vlsave, &vdwork[0], &vldwork, &ierror);
	if (ierror != 0) {
		fprintf(stderr, "vhseci_ error\n");
		exit(1);
	}

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

SphereOperator::~SphereOperator()
{
	--d->links;
	if (!d->links) 
	{
		delete d;
	}
}

