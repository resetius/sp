#ifndef SPHERE_OP_H
#define SPHERE_OP_H

#include <vector>

class SphereOperator
{
protected:
	typedef std::vector < double > array_t;

	long nlat;
	long nlon;
	long mdab;
	long mdb;

	long slsave;
	array_t swsave;

	long sldwork;
	array_t sdwork;

	long vlsave;
	array_t vwsave;

	long vldwork;
	array_t vdwork;

	long islsave;
	array_t iswsave;

	long isldwork;
	array_t isdwork;

	long ivlsave;
	array_t ivwsave;

	long ivldwork;
	array_t ivdwork;

	long lwork;
	array_t work;

public:
	SphereOperator(long nlat, long lon);
	~SphereOperator();
};

#endif /* SPHERE_OP_H */

