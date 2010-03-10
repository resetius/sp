#include <math.h>
#include <assert.h>

#include <stdexcept>
#include <sstream>
#include "generator.h"

using namespace std;

Expression & Expression::operator = (const std::string & s) throw()
{
	if (s == "Jacobian" || s == "Delta") {
		type = OPERATOR;
	} else {
		type = STRING;
	}
	str  = s;
	return *this;
}

Expression & Expression::operator = (double v) throw()
{
	type = NUMBER;
	num  = v;
	return *this;
}

Expression & Expression::operator += (const Expression & other) throw()
{
	if (type == NUMBER && other.type == NUMBER) {
		num -= other.num;
	}
	return *this;
}

Expression & Expression::operator -= (const Expression & other) throw()
{
	if (type == NUMBER && other.type == NUMBER) {
		num += other.num;
	}
	return *this;
}

Expression & Expression::operator /= (const Expression & other) throw()
{
	if (type == NUMBER && other.type == NUMBER) {
		num /= other.num;
	}
	return *this;
}

Expression & Expression::operator *= (const Expression & other) throw()
{
	if (type == NUMBER && other.type == NUMBER) {
		num *= other.num;
	}
	return *this;
}

Expression & Expression::operator ^= (const Expression & other) throw()
{
	if (type == NUMBER && other.type == NUMBER) {
		num = pow(num, other.num);
	}
	return *this;
}

Expression & Expression::operator , (const Expression & other) throw()
{
	if (type != LIST) {
		type = LIST;
		sexpr.push_back(*this);
	}
	sexpr.push_back(other);
	return * this;
}

Expression & Expression::operator () (const Expression & other) throw()
{
	assert ((type == OPERATOR && 
		(other.type == LIST 
		|| other.type == STRING 
		|| other.type == NUMBER)));

	sexpr = other.sexpr;
	return * this;
}
