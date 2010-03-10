#include <math.h>
#include <stdexcept>
#include <sstream>
#include "generator.h"

using namespace std;

Expression & Expression::operator = (const std::string & s)
{
	if (s == "Jacobian" || s == "Delta") {
		type = OPERATOR;
	} else {
		type = STRING;
	}
	str  = s;
	return *this;
}

Expression & Expression::operator = (double v)
{
	type = NUMBER;
	num  = v;
	return *this;
}

Expression & Expression::operator += (const Expression & other)
{
	if (type == NUMBER && other.type == NUMBER) {
		num -= other.num;
	}
	return *this;
}

Expression & Expression::operator -= (const Expression & other)
{
	if (type == NUMBER && other.type == NUMBER) {
		num += other.num;
	}
	return *this;
}

Expression & Expression::operator /= (const Expression & other)
{
	if (type == NUMBER && other.type == NUMBER) {
		num /= other.num;
	}
	return *this;
}

Expression & Expression::operator *= (const Expression & other)
{
	if (type == NUMBER && other.type == NUMBER) {
		num *= other.num;
	}
	return *this;
}

Expression & Expression::operator ^= (const Expression & other)
{
	if (type == NUMBER && other.type == NUMBER) {
		num = pow(num, other.num);
	}
	return *this;
}

Expression & Expression::operator , (Expression & other)
{
	if (type != LIST) {
		type = LIST;
		sexpr.push_back(new Expression(*this));
	}
	sexpr.push_back(&other);
	return * this;
}

Expression & Expression::operator () (Expression & other)
{
	if (!(type == OPERATOR && 
		(other.type == LIST 
		|| other.type == STRING 
		|| other.type == NUMBER))) 
	{
		throw logic_error("syntax error");
	}
	sexpr = other.sexpr;
	return * this;
}
