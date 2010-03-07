#include <stdexcept>
#include <sstream>
#include "generator.h"

using namespace std;

Parser::Parser(): 
	error(false), 
	precission(0),
	method(0)
{
}

Parser::~Parser()
{
}

void Parser::set_name(const std::string & n)
{
	name = n;
}

void Parser::set_precission(const std::string & p)
{
	if (p == "float") {
		precission = FLOAT;
	} else if (p == "double") {
		precission = DOUBLE;
	} else if (p == "both") {
		precission = BOTH;
	} else {
		std::stringstream str;
		str << "unknown precission " << p;
		throw std::logic_error(str.str());
	}
}

void Parser::set_method(const std::string & p)
{
	if (p == "harmonic") {
		method = HARMONIC;
	} else {
		std::stringstream str;
		str << "unknown method " << p;
		throw std::logic_error(str.str());
	}
}

void Parser::set_surface(const std::string & p)
{
}

void Parser::set_approximation(const std::string &p)
{
}

void Parser::add_map(const std::string & from, const std::string & to)
{
}

void Parser::add_scalar(const std::string &p)
{
}

void Parser::add_function(const std::string &p, int args)
{
}

void Parser::new_equation()
{
}
