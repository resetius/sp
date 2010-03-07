#include <stdexcept>
#include <sstream>
#include <map>
#include "generator.h"

using namespace std;

class Registrator
{
	Registrator() {}
	~Registrator() {}

public:
	map < string, generator_creator_t > generators;

	static Registrator & instance()
	{
		static Registrator r;
		return r;
	}
};

Generator::Generator(): precission(0)
{
}

Generator::~Generator()
{
}

int register_generator(const std::string & name, generator_creator_t creator)
{
	Registrator::instance().generators.insert(make_pair(name, creator));
	return 1;
}

Generator * make_generator(const std::string & name)
{
	map < string, generator_creator_t > & generators = Registrator::instance().generators;
	if (generators.find(name) != generators.end()) {
		return generators[name]();
	}
	return 0;
}

void Generator::set_precission(const std::string & p)
{
	if (p == "double") {
		precission = DOUBLE;
	} else if (p == "float") {
		precission = FLOAT;
	} else if (p == "both") {
		precission = BOTH;
	} else {
		stringstream str;
		str << "unknown precission " << p;
		throw std::logic_error(str.str());
	}
}

void Generator::add_map(const std::string & from, const std::string & to)
{
}

void Generator::add_scalar(const std::string &p)
{
}

void Generator::add_function(const std::string &p, int args)
{
}

void Generator::new_equation()
{
}

Parser::Parser()
{
}

Parser::~Parser()
{
}

void Parser::set_name(const std::string & n)
{
	name = n;
}

void Parser::check(const string &p)
{
	if (!generator.get())
	{
		stringstream str;
		str << p << ": declare method first ! ";
		throw logic_error(str.str());
	}
}

void Parser::set_precission(const std::string & p)
{
	check("precission");
	generator->set_precission(p);
}

void Parser::set_method(const std::string & p)
{
	if (generator.get()) {
		throw std::logic_error("method already set");
	}

	generator.reset(make_generator(p));

	if (!generator.get()) {
		std::stringstream str;
		str << "unknown method " << p;
		throw std::logic_error(str.str());
	}
}

void Parser::set_surface(const std::string & p)
{
	check("surface");
	generator->set_surface(p);
}

void Parser::set_approximation(const std::string &p)
{
	check("approximation");
	generator->set_approximation(p);
}

void Parser::add_map(const std::string & from, const std::string & to)
{
	check("map");
	generator->add_map(from, to);
}

void Parser::add_scalar(const std::string &p)
{
	check("scalar");
	generator->add_scalar(p);
}

void Parser::add_function(const std::string &p, int args)
{
	check("function");
	generator->add_function(p, args);
}

void Parser::new_equation()
{
	check("equation");
	generator->new_equation();
}
