#include <stdio.h>

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
	for (declarations_t::iterator it = declarations.begin();
		it != declarations.end(); ++it)
	{
		delete it->second;
	}
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
		throw logic_error(str.str());
	}
}

void Generator::add_initial(const std::string & p)
{
	fprintf(stderr, "new initial %s\n", p.c_str());

	initials.insert(p);
}

void Generator::add_scalar(const std::string &p)
{
	scalars.insert(p);
}

void Generator::add_function(const std::string &p, int args)
{
	if (args > 3 || args < 2) {
		stringstream str;
		str << p << ": args < 2 || args > 3 ";
		throw logic_error(str.str());
	}
	fprintf(stderr, "new function %s\n", p.c_str());

	functions[p] = args;
}

void Generator::add_declaration(const std::string &p, Expression * e)
{
	if (declarations.find(p) != declarations.end()) {
		stringstream str;
		str << p << " already declared";
		throw logic_error(str.str());
	}

	declarations[p] = e;
}

void Generator::check_var(const std::string &p)
{
	bool flag = false;
	flag |= (declarations.find(p) != declarations.end());
	flag |= (initials.find(p) != initials.end());

	if (!flag) {
		stringstream str;
		str << "symbol " << p << " is not defined";
		throw logic_error(str.str());
	}
}

void Generator::add_equation(Expression * lp, Expression * rp)
{
	fprintf(stderr, "new equation ");
}

Expression *  Generator::new_expression()
{
	Expression * exp = new Expression(this);
	return exp;
}

void Generator::make(const string & name, const string & h_name, const string & cpp_name)
{
	make_header(name, h_name);
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

void Parser::add_initial(const std::string & p)
{
	check("initial");
	generator->add_initial(p);
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

void Parser::add_equation(Expression * lp, Expression * rp)
{
	check("equation");
	generator->add_equation(lp, rp);
}

void Parser::add_declaration(const std::string &p, Expression * e)
{
	check("declaration");
	generator->add_declaration(p, e);
}

void Parser::make(const std::string & hname, const std::string & cppname)
{
	check("make");
	if (name.empty()) {
		throw logic_error("name not set !");
	}
	generator->make(name, hname, cppname);
}

void Parser::check_var(const std::string &p)
{
	generator->check_var(p);
}

Expression::~Expression()
{
}

Expression * Parser::new_expression()
{
	check("new expression");
	Expression * exp = generator->new_expression();
	return exp;
}


