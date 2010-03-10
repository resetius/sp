#ifndef GENERATOR_H
#define GENERATOR_H

#include <string>
#include <map>
#include <set>
#include <list>
#include <memory>

struct Generator;
typedef Generator * (*generator_creator_t)(void);
int register_generator(const std::string & name, generator_creator_t);
Generator * make_generator(const std::string & name);

struct Expression
{
	Generator * g;

	enum {
		STRING   = 0,
		NUMBER   = 1,
		OPERATOR = 2,
		LIST     = 3,
		COMPLEX  = 4,
	};
	int type;
	std::string str; // var or operator name
	double num;      // number

	typedef std::list < Expression > stack_t;
	stack_t stack;

	typedef std::list < Expression > list_t;
	list_t sexpr;

	Expression & operator = (const std::string & s) throw();
	Expression & operator = (double v) throw();
	Expression & operator += (const Expression & other) throw();
	Expression & operator -= (const Expression & other) throw();
	Expression & operator /= (const Expression & other) throw();
	Expression & operator *= (const Expression & other) throw();
	Expression & operator ^= (const Expression & other) throw();
	Expression & operator , (const Expression & other) throw();
	Expression & operator () (const Expression & other) throw();
	~Expression();

protected:
	Expression(Generator * p): g(p) {};
	friend struct Generator;
	friend class GC;
};

struct Generator
{
	enum {
		FLOAT = 1,
		DOUBLE = 2,
		BOTH = 3,
	};

	int precission;

	typedef std::set < std::string > scalars_t;
	scalars_t scalars;
	typedef std::set < std::string > initials_t;
	initials_t initials;

	typedef std::map < std::string, int > functions_t;
	functions_t functions;

	Generator();
	virtual ~Generator();

	virtual void set_surface(const std::string &p) = 0;
	virtual void set_approximation(const std::string &p) = 0;
	virtual void set_precission(const std::string & p);

	virtual void add_initial(const std::string &var);
	virtual void add_scalar(const std::string &p);
	virtual void add_function(const std::string &p, int args);

	typedef std::map < std::string , Expression * > declarations_t;
	declarations_t declarations;

	virtual void add_declaration(const std::string &p, Expression * e);
	virtual void add_equation(Expression * lp, Expression * rp);

	virtual void make_header(const std::string & name, const std::string & h_name);
	virtual void make(const std::string & class_name, const std::string & h_name, const std::string & cpp_name);
	virtual void check_var(const std::string &p);

	typedef std::set < std::string > methods_declrs_t;
	methods_declrs_t methods_declrs;
};

class GC
{
	std::list < char * > strings;
	std::list < Expression * > exprs;

public:

	char * new_string(const char * str);
	Expression * new_expression(Generator * gen);
	void collect_all();
};

struct Parser
{
	GC gc;
	std::string name;

	std::auto_ptr < Generator > generator;

	Parser();
	~Parser();

	void set_method(const std::string &p);
	void set_name(const std::string & n);

	void set_precission(const std::string & p);
	void set_surface(const std::string &p);
	void set_approximation(const std::string &p);

	void add_initial(const std::string &p);
	void add_scalar(const std::string &p);
	void add_function(const std::string &p, int args);
	void add_declaration(const std::string &p, Expression * e);

	void add_equation(Expression * lp, Expression * rp);

	void check(const std::string &p);
	void check_var(const std::string &p);
	void make(const std::string & hname, const std::string & cppname);

	Expression * new_expression();
};

#include "generator_harmonic.h"

#endif /* GENERATOR_H */

