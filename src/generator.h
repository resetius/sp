#ifndef GENERATOR_H
#define GENERATOR_H

#include <string>
#include <map>
#include <set>
#include <memory>

struct Generator;
typedef Generator * (*generator_creator_t)(void);
int register_generator(const std::string & name, generator_creator_t);
Generator * make_generator(const std::string & name);

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

	typedef std::map < std::string, int > functions_t;
	functions_t functions;

	Generator();
	virtual ~Generator();

	virtual void set_surface(const std::string &p) = 0;
	virtual void set_approximation(const std::string &p) = 0;
	virtual void set_precission(const std::string & p);

	virtual void add_map(const std::string & from, const std::string & to);
	virtual void add_scalar(const std::string &p);
	virtual void add_function(const std::string &p, int args);

	virtual void new_equation();

	virtual void make_header(const std::string & name, const std::string & h_name);
	virtual void make(const std::string & class_name, const std::string & h_name, const std::string & cpp_name);
};

struct Parser
{
	std::string name;

	std::auto_ptr < Generator > generator;

	Parser();
	~Parser();

	void set_method(const std::string &p);
	void set_name(const std::string & n);

	void set_precission(const std::string & p);
	void set_surface(const std::string &p);
	void set_approximation(const std::string &p);

	void add_map(const std::string & from, const std::string & to);
	void add_scalar(const std::string &p);
	void add_function(const std::string &p, int args);

	void new_equation();

	void check(const std::string &p);
	void make(const std::string & hname, const std::string & cppname);
};

#include "generator_harmonic.h"

#endif /* GENERATOR_H */
