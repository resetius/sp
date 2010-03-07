#ifndef GENERATOR_H
#define GENERATOR_H

#include <string>
struct Parser
{
	bool error;
	std::string name;

	enum {
		FLOAT = 1,
		DOUBLE = 2,
		BOTH = 3,
	};

	int precission;

	enum {
		FEM = 1,
		FDM = 2,
		HARMONIC = 3,
	};

	int method;

	enum {
		FLAT = 1,
		SPHERE = 2,
	};

	int surface;

	Parser();
	~Parser();

	void set_name(const std::string & n);
	void set_precission(const std::string & p);
	void set_method(const std::string &p);
	void set_surface(const std::string &p);
	void set_approximation(const std::string &p);

	void add_map(const std::string & from, const std::string & to);
	void add_scalar(const std::string &p);
	void add_function(const std::string &p, int args);

	void new_equation();
};

#endif /* GENERATOR_H */
