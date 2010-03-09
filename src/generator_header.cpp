#include <stdexcept>
#include <sstream>

#include <ctype.h>
#include <stdio.h>

#include "generator.h"

using namespace std;

static string upper(const string & p)
{
	string r;
	for (size_t i = 0; i < p.size(); ++i) {
		r += (char)toupper(p[i]);
	}
	return r;
}

void Generator::make_header(const string & name, const string & h_name)
{
	FILE * f = fopen(h_name.c_str(), "wb");
	if (!f) {
		stringstream str;
		str << "cannot open " << h_name << " for wrighting";
		throw runtime_error(str.str());
	}

	fprintf(f, 
"#ifndef %s\n"
"#define %s\n"
"/* generated file ! do not edit !*/ \n\n"
"#include <string.h>\n\n"
"#include <memory>\n\n"
, upper(name).c_str(), upper(name).c_str());

	fprintf(f,
"struct %sConf\n"
"{\n"
"	typedef vector < double > array_t; \n"
"	typedef double (*fun_xy_t ) (double x, double y, %sConf * conf ); \n"
"	typedef double (*fun_xyt_t ) (double x, double y, double t, %sConf * conf ); \n\n"
, name.c_str(), name.c_str(), name.c_str()
	);

	for (scalars_t::iterator it = scalars.begin(); it != scalars.end(); ++it)
	{
fprintf(f, "	double %s;\n", it->c_str());
	}

fprintf(f, "\n");

	for (functions_t::iterator it = functions.begin(); it != functions.end(); ++it)
	{
		if (it->second == 3) {
fprintf(f, "	fun_xyt_t %s; \n", it->first.c_str());
		} else if (it->second == 2) {
fprintf(f, "	fun_xy_t f_%s; \n", it->first.c_str());
fprintf(f, "	array_t a_%s; \n", it->first.c_str());
		}
	}

fprintf(f, "\n");
fprintf(f, "	%sConf() { memset(this, 0, sizeof(*this)); }\n", name.c_str());

fprintf(f,
"};\n\n");

fprintf(f,
"class %s\n"
"{\n"
"	class PImpl;\n"
"	std::auto_ptr < PImpl > impl;\n\n"
"public:\n"
"	%s(const %sConf & conf);\n\n"
, name.c_str(), name.c_str(), name.c_str()
);

	for (methods_declrs_t::iterator it = methods_declrs.begin(); it != methods_declrs.end(); ++it)
	{
fprintf(f, "	%s;", it->c_str());
	}

fprintf(f, "};\n\n");
	fprintf(f,
"#endif /*%s*/\n\n", upper(name).c_str());
	fclose(f);
}

