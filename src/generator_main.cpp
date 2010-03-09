#include <stdlib.h>
#include <stdexcept>
#include <string>
#include "generator.h"

#include "generator_scanner.h"
#include "generator_parser.hpp"

using namespace std;

int yyparse(Parser * ptr);

void usage(char * argv)
{
	fprintf(stderr, "usage: %s \n"
			"-in filename\n"
			"-header filename.h\n"
			"-cpp filename.cpp\n"
			, argv);
	exit(-1);
}

int yywrap()
{
	return 1;
}

void do_all(FILE * f, const string & hname, const string & cppname)
{
	Parser p;
	yyset_in(f);

	try {
		while (yyparse(&p));
	} catch (exception & e) {
		fprintf(stderr, "%s\n", e.what());
	}

	try {
		p.make(hname, cppname);
	} catch (exception & e) {
		fprintf(stderr, "%s\n", e.what());
	}
}

int main(int argc, char * argv[])
{
	string input;
	string hname;
	string cppname;

	for (int i = 0; i < argc; ++i) {
		if (!strcmp(argv[i], "-h")) {
			usage(argv[0]);
		} else if (!strcmp(argv[i], "-in")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}

			input = argv[i + 1];
		} else if (!strcmp(argv[i], "-header")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}

			hname = argv[i + 1];
		} else if (!strcmp(argv[i], "-cpp")) {
			if (i == argc - 1) {
				usage(argv[0]);
			}
			cppname = argv[i + 1];
		}
	}

	FILE * f = fopen(input.c_str(), "rb");
	if (!f) {
		fprintf(stderr, "cannot open '%s'\n", input.c_str());
		usage(argv[0]);
	}

	if (hname.empty()) {
		hname = input + ".h";
	}

	if (cppname.empty()) {
		cppname = input + ".cpp";
	}

	do_all(f, hname, cppname);
	fclose(f);

	return 0;
}

