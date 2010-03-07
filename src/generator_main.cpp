#include <stdlib.h>
#include <stdexcept>
#include "generator.h"

#include "generator_scanner.h"
#include "generator_parser.hpp"

using namespace std;

int yyparse(Parser * ptr);

void usage(char * argv)
{
	fprintf(stderr, "usage: %s filename\n", argv);
	exit(-1);
}

int yywrap()
{
	return 1;
}

void do_all(FILE * f)
{
	Parser p;
	yyset_in(f);

	try {
		while (yyparse(&p));
	} catch (exception & e) {
		fprintf(stderr, "%s\n", e.what());
	}
}

int main(int argc, char * argv[])
{
	if (argc < 2) {
		usage(argv[0]);
	}

	FILE * f = fopen(argv[1], "rb");
	if (!f) {
		usage(argv[0]);
	}

	do_all(f);
	fclose(f);

	return 0;
}
