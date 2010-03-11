%{
struct Parser;
#include "generator_parser.hpp"

#define YY_DECL int yylex( yyscan_t yyscanner, YYSTYPE & yylval )

int isatty (int )
{
	return 0;
}

%}
%option 8bit
%option outfile="generator_scanner.cpp" header-file="generator_scanner.h" case-insensitive
%option reentrant
%%
[ \t]+ { continue; }
\\.*\n { ++yylineno; continue; }
"name" { return NAME; }
"precission" { return PRECISSION; }
"method" { return METHOD; }
"surface" { return SURFACE; }
"approximation" { return APPROXIMATION; }
"initial" { return START; }
"declare" { return DECLARE; }
"scalar" { return SCALAR; }
"function" { return FUNCTION; }
"equation" { return EQUATION; }
"Delta" { return DELTA; }
"Jacobian" { return JACOBIAN; }
"->" { return ARROW; }
[0-9]+\.?[0-9]* { yylval.num = atof(yytext); return NUMBER; }
[a-zA-Z]+[0-9a-zA-Z]*'? { yylval.str = strdup(yytext); return VAR; }
"\n"    { ++yylineno; return '\n'; }
. { /*fprintf(stderr, " %c ", yytext[0]);*/return yytext[0]; }
%%
/**/
