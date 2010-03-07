%{
#include "generator_parser.hpp"
%}
%option 8bit
%option outfile="generator_scanner.c" header-file="generator_scanner.h" case-insensitive
%%
[ \t]+ { continue; }
\\.*\n { continue; }
"name" { return NAME; }
"precission" { return PRECISSION; }
"method" { return METHOD; }
"surface" { return SURFACE; }
"approximation" { return APPROXIMATION; }
"map" { return MAP; }
"scalar" { return SCALAR; }
"function" { return FUNCTION; }
"equation" { return EQUATION; }
"Delta" { return DELTA; }
"Jacobian" { return JACOBIAN; }
"->" { return ARROW; }
[0-9]+\.?[0-9]* { return NUMBER; }
[a-zA-Z]+[0-9a-zA-Z]*'? { return VAR; }
"\n"    { ++yylineno; return '\n'; }
. { return yytext[0]; }
%%
/**/
