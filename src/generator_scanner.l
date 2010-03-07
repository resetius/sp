%{
#include "generator_parser.hpp"
%}
%option 8bit
%option outfile="generator_scanner.c" header-file="generator_scanner.h" case-insensitive
WHITESPACE [ \t]
%%
"name" { return NAME; }
"precission" { return PRECISSION; }
"method" { return METHOD; }
"surface" { return SURFACE; }
"approximation" { return APPROXIMATION; }
"map" { return MAP; }
"scalar" { return SCALAR; }
"function" { return FUNCTION; }
"Delta" { return DELTA; }
"Jacobian" { return JACOBIAN; }
"->" { return ARROW; }
[0-9a-bA-B]+ { return VAR; }
WHITESPACE {}
"\n"    { ++yylineno; return '\n'; }
. { return yytext[0]; }
%%
/**/
