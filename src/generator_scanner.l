%{
#include "generator.h"
#include "generator_parser.hpp"

#define YY_DECL int yylex( YYSTYPE * yylval, yyscan_t yyscanner, Parser * parser )

int isatty (int )
{
	return 0;
}

char * my_strcpy(char *dest, const char *src){
               size_t i;

               for (i = 0 ; src[i] != '\0' ; i++)
                   dest[i] = src[i];
               dest[i] = '\0';

               return dest;
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
[0-9]+\.?[0-9]* { yylval->num = atof(yytext); return NUMBER; }
[a-zA-Z]+[0-9a-zA-Z]*'? { 
	yylval->str = (char*)calloc(128, 1); 
fprintf(stderr, "%s\n", yytext);
char * p = yytext;
//	my_strcpy(yylval->str, p);/*"harmonic";*/ /*strdup(yytext); */ /*parser->gc.new_string(yytext);*/ 
yylval->str = parser->gc.new_string(p);
	return VAR; 
}
"\n"    { ++yylineno; return '\n'; }
. { /*fprintf(stderr, " %c ", yytext[0]);*/return yytext[0]; }
%%
/**/
