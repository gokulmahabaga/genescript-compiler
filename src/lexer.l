/*
 * lexer.l -- Flex scanner for GeneScript.
 *
 * Satisfies Lab Experiment 1 ("Implementation of LEXR [lexer]") and,
 * together with parser.y, Experiment 7 ("Intro to Flex and Bison").
 *
 * Experiment 7's specific exercise -- terminating a statement with
 * "; b" instead of ";" to print output in binary -- is implemented
 * below as the SEMICOLON_BIN rule. Flex's maximal-munch rule means
 * ";b" is preferred over ";" automatically whenever the 'b' is present
 * with no whitespace in between.
 */
%option noyywrap
%option yylineno

%{
#include <string.h>
#include "ast.h"
#include "parser.tab.h"

static char *copy_unquoted(const char *text, int len) {
    /* strip the surrounding double quotes */
    char *s = (char *)malloc(len - 1);
    memcpy(s, text + 1, len - 2);
    s[len - 2] = '\0';
    return s;
}
%}

%%

[ \t\r\n]+              { /* skip whitespace */ }
"//"[^\n]*              { /* skip line comments */ }

"SEQUENCE"              { return SEQUENCE; }
"PRINT"                 { return PRINT; }
"COMPARE"               { return COMPARE; }
"WITH"                  { return WITH; }
"GC_CONTENT"            { return GC_CONTENT; }
"LENGTH"                { return LENGTH; }
"FIND_MOTIF"            { return FIND_MOTIF; }
"REVERSE_COMPLEMENT"    { return REVERSE_COMPLEMENT; }
"REVERSE"               { return REVERSE; }
"COMPLEMENT"            { return COMPLEMENT; }
"TRANSLATE"             { return TRANSLATE; }

\"[^\"]*\"              { yylval.str = copy_unquoted(yytext, yyleng); return STRING; }

[A-Za-z_][A-Za-z0-9_]*  { yylval.str = strdup(yytext); return IDENTIFIER; }

"="                     { return ASSIGN; }
"("                     { return LPAREN; }
")"                     { return RPAREN; }
","                     { return COMMA; }
";b"                    { return SEMICOLON_BIN; }
";"                     { return SEMICOLON; }

.                       { fprintf(stderr, "Lex Error (line %d): unexpected character '%s'\n", yylineno, yytext); }

%%
