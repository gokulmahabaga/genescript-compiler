/*
 * parser.y -- LALR(1) grammar for GeneScript, built with Bison.
 *
 * Satisfies Lab Experiment 6 ("Write a LR parser for the CFG language
 * ... implement it using LLVM") and Experiment 7 ("Intro to Flex and
 * Bison"). The formal CFG this implements is written out in GRAMMAR.md.
 *
 * Experiment 7 also asks: "Modify the scanner and parser so that
 * terminating a statement with '; b' instead of ';' results in the
 * output being printed in binary." That's implemented here as the
 * SEMICOLON_BIN token (see lexer.l) -- statements ending in ";b" set
 * print_binary=1 on the resulting AST node, and codegen/runtime print
 * the numeric result in binary instead of decimal.
 */
%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ast.h"

extern int yylex(void);
extern int yylineno;
void yyerror(const char *msg);

ASTNode *gs_ast_root = NULL;
%}

%union {
    ASTNode *node;
    ASTList *list;
    char *str;
}

%token <str> IDENTIFIER STRING
%token SEQUENCE PRINT COMPARE WITH
%token GC_CONTENT LENGTH FIND_MOTIF REVERSE_COMPLEMENT REVERSE COMPLEMENT TRANSLATE
%token ASSIGN LPAREN RPAREN COMMA SEMICOLON SEMICOLON_BIN

%type <node> statement expr
%type <list> stmt_list

%start program

%%

program
    : stmt_list { gs_ast_root = ast_new_program($1); }
    ;

stmt_list
    : /* empty */          { $$ = NULL; }
    | statement stmt_list  { $$ = ast_list_prepend($1, $2); }
    ;

statement
    : SEQUENCE IDENTIFIER ASSIGN STRING SEMICOLON {
          $$ = ast_new_sequence_decl($2, $4, yylineno);
          free($2); free($4);
      }
    | IDENTIFIER ASSIGN expr SEMICOLON {
          $$ = ast_new_assignment($1, $3, yylineno);
          free($1);
      }
    | PRINT expr SEMICOLON      { $$ = ast_new_print($2, 0, yylineno); }
    | PRINT expr SEMICOLON_BIN  { $$ = ast_new_print($2, 1, yylineno); }
    | COMPARE IDENTIFIER WITH IDENTIFIER SEMICOLON {
          $$ = ast_new_compare($2, $4, yylineno);
          free($2); free($4);
      }
    | expr SEMICOLON {
          $$ = ast_new_expr_stmt($1, yylineno);
      }
    | expr SEMICOLON_BIN {
          $$ = ast_new_expr_stmt($1, yylineno);
          $$->print_binary = 1;
      }
    ;

expr
    : IDENTIFIER { $$ = ast_new_identifier($1, yylineno); free($1); }
    | STRING     { $$ = ast_new_string_literal($1, yylineno); free($1); }
    | GC_CONTENT LPAREN expr RPAREN {
          $$ = ast_new_func_call(FUNC_GC_CONTENT, ast_list_append_single($3), yylineno);
      }
    | LENGTH LPAREN expr RPAREN {
          $$ = ast_new_func_call(FUNC_LENGTH, ast_list_append_single($3), yylineno);
      }
    | REVERSE LPAREN expr RPAREN {
          $$ = ast_new_func_call(FUNC_REVERSE, ast_list_append_single($3), yylineno);
      }
    | COMPLEMENT LPAREN expr RPAREN {
          $$ = ast_new_func_call(FUNC_COMPLEMENT, ast_list_append_single($3), yylineno);
      }
    | REVERSE_COMPLEMENT LPAREN expr RPAREN {
          $$ = ast_new_func_call(FUNC_REVERSE_COMPLEMENT, ast_list_append_single($3), yylineno);
      }
    | TRANSLATE LPAREN expr RPAREN {
          $$ = ast_new_func_call(FUNC_TRANSLATE, ast_list_append_single($3), yylineno);
      }
    | FIND_MOTIF LPAREN expr COMMA expr RPAREN {
          $$ = ast_new_func_call(FUNC_FIND_MOTIF,
                                  ast_list_prepend($3, ast_list_append_single($5)),
                                  yylineno);
      }
    ;

%%

void yyerror(const char *msg) {
    fprintf(stderr, "Parse Error (Bison, line %d): %s\n", yylineno, msg);
}
