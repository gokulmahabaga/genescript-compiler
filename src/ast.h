/*
 * ast.h -- GeneScript Abstract Syntax Tree
 *
 * Follows the classic LLVM RTTI idiom (see "Using RTTI in LLVM" in the
 * LLVM Programmer's Manual): every node carries a `Kind` tag set once at
 * construction, and callers use isa()/cast() instead of dynamic_cast /
 * C++ RTTI. This header is valid C and C++ so it can be included
 * unmodified by the Bison grammar (parser.y, plain C actions), the
 * hand-written recursive-descent parser (rdparser.cpp), and the LLVM
 * IR codegen (codegen.cpp).
 *
 * Satisfies Lab Experiment 8: "Using LLVM-style RTTI for the AST".
 */
#ifndef GENESCRIPT_AST_H
#define GENESCRIPT_AST_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    NODE_PROGRAM,
    NODE_SEQUENCE_DECL,   /* SEQUENCE name = "STRING"; */
    NODE_ASSIGNMENT,      /* name = expr; (expr's value is stored for later use) */
    NODE_PRINT,           /* PRINT expr ; or PRINT expr ;b (binary form) */
    NODE_EXPR_STMT,       /* bare function-call statement, auto-printed */
    NODE_COMPARE,         /* COMPARE a WITH b; */
    NODE_IDENTIFIER,
    NODE_STRING_LITERAL,
    NODE_FUNC_CALL,       /* GC_CONTENT(x), LENGTH(x), FIND_MOTIF(x,y), ... */
} NodeKind;

typedef enum {
    FUNC_GC_CONTENT,
    FUNC_LENGTH,
    FUNC_FIND_MOTIF,
    FUNC_REVERSE_COMPLEMENT,
    FUNC_REVERSE,
    FUNC_COMPLEMENT,
    FUNC_TRANSLATE,
} FuncKind;

typedef struct ASTNode ASTNode;

/* Simple intrusive singly-linked list used for statement lists / arg lists. */
typedef struct ASTList {
    ASTNode *node;
    struct ASTList *next;
} ASTList;

struct ASTNode {
    NodeKind kind;
    int line;

    /* NODE_PROGRAM */
    ASTList *statements;

    /* NODE_SEQUENCE_DECL: name = value (validated A/T/G/C only) */
    /* NODE_IDENTIFIER: name */
    char *name;

    /* NODE_STRING_LITERAL: value */
    char *strval;

    /* NODE_PRINT / NODE_EXPR_STMT: expr (+ binary output flag, exercise 7) */
    ASTNode *expr;
    int print_binary; /* set when the statement was terminated with ";b" */

    /* NODE_COMPARE: left WITH right */
    char *left;
    char *right;

    /* NODE_FUNC_CALL */
    FuncKind func;
    ASTList *args;
};

/* ---- construction helpers ---- */
ASTNode *ast_new_program(ASTList *statements);
ASTNode *ast_new_sequence_decl(const char *name, const char *value, int line);
ASTNode *ast_new_assignment(const char *name, ASTNode *expr, int line);
ASTNode *ast_new_print(ASTNode *expr, int print_binary, int line);
ASTNode *ast_new_expr_stmt(ASTNode *expr, int line);
ASTNode *ast_new_compare(const char *left, const char *right, int line);
ASTNode *ast_new_identifier(const char *name, int line);
ASTNode *ast_new_string_literal(const char *value, int line);
ASTNode *ast_new_func_call(FuncKind func, ASTList *args, int line);

ASTList *ast_list_prepend(ASTNode *node, ASTList *rest);
ASTList *ast_list_append_single(ASTNode *node); /* start a new list */
ASTList *ast_list_reverse(ASTList *list);

/* ---- LLVM-style RTTI: isa<T>-equivalent check ---- */
static inline int ast_isa(const ASTNode *n, NodeKind k) {
    return n != NULL && n->kind == k;
}

const char *ast_kind_name(NodeKind k);
const char *ast_func_name(FuncKind f);

#ifdef __cplusplus
}
#endif

#endif /* GENESCRIPT_AST_H */
