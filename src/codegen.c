/*
 * codegen.c -- AST -> LLVM IR, written entirely in C using LLVM's C API
 * (llvm-c/Core.h) rather than the C++ IRBuilder. The whole compiler
 * (lexer, both parsers, AST, codegen, driver, runtime) is now pure C.
 *
 * Satisfies:
 *   Experiment 3 ("Generating code with the LLVM backend") together
 *     with main.c, which takes the LLVMModuleRef produced here through
 *     bitcode -> llc -> object code -> linked executable.
 *   Experiment 8 ("... Generating IR from the AST") -- this file.
 *   Experiment 9 ("Converting types from an AST description to LLVM
 *     types") -- see gs_type_to_llvm() below.
 *
 * Design note: as in the C++ version this replaces, sequence analysis
 * operations lower to `call` instructions into the C runtime library
 * (runtime.c) rather than raw IR arithmetic. Variables (`x =
 * GC_CONTENT(dna);`) get a real `alloca` in the entry block, `store`
 * their computed value, and are `load`ed back out on use -- the
 * standard alloca/load/store pattern for locals before an optimizer's
 * mem2reg pass would promote them to SSA registers.
 */
#include "codegen.h"

#include <llvm-c/Core.h>
#include <llvm-c/Analysis.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* ---- Experiment 9: AST type -> LLVM type ---- */
typedef enum { GS_SEQUENCE, GS_INT32, GS_INT64, GS_DOUBLE } GSType;

static LLVMTypeRef gs_type_to_llvm(GSType t, LLVMContextRef ctx) {
    switch (t) {
        case GS_SEQUENCE: return LLVMPointerTypeInContext(ctx, 0); /* opaque i8*-equivalent */
        case GS_INT32:    return LLVMInt32TypeInContext(ctx);
        case GS_INT64:    return LLVMInt64TypeInContext(ctx);
        case GS_DOUBLE:   return LLVMDoubleTypeInContext(ctx);
    }
    return LLVMVoidTypeInContext(ctx);
}

/* The kind of value a GeneScript variable/expression can hold, tracked
   at compile time -- this is a flat, statically-typed-per-binding
   language, so no runtime type tags are needed. */
typedef enum {
    VT_SEQ_RAW,             /* a declared SEQUENCE, or an aliased copy of one */
    VT_REVERSE,             /* result of REVERSE(x) */
    VT_COMPLEMENT,          /* result of COMPLEMENT(x) */
    VT_REVERSE_COMPLEMENT,  /* result of REVERSE_COMPLEMENT(x) */
    VT_PROTEIN,             /* result of TRANSLATE(x) */
    VT_LENGTH,              /* result of LENGTH(x) -- i64 */
    VT_GC_PERCENT,          /* result of GC_CONTENT(x) -- double */
} GSValType;

static int is_string_like(GSValType t) { return t != VT_LENGTH && t != VT_GC_PERCENT; }

typedef struct VarInfo {
    char *name;
    LLVMValueRef storage; /* global string constant (VT_SEQ_RAW, non-alloca) or an alloca to load from */
    GSValType type;
    int is_alloca;
    struct VarInfo *next; /* newest binding for a name is found first (shadow/reassign) */
} VarInfo;

typedef struct {
    LLVMValueRef val;
    GSValType type;
} ComputedValue;

typedef struct {
    LLVMContextRef ctx;
    LLVMModuleRef module;
    LLVMBuilderRef builder;
    LLVMValueRef main_fn;
    LLVMBasicBlockRef entry_block;
    VarInfo *symbols;

    /* runtime function handles: value + its LLVM function type (needed
       by LLVMBuildCall2, which -- unlike the old getelementptr-typed
       call in the legacy API -- always wants the callee's type explicitly) */
    LLVMValueRef rt_fn[16];
    LLVMTypeRef rt_ty[16];
} CodeGen;

enum {
    RT_CALC_GC_CONTENT, RT_CALC_LENGTH, RT_CALC_REVERSE, RT_CALC_COMPLEMENT,
    RT_CALC_REVERSE_COMPLEMENT, RT_CALC_TRANSLATE,
    RT_PRINT_GC_VALUE, RT_PRINT_LENGTH_VALUE, RT_PRINT_REVERSE_VALUE,
    RT_PRINT_COMPLEMENT_VALUE, RT_PRINT_REVERSE_COMPLEMENT_VALUE,
    RT_PRINT_PROTEIN_VALUE, RT_PRINT_SEQ_RAW,
    RT_PRINT_FIND_MOTIF, RT_PRINT_COMPARE,
    RT_COUNT
};

static void semantic_error(int line, const char *msg) {
    fprintf(stderr, "Semantic Error (line %d): %s\n", line, msg);
    exit(1);
}

static void validate_nucleotides(const ASTNode *decl) {
    const char *seq = decl->strval;
    for (int i = 0; seq[i] != '\0'; i++) {
        char c = (char)toupper((unsigned char)seq[i]);
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
            char msg[160];
            snprintf(msg, sizeof(msg),
                     "Invalid nucleotide '%c' at position %d. Valid symbols: A, T, G, C",
                     seq[i], i + 1);
            semantic_error(decl->line, msg);
        }
    }
}

static void declare_runtime(CodeGen *cg) {
    LLVMTypeRef i8p = gs_type_to_llvm(GS_SEQUENCE, cg->ctx);
    LLVMTypeRef i32 = gs_type_to_llvm(GS_INT32, cg->ctx);
    LLVMTypeRef i64 = gs_type_to_llvm(GS_INT64, cg->ctx);
    LLVMTypeRef dbl = gs_type_to_llvm(GS_DOUBLE, cg->ctx);
    LLVMTypeRef voidTy = LLVMVoidTypeInContext(cg->ctx);

    struct { int id; const char *name; LLVMTypeRef ret; LLVMTypeRef args[2]; unsigned nargs; } specs[] = {
        { RT_CALC_GC_CONTENT,             "gs_calc_gc_content",             dbl,    {i8p},      1 },
        { RT_CALC_LENGTH,                 "gs_calc_length",                 i64,    {i8p},      1 },
        { RT_CALC_REVERSE,                "gs_calc_reverse",                i8p,    {i8p},      1 },
        { RT_CALC_COMPLEMENT,             "gs_calc_complement",             i8p,    {i8p},      1 },
        { RT_CALC_REVERSE_COMPLEMENT,     "gs_calc_reverse_complement",     i8p,    {i8p},      1 },
        { RT_CALC_TRANSLATE,              "gs_calc_translate",              i8p,    {i8p},      1 },
        { RT_PRINT_GC_VALUE,              "gs_print_gc_value",              voidTy, {dbl},      1 },
        { RT_PRINT_LENGTH_VALUE,          "gs_print_length_value",          voidTy, {i64, i32}, 2 },
        { RT_PRINT_REVERSE_VALUE,         "gs_print_reverse_value",         voidTy, {i8p},      1 },
        { RT_PRINT_COMPLEMENT_VALUE,      "gs_print_complement_value",      voidTy, {i8p},      1 },
        { RT_PRINT_REVERSE_COMPLEMENT_VALUE, "gs_print_reverse_complement_value", voidTy, {i8p}, 1 },
        { RT_PRINT_PROTEIN_VALUE,         "gs_print_protein_value",         voidTy, {i8p},      1 },
        { RT_PRINT_SEQ_RAW,               "gs_print_seq_raw",               voidTy, {i8p},      1 },
        { RT_PRINT_FIND_MOTIF,            "gs_print_find_motif",            voidTy, {i8p, i8p}, 2 },
        { RT_PRINT_COMPARE,               "gs_print_compare",               voidTy, {i8p, i8p}, 2 },
    };

    for (size_t i = 0; i < sizeof(specs) / sizeof(specs[0]); i++) {
        LLVMTypeRef fnTy = LLVMFunctionType(specs[i].ret, specs[i].args, specs[i].nargs, 0);
        LLVMValueRef fn = LLVMAddFunction(cg->module, specs[i].name, fnTy);
        cg->rt_fn[specs[i].id] = fn;
        cg->rt_ty[specs[i].id] = fnTy;
    }
}

static LLVMValueRef call1(CodeGen *cg, int rtId, LLVMValueRef arg0) {
    LLVMValueRef args[1] = { arg0 };
    return LLVMBuildCall2(cg->builder, cg->rt_ty[rtId], cg->rt_fn[rtId], args, 1, "");
}

static LLVMValueRef call2(CodeGen *cg, int rtId, LLVMValueRef arg0, LLVMValueRef arg1) {
    LLVMValueRef args[2] = { arg0, arg1 };
    return LLVMBuildCall2(cg->builder, cg->rt_ty[rtId], cg->rt_fn[rtId], args, 2, "");
}

/* allocas go at the top of the entry block, matching -O0 clang codegen
   for a local variable (and making them trivial for a later mem2reg
   pass to promote to SSA registers, if one is ever run). */
static LLVMValueRef create_entry_alloca(CodeGen *cg, LLVMTypeRef ty, const char *name) {
    LLVMBuilderRef tmp = LLVMCreateBuilderInContext(cg->ctx);
    LLVMValueRef firstInstr = LLVMGetFirstInstruction(cg->entry_block);
    if (firstInstr) LLVMPositionBuilderBefore(tmp, firstInstr);
    else LLVMPositionBuilderAtEnd(tmp, cg->entry_block);
    LLVMValueRef a = LLVMBuildAlloca(tmp, ty, name);
    LLVMDisposeBuilder(tmp);
    return a;
}

static VarInfo *symtab_find(CodeGen *cg, const char *name) {
    for (VarInfo *v = cg->symbols; v; v = v->next) {
        if (strcmp(v->name, name) == 0) return v;
    }
    return NULL;
}

static VarInfo *lookup(CodeGen *cg, const char *name, int line) {
    VarInfo *v = symtab_find(cg, name);
    if (!v) {
        char msg[160];
        snprintf(msg, sizeof(msg), "'%s' is not declared.", name);
        semantic_error(line, msg);
    }
    return v;
}

static void symtab_bind(CodeGen *cg, const char *name, LLVMValueRef storage, GSValType type, int is_alloca) {
    VarInfo *v = (VarInfo *)malloc(sizeof(VarInfo));
    v->name = strdup(name);
    v->storage = storage;
    v->type = type;
    v->is_alloca = is_alloca;
    v->next = cg->symbols;
    cg->symbols = v;
}

static GSType llvm_type_for(GSValType t) {
    switch (t) {
        case VT_LENGTH: return GS_INT64;
        case VT_GC_PERCENT: return GS_DOUBLE;
        default: return GS_SEQUENCE; /* all string-like results are pointers */
    }
}

/* Load a variable's usable value (dereferencing its alloca, if it has one). */
static LLVMValueRef gen_load_var(CodeGen *cg, const VarInfo *v, GSType llvmTy) {
    if (!v->is_alloca) return v->storage; /* VT_SEQ_RAW global constant: already usable */
    return LLVMBuildLoad2(cg->builder, gs_type_to_llvm(llvmTy, cg->ctx), v->storage, "loadtmp");
}

static LLVMValueRef gen_string_operand(CodeGen *cg, ASTNode *expr) {
    if (expr->kind == NODE_IDENTIFIER) {
        VarInfo *v = lookup(cg, expr->name, expr->line);
        if (!is_string_like(v->type)) {
            char msg[160];
            snprintf(msg, sizeof(msg), "'%s' does not hold a sequence/string value here.", expr->name);
            semantic_error(expr->line, msg);
        }
        return gen_load_var(cg, v, GS_SEQUENCE);
    }
    if (expr->kind == NODE_STRING_LITERAL) {
        return LLVMBuildGlobalStringPtr(cg->builder, expr->strval, "motif");
    }
    semantic_error(expr->line, "expected a sequence or string here");
    return NULL; /* unreachable */
}

/* Evaluate a FuncCall into a runtime value + its GeneScript type.
   FIND_MOTIF is intentionally excluded -- it's statement-only. */
static ComputedValue gen_compute(CodeGen *cg, ASTNode *call) {
    if (call->kind != NODE_FUNC_CALL || call->func == FUNC_FIND_MOTIF) {
        semantic_error(call->line,
            "FIND_MOTIF (and COMPARE) can only be used as their own statement, "
            "not assigned to a variable or composed into another expression.");
    }
    LLVMValueRef arg0 = gen_string_operand(cg, call->args->node);
    ComputedValue cv;
    switch (call->func) {
        case FUNC_GC_CONTENT:
            cv.val = call1(cg, RT_CALC_GC_CONTENT, arg0); cv.type = VT_GC_PERCENT; return cv;
        case FUNC_LENGTH:
            cv.val = call1(cg, RT_CALC_LENGTH, arg0); cv.type = VT_LENGTH; return cv;
        case FUNC_REVERSE:
            cv.val = call1(cg, RT_CALC_REVERSE, arg0); cv.type = VT_REVERSE; return cv;
        case FUNC_COMPLEMENT:
            cv.val = call1(cg, RT_CALC_COMPLEMENT, arg0); cv.type = VT_COMPLEMENT; return cv;
        case FUNC_REVERSE_COMPLEMENT:
            cv.val = call1(cg, RT_CALC_REVERSE_COMPLEMENT, arg0); cv.type = VT_REVERSE_COMPLEMENT; return cv;
        case FUNC_TRANSLATE:
            cv.val = call1(cg, RT_CALC_TRANSLATE, arg0); cv.type = VT_PROTEIN; return cv;
        default:
            semantic_error(call->line, "unsupported expression");
    }
    cv.val = NULL; cv.type = VT_SEQ_RAW;
    return cv; /* unreachable */
}

static void gen_print_value(CodeGen *cg, ComputedValue cv, int printBinary) {
    switch (cv.type) {
        case VT_GC_PERCENT:
            call1(cg, RT_PRINT_GC_VALUE, cv.val);
            return;
        case VT_LENGTH: {
            LLVMValueRef bin = LLVMConstInt(LLVMInt32TypeInContext(cg->ctx), printBinary ? 1 : 0, 0);
            call2(cg, RT_PRINT_LENGTH_VALUE, cv.val, bin);
            return;
        }
        case VT_REVERSE: call1(cg, RT_PRINT_REVERSE_VALUE, cv.val); return;
        case VT_COMPLEMENT: call1(cg, RT_PRINT_COMPLEMENT_VALUE, cv.val); return;
        case VT_REVERSE_COMPLEMENT: call1(cg, RT_PRINT_REVERSE_COMPLEMENT_VALUE, cv.val); return;
        case VT_PROTEIN: call1(cg, RT_PRINT_PROTEIN_VALUE, cv.val); return;
        case VT_SEQ_RAW: call1(cg, RT_PRINT_SEQ_RAW, cv.val); return;
    }
}

static void gen_assignment(CodeGen *cg, ASTNode *stmt) {
    ASTNode *rhs = stmt->expr;

    if (rhs->kind == NODE_IDENTIFIER) {
        /* x = y;  -- alias: share the existing binding's descriptor.
           Values here are all immutable once produced, so sharing the
           same storage is safe. */
        VarInfo *src = lookup(cg, rhs->name, stmt->line);
        symtab_bind(cg, stmt->name, src->storage, src->type, src->is_alloca);
        return;
    }

    if (rhs->kind != NODE_FUNC_CALL) {
        semantic_error(stmt->line,
            "assignment right-hand side must be a function call (GC_CONTENT, LENGTH, "
            "REVERSE, COMPLEMENT, REVERSE_COMPLEMENT, or TRANSLATE) or another variable.");
    }

    ComputedValue cv = gen_compute(cg, rhs);
    GSType llvmTy = llvm_type_for(cv.type);
    LLVMValueRef slot = create_entry_alloca(cg, gs_type_to_llvm(llvmTy, cg->ctx), stmt->name);
    LLVMBuildStore(cg->builder, cv.val, slot);
    symtab_bind(cg, stmt->name, slot, cv.type, 1);
}

static void gen_print_like(CodeGen *cg, ASTNode *expr, int printBinary) {
    if (expr->kind == NODE_IDENTIFIER) {
        VarInfo *v = lookup(cg, expr->name, expr->line);
        ComputedValue cv;
        cv.val = gen_load_var(cg, v, llvm_type_for(v->type));
        cv.type = v->type;
        gen_print_value(cg, cv, printBinary);
        return;
    }
    if (expr->kind == NODE_STRING_LITERAL) {
        LLVMValueRef s = LLVMBuildGlobalStringPtr(cg->builder, expr->strval, "strlit");
        call1(cg, RT_PRINT_SEQ_RAW, s);
        return;
    }
    if (expr->kind != NODE_FUNC_CALL) {
        semantic_error(expr->line, "expression cannot be printed");
    }
    if (expr->func == FUNC_FIND_MOTIF) {
        LLVMValueRef arg0 = gen_string_operand(cg, expr->args->node);
        LLVMValueRef arg1 = gen_string_operand(cg, expr->args->next->node);
        call2(cg, RT_PRINT_FIND_MOTIF, arg0, arg1);
        return;
    }
    gen_print_value(cg, gen_compute(cg, expr), printBinary);
}

static void gen_statement(CodeGen *cg, ASTNode *stmt) {
    switch (stmt->kind) {
        case NODE_SEQUENCE_DECL: {
            validate_nucleotides(stmt);
            LLVMValueRef g = LLVMBuildGlobalStringPtr(cg->builder, stmt->strval, stmt->name);
            symtab_bind(cg, stmt->name, g, VT_SEQ_RAW, 0);
            return;
        }
        case NODE_ASSIGNMENT:
            gen_assignment(cg, stmt);
            return;
        case NODE_PRINT:
            gen_print_like(cg, stmt->expr, stmt->print_binary);
            return;
        case NODE_EXPR_STMT:
            gen_print_like(cg, stmt->expr, stmt->print_binary);
            return;
        case NODE_COMPARE: {
            VarInfo *l = lookup(cg, stmt->left, stmt->line);
            VarInfo *r = lookup(cg, stmt->right, stmt->line);
            call2(cg, RT_PRINT_COMPARE,
                  gen_load_var(cg, l, GS_SEQUENCE), gen_load_var(cg, r, GS_SEQUENCE));
            return;
        }
        default: {
            char msg[128];
            snprintf(msg, sizeof(msg), "cannot generate code for a top-level %s", ast_kind_name(stmt->kind));
            semantic_error(stmt->line, msg);
        }
    }
}

LLVMModuleRef codegen_program(ASTNode *root, LLVMContextRef ctx, const char *moduleName) {
    CodeGen cg;
    memset(&cg, 0, sizeof(cg));
    cg.ctx = ctx;
    cg.module = LLVMModuleCreateWithNameInContext(moduleName, ctx);
    cg.builder = LLVMCreateBuilderInContext(ctx);
    declare_runtime(&cg);

    LLVMTypeRef mainTy = LLVMFunctionType(LLVMInt32TypeInContext(ctx), NULL, 0, 0);
    cg.main_fn = LLVMAddFunction(cg.module, "main", mainTy);
    cg.entry_block = LLVMAppendBasicBlockInContext(ctx, cg.main_fn, "entry");
    LLVMPositionBuilderAtEnd(cg.builder, cg.entry_block);

    for (ASTList *cur = root->statements; cur; cur = cur->next) {
        gen_statement(&cg, cur->node);
    }

    LLVMBuildRet(cg.builder, LLVMConstInt(LLVMInt32TypeInContext(ctx), 0, 0));

    char *err = NULL;
    if (LLVMVerifyModule(cg.module, LLVMReturnStatusAction, &err)) {
        fprintf(stderr, "Internal error: generated invalid LLVM IR:\n%s\n", err);
        LLVMDisposeMessage(err);
        exit(1);
    }
    if (err) LLVMDisposeMessage(err);

    LLVMDisposeBuilder(cg.builder);
    /* symbol table nodes are process-lifetime; freed on exit */
    return cg.module;
}
