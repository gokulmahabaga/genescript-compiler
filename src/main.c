/*
 * main.c -- GeneScript compiler driver ("gsc"), pure C.
 *
 * Pipeline:
 *   source.gs --[flex+bison OR hand-written recursive-descent]--> AST
 *            --[codegen.c, LLVM-C API]--> LLVM Module
 *            --[this file]--> .ll (text IR) + .bc (bitcode)
 *            --[llc, LLVM's backend]--> .s (assembly) + .o (object code)
 *            --[clang]--> linked executable (linked against runtime.o)
 *
 * This ties together Lab Experiments 1-3, 5-7, 9, 10 in one runnable
 * tool. See README.md for the full experiment-by-experiment mapping.
 *
 * Usage:
 *   gsc <file.gs> [--frontend=bison|rd] [--emit-llvm] [--no-run] [-o NAME]
 */
#include "ast.h"
#include "codegen.h"

#include <llvm-c/Core.h>
#include <llvm-c/BitWriter.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern FILE *yyin;
extern int yyparse(void);
extern ASTNode *gs_ast_root;
extern ASTNode *rd_parse_program(void);

static void strip_ext(const char *path, char *out, size_t outsz) {
    const char *slash = strrchr(path, '/');
    const char *base = slash ? slash + 1 : path;
    const char *dot = strrchr(base, '.');
    size_t len = dot ? (size_t)(dot - base) : strlen(base);
    if (len >= outsz) len = outsz - 1;
    memcpy(out, base, len);
    out[len] = '\0';
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <file.gs> [--frontend=bison|rd] [--emit-llvm] [--no-run] [-o NAME]\n", argv[0]);
        return 1;
    }

    const char *inputFile = NULL;
    const char *frontend = "bison";
    int emitLLVM = 0;
    int runAfter = 1;
    char outName[512] = "";

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (strncmp(a, "--frontend=", 11) == 0) frontend = a + 11;
        else if (strcmp(a, "--emit-llvm") == 0) emitLLVM = 1;
        else if (strcmp(a, "--no-run") == 0) runAfter = 0;
        else if (strcmp(a, "-o") == 0 && i + 1 < argc) { strncpy(outName, argv[++i], sizeof(outName) - 1); }
        else if (a[0] != '-') inputFile = a;
    }

    if (!inputFile) {
        fprintf(stderr, "Error: no input file given\n");
        return 1;
    }

    yyin = fopen(inputFile, "r");
    if (!yyin) {
        fprintf(stderr, "Error: could not open '%s'\n", inputFile);
        return 1;
    }

    ASTNode *root = NULL;
    if (strcmp(frontend, "bison") == 0) {
        if (yyparse() != 0) {
            fprintf(stderr, "Compilation failed (Bison/LALR front end).\n");
            return 1;
        }
        root = gs_ast_root;
    } else if (strcmp(frontend, "rd") == 0) {
        root = rd_parse_program();
    } else {
        fprintf(stderr, "Unknown --frontend value '%s' (use bison or rd)\n", frontend);
        return 1;
    }
    fclose(yyin);

    if (!root) {
        fprintf(stderr, "Compilation failed: no AST produced.\n");
        return 1;
    }

    char base[512];
    strip_ext(inputFile, base, sizeof(base));
    if (outName[0] == '\0') strncpy(outName, base, sizeof(outName) - 1);

    LLVMContextRef ctx = LLVMContextCreate();
    char moduleName[560];
    snprintf(moduleName, sizeof(moduleName), "%s_module", base);
    LLVMModuleRef module = codegen_program(root, ctx, moduleName);

    if (emitLLVM) {
        char *ir = LLVMPrintModuleToString(module);
        fputs(ir, stdout);
        LLVMDisposeMessage(ir);
    }

    char llPath[600], bcPath[600], sPath[600], oPath[600], exePath[600];
    snprintf(llPath, sizeof(llPath), "%s.ll", base);
    snprintf(bcPath, sizeof(bcPath), "%s.bc", base);
    snprintf(sPath, sizeof(sPath), "%s.s", base);
    snprintf(oPath, sizeof(oPath), "%s.o", base);
    snprintf(exePath, sizeof(exePath), "%s.exe", outName);

    char *err = NULL;
    if (LLVMPrintModuleToFile(module, llPath, &err)) {
        fprintf(stderr, "Could not write %s: %s\n", llPath, err ? err : "unknown error");
        if (err) LLVMDisposeMessage(err);
        return 1;
    }
    if (LLVMWriteBitcodeToFile(module, bcPath) != 0) {
        fprintf(stderr, "Could not write %s\n", bcPath);
        return 1;
    }

    fprintf(stderr, "[gsc] wrote LLVM IR:      %s\n", llPath);
    fprintf(stderr, "[gsc] wrote LLVM bitcode: %s\n", bcPath);

    /* Hand the bitcode to LLVM's own backend (llc) to emit real
       assembly text and a real object file -- Experiment 10. */
    char cmd[2048];
    snprintf(cmd, sizeof(cmd), "llc %s -o %s", bcPath, sPath);
    if (system(cmd) != 0) { fprintf(stderr, "llc (assembly) failed\n"); return 1; }
    snprintf(cmd, sizeof(cmd), "llc -filetype=obj %s -o %s", bcPath, oPath);
    if (system(cmd) != 0) { fprintf(stderr, "llc (object) failed\n"); return 1; }
    fprintf(stderr, "[gsc] wrote assembly:     %s\n", sPath);
    fprintf(stderr, "[gsc] wrote object code:  %s\n", oPath);

    /* Link against the runtime library and (optionally) run it. */
    snprintf(cmd, sizeof(cmd), "clang -no-pie %s runtime.o -o %s", oPath, exePath);
    if (system(cmd) != 0) {
        fprintf(stderr, "Linking failed (expected runtime.o to already be built -- see Makefile)\n");
        return 1;
    }
    fprintf(stderr, "[gsc] linked executable:  %s\n\n", exePath);

    LLVMDisposeModule(module);
    LLVMContextDispose(ctx);

    if (runAfter) {
        snprintf(cmd, sizeof(cmd), "./%s", exePath);
        system(cmd);
    }

    return 0;
}
