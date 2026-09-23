# GeneScript Compiler — LLVM / Flex / Bison edition

Built specifically for **BCSE307P (Compiler Design Lab)**, whose indicative
experiments require an LLVM-, Flex-, and Bison-based toolchain rather than a
hand-rolled interpreter. (There's also a pure-Python version of this same
GeneScript DSL — same language, same group project — built for the
lecture-style J-component; this folder is the separate lab deliverable.)

## What this is

A real compiler for GeneScript (a DSL for DNA sequence analysis: GC content,
motif search, reverse-complement, translation, mutation comparison), written
entirely in **C** (no C++ anywhere in the authored source — LLVM is driven
through its C API, `llvm-c/Core.h`, not the C++ `IRBuilder`). It has:

- **Two independent front ends** that parse the same grammar into the same
  AST: a Flex+Bison (LALR) front end, and a hand-written recursive-descent
  front end.
- **Real LLVM IR codegen** using LLVM's C API (`llvm-c/Core.h`), not the C++ IRBuilder.
- **Real machine code output**: LLVM bitcode → assembly (`.s`) → object
  code (`.o`) → a linked, runnable native executable — via LLVM's own `llc`
  backend and `clang` as the linker, not a toy VM.

## Build & run

Requires `flex`, `bison`, `clang`, and `llvm-18` (specifically
`llvm-config-18`, `llc`). All installed via `apt-get install flex bison
clang llvm llvm-dev`.

```bash
make              # builds ./gsc (the compiler driver) and runtime.o
./gsc examples/basic.gs                      # Bison front end (default)
./gsc examples/basic.gs --frontend=rd        # recursive-descent front end
./gsc examples/basic.gs --emit-llvm          # also print IR to stdout
./gsc examples/basic.gs --no-run             # compile/link but don't execute
make test         # runs every examples/*.gs through both front ends
make clean
```

Each run of `gsc` leaves behind `<name>.ll` (IR text), `<name>.bc`
(bitcode), `<name>.s` (assembly), `<name>.o` (object code), and
`<name>.exe` (linked binary) in the current directory — inspect any of
them to see the real intermediate artifacts.

## Mapping to the 10 lab experiments

| # | Experiment (from the syllabus) | Where it's satisfied |
|---|---|---|
| 1 | Implementation of LEXR using LLVM | `src/lexer.l` — Flex-generated scanner; tokens feed both front ends, which drive LLVM codegen |
| 2 | Handwritten parser using LLVM | `src/rdparser.c` — hand-written top-down parser, calls the same `codegen_program()` |
| 3 | Generating code with the LLVM backend | `src/main.c` — takes the `LLVMModuleRef` from codegen through bitcode → `llc` → object code → linked executable |
| 4 | Defining a real programming language | `GRAMMAR.md` — full CFG + lexical grammar + static semantic rules |
| 5 | Recursive descent parser for the CFG, implemented using LLVM | `src/rdparser.c` implements exactly the grammar in `GRAMMAR.md` and feeds LLVM codegen |
| 6 | LR parser for the CFG, implemented using LLVM | `src/parser.y` — Bison LALR(1) grammar, same AST/codegen path |
| 7 | Intro to Flex and Bison; the `"; b"` binary-output exercise | `src/lexer.l` + `src/parser.y`; the exercise itself is documented in `GRAMMAR.md` and demoed in `examples/binary_exercise.gs` |
| 8 | LLVM-style RTTI for the AST; generating IR from the AST | `src/ast.h` (`NodeKind` tag + `ast_isa()`, the classic LLVM RTTI idiom) and `src/codegen.c` |
| 9 | Converting AST types to LLVM types | `src/codegen.c::gs_type_to_llvm()` |
| 10 | Emitting assembler text and object code | `src/main.c` — `llc <bc> -o .s` and `llc -filetype=obj <bc> -o .o` |

## For teammates building on top of this

This repo is **just the compiler** (the BCSE307P lab deliverable) — the
group's ML classification and any Python/UI layer are meant to sit on top
of it, not inside it. A few things worth knowing before you start:

- **No `PREDICT`/`ANALYZE` keyword yet.** The Python interpreter version of
  GeneScript (a separate part of this project, not in this repo) has a
  placeholder GC-content-heuristic classifier wired in as `PREDICT(seq)`.
  This LLVM compiler doesn't have that keyword at all yet — it only
  supports `GC_CONTENT`, `LENGTH`, `FIND_MOTIF`, `REVERSE`, `COMPLEMENT`,
  `REVERSE_COMPLEMENT`, `TRANSLATE`. If the ML model needs to be callable
  from `.gs` source directly, that means: add a `PREDICT` token
  (`src/lexer.l`), a grammar rule (`src/parser.y` + `src/rdparser.c`), a
  `FUNC_PREDICT` case (`src/ast.h`, `src/codegen.c`), and a runtime function
  it calls into (`src/runtime.c`) — follow the same pattern any of the
  existing functions use, they're all one shape.
- **Easiest ML integration point**: `src/runtime.c` is a plain C file that
  gets linked into every compiled `.gs` program's executable. The
  lowest-friction way to plug in a model is a C function there (or one that
  calls out to a Python process / a C-exported model via a shared library)
  — no LLVM/codegen knowledge required to add it, just the same
  `PREDICT`-keyword plumbing above so the language can call it.
- **Python integration**: if the "python part" means calling this compiler
  from Python (e.g. a wrapper, a UI, or glue code) rather than extending the
  language, the cleanest boundary is: `gsc` takes a `.gs` file and produces
  a `.exe`; run that as a subprocess (`subprocess.run(["./gsc", path])`) and
  capture its stdout. No need to touch the C code at all for that.
- **`FIND_MOTIF` still isn't assignable to a variable** (see below) — if the
  ML model needs motif-search results as structured data rather than
  printed text, that's the other likely place someone will need to extend
  `codegen.c`.

## Scope note (read before extending)

Variables are supported: `x = GC_CONTENT(dna);` allocates a real stack slot
(`alloca`) in the generated LLVM IR, stores the computed value into it, and
later reads (`PRINT x;`, or using `x` as another function's argument) `load`
it back out — the same alloca/load/store pattern clang itself emits for a
local variable at `-O0`. A variable's type (int, double, or sequence
pointer) is tracked at compile time from how it was produced, so:

```
SEQUENCE dna = "ATGCGATCGATCG";
rc = REVERSE_COMPLEMENT(dna);   -- rc is a sequence, stored via alloca
FIND_MOTIF(rc, "CGA");          -- variables can feed into other functions
alias = dna;                    -- x = y aliases an existing binding
PRINT alias;
```

**Deliberately not supported**: assigning `FIND_MOTIF`'s result to a
variable (`m = FIND_MOTIF(...)`) — it's a composite result (a motif string
plus a variable-length list of positions), and this backend doesn't have a
struct/array runtime type to hold it in. It's still fully usable as its own
statement (`FIND_MOTIF(dna, "ATG");`). `COMPARE` is similarly statement-only,
matching the grammar (it was never an `expr` to begin with). Both give a
clear semantic error if you try. If you need `FIND_MOTIF` results
composable, that's a reasonable next step — it'd need a small runtime
struct (`{char* motif; long* positions; long count;}`) returned by pointer.

## File-by-file

```
src/ast.h, ast.c        AST node definitions + constructors (LLVM-style RTTI)
src/lexer.l             Flex scanner
src/parser.y            Bison LALR(1) grammar
src/rdparser.c        Hand-written recursive-descent parser (2nd front end)
src/codegen.h/.c       AST -> LLVM IR
src/runtime.c           C runtime: GC content, motif search, translate,
                         reverse-complement, mutation comparison, printing
src/main.c            Driver: parse -> codegen -> IR/bitcode -> llc -> clang
GRAMMAR.md              Formal CFG + lexical grammar + semantic rules
Makefile                Builds gsc + runtime.o; `make test` runs all examples
examples/*.gs           Sample programs
```
