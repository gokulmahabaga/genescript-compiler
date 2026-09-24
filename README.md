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

- **An optimizer** that runs on the AST between parsing and codegen:
  common-subexpression elimination by local value numbering (the DAG
  of a basic block) plus DNA-specific algebraic rewrites like
  `REVERSE(COMPLEMENT(x))` → `REVERSE_COMPLEMENT(x)`. See below.
- **Phase-by-phase dumps** (`--dump-tokens`, `--dump-ast`, `--dump-tac`,
  `--dump-symtab`, `--show-opt`) showing the output of every stage.
- **Two independent front ends** that parse the same grammar into the same
  AST: a Flex+Bison (LALR) front end, and a hand-written recursive-descent
  front end.
- **Real LLVM IR codegen** using LLVM's C API (`llvm-c/Core.h`), not the C++ IRBuilder.
- **Real machine code output**: LLVM bitcode → assembly (`.s`) → object
  code (`.o`) → a linked, runnable native executable — via LLVM's own `llc`
  backend and `clang` as the linker, not a toy VM.

## Build & run

Requires `flex`, `bison`, `clang`, and LLVM 18 **including its
development headers** (`llvm-config-18`, `llc`, `llvm-c/*.h`). On
Ubuntu / WSL:

```bash
sudo apt install -y build-essential flex bison clang llvm-18 llvm-18-dev
```

(Without `llvm-18-dev` the build fails with `'llvm-c/Types.h' file not found`.)

```bash
make              # builds ./gsc (the compiler driver) and runtime.o
./gsc examples/basic.gs                      # Bison front end (default)
./gsc examples/basic.gs --frontend=rd        # recursive-descent front end
./gsc examples/basic.gs --emit-llvm          # also print IR to stdout
./gsc examples/basic.gs --no-run             # compile/link but don't execute
./gsc examples/basic.gs --no-opt             # skip the optimizer (also -O0)
make test         # checks both front ends against saved expected output,
                  # and checks that every tests/errors/*.gs is rejected
make clean
```

`gsc` finds `runtime.o` next to its own executable, so it can be run
from any directory. It exits non-zero on any lexical, parse, or
semantic error. Each run of `gsc` leaves behind `<name>.ll` (IR text), `<name>.bc`
(bitcode), `<name>.s` (assembly), `<name>.o` (object code), and
`<name>.exe` (linked binary) in the current directory — inspect any of
them to see the real intermediate artifacts.

## Web UI (Streamlit)

A browser UI that shows every phase of a compile side by side: token
table, AST, symbol table, TAC before/after optimization with the list of
rewrites, LLVM IR, assembly, and the program's output. If compilation
fails, the pipeline bar marks the phase that failed (lexer, parser, or
semantic analysis) and shows the error.

```bash
make                                   # the UI runs the real ./gsc
python3 -m venv .venv && source .venv/bin/activate
pip install -r ui/requirements.txt
streamlit run ui/app.py                # then open http://localhost:8501
```

`ui/app.py` only runs `gsc` with its `--dump-*` flags and lays out the
result -- it doesn't reimplement any part of the compiler.

## Seeing every phase

```bash
./gsc tests/cases/optimizer_demo.gs --dump-tokens --dump-ast \
      --dump-tac --show-opt --dump-symtab
```

| Flag | Phase | Prints |
|---|---|---|
| `--dump-tokens` | Lexical analysis | every token: line, token type, lexeme |
| `--dump-ast` | Parsing | the syntax tree, drawn as a tree |
| `--dump-tac` | Intermediate code | three-address code ("Bio-IR"), before **and** after optimization |
| `--show-opt` | Optimization | each rewrite applied, with its line number |
| `--dump-symtab` | Semantic analysis | name, type, length, GC%, line, and source of every binding |
| `--emit-llvm` | Code generation | the LLVM IR |

All of these go to stdout, before the program's own output. With
`--no-run` you get only the dumps -- this is what a UI should call.

## The optimizer (`src/optimizer.c`)

Runs on the AST, after parsing and before LLVM codegen. Disable it with
`--no-opt` / `-O0`.

**Pass 1 -- domain-specific algebraic simplification.** DNA identities
that let the compiler do less work:

| Before | After | Why it's valid |
|---|---|---|
| `REVERSE(COMPLEMENT(x))`, `COMPLEMENT(REVERSE(x))` | `REVERSE_COMPLEMENT(x)` | one specialised runtime call instead of two passes + an intermediate copy |
| `LENGTH(REVERSE(x))` (or `COMPLEMENT`, `REVERSE_COMPLEMENT`) | `LENGTH(x)` | these transforms never change the length |
| `GC_CONTENT(REVERSE(x))` (same three) | `GC_CONTENT(x)` | complement swaps G↔C, so the G+C count is unchanged |
| `REVERSE(REVERSE(x))`, `COMPLEMENT(COMPLEMENT(x))`, `REVERSE_COMPLEMENT(REVERSE_COMPLEMENT(x))` | `x` | each is its own inverse (applied only when the value feeds another call, so printed labels don't change) |

**Pass 2 -- common-subexpression elimination (local value numbering).**
A GeneScript program is one basic block, so this is the classic DAG
construction for a block: every expression gets a *value number* built
from its operator and its operands' value numbers. If the same value is
computed more than once:

- and a variable already holds it, later uses read that variable
  (`x = GC_CONTENT(dna); y = GC_CONTENT(dna);` → `y = x;`);
- otherwise the compiler inserts a temporary (`$cse1`, `$cse2`, …),
  computes it once there, and every occurrence reads the temporary.

Reassigning a variable gives it a new value number, so a value is never
reused after the variable holding it changes
(`tests/cases/optimizer_safety.gs` checks this). Aliases (`a = dna;`)
share a value number, so `LENGTH(a)` and `LENGTH(dna)` are recognised
as the same computation.

The optimizer never changes a program's results. The one visible
difference: `PRINT REVERSE(COMPLEMENT(dna));` is labelled
"Reverse Complement:" instead of "Reverse:", because that is what it
now computes.

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
PRINT LENGTH(REVERSE(COMPLEMENT(dna)));   -- calls can be nested
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
src/optimizer.h/.c     AST optimizer: domain rewrites + CSE (value numbering)
src/dump.h/.c          --dump-tokens / --dump-ast / --dump-tac / --dump-symtab
src/codegen.h/.c       AST -> LLVM IR
src/runtime.c           C runtime: GC content, motif search, translate,
                         reverse-complement, mutation comparison, printing
src/main.c            Driver: parse -> optimize -> codegen -> IR/bitcode -> llc -> clang
GRAMMAR.md              Formal CFG + lexical grammar + semantic rules
Makefile                Builds gsc + runtime.o; `make test` runs the test suite
examples/*.gs           Sample programs (+ .expected output for each)
tests/run_tests.sh      Test runner: both front ends vs. expected output
tests/cases/*.gs        Extra programs that must compile and run correctly
                        (<name>.flags adds gsc options, e.g. --dump-tac)
ui/app.py               Streamlit web UI (runs gsc, shows every phase)
ui/requirements.txt     Python packages for the UI
tests/errors/*.gs       Programs that must be rejected with an error
```
