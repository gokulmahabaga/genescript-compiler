# GeneScript — Language Specification (Lab Experiment 4)

"Defining a real programming language" — this is GeneScript's formal
context-free grammar, as implemented in `src/parser.y` (Bison) and
`src/rdparser.cpp` (hand-written recursive descent). Both front ends
parse exactly this grammar and build the same AST (`src/ast.h`).

## Lexical grammar (tokens, see `src/lexer.l`)

```
IDENTIFIER          [A-Za-z_][A-Za-z0-9_]*        (excluding keywords below)
STRING              "\"" [^"]* "\""
KEYWORDS            SEQUENCE PRINT COMPARE WITH
                    GC_CONTENT LENGTH FIND_MOTIF
                    REVERSE_COMPLEMENT REVERSE COMPLEMENT TRANSLATE
SYMBOLS             = ( ) , ;
SEMICOLON_BIN       ";b"      -- see "Binary output" exercise below
COMMENT             "//" .* (to end of line, skipped)
WHITESPACE          [ \t\r\n]+ (skipped)
```

## Context-free grammar

```
program          -> statement*

statement        -> SEQUENCE IDENTIFIER '=' STRING terminator
                   | PRINT expr terminator
                   | COMPARE IDENTIFIER WITH IDENTIFIER ';'
                   | expr terminator

terminator       -> ';' | ';b'

expr             -> IDENTIFIER
                   | STRING
                   | GC_CONTENT           '(' expr ')'
                   | LENGTH               '(' expr ')'
                   | REVERSE              '(' expr ')'
                   | COMPLEMENT           '(' expr ')'
                   | REVERSE_COMPLEMENT   '(' expr ')'
                   | TRANSLATE            '(' expr ')'
                   | FIND_MOTIF           '(' expr ',' expr ')'
```

## Static (semantic) rules

- A `SEQUENCE` declaration's string literal must contain only the
  characters `A`, `T`, `G`, `C` (case-insensitive). Any other character
  is a semantic error reporting the exact 1-indexed position.
- Every `IDENTIFIER` used as a sequence operand must have been declared
  by an earlier `SEQUENCE` statement in the same program (no forward
  references, no scoping — this is a flat, single-pass language).
- `FIND_MOTIF`'s second argument is expected to be a string literal
  (the motif to search for); the grammar doesn't prevent passing an
  identifier there, but codegen will only accept it if that identifier
  is itself a declared sequence (used as a literal search string).

## "Binary output" exercise (Lab Experiment 7)

> Modify the scanner and parser so that terminating a statement with
> "; b" instead of ";" results in the output being printed in binary.

Implemented exactly as specified:

- `src/lexer.l` adds a `";b"` rule (`SEMICOLON_BIN` token). Flex's
  maximal-munch rule means `;b` wins over plain `;` automatically
  whenever the `b` immediately follows the semicolon with no space.
- `src/parser.y` and `src/rdparser.cpp` both accept either terminator
  on `PRINT` statements and bare expression statements, setting a
  `print_binary` flag on the AST node.
- `src/codegen.cpp` passes that flag through to `gs_print_length` in
  the runtime (`src/runtime.c`), which prints the integer in binary
  instead of decimal when the flag is set.

Example (`examples/binary_exercise.gs`):

```
SEQUENCE dna = "ATGCGATCGATCG";
PRINT LENGTH(dna);b
```

Output: `Length: 1101 (binary)` (13 in binary).

The flag is only wired up for `LENGTH` in this implementation (the
only integer-valued result in the language); it's accepted but has no
effect on sequence-valued or report-valued prints, since "binary" only
has an obvious meaning for a number. This is documented in `README.md`.
