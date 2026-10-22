# GeneScript Compiler

A DNA-focused programming language and compiler that translates GeneScript programs into executable code using LLVM, with DNA-specific optimization techniques and a multithreaded batch-processing subsystem.

## Overview

GeneScript is a compiler project designed to explore programming language implementation, compiler optimization, and DNA sequence processing. It combines a custom language front end with LLVM-based code generation and a runtime for DNA-related operations.

The project includes a main compiler pipeline and a separate batch-processing executable. The batch executable uses POSIX threads to process multiple DNA sequences concurrently.

## Key Features

- **Custom language front end:** Supports parsing GeneScript source code using the project's parser implementations.
- **DNA-specific operations:** Supports operations such as sequence length, GC content, reverse, complement, reverse complement, and translation in the batch-processing subsystem.
- **Compiler optimization:** Includes DNA-oriented algebraic rewrites and local value numbering/common subexpression elimination (CSE).
- **LLVM code generation:** Uses LLVM infrastructure to generate machine code through the compiler toolchain.
- **Parallel batch processing:** Uses POSIX threads (`pthreads`) to distribute sequence-processing tasks among worker threads.
- **Testing:** Includes batch-processing tests to help validate the executable's behavior.

## Architecture

The project can be understood as two related workflows.

### 1. Main compiler pipeline

1. **Source input:** Read a GeneScript program.
2. **Parsing:** Convert source code into a structured representation.
3. **Optimization:** Apply supported compiler optimization passes.
4. **Code generation:** Generate LLVM-based intermediate representation or code through the project's code-generation components.
5. **Machine-code generation:** Use the configured LLVM and compiler tools to produce an executable.

The main compiler pipeline currently performs its parsing, optimization, and code-generation stages sequentially.

### 2. DNA batch-processing pipeline

1. Read a batch of DNA sequences and requested operations.
2. Create a shared work queue.
3. Distribute sequence-processing tasks among worker threads.
4. Execute the requested operations using the runtime.
5. Emit results in the original sequence order.

The batch runner uses a mutex to coordinate access to the shared work queue.

## Parallelism

GeneScript currently implements **data-level parallelism in its DNA batch-processing subsystem**.

The batch runner in `src/batch_main.c` uses POSIX threads:

- It defines a maximum of four participating worker threads, including the main thread.
- Worker threads claim sequence-processing tasks from a shared queue.
- A mutex protects access to the shared task index.
- Results are stored by sequence index so they can be emitted in input order.
- If mutex initialization fails, the implementation falls back to serial processing.

**Important distinction:** This is parallel execution of DNA sequence-processing tasks, not parallel compilation of multiple GeneScript source programs. The main compiler pipeline does not currently implement parallel parsing, optimization, or code generation.

## DNA Operations

The batch-processing subsystem supports the following operation names:

| Operation | Description |
|---|---|
| `LENGTH` | Calculate the length of a DNA sequence. |
| `GC_CONTENT` | Calculate the proportion of G and C bases. |
| `REVERSE` | Reverse the sequence. |
| `COMPLEMENT` | Generate the complementary DNA sequence. |
| `REVERSE_COMPLEMENT` | Generate the reverse complement. |
| `TRANSLATE` | Translate a DNA sequence into its corresponding protein representation, as implemented by the runtime. |

The operations available to a complete GeneScript program depend on the language and runtime implementation.

## Compiler Optimizations

The project includes optimization techniques intended to reduce unnecessary computation:

- **DNA-specific algebraic rewrites:** Apply supported simplifications involving DNA-related expressions.
- **Local value numbering:** Identify equivalent computations within the relevant scope.
- **Common subexpression elimination:** Reuse the results of equivalent expressions where the optimizer can establish that this is valid.

Optimization behavior and effectiveness depend on the program and the supported optimization rules.

## Project Structure

```text
genescript-compiler/
├── src/
│   ├── main.c              # Main compiler entry point
│   ├── batch_main.c        # Multithreaded DNA batch runner
│   └── codegen.c           # Code-generation implementation
├── tests/
│   └── test_batch.py       # Batch-processing tests
├── Makefile                # Build and test targets
├── README.md               # Project documentation
└── ...
```

Other source files, parser components, runtime files, and optimization modules may be present in the repository.

## Requirements

Install the dependencies required by your environment and project configuration:

- GCC or another compatible C compiler
- GNU Make
- POSIX threads support (`pthread`)
- LLVM and the configured LLVM toolchain
- Clang
- Flex and Bison, if building the corresponding parser components
- Python 3 for the batch test script

Exact LLVM versions and additional dependencies should match the project's build configuration.

## Build Instructions

Clone the repository if you do not already have it:

```bash
git clone https://github.com/gokulmahabaga/genescript-compiler.git
cd genescript-compiler
```

Build the project using Make:

```bash
make
```

To build the batch-processing executable separately:

```bash
make genescript-batch
```

The batch executable is expected to be generated as `genescript-batch` in the project root when the target succeeds.

## Run the Batch Processor

The batch runner accepts JSON input containing an `operations` array and a `sequences` array.

For example:

```bash
printf '%s\n' '{"operations":["LENGTH","GC_CONTENT","REVERSE_COMPLEMENT"],"sequences":["ATGCGT","ACGT","TTT","GGCATC"]}' | ./genescript-batch
```

This example requests three operations for each of four DNA sequences.

The exact output format depends on the implementation in `src/batch_main.c`.

## Run Tests

Build the batch executable and run its test target:

```bash
make test-batch
```

You can also run the Python test script directly after building:

```bash
python3 tests/test_batch.py
```

A passing test suite helps verify functional correctness. To evaluate parallel speedup, benchmark the threaded implementation against a serial implementation using the same workloads and repeated measurements.

## Current Limitations

- The main compiler pipeline is sequential rather than a parallel compiler.
- Parallel processing is implemented in the DNA batch runner, not across compiler stages.
- The available DNA operations and language features are limited to those implemented in the source code and runtime.
- Performance improvements from multithreading depend on the workload, available CPU cores, and thread-management overhead.
- Performance claims should be supported by reproducible benchmarks.

## Future Improvements

- Investigate parallel compilation of independent source files or compilation units.
- Explore parallel optimization for independent functions or intermediate-representation units.
- Add reproducible serial-versus-parallel performance benchmarks.
- Expand compiler diagnostics and error reporting.
- Add more language examples and regression tests.
- Document supported language syntax and optimization rules in greater detail.

## Learning Outcomes

This project provides practical experience with:

- Compiler architecture and language parsing
- Intermediate representations and LLVM-based code generation
- Compiler optimization techniques
- DNA sequence processing
- POSIX threads and synchronization
- Build automation and test-driven development

## Repository

GitHub: https://github.com/gokulmahabaga/genescript-compiler

## License

Refer to the repository's `LICENSE` file for the applicable license terms.
