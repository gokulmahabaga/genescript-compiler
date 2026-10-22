LLVM_CONFIG := llvm-config-18
CC := clang

CFLAGS  := $(shell $(LLVM_CONFIG) --cflags) -std=c11 -Wno-unused-parameter -g
LDFLAGS := $(shell $(LLVM_CONFIG) --ldflags)
LIBS    := $(shell $(LLVM_CONFIG) --libs core bitwriter analysis)
SYSLIBS := $(shell $(LLVM_CONFIG) --system-libs)

SRC := src

.PHONY: all clean test test-batch

all: gsc runtime.o genescript-batch

# ---- Flex / Bison generated files (Lab Experiments 1, 6, 7) ----
$(SRC)/parser.tab.c $(SRC)/parser.tab.h: $(SRC)/parser.y
	bison -d -o $(SRC)/parser.tab.c $(SRC)/parser.y

$(SRC)/lex.yy.c: $(SRC)/lexer.l $(SRC)/parser.tab.h
	flex -o $(SRC)/lex.yy.c $(SRC)/lexer.l

# ---- Every compiled unit here is C (this whole compiler is pure C) ----
ast.o: $(SRC)/ast.c $(SRC)/ast.h
	$(CC) -c $(SRC)/ast.c -o ast.o -I$(SRC) -std=c11 -D_GNU_SOURCE

parser.tab.o: $(SRC)/parser.tab.c $(SRC)/ast.h
	$(CC) -c $(SRC)/parser.tab.c -o parser.tab.o -I$(SRC) -std=c11 -D_GNU_SOURCE

lex.yy.o: $(SRC)/lex.yy.c $(SRC)/parser.tab.h $(SRC)/ast.h
	$(CC) -c $(SRC)/lex.yy.c -o lex.yy.o -I$(SRC) -std=c11 -D_GNU_SOURCE -Wno-unused-function

runtime.o: $(SRC)/runtime.c
	$(CC) -c $(SRC)/runtime.c -o runtime.o -O2 -std=c11 -D_GNU_SOURCE

# Standalone JSON batch runner; does not alter the gsc compiler target.
genescript-batch: $(SRC)/batch_main.c runtime.o
	$(CC) $(SRC)/batch_main.c runtime.o -o genescript-batch -std=c11 -O2 -pthread

rdparser.o: $(SRC)/rdparser.c $(SRC)/ast.h $(SRC)/parser.tab.h
	$(CC) -c $(SRC)/rdparser.c -o rdparser.o -I$(SRC) -std=c11 -D_GNU_SOURCE

optimizer.o: $(SRC)/optimizer.c $(SRC)/optimizer.h $(SRC)/dump.h $(SRC)/ast.h
	$(CC) -c $(SRC)/optimizer.c -o optimizer.o -I$(SRC) -std=c11 -D_GNU_SOURCE

dump.o: $(SRC)/dump.c $(SRC)/dump.h $(SRC)/ast.h $(SRC)/parser.tab.h
	$(CC) -c $(SRC)/dump.c -o dump.o -I$(SRC) -std=c11 -D_GNU_SOURCE

codegen.o: $(SRC)/codegen.c $(SRC)/codegen.h $(SRC)/ast.h
	$(CC) -c $(SRC)/codegen.c -o codegen.o -I$(SRC) $(CFLAGS)

main.o: $(SRC)/main.c $(SRC)/ast.h $(SRC)/codegen.h $(SRC)/optimizer.h $(SRC)/dump.h
	$(CC) -c $(SRC)/main.c -o main.o -I$(SRC) $(CFLAGS)

# ---- Link the compiler driver itself (this is "gsc", the GeneScript compiler) ----
gsc: ast.o parser.tab.o lex.yy.o rdparser.o optimizer.o dump.o codegen.o main.o
	$(CC) $^ -o gsc $(LDFLAGS) $(LIBS) $(SYSLIBS)

clean:
	rm -f *.o gsc genescript-batch $(SRC)/parser.tab.c $(SRC)/parser.tab.h $(SRC)/lex.yy.c
	rm -f examples/*.ll examples/*.bc examples/*.s examples/*.o examples/*.exe
	rm -f *.ll *.bc *.s *.exe

# Run every example + test case through BOTH front ends and compare
# against the saved .expected output; check that tests/errors/ programs
# are rejected. See tests/run_tests.sh.
test: all
	@tests/run_tests.sh

test-batch: genescript-batch
	python3 tests/test_batch.py
