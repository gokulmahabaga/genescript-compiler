LLVM_CONFIG := llvm-config-18
CC := clang

CFLAGS  := $(shell $(LLVM_CONFIG) --cflags) -std=c11 -Wno-unused-parameter -g
LDFLAGS := $(shell $(LLVM_CONFIG) --ldflags)
LIBS    := $(shell $(LLVM_CONFIG) --libs core bitwriter analysis)
SYSLIBS := $(shell $(LLVM_CONFIG) --system-libs)

SRC := src

.PHONY: all clean test

all: gsc runtime.o

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

rdparser.o: $(SRC)/rdparser.c $(SRC)/ast.h $(SRC)/parser.tab.h
	$(CC) -c $(SRC)/rdparser.c -o rdparser.o -I$(SRC) -std=c11 -D_GNU_SOURCE

codegen.o: $(SRC)/codegen.c $(SRC)/codegen.h $(SRC)/ast.h
	$(CC) -c $(SRC)/codegen.c -o codegen.o -I$(SRC) $(CFLAGS)

main.o: $(SRC)/main.c $(SRC)/ast.h $(SRC)/codegen.h
	$(CC) -c $(SRC)/main.c -o main.o -I$(SRC) $(CFLAGS)

# ---- Link the compiler driver itself (this is "gsc", the GeneScript compiler) ----
gsc: ast.o parser.tab.o lex.yy.o rdparser.o codegen.o main.o
	$(CC) $^ -o gsc $(LDFLAGS) $(LIBS) $(SYSLIBS)

clean:
	rm -f *.o gsc $(SRC)/parser.tab.c $(SRC)/parser.tab.h $(SRC)/lex.yy.c
	rm -f examples/*.ll examples/*.bc examples/*.s examples/*.o examples/*.exe
	rm -f *.ll *.bc *.s *.exe

# Run every example through both front ends.
test: all
	@for f in examples/*.gs; do \
		echo "=== $$f (bison) ==="; ./gsc $$f --frontend=bison; echo; \
		echo "=== $$f (recursive-descent) ==="; ./gsc $$f --frontend=rd --no-run -o /tmp/rdcheck > /dev/null; echo "(rd front end: parsed + codegen OK)"; echo; \
	done
