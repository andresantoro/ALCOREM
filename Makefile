CC=gcc
CFLAGS=-O2 
LDLIBS=-lm -lgmp -lz
NAME=$(shell uname -s)
BIN_DIR1= Complexity_Multiplex
BIN_DIR2= Reducibility_Multiplex


CC_RUNS=$(shell \$($(CC) --version >/dev/null ))

##$(info $$CC_RUNS is [${CC_RUNS}])

#ifeq ($(CC_RUNS), "")
ifdef CC_RUNS
IS_CC_GCC=$(shell $(CC) --version 2>/dev/null | grep -Eic '(GCC)')
##$(info $$IS_CC_GCC is [${IS_CC_GCC}])
ifneq ($(IS_CC_GCC), 0)
# compiler is GCC
CFLAGS+=-fopenmp
else
IS_CC_ICC=$(shell $(CC) --version 2>/dev/null | grep -Eic 'Intel')
##$(info $$IS_CC_ICC is [${IS_CC_ICC}])
ifneq ($(IS_CC_ICC), 0)
# compiler is Intel icc
CFLAGS+=-qopenmp
else
IS_CC_CLANG=$(shell $(CC) --version 2>/dev/null | grep -Eic '(llvm|clang)')
##$(info $$IS_CC_CLANG is [${IS_CC_CLANG}])
ifeq ($(IS_CC_CLANG), 1)
# compiler is CLANG
LDLIBS+=-lomp
CFLAGS+=-Xpreprocessor -fopenmp
# ifeq ($(NAME),Darwin)
# 	CFLAGS+=-Xpreprocessor -fopenmp
# endif
endif
endif
endif
endif
all:  complexity reducibility

clean:
	rm -f $(BIN_DIR1)/compute_complexity_multiplex $(BIN_DIR2)/reducibility_complexity

complexity: $(BIN_DIR1)/Compute_Complexity_Multiplex.c
	$(CC) $(CFLAGS) $^ $(LDLIBS) -o $(BIN_DIR1)/compute_complexity_multiplex


reducibility: $(BIN_DIR2)/reducibility_complexity.c
	$(CC) $(CFLAGS) $^ $(LDLIBS) -o $(BIN_DIR2)/reducibility_complexity



