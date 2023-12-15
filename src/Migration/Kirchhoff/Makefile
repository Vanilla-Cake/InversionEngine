
# GCC=gcc

# MPICC=mpicc

# FLAGS = -ld_classic -lm #if not macOS, remove -ld_classic

INCLUDE = -I../../include

LIB = -L../../lib -lcore


migrationKirch : fmm.o kirchhoff.o
	$(GCC) -o migrationKirch fmm.o kirchhoff.o $(FLAGS) $(INCLUDE) $(LIB)
fmm.o : fmm.c
	$(GCC) -c fmm.c $(INCLUDE)

kirchhoff.o : kirchhoff.c
	$(GCC) -c kirchhoff.c $(INCLUDE)


clean:
	rm migrationKirch *.o

all: migrationKirch

install: migrationKirch
	cp  migrationKirch ../../../bin
