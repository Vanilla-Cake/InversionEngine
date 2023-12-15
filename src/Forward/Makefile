
# GCC=gcc

# MPICC=mpicc

# FLAGS = -ld_classic -lm #if not macOS, remove -ld_classic

INCLUDE = -I../include

LIB = -L../lib -lcore


forwardGBM : forwardGBM.o 
	$(GCC) -o forwardGBM forwardGBM.o $(FLAGS) $(INCLUDE) $(LIB)

forwardGBM.o : forwardGBM.c
	$(GCC) -c forwardGBM.c $(INCLUDE)

clean:
	rm forwardGBM *.o

all: forwardGBM

install: forwardGBM
	cp forwardGBM ../../bin
