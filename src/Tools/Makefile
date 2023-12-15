
# GCC=gcc

# MPICC=mpicc

# FLAGS = -ld_classic -lm #if not macOS, remove -ld_classic

INCLUDE = -I../include

LIB = -L../lib -lcore


makehead : makehead.o 
	$(GCC) -o makehead makehead.o $(FLAGS) $(INCLUDE) $(LIB)
subtract : subtract.o 
	$(GCC) -o subtract subtract.o $(FLAGS) $(INCLUDE) $(LIB)
generate : generate.o 
	$(GCC) -o generate generate.o $(FLAGS) $(INCLUDE) $(LIB)

# head
makehead.o : makehead.c
	$(GCC) -c makehead.c $(INCLUDE) 

# lsm
subtract.o : subtract.c
	$(GCC) -c subtract.c $(INCLUDE) 


# tools
generate.o : generate.c
	$(GCC) -c generate.c $(INCLUDE) 


clean:
	rm makehead subtract generate  *.o

all: makehead subtract generate

install: makehead subtract generate
	cp  makehead subtract generate ../../bin
