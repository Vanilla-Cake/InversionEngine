
# GCC=gcc

# MPICC=mpicc

# FLAGS = -ld_classic -lm #if not macOS, remove -ld_classic

INCLUDE = -I../../include

LIB = -L../../lib -lcore


migrationGBM : migrationGBM.o 
	$(MPICC) -o migrationGBM migrationGBM.o $(FLAGS) $(INCLUDE) $(LIB)

lsm : lsm.o 
	$(GCC) -o lsm lsm.o $(FLAGS) $(INCLUDE) $(LIB)

migrationGBS : migrationGBS.o 
	$(MPICC) -o migrationGBS migrationGBS.o $(FLAGS) $(INCLUDE) $(LIB)

migrationGBM.o : migrationGBM.c
	$(MPICC) -c migrationGBM.c $(INCLUDE)
migrationGBS.o : migrationGBS.c
	$(MPICC) -c migrationGBS.c $(INCLUDE) 

lsm.o : lsm.c
	$(GCC) -c lsm.c $(INCLUDE) 




clean:
	rm migrationGBM lsm migrationGBS *.o

all: migrationGBM lsm migrationGBS

install: migrationGBM lsm migrationGBS
	cp  migrationGBM lsm migrationGBS ../../../bin
