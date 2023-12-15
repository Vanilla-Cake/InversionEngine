
# GCC=gcc

# MPICC=mpicc

# FLAGS = -ld_classic -lm #if not macOS, remove -ld_classic

PAR_OBJ =  atopkge.o errpkge.o getpars.o subcalls.o syscalls.o

PAR = atopkge.c errpkge.c getpars.c subcalls.c syscalls.c

CORE_OBJ = alloc.o complex.o franuni.o fft.o sinc.o

CORE = alloc.c complex.c franuni.c fft.c sinc.c

INCLUDE = -I../include

LIB = -L../lib


libcore.a : $(CORE_OBJ) $(PAR_OBJ)
	ar r libcore.a $(CORE_OBJ) $(PAR_OBJ)
	ranlib libcore.a

# core
$(CORE_OBJ) : $(CORE)
	$(GCC) -c $(CORE) $(INCLUDE)

# par
$(PAR_OBJ) : $(PAR)
	$(GCC) -c $(PAR) $(INCLUDE)


clean:
	rm *.o *.a

install: libcore.a
	cp libcore.a ../lib


