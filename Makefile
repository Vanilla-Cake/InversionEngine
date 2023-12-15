GCC=gcc

MPICC=mpicc

FLAGS = -ld_classic -lm -O3 #if not macOS, remove -ld_classic

MODULES = \
	src/Core/ \
	src/Migration/GaussianBeam \
	src/Migration/Kirchhoff \
	src/Tools \
	src/Forward


export GCC MPICC FLAGS

all:
	$(foreach module, $(MODULES), $(MAKE) -C $(module);)

install:
	$(foreach module, $(MODULES), $(MAKE) -C $(module) install;)

clean:
	$(foreach module, $(MODULES), $(MAKE) -C $(module) clean;)



	


