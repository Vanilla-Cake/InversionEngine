# Makefile 

#--------------------------------------------------------
# edit here:

# source code for model generation
# MODEL = model_ainos_visc.c
# MODEL_EL = model_ainos.c
# MODEL_AC = model_acoustic.c
# MODEL_VTI = zinc_vti.c
# MODEL_TTI = zinc_tti.c

# Compiler (LAM: CC=hcc, CRAY T3E: CC=cc)

# ON Linux cluster running LAM
#CC=hcc
#LFLAGS=-lm -lmpi 
#CFLAGS=-Wall -O4 

# On CRAY T3E
# CC=cc

# On Linux NEC-cluster with Intel-MPI
#CC=mpiicc
#LFLAGS=-lm -lcseife -lfftw3 -ftz -lstdc++
#CFLAGS=-O3 -xAVX -ipo -fno-fnalias -restrict
#SFLAGS=-L./../libcseife -L$//sfs/fs2/work-sh2/sungw331/fftw-3.3.4/lib
#IFLAGS=-I./../libcseife -I$//sfs/fs2/work-sh2/sungw331/fftw-3.3.4/include -I./../include

# On Desktop computer with LinuxMint 17, OpenMPI and gcc 4.8.2
#  CC=mpiicc  #tesla:cc60,managed
#  LFLAGS=-lm -lcseife -lfftw3 -ftz -lstdc++
#  LFLAGS=-lm -ftz -lstdc++
#  CFLAGS=-O3 -xAVX -iop -fno-fnalias -restrict
#  SFLAGS=-L./../libcseife
#  IFLAGS=-I./../libcseife -I./../include
#  IFLAGS= -I./../include

# On Desktop computer with LinuxMint 17, OpenMPI and gcc 4.8.2
CFLAGS=-w -fno-stack-protector -D_FORTIFY_SOURCE=0
# SFLAGS=-L./../libcseife

IFLAGS=-I../../include
# On HLRN system
#CC=mpcc
#LFLAGS=-lm

# ALTIX
#CC=icc
#CFLAGS=-mp -O3 -ipo
#LFLAGS=-lmpi -lm -i-static

# after this line, no further editing should be necessary
# --------------------------------------------------------

.c.o:
	$(MPICC) $(CFLAGS) -c $< $(IFLAGS)

VPATH = snapmerge:PSV:VTI:TTI:AC
        
SNAPMERGE_SCR = \
	merge.c \
	read_par.c \
	readdsk.c \
	snapmerge.c \
	writedsk.c \
	util.c

SRC := $(wildcard *.c)

# add PSV source files
SRC += alloc_matPSV.c \
	  alloc_mpiPSV.c \
	  alloc_PSV.c \
 	  alloc_seisPSV.c \
	  alloc_seisPSVfull.c \
	  dealloc_PSV.c \
	  exchange_s_PSV.c \
	  exchange_v_PSV.c \
	  FD_PSV.c \
	  matcopy_PSV.c \
	  matcopy_elastic_PSV.c \
	  mem_PSV.c \
	  outseis_PSVfor.c \
	  psv.c \
	  readmod_elastic_PSV.c \
	  surface_elastic_PML_PSV.c \
	  update_s_elastic_PML_PSV.c \
	  update_v_PML_PSV.c \
	  zero_denise_elast_PSV.c

# add VTI source files
SRC += alloc_matVTI.c \
	  checkfd_ssg_VTI.c \
	  FD_VTI.c \
	  matcopy_elastic_VTI.c \
	  readmod_elastic_VTI.c \
	  seismo_ssg_VTI.c \
	  snap_ssg_VTI.c \
	  update_s_elastic_PML_VTI.c \
	  VTI.c

# add TTI source files
SRC += alloc_matTTI.c \
	  checkfd_ssg_TTI.c \
	  FD_TTI.c \
	  readmod_elastic_TTI.c \
	  rot_el_tensor_TTI.c \
	  update_s_elastic_PML_TTI.c \
	  TTI.c

# add AC source files
SRC += ac.c \
	  alloc_AC.c \
	  alloc_matAC.c \
	  checkfd_acoustic.c \
	  dealloc_AC.c \
	  exchange_p_AC.c \
	  exchange_v_AC.c \
	  FD_AC.c \
	  matcopy_acoustic_AC.c \
	  psource_AC.c \
	  readmod_AC.c \
	  seismo_AC.c \
	  snap_AC.c \
	  surface_acoustic_PML_AC.c \
	  update_s_acoustic_PML_AC.c \
	  update_v_PML_AC.c \
	  zero_denise_acoustic_AC.c


SNAPMERGE_OBJ = $(SNAPMERGE_SCR:.c=.o)
SRC_OBJ = $(SRC:.c=.o)      

all: seisM2D snapmerge

libfd.a: $(SRC_OBJ)
	ar rcs libfd.a $(SRC_OBJ)
	ranlib libfd.a

snapmerge2D:	$(SNAPMERGE_OBJ)
	$(MPICC) $(SNAPMERGE_OBJ) -o snapmerge2D $(FLAGS)

forwardFD2D:		$(SRC_OBJ)
	$(MPICC) $(SRC_OBJ) -o forwardFD2D $(FLAGS)
 
install : forwardFD2D snapmerge2D  libfd.a
	cp forwardFD2D snapmerge2D ../../../bin
	cp libfd.a ../../lib

clean:
	find . -name "*.o" -exec rm {} \; 
	find . -name "*.c%" -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	rm forwardFD2D snapmerge2D 

