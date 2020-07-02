MF=	Makefile

#COMPILER

# On ARCHER, ftn is a wrapper for whichever Fortran compiler - GNU
# (gfortran), Intel (ifort), or Cray (crayftn) - has been chosen by
# loading the appropriate PrgEnv module prior to compilation.

#FC=	ftn
#CC=	cc
#FFLAGS=	-g -DF_MPI_PRACTICAL
#CFLAGS= $(FFLAGS)
#LFLAGS=

#For Cirrus use mpif90 and mpicc
FC=	mpif90
CC=	mpicc
FFLAGS=	-g -DF_MPI_PRACTICAL
CFLAGS=	-cc=icc $(FFLAGS)
LFLAGS=

EXE=	sharpen

FSRC= \
	sharpen.f90	\
	dosharpen.f90 	\
	filter.f90	\
	fio.f90

CSRC=	utilities.c

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .f90 .c .o

OBJ=	$(FSRC:.f90=.o) $(CSRC:.c=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC) fuzzy.pgm submit.sh

clean:
	rm -f $(OBJ) $(EXE) core
