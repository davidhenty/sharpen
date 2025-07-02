MF=	Makefile

#COMPILER

# On ARCHER, cc is a wrapper for whichever C compiler - GNU (gcc),
# Intel (icc), or Cray (craycc) - has been chosen by loading the
# appropriate PrgEnv module prior to compilation.

#CC=	cc
#CFLAGS=	-g -DC_OPENMP_PRACTICAL
#LFLAGS=	-lm

#For Cirrus GPU use nvidia
CC=	nvcc
CFLAGS=	-arch=sm_70 -g -DC_SERIAL_PRACTICAL -dc
LFLAGS=	-arch=sm_70 -g -lm

EXE=	sharpen

SRC= \
	dosharpen.cu \
	dosharpenpixel.cu \
	filter.cu \
	sharpen.cu \
	cio.cu \
	utilities.cu


INC = \
	sharpen.h \
	utilities.h

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .cu .o

OBJ=	$(SRC:.cu=.o)

.cu.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF) $(INC)

clean:
	rm -f $(OBJ) $(EXE) core
