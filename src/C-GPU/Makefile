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
#CFLAGS=	-mp=gpu -g -DC_OPENMP_PRACTICAL
CFLAGS=	-g -DC_OPENMP_PRACTICAL
LFLAGS=	-lm

EXE=	sharpen

CSRC= \
	sharpen.c \
	filter.c \
	cio.c \
	utilities.c

CUSRC= \
	dosharpen.cu \
	dosharpenpixel.cu

INC = \
	sharpen.h \
	utilities.h

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .cu .o

COBJ=	$(CSRC:.c=.o)
CUOBJ=	$(CUSRC:.cu=.o)

OBJ=$(COBJ) $(CUOBJ)

.c.o:
	$(CC) $(CFLAGS) -c $<

.cu.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF) $(INC)

clean:
	rm -f $(OBJ) $(EXE) core
