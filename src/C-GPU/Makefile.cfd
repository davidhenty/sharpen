# Set the target system
# Can use one of this list: archer, cirrus or linux

SYSTEM=cirrusgpu

# System dependent definitions

# For ARCHER - default

ifeq ($(SYSTEM), archer)
CC=	cc
CFLAGS=	-fopenmp
LFLAGS=
endif

# For Cirrus using Intel compilers

ifeq ($(SYSTEM), cirrus)
CC=	icc
CFLAGS=	-qopenmp
LFLAGS=	-lm
endif

# For Cirrus using GPU offload

ifeq ($(SYSTEM), cirrusgpu)
CC=	nvc
CFLAGS=	-mp=gpu
LFLAGS=	-lm
endif

# For a standard Linux cluster using GNU compilers

ifeq ($(SYSTEM), linux)
CC=	gcc
CFLAGS=	-fopenmp
LFLAGS=	-lm
endif

# System independent definitions

MF=	Makefile

EXE=	cfd

INC= \
	arraymalloc.h \
	boundary.h \
	cfdio.h \
	jacobi.h

SRC= \
	arraymalloc.c \
	boundary.c \
	cfd.c \
	cfdio.c \
	jacobi.c

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(OBJ):	$(INC)

$(EXE):	$(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF)

tar:
	tar cvf cfd.tar $(MF) $(INC) $(SRC)

clean:
	rm -f $(OBJ) $(EXE) velocity.dat colourmap.dat cfd.plt core
