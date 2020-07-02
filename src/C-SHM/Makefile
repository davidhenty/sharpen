MF=	Makefile

#COMPILER
CC=	cc
# On ARCHER, cc is a wrapper for whichever C compiler - GNU (gcc), Intel (icc), or Cray (craycc) - has been chosen by loading the appropriate PrgEnv module prior to compilation.

CFLAGS=	-g -DC_OPENSHMEM_PRACTICAL
LFLAGS=	-lm

EXE=	sharpen

SRC= \
	sharpen.c \
	dosharpen.c \
	filter.c \
	cio.c \
	utilities.c

INC = \
	sharpen.h \
	utilities.h

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF) $(INC)

clean:
	rm -f $(OBJ) $(EXE) core
