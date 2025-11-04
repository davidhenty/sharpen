MF=	Makefile

#COMPILER

# On ARCHER, cc is a wrapper for whichever C compiler - GNU (gcc),
# Intel (icc), or Cray (craycc) - has been chosen by loading the
# appropriate PrgEnv module prior to compilation.

# For GPUs on ARCHER2 you need these modules:
# module load rocm
# module load craype-accel-amd-gfx90a
# module load craype-x86-milan

CC=	hipcc
CFLAGS= -D__HIP_PLATFORM_AMD__ --offload-arch=gfx90a -x hip -DC_SERIAL_PRACTICAL -g -fgpu-rdc
LFLAGS=$(CFLAGS) --hip-link

#CC=	cc
#CFLAGS= --offload-arch=gfx90a -x hip -DC_SERIAL_PRACTICAL -g -fgpu-rdc
#LFLAGS=$(CFLAGS) --hip-link

#For Cirrus GPU use nvidia
#CC=	nvcc
#CFLAGS=	-arch=sm_70 -g -DC_SERIAL_PRACTICAL -dc
#LFLAGS=	-arch=sm_70 -g -lm

EXE=	sharpen

SRC= \
	dosharpen.c \
	dosharpenpixel.c \
	filter.c \
	sharpen.c \
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
	$(CC) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF) $(INC)

clean:
	rm -f $(OBJ) $(EXE) core
