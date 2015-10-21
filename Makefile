# -----------------------------------------------------------------
# Version: 2.0
# Date: "Nov, 2007" 
# -----------------------------------------------------------------
# Programmer: Mukesh Kumar (muk139@psu.edu)@ PSU
# -----------------------------------------------------------------
# Makefile for PIHM 
#
# cvode/pihm/Makefile.  
# -----------------------------------------------------------------

SHELL = /bin/sh

srcdir       = ./src
rtdir        = ./src
builddir     = .
top_builddir = ../../
top_builddir = ../../
prefix       = ../sundials
exec_prefix  = ${prefix}
includedir   = ${prefix}/include
libdir       = ${exec_prefix}/lib



CPP      = /usr/bin/cc -E
CPPFLAGS = 
CC       = /usr/bin/gcc
CFLAGS   = -g
#CFLAGS   = -O2
LDFLAGS  = 
OPENMP   = -fopenmp
LIBS     = -lm
SRC    = ${srcdir}/pihm.c ${srcdir}/f.c ${srcdir}/read_alloc.c ${srcdir}/initialize.c ${srcdir}/seb.c ${srcdir}/update.c ${srcdir}/print.c ${srcdir}/meta.c ${srcdir}/swc.c   ${rtdir}/rt.c  ${rtdir}/os3d.c  ${rtdir}/react.c ${includedir}/smalldense.h

FORMATION    = ${srcdir}/flow.c  ${srcdir}/rt.h ${srcdir}/os3d.c  ${srcdir}/react.c  ${srcdir}/file_rt.h

COMPILER_PREFIX = 
LINKER_PREFIX   = 

SUNDIALS_INC_DIR = $(includedir)
SUNDIALS_LIB_DIR = $(libdir)
SUNDIALS_LIBS    = -lsundials_cvode -lsundials_nvecserial -lsundials_shared

# EXEC_FILES = cvdx cvdxe cvbx cvkx cvkxb cvdemd cvdemk

all:
	@(echo)
	@(echo '       make pihm     - make pihm                  ')
	@(echo '       make fert     - make 2d formation model    ')
	@(echo '       make clean    - remove all executable files')
	@(echo)

pihm:   $(SRC) Makefile 
	@echo '...Compiling PIHM ...'
	@$(CC) $(CFLAGS) $(OPENMP) -I$(SUNDIALS_INC_DIR) -L$(SUNDIALS_LIB_DIR) -o $(builddir)/pihm $(SRC) $(SUNDIALS_LIBS) $(LIBS)


fert:   $(FORMATION) Makefile
	@echo '...Compiling 2D Formation model ...'
	@$(CC) $(CFLAGS) -I$(SUNDIALS_INC_DIR) -L$(SUNDIALS_LIB_DIR) -o $(builddir)/fert $(FORMATION) $(SUNDIALS_LIBS) $(LIBS)

clean:
	@rm -f *.o
	@rm -f pihm

