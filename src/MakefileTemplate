#=============================================================================
# these preprocessor macros determine the data structure of a particle
# and can be used to determine what sections of the code need not be 
# compiled.  The default (and simplest) case is to have all of them
# commented out.  
#=============================================================================

#OPT += -DincHe    # if you want to include Helium
#OPT += -DincHrec  # if you want seperate rays traced for H recomb. radiation
#OPT += -DincHerec # if you want seperate rays traced for He recomb. radiation

OPT += -DoutGamma # if you want to output HI photoionization rate

#=============================================================================
# the following system definitions can be usefull in setting 
# the variables FC (Fortran Compiler) and FFLAGS (Compiler Flags)
# The first two generic definitions should be fine for most systems
# If you want architecture dependent compilation with the intel compiler
# use ifort --help | less to read about the -x? flags.
#=============================================================================

#SYSTEM=generic_g95_opt
#SYSTEM=generic_gfortran_opt
#SYSTEM=generic_ift_opt
#SYSTEM=eddington_ift_opt
#SYSTEM=catbus_ift_opt
#SYSTEM=silky_ift_opt
#SYSTEM=toto_ift_opt
#SYSTEM=profile
SYSTEM=debug

#=============================================================================
# these only need to be set if you are using the OpenGL visiualization tool
# Note that you need to compile f90gl with the same Fortran 90 compiler you 
# use to compile SPHRAY.  A make will produce the executable sphray while
# make glsphray will produce the executable glsphray
#=============================================================================

F90GLDIR=/home/galtay/usr/f90gl-1.2.11
GLUTDIR=/home/galtay/usr/glut-3.7.1

#=============================================================================
#=============================================================================

IFORTFLAGS= -ftz -fpe0 -g -traceback -heap-arrays -warn all -fpp 

ifeq ($(SYSTEM),debug)
   FC= ifort
   FFLAGS= $(IFORTFLAGS) 
endif

ifeq ($(SYSTEM),profile)
   FC= ifort
   FFLAGS= $(IFORTFLAGS) -p
endif

ifeq ($(SYSTEM),generic_ift_opt)
   FC= ifort
   FFLAGS= $(IFORTFLAGS) -O3 -ipo 
endif


ifeq ($(SYSTEM),catbus_ift_opt)
   FC= ifort
   FFLAGS= $(IFORTFLAGS) -fast 
endif


ifeq ($(SYSTEM),eddington_ift_opt)
   FC= ifort
   FFLAGS= $(IFORTFLAGS) -O3 -xP -ipo  
endif


ifeq ($(SYSTEM),toto_ift_opt)
   FC= ifort
   FFLAGS= $(IFORTFLAGS) -O3 -xN
endif


ifeq ($(SYSTEM),psycho_ift_opt)
   FC= ifort             
   FFLAGS= $(IFORTFLAGS) -O3 -xN -static
endif                                                                

ifeq ($(SYSTEM),silky_ift_opt)
   FC= ifort
   FFLAGS= $(IFORTFLAGS) -O3  
endif

ifeq ($(SYSTEM),generic_g95_opt)
   FC= g95
   FFLAGS= -O3 -ffast-math -ip -cpp
endif

ifeq ($(SYSTEM),generic_gfortran_opt)
   FC= gfortran
   FFLAGS= -O3 -Wall -frecord-marker=4 -g -frange-check
endif


APPS= sphray 

#--------------------

FILES= myf90.o physical_constants.o mt19937.o m_mrgrnk.o \
       pluecker.o cosmology.o particle_system.o b2cd.o \
       HuiGnedinAtomicRates.o CenAtomicRates.o \
       HummerAtomicRates.o VoronovAtomicRates.o \
       spectra.o sphpar.o octtree3.o ray.o raylist.o global.o \
       atomic_rates.o \
       source_input.o gadget_input.o main_input.o \
       ionpar.o euler.o bdf.o \
       iliev_comparison_project.o output.o \
       ion_temperature_update.o \
       initialize.o mainloop.o 

INCLUDE = 

LIBS = 

GLINCLUDE = -I/usr/include/GL -I$(GLUTDIR)/include \
 -I$(F90GLDIR)/include/GL

GLLIB = -L$(F90GLDIR)/lib -lf90GLU -lf90GL -lf90glut  -lGLU -lGL \
$(GLUTDIR)/lib/glut/libglut.a

X11LIB = -L/usr/X11R6/lib -lXaw -lXt -lXmu -lXi -lXext -lX11

THREAD_OPTIONS = -L. -lpthread -lpt

#=============================================================================

all:$(APPS)

sphray: $(FILES) sphray.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

glsphray: $(FILES) glsphray.o fthread.o viewer.o libpt.a
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) $(GLLIB) $(X11_OPTIONS) $(THREAD_OPTIONS) $(X11LIB)

%.o: %.f90 Makefile
	$(FC) $(FFLAGS) $(OPT) $(INCLUDE) -c -o $@ $<

density_test: $(FILES) density_test.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

viewer.o: viewer.f90
	$(FC) -cpp $(OPT)  $(GLINCLUDE) -c -o $@ $<

clean :
	rm -f libpt.a *.o *.mod 

cleanall :
	rm -f libpt.a *.o *.mod density_test glsphray $(APPS) 

tidy :
	rm -f *~ 

libpt.a: pt.c ptf77.c pt.h
	cc -c pt.c
	cc -c ptf77.c
	ar crv libpt.a pt.o ptf77.o
	ranlib libpt.a
