#=============================================================================
# these preprocessor macros determine the data structure of a particle
# and can be used to determine what sections of the code need not be 
# compiled.  The default (and simplest) case is to have all of them
# commented out.  
#=============================================================================

OPT += -DincHe    # if you want to include Helium
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
#SYSTEM=cdm
SYSTEM=profile
#SYSTEM=debug



#=============================================================================
#=============================================================================

IFORTFLAGS= -ftz -fpe0 -g -traceback -heap-arrays -warn all -fpp 

ifeq ($(SYSTEM),cdm)
   FC=f95
   FFLAGS= 
endif


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


ifeq ($(SYSTEM),generic_g95_opt)
   FC= g95
   FFLAGS= -O3 -ffast-math -cpp
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
       config.o atomic_rates.o \
       source_input.o gadget_input.o main_input.o \
       ionpar.o euler.o bdf.o \
       iliev_comparison_project.o output.o \
       ion_temperature_update.o \
       initialize.o mainloop.o 


#=============================================================================
# these only need to be set if you are using the OpenGL visiualization tool
# Note that you need to compile f90gl with the same Fortran 90 compiler you 
# use to compile SPHRAY.  A make will produce the executable sphray while
# make glsphray will produce the executable glsphray
#=============================================================================

#  Comment out for 32-bit architecture
BITS      := 64

#  The directories containing the OpenGL libraries and GLUT libraries
OGLLIBDIR := -L/usr/lib${BITS} -L/gal/galtay1/galtay/usr/f03gl

#  The directory containing the X11 libraries
X11LIBDIR := -L/usr/X11R6/lib${BITS}

#  Libraries for OpenGL, including GLUT, GLU and OpenGL
F90GLUTLIB := -lglut -lGL -lGLU

#  The X11 libraries
X11LIB := -lXaw -lXt -lXmu -lXi -lXext -lX11 -lm -lXxf86vm





#  If you use OpenGLUT or FreeGlut, change this variable and adjust
#  the libraries appropriately.
 GLUT      := glut
# GLUT      := openglut
# GLUT      := freeglut

LIBRARIES := ${OGLLIBDIR} ${X11LIBDIR} -L/usr/local/lib
GLLIBS    := ${F90GLUTLIB} ${X11LIB} \
             -lg2c -lpthread -ldl
MODULES    = GLUT_fonts.o OpenGL_${GLUT}.o OpenGL_glu.o OpenGL_gl.o 


GLINCLUDE = -I/home/galtay/f03gl

THREAD_OPTIONS = -L. -lpthread -lpt

#=============================================================================

all:$(APPS)

sphray: $(FILES) sphray.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)


%.o: %.F90 
	$(FC) $(FFLAGS) $(OPT) $(INCLUDE) -c -o $@ $<

%.o: %.f90 
	$(FC) $(FFLAGS) $(OPT) $(INCLUDE) -c -o $@ $<


density_test: $(FILES) density_test.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)



clean :
	rm -f libpt.a *.o *.mod *__genmod.f90

cleanall :
	rm -f libpt.a *.o *.mod density_test glsphray $(APPS) 

tidy :
	rm -f *~ 

libpt.a: pt.c ptf77.c pt.h
	cc -c pt.c
	cc -c ptf77.c
	ar crv libpt.a pt.o ptf77.o
	ranlib libpt.a


GLUT_fonts.o: GLUT_fonts.c
	cc -c -o GLUT_fonts.o GLUT_fonts.c

OpenGL_${GLUT}.o OpenGL_glu.o: OpenGL_gl.o



glsphray: %: %.F90 $(MODULES) $(FILES) viewer.o fthread.o libpt.a force
	${FC} ${FFLAGS} -c $<
	${FC} $@.o ${FILES} fthread.o viewer.o libpt.a \
	$(THREAD_OPTIONS) ${MODULES} ${GLINCLUDE} ${LIBRARIES} ${GLLIBS} -o $@

force:
