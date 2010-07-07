FC=ifort
FFLAGS= -traceback -fpe0 -C

APPS = create_spectra_file 

#--------------------

NRDIR=NR
SRCDIR=../../src

NRFILES = $(NRDIR)/nrtype.o $(NRDIR)/nrutil.o $(NRDIR)/nr.o \
          $(NRDIR)/trapzd.o $(NRDIR)/polint.o $(NRDIR)/qromb.o 



FILES = $(NRFILES) $(SRCDIR)/myf90.o $(SRCDIR)/physical_constants.o

INCLUDE= -I$(SRCDIR) 
OPT=
LIBS= 

#=============================================================================

all:$(APPS)

create_spectra_file: $(FILES) create_spectra_file.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) $(INCLUDE)

%.o: %.F90 
	$(FC) $(FFLAGS) $(OPT) $(INCLUDE) $(LIBS) -c -o $@ $<

%.o: %.f90 
	$(FC) $(FFLAGS) $(OPT) $(INCLUDE) $(LIBS) -c -o $@ $<

clean:
	rm -f *.o $(NRDIR)/*.o *.mod $(NRDIR)/*.mod

cleanall:
	rm -f *.o $(NRDIR)/*.o *.mod $(NRDIR)/*.mod $(APPS) 

tidy:
	rm -f *~