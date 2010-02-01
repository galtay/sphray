FC=ifort
FFLAGS= -O2 -ip -heap-arrays -traceback -fpe0 -C

APPS= write_rates_Cen write_rates_Hui write_rates_Hum write_rates_Vor write_rates

#--------------------

SRCDIR=../../src

FILES= $(SRCDIR)/myf90.o $(SRCDIR)/hui_gnedin_atomic_rates.o \
       $(SRCDIR)/cen_atomic_rates.o $(SRCDIR)/hummer_atomic_rates.o \
       $(SRCDIR)/voronov_atomic_rates.o 

INCLUDE = -I$(SRCDIR)
OPT=
LIBS=

#=============================================================================

all:$(APPS)

write_rates: $(FILES) write_rates.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) $(INCLUDE)

write_rates_Vor: $(FILES) write_rates_Vor.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) $(INCLUDE)

write_rates_Cen: $(FILES) write_rates_Cen.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) $(INCLUDE)

write_rates_Hui: $(FILES) write_rates_Hui.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) $(INCLUDE)

write_rates_Hum: $(FILES) write_rates_Hum.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS) $(INCLUDE)


%.o: %.F90 
	$(FC) $(FFLAGS) $(OPT) $(INCLUDE) -c -o $@ $<

%.o: %.f90 
	$(FC) $(FFLAGS) $(OPT) $(INCLUDE) -c -o $@ $<

clean:
	rm -f *.o *.mod $(APPS) $(SRCDIR)/*.o

cleanall:
	rm -f *.o *.mod $(APPS) $(SRCDIR)/*.o

tidy:
	rm -f *~
