OPT += -Dhdf5

include ../../src/makes/make.santabuntu.serial.debug



APPS= write_rates_Cen write_rates_Hui write_rates_Hum \
      write_rates test_rates_vs_iontables

#--------------------

SRCDIR=../../src

FILES= $(SRCDIR)/myf90.o $(SRCDIR)/hui_gnedin_atomic_rates.o \
       $(SRCDIR)/cen_atomic_rates.o $(SRCDIR)/hummer_atomic_rates.o \
       $(SRCDIR)/voronov_atomic_rates.o $(SRCDIR)/atomic_rates.o \
       $(SRCDIR)/ion_table_class.o

INCLUDE = -I$(SRCDIR)

ifeq (hdf5,$(findstring hdf5, $(OPT)))
   INCHDF= $(ISINC)$(DIRHDF)/include $(ISINC)$(DIRHDF)/lib
   MODHDF= $(ISMOD)$(DIRHDF)/lib
   RLIBHDF= $(ISRLIB)$(DIRHDF)/lib
   LIBHDF= $(ISLIB)$(DIRHDF)/lib 
   OPTHDF= $(INCHDF) $(MODHDF) $(RLIBHDF) $(LIBHDF)  -lhdf5 -lhdfwrapper
endif

OPT += $(OPTHDF)

#=============================================================================

all:$(APPS)

write_rates: $(FILES) write_rates.o
	$(FC) $(FFLAGS) -o $@ $^ $(OPT) $(INCLUDE)

write_rates_Vor: $(FILES) write_rates_Vor.o
	$(FC) $(FFLAGS) -o $@ $^ $(OPT) $(INCLUDE)

write_rates_Cen: $(FILES) write_rates_Cen.o
	$(FC) $(FFLAGS) -o $@ $^ $(OPT) $(INCLUDE)

write_rates_Hui: $(FILES) write_rates_Hui.o
	$(FC) $(FFLAGS) -o $@ $^ $(OPT) $(INCLUDE)

write_rates_Hum: $(FILES) write_rates_Hum.o
	$(FC) $(FFLAGS) -o $@ $^ $(OPT) $(INCLUDE)

test_rates_vs_iontables: $(FILES) test_rates_vs_iontables.o
	$(FC) $(FFLAGS) -o $@ $^ $(OPT) $(INCLUDE)


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
