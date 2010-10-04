
# Makefile for HDF5 wrapper library
#

# Configuration file with compiler names and flags and location of HDF5

CONFIG_FILE = ./config

# -----------------------------------------------------------------------

include $(CONFIG_FILE)

LIBDIR = `pwd`/../lib

CHDF_FLAGS = $(INCFLAG)$(HDF5)/include
FHDF_FLAGS = $(LIBFLAG)$(HDF5)/lib $(RFLAG)$(HDF5)/lib $(MODFLAG)$(HDF5)/lib $(MODFLAG)$(HDF5)/include -lhdf5 -lhdf5_fortran

ALLCFLAGS = $(CHDF_FLAGS) $(CFLAGS)
ALLFFLAGS = $(FHDF_FLAGS) $(FFLAGS)

.SUFFIXES:
.SUFFIXES: .f .c .F90 .f90 .o

.f.o:
	$(FC) $(ALLFFLAGS) -c $< -o $*.o
.f90.o:
	$(FC) $(ALLFFLAGS) -c $< -o $*.o
.c.o:
	$(CC) $(ALLCFLAGS) -c $< -o $*.o
.F90.o:
	$(FC) $(ALLFFLAGS) -c $< -o $*.o

all:	libhdfwrapper.so libhdfwrapper.a Test-Makefile

libhdfwrapper.so:	hdf5_wrapper.F90 read_dataset_c.o write_dataset_c.o read_attribute_c.o write_attribute_c.o byte_order.o
	$(FC) $(SHARED) hdf5_wrapper.F90 read_dataset_c.o write_dataset_c.o read_attribute_c.o write_attribute_c.o byte_order.o $(ALLFFLAGS) -o libhdfwrapper.so
	mkdir -p ../lib
	mv libhdfwrapper.so ../lib/	
	cp hdf5_wrapper.mod ../lib/


libhdfwrapper.a:	hdf5_wrapper.F90 hdf5_wrapper.o read_dataset_c.o write_dataset_c.o read_attribute_c.o write_attribute_c.o byte_order.o
	ar -r libhdfwrapper.a hdf5_wrapper.o read_dataset_c.o write_dataset_c.o read_attribute_c.o write_attribute_c.o byte_order.o
	mkdir -p ../lib
	mv libhdfwrapper.a ../lib/	
	cp hdf5_wrapper.mod ../lib/

Test-Makefile:	$(CONFIG_FILE) Makefile libhdfwrapper.so
	rm -f ../test/Makefile
	echo "FC = $(FC)" >> ../test/Makefile 
	echo "CC = $(CC)" >> ../test/Makefile 
	echo "FFLAGS = $(FFLAGS)" >> ../test/Makefile 
	echo "CFLAGS = $(CFLAGS)" >> ../test/Makefile 
	echo "LINKFLAGS = $(MODFLAG)$(HDF5)/lib $(MODFLAG)$(HDF5)/include/ $(INCFLAG)$(LIBDIR) $(LIBFLAG)$(LIBDIR) $(RFLAG)$(LIBDIR) $(MODFLAG)$(LIBDIR) -lhdfwrapper" >> ../test/Makefile
	cat Test-Makefile >> ../test/Makefile 

clean:
	rm -f *.o *.mod
	rm -f ../test/Makefile
	rm -f ../test/*.o ../test/*.exe ../test/test.hdf5
	rm -rf ../lib
	rm -f *~
