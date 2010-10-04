all:
	@\cd src && $(MAKE)
	@\cd test && $(MAKE)

clean:
	@\cd src && $(MAKE) clean
	\rm -f hdf5_wrapper.tar.gz

dist:
	@\mkdir hdf5_wrapper
	@\cp Makefile hdf5_wrapper
	@\mkdir hdf5_wrapper/src
	@\cp src/*Makefile hdf5_wrapper/src
	@\cp src/config* hdf5_wrapper/src
	@\cp src/*.inc hdf5_wrapper/src
	@\cp src/*.F90 hdf5_wrapper/src
	@\cp src/*.c   hdf5_wrapper/src
	@\mkdir hdf5_wrapper/test
	@\cp test/*.f90 hdf5_wrapper/test
	@\tar cvf hdf5_wrapper.tar hdf5_wrapper > /dev/null
	@\gzip hdf5_wrapper.tar
	@\rm -rf hdf5_wrapper
