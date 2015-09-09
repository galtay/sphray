# SPHRAY #

SPHRAY is a program for ray-tracing smoothed particle hydrodynamics (SPH) distributions to solve the radiative transfer problem.  Currently, the radiative transfer is not coupled to the hydrodynamics and so SPHRAY can be termed a post-processing radiative transfer code.  It accomplishes this goal using monte carlo ray tracing.  I've included a quick start guide on this page to help users get started.  In the external links section I've included links to public code I have made use of in SPHRAY.

**WARNING**
The results of this code are by no means guaranteed.  I have performed standard cosmological radiative transfer tests using this code, however the results depend strongly on the particulars of the density field, the number of sources used, and the number of rays traced.  I am happy to consult with users about projects, but please exercise caution in the interpretation of results and perform convergence studies.

## Optional Libraries ##

### HDF5 ###
If you would like to use the HDF5 I/O features of SPHRAY you will need a working installation of the HDF5 library.  In addition, it may be necessary to compile the HDF5 library, the fortran wrapper included with SPHRAY, and SPHRAY itself with the same Fortran compiler.

An example of how to do this with the intel fortran compiler ifort is given below.  Replace ifort with the compiler you will be using (gfortran for example).
  * Download the source for the latest official release of HDF5 (1.8.5-patch1 as of this writing) http://www.hdfgroup.org/HDF5/release/obtain5.html
  * untar the file creating the source directory hdf5-1.8.5-patch1
```
tar xvf hdf5-1.8.5-patch1.tar
```
  * rename the directory hdf5-1.8.5-patch1-source
```
mv hdf5-1.8.5-patch1 hdf5-1.8.5-patch1-source
```
  * create the build directory and move into it
```
mkdir hdf5-1.8.5-patch1-ifort-11.1
cd hdf5-1.8.5-patch1-ifort-11.1
```
  * set the environment variable FC to ifort (or whichever Fortran compiler you are using)
```
export FC=ifort
```
  * add the intel compiler install directory to your library search path.  I didn't have to do this for other compilers.
```
export LD_LIBRARY_PATH=/opt/intel/Compiler/11.1/073/lib/intel64:$LD_LIBRARY_PATH
```
  * run the configure script in the source directory with the enable fortran option
```
../hdf5-1.8.5-patch1-source/configure --enable-fortran
```
  * now run several makes and checks
```
make
make check
make install
make check-install
```

The end result should be a copy of the HDF5 library installed in the build directory.


## Getting Started ##
The quickest way to get a test problem running is to follow these step by step instructions.

  * Create a working directory where the code and the RTCP snapshots can be stored.  For example,
```
 mkdir sphray_download 
```

  * Make the working directory your current directory and enter the svn commands below to download the SPHRAY code and the snapshots to this working directory.
```
  cd sphray_download
  svn checkout http://sphray.googlecode.com/svn/trunk/ sphray
  svn checkout http://rt-comparison-project-snapshots.googlecode.com/svn/trunk/ rtcp_snapshots
```

  * Create directories for the output of the RTCP tests.
```
  cd sphray
  chmod u+x make_output_dirs.sh
  ./make_output_dirs.sh 
```
  * Make the sphray source directory your current directory, copy the Makefile template into a Makefile and compile the source.
```
  cd src
  cp MakefileTemplate Makefile
  make
```
  * Note that if you are using the intel compiler, you may have to add directories to your runtime search path.  Examples of how to do this can be found in the file sphray/src/makes/make.ferrari.serial.opt.
  * Run the first or second test using the configuration files in the data directory.
```
  ./sphray ../data/config/iliev_small_tests/iliev_test1_N64_R6.config
  ./sphray ../data/config/iliev_small_tests/iliev_test2_N64_R6.config
```
  * Analyze the results using the idl scripts in the analysis directory
```
  sphray_download/sphray/analysis/plot_IT1.pro
  sphray_download/sphray/analysis/plot_IT2x.pro
  sphray_download/sphray/analysis/plot_IT2t.pro
```
  * In addition, I have arranged a config file that uses SPHRAY's ability to solve for the Helium ionization fraction.  This file uses the same density field as the tests above, keeps the temperature fixed at 10,000 K like the first test and uses the thermal black body source from the second test.  The amount of Hydrogen is the same as in both tests, but on top of that I've added enough Helium to make the Helium mass fraction = 0.286.  I have also included an IDL analysis script for this test.
```
  ./sphray ../data/config/iliev_small_tests/iliev_test1_He_N64_R6.config
  sphray_download/sphray/analysis/plot_IT1_He.pro
```