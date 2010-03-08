#!/bin/tcsh
#
#$ -S /bin/tcsh
#$ -cwd
#$ -pe cordelia 1
#$ -l mem_free=1.5G

set EXEDIR=/data/rw14/galtay/sphray_update/sphray/src
set EXE=${EXEDIR}/sphray
set INDIR=/data/rw14/galtay/sphray_update/sphray/data/config/iliev_small_tests
set INFILE=${INDIR}/iliev_test2_N64_R6.config


${EXE} ${INFILE} 

