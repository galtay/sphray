#!/bin/tcsh
#
#$ -S /bin/tcsh
#$ -cwd
#$ -pe cordelia 1
#$ -l mem_free=1.5G

set EXEDIR=/data/rw14/galtay/sphray_update/sphray/src
set EXE=${EXEDIR}/sphray
set INDIR=/data/rw14/galtay/owls-gimic-dlas/sphray_config
set INFILE=${INDIR}/REF_L025N512_022.config


${EXE} ${INFILE} 

