#!/bin/tcsh
#
#$ -S /bin/tcsh

set OUTDIR=/data/rw14/galtay/sphray_update/sphray_output/IT1_He/r7
qsub -o ${OUTDIR} -e ${OUTDIR} it1_he_r7.sh 
