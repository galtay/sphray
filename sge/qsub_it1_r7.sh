#!/bin/tcsh
#
#$ -S /bin/tcsh

set OUTDIR=/data/rw14/galtay/sphray_update/sphray_output/IT1/r7
qsub -o ${OUTDIR} -e ${OUTDIR} it1_r7.sh 
