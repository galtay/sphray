#!/bin/tcsh
#
#$ -S /bin/tcsh

set OUTDIR=/data/rw14/galtay/sphray_update/sphray_output/IT2/r6
qsub -o ${OUTDIR} -e ${OUTDIR} it2_r6.sh 
