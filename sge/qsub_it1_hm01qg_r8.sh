#!/bin/tcsh
#
#$ -S /bin/tcsh

set OUTDIR=/data/rw14/galtay/sphray_update/sphray_output/IT1_HM01QG/r8
qsub -o ${OUTDIR} -e ${OUTDIR} it1_hm01qg_r8.sh 
