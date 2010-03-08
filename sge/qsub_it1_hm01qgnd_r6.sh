#!/bin/tcsh
#
#$ -S /bin/tcsh

set OUTDIR=/data/rw14/galtay/sphray_update/sphray_output/IT1_HM01QGnd/r6
qsub -o ${OUTDIR} -e ${OUTDIR} it1_hm01qgnd_r6.sh 
