#!/bin/tcsh
#
#$ -S /bin/tcsh

set OUTDIR=/data/rw14/galtay/owls-gimic-dlas/output/owls/REF_L025N512/r6
qsub -o ${OUTDIR} -e ${OUTDIR} ref_l025n512.sh 
