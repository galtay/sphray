#!/bin/tcsh
#$ -S /bin/tcsh
#$ -pe cordelia 64
#$ -cwd
#$ -l mem_free=1.0G
#$ -t 1:10

set EXEDIR=/data/rw14/galtay/idle-hands/codes/sph-to-grid/src
set EXE=${EXEDIR}/grid_example_mpi
set INDIR=/data/rw14/galtay/idle-hands/codes/sph-to-grid/infiles/cosma

set INFILES= ("" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "" "")


set INFILES[1]=${INDIR}/gimic/low_Sigmam2_023_z_H1_eosall_g16384.in
set INFILES[2]=${INDIR}/gimic/low_Sigmam1_023_z_H1_eosall_g16384.in
set INFILES[3]=${INDIR}/gimic/low_Sigma0_023_z_H1_eosall_g16384.in
set INFILES[4]=${INDIR}/gimic/low_Sigmap1_023_z_H1_eosall_g16384.in
set INFILES[5]=${INDIR}/gimic/low_Sigmap2_023_z_H1_eosall_g16384.in

set INFILES[6]=${INDIR}/gimic/low_Sigmam2_023_z_H1_eos1.0_g16384.in
set INFILES[7]=${INDIR}/gimic/low_Sigmam1_023_z_H1_eos1.0_g16384.in
set INFILES[8]=${INDIR}/gimic/low_Sigma0_023_z_H1_eos1.0_g16384.in
set INFILES[9]=${INDIR}/gimic/low_Sigmap1_023_z_H1_eos1.0_g16384.in
set INFILES[10]=${INDIR}/gimic/low_Sigmap2_023_z_H1_eos1.0_g16384.in


# set INFILES[11]=${INDIR}/gimic/low_Sigmam2_023_z_H1_eos0.0_g16384.in
# set INFILES[12]=${INDIR}/gimic/low_Sigmam1_023_z_H1_eos0.0_g16384.in
# set INFILES[13]=${INDIR}/gimic/low_Sigma0_023_z_H1_eos0.0_g16384.in
# set INFILES[14]=${INDIR}/gimic/low_Sigmap1_023_z_H1_eos0.0_g16384.in
# set INFILES[15]=${INDIR}/gimic/low_Sigmap2_023_z_H1_eos0.0_g16384.in



echo $INFILES[${SGE_TASK_ID}]

${MPIRUN} -np ${NSLOTS} ${CORDELIANGE} ${EXE} $INFILES[${SGE_TASK_ID}]






