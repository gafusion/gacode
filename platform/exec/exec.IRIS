#!/bin/sh
# GACODE Parallel execution script (MINT)

simdir=${1}
nmpi=${2}
exec=${3}
nomp=${4}
numa=${5}
mpinuma=${6}

cd $simdir
mpiexec -env OMP_NUM_THREADS $nomp -n $nmpi $exec

