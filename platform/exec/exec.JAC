#!/bin/sh
#
# SCRIPT:
#  gyro.DELPHI
#
# FUNCTION:
#  Parallel execution script

simdir=${1}
nmpi=${2}
exec=${3}
nomp=${4}
numa=${5}
mpinuma=${6}

cd $simdir
export OMP_NUM_THREADS=$nomp ; mpiexec -n $nmpi $exec

