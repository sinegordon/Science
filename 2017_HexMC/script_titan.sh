#!/bin/bash -login
#PBS -V
#PBS -A BIP151
#PBS -N hexaflake.test
#PBS -j oe
#PBS -l walltime=720:00:00,nodes=1

module unload PrgEnv-pgi
module load PrgEnv-gnu

cd $PBS_O_WORKDIR
CC -fopenmp hexaflake.c -o hexaflake.x
cp ./hexaflake.x $MEMBERWORK/BIP151/hexaflake.x
cd $MEMBERWORK/BIP151/
setenv OMP_NUM_THREADS 8
#Arguments - flake radius, flake iterations, threads count, nx and ny, steps count
aprun -n1 -d16 ./hexaflake.x 0.8 2 8 256 1000

wait #wait for all processes to finish before ending the job.
