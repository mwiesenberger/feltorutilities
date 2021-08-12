#!/bin/bash

#SBATCH -N 4 -n 16 --ntasks-per-node=4 --exclusive
#SBATCH --gres=gpu:4
#SBATCH --partition=m100_fua_prod
#SBATCH --account=fuac5_tsvv3
#SBATCH --time=24:00:00 # 24 hours is maximum

echo "Marconi 100 cluster with Tesla V100 GPUs"

hostname
date
module list
echo "$@"

: ${FELTOR_PATH:="../feltor"}

# $@ forwards all arguments
echo "1 1 16" | mpirun -gpu -n 16 $FELTOR_PATH/src/feltor/feltor_mpi "$@"

date
