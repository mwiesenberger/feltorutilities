#!/bin/bash

#SBATCH -N 4 -n 16 --ntasks-per-node=4 --exclusive
#####SBATCH --gres=gpu:4
#SBATCH --gpus-per-node=4
#SBATCH --partition=boost_fua_prod
#SBATCH --account=fual8_feltor
#SBATCH --time=24:00:00 # 24 hours is maximum

echo "Leonardo cluster with Tesla A100 GPUs"

hostname
date
module list
echo "$@"

: ${FELTOR_PATH:="../../feltor/build/mpi-gpu"}

# $@ forwards all arguments
echo "1 1 16" | mpirun -n 16 $FELTOR_PATH/src/feltor/feltor "$@"

date
