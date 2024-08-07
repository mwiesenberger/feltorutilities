#!/bin/bash
# submit (several times) with 
# sbatch  -o diag.out -J diag job_diag.sh

#SBATCH -N 1 
#SBATCH --partition=m100_all_serial
#SBATCH --account=fuac7_tsvv3
#SBATCH --time=4:00:00 # 24 hours is maximum

hostname
date
module load python
module list

python run_diag.py

date
