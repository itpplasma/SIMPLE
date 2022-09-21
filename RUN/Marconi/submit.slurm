#!/bin/bash
# Standard output and error:
#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J simple
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-core=2
#

export OMP_NUM_THREADS=96
export OMP_PLACES=threads

srun $HOME/src/SIMPLE/BUILD/simple.x
