#!/bin/bash
# Standard output and error:
#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J simple
#SBATCH --ntasks-per-node=1        
#SBATCH --ntasks-per-core=2
#SBATCH --cpus-per-task=144
#SBATCH --mail-user=calbert@ipp.mpg.de
#
# Wall clock Limit (max. is 24 hours):

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=threads

srun ../../BUILD/simple_mpi.x
