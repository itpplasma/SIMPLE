#!/bin/bash -l
# Submit from the SIMPLE root: sbatch examples/hpc/cineca/leonardo-serial.sbatch
# lrd_all_serial is budget-free: deliberately no account or GPU request.
#SBATCH --job-name=simple-smoke
#SBATCH --partition=lrd_all_serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:03:00
#SBATCH --output=cineca-serial-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the SIMPLE root}"
source examples/hpc/cineca/leonardo-gcc-env.sh
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores OMP_PROC_BIND=close
test -x build/simple.x
test -s run/cineca-smoke/wout.nc
cd run/cineca-smoke
srun --cpu-bind=cores ../../build/simple.x simple.in
