#!/bin/bash -l
# Derived template; not yet run on Pitagora. Replace PROJECT_ACCOUNT first.
# Export SIMPLE_ENV as an absolute path to your validated Pitagora module script.
#SBATCH --job-name=simple-cpu
#SBATCH --account=PROJECT_ACCOUNT
#SBATCH --partition=dcgp_fua_dbg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:03:00
#SBATCH --output=cineca-cpu-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the SIMPLE root}"
# shellcheck source=/dev/null
source "${SIMPLE_ENV:?Set an absolute path to the validated Pitagora environment}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_THREAD_LIMIT="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS=1 OMP_PLACES=cores OMP_PROC_BIND=close
test -x build/simple.x
test -s run/cineca-smoke/wout.nc
cd run/cineca-smoke
srun --cpu-bind=cores ../../build/simple.x simple.in
