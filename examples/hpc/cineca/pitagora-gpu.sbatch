#!/bin/bash -l
# Derived template; not yet run on Pitagora. Replace PROJECT_ACCOUNT first.
# SIMPLE_ENV must name a validated module script; prepare a CUDA Boozer input.
#SBATCH --job-name=simple-gpu
#SBATCH --account=PROJECT_ACCOUNT
#SBATCH --partition=boost_fua_dbg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --time=00:03:00
#SBATCH --output=cineca-gpu-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the SIMPLE root}"
# shellcheck source=/dev/null
source "${SIMPLE_ENV:?Set an absolute path to the validated Pitagora environment}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OMP_THREAD_LIMIT="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS=1 OMP_PLACES=cores OMP_PROC_BIND=close
test -x build-cuda/simple.x
test -s run/cineca-gpu/wout.nc
cd run/cineca-gpu
srun --cpu-bind=cores ../../build-cuda/simple.x simple.in
