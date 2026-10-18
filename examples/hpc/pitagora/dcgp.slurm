#!/usr/bin/env bash
#SBATCH --job-name=simple-cpu
#SBATCH --account=<project>
#SBATCH --partition=dcgp_fua_dbg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --hint=nomultithread
#SBATCH --time=00:30:00
#SBATCH --output=simple-cpu-%j.log
# UNTESTED TEMPLATE. Production: --partition=dcgp_fua_prod --time=24:00:00.
# Submit from a run directory containing simple.in and its public input files.
# Load the exact module stack used to build SIMPLE before sbatch (see docs).
set -euo pipefail
: "${SIMPLE_EXE:?Set SIMPLE_EXE to an absolute path to build/simple.x}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores OMP_PROC_BIND=spread OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
srun --cpu-bind=cores "$SIMPLE_EXE" simple.in
