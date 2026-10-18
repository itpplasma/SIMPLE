#!/usr/bin/env bash
#SBATCH --job-name=simple-gpu
#SBATCH --account=<project>
#SBATCH --partition=boost_fua_dbg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --time=00:30:00
#SBATCH --output=simple-gpu-%j.log
# UNTESTED TEMPLATE. Production: --partition=boost_fua_prod --time=24:00:00.
# One process uses one H100; the node has four. Multi-device tracing is experimental.
# Submit from a run directory containing simple.in and its public input files.
# Load the exact module stack used to build SIMPLE before sbatch (see docs).
set -euo pipefail
: "${SIMPLE_EXE:?Set SIMPLE_EXE to an absolute path to build-cuda/simple.x}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores OMP_PROC_BIND=spread OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# Keep Slurm's CUDA_VISIBLE_DEVICES mapping; never replace it with a physical ID.
srun --cpu-bind=cores --gpu-bind=single:1 "$SIMPLE_EXE" simple.in
