#!/bin/bash
# One SIMPLE process on one Perlmutter A100.
# Submit from a run directory containing simple.in and its input files.
#
#SBATCH --job-name=simple-cuda
#SBATCH --account=<account>
#SBATCH --constraint=gpu
#SBATCH --qos=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=job.out.%j
#SBATCH --error=job.err.%j

set -euo pipefail

module load PrgEnv-gnu cudatoolkit/12.9 craype-accel-nvidia80
module load cray-hdf5 cray-netcdf

RUN_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
SIMPLE_BUILD="${SIMPLE_BUILD:-$HOME/code/SIMPLE/build-cuda}"
SIMPLE_INPUT="${SIMPLE_INPUT:-$RUN_DIR/simple.in}"

test -x "$SIMPLE_BUILD/simple.x" || {
    echo "missing executable: $SIMPLE_BUILD/simple.x" >&2
    exit 1
}
test -r "$SIMPLE_INPUT" || {
    echo "missing input: $SIMPLE_INPUT" >&2
    exit 1
}

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

nvidia-smi --query-gpu=name,compute_cap,memory.total --format=csv,noheader

# --gpus-per-task=1 supplies one GPU and Slurm's per-task GPU binding.
srun --cpu-bind=cores "$SIMPLE_BUILD/simple.x" "$SIMPLE_INPUT"
