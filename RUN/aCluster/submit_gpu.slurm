#!/bin/bash
# One SIMPLE process on one aCluster Tesla T4.
# Submit from a run directory containing simple.in and its input files.
#
#SBATCH --job-name=simple-cuda
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=job.out.%j
#SBATCH --error=job.err.%j

set -euo pipefail

# Optional user-local NVHPC runtime, matching the current aCluster install.
NVHPC_ROOT="${NVHPC_ROOT:-$HOME/opt/nvidia/hpc_sdk/Linux_x86_64/26.3}"
if test -d "$NVHPC_ROOT/compilers/bin"; then
    export PATH="$NVHPC_ROOT/compilers/bin:$PATH"
    export LD_LIBRARY_PATH="$NVHPC_ROOT/compilers/lib:${LD_LIBRARY_PATH:-}"
fi

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
srun --cpu-bind=cores "$SIMPLE_BUILD/simple.x" "$SIMPLE_INPUT"
