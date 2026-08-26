# GPU Slurm examples

These examples run one SIMPLE process on one NVIDIA GPU. They assume that the
working directory contains `simple.in` and any referenced field or chartmap
files. Override the defaults when submitting from another location:

```bash
SIMPLE_BUILD=$HOME/code/SIMPLE/build-cuda \
SIMPLE_INPUT=$PWD/simple.in \
sbatch /path/to/SIMPLE/RUN/Perlmutter/submit_gpu.slurm
```

The production GPU path is native CUDA (`gpu_backend = 'cuda-native'`), not
the legacy OpenACC path. Build once before submitting a run. The CUDA
architecture must match the GPU:

| system | GPU | CUDA architecture |
| --- | --- | --- |
| Perlmutter | NVIDIA A100 | `80` |
| aCluster | NVIDIA Tesla T4 | `75` |
| sCluster | NVIDIA RTX PRO 6000 Blackwell Max-Q | `120` |

For Perlmutter, the normal GNU toolchain is sufficient for SIMPLE's Fortran
host code. The CUDA kernels are compiled by `nvcc`. A GNU build can use the
Cray wrappers and the current CPE/CUDA modules:

```bash
SIMPLE_ROOT=${SIMPLE_ROOT:-$HOME/code/SIMPLE}
SIMPLE_BUILD=${SIMPLE_BUILD:-$SIMPLE_ROOT/build-cuda}
module load PrgEnv-gnu cudatoolkit/12.9 craype-accel-nvidia80
module load cray-hdf5 cray-netcdf
cmake -S "$SIMPLE_ROOT" -B "$SIMPLE_BUILD" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=cc \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_Fortran_COMPILER=ftn \
  -DSIMPLE_ENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=80
cmake --build "$SIMPLE_BUILD" -j
```

If you want the NVIDIA Fortran toolchain, replace `PrgEnv-gnu` with
`PrgEnv-nvidia`. Perlmutter's `ftn` wrapper then invokes `nvfortran`. Check the
exact available version with `module spider nvidia` or `module spider nvhpc`.

The local cluster examples use the user-local NVHPC installations currently
available to the SIMPLE account. Set `NVHPC_ROOT` if the installation is in a
different location. `--gres=gpu:1` is the local-cluster request. Perlmutter
uses `--gpus-per-task=1`. For multiple Ax workers on Perlmutter, request one
GPU per task and use an explicit mapping such as
`--gpu-bind=map_gpu:0,1,2,3`.
