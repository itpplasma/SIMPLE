# Job scripts

Usage, CUDA architectures and the CINECA recipe:
https://itpplasma.github.io/SIMPLE/latest/hpc/index.html

- `Perlmutter/`, `aCluster/`, `sCluster/`: `submit_gpu.slurm`, one SIMPLE process on one GPU.
- `submit.slurm`, `submit_mpi.slurm`, `Marconi/`: legacy OpenMP templates for the old
  `BUILD/` layout; `simple_mpi.x` is no longer built.
- `alpha_lifetime_m.inp`, `datamodifier.f90`: legacy input and post-processing helpers.
