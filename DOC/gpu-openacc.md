# OpenACC GPU tracing

SIMPLE can offload the default orbit-tracing hot path to NVIDIA GPUs with
OpenACC. The offloaded configuration is the Boozer field (`isw_field_type=2`)
with the explicit-implicit symplectic Euler integrator
(`integmode=EXPL_IMPL_EULER`), no collisions, wall, or classifiers. Every other
configuration keeps the existing CPU path.

One particle runs per GPU thread. The integration is double precision, as the
symplectic scheme requires; the parallelism over particles is what the GPU
exploits, not single-orbit vectorization.

## Requirements

- NVIDIA HPC SDK (`nvfortran`), tested with 26.3.
- A CUDA-capable GPU. Tested on RTX 5060 Ti (cc120, Blackwell).
- NetCDF-Fortran and HDF5-Fortran built with `nvfortran` (compiler-specific
  `.mod` files; the distro packages are built with gfortran and are unusable).
  `nf-config` on `PATH` must point at the nvfortran build.

libneo must be built from a branch carrying the matching OpenACC fixes (managed
memory for batch splines, `acc routine seq` on the single-point evaluators, the
hdf5_tools nvfortran compile fix).

## Build

```bash
cmake -S . -B build-gpu -G Ninja \
  -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_C_COMPILER=nvc \
  -DCMAKE_BUILD_TYPE=Release -DSIMPLE_DETERMINISTIC_FP=ON \
  -DSIMPLE_ENABLE_OPENACC=ON -DENABLE_OPENACC=ON
cmake --build build-gpu -j
```

`SIMPLE_DETERMINISTIC_FP=ON` is required: the default `-fast -Mfprelaxed`
relaxes reciprocals and square roots, which breaks convergence of the
symplectic Newton iteration.

`SIMPLE_OPENACC_MEM` (default `unified`) selects the memory model. `unified`
needs an HMM-capable driver and keeps module statics plus allocatable spline
coefficients coherent across host and devices. `LIBNEO_OPENACC_MEM` must match.

The memory model matters for correctness, not only speed: without managed or
unified memory, `nvfortran` does not copy the allocatable `coeff` component of a
batch spline back from the device, so every spline evaluation returns zero.

## Multi-GPU

`trace_orbits_gpu` uses one device by default. `SIMPLE_GPU_NUM_DEVICES=N`
splits the particle batch across N devices, one host thread per device.

Splitting across cards is not yet profitable. With `mem:unified` the spline
coefficient array is shared, so two devices fault it back and forth on every
evaluation and throughput collapses. Profitable multi-GPU needs a per-device
resident copy of the read-only splines (construct or copy the coefficients once
per device), which the current batch spline construction does not do. Until
then, run one process per GPU over disjoint particle sets.

## Validation and benchmarking

Set `SIMPLE_GPU_BENCH=1` to run the GPU kernel and the CPU integrator on
identical per-particle initial states and report the agreement and speedup:

```bash
SIMPLE_GPU_BENCH=1 ./build-gpu/simple.x simple.in
```

`test_batch_splines_device` checks that batch spline evaluation inside an
OpenACC kernel matches the host result; it is registered when
`SIMPLE_ENABLE_OPENACC=ON`.

## Measured performance

RTX 5060 Ti (single card) versus a 32-thread Ryzen CPU, Boozer + Euler,
`trace_time=3e-4`:

| particles | CPU (OpenMP) | GPU    | speedup |
|-----------|--------------|--------|---------|
| 1024      | 1.74 s       | 5.92 s | 0.29x   |
| 4096      | 7.03 s       | 5.91 s | 1.19x   |
| 16384     | 27.6 s       | 23.0 s | 1.20x   |

The GPU saturates near 4096 particles and holds about 1.2x against the
32-thread CPU. Final-state positions agree to about 1e-5 and loss-step
classification matches for every particle; the small position difference is
floating-point reassociation in the parallel reductions, within the spread
expected for chaotic orbits.

The 1.2x ceiling on one card is the double-precision throughput of consumer
Blackwell (FP64 at 1/64 of FP32). Splitting across both cards scales the
particle batch further.
