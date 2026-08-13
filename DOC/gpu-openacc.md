# GPU tracing

SIMPLE's GPU production path is configured in the `config` namelist in
`simple.in`. The environment does not select the executable path or integrator.

```fortran
gpu_mode = 'production'       ! off, production, or benchmark
gpu_backend = 'auto'          ! native CUDA if built, otherwise OpenACC
gpu_method = 'auto'           ! tuned DOPRI on native CUDA, Euler on OpenACC
gpu_landreman = .True.        ! segmented alpha-loss objective
gpu_start_coordinates = 'reference'
```

The production input must use Boozer guiding-center tracing
(`isw_field_type = 2`) without collisions, walls, orbit output, or classifier
options. `relerr` is used as the integration tolerance. The supported methods
are:

| Backend | Methods | `gpu_method = 'auto'` |
| --- | --- | --- |
| OpenACC | Euler, midpoint, Cash-Karp, tuned DOPRI | explicit-implicit Euler |
| Native CUDA | Cash-Karp, tuned DOPRI | tuned DOPRI |

## OpenACC build

Requirements:

- NVIDIA HPC SDK (`nvfortran`, `nvc`), tested with 26.3.
- A CUDA-capable GPU.
- NetCDF-Fortran and HDF5-Fortran built with `nvfortran`.

The one feature switch propagates matching OpenACC settings to libneo and
Fortnum. The default memory model is `managed`.

```bash
cmake -S . -B build-gpu -G Ninja \
  -DCMAKE_Fortran_COMPILER=nvfortran \
  -DCMAKE_C_COMPILER=nvc \
  -DSIMPLE_ENABLE_OPENACC=ON
cmake --build build-gpu -j
(cd run && ../build-gpu/simple.x simple.in)
```

Set `SIMPLE_OPENACC_MEM` at configure time only when the target system needs a
different memory model. SIMPLE passes the same value to libneo.

## Native CUDA build

```bash
cmake -S . -B build-cuda -G Ninja -DSIMPLE_ENABLE_CUDA=ON
cmake --build build-cuda -j
(cd run && ../build-cuda/simple.x simple.in)
```

The native path uses Boozer RK field tables and reports allocation, upload,
kernel, download, and C-ABI timings. `gpu_compact_init = .True.` enables its
compact Boozer initializer. Use it only with native CUDA and Boozer starts.

## Production controls

The following namelist values have defaults suitable for the production path:

```fortran
gpu_num_devices = 1
gpu_t_block = 1d-3
gpu_loss_tau = 0.1d0
gpu_maxloss = 0.02d0
gpu_min_timestep = 1d-10
gpu_pilot_fraction = 0.04d0
gpu_work_order = .True.
gpu_particle_profile = ''
```

`gpu_num_devices > 1` is experimental. With managed memory, a shared spline
table can migrate between devices, so one process per GPU is usually better.
`gpu_particle_profile` writes initial and final states, field-evaluation counts,
loss steps, and loss times.

## Benchmark mode

Set `gpu_mode = 'benchmark'` in `simple.in` with `integmode = 1` (Euler) or
`integmode = 3` (midpoint). The benchmark runs the OpenACC GPU and CPU paths on
identical particles and reports agreement and speedup. It does not support
walls, collisions, classifier options, or boundary-reaching traces.

## Validation

`test_batch_splines_device` checks device spline evaluation against the host
result. `test_gpu_orbit_bench` and `test_gpu_midpoint_bench` exercise benchmark
mode. `test_gpu_landreman_segments` checks the segmented accounting and the
native CUDA bridge when the corresponding GPU build is available.

On an RTX 5060 Ti against a 32-thread Ryzen CPU, the OpenACC Euler path reached
about 1.2x at 4,096 and 16,384 particles for the documented short benchmark.
The result depends on particle count, GPU, compiler, and memory model.
