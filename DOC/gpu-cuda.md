# GPU tracing

SIMPLE's supported production GPU path is native CUDA. Configure it in the
`config` namelist in `simple.in`; no environment variable selects the backend
or integrator.

```fortran
gpu_mode = 'production'
gpu_backend = 'cuda-native'
gpu_method = 'dopri'
gpu_landreman = .True.
```

Use Boozer guiding-center tracing (`isw_field_type = 2`) without collisions,
walls, orbit output, or classifier options. Native CUDA supports tuned DOPRI
and Cash-Karp. DOPRI is the default. Both use `relerr` as their tolerance.

## Build and run

Requirements are the NVIDIA CUDA toolkit, a C++17 compiler, GNU or NVIDIA
Fortran, and the normal SIMPLE dependencies.

```bash
cmake -S . -B build-cuda -G Ninja -DSIMPLE_ENABLE_CUDA=ON
cmake --build build-cuda -j
(cd run && ../build-cuda/simple.x simple.in)
```

`gpu_compact_init = .True.` enables the compact Boozer initializer. Use it only
with native CUDA and Boozer starts. `gpu_start_coordinates = 'boozer'` means
that `start.dat` already contains Boozer coordinates.

## Production controls

The defaults are:

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

`gpu_particle_profile` writes initial and final states, field-evaluation
counts, loss steps, and loss times. Multiple devices are experimental; use one
process per GPU unless device-local field tables have been arranged.

The native CUDA tests are `test_cuda_native_rk54` and
`test_cuda_native_bridge`; the latter checks segmented loss accounting through
the Fortran-to-CUDA production bridge.
