# SPECTRE orbit profile and spline/exit-check speedups, 2026-10-16

macOS `sample` (1 ms) attached to the `study_spectre_orbits` executable
during tracing; NCSX, 3.5 keV, 256 matched starts, 16 threads, deterministic
Fast build. `leaf_*.txt` are the "top of stack" tables for the hybrid (entry
1.5, `npoiper2=512`), FO everywhere, and an almost pure GC run (entry 0.01).

## Profile before the changes

| Leaf | Hybrid | FO | Near-pure GC |
| --- | ---: | ---: | ---: |
| `evaluate_batch_splines_3d_der_core_general` | 74.1% | 80.7% | 36.0% |
| `evaluate_batch_splines_3d_der2_core_o555` | 6.5% | – | 16.4% |
| `memmove` / `memcpy` stub inside the spline kernel | 12.3% | 13.7% | 6.0% |
| idle threads (`psynch_cvwait`) | – | – | 34.6% |

Inclusive shares of marker-loop time:

| Routine | Hybrid | FO | Near-pure GC |
| --- | ---: | ---: | ---: |
| FO total (`continue_fo`) | 90% | 100% | 61% |
| RK45 stages (`dp_trial`) | 53% | 59% | 36% |
| exit check (`gc_map`) | 37% | 41% | 25% |
| of which `spectre_fo_field_locate` | 28% | 31% | 19% |
| GC steps (`orbit_timestep_sympl`) | 8% | – | 31% |
| proximity check (`spectre_fo_radial_larmor`) | 1.9% | – | 6.8% |

## Changes

1. **Exit check once per radian of gyration** (`spectre_fo_hybrid`,
   `EXIT_CHECK_ANGLE`). The guiding-centre map with its Newton coordinate
   inversion ran after every RK45 substep (at most 0.1 rad of gyration) only to
   test the exit distance. It now runs once per radian and at the end of every
   call, so the returned state and `last_y` remain current at microstep
   boundaries. An exit is detected at most 1 rad late, small against a
   one-gyroradius hysteresis band. FO-everywhere orbits are unchanged.
2. **Quintic first-derivative spline kernel** in libneo
   (`evaluate_batch_splines_3d_der_core_o555`, libneo branch
   `perf/spline3d-der-o555`). The general kernel's runtime loop bounds
   produced coefficient copies that gfortran emitted as library `memmove`
   calls. The new kernel uses fixed degree and register recurrences and is
   bitwise identical to the general kernel with FMA contraction off.

## Results (256 markers, 2 ms)

| Case | Losses | Status diff vs FO | Core-s / particle-s before | after | Speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| NCSX FO everywhere | 101 | – | 11,629 | 2,280 | 5.1 |
| NCSX hybrid, entry 1.5 | 99 | 6 (was 12) | 2,178 | 657 | 3.3 |
| G3V3L3Fi hybrid, entry 1.5 | 53 | 0 (was 0) | 632 | 162 | 3.9 |

FO-everywhere marker outcomes are bitwise identical before and after. The
hybrid changes only through the exit timing; its loss-curve maximum difference
from FO is 0.020 (NCSX) and 0.012 (G3V3L3Fi), within the 256-marker noise.
The tok2vol p_zeta invariant stays at 8.3e-8 (tolerance 1e-5), and the
affected native tests pass 24/25 (`test_spectre_construction_grid` needs
Linux `/proc`). Quick 64-marker runs isolate the two changes: the exit check
alone gives 1.57x (FO) and 1.44x (hybrid); the kernel then adds 3.5x and 3.0x.

Against production Boozer GC (19.7 core-s per particle-s, `../boozer-cost/`)
the NCSX hybrid is now about 33x and FO everywhere about 116x.
`summary_fast_fo.json` and `summary_g3_fo.json` hold the cases; the G3V3L3Fi
FO reference timing predates both changes.
