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
   Submitted as [libneo #427](https://github.com/itpplasma/libneo/pull/427).

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

## Boozer GC

A Boozer GC profile (production default) showed 27.5% in the generic-degree
single-quantity `der2` kernel although libneo has a quintic `nq1_o555`
kernel, and about 19% in heap allocation and copies from the single-point 1D
`der2`/`der3` wrappers. libneo #427 routes both to direct paths. With 2,560
markers and 2 ms, Boozer GC goes from 19.7 to 16.7 core-s per particle-s
(1.18x; midpoint with `npoiper2=1024`: 142.8 to 124.5, 1.15x). `times_lost.dat`
and `confined_fraction.dat` are byte-identical. Against the faster Boozer
baseline the NCSX hybrid is about 39x and FO everywhere about 137x.
`summary_fast_fo.json` and `summary_g3_fo.json` hold the cases; the G3V3L3Fi
FO reference timing predates both changes.

## Other field types (VMEC NCSX, 201 surfaces)

`test/test_data/wout_ncsx.nc`, 3.5 keV alphas, 2,560 markers, 10 ms, production
defaults (Euler1, `npoiper2=256`), sampled for 12 s during tracing. Inclusive
shares of `trace_orbit`:

| Field | Wall [s] | RK recoveries | Symplectic steps | RK recovery | Spline evaluation |
| --- | ---: | ---: | ---: | ---: | ---: |
| Canonical flux (0) | 71 | 7,349 | 83% | 9% | 70% |
| Boozer (2) | 47 | 552 | 90% | 1% | 59% |
| Meiss (3) | 147 | 57,491 | 19% | 75% | 76% |
| Albert (4) | 109 | 36,918 | 6% | 79% | 62% |

Meiss and Albert spend three quarters of their tracing time in the RK45
recovery that follows failed symplectic steps (`orbit_timestep_recovery`),
about 100x more recoveries than Boozer. On the 16-surface
`wout_li383_low_res_reference.nc` the counts are higher still (163,061 and
144,691 for 2,560 markers and 20 ms). This is the largest remaining cost for
these coordinates and is not a spline problem.

Cause: the recoveries were not Newton failures. The 0.5% momentum-drift gate
for RK recovery states also ran on converged symplectic steps, where z(4) is
reconstructed from H, which Euler1 conserves only to O(dt) (up to 2% at
`npoiper2=256`; rejections scale away at 1024). Symplectic steps now keep the
phase-space checks without that gate, and recovery starts from the reference
momentum so RK, validation and LCFS interpolation share one state. NCSX,
2,560 markers, 10 ms, 16 threads, matched starts: recoveries 34,084 to 1
(Meiss) and 33,471 to 1 (Albert); tracing 25-29 s to 18 s and 23-25 s to
17-19 s (four runs each, libneo f6090e8).
Losses 674/691 become 653/664 against 650 for Boozer, and paired status
differences against Boozer fall from 112/95 to 67/52.

Spline kernels: canonical flux and Boozer use the quintic single-quantity
`der2` kernel, Meiss and Albert the new quintic first-derivative kernel.

## gfortran IEEE save/restore

Every profile also showed `_gfortran_ieee_procedure_entry`/`_exit` with
`get/set_fpu_except_flags` (4% Boozer to 17% Albert of samples). gfortran
adds this save/restore to each procedure that accesses an IEEE intrinsic
module, including through a procedure-local `USE` of any module whose
dependencies import one at module scope (minimal reproduction: a procedure
with `use m, only: helper`, where `m` imports `ieee_arithmetic`). Module-scope
`USE` does not trigger it. Through fortio, libneo's `boozer_chartmap_io` and
SIMPLE's own modules, 96 modules carried an IEEE import, and `velo_can`,
`eval_field_can`, `validate_rk_state` and `integ_to_ref_albert` paid it on
every call.

SIMPLE now uses `finite_check` (`is_finite` by exponent bits, `quiet_nan`)
instead of the IEEE module; fortio#75 and libneo#428 move their imports into
the procedures that need them. Afterwards no tracing procedure carries the
hooks. Unprofiled A/B (16-surface equilibrium, 2,560 markers, 2 ms, three runs
each on a loaded host): Albert wall 31.4 to 27.1 s (about 1.2x in tracing),
canonical flux, Boozer and Meiss within noise. The profiler attribution
overstated the cost except for Albert, whose coordinate transform runs in
the inner loop. All outputs are byte-identical.
