# History audit 2022–2026

This audit rebuilt every first-parent commit of `main` since 2024 that
compiles, plus a sparse sample back to 2021, and ran each build on the same
fifteen physics cases. It asked two questions: did any change move results
without a reason, and how did run time develop. The answers, with the
measurements behind them, are below; raw tables are next to this page.

## Summary

- **Five regressions were still in `main` and are fixed.** Classification
  traced Meiss, Albert and Boozer-chartmap markers from untransformed start
  points (#681, fixed by #682); ANTS starts gave every marker the first input
  line (#685, #696); the 0.5 % momentum gate on symplectic endpoints returned
  and caused storms of Runge–Kutta recoveries and extra losses (#697, #707);
  fusion births accepted any particle energy while reporting 3.5 MeV alphas
  (#700, #701); and the starting field line cost seven times more than needed
  (#683, #684).
- **Every other result change in the audited cases has a recorded reason.**
  The largest, the drop of the QA loss fraction from 10.5 % to 6.9 % in 2022,
  is the removal of `bmod_ref` (the field rescaling to a reference value),
  documented in the changelog. A few round-off and loss-time-only changes
  were missing from the changelog and are added.
- **Symplectic Euler in Meiss and Albert coordinates loses more markers than
  in Boozer coordinates** (7.1 % and 11.9 % against 6.9 % on the QA case).
  This is step-size error of the first-order scheme, not a regression:
  midpoint, Gauss and finer steps remove the difference.
- **CPU throughput of the Boozer symplectic trace is about 1 million
  particle steps per second per core** and has stayed within a factor of
  two since 2024. Preprocessing grew from 1.2 s (2024) to about 10 s now; the
  last step, from 4.4 s to 9.8 s, comes from the unpinned libneo dependency
  (libneo#489).
- **The native CUDA path traces 16 384 markers for 10 ms in 47 s** on one
  RTX PRO 6000 (about 350 markers per second), three times faster than 32
  CPU threads; it has not changed since it was added in August 2026.

## Method

**Builds.** Each commit was built in Release mode as `make` builds it, with
the libneo revision that was current on `main` at the commit time (production
builds fetch libneo unpinned, so this is what a user got that day). Pinned
dependencies (fortnum, fortplot) were served from local mirrors. Commits
whose own sources do not compile are marked as not built.

**Cases.** All cases use the public QA equilibrium of the golden records
scaled to reactor size, or the NCSX equilibrium with its coils, or the
analytic TEST tokamak. The base case traces 1024 alpha markers started at
normalized toroidal flux $s = 0.5$ for 10 ms with symplectic Euler, 256 steps
per field period and 1000 output times. The variants change one thing:
canonical flux, Boozer, Meiss or Albert coordinates; RK45 (adaptive
Runge–Kutta) in Boozer and VMEC coordinates; midpoint, Gauss and Lobatto
integrators; collisions; classification (`tcut`, `fast_class`); the TEST
tokamak; start sampling; and Meiss coordinates on a coil field. For revisions
before the namelist existed the same physics was set positionally. Before
2022-06-29, SIMPLE rescaled the field to `bmod_ref`; for those revisions
`bmod_ref` was set to the field the later code uses, so that the cases stay
comparable.

**Comparison.** For each pair of consecutive builds the audit compares the
confined fraction, the fate of each marker (lost or confined) and the loss
times. A change is called *round-off* if loss times differ by less than
$10^{-6}$ relative, *loss-time* if no marker changes fate, and *changed* if at
least one does.

**Timing.** CPU timings use the paired benchmark of
`test/performance/benchmark.py` (#679): each revision and the reference (`main` at the
start of the audit) start at the same moment on one node, one OpenMP thread
each, so both see the same load; 32 markers for 1 ms, five repeats. A
revision's reported time is its median time ratio to the reference times
the reference's median time over all pairs, which removes slow drifts in
node speed. Preprocessing is the time to the
start of tracing; trace time is the rest. Throughput is particle steps per
second, with $(n_\mathrm{timstep}-1)\,n_\tau$ steps per marker.

## Results over time

```{figure} losses_timeline.png
:alt: Lost marker fraction per case against commit date
Lost fraction after 10 ms per case and commit (QA reactor-scale case).
```

The change sets that moved results, oldest first:

CHANGESETS

## Regressions found

| Issue | Introduced | Effect | Fix |
| ----- | ---------- | ------ | --- |
| #681 | Meiss/Albert coordinates (2024-09); Boozer chartmaps (2026) | Classification (`tcut > 0`, `fast_class`, `class_plot`) traced Meiss, Albert and Boozer-chartmap markers from untransformed reference coordinates and wrote integrator coordinates as end points. | #682 |
| #683 | 5a8f9d76 (2024-09) | The starting field line evaluated the full VMEC field with derivatives at every point: 3.9 s of a 10 s Boozer preprocessing on the QA case. | #684 (0.6 s) |
| #685 | c6de2792 (#77, 2025-03) | ANTS starts (`startmode = 3`) reopened the input for every marker, so all markers started from the first line. | #696 |
| #697 | 9eca50ed (#560, 2026-10-21) | A 0.5 % bound on the momentum of accepted symplectic steps, below the O(Δt) energy error of symplectic Euler, rejected ordinary steps: 68 000 to 554 000 Runge–Kutta recoveries per run in canonical-flux, Meiss and Albert coordinates, Meiss losses 7.0 % → 7.7 %, Albert 11.7 % → 13.6 %, TEST tokamak 0.8 % → 0 %, Boozer trace 47 % slower. | #707: bound 20 %, Meiss 7.1 %, Albert 11.9 % |
| #700 | #522 with #544 (2026-10) | Fusion-birth starts checked the default energy factor instead of the effective particle energy, so 10 keV markers passed as 3.5 MeV alphas. | #701 |

## Meiss and Albert coordinates against canonical flux

Canonical-flux and Boozer coordinates lose the same 71 markers (6.93 %) at
every step size tested. On `main` after #707, Meiss coordinates lose two more
markers at 256 steps per field period (7.13 %); with 1024 or 4096 steps, or
with the midpoint integrator, Meiss loses exactly the canonical-flux set. The
two extra markers are step-size error of symplectic Euler, not a defect. The
higher Meiss values of late October 2026 (8.98 % and 7.71 %) were the recovery
storms of #697.

Albert coordinates lost 11.4–11.7 % at every step size and integrator since
they were added: 48 markers more than canonical flux, none fewer, most lost
within 0.1 ms. The cause (#728) is the radial grid of the Albert field, which
ended at the smallest edge value of its radial coordinate over all angles,
so beyond about $s = 0.73$ at some angles the field was spline
extrapolation. With the grid extended to the largest edge value, Albert loses
7.03 % at 256 steps per field period and exactly the canonical-flux set at
1024 steps. ALBERTPR

## Performance

```{figure} perf_timeline.png
:alt: Preprocessing time, trace time and throughput per case against commit date
Paired single-thread timings, 32 markers for 1 ms.
```

PERF

### GPU

```{figure} gpu_throughput.png
:alt: Traced markers per second against number of markers, GPU and CPU
Native CUDA production mode on one GPU against the CPU on 32 threads.
```

The native CUDA production mode (`gpu_mode = 'production'`, Dormand–Prince
Runge–Kutta on a Boozer chartmap) traced 16 384 markers for 10 ms in 47 s on
one RTX PRO 6000 GPU at every commit since it was added (August 2026 to
November 2026, six revisions measured, all within 0.3 %). The kernel time
hardly depends on the marker count up to 16 384 (47 s for 1024 markers), so
small runs do not fill the GPU; 65 536 markers take 91 s, about 720 markers
per second. The CPU traces the same markers for the same time with
symplectic Euler on 32 threads in 144 s. The two use different integrators,
so this compares time to solution, not cost per step.

The OpenACC benchmark mode (`gpu_mode = 'benchmark'`, which runs the CPU
and the GPU on the same markers and compares them) stopped with an illegal
memory access on `main` because two targets were compiled with different
GPU memory models; #710 fixed this. Older OpenACC revisions do not build with
the NVHPC 26.3 compiler available for this audit, so no OpenACC timeline
exists.

### Golden-record run times

```{figure} golden_times.png
:alt: Run time of the golden-record inputs per commit, and CI duration of the golden records
Top: golden-record inputs run with one thread, paired against `main`.
Bottom: duration of the golden-record CI step.
```

Run single-threaded and paired against `main`, the golden-record inputs take
6–18 s each today (Boozer 9.6 s, canonical flux 6.0 s, Meiss 16 s, Albert
18 s, classifier 16 s, collisions 9.4 s, TEST tokamak 0.1 s). Most of this is
preprocessing, which is why the Boozer and collision inputs went from 1.3 s
in May 2024 to 7.3 s in August 2025 and back up to about 10 s with the
libneo update of October 2026. In CI the golden-record step took about 20 min
on GitHub runners in spring 2026, 28 min by August, and 40–70 min in mid
October 2026; since the gate moved to sCluster the whole build and test gate,
golden records included, typically takes 6–10 min and up to 45 min under load. The CI numbers depend on runner
load and on how many cases existed, so they show the order of magnitude, not
changes of individual commits.

## Code-quality issues

Each issue names the files and lines concerned.

| Issue | Topic | State |
| ----- | ----- | ----- |
| #686 | `macrostep` and `macrostep_with_wall_check` duplicate the microstep loop, so every recovery fix had to be made twice | ST686 |
| #687 | Code that `simple.x` never reaches (quasi-symplectic implicit steppers, `rk_gauss_newton`, `array_utils`, `test_collis`) | ST687 |
| #688 | The damped Newton iteration is written out separately for each implicit scheme, CPU and GPU | ST688 |
| #689 | Hand-written field derivative chains that the other kernels already generate with fortsym | ST689 |
| #690 | REFCOORDS field calls build the coordinate basis twice | ST690 |
| #691 | The Albert boundary check maps both step endpoints fully | ST691 |
| #692 | Per-marker integrator setup wrote shared procedure pointers inside the OpenMP region | ST692 |
| #693 | Classification loop: static schedule, critical sections and a counter race | ST693 |
| #694 | `simple_main.f90` (3700 lines) mixes CPU and GPU drivers and four marker loops | ST694 |
| #695 | Collision microsteps recompute radius and plasma state several times | ST695 |
| libneo#489 | Boozer setup about 20 times slower since libneo #461 | ST489 |

## Open items

OPEN

## Data

- `results.tsv`: one row per case and commit (status, wall time, confined
  fraction, lost markers, checksums of the output files).
- `changes.tsv`: one row per case and consecutive commit pair with the
  comparison class.
- `perf.tsv`: paired CPU timings.
- `gpu.tsv`: GPU and matching CPU timings.
- `ci_golden.tsv`: duration of the golden-record CI step.
