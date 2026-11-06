# History audit 2022–2026

This audit rebuilt every first-parent commit of `main` since 2024 that
compiles, plus a sparse sample back to 2021, and ran each build on the same
fifteen physics cases, mostly alpha-particle losses in a quasi-axisymmetric
(QA) stellarator. It asked two questions: did any change move results
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
  (#683, #684). Ten further defects that the audit cases exposed are fixed
  as well (see [Other defects](#other-defects)).
- **Every result change in the audited cases since 2024 has a changelog
  entry except two in 2024**, which moved only fast-classification results
  (see the list below). The largest change before that, the drop of the QA
  loss fraction from 10.5 % to 6.9 % in 2022, is the removal of `bmod_ref`
  (the field rescaling to a reference value).
- **All coordinate sets now agree on the QA case within two markers.** Boozer
  and canonical-flux coordinates lose 71 of 1024 markers (6.93 %), Albert 72
  (7.03 %) and Meiss 73 (7.13 %). Albert lost 11.6 % until its radial grid
  was extended to the whole plasma (#728, #732). With 1024 instead of 256
  steps per field period, Meiss and Albert lose exactly the canonical-flux
  set, so the extra markers are step-size error of symplectic Euler.
- **CPU throughput of the Boozer symplectic trace is 1.2 million particle
  steps per second per core** (1.0 million in 2022). Boozer preprocessing grew
  from 1.1 s in 2022 to 10.3 s in early November 2026 and is 6.4 s now. Most
  of it is the Boozer setup in the libneo dependency, which production builds
  fetch unpinned and which takes about 20 times longer since libneo #461
  (libneo#489, open).
- **One GPU traces 16 384 markers for 10 ms in 47 s (RTX PRO 6000) or 19 s
  (H100)** with an adaptive Runge–Kutta integrator, against 144 s with
  symplectic Euler on 32 CPU threads. The GPU time was the same at all six
  commits measured since the native CUDA path was added in August 2026.

## Method

**Builds.** Each commit was built in Release mode as `make` builds it, with
the libneo revision that was current on `main` at the commit time (production
builds fetch libneo unpinned, so this is what a user got that day). Pinned
dependencies (fortnum, fortplot) were served from local mirrors. Commits
whose own sources do not compile are marked as not built.

**Cases.** All cases use the public QA equilibrium of the golden records (the
reference outputs of SIMPLE's regression tests) scaled to reactor size, or the NCSX equilibrium with its coils, or the
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

The change sets that moved the fate of at least one marker, oldest first.
Each has a changelog row unless stated otherwise.

- **2021-08 to 2021-12** (77722c24, abfd48ca, 6a6699ab): ANTS starts and the
  first Boozer symplectic integrator; canonical-flux confinement 93.8 % to
  93.0 %. This predates the changelog.
- **2024-03-11** (6c5eb2ca): `bmin` and `bmax` from a local initialization
  change the trapped–passing parameter; collisions and fast classification
  change.
- **2024-05-25 and 2024-07-18** (579874be "Shorten module name", d582580f
  "Configurable collisions"): fast classification only, 60 to 99 markers
  change class. *No recorded reason*; neither commit intends a change of
  classification, and the audit did not trace the cause.
- **2024-08** (87a2afc7, #49; 74e871f4, #55): end states zeroed and OpenMP
  counters reorganized; `fast_class` alone ignored from #55 on (restored by
  #512 in 2026-10).
- **2025-03 and 2025-08** (c6de2792, #77; c7ca093a): sampler rewrite, after
  which fast classification reported every marker confined until c7ca093a
  initialized the starting surface again.
- **2025-08-22** (0016bff6, #172): RK45 switched to the libneo Cash–Karp
  integrator and Release builds to `-ffast-math`; Boozer RK45, Albert and
  coil-field Meiss change.
- **2025-12-17** (bc67a39b, #257): unified VMEC field convention; Albert
  changes by two markers.
- **2026-05-22** (e404496e, #323): `fast_class` stops regular orbits early
  again.
- **2026-06** (1eaa91cc, #354; 78bd4814, #384): Boozer axis handling and
  symplectic orbits continued through the axis; TEST tokamak only.
- **2026-07-13 to 2026-07-25** (714c1686, #454; cc5872a8, #459; 6eb1542f,
  #487; ba581b6f, #488; f60c57c4, #494; d8b0e019, #499; 69fe8c36, #503):
  failure reporting, located boundary crossings and warning-mode recovery.
  Midpoint, Gauss and Lobatto confinement fell to 57–82 % for one day and was
  restored by #487; every case changes at least in loss times.
- **2026-10-16 to 2026-10-21** (ffd7f4a1, #540; 27e5d3b4, #558; aa1d30bd,
  #512; 9eca50ed, #560): the momentum gate on symplectic endpoints removed and
  restored, collision streams per marker, `fast_class` alone enabled again.
- **2026-10-27 to 2026-10-31** (0a176137, #629; 094b7d25, #638; 4feb1bed,
  #640; 00a2b30f, #642): RK-recovery seeding, collision kicks with the carried
  field, Newton iterates through the axis.
- **2026-11-02** (b6dc7f9b, #643): the commit itself reproduces its parent;
  the change is libneo #461 of 2026-10-31, which production builds pick up
  unpinned (collisional losses 8.9 % to 8.4 %).
- **2026-11-03 to 2026-11-05** (b3a5aafb, #657; af9c3346, #647; 359ec8c4,
  #674; 619c217a, #675; 74eba10a, #707): round-off series, species-resolved
  Coulomb logarithms, CODATA 2018 constants, the 20 % momentum gate.
- **2026-11-05** (bcfac3c3, #732): Albert radial grid over the whole plasma,
  Albert losses 11.6 % to 7.0 %.

The TEST tokamak case moves between 0 and 63 lost markers of 8192 (0.8 %)
with each change of the recovery rules (#503, #540, #560 and the paired
benchmark commit 850ff9a2), and since #708 358 of its markers end
unresolved. This is a property of the test case, not of the code: its
normalized Larmor radius $\rho = 1$ is twice the minor radius, so
$B^*_\parallel = B + \rho\, v_\parallel\, \mathbf{h}\cdot(\nabla\times\mathbf{h})$
(with $\mathbf{h}$ the unit vector along the field) reaches about $-0.06$ for
some counter-passing markers and the guiding-centre
equations become singular there. Before #708 the RK field had no curl term,
so recovery did not meet the singularity. Which markers reach it depends on
rounding, so builds with different instruction-set flags differ.

## Regressions found

| Issue | Introduced | Effect | Fix |
| ----- | ---------- | ------ | --- |
| #681 | Meiss/Albert coordinates (2024-09); Boozer chartmaps (2026) | Classification (`tcut > 0`, `fast_class`, `class_plot`) traced Meiss, Albert and Boozer-chartmap markers from untransformed reference coordinates and wrote integrator coordinates as end points. | #682 |
| #683 | 5a8f9d76 (2024-09) | The starting field line evaluated the full VMEC field with derivatives at every point: 3.9 s of a 10 s Boozer preprocessing on the QA case. | #684 (0.6 s) |
| #685 | c6de2792 (#77, 2025-03) | ANTS starts (`startmode = 3`) reopened the input for every marker, so all markers started from the first line. | #696 |
| #697 | 9eca50ed (#560, 2026-10-21) | A 0.5 % bound on the momentum of accepted symplectic steps, below the O(Δt) energy error of symplectic Euler, rejected ordinary steps: 68 000 to 554 000 Runge–Kutta recoveries per run in canonical-flux, Meiss and Albert coordinates, Meiss losses 7.0 % → 7.7 %, Albert 11.7 % → 13.6 %, TEST tokamak 0.8 % → 0 %, Boozer trace 47 % slower. | #707: bound 20 %, Meiss 7.1 %, Albert 11.9 % |
| #700 | #522 with #544 (2026-10) | Fusion-birth starts checked the default energy factor instead of the effective particle energy, so 10 keV markers passed as 3.5 MeV alphas. | #701 |

### Other defects

Defects that the audit cases exposed and that no single recent change
introduced:

| Issue | Effect | Fix |
| ----- | ------ | --- |
| #728 | The Albert radial grid ended at the smallest edge value of its radius, so the field was extrapolated beyond $s \approx 0.73$ at some angles: 48 extra lost markers on the QA case. | #732 |
| #744 | REFCOORDS on a VMEC reference traced a marker given at $s$ from the surface $s^2$. | #750 |
| — | RK in VMEC coordinates collapsed its step when an angle passed through zero. | #747 |
| #754 | `magfie_vmec` radial derivatives returned NaN within $s < 10^{-3}$ of the axis and were too small near it; RK markers there ended unresolved. | #752, #755 |
| #508 | Markers carried by thousands of RK recoveries were reported as physical losses. | #753 (`max_rk_recoveries`) |
| — | Classification reported markers with a non-finite boundary state as LCFS losses. | #746 |
| — | Classification in Albert coordinates did not stop markers at $s = 1$. | #762 |
| — | Gauss and Lobatto steps across the axis kept $r < 0$ (NaN $R$, $Z$ in `orbits.nc`). | #763 |
| — | Boozer starts were moved to the first field period. | #733 |
| — | Meiss and Albert volume starts were not uniform in $s$. | #706 |

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
extrapolation. These numbers are from before #732. With the grid extended to
the largest edge value (#732),
Albert loses 72 markers (7.03 %) at 256 steps per field period and exactly
the canonical-flux set at 1024 steps. Every commit merged since then gives
the same 72 markers in this case.

## Performance

```{figure} perf_timeline.png
:alt: Preprocessing time, trace time and throughput per case against commit date
Paired single-thread timings, 32 markers for 1 ms.
```

| Case | Preprocessing 2022-06 | Preprocessing now | Trace throughput 2022-06 | Trace throughput now |
| ---- | --------------------- | ----------------- | ------------------------ | -------------------- |
| Boozer | 1.1 s | 6.4 s | 1.00 M steps/s | 1.22 M steps/s |
| Canonical flux | 2.4 s | 2.6 s | 0.54 M steps/s | 0.70 M steps/s |
| Collisions (Boozer) | 1.1 s | 6.3 s | 0.91 M steps/s | 0.81 M steps/s |
| Meiss | — | 13.0 s | — | 0.56 M steps/s |
| Albert | — | 17.5 s | — | — |

"Now" is 30431224 (2026-11-06); Meiss and Albert did not exist in 2022.

**Preprocessing.** Boozer preprocessing went from 4.5 s to 10.3 s with
libneo #461 (2026-10-31), which production builds picked up unpinned
(libneo#489, still open). It fell from 10.1 s to 6.4 s between b420006b and
7460ae9d, a range of ten commits without timings that contains the
starting-field-line fix #684 (3.9 s to 0.6 s when measured alone).

**Trace.** The trace got faster in three steps in November 2026, each larger than
the scatter of about 10 % between neighbouring commits: Newton
started from extrapolated states (#646; Boozer 3.2 s → 2.5 s, canonical flux
6.2 s → 4.8 s, Meiss 7.8 s → 5.7 s), the canonical-flux field stored as one
spline (#647; canonical flux 4.8 s → 3.3 s, Boozer 2.7 s → 2.0 s) and the
20 % momentum gate without the extra field evaluation of #560 (#707; Boozer
−25 %, canonical flux −23 %, Meiss −40 %). In Albert coordinates, steps that
end inside the plasma skip the radius map (#760): 8.2 % fewer instructions
on the QA case, counted on Pitagora with both builds on one libneo
revision, with identical loss times; its effect on the trace time is within the scatter.

**RK45 in VMEC coordinates** takes 35–45 min for the 1024-marker case and
exceeded the audit's time limit (30 min in the first runs, 90 min in later
ones) on several commits (status `FAIL 124` in `results.tsv`). These are runs that
did not finish in time, not failed markers.

### GPU

```{figure} gpu_throughput.png
:alt: Traced markers per second against number of markers, two GPUs and CPU
Native CUDA production mode on one GPU against the CPU on 32 threads.
```

The native CUDA production mode (`gpu_mode = 'production'`, Dormand–Prince
Runge–Kutta on a Boozer chartmap) traced 16 384 markers for 10 ms in 47 s on
one RTX PRO 6000 GPU at each of the six commits measured since it was added
(August 2026 to November 2026, all within 0.3 %). The kernel time
hardly depends on the marker count up to 16 384 (47.6 s for 1024 markers), so
small runs do not fill the GPU; 65 536 markers take 91 s, about 720 markers
per second. One H100 of the Pitagora cluster at CINECA (Booster partition) runs the same input in
19.4 s for 1024 markers, 19.1 s for 16 384 and 25.7 s for 65 536 (about
2500 markers per second), 2.4 to 3.6 times faster than the RTX PRO 6000. The
H100 ran 29eaf52d and the RTX PRO 6000 b420006b; given how little the GPU
time changed over six revisions, the difference is the hardware. The CPU
traces the same markers for the same time with symplectic Euler on 32
threads in 144 s. The two use different integrators, so this compares time
to solution, not cost per step.

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

The issues on GitHub name the files and lines concerned.

| Issue | Topic | State |
| ----- | ----- | ----- |
| #686 | `macrostep` and `macrostep_with_wall_check` duplicate the microstep loop, so every recovery fix had to be made twice | fixed by #729 |
| #687 | Code that `simple.x` never reaches (quasi-symplectic implicit steppers, `rk_gauss_newton`, `array_utils`, `test_collis`) | fixed by #709, #759 |
| #688 | The damped Newton iteration is written out separately for each implicit scheme, CPU and GPU | fixed by #731, #745 |
| #689 | Hand-written field derivative chains that the other kernels already generate with fortsym | fixed by #730 |
| #690 | REFCOORDS field calls build the coordinate basis twice | fixed by #711, #758 |
| #691 | The Albert boundary check maps both step endpoints fully | fixed by #712, #760 |
| #692 | Per-marker integrator setup wrote shared procedure pointers inside the OpenMP region | fixed by #698 |
| #693 | Classification loop: static schedule, critical sections and a counter race | fixed by #699, #757 |
| #694 | `simple_main.f90` (3700 lines) mixes CPU and GPU drivers and four marker loops | open; #762 removed one divergence between the classifier and the plain loop |
| #695 | Collision microsteps recompute radius and plasma state several times | open; #761 maps the radius once per microstep |
| libneo#489 | Boozer setup about 20 times slower since libneo #461 | open |

## Open items

- Split `simple_main.f90` and share one marker loop between plain tracing
  and classification (#694).
- Evaluate the plasma state once per collisional microstep (#695).
- Restore the Boozer setup time in libneo (libneo#489), or pin libneo in
  production builds as the golden records already do.
- Give the TEST tokamak case a normalized Larmor radius well below one, so
  that it stays inside guiding-centre validity.
- Find the cause of the 2024 fast-classification changes at 579874be and
  d582580f.

## Data

- `results.tsv`: one row per case and commit (status, wall time, confined
  fraction, lost markers, checksums of the output files).
- `changes.tsv`: one row per case and consecutive commit pair with the
  comparison class.
- `perf.tsv`: paired CPU timings.
- `gpu.tsv`: GPU and matching CPU timings; `gpu_h100.tsv`: the H100 runs.
- `ci_golden.tsv`: duration of the golden-record CI step.
