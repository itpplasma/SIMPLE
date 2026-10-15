# SPECTRE loss diagnosis, 2026-10-15

The 48-versus-53 G3V3L3Fi discrepancy contains a confirmed absorbing-boundary
inconsistency. Resolving the particle wall in the multi-volume hybrid recovers
four missed losses and gives **52/256 on both FO spline grids**. The remaining
marker is sensitive to the current switch and phase model. It is not accepted
as a validated physical difference or as evidence that the hybrid is converged.

The focused production changes are merged on main:
[FO performance #534](https://github.com/itpplasma/SIMPLE/pull/534),
[GC extrapolation #535](https://github.com/itpplasma/SIMPLE/pull/535), and
[particle-wall repair #536](https://github.com/itpplasma/SIMPLE/pull/536).
They merged as `934d964`, `71ad9ae` and `27e963c`, respectively. The frozen
measurement baseline remains `fe002a3`. Chris explicitly requested promotion;
independent human domain/senior review was not performed and CI was not awaited.
The requested [Sonnet reviews and adjudication](REVIEW.md) are advisory. They
found no confirmed new regression; the pre-existing midpoint axis-chart concern
remains open and is not excluded by the retained loss controls.

## Boundary repair

The previous proximity screen excluded the outer boundary. GC segments therefore
absorbed their guiding centre at `rho_g = Mvol`, while FO segments absorbed the
particle. A Larmor orbit can hit that same surface with its centre still inside.
Three of the five differing markers never entered FO in the old hybrid.

The repair applies the existing three-gyroradius entry and four-gyroradius exit
distances to the outer edge. The existing FO event solver localizes the particle
hit. A particle already outside after launch is lost immediately, rather than
counted as a failed coordinate inversion. No new production setting selects the
old wall rule. The separate single-volume GC model retains its GC-edge contract.

| Original marker | Refined FO loss [ms] | Repaired hybrid [ms] |
| --- | ---: | ---: |
| 62 | 0.922633 | 0.919965 |
| 124 | 0.502884 | 0.503914 |
| 126 | 0.656283 | 0.656637 |
| 129 | 1.988529 | survives at 2.000000 |
| 247 | 0.427111 | 0.428494 |

A direct-FO regression independently checks the hybrid's particle position,
velocity and event time at a wall hit with an interior guiding centre. A second
regression verifies seven exterior and nine interior launch particles against
the raw boundary geometry. Restoring the old outer exclusion makes both fail.
Ledger F8 records the removed multi-volume GC-only absorber and its evidence.

## Marker 129

Sixteen evenly spaced FO launch phases at the same GC start lose the particle
between 1.909902 and 2.055457 ms. Fifteen lose it before the 2 ms cutoff; one
survives that cutoff. This is a launch-phase sensitivity control, not a phase
average of the production hybrid. The latter discards phase between FO episodes.

The default hybrid survives through 2.125 ms. Refining its FO table to
48/96/48 or tightening the Meiss construction ODE tolerance from 1e-2 to 1e-6
does not change that outcome. The final reported particle radius changes by
about 4.1e-5 and 6.0e-7, respectively, from the default 2.879180.

| Microstep control | Steps through 2.125 ms | Final particle radius | Loss |
| --- | ---: | ---: | --- |
| Original | 14,212 | 2.879180 | no |
| approximately h/2 | 28,424 | 2.795599 | no |
| approximately h/4 | 56,780 | 2.759353 | no |
| approximately h/8 | 113,492 | 2.843292 | no |

These controls preserve survival but do **not** establish trajectory convergence.
The requested microstep changes both GC steps and FO call budgets; the internal
adaptive FO error control and gyro cap remain unchanged. The integer macrostep
schedule makes the finer nominal factors approximate. Equal loss counts do not
hide the non-monotone position differences.

Using six entry gyroradii instead of three changes this marker to a loss at
1.980922 ms. Keeping FO active throughout the hybrid pipeline reproduces the
same-grid phase-zero FO reference at 1.988529 ms, with a loss-time difference
below 1e-16 s. The discrepancy depends on the GC/FO decomposition. This evidence
supports the existing switch/phase limitation; it does not establish that every
remaining implementation path is bug-free. Matching 53 by widening the entry
region is not a repair and the production distance remains three gyroradii.

The next scientific change is the coherent action/phase-carrying GC model and
its entry/exit map, followed by joint step, entry-distance and phase controls.
The local first-order switch prototype alone does not supply that model.

## Full-duration checkpoint

All runs use 3.5 keV alphas, 256 matched volume-sampled GC starts and 2 ms.
The latest hybrid replays use 16 local threads and have zero numerical stops,
FO failures and warning holds. Unchanged FO references are reused, with their
original source hashes and run provenance kept separately.

| Field/model | Losses | Paired status differences from refined FO | Trace / end-to-end [s] |
| --- | ---: | ---: | ---: |
| NCSX refined FO | 99 | reference | 440.55 / 458.70 |
| NCSX repaired hybrid | 97 | 8 | 238.37 / 258.19 |
| G3V3L3Fi refined FO | 53 | reference | 109.02 / 128.60 |
| G3V3L3Fi repaired hybrid | 52 | 1 | 52.65 / 67.09 |

G3V3L3Fi's refined hybrid also loses 52, with identical final statuses to its
default hybrid. NCSX changes only marker 50 from loss to survival relative to
the pre-wall-repair hybrid (98 losses). Its eight disagreements remain open;
the wall repair is not a whole-hybrid acceptance. These shared-host accuracy
replays overlapped local verification/build work, so their timings are not a
new controlled speedup measurement. The earlier paired 5.8–7.4x tracing result
retains its original scope.

## Reproduction and verification

`results/loss-diagnosis/` retains all starts and marker outcomes, compact traces
including radial extrema and switch records, per-case configurations and source
hashes, and the native verification record. `results/full-wall-fixed/` retains
the latest six-control snapshot; `results/full-gc-fixed/` preserves the earlier
48-versus-53 result. Raw dense CSVs remain under ignored `build/study/`.

Replot without rerunning orbits:

```sh
python concept/spectre_loss_study/diagnose_report.py \
  concept/spectre_loss_study/results/loss-diagnosis build/study/diagnosis-replot
python concept/spectre_loss_study/full_report.py \
  concept/spectre_loss_study/results/full-wall-fixed build/study/full-wall-replot
```

Both retained-data round trips reproduce their numeric records and metadata.
The integrated affected sweep passed 24/25 checks; its sole failure exposed the
unintended change to single-volume GC. After restricting the repair to the
multi-volume hybrid, the unchanged single-volume test and both wall regressions
pass, 3/3. The other 22 results are reused because their multi-volume equations
are unchanged. The isolated performance-based wall branch passes nine native
checks, followed by the same three checks after the scope correction. No test
starts, energy settings or acceptance bounds were changed. Full-suite repetition
and CI were not awaited. Required domain/senior review still governs promotion.

Plots uploaded 2026-10-15 expire after three days; the sources and numeric data
above remain available for regeneration.

Five-marker diagnosis and launch-phase control:
https://box.sloppy.at/5cbec.png

Latest full-duration loss curves:
https://box.sloppy.at/5bf93.png

Latest measured cost (shared-host accuracy replay):
https://box.sloppy.at/a24e1.png
