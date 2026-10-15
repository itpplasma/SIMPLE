# SPECTRE orbits: plan

Updated 2026-10-15. Production repairs #534–#536 and research PR #537 are on
`main`, followed by the measured cost and residence records below. The frozen
measurement baseline of the earlier loss study remains `fe002a3`. Recovery and
archive history and all retained failed results remain available (see the end).

## Position

Spline-based, event-resolved FO is the SPECTRE reference and a usable
production path. The GC/FO hybrid is an accelerator whose value is measured by
cost against FO at matched loss statistics. Model work is justified only by a
measured cost or accuracy gain on that scale. Hamiltonian, canonical momenta,
phase evolution and particle reconstruction must still belong to one consistent
formulation; a frequency counter added to the fixed-mu GC state does not
provide one. Matching final loss counts alone is insufficient.

## Measured state (NCSX, 3.5 keV alphas, 2 ms, M3 Ultra, 16 threads)

Cost is core-seconds per simulated particle-second.

| Model | Cost | vs Boozer GC default | Source |
| --- | ---: | ---: | --- |
| Boozer GC on the source VMEC, production default | 19.7 | 1 | [boozer-cost](concept/spectre_loss_study/results/boozer-cost/README.md) |
| SPECTRE GC, midpoint `npoiper2=256` / 1024 | 106 / 390 | 5 / 20 | [fo-residence](concept/spectre_loss_study/results/fo-residence/README.md) |
| Hybrid, entry 1 radial gyroradius | 1,591 | 81 | fo-residence |
| Hybrid, entry 3 (default) | 3,862 | 196 | both |
| FO everywhere | 11,629 | 591 | both |

- FO costs about 30x GC on the same field. At entry 3 the hybrid spends 30%
  of particle time and about 90% of compute in FO.
- 63% of FO time is in 35 episodes longer than 1 ms: particles residing in the
  entry/exit band of a sheet. Crossings under 10 us hold about 1% and the outer
  wall 0.1%.
- FO time scales roughly linearly with the entry distance (2.5% at 0.01, 10% at
  1, 30% at 3). At 256 markers, loss curves for entries 0.5–4 are statistically
  indistinguishable from FO (binomial uncertainty about 0.03).
- Midpoint at `npoiper2=256` (the example input) cuts GC cost 3.7x without a
  resolvable loss change but makes the switch chatter at interface 2.
  Euler1 chatters far more and doubles total cost.
- The earlier 42–45% hybrid thread efficiency was the study driver's static
  schedule; production already uses dynamic scheduling, as the driver now does.

Earlier checkpoint (Linux host, refined FO references): NCSX hybrid 97/256
versus FO 99, eight paired status differences; G3V3L3Fi hybrid 52 versus 53.
FO default and refined grids differ by four NCSX statuses. Marker 129's 16 FO
launch phases span loss times 1.91–2.06 ms across the 2 ms cutoff. See the
[R&D report](concept/spectre_loss_study/RND.md).

## Working rule for iteration

Iterate with runs under about a minute (for example 64 markers, 200 us, matched
starts). Reserve 256-marker or larger 2 ms runs for decision gates, announce
them beforehand and run them in the background. Use matched starts and the
same build for any comparison, and normalize cost per particle-second.

## Next steps

1. **Resolve accuracy versus entry distance.** A 1,024-marker comparison of FO,
   entry 1 and entry 3 is running. If entry 1 agrees with FO within
   uncertainty, adopt a lower default entry with the present switch and record
   the gain (about 2.4x). Use bootstrap confidence bands for loss-time curves
   and keep phase-sampled ensembles separate from paired-phase runs.
2. **Stop the switch chatter at the coarser GC step.** Characterize the
   interface-2 episodes at `npoiper2=256` with short runs, then fix the
   hysteresis or landing so that midpoint 256 keeps its 3.7x GC saving without
   extra switching.
3. **Phase-carrying GC switch (GC1+) only if step 1 shows a gap.** Its purpose
   is to permit entry near one gyroradius with FO-level loss agreement. Build it
   on the retained local Darboux completion in `concept/gc_higher_order/` with
   per-volume spline tables; gate on scaled round trips at 1e-12, reconstruction
   orders with temporal error separated, and the loss-curve comparison of step 1.
4. **Longer and larger validation** at the chosen settings: longer durations,
   G3V3L3Fi alongside NCSX, and more markers once the short controls agree.

## Deferred, not rejected

Second-order GC (GC2), phase-carrying encounter maps and cycle compression,
resolved smooth-layer GL2, variational symplectic FO, NCSX equilibria with
2/4/8 volumes and further proof formalization stay retained in their concept
directories. Reopen one when a measured cost or accuracy gap requires it.
Residual FO at sheets (2.5–10% of particle time) is intrinsic; encounter maps
are the only listed route to reduce it further.

## Open concerns

- Marker 129 and NCSX paired differences: retained diagnostics, not tuning
  targets; do not adjust the entry distance to manufacture a matching count.
- Pre-existing midpoint axis-chart concern from the Sonnet review: needs a
  reproducible accepted crossing before a targeted repair
  ([REVIEW.md](concept/spectre_loss_study/REVIEW.md)).
- Independent human domain/senior review of #534–#537 was not performed; CI was
  not awaited.

## Proof and verification obligations

Proof work accompanies model development. It must distinguish source
correspondence, exact algebra, ideal mathematical maps, finite solves/arithmetic
and independent numerical behavior. One layer does not discharge the others.

| Retained evidence | Remaining obligation |
| --- | --- |
| GC1+/GC2 analytic prototypes: reconstruction-order and canonical midpoint/Gauss checks | Real SPECTRE tables, units, chart/axis patches, production integration and comparable-cost evidence |
| Local switch: 132 exact identities, regeneration, real-field Hessian/pullback checks, 10-case/seven-epsilon reconstruction study | Generic field hypotheses, extended canonical action/phase evolution and complete production entry/exit |
| Three GC source repairs with independent closed-Hamiltonian RED/GREEN witnesses | Uniform finite-Newton, field/domain and source-to-theorem bounds |
| Eight GC methods pass retained one-/100-step stellarator map checks | Long energy/domain acceptance; whole implemented map |
| FO Gauss2/entry-frozen/Lie and GC tableau Lean results with standard axioms | State-selected parameters, landing/exit time, chart changes, terminated Newton and whole hybrid |
| Provider-jet baseline: 1149 exact identities | Complete the recorded mutation/source-wiring and domain obligations |

The long-energy acceptance remains open. QH Boozer passes all eight methods
at both starts for 100,000 steps. NCSX trapped starts reach the chart boundary,
including unresolved status 5 near it; do not relabel this physical loss.
NCSX passing Gauss2 fails the original coarse-step envelope gate and passes
an equal-time half-step control. Preserve both outcomes and original bounds.
The known RED GC energy acceptance remains outside the default registered suite,
with its sources, data and registration patch retained explicitly.

Baseline verification is retained, not presented as a fresh full-suite pass:
the consolidated Fast sweep passed 98/99, followed by the repaired diagnostic's
passing rerun. The wall affected sweep passed 24/25; its single-volume scope
failure was repaired and rerun with both wall regressions, 3/3. The isolated
wall branch passed nine native checks plus that three-case scope rerun.
The old wall exclusion makes both wall witnesses fail. See the R&D report
and `results/loss-diagnosis/verification.json` for exact revisions and limits.

Use native Make/CMake authorities for affected behavior. Gremlin previously
returned `Gremlin state provider error: 2`; this is not a scientific result.
Do not relax thresholds, erase failures or add tests merely to reject obsolete
parameters. Do not wait for CI or repeat unaffected checks without a reason.

## Rejected paths and retained provenance

Ledger [F1–F8](DOC/spectre-rejected-interface-models.md) records fixed-mu sheet
resets, phase-discarding reversible/symplectic claims, refreshed frozen patches,
unresolved sharp stepping, energy-kick landing, Boris as a canonical guarantee,
the polar-axis spline defect and the multi-volume GC-only wall absorber.
Keep their failed implementations/settings retired. Higher-order smooth-volume
GC, carried-phase switching and resolved FO are not rejected paths.

| Recovered work | Authoritative location |
| --- | --- |
| Physical throughput/loss controls and current diagnosis | `concept/spectre_loss_study/` |
| GC1+/GC2 derivation, canonical completion and measured analytic results | `concept/gc_higher_order/` |
| Switch/provider/inexact-Newton/GC2 prerequisite work | `concept/fo_gc_derivations/`; local switch in `D_canonical_switch/` |
| GC source audit, map tests and open energy acceptance | `concept/gc_symplectic_proof/` |
| FO proofs and test-only Gauss2 | `concept/fo_symplectic_proof/` |
| FLR/planar encounter derivations | `concept/flr_options/` |
| Surviving candidate comparison and sharp-sheet controls | `concept/interface_candidates/`, `concept/symplectic_sheet/` |
| Superseded W4 implementation | `concept/fo_logical_spline/README.md` links its Git archive; unused code/patch removed during integration cleanup |

The [previous full recovery ledger](https://github.com/itpplasma/SIMPLE/blob/fa913186ab29a0cfba3175d6a6b2c7a515f77a1f/PLAN.md)
retains historical commands, timings, worker commits and verification details.
Recovery is `recovery/spectre-session-20261014` at `082bea4`; the pre-cleanup
archive is `archive/spectre-session-20261014-pre-cleanup` at `87aa81f`.
The research integration preserves recovery (`8afaecb`), GC repair (`0993478`),
switch work (`97fde90`, `060d17f`, `4d963cc`), wall repair (`535bc9e`, `5e88fe5`)
and ledger F8 (`c957615`). Separate audit/switch branches retain provenance.
Do not resume old agents or merge old worker tips over the cleaned sources.
Generic Lie/code-generation support stays in FortSym, including the retained
`feat/lie-pullback` dependency work; application physics stays in SIMPLE.

All new measured runs were local. The last retained read-only cluster snapshot
is 2026-10-14 20:18 CEST, not live monitoring. No protected-host work is
authorized: do not run/schedule/install/resume on faepop* or faepcr* without
an explicit request naming the host.

Keep regular commits/pushes and update this plan at completed decision gates.
Publish reproducible plots through slopbox and record expiring links in the
R&D report. The external R&D destination remains unspecified; no chat post
is claimed. Scratch data/build products stay outside product diffs.
