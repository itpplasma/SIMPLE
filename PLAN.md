# SPECTRE orbits: model comparison and execution plan

Updated 2026-10-15. Integration target: `main`; the research branch is being
retired after integration. Production PRs #534–#536 are merged, through
`27e963c`. The frozen measurement baseline remains `fe002a3` (including #533);
it is not the current main revision. Recovery/archive history and all retained
failed results remain available.

## Strategy and active choices

Keep the verified source repairs, particle-wall repair and cheap spline FO
provider. They improve the reference implementation and comparison conditions;
they do not select the final hybrid model. The present phase-discarding hybrid
remains a measured baseline with unresolved model and convergence errors.

Compare phase-complete first-order GC (`1+`), second-order GC (`GC2`) and FO.
Higher-order GC is an active comparison track. Its production integration
depends on real-stellarator accuracy and cost evidence. Build the shared
canonical action/phase infrastructure for both GC orders before choosing one.

| Approach | Role and next decision | Status |
| --- | --- | --- |
| Spline-based, event-resolved FO | Reference; resolve sharp-sheet and particle-wall encounters | Implemented and measured; no new symplectic guarantee |
| Present GC/FO hybrid | Baseline for cost, losses and switch error | Phase is discarded; trajectory convergence remains unresolved |
| Phase-complete first-order GC, `1+` | Essential GC baseline with consistent phase evolution and reconstruction | Local prototype exists; real SPECTRE implementation required |
| Second-order GC, `GC2` | Active contender using the same canonical infrastructure | Local prototype exists; stellarator benefit and production cost unmeasured |
| Encounter maps and cycle compression | Conditional acceleration if resolved FO encounters remain expensive | Exact planar controls exist; curved-field extension unvalidated |
| Resolved smooth-layer GL2 and variational FO | Secondary physics/structure controls; later alternatives | Concept-only; no production option or complete hybrid proof |

Higher order addresses smooth-volume dynamics. It supplies no GC crossing rule
for a sharp sheet. Keep resolved FO at sheets while comparing GC orders.
Hamiltonian, canonical momenta, phase evolution and particle reconstruction
must belong to one consistent formulation. A frequency counter added to the
current fixed-mu GC state does not provide that formulation.

Choose by physical error versus total runtime, including table construction,
memory and FO residence. Matching final loss counts alone is insufficient.
Retain the eventual bounded-energy and symplectic objectives, with separate
evidence for continuous physics, the numerical map and its implementation.

## Current implementation and measured checkpoint

| Production change | Review branch / revision | State |
| --- | --- | --- |
| FO logical splines, endpoint reuse and regular axis chart | [#534](https://github.com/itpplasma/SIMPLE/pull/534), `perf/spectre-fo-logical-spline`, `47d3fb0` | Merged as `934d964`; owner-directed promotion |
| Three GC final-Newton extrapolation repairs | [#535](https://github.com/itpplasma/SIMPLE/pull/535), `fix/gc-field-extrapolation`, `6970eff` | Merged as `71ad9ae`; owner-directed promotion |
| Multi-volume particle-wall loss repair | [#536](https://github.com/itpplasma/SIMPLE/pull/536), `fix/spectre-particle-edge-loss`, `914544f` | Merged as `27e963c`; owner-directed promotion |

All three changes are integrated on main at Chris's explicit merge request.
Independent human domain/senior review under `docs/review-policy.md` was not
performed; the owner-directed promotion is recorded in the merge descriptions.
Three requested read-only Claude Sonnet reviews are advisory, not that approval.
Their [findings and adjudication](concept/spectre_loss_study/REVIEW.md) retain
one unconfirmed pre-existing midpoint axis-chart concern. CI is not awaited.

Latest physical controls use 3.5 keV alphas, 256 volume-sampled GC starts,
2 ms and 16 local OpenMP threads. The hybrid replays have zero numerical stops,
FO failures and warning holds. The unchanged FO references retain their
original case-specific provenance.

| Field/model | Losses / 256 | Paired status differences from refined FO | Trace / total [s] |
| --- | ---: | ---: | ---: |
| NCSX refined FO | 99 | reference | 440.55 / 458.70 |
| NCSX repaired hybrid | 97 | 8 | 238.37 / 258.19 |
| G3V3L3Fi refined FO | 53 | reference | 109.02 / 128.60 |
| G3V3L3Fi repaired hybrid | 52 | 1 | 52.65 / 67.09 |

G3V3L3Fi's refined hybrid also loses 52 with identical final statuses.
NCSX default/refined FO disagree on four statuses and some loss times by
1.08 ms. G3V3L3Fi FO grids agree on every status, with maximum loss-time
difference 2.14e-12 s. These finite ensembles do not establish population
agreement or whole-hybrid accuracy.

The latest hybrid timings are shared-host accuracy replays that overlapped
other checks/builds. The separate 24-case, 32-marker, 200 us paired study
measured 5.8–7.4x tracing speedups against main, excluding setup. All health
counters were zero and paired loss statuses matched main in that study.
NCSX sheet FO separation from main decreased from 26.3 cm at the default grid
to 0.623 cm at the refined grid. No path-identity or >=10x claim is established.
The earlier 2–4 h estimate is superseded for this measured workload.

## What the 48-versus-53 discrepancy established

The multi-volume hybrid excluded the outer edge from FO proximity entry.
It absorbed the guiding centre there while the FO reference absorbed the
particle. Applying the existing FO region to the outer edge recovers markers
62, 124, 126 and 247: 48 becomes 52 on both spatial grids. An initially exterior
launch particle is now an immediate physical loss, not an inversion failure.
The separate single-volume GC model retains its GC-edge contract.

Marker 129 remains unresolved. Sixteen FO launch phases at its fixed GC start
span loss times 1.910–2.055 ms across the 2 ms cutoff. Finer FO tables and tighter
Meiss construction leave the hybrid surviving through 2.125 ms. Microstep
refinement preserves survival but does not establish trajectory convergence.
Changing the entry distance changes its outcome; FO throughout the hybrid
pipeline reproduces the same-grid phase-zero FO loss time to about 1e-16 s.
Do not tune the entry distance to manufacture 53 or accept the residual as
validated physical disagreement. NCSX's eight paired disagreements also remain
open; the wall repair changes only marker 50 from loss to survival there.

The [R&D report](concept/spectre_loss_study/RND.md) contains the marker controls,
verification record, reproduction commands and plots. Results are retained in
`concept/spectre_loss_study/results/full-wall-fixed/` and `results/loss-diagnosis/`.
The earlier `full/` and `full-gc-fixed/` snapshots remain intact.

## Ordered execution and decision gates

### 1. Preserve and promote the verified baseline

Production PRs #534/#535/#536 are merged. Integrate the retained research,
then retire the completed integration branches. Keep new model development
separate. Preserve the add/remove cleanup history,
independent axis-field oracle, GC source regressions and physical wall witnesses.
Use the repaired implementation for subsequent comparisons without treating
its current hybrid as the selected solution.

### 2. Build the shared canonical GC infrastructure

Use the retained local Darboux completion in `concept/gc_higher_order/`.
Carry gyroaction and phase, the radial phase primitive K and angular connection
tables C_a. Audit charge/orientation, normalization, gyrogauge patches and
matching to SIMPLE's canonical variables. Entry, evolution and exit must use
the same representation. Phase-complete `1+` needs the consistent H2 frequency
correction; `1+` and `GC2` therefore share substantial infrastructure.

Construct connection and H2 coefficient fields offline and use per-volume
splines in the tracing loop. Audit required derivatives before selecting the
table representation. The GC2 prerequisites note's third-/fourth-derivative
requirements are analysis, not completed machine proofs. Bound and measure
canonical momentum inversion, axis/chart handling, table memory and setup cost.
Do not restore raw spectral geometry or physical-coordinate inversion at
every orbit stage. Do not add production selectors before a viable model exists.

Gate: source-linked identities and independent field/derivative oracles,
consistent physical units and local chart handling, scaled round trips at
1e-12, and the expected reconstruction orders with temporal error separated.
Local prototype checks do not discharge the real-field obligations.

### 3. Compare 1+, GC2 and FO inside smooth stellarator volumes

Implement both GC orders on the shared provider and test NCSX/G3V3L3Fi volume
interiors before introducing sheet transitions. Start each model from the same
physical particle position and velocity, using its own consistent inverse map.
Equal numerical GC coordinates across orders are not matched particle states.

Use an epsilon scan and step/grid refinement to distinguish GC truncation,
phase error, field interpolation and time discretization. Check reconstruction,
phase/action, invariants and nontrivial energy envelopes against tighter FO.
Verify field-specific expected orders; do not infer a universal order from
selected analytic examples. Measure total cost and memory at comparable error.

Gate: real-field evidence that each implementation has its claimed accuracy
and a usable cost. GC2 remains an active contender even while its production
integration is deferred. If its extra order brings no measured benefit at the
target parameters, retain the result and prefer the cheaper adequate model.

### 4. Couple each viable GC order to the same resolved FO encounters

Replace entry and exit together with consistent particle/GC maps and carried
phase. Retain event-resolved FO at sharp sheets and the same particle absorbing
outer edge for every multi-volume candidate. Resolve owner/gauge changes and
phase-chart transitions coherently. A changed mu at a sharp sheet is physical;
a fixed-mu reset is not an alternative crossing model.

Gate: independent reconstruction through finite FO excursions, expected
epsilon scaling, phase/gauge consistency, clean event/failure classification
and step/entry-distance convergence. Marker 129 is a retained diagnostic,
not a special case to tune or an acceptance criterion by itself.

### 5. Run the matched loss and cost comparison, then select a model

Compare FO everywhere, the present hybrid control, `1+` hybrid and `GC2` hybrid
on NCSX/G3V3L3Fi. Begin with matched physical particle ensembles, at least
256 markers and 2 ms. Keep deterministic paired phases and sampled phase
ensembles distinct. Use the same physical boundary and equilibrium bytes.

First compare the default settings with one step/grid/tolerance refinement.
Then cover entry distances 1.5/2/3/4 radial gyroradii and FO tolerances
1e-7/1e-9/1e-11, with FO-everywhere reference refinement. Declare accuracy
criteria before looking at the results. Extend duration/ensemble size after
the initial controls converge; a 2 ms result is not long-time validation.

Publish loss-time curves, paired differences and uncertainty, reconstructed
orbit/invariant errors, FO residence/excursion distributions, memory, setup
and end-to-end cost. Keep numerical stops, recovered failures and warning
holds separate from physical losses. Control host load for timing comparisons.
Select on accuracy versus total cost; extra formal order or matching counts
does not win automatically. No complete model-selection matrix exists yet.

### 6. Pursue conditional acceleration and equilibrium checks

If resolved FO encounters still dominate after the comparison, evaluate
phase-carrying encounter maps/cycle compression against the same reference.
Exact planar controls do not justify indefinite frozen curved-field patches.
Stochastic phase models require a separate ensemble validation and do not
supply deterministic reversibility or symplecticity.

Obtain valid NCSX equilibria with 2/4/8 volumes and compare matched physical
starts to assess equilibrium-volume sensitivity. Clipping to one volume creates
an artificial absorber and is not an equivalent wall-loss comparison.
Keep resolved smooth-layer GL2 as a distinct physical-limit control.
Revisit variational symplectic FO if the accuracy, energy and cost evidence
requires it; its crossing-time solve and whole-map proof remain unfinished.

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
