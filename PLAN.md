# SPECTRE orbit programme: recovered status and plan

Updated 2026-10-14. Active study: `research/spectre-loss-study-20261014`.
Recovery branch: `recovery/spectre-session-20261014`.
Production baseline: `origin/main` at `fe002a3`, including cleanup PR #533.
Measurement baseline remains pinned to this revision even if main advances.
This consolidates the stopped Claude controller and its GPT Sol, strategy,
Opus and nested derivation workers. It records unfinished research rather
than claiming the scientific programme is complete. No previous worker is
running or being resumed by this recovery.

## Current measured checkpoint

The repaired 24-case controls are complete: 32 matched markers, 200 us,
16 local threads, NCSX/G3V3L3Fi, sheet/volume starts, FO/current hybrid,
and main/default/refined fields. Every case completed with zero numerical
stops, FO failures or warning-held intervals. All paired loss statuses agree
with main. Retained starts, sparse traces and source/dependency hashes are in
`concept/spectre_loss_study/results/long/{summary.json,traces.npz}`; reproduce
plots with `analyze.py` using that directory as input.

| Case | FO tracing speedup | Hybrid tracing speedup | FO / hybrid losses |
| --- | ---: | ---: | ---: |
| NCSX sheet | 7.43x | 7.44x | 4/32 / 3/32 |
| NCSX volume | 6.63x | 7.43x | 9/32 / 9/32 |
| G3V3L3Fi sheet | 5.79x | 6.47x | 0/32 / 0/32 |
| G3V3L3Fi volume | 5.87x | 6.11x | 2/32 / 2/32 |

These are tracing speedups on a shared host, excluding setup. They do not
establish trajectory identity: the largest NCSX sheet FO separation from
main is 26.3 cm at the default grid and 0.623 cm at the refined grid after
200 us. Hybrid sheet runs have up to four sampled GC/FO representation
mismatches; paired positions exclude those records. The repaired near-axis
marker completes, but its raw-main trajectory differs from both converging
spline grids. Independent axis-limit tests, rather than exact agreement with
a singular raw pole evaluation, check the field there. Aggregate loss
agreement alone is insufficient for scientific promotion.

The six 256-marker, 2 ms physical-volume controls are complete, each on
16 bound local OpenMP threads. Both hybrids were replayed after the GC repair
with zero stops/failures/holds. The table uses those latest hybrid rows and
the unchanged original FO references (FO throughout does not execute the
repaired GC stepper). Retained starts, sparse traces, case-specific provenance
and paired statistics are in `concept/spectre_loss_study/results/full-gc-fixed/`.
The original six-run snapshot remains in `results/full/`. Reproduce plots with
`full_report.py RESULTS_DIR OUTPUT_DIR`; the replay is an accuracy check on a
shared host, rather than a fresh controlled throughput benchmark.

| Field / model | Tracing [s] | End-to-end [s] | Losses / 256 | FO residence |
| --- | ---: | ---: | ---: | ---: |
| NCSX FO default | 436.59 | 451.36 | 101 | 100% |
| NCSX FO refined | 440.55 | 458.70 | 99 | 100% |
| NCSX current hybrid | 228.82 | 249.20 | 98 | 29.9% |
| G3V3L3Fi FO default | 105.17 | 117.70 | 53 | 100% |
| G3V3L3Fi FO refined | 109.02 | 128.60 | 53 | 100% |
| G3V3L3Fi current hybrid | 49.19 | 65.98 | 48 | 38.8% |

Relative to refined FO, NCSX default FO has four paired status disagreements
(net +0.781 percentage points, paired SE 0.781 pp); the repaired hybrid has seven
(net -0.391 pp, SE 1.035 pp). Both-lost NCSX markers can differ in loss time by
1.08 ms (default FO) / 1.24 ms (hybrid), so equal final counts are not pathwise
convergence. G3V3L3Fi default/refined FO have no status disagreements and a
maximum paired loss-time difference of 2.14e-12 s. Zero disagreements in 256
markers still permits an unobserved disagreement rate up to 1.16% at one-sided
95% confidence under an independent-marker model. G3V3L3Fi hybrid loses five
fewer markers (net -1.953 pp, paired SE 0.867 pp), all five status differences
in that direction. This discrepancy remains to be explained by the consistent
switch/entry and boundary study; it is not removed by the FO spatial grid.

Maximum sampled energy-shell errors are 5.97e-5 for NCSX FO, 5.00e-4 for its
hybrid, and 1.30–1.34e-5 for the G3V3L3Fi controls. These are measured errors,
not a bounded-energy or whole-map symplectic proof. The fixed production entry
phase and GC-centre/FO-particle absorbing-edge conventions still apply.

These direct costs replace the earlier hours-scale estimate for this workload;
thread counts, host load and implementation differed in that estimate. The
complete six-control batch took about 24 minutes sequentially including native
launch/setup/output. No direct 256-marker/2 ms raw-main run was measured here.
No protected-host job was started, resumed or cancelled.

Published plots (uploaded 2026-10-14; links expire after three days):

- Latest FO-reference / GC-fixed hybrid end-to-end cost: https://box.sloppy.at/2c9db.png
- Latest loss curves and pointwise uncertainty: https://box.sloppy.at/8af3d.png
- Paired main/spline tracing cost at 200 us: https://box.sloppy.at/17ce5.png
- GC energy/domain limits and equal-time refinement: https://box.sloppy.at/a4710.png
- Local first-order switch accuracy / pullback scaling: https://box.sloppy.at/667cc.png

The GC source/proof audit is retained on
`research/gc-source-audit-20261014` through `d9df7a8`. The independent closed-form
Hamiltonian witness reproduces three production bugs (RED `810f1e5`): an extra
midpoint-radius Taylor increment in both midpoint branches, a wrong Gauss
r/phi Hessian index, and omitted Gauss parallel-speed p_phi increment. The
fix is `6a4fc33`; each individual regression is rejected when reintroduced.
The dedicated main-based branch `fix/gc-field-extrapolation` contains only
those source repairs, the independent regression and affected documentation;
it is pushed as PR #535, commit `6970eff`. Native main-based checks pass
85/85 in 478.83 s; after fixing local interpreter discovery, the rebuilt
Python API suite passes all 17 cases in 99.52 s. Other backend-dependent Python
suites were not rerun. The repair is integrated here as `0993478`; affected
SPECTRE checks pass, 22/22 in 137.19 s. Both 256-marker hybrid controls
were replayed against retained FO references, giving 98 NCSX / 48 G3V3L3Fi
losses with no stops/failures/holds. The original NCSX hybrid lost 99; one
additional paired final-status difference is retained rather than suppressed.

All eight methods pass the one-step/100-step NCSX canonical-flux/Boozer and QH
Boozer map checks with unchanged bounds (1e-7/1e-5) and original spacings
(1e-4/5e-5). The original fourth-order finite-difference stencil failed two
NCSX 100-step checks before and after the source repair; an eight-point stencil
improves oracle resolution and passes. Both outcomes remain retained.

The long energy acceptance is **not closed**. QH Boozer passes all methods at
both starts for 100,000 steps. NCSX trapped starts reach the chart boundary;
some landing solves return status 5 near s=1, which is not relabelled physical
loss. NCSX Boozer passing Gauss2 fails the original coarse-step envelope gate;
halving the step at equal elapsed interval passes and reduces extrema from
roughly 1.6e-4 to 2.6e-7. Original failures, source hashes, samples and plots are
retained under the audit branch's `concept/gc_symplectic_proof/`. No thresholds
or starts were changed to obtain acceptance. Ideal fixed-step tableau theorems
survive; finite-Newton/field/domain and whole-hybrid proof bridges remain open.


The consistent-switch concept continues on
`research/canonical-switch-20261014` (`7626cd8` for the portable build and real
Hessians). Its original zero magnetic Hessians fail independent Cartesian
checks on both stellarators. Correct logical second jets plus geometry
chain-rule terms pass all six volume samples (max scaled error 2.62e-6,
fixed acceptance bound 2e-4). The build/runner now use native Make/CTest.
New independent two-form pullback and analytic-Jacobian checks pass on analytic
fields and all three volumes of each stellarator: selector 1 defects scale
with epsilon squared; selector 0 remains first order. These are numerical
checks, not exact emitted-map identities or production integration. The full
10-case, seven-epsilon reconstruction study is complete: finest-four mu-defect
orders are 1.995–2.001, worst scaled round trips below 4e-13 and temporal
refinement changes below 0.033%. The source-linked regeneration check reproduces
the committed kernel exactly. The repaired checker proves 132/132 identities
without probes, for its stated local field family and generic Lie-series
claims; this is not an exact general nonlinear SPECTRE-field identity. These
results are pushed as `c3a3f31`. Production still uses its existing switch and
phase. Actual extended canonical GC action/phase, connection tables and
production entry/exit remain to be implemented and validated.

The recovered branch has now been merged into this checkout (`8afaecb`), and
the completed switch work is integrated as `97fde90`, `060d17f`, `4d963cc`.
Latest GC audit sources, results and exact native registration patch are also
collected here; their known RED energy acceptance stays outside the default
suite. Historical Lean caches are preserved and ignored. Production changes
remain reviewable as the two focused PRs #534 and #535; this research branch
also contains the recovered concept inventory and all retained new evidence.
A fresh native Fast suite is running on the consolidated checkout.

## Goal and constraints

Physically faithful long stellarator orbits, without artificial loss or
confinement, with bounded oscillating energy error. The eventual symplectic
claim must cover the actual implementation through source-tied generation,
exact symbolic identities, mathematical theorems and executable oracles.
Continuous physics and numerical discretization must be validated separately.

NCSX and G3V3L3Fi are acceptance cases. Tokamak/axisymmetric examples provide
smoke checks and derivation controls only. Retain genuine failures and their
original thresholds. Rejected SPECTRE modes remain in history, without shims
or compatibility switches. Higher-order GC remains a standalone concept until
its field requirements and reconstruction are established in production.

The previous controller narrowed the immediate proof target to smooth-volume
GC and kept RK45 at sheets. This is a practical intermediate decision, not
fulfilment of the requested eventual symplectic hybrid. Neither a proof of the
whole current hybrid nor a general impossibility theorem has been established.
The last instruction in that session was to stop agents, record state and push.
This recovery publishes their work and an actionable continuation plan.

## What already shipped

- PR #533: collect the comparison in `4a1f842`, then remove demonstrated failed
  paths in `3a6fd37`, with the exact provenance pinned in the follow-up ledger.
  The merge preserves add/remove history. Removed executable candidates:
  fixed-mu C1L0/L1, phase-discarding C2fixed/pphi, curved-field frozen-patch C4,
  Boris C5 as a symplectic candidate, and the old unresolved branch/energy-kick
  landing methods. Their implementation/settings/scans and dedicated programs
  are gone. No new checks or tests reject outdated parameters. The scope and
  evidence of each rejection are in the
  [implementation ledger](DOC/spectre-rejected-interface-models.md#implementation-removal-ledger).

- PR #526: Fortran configuration-file API with a thin Python layer, matching
  executable initialization, API/docs cleanup and test pruning. The earlier
  API session was superseded by this merged work; do not restart its cleanup.
- PR #527: SPECTRE drift sign, relocated transmission ownership and fixture
  corrections. PR #529: non-folding NCSX free-boundary fixture.
- PR #531: remove rejected fixed-mu GC interface maps. This is included in
  the ancestry of #530 and the current main baseline.
- PR #530: symplectic GC in smooth volumes, RK45 full orbit near interior
  sheets, event resolution, gyromotion bounds, gauge-consistent exit and
  multi-volume Newton domain. SPECTRE Boris/layer controls and GC interface
  maps were removed. Multi-volume levels 0/1 and GC-only RK45 are rejected.
- `cbb4a7a`: trapping-parameter documentation. The separate session closed
  issue #532 with a checked reconstruction of confined fractions from NetCDF;
  the correct endpoint comparison is `times_lost >= t`.

The rejection record and its symbolic verifier already live in
`DOC/spectre-rejected-interface-models.md` and
`concept/symplectic_sheet/proofs/verify_rejected_models.py`. Old symplectic
crossing/reset constructions do not regain production status in this recovery.

## Scientific interpretation, evidence and limits

**Sharp sheet reference.** Integrate the Lorentz equation through the bounded
field discontinuity, resolving events; particle position and velocity are
continuous. A gyrophase-dependent decomposition into parallel speed and mu
can jump. Conserving mu across a sheet narrower than a gyroradius imposes a
false mirror barrier. Smooth-layer/GC and sharp-sheet/FO are different physical
limits; an integrator comparison cannot decide which equilibrium model is
appropriate for an unresolved physical current sheet.

**NCSX losses.** Historical strong-interface runs lost 39/64 with full orbit
everywhere and about 39–41/64 with the hybrid, versus 7/64 with the removed
fixed-mu GC map. The final stopped-session RK45 repeat reported 39/64 at
rho=1.95 (2359.1 s on a loaded host). This supports the sharp-sheet explanation
of the excess relative to GC, not a universal accuracy guarantee. The old
symplectic level-1 row also contained FO recovery excursions and is not a
clean GC baseline. Raw rerun directories under the deleted jobs scratch area
are unavailable; retained concept data and transcript summaries are identified
as historical evidence. The 64 interface-seeded markers do not establish
volume-averaged throughput or convergence with the number of equilibrium volumes.

**Switch diagnosis revised by W1.** The strategy estimated an O(epsilon) kick
for a zero-step round trip and inferred long-time random-walk drift. The actual
diagnostic measured a zero-step rms relative mu defect about 2.4e-5 with order
2, but an exit-phase-dependent amplitude about 6.9e-3 with order 1 (reference
rho_L/R0 about 3.4e-3). Thus immediate round-trip error and finite-excursion
phase error are distinct. The test asserts finite completion, not these fitted
orders; its log is evidence, not a scaling acceptance gate. Dominance of switch
noise over long times remains unproved. The proposed 20 ms GC comparison made
rho=1 absorbing by forcing one volume: 50/64 reached it. That is not comparable
to physical hybrid wall loss. The hybrid long run has no recovered final result.

**Gauss2 decision.** A recorded comparison was 811 s versus 109 s for RK45,
with roughly 12 times as many steps. Refined starts removed retained Newton
stops and showed fourth-order reduction of the marker-20 energy envelope;
these are useful concept results. Frozen entry parameters give a fixed-family
symplectic theorem, but differentiating the state-selected entry parameters,
exit time, GC conversion and re-entry adds terms outside that theorem.
Gauss2 therefore remains concept-only, without a namelist option.

**Higher-order GC.** The phase-complete first-order and second-order prototypes
have retained finite-time reconstruction and canonical midpoint/Gauss tests.
They preserve the gyrogauge connection using a radial primitive and angular
connection tables, rather than changing H alone. This establishes a local
canonical representation of a truncated model, not an arbitrary-order method,
not production stellarator validation and not a sheet-crossing rule.

## Collected work and provenance

| Workstream / source | Collected location | Status at recovery |
| --- | --- | --- |
| GPT FLR research worker | `concept/flr_options/` | Original 23-page note source and algebra/local-map/citation scripts reconstructed from the saved GPT session. Exact planar sheet maps and cycle compression; curved-sheet and production claims remain open. |
| GPT candidate comparison, `345bfa8`, SIMPLE-hybrid | `concept/interface_candidates/` | Historical marker data/manifests remain evidence. Main now retains only FO, C3 and resolved GL2 executable paths; rejected implementations/settings are removed by #533. C3 has no established speed gain or complete hybrid symplectic proof. |
| GPT GC orders + canonical completion, SIMPLE-gcho `c9f23b1`, `957d5a8`, `598f041`; Opus rebase/trim, SIMPLE-gclie `24c625a` | `concept/gc_higher_order/` | Final rebased sources and compact results integrated. Kernels generated in build with pinned output hashes. Production untouched. |
| GPT proof, SIMPLE-sproof/SIMPLE-derive `a27d1ec` | `concept/fo_symplectic_proof/` | Original Lean project, pinned predecessor, symbolic results and measurements retained; original explanatory text in `proof-foundation.md`. |
| GPT Gauss2 task I and Opus demotion, SIMPLE-fo-gauss2 `12aaf7f`, `286dc20`, `d956968`, `3ddc30b` plus stopped edits | Same FO concept directory; registered test-only library and probes | Final concept-only layout integrated. Production files matched main after demotion. Original radial hybrid test was removed by that demotion; its prior form remains in worker history. |
| Opus W1, SIMPLE-w1 `9b28f29` | `test/tests/field_can/test_spectre_switch_roundtrip.f90` | Diagnostic integrated and registered. Long-run comparison remains incomplete/non-equivalent. |
| Opus W4, SIMPLE-w4 uncommitted | `concept/fo_logical_spline/` | Generator, emitted RHS, spline provider and benchmark source retained. All tracked production edits in `pending-integration.patch`; not enabled. The historical prototype is superseded by the validated PR #534 implementation; preserve the draft for provenance, not promotion. |
| Opus GC proof, SIMPLE-gcproof uncommitted | `concept/gc_symplectic_proof/` | Lean/ideal-tableau results retained; numerical oracle and CMake delta preserved as drafts outside default tests. Stellarator L5 and source-to-theorem bridge unfinished. |
| Opus derivations, SIMPLE-derive uncommitted | `concept/fo_gc_derivations/`; two additional FO Lean files | All nonignored A–E sources and available results recovered. See below. |
| Cleanup worker, SIMPLE-remove `d2f5b44`, SIMPLE-mainclean `fc83a75` | Existing rejection docs/verifier and merged history | Already integrated upstream; do not reapply removed maps. |
| SIMPLE-sympl `c91970f`; task-i-pre-rebase `89babb3` | Original history | Superseded/rejected implementations, not active candidates. |
| SIMPLE-dependency-versions `aab1b38` | Merged PR #523 history | Shared dependency update already shipped. |

The original worktrees are retained locally, including their staging state, as
read-only recovery sources. Do not resume or merge their old tips over the
cleanup: the authoritative continuation is this rebased recovery branch.
The pre-cleanup integrated history is also published at
`archive/spectre-session-20261014-pre-cleanup` (`87aa81f`), including the original
worker commits. The first-parent recovery changes were replayed onto main;
ancestry-only experimental production branches were not replayed.
FortSym's consumer-required Lie/DAG/array-CSE work is `de7074c`, `2894087`,
`1372fbd` on `lazy-fortran/fortsym:feat/lie-pullback`. It stays in its owning
library and is published with this recovery; the SIMPLE concepts require that
revision rather than assuming upstream main already contains those APIs.
Build directories, Lean caches, virtual environments and temporary binaries are
excluded. The deleted `~/.claude/jobs/ac5435ab/tmp` area was reconstructed where
possible from the controller/agent JSONL records, including the strategy and
last progress updates. Personal transcripts and raw process logs are not product
files. This plan carries the durable decisions and recovery limits.

## Proof and derivation ledger

The proof chain has five separate obligations: L1 source correspondence,
L2 exact symbolic identities, L3 exact mathematical map, L4 finite solve/
arithmetic/domain bounds, L5 independent numerical behaviour. Completing one
layer does not discharge the others.

| Item | Result retained | Remaining obligation |
| --- | --- | --- |
| FO L2/L3 | Scalar-generated H and monitor; exact symbolic/readback checks; `Gauss2.lean` exact-stage theorem | Provider/domain hypotheses, regular stage branch, floating arithmetic and whole hybrid. |
| A: entry-frozen parameters | `EntryFrozen.lean`: composition, leaf/Poisson and selection-defect results with standard axioms | Chain-rule/free-particle/switch examples remain in the explicitly unchecked draft; state-selected full map not certified. |
| B: provider jets | Baseline 1149 exact identities; copied-source anchors and hashes | Final mutation did not apply; mutation outputs left one failed identity in the original JSON. Complete harness and W4 wiring/domain checks; baseline is rerun separately. |
| C: inexact Newton | Written derivative/defect bounds, audited full-Newton stopping rules | Constants, domain bounds and machine/executable verification. A small accepted residual alone does not bound the derivative of the terminated map; previous rejected residual matters. |
| D: canonical switch | First-order transform plus Newton inverse; seven-case retained study including three NCSX volumes | Portable dependency build, exact emitted-map identity checks, phase/frame and normalization audit, production wiring and long-run verification. |
| D study | First-order mu-defect last-four exponents about 1.995–2.001 on NCSX, versus order 1 controls; round-trip state defect about 2.1e-13 | Absolute chart round trip about 1.0e-11 (relative 6.5e-14): define scaled acceptance, do not claim unqualified absolute 1e-12 for every representation. |
| Generic Lie truncation | `LieTruncation.lean` generic polynomial/remainder statements | Establish its hypotheses for the concrete generated switch, including field derivatives and frame. |
| E: GC2 prerequisites | Written connection, H2 and derivative-order analysis | Charge/orientation/units checks and exact derivation. Fourth magnetic derivatives for H2 Hessian and connection/metric table requirements are analysis, not completed machine proofs. |
| All GC variants | Euler pair, midpoint/Gauss1, Gauss2–4, Lobatto3 exact PRK theorems; 22 standard-axiom results; coefficient audit | Source tie, field derivatives, radial inverse, extrapolation, terminated Newton and landing selections; stellarator map/energy oracle. |

At the stop, the midpoint increment, Gauss Hessian index and omitted vpar
term were source-audit leads. All three have now been independently reproduced
and repaired as recorded above. A coefficient theorem alone still does not
prove the remaining finite-solve/field/domain code paths.

## Performance and run status, 2026-10-14 20:18 CEST

Read-only Slurm queries found no jobs under the account on `acluster` and no
SIMPLE/SPECTRE loss jobs on `scluster`. The latter still runs the separate ITER
GVEC, VMEC, MEPHIT and NTV campaign; its launch commands/work directories were
checked. No local SIMPLE loss tracer was running. No jobs were started or
cancelled. Recovered Claude/Codex command logs contain no remote submissions
for this SPECTRE study. These are status observations, not a monitoring service.

The old eight-thread, 64-marker, 2 ms hybrid cases took 1848–3343 s on a loaded
host (`concept/symplectic_sheet/README.md`). Linear marker-count scaling gives
2.05–3.71 h for 256 markers with the same launch distribution and settings.
The controller initially budgeted 2–4 h per configuration, or 52–104 h for
26 configurations. This historical allowance is superseded by the direct
16-thread volume-sampled measurements above and must not be used as the
current ETA. The full entry/tolerance matrix and cluster capacity remain
unmeasured; longer times and larger ensembles increase work, subject to loss.
The 256-marker, 2 ms matrix is an initial comparison, not long-time validation.

At the stopped-session checkpoint, W4 had no measured speedup. Its old
pending patch calls `gc_map` after every
accepted FO substep, and that routine calls `spectre_fo_field_locate` before
checking exit. Thus it removes Cartesian inversion from RK stages but retains
an iterative coordinate inversion in every accepted-step exit check. A fast
implementation must also avoid this repeated reconstruction: retain logical
FO state, use a conservative inexpensive exit screen, and reconstruct the GC
only when needed for an actual transition or requested output. Small local
metric solves are distinct from the expensive nonlinear coordinate inversion.

## Focused performance implementation, 2026-10-14

A separate branch, `perf/spectre-fo-logical-spline`, now implements the focused
performance change against main `fe002a3`: [PR #534](https://github.com/itpplasma/SIMPLE/pull/534), commit `2e243d4`. It reuses W4's
per-volume provider, without importing the recovery branch's other concepts or
its draft canonical-momentum rewrite. RK45 advances logical position and
Cartesian velocity; regular interface stages need no position inversion or raw
spectral geometry/field evaluation. GC reconstruction retains a spline-only
inverse. A Cartesian position chart remains below rho=0.05 in volume 1.

The local short NCSX witness (32 markers per surface, RO=100 internally, ten
reference 1 T gyroperiods, tolerance 1e-9) measured 5.1–5.5x faster tracing against
main with the same workload. One/four-thread endpoints agree exactly as printed.
These timings exclude setup and are not a measured 2 ms loss-study ETA. Default
three-volume spline storage is about 1.27 GiB; measured setup was 1.92 s with one
thread and 1.46 s with four. Shared field reinitialization releases old tables.
Maximum sampled relative B-vector error was 2.42e-6 on the default 32/64/32 grid
and 6.60e-5 on the coarse grid. Geometry/A errors were 6.17e-9/1.63e-7.

Full native `make test CONFIG=Fast VERBOSE=1` passed 99/99. After final
terminal-boundary/cache-lifecycle adjustments, the affected SPECTRE and Python
API selection passed 19/19. Existing energy/momentum/event thresholds and stored
golden records were not changed. PR #534 is a draft; CI was not awaited. The old 2–4 h budget above concerns the Cartesian
implementation, not this candidate. Long-time stellarator loss agreement,
spline-grid convergence of losses and the >=10x target remain open. The
archived W4 patch is historical; do not apply it over this implementation.

## Ordered continuation plan

Execution is authorized through the stages below, with regular commits/pushes
and final research plots published through slopbox. Do not wait for CI. Keep
PR #534 focused; study sources and recovered research remain on this branch.

Initial milestone: the native physical study driver builds and passes its
8-marker, 0.2 microsecond NCSX smoke on main and the spline candidate. A
detached baseline at `fe002a3` uses the identical measurement source and an
accepted-FO-time diagnostic; no baseline solver equations are changed. Local
16-marker, 2 microsecond, 3.5 keV paired sheet/volume pilots come first, followed
by grid convergence, controlled throughput/profile measurements and a modest
loss comparison. Keep physical loss, numerical stop and recovered holds
separate. No production switch or symplectic proof is promoted by these runs.

The first 32-case pilot exposed a driver defect: it omitted production's
per-marker `set_spectre_volume_lock(0)` cleanup. The resulting hybrid starts
depended on the preceding marker on a thread; marker 12's NCSX toroidal start
shifted by 8.3e-4 rad and acquired a false FO excursion. A single-marker replay
isolated this. Those results/plots are superseded, kept only in ignored build
scratch, and must not be used as acceptance evidence. The driver now resets
the lock before initialization and after each marker. The corrected 32-case
matrix completed with no numerical stops, FO failures or warning holds.
All paired loss statuses and sampled GC/FO representations agree with main.
Initial GC angle round trips are at roundoff. Default-grid maximum paired FO
position errors are 0.90 micrometers on NCSX and 0.0014 micrometers on G3V3L3Fi;
the NCSX hybrid volume error is 0.20 micrometers. NCSX sheet FO errors decrease
18.55 -> 0.219 -> 0.0319 micrometers on the 16/32/48 radial grids. NCSX volume
losses are 3/16 for FO and 2/16 for the current hybrid; other pilot ensembles
have no loss by 2 microseconds. These are short controls, not long-time model
validation or a 256-marker acceptance study.

Retained inputs, sparse traces, marker outcomes, source provenance, summaries
and figures are in `concept/spectre_loss_study/results/pilot/`; reproduce plots
with `python concept/spectre_loss_study/analyze.py
concept/spectre_loss_study/results/pilot build/study/replot`. The default-grid
short tracing improvement spans 3.7–6.1x across fields and modes; setup is
excluded. The >=10x target and actual 2 ms throughput remain unestablished.

A 16-marker NCSX sheet FO profile at 200 microseconds completed without
stops/holds/failures. Profiling attributed 67.7% of all sampled user CPU time
(including setup) to the general batch first-derivative spline kernel; the
remaining coordinate inversion is a smaller contributor. Filtering call
chains to FO gives 8012/9264 samples (86.5%) in that kernel, with 2072 samples
(22.4%, inclusive) passing through the spline inverse and 4988 (53.8%) through
RK trials. Do not optimize an approximate exit screen before addressing this
measured cost. Its profiled 7.02 s trace is not an uninstrumented throughput
result. The endpoint-sample reuse change passed all 20 affected native tests in
226.64 s, and all four NCSX/G3V3L3Fi physical volume replays (FO/hybrid) retain
bitwise-identical sparse trajectory arrays. All 130 printed throughput endpoint
records also match the pre-reuse implementation. The RO=100 witness now takes
0.433/0.454 s on one thread and 0.117/0.123 s on four threads: about 1.35x faster
than the previous spline implementation, or 6.8–7.4x faster than main in this
short nonphysical witness. No field-order or transition-threshold change was
introduced. This is pushed as `6bbdef0` on PR #534 and `f876d28` on this branch.

The initial longer stage (before the axis repair) ran locally: 32 matched markers, 200 microseconds, 16
physical cores, both equilibria and sampling sets, FO/current hybrid, main vs
default/refined splines (24 cases). No cluster jobs were launched. FO tracing
improves 5.5–7.3x; all paired loss statuses match. However, the NCSX volume
hybrid is invalid: marker 31 has one FO locate failure and 4,666 held intervals
on both spline grids, while main has none. Do not promote the matching loss
counts as acceptance evidence or run the large hybrid matrix over this defect.

The cause is established independently: evaluating polar spline derivatives
over the vanishing axis Jacobian gives relative B error 49.0 at rho=1e-12.
See ledger F7 and `concept/spectre_loss_study/results/axis-failure.json`.
A regular xi/eta axis table now gives 4.21e-9 relative B error and resolves the
pole-seeded inverse. The exact marker replay has one entry/exit, zero holds or
failures, and advances to rho=0.2744 at 200 microseconds. All 21 affected native
tests pass in 229.96 s, including the G3V3L3Fi axis limit and transition-region
checks. This is pushed as `47d3fb0` on PR #534 / `9c62300` on this branch;
the repaired longer controls are complete as recorded above. This fix
retains cheap field evaluation without a raw fallback or a new setting.
The study controller now includes the provider in its source hashes and refuses
unresolved intervals as scientific acceptance. The later GC audit is complete to the stated source/numerical scope and
collected here; its three verified repairs are isolated in PR #535. Builds
were kept separate from the primary throughput measurements.
`full.py` completed the 256-marker, 2 ms volume controls on 16 local cores,
followed by both GC-fixed hybrid replays; see the current checkpoint above.
This answers the actual cost question before the full entry/tolerance matrix.
It is conditional on the current deterministic entry phase. FO uses particle
position at the outer edge; GC segments use the guiding-centre edge convention.
Their loss difference therefore does not isolate interface physics or establish
a coherent first-order switch. The later W2 acceptance matrix remains open.
Generated PNGs follow the repository's ignore rule; their versioned NumPy data
and plotting source reproduce them. The convergence plot was uploaded to
`https://box.sloppy.at/1dcca.png` on 2026-10-14 (slopbox links expire in 3 days).

1. **W4 first: review and validate the focused performance PR above.** The
   regular logical-position/Cartesian-velocity path is implemented and locally
   checked, including longer NCSX/G3V3L3Fi trajectory/loss and two-grid controls.
   Review the retained spatial/loss differences against the required accuracy
   before promotion. Keep the existing physical transition criteria.
   Profile the remaining spline-only GC reconstruction before introducing a
   conservative exit screen. Measure end-to-end cost, including setup and FO
   residence, under controlled load. The repaired 200 us tracing improvement is
   5.8–7.4x; the >=10x per-step target is not established. Current main retains
   its Cartesian path until a reviewed performance change is promoted.
2. **Finish the production GC proof bridge (three source repairs verified; uniform proof still open).** Map each of the eight variants'
   actual equations and field/extrapolation operations to its theorem. The
   stellarator one-step/100-step numerical checks pass; preserve their bounds
   (1e-7/1e-5), spacing and negative controls. Quantify finite-Newton and field
   errors and establish a suitable smooth-domain/time-step hypothesis for the
   unresolved long-energy claim. Preserve the failed energy/domain evidence.
3. **W1/W3: validate and integrate the consistent switch.** Item D's native build, local identity checks,
   regeneration and independent numerical pullback are complete. Implement
   the extended canonical action/phase and gyrogauge connection tables, then
   audit physical units and matching to production GC variables. Reuse the
   recovered local Darboux completion; a gyrofrequency counter on fixed-mu
   GC is insufficient. Phase-complete first-order evolution also needs the
   consistent H2 frequency correction; investigate splining its coefficient
   fields to keep tracing cheap, rather than evaluating high derivatives per
   stage. These tables and the real-field six-dimensional model are not yet
   implemented. Require scaled round-trip identity
   at 1e-12, O(epsilon^2) mu defect and independent FO reconstruction. Replace
   both entry and exit together with carried phase. Redesign the long-run
   comparison so it distinguishes sheet physics, switch error and artificial
   absorbing boundaries. No production switch change has yet been made.
4. **W2: paired volume-sampled benchmark.** At least 256 matched markers on
   NCSX and G3V3L3Fi, 2 ms initially. Compare entry distances 1.5/2/3/4 rho_L,
   tolerances 1e-7/1e-9/1e-11 and FO everywhere at 1e-9. Publish loss-time curves,
   paired uncertainty, FO residence, excursion-length and cost distributions.
   Control shared-host load before making timing decisions. No such acceptance
   matrix was completed by the stopped session.
5. **W5: encounter map concept, conditional on W2.** Exact velocity-continuity
   reset, carried/flux-weighted phase and short planar return cycles. Stochastic
   phase sampling is a modelling choice requiring ensemble validation, not a
   proof of deterministic symplecticity. Require NCSX/G3V3L3Fi loss agreement
   and a source-tied fixed-model energy proof before considering production.
6. **W6: equilibrium-volume convergence.** Obtain valid NCSX equilibria with
   2/4/8 volumes and repeat matched physical initial states. This determines
   sensitivity to the two-volume jump; more integrator refinement cannot do it.
7. **Eventual symplectic FO route.** A variational event integrator with crossing
   time as an unknown was proposed, not implemented or proved. Reassess after
   switch and cost evidence. Fixed-family Gauss2, unresolved splitting and
   state-dependent landing do not close the current whole-map proof gap.
8. **Parked concepts.** Keep higher-order GC, frozen-step and inexact-Newton
   investigations available. Revisit GC2 production only after the derivative
   provider, gyrogauge tables and actual stellarator benefit are established.

## Recovery verification and publication (historical checkpoint)

At the recovered `082bea4` checkpoint, production `src/` and `examples/`
were identical to main `fe002a3`; this paragraph is historical, superseded by
the active performance and GC changes above. Gauss2 is
linked only into a test-only concept library. W4's replacement path and GC
oracle registrations are preserved in patches, with explicit draft READMEs.

Checks on this consolidated checkout are being recorded below before push.
Historical tests are not promoted into fresh verification. The old final
Gauss2 worktree full suite had a chartmap-pipeline timeout; the GC prototype
worktree had Python-environment failures and timeouts (only some reruns passed).
A local Gremlin status query returned a state-provider error; this recovery uses
the repository's native make/CMake test authorities and starts no resident
workers. No protected hosts were used. Full long-orbit research sweeps remain
future work, not a condition for preserving this stopped state.

- Native deterministic Fast build and `make test-fast CONFIG=Fast`: **87/87
  passed** (286 s), including all fast FO stellarator map probes, the Gauss2
  concept field/gauge tests, production event tests and W1 diagnostic. The
  old chartmap-pipeline timeout passed here with its unchanged timeout.
- Recovered Lean checks: **passed**, standard axioms only. GC checker requires
  22 results; FO Gauss2, EntryFrozen and LieTruncation also checked. This does
  not certify the separated unchecked draft or the actual hybrid map.
- Provider-jet baseline: **1149/1149 exact**, no probes. Original rejection
  verifier: **passed**. The unfinished mutation harness is still a recorded
  limitation. Registered FO generated-source/readback check via `make test`: **1/1 passed**
  (160 s), including regenerated leaves and exact readback identities.
- Standalone candidate behavioural/symbolic checks: **passed**. GC concept
  `test`, `canonical-test`, `canonical-check-results`: **passed**. FLR note:
  **13 exact checks passed**; local-map checks passed (200-cycle discrepancy
  8.0e-13, FD symplectic defect 1.9e-10). Existing generated-GC continuation
  count warnings remain; no thresholds were relaxed.
- FortSym clean build and registered Lie/products/CSE-DAG tests: **3/3 passed**.
  Consumer generators and concept behavioural tests above also passed with it.
- Source inventory: all surviving nonignored untracked files from derive,
  gcproof and W4 are collected, with draft integrations preserved explicitly.
  Production `src/` and `examples/` have no delta. Diff inspection retained
  whitespace in literal patch context and unchanged generated EOF blanks so
  archived patches and exact regeneration output stay intact. No active-source
  whitespace finding remains. The full slow suite and long research sweeps
  were not repeated; fresh coverage is the fast suite plus the focused native
  generator test and standalone checks above.

The cleanup is on main via #533 and this recovery branch is rebased onto it.
The deleted implementations were not restored by replaying recovery commits.
The original worker histories remain in the published pre-cleanup archive;
the current branch carries the surviving concepts, proofs and unfinished drafts. The clean
FortSym `feat/lie-pullback` dependency branch is pushed separately in its owning
repository. PR #533 publishes the bounded failed-path cleanup on main; the remaining
research concepts/drafts stay on this recovery branch.
For continuation, start with item 1 above and update this ledger with results
for the actual next revision; do not infer closure from a worker's last message.

### Cleanup and rebase verification

- Cleanup: full existing native `make test CONFIG=Fast`: **95/95 passed**,
  including 15 slow tests (525 s). Existing standalone behavioural and symbolic
  checks passed; the remaining original sheet studies built without warnings.
- Five paired before/after ensembles, 256 markers each: FO A/B, carried-phase
  C3 A/B (switches exercised), resolved GL2 A. Every recorded numerical value
  matched exactly; wall-time metadata was excluded.
- No tests or runtime checks were added for rejecting retired parameters.
- Merged on explicit instruction without waiting for CI. GitHub required the
  Build and Test status, so administrator merge was used; that check was still
  running. Local verification above passed. Add/remove commits remain reachable
  from main through a merge commit.
- Rebased recovery focused native checks: the existing switch round-trip,
  FO symplectic map and Gauss2 concept-field tests passed (3/3). Main's cleaned
  concept trees, ledger, production sources and settings are retained verbatim.
