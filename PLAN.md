# SPECTRE orbit programme: recovered status and plan

Updated 2026-10-14. Recovery branch: `recovery/spectre-session-20261014`.
Production baseline: `origin/main` at `fe002a3`, including cleanup PR #533.
This consolidates the stopped Claude controller and its GPT Sol, strategy,
Opus and nested derivation workers. It records unfinished research rather
than claiming the scientific programme is complete. No previous worker is
running or being resumed by this recovery.

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
  the ancestry of #530; the live main ref remains `45c0598`.
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
| Opus W4, SIMPLE-w4 uncommitted | `concept/fo_logical_spline/` | Generator, emitted RHS, spline provider and benchmark source retained. All tracked production edits in `pending-integration.patch`; not enabled. Build, interface fidelity and speedup validation unfinished. |
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

GC worker source-audit leads: possible spurious midpoint extrapolation term,
Gauss derivative-index and omitted vpar extrapolation term. These were not
independently established before the stop. Reproduce each against production
before fixing it. A coefficient theorem is not a proof of those code paths.

## Ordered continuation plan

1. **Finish the production GC proof bridge.** Map each of the eight variants'
   actual equations and field/extrapolation operations to its theorem. Restore
   the draft oracle on a development branch; run NCSX canonical-flux/Boozer and
   QH cases. Preserve one-step/100-step defect gates (1e-7/1e-5), spacing controls
   and the negative control. Repair only reproduced implementation defects.
2. **W1/W3: validate and integrate the consistent switch.** Complete item D's
   generated identity checks and build recipe; audit units, gyrophase, gauge and
   the relation to production GC variables. Require scaled round-trip identity
   at 1e-12, O(epsilon^2) mu defect and independent FO reconstruction. Replace
   both entry and exit together with carried phase. Redesign the long-run
   comparison so it distinguishes sheet physics, switch error and artificial
   absorbing boundaries. No production switch change has yet been made.
3. **W4: finish splined logical-coordinate FO.** Resume the preserved sources
   and patch; validate derivatives and metric/gauge continuity per volume.
   Resolve rho-root events and crossing-and-return pairs without Cartesian
   inversion. Test axis/exterior ownership and compare energy/trajectory/loss
   results to current sharp-field Cartesian RK45. Target at least 10x per-step
   speedup. Retire the old path only after equivalent physical results and
   throughput are measured; the speedup is currently a target, not a result.
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

## Recovery verification and publication

Production `src/` and `examples/` are identical to main `fe002a3`. Gauss2 is
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
