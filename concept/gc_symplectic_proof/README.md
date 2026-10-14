# GC symplectic proof work recovered after the stopped session

`lean/GcSymplectic.lean` contains the generic PRK theorem and exact-tableau
instantiations for both Euler variants, midpoint/Gauss1, Gauss2–4 and Lobatto3.
`check_lean.sh` requires 22 theorem axiom lists on the standard Mathlib axioms.
Run it from the SIMPLE root. `check_tableaux.py` compares the production
coefficient literals with the ideal tableaux and writes `tableaux.json`.

This proves the ideal exact-stage methods under the stated hypotheses. The
bridge to the actual field evaluation, extrapolation, terminated Newton solve,
radial inversion and state-dependent SPECTRE landing is unfinished.

`test_gc_symplectic.f90` is the worker's numerical oracle draft, now registered
from its concept source in the native CMake suite on this audit branch. The
original registration proposal remains in `pending-integration.patch` for
provenance. Only its analytic tokamak smoke test completed before the stopped
session; fresh stellarator acceptance is pending below.
The predeclared defect thresholds are 1e-7 for one step and 1e-5 for 100 steps,
with spacing disagreement below the same threshold. Do not change these bounds
to conceal a failing result.

The worker reported possible midpoint/Gauss extrapolation defects in its source
audit. Those require independent reproduction before a production correction;
no correction is included in this recovery.

## Consolidated checkout

The native audit registrations are retained in `pending-integration.patch`,
rebased to this consolidated checkout. They are active and reproducible on
`research/gc-source-audit-20261014`; apply the patch in a research worktree to
resume them here. The current default suite keeps the independent production
regression but does not enable the known RED long-energy acceptance. The
original energy bounds, starts and failure records remain intact. This is a
scientific acceptance gap, not a reason to remove valid GC integrators.

## Source audit resumed, 2026-10-14

The closed-form polynomial canonical Hamiltonian witness in
`test/tests/test_gc_field_extrapolation.f90` fails against unchanged production
source (`9c62300`). It compares the actual Gauss update with independent
Hamiltonian quadrature at the solved stages and checks the midpoint endpoint
momentum against its closed expression. It exercises Gauss1–4 with/without a
cyclic toroidal coordinate, and both midpoint volume branches. The error bound
is quadratic in the last Newton correction plus binary64 roundoff. The failing
records are retained in `extrapolation_results.json`; bounds are frozen before
repair. The existing source-tableau exact-arithmetic check passes.

Candidate corrections are the endpoint-only midpoint Taylor increment, the
Gauss r/phi mixed derivative index and Gauss parallel-speed p_phi increment.
The three corrections now pass the frozen witness. Restoring each defect
individually fails it again; those isolated records are retained alongside RED
and GREEN. The standalone witness is numerical implementation evidence, not a
proof of the floating-point whole orbit map.

The first stellarator map run passes the NCSX canonical-flux and QH Boozer
cases and all one-step checks, but the original fourth-order FD stencil fails
two 100-step NCSX Boozer checks: midpoint/Gauss1 give 2.028e-5 at spacing 1e-4
and 2.280e-6 at 5e-5. Unchanged pre-fix production gives the same values.
The original failure is retained; oracle-resolution work is separate from the
three source repairs. No defect or spacing-disagreement bound is relaxed.

An eight-point central derivative now resolves all registered map cases at the
same spacings (1e-4, 5e-5) and unchanged bounds (1e-7 one step / 1e-5 for 100
steps, also for spacing disagreement). All five native selections pass in
44.54 s, including the independent Hamiltonian witness and analytic smoke.
The previously failing NCSX midpoint/Gauss1 case gives 6.895e-6 / 6.264e-7,
with spacing disagreement 2.583e-6. This improves the numerical oracle's
resolution; it is not an exact floating-point symplecticity proof. Both
stencils' records and the pre-fix failure are in `map_results.json`.
The 100,000-step stellarator energy checks and NCSX step-halving control are complete; their original acceptance failures remain open as detailed below.

## Long energy/domain result

QH Boozer passes all eight methods at both starts for 100,000 steps. NCSX
canonical-flux passes all methods at the passing start, but the trapped start
reaches the chart edge near step 9,800. NCSX Boozer's trapped start also reaches
the chart edge, with boundary status 4 or unresolved landing status 5 near
s=1. Neither is a completed 100,000-step smooth-volume orbit; status 5 is not
relabelled a physical loss.

At the NCSX Boozer passing start the coarse Gauss2 run fails the original
second-half envelope gate (max2/max1=2.079, energy extrema -1.576e-4/7.582e-5).
Halving h and doubling the steps to keep the same elapsed interval passes all
eight passing-start checks; Gauss2 extrema become -2.645e-7/2.275e-7 and
max2/max1=0.987. The trapped start still reaches the domain edge. The original
energy acceptance selections therefore remain RED: no threshold is relaxed,
start removed, domain exit turned into completion or proof closure claimed.
Sparse samples, all-step summary metrics and the plotting source are retained
under `results/energy/` and `plot_energy.py`.

The first unsupported proof bridge remains the identification of the ideal
exact-stage, fixed-step canonical theorem with the terminated-Newton,
extrapolated-field, domain-limited production map. The three verified source
repairs remove linear extrapolation defects; finite numerical oracles do not
supply missing uniform derivative/domain bounds. State-selected landing is
also outside the fixed-step theorem. The ideal tableau theorems survive with
their hypotheses; a whole-production/hybrid theorem remains open.
