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
They are not applied in this checkpoint. The standalone witness is numerical
implementation evidence, not a proof of the floating-point whole orbit map.
