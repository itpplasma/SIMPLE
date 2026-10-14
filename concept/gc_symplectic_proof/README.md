# GC symplectic proof work recovered after the stopped session

`lean/GcSymplectic.lean` contains the generic PRK theorem and exact-tableau
instantiations for both Euler variants, midpoint/Gauss1, Gauss2–4 and Lobatto3.
`check_lean.sh` requires 22 theorem axiom lists on the standard Mathlib axioms.
Run it from the SIMPLE root. `check_tableaux.py` compares the production
coefficient literals with the ideal tableaux and writes `tableaux.json`.

This proves the ideal exact-stage methods under the stated hypotheses. The
bridge to the actual field evaluation, extrapolation, terminated Newton solve,
radial inversion and state-dependent SPECTRE landing is unfinished.

`test_gc_symplectic.f90` is the worker's numerical oracle draft. Its proposed
CMake registration is preserved in `pending-integration.patch`, outside the
default test suite. To resume, copy the Fortran file to `test/tests/` and
integrate those registrations on a development branch. Only its analytic
tokamak smoke test completed before the stop; stellarator acceptance is open.
The predeclared defect thresholds are 1e-7 for one step and 1e-5 for 100 steps,
with spacing disagreement below the same threshold. Do not change these bounds
to conceal a failing result.

The worker reported possible midpoint/Gauss extrapolation defects in its source
audit. Those require independent reproduction before a production correction;
no correction is included in this recovery.
