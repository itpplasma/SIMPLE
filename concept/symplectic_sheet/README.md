# Symplectic full orbits through a current sheet

Standalone Fortran prototype, not part of the SIMPLE build: `make`, then run
the programs in `build/`. It studies fixed-step symplectic integration of
charged-particle full orbits through an MRxMHD-type current sheet, against an
adaptive RK45 reference.

## Model

`sheet_field.f90`: axisymmetric two-volume torus with continuous vector
potential `A_Z = -F ln R`, `A_phi = psi(r)/R`, `psi` continuous with a jump
of `psi'` at `r0`. The field direction jumps across the sheet (angle set by
`a_out`), `p_phi` is conserved, and the full orbit reduces to two degrees of
freedom with a Hamiltonian that is continuous with a jump in its gradient.
`width > 0` replaces the sheet by a C2 layer of that width.

The force jump at the sheet is `v_phi x [[B_p]]`. With
`p_phi = psi(r0)` the crossing markers have `v_phi = 0` there and barely feel
the sheet; the studies use `p_phi = 0.45` (`v_phi ~ 0.4 v`).

## Integrators

- `rk_reference.f90`: Dormand-Prince 5(4) with exact branch switching at the
  sheet (false-position landing, grazing touches resolved).
- `fo_sheet.f90`, fixed step `h`:
  - `SCHEME_MIDPOINT_REGION`: implicit midpoint on the piecewise `H`, branch
    from the midpoint; on a smooth layer this is the plain midpoint.
  - `SCHEME_NAIVE`: branch of the step start for the whole step.
  - `SCHEME_LANDING`: exact landing on the sheet and energy-matched normal
    kick (symplectic collision construction); exactly symplectic for
    transversal crossings (`tangent_map`, analytic), ill-defined at grazing
    incidence.
  - `poincare_step`: midpoint for `K = g(z)(H - E)`, a Poincare time
    transformation that shortens physical steps near the layer.

## Findings (strong sheet `a_out = 2`, `p_phi = 0.45`, ~8000 crossings per marker)

| scheme | h | width | rms dmu/mu | steps per unit time |
|---|---|---|---|---|
| RK45 reference, sharp | adaptive | 0 | 1.69e-2 | - |
| RK45, smooth layer | adaptive | 1e-2, 2.5e-3 | 1.68e-2, 1.71e-2 | - |
| midpoint, sharp | 0.1 / 0.05 / 0.025 / 0.0125 | 0 | 0.186 / 0.093 / 0.056 / 0.028 | 10-80 |
| naive, sharp | 0.1, 0.025 | 0 | 0.85 (mu destroyed) | 10, 40 |
| landing, sharp | 0.1 / 0.05 / 0.025 | 0 | 1.9 / 0.24 / 1.9e-2 | 10-40 |
| midpoint, smooth, resolved | 0.025 | 1e-2 | 1.68e-2 | 40 |
| Poincare, smooth | 0.1 | 2.5e-3 | 1.9e-2 to 2.4e-2 | 12-25 |

1. A finite layer of width `rho/20` to `rho/5` reproduces the sharp-sheet
   physics: the full-orbit zero-width limit is benign.
2. A fixed step that does not resolve the force jump makes a local error
   `~ h [[F]]` per crossing. Taking the branch of the midpoint makes its mean
   zero (sawtooth in the crossing fraction), taking the branch of the step
   start biases it. The zero-mean error still random-walks: spurious mu
   scattering grows like `h sqrt(N)`. Exact or piecewise symplecticity does
   not prevent it; resolving the jump does.
3. The energy-matched landing scheme fails at grazing incidence, where the
   matching kick exceeds the normal velocity.
4. The Poincare refinement is exactly symplectic and cheap but leaves an
   excess scattering that grows with the refinement length; the suspected
   cause is the term `(H - E) grad g` acting on the midpoint energy error.

`test_fo_symplectic` checks the analytic tangent map of landing steps
(`M^T J M = J` to 1e-13); `test_fo_stationary` checks long-run stationarity
of the microcanonical ensemble; `study_mu_scattering` produces the table.
