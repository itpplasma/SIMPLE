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
| Poincare (midpoint), smooth | 0.1 | 2.5e-3 | 1.9e-2 to 2.4e-2 | 12-25 |
| Poincare (Gauss-Legendre 2), smooth | 0.1 | 2.5e-3 | 1.70e-2 | 12.5 |

Cost and energy (256 markers, 32 threads, wall time; rms relative energy
error at the end):

| scheme | T = 4e4: wall s, dH/H | T = 4e5: wall s, dH/H, rms dmu/mu |
|---|---|---|
| RK45, tol 1e-6 | 1.8, 2.0e-4 | 16.9, 2.0e-3, 1.72e-2 |
| RK45, tol 1e-8 | 2.6, 1.0e-4 | 24.6, 1.0e-3, 1.68e-2 |
| Poincare Gauss-Legendre 2, h = 0.1 | 13.4, 1.1e-4 | 128, 4.2e-4, 1.50e-2 |

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
4. The Poincare refinement with the midpoint leaves an excess scattering
   from the residual `K = k` acting as the potential `k/g` (below); with the
   fourth-order Gauss-Legendre base method it matches the reference.
5. In this two-dimensional prototype with a cheap field, RK45 is 5 to 7 times
   faster than the Gauss-Legendre Poincare scheme. Its energy error grows
   linearly in time; the symplectic scheme's grows more slowly (factor 3.8
   for ten times the time) but is not bounded, plausibly because the C2
   layer limits backward error analysis. The 10 % lower long-time mu spread
   of the symplectic run is unexplained at this ensemble size.

## Error of a fixed step across a force jump

Let a step of length `h` cross the sheet at the fraction `theta` of the step,
and let `[[f]] = f+ - f-` be the jump of the Hamiltonian vector field (for the
full orbit a pure force jump, `v x [[B]]`). A one-step method that evaluates a
single branch for the whole step misses the impulse of the other branch on
its share of the step. To leading order:

- branch of the step start: `e = -(1 - theta) h [[f]]`; mean over `theta`
  `-h [[f]]/2`. Out- and inbound errors have opposite signs only if the
  branches swap roles, so in general the bias does not cancel.
- branch of the midpoint: `e = theta h [[f]]` for `theta < 1/2` and
  `-(1 - theta) h [[f]]` for `theta > 1/2`; zero mean, rms `h |[[f]]|/sqrt(12)`.

`theta` is equidistributed over many crossings, so the midpoint-branch errors
act as zero-mean kicks of size `~ h [[f]]` and random-walk: after `N`
crossings the spurious change of an invariant is `~ sqrt(N) h [[f]]`. The
map is piecewise symplectic (exactly symplectic on a smooth layer), but that
does not bound the walk: near the jump the map is not a small smooth
perturbation of the exact one, so no modified invariant survives. Only
resolving the jump removes the kicks: exact crossing times, or a smooth layer
with physical steps `v dt << w` inside it.

A Poincare time transformation `K = g(z) (H - E)` with a fixed step in `tau`
resolves the layer while keeping exact symplecticity. On the exact flow
`K = 0`. The discrete flow keeps `K = k` bounded but nonzero, and on `K = k`
the motion is the physical one on the shell `H = E + k/g(q)`: `k/g` acts as
an extra potential across the transition zone of `g`, with an impulse
`~ k (1 - g_min)/v_n` independent of the zone length. `k` is the energy error
of the base method, `O(h^2)` for the midpoint and `O(h^4)` for the two-stage
Gauss-Legendre method, which removes the residual (table below).

`test_fo_symplectic` checks the analytic tangent map of landing steps
(`M^T J M = J` to 1e-13); `test_fo_stationary` checks long-run stationarity
of the microcanonical ensemble; `study_mu_scattering` produces the table.
