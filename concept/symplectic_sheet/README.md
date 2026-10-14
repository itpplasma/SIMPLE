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
| RK45, tol 1e-10 | - | 50.8, 3.4e-5, 1.66e-2 |
| Poincare Gauss-Legendre 2, h = 0.1, L = 5e-3 | 13.4, 1.1e-4 | 128, 4.2e-4, 1.50e-2 |
| same, analytic layer, L = 1e-2 | 19.3, 5.4e-7 | 193, 5.1e-7, 1.61e-2 |

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
5. The energy error of the symplectic scheme stays bounded (5e-7 at both
   times) once every time scale of the transformed Hamiltonian is resolved:
   with `L = 5e-3` the transition of `g` was crossed in two fictitious steps
   and the energy random-walked (1e-4 to 4e-4); `L = 1e-2` fixes it at 15 %
   more steps. RK45's energy error grows linearly in time at any tolerance.
6. In this two-dimensional prototype with a cheap field, RK45 at tol 1e-10
   is about 4 times faster than the bounded-energy symplectic scheme over
   4e5 time units, with an energy error 70 times larger and growing.

## Boris

`boris_step` (Cartesian, field of the local region) keeps `|v|` and phase-space
volume exactly but is not symplectic. Strong sheet (`a_out = 2`,
`p_phi = 0.45`), rms dmu/mu against the RK45 reference 1.68e-2:

| field | h | rms dmu/mu | wall s (T = 4e4) |
|---|---|---|---|
| sharp | 0.1 / 0.05 / 0.025 | 0.166 / 0.090 / 0.049 | 1.1 / 2.2 / 6.4 |
| layer w = 1e-2 | 0.1 / 0.05 / 0.025 / 0.0125 | 2.44e-2 / 1.73e-2 / 1.66e-2 / 1.67e-2 | 1.0 / 2.4 / 4.8 / 8.4 |

On the sharp sheet it has the same `h sqrt(N)` scattering as every unresolved
fixed step. On a resolved layer it is accurate and as fast as RK45 with an
exact energy; at h = 0.05 over T = 4e5 the spread reaches 1.84e-2 against
1.66e-2, so h = 0.025 is the safe choice there.

## GC away from the sheet (`gc_hybrid.f90`, `study_hybrid.f90`)

GC drift equations (RK45, tol 1e-10) where the GC is farther than `D_OUT`
gyroradii from the sheet, FO (RK45 on the sharp sheet) within `D_IN`.
FO -> GC is the first-order transform with mu chosen to conserve the energy
exactly; GC -> FO lifts with the advanced gyrophase. Energy is conserved
through every switch (round trip error O(rho^2) in position, 1e-13 in
energy). Weak poloidal field (`a_in = 0.25`, `a_out = 0.75`, q ~ 2), D_IN = 3,
D_OUT = 4, T = 4e4:

| v (rho) | switches per marker | rms dp_phi (FO everywhere: 0) | rms dmu/mu hybrid / FO everywhere |
|---|---|---|---|
| 0.01 | 1.6 | 3.2e-5 | 1.8e-3 / 0.5e-3 |
| 0.02 | 2.4 | 2.1e-4 | 1.6e-3 / 0.9e-3 |
| 0.04 | 3.0 | 9.8e-4 | 3.4e-3 / 1.05e-3 |

Each switch makes an error set by the order of the GC transform: the
canonical momentum p_phi, exactly conserved physically, jumps by roughly
`rho^2` per switch, and these errors random-walk with the number of
switches. Symplecticity does not remove this; a higher-order transform or
fewer switches (hysteresis) does. A lift that picks the gyrophase keeping
p_phi exact (axisymmetric only) reduces but does not remove it.

## Real SPECTRE interfaces

`tools/spectre_interface_jumps.f90` (built against a SIMPLE build, see its
header) evaluates the true field at `rho_g = k -+ 1e-10` on a 128 x 64 grid:

| fixture, interface | direction jump max (rms) | `|B+|/|B-| - 1` |
|---|---|---|
| tok2vol_sheet | 0 | +1.1 % |
| tok2vol | 0.045 rad (0.029) | +1.0 to +1.2 % |
| NCSX free boundary, 1 | 0.054 rad (0.020) | +1 to +8 % |
| NCSX free boundary, 2 | 0.32 rad (0.11) | -25 to +38 % |

The outer NCSX interface is a strong sheet, comparable to or stronger than
the prototype's `a_out = 2` case. The tangential vector potential jumps at
every interface (per-volume gauge), so a full orbit in a global chart needs
the gauges of neighbouring volumes matched first.

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
Gauss-Legendre method, which removes the residual (table above). The
transformed Hamiltonian has the extra time scale `L/(g v)` of the transition
zone; a fixed step in `tau` keeps the energy error bounded only when it
resolves that scale too, as for any symplectic method.

`test_fo_symplectic` checks the analytic tangent map of landing steps
(`M^T J M = J` to 1e-13); `test_fo_stationary` checks long-run stationarity
of the microcanonical ensemble; `study_mu_scattering` produces the table.

## SIMPLE on NCSX: loss fractions against full orbit everywhere

`ncsx_freeboundary.h5` (Mvol = 3), 3.5 keV alphas (`facE_al = 1000`), 64
markers started next to an interface, `trace_time = 2e-3`, 32 threads. The
reference is full orbit everywhere: crossing_level 3, RK45 at tolerance 1e-9,
`spectre_fo_entry_larmor = 1e6`.

| model | s = 0.95 (interface 1): lost, wall s | s = 1.95 (interface 2): lost, wall s |
|---|---|---|
| full orbit everywhere, RK45 (reference) | 6, 920 | 39, 565 |
| level 3, RK45 | 7, 929 | 41, 409 |
| level 3, Boris (layer 0.2, rotation 0.0125) | 5, 1399 | 48, 480 |
| level 1, symplectic GC | 7, 15 | 29, 25 |
| level 1, RK45 GC | 1, 17 | 7, 30 |

Counting noise is about +-2 markers at interface 1 and +-4 at interface 2. At
the strong outer interface both guiding-centre maps confine particles
spuriously (RK45 GC: 32 of the 39 reference losses are missed). Level 3 with
RK45 matches the reference at both interfaces; Boris overestimates the losses
at interface 2 by about 2 sigma. Individual marker fates differ even between
full-orbit runs (chaotic orbits near the interface); only the fractions
compare. The hybrid is not cheaper than full orbit everywhere here: the
markers stay near the interfaces.
