# Symplectic full orbits through a current sheet

Standalone Fortran reference/convergence experiment, outside SIMPLE's CMake
build. Run from the SIMPLE root:

```sh
make -f concept/symplectic_sheet/Makefile
```

The executable study selectors are `rk45`, `poincare` and `poincare_gauss2`;
`study_crossing` compares single crossings (below).
The latter two require a resolved smooth layer and compare that explicitly
different field model. The `naive`, `midpoint_region` and `landing` selectors,
sharp-branch/energy-kick implementation and their dedicated tests were removed.
See the [removal ledger](../../DOC/spectre-rejected-interface-models.md#implementation-removal-ledger).
The following measurements and descriptions of removed schemes are historical.

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

## Integrators in the historical comparison

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

## Historical Boris comparison (removed)

The Boris implementation is retained only in git history; see
[the removal record](../../DOC/spectre-interface-crossing.md#removal-of-the-smoothed-layer-boris-option).
These are recorded measurements. The former `boris_step` (Cartesian, field of the local region) keeps `|v|` and phase-space
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

## Single crossings with SIMPLE's entry lift (`study_crossing.f90`)

SIMPLE's hybrid enters the full orbit with a uniform gyrophase hashed from the
GC state and, in axisymmetry, re-solves `v_par` at fixed speed so that the
particle's `p_phi` equals the GC's `psi(X) + R v_par b_phi` (`LIFT_UNIFORM`,
the rule of `entry_gyrophase` and `match_pzeta` in `src/spectre_fo_hybrid.f90`);
the exit is the first-order transform with mu from the energy.
`LIFT_UNIFORM_AVG` draws the same phase but targets the gyrophase average
(trapezoid, 16 phases) of the lifted particle's `p_phi` at the unshifted
`v_par`. `gyrophase` is the generic lift without re-solve. The study starts GCs
4 gyroradii inside the sheet (`a_in = 0.25`, `a_out = 0.75`, poloidal angle
`0 +- 0.1`, `v_par/v = 0.3 +- 0.05`) and follows each until its GC is
`D_OUT = 4` gyroradii from the sheet on either side (`D_IN = 3`). The
reference is RK45 full orbit everywhere from an independent draw of the same
GC distribution, lifted with a uniform gyrophase. mu and the GC position are
averaged over one gyroperiod after entry and after the exit; the displacement
is the poloidal arc of the GC from entry to exit. The `p_phi` jump is that of
one further FO -> GC -> FO round trip at entry, with a new gyrophase. Run as

```sh
study_crossing.x 0.02 4096 0 0.2 0.3 0.1 \
    [uniform|uniform_avg|gyrophase|avgentry|avgentry_p] \
    [pinned|plain|avgexit|avgexit_p|avgexit_x] [FILE]
```

The table below measures mu at both ends against the starting GC's mu and
uses the plain first-order exit (`plain`); its `p_phi` jump is that of the
plain exit. Hybrid against full orbit, N = 4096 markers per model
(transmitted fraction +-0.008; KS: two-sample Kolmogorov-Smirnov `D` and
p-value):

| lift | v (rho) | transmitted, hybrid / full | `mu1/mu0 - 1` mean, rms: hybrid / full | KS mu D (p) | KS displacement D (p) | `p_phi` jump rms (mean) |
|---|---|---|---|---|---|---|
| uniform (SIMPLE) | 0.01 | 0.615 / 0.619 | 1.37e-4, 1.6e-5 / -0.2e-5, 1.4e-4 | 0.70 (< 1e-300) | 0.026 (0.12) | 3.1e-5 (-2.4e-5) |
| | 0.02 | 0.609 / 0.610 | 2.81e-4, 3.5e-5 / 0.2e-5, 2.8e-4 | 0.69 (< 1e-300) | 0.043 (1.2e-3) | 1.25e-4 (-0.97e-4) |
| | 0.04 | 0.588 / 0.585 | 5.64e-4, 6.7e-5 / 0.3e-5, 5.6e-4 | 0.70 (< 1e-300) | 0.062 (2e-7) | 5.0e-4 (-3.9e-4) |
| uniform_avg | 0.01 | 0.613 / 0.619 | -0.3e-5, 0.3e-5 / -0.2e-5, 1.4e-4 | 0.49 (< 1e-300) | 0.011 (0.98) | 1.9e-5 (-0.01e-5) |
| | 0.02 | 0.604 / 0.610 | -0.1e-5, 0.8e-5 / 0.2e-5, 2.8e-4 | 0.49 (< 1e-300) | 0.016 (0.68) | 7.9e-5 (-0.3e-5) |
| | 0.04 | 0.581 / 0.585 | -0.6e-5, 1.5e-5 / 0.3e-5, 5.6e-4 | 0.49 (< 1e-300) | 0.020 (0.41) | 3.2e-4 (-0.9e-5) |
| gyrophase | 0.01 | 0.613 / 0.619 | -0.2e-5, 1.2e-4 / -0.2e-5, 1.4e-4 | 0.067 (2e-8) | 0.011 (0.96) | 2.8e-5 (-0.08e-5) |
| | 0.02 | 0.604 / 0.610 | -0.5e-5, 2.5e-4 / 0.2e-5, 2.8e-4 | 0.063 (1e-7) | 0.016 (0.63) | 1.1e-4 (-0.5e-5) |
| | 0.04 | 0.579 / 0.585 | -0.4e-5, 5.0e-4 / 0.3e-5, 5.6e-4 | 0.057 (4e-6) | 0.030 (0.05) | 4.5e-4 (-1.1e-5) |

KS mu `D` at N = 256 / 1024 / 4096 (v = 0.02; 5 % critical `D` 0.120 /
0.060 / 0.030): uniform 0.68 / 0.68 / 0.69, uniform_avg 0.49 / 0.50 / 0.49,
gyrophase 0.118 / 0.082 / 0.063. The other speeds behave the same.

The transmitted fraction agrees for every lift. With the SIMPLE lift, mu is
offset by `+1.4e-4 v/0.01` with a narrow spread, and the disagreement does not
shrink with N. The offset is present right after the lift: the particle's
`p_phi` averaged over the gyrophase differs from `psi(X) + R v_par b_phi` by
`O(rho^2)`, and matching every particle to the latter shifts `v_par` by
`O(rho^2/R)` in one direction. The same offset gives a systematic `p_phi`
drift of `-1e-4` per round trip at v = 0.02 (mean close to the rms), and
shortens the reflected markers' displacement (134.42 against 134.55 rho at
v = 0.02, KS p 7e-9; transmitted p 0.57).

Targeting the gyroaverage removes the offset and the `p_phi` drift and makes
the displacement agree at all N, but it narrows the mu distribution to rms
0.8e-5 against 2.8e-4 (v = 0.02), so KS mu fails with `D` about 0.49 at
every N. The round-trip `p_phi` jump keeps a zero-mean random part of
7.9e-5 rms, 30 % less than without any re-solve (1.1e-4): the first-order
FO -> GC transform does not return a GC whose gyroaveraged lift has the
particle's `p_phi`, so the next switch still changes it. Without the re-solve
the mean of mu agrees, but the width (2.5e-4 against 2.8e-4) still gives KS
`D` 0.063, p 1e-7 at N = 4096.

The spread of the reference is itself set by its first-order lift at the
start, while one crossing of this sheet changes mu by only about 5e-6 rms at
v = 0.02 in both models. mu statistics at the 1e-4 level here therefore test
the consistency of the first-order lifts more than the sheet physics.

### Full passage with SIMPLE's exit, measured GC to GC

SIMPLE's exit in axisymmetry is not the plain transform: it takes the
first-order GC position, sets the GC `v_par` so that the GC `p_phi` equals the
entry GC value, keeps the speed and takes mu from the energy (`pinned`,
`fo_to_gc_pinned`). With it the GC `p_phi` is conserved to round-off over a
round trip (rms 2e-13). Here the hybrid's statistic is the exit GC's mu over
the entry GC's mu, GC to GC, with no lift in the measurement. The full orbit
is measured with the same gyroaveraged first-order estimator at both ends:
mu over the gyroperiod after the exit over mu over the gyroperiod after the
start. One passage of the full orbit changes mu by -5.2e-6 / -5.1e-6 /
-3.9e-6 in the mean (standard error below 1.3e-7) with rms 1.5e-6 / 4.7e-6 /
8.4e-6 at v = 0.01 / 0.02 / 0.04; measured from entry instead of the start the
mean is -2.6e-6 / 0.5e-6 / 1.6e-6, so part of it is the estimator's drift on
the approach.

| entry + exit | v (rho) | `mu_after/mu_before - 1` mean, rms | KS mu D at N = 256 / 1024 / 4096 (p) | KS displacement D (p), N = 4096 |
|---|---|---|---|---|
| SIMPLE (GC `p_phi` + pinned) | 0.01 | 1.72e-4, 1.6e-4 | 0.64 / 0.61 / 0.62 (< 1e-48) | 0.026 (0.12) |
| | 0.02 | 3.63e-4, 3.3e-4 | 0.63 / 0.60 / 0.61 (< 1e-47) | 0.043 (1.2e-3) |
| | 0.04 | 7.70e-4, 7.2e-4 | 0.62 / 0.58 / 0.59 (< 1e-45) | 0.062 (2e-7) |
| averaged `p_phi` + pinned | 0.01 | 0.59e-4, 0.6e-4 | 0.63 / 0.60 / 0.61 (< 1e-47) | 0.011 (0.98) |
| | 0.02 | 1.30e-4, 1.3e-4 | 0.63 / 0.60 / 0.60 (< 1e-47) | 0.016 (0.68) |
| | 0.04 | 2.92e-4, 3.0e-4 | 0.62 / 0.57 / 0.58 (< 1e-45) | 0.020 (0.41) |
| averaged `p_phi` + plain | 0.01 | -0.93e-4, 1.0e-4 | 0.70 / 0.71 / 0.70 (< 1e-59) | 0.011 (0.98) |
| | 0.02 | -1.87e-4, 2.1e-4 | 0.68 / 0.72 / 0.72 (< 1e-57) | 0.016 (0.68) |
| | 0.04 | -3.98e-4, 4.0e-4 | 0.75 / 0.78 / 0.78 (< 1e-69) | 0.020 (0.41) |
| GC phase + plain (control) | 0.01 | -0.92e-4, 1.6e-4 | 0.67 / 0.67 / 0.67 (< 1e-53) | 0.011 (0.96) |
| | 0.02 | -1.92e-4, 3.3e-4 | 0.67 / 0.69 / 0.69 (< 1e-53) | 0.016 (0.63) |
| | 0.04 | -3.98e-4, 6.5e-4 | 0.71 / 0.69 / 0.69 (< 1e-61) | 0.030 (0.05) |

Standard errors of the means are below 1.1e-5. Transmitted fractions are as
in the previous table. In SIMPLE's rule the entry and exit shifts of mu do
not cancel but add: +2.8e-4 from the entry alone, +3.6e-4 GC to GC at
v = 0.02. Every rule changes mu by 1e-4 to 8e-4 per passage, proportional to
v, against a few 1e-6 for the full orbit, and the KS distance stays near 0.6
to 0.8 at every N. The GC-to-GC mu is that of one instant of the particle,
as the hybrid uses it, so its spread includes the gyrophase dependence of the
first-order exit, which the gyroaveraged full-orbit estimator does not have.
At this sheet the switching rules, not the sheet, set the change of mu per
passage.

### Gyroaveraged switches

`avgexit` continues the FO for one gyroperiod `2 pi/|B(X)|` (32 substeps),
averages the first-order X, `v_par` and `v_perp^2`, and takes
`mu = <v_perp^2>/(2 B(<X>))`, `v_par` from the energy with the sign of
`<v_par>`. `avgentry` lifts at the uniform phase, runs one gyroperiod of FO
through the same estimator, rescales the lift's mu by the ratio to the GC mu
and repeats until it agrees to 1e-9 (at most 8 trials); the trial orbits are
discarded. The `_p` variants also pin `p_phi`: `avgexit_p` takes `v_par` from
the entry GC `p_phi` at `<X>`, and `avgentry_p` shifts the lift position
along `grad psi` until the averaged GC has the GC `p_phi`. The full-orbit
change over the same passage, measured from entry with the gyroaveraged
estimator, is -2.6e-6 / 0.5e-6 / 1.6e-6 in the mean (standard error 0.02e-6 /
0.08e-6 / 0.15e-6) with rms 1.6e-6 / 5.4e-6 / 9.7e-6 at v = 0.01 / 0.02 /
0.04. KS is against that; against exit over start it is larger (0.47 to
0.70 for `avgentry` + `avgexit`).

| switches | v (rho) | GC to GC `mu_after/mu_before - 1` mean, rms | KS D at N = 256 / 1024 / 4096 (p at 4096) | `p_phi` jump rms (mean) | extra gyroperiods per entry + exit |
|---|---|---|---|---|---|
| SIMPLE | 0.01 | 1.72e-4, 1.6e-4 | 0.64 / 0.61 / 0.62 (< 1e-300) | 6e-15 | 0 |
| | 0.02 | 3.63e-4, 3.3e-4 | 0.63 / 0.60 / 0.61 (< 1e-300) | 2e-13 | 0 |
| | 0.04 | 7.70e-4, 7.2e-4 | 0.62 / 0.58 / 0.59 (< 1e-300) | 2e-13 | 0 |
| `avgentry` + `avgexit` | 0.01 | -2.1e-6, 1.9e-6 | 0.26 / 0.26 / 0.26 (1e-123) | 1.0e-6 (-0.9e-6) | 4.0 |
| | 0.02 | -0.2e-6, 4.4e-6 | 0.31 / 0.26 / 0.25 (4e-116) | 4.1e-6 (-3.6e-6) | 4.6 |
| | 0.04 | 3.2e-6, 9.2e-6 | 0.13 / 0.09 / 0.13 (4e-28) | 1.7e-5 (-1.5e-5) | 4.8 |
| `avgentry_p` + `avgexit_p` | 0.01 | 6.5e-5, 5.5e-5 | 0.63 / 0.60 / 0.61 (< 1e-300) | 4e-9 | 6.9 |
| | 0.02 | 1.42e-4, 1.2e-4 | 0.63 / 0.60 / 0.60 (< 1e-300) | 3e-8 | 8.7 |
| | 0.04 | 3.22e-4, 2.8e-4 | 0.62 / 0.57 / 0.58 (< 1e-300) | 4e-7 | 9.0 |

Gyroaveraging both switches without the `p_phi` constraint brings the mean
and rms of the change per passage down to the full-orbit level, about 80
times below SIMPLE's rule in rms. It is still distinguishable: the means
differ by 0.5e-6 to 1.7e-6 (7 to 12 standard errors), and the KS distance
does not shrink with N. The `p_phi` round-trip jump is 1e-6 to 2e-5 and mostly systematic.
Adding the `p_phi` constraint conserves `p_phi` to 4e-9 to 4e-7 but brings
back a bias of about half SIMPLE's. Its entry iteration often fails to
converge in 8 trials (73 of 4096 markers at v = 0.02, nearly all at
v = 0.04).

`avgexit_x` keeps the averaged exit's mu and instead moves the averaged GC
along `grad psi` until its `p_phi` equals the entry GC value, with `v_par`
from the energy at the new position (3 iterations). With `avgentry` the mu
statistics are unchanged, as constructed. The round-trip `p_phi` jump loses
its systematic part (rms 0.7e-6 / 2.7e-6 / 1.1e-5, means below 1.5e-7, against
1.0e-6 / 4.1e-6 / 1.7e-5 with means -0.9e-6 / -3.6e-6 / -1.5e-5 for
`avgexit`). The remaining jump comes from the unconstrained averaged entry.
The mean change of `psi` at the GC per passage (entry GC to exit GC; for the
full orbit, gyroaveraged X after entry to after the exit) agrees with the
full orbit within one standard error for SIMPLE's rule, `avgexit` and
`avgexit_x` alike: 4.61e-2 against 4.65e-2 (+-0.06e-2) at v = 0.02. At this
precision it does not separate the rules; its spread (rms 4.1e-2) is set by
the mix of transmitted and reflected markers.

## Real SPECTRE interfaces

`tools/spectre_interface_jumps.f90` (built against a SIMPLE build, see its
header) evaluates the true field at `rho_g = k -+ 1e-10` on a 128 x 64 grid:

| fixture, interface | direction jump max (rms) | `\|B+\|/\|B-\| - 1` |
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

The removed `test_fo_symplectic` checked an analytic tangent for the landing
construction on the crossings it completed. That did not validate grazing
coverage or physical long-time accuracy. Its source and the old stationary
comparison remain in history. Current smooth-layer methods retain their original
Hamiltonian equations; no exact symplectic claim for a terminated floating-point
solve or for the production hybrid follows from these historical tests.

## SIMPLE on NCSX: loss fractions against full orbit everywhere

`ncsx_freeboundary.h5` (Mvol = 3), 3.5 keV alphas (`facE_al = 1000`), 64
markers started next to an interface, `trace_time = 2e-3`. The original
study used 32 threads; the corrected-event hybrid row below uses eight threads
per run, with the two interfaces run concurrently, deterministic seed 12345,
and `relerr = 1e-12`. Other rows retain historical results and timings. The
reference is full orbit everywhere: crossing_level 3, RK45 at tolerance 1e-9,
`spectre_fo_entry_larmor = 1e6`.

| model | s = 0.95 (interface 1): lost, wall s | s = 1.95 (interface 2): lost, wall s |
|---|---|---|
| full orbit everywhere, RK45 (reference) | 6, 920 | 39, 565 |
| level 3, RK45 (corrected events, 8 threads) | 5, 3343 | 39, 1848 |
| level 3, Boris (layer 0.2, rotation 0.0125) | 5, 1399 | 48, 480 |
| level 1, symplectic GC | 7, 15 | 29, 25 |
| level 1, RK45 GC | 1, 17 | 7, 30 |

Counting noise is about +-2 markers at interface 1 and +-4 at interface 2. At
the strong outer interface both guiding-centre maps confine particles
spuriously (RK45 GC: 32 of the 39 reference losses are missed). Level 3 with
RK45 matches the reference at both interfaces; Boris overestimates the losses
at interface 2 by about 2 sigma. Individual marker fates differ even between
full-orbit runs (chaotic orbits near the interface); only the fractions
compare. The markers stay near the interfaces, so the hybrid spends substantial
time in full orbit. Corrected and historical timings have different thread counts
and are not a direct speed comparison.

For strong interfaces use `crossing_level = 3` with `integmode > 0`;
sharp-field RK45 with exact crossings is the full-orbit integrator. The GC direction-jump ordering and
operational warning are documented in
[`DOC/spectre-interface-crossing.md`](../../DOC/spectre-interface-crossing.md#warning-guiding-centre-maps-at-strong-interface-jumps).
