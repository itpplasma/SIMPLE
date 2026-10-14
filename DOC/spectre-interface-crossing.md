# SPECTRE orbit interface status

Multi-volume SPECTRE orbit tracing requires `integmode > 0` and
`crossing_level = 3` (the default). Levels 0/1 were removed without
compatibility aliases; requesting either is a configuration error. Single-volume
RK45 guiding-centre tracing and sampling-only initialization remain supported.

The conserved-magnetic-moment guiding-centre crossing maps were removed because
a sharp sheet changes the physical magnetic moment by a gyrophase-dependent
amount and the maps produce spurious confinement. The
[rejection record](spectre-rejected-interface-models.md) retains the mathematical
arguments, symbolic verifier, recorded measurements, and code provenance.

## Full orbit near interfaces (`crossing_level = 3`)

With `crossing_level = 3` and `integmode > 0`, the symplectic guiding-centre
path switches to the full orbit (`spectre_fo_hybrid`) when the guiding centre
comes within `spectre_fo_entry_larmor` radial gyroradii of an interior
interface, and back when it is one gyroradius farther out. RK45 uses the sharp
per-volume field and localizes interface crossings. Every trial stage's
geometric owner is checked, including shortened
bracketing and landing trials and ownership changes in both directions. A stage
that crosses and returns before the endpoint causes rejection and step halving;
a changed endpoint is localized on the first sheet leaving the initial
owner (`owner - 1` inward, `owner` outward), including trials that span a thin
intermediate volume. The remaining time is integrated in the new owner.
Starts on a sheet use the outgoing normal velocity; a later return is
bracketed at positive time. Repeated grazing reduces/retries and eventually
reports `SPECTRE_FO_STEP_LIMIT`, rather than committing an unlocalized step.
Away from the polar axis, RK45 advances logical position and Cartesian velocity.
The sheet is directly `rho_g = integer`, so stages and landing trials require no
Cartesian-to-logical coordinate inversion. Within `rho_g < 0.05` in volume 1,
Cartesian position provides a regular axis chart; this recovery region still
uses spline coordinate inversions.

`spectre_fo_field` constructs a separate quintic spline of R, Z, A_theta and
A_zeta for each volume. The embedding derivatives give the basis and its dual;
derivatives of A give B with libneo's Jacobian convention. Shared sheet geometry
uses the same interface samples, while A retains the two one-sided volume
gauges. Raw spectral geometry and field evaluation occur during construction,
not during FO integration. The grid is `spectre_fo_grid = 32,64,32` points in
rho per volume, theta and zeta per field period; axisymmetric fields use seven
identical toroidal planes. This interpolates the equilibrium: reducing the RK
tolerance does not remove field interpolation error. Refine this grid as well
as the RK tolerance for loss convergence.

GC reconstruction retains its nonlinear coordinate solve, now using only the
spline geometry and a local displacement prediction. Particle and GC field
samples are reused in the existing exit test. The switch transform and exit
criteria are unchanged; no whole-hybrid symplectic claim follows from this
performance change.

RK45 also caps each trial at `h |B| / ro0_bar <= 0.1` rad, about 1/63 of a
gyro-period. This cap is rechecked against each trial's initial field, including
zero-time owner changes. Stage checks cannot exclude arbitrarily shallow
excursions between samples. In the local uniform-field, planar-sheet approximation,
a circular cross-and-return within this entire step has depth at most
`rho_L (1 - cos(0.05))`, approximately `1.25e-3 rho_L`; interface curvature
and field variation limit this estimate. This bounds the geometric scale of
the remaining unresolved-grazing risk, rather than claiming that the RK error
estimate detects a force jump. The optional `h_scale` reduces this cap for
convergence studies. Failed trial stages halve the step without reading an
uncomputed error estimate. A location failure near the outer edge is numerical
failure, not by itself a physical loss; loss requires a located terminal point,
whose state and elapsed substep time are retained.

Each volume's vector potential has its own gauge (it vanishes on the volume's
inner boundary). Both the axisymmetric entry residual and exit mapping convert
momenta with the same summed interface jumps of `A_zeta`. A trial particle in
another volume is translated back to the entry gauge before comparison. If
there is no gyrophase root, entry accepts the closest sampled phase only when
`abs(delta P_zeta) <= 1e-3 max(abs(P_zeta_entry), 1)` in internal momentum
units; otherwise it reports `SPECTRE_FO_GYROPHASE_FAIL`. The signed accepted
defect is exposed as `spectre_fo_state_t%entry_pzeta_defect`. No gyroradius-squared
error estimate is asserted for this fallback. RK45 can drift in energy.
The unchanged long-trace bounds remain shell defect below `1e-4` and relative
momentum spread below `1e-5`. The integration tests use public `simple_api`
initialization and require every marker to complete and enter FO, with
round-trip witnesses and invariant samples in every visited volume. A focused
near-tangent test compares against independent geometric-field RK4 and the
axisymmetric momentum
invariant. The prototype in `concept/symplectic_sheet` documents the analysis.

The GC crossing-event log was removed with the rejected maps. For level-3
loss accounting use `times_lost.dat` (column 2: `0 < t < trace_time` is a loss);
FO loss counters cover only FO phases.
FO configuration is applied in shared SPECTRE field initialization, before
parallel tracing, identically for `simple.x` and `simple_api`/`pysimple`.

## Warning: guiding-centre maps at strong interface jumps

Guiding-centre (GC) crossing maps (`crossing_level = 0` or `1`) are not
valid where the field direction or `|B|` jumps strongly. Conservation of
energy and `mu` by a map does not establish its physical validity. For a
local sheet `B = B(z)(cos alpha(z), sin alpha(z), 0)`,
`h . curl h = -alpha'(z)`, so the GC parallel Jacobian is

    B*_par = B(1 - rho_par alpha'),   rho_par = m c v_par / (e B).

At a finite direction jump compressed into width `w`, `alpha' ~ [[alpha]]/w`;
for nonzero parallel gyroradius the GC ordering fails as `w -> 0`, and the
GC equations have no regular direction-jump limit. A strong magnitude jump
also violates the slow-field-variation ordering and can scatter `mu`. The
full-orbit equations remain well defined through the sharp sheet.

The NCSX free-boundary outer interface has direction jumps up to 0.32 rad
and magnitude jumps of -25 to +38 %. Recorded losses for 64 markers at
3.5 keV, traced for 2 ms from each interface, are:

| model | interface 1 (`rho_g = 0.95`) | interface 2 (`rho_g = 1.95`) |
|---|---:|---:|
| full orbit everywhere, RK45 (reference) | 6 | 39 |
| level 3, RK45 (corrected events) | 5 | 39 |
| level 3, Boris (layer 0.2, rotation 0.0125) | 5 | 48 |
| level 1, symplectic GC | 7 | 29 |
| level 1, RK45 GC | 1 | 7 |

At the strong outer interface the GC maps give spurious confinement;
level 3 with RK45 agrees with the reference loss fractions within counting
uncertainty (about 2 and 4 markers at the two interfaces). Boris on the
0.2-gyroradius layer overestimates the outer-interface losses. These are
ensemble comparisons; individual fates differ even between full-orbit runs.

Use `crossing_level = 3`, `integmode > 0`, and
RK45 at such interfaces. This resolves gyromotion
and costs far more than the GC maps. The corrected-event hybrid RK45 reruns
took 3343/1848 s at interfaces 1/2,
using two concurrent runs with eight OpenMP threads each, deterministic seed
12345, and `relerr = 1e-12`. Both loss counts agree with the retained full-orbit
references within the stated counting uncertainty. The earlier hybrid counts
were 7/41 with wall times 929/409 s; the retained reference and GC timings were
920/565 s and 15/25 s. These historical timings use different thread counts,
so they are not a direct speed comparison with the corrected runs. Setup and
full timing evidence are retained in
[`concept/symplectic_sheet/README.md`](../concept/symplectic_sheet/README.md#simple-on-ncsx-loss-fractions-against-full-orbit-everywhere).

`mu` conservation through a zero-width sheet is a modeling convention: the
Buechner-Zelenyi adiabaticity parameter vanishes for an ideal
discontinuity, and real crossings scatter `mu` (quasi-adiabatic regime).
The MRxMHD sheet stands in for a thin finite relaxed layer; full-orbit
comparisons quantify the scatter and bound the domain of validity of the
conserved-mu map.

## Removal of the smoothed-layer Boris option

The SPECTRE Boris integrator and its integrator/layer/rotation configuration
keys were removed without aliases. A finite smoothed layer changes the magnetic
field and therefore the physics; it does not integrate exactly the sharp-sheet
Hamiltonian. RK45 with exact crossings reproduces the retained full-orbit-everywhere
reference loss fractions within counting uncertainty (table above).
See [the symbolic rejection record](spectre-rejected-interface-models.md#t2-a-boris-smoothing-layer-changes-the-target-problem)
for the model distinction. Historical measurements above are evidence, not
available configurations. The last revision containing the removed code is
[`afec0c9`](https://github.com/itpplasma/SIMPLE/commit/afec0c9ccc5a90350ec71755a611450c23623083).

`spectre_fo_entry_larmor`, radial-gradient evaluation and exit-distance state
remain necessary for GC/RK45 switching. The fixed 0.1-rad RK45 gyro-angle bound
and `h_scale` still bound trial excursions; neither is a Boris rotation option.
The standalone study's analytic layer width also remains for its non-Boris
Gauss/Poincare convergence experiments.
