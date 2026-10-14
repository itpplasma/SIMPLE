# SPECTRE orbit interface status

Multi-volume SPECTRE orbit tracing (`Mvol > 1`) is disabled on main until the
full-orbit interface treatment of PR #530 lands. Configuration fails with:

> no valid interface crossing on this version; see DOC/spectre-rejected-interface-models.md

The conserved-magnetic-moment guiding-centre crossing maps were removed because
a sharp sheet changes the physical magnetic moment by a gyrophase-dependent
amount and the maps produce spurious confinement. The
[rejection record](spectre-rejected-interface-models.md) retains the mathematical
arguments, symbolic verifier, recorded measurements, and code provenance.

## State and invariants (SIMPLE RK45 path)

The tracer state is `z = (rho_g, theta, zeta, p, lambda_p)` with `p` the
momentum modulus (normalized to `sqrt(2 T/m)`) and `lambda_p = v_par/p` the
pitch cosine, so `v_par = p lambda_p`. The stored perpendicular invariant is

    perp_inv = p^2 (1 - lambda_p^2)/|B| = v_perp^2/|B|,   mu = perp_inv/2,

with `mu` held fixed across the sheet. Because `H = v_par^2/2 + mu |B| = p^2/2`,
the momentum modulus `p` is conserved by both map levels, and the pitch after
the map is `lambda_p' = v_par'/p` with `|v_par'| <= p`.

## Level 0 (implemented by SIMPLE#443)

Hold `(theta, zeta, mu)`; flip charts (`s: +1 -> -1`); update

    v_par'^2 = v_par^2 - 2 mu (B_plus - B_minus),

with both `B` evaluated at the same interface point from the respective
volumes. Sign of `v_par` carries over. If the radicand is negative the
crossing is forbidden: the particle stays in its volume (sheet
reflection) and the event is logged. Implementation regularization: the
sign of `v_par` is reversed on reflection; with a strictly unchanged
state the persistent outward drift re-triggers the event indefinitely,
while the sign flip conserves `H` and `mu` exactly, perturbs only the
small-`v_par` reflection class, and lets the orbit leave the sheet. The
Level-1 tangential kick resolves the sheet-skimming dynamics properly.
Energy `H = v_par^2/2 + mu B` is conserved exactly in both branches.

## Level 1 (the thin-layer limit; SIMPLE#440)

The regularized sheet dynamics converges, as the layer width goes to zero,
to an impulse along the Hamiltonian vector field of the interface function
`s` under the guiding-center Poisson bracket:

    Delta z = lambda {z, s},        H(z + Delta z)|_plus = H(z)|_minus,

one scalar `lambda` per crossing. Component formulas in a chart with
covariant `h_i`, Jacobian `sqrtg`, and `Bstar_par = B + (m c v_par / e)
h . curl h` (SIMPLE normalization: replace `m c/e` by its internal
constant):

    {theta, s} = - (c/e) h_zeta / (sqrtg Bstar_par)
    {zeta,  s} = + (c/e) h_theta / (sqrtg Bstar_par)
    {v_par, s} = (c v_par / e) (curl h)^rho / Bstar_par
    {s, s} = 0            (the map stays on the surface)

The tangential kick direction `(-h_zeta, +h_theta)` is perpendicular to B
within the surface: the sheet drift. Evaluate the bracket components in
the symmetrized state (average of the incoming state and the kicked
state, solved together with lambda by the same Newton loop) so the
discrete map is symplectic, not only its continuum limit. The reflection
branch takes the other root of the energy condition: same-side exit with
the tangential kick applied, energy exact.

## Level 1 in SIMPLE units and code (SIMPLE#440)

SIMPLE traces the per-volume drift with `velo_can` in internal Gaussian
units. Its guiding-center coefficient is the reference Larmor radius
`ro0 = m c v0 / (e B_ref)` (module `parmot_mod`), and its parallel Jacobian
is `hpstar = 1 + p_par (ro0/|B|) (h . curl h) = Bstar_par/|B|` with
`p_par = v_par` in the nonrelativistic limit SIMPLE runs. The crossing map
must sit on the same drift, so the generator uses the same `ro0` and the
same `Bstar_par`, not a separately reconstructed `c m/e`:

    Bstar_par = |B| + ro0 v_par (h . curl h),   h . curl h = hcovar . hcurl.

In coordinate order `(1, 2, 3) = (rho_g, theta, zeta)` the bracket
components `X = {z, rho_g}` that `interface_crossing.f90` applies are

    X_theta = - ro0 hcovar_zeta  / (sqrtg Bstar_par)
    X_zeta  = + ro0 hcovar_theta / (sqrtg Bstar_par)
    X_vpar  = + ro0 v_par (curl h)^rho / Bstar_par

with `hcovar`, `sqrtg`, and `hcurl` taken from `magfie`. These are the
`(c/e)`-scaled components of the previous section with `c/e -> ro0`: the
constant that scales the sheet drift is exactly the one `velo_can` uses for
its own `a_c(i) = hcurl(i) ro0/|B|` and for `hpstar`.

**Crossing vs mirror.** The energy criterion at the landing point,
`radicand0 = v_par^2 - 2 mu [[B]]` with `[[B]] = B_target - B_home` both at
the landing angles, decides the branch. It is the Level-0 criterion, so the
event split matches Level 0 and Level 1 only refines the crossing state. A
forbidden crossing (`radicand0 < 0`) reflects as the Level-0 mirror
(`v_par -> -v_par`, `mu`, `|B|` unchanged, stays home), which keeps the
trapped sheet-skimming class identical to Level 0.

**Refraction solve.** For a crossing the impulse `Delta z = lambda_k X`
acts on `(theta, zeta, v_par)`; `rho_g` stays on the interface and then
switches volume with the Level-0 nudge. The single scalar `lambda_k` solves
the exact energy condition

    H_target(z + lambda_k X) = H_home(z),   H = v_par^2/2 + mu |B|,

by a damped Newton on the residual `F(lambda_k) = 0.5 v_par'^2 +
mu B_target(kick) - H_home`, with `v_par' = v_par + lambda_k X_vpar` and the
target `|B|` at the kicked angles. Backtracking halves the step whenever it
fails to reduce `|F|`, because `B_target` varies poloidally and `F` is not
monotone in the kick. The generator is refreshed each residual at the
midpoint state `z + (lambda_k/2) X` from both-side field quantities
averaged, so the discrete map is symplectic, not only its continuum limit.
`mu` is held fixed and `p^2 = 2 H` is conserved, so the refracted crossing
keeps `v_par` to drift order and lands the angles on an equal-`|B|` target
point, where Level 0 instead rescales `v_par` at fixed angles.

**Thin-layer bound and fallback.** The equal-`|B|` point sits a tangential
distance `~ [[B]] / |grad_perp B|` away, that is `~ mu [[B]] / v_drift_normal`
in drift time. Where the poloidal `|B|` gradient along the sheet drift is
weak this exceeds the thin-layer picture, so the solve is capped at a kick
`|lambda_k X_tang| <= max_kick` (0.3 rad); beyond it the map falls back to
the energy-exact Level-0 rescale (`v_par' = sign(v_par) sqrt(radicand0)`, no
kick). Every energetically allowed marker therefore crosses, with the
tangential refraction applied where it is small. On the two-volume tok2vol
fixture most crossings refract, with kicks up to `~0.3` rad.

## Impulse-map validity domain (finding F5)

The single-impulse form is the thin-layer limit only while the energy
balance closes within a small tangential displacement:
`mu [[B]] < v_par * (tangential B variation over the kick scale)`.
Outside that domain the ideal sheet supplies no parallel mirror channel
(`B . n = 0` holds inside a tangential-discontinuity layer, so
`h . grad B` has no normal contribution), and the true limit is a long
tangential slide along the sheet toward an equal-|B| point, which leaves
the impulse map's ordering. Implementation: crossings whose refraction
root exceeds a 0.3 rad tangential kick fall back to the energy-exact
Level-0 rescale (the adiabatic-mirror convention). Mirror-vs-slide in this
regime is an open modeling question to be adjudicated by the full-orbit
comparison (SIMPLE#442); both conventions conserve H and mu exactly.

## Mixed-sheet ambiguity (finding F4) and the adopted convention

For pure `[[B]]` or pure `[[h]]` sheets the thin-layer limit is unique
(the kick integrands telescope). For simultaneous jumps the tangential
kick along `h` depends on the internal ordering of compression and
rotation; force balance does not fix it. Adopted convention:
**proportional mixing**: both jumps progress with a single interior
variable, which is what the single-impulse form above realizes. The
validation suite (SIMPLE#442) must quantify the spread between the two
extreme orderings (rotation-first vs compression-first, available from
the regularized slab integrator in `python/fig_crossing_convergence.py`)
and report it as model uncertainty next to the measured full-orbit
mu-scatter.

## Slab test oracles (exact, for the two-volume test equilibria)

Driven sheared slab (`B = B(z)(cos a(z), sin a(z), 0)`, force `F xhat`):

    v_par' = v_par sin(a_minus) / sin(a_plus)
    F Dx   = mu [[B]] + (v_par'^2 - v_par^2)/2
    Dy     = (v_par^2 sin^2(a_minus)/F) (cot a_plus - cot a_minus)
             - (mu/F) cot(a) [[B]]        (second term: constant-a sheets)

Energy `E = v_par^2/2 + mu B - F x` is conserved exactly by the model
equations and by the map. Level 0 in this model misses `F Dx` per
crossing: a width-independent defect (the "no step refinement fixes a
missing force" witness).

## Algorithm (per crossing)

1. Event: integrator lands exactly on the interface (dense-output root
   for RK45; Newton on the substep length for the implicit symplectic
   schemes (same scheme, shorter step); polynomial basis extension makes
   Newton overshoot past the surface safe).
2. Evaluate both-side `B`, `h`, metric at the landing point; solve the
   energy condition for lambda (Level 1) or the `v_par` rescale
   (Level 0); reflection branch on no real root.
3. Apply the kick in physical variables; switch chart.
4. Re-canonicalize in the target volume (same code path as orbit start);
   resume with the standard step size.

Crossings are drift-rare (once per radial volume transit), so per-crossing
cost is negligible.

### Symplectic path specifics (SIMPLE#441, `spectre_sympl_orbit.f90`)

The accepted step that leaves the home volume is never committed: the
substep length `h*` with `|rho_g(h*) - k| < 1e-10` is solved by Illinois
false position, each evaluation one full implicit step of the same scheme
from the same pre-step state, so the landed state lies on the exact
discrete orbit of a symplectic map. Per-volume splines are evaluated with
a home-volume lock and a C1 radial extension from an `EDGE_BAND` inside
the volume edge (`field_can_spectre`): `A` and `Bmod` continue linearly,
the radial slope of `h_theta`, `h_phi` blends to zero, which keeps
`dpth/dr = sqrtg*Bmod/(ro0*h_phi) > 0` in the whole sheet-adjacent zone
where the raw guiding-center 1-form degenerates (`Bstar_par -> 0`, the
canonical analogue of the RK45 drift band). Where the landing root falls
in a step-length regime with a multi-branched implicit solve (grazing
events) the substep is halved and re-solved; committed substeps that jump
more than half a volume width in `rho_g`, more than O(1) in energy, or
flip theta by pi away from the axis are unconverged-Newton artifacts and
are likewise halved. After `apply_crossing` the integrator is
re-initialized from the mapped physical state in the target volume's own
gauge (the orbit-start code path, so no gauge matching), and the remaining
microstep budget is completed there, preserving the scheme's order under
step halving. Only one-step schemes are admitted (`ntau = 1` asserted); a
persistent non-convergence terminates the orbit through the `CROSS_STOP`
fallback, which the tok2vol tests assert does not occur.

## Full orbit near interfaces (`crossing_level = 3`)

With `crossing_level = 3` and `integmode > 0`, the symplectic guiding-centre
path switches to the full orbit (`spectre_fo_hybrid`) when the guiding centre
comes within `spectre_fo_entry_larmor` radial gyroradii of an interior
interface, and back when it is one gyroradius farther out. The default
`spectre_fo_integrator = 'rk45'` uses the sharp per-volume field and localizes
detected interface crossings; `'boris'` remains an option for a resolved smooth layer.

A fixed step that does not resolve the force jump at the sheet makes a
zero-mean impulse error `~ h [[F]]` per crossing that random-walks the
invariants. `spectre_fo_layer_larmor` replaces each interior interface by a
smooth layer of that many gyroradii for the full orbit, and
`spectre_fo_rotation` bounds the Boris step. On tok2vol at 3.5 keV
(`test_spectre_fo_interfaces`) the relative spread of the exactly conserved
canonical toroidal momentum over all full-orbit phases is 14 % with the
sharp field at the default bound 0.05, 3.2 % at 0.0125, 4.0 % with a layer of
0.2 gyroradii at 0.05, and 0.93 % with the layer at 0.0125.
`spectre_fo_integrator = 'rk45'` integrates the full orbit adaptively on the
sharp field. Every trial stage's geometric owner is checked, including shortened
bracketing and landing trials and ownership changes in both directions. A stage
that crosses and returns before the endpoint causes rejection and step halving;
a changed endpoint is localized on the first sheet leaving the initial
owner (`owner - 1` inward, `owner` outward), including trials that span a thin
intermediate volume. The remaining time is integrated in the new owner.
Starts on a sheet use the outgoing normal velocity; a later return is
bracketed at positive time. Repeated grazing reduces/retries and eventually
reports `SPECTRE_FO_STEP_LIMIT`, rather than committing an unlocalized step.
The Cartesian radial residual refines near-sheet inverse coordinates that
would otherwise snap both bracket ends onto the sheet.

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
error estimate is asserted for this fallback. RK45 can drift in energy;
Boris preserves speed. The existing long-trace bounds remain shell defect
below `1e-4` and relative momentum spread below `1e-5` for RK45, and `1e-9`
and `2e-2` for Boris at layer 0.2 and rotation 0.0125. The integration tests use
the public `simple_api` initialization, exercising default RK45 and nondefault
Boris configuration, and require every marker to complete and enter FO, with
round-trip witnesses and invariant samples in every visited volume. A focused
near-tangent test compares against independent geometric-field RK4 and the
axisymmetric momentum
invariant. The prototype in `concept/symplectic_sheet` documents the analysis.

FO entries, sheet transits, exits, and FO losses are excluded from
`spectre_crossing_events.dat`, which records GC crossing-map events. For level-3
loss accounting use `times_lost.dat` (column 2: `0 < t < trace_time` is a loss)
and the FO loss counters, rather than counting type-3 crossing-log events.
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
`spectre_fo_integrator = 'rk45'` at such interfaces. This resolves gyromotion
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
