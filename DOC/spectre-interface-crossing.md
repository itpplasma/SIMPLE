# Guiding-center crossing of MRxMHD interfaces

Staging draft for `DOC/spectre-interface-crossing.md` in SIMPLE (issue #440).
Every formula here is machine-verified in `ca/06_crossing_map.wl` of
plasma/proj/spectre-orbits or derived in its `tex/spectre-orbits.pdf`,
section 6.

## Setting

A SPECTRE interface at integer `rho_g` is a flux surface from both sides
(`B^rho = 0`), carrying jumps `[[B]]` (from `[[p + B^2/2]] = 0`) and
`[[h]]` (tangential rotation of `h = B/|B|`). Guiding centers reach it by
drift only. The per-volume equations of motion contain no sheet force, so
the crossing is an explicit jump map applied between the exact-landing
event (arrival state on the surface) and re-canonicalization in the target
volume. The map acts on physical variables `z = (rho_g, theta, zeta,
v_par)` at fixed magnetic moment `mu`; canonical momenta are rebuilt from
the target volume's own gauge afterwards, so no gauge matching is needed.

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
fixture the poloidal variation is weak, so most crossings take the Level-0
rescale and a minority refract with kicks up to `~0.2` rad.

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

## Physical caveat

`mu` conservation through a zero-width sheet is a modeling convention: the
Buechner-Zelenyi adiabaticity parameter vanishes for an ideal
discontinuity, and real crossings scatter `mu` (quasi-adiabatic regime).
The MRxMHD sheet stands in for a thin finite relaxed layer; the full-orbit
Boris pusher quantifies the scatter (SIMPLE#442) and bounds the domain of
validity of the conserved-mu map.
