# SPECTRE interface crossing

How SIMPLE traces markers through the current sheets of a multi-volume
SPECTRE equilibrium: when it leaves the guiding-centre (GC) model, how the
full orbit (FO) is integrated and mapped back, and which settings control
accuracy. The model overview is in `docs/model/extensions.rst`; why GC maps
across a sheet are not used is derived in
`DOC/spectre-rejected-interface-models.md`.

## Configuration

Multi-volume tracing requires `integmode > 0` and `crossing_level = 3`
(the only accepted value); otherwise initialization stops with a message
pointing to the rejected-model note. Single-volume RK45 GC tracing and
sampling-only runs are unrestricted.

| Key | Default | Meaning |
|---|---|---|
| `spectre_fo_entry_larmor` | `1.5` | enter FO closer than this many radial gyroradii to a sheet (and to the outer boundary with `spectre_outer_loss = 'particle'`) |
| `spectre_fo_jump_tol` | `1d-3` | interfaces whose largest jump of $\lvert B\rvert/B$ or of $\mathbf b$ is smaller are crossed by the GC directly |
| `spectre_outer_loss` | `'guiding_centre'` | outer-boundary loss: when the guiding centre crosses it, as for VMEC and GVEC fields, or `'particle'` when the full orbit touches it |
| `spectre_fo_tolerance` | `1d-9` | relative tolerance of the FO RK45 |
| `spectre_fo_grid` | `32, 64, 32` | FO spline points per volume in $\rho_g$, $\vartheta$, $\zeta$ (per field period) |

## Switching rule

In volume $k$ the GC is advanced by the symplectic integrator in that
volume's Meiss chart (`src/spectre_sympl_orbit.f90`). After each microstep,
`interface_proximity` takes the nearest integer $k' = \mathrm{nint}(\rho_g)$
(an interior sheet, or with `spectre_outer_loss = 'particle'` also the outer
boundary $M_\text{vol}$) and the radial gyroradius

$$
\rho_r = \frac{\bar\rho_0\,\bar v_\perp}{B}\,|\nabla\rho_g|,
$$

and enters FO if $|\rho_g - k'| < \ell_e\rho_r$, unless $k'$ is a smooth
interface (below), with
$\ell_e = $ `spectre_fo_entry_larmor`. It leaves FO when the reconstructed
GC satisfies $|\rho_g - k'| > (\ell_e + 1)\rho_r$ and the field at the
particle and at the GC agree, $|B_p - B_{gc}| \le 0.05\,B_{gc}$ and
$|\mathbf b_p - \mathbf b_{gc}| \le 0.05$. The exit test runs once per radian
of gyration and at the end of each call.

## Smooth interfaces

At configuration SIMPLE samples every interior interface on a 64 by 32 grid
of one field period and records the largest $|B_+ - B_-|/B_-$ and
$|\mathbf b_+ - \mathbf b_-|$ (printed as `SPECTRE interface field jumps`).
An interface whose jump is below `spectre_fo_jump_tol` starts no FO. When
the GC lands on it, it continues in the neighbouring volume at the same
position, speed and pitch, re-canonicalized in that volume's chart. The
reset is exact for a continuous field (T1b in
`DOC/spectre-rejected-interface-models.md`); its error per crossing is of
the order of the jump. Interfaces of SPECTRE representations of smooth
equilibria (eight volumes on the VMEC surfaces of the Landreman-Paul QA
case: $10^{-4}$ to $3\times10^{-4}$ at Fourier order 10, up to
$9\times10^{-4}$ at order 6) fall below the default, current sheets
(tok2vol $4.5\times10^{-2}$, NCSX free boundary $7.7\times10^{-2}$ and
$0.38$) above it.

## GC to FO

`spectre_fo_enter` (`src/spectre_fo_hybrid.f90`) places the particle at

$$
\mathbf x = \mathbf X + \frac{\bar\rho_0}{B}\,\mathbf b\times\bar{\mathbf v}_\perp(\alpha),
\qquad \bar{\mathbf v} = \bar v_\parallel\mathbf b + \bar{\mathbf v}_\perp(\alpha).
$$

The gyrophase $\alpha$ is drawn uniformly, as a hash of the guiding-centre
state rounded to $10^{-9}$, so it is reproducible for any thread count.

A lift of this kind is correct only at zeroth order in the Larmor radius:
the gyroaverage of the lifted particle's magnetic moment differs from the
GC value $J_\perp$ by $O(\rho_*)$, with a mean over gyrophase that does not
vanish near a current sheet. The entry therefore calibrates the lift. A
scratch copy of the lifted particle is integrated for one gyroperiod
$2\pi\bar\rho_0/B(\mathbf X)$, and the gyroaveraged moment
$\langle J_\perp\rangle$ is computed with the estimator of the exit (below).
$\bar v_\perp$ is rescaled by $\sqrt{J_\perp/\langle J_\perp\rangle}$ at fixed
speed and sign of $\bar v_\parallel$, and the particle is lifted again at the
same $\alpha$, until $|\langle J_\perp\rangle/J_\perp - 1| < 10^{-9}$, with at most
five lifts. The scratch orbits are discarded and the FO starts from the last
lift; no time passes. The particle's $\bar p_\zeta$ is not constrained; its
difference from the GC value is kept in
`spectre_fo_state_t%entry_pzeta_defect` in axisymmetry.

If a lift lies outside the outer boundary, up to eight phases in all are
tried in golden-angle steps (with `spectre_outer_loss = 'particle'` the
marker is lost instead). If the iteration does not converge or its scratch
orbit fails, the entry falls back to the instantaneous lift: in an
axisymmetric field the particle's toroidal canonical momentum
$\bar p_\zeta = \bar{\mathbf v}\cdot\mathbf e_\zeta + A_\zeta/\bar\rho_0$ is
made equal to the GC value by a secant solve for $\bar v_\parallel$ at fixed
speed, with the same phase retries and `SPECTRE_FO_GYROPHASE_FAIL` when no
phase works. `simple.x` prints the number of such instantaneous entries and
exits after tracing.

## FO integration

The FO is advanced by RK45 in logical position $(\rho_g,\vartheta,\zeta)$
and Cartesian velocity in the sharp field of the owning volume. Each trial
step is capped at a gyro-angle $h B/\bar\rho_0 \le 0.1$, about 1/63 of a
gyro-period. In a uniform field near a planar sheet, a cross-and-return
excursion hidden inside such a step has depth at most
$\rho_L(1 - \cos 0.05) \approx 1.25\times10^{-3}\rho_L$.

Every stage's owning volume is checked. A stage that crosses a sheet and
returns before the step end rejects the step and halves it. A step whose
end lies in another volume is shortened until it lands on the first sheet
left by the initial owner ($k-1$ inward, $k$ outward); the rest of the step
continues in the new owner. Starts on a sheet use the sign of the normal
velocity. Repeated grazing ends with `SPECTRE_FO_STEP_LIMIT`. Because a
sheet is $\rho_g = $ integer, stage locations need no Cartesian-to-logical
inversion, except in $\rho_g < 0.05$ of volume 1, where Cartesian position
and a regular $(\xi,\eta)$ axis chart are used.

## Field for FO

`src/spectre_fo_field.f90` builds per volume a quintic spline of $R$, $Z$,
$A_\vartheta$ and $A_\zeta$ on `spectre_fo_grid`, from which the basis, its
dual and $\mathbf B = \nabla\times\mathbf A$ follow. Volumes share the sheet
geometry samples; $\mathbf A$ keeps each volume's gauge (zero on its inner
boundary). Axisymmetric fields use seven identical toroidal planes. Volume 1
has an additional axis table of $R, Z, A_\xi, A_\eta, A_\zeta$, blended with
the polar representation over $0.05 < \rho_g < 0.1$. This is an
interpolation of the equilibrium: refine `spectre_fo_grid` together with
`spectre_fo_tolerance` in loss-convergence studies.

## FO to GC

The first-order map removes the Larmor vector,
$\mathbf X = \mathbf x - (\bar\rho_0/B)\,\mathbf b\times\bar{\mathbf v}_\perp$,
locates $\mathbf X$, and gives $\bar v_\parallel = \bar{\mathbf v}\cdot\mathbf b(\mathbf X)$
and $\bar v_\perp^2 = \bar v^2 - \bar v_\parallel^2$. The exit test uses it at one
instant. When the test passes, a scratch copy of the FO is integrated for one
gyroperiod $2\pi\bar\rho_0/B(\mathbf X)$ and the map is applied at 32 equally
spaced times. $\mathbf X$ (Cartesian), $\bar v_\parallel$ and $\bar v_\perp^2$
are averaged, and the GC is placed at $\langle\mathbf X\rangle$ with
$J_\perp = \langle\bar v_\perp^2\rangle/(2B(\langle\mathbf X\rangle))$,
$v/v_0 = |\bar{\mathbf v}|/\sqrt2$ and $\bar v_\parallel$ from the energy with
the sign of $\langle\bar v_\parallel\rangle$. The scratch orbit is discarded, so
the GC resumes at the exit time with the state centred half a gyroperiod
later.

In an axisymmetric field the GC $\bar p_\zeta = \bar v_\parallel h_\zeta +
A_\zeta/\bar\rho_0$ is then restored to the entry value, moved to the exit
volume's gauge by the summed sheet jumps of $A_\zeta$. $\mathbf X$ is moved
along $\nabla A_\zeta/|\nabla A_\zeta|^2$ at fixed $J_\perp$ and energy, with
$\bar v_\parallel$ from the energy, until $\bar p_\zeta$ agrees to $10^{-12}$
(at most 12 steps). Fixing $\bar v_\parallel$ from $\bar p_\zeta$ instead would
bring back the $O(\rho_*)$ bias in $J_\perp$.

If the scratch orbit crosses back over the interface of this FO phase, the
particle still straddles the sheet and the FO continues; the next exit test
follows a gyroperiod later. If it reaches another interface or the outer
boundary, or the average fails, the exit falls back to the instantaneous
map at the exit time, with $\bar v_\parallel$ from the entry $\bar p_\zeta$ in
axisymmetry. `spectre_fo_to_gc` is the instantaneous map.

On the axisymmetric tok2vol fixture at $\rho_r = 0.02$
(`test_spectre_switch_mu_bias`), the mean change of $J_\perp$ per passage
that the switches add to that of the full orbit is $-3.2\times10^{-4}$
with the instantaneous maps and below $10^{-10}$ with the averaged ones,
and $\bar p_\zeta$ is conserved to $10^{-12}$. Each averaged switch costs
about one extra gyroperiod of FO, an entry one per lift.

## Losses and limits

With `spectre_outer_loss = 'guiding_centre'` (default) a marker is lost when
its guiding centre reaches $M_\text{vol}$, the definition used for VMEC and
GVEC fields. A full orbit whose particle reaches $M_\text{vol}$ while its
guiding centre is inside hands over to the guiding centre, since the
equilibrium defines no field outside; an entry tries further gyrophases
before leaving such a marker in the guiding-centre model. With `'particle'`
a marker is lost when the full orbit touches $M_\text{vol}$, a
finite-Larmor-radius loss that counts more losses near a compressed edge. A
location failure near the edge is a numerical failure. RK45 is not
symplectic and can drift in energy over long FO phases. Losses are read
from `times_lost.dat` (column 2, $0 < t <$ `trace_time`); FO loss counters
cover only FO phases. The FO settings are applied during shared SPECTRE
field initialization, identically for `simple.x` and `simple_api`/`pysimple`.

Symplectic SPECTRE warning-mode recovery can retain the last valid state for an unresolved time
interval and continue tracing. Such holds do not establish a resolved physical
trajectory. From the first affected macrostep, the marker is excluded from the
resolved confined fraction and included in `unresolved_fraction.dat`. Its public
exit code is 104 (numerical event) and its `times_lost` value is NaN, even if it
later reaches a physical boundary. A sparse `orbit_warning_holds.dat` preserves
one row per affected marker: marker index, cumulative held-interval count, first
affected macrostep (one-based), terminal native exit code and terminal native
time in seconds. The terminal label is evidence of the final event, not
acceptance of its prior trajectory. No-hold completion and physical-loss
reporting are unchanged. A hold count measures unresolved intervals, not Newton
iterations, output samples or full-orbit attempts; repeated identical failures
can bypass the complete recovery cascade. The integration and recovery policy
is unchanged by this reporting distinction.
The legacy SPECTRE RK45 warning-hold path is outside this reporting change.
