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
| `spectre_fo_entry_larmor` | `1.5` | enter FO closer than this many radial gyroradii to a sheet or the outer boundary |
| `spectre_fo_tolerance` | `1d-9` | relative tolerance of the FO RK45 |
| `spectre_fo_grid` | `32, 64, 32` | FO spline points per volume in $\rho_g$, $\vartheta$, $\zeta$ (per field period) |

## Switching rule

In volume $k$ the GC is advanced by the symplectic integrator in that
volume's Meiss chart (`src/spectre_sympl_orbit.f90`). After each microstep,
`interface_proximity` takes the nearest integer $k' = \mathrm{nint}(\rho_g)$
(an interior sheet or the outer boundary $M_\text{vol}$) and the radial
gyroradius

$$
\rho_r = \frac{\bar\rho_0\,\bar v_\perp}{B}\,|\nabla\rho_g|,
$$

and enters FO if $|\rho_g - k'| < \ell_e\rho_r$ with
$\ell_e = $ `spectre_fo_entry_larmor`. It leaves FO when the reconstructed
GC satisfies $|\rho_g - k'| > (\ell_e + 1)\rho_r$ and the field at the
particle and at the GC agree, $|B_p - B_{gc}| \le 0.05\,B_{gc}$ and
$|\mathbf b_p - \mathbf b_{gc}| \le 0.05$. The exit test runs once per radian
of gyration and at the end of each call.

## GC to FO

`spectre_fo_enter` (`src/spectre_fo_hybrid.f90`) places the particle at

$$
\mathbf x = \mathbf X + \frac{\bar\rho_0}{B}\,\mathbf b\times\bar{\mathbf v}_\perp(\alpha),
\qquad \bar{\mathbf v} = \bar v_\parallel\mathbf b + \bar{\mathbf v}_\perp(\alpha).
$$

In an axisymmetric field the gyrophase $\alpha$ is chosen so that the
particle's toroidal canonical momentum
$\bar p_\zeta = \bar{\mathbf v}\cdot\mathbf e_\zeta + A_\zeta/\bar\rho_0$
equals the GC value: 64 phases are sampled and a sign change is refined by
bisection. If no root exists (strongly passing markers), the best phase is
accepted only if $|\Delta\bar p_\zeta| \le 10^{-3}\max(|\bar p_\zeta|, 1)$;
otherwise the step fails with `SPECTRE_FO_GYROPHASE_FAIL`. The accepted
defect is kept in `spectre_fo_state_t%entry_pzeta_defect`.

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

`spectre_fo_to_gc` removes the Larmor vector,
$\mathbf X = \mathbf x - (\bar\rho_0/B)\,\mathbf b\times\bar{\mathbf v}_\perp$,
locates $\mathbf X$, and sets $v/v_0 = |\bar{\mathbf v}|/\sqrt2$ and
$J_\perp = \bar v_\perp^2/(2B_{gc})$. In an axisymmetric field
$\bar v_\parallel = (\bar p_\zeta - A_\zeta/\bar\rho_0)/h_\zeta$ with the
entry $\bar p_\zeta$ moved to the exit volume's gauge by the summed sheet
jumps of $A_\zeta$, so FO excursions conserve $\bar p_\zeta$ exactly;
otherwise $\bar v_\parallel = \bar{\mathbf v}\cdot\mathbf b_{gc}$.

## Losses and limits

A marker is lost only at a located terminal point beyond $M_\text{vol}$; a
location failure near the edge is a numerical failure. RK45 is not
symplectic and can drift in energy over long FO phases. Losses are read
from `times_lost.dat` (column 2, $0 < t <$ `trace_time`); FO loss counters
cover only FO phases. The FO settings are applied during shared SPECTRE
field initialization, identically for `simple.x` and `simple_api`/`pysimple`.

Warning-mode recovery can retain the last valid state for an unresolved time
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
