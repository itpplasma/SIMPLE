# Full-orbit integration in SIMPLE: feasibility analysis

Working notes on adding a Boris/VPA or symplectic full-orbit (Lorentz) pusher
alongside the existing symplectic guiding-center tracer. Records what the switch
costs, where it breaks across coordinate systems, and how ASCOT5 and VENUS-LEVIS
handle the stellarator-boundary problem. Source references are to the local trees
`/Users/ert/code/SIMPLE` and `/Users/ert/code/ascot5` and to the cited papers.

## 1. What SIMPLE integrates now

A 4D canonical guiding-center system, state `z = (r, theta, phi, p_phi)`, with the
magnetic moment `mu` carried as an adiabatic-invariant parameter, not a state
variable.

- `H = v_par^2/2 + mu*Bmod`
- `p_theta = h_theta*v_par + A_theta/ro0`
- `v_par = (p_phi - A_phi/ro0)/h_phi`

State type: `symplectic_integrator_t%z(4)` in
`src/orbit_symplectic_base.f90`. The integrators read only `Ath, Aph, hth, hph,
Bmod` and their first and second derivatives from `field_can_t`
(`src/field/field_can_base.f90:12-35`). The Hamiltonian and canonical momentum are
assembled in `src/field_can.f90` (`get_val`).

Methods (`src/orbit_symplectic_base.f90:12-14`, selected by `integmode`):
`RK45=0` (quasi, non-symplectic), `EXPL_IMPL_EULER=1`, `IMPL_EXPL_EULER=2`,
`MIDPOINT=3`, `GAUSS1..4 = 4..7`, `LOBATTO3=15`. All except RK45 are
implicit-symplectic, solved by Newton iteration.

The symplecticity comes from the canonical Hamiltonian structure of guiding-center
theory, not from the integrator alone. The canonical coordinate systems (Boozer,
flux, Meiss, Albert) are constructed so that `A_r = 0` and `h_r = 0`; only
`A_theta, A_phi, h_theta, h_phi, |B|` are splined. That radial-gauge construction
is exactly what makes the 4D motion canonical.

## 2. Two routes to full orbit

### Route A: Boris / volume-preserving (VPA)

The method ASCOT5 uses. Pushes the 6D Lorentz system in Cartesian, where the
magnetic substep is an exact rotation in SO(3). Volume-preserving and
energy-bounded, **not symplectic** (Qin 2013, proven; Ellison 2015 shows no
variational form exists). Timestep at gyro-resolution, roughly 1/10 to 1/100 of
the gyroperiod.

Code requirements in SIMPLE:
- 6D state: drop `mu` as a parameter, carry the full velocity vector. The
  `p_theta`/`v_par` reduction is gone.
- The real magnetic field vector `B(x)` at the particle, in a metric-flat frame.
  The canonical layer gives `h_theta, h_phi, |B|` with `h_r = 0`, which is
  insufficient. The splined reference field (`splined_field_t`) does store all
  seven components `A_r, A_theta, A_phi, h_r, h_theta, h_phi, |B|`
  (`src/field/field_splined.f90`), but only to first derivative order.
- Roughly 100 to 1000 times more steps per orbit than guiding-center.

### Route B: canonical symplectic full orbit

Needs `H = (p - qA)^2/2m` with `p = mv + qA`. Non-separable in `(x,p)`, hence
intrinsically implicit, and `p` is gauge-dependent. Explicit variants exist
(Zhang 2018; Xiao & Qin 2019) through generating-function or Hamiltonian-splitting
machinery that shares nothing with SIMPLE's current Newton-on-implicit-RK scheme.
Requires the vector potential `A(x)` in the chosen gauge. This is a new integrator
class, not an extension of `orbit_timestep_sympl`, and it does not reuse the
property that makes the guiding-center scheme symplectic.

## 3. Field and coordinate gaps

From `src/magfie.f90:15-32`, the `magfie` interface returns `bmod, sqrtg, bder,
hcovar, hctrvr, hcurl`. It does not return the vector potential `A`, the separate
covariant `B` components, the metric tensor `g_ij`, or Christoffel symbols.

| Quantity | Guiding-center needs | Full orbit needs | Available in SIMPLE |
| --- | --- | --- | --- |
| `h = B/|B|`, `|B|`, `grad ln|B|` | yes | yes | yes |
| full `B` vector (3 comp) | no | yes | only in reference/splined coords |
| vector potential `A` | partial | yes (Route B) | reference coords; `A_r=0` in canonical |
| metric `g_ij` | no | yes (curvilinear push) | VMEC/chartmap only, not canonical |
| Christoffel `Gamma^k_ij` | no | yes (curvilinear push) | not available; numeric diff |
| second derivatives | yes (some modes) | partial | canonical yes, splined no |

Reading the table by coordinate system:

| Coordinate system | Boris/VPA full orbit | Canonical symplectic |
| --- | --- | --- |
| Cylindrical / Cartesian (reference source) | natural: flat metric, exact rotation, full field available | possible; `A(x)` available, still implicit |
| VMEC `(s,theta,phi)` | needs `g_ij, sqrtg, Gamma` (libneo has them); velocity update no longer a single rotation | non-canonical; not the natural home |
| Boozer / flux / Meiss / Albert | fights the construction: `A_r, h_r` are zero by design, no explicit `g_ij`, no `Gamma` | defeats the purpose; these exist to make GC canonical |

Conclusion: a full-orbit pusher does not cleanly cover all current coordinate
systems. The canonical flux coordinates are the wrong representation for it, both
numerically and geometrically.

## 4. The stellarator-boundary problem

A cylindrical `(R,phi,z)` spline grid does not conform to a 3D plasma boundary, and
VMEC supplies a field only inside the last closed flux surface (LCFS). ASCOT5 and
VENUS-LEVIS resolve this in opposite directions.

### ASCOT5: decouple field grid from boundary

- The field is splined on a plain rectangular `(R,phi,z)` box, tricubic, periodic
  in `phi` over one field period (`B_STS`, `B_3DS`;
  `src/Bfield/B_STS.h:15-24`, `src/spline/interp.h:85-102`). The grid conforms to
  nothing. `psi0/psi1` are normalization constants for `rho`, not a grid boundary.
- The plasma shape lives in a separate 3D triangular wall mesh
  (`src/wall/wall_3d.c`). A particle is lost only when its trajectory intersects a
  wall triangle (ray-triangle, octree-accelerated). Field grid and geometry are
  independent.
- The cost moves to the field input. The box must hold a valid `B` everywhere,
  including the vacuum region between edge and wall, and ASCOT5 does no
  extrapolation (out-of-grid query throws; `src/Bfield/B_STS.c`). VMEC alone cannot
  supply this. The field must come from coils via Biot-Savart, or from VMEC
  extended into the vacuum by virtual casing plus a vacuum solution.

The boundary complexity never touches the field layer: spline a box, ignore the
plasma shape there, push the geometry into a mesh, and pay with a globally-defined
field.

### VENUS-LEVIS: flux coordinates conform to the boundary

- The field stays in VMEC flux coordinates `(s,u,v)`, Fourier in the angles, cubic
  spline in `s` (Pfefferle 2014, CPC 185, 3127). No `(R,phi,z)` box. RMP
  perturbations are computed on a cylindrical mesh but transformed into flux
  coordinates before use.
- The full orbit is integrated directly in curvilinear coordinates, with the
  Lorentz force plus the Christoffel term `-v^m v^n Gamma^i_mn` in the equation of
  motion (geodesic form, thesis Eq. 3.7). The state is the flux coordinates, so `B`
  is a forward Fourier-spline evaluation. There is no per-step inversion.
- The stellarator shape is captured exactly by the VMEC map `R(s,u,v), Z(s,u,v)`;
  flux surfaces fit the plasma by construction. No boundary problem in the field
  layer.
- The cost is the mirror of ASCOT5's. The domain is `s <= 1`; no field outside the
  LCFS, no scrape-off layer, no wall mesh. A particle reaching `s=1` is lost there.
  The pusher carries metric and Christoffel terms, and the axis singularity needs
  explicit handling: half-mesh boundary conditions, a `sqrt(s)` grid substitution
  to regularize the `1/sqrt(s)` inertial terms, and adaptive `dt` near the axis.

### Summary of the trade

ASCOT5: real-space box plus globally-defined field plus wall mesh. Gets vacuum/SOL
physics and true wall loads. Pays with field input and mesh infrastructure.

VENUS-LEVIS: flux coordinates. Boundary is free, no inversion, reuses the
equilibrium representation. Locked to `s <= 1`, pays with Christoffel terms in the
pusher and axis-singularity handling.

## 5. Implications for SIMPLE

The fit depends on what is needed at the edge.

- If the goal stays fast-ion confinement to the LCFS (alpha loss fraction, what
  SIMPLE already computes, lost = reached `s=1`), the VENUS-LEVIS route fits best.
  It avoids the cylindrical grid, the vacuum-field-extension problem, and the
  inversion, and it reuses SIMPLE's flux/canonical field machinery. The honest new
  costs are: full velocity vector, Christoffel terms in the pusher, and
  axis-singularity handling.
- If orbits must reach a physical wall through the vacuum region, only the ASCOT5
  route works. The real work item there is not the spline; it is producing a
  vacuum-extended `B(R,phi,z)` (coils or virtual casing) plus a wall mesh.

Tokamaks are a special case in either route. Axisymmetry collapses the field to
`B(R,Z)`, a 2D spline, and the vacuum field is easy (free-boundary equilibrium or
coils). Cylindrical full orbit is clean and fast for tokamaks. That cleanliness is
a tokamak property, not a general one; for stellarators the cylindrical route's
cost is the global field, not the grid.

Recommended first step if the purpose stays confinement in the closed-flux volume:
prototype a VENUS-LEVIS-style full orbit in flux coordinates, as a second
integrator family selected by a new mode, leaving the canonical guiding-center
integrators untouched. Reserve the cylindrical-box plus wall-mesh route for when
wall loads become the actual goal.

## 6. References

Numerical methods:
- H. Qin et al., "Why is Boris algorithm so good?", Phys. Plasmas 20, 084503
  (2013). DOI 10.1063/1.4818428. Boris preserves phase-space volume, is not
  symplectic, energy error stays bounded.
- N. Ellison, J. Burby, H. Qin, J. Comput. Phys. 301, 489 (2015).
  DOI 10.1016/j.jcp.2015.09.007. Proof that Boris is not variational.
- Y. He, Y. Sun, J. Liu, H. Qin, J. Comput. Phys. 281, 135 (2015).
  DOI 10.1016/j.jcp.2014.10.032. Volume-preserving algorithms; Boris as a special
  case.
- Y. He et al., Phys. Lett. A 381, 568 (2017). K-symplectic (non-canonical
  symplectic) explicit algorithms.
- H. Qin et al., Nucl. Fusion 56, 014001 (2016). DOI 10.1088/0029-5515/56/1/014001.
  Canonical symplectic PIC; implicit, fast local solve.
- R. Zhang et al., Phys. Plasmas (2018). DOI 10.1063/1.5012767. Explicit canonical
  symplectic via generating functions.
- Y. Xiao, H. Qin, Comput. Phys. Commun. (2019). DOI 10.1016/j.cpc.2019.04.003.
  Gauge-independent explicit symplectic.

Codes:
- E. Hirvijoki et al., ASCOT5; full orbit via VPA. Source
  `/Users/ert/code/ascot5/src/simulate/step/step_fo_vpa.c`, field
  `src/Bfield/B_STS.*`, wall `src/wall/wall_3d.*`. Overview arXiv:1908.02482.
- D. Pfefferle, W.A. Cooper, J.P. Graves, C. Misev, Comput. Phys. Commun. 185, 3127
  (2014). DOI 10.1016/j.cpc.2014.08.007. VENUS-LEVIS spline-Fourier field, full
  orbit in flux coordinates.
- D. Pfefferle, PhD thesis EPFL 6561 (2015). DOI 10.5075/epfl-thesis-6561.
  Axis-singularity handling, Eq. 3.7 curvilinear equations of motion.
- D. Muir et al., Comput. Phys. Commun. 267, 108078 (2021). VENUS-LEVIS with SPEC
  discontinuous fields.
