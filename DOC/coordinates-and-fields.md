# Coordinate Systems and Magnetic Fields in SIMPLE

This document describes the coordinate systems and magnetic-field paths used by
SIMPLE. It focuses on what exists in the code and how the main runtime paths
fit together.

Update this file when coordinates or field-loading paths change.

The optional IMAS actor imports axisymmetric rectangular equilibrium IDS through
libneo's G-EQDSK-to-Boozer converter. It maps DD3 COCOS11 and DD4 COCOS17 to
COCOS3 explicitly, preserving SI lengths and flux normalization. Cartesian
chartmap geometry covers the full toroidal period even for an axisymmetric field;
the covariant toroidal field is G, not the transform normalization G + iota*I.
Uniform rho and the independent uniform s abscissa for A_phi are distinct grids.

`simple_cylindrical_to_reference` and `simple_reference_to_cylindrical` expose
the loaded libneo reference coordinate map with `(R [m], phi [rad], Z [m])` at
the boundary. VMEC reference radius is s; chartmap reference radius is rho.
Failed inverses return a per-marker status and Python raises an error. Normalized
TEST and SPECTRE references are rejected by this physical-coordinate API.
These transforms do not change the integrator-coordinate conventions.
In particular, the tabulated collision profile abscissa remains s for both VMEC
and Boozer integration, so IDS rho_tor_norm is squared before loading profiles.
See the IMAS workflow in `python/README.rst` for supported IDS and weight semantics.

Fusion birth sampling (`startmode = 7`) uses the VMEC volume Jacobian in
`(s,theta,phi)` for VMEC input, or the field-derived Boozer Jacobian
`torflux*(iota*B_theta+B_zeta)/B^2` in integration `(s,theta_B,zeta_B)` for
Boozer chartmaps. The latter is a physical volume element only when the
chart geometry and field data are Boozer-consistent; this consistency is
not checked at load time. It converts chartmap starts to
reference coordinates only after sampling; profile arguments remain toroidal
flux `s`, not reference `sqrt(s)`. See [alpha births](https://github.com/itpplasma/SIMPLE/blob/main/docs/alpha-birth.md)
for source-grid resolution, units, supported inputs and marker normalization.

---

The equations behind canonical, Boozer, Meiss and Albert coordinates are in
{doc}`/model/coordinates`; the module map and the steps to add a field source
are in {doc}`/development/architecture`. This note documents the code paths:
types, routines, grids, units and conventions at each interface.

## 1. Coordinate systems

### 1.1 Interface (libneo)

Reference coordinates are a libneo `coordinate_system_t`
(`libneo_coordinates`). Every implementation provides

| Routine | Result |
|---|---|
| `evaluate_cart(u, x)` | Cartesian position $(x, y, z)$ |
| `evaluate_cyl(u, xcyl)` | cylindrical $(R, \varphi, Z)$ |
| `covariant_basis(u, e_cov)` | $3\times 3$ covariant basis $\partial\mathbf x/\partial u^i$ |
| `metric_tensor(u, g, ginv, sqrtg)` | $g_{ij}$, $g^{ij}$, $\sqrt g$ |
| `from_cyl(xcyl, u, ierr)` | inverse map from cylindrical coordinates |
| `cov_to_cart`, `ctr_to_cart` | vector components to Cartesian |

The run-wide instance is `reference_coordinates::ref_coords`, set by
`init_reference_coordinates(coord_input)` from the file type: a Boozer
chartmap or chartmap file gives `chartmap_coordinate_system_t` wrapped with
`vmec_RZ_scale` (`scaled_chartmap_coordinates`), a VMEC `wout` gives
`vmec_coordinate_system_t`.

### 1.2 VMEC coordinates

`vmec_coordinate_system_t` uses $(s, \vartheta, \varphi)$: normalised
toroidal flux, VMEC poloidal angle and geometric toroidal angle. $R$ and $Z$
come from `splint_vmec_data`; $x = R\cos\varphi$, $y = R\sin\varphi$. The
basis vectors are

$$
\mathbf e_{s,\vartheta} = \partial_{s,\vartheta} R\,\hat{\mathbf R}
  + \partial_{s,\vartheta} Z\,\hat{\mathbf z},\qquad
\mathbf e_\varphi = R\,\hat{\boldsymbol\varphi} + \partial_\varphi R\,\hat{\mathbf R}
  + \partial_\varphi Z\,\hat{\mathbf z}.
$$

`from_cyl` solves $(R, Z)(s, \vartheta) = (R_*, Z_*)$ at fixed $\varphi$ by
Newton iteration from $s \approx (d_\text{axis}/d_\text{boundary})^2$: at most
50 iterations, residual $10^{-10}$, step $10^{-12}$
(`libneo_coordinates_vmec.f90`). It can fail outside the plasma and near the
boundary ([#288](https://github.com/itpplasma/SIMPLE/issues/288)).

### 1.3 Chartmap coordinates

`chartmap_coordinate_system_t` reads a NetCDF file with Cartesian geometry
$x, y, z(\zeta, \vartheta, \rho)$ on a uniform grid over one field period,
dimensions `rho`, `theta`, `zeta` and attributes `num_field_periods`,
`zeta_convention`, `rho_convention`:

| `rho_convention` | Meaning of `rho` |
|---|---|
| `rho_tor`, `rho_pol` | $\sqrt{s}$ of toroidal or poloidal flux |
| `psi_tor_norm`, `psi_pol_norm` | normalised flux $s$ itself |
| `unknown` | geometric radius, e.g. from `map2disc` |

Two kinds of chartmap describe the same boundary with different interiors.
`vmec_to_chartmap.x` follows the VMEC flux surfaces (`rho_tor`, magnetic
axis at $\rho = 0$). `generate_chartmap_map2disc.py` fills the interior with a
conformal disc map (`unknown`, geometric centre at $\rho = 0$); it serves as a
smooth reference for Meiss coordinates. At equal nominal
$(\rho, \vartheta, \zeta)$ the two give different points and metrics; a
particle started at the same physical position follows the same orbit.

**Extended Boozer chartmap.** A chartmap with `boozer_field = 1` also stores
the field in Boozer coordinates; the format is specified in
{doc}`/boozer-chartmap-schema`. Geometry, `Bmod`, `B_theta` and `B_phi` use
the uniform `rho` grid with $s = \rho^2$; `A_phi` uses a uniform `s` grid
(`radial_abscissa = "s"`); $A_\vartheta = \psi_\text{tor}\, s$ is linear in
$s$. `theta` and `zeta` exclude the periodic endpoint; the reader appends it
for the splines. The major radius is not stored: the reader takes the
$(\vartheta, \zeta)$ average of $\sqrt{x^2 + y^2}$ on the innermost surface
(cm in the file, metres in `new_vmec_stuff_mod::rmajor`). Both consumers use
`boozer_chartmap_io::read_boozer_chartmap`. Writers are
`export_boozer_chartmap` (libneo `boozer_converter.F90`, app
`export_boozer_chartmap.x`), `tools/gvec_to_boozer_chartmap.py` and
`tools/booz_xform_to_boozer_chartmap.py` (from `boozmn*.nc`; recovers
`torflux` from the geometry when `phi_b` is absent).

### 1.4 Cartesian coordinates

`cartesian_coordinate_system_t` (`src/coordinates/cartesian_coordinates.f90`)
is the identity in cm with unit metric. It is the native system of
`coils_field_t`.

### 1.5 Radial scaling

Splines in the integrator radius $r$ resolve the axis better than splines in
$s$. `src/coordinates/coordinate_scaling.f90` provides

| Type | Map | $dr/ds$ |
|---|---|---|
| `sqrt_s_scaling_t` | $r = \sqrt{s}$ | $1/(2r)$ |
| `identity_scaling_t` | $r$ = reference radius | 1 |

`get_default_scaling` (`src/field/field_splined.f90`) selects
`sqrt_s_scaling_t` for VMEC references and for chartmaps with
`psi_tor_norm`/`psi_pol_norm`, and `identity_scaling_t` for chartmaps whose
`rho` is already $\sqrt s$ or geometric. Covariant radial components are
multiplied by $ds/dr = 2r$ under `sqrt_s_scaling_t`
([#291](https://github.com/itpplasma/SIMPLE/issues/291)).

## 2. Field objects

### 2.1 Interface

Fields extend libneo `magnetic_field_t` (module `field_base`) and carry their
native coordinate system in `coords`:

```fortran
subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    real(dp), intent(in) :: x(3)                   ! native coordinates
    real(dp), intent(out) :: Acov(3)               ! covariant A_i
    real(dp), intent(out) :: hcov(3)               ! covariant h_i = B_i/B
    real(dp), intent(out) :: Bmod                  ! B
    real(dp), intent(out), optional :: sqgBctr(3)  ! sqrt(g) B^i
end subroutine
```

`field_from_file` (`src/field.F90`) chooses the class: SPECTRE HDF5,
Boozer chartmap, VMEC `wout`, or a coils file (name starting with `coils` or
ending in `.coils`), which it wraps in a splined field on `ref_coords`.

| Type | File | Native coordinates | `sqgBctr` |
|---|---|---|---|
| `vmec_field_t` | libneo `field_vmec.f90` | VMEC $(s,\vartheta,\varphi)$ | not implemented |
| `coils_field_t` | `src/field/field_coils.f90` | Cartesian | Cartesian $\mathbf B$ |
| `boozer_chartmap_field_t` | `src/field/field_boozer_chartmap.f90` | $(\rho, \vartheta_B, \varphi_B)$ | not implemented |
| `splined_field_t` | `src/field/field_splined.f90` | reference coordinates | not implemented |
| `spectre_field_t` | libneo `field_spectre.f90` | $(\rho_g, \vartheta, \zeta)$ | exact (section 6) |

`sqgBctr` is inconsistent across types
([#290](https://github.com/itpplasma/SIMPLE/issues/290)).

`coils_field_t` evaluates Biot–Savart from libneo `neo_biotsavart`;
`create_coils_field` converts coil positions from m to cm and currents from A
to statA ($\times 2.99792458\times 10^{9}$).

`boozer_chartmap_field_t` needs no VMEC initialisation;
`is_boozer_chartmap` selects this path. On load it scales `Bmod` by
`vmec_B_scale`, `B_theta` and `B_phi` by `vmec_B_scale*vmec_RZ_scale`, and
`A_phi` and `torflux` by `vmec_B_scale*vmec_RZ_scale**2`. GVEC equilibria
enter only through this path after conversion with
`tools/gvec_to_boozer_chartmap.py`; SIMPLE does not link GVEC.

### 2.2 Splined field

`splined_field_t` resamples any field on the reference coordinates and then
forgets the source. `create_splined_field` (`src/field/field_splined.f90`):

1. `configure_spline_grid`: $62 \times 63 \times 64$ points in
   $(r, \vartheta, \varphi)$ over $r \in [10^{-12}, 1]$,
   $\vartheta \in [0, 2\pi]$, $\varphi \in [0, 2\pi/N_\text{fp}]$. For a
   chartmap reference it takes the chartmap grid size and $\rho$ range,
   pulled $0.1\,\Delta\rho$ inside the outermost surface, where the VMEC
   inversion can fail.
2. `get_default_scaling` (section 1.5).
3. `sample_spline_arrays` evaluates $A_i$, $h_i$ and $B$ at every node
   (OpenMP over $\varphi$), see section 4.
4. Gauge: subtract the value at the first node from each $A_i$.
5. `construct_batch_splines_3d`: quintic in all directions, periodic in the
   angles. Seven quantities: $A_r, A_\vartheta, A_\varphi, h_r, h_\vartheta,
   h_\varphi, B$.

At run time `evaluate_with_der` returns the seven values and their first
derivatives; `magfie_refcoords` adds the metric from `ref_coords`.

## 3. The `magfie` interface

`magfie_sub` (`src/magfie.f90`) evaluates the non-canonical field for RK45,
sampling, classification and diagnostics through a procedure pointer:

```fortran
subroutine magfie_base(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: bmod, sqrtg
    real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)
end subroutine
```

| Output | Meaning |
|---|---|
| `bmod` | $B$ (G) |
| `sqrtg` | signed Jacobian $\sqrt g$ |
| `bder` | $\partial_i \ln B$ |
| `hcovar`, `hctrvr` | $h_i$, $h^i$ |
| `hcurl` | $(\nabla\times\mathbf h)^i$ |

`init_magfie(id)` sets the pointer; the ids `TEST = -1` to `SPECTRE = 6` are
the `integ_coords` values ({doc}`/user-guide/field-sources`). The pointer is
global ([#292](https://github.com/itpplasma/SIMPLE/issues/292)).

`magfie_refcoords` (id `REFCOORDS = 5`) works with any field installed by
`set_magfie_refcoords_field`:

1. `evaluate_with_der` gives $A_i$, $h_i$, $B$ and $\partial_j h_k$,
   $\partial_j B$ in scaled coordinates $x$.
2. `scaled_to_ref_coords` maps $x$ to reference $u$ with Jacobian
   $J = du/dx$.
3. `metric_tensor` and `covariant_basis` at $u$ give $g^{ij}$ and the
   orientation; $\sqrt g_x = \operatorname{sign}(\det \mathbf e)\,
   \sqrt{|g_u|}\, J$ (`signed_jacobian`).
4. $h^i = g_x^{ij} h_j$ with the scaled inverse metric
   (`inverse_metric_scaled`).
5. `compute_hcurl`: $(\nabla\times\mathbf h)^i = \epsilon^{ijk}\partial_j
   h_k/\sqrt g$, from `dh(j, k)` $= \partial h_k/\partial x_j$ (derivative
   index first). Fields without analytic derivatives go through
   `magfie_generic_refcoords`, which builds the same array by central
   differences.

## 4. Sampling a field on reference coordinates

`sample_spline_arrays` covers three cases.

- **VMEC source, VMEC reference**: evaluate at $(s, \vartheta, \varphi)$
  directly; apply the radial scaling to the covariant $r$ components.
- **VMEC source, other reference** (e.g. chartmap): $u_\text{ref} \to
  \mathbf x$ (`evaluate_cart`), then `find_vmec_coords_for_cart_point`
  converts $\mathbf x$ to cylindrical coordinates and calls the VMEC
  `from_cyl`; evaluate there, convert $A$ and $h$ to Cartesian vectors
  (`cov_to_cart`) and project on the reference basis,
  $A_i = \mathbf A\cdot\mathbf e_i$.
- **Cartesian source** (coils): evaluate at $\mathbf x$ and project as above.

With a chartmap reference, symplectic and RK45 orbits can differ already in
$(\rho, \vartheta, \zeta)$; a difference seen in $(R, Z)$ is not
necessarily caused by the chartmap inversion.

## 5. Canonical coordinates

Canonical fields fill `field_can_t` (`src/field/field_can_base.f90`) with
$A_\vartheta, A_\varphi, h_\vartheta, h_\varphi, B$ and their first and
second derivatives in $(r, \vartheta, \varphi)$; the radial components of
$\mathbf A$ and $\mathbf h$ vanish by construction. The constructions are
module-level singletons: one canonical field per process
([#293](https://github.com/itpplasma/SIMPLE/issues/293)). They are written
during initialisation and read-only during tracing.

| Coordinates | File | $r$ | Construction |
|---|---|---|---|
| canonical flux | `src/field/field_can_flux.f90` | $s$ | libneo `vmec_to_can`, `can_to_vmec`, `splint_can_coord` |
| Boozer | `src/field/field_can_boozer.f90` | $s$ | libneo `splint_boozer_coord`, `vmec_to_boozer`, `boozer_to_vmec` |
| Meiss | `src/field/field_can_meiss.f90` | reference radius after scaling | radial ODE on a grid |
| Albert | `src/field/field_can_albert.f90` | $\psi = A_\vartheta/A_\text{norm}$ | Meiss, then regrid |
| SPECTRE | `src/field/field_can_spectre.f90` | $\rho_g$ | Meiss per volume (section 6) |
| test | `src/field/field_can_test.f90` | $r$ | analytic circular tokamak |

For Boozer chartmaps the Boozer splines come from the chartmap
(`load_boozer_from_chartmap`) and `integ_to_ref_boozer_chartmap` maps
$s \to \rho = \sqrt s$. The RK45 GPU path evaluates the same Boozer splines
on the device; `test_gpu_orbit_bench` compares it with the host.

**Meiss.** The construction needs only a toroidal covariant field component
that does not vanish, not flux surfaces, so VMEC, chartmap and coil fields on
any reference work. Phase 1 integrates the transformation ODE radially on the
`canonical_grid_*` grid (tolerance `canonical_ode_relerr`) and stores
$\lambda_\varphi$ and the gauge $\chi$ in batch splines; phase 2 samples
$A_\vartheta, A_\varphi, h_\vartheta, h_\varphi, B$ on the canonical grid.
The transforms are

$$
r = \text{scaled reference radius},\qquad \vartheta = \vartheta_c,\qquad
\varphi = \varphi_c + \lambda_\varphi(r, \vartheta_c, \varphi_c),
$$

`integ_to_ref_meiss` evaluates this directly; `ref_to_integ_meiss` solves
for $\varphi_c$ by Newton iteration (tolerance $10^{-12}$, at most 16
iterations). `diag_meiss.x` plots the stored components.

**Albert.** `init_albert` builds Meiss coordinates, maps $r \to \psi =
A_\vartheta(r)/A_\text{norm}$ with $A_\text{norm} = \pm\max|A_\vartheta|$
(`src/field/psi_transform.f90`) and regrids all components to $\psi$. Then
$A_\vartheta = A_\text{norm}\psi$ and $\partial_\psi A_\vartheta = A_\text{norm}$
are exact. `diag_albert.x` plots the four stored components at three radii
with fortplot (disabled without fortplot).

## 6. SPECTRE

SPECTRE/SPEC MRxMHD equilibria (HDF5) are read by libneo `spectre_reader`
(`load_spectre`); `field_input` is recognised by the `.h5` extension or the
HDF5 signature (`field::is_spectre_file`) and loaded by `create_spectre_field`,
which builds the stacked chart $\rho_g \in [0, M_\text{vol}]$, volume
$l = \min(\lfloor\rho_g\rfloor + 1, M_\text{vol})$. libneo `spectre_basis`
returns $(A_\vartheta, A_\zeta)$ and derivatives in the local coordinate
$s_l \in [-1, 1]$; the polynomial basis is valid beyond $|s_l| = 1$. The run
does not initialise VMEC. `init_spectre_field` (`src/simple_main.f90`) sets
$N_\text{fp}$, `fper` $= 2\pi/N_\text{fp}$ and `rmajor`, the $m = n = 0$
harmonic of $R$ on the outermost interface (m), so that `params_init` derives
$\Delta\tau$ as for VMEC.

### 6.1 Field lines

`src/field_line_spectre.f90` traces field lines inside one volume as a
one-degree-of-freedom Hamiltonian flow with time $\zeta$, coordinate
$\vartheta$ and momentum $\psi = A_\vartheta$; $s$ follows implicitly from
$\psi$:

$$
\frac{d\vartheta}{d\zeta} = -\frac{\partial_s A_\zeta}{\partial_s A_\vartheta},
\qquad
\frac{d\psi}{d\zeta} = \partial_\vartheta A_\zeta
 - \partial_\vartheta A_\vartheta\,\frac{\partial_s A_\zeta}{\partial_s A_\vartheta}.
$$

The map is semi-implicit symplectic Euler: Newton on $s$ for the momentum
update (analytic Jacobian, $|F| < 10^{-13}$), then explicit $\vartheta$ and
$\zeta$ updates. `spectre_poincare.x` (`app/spectre_poincare.f90`) reads
`spectre_poincare.in` and writes `poincare_vol<N>.dat` (seed, section, $s$,
$\vartheta$, $R$, $Z$) and `iota_vol<N>.dat` (least-squares rotational
transform) at $\zeta_0 + 2\pi k/N_\text{fp}$. $R, Z$ come from the interface
Fourier blend `spectre_rz`: linear in $s$ between interfaces, $\bar s^{m}$
in the axis volume.

### 6.2 Guiding centre, RK45 (`integmode = 0`)

`magfie_spectre` (`src/magfie.f90`) evaluates on $(\rho_g, \vartheta,
\zeta)$. One call of libneo `spectre_evaluate_der` gives $B$, $h_i$, the exact
$\sqrt g\,B^i$ and the derivatives of $B$ and $h_i$; the chart supplies
$g_{ij}$ and $\partial_k g_{ij}$. `bder` and `hcurl` follow analytically;
libneo returns $\partial h_i/\partial x_j$ with the component index first,
so `magfie_spectre` transposes it for `compute_hcurl`. Derivatives are taken
inside volume $l$, so the jump of $B$ at an interface never enters a drift.

The field is SI; the conversion to Gaussian units is in one block:

| Quantity | Factor | Reason |
|---|---|---|
| `bmod` | $10^4$ | T to G |
| `sqrtg` | $10^6$ | $L^3$, m to cm |
| `hcovar` | $10^2$ | $L$ |
| `hctrvr` | $10^{-2}$ | $1/L$ |
| `bder` | 1 | chart is dimensionless |

`sample_spectre_surface` (`startmode = 1`) samples the surface
$\rho_g$ = `sbeg(1)` with probability $\propto |\mathbf e_\vartheta \times
\mathbf e_\zeta|$ on the $(\vartheta, \zeta)$ nodes, uniform speed and random
pitch; `start.dat` holds $(\rho_g, \vartheta, \zeta)$ directly. It also sets
`bmin`, `bmax`, `bmod00`. The single-volume RK45 path (`src/spectre_orbit.f90`)
clamps evaluation at the boundary, reflects at the axis and absorbs at
$\rho_g = 1$. Multi-volume fields require `integmode > 0`.

### 6.3 Guiding centre, symplectic per volume (`integmode > 0`)

In SPECTRE coordinates the gauge $A_s = 0$ holds, but $h_s \neq 0$ (the
metric is not diagonal), so the guiding-centre one-form keeps a
$\rho_\parallel B_s\,ds$ term. The Meiss angle $\zeta_c$ with
$\partial_s\lambda = -B_s/B_\zeta$ and the gauge $\chi$ with
$\partial_s\chi = (\partial_s\lambda) A_\zeta$ remove both $A'_s$ and $B'_s$,
so the flow is canonical in $(\rho_g, \vartheta_c, \zeta_c)$.

`construct_spectre_coordinates` runs `init_meiss` and
`get_meiss_coordinates` once per volume on the smooth slab
$r \in [l-1, l]$ (volume 1 from $r = 10^{-3}$, avoiding the axis) with
identity radial scaling, on the grid `spectre_ncon_r` $\times$
`spectre_ncon_th` $\times n_\varphi$ (default $48 \times 48$).
`harvest_meiss_volume` moves the splines (`move_alloc`) into a
`meiss_volume_t` slot. `spectre_ncon_phi = -1` (default) uses
$n_\varphi = 32$ for fields with toroidal harmonics and `AXISYM_NPHI` $= 8$
for axisymmetric fields; a positive value is used as given. The quintic
coefficients take $n_q \cdot 6^3 \cdot n_r n_\vartheta n_\varphi$ doubles per
volume, which dominates peak memory; the axisymmetric clamp saves about a
factor 4. `set_spectre_construction_grid` passes the grid because
`field_can_*` cannot use `params` without a dependency cycle.

The construction stays in SI; `evaluate_spectre_can` converts values and
derivatives: $B$ by $10^4$, $h_\vartheta, h_\varphi$ by $10^2$,
$A_\vartheta, A_\varphi$ by $10^8$ (T m$^2$ to G cm$^2$), so that
$\sqrt g = (h_\varphi A_\vartheta' - h_\vartheta A_\varphi')/B$ is in
cm$^3$ as in `magfie_spectre`. It evaluates volume $l$; `integ_to_ref` and
`ref_to_integ` apply that volume's $\lambda_\varphi$ with identity radius.
Field evaluation keeps a $C^1$ radial extension past each volume edge for
Newton trials (bound $\rho_g = M_\text{vol} + 1$); this is not a crossing
model. Loss is located at $\rho_g = M_\text{vol}$.

### 6.4 Full orbit near interfaces

With `crossing_level = 3` (default), multi-volume tracing switches to full
orbit near interior interfaces and the outer boundary; the switch rules are
in {doc}`/development/design/spectre-interface-crossing`. A launch whose
gyro-offset places the particle outside is an immediate loss at zero time.

`src/spectre_fo_field.f90` splines $R$, $Z$, $A_\vartheta$ and $A_\zeta$ per
volume on `spectre_fo_grid` (default $32 \times 64 \times 32$ in
$\rho, \vartheta, \zeta$ per period; seven $\zeta$ planes for axisymmetric
fields). Basis vectors and $\mathbf B$ follow from the spline derivatives
with libneo's positive metric Jacobian. The orbit advances logical position
by $d\mathbf u/dt = \mathbf E^{-1}\mathbf v$ and Cartesian velocity by the
Lorentz force, so RK stages need no coordinate inversion. Near the volume-1
axis a regular table in $(\xi, \eta) = \rho_g(\cos\vartheta, \sin\vartheta)$
of $R$, $Z$ and three potential components, with
$2\max(3, n_\rho/4) + 1$ points per transverse direction, replaces the polar
splines; the two join over $0.05 < \rho_g < 0.1$, including the derivative of
the blending weight. The tables are built before the parallel region, shared
read-only and released on reinitialisation. Results must be converged in
the spline grid as well as in the time step.

## 7. libneo entry points

| Routine | Use |
|---|---|
| `spline_vmec_data`, `splint_vmec_data` | VMEC splines and their evaluation |
| `vmec_to_boozer`, `boozer_to_vmec`, `splint_boozer_coord` | Boozer coordinates |
| `vmec_to_can`, `can_to_vmec`, `splint_can_coord` | canonical flux coordinates |
| `make_vmec_coordinate_system`, `make_chartmap_coordinate_system`, `detect_refcoords_file_type` | reference coordinates |
| `create_spectre_field`, `spectre_evaluate_der` | SPECTRE fields |

`init_vmec` (`src/simple.f90`) calls `spline_vmec_data`, then `stevvo` for
$R_0$ and $N_\text{fp}$ (`fper` $= 2\pi/N_\text{fp}$), then
`volume_and_B00`.

## 8. Orbit-level interfaces

### 8.1 Canonical frequencies

`orbit_frequencies::compute_canonical_frequencies(tracer, initial_state,
options, result)` takes an initialised `tracer_t` (`dtaumin`, `v0`, `relerr`,
`integmode` set) and a state $(s, \vartheta, \varphi, v/v_0, \lambda)$ in
integration coordinates. It copies the tracer, so the caller's orbit state
does not advance. The result holds the bounce or transit frequency
$\omega_b = 2\pi/T > 0$ and the mean toroidal frequency
$\omega_\varphi = \Delta\varphi/T$ (rad/s, period in s). The flat wrapper
`compute_canonical_frequencies_flat` in `simple_main` takes a reference
state and is what f90wrap and `pysimple.compute_canonical_frequencies` call.

### 8.2 Initial states from invariants

`src/orbit_invariants.f90` maps the axisymmetric guiding-centre invariants
$(H_0, J_\perp, p_\varphi)$ to starting states and back, outside the
tracing loop. With $p = v/v_0$ and $\lambda = v_\parallel/v$,

$$
H_0 = p^2,\qquad J_\perp = \frac{p^2(1-\lambda^2)}{B},\qquad
\psi_*^\text{native} = A_\varphi + \rho_0\, p\lambda\, h_\varphi
 = \frac{c}{q} P_\varphi .
$$

No electrostatic potential enters. $\psi_*$ has magnetic-flux units; for VMEC
$A_\vartheta$ is the toroidal and $-A_\varphi$ the poloidal flux. The public
$p_\varphi$ uses a portable gauge,

$$
p_\varphi = \epsilon_\psi\,(\psi_*^\text{native} - A_{\varphi,\text{axis}}),
\qquad \epsilon_\psi = \operatorname{sign}(A_{\varphi,\text{edge}} - A_{\varphi,\text{axis}}),
$$

so its flux part is zero on axis and grows outward. `axis_flux_native` and
`edge_flux_native` are $A_\varphi$ on the axis and at the LCFS (for an
extended chartmap at `rho_lcfs`). POTATO values convert as
$p_\varphi = \operatorname{sign}(\psi_\text{edge} - \psi_\text{axis})
(\psi_*^\text{POTATO} - \psi_\text{axis})$; with $a = v_0^\text{POTATO} /
v_0^\text{SIMPLE}$, $H_0^\text{SIMPLE} = a^2 H_0^\text{POTATO}$ and
$J_\perp^\text{SIMPLE} = a^2 J_\perp^\text{POTATO}$. The symplectic state is
encoded by `encode_symplectic_state` ({doc}`/model/guiding-centre`).

The inverse searches the section $\nabla\psi \times \nabla B = 0$ (the $B$
extrema on each flux surface in axisymmetry) on both branches
$\sigma = \operatorname{sign} v_\parallel$ and returns every regular root,
ordered by major radius from the high- to the low-field side (POTATO's
$R_c$), with branch, sign and minimum/maximum metadata. No root is singled
out: several section points can belong to one orbit, and "inner"/"outer"
are topological classes, not field sides (Buchholz et al., J. Phys.: Conf.
Ser. 2397, 012012 (2022),
[doi:10.1088/1742-6596/2397/1/012012](https://doi.org/10.1088/1742-6596/2397/1/012012)).
Entry points: `invariants_from_state`, `states_from_invariants`,
`potato_to_simple_invariants`; the flat wrappers in `simple_main` convert
coordinates and add $R_c$ for Python.
