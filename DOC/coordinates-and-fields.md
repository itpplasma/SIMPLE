# Coordinate Systems and Magnetic Fields in SIMPLE

This document describes the coordinate systems and magnetic-field paths used by
SIMPLE. It focuses on what exists in the code and how the main runtime paths
fit together.

Update this file when coordinates or field-loading paths change.

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [Coordinate Systems](#2-coordinate-systems)
3. [Magnetic Field System](#3-magnetic-field-system)
4. [The magfie Interface](#4-the-magfie-interface)
5. [Splined Field Architecture](#5-splined-field-architecture)
6. [Canonical Coordinate Systems](#6-canonical-coordinate-systems)
7. [libneo Integration](#7-libneo-integration)
8. [Configuration Options](#8-configuration-options)
9. [Coordinate Transformation Workflows](#9-coordinate-transformation-workflows)
10. [Key Files Reference](#10-key-files-reference)

---

## 1. Architecture Overview

SIMPLE uses a layered architecture for coordinates and fields:

```
┌─────────────────────────────────────────────────────────────────┐
│                     Orbit Integrator Layer                       │
│  (orbit_symplectic.f90, orbit_symplectic_quasi.f90)             │
├─────────────────────────────────────────────────────────────────┤
│                      magfie Interface                            │
│  Unified field evaluation: magfie(x, bmod, sqrtg, bder, ...)    │
├─────────────────────────────────────────────────────────────────┤
│              Canonical Field Layer (optional)                    │
│  field_can_boozer, field_can_meiss, field_can_albert            │
├─────────────────────────────────────────────────────────────────┤
│                 Splined Field Decorator                          │
│  splined_field_t: Fast evaluation via batch B-splines           │
├─────────────────────────────────────────────────────────────────┤
│                 Magnetic Field Implementations                   │
│  vmec_field_t, coils_field_t, gvec_field_t                      │
├─────────────────────────────────────────────────────────────────┤
│                  Coordinate Systems (libneo)                     │
│  vmec_coordinate_system_t, chartmap_coordinate_system_t         │
└─────────────────────────────────────────────────────────────────┘
```

### Design Principles

1. **Polymorphism**: All fields inherit from abstract `magnetic_field_t`
2. **Decorator Pattern**: `splined_field_t` wraps any field for fast evaluation
3. **Strategy Pattern**: Runtime field selection via procedure pointer (`magfie`)
4. **Coordinate Independence**: Physics is independent of coordinate choice

---

## 2. Coordinate Systems

### 2.1 Abstract Base (libneo)

All coordinate systems derive from `coordinate_system_t` in libneo:

```fortran
type, abstract :: coordinate_system_t
contains
    procedure(evaluate_cart_if), deferred :: evaluate_cart
    procedure(evaluate_cyl_if), deferred :: evaluate_cyl
    procedure(covariant_basis_if), deferred :: covariant_basis
    procedure(metric_tensor_if), deferred :: metric_tensor
    procedure(from_cyl_if), deferred :: from_cyl
    procedure :: cov_to_cart => coordinate_system_cov_to_cart
    procedure :: ctr_to_cart => coordinate_system_ctr_to_cart
end type
```

**Core Operations**:
- `evaluate_cart(u, x)`: Convert coordinate u to Cartesian (x, y, z)
- `evaluate_cyl(u, xcyl)`: Convert to cylindrical (R, phi, Z)
- `covariant_basis(u, e_cov)`: Get 3x3 covariant basis vectors
- `metric_tensor(u, g, ginv, sqrtg)`: Compute metric tensor
- `from_cyl(xcyl, u, ierr)`: Invert cylindrical to native coordinates

### 2.2 VMEC Coordinates

**Type**: `vmec_coordinate_system_t`

**Coordinates**: (s, theta, varphi)
- `s`: Normalized toroidal flux, 0 <= s <= 1
- `theta`: VMEC poloidal angle, 0 to 2pi
- `varphi`: Geometrical toroidal angle, 0 to 2pi

**Position Evaluation**:
```
R, Z from splint_vmec_data(s, theta, varphi, ...)
x = R * cos(varphi)
y = R * sin(varphi)
z = Z
```

**Covariant Basis**:
```
e_1 = dR/ds * [cos(phi), sin(phi), 0] + dZ/ds * [0, 0, 1]
e_2 = dR/dtheta * [cos(phi), sin(phi), 0] + dZ/dtheta * [0, 0, 1]
e_3 = R * [-sin(phi), cos(phi), 0] + dR/dphi * [cos(phi), sin(phi), 0] + dZ/dphi * [0, 0, 1]
```

**Inversion (from_cyl)**: Newton iteration on (R, Z) -> (s, theta)
- Initial guess: s ~ (distance_from_axis / distance_to_boundary)^2
- Tolerance: 1e-10, max 50 iterations

### 2.3 Chartmap Coordinates

**Type**: `chartmap_coordinate_system_t`

**Coordinates**: `(rho, theta, zeta)`
- `rho`: normalized radial coordinate
- `theta`: poloidal angle
- `zeta`: toroidal angle over one field period

**Source**: NetCDF file with Cartesian geometry on a structured grid

**Required geometry layout**:
```
Dimensions: rho, theta, zeta
Variables: x(zeta, theta, rho), y(...), z(...)
Attributes: num_field_periods, zeta_convention, rho_convention
```

**Rho conventions**:
- `RHO_TOR`: rho = sqrt(s_toroidal)
- `RHO_POL`: rho = sqrt(s_poloidal)
- `PSI_TOR_NORM`: rho = s (normalized toroidal flux)
- `PSI_POL_NORM`: rho = s (normalized poloidal flux)
- `UNKNOWN`: Geometric disc mapping (e.g., map2disc)

**Common chartmap sources**:

1. `vmec_to_chartmap.x`
   - `rho_convention = rho_tor`
   - interior follows VMEC flux surfaces

2. `generate_chartmap_map2disc.py`
   - `rho_convention = unknown`
   - interior is geometric, not flux-based
   - mainly useful as a smooth reference map for Meiss coordinates

#### Extended Boozer chartmap (field profiles)

A chartmap with `boozer_field = 1` stores Boozer field data next to the
geometry. The radial contract is split by quantity:

- `rho` is uniform and increasing. Geometry, `Bmod`, `B_theta`, and `B_phi`
  use this grid.
- `s` is uniform and increasing. It spans `rho(1)^2` to `rho(n_rho)^2`.
  `A_phi` uses this grid and must carry `radial_abscissa = "s"`.

`rho_tor` means `s = rho^2`. `A_theta = torflux*s` stays linear in `s`.
`theta` and `zeta` are endpoint-excluded and shared by geometry and `Bmod`.
The reader appends exact periodic endpoint planes for spline construction. The
major radius is not stored: the reader derives it as the `(theta,
zeta)`-average of `sqrt(x^2 + y^2)` on the innermost `rho` surface (cm in the
file, metres in `new_vmec_stuff_mod::rmajor`). Both chartmap consumers use
`boozer_chartmap_io.read_boozer_chartmap`, so metadata, scaling, and grid layout
stay shared.

The required NetCDF variables and attributes are listed in
`docs/boozer-chartmap-schema.rst`.

Writers of the extended format: `export_boozer_chartmap`
(`boozer_converter.F90`), `tools/gvec_to_boozer_chartmap.py`, and
`tools/booz_xform_to_boozer_chartmap.py` (from booz_xform `boozmn*.nc`
harmonics; recovers `torflux` from the surface geometry when `phi_b` is
absent).

### 2.4 Cartesian Coordinates

**Type**: `cartesian_coordinate_system_t`

**Coordinates**: (x, y, z) in cm

- Identity coordinate system
- Used by `coils_field_t` (Biot-Savart evaluation)
- Trivial metric (identity), trivial basis (identity)

### 2.5 Coordinate Scaling

**Purpose**: Transform between flux coordinate (s) and radial coordinate (r)
for better spline resolution near magnetic axis.

**Types** (in `coordinate_scaling.f90`):

```fortran
type :: sqrt_s_scaling_t    ! Default
    ! Transform: s -> r = sqrt(s)
    ! Jacobian: dr/ds = 1/(2r)
end type

type :: identity_scaling_t
    ! No transformation: r = s
    ! Used when coordinates already use rho = sqrt(s)
end type
```

**Selection Logic**:
- PSI_TOR_NORM, PSI_POL_NORM: Use sqrt_s_scaling
- RHO_TOR, RHO_POL, UNKNOWN: Use identity_scaling

---

## 3. Magnetic Field System

### 3.1 Abstract Base Class

```fortran
type, abstract :: magnetic_field_t
    class(coordinate_system_t), allocatable :: coords
contains
    procedure(evaluate_interface), deferred :: evaluate
end type
```

**Evaluate Interface**:
```fortran
subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    real(dp), intent(in) :: x(3)           ! Position in field coords
    real(dp), intent(out) :: Acov(3)       ! Covariant vector potential
    real(dp), intent(out) :: hcov(3)       ! Covariant h = B/|B|
    real(dp), intent(out) :: Bmod          ! Field magnitude
    real(dp), intent(out), optional :: sqgBctr(3)  ! sqrt(g) * B^i
end subroutine
```

### 3.2 VMEC Field

**Type**: `vmec_field_t`
**File**: `libneo/src/field/field_vmec.f90`

- Native coordinates: VMEC (s, theta, varphi)
- Evaluation: Calls `splint_vmec_data()` from libneo
- sqgBctr: Not implemented

**Creation**:
```fortran
call create_vmec_field(field)
! Allocates vmec_coordinate_system_t into field%coords
```

### 3.3 Coils Field

**Type**: `coils_field_t`
**File**: `src/field/field_coils.f90`

- Native coordinates: Cartesian (x, y, z)
- Evaluation: Biot-Savart via `neo_biotsavart` library
- Unit conversions: meters -> cm, Amperes -> STATA

**Creation**:
```fortran
call create_coils_field('coils.simple', field)
```

### 3.4 Boozer Chartmap Field

**Type**: `boozer_chartmap_field_t`
**File**: `src/field/field_boozer_chartmap.f90`

- Reads an extended Boozer chartmap NetCDF (`boozer_field = 1`, see section 2)
- Native coordinates: chartmap (rho = sqrt(s), theta_B, phi_B)
- No VMEC initialization; runtime path selected by `is_boozer_chartmap()`

GVEC equilibria enter through this path. SIMPLE reads no raw GVEC state files
and links no GVEC library; convert with `tools/gvec_to_boozer_chartmap.py`
(section 8.2.1).

### 3.5 Splined Field (Decorator)

**Type**: `splined_field_t`
**File**: `src/field/field_splined.f90`

Wraps any magnetic field for fast evaluation via batch B-splines.

**Architecture**:
```
┌─────────────────────────────────────────────────────────────────┐
│  SPLINE CONSTRUCTION (one-time)                                  │
│                                                                  │
│  For each grid point (r, theta, phi) in reference coordinates:  │
│    1. ref_coords%evaluate_cart([r,theta,phi]) -> (x,y,z)        │
│    2. find_vmec_coords_for_cart_point -> (s, theta_v, phi_v)    │
│    3. source%evaluate(s, theta_v, phi_v) -> B_vmec              │
│    4. Transform to reference coords via basis vectors           │
│    5. Store in spline: spline(r, theta, phi) = B_ref            │
│                                                                  │
│  Result: Spline knows B(r,theta,phi) but NOT source origin      │
└─────────────────────────────────────────────────────────────────┘
                            |
                            v
┌─────────────────────────────────────────────────────────────────┐
│  RUNTIME EVALUATION (many calls)                                 │
│                                                                  │
│  Given position (r, theta, phi):                                 │
│    1. evaluate_batch_splines_3d -> (A_cov, h_cov, |B|)          │
│    2. ref_coords%metric_tensor -> g_ij, sqrt(g)                 │
│    3. Compute h^i = g^{ij} h_j                                  │
│                                                                  │
│  Source field is COMPLETELY FORGOTTEN after construction        │
└─────────────────────────────────────────────────────────────────┘
```

**Grid Configuration**:
- Default: 62 x 63 x 64 points (r x theta x phi)
- Spline order: 5 (quintic) in all dimensions
- Periodicity: [False, True, True] (r, theta, phi)

**Chartmap Special Case**:
- Range limited to rho in [0.1, 0.99]
- Avoids axis singularity (0.1) and boundary inversion failure (0.99)
- Covers 99% of plasma volume

**Stored Components** (7 total):
1. A_r (covariant vector potential)
2. A_theta
3. A_phi
4. h_r (covariant unit field direction)
5. h_theta
6. h_phi
7. |B| (field magnitude)

---

## 4. The magfie Interface

**File**: `src/magfie.f90`

Unified field evaluation interface using procedure pointer for runtime dispatch.

### 4.1 Interface Definition

```fortran
abstract interface
    subroutine magfie_base(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg
        real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)
    end subroutine
end interface

procedure(magfie_base), pointer :: magfie => null()
```

### 4.2 Output Quantities

| Variable | Type | Physical Meaning |
|----------|------|------------------|
| bmod | scalar | Magnetic field magnitude B/B_ref |
| sqrtg | scalar | Jacobian sqrt(det(g_ij)) |
| bder(3) | covariant | Derivatives d(ln B)/dx_i |
| hcovar(3) | covariant | Unit field direction h_i = B_i/|B| |
| hctrvr(3) | contravariant | Unit field direction h^i = B^i/|B| |
| hcurl(3) | contravariant | Curl of unit vector (curl h)^i |

### 4.3 Available Modes

| ID | Name | Coordinates | Description |
|----|------|-------------|-------------|
| -1 | TEST | (r, theta, phi) | Analytic circular tokamak |
| 0 | CANFLUX | (s, theta_c, phi_c) | Canonical flux coordinates |
| 1 | VMEC | (s, theta, phi) | Direct VMEC (non-canonical) |
| 2 | BOOZER | (s, theta_B, phi_B) | Boozer coordinates |
| 3 | MEISS | (r, theta_c, phi_c) | Meiss canonical |
| 4 | ALBERT | (psi, theta_c, phi_c) | Albert canonical |
| 5 | REFCOORDS | (varies) | Reference splined field |

### 4.4 Mode Initialization

```fortran
! Select mode at startup
call init_magfie(BOOZER)

! For REFCOORDS, first install the splined field
call set_magfie_refcoords_field(my_splined_field)
call init_magfie(REFCOORDS)
```

### 4.5 REFCOORDS Implementation Details

The REFCOORDS mode (`magfie_refcoords`) is the most flexible:

```fortran
subroutine magfie_refcoords(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    ! 1. Evaluate splined field with derivatives
    call refcoords_field%evaluate_with_der(x, Acov, hcovar, Bmod,
                                           dAcov, dhcov, dBmod)

    ! 2. Handle coordinate scaling (sqrt_s vs identity)
    call scaled_to_ref_coords(refcoords_field%coords, x, u_ref, J)

    ! 3. Get metric from reference coordinates
    call refcoords_field%coords%metric_tensor(u_ref, g, ginv_u, sqrtg_abs)
    call refcoords_field%coords%covariant_basis(u_ref, e_cov)

    ! 4. Compute signed Jacobian
    sqrtg = signed_jacobian(e_cov) * J

    ! 5. Transform to contravariant components
    call inverse_metric_scaled(J, ginv_u, ginv_x)
    hctrvr = matmul(ginv_x, hcovar)

    ! 6. Compute curl
    call compute_hcurl(sqrtg, dhcov, hcurl)
end subroutine
```

---

## 5. Splined Field Architecture

### 5.1 Construction Workflow

```fortran
subroutine create_splined_field(source, ref_coords, field, scaling, ...)
    ! 1. Configure grid
    call configure_spline_grid(ref_coords, ..., n_r, n_th, n_phi, xmin, xmax)

    ! 2. Select coordinate scaling
    scaling_ptr => get_default_scaling(ref_coords, ...)

    ! 3. Sample source field on grid
    call sample_spline_arrays(source, ref_coords, scaling_ptr, ...,
                              Ar, Ath, Aphi, hr, hth, hphi, Bmod_arr)

    ! 4. Build batch splines
    call construct_batch_splines_3d(xmin, xmax, y_batch, order, periodic, spl)
end subroutine
```

### 5.2 Coordinate Transformation During Sampling

For each grid point (r, theta, phi) in scaled reference coordinates:

**Case 1: VMEC source + VMEC reference**
- Direct evaluation, no transformation needed
- Scaling Jacobian applied to covariant components

**Case 2: VMEC source + Different reference (e.g., chartmap)**
```
1. ref_coords%evaluate_cart(x_ref) -> x_cart (Cartesian position)
2. find_vmec_coords_for_cart_point -> x_vmec (VMEC coords for same point)
3. source%evaluate(x_vmec) -> A_cov, h_cov in VMEC basis
4. source%coords%cov_to_cart -> A_cart, h_cart (Cartesian vectors)
5. ref_coords%covariant_basis -> e_cov (reference basis)
6. A_ref = A_cart * e_cov (project to reference coords)
```
When using a chartmap reference, differences between symplectic and RK45
trajectories can already appear in chartmap $(\rho,\theta,\zeta)$ space, so a
discrepancy seen in $R,Z$ is not automatically attributable to the chartmap
inversion alone.

**Case 3: Cartesian source (coils)**
```
1. ref_coords%evaluate_cart(x_ref) -> x_cart
2. source%evaluate(x_cart) -> A_cart, h_cart directly
3. ref_coords%covariant_basis -> e_cov
4. A_ref = A_cart * e_cov
```

### 5.3 The find_vmec_coords_for_cart_point Function

Critical for coordinate system interoperability:

```fortran
subroutine find_vmec_coords_for_cart_point(vmec_coords, ref_coords,
                                           u_ref, x_target, u_vmec)
    ! Convert Cartesian to cylindrical
    xcyl(1) = sqrt(x_target(1)**2 + x_target(2)**2)  ! R
    xcyl(2) = modulo(atan2(x_target(2), x_target(1)), twopi)  ! phi
    xcyl(3) = x_target(3)  ! Z

    ! Use VMEC native inversion (Newton iteration)
    call vcs%from_cyl(xcyl, u_vmec, ierr)
end subroutine
```

### 5.4 Gauge Fixing

Vector potential is defined up to gradient of scalar. Splines remove this freedom:

```fortran
Ar = Ar - Ar(1, 1, 1)      ! A_r(r_min, 0, 0) = 0
Ath = Ath - Ath(1, 1, 1)   ! A_theta(r_min, 0, 0) = 0
Aphi = Aphi - Aphi(1, 1, 1) ! A_phi(r_min, 0, 0) = 0
```

---

## 6. Canonical Coordinate Systems

Canonical coordinates simplify symplectic integration by giving the magnetic
field a special structure.

### 6.1 Common Structure

All canonical coordinate systems transform to a form where:
- Vector potential: A = (0, A_theta, A_phi) or A = (0, psi, A_phi)
- Unit field direction: h = (0, h_theta, h_phi)
- First (radial) component vanishes

### 6.2 Boozer Coordinates

**File**: `src/field/field_can_boozer.f90`
**Coordinates**: (r, theta_B, phi_B)

- Straight field line coordinates
- A_theta, A_phi independent of angles
- Uses `splint_boozer_coord()` from libneo

**Transformation**:
```fortran
call vmec_to_boozer(r, theta_v, phi_v, theta_B, phi_B)
call boozer_to_vmec(r, theta_B, phi_B, theta_v, phi_v)
```

### 6.3 Meiss Coordinates

**File**: `src/field/field_can_meiss.f90`
**Coordinates**: (r, theta_c, phi_c) where r is a generic radial coordinate
in the chosen reference system (often r = sqrt(s) for VMEC-based references).

**Applicability**:
- Meiss canonicalization is coordinate-agnostic: it does not require flux
  coordinates or nested flux surfaces.
- The only structural requirement is a non-vanishing toroidal field component
  in the domain, so VMEC, chartmap, or cylindrical reference coordinates are
  all valid inputs.

**Construction** (two-phase):

**Phase 1: Transformation Computation**
- Integrate ODEs radially to find lambda(r, theta, phi)
- lambda = phi_canonical - phi (angle transformation)
- Store in batch splines

**Phase 2: Field Sampling**
- Sample field components on canonical grid
- Store: A_theta(r), A_phi(r), h_theta, h_phi, B_mod

**Coordinate Transforms**:
```fortran
subroutine integ_to_ref_meiss(xinteg, xref)
    ! (r, theta_c, phi_c) -> (s, theta, phi)
    call evaluate_batch_splines_3d(spl_transform_batch, xinteg, y_batch)
    lambda = y_batch(1)
    xref(1) = xinteg(1)**2   ! s = r^2
    xref(2) = xinteg(2)      ! theta = theta_c
    xref(3) = xinteg(3) + lambda  ! phi = phi_c + lambda
end subroutine
```

### 6.4 Albert Coordinates

**File**: `src/field/field_can_albert.f90`
**Coordinates**: (psi, theta_c, phi_c) where psi ~ A_theta

**Key Property**: dA_theta/dpsi = constant (simplifies symplectic integrator)

**Construction**:
1. Build Meiss coordinates first
2. Transform radial coordinate: r -> psi = A_theta(r)
3. Regrid all field components from r to psi

**Evaluation**:
```fortran
f%Ath = Ath_norm * psi    ! Linear in psi
f%dAth = [Ath_norm, 0, 0] ! Constant derivative
```

### 6.5 Canonical Flux Coordinates

**File**: `src/field/field_can_flux.f90`
**Coordinates**: (s, theta_c, phi_c)

- Uses libneo functions: `vmec_to_can`, `can_to_vmec`
- Simpler than Meiss/Albert but less optimized

### 6.6 SPECTRE Field-Line Coordinates

**Files**: `src/field_line_spectre.f90`, `app/spectre_poincare.f90`
**Input**: SPECTRE/SPEC MRxMHD equilibrium (HDF5), read through libneo
`spectre_reader` (`spectre_data_t`, `load_spectre`); per-volume covariant
vector potential `(A_theta, A_zeta)` and its s/theta/zeta derivatives from
libneo `spectre_basis` (`eval_spectre_vector_potential`). Derivatives are with
respect to the local volume coordinate `s in [-1, 1]`; the polynomial basis is
valid past `|s| = 1`.

**Canonical structure**: inside a single volume the magnetic field line is the
flow of a one-degree-of-freedom canonical system with the toroidal angle `zeta`
as time, the poloidal angle `theta` as coordinate, and `psi = A_theta` as the
conjugate momentum. The radial coordinate `s` is recovered implicitly from `psi`
at fixed `(theta, zeta)`. The derived flow is

```
dtheta/dzeta = -dA_zeta/ds / dA_theta/ds
dpsi/dzeta   =  dA_zeta/dtheta - dA_theta/dtheta * dA_zeta/ds / dA_theta/ds
```

with the explicit `d/dzeta` terms cancelling. `field_line_spectre.f90` advances
`(theta, psi)` with a semi-implicit symplectic Euler map: a Newton solve on `s`
for the implicit momentum update (analytic Jacobian from the `ss` and `st`
second derivatives, `|F|` tolerance `1e-13`), then explicit `theta` and `zeta`
updates. Both update signs set the field-line orientation (the sign of iota).

**App `spectre_poincare.x`**: traces field lines per volume from a namelist
(`spectre_poincare.in`) and writes `poincare_vol<N>.dat` (seed, section, s,
theta, R, Z) and `iota_vol<N>.dat` (least-squares rotational transform). Field
periodicity places every section exactly on `zeta0 + k*2*pi/Nfp`, so no
interpolation is needed. R, Z come from the SPEC interface Fourier blend
(`spectre_rz`, mirroring libneo `libneo_coordinates_spectre` for output only:
linear-in-s between interfaces for interior volumes, the `sbar**m` power law for
the axis volume). Field lines never cross interfaces (`B.n = 0` on both sides),
so no crossing logic is involved. This app is a standalone diagnostic separate
from the guiding-center stack below.

#### Guiding-center RK45 per volume (`integ_coords = 6`)

**Files**: `src/magfie.f90` (`magfie_spectre`), `src/spectre_orbit.f90`,
`src/interface_crossing.f90`, `src/simple_main.f90` (`init_spectre_field`,
`trace_orbit_spectre`), `src/samplers.f90` (`sample_spectre_surface`).

`simple.x` traces guiding centers in a SPECTRE equilibrium on the non-canonical
RK45 path, one volume at a time. On reaching an interface the crossing map
(`crossing_level`, Level 1 by default) switches volume, reflects, or loses the
marker (see interface crossing below); the outermost interface is the loss
surface.

**Detection and loading**. `field_input` is detected as SPECTRE by a `.h5`
extension or the HDF5 magic bytes (`field.is_spectre_file`), and loaded through
libneo `create_spectre_field` (which builds the stacked-rho
`spectre_coordinate_system_t`). The run is VMEC-free: no `init_vmec`, no `wout`.
`init_spectre_field` sets the equilibrium globals the same way the Boozer
chartmap path does, so `stevvo` and `params_init` produce a consistent
`dphi`/`dtaumin`/`fper`: `nper = Nfp`, and `rmajor` is the `m=0, n=0` R harmonic
of the outermost interface (SI meters; `stevvo` scales it to cm). `fper` is
`2*pi/Nfp`. `integmode > 0` builds the per-volume canonical coordinates below.

**`magfie_spectre`** evaluates on the chart `x = (rho_g, theta, zeta)`. The
field supplies `|B|`, covariant `h = B/|B|`, and the exact `sqrt(g) B^i` (the
metric-free 2-form, no finite differences); the coordinate system supplies the
analytic metric. `grad log|B|` (`bder`) and `curl(h)` (`hcurl`), the small drift
corrections, come from central differences of `|B|` and `h`. `|B|` jumps across
SPEC interfaces (the pressure jump), so the radial stencil is kept inside the
current volume `[floor(rho_g), floor(rho_g)+1]`; straddling an integer would
inject the discontinuity into the drift and stall the ODE solver.

**Units**. The SPECTRE field returns SI (Tesla, meter); SIMPLE integrates in
Gaussian-CGS. The conversion lives in one block in `magfie_spectre`:

| quantity | factor | reason |
|----------|--------|--------|
| `bmod` = `\|B\|` | `1e4` | Tesla -> Gauss |
| `sqrtg` (Jacobian) | `1e6` | `L^3`, meter -> cm |
| `hcovar` (covariant `h_i`) | `1e2` | `L` |
| `hctrvr` (contravariant `h^i`) | `1e-2` | `1/L` |
| `bder` = `d log\|B\|/dx` | `1` | dimensionless (the chart is dimensionless) |

The covariant vector potential (`T*m -> G*cm`, factor `1e6`) does not enter the
non-canonical RK45 path, which uses only `bmod`/`sqrtg`/`bder`/`hcovar`/
`hctrvr`/`hcurl`. The alpha normalization then follows the standard SIMPLE chain
unchanged: `v0`, `ro0 = v0 m c / e` (cm*Gauss), and `dtaumin = 2*pi*rmajor/
npoiper2` (cm).

**Sampler**. `startmode = 1` samples the control surface `rho_g = sbeg(1)`:
`(theta, zeta)` grid nodes drawn with probability proportional to the surface
Jacobian `|e_theta x e_zeta|`, uniform velocity magnitude, random pitch. The
chart is the integration coordinate, so `start.dat` holds `(rho_g, theta, zeta)`
directly and every marker sits on the requested interface to round-off. The
sampler also sets `bmin`/`bmax`/`bmod00` from a scan of the surface.

**Interface crossing (Level 0, #443)**. A marker owns a home volume, fixed after
the first step as the pair of integer interfaces bracketing `rho_g`. When a
microstep leaves the home volume, `orbit_timestep_spectre` locates the interface
`rho_g = k` (Illinois false position on the substep to `< 1e-10`) and hands the
landing state to `apply_crossing(y_iface, iface, direction, mvol, level, y_out,
info)` in `interface_crossing.f90`. The `crossing_level` namelist knob selects
the map: `1` (default) is the Level-1 refraction map (#440), `0` is the Level-0
energy rescale (#443) kept for regression comparison.

The Level-0 map holds `(theta, zeta)` and the perpendicular invariant, evaluates
`|B|` from *both* volumes at the same interface point (`rho_g = k -+ 1e-12`, the
per-volume polynomial basis), and rescales the parallel velocity. Writing the
stored invariant `perp_inv = z(4)^2 (1 - z(5)^2)/|B| = 2 mu` (`v_perp^2/|B|`), so
`mu = perp_inv/2` is the magnetic moment:

    v_par'^2 = v_par^2 - 2 mu (B_target - B_home),   v_par = z(4) z(5).

- **Crossing** (`v_par'^2 >= 0`): keep the sign of `v_par`, switch to the
  neighbour volume. `H = v_par^2/2 + mu |B|` (each side in its own volume) is
  conserved to round-off, and `z(4)` (hence energy) is unchanged.
- **Reflection** (`v_par'^2 < 0`): the forbidden crossing into higher `|B|` is a
  magnetic mirror. The marker stays home with `v_par -> -v_par` (the same-side
  energy root), so `|v_par|`, `mu`, and `|B|` are unchanged and `H` is exact. A
  deeply trapped marker whose banana tip sits on the sheet stays pinned there,
  the crude Level-0 behaviour absent a tangential sheet-drift kick.
- **Loss**: crossing outward through the outermost interface (`direction = +1`,
  `iface = Mvol`) removes the marker; its loss time is recorded in
  `times_lost.dat` and `confined_fraction.dat`.
- **Axis**: reaching `rho_g = 0` (the coordinate axis, `sqrt(g) = 0`) reflects
  trivially before the singularity.

After a crossing `rho_g` is nudged `1e-6` into the resolved volume so the next
substep does not re-trigger the event within the odeint tolerance.

**Interface crossing (Level 1, #440)**. `crossing_level = 1` (default) replaces
the rescale with the thin-current-sheet limit of the drift: an impulse `Delta z =
lambda_k X` along the interface Hamiltonian vector field `X = {z, rho_g}`, with
the single scalar `lambda_k` fixed by exact energy conservation. The generator
uses the same `parmot_mod` constant `ro0` and the same
`Bstar_par = |B| + ro0 v_par (h.curl h)` as `velo_can`:

    X_theta = - ro0 hcovar_zeta  / (sqrtg Bstar_par)
    X_zeta  = + ro0 hcovar_theta / (sqrtg Bstar_par)
    X_vpar  = + ro0 v_par (curl h)^rho / Bstar_par.

The energy criterion at the landing point (`radicand0 = v_par^2 - 2 mu [[B]]`, the
Level-0 condition) decides crossing vs mirror, so the event split matches Level 0.
A forbidden crossing reflects as the Level-0 mirror (`v_par -> -v_par`, no kick).
For a crossing the angles `(theta, zeta)` and `v_par` receive `lambda_k X`;
`rho_g` stays on the interface and then switches volume with the Level-0 nudge.
`lambda_k` solves `H_target(z + lambda_k X) = H_home(z)` by a damped (backtracking)
Newton, with the generator refreshed each residual at the midpoint state
`z + (lambda_k/2) X` from both-side quantities averaged, so the discrete map is
symplectic and keeps `mu` and `p` fixed. The crossing lands the angles on an
equal-`|B|` target point instead of rescaling `v_par`. That point sits a
tangential distance `~ mu [[B]] / v_drift_normal` away; where the poloidal `|B|`
gradient along the sheet drift is weak this exceeds the `0.3` rad thin-layer
bound and the map falls back to the energy-exact Level-0 rescale (no kick), so
every energetically allowed marker still crosses. The full derivation is in
`DOC/spectre-interface-crossing.md`.

**Drift regularization**. The guiding-center drift diverges where
`Bstar_par = 1 + p_par (c/e) h.curl h -> 0`, a pitch-dependent phase-space
surface that the interface sheet current pushes into a thin layer beside each
interface. `velo_can_clamped` in `spectre_orbit.f90` freezes the field at the
volume edge minus a `2e-2` band and caps the RHS magnitude so the RK45 substep
integrates a finite, physical drift instead of collapsing to zero step. This
regularizes only the integrated orbit near the layer; the crossing map itself
uses the true both-side `|B|` at the exact interface, so energy conservation is
untouched.

**Crossing log**. Every crossing and reflection appends one line to
`spectre_crossing_events.dat`, grouped by particle in time order (a stable
counting sort, so the file is reproducible under OpenMP): particle, time,
interface, type (1 crossing / 2 reflection / 3 loss / 4 stop), exit volume,
entry volume, theta, zeta, `v_par` before, `v_par` after, mu, `|B|` home,
`|B|` target. Both the RK45 and the symplectic path feed this log; type 4 is
the symplectic path's pathological non-convergence fallback (#441).

#### Guiding-center symplectic per volume (`integmode > 0`)

**Files**: `src/field/field_can_spectre.f90`, `src/field/field_can_meiss.f90`
(`spectre_field_t` arms in `ah_cov_on_slice` / `init_canonical_field_components`,
the `meiss_volume_t` slot, and the `harvest_meiss_volume` mover),
`src/field_can.f90` (id `SPECTRE` dispatch), `src/spectre_sympl_orbit.f90`
(multi-volume microstepping with exact interface landing, #441),
`src/simple_main.f90` (`init_spectre_field`, `trace_orbit_spectre_sympl`).

**Corrected premise** (computer-algebra verified in spectre-orbits
`ca/05_gc_canonical.wl`): the SPECTRE `A_s = 0` gauge supplies only half of the
canonical requirement. The covariant `h_s` does not vanish in SPECTRE
coordinates (metric off-diagonals), so the guiding-center phase-space 1-form
keeps a `rho_par B_s ds` term. The Meiss angle transform `zeta_c` with
`d_s lambda = -B_s/B_zeta` plus gauge `chi` with `d_s chi = (d_s lambda) A_zeta`
restores both `A'_s = 0` and `B'_s = 0`, so the guiding-center flow is canonical
in `(rho_g, theta_c, zeta_c)`.

**Per-volume construction**. The construction ODE would integrate through the
interface current sheets in a global chart, so it is run once per volume on the
smooth slab `r in [lvol-1, lvol]` (`lvol = 1..Mvol`) using the existing Meiss
machinery (`init_meiss` + `get_meiss_coordinates`). The axis volume starts at
`r = 1e-3` so the construction never touches the `r = 0` coordinate singularity.
The SPECTRE chart is its own radial coordinate, so the identity radial scaling is
used (no `r = sqrt(s)`). Each volume is built on the `(r, theta, zeta)` grid set
by the `spectre_ncon_r/th/phi` namelist knob (default `48 x 48 x 32`); after
`get_meiss_coordinates`, `harvest_meiss_volume` moves the field and transform
batch splines (coefficient arrays via `move_alloc`, not copied) into a per-volume
`meiss_volume_t` slot in `field_can_spectre`.

**Construction-grid knob and memory** (#442). The quintic batch splines dominate
peak memory: the coefficient array is `n_quantities * 6^3 * n_r * n_th * n_phi`
per volume, ~1.8 GB at the `48 x 48 x 32` default on the two-volume tok2vol
fixture. `spectre_ncon_r/th/phi` (threaded in by
`set_spectre_construction_grid` before the build, kept as a setter because
`field_can_*` cannot `use params` without a dependency cycle) let a caller trade
resolution for memory. `spectre_ncon_phi = -1` (the default) is automatic: fields
carrying toroidal harmonics keep the historical `n_phi = 32` (bit-identical),
while an **axisymmetric** field (all `in == 0`, e.g. any tokamak) is phi-invariant
and auto-clamps to `AXISYM_NPHI = 8` -- the dominant memory saver, ~3.8x on
tok2vol (1.77 GB -> 0.46 GB), and bit-identical to `n_phi = 32` for well-confined
orbit energy conservation (verified to `2.9e-12`). A positive `spectre_ncon_phi`
overrides the clamp verbatim: raise it to 32 when high-order symplectic
convergence needs the full phi grid, or lower it for RK45-only runs (the RK45
path, `integmode = 0`, builds no construction at all). `AXISYM_NPHI` stays above
the quintic `order + 1 = 6` so the periodic spline is well conditioned.

**Units**. The construction stays in the field's SI covariant units (numerically
balanced, `O(1)`). The Gaussian-CGS conversion is applied only at
`evaluate_spectre_can`, uniformly to values and derivatives because
`(rho_g, theta, zeta)` are all dimensionless: `Bmod` `1e4` (Tesla -> Gauss),
covariant `h_theta`/`h_phi` `1e2` (`L`), covariant `A_theta`/`A_phi`
`1e8 = 1e4 * (1e2)^2` (flux-like `[B] L^2`, Tesla*m^2 -> Gauss*cm^2). The last
factor makes `sqrt(g) = (h_phi A_theta' - h_theta A_phi')/Bmod` come out in
`cm^3`, matching `magfie_spectre`.

**Evaluation**. `evaluate_spectre_can` picks `lvol = min(int(r)+1, Mvol)` --
or the volume pinned by `set_spectre_volume_lock` during multi-volume
stepping, so Newton iterates probing past an interior interface never read
the neighbour volume's discontinuous field -- and evaluates that volume's
splines; `integ_to_ref` / `ref_to_integ` apply the active volume's
`lambda_phi` transform (radial map is the identity). Within `EDGE_BAND` of a
volume edge and beyond it the field is the C1 radial extension from the band
edge (`extend_band`): `A` and `Bmod` linear, the radial slope of `h_theta`,
`h_phi` blended to zero, which removes the sheet-adjacent `Bstar_par -> 0`
degeneracy of the guiding-center 1-form from the zone (the canonical analogue
of the RK45 path's `drift_band`).

**Multi-volume orbits (#441)**. `trace_orbit_spectre_sympl` drives
`spectre_sympl_orbit`: an accepted step leaving `[home_lo, home_hi]` is
rewound and re-solved for the substep length landing on the interface to
`|rho_g - k| < 1e-10` (Illinois false position; every trial is one full
implicit step of the same scheme, so the landed state is on the discrete
orbit of a symplectic map). The landed state is converted to physical
variables, `apply_crossing` maps it across (or reflects), and the integrator
is re-initialized from the mapped state in the target volume's gauge exactly
like an orbit start; the remaining microstep budget is completed there.
Markers are lost only at the outermost interface, matching the RK45 path;
`sympl_rmax` sits at `Mvol + 1` so the internal Newton guards never fire on
iterates inside the domain or on the landing solve at `rho_g = Mvol`.
Grazing events whose landing solve stalls, and committed steps that jump in
`rho_g`, energy, or theta beyond physical bounds (silently unconverged
Newton), are retried with a halved substep; persistent non-convergence
terminates the orbit through the logged `CROSS_STOP` fallback. On tok2vol
this occurs only at full 3.5 MeV alpha energy, where the interior
`Bstar_par = 0` surface invalidates the guiding-center premise for a subset
of pitches.

**Optional follow-up** (documented, not implemented): a `no_K`-style variant that
drops `rho_par h_s` in SPECTRE coordinates directly -- the same
`O(Larmor-radius)` approximation as SIMPLE's Boozer `use_B_r = .false.`. The
exact per-volume Meiss construction here is the reference path.

---

## 7. libneo Integration

SIMPLE depends heavily on libneo for coordinate systems and VMEC data.

### 7.1 Key libneo Imports

```fortran
use libneo_coordinates, only: &
    coordinate_system_t,           &  ! Abstract base
    make_vmec_coordinate_system,   &  ! Create VMEC coords
    make_chartmap_coordinate_system, & ! Create chartmap coords
    detect_refcoords_file_type,    &  ! File type detection
    chartmap_coordinate_system_t,  &  ! Chartmap type
    PSI_TOR_NORM, PSI_POL_NORM,   &  ! Rho conventions
    RHO_TOR, RHO_POL, UNKNOWN
```

### 7.2 Key libneo Functions

| Function | Purpose |
|----------|---------|
| `splint_vmec_data` | Interpolate VMEC spline data |
| `vmec_to_boozer` | VMEC -> Boozer coordinate transform |
| `boozer_to_vmec` | Boozer -> VMEC inverse |
| `vmec_to_can` | VMEC -> Canonical flux |
| `can_to_vmec` | Canonical flux -> VMEC |
| `splint_boozer_coord` | Field in Boozer coordinates |
| `splint_can_coord` | Field in canonical flux coordinates |

### 7.3 VMEC Initialization

```fortran
subroutine init_vmec(vmec_file, ans_s, ans_tp, amultharm, fper)
    call spline_vmec_data()           ! Initialize libneo splines
    call stevvo(RT0, R0i, L1i, ...)   ! Get geometry parameters
    fper = twopi / dble(L1i)          ! Field period
    call volume_and_B00(volume, B00)  ! Get volume and B_00
end subroutine
```

### 7.4 Reference Coordinate Initialization

```fortran
subroutine init_reference_coordinates(coord_input)
    call detect_refcoords_file_type(coord_input, file_type, ierr, message)

    select case (file_type)
    case (refcoords_file_chartmap)
        call make_chartmap_coordinate_system(ref_coords, coord_input)
    case (refcoords_file_vmec_wout)
        call make_vmec_coordinate_system(ref_coords)
    case (refcoords_file_unknown)
        error stop "Unknown coordinate file type"
    end select
end subroutine
```

---

## 8. Configuration Options

### 8.1 Field Type Selection (`isw_field_type` / `integ_coords`)

| Value | Name | Coordinates | Integrators |
|-------|------|-------------|-------------|
| -1 | TEST | (r, theta, phi) | Symplectic only |
| 0 | CANFLUX | (s, theta_c, phi_c) | Symplectic |
| 1 | VMEC | (s, theta, phi) | RK45 only |
| 2 | BOOZER | (s, theta_B, phi_B) | Symplectic |
| 3 | MEISS | (r, theta_c, phi_c) | Symplectic |
| 4 | ALBERT | (psi, theta_c, phi_c) | Symplectic |
| 5 | REFCOORDS | (varies) | RK45 only |
| 6 | SPECTRE | (rho_g, theta, zeta) | RK45 (integmode = 0) or per-volume symplectic (integmode > 0) |

`integ_coords` is the integer mode id above; `6` selects SPECTRE. A string
alias `'spectre'` is not read from the namelist (the slot is integer).

### 8.2 Input File Parameters

| Parameter | Type | Purpose |
|-----------|------|---------|
| `netcdffile` | string | VMEC wout.nc file path |
| `field_input` | string | Field data file (VMEC wout, coils file, or Boozer chartmap NetCDF) |
| `coord_input` | string | Reference coordinate file (VMEC or chartmap) |

### 8.2.1 GVEC -> Boozer chartmap -> SIMPLE

This is the supported user workflow for using GVEC output with the new Boozer
chartmap path.

SIMPLE does **not** consume raw GVEC `parameter/state` files at runtime. First
convert them to an extended Boozer chartmap NetCDF:

```bash
python tools/gvec_to_boozer_chartmap.py \
    parameter_final.ini \
    State_final.dat \
    boozer_chartmap.nc
```

Then use that file as both field source and reference coordinates:

```fortran
&config
field_input = 'boozer_chartmap.nc'
coord_input = 'boozer_chartmap.nc'
isw_field_type = 2
startmode = 2
/
```

Notes:
- `startmode = 2` means particle starts are already given in the chartmap
  Boozer coordinates.
- If `field_input` has the `boozer_field = 1` metadata, SIMPLE takes the
  chartmap runtime path and does not initialize VMEC.
- The tested reference workflow is `test/tests/run_gvec_qa_roundtrip.py`.

### 8.2.2 STL Wall Intersection (`wall_input`, `wall_units`)

SIMPLE can optionally detect wall losses by intersecting the orbit segment in
Cartesian space with a triangulated STL surface (using CGAL).

| Parameter | Type | Purpose |
|-----------|------|---------|
| `wall_input` | string | STL file path. If empty, wall intersections are disabled. |
| `wall_units` | string | STL units: `m` or `mm` (STL has no unit metadata). |

**Requirements**:
- Must be built with `-DSIMPLE_ENABLE_CGAL=ON`.
- Requires `coord_input` to be a chartmap reference coordinate file.

**Performance gating**:
- The expensive STL intersection is only evaluated for `rho > rho_lcfs`, where
  `rho_lcfs` is read from the chartmap NetCDF metadata.

**Outputs** (when `output_results_netcdf = .True.`):
- `results.nc:wall_hit(particle)` (byte)
- `results.nc:wall_hit_cart(xyz, particle)` (cm, consistent with `xstart_cart`)

**Fallback Chain**:
1. coord_input (if set)
2. netcdffile (if coord_input empty)
3. field_input (if both empty)

### 8.2.3 SPECTRE Equilibria (`integ_coords = 6`)

VMEC-free multi-volume guiding-center tracing on a SPECTRE/SPEC MRxMHD
equilibrium (HDF5, `.h5`), detected automatically from `field_input`. The
end-to-end pipeline (issues #437-#442): symplectic field-line Poincare app
(#437), per-volume RK45 guiding centers (#438), per-volume Meiss canonical
coordinates (#439), Level-1 refraction (#440) and Level-0 rescale (#443)
interface crossing maps, symplectic multi-volume crossing (#441), and the
validation suite (#442). Example input: `examples/simple_spectre.in`.

| Parameter | Type | Default | Purpose |
|-----------|------|---------|---------|
| `field_input` | string | `''` | SPECTRE HDF5 equilibrium (`.h5`); VMEC-free |
| `integ_coords` | int | -1000 | `6` selects the SPECTRE per-volume stacked-rho chart |
| `integmode` | int | 0 | `0` = RK45 per volume; `> 0` = per-volume symplectic (builds the canonical construction) |
| `crossing_level` | int | 1 | interface crossing map: `1` = Level-1 refraction, `0` = Level-0 energy rescale |
| `spectre_ncon_r` | int | 48 | Meiss construction grid points in `r` (per volume) |
| `spectre_ncon_th` | int | 48 | construction grid points in `theta` |
| `spectre_ncon_phi` | int | -1 | construction grid points in `phi`; `-1` = automatic (32 with toroidal harmonics, `AXISYM_NPHI` for axisymmetric fields); a positive value forces that count |

The radial start `sbeg` is the SPECTRE label `rho_g in [0, Mvol]`; markers are
lost only at the outermost interface `rho_g = Mvol`. `crossing_level` and the
construction grid apply to `integmode > 0`; the RK45 path shares
`crossing_level` but builds no construction. On the axisymmetric two-volume
tok2vol fixture Level 0 and Level 1 give equal loss fractions (finding F5); the
committed reference behavior (loss accounting, L0-vs-L1, trapped/passing mirror
split, fixed-seed determinism) is locked by
`test/tests/test_spectre_validation.py`, and the construction-grid knob and
axisymmetric clamp by `test/tests/test_spectre_construction_grid.py`.

### 8.3 Integration Mode (`integmode`)

| Value | Name | Order | Speed |
|-------|------|-------|-------|
| 0 | RK45 | Adaptive | Slow |
| 1 | EXPL_IMPL_EULER | 2 | Fast |
| 2 | IMPL_EXPL_EULER | 2 | Fast |
| 3 | MIDPOINT | 2 | Medium |
| 4 | GAUSS1 | 2 | Slow |
| 5 | GAUSS2 | 4 | Medium |
| 6 | GAUSS3 | 6 | Medium-Slow |
| 7 | GAUSS4 | 8 | Slow |
| 15 | LOBATTO3 | 6 | Slow |

**Constraints**:
- integmode >= 1: Requires canonical field type (TEST, CANFLUX, BOOZER, MEISS, ALBERT, SPECTRE)
- integmode = 0 (RK45): works with any initialized field type (VMEC, REFCOORDS, or canonical modes); `simple_main.f90` guards only integmode >= 1
- SPECTRE multi-volume crossing (#441) admits one-step schemes only (all of
  the above, driven with `ntau = 1`, asserted in `spectre_sympl_orbit`):
  re-canonicalization after each interface crossing restarts the scheme from
  a single point, so no multistep history may carry over

### 8.4 Spline Parameters

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `ns_s` | 5 | Spline order for s dimension |
| `ns_tp` | 5 | Spline order for theta/phi |
| `multharm` | 5 | Angular grid multiplier |

### 8.5 Valid Configurations

**Fast Testing**:
```fortran
integmode = 1, isw_field_type = 2  ! Euler + Boozer
multharm = 3, ns_s = 3, ns_tp = 3
```

**Production**:
```fortran
integmode = 6, isw_field_type = 2  ! Gauss3 + Boozer
multharm = 5, ns_s = 5, ns_tp = 5
npoiper2 = 256
```

**Coil-Based**:
```fortran
integmode = 1, isw_field_type = 3  ! Euler + Meiss
field_input = 'coils.simple'
netcdffile = 'wout.nc'  ! Reference VMEC
```

---

## 9. Coordinate Transformation Workflows

### 9.1 Particle Initialization (VMEC to Chartmap)

When particles start in VMEC coordinates but integration uses chartmap:

```fortran
subroutine vmec_to_chartmap(z_vmec, z_chart)
    ! 1. Get Cartesian position from VMEC coords
    call ref_coords%evaluate_cart(z_vmec(1:3), x_target)

    ! 2. Initial guess (same as VMEC)
    z_chart(1:3) = z_vmec(1:3)

    ! 3. Newton iteration to find chartmap coords
    do iter = 1, max_iter
        call chartmap_coords%evaluate_cart(z_chart(1:3), x_test)
        err_vec = x_test - x_target
        if (norm(err_vec) < tol) exit

        ! Compute Jacobian and update
        call compute_jacobian(z_chart(1:3), J)
        z_chart(1:3) = z_chart(1:3) - J_inv * err_vec
    end do
end subroutine
```

### 9.2 Direct VMEC vs Map2disc Chartmaps

Both chartmap types represent the same boundary but use different interiors:

| Aspect | Direct VMEC | Map2disc |
|--------|-------------|----------|
| Interior | VMEC flux surfaces | Conformal disc |
| rho meaning | sqrt(s) | Geometric radius |
| rho at s=0.8 | 0.894 | ~0.874 |
| Axis position | Magnetic axis | Geometric center |
| rho_convention | rho_tor | unknown |

If particles are initialized at the same physical position and the field is
sampled from the same underlying equilibrium, the orbit physics should agree.

### 9.3 Metric Tensor Differences

The metric `g_ij` can differ between chartmaps at the same nominal
`(rho, theta, zeta)`:

```
At rho=0.5, theta=0, zeta=0:
  Direct VMEC: g_rr = 7820, sqrt(g) = 3.17
  Map2disc:    g_rr = 6695, sqrt(g) = 3.02  (~5% difference)
```

This changes coordinate-specific quantities, but not the intended physical
trajectory.

### 9.4 Orbit Trajectory Output

With `output_orbits_macrostep = .true.`, `orbits.nc` stores the trajectory in
reference coordinates as `s`, `theta`, and `phi`. The global
`radial_coordinate` attribute gives the meaning of the first variable: `s` for
VMEC and `rho` for chartmaps. The variable name `s` remains for file-format
compatibility.

When a reference coordinate system is available, the file also contains `R`
and `Z`. SIMPLE computes them through the active reference geometry, so VMEC
and chartmap runs use the same runtime map as the orbit calculation. Each
variable carries its native length unit (`cm` or `m`). Points without a finite
reference position remain NaN. Analytic test fields do not write `R` and `Z`
because they have no reference-coordinate geometry.

---

## 10. Key Files Reference

### 10.1 Coordinate Systems

| File | Purpose |
|------|---------|
| `libneo/src/coordinates/libneo_coordinates.f90` | Abstract base + VMEC + chartmap |
| `src/coordinates/reference_coordinates.f90` | Runtime coord system selection |
| `src/coordinates/coordinate_scaling.f90` | sqrt_s and identity scaling |
| `src/coordinates/cartesian_coordinates.f90` | Cartesian coordinate system |
| `src/coordinates/coordinates.f90` | Transformation function pointers |

### 10.2 Field Implementations

| File | Purpose |
|------|---------|
| `libneo/src/field/magnetic_field_base.f90` | Abstract magnetic_field_t |
| `libneo/src/field/field_vmec.f90` | VMEC field |
| `src/field/field_coils.f90` | Coils Biot-Savart (62 lines) |
| `src/field/field_boozer_chartmap.f90` | Extended Boozer chartmap field |
| `src/field/field_splined.f90` | Splined decorator (537 lines) |
| `src/field.F90` | Field factory (182 lines) |

### 10.3 Canonical Coordinates

| File | Purpose |
|------|---------|
| `src/field/field_can_base.f90` | field_can_t data structure (65 lines) |
| `src/field/field_can_boozer.f90` | Boozer canonical (129 lines) |
| `src/field/field_can_meiss.f90` | Meiss canonical (659 lines) |
| `src/field/field_can_albert.f90` | Albert canonical (334 lines) |
| `src/field/field_can_flux.f90` | Flux canonical (184 lines) |
| `src/field/psi_transform.f90` | r->psi regridding |

### 10.4 magfie Interface

| File | Purpose |
|------|---------|
| `src/magfie.f90` | Unified field evaluation interface |
| `src/magfie_can_boozer.f90` | Boozer/Canflux implementations |

### 10.5 Integration

| File | Purpose |
|------|---------|
| `src/orbit_symplectic_base.f90` | Integrator types and RK coefficients |
| `src/orbit_symplectic.f90` | Symplectic methods |
| `src/orbit_symplectic_quasi.f90` | Quasi-symplectic and RK45 |
| `src/alpha_lifetime_sub.f90` | orbit_timestep_axis |

---

## Appendix: Known Issues and Complications

### A.1 Chartmap Boundary Range

**Issue**: VMEC coordinate inversion fails at exact boundary (rho=1.0).
**Solution**: Spline range limited to [0.1, 0.99].
**Impact**: Particles beyond rho=0.99 use extrapolated field values.

### A.2 Coordinate Scaling Complexity

**Issue**: Multiple coordinate conventions (s vs rho vs r) require careful
handling of Jacobians and component transformations.
**Complication**: Easy to miss scaling factors when adding new coordinate
systems or field types.

### A.3 Single-Instance Canonical Coordinates

**Issue**: Meiss/Albert use module-level variables for splines and coordinate data.
**Note**: This is not a thread-safety issue since variables are written only at
initialization, then read-only during parallel particle tracing.
**Limitation**: Cannot have multiple canonical field instances with different
configurations simultaneously.

### A.4 VMEC from_cyl Newton Iteration

**Issue**: Newton iteration can fail near separatrix or for points outside
plasma.
**Workaround**: Error checking and fallback to cylindrical coordinates.

### A.5 Inconsistent sqgBctr Support

**Issue**: Some field types support sqgBctr output, others do not.
- vmec_field_t: Not implemented
- coils_field_t: Returns raw Cartesian B
- gvec_field_t: Fully implemented
- splined_field_t: Not implemented

### A.6 Mode Naming Confusion

**Issue**: Parameter `isw_field_type` is deprecated but still widely used.
New name `integ_coords` is more descriptive but less documented.
**Recommendation**: Use `integ_coords` in new code.

---

*Document last updated: 2025-12-29*
*Corresponding code version: Post v1.5.1*
