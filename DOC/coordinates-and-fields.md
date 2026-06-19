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
**File**: `src/field/field_vmec.f90`

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

### 3.4 GVEC Field

**Type**: `gvec_field_t`
**File**: `src/field/field_gvec.f90`

- Reads `.dat` files with B-spline magnetic field data
- Requires `GVEC_AVAILABLE` compile flag
- This is a separate runtime path from the chartmap-based Boozer path

For the workflow added in this PR, SIMPLE does **not** read raw GVEC state
files directly. Instead, GVEC output is converted to an extended Boozer
chartmap NetCDF and SIMPLE runs from that chartmap file.

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

### 6.6 The Curvilinear 6D Canonical Integrator

**Files**: `src/orbit_cpp_canonical.f90`, `src/orbit_cpp_vmec_metric.f90`,
`src/orbit_cpp_chartmap_metric.f90`, `src/field/field_can_test.f90`
(`eval_field_correct_test`)

The guiding-center integrators reduce the perpendicular motion to the magnetic
moment. The 6D canonical integrator in `orbit_cpp_canonical` keeps the full
phase space `(q, p)` and resolves (or, for the Pauli models, represents) that
motion directly. It is the SIMPLE port of the Egger-Feiel thesis
discrete-variational integrators, generalized to arbitrary curvilinear
coordinates with a full (non-diagonal) metric.

The Hamiltonian is `H = (1/2m)(p_i - qc A_i) g^ij (p_j - qc A_j) [+ mu|B|]`, so
`q_dot^k = (1/m) g^kj (p_j - qc A_j)` and
`p_dot_k = qc A_{j,k} v^j + (m/2) g_{ij,k} v^i v^j [- mu |B|_{,k}]`. Every term
carries the full metric `g_ij`, its inverse `g^ij`, and the direction
derivatives `g_{ij,k}`. The integrator reads them from a `block_t`: metric,
metric derivatives, covariant `A_i` with gradient, `|B|` with gradient, and the
covariant unit field `h_i`.

Three coordinate blocks fill that structure.

`COORD_TOK` is the analytic tokamak, inline and GPU-portable. The metric is
diagonal, `g = diag(1, r^2, (R0 + r cos theta)^2)`,
`sqrt(g) = r (R0 + r cos theta)`. The covariant vector potential is `A_r = 0`,
`A_theta = B0 (r^2/2 - r^3 cos(theta)/(3 R0))`,
`A_phi = -B0 iota0 (r^2/2 - r^4/(4 a^2))`, with `B0 = iota0 = R0 = 1`, `a = 0.5`.
The 6D path needs the exact field from the curl of `A`:
`B^k = eps^ijk A_{j,i} / sqrt(g)`, `|B| = sqrt(g_ij B^i B^j)`, so with `A_r = 0`,
`|B|^2 = A_{phi,r}^2 / (R0 + r cos theta)^2 + A_{theta,r}^2 / r^2`
(`eval_field_correct_test`). The guiding-center path keeps the linearized
`eval_field_test` (`|B| = B0(1 - r/R0 cos theta)`); at the reference start
`(r,theta) = (0.1, 1.5)` the exact `|B| = 0.99749` against the linearized
`0.99293`. The diagonal metric is the special case of the general arithmetic
(off-diagonals zero), so `COORD_TOK` reproduces the python oracle bit-for-bit
while the same residual runs on a stellarator metric.

`COORD_VMEC` runs on real VMEC equilibria in native flux coordinates
`(s, vartheta, varphi)`, wired through `orbit_cpp_vmec_metric`. It is the
production 6D loss chart (see below). The full metric `g_ij`, `g^ij` and
Christoffel symbols `Gamma^l_jk` come from libneo's `coordinate_system_t`
(issue #322, branch `feature/metric-christoffel`); the metric derivatives follow
from metric compatibility, `g_{ij,k} = g_il Gamma^l_jk + g_jl Gamma^l_ik`. The
covariant `A_i` and `|B|` come from SIMPLE's native VMEC field
(`vmec_field_evaluate`), with `dA` and `d|B|` by central difference. The metric
is consistent with the field: the covariant unit field obeys
`h_i g^ij h_j = |h|^2 = 1` to central-difference accuracy (the libneo metric
gives `g g^-1 = I` to `~1e-15`; the FD Christoffel sets the residual `~1e-2`).
This block is host-side: libneo's metric is `class()`-dispatched and reads 3D
splines, so it cannot run under `!$acc routine seq`.

`COORD_CHARTMAP` runs the 6D state in `(rho, theta_B, phi_B)` with
`rho = sqrt(s)`, wired through `orbit_cpp_chartmap_metric`: the chartmap metric
and Christoffel from `reference_coordinates%ref_coords`, the field reparametrized
from `s = rho^2` with the radial chain rule `dF/drho = 2 rho dF/ds`, and the
covariant `A_i`, `h_i`, `|B|` from the active `field_can_mod%evaluate` pointer.
This chart is NOT consistent and is not the production route. libneo splines the
chartmap Cartesian `x/y/z` with a periodic fit over one field period, but for
`nfp > 1` the Cartesian `x,y` are not field-period-periodic (they rotate by
`2pi/nfp`), so the periodic spline destroys the secular toroidal rotation: the
analytic spline `e_phi` loses its `~R` magnitude and the geometric metric gives
`h_i g^ij h_j ~ nfp^2` instead of 1. The defect is upstream in libneo's
Cartesian-storage path and in the storage convention itself, so it cannot be
repaired in the SIMPLE metric post-processor. A consistent chartmap route needs
an R,Z (cylindrical) Boozer-chartmap representation in libneo: R and Z are
field-period-periodic, and the reader's `has_spl_rz` path already adds the
toroidal rotation analytically. Until then the production 6D loss path uses
`COORD_VMEC`.

The magnetic coupling carries a length normalization. The Hamiltonian uses
`qc = charge/(c rho0)` with `rho0 = 1` by default, so `COORD_TOK` keeps the
thesis `qc = charge/c` (`c = 1`). The production wire threads
`rho0 = ro0_bar = ro0/sqrt(2)` so `qc = sqrt(2)/ro0`, which makes the canonical
momentum `p_i = vpar h_i + A_i/ro0_bar` match the guiding-center `pphi` seed of
`init_sympl`.

Three models share one integer-dispatched residual: `MODEL_CP` (full charged
particle), `MODEL_CPP_SYM` (Pauli symplectic midpoint, `H + mu|B|`),
`MODEL_CPP_VAR` (Pauli variational midpoint, discrete Euler-Lagrange). `MODEL_CP`
omits the `mu|B|` term: the full charged particle resolves the perpendicular
gyromotion directly, so its kinetic energy already holds the perpendicular
energy; the Pauli models drop the resolved gyromotion and reinstate it as the
guiding-center `mu|B|`. The state is fixed-size 6,
`z = (q1, q2, q3, p1, p2, p3)`; the position rows solve the canonical midpoint
and the momentum rows carry `p`, so the Jacobian is square `6x6` and solved with
the device LU `rk_solve` from `orbit_rk_core`. For `COORD_TOK` Newton uses the
analytic Jacobian: `d2g` and `d2A` come from central differences of the block's
own `dg`/`dA`, and the `O(mu)` `|B|` force gradient uses the block's analytic
Hessian `d2Bmod`, the closed-form second derivative of the corrected `|B|`. The
analytic-vs-finite-difference self-check passes for all three models. For
`COORD_VMEC` the Jacobian is a central difference of the whole residual,
consistent with the spline-based block; the FD step is taken relative to each
variable's own scale (angles `O(1)`, momenta their own magnitude) so the
physical-CGS momenta `~1e-8` keep an accurate column. A central-difference
Jacobian is accurate to `~1e-7`, so the Newton step cannot reach the
analytic-path floor `rtol = 1e-12`; the `COORD_VMEC` path converges on a
step tolerance `rtol_fd = 1e-8` while `COORD_TOK`/`COORD_CHARTMAP` keep the
tight `rtol`.

The field `|B|` and its derivatives are the exact closed form of
`|B| = sqrt(W)`, `W = A_phi,r^2/(R0 + r cos theta)^2 + A_theta,r^2/r^2`:
`d_k|B| = dW_k/(2|B|)`, `d2_kl|B| = d2W_kl/(2|B|) - dW_k dW_l/(4|B|^3)`,
finite-difference verified to `~1e-9` against `|B|`. An earlier port carried over
two errata from the python listing. The metric theta-derivative now carries the
factor `r`: `d g_33/d theta = -2 r (R0 + r cos theta) sin theta`. The field
`d|B|/d theta` now keeps the `A_theta,rth` chain-rule term the listing dropped.
With both corrected, the `CPP-sym` energy oscillation over a fixed time window
converges as `dt^2`: `max|dE/E0| = 4.2e-5, 1.0e-5, 2.2e-6, 4.8e-7` for
`dt = 80, 40, 20, 10`, each halving reducing the bound by about four. The buggy
`d|B|` instead produced a flat `~1e-3` plateau that did not converge; the test
asserts the `dt^2` behavior. The Fortran reproduces the regenerated python
oracle for all three models to solver tolerance.

`COORD_TOK` runs on the OpenACC device. `cpp_canon_step_tok` is the device entry:
the whole chain (`residual_tok` -> `eval_block_tok` /
`eval_field_correct_test` / `dLdq` / `raise` / `residual_blk`,
`jacobian_analytic` -> `grad_jacobian_tok`, `rk_solve`) is `!$acc routine seq`
with fixed-size state, integer model dispatch, the analytic Jacobian, and no
`class()` or procedure pointer, so one particle runs per GPU thread. The host
`cpp_canon_step` keeps the coordinate dispatcher; `COORD_VMEC` is host-only
because libneo's metric is `class()`-dispatched and reads 3D splines, which
cannot run under `!$acc routine seq`. `test/tests/test_cpp_canonical_device.f90`
(built only with `SIMPLE_ENABLE_OPENACC=ON`) runs all three models for a batch of
particles inside an `!$acc parallel loop` and checks the device result against
the host step.

`test/tests/test_cpp_canonical.f90` validates the analytic block against the
regenerated python oracle. `test/tests/test_cpp_vmec.f90` runs the same
integrator on `test/test_data/wout.nc` (a 2-field-period stellarator): the
libneo metric satisfies `g g^-1 = I` to machine precision, `CP` energy stays
bounded with no secular drift, and the big-step `CPP` orbit stays on a bounded
radial band with radial bounce points, the guiding-center confinement signature.
The stellarator is not axisymmetric, so the toroidal canonical momentum is not
conserved and is not asserted; near the axis `s -> 0` the flux metric is singular
and the central-difference gradients lose accuracy, so the test starts at
mid-radius.

The genuine 6D canonical CPP is wired into the production alpha-loss pipeline as
`orbit_model = ORBIT_CPP6D` (5), through `COORD_VMEC`. `init_cpp` in `simple.f90`
replicates the `init_sympl` sqrt(2) block (`mu` by factor 2,
`ro0_bar = ro0/sqrt(2)`, `vpar_bar = vpar sqrt(2)`), reading `|B|` and `h_i` from
the native VMEC field at the start, then seeds the state at `(s, theta, phi)`
(s direct, no rho) with `MODEL_CPP_SYM`, `mass = 1`, `dt = dtaumin/sqrt(2)`, and
`rho0 = ro0_bar`. Keeping `mass = 1` is what makes the wire well-posed: with the
consistent `|h|^2 = 1` metric the kinetic term `(1/2m)(p-qcA)g^ij(p-qcA)` along
the field reduces to `vpar_bar^2/2`, so the 6D Hamiltonian equals the GC one and
the velocities stay `O(vpar_bar)` (physical-CGS `mass ~ 1e-24` would blow up
`v^i = g^ij(...)/m` and wreck the Newton). The covariant momenta
`p_theta = vpar h_theta + A_theta/ro0_bar`, `p_phi = vpar h_phi + A_phi/ro0_bar`
match the GC seed; `p_s` carries the `g_si v^i` metric term, the genuine 6D
start. `orbit_timestep_cpp_canonical` advances one `dtaumin/sqrt(2)` step and
writes `z(1:5)` itself (`z(1) = s` for `COORD_VMEC`, `z(1) = rho^2` for
`COORD_CHARTMAP`, angles direct, `z(4) = pabs`,
`z(5) = vpar_bar/(pabs sqrt(2))`), so `times_lost`, `confined_fraction`, and the
trajectory output read `z(1:5)` exactly as on the GC path. The macrostep in
`simple_main` dispatches on `orbit_model`; the GC default and `ORBIT_PAULI` keep
`to_standard_z_coordinates`, `ORBIT_CPP6D` routes around it. `ORBIT_CPP6D`
requires a VMEC-backed canonical field (checked once in `main`, rejecting a
standalone Boozer-chartmap input), needs `swcoll = .false.` and no wall, and
error-stops in classification (`ntcut > 0` / `class_plot`); collisions, the wall
path, and the classifier stencil are the documented follow-ups. The COORD_VMEC
metric is attached once before the parallel region so per-thread `init_cpp` never
races on the allocation.

`test/tests/test_cpp6d_vs_gc.f90` drives the production `init_field` (BOOZER on
the real `test/test_data/wout.nc`), then runs the production GC
(`orbit_timestep_sympl`) and the production 6D wrapper (`COORD_VMEC`) from the
same start over 2000 bare GC macrosteps. It asserts the chart consistency
(`h_i g^ij h_j ~ 1`, the check the chartmap fails), `mu` held fixed and energy
bounded, both orbits confined with overlapping radial (`s`) bands, and the
`s >= 1` loss propagating through the wrapper. End to end, a 32-particle loss run
on the same equilibrium gives a CPP6D confined fraction within `O(rho*)` of the
GC one (`0.91` vs `0.97` at `trace_time = 1e-3 s`), the genuine-6D
finite-Larmor difference over the GC.

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
| 5 | REFCOORDS | (varies) | Both |

### 8.2 Input File Parameters

| Parameter | Type | Purpose |
|-----------|------|---------|
| `netcdffile` | string | VMEC wout.nc file path |
| `field_input` | string | Field data file (VMEC, coils, or GVEC) |
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
- integmode >= 1: Requires canonical field type (TEST, CANFLUX, BOOZER, MEISS, ALBERT)
- integmode = 0 (RK45): Requires VMEC field type

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
| `src/field/field_base.f90` | Abstract magnetic_field_t (28 lines) |
| `src/field/field_vmec.f90` | VMEC field (72 lines) |
| `src/field/field_coils.f90` | Coils Biot-Savart (62 lines) |
| `src/field/field_gvec.f90` | GVEC B-spline (214 lines) |
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
| `src/orbit_rk_core.f90` | Shared device LU and Newton shell |
| `src/orbit_cpp_canonical.f90` | Curvilinear 6D canonical-midpoint integrator (cp/cpp_sym/cpp_var) |
| `src/orbit_cpp_vmec_metric.f90` | VMEC metric/Christoffel + native field provider for the 6D integrator |
| `src/orbit_cpp_chartmap_metric.f90` | Production Boozer/chartmap metric + field_can provider for the 6D integrator (ORBIT_CPP6D) |
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
