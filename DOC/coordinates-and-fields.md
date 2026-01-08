# Coordinate Systems and Magnetic Fields in SIMPLE

This document provides comprehensive documentation of the coordinate systems,
magnetic field representations, and their interactions in SIMPLE. It covers
the complete architecture from abstract interfaces to concrete implementations.

**IMPORTANT**: This document must be revised whenever coordinates or fields
are changed or refactored. See CLAUDE.md for the maintenance requirement.

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

**Coordinates**: (rho, theta, zeta)
- `rho`: Normalized radial, 0 <= rho <= 1
- `theta`: Poloidal angle, 0 to 2pi
- `zeta`: Toroidal angle, 0 to 2pi/nfp

**Source**: NetCDF file with Cartesian mapping on structured grid

**File Format**:
```
Dimensions: rho, theta, zeta
Variables: x(zeta, theta, rho), y(...), z(...)
Attributes: num_field_periods, zeta_convention, rho_convention
```

**Rho Conventions**:
- `RHO_TOR`: rho = sqrt(s_toroidal)
- `RHO_POL`: rho = sqrt(s_poloidal)
- `PSI_TOR_NORM`: rho = s (normalized toroidal flux)
- `PSI_POL_NORM`: rho = s (normalized poloidal flux)
- `UNKNOWN`: Geometric disc mapping (e.g., map2disc)

**Two Types of Chartmaps**:

1. **Direct VMEC Chartmap**: Generated by `vmec_to_chartmap.x`
   - rho = sqrt(s), directly samples VMEC coordinates
   - Interior follows VMEC flux surfaces
   - `rho_convention = rho_tor`

2. **Map2disc Chartmap**: Generated by `generate_chartmap_map2disc.py`
   - Uses conformal disc mapping from boundary only
   - Interior is geometric (not flux-based)
   - `rho_convention = unknown`
   - Same boundary as VMEC, different interior
   - Suitable as a reference coordinate system for Meiss canonicalization;
     Meiss does not require flux coordinates, only a smooth reference mapping
     with non-vanishing toroidal field component in the domain.
   - For orbit integration stability, default splined-field resolution is
     increased for `rho_convention=unknown` chartmaps (currently
     96×97×96 unless overridden), and the sampling range remains
     `rho ∈ [0.1, 0.99]` to avoid VMEC inversion failures at the boundary.

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
- Complex Hamiltonian map transformation for evaluation
- Requires `GVEC_AVAILABLE` compile flag

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
    ! 1. Wrap chartmap angles to their periodic domains and evaluate splined field
    call refcoords_field%evaluate_with_der(x_eval, Acov, hcovar, Bmod,
                                           dAcov, dhcov, dBmod)

    ! 2. Handle coordinate scaling (sqrt_s vs identity)
    call scaled_to_ref_coords(refcoords_field%coords, x_eval, u_ref, J)

    ! 3. Get metric from reference coordinates
    call refcoords_field%coords%metric_tensor(u_ref, g, ginv_u, sqrtg_abs)
    call refcoords_field%coords%covariant_basis(u_ref, e_cov)

    ! 4. Compute signed Jacobian
    sqrtg = signed_jacobian(e_cov) * J

    ! 4b. Apply safety floors to avoid divide-by-zero in downstream steps
    bmod = max(abs(bmod), 1.0d-30)
    sqrtg = sign(max(abs(sqrtg), 1.0d-14), sqrtg)

    ! 5. Transform to contravariant components
    call inverse_metric_scaled(J, ginv_u, ginv_x)
    hctrvr = matmul(ginv_x, hcovar)

    ! 6. Compute curl
    call compute_hcurl(sqrtg, dhcov, hcurl)
end subroutine
```
For chartmap reference coordinates, `magfie_refcoords` wraps `theta` to
`[0, 2*pi)` and `zeta` to `[0, 2*pi/nfp)` before evaluating splines, clamps
`rho` to the spline domain, and floors `bmod`/`sqrtg` to avoid divide-by-zero
in `compute_hcurl`.

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
- The toroidal angle is treated with the field-period range
  `phi in [0, 2*pi/nfp)` during spline construction and wrapping.

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

Both chartmap types cover the same physical boundary but differ inside:

| Aspect | Direct VMEC | Map2disc |
|--------|-------------|----------|
| Interior | VMEC flux surfaces | Conformal disc |
| rho meaning | sqrt(s) | Geometric radius |
| rho at s=0.8 | 0.894 | ~0.874 |
| Axis position | Magnetic axis | Geometric center |
| rho_convention | rho_tor | unknown |

**Key Insight**: Despite different interior mappings, both give identical
physics because:
1. Particles start at same physical location (via Newton inversion)
2. B field sampled from VMEC at physical locations
3. Loss boundary (rho=1) is same physical surface
4. Orbit equations are coordinate-covariant

### 9.3 Metric Tensor Differences

The metric g_ij differs between chartmaps at same (rho, theta, zeta):

```
At rho=0.5, theta=0, zeta=0:
  Direct VMEC: g_rr = 7820, sqrt(g) = 3.17
  Map2disc:    g_rr = 6695, sqrt(g) = 3.02  (~5% difference)
```

This affects coordinate-specific quantities but not physics.

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
