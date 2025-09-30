# TODO: GEQDSK Support with Geoflux Reference Coordinates

## Overview

Add tokamak field support from GEQDSK files using new "geoflux" reference coordinates:
- **Radial**: s_tor = ψ_tor/ψ_tor,edge (normalized toroidal flux computed from ψ_pol and q-profile)
- **Poloidal**: θ_geo = atan2(Z - Z_axis, R - R_axis) (geometric angle from magnetic axis)
- **Toroidal**: φ = φ_cyl (cylindrical toroidal angle)

**Architecture**: Follow VMEC pattern with clean separation between libneo (reference coordinates) and SIMPLE (canonical coordinates).

**Key Simplification**: Reuse libneo's existing `field_divB0` module for all GEQDSK I/O and field evaluation. We only need to add:
1. Coordinate transformations: (s_tor, θ_geo, φ) ↔ (R, φ, Z)
2. Toroidal flux computation: s_tor(ψ_pol) from q-profile integration
3. Flux surface cache: Pre-compute (R, Z) on (s_tor, θ_geo) grid for performance

This significantly reduces implementation effort compared to reimplementing GEQDSK handling from scratch.

---

## Summary: 3D Field Superposition Capability

### What field_divB0 Provides

The `field_divB0` module in libneo supports **3D non-axisymmetric perturbations** on top of axisymmetric equilibria through the `ipert` parameter:

| Mode | Description | Use Case |
|------|-------------|----------|
| `ipert=0` | Axisymmetric equilibrium only | Basic tokamak physics, initial GEQDSK implementation |
| `ipert=1` | Vacuum perturbation (Biot-Savart) | RMPs, error fields, test coils (no plasma response) |
| `ipert=2` | Vacuum + plasma response | Linear MHD response (fast, no derivatives) |
| `ipert=3` | Vacuum + plasma response + derivatives | Full linear MHD (7× slower) |

**Key insight**: Reference coordinates remain axisymmetric (based on equilibrium), while the magnetic field includes 3D perturbations.

### How It Works

**Vacuum Perturbation Mode (ipert=1)**:
1. Pre-compute 3D coil field on cylindrical (R, φ, Z) grid using **vacfield.x**
2. Store in `pfile` (format controlled by `icftype`)
3. `field_divB0` converts to divergence-free Fourier representation (up to `ntor` harmonics)
4. At runtime: `B_total = B_equilibrium + ampl * B_perturbation`

**Generating Coil Fields**:
```bash
# Use libneo's vacfield.x program
vacfield.x AUG 2 coil_file1.dat coil_file2.dat Bvac grid.inp output.h5
```

Supported coil formats:
- **AUG**: ASDEX Upgrade format
- **GPEC**: General Perturbation Equilibrium Code
- **Nemov**: Wendelstein 7-X format
- **STELLOPT**: Filament format (coils.c09r00)

**Plasma Response Mode (ipert=2,3)**:
- Requires pre-computed linear MHD response in flux coordinates
- Typically from GPEC, MARS-F, or similar codes
- Reads Fourier harmonics from `fluxdatapath`
- Uses `field_fourier` / `field_fourier_derivs` for evaluation

### Tools Available

**Fortran**:
- `coil_tools.f90`: Coil I/O, Biot-Savart evaluation
- `coil_convert.x`: Convert between coil file formats
- `vacfield.x`: Generate field files from coil geometries
- `bdivfree.f90`: Divergence-free field representation

**Python**:
- `libneo.coils`: Read/write coil files, metadata extraction
- `libneo.mgrid`: STELLOPT M-grid format handling

### Current Status

**✅ Fully Implemented**:
- All 4 modes working in production codes
- Coil geometry readers for multiple formats
- Biot-Savart field computation
- Fourier decomposition and evaluation

**❌ Missing**:
- No automated regression tests (all test files have `ipert=0`)
- No example workflows documented
- No validation suite for 3D modes

**📋 For SIMPLE+GEQDSK**:
1. **Phase 1 (Priority 1)**: Implement pure equilibrium (`ipert=0`)
   - Focus on coordinate transformations and canonical coordinates
   - Get basic tokamak orbit tracing working
2. **Phase 2 (Future)**: Add RMP capability (`ipert=1`)
   - Interface SIMPLE with pre-computed coil fields
   - Enable RMP and error field studies
   - Requires: grid generation, `pfile` format support
3. **Phase 3 (Advanced)**: Full plasma response (`ipert=2,3`)
   - Couple to GPEC or MARS-F outputs
   - Study self-consistent perturbed equilibria
   - Out of scope for initial implementation

### Integration References

**In field_divB0.inp**:
```fortran
0                  ipert        ! 0=eq only, 1=vac, 2,3=vac+plas
1                  iequil       ! 0=pert. alone, 1=with equil.
1.00               ampl         ! amplitude of perturbation, a.u.
72                 ntor         ! number of toroidal harmonics
0.99               cutoff       ! inner cutoff in psi/psi_a units
4                  icftype      ! type of coil file
'test.geqdsk'      gfile        ! equilibrium file
'coils.dat'        pfile        ! coil field file
'convexwall.dat'   convexfile   ! convex file for stretchcoords
'fluxdata/'        fluxdatapath ! directory with plasma response
```

**Simpson rule for flux integration**: Reference implementation in `fix-integration` branch of itpplasma/KAMEL repo. Plan: extract to libneo utility module, then update KAMEL to use as library routine.

---

## Phase 1: libneo - Geoflux Reference Coordinate System

### 1.1 Create Geoflux Coordinate Module in libneo
**File**: `../libneo/src/coordinates/geoflux_coordinates.f90`

**Purpose**: Provide coordinate transformations between geoflux (s_tor, θ_geo, φ) and cylindrical (R, φ, Z), similar to `vmec_coordinates.f90`.

**Implementation details**:

```fortran
module geoflux_coordinates
  use iso_c_binding
  implicit none

  ! Geoflux coordinate data structure
  type :: geoflux_t
    ! GEQDSK data
    type(geqdsk_t) :: geqdsk

    ! Magnetic axis location
    real(8) :: R_axis, Z_axis

    ! Flux conversion data
    real(8) :: psi_pol_axis, psi_pol_edge
    real(8), allocatable :: s_tor_of_psi_pol(:)  ! 1D spline for s_tor(psi_pol)
    real(8), allocatable :: psi_pol_of_s_tor(:)  ! 1D spline for inverse

    ! Flux surface geometry cache (optional, for performance)
    integer :: ns_cache, ntheta_cache
    real(8), allocatable :: R_cache(:,:)  ! R(s_tor, theta_geo)
    real(8), allocatable :: Z_cache(:,:)  ! Z(s_tor, theta_geo)
    logical :: cache_built
  end type geoflux_t

contains

  ! Initialize from GEQDSK file
  subroutine init_geoflux(self, geqdsk_file)

  ! Compute toroidal flux from poloidal flux using q-profile
  ! Integral: s_tor = integral_0^psi_pol q(psi') dpsi' / (2*pi*psi_tor_edge)
  subroutine compute_toroidal_flux_mapping(self)

  ! Build flux surface geometry cache on (s_tor, theta_geo) grid
  subroutine build_flux_surface_cache(self, ns, ntheta)

  ! Core transformations (similar to vmec_coordinates)
  subroutine geoflux_to_cyl(self, s_tor, theta_geo, phi, R, Z, &
                            dRds, dRdtheta, dZds, dZdtheta)

  subroutine cyl_to_geoflux(self, R, phi, Z, s_tor, theta_geo, &
                            dsdr, dsdz, dthetadr, dthetadz)

  ! Get field components in geoflux coordinates
  subroutine get_field_geoflux(self, s_tor, theta_geo, phi, &
                               B_s, B_theta, B_phi, Bmod)

end module geoflux_coordinates
```

**Key algorithms**:

1. **Toroidal flux computation** (using q-profile from GEQDSK):
   ```
   Note: read_eqfile1 already reads qpsi(i) array (line 418 in field_divB0.f90)!

   For each psi_pol from axis to edge:
     s_tor(psi_pol) = (1/psi_tor_edge) * integral[q(psi') * dpsi', psi'=0..psi_pol]
   where psi_tor_edge can be computed from:
     psi_tor_edge = integral[q(psi) * dpsi, psi=axis..edge]

   Use trapezoidal rule or Simpson's rule on qpsi array from GEQDSK.

   Reference implementation: Simpson rule in branch `fix-integration` of
   itpplasma/KAMEL repo (private, github.com/itpplasma/KAMEL).
   Plan: Extract to utility module in libneo, then update KAMEL to use as library routine.
   ```

2. **Flux surface tracing** (for given s_tor, theta_geo → R, Z):
   - Convert s_tor → psi_pol using spline
   - Starting from geometric ray: R₀ = R_axis + r₀*cos(θ_geo), Z₀ = Z_axis + r₀*sin(θ_geo)
   - Newton iteration to find (R, Z) where:
     - ψ_interp(R, Z) = psi_pol (on flux surface)
     - atan2(Z - Z_axis, R - R_axis) = theta_geo (correct angle)
   - Cache results on grid for performance

3. **Inverse transformation** (R, Z → s_tor, theta_geo):
   - Interpolate ψ_pol(R, Z) from GEQDSK grid
   - Compute s_tor from psi_pol using spline
   - Compute theta_geo = atan2(Z - Z_axis, R - R_axis)

**Dependencies**:
- `use geqdsk_tools, only: geqdsk_t, geqdsk_read, geqdsk_standardise`
- `use interpolate, only: ...` (for spline construction/evaluation)
- `use binsrc_sub, only: ...` (for table lookup)

**Files to reference**:
- `../libneo/src/coordinates/vmec_coordinates.f90` (template structure)
- `../libneo/src/magfie/geqdsk_tools.f90` (GEQDSK I/O)
- `../libneo/src/magfie/field_divB0.f90` (field evaluation from GEQDSK)

---

### 1.2 Create Geoflux Spline Data Module
**File**: `../libneo/src/magfie/geoflux_field.f90`

**Purpose**: Provide spline-based field evaluation in geoflux coordinates, analogous to `splint_vmec_data` for VMEC.

**Key insight**: **Reuse existing `field_divB0` module!** It already has:
- `read_eqfile1`: Reads GEQDSK format (lines 377-443 in `field_divB0.f90`)
- `field_eq`: Evaluates field in cylindrical (R,φ,Z) with derivatives (lines 127-274)
- 2D bicubic spline interpolation on ψ(R,Z)
- F_pol(ψ) handling for toroidal field

```fortran
module geoflux_field
  use geoflux_coordinates
  use field_sub, only: field_eq  ! Reuse existing field evaluation!
  implicit none

  ! Global geoflux data instance (similar to VMEC pattern)
  type(geoflux_t), save :: geoflux_data

contains

  ! Initialize geoflux data (analogous to spline_vmec_data)
  subroutine spline_geoflux_data(geqdsk_file, ns, ntheta)
    character(len=*), intent(in) :: geqdsk_file
    integer, intent(in) :: ns, ntheta

    ! Read GEQDSK using existing field_divB0 infrastructure
    ! This initializes the module-level splines automatically
    call read_eqfile1(...)  ! from field_divB0

    ! Find magnetic axis from GEQDSK data
    call find_magnetic_axis(geoflux_data)

    ! Compute s_tor from psi_pol using q-profile
    call compute_toroidal_flux_mapping(geoflux_data)

    ! Build flux surface cache
    call build_flux_surface_cache(geoflux_data, ns, ntheta)
  end subroutine

  ! Evaluate field in geoflux coordinates (analogous to splint_can_sub)
  subroutine splint_geoflux_field(s_tor, theta_geo, phi, &
                                  Acov, hcov, Bmod, sqgBctr)
    real(8), intent(in) :: s_tor, theta_geo, phi
    real(8), intent(out) :: Acov(3), hcov(3), Bmod, sqgBctr
    real(8) :: R, Z, B_R, B_phi, B_Z
    real(8) :: dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ

    ! Transform to cylindrical using flux surface cache
    call geoflux_to_cyl(geoflux_data, s_tor, theta_geo, phi, R, Z, &
                        dRds, dRdtheta, dZds, dZdtheta)

    ! Call existing field_eq from field_divB0 module
    call field_eq(R, phi, Z, B_R, B_phi, B_Z, &
                  dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)

    ! Transform B from cylindrical to geoflux covariant components
    ! Using Jacobian ∂(R,Z)/∂(s_tor,θ_geo)
    ! Compute vector potential A and normalized field h
    ! Compute Jacobian sqrt(g) and metric tensor

  end subroutine

end module geoflux_field
```

**Key implementation notes**:
- **Reuse `field_divB0.f90`** for all GEQDSK I/O and field evaluation in cylindrical
- Only add coordinate transformation layer: (s_tor, θ_geo, φ) ↔ (R, φ, Z)
- `field_eq` already handles: 2D splines, F_pol(ψ), derivatives
- Toroidal symmetry ⟹ no φ-dependence (simplifies to 2D problem)
- Compute Jacobian: √g = R * |∂(R,Z)/∂(s_tor,θ_geo)|

**3D field superposition support**:
- The `field_divB0` module supports adding 3D perturbations on top of the axisymmetric equilibrium
- 3D fields (e.g., RMPs, error fields, test coils) affect the magnetic field evaluation
- Reference coordinates (s_tor, θ_geo, φ) remain based on the axisymmetric equilibrium
- This allows studying non-axisymmetric effects without losing the clean flux coordinate system

**How 3D perturbations work in field_divB0**:

The `ipert` switch controls perturbation mode:
- `ipert=0`: Axisymmetric equilibrium only (what we need for basic geoflux implementation)
- `ipert=1`: Vacuum perturbation (cylindrical coil field via Biot-Savart)
- `ipert=2`: Vacuum + plasma response (no derivatives, uses Fourier representation)
- `ipert=3`: Vacuum + plasma response with full derivatives (7× slower)

The `iequil` switch controls whether equilibrium is included:
- `iequil=0`: Perturbation field alone (useful for debugging)
- `iequil=1`: Total field = equilibrium + perturbation (normal mode)

**Vacuum perturbation workflow** (`ipert=1`):
1. Read 3D coil field from `pfile` on cylindrical (R, φ, Z) grid
2. Supported file formats (`icftype`):
   - Type 1-3: Legacy formats with fixed grid sizes
   - Type 4: Simple format with header: `nr np nz`, `rmin rmax`, `pmin pmax`, `zmin zmax`, then `Br Bp Bz` values
3. Convert field to divergence-free representation via vector potentials (`vector_potentials`)
   - Decomposes into Fourier harmonics (up to `ntor` modes)
   - Uses stretch coordinates to handle complex geometry
4. Evaluate perturbation field at any point via `field_divfree`
5. Add to equilibrium field with amplitude scaling: `B_total = B_eq + ampl * B_pert`

**Plasma response workflow** (`ipert=2,3`):
- Reads pre-computed plasma response in flux coordinates (`fluxdatapath`)
- Uses `field_fourier` / `field_fourier_derivs` for evaluation
- Requires `inthecore` cutoff to define region where plasma response is valid

**Generating vacuum field files**:
- Use `vacfield.x` program from libneo to compute coil fields
- Reads coil geometries in various formats:
  - AUG format (ASDEX Upgrade convention)
  - GPEC format (General Perturbation Equilibrium Code)
  - Nemov format (Wendelstein 7-X convention)
  - STELLOPT/MAKEGRID filament format (coils.c09r00)
- Computes Biot-Savart field on specified (R, φ, Z) grid
- Can output either real-space field or Fourier representation

**Coil geometry tools**:
- `coil_tools.f90` module provides coil I/O and Biot-Savart routines
- `coil_convert.x` program converts between coil file formats
- Python interface via `libneo.coils` module

**Current status in libneo tests**:
- All `field_divB0.inp` test files have `ipert=0` (equilibrium only)
- No automated tests for 3D field superposition exist yet
- The infrastructure is implemented and used in production codes, but not regression-tested

**For SIMPLE+GEQDSK implementation**:
1. Phase 1 (basic): Use `ipert=0` for pure axisymmetric equilibrium
2. Phase 2 (advanced): Enable `ipert=1` for RMP/error field studies
   - Requires pre-computing coil fields on appropriate grid
   - Use `vacfield.x` to generate field files from coil geometries
3. Phase 3 (full): Add `ipert=2,3` for self-consistent plasma response
   - Requires coupling to GPEC, MARS-F, or similar MHD codes
   - Out of scope for initial implementation

**Trade-offs on GEQDSK infrastructure**:

Option A: **Use `field_divB0` entirely** (RECOMMENDED)
- ✅ Battle-tested field evaluation with derivatives
- ✅ 2D bicubic splines already implemented
- ✅ Handles F_pol(ψ) correctly
- ⚠️ No COCOS standardization (assumes specific convention)
- ⚠️ Module-level variables (but threadprivate for OpenMP)

Option B: Use `geqdsk_tools` for I/O + implement field evaluation
- ✅ Modern interface with COCOS awareness
- ✅ Clean data structures (geqdsk_t type)
- ❌ Need to implement field evaluation from scratch
- ❌ More work, potential for bugs

**Decision**: Use `field_divB0` for field evaluation. Optionally wrap reading in `geqdsk_tools` for COCOS standardization, then pass data to `field_divB0` structures. Can refactor later if needed.

---

### 1.3 Update libneo Build System
**Files to modify**:
- `../libneo/CMakeLists.txt`: Add new source files
- `../libneo/fpm.toml`: Add to source list (if using fpm)

**Changes**:
```cmake
# Add to libneo sources
set(LIBNEO_SOURCES
  ...
  src/coordinates/geoflux_coordinates.f90
  src/magfie/geoflux_field.f90
  ...
)
```

---

### 1.4 Create libneo Tests for Geoflux
**File**: `../libneo/test/source/test_geoflux.f90`

**Test cases**:
1. **Load GEQDSK**: Read and standardize GEQDSK file
2. **Toroidal flux mapping**: Verify s_tor computation from q-profile
3. **Coordinate transformations**:
   - Round-trip: (s_tor, θ, φ) → (R, φ, Z) → (s_tor, θ, φ)
   - Check Jacobians and metric tensors
4. **Flux surface accuracy**: Verify ψ(R(s,θ), Z(s,θ)) = ψ(s) to tolerance
5. **Field evaluation**: Check ∇·B = 0, compare with known tokamak solutions

**Test data**: Use `../libneo/python/tests/test.geqdsk` (MAST equilibrium)

---

## Phase 2: SIMPLE - GEQDSK Field Class (Boozer-like representation)

### 2.1 Create GEQDSK Field Type
**File**: `src/field/field_geqdsk.f90`

**Purpose**: Extend `MagneticField` base class to handle GEQDSK files with geoflux reference coordinates.

```fortran
module field_geqdsk
  use field_base, only: MagneticField
  use geoflux_field, only: spline_geoflux_data, splint_geoflux_field
  implicit none

  type, extends(MagneticField) :: GeqdskField
    character(len=256) :: filename
  contains
    procedure :: init => init_geqdsk
    procedure :: evaluate => evaluate_geqdsk
  end type GeqdskField

contains

  subroutine init_geqdsk(self, filename)
    class(GeqdskField), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer :: ns, ntheta

    self%filename = filename

    ! Use reasonable grid resolution (adjust as needed)
    ns = 128
    ntheta = 256

    ! Initialize libneo geoflux data
    call spline_geoflux_data(filename, ns, ntheta)
  end subroutine

  subroutine evaluate_geqdsk(self, x, Acov, hcov, Bmod, sqgBctr)
    class(GeqdskField), intent(in) :: self
    real(8), intent(in) :: x(3)  ! (r=√s_tor, theta_geo, phi)
    real(8), intent(out) :: Acov(3), hcov(3), Bmod, sqgBctr
    real(8) :: s_tor

    ! Convert r to s_tor (SIMPLE convention: x(1) = √s)
    s_tor = x(1)**2

    ! Call libneo evaluation
    call splint_geoflux_field(s_tor, x(2), x(3), Acov, hcov, Bmod, sqgBctr)
  end subroutine

end module field_geqdsk
```

**Note**: This provides geoflux coordinates as reference system, similar to how `field_vmec.f90` provides VMEC coordinates as reference.

---

### 2.2 Update Field Dispatcher
**File**: `src/field.F90`

**Modify**: `field_from_file` subroutine (around line 15-44)

**Add**:
```fortran
use field_geqdsk, only: GeqdskField

! In field_from_file routine:
else if (endswith(filename, '.geqdsk') .or. endswith(filename, '.eqdsk')) then
  allocate(GeqdskField :: field)
  call field%init(filename)
```

---

### 2.3 Update Coordinate Module
**File**: `src/coordinates/coordinates.f90`

**Add**: Geoflux ↔ Cylindrical ↔ Cartesian transformations

```fortran
! Add to get_transform function
case('geoflux')
  select case(trim(to))
    case('cyl')
      get_transform => geoflux_to_cyl_wrapper
    case('cart')
      get_transform => geoflux_to_cart_wrapper
  end select
case('cyl')
  if (trim(to) == 'geoflux') then
    get_transform => cyl_to_geoflux_wrapper
  end if
```

**Implementation**: Wrappers call libneo's `geoflux_coordinates` module.

---

## Phase 3: SIMPLE - Canonical Coordinates on Geoflux

### 3.1 Meiss Canonical Coordinates for Geoflux
**File**: `src/field/field_can_meiss_geoflux.f90`

**Purpose**: Implement symmetry-flux canonical coordinates using geoflux as reference (instead of VMEC).

**Strategy**:
1. Copy `src/field/field_can_meiss.f90` as template
2. Replace all VMEC-specific calls with geoflux equivalents:
   - `vmec_to_cyl` → `geoflux_to_cyl`
   - VMEC spline evaluations → geoflux spline evaluations
3. Simplify for axisymmetry:
   - No φ-grid needed (∂/∂φ = 0 for equilibrium)
   - Can use 2D batch splines instead of 3D (save memory)
   - ODE integration only in (s_tor, θ_geo) plane

**Key changes from VMEC-based Meiss**:
```fortran
module field_can_meiss_geoflux
  use field_base
  use geoflux_field
  use interpolate  ! For batch splines

  type, extends(MagneticField) :: MeissGeofluxField
    ! Reference field (geoflux)
    character(len=256) :: geqdsk_file

    ! Canonical transformation batch splines (2D in s_tor, theta_geo)
    ! Axisymmetry ⟹ no phi dependence for equilibrium quantities
    type(batch_spline_2d_t) :: spline_transformation  ! [lambda_phi, chi_gauge]
    type(batch_spline_2d_t) :: spline_field           ! [A_*, h_*, Bmod, sqgBctr]
  contains
    procedure :: init => init_meiss_geoflux
    procedure :: evaluate => evaluate_meiss_geoflux
  end type

contains

  subroutine init_meiss_geoflux(self, geqdsk_file, ns, ntheta)
    ! 1. Initialize geoflux reference field
    call spline_geoflux_data(geqdsk_file, ns, ntheta)

    ! 2. Integrate ODEs for canonical transformation (simplified for axisymmetry)
    ! Grid: 2D (s_tor, theta_geo), no phi needed
    ! ODE: dlambda_phi/ds = -h_s / h_phi, dchi/ds = A_s + A_phi * dlambda/ds

    ! 3. Build 2D batch splines for transformation and field
  end subroutine

end module
```

**Testing**: Verify orbits match between Meiss-geoflux and Meiss-VMEC for same tokamak equilibrium (if VMEC file available).

---

### 3.2 Albert Canonical Coordinates for Geoflux
**File**: `src/field/field_can_albert_geoflux.f90`

**Purpose**: Implement poloidal-flux canonical coordinates on top of Meiss-geoflux (replacing s_tor with ψ_pol).

**Strategy**:
1. Copy `src/field/field_can_albert.f90` as template
2. Replace VMEC references with geoflux
3. Use existing `psi_transform` module to switch from s_tor to ψ_pol
4. This is the **primary coordinate system for GEQDSK**: matches natural poloidal-flux representation

**Key changes**:
```fortran
module field_can_albert_geoflux
  use field_can_meiss_geoflux
  use psi_transform

  type, extends(MeissGeofluxField) :: AlbertGeofluxField
    ! Add PSI transformation on top of Meiss
    type(psi_grid_t) :: psi_grid
  contains
    procedure :: init => init_albert_geoflux
  end type

contains

  subroutine init_albert_geoflux(self, geqdsk_file, ns, ntheta)
    ! 1. Initialize Meiss-geoflux base
    call self%MeissGeofluxField%init(geqdsk_file, ns, ntheta)

    ! 2. Apply PSI transformation: s_tor → psi_pol
    ! This makes radial coordinate proportional to poloidal flux
    call init_psi_transform(self%psi_grid, ...)

    ! 3. Rebuild splines with new radial coordinate
    call self%rebuild_splines_with_psi()
  end subroutine

end module
```

**Note**: This matches the "geoflux" name - Albert coordinates naturally use geometric flux (ψ_pol) as radial label.

---

### 3.3 Update Field Initialization in Main Code
**File**: `src/simple.f90`

**Modify**: Initialization section (around lines 47-75)

**Add**:
```fortran
! Detect field type from filename extension
if (endswith(equil_file, '.nc')) then
  ! VMEC file - existing logic
  call spline_vmec_data(...)
  ! Initialize VMEC-based canonical fields

else if (endswith(equil_file, '.geqdsk') .or. endswith(equil_file, '.eqdsk')) then
  ! GEQDSK file - new logic
  select case(trim(field_type))
    case('geoflux', 'boozer')
      ! Direct geoflux reference coordinates
      allocate(GeqdskField :: mag_field)

    case('meiss')
      ! Meiss canonical on geoflux reference
      allocate(MeissGeofluxField :: mag_field)

    case('albert')
      ! Albert canonical on geoflux reference
      allocate(AlbertGeofluxField :: mag_field)

    case default
      print *, 'For GEQDSK files, field_type must be: geoflux, meiss, or albert'
      stop
  end select

else
  print *, 'Unknown equilibrium file format'
  stop
end if
```

---

## Phase 4: Testing and Validation

### 4.1 Unit Tests

**Files to create**:

1. **libneo tests** (`../libneo/test/source/test_geoflux.f90`):
   - Load MAST GEQDSK: `../libneo/python/tests/test.geqdsk`
   - Test coordinate transformations
   - Verify field evaluation (∇·B = 0)
   - Check flux surface mapping accuracy

2. **SIMPLE integration test** (`test/fortran/test_geqdsk_field.f90`):
   - Load GEQDSK field
   - Initialize particle on flux surface
   - Trace single orbit
   - Verify energy conservation

3. **Canonical coordinate test** (`test/fortran/test_can_geoflux.f90`):
   - Compare Meiss-geoflux vs Meiss-VMEC (if VMEC available for same tokamak)
   - Verify Hamiltonian conservation
   - Check symplecticity of transformations

---

### 4.2 Integration Tests

**Test cases** (use `examples/` directory):

1. **Simple tokamak orbit**:
   - Input: `examples/simple_geqdsk.in` (new file)
   - GEQDSK: Use MAST equilibrium from libneo
   - Field type: `geoflux` (reference coordinates)
   - Trace 100 particles, verify confinement statistics

2. **Meiss canonical tokamak**:
   - Field type: `meiss`
   - Compare with `geoflux` results (should be equivalent modulo gauge)

3. **Albert canonical tokamak**:
   - Field type: `albert`
   - Test with various tokamak profiles (circular, shaped, etc.)

---

### 4.3 Comparison with Existing Codes

**Validation strategy**:
1. Generate GEQDSK from known VMEC tokamak equilibrium (if converter available)
2. Trace identical particle sets with both:
   - SIMPLE + VMEC file (existing)
   - SIMPLE + GEQDSK file (new)
3. Compare:
   - Orbit trajectories (should match to tolerance)
   - Confinement times
   - Poincaré sections

**Alternative**: Compare with established tokamak orbit codes (e.g., ORBIT, VENUS, etc.) using same GEQDSK.

---

### 4.4 Test Data Organization

**Location**: `test/data/geqdsk/`

**Files to include**:
- `test.geqdsk` (symlink to `../libneo/python/tests/test.geqdsk`)
- Additional GEQDSK files from literature or databases:
  - ITER equilibrium
  - NSTX/MAST (spherical tokamak)
  - DIII-D (shaped tokamak)
  - Circular tokamak (for analytical validation)

---

## Phase 5: Documentation and Examples

### 5.1 Update Documentation

**Files to modify**:

1. **README.md**: Add GEQDSK support to features list
2. **DESIGN.md**: Document geoflux coordinate system architecture
3. **examples/README.md**: Add GEQDSK usage examples

---

### 5.2 Example Input Files

**Create**: `examples/simple_geqdsk.in`

```fortran
&simple_config
  equil_file = 'test.geqdsk'
  field_type = 'albert'  ! Options: geoflux, meiss, albert

  ! Particle initialization
  nparticles = 1000
  sample_mode = 'surface'  ! Start on flux surface
  s_sample = 0.5

  ! Integration parameters
  integrator = 'lobatto4'
  dt = 1.0d-2
  nsteps = 100000

  ! Output
  output_step = 100
  write_orbits = .true.
/
```

---

### 5.3 Usage Examples

**Create**: `examples/example_geqdsk.f90`

Minimal example demonstrating:
1. Load GEQDSK file
2. Initialize geoflux field
3. Trace single particle orbit
4. Output to file

---

## Phase 6: Performance Optimization

### 6.1 Flux Surface Cache Optimization

**Approach**: Pre-compute (R, Z) on dense (s_tor, theta_geo) grid

**Trade-offs**:
- Memory: ~O(ns × ntheta × 2) doubles (e.g., 256×512×2 = 2.6 MB, negligible)
- Speed: Replace Newton iteration with 2D spline interpolation
- Accuracy: Control via grid resolution

**Implementation**: Already included in `geoflux_t%R_cache`, `geoflux_t%Z_cache`.

---

### 6.2 Profile Spline Optimization

**Approach**: Pre-spline all 1D profiles (fpol, q, etc.) during initialization

**Benefit**: Avoid repeated spline construction during orbit integration

---

### 6.3 Axisymmetry Exploitation

**Simplifications** for canonical coordinates:
- 2D batch splines instead of 3D (factor of ~nφ memory reduction)
- Skip φ-derivatives in ODE integration
- Simplified Jacobian calculations

**Expected speedup**: ~2-3× for canonical coordinate initialization compared to 3D stellarator case.

---

## Implementation Priority and Timeline

### Minimal Working Implementation (Priority 1)
**Goal**: Trace particles in GEQDSK file with geoflux reference coordinates

**Tasks**:
- [ ] Phase 1.1: `geoflux_coordinates.f90` in libneo (coordinate transformations)
- [ ] Phase 1.2: `geoflux_field.f90` in libneo (thin wrapper around `field_divB0`)
- [ ] Phase 2.1: `field_geqdsk.f90` in SIMPLE
- [ ] Phase 2.2: Update field dispatcher
- [ ] Phase 4.1: Basic unit tests

**Estimated effort**: 1-2 weeks (reduced from 2-3 weeks due to reusing `field_divB0`)

---

### Canonical Coordinates (Priority 2)
**Goal**: Support Meiss and Albert canonical coordinates for tokamaks

**Tasks**:
- [ ] Phase 3.1: `field_can_meiss_geoflux.f90`
- [ ] Phase 3.2: `field_can_albert_geoflux.f90`
- [ ] Phase 3.3: Update `simple.f90` initialization
- [ ] Phase 4.2-4.3: Integration and validation tests

**Estimated effort**: 2-3 weeks (after Priority 1 complete)

---

### Polish and Documentation (Priority 3)
**Goal**: Production-ready feature

**Tasks**:
- [ ] Phase 4.4: Comprehensive test suite
- [ ] Phase 5: Documentation and examples
- [ ] Phase 6: Performance optimization

**Estimated effort**: 1-2 weeks

---

## Technical Challenges and Mitigation

### Challenge 1: Toroidal Flux Computation
**Issue**: Need to integrate q(ψ_pol) from GEQDSK to get ψ_tor(ψ_pol)

**Solution**:
- Use trapezoidal or Simpson's rule on q-profile
- Verify: ι = 1/q should match rotational transform
- Test with known analytical profiles (e.g., circular tokamak)

---

### Challenge 2: Magnetic Axis Location
**Issue**: Need magnetic axis location (R_axis, Z_axis) for geometric angle θ_geo

**Solution**:
- **Already provided by GEQDSK!** `read_eqfile1` returns `rmaxis, zmaxis` (line 407)
- Validate: ∇ψ should vanish at axis
- Alternative: Search for minimum on ψ grid if needed: `minloc(geqdsk%psirz)`

---

### Challenge 3: Flux Surface Singular Behavior
**Issue**: At magnetic axis, θ_geo becomes ill-defined (polar singularity)

**Solution**:
- Exclude innermost surfaces: Start grid at s_tor_min = 0.01 (not 0)
- Use L'Hôpital's rule or Taylor expansion for near-axis derivatives
- Alternative: Use Cartesian representation near axis

---

### Challenge 4: Separatrix and SOL
**Issue**: Particles may reach separatrix (ψ = ψ_edge) or beyond

**Solution**:
- Detect when s_tor > 1.0 (outside separatrix)
- Terminate orbit or extend field to SOL with simple model
- Document limitation in user guide

---

### Challenge 5: COCOS Convention Complications
**Issue**: Different GEQDSK files may use different sign conventions

**Solution**:
- **Already handled by libneo**: `geqdsk_standardise` converts to COCOS 3
- Always call standardization after reading
- Include COCOS info in log output for debugging

---

## Testing Strategy

### Unit Tests (libneo)
- [x] GEQDSK read and standardize
- [ ] Toroidal flux computation from q-profile
- [ ] Flux surface tracing accuracy
- [ ] Coordinate transformation round-trips
- [ ] Field evaluation: verify ∇·B = 0

### Integration Tests (SIMPLE)
- [ ] Load GEQDSK and initialize field
- [ ] Single particle orbit (energy conservation)
- [ ] Multiple particles (statistics)
- [ ] Compare geoflux vs Meiss vs Albert

### Validation Tests
- [ ] Compare with VMEC-based results (if available)
- [ ] Compare with other tokamak codes
- [ ] Analytical test cases (circular tokamak)

---

## Dependencies and Prerequisites

### External Libraries (already available)
- ✓ LAPACK/BLAS (linear algebra)
- ✓ NetCDF (for GEQDSK? No, GEQDSK is ASCII)
- ✓ libneo (GEQDSK tools, splines, field evaluation)

### Internal Modules (reuse from VMEC implementation)
- ✓ `interpolate` module (batch splines)
- ✓ `psi_transform` module (Albert coordinates)
- ✓ `binsrc_sub` (binary search)
- ✓ Integration infrastructure (ODE solvers, symplectic integrators)

### New Modules Required
- [ ] `geoflux_coordinates` (Phase 1.1)
- [ ] `geoflux_field` (Phase 1.2)
- [ ] `field_geqdsk` (Phase 2.1)
- [ ] `field_can_meiss_geoflux` (Phase 3.1)
- [ ] `field_can_albert_geoflux` (Phase 3.2)

---

## File Structure Summary

### libneo (new files)
```
../libneo/
├── src/
│   ├── coordinates/
│   │   └── geoflux_coordinates.f90       [NEW - Phase 1.1]
│   └── magfie/
│       └── geoflux_field.f90              [NEW - Phase 1.2]
└── test/
    └── source/
        └── test_geoflux.f90               [NEW - Phase 4.1]
```

### SIMPLE (new files)
```
SIMPLE/
├── src/
│   ├── field/
│   │   ├── field_geqdsk.f90               [NEW - Phase 2.1]
│   │   ├── field_can_meiss_geoflux.f90    [NEW - Phase 3.1]
│   │   └── field_can_albert_geoflux.f90   [NEW - Phase 3.2]
│   ├── field.F90                          [MODIFY - Phase 2.2]
│   ├── coordinates/
│   │   └── coordinates.f90                [MODIFY - Phase 2.3]
│   └── simple.f90                         [MODIFY - Phase 3.3]
├── test/
│   ├── fortran/
│   │   ├── test_geqdsk_field.f90          [NEW - Phase 4.2]
│   │   └── test_can_geoflux.f90           [NEW - Phase 4.2]
│   └── data/
│       └── geqdsk/
│           └── test.geqdsk                [SYMLINK to libneo]
└── examples/
    ├── simple_geqdsk.in                   [NEW - Phase 5.2]
    └── example_geqdsk.f90                 [NEW - Phase 5.3]
```

---

## Key Design Decisions

### 1. Why implement in libneo?
**Rationale**: Follow VMEC pattern where coordinate transformations live in libneo, allowing reuse across multiple codes.

### 2. Why "geoflux" name?
**Rationale**:
- Emphasizes **geometric** poloidal angle (not Boozer/Hamada)
- Combines with toroidal **flux** as radial coordinate
- Distinguishes from VMEC's Fourier-based coordinates

### 3. Why support both Meiss and Albert?
**Rationale**:
- **Meiss**: Uses toroidal flux s_tor (natural for energy conservation)
- **Albert**: Uses poloidal flux ψ_pol (natural for GEQDSK representation)
- Gives users flexibility depending on application

### 4. Why cache flux surfaces?
**Rationale**:
- Newton iteration for (s_tor, θ_geo) → (R, Z) is expensive (~10-20 iterations)
- Pre-computing on grid reduces to 2D spline evaluation (~10× faster)
- Memory cost is negligible (<10 MB for typical grids)

---

## Success Criteria

### Milestone 1: Basic Functionality
- [ ] Load GEQDSK file and trace particle orbits
- [ ] Output matches expected tokamak behavior (banana orbits, passing orbits)
- [ ] Passes unit tests for coordinate transformations

### Milestone 2: Canonical Coordinates
- [ ] Meiss and Albert coordinates implemented
- [ ] Hamiltonian conserved to machine precision
- [ ] Symplectic structure preserved (Poincaré return map)

### Milestone 3: Production Ready
- [ ] Full test suite passing (unit + integration + validation)
- [ ] Documentation complete
- [ ] Performance comparable to VMEC-based workflow

---

## Future Extensions

### Included via field_divB0:
- **Non-axisymmetric perturbations**: 3D fields (RMPs, TBMs, error fields) on top of GEQDSK
  - Supported through `field_divB0` module's superposition capability
  - Reference coordinates remain axisymmetric (based on equilibrium)
  - Perturbations affect field evaluation but not coordinate system

### Out of Scope:
- **Kinetic MHD equilibria**: Pressure anisotropy, flows
- **Real-time equilibrium reconstruction**: EFIT/LIUQE interface
- **Scrape-off layer**: Extend field beyond separatrix

These could be added later as separate features building on the geoflux foundation.

---

## References

### Code References
- SIMPLE field system: `src/field/field_base.f90`, `src/field/field_vmec.f90`
- SIMPLE canonical coordinates: `src/field/field_can_*.f90`
- libneo VMEC: `../libneo/src/coordinates/vmec_coordinates.f90`
- libneo GEQDSK: `../libneo/src/magfie/geqdsk_tools.f90`

### Scientific References
- Meiss coordinates: [cite canonical coordinate papers]
- Albert coordinates: [cite geoflux canonical coordinate papers]
- GEQDSK format: EFIT documentation
- COCOS conventions: O. Sauter & S.Yu. Medvedev, Comput. Phys. Commun. 184, 293 (2013)

---

## Notes

- This TODO assumes familiarity with the existing SIMPLE and libneo architecture
- Adjust grid resolutions (ns, ntheta) based on performance testing
- Consider parallelizing flux surface tracing if it becomes a bottleneck
- Keep GEQDSK-specific code isolated for maintainability