# GVEC Magnetic Field Analysis

Analysis of the GVEC code in `thirdparty/gvec` for reading `.dat` files and constructing magnetic fields via splines.

## Overview

GVEC (Grad-Variational Equilibrium Code) provides comprehensive infrastructure for reading magnetic field data from various formats and constructing high-quality spline representations for field evaluation. The system supports both ASCII `.dat` files and binary/NetCDF formats.

## .dat File Format and Reading

### File Format Structure

GVEC uses ASCII `.dat` files to store MHD3D solution data. The format includes:

```
## MHD3D Solution... outputLevel and fileID:
0000,00002479
## grid: nElems, gridType
       2,       0
## grid: sp(0:nElems)  
  0.000000000000000E+00,  0.500000000000000E+00,  0.100000000000000E+01
## global: nfp,degGP,mn_nyq(2),hmap
       1,       7,      13,       1,       1
## X1_base: s%nbase,s%deg,s%continuity,f%modes,f%sin_cos,f%excl_mn_zero
       7,       5,       4,       4,       2,       0
## X2_base: s%nbase,s%deg,s%continuity,f%modes,f%sin_cos,f%excl_mn_zero
       7,       5,       4,       3,       1,       0
## LA_base: s%nbase,s%deg,s%continuity,f%modes,f%sin_cos,f%excl_mn_zero
       7,       5,       4,       3,       1,       1
## X1: m,n,X1(1:nbase,iMode)
       0,       0,  0.505481389082735E+01,  0.505481389082735E+01, ...
## X2: m,n,X2(1:nbase,iMode)
...
## LA: m,n,LA(1:nbase,iMode)
...
## profiles at X1_base IP points : spos,phi,chi,iota,pressure
...
## a_minor,r_major,volume
...
```

The file contains:
- **Grid information**: Elements, grid type, and flux surface coordinates (`sp`)
- **Basis specifications**: Spline degree, continuity, Fourier modes for X1, X2, LA components
- **Fourier coefficients**: Mode numbers (m,n) and coefficients for each magnetic field component
- **Profiles**: Radial profiles of φ (toroidal flux), χ (poloidal flux), ι (rotational transform), pressure
- **Geometry**: Minor radius, major radius, plasma volume

### Key Reading Function

**Location**: `/Users/ert/code/SIMPLE/thirdparty/gvec/src/readstate/readstate.f90`

**Function**: `ReadStateFileFromASCII(fileString, hmap_in)`

```fortran
SUBROUTINE ReadStateFileFromASCII(fileString,hmap_in)
```

**Key operations**:
1. **File validation**: Checks file existence
2. **Header reading**: Reads grid parameters, basis configurations, global parameters
3. **Fourier mode data**: Reads mode numbers and coefficients for X1, X2, LA components
4. **Profile data**: Reads radial profiles (φ, χ, ι, pressure) at interpolation points
5. **Spline construction**: Converts profile data to B-spline DOF (Degrees of Freedom)

**Key code sections**:
```fortran
! Read Fourier coefficients for each component
DO iMode=1,X1_modes_r
  read(ioUnit,*)X1_mn_r(:,iMode),X1_r(:,iMode)
END DO

! Read radial profiles at interpolation points
DO is=1,X1_nbase_r
  read(ioUnit,*)profiles_IP(is,:)  ! spos,phi,chi,iota,pressure
END DO

! Convert to spline DOF
profiles_1d(:,1) = X1_base_r%s%initDOF( profiles_IP(:,1+1) ) !phi
profiles_1d(:,3) = X1_base_r%s%initDOF( profiles_IP(:,1+3) ) !iota  
profiles_1d(:,4) = X1_base_r%s%initDOF( profiles_IP(:,1+4) ) !pressure
```

## Spline Construction System

### Core Spline Types

#### 1. Cubic Spline (`t_cubspl`)
**Location**: `/Users/ert/code/SIMPLE/thirdparty/gvec/src/base/cubic_spline.f90`

```fortran
TYPE :: t_cubspl
  REAL(wp),ALLOCATABLE  :: coefs(:)   ! B-Spline coefficients
  REAL(wp),ALLOCATABLE  :: knots(:)   ! knots (break points)
  CLASS(sll_c_bsplines),ALLOCATABLE :: bspl ! b-spline class
END TYPE t_cubspl
```

**Constructor**: `cubspl_new(x, f, BC, BC_val)`
- **x**: Position coordinates
- **f**: Function values at positions
- **BC**: Boundary conditions (0=not-a-knot, 1=first derivative, 2=second derivative)
- **BC_val**: Boundary values for BC types 1,2

**Key boundary condition handling**:
```fortran
SELECT CASE(BC(1))
CASE(0) !not-a-knot: leave out second knot
  nbreaks=nbreaks-1
  innerpts(1)=3
CASE(1,2) !derivative boundary conditions
  ncoefs=ncoefs+1
  nbc_left=1
END SELECT
```

#### 2. B-spline Profiles (`t_rProfile_bspl`)
**Location**: `/Users/ert/code/SIMPLE/thirdparty/gvec/src/profiles/rprofile_bspline.f90`

```fortran
TYPE, EXTENDS(c_rProfile) :: t_rProfile_bspl
  INTEGER               :: n_knots, deg
  REAL(wp), ALLOCATABLE :: knots(:)   ! knot values
  REAL(wp), ALLOCATABLE :: coefs(:)   ! B-Spline coefficients  
  CLASS(sll_c_bsplines),ALLOCATABLE :: bspl
END TYPE t_rProfile_bspl
```

**Key methods**:
- `eval_at_rho2(rho2, deriv)`: Evaluates profile at flux coordinate ρ²
- `antiderivative()`: Computes exact antiderivatives for integration

### B-spline Foundation Library

**Location**: `/Users/ert/code/SIMPLE/thirdparty/gvec/src/base/bsplines/`

Key modules:
- **`sll_m_bsplines.f90`**: Core B-spline basis functions
- **`sll_m_spline_1d.f90`**: 1D spline evaluation and coefficient management
- **`sll_m_spline_interpolator_1d.f90`**: 1D spline interpolation algorithms
- **`sll_m_spline_matrix.f90`**: Matrix operations for spline coefficient solving

## Magnetic Field Evaluation Workflow

### 1. Data Loading
```fortran
CALL ReadStateFileFromASCII(filename, hmap)
```
- Reads `.dat` file containing Fourier coefficients and profile data
- Validates grid and basis specifications

### 2. Spline Construction
```fortran
! Convert profiles to spline DOF
profiles_1d(:,1) = X1_base_r%s%initDOF( profiles_IP(:,1+1) ) !phi
profiles_1d(:,3) = X1_base_r%s%initDOF( profiles_IP(:,1+3) ) !iota
profiles_1d(:,4) = X1_base_r%s%initDOF( profiles_IP(:,1+4) ) !pressure
```
- Constructs B-spline representations of radial profiles
- Sets up basis functions for field components (X1, X2, LA)

### 3. Field Evaluation
```fortran
! Evaluate profiles at flux coordinate
phi_val = eval_phi_r(rho2)
iota_val = eval_iota_r(rho2) 
pres_val = eval_pres_r(rho2)
```
- Uses spline evaluation for smooth field interpolation
- Supports derivative computation for field gradients

### 4. Magnetic Field Components
The system handles three main magnetic field components:
- **X1**: Related to radial field component
- **X2**: Related to poloidal field component  
- **LA**: Lambda (magnetic stream function)

Each component is represented as Fourier series with spline-interpolated coefficients.

## VMEC Integration

**Location**: `/Users/ert/code/SIMPLE/thirdparty/gvec/src/vmec/vmec_readin.f90`

GVEC also supports reading VMEC equilibrium files:

```fortran
SUBROUTINE ReadVMEC(fileName, file_Format)
```

Key VMEC data structures:
- **Fourier coefficients**: `rmnc`, `zmns`, `lmns` (and asymmetric components if `lasym=.TRUE.`)
- **Profiles**: `iotaf` (iota), `presf` (pressure), `phi` (toroidal flux)
- **Grid**: `nFluxVMEC` flux surfaces, `mn_mode` Fourier modes

## Implementation Notes for SIMPLE Integration

### Recommended Integration Approach

1. **File Reading**: Use `ReadStateFileFromASCII()` pattern to parse `.dat` files
2. **Spline Setup**: Create `t_cubspl` objects for radial profiles (φ, ι, pressure)
3. **Field Evaluation**: Implement field evaluation using spline interpolation
4. **Fourier Reconstruction**: Use Fourier coefficients to reconstruct field components

### Key Data Structures to Implement

```fortran
type, extends(MagneticField) :: GvecField
    ! Grid information
    integer :: nElems, nfp
    real(dp), allocatable :: sp(:)  ! flux surface coordinates
    
    ! Fourier mode data
    integer :: X1_modes, X2_modes, LA_modes
    integer, allocatable :: X1_mn(:,:), X2_mn(:,:), LA_mn(:,:)
    real(dp), allocatable :: X1_coefs(:,:), X2_coefs(:,:), LA_coefs(:,:)
    
    ! Profile splines
    type(cubic_spline) :: phi_spline, iota_spline, pres_spline
    
    ! Geometry
    real(dp) :: a_minor, r_major, volume
end type GvecField
```

### Field Evaluation Strategy

```fortran
subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    ! x(1) = r = sqrt(ρ²), x(2) = θ, x(3) = φ
    
    ! 1. Evaluate profiles at flux coordinate
    rho2 = x(1)**2
    phi_val = self%phi_spline%eval(rho2)
    iota_val = self%iota_spline%eval(rho2)
    
    ! 2. Reconstruct field components using Fourier series
    ! Use X1, X2, LA coefficients with mode numbers
    
    ! 3. Transform to required coordinate system
    ! Compute Acov, hcov, Bmod from field components
end subroutine
```

This analysis provides the foundation for implementing GVEC `.dat` file reading and magnetic field construction in the SIMPLE field_gvec module.