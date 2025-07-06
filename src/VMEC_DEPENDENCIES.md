# VMEC Dependencies in SIMPLE Code

## Direct VMEC Calls

### 1. vmec_field (line 280)
- **Location**: `rhs_cancoord` subroutine
- **Purpose**: Evaluate magnetic field and related quantities
- **Inputs**: 
  - s (radial coordinate)
  - theta (poloidal angle)
  - varphi (toroidal angle)
- **Outputs**:
  - A_theta, A_phi (vector potential components)
  - dA_theta_ds, dA_phi_ds (radial derivatives)
  - aiota (rotational transform)
  - sqg (sqrt(g) - Jacobian)
  - alam (lambda function)
  - dl_ds, dl_dt, dl_dp (lambda derivatives)
  - Bctrvr_vartheta, Bctrvr_varphi (contravariant B components)
  - Bcovar_r, Bcovar_vartheta, Bcovar_varphi (covariant B components)

### 2. splint_lambda (lines 269, 997)
- **Location**: Newton iteration in `rhs_cancoord` and `get_torques` 
- **Purpose**: Interpolate lambda function for coordinate transformation
- **Inputs**: r, theta, varphi
- **Outputs**: alam, dl_dt (lambda and its theta derivative)

## Required Field Quantities for Abstraction

To replace direct VMEC calls with abstract field interface, we need:

1. **Magnetic Field Components**
   - Covariant components: B_r, B_theta, B_phi
   - Contravariant components: B^theta, B^phi

2. **Vector Potential**
   - A_theta, A_phi
   - Radial derivatives: dA_theta/ds, dA_phi/ds

3. **Geometric Quantities**
   - sqrt(g) (Jacobian)
   - Rotational transform (iota)
   - Lambda function and derivatives (for VMEC theta transformation)

4. **Coordinate Transformation**
   - VMEC theta <-> canonical theta conversion
   - Lambda function interpolation

## boozer_converter.f90 Dependencies

### 1. vmec_field (line 139)
- **Location**: `get_boozer_coordinates` subroutine
- **Purpose**: Same as in get_canonical_coordinates.f90
- **Context**: Computing Boozer coordinate transformation

### 2. Module imports
- Uses `spline_vmec_sub` module for spline interpolation
- Uses VMEC-specific data structures from various modules

## Refactoring Strategy

1. Create abstract interface for field evaluation that includes all required quantities
2. Implement adapter for VMEC that calls existing routines
3. Replace direct calls with interface calls
4. Test extensively to ensure numerical equivalence