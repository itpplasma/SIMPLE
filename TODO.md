# Spline Refactoring TODO

## Overview
Replace custom spline implementations with libneo's `interpolate` module to reduce code duplication and improve maintainability.

## Current State Analysis

### 1. Implementations Already Using libneo interpolate ✓
- `field_can_albert.f90` - Uses `SplineData3D` from libneo
- `field_can_meiss.f90` - Uses `SplineData3D` from libneo

### 2. Custom Implementations to Replace

#### 2.1 `splint_can_coord` in `get_canonical_coordinates.F90`
- **Location**: Lines 463-925
- **Purpose**: Interpolates canonical coordinate transformations
- **Features**:
  - 3D tensor product splines (r, theta, phi)
  - Computes up to 3rd order derivatives
  - Handles multiple quantities simultaneously
  - Uses custom `derf1`, `derf2`, `derf3` arrays for derivatives
- **Dependencies**: 
  - `canonical_coordinates_mod`
  - `vector_potential_mod`
  - `new_vmec_stuff_mod`

#### 2.2 `splint_boozer_coord` in `boozer_converter.F90`
- **Purpose**: Interpolates Boozer coordinate transformations
- **Features**:
  - Similar structure to `splint_can_coord`
  - 3D tensor product splines
  - Computes 1st and 2nd order derivatives
  - Uses same derivative arrays pattern

#### 2.3 `splint_vmec_data` in `spline_vmec_data.f90`
- **Purpose**: Interpolates VMEC equilibrium data
- **Features**:
  - 3D tensor product splines
  - Computes 1st order derivatives
  - Simpler than coordinate transformation splines

## Refactoring Plan

### Phase 1: Extend libneo interpolate Module
1. **Add derivative support** to libneo's interpolate module
   - [ ] Add derivative computation up to 3rd order
   - [ ] Support multiple derivative arrays (derf1, derf2, derf3)
   - [ ] Ensure backward compatibility with existing SplineData3D usage

2. **Create helper functions** for common operations
   - [ ] Batch evaluation for multiple quantities
   - [ ] Derivative array initialization
   - [ ] Periodic boundary handling

### Phase 2: Create Adapter Layer
1. **Design adapter interface** that mimics current calling conventions
   - [ ] Create `canonical_spline_adapter` module
   - [ ] Create `boozer_spline_adapter` module
   - [ ] Create `vmec_spline_adapter` module

2. **Implement adapters** that:
   - [ ] Convert between old data structures and SplineData3D
   - [ ] Handle derivative array passing
   - [ ] Maintain exact numerical compatibility

### Phase 3: Gradual Migration
1. **Start with simplest case**: `splint_vmec_data`
   - [ ] Create vmec_spline_adapter
   - [ ] Replace splint_vmec_data calls
   - [ ] Verify numerical results match exactly
   - [ ] Run full test suite

2. **Move to Boozer coordinates**: `splint_boozer_coord`
   - [ ] Create boozer_spline_adapter
   - [ ] Replace splint_boozer_coord calls
   - [ ] Verify field_can_boozer tests pass
   - [ ] Check conservation properties

3. **Finally canonical coordinates**: `splint_can_coord`
   - [ ] Create canonical_spline_adapter
   - [ ] Replace splint_can_coord calls
   - [ ] Verify field_can_flux tests pass
   - [ ] Validate physics results

### Phase 4: Cleanup
1. **Remove old implementations**
   - [ ] Delete splint_can_coord
   - [ ] Delete splint_boozer_coord
   - [ ] Delete splint_vmec_data
   - [ ] Remove associated helper routines

2. **Optimize and document**
   - [ ] Profile performance vs old implementation
   - [ ] Document new interfaces
   - [ ] Update examples

## Technical Challenges

### 1. Derivative Arrays
The current implementations use module-level derivative arrays (`derf1`, `derf2`, `derf3`) that would create circular dependencies. Solution: Pass these as arguments to the interpolation routines.

### 2. Numerical Compatibility
The existing code has been validated extensively. Any changes must maintain bit-for-bit compatibility or demonstrate improved accuracy.

### 3. Performance
Current implementations use hand-optimized nested loops. Need to ensure libneo's implementation is equally efficient.

## Success Criteria
- [ ] All tests pass with numerical tolerance < 1e-14
- [ ] No performance regression (< 5% slower)
- [ ] Reduced code size by at least 500 lines
- [ ] Improved code clarity and maintainability
- [ ] No changes to public APIs (drop-in replacement)

## Testing Strategy
1. Create golden reference outputs before changes
2. Compare all outputs after each phase
3. Run physics validation tests
4. Performance benchmarks at each stage

## Timeline Estimate
- Phase 1: 2-3 days (libneo extensions)
- Phase 2: 2 days (adapter design)
- Phase 3: 3-4 days (migration and testing)
- Phase 4: 1 day (cleanup)

Total: ~10 days of focused work