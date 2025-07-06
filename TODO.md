# Refactoring Plan: Abstract Field Support for Flux and Boozer Coordinates

## Goal
Refactor Flux and Boozer coordinate implementations to use the abstract `MagneticField` interface, enabling support for GVEC and other field representations beyond VMEC.

## Current State Analysis
- **Flux/Boozer**: Directly call `vmec_field` from `spline_vmec_sub`, tightly coupled to VMEC
- **Meiss/Albert**: Use abstract `MagneticField` interface, support any field type
- **Risk**: High - these are core routines used extensively in production

## Immediate Refactoring Steps (Phase 0)

### 0.1 Analysis of get_canonical_coordinates.f90
**File size**: 1141 lines (too large!)
**Key issues identified**:
1. Monolithic `get_canonical_coordinates` subroutine (lines 23-228, 205 lines!)
2. Mixed concerns: ODE integration, spline interpolation, coordinate transforms
3. Global state via threadprivate variables in `exchange_get_cancoord_mod`
4. Direct VMEC coupling in `rhs_cancoord` via `vmec_field` call
5. Complex nested interpolation in `splint_can_coord` (lines 442-904, 462 lines!)

### 0.2 Extractable Pure Functions
These can be extracted immediately with minimal risk:
1. **Stencil initialization** (lines 60-82): Pure computation, no dependencies
2. **Progress printing** (lines 284-295): Simple I/O utility
3. **Derivative array initialization** (lines 432-436): Pure computation
4. **Index boundary handling** (lines 87-100): Array indexing logic

### 0.3 Refactoring Strategy
1. **Start with pure functions**: Extract without changing behavior
2. **Add unit tests immediately**: Test each extracted function
3. **Gradual decoupling**: Replace direct VMEC calls with interfaces
4. **Preserve exact numerics**: Use bit-for-bit comparison tests

## Phase 1: Test Infrastructure (Week 1-2)

### 1.1 Create Golden Record Tests
- [ ] Generate reference outputs for standard VMEC test cases
  - W7-X standard configuration
  - NCSX configuration  
  - Simple axisymmetric case
- [ ] Store particle trajectories, loss times, confinement fractions
- [ ] Create automated comparison scripts with tolerance checks

### 1.2 Unit Tests for Coordinate Transformations
- [ ] Test `vmec_to_can` and `can_to_vmec` for Flux coordinates
- [ ] Test `vmec_to_boozer` and `boozer_to_vmec` transformations
- [ ] Verify Jacobian calculations and metric tensor components
- [ ] Test edge cases (axis, boundary, high aspect ratio)

### 1.3 Field Evaluation Tests
- [ ] Create mock field implementations for testing
- [ ] Test field component interpolation accuracy
- [ ] Verify derivative calculations (first and second order)

## Phase 2: Preparatory Refactoring (Week 3-4)

### 2.1 Extract VMEC-Specific Code
- [ ] Identify all VMEC-specific calls in `get_canonical_coordinates.f90`
- [ ] Identify all VMEC-specific calls in `boozer_converter.f90`
- [ ] Create inventory of required field quantities:
  - B components in various coordinate systems
  - Metric tensor elements
  - Jacobian
  - Flux surface geometry

### 2.2 Create Adapter Layer
- [ ] Design `VmecFieldAdapter` module that provides:
  ```fortran
  module vmec_field_adapter
    ! Provides high-level VMEC-specific operations
    ! built on top of MagneticField interface
    subroutine get_flux_surface_average(field, s, quantity, result)
    subroutine get_metric_tensor(field, x, g_ij)
    subroutine get_jacobian(field, x, jac)
  end module
  ```
- [ ] Implement using existing `VmecField%evaluate` calls
- [ ] Verify adapter produces identical results to direct calls

### 2.3 Refactor Without Changing Functionality
- [ ] Replace direct `vmec_field` calls with adapter calls
- [ ] Keep VMEC-only implementation for now
- [ ] Run full test suite after each change
- [ ] Ensure bit-for-bit reproducibility

## Phase 3: Abstract Interface Implementation (Week 5-6)

### 3.1 Modify Initialization Signatures
- [ ] Update `get_canonical_coordinates` to accept `class(MagneticField)`
- [ ] Update `get_boozer_coordinates` to accept `class(MagneticField)`
- [ ] Add backward compatibility wrappers that default to `VmecField`

### 3.2 Implement Field-Agnostic Algorithms
- [ ] Replace VMEC-specific spline calls with field%evaluate
- [ ] Handle coordinate system differences:
  - VMEC: (s, theta_vmec, phi)
  - GVEC: (s, theta*, phi)
  - Need coordinate transformation layer

### 3.3 Newton Iteration Refactoring
- [ ] Abstract the Newton iteration for finding VMEC theta
- [ ] Make it work for any periodic poloidal angle
- [ ] Test convergence for different coordinate systems

## Phase 4: GVEC Support Implementation (Week 7-8)

### 4.1 Extend GvecField Implementation
- [ ] Add missing methods required by canonical coordinates:
  - Flux surface averaging
  - Metric tensor computation
  - Jacobian calculation
- [ ] Implement theta* ↔ theta_vmec transformations

### 4.2 Integration Testing
- [ ] Test Flux coordinates with GVEC fields
- [ ] Test Boozer coordinates with GVEC fields
- [ ] Compare results against VMEC for same equilibrium
- [ ] Verify conservation properties

### 4.3 Performance Optimization
- [ ] Profile field evaluation calls
- [ ] Cache frequently used quantities
- [ ] Optimize coordinate transformations

## Phase 5: Validation and Documentation (Week 9-10)

### 5.1 Comprehensive Testing
- [ ] Run full test suite with both VMEC and GVEC
- [ ] Verify particle confinement statistics
- [ ] Check conservation of invariants
- [ ] Test with realistic stellarator configurations

### 5.2 Backwards Compatibility
- [ ] Ensure old input files still work
- [ ] Verify performance hasn't degraded
- [ ] Check memory usage hasn't increased significantly

### 5.3 Documentation Update
- [ ] Update code comments and docstrings
- [ ] Document new interfaces in CLAUDE.md
- [ ] Create migration guide for users
- [ ] Add examples using GVEC with canonical coordinates

## Technical Challenges to Address

### 1. Coordinate System Mismatch
- VMEC uses straight field line coordinates
- GVEC uses theta* (geometric poloidal angle)
- Need transformation layer or unified coordinate system

### 2. Flux Surface Averaging
- Currently relies on VMEC's internal representations
- Need field-agnostic algorithm for surface integrals
- May require additional methods in MagneticField interface

### 3. Performance Considerations
- Direct VMEC calls are highly optimized
- Abstract interface adds overhead
- Need careful optimization to maintain performance

### 4. Numerical Precision
- Canonical coordinate construction involves sensitive ODEs
- Small differences can accumulate
- Need robust error control and monitoring

## Success Criteria
- [ ] All existing tests pass with < 1e-10 relative error
- [ ] Can use GVEC fields with Flux and Boozer coordinates
- [ ] Performance degradation < 10%
- [ ] Code is more modular and maintainable
- [ ] No breaking changes for existing users

## Risk Mitigation
1. **Feature flags**: Add `use_legacy_vmec = .true.` option to preserve old behavior
2. **Incremental rollout**: Test thoroughly with volunteer users before release
3. **Rollback plan**: Git tags at each phase for easy reversion
4. **Parallel development**: Keep changes in feature branch until fully validated

## Notes
- Consider reaching out to original authors for historical context
- May discover additional VMEC dependencies during implementation
- Be prepared to extend MagneticField interface if needed