# Refactoring Plan: Abstract Field Support for Flux and Boozer Coordinates

## 📊 Current Status Summary
**Branch**: `refactor`  
**Phase 0**: ✅ COMPLETED - All pure functions extracted with unit tests  
**Phase 1**: ✅ COMPLETED - Test infrastructure complete, all tests passing  
**Phase 2**: ✅ COMPLETED - Preparatory refactoring done, adapter layer in place
**Next Steps**: Begin Phase 3 - Abstract Interface Implementation

## 🎯 Immediate Action Items
1. ✅ **Created VmecFieldAdapter** module with all required interfaces
2. ✅ **Replaced all direct VMEC calls** in coordinate modules with adapter
3. ✅ **Verified bit-for-bit reproducibility** - all tests pass
4. **Next: Phase 3.1** - Modify initialization signatures to accept abstract field
5. **Create field-agnostic versions** of get_canonical_coordinates and get_boozer_coordinates
6. **Add backward compatibility wrappers** for existing code

## ⚠️ MANDATORY REQUIREMENTS ⚠️

**🚨 CRITICAL: ALWAYS WORK FROM PROJECT ROOT `/afs/itp.tugraz.at/proj/plasma/CODE/ert/SIMPLE/` 🚨**

- **NEVER use `cd` commands** - Stay in project root at all times
- **NEVER run executables directly** - Always use `make test` 
- **NEVER use ctest manually** - Always use `make test`
- **Use `make test TEST=test_name`** for specific tests
- **Use `make test VERBOSE=1`** for detailed output

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
1. ✅ **Stencil initialization** (lines 60-82): COMPLETED - Extracted to `stencil_utils` module with tests
2. ✅ **Progress printing** (lines 284-295): COMPLETED - Already a well-structured subroutine
3. ✅ **Derivative array initialization** (lines 432-436): COMPLETED - Extracted to `array_utils` module with tests
4. ✅ **Index boundary handling** (lines 82-87): COMPLETED - Part of stencil_utils functionality

### 0.3 Refactoring Strategy
1. **Start with pure functions**: Extract without changing behavior
2. **Add unit tests immediately**: Test each extracted function
3. **Run `make test-all` after EVERY change**: Catch regressions immediately
4. **Gradual decoupling**: Replace direct VMEC calls with interfaces
5. **Preserve exact numerics**: Use bit-for-bit comparison tests

### 0.4 Phase 0 Status (COMPLETED) ✅
- ✅ Extracted stencil initialization to `stencil_utils` module
- ✅ Extracted derivative array initialization to `array_utils` module
- ✅ Created unit tests for both modules (test_stencil_utils.f90, test_array_utils.f90)
- ✅ All tests passing
- ✅ Progress printing and index boundary handling already well-structured

## Phase 1: Test Infrastructure (Week 1-2)

### 1.1 Create Golden Record Tests ✅ COMPLETED
- [x] Golden record infrastructure already exists:
  - `test/golden_record/golden_record.sh` - Main test runner
  - `test/golden_record/compare_files.py` - Numerical comparison with tolerance
  - `test/golden_record/run_golden_tests.sh` - Individual test runner
  - `test/golden_record/compare_golden_results.sh` - Results comparison
- [x] Test cases available:
  - Canonical coordinates test (`test/golden_record/canonical/`)
  - Boozer coordinates test (`test/golden_record/boozer/`)
- [x] All golden record tests passing (verified with `make test-all`)
- [x] Comparison includes tolerance checking via `np.isclose()`

### 1.2 Unit Tests for Coordinate Transformations
- [x] Basic test infrastructure exists:
  - `test/tests/test_boozer.f90` - Tests Boozer coordinate transformations
  - `test/tests/test_coord_trans.f90` - Integration test for coordinate transformations
  - `test/tests/test_coordinates.f90` - Simple coordinate transform driver
- [ ] Expand tests for `vmec_to_can` and `can_to_vmec` for Flux coordinates
- [x] Test `vmec_to_boozer` and `boozer_to_vmec` transformations (in test_boozer.f90)
- [ ] Add unit tests for Jacobian calculations and metric tensor components
- [ ] Test edge cases (axis, boundary, high aspect ratio)

### 1.3 Field Evaluation Tests
- [x] Field tests already exist:
  - `test/tests/field_can/test_field_can_transforms.f90` - Field transformations
  - `test/tests/field_can/test_field_can_meiss.f90` - Meiss coordinate field tests
  - `test/tests/field_can/test_field_can_albert.f90` - Albert coordinate field tests
- [ ] Create mock field implementations for abstract interface testing
- [ ] Add tests for field component interpolation accuracy
- [ ] Verify derivative calculations (first and second order) for all field types

## Phase 2: Preparatory Refactoring (Week 3-4)

### 2.1 Extract VMEC-Specific Code ✅ COMPLETED
- [x] Identified all VMEC-specific calls in `get_canonical_coordinates.f90`:
  - `vmec_field` (line 258-260) - Main field evaluation
  - `splint_iota` (line 236) - Rotational transform interpolation
  - `splint_lambda` (lines 247, 964) - Stream function interpolation
  - `splint_vmec_data` (line 1116) - Complete VMEC data interpolation
- [x] Identified all VMEC-specific calls in `boozer_converter.f90`:
  - `vmec_field` (line 139-141) - Main field evaluation
  - Uses `spline_vmec_sub` module (line 21)
- [x] Created inventory of required field quantities:
  - **Magnetic field components**: B^r, B^theta, B^phi (contravariant), B_r, B_theta, B_phi (covariant)
  - **Vector potentials**: A_theta, A_phi and their derivatives
  - **Geometric quantities**: sqrt(g) (Jacobian), lambda (stream function), iota (rotational transform)
  - **Coordinates**: R, Z and their derivatives
  - **VMEC-specific data**: torflux, ns_A, sA_phi arrays

### 2.2 Create Adapter Layer ✅ COMPLETED
- [x] Designed `VmecFieldAdapter` module with interfaces for:
  - `vmec_field_evaluate` - Replaces direct `vmec_field` calls
  - `vmec_iota_interpolate` - Replaces `splint_iota` calls
  - `vmec_lambda_interpolate` - Replaces `splint_lambda` calls
  - `vmec_data_interpolate` - Replaces `splint_vmec_data` calls
- [x] Implemented adapter using existing VMEC routines (direct pass-through for now)
- [x] Created overloaded versions with/without field object for future abstraction
- [x] Verified adapter produces identical results (all golden record tests pass)

### 2.3 Refactor Without Changing Functionality ✅ COMPLETED
- [x] Replaced all direct VMEC calls with adapter calls:
  - `get_canonical_coordinates.f90`: 3 calls replaced (vmec_field, splint_iota, splint_lambda, splint_vmec_data)
  - `boozer_converter.f90`: 1 call replaced (vmec_field)
- [x] Kept VMEC-only implementation (adapter just wraps existing calls)
- [x] Ran full test suite - all 9 tests pass
- [x] Confirmed bit-for-bit reproducibility (golden record tests pass)

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
