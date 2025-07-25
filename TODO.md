# TODO: Booz_xform Integration for SIMPLE

## Overview
Add support for reading Boozer coordinate magnetic fields directly from booz_xform output files as an alternative to the current internal VMEC-to-Boozer conversion.

## Background
- Currently: SIMPLE reads VMEC files and internally converts to Boozer coordinates using `boozer_converter.F90`
- Goal: Support direct input of pre-computed Boozer fields from [booz_xform](https://github.com/hiddenSymmetries/booz_xform)
- Benefit: Potentially faster initialization and standardized Boozer transformation

## Implementation Tasks

### High Priority

1. **Review booz_xform file format** ✓
   - [x] Study booz_xform NetCDF output structure
   - [x] Document required fields and data layout
   - [x] Identify mapping to SIMPLE's internal structures

2. **Create field_booz_xform.f90 module** ✓
   - [x] Implement NetCDF reader for booz_xform files
   - [x] Map booz_xform data to SIMPLE's field structures
   - [x] Handle coordinate conventions and normalizations

3. **Add booz_xform input option** ✓
   - [x] Add new field type to `params.f90` (e.g., `isw_field_type = 5`)
   - [x] Add `booz_xform_file` parameter to namelist
   - [x] Update field type selection logic

4. **Update field initialization** ✓
   - [x] Modify `simple_main.f90` to handle booz_xform input
   - [x] Ensure compatibility with existing Boozer field evaluation
   - [x] Handle field normalization and units

5. **Integration with existing Boozer evaluation** ✓
   - [x] Connect booz_xform reader to field evaluation through magfie interface
   - [x] Implement evaluate method with VMEC→Boozer coordinate transformation
   - [x] Set up proper field component calculations

### Medium Priority

6. **Validation and testing**
   - [ ] Complete `test_booz_xform.f90` implementation
   - [ ] Compare results: VMEC→internal Boozer vs VMEC→booz_xform→SIMPLE
   - [ ] Benchmark performance differences
   - [ ] Test with various stellarator configurations

7. **Example workflow**
   - [ ] Create example: VMEC wout.nc → booz_xform → SIMPLE
   - [ ] Document command sequences
   - [ ] Provide sample input files

### Low Priority

8. **Documentation**
   - [ ] Update user manual with booz_xform option
   - [ ] Add to CLAUDE.md for AI assistance
   - [ ] Create tutorial for booz_xform workflow

## Technical Considerations

### Data Structures
- Booz_xform outputs: `|B|`, `G`, `I`, `K`, `p`, `q`, derivatives
- SIMPLE needs: Magnetic field components, metric tensor, Jacobian
- Coordinate mapping: (s, θ_B, ζ_B) in both systems

### Compatibility
- Ensure backward compatibility with existing VMEC input
- Support switching between internal and external Boozer conversion
- Handle different normalization conventions

### Performance
- Pre-computed Boozer coordinates may reduce initialization time
- Consider memory vs computation tradeoffs
- Profile both approaches for typical use cases

## Current Status
### Completed ✓
- Basic test file created: `test/tests/test_booz_xform.f90`
- CMake integration ready
- Existing Boozer infrastructure in place: `field_can_boozer.f90`, `boozer_converter.F90`
- Created test framework with setup script and comparison plots
- Implemented `field_booz_xform.f90` module to read BOOZXFORM files
- Fixed NetCDF interface issues (using nctools_module from libneo)
- Test passes and generates comparison plots
- Added BOOZXFORM=5 constant to magfie.f90
- Updated simple_main.f90 to handle the new field type (sets boozxform_filename)
- Implemented proper evaluate method in field_booz_xform.f90 with coordinate transformation
- Added boozxform_file parameter to params.f90 namelist
- Created benchmark script: `benchmark/benchmark_boozxform.py`

### In Progress
- Fixed NetCDF loading for scalar and 1D arrays
- 2D array loading still has issues ("Start+count exceeds dimension bound")
- Need to implement proper coordinate transformation (pmns/gmn interpretation)

### Known Issues
- NetCDF 2D array reading fails for Fourier coefficients (rmnc_b, zmns_b, etc.)
- Coordinate transformation between VMEC and Boozer angles incomplete
- RKF45 integrator errors due to missing field data

## Next Steps
1. Fix 2D array NetCDF loading (possibly transpose issue)
2. Implement proper VMEC↔Boozer coordinate transformation
3. Complete validation with benchmark files