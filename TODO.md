# TODO: Booz_xform Integration for SIMPLE

## Overview
Add support for reading Boozer coordinate magnetic fields directly from booz_xform output files as an alternative to the current internal VMEC-to-Boozer conversion.

## Background
- Currently: SIMPLE reads VMEC files and internally converts to Boozer coordinates using `boozer_converter.F90`
- Goal: Support direct input of pre-computed Boozer fields from [booz_xform](https://github.com/hiddenSymmetries/booz_xform)
- Benefit: Potentially faster initialization and standardized Boozer transformation

## Implementation Tasks

### High Priority

1. **Review booz_xform file format** ⏳
   - [ ] Study booz_xform NetCDF output structure
   - [ ] Document required fields and data layout
   - [ ] Identify mapping to SIMPLE's internal structures

2. **Create field_booz_xform.f90 module**
   - [ ] Implement NetCDF reader for booz_xform files
   - [ ] Map booz_xform data to SIMPLE's field structures
   - [ ] Handle coordinate conventions and normalizations

3. **Add booz_xform input option**
   - [ ] Add new field type to `params.f90` (e.g., `isw_field_type = 5`)
   - [ ] Add `booz_xform_file` parameter to namelist
   - [ ] Update field type selection logic

4. **Update field initialization**
   - [ ] Modify `simple_main.f90` to handle booz_xform input
   - [ ] Ensure compatibility with existing Boozer field evaluation
   - [ ] Handle field normalization and units

5. **Integration with existing Boozer evaluation**
   - [ ] Connect booz_xform reader to `field_can_boozer.f90`
   - [ ] Ensure `eval_field_booz` works with external data
   - [ ] Verify derivatives and metric calculations

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
- Basic test file created: `test/tests/test_booz_xform.f90`
- CMake integration ready
- Existing Boozer infrastructure in place: `field_can_boozer.f90`, `boozer_converter.F90`

## Next Steps
1. Start with reviewing booz_xform NetCDF format
2. Design data structures for booz_xform input
3. Implement minimal reader for proof of concept