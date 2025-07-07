# TODO: Fix Compiler Warnings

## High Priority (Low Hanging Fruit)

### 1. Fix real-to-integer conversion warnings in samplers.f90
- [ ] Line 176: Fix conversion from REAL(8) to INTEGER(4)
- [ ] Line 183: Fix conversion from REAL(8) to INTEGER(4)
- **Solution**: Use explicit `int()` or `nint()` functions for proper conversion

## Medium Priority

### 2. Fix real equality/inequality comparisons in minpack.f90
- [ ] Line 139: Replace equality comparison for REAL(8)
- [ ] Line 149: Replace equality comparison for REAL(8) (2 instances)
- [ ] Line 280: Replace inequality comparison for REAL(8)
- [ ] Line 288: Replace equality comparison for REAL(8)
- [ ] Line 304: Replace equality comparison for REAL(8)
- [ ] Line 330: Replace inequality comparison for REAL(8)
- [ ] Line 346: Replace inequality comparison for REAL(8)
- [ ] Line 358: Replace inequality comparison for REAL(8)
- [ ] Line 517: Replace inequality comparison for REAL(8)
- [ ] Line 540: Replace inequality comparison for REAL(8)
- [ ] Line 544: Replace inequality comparison for REAL(8)
- [ ] Line 673: Replace equality comparison for REAL(8)
- [ ] Line 679: Replace equality comparison for REAL(8)
- [ ] Line 699: Replace equality comparison for REAL(8)
- [ ] Line 706: Replace equality comparison for REAL(8)
- [ ] Line 717: Replace equality comparison for REAL(8)
- [ ] Line 839: Replace equality comparison for REAL(8)
- [ ] Line 845: Replace equality comparison for REAL(8)
- [ ] Lines 1080-1535: Multiple real comparisons to fix
- [ ] Lines 1732-1851: Multiple real comparisons to fix
- [ ] Lines 4707-4901: Multiple real comparisons to fix
- [ ] Lines 5364-5474: Multiple real comparisons to fix
- **Solution**: Use tolerance-based comparisons with `abs(a - b) < epsilon`

### 3. Fix implicit interface warnings in minpack.f90
- [ ] Line 326: Add explicit interface for 'enorm' procedure
- [ ] Line 330: Add explicit interface for 'enorm' procedure
- [ ] Line 1080: Add explicit interface for 'fcn' procedure
- [ ] Line 1087: Add explicit interface for 'fcn' procedure
- [ ] Line 1110: Add explicit interface for 'fcn' procedure
- [ ] Line 1121: Add explicit interface for 'fcn' procedure
- [ ] Line 1210: Add explicit interface for 'enorm' procedure
- [ ] Line 1237: Add explicit interface for 'fdjac1' procedure
- **Solution**: Add interface blocks or use modules with explicit interfaces

## Low Priority

### 4. Fix missing terminating character warnings in GVEC library
- [ ] These are in third-party dependency code (build/_deps/gvec-src/)
- [ ] Multiple instances in comments with unterminated quotes
- **Note**: These are in external dependency, may not need fixing

## Notes
- Total warnings to fix: ~50+
- Most warnings are in `minpack.f90` (legacy MINPACK library)
- Real comparison warnings can be fixed using epsilon-based comparisons
- Implicit interface warnings require adding interface blocks
- Consider using compiler flags to suppress warnings in legacy code if needed
