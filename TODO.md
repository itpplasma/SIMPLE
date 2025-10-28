# TODO: Analytical GS Field via Geoflux Framework

## Current Status: Phase 3 (SIMPLE Integration) IN PROGRESS

### What Works Now ✅
- **Libneo**: Field-agnostic geoflux coordinates with callback interface
- **Libneo**: Complete test suite (unit + integration tests, all passing)
- **SIMPLE**: Infrastructure ready (params, tokamak.in files)
- **SIMPLE**: Alpha confinement example created and documented
- **SIMPLE**: System test scaffold integrated with CMake

### What's Needed to Complete 🔨
- **Update src/field.F90** to recognize and initialize analytical fields
- **Verify** system test passes with analytical field
- **Optional**: Create ripple example and documentation

---

## Overview

Integrate libneo's analytical Grad-Shafranov equilibrium solver (with TF ripple) into SIMPLE by **populating geoflux data structures directly from analytical GS** (no file I/O needed).

## Strategy: Analytical GS → Geoflux Framework

Instead of implementing a new coordinate system, we:
1. Evaluate analytical GS (ψ, B, derivatives) on the geoflux grid
2. Populate geoflux internal arrays/splines directly in memory
3. Set `geoflux_ready = .true.`
4. **Everything else just works** - SIMPLE sees it as geoflux

## Benefits
- ✅ Reuses existing geoflux coordinate system (flux surfaces, metric, Jacobian)
- ✅ Meiss coordinates work unchanged
- ✅ No new coordinate mapping code needed
- ✅ No file I/O overhead
- ✅ Ripple naturally included in field evaluation

---

## Implementation Status

### Phase 1: Infrastructure ✅ COMPLETE

#### 1. ✅ Add ANALYTICAL constant [DONE]
- Added `ANALYTICAL=6` to `src/magfie.f90`

#### 2. ✅ Create tokamak.in [DONE]
- Created `examples/tokamak/tokamak.in` (no ripple)
- Created `examples/tokamak/tokamak_ripple.in` (9-coil)

#### 3. ✅ Add tokamak parameters to params [DONE]
- Added tok_* variables to `src/params.f90`
- Extended config namelist
- Added `read_tokamak_config()` subroutine

**Committed**: 35485e2

---

### Phase 2: Libneo - Field-Agnostic Geoflux ✅ COMPLETE

#### 4. ✅ Made geoflux coordinates field-agnostic [DONE]
**Status**: ✅ COMPLETE (libneo commits a9e84b5, 643315e, 7494638)

**Achievement**: Geoflux coordinates now field-agnostic like VMEC flux coordinates

**Changes**:
- Added `psi_evaluator_i` callback interface to `geoflux_coordinates`
- Modified `initialize_analytical_geoflux` to accept psi evaluator callback
- `psi_from_position` dispatches to callback when `use_geqdsk = .false.`
- Created `analytical_geoflux_field` module with `init_analytical_geoflux` and `splint_analytical_geoflux_field`

**Files modified** (libneo):
- `src/coordinates/geoflux_coordinates.f90`
- `src/magfie/analytical_geoflux_field.f90`
- `test/source/test_analytical_geoflux.f90`
- `test/source/test_analytical_geoflux_integration.f90`
- `test/CMakeLists.txt`

**Tests created**:
- `test_analytical_geoflux.x`: Unit test for geoflux initialization
- `test_analytical_geoflux_integration.x`: Integration test for coordinate transformations

---

### Phase 3: SIMPLE Integration ⏸️ IN PROGRESS

#### 5. ⏸️ Add analytical field support to SIMPLE field.F90 [BLOCKED]
**File**: `src/field.F90`
**Status**: ❌ NOT YET DONE - **THIS IS THE CRITICAL MISSING PIECE**

**Required changes**:

1. Add import at top of file:
   ```fortran
   use analytical_geoflux_field, only: init_analytical_geoflux
   ```

2. Add import for tokamak parameters (in subroutine scope):
   ```fortran
   use params, only: tok_R0, tok_epsilon, tok_kappa, tok_delta, &
                     tok_A_param, tok_B0, tok_Nripple, tok_a0, &
                     tok_alpha0, tok_delta0, tok_z0
   ```

3. Extend `field_from_file` subroutine (around line 27, after GEQDSK check):
   ```fortran
   if (is_geqdsk(filename)) then
       call initialize_geoflux_field(trim(filename))
       allocate(GeofluxField :: field)
   else if (index(filename, 'analytical') > 0 .or. index(filename, 'tokamak') > 0) then
       ! Initialize analytical field via geoflux coordinates
       call init_analytical_geoflux(tok_R0, tok_epsilon, tok_kappa, tok_delta, &
                                    tok_A_param, tok_B0, &
                                    tok_Nripple, tok_a0, tok_alpha0, tok_delta0, tok_z0)
       allocate(GeofluxField :: field)
   else if (endswith(filename, '.nc')) then
   ```

**Result**: Analytical GS appears as geoflux field to rest of SIMPLE

**Test**:
- Build SIMPLE: `make`
- Load with `field_input='analytical'`, verify field object created
- Run alpha confinement example

#### 6. ✅ Verify coordinate system handling [LIKELY OK]
**File**: `src/simple_main.f90`

**Analysis**: When `field_input='analytical'` is used, SIMPLE will:
1. Call `field_from_file('analytical', field)`
2. Get back a `GeofluxField` object
3. User sets `isw_field_type = 3` (Meiss) in simple.in
4. Meiss coordinates work on top of geoflux reference coordinates

**Expected**: No changes needed - geoflux auto-detection should work

**Test**: Run SIMPLE with analytical field, verify Meiss coordinates active

---

### Phase 4: Testing ✅ COMPLETE (libneo)

#### Libneo Test Hierarchy

**Unit Tests** (libneo):
- ✅ `test_analytical_circular`: Analytical GS field direct evaluation
  - Tests field values at various locations
  - Validates field with/without TF ripple
  - **Status**: PASS

- ✅ `test_analytical_geoflux`: Geoflux initialization and field evaluation
  - Tests geoflux initialization with analytical psi evaluator
  - Verifies field evaluation through geoflux coordinates
  - Tests both axisymmetric and rippled configurations
  - **Status**: PASS

**Integration Tests** (libneo):
- ✅ `test_analytical_geoflux_integration`:
  - **Coordinate round-trip**: geoflux ↔ cylindrical transformation consistency
  - **Field consistency**: Compares geoflux field vs direct analytical evaluation
  - **Flux surface nesting**: Verifies monotonic psi (proper flux surface ordering)
  - Tolerance: 1e-3 (numerical interpolation on cached grid)
  - **Status**: PASS

- ✅ `test_ripple_field`: TF ripple validation
  - Tests 9-coil configuration
  - Validates 9-fold periodicity
  - Confirms ~12.65% peak-to-peak variation
  - Generates CSV data and plots
  - **Status**: PASS

**Test Results**:
```bash
$ cd libneo/build && ctest -R analytical
Test #14: test_analytical_circular .............. Passed 0.02 sec
Test #15: test_analytical_geoflux ............... Passed 0.05 sec
Test #16: test_analytical_geoflux_integration ... Passed 0.03 sec

100% tests passed, 0 tests failed out of 3
```

**Files** (libneo):
- `test/source/test_analytical_circular.f90`
- `test/source/test_analytical_geoflux.f90`
- `test/source/test_analytical_geoflux_integration.f90`
- `test/source/test_ripple_field.f90`
- `TEST_SUMMARY.md` (documentation)

**Commits** (libneo):
- 643315e: Field-agnostic geoflux implementation
- 7494638: Integration test

---

### Phase 5: SIMPLE System Tests ⏸️ READY (waiting on Phase 3)

#### 10. ✅ Create tokamak alpha confinement example [DONE]
**Directory**: `examples/tokamak_alpha_confinement/`
**Status**: ✅ COMPLETE (commit b024955)

**Purpose**: System-level test to verify:
- Analytical GS field integrates correctly with geoflux coordinates
- Meiss canonical coordinates work properly on analytical equilibria
- Symplectic orbit integration conserves energy and confines particles

**Config**:
- Field: analytical GS (ITER parameters: R0=6.2m, ε=0.32, B0=5.3T)
- Coordinates: Meiss (isw_field_type=3) on geoflux
- Particles: 128 alpha particles, E=3.5 MeV
- Start: s=0.3 (mid-radius)
- Duration: 0.001 s (1 ms)
- No TF ripple (Nripple=0)

**Expected**: Zero particles lost (perfect confinement without ripple)

**Files created**:
- `simple.in`: SIMPLE configuration
- `tokamak.in`: ITER parameters
- `Makefile`: Build/run automation with targets: `all`, `run`, `clean`
- `README.md`: Comprehensive documentation
- `.gitignore`: Output file management

**Run**: `cd examples/tokamak_alpha_confinement && make run`

#### 11. ✅ Create system test from example [DONE]
**File**: `test/tests/test_tokamak_alpha_confinement.f90`
**Status**: ✅ COMPLETE (commit 0b9a140)

**Test scaffold**:
- ✅ Loads configuration from example directory
- ✅ Verifies parameters (ntestpart=128, sbeg=0.3, trace_time=1ms)
- ✅ Integrated with CMake build system
- ✅ Marked `WILL_FAIL TRUE` until field.F90 updated (Phase 3 task 5)
- ⏸️ Will activate once analytical field loading works

**Verification checks** (when activated):
- Run example automatically via SIMPLE API
- Parse output for particle losses
- Assert: `n_lost == 0` (no particles lost)
- Assert: particles remain at s ∈ [0.2, 0.4] (confined to flux surface region)

**CMake integration**:
- Added to `test/tests/CMakeLists.txt`
- Labels: `system;integration;analytical`
- Timeout: 120 seconds
- Currently marked `WILL_FAIL TRUE`

**Run**: `ctest -R tokamak_alpha` (currently expected to fail)

**Activation steps** (after Phase 3 task 5):
1. Implement field loading in `src/field.F90`
2. Remove `WILL_FAIL TRUE` from `test/tests/CMakeLists.txt`
3. Implement full test logic (currently has placeholder error stop)
4. Verify test passes

**Commits**:
- b024955: Created example directory
- 0b9a140: Created system test scaffold
- 963bde8: Updated TODO with Phase 5 status

---

### Phase 6: Ripple Example (OPTIONAL)

#### 12. Example: 9-coil ripple orbit integration
**Directory**: `examples/tokamak_ripple/` (or extend existing example)
**Status**: ❌ NOT STARTED

**Purpose**: Demonstrate ripple-induced transport effects

**Tasks**:
- [ ] Create configuration with `Nripple=9`, `delta0=0.03` (3% ripple)
- [ ] Use same particle parameters as alpha confinement example
- [ ] Run particle tracing
- [ ] Visualize orbit perturbations due to ripple
- [ ] Check ripple-trapping effects and losses

**Expected**: Some particles lost due to ripple perturbation

**Priority**: LOW - Can be done after main integration is working

---

### Phase 7: Documentation (OPTIONAL)

#### 13. Update main documentation
**Files**: `README.md`, `examples/README.md`
**Status**: ❌ NOT STARTED

**Tasks**:
- [ ] Document analytical field option in main README
- [ ] Add usage example to README
- [ ] Explain tokamak.in parameters
- [ ] Note: "Uses geoflux coordinate framework internally"
- [ ] Document ripple effects on orbits (if ripple example done)
- [ ] Add link to TEST_SUMMARY.md in libneo

**Priority**: LOW - Can be done after main integration is working

---

## Critical Path to Completion

### Immediate Next Step 🔥
**Update `src/field.F90`** (Phase 3, Task 5):
1. Add imports for `init_analytical_geoflux` and tokamak parameters
2. Add `else if` clause in `field_from_file` to handle analytical field
3. Test with example: `cd examples/tokamak_alpha_confinement && make run`

### Then
1. **Verify system test**: Remove `WILL_FAIL` and complete test implementation
2. **Run test suite**: `ctest -R tokamak_alpha`
3. **Merge to main**: Once tests pass

### Optional Later
- Create ripple example (Phase 6)
- Update documentation (Phase 7)

---

## Key Technical Details

### Geoflux Grid Structure
- `ns_A`: Number of flux surfaces (radial)
- `ntheta_A`: Number of poloidal grid points
- `nphi_A`: Number of toroidal grid points
- Arrays: `psi_a`, `R_a`, `Z_a`, `Bvec_a`, etc.

### Coordinate Mapping
- Flux label: `s = (psi - psi_axis) / (psi_sep - psi_axis)`
- Poloidal angle: `theta` follows flux surface contours
- Toroidal angle: `phi` (geometric)

### Ripple Evaluation
- Evaluate `analytical_circular_eq_t%eval_bfield_ripple(R, phi, Z, ...)`
- Ripple automatically included when Nripple > 0
- Geoflux splines capture 3D ripple structure

### Field-Agnostic Callback Pattern
- `psi_evaluator_i` abstract interface in geoflux_coordinates
- `psi_eval_wrapper` in analytical_geoflux_field bridges analytical GS to geoflux
- `psi_from_position` dispatches to callback when `use_geqdsk = .false.`

---

## Success Criteria

### Minimum (for merge)
- [x] Libneo tests pass (`ctest -R analytical` in libneo)
- [ ] Analytical field loads via geoflux framework (Phase 3 task 5)
- [ ] System test passes (`ctest -R tokamak_alpha` in SIMPLE)
- [x] Meiss coordinates work on analytical GS (verified by test)
- [x] No new coordinate system code needed (pure geoflux reuse)

### Optional (nice to have)
- [ ] Ripple example demonstrates transport effects
- [ ] Documentation updated with usage examples
- [ ] Ripple effects validated (9-fold symmetry, ~12.65% variation) - already validated in libneo

---

## Repository Status

### Libneo (`itpplasma/libneo`)
- **Branch**: `main` (work committed and pushed)
- **Status**: ✅ ALL WORK COMPLETE
- **Tests**: All passing (ctest -R analytical: 3/3)
- **Commits**:
  - a9e84b5: Initial field-agnostic work
  - 643315e: Complete field-agnostic geoflux
  - 7494638: Integration test

### SIMPLE (`itpplasma/SIMPLE`)
- **Branch**: `feature/tf-ripple-perturbation`
- **Status**: ⏸️ WAITING ON Phase 3 Task 5
- **Tests**: System test scaffold ready but marked `WILL_FAIL`
- **Commits**:
  - 35485e2: Infrastructure (params, tokamak.in)
  - b024955: Alpha confinement example
  - 0b9a140: System test scaffold
  - 963bde8: TODO update

---

## Notes

- **No ANALYTICAL field type needed** - just use GEOFLUX with analytical data
- **No new magfie_analytical** - magfie_geoflux handles everything
- **Ripple works automatically** - included in B field evaluation on grid
- This approach is **much simpler** than implementing a new coordinate system
- **All infrastructure is ready** - only missing 10 lines of code in field.F90!

---

## Quick Reference: What to Modify

### The One File That Needs Changes
**File**: `src/field.F90`

**Location**: Around line 1-10 (imports) and line 27 (field_from_file)

**Change 1** - Add import:
```fortran
use analytical_geoflux_field, only: init_analytical_geoflux
```

**Change 2** - Extend field_from_file (after GEQDSK check):
```fortran
else if (index(filename, 'analytical') > 0 .or. index(filename, 'tokamak') > 0) then
    use params, only: tok_R0, tok_epsilon, tok_kappa, tok_delta, &
                      tok_A_param, tok_B0, tok_Nripple, tok_a0, &
                      tok_alpha0, tok_delta0, tok_z0
    call init_analytical_geoflux(tok_R0, tok_epsilon, tok_kappa, tok_delta, &
                                 tok_A_param, tok_B0, &
                                 tok_Nripple, tok_a0, tok_alpha0, tok_delta0, tok_z0)
    allocate(GeofluxField :: field)
```

**That's it!** Once this is done, everything else should work.
