# TODO: Analytical GS Field via Geoflux Framework

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

## Tasks

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
**Status**: ✅ COMPLETE (libneo commits a9e84b5, 643315e)

**Achievement**: Geoflux coordinates now field-agnostic like VMEC flux coordinates

**Changes**:
- Added `psi_evaluator_i` callback interface to `geoflux_coordinates`
- Modified `initialize_analytical_geoflux` to accept psi evaluator callback
- `psi_from_position` dispatches to callback when `use_geqdsk = .false.`
- Created `analytical_geoflux_field` module with `init_analytical_geoflux` and `splint_analytical_geoflux_field`

**Test**: `test_analytical_geoflux.x` passes - analytical GS works through geoflux coordinates

---

### Phase 3: SIMPLE Integration (IN SIMPLE REPO)

#### 5. Add initialize_analytical_geoflux to SIMPLE field.F90
**File**: `src/field.F90`

- [ ] Add import:
  ```fortran
  use geoflux_field, only: initialize_analytical_geoflux
  ```

- [ ] Extend `field_from_file`:
  ```fortran
  else if (index(filename, 'analytical') > 0 .or. index(filename, 'tokamak') > 0) then
      use params, only: tok_R0, tok_epsilon, tok_kappa, tok_delta, &
                        tok_A_param, tok_B0, tok_Nripple, tok_a0, &
                        tok_alpha0, tok_delta0, tok_z0

      call initialize_analytical_geoflux(tok_R0, tok_epsilon, tok_kappa, tok_delta, &
                                         tok_A_param, tok_B0, &
                                         tok_Nripple, tok_a0, tok_alpha0, tok_delta0, tok_z0)
      allocate(GeofluxField :: field)
  ```

**Result**: Analytical GS appears as geoflux field to rest of SIMPLE

**Test**: Load with `field_input='analytical'`, verify field object created

#### 6. Update simple_main.f90 to use GEOFLUX mode
**File**: `src/simple_main.f90`

When analytical field detected, ensure `isw_field_type` uses geoflux coordinate system:

- [ ] Check if special handling needed, or if geoflux auto-detection suffices

**Test**: Run SIMPLE with analytical field, verify geoflux coordinate system active

---

### Phase 4: Testing ✅ COMPLETE (libneo)

#### Libneo Test Hierarchy

**Unit Tests** (libneo):
- ✅ `test_analytical_circular`: Analytical GS field direct evaluation
- ✅ `test_analytical_geoflux`: Geoflux initialization and field evaluation

**Integration Tests** (libneo):
- ✅ `test_analytical_geoflux_integration`:
  * Coordinate round-trip (geoflux ↔ cylindrical)
  * Field consistency (geoflux vs direct)
  * Flux surface nesting
- ✅ `test_ripple_field`: TF ripple (9-coil, 9-fold symmetry, ~12.65% variation)

**All libneo tests pass**: `ctest -R analytical` (3/3 passed)

---

### Phase 5: SIMPLE System Tests ⏸️ READY (waiting on Phase 3)

#### 10. ✅ Create tokamak example (ITER-size, no ripple) [DONE]
**Directory**: `examples/tokamak_alpha_confinement/`
**Status**: ✅ COMPLETE (commit b024955)

**Config**:
- Field: analytical GS (ITER parameters: R0=6.2m, ε=0.32, B0=5.3T)
- Coordinates: Meiss (isw_field_type=3) on geoflux
- Particles: 128 alpha particles, E=3.5 MeV
- Start: s=0.3 (mid-radius)
- Duration: 0.001 s

**Expected**: Zero particles lost (perfect confinement without ripple)

**Files created**:
- `simple.in`: SIMPLE configuration
- `tokamak.in`: ITER parameters
- `Makefile`: Build/run automation
- `README.md`: Documentation
- `.gitignore`: Output files

#### 11. ✅ Create system test from example [DONE]
**File**: `test/tests/test_tokamak_alpha_confinement.f90`
**Status**: ✅ COMPLETE (commit 0b9a140)

**Test**:
- ✅ Scaffold created and integrated with CMake
- ✅ Marked `WILL_FAIL TRUE` until field integration complete
- ⏸️ Will activate once field.F90 updated (Phase 3 task 5)

**Checks** (when activated):
- [ ] Run example automatically
- [ ] Parse output
- [ ] Assert: `n_lost == 0`
- [ ] Assert: all particles remain at s ∈ [0.2, 0.4]

**Run**: `ctest -R tokamak_alpha` (currently expected to fail)

---

### Phase 5: Examples

#### 10. Example: Circular tokamak orbit integration
**Directory**: `examples/tokamak/`

- [ ] Create `simple.in`:
  ```
  &config
    field_input = 'analytical'
    tokamak_input = 'tokamak.in'
    isw_field_type = 3  ! Meiss (will auto-use geoflux)
    nper = 100
    ntimstep = 1000
    notrace_passing = 0
    ...
  /
  ```

- [ ] Create `Makefile` with targets: `all`, `run`, `clean`
- [ ] Run and verify orbits physically reasonable
- [ ] Check particle confinement

**Test**: Example runs to completion, produces output

#### 11. Example: 9-coil ripple orbit integration
**Directory**: `examples/tokamak_ripple/`

- [ ] Use `tokamak_ripple.in` with Nripple=9
- [ ] Run particle tracing
- [ ] Visualize orbit perturbations due to ripple
- [ ] Check ripple-trapping effects

**Test**: Ripple effects visible in trajectories

---

### Phase 6: Documentation

#### 12. Update documentation
**Files**: `README.md`, `examples/tokamak/README.md`

- [ ] Document analytical field option
- [ ] Explain tokamak.in parameters
- [ ] Add usage examples
- [ ] Note: "Uses geoflux coordinate framework internally"
- [ ] Document ripple effects on orbits

**Test**: Documentation clear and accurate

---

## Implementation Order

1. **Phase 2 in libneo** (task 4)
2. **Phase 3 in SIMPLE** (tasks 5-6)
3. **Phase 4 tests** (tasks 7-9)
4. **Phase 5 examples** (tasks 10-11)
5. **Phase 6 docs** (task 12)

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

## Success Criteria

- [ ] All tests pass (`ctest --output-on-failure`)
- [ ] Analytical field loads via geoflux framework
- [ ] Meiss coordinates work on analytical GS
- [ ] Ripple effects validated (9-fold symmetry, ~12-13% variation)
- [ ] Examples run and produce physical orbits
- [ ] No new coordinate system code needed (pure geoflux reuse)

## Notes

- **No ANALYTICAL field type needed** - just use GEOFLUX with analytical data
- **No new magfie_analytical** - magfie_geoflux handles everything
- **Ripple works automatically** - included in B field evaluation on grid
- This approach is **much simpler** than implementing a new coordinate system
