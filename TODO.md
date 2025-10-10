# TODO: Analytical GS Field Integration for SIMPLE

## Overview
Integrate libneo's analytical Grad-Shafranov equilibrium solver (with TF ripple) into SIMPLE's three-level coordinate/field system. Focus on Meiss canonical coordinates on top of analytical GS.

## Three Levels in SIMPLE
1. **Reference Coordinates**: Define flux surfaces (VMEC, GEOFLUX, or **ANALYTICAL GS**)
2. **Evaluation Field**: Compute B(x) at points (VMEC, GEOFLUX, coils, GVEC, or **ANALYTICAL GS**)
3. **Canonicalized Coordinates**: Transform for integration (focus on **MEISS = 3**)

## Context from Existing Work
- Meiss coordinates already work on VMEC and GEOFLUX
- `coordinate_system_t` from libneo provides unified geometry interface
- `geoflux_ready()` pattern distinguishes VMEC vs GEOFLUX
- Golden record tests validate system-level behavior

## Tasks

### 1. Add ANALYTICAL constant to magfie module
**File**: `src/magfie.f90`

- [ ] Add constant after GEOFLUX:
  ```fortran
  integer, parameter :: TEST=-1, CANFLUX=0, VMEC=1, BOOZER=2, MEISS=3, ALBERT=4, GEOFLUX=5, ANALYTICAL=6
  ```
- [ ] Verify export in module

**Test**: `grep "ANALYTICAL=6" src/magfie.f90`

### 2. Create tokamak.in namelist file
**File**: `examples/tokamak/tokamak.in`

- [ ] Define namelist structure:
  ```fortran
  &tokamak_params
    ! Equilibrium parameters (Cerfon-Freidberg)
    R0 = 6.2d0        ! Major radius [m]
    epsilon = 0.32d0  ! Inverse aspect ratio
    kappa = 1.0d0     ! Elongation (1.0 = circular)
    delta = 0.0d0     ! Triangularity (0.0 = no triangularity)
    A_param = -0.142d0  ! Shafranov parameter (pressure)
    B0 = 5.3d0        ! Toroidal field on axis [T]

    ! TF ripple parameters (optional)
    Nripple = 0       ! Number of TF coils (0 = axisymmetric)
    delta0 = 0.0d0    ! Ripple amplitude
    alpha0 = 2.0d0    ! Ripple radial exponent
    a0 = 1.984d0      ! Ripple reference radius [m] (= epsilon*R0)
    z0 = 0.0d0        ! Ripple vertical center [m]
  /
  ```

- [ ] Add example with 9-coil ripple (delta0=0.10, Nripple=9)
- [ ] Document parameters in comments

**Test**: `cat examples/tokamak/tokamak.in`

### 3. Add tokamak parameters to params module
**File**: `src/params.f90`

- [ ] Add module variables after existing field params:
  ```fortran
  ! Analytical tokamak parameters
  real(dp) :: tok_R0 = 6.2d0, tok_epsilon = 0.32d0
  real(dp) :: tok_kappa = 1.0d0, tok_delta = 0.0d0
  real(dp) :: tok_A_param = -0.142d0, tok_B0 = 5.3d0
  integer :: tok_Nripple = 0
  real(dp) :: tok_a0 = 1.984d0, tok_alpha0 = 2.0d0
  real(dp) :: tok_delta0 = 0.0d0, tok_z0 = 0.0d0
  character(1000) :: tokamak_input = 'tokamak.in'
  ```

- [ ] Add to config namelist:
  ```fortran
  namelist /config/ ..., tokamak_input, &
    tok_R0, tok_epsilon, tok_kappa, tok_delta, tok_A_param, tok_B0, &
    tok_Nripple, tok_a0, tok_alpha0, tok_delta0, tok_z0
  ```

- [ ] Add read_tokamak_config subroutine:
  ```fortran
  subroutine read_tokamak_config
    logical :: exists
    namelist /tokamak_params/ tok_R0, tok_epsilon, tok_kappa, tok_delta, &
      tok_A_param, tok_B0, tok_Nripple, tok_a0, tok_alpha0, tok_delta0, tok_z0

    inquire(file=trim(tokamak_input), exist=exists)
    if (exists) then
      open(1, file=trim(tokamak_input), status='old', action='read')
      read(1, nml=tokamak_params)
      close(1)
      print *, 'Loaded tokamak parameters from ', trim(tokamak_input)
    else
      print *, 'Using default tokamak parameters (no tokamak.in found)'
    end if
  end subroutine read_tokamak_config
  ```

- [ ] Call from read_config when `isw_field_type == ANALYTICAL` or `field_input` contains "analytical"

**Test**: Compile, check namelist reads correctly

### 4. Rewrite field_analytical_gs for coordinate system interface
**File**: `src/field/field_analytical_gs.f90`

Currently uses simple circular mapping. Need to:

- [ ] Import libneo coordinate utilities:
  ```fortran
  use analytical_tokamak_field, only: analytical_circular_eq_t
  use libneo_coordinates, only: coordinate_system_t
  ```

- [ ] Add coordinate system to type:
  ```fortran
  type, extends(MagneticField) :: AnalyticalGSField
      type(analytical_circular_eq_t) :: eq
      class(coordinate_system_t), allocatable :: coords
      logical :: initialized = .false.
  contains
      procedure :: evaluate
      procedure :: init_coordinates
  end type
  ```

- [ ] Implement proper flux surface mapping:
  - Define flux label s = ψ_normalized (0 at axis, 1 at separatrix)
  - Map (s, theta, phi) → (R, phi, Z) following flux surfaces
  - Compute metric tensor from coordinate Jacobian

- [ ] Keep eval_bfield_ripple for evaluation but transform coordinates first

**Note**: This may require extending libneo's analytical_tokamak_field to provide coordinate_system_t interface. Check if `make_analytical_coordinate_system` exists in libneo first.

**Test**: Check libneo has analytical coordinate system support

### 5. Add analytical field initialization to field.F90
**File**: `src/field.F90`

- [ ] Add import:
  ```fortran
  use field_analytical_gs, only: AnalyticalGSField, create_analytical_gs_field
  ```

- [ ] Extend `field_from_file` detection:
  ```fortran
  else if (index(filename, 'analytical') > 0 .or. index(filename, 'tokamak') > 0) then
      call create_analytical_gs_field_from_params(field)
  ```

- [ ] Implement `create_analytical_gs_field_from_params`:
  ```fortran
  subroutine create_analytical_gs_field_from_params(field)
      use params, only: tok_R0, tok_epsilon, tok_kappa, tok_delta, &
                        tok_A_param, tok_B0, tok_Nripple, tok_a0, &
                        tok_alpha0, tok_delta0, tok_z0

      class(MagneticField), allocatable, intent(out) :: field
      class(AnalyticalGSField), allocatable :: gs_temp

      call create_analytical_gs_field(tok_R0, tok_epsilon, &
          kappa=tok_kappa, delta=tok_delta, &
          A_param=tok_A_param, B0=tok_B0, &
          Nripple=tok_Nripple, a0_ripple=tok_a0, &
          alpha0=tok_alpha0, delta0=tok_delta0, z0=tok_z0, &
          gs_field=gs_temp)

      call move_alloc(gs_temp, field)
  end subroutine
  ```

**Test**: Compile and check factory works

### 6. Add magfie_analytical to magfie.f90
**File**: `src/magfie.f90`

Follow `magfie_geoflux` pattern closely:

- [ ] Add subroutine after `magfie_geoflux`:
  ```fortran
  subroutine magfie_analytical(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    real(dp), intent(in) :: x(3)  ! (s, theta, phi)
    real(dp), intent(out) :: bmod, sqrtg
    real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)

    ! Similar structure to magfie_geoflux:
    ! 1. Clip/validate coordinates
    ! 2. Map (s,theta,phi) → (R,phi,Z) via analytical equilibrium
    ! 3. Evaluate B field with ripple at (R,phi,Z)
    ! 4. Compute metric tensor from coordinate Jacobian
    ! 5. Transform B to covariant components in (s,theta,phi)
    ! 6. Compute derivatives via finite differences
    ! 7. Compute curl from derivatives

    ! TODO: Implement following geoflux_eval_point pattern
    error stop 'magfie_analytical: not yet implemented'
  end subroutine magfie_analytical
  ```

- [ ] Add case to `init_magfie`:
  ```fortran
  case(ANALYTICAL)
    magfie => magfie_analytical
  ```

- [ ] Add helper functions following geoflux pattern:
  - `analytical_eval_point` (full evaluation with metric)
  - `analytical_eval_basic` (just B field)

**Test**: Compile, check init_magfie accepts ANALYTICAL

### 7. Add analytical_ready() function
**File**: `src/field/field_analytical_gs.f90`

Follow `geoflux_ready()` pattern:

- [ ] Add module variable:
  ```fortran
  logical, save :: analytical_initialized = .false.
  ```

- [ ] Add public function:
  ```fortran
  function analytical_ready()
    logical :: analytical_ready
    analytical_ready = analytical_initialized
  end function
  ```

- [ ] Set flag in `create_analytical_gs_field`

**Test**: Function callable from other modules

### 8. Update simple_main.f90 for analytical field
**File**: `src/simple_main.f90`

- [ ] Check if analytical case needs special handling in `init_field_can`
- [ ] Ensure analytical field loaded when `field_input` contains "analytical"
- [ ] Verify `init_magfie(ANALYTICAL)` called correctly

**Test**: Run with analytical field

### 9. Test: Circular tokamak without ripple
**File**: `test/tests/test_field_analytical.f90`

- [ ] Create test similar to `test_field_geoflux.f90`:
  ```fortran
  program test_field_analytical
    use field_analytical_gs
    use params, only: tok_R0, tok_epsilon, tok_kappa, tok_delta, &
                      tok_A_param, tok_B0, tok_Nripple

    ! Set circular parameters
    tok_R0 = 6.2d0
    tok_epsilon = 0.32d0
    tok_kappa = 1.0d0
    tok_delta = 0.0d0
    tok_A_param = -0.142d0
    tok_B0 = 5.3d0
    tok_Nripple = 0  ! No ripple

    ! Create field
    ! Evaluate at test points
    ! Check B field values are reasonable
    ! Check flux surfaces are nested

    print *, 'Analytical field test PASSED'
  end program
  ```

- [ ] Add to `test/tests/CMakeLists.txt`
- [ ] Run: `ctest -R test_field_analytical`

**Test**: Must pass, <1s runtime

### 10. Test: 9-coil ripple
**File**: `test/tests/test_field_analytical_ripple.f90`

- [ ] Same as test 9 but with:
  ```fortran
  tok_Nripple = 9
  tok_delta0 = 0.10d0
  ```

- [ ] Verify ripple pattern has 9-fold symmetry
- [ ] Check peak-to-peak variation ~12-13%
- [ ] Scan toroidal angle, check periodicity

**Test**: Must pass

### 11. Test: Meiss coordinates on analytical GS
**File**: `test/tests/test_field_can_meiss_analytical.f90`

Follow `test_field_can_meiss_vmec.f90` / `test_field_can_meiss_eqdsk.f90` pattern:

- [ ] Initialize analytical field
- [ ] Initialize Meiss coordinates via `init_field_can(MEISS, field)`
- [ ] Evaluate at test points
- [ ] Check invariants:
  - Energy conservation: `H = 0.5*vpar^2 + mu*Bmod`
  - Orthogonality: `dot_product(hcov, dhth) ≈ 0`
- [ ] Compare with reference values (tolerance 1e-10)

**Test**: Must pass

### 12. Example: Circular tokamak orbit integration
**Directory**: `examples/tokamak/`

- [ ] Create `simple.in` with:
  ```
  &config
    field_input = 'analytical'
    tokamak_input = 'tokamak.in'
    isw_field_type = 3  ! Meiss
    nper = 100
    ntimstep = 1000
    ...
  /
  ```

- [ ] Create `tokamak.in` (circular, no ripple)
- [ ] Create `Makefile` with targets: `all`, `run`, `clean`
- [ ] Run and verify orbits look reasonable

**Test**: Example runs to completion

### 13. Example: 9-coil ripple orbit integration
**Directory**: `examples/tokamak_ripple/`

- [ ] Same as example 12 but with ripple enabled
- [ ] Visualize orbits showing ripple perturbation effects
- [ ] Document expected behavior

**Test**: Example runs, ripple visible in orbit traces

### 14. Documentation
**Files**: `README.md`, `examples/tokamak/README.md`

- [ ] Update main README with analytical field option
- [ ] Document tokamak.in format and parameters
- [ ] Add usage examples
- [ ] Explain coordinate system conventions
- [ ] Note compatibility with Meiss/Albert canonical coordinates

**Test**: Documentation is clear and complete

## Implementation Order

Execute tasks 1-14 sequentially. After each task:
1. Commit changes with clear message
2. Run relevant tests
3. Verify no regressions

## Key Dependencies

- libneo must have analytical_tokamak_field with field_t interface ✓ (PR #149)
- SIMPLE must have working Meiss coordinates on GEOFLUX (current TODO.md context suggests this exists)
- Coordinate system interface must support analytical equilibrium (may need libneo extension)

## Success Criteria

- [ ] All unit tests pass (`ctest --output-on-failure`)
- [ ] Analytical field works with Meiss canonical coordinates
- [ ] Ripple perturbation validated against libneo tests
- [ ] Examples run and produce physical results
- [ ] Documentation complete

## Notes

- Analytical field serves BOTH as reference (flux surfaces) AND evaluation (B field)
- Coordinate mapping must be consistent with VMEC/GEOFLUX conventions
- Focus on Meiss coordinates (MEISS=3) for canonical transforms
- Ripple is optional: Nripple=0 gives axisymmetric equilibrium
