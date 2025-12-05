# Coordinate Systems and Field Formats in SIMPLE

## Field Formats Supported

### 1. VMEC (.nc files)
- **Primary stellarator equilibrium format**
- NetCDF files from VMEC code
- Contains magnetic field data in flux coordinates (s, θ, φ)
- Implemented in `src/field/field_vmec.f90` via `VmecField` class
- Uses spline interpolation from `spline_vmec_sub`

### 2. Coil Files (.coils or files starting with "coils")
- **External coil magnetic fields**
- Uses Biot-Savart calculations via `neo_biotsavart` library
- Handled by `src/field/field_coils.f90` with 3D spline interpolation
- Supports direct evaluation and pre-computed spline grids

### 3. GVEC Files (.dat)
- **Magnetic field data using GVEC B-spline library**
- Requires `GVEC_AVAILABLE` compilation flag (`-DENABLE_GVEC=ON`)
- Processed by `src/field/field_gvec.f90` using minimal GVEC from `thirdparty/gvec/`
- Provides B-spline and cubic spline functionality for field interpolation

## Coordinate Systems

### Base Coordinates (Reference System)
**VMEC coordinates**: (r=√s, θ, φ)
- `r = sqrt(s)` where `s` is normalized toroidal flux
- `θ` is poloidal angle (VMEC poloidal angle, not θ*)
- `φ` is toroidal angle
- Used as reference system for all coordinate transformations

### Canonical Coordinates
Five canonical coordinate variants selected via `isw_field_type` parameter:

1. **Test** (`field_can_test.f90`)
   - Simple analytical test field for validation

2. **Flux** (`field_can_flux.f90`) 
   - Canonical flux coordinates
   - Direct canonical transformation

3. **Boozer** (`field_can_boozer.f90`)
   - Boozer coordinates (θ_B, ζ_B)
   - Straight field line coordinates
   - Uses `boozer_sub` module for coordinate conversion

4. **Meiss** (`field_can_meiss.f90`)
   - Meiss canonical coordinates with gauge transformation
   - Requires non-canonical field input for initialization
   - Uses 3D spline interpolation for field components
   - Includes lambda (angle difference) and chi (gauge) transformations

5. **Albert** (`field_can_albert.f90`)
   - Albert canonical coordinates with ψ transformation  
   - Similar structure to Meiss but with psi coordinate transformation
   - Uses `psi_transform` module for coordinate mapping

### Physical Coordinates
**Cylindrical coordinates**: (R, φ, Z)
- Standard cylindrical coordinates
- Intermediate step for Cartesian conversion

**Cartesian coordinates**: (X, Y, Z)
- Standard Cartesian coordinates
- Final coordinate system for visualization

**Transformation Chain**: VMEC → Cylindrical → Cartesian
- Implemented in `src/coordinates/coordinates.f90`
- Supports both forward transforms and Jacobian computation
- Function pointer system for runtime coordinate selection

## Implementation Architecture

### Abstract Field Interface
- Base class: `magnetic_field_t` in `src/field/field_base.f90`
- Abstract `evaluate()` method for polymorphic field evaluation
- Common interface: `evaluate(x, Acov, hcov, Bmod, sqgBctr)`

### Field Selection Logic
Runtime field format detection in `field_from_file()` (`src/field.F90`):
```fortran
if (endswith(filename, '.nc')) then
    allocate(VmecField :: field_from_file)
else if (startswidth(stripped_name, 'coils') .or. endswith(filename, '.coils')) then
    field_from_file = create_coils_field(filename)
else if (endswith(filename, '.dat')) then
    field_from_file = create_gvec_field(filename)
```

### Canonical Coordinate Initialization
- Function: `init_field_can(field_id, field_noncan)` in `src/field_can.f90`
- Field-agnostic initialization via polymorphic `magnetic_field_t` interface
- Meiss and Albert coordinates require non-canonical field input for spline construction

### Coordinate Transformation System
- Function pointers for runtime coordinate selection
- `integ_to_ref` and `ref_to_integ` transformation procedures
- Identity transforms for test/flux coordinates
- Specialized transforms for Boozer/Meiss/Albert coordinates

## Usage Pattern

1. **Field Loading**: File extension determines field type via `field_from_file()`
2. **Canonical Setup**: `init_field_can()` initializes chosen canonical system
3. **Runtime Evaluation**: Polymorphic `evaluate()` calls appropriate field implementation
4. **Coordinate Conversion**: Function pointers enable runtime coordinate transformations
5. **Integration**: Symplectic integrators work in canonical coordinates with field evaluation

## Key Files

### Core Field Modules
- `src/field/field_base.f90` - Abstract base class
- `src/field.F90` - Field factory and format detection
- `src/field_can.f90` - Canonical coordinate system manager

### Field Implementations  
- `src/field/field_vmec.f90` - VMEC equilibria
- `src/field/field_coils.f90` - External coil fields
- `src/field/field_gvec.f90` - GVEC B-spline fields

### Canonical Coordinates
- `src/field/field_can_*.f90` - Individual canonical coordinate implementations
- `src/get_canonical_coordinates.F90` - Canonical coordinate computation

### Coordinate Transformations
- `src/coordinates/coordinates.f90` - Physical coordinate transformations
- `src/field/psi_transform.f90` - ψ coordinate transformation utilities
