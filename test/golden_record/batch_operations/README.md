# Batch Operations Golden Record Test

This test case validates that the Python batch API produces identical results to individual particle processing.

## Test Configuration

- **Particles**: 1000 test particles
- **Integration**: Symplectic method (integmode=3)
- **Time steps**: 200 steps with dtau=1e-3
- **Classification**: Enabled (ntcut=10)
- **Coordinate system**: Canonical coordinates (isw_field_type=1)

## Validation Criteria

1. **Deterministic Results**: Same particle initial conditions must produce identical final states
2. **Batch vs Individual**: Batch processing must match sequential individual particle processing
3. **Statistical Equivalence**: Confinement statistics must be identical within numerical precision
4. **Memory Layout**: SoA access must not affect numerical results

## Expected Outputs

- `confined_fraction.dat`: Time evolution of confined particle fraction
- `times_lost.dat`: Individual particle loss times and classification
- `class_parts.dat`: Particle classification results
- `start.dat`: Initial particle conditions (for reproducibility)

## Performance Requirements

- Python batch API overhead < 5% vs direct Fortran execution
- Memory access must be zero-copy from SoA arrays
- Thread safety for OpenMP parallelization

## Usage

This test case is used by the golden record comparison system to validate:

1. **Fortran baseline**: Run with standard SIMPLE executable
2. **Python batch API**: Run with new batch-oriented Python interface
3. **Comparison**: Verify identical results within numerical tolerance