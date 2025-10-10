# Alpha Particle Confinement Test

## Purpose

This example tests alpha particle confinement in an ITER-size analytical tokamak without TF ripple. It serves as a baseline system test to verify that:

1. The analytical Grad-Shafranov field integrates correctly with geoflux coordinates
2. Meiss canonical coordinates work properly on analytical equilibria
3. Symplectic orbit integration conserves energy and confines particles

## Configuration

**Plasma parameters** (ITER-like):
- Major radius: R0 = 6.2 m
- Minor radius: a = 2.0 m (ε = 0.32)
- Toroidal field: B0 = 5.3 T
- Circular cross-section (κ = 1.0, δ = 0.0)
- No TF ripple (Nripple = 0)

**Particle parameters**:
- Species: Alpha particles (He-4, Z=2, A=4)
- Energy: 3.5 MeV (fusion alpha birth energy)
- Number: 128 particles
- Starting position: s = 0.25 (inner mid-radius flux surface)
- Duration: 0.5 ms (5×10⁻⁴ s)

**Coordinate system**:
- Field: Analytical GS via field-agnostic geoflux
- Orbit coordinates: Meiss canonical (isw_field_type=3)

## Expected Result

**Zero particles lost** - In an axisymmetric tokamak without ripple, all well-confined alpha particles starting at mid-radius remain confined during the 0.5 ms integration window.

Particles should:
- Remain on or near their initial flux surface (s ∈ [0.2, 0.4])
- Exhibit regular passing or trapped orbits
- Conserve energy to numerical precision

## Running the Example

```bash
make run
```

This will:
1. Execute SIMPLE with the analytical field
2. Trace 128 alpha particles for 0.5 ms
3. Output results to `fort.*` files

## Verification

Check the output for:
- `n_lost = 0` (no particles lost)
- Particle trajectories remain near s = 0.3
- Energy conservation within tolerance

## Notes

This test establishes the baseline for comparison with:
- Tokamak with TF ripple (expecting some losses)
- Different starting positions
- Different particle energies
