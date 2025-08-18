#!/usr/bin/env python3
"""
Basic SIMPLE Python API Usage Example

Demonstrates the essential interface to existing Fortran functionality:
- Loading particles using existing samplers.f90 functions
- Running simulations via existing simple_main.f90 execution
- Accessing results through zero-copy Fortran array wrappers

This is a minimal example showing the pure interface layer.
"""

import sys
from pathlib import Path

# Add SIMPLE Python API to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import simple


def download_test_vmec():
    """Download test VMEC file if needed"""
    vmec_file = "wout.nc"
    if not Path(vmec_file).exists():
        print("Downloading test VMEC file...")
        import urllib.request
        url = "https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
        try:
            urllib.request.urlretrieve(url, vmec_file)
            print(f"Downloaded {vmec_file}")
        except Exception as e:
            print(f"Failed to download: {e}")
            return None
    return vmec_file


def basic_simulation_example():
    """Basic simulation using existing Fortran functionality"""
    print("=== Basic SIMPLE Python API Example ===")
    
    # 1. Initialize particles using existing samplers.f90
    print("1. Initializing particles using existing surface sampling...")
    vmec_file = "wout.nc"
    
    # Use existing SurfaceSampler (wraps samplers.f90)
    sampler = simple.SurfaceSampler(vmec_file)
    particle_data = sampler.sample_surface_fieldline(n_particles=1000)
    
    # Create batch wrapper around Fortran arrays
    particles = simple.ParticleBatch.from_fortran_arrays(particle_data)
    print(f"   Initialized {particles.n_particles} particles")
    print(f"   s range: [{particles.coordinates.s.min():.3f}, {particles.coordinates.s.max():.3f}]")
    
    # 2. Run simulation using existing simple_main.f90
    print("2. Running simulation via existing Fortran execution...")
    results = simple.trace_orbits(
        particles,
        tmax=500.0,
        integrator='symplectic_midpoint',
        verbose=True
    )
    
    # 3. Access results through zero-copy wrappers
    print("3. Accessing results through Fortran array wrappers...")
    stats = results.confinement_statistics()
    
    print(f"Results:")
    print(f"   Total particles: {stats.n_total}")
    print(f"   Confined: {stats.n_confined} ({stats.confined_fraction:.1%})")
    print(f"   Lost: {stats.n_lost}")
    if stats.n_lost > 0:
        print(f"   Mean loss time: {stats.mean_loss_time:.3f}")
    
    return results


def volume_sampling_example():
    """Example using existing volume sampling"""
    print("\n=== Volume Sampling Example ===")
    
    # Use existing VolumeSampler (wraps samplers.f90)
    vmec_file = "wout.nc"
    sampler = simple.VolumeSampler(vmec_file)
    
    # Sample between flux surfaces using existing sample_volume_single
    particle_data = sampler.sample_volume_single(
        n_particles=500, 
        s_inner=0.1, 
        s_outer=0.9
    )
    
    particles = simple.ParticleBatch.from_fortran_arrays(particle_data)
    print(f"Volume sampled {particles.n_particles} particles")
    print(f"s range: [{particles.coordinates.s.min():.3f}, {particles.coordinates.s.max():.3f}]")
    
    # Quick simulation
    results = simple.trace_orbits(particles, tmax=200.0, verbose=False)
    stats = results.confinement_statistics()
    print(f"Confined fraction: {stats.confined_fraction:.3f}")
    
    return results


def file_loading_example():
    """Example using existing file loading"""
    print("\n=== File Loading Example ===")
    
    # First save some particles to demonstrate file loading
    vmec_file = "wout.nc"
    sampler = simple.SurfaceSampler(vmec_file)
    particle_data = sampler.sample_surface_fieldline(n_particles=100)
    
    # Save to file (uses existing save_starting_points)
    import numpy as np
    with open("example_particles.dat", "w") as f:
        for i in range(particle_data.shape[1]):
            coords = particle_data[:, i]
            f.write(f"{coords[0]:.6e} {coords[1]:.6e} {coords[2]:.6e} {coords[3]:.6e} {coords[4]:.6e}\n")
    
    # Load using existing file sampler (wraps sample_read)
    file_sampler = simple.FileSampler(vmec_file)
    loaded_data = file_sampler.load_from_file("example_particles.dat")
    
    particles = simple.ParticleBatch.from_fortran_arrays(loaded_data)
    print(f"Loaded {particles.n_particles} particles from file")
    
    # Verify loaded data matches
    print(f"Data integrity: {np.allclose(particle_data, loaded_data, rtol=1e-10)}")
    
    # Cleanup
    Path("example_particles.dat").unlink()
    
    return particles


def integrator_comparison():
    """Compare different integrators using existing methods"""
    print("\n=== Integrator Comparison ===")
    
    # Create consistent initial conditions
    vmec_file = "wout.nc"
    sampler = simple.SurfaceSampler(vmec_file)
    particle_data = sampler.sample_surface_fieldline(n_particles=500)
    
    integrators = ['symplectic_euler', 'symplectic_midpoint', 'symplectic_gauss']
    
    for integrator in integrators:
        # Fresh copy for each test
        particles = simple.ParticleBatch.from_fortran_arrays(particle_data.copy())
        
        # Run using existing integrator implementation
        results = simple.trace_orbits(
            particles,
            tmax=300.0,
            integrator=integrator,
            verbose=False
        )
        
        stats = results.confinement_statistics()
        print(f"{integrator:20}: {stats.confined_fraction:.4f} confined")


def main():
    """Main example execution"""
    print("SIMPLE Python API - Pure Interface Examples")
    print("All functionality delegates to existing Fortran implementations")
    print("=" * 60)
    
    # Download test data
    vmec_file = download_test_vmec()
    if vmec_file is None:
        print("Cannot proceed without VMEC file")
        return
    
    try:
        # Run examples
        basic_simulation_example()
        volume_sampling_example()
        file_loading_example()
        integrator_comparison()
        
        print("\n" + "=" * 60)
        print("All examples completed successfully!")
        print("The Python API provides zero-copy access to existing Fortran functionality.")
        
    except Exception as e:
        print(f"\nExample failed with error: {e}")
        print("Please check:")
        print("1. SIMPLE was built with `make` (includes Python support)")
        print("2. pysimple module is available in build/ directory")
        print("3. VMEC file exists and is valid")
        raise


if __name__ == "__main__":
    main()