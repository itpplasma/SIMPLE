#!/usr/bin/env python3
"""
Basic Batch Processing Example

Demonstrates fundamental usage of the SIMPLE Python API for batch-oriented
particle orbit tracing with performance validation.

Requirements:
- SIMPLE built with Python support (run `make` in project root)
- VMEC equilibrium file (download example below)
- Python packages: numpy, matplotlib (optional)
"""

import sys
import time
import numpy as np
from pathlib import Path

# Add SIMPLE Python API to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import simple


def download_test_vmec():
    """Download test VMEC file if not present"""
    vmec_file = "wout.nc"
    if not Path(vmec_file).exists():
        print("Downloading test VMEC file...")
        import urllib.request
        url = "https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
        try:
            urllib.request.urlretrieve(url, vmec_file)
            print(f"Downloaded {vmec_file}")
        except Exception as e:
            print(f"Failed to download VMEC file: {e}")
            print("Please download manually or provide your own wout.nc file")
            return None
    return vmec_file


def basic_surface_simulation():
    """Basic simulation with surface sampling"""
    print("=== Basic Surface Simulation ===")
    
    # Create particle batch
    n_particles = 10000
    particles = simple.ParticleBatch(n_particles)
    print(f"Created batch with {particles.n_particles:,} particles")
    
    # Initialize particles on flux surface
    vmec_file = "wout.nc"
    s_surface = 0.9  # Near edge for interesting dynamics
    particles.initialize_surface(vmec_file, s=s_surface)
    
    # Verify initialization
    coords = particles.coordinates
    print(f"Initialization complete:")
    print(f"  s range: [{coords.s.min():.3f}, {coords.s.max():.3f}]")
    print(f"  theta range: [{coords.theta.min():.3f}, {coords.theta.max():.3f}]")
    print(f"  phi range: [{coords.phi.min():.3f}, {coords.phi.max():.3f}]")
    
    # Run simulation
    print("\nRunning simulation...")
    start_time = time.time()
    
    results = simple.trace_orbits(
        particles,
        tmax=1000.0,
        integrator='symplectic_midpoint',
        openmp_threads=4,  # Adjust for your system
        verbose=True
    )
    
    simulation_time = time.time() - start_time
    print(f"Simulation completed in {simulation_time:.2f} seconds")
    print(f"Performance: {n_particles / simulation_time:.0f} particles/second")
    
    # Analyze results
    stats = results.confinement_statistics()
    print(f"\nResults Summary:")
    print(f"  Total particles: {stats.n_total:,}")
    print(f"  Confined: {stats.n_confined:,} ({stats.confined_fraction:.1%})")
    print(f"  Lost: {stats.n_lost:,}")
    if stats.n_lost > 0:
        print(f"  Mean loss time: {stats.mean_loss_time:.3f}")
    
    return results


def validate_soa_performance():
    """Validate Structure of Arrays performance characteristics"""
    print("\n=== SoA Performance Validation ===")
    
    # Create test batch
    particles = simple.ParticleBatch(100000)
    particles.initialize_surface("wout.nc", s=0.9)
    
    # Run performance validation
    from simple.core.batch import validate_soa_performance
    perf_results = validate_soa_performance(particles, verbose=True)
    
    # Check performance characteristics
    if perf_results['cache_friendly']:
        print("✓ Cache-friendly memory layout confirmed")
    else:
        print("⚠ Memory layout may not be optimal")
    
    if perf_results['c_contiguous']:
        print("✓ C-contiguous memory layout confirmed")
    else:
        print("⚠ Memory is not C-contiguous")


def demonstrate_data_access():
    """Demonstrate different data access patterns"""
    print("\n=== Data Access Patterns ===")
    
    # Small batch for demonstration
    particles = simple.ParticleBatch(1000)
    particles.initialize_surface("wout.nc", s=0.9)
    
    # Run quick simulation
    results = simple.trace_orbits(particles, tmax=100.0, verbose=False)
    
    print("1. Raw array access:")
    print(f"   positions shape: {particles.positions.shape}")
    print(f"   loss_times shape: {results.loss_times.shape}")
    print(f"   final_positions shape: {results.final_positions.shape}")
    
    print("\n2. Structured coordinate access:")
    coords = particles.coordinates
    print(f"   s coordinates: {coords.s.shape}, mean = {coords.s.mean():.3f}")
    print(f"   theta coordinates: {coords.theta.shape}, range = [{coords.theta.min():.2f}, {coords.theta.max():.2f}]")
    
    print("\n3. Boolean mask access:")
    confined_mask = results.confined_mask
    lost_mask = results.lost_mask
    print(f"   confined_mask: {confined_mask.sum()} True values")
    print(f"   lost_mask: {lost_mask.sum()} True values")
    
    print("\n4. Single particle access:")
    particle_0 = particles.get_particle(0)
    print(f"   particle 0: s={particle_0[0]:.3f}, theta={particle_0[1]:.3f}")
    
    print("\n5. Particle subsets:")
    confined_final = results.get_confined_particles()
    lost_final, lost_times, lost_trap = results.get_lost_particles()
    print(f"   confined subset: {confined_final.shape}")
    print(f"   lost subset: {lost_final.shape}, {len(lost_times)} loss times")


def compare_integrators():
    """Compare different symplectic integrators"""
    print("\n=== Integrator Comparison ===")
    
    # Create consistent initial conditions
    particles = simple.ParticleBatch(5000)
    particles.initialize_surface("wout.nc", s=0.9)
    
    # Test different integrators
    integrators = [
        'symplectic_euler',
        'symplectic_midpoint', 
        'symplectic_gauss',
        'rk45'
    ]
    
    results_dict = {}
    
    for integrator in integrators:
        print(f"\nTesting {integrator}...")
        
        # Create fresh copy for each test
        test_particles = particles.copy()
        
        start_time = time.time()
        try:
            results = simple.trace_orbits(
                test_particles,
                tmax=500.0,
                integrator=integrator,
                verbose=False
            )
            
            exec_time = time.time() - start_time
            stats = results.confinement_statistics()
            
            results_dict[integrator] = {
                'results': results,
                'time': exec_time,
                'confined_fraction': stats.confined_fraction,
                'mean_loss_time': stats.mean_loss_time if stats.n_lost > 0 else np.inf
            }
            
            print(f"  Execution time: {exec_time:.3f}s")
            print(f"  Confined fraction: {stats.confined_fraction:.4f}")
            if stats.n_lost > 0:
                print(f"  Mean loss time: {stats.mean_loss_time:.3f}")
                
        except Exception as e:
            print(f"  Failed: {e}")
    
    # Summary comparison
    print(f"\nIntegrator Comparison Summary:")
    print(f"{'Integrator':<20} {'Time (s)':<10} {'Confined':<10} {'Mean Loss':<12}")
    print("-" * 55)
    
    for name, data in results_dict.items():
        print(f"{name:<20} {data['time']:<10.3f} {data['confined_fraction']:<10.4f} {data['mean_loss_time']:<12.3f}")


def demonstrate_batch_operations():
    """Demonstrate batch manipulation operations"""
    print("\n=== Batch Operations ===")
    
    # Create original batch
    original = simple.ParticleBatch(10000)
    original.initialize_surface("wout.nc", s=0.9)
    print(f"Original batch: {original.n_particles:,} particles")
    
    # Copy operation
    backup = original.copy()
    print(f"Copied batch: {backup.n_particles:,} particles")
    print(f"Memory shared: {np.shares_memory(original.positions, backup.positions)}")
    
    # Slice operation
    subset = original.slice(0, 1000)
    print(f"Sliced batch: {subset.n_particles:,} particles")
    
    # Verify slice data
    orig_first_1000 = original.positions[:, :1000]
    subset_data = subset.positions
    print(f"Slice data matches: {np.allclose(orig_first_1000, subset_data)}")
    
    # Serialization
    data_dict = original.to_dict()
    restored = simple.ParticleBatch.from_dict(data_dict)
    print(f"Serialization successful: {restored.n_particles:,} particles restored")
    
    # Custom initialization
    custom_data = np.random.rand(5, 1000)
    custom_data[0, :] = 0.8  # Set s = 0.8
    custom_data[1, :] = np.random.uniform(0, 2*np.pi, 1000)  # theta
    custom_data[2, :] = np.random.uniform(0, 2*np.pi, 1000)  # phi
    
    custom_batch = simple.ParticleBatch(1000)
    custom_batch.initialize_from_array(custom_data)
    print(f"Custom initialization: s range [{custom_batch.coordinates.s.min():.3f}, {custom_batch.coordinates.s.max():.3f}]")


def export_results_example():
    """Demonstrate result export capabilities"""
    print("\n=== Results Export ===")
    
    # Run simulation
    particles = simple.ParticleBatch(5000)
    particles.initialize_surface("wout.nc", s=0.9)
    results = simple.trace_orbits(particles, tmax=500.0, verbose=False)
    
    # NumPy export
    numpy_file = "example_results.npz"
    results.save_numpy(numpy_file)
    print(f"Saved to NumPy format: {numpy_file}")
    
    # Load and verify
    loaded_results = simple.BatchResults.load_numpy(numpy_file)
    print(f"Loaded results: {loaded_results.n_particles:,} particles")
    
    # Verify data integrity
    original_stats = results.confinement_statistics()
    loaded_stats = loaded_results.confinement_statistics()
    
    print(f"Data integrity check:")
    print(f"  Original confined fraction: {original_stats.confined_fraction:.6f}")
    print(f"  Loaded confined fraction: {loaded_stats.confined_fraction:.6f}")
    print(f"  Match: {abs(original_stats.confined_fraction - loaded_stats.confined_fraction) < 1e-10}")
    
    # HDF5 export (if h5py available)
    try:
        import h5py
        hdf5_file = "example_results.h5"
        results.save_hdf5(hdf5_file)
        print(f"Saved to HDF5 format: {hdf5_file}")
        
        # Quick HDF5 inspection
        with h5py.File(hdf5_file, 'r') as f:
            print(f"HDF5 groups: {list(f.keys())}")
            if 'results' in f:
                grp = f['results']
                print(f"HDF5 datasets: {list(grp.keys())}")
                print(f"HDF5 attributes: {dict(grp.attrs)}")
                
    except ImportError:
        print("h5py not available - skipping HDF5 export")
    
    # Clean up files
    import os
    for filename in [numpy_file, "example_results.h5"]:
        if os.path.exists(filename):
            os.remove(filename)
            print(f"Cleaned up: {filename}")


def main():
    """Main example execution"""
    print("SIMPLE Python API - Basic Batch Processing Example")
    print("=" * 55)
    
    # Download test data
    vmec_file = download_test_vmec()
    if vmec_file is None:
        print("Cannot proceed without VMEC file")
        return
    
    try:
        # Run examples
        basic_surface_simulation()
        validate_soa_performance()
        demonstrate_data_access()
        compare_integrators()
        demonstrate_batch_operations()
        export_results_example()
        
        print("\n" + "=" * 55)
        print("All examples completed successfully!")
        print("\nNext steps:")
        print("- Try the large_scale_streaming.py example")
        print("- Explore the memory_efficient_processing.py example")
        print("- Read the API reference documentation")
        
    except Exception as e:
        print(f"\nExample failed with error: {e}")
        print("Please check:")
        print("1. SIMPLE was built with `make` (includes Python support)")
        print("2. pysimple module is available in build/ directory")
        print("3. VMEC file wout.nc exists and is valid")
        raise


if __name__ == "__main__":
    main()