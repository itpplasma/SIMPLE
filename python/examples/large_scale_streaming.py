#!/usr/bin/env python3
"""
Large-Scale Streaming Example

Demonstrates memory-efficient processing of millions of particles using
the SIMPLE Python API streaming capabilities with constant memory usage.

This example shows how to:
- Process millions of particles with limited memory
- Stream results to HDF5 for large dataset handling
- Monitor memory usage and optimize batch sizes
- Analyze streaming results without loading all data

Requirements:
- SIMPLE built with Python support
- h5py for HDF5 streaming (pip install h5py)
- psutil for memory monitoring (pip install psutil)
- At least 4GB RAM recommended for this example
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
            return None
    return vmec_file


def memory_estimation_example():
    """Demonstrate memory estimation and batch size optimization"""
    print("=== Memory Estimation and Optimization ===")
    
    # Estimate memory for different particle counts
    particle_counts = [100_000, 1_000_000, 10_000_000, 100_000_000]
    
    print("Memory estimates:")
    print(f"{'Particles':<12} {'Particles MB':<15} {'Results MB':<12} {'Total GB':<10}")
    print("-" * 55)
    
    for n_particles in particle_counts:
        memory_est = simple.estimate_memory_usage(n_particles, include_results=True)
        print(f"{n_particles:<12,} {memory_est['particles_mb']:<15.1f} "
              f"{memory_est['results_mb']:<12.1f} {memory_est['total_gb']:<10.2f}")
    
    # Optimize batch size for available memory
    available_memory_gb = 4.0  # Assume 4GB available
    total_particles = 5_000_000
    
    optimal_batch = simple.optimize_batch_size(
        target_memory_gb=available_memory_gb,
        n_total=total_particles,
        safety_factor=0.7  # Conservative
    )
    
    print(f"\nBatch size optimization:")
    print(f"  Available memory: {available_memory_gb:.1f} GB")
    print(f"  Total particles: {total_particles:,}")
    print(f"  Optimal batch size: {optimal_batch:,}")
    print(f"  Number of batches: {(total_particles + optimal_batch - 1) // optimal_batch}")
    
    # Verify memory estimate for optimal batch
    optimal_memory = simple.estimate_memory_usage(optimal_batch)
    print(f"  Estimated memory per batch: {optimal_memory['total_gb']:.2f} GB")
    
    return optimal_batch


def streaming_simulation_example():
    """Demonstrate large-scale streaming simulation"""
    print("\n=== Large-Scale Streaming Simulation ===")
    
    # Parameters for streaming simulation
    n_total = 500_000  # Reduced for example (use millions in practice)
    batch_size = 50_000  # Smaller batches for demonstration
    output_file = "streaming_results.h5"
    
    print(f"Streaming simulation parameters:")
    print(f"  Total particles: {n_total:,}")
    print(f"  Batch size: {batch_size:,}")
    print(f"  Number of batches: {(n_total + batch_size - 1) // batch_size}")
    print(f"  Output file: {output_file}")
    
    # Clean up any existing output file
    if Path(output_file).exists():
        Path(output_file).unlink()
        print(f"  Removed existing {output_file}")
    
    # Run streaming simulation with memory monitoring
    print(f"\nStarting streaming simulation...")
    start_time = time.time()
    
    try:
        # Check dependencies
        import h5py
        import psutil
        
        # Run streaming simulation
        stream_results = simple.process_large_simulation(
            vmec_file="wout.nc",
            n_total=n_total,
            tmax=500.0,  # Shorter simulation for example
            batch_size=batch_size,
            output_file=output_file,
            s_surface=0.9,
            integrator='symplectic_midpoint',
            verbose=True,
            memory_limit_gb=8.0
        )
        
        total_time = time.time() - start_time
        
        print(f"\nStreaming simulation completed!")
        print(f"  Total time: {total_time:.1f} seconds")
        print(f"  Processing rate: {n_total / total_time:.0f} particles/second")
        
        return stream_results
        
    except ImportError as e:
        print(f"Missing dependency: {e}")
        print("Please install: pip install h5py psutil")
        return None
    except Exception as e:
        print(f"Streaming simulation failed: {e}")
        return None


def analyze_streaming_results(stream_results):
    """Analyze streaming results without loading all data"""
    print("\n=== Streaming Results Analysis ===")
    
    if stream_results is None:
        print("No streaming results to analyze")
        return
    
    # Get overall summary statistics
    summary = stream_results.get_summary_statistics()
    
    print(f"Overall simulation summary:")
    print(f"  Total particles: {summary['total_particles']:,}")
    print(f"  Total confined: {summary['total_confined']:,}")
    print(f"  Total lost: {summary['total_lost']:,}")
    print(f"  Confined fraction: {summary['confined_fraction']:.4f}")
    
    if summary['total_lost'] > 0:
        print(f"  Mean loss time: {summary['mean_loss_time']:.3f}")
        print(f"  Median loss time: {summary['median_loss_time']:.3f}")
    
    print(f"  Output file: {summary['output_file']}")
    
    # Analyze batch-by-batch statistics
    print(f"\nBatch-by-batch analysis:")
    print(f"{'Batch':<8} {'Particles':<12} {'Confined':<10} {'Lost':<8} {'Fraction':<10}")
    print("-" * 55)
    
    batch_stats = []
    for i, batch_results in enumerate(stream_results.iter_batch_results()):
        stats = batch_results.confinement_statistics()
        batch_stats.append(stats)
        
        print(f"{i:<8} {stats.n_total:<12,} {stats.n_confined:<10,} "
              f"{stats.n_lost:<8,} {stats.confined_fraction:<10.4f}")
        
        # Only show first few batches in example
        if i >= 5:
            print("  ... (showing first 6 batches)")
            break
    
    # Statistical analysis across batches
    if batch_stats:
        confined_fractions = [stats.confined_fraction for stats in batch_stats]
        
        print(f"\nBatch statistics:")
        print(f"  Mean confined fraction: {np.mean(confined_fractions):.4f}")
        print(f"  Std confined fraction: {np.std(confined_fractions):.4f}")
        print(f"  Min confined fraction: {np.min(confined_fractions):.4f}")
        print(f"  Max confined fraction: {np.max(confined_fractions):.4f}")


def memory_efficient_analysis():
    """Demonstrate memory-efficient analysis of large results"""
    print("\n=== Memory-Efficient Analysis ===")
    
    output_file = "streaming_results.h5"
    if not Path(output_file).exists():
        print(f"Results file {output_file} not found - skipping analysis")
        return
    
    try:
        import h5py
        
        # Manual HDF5 analysis for very large datasets
        print("Manual HDF5 analysis (for very large datasets):")
        
        total_particles = 0
        total_confined = 0
        all_loss_times = []
        
        with h5py.File(output_file, 'r') as f:
            print(f"  HDF5 file groups: {list(f.keys())}")
            
            # Process each batch group
            for group_name in f.keys():
                if group_name.startswith('batch_'):
                    grp = f[group_name]
                    
                    # Get batch metadata
                    n_particles = grp.attrs['n_particles']
                    n_confined = grp.attrs['n_confined']
                    n_lost = grp.attrs['n_lost']
                    
                    total_particles += n_particles
                    total_confined += n_confined
                    
                    # Optionally load loss times for detailed analysis
                    if len(all_loss_times) < 10000:  # Limit for example
                        loss_times = grp['loss_times'][:]
                        finite_times = loss_times[loss_times < np.inf]
                        all_loss_times.extend(finite_times[:1000])  # Sample for demo
        
        print(f"  Manual count - Total particles: {total_particles:,}")
        print(f"  Manual count - Total confined: {total_confined:,}")
        print(f"  Manual count - Confined fraction: {total_confined / total_particles:.4f}")
        
        if all_loss_times:
            all_loss_times = np.array(all_loss_times)
            print(f"  Sample loss time statistics:")
            print(f"    Mean: {all_loss_times.mean():.3f}")
            print(f"    Median: {np.median(all_loss_times):.3f}")
            print(f"    Range: [{all_loss_times.min():.3f}, {all_loss_times.max():.3f}]")
        
    except ImportError:
        print("h5py not available - skipping HDF5 analysis")


def streaming_workflow_example():
    """Complete workflow for streaming large-scale simulations"""
    print("\n=== Complete Streaming Workflow ===")
    
    # Step 1: Define simulation parameters
    simulation_params = {
        'vmec_file': 'wout.nc',
        'n_total': 200_000,  # Reduced for example
        'tmax': 300.0,
        's_surface': 0.9,
        'integrator': 'symplectic_midpoint'
    }
    
    # Step 2: Optimize memory usage
    available_memory = 2.0  # GB
    optimal_batch = simple.optimize_batch_size(
        target_memory_gb=available_memory,
        n_total=simulation_params['n_total'],
        safety_factor=0.8
    )
    
    print(f"Workflow parameters:")
    print(f"  Total particles: {simulation_params['n_total']:,}")
    print(f"  Available memory: {available_memory:.1f} GB")
    print(f"  Optimal batch size: {optimal_batch:,}")
    
    # Step 3: Run streaming simulation with memory monitoring
    output_file = "workflow_results.h5"
    
    try:
        with simple.MemoryMonitor("streaming_workflow", verbose=True) as monitor:
            stream_results = simple.process_large_simulation(
                **simulation_params,
                batch_size=optimal_batch,
                output_file=output_file,
                verbose=True
            )
            
            # Monitor peak memory
            peak_memory = monitor.sample()
            
        # Step 4: Analyze results efficiently
        summary = stream_results.get_summary_statistics()
        
        print(f"\nWorkflow results:")
        print(f"  Peak memory usage: {peak_memory:.2f} GB")
        print(f"  Confined fraction: {summary['confined_fraction']:.4f}")
        print(f"  Output file size: {Path(output_file).stat().st_size / (1024**2):.1f} MB")
        
        # Step 5: Selective data loading for visualization
        print(f"\nSelective analysis example:")
        
        # Load only first batch for detailed analysis
        first_batch = stream_results.load_batch_results(0)
        first_stats = first_batch.confinement_statistics()
        
        print(f"  First batch analysis:")
        print(f"    Particles: {first_stats.n_total:,}")
        print(f"    Confined: {first_stats.confined_fraction:.4f}")
        
        # Access specific arrays for plotting
        first_loss_times = first_batch.loss_times
        confined_times = first_loss_times[first_loss_times == np.inf]
        lost_times = first_loss_times[first_loss_times < np.inf]
        
        print(f"    Arrays loaded: {len(confined_times):,} confined, {len(lost_times):,} lost")
        
        # Clean up
        if Path(output_file).exists():
            Path(output_file).unlink()
            print(f"  Cleaned up: {output_file}")
            
    except Exception as e:
        print(f"Workflow failed: {e}")


def performance_scaling_analysis():
    """Analyze performance scaling with particle count and batch size"""
    print("\n=== Performance Scaling Analysis ===")
    
    # Test different batch sizes
    test_particles = [10_000, 25_000, 50_000, 100_000]
    
    print(f"Performance scaling test:")
    print(f"{'Particles':<12} {'Time (s)':<10} {'Rate (p/s)':<12} {'Memory (MB)':<12}")
    print("-" * 50)
    
    for n_particles in test_particles:
        try:
            # Create and run simulation
            particles = simple.ParticleBatch(n_particles)
            particles.initialize_surface("wout.nc", s=0.9)
            
            start_time = time.time()
            results = simple.trace_orbits(particles, tmax=100.0, verbose=False)
            exec_time = time.time() - start_time
            
            # Estimate memory usage
            memory_est = simple.estimate_memory_usage(n_particles)
            
            rate = n_particles / exec_time
            
            print(f"{n_particles:<12,} {exec_time:<10.3f} {rate:<12.0f} {memory_est['total_mb']:<12.1f}")
            
        except Exception as e:
            print(f"{n_particles:<12,} {'FAILED':<10} {'N/A':<12} {'N/A':<12}")
    
    # Benchmark with different thread counts
    print(f"\nOpenMP thread scaling (50k particles):")
    print(f"{'Threads':<10} {'Time (s)':<10} {'Speedup':<10} {'Efficiency':<12}")
    print("-" * 45)
    
    n_particles = 50_000
    baseline_time = None
    
    for n_threads in [1, 2, 4, 8]:
        try:
            particles = simple.ParticleBatch(n_particles)
            particles.initialize_surface("wout.nc", s=0.9)
            
            start_time = time.time()
            results = simple.trace_orbits(
                particles, 
                tmax=100.0, 
                openmp_threads=n_threads,
                verbose=False
            )
            exec_time = time.time() - start_time
            
            if baseline_time is None:
                baseline_time = exec_time
                speedup = 1.0
                efficiency = 1.0
            else:
                speedup = baseline_time / exec_time
                efficiency = speedup / n_threads
            
            print(f"{n_threads:<10} {exec_time:<10.3f} {speedup:<10.2f} {efficiency:<12.3f}")
            
        except Exception as e:
            print(f"{n_threads:<10} {'FAILED':<10} {'N/A':<10} {'N/A':<12}")


def main():
    """Main example execution"""
    print("SIMPLE Python API - Large-Scale Streaming Example")
    print("=" * 55)
    
    # Check dependencies
    try:
        import h5py
        import psutil
        print("✓ Required dependencies available (h5py, psutil)")
    except ImportError as e:
        print(f"✗ Missing dependency: {e}")
        print("Please install: pip install h5py psutil")
        return
    
    # Download test data
    vmec_file = download_test_vmec()
    if vmec_file is None:
        print("Cannot proceed without VMEC file")
        return
    
    try:
        # Run examples
        optimal_batch = memory_estimation_example()
        stream_results = streaming_simulation_example()
        analyze_streaming_results(stream_results)
        memory_efficient_analysis()
        streaming_workflow_example()
        performance_scaling_analysis()
        
        print("\n" + "=" * 55)
        print("Large-scale streaming examples completed successfully!")
        print("\nKey takeaways:")
        print("- Memory usage is constant regardless of total particle count")
        print("- Batch size should be optimized for available memory")
        print("- Streaming enables analysis of arbitrarily large datasets")
        print("- HDF5 format provides efficient storage and access patterns")
        
        # Clean up
        for filename in ["streaming_results.h5", "workflow_results.h5"]:
            if Path(filename).exists():
                Path(filename).unlink()
                print(f"Cleaned up: {filename}")
        
    except Exception as e:
        print(f"\nExample failed with error: {e}")
        print("Please check:")
        print("1. SIMPLE was built with `make` (includes Python support)")
        print("2. Dependencies installed: pip install h5py psutil")
        print("3. Sufficient memory available (at least 2GB recommended)")
        raise


if __name__ == "__main__":
    main()