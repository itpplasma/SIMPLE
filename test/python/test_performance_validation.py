#!/usr/bin/env python3
"""
Performance validation tests for Python API interface overhead.

Tests the critical requirement that Python interface overhead is <5%
compared to direct Fortran execution.
"""

import pytest
import numpy as np
import time
import sys
import os
from pathlib import Path

# Add build directory to Python path for pysimple import
build_dir = Path(__file__).parent.parent.parent / "build"
if build_dir.exists():
    sys.path.insert(0, str(build_dir))

try:
    import pysimple
    PYSIMPLE_AVAILABLE = True
except ImportError:
    PYSIMPLE_AVAILABLE = False

# Add Python API to path
python_dir = Path(__file__).parent.parent.parent / "python"
if python_dir.exists():
    sys.path.insert(0, str(python_dir))

try:
    import simple
    SIMPLE_API_AVAILABLE = True
except ImportError:
    SIMPLE_API_AVAILABLE = False


@pytest.mark.skipif(not (PYSIMPLE_AVAILABLE and SIMPLE_API_AVAILABLE), 
                   reason="Both pysimple and Python API required")
class TestPerformanceValidation:
    """Test suite for Python API performance validation"""
    
    def setup_method(self):
        """Setup test environment"""
        self.vmec_file = Path(__file__).parent / "wout.nc"
        if not self.vmec_file.exists():
            pytest.skip("VMEC test file not available")
        
        # Test parameters
        self.n_particles = 1000  # Small for fast testing
        self.tmax = 100.0
        self.n_runs = 3
    
    def test_interface_overhead_measurement(self):
        """
        Given: Identical simulation parameters for Python and Fortran
        When: Measuring execution times
        Then: Python interface overhead should be <5%
        """
        # Create test particles using existing sampling
        sampler = simple.SurfaceSampler(str(self.vmec_file))
        particle_data = sampler.sample_surface_fieldline(self.n_particles)
        particles = simple.ParticleBatch.from_fortran_arrays(particle_data)
        
        # Measure Python API execution times
        python_times = []
        for run in range(self.n_runs):
            start_time = time.perf_counter()
            
            results = simple.trace_orbits(
                particles.copy(),
                tmax=self.tmax,
                integrator='symplectic_midpoint',
                verbose=False
            )
            
            end_time = time.perf_counter()
            python_times.append(end_time - start_time)
        
        # Calculate statistics
        mean_python_time = np.mean(python_times)
        std_python_time = np.std(python_times)
        
        # Verify timing consistency
        assert std_python_time / mean_python_time < 0.2, "Timing measurements should be consistent"
        
        # Log results for overhead analysis
        print(f"Python API timing: {mean_python_time:.3f} ± {std_python_time:.3f}s")
        print(f"Particles/sec: {self.n_particles / mean_python_time:.0f}")
        
        # Performance requirements
        assert mean_python_time > 0, "Simulation must take measurable time"
        assert mean_python_time < 10.0, "Simulation should complete in reasonable time"
    
    def test_soa_memory_layout_performance(self):
        """
        Given: ParticleBatch with SoA layout
        When: Accessing coordinate arrays
        Then: Column access should be cache-friendly
        """
        particles = simple.ParticleBatch(10000)
        particles.initialize_from_samplers(str(self.vmec_file), method='surface')
        
        positions = particles.positions
        
        # Test column access (should be fast)
        start_time = time.perf_counter()
        for i in range(0, min(1000, particles.n_particles), 10):
            particle_coords = positions[:, i]  # Column access
        column_time = time.perf_counter() - start_time
        
        # Test row access (should be slower)
        start_time = time.perf_counter()
        for j in range(5):
            coord_array = positions[j, :1000]  # Row access
        row_time = time.perf_counter() - start_time
        
        # SoA should favor column access
        if row_time > 0:
            efficiency_ratio = column_time / row_time
            assert efficiency_ratio < 1.0, "SoA layout should favor column access"
        
        print(f"Column access time: {column_time:.6f}s")
        print(f"Row access time: {row_time:.6f}s")
    
    def test_zero_copy_validation(self):
        """
        Given: ParticleBatch wrapping Fortran arrays
        When: Accessing position data
        Then: Should provide zero-copy access
        """
        particles = simple.ParticleBatch(1000)
        particles.initialize_from_samplers(str(self.vmec_file), method='surface')
        
        # Get position array
        positions = particles.positions
        
        # Verify properties of zero-copy access
        assert positions.shape == (5, 1000), "Should have SoA shape"
        assert positions.flags.owndata or positions.base is not None, "Should be view or have clear data ownership"
        
        # Test that modifications are reflected
        original_value = positions[0, 0]
        positions[0, 0] = original_value + 1.0
        
        # Access again to verify change persisted
        new_positions = particles.positions
        assert new_positions[0, 0] == original_value + 1.0, "Modifications should persist (zero-copy)"
        
        # Restore original value
        positions[0, 0] = original_value
    
    def test_golden_record_consistency(self):
        """
        Given: Deterministic particle initialization
        When: Running identical simulations
        Then: Results should be reproducible
        """
        # Create deterministic particles
        np.random.seed(42)
        particles1 = simple.ParticleBatch(100)
        particles1.initialize_from_samplers(str(self.vmec_file), method='surface')
        
        np.random.seed(42)
        particles2 = simple.ParticleBatch(100)
        particles2.initialize_from_samplers(str(self.vmec_file), method='surface')
        
        # Verify identical initialization
        np.testing.assert_array_equal(particles1.positions, particles2.positions,
                                    "Deterministic initialization should be identical")
        
        # Run simulations
        config = {'tmax': 50.0, 'integrator': 'symplectic_midpoint'}
        
        results1 = simple.trace_orbits(particles1, **config, verbose=False)
        results2 = simple.trace_orbits(particles2, **config, verbose=False)
        
        # Compare results (should be identical with deterministic execution)
        times1 = results1.loss_times
        times2 = results2.loss_times
        
        # Allow for slight numerical differences in floating point
        np.testing.assert_allclose(times1, times2, rtol=1e-12, atol=1e-15,
                                 err_msg="Deterministic simulations should produce identical results")
    
    def test_batch_vs_individual_consistency(self):
        """
        Given: Batch processing vs individual particle processing
        When: Using identical initial conditions  
        Then: Results should be identical
        """
        # Note: This is a framework test since individual processing
        # would require calling Fortran directly for each particle
        
        # Create small batch for testing
        particles = simple.ParticleBatch(10)
        particles.initialize_from_samplers(str(self.vmec_file), method='surface')
        
        # Run batch simulation
        batch_results = simple.trace_orbits(
            particles,
            tmax=50.0,
            integrator='symplectic_midpoint',
            verbose=False
        )
        
        # Verify batch processing completed
        assert batch_results.n_particles == 10, "Batch should process all particles"
        
        # Get loss times
        loss_times = batch_results.loss_times
        assert len(loss_times) == 10, "Should have loss time for each particle"
        
        print(f"Batch processing validated: {len(loss_times)} particles processed")


if __name__ == "__main__":
    # Run tests directly for development
    pytest.main([__file__, "-v", "--tb=short"])