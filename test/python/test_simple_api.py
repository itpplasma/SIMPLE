#!/usr/bin/env python3
"""
Test suite for simple Python API (Issue #90)

Tests the simple functional interface:
- simple.load_field()
- simple.sample_surface() / simple.sample_volume()
- simple.trace()
- simple.get_confined() / simple.get_lost()
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add python directory to path for simple import
python_dir = Path(__file__).parent.parent.parent / "python"
sys.path.insert(0, str(python_dir))

try:
    import simple
    SIMPLE_AVAILABLE = True
except ImportError as e:
    SIMPLE_AVAILABLE = False
    print(f"Warning: simple module not available: {e}")


@pytest.mark.skipif(not SIMPLE_AVAILABLE, reason="simple module not available")
class TestSimpleAPI:
    """Test suite for simple functional API"""
    
    def setup_method(self):
        """Setup test environment"""
        self.vmec_file = Path(__file__).parent / "wout.nc"
        if not self.vmec_file.exists():
            pytest.skip("VMEC test file not available")
    
    def test_load_field(self):
        """
        Given: VMEC equilibrium file
        When: Loading field with simple.load_field()
        Then: Field loads without error
        """
        # This will call pysimple.init_field() internally
        simple.load_field(str(self.vmec_file))
        
        # If we get here without exception, field loading succeeded
        assert True
    
    def test_sample_surface(self):
        """
        Given: Loaded field
        When: Sampling particles on surface
        Then: Returns proper SoA array format
        """
        simple.load_field(str(self.vmec_file))
        
        n_particles = 100
        particles = simple.sample_surface(n_particles=n_particles, s=0.9)
        
        assert particles.shape == (5, n_particles), "Must return (5, n_particles) SoA format"
        assert isinstance(particles, np.ndarray), "Must return numpy array"
    
    def test_sample_volume(self):
        """
        Given: Loaded field
        When: Sampling particles in volume
        Then: Returns proper SoA array format
        """
        simple.load_field(str(self.vmec_file))
        
        n_particles = 50
        particles = simple.sample_volume(n_particles=n_particles, s_inner=0.1, s_outer=0.9)
        
        assert particles.shape == (5, n_particles), "Must return (5, n_particles) SoA format"
        assert isinstance(particles, np.ndarray), "Must return numpy array"
    
    def test_trace_particles(self):
        """
        Given: Particles and loaded field
        When: Tracing with simple.trace()
        Then: Returns results dictionary
        """
        simple.load_field(str(self.vmec_file))
        
        # Sample some particles
        particles = simple.sample_surface(n_particles=10, s=0.9)
        
        # Trace them
        results = simple.trace(particles, tmax=10.0)
        
        assert isinstance(results, dict), "Must return results dictionary"
        # Note: Results content depends on actual Fortran integration
    
    def test_get_confined_and_lost(self):
        """
        Given: Simulation results
        When: Analyzing confinement
        Then: Can separate confined and lost particles
        """
        simple.load_field(str(self.vmec_file))
        
        particles = simple.sample_surface(n_particles=20, s=0.9)
        results = simple.trace(particles, tmax=50.0)
        
        # Analyze results
        confined = simple.get_confined(results)
        lost = simple.get_lost(results)
        
        assert isinstance(confined, np.ndarray), "get_confined must return array"
        assert isinstance(lost, dict), "get_lost must return dict with positions and times"
    
    def test_load_particles_from_file(self):
        """
        Given: Particle file and loaded field
        When: Loading particles from file
        Then: Returns proper SoA array format
        """
        # Create test particle file
        test_file = Path(__file__).parent / "test_particles.dat"
        test_particles = np.random.randn(10, 5)  # 10 particles, 5 coordinates each
        np.savetxt(test_file, test_particles)
        
        try:
            simple.load_field(str(self.vmec_file))
            particles = simple.load_particles(str(test_file))
            
            assert particles.shape[0] == 5, "Must return (5, n_particles) SoA format"
            assert particles.shape[1] == 10, "Must load correct number of particles"
        finally:
            # Clean up test file
            if test_file.exists():
                test_file.unlink()


@pytest.mark.skipif(not SIMPLE_AVAILABLE, reason="simple module not available")
class TestSimpleAPIWorkflow:
    """Test complete workflow with simple API"""
    
    def setup_method(self):
        """Setup test environment"""
        self.vmec_file = Path(__file__).parent / "wout.nc"
        if not self.vmec_file.exists():
            pytest.skip("VMEC test file not available")
    
    def test_complete_workflow(self):
        """
        Given: Simple API
        When: Running complete simulation workflow
        Then: All steps complete without error
        """
        # Load field
        simple.load_field(str(self.vmec_file))
        
        # Sample particles
        particles = simple.sample_surface(n_particles=50, s=0.8)
        
        # Trace orbits
        results = simple.trace(particles, tmax=20.0)
        
        # Analyze results
        confined = simple.get_confined(results)
        lost = simple.get_lost(results)
        
        # Basic sanity checks
        assert particles.shape[0] == 5, "Initial particles in SoA format"
        assert isinstance(results, dict), "Results is dictionary"
        assert isinstance(confined, np.ndarray), "Confined particles array"
        assert isinstance(lost, dict), "Lost particles dictionary"
        
        print(f"Workflow completed: {particles.shape[1]} initial particles")
        print(f"Confined: {confined.shape[1] if confined.size > 0 else 0}")
        print(f"Lost: {lost['positions'].shape[1] if 'positions' in lost and lost['positions'].size > 0 else 0}")


if __name__ == "__main__":
    # Run tests directly for development
    pytest.main([__file__, "-v", "--tb=short"])