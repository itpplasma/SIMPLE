#!/usr/bin/env python3
"""
Real unit tests for simple Python API using pytest and TDD.
Tests actual Fortran behavior without defensive programming.
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

# Note: vmec_file fixture is provided by test/conftest.py


@pytest.mark.skipif(not SIMPLE_AVAILABLE, reason="simple module not available")
class TestSimpleConstants:
    """Test integration constants match Fortran."""
    
    def test_integration_constants(self):
        """
        Given: orbit_symplectic_base.f90 constants
        When: Importing simple module
        Then: Constants match Fortran values exactly
        """
        # From orbit_symplectic_base.f90:
        # integer, parameter :: RK45 = 0, EXPL_IMPL_EULER = 1, IMPL_EXPL_EULER = 2, &
        #   MIDPOINT = 3, GAUSS1 = 4, GAUSS2 = 5, GAUSS3 = 6, GAUSS4 = 7, LOBATTO3 = 15
        
        assert simple.RK45 == 0
        assert simple.EXPL_IMPL_EULER == 1
        assert simple.IMPL_EXPL_EULER == 2
        assert simple.MIDPOINT == 3
        assert simple.GAUSS1 == 4
        assert simple.GAUSS2 == 5
        assert simple.GAUSS3 == 6
        assert simple.GAUSS4 == 7
        assert simple.LOBATTO3 == 15
    
    def test_default_values(self):
        """
        Given: params.f90 default values
        When: Importing simple module  
        Then: Defaults match Fortran
        """
        # From params.f90: sbeg=0.5d0, trace_time=1d-1, integmode = EXPL_IMPL_EULER
        assert simple.DEFAULT_SURFACE == 0.5
        assert simple.DEFAULT_TMAX == 0.1
        assert simple.DEFAULT_INTEGRATOR == simple.EXPL_IMPL_EULER


@pytest.mark.skipif(not SIMPLE_AVAILABLE, reason="simple module not available")  
class TestFieldLoading:
    """Test field loading behavior."""
    
    def test_load_field_basic(self, vmec_file):
        """
        Given: Valid VMEC file
        When: Loading field
        Then: No exceptions raised
        """
        simple.load_field(vmec_file)
        # If we get here, loading succeeded

    def test_load_field_double_init(self, vmec_file):
        """
        Given: Field already loaded  
        When: Loading field again
        Then: Reinitializes without crash (TDD: test actual behavior)
        """
        simple.load_field(vmec_file)
        # Load again - should reinitialize, not crash
        simple.load_field(vmec_file)


@pytest.mark.skipif(not SIMPLE_AVAILABLE, reason="simple module not available")
class TestSampling:
    """Test particle sampling functions."""
    
    def test_sample_surface_basic(self, vmec_file):
        """
        Given: Loaded field
        When: Sampling particles on surface
        Then: Returns proper (5, n_particles) array
        """
        simple.load_field(vmec_file)
        
        n_particles = 10
        particles = simple.sample_surface(n_particles, s=0.8)
        
        assert particles.shape == (5, n_particles)
        assert isinstance(particles, np.ndarray)
        # Check that we got actual particle data (not zeros)
        assert not np.allclose(particles, 0)
    
    def test_sample_volume_basic(self, vmec_file):
        """
        Given: Loaded field
        When: Sampling particles in volume
        Then: Returns proper (5, n_particles) array  
        """
        simple.load_field(vmec_file)
        
        n_particles = 5
        particles = simple.sample_volume(n_particles, s_inner=0.2, s_outer=0.8)
        
        assert particles.shape == (5, n_particles)
        assert isinstance(particles, np.ndarray)
        assert not np.allclose(particles, 0)
    
    def test_sample_different_sizes(self, vmec_file):
        """
        Given: Loaded field
        When: Sampling different numbers of particles
        Then: Arrays have correct sizes
        """
        simple.load_field(vmec_file)
        
        for n in [1, 5, 20, 100]:
            particles = simple.sample_surface(n)
            assert particles.shape == (5, n)


@pytest.mark.skipif(not SIMPLE_AVAILABLE, reason="simple module not available")
class TestTracing:
    """Test orbit tracing."""
    
    def test_trace_basic(self, vmec_file):
        """
        Given: Particles and loaded field
        When: Tracing orbits with short time
        Then: Returns results dictionary
        """
        simple.load_field(vmec_file)
        
        # Use short integration time for tests
        particles = simple.sample_surface(5, s=0.8)
        results = simple.trace(particles, tmax=1e-4, integrator=simple.EXPL_IMPL_EULER)
        
        assert isinstance(results, dict)
        assert 'final_positions' in results
        assert 'loss_times' in results
        assert 'tmax' in results
        
        # Check shapes
        assert results['final_positions'].shape == particles.shape
        assert results['loss_times'].shape == (particles.shape[1],)
        assert results['tmax'] == 1e-4
    
    def test_trace_different_integrators(self, vmec_file):
        """
        Given: Particles and loaded field
        When: Using different integrators
        Then: All complete without crash
        """
        simple.load_field(vmec_file)
        particles = simple.sample_surface(3, s=0.8)
        
        integrators = [simple.RK45, simple.EXPL_IMPL_EULER, simple.MIDPOINT]
        
        for integrator in integrators:
            results = simple.trace(particles, tmax=1e-4, integrator=integrator)
            assert isinstance(results, dict)
    
    def test_trace_array_transpose(self, vmec_file):
        """
        Given: Particles in wrong orientation 
        When: Tracing orbits
        Then: Automatically transposes and works
        """
        simple.load_field(vmec_file)
        
        # Create particles in (n_particles, 5) format instead of (5, n_particles)
        particles_wrong = np.random.rand(3, 5) * 0.1 + [0.5, 0, 0, 1, 0.5]
        
        # Should automatically transpose
        results = simple.trace(particles_wrong, tmax=1e-4)
        assert isinstance(results, dict)


@pytest.mark.skipif(not SIMPLE_AVAILABLE, reason="simple module not available")
class TestConfinementAnalysis:
    """Test confinement analysis functions."""
    
    def test_get_confined_and_lost(self, vmec_file):
        """
        Given: Simulation results
        When: Analyzing confinement
        Then: Can separate confined and lost particles
        """
        simple.load_field(vmec_file)
        
        particles = simple.sample_surface(10, s=0.9)
        results = simple.trace(particles, tmax=1e-4)
        
        confined = simple.get_confined(results)
        lost = simple.get_lost(results)
        
        assert isinstance(confined, np.ndarray)
        assert isinstance(lost, dict)
        assert 'positions' in lost
        assert 'loss_times' in lost
        
        # Check that total particles = confined + lost
        n_confined = confined.shape[1] if confined.size > 0 else 0
        n_lost = lost['positions'].shape[1] if lost['positions'].size > 0 else 0
        assert n_confined + n_lost == particles.shape[1]


@pytest.mark.skipif(not SIMPLE_AVAILABLE, reason="simple module not available")
class TestParameterAccess:
    """Test parameter get/set functions."""
    
    def test_set_get_parameters(self):
        """
        Given: Parameter functions
        When: Setting and getting parameters
        Then: Values are correctly stored and retrieved
        """
        # Test setting parameters
        simple.set_parameters(ntimstep=5000, relerr=1e-12)
        
        # Test getting specific parameters
        ntimstep = simple.get_parameters('ntimstep')
        assert ntimstep == 5000
        
        # Test getting multiple parameters
        params = simple.get_parameters('ntimstep', 'relerr')
        assert params['ntimstep'] == 5000
        assert params['relerr'] == 1e-12
        
        # Test getting all common parameters
        all_params = simple.get_parameters()
        assert isinstance(all_params, dict)
        assert 'ntimstep' in all_params


@pytest.mark.skipif(not SIMPLE_AVAILABLE, reason="simple module not available")
class TestRobustness:
    """Test robustness and edge cases."""
    
    def test_multiple_reinitializations(self, vmec_file):
        """
        Given: Loaded field
        When: Multiple reinitializations
        Then: No crashes, last initialization takes effect
        """
        # Test multiple calls to load_field
        for _ in range(3):
            simple.load_field(vmec_file)
        
        # Should still work
        particles = simple.sample_surface(5)
        results = simple.trace(particles, tmax=1e-4)
        assert isinstance(results, dict)
    
    def test_parameter_changes_between_runs(self, vmec_file):
        """
        Given: Loaded field
        When: Changing parameters between runs
        Then: New parameters take effect
        """
        simple.load_field(vmec_file)
        particles = simple.sample_surface(5)
        
        # First run with default integrator
        results1 = simple.trace(particles, tmax=1e-4, integrator=simple.EXPL_IMPL_EULER)
        
        # Second run with different integrator  
        results2 = simple.trace(particles, tmax=1e-4, integrator=simple.MIDPOINT)
        
        # Both should succeed
        assert isinstance(results1, dict)
        assert isinstance(results2, dict)


if __name__ == "__main__":
    # Run tests directly for development
    pytest.main([__file__, "-v", "--tb=short"])