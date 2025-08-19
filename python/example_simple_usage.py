#!/usr/bin/env python3
"""
Simple usage example for SIMPLE Python interface.

Demonstrates the clean functional API.
"""

import simple
import numpy as np

def basic_example():
    """Basic workflow example."""
    print("=== Basic Example ===")
    
    # Load field once for whole simulation
    print("Loading field...")
    simple.load_field('wout.nc')
    
    # Sample particles on flux surface
    print("Sampling particles...")
    particles = simple.sample_surface(n_particles=100, s=0.9)
    print(f"Particles shape: {particles.shape}")
    
    # Trace orbits
    print("Tracing orbits...")
    results = simple.trace(particles, tmax=50.0)
    
    # Analyze results
    confined = simple.get_confined(results)
    lost = simple.get_lost(results)
    
    print(f"Total: {particles.shape[1]}")
    print(f"Confined: {confined.shape[1]}")
    print(f"Lost: {lost['positions'].shape[1]}")
    print(f"Confinement: {confined.shape[1] / particles.shape[1]:.3f}")

def confinement_study():
    """Study confinement vs flux surface."""
    print("\n=== Confinement Study ===")
    
    simple.load_field('wout.nc')
    
    s_values = [0.3, 0.6, 0.9]
    for s in s_values:
        particles = simple.sample_surface(n_particles=200, s=s)
        results = simple.trace(particles, tmax=100.0)
        confined = simple.get_confined(results)
        fraction = confined.shape[1] / particles.shape[1]
        print(f"s={s}: confinement = {fraction:.3f}")

def volume_example():
    """Volume sampling example."""
    print("\n=== Volume Sampling ===")
    
    simple.load_field('wout.nc')
    
    # Sample in volume
    particles = simple.sample_volume(n_particles=150, s_inner=0.1, s_outer=0.8)
    results = simple.trace(particles, tmax=75.0, integrator='symplectic')
    
    confined = simple.get_confined(results)
    lost = simple.get_lost(results)
    
    print(f"Volume particles: {particles.shape[1]}")
    print(f"Confined: {confined.shape[1]}")
    if lost['loss_times'].size > 0:
        print(f"Mean loss time: {np.mean(lost['loss_times']):.2f}")

def main():
    """Run all examples."""
    # Check availability
    try:
        simple.info()
        print()
    except:
        print("SIMPLE not available - check build")
        return
    
    try:
        basic_example()
        confinement_study() 
        volume_example()
    except FileNotFoundError:
        print("Note: wout.nc not found - download example VMEC file")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == '__main__':
    main()