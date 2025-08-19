#!/usr/bin/env python3
"""
Simple usage example for SIMPLE Python interface.

Direct functional API - no OOP complexity.
"""

import simple
import numpy as np

def main():
    """Demonstrate simple SIMPLE usage."""
    
    # Check if SIMPLE is available
    simple.info()
    
    # Load field once for whole simulation
    print("Loading field from VMEC equilibrium...")
    simple.load_field('wout.nc')
    
    # Sample particles on flux surface
    print("Sampling 1000 particles on s=0.9 surface...")
    particles = simple.sample_surface(n_particles=1000, s=0.9)
    print(f"Sampled particles shape: {particles.shape}")
    
    # Trace orbits
    print("Tracing particle orbits...")
    results = simple.trace(particles, tmax=100.0)
    
    # Analyze results
    confined = simple.get_confined(results)
    lost = simple.get_lost(results)
    
    print(f"Total particles: {particles.shape[1]}")
    print(f"Confined particles: {confined.shape[1]}")
    print(f"Lost particles: {lost['positions'].shape[1]}")
    print(f"Confinement fraction: {confined.shape[1] / particles.shape[1]:.3f}")
    
    # Volume sampling
    print("\nVolume sampling example...")
    volume_particles = simple.sample_volume(n_particles=500)
    print(f"Volume particles shape: {volume_particles.shape}")


if __name__ == '__main__':
    main()