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
    
    simple.load_field('wout.nc')
    particles = simple.sample_surface(100, s=0.9)
    results = simple.trace(particles, tmax=50.0, integrator=simple.SYMPLECTIC)
    
    confined = simple.get_confined(results)
    lost = simple.get_lost(results)
    
    print(f"Total: {particles.shape[1]}")
    print(f"Confined: {confined.shape[1]}")
    print(f"Confinement: {confined.shape[1] / particles.shape[1]:.3f}")

def confinement_study():
    """Study confinement vs flux surface."""
    print("\n=== Confinement Study ===")
    
    simple.load_field('wout.nc')
    
    for s in [0.3, 0.6, 0.9]:
        particles = simple.sample_surface(200, s=s)
        results = simple.trace(particles, tmax=100.0, integrator=simple.MIDPOINT)
        confined = simple.get_confined(results)
        print(f"s={s}: {confined.shape[1]/particles.shape[1]:.3f}")

def volume_example():
    """Volume sampling example."""
    print("\n=== Volume Sampling ===")
    
    simple.load_field('wout.nc')
    particles = simple.sample_volume(150, s_inner=0.1, s_outer=0.8)
    results = simple.trace(particles, tmax=75.0, integrator=simple.RK4)
    
    confined = simple.get_confined(results)
    lost = simple.get_lost(results)
    
    print(f"Volume: {particles.shape[1]}, confined: {confined.shape[1]}")
    if lost['loss_times'].size > 0:
        print(f"Mean loss time: {np.mean(lost['loss_times']):.2f}")

def main():
    """Run all examples."""
    try:
        simple.info()
        basic_example()
        confinement_study() 
        volume_example()
    except FileNotFoundError:
        print("wout.nc not found")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == '__main__':
    main()