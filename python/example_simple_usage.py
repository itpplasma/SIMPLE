#!/usr/bin/env python3
"""Simple usage example for SIMPLE Python interface."""

import simple

def main():
    simple.load_field('wout.nc')
    
    # Basic example
    particles = simple.sample_surface(100, s=0.3)
    results = simple.trace(particles, tmax=0.1)
    confined = simple.get_confined(results)
    print(f"Confined: {confined.shape[1]}/{particles.shape[1]}")
    
    # Compare surfaces
    for s in [0.3, 0.6, 0.9]:
        particles = simple.sample_surface(50, s=s)
        results = simple.trace(particles, tmax=0.1)
        confined = simple.get_confined(results)
        print(f"s={s}: {confined.shape[1]/particles.shape[1]:.3f}")

if __name__ == '__main__':
    main()