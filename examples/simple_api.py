#!/usr/bin/env python3
"""
Example using the clean SIMPLE Python API.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
sys.path.append('../python')
import simple

def main():
    # Download test VMEC file if needed
    vmec_file = "wout.nc"
    if not Path(vmec_file).exists():
        print(f"Download VMEC file: wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O {vmec_file}")
        return

    # Load field from VMEC equilibrium
    print("Loading VMEC field...")
    simple.load_field(vmec_file)
    
    # Sample particles on flux surface s=0.3
    print("Sampling particles...")
    particles = simple.sample_surface(n_particles=100, s=0.3)
    print(f"Sampled {particles.shape[1]} particles, shape: {particles.shape}")
    
    # Trace particle orbits
    print("Tracing orbits...")
    results = simple.trace(particles, tmax=0.1, integrator=simple.MIDPOINT)
    
    # Analyze results  
    confined = simple.get_confined(results)
    lost = simple.get_lost(results)
    
    print(f"Results: {confined.shape[1]} confined, {len(lost['loss_times'])} lost")
    
    # Plot confined fraction vs time
    loss_times = results['loss_times']
    time_points = np.linspace(0, results['tmax'], 100)
    confined_fraction = []
    
    for t in time_points:
        n_confined = np.sum(loss_times >= t)
        confined_fraction.append(n_confined / len(loss_times))
    
    plt.figure(figsize=(8, 6))
    plt.plot(time_points, confined_fraction)
    plt.xlabel('Time [normalized]')
    plt.ylabel('Confined fraction')
    plt.title('Particle confinement vs time')
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    main()