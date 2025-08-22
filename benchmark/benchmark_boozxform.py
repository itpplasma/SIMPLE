#!/usr/bin/env python3
"""
Benchmark script to compare orbit integration using VMEC vs BOOZXFORM fields
"""

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import shutil
import time

# Set paths
SIMPLE_EXEC = "../build/simple.x"
VMEC_FILE = "../../benchmark_orbit/booz_xform/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
BOOZ_FILE = "../../benchmark_orbit/booz_xform/boozmn_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"

# Create two input files - one for VMEC, one for BOOZXFORM
simple_in_vmec = """&simple
  vmec_file = 'wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc'
  trace_time = 1.0e-4  ! short integration time
  dtau = 5.0e-7
  ntau = 10000
  ntestpart = 5
  ntimstep = 1000
  sbeg = 0.5
  phibeg = 0.0
  thetabeg = 0.0
  bmod00 = 5.0
  isw_field_type = 1  ! VMEC
/
"""

simple_in_booz = """&simple
  vmec_file = 'wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc'
  boozxform_file = 'boozmn_LandremanPaul2021_QA_reactorScale_lowres_reference.nc'
  trace_time = 1.0e-4  ! short integration time
  dtau = 5.0e-7
  ntau = 10000
  ntestpart = 5
  ntimstep = 1000
  sbeg = 0.5
  phibeg = 0.0
  thetabeg = 0.0
  bmod00 = 5.0
  isw_field_type = 5  ! BOOZXFORM
/
"""

def run_benchmark(case_name, input_content):
    """Run SIMPLE with given input and return results"""
    print(f"\nRunning {case_name} case...")
    
    # Create directory for this case
    case_dir = f"benchmark_{case_name}"
    os.makedirs(case_dir, exist_ok=True)
    os.chdir(case_dir)
    
    # Copy necessary files
    shutil.copy2(f"../{VMEC_FILE}", ".")
    if case_name == "boozxform":
        shutil.copy2(f"../{BOOZ_FILE}", ".")
    
    # Write input file
    with open("simple.in", "w") as f:
        f.write(input_content)
    
    # Run SIMPLE
    start_time = time.time()
    result = subprocess.run([f"../{SIMPLE_EXEC}"], capture_output=True, text=True)
    end_time = time.time()
    
    if result.returncode != 0:
        print(f"Error running {case_name}:")
        print(result.stderr)
        os.chdir("..")
        return None, None
    
    print(f"Execution time: {end_time - start_time:.2f} seconds")
    
    # Read results
    if os.path.exists("confined_fraction.dat"):
        confined_data = np.loadtxt("confined_fraction.dat")
    else:
        print(f"No confined_fraction.dat found for {case_name}")
        confined_data = None
    
    if os.path.exists("times_lost.dat"):
        times_lost = np.loadtxt("times_lost.dat")
    else:
        times_lost = None
    
    os.chdir("..")
    return confined_data, times_lost

def plot_comparison(vmec_data, booz_data):
    """Plot comparison of results"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Plot confined fraction vs time
    if vmec_data[0] is not None and booz_data[0] is not None:
        ax1.plot(vmec_data[0][:, 0], vmec_data[0][:, 1], 'b-', label='VMEC', linewidth=2)
        ax1.plot(booz_data[0][:, 0], booz_data[0][:, 1], 'r--', label='BOOZXFORM', linewidth=2)
        ax1.set_xlabel('Time')
        ax1.set_ylabel('Confined Fraction')
        ax1.set_title('Confinement Comparison: VMEC vs BOOZXFORM')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
    
    # Plot difference
    if vmec_data[0] is not None and booz_data[0] is not None:
        # Interpolate to common time grid if needed
        t_common = vmec_data[0][:, 0]
        vmec_conf = vmec_data[0][:, 1]
        booz_conf = np.interp(t_common, booz_data[0][:, 0], booz_data[0][:, 1])
        
        diff = abs(vmec_conf - booz_conf)
        ax2.semilogy(t_common, diff, 'g-', linewidth=2)
        ax2.set_xlabel('Time')
        ax2.set_ylabel('|Difference|')
        ax2.set_title('Absolute Difference in Confined Fraction')
        ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('benchmark_comparison.png', dpi=150)
    print("\nSaved comparison plot to benchmark_comparison.png")

def main():
    print("=== SIMPLE BOOZXFORM Benchmark ===")
    print(f"Comparing VMEC field evaluation vs BOOZXFORM field evaluation")
    print(f"VMEC file: {VMEC_FILE}")
    print(f"BOOZ file: {BOOZ_FILE}")
    
    # Check if files exist
    if not os.path.exists(SIMPLE_EXEC):
        print(f"Error: SIMPLE executable not found at {SIMPLE_EXEC}")
        print("Please build SIMPLE first with 'make' in the project root")
        return
    
    if not os.path.exists(VMEC_FILE):
        print(f"Error: VMEC file not found at {VMEC_FILE}")
        return
    
    if not os.path.exists(BOOZ_FILE):
        print(f"Error: BOOZXFORM file not found at {BOOZ_FILE}")
        return
    
    # Run benchmarks
    vmec_data = run_benchmark("vmec", simple_in_vmec)
    booz_data = run_benchmark("boozxform", simple_in_booz)
    
    # Compare results
    if vmec_data[0] is not None and booz_data[0] is not None:
        print("\n=== Results Summary ===")
        print(f"VMEC final confined fraction: {vmec_data[0][-1, 1]:.6f}")
        print(f"BOOZ final confined fraction: {booz_data[0][-1, 1]:.6f}")
        print(f"Difference: {abs(vmec_data[0][-1, 1] - booz_data[0][-1, 1]):.6e}")
        
        # Plot comparison
        plot_comparison(vmec_data, booz_data)
    else:
        print("\nBenchmark failed - check error messages above")

if __name__ == "__main__":
    main()