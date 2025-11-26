#!/usr/bin/env python3
"""Compare guiding-center and Pauli particle orbit models.

This script runs SIMPLE with different orbit_model settings (GC vs Pauli)
and plots the trajectories in R-Z and R-phi projections for comparison.

Usage:
    python compare_gc_pauli.py [working_directory]

The working directory must contain:
    - wout.nc: VMEC equilibrium file
    - coils file: Coils file for Biot-Savart field
    - simple.in: Base configuration file

Output:
    - gc_pauli_comparison.png: R-Z and R-phi trajectory comparison
"""

import os
import sys
import subprocess
import tempfile
import shutil
import urllib.request
import numpy as np

# Try to import matplotlib; skip plotting if not available
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Warning: matplotlib not available, skipping plot generation")


def download_test_data(dest_dir):
    """Download VMEC test file if not present."""
    wout_file = os.path.join(dest_dir, 'wout.nc')
    if not os.path.exists(wout_file):
        url = ("https://github.com/hiddenSymmetries/simsopt/raw/master/"
               "tests/test_files/wout_LandremanPaul2021_QA_reactorScale_"
               "lowres_reference.nc")
        print(f"Downloading VMEC test file to {wout_file}...")
        urllib.request.urlretrieve(url, wout_file)
        print("Download complete.")
    return wout_file


def run_simple(work_dir, orbit_model, trace_time=1e-6, ntimstep=100):
    """Run SIMPLE with specified orbit model and return trajectory."""
    # Create config file
    config = f"""&config
trace_time = {trace_time}
sbeg = 0.6d0
ntestpart = 1
netcdffile = 'wout.nc'
isw_field_type = 1
integmode = 3
npoiper2 = 128
deterministic = .True.
field_input = 'coils.5C'
orbit_model = {orbit_model}
ntimstep = {ntimstep}
/
"""

    config_path = os.path.join(work_dir, f'simple_orbit{orbit_model}.in')
    with open(config_path, 'w') as f:
        f.write(config)

    # Run SIMPLE
    simple_exe = os.path.join(os.path.dirname(__file__), '..', 'build', 'simple.x')
    simple_exe = os.path.abspath(simple_exe)

    print(f"Running SIMPLE with orbit_model={orbit_model}...")
    result = subprocess.run(
        [simple_exe, config_path],
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=300
    )

    if result.returncode != 0:
        print(f"SIMPLE failed with orbit_model={orbit_model}")
        print("STDOUT:", result.stdout[-1000:] if len(result.stdout) > 1000 else result.stdout)
        print("STDERR:", result.stderr[-1000:] if len(result.stderr) > 1000 else result.stderr)
        return None

    # Read orbit output
    try:
        import netCDF4 as nc
        orbits_path = os.path.join(work_dir, 'orbits.nc')
        if os.path.exists(orbits_path):
            with nc.Dataset(orbits_path, 'r') as ds:
                # Get trajectory for first particle
                s = ds.variables['s'][0, :]
                theta = ds.variables['theta'][0, :]
                phi = ds.variables['phi'][0, :]
                time = ds.variables['time'][0, :]

                # Convert to cylindrical for plotting
                # (would need VMEC splines, so use s,theta,phi directly)
                return {
                    's': s,
                    'theta': theta,
                    'phi': phi,
                    'time': time
                }
    except ImportError:
        print("netCDF4 not available, cannot read orbit output")
    except Exception as e:
        print(f"Error reading orbit output: {e}")

    return None


def plot_comparison(gc_traj, pauli_traj, output_file):
    """Generate comparison plot of GC and Pauli trajectories."""
    if not HAS_MATPLOTLIB:
        print("Skipping plot generation (matplotlib not available)")
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # s vs phi (similar to R-phi projection in flux coordinates)
    ax = axes[0, 0]
    if gc_traj is not None:
        ax.plot(gc_traj['phi'], gc_traj['s'], 'b-', label='Guiding-center', alpha=0.7)
    if pauli_traj is not None:
        ax.plot(pauli_traj['phi'], pauli_traj['s'], 'r-', label='Pauli particle', alpha=0.7)
    ax.set_xlabel('phi (rad)')
    ax.set_ylabel('s (normalized flux)')
    ax.set_title('s vs phi')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # theta vs phi
    ax = axes[0, 1]
    if gc_traj is not None:
        ax.plot(gc_traj['phi'], gc_traj['theta'], 'b-', label='Guiding-center', alpha=0.7)
    if pauli_traj is not None:
        ax.plot(pauli_traj['phi'], pauli_traj['theta'], 'r-', label='Pauli particle', alpha=0.7)
    ax.set_xlabel('phi (rad)')
    ax.set_ylabel('theta (rad)')
    ax.set_title('theta vs phi')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # s vs theta (similar to R-Z projection in flux coordinates)
    ax = axes[1, 0]
    if gc_traj is not None:
        ax.plot(gc_traj['theta'], gc_traj['s'], 'b-', label='Guiding-center', alpha=0.7)
    if pauli_traj is not None:
        ax.plot(pauli_traj['theta'], pauli_traj['s'], 'r-', label='Pauli particle', alpha=0.7)
    ax.set_xlabel('theta (rad)')
    ax.set_ylabel('s (normalized flux)')
    ax.set_title('s vs theta')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Time evolution of s
    ax = axes[1, 1]
    if gc_traj is not None:
        ax.plot(gc_traj['time'], gc_traj['s'], 'b-', label='Guiding-center', alpha=0.7)
    if pauli_traj is not None:
        ax.plot(pauli_traj['time'], pauli_traj['s'], 'r-', label='Pauli particle', alpha=0.7)
    ax.set_xlabel('time (s)')
    ax.set_ylabel('s (normalized flux)')
    ax.set_title('s vs time')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Plot saved to: {output_file}")
    plt.close()


def main():
    if len(sys.argv) < 2:
        # Use default test directory
        work_dir = '/home/ert/data/SQUID/SIMPLE/losses/meiss_coils'
    else:
        work_dir = sys.argv[1]

    if not os.path.exists(work_dir):
        print(f"Error: Directory {work_dir} does not exist")
        sys.exit(1)

    print(f"Working directory: {work_dir}")

    # Check for required files
    wout_file = os.path.join(work_dir, 'wout.nc')
    coils_file = os.path.join(work_dir, 'coils.5C')

    if not os.path.exists(wout_file):
        download_test_data(work_dir)

    if not os.path.exists(coils_file):
        print(f"Error: Coils file {coils_file} not found")
        print("Please provide a coils file for your configuration")
        sys.exit(1)

    # Run with guiding-center (orbit_model=0)
    print("\n" + "="*60)
    print("Running guiding-center simulation (orbit_model=0)")
    print("="*60)
    gc_traj = run_simple(work_dir, orbit_model=0, trace_time=1e-6, ntimstep=100)

    # Run with Pauli particle (orbit_model=1)
    print("\n" + "="*60)
    print("Running Pauli particle simulation (orbit_model=1)")
    print("="*60)
    pauli_traj = run_simple(work_dir, orbit_model=1, trace_time=1e-6, ntimstep=100)

    # Generate comparison plot
    print("\n" + "="*60)
    print("Generating comparison plot")
    print("="*60)

    output_file = os.path.join(work_dir, 'gc_pauli_comparison.png')
    plot_comparison(gc_traj, pauli_traj, output_file)

    print("\n" + "="*60)
    print("Results:")
    print("="*60)
    print(f"  Guiding-center: {'OK' if gc_traj is not None else 'FAILED'}")
    print(f"  Pauli particle: {'OK' if pauli_traj is not None else 'FAILED'}")
    if HAS_MATPLOTLIB:
        print(f"  Plot: {output_file}")


if __name__ == '__main__':
    main()
