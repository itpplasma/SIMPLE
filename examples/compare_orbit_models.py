#!/usr/bin/env python3
"""Compare guiding-center, Pauli, and full orbit models.

This script runs SIMPLE with different orbit_model settings and plots
the trajectories in R-Z and R-phi projections for comparison.

Usage:
    python compare_orbit_models.py [working_directory]

The working directory must contain:
    - wout.nc: VMEC equilibrium file
    - coils.5C or similar: Coils file for Biot-Savart field
    - simple.in: Configuration file (will be modified temporarily)
"""

import os
import sys
import subprocess
import tempfile
import shutil
import numpy as np
import matplotlib.pyplot as plt


def run_simple(work_dir, orbit_model, ntimstep=100, trace_time=1e-4):
    """Run SIMPLE with specified orbit model and return trajectory."""

    # Read original config
    config_path = os.path.join(work_dir, 'simple.in')
    with open(config_path, 'r') as f:
        original_config = f.read()

    # Create modified config with orbit_model
    # Parse and modify
    lines = original_config.strip().split('\n')
    new_lines = []
    has_orbit_model = False
    has_ntimstep = False
    has_trace_time = False
    has_ntestpart = False

    for line in lines:
        # Skip existing orbit_model, ntimstep, trace_time lines
        if 'orbit_model' in line.lower() and '=' in line:
            has_orbit_model = True
            new_lines.append(f'orbit_model = {orbit_model}')
            continue
        if 'ntimstep' in line.lower() and '=' in line:
            has_ntimstep = True
            new_lines.append(f'ntimstep = {ntimstep}')
            continue
        if 'trace_time' in line.lower() and '=' in line:
            has_trace_time = True
            new_lines.append(f'trace_time = {trace_time}')
            continue
        if 'ntestpart' in line.lower() and '=' in line:
            has_ntestpart = True
            new_lines.append('ntestpart = 1')
            continue
        if line.strip() == '/':
            # Add missing parameters before closing
            if not has_orbit_model:
                new_lines.append(f'orbit_model = {orbit_model}')
            if not has_ntimstep:
                new_lines.append(f'ntimstep = {ntimstep}')
            if not has_trace_time:
                new_lines.append(f'trace_time = {trace_time}')
            if not has_ntestpart:
                new_lines.append('ntestpart = 1')
        new_lines.append(line)

    modified_config = '\n'.join(new_lines)

    # Write modified config
    with open(config_path, 'w') as f:
        f.write(modified_config)

    try:
        # Run SIMPLE
        simple_exe = os.path.join(os.path.dirname(__file__), '..', 'build', 'simple.x')
        result = subprocess.run(
            [simple_exe],
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=300
        )

        if result.returncode != 0:
            print(f"SIMPLE failed with orbit_model={orbit_model}")
            print("STDOUT:", result.stdout[-2000:] if len(result.stdout) > 2000 else result.stdout)
            print("STDERR:", result.stderr[-2000:] if len(result.stderr) > 2000 else result.stderr)
            return None

        # Read trajectory from times_lost.dat or start.dat
        # For now, just check if it ran successfully
        times_lost_path = os.path.join(work_dir, 'times_lost.dat')
        if os.path.exists(times_lost_path):
            data = np.loadtxt(times_lost_path)
            return data

        return None

    finally:
        # Restore original config
        with open(config_path, 'w') as f:
            f.write(original_config)


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
    print()

    # Check required files
    required_files = ['wout.nc', 'simple.in']
    for f in required_files:
        if not os.path.exists(os.path.join(work_dir, f)):
            print(f"Error: Required file {f} not found in {work_dir}")
            sys.exit(1)

    # Run with guiding-center (orbit_model=0)
    print("Running with guiding-center orbit model (orbit_model=0)...")
    gc_result = run_simple(work_dir, orbit_model=0, ntimstep=1000, trace_time=1e-3)

    # Run with Pauli particle (orbit_model=1)
    print("Running with Pauli particle orbit model (orbit_model=1)...")
    pauli_result = run_simple(work_dir, orbit_model=1, ntimstep=1000, trace_time=1e-3)

    # Run with full orbit (orbit_model=2)
    print("Running with full orbit model (orbit_model=2)...")
    full_result = run_simple(work_dir, orbit_model=2, ntimstep=1000, trace_time=1e-3)

    print()
    print("Results:")
    print(f"  Guiding-center: {'OK' if gc_result is not None else 'FAILED'}")
    print(f"  Pauli particle: {'OK' if pauli_result is not None else 'FAILED'}")
    print(f"  Full orbit: {'OK' if full_result is not None else 'FAILED'}")


if __name__ == '__main__':
    main()
