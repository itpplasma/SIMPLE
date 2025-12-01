#!/usr/bin/env python3
"""Plot a single VMEC RK guiding-center orbit.

This script runs SIMPLE in VMEC Runge–Kutta (integmode = -1) mode
for a single particle and produces a simple diagnostic plot.

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
import numpy as np

try:
    import netCDF4 as nc
except ImportError:
    nc = None

try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def run_simple_vmec_rk(work_dir, ntimstep=200, trace_time=1e-4):
    """Run SIMPLE in VMEC RK mode (integmode = -1, orbit_model = 0)."""

    config = f"""&config
trace_time = {trace_time}
sbeg = 0.6d0
ntestpart = 1
ntimstep = {ntimstep}
netcdffile = 'wout.nc'
isw_field_type = 1
integmode = -1
npoiper2 = 128
deterministic = .True.
orbit_model = 0
n_e = 2
n_d = 4
contr_pp = -1000d0
output_orbits_macrostep = .True.
/
"""

    config_path = os.path.join(work_dir, 'simple_vmec_rk.in')
    with open(config_path, 'w') as f:
        f.write(config)

    simple_exe = os.path.join(os.path.dirname(__file__), '..', 'build', 'simple.x')
    result = subprocess.run(
        [simple_exe, config_path],
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=300
    )

    if result.returncode != 0:
        print("SIMPLE VMEC RK run failed")
        print("STDOUT:", result.stdout[-2000:] if len(result.stdout) > 2000 else result.stdout)
        print("STDERR:", result.stderr[-2000:] if len(result.stderr) > 2000 else result.stderr)
        return None

    orbits_path = os.path.join(work_dir, 'orbits.nc')
    if nc is None or not os.path.exists(orbits_path):
        print("orbits.nc not found or netCDF4 not available")
        return None

    with nc.Dataset(orbits_path, 'r') as ds:
        # Variables are stored as (timestep, particle)
        s = ds.variables['s'][:, 0]
        theta = ds.variables['theta'][:, 0]
        phi = ds.variables['phi'][:, 0]
        time = ds.variables['time'][:, 0]

    return {
        's': s,
        'theta': theta,
        'phi': phi,
        'time': time,
    }


def plot_vmec_rk(traj, output_file):
    """Plot VMEC RK orbit (s, theta, phi vs time)."""
    if not HAS_MATPLOTLIB:
        print("matplotlib not available; skipping plot generation")
        return

    s = traj['s']
    theta = traj['theta']
    phi = traj['phi']
    time = traj['time']

    mask = ~np.isnan(s)

    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    axes[0].plot(time[mask], s[mask], 'b-')
    axes[0].set_ylabel('s')
    axes[0].set_title('VMEC RK orbit: s, theta, phi vs time')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(time[mask], theta[mask], 'r-')
    axes[1].set_ylabel('theta')
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(time[mask], phi[mask], 'g-')
    axes[2].set_ylabel('phi')
    axes[2].set_xlabel('time')
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
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
    print()

    # Check required files
    required_files = ['wout.nc', 'simple.in']
    for f in required_files:
        if not os.path.exists(os.path.join(work_dir, f)):
            print(f"Error: Required file {f} not found in {work_dir}")
            sys.exit(1)

    trace_time = 1e-4
    ntimstep = 1000

    print("Running VMEC RK guiding-center orbit (integmode = -1, orbit_model = 0)...")
    traj = run_simple_vmec_rk(work_dir, ntimstep=ntimstep, trace_time=trace_time)

    if traj is None:
        print("VMEC RK run failed; no plot generated")
        sys.exit(1)

    output_file = os.path.join(work_dir, 'vmec_rk_orbit.png')
    plot_vmec_rk(traj, output_file)
    if HAS_MATPLOTLIB:
        print(f"  Plot written to: {output_file}")


if __name__ == '__main__':
    main()
