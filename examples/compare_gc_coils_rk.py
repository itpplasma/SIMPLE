#!/usr/bin/env python3
"""Compare RK guiding-center orbits: VMEC vs coils field.

This script runs SIMPLE twice:
  1) VMEC RK guiding-center (integmode = -1)
  2) Coil-based RK guiding-center using coils magfie backend
Both runs use VMEC reference coordinates (s, theta, phi) via NetCDF output.
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


def run_simple(config_text: str, work_dir: str, tag: str):
    """Write a config file, run SIMPLE, and return trajectory from orbits.nc."""
    if nc is None:
        print("netCDF4 not available; cannot read orbits.nc")
        return None

    config_path = os.path.join(work_dir, f'simple_{tag}.in')
    with open(config_path, 'w') as f:
        f.write(config_text)

    simple_exe = os.path.join(os.path.dirname(__file__), '..', 'build', 'simple.x')
    result = subprocess.run(
        [simple_exe, config_path],
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=300,
    )

    if result.returncode != 0:
        print(f"SIMPLE run '{tag}' failed")
        print("STDOUT:", result.stdout[-2000:] if len(result.stdout) > 2000 else result.stdout)
        print("STDERR:", result.stderr[-2000:] if len(result.stderr) > 2000 else result.stderr)
        return None

    orbits_path = os.path.join(work_dir, 'orbits.nc')
    if not os.path.exists(orbits_path):
        print(f"orbits.nc not found for run '{tag}'")
        return None
    with nc.Dataset(orbits_path, 'r') as ds:
        # Variables are stored as (timestep, particle)
        s = ds.variables['s'][:, 0]
        theta = ds.variables['theta'][:, 0]
        phi = ds.variables['phi'][:, 0]

    return {"s": s, "theta": theta, "phi": phi}


def main():
    if len(sys.argv) < 2:
        work_dir = os.path.join(os.getcwd(), "examples")
    else:
        work_dir = sys.argv[1]

    work_dir = os.path.abspath(work_dir)
    print(f"Working directory: {work_dir}")

    for f in ("wout.nc", "coils.simple"):
        if not os.path.exists(os.path.join(work_dir, f)):
            print(f"Error: Required file {f} not found in {work_dir}")
            sys.exit(1)

    trace_time = 1e-4
    ntimstep = 1000
    facE_al = 1000.0  # Reduce orbit energy by factor 1000

    # 1) VMEC RK guiding-center: VMEC field, no field_input
    vmec_config = f"""&config
trace_time = {trace_time}
sbeg = 0.6d0
ntestpart = 1
ntimstep = {ntimstep}
netcdffile = 'wout.nc'
isw_field_type = 1        ! VMEC
integmode = -2            ! Fixed-step RK4 in VMEC field
npoiper2 = 128
deterministic = .True.
facE_al = {facE_al}
contr_pp = -1000d0
output_orbits_macrostep = .True.
/
"""

    # 2) Coils RK guiding-center: coils magfie backend in VMEC reference frame
    coils_config = f"""&config
trace_time = {trace_time}
sbeg = 0.6d0
ntestpart = 1
ntimstep = {ntimstep}
netcdffile = 'wout.nc'
isw_field_type = 5        ! Coils magfie backend (VMEC reference coordinates)
integmode = -2
npoiper2 = 128
deterministic = .True.
facE_al = {facE_al}
field_input = 'coils.simple'
contr_pp = -1000d0
output_orbits_macrostep = .True.
/
"""

    print("Running VMEC RK guiding-center orbit...")
    vmec_traj = run_simple(vmec_config, work_dir, "vmec_rk")

    print("Running coils RK guiding-center orbit (coils backend)...")
    coils_traj = run_simple(coils_config, work_dir, "coils_rk")

    if vmec_traj is None or coils_traj is None:
        print("One of the runs failed; no plot generated")
        sys.exit(1)

    if not HAS_MATPLOTLIB:
        print("matplotlib not available; not plotting")
        sys.exit(0)

    # Mask out NaNs and wrap angles
    s_vmec = vmec_traj["s"]
    theta_vmec = np.mod(vmec_traj["theta"], 2.0 * np.pi)
    phi_vmec = np.mod(vmec_traj["phi"], 2.0 * np.pi)
    s_coils = coils_traj["s"]
    theta_coils = np.mod(coils_traj["theta"], 2.0 * np.pi)
    phi_coils = np.mod(coils_traj["phi"], 2.0 * np.pi)

    # Reconstruct physical time from trace_time and number of steps
    n_steps = s_vmec.shape[0]
    t_axis = np.linspace(0.0, trace_time, n_steps)

    mask_vmec = ~np.isnan(s_vmec)
    mask_coils = ~np.isnan(s_coils)

    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    axes[0].plot(t_axis[mask_vmec], s_vmec[mask_vmec], "b-", label="VMEC RK")
    axes[0].plot(t_axis[mask_coils], s_coils[mask_coils], "r--", label="Coils RK")
    axes[0].set_ylabel("s")
    axes[0].set_title("Guiding-center RK: VMEC vs coils (VMEC reference coords)")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t_axis[mask_vmec], theta_vmec[mask_vmec], "b-")
    axes[1].plot(t_axis[mask_coils], theta_coils[mask_coils], "r--")
    axes[1].set_ylabel("theta mod 2π")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t_axis[mask_vmec], phi_vmec[mask_vmec], "b-")
    axes[2].plot(t_axis[mask_coils], phi_coils[mask_coils], "r--")
    axes[2].set_ylabel("phi mod 2π")
    axes[2].set_xlabel("time")
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    output_png = os.path.join(work_dir, "gc_vmec_vs_coils_rk.png")
    plt.savefig(output_png, dpi=150)
    plt.close()

    print(f"Comparison plot written to: {output_png}")


if __name__ == "__main__":
    main()
