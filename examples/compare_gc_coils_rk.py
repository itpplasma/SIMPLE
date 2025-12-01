#!/usr/bin/env python3
"""Compare guiding-center orbits: VMEC vs coils field, RK4 vs Cash–Karp RK.

This script runs SIMPLE four times:
  1) VMEC GC, fixed-step RK4           (integmode = -2)
  2) Coils GC, fixed-step RK4          (integmode = -2, coils magfie backend)
  3) VMEC GC, adaptive Cash–Karp RK5(4) (integmode = -1)
  4) Coils GC, adaptive Cash–Karp RK5(4) (integmode = -1, coils magfie backend)
All runs use VMEC reference coordinates (s, theta, phi) via NetCDF output.
"""

import os
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

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

try:
    # Reuse NCSX data downloader from plot script
    from plot_ncsx_vmec_coils import ensure_ncsx_files  # type: ignore
except Exception:
    ensure_ncsx_files = None


def run_simple_subprocess(simple_exe: str, base_dir: str, tag: str, config_text: str) -> None:
    """
    Run SIMPLE in a dedicated subdirectory (base_dir/tag) to avoid file clashes.

    This function only runs the code; reading orbits.nc happens separately.
    """
    work_dir = os.path.join(base_dir, tag)
    os.makedirs(work_dir, exist_ok=True)

    config_path = os.path.join(work_dir, "simple.in")
    with open(config_path, "w") as f:
        f.write(config_text)

    # Symlink field files from base_dir (they must exist there already)
    for fname in ("wout.nc", "coils.simple"):
        src = os.path.join(base_dir, fname)
        dst = os.path.join(work_dir, fname)
        if os.path.exists(src) and not os.path.exists(dst):
            os.symlink(src, dst)

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
        raise RuntimeError(f"SIMPLE run '{tag}' failed with exit code {result.returncode}")


def read_orbit(work_dir: str):
    """Read s, theta, phi from orbits.nc in work_dir."""
    if nc is None:
        print("netCDF4 not available; cannot read orbits.nc")
        return None

    orbits_path = os.path.join(work_dir, "orbits.nc")
    if not os.path.exists(orbits_path):
        print(f"orbits.nc not found in {work_dir}")
        return None

    with nc.Dataset(orbits_path, "r") as ds:
        s = ds.variables["s"][:, 0]
        theta = ds.variables["theta"][:, 0]
        phi = ds.variables["phi"][:, 0]

    return {"s": s, "theta": theta, "phi": phi}


def main():
    if len(sys.argv) < 2:
        work_dir = os.path.join(os.getcwd(), "examples")
    else:
        work_dir = sys.argv[1]

    work_dir = os.path.abspath(work_dir)
    os.makedirs(work_dir, exist_ok=True)
    print(f"Working directory: {work_dir}")

    # Ensure required field files exist. Prefer local files; otherwise,
    # for NCSX runs, lazily download from STELLOPT via ensure_ncsx_files.
    wout_path = os.path.join(work_dir, "wout.nc")
    coils_path = os.path.join(work_dir, "coils.simple")

    missing_wout = not os.path.exists(wout_path)
    missing_coils = not os.path.exists(coils_path)

    if (missing_wout or missing_coils) and ensure_ncsx_files is not None:
        script_dir = os.path.abspath(os.path.dirname(__file__))
        test_data_dir = os.path.join(script_dir, "..", "golden_record", "test_data")
        test_data_dir = os.path.abspath(test_data_dir)
        os.makedirs(test_data_dir, exist_ok=True)

        ncsx_wout, coils_simple = ensure_ncsx_files(test_data_dir)

        if missing_wout:
            if os.path.islink(wout_path) or os.path.exists(wout_path):
                os.remove(wout_path)
            os.symlink(ncsx_wout, wout_path)
        if missing_coils:
            if os.path.islink(coils_path) or os.path.exists(coils_path):
                os.remove(coils_path)
            os.symlink(coils_simple, coils_path)

    for f in ("wout.nc", "coils.simple"):
        if not os.path.exists(os.path.join(work_dir, f)):
            print(f"Error: Required file {f} not found in {work_dir}")
            sys.exit(1)

    trace_time = 1e-4
    ntimstep = 1000
    facE_al = 1000.0  # Reduce orbit energy by factor 1000

    # 1) VMEC GC, fixed-step RK4 (no field_input)
    vmec_config_rk4 = f"""&config
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

    # 2) Coils GC, fixed-step RK4: coils magfie backend in VMEC reference frame
    coils_config_rk4 = f"""&config
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

    # 3) VMEC GC, adaptive Cash–Karp RK5(4)
    vmec_config_ck = f"""&config
trace_time = {trace_time}
sbeg = 0.6d0
ntestpart = 1
ntimstep = {ntimstep}
netcdffile = 'wout.nc'
isw_field_type = 1        ! VMEC
integmode = -1            ! Adaptive Cash–Karp RK5(4) in VMEC field
npoiper2 = 128
deterministic = .True.
facE_al = {facE_al}
contr_pp = -1000d0
output_orbits_macrostep = .True.
/
"""

    # 4) Coils GC, adaptive Cash–Karp RK5(4)
    coils_config_ck = f"""&config
trace_time = {trace_time}
sbeg = 0.6d0
ntestpart = 1
ntimstep = {ntimstep}
netcdffile = 'wout.nc'
isw_field_type = 5        ! Coils magfie backend (VMEC reference coordinates)
integmode = -1
npoiper2 = 128
deterministic = .True.
facE_al = {facE_al}
field_input = 'coils.simple'
contr_pp = -1000d0
output_orbits_macrostep = .True.
/
"""
    simple_exe = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "build", "simple.x"))
    if not os.path.exists(simple_exe):
        print(f"Error: SIMPLE executable not found at {simple_exe}")
        sys.exit(1)

    # Run all four cases in parallel, each in its own subdirectory
    configs = {
        "vmec_rk4": vmec_config_rk4,
        "coils_rk4": coils_config_rk4,
        "vmec_ck": vmec_config_ck,
        "coils_ck": coils_config_ck,
    }

    print("Running GC orbits (RK4 and Cash–Karp) in parallel...")
    with ThreadPoolExecutor(max_workers=len(configs)) as executor:
        futures = {
            executor.submit(run_simple_subprocess, simple_exe, work_dir, tag, cfg): tag
            for tag, cfg in configs.items()
        }
        for fut in as_completed(futures):
            tag = futures[fut]
            try:
                fut.result()
            except Exception as exc:
                print(f"Run '{tag}' failed: {exc}")
                sys.exit(1)

    # Read back trajectories
    vmec_rk4_traj = read_orbit(os.path.join(work_dir, "vmec_rk4"))
    coils_rk4_traj = read_orbit(os.path.join(work_dir, "coils_rk4"))
    vmec_ck_traj = read_orbit(os.path.join(work_dir, "vmec_ck"))
    coils_ck_traj = read_orbit(os.path.join(work_dir, "coils_ck"))

    if any(traj is None for traj in (vmec_rk4_traj, coils_rk4_traj, vmec_ck_traj, coils_ck_traj)):
        print("One of the runs failed or produced no orbits; no plot generated")
        sys.exit(1)

    if not HAS_MATPLOTLIB:
        print("matplotlib not available; not plotting")
        sys.exit(0)

    # Mask out NaNs and wrap angles
    s_vmec_rk4 = vmec_rk4_traj["s"]
    theta_vmec_rk4 = np.mod(vmec_rk4_traj["theta"], 2.0 * np.pi)
    phi_vmec_rk4 = np.mod(vmec_rk4_traj["phi"], 2.0 * np.pi)

    s_coils_rk4 = coils_rk4_traj["s"]
    theta_coils_rk4 = np.mod(coils_rk4_traj["theta"], 2.0 * np.pi)
    phi_coils_rk4 = np.mod(coils_rk4_traj["phi"], 2.0 * np.pi)

    s_vmec_ck = vmec_ck_traj["s"]
    theta_vmec_ck = np.mod(vmec_ck_traj["theta"], 2.0 * np.pi)
    phi_vmec_ck = np.mod(vmec_ck_traj["phi"], 2.0 * np.pi)

    s_coils_ck = coils_ck_traj["s"]
    theta_coils_ck = np.mod(coils_ck_traj["theta"], 2.0 * np.pi)
    phi_coils_ck = np.mod(coils_ck_traj["phi"], 2.0 * np.pi)

    # Reconstruct physical time from trace_time and number of steps
    n_steps = s_vmec_rk4.shape[0]
    t_axis = np.linspace(0.0, trace_time, n_steps)

    mask_vmec_rk4 = ~np.isnan(s_vmec_rk4)
    mask_coils_rk4 = ~np.isnan(s_coils_rk4)
    mask_vmec_ck = ~np.isnan(s_vmec_ck)
    mask_coils_ck = ~np.isnan(s_coils_ck)

    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    axes[0].plot(t_axis[mask_vmec_rk4], s_vmec_rk4[mask_vmec_rk4], "b-", label="VMEC RK4")
    axes[0].plot(t_axis[mask_coils_rk4], s_coils_rk4[mask_coils_rk4], "r--", label="Coils RK4")
    axes[0].plot(t_axis[mask_vmec_ck], s_vmec_ck[mask_vmec_ck], "b:", label="VMEC Cash–Karp")
    axes[0].plot(t_axis[mask_coils_ck], s_coils_ck[mask_coils_ck], "r-.", label="Coils Cash–Karp")
    axes[0].set_ylabel("s")
    axes[0].set_title("Guiding-center RK: VMEC vs coils (VMEC reference coords)")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t_axis[mask_vmec_rk4], theta_vmec_rk4[mask_vmec_rk4], "b-")
    axes[1].plot(t_axis[mask_coils_rk4], theta_coils_rk4[mask_coils_rk4], "r--")
    axes[1].plot(t_axis[mask_vmec_ck], theta_vmec_ck[mask_vmec_ck], "b:")
    axes[1].plot(t_axis[mask_coils_ck], theta_coils_ck[mask_coils_ck], "r-.")
    axes[1].set_ylabel("theta mod 2π")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t_axis[mask_vmec_rk4], phi_vmec_rk4[mask_vmec_rk4], "b-")
    axes[2].plot(t_axis[mask_coils_rk4], phi_coils_rk4[mask_coils_rk4], "r--")
    axes[2].plot(t_axis[mask_vmec_ck], phi_vmec_ck[mask_vmec_ck], "b:")
    axes[2].plot(t_axis[mask_coils_ck], phi_coils_ck[mask_coils_ck], "r-.")
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
