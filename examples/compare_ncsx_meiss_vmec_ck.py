#!/usr/bin/env python3
"""
Compare NCSX guiding-center orbits: Cash–Karp VMEC vs Meiss (coils).

This script:
  - Ensures the NCSX free-boundary VMEC equilibrium and coils file from
    the STELLOPT tutorial are available locally (same URLs as in
    test/golden_record/run_golden_tests.sh).
  - Runs SIMPLE twice in a dedicated working directory:
      * VMEC GC, adaptive Cash–Karp RK5(4) (integmode = -1, isw_field_type = 1)
      * Meiss GC from coils, adaptive Cash–Karp RK5(4)
        (integmode = -1, isw_field_type = 3, field_input = 'coils.simple')
  - Both runs use the same NCSX wout and coils files in VMEC reference
    coordinates.
  - Reads orbits.nc from both runs and plots s, theta mod 2π, phi mod 2π
    vs time for visual comparison.
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


NCSX_WOUT_URL = (
    "https://github.com/hiddenSymmetries/simsopt/raw/master/"
    "tests/test_files/wout_c09r00_fixedBoundary_0.5T_vacuum_ns201.nc"
)

NCSX_COILS_URL = (
    "https://princetonuniversity.github.io/STELLOPT/examples/coils.c09r00"
)


def ensure_ncsx_files(test_data_dir: str):
    """
    Ensure NCSX wout_ncsx.nc and coils.c09r00.simple exist in test_data_dir.

    Returns (wout_path, coils_simple_path).
    """
    from urllib.request import urlretrieve

    os.makedirs(test_data_dir, exist_ok=True)

    wout_ncsx = os.path.join(test_data_dir, "wout_ncsx.nc")
    coils_stellopt = os.path.join(test_data_dir, "coils.c09r00")
    coils_simple = os.path.join(test_data_dir, "coils.c09r00.simple")

    if not os.path.exists(wout_ncsx):
        print(f"Downloading NCSX wout: {NCSX_WOUT_URL}")
        urlretrieve(NCSX_WOUT_URL, wout_ncsx)

    if not os.path.exists(coils_simple):
        if not os.path.exists(coils_stellopt):
            print(f"Downloading NCSX coils: {NCSX_COILS_URL}")
            urlretrieve(NCSX_COILS_URL, coils_stellopt)

        # Convert STELLOPT coils file to SIMPLE's coils.simple format
        print("Converting NCSX coils to SIMPLE format (coils.c09r00.simple)")
        coords = []
        currents = []
        in_filament = False
        with open(coils_stellopt, "r") as f:
            for line in f:
                line = line.strip()
                lower = line.lower()
                if "begin filament" in lower:
                    in_filament = True
                    continue
                if not in_filament or not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                try:
                    x, y, z, I = map(float, parts[:4])
                except ValueError:
                    continue
                coords.append((x, y, z))
                currents.append(I)
        if not coords:
            raise RuntimeError("No filament coordinates found in coils.c09r00")
        with open(coils_simple, "w") as f:
            f.write(f"{len(coords)}\n")
            for (x, y, z), I in zip(coords, currents):
                f.write(f"{x:.14E}   {y:.14E}   {z:.14E}   {I:.14E}\n")

    return wout_ncsx, coils_simple


def run_simple(simple_exe: str, work_dir: str, tag: str, config_text: str) -> None:
    """
    Run SIMPLE in work_dir/tag with the given namelist config.
    """
    run_dir = os.path.join(work_dir, tag)
    os.makedirs(run_dir, exist_ok=True)

    cfg_path = os.path.join(run_dir, "simple.in")
    with open(cfg_path, "w") as f:
        f.write(config_text)

    # Link field files from work_dir
    for fname in ("wout.nc", "coils.simple"):
        src = os.path.join(work_dir, fname)
        dst = os.path.join(run_dir, fname)
        if os.path.exists(src) and not os.path.exists(dst):
            os.symlink(src, dst)

    print(f"Running SIMPLE [{tag}] in {run_dir}")
    res = subprocess.run(
        [simple_exe, cfg_path],
        cwd=run_dir,
        capture_output=True,
        text=True,
        timeout=1800,
    )
    if res.returncode != 0:
        print(f"SIMPLE run '{tag}' failed")
        print("STDOUT:", res.stdout[-2000:] if len(res.stdout) > 2000 else res.stdout)
        print("STDERR:", res.stderr[-2000:] if len(res.stderr) > 2000 else res.stderr)
        raise RuntimeError(f"SIMPLE run '{tag}' failed with exit code {res.returncode}")


def read_orbit(orbits_path: str):
    """
    Read s, theta, phi, time from orbits.nc (first particle).
    """
    if nc is None:
        raise RuntimeError("netCDF4 is not available")
    if not os.path.exists(orbits_path):
        raise FileNotFoundError(orbits_path)

    with nc.Dataset(orbits_path, "r") as ds:
        s = ds.variables["s"][:, 0]
        theta = ds.variables["theta"][:, 0]
        phi = ds.variables["phi"][:, 0]
        t = ds.variables["time"][:, 0]
    return t, s, theta, phi


def main():
    script_dir = os.path.abspath(os.path.dirname(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, ".."))

    if len(sys.argv) > 1:
        work_dir = os.path.abspath(sys.argv[1])
    else:
        work_dir = os.path.join(repo_root, "build", "test", "ncsx_meiss_vmec_ck")
    os.makedirs(work_dir, exist_ok=True)
    print(f"Working directory: {work_dir}")

    # Ensure NCSX wout and coils from STELLOPT exist under golden_record test_data
    test_data_dir = os.path.join(repo_root, "golden_record", "test_data")
    wout_ncsx, coils_simple = ensure_ncsx_files(test_data_dir)

    # Link them into work_dir as wout.nc and coils.simple
    wout_link = os.path.join(work_dir, "wout.nc")
    coils_link = os.path.join(work_dir, "coils.simple")
    for src, dst in ((wout_ncsx, wout_link), (coils_simple, coils_link)):
        if os.path.islink(dst) or os.path.exists(dst):
            os.remove(dst)
        os.symlink(src, dst)

    simple_exe = os.path.join(repo_root, "build", "simple.x")
    if not os.path.exists(simple_exe):
        print(f"Error: SIMPLE executable not found at {simple_exe}")
        sys.exit(1)

    trace_time = 1e-4
    ntimstep = 1000
    facE_al = 1000.0

    vmec_ck_cfg = f"""&config
trace_time = {trace_time}
sbeg = 0.6d0
ntestpart = 1
ntimstep = {ntimstep}
netcdffile = 'wout.nc'
isw_field_type = 1        ! VMEC
integmode = -1            ! Cash–Karp RK5(4)
npoiper2 = 128
deterministic = .True.
facE_al = {facE_al}
contr_pp = -1000d0
output_orbits_macrostep = .True.
/
"""

    meiss_ck_cfg = f"""&config
trace_time = {trace_time}
sbeg = 0.6d0
ntestpart = 1
ntimstep = {ntimstep}
netcdffile = 'wout.nc'
isw_field_type = 3        ! Meiss canonical coordinates
integmode = -1            ! Cash–Karp RK5(4)
npoiper2 = 128
deterministic = .True.
facE_al = {facE_al}
field_input = 'coils.simple'
contr_pp = -1000d0
output_orbits_macrostep = .True.
/
"""

    # Run both cases
    run_simple(simple_exe, work_dir, "vmec_ck", vmec_ck_cfg)
    run_simple(simple_exe, work_dir, "meiss_ck", meiss_ck_cfg)

    # Read trajectories
    t_vmec, s_vmec, th_vmec, ph_vmec = read_orbit(os.path.join(work_dir, "vmec_ck", "orbits.nc"))
    t_meiss, s_meiss, th_meiss, ph_meiss = read_orbit(os.path.join(work_dir, "meiss_ck", "orbits.nc"))

    if not np.allclose(t_vmec, t_meiss):
        print("Warning: time axes differ between VMEC and Meiss; plotting VMEC times only")
    t_axis = t_vmec

    if not HAS_MATPLOTLIB:
        print("matplotlib not available; cannot plot")
        sys.exit(0)

    # Wrap angles
    th_vmec_mod = np.mod(th_vmec, 2.0 * np.pi)
    ph_vmec_mod = np.mod(ph_vmec, 2.0 * np.pi)
    th_meiss_mod = np.mod(th_meiss, 2.0 * np.pi)
    ph_meiss_mod = np.mod(ph_meiss, 2.0 * np.pi)

    mask_vmec = ~np.isnan(s_vmec)
    mask_meiss = ~np.isnan(s_meiss)

    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    axes[0].plot(t_axis[mask_vmec], s_vmec[mask_vmec],
                 color="C0", linestyle="-", linewidth=1.8, label="VMEC Cash–Karp")
    axes[0].plot(t_axis[mask_meiss], s_meiss[mask_meiss],
                 color="C1", linestyle="--", linewidth=1.8, label="Meiss(coils) Cash–Karp")
    axes[0].set_ylabel("s")
    axes[0].set_title("NCSX: GC Cash–Karp, VMEC vs Meiss (coils)")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t_axis[mask_vmec], th_vmec_mod[mask_vmec],
                 color="C0", linestyle="-", linewidth=1.8)
    axes[1].plot(t_axis[mask_meiss], th_meiss_mod[mask_meiss],
                 color="C1", linestyle="--", linewidth=1.8)
    axes[1].set_ylabel("theta mod 2π")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t_axis[mask_vmec], ph_vmec_mod[mask_vmec],
                 color="C0", linestyle="-", linewidth=1.8)
    axes[2].plot(t_axis[mask_meiss], ph_meiss_mod[mask_meiss],
                 color="C1", linestyle="--", linewidth=1.8)
    axes[2].set_ylabel("phi mod 2π")
    axes[2].set_xlabel("time")
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = os.path.join(work_dir, "gc_meiss_vs_vmec_ncsx_ck.png")
    plt.savefig(out_png, dpi=150)
    plt.close()

    print(f"Comparison plot written to: {out_png}")


if __name__ == "__main__":
    main()

