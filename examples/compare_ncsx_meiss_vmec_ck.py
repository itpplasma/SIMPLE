#!/usr/bin/env python3
"""Compare NCSX GC orbits integrated in VMEC vs Meiss canonical coordinates.

Both runs use the same VMEC Fourier B-field with Cash-Karp RK5(4). The VMEC
leg integrates in (s, theta, phi); the Meiss leg integrates in the canonical
(r, theta_c, phi_c) frame whose construction (Meiss ODE + spline of the
canonical field components and lambda/chi gauge) is built from the same
VMEC field. Agreement at the macrostep grid guards against regressions in
the Meiss coordinate construction, batch-spline accuracy, and the
integ_to_ref / ref_to_integ reference-coordinate mapping.

This is NOT a coils-vs-VMEC field cross-check: with integmode=-1 SIMPLE
does not load field_input, so both legs share the VMEC Fourier B-field.

Fails (non-zero exit) if the final-position deviation exceeds tolerance.
"""

import os
import sys
import subprocess
from urllib.request import urlretrieve

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

WOUT_URL = ("https://github.com/hiddenSymmetries/simsopt/raw/master/"
            "tests/test_files/wout_c09r00_fixedBoundary_0.5T_vacuum_ns201.nc")

CFG = """&config
trace_time = 1.0d-4
sbeg = 0.6d0
ntestpart = 1
ntimstep = 1000
netcdffile = 'wout.nc'
isw_field_type = {field_type}
integmode = -1
npoiper2 = 128
deterministic = .True.
facE_al = 1000.0
contr_pp = -1000d0
output_orbits_macrostep = .True.
/
"""

TOL_S_ABS = 5.0e-4
TOL_ANG_ABS = 5.0e-3


def ensure_wout(data_dir):
    os.makedirs(data_dir, exist_ok=True)
    wout = os.path.join(data_dir, "wout_ncsx.nc")
    if not os.path.exists(wout):
        urlretrieve(WOUT_URL, wout)
    return wout


def run_simple(exe, work_dir, tag, cfg):
    run_dir = os.path.join(work_dir, tag)
    os.makedirs(run_dir, exist_ok=True)
    cfg_path = os.path.join(run_dir, "simple.in")
    with open(cfg_path, "w") as f:
        f.write(cfg)
    wout_link = os.path.join(run_dir, "wout.nc")
    src = os.path.join(work_dir, "wout.nc")
    if os.path.exists(src) and not os.path.exists(wout_link):
        os.symlink(src, wout_link)
    res = subprocess.run([exe, cfg_path], cwd=run_dir, capture_output=True,
                         text=True, timeout=1800)
    if res.returncode != 0:
        sys.stderr.write(res.stdout[-2000:] + res.stderr[-2000:])
        raise RuntimeError(f"SIMPLE run '{tag}' failed ({res.returncode})")
    with nc.Dataset(os.path.join(run_dir, "orbits.nc")) as ds:
        return {k: ds.variables[k][:, 0] for k in ("time", "s", "theta", "phi")}


def angular_diff(a, b):
    """Minimal signed angle a - b wrapped to [-pi, pi]."""
    d = (a - b + np.pi) % (2 * np.pi) - np.pi
    return d


def compare(vmec, meiss):
    mask = ~(np.isnan(vmec["s"]) | np.isnan(meiss["s"]))
    if mask.sum() < 2:
        raise RuntimeError("Both orbits failed almost immediately")

    last = np.where(mask)[0][-1]
    ds = abs(vmec["s"][last] - meiss["s"][last])
    dth = abs(angular_diff(vmec["theta"][last], meiss["theta"][last]))
    dph = abs(angular_diff(vmec["phi"][last], meiss["phi"][last]))

    print(f"Final-position deviation at t={vmec['time'][last]:.4e}:")
    print(f"  |Δs|     = {ds:.3e}  (tol {TOL_S_ABS:.1e})")
    print(f"  |Δθ|    = {dth:.3e}  (tol {TOL_ANG_ABS:.1e})")
    print(f"  |Δφ|    = {dph:.3e}  (tol {TOL_ANG_ABS:.1e})")

    failures = []
    if ds > TOL_S_ABS:
        failures.append(f"s deviation {ds:.3e} > {TOL_S_ABS:.1e}")
    if dth > TOL_ANG_ABS:
        failures.append(f"theta deviation {dth:.3e} > {TOL_ANG_ABS:.1e}")
    if dph > TOL_ANG_ABS:
        failures.append(f"phi deviation {dph:.3e} > {TOL_ANG_ABS:.1e}")
    return failures


def plot(vmec, meiss, out_png):
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    axes[0].set_title("NCSX RK45: VMEC coords vs Meiss canonical coords")
    for ax, key, ylabel, wrap in (
        (axes[0], "s", "s", False),
        (axes[1], "theta", "theta mod 2π", True),
        (axes[2], "phi", "phi mod 2π", True),
    ):
        for run, color, ls, label in ((vmec, "C0", "-", "VMEC coords"),
                                      (meiss, "C1", "--", "Meiss coords")):
            y = np.mod(run[key], 2 * np.pi) if wrap else run[key]
            mask = ~np.isnan(run["s"])
            ax.plot(run["time"][mask], y[mask], color=color, linestyle=ls,
                    linewidth=1.8, label=label if key == "s" else None)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
    axes[0].legend()
    axes[2].set_xlabel("time")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)


def main():
    repo = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    work_dir = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 \
        else os.path.join(repo, "build", "test", "ncsx_meiss_vmec_ck")
    os.makedirs(work_dir, exist_ok=True)

    wout = ensure_wout(os.path.join(repo, "golden_record", "test_data"))
    link = os.path.join(work_dir, "wout.nc")
    if os.path.lexists(link):
        os.remove(link)
    os.symlink(wout, link)

    exe = os.path.join(repo, "build", "simple.x")
    if not os.path.exists(exe):
        sys.exit(f"SIMPLE not built at {exe}")

    vmec = run_simple(exe, work_dir, "vmec_ck", CFG.format(field_type=1))
    meiss = run_simple(exe, work_dir, "meiss_ck", CFG.format(field_type=3))

    out_png = os.path.join(work_dir, "gc_meiss_vs_vmec_ncsx_ck.png")
    plot(vmec, meiss, out_png)
    print(f"Comparison plot: {out_png}")

    failures = compare(vmec, meiss)
    if failures:
        for f in failures:
            print(f"FAIL: {f}")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
