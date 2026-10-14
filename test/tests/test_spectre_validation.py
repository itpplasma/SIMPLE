#!/usr/bin/env python3
"""SPECTRE configuration gate, single-volume classification and determinism."""
import os
import re
import subprocess
import sys
import tempfile
import numpy as np
from spectre_test_utils import axis_volume

TRACE_TIME = 1.0e-4
NPART = 64
FACE_AL = 8.0
NPOIPER2 = 128
CLASS_MATCH_MIN = 1.0

def write_input(path, h5, ntestpart=NPART,
                trace_time=TRACE_TIME, sbeg=0.5, ntimstep=100,
                npoiper2=NPOIPER2, relerr="1d-11", face_al=FACE_AL,
                ran_seed=12345):
    lines = [
        "&config",
        f"  trace_time = {trace_time}",
        f"  sbeg = {sbeg}",
        f"  ntestpart = {ntestpart}",
        f"  ntimstep = {ntimstep}",
        f"  npoiper2 = {npoiper2}",
        f"  relerr = {relerr}",
        f"  facE_al = {face_al}",
        f"  field_input = '{h5}'",
        "  integ_coords = 6",
        "  integmode = 3",
        "  spectre_ncon_phi = 32",
        "  deterministic = .True.",
        f"  ran_seed = {ran_seed}",
        "/",
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def run(binary, workdir, h5, **kwargs):
    write_input(os.path.join(workdir, "simple.in"), h5, **kwargs)
    proc = subprocess.run([binary, "simple.in"], cwd=workdir,
                          capture_output=True, text=True, timeout=280)
    if proc.returncode != 0:
        raise RuntimeError(f"simple.x failed:\n{proc.stdout[-3000:]}\n"
                           f"{proc.stderr[-1500:]}")
    return proc.stdout


def classification(binary, h5, failures):
    with tempfile.TemporaryDirectory() as work:
        out = run(binary, work, h5, ntestpart=128, trace_time=1.0e-6, ntimstep=5,
                  relerr="1d-8")
        tl = np.loadtxt(os.path.join(work, "times_lost.dat"))
        start = np.loadtxt(os.path.join(work, "start.dat"))

    m = re.search(r"bmin = *([0-9.E+-]+) *bmax = *([0-9.E+-]+)", out)
    if not m:
        failures.append("classification: surface bmin/bmax line missing")
        return
    bmin, bmax = float(m.group(1)), float(m.group(2))

    trap_par, perp_inv, lam = tl[:, 2], tl[:, 4], start[:, 4]
    code_trapped = trap_par > 0.0
    mirror_trapped = perp_inv > 1.0 / bmax   # v_par^2 < 2 mu (Bmax - B), v = v0
    match = float(np.mean(code_trapped == mirror_trapped))
    if match < CLASS_MATCH_MIN:
        failures.append(f"classification: mirror criterion matches trap_par on "
                        f"{match:.1%} of markers (< {CLASS_MATCH_MIN:.0%})")
    if not (code_trapped.any() and (~code_trapped).any()):
        failures.append("classification: split degenerate (one class empty)")

    # perp_inv encodes a real on-surface |B| = (1 - lambda^2)/perp_inv in
    # [Bmin, Bmax] (independent physical consistency of the sampled markers).
    b_marker = (1.0 - lam ** 2) / perp_inv
    if not np.all(b_marker <= bmax * (1.0 + 1e-6)):
        failures.append("classification: reconstructed |B| exceeds surface Bmax")
    if not np.all(b_marker >= bmin * (1.0 - 1e-6)):
        failures.append("classification: reconstructed |B| below surface Bmin")

    print(f"classification: trapped={int(code_trapped.sum())}/{len(tl)} "
          f"(frac={code_trapped.mean():.3f}) mirror_match={match:.1%} "
          f"Bmax={bmax:.3e} Bmin={bmin:.3e}")


def main():
    binary, source, field_driver = sys.argv[1:]
    for mode, level in ((0, 3), (1, 0), (3, 1), (5, 1)):
        with tempfile.TemporaryDirectory() as work:
            write_input(os.path.join(work, "simple.in"), source)
            path = os.path.join(work, "simple.in")
            with open(path) as f:
                config = f.read().replace("integmode = 3", f"integmode = {mode}")
                config = config.replace("&config", f"&config\n crossing_level = {level}")
            with open(path, "w") as f:
                f.write(config)
            proc = subprocess.run([binary, "simple.in"], cwd=work,
                                  capture_output=True, text=True, timeout=30)
            output = proc.stdout + proc.stderr
            assert proc.returncode != 0, f"mode {mode} traced a multi-volume field"
            assert ("multi-volume SPECTRE tracing requires" if level == 3 else
                    "crossing_level must be 3") in output
            if level == 3:
                assert "DOC/spectre-rejected-interface-models.md" in output
            assert "spectre_construction_grid:" not in output
            assert not os.path.exists(os.path.join(work, "times_lost.dat"))
    # Removed controls must fail during namelist reading, even at old defaults.
    for key, value in (("spectre_fo_integrator", "'boris'"),
                       ("spectre_fo_integrator", "'rk45'"),
                       ("spectre_fo_layer_larmor", "0d0"),
                       ("spectre_fo_rotation", "0.05d0")):
        with tempfile.TemporaryDirectory() as work:
            path = os.path.join(work, "simple.in")
            write_input(path, source)
            with open(path) as f:
                config = f.read().replace("&config", f"&config\n {key} = {value}")
            with open(path, "w") as f:
                f.write(config)
            proc = subprocess.run([binary, "simple.in"], cwd=work,
                                  capture_output=True, text=True, timeout=30)
            assert proc.returncode != 0, f"removed control {key} was accepted"
            assert key in (proc.stdout + proc.stderr).lower()
            assert not os.path.exists(os.path.join(work, "times_lost.dat"))
    # Exercise the executable's restored multi-volume level-3 tracing path.
    with tempfile.TemporaryDirectory() as work:
        run(binary, work, source, ntestpart=8, trace_time=1.0e-6,
            sbeg=0.97, ntimstep=5, face_al=1000.0, npoiper2=1024,
            relerr="1d-13")
        records = np.loadtxt(os.path.join(work, "times_lost.dat"))
        assert records.shape[0] == 8
        assert np.all(np.isfinite(records[:, 1]))
        assert np.all(records[:, 1] > 0.0), "level-3 trace did not resolve markers"
        assert not os.path.exists(os.path.join(work, "spectre_crossing_events.dat"))
    # A library caller must not turn a sampling-only field into an orbit run.
    with tempfile.TemporaryDirectory() as work:
        proc = subprocess.run([field_driver, "reject", source], cwd=work,
                              capture_output=True, text=True, timeout=30)
        output = proc.stdout + proc.stderr
        assert "SPECTRE sampling-only configuration accepted" in output
        assert proc.returncode != 0
        assert "multi-volume SPECTRE tracing requires" in output
        assert "DOC/spectre-rejected-interface-models.md" in output
    failures = []
    with axis_volume(source) as field:
        classification(binary, field, failures)
        records = []
        for _ in range(2):
            with tempfile.TemporaryDirectory() as work:
                run(binary, work, field, trace_time=2.0e-5)
                with open(os.path.join(work, "times_lost.dat"), "rb") as f:
                    records.append(f.read())
        assert records[0] == records[1], "fixed-seed tracing changed between runs"
    assert not failures, failures
    print("SPECTRE config gate, classification and deterministic tracing PASS")


if __name__ == "__main__":
    main()
