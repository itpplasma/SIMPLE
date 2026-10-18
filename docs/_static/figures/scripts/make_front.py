#!/usr/bin/env python3
"""Generate docs/_static/figures/front.png from a real SIMPLE run.

Configuration: the public Landreman & Paul (2022) precise quasi-axisymmetric
equilibrium scaled to reactor size, as shipped in the simsopt test data and
used by the SIMPLE test suite. The script

1. downloads the VMEC file at a pinned simsopt commit and checks its SHA-256,
2. writes ``simple.in`` and runs ``simple.x`` (3.5 MeV alphas, CPU/OpenMP),
3. plots the guiding-centre loss positions on the last closed flux surface
   (LCFS) and the confined fraction versus time.

Usage (from the repository root, after ``make``)::

    python docs/_static/figures/scripts/make_front.py            # run + plot
    python docs/_static/figures/scripts/make_front.py --plot-only

The run directory defaults to ``build/figures/front``.
"""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import simple_style  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
SIMSOPT_COMMIT = "205ba8e5053032b05b91020be4592d00d1c15a52"
WOUT_NAME = "wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
WOUT_URL = (
    f"https://github.com/hiddenSymmetries/simsopt/raw/{SIMSOPT_COMMIT}"
    f"/tests/test_files/{WOUT_NAME}"
)
WOUT_SHA256 = "333ab8cef64c4b4f406e76209266fc3bf4d7976e9c28c162983c2120e112e771"

SIMPLE_IN = """\
&config
netcdffile = 'wout.nc'   ! public QA equilibrium (see make_front.py)
isw_field_type = 2       ! Boozer canonical coordinates
integmode = 1            ! symplectic Euler
trace_time = 5d-2        ! s
sbeg = 0.5d0             ! start surface, normalized toroidal flux
ntestpart = 3000
deterministic = .true.
output_results_netcdf = .true.
/
"""


def fetch_wout(path: Path) -> None:
    if not path.exists():
        print(f"downloading {WOUT_URL}")
        urllib.request.urlretrieve(WOUT_URL, path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != WOUT_SHA256:
        raise SystemExit(f"{path}: SHA-256 {digest} != {WOUT_SHA256}")


def run_simple(workdir: Path, exe: Path) -> None:
    workdir.mkdir(parents=True, exist_ok=True)
    fetch_wout(workdir / "wout.nc")
    (workdir / "simple.in").write_text(SIMPLE_IN)
    with open(workdir / "simple.log", "w") as log:
        subprocess.run(
            [str(exe), "simple.in"], cwd=workdir, stdout=log, stderr=log, check=True
        )


def lcfs(wout_path: Path, ntheta: int, nphi: int):
    """Return a full-torus angle grid and VMEC LCFS Fourier coefficients."""
    import netCDF4

    with netCDF4.Dataset(wout_path) as w:
        if int(w["lasym__logical__"][:]) != 0:
            raise SystemExit("stellarator-asymmetric wout not supported here")
        xm = np.asarray(w["xm"][:])
        xn = np.asarray(w["xn"][:])
        rmnc = np.asarray(w["rmnc"][-1, :])
        zmns = np.asarray(w["zmns"][-1, :])
    theta = np.linspace(0.0, 2.0 * np.pi, ntheta)
    phi = np.linspace(0.0, 2.0 * np.pi, nphi)
    return theta, phi, xm, xn, rmnc, zmns


def eval_rz(theta, phi, xm, xn, rmnc, zmns):
    arg = np.multiply.outer(theta, xm) - np.multiply.outer(phi, xn)
    return np.cos(arg) @ rmnc, np.sin(arg) @ zmns


def load_results(workdir: Path):
    import netCDF4

    with netCDF4.Dataset(workdir / "results.nc") as d:
        res = {
            k: np.asarray(d[k][:])
            for k in ("zend", "xend_cart", "times_lost", "orbit_exit_code", "trap_par")
        }
        res["version"] = d.getncattr("simple_version")
    cf = np.loadtxt(workdir / "confined_fraction.dat")
    return res, cf


def view_matrix(azim_deg: float, elev_deg: float) -> np.ndarray:
    a, e = np.radians(azim_deg), np.radians(elev_deg)
    rz = np.array([[np.cos(a), np.sin(a), 0], [-np.sin(a), np.cos(a), 0], [0, 0, 1]])
    rx = np.array([[1, 0, 0], [0, np.sin(e), np.cos(e)], [0, -np.cos(e), np.sin(e)]])
    # Rows: screen x, screen y, depth (towards the viewer is positive).
    return rx @ rz


_ang = np.linspace(0.0, 2.0 * np.pi, 13)
MARKER_SHAPES = {
    # Same visual area: disc for trapped, diamond for passing.
    "trapped": np.stack([np.cos(_ang), np.sin(_ang)], axis=-1),
    "passing": 1.25 * np.array([[0, 1], [1, 0], [0, -1], [-1, 0], [0, 1.0]]),
}


def confined_curves(times_lost, exit_code, trap_par, t):
    """Confined fractions on times ``t`` from the per-marker loss times."""
    lost = exit_code == 1
    if np.any(exit_code >= 101):
        raise SystemExit("numerically failed markers present; curves need care")
    n = exit_code.size
    passing = trap_par < 0
    curves = {}
    for key, sel in (("passing", passing), ("trapped", ~passing)):
        tl = np.sort(times_lost[lost & sel])
        curves[key] = (np.count_nonzero(sel) - np.searchsorted(tl, t, side="right")) / n
    curves["total"] = curves["passing"] + curves["trapped"]
    return curves


def draw_torus(ax, wout_path: Path, theta_l, phi_l, trapped):
    """Painter's-algorithm render of the LCFS with loss markers on it."""
    ntheta, nphi = 64, 361
    theta, phi, xm, xn, rmnc, zmns = lcfs(wout_path, ntheta, nphi)
    T, P = np.meshgrid(theta, phi, indexing="ij")
    R, Z = eval_rz(T.ravel(), P.ravel(), xm, xn, rmnc, zmns)
    R, Z = R.reshape(T.shape), Z.reshape(T.shape)
    X, Y = R * np.cos(P), R * np.sin(P)
    M = view_matrix(azim_deg=-70.0, elev_deg=24.0)
    pts = np.stack([X, Y, Z], axis=-1) @ M.T

    quads = np.stack(
        [pts[:-1, :-1], pts[1:, :-1], pts[1:, 1:], pts[:-1, 1:]], axis=2
    )  # (nt-1, np-1, 4, 3)
    quads = quads.reshape(-1, 4, 3)
    e1 = quads[:, 2] - quads[:, 0]
    e2 = quads[:, 3] - quads[:, 1]
    normal = np.cross(e1, e2)
    normal /= np.linalg.norm(normal, axis=1, keepdims=True)
    light = np.array([-0.4, 0.55, 0.75])
    light /= np.linalg.norm(light)
    lam = np.abs(normal @ light)
    base = np.array(matplotlib.colors.to_rgb(simple_style.COLORS["surface"]))
    shade = 0.55 + 0.45 * lam
    face = np.clip(base[None, :] * shade[:, None] + 0.08 * lam[:, None] ** 8, 0, 1)
    depth = quads[:, :, 2].mean(axis=1)

    # Markers: small screen-space discs at the loss points, slightly in front.
    Rl, Zl = eval_rz(theta_l, phi_l, xm, xn, rmnc, zmns)
    pl = np.stack([Rl * np.cos(phi_l), Rl * np.sin(phi_l), Zl], axis=-1) @ M.T
    rad = 0.0068 * float(np.ptp(pts[..., 0]))
    marker_polys = [
        p[:2] + rad * MARKER_SHAPES["trapped" if t else "passing"]
        for p, t in zip(pl, trapped)
    ]
    marker_depth = pl[:, 2] + 0.02 * float(np.ptp(pts[..., 2]))
    rgba = {k: matplotlib.colors.to_rgba(simple_style.COLORS[k])
            for k in ("trapped", "passing")}
    marker_face = np.array([rgba["trapped"] if t else rgba["passing"] for t in trapped])

    polys = list(quads[:, :, :2]) + list(marker_polys)
    faces = np.vstack([np.column_stack([face, np.ones(len(face))]), marker_face])
    edges = np.vstack(
        [np.column_stack([face, np.ones(len(face))]),
         np.tile([0.1, 0.1, 0.1, 1.0], (len(marker_polys), 1))]
    )
    widths = np.concatenate([np.full(len(face), 0.3), np.full(len(marker_polys), 0.4)])
    order = np.argsort(np.concatenate([depth, marker_depth]))
    coll = PolyCollection(
        [polys[i] for i in order],
        facecolors=faces[order],
        edgecolors=edges[order],
        linewidths=widths[order],
    )
    ax.add_collection(coll)
    ax.set_xlim(pts[..., 0].min() - 0.3, pts[..., 0].max() + 0.3)
    ax.set_ylim(pts[..., 1].min() - 0.3, pts[..., 1].max() + 0.3)
    ax.set_aspect("equal")
    ax.axis("off")


def plot(workdir: Path, out: Path) -> None:
    simple_style.use()
    res, cf = load_results(workdir)
    lost = res["orbit_exit_code"] == 1
    zl = res["zend"][lost]
    if not np.allclose(zl[:, 0], 1.0, atol=1e-6):
        raise SystemExit("lost markers do not end on s = 1")
    theta_l, phi_l = zl[:, 1], zl[:, 2]

    # Consistency check: VMEC LCFS at (theta, phi) must reproduce xend_cart.
    _, _, xm, xn, rmnc, zmns = lcfs(workdir / "wout.nc", 2, 2)
    Rl, Zl = eval_rz(theta_l, phi_l, xm, xn, rmnc, zmns)
    xyz = np.stack([Rl * np.cos(phi_l), Rl * np.sin(phi_l), Zl], axis=-1)
    err = np.abs(xyz * 100.0 - res["xend_cart"][lost]).max()
    print(f"max |LCFS(theta, phi) - xend_cart| = {err:.3g} cm")
    if err > 1.0:
        raise SystemExit("loss positions are inconsistent with xend_cart")

    t_end = cf[-1, 0]
    # Check the event-resolved curves against SIMPLE's confined_fraction.dat.
    grid = confined_curves(res["times_lost"], res["orbit_exit_code"],
                           res["trap_par"], cf[:, 0])
    dev = max(np.abs(grid["passing"] - cf[:, 1]).max(),
              np.abs(grid["trapped"] - cf[:, 2]).max())
    print(f"max deviation from confined_fraction.dat = {dev:.3g}")
    if dev > 1.5 / lost.size:
        raise SystemExit("event-resolved curves disagree with confined_fraction.dat")
    trapped = res["trap_par"][lost] >= 0

    fig = plt.figure(figsize=(10.0, 4.6))
    ax3 = fig.add_axes([0.0, 0.10, 0.60, 0.77])
    draw_torus(ax3, workdir / "wout.nc", theta_l, phi_l, trapped)
    fig.text(0.015, 0.975, "Alpha-particle loss positions on the last closed flux surface",
             fontsize=12, fontweight="bold", va="top")
    fig.text(0.015, 0.915,
             "Landreman-Paul QA, reactor scale  ·  3.5 MeV alphas started at s = 0.5  "
             "·  guiding centre, collisionless",
             fontsize=9, va="top", color=simple_style.COLORS["neutral"])
    handles = []
    for key, mk in (("trapped", "o"), ("passing", "D")):
        n_k = int(np.count_nonzero(trapped if key == "trapped" else ~trapped))
        handles.append(plt.Line2D([], [], ls="none", marker=mk, ms=6 if mk == "o" else 5,
                                  mfc=simple_style.COLORS[key], mec="#1A1A1A", mew=0.5,
                                  label=f"lost, {key} at start ({n_k})"))
    ax3.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.03, -0.1),
               ncol=2, handletextpad=0.3, columnspacing=1.2)

    ax = fig.add_axes([0.67, 0.17, 0.27, 0.66])
    # Include every loss event so the step locations retain their actual times.
    loss_times = res["times_lost"][lost]
    t_start = min(1e-6, loss_times.min()) if loss_times.size else 1e-6
    t = np.unique(np.concatenate(([t_start, t_end], loss_times)))
    curves = confined_curves(res["times_lost"], res["orbit_exit_code"],
                             res["trap_par"], t)
    for key in ("total", "passing", "trapped"):
        y = curves[key]
        ax.step(t, y, where="post", color=simple_style.COLORS[key],
                ls=simple_style.LINESTYLES[key])
        ax.annotate(key, (t[-1], y[-1]), xytext=(4, 0), textcoords="offset points",
                    va="center", fontsize=9, color=simple_style.COLORS[key],
                    annotation_clip=False)
    ax.set_xscale("log")
    ax.set_xlim(t[0], t_end)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("confined fraction of all markers")
    ax.set_title(f"{lost.sum()} of {lost.size} markers lost by {t_end * 1e3:g} ms",
                 fontsize=10, loc="left")
    fig.text(0.995, 0.005, f"SIMPLE {res['version']}", ha="right", va="bottom",
             fontsize=7, color=simple_style.COLORS["neutral"])
    simple_style.save(fig, out, dpi=130)
    print(f"wrote {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--simple-exe", type=Path, default=REPO / "build" / "simple.x")
    parser.add_argument("--workdir", type=Path, default=REPO / "build" / "figures" / "front")
    parser.add_argument("--out", type=Path,
                        default=REPO / "docs" / "_static" / "figures" / "front.png")
    parser.add_argument("--plot-only", action="store_true",
                        help="reuse an existing run in --workdir")
    args = parser.parse_args()
    if not args.plot_only:
        run_simple(args.workdir.resolve(), args.simple_exe.resolve())
    plot(args.workdir.resolve(), args.out)


if __name__ == "__main__":
    main()
