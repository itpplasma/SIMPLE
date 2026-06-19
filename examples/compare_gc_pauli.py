#!/usr/bin/env python3
"""Banana-orbit cross-validation figures: GC vs CPP vs full orbit.

Plotting-only harness for the figures of the full-orbit / CPP cross-validation
suite. It runs the production trajectory dumper diag_traj_models.x (which traces
ONE trapped banana with the real symplectic GC integrator, the CPP pusher, and
the gyro-resolved Boris full orbit), reads its .dat output, and produces, per
field case:

  Panel A  R-Z poloidal projection of the banana, GC vs CPP vs full orbit
           overlaid. GC traces the smooth banana, CPP overlaps GC, the full
           orbit gyrates around it.
  Panel B  conservation vs time: H(t)/H(0)-1, p_phi(t)-p_phi(0), mu(t)/mu(0).

Field cases:
  analytic  circular-tokamak "test" canonical chart (GC/CPP) and the matching
            analytic cylindrical tokamak field (full orbit).
  vmec      BOOZER chart of the QA wout (GC/CPP). The full-orbit VMEC
            curvilinear path depends on libneo metric/Christoffel (PR #322) and
            is not wired here, so the vmec full-orbit curve is omitted.

PNGs are written to the output directory (default /tmp/crossval), four files:
  {analytic,vmec}_rz_overlay.png and {analytic,vmec}_conservation.png

Usage:
  python examples/compare_gc_pauli.py [--build-dir build] [--outdir /tmp/crossval]
                                      [--nstep-analytic 1500] [--nstep-vmec 2500]
                                      [--no-run]
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

MODELS = ("gc", "cpp", "fullorbit")
STYLE = {
    "gc":        dict(color="tab:blue",   lw=2.2, ls="-",  label="GC (symplectic)"),
    "cpp":       dict(color="tab:orange", lw=1.3, ls="--", label="CPP (Pauli)"),
    "fullorbit": dict(color="tab:green",  lw=0.6, ls="-",  alpha=0.7,
                      label="full orbit (Boris)"),
}


def load(outdir, case, model):
    path = os.path.join(outdir, f"{case}_{model}.dat")
    if not os.path.isfile(path):
        return None
    data = np.loadtxt(path)
    if data.ndim == 1 or data.shape[0] < 2:
        return None
    return data  # columns: t R Z Hrel pphi mu


def run_dumper(exe, case, outdir, nstep):
    cmd = [exe, case, outdir, str(nstep)]
    print("run:", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    sys.stdout.write(res.stdout)
    if res.returncode != 0:
        sys.stderr.write(res.stderr)
        print(f"  WARNING: dumper exited {res.returncode} for case {case}")


def plot_rz(outdir, case, title):
    fig, ax = plt.subplots(figsize=(5.2, 6.0))
    any_curve = False
    for model in MODELS:
        d = load(outdir, case, model)
        if d is None:
            continue
        ax.plot(d[:, 1], d[:, 2], **STYLE[model])
        any_curve = True
    if not any_curve:
        plt.close(fig)
        return None
    ax.set_xlabel("R")
    ax.set_ylabel("Z")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_title(f"{title}\nbanana poloidal projection")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = os.path.join(outdir, f"{case}_rz_overlay.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("wrote", out)
    return out


def plot_conservation(outdir, case, title):
    fig, axes = plt.subplots(3, 1, figsize=(7.0, 7.5), sharex=True)
    hax, pax, max_ = axes
    any_curve = False
    for model in MODELS:
        d = load(outdir, case, model)
        if d is None:
            continue
        t = d[:, 0]
        hax.plot(t, d[:, 3], **STYLE[model])
        # p_phi is canonical only for GC/CPP; full orbit reports 0 -> skip.
        if model != "fullorbit":
            pax.plot(t, d[:, 4] - d[0, 4], **STYLE[model])
        mu0 = d[0, 5]
        if mu0 != 0.0:
            max_.plot(t, d[:, 5] / mu0 - 1.0, **STYLE[model])
        any_curve = True
    if not any_curve:
        plt.close(fig)
        return None
    hax.set_ylabel(r"$H/H_0 - 1$")
    hax.set_title(f"{title}\ninvariant conservation vs time")
    pax.set_ylabel(r"$p_\varphi(t) - p_\varphi(0)$")
    max_.set_ylabel(r"$\mu/\mu_0 - 1$")
    max_.set_xlabel("t")
    for a in axes:
        a.grid(True, alpha=0.3)
        a.legend(loc="best", fontsize=8)
    fig.tight_layout()
    out = os.path.join(outdir, f"{case}_conservation.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("wrote", out)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-dir", default="build")
    ap.add_argument("--outdir", default="/tmp/crossval")
    ap.add_argument("--nstep-analytic", type=int, default=1500)
    ap.add_argument("--nstep-vmec", type=int, default=2500)
    ap.add_argument("--no-run", action="store_true",
                    help="reuse existing .dat files, do not run the dumper")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    exe = os.path.join(args.build_dir, "diag_traj_models.x")

    if not args.no_run:
        if not os.path.isfile(exe):
            sys.exit(f"dumper not found: {exe} (build with make CONFIG=Fast)")
        run_dumper(exe, "analytic", args.outdir, args.nstep_analytic)
        run_dumper(exe, "vmec", args.outdir, args.nstep_vmec)

    produced = []
    for case, title in (("analytic", "Analytic circular tokamak"),
                        ("vmec", "VMEC QA stellarator")):
        rz = plot_rz(args.outdir, case, title)
        co = plot_conservation(args.outdir, case, title)
        produced += [p for p in (rz, co) if p]

    print("\nfigures:")
    for p in produced:
        print(" ", p)


if __name__ == "__main__":
    main()
