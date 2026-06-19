#!/usr/bin/env python3
"""Regenerate the three-system cross-validation figures from diag_crossval_figures.x.

Run the Fortran dump first (from a directory containing wout.nc for the VMEC
panels):

    ./build/diag_crossval_figures.x <data_dir>
    python examples/plot_crossval_figures.py <data_dir> <out_dir>

Figures, all from real integrator output (no fabricated curves):
  banana_analytic.png      GC drift vs 6D Pauli (shared field), CP/CPP-sym/CPP-var
  conservation_analytic.png energy/mu/p_phi vs time + symplectic dt plateau
  banana_vmec.png          VMEC CP (gyro) vs CPP-sym (big-step) flux poloidal
  conservation_vmec.png    VMEC Hamiltonian energy vs time, CP and CPP-sym
"""
import sys
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(d, name):
    return np.loadtxt(os.path.join(d, name))


def banana_analytic(d, out):
    gc = load(d, "banana_gc_drift.dat")
    pa = load(d, "banana_pauli.dat")
    cp = load(d, "analytic_cp.dat")
    cs = load(d, "analytic_cpp_sym.dat")
    cv = load(d, "analytic_cpp_var.dat")

    fig, ax = plt.subplots(1, 2, figsize=(11, 5.2))

    # Shared trapped launch: the full 6D Pauli orbit wraps the smooth GC banana
    # with gyration of order rho*.
    ax[0].plot(pa[:, 1], pa[:, 2], ",", color="0.6", label="6D Pauli (full orbit)")
    ax[0].plot(gc[:, 1], gc[:, 2], "-", color="C3", lw=1.5, label="GC drift")
    ax[0].set_title("Shared-field banana: full orbit vs GC")
    ax[0].set_xlabel(r"$R-R_0$")
    ax[0].set_ylabel(r"$Z$")
    ax[0].set_aspect("equal")
    ax[0].legend(loc="upper right", fontsize=8)

    # Thesis canonical banana (r0=0.1, mu=1e-5): CP gyro-resolved, CPP-sym and
    # CPP-var at large steps trace the same trapped section. Decimate the dense
    # CP run so the banana stays legible.
    cpd = cp[:: max(1, len(cp) // 4000)]
    ax[1].plot(cpd[:, 4], cpd[:, 5], ".", ms=1.5, color="C0",
               label=r"CP $\Delta t{=}1$ (gyro-resolved)")
    ax[1].plot(cs[:, 4], cs[:, 5], "o", ms=2.5, color="C1",
               label=r"CPP-sym $\Delta t{=}80$")
    ax[1].plot(cv[:, 4], cv[:, 5], "x", ms=4, color="C2",
               label=r"CPP-var $\Delta t{=}800$")
    ax[1].set_title("Canonical-midpoint banana, increasing step")
    ax[1].set_xlabel(r"$R-R_0$")
    ax[1].set_ylabel(r"$Z$")
    ax[1].set_aspect("equal")
    ax[1].legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    p = os.path.join(out, "banana_analytic.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def conservation_analytic(d, out):
    cp = load(d, "analytic_cp.dat")
    cs = load(d, "analytic_cpp_sym.dat")
    sw = load(d, "analytic_cpp_sym_dtsweep.dat")

    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))

    # Relative energy: bounded oscillation, no secular drift (symplectic). Each
    # model spans a different physical time, so plot vs total elapsed time.
    tcp = cp[:, 0] * 1.0
    tcs = cs[:, 0] * 80.0
    ax[0, 0].plot(tcp, cp[:, 6], color="C0", lw=0.8, label=r"CP $\Delta t{=}1$")
    ax[0, 0].plot(tcs, cs[:, 6], color="C1", lw=0.8, label=r"CPP-sym $\Delta t{=}80$")
    ax[0, 0].set_xlabel("time")
    ax[0, 0].set_ylabel(r"$\Delta E / E_0$")
    ax[0, 0].set_title("Relative energy error (bounded, no drift)")
    ax[0, 0].legend(fontsize=8)

    # mu: CP keeps it to machine zero (mu is the radial gyration energy / |B|).
    ax[0, 1].plot(tcp, (cp[:, 7] - cp[0, 7]) / cp[0, 7], color="C0")
    ax[0, 1].set_xlabel("time")
    ax[0, 1].set_ylabel(r"$\Delta \mu / \mu_0$")
    ax[0, 1].set_title("Relative magnetic-moment error (CP)")
    ax[0, 1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    # Toroidal canonical momentum: axisymmetric field, conserved to machine zero.
    ax[1, 0].plot(tcp, cp[:, 8] - cp[0, 8], color="C0", label="CP")
    ax[1, 0].plot(tcs, cs[:, 8] - cs[0, 8], color="C1", label="CPP-sym")
    ax[1, 0].set_xlabel("time")
    ax[1, 0].set_ylabel(r"$p_\varphi - p_{\varphi,0}$")
    ax[1, 0].set_title(r"Toroidal canonical momentum $p_\varphi$ (machine zero)")
    ax[1, 0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax[1, 0].legend(fontsize=8)

    # Energy band vs dt: flat plateau is the symplectic signature.
    ax[1, 1].loglog(sw[:, 0], sw[:, 1], "o-", color="C1")
    ax[1, 1].set_xlabel(r"$\Delta t$")
    ax[1, 1].set_ylabel(r"max $|\Delta E / E_0|$")
    ax[1, 1].set_title("Energy band vs step (CPP-sym plateau)")
    ax[1, 1].grid(True, which="both", ls=":", alpha=0.5)

    fig.tight_layout()
    p = os.path.join(out, "conservation_analytic.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def banana_vmec(d, out):
    cp = load(d, "vmec_cp.dat")
    cs = load(d, "vmec_cpp_sym.dat")

    fig, ax = plt.subplots(1, 2, figsize=(11, 5.2))
    ax[0].plot(cp[:, 4], cp[:, 5], ".", ms=1.0, color="C0")
    ax[0].set_title(r"VMEC CP (gyro-resolved, small $\Delta t$)")
    ax[1].plot(cs[:, 4], cs[:, 5], "o", ms=2.5, color="C1")
    ax[1].set_title(r"VMEC CPP-sym (GC band, big $\Delta t$)")
    for a in ax:
        a.set_xlabel(r"$\sqrt{s}\,\cos\vartheta$")
        a.set_ylabel(r"$\sqrt{s}\,\sin\vartheta$")
        a.set_aspect("equal")
    fig.tight_layout()
    p = os.path.join(out, "banana_vmec.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def conservation_vmec(d, out):
    cp = load(d, "vmec_cp.dat")
    cs = load(d, "vmec_cpp_sym.dat")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].plot(cp[:, 0], cp[:, 6], color="C0", lw=0.8)
    ax[0].set_title(r"VMEC CP energy (gyro-resolved)")
    ax[1].plot(cs[:, 0], cs[:, 6], color="C1", lw=0.8)
    ax[1].set_title(r"VMEC CPP-sym energy (big step)")
    for a in ax:
        a.set_xlabel("step")
        a.set_ylabel(r"$\Delta E / E_0$")
    fig.tight_layout()
    p = os.path.join(out, "conservation_vmec.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "crossval_data"
    out = sys.argv[2] if len(sys.argv) > 2 else d
    os.makedirs(out, exist_ok=True)

    made = [banana_analytic(d, out), conservation_analytic(d, out)]
    if os.path.exists(os.path.join(d, "vmec_cp.dat")):
        made += [banana_vmec(d, out), conservation_vmec(d, out)]
    for p in made:
        print("wrote", p)


if __name__ == "__main__":
    main()
