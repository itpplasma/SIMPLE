#!/usr/bin/env python3
"""Evidence for itpplasma/SIMPLE#398: the symplectic step's Newton Jacobian
consumes the Boozer field's second radial derivatives, which diverge at the
axis in s = rho^2 coordinates; RK never uses them. Three figures:

  fig1_d2_divergence.png  - d2Bmod/ds2, d2hth/ds2 vs s (internal vs booz_xform):
                            both diverge ~ s^(-3/2); the symplectic Jacobian uses them.
  fig2_orbit_energy.png   - one axis-grazing orbit (idx 12): radius s(t) and energy
                            |dH/H|(t) for RK vs symplectic Euler1/midpoint.
  fig3_values_match.png   - iota(s) and mirror ratio Bmax/Bmin(s) agree to <0.2%
                            between the two fields: it is not the field values.

Inputs (paths as argv or defaults): internal field-derivative dump from
diag_field_deriv.x, booz_xform and SIMPLE chartmap .nc, and idx-12 orbit traces
from diag_traj.x (integmode 0/1/3).
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.interpolate import CubicSpline

FD_INT = sys.argv[1] if len(sys.argv) > 1 else "/tmp/fd_int/field_internal.dat"
BX_NC = sys.argv[2] if len(sys.argv) > 2 else "boozer_chartmap_w7x.nc"
SM_NC = sys.argv[3] if len(sys.argv) > 3 else "122-chartmap_simple_axispow_t1s/chartmap_simple_axispow.nc"
ORBDIR = sys.argv[4] if len(sys.argv) > 4 else "/tmp/traj_diag"
THETA0, PHI0 = 0.3, 0.2


def bx_bmod_line(nc, theta0, phi0):
    """Bmod(s) at the angle closest to (theta0,phi0) for a chartmap .nc."""
    d = Dataset(nc)
    rho = d["rho"][:]
    B = d["Bmod"][:]
    th = d["theta"][:] if "theta" in d.variables else np.linspace(0, 2*np.pi, B.shape[0], endpoint=False)
    ze = d["zeta"][:] if "zeta" in d.variables else np.linspace(0, 2*np.pi, B.shape[1], endpoint=False)
    tf = float(d.torflux)
    Aphi = d["A_phi"][:]
    s_a = d["s"][:] if "s" in d.variables else rho**2
    d.close()
    rax = [i for i, n in enumerate(B.shape) if n == len(rho)][0]
    B = np.moveaxis(B, rax, 0)            # (rho, A, B)
    ia = int(np.argmin(np.abs(th[:B.shape[1]] - theta0)))
    ib = int(np.argmin(np.abs(ze[:B.shape[2]] - phi0)))
    return rho**2, B[:, ia, ib], B, rho, tf, Aphi, s_a


# ---- fig1: second-derivative divergence ----
fd = np.loadtxt(FD_INT)
s_i = fd[:, 0]
d2Bmod_i, d2hth_i = fd[:, 10], fd[:, 11]
s_bx, Bline_bx, Bfull_bx, rho_bx, tf_bx, Aphi_bx, sa_bx = bx_bmod_line(BX_NC, THETA0, PHI0)
cs = CubicSpline(s_bx, Bline_bx)
d2Bmod_bx = np.array([cs(x, 2) for x in s_i])
# hth_bx = B_theta/Bmod; B_theta tiny near axis -> compare Bmod second deriv (dominant)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
ax[0].loglog(s_i, np.abs(d2Bmod_i), "b.-", label="internal transform")
ax[0].loglog(s_i, np.abs(d2Bmod_bx), "r.-", label="booz_xform chartmap")
sref = np.array([1e-6, 1e-2]); ax[0].loglog(sref, 5e2*sref**-1.5, "k--", lw=1, label=r"$\propto s^{-3/2}$")
ax[0].axvspan(s_i.min(), 1e-4, color="orange", alpha=0.15, label="orbit grazing band")
ax[0].set_xlabel("s"); ax[0].set_ylabel(r"$|d^2 B_{mod}/ds^2|$")
ax[0].set_title("Second radial derivative of |B| diverges at axis"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, which="both")
ax[1].loglog(s_i, np.abs(d2hth_i), "b.-", label="internal $d^2 h_\\theta/ds^2$")
ax[1].axvspan(s_i.min(), 1e-4, color="orange", alpha=0.15)
ax[1].set_xlabel("s"); ax[1].set_ylabel(r"$|d^2 h_\theta/ds^2|$")
ax[1].set_title("These enter the symplectic Newton Jacobian (RK never uses them)")
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, which="both")
fig.tight_layout(); fig.savefig("fig1_d2_divergence.png", dpi=140); print("fig1 saved")

# ---- fig2: orbit radius and energy ----
fig, ax = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)
styles = {0: ("k", "RK (relerr 1e-12)"), 1: ("tab:blue", "symplectic Euler1"), 3: ("tab:red", "symplectic midpoint")}
for im, (c, lab) in styles.items():
    try:
        d = np.genfromtxt(f"{ORBDIR}/orbfig_im{im}.dat"); d = d[np.isfinite(d[:, 1])]
    except Exception:
        continue
    t, s, H = d[:, 0], d[:, 1], d[:, 5]
    ax[0].plot(t, s, c=c, lw=1.0, label=lab)
    ax[1].plot(t, np.abs(H - H[0]) / abs(H[0]), c=c, lw=1.2, label=lab)
ax[0].set_ylabel("s (flux radius)"); ax[0].set_yscale("log")
ax[0].set_title("Axis-grazing trapped orbit (idx 12): RK stays bounded, symplectic diffuses"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
ax[1].set_ylabel(r"$|\Delta H/H|$"); ax[1].set_xlabel("t [s]")
ax[1].set_title("Energy: RK conserved (~1e-12); symplectic random-walks at each near-axis pass"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
fig.tight_layout(); fig.savefig("fig2_orbit_energy.png", dpi=140); print("fig2 saved")

# ---- fig3: field values agree ----
def profile(nc, aphi_s):
    s2, Bl, Bf, rho, tf, Aphi, s_a = bx_bmod_line(nc, THETA0, PHI0)
    rax0 = 0
    Bmin = Bf.min(axis=(1, 2)); Bmax = Bf.max(axis=(1, 2))
    dA = np.gradient(Aphi, s_a if aphi_s else rho**2)
    iota = -dA / tf
    return rho**2, Bmin, Bmax, (s_a if aphi_s else rho**2), iota
sB, BminB, BmaxB, saB, iotaB = profile(BX_NC, False)
sS, BminS, BmaxS, saS, iotaS = profile(SM_NC, True)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
ax[0].plot(saB, iotaB, "r-", label="booz_xform"); ax[0].plot(saS, iotaS, "b--", label="SIMPLE internal")
ax[0].set_xlim(0, 0.5); ax[0].set_ylim(0.84, 0.92); ax[0].set_xlabel("s"); ax[0].set_ylabel(r"$\iota$")
ax[0].set_title(r"$\iota(s)$ agree (axis value 0.856)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
ax[1].plot(sB, BmaxB/BminB, "r-", label="booz_xform"); ax[1].plot(sS, BmaxS/BminS, "b--", label="SIMPLE internal")
ax[1].set_xlim(0, 0.5); ax[1].set_xlabel("s"); ax[1].set_ylabel(r"$B_{max}/B_{min}$")
ax[1].set_title("Mirror ratio agrees to <0.2%: not the field values"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
fig.tight_layout(); fig.savefig("fig3_values_match.png", dpi=140); print("fig3 saved")
