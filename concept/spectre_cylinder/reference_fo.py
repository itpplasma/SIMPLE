"""Independent full-orbit reference in the analytic Lundquist cylinder.

Straight periodic cylinder (Cartesian X, Y in the cross-section, s along the
axis), field of the volume that owns the particle, DOP853 with the interface
r = r1 as an exact terminal event that switches the owner. Units: lengths in
m, |v| = 1, time as path length in m, dv/dt = v x B / rho0 with
rho0 = m v0 / q in T m.

For each SIMPLE marker of cylinder_markers.x two reference orbits are run:
  independent: same guiding-centre distance from the sheet, side and pitch,
               NPHASE uniform random gyrophases, plain Larmor lift;
  paired:      SIMPLE's entry particle mapped to the cylinder.
Both are sampled on SIMPLE's microstep grid with the same far-side recording
rule. Usage: reference_fo.py MARKERS.dat FIXTURE.h5 OUT.npz
"""

import sys

import h5py
import numpy as np
from scipy.integrate import solve_ivp

from lundquist import from_fluxes, spectre_fluxes

NPER, NGYRO, NPHASE = 32, 3, 4
RTOL, ATOL = 1e-11, 1e-13


def load_field(h5file):
    with h5py.File(h5file, "r") as f:
        ph = f["input/physics"]
        phiedge = float(ph["phiedge"][()].ravel()[0])
        tflux, pflux, mu = (np.ravel(ph[k][()]) for k in ("tflux", "pflux", "mu"))
        im = np.ravel(f["output/im"][()])
        rbc = np.atleast_2d(f["output/Rbc"][()])
    rbc = rbc if rbc.shape[-1] == len(im) else rbc.T
    m0, m1 = np.flatnonzero(im == 0)[0], np.flatnonzero(im == 1)[0]
    r0, a, r1 = rbc[-1, m0], rbc[-1, m1], rbc[1, m1]
    tor, pol = spectre_fluxes(phiedge, tflux, pflux)
    return from_fluxes([r1], a, mu, tor, pol, 2*np.pi*r0), r0, r1


class Orbit:
    def __init__(self, field, r1, rho0):
        self.f, self.r1, self.rho0 = field, r1, rho0

    def bvec(self, owner, x, y):
        r = np.hypot(x, y)
        bt, bs = self.f.field_in(owner, r)
        return np.array([-bt*y/r, bt*x/r, bs])

    def rhs(self, owner):
        def f(t, z):
            b = self.bvec(owner, z[0], z[1])
            v = z[3:]
            return np.concatenate([v, np.cross(v, b)/self.rho0])
        return f

    def event(self, owner):
        def g(t, z):
            return np.hypot(z[0], z[1]) - self.r1
        g.terminal = True
        g.direction = 1 if owner == 0 else -1
        return g

    def run(self, z0, owner0, dt, nstep):
        """First grid sample in the other volume: (k, z) or (0, None)."""
        t_grid = dt*np.arange(1, nstep + 1)
        t, z, owner = 0.0, np.array(z0, float), owner0
        b_scale = np.linalg.norm(self.bvec(owner, z[0], z[1]))
        hmax = 2*np.pi*self.rho0/b_scale/64
        while t < t_grid[-1]:
            sol = solve_ivp(self.rhs(owner), (t, t_grid[-1]), z, method="DOP853",
                            rtol=RTOL, atol=ATOL, max_step=hmax,
                            events=self.event(owner), dense_output=True)
            t_end = sol.t[-1]
            if owner != owner0:
                inside = (t_grid > t) & (t_grid <= t_end)
                if inside.any():
                    k = int(np.flatnonzero(inside)[0])
                    return k + 1, sol.sol(t_grid[k]), owner
            if sol.status != 1:
                return 0, None, owner
            t, z = t_end, sol.y_events[0][0]
            owner = 1 - owner
        return 0, None, owner

    def far_side(self, z, owner):
        """mu (|v|^2 - v_par^2)/B and guiding-centre radius at state z."""
        b = self.bvec(owner, z[0], z[1])
        bmod = np.linalg.norm(b)
        b /= bmod
        v = z[3:]
        vpar = v @ b
        xgc = z[:3] - self.rho0/bmod*np.cross(b, v)
        return (v @ v - vpar**2)/bmod, np.hypot(xgc[0], xgc[1])

    def lift(self, r_gc, owner, lam, alpha):
        b = self.bvec(owner, r_gc, 0.0)
        bmod = np.linalg.norm(b)
        b /= bmod
        e1 = np.array([1.0, 0.0, 0.0])
        e2 = np.cross(b, e1)
        vperp = np.sqrt(1 - lam**2)
        vp = vperp*(np.cos(alpha)*e1 + np.sin(alpha)*e2)
        x = np.array([r_gc, 0.0, 0.0]) + self.rho0/bmod*np.cross(b, vp)
        return np.concatenate([x, lam*b + vp]), bmod


def torus_to_cylinder(x_cm, v, r0, r_cyl):
    """Particle from the torus to the straight cylinder at radius r_cyl.

    Velocity components are taken in the local (radial, poloidal, toroidal)
    frame about the torus centre line; r_cyl keeps the distance to the sheet.
    """
    x = np.asarray(x_cm)/100.0
    big_r = np.hypot(x[0], x[1])
    e_r = np.array([x[0], x[1], 0.0])/big_r
    e_phi = np.array([-x[1], x[0], 0.0])/big_r
    dr = np.array([big_r - r0, x[2]])
    r = np.hypot(*dr)
    rhat = (dr[0]*e_r + dr[1]*np.array([0, 0, 1.0]))/r
    that = np.cross(e_phi, rhat)
    v = np.asarray(v)/np.sqrt(2.0)
    vr, vt, vs = v @ rhat, v @ that, v @ e_phi
    return np.array([r_cyl, 0.0, 0.0, vr, vt, vs])


def read_markers(path):
    head = open(path).readline().split()
    meta = dict(zip(head[1::2], map(float, head[2::2])))
    d = np.loadtxt(path)
    return meta, d


def main():
    markers, h5file, out = sys.argv[1:4]
    meta, d = read_markers(markers)
    field, r0, r1 = load_field(h5file)
    rho0 = meta["ro0_cmG"]*1e-6
    dt = meta["dt_cm"]/100.0
    orb = Orbit(field, r1, rho0)
    rng = np.random.default_rng(20261107)
    nstep = NGYRO*NPER
    ind, pair = [], []
    for row in d:
        side, lam, d_gc = int(row[1]) - 1, row[2], row[7]
        for alpha in rng.uniform(0, 2*np.pi, NPHASE):
            z0, bmod = orb.lift(r1 + d_gc, side, lam, alpha)
            mu0 = (1 - lam**2)/bmod
            owner0 = int(np.hypot(z0[0], z0[1]) > r1)
            k, z, owner = orb.run(z0, owner0, dt, nstep)
            if k:
                mu1, r_after = orb.far_side(z, owner)
                ind.append([d_gc, lam, 1, mu1/mu0, r_after - r1])
            else:
                ind.append([d_gc, lam, 0, np.nan, np.nan])
        z0 = torus_to_cylinder(row[10:13], row[13:16], r0, r1 + row[23])
        owner0 = int(row[16]) - 1
        bt, bs = field.field_in(side, r1 + d_gc)
        mu0 = (1 - lam**2)/np.hypot(bt, bs)
        k, z, owner = orb.run(z0, owner0, dt, nstep)
        if k:
            mu1, r_after = orb.far_side(z, owner)
            pair.append([1, mu1/mu0, r_after - r1, k])
        else:
            pair.append([0, np.nan, np.nan, 0])
    np.savez(out, independent=np.array(ind), paired=np.array(pair),
             r1=r1, rho0=rho0, dt=dt)


if __name__ == "__main__":
    main()
