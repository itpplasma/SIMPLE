"""Two-volume Lundquist (constant-mu Beltrami) field of a straight cylinder.

Right-handed cylinder coordinates (r, theta, s): B = B_theta e_theta + B_s e_s,
curl B = mu_k B in volume k. With k = |mu| and sigma = sign(mu),

    B_s = c_J J0(k r) + c_Y Y0(k r),  B_theta = sigma (c_J J1(k r) + c_Y Y1(k r)),

and c_Y = 0 in the core volume. Lengths in m, B in T, fluxes in Wb. The
toroidal flux of a volume is int B_s 2 pi r dr; its poloidal flux is the flux
of B_theta through a ribbon theta = const of axial length L = 2 pi R0.

SPECTRE conventions (igeometry 3, boundary R = R0 + a cos t, Z = -a sin t, so
that (rho, t, phi) is right-handed and t is the cylinder angle about e_phi):
the toroidal flux of volume k is phiedge (tflux_k - tflux_{k-1}) and its
poloidal flux in the +theta sense is phiedge (pflux_k - pflux_{k-1}). The sign
of the latter was fixed empirically: only it gives the designed pressure
balance in SPECTRE's interface force.
"""

import numpy as np
from scipy.special import j0, j1, y0, y1


class Lundquist:
    def __init__(self, radii, mu, cj, cy):
        self.radii = np.asarray(radii, float)  # outer radius of each volume
        self.mu = np.asarray(mu, float)
        self.cj = np.asarray(cj, float)
        self.cy = np.asarray(cy, float)

    def volume(self, r):
        return np.minimum(np.searchsorted(self.radii, r, side="left"),
                          len(self.radii) - 1)

    def field_in(self, k, r):
        """(B_theta, B_s) of volume k's solution at radius r."""
        r = np.asarray(r, float)
        mu = self.mu[k]
        if mu == 0.0:
            return self.cy[k]/r, self.cj[k] + 0*r
        q = abs(mu)*r
        bs = self.cj[k]*j0(q)
        bt = self.cj[k]*j1(q)
        if self.cy[k] != 0.0:
            bs = bs + self.cy[k]*y0(q)
            bt = bt + self.cy[k]*y1(q)
        return np.sign(mu)*bt, bs

    def field(self, r):
        r = np.asarray(r, float)
        k = self.volume(r)
        bt = np.empty_like(r)
        bs = np.empty_like(r)
        for kk in range(len(self.radii)):
            m = k == kk
            bt[m], bs[m] = self.field_in(kk, r[m])
        return bt, bs

    def vector_potential_in(self, k, r):
        """(A_theta, A_s) with B = curl A in volume k; A = B/mu up to gauge."""
        bt, bs = self.field_in(k, r)
        return bt/self.mu[k], bs/self.mu[k]


def flux_rows(mu, r_in, r_out, length):
    """Toroidal and poloidal flux of the J and Y solutions on [r_in, r_out]."""
    k = abs(mu)
    sg = np.sign(mu)
    tor_j = 2*np.pi*(r_out*j1(k*r_out) - r_in*j1(k*r_in))/k
    pol_j = sg*length*(j0(k*r_in) - j0(k*r_out))/k
    if r_in == 0.0:
        return np.array([tor_j]), np.array([pol_j])
    tor_y = 2*np.pi*(r_out*y1(k*r_out) - r_in*y1(k*r_in))/k
    pol_y = sg*length*(y0(k*r_in) - y0(k*r_out))/k
    return np.array([tor_j, tor_y]), np.array([pol_j, pol_y])


def from_fluxes(r_iface, a, mu, tor_flux, pol_flux, length):
    """Lundquist field with the given per-volume fluxes (core: toroidal only)."""
    radii = list(r_iface) + [a]
    cj, cy = [], []
    for k, m in enumerate(mu):
        r_in = 0.0 if k == 0 else radii[k - 1]
        tor, pol = flux_rows(m, r_in, radii[k], length)
        if k == 0:
            cj.append(tor_flux[0]/tor[0])
            cy.append(0.0)
        else:
            c = np.linalg.solve(np.array([tor, pol]),
                                np.array([tor_flux[k], pol_flux[k]]))
            cj.append(c[0])
            cy.append(c[1])
    return Lundquist(radii, mu, cj, cy)


def spectre_fluxes(phiedge, tflux, pflux):
    """Physical per-volume (toroidal, +theta poloidal) fluxes from SPECTRE input."""
    t = phiedge*np.diff(np.concatenate([[0.0], tflux]))
    p = phiedge*np.diff(np.concatenate([[0.0], pflux]))
    return t, p


def pitch(bt, bs):
    return np.arctan2(bt, bs)
