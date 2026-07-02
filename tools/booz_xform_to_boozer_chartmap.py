#!/usr/bin/env python3
"""Convert booz_xform output (boozmn*.nc) to a Boozer chartmap NetCDF file.

Reads the Boozer-coordinate Fourier harmonics produced by booz_xform and
writes the extended chartmap format that SIMPLE reads without any VMEC or
booz_xform library at runtime. The boozmn file alone is sufficient: the
toroidal flux is taken from phi_b when present and otherwise recovered from
the surface geometry via the Boozer identity sqrt(g) B^2 = psi' (G + iota I).

Units: booz_xform is SI (T, m), the chartmap stores Gauss and cm.

Radial staggering: bmnc_b/rmnc_b/zmns_b/pmns_b live on the half grid
(pack_rad, surfaces jlist), iota_b/buco_b/bvco_b/phi_b on the full grid
(ns_b surfaces, index jlist[k]-1 holds the matching half-grid value for
surface k).

Half-grid s: s_half[k] = (jlist[k] - 1.5) / (ns_b - 1), the midpoint of the
interval between full-grid surfaces jlist[k]-1 and jlist[k] (1-based); this
matches booz_xform's own s_in for jlist[0] = 2.

Usage:
    python tools/booz_xform_to_boozer_chartmap.py <boozmn.nc> <output.nc>
        [--nrho N] [--ntheta N] [--nzeta N]
"""

import argparse
import sys

import numpy as np

try:
    import netCDF4
except ImportError:
    print("ERROR: netCDF4 required. pip install netCDF4")
    sys.exit(1)

try:
    from scipy.interpolate import CubicSpline
except ImportError:
    print("ERROR: scipy required. pip install scipy")
    sys.exit(1)

TWOPI = 2.0 * np.pi


def read_boozmn(filename):
    """Read a boozmn file into a plain dict of numpy arrays."""
    d = {}
    with netCDF4.Dataset(filename) as ds:
        def var(n):
            return np.asarray(ds.variables[n][:])

        d["ns"] = int(var("ns_b"))
        d["nfp"] = int(var("nfp_b"))
        d["lasym"] = bool(var("lasym__logical__"))
        d["jlist"] = var("jlist").astype(int)
        d["ixm"] = var("ixm_b").astype(int)
        d["ixn"] = var("ixn_b").astype(int)
        d["iota"] = var("iota_b")
        d["buco"] = var("buco_b")
        d["bvco"] = var("bvco_b")
        d["phi"] = var("phi_b")
        d["bmnc"] = var("bmnc_b")
        d["rmnc"] = var("rmnc_b")
        d["zmns"] = var("zmns_b")
        d["pmns"] = var("pmns_b")
        if d["lasym"]:
            d["bmns"] = var("bmns_b")
            d["rmns"] = var("rmns_b")
            d["zmnc"] = var("zmnc_b")
            d["pmnc"] = var("pmnc_b")
    # VMEC half-grid convention: s at the midpoint between full-grid surfaces
    # jlist[k]-1 and jlist[k] (1-based) is (jlist[k] - 1.5) / (ns_b - 1).
    d["s_half"] = (d["jlist"] - 1.5) / (d["ns"] - 1)
    return d


def interp_coeffs(coeffs, ixm, rho_half, rho_out):
    """Interpolate packed Fourier coefficients from the half grid to rho_out.

    Cubic spline in rho per mode; below the innermost half-grid surface the
    poloidal mode number sets the near-axis behaviour c ~ rho^m, so modes are
    continued with a power law instead of the spline's polynomial tail.
    """
    out = CubicSpline(rho_half, coeffs, axis=0)(rho_out)
    inner = rho_out < rho_half[0]
    if np.any(inner):
        ratio = rho_out[inner, None] / rho_half[0]
        m = np.minimum(ixm[None, :], 50)
        out[inner, :] = coeffs[0][None, :] * ratio**m
    return out


def fourier_eval(coeffs_grid, ixm, ixn, theta, zeta, parity):
    """Sum coeffs(rho, mn) * trig(m*theta - n*zeta) on a tensor angle grid.

    parity is 'cos' or 'sin'. Returns an (n_rho, n_theta, n_zeta) array.
    """
    angle = (ixm[:, None, None] * theta[None, :, None]
             - ixn[:, None, None] * zeta[None, None, :])
    basis = np.cos(angle) if parity == "cos" else np.sin(angle)
    return np.tensordot(coeffs_grid, basis, axes=([1], [0]))


def derive_psi_prime(d):
    """Recover psi' = Phi'(s)/(2*pi) from geometry when phi_b is absent.

    On each surface sqrt(g) B^2 = psi' (G + iota I) with sqrt(g) the Jacobian
    of (s, theta_B, zeta_B), so the angle average of J B^2 / (G + iota I)
    gives psi'. J comes from spectral angle derivatives plus radial spline
    derivatives of the geometry harmonics. s is the normalized toroidal
    flux, so psi' is one constant (carrying the VMEC handedness sign, equal
    to -phi_wout(ns)/(2*pi)); the spread over mid-radius surfaces, where the
    radial splines are accurate, is a consistency check.
    """
    ixm, ixn = d["ixm"], d["ixn"]
    s_half = d["s_half"]
    mid = np.where((s_half > 0.25) & (s_half < 0.95))[0]
    if len(mid) < 3:
        mid = np.arange(1, len(s_half) - 1)
    nth, nze = 73, 73
    theta = np.linspace(0.0, TWOPI, nth, endpoint=False)
    zeta = np.linspace(0.0, TWOPI / d["nfp"], nze, endpoint=False)
    angle = (ixm[:, None, None] * theta[None, :, None]
             - ixn[:, None, None] * zeta[None, None, :])
    cosg, sing = np.cos(angle), np.sin(angle)

    def ev(c, basis):
        return np.tensordot(c, basis, axes=([0], [0]))

    names = ["rmnc", "zmns", "pmns"]
    if d["lasym"]:
        names += ["rmns", "zmnc", "pmnc"]
    ds_coeffs = {name: CubicSpline(s_half, d[name], axis=0)(s_half[mid], nu=1)
                 for name in names}

    psi_primes = []
    for i, k in enumerate(mid):
        R = ev(d["rmnc"][k], cosg)
        B = ev(d["bmnc"][k], cosg)
        R_s = ev(ds_coeffs["rmnc"][i], cosg)
        Z_s = ev(ds_coeffs["zmns"][i], sing)
        phi_s = ev(ds_coeffs["pmns"][i], sing)
        R_t = ev(-ixm * d["rmnc"][k], sing)
        Z_t = ev(ixm * d["zmns"][k], cosg)
        phi_t = ev(ixm * d["pmns"][k], cosg)
        R_z = ev(ixn * d["rmnc"][k], sing)
        Z_z = ev(-ixn * d["zmns"][k], cosg)
        phi_z = 1.0 + ev(-ixn * d["pmns"][k], cosg)
        if d["lasym"]:
            R += ev(d["rmns"][k], sing)
            B += ev(d["bmns"][k], sing)
            R_s += ev(ds_coeffs["rmns"][i], sing)
            Z_s += ev(ds_coeffs["zmnc"][i], cosg)
            phi_s += ev(ds_coeffs["pmnc"][i], cosg)
            R_t += ev(ixm * d["rmns"][k], cosg)
            Z_t += ev(-ixm * d["zmnc"][k], sing)
            phi_t += ev(-ixm * d["pmnc"][k], sing)
            R_z += ev(-ixn * d["rmns"][k], cosg)
            Z_z += ev(ixn * d["zmnc"][k], sing)
            phi_z += ev(ixn * d["pmnc"][k], sing)

        # Jacobian in cylindrical components (dR, R dphi, dZ) of the three
        # tangent vectors; the e_phi factor contributes the leading R.
        J = R * (R_s * (phi_t * Z_z - phi_z * Z_t)
                 + phi_s * (Z_t * R_z - Z_z * R_t)
                 + Z_s * (R_t * phi_z - R_z * phi_t))

        jidx = d["jlist"][k] - 1
        G, I, iota = d["bvco"][jidx], d["buco"][jidx], d["iota"][jidx]
        psi_primes.append(np.mean(J * B**2) / (G + iota * I))

    psi_primes = np.array(psi_primes)
    spread = np.ptp(psi_primes) / np.abs(np.median(psi_primes))
    if spread > 1.0e-3:
        raise RuntimeError(
            f"psi' from geometry varies by {spread:.2e} across surfaces; "
            "boozmn file is inconsistent"
        )
    return float(np.median(psi_primes))


def main():
    parser = argparse.ArgumentParser(
        description="Convert booz_xform boozmn output to a Boozer chartmap "
                    "NetCDF file"
    )
    parser.add_argument("boozmn", help="booz_xform output file (boozmn*.nc)")
    parser.add_argument("output", help="Output Boozer chartmap .nc file")
    parser.add_argument("--nrho", type=int, default=50)
    parser.add_argument("--ntheta", type=int, default=48,
                        help="poloidal points, endpoint-excluded geometry grid")
    parser.add_argument("--nzeta", type=int, default=96,
                        help="toroidal points per field period, endpoint-excluded")
    args = parser.parse_args()

    print(f"Loading boozmn file: {args.boozmn}")
    d = read_boozmn(args.boozmn)
    nfp = d["nfp"]
    j = d["jlist"] - 1  # full-grid indices matching each half-grid surface
    s_half = d["s_half"]
    rho_half = np.sqrt(s_half)

    # Extract half-grid values of surface functions via full-grid index j.
    # j[0] = jlist[0]-1 = 1 for jlist[0]=2; skip the axis point (index 0)
    # where buco/bvco are set to zero as a boundary artefact.
    iota_h = d["iota"][j]
    buco_h = d["buco"][j]
    bvco_h = d["bvco"][j]

    # Convert booz_xform's cumulative phi_b to the chartmap A_theta coefficient.
    # The golden export test below fixes the sign against SIMPLE's VMEC path.
    if np.any(d["phi"] != 0.0):
        torflux_si = -float(d["phi"][-1]) / TWOPI
        print(f"Toroidal flux from phi_b: torflux = {torflux_si:.6e} T m^2")
    else:
        torflux_si = derive_psi_prime(d)
        print(f"Toroidal flux from geometry: torflux = {torflux_si:.6e} T m^2")

    n_rho = args.nrho
    n_theta_geom = args.ntheta
    n_phi_geom = args.nzeta

    # Uniform from rho_min: the chartmap reader assumes a strictly uniform rho
    # grid, so clipping the first point of a [0,1] grid
    # to rho_min would misplace every interior sample by up to half a cell.
    rho_grid = np.linspace(1.0e-3, 1.0, n_rho)
    s_grid = rho_grid**2
    s = np.linspace(rho_grid[0] ** 2, 1.0, n_rho)
    theta_geom = np.linspace(0.0, TWOPI, n_theta_geom, endpoint=False)
    zeta_geom = np.linspace(0.0, TWOPI / nfp, n_phi_geom, endpoint=False)

    # Surface functions on uniform rho grid.
    # buco/bvco start at jlist[0]-1 = 1 (not 0); the axis value (index 0) is
    # a boundary-condition zero and must not be included in the spline.
    iota_spline = CubicSpline(s_half, iota_h)
    B_theta = CubicSpline(s_half, buco_h)(s_grid)
    B_phi = CubicSpline(s_half, bvco_h)(s_grid)
    # A_phi(s) = -torflux * integral_0^s iota(s') ds'; sign follows SIMPLE's
    # chi = -A_phi convention (libneo spline_vmec_data, compute_boozer_data).
    iota_int = iota_spline.antiderivative()
    A_phi = -torflux_si * (iota_int(s) - iota_int(0.0))

    print(f"Evaluating Bmod on {n_rho} x {n_theta_geom} x {n_phi_geom}...")
    bmnc_g = interp_coeffs(d["bmnc"], d["ixm"], rho_half, rho_grid)
    Bmod = fourier_eval(bmnc_g, d["ixm"], d["ixn"], theta_geom, zeta_geom,
                        "cos")
    if d["lasym"]:
        bmns_g = interp_coeffs(d["bmns"], d["ixm"], rho_half, rho_grid)
        Bmod += fourier_eval(bmns_g, d["ixm"], d["ixn"], theta_geom,
                             zeta_geom, "sin")

    print(f"Evaluating geometry on {n_rho} x {n_theta_geom} x {n_phi_geom}...")
    rmnc_g = interp_coeffs(d["rmnc"], d["ixm"], rho_half, rho_grid)
    zmns_g = interp_coeffs(d["zmns"], d["ixm"], rho_half, rho_grid)
    pmns_g = interp_coeffs(d["pmns"], d["ixm"], rho_half, rho_grid)
    R = fourier_eval(rmnc_g, d["ixm"], d["ixn"], theta_geom, zeta_geom, "cos")
    Z = fourier_eval(zmns_g, d["ixm"], d["ixn"], theta_geom, zeta_geom, "sin")
    p = fourier_eval(pmns_g, d["ixm"], d["ixn"], theta_geom, zeta_geom, "sin")
    if d["lasym"]:
        rmns_g = interp_coeffs(d["rmns"], d["ixm"], rho_half, rho_grid)
        zmnc_g = interp_coeffs(d["zmnc"], d["ixm"], rho_half, rho_grid)
        pmnc_g = interp_coeffs(d["pmnc"], d["ixm"], rho_half, rho_grid)
        R += fourier_eval(rmns_g, d["ixm"], d["ixn"], theta_geom, zeta_geom,
                          "sin")
        Z += fourier_eval(zmnc_g, d["ixm"], d["ixn"], theta_geom, zeta_geom,
                          "cos")
        p += fourier_eval(pmnc_g, d["ixm"], d["ixn"], theta_geom, zeta_geom,
                          "cos")
    # Cylindrical angle: phi_cyl = zeta_B + p.
    phi_cyl = zeta_geom[None, None, :] + p
    X = R * np.cos(phi_cyl)
    Y = R * np.sin(phi_cyl)

    # SI -> CGS (chartmap stores Gauss and cm).
    Bmod = Bmod * 1.0e4              # T -> G
    B_theta = B_theta * 1.0e6        # T*m -> G*cm
    B_phi = B_phi * 1.0e6            # T*m -> G*cm
    A_phi = A_phi * 1.0e8            # T*m^2 -> G*cm^2
    torflux = torflux_si * 1.0e8     # T*m^2 -> G*cm^2
    X, Y, Z = X * 1.0e2, Y * 1.0e2, Z * 1.0e2  # m -> cm

    print(f"Writing {args.output}")
    ds = netCDF4.Dataset(args.output, "w", format="NETCDF4")
    ds.createDimension("rho", n_rho)
    ds.createDimension("s", n_rho)
    ds.createDimension("theta", n_theta_geom)
    ds.createDimension("zeta", n_phi_geom)

    v = ds.createVariable("rho", "f8", ("rho",)); v[:] = rho_grid
    v = ds.createVariable("s", "f8", ("s",)); v[:] = s
    v = ds.createVariable("theta", "f8", ("theta",)); v[:] = theta_geom
    v = ds.createVariable("zeta", "f8", ("zeta",)); v[:] = zeta_geom

    # NetCDF dimensions are (zeta, theta, rho) so NF90 reads as (rho, theta, zeta).
    for name, arr in [("x", X), ("y", Y), ("z", Z)]:
        v = ds.createVariable(name, "f8", ("zeta", "theta", "rho"))
        v[:] = np.transpose(arr, (2, 1, 0))
        v.units = "cm"

    v = ds.createVariable("A_phi", "f8", ("s",))
    v[:] = A_phi
    v.radial_abscissa = "s"
    v = ds.createVariable("B_theta", "f8", ("rho",)); v[:] = B_theta
    v = ds.createVariable("B_phi", "f8", ("rho",)); v[:] = B_phi
    v = ds.createVariable("Bmod", "f8", ("zeta", "theta", "rho"))
    v[:] = np.transpose(Bmod, (2, 1, 0))

    v = ds.createVariable("num_field_periods", "i4"); v[:] = np.int32(nfp)

    ds.rho_convention = "rho_tor"
    ds.zeta_convention = "boozer"
    ds.rho_lcfs = float(rho_grid[-1])
    ds.boozer_field = np.int32(1)
    ds.torflux = torflux
    # No rmajor attribute: the chartmap reader derives the major radius from
    # the innermost-surface geometry (see boozer_chartmap_io.f90).
    ds.booz2chartmap_source = args.boozmn

    ds.close()
    print(f"Done. nfp={nfp}, torflux={torflux:.6e} G cm^2")


if __name__ == "__main__":
    main()
