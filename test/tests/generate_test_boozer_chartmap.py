#!/usr/bin/env python3
"""Generate a minimal Boozer chartmap NetCDF test file with known analytic values."""

import numpy as np
import netCDF4

def generate_test_file(filename):
    """Create a test Boozer chartmap with simple analytic field.

    Uses a tokamak-like field: Bmod ~ B0 * (1 - epsilon * cos(theta))
    with B_theta ~ constant, B_phi ~ constant, A_phi ~ -torflux * iota * s.
    """
    nfp = 2
    n_rho = 17
    n_theta = 33
    n_phi = 33
    B0 = 2.0e4  # 2 T in Gauss
    R0 = 150.0  # 1.5 m in cm
    a = 50.0    # 0.5 m in cm
    epsilon = a / R0
    iota = 0.5
    torflux = B0 * np.pi * a**2 / (2.0 * np.pi)  # approximate

    twopi = 2.0 * np.pi

    rho_grid = np.linspace(0.0, 1.0, n_rho)
    rho_grid[0] = 1e-6
    theta_grid = np.linspace(0.0, twopi, n_theta, endpoint=False)
    zeta_grid = np.linspace(0.0, twopi / nfp, n_phi, endpoint=False)

    # Geometry: concentric circular cross-section tokamak
    X = np.zeros((n_rho, n_theta, n_phi))
    Y = np.zeros((n_rho, n_theta, n_phi))
    Z = np.zeros((n_rho, n_theta, n_phi))
    Bmod = np.zeros((n_rho, n_theta, n_phi))

    for ir in range(n_rho):
        rho = rho_grid[ir]
        r = a * rho  # minor radius
        for it in range(n_theta):
            theta = theta_grid[it]
            R = R0 + r * np.cos(theta)
            Zval = r * np.sin(theta)
            for ip in range(n_phi):
                phi = zeta_grid[ip]
                X[ir, it, ip] = R * np.cos(phi)
                Y[ir, it, ip] = R * np.sin(phi)
                Z[ir, it, ip] = Zval
                Bmod[ir, it, ip] = B0 / (1.0 + rho * epsilon * np.cos(theta))

    # Profiles
    A_phi = -torflux * iota * rho_grid**2
    B_theta = np.full(n_rho, B0 * epsilon * iota)
    B_phi = np.full(n_rho, B0 * R0)

    # Write
    ds = netCDF4.Dataset(filename, "w", format="NETCDF4")

    ds.createDimension("rho", n_rho)
    ds.createDimension("theta", n_theta)
    ds.createDimension("zeta", n_phi)

    v = ds.createVariable("rho", "f8", ("rho",))
    v[:] = rho_grid

    v = ds.createVariable("theta", "f8", ("theta",))
    v[:] = theta_grid

    v = ds.createVariable("zeta", "f8", ("zeta",))
    v[:] = zeta_grid

    for name, arr in [("x", X), ("y", Y), ("z", Z)]:
        v = ds.createVariable(name, "f8", ("zeta", "theta", "rho"))
        v[:] = np.transpose(arr, (2, 1, 0))
        v.units = "cm"

    v = ds.createVariable("A_phi", "f8", ("rho",))
    v[:] = A_phi

    v = ds.createVariable("B_theta", "f8", ("rho",))
    v[:] = B_theta

    v = ds.createVariable("B_phi", "f8", ("rho",))
    v[:] = B_phi

    v = ds.createVariable("Bmod", "f8", ("zeta", "theta", "rho"))
    v[:] = np.transpose(Bmod, (2, 1, 0))

    # num_field_periods must be a scalar variable (not attribute) for libneo
    v = ds.createVariable("num_field_periods", "i4")
    v[:] = np.int32(nfp)

    ds.rho_convention = "rho_tor"
    ds.zeta_convention = "cyl"
    ds.rho_lcfs = 1.0
    ds.boozer_field = np.int32(1)
    ds.torflux = torflux

    ds.close()
    print(f"Generated {filename}")
    print(f"  torflux={torflux:.6e}, B0={B0:.1f}, R0={R0:.1f}, a={a:.1f}")
    print(f"  nfp={nfp}, nrho={n_rho}, ntheta={n_theta}, nphi={n_phi}")


if __name__ == "__main__":
    generate_test_file("test_boozer_chartmap.nc")
