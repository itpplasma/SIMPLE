#!/usr/bin/env python3
"""Convert a GVEC state to a Boozer chartmap NetCDF file.

Uses GVEC Python library to evaluate fields in Boozer coordinates
and writes the result in the extended chartmap format that SIMPLE
can read without any GVEC or VMEC library at runtime.

Usage:
    .venv/bin/python gvec_to_boozer_chartmap.py <parameter.ini> <state.dat> <output.nc>
"""

import argparse
import sys
import numpy as np

try:
    import gvec
except ImportError:
    print("ERROR: gvec Python package required. Install via pip or use the .venv")
    sys.exit(1)

try:
    import netCDF4
except ImportError:
    print("ERROR: netCDF4 required. pip install netCDF4")
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Convert GVEC state to Boozer chartmap NetCDF"
    )
    parser.add_argument("paramfile", help="GVEC parameter file (.ini)")
    parser.add_argument("statefile", help="GVEC state file (.dat)")
    parser.add_argument("output", help="Output Boozer chartmap .nc file")
    parser.add_argument("--nrho", type=int, default=50)
    parser.add_argument("--ntheta", type=int, default=36)
    parser.add_argument("--nphi", type=int, default=81)
    args = parser.parse_args()

    print(f"Loading GVEC state: {args.paramfile} + {args.statefile}")
    state = gvec.load_state(args.paramfile, args.statefile)
    nfp = int(state.nfp)
    twopi = 2.0 * np.pi
    phi_period = twopi / nfp

    n_rho = args.nrho
    # Geometry: endpoint-excluded (chartmap validator)
    n_theta_geom = args.ntheta
    n_phi_geom = args.nphi
    # Field: endpoint-included (exact spline reproduction)
    n_theta_field = n_theta_geom + 1
    n_phi_field = n_phi_geom + 1

    rho_grid = np.linspace(0, 1, n_rho)
    rho_grid[0] = 1e-3
    theta_geom = np.linspace(0, twopi, n_theta_geom, endpoint=False)
    zeta_geom = np.linspace(0, phi_period, n_phi_geom, endpoint=False)
    theta_field = np.linspace(0, twopi, n_theta_field)
    zeta_field = np.linspace(0, phi_period, n_phi_field)

    # Evaluate field quantities in Boozer coordinates (endpoint-included grid)
    print(f"Evaluating field on {n_rho} x {n_theta_field} x {n_phi_field}...")
    ev = gvec.EvaluationsBoozer(
        rho=rho_grid, theta_B=theta_field, zeta_B=zeta_field, state=state)
    gvec.compute(ev, 'mod_B', state=state)
    gvec.compute(ev, 'B_theta_B', state=state)
    gvec.compute(ev, 'B_zeta_B', state=state)
    gvec.compute(ev, 'iota', state=state)

    Bmod_3d = ev['mod_B'].values   # (n_rho, n_theta_field, n_phi_field)
    B_theta_3d = ev['B_theta_B'].values
    B_phi_3d = ev['B_zeta_B'].values
    iota_vals = ev['iota'].values
    # Surface functions: take mean over angles (should be constant)
    if B_theta_3d.ndim == 3:
        B_theta = B_theta_3d[:, 0, 0]
        B_phi = B_phi_3d[:, 0, 0]
    else:
        B_theta = B_theta_3d
        B_phi = B_phi_3d
    iota_prof = iota_vals.ravel()[:n_rho] if iota_vals.ndim > 1 else iota_vals
    s_grid = rho_grid**2

    # Toroidal flux: Phi_edge from evaluations
    ev_edge = gvec.Evaluations(rho=np.array([1.0]),
                                theta=np.array([0.0]),
                                zeta=np.array([0.0]), state=state)
    gvec.compute(ev_edge, 'Phi_edge', state=state)
    # GVEC Phi_edge is already Phi/(2*pi) in SI (Wb/(2*pi) = T*m^2/(2*pi))
    torflux_SI = float(ev_edge['Phi_edge'].values.flat[0])

    # CGS conversion (GVEC: SI, SIMPLE: CGS)
    Bmod_3d = Bmod_3d * 1e4         # T -> G
    B_theta = B_theta * 1e4 * 1e2   # T*m -> G*cm
    B_phi = B_phi * 1e4 * 1e2
    torflux = torflux_SI * 1e8       # T*m^2/(2pi) -> G*cm^2

    # A_phi: integrate -torflux * iota over s
    A_phi = np.zeros(n_rho)
    for i in range(1, n_rho):
        ds = s_grid[i] - s_grid[i - 1]
        A_phi[i] = A_phi[i - 1] - torflux * iota_prof[i] * ds

    # Geometry on endpoint-excluded grid (Boozer angles)
    print(f"Evaluating geometry on {n_rho} x {n_theta_geom} x {n_phi_geom}...")
    ev_geom = gvec.EvaluationsBoozer(
        rho=rho_grid, theta_B=theta_geom, zeta_B=zeta_geom, state=state)
    gvec.compute(ev_geom, 'pos', state=state)
    pos = ev_geom['pos'].values  # (3, n_rho, n_theta_geom, n_phi_geom)
    # pos[0]=X, pos[1]=Y, pos[2]=Z in meters

    # Use pseudo-Cartesian: R from true geometry, phi_B as azimuthal angle
    X_true = pos[0] * 1e2  # m -> cm
    Y_true = pos[1] * 1e2
    Z_geom = pos[2] * 1e2
    R_geom = np.sqrt(X_true**2 + Y_true**2)

    zeta_3d = np.broadcast_to(
        zeta_geom[np.newaxis, np.newaxis, :],
        (n_rho, n_theta_geom, n_phi_geom)
    )
    X = R_geom * np.cos(zeta_3d)
    Y = R_geom * np.sin(zeta_3d)

    # Write NetCDF
    print(f"Writing {args.output}")
    ds = netCDF4.Dataset(args.output, "w", format="NETCDF4")

    ds.createDimension("rho", n_rho)
    ds.createDimension("theta", n_theta_geom)
    ds.createDimension("zeta", n_phi_geom)
    ds.createDimension("theta_field", n_theta_field)
    ds.createDimension("zeta_field", n_phi_field)

    v = ds.createVariable("rho", "f8", ("rho",)); v[:] = rho_grid
    v = ds.createVariable("theta", "f8", ("theta",)); v[:] = theta_geom
    v = ds.createVariable("zeta", "f8", ("zeta",)); v[:] = zeta_geom

    for name, arr in [("x", X), ("y", Y), ("z", Z_geom)]:
        v = ds.createVariable(name, "f8", ("rho", "theta", "zeta"))
        v[:] = arr
        v.units = "cm"

    v = ds.createVariable("A_phi", "f8", ("rho",)); v[:] = A_phi
    v = ds.createVariable("B_theta", "f8", ("rho",)); v[:] = B_theta
    v = ds.createVariable("B_phi", "f8", ("rho",)); v[:] = B_phi
    v = ds.createVariable("Bmod", "f8", ("rho", "theta_field", "zeta_field"))
    v[:] = Bmod_3d

    v = ds.createVariable("num_field_periods", "i4"); v[:] = np.int32(nfp)

    ds.rho_convention = "rho_tor"
    ds.zeta_convention = "cyl"
    ds.rho_lcfs = float(rho_grid[-1])
    ds.boozer_field = np.int32(1)
    ds.torflux = torflux

    ds.close()
    print(f"Done. nfp={nfp}, torflux={torflux:.6e}")


if __name__ == "__main__":
    main()
