#!/usr/bin/env python3
"""Convert a GVEC state to a Boozer chartmap NetCDF file.

Uses GVEC Python library to evaluate fields in Boozer coordinates
and writes the result in the extended chartmap format that SIMPLE
can read without any GVEC or VMEC library at runtime.

Usage:
    python tools/gvec_to_boozer_chartmap.py <parameter.ini> <state.dat> <output.nc>
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
    parser.add_argument("--boozer-factor", type=int, default=2)
    parser.add_argument("--Aphi-coord", choices=["rho", "rho2"], default="rho2", help="Coordinate for A_phi: rho or rho^2")
    parser.add_argument("--Bcov", choices=["avg", "boozer-avg", "boozer-0"], default="avg", help="Method for computing B_theta and B_phi surface functions.")
    parser.add_argument("--xyz-in-logical", action="store_true", help="Output XYZ in terms of the logical zeta, instead of the boozer-zeta.")
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
    print(f"Evaluating field & geometry on {n_rho} x {n_theta_field} x {n_phi_field}...")
    ev = state.evaluate_sfl(
        "mod_B", "B_theta_B", "B_zeta_B", "chi", "Phi_edge", "pos",
        rho=rho_grid,
        theta=theta_field,
        zeta=zeta_field,
        sfl="boozer",
        MNfactor=args.boozer_factor,
    )
    ev = ev.transpose("xyz", 'rad', 'pol', 'tor')
    ev_geom = ev.isel(pol=slice(0, -1), tor=slice(0, -1))  # Exclude endpoints
    ev_s = state.evaluate("chi", rho=np.sqrt(rho_grid))

    Bmod_3d = ev['mod_B'].values   # (n_rho, n_theta_field, n_phi_field)
    # Surface functions: (should be constant)
    if args.Bcov == "avg":
        ev_Bcov = state.evaluate("B_theta_avg", "B_zeta_avg", rho=rho_grid)
        B_theta = ev_Bcov["B_theta_avg"].values
        B_phi = ev_Bcov["B_zeta_avg"].values
    elif args.Bcov == "boozer-avg":
        B_theta = ev_geom['B_theta_B'].mean(["pol", "tor"]).values
        B_phi = ev_geom['B_zeta_B'].mean(["pol", "tor"]).values
    elif args.Bcov == "boozer-0":
        B_theta = ev_geom['B_theta_B'].values[:, 0, 0]
        B_phi = ev_geom['B_zeta_B'].values[:, 0, 0]
    else:
        raise ValueError(f"Invalid Bcov option: {args.Bcov}")

    if args.Aphi_coord == "rho":
        A_phi = -ev_s['chi'].values
    elif args.Aphi_coord == "rho2":
        A_phi = -ev["chi"].values
    else:
        raise ValueError(f"Invalid Aphi-coord option: {args.Aphi_coord}")

    # Toroidal flux: Phi_edge from evaluations
    # GVEC Phi_edge is already Phi/(2*pi) in SI (Wb/(2*pi) = T*m^2/(2*pi))
    torflux_SI = float(ev['Phi_edge'].item())

    # CGS conversion (GVEC: SI, SIMPLE: CGS)
    Bmod_3d = Bmod_3d * 1e4         # T -> G
    B_theta = B_theta * 1e4 * 1e2   # T*m -> G*cm
    B_phi = B_phi * 1e4 * 1e2
    torflux = torflux_SI * 1e8       # T*m^2/(2pi) -> G*cm^2
    A_phi = A_phi * 1e8              # T*m^2 -> G*cm^2

    if args.xyz_in_logical:
        ev_logical = state.evaluate("pos", rho=rho_grid, theta=theta_geom, zeta=-zeta_geom)
        ev_logical = ev_logical.transpose("xyz", 'rad', 'pol', 'tor')
        pos = ev_logical['pos'].values  # (3, n_rho, n_theta_geom, n_phi_geom)
    else:
        # Geometry on endpoint-excluded grid (Boozer angles)
        pos = ev_geom['pos'].values  # (3, n_rho, n_theta_geom, n_phi_geom)
        # pos[0]=X, pos[1]=Y, pos[2]=Z in meters

    X = pos[0] * 1e2  # m -> cm
    Y = pos[1] * 1e2
    Z_geom = pos[2] * 1e2

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

    # NetCDF dimensions are (zeta, theta, rho) so NF90 reads as (rho, theta, zeta)
    for name, arr in [("x", X), ("y", Y), ("z", Z_geom)]:
        v = ds.createVariable(name, "f8", ("zeta", "theta", "rho"))
        v[:] = np.transpose(arr, (2, 1, 0))
        v.units = "cm"

    v = ds.createVariable("A_phi", "f8", ("rho",)); v[:] = A_phi
    v = ds.createVariable("B_theta", "f8", ("rho",)); v[:] = B_theta
    v = ds.createVariable("B_phi", "f8", ("rho",)); v[:] = B_phi
    # Bmod uses field dimensions (also reversed for NF90)
    v = ds.createVariable("Bmod", "f8", ("zeta_field", "theta_field", "rho"))
    v[:] = np.transpose(Bmod_3d, (2, 1, 0))

    v = ds.createVariable("num_field_periods", "i4"); v[:] = np.int32(nfp)

    ds.rho_convention = "rho_tor"
    if args.xyz_in_logical:
        ds.zeta_convention = "cyl"
    else:
        ds.zeta_convention = "boozer"
    ds.rho_lcfs = float(rho_grid[-1])
    ds.boozer_field = np.int32(1)
    ds.torflux = torflux

    ds.close()
    print(f"Done. nfp={nfp}, torflux={torflux:.6e}")


if __name__ == "__main__":
    main()
