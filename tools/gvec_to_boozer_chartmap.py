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
    parser.add_argument("--boozer-factor", type=int, default=1)
    parser.add_argument("--Bcov", choices=["avg", "boozer-avg", "boozer-0"], default="boozer-avg", help="Method for computing B_theta and B_phi surface functions.")
    parser.add_argument("--flip", choices=["pol", "tor"], default="tor", help="Flip the sign of the poloidal or toroidal angle (to obtain left-handed coordinates).")
    args = parser.parse_args()

    print(f"Loading GVEC state: {args.paramfile} + {args.statefile}")
    state = gvec.load_state(args.paramfile, args.statefile)
    nfp = int(state.nfp)
    twopi = 2.0 * np.pi
    phi_period = twopi / nfp

    n_rho = args.nrho
    n_theta_geom = args.ntheta
    n_phi_geom = args.nphi

    rho_grid = np.linspace(1e-3, 1, n_rho)
    s = np.linspace(rho_grid[0] ** 2, 1.0, n_rho)
    rho_profile = np.sqrt(s)
    theta_geom = np.linspace(0, twopi, n_theta_geom, endpoint=False)
    zeta_geom = np.linspace(0, phi_period, n_phi_geom, endpoint=False)

    print(f"Evaluating field & geometry on {n_rho} x {n_theta_geom} x {n_phi_geom}...")
    ev = state.evaluate_sfl(
        "mod_B", "B_theta_B", "B_zeta_B", "Phi_edge", "pos",
        rho=rho_grid,
        theta=-theta_geom if args.flip == "pol" else theta_geom,
        zeta=-zeta_geom if args.flip == "tor" else zeta_geom,
        sfl="boozer",
        MNfactor=args.boozer_factor,
    )
    ev = ev.transpose("xyz", 'rad', 'pol', 'tor')

    Bmod_3d = ev.mod_B.values   # (n_rho, n_theta_geom, n_phi_geom)

    # Surface functions: (should be constant)
    if args.Bcov == "avg":
        ev_Bcov = state.evaluate("B_theta_avg", "B_zeta_avg", rho=rho_grid)
        B_theta = ev_Bcov.B_theta_avg.values
        B_phi = ev_Bcov.B_zeta_avg.values
    elif args.Bcov == "boozer-avg":
        B_theta = ev.B_theta_B.mean(["pol", "tor"]).values
        B_phi = ev.B_zeta_B.mean(["pol", "tor"]).values
    elif args.Bcov == "boozer-0":
        B_theta = ev.B_theta_B.values[:, 0, 0]
        B_phi = ev.B_zeta_B.values[:, 0, 0]
    else:
        raise ValueError(f"Invalid Bcov option: {args.Bcov}")
    ev_aphi = state.evaluate_sfl(
        "chi",
        rho=rho_profile,
        theta=np.array([0.0]),
        zeta=np.array([0.0]),
        sfl="boozer",
        MNfactor=args.boozer_factor,
    )
    A_phi = -np.asarray(ev_aphi.chi.values).reshape(n_rho, -1)[:, 0]

    # A_theta (on edge) = toroidal flux
    # GVEC Phi_edge is already Phi/(2*pi) in SI (Wb/(2*pi) = T*m^2/(2*pi))
    # the 2*pi factor is used to convert between vector potential components and integral fluxes
    # GVEC's profiles correspond to the vector potential components
    A_theta_edge = ev.Phi_edge.item()

    pos = ev.pos.values  # (3, n_rho, n_theta_geom, n_phi_geom)
    X, Y, Z = pos[0], pos[1], pos[2]

    # GVEC right-handed coordinates -> SIMPLE left-handed coordinates
    if args.flip == "tor":
        B_phi = -B_phi
        A_phi = -A_phi
    elif args.flip == "pol":
        B_theta = -B_theta
        A_theta_edge = -A_theta_edge
    else:
        raise ValueError(f"Invalid flip option: {args.flip}")

    # CGS conversion (GVEC: SI, SIMPLE: CGS)
    Bmod_3d = Bmod_3d * 1e4            # T -> G
    B_theta = B_theta * 1e4 * 1e2      # T*m -> G*cm
    B_phi = B_phi * 1e4 * 1e2          # T*m -> G*cm 
    A_theta_edge = A_theta_edge * 1e8  # T*m^2 -> G*cm^2
    A_phi = A_phi * 1e8                # T*m^2 -> G*cm^2
    X = X * 1e2  # m -> cm
    Y = Y * 1e2
    Z = Z * 1e2

    # Write NetCDF
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

    # NetCDF dimensions are (zeta, theta, rho) so NF90 reads as (rho, theta, zeta)
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
    v[:] = np.transpose(Bmod_3d, (2, 1, 0))

    v = ds.createVariable("num_field_periods", "i4"); v[:] = np.int32(nfp)

    ds.rho_convention = "rho_tor"
    ds.zeta_convention = "boozer"
    ds.rho_lcfs = float(rho_grid[-1])
    ds.boozer_field = np.int32(1)
    ds.torflux = A_theta_edge
    ds.gvec2chartmap_boozer_factor = args.boozer_factor
    ds.gvec2chartmap_Bcov_method = args.Bcov
    ds.gvec2chartmap_flip = args.flip

    ds.close()
    print(f"Done. nfp={nfp}, torflux={A_theta_edge:.6e}")


if __name__ == "__main__":
    main()
