#!/usr/bin/env python3
"""
Convert VMEC equilibrium to chartmap coordinate file (passthrough mode).

Creates a tabulated (rho, theta, zeta) -> (X, Y, Z) mapping from VMEC data
by directly evaluating VMEC Fourier series on a uniform grid.

Usage:
    python vmec_to_chartmap.py wout.nc output.chartmap.nc [--nrho 64] [--ntheta 65] [--nzeta 33]

Requirements:
    pip install netCDF4 numpy
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


# Conversion factor: VMEC R,Z are in meters, libneo uses cm
M_TO_CM = 100.0


def vmec_to_chartmap_passthrough(
    vmec_file: str,
    output_file: str,
    nrho: int = 64,
    ntheta: int = 65,
    nzeta: int = 33,
) -> None:
    """Generate chartmap file by direct VMEC Fourier evaluation (passthrough mode).

    Args:
        vmec_file: Input VMEC equilibrium file (wout.nc)
        output_file: Output chartmap file (.chartmap.nc)
        nrho: Number of radial points
        ntheta: Number of poloidal points
        nzeta: Number of toroidal points
    """
    print(f"Reading VMEC file: {vmec_file}")

    # Read VMEC data
    with Dataset(vmec_file, "r") as ds:
        xm = ds.variables["xm"][:]  # Poloidal mode numbers
        xn = ds.variables["xn"][:]  # Toroidal mode numbers
        rmnc = ds.variables["rmnc"][:, :]  # R Fourier coefficients [ns, nmodes]
        zmns = ds.variables["zmns"][:, :]  # Z Fourier coefficients
        nfp = int(ds.variables["nfp"][:])
        ns_vmec = rmnc.shape[0]

        # Check for stellarator asymmetry
        lasym = ds.variables.get("lasym__logical__")
        if lasym is not None:
            lasym = bool(lasym[:])
        else:
            lasym = False

        if lasym:
            rmns = ds.variables["rmns"][:, :]
            zmnc = ds.variables["zmnc"][:, :]
        else:
            rmns = None
            zmnc = None

    print(f"VMEC data: ns={ns_vmec}, nmodes={len(xm)}, nfp={nfp}, lasym={lasym}")

    # Create uniform grid in (rho, theta, zeta)
    rho = np.linspace(0, 1, nrho)
    theta = np.linspace(0, 2*np.pi, ntheta)
    zeta = np.linspace(0, 2*np.pi/nfp, nzeta)

    print(f"Generating chartmap grid: {nrho} x {ntheta} x {nzeta}")

    # Allocate arrays
    R = np.zeros((nrho, ntheta, nzeta))
    Z = np.zeros((nrho, ntheta, nzeta))

    # Evaluate VMEC Fourier series on grid
    nmodes = len(xm)
    for js, rho_val in enumerate(rho):
        for ith, theta_val in enumerate(theta):
            for iz, zeta_val in enumerate(zeta):
                # Convert rho to VMEC s (normalized toroidal flux)
                s_val = rho_val**2

                # Linear interpolation in s (radial coordinate)
                s_idx = s_val * (ns_vmec - 1)
                js_lo = int(np.floor(s_idx))
                js_hi = min(js_lo + 1, ns_vmec - 1)
                ds_frac = s_idx - js_lo

                R_val = 0.0
                Z_val = 0.0

                # Sum Fourier series
                for im in range(nmodes):
                    m = xm[im]
                    n = xn[im] / nfp
                    angle = m * theta_val - n * zeta_val

                    # Interpolate Fourier coefficients
                    rmnc_interp = (1 - ds_frac) * rmnc[js_lo, im] + ds_frac * rmnc[js_hi, im]
                    zmns_interp = (1 - ds_frac) * zmns[js_lo, im] + ds_frac * zmns[js_hi, im]

                    R_val += rmnc_interp * np.cos(angle)
                    Z_val += zmns_interp * np.sin(angle)

                    if lasym and rmns is not None and zmnc is not None:
                        rmns_interp = (1 - ds_frac) * rmns[js_lo, im] + ds_frac * rmns[js_hi, im]
                        zmnc_interp = (1 - ds_frac) * zmnc[js_lo, im] + ds_frac * zmnc[js_hi, im]
                        R_val += rmns_interp * np.sin(angle)
                        Z_val += zmnc_interp * np.cos(angle)

                R[js, ith, iz] = R_val
                Z[js, ith, iz] = Z_val

    # Convert cylindrical (R, phi, Z) to Cartesian (X, Y, Z)
    # Note: VMEC stores in meters, libneo uses cm
    R = R * M_TO_CM
    Z = Z * M_TO_CM

    X = np.zeros((nrho, ntheta, nzeta))
    Y = np.zeros((nrho, ntheta, nzeta))
    Zcart = np.zeros((nrho, ntheta, nzeta))

    for iz, zeta_val in enumerate(zeta):
        X[:, :, iz] = R[:, :, iz] * np.cos(zeta_val)
        Y[:, :, iz] = R[:, :, iz] * np.sin(zeta_val)
        Zcart[:, :, iz] = Z[:, :, iz]

    # Write chartmap NetCDF file
    write_chartmap(output_file, rho, theta, zeta, X, Y, Zcart, nfp)

    print(f"Chartmap file written: {output_file}")
    print(f"Grid dimensions: rho={nrho}, theta={ntheta}, zeta={nzeta}")


def write_chartmap(
    outfile: str | Path,
    rho: np.ndarray,
    theta: np.ndarray,
    zeta: np.ndarray,
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
    nfp: int = 1,
) -> None:
    """Write chartmap NetCDF file in libneo format.

    NetCDF uses C (row-major) order, Fortran uses column-major order.
    Dimensions appear reversed: Python stores as (rho, theta, zeta) but
    NetCDF dimensions are declared as (zeta, theta, rho) so Fortran reads
    them in the correct order with rho varying fastest.
    """
    with Dataset(outfile, "w") as ds:
        nrho = len(rho)
        ntheta = len(theta)
        nzeta = len(zeta)

        # Create dimensions
        ds.createDimension("rho", nrho)
        ds.createDimension("theta", ntheta)
        ds.createDimension("zeta", nzeta)

        # Create coordinate variables
        v_rho = ds.createVariable("rho", "f8", ("rho",))
        v_theta = ds.createVariable("theta", "f8", ("theta",))
        v_zeta = ds.createVariable("zeta", "f8", ("zeta",))

        # Create Cartesian coordinate variables with reversed dimension order
        # for Fortran column-major reading
        v_x = ds.createVariable("x", "f8", ("zeta", "theta", "rho"))
        v_y = ds.createVariable("y", "f8", ("zeta", "theta", "rho"))
        v_z = ds.createVariable("z", "f8", ("zeta", "theta", "rho"))
        v_nfp = ds.createVariable("nfp", "i4")

        # Write coordinate arrays
        v_rho[:] = rho
        v_theta[:] = theta
        v_zeta[:] = zeta

        # Write Cartesian coordinates with transposed order
        # X has shape (nrho, ntheta, nzeta), transpose to (nzeta, ntheta, nrho)
        v_x[:, :, :] = np.transpose(X, (2, 1, 0))
        v_y[:, :, :] = np.transpose(Y, (2, 1, 0))
        v_z[:, :, :] = np.transpose(Z, (2, 1, 0))
        v_nfp.assignValue(nfp)

        # Add attributes
        ds.description = "Chartmap coordinate file generated from VMEC (passthrough mode)"
        ds.source = "vmec_to_chartmap.py"


def main():
    parser = argparse.ArgumentParser(
        description="Convert VMEC equilibrium to chartmap coordinate file (passthrough mode)"
    )
    parser.add_argument("vmec_file", help="Input VMEC equilibrium file (wout.nc)")
    parser.add_argument("output_file", help="Output chartmap file (.chartmap.nc)")
    parser.add_argument("--nrho", type=int, default=64, help="Number of radial points (default: 64)")
    parser.add_argument("--ntheta", type=int, default=65, help="Number of poloidal points (default: 65)")
    parser.add_argument("--nzeta", type=int, default=33, help="Number of toroidal points (default: 33)")

    args = parser.parse_args()

    vmec_to_chartmap_passthrough(
        args.vmec_file,
        args.output_file,
        nrho=args.nrho,
        ntheta=args.ntheta,
        nzeta=args.nzeta,
    )


if __name__ == "__main__":
    main()
