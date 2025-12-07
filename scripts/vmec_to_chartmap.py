#!/usr/bin/env python3
"""
Convert VMEC equilibrium to chartmap coordinate file.

Creates a tabulated (rho, theta, zeta) -> (X, Y, Z) mapping from VMEC data
that can be used as reference coordinates in SIMPLE. The output is a NetCDF
file compatible with libneo chartmap_coordinate_system_t.

Usage:
    python vmec_to_chartmap.py wout.nc chartmap.nc --nrho 64 --ntheta 65 --nzeta 66
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def vmec_to_cartesian(
    vmec_file: str,
    s: np.ndarray,
    theta: np.ndarray,
    zeta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Evaluate VMEC coordinates at given (s, theta, zeta) grid points.

    Returns Cartesian (X, Y, Z) arrays with shape (ns, ntheta, nzeta).
    Units are in cm (same as libneo VMEC splines).
    """
    # Conversion factor: VMEC R,Z are in meters, libneo uses cm
    M_TO_CM = 100.0

    with Dataset(vmec_file, "r") as ds:
        xm = ds.variables["xm"][:]
        xn = ds.variables["xn"][:]
        rmnc = ds.variables["rmnc"][:, :]
        zmns = ds.variables["zmns"][:, :]
        nfp = int(ds.variables["nfp"][:])

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

    ns_vmec = rmnc.shape[0]
    nmodes = len(xm)

    ns = len(s)
    nth = len(theta)
    nz = len(zeta)

    R = np.zeros((ns, nth, nz))
    Z = np.zeros((ns, nth, nz))

    for iz, zeta_val in enumerate(zeta):
        for ith, theta_val in enumerate(theta):
            for js, s_val in enumerate(s):
                s_idx = s_val * (ns_vmec - 1)
                js_lo = int(np.floor(s_idx))
                js_hi = min(js_lo + 1, ns_vmec - 1)
                ds_frac = s_idx - js_lo

                R_val = 0.0
                Z_val = 0.0

                for im in range(nmodes):
                    m = xm[im]
                    n = xn[im] / nfp
                    angle = m * theta_val - n * zeta_val

                    rmnc_interp = (1 - ds_frac) * rmnc[js_lo, im] + ds_frac * rmnc[
                        js_hi, im
                    ]
                    zmns_interp = (1 - ds_frac) * zmns[js_lo, im] + ds_frac * zmns[
                        js_hi, im
                    ]

                    R_val += rmnc_interp * np.cos(angle)
                    Z_val += zmns_interp * np.sin(angle)

                    if lasym and rmns is not None and zmnc is not None:
                        rmns_interp = (1 - ds_frac) * rmns[js_lo, im] + ds_frac * rmns[
                            js_hi, im
                        ]
                        zmnc_interp = (1 - ds_frac) * zmnc[js_lo, im] + ds_frac * zmnc[
                            js_hi, im
                        ]
                        R_val += rmns_interp * np.sin(angle)
                        Z_val += zmnc_interp * np.cos(angle)

                R[js, ith, iz] = R_val
                Z[js, ith, iz] = Z_val

    # Convert from meters to cm (VMEC stores in meters, libneo uses cm)
    R = R * M_TO_CM
    Z = Z * M_TO_CM

    X = np.zeros((ns, nth, nz))
    Y = np.zeros((ns, nth, nz))
    Zcart = np.zeros((ns, nth, nz))

    for iz, zeta_val in enumerate(zeta):
        X[:, :, iz] = R[:, :, iz] * np.cos(zeta_val)
        Y[:, :, iz] = R[:, :, iz] * np.sin(zeta_val)
        Zcart[:, :, iz] = Z[:, :, iz]

    return X, Y, Zcart


def write_chartmap(
    outfile: Path,
    rho: np.ndarray,
    theta: np.ndarray,
    zeta: np.ndarray,
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
    nfp: int = 1,
) -> None:
    """Write chartmap NetCDF file in libneo format."""
    with Dataset(outfile, "w") as ds:
        nrho = len(rho)
        ntheta = len(theta)
        nzeta = len(zeta)

        ds.createDimension("rho", nrho)
        ds.createDimension("theta", ntheta)
        ds.createDimension("zeta", nzeta)

        v_rho = ds.createVariable("rho", "f8", ("rho",))
        v_theta = ds.createVariable("theta", "f8", ("theta",))
        v_zeta = ds.createVariable("zeta", "f8", ("zeta",))

        v_x = ds.createVariable("x", "f8", ("zeta", "theta", "rho"))
        v_y = ds.createVariable("y", "f8", ("zeta", "theta", "rho"))
        v_z = ds.createVariable("z", "f8", ("zeta", "theta", "rho"))
        v_nfp = ds.createVariable("nfp", "i4")

        v_rho[:] = rho
        v_theta[:] = theta
        v_zeta[:] = zeta

        v_x[:, :, :] = np.transpose(X, (2, 1, 0))
        v_y[:, :, :] = np.transpose(Y, (2, 1, 0))
        v_z[:, :, :] = np.transpose(Z, (2, 1, 0))
        v_nfp.assignValue(nfp)


def vmec_to_chartmap(
    vmec_file: str,
    output_file: str,
    nrho: int = 64,
    ntheta: int = 65,
    nzeta: int = 66,
) -> None:
    """Convert VMEC file to chartmap format.

    Note: Although the chartmap format uses 'rho' as the variable name,
    we store the data at VMEC s values (s = psi/psi_edge, the normalized
    toroidal flux) for direct compatibility with VMEC coordinate system.
    This means rho = s, not rho = sqrt(s).
    """
    with Dataset(vmec_file, "r") as ds:
        nfp = int(ds.variables["nfp"][:])

    # Use s directly (not sqrt(s)) for compatibility with VMEC coordinate system
    # The variable is called 'rho' in chartmap format but represents s here
    rho = np.linspace(0.0, 1.0, nrho)
    theta = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi / nfp, nzeta, endpoint=False)

    s = rho  # Use s directly, not s = rho**2

    print(f"Reading VMEC file: {vmec_file}")
    print(f"Grid: {nrho} x {ntheta} x {nzeta} (rho x theta x zeta)")
    print(f"Number of field periods: {nfp}")

    X, Y, Z = vmec_to_cartesian(vmec_file, s, theta, zeta)

    print(f"Writing chartmap file: {output_file}")
    write_chartmap(Path(output_file), rho, theta, zeta, X, Y, Z, nfp)

    print("Done.")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="vmec_to_chartmap",
        description="Convert VMEC equilibrium to chartmap coordinate file",
    )
    parser.add_argument("vmec_file", help="Input VMEC wout file (.nc)")
    parser.add_argument("output_file", help="Output chartmap file (.nc or .h5)")
    parser.add_argument("--nrho", type=int, default=64, help="Number of rho points")
    parser.add_argument(
        "--ntheta", type=int, default=65, help="Number of theta points"
    )
    parser.add_argument("--nzeta", type=int, default=66, help="Number of zeta points")

    args = parser.parse_args(argv)

    vmec_to_chartmap(
        args.vmec_file,
        args.output_file,
        nrho=args.nrho,
        ntheta=args.ntheta,
        nzeta=args.nzeta,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
