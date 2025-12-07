#!/usr/bin/env python3
"""
Convert VMEC equilibrium to chartmap coordinate file.

Creates a tabulated (rho, theta, zeta) -> (X, Y, Z) mapping from VMEC data.
The output is a NetCDF file compatible with libneo chartmap_coordinate_system_t.

Two modes are available (one must be specified):

  --map2disc: Uses boundary-conforming conformal mappings via map2disc.
              Creates conformal mappings from a unit disk to the flux surface
              cross-section, ensuring good coordinate properties (orthogonality
              near the boundary, smooth behavior at the magnetic axis).
              Requires: pip install map2disc

  --passthrough: Direct VMEC Fourier evaluation without conformal mapping.
                 Simply evaluates VMEC coordinates on a regular grid.

Usage:
    python vmec_to_chartmap.py wout.nc chartmap.nc --map2disc
    python vmec_to_chartmap.py wout.nc chartmap.nc --passthrough

Requirements:
    pip install netCDF4 numpy
    pip install map2disc  # only for --map2disc mode
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

try:
    from map2disc import map as m2d
    HAS_MAP2DISC = True
except ImportError:
    HAS_MAP2DISC = False


# Conversion factor: VMEC R,Z are in meters, libneo uses cm
M_TO_CM = 100.0


def get_vmec_boundary(vmec_file: str, zeta: float = 0.0) -> tuple[np.ndarray, int]:
    """Extract VMEC boundary (R, Z) at given toroidal angle.

    Returns boundary curve function and nfp.
    """
    with Dataset(vmec_file, "r") as ds:
        xm = ds.variables["xm"][:]
        xn = ds.variables["xn"][:]
        rmnc = ds.variables["rmnc"][-1, :]  # Last flux surface = boundary
        zmns = ds.variables["zmns"][-1, :]
        nfp = int(ds.variables["nfp"][:])

        lasym = ds.variables.get("lasym__logical__")
        if lasym is not None:
            lasym = bool(lasym[:])
        else:
            lasym = False

        if lasym:
            rmns = ds.variables["rmns"][-1, :]
            zmnc = ds.variables["zmnc"][-1, :]
        else:
            rmns = None
            zmnc = None

    def boundary_curve(theta: np.ndarray) -> np.ndarray:
        """Return (R, Z) boundary points for given theta array."""
        theta = np.atleast_1d(np.asarray(theta, dtype=np.float64))
        R = np.zeros_like(theta, dtype=np.float64)
        Z = np.zeros_like(theta, dtype=np.float64)
        for im in range(len(xm)):
            m = xm[im]
            n = xn[im] / nfp
            angle = m * theta - n * zeta
            R += rmnc[im] * np.cos(angle)
            Z += zmns[im] * np.sin(angle)
            if lasym and rmns is not None and zmnc is not None:
                R += rmns[im] * np.sin(angle)
                Z += zmnc[im] * np.cos(angle)
        return np.array([R * M_TO_CM, Z * M_TO_CM])

    return boundary_curve, nfp


def vmec_to_chartmap_map2disc(
    vmec_file: str,
    output_file: str,
    nrho: int = 64,
    ntheta: int = 65,
    nzeta: int = 66,
    M: int = 16,
    Nt: int = 256,
    Ng: tuple[int, int] = (256, 256),
) -> None:
    """Convert VMEC to chartmap using map2disc conformal mapping.

    Args:
        vmec_file: Input VMEC wout file
        output_file: Output chartmap NetCDF file
        nrho: Number of radial points
        ntheta: Number of poloidal points
        nzeta: Number of toroidal points
        M: map2disc Fourier truncation parameter
        Nt: map2disc boundary discretization
        Ng: map2disc grid size for solver
    """
    if not HAS_MAP2DISC:
        raise ImportError(
            "map2disc not installed. Install with: pip install map2disc\n"
            "Or use --simple mode for direct VMEC evaluation."
        )

    # Get VMEC boundary and nfp
    _, nfp = get_vmec_boundary(vmec_file, zeta=0.0)

    rho = np.linspace(0.0, 1.0, nrho)
    theta = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
    zeta_grid = np.linspace(0.0, 2.0 * np.pi / nfp, nzeta, endpoint=False)

    print(f"Reading VMEC file: {vmec_file}")
    print(f"Grid: {nrho} x {ntheta} x {nzeta} (rho x theta x zeta)")
    print(f"Number of field periods: {nfp}")
    print(f"Using map2disc conformal mapping (M={M}, Nt={Nt})")

    X = np.zeros((nrho, ntheta, nzeta))
    Y = np.zeros((nrho, ntheta, nzeta))
    Z = np.zeros((nrho, ntheta, nzeta))

    for iz, zeta_val in enumerate(zeta_grid):
        print(f"  Processing zeta slice {iz + 1}/{nzeta} (zeta={zeta_val:.4f})")

        # Get boundary curve at this toroidal angle
        boundary_curve, _ = get_vmec_boundary(vmec_file, zeta=zeta_val)

        # Create map2disc conformal mapping for this slice
        bcm = m2d.BoundaryConformingMapping(
            curve=boundary_curve,
            M=M,
            Nt=Nt,
            Ng=Ng,
        )
        bcm.solve_domain2disk()
        bcm.solve_disk2domain()

        # Evaluate mapping: (rho, theta) -> (R, Z) in poloidal plane
        # eval_rt_1d returns shape (2, nrho, ntheta)
        rz = bcm.eval_rt_1d(rho, theta)
        R_2d = rz[0]  # shape (nrho, ntheta), already in cm
        Z_2d = rz[1]  # shape (nrho, ntheta), already in cm

        # Convert to Cartesian (X, Y, Z)
        X[:, :, iz] = R_2d * np.cos(zeta_val)
        Y[:, :, iz] = R_2d * np.sin(zeta_val)
        Z[:, :, iz] = Z_2d

    print(f"Writing chartmap file: {output_file}")
    write_chartmap(Path(output_file), rho, theta, zeta_grid, X, Y, Z, nfp)
    print("Done.")


def vmec_to_cartesian_simple(
    vmec_file: str,
    s: np.ndarray,
    theta: np.ndarray,
    zeta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Evaluate VMEC coordinates at given (s, theta, zeta) grid points.

    This is the simple/direct method without map2disc conformal mapping.
    Returns Cartesian (X, Y, Z) arrays with shape (ns, ntheta, nzeta).
    Units are in cm (same as libneo VMEC splines).
    """
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
    """Write chartmap NetCDF file in libneo format.

    NetCDF uses C (row-major) order, Fortran uses column-major order.
    Dimensions appear reversed: Python ("zeta", "theta", "rho") is read by
    Fortran as (rho, theta, zeta) with rho varying fastest.
    """
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

        # Reversed dimension order for Fortran column-major reading
        v_x = ds.createVariable("x", "f8", ("zeta", "theta", "rho"))
        v_y = ds.createVariable("y", "f8", ("zeta", "theta", "rho"))
        v_z = ds.createVariable("z", "f8", ("zeta", "theta", "rho"))
        v_nfp = ds.createVariable("nfp", "i4")

        v_rho[:] = rho
        v_theta[:] = theta
        v_zeta[:] = zeta

        # X has shape (nrho, ntheta, nzeta), transpose to (nzeta, ntheta, nrho)
        v_x[:, :, :] = np.transpose(X, (2, 1, 0))
        v_y[:, :, :] = np.transpose(Y, (2, 1, 0))
        v_z[:, :, :] = np.transpose(Z, (2, 1, 0))
        v_nfp.assignValue(nfp)


def vmec_to_chartmap_simple(
    vmec_file: str,
    output_file: str,
    nrho: int = 64,
    ntheta: int = 65,
    nzeta: int = 66,
) -> None:
    """Convert VMEC file to chartmap format using direct evaluation.

    This is a simpler method that directly evaluates VMEC coordinates without
    the conformal mapping from map2disc. Use --simple flag to enable this mode.

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
    print("Using simple/direct VMEC evaluation (no map2disc)")

    X, Y, Z = vmec_to_cartesian_simple(vmec_file, s, theta, zeta)

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

    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        "--map2disc",
        action="store_true",
        help="Use map2disc boundary-conforming conformal mapping (requires map2disc)",
    )
    mode_group.add_argument(
        "--passthrough",
        action="store_true",
        help="Use direct VMEC Fourier evaluation (no conformal mapping)",
    )

    parser.add_argument(
        "--M", type=int, default=16, help="map2disc Fourier truncation (default: 16)"
    )
    parser.add_argument(
        "--Nt", type=int, default=256, help="map2disc boundary discretization"
    )

    args = parser.parse_args(argv)

    if args.map2disc:
        if not HAS_MAP2DISC:
            print("Error: map2disc not installed. Install with: pip install map2disc")
            return 1
        vmec_to_chartmap_map2disc(
            args.vmec_file,
            args.output_file,
            nrho=args.nrho,
            ntheta=args.ntheta,
            nzeta=args.nzeta,
            M=args.M,
            Nt=args.Nt,
        )
    else:  # args.passthrough
        vmec_to_chartmap_simple(
            args.vmec_file,
            args.output_file,
            nrho=args.nrho,
            ntheta=args.ntheta,
            nzeta=args.nzeta,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
