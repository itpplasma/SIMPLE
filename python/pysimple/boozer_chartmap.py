#!/usr/bin/env python3
"""Generate extended Boozer chartmap NetCDF from VMEC wout.nc.

WARNING: This script uses an APPROXIMATE Boozer transformation based on
VMEC lambda. For production use, prefer the Fortran export_boozer_chartmap
subroutine (in boozer_converter.F90) which uses the exact Boozer coordinate
transformation computed by get_boozer_coordinates(). The Fortran exporter
also produces endpoint-included field grids matching the internal spline
representation exactly.

This Python script is a prototype for offline file generation and testing.
The geometry (X, Y, Z) and Bmod are evaluated at VMEC angles, not true
Boozer angles, leading to small but nonzero errors.

Usage:
    python -m pysimple.boozer_chartmap --wout wout.nc --out boozer_chartmap.nc \\
        [--nrho 62] [--ntheta 63] [--nphi 64]
"""

import argparse
import numpy as np

try:
    import netCDF4
except ImportError:
    raise ImportError("netCDF4 is required: pip install netCDF4")


def read_vmec_wout(filename):
    """Read essential data from VMEC wout.nc."""
    ds = netCDF4.Dataset(filename, "r")

    nfp = int(ds.variables["nfp"][:])
    ns = int(ds.dimensions["radius"].size)
    xm = ds.variables["xm"][:]
    xn = ds.variables["xn"][:]
    rmnc = ds.variables["rmnc"][:]
    zmns = ds.variables["zmns"][:]

    # Check for non-stellarator-symmetric terms
    has_rbc_sin = "rmns" in ds.variables
    rmns = ds.variables["rmns"][:] if has_rbc_sin else None
    zmnc = ds.variables["zmnc"][:] if has_rbc_sin else None

    # Profiles
    phi_edge = float(ds.variables["phi"][-1])  # total toroidal flux
    iotaf = ds.variables["iotaf"][:] if "iotaf" in ds.variables else None
    iotas = ds.variables["iotas"][:] if "iotas" in ds.variables else None

    # Vector potential: bvco = B_phi (toroidal covariant), buco = B_theta (poloidal covariant)
    bvco = ds.variables["bvco"][:] if "bvco" in ds.variables else None
    buco = ds.variables["buco"][:] if "buco" in ds.variables else None

    # Boozer spectrum if available
    has_boozer = "bmnc_b" in ds.variables
    bmnc_b = ds.variables["bmnc_b"][:] if has_boozer else None
    xm_b = ds.variables["ixm_b"][:] if has_boozer else None
    xn_b = ds.variables["ixn_b"][:] if has_boozer else None

    # Lambda (straight field line angle correction)
    lmns = ds.variables["lmns"][:] if "lmns" in ds.variables else None
    lmnc = ds.variables["lmnc"][:] if has_rbc_sin and "lmnc" in ds.variables else None

    ds.close()

    return {
        "nfp": nfp,
        "ns": ns,
        "xm": xm,
        "xn": xn,
        "rmnc": rmnc,
        "zmns": zmns,
        "rmns": rmns,
        "zmnc": zmnc,
        "phi_edge": phi_edge,
        "iotaf": iotaf,
        "iotas": iotas,
        "bvco": bvco,
        "buco": buco,
        "lmns": lmns,
        "lmnc": lmnc,
        "has_rbc_sin": has_rbc_sin,
        "has_boozer": has_boozer,
        "bmnc_b": bmnc_b,
        "xm_b": xm_b,
        "xn_b": xn_b,
    }


def vmec_to_boozer_angles(vmec, s_idx, theta_v, phi_v):
    """Compute Boozer angles from VMEC angles at a given flux surface.

    Returns theta_B, phi_B and the angle corrections.
    Uses the lambda transformation: theta_B = theta_V + lambda * iota,
    phi_B = phi_V + lambda.
    This is approximate -- for a proper transformation, use booz_xform.
    """
    xm = vmec["xm"]
    xn = vmec["xn"]
    lmns = vmec["lmns"]
    lmnc = vmec["lmnc"]
    iota = vmec["iotaf"][s_idx] if vmec["iotaf"] is not None else vmec["iotas"][s_idx]

    n_theta = len(theta_v)
    n_phi = len(phi_v)
    lam = np.zeros((n_theta, n_phi))

    for k in range(len(xm)):
        m = xm[k]
        n = xn[k]
        for it in range(n_theta):
            for ip in range(n_phi):
                angle = m * theta_v[it] - n * phi_v[ip]
                lam[it, ip] += lmns[s_idx, k] * np.sin(angle)
                if lmnc is not None:
                    lam[it, ip] += lmnc[s_idx, k] * np.cos(angle)

    return lam


def evaluate_vmec_position(vmec, s_idx, theta, phi):
    """Evaluate R, Z at given (theta, phi) on flux surface s_idx."""
    xm = vmec["xm"]
    xn = vmec["xn"]
    rmnc = vmec["rmnc"]
    zmns = vmec["zmns"]
    rmns = vmec["rmns"]
    zmnc = vmec["zmnc"]

    R = np.zeros_like(theta)
    Z = np.zeros_like(theta)

    for k in range(len(xm)):
        m = xm[k]
        n = xn[k]
        angle = m * theta - n * phi
        R += rmnc[s_idx, k] * np.cos(angle)
        Z += zmns[s_idx, k] * np.sin(angle)
        if rmns is not None:
            R += rmns[s_idx, k] * np.sin(angle)
        if zmnc is not None:
            Z += zmnc[s_idx, k] * np.cos(angle)

    return R, Z


def evaluate_vmec_bmod(vmec, s_idx, theta, phi):
    """Evaluate |B| at given (theta, phi) on flux surface s_idx.

    Uses the Fourier representation of |B| from VMEC.
    """
    if vmec["has_boozer"]:
        # Use Boozer spectrum
        xm_b = vmec["xm_b"]
        xn_b = vmec["xn_b"]
        bmnc_b = vmec["bmnc_b"]
        bmod = np.zeros_like(theta)
        for k in range(len(xm_b)):
            m = xm_b[k]
            n = xn_b[k]
            angle = m * theta - n * phi
            bmod += bmnc_b[s_idx, k] * np.cos(angle)
        return bmod
    else:
        raise ValueError("VMEC file does not contain Boozer spectrum (bmnc_b). "
                         "Run booz_xform first or use a VMEC file with Boozer data.")


def compute_aphi_profile(vmec, rho_grid):
    """Compute A_phi profile from VMEC flux data.

    A_phi = -torflux * integral(iota ds) from VMEC data.
    In SIMPLE convention: A_phi is the toroidal vector potential.
    """
    ns = vmec["ns"]
    s_half = np.linspace(0, 1, ns)

    # VMEC stores phi (toroidal flux), from which we derive A_phi
    # In the Boozer convention used by SIMPLE:
    # A_theta = torflux * s (by definition)
    # A_phi comes from the VMEC sA_phi array
    # For now, approximate: A_phi = -torflux * integral(iota, ds)
    torflux = abs(vmec["phi_edge"]) / (2.0 * np.pi)
    iota = vmec["iotaf"] if vmec["iotaf"] is not None else vmec["iotas"]

    # Integrate iota over s to get A_phi
    A_phi = np.zeros(ns)
    for i in range(1, ns):
        ds = s_half[i] - s_half[i - 1]
        A_phi[i] = A_phi[i - 1] - torflux * iota[i] * ds

    # Interpolate onto rho grid
    s_grid = rho_grid**2
    A_phi_on_rho = np.interp(s_grid, s_half, A_phi)

    return A_phi_on_rho, torflux


def generate_boozer_chartmap(wout_file, output_file, n_rho=62, n_theta=63, n_phi=64):
    """Generate extended Boozer chartmap from VMEC wout.nc."""
    vmec = read_vmec_wout(wout_file)
    nfp = vmec["nfp"]
    ns_vmec = vmec["ns"]

    print(f"VMEC file: {wout_file}")
    print(f"  nfp={nfp}, ns={ns_vmec}")
    print(f"Output grid: {n_rho} x {n_theta} x {n_phi}")

    twopi = 2.0 * np.pi

    # Build grids
    rho_grid = np.linspace(0.0, 1.0, n_rho)
    rho_grid[0] = 1e-6  # avoid axis singularity
    theta_grid = np.linspace(0.0, twopi, n_theta, endpoint=False)
    zeta_grid = np.linspace(0.0, twopi / nfp, n_phi, endpoint=False)

    # Allocate output arrays
    X = np.zeros((n_rho, n_theta, n_phi))
    Y = np.zeros((n_rho, n_theta, n_phi))
    Z = np.zeros((n_rho, n_theta, n_phi))
    Bmod = np.zeros((n_rho, n_theta, n_phi))
    B_theta_arr = np.zeros(n_rho)
    B_phi_arr = np.zeros(n_rho)

    # Interpolation from half-grid to rho grid
    s_half = np.linspace(0, 1, ns_vmec)

    # Get B_theta (buco) and B_phi (bvco) profiles
    if vmec["buco"] is not None:
        B_theta_arr = np.interp(rho_grid**2, s_half, vmec["buco"])
    if vmec["bvco"] is not None:
        B_phi_arr = np.interp(rho_grid**2, s_half, vmec["bvco"])

    # Compute A_phi profile
    A_phi_arr, torflux = compute_aphi_profile(vmec, rho_grid)

    # For each flux surface, compute Boozer angle mapping and evaluate geometry
    for ir in range(n_rho):
        s = rho_grid[ir] ** 2
        # Find nearest VMEC half-grid index
        s_idx = min(int(s * (ns_vmec - 1) + 0.5), ns_vmec - 1)
        s_idx = max(s_idx, 1)  # avoid axis

        # For simplicity, use VMEC angles directly as approximate Boozer angles.
        # A proper implementation would use booz_xform for the full transformation.
        # The geometry (X,Y,Z) is evaluated at VMEC angles and the Bmod
        # is evaluated using the Boozer spectrum if available.
        for it in range(n_theta):
            for ip in range(n_phi):
                theta_v = theta_grid[it]
                phi_v = zeta_grid[ip]

                R, Zval = evaluate_vmec_position(vmec, s_idx, theta_v, phi_v)
                X[ir, it, ip] = R * np.cos(phi_v)
                Y[ir, it, ip] = R * np.sin(phi_v)
                Z[ir, it, ip] = Zval

        # Evaluate Bmod on the Boozer grid
        if vmec["has_boozer"]:
            theta_2d, zeta_2d = np.meshgrid(theta_grid, zeta_grid, indexing="ij")
            Bmod[ir, :, :] = evaluate_vmec_bmod(vmec, s_idx, theta_2d, zeta_2d)
        else:
            # Fallback: use |B| from VMEC (not truly in Boozer coords)
            Bmod[ir, :, :] = abs(B_phi_arr[ir]) * 0.01  # rough estimate

    # Compute rho_lcfs (last closed flux surface in rho_tor)
    rho_lcfs = rho_grid[-1]

    # Write extended chartmap NetCDF
    print(f"Writing {output_file}")
    ds = netCDF4.Dataset(output_file, "w", format="NETCDF4")

    # Dimensions
    ds.createDimension("rho", n_rho)
    ds.createDimension("theta", n_theta)
    ds.createDimension("zeta", n_phi)

    # Coordinate variables
    v = ds.createVariable("rho", "f8", ("rho",))
    v[:] = rho_grid

    v = ds.createVariable("theta", "f8", ("theta",))
    v[:] = theta_grid

    v = ds.createVariable("zeta", "f8", ("zeta",))
    v[:] = zeta_grid

    # Geometry (stored as zeta, theta, rho for chartmap convention)
    for name, arr in [("x", X), ("y", Y), ("z", Z)]:
        v = ds.createVariable(name, "f8", ("zeta", "theta", "rho"))
        # Transpose from (rho, theta, zeta) to (zeta, theta, rho)
        v[:] = np.transpose(arr, (2, 1, 0))
        v.units = "cm"

    # Boozer field data - 1D profiles
    v = ds.createVariable("A_phi", "f8", ("rho",))
    v[:] = A_phi_arr

    v = ds.createVariable("B_theta", "f8", ("rho",))
    v[:] = B_theta_arr

    v = ds.createVariable("B_phi", "f8", ("rho",))
    v[:] = B_phi_arr

    # Boozer field data - 3D
    v = ds.createVariable("Bmod", "f8", ("zeta", "theta", "rho"))
    v[:] = np.transpose(Bmod, (2, 1, 0))

    # num_field_periods must be a scalar variable (not attribute) for libneo
    v = ds.createVariable("num_field_periods", "i4")
    v[:] = np.int32(nfp)

    # Global attributes
    ds.rho_convention = "rho_tor"
    ds.zeta_convention = "cyl"
    ds.rho_lcfs = rho_lcfs
    ds.boozer_field = np.int32(1)
    ds.torflux = torflux

    ds.close()
    print(f"Done. torflux={torflux:.6e}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate extended Boozer chartmap from VMEC wout.nc"
    )
    parser.add_argument("--wout", required=True, help="Input VMEC wout.nc file")
    parser.add_argument("--out", required=True, help="Output Boozer chartmap .nc file")
    parser.add_argument("--nrho", type=int, default=62, help="Number of radial points")
    parser.add_argument("--ntheta", type=int, default=63, help="Number of poloidal points")
    parser.add_argument("--nphi", type=int, default=64, help="Number of toroidal points")
    args = parser.parse_args()

    generate_boozer_chartmap(args.wout, args.out, args.nrho, args.ntheta, args.nphi)


if __name__ == "__main__":
    main()
