#!/usr/bin/env python3
"""Convert a QUASR simsopt JSON surface to the GVEC boundary NetCDF format."""

import argparse
import json
import math
from pathlib import Path

import netCDF4
import numpy as np


def find_surface(obj):
    for key, value in obj["simsopt_objs"].items():
        if key.startswith("SurfaceXYZTensorFourier"):
            return value
    raise ValueError("No SurfaceXYZTensorFourier object found")


def skip_mode(stellsym, dim, m, n, mpol, ntor):
    if not stellsym:
        return False
    if dim == 0:
        return (n <= ntor and m > mpol) or (n > ntor and m <= mpol)
    return (n <= ntor and m <= mpol) or (n > ntor and m > mpol)


def unpack_coefficients(surface, simsopt_objs):
    mpol = int(surface["mpol"])
    ntor = int(surface["ntor"])
    stellsym = bool(surface["stellsym"])
    dofs_ref = surface["dofs"]["value"]
    dofs = np.array(simsopt_objs[dofs_ref]["x"]["data"], dtype=float)

    coeffs = [np.zeros((2 * mpol + 1, 2 * ntor + 1)) for _ in range(3)]
    counter = 0
    for dim in range(3):
        for m in range(2 * mpol + 1):
            for n in range(2 * ntor + 1):
                if skip_mode(stellsym, dim, m, n, mpol, ntor):
                    continue
                coeffs[dim][m, n] = dofs[counter]
                counter += 1

    if counter != len(dofs):
        raise ValueError(f"Consumed {counter} dofs, expected {len(dofs)}")
    return coeffs


def basis_theta(m, theta, mpol):
    if m <= mpol:
        return math.cos(m * theta)
    return math.sin((m - mpol) * theta)


def basis_phi(n, phi, ntor, nfp):
    if n <= ntor:
        return math.cos(nfp * n * phi)
    return math.sin(nfp * (n - ntor) * phi)


def boundary_enforcer(enabled, phi, theta, nfp):
    if not enabled:
        return 1.0
    return math.sin(nfp * phi / 2.0) ** 2 + math.sin(theta / 2.0) ** 2


def evaluate_surface(surface, coeffs, ntheta, nzeta):
    mpol = int(surface["mpol"])
    ntor = int(surface["ntor"])
    nfp = int(surface["nfp"])
    clamped_dims = surface["clamped_dims"]

    theta_vals = np.linspace(0.0, 1.0, ntheta, endpoint=False)
    phi_vals = np.linspace(0.0, 1.0, nzeta * nfp, endpoint=False)
    xyz = np.zeros((nzeta * nfp, ntheta, 3))

    for iz, phi_norm in enumerate(phi_vals):
        phi = 2.0 * math.pi * phi_norm
        cosphi = math.cos(phi)
        sinphi = math.sin(phi)
        phi_basis = [basis_phi(n, phi, ntor, nfp) for n in range(2 * ntor + 1)]
        for it, theta_norm in enumerate(theta_vals):
            theta = 2.0 * math.pi * theta_norm
            theta_basis = [basis_theta(m, theta, mpol) for m in range(2 * mpol + 1)]
            xhat = 0.0
            yhat = 0.0
            z = 0.0
            for m in range(2 * mpol + 1):
                w = theta_basis[m]
                for n in range(2 * ntor + 1):
                    xhat += coeffs[0][m, n] * w * phi_basis[n] * boundary_enforcer(
                        clamped_dims[0] and n <= ntor and m <= mpol, phi, theta, nfp
                    )
                    yhat += coeffs[1][m, n] * w * phi_basis[n] * boundary_enforcer(
                        clamped_dims[1] and n <= ntor and m <= mpol, phi, theta, nfp
                    )
                    z += coeffs[2][m, n] * w * phi_basis[n] * boundary_enforcer(
                        clamped_dims[2] and n <= ntor and m <= mpol, phi, theta, nfp
                    )
            xyz[iz, it, 0] = xhat * cosphi - yhat * sinphi
            xyz[iz, it, 1] = xhat * sinphi + yhat * cosphi
            xyz[iz, it, 2] = z

    return theta_vals, phi_vals, xyz


def write_boundary(out_path, theta_vals, phi_vals, xyz, nfp):
    with netCDF4.Dataset(out_path, "w") as ds:
        ds.createDimension("zeta", xyz.shape[0])
        ds.createDimension("theta", xyz.shape[1])
        ds.createDimension("xyz", 3)
        ds.createVariable("pos", "f8", ("zeta", "theta", "xyz"))[:] = xyz
        ds.createVariable("nfp", "i4")[:] = np.int32(nfp)
        ds.createVariable("theta", "f8", ("theta",))[:] = 2.0 * np.pi * theta_vals
        ds.createVariable("zeta", "f8", ("zeta",))[:] = 2.0 * np.pi * phi_vals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_json")
    parser.add_argument("output_nc")
    parser.add_argument("--ntheta", type=int, default=81)
    parser.add_argument("--nzeta", type=int, default=81)
    args = parser.parse_args()

    obj = json.loads(Path(args.input_json).read_text())
    surface = find_surface(obj)
    coeffs = unpack_coefficients(surface, obj["simsopt_objs"])
    theta_vals, phi_vals, xyz = evaluate_surface(
        surface, coeffs, args.ntheta, args.nzeta
    )
    write_boundary(args.output_nc, theta_vals, phi_vals, xyz, int(surface["nfp"]))

    print(args.output_nc)
    print("shape", xyz.shape)


if __name__ == "__main__":
    main()
