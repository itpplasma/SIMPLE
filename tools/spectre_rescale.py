#!/usr/bin/env python3
"""Rescale a SPEC or SPECTRE output file in size and field strength.

Lengths are multiplied by LENGTH and the magnetic field by FIELD, so the
equilibrium keeps its shape, rotational transform and flux fractions:

    R, Z harmonics                         -> LENGTH * R, Z
    A_theta, A_zeta (covariant, B L^2)     -> FIELD * LENGTH**2 * A
    phiedge (B L^2)                        -> FIELD * LENGTH**2 * phiedge
    mu (1/L)                               -> mu / LENGTH
    pressure (B^2)                         -> FIELD**2 * pressure

The output holds only the datasets that SIMPLE's SPECTRE reader uses, so it
carries no unscaled data.

Usage: spectre_rescale.py IN.h5 OUT.h5 --length L --field B
       spectre_rescale.py IN.h5 OUT.h5 --minor-radius A --axis-field B
"""
import argparse

import numpy as np
from netCDF4 import Dataset

SCALING = {
    "input/physics/Igeometry": (0, 0), "input/physics/Nfp": (0, 0),
    "input/physics/Mpol": (0, 0), "input/physics/Ntor": (0, 0),
    "input/physics/Nvol": (0, 0), "input/physics/Lrad": (0, 0),
    "input/physics/tflux": (0, 0), "input/physics/pflux": (0, 0),
    "input/physics/phiedge": (2, 1), "input/physics/pressure": (0, 2),
    "output/Mvol": (0, 0), "output/mn": (0, 0),
    "output/im": (0, 0), "output/in": (0, 0), "output/mu": (-1, 0),
    "output/Rbc": (1, 0), "output/Rbs": (1, 0),
    "output/Zbc": (1, 0), "output/Zbs": (1, 0),
    "vector_potential/Ate": (2, 1), "vector_potential/Ato": (2, 1),
    "vector_potential/Aze": (2, 1), "vector_potential/Azo": (2, 1),
}
"""Dataset -> (power of LENGTH, power of FIELD)."""


def outer_minor_radius(src):
    """Half the radial width of the outermost interface at zeta = 0."""
    m = np.ravel(src["output/im"][:])
    rbc = np.asarray(src["output/Rbc"][:])[-1]
    r_out, r_in = np.sum(rbc), np.sum(rbc*np.cos(m*np.pi))
    return 0.5*abs(r_out - r_in)


def mean_field(src):
    """phiedge / (pi a^2), the mean toroidal field through the boundary."""
    a = outer_minor_radius(src)
    return float(np.ravel(src["input/physics/phiedge"][:])[0])/(np.pi*a*a)


def rescale(source, target, length, field):
    with Dataset(source, "r") as src, Dataset(target, "w") as dst:
        src.set_auto_mask(False)
        for key, (p_length, p_field) in SCALING.items():
            values = np.asarray(src[key][:])
            if p_length or p_field:
                values = values*length**p_length*field**p_field
            parent, name = key.rsplit("/", 1)
            group = dst
            for part in parent.split("/"):
                group = group.groups[part] if part in group.groups \
                    else group.createGroup(part)
            dims = tuple(f"{name}_d{i}" for i in range(values.ndim))
            for dim, size in zip(dims, values.shape):
                group.createDimension(dim, size)
            var = group.createVariable(name, values.dtype, dims)
            if values.ndim:
                var[:] = values
            else:
                var.assignValue(values)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("source")
    parser.add_argument("target")
    parser.add_argument("--length", type=float, help="length factor")
    parser.add_argument("--field", type=float, help="field factor")
    parser.add_argument("--minor-radius", type=float,
                        help="target minor radius of the outer interface")
    parser.add_argument("--axis-field", type=float,
                        help="target phiedge / (pi a^2)")
    args = parser.parse_args()
    with Dataset(args.source, "r") as src:
        a0, b0 = outer_minor_radius(src), mean_field(src)
    length = args.length or (args.minor_radius/a0 if args.minor_radius else 1.0)
    field = args.field or (args.axis_field/b0 if args.axis_field else 1.0)
    rescale(args.source, args.target, length, field)
    print(f"length x {length:.6g} (a {a0:.6g} -> {a0*length:.6g}), "
          f"field x {field:.6g} (phiedge/(pi a^2) {b0:.6g} -> {b0*field:.6g})")


if __name__ == "__main__":
    main()
