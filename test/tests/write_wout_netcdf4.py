"""Copy a wout file to netCDF-4 (HDF5) with netCDF4-python's default layout."""
import sys

import netCDF4


def main(source, target):
    with netCDF4.Dataset(source) as src, netCDF4.Dataset(target, "w", format="NETCDF4") as dst:
        for name, dim in src.dimensions.items():
            dst.createDimension(name, len(dim))
        for name, var in src.variables.items():
            out = dst.createVariable(name, var.datatype, var.dimensions)
            out[...] = var[...]


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
