"""Write a closed STL wall from one radial surface of a chartmap file.

Usage: python make_wall_stl.py CHARTMAP RHO_WALL OUTPUT.stl
"""
import sys

import numpy as np
from netCDF4 import Dataset

chartmap, rho_wall, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]

with Dataset(chartmap) as ds:
    rho = np.asarray(ds["rho"][:])
    nfp = int(ds["num_field_periods"][...])
    ir = int(np.argmin(abs(rho - rho_wall)))
    # x, y, z have shape (zeta, theta, rho) and unit cm; one field period.
    xyz = np.stack([np.asarray(ds[c][:, :, ir]) for c in ("x", "y", "z")], -1)

# Rotate the field period to cover the full torus; grids exclude endpoints.
period = 2 * np.pi / nfp
parts = []
for k in range(nfp):
    c, s = np.cos(k * period), np.sin(k * period)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    parts.append(xyz @ rot.T)
points = np.concatenate(parts, axis=0) * 1e-2  # cm -> m
nz, nt, _ = points.shape

with open(out, "w") as f:
    f.write("solid wall\n")
    for i in range(nz):
        for j in range(nt):
            p00 = points[i, j]
            p10 = points[(i + 1) % nz, j]
            p01 = points[i, (j + 1) % nt]
            p11 = points[(i + 1) % nz, (j + 1) % nt]
            for a, b, c in ((p00, p10, p11), (p00, p11, p01)):
                n = np.cross(b - a, c - a)
                n /= np.linalg.norm(n)
                f.write(f"facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n")
                f.write(" outer loop\n")
                for p in (a, b, c):
                    f.write(f"  vertex {p[0]:.9e} {p[1]:.9e} {p[2]:.9e}\n")
                f.write(" endloop\nendfacet\n")
    f.write("endsolid wall\n")
print(f"rho = {rho[ir]:.4f}: {2 * nz * nt} triangles written to {out}")
