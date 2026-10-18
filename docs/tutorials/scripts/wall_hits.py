"""Summarise STL wall hits in a SIMPLE results.nc file.

Usage: python wall_hits.py RESULTS.nc
"""
import sys

import numpy as np
from netCDF4 import Dataset

with Dataset(sys.argv[1]) as ds:
    exit_code = np.asarray(ds["orbit_exit_code"][:])
    hit = np.asarray(ds["wall_hit"][:]).astype(bool)
    x_hit = np.asarray(ds["wall_hit_cart"][:])  # (particle, xyz)
    to_m = 0.01 if ds["wall_hit_cart"].units == "cm" else 1.0
    angle = np.asarray(ds["wall_hit_angle_rad"][:])  # to the wall normal
    times = np.asarray(ds["times_lost"][:])

codes, counts = np.unique(exit_code, return_counts=True)
print("exit codes:", dict(zip(codes.tolist(), counts.tolist())))
print(f"wall hits: {hit.sum()} of {hit.size}")
if hit.any():
    r = np.hypot(x_hit[hit, 0], x_hit[hit, 1]) * to_m
    z = x_hit[hit, 2] * to_m
    print(f"hit points: R = {r.min():.2f}..{r.max():.2f} m, "
          f"Z = {z.min():.2f}..{z.max():.2f} m")
    print(f"median incidence angle: {np.degrees(np.median(angle[hit])):.0f} deg")
    print(f"median loss time: {np.median(times[hit]):.2e} s")
