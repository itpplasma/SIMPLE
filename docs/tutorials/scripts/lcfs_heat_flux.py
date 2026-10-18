"""Heat flux of lost alphas on the last closed flux surface.

Usage: python lcfs_heat_flux.py RESULTS.nc CHARTMAP.nc

Bins the guiding-centre exit points of lost markers in the poloidal and
toroidal angle and divides the carried power by the area of each bin.
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from pysimple.engineering import compute_bin_areas_from_chartmap

P_ALPHA = 600.0  # MW, alpha power of a 3 GW D-T reactor

results, chartmap = sys.argv[1], sys.argv[2]
with Dataset(results) as ds:
    exit_code = np.asarray(ds["orbit_exit_code"][:])
    zend = np.asarray(ds["zend"][:]).T  # rows: s, theta, phi, v/v0, v_par/v
    n_total = int(ds.ntestpart)

lost = (exit_code == 1) | (exit_code == 2)  # LCFS or wall
theta = (zend[1, lost] + np.pi) % (2 * np.pi) - np.pi
phi = zend[2, lost] % (2 * np.pi)
weight = zend[3, lost] ** 2  # fraction of the birth energy at the exit

theta_edges = np.linspace(-np.pi, np.pi, 9)
phi_edges = np.linspace(0, 2 * np.pi, 17)
energy, _, _ = np.histogram2d(theta, phi, [theta_edges, phi_edges], weights=weight)
area = compute_bin_areas_from_chartmap(chartmap, theta_edges, phi_edges)
flux = energy * (P_ALPHA / n_total) / area  # MW/m^2

print(f"lost markers: {lost.sum()} of {n_total} ({lost.mean():.1%})")
print(f"lost power: {weight.sum() / n_total * P_ALPHA:.0f} MW")
print(f"LCFS area: {area.sum():.0f} m^2")
print(f"peak flux: {flux.max():.2f} MW/m^2, mean flux: {flux.mean():.3f} MW/m^2")

plt.pcolormesh(phi_edges, theta_edges, flux, cmap="inferno")
plt.colorbar(label="heat flux (MW/m$^2$)")
plt.xlabel("toroidal angle $\\varphi$ (rad)")
plt.ylabel("poloidal angle $\\vartheta$ (rad)")
plt.savefig("lcfs_heat_flux.png", dpi=150)
