#!/usr/bin/env python3
"""
Plot a single orbit from orbits.nc NetCDF file.

Usage: python plot_orbit.py <particle_id>

Creates four subplots:
  1. s vs theta (poloidal cross-section)
  2. s vs time
  3. theta vs time
  4. phi vs time
"""

import sys

try:
    import xarray as xr
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError as exc:
    # Allow this script to act as an optional plotting helper in
    # environments where xarray/matplotlib/numpy are not installed.
    print(f"Skipping orbit plot (missing dependency: {exc})")
    sys.exit(0)


def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_orbit.py <particle_id>")
        sys.exit(1)

    try:
        particle_id = int(sys.argv[1])
    except ValueError:
        print(f"Error: particle_id must be an integer, got '{sys.argv[1]}'")
        sys.exit(1)

    # Load NetCDF file (relative to build/test/tests when run from ctest)
    orbits_path = "../../../examples/orbit_output_test/orbits.nc"
    try:
        ds = xr.open_dataset(orbits_path)
    except FileNotFoundError:
        print(f"Error: {orbits_path} not found.")
        sys.exit(1)

    # Check if particle exists
    if particle_id not in ds.particle.values:
        print(f"Error: Particle {particle_id} not found in orbits.nc")
        print(f"Available particles: {sorted(ds.particle.values)}")
        sys.exit(1)

    # Extract orbit data
    orbit = ds.sel(particle=particle_id)
    time = orbit["time"].values
    s = orbit["s"].values
    theta = orbit["theta"].values
    phi = orbit["phi"].values

    # Remove NaN values
    valid = ~np.isnan(time)
    time = time[valid]
    s = s[valid]
    theta = theta[valid]
    phi = phi[valid]

    if len(time) == 0:
        print(f"Error: No valid data for particle {particle_id}")
        sys.exit(1)

    # Create figure with 4 subplots
    try:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f"Orbit for Particle {particle_id}", fontsize=14, fontweight="bold")

        # Subplot 1: s vs theta (poloidal plane)
        ax = axes[0, 0]
        ax.plot(theta, s, "b,", markersize=2)
        ax.scatter(theta[0], s[0], color="green", s=100, marker="o", label="Start", zorder=5)
        ax.scatter(theta[-1], s[-1], color="red", s=100, marker="x", label="End", zorder=5)
        ax.set_xlabel(r"$\theta$ (rad)")
        ax.set_ylabel(r"$s$ (flux)")
        ax.set_title("Poloidal Cross-Section (s-θ)")
        ax.grid(True, alpha=0.3)
        ax.legend()

        # Subplot 2: s vs time
        ax = axes[0, 1]
        ax.plot(time, s, "r,", markersize=2)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(r"$s$ (flux)")
        ax.set_title(r"Flux Surface Evolution")
        ax.grid(True, alpha=0.3)

        # Subplot 3: theta vs time
        ax = axes[1, 0]
        ax.plot(time, theta, "g,", markersize=2)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(r"$\theta$ (rad)")
        ax.set_title("Poloidal Angle Evolution")
        ax.grid(True, alpha=0.3)

        # Subplot 4: phi vs time
        ax = axes[1, 1]
        ax.plot(time, phi, "m,", markersize=2)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(r"$\phi$ (rad)")
        ax.set_title("Toroidal Angle Evolution")
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        # Save as PNG
        png_filename = f"orbit_{particle_id}.png"
        plt.savefig(png_filename, dpi=150, bbox_inches="tight")
        print(f"Saved plot to {png_filename}")
    except Exception as exc:
        print(f"Skipping orbit plot (matplotlib failed: {exc})")
        ds.close()
        sys.exit(0)

    ds.close()


if __name__ == "__main__":
    main()
