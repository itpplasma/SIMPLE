#!/usr/bin/env python3
"""
Plot orbit comparison between magfie_vmec and magfie_refcoords.

Reads orbit_refcoords_comparison.nc and creates a multi-panel comparison plot:
  1. s-theta poloidal cross-section (both trajectories overlaid)
  2. s vs time comparison
  3. theta vs time comparison
  4. phi vs time comparison
  5. mu vs time comparison (invariant conservation)
  6. Trajectory deviations vs time
"""

import sys

try:
    import xarray as xr
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError as exc:
    print(f"Skipping orbit comparison plot (missing dependency: {exc})")
    sys.exit(0)


def main():
    nc_file = "orbit_refcoords_comparison.nc"
    try:
        ds = xr.open_dataset(nc_file)
    except FileNotFoundError:
        print(f"Error: {nc_file} not found.")
        print("Run test_orbit_refcoords_rk45.x first to generate the data file.")
        sys.exit(1)

    time = ds["time"].values
    s_vmec = ds["s_vmec"].values
    theta_vmec = ds["theta_vmec"].values
    phi_vmec = ds["phi_vmec"].values
    mu_vmec = ds["mu_vmec"].values

    s_ref = ds["s_refcoords"].values
    theta_ref = ds["theta_refcoords"].values
    phi_ref = ds["phi_refcoords"].values
    mu_ref = ds["mu_refcoords"].values

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("RK45 Orbit Comparison: magfie_vmec vs magfie_refcoords",
                 fontsize=14, fontweight="bold")

    ax = axes[0, 0]
    ax.plot(theta_vmec, s_vmec, "b-", linewidth=1.5, label="vmec", alpha=0.8)
    ax.plot(theta_ref, s_ref, "r--", linewidth=1.5, label="refcoords", alpha=0.8)
    ax.scatter(theta_vmec[0], s_vmec[0], color="green", s=100, marker="o",
               label="Start", zorder=5)
    ax.scatter(theta_vmec[-1], s_vmec[-1], color="blue", s=80, marker="x",
               label="End (vmec)", zorder=5)
    ax.scatter(theta_ref[-1], s_ref[-1], color="red", s=80, marker="+",
               label="End (refcoords)", zorder=5)
    ax.set_xlabel(r"$\theta$ (rad)")
    ax.set_ylabel(r"$s$ (flux)")
    ax.set_title("Poloidal Cross-Section")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.plot(time, s_vmec, "b-", linewidth=1.5, label="vmec")
    ax.plot(time, s_ref, "r--", linewidth=1.5, label="refcoords")
    ax.set_xlabel("Time (normalized)")
    ax.set_ylabel(r"$s$ (flux)")
    ax.set_title("Flux Surface Evolution")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = axes[0, 2]
    ax.plot(time, theta_vmec, "b-", linewidth=1.5, label="vmec")
    ax.plot(time, theta_ref, "r--", linewidth=1.5, label="refcoords")
    ax.set_xlabel("Time (normalized)")
    ax.set_ylabel(r"$\theta$ (rad)")
    ax.set_title("Poloidal Angle Evolution")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = axes[1, 0]
    ax.plot(time, phi_vmec, "b-", linewidth=1.5, label="vmec")
    ax.plot(time, phi_ref, "r--", linewidth=1.5, label="refcoords")
    ax.set_xlabel("Time (normalized)")
    ax.set_ylabel(r"$\phi$ (rad)")
    ax.set_title("Toroidal Angle Evolution")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = axes[1, 1]
    ax.plot(time, mu_vmec, "b-", linewidth=1.5, label="vmec")
    ax.plot(time, mu_ref, "r--", linewidth=1.5, label="refcoords")
    ax.set_xlabel("Time (normalized)")
    ax.set_ylabel(r"$\mu$ (magnetic moment)")
    ax.set_title("Invariant Conservation")
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.ticklabel_format(style="scientific", axis="y", scilimits=(0, 0))

    ax = axes[1, 2]
    dev_s = np.abs(s_vmec - s_ref)
    dev_theta = np.abs(theta_vmec - theta_ref)
    dev_phi = np.abs(phi_vmec - phi_ref)
    ax.semilogy(time, dev_s, "g-", linewidth=1.5, label=r"$|\Delta s|$")
    ax.semilogy(time, dev_theta, "m-", linewidth=1.5, label=r"$|\Delta\theta|$")
    ax.semilogy(time, dev_phi, "c-", linewidth=1.5, label=r"$|\Delta\phi|$")
    ax.set_xlabel("Time (normalized)")
    ax.set_ylabel("Absolute Deviation")
    ax.set_title("Trajectory Deviations")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend()

    plt.tight_layout()

    png_filename = "orbit_refcoords_comparison.png"
    try:
        plt.savefig(png_filename, dpi=150, bbox_inches="tight")
    except Exception as exc:
        print(f"Skipping orbit comparison plot (savefig failed: {exc})")
        ds.close()
        sys.exit(0)
    print(f"Saved plot to {png_filename}")

    ds.close()


if __name__ == "__main__":
    main()
