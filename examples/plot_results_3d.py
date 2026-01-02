#!/usr/bin/env python3
"""
Visualize SIMPLE particle tracing results from NetCDF output.

Creates side-by-side 3D scatter plots showing start and end positions
of particles, with lost particles highlighted.

Usage:
    python plot_results_3d.py results.nc
    python plot_results_3d.py results.nc --lost-only
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def load_results(results_path: str | Path) -> dict:
    """Load results from NetCDF file."""
    try:
        import netCDF4 as nc
    except ImportError:
        raise ImportError("netCDF4 required: pip install netCDF4")

    results_path = Path(results_path)
    if not results_path.exists():
        raise FileNotFoundError(f"Results file not found: {results_path}")

    with nc.Dataset(results_path, 'r') as ds:
        results = {
            'xstart_cart': ds.variables['xstart_cart'][:],
            'xend_cart': ds.variables['xend_cart'][:],
            'class_lost': ds.variables['class_lost'][:].astype(bool),
            'times_lost': ds.variables['times_lost'][:],
            'ntestpart': ds.ntestpart,
            'trace_time': ds.trace_time,
        }
    return results


def plot_3d_positions(results: dict, lost_only: bool = False,
                      output_path: str | None = None) -> None:
    """Create side-by-side 3D scatter plots of start/end positions."""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    # NetCDF stores (xyz, particle) in Fortran order, which is
    # (particle, xyz) in Python (C order)
    xstart = results['xstart_cart']  # (n_particles, 3)
    xend = results['xend_cart']
    lost = results['class_lost']

    if lost_only:
        xstart = xstart[lost, :]
        xend = xend[lost, :]
        title_suffix = " (lost particles only)"
        n_shown = np.sum(lost)
    else:
        title_suffix = ""
        n_shown = len(lost)

    if n_shown == 0:
        print("No particles to plot.")
        return

    fig = plt.figure(figsize=(14, 6))

    # Start positions
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.scatter(xstart[:, 0], xstart[:, 1], xstart[:, 2],
                c='blue', alpha=0.6, s=20, label='Start')
    ax1.set_xlabel('X (cm)')
    ax1.set_ylabel('Y (cm)')
    ax1.set_zlabel('Z (cm)')
    ax1.set_title(f'Start Positions{title_suffix}')

    # End positions
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.scatter(xend[:, 0], xend[:, 1], xend[:, 2],
                c='red', alpha=0.6, s=20, label='End')
    ax2.set_xlabel('X (cm)')
    ax2.set_ylabel('Y (cm)')
    ax2.set_zlabel('Z (cm)')
    ax2.set_title(f'End Positions{title_suffix}')

    # Match axis limits
    all_coords = np.vstack([xstart, xend])
    max_range = np.max(np.abs(all_coords)) * 1.1
    for ax in [ax1, ax2]:
        ax.set_xlim(-max_range, max_range)
        ax.set_ylim(-max_range, max_range)
        ax.set_zlim(-max_range, max_range)

    n_lost = np.sum(results['class_lost'])
    n_total = len(results['class_lost'])
    fig.suptitle(f'SIMPLE Results: {n_lost}/{n_total} particles lost '
                 f'({100*n_lost/n_total:.1f}%)',
                 fontsize=12)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved to {output_path}")
    else:
        plt.show()


def plot_wall_heatmap(results: dict, output_path: str | None = None) -> None:
    """Create 2D heatmap of particle end positions on wall (R-Z plane)."""
    import matplotlib.pyplot as plt

    xend = results['xend_cart']
    lost = results['class_lost']

    # Only plot lost particles (those that hit the wall)
    xend_lost = xend[lost, :]
    if len(xend_lost) == 0:
        print("No lost particles to plot.")
        return

    # Convert to R-Z cylindrical coordinates
    x, y, z = xend_lost[:, 0], xend_lost[:, 1], xend_lost[:, 2]
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # R-Z scatter plot
    ax1 = axes[0]
    ax1.scatter(r, z, alpha=0.5, s=10, c='red')
    ax1.set_xlabel('R (cm)')
    ax1.set_ylabel('Z (cm)')
    ax1.set_title('Wall Hit Locations (R-Z)')
    ax1.set_aspect('equal')

    # Phi-Z scatter plot (toroidal angle vs Z)
    ax2 = axes[1]
    ax2.scatter(np.degrees(phi), z, alpha=0.5, s=10, c='red')
    ax2.set_xlabel('Phi (degrees)')
    ax2.set_ylabel('Z (cm)')
    ax2.set_title('Wall Hit Locations (Phi-Z)')

    n_lost = np.sum(lost)
    n_total = len(lost)
    fig.suptitle(f'Wall Heat Load: {n_lost} particles lost '
                 f'({100*n_lost/n_total:.1f}%)',
                 fontsize=12)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved to {output_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Visualize SIMPLE particle tracing results')
    parser.add_argument('results_file', type=str,
                        help='Path to results.nc file')
    parser.add_argument('--lost-only', action='store_true',
                        help='Show only lost particles')
    parser.add_argument('--heatmap', action='store_true',
                        help='Show wall heatmap instead of 3D scatter')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output file path (default: display)')

    args = parser.parse_args()

    results = load_results(args.results_file)

    print(f"Loaded {results['ntestpart']} particles")
    print(f"  Lost: {np.sum(results['class_lost'])}")
    print(f"  Confined: {np.sum(~results['class_lost'])}")
    print(f"  Trace time: {results['trace_time']:.2e} s")

    if args.heatmap:
        plot_wall_heatmap(results, output_path=args.output)
    else:
        plot_3d_positions(results, lost_only=args.lost_only,
                          output_path=args.output)


if __name__ == '__main__':
    main()
