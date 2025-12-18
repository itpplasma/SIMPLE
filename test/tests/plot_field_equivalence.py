#!/usr/bin/env python3
"""Plot field equivalence test results from binary data files."""

import numpy as np
import matplotlib.pyplot as plt
import sys


def read_grid_data(prefix):
    """Read grid data from binary files written by Fortran."""
    r_grid = np.fromfile(f'{prefix}_r_grid.bin', dtype=np.float64)
    theta_grid = np.fromfile(f'{prefix}_theta_grid.bin', dtype=np.float64)

    n_r = len(r_grid) - 1
    n_th = len(theta_grid) - 1

    Bmod_direct = np.fromfile(f'{prefix}_Bmod_direct.bin', dtype=np.float64)
    Bmod_vmec = np.fromfile(f'{prefix}_Bmod_vmec.bin', dtype=np.float64)
    Bmod_chart = np.fromfile(f'{prefix}_Bmod_chart.bin', dtype=np.float64)
    err_vmec = np.fromfile(f'{prefix}_err_vmec.bin', dtype=np.float64)
    err_chart = np.fromfile(f'{prefix}_err_chart.bin', dtype=np.float64)
    R_vmec = np.fromfile(f'{prefix}_R_vmec.bin', dtype=np.float64)
    Z_vmec = np.fromfile(f'{prefix}_Z_vmec.bin', dtype=np.float64)
    R_chart = np.fromfile(f'{prefix}_R_chart.bin', dtype=np.float64)
    Z_chart = np.fromfile(f'{prefix}_Z_chart.bin', dtype=np.float64)

    # Reshape to 2D (Fortran column-major order)
    Bmod_direct = Bmod_direct.reshape((n_th, n_r), order='F')
    Bmod_vmec = Bmod_vmec.reshape((n_th, n_r), order='F')
    Bmod_chart = Bmod_chart.reshape((n_th, n_r), order='F')
    err_vmec = err_vmec.reshape((n_th, n_r), order='F')
    err_chart = err_chart.reshape((n_th, n_r), order='F')
    R_vmec = R_vmec.reshape((n_th, n_r), order='F')
    Z_vmec = Z_vmec.reshape((n_th, n_r), order='F')
    R_chart = R_chart.reshape((n_th, n_r), order='F')
    Z_chart = Z_chart.reshape((n_th, n_r), order='F')

    return (r_grid, theta_grid, Bmod_direct, Bmod_vmec, Bmod_chart,
            err_vmec, err_chart, R_vmec, Z_vmec, R_chart, Z_chart)


def plot_comparison(prefix):
    """Generate comparison plots."""
    (r_grid, theta_grid, Bmod_direct, Bmod_vmec, Bmod_chart,
     err_vmec, err_chart, R_vmec, Z_vmec, R_chart, Z_chart) = read_grid_data(prefix)

    Bmin = min(Bmod_direct.min(), Bmod_vmec.min(), Bmod_chart.min())
    Bmax = max(Bmod_direct.max(), Bmod_vmec.max(), Bmod_chart.max())
    err_min = min(err_vmec.min(), err_chart.min())
    err_max = max(err_vmec.max(), err_chart.max())

    fig, axes = plt.subplots(2, 3, figsize=(12, 8))

    # Row 1: Bmod values
    im0 = axes[0, 0].pcolormesh(r_grid, theta_grid, Bmod_direct,
                                 cmap='viridis', vmin=Bmin, vmax=Bmax)
    axes[0, 0].set_xlabel('r = sqrt(s)')
    axes[0, 0].set_ylabel('theta')
    axes[0, 0].set_title('Bmod direct (coils)')
    plt.colorbar(im0, ax=axes[0, 0])

    im1 = axes[0, 1].pcolormesh(r_grid, theta_grid, Bmod_vmec,
                                 cmap='viridis', vmin=Bmin, vmax=Bmax)
    axes[0, 1].set_xlabel('r = sqrt(s)')
    axes[0, 1].set_ylabel('theta')
    axes[0, 1].set_title('Bmod VMEC-ref spline')
    plt.colorbar(im1, ax=axes[0, 1])

    im2 = axes[0, 2].pcolormesh(r_grid, theta_grid, Bmod_chart,
                                 cmap='viridis', vmin=Bmin, vmax=Bmax)
    axes[0, 2].set_xlabel('r = sqrt(s)')
    axes[0, 2].set_ylabel('theta')
    axes[0, 2].set_title('Bmod chartmap-ref spline')
    plt.colorbar(im2, ax=axes[0, 2])

    # Row 2: Errors
    im3 = axes[1, 0].pcolormesh(r_grid, theta_grid, err_vmec,
                                 cmap='RdBu_r', vmin=err_min, vmax=err_max)
    axes[1, 0].set_xlabel('r = sqrt(s)')
    axes[1, 0].set_ylabel('theta')
    axes[1, 0].set_title('VMEC error (log10)')
    plt.colorbar(im3, ax=axes[1, 0])

    im4 = axes[1, 1].pcolormesh(r_grid, theta_grid, err_chart,
                                 cmap='RdBu_r', vmin=err_min, vmax=err_max)
    axes[1, 1].set_xlabel('r = sqrt(s)')
    axes[1, 1].set_ylabel('theta')
    axes[1, 1].set_title('Chartmap error (log10)')
    plt.colorbar(im4, ax=axes[1, 1])

    im5 = axes[1, 2].pcolormesh(r_grid, theta_grid, Bmod_vmec - Bmod_chart,
                                 cmap='RdBu_r')
    axes[1, 2].set_xlabel('r = sqrt(s)')
    axes[1, 2].set_ylabel('theta')
    axes[1, 2].set_title('VMEC - Chartmap diff')
    plt.colorbar(im5, ax=axes[1, 2])

    plt.tight_layout()
    plt.savefig(f'{prefix}_comparison.png', dpi=150)
    plt.close()

    # Flux surface plot: VMEC vs Chartmap overlaid
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot VMEC flux surfaces (solid lines, colored by Bmod)
    n_th, n_r = R_vmec.shape
    for i_r in range(n_r):
        # Close the flux surface by appending the first point
        R_line = np.append(R_vmec[:, i_r], R_vmec[0, i_r])
        Z_line = np.append(Z_vmec[:, i_r], Z_vmec[0, i_r])
        ax.plot(R_line, Z_line, 'b-', linewidth=1.5, alpha=0.8,
                label='VMEC' if i_r == 0 else None)

    # Plot chartmap flux surfaces (dashed lines)
    for i_r in range(n_r):
        R_line = np.append(R_chart[:, i_r], R_chart[0, i_r])
        Z_line = np.append(Z_chart[:, i_r], Z_chart[0, i_r])
        ax.plot(R_line, Z_line, 'r--', linewidth=1.5, alpha=0.8,
                label='Chartmap' if i_r == 0 else None)

    ax.set_xlabel('R [cm]')
    ax.set_ylabel('Z [cm]')
    ax.set_title('Flux surface comparison (phi=0): VMEC (solid) vs Chartmap (dashed)')
    ax.set_aspect('equal')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{prefix}_flux_surface.png', dpi=150)
    plt.close()

    # Additional plot: Position difference between VMEC and chartmap
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    dR = R_vmec - R_chart
    dZ = Z_vmec - Z_chart
    pos_diff = np.sqrt(dR**2 + dZ**2)

    im0 = axes[0].pcolormesh(r_grid, theta_grid, pos_diff, cmap='hot')
    axes[0].set_xlabel('r = sqrt(s)')
    axes[0].set_ylabel('theta')
    axes[0].set_title('|pos_vmec - pos_chart| [cm]')
    plt.colorbar(im0, ax=axes[0])

    im1 = axes[1].pcolormesh(r_grid, theta_grid,
                              np.log10(np.maximum(pos_diff, 1e-10)), cmap='hot')
    axes[1].set_xlabel('r = sqrt(s)')
    axes[1].set_ylabel('theta')
    axes[1].set_title('log10(|pos_vmec - pos_chart|) [cm]')
    plt.colorbar(im1, ax=axes[1])

    plt.tight_layout()
    plt.savefig(f'{prefix}_position_diff.png', dpi=150)
    plt.close()

    print(
        f'Plots saved: {prefix}_comparison.png, {prefix}_flux_surface.png, '
        f'{prefix}_position_diff.png'
    )


if __name__ == '__main__':
    prefix = sys.argv[1] if len(sys.argv) > 1 else 'field_equiv'
    plot_comparison(prefix)
