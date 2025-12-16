#!/usr/bin/env python3
"""Plot field equivalence test results from binary data files."""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os


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
    R_2d = np.fromfile(f'{prefix}_R_2d.bin', dtype=np.float64)
    Z_2d = np.fromfile(f'{prefix}_Z_2d.bin', dtype=np.float64)

    # Reshape to 2D (Fortran column-major order)
    Bmod_direct = Bmod_direct.reshape((n_th, n_r), order='F')
    Bmod_vmec = Bmod_vmec.reshape((n_th, n_r), order='F')
    Bmod_chart = Bmod_chart.reshape((n_th, n_r), order='F')
    err_vmec = err_vmec.reshape((n_th, n_r), order='F')
    err_chart = err_chart.reshape((n_th, n_r), order='F')
    R_2d = R_2d.reshape((n_th, n_r), order='F')
    Z_2d = Z_2d.reshape((n_th, n_r), order='F')

    return (r_grid, theta_grid, Bmod_direct, Bmod_vmec, Bmod_chart,
            err_vmec, err_chart, R_2d, Z_2d)


def plot_comparison(prefix):
    """Generate comparison plots."""
    (r_grid, theta_grid, Bmod_direct, Bmod_vmec, Bmod_chart,
     err_vmec, err_chart, R_2d, Z_2d) = read_grid_data(prefix)

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
    plt.savefig('field_equiv_comparison.png', dpi=150)
    plt.close()

    # Flux surface plot
    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(R_2d.flatten(), Z_2d.flatten(),
                    c=Bmod_direct.flatten(), cmap='viridis', s=20)
    ax.set_xlabel('R [cm]')
    ax.set_ylabel('Z [cm]')
    ax.set_title('Flux surface cross-section (phi=0), colored by Bmod')
    ax.set_aspect('equal')
    plt.colorbar(sc, ax=ax, label='Bmod')
    plt.tight_layout()
    plt.savefig('field_equiv_flux_surface.png', dpi=150)
    plt.close()

    print('Plots saved: field_equiv_comparison.png, field_equiv_flux_surface.png')


if __name__ == '__main__':
    prefix = sys.argv[1] if len(sys.argv) > 1 else 'field_equiv'
    plot_comparison(prefix)
