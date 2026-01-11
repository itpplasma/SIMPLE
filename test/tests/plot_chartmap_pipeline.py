#!/usr/bin/env python3
"""Plot chartmap pipeline unit test results."""

import sys

try:
    import numpy as np
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
except Exception as exc:
    print(f"Skipping chartmap pipeline plots (matplotlib/numpy unavailable: {exc})")
    raise SystemExit(0)


def _save_or_skip(path: str) -> bool:
    try:
        plt.savefig(path, dpi=150)
        return True
    except Exception as exc:
        print(f"Skipping chartmap pipeline plots (savefig failed: {exc})")
        return False


def read_test_data(prefix):
    """Read test data from binary files."""
    r_grid = np.fromfile(f'{prefix}_r_grid.bin', dtype=np.float64)
    theta_grid = np.fromfile(f'{prefix}_theta_grid.bin', dtype=np.float64)
    err_vmec = np.fromfile(f'{prefix}_err_vmec.bin', dtype=np.float64)
    err_chart = np.fromfile(f'{prefix}_err_chart.bin', dtype=np.float64)

    n_r = len(r_grid) - 1
    n_th = len(theta_grid) - 1

    err_vmec = err_vmec.reshape((n_th, n_r), order='F')
    err_chart = err_chart.reshape((n_th, n_r), order='F')

    return r_grid, theta_grid, err_vmec, err_chart


def plot_test_A():
    """Plot Test A: roundtrip error."""
    r_grid, theta_grid, err_vmec, err_chart = read_test_data('pipeline_A')

    try:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        vmin = min(
            err_vmec[err_vmec > -15].min() if np.any(err_vmec > -15) else -15,
            err_chart[err_chart > -15].min() if np.any(err_chart > -15) else -15,
        )
        vmax = max(err_vmec.max(), err_chart.max())

        im0 = axes[0].pcolormesh(
            r_grid, theta_grid, err_vmec, cmap='hot', vmin=vmin, vmax=vmax
        )
        axes[0].set_xlabel('r = sqrt(s)')
        axes[0].set_ylabel('theta')
        axes[0].set_title('VMEC->chart: log10(cyl roundtrip error)')
        plt.colorbar(im0, ax=axes[0])

        im1 = axes[1].pcolormesh(
            r_grid, theta_grid, err_chart, cmap='hot', vmin=vmin, vmax=vmax
        )
        axes[1].set_xlabel('r = sqrt(s)')
        axes[1].set_ylabel('theta')
        axes[1].set_title('Chartmap: log10(u roundtrip error)')
        plt.colorbar(im1, ax=axes[1])

        plt.suptitle('Test A: Coordinate mapping roundtrip')
        plt.tight_layout()
        if not _save_or_skip('pipeline_A_roundtrip.png'):
            return
        plt.close()
        print('Saved: pipeline_A_roundtrip.png')
    except Exception as exc:
        print(f"Skipping chartmap pipeline plots (matplotlib failed: {exc})")
        return


def plot_test_B():
    """Plot Test B: covariant basis FD error."""
    r_grid, theta_grid, err_vmec, err_chart = read_test_data('pipeline_B')

    try:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        vmin = min(err_vmec.min(), err_chart.min())
        vmax = max(err_vmec.max(), err_chart.max())

        im0 = axes[0].pcolormesh(
            r_grid, theta_grid, err_vmec, cmap='hot', vmin=vmin, vmax=vmax
        )
        axes[0].set_xlabel('r = sqrt(s)')
        axes[0].set_ylabel('theta')
        axes[0].set_title('VMEC: log10(basis FD error)')
        plt.colorbar(im0, ax=axes[0])

        im1 = axes[1].pcolormesh(
            r_grid, theta_grid, err_chart, cmap='hot', vmin=vmin, vmax=vmax
        )
        axes[1].set_xlabel('r = sqrt(s)')
        axes[1].set_ylabel('theta')
        axes[1].set_title('Chartmap: log10(basis FD error)')
        plt.colorbar(im1, ax=axes[1])

        plt.suptitle('Test B: covariant_basis vs finite-difference')
        plt.tight_layout()
        if not _save_or_skip('pipeline_B_basis_fd.png'):
            return
        plt.close()
        print('Saved: pipeline_B_basis_fd.png')
    except Exception as exc:
        print(f"Skipping chartmap pipeline plots (matplotlib failed: {exc})")
        return


def plot_test_C():
    """Plot Test C: tensor transform identity error."""
    r_grid, theta_grid, err_vmec, err_chart = read_test_data('pipeline_C')

    try:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        vmin = min(err_vmec.min(), err_chart.min())
        vmax = max(err_vmec.max(), err_chart.max())

        im0 = axes[0].pcolormesh(
            r_grid, theta_grid, err_vmec, cmap='hot', vmin=vmin, vmax=vmax
        )
        axes[0].set_xlabel('r = sqrt(s)')
        axes[0].set_ylabel('theta')
        axes[0].set_title('VMEC: log10(transform identity error)')
        plt.colorbar(im0, ax=axes[0])

        im1 = axes[1].pcolormesh(
            r_grid, theta_grid, err_chart, cmap='hot', vmin=vmin, vmax=vmax
        )
        axes[1].set_xlabel('r = sqrt(s)')
        axes[1].set_ylabel('theta')
        axes[1].set_title('Chartmap: log10(transform identity error)')
        plt.colorbar(im1, ax=axes[1])

        plt.suptitle('Test C: Cartesian <-> ref basis transform identity')
        plt.tight_layout()
        if not _save_or_skip('pipeline_C_transform.png'):
            return
        plt.close()
        print('Saved: pipeline_C_transform.png')
    except Exception as exc:
        print(f"Skipping chartmap pipeline plots (matplotlib failed: {exc})")
        return


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: plot_chartmap_pipeline.py [A|B|C]')
        sys.exit(1)

    test = sys.argv[1].upper()
    if test == 'A':
        plot_test_A()
    elif test == 'B':
        plot_test_B()
    elif test == 'C':
        plot_test_C()
    else:
        print(f'Unknown test: {test}')
        sys.exit(1)
