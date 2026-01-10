#!/usr/bin/env python3
"""
Plot chartmap wall losses comparison.

Reads chartmap_wall_losses.nc and creates:
1. Loss fraction over time comparison (VMEC vs Chartmap vs Meiss)
2. Loss positions on wall in (theta, phi) space
3. Summary bar chart of confined fractions
"""

import numpy as np
import netCDF4 as nc
import sys
from pathlib import Path


def _maybe_import_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        return plt
    except Exception as exc:
        print(f"Skipping plot_chartmap_wall_losses (matplotlib unavailable: {exc})")
        return None


def load_data(filename):
    """Load loss data from NetCDF file."""
    ds = nc.Dataset(filename, 'r')

    data = {
        'trace_time': ds.getncattr('trace_time'),
        'n_steps': ds.getncattr('n_steps'),
        'z0': ds.variables['z0'][:],
        'vmec': {
            'loss_time': ds.variables['loss_time_vmec'][:],
            'loss_pos': ds.variables['loss_pos_vmec'][:],
            'lost_step': ds.variables['lost_step_vmec'][:]
        },
        'chartmap': {
            'loss_time': ds.variables['loss_time_chartmap'][:],
            'loss_pos': ds.variables['loss_pos_chartmap'][:],
            'lost_step': ds.variables['lost_step_chartmap'][:]
        },
        'boozer': {
            'loss_time': ds.variables['loss_time_boozer'][:],
            'loss_pos': ds.variables['loss_pos_boozer'][:],
            'lost_step': ds.variables['lost_step_boozer'][:]
        },
        'meiss_vmec': {
            'loss_time': ds.variables['loss_time_meiss_vmec'][:],
            'loss_pos': ds.variables['loss_pos_meiss_vmec'][:],
            'lost_step': ds.variables['lost_step_meiss_vmec'][:]
        },
        'meiss_chart': {
            'loss_time': ds.variables['loss_time_meiss_chart'][:],
            'loss_pos': ds.variables['loss_pos_meiss_chart'][:],
            'lost_step': ds.variables['lost_step_meiss_chart'][:]
        }
    }

    ds.close()
    return data


def compute_loss_curves(data):
    """Compute cumulative loss fraction over time."""
    n_particles = len(data['vmec']['loss_time'])
    trace_time = data['trace_time']
    n_steps = data['n_steps']

    time = np.linspace(0, trace_time, n_steps)
    curves = {}

    for name in ['vmec', 'chartmap', 'boozer', 'meiss_vmec', 'meiss_chart']:
        loss_times = data[name]['loss_time']
        lost_fraction = np.zeros(n_steps)
        for i, t in enumerate(time):
            lost_fraction[i] = np.sum(loss_times <= t) / n_particles
        curves[name] = {
            'time': time,
            'lost': lost_fraction,
            'confined': 1.0 - lost_fraction
        }

    return curves


def plot_loss_over_time(curves, output_file):
    """Plot loss and confined fractions over time."""
    plt = _maybe_import_matplotlib()
    if plt is None:
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    colors = {
        'vmec': 'black', 'chartmap': 'red', 'boozer': 'orange',
        'meiss_vmec': 'blue', 'meiss_chart': 'green'
    }
    labels = {
        'vmec': 'VMEC (RK45)',
        'chartmap': 'Chartmap (RK45)',
        'boozer': 'Boozer (RK45)',
        'meiss_vmec': 'Meiss from VMEC',
        'meiss_chart': 'Meiss from Chartmap'
    }

    ax = axes[0]
    for name in ['vmec', 'chartmap', 'boozer', 'meiss_vmec', 'meiss_chart']:
        curve = curves[name]
        ax.plot(curve['time'] * 1e6, curve['lost'] * 100,
                color=colors[name], linewidth=2, label=labels[name])
    ax.set_xlabel('Time [microseconds]')
    ax.set_ylabel('Lost Fraction [%]')
    ax.set_title('Particle Losses Over Time')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 100])

    ax = axes[1]
    for name in ['vmec', 'chartmap', 'boozer', 'meiss_vmec', 'meiss_chart']:
        curve = curves[name]
        ax.plot(curve['time'] * 1e6, curve['confined'] * 100,
                color=colors[name], linewidth=2, label=labels[name])
    ax.set_xlabel('Time [microseconds]')
    ax.set_ylabel('Confined Fraction [%]')
    ax.set_title('Particle Confinement Over Time')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 100])

    try:
        plt.tight_layout()
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
    except Exception as exc:
        print(f"Skipping plot_chartmap_wall_losses (savefig failed: {exc})")
        return
    print(f'Saved loss time plot: {output_file}')
    plt.close()


def plot_loss_positions(data, output_file):
    """Plot loss positions on wall in (theta, phi) coordinates."""
    plt = _maybe_import_matplotlib()
    if plt is None:
        return

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    titles = {
        'vmec': 'VMEC (RK45)',
        'chartmap': 'Chartmap (RK45)',
        'boozer': 'Boozer (RK45)',
        'meiss_vmec': 'Meiss from VMEC',
        'meiss_chart': 'Meiss from Chartmap'
    }

    trace_time = data['trace_time']
    modes = ['vmec', 'chartmap', 'boozer', 'meiss_vmec', 'meiss_chart']

    for i, name in enumerate(modes):
        ax = axes[i]
        loss_pos = data[name]['loss_pos']
        loss_time = data[name]['loss_time']

        lost_mask = loss_time < trace_time * 0.999
        if np.any(lost_mask):
            # NetCDF stores as (dim3, n_particles), Python reads as (n_particles, dim3)
            s_lost = loss_pos[lost_mask, 0]
            theta_lost = loss_pos[lost_mask, 1]
            phi_lost = loss_pos[lost_mask, 2]
            time_lost = loss_time[lost_mask]

            sc = ax.scatter(np.degrees(phi_lost), np.degrees(theta_lost),
                            c=time_lost * 1e6, cmap='viridis', s=50, alpha=0.8)
            plt.colorbar(sc, ax=ax, label='Loss time [us]')
        else:
            ax.text(0.5, 0.5, 'No losses', transform=ax.transAxes,
                    ha='center', va='center', fontsize=14)

        ax.set_xlabel('Toroidal angle phi [degrees]')
        ax.set_ylabel('Poloidal angle theta [degrees]')
        ax.set_title(f'{titles[name]}\n(n_lost = {np.sum(lost_mask)})')
        ax.grid(True, alpha=0.3)

    # Hide the 6th subplot (unused)
    axes[5].set_visible(False)

    try:
        plt.tight_layout()
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
    except Exception as exc:
        print(f"Skipping plot_chartmap_wall_losses (savefig failed: {exc})")
        return
    print(f'Saved loss positions plot: {output_file}')
    plt.close()


def plot_summary(data, curves, output_file):
    """Create summary comparison plot."""
    plt = _maybe_import_matplotlib()
    if plt is None:
        return

    fig, axes = plt.subplots(1, 2, figsize=(16, 5))

    names = ['vmec', 'chartmap', 'boozer', 'meiss_vmec', 'meiss_chart']
    labels = ['VMEC\n(RK45)', 'Chartmap\n(RK45)', 'Boozer\n(RK45)',
              'Meiss\n(VMEC)', 'Meiss\n(Chart)']
    colors = ['black', 'red', 'orange', 'blue', 'green']

    confined_frac = [curves[n]['confined'][-1] * 100 for n in names]

    ax = axes[0]
    bars = ax.bar(labels, confined_frac, color=colors, alpha=0.7)
    ax.set_ylabel('Confined Fraction [%]')
    ax.set_title('Final Confinement (t = {:.0f} us)'.format(
        data['trace_time'] * 1e6))
    ax.set_ylim([0, 100])
    ax.grid(True, alpha=0.3, axis='y')

    for bar, val in zip(bars, confined_frac):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{val:.1f}%', ha='center', va='bottom', fontsize=10)

    ax = axes[1]
    n_particles = len(data['vmec']['loss_time'])
    trace_time = data['trace_time']

    for i, name in enumerate(names):
        loss_times = data[name]['loss_time']
        lost_mask = loss_times < trace_time * 0.999
        if np.any(lost_mask):
            ax.hist(loss_times[lost_mask] * 1e6, bins=20, alpha=0.5,
                    color=colors[i], label=f'{labels[i].replace(chr(10), " ")}',
                    density=True)

    ax.set_xlabel('Loss Time [microseconds]')
    ax.set_ylabel('Density')
    ax.set_title('Distribution of Loss Times')
    ax.legend()
    ax.grid(True, alpha=0.3)

    try:
        plt.tight_layout()
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
    except Exception as exc:
        print(f"Skipping plot_chartmap_wall_losses (savefig failed: {exc})")
        return
    print(f'Saved summary plot: {output_file}')
    plt.close()


def print_statistics(data, curves):
    """Print loss statistics."""
    print('\nLoss Statistics:')
    print('=' * 60)

    trace_time = data['trace_time']
    n_particles = len(data['vmec']['loss_time'])

    labels = {
        'vmec': 'VMEC (RK45)',
        'chartmap': 'Chartmap',
        'boozer': 'Boozer',
        'meiss_vmec': 'Meiss(VMEC)',
        'meiss_chart': 'Meiss(Chart)'
    }

    for name in ['vmec', 'chartmap', 'boozer', 'meiss_vmec', 'meiss_chart']:
        loss_times = data[name]['loss_time']
        lost_mask = loss_times < trace_time * 0.999
        n_lost = np.sum(lost_mask)
        confined = curves[name]['confined'][-1] * 100

        print(f'  {labels[name]:14s}: {n_lost:3d}/{n_particles} lost '
              f'({100 - confined:.1f}%), {confined:.1f}% confined')

        if n_lost > 0:
            mean_loss = np.mean(loss_times[lost_mask])
            min_loss = np.min(loss_times[lost_mask])
            max_loss = np.max(loss_times[lost_mask])
            print(f'                 Loss time: '
                  f'min={min_loss*1e6:.2f}us, '
                  f'mean={mean_loss*1e6:.2f}us, '
                  f'max={max_loss*1e6:.2f}us')

    print('=' * 60)


def main():
    nc_file = 'chartmap_wall_losses.nc'
    if len(sys.argv) > 1:
        nc_file = sys.argv[1]

    output_prefix = 'chartmap_wall_losses'
    if nc_file != 'chartmap_wall_losses.nc':
        output_prefix = Path(nc_file).with_suffix('').name

    if not Path(nc_file).exists():
        print(f'Error: {nc_file} not found')
        print('Run test_chartmap_wall_losses first')
        sys.exit(1)

    print(f'Loading data from {nc_file}...')
    data = load_data(nc_file)

    print('Computing loss curves...')
    curves = compute_loss_curves(data)

    print_statistics(data, curves)

    print('\nGenerating plots...')
    plot_loss_over_time(curves, f'{output_prefix}_time.png')
    plot_loss_positions(data, f'{output_prefix}_positions.png')
    plot_summary(data, curves, f'{output_prefix}_summary.png')

    print('\nDone!')


if __name__ == '__main__':
    main()
