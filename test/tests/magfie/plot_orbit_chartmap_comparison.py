#!/usr/bin/env python3
"""
Plot orbit comparison: VMEC vs Meiss-VMEC vs Chartmap vs Meiss-Chartmap.

Reads orbit_chartmap_comparison.nc and creates 3D trajectory plots
comparing all four integration paths in Cartesian coordinates.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import netCDF4 as nc
import sys
from pathlib import Path


def load_data(filename):
    """Load trajectory data from NetCDF file."""
    ds = nc.Dataset(filename, 'r')

    data = {
        'time': ds.variables['time'][:],
        'vmec': {
            'x': ds.variables['x_vmec'][:],
            'y': ds.variables['y_vmec'][:],
            'z': ds.variables['z_vmec'][:]
        },
        'meiss_vmec': {
            'x': ds.variables['x_meiss_vmec'][:],
            'y': ds.variables['y_meiss_vmec'][:],
            'z': ds.variables['z_meiss_vmec'][:]
        },
        'chartmap': {
            'x': ds.variables['x_chartmap'][:],
            'y': ds.variables['y_chartmap'][:],
            'z': ds.variables['z_chartmap'][:]
        },
        'meiss_chartmap': {
            'x': ds.variables['x_meiss_chartmap'][:],
            'y': ds.variables['y_meiss_chartmap'][:],
            'z': ds.variables['z_meiss_chartmap'][:]
        }
    }

    ds.close()
    return data


def compute_deviations(data, reference='vmec'):
    """Compute Euclidean distance from reference trajectory."""
    ref = data[reference]
    deviations = {}

    for name in ['meiss_vmec', 'chartmap', 'meiss_chartmap']:
        traj = data[name]
        dev = np.sqrt(
            (traj['x'] - ref['x'])**2 +
            (traj['y'] - ref['y'])**2 +
            (traj['z'] - ref['z'])**2
        )
        deviations[name] = dev

    return deviations


def plot_3d_trajectories(data, output_file):
    """Create 3D plot of all trajectories."""
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')

    labels = {
        'vmec': 'VMEC (reference)',
        'meiss_vmec': 'Meiss from VMEC',
        'chartmap': 'Chartmap',
        'meiss_chartmap': 'Meiss from Chartmap'
    }

    colors = {
        'vmec': 'black',
        'meiss_vmec': 'blue',
        'chartmap': 'red',
        'meiss_chartmap': 'green'
    }

    linewidths = {
        'vmec': 2.5,
        'meiss_vmec': 1.5,
        'chartmap': 1.5,
        'meiss_chartmap': 1.5
    }

    for name in ['vmec', 'meiss_vmec', 'chartmap', 'meiss_chartmap']:
        traj = data[name]
        ax.plot(traj['x'], traj['y'], traj['z'],
               label=labels[name], color=colors[name],
               linewidth=linewidths[name], alpha=0.8)

    ax.scatter([data['vmec']['x'][0]], [data['vmec']['y'][0]], [data['vmec']['z'][0]],
              color='black', s=100, marker='o', label='Start')
    ax.scatter([data['vmec']['x'][-1]], [data['vmec']['y'][-1]], [data['vmec']['z'][-1]],
              color='black', s=100, marker='s', label='End (VMEC)')

    ax.set_xlabel('X [cm]')
    ax.set_ylabel('Y [cm]')
    ax.set_zlabel('Z [cm]')
    ax.set_title('Orbit Comparison: 4 Coordinate Systems\n(Cartesian Space)')
    ax.legend(loc='upper left')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f'Saved 3D trajectory plot: {output_file}')
    plt.close()


def plot_deviations(data, deviations, output_file):
    """Plot deviations from reference trajectory over time."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    time = data['time']

    ax = axes[0, 0]
    for name, dev in deviations.items():
        ax.semilogy(time, dev, label=name.replace('_', ' ').title())
    ax.set_xlabel('Time [normalized]')
    ax.set_ylabel('Deviation from VMEC [cm]')
    ax.set_title('Deviation in Cartesian Space')
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(time, data['vmec']['x'], 'k-', label='VMEC', linewidth=2)
    ax.plot(time, data['meiss_vmec']['x'], 'b--', label='Meiss-VMEC', alpha=0.7)
    ax.plot(time, data['chartmap']['x'], 'r--', label='Chartmap', alpha=0.7)
    ax.plot(time, data['meiss_chartmap']['x'], 'g--', label='Meiss-Chart', alpha=0.7)
    ax.set_xlabel('Time')
    ax.set_ylabel('X [cm]')
    ax.set_title('X Coordinate')
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(time, data['vmec']['y'], 'k-', label='VMEC', linewidth=2)
    ax.plot(time, data['meiss_vmec']['y'], 'b--', label='Meiss-VMEC', alpha=0.7)
    ax.plot(time, data['chartmap']['y'], 'r--', label='Chartmap', alpha=0.7)
    ax.plot(time, data['meiss_chartmap']['y'], 'g--', label='Meiss-Chart', alpha=0.7)
    ax.set_xlabel('Time')
    ax.set_ylabel('Y [cm]')
    ax.set_title('Y Coordinate')
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(time, data['vmec']['z'], 'k-', label='VMEC', linewidth=2)
    ax.plot(time, data['meiss_vmec']['z'], 'b--', label='Meiss-VMEC', alpha=0.7)
    ax.plot(time, data['chartmap']['z'], 'r--', label='Chartmap', alpha=0.7)
    ax.plot(time, data['meiss_chartmap']['z'], 'g--', label='Meiss-Chart', alpha=0.7)
    ax.set_xlabel('Time')
    ax.set_ylabel('Z [cm]')
    ax.set_title('Z Coordinate')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f'Saved deviation plot: {output_file}')
    plt.close()


def plot_final_comparison(data, deviations, output_file):
    """Create summary comparison plot."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    names = list(deviations.keys())
    max_devs = [np.max(deviations[n]) for n in names]
    mean_devs = [np.mean(deviations[n]) for n in names]

    x = np.arange(len(names))
    width = 0.35

    bars1 = ax.bar(x - width/2, max_devs, width, label='Max deviation', color='coral')
    bars2 = ax.bar(x + width/2, mean_devs, width, label='Mean deviation', color='steelblue')

    ax.set_ylabel('Deviation [cm]')
    ax.set_title('Trajectory Deviations from VMEC Reference')
    ax.set_xticks(x)
    ax.set_xticklabels([n.replace('_', '\n') for n in names])
    ax.legend()
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, axis='y')

    for bar, val in zip(bars1, max_devs):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
               f'{val:.2e}', ha='center', va='bottom', fontsize=8, rotation=45)

    ax = axes[1]
    ref = data['vmec']
    R_ref = np.sqrt(ref['x']**2 + ref['y']**2)
    phi_ref = np.arctan2(ref['y'], ref['x'])

    ax.plot(R_ref * np.cos(phi_ref), R_ref * np.sin(phi_ref), 'k-',
           linewidth=2, label='VMEC')

    for name, color in [('meiss_vmec', 'blue'), ('chartmap', 'red'),
                        ('meiss_chartmap', 'green')]:
        traj = data[name]
        R = np.sqrt(traj['x']**2 + traj['y']**2)
        phi = np.arctan2(traj['y'], traj['x'])
        ax.plot(R * np.cos(phi), R * np.sin(phi), '--',
               color=color, alpha=0.7, label=name.replace('_', ' '))

    ax.set_xlabel('X [cm]')
    ax.set_ylabel('Y [cm]')
    ax.set_title('Top View (X-Y Plane)')
    ax.legend()
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f'Saved summary plot: {output_file}')
    plt.close()


def main():
    nc_file = 'orbit_chartmap_comparison.nc'

    if not Path(nc_file).exists():
        print(f'Error: {nc_file} not found')
        print('Run test_orbit_chartmap_comparison first')
        sys.exit(1)

    print(f'Loading data from {nc_file}...')
    data = load_data(nc_file)

    print('Computing deviations...')
    deviations = compute_deviations(data)

    print('\nDeviation statistics (from VMEC reference):')
    for name, dev in deviations.items():
        print(f'  {name:20s}: max={np.max(dev):.4e} cm, mean={np.mean(dev):.4e} cm')

    print('\nGenerating plots...')
    plot_3d_trajectories(data, 'orbit_chartmap_3d.png')
    plot_deviations(data, deviations, 'orbit_chartmap_deviations.png')
    plot_final_comparison(data, deviations, 'orbit_chartmap_summary.png')

    print('\nDone!')


if __name__ == '__main__':
    main()
