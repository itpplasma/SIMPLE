#!/usr/bin/env python3
"""Plot measured mean wall time and derived strong scaling, with min/max ranges."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('measurements', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    rows = np.genfromtxt(args.measurements, delimiter=',', names=True, dtype=None,
                         encoding='utf-8')
    threads = np.unique(rows['threads'])
    walls = [rows['wall_s'][rows['threads'] == t] for t in threads]
    means = np.array([np.mean(values) for values in walls])
    bounds = np.array([[np.min(values), np.max(values)] for values in walls]).T
    speedup = means[0] / means
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.3), layout='constrained')
    # Okabe-Ito blue and black, with markers/dashes as redundant encoding.
    axes[0].errorbar(threads, means, yerr=[means - bounds[0], bounds[1] - means],
                     color='#0072B2', marker='o', capsize=3)
    axes[0].set_ylabel('Wall time (s)')
    axes[1].plot(threads, speedup, 'o-', color='#0072B2', label='Measured mean ratio')
    axes[1].plot(threads, threads, '--', color='black', label='Ideal')
    axes[1].set_ylabel('Speedup vs 1 thread')
    axes[1].legend(fontsize=8)
    axes[2].plot(threads, 100 * speedup / threads, 'o-', color='#0072B2')
    axes[2].set_ylabel('Parallel efficiency (%)')
    axes[2].set_ylim(0, max(105, 1.05 * np.max(100 * speedup / threads)))
    for axis in axes:
        axis.set_xscale('log', base=2)
        axis.set_xticks(threads, labels=[str(t) for t in threads])
        axis.set_xlabel('OpenMP threads')
        axis.grid(alpha=.25)
    fig.suptitle(f'SIMPLE: {int(rows["markers"][0])} markers, whole-process wall time')
    fig.savefig(args.output, dpi=130)


if __name__ == '__main__':
    main()
