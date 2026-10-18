#!/usr/bin/env python3
"""Plot marker scaling and paired CPU/GPU throughput from raw repeated timings."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu', nargs='+', type=Path, required=True)
    parser.add_argument('--gpu', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rows = np.concatenate([np.atleast_1d(np.genfromtxt(
        path, delimiter=',', names=True, dtype=None, encoding='utf-8'))
        for path in args.cpu])
    gpu = np.atleast_1d(np.genfromtxt(
        args.gpu, delimiter=',', names=True, dtype=None, encoding='utf-8'))
    rows = rows[rows['threads'] == 16]
    counts = np.unique(rows['markers'])
    groups = [rows[rows['markers'] == n] for n in counts]
    wall = np.array([np.mean(group['wall_s']) for group in groups])
    trace = np.array([np.mean(group['trace_s']) for group in groups])
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.3), layout='constrained')
    axes[0].plot(counts, wall, 'o-', color='#0072B2', label='Total wall time')
    axes[0].plot(counts, trace, 's--', color='black', label='Tracing phase')
    axes[0].set(xlabel='Markers', ylabel='Time (s)', title='CPU, 16 OpenMP threads')
    axes[0].legend(fontsize=8)
    axes[0].set_xscale('log', base=2)
    axes[0].set_xticks(counts, labels=[str(n) for n in counts])
    largest = counts[-1]
    if not np.all(gpu['markers'] == largest):
        raise ValueError('GPU and largest CPU marker counts must match')
    gpu_wall = np.mean(gpu['wall_s'])
    axes[1].bar(['CPU 16 threads', 'CUDA 1 GPU'],
                [largest / wall[-1], largest / gpu_wall],
                color=['#0072B2', '#E69F00'], hatch=['', '//'], edgecolor='black')
    axes[1].set(ylabel='Markers / total wall second', title=f'{largest} identical starts')
    for axis in axes:
        axis.grid(axis='y', alpha=.25)
        axis.set_axisbelow(True)
    fig.suptitle('SIMPLE: repeated means, including setup and output')
    fig.savefig(args.output, dpi=130)


if __name__ == '__main__':
    main()
