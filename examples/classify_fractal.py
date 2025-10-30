#!/usr/bin/env python3
"""Fractal classification using fractal dimension."""

import pysimple

vmec_file = '../test/test_data/wout.nc'

pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, ntestpart=16)

particles = pysimple.sample_surface(16, s=0.4)

results = pysimple.classify_fractal(particles, tcut=0.1)

n_passing = results['passing'].sum()
n_trapped = (~results['passing']).sum()
n_lost = results['lost'].sum()

fractal_regular = (results['fractal'] == 1).sum()
fractal_chaotic = (results['fractal'] == 2).sum()

print(f"Fractal classification of {16} particles:")
print(f"  Passing: {n_passing}, Trapped: {n_trapped}, Lost: {n_lost}")
print(f"  Fractal: {fractal_regular} regular, {fractal_chaotic} chaotic")
