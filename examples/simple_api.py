#!/usr/bin/env python3
"""Simple examples using pysimple API."""

import pysimple

vmec_file = '../test/test_data/wout.nc'

pysimple.init(vmec_file, deterministic=True, trace_time=5e-5, ntestpart=32)

particles = pysimple.sample_surface(32, s=0.3)

results = pysimple.trace_parallel(particles, integrator='midpoint')

n_lost = (results['loss_times'] < 5e-5).sum()
print(f"Traced 32 particles: {n_lost} lost")

pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, ntestpart=16)

particles = pysimple.sample_surface(16, s=0.4)

results = pysimple.classify_fast(particles, nturns=8)

n_passing = results['passing'].sum()
n_trapped = (~results['passing']).sum()
print(f"Fast classification of 16 particles: {n_passing} passing, {n_trapped} trapped")
