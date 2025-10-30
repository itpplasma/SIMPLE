#!/usr/bin/env python3
"""Simple particle tracing example using pysimple API."""

import pysimple

vmec_file = '../test/test_data/wout.nc'

pysimple.init(vmec_file, deterministic=True, trace_time=5e-5, ntestpart=32)

particles = pysimple.sample_surface(32, s=0.3)

results = pysimple.trace_parallel(particles, integrator='midpoint')

n_lost = (results['loss_times'] < 5e-5).sum()
n_confined = (results['loss_times'] >= 5e-5).sum()
print(f"Traced 32 particles: {n_confined} confined, {n_lost} lost")
