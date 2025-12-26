#!/usr/bin/env python3
"""Trace orbit with trajectory output."""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "python"))
sys.modules.pop("pysimple", None)

import pysimple

vmec_file = '../test/test_data/wout.nc'

pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, ntestpart=1)

particles = pysimple.sample_surface(1, s=0.35)
particle = particles[:, 0]

result = pysimple.trace_orbit(particle, integrator='midpoint', return_trajectory=True)

n_timesteps = result['times'].shape[0]
print(f"Traced 1 particle with {n_timesteps} timesteps")
print(f"Loss time: {result['loss_time']:.3e}s")
