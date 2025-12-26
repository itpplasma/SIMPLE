#!/usr/bin/env python3
"""Fast classification using jpar and topology."""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "python"))
sys.modules.pop("pysimple", None)

import pysimple

vmec_file = '../test/test_data/wout.nc'

pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, ntestpart=16)

particles = pysimple.sample_surface(16, s=0.4)

results = pysimple.classify_fast(particles, nturns=8)

n_passing = results['passing'].sum()
n_trapped = (~results['passing']).sum()

jpar_regular = (results['jpar'] == 1).sum()
jpar_stochastic = (results['jpar'] == 2).sum()

topology_ideal = (results['topology'] == 1).sum()
topology_non_ideal = (results['topology'] == 2).sum()

print(f"Fast classification of 16 particles: {n_passing} passing, {n_trapped} trapped")
print(f"  J_parallel: {jpar_regular} regular, {jpar_stochastic} stochastic")
print(f"  Topology: {topology_ideal} ideal, {topology_non_ideal} non-ideal")
