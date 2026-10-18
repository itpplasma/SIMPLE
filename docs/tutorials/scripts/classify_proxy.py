"""Fast regular/chaotic classification of trapped orbits with pysimple.

Usage: python classify_proxy.py WOUT
"""
import sys
import time

import numpy as np
import pysimple

pysimple.init(sys.argv[1], isw_field_type=2, trace_time=1e-2, deterministic=True)
start = pysimple.sample_surface(128, s=0.6)
t0 = time.perf_counter()
c = pysimple.classify_fast(start, nturns=8)
trapped = ~c["passing"]
print(f"classified in {time.perf_counter() - t0:.0f} s")
print(f"trapped: {trapped.sum()}, J_par regular: {(c['jpar'][trapped] == 1).sum()}, "
      f"chaotic: {(c['jpar'][trapped] == 2).sum()}, "
      f"unclassified: {(c['jpar'][trapped] == 0).sum()}")
print(f"lost during classification: {c['lost'].sum()}")
