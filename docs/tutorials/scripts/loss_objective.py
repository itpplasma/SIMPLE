"""Alpha-particle loss fraction as an objective function with pysimple.

Usage: python loss_objective.py WOUT_QA WOUT_QH
"""
import sys
import time

import numpy as np
import pysimple

TRACE_TIME = 1e-2  # s
N_PARTICLES = 128
S_START = 0.6


def loss_fraction(wout, **settings):
    """Return the loss fraction and its standard error for one equilibrium.

    ``settings`` are further &config keys. The same seed gives the same start
    points for every call (common random numbers).
    """
    pysimple.init(wout, isw_field_type=2, trace_time=TRACE_TIME,
                  deterministic=True, **settings)
    start = pysimple.sample_surface(N_PARTICLES, s=S_START)
    assert np.allclose(start[0], S_START), "start points are not on the surface"
    result = pysimple.trace_parallel(start)
    t = result["loss_times"]
    lost = np.isfinite(t) & (t > 0) & (t < TRACE_TIME)
    n_resolved = np.isfinite(t).sum()
    if n_resolved == 0:
        raise RuntimeError("no numerically resolved particles")
    f = lost.sum() / n_resolved
    return f, np.sqrt(f * (1 - f) / n_resolved)


if __name__ == "__main__":
    wout_qa, wout_qh = sys.argv[1], sys.argv[2]
    for scale in (0.5, 0.75, 1.0):
        t0 = time.perf_counter()
        f, err = loss_fraction(wout_qa, vmec_B_scale=scale)
        print(f"QA, B x {scale:4.2f}: loss fraction {f:.3f} +- {err:.3f}"
              f"  ({time.perf_counter() - t0:.0f} s)")
    f, err = loss_fraction(wout_qh)
    print(f"QH, B x 1.00: loss fraction {f:.3f} +- {err:.3f}")
