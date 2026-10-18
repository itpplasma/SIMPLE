#!/usr/bin/env python3
"""Compare seeded CPU/GPU fraction curves without treating unresolved as losses."""
import argparse
import json
from pathlib import Path

import numpy as np

from efss_scaling import check_traced, physics


def compare(cpu, gpu):
    a = np.loadtxt(cpu / 'start.dat', ndmin=2)
    b = np.loadtxt(gpu / 'start.dat', ndmin=2)
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Nonfinite starting markers')
    np.testing.assert_array_equal(a, b)
    n = len(a)
    check_traced(cpu, n)
    check_traced(gpu, n)
    ta, ca, ua = physics(cpu)
    tb, cb, ub = physics(gpu)
    if ta.shape != tb.shape:
        raise ValueError('CPU/GPU time grid lengths differ')
    np.testing.assert_allclose(ta, tb, rtol=0, atol=1e-14)
    # Engineering threshold for aggregate fraction curves, not a bound on the
    # number of paired marker-status changes or a statistical confidence bound.
    tolerance = 2 / n
    confined_delta = float(np.max(np.abs(ca - cb)))
    unresolved_delta = float(np.max(np.abs(ua - ub)))
    result = dict(markers=n, identical_starts=True,
                  cpu_confined_final=float(ca[-1]), gpu_confined_final=float(cb[-1]),
                  cpu_unresolved_final=float(ua[-1]), gpu_unresolved_final=float(ub[-1]),
                  max_confined_delta=confined_delta,
                  max_unresolved_delta=unresolved_delta,
                  absolute_tolerance=tolerance,
                  passed=max(confined_delta, unresolved_delta) <= tolerance)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cpu', type=Path)
    parser.add_argument('gpu', type=Path)
    args = parser.parse_args()
    result = compare(args.cpu, args.gpu)
    print(json.dumps(result, indent=2))
    if not result['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
