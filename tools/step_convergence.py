#!/usr/bin/env python3
"""Statistical step-size convergence check for the symplectic integrator.

Symplectic integrators have no per-step error estimator, and individual chaotic
orbits diverge exponentially between step sizes, so a trajectory-level comparison
is meaningless. The loss statistics, however, are a property of the coarse-grained
transport topology and are robust to local chaos. This tool therefore compares two
runs at different ``npoiper2`` (same initial conditions, i.e. the same start.dat)
at the level of the loss statistic:

1. Confined-fraction bias from the paired outcome flips (McNemar). Only orbits
   whose lost/confined outcome actually changes between the two step sizes
   contribute, so local chaos that does not change the outcome is ignored.
2. Loss-time distribution distance (Kolmogorov-Smirnov and 1-Wasserstein) over
   the lost particles, compared unpaired (tolerant to chaotic loss-time jitter).

The coarse step is adequate when the bias and the distribution distance are below
the Monte-Carlo error of the confined fraction (~1/sqrt(N)). Evaluate this at the
full trace time only: a coarse step can look converged early and diffuse later.

times_lost.dat convention: column 2 is the loss time. A particle is confined if
it is either never traced (sentinel < 0) or survives to the trace time
(loss_time >= trace_time); it is lost only if 0 < loss_time < trace_time. The
trace time is read from the simple.in next to each times_lost.dat.
"""
import argparse
import os
import re
import numpy as np


def read_times_lost(path):
    """Return (index, loss_time) arrays from a times_lost.dat file."""
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[None, :]
    return data[:, 0].astype(int), data[:, 1]


def read_trace_time(times_lost_path):
    """trace_time from the simple.in sibling of a times_lost.dat path."""
    cfg = os.path.join(os.path.dirname(os.path.abspath(times_lost_path)), "simple.in")
    m = re.search(r"^\s*trace_time\s*=\s*([0-9.eEdD+-]+)", open(cfg).read(), re.M)
    if not m:
        raise SystemExit(f"trace_time not found in {cfg}")
    return float(m.group(1).replace("d", "e").replace("D", "e"))


def confined_mask(loss_time, trace_time):
    """Confined = never traced (loss_time < 0) or survived to trace_time.
    Lost = 0 < loss_time < trace_time."""
    return (loss_time < 0.0) | (loss_time >= trace_time * (1.0 - 1e-9))


def ks_statistic(a, b):
    """Two-sample Kolmogorov-Smirnov statistic, no SciPy dependency."""
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    grid = np.sort(np.concatenate([a, b]))
    ca = np.searchsorted(np.sort(a), grid, side="right") / len(a)
    cb = np.searchsorted(np.sort(b), grid, side="right") / len(b)
    return float(np.max(np.abs(ca - cb)))


def wasserstein1(a, b):
    """1-Wasserstein distance between two empirical 1-D distributions."""
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    qs = np.linspace(0.0, 1.0, 1024)
    qa = np.quantile(a, qs)
    qb = np.quantile(b, qs)
    return float(np.mean(np.abs(qa - qb)))


def compare(coarse_path, fine_path):
    ic, tc = read_times_lost(coarse_path)
    iff, tf = read_times_lost(fine_path)
    # Align on the shared particle indices (same initial conditions).
    common = np.intersect1d(ic, iff)
    if len(common) == 0:
        raise SystemExit("no shared particle indices; runs must share start.dat")
    tc = tc[np.searchsorted(ic, common)]
    tf = tf[np.searchsorted(iff, common)]
    n = len(common)

    conf_c = confined_mask(tc, read_trace_time(coarse_path))
    conf_f = confined_mask(tf, read_trace_time(fine_path))
    cf_c = conf_c.mean()
    cf_f = conf_f.mean()

    # Paired outcome flips (McNemar): discordant pairs only.
    n_cf = int(np.sum(conf_c & ~conf_f))   # confined coarse, lost fine
    n_fc = int(np.sum(~conf_c & conf_f))   # lost coarse, confined fine
    bias = (n_cf - n_fc) / n               # cf_c - cf_f
    flip_rate = (n_cf + n_fc) / n
    se_bias = np.sqrt(n_cf + n_fc) / n      # McNemar paired standard error

    # Monte-Carlo error of the confined fraction itself.
    p = 0.5 * (cf_c + cf_f)
    mc_err = np.sqrt(max(p * (1.0 - p), 0.0) / n)

    # Loss-time distribution distance over lost particles (unpaired).
    lost_c = tc[~conf_c]
    lost_f = tf[~conf_f]
    ks = ks_statistic(lost_c, lost_f)
    w1 = wasserstein1(lost_c, lost_f)

    converged = abs(bias) <= mc_err and flip_rate <= 3.0 * mc_err
    return dict(n=n, cf_coarse=cf_c, cf_fine=cf_f, bias=bias, se_bias=se_bias,
                flip_rate=flip_rate, mc_err=mc_err, ks=ks, wasserstein=w1,
                n_cf=n_cf, n_fc=n_fc, converged=converged)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("coarse", help="times_lost.dat of the coarser (larger step) run")
    ap.add_argument("fine", help="times_lost.dat of the finer (smaller step) run")
    a = ap.parse_args()
    r = compare(a.coarse, a.fine)
    print(f"particles compared        : {r['n']}")
    print(f"confined fraction coarse  : {r['cf_coarse']:.4f}")
    print(f"confined fraction fine    : {r['cf_fine']:.4f}")
    print(f"bias (coarse - fine)      : {r['bias']:+.4f} +- {r['se_bias']:.4f} (paired)")
    print(f"outcome flip rate         : {r['flip_rate']:.4f}  ({r['n_cf']} conf->lost, {r['n_fc']} lost->conf)")
    print(f"Monte-Carlo error of c.f. : {r['mc_err']:.4f}")
    print(f"loss-time KS distance     : {r['ks']:.4f}")
    print(f"loss-time 1-Wasserstein   : {r['wasserstein']:.4g}")
    print(f"VERDICT                   : {'CONVERGED' if r['converged'] else 'NOT CONVERGED'} "
          f"(bias and flips within Monte-Carlo error)")
    raise SystemExit(0 if r["converged"] else 2)


if __name__ == "__main__":
    main()
