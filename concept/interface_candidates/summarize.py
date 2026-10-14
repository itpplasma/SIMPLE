"""Loss Wilson intervals, paired bootstrap differences, quantile intervals."""
from pathlib import Path
import csv
import gzip
import json
import math
import numpy as np

ROOT = Path("concept/interface_candidates")
RESULTS = ROOT / "results"

def wilson(k, n, z=1.959963984540054):
    p = k/n
    den = 1+z*z/n
    center = (p+z*z/(2*n))/den
    half = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [max(0., center-half), min(1., center+half)]

def binomial_cdf(k, n, p):
    if k < 0:
        return 0.
    if k >= n or p <= 0:
        return 1.
    if p >= 1:
        return 0.
    if p > .5:
        return 1-binomial_cdf(n-k-1, n, 1-p)
    term = (1-p)**n
    value = term
    for j in range(1, k+1):
        term *= (n-j+1)/j*p/(1-p)
        value += term
    return min(1., max(0., value))

def clopper_pearson(k, n, alpha=.025):
    def solve(k, target):
        lo, hi = 0., 1.
        for _ in range(56):
            mid = (lo+hi)/2
            if binomial_cdf(k, n, mid) > target:
                lo = mid
            else:
                hi = mid
        return (lo+hi)/2
    lower = 0. if k == 0 else solve(k-1, 1-alpha/2)
    upper = 1. if k == n else solve(k, alpha/2)
    return lower, upper

def paired_loss_interval(delta):
    # Each discordance probability has an exact 97.5% binomial interval.
    # Bonferroni gives >=95% simultaneous coverage without independence.
    plus = clopper_pearson(int(sum(delta > 0)), len(delta))
    minus = clopper_pearson(int(sum(delta < 0)), len(delta))
    return [plus[0]-minus[1], plus[1]-minus[0]]

def load(record):
    opener = gzip.open if record["file"].endswith(".gz") else open
    with opener(record["file"], "rt") as f:
        meta = f.readline().strip().lstrip("# ")
    times = dict(q.split("=") for q in meta.split(","))
    array = np.genfromtxt(record["file"], delimiter=",", names=True, skip_header=1)
    return array, float(times["wall"])

def quantile_ci(x, quantile, samples):
    estimates = np.quantile(x[samples], quantile, axis=1)
    return np.quantile(estimates, [.025, .975]).tolist()

def key(c):
    return c["case"], c["aout"], c["speed"], c["tend"], c["sample"]

def rms(x):
    return float(np.sqrt(np.mean(x*x)))

def summarize(prefix="main"):
    manifest = json.loads((RESULTS / f"{prefix}_manifest.json").read_text())
    if prefix == "main":
        for extra in ("long", "sampling"):
            if (RESULTS / f"{extra}_manifest.json").exists():
                manifest += json.loads((RESULTS / f"{extra}_manifest.json").read_text())
    references = {}
    arrays = {}
    for rec in manifest:
        data, wall = load(rec)
        arrays[rec["file"]] = (data, wall)
        if rec["config"]["method"] == "FO":
            k = key(rec["config"])
            if k not in references or rec["config"]["tol"] < references[k]["config"]["tol"]:
                references[k] = rec
    rng = np.random.default_rng(20261014)
    boot = rng.integers(0, 256, (4000, 256))
    summaries = []
    for rec in manifest:
        c = rec["config"]
        a, wall = arrays[rec["file"]]
        valid = a["failed"] == 0
        n, failed = len(a), int(sum(~valid))
        k = int(sum(a["lost"]))
        pnorm = 10*c["speed"]
        rev = (a["rev_position"] >= 0) & valid
        q = dict(config=c, file=rec["file"], lost=k, failed=failed,
                 loss=k/n, loss_ci=wilson(k, n),
                 loss_with_failed_bounds=[k/n, (k+failed)/n], wall=wall,
                 energy_rms=rms(a["emax"][valid]),
                 energy_max=float(max(a["emax"][valid], default=math.nan)),
                 pphi_rms=rms(a["pmax"][valid])/pnorm if c["case"] == "A" else None,
                 pphi_max=float(max(a["pmax"][valid], default=math.nan))/pnorm
                    if c["case"] == "A" else None,
                 reverse_position=rms(a["rev_position"][rev]) if sum(rev) else None,
                 reverse_velocity=rms(a["rev_velocity"][rev]) if sum(rev) else None,
                 reverse_count=int(sum(rev)), switches=float(np.mean(a["switches"])),
                 entries=int(sum(np.floor(a["switches"]/2))),
                 exits=int(sum(np.ceil(a["switches"]/2))),
                 returned_markers=int(sum(a["switches"] >= 2)),
                 crossings=float(np.mean(a["crossings"])),
                 fo_fraction=float(np.mean(a["fo_fraction"])),
                 root_failures=int(sum(a["root_failures"])),
                 switch_energy=float(max(a["switch_energy"])),
                 cost_per_time=wall/(n*c["tend"]))
        # Incomplete trajectories are not treated as a converged distribution.
        if not failed:
            q["excursion_quantiles"] = np.quantile(a["excursion"], [.1, .5, .9]).tolist()
            q["excursion_ci"] = [quantile_ci(a["excursion"], p, boot) for p in (.1, .5, .9)]
        ref = references.get(key(c))
        if ref is None:
            # Sampling refinements compare with default sampled reference too.
            ref = references.get((c["case"], c["aout"], c["speed"], c["tend"], .5))
        if ref is not None:
            ar, wr = arrays[ref["file"]]
            q["reference"] = ref["file"]
            # Costs use the same RK tolerance as the formulation comparison;
            # the tighter reference is used only for physical error estimates.
            cost_ref = next((r for r in manifest if key(r["config"]) == key(c)
                             and r["config"]["method"] == "FO"
                             and r["config"]["tol"] == 1e-9), ref)
            q["cost_reference"] = cost_ref["file"]
            q["cost_relative"] = wall/arrays[cost_ref["file"]][1]
            if not failed:
                losses = a["lost"]-ar["lost"]
                q["loss_delta"] = float(np.mean(losses))
                q["loss_delta_ci"] = paired_loss_interval(losses)
                q["loss_delta_bootstrap_ci"] = np.quantile(np.mean(losses[boot], axis=1), [.025, .975]).tolist()
                q["fate_disagreement"] = float(np.mean(losses != 0))
                q["fate_disagreement_ci"] = wilson(int(sum(losses != 0)), n)
                excursion_delta = a["excursion"]-ar["excursion"]
                q["excursion_delta"] = float(np.mean(excursion_delta))
                q["excursion_delta_ci"] = np.quantile(np.mean(excursion_delta[boot], axis=1), [.025, .975]).tolist()
                q["excursion_rmse"] = rms(excursion_delta)
                q["excursion_wasserstein"] = float(np.mean(abs(np.sort(a["excursion"])-np.sort(ar["excursion"]))))
                grid = np.unique(np.r_[a["excursion"], ar["excursion"]])
                ca = np.searchsorted(np.sort(a["excursion"]), grid, side="right")/n
                cb = np.searchsorted(np.sort(ar["excursion"]), grid, side="right")/n
                q["excursion_ks"] = float(max(abs(ca-cb)))
        summaries.append(q)
    (RESULTS / f"{prefix}_summary.json").write_text(json.dumps(summaries, indent=2, allow_nan=False)+"\n")
    columns = ["case", "aout", "speed", "method", "order", "din", "tol", "knob", "width", "sample",
               "tend", "lost", "failed", "loss", "ci_low", "ci_high", "energy_rms", "pphi_rms",
               "reverse_position", "reverse_velocity", "switches", "fo_fraction", "cost_relative",
               "excursion_wasserstein", "excursion_ks", "loss_delta", "delta_ci_low", "delta_ci_high"]
    with (RESULTS / f"{prefix}_summary.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for q in summaries:
            out = {k: q.get(k) for k in columns}
            out.update({k: q["config"].get(k) for k in q["config"]})
            out.update(ci_low=q["loss_ci"][0], ci_high=q["loss_ci"][1],
                       delta_ci_low=q.get("loss_delta_ci", [None, None])[0],
                       delta_ci_high=q.get("loss_delta_ci", [None, None])[1])
            writer.writerow({k: out.get(k) for k in columns})
    return summaries

if __name__ == "__main__":
    s = summarize()
    print(f"Summarized {len(s)} configurations; 95% Wilson/exact paired loss intervals; bootstrap excursion intervals.")
