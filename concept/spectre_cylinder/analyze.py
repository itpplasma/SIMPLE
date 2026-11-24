"""Statistics and figure: SIMPLE against the analytic Lundquist cylinder.

Usage: analyze.py FIELD100.dat FIELD50.dat RUN100 RUN50 OUT.png
FIELD*.dat are written by test_spectre_cylinder_field.x; RUN* directories
hold markers.dat (cylinder_markers.x) and ref.npz (reference_fo.py).
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp

from reference_fo import read_markers


def field_table(path):
    d = np.loadtxt(path)
    return dict(theta=d[:, 0], r=d[:, 1], owner=d[:, 2], b=d[:, 3], b_an=d[:, 5],
                p=d[:, 6], p_an=d[:, 8])


def mean_se(x):
    x = x[np.isfinite(x)]
    return x.mean(), x.std(ddof=1)/np.sqrt(len(x))


def runs(path):
    meta, d = read_markers(f"{path}/markers.dat")
    ref = np.load(f"{path}/ref.npz")
    lam, b_gc = d[:, 2], d[:, 8]
    rho_l = meta["ro0_cmG"]*1e-6*np.sqrt(1 - lam**2)/b_gc
    sim = dict(dist=np.abs(d[:, 7])/rho_l, crossed=d[:, 17] > 0,
               mu=np.where(d[:, 17] > 0, d[:, 19]/d[:, 9], np.nan),
               d_after=np.where(d[:, 17] > 0, d[:, 21], np.nan))
    ind = ref["independent"]
    nph = len(ind)//len(d)
    ref_ind = dict(dist=np.repeat(sim["dist"], nph), crossed=ind[:, 2] > 0,
                   mu=ind[:, 3], d_after=ind[:, 4])
    pair = ref["paired"]
    ref_pair = dict(crossed=pair[:, 0] > 0, mu=pair[:, 1], d_after=pair[:, 2])
    return sim, ref_ind, ref_pair, rho_l


def report(name, sim, ref, pair, rho_l):
    out = [f"== {name}: {len(sim['crossed'])} SIMPLE markers, "
           f"{len(ref['crossed'])} independent reference orbits"]
    for lab, s in (("SIMPLE", sim), ("reference", ref)):
        p = s["crossed"].mean()
        out.append(f"crossing fraction {lab:9s} {p:.4f} +- "
                   f"{np.sqrt(p*(1 - p)/len(s['crossed'])):.4f}")
    for key, unit in (("mu", "mu_after/mu_before"), ("d_after", "d_after [m]")):
        ms, es = mean_se(sim[key])
        mr, er = mean_se(ref[key])
        ks = ks_2samp(sim[key][np.isfinite(sim[key])], ref[key][np.isfinite(ref[key])])
        out.append(f"{unit:19s} mean SIMPLE {ms:.4f}+-{es:.4f}  reference "
                   f"{mr:.4f}+-{er:.4f}  diff {ms - mr:+.4f}+-{np.hypot(es, er):.4f}"
                   f"  KS p={ks.pvalue:.3f}")
    both = sim["crossed"] & pair["crossed"]
    out.append(f"paired: crossing flag agrees {np.mean(sim['crossed'] == pair['crossed']):.4f}; "
               f"on {both.sum()} common crossings median |d mu|/mu "
               f"{np.median(np.abs(sim['mu'][both]/pair['mu'][both] - 1)):.2e}, "
               f"median |d d_after|/rho_L "
               f"{np.median(np.abs(sim['d_after'][both] - pair['d_after'][both])/rho_l[both]):.2e}")
    out.append(f"rho_L/a range {rho_l.min():.4f} .. {rho_l.max():.4f}, median {np.median(rho_l):.4f}")
    return out


def main():
    f100, f50, run100, run50, png = sys.argv[1:6]
    t100, t50 = field_table(f100), field_table(f50)
    lines = []
    for name, t in (("R0=100", t100), ("R0=50", t50)):
        lines.append(f"{name}: max |B| rel dev {np.max(np.abs(t['b']/t['b_an'] - 1)):.3e}, "
                     f"max pitch dev {np.max(np.abs(t['p'] - t['p_an'])):.3e} rad")
    res = {}
    for name, path in (("R0=100", run100), ("R0=50", run50)):
        res[name] = runs(path)
        lines += report(name, *res[name])
    print("\n".join(lines))

    fig, ax = plt.subplots(2, 3, figsize=(13, 7.2), constrained_layout=True)
    colors = {0.0: "C0", 1.5708: "C1", 3.1416: "C2"}
    labels = {0.0: r"$\theta=0$ (outboard)", 1.5708: r"$\theta=\pi/2$",
              3.1416: r"$\theta=\pi$ (inboard)"}
    for th, c in colors.items():
        m = np.abs(t100["theta"] - th) < 1e-3
        for own in (1, 2):
            mm = m & (t100["owner"] == own)
            o = np.argsort(t100["r"][mm])
            kw = dict(label=labels[th]) if own == 1 else {}
            ax[0, 0].plot(t100["r"][mm][o], t100["b"][mm][o], c, lw=1.2, **kw)
            ax[0, 1].plot(t100["r"][mm][o], t100["p"][mm][o], c, lw=1.2)
    for own in (1, 2):
        mm = (np.abs(t100["theta"]) < 1e-3) & (t100["owner"] == own)
        o = np.argsort(t100["r"][mm])
        kw = dict(label="Lundquist (cylinder)") if own == 1 else {}
        ax[0, 0].plot(t100["r"][mm][o], t100["b_an"][mm][o], "k--", lw=1, **kw)
        ax[0, 1].plot(t100["r"][mm][o], t100["p_an"][mm][o], "k--", lw=1)
    ax[0, 0].set(xlabel="minor radius r [m]", ylabel="|B| [T]",
                 title="(a) |B|, SIMPLE (R0 = 100 m) vs analytic")
    ax[0, 0].legend(fontsize=8)
    ax[0, 1].set(xlabel="minor radius r [m]",
                 ylabel=r"pitch angle atan($B_\theta/B_z$) [rad]",
                 title="(b) field direction: jump at the sheet")
    for name, t, c in (("R0 = 100 m", t100, "C3"), ("R0 = 50 m", t50, "C4")):
        dev = t["b"]/t["b_an"] - 1
        ax[0, 2].plot(t["r"], dev*100, ".", ms=1.5, color=c,
                      label=f"{name}: max {np.max(np.abs(dev))*100:.2f} %")
    ax[0, 2].set(xlabel="minor radius r [m]", ylabel=r"$|B|_{SIMPLE}/|B|_{an} - 1$ [%]",
                 title="(c) deviation, all 16 poloidal angles")
    ax[0, 2].legend(fontsize=8, markerscale=6)
    for a_ in ax[0]:
        a_.axvline(0.5, color="0.6", lw=0.8, ls=":")

    sim, ref, pair, rho_l = res["R0=100"]
    bins = np.linspace(0, 1.4, 8)
    mid = 0.5*(bins[1:] + bins[:-1])
    for lab, s, c, off in (("SIMPLE hybrid", sim, "C0", -0.01), ("reference FO", ref, "k", 0.01)):
        idx = np.digitize(s["dist"], bins) - 1
        p = np.array([s["crossed"][idx == k].mean() for k in range(len(mid))])
        n = np.array([(idx == k).sum() for k in range(len(mid))])
        ax[1, 0].errorbar(mid + off, p, np.sqrt(p*(1 - p)/n), fmt="o", color=c,
                          ms=4, capsize=2, label=f"{lab}: {s['crossed'].mean():.3f}")
    ax[1, 0].set(xlabel=r"launch distance of the guiding centre $|d|/\rho_L$",
                 ylabel="fraction reaching the other volume",
                 title="(d) crossing fraction (R0 = 100 m)")
    ax[1, 0].legend(fontsize=8)
    hb = np.linspace(0, 3, 31)
    ax[1, 1].hist(ref["mu"][np.isfinite(ref["mu"])], hb, density=True, histtype="step",
                  color="k", label="reference FO")
    ax[1, 1].hist(sim["mu"][np.isfinite(sim["mu"])], hb, density=True, histtype="step",
                  color="C0", lw=1.5, label="SIMPLE hybrid")
    ax[1, 1].set(xlabel=r"$\mu_{after}/\mu_{before}$",
                 ylabel="probability density",
                 title=r"(e) magnetic moment after one passage")
    ax[1, 1].legend(fontsize=8)
    hd = np.linspace(-8, 8, 33)
    ax[1, 2].hist(100*ref["d_after"][np.isfinite(ref["d_after"])], hd, density=True,
                  histtype="step", color="k", label="reference FO")
    ax[1, 2].hist(100*sim["d_after"][np.isfinite(sim["d_after"])], hd, density=True,
                  histtype="step", color="C0", lw=1.5, label="SIMPLE hybrid")
    ax[1, 2].set(xlabel="guiding centre to sheet distance [cm]",
                 ylabel=r"probability density [cm$^{-1}$]",
                 title="(f) guiding-centre position after one passage")
    ax[1, 2].legend(fontsize=8)
    fig.suptitle("Two-volume Lundquist cylinder limit (a = 1 m, sheet at r = 0.5 m, "
                 r"3.5 MeV alphas, $\rho_L/a$ = 2.5-6 %; (d)-(f): 800 SIMPLE markers, 3200 reference orbits)", fontsize=11)
    fig.savefig(png, dpi=150)
    with open(png.rsplit(".", 1)[0] + ".txt", "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
