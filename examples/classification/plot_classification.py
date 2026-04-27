import numpy as np
import matplotlib.pyplot as plt


def bin_classification(s, perp_inv, iclass, ns=50, nperp=100):
    """Bin particles into (s, perp_inv) grid by classification type.

    Returns prompt, regular, stochastic arrays of shape (nperp, ns).
    """
    hs = 1.0 / ns
    pmin = 0.0
    pmax = np.max(perp_inv)
    hp = (pmax - pmin) / nperp

    prompt = np.zeros((nperp, ns))
    regular = np.zeros((nperp, ns))
    stochastic = np.zeros((nperp, ns))

    for ipart in range(len(s)):
        i = min(ns, max(1, int(np.ceil(s[ipart] / hs)))) - 1
        k = min(nperp, max(1, int(np.ceil(perp_inv[ipart] / hp)))) - 1

        if iclass[ipart] == 0:
            prompt[k, i] += 1.0
        elif iclass[ipart] == 1:
            regular[k, i] += 1.0
        elif iclass[ipart] == 2:
            stochastic[k, i] += 1.0

    return prompt, regular, stochastic


def doplot_inner(prompt, regular, stochastic, bminmax):
    tot = prompt + regular + stochastic
    wpr = -1.0
    wst = 0.0
    wre = 1.0
    cla = wre * regular + wst * stochastic + wpr * prompt
    cla[tot == 0] = np.nan
    tot[tot == 0] = np.nan
    cla = cla / tot

    bmin_global = np.min(bminmax[:, 1])
    bmax_global = np.max(bminmax[:, 2])
    jpmin = bmin_global / bmax_global

    plt.figure(figsize=(5, 4))
    plt.imshow(cla, origin='lower', extent=[0, 1, 0, 1], aspect='auto',
               vmin=-1.0, vmax=1.0, clim=(-1.0, 1.0))
    plt.plot(bminmax[:, 0], bmin_global / bminmax[:, 1], 'k')
    plt.plot(bminmax[:, 0], bmin_global / bminmax[:, 2], 'k--')

    # Vertical line at s=0.25
    idx_025 = np.argmin(np.abs(bminmax[:, 0] - 0.25))
    plt.plot([0.25, 0.25], [0.0, bmin_global / bminmax[idx_025, 1]],
             'k', lw=0.75)

    plt.xticks([0, 0.25, 0.5, 0.75, 1])
    plt.ylim(jpmin, 1)
    plt.xlabel(r'Normalized toroidal flux $s$')
    plt.ylabel(r'Perpendicular invariant $J_\perp$')
    plt.tight_layout()


def main():
    # Read class_parts.dat: columns are
    # ipart, s, perp_inv, iclass_jpar, iclass_ideal, iclass_fractal
    data = np.loadtxt('class_parts.dat')
    s = data[:, 1]
    perp_inv = data[:, 2]
    icl_jpar = data[:, 3].astype(int)
    icl_ideal = data[:, 4].astype(int)

    bminmax = np.loadtxt('bminmax.dat')

    # J_parallel classification
    prompt1, regular1, stochastic1 = bin_classification(
        s, perp_inv, icl_jpar)
    doplot_inner(prompt1, regular1, stochastic1, bminmax)
    plt.title(r'$J_\parallel$ classifier')
    cb = plt.colorbar()
    cb.set_ticks([-1, 0, 1])
    cb.set_ticklabels(['Prompt losses', 'Non-ideal', 'Regular'])
    plt.savefig('class_jpar.pdf')

    # Topological (ideal/non-ideal) classification
    prompt2, regular2, stochastic2 = bin_classification(
        s, perp_inv, icl_ideal)
    doplot_inner(prompt2, regular2, stochastic2, bminmax)
    plt.title('Topological classifier')
    cb = plt.colorbar()
    cb.set_ticks([-1, 0, 1])
    cb.set_ticklabels(['Prompt losses', 'Non-ideal', 'Ideal'])
    plt.savefig('class_ideal.pdf')


if __name__ == '__main__':
    main()
