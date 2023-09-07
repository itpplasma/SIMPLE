#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 14:13:24 2019

@author: ert
"""
# %%
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.stats import gaussian_kde

#basedir = '/home/calbert/cobraruns/NEO-ORB/QA/2019-09-17_10k'

#basedir = '/home/calbert/net/cobra/run/SIMPLE/QA/2020-01-22_s0.3_withcut'
#basedir = '/home/calbert/net/cobra/run/SIMPLE/QH/2020-02-10_s0.3_withcut'
basedir = '/home/calbert/net/cobra/run/SIMPLE/QA/2020-02-10_s0.6_withcut'
#basedir = '/home/calbert/net/cobra/run/NEO-ORB/QH/2019-09-17_1k'

data = np.loadtxt(os.path.join(basedir, 'start.dat'))
r0 = data[:, 0]
th0 = np.mod(data[:, 1], 2*np.pi)
ph0 = np.mod(data[:, 2], 2*np.pi)
p0 = data[:, 3]
lam0 = data[:, 4]

# plt.figure()
#plt.plot(th0, ph0, ',')

# %%

npart = data.shape[0]

data = np.loadtxt(os.path.join(basedir, 'times_lost.dat'))
tlost = np.abs(data[:, 1])
trap_par = data[:, 2]

data = np.loadtxt(os.path.join(basedir, 'confined_fraction.dat'))
tconf = np.abs(data[1:, 0])
conffrac = data[1:, 1]+data[1:, 2]

# %%
plt.figure()
plt.scatter(th0, lam0, c=np.log10(tlost), s=0.5)
plt.colorbar()
# %%

nbins = 25
bins_lam = np.linspace(-0.9, 0.9, nbins-1)
bins_lam_center = np.linspace(-0.95, 0.95, nbins)
bin_idx = np.digitize(lam0, bins_lam)

t1 = 0.0
t2 = 0.99
lossfrac = np.zeros(nbins)
npart_bins = np.zeros(nbins)
for k in range(npart):
    if t1 < np.abs(tlost[k]):
        npart_bins[bin_idx[k]] += 1.0
    if t1 < np.abs(tlost[k]) and np.abs(tlost[k]) < t2:
        lossfrac[bin_idx[k]] += 1.0

lossfrac = lossfrac/npart_bins


kde2 = gaussian_kde([np.log10(tlost[tlost < t2]), trap_par[tlost < t2]])
X, Y = np.mgrid[-5:0:100j, -2:1:200j]
positions = np.vstack([X.ravel(), Y.ravel()])
Z = np.reshape(kde2(positions).T, X.shape)

plt.figure()
plt.imshow(np.rot90(Z), cmap=plt.cm.Blues, extent=[-5, 0, -2, 1])
plt.plot(np.log10(tlost), trap_par, 'kx', markersize=0.3, alpha=0.5)
#plt.plot(np.log10(tlost[wrong]), trap_par[wrong], 'rx', markersize=2, alpha=1)
plt.plot([-5, 0], [0, 0], 'k--', linewidth=.8)
plt.ylim([-2.1, 1.0])
plt.xlim(-5.0, 0.1)
plt.ylabel(r'trapping parameter')
plt.xlabel(r'$t\,/\,\mathrm{s}$')
plt.gca().set_xticklabels([r'$10^{-5}$', r'$10^{-4}$',
                           r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$1$'])
ax2 = plt.gca().twinx()
ax2.plot(np.log10(tconf), conffrac, color='tab:red')
ax2.set_ylim([0.5, 1.25])
ax2.set_ylabel('confined fraction                    ', color='tab:red')
ax2.set_yticks(np.linspace(0.5, 1, 6))
ax2.tick_params('y', colors='tab:red')
# plt.savefig('qi.png',dpi=150)
