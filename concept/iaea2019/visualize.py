#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 14:13:24 2019

@author: ert
"""
# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import glob
import os
from scipy.stats import gaussian_kde
import exportfig

# basedir_cut = '/home/calbert/cobraruns/NEO-ORB/QA/2019-09-18_cut'
#basedir = '/home/calbert/net/cobra/run/NEO-ORB/QI/2019-09-17_1k'
#basedir_cut = '/home/calbert/net/cobra/run/NEO-ORB/QI/2019-09-18_withcut'

basedir = '/home/calbert/net/cobra/run/SIMPLE/QA/2020-02-10_s0.3'
basedir_cut = basedir+'_withcut'
basedir_plot = basedir+'_cut'
exportname = 'qa_s03'


data = np.loadtxt(os.path.join(basedir, 'start.dat'))
r0 = data[:, 0]
th0 = np.mod(data[:, 1], 2*np.pi)
ph0 = np.mod(data[:, 2], 2*np.pi)
p0 = data[:, 3]
lam0 = data[:, 4]

npart = data.shape[0]

data = np.loadtxt(os.path.join(basedir, 'times_lost.dat'))
tlost = np.abs(data[:, 1])
trap_par = data[:, 2]

data = np.loadtxt(os.path.join(basedir, 'confined_fraction.dat'))
tconf = np.abs(data[1:, 0])
conffrac = data[1:, 1]+data[1:, 2]


filepath = os.path.join(basedir_cut, 'tjob.out*')
outfile = glob.glob(filepath)[0]
with open(os.path.join(basedir_cut, outfile)) as f:
    outtxt = f.readlines()

ntip = np.zeros(npart, dtype=int)
nper = np.zeros(npart, dtype=int)
regtip = np.zeros(npart, dtype=bool)
regper = np.zeros(npart, dtype=bool)

for line in outtxt:
    linespl = line.split()
    if len(linespl) < 2:
        continue
    if linespl[-2] == 'tip':
        ntip[int(linespl[0])-1] = int(linespl[-1])
        if 'regular' in line:
            regtip[int(linespl[0])-1] = True
    elif linespl[-2] == 'per':
        nper[int(linespl[0])-1] = int(linespl[-1])
        if 'regular' in line:
            regper[int(linespl[0])-1] = True

kregtip = np.where(regtip)[0]
kregper = np.where(regper)[0]
kreg = np.where(np.logical_and(regtip, regper))[0]

wrongtip = kregtip[np.logical_and(tlost[kregtip] < 1.0, tlost[kregtip] > 1e-1)]
wrongper = kregper[np.logical_and(tlost[kregper] < 1.0, tlost[kregper] > 1e-1)]
wrong = kreg[np.logical_and(tlost[kreg] < 1.0, tlost[kreg] > 1e-1)]

print(basedir)
print('wrong: {}'.format(wrong+1))


# plt.figure()
# plt.scatter(th0, ph0, c=lam0)
# plt.colorbar()

# plt.figure()
# plt.scatter(th0, ph0, c=np.log10(tlost))
# plt.colorbar()

plt.figure()
plt.scatter(th0, lam0, c=tlost)
plt.colorbar()

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


# plt.figure()
# plt.bar(bins_lam_center, lossfrac, width=2.0/nbins)
# plt.xlabel(r'pitch angle $v_\parallel/v$')
# plt.ylabel(r'loss fraction')
# plt.title(r'total losses')

kde = gaussian_kde(np.log10(tlost[tlost < t2]))
xpl = np.linspace(-5, 0, 100)

# plt.figure()
# plt.hist(np.log10(tlost[tlost < t2]), bins=12)
# plt.plot(xpl, kde(xpl)*len(data)/12)
# plt.xlabel(r'$\log_{10}(t_{\mathrm{\,lost}}\,/\,\mathrm{s})$')
# plt.ylabel('number of particles')
# plt.title('distribution of loss times')

# plt.figure()
# plt.hist2d(np.log10(tlost[tlost < t2]), trap_par[tlost < t2],
#            range=[[-5, 0], [-1, 1]], bins=[15, 15])
# plt.colorbar()

kde2 = gaussian_kde([np.log10(tlost[tlost < t2]), trap_par[tlost < t2]])
X, Y = np.mgrid[-5:0:100j, -2:1:200j]
positions = np.vstack([X.ravel(), Y.ravel()])
Z = np.reshape(kde2(positions).T, X.shape)

plt.figure(figsize=(2.9, 2.4))
plt.imshow(np.rot90(Z), cmap=plt.cm.Blues,
           extent=[-5, 0, -2, 1], aspect='auto')
plt.plot(np.log10(tlost), trap_par, 'k,', marker=',', color='k', markersize=0.1, alpha=0.5)
#plt.plot(np.log10(tlost[wrong]), trap_par[wrong], 'rx', markersize=2, alpha=1)
plt.plot([-5, 0], [0, 0], 'k--', linewidth=.8)
plt.ylim([-2.0, 1.0])
plt.xlim(-5.0, 0.05)
plt.ylabel(r'trapping parameter')
plt.xlabel(r'$t\,/\,\mathrm{s}$')
ax = plt.gca()
ax.set_xticks([-5, -4, -3, -2, -1, 0])
ax.set_xticklabels([r'$10^{-5}$', r'$10^{-4}$',
                    r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$1$'])
ax.set_xticks(np.concatenate([np.log10(np.linspace(1e-5, 1e-4, 10)),
                              np.log10(np.linspace(1e-4, 1e-3, 10)),
                              np.log10(np.linspace(1e-3, 1e-2, 10)),
                              np.log10(np.linspace(1e-2, 1e-1, 10)),
                              np.log10(np.linspace(1e-1, 1, 10))]),
              minor=True)

sig = np.sqrt(conffrac*(1.0-conffrac)/1000.0)

ax2 = plt.gca().twinx()
shift=-1000
ax2.plot(np.log10(tconf[:shift]), conffrac[:shift], color='tab:red')
ax2.fill_between(np.log10(tconf[:shift]), conffrac[:shift]-1.96*sig[:shift], 
                 conffrac[:shift]+1.96*sig[:shift], color='#FFBBBB')



ax2.set_xlim(-5.0, 0.05)
ax2.set_ylim([0.6, 1.2])
ax2.set_ylabel(r'confined fraction                    ', color='tab:red')
ax2.set_yticks(np.linspace(0.6, 1, 5))
ax2.tick_params('y', colors='tab:red')
plt.tight_layout()

exportfig.exportfig(exportname)
exportfig.exporteps(exportname)

# %%

# th = np.linspace(0, 2*np.pi, 100)


# for wt in wrong:
#     plt.figure(figsize=(4, 4))
#     data = np.loadtxt(os.path.join(basedir_cut, 'fort.{}'.format(10000+wt+1)))
#     plt.plot(np.cos(th), np.sin(th), 'k')
#     plt.plot(np.sqrt(data[:, 0])*np.cos(data[:, 1]),
#              np.sqrt(data[:, 0])*np.sin(data[:, 1]), ',', color='tab:red')
#     plt.xlabel(r'$R-R_0$')
#     plt.ylabel(r'$Z$')
#     plt.title('tip cut')
#     plt.axis('equal')
#     plt.figure(figsize=(4, 4))
#     data = np.loadtxt(os.path.join(basedir_cut, 'fort.{}'.format(20000+wt+1)))
#     plt.plot(np.cos(th), np.sin(th), 'k')
#     plt.plot(np.sqrt(data[:, 0])*np.cos(data[:, 1]),
#              np.sqrt(data[:, 0])*np.sin(data[:, 1]), ',', color='tab:red')
#     plt.xlabel(r'$R-R_0$')
#     plt.ylabel(r'$Z$')
#     plt.title('periodic cut')
#     plt.axis('equal')
#     plt.tight_layout()

# %%
wt = 24
#wt = 69
#wt = 882
#wt = 982

data = np.loadtxt(os.path.join(basedir_plot, 'fort.{}'.format(20000+wt+1)))

x = data[:, 0]*np.cos(data[:, 1])
y = data[:, 0]*np.sin(data[:, 1])

x = 0.05+0.9*(x - np.min(x))/(np.max(x) - np.min(x))
y = 0.05+0.9*(y - np.min(y))/(np.max(y) - np.min(y))

#plt.plot(np.cos(th), np.sin(th), 'k')
# plt.plot(np.sqrt(data[:, 0])*np.cos(data[:, 1]),
#         np.sqrt(data[:, 0])*np.sin(data[:, 1]), ',', color='tab:red')

#plt.title('periodic cut')

nstep = 5
eps = 1.0/(nstep-1)

phase = 0.0

plt.figure(figsize=(3.2, 3.2))
for k in range(nstep-1):
    for l in range(nstep-1):
        rn = x
        tn = y
        t = np.where((rn >= k*eps) & (rn < (k+1)*eps)
                     & (tn >= l*eps) & (tn < (l+1)*eps))
        if t[0].size > 0:
            plt.fill_between(np.array([k*eps, (k+1)*eps]),
                             np.array([l*eps, l*eps]),
                             np.array([(l+1)*eps, (l+1)*eps]),
                             color='k',
                             alpha=0.2)
        boxx = np.array([k*eps, (k+1)*eps, (k+1)*eps, k*eps, k*eps])
        boxy = np.array([l*eps, l*eps, (l+1)*eps, (l+1)*eps, l*eps])
        plt.plot(boxx, boxy, 'k-')

plt.plot(x, y, ',', color='tab:red')

# plt.plot(data[:, 0],
#        (np.mod(data[:, 1] + phase, 2*np.pi)-phase)/(np.pi), ',', color='tab:red')

# plt.axis('equal')
# plt.xlim([0,1.0])
# plt.ylim([-.6,.6])
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
# plt.ylim([0,2])
# plt.xlabel(r'$r$')
#plt.ylabel(r'$\vartheta / \pi$')
#plt.gca().set_yticks([0, .5, 1.0, 1.5, 2.0])
# plt.tight_layout()
# %%


# %%
