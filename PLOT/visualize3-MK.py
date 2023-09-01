#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 14:13:24 2019

@author: Christopher Albert, Majid Khan, Daniele Corrias
"""
# %%
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.stats import gaussian_kde


"""
Input file simple.in
* tcut: Cut time for Minkowski classification
* fast_class: Use new classifiers (2023), don't need tcut
* class_plot: exit after Minkowski classification to generate
"Brazilian flag" plot (Fig. 8 of Accelerated Methods paper)

Data columns in files

* start.dat
1) zstart(1,ipart)=r (normalized toroidal flux s)
2) zstart(2,ipart)=theta_vmec
3) zstart(3,ipart)=varphi_vmec
4) normalized velocity module z(4) = v / v_0:
zstart(4,ipart)=1.d0  Means all particles start with 3.5 Mev energies
5) starting pitch z(5)=v_\parallel / v:

* times_lost.dat
1) Particle index. Corresponds to line number in start.dat .
2) Time t_loss [s] when the particle is lost. [It can have three values, -1, +1, any other value]
  * If never lost or classified as regular, maximum tracing time trace_time is written.
  * If ignored due to contr_pp, which defines deep passing region as confined (we don consider them anymore), -1 is written.
3) Trapping parameter trap_par that is 1 for deeply trapped, 0 for tp boundary
   and negative for passing. Eq. (3.1) in Accelerated Methods paper.
   Whenever trap_par < contr_pp, particle is not traced and counted as confined.

* confined_fraction.dat
1) Time in s, according to number of recording timesteps ntimstep.
2) confpart_pass: Number of confined passing particles / total number of particles (confined fraction of passing)
3) confpart_trap: Number of confined trapped particles / total number of particles (confined fraction of trapped)
4) total number of particles

Col 2 + Col 3 gives Number of confined particles / total number of particles  (total confined fraction)
In the first line, the total fractions of trapped and passing particles are written (all confined in 1st step)
"""

#basedir = '/home/calbert/cobraruns/NEO-ORB/QA/2019-09-17_10k'

#basedir = '/home/calbert/net/cobra/run/SIMPLE/QA/2020-01-22_s0.3_withcut'
#basedir = '/home/calbert/net/cobra/run/SIMPLE/QH/2020-02-10_s0.3_withcut'
#basedir = '/home/calbert/net/cobra/run/SIMPLE/QA/2020-02-10_s0.6_withcut'
#basedir = '/home/calbert/net/cobra/run/NEO-ORB/QH/2019-09-17_1k'
basedir = '/afs/itp.tugraz.at/user/khanm/Downloads/ISHW_2019/data/QI'
#basedir = '/afs/itp.tugraz.at/user/khanm/Downloads/SIMPLE-develop-2023/build/'
#basedir = '/afs/itp.tugraz.at/user/khanm/Downloads/ISHW_2019/data/QA/10k'
# %%
data = np.loadtxt(os.path.join(basedir, 'times_lost.dat'))
tlost = np.abs(data[:, 1])
trap_par = data[:, 2]

npart = data.shape[0]    # the number of rows in the data array is stored in the npart, which is actually total particles.

data = np.loadtxt(os.path.join(basedir, 'start.dat'))
r0 = data[:npart, 0]
th0 = np.mod(data[:npart, 1], 2*np.pi)  #to ensure that values in th0 and ph0 are within the range of 0 to 2π.
ph0 = np.mod(data[:npart, 2], 2*np.pi)
p0 = data[:npart, 3]    # Normalized initial velocity of alphas, =1
lam0 = data[:npart, 4]  # Pitch angle cosine

plt.figure()
plt.plot(th0, ph0, '.')    # plots values of th0 on the x-axis and ph0 on the y-axis using dots ('.') as markers.
# %%

data = np.loadtxt(os.path.join(basedir, 'confined_fraction.dat'))
tconf = np.abs(data[1:, 0])        #  Starting from line 2 till the end in column 0
conffrac = data[1:, 1]+data[1:, 2] #  total confined fraction is sume of two columns (passing + trapp)
plt.figure()
plt.plot(np.log10(tconf), conffrac, color='tab:red')
plt.xlabel(r'Log(t)')   # r stands for raw string, backslashes within the string should be treated as literal characters rather than escape characters.

# %%
plt.figure()
plt.scatter(th0, lam0, c=np.log10(tlost), s=0.6)  # The color of each point is determined by np.log10(tlost), and the size of the points= 0.5.
plt.colorbar() #The color bar is added to the plot. Note that tlost=1 for confined particles and log(1)=0 for all such particles
               # for example the yellow color in figure represents confined. All other are lost particles having different colours.
# %%

nbins = 25
bins_lam = np.linspace(-0.9, 0.9, nbins-1)


bins_lam_center = np.linspace(-0.95, 0.95, nbins)  # Not used here in this code
bin_idx = np.digitize(lam0, bins_lam)

t1 = 0.0
t2 = 0.99
lossfrac = np.zeros(nbins)
npart_bins = np.zeros(nbins)
for k in range(npart):
    if t1 < np.abs(tlost[k]):   # I think this will always be satisfied bcz t1=0, and abs(tlost) >0
        npart_bins[bin_idx[k]] += 1.0
    if t1 < np.abs(tlost[k]) and np.abs(tlost[k]) < t2:  # Which means tlost[k] is between t1 and t2, means a lost particle
        lossfrac[bin_idx[k]] += 1.0

lossfrac = lossfrac/npart_bins


kde2 = gaussian_kde([np.log10(tlost[tlost < t2]), trap_par[tlost < t2]])
X, Y = np.mgrid[-5:0:100j, -2:1:200j]       # creates an array with 100 equally spaced values from -5 to 0 in the first dimension (X-axis).
positions = np.vstack([X.ravel(), Y.ravel()]) #that vertically stacks arrays. It concatenates the flattened X and Y arrays together to create
                                              #  a new array positions with shape (2, num_points), where num_points is the total number of points in the grid.
Z = np.reshape(kde2(positions).T, X.shape) #It takes the positions array as input and returns a 1D array of estimated densities.
                       #After executing these lines of code, Z will hold the estimated densities evaluated on the grid defined by X and Y
plt.figure()
plt.imshow(np.rot90(Z), cmap=plt.cm.Blues, extent=[-5, 0, -2, 1])    #extent controls the blue color
plt.plot(np.log10(tlost), trap_par, 'kx', markersize=0.5, alpha=0.6)  # k for black, x for marker, alpha sets the transparency of the markers to 50%.
#plt.plot(np.log10(tlost[wrong]), trap_par[wrong], 'rx', markersize=2, alpha=1)
plt.plot([-5, 0], [0, 0], 'k--', linewidth=.8)   #is used to create a dashed black line at y=0
plt.ylim([-2.1, 1.0])
plt.xlim(-5.0, 0.1)
plt.ylabel(r'trapping parameter')
plt.xlabel(r'$t\,/\,\mathrm{s}$')
plt.gca().set_xticklabels([r'$10^{-5}$', r'$10^{-4}$',
                            r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$1$'])
ax2 = plt.gca().twinx()     #his section creates a twin axes object (ax2) sharing the x-axis with the previous plot.
ax2.plot(np.log10(tconf), conffrac, color='tab:red')   #
ax2.set_ylim([0.5, 1.25])
ax2.set_ylabel('confined fraction                    ', color='tab:red')
ax2.set_yticks(np.linspace(0.5, 1, 6))
ax2.tick_params('y', colors='tab:red')
# plt.savefig('qi.png',dpi=150)
plt.show()
# %%
