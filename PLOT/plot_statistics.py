#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 14:51:30 2019

@author: calbert
"""
#%%
import os
import numpy as np
import matplotlib.pyplot as plt

inpath = ['/home/calbert/mnt/draco/run/NEO-ORB/2019-01-08_alphas/',
          '/home/calbert/mnt/draco/run/NEO-ORB/2019-01-08_alphas_sympl/',
          '/home/calbert/mnt/draco/run/NEO-ORB/2019-01-08_alphas_sympl2/',
          '/home/calbert/mnt/draco/run/NEO-ORB/2019-01-08_alphas_sympl3/',
          '/home/calbert/mnt/draco/run/NEO-ORB/2019-01-08_alphas_sympl4/',
          '/home/calbert/run/NEO-ORB/2019-01-08_alphas_sympl/']

def plotdata(infile):
    data = np.loadtxt(infile)
    t = data[:,0]
    fracpass = data[:,1]
    fractrap = data[:,2]
    
    plt.figure(1)
    plt.semilogx(t+1e-3, fracpass+fractrap)
    plt.xlabel('t')
    plt.ylabel('% confined (total)')
    plt.grid(True, which="both")
    plt.tight_layout()
    
    plt.figure(2)
    plt.semilogx(t+1e-3, fractrap)
    plt.xlabel('t')
    plt.ylabel('% confined (trapped, relative)')
    plt.grid(True, which="both")
    plt.tight_layout()
    
#cases = ['current', 'RKV', 'RK', '0016', '0032', '0064', '0128']
#cases = ['reference', '64', '32', '16', r'RK $10^{-6}$', r'RK $10^{-5}$', r'RK $10^{-4}$']
cases = ['reference', '64', '32', 'test']

plotdata(os.path.join(inpath[1], 'confined_fraction_RK.dat'))
#plotdata(os.path.join(inpath[1], 'confined_fraction_1024.dat'))
plotdata(os.path.join(inpath[1], 'confined_fraction_0064.dat'))
plotdata(os.path.join(inpath[1], 'confined_fraction_0032.dat'))
#plotdata(os.path.join(inpath[1], 'confined_fraction_0016.dat'))
#plotdata(os.path.join(inpath[1], 'confined_fraction_RK1d-6.dat'))
#plotdata(os.path.join(inpath[1], 'confined_fraction_RK1d-5.dat'))
#plotdata(os.path.join(inpath[1], 'confined_fraction_RK1d-4.dat'))
#plotdata(os.path.join(inpath[1], 'confined_fraction_RK1d-3.dat'))
plotdata(os.path.join(inpath[1], 'confined_fraction.dat'))
#for run in cases[1:3]:
#  plotdata(os.path.join(inpath[5], 'confined_fraction_{}.dat'.format(run)))

plt.figure(1)
plt.legend(cases)
plt.figure(2)
plt.legend(cases)

plt.show()

#%%
