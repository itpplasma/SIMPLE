# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 18:12:14 2018

@author: chral
"""

#import matplotlib
#matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt


def doplot(z, marker):
    plt.figure(1)
    plt.plot(z[:,2], z[:,0], marker)
    
#    plt.figure(2)
#    plt.plot(z[:,2], z[:,0], marker)
#    plt.xlim([.9985,1.001])
#    plt.ylim([.4875,.491])
#%%

#z = np.loadtxt('/home/calbert/run/NEO-ORB/poiplot_euler16.dat')
#doplot(z, 'g,')
#
#z = np.loadtxt('/home/calbert/run/NEO-ORB/poiplot_verlet16.dat')
#doplot(z, 'b,')
#
#z = np.loadtxt('/home/calbert/run/NEO-ORB/poiplot_verlet8.dat')
#doplot(z, 'c,')
#
#z = np.loadtxt('/home/calbert/run/NEO-ORB/poiplot_rk16.dat')
#doplot(z, 'k,')
#    
#z = np.loadtxt('/home/calbert/run/NEO-ORB/poiplot.dat')
#doplot(z, 'r,')

#plt.legend(['Euler16 w Taylor', 'Euler16', 'Verlet16', 'Verlet8 old', 'RK16'])
    
z = np.loadtxt('L:/run/NEO-ORB/fort.10003')
doplot(z, 'b,')

z = np.loadtxt('L:/run/NEO-ORB/fort.11003')
doplot(z, 'r,')
