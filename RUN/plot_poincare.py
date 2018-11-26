# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 18:12:14 2018

@author: chral
"""

#import matplotlib
#matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt

z = np.loadtxt('/home/calbert/run/NEO-ORB/poiplot.dat')

plt.figure()
plt.plot(z[:,2], z[:,0], ',')

plt.figure()
plt.plot(z[:,2], z[:,0], '.')
plt.xlim([.9005,.903])
plt.ylim([.52793,.52796])

#%%