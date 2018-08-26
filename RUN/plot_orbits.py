# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 18:12:14 2018

@author: chral
"""

from numpy import *

import matplotlib
#matplotlib.use('Agg')
from matplotlib.pyplot import *

data_rk = loadtxt('fort.3001') # Runge Kutta
data_can = loadtxt('fort.3003') # Runge Kutta on canonicalized coordinates
data_sym = loadtxt('fort.3004') # Symplectic Euler r, vpar

nrange = range(1,len(data_sym))
#nrange = range(1,1000)

data_rk = data_rk[nrange,:]
data_can = data_can[nrange,:]
data_sym = data_sym[nrange,:]

figure()
plot(data_can[:,0], data_can[:,1],'b-')
plot(data_sym[:,0], data_sym[:,1],'r-')
plot(data_rk[:,0], data_rk[:,1],'g-')
savefig('orbit_stell.png')

figure()
plot(data_can[:,0], data_can[:,-2],'b-')
plot(data_sym[:,0], data_sym[:,-2],'r-')
savefig('orbit_stell_thcan.png')

figure()
plot(data_can[:,0], data_can[:,-1],'b-')
plot(data_sym[:,0], data_sym[:,-1],'r-')
savefig('orbit_stell_phcan.png')

