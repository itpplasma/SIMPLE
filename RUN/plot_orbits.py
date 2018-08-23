# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 18:12:14 2018

@author: chral
"""

from numpy import *
from matplotlib.pyplot import *

data_rk = loadtxt('fort.3001') # Runge Kutta
data_can = loadtxt('fort.3003') # Runge Kutta on canonicalized coordinates
data_sym = loadtxt('fort.3004') # Symplectic Euler r, vpar

nrange = range(1,1001)
data_rk = data_rk[nrange,:]
data_can = data_can[nrange,:]
data_sym = data_sym[nrange,:]

figure()
plot(data_can[:,0], data_can[:,1],'.-')
plot(data_sym[:,0], data_sym[:,1],'.-')
plot(data_rk[:,0], data_rk[:,1],'.-')

figure()
plot(data_can[:,0], data_can[:,-2],'.-')
plot(data_sym[:,0], data_sym[:,-2],'.-')

figure()
plot(data_can[:,0], data_can[:,-1],'.-')
plot(data_sym[:,0], data_sym[:,-1],'.-')
