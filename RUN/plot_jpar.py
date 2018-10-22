# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 18:12:14 2018

@author: chral
"""

#import matplotlib
#matplotlib.use('Agg')

from numpy import *
from matplotlib.pyplot import *

data_rk = loadtxt('fort.100') # Runge Kutta on VMEC coordinates
data_can = loadtxt('fort.101') # Runge Kutta on canonicalized coordinates
data_sym = loadtxt('fort.102') # Symplectic Euler r, vpar

figure()
plot(data_can[1:,0], data_can[1:,1],'b,-')
plot(data_sym[1:,0], data_sym[1:,1],'r,-')
xlabel('time step')
ylabel('J_par (approx)')
legend(['Verlet','RK45'])
#plot(data_rk[1:,0], data_rk[1:,1],'g,')

