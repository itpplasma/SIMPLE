# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 18:12:14 2018

@author: chral
"""

#import matplotlib
#matplotlib.use('Agg')

from numpy import *
from matplotlib.pyplot import *

#data_rk = loadtxt('fort.100') # Runge Kutta on VMEC coordinates
data_axi = loadtxt('fort.103') # Runge Kutta on canonicalized coordinates
data_sym = loadtxt('fort.102') # Symplectic Euler r, vpar


figure()
plot(data_axi[1:,0], data_axi[1:,1],'b,-')
plot(data_sym[1:,0], data_sym[1:,1],'r,-')
xlabel('time step')
ylabel('J_par (approx)')
legend(['RK45','symplectic'])
#plot(data_rk[1:,0], data_rk[1:,1],'g,')

nsmooth = 1000
J_axi_sm = convolve(data_axi[1:,1], ones((nsmooth,))/nsmooth, mode='valid')
J_sym_sm = convolve(data_sym[1:,1], ones((nsmooth,))/nsmooth, mode='valid')
figure()
plot(arange(len(J_axi_sm))*nsmooth, J_axi_sm,'b,-')
plot(arange(len(J_sym_sm))*nsmooth, J_sym_sm,'r,-')
#ylim([0,500])
xlabel('time step')
ylabel(r'$J_{\parallel}$')
legend(['RK45','symplectic'])
