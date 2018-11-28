# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 18:12:14 2018

@author: chral
"""

#import matplotlib
#matplotlib.use('Agg')

from numpy import *
from matplotlib.pyplot import *

parent = '/home/calbert/run/NEO-ORB/'
#parent = '/home/calbert/code/NEO-ORB/RUN'

#data_rk = loadtxt('orbit_vmec.out') # Runge Kutta
#data_can = loadtxt('orbit_can.out') # Runge Kutta on canonicalized coordinates
data_axi = loadtxt(parent+'/orbit_axis.out') # Runge Kutta canonical with axis coorection
data_sym = loadtxt(parent+'/orbit_sympl.out') # Symplectic Euler r, vpar

nrange = range(1,len(data_sym))
#nrange = range(1,146)

#data_rk = data_rk[nrange,:]
#data_can = data_can[nrange,:]
#data_sym = data_sym[nrange,:]


figure()
plot(data_axi[2:,0], data_axi[2:,1],'b,')
plot(data_sym[2:,0], data_sym[2:,1],'r,')
#plot(data_can[:,0], data_can[:,1]+.02,'b,')
#plot(data_rk[:,0], data_rk[:,1]+.03,'k,')
legend(['can','sympl'])
#plot(data_rk[:,0], data_rk[:,1],'g,')
#savefig('orbit_stell.png')

#figure()
#plot(data_can[:,0], data_can[:,-2],'b,')
#plot(data_sym[:,0], data_sym[:,-2],'r,')
#plot(data_axi[:,0], data_axi[:,-2],'g,')
#legend(['can','sympl','axi'])
#savefig('orbit_stell_thcan.png')
#
#figure()
#plot(data_can[:,0], data_can[:,-1],'b,')
#plot(data_sym[:,0], data_sym[:,-1],'r,')
#plot(data_axi[:,0], data_axi[:,-1],'g,')
#legend(['can','sympl','axi'])
#savefig('orbit_stell_phcan.png')

#%%
nsmooth = 100
H_sm = convolve(data_sym[:,-4], ones((nsmooth,))/nsmooth, mode='valid')
#
figure()
plot(H_sm,'b,')

#%%

figure()
plot(np.mod(data_axi[:,-1], 2*np.pi), data_axi[:,1],'b,')
plot(np.mod(data_sym[:,-1], 2*np.pi), data_sym[:,1],'r,')
xlabel('phi')
ylabel('r')

figure()
plot(np.mod(data_axi[:,-1], 2*np.pi), np.mod(data_axi[:,-2], 2*np.pi),'b,')
plot(np.mod(data_sym[:,-1], 2*np.pi), np.mod(data_sym[:,-2], 2*np.pi),'r,')
xlabel('phi')
ylabel('theta')
