# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 18:12:14 2018

@author: chral
"""

#import matplotlib
#matplotlib.use('Agg')

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def doplot(z, z2, marker):
    plt.plot(z[:,0]*np.cos(z[:,1]), z[:,0]*np.sin(z[:,1]), marker)
    plt.plot(z2[:,0]*np.cos(z2[:,1]), z2[:,0]*np.sin(z2[:,1]), marker)
    plt.xlabel('R-R0')
    plt.ylabel('Z')
    plt.grid(True)
    
def doplot3(ax, z, z2, marker):
    ax.plot(z[:,0]*np.cos(z[:,2])*np.cos(z[:,1]), 
            z[:,0]*np.sin(z[:,2])*np.cos(z[:,1]),
            z[:,0]*np.sin(z[:,1]), marker)
    ax.plot(z2[:,0]*np.cos(z2[:,2])*np.cos(z2[:,1]), 
            z2[:,0]*np.sin(z2[:,2])*np.cos(z2[:,1]),
            z2[:,0]*np.sin(z2[:,1]), marker)
    plt.xlabel('X')
    plt.ylabel('Y')
    ax.set_zlabel('Z')
    ax.grid(False)
    
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

#prefix = '/home/calbert/run/NEO-ORB/RK_1em10_Euler1_32/'
prefix = '/home/calbert/run/NEO-ORB/RK_1em8_Euler1_32/'
#prefix = '/home/calbert/run/NEO-ORB/RK_1em10_Euler1_128/'
#prefix = 'L:/run/NEO-ORB/'
#prefix = '/home/calbert/mnt/marconi_scratch/NEO-ORB/RK_1em10_Verlet_32/'
#prefix = '/home/calbert/run/NEO-ORB/'

data = np.loadtxt(os.path.join(prefix,'orbit_kinds.out'))
difffilter = abs(data[:,6] - data[:,7]) > 1e-5
datadiff = data[difffilter,:]
ndiff = len(datadiff)
partdiff = datadiff[:,0].astype(int)

#%%
ipart = data[10,0].astype(int)
#ipart = 601
#ipart = partdiff[3]
#ipart = 1
num = '{:03d}'.format(ipart)
    
z = np.empty((1,3)); z[:] = np.nan
z2 = np.empty((1,3)); z2[:] = np.nan
try:
    z = np.loadtxt(prefix+'fort.{}{}'.format('10',num))
    if(len(z) == 0) : z = np.empty((1,3)); z[:] = np.nan
except:
    pass
try:
    z2 = np.loadtxt(prefix+'fort.{}{}'.format('11',num))
    if(len(z2) == 0) : z2 = np.empty((1,3)); z2[:] = np.nan
except:
    pass

plt.figure()
nsmooth = 1000
par_inv_sy = np.convolve(z2[1:,-1], np.ones((nsmooth,))/nsmooth, mode='valid')
par_inv_rk = np.convolve(z[1:,-1], np.ones((nsmooth,))/nsmooth, mode='valid')

nb_sy = z.shape[0]/par_inv_sy.shape[0]*np.arange(par_inv_sy.shape[0])
nb_rk = z.shape[0]/par_inv_rk.shape[0]*np.arange(par_inv_rk.shape[0])

plt.plot(nb_sy, par_inv_sy/par_inv_sy[0])
plt.plot(nb_rk, par_inv_rk/par_inv_rk[0])
#%%
plt.figure(figsize=(2.4,2.4))
doplot(z, z2, ',')
frame1 = plt.gca()
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])
plt.axis('equal')
plt.title('tip cut, ipart={:3d}'.format(ipart))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
doplot3(ax, z, z2, ',')
plt.title('tip cut, ipart={:3d}'.format(ipart))

#%%
plt.figure()
z = np.zeros([1,2]); z[:,0] = 0.5
z2 = np.zeros([1,2]); z2[:,0] = 0.5
try:
    z = np.loadtxt(prefix+'fort.{}{}'.format('20',num))
    if(len(z) == 0) : z = np.empty((1,3)); z[:] = np.nan
except:
    pass
try:
    z2 = np.loadtxt(prefix+'fort.{}{}'.format('21',num))
    if(len(z2) == 0) : z2 = np.empty((1,3)); z2[:] = np.nan
except:
    pass
doplot(z, z2, ',')
plt.title('period cut, ipart={:3d}'.format(ipart))

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#doplot3(ax, z, z2, ',')
#plt.title('period cut, ipart={:3d}'.format(ipart))
plt.show()
print('different classifications: {}'.format(ndiff))
print(sorted(partdiff))

