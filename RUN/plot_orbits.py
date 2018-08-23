# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 18:12:14 2018

@author: chral
"""

from numpy import *
from matplotlib.pyplot import *

data0 = loadtxt('fort.3001')
data1 = loadtxt('fort.3003')
data2 = loadtxt('fort.3004')

data2 = data2[1:100,:]

figure()
plot(data0[:len(data2),1])
plot(data1[:len(data2),1])
plot(data2[:,1])

figure()
plot(data0[:len(data2),2])
plot(data1[:len(data2),2])
plot(data2[:,2])

figure()
plot(data0[:len(data2),3])
plot(data1[:len(data2),3])
plot(data2[:len(data2),3])

#figure()
#plot(data0[:,0],data0[:,1])
#plot(data1[:,0],data1[:,1])
#plot(data2[:,0],data2[:,1])
figure()
plot(data2[:,4])
figure()
plot(data2[:,5])
figure()
plot(data0[:,2],data0[:,3])
plot(data1[:,2],data1[:,3])
plot(data2[:,2],data2[:,3])
plot()

#%%

figure()
plot(data1[:len(data2),0],data1[:len(data2),1],'b,')
plot(data2[:len(data2),0],data2[:len(data2),1],'r,')
xlabel('t/tau0')
ylabel('r(t)')

figure()
plot(data1[:len(data2),-2],data1[:len(data2),-1],'b,')
plot(data2[:len(data2),-2],data2[:len(data2),-1],'r,')
xlabel('th_c(t)')
ylabel('ph_c(t)')
