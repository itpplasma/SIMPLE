# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 18:12:14 2018

@author: chral
"""
import matplotlib
#matplotlib.use('Agg')

from numpy import *
from matplotlib.pyplot import *

z = loadtxt('fort.4001')

figure()
plot(1+z[:,0]*cos(z[:,1]),1+z[:,0]*sin(z[:,1]),',')
savefig('orbit.png')

figure()
plot(z[:,3],',-')
savefig('ham.png')