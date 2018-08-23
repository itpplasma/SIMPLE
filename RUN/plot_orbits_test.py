# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 18:12:14 2018

@author: chral
"""

from numpy import *
from matplotlib.pyplot import *

z = loadtxt('fort.4001')

figure()
plot(1+z[:,0]*cos(z[:,1]),1+z[:,0]*sin(z[:,1]),',')

figure()
plot(z[:,3],',-')
