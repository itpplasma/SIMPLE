#%% 
import matplotlib.pyplot as plt
import numpy as np
import common
from cpp_sym_midpoint import cpp_sym_mid
from common import (r0,th0,ph0,vpar0,timesteps,get_val,get_der,mu,qe,m,c,)
from field_correct_test import field
from plotale import plot_orbit, plot_mani, plot_cost_function
#from manifolds import surf_R

f = field() 

dt = 600
nt = 1000

# Initial Conditions
z = np.zeros([3, nt + 1])
z[:, 0] = [r0, th0, ph0]

f.evaluate(r0, th0, ph0)
v = np.zeros(3)
p = np.zeros(3)

v[0] = vpar0*f.co_hr
v[1] = vpar0*f.co_hth
v[2] = vpar0*f.co_hph
p[0] = v[0] 
p[1] = v[1] + qe/c*f.co_Ath
p[2] = v[2] + qe/c*f.co_Aph

P = np.zeros([3,nt+1])
P[:,0] = p

print(P)

# Symplectic midpoint
for kt in range(nt):

    z[:,kt+1] = cpp_sym_mid(z[:,kt])


plot_orbit(z)
plt.plot(z[0,:5]*np.cos(z[1,:5]), z[0,:5]*np.sin(z[1,:5]), "o")
plt.show()