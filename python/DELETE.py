#%% 
import matplotlib.pyplot as plt
import numpy as np
import common
from common import (r0,th0,ph0,pph0,timesteps,get_val,get_der,mu,qe,m,c,)
from field_correct_test import field
from plotale import plot_orbit, plot_mani, plot_cost_function
from cpp_sym_midpoint import z
from manifolds import surf_R

print(z.shape)
print(surf_R[:,0])

R0 = 1
z0 = 1


# Plotting
#fig = plt.figure(figsize=(8, 6))
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(X, Y, Z, cmap='plasma', edgecolor='none')

x = (R0 + z[0,:]*np.cos(z[1,:]))*np.cos(z[2,:])
y = (R0 + z[0,:]*np.cos(z[1,:]))*np.sin(z[2,:])
z = z0 + z[0,:]*np.sin(z[1,:])
#ax.set_zlim([-1,1])


fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter(surf_R[:,0], surf_R[:,1], surf_R[:,2], c='b', s=0.1)

sc = ax.scatter(x,y,z)

plt.show()