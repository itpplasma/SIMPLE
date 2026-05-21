#%%
'''
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
'''

import numpy as np
import matplotlib.pyplot as plt

# ----- define plane Q -----
x = np.linspace(-2, 2, 20)
y = np.linspace(-2, 2, 20)
X, Y = np.meshgrid(x, y)
Z = 0.5 * X + 0.2 * Y   # plane Q

# ----- normal vector (perpendicular axis) -----
n = np.array([-0.5, -0.2, 1.0])
n = n / np.linalg.norm(n)

# axis line
t = np.linspace(-2, 2, 100)
Lx = t * n[0]
Ly = t * n[1]
Lz = t * n[2]

# ----- plot -----
fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111, projection="3d")

# plane Q
ax.plot_surface(X, Y, Z, alpha=0.4)

# perpendicular axis
ax.plot(Lx, Ly, Lz, linewidth=3)
ax.scatter(0, 0, 0, s=50)  # origin

# labels
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Plane Q and Perpendicular Axis", fontsize=14)

# equal scaling
ax.set_box_aspect((1, 1, 1))

# clean look
ax.grid(False)

plt.show()
