from numba import njit
from field import afield_simple, bfield_simple, bmod_simple
from scipy.integrate import solve_ivp
from scipy.optimize import root
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

taub = 7800 # estimated bounce time

def timesteps(steps_per_bounce, nbounce):
    """ compute step size and number of timesteps """
    return taub/steps_per_bounce, nbounce*steps_per_bounce


qe = 1.0; m = 1.0; c = 1.0 # particle charge, mass and speed of light
dt, nt = timesteps(steps_per_bounce = 8, nbounce = 100)

dt = 0.104
nt = 5000
z0_cpp = np.array([1.04, 2.04, 2.94, 0.1, 0.1, 0.1])  # Guessed guiding center
p = np.zeros(3)

# Set initial conditions for guiding center CPP
B = np.empty(3)
Bmod = np.empty(1)
dBmod = np.empty(3)
A = np.empty(3)
dA3 = np.empty(3)
bfield_simple(z0_cpp[0], z0_cpp[1], z0_cpp[2], B)
bmod_simple(z0_cpp[0], z0_cpp[1], z0_cpp[2], Bmod, dBmod)
afield_simple(z0_cpp[0], z0_cpp[1], z0_cpp[2], A, dA3)

vpar = np.dot(z0_cpp[3:], B/Bmod)
vperp = np.sqrt(np.sum(z0_cpp[3:]**2) - vpar**2)
mu = 0.5 * vperp**2 / Bmod
mu0 = mu[0]
print(mu0)

z = np.zeros([3,nt+1])
z[:,0] = z0_cpp[0:3]
vel = np.zeros([3,nt+1])
vel[:,0] = vpar*B/Bmod
v = vel[:,0]
vold = v


def implicit_v(v, vold, B, dBmod, mu0):
    ret = np.zeros(3)

    ret[0] = v[0] - vold[0] - dt*(v[1]*B[2] - v[2]*B[1] - mu0*dBmod[0]) 
    ret[1] = v[1] - vold[1] - dt*(v[2]*B[0] - v[0]*B[2] - mu0*dBmod[1])
    ret[2] = v[2] - vold[2] - dt*(v[0]*B[1] - v[1]*B[0] - mu0*dBmod[2])

    return ret


for kt in range(nt):

    B = np.empty(3)
    Bmod = np.empty(1)
    dBmod = np.empty(3)

    bmod_simple(z[0,kt], z[1,kt], z[2,kt], Bmod, dBmod)
    bfield_simple(z[0,kt], z[1,kt], z[2,kt], B)
    

    sol = root(implicit_v, v, method='hybr',tol=1e-12,args=(vold, B, dBmod, mu0))
    vold = v
    v = sol.x


    z[0,kt+1] = z[0,kt] + dt*v[0]
    z[1,kt+1] = z[1,kt] + dt*v[1]
    z[2,kt+1] = z[2,kt] + dt*v[2]



plt.figure()
plt.plot(z[0,:], z[1,:])
plt.plot(z0_cpp[0], z0_cpp[1], "ro")
plt.show()


# %%
