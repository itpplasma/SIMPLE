# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import (r0,th0,ph0,vpar0,timesteps,get_val,get_der,mu,qe,m,c,)
from field_correct_test import field
from plotale import plot_orbit, plot_mani, plot_cost_function

f = field()

dt = 80
nt = 1000

metric = lambda z: {
    "_ii":  np.array([1, z[0]**2, (1 + z[0]*np.cos(z[1]))**2]),
    "^ii":  np.array([1, 1/z[0]**2, 1/(1 + z[0]*np.cos(z[1]))**2]),
    "d_11": np.array([0,0,0]),
    "d_22": np.array([2*z[0], 0,0]),
    "d_33": np.array([2*(1+z[0]*np.cos(z[1]))*np.cos(z[1]), -2*(1+z[0]*np.cos(z[1]))*np.sin(z[1]), 0]),
}

def dLdx(v, dAth, dAph, dB, g):
    ret = np.zeros(3)
    ret = (m/2*(g['d_11']*v[0]**2 + g['d_22']*v[1]**2 + g['d_33']*v[2]**2) + qe/c*(dAth*v[1] + dAph*v[2]) - mu*dB)
    return ret


def F(x, xold, pold):
    global p
    ret = np.zeros(3)

    xmid = (x + xold)/2
    vmid = (x - xold)/dt

    f.evaluate(xmid[0], xmid[1], xmid[2])
    Amid = np.array([0, f.co_Ath, f.co_Aph])
    g = metric(xmid)
    pmid = pold + dt/2 * dLdx(vmid, f.co_dAth, f.co_dAph, f.dB, g)
    p = pold + dt * dLdx(vmid, f.co_dAth, f.co_dAph, f.dB, g)
    # d/dt p = d/dt dL/dq. = dL/dq
    ret = x - xold - dt/m*g['^ii']*(pmid - qe/c*Amid)
    return ret


# Initial Conditions
z = np.zeros([3, nt + 1])
z[:, 0] = [r0, th0, ph0]

f.evaluate(r0, th0, ph0)
g = metric(z[:,0])
v = np.zeros(3)
p = np.zeros(3)
vel = np.zeros([3, nt + 1])
vel[:,0] = [vpar0*f.co_hr, vpar0*f.co_hth, vpar0*f.co_hph]

p[0] = vel[0,0] 
p[1] = vel[1,0] + qe/c*f.co_Ath
p[2] = vel[2,0] + qe/c*f.co_Aph

vper = np.zeros(nt)
vpar = np.zeros(nt)
v = np.zeros(nt)

from time import time
tic = time()
for kt in range(nt):
    pold = p

    sol = root(F, z[:,kt], method='hybr',tol=1e-12,args=(z[:,kt], pold))
    z[:,kt+1] = sol.x

    f.evaluate(z[0,kt+1], z[1,kt+1], z[2,kt+1])
    g = metric(z[:,kt+1])
    vel[:,kt+1] = 1/m * g['^ii']*(p - qe/c*np.array([0, f.co_Ath, f.co_Aph]))


    vper[kt] = np.sqrt(g['^ii'][0]*mu*2*f.B) 
    vpar[kt] = np.sqrt(np.sum(vel[:,kt+1]**2)) - vper[kt]
    v[kt] = np.sum(vel[:,kt+1])






'''
t = np.linspace(0, nt, nt)

# k = nt
k = 100

plt.plot(t[:k], vpar[:k], linewidth=0.5, c='r', label='vpar')

plt.plot(t[:k], vper[:k], linewidth=0.5, c='g', label='vper')

plt.plot(t[:k], v[:k], linewidth=0.5, c='b', label='v')
plt.show()


plt.hist(v, bins=30, color='blue', edgecolor='black', alpha=0.7)
plt.hist(vpar, bins=30, color='red', edgecolor='black', alpha=0.7)
plt.hist(vper, bins=30, color='green', edgecolor='black', alpha=0.7)

plt.xlabel("Value")
plt.ylabel("Frequency")
plt.title("Histogram of velocity")
plt.grid(True, linestyle="--", alpha=0.5)
plt.show()

plt.plot
plot_mani(z)
#plt.plot(z[0,:]*np.cos(z[1,:]), z[0,:]*np.sin(z[1,:]), ".")
plt.show()

fig, ax = plt.subplots()
plot_orbit(z, ax = ax)
plt.plot(z[0,:5]*np.cos(z[1,:5]), z[0,:5]*np.sin(z[1,:5]), "o")
plt.show()
'''