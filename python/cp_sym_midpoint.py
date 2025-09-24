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

dt = 1
nt = 1000
print(dt*nt)
metric = lambda z: {
    "_ii":  np.array([1, z[0]**2, (1 + z[0]*np.cos(z[1]))**2]),
    "^ii":  np.array([1, 1/z[0]**2, 1/(1 + z[0]*np.cos(z[1]))**2]),
    "d_11": np.array([0,0,0]),
    "d_22": np.array([2*z[0], 0,0]),
    "d_33": np.array([2*(1+z[0]*np.cos(z[1]))*np.cos(z[1]), -2*(1+z[0]*np.cos(z[1]))*np.sin(z[1]), 0]),
}


def dLdx(v, dAth, dAph, g):
    ret = np.zeros(3)
    ret = (m/2*(g['d_11']*v[0]**2 + g['d_22']*v[1]**2 + g['d_33']*v[2]**2) + qe/c*(dAth*v[1] + dAph*v[2]))
    return ret


def F(x, xold, pold):
    global p
    ret = np.zeros(3)

    xmid = (x + xold)/2
    vmid = (x - xold)/dt

    f.evaluate(xmid[0], xmid[1], xmid[2])
    Amid = np.array([0, f.co_Ath, f.co_Aph])
    g = metric(xmid)
    pmid = pold + dt/2 * dLdx(vmid, f.co_dAth, f.co_dAph, g)
    p = pold + dt * dLdx(vmid, f.co_dAth, f.co_dAph, g)
    # d/dt p = d/dt dL/dq. = dL/dq
    ret = x - xold - dt/m*g['^ii']*(pmid - qe/c*Amid)
    return ret


# Initial Conditions
z = np.zeros([3, nt + 1])
z[:, 0] = [r0, th0, ph0]

f.evaluate(r0, th0, ph0)
g = metric(z[:,0])

vel = np.zeros([3, nt + 1])
vel[0,0] = np.sqrt(g['^ii'][0]*mu*2*f.B)
vel[1,0] = 0 #-np.sqrt(g['^22']*mu*2*f.B/3)
vel[2,0] = 0 #-np.sqrt(g['^33']*mu*2*f.B/3)
p = np.zeros(3)
p[0] = g['_ii'][0]*vel[0,0]
p[1] = g['_ii'][1]*vel[1,0] + qe/c * f.co_Ath
p[2] = g['_ii'][2]*vel[2,0] + qe/c * f.co_Aph

vper = np.zeros(nt)
vpar = np.zeros(nt)
v = np.zeros(nt)

'''
v = np.zeros(3)
p = np.zeros(3)
vel = np.zeros([3, nt + 1])

p[0] = vel[0,0] 
p[1] = vel[1,0] + qe/c*f.co_Ath
p[2] = vel[2,0] + qe/c*f.co_Aph
'''

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


q_cp = z
v_cp = vel


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



plot_mani(q_cp)
plt.show()

fig, ax = plt.subplots()
plot_orbit(q_cp, ax = ax)
plt.plot(q_cp[0,:100]*np.cos(q_cp[1,:100]), q_cp[0,:100]*np.sin(q_cp[1,:100]), "o")
plt.show()
'''






'''
dt, nt = timesteps(steps_per_bounce=8, nbounce=100)
print(dt)
print(nt)

f = field()
dt = 1
nt = 10000

z = np.zeros([3, nt + 1])
z[:, 0] = [r0, th0, ph0]

metric = lambda z: {
    "_11": 1,
    "_22": z[0]**2,
    "_33": (1 + z[0]*np.cos(z[1]))**2,
    "^11": 1,
    "^22": 1/z[0]**2,







    "^33": 1/(1 + z[0]*np.cos(z[1]))**2,
    "d_11": [0,0,0],
    "d_22": [2*z[0], 0,0],
    "d_33": [2*(1+z[0]*np.cos(z[1]))*np.cos(z[1]), -2*(1+z[0]*np.cos(z[1]))*np.sin(z[1]), 0],
}

def implicit_p(p, pold, Ath, Aph, dAth, dAph, g):
    ret = np.zeros(3)

    ctrv = np.zeros(3)
    ctrv[0] = 1/m *g['^11'] *(p[0])
    ctrv[1] = 1/m *g['^22'] *(p[1] - qe/c*Ath)
    ctrv[2] = 1/m *g['^33'] *(p[2] - qe/c*Aph)

    ret[0] = p[0] - pold[0] - dt*(qe/c*(ctrv[1]*dAth[0] + ctrv[2]*dAph[0]) + m/2*(g['d_11'][0]*ctrv[0]**2 + g['d_22'][0]*ctrv[1]**2 + g['d_33'][0]*ctrv[2]**2))
    ret[1] = p[1] - pold[1] - dt*(qe/c*(ctrv[1]*dAth[1] + ctrv[2]*dAph[1]) + m/2*(g['d_11'][1]*ctrv[0]**2 + g['d_22'][1]*ctrv[1]**2 + g['d_33'][1]*ctrv[2]**2))
    ret[2] = p[2] - pold[2] - dt*(qe/c*(ctrv[1]*dAth[2] + ctrv[2]*dAph[2]) + m/2*(g['d_11'][2]*ctrv[0]**2 + g['d_22'][2]*ctrv[1]**2 + g['d_33'][2]*ctrv[2]**2))

    return ret

#Initial Conditions 
f.evaluate(r0, th0, ph0)
g = metric(z[:,0])
ctrv = np.zeros(3)
ctrv[0] = np.sqrt(g['^11']*mu*2*f.B)
ctrv[1] = 0 #-np.sqrt(g['^22']*mu*2*f.B/3)
ctrv[2] = 0 #-np.sqrt(g['^33']*mu*2*f.B/3)
p = np.zeros(3)
p[0] = g['_11']*ctrv[0]
p[1] = g['_22']*ctrv[1] + qe/c * f.co_Ath
p[2] = g['_33']*ctrv[2] + qe/c * f.co_Aph
pold = p

from time import time
tic = time()
for kt in range(nt):

    sol = root(implicit_p, p, method='hybr',tol=1e-12,args=(pold, f.co_Ath, f.co_Aph, f.co_dAth, f.co_dAph, g))
    p = sol.x
    pold = p 

    z[0,kt+1] = z[0,kt] + dt/m * g['^11']*(p[0])
    z[1,kt+1] = z[1,kt] + dt/m * g['^22']*(p[1] - qe/c*f.co_Ath)
    z[2,kt+1] = z[2,kt] + dt/m * g['^33']*(p[2] - qe/c*f.co_Aph)

    f.evaluate(z[0, kt+1], z[1, kt+1], z[2, kt+1])
    g = metric(z[:,kt+1])

    p = np.zeros(3)
    p[0] = g['_11']*m*(z[0,kt+1]-z[0,kt])/dt
    p[1] = g['_22']*m*(z[1,kt+1]-z[1,kt])/dt + qe/c*f.co_Ath
    p[2] = g['_33']*m*(z[2,kt+1]-z[2,kt])/dt + qe/c*f.co_Aph
'''



# %%
