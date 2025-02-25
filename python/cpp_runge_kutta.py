# %%
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import (f,r0,th0,ph0,pph0,timesteps,get_val,get_der,mu,qe,m,c,)
from plotting import plot_orbit, plot_cost_function

B0 = 1.0    # magnetic field modulus normalization
iota0 = 1.0 # constant part of rotational transform
a = 0.5     # (equivalent) minor radius
R0 = 1.0    # (equivalent) major radius

z0_cpp = np.array([r0, th0, ph0])

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

def field_Ath(r, th, ph):
    Ath = B0*(r**2/2 - r**3/(3*R0)*np.cos(th))
    dAth = np.array((B0*(r - r**2/R0*np.cos(th)), B0*r**3*np.sin(th)/(3*R0), 0))
    return Ath, dAth

def field_Aph(r, th, ph):
    Aph = -B0*iota0*(r**2/2-r**4/(4*a**2))
    dAph = np.array((-B0*iota0*(r-r**3/a**2), 0, 0))
    return Aph, dAph


def field_dB(r, th, ph):
    dB = np.array((iota0*(1-3*r**2/a**2)/(R0+r*np.cos(th)) - iota0*(r-r**3/a**2)/(R0+r*np.cos(th))**2 * np.cos(th) - np.cos(th)/R0, 
                    iota0*(r-r**3/a**2)/(R0+r*np.cos(th))**2 *r*np.sin(th)  + r*np.sin(th), 0))
    return dB


#Initial Conditions 
Ath = np.empty(1)
Aph = np.empty(1)
dAth = np.empty(3)
dAph = np.empty(3)
dB = np.empty(3)
Ath, dAth = field_Ath(r0, th0, ph0)
Aph, dAph = field_Aph(r0, th0, ph0)
dB = field_dB(r0, th0, ph0)

g = metric(z0_cpp[:])
p = np.zeros(3)
p[0] = 0
p[1] = qe/(m*c) * Ath
p[2] = qe/(m*c) * Aph
pold = p

#Initial Conditions Runge-Kutta
z0_cpp = np.array([r0, th0, ph0, p[0], p[1], p[2]])
print(z0_cpp)

def momenta_cpp(t, x):
    Ath = np.empty(1)
    Aph = np.empty(1)
    dAth = np.empty(3)
    dAph = np.empty(3)
    dB = np.empty(3)
    Ath, dAth = field_Ath(x[0], x[1], x[2])
    Aph, dAph = field_Aph(x[0], x[1], x[2])
    dB = field_dB(x[0], x[1], x[2])
    g = metric(x[:3])

    ctrv = np.empty(3)
    ctrv[0] = 1/m *g['^11'] *(x[3])
    ctrv[1] = 1/m *g['^22'] *(x[4] - qe/c*Ath)
    ctrv[2] = 1/m *g['^33'] *(x[5] - qe/c*Aph)

    xdot = np.empty(6)
    xdot[0] = 1/m * g['^11']*x[3]
    xdot[1] = 1/m * g['^22']*(x[4] - qe/c*Ath)
    xdot[2] = 1/m * g['^33']*(x[5] - qe/c*Aph)
    xdot[3] = qe/c*(ctrv[1]*dAth[0] + ctrv[2]*dAph[0]) - mu*dB[0] + m/2*(g['d_11'][0]*ctrv[0]**2 + g['d_22'][0]*ctrv[1]**2 + g['d_33'][0]*ctrv[2]**2)
    xdot[4] = qe/c*(ctrv[1]*dAth[1] + ctrv[2]*dAph[1]) - mu*dB[1] + m/2*(g['d_11'][1]*ctrv[0]**2 + g['d_22'][1]*ctrv[1]**2 + g['d_33'][1]*ctrv[2]**2)
    xdot[5] = qe/c*(ctrv[1]*dAth[2] + ctrv[2]*dAph[2]) - mu*dB[2] + m/2*(g['d_11'][2]*ctrv[0]**2 + g['d_22'][2]*ctrv[1]**2 + g['d_33'][2]*ctrv[2]**2)

    return xdot

sol_cpp = solve_ivp(
    momenta_cpp, [0, 7000], z0_cpp, method="RK45", rtol=1e-9, atol=1e-9,
)


nt = len(sol_cpp.t)
print(nt)  # This will print the number of time steps the solver used
print(7000/nt)
z = np.zeros([3, nt])
z[0,:] = sol_cpp.y[0]
z[1,:] = sol_cpp.y[1]
z[2,:] = sol_cpp.y[2]


plot_orbit(z)
plt.plot(z[0,:5]*np.cos(z[1,:5]), z[0,:5]*np.sin(z[1,:5]), "o")
plt.show()



# %%
