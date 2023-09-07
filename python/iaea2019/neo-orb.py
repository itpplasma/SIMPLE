#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 12:35:03 2019

@author: calbert
"""
# %%
import time
import numpy as np
import matplotlib.pyplot as plt
from pysimple import simple, params, orbit_symplectic, vmec_to_can,\
     can_to_vmec, vmec_to_cyl

from mpl_toolkits.mplot3d import Axes3D

s0 = 0.6
th0 = 0.0
ph0 = 0.25*np.pi
lam = -0.4  # trapped regular
#lam = 0.8   # passing regular

pi = 3.14159265358979
c = 2.9979e10
e_charge = 4.8032e-10
e_mass = 9.1094e-28
p_mass = 1.6726e-24
ev = 1.6022e-12

bmod_ref = 5.0e4
Z_charge = 2
m_mass = 4
E_kin = 3.5e6
trace_time = 1.0e-2

v0 = np.sqrt(2.0*E_kin*ev/(m_mass*p_mass))    # alpha particle velocity, cm/s
rlarm = v0*m_mass*p_mass*c/(Z_charge*e_charge*bmod_ref)  # reference gyroradius
tau = trace_time*v0

# %%
tracy = params.Tracer()

simple.init_field(tracy, "wout.nc", 5, 5, 3, 1)
simple.init_params(tracy, Z_charge, m_mass, E_kin, 256, 1, 1e-13)

# %%
s = np.array(0.7)
theta = np.array(0.0)
varphi = np.array(0.0)
A_phi = np.array(0.0)
A_theta = np.array(0.0)
dA_phi_ds = np.array(0.0)
dA_theta_ds = np.array(0.0)
aiota = np.array(0.0)
alam = np.array(0.0)
dR_ds = np.array(0.0)
dR_dt = np.array(0.0)
dR_dp = np.array(0.0)
dZ_ds = np.array(0.0)
dZ_dt = np.array(0.0)
dZ_dp = np.array(0.0)
dl_ds = np.array(0.0)
dl_dt = np.array(0.0)
dl_dp = np.array(0.0)
Bctrvr_vartheta = np.array(0.0)
Bctrvr_varphi = np.array(0.0)
Bcovar_r = np.array(0.0)
Bcovar_vartheta = np.array(0.0)
Bcovar_varphi = np.array(0.0)
sqg = np.array(0.0)
# %%
theta_vmec, varphi_vmec = can_to_vmec(s, theta, varphi)
R, Z = vmec_to_cyl(s, theta_vmec, varphi_vmec)
print("Initial canonical coordinates (s,th_c,ph_c) = ", s, theta, varphi)
print("Initial VMEC coordinates (s,th,ph) = ", s, theta_vmec, varphi_vmec)
print("Initial cylindrical coordinates (R,PH,Z) = ", R, varphi_vmec, Z)

# %%
npoiper2 = 128                       # interation points per field period
rbig = new_vmec_stuff.rmajor*1.0e2  # major radius in cm
dtaumax = 2.0*pi*rbig/npoiper2      # maximum time step for integrator
dtau = dtaumax                   # time step for output
neo_orb.init_params(Z_charge, m_mass, E_kin, dtau, dtaumax, 1e-8)

# %%
ntimstep = 10000
neo_orb.firstrun = 1  # to re-initialize symplectic integrator in case

z = np.array([s0, th0, ph0, 1.0, lam])
neo_orb.init_integrator(z)

zs = np.empty([4, ntimstep+1])
zs[:, 0] = z[[0, 1, 2, 4]]

t = time.time()
ierr = 0  # doesn't work yet with pass-by-reference
for kt in range(ntimstep):
    neo_orb.timestep_sympl_z(z, ierr)
    zs[:, kt+1] = z[[0, 1, 2, 4]]
print(z)
print('time elapsed: {} s'.format(time.time() - t))
kmax = ntimstep


# %%
neo_orb.firstrun = 1  # to re-initialize symplectic integrator in case
z = np.array([s0, th0, ph0, 1.0, lam])
neo_orb.init_integrator(z)
cut_detector.init(z)
ierr = 0
var_cut = np.zeros(6)
ncut = 1000
# variables to evaluate at tip: z(1..5), par_inv
var_tips = np.zeros([6, ncut])
#
t = time.time()
for k in np.arange(ncut):
    cut_type = 0
    cut_detector.trace_to_cut(z, var_cut, cut_type, ierr)
    # print(z)
    var_tips[:, k] = var_cut
print('time elapsed: {} s'.format(time.time() - t))

#%%

z_vmec = np.zeros_like(zs)
tip_vmec = np.zeros_like(var_tips)

for k in np.arange(zs.shape[1]):
    s[0] = zs[0, k]
    theta[0] = zs[1, k]
    varphi[0] = zs[2, k]
    if s[0] <= 0.0 or s[0] >= 1.0:
        print(s[0])
        break
    libneo_orb._lib.can_to_vmec_(s, theta, varphi, theta_vmec, varphi_vmec)
    z_vmec[0, k] = s[0]
    z_vmec[1, k] = theta_vmec[0]
    z_vmec[2, k] = varphi_vmec[0]


for k in np.arange(var_tips.shape[1]):
    s[0] = var_tips[0, k]
    theta[0] = var_tips[1, k]
    varphi[0] = var_tips[2, k]
    if s[0] <= 0.0 or s[0] >= 1.0:
        print(s[0])
        break
    libneo_orb._lib.can_to_vmec_(s, theta, varphi, theta_vmec, varphi_vmec)
    tip_vmec[0, k] = s[0]
    tip_vmec[1, k] = theta_vmec[0]
    tip_vmec[2, k] = varphi_vmec[0]

#%%
#
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111, projection='3d')
ax._axis3don = False
R0 = 10
phi = np.linspace(0, 2 * np.pi, endpoint=True, num=75)
theta = np.linspace(0, 2 * np.pi, endpoint=True, num=75)
theta, phi = np.meshgrid(theta, phi)
x = np.cos(phi) * (R0 + np.cos(theta))
y = np.sin(phi) * (R0 + np.cos(theta))
z = np.sin(theta)
ax.plot_surface(x, y, z, rstride=2, cstride=2, alpha=0.1)
phase = -np.pi
R = np.sqrt(zs[0, :])*np.cos(zs[1, :])
Z = np.sqrt(zs[0, :])*np.sin(zs[1, :])
P = zs[2, :] + phase
ax.plot((R0+R)*np.cos(P), -(R0+R)*np.sin(P),
        Z, linewidth=0.5, alpha=0.7, color='k')
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_zlim(-3, 3)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.tight_layout()

#%%    
np.save('z.npy', z_vmec)
np.save('ztip.npy', tip_vmec)


# %%
