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
from pysimple import simple, params, orbit_symplectic, cut_detector, \
    vmec_to_can, can_to_vmec, vmec_to_cyl, spline_vmec_data

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
s = 0.7
theta = 0.0
varphi = 0.0
# %%
theta_vmec, varphi_vmec = can_to_vmec(s, theta, varphi)
R, Z = vmec_to_cyl(s, theta_vmec, varphi_vmec)
print("Initial canonical coordinates (s,th_c,ph_c) = ", s, theta, varphi)
print("Initial VMEC coordinates (s,th,ph) = ", s, theta_vmec, varphi_vmec)
print("Initial cylindrical coordinates (R,PH,Z) = ", R, varphi_vmec, Z)

# %%
ntimstep = 10000

z = np.array([s0, th0, ph0, 1.0, lam])
simple.init_integrator(tracy, z)

zs = np.empty([4, ntimstep+1])
zs[:, 0] = z[[0, 1, 2, 4]]

t = time.time()
for kt in range(ntimstep):
    ierr = simple._timestep_sympl_z(tracy.si, tracy.f, z)
    zs[:, kt+1] = z[[0, 1, 2, 4]]
print(z)
print('time elapsed: {} s'.format(time.time() - t))
kmax = ntimstep

fig, ax = plt.subplots()
ax.plot(zs[0,:]*np.sin(zs[1,:]), zs[0,:]*np.cos(zs[1,:]))


# %%
z = np.array([s0, th0, ph0, 1.0, lam])
simple.init_integrator(tracy, z)

cutty = cut_detector.CutDetector()
cut_detector.init(cutty, tracy.fper, z)
ierr = 0
var_cut = np.zeros(6)
ncut = 500
# variables to evaluate at tip: z(1..5), par_inv
var_tips = np.zeros([6, ncut])
#
t = time.time()
for k in np.arange(ncut):
    cut_type, ierr = cut_detector.trace_to_cut(cutty, tracy.si, tracy.f, z, var_cut)
    var_tips[:, k] = var_cut
print('time elapsed: {} s'.format(time.time() - t))

fig, ax = plt.subplots()
ax.plot(var_tips[0,:]*np.sin(var_tips[1,:]), var_tips[0,:]*np.cos(var_tips[1,:]), ',')

#%%

z_vmec = np.zeros_like(zs)
tip_vmec = np.zeros_like(var_tips)

for k in np.arange(zs.shape[1]):
    s = zs[0, k]
    theta = zs[1, k]
    varphi = zs[2, k]
    if s <= 0.0 or s >= 1.0:
        print(s)
        break

    z_vmec[0, k] = s
    z_vmec[1, k], z_vmec[2, k] = can_to_vmec(s, theta, varphi)


for k in np.arange(var_tips.shape[1]):
    s = var_tips[0, k]
    theta = var_tips[1, k]
    varphi = var_tips[2, k]
    if s <= 0.0 or s >= 1.0:
        print(s)
        break
    tip_vmec[1, k], tip_vmec[2, k] = can_to_vmec(s, theta, varphi)

#%% Plot orbit in canonical coordinates stretched to a torus
#   This is not the actual geometry, but the topology is the same.
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

# %% Same for Poincare cut
plt.figure()
plt.plot(np.cos(thring), np.sin(thring), '-', color='tab:blue')

plt.plot(np.sqrt(var_tips[0, :])*np.cos(var_tips[1, :]), 
         np.sqrt(var_tips[0, :])*np.sin(var_tips[1, :]),
         'o', markersize=0.5, color='tab:red')
plt.xlabel(r'$R$')
plt.ylabel(r'$Z$')
plt.axis('equal')
plt.tight_layout()
#plt.savefig('orbit2_proj_topo.png', dpi=300)

#exportfig.exporteps('orbit2_proj_topo')

#%% Compute orbit in cylindrical coordinates   

phase = -np.pi

RR = np.zeros_like(zs[0, :])
ZZ = np.zeros_like(RR)
PP = np.zeros_like(RR)

for k in np.arange(len(RR)):
    s = zs[0, k]
    theta = zs[1, k]
    varphi = zs[2, k]
    if s <= 0.0 or s >= 1.0:
        break

    RR[k], ZZ[k] = vmec_to_cyl(s, theta, varphi)
    RR[k] = RR[k]/100.0  # cm to m
    ZZ[k] = ZZ[k]/100.0
    PP[k] = varphi



## Compute outer flux surface for 3D plot
ph = np.linspace(0, 2 * np.pi, endpoint=True, num=75)
th = np.linspace(0, 2 * np.pi, endpoint=True, num=75)
th, ph = np.meshgrid(th, ph)
thflat = th.flatten()
phflat = ph.flatten()
x = np.empty_like(thflat)
y = np.empty_like(thflat)
z = np.empty_like(thflat)


for k in np.arange(len(thflat)):
    s = 1.0
    theta = thflat[k]
    varphi = phflat[k]
    R, Z = vmec_to_cyl(s, theta, varphi)
    x[k] = np.cos(-varphi) * R/100.0
    y[k] = np.sin(-varphi) * R/100.0
    z[k] = Z/100.0

x = x.reshape(th.shape)
y = y.reshape(th.shape)
z = z.reshape(th.shape)

#%% Plot orbits and flux surfaces in 3D
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111, projection='3d')
ax._axis3don = False

ax.plot_surface(x, y, z, rstride=2, cstride=2, color='tab:blue', alpha=0.1)

ax.plot(RR*np.cos(-PP), RR*np.sin(-PP),
        ZZ, linewidth=0.5, alpha=0.7, color='k')
ax.set_xlim(-16, 16)
ax.set_ylim(-16, 16)
ax.set_zlim(-5, 5)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.tight_layout()



## Compute and plot Poincare cut (tips of banana)
RR = np.zeros_like(var_tips[0, :])
ZZ = np.zeros_like(RR)
PP = np.zeros_like(RR)

for k in np.arange(len(RR)):
    s = var_tips[0, k]
    theta = var_tips[1, k]
    varphi = var_tips[2, k]
    if s <= 0.0 or s >= 1.0:
        break

    RR[k], ZZ[k] = vmec_to_cyl(s, theta, varphi)
    RR[k] = RR[k]/100.0  # cm to m
    ZZ[k] = ZZ[k]/100.0
    PP[k] = varphi


ax.plot(RR*np.cos(-PP), RR*np.sin(-PP),
        ZZ, 'o', markersize=1, color='tab:red')

# plt.savefig('orbit2.png', dpi=150)
# exportfig.exportfig('orbit2')
# exportfig.exporteps('orbit2')

# %% Plot cross-section. TODO: Make number of device field periods a parameter
fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

#ax._axis3don = False
th = np.linspace(0, 2*np.pi, 100)
#plt.plot(RR, ZZ, ',')


thring = np.linspace(0, 2 * np.pi, endpoint=True, num=len(PP))
x = np.empty_like(thring)
y = np.empty_like(thring)
z = np.empty_like(thring)

for k in np.arange(len(RR)):
    s = 1.0
    theta = var_tips[1, k]
    varphi = np.mod(var_tips[2, k], 2.0*np.pi)
    R, Z = vmec_to_cyl(s, theta, varphi)
    x[k] = R/100.0
    #y[k] = np.sin(-varphi[0]) * R[0]/100.0
    z[k] = Z/100.0
 
plt.plot(x, z, '.', markersize=1, color='tab:blue')
plt.plot(RR, ZZ, '.', markersize=0.5, color='tab:red')

plt.xlabel(r'$R$ / m', labelpad=2)
plt.ylabel(r'$Z$ / m', labelpad=-2)
plt.axis('equal')
plt.tight_layout()

#plt.savefig('orbit2_proj.png', dpi=300)
#exportfig.exporteps('orbit2_proj')



# %%# %%
