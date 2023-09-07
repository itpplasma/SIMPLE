#%%
import time
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from neo_orb_fffi import *
from os import path

import exportfig

#basedir = '/home/calbert/net/cobra/ipp/ishw19/scripts'
basedir = '.'

zs = np.load(path.join(basedir, 'z.npy'))
var_tips = np.load(path.join(basedir, 'ztip.npy'))

#%%
libneo_orb.compile(verbose=1)
neo_orb.load()
new_vmec_stuff.load()
neo_orb.init_field(5, 5, 3, -1)
#%%
s = libneo_orb._ffi.new('double*', 0.7)
theta = libneo_orb._ffi.new('double*', 0.0)
varphi = libneo_orb._ffi.new('double*', 0.0)
theta_vmec = libneo_orb._ffi.new('double*', 0.0)
varphi_vmec = libneo_orb._ffi.new('double*', 0.0)
A_phi = libneo_orb._ffi.new('double*', 0.0)
A_theta = libneo_orb._ffi.new('double*', 0.0)
dA_phi_ds = libneo_orb._ffi.new('double*', 0.0)
dA_theta_ds = libneo_orb._ffi.new('double*', 0.0)
aiota = libneo_orb._ffi.new('double*', 0.0)
R = libneo_orb._ffi.new('double*', 0.0)
Z = libneo_orb._ffi.new('double*', 0.0)
alam = libneo_orb._ffi.new('double*', 0.0)
dR_ds = libneo_orb._ffi.new('double*', 0.0)
dR_dt = libneo_orb._ffi.new('double*', 0.0)
dR_dp = libneo_orb._ffi.new('double*', 0.0)
dZ_ds = libneo_orb._ffi.new('double*', 0.0)
dZ_dt = libneo_orb._ffi.new('double*', 0.0)
dZ_dp = libneo_orb._ffi.new('double*', 0.0)
dl_ds = libneo_orb._ffi.new('double*', 0.0)
dl_dt = libneo_orb._ffi.new('double*', 0.0)
dl_dp = libneo_orb._ffi.new('double*', 0.0)
Bctrvr_vartheta = libneo_orb._ffi.new('double*', 0.0)
Bctrvr_varphi = libneo_orb._ffi.new('double*', 0.0)
Bcovar_r = libneo_orb._ffi.new('double*', 0.0)
Bcovar_vartheta = libneo_orb._ffi.new('double*', 0.0)
Bcovar_varphi = libneo_orb._ffi.new('double*', 0.0)
sqg = libneo_orb._ffi.new('double*', 0.0)

# %%

phase = -np.pi

RR = np.zeros_like(zs[0, :])
ZZ = np.zeros_like(RR)
PP = np.zeros_like(RR)

for k in np.arange(len(RR)):
    s[0] = zs[0, k]
    theta[0] = zs[1, k]
    varphi[0] = zs[2, k]
    if s[0] <= 0.0 or s[0] >= 1.0:
        break
    neo_orb.spline_vmec(
        s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, R, Z,
        alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp
    )
    RR[k] = R[0]/100.0
    ZZ[k] = Z[0]/100.0
    PP[k] = varphi[0]

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111, projection='3d')
ax._axis3don = False

ph = np.linspace(0, 2 * np.pi, endpoint=True, num=75)
th = np.linspace(0, 2 * np.pi, endpoint=True, num=75)
th, ph = np.meshgrid(th, ph)
thflat = th.flatten()
phflat = ph.flatten()
x = np.empty_like(thflat)
y = np.empty_like(thflat)
z = np.empty_like(thflat)

for k in np.arange(len(thflat)):
    s[0] = 1.0
    theta[0] = thflat[k]
    varphi[0] = phflat[k]
    neo_orb.spline_vmec(
        s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, R, Z,
        alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp
    )
    x[k] = np.cos(-varphi[0]) * R[0]/100.0
    y[k] = np.sin(-varphi[0]) * R[0]/100.0
    z[k] = Z[0]/100.0

x = x.reshape(th.shape)
y = y.reshape(th.shape)
z = z.reshape(th.shape)

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

RR = np.zeros_like(var_tips[0, :])
ZZ = np.zeros_like(RR)
PP = np.zeros_like(RR)

for k in np.arange(len(RR)):
    s[0] = var_tips[0, k]
    theta[0] = var_tips[1, k]
    varphi[0] = var_tips[2, k]
    if s[0] <= 0.0 or s[0] >= 1.0:
        break
    neo_orb.spline_vmec(
        s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, R, Z,
        alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp
    )
    RR[k] = R[0]/100.0
    ZZ[k] = Z[0]/100.0
    PP[k] = varphi[0]


ax.plot(RR*np.cos(-PP), RR*np.sin(-PP),
        ZZ, 'o', markersize=1, color='tab:red')
plt.savefig('orbit2.png', dpi=150)
exportfig.exportfig('orbit2')
exportfig.exporteps('orbit2')
# %%
fig = plt.figure(figsize=(2.0, 2.0))
#ax = fig.add_subplot(111, projection='3d')

#ax._axis3don = False
th = np.linspace(0, 2*np.pi, 100)
#plt.plot(RR, ZZ, ',')


thring = np.linspace(0, 2 * np.pi, endpoint=True, num=len(PP))
x = np.empty_like(thring)
y = np.empty_like(thring)
z = np.empty_like(thring)

for k in np.arange(len(RR)):
    s[0] = 1.0
    theta[0] = var_tips[1, k]
    varphi[0] = np.mod(var_tips[2, k]+2.0*np.pi/10.0, 2.0*np.pi/5.0)-2.0*np.pi/10.0
    neo_orb.spline_vmec(
        s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, R, Z,
        alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp
    )
    x[k] = R[0]/100.0
    #y[k] = np.sin(-varphi[0]) * R[0]/100.0
    z[k] = Z[0]/100.0
 
plt.plot(x, z, '.', markersize=1, color='tab:blue')
plt.plot(RR, ZZ, '.', markersize=0.5, color='tab:red')

plt.xlabel(r'$R$ / m', labelpad=2)
plt.ylabel(r'$Z$ / m', labelpad=-2)
plt.axis('equal')
plt.tight_layout()
plt.savefig('orbit2_proj.png', dpi=300)
exportfig.exporteps('orbit2_proj')

# %%
plt.figure(figsize=(1.8, 1.8))
plt.plot(np.cos(thring), np.sin(thring), '-', color='tab:blue')

plt.plot(np.sqrt(var_tips[0, :])*np.cos(var_tips[1, :]), 
         np.sqrt(var_tips[0, :])*np.sin(var_tips[1, :]),
         'o', markersize=0.5, color='tab:red')
plt.xlabel(r'$R$')
plt.ylabel(r'$Z$')
plt.axis('equal')
plt.tight_layout()
plt.savefig('orbit2_proj_topo.png', dpi=300)

exportfig.exporteps('orbit2_proj_topo')


# %%
