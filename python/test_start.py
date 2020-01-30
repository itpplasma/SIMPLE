"""
Created: Tue Aug 20 15:07:28 2019
@author: Christopher Albert <albert@alumni.tugraz.at>
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from fffi import fortran_library, fortran_module

libneo_orb = fortran_library('simple')

# %%
neoorb = fortran_module(libneo_orb, 'neo_orb_global')

neoorb.fdef(
    """
    subroutine init_field(ans_s, ans_tp, amultharm, aintegmode)
        integer, intent(in) :: ans_s, ans_tp, amultharm, aintegmode
    end

    subroutine spline_vmec(s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)
  
    double precision, intent(in) :: s, theta, varphi
    double precision, intent(out) :: A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp

    end
    """
)

# %%
libneo_orb.compile(verbose=10)
neoorb.load()

# %%
# %%
neoorb.init_field(5, 5, 3, -1)
# %%
s = libneo_orb._ffi.new('double*', 0.25)
theta = libneo_orb._ffi.new('double*', 0.0)
varphi = libneo_orb._ffi.new('double*', 0.0)
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
# %%
data = np.loadtxt(
    '/home/calbert/net/cobra/run/SIMPLE/QA/2020-01-22_s0.25_vmec/start.dat')
ss = data[:, 0]
th = data[:, 1]
ph = data[:, 2]
PITCH = data[:, 4]

RR = np.zeros_like(ss)
ZZ = np.zeros_like(ss)
PP = np.zeros_like(ss)

for k in np.arange(len(ss)):
    s = ss[k]
    theta = th[k]
    varphi = ph[k]
    neoorb.spline_vmec(s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota,
                       R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)
    RR[k] = R[0]/100.0
    ZZ[k] = Z[0]/100.0
    PP[k] = varphi

# for pseudo-cylinder coordinates
# RR = 1.0 + s*np.cos(th)
# ZZ = s*np.sin(th)
# PP = ph

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.plot(RR*np.cos(PP), RR*np.sin(PP), ZZ, ',')

plt.figure()
plt.plot(PITCH, ',')


# %%
data = np.loadtxt('launch_0.25', skiprows=1)
#data = np.loadtxt('launch_hydra_a12_6_hi_7_bf_0.9.txt', skiprows=1)
i = data[:, 0]
E = data[:, 1]
r = data[:, 3]
phi = data[:, 4]
z = data[:, 5]
vr = data[:, 6]
vphi = data[:, 7]
vz = data[:, 8]
pitch = data[:, 12]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.plot(r*np.cos(phi), r*np.sin(phi), z, 'bo', markersize=1)
plt.plot(RR*np.cos(PP), RR*np.sin(PP), ZZ, 'rs', markersize=1)

plt.figure()
plt.plot(pitch, ',')

# %%
plt.show()
