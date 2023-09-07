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
from fffi import fortran_library, fortran_module
from netCDF4 import Dataset

from mpl_toolkits.mplot3d import Axes3D

pi = np.pi

libneo_orb = fortran_library(
    'simple', compiler={'name': 'gfortran', 'version': 9})
# compiler={'name': 'ifort', 'version': 18})

neo_orb = fortran_module(libneo_orb, 'neo_orb_global')
neo_orb.fdef("""
  subroutine init_field(ans_s, ans_tp, amultharm, aintegmode)
    integer :: ans_s, ans_tp, amultharm, aintegmode
  end

  subroutine init_params(Z_charge, m_mass, E_kin, adtau, adtaumax, arelerr)
    integer :: Z_charge, m_mass
    double precision :: E_kin, adtau, adtaumax
    double precision :: arelerr
  end

  subroutine init_integrator(z0)
    double precision, dimension(:), intent(in) :: z0
  end

  subroutine timestep_z(z, ierr)
    double precision, dimension(:) :: z
    integer :: ierr
  end

  subroutine timestep_sympl_z(z, ierr)
    double precision, dimension(:) :: z
    integer :: ierr
  end

  subroutine spline_vmec(s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)
    double precision, intent(in) :: s, theta, varphi
    double precision, intent(out) :: A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
  end

  subroutine field_vmec(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds, aiota, sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi, Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
    double precision :: s,theta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota, R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
    double precision :: Bctrvr_vartheta,Bctrvr_varphi,Bcovar_r,Bcovar_vartheta,Bcovar_varphi,sqg
  end
""")

new_vmec_stuff = fortran_module(libneo_orb, 'new_vmec_stuff_mod')
new_vmec_stuff.fdef("""
    double precision :: rmajor, h_theta, h_phi
    integer :: nsurfm,nstrm,nper,kpar
""")

libneo_orb.compile(verbose=1)
neo_orb.load()
new_vmec_stuff.load()
# %%
neo_orb.init_field(5, 5, 3, -1)

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
Bctrvr_vartheta = libneo_orb._ffi.new('double*', 0.0)
Bctrvr_varphi = libneo_orb._ffi.new('double*', 0.0)
Bcovar_r = libneo_orb._ffi.new('double*', 0.0)
Bcovar_vartheta = libneo_orb._ffi.new('double*', 0.0)
Bcovar_varphi = libneo_orb._ffi.new('double*', 0.0)
sqg = libneo_orb._ffi.new('double*', 0.0)

# %%
rootgrp = Dataset("wout.nc", "r", format="NETCDF4")
ks = 20

sdat = rootgrp.variables['phi']/max(rootgrp.variables['phi'])
s0 = (sdat[ks] + sdat[ks-1])/2.0
#s0 = sdat[ks]
print('s = {}'.format(s0))

nth = 144+1
nph = 72+1
ss = s0*np.ones(nth*nph).copy()
th = np.linspace(-pi/2, 3*pi/2, nth, endpoint=True)
ph = np.linspace(0, pi, nph, endpoint=True)
[PH, TH] = np.meshgrid(ph, th)
TH = TH.flatten()
PH = PH.flatten()
TH2 = np.zeros_like(TH)

RR = np.zeros_like(ss)
ZZ = np.zeros_like(ss)
PP = np.zeros_like(ss)
BB = np.zeros_like(ss)
# %%
for k in np.arange(len(ss)):
    s[0] = ss[k]
    theta[0] = TH[k]
    varphi[0] = PH[k]
    #print(s[0], theta[0], varphi[0])
    neo_orb.spline_vmec(
        s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, R, Z,
        alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp
    )
    neo_orb.field_vmec(
        s, theta, varphi, A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, sqg,
        alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, Bcovar_r,
        Bcovar_vartheta, Bcovar_varphi
    )

    #print(R[0], Z[0])
    RR[k] = R[0]/100.0
    ZZ[k] = Z[0]/100.0
    PP[k] = varphi[0]
    BB[k] = np.sqrt(Bctrvr_vartheta[0]*Bcovar_vartheta[0]
                    + Bctrvr_varphi[0]*Bcovar_varphi[0])
    TH2[k] = TH[k] + alam[0]
    #BB[k] = Bcovar_vartheta[0]

fig = plt.figure(figsize=(3.2, 3.2))
ax = fig.add_subplot(111, projection='3d')
plt.plot(RR*np.cos(PP), RR*np.sin(PP), ZZ, 'k,')
max_range = np.array([2*RR.max(), 2*RR.max(), ZZ.max()-ZZ.min()]).max() / 2.0

ax.set_xlim(- max_range, max_range)
ax.set_ylim(- max_range, max_range)
ax.set_zlim(- max_range, max_range)


# %%
plt.figure(figsize=(3.2, 3.2))
plt.contour(PH.reshape(nth, nph)/pi, TH.reshape(nth, nph)/pi,
            BB.reshape(nth, nph)/min(BB), cmap='jet', levels=100)
plt.xlabel(r'$\varphi/\pi$')
plt.ylabel(r'$\vartheta/\pi$')

plt.xticks(np.arange(0.0, 1.01, 0.5))
plt.yticks(np.arange(-0.5, 1.51, 0.5))

print(new_vmec_stuff.nsurfm)


# %%

mn_mode_nyq = rootgrp.dimensions['mn_mode_nyq']
xm_nyq = rootgrp.variables['xm_nyq']
xn_nyq = rootgrp.variables['xn_nyq']
xm = rootgrp.variables['xm']
xn = rootgrp.variables['xn']
bmnc = rootgrp.variables['bmnc']
bsubumnc = rootgrp.variables['bsubumnc']
bsubvmnc = rootgrp.variables['bsubvmnc']
bsupumnc = rootgrp.variables['bsupumnc']
bsupvmnc = rootgrp.variables['bsupvmnc']


b = np.zeros(TH.shape)
bsubu = np.zeros(TH.shape)
bsubv = np.zeros(TH.shape)
bsupu = np.zeros(TH.shape)
bsupv = np.zeros(TH.shape)

mmin = min(xm)-0.5
mmax = max(xm)+0.5
nmin = min(xn)-0.5
nmax = max(xn)+0.5

for k in np.arange(mn_mode_nyq.size):
    m = xm_nyq[k]
    n = xn_nyq[k]
    if m<mmin or m>mmax or n<nmin or n>nmax:
        continue

    fac = np.cos(m*TH - n*PH)
    b += bmnc[ks, k]*fac
    bsubu += bsubumnc[ks, k]*fac
    bsubv += bsubvmnc[ks, k]*fac
    bsupu += bsupumnc[ks, k]*fac
    bsupv += bsupvmnc[ks, k]*fac

# %%

plt.figure(figsize=(3.2, 3.2))
plt.contour(PH.reshape(nth, nph)/pi, TH.reshape(nth, nph)/pi,
            b.reshape(nth, nph)/min(b), cmap='jet', levels=100)
plt.xlabel(r'$\varphi/\pi$')
plt.ylabel(r'$\vartheta/\pi$')

plt.xticks(np.arange(0.0, 1.01, 0.5))
plt.yticks(np.arange(-0.5, 1.51, 0.5))

# %%
b2 = np.sqrt(bsubu*bsupu + bsubv*bsupv)

plt.figure(figsize=(3.2, 3.2))
plt.contour(PH.reshape(nth, nph)/pi, TH.reshape(nth, nph)/pi,
            b2.reshape(nth, nph)/min(b), cmap='jet', levels=100)
plt.xlabel(r'$\varphi/\pi$')
plt.ylabel(r'$\vartheta/\pi$')

plt.xticks(np.arange(0.0, 1.01, 0.5))
plt.yticks(np.arange(-0.5, 1.51, 0.5))

rootgrp.close()


# %%
