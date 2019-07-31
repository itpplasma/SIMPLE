#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 12:35:03 2019

@author: calbert
"""
import time
import numpy as np
import matplotlib.pyplot as plt
from fffi import fortran_library, fortran_module

pi=3.14159265358979
c=2.9979e10
e_charge=4.8032e-10
e_mass=9.1094e-28
p_mass=1.6726e-24
ev=1.6022e-12

bmod_ref = 5.0e4
Z_charge = 2
m_mass = 4
E_kin = 3.5e6
trace_time = 1.0e-2
#rbig=rlarm*bmod00
#dphi=2.d0*pi/(L1i*npoiper)

v0 = np.sqrt(2.0*E_kin*ev/(m_mass*p_mass))    # alpha particle velocity, cm/s
rlarm=v0*m_mass*p_mass*c/(Z_charge*e_charge*bmod_ref) # reference gyroradius
tau=trace_time*v0

libneo_orb = fortran_library('neo_orb', path='../lib')

neo_orb = fortran_module(libneo_orb, 'neo_orb')
neo_orb.fdef("""
  double precision :: fper
  double precision :: dtau, dtaumax, v0
  double precision, dimension(5) :: z
  integer          :: n_e, n_d
  integer          :: firstrun
  
  subroutine init_field(ans_s, ans_tp, amultharm, aintegmode)
    integer :: ans_s, ans_tp, amultharm, aintegmode
  end
  
  subroutine init_params(Z_charge, m_mass, E_kin, adtau, adtaumax, arelerr)
    integer :: Z_charge, m_mass
    double precision :: E_kin, adtau, adtaumax
    double precision :: arelerr
  end
  
  subroutine timestep_global(ierr)
    integer :: ierr
  end
""")

new_vmec_stuff = fortran_module(libneo_orb, 'new_vmec_stuff_mod')
new_vmec_stuff.fdef("""
    double precision :: rmajor, h_theta, h_phi
""")


cut_detector = fortran_module(libneo_orb, 'cut_detector')
cut_detector.fdef("""
    subroutine init()
    end
    
    subroutine trace_to_cut(t, var_tip, cut_type, ierr)
      double precision, intent(inout) :: t
      double precision, dimension(:), intent(inout) :: var_tip
      integer, intent(out) :: cut_type
      integer, intent(out) :: ierr
    end
""")

libneo_orb.compile()
neo_orb.load()
new_vmec_stuff.load()

#%%
neo_orb.init_field(5, 5, 3, 0)

#%%
npoiper2 = 64                       # interation points per field period
rbig = new_vmec_stuff.rmajor*1.0e2  # major radius in cm
dtaumax = 2.0*pi*rbig/npoiper2      # maximum time step for integrator
dtau = dtaumax                   # time step for output
neo_orb.init_params(Z_charge, m_mass, E_kin, dtau, dtaumax, 1e-8)
#%%
cut_detector.load()
cut_detector.init()

#%%
ntimstep = 10000
neo_orb.firstrun = 1  # to re-initialize symplectic integrator in case

s = 0.5
th = 0.0
ph = 0.314
lam = 0.22

z = neo_orb.z
z[:] = [s,th,ph,1.0,lam]

zs = np.empty([4, ntimstep])
zs[:,0] = z[[0,1,2,4]]

t = time.time()
ierr = 0  # doesn't work yet with pass-by-reference
for kt in range(1, ntimstep):
    neo_orb.timestep_global(ierr)
    zs[:, kt] = z[[0,1,2,4]]
print(z)
print('time elapsed: {} s'.format(time.time() - t))

t = 0.0
ierr = 0
var_tip = np.zeros(6)
ncut = 1000
# variables to evaluate at tip: z(1..5), par_inv
var_tips = np.empty([6, ncut])

t = time.time()
for k in np.arange(ncut):
    cut_type = 0
    cut_detector.trace_to_cut(0.0, var_tip, cut_type, ierr)
    var_tips[:,k] = var_tip
print(z)
print('time elapsed: {} s'.format(time.time() - t))

plt.plot(zs[0,:]*np.cos(zs[1,:]), zs[0,:]*np.sin(zs[1,:]))
plt.plot(var_tips[0,:]*np.cos(var_tips[1,:]), var_tips[0,:]*np.sin(var_tips[1,:]), ',')
