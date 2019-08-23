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

field_can = fortran_module(libneo_orb, 'field_can_mod')
field_can.fdef("""
type :: FieldCan
  integer :: field_type

  double precision :: Ath, Aph
  double precision :: hth, hph
  double precision :: Bmod

  double precision, dimension(3) :: dAth, dAph
  double precision, dimension(3) :: dhth, dhph
  double precision, dimension(3) :: dBmod

  double precision, dimension(6) :: d2Ath, d2Aph
  double precision, dimension(6) :: d2hth, d2hph
  double precision, dimension(6) :: d2Bmod

  double precision :: H, pth, vpar
  double precision, dimension(4) :: dvpar, dH, dpth
  
  double precision, dimension(10) :: d2vpar, d2H, d2pth

  double precision :: mu, ro0
end type               

subroutine fieldcan_init(f, mu, ro0, vpar, field_type)
  type(FieldCan), intent(inout) :: f
  double precision, intent(in) :: mu, ro0, vpar
  integer, intent(in) :: field_type
end

subroutine eval_field(f, r, th_c, ph_c, mode_secders)
  type(FieldCan), intent(inout) :: f
  double precision, intent(in) :: r, th_c, ph_c
  integer, intent(in) :: mode_secders
end

subroutine get_derivatives(f, pphi)
  type(FieldCan), intent(inout) :: f
  double precision, intent(in) :: pphi
end
""")

orbit_symplectic = fortran_module(libneo_orb, 'orbit_symplectic')
orbit_symplectic.fdef("""  
type :: SymplecticIntegrator
  integer :: nlag
  integer :: nbuf
  logical :: extrap_field

  double precision :: atol
  double precision :: rtol

  double precision, dimension(4) :: z
  double precision :: pthold

  integer :: kbuf
  integer :: kt
  integer :: k
  integer :: bufind1, bufind2, bufind3
  double precision, dimension(264) :: zbuf
  double precision, dimension(3) :: coef

  integer :: ntau
  double precision :: dt
  double precision :: pabs

  integer :: mode
end type

subroutine orbit_sympl_init(si, f, z, dt, ntau, rtol_init, mode_init, nlag)
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f
  double precision, intent(in) :: z(:)
  double precision, intent(in) :: dt
  integer, intent(in) :: ntau
  double precision, intent(in) :: rtol_init
  integer, intent(in) :: mode_init
  integer, intent(in) :: nlag
end

subroutine orbit_timestep_sympl(si, f, ierr)
  type(SymplecticIntegrator), intent(inout) :: si
  type(FieldCan), intent(inout) :: f
  integer, intent(out) :: ierr
end
""")

libneo_orb.compile()
field_can.load()
orbit_symplectic.load()


# Initial conditions
z0 = np.array([0.1, 1.5, 0.0, 0.0])
vpar0 = 0.0

f = field_can.new('FieldCan')
field_can.fieldcan_init(f, 1.0e-5, 1.0, vpar0, -1)
field_can.eval_field(f, z0[0], z0[1], z0[2], 0)
z0[3] = vpar0*f.hph + f.Aph/f.ro0  # pphi

taub = 1400.0  # estimated bounce time
npoiper = 16
dt = taub/npoiper
integ_mode = 1
rtol = 1e-13
nlag = 0

si = orbit_symplectic.new('SymplecticIntegrator')

#%%
# Grid in th,r space
ths = np.linspace(0.5, 1.5, 30)
rs = np.linspace(0.05, 0.2, 30)
[TH, R] = np.meshgrid(ths, rs); TH=TH.flatten(); R=R.flatten()
PTH  = np.empty(TH.shape)

nt = 128
zt = np.empty([4, len(TH), nt]) 
PTHt = np.empty([len(TH), nt]) 

# Grid in th,pth space
for k in range(len(TH)):
    z0[0] = R[k]
    z0[1] = TH[k]
    z0[2] = 0.0
    field_can.eval_field(f, z0[0], z0[1], z0[2], 0)
    # z0[3] = vpar0*f.hph + f.Aph/f.ro0  # keep pphi
    field_can.get_derivatives(f, z0[3])
    PTH[k] = f.pth
    
    zt[:,k,0] = z0
    orbit_symplectic.orbit_sympl_init(si, f, z0, dt, 1, rtol, integ_mode, nlag)
    for kt in range(nt):
        ierr = 0
        orbit_symplectic.orbit_timestep_sympl(si, f, ierr)
        zt[0, k, kt] = si.z[0]
        zt[1, k, kt] = si.z[1]
        zt[2, k, kt] = si.z[2]
        zt[3, k, kt] = si.z[3]
        PTHt[k, kt] = f.pth
    
plt.figure()
plt.plot(TH, PTH, '.')
t1 = 4
plt.plot(zt[1,:,t1], PTHt[:,t1], '.')
t1 = 8
plt.plot(zt[1,:,t1], PTHt[:,t1], '.')
t1 = 16
plt.plot(zt[1,:,t1], PTHt[:,t1], '.')
plt.xlabel('th')
plt.ylabel('pth')
plt.legend([0,4,8,16])

#%%
plt.figure()
plt.plot(zt[0,:,:].T*np.cos(zt[1,:,:]).T, zt[0,:,:].T*np.sin(zt[1,:,:]).T, '-')
plt.show()
#%%

