# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 16:10:06 2019

@author: chral
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

orbit_symplectic = fortran_module(libneo_orb, 'orbit_symplectic')
orbit_symplectic.fdef("""
  type :: SymplecticIntegrator
    double precision :: atol
    double precision :: rtol
  
  ! Current phase-space coordinates z and old pth
    double precision, dimension(4) :: z  ! z = (r, th, ph, pphi)
    double precision :: pthold
  
    ! Buffer for Lagrange polynomial interpolation
    integer :: kbuf
    integer :: kt
    integer :: k
    integer :: bufind(0:nlag)
    double precision, dimension(4, nbuf) :: zbuf
    double precision, dimension(0:0, nlag+1) :: coef
  
    ! Timestep and variables from z0
    integer :: ntau
    double precision :: dt
    double precision :: pabs
  
    ! Integrator mode
    integer :: mode  ! 1 = euler1, 2 = euler2, 3 = verlet
  end type SymplecticIntegrator                      
""")
