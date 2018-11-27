!
  use new_vmec_stuff_mod, only : netcdffile,multharm,ns_A,ns_s,ns_tp
  use vmec_stuff_mod, only : kpar
  use chamb_mod,  only : rbig,rcham2
  use parmot_mod, only : rmu,ro0,eeff
  use velo_mod,   only : isw_field_type
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
  double precision,parameter  :: c=2.9979d10
  double precision,parameter  :: e_charge=4.8032d-10
  double precision,parameter  :: e_mass=9.1094d-28
  double precision,parameter  :: p_mass=1.6726d-24
  double precision,parameter  :: ev=1.6022d-12
!
  integer          :: npoi,ierr,L1i,nper,npoiper,i,ntimstep,ntestpart
  integer          :: ipart,notrace_passing,loopskip,iskip,ilost
  real             :: zzg
  double precision :: dphi,rbeg,phibeg,zbeg,bmod00,rcham,rlarm,bmax,bmin
  double precision :: tau,dtau,dtaumin,xi,v0,bmod_ref,E_alpha,trace_time
  double precision :: RT0,R0i,cbfi,bz0i,bf0,trap_par
  double precision :: sbeg,thetabeg
  double precision, dimension(5) :: z
! 24.03.2016
  integer          :: npoiper2
! 24.03.2016 end
! 13.02.2013
  double precision :: contr_pp
! 14.04.2013
  double precision :: facE_al
  integer          :: ibins
! 14.04.2013 end
! 22.09.2013
  integer          :: n_e,n_d,n_b
! 22.09.2013 end
! 13.02.2013 end
  double precision :: r,vartheta_c,varphi_c,theta_vmec,varphi_vmec,alam0
!
  open(1,file='alpha_lifetime_m.inp')
  read (1,*) notrace_passing   !skip tracing passing prts if notrace_passing=1
  read (1,*) nper              !number of periods for initial field line
  read (1,*) npoiper           !number of points per period on this field line
  read (1,*) ntimstep          !number of time steps per slowing down time
  read (1,*) ntestpart         !number of test particles
!<=2017  read (1,*) rcham             !vacuum chamber radius, cm
  read (1,*) bmod_ref          !reference field, G, for Boozer $B_{00}$
  read (1,*) trace_time        !slowing down time, s
  read (1,*) sbeg              !starting s for field line                       !<=2017
  read (1,*) phibeg            !starting phi for field line                     !<=2017
  read (1,*) thetabeg          !starting theta for field line                   !<=2017
  read (1,*) loopskip          !how many loops to skip to shift random numbers
! 13.02.2013
  read (1,*) contr_pp          !control of passing particle fraction
! 13.02.2013 end
! 14.04.2013
  read (1,*) facE_al           !facE_al test particle energy reduction factor
! 14.04.2013
! 24.03.2016
  read (1,*) npoiper2          !additional integration step split factor
! 24.03.2016 end
! 22.09.2013
  read (1,*) n_e               !test particle charge number (the same as Z)
  read (1,*) n_d               !test particle mass number (the same as A)
! 22.09.2013 end
! 26.09.2013
!<=2017  read (1,*) n_b               !determination bmax
  read (1,*) netcdffile        !name of VMEC file in NETCDF format <=2017 NEW
  close(1)
!
  multharm=7
  ns_A=5
  ns_s=5
  ns_tp=5
!
  call spline_vmec_data
!
  call get_canonical_coordinates
!
  call testing
!
! inverse relativistic temperature
  rmu=1d8
!
! 14.04.2013
! alpha particle energy, eV:
!  E_alpha=3.5d6/16d0
  E_alpha=3.5d6/facE_al
! alpha particle velocity, cm/s
! 22.09.2013  v0=sqrt(2.d0*E_alpha*ev/(4.d0*p_mass))
  v0=sqrt(2.d0*E_alpha*ev/(n_d*p_mass))
! 14.04.2013 end
!
! Larmor radius:
!  rlarm=v0*4.d0*p_mass*c/(e_charge*bmod_ref)
! 22.09.2013  rlarm=v0*2.d0*p_mass*c/(e_charge*bmod_ref)
  rlarm=v0*n_d*p_mass*c/(n_e*e_charge*bmod_ref)
! normalized slowing down time:
  tau=trace_time*v0
! normalized time step:
  dtau=tau/dfloat(ntimstep-1)
!
! 14.11.2011  call stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0)
  call stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0)         !<=2017
!<=2017  call tj2vvo(RT0,R0i,L1i,cbfi,bz0i,bf0)
! 05.03.2014 for w7x
! 13.07.2016  rt0=145d0
! 13.07.2016  L1i=3
! 05.03.2014 end
!
! parameters for the vacuum chamber:
  rbig=rt0
!<=2017  rcham2=rcham**2
! field line integration step step over phi (to check chamber wall crossing)
  dphi=2.d0*pi/(L1i*npoiper)
! orbit integration time step (to check chamber wall crossing)
! 24.03.2016  dtaumin=dphi*rbig
  dtaumin=dphi*rbig/npoiper2!
dtau=2*dtaumin
print *,dtau
bmod00=281679.46317784750d0
!
! Larmor raidus corresponds to the field stregth egual to $B_{00}$ harmonic
! in Boozer coordinates:
! 14.11.2011  bmod00=bmod_ref  !<=deactivated, use value from the 'alpha_lifetime.inp'
  ro0=rlarm*bmod00  ! 23.09.2013
!
  isw_field_type=1
  r=0.7d0
  vartheta_c=0.5d0
  varphi_c=0.5d0
  alam0=0.5d0
!
  call can_to_vmec(r,vartheta_c,varphi_c,theta_vmec,varphi_vmec)
!
  print *,'can : ',r,vartheta_c,varphi_c
  print *,'VMEC: ',r,theta_vmec,varphi_vmec
!
  z(1)=r
  z(2)=theta_vmec
  z(3)=varphi_vmec
  z(4)=1.d0
  z(5)=alam0
!
print *,'VMEC, splines'
  do i=1,L1i*npoiper*npoiper2*10
!
    call orbit_timestep_can(z,dtau,dtaumin,ierr)
!
    write (3001,*) dtau*dfloat(i),z
  enddo
print *,'done'
!
  z(1)=r*98.d0*cos(theta_vmec)
  z(2)=varphi_vmec
  z(3)=r*98.d0*sin(theta_vmec)
  z(4)=1.d0
  z(5)=alam0
!
print *,'VMEC, orig'
  do i=1,L1i*npoiper*npoiper2*10
!
    call orbit_timestep(z,dtau,dtaumin,ierr)
!
    write (3002,*) dtau*dfloat(i),sqrt(z(1)**2+z(3)**2)/98.d0,atan2(z(3),z(1)),z(2),z(4:5)
  enddo
print *,'done'
!
  isw_field_type=0
  z(1)=r
  z(2)=vartheta_c
  z(3)=varphi_c
  z(4)=1.d0
  z(5)=alam0
!
print *,'canonical'
  do i=1,L1i*npoiper*npoiper2*10
!
    call orbit_timestep_can(z,dtau,dtaumin,ierr)
!
    call can_to_vmec(z(1),z(2),z(3),theta_vmec,varphi_vmec)
!
!
    write (3003,*) dtau*dfloat(i),z(1),theta_vmec,varphi_vmec,z(4:5),z(2:3)
  enddo
print *,'done'
!
  end
!
