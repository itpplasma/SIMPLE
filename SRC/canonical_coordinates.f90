!
  use new_vmec_stuff_mod, only : netcdffile,multharm,ns_A,ns_s,ns_tp
  use chamb_mod,  only : rbig,rcham2
  use parmot_mod, only : rmu,ro0,eeff
  use velo_mod,   only : isw_field_type
  use orbit_symplectic, only : orbit_sympl_init, orbit_timestep_sympl
  use field_can_mod, only : neval
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
  integer          :: npoiper2
  double precision :: contr_pp
  double precision :: facE_al
  integer          :: ibins
  integer          :: n_e,n_d,n_b
  double precision :: r,vartheta_c,varphi_c,theta_vmec,varphi_vmec,alam0
!
  open(1,file='alpha_lifetime_m.inp')
  read (1,*) notrace_passing   !skip tracing passing prts if notrace_passing=1
  read (1,*) nper              !number of periods for initial field line
  read (1,*) npoiper           !number of points per period on this field line
  read (1,*) ntimstep          !number of time steps per slowing down time
  read (1,*) ntestpart         !number of test particles
  read (1,*) bmod_ref          !reference field, G, for Boozer $B_{00}$
  read (1,*) trace_time        !slowing down time, s
  read (1,*) sbeg              !starting s for field line                       !<=2017
  read (1,*) phibeg            !starting phi for field line                     !<=2017
  read (1,*) thetabeg          !starting theta for field line                   !<=2017
  read (1,*) loopskip          !how many loops to skip to shift random numbers
  read (1,*) contr_pp          !control of passing particle fraction
  read (1,*) facE_al           !facE_al test particle energy reduction factor
  read (1,*) npoiper2          !additional integration step split factor
  read (1,*) n_e               !test particle charge number (the same as Z)
  read (1,*) n_d               !test particle mass number (the same as A)
  read (1,*) netcdffile        !name of VMEC file in NETCDF format <=2017 NEW
  close(1)
!
! inverse relativistic temperature
  rmu=1d8
!
! alpha particle energy, eV:
  E_alpha=3.5d6/facE_al
! alpha particle velocity, cm/s
  v0=sqrt(2.d0*E_alpha*ev/(n_d*p_mass))
! 14.04.2013 end
!
! Larmor radius:
  rlarm=v0*n_d*p_mass*c/(n_e*e_charge*bmod_ref)
! normalized slowing down time:
  tau=trace_time*v0
! normalized time step:
  dtau=tau/dfloat(ntimstep-1)
!
bmod00=281679.46317784750d0
! Larmor raidus corresponds to the field stregth egual to $B_{00}$ harmonic
! in Boozer coordinates:
! 14.11.2011  bmod00=bmod_ref  !<=deactivated, use value from the 'alpha_lifetime.inp'
  ro0=rlarm*bmod00  ! 23.09.2013
!
  multharm=3 !7
  ns_A=5
  ns_s=5
  ns_tp=5
!
  call spline_vmec_data
!
  call stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0)         !<=2017
!
  rbig=rt0
! field line integration step step over phi (to check chamber wall crossing)
  dphi=2.d0*pi/(L1i*npoiper)
! orbit integration time step (to check chamber wall crossing)
  dtaumin=dphi*rbig/npoiper2!
dtau=2*dtaumin
print *,dtau
!
  call get_canonical_coordinates
! call testing
!
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
  print *,'can : ',r,vartheta_c,varphi_c
!
  isw_field_type=0
  z(1)=r
  z(2)=vartheta_c
  z(3)=varphi_c
  z(4)=1.d0
  z(5)=alam0
!
print *,'symplectic'
  call orbit_sympl_init(z) 
  do i=1,L1i*npoiper*npoiper2*10
!
    call orbit_timestep_sympl(z,dtau,dtaumin,ierr)
!
    call can_to_vmec(z(1),z(2),z(3),theta_vmec,varphi_vmec)
!
!
    write (3004,*) dtau*dfloat(i),z(1),theta_vmec,varphi_vmec,z(4:5),z(2:3)
  enddo
print *,'done. Evaluations: ', neval
!
  call can_to_vmec(r,vartheta_c,varphi_c,theta_vmec,varphi_vmec)
!
!  call deallocate_can_coord
!
  call spline_vmec_data
!
  isw_field_type=1
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
!call testing
!
  end
!
