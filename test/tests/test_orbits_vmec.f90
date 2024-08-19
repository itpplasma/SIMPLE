!
  use new_vmec_stuff_mod, only : netcdffile,multharm,ns_A,ns_s,ns_tp
!  use chamb_mod,  only : rbig,rcham2
  use parmot_mod, only : rmu,ro0,eeff
  use velo_mod,   only : isw_field_type
  use boozer_coordinates_mod, only : use_B_r
  use get_can_sub
  use boozer_sub
  use alpha_lifetime_sub
  use spline_vmec_sub
use diag_mod, only : icounter
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
  double precision,parameter  :: c=2.9979d10
  double precision,parameter  :: e_charge=4.8032d-10
  double precision,parameter  :: e_mass=9.1094d-28
  double precision,parameter  :: p_mass=1.6726d-24
  double precision,parameter  :: ev=1.6022d-12
  double precision,parameter  :: snear_axis=0.05d0
!
  logical :: near_axis
  integer          :: npoi,ierr,L1i,nper,npoiper,i,ntimstep,ntestpart
  integer          :: ipart,notrace_passing,loopskip,iskip,ilost,it
  real             :: zzg
  double precision :: dphi,rbeg,phibeg,zbeg,bmod00,rcham,rlarm,bmax,bmin
  double precision :: tau,dtau,dtaumin,xi,v0,bmod_ref,E_alpha,trace_time
  double precision :: RT0,R0i,cbfi,bz0i,bf0,trap_par
  double precision :: sbeg,thetabeg
  double precision :: z1,z2
  double precision, dimension(5) :: z
  integer          :: npoiper2
  double precision :: contr_pp
  double precision :: facE_al
  integer          :: ibins
  integer          :: n_e,n_d,n_b
  double precision :: r,vartheta_c,varphi_c,theta_vmec,varphi_vmec,alam0
  double precision :: vartheta,varphi
  double precision :: alam,alam_prev,par_inv
  real :: tstart, tend
  integer, parameter :: runlen = 1
  logical, parameter :: jparmode = .false. ! if true output orbit points only when vpar changes sign from positive to negative
  integer, parameter :: mode_sympl = 0 ! 0=NEO-ORB in vmec (sergei); 1=symplectic in canonical (christopher; euler1)
  integer :: ntau

  open(1,file='alpha_lifetime_m.inp', recl=1024)
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

isw_field_type=2
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
  dtau=tau/dble(ntimstep-1)
!
bmod00=281679.46317784750d0
! Larmor raidus corresponds to the field stregth egual to $B_{00}$ harmonic
! in Boozer coordinates:
! 14.11.2011  bmod00=bmod_ref  !<=deactivated, use value from the 'alpha_lifetime.inp'
  ro0=rlarm*bmod00  ! 23.09.2013
!
  multharm=7 !3 ! resolution of interpolation in angles
  ns_A=5
  ns_s=5
  ns_tp=5
!
call spline_vmec_data
!call testing
!
  call stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0)         !<=2017
! field line integration step step over phi (to check chamber wall crossing)
  dphi=2.d0*pi/(L1i*npoiper)
! orbit integration time step (to check chamber wall crossing)
  dtaumin=dphi*rt0/npoiper2!
!dtau=2*dtaumin
!ntau = ceiling(dtau/dtaumin)
!dtaumin = dtau/ntau
dtau = dtaumin
ntimstep = L1i*npoiper*npoiper2*100
print *, 'dtau = ', dtau, ' dtau/dtaumin = ', dtau/dtaumin, 'tau = ', tau
!
! initial conditions
r = 0.35
theta_vmec = 0.33
varphi_vmec = 0.97
alam0 = 0.1

if (isw_field_type == 0) then
    call get_canonical_coordinates
    call vmec_to_can(r,theta_vmec,varphi_vmec,vartheta,varphi)
elseif (isw_field_type == 2) then
    use_B_r=.false.
!    use_B_r=.true.
    call boozer_converter
    call vmec_to_boozer(r,theta_vmec,varphi_vmec,vartheta,varphi)
else
    vartheta=theta_vmec
    varphi=varphi_vmec
endif

! !
z(1)=r           ! radius ~psi_pol
z(2)=vartheta    ! poloidal angle
z(3)=varphi      ! toroidal angle
z(4)=1.d0        ! normalized momentum
z(5)=alam0       ! parameter for pitch angle between parallel and perp velocity
! !

par_inv=0.d0
alam=alam0
alam_prev=alam

print *,'VMEC, splines'
call cpu_time(tstart)
open(3001, file='orbit_rk.out', recl=1024)
do i=1,L1i*npoiper*npoiper2*runlen
! !
    if (mode_sympl==0) call orbit_timestep_axis(z,dtau,dtaumin,1d-10,ierr)
    r=z(1)
    vartheta=z(2)
    varphi=z(3)
    if (isw_field_type == 0) then
      call can_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
    elseif (isw_field_type == 1) then
      theta_vmec=vartheta
      varphi_vmec=varphi
    elseif (isw_field_type == 2) then
      call boozer_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
    endif

    if (.not. jparmode) then
      write (3001,*) dtau*dble(i),z,theta_vmec,varphi_vmec
    else
      alam=z(5)
      par_inv=par_inv+alam**2*dtau
      if(alam_prev.lt.0.d0.and.alam.gt.0.d0) then
        write (100,*) i,par_inv
        write (3001,*) dtau*dble(i),z,theta_vmec,varphi_vmec
        par_inv=0.d0
      endif
      alam_prev=alam
    endif
enddo
close(3001)
call cpu_time(tend)
print *, 'CPU time (s): ', tend - tstart
end
