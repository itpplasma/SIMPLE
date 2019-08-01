!
  use new_vmec_stuff_mod, only : netcdffile,multharm,ns_A,ns_s,ns_tp
!  use chamb_mod,  only : rbig,rcham2
  use parmot_mod, only : rmu,ro0,eeff
  use velo_mod,   only : isw_field_type
use diag_mod, only : icounter
  use orbit_symplectic, only : SymplecticIntegrator, orbit_sympl_init, orbit_timestep_sympl
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
  double precision :: rbig,z1,z2
  double precision, dimension(5) :: z
  integer          :: npoiper2
  double precision :: contr_pp
  double precision :: facE_al
  integer          :: ibins
  integer          :: n_e,n_d,n_b
  double precision :: r,vartheta_c,varphi_c,theta_vmec,varphi_vmec,alam0

  integer, parameter :: mode_sympl = 0 ! 0 = Euler1, 1 = Euler2, 2 = Verlet
!
!---------------------------------------------------------------------------
! Prepare calculation of orbit tip by interpolation
!
  integer                                       :: nplagr,nder,itip,npl_half
  double precision                              :: alam_prev,zerolam,twopi,fraction
  double precision, dimension(5)                :: z_tip
  integer,          dimension(:),   allocatable :: ipoi
  double precision, dimension(:),   allocatable :: xp
  double precision, dimension(:,:), allocatable :: coef,orb_sten
!
  type(SymplecticIntegrator) :: si
  
  zerolam=0.d0
  twopi=2.d0*pi
  nplagr=4
  nder=0
  npl_half=nplagr/2
  allocate(ipoi(nplagr),coef(0:nder,nplagr),orb_sten(5,nplagr),xp(nplagr))
  do i=1,nplagr
    ipoi(i)=i
  enddo
!
! End prepare calculation of orbit tip by interpolation
!--------------------------------------------------------------------------
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
  dtau=tau/dble(ntimstep-1)
!
bmod00=281679.46317784750d0
! Larmor raidus corresponds to the field stregth egual to $B_{00}$ harmonic
! in Boozer coordinates:
! 14.11.2011  bmod00=bmod_ref  !<=deactivated, use value from the 'alpha_lifetime.inp'
  ro0=rlarm*bmod00  ! 23.09.2013
!
  multharm=3 !3 !7
  ns_A=5
  ns_s=5
  ns_tp=5
!
  call spline_vmec_data
!call testing
!
  call stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0)         !<=2017
!
  rbig=rt0
! field line integration step step over phi (to check chamber wall crossing)
  dphi=2.d0*pi/(L1i*npoiper)
! orbit integration time step (to check chamber wall crossing)
  dtaumin=dphi*rbig/npoiper2!
!dtau=2*dtaumin
dtau=dtaumin
ntimstep = L1i*npoiper*npoiper2*10000
print *, 'dtau = ', dtau, ' dtau/dtaumin = ', dtau/dtaumin
print *, 'ttrace = ', ntimstep*dtau/v0, 'nstep = ', ntimstep
!
  call get_canonical_coordinates
!call testing
!
!do 
  !print *, 'Enter r, theta, phi, lambda: '
  !read *,r,vartheta_c,varphi_c,alam0
  r = 0.5
  vartheta_c = 0.0
  varphi_c = 0.314
  alam0 = 0.0
!
  isw_field_type=0
  z(1)=r
  z(2)=vartheta_c
  z(3)=varphi_c
  z(4)=1.d0
  z(5)=alam0
!
icounter=0
  call orbit_sympl_init(si, z, dtau, dtaumin, 1d-12, mode_sympl) 
!
!--------------------------------
! Initialize tip detector
!
  itip=3
  alam_prev=z(5)
!
! End initialize tip detector
!--------------------------------
!
  open(101,file='poiplot.dat')
!
  do i=1,ntimstep !300 !10
!
!    call orbit_timestep_axis(z,dtau,dtaumin,ierr)
    call orbit_timestep_sympl(si, z,ierr)
!
    if(ierr.ne.0) exit
!
!-------------------------------------------------------------------------
! Tip detection and interpolation
!
    if(alam_prev.lt.0.d0.and.z(5).gt.0.d0) itip=0   !<=tip has been passed
    itip=itip+1
    alam_prev=z(5)
    if(i.le.nplagr) then          !<=first nplagr points to initialize stencil
      orb_sten(:,i)=z
    else                          !<=normal case, shift stencil
      orb_sten(:,ipoi(1))=z
      ipoi=cshift(ipoi,1)
      if(itip.eq.npl_half) then   !<=stencil around tip is complete, interpolate
        xp=orb_sten(5,ipoi)
!
        call plag_coeff(nplagr,nder,zerolam,xp,coef)
!
        z_tip=matmul(orb_sten(:,ipoi),coef(0,:))
        z_tip(2)=modulo(z_tip(2),twopi)
        z_tip(3)=modulo(z_tip(3),twopi)
        write(101,*) z_tip
      endif
    endif
!
! End tip detection and interpolation
!-------------------------------------------------------------------------
!
  enddo
  close(101)
!
print *,'done  ',icounter,'  field calls', icounter*1.0d0/ntimstep, 'per step'
!
  call fract_dimension(fraction)
!
  if(fraction.gt.0.3d0) then
    print *,'chaotic orbit'
  else
    print *,'regular orbit'
  endif
!enddo
!
  call deallocate_can_coord
!
  end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine fract_dimension(fraction)
!
  implicit none
!
  integer, parameter :: iunit=171
  integer :: itr,ntr,ir,it,ngrid,nrefine,irefine,kr,kt,nboxes
  double precision :: fraction,r,rmax,rmin,tmax,tmin,hr,ht
  logical,          dimension(:,:), allocatable :: free
  double precision, dimension(:,:), allocatable :: rt
!
  ntr=0
!
  open(iunit,file='poiplot.dat')
  do
    read(iunit,*,end=1) r
    ntr=ntr+1
  enddo
1 close(iunit)
!
  allocate(rt(2,ntr))
  open(iunit,file='poiplot.dat')
  do itr=1,ntr
    read(iunit,*) rt(:,itr)
  enddo
  close(iunit)
!
  rmin=minval(rt(1,:))
  rmax=maxval(rt(1,:))
  tmin=minval(rt(2,:))
  tmax=maxval(rt(2,:))
!
  nrefine=int(log(dble(ntr))/log(4.d0))
!
  open(iunit,file='boxcount.dat')
  ngrid=1
  nrefine=nrefine+3       !<=add 3 for curiousity
  do irefine=1,nrefine
    ngrid=ngrid*2
    allocate(free(0:ngrid,0:ngrid))
    free=.true.
    hr=(rmax-rmin)/dble(ngrid)
    ht=(tmax-tmin)/dble(ngrid)
    nboxes=0
    do itr=1,ntr
      kr=int((rt(1,itr)-rmin)/hr)
      kr=min(ngrid-1,max(0,kr))
      kt=int((rt(2,itr)-tmin)/ht)
      kt=min(ngrid-1,max(0,kt))
      if(free(kr,kt)) then
        free(kr,kt)=.false.
        nboxes=nboxes+1
      endif
    enddo
    deallocate(free)
    write(iunit,*) dble(irefine),dble(nboxes)/dble(ngrid**2)
    if(irefine.eq.nrefine-3) fraction=dble(nboxes)/dble(ngrid**2)
  enddo
  close(iunit)
  deallocate(rt)
!
  end subroutine fract_dimension
