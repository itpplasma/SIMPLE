!
  use new_vmec_stuff_mod, only : netcdffile,multharm,ns_A,ns_s,ns_tp
!  use chamb_mod,  only : rbig,rcham2
  use parmot_mod, only : rmu,ro0,eeff
  use velo_mod,   only : isw_field_type
use diag_mod, only : icounter
  use orbit_symplectic, only : orbit_sympl_init, orbit_timestep_sympl
  use common, only: pi,c,e_charge,e_mass,p_mass,ev
!
#ifdef _OPENMP
  use omp_lib
#endif

!
  implicit none
!

  double precision,parameter  :: snear_axis=0.05d0
!
  logical          :: near_axis
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
  integer, parameter :: npart = 984 !960
  double precision :: r,vartheta_c(npart),varphi_c(npart),theta_vmec,varphi_vmec,alam0(npart)
!
  integer :: i_ctr ! for nice counting in parallel
!
!---------------------------------------------------------------------------  
! buffer for Poincare plot:
!
!---------------------------------------------------------------------------
! Prepare calculation of orbit tip by interpolation and buffer for Poincare plot:
!
  integer                                       :: nplagr,nder,npl_half,itip,iper
  integer                                       :: npass_regular,npass_chaotic,ntr_regular,ntr_chaotic
  integer                                       :: nstep_tot,norbper,ifp_tip,ifp_per
  integer                                       :: nfp,nfp_tip,nfp_per,kper
  double precision                              :: alam_prev,zerolam,twopi,fraction,fper,phiper
  double precision, dimension(5)                :: z_tip
  integer,          dimension(:),   allocatable :: ipoi
  double precision, dimension(:),   allocatable :: xp
  double precision, dimension(:,:), allocatable :: coef,orb_sten
  double precision, dimension(:,:), allocatable :: zpoipl_tip,zpoipl_per,dummy2d

  zerolam=0.d0
  twopi=2.d0*pi
  nplagr=6
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
print *,dtau
!
  fper=twopi/dble(L1i)   !<= field period
!
  call get_canonical_coordinates
!
  npass_regular=0
  npass_chaotic=0
  ntr_regular=0
  ntr_chaotic=0
!
print *, 'generating random initial conditions'
r=0.5d0
do ipart=1,npart
!
  call random_number(zzg)
!
  vartheta_c(ipart)=twopi*zzg
!
  call random_number(zzg)
!
  varphi_c(ipart)=twopi*zzg
!
  call random_number(zzg)
!
  alam0(ipart)=2.d0*zzg-1.d0
enddo

isw_field_type=0
i_ctr=0
!$omp parallel private(z,ifp,alam_prev,itip,ierr,orb_sten,xp,z_tip,i) &
!$omp& firstprivate(ipoi)
!print *, 'run started on thread ', omp_get_thread_num()

!$omp do
do ipart=1,npart  
!print *, 'particle ', ipart, '/', npart
!
  z(1)=r
  z(2)=vartheta_c(ipart)
  z(3)=varphi_c(ipart)
  z(4)=1.d0
  z(5)=alam0(ipart)
!  
!print *, 'z=', z
!
  norbper=10000 !300 !10
  nstep_tot=L1i*npoiper*norbper
  nfp=L1i*norbper         !<= guess for footprint number
  nfp_tip=nfp             !<= initial array dimension for tips
  nfp_per=nfp             !<= initial array dimension for periods
!
  allocate(zpoipl_tip(2,nfp_tip),zpoipl_per(2,nfp_per))
!
  ifp_tip=0               !<= initialize footprint counter on tips
  ifp_per=0               !<= initialize footprint counter on periods
!
icounter=0
  call orbit_sympl_init(z, dtau, dtaumin, 0) 
!
!--------------------------------
! Initialize tip detector
!
  itip=npl_half+1
  alam_prev=z(5)
!
! End initialize tip detector
!--------------------------------
! Initialize period crossing detector
!
  iper=npl_half+1
  kper=int(z(3)/fper)
!
! End initialize period crossing detector
!--------------------------------
!
!
  do i=1,nstep_tot
!
!    call orbit_timestep_axis(z,dtau,dtaumin,ierr)    
    call orbit_timestep_sympl(z, ierr)
!
    if(ierr.ne.0) exit
!
    if(i.le.nplagr) then          !<=first nplagr points to initialize stencil
      orb_sten(:,i)=z
    else                          !<=normal case, shift stencil
      orb_sten(:,ipoi(1))=z
      ipoi=cshift(ipoi,1)
    endif
!
!-------------------------------------------------------------------------
! Tip detection and interpolation
!
    if(alam_prev.lt.0.d0.and.z(5).gt.0.d0) itip=0   !<=tip has been passed
    itip=itip+1
    alam_prev=z(5)
    if(i.gt.nplagr) then          !<=use only initialized stencil
      if(itip.eq.npl_half) then   !<=stencil around tip is complete, interpolate
        xp=orb_sten(5,ipoi)
!
        call plag_coeff(nplagr,nder,zerolam,xp,coef)
!
        z_tip=matmul(orb_sten(:,ipoi),coef(0,:))
        z_tip(2)=modulo(z_tip(2),twopi)
        z_tip(3)=modulo(z_tip(3),twopi)
write(1001,*) z_tip
        ifp_tip=ifp_tip+1
        if(ifp_tip.gt.nfp_tip) then   !<=increase the buffer for banana tips
          allocate(dummy2d(2,ifp_tip-1))
          dummy2d=zpoipl_tip(:,1:ifp_tip-1)
          deallocate(zpoipl_tip)
          nfp_tip=nfp_tip+nfp
          allocate(zpoipl_tip(2,nfp_tip))
          zpoipl_tip(:,1:ifp_tip-1)=dummy2d
          deallocate(dummy2d)
        endif
        zpoipl_tip(:,ifp_tip)=z_tip(1:2)
      endif
    endif
!
! End tip detection and interpolation
!-------------------------------------------------------------------------
! Periodic boundary footprint detection and interpolation
!
    if(z(3).gt.dble(kper+1)*fper) then
      iper=0   !<=periodic boundary has been passed
      phiper=dble(kper+1)*fper
      kper=kper+1
    elseif(z(3).lt.dble(kper)*fper) then
      iper=0   !<=periodic boundary has been passed
      phiper=dble(kper)*fper
      kper=kper-1
    endif
    iper=iper+1
    if(i.gt.nplagr) then          !<=use only initialized stencil 
      if(iper.eq.npl_half) then   !<=stencil around periodic boundary is complete, interpolate
        xp=orb_sten(3,ipoi)-phiper
!
        call plag_coeff(nplagr,nder,zerolam,xp,coef)
!
        z_tip=matmul(orb_sten(:,ipoi),coef(0,:))
        z_tip(2)=modulo(z_tip(2),twopi)
        z_tip(3)=modulo(z_tip(3),twopi)
write(1002,*) z_tip
        ifp_per=ifp_per+1
        if(ifp_per.gt.nfp_per) then   !<=increase the buffer for periodic boundary footprints
          allocate(dummy2d(2,ifp_per-1))
          dummy2d=zpoipl_per(:,1:ifp_per-1)
          deallocate(zpoipl_per)
          nfp_per=nfp_per+nfp
          allocate(zpoipl_per(2,nfp_per))
          zpoipl_per(:,1:ifp_per-1)=dummy2d
          deallocate(dummy2d)
        endif
        zpoipl_per(:,ifp_per)=z_tip(1:2)
      endif
    endif
!
! End periodic boundary footprint detection and interpolation
!-------------------------------------------------------------------------
!
  enddo
close(1001)
close(1002)
!
!$omp critical
i_ctr = i_ctr+1
print *, 'particle ', i_ctr, '/', npart, ' done: ',icounter,'  field calls'
!$omp end critical
!
  if(ifp_tip.eq.0) then
!
    call fract_dimension(ifp_per,zpoipl_per(:,1:ifp_per),fraction)
!
    if(fraction.gt.0.2d0) then
      print *,'chaotic passing orbit'
!$omp atomic
      npass_chaotic=npass_chaotic+1
    else
      print *,'regular passing orbit'
!$omp atomic
      npass_regular=npass_regular+1
    endif
  else
!
    call fract_dimension(ifp_tip,zpoipl_tip(:,1:ifp_tip),fraction)
!
    if(fraction.gt.0.2d0) then
      print *,'chaotic trapped orbit'
!$omp atomic
      ntr_chaotic=ntr_chaotic+1
    else
      print *,'regular trapped orbit'
!$omp atomic
      ntr_regular=ntr_regular+1
    endif
  endif
!
  deallocate(zpoipl_tip,zpoipl_per)
print *,npass_regular,' passing regular ',npass_chaotic,' passing chaotic ', &
        ntr_regular,' trapped regular ',ntr_chaotic,' trapped chaotic'
enddo
!$omp end do
!$omp end parallel
!
  call deallocate_can_coord
!
  end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine fract_dimension(ntr,rt,fraction)
!
#ifdef _OPENMP
  use omp_lib
#endif
!
  implicit none
!
integer, parameter :: iunit=1003
  integer :: itr,ntr,ir,it,ngrid,nrefine,irefine,kr,kt,nboxes
  double precision :: fraction,r,rmax,rmin,tmax,tmin,hr,ht
  double precision, dimension(2,ntr)              :: rt
  logical,          dimension(:,:),   allocatable :: free
!
  rmin=minval(rt(1,:))
  rmax=maxval(rt(1,:))
  tmin=minval(rt(2,:))
  tmax=maxval(rt(2,:))
!
  nrefine=int(log(dble(ntr))/log(4.d0))
!
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
!
  end subroutine fract_dimension
