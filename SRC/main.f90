program neo_orb_main
  use omp_lib
  use util, only: pi, twopi, c, e_charge, e_mass, p_mass, ev, sqrt2
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp, &
    vmec_B_scale, vmec_RZ_scale

  use velo_mod,   only : isw_field_type
  use orbit_symplectic, only : orbit_timestep_sympl, get_val
  use simple, only : init_field, init_sympl, Tracer, debug, eval_field
  use cut_detector, only : fract_dimension
  use diag_mod, only : icounter
  use collis_alp, only : loacol_alpha, stost
  use params

  implicit none

  double precision :: facE_al, bmod_ref, trace_time
  integer :: ntimstep, npoiper, npoiper2, n_e, n_d
  type(Tracer) :: norb

! read config file
  call read_config

! initialize field geometry
  call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
  call params_init(3.5d6/facE_al, bmod_ref, trace_time, ntimstep, npoiper, &
    npoiper2, n_e, n_d)
  print *, 'tau: ', dtau, dtaumin, min(dabs(mod(dtau, dtaumin)), &
                    dabs(mod(dtau, dtaumin)-dtaumin))/dtaumin, ntau
  print *, 'v0 = ', v0


! init collisions
  if (swcoll) then
    call loacol_alpha(am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe, &
      3.5d6/facE_al,v0,dchichi,slowrate,dchichi_norm,slowrate_norm)
  endif
  print *, 'v0 = ', v0

! pre-compute starting flux surface
  npoi=nper*npoiper ! total number of starting points
  allocate(xstart(3,npoi),bstart(npoi),volstart(npoi))
  call init_starting_surf

! initialize array of confined particle percentage
  allocate(confpart_trap(ntimstep),confpart_pass(ntimstep))
  confpart_trap=0.d0
  confpart_pass=0.d0

! initialize lost times when particles get lost
  allocate(times_lost(ntestpart), trap_par(ntestpart), perp_inv(ntestpart))
  allocate(iclass(3,ntestpart))
  times_lost = -1.d0

  allocate(zstart(5,ntestpart), zend(5,ntestpart))
  if(local) then
    call init_starting_points
  else
    call init_starting_points_global
  endif
  if (startmode == 0) stop

  icounter=0 ! evaluation counter

! do particle tracing in parallel

  !$omp parallel firstprivate(norb)
  !$omp do
  do i=1,ntestpart
    !$omp critical
    kpart = kpart+1
    print *, kpart, ' / ', ntestpart, 'particle: ', i, 'thread: ', omp_get_thread_num()
    !$omp end critical
    call trace_orbit(norb, i)
  end do
  !$omp end do
  !$omp end parallel

  confpart_pass=confpart_pass/ntestpart
  confpart_trap=confpart_trap/ntestpart

  open(1,file='confined_fraction.dat',recl=1024)
  do i=1,ntimstep
    write(1,*) dble(i-1)*dtau/v0,confpart_pass(i),confpart_trap(i),ntestpart
  enddo
  close(1)

  open(1,file='times_lost.dat',recl=1024)
  do i=1,ntestpart
    write(1,*) i, times_lost(i), trap_par(i), zstart(1,i), perp_inv(i), zend(:,i)
  enddo
  close(1)

  open(1,file='class_parts.dat',recl=1024)
  do i=1,ntestpart
    write(1,*) i, zstart(1,i), perp_inv(i), iclass(:,i)
  enddo
  close(1)

  if (integmode >= 0) call deallocate_can_coord

  deallocate(times_lost, confpart_trap, confpart_pass, trap_par, perp_inv, iclass)
  deallocate(xstart, bstart, volstart, zstart)

contains

subroutine read_config
  open(1,file='simple.in',recl=1024)
  read (1,*) notrace_passing   !skip tracing passing prts if notrace_passing=1
  read (1,*) nper              !number of periods for initial field line        ! TODO: increase
  read (1,*) npoiper           !number of points per period on this field line  ! TODO: increase
  read (1,*) ntimstep          !number of time steps per slowing down time
  read (1,*) ntestpart         !number of test particles
  read (1,*) bmod_ref          !reference field, G, for Boozer $B_{00}$
  read (1,*) trace_time        !slowing down time, s
  read (1,*) sbeg              !starting s for field line                       !<=2017
  read (1,*) phibeg            !starting phi for field line                     !<=2017
  read (1,*) thetabeg          !starting theta for field line                   !<=2017
  read (1,*) loopskip          !how many loops to skip to shift random numbers
  read (1,*) contr_pp          !control of passing particle fraction            ! UNUSED (2019)
  read (1,*) facE_al           !facE_al test particle energy reduction factor
  read (1,*) npoiper2          !additional integration step split factor
  read (1,*) n_e               !test particle charge number (the same as Z)
  read (1,*) n_d               !test particle mass number (the same as A)
  read (1,*) netcdffile        !name of VMEC file in NETCDF format <=2017 NEW
  read (1,*) ns_s              !spline order for 3D quantities over s variable
  read (1,*) ns_tp             !spline order for 3D quantities over theta and phi
  read (1,*) multharm          !angular grid factor (n_grid=multharm*n_harm_max where n_harm_max - maximum Fourier index)
  read (1,*) isw_field_type    !field type: -1 - Testing, 0 - Canonical, 1 - VMEC, 2 - Boozer
  read (1,*) startmode         !mode for initial conditions: 0=generate and store, 1=generate, store, and run, 2=read and run, 3=read ANTS and run
  read (1,*) integmode         !mode for integrator: -1 = RK VMEC, 0 = RK CAN, 1 = Euler1, 2 = Euler2, 3 = Verlet
  read (1,*) relerr            !relative error for RK integrator
  read (1,*) tcut              !time when to do cut for classification, usually 1d-1, or -1 if no cuts desired
  read (1,*) debug             !produce debugging output (.True./.False.). Use only in non-parallel mode!
  read (1,*) class_plot        !write starting points at phi=const cut for classification plot (.True./.False.).  !<=AAA
  read (1,*) cut_in_per        !normalized phi-cut position within field period, [0:1], used if class_plot=.True. !<=AAA
  read (1,*) fast_class        !if .True. quit immeadiately after fast classification and don't trace orbits to the end
  read (1,*) local             !if .True. orbits are started on a single flux surface, if .False. (not yet equally distributed) in volume
  read (1,*) vmec_B_scale      !factor to scale the B field from VMEC
  read (1,*) vmec_RZ_scale     !factor to scale the device size from VMEC
  read (1,*) swcoll            !if .True. enables collisions. This is incompatible with classification.
  read (1,*) deterministic     !if .True. put seed for the same random walk always
! TODO: add new config options for collisions
!  read (1,*) am1           !mass number of the 1-st field ion species
!  read (1,*) am2           !mass number of the 2-nd field ion species
!  read (1,*) Z1            !charge number of the 1-st field ion species
!  read (1,*) Z2            !charge number of the 2-nd field ion species
!  read (1,*) densi1        !density of the 1-st field ion species, 1/cm^3
!  read (1,*) densi2        !density of the 2-nd field ion species, 1/cm^3
!  read (1,*) tempi1        !temperature of the 1-st field ion species, eV
!  read (1,*) tempi2        !temperature of the 2-nd field ion species, eV
!  read (1,*) tempe         !temperature of electrons, eV
  close(1)

  if (swcoll .and. (tcut > 0.0d0 .or. class_plot .or. fast_class)) then
    stop 'Collisions are incompatible with classification'
  endif

end subroutine read_config

subroutine init_starting_surf
  integer :: ierr

  xstart=0.d0
  bstart=0.d0
  volstart=0.d0

  call integrate_mfl_can( &
    npoi,dphi,sbeg,phibeg,thetabeg, &
    xstart,bstart,volstart,bmod00,ierr)

  if(ierr.ne.0) then
    print *,'starting field line has points outside the chamber'
    stop
  endif

! Larmor radius corresponds to the field stregth egual to $B_{00}$ harmonic
! in Boozer coordinates:
  ro0=rlarm*bmod00
! maximum value of B module:
  bmax=maxval(bstart)
  bmin=minval(bstart)

  print *, 'bmod00 = ', bmod00, 'bmin = ', bmin, 'bmax = ', bmax
end subroutine init_starting_surf

subroutine init_starting_points_ants(unit)
  use parse_ants, only: process_line
  integer, intent(in) :: unit

  integer, parameter :: maxlen = 4096
  character(len=maxlen) :: line
  real(8) :: v_par, v_perp, u, v, s
  real(8) :: th, ph, th_c, ph_c  ! Canonical flux coordinate angles
  integer :: ipart

  do ipart=1,ntestpart
    read(unit, '(A)') line
    call process_line(line, v_par, v_perp, u, v, s)
    ! In the test case, u runs from 0 to 1 and v from 0 to 4
    th = 2d0*pi*u
    ph = 2d0*pi*v/4d0
    call vmec_to_can(s, th, ph, th_c, ph_c)
    zstart(1, ipart) = s
    zstart(2, ipart) = ph_c
    zstart(3, ipart) = th_c
    zstart(4, ipart) = 1.d0
    zstart(5, ipart) = v_par / sqrt(v_par**2 + v_perp**2)
  enddo
end subroutine

subroutine init_starting_points
  integer :: ipart
  real :: zzg
  double precision :: r,vartheta,varphi,theta_vmec,varphi_vmec

  ! skip random numbers according to configuration
    do iskip=1,loopskip
      do ipart=1,ntestpart
        xi=zzg()
        xi=zzg()
      enddo
    enddo

  ! files for storing starting coords
  open(1,file='start.dat',recl=1024)
  ! determine the starting point:
  if (startmode == 0 .or. startmode == 1) then
    do ipart=1,ntestpart
      xi=zzg()
      call binsrc(volstart,1,npoi,xi,i)
      ibins=i
      ! coordinates: z(1) = r, z(2) = vartheta, z(3) = varphi
      r=xstart(1,i)
      vartheta=xstart(2,i)
      varphi=xstart(3,i)
!
! we store starting points in VMEC coordinates:
      if(isw_field_type.eq.0) then
        call can_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
      elseif(isw_field_type.eq.1) then
        theta_vmec=vartheta
        varphi_vmec=varphi
      elseif(isw_field_type.eq.2) then
        call boozer_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
      else
        print *,'init_starting_points: unknown field type'
      endif
!
      zstart(1,ipart)=r
      zstart(2,ipart)=theta_vmec
      zstart(3,ipart)=varphi_vmec
      ! normalized velocity module z(4) = v / v_0:
      zstart(4,ipart)=1.d0
      ! starting pitch z(5)=v_\parallel / v:
      xi=zzg()
      zstart(5,ipart)=2.d0*(xi-0.5d0)
      write(1,*) zstart(:,ipart)
    enddo
  else if (startmode == 2) then
    do ipart=1,ntestpart
      read(1,*) zstart(:,ipart)
    enddo
  else if (startmode == 3) then  ! ANTS input mode
    call init_starting_points_ants(1)
  endif

  close(1)
end subroutine init_starting_points

subroutine init_starting_points_global
  integer, parameter :: ns=1000
  integer :: ipart,is
  real :: zzg
  double precision :: r,vartheta,varphi,theta_vmec,varphi_vmec
  double precision :: s,bmin,bmax
!
  open(1,file='bminmax.dat',recl=1024)
  do is=0,ns
    s=dble(is)/dble(ns)
!
    call get_bminmax(s,bmin,bmax)
!
    write(1,*) s,bmin,bmax
  enddo
  close(1)

  ! skip random numbers according to configuration
    do iskip=1,loopskip
      do ipart=1,ntestpart
        xi=zzg()
        xi=zzg()
        xi=zzg()
        xi=zzg()
      enddo
    enddo

  ! files for storing starting coords
  open(1,file='start.dat',recl=1024)
  ! determine the starting point:
  if (startmode == 0 .or. startmode == 1) then
    do ipart=1,ntestpart
      xi=zzg()
      r=xi
      xi=zzg()
      vartheta=twopi*xi
      xi=zzg()
      varphi=twopi*xi
!
! we store starting points in VMEC coordinates:
      if(isw_field_type.eq.0) then
        call can_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
      elseif(isw_field_type.eq.1) then
        theta_vmec=vartheta
        varphi_vmec=varphi
      elseif(isw_field_type.eq.2) then
        call boozer_to_vmec(r,vartheta,varphi,theta_vmec,varphi_vmec)
      else
        print *,'init_starting_points: unknown field type'
      endif
!
      zstart(1,ipart)=r
      zstart(2,ipart)=theta_vmec
      zstart(3,ipart)=varphi_vmec
      ! normalized velocity module z(4) = v / v_0:
      zstart(4,ipart)=1.d0
      ! starting pitch z(5)=v_\parallel / v:
      xi=zzg()
      zstart(5,ipart)=2.d0*(xi-0.5d0)
      write(1,*) zstart(:,ipart)
    enddo
  else if (startmode == 2) then
    do ipart=1,ntestpart
      read(1,*) zstart(:,ipart)
    enddo
  else if (startmode == 3) then  ! ANTS input mode
    call init_starting_points_ants(1)
  endif

  close(1)
end subroutine init_starting_points_global

subroutine trace_orbit(anorb, ipart)
  type(Tracer), intent(inout) :: anorb
  integer, intent(in) :: ipart
  integer :: ierr, ierr_coll
  double precision, dimension(5) :: z
  double precision :: bmod,sqrtg
  double precision, dimension(3) :: bder, hcovar, hctrvr, hcurl
  integer :: it, ktau
  integer(8) :: kt
  logical :: passing

  integer                                       :: ifp_tip,ifp_per
  integer,          dimension(:),   allocatable :: ipoi
  double precision, dimension(:),   allocatable :: xp
  double precision, dimension(:,:), allocatable :: coef,orb_sten
  double precision, dimension(:,:), allocatable :: zpoipl_tip,zpoipl_per,dummy2d
  double precision, dimension(n_tip_vars)       :: var_tip
  double precision :: phiper, alam_prev, par_inv
  integer :: iper, itip, kper, nfp_tip, nfp_per

  double precision :: fraction
  double precision :: r,theta_vmec,varphi_vmec
  logical :: regular
! output files:
! iaaa_bou - trapped-passing boundary
! iaaa_pnt - forced regular passing
! iaaa_prp - lossed passing
! iaaa_prt - lossed trapped
! iaaa_rep - regular passing
! iaaa_ret - regular trapped
! iaaa_stp - stochastic passing
! iaaa_stt - stochastic trapped
  integer, parameter :: iaaa_bou=20000, iaaa_pnt=10000, iaaa_prp=10001, iaaa_prt=10002, &  !<=AAA
                        iaaa_rep=10011, iaaa_ret=10012, iaaa_stp=10021, iaaa_stt=10022     !<=AAA

! Variables and settings for classification by J_parallel and ideal orbit condition:
  integer, parameter :: nfp_dim=3, nturns=8
  integer :: nfp_cot,ideal,ijpar,ierr_cot,iangvar
  double precision, dimension(nfp_dim) :: fpr_in
! output files:
! iaaa_jre - regular trapped by J_parallel
! iaaa_jst - stochastic trapped by J_parallel
! iaaa_jer - non-classified trapped by J_parallel
! iaaa_ire - ideal trapped by recurrences and monotonicity
! iaaa_ist - non-ideal trapped by recurrences and monotonicity
! iaaa_ier - non-classified trapped by recurrences and monotonicity
  integer, parameter :: iaaa_jre=40012, iaaa_jst=40022, iaaa_jer=40032, &
                        iaaa_ire=50012, iaaa_ist=50022, iaaa_ier=50032

! for run with fixed random seed
  integer :: seedsize
  integer, allocatable :: seed(:)

  if (deterministic) then
    call random_seed(size = seedsize)
    allocate(seed(seedsize))
    seed = 0
    call random_seed(put=seed)
  endif
!
  iangvar=2
! End variables and settings for classification by J_parallel and ideal orbit condition
!

!  open(unit=10000+ipart, iostat=stat, status='old')
!  if (stat == 0) close(10000+ipart, status='delete')
!  open(unit=20000+ipart, iostat=stat, status='old')
!  if (stat == 0) close(20000+ipart, status='delete')

! Write out trapped-passing boundary at the classification cut:
  if(class_plot) then
    if(ipart.eq.1) then
      z(1)=zstart(1,ipart)
      z(3)=cut_in_per*fper
      do kt=0,1000
        z(2)=1d-3*twopi*dble(kt)
        if(isw_field_type.eq.0) then
          call magfie_can(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
        elseif(isw_field_type.eq.1) then
          call magfie_vmec(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
        elseif(isw_field_type.eq.2) then
          call magfie_boozer(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
        else
          print *,'unknown field type'
        endif
        write(iaaa_bou,*) z(2),sqrt(1.d0-bmod/bmax)
      enddo
    endif
  endif
! End write out trapped-passing boundary at the classification cut
!
  z = zstart(:, ipart)
  r=z(1)
  theta_vmec=z(2)
  varphi_vmec=z(3)
!
  if(isw_field_type.eq.0) then
      call vmec_to_can(r,theta_vmec,varphi_vmec,z(2),z(3))
  elseif(isw_field_type.eq.2) then
      call vmec_to_boozer(r,theta_vmec,varphi_vmec,z(2),z(3))
  endif

! In case of classification plot all starting points are moved to the classification cut:
  if(class_plot) then
    z(3)=cut_in_per*fper
    zstart(2,ipart)=modulo(zstart(2,ipart),twopi)
  endif
! End moving starting points to the classification cut

  if (integmode>0) call init_sympl(anorb%si, anorb%f, z, dtaumin, dtaumin, relerr, integmode)

  if(isw_field_type.eq.0) then
      call magfie_can(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
  elseif(isw_field_type.eq.1) then
      call magfie_vmec(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
  elseif(isw_field_type.eq.2) then
      call magfie_boozer(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
  else
      print *,'unknown field type'
  endif
!
  if(.not.local) then
!
    call get_bminmax(z(1),bmin,bmax)
!
  endif
!
  passing = z(5)**2.gt.1.d0-bmod/bmax
  trap_par(ipart) = ((1.d0-z(5)**2)*bmax/bmod-1.d0)*bmin/(bmax-bmin)
  perp_inv(ipart) = z(4)**2*(1.d0-z(5)**2)/bmod
  iclass(:,ipart) = 0

! Forced classification of passing as regular:
  if(passing.and.(notrace_passing.eq.1 .or. trap_par(ipart).le.contr_pp)) then
    ! passing particle
    ! no tracing of passing particles, assume that all are confined
    ! or: strongly passing particles that are certainly confined
    !$omp critical
    confpart_pass=confpart_pass+1.d0
    !$omp end critical
    if(class_plot) then
!$omp critical
      write (iaaa_pnt,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
!$omp end critical
    endif
    iclass(:,ipart) = 1
    return
  endif
! End forced classification of passing as regular

!$omp critical
  allocate(ipoi(nplagr),coef(0:nder,nplagr),orb_sten(6,nplagr),xp(nplagr))
!$omp end critical
  do it=1,nplagr
    ipoi(it)=it
  enddo

  nfp_tip=nfp             !<= initial array dimension for tips
  nfp_per=nfp             !<= initial array dimension for periods
!$omp critical
  allocate(zpoipl_tip(2,nfp_tip),zpoipl_per(2,nfp_per))
!$omp end critical

!  open(unit=10000+ipart, recl=1024, position='append')
!  open(unit=20000+ipart, recl=1024, position='append')

  ifp_tip=0               !<= initialize footprint counter on tips
  ifp_per=0               !<= initialize footprint counter on periods

  icounter=0
  phiper=0.0d0


  kt = 0
  !$omp atomic
  confpart_pass(1)=confpart_pass(1)+1.d0
  !$omp atomic
  confpart_trap(1)=confpart_trap(1)+1.d0

  !--------------------------------
      ! Initialize tip detector

  itip=npl_half+1
  alam_prev=z(5)

  ! End initialize tip detector
  !--------------------------------
  ! Initialize period crossing detector

  iper=npl_half+1
  kper=int(z(3)/fper)

  ! End initialize period crossing detector
  !--------------------------------
!
! Initialize classification by J_parallel and ideal orbit condition:
  nfp_cot=0
! End Initialize classification by J_parallel and ideal orbit condition
!
  par_inv = 0d0
  regular = .False.
  do it=2,ntimstep
    if (regular) then  ! regular orbit, will not be lost
      if(passing) then
        !$omp atomic
        confpart_pass(it)=confpart_pass(it)+1.d0
      else
        !$omp atomic
        confpart_trap(it)=confpart_trap(it)+1.d0
      endif
      kt = kt+ntau
      cycle
    endif
    do ktau=1,ntau
      if (integmode <= 0) then
        call orbit_timestep_axis(z, dtaumin, dtaumin, relerr, ierr)
      else
        if (swcoll) then  ! TODO: move this inside modules
          anorb%si%pabs = z(4)
          anorb%f%vpar = z(4)*z(5)*sqrt2
          anorb%f%mu = z(4)**2*(1.d0-z(5)**2)/anorb%f%Bmod
          anorb%si%z(4) = anorb%f%vpar*anorb%f%hph + anorb%f%Aph/anorb%f%ro0
          call get_val(anorb%f, anorb%si%z(4)) ! for pth
          anorb%si%pthold = anorb%f%pth
        endif
        call orbit_timestep_sympl(anorb%si, anorb%f, ierr)
        z(1:3) = anorb%si%z(1:3)
        z(4) = dsqrt(anorb%f%mu*anorb%f%Bmod+0.5d0*anorb%f%vpar**2)
        z(5) = anorb%f%vpar/(z(4)*sqrt2)
      endif

      ! Collisions
      if (swcoll) then
        call stost(z, dtaumin, 1, ierr_coll)
        if (ierr_coll /= 0) then
          print *, 'Error in stost: ', ierr_coll, 'z = ', z, 'dtaumin = ', dtaumin
        endif
      endif

      ! Write starting data for orbits which were lost in case of classification plot
      if(class_plot) then
        if(ierr.ne.0) then
          !$omp critical
          if(passing) then
            write (iaaa_prp,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
          else
            write (iaaa_prt,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
          endif
          !$omp end critical
        endif
      endif
      ! End write starting data for orbits which were lost in case of classification plot

      if(ierr.ne.0) exit
      kt = kt+1

      par_inv = par_inv+z(5)**2*dtaumin ! parallel adiabatic invariant
      if(kt.le.nplagr) then          !<=first nplagr points to initialize stencil
        orb_sten(1:5,kt)=z
        orb_sten(6,kt)=par_inv
      else                          !<=normal case, shift stencil
        orb_sten(1:5,ipoi(1))=z
        orb_sten(6,ipoi(1))=par_inv
        ipoi=cshift(ipoi,1)
      endif

      ! Tip detection and interpolation
      if(alam_prev.lt.0.d0.and.z(5).gt.0.d0) itip=0   !<=tip has been passed
      itip=itip+1
      alam_prev=z(5)
      if(kt.gt.nplagr) then          !<=use only initialized stencil
        if(itip.eq.npl_half) then   !<=stencil around tip is complete, interpolate
          xp=orb_sten(5,ipoi)

          call plag_coeff(nplagr,nder,zerolam,xp,coef)

          var_tip=matmul(orb_sten(:,ipoi),coef(0,:))
          var_tip(2)=modulo(var_tip(2),twopi)
          var_tip(3)=modulo(var_tip(3),twopi)

!          write(10000+ipart,*) var_tip

          ifp_tip=ifp_tip+1
          if(ifp_tip.gt.nfp_tip) then   !<=increase the buffer for banana tips
            !$omp critical
            allocate(dummy2d(2,ifp_tip-1))
            !$omp end critical
            dummy2d=zpoipl_tip(:,1:ifp_tip-1)
            !$omp critical
            deallocate(zpoipl_tip)
            !$omp end critical
            nfp_tip=nfp_tip+nfp
            !$omp critical
            allocate(zpoipl_tip(2,nfp_tip))
            !$omp end critical
            zpoipl_tip(:,1:ifp_tip-1)=dummy2d
            !$omp critical
            deallocate(dummy2d)
            !$omp end critical
          endif
          zpoipl_tip(:,ifp_tip)=var_tip(1:2)
          par_inv = par_inv - var_tip(6)
!
! Classification by J_parallel and ideal orbit conditions:
          fpr_in(1)=var_tip(1)
          fpr_in(2)=var_tip(iangvar)
          fpr_in(3)=var_tip(6)
!
          call check_orbit_type(nturns,nfp_cot,fpr_in,ideal,ijpar,ierr_cot)
!
          iclass(1,ipart) = ijpar
          iclass(2,ipart) = ideal
          if(fast_class) ierr=ierr_cot
!
! End classification by J_parallel and ideal orbit conditions
        endif
      endif
      ! End tip detection and interpolation

      ! Periodic boundary footprint detection and interpolation
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
      if(kt.gt.nplagr) then          !<=use only initialized stencil
        if(iper.eq.npl_half) then   !<=stencil around periodic boundary is complete, interpolate
          xp=orb_sten(3,ipoi)-phiper

          call plag_coeff(nplagr,nder,zerolam,xp,coef)

          var_tip=matmul(orb_sten(:,ipoi),coef(0,:))
          var_tip(2)=modulo(var_tip(2),twopi)
          var_tip(3)=modulo(var_tip(3),twopi)
! write(20000+ipart,*) var_tip
          ifp_per=ifp_per+1
          if(ifp_per.gt.nfp_per) then   !<=increase the buffer for periodic boundary footprints
            !$omp critical
            allocate(dummy2d(2,ifp_per-1))
            !$omp end critical
            dummy2d=zpoipl_per(:,1:ifp_per-1)
            !$omp critical
            deallocate(zpoipl_per)
            !$omp end critical
            nfp_per=nfp_per+nfp
            !$omp critical
            allocate(zpoipl_per(2,nfp_per))
            !$omp end critical
            zpoipl_per(:,1:ifp_per-1)=dummy2d
            !$omp critical
            deallocate(dummy2d)
            !$omp end critical
          endif
          zpoipl_per(:,ifp_per)=var_tip(1:2)
        endif
      endif
      ! End periodic boundary footprint detection and interpolation

      ! Cut classification into regular or chaotic
      if (kt == ntcut) then
        regular = .True.

        if(ifp_per > 0) then

          call fract_dimension(ifp_per,zpoipl_per(:,1:ifp_per),fraction)

          if(fraction.gt.0.2d0) then
            print *, ipart, ' chaotic per ', ifp_per
            regular = .False.
          else
            print *, ipart, ' regular per', ifp_per
          endif
        endif

        if(ifp_tip > 0) then

          call fract_dimension(ifp_tip,zpoipl_tip(:,1:ifp_tip),fraction)

          if(fraction.gt.0.2d0) then
            print *, ipart, ' chaotic tip ', ifp_tip
            regular = .False.
            iclass(3,ipart) = 2
          else
            print *, ipart, ' regular tip ', ifp_tip
            iclass(3,ipart) = 1
          endif
        endif

        if(class_plot) then
!$omp critical
! Output of classification by Minkowsky dimension:
          if(regular) then
            if(passing) then
              write (iaaa_rep,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
            else
              write (iaaa_ret,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
            endif
          else
            if(passing) then
              write (iaaa_stp,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
            else
              write (iaaa_stt,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
            endif
          endif
!End output of classification by Minkowsky dimension
!$omp end critical
          ierr=1
        endif
      endif
!
      if(ierr.ne.0) then
        if(class_plot) then
! Output of classification by J_parallel and ideal orbit condition:
          if(.not.passing) then
!$omp critical
            select case(ijpar)
            case(0)
              write (iaaa_jer,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
            case(1)
              write (iaaa_jre,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
            case(2)
              write (iaaa_jst,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
            end select
!
            select case(ideal)
            case(0)
              write (iaaa_ier,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
            case(1)
              write (iaaa_ire,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
            case(2)
              write (iaaa_ist,*) zstart(2,ipart),zstart(5,ipart),trap_par(ipart)
            end select
!$omp end critical
          endif
        endif
! End output of classification by J_parallel and ideal orbit condition
        exit
      endif
!    write(999, *) kt*dtaumin/v0, z
    enddo
    if(ierr.ne.0) exit
    if(passing) then
      !$omp atomic
      confpart_pass(it)=confpart_pass(it)+1.d0
    else
      !$omp atomic
      confpart_trap(it)=confpart_trap(it)+1.d0
    endif
  enddo
  zend(:,ipart) = z
  times_lost(ipart) = kt*dtaumin/v0
  !$omp critical
  deallocate(zpoipl_tip, zpoipl_per)
  !$omp end critical
!  close(unit=10000+ipart)
!  close(unit=10000+ipart)
end subroutine trace_orbit

end program neo_orb_main
