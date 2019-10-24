program neo_orb_main
  use omp_lib
  use common, only: pi, twopi, c, e_charge, e_mass, p_mass, ev
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp

  use parmot_mod, only : ro0
  use velo_mod,   only : isw_field_type
  use field_can_mod, only : FieldCan
  use orbit_symplectic, only : SymplecticIntegrator, orbit_timestep_sympl
  use neo_orb, only : init_field, init_sympl, NeoOrb
  use cut_detector, only : fract_dimension
  use diag_mod, only : icounter
  
  implicit none

  integer          :: npoi,L1i,nper,npoiper,i,ntimstep,ntestpart
  integer          :: notrace_passing,loopskip,iskip
  double precision :: dphi,phibeg,bmod00,rlarm,bmax,bmin
  double precision :: tau,dtau,dtaumin,xi,v0,bmod_ref,E_alpha,trace_time
  double precision :: RT0,R0i,cbfi,bz0i,bf0,rbig
  double precision :: sbeg,thetabeg
  double precision, dimension(:),   allocatable :: bstart,volstart
  double precision, dimension(:,:), allocatable :: xstart
  double precision, dimension(:,:), allocatable :: zstart
  double precision, dimension(:), allocatable :: confpart_trap,confpart_pass
  double precision, dimension(:), allocatable :: times_lost
  integer          :: npoiper2
  double precision :: contr_pp
  double precision :: facE_al
  integer          :: ibins
  integer          :: n_e,n_d
  integer          :: startmode

  integer :: ntau ! number of dtaumin in dtau
  integer :: integmode = 0 ! 0 = RK, 1 = Euler1, 2 = Euler2, 3 = Verlet

  integer :: kpart = 0 ! progress counter for particles

  double precision :: relerr

  type(NeoOrb) :: norb
  double precision, allocatable :: trap_par(:)

  integer, parameter :: n_tip_vars = 6  ! variables to evaluate at tip: z(1..5), par_inv
  integer :: nplagr,nder,npl_half
  integer :: norbper,nfp
  double precision :: fper, zerolam = 0d0

! read config file
  call read_config

! initialize field geometry
  call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
  call init_params
  print *, 'tau: ', dtau, dtaumin, min(dabs(mod(dtau, dtaumin)), &
                    dabs(mod(dtau, dtaumin)-dtaumin))/dtaumin, ntau

! pre-compute starting flux surface
  npoi=nper*npoiper ! total number of starting points
  allocate(xstart(3,npoi),bstart(npoi),volstart(npoi))
  call init_starting_surf

! initialize array of confined particle percentage
  allocate(confpart_trap(ntimstep),confpart_pass(ntimstep))
  confpart_trap=0.d0
  confpart_pass=0.d0

! initialize lost times when particles get lost
  allocate(times_lost(ntestpart), trap_par(ntestpart))
  times_lost = -1.d0

  allocate(zstart(5,ntestpart))
  call init_starting_points
  if (startmode == 0) stop

  icounter=0 ! evaluation counter

! do particle tracing in parallel

  !$omp parallel private(norb)
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
    write(1,*) i, times_lost(i), trap_par(i)
  enddo
  close(1)

  if (integmode >= 0) call deallocate_can_coord
  deallocate(times_lost, confpart_trap, confpart_pass, trap_par)

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
  read (1,*) startmode         !mode for initial conditions: 0=generate and store, 1=generate, store, and run, 2=read and run
  read (1,*) integmode         !mode for integrator: -1 = RK VMEC, 0 = RK CAN, 1 = Euler1, 2 = Euler2, 3 = Verlet
  read (1,*) relerr            !relative error for RK integrator
  close(1)
end subroutine read_config

subroutine init_params
! set alpha energy, velocity, and Larmor radius
  E_alpha=3.5d6/facE_al
  v0=sqrt(2.d0*E_alpha*ev/(n_d*p_mass))
  rlarm=v0*n_d*p_mass*c/(n_e*e_charge*bmod_ref)

! normalized slowing down time:
  tau=trace_time*v0
! normalized time step:
  dtau=tau/dble(ntimstep-1)
! parameters for the vacuum chamber:
  call stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0) ! TODO: why again?
  rbig=rt0
! field line integration step step over phi (to check chamber wall crossing)
  dphi=2.d0*pi/(L1i*npoiper)
! orbit integration time step (to check chamber wall crossing)
  dtaumin=2.d0*pi*rbig/npoiper2
  ntau=ceiling(dtau/dtaumin)
  dtaumin=dtau/ntau

  norbper=ceiling(1d0*ntau*ntimstep/(L1i*npoiper2))
  nfp=L1i*norbper         !<= guess for footprint number

  zerolam=0.d0
  nplagr=4
  nder=0
  npl_half=nplagr/2

  fper = 2d0*pi/dble(L1i)   !<= field period
end subroutine init_params

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

subroutine init_starting_points
  integer :: ipart
  real :: zzg
  
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
      ! coordinates: z(1) = R, z(2) = phi, z(3) = Z
      zstart(1:3,ipart)=xstart(:,i)
      ! normalized velocity module z(4) = v / v_0:
      zstart(4,ipart)=1.d0
      ! starting pitch z(5)=v_\parallel / v:
      xi=zzg()
      zstart(5,ipart)=2.d0*(xi-0.5d0)
      write(1,*) zstart(:,ipart)
    enddo
  else
    do ipart=1,ntestpart
      read(1,*) zstart(:,ipart)
    enddo
  endif
  
  close(1)
end subroutine init_starting_points

subroutine trace_orbit(anorb, ipart)
  type(NeoOrb), intent(inout) :: anorb
  integer, intent(in) :: ipart
  integer :: ierr
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
  integer :: stat
  double precision :: phiper, alam_prev, par_inv
  integer :: iper, itip, kper, nfp_tip, nfp_per

  double precision :: fraction
  logical :: regular
  integer :: ntcut


  !open(unit=10000+ipart, iostat=stat, status='old')
  !if (stat == 0) close(10000+ipart, status='delete')
  !open(unit=20000+ipart, iostat=stat, status='old')
  !if (stat == 0) close(20000+ipart, status='delete')

  z = zstart(:, ipart)
  if (integmode>0) call init_sympl(anorb%si, anorb%f, z, dtaumin, dtaumin, relerr, integmode)

  if(isw_field_type.eq.0) then
      call magfie_can(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
  else
      call magfie_vmec(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
  endif

  passing = z(5)**2.gt.1.d0-bmod/bmax
  trap_par(ipart) = ((1.d0-z(5)**2)*bmax/bmod-1.d0)*bmin/(bmax-bmin)

  if(passing.and.(notrace_passing.eq.1 .or. trap_par(ipart).le.contr_pp)) then
    ! passing particle
    ! no tracing of passing particles, assume that all are confined
    ! or: strongly passing particles that are certainly confined
    !$omp critical
    confpart_pass=confpart_pass+1.d0
    !$omp end critical
    return
  endif

  allocate(ipoi(nplagr),coef(0:nder,nplagr),orb_sten(6,nplagr),xp(nplagr))
  do it=1,nplagr
    ipoi(it)=it
  enddo

  nfp_tip=nfp             !<= initial array dimension for tips
  nfp_per=nfp             !<= initial array dimension for periods
  allocate(zpoipl_tip(2,nfp_tip),zpoipl_per(2,nfp_per))

  !open(unit=10000+ipart, recl=1024, position='append')
  !open(unit=20000+ipart, recl=1024, position='append')

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

  par_inv = 0d0
  regular = .False.
  ntcut = ntimstep*ntau/10
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
        call orbit_timestep_sympl(anorb%si, anorb%f, ierr)
        z(1:3) = anorb%si%z(1:3)
        z(4) = 1d0
        z(5) = anorb%f%vpar/dsqrt(2d0)
      endif
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

          !write(10000+ipart,*) var_tip

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
          zpoipl_tip(:,ifp_tip)=var_tip(1:2)
          par_inv = par_inv - var_tip(6)
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
          !write(20000+ipart,*) var_tip
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
          else
            print *, ipart, ' regular tip ', ifp_tip
          endif
        endif
      endif
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

  times_lost(ipart) = kt*dtaumin/v0
  deallocate(zpoipl_tip, zpoipl_per)
end subroutine trace_orbit

end program neo_orb_main
