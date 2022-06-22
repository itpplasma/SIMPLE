program alpha_lifetime
  use omp_lib
  use util, only: pi, c, e_charge, e_mass, p_mass, ev
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp

  use parmot_mod, only : rmu, ro0
  use velo_mod,   only : isw_field_type
  use field_can_mod, only : FieldCan
  use orbit_symplectic, only : SymplecticIntegrator, orbit_timestep_sympl
  use simple, only : init_sympl
use diag_mod, only : icounter

  implicit none

  integer          :: npoi,ierr,L1i,nper,npoiper,i,ntimstep,ntestpart
  integer          :: ipart,notrace_passing,loopskip,iskip
  real             :: zzg
  double precision :: dphi,phibeg,bmod00,rlarm,bmax,bmin
  double precision :: tau,dtau,dtaumin,xi,v0,bmod_ref,E_alpha,trace_time
  double precision :: RT0,R0i,cbfi,bz0i,bf0,rbig
  double precision :: sbeg,thetabeg
  double precision, dimension(5) :: z
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

  double precision :: bmod,sqrtg
  double precision, dimension(3) :: bder,hcovar,hctrvr,hcurl

  double precision :: relerr

  type(FieldCan) :: f
  type(SymplecticIntegrator) :: si

  rmu=1d5 ! inverse relativistic temperature

  open(1,file='alpha_lifetime.inp',recl=1024)
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

! initialize field geometry
  call spline_vmec_data ! initialize splines for VMEC field
  call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0) ! initialize periods and major radius
  isw_field_type = 1 ! evaluate fields in VMEC coords (0 = CAN, 1 = VMEC)
  if (integmode>=0) then
    call get_canonical_coordinates ! pre-compute transformation to canonical coords
    isw_field_type = 0 ! evaluate fields in canonical coords (0 = CAN, 1 = VMEC)
  end if

! initialize position and do first check if z is inside vacuum chamber
  z = 0.0d0
  call chamb_can(z(1:2), z(3), ierr)

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
  print *, 'tau: ', dtau, dtaumin, min(dabs(mod(dtau, dtaumin)), dabs(mod(dtau, dtaumin)-dtaumin))/dtaumin, ntau

! pre-compute starting points
  npoi=nper*npoiper ! total number of starting points
  allocate(xstart(3,npoi),bstart(npoi),volstart(npoi),zstart(5,ntestpart))
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

  open(1,file='bminmax.dat',recl=1024)
  write(1,*)bmin,bmax,bmod00
  close(1)
  open(1,file='volstart.dat',recl=1024)
  do i=1,npoi
  write(1,*)i,volstart(i)
  end do
  close(1)

  open(1,file='starting_surface.dat',recl=1024)
  do i=1,npoi,npoiper
    write (1,*) xstart(1,i),xstart(3,i)
  enddo
  close(1)

! initialize array of confined particle percentage
  allocate(confpart_trap(ntimstep),confpart_pass(ntimstep))
  confpart_trap=0.d0
  confpart_pass=0.d0

! initialize lost times when particles get lost
  allocate(times_lost(ntestpart))
  times_lost = -1.d0

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

  if (startmode == 0) stop

  icounter=0 ! evaluation counter

! do particle tracing in parallel
!$omp parallel private(ibins, xi, i, z, ierr, bmod, sqrtg, bder, hcovar, hctrvr, hcurl, si, f)
!$omp do
  do ipart=1,ntestpart
!$omp atomic
    kpart = kpart+1
    print *, kpart, ' / ', ntestpart, 'particle: ', ipart, 'thread: ', omp_get_thread_num()
    z = zstart(:,ipart)
    if (integmode>0) call init_sympl(si, f, z, dtau, dtaumin, relerr, integmode)

    if(isw_field_type.eq.0) then
        call magfie_can(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
    else
        call magfie_vmec(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
    endif

    if(z(5)**2.gt.1.d0-bmod/bmax) then
    ! passing particle
      if(notrace_passing.eq.1) then
      ! no tracing of passing particles, assume that all are confined
!$omp critical
        confpart_pass=confpart_pass+1.d0
!$omp end critical
        cycle
      endif
      ! trace passing particle
!$omp atomic
      confpart_pass(1)=confpart_pass(1)+1.d0
      do i=2,ntimstep
        if (integmode <= 0) then
          call orbit_timestep_axis(z, dtau, dtaumin, relerr, ierr)
        else
          call orbit_timestep_sympl(si, f, ierr)
        endif
        if(ierr.ne.0) exit
!$omp atomic
        confpart_pass(i)=confpart_pass(i)+1.d0
      enddo
      if (integmode <= 0) then
        times_lost(ipart) = dble(i-1)*dtau/v0
      else
        times_lost(ipart) = si%kt*dtaumin/v0  ! more accurate sub-steps
      endif
    else
! trapped particle (traced always)
!$omp atomic
      confpart_trap(1)=confpart_trap(1)+1.d0
      do i=2,ntimstep
        if (integmode <= 0) then
          call orbit_timestep_axis(z, dtau, dtaumin, relerr, ierr)
        else
          call orbit_timestep_sympl(si, f, ierr)
        endif
        if(ierr.ne.0) exit
!$omp atomic
        confpart_trap(i)=confpart_trap(i)+1.d0
      enddo
      if (integmode <= 0) then
        times_lost(ipart) = dble(i-1)*dtau/v0
      else
        times_lost(ipart) = si%kt*dtaumin/v0  ! more accurate sub-steps
      endif
    endif
  enddo
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
    write(1,*) i, times_lost(i)
  enddo
  close(1)

  if (integmode >= 0) call deallocate_can_coord
  deallocate(times_lost, confpart_trap, confpart_pass)

! log initial configuration
  open(1,file='alpha_lifetime.log',recl=1024)
  write (1,*) 'notrace_passing = ',notrace_passing
  write (1,*) 'nper = ',nper
  write (1,*) 'npoiper = ',npoiper
  write (1,*) 'ntimestep = ',ntimstep
  write (1,*) 'ntestpart = ',ntestpart
  write (1,*) 'bmod_ref = ',bmod_ref
  write (1,*) 'trace_time = ',trace_time
  write (1,*) 'sbeg = ',sbeg
  write (1,*) 'phibeg = ',phibeg
  write (1,*) 'thetabeg = ',thetabeg
  write (1,*) 'Rbig = ',rbig
  write (1,*) 'dphi = ',dphi
  write (1,*) 'v0 = ',v0
  write (1,*) 'rlarm = ',rlarm
  write (1,*) 'dtau = ',dtau
  write (1,*) 'dtaumin = ',dtaumin
  write (1,*) 'E_alpha = ',E_alpha
  write (1,*) 'contr_pp =', contr_pp
  write (1,*) 'n_e = ',n_e
  write (1,*) 'n_d = ',n_d
  write (1,*) 'ns_s = ',ns_s
  write (1,*) 'ns_tp = ',ns_tp
  write (1,*) 'multharm = ',multharm
  write (1,*) 'startmode = ',startmode
  write (1,*) 'npoiper2 = ',npoiper2
  write (1,*) 'integmode = ',integmode
  write (1,*) 'neval = ', icounter
  close(1)


  print *, 'integ_mode'
end program alpha_lifetime
