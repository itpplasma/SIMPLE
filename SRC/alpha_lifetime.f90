program alpha_lifetime
  use omp_lib
  use common, only: pi, c, e_charge, e_mass, p_mass, ev
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_A, ns_s, ns_tp
  use chamb_mod,  only : rnegflag
  use parmot_mod, only : rmu, ro0, eeff
  use velo_mod,   only : isw_field_type
use diag_mod, only : dodiag

  implicit none

  integer          :: npoi,ierr,L1i,nper,npoiper,i,ntimstep,ntestpart
  integer          :: ipart,notrace_passing,loopskip,iskip,ilost
  real             :: zzg
  double precision :: dphi,rbeg,phibeg,zbeg,bmod00,rcham,rlarm,bmax,bmin
  double precision :: tau,dtau,dtaumin,xi,v0,bmod_ref,E_alpha,trace_time
  double precision :: RT0,R0i,cbfi,bz0i,bf0,trap_par,rbig
  double precision :: sbeg,thetabeg
  double precision, dimension(5) :: z
  double precision, dimension(:),   allocatable :: bstart,volstart,confpart
  double precision, dimension(:,:), allocatable :: xstart
  double precision, dimension(:,:), allocatable :: zstart
  double precision, dimension(:), allocatable :: confpart_trap,confpart_pass
  integer          :: npoiper2
  double precision :: contr_pp
  double precision :: facE_al
  integer          :: ibins
  integer          :: n_e,n_d,n_b
  integer          :: startmode

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
  read (1,*) contr_pp          !control of passing particle fraction
  read (1,*) facE_al           !facE_al test particle energy reduction factor
  read (1,*) npoiper2          !additional integration step split factor
  read (1,*) n_e               !test particle charge number (the same as Z)
  read (1,*) n_d               !test particle mass number (the same as A)
  read (1,*) netcdffile        !name of VMEC file in NETCDF format <=2017 NEW
  read (1,*) ns_s              !spline order for 3D quantities over s variable
  read (1,*) ns_tp             !spline order for 3D quantities over theta and phi
  read (1,*) multharm          !angular grid factor (n_grid=multharm*n_harm_max where n_harm_max - maximum Fourier index)
  read (1,*) startmode         !mode for initial conditions: 0=generate and store, 1=generate, store, and run, 2=read and run
  close(1)

! initialize field geometry
  call spline_vmec_data ! initialize splines for VMEC field
  call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0) ! initialize periods and major radius 
  call get_canonical_coordinates ! pre-compute transformation to canonical coords
  isw_field_type = 0 ! evaluate fields in canonical coords (0 = CAN, 1 = VMEC)

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
  dtau=tau/dfloat(ntimstep-1)
! parameters for the vacuum chamber:
  call stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0) ! TODO: why again?
  rbig=rt0
! field line integration step step over phi (to check chamber wall crossing)
  dphi=2.d0*pi/(L1i*npoiper)
! orbit integration time step (to check chamber wall crossing)
  dtaumin=dphi*rbig/npoiper2

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
  close(1)

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

! do particle tracing in parallel
!$omp parallel private(ibins, xi, i, z, trap_par)
!$omp do
  do ipart=1,ntestpart
    print *, ipart, ' / ', ntestpart, 'thread: ', omp_get_thread_num()
    z = zstart(:,ipart)    

    if(z(5)**2.gt.1.d0-bstart(i)/bmax) then
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
      ilost=ntimstep-1
      do i=2,ntimstep
      if(trap_par.le.contr_pp) go to 111
        call orbit_timestep_axis(z,dtau,dtaumin,ierr)
        if(ierr.ne.0) exit
  111  continue
        ilost=ntimstep-i
print *,'passing particle ',ipart,' step ',i,' of ',ntimstep
!$omp atomic
        confpart_pass(i)=confpart_pass(i)+1.d0
      enddo
    else
! trapped particle (traced always)
      confpart_trap(1)=confpart_trap(1)+1.d0
      ilost=ntimstep-1
      do i=2,ntimstep

        call orbit_timestep_axis(z,dtau,dtaumin,ierr)

        if(ierr.ne.0) exit
        ilost=ntimstep-i

        confpart_trap(i)=confpart_trap(i)+1.d0
      enddo
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

  call deallocate_can_coord
end program alpha_lifetime
