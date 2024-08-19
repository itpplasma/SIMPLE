program test_coord_trans
  use omp_lib
  use util, only: pi, twopi, c, e_charge, e_mass, p_mass, ev
  use new_vmec_stuff_mod, only : netcdffile, multharm, ns_s, ns_tp

  use parmot_mod, only : ro0, rmu
  use field_can_mod, only : FieldCan
  use orbit_symplectic, only : SymplecticIntegrator, orbit_timestep_sympl
  use simple, only : init_field, init_sympl
  use cut_detector, only : fract_dimension
  use params, only : Tracer, debug
  use binsrc_sub, only: binsrc
  use get_can_sub, only: can_to_vmec, vmec_to_can
  use alpha_lifetime_sub, only: integrate_mfl_can

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

  type(Tracer) :: norb
  double precision, allocatable :: trap_par(:)

  integer, parameter :: n_tip_vars = 6  ! variables to evaluate at tip: z(1..5), par_inv
  integer :: nplagr,nder,npl_half
  integer :: norbper,nfp
  double precision :: fper, zerolam = 0d0

  double precision :: tcut
  integer :: ntcut
  logical          :: class_plot     !<=AAA
  double precision :: cut_in_per     !<=AAA

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

  allocate(zstart(5,ntestpart))
  call init_starting_points
  if (startmode == 0) stop

  call testcoordtrans

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
  read (1,*) startmode         !mode for initial conditions: 0=generate and store, 1=generate, store, and run, 2=read and run, 3=read ANTS and run
  read (1,*) integmode         !mode for integrator: -1 = RK VMEC, 0 = RK CAN, 1 = Euler1, 2 = Euler2, 3 = Verlet
  read (1,*) relerr            !relative error for RK integrator
  read (1,*) tcut              !time when to do cut for classification, usually 1d-1, or -1 if no cuts desired
  read (1,*) debug             !produce debugging output (.True./.False.). Use only in non-parallel mode!
  read (1,*) class_plot        !write starting points at phi=const cut for classification plot (.True./.False.).  !<=AAA
  read (1,*) cut_in_per        !normalized phi-cut position within field period, [0:1], used if class_plot=.True. !<=AAA
  close(1)
end subroutine read_config

subroutine init_params
! set alpha energy, velocity, and Larmor radius
  E_alpha=3.5d6/facE_al
  v0=sqrt(2.d0*E_alpha*ev/(n_d*p_mass))
  rlarm=v0*n_d*p_mass*c/(n_e*e_charge*bmod_ref)

! Neglect relativistic effects by large inverse relativistic temperature
  rmu=1d8

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

  ntcut = ceiling(ntimstep*ntau*tcut/trace_time)

  norbper=ceiling(1d0*ntau*ntimstep/(L1i*npoiper2))
  nfp=L1i*norbper         !<= guess for footprint number

  zerolam=0.d0
  nplagr=4
  nder=0
  npl_half=nplagr/2

  fper = 2d0*pi/dble(L1i)   !<= field period
end subroutine init_params

subroutine init_starting_surf
  use alpha_lifetime_sub, only: integrate_mfl_can
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
  integer :: ipart

  do ipart=1,ntestpart
    read(unit, '(A)') line
    call process_line(line, v_par, v_perp, u, v, s)
    print *, v_par, v_perp, u, v, s
  enddo
end subroutine

subroutine init_starting_points
  use binsrc_sub

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
  else if (startmode == 2) then
    do ipart=1,ntestpart
      read(1,*) zstart(:,ipart)
    enddo
  else if (startmode == 3) then  ! ANTS input mode
    call init_starting_points_ants(1)
  endif

  close(1)
end subroutine init_starting_points

!--------------------------------------------------

subroutine testcoordtrans
use get_can_sub
!
integer :: it,ip,nt,np,nper
double precision :: r,vartheta_c,varphi_c,theta_vmec,varphi_vmec,ht,hp
double precision :: vartheta_c_back,varphi_c_back
!
nper=5
nt=100
np=101
ht=2.d0*pi/dble(nt)
hp=2.d0*pi/dble(np*nper)
r=0.7
do it=0,nt
  do ip=0,np
    vartheta_c=ht*dble(it)
    varphi_c=hp*dble(ip)
!
    call can_to_vmec(r,vartheta_c,varphi_c,theta_vmec,varphi_vmec)
    call vmec_to_can(r,theta_vmec,varphi_vmec,vartheta_c_back,varphi_c_back)
!
    write(20001,*) vartheta_c,varphi_c,vartheta_c_back-vartheta_c
    write(20002,*) vartheta_c,varphi_c,varphi_c_back-varphi_c
  enddo
  write(20001,*) ' '
  write(20002,*) ' '
enddo

stop

end subroutine testcoordtrans

!-----------------------------------------

end program test_coord_trans
