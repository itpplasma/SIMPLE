module params
#ifdef MPI
    use mpi
#endif
  use util
  use parmot_mod, only : ro0, rmu
  use new_vmec_stuff_mod, only : old_axis_healing, old_axis_healing_boundary, &
    netcdffile, ns_s, ns_tp, multharm, vmec_B_scale, vmec_RZ_scale
  use velo_mod,   only : isw_field_type
  use field_can_mod, only : FieldCan
  use orbit_symplectic, only : SymplecticIntegrator, MultistageIntegrator

  implicit none
  save

  type :: Tracer
    double precision :: fper
    double precision :: dtau, dtaumin, v0
    integer          :: n_e, n_d

    integer :: integmode = 0 ! 0 = RK, 1 = Euler1, 2 = Euler2, 3 = Verlet
    double precision :: relerr

    type(FieldCan) :: f
    type(SymplecticIntegrator) :: si
    type(MultistageIntegrator) :: mi
  end type Tracer

  integer          :: npoi, L1i, nper=1000, i, ntestpart=1024
  integer          :: loopskip=0,iskip
  double precision :: dphi,phibeg,bmod00,rlarm,bmax,bmin
  double precision :: tau,dtau,dtaumin,xi
  double precision :: RT0,R0i,cbfi,bz0i,bf0,rbig
  double precision :: sbeg=0.5d0, thetabeg=0.0d0
  double precision, dimension(:), allocatable :: sbeg_arr
  double precision, dimension(:),   allocatable :: bstart,volstart
  double precision, dimension(:,:), allocatable :: xstart
  double precision, dimension(:,:), allocatable :: zstart, zend
  double precision, dimension(:), allocatable :: confpart_trap,confpart_pass
  double precision, dimension(:), allocatable :: times_lost
  double precision :: contr_pp=-1d0
  integer          :: ibins
  integer          :: startmode=1

  integer :: ntau ! number of dtaumin in dtau
  integer :: integmode = 1 ! 0 = RK, 1 = Euler1, 2 = Euler2, 3 = Verlet

  integer :: kpart = 0 ! progress counter for particles

  double precision :: relerr = 1d-13

  double precision, allocatable :: trap_par(:), perp_inv(:)
  integer,          allocatable :: iclass(:,:)

  integer, parameter :: n_tip_vars = 6  ! variables to evaluate at tip: z(1..5), par_inv
  integer :: nplagr,nder,npl_half
  integer :: norbper,nfp
  double precision :: fper, zerolam = 0d0

  double precision :: tcut = -1d0
  integer :: ntcut
  logical          :: class_plot = .False.    !<=AAA
  double precision :: cut_in_per = 0d0        !<=AAA
  logical          :: local=.True.

  logical :: fast_class=.False.  !if .true. quit immeadiately after fast classification

  ! Colliding with D-T reactor plasma. TODO: Make configurable
  logical :: swcoll = .False.
  double precision, parameter :: am1=2.0d0, am2=3.0d0, Z1=1.0d0, Z2=1.0d0, &
    densi1=0.5d14, densi2=0.5d14, tempi1=1.0d4, tempi2=1.0d4, tempe=1.0d4
  double precision :: dchichi,slowrate,dchichi_norm,slowrate_norm
  logical :: deterministic = .False.

  ! Further configuration parameters
  integer          :: notrace_passing = 0
  double precision :: facE_al=1d0, trace_time=1d-1
  integer :: ntimstep=10000, npoiper=100, npoiper2=256, n_e=2, n_d=4

  double precision :: v0

  integer :: nfirstpart, nlastpart

  logical :: debug = .False.
  integer :: ierr
#ifdef MPI
  integer :: mpirank, mpisize

  double precision, dimension(:), allocatable :: &
    confpart_trap_glob, confpart_pass_glob
#endif
  integer :: batch_size=10000
  integer :: ran_seed=12345
  logical :: reuse_batch =.False.
  logical :: sbeg_multi = .False.
  integer, dimension (:), allocatable :: idx

  namelist /config/ notrace_passing, nper, npoiper, ntimstep, ntestpart, &
    trace_time, sbeg_multi, sbeg_arr, sbeg, phibeg, thetabeg, loopskip, contr_pp,              &
    facE_al, npoiper2, n_e, n_d, netcdffile, ns_s, ns_tp, multharm,      &
    isw_field_type, startmode, integmode, relerr, tcut, debug,           &
    class_plot, cut_in_per, fast_class, local, vmec_B_scale,             &
    vmec_RZ_scale, swcoll, deterministic, old_axis_healing,              &
    old_axis_healing_boundary, &
    batch_size, ran_seed, reuse_batch

contains

  subroutine read_config
    integer :: iostat
    character(1024) :: iomsg

    open(1, file='simple.in', recl=1024, iostat=iostat, iomsg=iomsg)
    if (iostat /= 0) goto 666

    read(1, nml=config, iostat=iostat, iomsg=iomsg)
    if (iostat /= 0) goto 666

    close(1)

    if (swcoll .and. (tcut > 0.0d0 .or. class_plot .or. fast_class)) then
      stop 'Collisions are incompatible with classification'
    endif

    return

    666 stop iomsg
  end subroutine read_config


  subroutine params_init
    double precision :: E_alpha
    integer :: iostat, i, n
    character, dimension (:), allocatable :: batch_file
    logical :: old_batch
    real :: ran_tmp

    E_alpha = 3.5d6/facE_al

  ! set alpha energy, velocity, and Larmor radius
    v0=sqrt(2.d0*E_alpha*ev/(n_d*p_mass))
    rlarm=v0*n_d*p_mass*c/(n_e*e_charge)
    ro0=rlarm

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
    dtaumin=2.d0*pi*rbig/npoiper2  ! ntimstep =
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

    npoi=nper*npoiper ! total number of starting points

    !See if batch is wanted, if yes find random or re-use from previous file.
    if (ntestpart > batch_size) then
      if (reuse_batch) then
        INQUIRE(FILE="batch.dat", EXIST=old_batch)

        if (old_batch) then
          allocate(batch_file(batch_size))
          allocate(idx(batch_size))
          open(1, file='batch.dat', recl=batch_size, iostat=iostat)
          if (iostat /= 0) goto 666

          do i=1,batch_size
            read(1,iostat=iostat) batch_file(i)
            if (iostat /= 0) goto 667
            idx (i) = ICHAR(batch_file(i)) ! TODO batch_file list is pointless
          end do
          deallocate(batch_file)
          close(1)

        else !old_batch
          !Create a random list(batch_size) of indices using ran_seed.
          allocate(idx(batch_size))
          call SRAND(ran_seed)
          do i=0,batch_size
            
            call random_number(ran_tmp)
            idx(i) = FLOOR(batch_size * ran_tmp)
          end do
          call sort_idx(idx, batch_size)
          
          open(1, file='batch.dat', recl=batch_size*2, iostat=iostat)
          do n=1, batch_size
            write(1, *) idx(i)
          end do
          close(1)
        endif !old_batch

      endif !reuse_batch

      !Set ntestpart to batch_size for rest of the run.
      ntestpart = batch_size
    endif !batches wanted

    allocate(zstart(5,ntestpart), zend(5,ntestpart))
    allocate(times_lost(ntestpart), trap_par(ntestpart), perp_inv(ntestpart))
    allocate(xstart(3,npoi),bstart(npoi),volstart(npoi))
    allocate(confpart_trap(ntimstep),confpart_pass(ntimstep))
#ifdef MPI
      allocate(confpart_trap_glob(ntimstep),confpart_pass_glob(ntimstep))
#endif
    allocate(iclass(3,ntestpart))

    666 stop iostat
    667 stop iostat
  end subroutine params_init
  
  
  ! TODO sort, then check for duplicates
  subroutine sort_idx(idx_arr, N)
    integer, dimension (N) :: idx_arr
    integer :: N
    ! TODO implement
  end subroutine sort_idx

end module params
