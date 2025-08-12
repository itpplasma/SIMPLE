module params
  use util, only: pi, c, e_charge, p_mass, ev
  use parmot_mod, only : ro0, rmu
  use new_vmec_stuff_mod, only : old_axis_healing, old_axis_healing_boundary, &
    netcdffile, ns_s, ns_tp, multharm, vmec_B_scale, vmec_RZ_scale
  use velo_mod,   only : isw_field_type
  use field_can_mod, only : eval_field => evaluate, FieldCan
  use orbit_symplectic_base, only : SymplecticIntegrator, MultistageIntegrator, &
    EXPL_IMPL_EULER
  use vmecin_sub, only : stevvo
  use callback, only : output_error, output_orbits_macrostep

  implicit none

  ! Define real(dp) kind parameter
  integer, parameter :: dp = kind(1.0d0)
  integer          :: nper=1000, ntestpart=1024
  integer          :: zstart_dim1 = 5
  real(dp) :: dphi,phibeg=0d0,bmod00,rlarm,bmax,bmin
  real(dp) :: tau,dtau,dtaumin,xi
  real(dp) :: RT0,R0i,cbfi,bz0i,bf0,rbig
  real(dp), dimension(1024) :: sbeg=0.5d0
  real(dp) :: thetabeg=0.0d0
  real(dp), dimension(:),   allocatable :: bstart,volstart !where npoi is used, add a dimension at the end for sbeg
  real(dp), dimension(:,:), allocatable :: xstart
  real(dp), dimension(:,:), allocatable :: zstart, zend
  real(dp), dimension(:), allocatable :: confpart_trap,confpart_pass
  real(dp), dimension(:), allocatable :: times_lost
  real(dp) :: contr_pp=-1d0
  integer          :: ibins
  logical          :: generate_start_only=.False.
  integer          :: startmode=1
  real(dp) :: grid_density=0d0
  logical          :: special_ants_file=.False.

  integer :: ntau ! number of dtaumin in dtau

  integer :: integmode = EXPL_IMPL_EULER

  integer :: kpart = 0 ! progress counter for particles

  real(dp) :: relerr = 1d-13

  real(dp), allocatable :: trap_par(:), perp_inv(:)
  integer,          allocatable :: iclass(:,:)

  integer, parameter :: n_tip_vars = 6  ! variables to evaluate at tip: z(1..5), par_inv
  integer :: nplagr,nder,npl_half
  integer :: norbper,nfp
  real(dp) :: fper, zerolam = 0d0

  real(dp) :: tcut = -1d0
  integer :: ntcut
  logical          :: class_plot = .False.    !<=AAA
  real(dp) :: cut_in_per = 0d0        !<=AAA

  logical :: fast_class=.False.  !if .true. quit immeadiately after fast classification

  ! Colliding with D-T reactor plasma. TODO: Make configurable
  logical :: swcoll = .False.
  real(dp) :: am1=2.0d0, am2=3.0d0, Z1=1.0d0, Z2=1.0d0, &
    densi1=0.5d14, densi2=0.5d14, tempi1=1.0d4, tempi2=1.0d4, tempe=1.0d4
  real(dp) :: dchichi,slowrate,dchichi_norm,slowrate_norm
  logical :: deterministic = .False.

  ! Further configuration parameters
  integer          :: notrace_passing = 0
  real(dp) :: facE_al=1d0, trace_time=1d-1
  integer :: ntimstep=10000, npoiper=100, npoiper2=256, n_e=2
  real(dp) :: n_d=4

  real(dp) :: v0

  logical :: debug = .False.
  integer :: ierr

  integer :: batch_size=2000000000  ! Initialize large so batch mode is not default
  integer :: ran_seed=12345
  integer :: num_surf=1
  logical :: reuse_batch =.False.
  integer, dimension (:), allocatable :: idx

  character(1000) :: field_input = ''

  namelist /config/ notrace_passing, nper, npoiper, ntimstep, ntestpart, &
    trace_time, num_surf, sbeg, phibeg, thetabeg, contr_pp,              &
    facE_al, npoiper2, n_e, n_d, netcdffile, ns_s, ns_tp, multharm,      &
    isw_field_type, generate_start_only, startmode, grid_density, special_ants_file, integmode, relerr, tcut, debug,           &
    class_plot, cut_in_per, fast_class, vmec_B_scale,             &
    vmec_RZ_scale, swcoll, deterministic, old_axis_healing,              &
    old_axis_healing_boundary, am1, am2, Z1, Z2, &
    densi1, densi2, tempi1, tempi2, tempe, &
    batch_size, ran_seed, reuse_batch, field_input, &
    output_error, output_orbits_macrostep  ! callback

contains

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Validate configuration parameters
  subroutine validate_configuration()
    ! Validates that configuration parameters are within acceptable ranges
    
    ! Check for incompatible settings
    if (swcoll .and. (tcut > 0.0d0 .or. class_plot .or. fast_class)) then
      error stop 'Collisions are incompatible with classification'
    endif
    
    ! Validate physical parameters
    if (ntestpart <= 0) then
      error stop 'Number of test particles must be positive'
    endif
    
    if (ntimstep <= 0) then
      error stop 'Number of time steps must be positive'
    endif
    
    if (trace_time <= 0.0d0) then
      error stop 'Trace time must be positive'
    endif
    
    if (relerr <= 0.0d0 .or. relerr >= 1.0d0) then
      error stop 'Relative error must be between 0 and 1'
    endif
    
    ! Validate batch parameters
    if (batch_size <= 0) then
      error stop 'Batch size must be positive'
    endif
    
    if (batch_size > ntestpart .and. batch_size /= 2000000000) then
      print *, 'Warning: batch_size exceeds ntestpart, batch mode disabled'
    endif
  end subroutine validate_configuration

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Calculate derived parameters
  subroutine calculate_derived_params(E_alpha, L1i)
    real(dp), intent(in) :: E_alpha
    integer, intent(in) :: L1i
    
    ! Calculate velocity and Larmor radius
    v0 = sqrt(2.d0*E_alpha*ev/(n_d*p_mass))
    rlarm = v0*n_d*p_mass*c/(n_e*e_charge)
    ro0 = rlarm
    
    ! Set inverse relativistic temperature (neglect relativistic effects)
    rmu = 1d8
    
    ! Calculate time parameters
    tau = trace_time*v0
    dtau = tau/dble(ntimstep-1)
    
    ! Field line integration parameters
    dphi = 2.d0*pi/(L1i*npoiper)
    dtaumin = 2.d0*pi*rbig/npoiper2
    ntau = ceiling(dtau/dtaumin)
    dtaumin = dtau/ntau
    
    ! Other derived parameters
    ntcut = ceiling(ntimstep*ntau*tcut/trace_time)
    norbper = ceiling(1d0*ntau*ntimstep/(L1i*npoiper2))
    nfp = L1i*norbper
    fper = 2d0*pi/dble(L1i)
    
    ! Initialize Lagrange interpolation parameters
    zerolam = 0.d0
    nplagr = 4
    nder = 0
    npl_half = nplagr/2
  end subroutine calculate_derived_params

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Load batch indices from file
  function load_batch_indices(batch_size) result(success)
    integer, intent(in) :: batch_size
    logical :: success
    integer :: iostat, i
    character, dimension(:), allocatable :: batch_file
    
    success = .false.
    
    allocate(batch_file(batch_size))
    allocate(idx(batch_size))
    
    open(1, file='batch.dat', recl=batch_size, iostat=iostat)
    if (iostat /= 0) then
      deallocate(batch_file)
      deallocate(idx)
      return
    endif
    
    do i = 1, batch_size
      read(1, iostat=iostat) batch_file(i)
      if (iostat /= 0) then
        close(1)
        deallocate(batch_file)
        deallocate(idx)
        return
      endif
      idx(i) = ichar(batch_file(i))
    end do
    
    close(1)
    deallocate(batch_file)
    success = .true.
  end function load_batch_indices

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Helper: Generate random batch indices
  subroutine generate_batch_indices(batch_size, ntestpart_in)
    integer, intent(in) :: batch_size, ntestpart_in
    integer :: i, n, iostat
    real :: ran_tmp
    
    allocate(idx(batch_size))
    
    call srand(ran_seed)
    
    ! Generate random indices, leaving out upper 1% as margin
    do i = 0, batch_size
      call random_number(ran_tmp)
      idx(i) = floor(ceiling((ntestpart_in - (ntestpart_in*0.01))) * ran_tmp)
    end do
    
    ! Sort and remove duplicates
    call sort_idx(idx, batch_size)
    
    ! Save to file
    open(1, file='batch.dat', recl=batch_size*2, iostat=iostat)
    if (iostat == 0) then
      do n = 1, batch_size
        write(1, *) idx(n)
      end do
      close(1)
    else
      print *, 'Warning: Could not save batch indices to file'
    endif
  end subroutine generate_batch_indices

  subroutine read_config(config_file)
    character(256), intent(in) :: config_file

    open(1, file=config_file, status='old', action='read')
    read(1, nml=config)
    close(1)

    call reset_seed_if_deterministic
    call validate_configuration
  end subroutine read_config


  subroutine params_init
    real(dp) :: E_alpha
    integer :: L1i

    ! Calculate alpha particle energy
    E_alpha = 3.5d6/facE_al

    ! Get vacuum chamber parameters
    call stevvo(RT0,R0i,L1i,cbfi,bz0i,bf0)
    rbig = rt0
    
    ! Calculate all derived parameters
    call calculate_derived_params(E_alpha, L1i)

    ! Initialize batch processing if needed
    call init_batch
    
    ! Allocate arrays based on configuration
    call reallocate_arrays
  end subroutine params_init

  subroutine init_batch
    logical :: old_batch, batch_loaded

    ! Check if batch processing is needed
    if (ntestpart <= batch_size) return
    
    if (reuse_batch) then
      ! Try to reuse existing batch file
      inquire(file="batch.dat", exist=old_batch)
      
      if (old_batch) then
        batch_loaded = load_batch_indices(batch_size)
        if (.not. batch_loaded) then
          print *, 'Warning: Could not load batch file, generating new indices'
          reuse_batch = .false.
        endif
      else
        ! No batch file found, generate new one
        reuse_batch = .false.
      endif
    endif
    
    if (.not. reuse_batch) then
      ! Generate new random batch indices
      call generate_batch_indices(batch_size, ntestpart)
    endif
    
    ! Update ntestpart to batch_size for the rest of the run
    ntestpart = batch_size
  end subroutine init_batch

  subroutine reallocate_arrays
    integer :: npoi
    npoi=nper*npoiper ! total number of starting points

    if(allocated(zstart))  deallocate(zstart)
    if(allocated(zend))  deallocate(zend)
    if(allocated(times_lost))  deallocate(times_lost)
    if(allocated(trap_par))  deallocate(trap_par)
    if(allocated(perp_inv))  deallocate(perp_inv)
    if(allocated(iclass))  deallocate(iclass)
    if(allocated(xstart))  deallocate(xstart)
    if(allocated(bstart))  deallocate(bstart)
    if(allocated(volstart))  deallocate(volstart)
    if(allocated(confpart_trap))  deallocate(confpart_trap)
    if(allocated(confpart_pass))  deallocate(confpart_pass)

    allocate(zstart(zstart_dim1,ntestpart), zend(zstart_dim1,ntestpart))
    allocate(times_lost(ntestpart), trap_par(ntestpart), perp_inv(ntestpart))
    allocate(xstart(3,npoi),bstart(npoi),volstart(npoi))
    allocate(confpart_trap(ntimstep),confpart_pass(ntimstep))
    allocate(iclass(3,ntestpart))
  end subroutine reallocate_arrays


  subroutine sort_idx(idx_arr, N)
    ! sort particle indices.
    integer :: N
    integer, dimension (N) :: idx_arr
    integer :: i, j, temp, o, r, num_removed, p
    logical :: swapped

    do j = n-1, 1, -1
       swapped = .false.
       do i = 1, j
          if (idx_arr(i) > idx_arr(i+1)) then
             temp = idx_arr(i)
             idx_arr(i) = idx_arr(i+1)
             idx_arr(i+1) = temp
             swapped = .true.
          end if
       end do
       if (.not. swapped) exit
    end do

    ! eliminate the duplicates, replace them by either appending from higher  indices, dependent on available indices.
    num_removed = 0
    do o=1, N
      if ((N-o) < (num_removed+1)) exit

      if (idx_arr(o) == idx_arr(o+1)) then
        num_removed = num_removed + 1
        do r=o+1, N-1
          idx_arr(r) = idx_arr(r+1)
        end do
      end if
    end do

    do p=N-num_removed+1, N
      idx_arr(p) = idx_arr(p-1) + 1
    end do

    if (idx_arr(N) > ntestpart) then
      print *,'ERROR - Invalid indices (Out of Range)!'
      error stop
    end if

  end subroutine sort_idx


  function should_skip(ipart)
    ! notrace_passing: no tracing of passing particles, assume that all are confined
    ! or skip strongly passing particles that are certainly confined
    logical :: should_skip
    integer, intent(in) :: ipart

    should_skip = (notrace_passing .eq. 1) .or. (trap_par(ipart) .le. contr_pp)
  end function should_skip


  subroutine reset_seed_if_deterministic
    ! for run with fixed random seed
    integer :: seedsize
    integer, allocatable :: seed(:)

    if (deterministic) then
      call random_seed(size = seedsize)
      if (.not. allocated(seed)) allocate(seed(seedsize))
      seed = 0
      call random_seed(put=seed)
    endif
  end subroutine reset_seed_if_deterministic

end module params
