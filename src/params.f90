module params
    use, intrinsic :: iso_fortran_env, only: int8, int64
    use util, only: pi, c, e_charge, p_mass, ev
    use parmot_mod, only: ro0, rmu
    use new_vmec_stuff_mod, only: old_axis_healing, old_axis_healing_boundary, &
                                  axis_healing_power_law, rho_axis_heal, &
                                  axis_healing, s_axis_heal, &
                                  axis_healing_polyfit_degree, &
                                  netcdffile, ns_s, ns_tp, multharm, vmec_B_scale, &
                                  vmec_RZ_scale
    use velo_mod, only: isw_field_type
    use magfie_sub, only: TEST
    use field_can_mod, only: eval_field => evaluate, field_can_t
    use orbit_symplectic_base, only: symplectic_integrator_t, multistage_integrator_t, &
                                     EXPL_IMPL_EULER, &
                                     boundary_event_fraction_tolerance, &
                                     boundary_event_radial_tolerance, &
                                     symplectic_newton_warning_mode
    use vmecin_sub, only: stevvo
    use callback, only: output_error, output_orbits_macrostep

    implicit none

    private :: config_value_is_finite

    ! Define real(dp) kind parameter
    integer, parameter :: dp = kind(1.0d0)
    integer          :: nper = 1000, ntestpart = 1024
    integer          :: zstart_dim1 = 5
    real(dp) :: dphi, phibeg = 0d0, bmod00, rlarm, bmax, bmin
    real(dp) :: tau, dtau, dtaumin, xi
    real(dp) :: RT0, R0i, cbfi, bz0i, bf0, rbig
    real(dp), dimension(1024) :: sbeg = 0.5d0
    real(dp) :: thetabeg = 0.0d0
    real(dp), dimension(:), allocatable :: bstart, volstart
    !where npoi is used, add a dimension at the end for sbeg
    real(dp), dimension(:, :), allocatable :: xstart
    real(dp), dimension(:, :), allocatable :: zstart, zend
    real(dp), dimension(:), allocatable :: confpart_trap, confpart_pass
    integer, dimension(:), allocatable :: unresolved_orbits
    real(dp), dimension(:), allocatable :: times_lost
    real(dp), dimension(:), allocatable :: boundary_event_radial_residual
    real(dp), dimension(:), allocatable :: boundary_event_time_width
    integer, dimension(:), allocatable :: orbit_exit_code
    integer, parameter :: ORBIT_EXIT_COMPLETED = 0
    integer, parameter :: ORBIT_EXIT_LCFS = 1
    integer, parameter :: ORBIT_EXIT_WALL = 2
    integer, parameter :: ORBIT_EXIT_SKIPPED = 3
    integer, parameter :: ORBIT_EXIT_NUMERICAL_DOMAIN = 101
    integer, parameter :: ORBIT_EXIT_NUMERICAL_MAXITER = 102
    integer, parameter :: ORBIT_EXIT_NUMERICAL_LINEAR = 103
    integer, parameter :: ORBIT_EXIT_NUMERICAL_EVENT = 104
    integer, parameter :: ORBIT_EXIT_NUMERICAL_FULL_ORBIT = 105
    real(dp) :: contr_pp = -1d0
    integer          :: ibins
    logical          :: generate_start_only = .False.
    integer          :: startmode = 1
    real(dp) :: grid_density = 0d0
    logical          :: special_ants_file = .False.

	    integer :: ntau ! number of dtaumin in dtau
	    integer(8) :: n_microsteps_total = 0_8
	    character(16) :: macrostep_time_grid = 'linear'
	    integer, allocatable :: ntau_macro(:)
	    integer(8), allocatable :: kt_macro(:)

    integer :: integmode = EXPL_IMPL_EULER

    ! Orbit model selector. 0 = guiding-center (GC), the default symplectic
    ! gyro-averaged path. 7 = full orbit (FO), the gyro-resolved Boris pusher in
    ! Cartesian on the Boozer/chartmap chart, the ASCOT-style counterpart to GC.
    ! Values 1-6 are reserved for other models.
    integer, parameter :: ORBIT_GC = 0
    integer, parameter :: ORBIT_FULL_ORBIT = 7
    integer :: orbit_model = ORBIT_GC

    ! Chart for the full-orbit field+geometry: full orbit currently supports only
    ! orbit_coord = 1 (Boozer/chartmap), which shares the production GC field.
    integer :: orbit_coord = 0

    integer :: kpart = 0 ! progress counter for particles

    real(dp) :: relerr = 1d-13

    ! Legacy progress/checkpoint interval, retained for input compatibility.
    ! Tracing performs no timed I/O; results and event totals are emitted after
    ! the OpenMP team joins.
    real(dp) :: checkpoint_interval = 10.0d0
    integer :: canonical_grid_nr = 62
    integer :: canonical_grid_ntheta = 63
    integer :: canonical_grid_nphi = 64
    real(dp) :: canonical_ode_relerr = 1d-11

    real(dp), allocatable :: trap_par(:), perp_inv(:)
    integer, allocatable :: iclass(:, :)
    logical, allocatable :: class_passing(:), class_lost(:)

    integer, parameter :: n_tip_vars = 6
    ! variables to evaluate at tip: z(1..5), par_inv
    integer :: nplagr, nder, npl_half
    integer :: norbper, nfp
    real(dp) :: fper, zerolam = 0d0

    real(dp) :: tcut = -1d0
    integer(8) :: ntcut
    integer :: nturns = 8
    logical          :: class_plot = .False.    !<=AAA
    real(dp) :: cut_in_per = 0d0        !<=AAA

    logical :: fast_class = .False.
    !if .true. quit immeadiately after fast classification

    ! Colliding with D-T reactor plasma. TODO: Make configurable
    logical :: swcoll = .False.
    real(dp) :: am1 = 2.0d0, am2 = 3.0d0, Z1 = 1.0d0, Z2 = 1.0d0, &
                densi1 = 0.5d14, densi2 = 0.5d14, tempi1 = 1.0d4, tempi2 = 1.0d4, &
                tempe = 1.0d4
    real(dp) :: dchichi, slowrate, dchichi_norm, slowrate_norm
    logical :: deterministic = .False.

    ! Further configuration parameters
    integer          :: notrace_passing = 0
    real(dp) :: facE_al = 1d0, trace_time = 1d-1
    integer :: ntimstep = 10000, npoiper = 100, npoiper2 = 256, n_e = 2
    real(dp) :: n_d = 4

    real(dp) :: v0

    logical :: debug = .False.
    logical :: output_results_netcdf = .False.
    logical :: restart = .False.
    integer :: ierr

    integer :: batch_size = 2000000000  ! Initialize large so batch mode is not default
    integer :: ran_seed = 42
    integer :: num_surf = 1
    logical :: reuse_batch = .False.
    integer, dimension(:), allocatable :: idx

    !> Warning mode may consume one isolated, fully exhausted numerical
    !> microstep. A second consecutive failure is terminal rather than freezing
    !> the remainder of the orbit at one phase-space point.
    integer, parameter :: max_consecutive_warning_holds = 1

    character(1000) :: field_input = ''
    character(1000) :: coord_input = ''
    character(1000) :: wall_input = ''
    character(16) :: wall_units = 'm'
    !> Radial threshold for querying an external STL wall. A negative value
    !> means that no external-wall query gate is active for this run.
    real(dp) :: wall_query_rho_lcfs = -1.0_dp
    character(16) :: chart_boundary_kind = 'auto'
    character(16) :: chart_boundary_kind_effective = 'lcfs'
    integer :: integ_coords = -1000  ! Sentinel: -1000 means user did not set it
    !> SPECTRE RK45 interface-crossing map: 1 = Level-1 refraction (default),
    !> 0 = Level-0 energy rescale (regression comparison).
    integer :: crossing_level = 1
    logical :: spectre_sbeg_is_toroidal_flux = .false.

    !> SPECTRE per-volume Meiss construction grid (n_r, n_th, n_phi). Each volume
    !> allocates rank-3 arrays plus batch splines of this size, so the
    !> total peak memory scales as Mvol*n_r*n_th*n_phi. spectre_ncon_phi = -1 means
    !> automatic: 32 for fields with toroidal harmonics and a minimal phi grid
    !> for axisymmetric fields (all n = 0),
    !> the dominant memory saver for tokamak cases. A positive value forces that
    !> phi count verbatim, bypassing the axisymmetric clamp.
    integer :: spectre_ncon_r = 48
    integer :: spectre_ncon_th = 48
    integer :: spectre_ncon_phi = -1
    integer :: spectre_ncon_order = 5
    integer :: spectre_ncon_ode_max_steps = 1000000
    real(dp) :: spectre_ncon_ode_relerr = 1.0e-2_dp

	    namelist /config/ notrace_passing, nper, npoiper, ntimstep, ntestpart, &
	        trace_time, num_surf, sbeg, phibeg, thetabeg, contr_pp, &
	        facE_al, npoiper2, n_e, n_d, netcdffile, ns_s, ns_tp, multharm, &
	        isw_field_type, generate_start_only, startmode, grid_density, &
	        special_ants_file, integmode, orbit_model, orbit_coord, relerr, &
	        symplectic_newton_warning_mode, &
	        boundary_event_fraction_tolerance, boundary_event_radial_tolerance, &
	        canonical_grid_nr, canonical_grid_ntheta, canonical_grid_nphi, &
	        canonical_ode_relerr, &
	        tcut, nturns, debug, &
	        class_plot, cut_in_per, fast_class, vmec_B_scale, &
	        vmec_RZ_scale, swcoll, deterministic, old_axis_healing, &
	        old_axis_healing_boundary, axis_healing_power_law, rho_axis_heal, &
	        axis_healing, s_axis_heal, axis_healing_polyfit_degree, &
	        am1, am2, Z1, Z2, &
	        densi1, densi2, tempi1, tempi2, tempe, &
	        batch_size, ran_seed, reuse_batch, field_input, coord_input, &
	        wall_input, wall_units, chart_boundary_kind, integ_coords, crossing_level, &
	        spectre_sbeg_is_toroidal_flux, spectre_ncon_r, spectre_ncon_th, &
	        spectre_ncon_phi, spectre_ncon_order, spectre_ncon_ode_max_steps, &
	        spectre_ncon_ode_relerr, &
        output_results_netcdf, &
	        output_error, output_orbits_macrostep, &  ! callback
	        macrostep_time_grid, checkpoint_interval, restart

    integer(int8), allocatable :: wall_hit(:)
    real(dp), allocatable :: wall_hit_cart(:, :)
    real(dp), allocatable :: wall_hit_normal_cart(:, :)
    real(dp), allocatable :: wall_hit_cos_incidence(:)
    real(dp), allocatable :: wall_hit_angle_rad(:)

contains

    subroutine read_config(config_file)
        character(256), intent(in) :: config_file

        open (1, file=config_file, status='old', action='read')
        read (1, nml=config)
        close (1)

        call apply_config_aliases

        call reset_seed_if_deterministic

        if (integmode > 0 .and. symplectic_newton_warning_mode) then
            print *, 'WARNING: symplectic integrators accept bounded Newton ', &
                'corrections after maxit; failed steps retry from the last ', &
                'valid state; see maxit/step-skip diagnostics'
        end if

        call validate_boundary_event_tolerances
        call validate_chart_boundary_kind
        if (min(canonical_grid_nr, canonical_grid_ntheta, &
            canonical_grid_nphi) < 6) then
            error stop 'canonical map grid dimensions must be at least 6'
        end if
        if (max(canonical_grid_nr, canonical_grid_ntheta, &
            canonical_grid_nphi) > 512) then
            error stop 'canonical map grid dimensions must not exceed 512'
        end if
        if (int(canonical_grid_nr, int64)*int(canonical_grid_ntheta, int64)* &
            int(canonical_grid_nphi, int64) > 2097152_int64) then
            error stop 'canonical map grid exceeds the 2097152-point limit'
        end if
        if (.not. config_value_is_finite(canonical_ode_relerr) .or. &
            canonical_ode_relerr <= 0d0) then
            error stop 'canonical_ode_relerr must be finite and positive'
        end if

        if (swcoll .and. (tcut > 0.0d0 .or. class_plot .or. fast_class)) then
            error stop 'Collisions are incompatible with classification'
        end if

    end subroutine read_config

    subroutine validate_chart_boundary_kind
        select case (trim(chart_boundary_kind))
        case ('auto', 'lcfs', 'wall', 'domain')
        case default
            error stop "chart_boundary_kind must be auto, lcfs, wall, or domain"
        end select
    end subroutine validate_chart_boundary_kind

    subroutine validate_boundary_event_tolerances
        if (.not. config_value_is_finite(boundary_event_fraction_tolerance) .or. &
            (boundary_event_fraction_tolerance /= -1d0 .and. &
             boundary_event_fraction_tolerance <= 0d0)) then
            error stop 'boundary_event_fraction_tolerance must be finite and positive or -1'
        end if
        if (.not. config_value_is_finite(boundary_event_radial_tolerance) .or. &
            (boundary_event_radial_tolerance /= -1d0 .and. &
             boundary_event_radial_tolerance <= 0d0)) then
            error stop 'boundary_event_radial_tolerance must be finite and positive or -1'
        end if
    end subroutine validate_boundary_event_tolerances

    pure logical function config_value_is_finite(value)
        real(dp), intent(in) :: value
        integer(int64), parameter :: exponent_mask = &
            int(z'7ff0000000000000', int64)
        integer(int64) :: bits

        bits = transfer(value, bits)
        config_value_is_finite = iand(bits, exponent_mask) /= exponent_mask
    end function config_value_is_finite

	    subroutine params_init
	        real(dp) :: E_alpha
	        integer :: L1i
	        real(dp) :: weight_sum, cumul_weight, w
	        integer :: i, nintv
	        integer(8) :: kt_target, kt_prev
	        character(16) :: grid_kind

	        if (isw_field_type == TEST) then
            ! TEST field uses normalized units: B0=1, R0=1, a=0.5
            ! Use normalized parameters for orbit integration
            v0 = 1.0d0          ! Normalized velocity
            rlarm = 1.0d0       ! Normalized Larmor radius (will be scaled by ro0)
            ro0 = 1.0d0         ! Normalized ro0 = mc/e * v0
            rmu = 1d8           ! Large inverse relativistic temperature
            tau = trace_time    ! Already in normalized units
	            dtau = tau/dble(ntimstep - 1)
            L1i = 1             ! One field period (full torus)
            rbig = 1.0d0        ! Major radius R0=1
            dphi = 2.d0*pi/(L1i*npoiper)
	            dtaumin = 2.d0*pi*rbig/npoiper2
	            ntau = ceiling(dtau/dtaumin)
	            dtaumin = dtau/ntau
	            fper = 2d0*pi       ! Full torus
	        else
            E_alpha = 3.5d6/facE_al
            ! set alpha energy, velocity, and Larmor radius
            v0 = sqrt(2.d0*E_alpha*ev/(n_d*p_mass))
            rlarm = v0*n_d*p_mass*c/(n_e*e_charge)
            ro0 = rlarm
            ! Neglect relativistic effects by large inverse relativistic temperature
            rmu = 1d8
            ! normalized slowing down time:
            tau = trace_time*v0
	            ! normalized time step:
	            dtau = tau/dble(ntimstep - 1)
            ! parameters for the vacuum chamber:
            call stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0)
            rbig = rt0
            ! field line integration step step over phi (to check chamber wall crossing)
            dphi = 2.d0*pi/(L1i*npoiper)
	            ! orbit integration time step (to check chamber wall crossing)
	            dtaumin = 2.d0*pi*rbig/npoiper2
	            ntau = ceiling(dtau/dtaumin)
	            dtaumin = dtau/ntau
	            fper = 2d0*pi/dble(L1i)
	        end if

	        ! Macrostep schedule (number of microsteps per macrostep).
	        ! - Default is linear: constant ntau per macrostep.
	        ! - Log schedule: macrosteps are distributed logarithmically in time
	        !   while keeping the microstep resolution dtaumin. The total number
	        !   of microsteps is preserved as (ntimstep-1)*ntau.
	        nintv = max(1, ntimstep - 1)
	        n_microsteps_total = int(nintv, kind=8) * int(ntau, kind=8)
	        if (allocated(ntau_macro)) deallocate(ntau_macro)
	        if (allocated(kt_macro)) deallocate(kt_macro)
	        allocate (ntau_macro(ntimstep))
	        allocate (kt_macro(ntimstep))
	        ntau_macro = 0
	        kt_macro = 0_8

	        grid_kind = to_lower(macrostep_time_grid)
	        if (trim(grid_kind) == 'log') then
	            ! Logarithmic macrostep grid: timesteps appear equally spaced
	            ! on a log10 plot of time. Weight for interval i is 10^((i-1)/(N-1)),
	            ! giving a 10:1 ratio between last and first intervals.
	            ! Uses cumulative rounding to guarantee exact total microsteps.

	            ! First pass: compute total weight
	            weight_sum = 0.0d0
	            do i = 1, nintv
	                weight_sum = weight_sum + 10.0d0 ** (dble(i - 1) / dble(max(1, nintv - 1)))
	            end do

	            ! Second pass: distribute microsteps via cumulative rounding
	            ! Guarantee at least 1 microstep per interval to avoid empty loops
	            cumul_weight = 0.0d0
	            kt_prev = 0_8
	            do i = 1, nintv
	                w = 10.0d0 ** (dble(i - 1) / dble(max(1, nintv - 1)))
	                cumul_weight = cumul_weight + w
	                kt_target = nint(cumul_weight / weight_sum * dble(n_microsteps_total), kind=8)
	                ntau_macro(i + 1) = max(1, int(kt_target - kt_prev))
	                kt_macro(i + 1) = kt_prev + int(ntau_macro(i + 1), kind=8)
	                kt_prev = kt_macro(i + 1)
	            end do
	        else
	            do i = 2, ntimstep
	                ntau_macro(i) = ntau
	                kt_macro(i) = kt_macro(i - 1) + int(ntau, kind=8)
	            end do
	        end if

	        ntcut = microstep_cut_index(ntimstep, ntau, tcut, trace_time)
	        norbper = ceiling(1d0*ntau*ntimstep/(L1i*npoiper2))
	        nfp = L1i*norbper

        zerolam = 0.d0
        nplagr = 4
        nder = 0
        npl_half = nplagr/2

        call init_batch
        call reallocate_arrays
	    end subroutine params_init

    pure function microstep_cut_index(ntimstep_in, ntau_in, tcut_in, &
            trace_time_in) result(ntcut_out)
        ! Microstep index of the classification cut time tcut; <=0 disables it.
        ! ntimstep*ntau reaches ~1e10 for second-scale traces and overflows a
        ! 32-bit product, flipping the sign and spuriously enabling classification,
        ! so evaluate in real(dp) and return int64.
        integer, intent(in) :: ntimstep_in, ntau_in
        real(dp), intent(in) :: tcut_in, trace_time_in
        integer(8) :: ntcut_out

        if (tcut_in > 0d0 .and. trace_time_in > 0d0) then
            ntcut_out = ceiling( &
                real(ntimstep_in, dp)*real(ntau_in, dp)*tcut_in/trace_time_in, &
                kind=8)
        else
            ntcut_out = -1_8
        end if
    end function microstep_cut_index

	    pure function to_lower(s) result(out)
	        character(*), intent(in) :: s
	        character(len(s)) :: out
	        integer :: i, c
	        out = s
	        do i = 1, len(s)
	            c = iachar(out(i:i))
	            if (c >= iachar('A') .and. c <= iachar('Z')) then
	                out(i:i) = achar(c + 32)
	            end if
	        end do
	    end function to_lower

    subroutine init_batch
        integer :: iostat, i, n
        real :: ran_tmp
        character, dimension(:), allocatable :: batch_file
        logical :: old_batch

        !See if batch is wanted, if yes find random or re-use from previous file.
        if (ntestpart > batch_size) then
            if (reuse_batch) then
                inquire (file="batch.dat", exist=old_batch)

                if (old_batch) then
                    allocate (batch_file(batch_size))
                    allocate (idx(batch_size))
                    open (1, file='batch.dat', recl=batch_size, iostat=iostat)
                    if (iostat /= 0) goto 666

                    do i = 1, batch_size
                        read (1, iostat=iostat) batch_file(i)
                        if (iostat /= 0) goto 667
                        idx(i) = ichar(batch_file(i))
                    end do
                    deallocate (batch_file)
                    close (1)
                else !old_batch
                    ! no old batch file found, so we pretend the user knew this and
                    ! set the flag to .False.
                    reuse_batch = .False.
                end if !old_batch
            end if !reuse_batch

            if (.not. reuse_batch) then
                !Create a random list(batch_size) of indices using ran_seed.
                allocate (idx(batch_size))

                call srand(ran_seed)

                do i = 1, batch_size
                    call random_number(ran_tmp)
                    ! Create randomized indices from the amount available, leaving out
                    ! the upper 1% as margin for later sorting and replacing of
                    ! duplicates (relevant for smaller batches).
                    idx(i) = floor(ceiling((ntestpart - (ntestpart*0.01)))*ran_tmp)
                end do

                call sort_idx(idx, batch_size)

                open (1, file='batch.dat', recl=batch_size*2, iostat=iostat)

                do n = 1, batch_size
                    write (1, *) idx(n)
                end do
                close (1)
            end if !reuse_batch again

            !Set ntestpart to batch_size for rest of the run.
            ntestpart = batch_size
        end if !batches wanted

        return

666     print *, iostat
        error stop
667     print *, iostat
        error stop
    end subroutine init_batch

    subroutine reallocate_arrays
        integer :: npoi
        npoi = nper*npoiper ! total number of starting points

        if (allocated(zstart)) deallocate (zstart)
        if (allocated(zend)) deallocate (zend)
        if (allocated(times_lost)) deallocate (times_lost)
        if (allocated(boundary_event_radial_residual)) &
            deallocate (boundary_event_radial_residual)
        if (allocated(boundary_event_time_width)) deallocate (boundary_event_time_width)
        if (allocated(orbit_exit_code)) deallocate (orbit_exit_code)
        if (allocated(trap_par)) deallocate (trap_par)
        if (allocated(perp_inv)) deallocate (perp_inv)
        if (allocated(iclass)) deallocate (iclass)
        if (allocated(class_passing)) deallocate (class_passing)
        if (allocated(class_lost)) deallocate (class_lost)
        if (allocated(xstart)) deallocate (xstart)
        if (allocated(bstart)) deallocate (bstart)
        if (allocated(volstart)) deallocate (volstart)
        if (allocated(confpart_trap)) deallocate (confpart_trap)
        if (allocated(confpart_pass)) deallocate (confpart_pass)
        if (allocated(unresolved_orbits)) deallocate (unresolved_orbits)
        if (allocated(wall_hit)) deallocate (wall_hit)
        if (allocated(wall_hit_cart)) deallocate (wall_hit_cart)
        if (allocated(wall_hit_normal_cart)) deallocate (wall_hit_normal_cart)
        if (allocated(wall_hit_cos_incidence)) deallocate (wall_hit_cos_incidence)
        if (allocated(wall_hit_angle_rad)) deallocate (wall_hit_angle_rad)

        allocate (zstart(zstart_dim1, ntestpart), zend(zstart_dim1, ntestpart))
        allocate (times_lost(ntestpart), trap_par(ntestpart), perp_inv(ntestpart))
        allocate (boundary_event_radial_residual(ntestpart), &
                  boundary_event_time_width(ntestpart))
        allocate (orbit_exit_code(ntestpart))
        allocate (xstart(3, npoi), bstart(npoi), volstart(npoi))
        allocate (confpart_trap(ntimstep), confpart_pass(ntimstep))
        allocate (unresolved_orbits(ntimstep))
        allocate (iclass(3, ntestpart))
        allocate (class_passing(ntestpart), class_lost(ntestpart))
        allocate (wall_hit(ntestpart))
        allocate (wall_hit_cart(3, ntestpart))
        allocate (wall_hit_normal_cart(3, ntestpart))
        allocate (wall_hit_cos_incidence(ntestpart))
        allocate (wall_hit_angle_rad(ntestpart))

        times_lost = 0.0d0
        boundary_event_radial_residual = -1d0
        boundary_event_time_width = -1d0
        orbit_exit_code = 0
        trap_par = 0.0d0
        perp_inv = 0.0d0
        iclass = 0
        class_passing = .false.
        class_lost = .false.
        wall_hit = 0_int8
        wall_hit_cart = 0.0d0
        wall_hit_normal_cart = 0.0d0
        wall_hit_cos_incidence = 0.0d0
        wall_hit_angle_rad = 0.0d0
    end subroutine reallocate_arrays

    subroutine sort_idx(idx_arr, N)
        ! sort particle indices.
        integer :: N
        integer, dimension(N) :: idx_arr
        integer :: i, j, temp, o, r, num_removed, p
        logical :: swapped

        do j = n - 1, 1, -1
            swapped = .false.
            do i = 1, j
                if (idx_arr(i) > idx_arr(i + 1)) then
                    temp = idx_arr(i)
                    idx_arr(i) = idx_arr(i + 1)
                    idx_arr(i + 1) = temp
                    swapped = .true.
                end if
            end do
            if (.not. swapped) exit
        end do

        ! eliminate the duplicates, replace them by either appending from higher
        ! indices, dependent on available indices.
        num_removed = 0
        do o = 1, N
            if ((N - o) < (num_removed + 1)) exit

            if (idx_arr(o) == idx_arr(o + 1)) then
                num_removed = num_removed + 1
                do r = o + 1, N - 1
                    idx_arr(r) = idx_arr(r + 1)
                end do
            end if
        end do

        do p = N - num_removed + 1, N
            idx_arr(p) = idx_arr(p - 1) + 1
        end do

        if (idx_arr(N) > ntestpart) then
            print *, 'ERROR - Invalid indices (Out of Range)!'
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

    subroutine apply_config_aliases
        ! Handle deprecated aliases and apply fallback logic for config parameters.
        ! Priority: new name > deprecated alias > default
        !
        ! Aliases:
        !   netcdffile -> fallback for field_input and coord_input (deprecated)
        !   isw_field_type -> integ_coords (deprecated)
        !
        ! Fallback chain:
        !   field_input: explicit > netcdffile > ''
        !   coord_input: explicit > netcdffile > field_input > ''

        ! netcdffile serves as fallback for both field_input and coord_input
        if (field_input == '' .and. len_trim(netcdffile) > 0) then
            field_input = netcdffile
        end if

        if (coord_input == '' .and. len_trim(netcdffile) > 0) then
            coord_input = netcdffile
        end if

        ! coord_input falls back to field_input if still not set
        if (coord_input == '') then
            coord_input = field_input
        end if

        ! Do not force netcdffile to follow coord_input:
        ! - coord_input may point to a chartmap file (reference coordinates)
        ! - netcdffile is still used as the VMEC equilibrium file by legacy code paths
        ! The VMEC file selection is handled in simple_main.init_field.

        ! isw_field_type is deprecated alias for integ_coords
        ! integ_coords == -1000 means user did not set it
        if (integ_coords == -1000) then
            integ_coords = isw_field_type
        else
            isw_field_type = integ_coords
        end if
    end subroutine apply_config_aliases

    subroutine reset_seed_if_deterministic
        ! for run with fixed random seed
        integer :: seedsize
        integer, allocatable :: seed(:)

        if (deterministic) then
            call random_seed(size=seedsize)
            if (.not. allocated(seed)) allocate (seed(seedsize))
            seed = ran_seed
            call random_seed(put=seed)
        end if
    end subroutine reset_seed_if_deterministic

end module params
