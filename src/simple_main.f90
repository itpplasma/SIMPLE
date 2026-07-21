module simple_main
    use, intrinsic :: iso_fortran_env, only: int8
    use omp_lib
    use util, only: sqrt2, twopi
    use simple, only: init_vmec, init_sympl, init_fo, orbit_timestep_fo, &
        orbit_timestep_fo_bridge, reseed_sympl, tracer_t, ORBIT_FO_LOSS, &
        ORBIT_FO_NUMERICAL
    use diag_mod, only: icounter
    use collis_alp, only: loacol_alpha, stost, init_collision_profiles
    use samplers, only: sample
    use field_can_mod, only: field_can_t, integ_to_ref, ref_to_integ, &
        init_field_can
    use callback, only: output_orbits_macrostep
    use params, only: swcoll, ntestpart, startmode, special_ants_file, num_surf, &
        grid_density, dtau, dtaumin, ntau, v0, kpart, confpart_pass, &
        confpart_trap, unresolved_orbits, times_lost, orbit_exit_code, &
        boundary_event_radial_residual, &
        boundary_event_time_width, integmode, relerr, trace_time, class_plot, &
        fast_class, generate_start_only, ntcut, iclass, bmin, bmax, zstart, &
        zend, trap_par, perp_inv, sbeg, ntimstep, should_skip, &
        reset_seed_if_deterministic, field_input, isw_field_type, reuse_batch, &
        max_consecutive_warning_holds, &
        coord_input, wall_input, wall_units, wall_hit, wall_hit_cart, &
        wall_query_rho_lcfs, &
        chart_boundary_kind, chart_boundary_kind_effective, &
        wall_hit_normal_cart, wall_hit_cos_incidence, wall_hit_angle_rad, &
        ntau_macro, kt_macro, checkpoint_interval, orbit_model, orbit_coord, &
        ORBIT_GC, ORBIT_FULL_ORBIT, ORBIT_EXIT_COMPLETED, ORBIT_EXIT_LCFS, &
        ORBIT_EXIT_WALL, ORBIT_EXIT_SKIPPED, ORBIT_EXIT_NUMERICAL_DOMAIN, &
        ORBIT_EXIT_NUMERICAL_MAXITER, ORBIT_EXIT_NUMERICAL_LINEAR, &
        ORBIT_EXIT_NUMERICAL_EVENT, ORBIT_EXIT_NUMERICAL_FULL_ORBIT
    use params, only: canonical_grid_nr, canonical_grid_ntheta, &
        canonical_grid_nphi, canonical_ode_relerr
    use diag_counters, only: diag_counters_init, count_event, &
        EVT_WARNING_STEP_SKIP, EVT_SYMPLECTIC_RK_RECOVERY, &
        EVT_SYMPLECTIC_RESUME, EVT_TOROIDAL_REGULARIZATION, &
        EVT_AXIS_FO_BRIDGE
    use progress_monitor, only: progress_init, progress_tick, progress_finalize
    use restart_mod, only: particle_done, read_restart_data, restore_confined_counts
    use chartmap_metadata, only: chartmap_metadata_t, read_chartmap_metadata
    use reference_coordinates, only: ref_coords
    use stl_wall_intersection, only: stl_wall_t, stl_wall_init, &
        stl_wall_finalize, &
        stl_wall_first_hit_segment_with_normal
    use libneo_coordinates, only: chartmap_coordinate_system_t
    use orbit_symplectic_base, only: symplectic_integrator_t, &
        SYMPLECTIC_STEP_BOUNDARY, &
        SYMPLECTIC_STEP_OUTSIDE_DOMAIN, SYMPLECTIC_STEP_MAXITER, &
        SYMPLECTIC_STEP_LINEAR_SOLVE, SYMPLECTIC_STEP_EVENT_NOT_CONVERGED, &
        SYMPLECTIC_STEP_BOUNDARY_LIMITED, symplectic_newton_warning_mode

    implicit none

    ! Define real(dp) kind parameter
    integer, parameter :: dp = kind(1.0d0)

    logical, save :: wall_enabled = .false.
    integer, save :: map_boundary_exit_code = ORBIT_EXIT_LCFS
    real(dp), save :: chartmap_cart_scale_to_m = -1.0d0
    type(stl_wall_t), save :: wall

    integer, parameter :: RK_RECOVERY_GENERIC = 1
    integer, parameter :: RK_RECOVERY_AXIS = 2
    integer, parameter :: RK_RECOVERY_EDGE = 3
    real(dp), parameter :: RK_AXIS_ENTER = 0.1_dp
    real(dp), parameter :: RK_AXIS_EXIT = 0.12_dp
    real(dp), parameter :: FO_AXIS_ENTER = 0.02_dp
    real(dp), parameter :: RK_EDGE_ENTER = 0.98_dp
    integer, parameter :: RK_RECOVERY_STEP_LIMIT = 10000
    integer, parameter :: TOROIDAL_RECOVERY_STEP_LIMIT = 100000
    integer, parameter :: RK_SAFE_STREAK_TO_RESUME = 256
    real(dp), parameter :: RK_RECOVERY_AXIS_ENTER = 0.1_dp
    ! Retry only a failed legacy axis solve with a smoother, still axis-local
    ! Cartesian regularization. Successful release-era RK paths stay unchanged.
    real(dp), parameter :: RK_RECOVERY_AXIS_FLOOR = 1.0e-3_dp
    real(dp), parameter :: FO_BRIDGE_BSTAR_ENTER = 0.2_dp
    real(dp), parameter :: RK_RECOVERY_BSTAR_EXIT = 0.3_dp
    integer, parameter :: RECOVERY_STEP_STANDARD = 0
    integer, parameter :: RECOVERY_STEP_AXIS = 1
    integer, parameter :: RECOVERY_STEP_TOROIDAL = 2
    integer, parameter :: RECOVERY_STEP_FULL_ORBIT = 3

    type :: rk_recovery_state_t
        logical :: active = .false.
        logical :: has_momentum_reference = .false.
        integer :: reason = RK_RECOVERY_GENERIC
        integer :: steps = 0
        integer :: safe_steps = 0
        integer :: failures = 0
        real(dp) :: momentum_reference = 0.0_dp
    end type rk_recovery_state_t

contains

    subroutine orbit_timestep_recovery(fo, z, interval, tolerance, ierr, method)
        use alpha_lifetime_sub, only: orbit_timestep_axis, &
            bstar_parallel_factor, orbit_timestep_toroidal_regularized
        use orbit_fo_boris, only: fo_state_t

        type(fo_state_t), intent(inout) :: fo
        real(dp), intent(inout) :: z(5)
        real(dp), intent(in) :: interval, tolerance
        integer, intent(out) :: ierr
        integer, intent(out) :: method
        real(dp) :: z_start(5)
        real(dp) :: bstar_factor
        real(dp) :: x_reference(3)

        method = RECOVERY_STEP_STANDARD
        z_start = z
        call orbit_timestep_axis(z, interval, interval, tolerance, ierr, &
            RK_RECOVERY_STEP_LIMIT)
        if (ierr /= 2) return

        call bstar_parallel_factor(z_start, bstar_factor)
        call integ_to_ref(z_start(1:3), x_reference)
        if (x_reference(1) < RK_RECOVERY_AXIS_ENTER .and. &
                abs(bstar_factor) >= FO_BRIDGE_BSTAR_ENTER) then
            z = z_start
            call orbit_timestep_axis(z, interval, interval, tolerance, ierr, &
                RK_RECOVERY_STEP_LIMIT, RK_RECOVERY_AXIS_FLOOR, &
                RK_RECOVERY_AXIS_ENTER)
            if (ierr /= 2) then
                method = RECOVERY_STEP_AXIS
                return
            end if
        end if

        z = z_start
        call orbit_timestep_toroidal_regularized(z, interval, &
            max(tolerance, 1.0e-10_dp), ierr, &
            TOROIDAL_RECOVERY_STEP_LIMIT)
        if (ierr == 0) then
            method = RECOVERY_STEP_TOROIDAL
            call count_event(EVT_TOROIDAL_REGULARIZATION)
            return
        end if
        if (ierr /= 2 .or. x_reference(1) >= FO_AXIS_ENTER) return

        z = z_start
        call orbit_timestep_fo_bridge(fo, z, interval, ierr)
        if (ierr == 0) then
            method = RECOVERY_STEP_FULL_ORBIT
            call count_event(EVT_AXIS_FO_BRIDGE)
        end if
    end subroutine orbit_timestep_recovery

    subroutine activate_rk_recovery(state, z, momentum_reference)
        type(rk_recovery_state_t), intent(inout) :: state
        real(dp), intent(in) :: z(5), momentum_reference
        real(dp) :: x_reference(3)

        call integ_to_ref(z(1:3), x_reference)
        state%active = .true.
        state%has_momentum_reference = .true.
        state%steps = 0
        state%safe_steps = 0
        state%momentum_reference = momentum_reference
        state%failures = min(state%failures + 1, 9)
        if (x_reference(1) < RK_AXIS_ENTER) then
            state%reason = RK_RECOVERY_AXIS
        else if (x_reference(1) > RK_EDGE_ENTER) then
            state%reason = RK_RECOVERY_EDGE
        else
            state%reason = RK_RECOVERY_GENERIC
        end if
    end subroutine activate_rk_recovery

    logical function should_resume_symplectic(state, z)
        type(rk_recovery_state_t), intent(in) :: state
        real(dp), intent(in) :: z(5)
        integer :: cooldown
        real(dp) :: x_reference(3)

        call integ_to_ref(z(1:3), x_reference)
        select case (state%reason)
        case (RK_RECOVERY_AXIS)
            should_resume_symplectic = x_reference(1) >= RK_AXIS_EXIT .and. &
                state%safe_steps >= RK_SAFE_STREAK_TO_RESUME
        case (RK_RECOVERY_EDGE)
            should_resume_symplectic = x_reference(1) <= RK_EDGE_ENTER .and. &
                state%safe_steps >= RK_SAFE_STREAK_TO_RESUME
        case default
            cooldown = 2**(7 + state%failures)
            should_resume_symplectic = state%safe_steps >= cooldown
        end select
    end function should_resume_symplectic

    subroutine update_recovery_streak(state, z, method)
        use alpha_lifetime_sub, only: bstar_parallel_factor

        type(rk_recovery_state_t), intent(inout) :: state
        real(dp), intent(in) :: z(5)
        integer, intent(in) :: method
        real(dp) :: bstar_factor

        state%steps = state%steps + 1
        if (method /= RECOVERY_STEP_STANDARD) then
            state%safe_steps = 0
            return
        end if
        call bstar_parallel_factor(z, bstar_factor)
        if (abs(bstar_factor) < RK_RECOVERY_BSTAR_EXIT) then
            state%safe_steps = 0
            return
        end if
        state%safe_steps = state%safe_steps + 1
    end subroutine update_recovery_streak

    subroutine main
        use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, &
            integmode, params_init, swcoll, generate_start_only, &
            isw_field_type, field_input, startmode, &
            ntestpart, ntimstep, coord_input, restart
        use timing, only: init_timer, print_phase_time
        use magfie_sub, only: TEST, VMEC, SPECTRE, init_magfie
        use samplers, only: init_starting_surf, sample_spectre_surface, &
            init_spectre_start_bounds
        use version, only: simple_version
        use field_boozer_chartmap, only: is_boozer_chartmap
        use field, only: is_spectre_file

        implicit none

        character(256) :: config_file
        type(tracer_t) :: norb
        logical :: chartmap_mode, spectre_mode

        ! Print version on startup
        print '(A,A)', 'SIMPLE version ', simple_version

        ! Initialize timing
        call init_timer()

        ! read configuration file name from command line arguments
        if (command_argument_count() == 0) then
            config_file = 'simple.in'
        else
            call get_command_argument(1, config_file)
        end if
        call print_phase_time('Command line parsing completed')

        ! Must be called in this order. TODO: Fix
        call read_config(config_file)
        call validate_orbit_model_config
        call print_phase_time('Configuration reading completed')

        call read_profiles_config(config_file)
        call print_phase_time('Profiles configuration reading completed')

        call init_field(norb, netcdffile, ns_s, ns_tp, multharm, integmode)
        call print_phase_time('Field initialization completed')

        call params_init
        call print_phase_time('Parameter initialization completed')

        call print_parameters
        call print_phase_time('Parameter printing completed')

        if (swcoll) then
            call init_collisions
            call print_phase_time('Collision initialization completed')
        end if

        chartmap_mode = .false.
        spectre_mode = .false.
        if (len_trim(field_input) > 0) then
            spectre_mode = is_spectre_file(field_input)
            if (.not. spectre_mode) chartmap_mode = is_boozer_chartmap(field_input)
        end if

        if (isw_field_type == TEST) then
            ! TEST field uses analytic tokamak - no VMEC needed for sampling
            call init_magfie(TEST)
            call print_phase_time('TEST field initialization completed')
            call sample_particles_test_field
            call print_phase_time('TEST field particle sampling completed')
        else if (spectre_mode) then
            ! SPECTRE: VMEC-free, per-volume RK45 guiding centers. Sampling and
            ! tracing use magfie_spectre directly on the stacked-rho chart.
            call init_magfie(SPECTRE)
            call print_phase_time('SPECTRE magfie initialization completed')

            if (startmode == 1) then
                call sample_spectre_surface(zstart)
            else
                call sample_particles(.true.)
                if (startmode == 2) call init_spectre_start_bounds(zstart)
            end if
            call print_phase_time('SPECTRE particle sampling completed')

            if (generate_start_only) stop 'stopping after generating start.dat'
        else if (chartmap_mode) then
            ! Boozer chartmap: no VMEC, use Boozer magfie for field line tracing.
            call init_magfie(isw_field_type)
            call print_phase_time('Boozer magfie initialization completed')

            call init_starting_surf
            call print_phase_time('Starting surface initialization completed')

            call sample_particles(.true.)
            call print_phase_time('Particle sampling completed')

            if (generate_start_only) stop 'stopping after generating start.dat'
        else
            call init_magfie(VMEC)
            call print_phase_time('VMEC magnetic field initialization completed')

            call init_starting_surf
            call print_phase_time('Starting surface initialization completed')

            call sample_particles
            call print_phase_time('Particle sampling completed')

            if (generate_start_only) stop 'stopping after generating start.dat'

            call init_magfie(isw_field_type)
            call print_phase_time('Field type initialization completed')
        end if

        if ((isw_field_type /= TEST) .and. needs_bminmax_cache()) then
            call init_bminmax
            call print_phase_time('Bmin/Bmax initialization completed')
        end if

        call init_counters
        call print_phase_time('Counter initialization completed')

        if (restart) then
            call read_restart_data()
            call restore_confined_counts()
            call print_phase_time('Restart data loaded')
        end if

        block
            character(32) :: gpu_bench_env
            integer :: gpu_bench_len, gpu_bench_stat
            call get_environment_variable('SIMPLE_GPU_BENCH', gpu_bench_env, &
                gpu_bench_len, gpu_bench_stat)
            if (gpu_bench_stat == 0 .and. gpu_bench_len > 0) then
                call trace_compare_gpu(norb)
                return
            end if
        end block

        call diag_counters_init
        if (isw_field_type == SPECTRE) then
            block
                use interface_crossing, only: crossing_log_reset
                use spectre_sympl_orbit, only: sympl_landing_stats_reset
                call crossing_log_reset(256)
                call sympl_landing_stats_reset
            end block
        end if
        call progress_init(checkpoint_interval, ntestpart)
        call trace_parallel(norb)
        call progress_finalize
        call print_phase_time('Parallel particle tracing completed')

        call write_output
        call print_phase_time('Output writing completed')

        call stl_wall_finalize(wall)
    end subroutine main

    ! Reject orbit_model values this build does not implement, with a clear message,
    ! before any tracing starts. Only guiding-center (the default symplectic path)
    ! and full orbit (the gyro-resolved Boris pusher) are available here.
    subroutine validate_orbit_model_config
        select case (orbit_model)
        case (ORBIT_GC)
            continue
        case (ORBIT_FULL_ORBIT)
            if (orbit_coord /= 1) error stop &
                'orbit_model=ORBIT_FULL_ORBIT supports only orbit_coord=1 (Boozer)'
            if (class_plot .or. fast_class) error stop &
                'orbit_model=ORBIT_FULL_ORBIT does not support orbit classification'
        case default
            error stop 'unsupported orbit_model (use 0 = guiding-center or '// &
                '7 = full orbit)'
        end select
    end subroutine validate_orbit_model_config

    subroutine init_field(self, vmec_file, ans_s, ans_tp, amultharm, aintegmode)
        use field_base, only: magnetic_field_t
        use field, only: field_from_file, is_spectre_file
        use field_boozer_chartmap, only: boozer_chartmap_field_t, is_boozer_chartmap
        use timing, only: print_phase_time
        use magfie_sub, only: TEST, CANFLUX, VMEC, BOOZER, MEISS, ALBERT, &
            REFCOORDS, SPECTRE, set_magfie_refcoords_field
        use field_splined, only: splined_field_t, create_splined_field
        use field_vmec, only: vmec_field_t
        use reference_coordinates, only: init_reference_coordinates, ref_coords
        use params, only: coord_input, field_input, wall_input, wall_units

        character(*), intent(in) :: vmec_file
        type(tracer_t), intent(inout) :: self
        integer, intent(in) :: ans_s, ans_tp, amultharm, aintegmode
        class(magnetic_field_t), allocatable :: field_temp
        character(:), allocatable :: vmec_equilibrium_file
        logical :: use_boozer_chartmap, use_spectre

        self%integmode = aintegmode
        map_boundary_exit_code = ORBIT_EXIT_LCFS
        chart_boundary_kind_effective = 'lcfs'

        ! Check if field_input is a Boozer chartmap or SPECTRE file (no VMEC needed)
        use_boozer_chartmap = .false.
        use_spectre = .false.
        if (len_trim(field_input) > 0) then
            use_spectre = is_spectre_file(field_input)
            if (.not. use_spectre) use_boozer_chartmap = is_boozer_chartmap(field_input)
        end if

        ! TEST field is analytic - no VMEC or field files needed
        if (isw_field_type == TEST) then
            self%fper = twopi ! Full torus for analytic tokamak
            call print_phase_time('TEST field mode - no input files required')
        else if (use_spectre) then
            call init_spectre_field(self)
            call print_phase_time('SPECTRE field loading completed')
        else if (use_boozer_chartmap) then
            ! Boozer chartmap: file-based, no VMEC initialization needed
            call init_reference_coordinates(coord_input)
            call print_phase_time('Reference coordinate system '// &
                'initialization completed')

            block
                use libneo_coordinates, only: chartmap_coordinate_system_t
                select type (cs => ref_coords)
                class is (chartmap_coordinate_system_t)
                    self%fper = twopi / real(cs%num_field_periods, dp)
                class default
                    self%fper = twopi
                end select
            end block

            call init_stl_wall_if_enabled(coord_input)
            call print_phase_time('STL wall initialization completed')

            if (self%integmode >= 0) then
                call field_from_file(field_input, field_temp)
                call print_phase_time('Boozer chartmap field loading completed')
            end if
        else
            vmec_equilibrium_file = select_vmec_equilibrium_file(vmec_file, &
                field_input, &
                coord_input)
            call init_vmec(vmec_equilibrium_file, ans_s, ans_tp, amultharm, &
                self%fper)
            call print_phase_time('VMEC initialization completed')

            call init_reference_coordinates(coord_input)
            call print_phase_time('Reference coordinate system '// &
                'initialization completed')

            call init_stl_wall_if_enabled(coord_input)
            call print_phase_time('STL wall initialization completed')

            if (self%integmode >= 0) then
                if (trim(field_input) == '') then
                    print *, 'simple_main.init_field: field_input must be set (see ', &
                        'params.apply_config_aliases)'
                    error stop
                end if

                call field_from_file(field_input, field_temp)
                call print_phase_time('Field from file loading completed')

                if (isw_field_type == REFCOORDS) then
                    select type (field_temp)
                        type is (splined_field_t)
                        call set_magfie_refcoords_field(field_temp)
                        type is (vmec_field_t)
                        block
                            type(splined_field_t) :: splined_vmec
                            call create_splined_field(field_temp, ref_coords, &
                                splined_vmec)
                            call set_magfie_refcoords_field(splined_vmec)
                        end block
                    class default
                        print *, &
                            'simple_main.init_field: REFCOORDS requires '// &
                            'a splined field'
                        print *, &
                            'Supported inputs: coils (auto-splined) or VMEC wout', &
                            ' (splined onto coord_input)'
                        error stop
                    end select
                end if
            end if
        end if

        if (self%integmode > 0) then
            select case (isw_field_type)
            case (VMEC)
                error stop &
                    'Symplectic guiding-center integrators require '// &
                    'canonical field representation (set isw_field_type to TEST, '// &
                    'CANFLUX, BOOZER, MEISS, ALBERT, or SPECTRE)'
            case (TEST, CANFLUX, BOOZER, MEISS, ALBERT, SPECTRE)
                continue
            case default
                error stop &
                    'Unknown canonical field type for symplectic guiding-center '// &
                    'integrator'
            end select
        end if

        if (isw_field_type == TEST) then
            ! TEST field is fully analytic - no field file needed
            call init_field_can(isw_field_type, n_r=canonical_grid_nr, &
                n_th=canonical_grid_ntheta, n_phi=canonical_grid_nphi, &
                transformation_relerr=canonical_ode_relerr)
            call print_phase_time('Canonical field initialization completed')
        else if (isw_field_type == CANFLUX .or. isw_field_type == BOOZER .or. &
                isw_field_type == MEISS .or. isw_field_type == ALBERT) then
            call init_field_can(isw_field_type, field_temp, canonical_grid_nr, &
                canonical_grid_ntheta, canonical_grid_nphi, &
                canonical_ode_relerr)
            if (isw_field_type == BOOZER .or. &
                    isw_field_type == CANFLUX) then
                call set_magfie_refcoords_field(field_temp)
            end if
            call print_phase_time('Canonical field initialization completed')
        end if
    end subroutine init_field

    subroutine init_spectre_field(self)
        !> VMEC-free SPECTRE setup. Mirrors the Boozer-chartmap precedent: no
        !> init_vmec; the equilibrium globals nper and rmajor are taken from the
        !> SPECTRE file so stevvo and params_init produce a consistent
        !> dphi/dtaumin/fper. rmajor is the m=0,n=0 R harmonic of the outermost
        !> interface (SI meters); stevvo scales it to cm. For integmode > 0 the
        !> per-volume Meiss canonical coordinates are built here (#439).
        use field_spectre, only: spectre_field_t, create_spectre_field
        use field_can_spectre, only: set_spectre_construction_grid
        use magfie_sub, only: set_magfie_spectre_field, SPECTRE
        use new_vmec_stuff_mod, only: nper, rmajor
        use util, only: twopi
        use params, only: field_input, spectre_ncon_r, spectre_ncon_th, &
            spectre_ncon_phi, spectre_ncon_order, &
            spectre_ncon_ode_max_steps
        use params, only: spectre_ncon_ode_relerr
        use timing, only: print_phase_time
        use orbit_symplectic_base, only: sympl_rmax

        type(tracer_t), intent(inout) :: self

        type(spectre_field_t) :: sf
        integer :: ierr, ii
        real(dp) :: r00

        call create_spectre_field(sf, field_input, ierr)
        if (ierr /= 0) then
            print *, 'init_spectre_field: create_spectre_field failed for ', &
                trim(field_input), ' (ierr = ', ierr, ')'
            error stop
        end if
        call set_magfie_spectre_field(sf)

        r00 = 0.0d0
        do ii = 1, sf%data%mn
            if (sf%data%im(ii) == 0) then
                if (sf%data%in(ii) == 0) then
                    r00 = sf%data%Rbc(ii, sf%data%Mvol)
                end if
            end if
        end do
        if (r00 <= 0.0d0) error stop &
            'init_spectre_field: no positive m=0,n=0 R harmonic on outer interface'

        nper = sf%data%Nfp
        rmajor = r00
        self%fper = twopi/real(sf%data%Nfp, dp)

        print *, 'SPECTRE field: Nfp = ', sf%data%Nfp, ' Mvol = ', sf%data%Mvol, &
            ' rmajor = ', rmajor, ' m'

        if (self%integmode > 0) then
            ! The symplectic Newton/loss guards treat r > sympl_rmax as "left the
            ! domain". SPECTRE integrates in rho_g in [0, Mvol] and loses markers
            ! through the crossing pipeline at rho_g = Mvol, so the guard sits one
            ! volume width beyond: at exactly Mvol it would abort Newton iterates
            ! mid-solve during the exact-landing substep onto the outermost
            ! interface (#441). Beyond the edge the locked per-volume field is
            ! extended linearly, so iterates out there stay finite.
            sympl_rmax = real(sf%data%Mvol + 1, dp)
            call set_spectre_construction_grid(spectre_ncon_r, spectre_ncon_th, &
                spectre_ncon_phi, spectre_ncon_order, &
                spectre_ncon_ode_max_steps, &
                spectre_ncon_ode_relerr)
            call init_field_can(SPECTRE, sf)
            call print_phase_time('SPECTRE per-volume canonical construction completed')
        end if
    end subroutine init_spectre_field

    subroutine init_stl_wall_if_enabled(coord_file)
        character(len=*), intent(in) :: coord_file

        type(chartmap_metadata_t) :: meta
        character(len=16) :: units

        wall_enabled = .false.
        map_boundary_exit_code = ORBIT_EXIT_LCFS
        chart_boundary_kind_effective = 'lcfs'
        chartmap_cart_scale_to_m = -1.0d0
        wall_query_rho_lcfs = -1.0d0

        if (.not. allocated(ref_coords)) then
            if (len_trim(wall_input) == 0) return
            error stop "wall_input set but reference coordinates not initialized"
        end if

        select type (ref_coords)
        class is (chartmap_coordinate_system_t)
            call read_chartmap_metadata(coord_file, meta)
            chartmap_cart_scale_to_m = meta%cart_scale_to_m
            select case (trim(chart_boundary_kind))
            case ('auto')
                if (abs(meta%rho_lcfs - 1d0) <= 1d-12) then
                    map_boundary_exit_code = ORBIT_EXIT_LCFS
                    chart_boundary_kind_effective = 'lcfs'
                else
                    map_boundary_exit_code = ORBIT_EXIT_NUMERICAL_DOMAIN
                    chart_boundary_kind_effective = 'domain'
                end if
            case ('lcfs')
                map_boundary_exit_code = ORBIT_EXIT_LCFS
                chart_boundary_kind_effective = 'lcfs'
            case ('wall')
                map_boundary_exit_code = ORBIT_EXIT_WALL
                chart_boundary_kind_effective = 'wall'
            case ('domain')
                map_boundary_exit_code = ORBIT_EXIT_NUMERICAL_DOMAIN
                chart_boundary_kind_effective = 'domain'
            end select
        class default
            if (len_trim(wall_input) == 0) return
            error stop "wall_input requires chartmap reference coordinates"
        end select

        if (len_trim(wall_input) == 0) return

        ! The STL is an external physical wall enclosing the LCFS. Querying a
        ! straight Cartesian chord while the accepted endpoint is still inside
        ! the LCFS can create a false intersection in a non-convex toroidal
        ! geometry. Preserve the historical SIMPLE rule: activate wall queries
        ! only once an accepted microstep endpoint has passed the recorded LCFS.
        wall_query_rho_lcfs = meta%rho_lcfs

        units = adjustl(wall_units)
        select case (trim(units))
        case ("m", "M")
            call stl_wall_init(wall, trim(wall_input), 1.0d0)
        case ("mm", "MM")
            call stl_wall_init(wall, trim(wall_input), 1.0d-3)
        case default
            print *, "init_stl_wall_if_enabled: invalid wall_units=", trim(wall_units)
            error stop "wall_units must be m or mm"
        end select

        wall_enabled = .true.
    end subroutine init_stl_wall_if_enabled

    function select_vmec_equilibrium_file(vmec_file_default, field_file, coord_file) &
            result(vmec_file)
        use libneo_coordinates, only: detect_refcoords_file_type, &
            refcoords_file_vmec_wout
        character(*), intent(in) :: vmec_file_default
        character(*), intent(in) :: field_file
        character(*), intent(in) :: coord_file
        character(:), allocatable :: vmec_file

        integer :: file_type, ierr
        character(len=2048) :: message

        vmec_file = trim(vmec_file_default)

        if (len_trim(field_file) > 0) then
            if (ends_with_nc(field_file)) then
                call detect_refcoords_file_type(trim(field_file), file_type, &
                    ierr, message)
                if (ierr == 0 .and. file_type == refcoords_file_vmec_wout) then
                    vmec_file = trim(field_file)
                    return
                end if
            end if
        end if

        if (len_trim(coord_file) > 0) then
            if (ends_with_nc(coord_file)) then
                call detect_refcoords_file_type(trim(coord_file), file_type, &
                    ierr, message)
                if (ierr == 0 .and. file_type == refcoords_file_vmec_wout) then
                    vmec_file = trim(coord_file)
                    return
                end if
            end if
        end if
    end function select_vmec_equilibrium_file

    function ends_with_nc(filename) result(is_nc)
        character(*), intent(in) :: filename
        logical :: is_nc
        integer :: n

        n = len_trim(filename)
        is_nc = .false.
        if (n >= 3) then
            if (filename(n - 2:n) == '.nc') then
                is_nc = .true.
            end if
        end if
    end function ends_with_nc

    subroutine trace_parallel(norb)
        use netcdf_orbit_output, only: init_orbit_netcdf, close_orbit_netcdf, &
            write_orbit_to_netcdf
#ifdef SIMPLE_ENABLE_DEBUG_OUTPUT
        use params, only: debug
#endif

        type(tracer_t), intent(inout) :: norb
        integer :: i
        real(dp), allocatable :: traj(:, :), times(:)

        if (output_orbits_macrostep) then
            call init_orbit_netcdf(ntestpart, ntimstep)
        end if

        !$omp parallel firstprivate(norb) private(traj, times, i)
        allocate (traj(5, ntimstep), times(ntimstep))

        ! Independent markers have highly nonuniform trace costs. One-marker
        ! dynamic chunks keep unfinished work available to every thread instead
        ! of queueing it behind a pathological marker in a static block.
        !$omp do schedule(dynamic, 1)
        do i = 1, ntestpart
            if (allocated(particle_done)) then
                if (particle_done(i)) then
                    call progress_tick
                    cycle
                end if
            end if

#ifdef SIMPLE_ENABLE_DEBUG_OUTPUT
            if (debug) then
                !$omp critical
                kpart = kpart + 1
                print *, kpart, ' / ', ntestpart, 'particle: ', i, 'thread: ', &
                    omp_get_thread_num()
                !$omp end critical
            end if
#endif

            call trace_orbit(norb, i, traj, times)

            if (output_orbits_macrostep) then
                !$omp critical
                call write_orbit_to_netcdf(i, traj, times)
                !$omp end critical
            end if

            call progress_tick
        end do
        !$omp end do
        !$omp end parallel

        if (output_orbits_macrostep) then
            call close_orbit_netcdf()
        end if
    end subroutine trace_parallel

    subroutine trace_compare_gpu(norb)
        ! Validate and benchmark the OpenACC GPU tracing kernel against the CPU
        ! symplectic integrator on identical per-particle initial states.
        ! Triggered by the SIMPLE_GPU_BENCH environment variable and restricted
        ! to the Boozer + EXPL_IMPL_EULER path without wall, collision, or
        ! classifier options.
        use orbit_symplectic, only: orbit_timestep_sympl
        use orbit_symplectic_base, only: symplectic_integrator_t, EXPL_IMPL_EULER
        use field_can_mod, only: field_can_t
        use magfie_sub, only: BOOZER
        use boozer_sub, only: sync_boozer_state
        use simple_gpu, only: trace_orbits_gpu

        type(tracer_t), intent(inout) :: norb

        type(symplectic_integrator_t), allocatable :: si_cpu(:), si_gpu(:)
        type(field_can_t), allocatable :: f_cpu(:), f_gpu(:)
        integer, allocatable :: cpu_loss(:), gpu_loss(:)
        real(dp), allocatable :: cpu_zend(:, :), gpu_zend(:, :)
        real(dp) :: z(5)
        integer :: i, it, ktau, ierr, loss_mismatch
        integer :: cpu_lost, gpu_lost, flip
        real(dp) :: t0, t1, t_cpu, t_gpu, maxz

        if (isw_field_type /= BOOZER .or. integmode /= EXPL_IMPL_EULER .or. swcoll .or. &
            len_trim(wall_input) > 0 .or. class_plot .or. fast_class .or. generate_start_only) then
            error stop "simple_main.trace_compare_gpu: SIMPLE_GPU_BENCH requires Boozer, " // &
                "EXPL_IMPL_EULER, and no wall/collision/classifier options"
        end if

        call sync_boozer_state

        allocate (si_cpu(ntestpart), si_gpu(ntestpart))
        allocate (f_cpu(ntestpart), f_gpu(ntestpart))
        allocate (cpu_loss(ntestpart), gpu_loss(ntestpart))
        allocate (cpu_zend(4, ntestpart), gpu_zend(4, ntestpart))

        ! Identical per-particle initialisation (host). init_sympl sets the
        ! orbit_timestep_sympl procedure pointer for the CPU reference.
        do i = 1, ntestpart
            call ref_to_integ(zstart(1:3, i), z(1:3))
            z(4:5) = zstart(4:5, i)
            call init_sympl(si_cpu(i), f_cpu(i), z, dtaumin, dtaumin, relerr, integmode)
            si_gpu(i) = si_cpu(i)
            f_gpu(i) = f_cpu(i)
        end do

        ! CPU reference (OpenMP over particles)
        t0 = omp_get_wtime()
        !$omp parallel do private(i, it, ktau, ierr) schedule(dynamic)
        do i = 1, ntestpart
            ierr = 0
            cpu_loss(i) = ntimstep
            do it = 2, ntimstep
                do ktau = 1, ntau_macro(it)
                    call orbit_timestep_sympl(si_cpu(i), f_cpu(i), ierr)
                    if (ierr /= 0) exit
                end do
                if (ierr /= 0) then
                    cpu_loss(i) = it
                    exit
                end if
            end do
            cpu_zend(:, i) = si_cpu(i)%z(1:4)
        end do
        !$omp end parallel do
        t1 = omp_get_wtime()
        t_cpu = t1 - t0

        if (any(cpu_loss < ntimstep)) then
            error stop 'SIMPLE_GPU_BENCH does not support boundary-reaching traces'
        end if

        ! GPU kernel
        t0 = omp_get_wtime()
        call trace_orbits_gpu(si_gpu, f_gpu, ntestpart, ntimstep, ntau_macro, &
            gpu_loss, gpu_zend)
        t1 = omp_get_wtime()
        t_gpu = t1 - t0

        ! Compare
        maxz = 0d0
        loss_mismatch = 0
        cpu_lost = 0
        gpu_lost = 0
        flip = 0
        do i = 1, ntestpart
            maxz = max(maxz, maxval(dabs(cpu_zend(:, i) - gpu_zend(:, i))))
            if (cpu_loss(i) /= gpu_loss(i)) loss_mismatch = loss_mismatch + 1
            if (cpu_loss(i) < ntimstep) cpu_lost = cpu_lost + 1
            if (gpu_loss(i) < ntimstep) gpu_lost = gpu_lost + 1
            if ((cpu_loss(i) < ntimstep) .neqv. (gpu_loss(i) < ntimstep)) flip = flip + 1
        end do

        print *, '==================== GPU vs CPU tracing ===================='
        print '(a,i0,a,i0)', ' particles = ', ntestpart, '   timesteps = ', ntimstep
        print '(a,es12.4)', ' max |z_cpu - z_gpu| (final state) = ', maxz
        print '(a,i0,a,i0)', ' loss-step mismatches = ', loss_mismatch, ' / ', ntestpart
        print '(a,i0,a,i0,a,f7.4)', ' CPU lost = ', cpu_lost, ' / ', ntestpart, &
            '   confined frac = ', 1d0 - real(cpu_lost, dp)/real(ntestpart, dp)
        print '(a,i0,a,i0,a,f7.4)', ' GPU lost = ', gpu_lost, ' / ', ntestpart, &
            '   confined frac = ', 1d0 - real(gpu_lost, dp)/real(ntestpart, dp)
        print '(a,i0,a,i0)', ' lost<->confined flips = ', flip, ' / ', ntestpart
        print '(a,f10.4,a)', ' CPU time (OpenMP) = ', t_cpu, ' s'
        print '(a,f10.4,a)', ' GPU time          = ', t_gpu, ' s'
        if (t_gpu > 0d0) print '(a,f8.2,a)', ' speedup (CPU/GPU) = ', t_cpu/t_gpu, ' x'
        print *, '============================================================'
    end subroutine trace_compare_gpu

    subroutine classify_parallel(norb)
        use classification, only: trace_orbit_with_classifiers, classification_result_t
        use params, only: class_passing, class_lost
#ifdef SIMPLE_ENABLE_DEBUG_OUTPUT
        use params, only: debug
#endif

        type(tracer_t), intent(inout) :: norb
        integer :: i
        type(classification_result_t) :: class_result

        !$omp parallel firstprivate(norb) private(class_result, i)
        !$omp do
        do i = 1, ntestpart
#ifdef SIMPLE_ENABLE_DEBUG_OUTPUT
            if (debug) then
                !$omp critical
                kpart = kpart + 1
                print *, kpart, ' / ', ntestpart, 'particle: ', i, 'thread: ', &
                    omp_get_thread_num()
                !$omp end critical
            end if
#endif

            if (swcoll) call reset_seed_if_deterministic
            call trace_orbit_with_classifiers(norb, i, class_result)

            ! Store classification flags in params arrays
            class_passing(i) = class_result%passing
            class_lost(i) = class_result%lost
            ! iclass already populated by trace_orbit_with_classifiers
            ! Other results (zend, times_lost, trap_par, perp_inv) also already stored
        end do
        !$omp end do
        !$omp end parallel
    end subroutine classify_parallel

    subroutine print_parameters
        print *, 'tau: ', dtau, dtaumin, min(dabs(mod(dtau, dtaumin)), &
            dabs(mod(dtau, dtaumin) - &
            dtaumin))/dtaumin, ntau
        print *, 'v0 = ', v0
    end subroutine print_parameters

    subroutine read_profiles_config(config_file)
        use simple_profiles, only: read_profiles

        character(256), intent(in) :: config_file

        call read_profiles(config_file)
    end subroutine read_profiles_config

    subroutine init_collisions
        use params, only: am1, am2, Z1, Z2, facE_al, dchichi, slowrate, &
            dchichi_norm, slowrate_norm, v0
        use simple_profiles, only: Te_scale, Ti1_scale, Ti2_scale, &
            ni1_scale, ni2_scale

        real(dp) :: v0_coll, ealpha
        real(dp) :: densi1, densi2, tempi1, tempi2, tempe

        ealpha = 3.5d6/facE_al
        densi1 = ni1_scale*1.0d-6
        densi2 = ni2_scale*1.0d-6
        tempi1 = Ti1_scale
        tempi2 = Ti2_scale
        tempe = Te_scale

        call loacol_alpha(am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, &
            ealpha, v0_coll, dchichi, slowrate, dchichi_norm, &
            slowrate_norm)

        if (abs(v0_coll - v0) > 1d-6) then
            error stop 'simple_main.init_collisions: v0_coll != v0'
        end if

        call init_collision_profiles(am1, am2, Z1, Z2, ealpha, v0)
    end subroutine init_collisions

    subroutine sample_particles(xstart_is_integ_coords)
        use samplers, only: sample, START_FILE, sample_grid, &
            sample_surface_fieldline, sample_surface_fieldline_from_integ

        logical, intent(in), optional :: xstart_is_integ_coords
        logical :: convert_surface_starts

        convert_surface_starts = .false.
        if (present(xstart_is_integ_coords)) then
            convert_surface_starts = xstart_is_integ_coords
        end if

        if (1 == startmode) then
            if ((0d0 < grid_density) .and. (1d0 > grid_density)) then
                call sample_grid(zstart, grid_density, convert_surface_starts)
            else
                if (convert_surface_starts) then
                    call sample_surface_fieldline_from_integ(zstart)
                else
                    call sample_surface_fieldline(zstart)
                end if
            end if

        elseif (2 == startmode) then
            call sample(zstart, START_FILE)

        elseif (3 == startmode) then
            call sample(special_ants_file)

        elseif (4 == startmode) then
            call sample(zstart, reuse_batch)

        elseif (5 == startmode) then
            if (0 == num_surf) then
                call sample(zstart, 0.0d0, 1.0d0)
            elseif (1 == num_surf) then
                call sample(zstart, 0.0d0, sbeg(1))
            elseif (2 == num_surf) then
                call sample(zstart, sbeg(1), sbeg(num_surf))
            else
                print *, 'Invalid surface range for volume sample defined ', &
                    '(2 < num_surf), stopping.'
                stop
            end if

        else
            print *, 'Unknown startmode: ', startmode
        end if
    end subroutine sample_particles

    subroutine sample_particles_test_field
        !> Sample particles for the analytic circular tokamak TEST field.
        !> TEST field uses (r, theta, phi) coordinates with B0=1, R0=1, a=0.5, iota=1.
        !> bmod = B0 * (1 - r/R0 * cos(theta))
        use params, only: ntestpart, sbeg, bmod00, bmin, bmax, zstart, &
            reset_seed_if_deterministic

        real(dp), parameter :: B0 = 1.0d0, R0 = 1.0d0, a = 0.5d0
        real(dp) :: r_start, tmp_rand
        integer :: ipart

        ! Use sbeg(1) as the starting minor radius (mapped to r for tokamak)
        ! sbeg(1) is input as a flux-like value; for TEST field, interpret as r/a
        r_start = sbeg(1)*a

        ! Set magnetic field bounds for this r value
        ! bmod = B0*(1 - r/R0*cos(theta))
        ! Maximum at theta=pi: bmax = B0*(1 + r/R0)
        ! Minimum at theta=0:  bmin = B0*(1 - r/R0)
        bmod00 = B0
        bmax = B0*(1.0d0 + r_start/R0)
        bmin = B0*(1.0d0 - r_start/R0)

        print *, 'TEST field: bmod00 = ', bmod00, 'bmin = ', bmin, 'bmax = ', bmax

        call reset_seed_if_deterministic

        ! Sample particles uniformly in (theta, phi) at fixed r
        do ipart = 1, ntestpart
            zstart(1, ipart) = r_start
            call random_number(tmp_rand)
            zstart(2, ipart) = twopi*tmp_rand
            call random_number(tmp_rand)
            zstart(3, ipart) = twopi*tmp_rand
            zstart(4, ipart) = 1.0d0
            call random_number(tmp_rand)
            zstart(5, ipart) = 2.0d0*(tmp_rand - 0.5d0)
        end do

        call save_starting_points_test(zstart)

    contains
        subroutine save_starting_points_test(zs)
            real(dp), intent(in) :: zs(:, :)
            integer :: i, unit

            open (newunit=unit, file='start.dat', recl=1024)
            do i = 1, size(zs, 2)
                write (unit, *) zs(:, i)
            end do
            close (unit)
        end subroutine save_starting_points_test
    end subroutine sample_particles_test_field

    subroutine init_bminmax
        use find_bminmax_sub, only: init_bminmax_arrays

        ! Populate bminmax arrays with the active tracing magfie backend.
        ! init_starting_surf supplies bmin/bmax for a single sampled surface;
        ! this cache is only needed when particles span multiple surfaces.
        call init_bminmax_arrays
    end subroutine init_bminmax

    logical function needs_bminmax_cache()
        ! Match the current readers. Classifier semantics are kept unchanged
        ! here; broad classifier changes belong in the classifier PR.
        if ((ntcut > 0) .or. class_plot) then
            needs_bminmax_cache = num_surf > 1
        else
            needs_bminmax_cache = num_surf /= 1
        end if
    end function needs_bminmax_cache

    subroutine init_counters
        icounter = 0 ! evaluation counter
        kpart = 0

        ! initialize array of confined particle percentage
        confpart_trap = 0.d0
        confpart_pass = 0.d0
        unresolved_orbits = 0

        ! initialize lost times when particles get lost
        times_lost = -1.d0
        orbit_exit_code = 0
        boundary_event_radial_residual = -1d0
        boundary_event_time_width = -1d0
    end subroutine init_counters

    subroutine compute_canonical_frequencies_flat(anorb, initial_state, &
            n_periods, max_steps, &
            values, metadata)
        use orbit_frequencies, only: frequency_options_t, frequency_result_t, &
            compute_canonical_frequencies

        type(tracer_t), intent(inout) :: anorb
        real(dp), intent(in) :: initial_state(5)
        integer, intent(in) :: n_periods, max_steps
        real(dp), intent(out) :: values(6)
        integer, intent(out) :: metadata(5)

        type(frequency_options_t) :: options
        type(frequency_result_t) :: result
        real(dp) :: z(5)

        call ref_to_integ(initial_state(1:3), z(1:3))
        z(4:5) = initial_state(4:5)

        anorb%dtaumin = dtaumin
        anorb%dtau = dtaumin
        anorb%v0 = v0
        anorb%relerr = relerr
        anorb%integmode = integmode
        options%n_periods = n_periods
        options%max_steps = max_steps

        call compute_canonical_frequencies(anorb, z, options, result)

        values(1) = result%period
        values(2) = result%period_std
        values(3) = result%delta_phi
        values(4) = result%delta_phi_std
        values(5) = result%omega_b
        values(6) = result%omega_phi
        metadata(1) = result%status
        metadata(2) = result%orbit_class
        metadata(3) = result%parallel_direction
        metadata(4) = result%n_periods
        metadata(5) = result%n_steps
    end subroutine compute_canonical_frequencies_flat

    subroutine invariants_from_state_flat(anorb, initial_state, values, status)
        use orbit_invariants, only: guiding_center_invariants_t, &
            invariants_from_state, invariant_flux_convention

        type(tracer_t), intent(in) :: anorb
        real(dp), intent(in) :: initial_state(5)
        real(dp), intent(out) :: values(5)
        integer, intent(out) :: status

        type(guiding_center_invariants_t) :: invariants
        real(dp) :: state(5)

        call ref_to_integ(initial_state(1:3), state(1:3))
        state(4:5) = initial_state(4:5)
        call invariants_from_state(anorb, state, invariants, status)
        values(1:3) = [invariants%h0, invariants%j_perp, invariants%p_phi]
        call invariant_flux_convention(anorb, values(4), values(5))
    end subroutine invariants_from_state_flat

    subroutine potato_invariants_to_simple_flat(potato_values, psi_axis, &
            psi_edge, v0_ratio, simple_values)
        use orbit_invariants, only: guiding_center_invariants_t, &
            potato_to_simple_invariants

        real(dp), intent(in) :: potato_values(3)
        real(dp), intent(in) :: psi_axis, psi_edge, v0_ratio
        real(dp), intent(out) :: simple_values(3)

        type(guiding_center_invariants_t) :: potato, converted

        potato%h0 = potato_values(1)
        potato%j_perp = potato_values(2)
        potato%p_phi = potato_values(3)
        call potato_to_simple_invariants(potato, psi_axis, psi_edge, v0_ratio, &
            converted)
        simple_values = [converted%h0, converted%j_perp, converted%p_phi]
    end subroutine potato_invariants_to_simple_flat

    subroutine states_from_invariants_flat(anorb, values, max_solutions, &
            states, metadata, residuals, cylindrical, n_solutions, status)
        use orbit_invariants, only: guiding_center_invariants_t, &
            invariant_start_result_t, states_from_invariants, &
            INVARIANT_CAPACITY

        type(tracer_t), intent(in) :: anorb
        real(dp), intent(in) :: values(3)
        integer, intent(in) :: max_solutions
        real(dp), intent(out) :: states(5, max_solutions)
        integer, intent(out) :: metadata(3, max_solutions)
        real(dp), intent(out) :: residuals(max_solutions)
        real(dp), intent(out) :: cylindrical(3, max_solutions)
        integer, intent(out) :: n_solutions, status

        type(guiding_center_invariants_t) :: invariants
        type(invariant_start_result_t) :: result
        integer :: i

        states = 0.0_dp
        metadata = 0
        residuals = 0.0_dp
        cylindrical = 0.0_dp
        n_solutions = 0
        invariants%h0 = values(1)
        invariants%j_perp = values(2)
        invariants%p_phi = values(3)
        call states_from_invariants(anorb, invariants, result)
        status = result%status
        n_solutions = min(size(result%starts), max_solutions)
        do i = 1, n_solutions
            states(:, i) = result%starts(i)%state
            metadata(:, i) = [result%starts(i)%sigma, &
                result%starts(i)%section_branch, result%starts(i)%section_kind]
            residuals(i) = result%starts(i)%residual
            if (allocated(ref_coords)) then
                call ref_coords%evaluate_cyl(states(1:3, i), cylindrical(:, i))
            end if
        end do
        call sort_invariant_starts(states, metadata, residuals, cylindrical, &
            n_solutions)
        if (size(result%starts) > max_solutions) status = INVARIANT_CAPACITY
    end subroutine states_from_invariants_flat

    subroutine sort_invariant_starts(states, metadata, residuals, cylindrical, &
            count)
        real(dp), intent(inout) :: states(:, :), residuals(:), cylindrical(:, :)
        integer, intent(inout) :: metadata(:, :)
        integer, intent(in) :: count

        real(dp) :: state_tmp(size(states, 1)), residual_tmp
        real(dp) :: cylindrical_tmp(size(cylindrical, 1))
        integer :: metadata_tmp(size(metadata, 1))
        integer :: i, j

        if (.not. allocated(ref_coords)) return
        do i = 2, count
            state_tmp = states(:, i)
            metadata_tmp = metadata(:, i)
            residual_tmp = residuals(i)
            cylindrical_tmp = cylindrical(:, i)
            j = i - 1
            do while (j >= 1)
                if (cylindrical(1, j) <= cylindrical_tmp(1)) exit
                states(:, j + 1) = states(:, j)
                metadata(:, j + 1) = metadata(:, j)
                residuals(j + 1) = residuals(j)
                cylindrical(:, j + 1) = cylindrical(:, j)
                j = j - 1
            end do
            states(:, j + 1) = state_tmp
            metadata(:, j + 1) = metadata_tmp
            residuals(j + 1) = residual_tmp
            cylindrical(:, j + 1) = cylindrical_tmp
        end do
    end subroutine sort_invariant_starts

    subroutine trace_to_cut_flat(anorb, initial_state, requested_cut_type, &
            max_events, cut_state, cut_type, status)
        use cut_detector, only: cut_detector_t, init_cut_detector => init, &
            trace_to_cut

        type(tracer_t), intent(inout) :: anorb
        real(dp), intent(in) :: initial_state(5)
        integer, intent(in) :: requested_cut_type, max_events
        real(dp), intent(out) :: cut_state(6)
        integer, intent(out) :: cut_type, status

        type(cut_detector_t) :: detector
        real(dp) :: z(5), reference_position(3)
        integer :: event_index, ierr

        cut_state = 0.0_dp
        cut_type = -1
        status = 1
        if (requested_cut_type < -1 .or. requested_cut_type > 1) return
        if (max_events < 1 .or. integmode <= 0) return

        call ref_to_integ(initial_state(1:3), z(1:3))
        z(4:5) = initial_state(4:5)
        anorb%dtaumin = dtaumin
        anorb%dtau = dtaumin
        anorb%v0 = v0
        anorb%relerr = relerr
        anorb%integmode = integmode

        call init_sympl(anorb%si, anorb%f, z, dtaumin, dtaumin, relerr, integmode)
        call init_cut_detector(detector, anorb%fper, z)

        do event_index = 1, max_events
            call trace_to_cut(detector, anorb%si, anorb%f, z, cut_state, &
                cut_type, ierr)
            if (ierr /= 0) then
                status = ierr
                return
            end if
            if (requested_cut_type == -1) exit
            if (cut_type == requested_cut_type) exit
        end do

        if (requested_cut_type /= -1) then
            if (cut_type /= requested_cut_type) then
                status = 4
                return
            end if
        end if

        call integ_to_ref(cut_state(1:3), reference_position)
        cut_state(1:3) = reference_position
        status = 0
    end subroutine trace_to_cut_flat
    pure integer function classify_orbit_exit(ierr_orbit, model, mode, &
            boundary_exit_code)
        integer, intent(in) :: ierr_orbit, model, mode
        integer, intent(in) :: boundary_exit_code

        if (ierr_orbit == 0) then
            classify_orbit_exit = ORBIT_EXIT_COMPLETED
        else if (ierr_orbit == 77) then
            classify_orbit_exit = ORBIT_EXIT_WALL
        else if (model == ORBIT_FULL_ORBIT) then
            if (ierr_orbit == ORBIT_FO_LOSS) then
                classify_orbit_exit = boundary_exit_code
            else
                classify_orbit_exit = ORBIT_EXIT_NUMERICAL_FULL_ORBIT
            end if
        else if (mode <= 0) then
            if (ierr_orbit == 1) then
                classify_orbit_exit = boundary_exit_code
            else
                classify_orbit_exit = ORBIT_EXIT_NUMERICAL_EVENT
            end if
        else
            select case (ierr_orbit)
            case (SYMPLECTIC_STEP_BOUNDARY)
                classify_orbit_exit = boundary_exit_code
            case (SYMPLECTIC_STEP_OUTSIDE_DOMAIN)
                classify_orbit_exit = ORBIT_EXIT_NUMERICAL_DOMAIN
            case (SYMPLECTIC_STEP_MAXITER)
                classify_orbit_exit = ORBIT_EXIT_NUMERICAL_MAXITER
            case (SYMPLECTIC_STEP_LINEAR_SOLVE)
                classify_orbit_exit = ORBIT_EXIT_NUMERICAL_LINEAR
            case (SYMPLECTIC_STEP_EVENT_NOT_CONVERGED)
                classify_orbit_exit = ORBIT_EXIT_NUMERICAL_EVENT
            case (SYMPLECTIC_STEP_BOUNDARY_LIMITED)
                classify_orbit_exit = ORBIT_EXIT_NUMERICAL_EVENT
            case default
                classify_orbit_exit = ORBIT_EXIT_NUMERICAL_EVENT
            end select
        end if
    end function classify_orbit_exit

    subroutine trace_orbit(anorb, ipart, orbit_traj, orbit_times)
        use classification, only: trace_orbit_with_classifiers, &
            classification_result_t, &
            write_classification_results
        use magfie_sub, only: SPECTRE
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

        type(tracer_t), intent(inout) :: anorb
        integer, intent(in) :: ipart
        real(dp), intent(out) :: orbit_traj(:, :) ! (5, ntimstep)
        real(dp), intent(out) :: orbit_times(:) ! (ntimstep)

        real(dp), dimension(5) :: z
        real(dp) :: u_ref_prev(3), x_prev(3), x_prev_m(3), exit_step
        integer :: it, ierr_orbit, it_final, it_f
        integer :: first_unresolved_it, hold_streak
        integer(8) :: kt
        logical :: passing, faulted, numerical_hold, physical_exit
        type(rk_recovery_state_t) :: rk_recovery
        type(classification_result_t) :: class_result

        ierr_orbit = 0
        faulted = .false.
        physical_exit = .false.
        first_unresolved_it = 0
        hold_streak = 0
        rk_recovery = rk_recovery_state_t()
        exit_step = -1d0
        orbit_traj = ieee_value(0.0d0, ieee_quiet_nan)
        orbit_times = ieee_value(0.0d0, ieee_quiet_nan)

        if (swcoll) call reset_seed_if_deterministic

        if (ntcut > 0 .or. class_plot) then
            call trace_orbit_with_classifiers(anorb, ipart, class_result)
            if (class_plot) then
                call write_classification_results(ipart, class_result)
            end if
            return
        end if

        call ref_to_integ(zstart(1:3, ipart), z(1:3))
        z(4:5) = zstart(4:5, ipart)
        orbit_traj(:, 1) = zstart(:, ipart)
        orbit_times(1) = 0.0_dp
        zend(:, ipart) = 0d0

        if (wall_enabled) then
            u_ref_prev = zstart(1:3, ipart)
            call ref_coords%evaluate_cart(u_ref_prev, x_prev)
            x_prev_m = x_prev*chartmap_cart_scale_to_m
        end if

        if (orbit_model == ORBIT_FULL_ORBIT) then
            if (wall_enabled) error stop &
                'orbit_model=ORBIT_FULL_ORBIT does not support wall loss yet'
            if (swcoll) error stop &
                'orbit_model=ORBIT_FULL_ORBIT does not support collisions yet'
            call init_fo(anorb%fo, z, dtaumin)
        else if (integmode > 0) then
            call init_sympl(anorb%si, anorb%f, z, dtaumin, dtaumin, relerr, integmode)
        end if

        call compute_pitch_angle_params(z, passing, trap_par(ipart), perp_inv(ipart))

        ! SPECTRE traces every marker per volume: the passing-skip optimisation
        ! (a well-confined-tokamak assumption) would wrongly mark interface-bound
        ! markers as confined, so it is bypassed here.
        if (isw_field_type == SPECTRE) then
            if (integmode > 0) then
                call trace_orbit_spectre_sympl(anorb, ipart, z, passing, &
                    orbit_traj, orbit_times)
            else
                call trace_orbit_spectre(ipart, z, passing, orbit_traj, orbit_times)
            end if
            return
        end if

        if (passing .and. should_skip(ipart)) then
            zend(:, ipart) = zstart(:, ipart)
            times_lost(ipart) = -1.d0
            orbit_exit_code(ipart) = ORBIT_EXIT_SKIPPED
            do it = 1, ntimstep
                !$omp atomic update
                confpart_pass(it) = confpart_pass(it) + 1.d0
            end do
            return
        end if

        kt = 0
        it_final = 0
        do it = 1, ntimstep
            if (it >= 2) then
                if (wall_enabled) then
                    ! Use microstep-level wall checking for STL walls
                    call macrostep_with_wall_check(anorb, z, kt, ierr_orbit, &
                        ntau_macro(it), ipart, x_prev_m, exit_step, &
                        hold_streak=hold_streak, &
                        numerical_hold_any=numerical_hold, &
                        rk_recovery=rk_recovery)
                else
                    call macrostep(anorb, z, kt, ierr_orbit, ntau_macro(it), &
                        exit_step, hold_streak=hold_streak, &
                        numerical_hold_any=numerical_hold, &
                        rk_recovery=rk_recovery)
                end if
            end if

            if (ierr_orbit .ne. 0) then
                orbit_exit_code(ipart) = classify_orbit_exit( &
                    ierr_orbit, orbit_model, integmode, map_boundary_exit_code)
                it_final = it
                physical_exit = orbit_exit_code(ipart) == ORBIT_EXIT_LCFS .or. &
                    orbit_exit_code(ipart) == ORBIT_EXIT_WALL
                if (physical_exit) then
                    orbit_traj(:, it) = z
                    orbit_times(it) = exit_step*dtaumin/v0
                    if (orbit_exit_code(ipart) == ORBIT_EXIT_LCFS .and. &
                        orbit_model == ORBIT_GC .and. integmode > 0 .and. &
                        ierr_orbit == SYMPLECTIC_STEP_BOUNDARY) then
                        boundary_event_radial_residual(ipart) = &
                            anorb%si%last_event_radial_residual
                        boundary_event_time_width(ipart) = &
                            anorb%si%last_event_fraction_width*dtaumin/v0
                    end if
                else
                    if (first_unresolved_it == 0) first_unresolved_it = it
                    faulted = .true.
                end if
                exit
            end if

            ! Store trajectory data (after macrostep so time is correct)
            orbit_traj(:, it) = z
            orbit_times(it) = kt*dtaumin/v0

            if (first_unresolved_it == 0) call increase_confined_count(it, passing)
            it_final = it
        end do

        if (first_unresolved_it > 0) then
            do it_f = first_unresolved_it, ntimstep
                !$omp atomic update
                unresolved_orbits(it_f) = unresolved_orbits(it_f) + 1
            end do
            physical_exit = .false.
            faulted = .true.
            if (orbit_exit_code(ipart) < ORBIT_EXIT_NUMERICAL_DOMAIN) &
                orbit_exit_code(ipart) = ORBIT_EXIT_NUMERICAL_EVENT
        end if

        ! Fill remaining timesteps with NaN if particle left domain early
        if (it_final < ntimstep) then
            do it = it_final + 1, ntimstep
                orbit_traj(:, it) = ieee_value(0.0d0, ieee_quiet_nan)
                orbit_times(it) = ieee_value(0.0d0, ieee_quiet_nan)
            end do
        end if

        !$omp critical
        call integ_to_ref(z(1:3), zend(1:3, ipart))
        zend(4:5, ipart) = z(4:5)
        if (physical_exit) then
            times_lost(ipart) = exit_step*dtaumin/v0
        else if (faulted) then
            times_lost(ipart) = ieee_value(0.0d0, ieee_quiet_nan)
        else
            times_lost(ipart) = kt*dtaumin/v0
        end if
        !$omp end critical
    end subroutine trace_orbit

    subroutine trace_orbit_spectre(ipart, z, passing, orbit_traj, orbit_times)
        !> Per-volume RK45 guiding-center trace for SPECTRE (integmode=0). The
        !> marker traverses volumes through the Level-0 crossing map and stays
        !> confined for the full trace, reflects at forbidden interfaces, or is
        !> lost at the outermost interface. Every crossing/reflection is logged.
        use spectre_orbit, only: spectre_orbit_state_t, spectre_event_t, &
            spectre_state_reset, orbit_timestep_spectre, &
            SPECTRE_OK, SPECTRE_BOUNDARY, SPECTRE_FAULT
        use magfie_sub, only: spectre_field
        use interface_crossing, only: crossing_log_record
        use params, only: crossing_level
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

        integer, intent(in) :: ipart
        real(dp), intent(inout) :: z(5)
        logical, intent(in) :: passing
        real(dp), intent(out) :: orbit_traj(:, :)
        real(dp), intent(out) :: orbit_times(:)

        type(spectre_orbit_state_t) :: state
        type(spectre_event_t) :: event
        integer :: it, it_f, ktau, ierr_orbit, it_final
        integer(8) :: kt
        real(dp) :: t_stop

        call spectre_state_reset(state, spectre_field%data%Mvol)

        kt = 0
        it_final = 0
        ierr_orbit = SPECTRE_OK
        do it = 1, ntimstep
            if (it >= 2) then
                do ktau = 1, ntau_macro(it)
                    call orbit_timestep_spectre(state, z, dtaumin, relerr, &
                        crossing_level, ierr_orbit, event)
                    if (event%occurred) then
                        call crossing_log_record(ipart, &
                            (real(kt, dp) + event%t_frac)*dtaumin/v0, event%info)
                    end if
                    if (ierr_orbit == SPECTRE_FAULT .and. &
                        symplectic_newton_warning_mode) then
                        ! orbit_timestep_spectre rolls back to the last accepted
                        ! position. Consume this held interval and retry the same
                        ! marker on the next microstep; warning mode never turns an
                        ! isolated RK failure into a terminal marker.
                        call count_event(EVT_WARNING_STEP_SKIP)
                        ierr_orbit = SPECTRE_OK
                        kt = kt + 1
                        cycle
                    end if
                    if (ierr_orbit /= SPECTRE_OK) exit
                    kt = kt + 1
                end do
            end if

            if (ierr_orbit /= SPECTRE_OK) then
                it_final = it
                if (ierr_orbit == SPECTRE_BOUNDARY) then
                    orbit_exit_code(ipart) = ORBIT_EXIT_LCFS
                else
                    orbit_exit_code(ipart) = ORBIT_EXIT_NUMERICAL_EVENT
                    do it_f = it, ntimstep
                        !$omp atomic update
                        unresolved_orbits(it_f) = unresolved_orbits(it_f) + 1
                    end do
                end if
                exit
            end if

            orbit_traj(:, it) = z
            orbit_times(it) = kt*dtaumin/v0
            call increase_confined_count(it, passing)
            it_final = it
        end do

        if (it_final < ntimstep) then
            do it = it_final + 1, ntimstep
                orbit_traj(:, it) = ieee_value(0.0d0, ieee_quiet_nan)
                orbit_times(it) = ieee_value(0.0d0, ieee_quiet_nan)
            end do
        end if

        if (ierr_orbit == SPECTRE_BOUNDARY) then
            t_stop = (real(kt, dp) + event%t_frac)*dtaumin/v0
        else if (ierr_orbit /= SPECTRE_OK) then
            t_stop = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            t_stop = real(kt, dp)*dtaumin/v0
            orbit_exit_code(ipart) = ORBIT_EXIT_COMPLETED
        end if

        !$omp critical
        call integ_to_ref(z(1:3), zend(1:3, ipart))
        zend(4:5, ipart) = z(4:5)
        times_lost(ipart) = t_stop
        !$omp end critical
    end subroutine trace_orbit_spectre

    subroutine trace_orbit_spectre_sympl(anorb, ipart, z, passing, orbit_traj, &
            orbit_times)
        !> Per-volume symplectic guiding-center trace for SPECTRE (integmode > 0).
        !> Each microstep resolves interface events by an exact-landing substep of
        !> the same implicit scheme, applies the crossing map, and re-canonicalizes
        !> in the resolved volume (#441). Markers are lost only at the outermost
        !> interface, matching the RK45 path; CROSS_STOP remains only as the
        !> pathological non-convergence fallback inside the microstepper.
        use spectre_sympl_orbit, only: sympl_spectre_state_t, sympl_spectre_reset, &
            orbit_microstep_sympl_spectre, &
            SYMPL_SPECTRE_OK, SYMPL_SPECTRE_LOSS, &
            SYMPL_SPECTRE_STOP, SYMPL_SPECTRE_SKIM
        use field_can_spectre, only: spectre_mvol, set_spectre_volume_lock
        use params, only: crossing_level
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

        type(tracer_t), intent(inout) :: anorb
        integer, intent(in) :: ipart
        real(dp), intent(inout) :: z(5)
        logical, intent(in) :: passing
        real(dp), intent(out) :: orbit_traj(:, :)
        real(dp), intent(out) :: orbit_times(:)

        type(sympl_spectre_state_t) :: state
        integer :: it, it_f, ktau, ierr_orbit, it_final
        integer :: first_unresolved_it
        integer(8) :: kt
        real(dp) :: t_stop, t_frac

        call sympl_spectre_reset(state, anorb%si, spectre_mvol, integmode, &
            crossing_level)

        kt = 0
        it_final = 0
        ierr_orbit = SYMPL_SPECTRE_OK
        t_frac = 1.0d0
        first_unresolved_it = 0

        do it = 1, ntimstep
            if (it >= 2) then
                do ktau = 1, ntau_macro(it)
                    call orbit_microstep_sympl_spectre(state, anorb%si, anorb%f, &
                        ipart, &
                        real(kt, dp)*dtaumin/v0, &
                        dtaumin/v0, ierr_orbit, &
                        t_frac)
                    if (ierr_orbit /= SYMPL_SPECTRE_OK) exit
                    kt = kt + 1
                end do
            end if

            if (ierr_orbit /= SYMPL_SPECTRE_OK) then
                it_final = it
                select case (ierr_orbit)
                case (SYMPL_SPECTRE_LOSS)
                    orbit_exit_code(ipart) = ORBIT_EXIT_LCFS
                case (SYMPL_SPECTRE_STOP)
                    orbit_exit_code(ipart) = ORBIT_EXIT_NUMERICAL_FULL_ORBIT
                case (SYMPL_SPECTRE_SKIM)
                    orbit_exit_code(ipart) = ORBIT_EXIT_COMPLETED
                case default
                    orbit_exit_code(ipart) = ORBIT_EXIT_NUMERICAL_EVENT
                end select
                if (ierr_orbit /= SYMPL_SPECTRE_LOSS .and. &
                    ierr_orbit /= SYMPL_SPECTRE_SKIM .and. &
                    first_unresolved_it == 0) first_unresolved_it = it
                exit
            end if

            if (state%fo%active .and. state%fo%has_y) then
                z = state%fo%last_y
            else
                call to_standard_z_coordinates(anorb, z)
            end if
            orbit_traj(:, it) = z
            orbit_times(it) = kt*dtaumin/v0
            if (first_unresolved_it == 0) call increase_confined_count(it, passing)
            it_final = it
        end do

        if (first_unresolved_it > 0) then
            do it_f = first_unresolved_it, ntimstep
                !$omp atomic update
                unresolved_orbits(it_f) = unresolved_orbits(it_f) + 1
            end do
        end if

        ! A mirror-confined (skimming) marker is confined for the rest of the
        ! trace, so it must count as confined at every remaining step for the
        ! confined_fraction series to stay consistent with times_lost.
        if (ierr_orbit == SYMPL_SPECTRE_SKIM) then
            do it = it_final, ntimstep
                call increase_confined_count(it, passing)
            end do
        end if

        if (it_final < ntimstep) then
            do it = it_final + 1, ntimstep
                orbit_traj(:, it) = ieee_value(0.0d0, ieee_quiet_nan)
                orbit_times(it) = ieee_value(0.0d0, ieee_quiet_nan)
            end do
        end if

        if (ierr_orbit == SYMPL_SPECTRE_OK) then
            t_stop = real(kt, dp)*dtaumin/v0
            orbit_exit_code(ipart) = ORBIT_EXIT_COMPLETED
        else if (ierr_orbit == SYMPL_SPECTRE_SKIM) then
            ! Mirror-confined at an interior interface: cannot be lost, so record
            ! it as confined (times_lost = trace_time) rather than at its stop.
            t_stop = trace_time
        else if (ierr_orbit /= SYMPL_SPECTRE_LOSS) then
            t_stop = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            t_stop = (real(kt, dp) + t_frac)*dtaumin/v0
        end if

        if (state%fo%active .and. state%fo%has_y) then
            z = state%fo%last_y
        else
            call to_standard_z_coordinates(anorb, z)
        end if
        ! The next marker on this thread starts with an unlocked field dispatch.
        call set_spectre_volume_lock(0)

        !$omp critical
        call integ_to_ref(z(1:3), zend(1:3, ipart))
        zend(4:5, ipart) = z(4:5)
        times_lost(ipart) = t_stop
        !$omp end critical
    end subroutine trace_orbit_spectre_sympl

    pure subroutine locate_linear_lcfs(z_start, z_end, field_period, z_event, &
            event_fraction)
        real(dp), intent(in) :: z_start(5), z_end(5), field_period
        real(dp), intent(out) :: z_event(5), event_fraction
        real(dp) :: theta_delta, phi_delta

        event_fraction = 1d0
        if (z_end(1) > z_start(1)) then
            event_fraction = min(1d0, max(0d0, &
                (1d0 - z_start(1))/(z_end(1) - z_start(1))))
        end if
        theta_delta = modulo(z_end(2) - z_start(2) + 0.5d0*twopi, twopi) - &
            0.5d0*twopi
        phi_delta = modulo(z_end(3) - z_start(3) + 0.5d0*field_period, &
            field_period) - 0.5d0*field_period
        z_event(1) = 1d0
        z_event(2) = z_start(2) + event_fraction*theta_delta
        z_event(3) = z_start(3) + event_fraction*phi_delta
        z_event(4:5) = z_start(4:5) + &
            event_fraction*(z_end(4:5) - z_start(4:5))
    end subroutine locate_linear_lcfs

    subroutine locate_validated_lcfs(z_start, z_end, field_period, z_event, &
            event_fraction, ierr)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        real(dp), intent(in) :: z_start(5), z_end(5), field_period
        real(dp), intent(out) :: z_event(5), event_fraction
        integer, intent(inout) :: ierr
        real(dp), parameter :: radial_sanity_band = 0.05_dp
        real(dp), parameter :: max_local_radial_step = 0.5_dp
        real(dp) :: u_event(3)

        if (.not. all(ieee_is_finite(z_start)) .or. &
            .not. all(ieee_is_finite(z_end)) .or. &
            abs(z_end(1) - z_start(1)) > max_local_radial_step) then
            z_event = z_start
            event_fraction = 0.0_dp
            ierr = 2
            return
        end if

        call locate_linear_lcfs(z_start, z_end, field_period, z_event, &
            event_fraction)
        call integ_to_ref(z_event(1:3), u_event)
        if (all(ieee_is_finite(z_event)) .and. &
            all(ieee_is_finite(u_event)) .and. &
            abs(u_event(1) - 1.0_dp) <= radial_sanity_band) return

        ! The legacy RK chamber flag is based on its first integration
        ! coordinate. In a canonical map, a gross numerical jump can therefore
        ! resemble an LCFS crossing. Do not turn that failure into a physical
        ! loss unless the candidate is local and maps back to a finite LCFS point.
        z_event = z_start
        event_fraction = 0.0_dp
        ierr = 2
    end subroutine locate_validated_lcfs

    subroutine validate_rk_state(z_start, z_end, ierr, expected_momentum)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        real(dp), intent(in) :: z_start(5)
        real(dp), intent(inout) :: z_end(5)
        integer, intent(inout) :: ierr
        real(dp), intent(in), optional :: expected_momentum
        real(dp), parameter :: radial_sanity_band = 0.05_dp
        real(dp), parameter :: max_local_radial_step = 0.5_dp
        real(dp), parameter :: momentum_sanity_fraction = 0.005_dp
        real(dp) :: momentum_reference
        real(dp) :: u_end(3)
        logical :: momentum_valid

        if (ierr /= 0) return
        momentum_reference = z_start(4)
        momentum_valid = .true.
        if (present(expected_momentum)) then
            momentum_reference = expected_momentum
            momentum_valid = abs(z_end(4) - momentum_reference) <= &
                momentum_sanity_fraction* &
                max(abs(momentum_reference), tiny(1.0_dp))
        end if
        if (all(ieee_is_finite(z_end)) .and. z_end(4) > 0.0_dp .and. &
            abs(z_end(5)) <= 1.0_dp + 100.0_dp*epsilon(1.0_dp) .and. &
            momentum_valid .and. &
            abs(z_end(1) - z_start(1)) <= max_local_radial_step) then
            call integ_to_ref(z_end(1:3), u_end)
            if (all(ieee_is_finite(u_end)) .and. &
                u_end(1) >= -radial_sanity_band .and. &
                u_end(1) <= 1.0_dp + radial_sanity_band) return
        end if

        z_end = z_start
        ierr = 2
    end subroutine validate_rk_state

    subroutine macrostep(anorb, z, kt, ierr_orbit, ntau_local, exit_step, &
            hold_streak, numerical_hold_any, rk_recovery)
        use alpha_lifetime_sub, only: orbit_timestep_axis
        use orbit_symplectic, only: advance_symplectic_with_retry, &
            orbit_timestep_sympl

        type(tracer_t), intent(inout) :: anorb
        real(dp), intent(inout) :: z(5)
        integer(8), intent(inout) :: kt
        integer, intent(out) :: ierr_orbit
        integer, intent(in) :: ntau_local
        real(dp), intent(out), optional :: exit_step
        integer, intent(inout), optional :: hold_streak
        logical, intent(out), optional :: numerical_hold_any
        type(rk_recovery_state_t), intent(inout), optional :: rk_recovery

        integer :: hold_streak_local, ktau, recovery_method
        real(dp) :: z_step_start(5), z_step_end(5), loss_fraction
        logical :: numerical_hold, numerical_hold_any_local
        type(rk_recovery_state_t) :: rk_recovery_local
        type(symplectic_integrator_t) :: si_step_start
        type(field_can_t) :: f_step_start

        if (present(exit_step)) exit_step = real(kt, dp)
        hold_streak_local = 0
        if (present(hold_streak)) hold_streak_local = hold_streak
        numerical_hold_any_local = .false.
        rk_recovery_local = rk_recovery_state_t()
        if (present(rk_recovery)) rk_recovery_local = rk_recovery

        do ktau = 1, ntau_local
            numerical_hold = .false.
            z_step_start = z
            if (orbit_model == ORBIT_FULL_ORBIT) then
                call orbit_timestep_fo(anorb%fo, z, ierr_orbit)
                if (ierr_orbit .ne. 0) then
                    if (ierr_orbit == ORBIT_FO_LOSS) then
                        z_step_end = z
                        call locate_linear_lcfs(z_step_start, z_step_end, &
                            anorb%fper, &
                            z, loss_fraction)
                        if (present(exit_step)) then
                            exit_step = real(kt, dp) + loss_fraction
                        end if
                    else if (symplectic_newton_warning_mode) then
                        call hold_isolated_warning_failure(z, z_step_start, &
                            ierr_orbit, hold_streak_local, numerical_hold, &
                            numerical_hold_any_local)
                    end if
                    if (ierr_orbit .ne. 0) exit
                end if
                if (numerical_hold) then
                    kt = kt + 1
                    cycle
                end if
                hold_streak_local = 0
                kt = kt + 1
                cycle
            end if
            if (integmode <= 0 .or. rk_recovery_local%active) then
                if (rk_recovery_local%active) then
                    call orbit_timestep_recovery(anorb%fo, z, dtaumin, relerr, &
                        ierr_orbit, recovery_method)
                else
                    call orbit_timestep_axis(z, dtaumin, dtaumin, relerr, &
                        ierr_orbit)
                end if
                if (rk_recovery_local%has_momentum_reference) then
                    call validate_rk_state(z_step_start, z, ierr_orbit, &
                        rk_recovery_local%momentum_reference)
                else
                    call validate_rk_state(z_step_start, z, ierr_orbit, &
                        z_step_start(4))
                end if
                if (ierr_orbit == 1) then
                    z_step_end = z
                    call locate_validated_lcfs(z_step_start, z_step_end, &
                        anorb%fper, z, loss_fraction, ierr_orbit)
                end if
                if (ierr_orbit == 1) then
                    if (present(exit_step)) then
                        exit_step = real(kt, dp) + loss_fraction
                    end if
                    if (rk_recovery_local%active) then
                        anorb%si%last_step_fraction = loss_fraction
                        anorb%si%last_event_radial_residual = abs(z(1) - 1.d0)
                        anorb%si%last_event_fraction_width = 1.d0
                        ierr_orbit = SYMPLECTIC_STEP_BOUNDARY
                    end if
                else if (ierr_orbit /= 0) then
                    call hold_isolated_warning_failure(z, z_step_start, &
                        ierr_orbit, hold_streak_local, numerical_hold, &
                        numerical_hold_any_local)
                end if
            else
                if (swcoll) call update_momentum(anorb, z)
                si_step_start = anorb%si
                f_step_start = anorb%f
                call advance_symplectic_with_retry(anorb%si, anorb%f, &
                    orbit_timestep_sympl, ierr_orbit)
                if (ierr_orbit == 0) then
                    call to_standard_z_coordinates(anorb, z)
                    call validate_rk_state(z_step_start, z, ierr_orbit, &
                        si_step_start%pabs)
                    if (ierr_orbit /= 0) then
                        anorb%si = si_step_start
                        anorb%f = f_step_start
                    end if
                end if
                if (ierr_orbit == SYMPLECTIC_STEP_BOUNDARY) then
                    call to_standard_z_coordinates(anorb, z)
                    if (present(exit_step)) exit_step = real(kt, dp) + &
                        anorb%si%last_step_fraction
                    exit
                end if
                if (ierr_orbit .ne. 0 .and. &
                    symplectic_newton_warning_mode) then
                    ! A failed symplectic map has restored its last accepted
                    ! state. Advance that same interval with the established
                    ! adaptive RK/axis path, then thread-locally reseed the
                    ! symplectic state. Keep this marker on RK while it remains
                    ! in the failing axis/edge region; generic failures use an
                    ! exponentially backed-off probe interval.
                    z = z_step_start
                    z(4) = si_step_start%pabs
                    call orbit_timestep_recovery(anorb%fo, z, dtaumin, relerr, &
                        ierr_orbit, recovery_method)
                    call validate_rk_state(z_step_start, z, ierr_orbit, &
                        si_step_start%pabs)
                    if (ierr_orbit == 1) then
                        z_step_end = z
                        call locate_validated_lcfs(z_step_start, z_step_end, &
                            anorb%fper, z, loss_fraction, ierr_orbit)
                    end if
                    if (ierr_orbit == 0) then
                        call reseed_sympl(anorb%si, anorb%f, z, &
                            si_step_start%pabs)
                        call count_event(EVT_SYMPLECTIC_RK_RECOVERY)
                        call activate_rk_recovery(rk_recovery_local, &
                            z_step_start, si_step_start%pabs)
                        hold_streak_local = 0
                    else if (ierr_orbit == 1) then
                        anorb%si%last_step_fraction = loss_fraction
                        anorb%si%last_event_radial_residual = abs(z(1) - 1.d0)
                        anorb%si%last_event_fraction_width = 1.d0
                        if (present(exit_step)) exit_step = real(kt, dp) + &
                            loss_fraction
                        ierr_orbit = SYMPLECTIC_STEP_BOUNDARY
                        exit
                    else
                        z = z_step_start
                        anorb%si = si_step_start
                        anorb%f = f_step_start
                        call hold_isolated_warning_failure(z, z_step_start, &
                            ierr_orbit, hold_streak_local, numerical_hold, &
                            numerical_hold_any_local)
                    end if
                else if (ierr_orbit == 0) then
                    call to_standard_z_coordinates(anorb, z)
                    hold_streak_local = 0
                    rk_recovery_local = rk_recovery_state_t()
                end if
                if (ierr_orbit .ne. 0) exit
            end if
            if (ierr_orbit .ne. 0) exit
            if (numerical_hold) then
                kt = kt + 1
                cycle
            end if
            if (rk_recovery_local%active .and. .not. numerical_hold) then
                hold_streak_local = 0
                call update_recovery_streak(rk_recovery_local, z, &
                    recovery_method)
                if (should_resume_symplectic(rk_recovery_local, z)) then
                    if (rk_recovery_local%has_momentum_reference) then
                        call reseed_sympl(anorb%si, anorb%f, z, &
                            rk_recovery_local%momentum_reference)
                    else
                        call reseed_sympl(anorb%si, anorb%f, z)
                    end if
                    rk_recovery_local%active = .false.
                    call count_event(EVT_SYMPLECTIC_RESUME)
                end if
            end if
            if (integmode <= 0 .or. rk_recovery_local%active) &
                hold_streak_local = 0
            if (swcoll) then
                call collide(z, dtaumin) ! Collisions
                if (rk_recovery_local%active) then
                    rk_recovery_local%momentum_reference = z(4)
                end if
            end if
            kt = kt + 1
        end do
        if (present(hold_streak)) hold_streak = hold_streak_local
        if (present(numerical_hold_any)) &
            numerical_hold_any = numerical_hold_any_local
        if (present(rk_recovery)) rk_recovery = rk_recovery_local
    end subroutine macrostep

    subroutine hold_isolated_warning_failure(z, z_step_start, ierr_orbit, &
            hold_streak, numerical_hold, numerical_hold_any)
        real(dp), intent(inout) :: z(5)
        real(dp), intent(in) :: z_step_start(5)
        integer, intent(inout) :: ierr_orbit, hold_streak
        logical, intent(inout) :: numerical_hold, numerical_hold_any

        z = z_step_start
        if (.not. symplectic_newton_warning_mode) return
        if (hold_streak >= max_consecutive_warning_holds) return

        ! Preserve the historical warning-level behavior for one isolated
        ! unresolved microstep. Unlike the former permanent latch, the next
        ! microstep retries both pushers; a successful step resets the streak.
        call count_event(EVT_WARNING_STEP_SKIP)
        hold_streak = hold_streak + 1
        numerical_hold = .true.
        numerical_hold_any = .true.
        ierr_orbit = 0
    end subroutine hold_isolated_warning_failure

    subroutine locate_wall_segment(z, z_start, z_end, ipart, x_start_m, &
            u_start, x_end_m, u_end, field_period, start_step, segment_steps, &
            hit, hit_step)
        real(dp), intent(inout) :: z(5)
        real(dp), intent(in) :: z_start(5), z_end(5)
        integer, intent(in) :: ipart
        real(dp), intent(in) :: x_start_m(3), u_start(3), x_end_m(3), u_end(3)
        real(dp), intent(in) :: field_period
        real(dp), intent(in) :: start_step, segment_steps
        logical, intent(out) :: hit
        real(dp), intent(out) :: hit_step

        real(dp) :: x_hit_m(3), x_hit(3), normal_m(3), vhat(3)
        real(dp) :: segment_length, hit_distance, hit_fraction, cos_inc
        real(dp) :: u_hit(3), theta_delta, phi_delta
        integer :: ierr_from_cart

        hit = .false.
        hit_step = start_step + segment_steps
        if (wall_query_rho_lcfs >= 0.0_dp .and. &
            u_end(1) <= wall_query_rho_lcfs) return
        call stl_wall_first_hit_segment_with_normal( &
            wall, x_start_m, x_end_m, hit, x_hit_m, normal_m)
        if (.not. hit) return

        x_hit = x_hit_m/chartmap_cart_scale_to_m
        wall_hit(ipart) = 1_int8
        wall_hit_cart(:, ipart) = x_hit
        wall_hit_normal_cart(:, ipart) = normal_m

        vhat = x_end_m - x_start_m
        segment_length = sqrt(sum(vhat*vhat))
        hit_fraction = 1d0
        if (segment_length > 0d0) then
            hit_distance = sqrt(sum((x_hit_m - x_start_m)**2))
            hit_fraction = min(1d0, max(0d0, hit_distance/segment_length))
            vhat = vhat/segment_length
            cos_inc = min(1d0, max(0d0, abs(sum(vhat*normal_m))))
            wall_hit_cos_incidence(ipart) = cos_inc
            wall_hit_angle_rad(ipart) = acos(cos_inc)
        end if

        ierr_from_cart = 0
        select type (ccs => ref_coords)
        class is (chartmap_coordinate_system_t)
            call ccs%from_cart(x_hit, u_hit, ierr_from_cart)
        class default
            ierr_from_cart = 1
        end select

        if (ierr_from_cart /= 0) then
            theta_delta = modulo(u_end(2) - u_start(2) + 0.5d0*twopi, &
                twopi) - 0.5d0*twopi
            phi_delta = modulo(u_end(3) - u_start(3) + 0.5d0*field_period, &
                field_period) - 0.5d0*field_period
            u_hit(1) = u_start(1) + hit_fraction*(u_end(1) - u_start(1))
            u_hit(2) = u_start(2) + hit_fraction*theta_delta
            u_hit(3) = u_start(3) + hit_fraction*phi_delta
        end if
        call ref_to_integ(u_hit, z(1:3))
        z(4:5) = z_start(4:5) + hit_fraction*(z_end(4:5) - z_start(4:5))
        hit_step = start_step + hit_fraction*segment_steps
    end subroutine locate_wall_segment

    subroutine macrostep_with_wall_check(anorb, z, kt, ierr_orbit, ntau_local, &
            ipart, x_prev_m, exit_step, hold_streak, numerical_hold_any, &
            rk_recovery)
        use alpha_lifetime_sub, only: orbit_timestep_axis
        use orbit_symplectic, only: advance_symplectic_with_retry, &
            orbit_timestep_sympl

        type(tracer_t), intent(inout) :: anorb
        real(dp), intent(inout) :: z(5)
        integer(8), intent(inout) :: kt
        integer, intent(out) :: ierr_orbit
        integer, intent(in) :: ntau_local
        integer, intent(in) :: ipart
        real(dp), intent(inout) :: x_prev_m(3)
        real(dp), intent(out), optional :: exit_step
        integer, intent(inout), optional :: hold_streak
        logical, intent(out), optional :: numerical_hold_any
        type(rk_recovery_state_t), intent(inout), optional :: rk_recovery

        integer :: hold_streak_local, ktau, recovery_method
        real(dp) :: u_ref_prev(3), u_ref_cur(3), x_cur(3), x_cur_m(3)
        real(dp) :: z_step_start(5), z_step_end(5), segment_duration
        real(dp) :: boundary_fraction
        real(dp) :: wall_exit_step
        logical :: hit, numerical_hold, numerical_hold_any_local
        type(rk_recovery_state_t) :: rk_recovery_local
        type(symplectic_integrator_t) :: si_step_start
        type(field_can_t) :: f_step_start

        call integ_to_ref(z(1:3), u_ref_prev)
        if (present(exit_step)) exit_step = real(kt, dp)
        hold_streak_local = 0
        if (present(hold_streak)) hold_streak_local = hold_streak
        numerical_hold_any_local = .false.
        rk_recovery_local = rk_recovery_state_t()
        if (present(rk_recovery)) rk_recovery_local = rk_recovery
        do ktau = 1, ntau_local
            numerical_hold = .false.
            segment_duration = 1.0_dp
            z_step_start = z
            if (integmode <= 0 .or. rk_recovery_local%active) then
                if (rk_recovery_local%active) then
                    call orbit_timestep_recovery(anorb%fo, z, dtaumin, relerr, &
                        ierr_orbit, recovery_method)
                else
                    call orbit_timestep_axis(z, dtaumin, dtaumin, relerr, &
                        ierr_orbit)
                end if
                if (rk_recovery_local%has_momentum_reference) then
                    call validate_rk_state(z_step_start, z, ierr_orbit, &
                        rk_recovery_local%momentum_reference)
                else
                    call validate_rk_state(z_step_start, z, ierr_orbit, &
                        z_step_start(4))
                end if
                if (ierr_orbit == 1) then
                    z_step_end = z
                    call locate_validated_lcfs(z_step_start, z_step_end, &
                        anorb%fper, z, boundary_fraction, ierr_orbit)
                end if
                if (ierr_orbit == 1) then
                    z_step_end = z
                    call integ_to_ref(z(1:3), u_ref_cur)
                    call ref_coords%evaluate_cart(u_ref_cur, x_cur)
                    x_cur_m = x_cur*chartmap_cart_scale_to_m
                    call locate_wall_segment(z, z_step_start, z_step_end, ipart, &
                        x_prev_m, u_ref_prev, x_cur_m, u_ref_cur, anorb%fper, &
                        real(kt, dp), boundary_fraction, hit, wall_exit_step)
                    if (hit) ierr_orbit = 77
                    if (.not. hit .and. rk_recovery_local%active) then
                        ierr_orbit = SYMPLECTIC_STEP_BOUNDARY
                        anorb%si%last_step_fraction = boundary_fraction
                        anorb%si%last_event_radial_residual = abs(z(1) - 1.d0)
                        anorb%si%last_event_fraction_width = 1.d0
                    end if
                    if (present(exit_step)) then
                        if (hit) then
                            exit_step = wall_exit_step
                        else
                            exit_step = real(kt, dp) + boundary_fraction
                        end if
                    end if
                    exit
                else if (ierr_orbit /= 0) then
                    call hold_isolated_warning_failure(z, z_step_start, &
                        ierr_orbit, hold_streak_local, numerical_hold, &
                        numerical_hold_any_local)
                end if
            else
                if (swcoll) call update_momentum(anorb, z)
                si_step_start = anorb%si
                f_step_start = anorb%f
                call advance_symplectic_with_retry(anorb%si, anorb%f, &
                    orbit_timestep_sympl, ierr_orbit, segment_duration)
                if (ierr_orbit == 0) then
                    call to_standard_z_coordinates(anorb, z)
                    call validate_rk_state(z_step_start, z, ierr_orbit, &
                        si_step_start%pabs)
                    if (ierr_orbit /= 0) then
                        anorb%si = si_step_start
                        anorb%f = f_step_start
                    end if
                end if
                if (ierr_orbit == SYMPLECTIC_STEP_BOUNDARY) then
                    call to_standard_z_coordinates(anorb, z)
                    call integ_to_ref(z(1:3), u_ref_cur)
                    call ref_coords%evaluate_cart(u_ref_cur, x_cur)
                    x_cur_m = x_cur*chartmap_cart_scale_to_m
                    z_step_end = z
                    segment_duration = anorb%si%last_step_fraction
                    call locate_wall_segment(z, z_step_start, z_step_end, ipart, &
                        x_prev_m, u_ref_prev, x_cur_m, u_ref_cur, anorb%fper, &
                        real(kt, dp), segment_duration, hit, wall_exit_step)
                    if (hit) ierr_orbit = 77
                    if (present(exit_step)) then
                        if (hit) then
                            exit_step = wall_exit_step
                        else
                            exit_step = real(kt, dp) + segment_duration
                        end if
                    end if
                    exit
                end if
                if (ierr_orbit .ne. 0 .and. &
                    symplectic_newton_warning_mode) then
                    z = z_step_start
                    z(4) = si_step_start%pabs
                    call orbit_timestep_recovery(anorb%fo, z, dtaumin, relerr, &
                        ierr_orbit, recovery_method)
                    call validate_rk_state(z_step_start, z, ierr_orbit, &
                        si_step_start%pabs)
                    if (ierr_orbit == 1) then
                        z_step_end = z
                        call locate_validated_lcfs(z_step_start, z_step_end, &
                            anorb%fper, z, boundary_fraction, ierr_orbit)
                    end if
                    if (ierr_orbit == 0) then
                        call reseed_sympl(anorb%si, anorb%f, z, &
                            si_step_start%pabs)
                        call count_event(EVT_SYMPLECTIC_RK_RECOVERY)
                        call activate_rk_recovery(rk_recovery_local, &
                            z_step_start, si_step_start%pabs)
                        hold_streak_local = 0
                    else if (ierr_orbit == 1) then
                        call integ_to_ref(z(1:3), u_ref_cur)
                        call ref_coords%evaluate_cart(u_ref_cur, x_cur)
                        x_cur_m = x_cur*chartmap_cart_scale_to_m
                        z_step_end = z
                        call locate_wall_segment(z, z_step_start, z_step_end, &
                            ipart, &
                            x_prev_m, u_ref_prev, x_cur_m, u_ref_cur, anorb%fper, &
                            real(kt, dp), boundary_fraction, hit, wall_exit_step)
                        if (hit) then
                            ierr_orbit = 77
                        else
                            ierr_orbit = SYMPLECTIC_STEP_BOUNDARY
                            anorb%si%last_step_fraction = boundary_fraction
                            anorb%si%last_event_radial_residual = abs(z(1) - 1.d0)
                            anorb%si%last_event_fraction_width = 1.d0
                        end if
                        if (present(exit_step)) then
                            if (hit) then
                                exit_step = wall_exit_step
                            else
                                exit_step = real(kt, dp) + boundary_fraction
                            end if
                        end if
                        exit
                    else
                        z = z_step_start
                        anorb%si = si_step_start
                        anorb%f = f_step_start
                        call hold_isolated_warning_failure(z, z_step_start, &
                            ierr_orbit, hold_streak_local, numerical_hold, &
                            numerical_hold_any_local)
                    end if
                else if (ierr_orbit == 0) then
                    call to_standard_z_coordinates(anorb, z)
                    hold_streak_local = 0
                    rk_recovery_local = rk_recovery_state_t()
                end if
                if (ierr_orbit .ne. 0) exit
            end if
            if (ierr_orbit .ne. 0) exit
            if (numerical_hold) then
                kt = kt + 1
                cycle
            end if
            z_step_end = z
            call integ_to_ref(z(1:3), u_ref_cur)
            call ref_coords%evaluate_cart(u_ref_cur, x_cur)
            x_cur_m = x_cur*chartmap_cart_scale_to_m
            call locate_wall_segment(z, z_step_start, z_step_end, ipart, &
                x_prev_m, u_ref_prev, x_cur_m, u_ref_cur, anorb%fper, &
                real(kt, dp), segment_duration, hit, wall_exit_step)
            if (hit) then
                ierr_orbit = 77
                if (present(exit_step)) exit_step = wall_exit_step
                exit
            end if

            if (swcoll) then
                call collide(z, dtaumin)
                if (rk_recovery_local%active) then
                    rk_recovery_local%momentum_reference = z(4)
                end if
            end if
            kt = kt + 1
            if (rk_recovery_local%active) then
                hold_streak_local = 0
                call update_recovery_streak(rk_recovery_local, z, &
                    recovery_method)
                if (should_resume_symplectic(rk_recovery_local, z)) then
                    if (rk_recovery_local%has_momentum_reference) then
                        call reseed_sympl(anorb%si, anorb%f, z, &
                            rk_recovery_local%momentum_reference)
                    else
                        call reseed_sympl(anorb%si, anorb%f, z)
                    end if
                    rk_recovery_local%active = .false.
                    call count_event(EVT_SYMPLECTIC_RESUME)
                end if
            end if
            if (integmode <= 0 .or. rk_recovery_local%active) &
                hold_streak_local = 0
            x_prev_m = x_cur_m
            u_ref_prev = u_ref_cur
        end do
        if (present(hold_streak)) hold_streak = hold_streak_local
        if (present(numerical_hold_any)) &
            numerical_hold_any = numerical_hold_any_local
        if (present(rk_recovery)) rk_recovery = rk_recovery_local
    end subroutine macrostep_with_wall_check

    subroutine to_standard_z_coordinates(anorb, z)
        type(tracer_t), intent(in) :: anorb
        real(dp), intent(inout) :: z(5)

        z(1:3) = anorb%si%z(1:3)
        z(4) = dsqrt(anorb%f%mu*anorb%f%Bmod + 0.5d0*anorb%f%vpar**2)
        z(5) = anorb%f%vpar/(z(4)*sqrt2)
    end subroutine to_standard_z_coordinates

    subroutine increase_confined_count(it, passing)
        integer, intent(in) :: it
        logical, intent(in) :: passing

        if (passing) then
            !$omp atomic update
            confpart_pass(it) = confpart_pass(it) + 1.d0
        else
            !$omp atomic update
            confpart_trap(it) = confpart_trap(it) + 1.d0
        end if
    end subroutine increase_confined_count

    subroutine compute_pitch_angle_params(z, passing, trap_par_, perp_inv_)
        use find_bminmax_sub, only: get_bminmax

        real(dp), intent(in) :: z(5)
        logical, intent(out) :: passing
        real(dp), intent(out) :: trap_par_, perp_inv_

        real(dp) :: bmod

        !$omp critical
        bmod = compute_bmod(z(1:3))
        if (num_surf /= 1) then
            call get_bminmax(z(1), bmin, bmax)
        end if
        passing = z(5)**2 .gt. 1.d0 - bmod/bmax
        trap_par_ = ((1.d0 - z(5)**2)*bmax/bmod - 1.d0)*bmin/(bmax - bmin)
        perp_inv_ = z(4)**2*(1.d0 - z(5)**2)/bmod
        !$omp end critical
    end subroutine compute_pitch_angle_params

    function compute_bmod(z) result(bmod)
        use magfie_sub, only: magfie

        real(dp) :: bmod
        real(dp), intent(in) :: z(3)

        real(dp) :: sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        call magfie(z(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    end function compute_bmod

    subroutine update_momentum(anorb, z)
        use orbit_symplectic, only: get_val

        type(tracer_t), intent(inout) :: anorb
        real(dp), intent(in) :: z(5)

        anorb%si%pabs = z(4)
        anorb%f%vpar = z(4)*z(5)*sqrt2
        anorb%f%mu = z(4)**2*(1.d0 - z(5)**2)/anorb%f%Bmod
        anorb%si%z(4) = anorb%f%vpar*anorb%f%hph + anorb%f%Aph/anorb%f%ro0
        call get_val(anorb%f, anorb%si%z(4)) ! for pth
        anorb%si%pthold = anorb%f%pth
    end subroutine update_momentum

    subroutine collide(z, dt)
        real(dp), intent(inout) :: z(5)
        real(dp), intent(in) :: dt
        integer :: ierr_coll

        call stost(z, dt, 1, ierr_coll)
        if (ierr_coll /= 0) then
            print *, 'Error in stost: ', ierr_coll, 'z = ', z, 'dtaumin = ', dtaumin
        end if
    end subroutine collide

    subroutine write_output
        use field_can_base, only: n_field_evaluations
        use params, only: output_results_netcdf
        use netcdf_results_output, only: write_results_netcdf

        integer(8) :: total_field_evaluations

        ! Sum field evaluations across all threads
        total_field_evaluations = 0
        !$omp parallel reduction(+:total_field_evaluations)
        total_field_evaluations = total_field_evaluations + n_field_evaluations
        !$omp end parallel

        print *, "Total field evaluations: ", total_field_evaluations

        call write_results

        if (output_results_netcdf) then
            call write_results_netcdf('results.nc')
        end if
    end subroutine write_output

    subroutine write_results
        !> Write the per-particle and confined-fraction result files from the
        !> shared result arrays after particle tracing has reached a team-safe
        !> point. Confined fractions are conditional on the numerically resolved
        !> population at each time. Particles not yet finished (for example in an
        !> explicitly requested partial write) keep times_lost = -1.
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

        integer :: i, num_lost, resolved_count, unit
        real(dp) :: inverse_times_lost_sum, norm, pass_fraction, trap_fraction

        norm = real(max(ntestpart, 1), dp)

        open (newunit=unit, file='orbit_exit_code.dat', status='replace', &
            action='write')
        close (unit)

        open (newunit=unit, file='times_lost.dat', status='replace', &
            action='write', recl=1024)
        num_lost = 0
        inverse_times_lost_sum = 0.0d0
        do i = 1, ntestpart
            if (orbit_exit_code(i) >= ORBIT_EXIT_NUMERICAL_DOMAIN) then
                times_lost(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
            write (unit, *) i, times_lost(i), trap_par(i), zstart(1, i), &
                perp_inv(i), zend(:, i)
            if (orbit_exit_code(i) == ORBIT_EXIT_LCFS .or. &
                orbit_exit_code(i) == ORBIT_EXIT_WALL) then
                if (times_lost(i) > 0d0) then
                    num_lost = num_lost + 1
                    inverse_times_lost_sum = inverse_times_lost_sum + &
                        1.0/times_lost(i)
                end if
            end if
        end do
        close (unit)

        open (newunit=unit, file='orbit_exit_code.dat', status='replace', &
            action='write', recl=1024)
        do i = 1, ntestpart
            write (unit, *) i, orbit_exit_code(i), &
                times_lost(i), boundary_event_radial_residual(i), &
                boundary_event_time_width(i)
        end do
        close (unit)

        open (newunit=unit, file='avg_inverse_t_lost.dat', status='replace', &
            action='write', recl=1024)
        if (num_lost > 0) then
            write (unit, *) inverse_times_lost_sum/num_lost
        else
            write (unit, *) 0d0
        end if
        close (unit)

        open (newunit=unit, file='confined_fraction.dat', status='replace', &
            action='write', recl=1024)
        do i = 1, ntimstep
            resolved_count = ntestpart - unresolved_orbits(i)
            if (resolved_count > 0) then
                pass_fraction = confpart_pass(i)/real(resolved_count, dp)
                trap_fraction = confpart_trap(i)/real(resolved_count, dp)
            else
                pass_fraction = ieee_value(0.0_dp, ieee_quiet_nan)
                trap_fraction = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
            write (unit, *) dble(kt_macro(i))*dtaumin/v0, pass_fraction, &
                trap_fraction, resolved_count
        end do
        close (unit)

        open (newunit=unit, file='unresolved_fraction.dat', status='replace', &
            action='write', recl=1024)
        do i = 1, ntimstep
            write (unit, *) dble(kt_macro(i))*dtaumin/v0, &
                real(unresolved_orbits(i), dp)/norm, ntestpart
        end do
        close (unit)

        call write_spectre_crossing_events

        if (ntcut > 0 .or. class_plot) then
            open (newunit=unit, file='class_parts.dat', recl=1024)
            do i = 1, ntestpart
                write (unit, *) i, zstart(1, i), perp_inv(i), iclass(:, i)
            end do
            close (unit)

            if (needs_bminmax_cache()) then
                block
                    use bminmax_mod, only: nsbmnx, hsbmnx, bmin_arr, bmax_arr

                    open (newunit=unit, file='bminmax.dat', recl=1024)
                    do i = 0, nsbmnx
                        write (unit, *) hsbmnx*dble(i), bmin_arr(i), bmax_arr(i)
                    end do
                    close (unit)
                end block
            end if
        end if
    end subroutine write_results

    subroutine write_spectre_crossing_events
        !> Interface crossing log (#443/#440), one line per crossing/reflection
        !> event across all SPECTRE markers, grouped by particle in time order.
        !> For symplectic runs the exact-landing statistics (#441) are printed so
        !> tests can assert the interface landing accuracy from stdout.
        use interface_crossing, only: crossing_log_write
        use spectre_sympl_orbit, only: sympl_landing_stats, sympl_sheet_stats, &
            sympl_fo_stats
        use magfie_sub, only: SPECTRE

        integer :: landings, stops, sheet_entries, sheet_exits
        integer :: sheet_init_failures, sheet_advance_failures
        integer :: sheet_failure_status(5)
        integer :: sheet_stop_reason(5)
        integer :: fo_entries, fo_exits, fo_losses, fo_failures
        integer :: fo_failure_status(5)
        real(dp) :: max_resid

        if (isw_field_type /= SPECTRE) return

        call crossing_log_write('spectre_crossing_events.dat')

        if (integmode > 0) then
            call sympl_landing_stats(landings, max_resid, stops)
            print '(A,I0,A,ES12.4,A,I0)', 'sympl_landing_stats: count= ', &
                landings, ' max_resid= ', max_resid, ' cross_stop= ', stops
            call sympl_sheet_stats(sheet_entries, sheet_exits, &
                sheet_init_failures, sheet_advance_failures, sheet_failure_status, &
                sheet_stop_reason)
            print '(A,I0,A,I0,A,I0,A,I0,A,5(I0,1X),A,5(I0,1X))', &
                'sympl_sheet_stats: entries= ', &
                sheet_entries, ' exits= ', sheet_exits, ' init_fail= ', &
                sheet_init_failures, ' advance_fail= ', sheet_advance_failures, &
                ' status= ', sheet_failure_status, ' stop_reason= ', &
                sheet_stop_reason
            call sympl_fo_stats(fo_entries, fo_exits, fo_losses, fo_failures, &
                fo_failure_status)
            print '(A,I0,A,I0,A,I0,A,I0,A,5(I0,1X))', &
                'sympl_fo_stats: entries= ', fo_entries, ' exits= ', fo_exits, &
                ' losses= ', fo_losses, ' failures= ', fo_failures, ' status= ', &
                fo_failure_status
        end if
    end subroutine write_spectre_crossing_events

end module simple_main
