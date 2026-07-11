module simple_main
    use, intrinsic :: iso_fortran_env, only: int8
    use omp_lib
    use util, only: sqrt2
    use simple, only: init_vmec, init_sympl, init_fo, orbit_timestep_fo, tracer_t
    use diag_mod, only: icounter
    use collis_alp, only: loacol_alpha, stost, init_collision_profiles
    use samplers, only: sample
    use field_can_mod, only: integ_to_ref, ref_to_integ, init_field_can
	    use callback, only: output_orbits_macrostep
	    use params, only: swcoll, ntestpart, startmode, special_ants_file, num_surf, &
	                      grid_density, dtau, dtaumin, ntau, v0, &
	                      kpart, confpart_pass, confpart_trap, times_lost, integmode, &
	                      relerr, trace_time, &
	                      class_plot, fast_class, generate_start_only, ntcut, iclass, &
	                      bmin, bmax, zstart, zend, trap_par, perp_inv, sbeg, &
	                      ntimstep, should_skip, reset_seed_if_deterministic, &
	                      field_input, isw_field_type, reuse_batch, coord_input, &
	                      wall_input, wall_units, wall_hit, wall_hit_cart, &
	                      wall_hit_normal_cart, wall_hit_cos_incidence, &
	                      wall_hit_angle_rad, ntau_macro, kt_macro, &
	                      checkpoint_interval, orbit_model, orbit_coord, &
	                      ORBIT_GC, ORBIT_FULL_ORBIT
    use diag_counters, only: diag_counters_init
    use progress_monitor, only: progress_init, progress_tick, progress_finalize
    use restart_mod, only: particle_done, read_restart_data, restore_confined_counts
    use chartmap_metadata, only: chartmap_metadata_t, read_chartmap_metadata
    use reference_coordinates, only: ref_coords
    use stl_wall_intersection, only: stl_wall_t, stl_wall_init, &
                                     stl_wall_finalize, &
                                     stl_wall_first_hit_segment_with_normal
    use libneo_coordinates, only: chartmap_coordinate_system_t

    implicit none

    ! Define real(dp) kind parameter
    integer, parameter :: dp = kind(1.0d0)

    logical, save :: wall_enabled = .false.
    real(dp), save :: wall_rho_lcfs = -1.0d0
    real(dp), save :: chartmap_cart_scale_to_m = -1.0d0
    type(stl_wall_t), save :: wall

contains

    subroutine main
        use params, only: read_config, netcdffile, ns_s, ns_tp, multharm, &
                          integmode, params_init, swcoll, generate_start_only, &
                          isw_field_type, field_input, startmode, &
                          ntestpart, ntimstep, coord_input, restart
        use timing, only: init_timer, print_phase_time
        use magfie_sub, only: TEST, VMEC, SPECTRE, init_magfie
        use samplers, only: init_starting_surf, sample_spectre_surface
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

            call sample_spectre_surface(zstart)
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
        call progress_init(checkpoint_interval, ntestpart, write_results)
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
        use util, only: twopi
        use reference_coordinates, only: init_reference_coordinates, ref_coords
        use params, only: coord_input, field_input, wall_input, wall_units

        character(*), intent(in) :: vmec_file
        type(tracer_t), intent(inout) :: self
        integer, intent(in) :: ans_s, ans_tp, amultharm, aintegmode
        class(magnetic_field_t), allocatable :: field_temp
        character(:), allocatable :: vmec_equilibrium_file
        logical :: use_boozer_chartmap, use_spectre

        self%integmode = aintegmode

        ! Check if field_input is a Boozer chartmap or SPECTRE file (no VMEC needed)
        use_boozer_chartmap = .false.
        use_spectre = .false.
        if (len_trim(field_input) > 0) then
            use_spectre = is_spectre_file(field_input)
            if (.not. use_spectre) use_boozer_chartmap = is_boozer_chartmap(field_input)
        end if

        ! TEST field is analytic - no VMEC or field files needed
        if (isw_field_type == TEST) then
            self%fper = twopi  ! Full torus for analytic tokamak
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
            call init_field_can(isw_field_type)
            call print_phase_time('Canonical field initialization completed')
        else if (isw_field_type == CANFLUX .or. isw_field_type == BOOZER .or. &
                 isw_field_type == MEISS .or. isw_field_type == ALBERT) then
            call init_field_can(isw_field_type, field_temp)
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
                          spectre_ncon_phi
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
                                               spectre_ncon_phi)
            call init_field_can(SPECTRE, sf)
            call print_phase_time('SPECTRE per-volume canonical construction completed')
        end if
    end subroutine init_spectre_field

    subroutine init_stl_wall_if_enabled(coord_file)
        character(len=*), intent(in) :: coord_file

        type(chartmap_metadata_t) :: meta
        character(len=16) :: units

        wall_enabled = .false.
        wall_rho_lcfs = -1.0d0
        chartmap_cart_scale_to_m = -1.0d0

        if (len_trim(wall_input) == 0) then
            return
        end if

        if (.not. allocated(ref_coords)) then
            error stop "wall_input set but reference coordinates not initialized"
        end if

        select type (ref_coords)
        class is (chartmap_coordinate_system_t)
            continue
        class default
            error stop "wall_input requires chartmap reference coordinates"
        end select

        call read_chartmap_metadata(coord_file, meta)
        wall_rho_lcfs = meta%rho_lcfs
        chartmap_cart_scale_to_m = meta%cart_scale_to_m

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

!$omp do
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
        use util, only: twopi
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

        ! initialize lost times when particles get lost
        times_lost = -1.d0
    end subroutine init_counters

    subroutine trace_orbit(anorb, ipart, orbit_traj, orbit_times)
        use classification, only: trace_orbit_with_classifiers, &
                                  classification_result_t, &
                                  write_classification_results
        use magfie_sub, only: SPECTRE
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

        type(tracer_t), intent(inout) :: anorb
        integer, intent(in) :: ipart
        real(dp), intent(out) :: orbit_traj(:, :)  ! (5, ntimstep)
        real(dp), intent(out) :: orbit_times(:)   ! (ntimstep)

        real(dp), dimension(5) :: z
        real(dp) :: u_ref_prev(3), u_ref_cur(3), u_ref_hit(3)
        real(dp) :: x_prev(3), x_cur(3)
        real(dp) :: x_prev_m(3), x_cur_m(3), x_hit_m(3), x_hit(3)
        real(dp) :: normal_m(3), vhat(3), vnorm, cos_inc
        real(dp) :: segment_length, hit_distance, t_frac
        integer :: it, ierr_orbit, it_final, it_f
        integer(8) :: kt
        logical :: passing, faulted
        type(classification_result_t) :: class_result

        ierr_orbit = 0
        faulted = .false.
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
                        ntau_macro(it), ipart, x_prev_m)
                else
                    call macrostep(anorb, z, kt, ierr_orbit, ntau_macro(it))
                end if
            end if

            if (ierr_orbit .ne. 0) then
                it_final = it
                if (orbit_model == ORBIT_FULL_ORBIT .and. ierr_orbit == 3) then
                    ! Last-resort full-orbit fallback: the Cartesian inversion could
                    ! not resolve the position (near-axis below chartmap resolution,
                    ! or a field-period seam) and orbit_timestep_fo already warned.
                    ! This is NOT a physical loss -- the state is at the last resolved
                    ! position. Count the marker confined for the rest of the trace
                    ! so an unresolved step never registers as a lost particle.
                    do it_f = it, ntimstep
                        call increase_confined_count(it_f, passing)
                    end do
                    faulted = .true.
                end if
                exit
            end if

            ! Store trajectory data (after macrostep so time is correct)
            orbit_traj(:, it) = z
            orbit_times(it) = kt*dtaumin/v0

            call increase_confined_count(it, passing)
            it_final = it
        end do

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
        times_lost(ipart) = kt*dtaumin/v0
        if (faulted) times_lost(ipart) = trace_time   ! unresolved: confined, not lost
!$omp end critical
    end subroutine trace_orbit

    subroutine trace_orbit_spectre(ipart, z, passing, orbit_traj, orbit_times)
        !> Per-volume RK45 guiding-center trace for SPECTRE (integmode=0). The
        !> marker traverses volumes through the Level-0 crossing map and stays
        !> confined for the full trace, reflects at forbidden interfaces, or is
        !> lost at the outermost interface. Every crossing/reflection is logged.
        use spectre_orbit, only: spectre_orbit_state_t, spectre_event_t, &
                                 spectre_state_reset, orbit_timestep_spectre, &
                                 SPECTRE_OK, SPECTRE_BOUNDARY
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
        integer :: it, ktau, ierr_orbit, it_final
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
                    if (ierr_orbit /= SPECTRE_OK) exit
                    kt = kt + 1
                end do
            end if

            if (ierr_orbit /= SPECTRE_OK) then
                it_final = it
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
        else
            t_stop = real(kt, dp)*dtaumin/v0
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
                                       SYMPL_SPECTRE_OK, SYMPL_SPECTRE_SKIM
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
        integer :: it, ktau, ierr_orbit, it_final
        integer(8) :: kt
        real(dp) :: t_stop, t_frac

        call sympl_spectre_reset(state, anorb%si, spectre_mvol, integmode, &
                                 crossing_level)

        kt = 0
        it_final = 0
        ierr_orbit = SYMPL_SPECTRE_OK
        t_frac = 1.0d0

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
                exit
            end if

            call to_standard_z_coordinates(anorb, z)
            orbit_traj(:, it) = z
            orbit_times(it) = kt*dtaumin/v0
            call increase_confined_count(it, passing)
            it_final = it
        end do

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
        else if (ierr_orbit == SYMPL_SPECTRE_SKIM) then
            ! Mirror-confined at an interior interface: cannot be lost, so record
            ! it as confined (times_lost = trace_time) rather than at its stop.
            t_stop = trace_time
        else
            t_stop = (real(kt, dp) + t_frac)*dtaumin/v0
        end if

        call to_standard_z_coordinates(anorb, z)
        ! The next marker on this thread starts with an unlocked field dispatch.
        call set_spectre_volume_lock(0)

!$omp critical
        call integ_to_ref(z(1:3), zend(1:3, ipart))
        zend(4:5, ipart) = z(4:5)
        times_lost(ipart) = t_stop
!$omp end critical
    end subroutine trace_orbit_spectre_sympl

	    subroutine macrostep(anorb, z, kt, ierr_orbit, ntau_local)
        use alpha_lifetime_sub, only: orbit_timestep_axis
        use orbit_symplectic, only: orbit_timestep_sympl

        type(tracer_t), intent(inout) :: anorb
        real(dp), intent(inout) :: z(5)
        integer(8), intent(inout) :: kt
        integer, intent(out) :: ierr_orbit
        integer, intent(in) :: ntau_local

        integer :: ktau

        do ktau = 1, ntau_local
            if (orbit_model == ORBIT_FULL_ORBIT) then
                call orbit_timestep_fo(anorb%fo, z, ierr_orbit)
                if (ierr_orbit .ne. 0) exit
                kt = kt + 1
                cycle
            end if
            if (integmode <= 0) then
                call orbit_timestep_axis(z, dtaumin, dtaumin, relerr, ierr_orbit)
            else
                if (swcoll) call update_momentum(anorb, z)
                call orbit_timestep_sympl(anorb%si, anorb%f, ierr_orbit)
                call to_standard_z_coordinates(anorb, z)
            end if
            if (swcoll) call collide(z, dtaumin) ! Collisions
            if (ierr_orbit .ne. 0) exit
            kt = kt + 1
        end do
    end subroutine macrostep

    subroutine macrostep_with_wall_check(anorb, z, kt, ierr_orbit, ntau_local, &
            ipart, x_prev_m)
        use alpha_lifetime_sub, only: orbit_timestep_axis
        use orbit_symplectic, only: orbit_timestep_sympl

        type(tracer_t), intent(inout) :: anorb
        real(dp), intent(inout) :: z(5)
        integer(8), intent(inout) :: kt
        integer, intent(out) :: ierr_orbit
        integer, intent(in) :: ntau_local
        integer, intent(in) :: ipart
        real(dp), intent(inout) :: x_prev_m(3)

        integer :: ktau, wall_check_interval
        real(dp) :: u_ref_prev(3), u_ref_cur(3), x_cur(3), x_cur_m(3)
        real(dp) :: x_hit_m(3), x_hit(3), normal_m(3)
        real(dp) :: vhat(3), vnorm, cos_inc
        real(dp) :: u_ref_hit(3)
        real(dp) :: segment_length, hit_distance, t_frac
        logical :: hit
        integer :: ierr_from_cart

        call integ_to_ref(z(1:3), u_ref_prev)

        ! Check wall every N microsteps to limit overhead
        ! For small macrosteps (ntau_local<=32), check at end only
        ! For larger macrosteps, check every 32 microsteps
        wall_check_interval = max(1, min(32, ntau_local))

        do ktau = 1, ntau_local
            if (integmode <= 0) then
                call orbit_timestep_axis(z, dtaumin, dtaumin, relerr, ierr_orbit)
            else
                if (swcoll) call update_momentum(anorb, z)
                call orbit_timestep_sympl(anorb%si, anorb%f, ierr_orbit)
                call to_standard_z_coordinates(anorb, z)
            end if
            if (swcoll) call collide(z, dtaumin) ! Collisions
            if (ierr_orbit .ne. 0) exit
            kt = kt + 1

            ! Check wall intersection periodically or at end of macrostep
            if (mod(ktau, wall_check_interval) == 0 .or. ktau == ntau_local) then
                call integ_to_ref(z(1:3), u_ref_cur)
                call ref_coords%evaluate_cart(u_ref_cur, x_cur)
                x_cur_m = x_cur*chartmap_cart_scale_to_m

                if (u_ref_cur(1) > wall_rho_lcfs) then
                    call stl_wall_first_hit_segment_with_normal( &
                        wall, x_prev_m, x_cur_m, hit, x_hit_m, normal_m)
                    if (hit) then
                        wall_hit(ipart) = 1_int8
                        x_hit = x_hit_m/chartmap_cart_scale_to_m
                        wall_hit_cart(:, ipart) = x_hit
                        wall_hit_normal_cart(:, ipart) = normal_m

                        vhat = x_cur_m - x_prev_m
                        vnorm = sqrt(sum(vhat*vhat))
                        if (vnorm > 0.0_dp) then
                            vhat = vhat/vnorm
                            cos_inc = abs(sum(vhat*normal_m))
                            cos_inc = min(1.0_dp, max(0.0_dp, cos_inc))
                            wall_hit_cos_incidence(ipart) = cos_inc
                            wall_hit_angle_rad(ipart) = acos(cos_inc)
                        end if

                        ierr_from_cart = 0
                        select type (ccs => ref_coords)
                        class is (chartmap_coordinate_system_t)
                            call ccs%from_cart(x_hit, u_ref_hit, ierr_from_cart)
                        class default
                            ierr_from_cart = 1
                        end select

                        if (ierr_from_cart == 0) then
                            call ref_to_integ(u_ref_hit, z(1:3))
                        else
                            ! Fallback: linear interpolation of reference coordinates
                            ! when from_cart fails (ill-conditioned regions of chartmap)
                            segment_length = sqrt(sum((x_cur_m - x_prev_m)**2))
                            if (segment_length > 0.0_dp) then
                                hit_distance = sqrt(sum((x_hit_m - x_prev_m)**2))
                                t_frac = hit_distance / segment_length
                                u_ref_hit = u_ref_prev + t_frac * (u_ref_cur - u_ref_prev)
                                call ref_to_integ(u_ref_hit, z(1:3))
                            end if
                        end if

                        ierr_orbit = 77
                        exit
                    end if
                end if

                x_prev_m = x_cur_m
                u_ref_prev = u_ref_cur
            end if
        end do
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
        !> shared result arrays. progress_monitor calls this every
        !> checkpoint_interval seconds, so a run killed mid-flight keeps its
        !> last flushed output. Confined fractions are normalised by ntestpart,
        !> as in the final write, so a partial file is a converging lower bound
        !> and never exceeds one; particles not yet finished keep their sentinel
        !> values (times_lost = -1).
        integer :: i, num_lost, unit
        real(dp) :: inverse_times_lost_sum, norm

        norm = real(max(ntestpart, 1), dp)

        open (newunit=unit, file='times_lost.dat', recl=1024)
        num_lost = 0
        inverse_times_lost_sum = 0.0d0
        do i = 1, ntestpart
            write (unit, *) i, times_lost(i), trap_par(i), zstart(1, i), &
                perp_inv(i), zend(:, i)
            if (times_lost(i) > 0.0d0 .and. times_lost(i) < trace_time) then
                num_lost = num_lost + 1
                inverse_times_lost_sum = inverse_times_lost_sum + 1.0/times_lost(i)
            end if
        end do
        close (unit)

        if (num_lost > 0) then
            ! Write average loss time.
            open (newunit=unit, file='avg_inverse_t_lost.dat', recl=1024)
            write (unit, *) inverse_times_lost_sum/num_lost
            close (unit)
        end if

        open (newunit=unit, file='confined_fraction.dat', recl=1024)
        do i = 1, ntimstep
            write (unit, *) dble(kt_macro(i))*dtaumin/v0, confpart_pass(i)/norm, &
                confpart_trap(i)/norm, ntestpart
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
        use spectre_sympl_orbit, only: sympl_landing_stats
        use magfie_sub, only: SPECTRE

        integer :: landings, stops
        real(dp) :: max_resid

        if (isw_field_type /= SPECTRE) return

        call crossing_log_write('spectre_crossing_events.dat')

        if (integmode > 0) then
            call sympl_landing_stats(landings, max_resid, stops)
            print '(A,I0,A,ES12.4,A,I0)', 'sympl_landing_stats: count= ', &
                landings, ' max_resid= ', max_resid, ' cross_stop= ', stops
        end if
    end subroutine write_spectre_crossing_events

end module simple_main
