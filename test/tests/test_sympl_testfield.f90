module failed_symplectic_step_backend
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_base, only: field_can_t
    use field_can_mod, only: evaluate
    use orbit_symplectic_base, only: symplectic_integrator_t, SYMPLECTIC_STEP_BOUNDARY

    implicit none

    integer :: successful_drift_calls = 0
    integer :: injected_status = 0
    real(dp) :: injected_momentum_ratio = 1d0
    logical :: injected_stale_cache = .false.

contains

    subroutine fail_symplectic_step(si, f, ierr)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        f%Bmod = 8.0d0
        f%vpar = 3.0d0
        f%mu = 1.0d0
        ierr = 1
    end subroutine fail_symplectic_step

    subroutine locate_lcfs_step(si, f, ierr)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        si%z = [1.0d0, 0.2d0, 0.3d0, 0.4d0]
        si%last_step_fraction = 0.25d0
        si%last_event_radial_residual = 1.0d-11
        si%last_event_fraction_width = 1.0d-12
        call evaluate(f, si%z(1), si%z(2), si%z(3), 0)
        f%ro0 = 1d0
        f%mu = si%pabs**2/f%Bmod
        f%vpar = 0d0
        si%z(4) = f%Aph/f%ro0
        ierr = SYMPLECTIC_STEP_BOUNDARY
    end subroutine locate_lcfs_step

    subroutine wrong_energy_step(si, f, ierr)
        use params, only: radial_electric_potential_slope
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        si%z(1) = 0.3d0
        if (injected_status == SYMPLECTIC_STEP_BOUNDARY) then
            si%z(1) = 1d0
            si%last_step_fraction = 0.25d0
        end if
        call evaluate(f, si%z(1), si%z(2), si%z(3), 0)
        f%vpar = sqrt(2d0*((injected_momentum_ratio*si%pabs)**2 - f%mu*f%Bmod))
        si%z(4) = f%vpar*f%hph + f%Aph/f%ro0
        ! Model Euler's stale cache without changing its canonical endpoint.
        if (injected_stale_cache) &
            f%vpar = sqrt(2d0*((1.1d0*si%pabs)**2 - f%mu*f%Bmod))
        if (radial_electric_potential_slope /= 0d0) &
            f%Phie = f%Phie + (1d0 - injected_momentum_ratio**2)*si%pabs**2
        f%H = 0.5d0*f%vpar**2 + f%mu*f%Bmod + f%Phie
        ierr = injected_status
    end subroutine wrong_energy_step

    subroutine successful_drift_step(si, f, ierr)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        successful_drift_calls = successful_drift_calls + 1
        si%z(1) = si%z(1) + 0.01d0
        call evaluate(f, si%z(1), si%z(2), si%z(3), 0)
        f%mu = si%pabs**2/f%Bmod
        f%vpar = 0d0
        si%z(4) = f%Aph/f%ro0
        f%H = si%pabs**2
        ierr = 0
    end subroutine successful_drift_step

end module failed_symplectic_step_backend

program test_sympl_testfield
    use, intrinsic :: iso_fortran_env, only : dp => real64, int64
    use failed_symplectic_step_backend, only: fail_symplectic_step, &
                         locate_lcfs_step, successful_drift_step, successful_drift_calls
    use failed_symplectic_step_backend, only: wrong_energy_step, injected_status, &
                                              injected_momentum_ratio, &
                                              injected_stale_cache
    use simple_main, only: update_momentum
    use classification, only: classify_classifier_exit
    use simple_main, only : classify_orbit_exit, init_field, locate_linear_lcfs, &
        locate_validated_lcfs, macrostep, macrostep_with_wall_check, &
        recovery_state_has_valid_phase_space, recovery_state_is_physical, &
        rk_recovery_state_t, validate_rk_state, &
        validate_rk_step_or_reference_boundary, &
        activate_rk_recovery, should_resume_symplectic, &
        normalized_flux_from_reference_radial, &
        should_classify_numerically_confined
    use simple, only : tracer_t, init_sympl, reseed_sympl, ORBIT_FO_LOSS, &
        ORBIT_FO_NUMERICAL
    use params, only : isw_field_type, field_input, coord_input, integmode, &
        swcoll, orbit_model, ORBIT_GC, ORBIT_FULL_ORBIT, ORBIT_EXIT_LCFS, &
        ORBIT_EXIT_WALL, ORBIT_EXIT_NUMERICAL_DOMAIN, &
        ORBIT_EXIT_NUMERICAL_MAXITER, ORBIT_EXIT_NUMERICAL_EVENT, &
        ORBIT_EXIT_NUMERICAL_FULL_ORBIT
    use magfie_sub, only : TEST, init_magfie
    use field_can_mod, only : evaluate
    use orbit_symplectic, only : orbit_timestep_sympl
    use orbit_symplectic_base, only: SYMPLECTIC_STEP_BOUNDARY, &
        SYMPLECTIC_STEP_OK, SYMPLECTIC_STEP_MAXITER, &
        symplectic_newton_warning_mode

    implicit none

    type(tracer_t) :: norb
    character(len=:), allocatable :: vmec_file
    integer, parameter :: ans_s = 5, ans_tp = 5, amultharm = 5
    integer :: ierr, path_length
    real(dp) :: initial_si_z(4)

    call get_command_argument(1, length=path_length)
    if (path_length < 1) error stop 'Expected VMEC fixture path'
    allocate (character(len=path_length) :: vmec_file)
    call get_command_argument(1, vmec_file)

    ! Configure symplectic GC with TEST field to ensure initialization succeeds
    isw_field_type = TEST
    field_input = vmec_file
    coord_input = vmec_file
    integmode = 1

    call init_field(norb, vmec_file, ans_s, ans_tp, amultharm, integmode)
    call init_magfie(TEST)

    if (.not. associated(evaluate)) then
        print *, 'evaluate pointer not associated for TEST field'
        stop 1
    end if

    ! Smoke step: one symplectic timestep with test field
    norb%relerr = 1.0e-12_dp
    norb%dtaumin = 1.0e-3_dp
    norb%dtau = 1.0e-3_dp

    call init_sympl(norb%si, norb%f, (/0.2_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/), &
        norb%dtau, norb%dtaumin, norb%relerr, 1)

    ierr = 0
    initial_si_z = norb%si%z
    call orbit_timestep_sympl(norb%si, norb%f, ierr)
    if (ierr == SYMPLECTIC_STEP_OK) then
        error stop 'unconverged Euler test-field step lost its failure status'
    end if
    if (any(norb%si%z /= initial_si_z)) then
        error stop 'unconverged Euler test-field step changed accepted position'
    end if

    call test_macrostep_lcfs_event
    call test_failed_step_preserves_state
    call test_collision_failure_propagation
    call test_wrong_energy_symplectic_endpoint
    call test_rk_recovery_regions
    call test_exit_classification
    call test_numerically_confined_classification
    call test_fo_lcfs_location
    call test_validated_rk_lcfs_location
    call test_rk_state_validation
    call test_symplectic_step_validation
    call test_recovery_state_validation
    call test_reseed_momentum_reference

    print *, 'TEST field symplectic step succeeded'

contains

    subroutine test_macrostep_lcfs_event
        real(dp) :: z(5), exit_step
        integer(int64) :: kt
        integer :: step_error

        swcoll = .false.
        orbit_model = ORBIT_GC
        orbit_timestep_sympl => locate_lcfs_step
        z = [0.9_dp, 0.1_dp, 0.2_dp, 1.0_dp, 0.0_dp]
        kt = 0_int64
        call macrostep(norb, z, kt, step_error, 1, exit_step)

        if (step_error /= SYMPLECTIC_STEP_BOUNDARY) then
            error stop 'located LCFS event lost its physical status'
        end if
        if (kt /= 0_int64) error stop 'fractional LCFS event advanced a full step'
        if (exit_step /= 0.25_dp) error stop 'LCFS event time lost its step fraction'
        if (z(1) /= 1.0_dp) error stop 'LCFS event endpoint was not committed'
    end subroutine test_macrostep_lcfs_event

    subroutine test_failed_step_preserves_state
        real(dp), parameter :: initial_state(5) = [0.2_dp, 1.0_dp, 2.0_dp, &
            1.0_dp, 0.25_dp]
        real(dp) :: z(5)
        real(dp) :: x_previous(3)
        integer(int64) :: kt
        integer :: hold_streak, step_error
        logical :: numerical_hold

        swcoll = .false.
        orbit_model = ORBIT_GC
        orbit_timestep_sympl => fail_symplectic_step
        z = initial_state
        kt = 0_int64
        symplectic_newton_warning_mode = .false.

        call macrostep(norb, z, kt, step_error, 1)

        if (step_error /= 1) error stop 'failed step status was not preserved'
        if (kt /= 0_int64) error stop 'failed step advanced the time index'
        if (any(z /= initial_state)) then
            error stop 'failed symplectic step changed the accepted state'
        end if

        z = initial_state
        x_previous = 0.0_dp
        kt = 0_int64
        call macrostep_with_wall_check(norb, z, kt, step_error, 1, 1, x_previous)
        if (step_error /= 1) error stop 'wall path lost failed step status'
        if (kt /= 0_int64) error stop 'failed wall path advanced the time index'
        if (any(z /= initial_state)) then
            error stop 'failed wall path changed the accepted state'
        end if

        z = initial_state
        kt = 0_int64
        hold_streak = 0
        symplectic_newton_warning_mode = .true.
        call macrostep(norb, z, kt, step_error, 1, hold_streak=hold_streak, &
            numerical_hold_any=numerical_hold)
        if (step_error /= 0) error stop 'isolated warning failure stopped orbit'
        if (kt /= 1_int64) error stop 'isolated warning hold lost its interval'
        if (.not. numerical_hold .or. hold_streak /= 1) &
            error stop 'isolated warning hold was not reported'
        if (any(z /= initial_state)) then
            error stop 'isolated warning hold changed the last accepted state'
        end if

        call macrostep(norb, z, kt, step_error, 1, hold_streak=hold_streak, &
            numerical_hold_any=numerical_hold)
        if (step_error == 0) error stop 'repeated warning failure was hidden'
        if (kt /= 1_int64) error stop 'repeated failure consumed another interval'
        if (numerical_hold .or. hold_streak /= 1) &
            error stop 'repeated warning failure was misreported'
        if (any(z /= initial_state)) then
            error stop 'repeated warning failure changed the last accepted state'
        end if

        z = initial_state
        x_previous = 0.0_dp
        kt = 0_int64
        hold_streak = 0
        call macrostep_with_wall_check(norb, z, kt, step_error, 1, 1, x_previous, &
            hold_streak=hold_streak, numerical_hold_any=numerical_hold)
        if (step_error /= 0) error stop 'isolated wall warning stopped orbit'
        if (kt /= 1_int64) error stop 'isolated wall hold lost its interval'
        if (.not. numerical_hold .or. hold_streak /= 1) &
            error stop 'isolated wall hold was not reported'
        if (any(z /= initial_state)) then
            error stop 'isolated wall hold changed the last accepted state'
        end if

        call macrostep_with_wall_check(norb, z, kt, step_error, 1, 1, x_previous, &
            hold_streak=hold_streak, numerical_hold_any=numerical_hold)
        if (step_error == 0) error stop 'repeated wall failure was hidden'
        if (kt /= 1_int64) error stop 'repeated wall failure consumed another interval'
        if (numerical_hold .or. hold_streak /= 1) &
            error stop 'repeated wall failure was misreported'
        if (any(z /= initial_state)) then
            error stop 'repeated wall failure changed the last accepted state'
        end if
    end subroutine test_failed_step_preserves_state

    subroutine test_collision_failure_propagation
        use collis_alp, only: efcolf_grid, velrat_grid, enrat_grid, stost
        use parmot_mod, only: ro0
        use finite_check, only: is_finite
        use reference_coordinates, only: ref_coords
        use libneo_coordinates, only: coordinate_system_t
        use cartesian_coordinates, only: cartesian_coordinate_system_t
        use field_can_mod, only: field_can_t
        use orbit_symplectic_base, only: symplectic_integrator_t, &
                                         SYMPLECTIC_STEP_EVENT_NOT_CONVERGED
        use params, only: dtaumin, wall_query_rho_lcfs
        type(tracer_t) :: local_orbit
        class(coordinate_system_t), allocatable :: saved_ref_coords
        type(symplectic_integrator_t) :: saved_si
        type(field_can_t) :: saved_field
        real(dp), parameter :: initial(5) = [0.2_dp, 1.0_dp, 0.2_dp, &
                                             1.0_dp, 0.0_dp]
        real(dp) :: z(5), x_previous(3), saved_wall_gate, saved_ro0, saved_dtaumin
        integer(int64) :: kt
        integer :: warning, wall_path, step_error, hold_streak
        logical :: numerical_hold

        ! A real stost failure, not a mock collision callback: negative
        ! diffusion coefficients must be rejected by the production operator.
        saved_ro0 = ro0
        saved_dtaumin = dtaumin
        ro0 = 1.0_dp ! TEST field uses explicit normalized finite units.
        efcolf_grid = -1.0_dp
        velrat_grid = 1.0_dp
        enrat_grid = 1.0_dp
        dtaumin = 1.0e-3_dp
        saved_wall_gate = wall_query_rho_lcfs
        ! TEST field initialization supplies no reference coordinate object.
        ! Use the production Cartesian implementation for the wall path.
        call move_alloc(ref_coords, saved_ref_coords)
        allocate (cartesian_coordinate_system_t :: ref_coords)
        wall_query_rho_lcfs = huge(1.0_dp)
        orbit_model = ORBIT_GC
        integmode = 3
        swcoll = .true.
        z = initial
        call stost(z, dtaumin, 1, step_error)
        if (step_error /= 1 .or. any(z /= initial)) &
            error stop 'Collision rejection fixture did not reject bad diffusion'
        do warning = 0, 1
            symplectic_newton_warning_mode = warning == 1
            do wall_path = 0, 1
                local_orbit = norb
                call init_sympl(local_orbit%si, local_orbit%f, initial, &
                                1.0e-3_dp, 1.0e-3_dp, 1.0e-12_dp, 3)
                if (.not. all(is_finite(local_orbit%si%z)) .or. &
                    .not. all(is_finite([local_orbit%f%H, &
                                         local_orbit%f%mu, local_orbit%f%vpar]))) &
                    error stop 'Collision caller fixture initialized nonfinite state'
                if (abs(local_orbit%f%H - initial(4)**2) > 1.0e-12_dp) &
                    error stop 'Collision caller fixture initialized wrong energy'
                saved_si = local_orbit%si
                saved_field = local_orbit%f
                orbit_timestep_sympl => successful_drift_step
                successful_drift_calls = 0
                z = initial
                x_previous = [1.0_dp, 2.0_dp, 3.0_dp]
                kt = 7_int64
                hold_streak = 0
                if (wall_path == 0) then
                    call macrostep(local_orbit, z, kt, step_error, 1, &
                             hold_streak=hold_streak, numerical_hold_any=numerical_hold)
                else
                    call macrostep_with_wall_check(local_orbit, z, kt, &
                                step_error, 1, 1, x_previous, hold_streak=hold_streak, &
                                                   numerical_hold_any=numerical_hold)
                end if
                if (successful_drift_calls /= 1) &
                    error stop 'Collision caller fixture did not complete its drift'
                if (step_error /= SYMPLECTIC_STEP_EVENT_NOT_CONVERGED) &
                    error stop 'Collision failure did not terminate microstep'
                if (classify_orbit_exit(step_error, ORBIT_GC, integmode, &
                                       ORBIT_EXIT_LCFS) /= ORBIT_EXIT_NUMERICAL_EVENT) &
                    error stop 'Collision failure was classified as physical'
                if (kt /= 7_int64 .or. any(z /= initial)) &
                    error stop 'Collision failure committed drift or elapsed time'
                if (hold_streak /= 0 .or. numerical_hold) &
                    error stop 'Collision failure was consumed as warning hold'
                if (any(local_orbit%si%z /= saved_si%z) .or. &
                    local_orbit%si%pabs /= saved_si%pabs .or. &
                    local_orbit%f%H /= saved_field%H .or. &
                    local_orbit%f%mu /= saved_field%mu .or. &
                    local_orbit%f%Bmod /= saved_field%Bmod) &
                    error stop 'Collision failure changed accepted hidden state'
                if (any(x_previous /= [1.0_dp, 2.0_dp, 3.0_dp])) &
                    error stop 'Collision failure committed wall segment origin'
            end do
        end do
        wall_query_rho_lcfs = saved_wall_gate
        deallocate (ref_coords)
        call move_alloc(saved_ref_coords, ref_coords)
        ro0 = saved_ro0
        dtaumin = saved_dtaumin
        efcolf_grid = 0.0_dp
        swcoll = .false.
    end subroutine test_collision_failure_propagation

    subroutine test_wrong_energy_symplectic_endpoint
        use orbit_symplectic_base, only: SYMPLECTIC_STEP_EVENT_NOT_CONVERGED
        use parmot_mod, only: ro0
        use reference_coordinates, only: ref_coords
        use libneo_coordinates, only: coordinate_system_t
        use cartesian_coordinates, only: cartesian_coordinate_system_t
        use params, only: wall_query_rho_lcfs, radial_electric_potential_slope
        type(tracer_t) :: local_orbit, expected_orbit, candidate_orbit
        class(coordinate_system_t), allocatable :: saved_ref_coords
        real(dp) :: z(5), initial(5), x_previous(3), saved_ro0, exit_step
        real(dp) :: saved_wall_gate, saved_potential_slope, cached_ratio
        integer(int64) :: kt
        integer :: collisions, endpoint, wall_path, step_error, excessive
        integer :: candidate_error

        call move_alloc(ref_coords, saved_ref_coords)
        allocate (cartesian_coordinate_system_t :: ref_coords)
        saved_potential_slope = radial_electric_potential_slope
        saved_wall_gate = wall_query_rho_lcfs
        wall_query_rho_lcfs = 1.0_dp
        local_orbit = norb
        saved_ro0 = ro0
        ro0 = 1.0_dp ! Normalized TEST-field charge-to-mass scale.
        symplectic_newton_warning_mode = .false.
        orbit_model = ORBIT_GC
        integmode = 3
        do excessive = 0, 3
            injected_stale_cache = excessive == 3
            radial_electric_potential_slope = 0.0_dp
            if (excessive == 2) radial_electric_potential_slope = -1.0_dp
            injected_momentum_ratio = 1.004_dp
            if (excessive == 1) injected_momentum_ratio = 1.006_dp
            if (excessive == 2) injected_momentum_ratio = 1.1_dp
            do collisions = 0, 1
                swcoll = collisions == 1
                ! A collision can change both momentum and pitch before this push.
                initial = [0.2_dp, 1.0_dp, 0.2_dp, 1.0_dp, 0.25_dp]
                do endpoint = 0, 1
                    ! Bounded boundary states return before another collision kick.
                    if (excessive /= 1 .and. excessive /= 3 .and. endpoint == 0) cycle
                    if (swcoll .and. endpoint == 0 .and. excessive == 3) cycle
                    do wall_path = 0, 1
                        initial(4:5) = [1.0_dp, 0.25_dp]
                        call init_sympl(local_orbit%si, local_orbit%f, initial, &
                                        1.0e-3_dp, 1.0e-3_dp, 1.0e-12_dp, 3)
                        ! Leave a stale pre-collision reference in the active orbit.
                        if (swcoll) initial(4:5) = [0.7_dp, -0.4_dp]
                        expected_orbit = local_orbit
                        if (swcoll) call update_momentum(expected_orbit, initial)
                        injected_status = SYMPLECTIC_STEP_OK
                        if (endpoint == 1) injected_status = SYMPLECTIC_STEP_BOUNDARY
                        orbit_timestep_sympl => wrong_energy_step
                        candidate_orbit = expected_orbit
                        call wrong_energy_step(candidate_orbit%si, &
                                               candidate_orbit%f, candidate_error)
                        z = initial
                        x_previous = 0.0_dp
                        kt = 0_int64
                        if (wall_path == 0) then
                            call macrostep(local_orbit, z, kt, step_error, 1, exit_step)
                        else
                            call macrostep_with_wall_check(local_orbit, z, kt, &
                                                step_error, 1, 1, x_previous, exit_step)
                        end if
                        if (excessive /= 1) then
                            cached_ratio = injected_momentum_ratio
                            if (injected_stale_cache) cached_ratio = 1.1_dp
                            if (step_error /= injected_status .or. &
                                abs(z(4)/initial(4) - cached_ratio) > 1.0e-12_dp) &
                                error stop 'Physical endpoint momentum was rejected'
                            if (any(local_orbit%si%z /= candidate_orbit%si%z) .or. &
                                local_orbit%f%Bmod /= candidate_orbit%f%Bmod .or. &
                                local_orbit%f%vpar /= candidate_orbit%f%vpar .or. &
                                local_orbit%f%hph /= candidate_orbit%f%hph .or. &
                                local_orbit%f%Aph /= candidate_orbit%f%Aph) &
                                error stop 'Validation mutated the carried endpoint'
                            if (excessive == 2) then
                                if (abs(local_orbit%f%H - expected_orbit%f%H) > &
                                    1.0e-12_dp) error stop 'Electric energy mismatch'
                            end if
                            if (endpoint == 1) then
                                if (kt /= 0_int64 .or. exit_step /= 0.25_dp) &
                                    error stop 'Bounded boundary time was changed'
                            else
                                if (kt /= 1_int64) &
                                    error stop 'Accepted ordinary step did not advance'
                            end if
                            cycle
                        end if
                        if (step_error /= SYMPLECTIC_STEP_EVENT_NOT_CONVERGED) &
                            error stop 'Wrong-energy endpoint accepted as a step or hit'
                        if (classify_orbit_exit(step_error, ORBIT_GC, integmode, &
                                                ORBIT_EXIT_LCFS) /= &
                            ORBIT_EXIT_NUMERICAL_EVENT) &
                        error stop 'Wrong-energy rejection classified as Newton failure'
                        if (kt /= 0_int64 .or. exit_step /= 0.0_dp .or. &
                            any(z /= initial) .or. any(x_previous /= 0.0_dp)) &
                            error stop 'Wrong-energy rejection changed particle or time'
                        if (any(local_orbit%si%z /= expected_orbit%si%z) .or. &
                            local_orbit%si%pabs /= expected_orbit%si%pabs .or. &
                            local_orbit%f%vpar /= expected_orbit%f%vpar .or. &
                            local_orbit%f%mu /= expected_orbit%f%mu .or. &
                            local_orbit%f%H /= expected_orbit%f%H) &
                        error stop 'Wrong-energy rejection did not restore hidden state'
                    end do
                end do
            end do
        end do
        injected_stale_cache = .false.
        swcoll = .false.
        ro0 = saved_ro0
        radial_electric_potential_slope = saved_potential_slope
        wall_query_rho_lcfs = saved_wall_gate
        deallocate (ref_coords)
        call move_alloc(saved_ref_coords, ref_coords)
    end subroutine test_wrong_energy_symplectic_endpoint

    subroutine test_rk_recovery_regions
        type(rk_recovery_state_t) :: recovery
        real(dp) :: z(5)

        z = [0.005_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
        call activate_rk_recovery(recovery, z, 1.0_dp)
        if (.not. recovery%active) error stop 'near-axis RK recovery did not start'
        if (.not. recovery%has_momentum_reference) then
            error stop 'RK recovery did not retain its momentum reference'
        end if
        if (recovery%momentum_reference /= 1.0_dp) then
            error stop 'RK recovery changed its momentum reference'
        end if
        recovery%safe_steps = 255
        z(1) = 0.12_dp
        if (should_resume_symplectic(recovery, z)) then
            error stop 'near-axis recovery resumed inside its hysteresis band'
        end if
        recovery%safe_steps = 256
        if (.not. should_resume_symplectic(recovery, z)) then
            error stop 'near-axis recovery did not resume after a safe streak'
        end if

        recovery = rk_recovery_state_t()
        z(1) = 0.5_dp
        call activate_rk_recovery(recovery, z, 1.0_dp)
        recovery%safe_steps = 255
        if (should_resume_symplectic(recovery, z)) then
            error stop 'generic recovery ignored its first cooldown'
        end if
        recovery%safe_steps = 256
        if (.not. should_resume_symplectic(recovery, z)) then
            error stop 'generic recovery did not probe after its cooldown'
        end if
        call activate_rk_recovery(recovery, z, 1.0_dp)
        recovery%safe_steps = 511
        if (should_resume_symplectic(recovery, z)) then
            error stop 'repeated generic recovery did not back off'
        end if
        recovery%safe_steps = 512
        if (.not. should_resume_symplectic(recovery, z)) then
            error stop 'backed-off generic recovery never resumed probing'
        end if

        recovery = rk_recovery_state_t()
        z(1) = 0.995_dp
        call activate_rk_recovery(recovery, z, 1.0_dp)
        recovery%safe_steps = 256
        z(1) = 0.99_dp
        if (should_resume_symplectic(recovery, z)) then
            error stop 'edge recovery resumed before moving safely inward'
        end if
        z(1) = 0.98_dp
        if (.not. should_resume_symplectic(recovery, z)) then
            error stop 'edge recovery did not resume after moving inward'
        end if
    end subroutine test_rk_recovery_regions

    subroutine test_reseed_momentum_reference
        real(dp) :: z(5)

        z = [0.2_dp, 1.0_dp, 0.0_dp, 1.02_dp, 0.25_dp]
        call reseed_sympl(norb%si, norb%f, z, 1.0_dp)
        if (norb%si%pabs /= 1.0_dp) then
            error stop 'reseed ratcheted the symplectic momentum reference'
        end if
        if (abs(norb%f%H - 1.0_dp) > 100.0_dp*epsilon(1.0_dp)) then
            error stop 'reseed did not project onto the reference energy shell'
        end if
        if (abs(norb%f%vpar/sqrt(2.0_dp) - 0.25_dp) > &
                100.0_dp*epsilon(1.0_dp)) then
            error stop 'momentum-shell projection changed the pitch coordinate'
        end if
    end subroutine test_reseed_momentum_reference

    subroutine test_exit_classification
        logical :: classifier_lost
        integer :: classifier_exit

        if (classify_orbit_exit(SYMPLECTIC_STEP_BOUNDARY, ORBIT_GC, 3, &
            ORBIT_EXIT_LCFS) /= &
            ORBIT_EXIT_LCFS) then
            error stop 'converged LCFS event was not classified as physical'
        end if
        if (classify_orbit_exit(SYMPLECTIC_STEP_BOUNDARY, ORBIT_GC, 3, &
            ORBIT_EXIT_NUMERICAL_DOMAIN) /= &
            ORBIT_EXIT_NUMERICAL_DOMAIN) then
            error stop 'extended map boundary was classified as physical LCFS'
        end if
        if (classify_orbit_exit(SYMPLECTIC_STEP_BOUNDARY, ORBIT_GC, 3, &
            ORBIT_EXIT_WALL) /= ORBIT_EXIT_WALL) then
            error stop 'coordinate wall was not classified as physical'
        end if
        if (classify_orbit_exit(77, ORBIT_GC, 3, ORBIT_EXIT_LCFS) /= &
            ORBIT_EXIT_WALL) then
            error stop 'wall event was not classified as physical'
        end if
        if (classify_orbit_exit(77, ORBIT_GC, 0, ORBIT_EXIT_LCFS) /= &
            ORBIT_EXIT_WALL) then
            error stop 'axis wall event was not classified as physical'
        end if
        if (classify_orbit_exit(1, ORBIT_GC, 3, ORBIT_EXIT_LCFS) /= &
            ORBIT_EXIT_NUMERICAL_DOMAIN) then
            error stop 'exterior Newton iterate was classified as physical'
        end if
        if (classify_orbit_exit(ORBIT_FO_LOSS, ORBIT_FULL_ORBIT, 3, &
            ORBIT_EXIT_LCFS) /= &
            ORBIT_EXIT_LCFS) then
            error stop 'full-orbit LCFS exit was not classified as physical'
        end if
        if (classify_orbit_exit(ORBIT_FO_LOSS, ORBIT_FULL_ORBIT, 3, &
            ORBIT_EXIT_NUMERICAL_DOMAIN) /= &
            ORBIT_EXIT_NUMERICAL_DOMAIN) then
            error stop 'extended full-orbit map boundary was classified as physical'
        end if
        if (classify_orbit_exit(ORBIT_FO_NUMERICAL, ORBIT_FULL_ORBIT, 3, &
            ORBIT_EXIT_LCFS) /= &
            ORBIT_EXIT_NUMERICAL_FULL_ORBIT) then
            error stop 'full-orbit locate failure was classified as physical'
        end if
        if (classify_orbit_exit(1, ORBIT_GC, 0, ORBIT_EXIT_LCFS) /= &
            ORBIT_EXIT_LCFS) then
            error stop 'RK LCFS exit was classified as numerical'
        end if
        if (classify_orbit_exit(1, ORBIT_GC, 0, &
            ORBIT_EXIT_NUMERICAL_DOMAIN) /= &
            ORBIT_EXIT_NUMERICAL_DOMAIN) then
            error stop 'RK extended-map boundary was classified as physical'
        end if
        if (classify_orbit_exit(2, ORBIT_GC, 0, ORBIT_EXIT_LCFS) /= &
            ORBIT_EXIT_NUMERICAL_EVENT) then
            error stop 'RK integration fault was classified as physical'
        end if

        call classify_classifier_exit(SYMPLECTIC_STEP_MAXITER, 3, &
            classifier_lost, classifier_exit)
        if (classifier_lost .or. classifier_exit /= ORBIT_EXIT_NUMERICAL_MAXITER) then
            error stop 'classifier mode promoted Newton maxiter to a physical loss'
        end if
        call classify_classifier_exit(SYMPLECTIC_STEP_BOUNDARY, 3, &
            classifier_lost, classifier_exit)
        if (.not. classifier_lost .or. classifier_exit /= ORBIT_EXIT_LCFS) then
            error stop 'classifier mode lost its physical boundary classification'
        end if
        call classify_classifier_exit(2, 0, classifier_lost, classifier_exit)
        if (classifier_lost .or. classifier_exit /= ORBIT_EXIT_NUMERICAL_EVENT) then
            error stop 'classifier mode promoted an RK fault to a physical loss'
        end if
    end subroutine test_exit_classification

    subroutine test_numerically_confined_classification
        if (abs(normalized_flux_from_reference_radial(0.1_dp, .true.) - &
                0.01_dp) > epsilon(1.0_dp)) then
            error stop 'chartmap rho was not converted to normalized flux s'
        end if
        if (normalized_flux_from_reference_radial(0.1_dp, .false.) /= 0.1_dp) then
            error stop 'VMEC normalized flux s was incorrectly rescaled'
        end if
        if (.not. should_classify_numerically_confined( &
                ORBIT_EXIT_NUMERICAL_MAXITER, ORBIT_GC, 1, .true., .false., &
                .true., 0.0099_dp)) then
            error stop 'eligible core recovery failure was not classified confined'
        end if
        if (should_classify_numerically_confined(ORBIT_EXIT_NUMERICAL_MAXITER, &
                ORBIT_GC, 1, .true., .false., .true., 0.01_dp)) then
            error stop 'axis fallback included its outer radial boundary'
        end if
        if (should_classify_numerically_confined(ORBIT_EXIT_NUMERICAL_MAXITER, &
                ORBIT_GC, 1, .true., .true., .true., 0.01_dp)) then
            error stop 'collisional marker was classified numerically confined'
        end if
        if (should_classify_numerically_confined(ORBIT_EXIT_NUMERICAL_MAXITER, &
                ORBIT_GC, 1, .false., .false., .true., 0.01_dp)) then
            error stop 'strict-mode marker was classified numerically confined'
        end if
        if (should_classify_numerically_confined(ORBIT_EXIT_NUMERICAL_MAXITER, &
                ORBIT_GC, 1, .true., .false., .false., 0.01_dp)) then
            error stop 'failure without exhausted recovery was classified confined'
        end if
        if (should_classify_numerically_confined(ORBIT_EXIT_NUMERICAL_MAXITER, &
                ORBIT_FULL_ORBIT, 1, .true., .false., .true., 0.01_dp)) then
            error stop 'full-orbit marker was classified numerically confined'
        end if
        if (should_classify_numerically_confined(ORBIT_EXIT_NUMERICAL_MAXITER, &
                ORBIT_GC, 0, .true., .false., .true., 0.01_dp)) then
            error stop 'non-symplectic marker was classified numerically confined'
        end if
        if (should_classify_numerically_confined(ORBIT_EXIT_NUMERICAL_DOMAIN, &
                ORBIT_GC, 1, .true., .false., .true., 0.01_dp)) then
            error stop 'domain failure was classified numerically confined'
        end if
    end subroutine test_numerically_confined_classification

    subroutine test_fo_lcfs_location
        real(dp), parameter :: z_before(5) = [0.9_dp, 6.2_dp, 3.0_dp, 1.0_dp, &
            -0.5_dp]
        real(dp), parameter :: z_after(5) = [1.1_dp, 0.1_dp, -3.0_dp, 2.0_dp, &
            0.5_dp]
        real(dp), parameter :: tolerance = 32.0_dp*epsilon(1.0_dp)
        real(dp) :: z_event(5), event_fraction

        call locate_linear_lcfs(z_before, z_after, 6.0_dp, z_event, event_fraction)
        if (abs(event_fraction - 0.5_dp) > tolerance) then
            error stop 'full-orbit LCFS fraction is incorrect'
        end if
        if (z_event(1) /= 1.0_dp) error stop 'full-orbit event is not on the LCFS'
        if (abs(z_event(4) - 1.5_dp) > tolerance .or. &
            abs(z_event(5)) > tolerance) then
            error stop 'full-orbit phase-space state is inconsistent with event time'
        end if
        if (abs(z_event(2) - (6.2_dp + 0.5_dp*(0.1_dp - 6.2_dp + &
            2.0_dp*acos(-1.0_dp)))) > tolerance) then
            error stop 'full-orbit poloidal seam was not interpolated periodically'
        end if
        if (abs(z_event(3) - 3.0_dp) > tolerance) then
            error stop 'full-orbit field-period seam was not interpolated periodically'
        end if
    end subroutine test_fo_lcfs_location

    subroutine test_validated_rk_lcfs_location
        real(dp), parameter :: z_before(5) = [0.9_dp, 0.2_dp, 0.3_dp, 1.0_dp, &
            -0.5_dp]
        real(dp), parameter :: z_after(5) = [1.1_dp, 0.4_dp, 0.5_dp, 1.004_dp, &
            0.5_dp]
        real(dp), parameter :: z_no_crossing(5) = [0.8_dp, 0.4_dp, 0.5_dp, &
            2.0_dp, 0.5_dp]
        real(dp), parameter :: z_invalid_event(5) = [1.1_dp, 0.4_dp, 0.5_dp, &
            1.0_dp, 4.0_dp]
        real(dp), parameter :: tolerance = 32.0_dp*epsilon(1.0_dp)
        real(dp) :: z_event(5), event_fraction
        integer :: step_error

        step_error = 1
        call locate_validated_lcfs(z_before, z_after, 6.0_dp, z_event, &
            event_fraction, step_error)
        if (step_error /= 1 .or. abs(event_fraction - 0.5_dp) > tolerance .or. &
            abs(z_event(1) - 1.0_dp) > tolerance) then
            error stop 'valid adaptive-RK LCFS crossing was rejected'
        end if

        step_error = 1
        call locate_validated_lcfs(z_before, z_no_crossing, 6.0_dp, z_event, &
            event_fraction, step_error)
        if (step_error /= 2 .or. event_fraction /= 0.0_dp .or. &
            any(z_event /= z_before)) then
            error stop 'adaptive-RK chamber flag without a crossing became a loss'
        end if

        step_error = 1
        call locate_validated_lcfs(z_before, z_invalid_event, 6.0_dp, z_event, &
            event_fraction, step_error)
        if (step_error /= 2 .or. event_fraction /= 0.0_dp .or. &
            any(z_event /= z_before)) then
            error stop 'unphysical adaptive-RK chamber event became a loss'
        end if
    end subroutine test_validated_rk_lcfs_location

    subroutine test_rk_state_validation
        real(dp), parameter :: z_before(5) = [0.2_dp, 0.2_dp, 0.3_dp, 1.0_dp, &
            -0.5_dp]
        real(dp) :: z_after(5)
        integer :: step_error

        z_after = [0.21_dp, 0.4_dp, 0.5_dp, 1.1_dp, 0.5_dp]
        step_error = 0
        call validate_rk_state(z_before, z_after, step_error)
        if (step_error /= 0) error stop 'valid adaptive-RK state was rejected'

        z_after = [0.21_dp, 0.4_dp, 0.5_dp, 1.1_dp, 0.5_dp]
        step_error = 0
        call validate_rk_state(z_before, z_after, step_error, 1.0_dp)
        if (step_error /= 2 .or. any(z_after /= z_before)) then
            error stop 'wrong-energy adaptive-RK branch was committed'
        endif

        z_after = [0.21_dp, 0.4_dp, 0.5_dp, 1.004_dp, 0.5_dp]
        step_error = 0
        call validate_rk_state(z_before, z_after, step_error, 1.0_dp)
        if (step_error /= 0) error stop 'bounded adaptive-RK momentum drift was rejected'

        z_after = [0.21_dp, 0.4_dp, 0.5_dp, 1.006_dp, 0.5_dp]
        step_error = 0
        call validate_rk_state(z_before, z_after, step_error, 1.0_dp)
        if (step_error /= 2 .or. any(z_after /= z_before)) then
            error stop 'excessive adaptive-RK momentum drift was committed'
        end if

        z_after = [0.21_dp, 0.4_dp, 0.5_dp, 1.1_dp, 2.0_dp]
        step_error = 0
        call validate_rk_state(z_before, z_after, step_error)
        if (step_error /= 2 .or. any(z_after /= z_before)) then
            error stop 'invalid adaptive-RK state was committed'
        end if
    end subroutine test_rk_state_validation

    subroutine test_symplectic_step_validation
        real(dp), parameter :: z_before(5) = [0.2_dp, 0.2_dp, 0.3_dp, 1.0_dp, &
            -0.5_dp]
        real(dp) :: z_after(5), z(5), event_fraction
        integer :: step_error
        logical :: reference_boundary

        ! Symplectic endpoints allow bounded truncation error, but share the
        ! recovery momentum sanity bound before committing a physical state.
        z_after = [0.21_dp, 0.4_dp, 0.5_dp, 1.004_dp, -0.49_dp]
        step_error = 0
        call validate_rk_step_or_reference_boundary(z_before, z_after, &
                      1.0_dp, 1.0_dp, z, event_fraction, step_error, reference_boundary)
        if (step_error /= 0 .or. any(z /= z_after)) &
            error stop 'bounded symplectic truncation error was rejected'

        z_after(4) = 1.015_dp
        step_error = 0
        call validate_rk_step_or_reference_boundary(z_before, z_after, &
            1.0_dp, 1.0_dp, z, event_fraction, step_error, reference_boundary)
        if (step_error /= 2) &
            error stop 'excessive symplectic momentum drift was accepted'

        z_after(5) = 2.0_dp
        step_error = 0
        call validate_rk_step_or_reference_boundary(z_before, z_after, &
            1.0_dp, z=z, event_fraction=event_fraction, ierr=step_error, &
            reference_boundary=reference_boundary)
        if (step_error /= 2 .or. any(z /= z_before)) &
            error stop 'invalid symplectic phase space was committed'
    end subroutine test_symplectic_step_validation

    subroutine test_recovery_state_validation
        real(dp) :: z(5)

        z = [0.2_dp, 0.3_dp, 0.4_dp, 1.0_dp, 0.5_dp]
        if (.not. recovery_state_has_valid_phase_space(z)) then
            error stop 'physical recovery state was rejected'
        end if
        if (.not. recovery_state_is_physical(z, 1.0_dp)) then
            error stop 'recovery state on the momentum shell was rejected'
        end if

        z(5) = 2.0_dp
        if (recovery_state_has_valid_phase_space(z)) then
            error stop 'recovery state with impossible pitch was accepted'
        end if

        z = [0.2_dp, 0.3_dp, 0.4_dp, 1.004_dp, 0.5_dp]
        if (.not. recovery_state_is_physical(z, 1.0_dp)) then
            error stop 'bounded recovery momentum drift was rejected'
        end if
        z(4) = 1.006_dp
        if (recovery_state_is_physical(z, 1.0_dp)) then
            error stop 'excessive recovery momentum drift was accepted'
        end if
    end subroutine test_recovery_state_validation

end program test_sympl_testfield
