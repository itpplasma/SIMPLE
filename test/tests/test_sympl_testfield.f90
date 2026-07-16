module failed_symplectic_step_backend
    use field_can_base, only: field_can_t
    use orbit_symplectic_base, only: symplectic_integrator_t, SYMPLECTIC_STEP_BOUNDARY

    implicit none

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
        f%Bmod = 1.0d0
        f%mu = 1.0d0
        f%vpar = 0.0d0
        ierr = SYMPLECTIC_STEP_BOUNDARY
    end subroutine locate_lcfs_step

end module failed_symplectic_step_backend

program test_sympl_testfield
    use, intrinsic :: iso_fortran_env, only : dp => real64, int64
    use failed_symplectic_step_backend, only: fail_symplectic_step, locate_lcfs_step
    use classification, only: classify_classifier_exit
    use simple_main, only : classify_orbit_exit, init_field, locate_linear_lcfs, &
        macrostep, macrostep_with_wall_check, rk_recovery_state_t, &
        activate_rk_recovery, should_resume_symplectic
    use simple, only : tracer_t, init_sympl, ORBIT_FO_LOSS, ORBIT_FO_NUMERICAL
    use params, only : isw_field_type, field_input, coord_input, integmode, &
        swcoll, orbit_model, ORBIT_GC, ORBIT_FULL_ORBIT, ORBIT_EXIT_LCFS, &
        ORBIT_EXIT_WALL, ORBIT_EXIT_NUMERICAL_DOMAIN, &
        ORBIT_EXIT_NUMERICAL_MAXITER, ORBIT_EXIT_NUMERICAL_EVENT, &
        ORBIT_EXIT_NUMERICAL_FULL_ORBIT
    use magfie_sub, only : TEST
    use field_can_mod, only : evaluate
    use orbit_symplectic, only : orbit_timestep_sympl
    use orbit_symplectic_base, only: SYMPLECTIC_STEP_BOUNDARY, &
        SYMPLECTIC_STEP_OK, SYMPLECTIC_STEP_MAXITER, &
        symplectic_newton_warning_mode

    implicit none

    type(tracer_t) :: norb
    character(*), parameter :: vmec_file = WOUT_FILE
    integer, parameter :: ans_s = 5, ans_tp = 5, amultharm = 5
    integer :: ierr
    real(dp) :: initial_si_z(4)

    ! Configure symplectic GC with TEST field to ensure initialization succeeds
    isw_field_type = TEST
    field_input = vmec_file
    coord_input = vmec_file
    integmode = 1

    call init_field(norb, vmec_file, ans_s, ans_tp, amultharm, integmode)

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
    call test_rk_recovery_regions
    call test_exit_classification
    call test_fo_lcfs_location

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
        if (step_error /= 0) error stop 'warning mode stopped a failed step'
        if (kt /= 1_int64) error stop 'warning mode did not consume the failed step'
        if (.not. numerical_hold .or. hold_streak == 0) &
            error stop 'warning mode did not report the held interval'
        if (any(z /= initial_state)) then
            error stop 'warning-mode skip changed the last accepted state'
        end if
        call macrostep(norb, z, kt, step_error, 1, hold_streak=hold_streak, &
            numerical_hold_any=numerical_hold)
        if (step_error /= 0) &
            error stop 'repeated warning failure stopped the orbit'
        if (kt /= 2_int64) &
            error stop 'repeated warning failure did not consume one interval'
        if (.not. numerical_hold .or. hold_streak == 0) &
            error stop 'repeated warning failure lost its diagnostic status'
        if (any(z /= initial_state)) &
            error stop 'repeated warning hold changed the accepted state'

        z = initial_state
        x_previous = 0.0_dp
        kt = 0_int64
        hold_streak = 0
        call macrostep_with_wall_check(norb, z, kt, step_error, 1, 1, x_previous, &
            hold_streak=hold_streak, numerical_hold_any=numerical_hold)
        if (step_error /= 0) error stop 'warning-mode wall path stopped a failed step'
        if (kt /= 1_int64) &
            error stop 'warning-mode wall path did not consume the failed step'
        if (.not. numerical_hold .or. hold_streak == 0) &
            error stop 'warning-mode wall path did not latch the failure'
        if (any(z /= initial_state)) then
            error stop 'warning-mode wall skip changed the last accepted state'
        end if
        call macrostep_with_wall_check(norb, z, kt, step_error, 1, 1, x_previous, &
            hold_streak=hold_streak, numerical_hold_any=numerical_hold)
        if (step_error /= 0) &
            error stop 'repeated wall warning failure stopped the orbit'
        if (kt /= 2_int64) &
            error stop 'repeated wall warning did not consume one interval'
        if (.not. numerical_hold .or. hold_streak == 0) &
            error stop 'repeated wall warning lost its diagnostic status'
        if (any(z /= initial_state)) &
            error stop 'repeated wall warning changed the accepted state'
    end subroutine test_failed_step_preserves_state

    subroutine test_rk_recovery_regions
        type(rk_recovery_state_t) :: recovery

        call activate_rk_recovery(recovery, 0.005_dp)
        if (.not. recovery%active) error stop 'near-axis RK recovery did not start'
        if (should_resume_symplectic(recovery, 0.015_dp)) then
            error stop 'near-axis recovery resumed inside its hysteresis band'
        end if
        if (.not. should_resume_symplectic(recovery, 0.02_dp)) then
            error stop 'near-axis recovery did not resume outside the axis band'
        end if

        recovery = rk_recovery_state_t()
        call activate_rk_recovery(recovery, 0.5_dp)
        recovery%steps = 255
        if (should_resume_symplectic(recovery, 0.5_dp)) then
            error stop 'generic recovery ignored its first cooldown'
        end if
        recovery%steps = 256
        if (.not. should_resume_symplectic(recovery, 0.5_dp)) then
            error stop 'generic recovery did not probe after its cooldown'
        end if
        call activate_rk_recovery(recovery, 0.5_dp)
        recovery%steps = 511
        if (should_resume_symplectic(recovery, 0.5_dp)) then
            error stop 'repeated generic recovery did not back off'
        end if
        recovery%steps = 512
        if (.not. should_resume_symplectic(recovery, 0.5_dp)) then
            error stop 'backed-off generic recovery never resumed probing'
        end if

        recovery = rk_recovery_state_t()
        call activate_rk_recovery(recovery, 0.995_dp)
        if (should_resume_symplectic(recovery, 0.99_dp)) then
            error stop 'edge recovery resumed before moving safely inward'
        end if
        if (.not. should_resume_symplectic(recovery, 0.98_dp)) then
            error stop 'edge recovery did not resume after moving inward'
        end if
    end subroutine test_rk_recovery_regions

    subroutine test_exit_classification
        logical :: classifier_lost
        integer :: classifier_exit

        if (classify_orbit_exit(SYMPLECTIC_STEP_BOUNDARY, ORBIT_GC, 3, .true.) /= &
            ORBIT_EXIT_LCFS) then
            error stop 'converged LCFS event was not classified as physical'
        end if
        if (classify_orbit_exit(SYMPLECTIC_STEP_BOUNDARY, ORBIT_GC, 3, .false.) /= &
            ORBIT_EXIT_NUMERICAL_DOMAIN) then
            error stop 'extended map boundary was classified as physical LCFS'
        end if
        if (classify_orbit_exit(77, ORBIT_GC, 3, .true.) /= ORBIT_EXIT_WALL) then
            error stop 'wall event was not classified as physical'
        end if
        if (classify_orbit_exit(77, ORBIT_GC, 0, .true.) /= ORBIT_EXIT_WALL) then
            error stop 'axis wall event was not classified as physical'
        end if
        if (classify_orbit_exit(1, ORBIT_GC, 3, .true.) /= &
            ORBIT_EXIT_NUMERICAL_DOMAIN) then
            error stop 'exterior Newton iterate was classified as physical'
        end if
        if (classify_orbit_exit(ORBIT_FO_LOSS, ORBIT_FULL_ORBIT, 3, .true.) /= &
            ORBIT_EXIT_LCFS) then
            error stop 'full-orbit LCFS exit was not classified as physical'
        end if
        if (classify_orbit_exit(ORBIT_FO_LOSS, ORBIT_FULL_ORBIT, 3, .false.) /= &
            ORBIT_EXIT_NUMERICAL_DOMAIN) then
            error stop 'extended full-orbit map boundary was classified as physical'
        end if
        if (classify_orbit_exit(ORBIT_FO_NUMERICAL, ORBIT_FULL_ORBIT, 3, .true.) /= &
            ORBIT_EXIT_NUMERICAL_FULL_ORBIT) then
            error stop 'full-orbit locate failure was classified as physical'
        end if
        if (classify_orbit_exit(1, ORBIT_GC, 0, .true.) /= ORBIT_EXIT_LCFS) then
            error stop 'RK LCFS exit was classified as numerical'
        end if
        if (classify_orbit_exit(1, ORBIT_GC, 0, .false.) /= &
            ORBIT_EXIT_NUMERICAL_DOMAIN) then
            error stop 'RK extended-map boundary was classified as physical'
        end if
        if (classify_orbit_exit(2, ORBIT_GC, 0, .true.) /= &
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

end program test_sympl_testfield
