program test_result_accounting
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use magfie_sub, only: TEST
    use params, only: ntestpart, ntimstep, kt_macro, dtaumin, v0, &
        confpart_pass, confpart_trap, unresolved_orbits, &
        times_lost, orbit_exit_code, trap_par, perp_inv, &
        zstart, zend, boundary_event_radial_residual, &
        boundary_event_time_width, isw_field_type, class_plot, &
        ntcut, ORBIT_EXIT_COMPLETED, ORBIT_EXIT_LCFS, &
        ORBIT_EXIT_NUMERICAL_MAXITER
    use simple_main, only: write_results
    implicit none

    integer :: nerr

    nerr = 0
    call setup_results()
    call test_resolved_denominator()
    call test_zero_resolved()
    call teardown_results()

    if (nerr > 0) error stop 1
    print *, 'All result accounting tests passed!'

contains

    subroutine setup_results()
        ntestpart = 3
        ntimstep = 1
        dtaumin = 1.0_dp
        v0 = 1.0_dp
        isw_field_type = TEST
        class_plot = .false.
        ntcut = 0

        allocate (kt_macro(1), confpart_pass(1), confpart_trap(1))
        allocate (unresolved_orbits(1), times_lost(3), orbit_exit_code(3))
        allocate (trap_par(3), perp_inv(3), zstart(5, 3), zend(5, 3))
        allocate (boundary_event_radial_residual(3))
        allocate (boundary_event_time_width(3))

        kt_macro = [10_int64]
        trap_par = 0.0_dp
        perp_inv = 0.0_dp
        zstart = 0.0_dp
        zend = 0.0_dp
        boundary_event_radial_residual = -1.0_dp
        boundary_event_time_width = -1.0_dp
    end subroutine setup_results

    subroutine test_resolved_denominator()
        integer :: idx, unit
        real(dp) :: time, pass_fraction, trap_fraction, loss_time
        real(dp) :: unresolved_fraction
        integer :: resolved_count, total_count

        confpart_pass = 1.0_dp
        confpart_trap = 0.0_dp
        unresolved_orbits = 1
        times_lost = [10.0_dp, 5.0_dp, 5.0_dp]
        orbit_exit_code = [ORBIT_EXIT_COMPLETED, ORBIT_EXIT_LCFS, &
            ORBIT_EXIT_NUMERICAL_MAXITER]

        call write_results()

        open (newunit=unit, file='confined_fraction.dat', status='old')
        read (unit, *) time, pass_fraction, trap_fraction, resolved_count
        close (unit)
        call assert_close(pass_fraction, 0.5_dp, 'resolved confined fraction')
        call assert_close(trap_fraction, 0.0_dp, 'resolved trapped fraction')
        call assert_true(resolved_count == 2, 'resolved denominator is two')

        open (newunit=unit, file='times_lost.dat', status='old')
        read (unit, *) idx, loss_time
        read (unit, *) idx, loss_time
        read (unit, *) idx, loss_time
        close (unit)
        call assert_true(idx == 3, 'read numerical particle row')
        call assert_true(ieee_is_nan(loss_time), 'numerical loss time is NaN')

        open (newunit=unit, file='unresolved_fraction.dat', status='old')
        read (unit, *) time, unresolved_fraction, total_count
        close (unit)
        call assert_close(unresolved_fraction, 1.0_dp/3.0_dp, &
            'unresolved fraction keeps total denominator')
        call assert_true(total_count == 3, 'unresolved denominator is total')
    end subroutine test_resolved_denominator

    subroutine test_zero_resolved()
        integer :: unit, resolved_count
        real(dp) :: time, pass_fraction, trap_fraction

        confpart_pass = 0.0_dp
        confpart_trap = 0.0_dp
        unresolved_orbits = ntestpart

        call write_results()

        open (newunit=unit, file='confined_fraction.dat', status='old')
        read (unit, *) time, pass_fraction, trap_fraction, resolved_count
        close (unit)
        call assert_true(resolved_count == 0, 'zero resolved denominator')
        call assert_true(ieee_is_nan(pass_fraction), 'zero denominator pass NaN')
        call assert_true(ieee_is_nan(trap_fraction), 'zero denominator trap NaN')
    end subroutine test_zero_resolved

    subroutine teardown_results()
        deallocate (kt_macro, confpart_pass, confpart_trap, unresolved_orbits)
        deallocate (times_lost, orbit_exit_code, trap_par, perp_inv)
        deallocate (zstart, zend, boundary_event_radial_residual)
        deallocate (boundary_event_time_width)
    end subroutine teardown_results

    subroutine assert_close(actual, expected, message)
        real(dp), intent(in) :: actual, expected
        character(*), intent(in) :: message

        if (abs(actual - expected) > 1.0e-12_dp) then
            print *, 'FAIL:', message, 'expected', expected, 'got', actual
            nerr = nerr + 1
        end if
    end subroutine assert_close

    subroutine assert_true(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) then
            print *, 'FAIL:', message
            nerr = nerr + 1
        end if
    end subroutine assert_true

end program test_result_accounting
