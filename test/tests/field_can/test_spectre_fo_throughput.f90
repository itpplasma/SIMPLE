program test_spectre_fo_throughput
    ! Short stellarator timing witness, shared by the Cartesian baseline and
    ! logical implementation. Timing is reported, never a noisy pass/fail gate.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use field_spectre, only: spectre_field_t, create_spectre_field
    use magfie_sub, only: set_magfie_spectre_field
    use spectre_fo_hybrid, only: spectre_fo_state_t, spectre_fo_configure, &
        spectre_fo_enter, spectre_fo_advance_until_exit, SPECTRE_FO_OK
    implicit none
    integer, parameter :: NMARK = 32, NCHUNKS = 100
    real(dp), parameter :: RO = 100.0_dp, PI = 4.0_dp*atan(1.0_dp)
    type(spectre_field_t) :: field
    type(spectre_fo_state_t) :: state
    character(len=1024) :: path
    real(dp) :: y(5), yo(5), used, dt, rho, started, finished
    real(dp) :: shell(NMARK), endpoint(6, NMARK)
    integer :: status, owner, m, k, volume, failures, clock_rate, clock_start
    logical :: exited

    call get_command_argument(1, path)
    call create_spectre_field(field, trim(path), status)
    if (status /= 0) error stop 'throughput field initialization'
    call set_magfie_spectre_field(field)
    call system_clock(clock_start, clock_rate)
    started = real(clock_start, dp)/real(clock_rate, dp)
    call spectre_fo_configure(1.0e-9_dp)
    call system_clock(clock_start)
    finished = real(clock_start, dp)/real(clock_rate, dp)
    print '(a,f10.4)', 'FO spline setup wall s ', finished - started
    dt = 2.0_dp*PI*RO/1.0e4_dp/10.0_dp
    do volume = 1, 2
        rho = real(volume, dp) - 0.05_dp
        failures = 0
        shell = 0.0_dp
        call system_clock(clock_start, clock_rate)
        started = real(clock_start, dp)/real(clock_rate, dp)
        !$omp parallel do private(m, k, y, yo, state, status, owner, used, exited) &
        !$omp reduction(+:failures)
        do m = 1, NMARK
            y = [rho, 0.7_dp*real(m, dp), 0.3_dp*real(m, dp), 1.0_dp, &
                -0.9_dp + 1.8_dp*real(m - 1, dp)/real(NMARK - 1, dp)]
            call spectre_fo_enter(state, y, volume, volume, RO, status)
            if (status /= SPECTRE_FO_OK) then
                failures = failures + 1
                cycle
            end if
            state%exit_rho = 1.0e6_dp
            do k = 1, NCHUNKS
                call spectre_fo_advance_until_exit(state, dt, RO, used, yo, &
                    owner, exited, status)
                if (status /= SPECTRE_FO_OK .or. exited) exit
                shell(m) = max(shell(m), abs(dot_product(state%v, state%v)/2 - 1))
            end do
            if (status /= SPECTRE_FO_OK .or. exited) failures = failures + 1
            endpoint(:, m) = [state%u, state%v]
        end do
        !$omp end parallel do
        call system_clock(clock_start)
        finished = real(clock_start, dp)/real(clock_rate, dp)
        print '(a,i1,a,f10.4,a,es12.4)', 'volume ', volume, ' wall s ', &
            finished - started, ' shell ', maxval(shell)
        if (failures /= 0) error stop 'throughput trace did not complete'
        if (.not. all(ieee_is_finite(endpoint))) error stop 'nonfinite orbit'
        if (maxval(shell) > 1.0e-6_dp) error stop 'short FO energy drift'
        do m = 1, NMARK
            print '(a,2i4,6es24.15)', 'endpoint ', volume, m, endpoint(:, m)
        end do
    end do
end program test_spectre_fo_throughput
