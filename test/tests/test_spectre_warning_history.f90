program test_spectre_warning_history
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use params, only: ntestpart, ntimstep, nper, npoiper, reallocate_arrays, &
               kt_macro, dtaumin, v0, confpart_pass, confpart_trap, unresolved_orbits, &
                      times_lost, orbit_exit_code, orbit_warning_hold_count, &
              orbit_first_warning_step, orbit_terminal_exit_code, orbit_terminal_time, &
                      zstart, zend, class_plot, ntcut, ORBIT_EXIT_NUMERICAL_EVENT
    use simple_main, only: record_spectre_warning_holds, &
                           finalize_spectre_warning_holds, write_results
    implicit none
    integer :: i, it, first_unresolved, count

    ntestpart = 4
    ntimstep = 3
    nper = 1
    npoiper = 1
    call reallocate_arrays
    allocate (kt_macro(3))
    kt_macro = [0_int64, 1_int64, 2_int64]
    dtaumin = 1.0_dp
    v0 = 1.0_dp
    class_plot = .false.
    ntcut = 0
    zstart = 0.0_dp
    zend = 0.0_dp
    confpart_pass = 0.0_dp
    confpart_trap = 0.0_dp
    unresolved_orbits = 0

    ! Synthetic native recovery histories: held completion, held then physical
    ! exit, ordinary completion and immediate physical exit. Each worker owns
    ! one marker, as in trace_parallel; shared time-bin counters are atomic.
    orbit_exit_code = [0, 1, 0, 1]
    times_lost = [2.0_dp, 1.75_dp, 2.0_dp, 0.5_dp]
    !$omp parallel do private(it, first_unresolved, count)
    do i = 1, ntestpart
        first_unresolved = 0
        do it = 1, ntimstep
            count = 0
            if (i == 1 .and. it == 2) count = 1
            if (i == 1 .and. it == 3) count = 3
            if (i == 2 .and. it == 3) count = 2
            call record_spectre_warning_holds(i, count, it, first_unresolved)
            if (first_unresolved == 0 .and. &
                real(it - 1, dp) <= times_lost(i)) then
                !$omp atomic update
                confpart_pass(it) = confpart_pass(it) + 1.0_dp
            end if
        end do
        if (first_unresolved > 0) then
            do it = first_unresolved, ntimstep
                !$omp atomic update
                unresolved_orbits(it) = unresolved_orbits(it) + 1
            end do
        end if
        call finalize_spectre_warning_holds(i)
    end do
    !$omp end parallel do

    if (any(orbit_warning_hold_count /= [3, 2, 0, 0])) &
        error stop 'hold histories mixed between markers'
    if (any(orbit_first_warning_step /= [2, 3, 0, 0])) &
        error stop 'first hold time was overwritten'
    if (any(orbit_terminal_exit_code /= [0, 1, 0, 1])) &
        error stop 'terminal physical labels were erased'
    if (maxval(abs(orbit_terminal_time - [2.0_dp, 1.75_dp, 2.0_dp, 0.5_dp])) &
        > 1.0e-14_dp) error stop 'terminal times were erased'
    if (any(orbit_exit_code(1:2) /= ORBIT_EXIT_NUMERICAL_EVENT)) &
        error stop 'held trajectories remained resolved'
    if (.not. all(ieee_is_nan(times_lost(1:2)))) &
        error stop 'held trajectory exposed a physical loss time'
    if (any(orbit_exit_code(3:4) /= [0, 1])) &
        error stop 'ordinary completion or physical loss changed'
    if (any(unresolved_orbits /= [0, 1, 2])) &
        error stop 'unresolved time accounting changed'
    call write_results
    print *, 'SPECTRE unresolved history and terminal events PASS'
end program test_spectre_warning_history
