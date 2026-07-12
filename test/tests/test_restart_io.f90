program test_restart_io
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use params, only: ntestpart, times_lost, trap_par, perp_inv, zend, zstart, &
                      orbit_exit_code, boundary_event_radial_residual, &
                      boundary_event_time_width, ORBIT_EXIT_COMPLETED, &
                      ORBIT_EXIT_LCFS, &
                      confpart_pass, confpart_trap, ntimstep, kt_macro, &
                      v0, dtaumin, trace_time
    use restart_mod, only: particle_done, read_restart_data, restore_confined_counts
    implicit none

    integer, parameter :: np = 8, nt = 5
    integer :: nerr

    nerr = 0
    call test_read_restart()
    call test_legacy_restart_retraces_early_exit()
    call test_restore_counts()

    if (nerr > 0) then
        print *, 'FAILED:', nerr, 'errors'
        error stop 1
    end if
    print *, 'All restart IO tests passed!'

contains

    subroutine test_read_restart()
        integer :: i, unit
        real(dp) :: expected_tl(np), expected_tp(np), expected_pi(np)
        real(dp) :: expected_ze(5, np)

        print *, '=== Test: read_restart_data ==='

        ntestpart = np
        trace_time = 1.0d0
        allocate (times_lost(np), trap_par(np), perp_inv(np))
        allocate (orbit_exit_code(np), boundary_event_radial_residual(np), &
                  boundary_event_time_width(np))
        allocate (zend(5, np), zstart(5, np))
        times_lost = -1.0d0
        orbit_exit_code = 0
        boundary_event_radial_residual = -1.0d0
        boundary_event_time_width = -1.0d0
        trap_par = 0.0d0
        perp_inv = 0.0d0
        zend = 0.0d0
        zstart = 0.0d0

        expected_tl = [-1.0d0, 0.05d0, -1.0d0, 0.10d0, &
                       -1.0d0, 0.02d0, -1.0d0, 0.08d0]
        expected_tp = [0.0d0, -0.3d0, 0.0d0, 0.5d0, &
                       0.0d0, -0.1d0, 0.0d0, 0.7d0]
        expected_pi = [0.0d0, 1.2d0, 0.0d0, 0.8d0, &
                       0.0d0, 1.5d0, 0.0d0, 0.9d0]
        expected_ze = 0.0d0
        do i = 1, np
            if (expected_tl(i) > 0.0d0) then
                expected_ze(1, i) = 0.1d0*i
                expected_ze(2, i) = 0.2d0*i
                expected_ze(3, i) = 0.3d0*i
                expected_ze(4, i) = 1.0d0
                expected_ze(5, i) = 0.4d0*i
            end if
        end do

        open (newunit=unit, file='times_lost.dat', recl=1024)
        do i = 1, np
            write (unit, *) i, expected_tl(i), expected_tp(i), zstart(1, i), &
                expected_pi(i), expected_ze(:, i)
        end do
        close (unit)

        open (newunit=unit, file='orbit_exit_code.dat', recl=1024)
        do i = 1, np
            write (unit, *) i, merge(ORBIT_EXIT_LCFS, 0, expected_tl(i) > 0d0), &
                expected_tl(i), &
                merge(1.0d-11, -1.0d0, expected_tl(i) > 0d0), &
                merge(1.0d-12, -1.0d0, expected_tl(i) > 0d0)
        end do
        close (unit)

        call read_restart_data()

        do i = 1, np
            if (expected_tl(i) < 0.0d0) then
                call assert(.not. particle_done(i), 'unfinished not done', i)
                call assert(times_lost(i) < 0.0d0, 'unfinished tl=-1', i)
            else
                call assert(particle_done(i), 'finished is done', i)
                call assert_eq(times_lost(i), expected_tl(i), 'times_lost', i)
                call assert_eq(trap_par(i), expected_tp(i), 'trap_par', i)
                call assert_eq(perp_inv(i), expected_pi(i), 'perp_inv', i)
                call assert_eq(zend(1, i), expected_ze(1, i), 'zend(1)', i)
                call assert_eq(zend(5, i), expected_ze(5, i), 'zend(5)', i)
                call assert(orbit_exit_code(i) == ORBIT_EXIT_LCFS, &
                    'physical exit code restored', i)
                call assert_eq(boundary_event_radial_residual(i), 1.0d-11, &
                    'event residual', i)
            end if
        end do

        deallocate (times_lost, orbit_exit_code, boundary_event_radial_residual, &
                    boundary_event_time_width, trap_par, perp_inv, zend, zstart, &
                    particle_done)
    end subroutine test_read_restart

    subroutine test_legacy_restart_retraces_early_exit()
        integer :: unit

        ntestpart = 2
        trace_time = 1d0
        allocate (times_lost(2), trap_par(2), perp_inv(2), orbit_exit_code(2))
        allocate (boundary_event_radial_residual(2), boundary_event_time_width(2))
        allocate (zend(5, 2), zstart(5, 2))
        times_lost = -1d0
        orbit_exit_code = ORBIT_EXIT_COMPLETED
        trap_par = 0d0
        perp_inv = 0d0
        zend = 0d0
        zstart = 0d0

        open (newunit=unit, file='times_lost.dat', recl=1024)
        write (unit, *) 1, 0.5d0, 0d0, 0d0, 0d0, zend(:, 1)
        write (unit, *) 2, 1.0d0, 0d0, 0d0, 0d0, zend(:, 2)
        close (unit)
        open (newunit=unit, file='orbit_exit_code.dat', status='old')
        close (unit, status='delete')

        call read_restart_data()
        call assert(.not. particle_done(1), &
            'ambiguous legacy early exit must be retraced', 1)
        call assert(particle_done(2), 'legacy completed trace is reusable', 2)

        deallocate (times_lost, orbit_exit_code, boundary_event_radial_residual, &
                    boundary_event_time_width, trap_par, perp_inv, zend, zstart, &
                    particle_done)
    end subroutine test_legacy_restart_retraces_early_exit

    subroutine test_restore_counts()
        integer :: it

        print *, '=== Test: restore_confined_counts ==='

        ntestpart = 4
        ntimstep = nt
        v0 = 1.0d0
        dtaumin = 1.0d-3

        allocate (times_lost(4), trap_par(4), perp_inv(4), orbit_exit_code(4))
        allocate (zend(5, 4), zstart(5, 4))
        allocate (confpart_pass(nt), confpart_trap(nt))
        allocate (kt_macro(nt))

        kt_macro = [0_int64, 100_int64, 200_int64, 300_int64, 400_int64]

        confpart_pass = 0.0d0
        confpart_trap = 0.0d0
        times_lost = -1.0d0
        orbit_exit_code = ORBIT_EXIT_LCFS
        trap_par = 0.0d0
        perp_inv = 0.0d0
        zend = 0.0d0
        zstart = 0.0d0

        if (allocated(particle_done)) deallocate(particle_done)
        allocate (particle_done(4))

        ! Particle 1: passing, exits exactly at step 200.
        particle_done(1) = .true.
        times_lost(1) = 200.0d0*dtaumin/v0
        trap_par(1) = -0.5d0

        ! Particle 2: trapped, survived all 400+ microsteps → counted at all it
        particle_done(2) = .true.
        times_lost(2) = 500.0d0*dtaumin/v0
        trap_par(2) = 0.3d0
        orbit_exit_code(2) = ORBIT_EXIT_COMPLETED

        ! Particle 3: not done
        particle_done(3) = .false.

        ! Particle 4: passing, exits just after step 200.
        particle_done(4) = .true.
        times_lost(4) = 200.25d0*dtaumin/v0
        trap_par(4) = -0.2d0

        call restore_confined_counts()

        ! it=1: kt_macro=0 → p1(pass)+p2(trap)+p4(pass) → pass=2, trap=1
        call assert_eq(confpart_pass(1), 2.0d0, 'pass(1)', 0)
        call assert_eq(confpart_trap(1), 1.0d0, 'trap(1)', 0)

        ! At step 100, both passing particles remain confined.
        call assert_eq(confpart_pass(2), 2.0d0, 'pass(2)', 0)
        call assert_eq(confpart_trap(2), 1.0d0, 'trap(2)', 0)

        ! At step 200, p1 is already lost; p4 exits just afterward.
        call assert_eq(confpart_pass(3), 1.0d0, 'pass(3)', 0)
        call assert_eq(confpart_trap(3), 1.0d0, 'trap(3)', 0)

        ! it=4: kt_macro=300 → p1(250<300,skip)+p2(500>=300,trap) → pass=0, trap=1
        call assert_eq(confpart_pass(4), 0.0d0, 'pass(4)', 0)
        call assert_eq(confpart_trap(4), 1.0d0, 'trap(4)', 0)

        ! it=5: kt_macro=400 → p2(500>=400,trap) → pass=0, trap=1
        call assert_eq(confpart_pass(5), 0.0d0, 'pass(5)', 0)
        call assert_eq(confpart_trap(5), 1.0d0, 'trap(5)', 0)

        deallocate (times_lost, orbit_exit_code, trap_par, perp_inv, zend, zstart)
        deallocate (confpart_pass, confpart_trap, kt_macro, particle_done)
    end subroutine test_restore_counts

    subroutine assert(cond, msg, idx)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        integer, intent(in) :: idx

        if (.not. cond) then
            print *, 'FAIL:', msg, 'particle', idx
            nerr = nerr + 1
        end if
    end subroutine assert

    subroutine assert_eq(actual, expected, msg, idx)
        real(dp), intent(in) :: actual, expected
        character(*), intent(in) :: msg
        integer, intent(in) :: idx

        if (abs(actual - expected) > 1.0d-12) then
            print *, 'FAIL:', msg, 'particle', idx, &
                     'expected', expected, 'got', actual
            nerr = nerr + 1
        end if
    end subroutine assert_eq

end program test_restart_io
