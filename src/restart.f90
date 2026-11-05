module restart_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use finite_check, only: is_finite
    use params, only: ntestpart, times_lost, orbit_exit_code, &
        boundary_event_radial_residual, boundary_event_time_width, &
        trap_par, perp_inv, zend, &
        confpart_pass, confpart_trap, ntimstep, kt_macro, &
        v0, dtaumin, trace_time, ORBIT_EXIT_COMPLETED, &
        ORBIT_EXIT_LCFS, ORBIT_EXIT_WALL, ORBIT_EXIT_NUMERICAL_CONFINED, &
        ORBIT_EXIT_THERMALIZED
    use params, only: zstart, iclass, class_scores, classification_enabled
    implicit none
    private

    public :: particle_done, read_restart_data, restore_confined_counts

    logical, allocatable :: particle_done(:)

contains

    subroutine read_restart_data()
        integer :: i, idx, unit, ios, exit_code
        real(dp) :: tl, tp, s1, pi_val, exit_time, radial_residual, time_width
        real(dp) :: ze(5)

        if (allocated(particle_done)) deallocate(particle_done)
        allocate (particle_done(ntestpart))
        particle_done = .false.

        open (newunit=unit, file='times_lost.dat', status='old', &
              action='read', iostat=ios, recl=1024)
        if (ios /= 0) then
            print *, 'restart: cannot open times_lost.dat, starting from scratch'
            return
        end if

        do i = 1, ntestpart
            read (unit, *, iostat=ios) idx, tl, tp, s1, pi_val, ze
            if (ios /= 0) exit
            if (idx < 1 .or. idx > ntestpart) cycle
            if (.not. is_finite(tl) .or. tl < 0.0d0) cycle
            times_lost(idx) = tl
            if (tl >= trace_time) then
                particle_done(idx) = .true.
                orbit_exit_code(idx) = ORBIT_EXIT_COMPLETED
            end if
            trap_par(idx) = tp
            perp_inv(idx) = pi_val
            zend(:, idx) = ze
        end do

        close (unit)

        open (newunit=unit, file='orbit_exit_code.dat', status='old', &
              action='read', iostat=ios, recl=1024)
        if (ios == 0) then
            do i = 1, ntestpart
                read (unit, *, iostat=ios) idx, exit_code, exit_time, &
                    radial_residual, time_width
                if (ios /= 0) exit
                if (idx < 1 .or. idx > ntestpart) cycle
                if (.not. is_finite(exit_time) .or. exit_time < 0d0) cycle
                if (abs(exit_time - times_lost(idx)) > &
                    16d0*epsilon(exit_time)*max(1d0, abs(exit_time))) cycle
                orbit_exit_code(idx) = exit_code
                boundary_event_radial_residual(idx) = radial_residual
                boundary_event_time_width(idx) = time_width
                particle_done(idx) = exit_code == ORBIT_EXIT_COMPLETED .or. &
                    exit_code == ORBIT_EXIT_LCFS .or. &
                    exit_code == ORBIT_EXIT_WALL .or. &
                    exit_code == ORBIT_EXIT_NUMERICAL_CONFINED .or. &
                    exit_code == ORBIT_EXIT_THERMALIZED
            end do
            close (unit)
        end if
        if (classification_enabled()) call read_restart_classes()
        print *, 'restart: loaded', count(particle_done), '/', ntestpart, &
                 'finished particles'
    end subroutine read_restart_data

    !> Restore the classes of finished markers from class_parts.dat and
    !> class_scores.dat. A finished marker without a class row for the same
    !> start radius is traced again, so class_parts.dat never loses classes.
    subroutine read_restart_classes()
        integer :: i, idx, unit, ios
        integer :: ic(3)
        real(dp) :: s_start, pi_val, tp
        real(dp) :: scores(size(class_scores, 1))
        logical, allocatable :: restored(:)

        allocate (restored(ntestpart), source=.false.)
        open (newunit=unit, file='class_parts.dat', status='old', &
              action='read', iostat=ios, recl=1024)
        if (ios == 0) then
            do i = 1, ntestpart
                read (unit, *, iostat=ios) idx, s_start, pi_val, ic
                if (ios /= 0) exit
                if (idx < 1 .or. idx > ntestpart) cycle
                if (abs(s_start - zstart(1, idx)) > &
                    16d0*epsilon(s_start)*max(1d0, abs(s_start))) cycle
                if (.not. particle_done(idx)) cycle
                iclass(:, idx) = ic
                restored(idx) = .true.
            end do
            close (unit)
        end if
        particle_done = particle_done .and. restored

        open (newunit=unit, file='class_scores.dat', status='old', &
              action='read', iostat=ios, recl=1024)
        if (ios /= 0) return
        do i = 1, ntestpart
            read (unit, *, iostat=ios) idx, scores, tp
            if (ios /= 0) exit
            if (idx < 1 .or. idx > ntestpart) cycle
            if (particle_done(idx)) class_scores(:, idx) = scores
        end do
        close (unit)
    end subroutine read_restart_classes

    subroutine restore_confined_counts()
        integer :: i, it
        real(dp) :: survived_steps
        logical :: passing

        do i = 1, ntestpart
            if (.not. particle_done(i)) cycle
            passing = trap_par(i) < 0.0d0
            survived_steps = times_lost(i)*v0/dtaumin
            do it = 1, ntimstep
                if (orbit_exit_code(i) == ORBIT_EXIT_COMPLETED) then
                    if (real(kt_macro(it), dp) > survived_steps) exit
                else
                    if (real(kt_macro(it), dp) >= survived_steps) exit
                end if
                if (passing) then
                    confpart_pass(it) = confpart_pass(it) + 1.0d0
                else
                    confpart_trap(it) = confpart_trap(it) + 1.0d0
                end if
            end do
        end do
    end subroutine restore_confined_counts

end module restart_mod
