module restart_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use params, only: ntestpart, times_lost, trap_par, perp_inv, zend, &
                      confpart_pass, confpart_trap, ntimstep, kt_macro, &
                      v0, dtaumin
    implicit none
    private

    public :: particle_done, read_restart_data, restore_confined_counts

    logical, allocatable :: particle_done(:)

contains

    subroutine read_restart_data()
        integer :: i, idx, unit, ios
        real(dp) :: tl, tp, s1, pi_val
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
            if (tl < 0.0d0) cycle
            particle_done(idx) = .true.
            times_lost(idx) = tl
            trap_par(idx) = tp
            perp_inv(idx) = pi_val
            zend(:, idx) = ze
        end do

        close (unit)
        print *, 'restart: loaded', count(particle_done), '/', ntestpart, &
                 'finished particles'
    end subroutine read_restart_data

    subroutine restore_confined_counts()
        integer :: i, it
        integer(int64) :: kt_survived
        logical :: passing

        do i = 1, ntestpart
            if (.not. particle_done(i)) cycle
            passing = trap_par(i) < 0.0d0
            kt_survived = nint(times_lost(i)*v0/dtaumin, kind=int64)
            do it = 1, ntimstep
                if (kt_macro(it) > kt_survived) exit
                if (passing) then
                    confpart_pass(it) = confpart_pass(it) + 1.0d0
                else
                    confpart_trap(it) = confpart_trap(it) + 1.0d0
                end if
            end do
        end do
    end subroutine restore_confined_counts

end module restart_mod
