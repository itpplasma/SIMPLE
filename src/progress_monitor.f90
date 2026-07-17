module progress_monitor
    !> Periodic progress reporting and team-safe final result output for the
    !> particle tracing loop.
    !>
    !> Once per finished particle a thread calls progress_tick. The OpenMP path
    !> is only one atomic increment: no lock, formatted I/O, file access, clock
    !> query, or callback runs while particles are being traced. The joined
    !> caller emits one final event summary; the owner then writes result files
    !> serially. The legacy interval argument is retained for input/API
    !> compatibility but no longer triggers mid-trace output.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, output_unit
    use diag_counters, only: N_EVENT, diag_counters_total, event_name
    implicit none
    private

    public :: progress_init, progress_tick, progress_finalize
    real(dp) :: t_start = 0.0d0
    integer :: n_total = 0
    integer(int64) :: n_done_shared = 0_int64
    logical :: active = .false.

contains

    subroutine progress_init(interval_seconds, n_particles)
        real(dp), intent(in) :: interval_seconds
        integer, intent(in) :: n_particles

        n_total = n_particles
        n_done_shared = 0_int64
        active = n_particles > 0
        t_start = wall_time()
    end subroutine progress_init

    subroutine progress_tick()
        !> Call once per finished particle from inside the parallel region.
        if (.not. active) return

!$omp atomic update
        n_done_shared = n_done_shared + 1_int64
!$omp end atomic
    end subroutine progress_tick

    subroutine progress_finalize()
        !> Emit one final summary so the event totals (e.g. fo_loss / fo_fault)
        !> are always reported, including for runs too short to cross a periodic
        !> interval. A silent fo_fault count would hide inversion faults.
        if (n_total > 0) call emit(int(n_done_shared), wall_time(), .true.)
        active = .false.
    end subroutine progress_finalize

    subroutine emit(n_done, now, final)
        integer, intent(in) :: n_done
        real(dp), intent(in) :: now
        logical, intent(in) :: final

        real(dp) :: elapsed, frac, eta
        integer :: id
        integer(int64) :: total
        character(len=512) :: events

        elapsed = now - t_start
        frac = real(n_done, dp)/real(max(n_total, 1), dp)
        eta = 0.0d0
        if (frac > 0.0d0) eta = elapsed*(1.0d0 - frac)/frac

        events = ''
        if (final) then
            do id = 1, N_EVENT
                total = diag_counters_total(id)
                if (total > 0_int64) then
                    events = trim(events)//' '//event_name(id)//'='// &
                             trim(int_to_str(total))
                end if
            end do
        end if

        write (output_unit, '(A, F6.2, A, I0, A, I0, A, F0.1, A, F0.1, A, A)') &
            ' [progress] ', 100.0d0*frac, '%  ', n_done, '/', n_total, &
            '  t=', elapsed, 's  eta=', eta, 's', trim(events)
        flush (output_unit)
    end subroutine emit

    function wall_time() result(t)
        !> Portable wall-clock seconds; works with or without OpenMP.
        real(dp) :: t
        integer(int64) :: c, rate

        call system_clock(c, rate)
        t = 0.0d0
        if (rate > 0_int64) t = real(c, dp)/real(rate, dp)
    end function wall_time

    function int_to_str(value) result(str)
        integer(int64), intent(in) :: value
        character(:), allocatable :: str
        character(len=32) :: buf

        write (buf, '(I0)') value
        str = trim(buf)
    end function int_to_str

end module progress_monitor
