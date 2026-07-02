module progress_monitor
    !> Periodic progress reporting and partial-result checkpointing for the
    !> particle tracing loop.
    !>
    !> Once per finished particle a thread calls progress_tick. The hot path is
    !> a single atomic increment of the completed counter and a wall-clock read,
    !> with no critical and no file access. When `interval` seconds have passed
    !> one thread enters a rarely taken critical, writes a status line and the
    !> partial results to disk, and advances the deadline. A run killed mid-flight
    !> therefore leaves its most recent output on disk instead of nothing.
    !>
    !> The monitor never touches physics arrays itself. The owner registers a
    !> dump callback at init, so this module depends on no tracing code and the
    !> dependency runs one way.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, output_unit
    use diag_counters, only: N_EVENT, diag_counters_total, event_name
    implicit none
    private

    public :: progress_init, progress_tick, progress_finalize, dump_proc_i

    abstract interface
        subroutine dump_proc_i()
            !> Write current results to disk from the shared result arrays.
        end subroutine dump_proc_i
    end interface

    procedure(dump_proc_i), pointer :: dump_results => null()
    real(dp) :: interval = 0.0d0   ! <= 0 disables periodic checkpointing
    real(dp) :: t_start = 0.0d0
    real(dp) :: t_next = 0.0d0
    integer :: n_total = 0
    integer(int64) :: n_done_shared = 0_int64
    logical :: active = .false.

contains

    subroutine progress_init(interval_seconds, n_particles, dump)
        real(dp), intent(in) :: interval_seconds
        integer, intent(in) :: n_particles
        procedure(dump_proc_i) :: dump

        interval = interval_seconds
        n_total = n_particles
        n_done_shared = 0_int64
        dump_results => dump
        active = interval > 0.0d0 .and. n_particles > 0
        t_start = wall_time()
        t_next = t_start + interval
    end subroutine progress_init

    subroutine progress_tick()
        !> Call once per finished particle from inside the parallel region.
        real(dp) :: now
        integer(int64) :: done

        if (.not. active) return

!$omp atomic capture
        n_done_shared = n_done_shared + 1_int64
        done = n_done_shared
!$omp end atomic

        now = wall_time()
        if (now < t_next) return

!$omp critical (progress_flush)
        if (now >= t_next) then   ! double-checked: only one thread per interval
            t_next = now + interval
            call emit(int(done), now)
        end if
!$omp end critical (progress_flush)
    end subroutine progress_tick

    subroutine progress_finalize()
        !> Emit one final summary so the event totals (e.g. fo_loss / fo_fault)
        !> are always reported, including for runs too short to cross a periodic
        !> interval. A silent fo_fault count would hide inversion faults.
        if (n_total > 0) call emit(int(n_done_shared), wall_time())
        active = .false.
        dump_results => null()
    end subroutine progress_finalize

    subroutine emit(n_done, now)
        integer, intent(in) :: n_done
        real(dp), intent(in) :: now

        real(dp) :: elapsed, frac, eta
        integer :: id
        integer(int64) :: total
        character(len=512) :: events

        if (associated(dump_results)) call dump_results()

        elapsed = now - t_start
        frac = real(n_done, dp)/real(max(n_total, 1), dp)
        eta = 0.0d0
        if (frac > 0.0d0) eta = elapsed*(1.0d0 - frac)/frac

        events = ''
        do id = 1, N_EVENT
            total = diag_counters_total(id)
            if (total > 0_int64) then
                events = trim(events)//' '//event_name(id)//'='// &
                         trim(int_to_str(total))
            end if
        end do

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
