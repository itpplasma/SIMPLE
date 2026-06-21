module diag_counters
    !> Lock-free per-thread event counters for hot-loop diagnostics.
    !>
    !> The symplectic integrators hit non-convergence and out-of-domain cases
    !> deep inside their inner Newton iterations. Printing there floods stdout
    !> and serializes threads. Instead each event bumps a counter; progress_monitor
    !> aggregates and reports them. Every thread owns its own padded column, so an
    !> increment touches only that thread's cache line and needs no atomic or
    !> critical. Totals are a plain reduction over the thread columns, taken
    !> outside the hot path.
    use, intrinsic :: iso_fortran_env, only: int64
!$  use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    implicit none
    private

    public :: EVT_NEWTON1_MAXIT, EVT_NEWTON2_MAXIT, EVT_RK_GAUSS_MAXIT, &
              EVT_RK_LOBATTO_MAXIT, EVT_FIXPOINT_MAXIT, EVT_R_NEGATIVE, &
              EVT_CPP_LU_FAIL, EVT_CPP_NONCONV, EVT_CPP_SBOUND, N_EVENT
    public :: diag_counters_init, count_event, diag_counters_total, &
              diag_counters_reset, event_name

    integer, parameter :: EVT_NEWTON1_MAXIT = 1
    integer, parameter :: EVT_NEWTON2_MAXIT = 2
    integer, parameter :: EVT_RK_GAUSS_MAXIT = 3
    integer, parameter :: EVT_RK_LOBATTO_MAXIT = 4
    integer, parameter :: EVT_FIXPOINT_MAXIT = 5
    integer, parameter :: EVT_R_NEGATIVE = 6
    ! 6D canonical CP/CPP midpoint Newton outcomes, kept separate from physical
    ! edge loss so a numerical failure is never silently counted as a lost
    ! particle. SBOUND = genuine s >= 1 crossing (physical loss); LU_FAIL =
    ! singular Newton matrix; NONCONV = residual not converged in maxit.
    integer, parameter :: EVT_CPP_LU_FAIL = 7
    integer, parameter :: EVT_CPP_NONCONV = 8
    integer, parameter :: EVT_CPP_SBOUND = 9
    integer, parameter :: N_EVENT = 9

    ! One cache line per thread column, so neighbouring threads never share a
    ! line. The event id indexes within a column; STRIDE >= N_EVENT.
    integer, parameter :: STRIDE = 16
    integer(int64), allocatable :: counts(:, :)  ! (STRIDE, 0:nthreads-1)

contains

    subroutine diag_counters_init()
        integer :: nthreads

        nthreads = 1
!$      nthreads = omp_get_max_threads()
        if (allocated(counts)) deallocate (counts)
        allocate (counts(STRIDE, 0:nthreads - 1))
        counts = 0_int64
    end subroutine diag_counters_init

    subroutine count_event(id)
        !> Record one occurrence of event `id`. Safe to call from inside a
        !> parallel region: writes only the calling thread's column.
        integer, intent(in) :: id
        integer :: tid

        if (.not. allocated(counts)) return
        tid = 0
!$      tid = omp_get_thread_num()
        counts(id, tid) = counts(id, tid) + 1_int64
    end subroutine count_event

    function diag_counters_total(id) result(total)
        !> Sum of event `id` across all threads. Call outside the hot path.
        integer, intent(in) :: id
        integer(int64) :: total

        total = 0_int64
        if (allocated(counts)) total = sum(counts(id, :))
    end function diag_counters_total

    subroutine diag_counters_reset()
        if (allocated(counts)) counts = 0_int64
    end subroutine diag_counters_reset

    function event_name(id) result(name)
        integer, intent(in) :: id
        character(:), allocatable :: name

        select case (id)
        case (EVT_NEWTON1_MAXIT)
            name = 'newton1_maxit'
        case (EVT_NEWTON2_MAXIT)
            name = 'newton2_maxit'
        case (EVT_RK_GAUSS_MAXIT)
            name = 'rk_gauss_maxit'
        case (EVT_RK_LOBATTO_MAXIT)
            name = 'rk_lobatto_maxit'
        case (EVT_FIXPOINT_MAXIT)
            name = 'fixpoint_maxit'
        case (EVT_R_NEGATIVE)
            name = 'r_negative'
        case (EVT_CPP_LU_FAIL)
            name = 'cpp_lu_fail'
        case (EVT_CPP_NONCONV)
            name = 'cpp_nonconv'
        case (EVT_CPP_SBOUND)
            name = 'cpp_sbound'
        case default
            name = 'unknown'
        end select
    end function event_name

end module diag_counters
