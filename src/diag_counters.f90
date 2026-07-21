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
        EVT_FO_LOSS, EVT_FO_FAULT, EVT_MIDPOINT_MAXIT, &
        EVT_WARNING_STEP_SKIP, EVT_SPECTRE_REF_INVERSE_MAXIT, &
        EVT_SPECTRE_INVALID_STATE, EVT_SYMPLECTIC_RK_RECOVERY, &
        EVT_SYMPLECTIC_RESUME, EVT_WARNING_MAXIT_ACCEPT, &
        EVT_WARNING_MAXIT_REJECT_100, EVT_WARNING_MAXIT_REJECT_1E4, &
        EVT_WARNING_MAXIT_REJECT_LARGE, EVT_WARNING_MAXIT_NONFINITE, &
        EVT_RETRY_EXHAUSTED, EVT_WARNING_AXIS_CROSSING_ACCEPT, &
        EVT_TOROIDAL_REGULARIZATION, EVT_AXIS_FO_BRIDGE, N_EVENT
    public :: diag_counters_init, count_event, diag_counters_total, &
        diag_counters_reset, event_name

    integer, parameter :: EVT_NEWTON1_MAXIT = 1
    integer, parameter :: EVT_NEWTON2_MAXIT = 2
    integer, parameter :: EVT_RK_GAUSS_MAXIT = 3
    integer, parameter :: EVT_RK_LOBATTO_MAXIT = 4
    integer, parameter :: EVT_FIXPOINT_MAXIT = 5
    integer, parameter :: EVT_R_NEGATIVE = 6
    ! Full-orbit (FO) outcomes, kept separate from physical edge loss so a numerical
    ! locate fault is never silently counted as a lost particle. LOSS = guiding-centre
    ! s >= 1 crossing (physical loss); FAULT = field-inversion non-convergence.
    integer, parameter :: EVT_FO_LOSS = 7
    integer, parameter :: EVT_FO_FAULT = 8
    integer, parameter :: EVT_MIDPOINT_MAXIT = 9
    integer, parameter :: EVT_WARNING_STEP_SKIP = 10
    integer, parameter :: EVT_SPECTRE_REF_INVERSE_MAXIT = 11
    integer, parameter :: EVT_SPECTRE_INVALID_STATE = 12
    integer, parameter :: EVT_SYMPLECTIC_RK_RECOVERY = 13
    integer, parameter :: EVT_SYMPLECTIC_RESUME = 14
    integer, parameter :: EVT_WARNING_MAXIT_ACCEPT = 15
    integer, parameter :: EVT_WARNING_MAXIT_REJECT_100 = 16
    integer, parameter :: EVT_WARNING_MAXIT_REJECT_1E4 = 17
    integer, parameter :: EVT_WARNING_MAXIT_REJECT_LARGE = 18
    integer, parameter :: EVT_WARNING_MAXIT_NONFINITE = 19
    integer, parameter :: EVT_RETRY_EXHAUSTED = 20
    integer, parameter :: EVT_WARNING_AXIS_CROSSING_ACCEPT = 21
    integer, parameter :: EVT_TOROIDAL_REGULARIZATION = 22
    integer, parameter :: EVT_AXIS_FO_BRIDGE = 23
    integer, parameter :: N_EVENT = 23

    ! Whole cache lines per thread column keep neighbouring threads from sharing
    ! a line. The event id indexes within a column; STRIDE >= N_EVENT.
    integer, parameter :: STRIDE = 32
    integer(int64), allocatable :: counts(:, :) ! (STRIDE, 0:nthreads-1)

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
        case (EVT_FO_LOSS)
            name = 'fo_loss'
        case (EVT_FO_FAULT)
            name = 'fo_fault'
        case (EVT_MIDPOINT_MAXIT)
            name = 'midpoint_maxit'
        case (EVT_WARNING_STEP_SKIP)
            name = 'warning_step_skip'
        case (EVT_SPECTRE_REF_INVERSE_MAXIT)
            name = 'spectre_ref_inverse_maxit'
        case (EVT_SPECTRE_INVALID_STATE)
            name = 'spectre_invalid_state'
        case (EVT_SYMPLECTIC_RK_RECOVERY)
            name = 'symplectic_rk_recovery'
        case (EVT_SYMPLECTIC_RESUME)
            name = 'symplectic_resume'
        case (EVT_WARNING_MAXIT_ACCEPT)
            name = 'warning_maxit_accept'
        case (EVT_WARNING_MAXIT_REJECT_100)
            name = 'warning_maxit_reject_100'
        case (EVT_WARNING_MAXIT_REJECT_1E4)
            name = 'warning_maxit_reject_1e4'
        case (EVT_WARNING_MAXIT_REJECT_LARGE)
            name = 'warning_maxit_reject_large'
        case (EVT_WARNING_MAXIT_NONFINITE)
            name = 'warning_maxit_nonfinite'
        case (EVT_RETRY_EXHAUSTED)
            name = 'retry_exhausted'
        case (EVT_WARNING_AXIS_CROSSING_ACCEPT)
            name = 'warning_axis_crossing_accept'
        case (EVT_TOROIDAL_REGULARIZATION)
            name = 'toroidal_regularization'
        case (EVT_AXIS_FO_BRIDGE)
            name = 'axis_fo_bridge'
        case default
            name = 'unknown'
        end select
    end function event_name

end module diag_counters
