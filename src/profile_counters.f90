module profile_counters
    use, intrinsic :: iso_fortran_env, only: int64, real64, output_unit

    implicit none
    private

    public :: prof_reset
    public :: prof_print
    public :: prof_add_old_newton_call
    public :: prof_add_old_newton_iter
    public :: prof_add_old_f_eval
    public :: prof_add_soa_newton_call
    public :: prof_add_soa_newton_iter
    public :: prof_add_soa_eval_pre
    public :: prof_add_soa_eval_newton

    integer(int64), save :: old_newton_calls = 0_int64
    integer(int64), save :: old_newton_iters = 0_int64
    integer(int64), save :: old_newton_iters_max = 0_int64
    integer(int64), save :: old_f_eval_calls = 0_int64

    integer(int64), save :: soa_newton_calls = 0_int64
    integer(int64), save :: soa_newton_iters = 0_int64
    integer(int64), save :: soa_newton_iters_max = 0_int64
    integer(int64), save :: soa_eval_pre_calls = 0_int64
    integer(int64), save :: soa_eval_pre_points = 0_int64
    integer(int64), save :: soa_eval_newton_calls = 0_int64
    integer(int64), save :: soa_eval_newton_points = 0_int64

contains

    subroutine prof_reset
        old_newton_calls = 0_int64
        old_newton_iters = 0_int64
        old_newton_iters_max = 0_int64
        old_f_eval_calls = 0_int64

        soa_newton_calls = 0_int64
        soa_newton_iters = 0_int64
        soa_newton_iters_max = 0_int64
        soa_eval_pre_calls = 0_int64
        soa_eval_pre_points = 0_int64
        soa_eval_newton_calls = 0_int64
        soa_eval_newton_points = 0_int64
    end subroutine prof_reset


    subroutine prof_add_old_newton_call
        old_newton_calls = old_newton_calls + 1_int64
    end subroutine prof_add_old_newton_call


    subroutine prof_add_old_newton_iter(niter)
        integer, intent(in) :: niter
        old_newton_iters = old_newton_iters + 1_int64
        old_newton_iters_max = max(old_newton_iters_max, int(niter, int64))
    end subroutine prof_add_old_newton_iter


    subroutine prof_add_old_f_eval
        old_f_eval_calls = old_f_eval_calls + 1_int64
    end subroutine prof_add_old_f_eval


    subroutine prof_add_soa_newton_call
        soa_newton_calls = soa_newton_calls + 1_int64
    end subroutine prof_add_soa_newton_call


    subroutine prof_add_soa_newton_iter(niter)
        integer, intent(in) :: niter
        soa_newton_iters = soa_newton_iters + 1_int64
        soa_newton_iters_max = max(soa_newton_iters_max, int(niter, int64))
    end subroutine prof_add_soa_newton_iter


    subroutine prof_add_soa_eval_pre(npts)
        integer, intent(in) :: npts
        soa_eval_pre_calls = soa_eval_pre_calls + 1_int64
        soa_eval_pre_points = soa_eval_pre_points + int(npts, int64)
    end subroutine prof_add_soa_eval_pre


    subroutine prof_add_soa_eval_newton(npts)
        integer, intent(in) :: npts
        soa_eval_newton_calls = soa_eval_newton_calls + 1_int64
        soa_eval_newton_points = soa_eval_newton_points + int(npts, int64)
    end subroutine prof_add_soa_eval_newton


    subroutine prof_print(unit)
        integer, intent(in), optional :: unit
        integer :: u
        real(real64) :: old_iters_per_call, soa_iters_per_call
        real(real64) :: soa_pre_pts_per_call, soa_newton_pts_per_call

        u = output_unit
        if (present(unit)) u = unit

        old_iters_per_call = 0.0_real64
        if (old_newton_calls > 0_int64) then
            old_iters_per_call = real(old_newton_iters, real64) / &
                                 real(old_newton_calls, real64)
        end if

        soa_iters_per_call = 0.0_real64
        if (soa_newton_calls > 0_int64) then
            soa_iters_per_call = real(soa_newton_iters, real64) / &
                                 real(soa_newton_calls, real64)
        end if

        soa_pre_pts_per_call = 0.0_real64
        if (soa_eval_pre_calls > 0_int64) then
            soa_pre_pts_per_call = real(soa_eval_pre_points, real64) / &
                                   real(soa_eval_pre_calls, real64)
        end if

        soa_newton_pts_per_call = 0.0_real64
        if (soa_eval_newton_calls > 0_int64) then
            soa_newton_pts_per_call = real(soa_eval_newton_points, real64) / &
                                      real(soa_eval_newton_calls, real64)
        end if

        write (u, *) 'PROFILE_COUNTERS'
        write (u, *) '  old_newton_calls:      ', old_newton_calls
        write (u, *) '  old_newton_iters:      ', old_newton_iters
        write (u, *) '  old_newton_iters_max:  ', old_newton_iters_max
        write (u, *) '  old_newton_iters_per_call: ', old_iters_per_call
        write (u, *) '  old_f_sympl_euler1_calls:  ', old_f_eval_calls
        write (u, *) '  soa_newton_calls:      ', soa_newton_calls
        write (u, *) '  soa_newton_iters:      ', soa_newton_iters
        write (u, *) '  soa_newton_iters_max:  ', soa_newton_iters_max
        write (u, *) '  soa_newton_iters_per_call: ', soa_iters_per_call
        write (u, *) '  soa_eval_pre_calls:    ', soa_eval_pre_calls
        write (u, *) '  soa_eval_pre_points:   ', soa_eval_pre_points
        write (u, *) '  soa_eval_pre_pts_per_call: ', soa_pre_pts_per_call
        write (u, *) '  soa_eval_newton_calls: ', soa_eval_newton_calls
        write (u, *) '  soa_eval_newton_points:', soa_eval_newton_points
        write (u, *) '  soa_eval_newton_pts_per_call: ', soa_newton_pts_per_call
    end subroutine prof_print

end module profile_counters
