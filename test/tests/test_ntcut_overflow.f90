program test_ntcut_overflow
    ! The classification cut index ntcut = ntimstep*ntau*tcut/trace_time must stay
    ! disabled (<=0) when tcut<=0, for any trace length. The product ntimstep*ntau
    ! reaches ~2.3e9 for second-scale traces and overflows a 32-bit integer; the
    ! old int32 product wrapped negative and, times tcut<0, flipped ntcut positive,
    ! silently routing long traces through the classifier path. Regression guard.
    use params, only: microstep_cut_index
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none

    integer :: errors
    integer(int64) :: nt

    errors = 0

    ! 1 s W7-X reactor case (npoiper2=16384): ntau ~ 2.3e6, ntimstep*ntau ~ 2.3e9
    ! overflows int32. tcut<0 -> classification disabled -> ntcut must be <= 0.
    nt = microstep_cut_index(1000, 2297297, -1.0d0, 1.0d0)
    if (nt > 0_int64) then
        print *, "FAIL: overflow case gave ntcut > 0:", nt
        errors = errors + 1
    end if

    ! Short trace, tcut<0: also disabled.
    nt = microstep_cut_index(1000, 2295, -1.0d0, 1.0d-3)
    if (nt > 0_int64) then
        print *, "FAIL: short-trace disabled case gave ntcut > 0:", nt
        errors = errors + 1
    end if

    ! tcut>0, no overflow: cut at half the trace -> ntimstep*ntau/2.
    nt = microstep_cut_index(1000, 1000, 0.5d0, 1.0d0)
    if (nt /= 500000_int64) then
        print *, "FAIL: enabled no-overflow case gave", nt, "expected 500000"
        errors = errors + 1
    end if

    ! tcut>0 with an overflow-scale product must stay correct in int64, not wrap.
    nt = microstep_cut_index(1000, 2297297, 1.0d0, 1.0d0)
    if (nt /= 2297297000_int64) then
        print *, "FAIL: enabled overflow-scale case gave", nt, &
            "expected 2297297000"
        errors = errors + 1
    end if

    if (errors == 0) then
        print *, "All ntcut overflow tests passed!"
    else
        print *, "ERROR:", errors, "test(s) failed!"
        stop 1
    end if
end program test_ntcut_overflow
