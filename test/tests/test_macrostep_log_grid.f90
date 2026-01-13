! Test macrostep log grid implementation
! Verifies that logarithmic macrostep time grid:
! 1. Preserves total microsteps exactly
! 2. Produces log10-spaced distribution (10:1 ratio last:first)
! 3. Runs in O(n) time, not O(n^2) or worse
program test_macrostep_log_grid
    use iso_fortran_env, only: int64
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    call test_log_grid_total_steps_preserved()
    call test_log_grid_distribution()
    call test_log_grid_performance()
    call test_log_grid_edge_case_min_microsteps()

    print *, "All macrostep log grid tests passed!"

contains

    subroutine test_log_grid_total_steps_preserved()
        ! Test that logarithmic macrostep grid preserves total microsteps exactly
        integer, parameter :: ntimstep = 20000
        integer, parameter :: ntau = 100

        integer :: nintv, i
        integer(int64) :: n_microsteps_total, kt_target, kt_prev
        integer, allocatable :: ntau_macro(:)
        real(dp) :: weight_sum, cumul_weight, w

        print *, "=== Test: Log grid preserves total microsteps ==="

        nintv = ntimstep - 1
        n_microsteps_total = int(nintv, kind=int64) * int(ntau, kind=int64)

        allocate(ntau_macro(ntimstep))
        ntau_macro = 0

        ! Compute log grid using the clean cumulative algorithm
        weight_sum = 0.0d0
        do i = 1, nintv
            weight_sum = weight_sum + 10.0d0 ** (dble(i - 1) / dble(max(1, nintv - 1)))
        end do

        cumul_weight = 0.0d0
        kt_prev = 0_int64
        do i = 1, nintv
            w = 10.0d0 ** (dble(i - 1) / dble(max(1, nintv - 1)))
            cumul_weight = cumul_weight + w
            kt_target = nint(cumul_weight / weight_sum * dble(n_microsteps_total), kind=int64)
            ntau_macro(i + 1) = int(kt_target - kt_prev)
            kt_prev = kt_target
        end do

        print *, "n_microsteps_total =", n_microsteps_total
        print *, "Final kt_target    =", kt_target
        print *, "First 5 ntau_macro:", ntau_macro(1:5)
        print *, "Last 5 ntau_macro:", ntau_macro(ntimstep-4:ntimstep)

        ! Total must match exactly
        if (kt_target /= n_microsteps_total) then
            print *, "FAIL: Total microsteps mismatch"
            print *, "  Expected:", n_microsteps_total
            print *, "  Got:", kt_target
            stop 1
        end if

        print *, "PASS: Total microsteps preserved exactly"
        deallocate(ntau_macro)
    end subroutine test_log_grid_total_steps_preserved

    subroutine test_log_grid_distribution()
        ! Test that the distribution has approximately 10:1 ratio last:first
        integer, parameter :: ntimstep = 1000
        integer, parameter :: ntau = 100

        integer :: nintv, i
        integer(int64) :: n_microsteps_total, kt_target, kt_prev
        integer, allocatable :: ntau_macro(:)
        real(dp) :: weight_sum, cumul_weight, w, ratio

        print *, ""
        print *, "=== Test: Log grid has correct distribution ==="

        nintv = ntimstep - 1
        n_microsteps_total = int(nintv, kind=int64) * int(ntau, kind=int64)

        allocate(ntau_macro(ntimstep))
        ntau_macro = 0

        weight_sum = 0.0d0
        do i = 1, nintv
            weight_sum = weight_sum + 10.0d0 ** (dble(i - 1) / dble(max(1, nintv - 1)))
        end do

        cumul_weight = 0.0d0
        kt_prev = 0_int64
        do i = 1, nintv
            w = 10.0d0 ** (dble(i - 1) / dble(max(1, nintv - 1)))
            cumul_weight = cumul_weight + w
            kt_target = nint(cumul_weight / weight_sum * dble(n_microsteps_total), kind=int64)
            ntau_macro(i + 1) = int(kt_target - kt_prev)
            kt_prev = kt_target
        end do

        ! Check ratio between last and first interval
        ratio = dble(ntau_macro(ntimstep)) / dble(ntau_macro(2))
        print *, "First interval ntau:", ntau_macro(2)
        print *, "Last interval ntau:", ntau_macro(ntimstep)
        print *, "Ratio last/first:", ratio

        ! Should be close to 10
        if (abs(ratio - 10.0d0) > 0.5d0) then
            print *, "FAIL: Ratio should be close to 10"
            stop 1
        end if

        print *, "PASS: Distribution has correct 10:1 ratio"
        deallocate(ntau_macro)
    end subroutine test_log_grid_distribution

    subroutine test_log_grid_performance()
        ! Test that log grid computation is O(n), not O(n^2) or worse
        integer, parameter :: ntimstep = 20000
        integer, parameter :: ntau = 100

        integer :: nintv, i
        integer(int64) :: n_microsteps_total, kt_target, kt_prev
        integer, allocatable :: ntau_macro(:)
        real(dp) :: weight_sum, cumul_weight, w
        real(dp) :: t_start, t_end

        print *, ""
        print *, "=== Test: Log grid performance ==="

        nintv = ntimstep - 1
        n_microsteps_total = int(nintv, kind=int64) * int(ntau, kind=int64)

        allocate(ntau_macro(ntimstep))

        call cpu_time(t_start)

        ! This should complete in milliseconds, not minutes
        weight_sum = 0.0d0
        do i = 1, nintv
            weight_sum = weight_sum + 10.0d0 ** (dble(i - 1) / dble(max(1, nintv - 1)))
        end do

        cumul_weight = 0.0d0
        kt_prev = 0_int64
        ntau_macro = 0
        do i = 1, nintv
            w = 10.0d0 ** (dble(i - 1) / dble(max(1, nintv - 1)))
            cumul_weight = cumul_weight + w
            kt_target = nint(cumul_weight / weight_sum * dble(n_microsteps_total), kind=int64)
            ntau_macro(i + 1) = int(kt_target - kt_prev)
            kt_prev = kt_target
        end do

        call cpu_time(t_end)

        print *, "Log grid computation time:", t_end - t_start, "seconds"
        print *, "Number of intervals:", nintv

        ! Should complete in under 1 second (O(n) algorithm)
        if (t_end - t_start > 1.0d0) then
            print *, "FAIL: Log grid computation took too long (possible O(n^2) or infinite loop)"
            stop 1
        end if

        print *, "PASS: Log grid computation is fast (O(n))"
        deallocate(ntau_macro)
    end subroutine test_log_grid_performance

    subroutine test_log_grid_edge_case_min_microsteps()
        ! Edge case: very few microsteps per interval
        ! This tests the max(1, ...) safeguard ensuring no empty intervals
        integer, parameter :: ntimstep = 100
        integer, parameter :: ntau = 1  ! Only 1 microstep per interval on average!

        integer :: nintv, i, min_ntau
        integer(int64) :: n_microsteps_total, kt_target, kt_prev, total
        integer, allocatable :: ntau_macro(:), kt_macro(:)
        real(dp) :: weight_sum, cumul_weight, w

        print *, ""
        print *, "=== Test: Edge case - minimum microsteps per interval ==="

        nintv = ntimstep - 1
        n_microsteps_total = int(nintv, kind=int64) * int(ntau, kind=int64)

        allocate(ntau_macro(ntimstep))
        allocate(kt_macro(ntimstep))
        ntau_macro = 0
        kt_macro = 0_int64

        ! Compute log grid with max(1, ...) safeguard
        weight_sum = 0.0d0
        do i = 1, nintv
            weight_sum = weight_sum + 10.0d0 ** (dble(i - 1) / dble(max(1, nintv - 1)))
        end do

        cumul_weight = 0.0d0
        kt_prev = 0_int64
        do i = 1, nintv
            w = 10.0d0 ** (dble(i - 1) / dble(max(1, nintv - 1)))
            cumul_weight = cumul_weight + w
            kt_target = nint(cumul_weight / weight_sum * dble(n_microsteps_total), kind=int64)
            ntau_macro(i + 1) = max(1, int(kt_target - kt_prev))
            kt_macro(i + 1) = kt_prev + int(ntau_macro(i + 1), kind=int64)
            kt_prev = kt_macro(i + 1)
        end do

        ! Check that no interval has 0 microsteps
        min_ntau = minval(ntau_macro(2:ntimstep))
        print *, "n_microsteps_total =", n_microsteps_total
        print *, "Actual total       =", kt_macro(ntimstep)
        print *, "Min ntau_macro     =", min_ntau
        print *, "First 5 ntau_macro:", ntau_macro(1:5)

        if (min_ntau < 1) then
            print *, "FAIL: Some intervals have 0 microsteps"
            stop 1
        end if

        ! Total may exceed original due to max(1,...) forcing
        total = kt_macro(ntimstep)
        if (total < n_microsteps_total) then
            print *, "FAIL: Total microsteps less than expected"
            stop 1
        end if

        print *, "PASS: All intervals have at least 1 microstep"
        deallocate(ntau_macro, kt_macro)
    end subroutine test_log_grid_edge_case_min_microsteps

end program test_macrostep_log_grid
