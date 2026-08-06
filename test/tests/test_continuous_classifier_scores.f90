program test_continuous_classifier_scores
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use check_orbit_type_sub, only: compute_continuous_tip_scores
    implicit none

    integer :: errors

    errors = 0
    call test_constant_drift_rate(errors)
    call test_rotation_frequency_drift(errors)
    call test_unresolved_reference(errors)
    call test_too_few_tips(errors)

    if (errors > 0) then
        print *, "ERROR:", errors, "continuous classifier score checks failed"
        stop 1
    end if
    print *, "Continuous classifier score checks passed"

contains

    subroutine test_constant_drift_rate(errors)
        integer, intent(inout) :: errors
        integer, parameter :: n = 6
        integer :: jpar_samples, k, rotation_samples
        real(dp) :: jpar_rate, rotation_drift, turns
        real(dp) :: tips(3, n)
        real(dp), parameter :: pi = acos(-1.0_dp)

        tips = 0.0_dp
        do k = 1, n
            tips(2, k) = modulo(real(k - 1, dp)*0.5_dp*pi, 2.0_dp*pi)
        end do
        tips(3, 1) = 5.0_dp
        do k = 2, n
            tips(3, k) = 10.0_dp + real(k - 2, dp)
        end do

        call compute_continuous_tip_scores(tips, n, jpar_rate, &
            rotation_drift, turns, jpar_samples, &
            rotation_samples)

        call assert_close("constant relative drift rate", jpar_rate, 0.4_dp, errors)
        call assert_close("one achieved turn", turns, 1.0_dp, errors)
        call assert_close("constant rotation number", rotation_drift, 0.0_dp, errors)
        call assert_equal("J_parallel sample count", jpar_samples, 4, errors)
        call assert_equal("rotation half sample count", rotation_samples, 2, errors)
    end subroutine test_constant_drift_rate

    subroutine test_rotation_frequency_drift(errors)
        integer, intent(inout) :: errors
        integer, parameter :: n = 9
        integer :: jpar_samples, k, rotation_samples
        real(dp) :: angle, jpar_rate, rotation_drift, turns
        real(dp) :: tips(3, n)
        real(dp), parameter :: pi = acos(-1.0_dp)

        tips = 0.0_dp
        tips(3, :) = 10.0_dp
        angle = 0.0_dp
        do k = 2, 5
            angle = angle + 0.2_dp*pi
            tips(2, k) = modulo(angle, 2.0_dp*pi)
        end do
        do k = 6, n
            angle = angle + 0.4_dp*pi
            tips(2, k) = modulo(angle, 2.0_dp*pi)
        end do

        call compute_continuous_tip_scores(tips, n, jpar_rate, &
            rotation_drift, turns, jpar_samples, &
            rotation_samples)

        call assert_close("rotation drift", rotation_drift, 0.1_dp, errors)
        call assert_close("wrapped-angle exposure", turns, 1.1_dp, errors)
        call assert_close("constant J_parallel", jpar_rate, 0.0_dp, errors)
        call assert_equal("rotation half sample count", rotation_samples, 4, errors)
        call assert_equal("J_parallel sample count", jpar_samples, 7, errors)
    end subroutine test_rotation_frequency_drift

    subroutine test_unresolved_reference(errors)
        integer, intent(inout) :: errors
        integer, parameter :: n = 5
        integer :: jpar_samples, k, rotation_samples
        real(dp) :: jpar_rate, rotation_drift, turns
        real(dp) :: tips(3, n)
        real(dp), parameter :: pi = acos(-1.0_dp)

        tips = 0.0_dp
        do k = 1, n
            tips(2, k) = real(k - 1, dp)*0.25_dp*pi
        end do
        call compute_continuous_tip_scores(tips, n, jpar_rate, &
            rotation_drift, turns, jpar_samples, &
            rotation_samples)

        call assert_close("zero reference score", jpar_rate, 0.0_dp, errors)
        call assert_equal("zero reference samples", jpar_samples, 0, errors)
        call assert_equal("rotation remains resolved", rotation_samples, 2, errors)
    end subroutine test_unresolved_reference

    subroutine test_too_few_tips(errors)
        integer, intent(inout) :: errors
        integer :: jpar_samples, rotation_samples
        real(dp) :: jpar_rate, rotation_drift, turns
        real(dp) :: tips(3, 2)

        tips = 0.0_dp
        call compute_continuous_tip_scores(tips, 2, jpar_rate, &
            rotation_drift, turns, jpar_samples, &
            rotation_samples)

        call assert_close("unresolved J_parallel", jpar_rate, 0.0_dp, errors)
        call assert_close("unresolved rotation", rotation_drift, 0.0_dp, errors)
        call assert_close("unresolved turns", turns, 0.0_dp, errors)
        call assert_equal("unresolved J_parallel samples", jpar_samples, 0, errors)
        call assert_equal("unresolved rotation samples", rotation_samples, 0, errors)
    end subroutine test_too_few_tips

    subroutine assert_close(name, actual, expected, errors)
        character(*), intent(in) :: name
        real(dp), intent(in) :: actual, expected
        integer, intent(inout) :: errors

        if (abs(actual - expected) > 1.0e-12_dp) then
            print *, "FAIL:", name, "actual", actual, "expected", expected
            errors = errors + 1
        end if
    end subroutine assert_close

    subroutine assert_equal(name, actual, expected, errors)
        character(*), intent(in) :: name
        integer, intent(in) :: actual, expected
        integer, intent(inout) :: errors

        if (actual /= expected) then
            print *, "FAIL:", name, "actual", actual, "expected", expected
            errors = errors + 1
        end if
    end subroutine assert_equal

end program test_continuous_classifier_scores
