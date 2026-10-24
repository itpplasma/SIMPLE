program test_continuous_classifier_scores
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use check_orbit_type_sub, only: check_orbit_type, compute_continuous_tip_scores
    use detect_oneline_mod, only: jpar_variation_rate, &
                                  rotation_number_drift, precession_turns, &
                                  jpar_sample_count, rotation_half_count, &
                                  tip_count, jpar_spread, jpar_ref, &
                                  topology_margin, score_status
    implicit none

    integer :: errors, k
    real(dp) :: legacy_scores(12, 3), legacy_parts(6, 3)
    character(len=32) :: argument

    errors = 0
    call test_constant_drift_rate(errors)
    call test_rotation_frequency_drift(errors)
    call test_unresolved_reference(errors)
    call test_nonfinite_tip(errors)
    call test_too_few_tips(errors)
    call test_regular_tip_maps(errors, legacy_scores, legacy_parts)

    if (errors > 0) then
        print *, "ERROR:", errors, "continuous classifier score checks failed"
        stop 1
    end if
    print *, "Continuous classifier score checks passed"
    if (command_argument_count() == 1) then
        call get_command_argument(1, argument)
        if (trim(argument) /= '--legacy-controls') error stop 'unknown argument'
        do k = 1, 3
            write (*, '(a,12(1x,es25.16e3))') 'legacy_score', legacy_scores(:, k)
            write (*, '(a,6(1x,es25.16e3))') 'legacy_parts', legacy_parts(:, k)
        end do
    end if

contains

    subroutine test_regular_tip_maps(errors, scores, parts)
        integer, intent(inout) :: errors
        real(dp), intent(out) :: scores(12, 3), parts(6, 3)
        real(dp), parameter :: step = 0.1_dp, twopi = 2.0_dp*acos(-1.0_dp)
        real(dp), parameter :: spreads(3) = [0.0_dp, 15.0_dp, 16.0_dp]
        integer, parameter :: expected_classes(3) = [1, 1, 2]
        integer, parameter :: turns = 4, max_tips = ceiling(turns*twopi/step) + 2
        real(dp) :: tip(3)
        integer :: control, nfp, ideal, jpar, status, k

        ! A rigid rotation on one invariant circle has return intervals of
        ! 62 or 63 tips, well above the early-stochastic cutoff of ten.
        ! Four recurrences preserve their angular order, so topology is regular.
        ! Constant J is regular. Exact spreads at and above the unchanged 15 cm
        ! threshold test its strict comparison without a sampled trajectory.
        do control = 1, 3
            nfp = 0
            ideal = 0
            jpar = 0
            do k = 1, max_tips
                tip = [0.3_dp, modulo(real(k - 1, dp)*step, twopi), 20.0_dp]
                if (k >= 3) tip(3) = tip(3) + spreads(control)
                call check_orbit_type(turns, nfp, tip, ideal, jpar, status)
                if (score_status /= 0) exit
            end do
            call assert_equal('controlled topology', ideal, 1, errors)
            call assert_equal('controlled J_parallel', jpar, &
                              expected_classes(control), errors)
            call assert_equal('controlled score status', score_status, 1, errors)
            call assert_equal('controlled resolved tips', tip_count, &
                              ceiling(turns*twopi/step) + 1, errors)
            call assert_close('controlled J_parallel spread', jpar_spread, &
                              spreads(control), errors)
            if (topology_margin <= 0.0_dp) then
                print *, 'FAIL: controlled topology margin is not positive'
                errors = errors + 1
            end if
            scores(:, control) = [real(control, dp), &
                                  jpar_variation_rate, rotation_number_drift, &
                                  precession_turns, real(jpar_sample_count, dp), &
                                  real(rotation_half_count, dp), &
                                  real(tip_count, dp), jpar_spread, jpar_ref, &
                                  topology_margin, real(score_status, dp), 1.0_dp]
            parts(:, control) = [real(control, dp), tip(1), 0.0_dp, &
                                 real(jpar, dp), real(ideal, dp), 0.0_dp]
        end do
    end subroutine test_regular_tip_maps

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

    subroutine test_nonfinite_tip(errors)
        integer, intent(inout) :: errors
        integer, parameter :: n = 5
        integer :: jpar_samples, rotation_samples
        real(dp) :: jpar_rate, rotation_drift, turns
        real(dp) :: tips(3, n)

        tips = 0.0_dp
        tips(2, :) = [0.0_dp, 0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp]
        tips(3, :) = 1.0_dp
        tips(3, 4) = ieee_value(0.0_dp, ieee_quiet_nan)

        call compute_continuous_tip_scores(tips, n, jpar_rate, &
                                  rotation_drift, turns, jpar_samples, rotation_samples)

        if (.not. all(ieee_is_finite([jpar_rate, rotation_drift, turns]))) then
            print *, "FAIL: nonfinite tip produced a nonfinite score"
            errors = errors + 1
        end if
        call assert_close("nonfinite J_parallel score", jpar_rate, 0.0_dp, errors)
        call assert_close("nonfinite rotation score", rotation_drift, 0.0_dp, errors)
        call assert_close("nonfinite turns", turns, 0.0_dp, errors)
        call assert_equal("nonfinite J_parallel samples", jpar_samples, 0, errors)
        call assert_equal("nonfinite rotation samples", rotation_samples, 0, errors)
    end subroutine test_nonfinite_tip

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
