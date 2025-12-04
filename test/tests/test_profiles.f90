program test_profiles
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use simple_profiles
    use collis_alp
    implicit none

    logical :: all_passed
    integer :: n_failed

    all_passed = .true.
    n_failed = 0

    call test_power_series_evaluation(all_passed, n_failed)
    call test_two_power_evaluation(all_passed, n_failed)
    call test_flat_profile_coefficients(all_passed, n_failed)
    call test_peaked_profile_collision_rates(all_passed, n_failed)
    call plot_collision_frequency_comparison()

    if (all_passed) then
        print *, 'All profile tests passed'
    else
        print *, 'FAILED: ', n_failed, ' profile tests failed'
        error stop 1
    end if

contains

    subroutine test_power_series_evaluation(passed, nfail)
        logical, intent(inout) :: passed
        integer, intent(inout) :: nfail

        real(dp) :: coef(MAX_POWER_SERIES), result, expected
        real(dp), parameter :: tol = 1.0d-12

        print *, 'Testing power series evaluation...'

        coef = 0.0d0
        coef(1) = 1.0d0
        coef(2) = -1.0d0
        result = eval_power_series(0.5d0, 1.0d0, coef)
        expected = 0.5d0
        if (abs(result - expected) > tol) then
            print *, '  FAIL: 1-s at s=0.5 expected', expected, 'got', result
            passed = .false.
            nfail = nfail + 1
        else
            print *, '  PASS: 1-s at s=0.5'
        end if

        coef = 0.0d0
        coef(1) = 1.0d0
        coef(2) = -2.0d0
        coef(3) = 1.0d0
        result = eval_power_series(0.5d0, 1.0d0, coef)
        expected = 0.25d0
        if (abs(result - expected) > tol) then
            print *, '  FAIL: (1-s)^2 at s=0.5 expected', expected, 'got', result
            passed = .false.
            nfail = nfail + 1
        else
            print *, '  PASS: (1-s)^2 at s=0.5'
        end if

        result = eval_power_series(0.0d0, 100.0d0, coef)
        expected = 100.0d0
        if (abs(result - expected) > tol) then
            print *, '  FAIL: scale at s=0 expected', expected, 'got', result
            passed = .false.
            nfail = nfail + 1
        else
            print *, '  PASS: scale at s=0'
        end if
    end subroutine test_power_series_evaluation

    subroutine test_two_power_evaluation(passed, nfail)
        logical, intent(inout) :: passed
        integer, intent(inout) :: nfail

        real(dp) :: result, expected
        real(dp), parameter :: tol = 1.0d-12

        print *, 'Testing two-power evaluation...'

        result = eval_two_power(0.5d0, 1.0d0, 1.0d0, 2.0d0)
        expected = 0.25d0
        if (abs(result - expected) > tol) then
            print *, '  FAIL: (1-s)^2 at s=0.5 expected', expected, 'got', result
            passed = .false.
            nfail = nfail + 1
        else
            print *, '  PASS: (1-s)^2 at s=0.5'
        end if

        result = eval_two_power(0.5d0, 1.0d0, 1.0d0, 3.5d0)
        expected = 0.5d0**3.5d0
        if (abs(result - expected) > tol) then
            print *, '  FAIL: (1-s)^3.5 at s=0.5 expected', expected, 'got', result
            passed = .false.
            nfail = nfail + 1
        else
            print *, '  PASS: (1-s)^3.5 at s=0.5'
        end if

        result = eval_two_power(0.0d0, 10010.0d0, 1.0d0, 2.0d0)
        expected = 10010.0d0
        if (abs(result - expected) > tol) then
            print *, '  FAIL: on-axis value expected', expected, 'got', result
            passed = .false.
            nfail = nfail + 1
        else
            print *, '  PASS: on-axis value'
        end if

        result = eval_two_power(1.0d0, 10010.0d0, 1.0d0, 2.0d0)
        expected = 0.0d0
        if (abs(result - expected) > tol) then
            print *, '  FAIL: edge value expected', expected, 'got', result
            passed = .false.
            nfail = nfail + 1
        else
            print *, '  PASS: edge value is zero'
        end if

        result = eval_two_power(1.0d0, 10010.0d0, 1.0d0, 0.0d0)
        expected = 10010.0d0
        if (abs(result - expected) > tol) then
            print *, '  FAIL: flat profile at edge expected', expected, 'got', result
            passed = .false.
            nfail = nfail + 1
        else
            print *, '  PASS: flat profile at edge gives scale'
        end if
    end subroutine test_two_power_evaluation

    subroutine test_flat_profile_coefficients(passed, nfail)
        logical, intent(inout) :: passed
        integer, intent(inout) :: nfail

        real(dp) :: am1, am2, Z1, Z2, ealpha, v0
        real(dp) :: densi1, densi2, tempi1, tempi2, tempe
        real(dp) :: dchichi, slowrate, dchichi_norm, slowrate_norm
        real(dp) :: efcolf_scalar(nsorts), velrat_scalar(nsorts), enrat_scalar(nsorts)
        real(dp) :: efcolf_profile(nsorts), velrat_profile(nsorts), enrat_profile(nsorts)
        real(dp), parameter :: tol = 1.0d-10
        integer :: i

        print *, 'Testing flat profile vs scalar coefficients...'

        am1 = 2.0d0
        am2 = 3.0d0
        Z1 = 1.0d0
        Z2 = 1.0d0
        densi1 = 0.5d14
        densi2 = 0.5d14
        tempi1 = 1.0d4
        tempi2 = 1.0d4
        tempe = 1.0d4
        ealpha = 3.5d6

        call loacol_alpha(am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, &
                          ealpha, v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
        efcolf_scalar = efcolf
        velrat_scalar = velrat
        enrat_scalar = enrat

        profile_type = "two_power"
        active_profile = TWO_POWER
        Te_scale = tempe
        Te_p1 = 1.0d0
        Te_p2 = 0.0d0
        Ti1_scale = tempi1
        Ti1_p1 = 1.0d0
        Ti1_p2 = 0.0d0
        Ti2_scale = tempi2
        Ti2_p1 = 1.0d0
        Ti2_p2 = 0.0d0
        ni1_scale = densi1*1.0d6
        ni1_p1 = 1.0d0
        ni1_p2 = 0.0d0
        ni2_scale = densi2*1.0d6
        ni2_p1 = 1.0d0
        ni2_p2 = 0.0d0

        call init_collision_profiles(am1, am2, Z1, Z2, ealpha, v0)

        call get_local_coeffs(0.5d0, efcolf_profile, velrat_profile, enrat_profile)

        do i = 1, nsorts
            if (abs(efcolf_profile(i) - efcolf_scalar(i))/abs(efcolf_scalar(i)) > tol) then
                print *, '  FAIL: efcolf mismatch for species', i
                print *, '    scalar:', efcolf_scalar(i), 'profile:', efcolf_profile(i)
                passed = .false.
                nfail = nfail + 1
            end if
            if (abs(velrat_profile(i) - velrat_scalar(i))/abs(velrat_scalar(i)) > tol) then
                print *, '  FAIL: velrat mismatch for species', i
                passed = .false.
                nfail = nfail + 1
            end if
            if (abs(enrat_profile(i) - enrat_scalar(i))/abs(enrat_scalar(i)) > tol) then
                print *, '  FAIL: enrat mismatch for species', i
                passed = .false.
                nfail = nfail + 1
            end if
        end do

        print *, '  PASS: flat profile coefficients match scalar'
    end subroutine test_flat_profile_coefficients

    subroutine test_peaked_profile_collision_rates(passed, nfail)
        logical, intent(inout) :: passed
        integer, intent(inout) :: nfail

        real(dp) :: am1, am2, Z1, Z2, ealpha, v0
        real(dp) :: densi1, densi2, tempi1, tempi2, tempe
        real(dp) :: dchichi, slowrate, dchichi_norm, slowrate_norm
        real(dp) :: efcolf_core(nsorts), efcolf_edge(nsorts)
        real(dp) :: velrat_core(nsorts), velrat_edge(nsorts)
        real(dp) :: enrat_core(nsorts), enrat_edge(nsorts)
        real(dp) :: s_core, s_edge

        print *, 'Testing peaked profile produces different collision rates...'

        am1 = 2.0d0
        am2 = 3.0d0
        Z1 = 1.0d0
        Z2 = 1.0d0
        densi1 = 0.5d14
        densi2 = 0.5d14
        tempi1 = 1.0d4
        tempi2 = 1.0d4
        tempe = 1.0d4
        ealpha = 3.5d6

        call loacol_alpha(am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, &
                          ealpha, v0, dchichi, slowrate, dchichi_norm, slowrate_norm)

        profile_type = "two_power"
        active_profile = TWO_POWER
        Te_scale = 10010.0d0
        Te_p1 = 1.0d0
        Te_p2 = 2.0d0
        Ti1_scale = 10010.0d0
        Ti1_p1 = 1.0d0
        Ti1_p2 = 2.0d0
        Ti2_scale = 0.0d0
        Ti2_p1 = 1.0d0
        Ti2_p2 = 0.0d0
        ni1_scale = 1.822d21
        ni1_p1 = 1.0d0
        ni1_p2 = 3.5d0
        ni2_scale = 0.0d0
        ni2_p1 = 1.0d0
        ni2_p2 = 0.0d0

        call init_collision_profiles(am1, am2, Z1, Z2, ealpha, v0)

        s_core = 0.04d0
        s_edge = 0.64d0

        call get_local_coeffs(s_core, efcolf_core, velrat_core, enrat_core)
        call get_local_coeffs(s_edge, efcolf_edge, velrat_edge, enrat_edge)

        if (efcolf_core(1) <= efcolf_edge(1)) then
            print *, '  FAIL: core collision rate should be higher than edge'
            print *, '    core efcolf(1):', efcolf_core(1)
            print *, '    edge efcolf(1):', efcolf_edge(1)
            passed = .false.
            nfail = nfail + 1
        else
            print *, '  PASS: core has higher collision rate than edge'
            print *, '    core/edge ratio:', efcolf_core(1)/efcolf_edge(1)
        end if
    end subroutine test_peaked_profile_collision_rates

    subroutine plot_collision_frequency_comparison()
        use fortplot, only: figure, plot, savefig, xlabel, ylabel, title, legend

        real(dp) :: am1, am2, Z1, Z2, ealpha, v0
        real(dp) :: densi1, densi2, tempi1, tempi2, tempe
        real(dp) :: dchichi, slowrate, dchichi_norm, slowrate_norm
        real(dp) :: efcolf_loc(nsorts), velrat_loc(nsorts), enrat_loc(nsorts)
        real(dp) :: efcolf_onaxis
        integer, parameter :: npts = 101
        real(dp) :: s_arr(npts), efcolf_profile(npts), efcolf_flat(npts)
        real(dp) :: s, ds
        integer :: i

        print *, 'Generating collision frequency comparison plot...'

        am1 = 2.0d0
        am2 = 3.0d0
        Z1 = 1.0d0
        Z2 = 1.0d0
        ealpha = 3.5d6

        profile_type = "two_power"
        active_profile = TWO_POWER
        Te_scale = 10010.0d0
        Te_p1 = 1.0d0
        Te_p2 = 2.0d0
        Ti1_scale = 10010.0d0
        Ti1_p1 = 1.0d0
        Ti1_p2 = 2.0d0
        Ti2_scale = 0.0d0
        Ti2_p1 = 1.0d0
        Ti2_p2 = 0.0d0
        ni1_scale = 1.822d21
        ni1_p1 = 1.0d0
        ni1_p2 = 3.5d0
        ni2_scale = 0.0d0
        ni2_p1 = 1.0d0
        ni2_p2 = 0.0d0

        densi1 = ni1_scale*1.0d-6
        densi2 = ni2_scale*1.0d-6
        tempi1 = Ti1_scale
        tempi2 = Ti2_scale
        tempe = Te_scale

        call loacol_alpha(am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, &
                          ealpha, v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
        efcolf_onaxis = efcolf(1)

        call init_collision_profiles(am1, am2, Z1, Z2, ealpha, v0)

        ds = 1.0d0/dble(npts - 1)
        do i = 1, npts
            s = dble(i - 1)*ds
            s_arr(i) = s
            call get_local_coeffs(s, efcolf_loc, velrat_loc, enrat_loc)
            efcolf_profile(i) = efcolf_loc(1)
            efcolf_flat(i) = efcolf_onaxis
        end do

        call figure()
        call plot(s_arr, efcolf_profile, label="Peaked profile", linestyle="b-")
        call plot(s_arr, efcolf_flat, label="Uniform (on-axis)", linestyle="r--")
        call xlabel("s (normalized toroidal flux)")
        call ylabel("Collision frequency [1/s]")
        call title("Collision Frequency: Profile vs Uniform")
        call legend()
        call savefig("collision_frequency_comparison.png")

        print *, '  Saved plot to collision_frequency_comparison.png'
    end subroutine plot_collision_frequency_comparison

end program test_profiles
