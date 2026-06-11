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
    call test_flat_intermediate_values(all_passed, n_failed)
    call test_slowing_down_distribution(all_passed, n_failed)
    call test_loss_fraction_flat_vs_scalar(all_passed, n_failed)
    call test_maxwellian_fixed_point(all_passed, n_failed)
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

    ! Flat radial profiles must reproduce the scalar (constant-coefficient)
    ! collision operator bit-for-bit. This checks every intermediate value the
    ! stochastic step depends on: the interpolated coefficients at off-node s
    ! and the coleff outputs dpp/dhh/fpeff over the whole momentum range.
    subroutine test_flat_intermediate_values(passed, nfail)
        logical, intent(inout) :: passed
        integer, intent(inout) :: nfail

        real(dp) :: am1, am2, Z1, Z2, ealpha, v0
        real(dp) :: densi1, densi2, tempi1, tempi2, tempe
        real(dp) :: dchichi, slowrate, dchichi_norm, slowrate_norm
        real(dp) :: efs(nsorts), vrs(nsorts), ers(nsorts)
        real(dp) :: efl(nsorts), vrl(nsorts), erl(nsorts)
        real(dp) :: dpp_s, dhh_s, fp_s, dpp_l, dhh_l, fp_l
        real(dp) :: s, p
        integer :: i, j, k
        logical :: coeff_ok, oper_ok

        print *, 'Testing flat profile reproduces scalar at all intermediate values...'

        am1 = 2.0d0; am2 = 3.0d0; Z1 = 1.0d0; Z2 = 1.0d0; ealpha = 3.5d6
        ni1_scale = 0.5d20; ni2_scale = 0.5d20
        Te_scale = 1.2d4; Ti1_scale = 1.0d4; Ti2_scale = 0.8d4
        densi1 = ni1_scale*1.0d-6; densi2 = ni2_scale*1.0d-6
        tempi1 = Ti1_scale; tempi2 = Ti2_scale; tempe = Te_scale

        call loacol_alpha(am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, &
                          ealpha, v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
        efs = efcolf; vrs = velrat; ers = enrat

        call set_flat_two_power()
        call init_collision_profiles(am1, am2, Z1, Z2, ealpha, v0)

        ! Interpolated coefficients must equal the scalar values exactly, including
        ! at s between grid nodes where the interpolation weight is nontrivial.
        coeff_ok = .true.
        do j = 0, 200
            s = (dble(j) + 0.37d0)/200.5d0
            if (s > 1.0d0) cycle
            call get_local_coeffs(s, efl, vrl, erl)
            do i = 1, nsorts
                if (efl(i) /= efs(i) .or. vrl(i) /= vrs(i) .or. erl(i) /= ers(i)) &
                    coeff_ok = .false.
            end do
        end do
        if (coeff_ok) then
            print *, '  PASS: efcolf/velrat/enrat exact at all s'
        else
            print *, '  FAIL: interpolated coefficients differ from scalar'
            passed = .false.; nfail = nfail + 1
        end if

        ! The diffusion/drag outputs must match between the constant operator
        ! (coleff, module coefficients) and the local operator (coleff_local,
        ! interpolated coefficients) over the full momentum range.
        oper_ok = .true.
        s = 0.413d0
        call get_local_coeffs(s, efl, vrl, erl)
        do k = 1, 200
            p = dble(k)/100.0d0
            call coleff(p, dpp_s, dhh_s, fp_s)
            call coleff_local(p, efl, vrl, erl, dpp_l, dhh_l, fp_l)
            if (dpp_l /= dpp_s .or. dhh_l /= dhh_s .or. fp_l /= fp_s) oper_ok = .false.
        end do
        if (oper_ok) then
            print *, '  PASS: dpp/dhh/fpeff exact at all p'
        else
            print *, '  FAIL: coleff_local differs from coleff'
            passed = .false.; nfail = nfail + 1
        end if
    end subroutine test_flat_intermediate_values

    ! The deterministic drag operator slows a particle from v0 toward the thermal
    ! bath. Residence time per momentum bin is the steady-state slowing-down
    ! distribution n(p) for a constant source. The classical (Spitzer) result is
    ! n(p) ~ p^2/(p^3 + p_c^3) with p_c the critical velocity. This checks that
    ! shape and that the flat-profile and constant operators give the same
    ! distribution bit-for-bit.
    subroutine test_slowing_down_distribution(passed, nfail)
        logical, intent(inout) :: passed
        integer, intent(inout) :: nfail

        integer, parameter :: nbin = 14
        real(dp), parameter :: plo = 0.30d0, phi = 0.92d0
        real(dp), parameter :: pstop = 0.25d0
        real(dp) :: am1, am2, Z1, Z2, ealpha, v0
        real(dp) :: densi1, densi2, tempi1, tempi2, tempe
        real(dp) :: dchichi, slowrate, dchichi_norm, slowrate_norm
        real(dp) :: hist_scalar(nbin), hist_flat(nbin), hist_ana(nbin)
        real(dp) :: pcen(nbin), dbin
        real(dp) :: dpp, dhh, fpeff, dtauc
        real(dp) :: f_hi, f_lo, p_hi, p_lo, acoef, xc3, xc, model, relmax
        integer :: ib
        logical :: same, shape_ok, form_ok

        print *, 'Testing slowing-down distribution vs analytic...'

        am1 = 2.0d0; am2 = 3.0d0; Z1 = 1.0d0; Z2 = 1.0d0; ealpha = 3.5d6
        ni1_scale = 0.5d20; ni2_scale = 0.5d20
        Te_scale = 1.0d4; Ti1_scale = 1.0d4; Ti2_scale = 1.0d4
        densi1 = ni1_scale*1.0d-6; densi2 = ni2_scale*1.0d-6
        tempi1 = Ti1_scale; tempi2 = Ti2_scale; tempe = Te_scale

        call loacol_alpha(am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, &
                          ealpha, v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
        call set_flat_two_power()

        dbin = (phi - plo)/dble(nbin)
        do ib = 1, nbin
            pcen(ib) = plo + (dble(ib) - 0.5d0)*dbin
        end do

        call coleff(0.5d0, dpp, dhh, fpeff)
        dtauc = 0.004d0*dbin/abs(fpeff)

        ! Variant "old constant mode": grid is the constant scalar coefficients.
        call fill_constant_grid()
        call drag_residence(dtauc, plo, phi, nbin, pstop, hist_scalar)

        ! Variant "new flat-profile mode": grid built from the flat profile.
        call init_collision_profiles(am1, am2, Z1, Z2, ealpha, v0)
        call drag_residence(dtauc, plo, phi, nbin, pstop, hist_flat)

        call normalize(hist_scalar, nbin)
        call normalize(hist_flat, nbin)

        ! Both variants must yield the same slowing-down distribution exactly.
        same = .true.
        do ib = 1, nbin
            if (hist_flat(ib) /= hist_scalar(ib)) same = .false.
        end do
        if (same) then
            print *, '  PASS: flat-profile and constant distributions identical'
        else
            print *, '  FAIL: slowing-down distribution differs between variants'
            passed = .false.; nfail = nfail + 1
        end if

        ! Critical velocity p_c from the analytic drag form |fpeff|=A(p^3+p_c^3)/p^2.
        p_hi = 0.90d0; p_lo = 0.35d0
        call coleff(p_hi, dpp, dhh, f_hi)
        f_hi = abs(f_hi)
        call coleff(p_lo, dpp, dhh, f_lo)
        f_lo = abs(f_lo)
        acoef = (f_hi*p_hi**2 - f_lo*p_lo**2)/(p_hi**3 - p_lo**3)
        xc3 = (f_hi*p_hi**2)/acoef - p_hi**3
        xc = sign(abs(xc3)**(1.0d0/3.0d0), xc3)

        ! Analytic slowing-down distribution n(p) ~ p^2/(p^3 + p_c^3).
        do ib = 1, nbin
            hist_ana(ib) = pcen(ib)**2/(pcen(ib)**3 + xc3)
        end do
        call normalize(hist_ana, nbin)

        shape_ok = .true.; relmax = 0.0d0
        do ib = 1, nbin
            relmax = max(relmax, abs(hist_scalar(ib) - hist_ana(ib))/hist_ana(ib))
        end do
        if (relmax > 0.04d0) shape_ok = .false.
        if (shape_ok) then
            print *, '  PASS: distribution matches p^2/(p^3+pc^3), max reldev', relmax
        else
            print *, '  FAIL: distribution deviates from analytic, max reldev', relmax
            passed = .false.; nfail = nfail + 1
        end if

        ! The drag operator must be of Spitzer slowing-down form with a physical
        ! critical velocity (between thermal and birth speed).
        form_ok = (xc > 0.2d0 .and. xc < 0.5d0)
        relmax = 0.0d0
        do ib = 1, nbin
            call coleff(pcen(ib), dpp, dhh, fpeff)
            model = acoef*(pcen(ib)**3 + xc3)/pcen(ib)**2
            relmax = max(relmax, abs(abs(fpeff) - model)/abs(fpeff))
        end do
        if (relmax > 0.03d0) form_ok = .false.
        if (form_ok) then
            print *, '  PASS: drag is Spitzer form, pc/v0 =', xc, ' max reldev', relmax
        else
            print *, '  FAIL: drag not Spitzer form, pc/v0 =', xc, ' max reldev', relmax
            passed = .false.; nfail = nfail + 1
        end if
    end subroutine test_slowing_down_distribution

    ! End-to-end stochastic check: an ensemble of alphas under the full collision
    ! operator (pitch scattering, energy scattering, drag). With a fixed RNG seed
    ! the flat-profile mode and the constant mode must give identical thermalized
    ! fraction and identical final (momentum, pitch) coordinates.
    subroutine test_loss_fraction_flat_vs_scalar(passed, nfail)
        logical, intent(inout) :: passed
        integer, intent(inout) :: nfail

        integer, parameter :: npart = 64, nstep = 400
        real(dp), parameter :: p_therm = 0.1d0
        real(dp) :: am1, am2, Z1, Z2, ealpha, v0
        real(dp) :: densi1, densi2, tempi1, tempi2, tempe
        real(dp) :: dchichi, slowrate, dchichi_norm, slowrate_norm
        real(dp) :: dpp, dhh, fpeff, dtauc
        real(dp) :: zend_s(2, npart), zend_f(2, npart)
        integer :: nlost_s, nlost_f, ip
        logical :: coords_ok

        print *, 'Testing thermalized fraction and end coordinates (flat vs scalar)...'

        am1 = 2.0d0; am2 = 3.0d0; Z1 = 1.0d0; Z2 = 1.0d0; ealpha = 3.5d6
        ni1_scale = 0.5d20; ni2_scale = 0.5d20
        Te_scale = 1.0d4; Ti1_scale = 1.0d4; Ti2_scale = 1.0d4
        densi1 = ni1_scale*1.0d-6; densi2 = ni2_scale*1.0d-6
        tempi1 = Ti1_scale; tempi2 = Ti2_scale; tempe = Te_scale

        call loacol_alpha(am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, &
                          ealpha, v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
        call set_flat_two_power()

        call coleff(1.0d0, dpp, dhh, fpeff)
        dtauc = 2.0d-3/abs(fpeff)

        call fill_constant_grid()
        call run_ensemble(npart, nstep, dtauc, p_therm, nlost_s, zend_s)

        call init_collision_profiles(am1, am2, Z1, Z2, ealpha, v0)
        call run_ensemble(npart, nstep, dtauc, p_therm, nlost_f, zend_f)

        coords_ok = .true.
        do ip = 1, npart
            if (zend_f(1, ip) /= zend_s(1, ip) .or. zend_f(2, ip) /= zend_s(2, ip)) &
                coords_ok = .false.
        end do

        if (nlost_f == nlost_s .and. coords_ok) then
            print *, '  PASS: thermalized', nlost_s, 'of', npart, &
                'and end coordinates identical'
        else
            if (nlost_f /= nlost_s) print *, '  FAIL: thermalized count differs', &
                nlost_s, nlost_f
            if (.not. coords_ok) print *, '  FAIL: end coordinates differ'
            passed = .false.; nfail = nfail + 1
        end if
    end subroutine test_loss_fraction_flat_vs_scalar

    ! Acceptance test (issue #363): the background Maxwellian is the fixed point
    ! of the collision operator. An ensemble sampled from the Maxwellian at the
    ! background temperature, evolved under energy scattering and drag, keeps its
    ! mean energy and its shape. enrat = E_alpha/T sets the equilibrium width:
    ! f(p) ~ p^2 exp(-enrat p^2), <p^2> = 3/(2 enrat).
    subroutine test_maxwellian_fixed_point(passed, nfail)
        logical, intent(inout) :: passed
        integer, intent(inout) :: nfail

        integer, parameter :: npart = 8000, nstep = 4000, nbin = 12
        real(dp) :: am1, am2, Z1, Z2, ealpha, v0
        real(dp) :: densi1, densi2, tempi1, tempi2, tempe
        real(dp) :: dchichi, slowrate, dchichi_norm, slowrate_norm
        real(dp) :: enr, sig, dpp, dhh, fpeff, dtauc, twopi
        real(dp) :: p0(npart), z(5)
        real(dp) :: h_init(nbin), h_final(nbin), h_ana(nbin)
        real(dp) :: plo, phi, db, pc, p2_init, p2_final, p2_eq
        real(dp) :: u1, u2, u3, u4, vx, vy, vz, reldev_shape, reldev_ana, ratio
        integer :: ip, it, ib, ierr
        logical :: moment_ok, shape_ok

        print *, 'Testing Maxwellian fixed point of the collision operator...'

        am1 = 2.0d0; am2 = 3.0d0; Z1 = 1.0d0; Z2 = 1.0d0; ealpha = 3.5d6
        tempe = 1.0d4; tempi1 = 1.0d4; tempi2 = 1.0d4
        ni1_scale = 0.5d20; ni2_scale = 0.5d20
        Te_scale = tempe; Ti1_scale = tempi1; Ti2_scale = tempi2
        densi1 = ni1_scale*1.0d-6; densi2 = ni2_scale*1.0d-6

        call loacol_alpha(am1, am2, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, &
                          ealpha, v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
        call set_flat_two_power()
        call init_collision_profiles(am1, am2, Z1, Z2, ealpha, v0)

        enr = enrat(1)
        sig = 1.0d0/sqrt(2.0d0*enr)
        p2_eq = 1.5d0/enr
        twopi = 8.0d0*atan(1.0d0)
        call coleff(1.0d0/sqrt(enr), dpp, dhh, fpeff)
        dtauc = 0.01d0/enr/max(dpp, abs(fpeff))

        plo = 0.0d0; phi = 4.0d0/sqrt(enr); db = (phi - plo)/dble(nbin)
        h_init = 0.0d0; h_final = 0.0d0

        call seed_rng()
        ! sample Maxwellian momenta p = |v|, each Cartesian component ~ N(0, sig^2)
        p2_init = 0.0d0
        do ip = 1, npart
            call random_number(u1); call random_number(u2)
            call random_number(u3); call random_number(u4)
            vx = sig*sqrt(-2.0d0*log(u1 + 1.0d-30))*cos(twopi*u2)
            vy = sig*sqrt(-2.0d0*log(u1 + 1.0d-30))*sin(twopi*u2)
            vz = sig*sqrt(-2.0d0*log(u3 + 1.0d-30))*cos(twopi*u4)
            p0(ip) = sqrt(vx*vx + vy*vy + vz*vz)
            p2_init = p2_init + p0(ip)**2
            ib = int((p0(ip) - plo)/db) + 1
            if (ib >= 1 .and. ib <= nbin) h_init(ib) = h_init(ib) + 1.0d0
        end do

        ! evolve with energy scattering + drag (iswmode = 2)
        p2_final = 0.0d0
        do ip = 1, npart
            z = 0.0d0; z(1) = 0.5d0; z(4) = p0(ip); z(5) = 0.0d0
            do it = 1, nstep
                call stost(z, dtauc, 2, ierr)
            end do
            p2_final = p2_final + z(4)**2
            ib = int((z(4) - plo)/db) + 1
            if (ib >= 1 .and. ib <= nbin) h_final(ib) = h_final(ib) + 1.0d0
        end do

        do ib = 1, nbin
            pc = plo + (dble(ib) - 0.5d0)*db
            h_ana(ib) = pc*pc*exp(-enr*pc*pc)
        end do
        call normalize(h_init, nbin)
        call normalize(h_final, nbin)
        call normalize(h_ana, nbin)

        ratio = (p2_final/dble(npart))/p2_eq
        moment_ok = (ratio > 0.95d0 .and. ratio < 1.05d0)
        if (moment_ok) then
            print *, '  PASS: mean energy preserved, <p^2> final/eq =', ratio
        else
            print *, '  FAIL: mean energy drifted, <p^2> final/eq =', ratio
            passed = .false.; nfail = nfail + 1
        end if

        reldev_shape = 0.0d0; reldev_ana = 0.0d0
        do ib = 1, nbin
            if (h_init(ib) > 0.05d0) then
                reldev_shape = max(reldev_shape, &
                                   abs(h_final(ib) - h_init(ib))/h_init(ib))
                reldev_ana = max(reldev_ana, abs(h_final(ib) - h_ana(ib))/h_ana(ib))
            end if
        end do
        shape_ok = (reldev_shape < 0.1d0 .and. reldev_ana < 0.1d0)
        if (shape_ok) then
            print *, '  PASS: distribution stays Maxwellian, max reldev', &
                max(reldev_shape, reldev_ana)
        else
            print *, '  FAIL: distribution distorted, reldev init/ana', &
                reldev_shape, reldev_ana
            passed = .false.; nfail = nfail + 1
        end if
    end subroutine test_maxwellian_fixed_point

    subroutine set_flat_two_power()
        profile_type = "two_power"
        active_profile = TWO_POWER
        Te_p1 = 1.0d0; Te_p2 = 0.0d0
        Ti1_p1 = 1.0d0; Ti1_p2 = 0.0d0
        Ti2_p1 = 1.0d0; Ti2_p2 = 0.0d0
        ni1_p1 = 1.0d0; ni1_p2 = 0.0d0
        ni2_p1 = 1.0d0; ni2_p2 = 0.0d0
    end subroutine set_flat_two_power

    subroutine fill_constant_grid()
        integer :: i
        do i = 1, N_S_GRID
            efcolf_grid(:, i) = efcolf
            velrat_grid(:, i) = velrat
            enrat_grid(:, i) = enrat
        end do
    end subroutine fill_constant_grid

    subroutine drag_residence(dtauc, plo, phi, nbin, pstop, hist)
        real(dp), intent(in) :: dtauc, plo, phi, pstop
        integer, intent(in) :: nbin
        real(dp), intent(out) :: hist(nbin)

        real(dp) :: z(5), p, dbin
        integer :: it, ib, ierr
        integer, parameter :: maxstep = 2000000

        dbin = (phi - plo)/dble(nbin)
        hist = 0.0d0
        z = 0.0d0; z(1) = 0.5d0; z(4) = 1.0d0; z(5) = 0.0d0
        do it = 1, maxstep
            call stost(z, dtauc, 3, ierr)
            p = z(4)
            if (p < pstop) exit
            if (p >= plo .and. p < phi) then
                ib = int((p - plo)/dbin) + 1
                if (ib >= 1 .and. ib <= nbin) hist(ib) = hist(ib) + 1.0d0
            end if
        end do
    end subroutine drag_residence

    subroutine run_ensemble(npart, nstep, dtauc, p_therm, nlost, zend)
        integer, intent(in) :: npart, nstep
        real(dp), intent(in) :: dtauc, p_therm
        integer, intent(out) :: nlost
        real(dp), intent(out) :: zend(2, npart)

        real(dp) :: z(5)
        integer :: ip, it, ierr

        call seed_rng()
        nlost = 0
        do ip = 1, npart
            z = 0.0d0; z(1) = 0.5d0
            z(4) = 0.2d0 + 0.8d0*(dble(ip) - 0.5d0)/dble(npart)
            z(5) = -0.99d0 + 1.98d0*(dble(ip) - 0.5d0)/dble(npart)
            do it = 1, nstep
                call stost(z, dtauc, 1, ierr)
                if (z(4) < p_therm) exit
            end do
            zend(1, ip) = z(4); zend(2, ip) = z(5)
            if (z(4) < p_therm) nlost = nlost + 1
        end do
    end subroutine run_ensemble

    subroutine seed_rng()
        integer :: n, i
        integer, allocatable :: seed(:)
        call random_seed(size=n)
        allocate (seed(n))
        do i = 1, n
            seed(i) = 4242 + i*7
        end do
        call random_seed(put=seed)
    end subroutine seed_rng

    subroutine normalize(h, n)
        integer, intent(in) :: n
        real(dp), intent(inout) :: h(n)
        real(dp) :: s
        s = sum(h)
        if (s > 0.0d0) h = h/s
    end subroutine normalize

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
