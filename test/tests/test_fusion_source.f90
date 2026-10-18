module fusion_source_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
contains
    subroutine analytic_magfie(x, bmod, jac, bder, hcov, hctr, hcurl)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, jac, bder(3), hcov(3), hctr(3), hcurl(3)
        ! An analytic, left-handed toroidal volume chart in cm^3.
        jac = -2.0e6_dp*(1.0_dp + x(1))*(1.0_dp + 0.4_dp*cos(x(2)))* &
            (1.0_dp + 0.3_dp*cos(x(3)))
        bmod = 1.0_dp
        bder = 0.0_dp
        hcov = 0.0_dp
        hctr = 0.0_dp
        hcurl = 0.0_dp
    end subroutine analytic_magfie
end module fusion_source_fixture

program test_fusion_source
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fusion_source, only: dt_reactivity, build_birth_cdf, draw_birth, &
        alpha_energy_joule
    use fusion_source_fixture, only: analytic_magfie
    use simple_profiles
    use magfie_sub, only: magfie
    implicit none
    integer, parameter :: n = 200000, nb = 20
    integer :: grid(3), counts(nb), i, j, b, unit, nseed
    integer, allocatable :: seed(:)
    real(dp), allocatable :: cdf(:)
    real(dp) :: power, z(5), expected(nb), a, h, chi2, moments(4)
    real(dp) :: reference, s, f, factor, power_coarse
    character(32) :: mode

    call get_command_argument(1, mode)
    call reset_profiles
    magfie => analytic_magfie
    grid = [200, 32, 32]
    select case (trim(mode))
    case ('negative_density')
        ni1_scale = -1.0_dp
    case ('zero_source')
        ni1_scale = 0.0_dp
    case ('invalid_bounds')
        call build_birth_cdf(0.5_dp, 0.2_dp, grid, cdf, power)
        stop ! WILL_FAIL test must fail if invalid inputs return normally.
    case ('invalid_grid')
        grid(1) = 0
    case ('hot_profile')
        Ti1_scale = 2.0e5_dp
        Ti2_scale = 2.0e5_dp
    case ('nonfinite')
        block
            use finite_check, only: quiet_nan
            ni1_scale = quiet_nan()
        end block
    case ('peaked')
        profile_type = 'two_power'
        ni1_p2 = 1.0_dp
        ni2_p2 = 1.0_dp
        call select_active_profile
    case ('power')
        Ti1_coef(2) = -0.8_dp
        Ti2_coef(2) = -0.3_dp
    end select

    ! Values from an independent evaluation of the published Table VII fit.
    if (abs(dt_reactivity(1.0e3_dp)/6.85688424959328e-27_dp - 1) > 1e-12_dp) &
        error stop '1 keV reactivity'
    if (abs(dt_reactivity(1.0e4_dp)/1.1361654705836232e-22_dp - 1) > 1e-12_dp) &
        error stop '10 keV reactivity'
    if (abs(dt_reactivity(1.0e5_dp)/8.447662458665835e-22_dp - 1) > 1e-12_dp) &
        error stop '100 keV reactivity'
    if (dt_reactivity(0.0_dp) /= 0.0_dp) error stop 'cold reactivity'
    call build_birth_cdf(0.0_dp, 1.0_dp, grid, cdf, power)
    select case (trim(mode))
    case ('negative_density', 'zero_source', 'invalid_grid', 'hot_profile', 'nonfinite')
        stop ! Returning normally makes the registered WILL_FAIL test fail.
    end select
    if (trim(mode) == 'power') then
        ! Independent composite Simpson quadrature (10000 radial intervals);
        ! angular integral is analytic, not the sampler's cell sum.
        reference = 0.0_dp
        do i = 0, 10000
            s = real(i, dp)/10000
            factor = 2.0_dp
            if (mod(i, 2) == 1) factor = 4.0_dp
            if (i == 0 .or. i == 10000) factor = 1.0_dp
            f = ni1_scale*ni2_scale*dt_reactivity(1.0e4_dp*(1 - 0.6_dp*s))
            reference = reference + factor*f*(1 + s)
        end do
        reference = reference/30000*(2*acos(-1.0_dp))**2*2*alpha_energy_joule
        print *, 'power, independent Simpson [W], relative error:', &
            power, reference, abs(power/reference - 1)
        if (abs(power/reference - 1) > 3e-5_dp) error stop 'power quadrature'
        power_coarse = power
        grid(1) = 400
        call build_birth_cdf(0.0_dp, 1.0_dp, grid, cdf, power)
        if (abs(power - reference) > 0.3_dp*abs(power_coarse - reference)) &
            error stop 'midpoint power does not converge quadratically'
        print *, 'PASS power quadrature and grid convergence'
        stop
    end if

    call random_seed(size=nseed)
    allocate (seed(nseed))
    seed = 543
    call random_seed(put=seed)
    counts = 0
    moments = 0.0_dp
    do i = 1, n
        call draw_birth(0.0_dp, 1.0_dp, grid, cdf, z)
        if (minval(z(1:3)) < 0 .or. z(1) >= 1) error stop 'position bounds'
        if (maxval(z(2:3)) >= 2*acos(-1.0_dp)) error stop 'angle bounds'
        if (z(4) /= 1 .or. abs(z(5)) > 1) error stop 'velocity bounds'
        b = min(nb, int(z(1)*nb) + 1)
        counts(b) = counts(b) + 1
        moments = moments + [cos(z(2)), cos(z(3)), z(5), z(5)**2]
    end do
    moments = moments/n
    if (maxval(abs(moments - [0.2_dp, 0.15_dp, 0.0_dp, 1.0_dp/3])) > &
        0.008_dp) error stop 'Jacobian angles or isotropic pitch'
    open (newunit=unit, file=trim(mode)//'_radial.csv', status='replace')
    write (unit, '(a)') 's_left,s_right,count,analytic_probability'
    do j = 1, nb
        a = real(j - 1, dp)/nb
        h = real(j, dp)/nb
        if (trim(mode) == 'peaked') then
            expected(j) = (primitive(h) - primitive(a))/primitive(1.0_dp)
        else
            expected(j) = (h + h*h/2 - a - a*a/2)/1.5_dp
        end if
        write (unit, '(2(f8.4,","),i8,",",es24.16)') a, h, counts(j), expected(j)
    end do
    close (unit)
    chi2 = sum((counts - n*expected)**2/(n*expected))
    ! 19 degrees of freedom; upper tail probability at 65 is < 1e-6.
    print *, trim(mode), ': chi2(19) = ', chi2, ' limit = 65'
    if (chi2 > 65.0_dp) error stop 'radial chi-squared'
    reference = 1.5_dp
    if (trim(mode) == 'peaked') reference = primitive(1.0_dp)
    reference = reference*2*(2*acos(-1.0_dp))**2*ni1_scale*ni2_scale* &
        dt_reactivity(1.0e4_dp)*alpha_energy_joule
    if (abs(power/reference - 1) > 2e-5_dp) error stop 'analytic total power'
    print *, 'PASS radial histogram, angular Jacobian, pitch, energy, total power'
contains
    pure function primitive(s) result(v)
        real(dp), intent(in) :: s
        real(dp) :: v
        ! Integral of (1+s)*(1-s)^2: peaked density profiles, flat T.
        v = s - s**2/2 - s**3/3 + s**4/4
    end function primitive
end program test_fusion_source
