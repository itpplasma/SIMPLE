program test_collision_coulomb_log
    ! Coulomb logarithms against an independent SI evaluation of
    ! ln(lambda_D/b_min), and zero Maxwellian flux of the speed operator.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use collis_alp, only: compute_collision_coeffs, loacol_alpha, coleff, &
                          log_factor, nsorts, enrat
    implicit none
    real(dp), parameter :: e = 1.602176634e-19_dp, eps0 = 8.8541878128e-12_dp
    real(dp), parameter :: hbar = 1.054571817e-34_dp, pi = acos(-1.0_dp)
    real(dp), parameter :: mp = 1.6726e-27_dp, me = 9.1094e-31_dp
    real(dp), parameter :: energy = 55000.0_dp, temp = 1000.0_dp
    real(dp), parameter :: n1 = 1.0e13_dp, n2 = 1.0e12_dp, z2 = 6.0_dp, a2 = 12.0_dp
    real(dp) :: coef(nsorts), velrat(nsorts), enr(nsorts), clog(nsorts), bratio(nsorts)
    real(dp) :: v0, ne, ref(nsorts), mass(nsorts), charge(nsorts), f, df
    real(dp) :: p, h, dpp, dhh, fpeff, dplus, dminus, expected, rates(4)
    real(dp), parameter :: speeds(3) = [0.05_dp, 0.4_dp, 1.5_dp]
    integer :: i, k

    call compute_collision_coeffs(1.0_dp, a2, 1.0_dp, z2, n1, n2, temp, temp, temp, &
                                  energy, v0, coef, velrat, enr, 1.0_dp, 1.0_dp, &
                                  clog, bratio)
    ne = n1 + z2*n2
    mass = [mp, a2*mp, me]
    charge = [1.0_dp, z2, 1.0_dp]
    do i = 1, nsorts
        ref(i) = si_log(i, 1.0_dp)
        if (abs(clog(i) - ref(i)) > 1.0e-5_dp*ref(i)) &
            error stop 'Coulomb logarithm differs from ln(lambda_D/b_min)'
    end do
    ! The coefficients carry the logarithms: main ion relative to electrons.
    expected = n1/ne*ref(1)/ref(3)*velrat(1)/velrat(3)
    if (abs(coef(1)/coef(3) - expected) > 1.0e-5_dp*expected) &
        error stop 'Collision coefficients do not carry the species logarithms'
    do k = 1, 2
        p = merge(0.3_dp, 2.0_dp, k == 1)
        do i = 1, nsorts
            call log_factor(i, p, velrat, f, df, clog, bratio)
            if (abs(f - si_log(i, p)/ref(i)) > 1.0e-5_dp) &
                error stop 'Speed dependence of the Coulomb logarithm is wrong'
        end do
    end do
    ! Equal temperatures: the Maxwellian p**2*exp(-enrat*p**2) has zero flux,
    ! so the Ito drift equals D' + D*(2/p - 2*enrat*p) for the speed-dependent D.
    call loacol_alpha(1.0_dp, a2, 1.0_dp, z2, n1, n2, temp, temp, temp, energy, v0, &
                      rates(1), rates(2), rates(3), rates(4), 1.0_dp, 1.0_dp)
    do k = 1, 3
        p = speeds(k)
        h = 1.0e-5_dp*p
        call coleff(p + h, dplus, dhh, fpeff)
        call coleff(p - h, dminus, dhh, fpeff)
        call coleff(p, dpp, dhh, fpeff)
        expected = (dplus - dminus)/(2*h) + dpp*(2/p - 2*enrat(1)*p)
        if (abs(fpeff - expected) > 1.0e-6_dp*abs(fpeff)) &
            error stop 'Speed operator does not preserve the Maxwellian'
    end do
    print *, 'PASS: species Coulomb logarithms and Maxwellian flux'

contains

    real(dp) function si_log(i, p)
        integer, intent(in) :: i
        real(dp), intent(in) :: p
        real(dp) :: debye, v2, mr, bcl, bqm, t
        t = temp*e
        debye = sqrt(eps0*t/(e**2*(n1*1e6 + z2**2*n2*1e6 + ne*1e6)))
        v2 = (p*v0*1e-2_dp)**2 + 2*t/mass(i)
        mr = mp*mass(i)/(mp + mass(i))
        bcl = charge(i)*e**2/(4*pi*eps0*mr*v2)
        bqm = hbar/(2*mr*sqrt(v2))
        si_log = log(debye/max(bcl, bqm))
    end function si_log

end program test_collision_coulomb_log
