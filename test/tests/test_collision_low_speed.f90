!> Low-speed behaviour of collis_alp::stost in a cold plasma. Below the thermal
!> speed of the slowest field species the radial speed process is Bessel-like,
!> and the Euler walk in p can collapse to p ~ 1e-10, where it cannot represent
!> further time progress (ierr = 4). The Cartesian substep is regular there.
program test_collision_low_speed
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use collis_alp, only: stost, coleff, efcolf, velrat, enrat, efcolf_grid, &
                          velrat_grid, enrat_grid, N_S_GRID, &
                          compute_collision_coeffs
    implicit none

    real(dp), parameter :: ealpha = 55.0e3_dp, temp = 20.0_dp
    real(dp) :: p_switch, rate, v0
    integer :: nseed
    integer, allocatable :: seed(:)

    call random_seed(size=nseed)
    allocate (seed(nseed))
    call setup_cold_plasma()
    p_switch = 1.0_dp/maxval(velrat)
    rate = sum(efcolf*enrat)*2.0_dp*0.75225278_dp
    print '(a,es10.3,a,es10.3)', 'p_switch = ', p_switch, ', rate = ', rate

    call check_tiny_speed_start()
    call check_generator_moments()
    call check_maxwellian()
    print *, 'PASS: low-speed collision substep'
contains

    subroutine setup_cold_plasma()
        ! H test ion in a 20 eV H plasma with 2.5% C6+, n_e = 2e19 m^-3.
        real(dp) :: ni1, ni2
        integer :: node
        ni1 = 2.0e13_dp/(1.0_dp + 6.0_dp*0.025_dp)
        ni2 = 0.025_dp*ni1
        call compute_collision_coeffs(1.0_dp, 12.0_dp, 1.0_dp, 6.0_dp, ni1, ni2, &
                                      temp, temp, temp, ealpha, v0, &
                                      efcolf, velrat, enrat)
        do node = 1, N_S_GRID
            efcolf_grid(:, node) = efcolf
            velrat_grid(:, node) = velrat
            enrat_grid(:, node) = enrat
        end do
    end subroutine setup_cold_plasma

    subroutine check_tiny_speed_start()
        ! Before the Cartesian substep, a start at p = 1e-9 or 1e-12 ended
        ! immediately with ierr = 4 (substep ~ p**2 below ulp of the interval).
        real(dp), parameter :: starts(2) = [1.0e-9_dp, 1.0e-12_dp]
        integer :: k, mode, ierr
        real(dp) :: z(5)
        do mode = 1, 2
            do k = 1, 2
                z = [0.5_dp, 0.2_dp, 0.3_dp, starts(k), 0.4_dp]
                call stost(z, 6.76_dp, mode, ierr)
                if (ierr /= 0) error stop 'Tiny-speed start returned an error'
                if (.not. all(z == z) .or. z(4) <= 0.0_dp .or. &
                    abs(z(5)) > 1.0_dp) error stop 'Tiny-speed start left bad state'
            end do
        end do
    end subroutine check_tiny_speed_start

    subroutine check_generator_moments()
        ! One substep from either side of p_switch: first moments of p and
        ! lambda must match fpeff*h and -2*dhh*lambda*h.
        integer, parameter :: n = 400000
        real(dp), parameter :: lam = 0.5_dp
        real(dp) :: fac(3) = [0.5_dp, 0.98_dp, 1.02_dp]
        real(dp) :: p, h, dpp, dhh, fpeff, z(5), dp1, dl1, m1, m2, s1, s2
        integer :: k, i, ierr
        seed = 77
        call random_seed(put=seed)
        do k = 1, size(fac)
            p = fac(k)*p_switch
            if (k == 3) p = 1.02_dp*p_switch
            call coleff(p, dpp, dhh, fpeff)
            h = 0.5_dp*min(0.01_dp/dhh, 0.1_dp*p/abs(fpeff), &
                           0.01_dp*p**2/(6.0_dp*dpp))
            m1 = 0.0_dp
            m2 = 0.0_dp
            s1 = 0.0_dp
            s2 = 0.0_dp
            do i = 1, n
                z = [0.5_dp, 0.2_dp, 0.3_dp, p, lam]
                call stost(z, h, 1, ierr)
                if (ierr /= 0) error stop 'Generator substep failed'
                dp1 = z(4) - p
                dl1 = z(5) - lam
                m1 = m1 + dp1
                s1 = s1 + dp1**2
                m2 = m2 + dl1
                s2 = s2 + dl1**2
            end do
            m1 = m1/n
            m2 = m2/n
            s1 = sqrt((s1/n - m1**2)/n)
            s2 = sqrt((s2/n - m2**2)/n)
            print '(a,f6.3,4es12.4)', 'moments p/p_switch, dp: ', p/p_switch, &
                m1, fpeff*h, m2, -2.0_dp*dhh*lam*h
            if (abs(m1 - fpeff*h) > 5.0_dp*s1) &
                error stop 'Speed drift disagrees with fpeff'
            if (abs(m2 + 2.0_dp*dhh*lam*h) > 5.0_dp*s2) &
                error stop 'Pitch drift disagrees with -2*dhh*lambda'
        end do
    end subroutine check_generator_moments

    subroutine check_maxwellian()
        integer, parameter :: n = 4096
        real(dp) :: z(5), sp2, sl, sl2, target, sig, duration
        integer :: i, ierr
        duration = 6.0_dp/rate
        seed = 20261023
        call random_seed(put=seed)
        sp2 = 0.0_dp
        sl = 0.0_dp
        sl2 = 0.0_dp
        do i = 1, n
            z = [0.5_dp, 0.2_dp, 0.3_dp, 1.0e-6_dp, 0.0_dp]
            call stost(z, duration, 1, ierr)
            if (ierr /= 0) error stop 'Thermalization ended in an error'
            sp2 = sp2 + z(4)**2/n
            sl = sl + z(5)/n
            sl2 = sl2 + z(5)**2/n
        end do
        target = 1.5_dp*temp/ealpha
        ! Maxwellian: var(E) = 1.5*T**2, var(lambda) = 1/3, var(lambda**2) = 4/45.
        sig = sqrt(1.5_dp)*temp/ealpha/sqrt(real(n, dp))
        print '(a,3es12.4)', '<p^2>, target, sigma: ', sp2, target, sig
        print '(a,2es12.4)', '<lambda>, <lambda^2>: ', sl, sl2
        if (abs(sp2 - target) > 5.0_dp*sig) error stop '<p^2> off Maxwellian'
        if (abs(sl) > 5.0_dp*sqrt(1.0_dp/3.0_dp/n)) error stop '<lambda> not 0'
        if (abs(sl2 - 1.0_dp/3.0_dp) > 5.0_dp*sqrt(4.0_dp/45.0_dp/n)) &
            error stop '<lambda^2> not 1/3'
        ! Subthermal starts never error, from many seeds.
        do i = 1, 2000
            seed = 1000 + 7*i
            call random_seed(put=seed)
            z = [0.5_dp, 0.2_dp, 0.3_dp, 0.5_dp*p_switch, 0.3_dp]
            call stost(z, 6.76_dp, 1, ierr)
            if (ierr /= 0) error stop 'Subthermal marker errored'
        end do
    end subroutine check_maxwellian
end program test_collision_low_speed
