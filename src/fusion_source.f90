module fusion_source
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use finite_check, only: is_finite
    implicit none
    private
    public :: dt_reactivity, build_birth_cdf, draw_birth, alpha_energy_joule

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: alpha_energy_joule = 3.5e6_dp*1.602176634e-19_dp

contains

    function dt_reactivity(temperature_ev) result(rate)
        ! H.-S. Bosch & G.M. Hale, Nuclear Fusion 32 (1992) 611-631,
        ! doi:10.1088/0029-5515/32/4/I07, Eq. (12), Table VII (D-T).
        ! Input eV; published fit uses keV and cm^3/s; return m^3/s.
        real(dp), intent(in) :: temperature_ev
        real(dp) :: rate, t, theta, xi
        real(dp), parameter :: bg = 34.3827_dp, mrc2 = 1124656.0_dp
        real(dp), parameter :: c(7) = [1.17302e-9_dp, 1.51361e-2_dp, &
            7.51886e-2_dp, 4.60643e-3_dp, 1.35000e-2_dp, &
            -1.06750e-4_dp, 1.36600e-5_dp]

        if (.not. is_finite(temperature_ev)) &
            error stop 'fusion source: nonfinite ion temperature'
        if (temperature_ev < 0.0_dp .or. temperature_ev > 1.0e5_dp) &
            error stop 'fusion source: ion temperature outside [0,100000] eV'
        rate = 0.0_dp
        t = temperature_ev*1.0e-3_dp
        ! Below this temperature exp(-3*xi) underflows to zero.
        if (t < 1.0e-6_dp) return
        theta = t/(1.0_dp - t*(c(2) + t*(c(4) + t*c(6)))/ &
            (1.0_dp + t*(c(3) + t*(c(5) + t*c(7)))))
        xi = (bg**2/(4.0_dp*theta))**(1.0_dp/3.0_dp)
        rate = 1.0e-6_dp*c(1)*theta*sqrt(xi/(mrc2*t**3))*exp(-3.0_dp*xi)
    end function dt_reactivity

    subroutine build_birth_cdf(slo, shi, grid, cdf, power)
        use simple_profiles, only: get_plasma_params
        use magfie_sub, only: magfie
        real(dp), intent(in) :: slo, shi
        integer, intent(in) :: grid(3)
        real(dp), allocatable, intent(out) :: cdf(:)
        real(dp), intent(out) :: power
        real(dp) :: x(3), te, td, tt, nd, nt, source, jac, bmod
        real(dp) :: bder(3), hcov(3), hctr(3), hcurl(3), cell_volume
        integer :: ir, it, iz, k, n

        if (.not. all(is_finite([slo, shi]))) &
            error stop 'fusion source: nonfinite radial bounds'
        if (slo < 0.0_dp .or. shi > 1.0_dp .or. shi <= slo) &
            error stop 'fusion source: require 0 <= s_inner < s_outer <= 1'
        if (any(grid < 1) .or. any(grid > 1024)) &
            error stop 'fusion source: grid dimensions must be in [1,1024]'
        if (product(int(grid, int64)) > 8388608_int64) &
            error stop 'fusion source: grid exceeds 8388608 cells'
        n = product(grid)
        allocate (cdf(0:n))
        cdf(0) = 0.0_dp
        k = 0
        ! Full torus, not one field period. magfie Jacobian is in cm^3.
        cell_volume = (shi - slo)*(2.0_dp*pi)**2/real(n, dp)*1.0e-6_dp
        do ir = 1, grid(1)
            x(1) = slo + (shi - slo)*(real(ir, dp) - 0.5_dp)/grid(1)
            call get_plasma_params(x(1), te, td, tt, nd, nt)
            if (.not. all(is_finite([td, tt, nd, nt]))) &
                error stop 'fusion source: nonfinite profiles'
            if (min(td, tt, nd, nt) < 0.0_dp) &
                error stop 'fusion source: negative profiles'
            ! Relative Maxwellian temperature: mu*(T_D/m_D + T_T/m_T).
            source = nd*nt*dt_reactivity((3.0_dp*td + 2.0_dp*tt)/5.0_dp)
            do it = 1, grid(2)
                x(2) = 2.0_dp*pi*(real(it, dp) - 0.5_dp)/grid(2)
                do iz = 1, grid(3)
                    x(3) = 2.0_dp*pi*(real(iz, dp) - 0.5_dp)/grid(3)
                    call magfie(x, bmod, jac, bder, hcov, hctr, hcurl)
                    if (.not. is_finite(jac)) &
                        error stop 'fusion source: nonfinite volume Jacobian'
                    k = k + 1
                    cdf(k) = cdf(k - 1) + source*abs(jac)*cell_volume
                end do
            end do
        end do
        if (.not. is_finite(cdf(n)) .or. cdf(n) <= 0.0_dp) &
            error stop 'fusion source: total birth rate must be finite and positive'
        power = cdf(n)*alpha_energy_joule
        cdf = cdf/cdf(n)
    end subroutine build_birth_cdf

    subroutine draw_birth(slo, shi, grid, cdf, z)
        real(dp), intent(in) :: slo, shi, cdf(0:)
        integer, intent(in) :: grid(3)
        real(dp), intent(out) :: z(5)
        real(dp) :: u(5)
        integer :: lo, hi, mid, ir, it, iz, cell

        call random_number(u)
        lo = 0
        hi = size(cdf) - 1
        do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (cdf(mid) <= u(1)) then
                lo = mid
            else
                hi = mid
            end if
        end do
        cell = hi - 1
        iz = mod(cell, grid(3))
        it = mod(cell/grid(3), grid(2))
        ir = cell/(grid(2)*grid(3))
        ! Continuous uniform positions inside the selected midpoint cell.
        z(1) = slo + (shi - slo)*(ir + u(2))/grid(1)
        z(2) = 2.0_dp*pi*(it + u(3))/grid(2)
        z(3) = 2.0_dp*pi*(iz + u(4))/grid(3)
        z(4) = 1.0_dp
        z(5) = 2.0_dp*u(5) - 1.0_dp
    end subroutine draw_birth

end module fusion_source
