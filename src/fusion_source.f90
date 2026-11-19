module fusion_source
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use finite_check, only: is_finite
    implicit none
    private
    public :: dt_reactivity, build_birth_cdf, draw_birth, alpha_energy_joule, &
             dd_neutron_cross_section, dd_neutron_reactivity, dd_beam_maxwellian_rate, &
              init_hermite

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: alpha_energy_joule = 3.5e6_dp*1.602176634e-19_dp
    integer, parameter :: nherm = 32
    real(dp) :: hermite_x(nherm) = 0.0_dp, hermite_w(nherm) = 0.0_dp
    logical :: hermite_ready = .false.

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

    function dd_neutron_cross_section(e_cm_kev) result(sigma)
        ! H.-S. Bosch & G.M. Hale, Nuclear Fusion 32 (1992) 611, Eqs. (8), (9),
        ! Table IV, D(d,n)3He, valid for 0.5 <= E_cm <= 4900 keV. Returns m^2.
        real(dp), intent(in) :: e_cm_kev
        real(dp) :: sigma, e, s
        real(dp), parameter :: bg = 31.3970_dp
        real(dp), parameter :: a(5) = [5.3701e4_dp, 3.3027e2_dp, -1.2706e-1_dp, &
                                       2.9327e-5_dp, -2.5151e-9_dp]

        if (.not. is_finite(e_cm_kev)) error stop 'dd cross section: nonfinite energy'
        if (e_cm_kev > 4900.0_dp) error stop 'dd cross section: energy above 4900 keV'
        sigma = 0.0_dp
        if (e_cm_kev < 0.5_dp) return
        e = e_cm_kev
        s = a(1) + e*(a(2) + e*(a(3) + e*(a(4) + e*a(5))))
        ! mb to m^2
        sigma = s/(e*exp(bg/sqrt(e)))*1.0e-31_dp
    end function dd_neutron_cross_section

    function dd_neutron_reactivity(temperature_ev) result(rate)
        ! Bosch & Hale (1992) Eq. (12), Table VII, D(d,n)3He, valid for
        ! 0.2-100 keV. Input eV, returns m^3/s.
        real(dp), intent(in) :: temperature_ev
        real(dp) :: rate, t, theta, xi
        real(dp), parameter :: bg = 31.3970_dp, mrc2 = 937814.0_dp
        real(dp), parameter :: c1 = 5.43360e-12_dp, c2 = 5.85778e-3_dp, &
                               c3 = 7.68222e-3_dp, c5 = -2.964e-6_dp

        if (.not. is_finite(temperature_ev)) &
            error stop 'dd reactivity: nonfinite ion temperature'
        if (temperature_ev < 0.0_dp .or. temperature_ev > 1.0e5_dp) &
            error stop 'dd reactivity: ion temperature outside [0,100000] eV'
        rate = 0.0_dp
        t = temperature_ev*1.0e-3_dp
        if (t < 1.0e-6_dp) return
        theta = t/(1.0_dp - t*c2/(1.0_dp + t*(c3 + t*c5)))
        xi = (bg**2/(4.0_dp*theta))**(1.0_dp/3.0_dp)
        rate = 1.0e-6_dp*c1*theta*sqrt(xi/(mrc2*t**3))*exp(-3.0_dp*xi)
    end function dd_neutron_reactivity

    function dd_beam_maxwellian_rate(speed, temperature_ev, mass_kg) result(rate)
        ! <sigma u> [m^3/s] of a deuteron of speed v [m/s] on Maxwellian deuterons
        ! of temperature T [eV] and mass m [kg], relative speed u, E_cm = m u^2/4,
        ! target thermal speed w = sqrt(2 T/m):
        ! <sigma u> = (sqrt(pi) v w)^-1 int_R sgn(u) u^2 sigma(|u|) exp(-(u - v)^2/w^2) du.
        real(dp), intent(in) :: speed, temperature_ev, mass_kg
        real(dp) :: rate, w, lo, hi, h, u, f, weight
        real(dp), parameter :: ev = 1.602176634e-19_dp, span = 10.0_dp, emax = 4900.0_dp
        integer, parameter :: n = 400
        integer :: k

        if (.not. all(is_finite([speed, temperature_ev, mass_kg])) .or. speed < 0.0_dp &
            .or. temperature_ev < 0.0_dp .or. mass_kg <= 0.0_dp) &
            error stop 'dd beam rate: invalid speed, temperature or mass'
        w = sqrt(2.0_dp*temperature_ev*ev/mass_kg)
        if (w <= 1.0e-6_dp*speed) then
            rate = speed*dd_neutron_cross_section(e_cm(speed))
            return
        end if
        rate = 0.0_dp
        if (speed >= 2.0_dp*w) then
            ! Gauss-Hermite in y = (u - v)/w; the integrand is smooth here.
            if (.not. hermite_ready) call init_hermite
            do k = 1, nherm
                u = speed + w*hermite_x(k)
                if (u > 0.0_dp .and. e_cm(u) <= emax) rate = rate + hermite_w(k)*u**2* &
                                                       dd_neutron_cross_section(e_cm(u))
            end do
            rate = rate/(sqrt(pi)*speed)
            return
        end if
        lo = 0.0_dp
        hi = speed + span*w
        ! The cross-section fit ends at emax, the Maxwellian tail there is negligible.
        hi = min(hi, sqrt(4.0e3_dp*emax*ev/mass_kg))
        h = (hi - lo)/n
        do k = 1, n
            u = lo + k*h
            weight = merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(k, 2) == 1), k == n)
            if (speed < 1.0e-3_dp*w) then
                ! Limit v -> 0 of the bracket divided by v.
                f = 4.0_dp*u**3*exp(-(u/w)**2)/w**2
            else
                f = u**2*(exp(-((u - speed)/w)**2) - exp(-((u + speed)/w)**2))/speed
            end if
            rate = rate + weight*f*dd_neutron_cross_section(e_cm(u))
        end do
        rate = rate*h/3.0_dp/(sqrt(pi)*w)
    contains
        pure function e_cm(v) result(e)
            real(dp), intent(in) :: v
            real(dp) :: e
            e = 0.25_dp*mass_kg*v**2/ev*1.0e-3_dp
        end function e_cm
    end function dd_beam_maxwellian_rate

    subroutine init_hermite
        ! Nodes and weights for int exp(-y^2) f(y) dy by Newton iteration on the
        ! orthonormal Hermite recurrence (Press et al., Numerical Recipes, gauher).
        real(dp) :: x(nherm), wt(nherm), z, z1, p1, p2, p3, pp
        integer :: i, j, its, m

        m = (nherm + 1)/2
        z = 0.0_dp
        do i = 1, m
            select case (i)
            case (1)
    z = sqrt(real(2*nherm + 1, dp)) - 1.85575_dp*real(2*nherm + 1, dp)**(-1.0_dp/6.0_dp)
            case (2)
                z = z - 1.14_dp*real(nherm, dp)**0.426_dp/z
            case (3)
                z = 1.86_dp*z - 0.86_dp*x(1)
            case (4)
                z = 1.91_dp*z - 0.91_dp*x(2)
            case default
                z = 2.0_dp*z - x(i - 2)
            end select
            do its = 1, 100
                p1 = pi**(-0.25_dp)
                p2 = 0.0_dp
                do j = 1, nherm
                    p3 = p2
                    p2 = p1
                    p1 = z*sqrt(2.0_dp/j)*p2 - sqrt(real(j - 1, dp)/j)*p3
                end do
                pp = sqrt(2.0_dp*nherm)*p2
                z1 = z
                z = z1 - p1/pp
                if (abs(z - z1) <= 1.0e-14_dp) exit
            end do
            x(i) = z
            x(nherm + 1 - i) = -z
            wt(i) = 2.0_dp/(pp*pp)
            wt(nherm + 1 - i) = wt(i)
        end do
        !$omp critical (fusion_source_hermite)
        if (.not. hermite_ready) then
            hermite_x = x
            hermite_w = wt
            hermite_ready = .true.
        end if
        !$omp end critical (fusion_source_hermite)
    end subroutine init_hermite

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
        ! Full torus, not one field period. magfie Jacobian is in cm^3.
        cell_volume = (shi - slo)*(2.0_dp*pi)**2/real(n, dp)*1.0e-6_dp
        ! Sampling backends read initialized VMEC/Boozer spline tables. Midpoints
        ! have s > 0, so Boozer's negative-radius flag is never written here.
        !$omp parallel do default(none) schedule(static) &
        !$omp shared(slo, shi, grid, cdf, cell_volume, magfie) &
        !$omp private(ir, it, iz, k, x, te, td, tt, nd, nt, source, &
        !$omp jac, bmod, bder, hcov, hctr, hcurl)
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
                    k = ((ir - 1)*grid(2) + it - 1)*grid(3) + iz
                    cdf(k) = source*abs(jac)*cell_volume
                end do
            end do
        end do
        !$omp end parallel do
        ! Keep the prefix sum ordered: thread count cannot change cell masses
        ! or their accumulation order, including for fixed-seed sampling.
        do k = 1, n
            cdf(k) = cdf(k - 1) + cdf(k)
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
