!> Normalized drift and diffusion of collis_alp::coleff against the
!> test-particle rates of the NRL Plasma Formulary (2019, p. 31) in physical
!> units: slowing down nu_s, perpendicular diffusion nu_perp and parallel
!> diffusion nu_par, written with the Maxwell integral psi(x) instead of the
!> Chandrasekhar G(x) of the code. For speed v = p v0 and path tau = v0 t,
!>   dpp   = nu_par p^2/(2 v0),
!>   dhh   = nu_perp/(4 v0),
!>   fpeff = (nu_perp/2 - nu_s) p/v0,
!> the last from d<|v|>/dt = -nu_s v + nu_perp v/2.
program test_collision_rates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use collis_alp, only: coleff, efcolf, velrat, enrat, nsorts, &
                          compute_collision_coeffs, clog_ref, clog_bratio
    use util, only: pi, e_charge, p_mass, e_mass, ev
    implicit none

    ! 93 keV D in D and C6+ at Z_eff = 1.6, n_e = 5e13 cm^-3, T = 2 keV.
    real(dp), parameter :: ebeam = 93.0e3_dp, temp = 2.0e3_dp, ne = 5.0e13_dp
    real(dp), parameter :: nc = ne*0.6_dp/30.0_dp, nd = ne - 6.0_dp*nc
    real(dp), parameter :: mass(nsorts) = [2.0_dp*p_mass, 12.0_dp*p_mass, e_mass]
    real(dp), parameter :: charge(nsorts) = [1.0_dp, 6.0_dp, 1.0_dp]
    real(dp), parameter :: dens(nsorts) = [nd, nc, ne]
    real(dp), parameter :: ps(6) = [0.05_dp, 0.2_dp, 0.5_dp, 1.0_dp, 1.5_dp, 3.0_dp]
    real(dp) :: v0, clog(nsorts), bratio(nsorts), dpp, dhh, fpeff
    real(dp) :: nus, nuperp, nupar, err
    integer :: k

    call compute_collision_coeffs(2.0_dp, 12.0_dp, 1.0_dp, 6.0_dp, nd, nc, &
                                  temp, temp, temp, ebeam, v0, efcolf, velrat, enrat, &
                                  test_mass=2.0_dp, test_charge=1.0_dp, &
                                  clog_out=clog, bratio_out=bratio)
    ! Speed-independent logarithms, as in the formulary.
    clog_ref = 0.0_dp
    clog_bratio = 0.0_dp
    err = 0.0_dp
    do k = 1, size(ps)
        call coleff(ps(k), dpp, dhh, fpeff)
        call nrl_rates(ps(k)*v0, nus, nuperp, nupar)
        err = max(err, abs(dpp/(nupar*ps(k)**2/(2.0_dp*v0)) - 1.0_dp), &
                  abs(dhh/(nuperp/(4.0_dp*v0)) - 1.0_dp), &
                  abs(fpeff/((0.5_dp*nuperp - nus)*ps(k)/v0) - 1.0_dp))
    end do
    print '(a,es9.2)', 'max relative difference from the formulary ', err
    ! The series below v = 0.01 truncates at O(v^4).
    if (err > 1.0e-9_dp) error stop 'collision coefficients differ from the NRL rates'
    print '(a)', 'PASS'

contains

    subroutine nrl_rates(v, nus, nuperp, nupar)
        !> Sums over D, C6+ and electrons of nu_s, nu_perp, nu_par [1/s], cgs.
        real(dp), intent(in) :: v
        real(dp), intent(out) :: nus, nuperp, nupar
        real(dp) :: x, psi, dpsi, nu0
        integer :: b
        nus = 0.0_dp
        nuperp = 0.0_dp
        nupar = 0.0_dp
        do b = 1, nsorts
            x = mass(b)*v**2/(2.0_dp*temp*ev)
            psi = erf(sqrt(x)) - 2.0_dp*sqrt(x/pi)*exp(-x)
            dpsi = 2.0_dp*sqrt(x/pi)*exp(-x)
            nu0 = 4.0_dp*pi*e_charge**4*charge(b)**2*clog(b)*dens(b) &
                  /((2.0_dp*p_mass)**2*v**3)
            nus = nus + (1.0_dp + 2.0_dp*p_mass/mass(b))*psi*nu0
            nuperp = nuperp + 2.0_dp*((1.0_dp - 0.5_dp/x)*psi + dpsi)*nu0
            nupar = nupar + psi/x*nu0
        end do
    end subroutine nrl_rates
end program test_collision_rates
