program test_collision_slowing_down
    ! Mean energy of a fast-ion ensemble under the full Coulomb operator
    ! against the test-particle energy-loss rate d<E>/dt = -sum_b nu_eps^b E
    ! of the NRL Plasma Formulary (2019, p. 31), which is written in terms of
    ! the Maxwell integral psi(x) rather than the Chandrasekhar G(x) used by
    ! the code. The ODE is integrated independently in physical time.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use collis_alp, only: stost, loacol_alpha, efcolf, velrat, enrat, &
        efcolf_grid, velrat_grid, enrat_grid, N_S_GRID
    implicit none

    integer, parameter :: npart = 10000, ncheck = 3
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: pmass = 1.6726d-24, emass = 9.1094d-28
    real(dp), parameter :: e_cgs = 4.8032d-10, ev_cgs = 1.6022d-12
    ! 55 keV H+ in an H+/C6+ plasma at 2 keV, ne = 5.6e13 cm^-3 (Zeff = 1.6).
    real(dp), parameter :: am1 = 1.0_dp, am2 = 12.0_dp, z1 = 1.0_dp, z2 = 6.0_dp
    real(dp), parameter :: densi1 = 5.0d13, densi2 = 1.0d12
    real(dp), parameter :: tempi = 2.0d3, tempe = 2.0d3, energy0 = 5.5d4
    real(dp), parameter :: test_mass = 1.0_dp, test_charge = 1.0_dp
    real(dp) :: v0, dchichi, slowrate, dchichi_norm, slowrate_norm
    real(dp) :: z(5), mean(ncheck), meansq(ncheck), oracle(ncheck), sigma
    real(dp) :: times(ncheck), tau_prev, t_prev, energy
    integer :: i, k, ierr, nseed
    integer, allocatable :: seed(:)

    call loacol_alpha(am1, am2, z1, z2, densi1, densi2, tempi, tempi, tempe, &
        energy0, v0, dchichi, slowrate, dchichi_norm, &
        slowrate_norm, test_mass, test_charge)
    do i = 1, N_S_GRID
        efcolf_grid(:, i) = efcolf
        velrat_grid(:, i) = velrat
        enrat_grid(:, i) = enrat
    end do

    ! Checkpoints at 10, 30 and 60 % of the initial energy-loss time.
    times = [0.1_dp, 0.3_dp, 0.6_dp]/loss_rate(energy0)
    energy = energy0
    t_prev = 0.0_dp
    do k = 1, ncheck
        call integrate_oracle(energy, times(k) - t_prev)
        oracle(k) = energy/energy0
        t_prev = times(k)
    end do

    call random_seed(size=nseed)
    allocate (seed(nseed))
    seed = 20261028
    call random_seed(put=seed)
    mean = 0.0_dp
    meansq = 0.0_dp
    do i = 1, npart
        z = [0.5_dp, 0.0_dp, 0.0_dp, 1.0_dp, 2.0_dp*(i - 0.5_dp)/npart - 1.0_dp]
        tau_prev = 0.0_dp
        do k = 1, ncheck
            call stost(z, v0*times(k) - tau_prev, 1, ierr)
            if (ierr /= 0) error stop 'stost failed during slowing down'
            tau_prev = v0*times(k)
            mean(k) = mean(k) + z(4)**2/npart
            meansq(k) = meansq(k) + z(4)**4/npart
        end do
    end do

    do k = 1, ncheck
        sigma = sqrt((meansq(k) - mean(k)**2)/npart)
        write (*, '(a,es11.4,a,3es14.6)') 't=', times(k), &
            ' s  <E>/E0, oracle, sigma: ', mean(k), oracle(k), sigma
        if (abs(mean(k) - oracle(k)) > 3.0_dp*sigma) &
            error stop 'Mean energy disagrees with the NRL energy-loss rate'
    end do

contains

    real(dp) function loss_rate(e) result(nu)
        ! sum_b nu_eps^b at test-particle energy e (eV), in 1/s.
        real(dp), intent(in) :: e
        real(dp) :: m, v, dense, lam(3), mb(3), qb(3), nb(3), tb(3), x, nu0
        integer :: b

        m = test_mass*pmass
        v = sqrt(2.0_dp*e*ev_cgs/m)
        dense = z1*densi1 + z2*densi2
        mb = [am1*pmass, am2*pmass, emass]
        qb = [z1, z2, -1.0_dp]*e_cgs
        nb = [densi1, densi2, dense]
        tb = [tempi, tempi, tempe]
        ! Coulomb logarithms are inputs of the model; same expressions as the code.
        lam(1) = 23.0_dp - log(sqrt(densi1*z1**2/tempi)*test_charge*z1* &
            (test_mass + am1)/(test_mass*tempi + am1*energy0))
        lam(2) = 23.0_dp - log(sqrt(densi2*z2**2/tempi)*test_charge*z2* &
            (test_mass + am2)/(test_mass*tempi + am2*energy0))
        lam(3) = 24.0_dp - log(sqrt(dense)/tempe)
        nu = 0.0_dp
        do b = 1, 3
            x = mb(b)*v**2/(2.0_dp*tb(b)*ev_cgs)
            nu0 = 4.0_dp*pi*(test_charge*e_cgs)**2*qb(b)**2*lam(b)*nb(b)/(m**2*v**3)
            nu = nu + 2.0_dp*(m/mb(b)*psi(x) - dpsi(x))*nu0
        end do
    end function loss_rate

    real(dp) function psi(x)
        real(dp), intent(in) :: x
        psi = erf(sqrt(x)) - 2.0_dp/sqrt(pi)*sqrt(x)*exp(-x)
    end function psi

    real(dp) function dpsi(x)
        real(dp), intent(in) :: x
        dpsi = 2.0_dp/sqrt(pi)*sqrt(x)*exp(-x)
    end function dpsi

    subroutine integrate_oracle(e, duration)
        ! Classical RK4 for dE/dt = -nu_eps(E) E; 4000 steps per interval.
        real(dp), intent(inout) :: e
        real(dp), intent(in) :: duration
        integer, parameter :: nstep = 4000
        real(dp) :: h, k1, k2, k3, k4
        integer :: n

        h = duration/nstep
        do n = 1, nstep
            k1 = -loss_rate(e)*e
            k2 = -loss_rate(e + 0.5_dp*h*k1)*(e + 0.5_dp*h*k1)
            k3 = -loss_rate(e + 0.5_dp*h*k2)*(e + 0.5_dp*h*k2)
            k4 = -loss_rate(e + h*k3)*(e + h*k3)
            e = e + h*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
        end do
    end subroutine integrate_oracle
end program test_collision_slowing_down
