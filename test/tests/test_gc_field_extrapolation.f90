module gc_extrapolation_oracle
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_mod, only: field_can_t
    implicit none
    real(dp), parameter :: MU0 = 0.6_dp
    logical :: axial = .false.
contains
    subroutine polynomial_field(f, r, theta, phi, mode)
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: r, theta, phi
        integer, intent(in) :: mode
        real(dp) :: c, s
        c = 1.0_dp
        s = 0.0_dp
        if (.not. axial) then
            c = cos(phi)
            s = sin(phi)
        end if
        ! Canonical Hamiltonian witness, not an equilibrium:
        ! p_theta = r + (0.1 sin(theta) + 0.2 p_phi) r^2,
        ! H = p_phi^2/2 + mu (1 + r*(0.2 cos(phi) + 0.15 sin(theta))).
        f%Ath = r + 0.1_dp*r*r*sin(theta)
        f%Aph = 0.0_dp
        f%hth = 0.2_dp*r*r
        f%hph = 1.0_dp
        f%Bmod = 1.0_dp + r*(0.2_dp*c + 0.15_dp*sin(theta))
        f%dAth = [1 + 0.2_dp*r*sin(theta), 0.1_dp*r*r*cos(theta), 0.0_dp]
        f%dAph = 0.0_dp
        f%dhth = [0.4_dp*r, 0.0_dp, 0.0_dp]
        f%dhph = 0.0_dp
        f%dBmod = [0.2_dp*c + 0.15_dp*sin(theta), &
            0.15_dp*r*cos(theta), -0.2_dp*r*s]
        if (mode <= 0) return
        f%d2Ath = [0.2_dp*sin(theta), 0.2_dp*r*cos(theta), 0.0_dp, &
            -0.1_dp*r*r*sin(theta), 0.0_dp, 0.0_dp]
        f%d2Aph = 0.0_dp
        f%d2hth = [0.4_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        f%d2hph = 0.0_dp
        f%d2Bmod = [0.0_dp, 0.15_dp*cos(theta), -0.2_dp*s, &
            -0.15_dp*r*sin(theta), 0.0_dp, &
            merge(0.0_dp, -0.2_dp*r*c, axial)]
    end subroutine polynomial_field

    pure function canonical_state(z) result(w)
        real(dp), intent(in) :: z(4)
        real(dp) :: w(4)
        w = [z(2), z(3), z(1) + &
            (0.1_dp*sin(z(2)) + 0.2_dp*z(4))*z(1)**2, z(4)]
    end function canonical_state

    function canonical_rhs(z) result(rhs)
        real(dp), intent(in) :: z(4)
        real(dp) :: rhs(4), d, g, c, s
        c = 1.0_dp
        s = 0.0_dp
        if (.not. axial) then
            c = cos(z(3))
            s = sin(z(3))
        end if
        d = 1 + 2*(0.1_dp*sin(z(2)) + 0.2_dp*z(4))*z(1)
        g = MU0*(0.2_dp*c + 0.15_dp*sin(z(2)))
        ! Closed derivatives of H after solving the quadratic p_theta(r),
        ! independent of the production field jets.
        rhs = [g/d, z(4) - g*0.2_dp*z(1)**2/d, &
            -MU0*0.15_dp*z(1)*cos(z(2)) + g*0.1_dp*z(1)**2*cos(z(2))/d, &
            MU0*0.2_dp*z(1)*s]
    end function canonical_rhs
end module gc_extrapolation_oracle

program test_gc_field_extrapolation
    use gc_extrapolation_oracle
    use field_can_mod, only: evaluate, field_can_t
    use orbit_symplectic, only: orbit_sympl_init, newton_midpoint, &
        attempt_midpoint_step, orbit_timestep_sympl_midpoint, &
        newton_rk_gauss, orbit_timestep_sympl_rk_gauss
    use orbit_symplectic_base, only: symplectic_integrator_t, &
        MIDPOINT, GAUSS1, coeff_rk_gauss, sympl_rmax
    use vector_potentail_mod, only: torflux
    implicit none
    real(dp), parameter :: HSTEP = 0.4_dp
    real(dp), parameter :: START(4) = [0.35_dp, 0.7_dp, 1.1_dp, 0.8_dp]
    integer :: s, model, failures
    evaluate => polynomial_field
    torflux = 1.0_dp
    failures = 0
    do model = 1, 2
        axial = model == 1
        do s = 1, 4
            call check_gauss(s)
        end do
    end do
    axial = .false.
    call check_midpoint(1.0e-3_dp, .false.)
    call check_midpoint(1.0e-5_dp, .false.)
    call check_midpoint(1.0e-3_dp, .true.)
    call check_midpoint(1.0e-5_dp, .true.)
    if (failures /= 0) error stop 'GC field extrapolation oracle failed'
contains
    subroutine initialize(si, f, tolerance, mode)
        type(symplectic_integrator_t), intent(out) :: si
        type(field_can_t), intent(out) :: f
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: mode
        f%mu = MU0
        f%ro0 = 1.0_dp
        si%pabs = 1.0_dp
        call orbit_sympl_init(si, f, START, HSTEP, 1, tolerance, mode)
    end subroutine initialize

    subroutine check_gauss(stages)
        integer, intent(in) :: stages
        type(symplectic_integrator_t) :: si, initial_si
        type(field_can_t) :: f, initial_f, fs(stages)
        real(dp) :: x(4*stages), previous(4*stages), a(stages, stages)
        real(dp) :: weights(stages), c(stages), expected(4), actual(4)
        real(dp) :: correction, error, bound
        integer :: k, status
        call initialize(si, f, 1.0e-4_dp, GAUSS1 + stages - 1)
        initial_si = si
        initial_f = f
        do k = 1, stages
            fs(k) = f
            x(4*k - 3:4*k) = si%z
        end do
        call newton_rk_gauss(si, fs, stages, x, si%atol, si%rtol, 32, &
            previous, status)
        if (status /= 0) error stop 'oracle Gauss stage solve'
        correction = maxval(abs(x - previous))
        call coeff_rk_gauss(stages, a, weights, c)
        expected = canonical_state(START)
        do k = 1, stages
            expected = expected + HSTEP*weights(k)*canonical_rhs(x(4*k - 3:4*k))
        end do
        si = initial_si
        f = initial_f
        call orbit_timestep_sympl_rk_gauss(si, f, stages, status)
        if (status /= 0) error stop 'oracle Gauss production step'
        actual = [si%z(2), si%z(3), f%pth, si%z(4)]
        error = maxval(abs(actual - expected))
        ! Quadratic Taylor remainder on this compact polynomial witness;
        ! the loose solve exposes a spurious linear term above roundoff.
        bound = 50*HSTEP*correction**2 + 200*epsilon(1.0_dp)
        print '(a,i0,a,l1,a,3es12.3)', 'Gauss ', stages, ' axial=', axial, &
            ' correction/error/bound ', correction, error, bound
        if (error > bound) failures = failures + 1
        if (axial .and. abs(si%z(4) - START(4)) > &
            200*epsilon(1.0_dp)) failures = failures + 1
    end subroutine check_gauss

    subroutine check_midpoint(tolerance, multiple_volume)
        real(dp), intent(in) :: tolerance
        logical, intent(in) :: multiple_volume
        type(symplectic_integrator_t) :: si, initial_si
        type(field_can_t) :: f, initial_f
        real(dp) :: x(5), previous(5), correction, expected(4), error, bound
        integer :: status
        sympl_rmax = merge(2.0_dp, 1.0_dp, multiple_volume)
        call initialize(si, f, tolerance, MIDPOINT)
        initial_si = si
        initial_f = f
        x = [si%z, si%z(1)]
        call newton_midpoint(si, f, x, si%atol, si%rtol, 32, previous, status)
        if (status /= 0) error stop 'oracle midpoint stage solve'
        correction = maxval(abs(x(1:4) - previous(1:4)))
        si = initial_si
        f = initial_f
        if (multiple_volume) then
            call orbit_timestep_sympl_midpoint(si, f, status)
        else
            call attempt_midpoint_step(si, f, status)
        end if
        if (status /= 0) error stop 'oracle midpoint production step'
        expected = canonical_state(si%z)
        error = abs(f%pth - expected(3))
        bound = 10*correction**2 + 200*epsilon(1.0_dp)
        print '(a,l1,a,4es12.3)', 'Midpoint multivolume=', multiple_volume, &
            ' tol/correction/error/bound ', tolerance, correction, error, bound
        if (error > bound) failures = failures + 1
        sympl_rmax = 1.0_dp
    end subroutine check_midpoint
end program test_gc_field_extrapolation
