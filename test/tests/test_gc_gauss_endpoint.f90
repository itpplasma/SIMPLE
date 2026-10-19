program test_gc_gauss_endpoint
    ! Gauss steps advance (theta, phi, p_theta, p_phi); the radius is defined
    ! implicitly by p_theta(r, theta, phi, p_phi). The step reports a Hermite
    ! extrapolation of the stage radii. Check against an independent closed-form
    ! inversion that it converges at the order 2s of the method (the last stage
    ! radius, reported before, converges only at first order), that the
    ! macrostep leaves the endpoint field cache, and that rejected steps leave
    ! the state unchanged.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_mod, only: field_can_t, evaluate
    use field_can_test, only: evaluate_test
    use orbit_symplectic_base, only: symplectic_integrator_t, sympl_rmax, &
                                     GAUSS1, SYMPLECTIC_STEP_OK
    use orbit_symplectic, only: orbit_sympl_init, orbit_timestep_sympl, &
                                orbit_timestep_sympl_rk_gauss
    implicit none
    real(dp), parameter :: mu = 0.2_dp, ro0 = 0.01_dp
    real(dp), parameter :: start(4) = [0.2_dp, 0.7_dp, 0.3_dp, &
                                       0.4_dp*(1 + 0.2_dp*cos(0.7_dp)) &
                                       - (0.2_dp**2/2 - 0.2_dp**4)/ro0]
    integer, parameter :: nlevel = 4
    real(dp), parameter :: dt0 = 0.4_dp
    integer :: stages, level, status, k, rejected
    type(field_can_t) :: f, initial, saved_field
    type(symplectic_integrator_t) :: si, saved
    real(dp) :: radius_error(nlevel), rate
    real(dp) :: exact_radius, physical_pth, physical_h

    evaluate => evaluate_test
    sympl_rmax = 1.0_dp
    initial%mu = mu
    initial%ro0 = ro0

    do stages = 1, 4
        do level = 1, nlevel
            f = initial
            ! Higher orders reach roundoff sooner; start them at longer steps.
            call orbit_sympl_init(si, f, start, &
                                  dt0*2**max(0, stages - 2)/2**(level - 1), 1, &
                                  1e-13_dp, GAUSS1 + stages - 1)
            call orbit_timestep_sympl_rk_gauss(si, f, stages, status)
            if (status /= SYMPLECTIC_STEP_OK) error stop 'Gauss interior step'
            exact_radius = inverted_radius(si%z, f%pth)
            radius_error(level) = abs(si%z(1) - exact_radius)
        end do
        rate = observed_rate(radius_error)
        print '(a,i0,a,4es10.2,a,f5.2)', 'Gauss', stages, &
            ' endpoint radius error', radius_error, '  order', rate
        ! Gauss4's eighth order reaches roundoff before it is asymptotic in
        ! double precision; require at least sixth order there.
        if (rate < min(2*stages, 6) - 0.5_dp) &
            error stop 'endpoint radius below method order'
    end do

    ! A full macrostep refreshes the field cache at the endpoint.
    do stages = 1, 4
        f = initial
        call orbit_sympl_init(si, f, start, 0.05_dp, 3, 1e-13_dp, &
                              GAUSS1 + stages - 1)
        call orbit_timestep_sympl(si, f, status)
        if (status /= SYMPLECTIC_STEP_OK) error stop 'Gauss macrostep'
        call analytic_values(si%z, physical_pth, physical_h)
        ! The cache is evaluated at the reported state itself; how close that
        ! state is to the canonical p_theta is the order check above.
        if (abs(physical_h - f%H) > 1e-12_dp) error stop 'stale field cache energy'
    end do

    ! Rejected steps near the boundary roll back the accepted state.
    rejected = 0
    do k = 1, 30
        sympl_rmax = 0.2_dp + k*1e-6_dp
        f = initial
        call orbit_sympl_init(si, f, start, -0.01_dp, 1, 1e-12_dp, GAUSS1 + 1)
        saved = si
        saved_field = f
        call orbit_timestep_sympl_rk_gauss(si, f, 2, status)
        if (status /= SYMPLECTIC_STEP_OK) then
            rejected = rejected + 1
            call unchanged()
        end if
    end do
    if (rejected == 0) error stop 'rollback cases not exercised'
    print *, 'Gauss endpoint order, cache and rollback checks passed'

contains

    real(dp) function observed_rate(errors)
        ! Least-squares slope of log2(error) over the halving levels, using only
        ! levels above the roundoff floor.
        real(dp), intent(in) :: errors(:)
        real(dp) :: x(size(errors)), y(size(errors)), xm, ym
        integer :: i, n

        n = 0
        do i = 1, size(errors)
            if (errors(i) > 1e-14_dp) then
                n = n + 1
                x(n) = real(i - 1, dp)
                y(n) = log(errors(i))/log(2.0_dp)
            end if
        end do
        if (n < 2) error stop 'too few levels above roundoff'
        xm = sum(x(1:n))/n
        ym = sum(y(1:n))/n
        observed_rate = -sum((x(1:n) - xm)*(y(1:n) - ym))/sum((x(1:n) - xm)**2)
    end function observed_rate

    real(dp) function inverted_radius(z, target_pth)
        ! Bisection on the closed-form p_theta; independent of the production
        ! field jets and of the extrapolation under test.
        real(dp), intent(in) :: z(4), target_pth
        real(dp) :: low, high, mid, pth, h
        integer :: iteration

        low = 0.05_dp
        high = 0.4_dp
        do iteration = 1, 200
            mid = 0.5_dp*(low + high)
            call analytic_values([mid, z(2:4)], pth, h)
            if (pth < target_pth) then
                low = mid
            else
                high = mid
            end if
        end do
        inverted_radius = 0.5_dp*(low + high)
    end function inverted_radius

    subroutine analytic_values(z, pth, h)
        real(dp), intent(in) :: z(4)
        real(dp), intent(out) :: pth, h
        real(dp) :: r, c, vpar

        r = z(1)
        c = cos(z(2))
        ! Closed expressions for field_can_test; no production field jets or
        ! inversion routine are used as the physical-state oracle.
        vpar = (z(4) + (r*r/2 - r**4)/ro0)/(1 + r*c)
        pth = (r*r/2 - r**3*c/3)/ro0 + vpar*(1 - 4*r*r)*r*r
        h = vpar*vpar/2 + mu*(1 - r*c)
    end subroutine analytic_values

    subroutine unchanged()
        if (any(si%z /= saved%z) .or. f%pth /= saved_field%pth .or. &
            f%H /= saved_field%H .or. f%Bmod /= saved_field%Bmod .or. &
            f%vpar /= saved_field%vpar) error stop 'failure mutated accepted state'
    end subroutine unchanged
end program test_gc_gauss_endpoint
