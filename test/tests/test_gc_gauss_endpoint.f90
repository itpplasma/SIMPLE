program test_gc_gauss_endpoint
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use field_can_mod, only: field_can_t, evaluate, get_val
    use field_can_test, only: evaluate_test
    use orbit_symplectic_base, only: symplectic_integrator_t, sympl_rmax, &
                              GAUSS1, SYMPLECTIC_STEP_OK, SYMPLECTIC_STEP_OUTSIDE_DOMAIN
    use orbit_symplectic, only: orbit_sympl_init, &
                               orbit_timestep_sympl_rk_gauss, reconstruct_gauss_endpoint
    implicit none
    real(dp), parameter :: mu = 0.2_dp, ro0 = 0.01_dp
    real(dp), parameter :: start(4) = [0.2_dp, 0.7_dp, 0.3_dp, &
                        0.4_dp*(1 + 0.2_dp*cos(0.7_dp)) - (0.2_dp**2/2 - 0.2_dp**4)/ro0]
    integer :: stages, substeps, level, status, k, rejected
    type(field_can_t) :: f, initial, target, saved_field
    type(symplectic_integrator_t) :: si, saved
    real(dp) :: h0, physical_pth, physical_h, errors(3), target_pth

    evaluate => evaluate_test
    sympl_rmax = 1.0_dp
    initial%mu = mu
    initial%ro0 = ro0
    call orbit_sympl_init(si, initial, start, 0.4_dp, 1, 1e-13_dp, GAUSS1)
    call analytic_values(start, physical_pth, h0)

    ! The built-in circular-tokamak field has nonzero h_theta. Check the
    ! returned physical radius and cache, including multiple internal steps.
    do stages = 1, 4
        do substeps = 1, 3, 2
            do level = 1, 3
                f = initial
                call orbit_sympl_init(si, f, start, 0.4_dp/2**(level - 1), &
                                      substeps, 1e-13_dp, GAUSS1 + stages - 1)
                call orbit_timestep_sympl_rk_gauss(si, f, stages, status)
                if (status /= SYMPLECTIC_STEP_OK) error stop 'Gauss interior step'
                call analytic_values(si%z, physical_pth, physical_h)
                if (abs(physical_pth - f%pth) > 3e-12_dp) &
                    error stop 'Gauss returned collocation radius instead of endpoint'
                if (abs(physical_h - f%H) > 3e-13_dp) &
                    error stop 'Gauss returned stale collocation field cache'
                if (stages == 2 .and. substeps == 1) &
                    errors(level) = abs((physical_h - h0)/h0)
            end do
        end do
    end do
    ! Fourth-order Gauss2 has fifth-order one-step error on this smooth field.
    if (any(errors(1:2)/errors(2:3) < 16.0_dp)) &
        error stop 'Gauss2 physical Hamiltonian local convergence'

    ! Invert a known finite-h_theta target from a deliberately displaced guess.
    f = initial
    target = initial
    call evaluate_test(target, 0.23_dp, start(2), start(3), 0)
    call get_val(target, start(4))
    target_pth = target%pth
    call orbit_sympl_init(si, f, start, 0.01_dp, 1, 1e-12_dp, GAUSS1)
    f%pth = target_pth
    call reconstruct_gauss_endpoint(si, f, status)
    if (status /= SYMPLECTIC_STEP_OK .or. abs(si%z(1) - 0.23_dp) > 2e-12_dp) &
        error stop 'finite-h_theta canonical inversion'
    if (f%pth /= target_pth) error stop 'canonical update changed during inversion'

    ! Reject genuine boundary crossings without changing an accepted state.
    sympl_rmax = 0.21_dp
    si%z = start
    f = initial
    f%pth = target_pth
    saved = si
    saved_field = f
    call reconstruct_gauss_endpoint(si, f, status)
    if (status /= SYMPLECTIC_STEP_OUTSIDE_DOMAIN) error stop 'outer crossing status'
    call unchanged()
    sympl_rmax = 1.0_dp
    si%z(1) = -0.01_dp
    saved = si
    saved_field = f
    call reconstruct_gauss_endpoint(si, f, status)
    if (status == SYMPLECTIC_STEP_OK) error stop 'negative radius accepted'
    call unchanged()
    si%z = start
    f = initial
    f%pth = ieee_value(0.0_dp, ieee_quiet_nan)
    saved = si
    call reconstruct_gauss_endpoint(si, f, status)
    if (status == SYMPLECTIC_STEP_OK .or. any(si%z /= saved%z)) &
        error stop 'nonfinite canonical target accepted'

    ! Exercise the raw step's rollback, not only the reconstruction helper.
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
    print *, 'Gauss endpoint, cache, local order and rollback checks passed'
contains
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
