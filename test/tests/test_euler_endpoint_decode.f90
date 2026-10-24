module euler_endpoint_witness
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use field_can_mod, only: field_can_t
    implicit none
    integer :: backend_kind = 0
contains
    subroutine witness_field(f, r, theta, phi, mode)
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: r, theta, phi
        integer, intent(in) :: mode

        ! Synthetic canonical Hamiltonian, not machine or equilibrium data:
        ! p_theta = r + .2 sin(theta), H = p_phi**2/2 + mu*(1 + r).
        f%Ath = r + 0.2_dp*sin(theta)
        f%Aph = 0.0_dp
        f%hth = 0.0_dp
        f%hph = 1.0_dp
        f%Bmod = 1.0_dp + r
        f%Phie = 0.0_dp
        f%dAth = [1.0_dp, 0.2_dp*cos(theta), 0.0_dp]
        f%dAph = 0.0_dp
        f%dhth = 0.0_dp
        f%dhph = 0.0_dp
        f%dBmod = [1.0_dp, 0.0_dp, 0.0_dp]
        f%dPhie = 0.0_dp
        f%d2Ath = [0.0_dp, 0.0_dp, 0.0_dp, &
                   -0.2_dp*sin(theta), 0.0_dp, 0.0_dp]
        f%d2Aph = 0.0_dp
        f%d2hth = 0.0_dp
        f%d2hph = 0.0_dp
        f%d2Bmod = 0.0_dp
        f%d2Phie = 0.0_dp
        if (backend_kind == 1) then
            f%Ath = 0.2_dp*sin(theta)
            f%dAth(1) = 0.0_dp
        else if (backend_kind == 2) then
            f%Ath = ieee_value(0.0_dp, ieee_quiet_nan)
        end if
    end subroutine witness_field

    function field_values(f) result(values)
        type(field_can_t), intent(in) :: f
        real(dp) :: values(107)
        ! Compare every observable component, without hashing object padding.
        values = [f%Ath, f%Aph, f%hth, f%hph, f%Bmod, f%Phie, &
                  f%dAth, f%dAph, f%dhth, f%dhph, f%dBmod, f%dPhie, &
                  f%d2Ath, f%d2Aph, f%d2hth, f%d2hph, f%d2Bmod, f%d2Phie, &
                  f%H, f%pth, f%vpar, f%dvpar, f%dH, f%dpth, &
                  f%d2vpar, f%d2H, f%d2pth, f%mu, f%ro0]
    end function field_values
end module euler_endpoint_witness

program test_euler_endpoint_decode
    use euler_endpoint_witness
    use simple, only: tracer_t, decode_symplectic_state
    use simple_main, only: to_standard_z_coordinates, validate_rk_state, &
                           validate_symplectic_endpoint, symplectic_recovery_start, &
                           validate_decoded_symplectic_endpoint
    use field_can_mod, only: evaluate
    use orbit_symplectic, only: orbit_sympl_init, orbit_timestep_sympl
    use orbit_symplectic_base, only: symplectic_integrator_t, &
                                     EXPL_IMPL_EULER, IMPL_EXPL_EULER, MIDPOINT
    use vector_potentail_mod, only: torflux
    use params, only: isw_field_type, radial_electric_potential_slope
    use magfie_sub, only: TEST
    implicit none
    real(dp), parameter :: START(4) = [0.3_dp, 0.0_dp, 0.0_dp, 0.8_dp]
    real(dp), parameter :: H0 = 1.62_dp, TOL = 2.0e-12_dp
    integer :: mode, level
    real(dp) :: error(2), zstart(5), candidate(5)

    evaluate => witness_field
    torflux = 1.0_dp
    isw_field_type = TEST
    radial_electric_potential_slope = 0.0_dp
    zstart = [START(1:3), sqrt(H0), START(4)/sqrt(2.0_dp*H0)]
    do mode = EXPL_IMPL_EULER, IMPL_EXPL_EULER
        do level = 1, 2
            call check_euler(mode, 0.1_dp/2**(level - 1), error(level))
        end do
        if (error(2) <= 0.0_dp .or. error(1)/error(2) < 3.9_dp) &
            error stop 'Reconstructed Euler endpoint failed timestep convergence'
    end do
    call check_validation_policy()
    call check_midpoint()
    call check_failures()

    ! Exercise both sides of the existing 0.5 percent momentum limit.
    candidate = zstart
    candidate(4) = 1.0049_dp*zstart(4)
    call check_momentum(candidate, .true.)
    candidate = zstart
    candidate(4) = 1.0051_dp*zstart(4)
    call check_momentum(candidate, .false.)
contains
    subroutine initialize(anorb, method, dt)
        type(tracer_t), intent(out) :: anorb
        integer, intent(in) :: method
        real(dp), intent(in) :: dt

        anorb%integmode = method
        anorb%fper = 8.0_dp*atan(1.0_dp)
        anorb%f%mu = 1.0_dp
        anorb%f%ro0 = 1.0_dp
        anorb%f%dvpar = 0.0_dp
        anorb%f%dH = 0.0_dp
        anorb%f%dpth = 0.0_dp
        anorb%f%d2vpar = 0.0_dp
        anorb%f%d2H = 0.0_dp
        anorb%f%d2pth = 0.0_dp
        anorb%si%pabs = sqrt(H0)
        call orbit_sympl_init(anorb%si, anorb%f, START, dt, 1, &
                              1.0e-13_dp, method)
    end subroutine initialize

    subroutine check_euler(method, dt, energy_error)
        integer, intent(in) :: method
        real(dp), intent(in) :: dt
        real(dp), intent(out) :: energy_error
        type(tracer_t) :: anorb
        type(symplectic_integrator_t) :: initial_si, saved_si
        real(dp) :: saved_f(107), decoded(5), wrapped(5), checked(5)
        real(dp) :: recovery(5), stage_p
        real(dp) :: expected_r, expected_pth, expected_stage, fraction
        integer :: ierr
        logical :: boundary

        call initialize(anorb, method, dt)
        initial_si = anorb%si
        call orbit_timestep_sympl(anorb%si, anorb%f, ierr)
        if (ierr /= 0) error stop 'Synthetic Euler map did not converge'
        if (method == EXPL_IMPL_EULER) then
            expected_stage = START(1) + 0.2_dp*dt
            expected_pth = expected_stage
        else
            expected_stage = START(1) - 0.2_dp*sin(dt)
            expected_pth = START(1) + 0.2_dp*dt*cos(dt)
        end if
        expected_r = expected_pth - 0.2_dp*sin(dt)
        if (abs(anorb%si%z(1) - expected_stage) > TOL .or. &
            abs(anorb%si%z(2) - dt) > TOL .or. &
            abs(anorb%si%z(3) - START(4)*dt) > TOL .or. &
            abs(anorb%f%pth - expected_pth) > TOL) &
            error stop 'Euler map differs from independent canonical witness'
        saved_si = anorb%si
        saved_f = field_values(anorb%f)
        call decode_symplectic_state(anorb%si, anorb%f, method, decoded, ierr)
        if (ierr /= 0) error stop 'Euler endpoint reconstruction failed'
        wrapped = zstart
        call to_standard_z_coordinates(anorb, wrapped)
        if (maxval(abs(wrapped - decoded)) > TOL) &
            error stop 'Public endpoint conversion disagrees with decoder'
        if (abs(decoded(1) - expected_r) > TOL .or. &
            abs(decoded(1) + 0.2_dp*sin(decoded(2)) - expected_pth) > TOL) &
            error stop 'Decoded endpoint violates carried canonical momentum'
        if (abs(decoded(4)**2 - (0.5_dp*START(4)**2 + 1 + expected_r)) > TOL .or. &
            abs(sqrt(2.0_dp)*decoded(4)*decoded(5) - START(4)) > TOL) &
            error stop 'Decoded physical energy or pitch differs from oracle'
        energy_error = abs(decoded(4)**2 - H0)
        checked = decoded
        call validate_symplectic_endpoint(initial_si, anorb, zstart, decoded, &
                                          checked, fraction, ierr, boundary)
        if (ierr /= 0 .or. boundary) &
            error stop 'Valid reconstructed Euler endpoint was rejected'
        call assert_unchanged(anorb, saved_si, saved_f)

        ! RK recovery and reseeding restart from the mixed-stage state the map
        ! carries (stage radius, advanced angles, pitch from the cached field),
        ! not from the decoded endpoint, and always with the reference momentum.
        call symplectic_recovery_start(anorb%si, anorb%f, decoded, .true., &
                                       recovery)
        stage_p = sqrt(anorb%f%mu*anorb%f%Bmod + 0.5_dp*anorb%f%vpar**2)
        if (any(recovery(1:3) /= saved_si%z(1:3)) .or. &
            abs(recovery(1) - expected_stage) > TOL .or. &
            abs(recovery(1) - decoded(1)) <= 1.0e3_dp*TOL .or. &
            recovery(4) /= anorb%si%pabs .or. &
            abs(recovery(5) - anorb%f%vpar/(sqrt(2.0_dp)*stage_p)) > TOL) &
            error stop 'RK recovery does not start from the mixed-stage state'
        call symplectic_recovery_start(anorb%si, anorb%f, decoded, .false., &
                                       recovery)
        if (any(recovery([1, 2, 3, 5]) /= decoded([1, 2, 3, 5])) .or. &
            recovery(4) /= anorb%si%pabs) &
            error stop 'RK recovery changed a state not decoded from the map'
        call assert_unchanged(anorb, saved_si, saved_f)

        ! An actual wrong-energy canonical state must still be rejected.
        anorb%si%z(4) = 1.1_dp*START(4)
        call decode_symplectic_state(anorb%si, anorb%f, method, decoded, ierr)
        if (ierr /= 0) error stop 'Wrong-energy witness could not be decoded'
        checked = decoded
        call validate_symplectic_endpoint(initial_si, anorb, zstart, decoded, &
                                          checked, fraction, ierr, boundary)
        if (ierr == 0 .or. any(checked /= zstart)) &
            error stop 'Endpoint validator accepted a true energy violation'
    end subroutine check_euler

    subroutine check_validation_policy()
        type(tracer_t) :: anorb
        type(symplectic_integrator_t) :: initial_si
        type(field_can_t) :: initial_f
        real(dp) :: decoded(5), checked(5), fraction, dt
        integer :: method, level, ierr
        logical :: boundary

        do method = EXPL_IMPL_EULER, IMPL_EXPL_EULER
            do level = 1, 2
                dt = 0.1_dp/2**(level - 1)
                call initialize(anorb, method, dt)
                initial_si = anorb%si
                initial_f = anorb%f
                call orbit_timestep_sympl(anorb%si, anorb%f, ierr)
                if (ierr /= 0) error stop 'Validation witness map failed'
                call decode_symplectic_state(anorb%si, anorb%f, method, &
                                             decoded, ierr)
                if (ierr /= 0) error stop 'Validation witness decode failed'
                checked = decoded
                call validate_decoded_symplectic_endpoint(initial_si, &
                                  initial_f, anorb, zstart, .false., decoded, checked, &
                                                          fraction, ierr, boundary)
                ! At dt=.1 the mixed-stage energy differs by over 0.5 percent
                ! in momentum, while the decoded energy is within that limit.
                if (level == 1) then
                    if (ierr == 0 .or. any(checked /= zstart)) &
                        error stop 'Decoded energy changed legacy acceptance'
                else
                    if (ierr /= 0 .or. boundary .or. &
                        any(checked /= decoded)) &
                        error stop 'Accepted step did not report decoded endpoint'
                end if
            end do
        end do
    end subroutine check_validation_policy

    subroutine check_midpoint()
        type(tracer_t) :: anorb
        type(symplectic_integrator_t) :: saved_si
        real(dp) :: decoded(5), saved_f(107)
        integer :: ierr

        call initialize(anorb, MIDPOINT, 0.1_dp)
        call orbit_timestep_sympl(anorb%si, anorb%f, ierr)
        if (ierr /= 0) error stop 'Synthetic midpoint map did not converge'
        saved_si = anorb%si
        saved_f = field_values(anorb%f)
        call decode_symplectic_state(anorb%si, anorb%f, MIDPOINT, decoded, ierr)
        if (ierr /= 0 .or. any(decoded(1:3) /= saved_si%z(1:3))) &
            error stop 'Midpoint endpoint coordinates changed during decoding'
        call assert_unchanged(anorb, saved_si, saved_f)
    end subroutine check_midpoint

    subroutine check_failures()
        type(tracer_t) :: anorb
        real(dp) :: decoded(5)
        integer :: kind, ierr

        call initialize(anorb, EXPL_IMPL_EULER, 0.1_dp)
        do kind = 1, 2
            backend_kind = kind
            call decode_symplectic_state(anorb%si, anorb%f, &
                                         EXPL_IMPL_EULER, decoded, ierr)
            if (ierr == 0) &
                error stop 'Singular or nonfinite endpoint constraint was accepted'
        end do
        backend_kind = 0
    end subroutine check_failures

    subroutine check_momentum(z, accepted)
        real(dp), intent(inout) :: z(5)
        logical, intent(in) :: accepted
        integer :: ierr

        ierr = 0
        call validate_rk_state(zstart, z, ierr, zstart(4))
        if ((ierr == 0) .neqv. accepted) &
            error stop 'The 0.5 percent momentum sanity limit changed'
    end subroutine check_momentum

    subroutine assert_unchanged(anorb, si, values)
        type(tracer_t), intent(in) :: anorb
        type(symplectic_integrator_t), intent(in) :: si
        real(dp), intent(in) :: values(107)

        if (any(field_values(anorb%f) /= values)) &
            error stop 'Endpoint conversion mutated the carried field cache'
        if (any(anorb%si%z /= si%z) .or. anorb%si%ntau /= si%ntau .or. &
            any([anorb%si%atol, anorb%si%rtol, anorb%si%warning_factor, &
                 anorb%si%pthold, anorb%si%dt, anorb%si%pabs, &
                 anorb%si%last_step_fraction, anorb%si%last_event_radial_residual, &
                 anorb%si%last_event_fraction_width] /= &
                [si%atol, si%rtol, si%warning_factor, si%pthold, si%dt, si%pabs, &
                 si%last_step_fraction, si%last_event_radial_residual, &
                 si%last_event_fraction_width])) &
            error stop 'Endpoint conversion mutated the canonical map state'
    end subroutine assert_unchanged
end program test_euler_endpoint_decode
