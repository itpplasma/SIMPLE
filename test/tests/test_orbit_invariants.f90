program test_orbit_invariants
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_mod, only: evaluate_field => evaluate
    use magfie_sub, only: TEST
    use orbit_invariants, only: guiding_center_invariants_t, &
        invariant_start_result_t, invariants_from_state, states_from_invariants, &
        potato_to_simple_invariants, invariant_flux_convention, INVARIANT_SUCCESS
    use orbit_symplectic_base, only: symplectic_integrator_t
    use params, only: params_init
    use parmot_mod, only: rho0 => ro0
    use simple, only: tracer_t, encode_symplectic_state, init_sympl
    use simple_main, only: init_field
    use util, only: sqrt2
    use velo_mod, only: isw_field_type

    implicit none

    type(tracer_t) :: tracer, work
    type(symplectic_integrator_t) :: symplectic
    type(guiding_center_invariants_t) :: invariants, reconstructed
    type(guiding_center_invariants_t) :: potato, converted
    type(invariant_start_result_t) :: result
    real(dp) :: state(5), zsympl(4), psi_native
    real(dp) :: axis_flux, outward_sign
    integer :: i, status

    isw_field_type = TEST
    call init_field(tracer, '', 5, 5, 3, 3)
    call params_init
    rho0 = 0.02_dp

    state = [0.2_dp, 0.7_dp, 0.2_dp, 0.8_dp, 0.6_dp]
    call invariants_from_state(tracer, state, invariants, status)
    call require(status == INVARIANT_SUCCESS, 'forward conversion status')

    work = tracer
    call evaluate_field(work%f, state(1), state(2), state(3), 0)
    call invariant_flux_convention(tracer, axis_flux, outward_sign)
    psi_native = work%f%Aph + rho0*state(4)*state(5)*work%f%hph
    call require_close(invariants%h0, state(4)**2, 1.0e-13_dp, 'H0')
    call require_close(invariants%j_perp, &
        state(4)**2*(1.0_dp - state(5)**2)/work%f%Bmod, 1.0e-13_dp, 'J_perp')
    call require_close(invariants%p_phi, &
        outward_sign*(psi_native - axis_flux), 1.0e-13_dp, 'p_phi')

    work = tracer
    call encode_symplectic_state(work%f, state, zsympl)
    call require_close(zsympl(4), sqrt2*psi_native/rho0, 1.0e-13_dp, &
        'private symplectic p_phi normalization')

    work = tracer
    call init_sympl(symplectic, work%f, state, 0.1_dp, 0.1_dp, 1.0e-10_dp, 3)
    call require_close(symplectic%z(4), zsympl(4), 1.0e-13_dp, &
        'shared symplectic encoding')
    call require_close(symplectic%dt*sqrt2, 0.1_dp, 1.0e-13_dp, &
        'public tau versus private symplectic time')

    call states_from_invariants(tracer, invariants, result)
    call require(result%status == INVARIANT_SUCCESS, 'inverse conversion status')
    call require(size(result%starts) >= 2, 'all cut intersections returned')
    do i = 1, size(result%starts)
        call invariants_from_state(tracer, result%starts(i)%state, reconstructed, &
            status)
        call require(status == INVARIANT_SUCCESS, 'reconstructed status')
        call require_close(reconstructed%h0, invariants%h0, 1.0e-11_dp, &
            'reconstructed H0')
        call require_close(reconstructed%j_perp, invariants%j_perp, 1.0e-10_dp, &
            'reconstructed J_perp')
        call require_close(reconstructed%p_phi, invariants%p_phi, 1.0e-9_dp, &
            'reconstructed p_phi')
    end do

    potato = guiding_center_invariants_t(1.0_dp, 2.0_dp, 12.0_dp)
    call potato_to_simple_invariants(potato, 10.0_dp, 20.0_dp, 0.5_dp, converted)
    call require_close(converted%h0, 0.25_dp, 1.0e-14_dp, 'POTATO H0 scale')
    call require_close(converted%j_perp, 0.5_dp, 1.0e-14_dp, &
        'POTATO J_perp scale')
    call require_close(converted%p_phi, 2.0_dp, 1.0e-14_dp, &
        'POTATO outward flux convention')
    potato%p_phi = 8.0_dp
    call potato_to_simple_invariants(potato, 10.0_dp, 0.0_dp, 1.0_dp, converted)
    call require_close(converted%p_phi, 2.0_dp, 1.0e-14_dp, &
        'POTATO reversed-flux convention')

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) then
            write (*, '(A)') 'FAILED: '//message
            error stop 1
        end if
    end subroutine require

    subroutine require_close(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual, expected, tolerance
        character(*), intent(in) :: message

        call require(abs(actual - expected) <= tolerance*max(1.0_dp, abs(expected)), &
            message)
    end subroutine require_close

end program test_orbit_invariants
