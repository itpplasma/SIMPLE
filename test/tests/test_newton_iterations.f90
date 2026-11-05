program test_newton_iterations
    ! Newton cost and accuracy of the implicit symplectic steps on the
    ! stellarator equilibrium of test_sympl_stell. With the extrapolated starting
    ! guess and the contraction-rate stop, each Newton solve needs two
    ! iterations, one field evaluation per stage and iteration. The oracle for
    ! accuracy is the same scheme run from z_n with the correction stop, which
    ! iterates until the last correction is below tolerance.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use field_can_base, only: n_field_evaluations_d2
    use field_can_mod, only: eval_field => evaluate, field_can_t
    use orbit_symplectic, only: orbit_sympl_init, orbit_timestep_sympl
    use orbit_symplectic_base, only: symplectic_integrator_t, sympl_predictor, &
                    sympl_stop_rule, PREDICT_NONE, PREDICT_QUADRATIC, STOP_CORRECTION, &
                                STOP_RATE, EXPL_IMPL_EULER, IMPL_EXPL_EULER, MIDPOINT, &
                                     GAUSS1, GAUSS2, GAUSS3, GAUSS4, LOBATTO3
    use simple, only: tracer_t, init_params
    use simple_main, only: init_field
    use params, only: coord_input, field_input
    use new_vmec_stuff_mod, only: rmajor
    use util, only: twopi
    use parmot_mod, only: ro0

    implicit none

    integer, parameter :: nsteps = 400
    real(dp), parameter :: rtol = 1d-12
    type(tracer_t) :: norb
    type(field_can_t) :: f0
    real(dp) :: z0(4), dt
    integer :: nfail

    nfail = 0
    call setup(z0, f0, dt)

    ! Field evaluations per step: stages times two Newton iterations, with a
    ! margin for the few steps that need a third.
    call check_scheme('euler1', EXPL_IMPL_EULER, 1, nfail)
    call check_scheme('midpoint', MIDPOINT, 2, nfail)
    call check_scheme('gauss4', GAUSS2, 2, nfail)
    call check_scheme('lobatto4', LOBATTO3, 3, nfail)
    call check_history_reset(nfail)
    call check_iteration_counts(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, ' check(s) failed'
        error stop 1
    end if
    write (*, '(a)') 'test_newton_iterations: all checks passed'

contains

    subroutine setup(z0, f, dt)
        real(dp), intent(out) :: z0(4), dt
        type(field_can_t), intent(out) :: f
        real(dp) :: vpar0

        field_input = 'wout.nc'
        coord_input = 'wout.nc'
        call init_field(norb, 'wout.nc', 5, 5, 3, 2)
        call init_params(norb, 2, 4, 3.5d6, 64, 1, 1d-8)
        z0(1:3) = [0.4d0, 0.7d0, 0.1d0]
        vpar0 = 0.8d0
        call eval_field(f, z0(1), z0(2), z0(3), 0)
        f%mu = 0.5d0**2*(1d0 - vpar0**2)/f%Bmod*2d0
        f%ro0 = ro0/sqrt(2d0)
        f%vpar = vpar0*sqrt(2d0)
        z0(4) = vpar0*f%hph + f%Aph/f%ro0
        dt = twopi*rmajor*1d2/48d0
    end subroutine setup

    subroutine run(mode, predictor, stop_rule, z, evaluations_per_step, evaluations)
        integer, intent(in) :: mode, predictor, stop_rule
        real(dp), intent(out) :: z(4), evaluations_per_step
        integer(8), intent(out), optional :: evaluations
        type(symplectic_integrator_t) :: si
        type(field_can_t) :: f
        integer(8) :: d2_start
        integer :: k, ierr

        sympl_predictor = predictor
        sympl_stop_rule = stop_rule
        f = f0
        call orbit_sympl_init(si, f, z0, dt, 1, rtol, mode)
        d2_start = n_field_evaluations_d2
        do k = 1, nsteps
            call orbit_timestep_sympl(si, f, ierr)
            if (ierr /= 0) error stop 'symplectic step failed'
        end do
        evaluations_per_step = real(n_field_evaluations_d2 - d2_start, dp)/nsteps
        if (present(evaluations)) evaluations = n_field_evaluations_d2 - d2_start
        z = si%z
        sympl_predictor = PREDICT_QUADRATIC
        sympl_stop_rule = STOP_RATE
    end subroutine run

    subroutine check_scheme(name, mode, stages, nfail)
        character(*), intent(in) :: name
        integer, intent(in) :: mode, stages
        integer, intent(inout) :: nfail
        real(dp) :: z(4), zref(4), cost, cost_ref, deviation

        call run(mode, PREDICT_QUADRATIC, STOP_RATE, z, cost)
        call run(mode, PREDICT_NONE, STOP_CORRECTION, zref, cost_ref)
        ! Tolerance-sized per-step differences add up over the run.
        deviation = maxval(abs(z - zref)/[1d0, twopi, twopi, abs(zref(4))])
        write (*, '(a10,a,f7.3,a,f7.3,a,es9.2)') name, ': evaluations/step ', &
            cost, ' (reference ', cost_ref, '), deviation ', deviation
        if (cost > 2.2d0*stages) then
            write (error_unit, '(a,a)') '  too many Newton iterations: ', name
            nfail = nfail + 1
        end if
        if (deviation > nsteps*rtol) then
           write (error_unit, '(a,a)') '  deviates from the converged reference: ', name
            nfail = nfail + 1
        end if
    end subroutine check_scheme

    subroutine check_history_reset(nfail)
        ! A state changed from outside (interface map, chart switch) must not
        ! reuse the old history: the next step then equals a fresh start.
        integer, intent(inout) :: nfail
        type(symplectic_integrator_t) :: si, fresh
        type(field_can_t) :: f, f_fresh
        real(dp) :: z(4)
        integer :: k, ierr

        f = f0
        call orbit_sympl_init(si, f, z0, dt, 1, rtol, MIDPOINT)
        do k = 1, 20
            call orbit_timestep_sympl(si, f, ierr)
        end do
        z = si%z
        z(1) = z(1) + 1d-3
        f_fresh = f0
        call orbit_sympl_init(fresh, f_fresh, z, dt, 1, rtol, MIDPOINT)
        si%z = z
        si%pthold = fresh%pthold
        f = f_fresh
        call orbit_timestep_sympl(si, f, ierr)
        call orbit_timestep_sympl(fresh, f_fresh, ierr)
        if (any(si%z /= fresh%z)) then
         write (error_unit, '(a)') '  stale history used after an external state change'
            nfail = nfail + 1
        end if
    end subroutine check_history_reset

    subroutine check_iteration_counts(nfail)
        ! Refactor guard for the Newton drivers: field evaluations, one per stage
        ! and Newton iteration, as recorded with the per-scheme drivers that
        ! preceded the shared one, for the rate stop and the correction stop.
        integer, intent(inout) :: nfail
        integer, parameter :: nschemes = 8
        character(8), parameter :: names(nschemes) = [character(8) :: 'euler1', &
               'euler2', 'midpoint', 'gauss1', 'gauss2', 'gauss3', 'gauss4', 'lobatto3']
        integer, parameter :: modes(nschemes) = [EXPL_IMPL_EULER, IMPL_EXPL_EULER, &
                                     MIDPOINT, GAUSS1, GAUSS2, GAUSS3, GAUSS4, LOBATTO3]
        integer, parameter :: rate_counts(nschemes) = &
                              [802, 964, 1678, 824, 1602, 2403, 3204, 2403]
        integer, parameter :: correction_counts(nschemes) = &
                              [1332, 1698, 3200, 1600, 3200, 4800, 6400, 4800]
        real(dp) :: z(4), cost
        integer(8) :: counts(2)
        integer :: k

        do k = 1, nschemes
            call run(modes(k), PREDICT_QUADRATIC, STOP_RATE, z, cost, counts(1))
            call run(modes(k), PREDICT_NONE, STOP_CORRECTION, z, cost, counts(2))
            write (*, '(a10,a,2i8)') names(k), ': field evaluations ', counts
            if (counts(1) /= rate_counts(k) .or. &
                counts(2) /= correction_counts(k)) then
                write (error_unit, '(a,a,a,2i8)') '  Newton iterations changed: ', &
                    trim(names(k)), ', expected', rate_counts(k), correction_counts(k)
                nfail = nfail + 1
            end if
        end do
    end subroutine check_iteration_counts

end program test_newton_iterations
