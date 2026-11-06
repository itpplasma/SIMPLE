module orbit_symplectic

use, intrinsic :: iso_fortran_env, only: dp => real64, int64
use util, only: pi, twopi
!$  use omp_lib, only: omp_in_parallel
use field_can_mod, only: field_can_t, get_val, get_derivatives, get_derivatives2, &
  eval_field => evaluate
use orbit_symplectic_base, only: symplectic_integrator_t, multistage_integrator_t, &
  RK45, EXPL_IMPL_EULER, IMPL_EXPL_EULER, MIDPOINT, GAUSS1, GAUSS2, GAUSS3, GAUSS4, &
  LOBATTO3, RADAU15, GBS16, CASHKARP45, TDRK24, TDRK24A, S_MAX, orbit_timestep_sympl_i, extrap_field, sympl_rmax, &
  coeff_rk_gauss, coeff_rk_lobatto, f_rk_lobatto, &
  SYMPLECTIC_STEP_OK, SYMPLECTIC_STEP_OUTSIDE_DOMAIN, &
  SYMPLECTIC_STEP_MAXITER, SYMPLECTIC_STEP_LINEAR_SOLVE, &
  SYMPLECTIC_STEP_BOUNDARY, SYMPLECTIC_STEP_EVENT_NOT_CONVERGED, &
  SYMPLECTIC_STEP_BOUNDARY_LIMITED, boundary_event_fraction_tolerance, &
                      boundary_event_radial_tolerance, symplectic_newton_warning_mode, &
                         history_valid, predict_state, record_history, rate_converged, &
                                     max_retry_depth, retry_split, retry_next
    use orbit_symplectic_quasi, only: &
                  orbit_timestep_radau15, orbit_timestep_gbs16, orbit_timestep_tdrk24, &
                                    orbit_timestep_tdrk24a, orbit_timestep_cashkarp45, &
                                      reset_explicit_step_carry, &
                                      si_quasi => si, f_quasi => f
use orbit_symplectic_euler1, only: sympl_euler1_residual, sympl_euler1_jacobian, &
                                   sympl_euler1_solve, sympl_euler1_extrapolate_field, &
  sympl_euler1_advance_angles
    use orbit_systems_generated, only: euler2_residual, euler2_jacobian, &
                                 midpoint_residual, midpoint_jacobian, lobatto3_jacobian
    use orbit_systems_generated, only: gc_geometry_gradient_transport
    use orbit_newton, only: newton_state_t, newton_begin, newton_advance, &
                            SYSTEM_EULER1, SYSTEM_EULER2, SYSTEM_MIDPOINT, &
                            SYSTEM_RK_GAUSS, SYSTEM_RK_LOBATTO, NEWTON_NMAX, &
                            WARNING_NONE, WARNING_ACCEPT, WARNING_REJECT_100, &
                            WARNING_REJECT_1E4, WARNING_REJECT_LARGE, &
                            WARNING_NONFINITE, WARNING_AXIS_ACCEPT, &
                            finite_newton_iterate, warning_maxiter_class, &
                            axis_crossing_converged, matrix3_near_singular, &
                            limit_radial_newton_step, radial_boundary_reached
use orbit_rk_core, only: rk_gauss_residual, rk_gauss_jacobian, rk_solve, MODEL_GC
use vector_potentail_mod, only: torflux
use diag_counters, only: count_event, EVT_NEWTON1_MAXIT, EVT_NEWTON2_MAXIT, &
  EVT_RK_GAUSS_MAXIT, EVT_RK_LOBATTO_MAXIT, EVT_FIXPOINT_MAXIT, EVT_R_NEGATIVE, &
  EVT_MIDPOINT_MAXIT, EVT_WARNING_MAXIT_ACCEPT, &
  EVT_WARNING_MAXIT_REJECT_100, EVT_WARNING_MAXIT_REJECT_1E4, &
  EVT_WARNING_MAXIT_REJECT_LARGE, EVT_WARNING_MAXIT_NONFINITE, &
  EVT_RETRY_EXHAUSTED, EVT_WARNING_AXIS_CROSSING_ACCEPT

implicit none

! Tests may substitute a stepper; production code never reassigns it.
    procedure(orbit_timestep_sympl_i), pointer :: orbit_timestep_sympl => &
        orbit_timestep_sympl_dispatch

contains

    logical function accept_warning_maxiter(x, xlast, tolref, rtol, warning_factor)
        ! Continue only when the final Newton correction is finite and inside this
        ! integrator's tolerance-scaled guard. Other iterates enter retry recovery.
        real(dp), intent(in) :: x(:), xlast(:), tolref(:), rtol, warning_factor
        integer :: warning

        warning = warning_maxiter_class(symplectic_newton_warning_mode, x, xlast, &
                                        tolref, rtol, warning_factor)
        call count_warning(warning)
        accept_warning_maxiter = warning == WARNING_ACCEPT
    end function accept_warning_maxiter

    subroutine count_warning(warning)
        integer, intent(in) :: warning

        select case (warning)
        case (WARNING_ACCEPT)
            call count_event(EVT_WARNING_MAXIT_ACCEPT)
        case (WARNING_REJECT_100)
            call count_event(EVT_WARNING_MAXIT_REJECT_100)
        case (WARNING_REJECT_1E4)
            call count_event(EVT_WARNING_MAXIT_REJECT_1E4)
        case (WARNING_REJECT_LARGE)
            call count_event(EVT_WARNING_MAXIT_REJECT_LARGE)
        case (WARNING_NONFINITE)
            call count_event(EVT_WARNING_MAXIT_NONFINITE)
        case (WARNING_AXIS_ACCEPT)
            call count_event(EVT_WARNING_AXIS_CROSSING_ACCEPT)
        end select
    end subroutine count_warning

    logical function accept_warning_axis_crossing(si, x, xlast, tolref, rtol)
        type(symplectic_integrator_t), intent(in) :: si
        real(dp), intent(in) :: x(5), xlast(5), tolref(5), rtol

        accept_warning_axis_crossing = axis_crossing_converged( &
                                       symplectic_newton_warning_mode, si%z(1), x, &
                                       xlast, tolref, rtol, si%warning_factor)
        if (accept_warning_axis_crossing) &
            call count_event(EVT_WARNING_AXIS_CROSSING_ACCEPT)
    end function accept_warning_axis_crossing

    subroutine advance_symplectic_with_retry(si, f, stepper, status, &
                                             accepted_fraction)
        ! For a fixed subdivision, the accepted maps remain symplectic and cover the
        ! original interval exactly. State-dependent subdivision is a warning-mode
        ! recovery path, not a globally fixed symplectic map.
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        procedure(orbit_timestep_sympl_i) :: stepper
        integer, intent(out) :: status
        real(dp), intent(out), optional :: accepted_fraction
        type(symplectic_integrator_t) :: initial_integrator, substep_integrator
        type(field_can_t) :: initial_field, substep_field
        real(dp) :: interval, fraction
        integer :: level, slot, k
        logical :: done, split

        initial_integrator = si
        initial_field = f
        interval = si%dt
        level = 0
        slot = 0
        fraction = 0.0_dp
        do
            if (level > 0) then
                substep_integrator = si
                substep_field = f
            end if
            si%dt = interval*0.5_dp**level
            call stepper(si, f, status)
            if (status == SYMPLECTIC_STEP_OK) then
                call retry_next(level, slot, done)
                if (.not. done) cycle
                fraction = 1.0_dp
                exit
            end if
            if (status == SYMPLECTIC_STEP_BOUNDARY) then
                ! Fold the sub-interval fraction outwards one halving at a time, in the
                ! order of the former recursive retry, so the arithmetic is unchanged.
                do k = 1, level
                    if (mod(slot, 2) == 0) then
                        si%last_step_fraction = 0.5_dp*si%last_step_fraction
                    else
                        si%last_step_fraction = 0.5_dp + 0.5_dp*si%last_step_fraction
                    end if
                    si%last_event_fraction_width = 0.5_dp*si%last_event_fraction_width
                    slot = slot/2
                end do
                fraction = si%last_step_fraction
                exit
            end if

            if (level > 0) then
                si = substep_integrator
                f = substep_field
            else
                si = initial_integrator
                f = initial_field
            end if
            call retry_split(status, symplectic_newton_warning_mode, level, slot, split)
            if (.not. split) then
                if (level >= max_retry_depth .and. &
                    (status == SYMPLECTIC_STEP_MAXITER .or. &
                     status == SYMPLECTIC_STEP_LINEAR_SOLVE)) &
                    call count_event(EVT_RETRY_EXHAUSTED)
                si = initial_integrator
                f = initial_field
                exit
            end if
        end do
        si%dt = interval
        if (present(accepted_fraction)) accepted_fraction = fraction
    end subroutine advance_symplectic_with_retry

subroutine apply_axis_chart_switch(radius, theta, crossed, status)
  real(dp), intent(inout) :: radius, theta
  logical, intent(out) :: crossed
  integer, intent(out) :: status

  crossed = radius < 0d0
  status = SYMPLECTIC_STEP_OK
  if (.not. crossed) return

  radius = -radius
  theta = theta + pi
  call count_event(EVT_R_NEGATIVE)
  if (radius > sympl_rmax) status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
end subroutine apply_axis_chart_switch

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_sympl_init(si, f, z, dt, ntau, rtol_init, mode_init)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  real(dp), intent(in) :: z(:)
  real(dp), intent(in) :: dt
  integer, intent(in) :: ntau
  real(dp), intent(in) :: rtol_init
  integer, intent(in) :: mode_init

  si%atol = 1d-15
  si%rtol = rtol_init
        call reset_explicit_step_carry

        si%mode = mode_init
  si%ntau = ntau
  si%dt = dt
  si%last_step_fraction = 1d0
  si%last_event_radial_residual = 0d0
  si%last_event_fraction_width = 0d0
        si%nhist = 0

  si%z = z

  call eval_field(f, z(1), z(2), z(3), 0)
  call get_val(f, si%z(4)) ! for pth
  si%pthold = f%pth

        if (.not. any(mode_init == [RK45, EXPL_IMPL_EULER, IMPL_EXPL_EULER, MIDPOINT, &
                                    GAUSS1, GAUSS2, GAUSS3, GAUSS4, LOBATTO3, RADAU15, &
                                    GBS16, TDRK24, TDRK24A, CASHKARP45])) then
            print *, 'invalid mode for orbit_timestep_sympl: ', mode_init
            error stop
        end if
end subroutine orbit_sympl_init

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
! Bridges the explicit high-order steppers, which live on the quasi path and
! read the module-level state there, into the symplectic stepper signature the
! driver calls. The quasi state is threadprivate, so this is thread safe.
!
! The copies are the price of reusing f_ode unchanged, which is what keeps the
! comparison honest: these methods integrate exactly the same right-hand side
! as the symplectic schemes, through the same field evaluations.
    recursive subroutine orbit_timestep_sympl_radau15(si, f, ierr)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        si_quasi = si
        f_quasi = f
        call orbit_timestep_radau15(ierr)
        si = si_quasi
        f = f_quasi
    end subroutine orbit_timestep_sympl_radau15

    recursive subroutine orbit_timestep_sympl_gbs16(si, f, ierr)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        si_quasi = si
        f_quasi = f
        call orbit_timestep_gbs16(ierr)
        si = si_quasi
        f = f_quasi
    end subroutine orbit_timestep_sympl_gbs16

    recursive subroutine orbit_timestep_sympl_tdrk24(si, f, ierr)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        si_quasi = si
        f_quasi = f
        call orbit_timestep_tdrk24(ierr)
        si = si_quasi
        f = f_quasi
    end subroutine orbit_timestep_sympl_tdrk24

    recursive subroutine orbit_timestep_sympl_tdrk24a(si, f, ierr)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        si_quasi = si
        f_quasi = f
        call orbit_timestep_tdrk24a(ierr)
        si = si_quasi
        f = f_quasi
    end subroutine orbit_timestep_sympl_tdrk24a

    recursive subroutine orbit_timestep_sympl_cashkarp45(si, f, ierr)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        si_quasi = si
        f_quasi = f
        call orbit_timestep_cashkarp45(ierr)
        si = si_quasi
        f = f_quasi
    end subroutine orbit_timestep_sympl_cashkarp45


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_timestep_sympl_gauss1(si, f, ierr)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  call orbit_timestep_sympl_rk_gauss(si, f, 1, ierr)
end subroutine orbit_timestep_sympl_gauss1

recursive subroutine orbit_timestep_sympl_gauss2(si, f, ierr)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  call orbit_timestep_sympl_rk_gauss(si, f, 2, ierr)
end subroutine orbit_timestep_sympl_gauss2

recursive subroutine orbit_timestep_sympl_gauss3(si, f, ierr)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  call orbit_timestep_sympl_rk_gauss(si, f, 3, ierr)
end subroutine orbit_timestep_sympl_gauss3

recursive subroutine orbit_timestep_sympl_gauss4(si, f, ierr)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  call orbit_timestep_sympl_rk_gauss(si, f, 4, ierr)
end subroutine orbit_timestep_sympl_gauss4

recursive subroutine orbit_timestep_sympl_lobatto3(si, f, ierr)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  call orbit_timestep_sympl_rk_lobatto(si, f, 3, ierr)
end subroutine orbit_timestep_sympl_lobatto3



  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine f_sympl_euler1(si, f, n, x, fvec, iflag)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(in) :: n
  real(dp), intent(in) :: x(n)
  real(dp), intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  call eval_field(f, x(1), si%z(2), si%z(3), 2)
  call get_derivatives2(f, x(2))
  call sympl_euler1_residual(si, f, x, fvec)

end subroutine f_sympl_euler1


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine jac_sympl_euler1(si, f, x, jac)
  !
  type(symplectic_integrator_t), intent(in) :: si
  type(field_can_t), intent(inout) :: f

  real(dp), intent(in)  :: x(2)
  real(dp), intent(out) :: jac(2, 2)

  call sympl_euler1_jacobian(si, f, x, jac)

end subroutine jac_sympl_euler1

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine f_sympl_euler2(si, f, n, x, fvec, iflag)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(in) :: n
  real(dp), intent(in) :: x(n)
  real(dp), intent(out) :: fvec(n)
  integer, intent(in) :: iflag

  call eval_field(f, x(1), x(2), x(3), 2)
  call get_derivatives2(f, si%z(4))
        call euler2_residual(si%dt, si%pthold, si%z, x, f, fvec)

end subroutine f_sympl_euler2


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine jac_sympl_euler2(si, f, x, jac)
  !
  type(symplectic_integrator_t), intent(in) :: si
  type(field_can_t), intent(inout) :: f
  real(dp), intent(in)  :: x(3)
  real(dp), intent(out) :: jac(3, 3)

        call euler2_jacobian(si%dt, si%z, x, f, jac)

end subroutine jac_sympl_euler2


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
    recursive subroutine midpoint_system(si, f, fmid, x, fvec, jac)
        ! Midpoint residual and Jacobian at x = (rend, thend, phend, pphend, rmid).
        ! On return f holds the field at the end point, fmid at the midpoint.
        type(symplectic_integrator_t), intent(in) :: si
        type(field_can_t), intent(inout) :: f, fmid
        real(dp), intent(in) :: x(5)
        real(dp), intent(out) :: fvec(5), jac(5, 5)

        call eval_field(f, x(5), 0.5d0*(x(2) + si%z(2)), 0.5d0*(x(3) + si%z(3)), 2)
        call get_derivatives2(f, 0.5d0*(x(4) + si%z(4)))
        fmid = f
        call eval_field(f, x(1), x(2), x(3), 2)
        call get_derivatives2(f, x(4))
        call midpoint_residual(si%dt, si%pthold, si%z, x, fmid, f, fvec)
        call midpoint_jacobian(si%dt, si%pthold, si%z, x, fmid, f, jac)
    end subroutine midpoint_system


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
    recursive subroutine newton1(si, f, x, maxit, xlast, status)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(inout) :: x(2)
        integer, intent(in) :: maxit
        real(dp), intent(out) :: xlast(2)
        integer, intent(out) :: status

        call newton_solve(SYSTEM_EULER1, 2, si, x, si%atol, si%rtol, maxit, xlast, &
                          status, f=f)
    end subroutine newton1

    ! The damped Newton iteration of all implicit schemes. The scheme selects the
    ! system, the linear solve, the domain guard after the step and the tolerance
    ! references.
    recursive subroutine newton_solve(scheme, n, si, x, atol, rtol, maxit, xlast, &
                                      status, f, fs, s)
        integer, intent(in) :: scheme, n
        type(symplectic_integrator_t), intent(inout) :: si
        real(dp), intent(inout) :: x(n)
        real(dp), intent(in) :: atol, rtol
        integer, intent(in) :: maxit
        real(dp), intent(out) :: xlast(n)
        integer, intent(out) :: status
        type(field_can_t), intent(inout), optional :: f
        type(field_can_t), intent(inout), optional :: fs(:)
        integer, intent(in), optional :: s

        type(newton_state_t) :: state
        ! Fixed size: arrays of size n are heap-allocated by some compilers on
        ! every call. The Jacobian is passed by sequence association.
        real(dp) :: fvec(NEWTON_NMAX), fjac(NEWTON_NMAX*NEWTON_NMAX)
        real(dp) :: flux_scale, aphi
        integer :: stages
        logical :: need_system

        if (n > NEWTON_NMAX) error stop 'newton_solve: system too large'
        stages = 0
        if (present(s)) stages = s
        flux_scale = 0d0
        aphi = 0d0
        if (present(f)) then
            flux_scale = 1d1*torflux/f%ro0
            aphi = f%Aph
        end if
        call newton_begin(state, scheme, n, stages, maxit, atol, rtol, &
                          symplectic_newton_warning_mode, si%warning_factor, &
                          si%z(1), flux_scale, aphi)
        call newton_advance(state, x, xlast, fvec, fjac, need_system)
        do while (need_system)
            select case (scheme)
            case (SYSTEM_EULER1)
                call eval_field(f, x(1), si%z(2), si%z(3), 2)
                call get_derivatives2(f, x(2))
                call sympl_euler1_residual(si, f, x, fvec)
                call sympl_euler1_jacobian(si, f, x, fjac)
            case (SYSTEM_EULER2)
                call f_sympl_euler2(si, f, 3, x, fvec, 1)
                call jac_sympl_euler2(si, f, x, fjac)
            case (SYSTEM_MIDPOINT)
                call midpoint_newton_system(si, f, x, fvec, fjac)
            case (SYSTEM_RK_GAUSS)
                call f_rk_gauss(si, fs, s, x, fvec)
                call jac_rk_gauss(si, fs, s, fjac)
            case (SYSTEM_RK_LOBATTO)
                call f_rk_lobatto(si, fs, s, x, fvec, 2)
                call jac_rk_lobatto(si, fs, s, fjac)
            end select
            call newton_advance(state, x, xlast, fvec, fjac, need_system)
        end do
        status = state%status

        if (.not. state%maxit_reached) return
        select case (scheme)
        case (SYSTEM_EULER1)
            call count_event(EVT_NEWTON1_MAXIT)
        case (SYSTEM_EULER2)
            call count_event(EVT_NEWTON2_MAXIT)
        case (SYSTEM_MIDPOINT)
            call count_event(EVT_MIDPOINT_MAXIT)
        case (SYSTEM_RK_GAUSS)
            call count_event(EVT_RK_GAUSS_MAXIT)
        case (SYSTEM_RK_LOBATTO)
            call count_event(EVT_RK_LOBATTO_MAXIT)
        end select
        call count_warning(state%warning)
    end subroutine newton_solve

    recursive subroutine newton2(si, f, x, atol, rtol, maxit, xlast, status)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(inout) :: x(3)
        real(dp), intent(in) :: atol, rtol
        integer, intent(in) :: maxit
        real(dp), intent(out) :: xlast(3)
        integer, intent(out) :: status

        call newton_solve(SYSTEM_EULER2, 3, si, x, atol, rtol, maxit, xlast, &
                          status, f=f)
    end subroutine newton2

! Small dense Newton solve (n = 2 .. 4*s) by in-house LU with partial pivoting.
! No library call and no allocation: LAPACK/OpenBLAS dgesv per step exhausted
! OpenBLAS memory regions and scaled badly with thread count.
subroutine solve_newton_system(matrix, residual, status)
  real(dp), intent(inout) :: matrix(:, :), residual(:)
  integer, intent(out) :: status

  integer :: info, n

  n = size(residual)
  call rk_solve(n, matrix, residual, info)
  if (info == 0) then
    status = SYMPLECTIC_STEP_OK
  else
    status = SYMPLECTIC_STEP_LINEAR_SOLVE
  end if
end subroutine solve_newton_system

pure logical function boundary_event_converged(fraction_width, radial_residual, &
    fraction_tolerance, radial_tolerance)
  real(dp), intent(in) :: fraction_width, radial_residual
  real(dp), intent(in) :: fraction_tolerance, radial_tolerance

  boundary_event_converged = fraction_width <= fraction_tolerance .and. &
    radial_residual <= radial_tolerance
end function boundary_event_converged

    recursive subroutine newton_midpoint(si, f, x, atol, rtol, maxit, xlast, status)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(inout) :: x(5)  ! = (rend, thend, phend, pphend, rmid)
        real(dp), intent(in) :: atol, rtol
        integer, intent(in) :: maxit
        real(dp), intent(out) :: xlast(5)
        integer, intent(out) :: status

        call newton_solve(SYSTEM_MIDPOINT, 5, si, x, atol, rtol, maxit, xlast, &
                          status, f=f)
    end subroutine newton_midpoint

    recursive subroutine midpoint_newton_system(si, f, x, fvec, jac)
        type(symplectic_integrator_t), intent(in) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: x(5)
        real(dp), intent(out) :: fvec(5), jac(5, 5)

        type(field_can_t) :: fmid

        call midpoint_system(si, f, fmid, x, fvec, jac)
    end subroutine midpoint_newton_system

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! Gauss-Legendre Runge-Kutta method with s internal stages (n=4*s variables)
  !
recursive subroutine f_rk_gauss(si, fs, s, x, fvec)
  ! GC Gauss residual: thin wrapper over the shared core (MODEL_GC). The
  ! arithmetic lives in orbit_rk_core::gauss_canfield_residual.
  type(symplectic_integrator_t), intent(inout) :: si
  integer, intent(in) :: s
  type(field_can_t), intent(inout) :: fs(:)
  real(dp), intent(in) :: x(4*s)  ! = (rend, thend, phend, pphend)
  real(dp), intent(out) :: fvec(4*s)

  call rk_gauss_residual(MODEL_GC, si, fs, s, x, fvec)
  end subroutine f_rk_gauss



  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine jac_rk_gauss(si, fs, s, jac)
  ! GC Gauss Jacobian: thin wrapper over the shared core (MODEL_GC). The
  ! arithmetic lives in orbit_rk_core::gauss_canfield_jacobian.
  type(symplectic_integrator_t), intent(in) :: si
  integer, intent(in) :: s
  type(field_can_t), intent(in) :: fs(:)
  real(dp), intent(out) :: jac(4*s, 4*s)

  call rk_gauss_jacobian(MODEL_GC, si, fs, s, jac)

end subroutine jac_rk_gauss

    recursive subroutine newton_rk_gauss(si, fs, s, x, atol, rtol, maxit, xlast, status)
        type(symplectic_integrator_t), intent(inout) :: si
        integer, intent(in) :: s
        type(field_can_t), intent(inout) :: fs(:)
        real(dp), intent(inout) :: x(4*s)
        real(dp), intent(in) :: atol, rtol
        integer, intent(in) :: maxit
        real(dp), intent(out) :: xlast(4*s)
        integer, intent(out) :: status

        call newton_solve(SYSTEM_RK_GAUSS, 4*s, si, x, atol, rtol, maxit, xlast, &
                          status, fs=fs, s=s)
    end subroutine newton_rk_gauss

recursive subroutine fixpoint_rk_gauss(si, fs, s, x, atol, rtol, maxit, xlast)
  ! TODO: this doesn't work well yet
  type(symplectic_integrator_t), intent(inout) :: si

  integer, intent(in) :: s
  type(field_can_t), intent(inout) :: fs(:)
  integer :: kit, ks

  real(dp), intent(inout) :: x(4*s)
  real(dp), intent(in) :: atol, rtol
  integer, intent(in) :: maxit
  real(dp), intent(out) :: xlast(4*s)

  real(dp) :: fvec(4*s)

  real(dp) :: xabs(4*s), tolref(4*s), fabs(4*s)

  real(dp) :: a(s,s), b(s), c(s), Hprime(s), dHprimedr(s)
  integer :: k, l

  real(dp) :: pthnew, dpthnewdr, damp

  call coeff_rk_gauss(s, a, b, c)  ! TODO: move this to preprocessing

  do kit = 1, 4096

    ! Check if radius left the boundary
    do ks = 1, s
                if (abs(x(4*ks - 3)) > sympl_rmax) return
                ! A negative stage radius is the axis continuation (#370); the
                ! field evaluation maps it to (|s|, theta+pi).
    end do

    call f_rk_gauss(si, fs, s, x, fvec)
    fabs = dabs(fvec)
    xlast = x

    do k = 1, s
      call eval_field(fs(k), x(4*k-3), x(4*k-2), x(4*k-1), 2)
      call get_derivatives2(fs(k), x(4*k))
      Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
      dHprimedr(k) = (fs(k)%d2H(1) - Hprime(k)*fs(k)%d2pth(1))/fs(k)%dpth(1)
    end do

    k = 1
    l = 1

    pthnew = si%pthold - si%dt*a(k,l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))
    dpthnewdr = -si%dt*a(k,l)*( &
    fs(l)%d2H(2) - Hprime(l)*fs(l)%d2pth(2) - dHprimedr(l)*fs(l)%dpth(2))
    !print *, pthnew, dpthnewdr, pthnew/dpthnewdr

    damp = 0.99d0

    x(1) = x(1) - (1d0-damp)*pthnew/dpthnewdr
    x(2) = damp*x(2) + (1d0-damp)*(si%z(2) + si%dt*a(k,l)*si%dt*a(k,l)*Hprime(l))
    x(3) = damp*x(3) + (1d0-damp)*(si%z(3) + si%dt*a(k,l)*(fs(l)%vpar &
                                    - Hprime(l)*fs(l)%hth)/fs(l)%hph)
    x(4) = damp*x(4) + (1d0-damp)*(si%z(4) - si%dt*a(k,l)*(fs(l)%dH(3) &
                                    - Hprime(l)*fs(l)%dpth(3)))

    print *, x

    xabs = dabs(x - xlast)

    do ks = 1, s
      xabs(4*ks-2) = modulo(xabs(4*ks-2), pi)
      xabs(4*ks-1) = modulo(xabs(4*ks-1), pi)
    end do

    do ks = 1, s
      tolref(4*ks-3) = 1d0
      tolref(4*ks-2) = twopi
      tolref(4*ks-1) = twopi
      tolref(4*ks) = dabs(x(4*ks))
    end do

    if (all(fabs < atol)) return
    if (all(xabs < rtol*tolref)) return
  enddo
  call count_event(EVT_FIXPOINT_MAXIT)
end subroutine fixpoint_rk_gauss


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine jac_rk_lobatto(si, fs, s, jac)
  !
  type(symplectic_integrator_t), intent(in) :: si
  integer, intent(in) :: s
  type(field_can_t), intent(in) :: fs(:)
  real(dp), intent(out) :: jac(4*s-2, 4*s-2)

        real(dp) :: a(s, s), ahat(s, s), b(s), c(s)

  call coeff_rk_lobatto(s, a, ahat, b, c)
        ! coeff_rk_lobatto implements s = 3 only.
        if (s /= 3) error stop 'jac_rk_lobatto: only s = 3 is implemented'
        call lobatto3_jacobian(si%dt, a, ahat, fs, jac)
end subroutine jac_rk_lobatto

subroutine guard_lobatto_stage_radii(x, s, status)
  integer, intent(in) :: s
  real(dp), intent(inout) :: x(4*s-2)
  integer, intent(out) :: status

  integer :: ks, radius_index

  status = SYMPLECTIC_STEP_OK
        if (abs(x(1)) > sympl_rmax) then
    status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
    return
  end if

  do ks = 2, s
    radius_index = 4*ks - 5
            if (abs(x(radius_index)) > sympl_rmax) then
      status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
      return
    end if
  end do
end subroutine guard_lobatto_stage_radii


    recursive subroutine newton_rk_lobatto(si, fs, s, x, atol, rtol, maxit, xlast, &
                                           status)
        type(symplectic_integrator_t), intent(inout) :: si
        integer, intent(in) :: s
        type(field_can_t), intent(inout) :: fs(:)
        real(dp), intent(inout) :: x(4*s - 2)
        real(dp), intent(in) :: atol, rtol
        integer, intent(in) :: maxit
        real(dp), intent(out) :: xlast(4*s - 2)
        integer, intent(out) :: status

        call newton_solve(SYSTEM_RK_LOBATTO, 4*s - 2, si, x, atol, rtol, maxit, &
                          xlast, status, fs=fs, s=s)
    end subroutine newton_rk_lobatto

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_timestep_sympl_multi(mi, f, ierr)
  !
  type(multistage_integrator_t), intent(inout) :: mi
  type(field_can_t), intent(inout) :: f

  integer, intent(out) :: ierr

  integer :: kstage
        logical :: chained
        real(dp) :: z_entry(4), pthold_entry
        type(field_can_t) :: f_entry

        ! Each stage starts where the previous one ended, so its history stays
        ! valid across the composition as long as nobody changed the state between
        ! composed steps.
        chained = all(mi%stages(1)%z == mi%stages(2*mi%s)%hist_z)
        do kstage = 1, 2*mi%s
            if (.not. chained) mi%stages(kstage)%nhist = 0
        end do
        if (chained) mi%stages(1)%hist_z = mi%stages(1)%z
        z_entry = mi%stages(1)%z
        pthold_entry = mi%stages(1)%pthold
        f_entry = f
  call orbit_timestep_sympl(mi%stages(1), f, ierr)

  do kstage = 2, 2*mi%s
            ! A failed stage fails the whole composed step, which then restarts
            ! from its entry state; later stages must not run from a stale one.
            if (ierr /= 0) exit
    mi%stages(kstage)%z = mi%stages(kstage-1)%z
            if (chained) mi%stages(kstage)%hist_z = mi%stages(kstage)%z
    mi%stages(kstage)%pthold = f%pth
    call orbit_timestep_sympl(mi%stages(kstage), f, ierr)
  end do
        if (ierr /= 0) then
            mi%stages(1)%z = z_entry
            mi%stages(1)%pthold = pthold_entry
            f = f_entry
            return
        end if

  mi%stages(1)%z = mi%stages(2*mi%s)%z
  mi%stages(1)%pthold = f%pth

end subroutine orbit_timestep_sympl_multi


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
        !
        ! Stage pair k runs euler1 over alpha(k)*dtau, then euler2 over beta(k)*dtau.
        ! This is Hairer, Lubich & Wanner (2006) V.(3.1) with Phi* = euler1, so
        ! alpha holds their beta and beta their alpha.
  !
  type(multistage_integrator_t), intent(inout) :: mi
  type(field_can_t), intent(inout) :: f

  real(dp), intent(in) :: z(:)
  real(dp), intent(in) :: dtau
  integer, intent(in) :: ntau
  real(dp), intent(in) :: rtol_init

  real(dp), intent(in) :: alpha(:), beta(:)

  integer :: ks

  mi%s = size(alpha)

  do ks = 1, mi%s
    call orbit_sympl_init(mi%stages(2*ks-1), f, z, &
      alpha(ks)*dtau, ntau, rtol_init, 1)
    call orbit_sympl_init(mi%stages(2*ks), f, z, &
      beta(ks)*dtau, ntau, rtol_init, 2)
  end do
end subroutine orbit_sympl_init_multi



  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_sympl_init_verlet(mi, f, z, dtau, ntau, rtol_init)
  !
  type(multistage_integrator_t), intent(inout) :: mi

  type(field_can_t), intent(inout) :: f

  real(dp), intent(in) :: z(:)
  real(dp), intent(in) :: dtau
  integer, intent(in) :: ntau
  real(dp), intent(in) :: rtol_init

  real(dp) :: alpha(1), beta(1)

  alpha(1) = 0.5d0
  beta(1)  = 0.5d0

  call orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
end subroutine orbit_sympl_init_verlet


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_sympl_init_order4(mi, f, z, dtau, ntau, rtol_init)
  !
  ! Composition method of order 4 with s=3
  !
  !
    type(multistage_integrator_t), intent(inout) :: mi
    type(field_can_t), intent(inout) :: f

    real(dp), intent(in) :: z(:)
    real(dp), intent(in) :: dtau
    integer, intent(in) :: ntau
    real(dp), intent(in) :: rtol_init

        real(dp) :: alpha(3), beta(3)

    alpha(1) = 1d0/(2d0*(2d0 - 2d0**(1d0/3d0)))
        alpha(2) = -2d0**(1d0/3d0)/(2d0*(2d0 - 2d0**(1d0/3d0)))
    alpha(3) = alpha(1)

    beta = alpha

    call orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
  end subroutine orbit_sympl_init_order4


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_sympl_init_mclachlan4(mi, f, z, dtau, ntau, rtol_init)
  !
  ! Composition method of order 4 with s=5 by McLachlan (1995)
  !
  !
    type(multistage_integrator_t), intent(inout) :: mi
    type(field_can_t), intent(inout) :: f

    real(dp), intent(in) :: z(:)
    real(dp), intent(in) :: dtau
    integer, intent(in) :: ntau
    real(dp), intent(in) :: rtol_init

    real(dp) :: alpha(5), beta(5)

        beta(1) = (146d0 + 5d0*dsqrt(19d0))/540d0
        beta(2) = (-2d0 + 10d0*dsqrt(19d0))/135d0
        beta(3) = 1d0/5d0
        beta(4) = (-23d0 - 20d0*dsqrt(19d0))/270d0
        beta(5) = (14d0 - dsqrt(19d0))/108d0

        alpha(5) = beta(1)
        alpha(4) = beta(2)
        alpha(3) = beta(3)
        alpha(2) = beta(4)
        alpha(1) = beta(5)

    call orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
  end subroutine orbit_sympl_init_mclachlan4


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_sympl_init_blanes4(mi, f, z, dtau, ntau, rtol_init)
  !
  ! Composition method of order 4 with s=6 by Blanes&Moan (2002)
  ! with coefficients in the form of Hairer (2002)
  !
  !
    type(multistage_integrator_t), intent(inout) :: mi
    type(field_can_t), intent(inout) :: f

    real(dp), intent(in) :: z(:)
    real(dp), intent(in) :: dtau
    integer, intent(in) :: ntau
    real(dp), intent(in) :: rtol_init

    real(dp) :: alpha(6), beta(6)

        beta(1) = 0.16231455076687d0
        beta(2) = 0.37087741497958d0
        beta(3) = 0.059762097006575d0
        beta(4) = -0.40993371990193d0
        beta(5) = 0.23399525073150d0
        beta(6) = 0.082984406417405d0

        alpha(6) = beta(1)
        alpha(5) = beta(2)
        alpha(4) = beta(3)
        alpha(3) = beta(4)
        alpha(2) = beta(5)
        alpha(1) = beta(6)

    call orbit_sympl_init_multi(mi, f, z, dtau, ntau, rtol_init, alpha, beta)
  end subroutine orbit_sympl_init_blanes4


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_sympl_init_kahan6(mi, f, z, dtau, ntau, rtol_init)
  !
  ! Composition method of order 6 with s=9 by Kahan&Li (1995)
  ! with coefficients in the form of Hairer (2002)
  !
  !
    type(multistage_integrator_t), intent(inout) :: mi
    type(field_can_t), intent(inout) :: f

    real(dp), intent(in) :: z(:)
    real(dp), intent(in) :: dtau
    integer, intent(in) :: ntau
    real(dp), intent(in) :: rtol_init

    real(dp) :: gam(9)

    gam(1) = 0.39216144400731413927925056d0
    gam(2) = 0.33259913678935943859974864d0
    gam(3) = -0.70624617255763935980996482d0
    gam(4) = 0.08221359629355080023149045d0
    gam(5) = 0.79854399093482996339895035d0
    gam(6) = gam(4)
    gam(7) = gam(3)
    gam(8) = gam(2)
    gam(9) = gam(1)

    call orbit_sympl_init_multi(mi,f,z,dtau,ntau,rtol_init,gam/2d0,gam/2d0)
  end subroutine orbit_sympl_init_kahan6


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  recursive subroutine orbit_sympl_init_kahan8(mi, f, z, dtau, ntau, rtol_init)
    !
    ! Composition method of order 8 with s=17 by Kahan&Li (1995)
    ! with coefficients in the form of Hairer (2002)
    !
    !
      type(multistage_integrator_t), intent(inout) :: mi
      type(field_can_t), intent(inout) :: f

      real(dp), intent(in) :: z(:)
      real(dp), intent(in) :: dtau
      integer, intent(in) :: ntau
      real(dp), intent(in) :: rtol_init

      real(dp) :: gam(17)

      gam(1) =  0.13020248308889008087881763d0
      gam(2) =  0.56116298177510838456196441d0
      gam(3) = -0.38947496264484728640807860d0
      gam(4) =  0.15884190655515560089621075d0
      gam(5) = -0.39590389413323757733623154d0
      gam(6) =  0.18453964097831570709183254d0
      gam(7) =  0.25837438768632204729397911d0
      gam(8) =  0.29501172360931029887096624d0
      gam(9) = -0.60550853383003451169892108d0
      gam(10) = gam(8)
      gam(11) = gam(7)
      gam(12) = gam(6)
      gam(13) = gam(5)
      gam(14) = gam(4)
      gam(15) = gam(3)
      gam(16) = gam(2)
      gam(17) = gam(1)

      call orbit_sympl_init_multi(mi,f,z,dtau,ntau,rtol_init,gam/2d0,gam/2d0)
    end subroutine orbit_sympl_init_kahan8

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_timestep_sympl_expl_impl_euler(si, f, ierr)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  integer, parameter :: n = 2
  integer, parameter :: maxit = 32

  real(dp), dimension(n) :: x, xlast
        real(dp) :: z_start(4), z_guess(4)
  integer :: ktau, newton_status
  logical :: crossed
  type(symplectic_integrator_t) :: accepted_integrator
  type(field_can_t) :: accepted_field

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    accepted_integrator = si
    accepted_field = f
    si%pthold = f%pth

            z_start = si%z
            z_guess = predict_state(si, 1d0)
            x(1) = z_guess(1)
            x(2) = z_guess(4)

    call newton1(si, f, x, maxit, xlast, newton_status)
    if (newton_status /= SYMPLECTIC_STEP_OK) then
      si = accepted_integrator
      f = accepted_field
      ierr = newton_status
      return
    end if

            if (abs(x(1)) > sympl_rmax) then
      si = accepted_integrator
      f = accepted_field
      ierr = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
      return
    end if

    call apply_axis_chart_switch(x(1), si%z(2), crossed, ierr)
    if (ierr /= SYMPLECTIC_STEP_OK) then
      si = accepted_integrator
      f = accepted_field
      return
    end if

    si%z(1) = x(1)
    si%z(4) = x(2)

    if (extrap_field .and. .not. crossed) then
      call sympl_euler1_extrapolate_field(si, f, x, xlast)
    else
      ! After a chart switch xlast lives in the other chart; extrapolation
      ! across the flip is invalid, evaluate the field fresh instead.
      call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
      call get_derivatives(f, si%z(4))
    endif

    call sympl_euler1_advance_angles(si, f)
            call record_history(si, z_start,.not. crossed)

    ktau = ktau+1
  enddo

end subroutine orbit_timestep_sympl_expl_impl_euler


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_timestep_sympl_impl_expl_euler(si, f, ierr)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f

  integer, intent(out) :: ierr

  integer, parameter :: n = 3
  integer, parameter :: maxit = 32

        real(dp), dimension(n) :: x, xlast, dz
        real(dp) :: z_start(4), z_guess(4)
  integer :: ktau, newton_status
  logical :: crossed
  type(symplectic_integrator_t) :: accepted_integrator
  type(field_can_t) :: accepted_field

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    accepted_integrator = si
    accepted_field = f
    si%pthold = f%pth

            z_start = si%z
            z_guess = predict_state(si, 1d0)
            x = z_guess(1:3)

    call newton2(si, f, x, si%atol, si%rtol, maxit, xlast, newton_status)
    if (newton_status /= SYMPLECTIC_STEP_OK) then
      si = accepted_integrator
      f = accepted_field
      ierr = newton_status
      return
    end if

            if (abs(x(1)) > sympl_rmax) then
      si = accepted_integrator
      f = accepted_field
      ierr = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
      return
    end if

    call apply_axis_chart_switch(x(1), x(2), crossed, ierr)
    if (ierr /= SYMPLECTIC_STEP_OK) then
      si = accepted_integrator
      f = accepted_field
      return
    end if

    si%z(1:3) = x

    if (extrap_field .and. .not. crossed) then
      dz(1) = x(1)-xlast(1)
      dz(2) = x(2)-xlast(2)
      dz(3) = x(3)-xlast(3)

      f%dH(1) = f%dH(1) + f%d2H(1)*dz(1) + f%d2H(2)*dz(2) + f%d2H(3)*dz(3)
      f%dH(2) = f%dH(2) + f%d2H(2)*dz(1) + f%d2H(4)*dz(2) + f%d2H(5)*dz(3)
      f%dH(3) = f%dH(3) + f%d2H(3)*dz(1) + f%d2H(5)*dz(2) + f%d2H(6)*dz(3)

      f%dpth(1) = f%dpth(1)+f%d2pth(1)*dz(1)+f%d2pth(2)*dz(2)+f%d2pth(3)*dz(3)
      f%dpth(2) = f%dpth(2)+f%d2pth(2)*dz(1)+f%d2pth(4)*dz(2)+f%d2pth(5)*dz(3)
      f%dpth(3) = f%dpth(3)+f%d2pth(3)*dz(1)+f%d2pth(5)*dz(2)+f%d2pth(6)*dz(3)
    else
      call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
      call get_derivatives(f, si%z(4))
    endif

    f%pth = si%pthold - si%dt*(f%dH(2) - f%dH(1)*f%dpth(2)/f%dpth(1))
    si%z(4) = si%z(4) - si%dt*(f%dH(3) - f%dH(1)*f%dpth(3)/f%dpth(1))
            call record_history(si, z_start,.not. crossed)

    ktau = ktau+1
  enddo

end subroutine orbit_timestep_sympl_impl_expl_euler


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
subroutine attempt_midpoint_step(si, f, status)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: status

        integer, parameter :: n = 5, maxit = 32
        real(dp) :: x(n), xlast(n), z_start(4), z_guess(4)
  logical :: crossed

  ! 32 as in the other schemes: the atol/rtol exits keep converged steps cheap,
  ! and the exact-landing substep solve (#441) needs fully converged trial
  ! steps, which 8 iterations do not guarantee near interfaces.
  si%pthold = f%pth
        z_start = si%z
        x(1:4) = predict_state(si, 1d0)
        z_guess = predict_state(si, 0.5d0)
        x(5) = z_guess(1)

  call newton_midpoint(si, f, x, si%atol, si%rtol, maxit, xlast, status)
  if (status /= SYMPLECTIC_STEP_OK) return

  call apply_axis_chart_switch(x(1), x(2), crossed, status)
  if (status /= SYMPLECTIC_STEP_OK) return

  si%z = x(1:4)
        if (extrap_field .and. .not. crossed) then
            ! The carried momentum uses endpoint coordinates, not x(5).
            ! Old-chart increments cannot extrapolate p_theta after a switch.
    f%pth = f%pth + f%dpth(1)*(x(1)-xlast(1)) &
      + f%dpth(2)*(x(2)-xlast(2)) + f%dpth(3)*(x(3)-xlast(3)) &
      + f%dpth(4)*(x(4)-xlast(4))
  else
    call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
    call get_val(f, si%z(4))
  end if
        call record_history(si, z_start,.not. crossed)
end subroutine attempt_midpoint_step

subroutine get_boundary_event_tolerances(rtol, fraction_tolerance, radial_tolerance)
  real(dp), intent(in) :: rtol
  real(dp), intent(out) :: fraction_tolerance, radial_tolerance

  fraction_tolerance = max(1d-12, 10d0*rtol)
  radial_tolerance = max(1d-10, 10d0*rtol)
  if (boundary_event_fraction_tolerance > 0d0) then
    fraction_tolerance = boundary_event_fraction_tolerance
  end if
  if (boundary_event_radial_tolerance > 0d0) then
    radial_tolerance = boundary_event_radial_tolerance
  end if
end subroutine get_boundary_event_tolerances

subroutine locate_symplectic_boundary(accepted_integrator, accepted_field, &
    raw_step, si, f, status, event_fraction, radial_residual, fraction_width)
  type(symplectic_integrator_t), intent(in) :: accepted_integrator
  type(field_can_t), intent(in) :: accepted_field
  procedure(orbit_timestep_sympl_i) :: raw_step
  type(symplectic_integrator_t), intent(out) :: si
  type(field_can_t), intent(out) :: f
  integer, intent(out) :: status
  real(dp), intent(out) :: event_fraction, radial_residual, fraction_width

  integer, parameter :: max_event_iterations = 64
  type(symplectic_integrator_t) :: trial_integrator, best_integrator
  type(field_can_t) :: trial_field, best_field
  real(dp) :: fraction_low, fraction_ceiling, trial_fraction
  real(dp) :: previous_fraction, previous_radius
  real(dp) :: root_estimate
  real(dp) :: fraction_tolerance, radial_tolerance
  integer :: iteration, trial_status
  logical :: root_estimate_available

  fraction_low = 0d0
  fraction_ceiling = 1d0
  trial_fraction = 0.5d0
  root_estimate_available = .false.
  best_integrator = accepted_integrator
  best_field = accepted_field
  call get_boundary_event_tolerances(accepted_integrator%rtol, &
    fraction_tolerance, radial_tolerance)

  do iteration = 1, max_event_iterations
    trial_integrator = accepted_integrator
    trial_field = accepted_field
    trial_integrator%dt = accepted_integrator%dt*trial_fraction
    trial_integrator%ntau = 1
    call raw_step(trial_integrator, trial_field, trial_status)

    select case (trial_status)
    case (SYMPLECTIC_STEP_OK)
      previous_fraction = fraction_low
      previous_radius = best_integrator%z(1)
      fraction_low = trial_fraction
      best_integrator = trial_integrator
      best_field = trial_field
      radial_residual = abs(sympl_rmax - best_integrator%z(1))

      if (radial_residual == 0d0) then
        fraction_ceiling = fraction_low
        root_estimate = fraction_low
        root_estimate_available = .true.
      else if (best_integrator%z(1) > previous_radius .and. &
          fraction_low > previous_fraction) then
        root_estimate = fraction_low + radial_residual* &
          (fraction_low - previous_fraction)/ &
          (best_integrator%z(1) - previous_radius)
        root_estimate = min(fraction_ceiling, max(fraction_low, root_estimate))
        root_estimate_available = .true.
      else
        root_estimate_available = .false.
      end if
      fraction_width = fraction_ceiling - fraction_low

      if (root_estimate_available .and. &
          boundary_event_converged(fraction_width, radial_residual, &
          fraction_tolerance, radial_tolerance)) then
        si = best_integrator
        f = best_field
        si%dt = accepted_integrator%dt
        si%ntau = accepted_integrator%ntau
        si%z(1) = sympl_rmax
        call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
        call get_val(f, si%z(4))
        event_fraction = fraction_low
        status = SYMPLECTIC_STEP_BOUNDARY
        return
      end if

      if (root_estimate_available .and. root_estimate > fraction_low) then
        trial_fraction = fraction_low + 0.8d0* &
          (root_estimate - fraction_low)
      else
        trial_fraction = 0.5d0*(fraction_low + fraction_ceiling)
      end if
    case (SYMPLECTIC_STEP_OUTSIDE_DOMAIN)
      fraction_ceiling = trial_fraction
      trial_fraction = 0.5d0*(fraction_low + fraction_ceiling)
    case (SYMPLECTIC_STEP_BOUNDARY_LIMITED)
      si = accepted_integrator
      f = accepted_field
      status = SYMPLECTIC_STEP_EVENT_NOT_CONVERGED
      event_fraction = 0d0
      radial_residual = abs(sympl_rmax - best_integrator%z(1))
      fraction_width = fraction_ceiling - fraction_low
      return
    case default
      si = accepted_integrator
      f = accepted_field
      status = trial_status
      event_fraction = 0d0
      radial_residual = abs(sympl_rmax - best_integrator%z(1))
      fraction_width = fraction_ceiling - fraction_low
      return
    end select

    if (trial_fraction <= fraction_low .or. &
        trial_fraction >= fraction_ceiling) then
      exit
    end if
  end do

  si = accepted_integrator
  f = accepted_field
  status = SYMPLECTIC_STEP_EVENT_NOT_CONVERGED
  event_fraction = 0d0
  radial_residual = abs(sympl_rmax - best_integrator%z(1))
  fraction_width = fraction_ceiling - fraction_low
end subroutine locate_symplectic_boundary

recursive subroutine orbit_timestep_sympl_midpoint(si, f, ierr)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  integer, parameter :: n = 5, maxit = 32
        real(dp) :: x(n), xlast(n), z_start(4), z_guess(4)
  integer :: ktau, step_status
  type(symplectic_integrator_t) :: accepted_integrator
  type(field_can_t) :: accepted_field
        logical :: crossed

  if (sympl_rmax > 1d0) then
    ierr = SYMPLECTIC_STEP_OK
    do ktau = 1, si%ntau
      accepted_integrator = si
      accepted_field = f
      si%pthold = f%pth
                z_start = si%z
                x(1:4) = predict_state(si, 1d0)
                z_guess = predict_state(si, 0.5d0)
                x(5) = z_guess(1)
      call newton_midpoint(si, f, x, si%atol, si%rtol, maxit, xlast, &
        step_status)
      if (step_status /= SYMPLECTIC_STEP_OK) then
        si = accepted_integrator
        f = accepted_field
        ierr = step_status
        return
      end if
                if (abs(x(1)) > sympl_rmax) then
        ierr = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
        return
      end if
                call apply_axis_chart_switch(x(1), x(2), crossed, ierr)
                if (ierr /= SYMPLECTIC_STEP_OK) return
      si%z = x(1:4)
                if (extrap_field .and. .not. crossed) then
        ! Only the endpoint coordinates enter the carried momentum.
                    ! Across a chart switch, evaluate p_theta at the new endpoint.
        f%pth = f%pth + f%dpth(1)*(x(1)-xlast(1)) &
          + f%dpth(2)*(x(2)-xlast(2)) + f%dpth(3)*(x(3)-xlast(3)) &
          + f%dpth(4)*(x(4)-xlast(4))
      else
        call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
        call get_val(f, si%z(4))
      end if
                call record_history(si, z_start,.not. crossed)
    end do
    return
  end if

  if (si%ntau /= 1) error stop &
    'orbit_timestep_sympl_midpoint raw step requires ntau=1'
  call attempt_midpoint_step(si, f, ierr)
end subroutine orbit_timestep_sympl_midpoint

    recursive subroutine orbit_timestep_sympl_with_boundary(si, f, ierr)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        call advance_symplectic_with_boundary(si, f, raw_timestep_sympl, ierr)
    end subroutine orbit_timestep_sympl_with_boundary

!> Advance si by one step of the scheme it was initialised with. Dispatching
!> on si%mode keeps integrators of different schemes independent and avoids
!> shared mutable state between OpenMP threads.
    recursive subroutine orbit_timestep_sympl_dispatch(si, f, ierr)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        if (sympl_rmax <= 1d0) then
            call advance_symplectic_with_boundary(si, f, raw_timestep_sympl, ierr)
        else
            call raw_timestep_sympl(si, f, ierr)
        end if
    end subroutine orbit_timestep_sympl_dispatch

    recursive subroutine raw_timestep_sympl(si, f, ierr)
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        select case (si%mode)
        case (EXPL_IMPL_EULER)
            call orbit_timestep_sympl_expl_impl_euler(si, f, ierr)
        case (IMPL_EXPL_EULER)
            call orbit_timestep_sympl_impl_expl_euler(si, f, ierr)
        case (MIDPOINT)
            call orbit_timestep_sympl_midpoint(si, f, ierr)
        case (GAUSS1)
            call orbit_timestep_sympl_gauss1(si, f, ierr)
        case (GAUSS2)
            call orbit_timestep_sympl_gauss2(si, f, ierr)
        case (GAUSS3)
            call orbit_timestep_sympl_gauss3(si, f, ierr)
        case (GAUSS4)
            call orbit_timestep_sympl_gauss4(si, f, ierr)
        case (LOBATTO3)
            call orbit_timestep_sympl_lobatto3(si, f, ierr)
        case (RADAU15)
            call orbit_timestep_sympl_radau15(si, f, ierr)
        case (GBS16)
            call orbit_timestep_sympl_gbs16(si, f, ierr)
        case (TDRK24)
            call orbit_timestep_sympl_tdrk24(si, f, ierr)
        case (TDRK24A)
            call orbit_timestep_sympl_tdrk24a(si, f, ierr)
        case (CASHKARP45)
            call orbit_timestep_sympl_cashkarp45(si, f, ierr)
        case default
            error stop 'symplectic timestep implementation is not initialized'
        end select
    end subroutine raw_timestep_sympl

subroutine advance_symplectic_with_boundary(si, f, raw_step, ierr)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  procedure(orbit_timestep_sympl_i) :: raw_step
  integer, intent(out) :: ierr

  type(symplectic_integrator_t) :: entry_integrator, accepted_integrator
  type(field_can_t) :: entry_field, accepted_field
  real(dp) :: event_fraction, radial_residual, fraction_width
  integer :: ktau, step_status, original_ntau

  entry_integrator = si
  entry_field = f
  original_ntau = si%ntau
  si%last_step_fraction = 0d0
  si%last_event_radial_residual = 0d0
  si%last_event_fraction_width = 0d0
  radial_residual = 0d0
  fraction_width = 0d0
  ierr = SYMPLECTIC_STEP_OK

  if (sympl_rmax > 1d0) then
    do ktau = 1, original_ntau
      accepted_integrator = si
      accepted_field = f
      si%ntau = 1
      call raw_step(si, f, step_status)
      si%ntau = original_ntau
      if (step_status /= SYMPLECTIC_STEP_OK) then
        si = accepted_integrator
        f = accepted_field
        si%last_step_fraction = real(ktau - 1, dp)/real(original_ntau, dp)
        ierr = step_status
        return
      end if
    end do
    si%last_step_fraction = 1d0
    return
  end if

  do ktau = 1, original_ntau
    accepted_integrator = si
    accepted_field = f
    si%ntau = 1
    call raw_step(si, f, step_status)
    si%ntau = original_ntau
    if (step_status == SYMPLECTIC_STEP_OK) then
      if (si%z(1) == sympl_rmax) then
        si%last_step_fraction = real(ktau, dp)/real(original_ntau, dp)
        si%last_event_radial_residual = 0d0
        si%last_event_fraction_width = 0d0
        ierr = SYMPLECTIC_STEP_BOUNDARY
        return
      end if
      cycle
    end if

    if (step_status == SYMPLECTIC_STEP_OUTSIDE_DOMAIN) then
      accepted_integrator%ntau = original_ntau
      call locate_symplectic_boundary(accepted_integrator, accepted_field, &
        raw_step, si, f, step_status, event_fraction, radial_residual, &
        fraction_width)
      if (step_status == SYMPLECTIC_STEP_BOUNDARY) then
        si%last_step_fraction = &
          (real(ktau - 1, dp) + event_fraction)/real(original_ntau, dp)
        si%last_event_radial_residual = radial_residual
        si%last_event_fraction_width = &
          fraction_width/real(original_ntau, dp)
        ierr = step_status
        return
      end if
    end if

    si = entry_integrator
    f = entry_field
    si%last_step_fraction = 0d0
    si%last_event_radial_residual = radial_residual
    si%last_event_fraction_width = fraction_width
    ierr = step_status
    return
  end do

  si%last_step_fraction = 1d0
  si%ntau = original_ntau
end subroutine advance_symplectic_with_boundary


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
    pure subroutine hermite_endpoint(s, nq, c, value, rate, value_end)
        ! Extrapolate smooth functions of the Gauss collocation polynomial u(tau)
        ! to the step end tau = 1 from their values and tau-derivatives at the s
        ! stages. By collocation u'(c_k) = f(u(c_k)), so the stage rates follow
        ! from the stage field caches, and the 2s Hermite conditions give the
        ! endpoint to O(dt**(2s)) per step, the order of the method, with no
        ! field evaluation. Columns of value and rate are the nq quantities.
        integer, intent(in) :: s, nq
        real(dp), intent(in) :: c(s), value(s, nq), rate(s, nq)
        real(dp), intent(out) :: value_end(nq)
        real(dp) :: node(2*s), table(2*s)
        integer :: i, j, n, q

        n = 2*s
        do i = 1, s
            node(2*i - 1) = c(i)
            node(2*i) = c(i)
        end do
        do q = 1, nq
            do i = 1, s
                table(2*i - 1) = value(i, q)
                table(2*i) = value(i, q)
            end do
            ! Divided differences in place; repeated nodes take the derivative.
            do j = 1, n - 1
                do i = n, j + 1, -1
                    if (j == 1 .and. mod(i, 2) == 0) then
                        table(i) = rate(i/2, q)
                    else
                        table(i) = (table(i) - table(i - 1))/(node(i) - node(i - j))
                    end if
                end do
            end do
            value_end(q) = table(n)
            do i = n - 1, 1, -1
                value_end(q) = table(i) + (1d0 - node(i))*value_end(q)
            end do
        end do
    end subroutine hermite_endpoint

recursive subroutine orbit_timestep_sympl_rk_gauss(si, f, s, ierr)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f

  integer, intent(out) :: ierr

  integer, parameter :: maxit = 32

  integer, intent(in) :: s
  real(dp), dimension(4*s) :: x, xlast
  integer :: k, l, ktau, newton_status

  type(field_can_t) :: fs(s)
  type(field_can_t) :: accepted_field
  type(symplectic_integrator_t) :: accepted_integrator
  real(dp) :: a(s,s), b(s), c(s), Hprime(s), dz(4*s)
        real(dp) :: stage_value(s, 6), stage_rate(s, 6), endpoint_value(6)
        real(dp) :: geometry_gradients(3, 5)
        real(dp) :: radius_end, r_rate, theta_rate, phi_rate, pphi_rate, canonical_pth
        real(dp) :: z_start(4)
        logical :: crossed

  do k = 1,s
    fs(k) = f
  end do
        call coeff_rk_gauss(s, a, b, c)

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    accepted_integrator = si
    accepted_field = f
    si%pthold = f%pth

            z_start = si%z
            do k = 1, s
                x(4*k - 3:4*k) = predict_state(si, c(k))
    end do

    call newton_rk_gauss(si, fs, s, x, si%atol, si%rtol, maxit, xlast, &
      newton_status)
    if (newton_status /= SYMPLECTIC_STEP_OK) then
      si = accepted_integrator
      f = accepted_field
      ierr = newton_status
      return
    end if
    !optionally try fixed point iterations, doesn't work yet
    !call fixpoint_rk_gauss(si, fs, s, x, si%atol, si%rtol, maxit, xlast)

    if (any(x(1:4*s:4) > sympl_rmax)) then
      si = accepted_integrator
      f = accepted_field
      ierr = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
      return
    end if

    if (extrap_field) then
      do k = 1, s
        dz(1) = x(4*k-3)-xlast(4*k-3)
        dz(2) = x(4*k-2)-xlast(4*k-2)
        dz(3) = x(4*k-1)-xlast(4*k-1)
        dz(4) = x(4*k)-xlast(4*k)

        fs(k)%pth = fs(k)%pth + fs(k)%dpth(1)*dz(1) &  ! d/dr
                      + fs(k)%dpth(2)*dz(2) &          ! d/dth
                      + fs(k)%dpth(3)*dz(3) &          ! d/dph
                      + fs(k)%dpth(4)*dz(4)            ! d/dpph

        fs(k)%dpth(1) = fs(k)%dpth(1) + fs(k)%d2pth(1)*dz(1) &
                      + fs(k)%d2pth(2)*dz(2) &
                      + fs(k)%d2pth(3)*dz(3) &
                      + fs(k)%d2pth(7)*dz(4)

        fs(k)%dpth(2) = fs(k)%dpth(2) + fs(k)%d2pth(2)*dz(1)&
                      + fs(k)%d2pth(4)*dz(2) &
                      + fs(k)%d2pth(5)*dz(3) &
                      + fs(k)%d2pth(8)*dz(4)

        fs(k)%dpth(3) = fs(k)%dpth(3) + fs(k)%d2pth(3)*dz(1) &
                      + fs(k)%d2pth(5)*dz(2) &
                      + fs(k)%d2pth(6)*dz(3) &
                      + fs(k)%d2pth(9)*dz(4)

        fs(k)%dH(1) = fs(k)%dH(1) + fs(k)%d2H(1)*dz(1) &
                      + fs(k)%d2H(2)*dz(2) &
                      + fs(k)%d2H(3)*dz(3) &
                      + fs(k)%d2H(7)*dz(4)

        fs(k)%dH(2) = fs(k)%dH(2) + fs(k)%d2H(2)*dz(1) &
                      + fs(k)%d2H(4)*dz(2) &
                      + fs(k)%d2H(5)*dz(3) &
                      + fs(k)%d2H(8)*dz(4)

        fs(k)%dH(3) = fs(k)%dH(3) + fs(k)%d2H(3)*dz(1) &
                      + fs(k)%d2H(5)*dz(2) &
                      + fs(k)%d2H(6)*dz(3) &
                      + fs(k)%d2H(9)*dz(4)

        fs(k)%vpar = fs(k)%vpar + fs(k)%dvpar(1)*dz(1) &  ! d/dr
                   + fs(k)%dvpar(2)*dz(2) &
                   + fs(k)%dvpar(3)*dz(3) &
                   + fs(k)%dvpar(4)*dz(4)

        fs(k)%hth = fs(k)%hth + fs(k)%dhth(1)*dz(1) &  ! d/dr
                      + fs(k)%dhth(2)*dz(2) &
                      + fs(k)%dhth(3)*dz(3)

        fs(k)%hph = fs(k)%hph + fs(k)%dhph(1)*dz(1) &  ! d/dr
                      + fs(k)%dhph(2)*dz(2) &
                      + fs(k)%dhph(3)*dz(3)

                    ! Not used by the update, but extrapolated to the endpoint cache.
                    fs(k)%Ath = fs(k)%Ath + fs(k)%dAth(1)*dz(1) + fs(k)%dAth(2)*dz(2) &
                                + fs(k)%dAth(3)*dz(3)
                    fs(k)%Aph = fs(k)%Aph + fs(k)%dAph(1)*dz(1) + fs(k)%dAph(2)*dz(2) &
                                + fs(k)%dAph(3)*dz(3)
                    fs(k)%Bmod = fs(k)%Bmod + fs(k)%dBmod(1)*dz(1) &
                                 + fs(k)%dBmod(2)*dz(2) + fs(k)%dBmod(3)*dz(3)
                    ! Hermite chain rates need gradients at the corrected stage.
                    call gc_geometry_gradient_transport(dz(1:3), fs(k), &
                                                        geometry_gradients)
                    fs(k)%dAth = geometry_gradients(:, 1)
                    fs(k)%dAph = geometry_gradients(:, 2)
                    fs(k)%dhth = geometry_gradients(:, 3)
                    fs(k)%dhph = geometry_gradients(:, 4)
                    fs(k)%dBmod = geometry_gradients(:, 5)
      end do
    else
      do k = 1, s
        call eval_field(fs(k), x(4*k-3), x(4*k-2), x(4*k-1), 0)
        call get_derivatives(fs(k), x(4*k))
      end do
    end if

    do k = 1, s
      Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
    end do

    f = fs(s)
    f%pth = si%pthold
    si%z(1) = x(4*s-3)

    do l = 1, s
      f%pth = f%pth - si%dt*b(l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))
      si%z(2) = si%z(2) + si%dt*b(l)*Hprime(l)
      si%z(3) = si%z(3)+si%dt*b(l)*(fs(l)%vpar-Hprime(l)*fs(l)%hth)/fs(l)%hph
      si%z(4) = si%z(4) - si%dt*b(l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))
    end do

            ! The update advances (theta, phi, p_theta, p_phi); the radius is
            ! defined only implicitly by p_theta. Report the endpoint radius and
            ! field cache by Hermite extrapolation of the stages. The error does not
            ! accumulate: each step solves its stage radii from the constraint
            ! again, and the endpoint only seeds the next Newton guess.
            ! Stage values and rates in tau = t/dt of r and of the cached field
            ! values callers read after the step (vpar and H follow from p_phi).
            ! dr/dt follows from d/dt p_theta(r, theta, phi, p_phi) = dp_theta/dt
            ! along the stage; the theta terms cancel against p_theta's own rate.
            do k = 1, s
                theta_rate = si%dt*Hprime(k)
                phi_rate = si%dt*(fs(k)%vpar - Hprime(k)*fs(k)%hth)/fs(k)%hph
                pphi_rate = -si%dt*(fs(k)%dH(3) - Hprime(k)*fs(k)%dpth(3))
                r_rate = (-si%dt*fs(k)%dH(2) - fs(k)%dpth(3)*phi_rate &
                          - fs(k)%hth/fs(k)%hph*pphi_rate)/fs(k)%dpth(1)
                stage_value(k, 1) = x(4*k - 3)
                stage_value(k, 2) = fs(k)%Ath
                stage_value(k, 3) = fs(k)%Aph
                stage_value(k, 4) = fs(k)%hth
                stage_value(k, 5) = fs(k)%hph
                stage_value(k, 6) = fs(k)%Bmod
                stage_rate(k, 1) = r_rate
                stage_rate(k, 2) = chain_rate(fs(k)%dAth)
                stage_rate(k, 3) = chain_rate(fs(k)%dAph)
                stage_rate(k, 4) = chain_rate(fs(k)%dhth)
                stage_rate(k, 5) = chain_rate(fs(k)%dhph)
                stage_rate(k, 6) = chain_rate(fs(k)%dBmod)
            end do
            call hermite_endpoint(s, 6, c, stage_value, stage_rate, endpoint_value)
            radius_end = endpoint_value(1)
            if (.not. finite_newton_iterate([radius_end])) then
                si = accepted_integrator
                f = accepted_field
                ierr = SYMPLECTIC_STEP_MAXITER
                return
            end if
            if (radius_end > sympl_rmax) then
                ! A crossing inside the last stage interval; the boundary event solve
                ! locates it on the step fraction.
                si = accepted_integrator
                f = accepted_field
                ierr = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
                return
            end if
            ! The polar chart cannot represent a negative radius. Near the axis keep
            ! the last stage radius and cache, as before; the next step solves its
            ! stages anew.
            if (radius_end >= 0d0) then
                si%z(1) = radius_end
                canonical_pth = f%pth
                f%Ath = endpoint_value(2)
                f%Aph = endpoint_value(3)
                f%hth = endpoint_value(4)
                f%hph = endpoint_value(5)
                f%Bmod = endpoint_value(6)
                call get_val(f, si%z(4))
                f%pth = canonical_pth
            end if
            ! A stage radius past the axis is the point (|r|, theta+pi); commit that
            ! chart switch as the Euler and midpoint steps do (#370).
            call apply_axis_chart_switch(si%z(1), si%z(2), crossed, ierr)
            if (ierr /= SYMPLECTIC_STEP_OK) then
                si = accepted_integrator
                f = accepted_field
                return
            end if
            if (crossed) then
                ! The stage cache lives in the other chart: evaluate the field
                ! and p_theta afresh at the new endpoint, as Euler and midpoint do.
                call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
                call get_val(f, si%z(4))
            end if
            call record_history(si, z_start,.not. crossed)

    ktau = ktau+1
  enddo

    contains

        pure real(dp) function chain_rate(gradient)
            ! d/dtau of a field value along the stage from its (r, theta, phi) gradient.
            real(dp), intent(in) :: gradient(3)

            chain_rate = gradient(1)*r_rate + gradient(2)*theta_rate &
                         + gradient(3)*phi_rate
        end function chain_rate

end subroutine orbit_timestep_sympl_rk_gauss


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
recursive subroutine orbit_timestep_sympl_rk_lobatto(si, f, s, ierr)
  !
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), intent(inout) :: f
  integer, intent(out) :: ierr

  integer, parameter :: maxit = 32

  integer, intent(in) :: s
  real(dp), dimension(4*s-2) :: x, xlast
  integer :: ktau, k, l, newton_status

  type(field_can_t) :: fs(s)
  type(field_can_t) :: accepted_field
  type(symplectic_integrator_t) :: accepted_integrator
  real(dp) :: a(s,s), ahat(s,s), b(s), c(s), Hprime(s)
        real(dp) :: z_start(4)
        logical :: crossed

  do k = 1,s
    fs(k) = f
  end do
        call coeff_rk_lobatto(s, a, ahat, b, c)

  ierr = 0
  ktau = 0
  do while(ktau .lt. si%ntau)
    accepted_integrator = si
    accepted_field = f
    si%pthold = f%pth

            z_start = si%z
    x(1) = si%z(1)
    x(2) = si%z(4)
    do k = 2,s
                x((4*k - 3 - 2):(4*k - 2)) = predict_state(si, c(k))
    end do

    call newton_rk_lobatto(si, fs, s, x, si%atol, si%rtol, maxit, xlast, &
      newton_status)
    if (newton_status /= SYMPLECTIC_STEP_OK) then
      si = accepted_integrator
      f = accepted_field
      ierr = newton_status
      return
    end if

    call eval_field(fs(1), x(1), si%z(2), si%z(3), 0)
    call get_derivatives(fs(1), x(2))
    Hprime(1) = fs(1)%dH(1)/fs(1)%dpth(1)

    do k = 2, s
      call eval_field(fs(k), x(4*k-3-2), x(4*k-2-2), x(4*k-1-2), 0)
      call get_derivatives(fs(k), x(4*k-2))
      Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
    end do

    f = fs(s)
    f%pth   = si%pthold
    si%z(1) = x(4*s-3-2)

    do l = 1, s
      f%pth   = f%pth   - si%dt*b(l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))
      si%z(2) = si%z(2) + si%dt*b(l)*Hprime(l)
      si%z(3) = si%z(3) + si%dt*b(l)*(fs(l)%vpar-Hprime(l)*fs(l)%hth)/fs(l)%hph
      si%z(4) = si%z(4) - si%dt*b(l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))
    end do
            ! A stage radius past the axis is the point (|r|, theta+pi); commit that
            ! chart switch as the Euler and midpoint steps do (#370).
            call apply_axis_chart_switch(si%z(1), si%z(2), crossed, ierr)
            if (ierr /= SYMPLECTIC_STEP_OK) then
                si = accepted_integrator
                f = accepted_field
                return
            end if
            if (crossed) then
                ! The stage cache lives in the other chart: evaluate the field
                ! and p_theta afresh at the new endpoint, as Euler and midpoint do.
                call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
                call get_val(f, si%z(4))
            end if
            call record_history(si, z_start,.not. crossed)

    ktau = ktau+1
  enddo

end subroutine orbit_timestep_sympl_rk_lobatto


recursive subroutine debug_root(si, f, x0)
  type(symplectic_integrator_t), intent(inout) :: si
  type(field_can_t), pointer :: f

  real(dp) :: x0(2), x(2)
  integer :: k, l, iflag
  integer, parameter :: n = 100
  real(dp), parameter :: eps = 1d-15

  real(dp) :: fvec(2)

  do k = -n,n
    do l = -n,n
      x = x0 + l*eps/n*(/x0(1),0d0/) + k*eps/n*(/0d0,x0(2)/)
      call f_sympl_euler1(si, f, 2, x, fvec, iflag)
      write(5001,*) x(1), x(2), fvec
    end do
  end do

end subroutine debug_root

end module orbit_symplectic
