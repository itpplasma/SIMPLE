module orbit_newton
    ! Damped Newton iteration of the implicit symplectic steps, shared by the CPU
    ! and the OpenACC steppers. The driver asks its caller for the residual and
    ! Jacobian at each iterate (reverse communication), so it needs no procedure
    ! argument and runs on the device; the caller evaluates the field there.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use util, only: twopi
    use orbit_symplectic_base, only: sympl_rmax, rate_converged, &
                                     SYMPLECTIC_STEP_OK, &
                                     SYMPLECTIC_STEP_OUTSIDE_DOMAIN, &
                                     SYMPLECTIC_STEP_MAXITER, &
                                     SYMPLECTIC_STEP_LINEAR_SOLVE, &
                                     SYMPLECTIC_STEP_BOUNDARY_LIMITED
    use orbit_symplectic_euler1, only: sympl_euler1_solve
    use orbit_rk_core, only: rk_solve

    implicit none
    private

    ! Implicit systems.
    integer, parameter, public :: SYSTEM_EULER1 = 1, SYSTEM_EULER2 = 2, &
                                  SYSTEM_MIDPOINT = 3, SYSTEM_RK_GAUSS = 4, &
                                  SYSTEM_RK_LOBATTO = 5
    ! Largest system: four Gauss stages of four unknowns each.
    integer, parameter, public :: NEWTON_NMAX = 16

    ! How warning mode decided after the last iteration; the CPU counts these.
    integer, parameter, public :: WARNING_NONE = 0, WARNING_ACCEPT = 1, &
                                  WARNING_REJECT_100 = 2, &
                                  WARNING_REJECT_1E4 = 3, &
                                  WARNING_REJECT_LARGE = 4, &
                                  WARNING_NONFINITE = 5, &
                                  WARNING_AXIS_ACCEPT = 6

    type, public :: newton_state_t
        integer :: system = 0, n = 0, s = 0, nradial = 0
        integer :: kit = 0, maxit = 0
        integer :: status = SYMPLECTIC_STEP_MAXITER
        integer :: warning = WARNING_NONE
        integer :: radial_indices(NEWTON_NMAX) = 0
        real(dp) :: atol = 0d0, rtol = 0d0, warning_factor = 0d0
        real(dp) :: radius_start = 0d0, dnorm_prev = 0d0
        real(dp) :: tolref(NEWTON_NMAX) = 1d0
        logical :: warning_mode = .false.
        logical :: boundary_limited = .false., step_limited = .false.
        logical :: maxit_reached = .false.
    end type newton_state_t

    public :: newton_begin, newton_advance
    public :: finite_newton_iterate, warning_maxiter_class, axis_crossing_converged
    public :: solve_euler2_system, matrix3_near_singular
    public :: limit_radial_newton_step, radial_boundary_reached

contains

    subroutine newton_begin(state, system, n, s, maxit, atol, rtol, &
                            warning_mode, warning_factor, radius_start, &
                            flux_scale, aphi)
        ! flux_scale = 10 torflux/ro0 and aphi = A_phi at the start set the
        ! tolerance reference of the canonical momentum.
        !$acc routine seq
        type(newton_state_t), intent(out) :: state
        integer, intent(in) :: system, n, s, maxit
        real(dp), intent(in) :: atol, rtol, warning_factor, radius_start
        real(dp), intent(in) :: flux_scale, aphi
        logical, intent(in) :: warning_mode

        integer :: ks

        state%system = system
        state%n = n
        state%s = s
        state%maxit = maxit
        state%atol = atol
        state%rtol = rtol
        state%warning_mode = warning_mode
        state%warning_factor = warning_factor
        state%radius_start = radius_start
        state%kit = 0
        state%status = SYMPLECTIC_STEP_MAXITER
        state%warning = WARNING_NONE
        state%dnorm_prev = 0d0
        state%boundary_limited = .false.
        state%step_limited = .false.
        state%maxit_reached = .false.
        state%tolref = 1d0

        select case (system)
        case (SYSTEM_EULER1)
            state%nradial = 1
            state%radial_indices(1) = 1
            state%tolref(2) = dabs(flux_scale)
        case (SYSTEM_EULER2)
            state%nradial = 1
            state%radial_indices(1) = 1
            state%tolref(2:3) = twopi
        case (SYSTEM_MIDPOINT)
            state%nradial = 2
            state%radial_indices(1) = 1
            state%radial_indices(2) = 5
            state%tolref(2:3) = twopi
            state%tolref(4) = max(dabs(aphi), dabs(flux_scale))
        case (SYSTEM_RK_GAUSS)
            state%nradial = s
            do ks = 1, s
                state%radial_indices(ks) = 4*ks - 3
            end do
        case (SYSTEM_RK_LOBATTO)
            state%nradial = s
            state%radial_indices(1) = 1
            do ks = 2, s
                state%radial_indices(ks) = 4*ks - 5
            end do
        end select
    end subroutine newton_begin

    subroutine newton_advance(state, x, xlast, fvec, fjac, need_system)
        ! Call once after newton_begin and then after each evaluation of the
        ! residual fvec and Jacobian fjac at x, while need_system is true. On
        ! return with need_system false, state%status holds the result.
        !$acc routine seq
        type(newton_state_t), intent(inout) :: state
        real(dp), intent(inout) :: x(state%n)
        real(dp), intent(inout) :: xlast(state%n)
        real(dp), intent(in) :: fvec(state%n)
        real(dp), intent(inout) :: fjac(state%n*state%n)
        logical, intent(out) :: need_system

        need_system = .false.
        if (state%kit == 0) then
            xlast = x
        else
            call newton_correct(state, x, xlast, fvec, fjac)
            if (state%status /= SYMPLECTIC_STEP_MAXITER) return
        end if
        if (state%kit >= state%maxit) then
            call newton_finish(state, x, xlast)
            return
        end if
        ! A negative radial label is the axis continuation (#370): the field
        ! evaluation maps it to (|s|, theta+pi), and the caller commits the chart
        ! switch on the converged solution.
        if (radius_beyond(x, state%radial_indices(1:state%nradial), .true.)) then
            state%status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
            return
        end if
        state%kit = state%kit + 1
        need_system = .true.
    end subroutine newton_advance

    subroutine newton_correct(state, x, xlast, fvec, fjac)
        ! One Newton correction from the residual and Jacobian at x.
        !$acc routine seq
        type(newton_state_t), intent(inout) :: state
        real(dp), intent(inout) :: x(state%n)
        real(dp), intent(out) :: xlast(state%n)
        real(dp), intent(in) :: fvec(state%n)
        real(dp), intent(inout) :: fjac(state%n*state%n)

        ! Fixed-size work arrays: arrays of size n are heap-allocated by some
        ! compilers on every call, which tripled the CPU cost of a step.
        real(dp) :: correction(NEWTON_NMAX), fabs(NEWTON_NMAX)
        real(dp) :: xabs(NEWTON_NMAX), tolerance(NEWTON_NMAX)
        real(dp) :: step_scale
        integer :: n, ks, info
        logical :: linear_failed, nonnegative

        n = state%n
        fabs(1:n) = dabs(fvec)
        xlast = x

        ! The Euler schemes keep their closed-form inverses: rk_solve would
        ! round differently and change the iterates.
        select case (state%system)
        case (SYSTEM_EULER1)
            call sympl_euler1_solve(fjac, fvec, correction, linear_failed)
        case (SYSTEM_EULER2)
            call solve_euler2_system(fjac, fvec, correction, linear_failed)
        case default
            correction(1:n) = fvec
            call rk_solve(n, fjac, correction, info)
            linear_failed = info /= 0
        end select
        if (linear_failed) then
            state%status = SYMPLECTIC_STEP_LINEAR_SOLVE
            return
        end if

        call limit_radial_newton_step(x, correction(1:n), &
                                      state%radial_indices(1:state%nradial), &
                                      step_scale, state%step_limited)
        state%boundary_limited = state%boundary_limited .or. state%step_limited
        x = x - step_scale*correction(1:n)

        ! Midpoint and Euler-1 bound |r| like the guard before the step, the
        ! Runge-Kutta schemes only r itself, and Euler-2 leaves it to the next
        ! guard.
        select case (state%system)
        case (SYSTEM_EULER1, SYSTEM_MIDPOINT)
            if (radius_beyond(x, state%radial_indices(1:state%nradial), .true.)) then
                state%status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
                return
            end if
        case (SYSTEM_RK_GAUSS, SYSTEM_RK_LOBATTO)
            if (radius_beyond(x, state%radial_indices(1:state%nradial), .false.)) then
                state%status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
                return
            end if
        end select
        xabs(1:n) = dabs(x - xlast)

        select case (state%system)
        case (SYSTEM_EULER1)
            state%tolref(2) = max(dabs(x(2)), state%tolref(2))
        case (SYSTEM_MIDPOINT)
            ! Don't take too small values in pphi as tolerance reference
            state%tolref(4) = max(dabs(x(4)), state%tolref(4))
        case (SYSTEM_RK_GAUSS)
            do ks = 1, state%s
                state%tolref(4*ks - 3) = 1d0
                state%tolref(4*ks - 2) = twopi
                state%tolref(4*ks - 1) = twopi
                state%tolref(4*ks) = dabs(xlast(4*ks))
            end do
        case (SYSTEM_RK_LOBATTO)
            state%tolref(1) = 1d0
            state%tolref(2) = dabs(xlast(2))
            do ks = 2, state%s
                state%tolref(4*ks - 2 - 3) = 1d0
                state%tolref(4*ks - 2 - 2) = twopi
                state%tolref(4*ks - 2 - 1) = twopi
                state%tolref(4*ks - 2) = dabs(xlast(4*ks - 2))
            end do
        end select
        tolerance(1:n) = state%rtol*state%tolref(1:n)

        ! Gauss stages accept only nonnegative radii.
        nonnegative = .true.
        if (state%system == SYSTEM_RK_GAUSS) then
            do ks = 1, state%nradial
                nonnegative = nonnegative .and. &
                              x(state%radial_indices(ks)) >= 0d0
            end do
        end if

        if (.not. state%step_limited .and. nonnegative .and. &
            all(fabs(1:n) < state%atol)) then
            state%status = SYMPLECTIC_STEP_OK
            return
        end if
        if (.not. state%step_limited .and. nonnegative .and. &
            all(xabs(1:n) < tolerance(1:n))) then
            state%status = SYMPLECTIC_STEP_OK
            return
        end if
        if (accept_by_rate(xabs, tolerance, n, state%step_limited, &
                           state%dnorm_prev) .and. nonnegative) then
            state%status = SYMPLECTIC_STEP_OK
            return
        end if
    end subroutine newton_correct

    subroutine newton_finish(state, x, xlast)
        ! Status after maxit iterations without convergence.
        !$acc routine seq
        type(newton_state_t), intent(inout) :: state
        real(dp), intent(in) :: x(state%n), xlast(state%n)

        integer :: nr
        logical :: at_boundary

        nr = state%nradial
        if (state%boundary_limited) then
            at_boundary = radial_boundary_reached(x, state%radial_indices(1:nr), &
                                                  state%rtol)
            if (state%step_limited .and. at_boundary) then
                state%status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
            else
                state%status = SYMPLECTIC_STEP_BOUNDARY_LIMITED
            end if
            return
        end if
        state%maxit_reached = .true.
        if (state%system == SYSTEM_MIDPOINT) then
            if (axis_crossing_converged(state%warning_mode, state%radius_start, &
                                        x, xlast, state%tolref, state%rtol, &
                                        state%warning_factor)) then
                state%warning = WARNING_AXIS_ACCEPT
                state%status = SYMPLECTIC_STEP_OK
                return
            end if
        end if
        state%warning = warning_maxiter_class(state%warning_mode, x, xlast, &
                                              state%tolref(1:state%n), &
                                              state%rtol, state%warning_factor)
        if (state%warning == WARNING_ACCEPT) state%status = SYMPLECTIC_STEP_OK
    end subroutine newton_finish

    logical function accept_by_rate(xabs, tol, n, step_limited, dnorm_prev)
        ! Rate-based acceptance of the current Newton iterate; dnorm_prev carries
        ! the scaled norm of the previous correction between iterations.
        !$acc routine seq
        integer, intent(in) :: n
        real(dp), intent(in) :: xabs(n), tol(n)
        logical, intent(in) :: step_limited
        real(dp), intent(inout) :: dnorm_prev
        real(dp) :: dnorm

        accept_by_rate = .false.
        ! A zero tolerance demands the correction rule; no estimate meets it.
        if (step_limited .or. .not. all(tol > 0d0)) then
            dnorm_prev = 0d0
            return
        end if
        dnorm = maxval(xabs/tol)
        accept_by_rate = rate_converged(dnorm, dnorm_prev)
        dnorm_prev = dnorm
    end function accept_by_rate

    pure logical function radius_beyond(x, radial_indices, absolute)
        !$acc routine seq
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: radial_indices(:)
        logical, intent(in) :: absolute

        integer :: k

        radius_beyond = .false.
        do k = 1, size(radial_indices)
            if (absolute) then
                radius_beyond = abs(x(radial_indices(k))) > sympl_rmax
            else
                radius_beyond = x(radial_indices(k)) > sympl_rmax
            end if
            if (radius_beyond) return
        end do
    end function radius_beyond

    pure subroutine limit_radial_newton_step(x, correction, radial_indices, &
                                             step_scale, limited)
        !$acc routine seq
        real(dp), intent(in) :: x(:), correction(:)
        integer, intent(in) :: radial_indices(:)
        real(dp), intent(out) :: step_scale
        logical, intent(out) :: limited

        real(dp) :: outward_step
        integer :: k, radius_index

        step_scale = 1d0
        if (sympl_rmax > 1d0) then
            limited = .false.
            return
        end if
        do k = 1, size(radial_indices)
            radius_index = radial_indices(k)
            outward_step = -correction(radius_index)
            if (outward_step > 0d0) then
                step_scale = min(step_scale, &
                                 0.8d0*max(0d0, sympl_rmax - x(radius_index))/ &
                                 outward_step)
            end if
        end do
        limited = step_scale < 1d0
    end subroutine limit_radial_newton_step

    pure logical function radial_boundary_reached(x, radial_indices, rtol)
        !$acc routine seq
        real(dp), intent(in) :: x(:), rtol
        integer, intent(in) :: radial_indices(:)

        integer :: k

        radial_boundary_reached = .false.
        do k = 1, size(radial_indices)
            if (abs(sympl_rmax - x(radial_indices(k))) <= max(1d-10, 10d0*rtol)) then
                radial_boundary_reached = .true.
                return
            end if
        end do
    end function radial_boundary_reached

    pure logical function finite_newton_iterate(x)
        !$acc routine seq
        real(dp), intent(in) :: x(:)

        integer(int64), parameter :: exponent_mask = &
                                     int(z'7FF0000000000000', int64)
        integer(int64) :: bits
        integer :: i

        finite_newton_iterate = .false.
        do i = 1, size(x)
            bits = transfer(x(i), bits)
            if (iand(bits, exponent_mask) == exponent_mask) return
        end do
        finite_newton_iterate = .true.
    end function finite_newton_iterate

    pure integer function warning_maxiter_class(warning_mode, x, xlast, tolref, &
                                                rtol, warning_factor)
        ! Warning mode continues only when the final Newton correction is finite
        ! and inside this integrator's tolerance-scaled guard. Other iterates
        ! enter retry recovery.
        !$acc routine seq
        logical, intent(in) :: warning_mode
        real(dp), intent(in) :: x(:), xlast(:), tolref(:), rtol, warning_factor
        real(dp) :: correction_ratio

        warning_maxiter_class = WARNING_NONE
        if (.not. warning_mode .or. rtol <= 0.0_dp) return
        if (.not. finite_newton_iterate(x) .or. &
            .not. finite_newton_iterate(xlast) .or. &
            .not. finite_newton_iterate(tolref)) then
            warning_maxiter_class = WARNING_NONFINITE
            return
        end if
        correction_ratio = maxval(abs(x - xlast)/ &
                                  max(rtol*abs(tolref), tiny(1.0_dp)))
        if (correction_ratio <= warning_factor) then
            warning_maxiter_class = WARNING_ACCEPT
        else if (correction_ratio <= 100.0_dp) then
            warning_maxiter_class = WARNING_REJECT_100
        else if (correction_ratio <= 1.0e4_dp) then
            warning_maxiter_class = WARNING_REJECT_1E4
        else
            warning_maxiter_class = WARNING_REJECT_LARGE
        end if
    end function warning_maxiter_class

    pure logical function axis_crossing_converged(warning_mode, radius_start, x, &
                                                  xlast, tolref, rtol, &
                                                  warning_factor)
        ! Near the polar-coordinate axis, the endpoint and midpoint radii can
        ! change sign while the physical step and its angular/momentum components
        ! have converged. Accept only that local chart crossing; the caller
        ! immediately maps the negative endpoint to the equivalent positive-radius
        ! chart.
        !$acc routine seq
        logical, intent(in) :: warning_mode
        real(dp), intent(in) :: radius_start, x(5), xlast(5), tolref(5)
        real(dp), intent(in) :: rtol, warning_factor
        real(dp), parameter :: axis_radius = 1.0e-6_dp

        axis_crossing_converged = .false.
        if (.not. warning_mode .or. rtol <= 0.0_dp) return
        if (radius_start < 0.0_dp .or. radius_start > axis_radius) return
        if (.not. finite_newton_iterate(x)) return
        if (.not. finite_newton_iterate(xlast)) return
        if (.not. finite_newton_iterate(tolref)) return
        if (x(1) >= 0.0_dp) return
        if (max(abs(x(1)), abs(x(5))) > axis_radius) return
        axis_crossing_converged = all(abs(x(2:4) - xlast(2:4)) <= &
                                      warning_factor*rtol*abs(tolref(2:4)))
    end function axis_crossing_converged

    pure subroutine solve_euler2_system(fjac, fvec, correction, failed)
        !$acc routine seq
        real(dp), intent(in) :: fjac(3, 3), fvec(3)
        real(dp), intent(out) :: correction(3)
        logical, intent(out) :: failed

        real(dp) :: jinv(3, 3), det

        det = fjac(1, 1)*fjac(2, 2)*fjac(3, 3) + fjac(1, 2)*fjac(2, 3)*fjac(3, 1) &
              + fjac(1, 3)*fjac(2, 1)*fjac(3, 2) - fjac(1, 3)*fjac(2, 2)*fjac(3, 1) &
              - fjac(1, 1)*fjac(2, 3)*fjac(3, 2) - fjac(1, 2)*fjac(2, 1)*fjac(3, 3)
        failed = matrix3_near_singular(fjac, det)
        if (failed) then
            correction = 0d0
            return
        end if

        jinv(1, 1) = fjac(2, 2)*fjac(3, 3) - fjac(2, 3)*fjac(3, 2)
        jinv(1, 2) = fjac(1, 3)*fjac(3, 2) - fjac(1, 2)*fjac(3, 3)
        jinv(1, 3) = fjac(1, 2)*fjac(2, 3) - fjac(1, 3)*fjac(2, 2)

        jinv(2, 1) = fjac(2, 3)*fjac(3, 1) - fjac(2, 1)*fjac(3, 3)
        jinv(2, 2) = fjac(1, 1)*fjac(3, 3) - fjac(1, 3)*fjac(3, 1)
        jinv(2, 3) = fjac(1, 3)*fjac(2, 1) - fjac(1, 1)*fjac(2, 3)

        jinv(3, 1) = fjac(2, 1)*fjac(3, 2) - fjac(3, 1)*fjac(2, 2)
        jinv(3, 2) = fjac(1, 2)*fjac(3, 1) - fjac(3, 2)*fjac(1, 1)
        jinv(3, 3) = fjac(1, 1)*fjac(2, 2) - fjac(1, 2)*fjac(2, 1)

        correction = matmul(jinv, fvec)/det
    end subroutine solve_euler2_system

    pure logical function matrix3_near_singular(matrix, determinant)
        !$acc routine seq
        real(dp), intent(in) :: matrix(3, 3), determinant
        real(dp) :: matrix_scale

        matrix_scale = maxval(abs(matrix))
        matrix3_near_singular = matrix_scale == 0d0 .or. &
                                abs(determinant) <= epsilon(1d0)*matrix_scale**3
    end function matrix3_near_singular

end module orbit_newton
