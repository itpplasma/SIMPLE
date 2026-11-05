module orbit_quasi_implicit
! Quasi-Newton (fortnum multiroot_hybrids) versions of the implicit symplectic
! steps. simple.x uses the Newton steppers of orbit_symplectic; these remain
! only as an independent root solver for the tests that compare against them.

    use util, only: pi, twopi
    use field_can_mod, only: eval_field => evaluate, field_can_t, get_derivatives
    use orbit_symplectic_base, only: symplectic_integrator_t, multistage_integrator_t, &
   orbit_timestep_quasi_i, coeff_rk_gauss, coeff_rk_lobatto, f_rk_lobatto, sympl_rmax, &
            SYMPLECTIC_STEP_OK, SYMPLECTIC_STEP_MAXITER, SYMPLECTIC_STEP_LINEAR_SOLVE, &
              symplectic_newton_warning_mode, EXPL_IMPL_EULER, IMPL_EXPL_EULER, S_MAX, &
             RK45, MIDPOINT, GAUSS1, GAUSS2, GAUSS3, GAUSS4, LOBATTO3, RADAU15, GBS16, &
                                     TDRK24, TDRK24A, CASHKARP45
 use orbit_symplectic_quasi, only: f, si, orbit_timestep_rk45, orbit_timestep_radau15, &
                  orbit_timestep_gbs16, orbit_timestep_tdrk24, orbit_timestep_tdrk24a, &
                                      orbit_timestep_cashkarp45
    use orbit_systems_generated, only: euler1_residual, euler2_residual, &
                 midpoint_residual, gauss1_residual, gauss2_residual, gauss3_residual, &
                                       gauss4_residual
    use fortnum_multiroot, only: multiroot_hybrids, multiroot_hybrids_warm_t
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use vector_potentail_mod, only: torflux
    use diag_counters, only: count_event, EVT_R_NEGATIVE, EVT_WARNING_MAXIT_ACCEPT

    implicit none
    integer, parameter :: dp = kind(1.0d0)
    save

    logical, parameter :: exact_steps = .False.
! Largest stage count coeff_rk_gauss has coefficients for.
    integer, parameter :: S_MAX_GAUSS = 4

! Stage fields of one implicit Runge-Kutta solve. They reach the residual
! through ctx and live on the solving thread's stack, so concurrent particles
! cannot overwrite each other's stages.
    type :: stage_fields_t
        type(field_can_t), pointer :: fs(:) => null()
    end type stage_fields_t

    procedure(orbit_timestep_quasi_i), pointer :: orbit_timestep_quasi => null()

! Jacobian estimate of the previous step's implicit system, so that
! multiroot_hybrids does not rebuild it by n finite differences every step.
! One per composition stage (kwarm), because each stage solves with its own dt.
! Like the carry above it belongs to one orbit; orbit_sympl_init resets it.
! Used where it measured faster: euler2 (all composition schemes) and midpoint.
! For euler1 and Gauss the secant iterations it adds cost more than the n
! difference evaluations it saves.
    type(multiroot_hybrids_warm_t) :: warm(2*S_MAX)
    integer :: kwarm = 1
    !$omp threadprivate(warm, kwarm)

contains

    ! Select the stepper of mode and forget the previous orbit's warm start.
    subroutine select_quasi_stepper(mode)
        integer, intent(in) :: mode

        call reset_quasi_warm_start
        select case (mode)
        case (RK45)
            orbit_timestep_quasi => orbit_timestep_rk45
        case (EXPL_IMPL_EULER)
            orbit_timestep_quasi => timestep_expl_impl_euler_quasi
        case (IMPL_EXPL_EULER)
            orbit_timestep_quasi => timestep_impl_expl_euler_quasi
        case (MIDPOINT)
            orbit_timestep_quasi => timestep_midpoint_quasi
        case (GAUSS1)
            orbit_timestep_quasi => orbit_timestep_quasi_gauss1
        case (GAUSS2)
            orbit_timestep_quasi => orbit_timestep_quasi_gauss2
        case (GAUSS3)
            orbit_timestep_quasi => orbit_timestep_quasi_gauss3
        case (GAUSS4)
            orbit_timestep_quasi => orbit_timestep_quasi_gauss4
        case (LOBATTO3)
            orbit_timestep_quasi => orbit_timestep_quasi_lobatto3
        case (RADAU15)
            orbit_timestep_quasi => orbit_timestep_radau15
        case (GBS16)
            orbit_timestep_quasi => orbit_timestep_gbs16
        case (TDRK24)
            orbit_timestep_quasi => orbit_timestep_tdrk24
        case (TDRK24A)
            orbit_timestep_quasi => orbit_timestep_tdrk24a
        case (CASHKARP45)
            orbit_timestep_quasi => orbit_timestep_cashkarp45
        case default
            print *, 'invalid mode for select_quasi_stepper: ', mode
            error stop
        end select
    end subroutine select_quasi_stepper

    subroutine reset_quasi_warm_start
        integer :: k
        do k = 1, size(warm)
            warm(k)%valid = .false.
        end do
        kwarm = 1
    end subroutine reset_quasi_warm_start

! The quasi solves stop like the Newton steppers in orbit_symplectic: on
! |F| < atol, or on a full quasi-Newton step below rtol*tolref with the same
! tolref per unknown. |F| stalls at its rounding floor far above atol, so the
! step test is what ends a converged solve. Where that floor keeps the step
! above rtol*tolref, warning mode accepts, as for the Newton steppers, a root
! whose step falls below warning_factor times that. The relaxed solve restarts
! from the last iterate with a fresh Jacobian, since a secant estimate updated
! on noise can point far off. Otherwise the solve failed, and its status is the
! one the Newton steppers report.
    integer function quasi_solve_status(status)
        type(fortnum_status_t), intent(in) :: status

        quasi_solve_status = SYMPLECTIC_STEP_OK
        if (status%code == FORTNUM_DOMAIN_ERROR) then
            quasi_solve_status = SYMPLECTIC_STEP_LINEAR_SOLVE
        else if (status%code /= FORTNUM_OK) then
            quasi_solve_status = SYMPLECTIC_STEP_MAXITER
        end if
    end function quasi_solve_status

    !
    ! Wrapper routines for ODEPACK
    !

    subroutine f_exact_quasi(x, fvec, ctx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fvec(:)
        class(*), intent(in), optional :: ctx

        type(field_can_t) :: f2

        f2 = f

        call eval_field(f2, x(1), si%z(2), si%z(3), 0)
        call get_derivatives(f2, si%z(4))

        fvec(1) = f%pth - f2%pth

        f%H = f2%H

    end subroutine f_exact_quasi

    subroutine f_euler1_quasi(x, fvec, ctx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fvec(:)
        class(*), intent(in), optional :: ctx

        call eval_field(f, x(1), si%z(2), si%z(3), 0)
        call get_derivatives(f, x(2))

        call euler1_residual(si%dt, si%pthold, si%z, x, f, fvec)

    end subroutine f_euler1_quasi

    subroutine f_euler2_quasi(x, fvec, ctx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fvec(:)
        class(*), intent(in), optional :: ctx

        call eval_field(f, x(1), x(2), x(3), 0)
        call get_derivatives(f, si%z(4))

        call euler2_residual(si%dt, si%pthold, si%z, x, f, fvec)

    end subroutine f_euler2_quasi

    subroutine f_midpoint_quasi(x, fvec, ctx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fvec(:)
        class(*), intent(in), optional :: ctx

        type(field_can_t) :: fmid

        fmid = f
        call eval_field(fmid, x(5), 0.5*(x(2) + si%z(2)), 0.5*(x(3) + si%z(3)), 0)
        call get_derivatives(fmid, 0.5*(x(4) + si%z(4)))
        call eval_field(f, x(1), x(2), x(3), 0)
        call get_derivatives(f, x(4))
        call midpoint_residual(si%dt, si%pthold, si%z, x, fmid, f, fvec)

    end subroutine f_midpoint_quasi

    subroutine f_rk_gauss_quasi(x, fvec, ctx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fvec(:)
        class(*), intent(in), optional :: ctx

        select type (ctx)
        type is (stage_fields_t)
            call rk_gauss_quasi_residual(ctx%fs, x, fvec)
        class default
            error stop 'f_rk_gauss_quasi: ctx must carry the stage fields'
        end select
    end subroutine f_rk_gauss_quasi

    subroutine rk_gauss_quasi_residual(fs, x, fvec)

        type(field_can_t), intent(inout) :: fs(:)
        real(dp), intent(in) :: x(:)  ! = (rend, thend, phend, pphend, rmid)
        real(dp), intent(out) :: fvec(:)

        integer :: n, k
        real(dp) :: a(size(x)/4, size(x)/4), b(size(x)/4), c(size(x)/4)

        n = size(x)
        call coeff_rk_gauss(n/4, a, b, c)  ! TODO: move this to preprocessing

        do k = 1, n/4
            call eval_field(fs(k), x(4*k - 3), x(4*k - 2), x(4*k - 1), 0)
            call get_derivatives(fs(k), x(4*k))
        end do

        select case (n/4)
        case (1)
            call gauss1_residual(si%dt, si%pthold, si%z, a, x, fs, fvec)
        case (2)
            call gauss2_residual(si%dt, si%pthold, si%z, a, x, fs, fvec)
        case (3)
            call gauss3_residual(si%dt, si%pthold, si%z, a, x, fs, fvec)
        case (4)
            call gauss4_residual(si%dt, si%pthold, si%z, a, x, fs, fvec)
        end select

    end subroutine rk_gauss_quasi_residual

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    ! Lobatto (IIIA)-(IIIB) Runge-Kutta method, s stages (n=4*s-2 variables)
    !
    subroutine f_rk_lobatto_quasi(x, fvec, ctx)
        !
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: fvec(:)
        class(*), intent(in), optional :: ctx

        select type (ctx)
        type is (stage_fields_t)
            call f_rk_lobatto(si, ctx%fs, (size(x) + 2)/4, x, fvec, 0)
        class default
            error stop 'f_rk_lobatto_quasi: ctx must carry the stage fields'
        end select

    end subroutine f_rk_lobatto_quasi

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine orbit_timestep_multi_quasi(mi, ierr)
        !
        type(multistage_integrator_t), intent(inout) :: mi

        integer, intent(out) :: ierr

        integer :: kstage

        si = mi%stages(1)
        kwarm = 1
        call timestep_stage_quasi(ierr)
        if (ierr /= 0) return

        do kstage = 2, 2*mi%s
            mi%stages(kstage)%z = si%z
            mi%stages(kstage)%pthold = f%pth
            si = mi%stages(kstage)
            kwarm = kstage
            call timestep_stage_quasi(ierr)
            if (ierr /= 0) then
                kwarm = 1
                return
            end if
        end do
        kwarm = 1
        mi%stages(1)%z = si%z
        mi%stages(1)%pthold = f%pth

    end subroutine orbit_timestep_multi_quasi

! Composition stages alternate the Euler step and its adjoint, so each one
! runs the method stored with it, not the global orbit_timestep_quasi.
    subroutine timestep_stage_quasi(ierr)
        integer, intent(out) :: ierr

        select case (si%mode)
        case (EXPL_IMPL_EULER)
            call timestep_expl_impl_euler_quasi(ierr)
        case (IMPL_EXPL_EULER)
            call timestep_impl_expl_euler_quasi(ierr)
        case default
            error stop 'composition stage needs an Euler method'
        end select
    end subroutine timestep_stage_quasi

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine timestep_midpoint_quasi(ierr)
        !
        integer, intent(out) :: ierr

        integer, parameter :: n = 5

        real(dp), dimension(n) :: x, x0, tolref
        type(fortnum_status_t) :: status
        integer :: ktau
        type(symplectic_integrator_t) :: accepted_integrator
        type(field_can_t) :: accepted_field

        ierr = 0
        ktau = 0
        do while (ktau .lt. si%ntau)
            accepted_integrator = si
            accepted_field = f
            si%pthold = f%pth

            x0(1:4) = si%z
            x0(5) = si%z(1)

            tolref(1) = 1d0
            tolref(2) = twopi
            tolref(3) = twopi
            tolref(4) = max(dabs(f%Aph), dabs(1d1*torflux/f%ro0), dabs(si%z(4)))
            tolref(5) = 1d0

            call multiroot_hybrids(f_midpoint_quasi, n, x0, x, status, &
                  xtol=si%rtol, ftol=si%atol, warm=warm(kwarm), step_tol=si%rtol*tolref)
            if (status%code /= FORTNUM_OK .and. symplectic_newton_warning_mode) then
                x0 = x
                call multiroot_hybrids(f_midpoint_quasi, n, x0, x, status, &
                  xtol=si%rtol, ftol=si%atol, step_tol=si%warning_factor*si%rtol*tolref)
               if (status%code == FORTNUM_OK) call count_event(EVT_WARNING_MAXIT_ACCEPT)
            end if
            ierr = quasi_solve_status(status)
            if (ierr /= SYMPLECTIC_STEP_OK) then
                si = accepted_integrator
                f = accepted_field
                return
            end if

            if (x(1) > sympl_rmax) then
                warm(kwarm)%valid = .false.
                ierr = 1
                return
            end if

            if (x(1) < 0.0) then
                ! Axis crossing: continue on the opposite ray, (r, theta) ->
                ! (|r|, theta + pi), instead of the former teleport to r = 0.01 (#370).
                call count_event(EVT_R_NEGATIVE)
                warm(kwarm)%valid = .false.
                x(1) = -x(1)
                x(2) = x(2) + pi
                if (x(1) > sympl_rmax) then
                    ierr = 1
                    return
                end if
            end if

            si%z = x(1:4)

            ktau = ktau + 1
        end do

    end subroutine timestep_midpoint_quasi

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine timestep_expl_impl_euler_quasi(ierr)
        !
        integer, intent(out) :: ierr

        integer, parameter :: n = 2

        real(dp), dimension(n) :: x, x0, tolref
        real(dp), dimension(1) :: xe0, xe
        type(fortnum_status_t) :: status
        integer :: ktau
        type(symplectic_integrator_t) :: accepted_integrator
        type(field_can_t) :: accepted_field

        tolref(1) = 1d0
        tolref(2) = dabs(1d1*torflux/f%ro0)

        ierr = 0
        ktau = 0
        do while (ktau .lt. si%ntau)
            accepted_integrator = si
            accepted_field = f
            si%pthold = f%pth

            x0(1) = si%z(1)
            x0(2) = si%z(4)

            call multiroot_hybrids(f_euler1_quasi, n, x0, x, status, &
                                   xtol=si%rtol, ftol=si%atol, step_tol=si%rtol*tolref)
            if (status%code /= FORTNUM_OK .and. symplectic_newton_warning_mode) then
                x0 = x
                call multiroot_hybrids(f_euler1_quasi, n, x0, x, status, &
                  xtol=si%rtol, ftol=si%atol, step_tol=si%warning_factor*si%rtol*tolref)
               if (status%code == FORTNUM_OK) call count_event(EVT_WARNING_MAXIT_ACCEPT)
            end if
            ierr = quasi_solve_status(status)
            if (ierr /= SYMPLECTIC_STEP_OK) then
                si = accepted_integrator
                f = accepted_field
                return
            end if

            if (x(1) > sympl_rmax) then
                ierr = 1
                return
            end if

            if (x(1) < 0.0) then
                ! Axis crossing: continue on the opposite ray, (r, theta) ->
                ! (|r|, theta + pi), instead of the former teleport to r = 0.01 (#370).
                call count_event(EVT_R_NEGATIVE)
                x(1) = -x(1)
                x(2) = x(2) + pi
                if (x(1) > sympl_rmax) then
                    ierr = 1
                    return
                end if
            end if

            si%z(1) = x(1)
            si%z(4) = x(2)

            call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
            call get_derivatives(f, si%z(4))

            si%z(2) = si%z(2) + si%dt*f%dH(1)/f%dpth(1)
            si%z(3) = si%z(3) + si%dt*(f%vpar - f%dH(1)/f%dpth(1)*f%hth)/f%hph

            if (exact_steps) then
                xe0(1) = si%z(1)
                call multiroot_hybrids(f_exact_quasi, 1, xe0, xe, status, &
                                       xtol=si%rtol, ftol=si%rtol)
                si%z(1) = xe(1)
                call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
                call get_derivatives(f, si%z(4))
            end if

            ktau = ktau + 1
        end do

    end subroutine timestep_expl_impl_euler_quasi

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    subroutine timestep_impl_expl_euler_quasi(ierr)
        !
        integer, intent(out) :: ierr

        integer, parameter :: n = 3

        real(dp), dimension(n) :: x, x0
        real(dp), dimension(1) :: xe0, xe
        type(fortnum_status_t) :: status
        integer :: ktau
        type(symplectic_integrator_t) :: accepted_integrator
        type(field_can_t) :: accepted_field
        real(dp), parameter :: tolref(n) = [1d0, twopi, twopi]

        ierr = 0
        ktau = 0
        do while (ktau .lt. si%ntau)
            accepted_integrator = si
            accepted_field = f
            si%pthold = f%pth

            x0 = si%z(1:3)

            call multiroot_hybrids(f_euler2_quasi, n, x0, x, status, &
                  xtol=si%rtol, ftol=si%atol, warm=warm(kwarm), step_tol=si%rtol*tolref)
            if (status%code /= FORTNUM_OK .and. symplectic_newton_warning_mode) then
                x0 = x
                call multiroot_hybrids(f_euler2_quasi, n, x0, x, status, &
                  xtol=si%rtol, ftol=si%atol, step_tol=si%warning_factor*si%rtol*tolref)
               if (status%code == FORTNUM_OK) call count_event(EVT_WARNING_MAXIT_ACCEPT)
            end if
            ierr = quasi_solve_status(status)
            if (ierr /= SYMPLECTIC_STEP_OK) then
                si = accepted_integrator
                f = accepted_field
                return
            end if

            if (x(1) > sympl_rmax) then
                warm(kwarm)%valid = .false.
                ierr = 1
                return
            end if

            if (x(1) < 0.0) then
                ! Axis crossing: continue on the opposite ray, (r, theta) ->
                ! (|r|, theta + pi), instead of the former teleport to r = 0.01 (#370).
                call count_event(EVT_R_NEGATIVE)
                warm(kwarm)%valid = .false.
                x(1) = -x(1)
                x(2) = x(2) + pi
                if (x(1) > sympl_rmax) then
                    ierr = 1
                    return
                end if
            end if

            si%z(1:3) = x

            call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
            call get_derivatives(f, si%z(4))

            f%pth = si%pthold - si%dt*(f%dH(2) - f%dH(1)*f%dpth(2)/f%dpth(1))
            si%z(4) = si%z(4) - si%dt*(f%dH(3) - f%dH(1)*f%dpth(3)/f%dpth(1))

            if (exact_steps) then
                xe0(1) = si%z(1)
                call multiroot_hybrids(f_exact_quasi, 1, xe0, xe, status, &
                                       xtol=si%rtol, ftol=si%rtol)
                si%z(1) = xe(1)
                call eval_field(f, si%z(1), si%z(2), si%z(3), 0)
                call get_derivatives(f, si%z(4))
            end if

            ktau = ktau + 1
        end do

    end subroutine timestep_impl_expl_euler_quasi

    subroutine timestep_rk_gauss_quasi(s, ierr)
        !
        integer, intent(out) :: ierr

        integer, intent(in) :: s
        real(dp), dimension(4*s) :: x, x0
        type(fortnum_status_t) :: status
        integer :: ktau, k, l

        real(dp) :: a(s, s), b(s), c(s), Hprime(s), tolref(4*s)
        type(field_can_t), target :: fs(s)
        type(stage_fields_t) :: stages
        type(symplectic_integrator_t) :: accepted_integrator
        type(field_can_t) :: accepted_field

        if (s < 1 .or. s > S_MAX_GAUSS) &
            error stop 'timestep_rk_gauss_quasi: no Gauss coefficients for this s'

        do k = 1, s
            fs(k) = f
        end do
        stages%fs => fs

        ierr = 0
        ktau = 0
        do while (ktau .lt. si%ntau)
            accepted_integrator = si
            accepted_field = f
            si%pthold = f%pth

            do k = 1, s
                x0((4*k - 3):(4*k)) = si%z
                tolref((4*k - 3):(4*k)) = [1d0, twopi, twopi, dabs(si%z(4))]
            end do

            call multiroot_hybrids(f_rk_gauss_quasi, 4*s, x0, x, status, &
                        xtol=si%rtol, ftol=si%atol, ctx=stages, step_tol=si%rtol*tolref)
            if (status%code /= FORTNUM_OK .and. symplectic_newton_warning_mode) then
                x0 = x
                call multiroot_hybrids(f_rk_gauss_quasi, 4*s, x0, x, status, &
      xtol=si%rtol, ftol=si%atol, ctx=stages, step_tol=si%warning_factor*si%rtol*tolref)
               if (status%code == FORTNUM_OK) call count_event(EVT_WARNING_MAXIT_ACCEPT)
            end if
            ierr = quasi_solve_status(status)
            if (ierr /= SYMPLECTIC_STEP_OK) then
                si = accepted_integrator
                f = accepted_field
                return
            end if

            if (x(1) > sympl_rmax) then
                ierr = 1
                return
            end if

            if (x(1) < 0.0) then
                ! Axis crossing: continue on the opposite ray, (r, theta) ->
                ! (|r|, theta + pi), instead of the former teleport to r = 0.01 (#370).
                call count_event(EVT_R_NEGATIVE)
                x(1) = -x(1)
                x(2) = x(2) + pi
                if (x(1) > sympl_rmax) then
                    ierr = 1
                    return
                end if
            end if

            call coeff_rk_gauss(s, a, b, c)

            f = fs(s)
            f%pth = si%pthold
            si%z(1) = x(4*s - 3)

            do l = 1, s
                Hprime(l) = fs(l)%dH(1)/fs(l)%dpth(1)
                f%pth = f%pth - si%dt*b(l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))
                si%z(2) = si%z(2) + si%dt*b(l)*Hprime(l)
             si%z(3) = si%z(3) + si%dt*b(l)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph
                si%z(4) = si%z(4) - si%dt*b(l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))
            end do

            ktau = ktau + 1
        end do

    end subroutine timestep_rk_gauss_quasi

    subroutine timestep_rk_lobatto_quasi(s, ierr)
        !
        integer, intent(out) :: ierr

        integer, intent(in) :: s
        real(dp), dimension(4*s - 2) :: x, x0
        type(fortnum_status_t) :: status
        integer :: ktau, k, l

        real(dp) :: a(s, s), ahat(s, s), b(s), c(s), Hprime(s), tolref(4*s - 2)
        type(field_can_t), target :: fs(s)
        type(stage_fields_t) :: stages
        type(symplectic_integrator_t) :: accepted_integrator
        type(field_can_t) :: accepted_field

        do k = 1, s
            fs(k) = f
        end do
        stages%fs => fs

        ierr = 0
        ktau = 0
        do while (ktau .lt. si%ntau)
            accepted_integrator = si
            accepted_field = f
            si%pthold = f%pth

            x0(1) = si%z(1)
            x0(2) = si%z(4)
            tolref(1:2) = [1d0, dabs(si%z(4))]
            do k = 2, s
                x0((4*k - 3 - 2):(4*k - 2)) = si%z
                tolref((4*k - 3 - 2):(4*k - 2)) = [1d0, twopi, twopi, dabs(si%z(4))]
            end do

            !
            ! q1 = q0
            ! q2 = q0 + h*(a21*f(w1) + a22*f(z2) + a23*f(z3)
            ! q3 = q0 + h*(a31*f(w1) + a32*f(z2) + a33*f(z3))
            !
            ! p1(w1) = p0 + h*(o11*f(z1) + o12*f(z2))
            ! p2(z2) = p0 + h*(o21*f(z1) + o22*f(z2))
            ! p3(z3) = p0 + h*(o31*f(z1) + o32*f(z2))
            !
            ! 0 0 0 0 0 0
            ! x x x x x x
            ! x x x x x x
            ! x x x x 0 0
            ! x x x x 0 0
            ! x x x x 0 0
            !
            ! effectively
            !
            ! x x x x x
            ! x x x x x
            ! x x x 0 0
            ! x x x 0 0
            ! x x x 0 0

            call multiroot_hybrids(f_rk_lobatto_quasi, 4*s - 2, x0, x, status, &
                        xtol=si%rtol, ftol=si%atol, ctx=stages, step_tol=si%rtol*tolref)
            if (status%code /= FORTNUM_OK .and. symplectic_newton_warning_mode) then
                x0 = x
                call multiroot_hybrids(f_rk_lobatto_quasi, 4*s - 2, x0, x, status, &
      xtol=si%rtol, ftol=si%atol, ctx=stages, step_tol=si%warning_factor*si%rtol*tolref)
               if (status%code == FORTNUM_OK) call count_event(EVT_WARNING_MAXIT_ACCEPT)
            end if
            ierr = quasi_solve_status(status)
            if (ierr /= SYMPLECTIC_STEP_OK) then
                si = accepted_integrator
                f = accepted_field
                return
            end if

            call coeff_rk_lobatto(s, a, ahat, b, c)

            call eval_field(fs(1), x(1), si%z(2), si%z(3), 0)
            call get_derivatives(fs(1), x(2))
            Hprime(1) = fs(1)%dH(1)/fs(1)%dpth(1)

            do k = 2, s
               call eval_field(fs(k), x(4*k - 3 - 2), x(4*k - 2 - 2), x(4*k - 1 - 2), 0)
                call get_derivatives(fs(k), x(4*k - 2))
                Hprime(k) = fs(k)%dH(1)/fs(k)%dpth(1)
            end do

            f = fs(s)
            f%pth = si%pthold
            si%z(1) = x(4*s - 3 - 2)

            do l = 1, s
                f%pth = f%pth - si%dt*b(l)*(fs(l)%dH(2) - Hprime(l)*fs(l)%dpth(2))
                si%z(2) = si%z(2) + si%dt*b(l)*Hprime(l)
             si%z(3) = si%z(3) + si%dt*b(l)*(fs(l)%vpar - Hprime(l)*fs(l)%hth)/fs(l)%hph
                si%z(4) = si%z(4) - si%dt*b(l)*(fs(l)%dH(3) - Hprime(l)*fs(l)%dpth(3))
            end do

            ktau = ktau + 1
        end do

    end subroutine timestep_rk_lobatto_quasi

    recursive subroutine orbit_timestep_quasi_gauss1(ierr)
        integer, intent(out) :: ierr
        call timestep_rk_gauss_quasi(1, ierr)
    end subroutine orbit_timestep_quasi_gauss1

    recursive subroutine orbit_timestep_quasi_gauss2(ierr)
        integer, intent(out) :: ierr
        call timestep_rk_gauss_quasi(2, ierr)
    end subroutine orbit_timestep_quasi_gauss2

    recursive subroutine orbit_timestep_quasi_gauss3(ierr)
        integer, intent(out) :: ierr
        call timestep_rk_gauss_quasi(3, ierr)
    end subroutine orbit_timestep_quasi_gauss3

    recursive subroutine orbit_timestep_quasi_gauss4(ierr)
        integer, intent(out) :: ierr
        call timestep_rk_gauss_quasi(4, ierr)
    end subroutine orbit_timestep_quasi_gauss4

    recursive subroutine orbit_timestep_quasi_lobatto3(ierr)
        integer, intent(out) :: ierr
        call timestep_rk_lobatto_quasi(3, ierr)
    end subroutine orbit_timestep_quasi_lobatto3

end module orbit_quasi_implicit
