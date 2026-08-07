module simple_gpu
    ! OpenACC GPU offload of the default orbit-tracing hot path:
    ! Boozer field (isw_field_type=2), one particle per GPU thread. The
    ! fixed-step symplectic Euler path and allocation-free adaptive RK5(4)
    ! paths share the same resident field representation.
    !
    ! The CPU integrator dispatches the field evaluation through a procedure
    ! pointer, which cannot be called inside an OpenACC compute region. This
    ! module keeps the device-side field evaluation local and reuses the shared
    ! symplectic-Euler algebra from orbit_symplectic_euler1. Equivalence is
    ! checked against the CPU path by test_gpu_orbit_bench.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use util, only: pi, twopi
    use field_can_mod, only: field_can_t, get_val, get_derivatives, get_derivatives2
    use field_can_boozer, only: eval_field_booz_device
    use orbit_symplectic_base, only: symplectic_integrator_t, &
        SYMPLECTIC_STEP_OK, SYMPLECTIC_STEP_OUTSIDE_DOMAIN, &
        SYMPLECTIC_STEP_MAXITER, SYMPLECTIC_STEP_LINEAR_SOLVE, &
        SYMPLECTIC_STEP_BOUNDARY_LIMITED, symplectic_newton_warning_mode
    use orbit_symplectic_euler1, only: sympl_euler1_newton_iter
    use orbit_symplectic_euler1, only: sympl_euler1_extrapolate_field, sympl_euler1_advance_angles
    use orbit_symplectic_base, only: EXPL_IMPL_EULER, MIDPOINT, CASH_KARP, &
        DORMAND_PRINCE
    use orbit_rk_core, only: rk_solve
    use fortnum_ode_rk54_device, only: rk54_controls4_t, rk54_state4_t, &
        rk54_initialize4, rk54_request4, rk54_supply4, &
        RK54_CASH_KARP, RK54_DORMAND_PRINCE, RK54_NEED_RHS, RK54_ACCEPTED, &
        RK54_REJECTED, RK54_FAILED
    use boozer_sub, only: boozer_state
    use omp_lib, only: omp_get_thread_num
#ifdef _OPENACC
    use openacc, only: acc_get_num_devices, acc_set_device_num, acc_device_nvidia
#endif

    implicit none
    private

    integer, parameter :: maxit = 32

    public :: trace_orbits_gpu, trace_orbits_gpu_range
    public :: trace_orbits_gpu_method, trace_orbits_gpu_rk54_range
    public :: trace_orbits_gpu_midpoint_range

contains

    logical function gpu_finite_iterate(x)
        !$acc routine seq
        real(dp), intent(in) :: x(:)

        integer(int64), parameter :: exponent_mask = &
            int(z'7FF0000000000000', int64)
        integer(int64) :: bits
        integer :: i

        gpu_finite_iterate = .false.
        do i = 1, size(x)
            bits = transfer(x(i), bits)
            if (iand(bits, exponent_mask) == exponent_mask) return
        end do
        gpu_finite_iterate = .true.
    end function gpu_finite_iterate

    subroutine gpu_newton1(si, f, x, xlast, warning_mode, status, nfev)
        !$acc routine seq
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(inout) :: x(2)
        real(dp), intent(out) :: xlast(2)
        logical, intent(in) :: warning_mode
        integer, intent(out) :: status
        integer, intent(out) :: nfev

        real(dp) :: tolref(2)
        integer :: kit
        logical :: converged, linear_failed, boundary_limited
        logical :: step_boundary_limited

        status = SYMPLECTIC_STEP_MAXITER
        nfev = 0
        boundary_limited = .false.

        tolref(1) = 1d0
        tolref(2) = dabs(1d1*boozer_state%torflux/f%ro0)

        do kit = 1, maxit
            if (x(1) > 1d0) then
                status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
                return
            end if
            ! Transient guard; converged-negative handled by the caller
            ! (mirrors newton1 in orbit_symplectic, #370).
            if (x(1) < 0d0) x(1) = 0.01d0

            call eval_field_booz_device(f, x(1), si%z(2), si%z(3), 2)
            nfev = nfev + 1
            call get_derivatives2(f, x(2))
            call sympl_euler1_newton_iter(si, f, x, tolref, xlast, converged, &
                linear_failed, step_boundary_limited)
            boundary_limited = boundary_limited .or. step_boundary_limited

            if (linear_failed) then
                status = SYMPLECTIC_STEP_LINEAR_SOLVE
                return
            end if
            if (x(1) > 1d0) then
                status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
                return
            end if
            if (converged) then
                status = SYMPLECTIC_STEP_OK
                return
            end if
        end do
        if (boundary_limited) then
            if (step_boundary_limited .and. &
                abs(1d0 - x(1)) <= max(1d-10, 10d0*si%rtol)) then
                status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
            else
                status = SYMPLECTIC_STEP_BOUNDARY_LIMITED
            end if
        else if (warning_mode .and. gpu_finite_iterate(x) .and. &
                 gpu_finite_iterate(xlast) .and. gpu_finite_iterate(tolref) .and. &
                 all(abs(x - xlast) <= &
                     si%warning_factor*si%rtol*abs(tolref))) then
            status = SYMPLECTIC_STEP_OK
        end if
        ! Non-convergence diagnostics (CPU writes fort.6601) are omitted on device.
    end subroutine gpu_newton1

    subroutine gpu_timestep_euler(si, f, warning_mode, ierr, nfev)
        !$acc routine seq
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        logical, intent(in) :: warning_mode
        integer, intent(out) :: ierr
        integer, intent(out) :: nfev

        real(dp) :: x(2), xlast(2)
        integer :: ktau, newton_status, newton_nfev
        logical :: crossed
        type(symplectic_integrator_t) :: accepted_integrator
        type(field_can_t) :: accepted_field

        ierr = 0
        nfev = 0
        ktau = 0
        do while (ktau < si%ntau)
            accepted_integrator = si
            accepted_field = f
            si%pthold = f%pth

            x(1) = si%z(1)
            x(2) = si%z(4)

            call gpu_newton1(si, f, x, xlast, warning_mode, newton_status, &
                newton_nfev)
            nfev = nfev + newton_nfev
            if (newton_status /= SYMPLECTIC_STEP_OK) then
                si = accepted_integrator
                f = accepted_field
                ierr = newton_status
                return
            end if

            if (x(1) > 1.0d0) then
                si = accepted_integrator
                f = accepted_field
                ierr = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
                return
            end if
            crossed = .false.
            if (x(1) < 0.0d0) then
                ! The converged radius lies beyond the axis: commit the
                ! chart switch (r, theta) -> (|r|, theta + pi) (#370).
                x(1) = -x(1)
                si%z(2) = si%z(2) + pi
                crossed = .true.
                if (x(1) > 1.0d0) then
                    si = accepted_integrator
                    f = accepted_field
                    ierr = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
                    return
                end if
            end if

            si%z(1) = x(1)
            si%z(4) = x(2)

            if (crossed) then
                ! xlast lives in the other chart; extrapolation across the
                ! flip is invalid, evaluate the field fresh instead.
                call eval_field_booz_device(f, si%z(1), si%z(2), si%z(3), 0)
                call get_derivatives(f, si%z(4))
                nfev = nfev + 1
            else
                call sympl_euler1_extrapolate_field(si, f, x, xlast)
            end if
            call sympl_euler1_advance_angles(si, f)

            ktau = ktau + 1
        end do
    end subroutine gpu_timestep_euler

    subroutine gpu_midpoint_system(si, f, x, residual, jacobian)
        ! Five-equation degenerate-variational midpoint system and analytic
        ! Jacobian, specialized to the resident Boozer field evaluator.
        !$acc routine seq
        type(symplectic_integrator_t), intent(in) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: x(5)
        real(dp), intent(out) :: residual(5), jacobian(5, 5)

        type(field_can_t) :: fmid
        real(dp) :: dpthmid, pthdotbar

        call eval_field_booz_device(f, x(5), 0.5_dp*(x(2) + si%z(2)), &
            0.5_dp*(x(3) + si%z(3)), 2)
        call get_derivatives2(f, 0.5_dp*(x(4) + si%z(4)))

        residual(2) = f%dpth(1)*(x(2) - si%z(2)) - si%dt*f%dH(1)
        residual(3) = f%dpth(1)*f%hph*(x(3) - si%z(3)) - &
            si%dt*(f%dpth(1)*f%vpar - f%dH(1)*f%hth)
        residual(4) = f%dpth(1)*(x(4) - si%z(4)) + &
            si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))
        residual(5) = f%dpth(1)*(f%pth - si%pthold) + 0.5_dp*si%dt*( &
            f%dpth(1)*f%dH(2) - f%dpth(2)*f%dH(1))

        jacobian(2, 1) = 0.0_dp
        jacobian(2, 5) = f%d2pth(1)*(x(2) - si%z(2)) - si%dt*f%d2H(1)
        jacobian(2, 2:3) = 0.5_dp*(f%d2pth(2:3)*(x(2) - si%z(2)) - &
            si%dt*f%d2H(2:3))
        jacobian(2, 2) = jacobian(2, 2) + f%dpth(1)
        jacobian(2, 4) = 0.5_dp*(f%d2pth(7)*(x(2) - si%z(2)) - &
            si%dt*f%d2H(7))

        jacobian(3, 1) = 0.0_dp
        jacobian(3, 5) = (f%d2pth(1)*f%hph + f%dpth(1)*f%dhph(1))* &
            (x(3) - si%z(3)) - si%dt*(f%d2pth(1)*f%vpar + &
            f%dpth(1)*f%dvpar(1) - f%d2H(1)*f%hth - f%dH(1)*f%dhth(1))
        jacobian(3, 2:3) = 0.5_dp*((f%d2pth(2:3)*f%hph + &
            f%dpth(1)*f%dhph(2:3))*(x(3) - si%z(3)) - si%dt*( &
            f%d2pth(2:3)*f%vpar + f%dpth(1)*f%dvpar(2:3) - &
            f%d2H(2:3)*f%hth - f%dH(1)*f%dhth(2:3)))
        jacobian(3, 3) = jacobian(3, 3) + f%dpth(1)*f%hph
        jacobian(3, 4) = 0.5_dp*(f%d2pth(7)*f%hph*(x(3) - si%z(3)) - &
            si%dt*(f%d2pth(7)*f%vpar + f%dpth(1)*f%dvpar(4) - &
            f%d2H(7)*f%hth))

        jacobian(4, 1) = 0.0_dp
        jacobian(4, 5) = f%d2pth(1)*(x(4) - si%z(4)) + si%dt*( &
            f%d2H(3)*f%dpth(1) + f%dH(3)*f%d2pth(1) - &
            f%d2H(1)*f%dpth(3) - f%dH(1)*f%d2pth(3))
        jacobian(4, 2) = 0.5_dp*(f%d2pth(2)*(x(4) - si%z(4)) + si%dt*( &
            f%d2H(5)*f%dpth(1) + f%dH(3)*f%d2pth(2) - &
            f%d2H(2)*f%dpth(3) - f%dH(1)*f%d2pth(5)))
        jacobian(4, 3) = 0.5_dp*(f%d2pth(3)*(x(4) - si%z(4)) + si%dt*( &
            f%d2H(6)*f%dpth(1) + f%dH(3)*f%d2pth(3) - &
            f%d2H(3)*f%dpth(3) - f%dH(1)*f%d2pth(6)))
        jacobian(4, 4) = 0.5_dp*(f%d2pth(7)*(x(4) - si%z(4)) + si%dt*( &
            f%dH(3)*f%d2pth(7) + f%d2H(9)*f%dpth(1) - &
            f%d2H(7)*f%dpth(3) - f%dH(1)*f%d2pth(9))) + f%dpth(1)

        jacobian(5, 1) = 0.0_dp
        jacobian(5, 5) = f%d2pth(1)*(f%pth - si%pthold) + &
            f%dpth(1)*f%dpth(1) + 0.5_dp*si%dt*(f%d2pth(1)*f%dH(2) + &
            f%dpth(1)*f%d2H(2) - f%d2pth(2)*f%dH(1) - f%dpth(2)*f%d2H(1))
        jacobian(5, 2) = 0.5_dp*(f%d2pth(2)*(f%pth - si%pthold) + &
            f%dpth(1)*f%dpth(2) + 0.5_dp*si%dt*(f%d2pth(2)*f%dH(2) + &
            f%dpth(1)*f%d2H(4) - f%d2pth(4)*f%dH(1) - f%dpth(2)*f%d2H(2)))
        jacobian(5, 3) = 0.5_dp*(f%d2pth(3)*(f%pth - si%pthold) + &
            f%dpth(1)*f%dpth(3) + 0.5_dp*si%dt*(f%d2pth(3)*f%dH(2) + &
            f%dpth(1)*f%d2H(5) - f%d2pth(5)*f%dH(1) - f%dpth(2)*f%d2H(3)))
        jacobian(5, 4) = 0.5_dp*(f%d2pth(7)*(f%pth - si%pthold) + &
            f%dpth(1)*f%dpth(4) + 0.5_dp*si%dt*(f%d2pth(7)*f%dH(2) + &
            f%dpth(1)*f%d2H(8) - f%d2pth(8)*f%dH(1) - f%dpth(2)*f%d2H(7)))

        fmid = f
        dpthmid = f%dpth(1)
        pthdotbar = f%dpth(1)*f%dH(2) - f%dpth(2)*f%dH(1)
        call eval_field_booz_device(f, x(1), x(2), x(3), 2)
        call get_derivatives2(f, x(4))
        residual(1) = dpthmid*(f%pth - si%pthold) + si%dt*pthdotbar

        jacobian(1, 1) = fmid%dpth(1)*f%dpth(1)
        jacobian(1, 2) = 0.5_dp*(fmid%d2pth(2)*(f%pth - si%pthold) + &
            si%dt*(fmid%d2pth(2)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(4) - &
            fmid%dpth(2)*fmid%d2H(2) - fmid%d2pth(4)*fmid%dH(1))) + &
            fmid%dpth(1)*f%dpth(2)
        jacobian(1, 3) = 0.5_dp*(fmid%d2pth(3)*(f%pth - si%pthold) + &
            si%dt*(fmid%d2pth(3)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(5) - &
            fmid%dpth(2)*fmid%d2H(3) - fmid%d2pth(5)*fmid%dH(1))) + &
            fmid%dpth(1)*f%dpth(3)
        jacobian(1, 4) = 0.5_dp*(fmid%d2pth(7)*(f%pth - si%pthold) + &
            si%dt*(fmid%d2pth(7)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(8) - &
            fmid%dpth(2)*fmid%d2H(7) - fmid%d2pth(8)*fmid%dH(1))) + &
            fmid%dpth(1)*f%dpth(4)
        jacobian(1, 5) = fmid%d2pth(1)*(f%pth - si%pthold) + si%dt*( &
            fmid%d2pth(1)*fmid%dH(2) + fmid%dpth(1)*fmid%d2H(2) - &
            fmid%dpth(2)*fmid%d2H(1) - fmid%d2pth(2)*fmid%dH(1))
    end subroutine gpu_midpoint_system

    subroutine gpu_newton_midpoint(si, f, x, warning_mode, status, nfev)
        !$acc routine seq
        type(symplectic_integrator_t), intent(in) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(inout) :: x(5)
        logical, intent(in) :: warning_mode
        integer, intent(out) :: status, nfev

        integer :: kit, info, ir
        integer, parameter :: radial_indices(2) = [1, 5]
        real(dp) :: residual(5), jacobian(5, 5), xlast(5), tolref(5)
        real(dp) :: residual_abs(5), correction_abs(5), scale, outward
        logical :: limited

        status = SYMPLECTIC_STEP_MAXITER
        nfev = 0
        tolref = [1.0_dp, twopi, twopi, &
            max(abs(f%Aph), abs(10.0_dp*boozer_state%torflux/f%ro0)), 1.0_dp]
        xlast = x

        do kit = 1, maxit
            if (x(1) > 1.0_dp .or. x(5) > 1.0_dp) then
                status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
                return
            end if
            if (x(1) < 0.0_dp) x(1) = 0.01_dp
            if (x(5) < 0.0_dp) x(5) = 0.01_dp

            call gpu_midpoint_system(si, f, x, residual, jacobian)
            nfev = nfev + 2
            residual_abs = abs(residual)
            xlast = x
            call rk_solve(5, jacobian, residual, info)
            if (info /= 0) then
                status = SYMPLECTIC_STEP_LINEAR_SOLVE
                return
            end if

            scale = 1.0_dp
            limited = .false.
            do ir = 1, 2
                outward = -residual(radial_indices(ir))
                if (outward > 0.0_dp) then
                    scale = min(scale, 0.8_dp*max(0.0_dp, &
                        1.0_dp - x(radial_indices(ir)))/outward)
                end if
            end do
            limited = scale < 1.0_dp
            x = x - scale*residual
            if (x(1) > 1.0_dp .or. x(5) > 1.0_dp) then
                status = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
                return
            end if
            correction_abs = abs(x - xlast)
            tolref(4) = max(abs(x(4)), tolref(4))
            if (.not. limited .and. (all(residual_abs < si%atol) .or. &
                    all(correction_abs < si%rtol*tolref))) then
                status = SYMPLECTIC_STEP_OK
                return
            end if
        end do

        if (warning_mode .and. gpu_finite_iterate(x) .and. &
                gpu_finite_iterate(xlast) .and. gpu_finite_iterate(tolref) .and. &
                all(abs(x - xlast) <= si%warning_factor*si%rtol*abs(tolref))) &
            status = SYMPLECTIC_STEP_OK
    end subroutine gpu_newton_midpoint

    subroutine gpu_timestep_midpoint(si, f, warning_mode, ierr, nfev)
        !$acc routine seq
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        logical, intent(in) :: warning_mode
        integer, intent(out) :: ierr, nfev

        type(symplectic_integrator_t) :: accepted_integrator
        type(field_can_t) :: accepted_field
        real(dp) :: x(5)
        integer :: ktau, step_status, step_nfev

        ierr = 0
        nfev = 0
        do ktau = 1, si%ntau
            accepted_integrator = si
            accepted_field = f
            si%pthold = f%pth
            x(1:4) = si%z
            x(5) = si%z(1)
            call gpu_newton_midpoint(si, f, x, warning_mode, step_status, step_nfev)
            nfev = nfev + step_nfev
            if (step_status /= SYMPLECTIC_STEP_OK) then
                si = accepted_integrator
                f = accepted_field
                ierr = step_status
                return
            end if
            if (x(1) < 0.0_dp) then
                x(1) = -x(1)
                x(2) = x(2) + pi
            end if
            if (x(1) > 1.0_dp) then
                si = accepted_integrator
                f = accepted_field
                ierr = SYMPLECTIC_STEP_OUTSIDE_DOMAIN
                return
            end if
            si%z = x(1:4)
            call eval_field_booz_device(f, si%z(1), si%z(2), si%z(3), 0)
            call get_val(f, si%z(4))
            nfev = nfev + 1
        end do
    end subroutine gpu_timestep_midpoint

    subroutine gpu_rhs_canonical(y, field_template, dydt)
        ! Canonical guiding-centre equations used by SIMPLE's CPU RK45 path,
        ! with the Boozer evaluator bound statically for device compilation.
        !$acc routine seq
        real(dp), intent(in) :: y(4)
        type(field_can_t), intent(in) :: field_template
        real(dp), intent(out) :: dydt(4)

        type(field_can_t) :: feval
        real(dp) :: hprime

        feval = field_template
        call eval_field_booz_device(feval, y(1), y(2), y(3), 0)
        call get_derivatives(feval, y(4))
        hprime = feval%dH(1)/feval%dpth(1)

        dydt(1) = -(feval%dH(2) - feval%hth/feval%hph*feval%dH(3))/ &
            feval%dpth(1)
        dydt(2) = hprime
        dydt(3) = (feval%vpar - hprime*feval%hth)/feval%hph
        dydt(4) = -(feval%dH(3) - hprime*feval%dpth(3))
    end subroutine gpu_rhs_canonical

    subroutine gpu_trace_rk54(field_template, y, duration, controls, loss_time, &
            nfev, ierr)
        !$acc routine seq
        type(field_can_t), intent(in) :: field_template
        real(dp), intent(inout) :: y(4)
        real(dp), intent(in) :: duration
        type(rk54_controls4_t), intent(in) :: controls
        real(dp), intent(out) :: loss_time
        integer, intent(out) :: nfev, ierr

        integer, parameter :: max_attempts = 100000000
        type(rk54_state4_t) :: state
        real(dp) :: t_eval, y_eval(4), dydt(4), h0
        integer :: request, attempt

        loss_time = duration
        nfev = 0
        ierr = 0
        if (duration <= 0.0_dp) return

        h0 = min(duration, max(controls%hmin, 1.0e-3_dp*controls%hmax))
        call rk54_initialize4(state, 0.0_dp, y, h0)

        do attempt = 1, max_attempts
            if (state%t >= duration) exit
            if (state%y(1) >= 1.0_dp) then
                loss_time = state%t
                exit
            end if

            state%h = min(state%h, duration - state%t)
            call rk54_request4(state, controls, t_eval, y_eval, request)
            do while (request == RK54_NEED_RHS)
                call gpu_rhs_canonical(y_eval, field_template, dydt)
                call rk54_supply4(state, controls, dydt, t_eval, y_eval, request)
            end do

            if (request == RK54_FAILED) then
                ierr = RK54_FAILED
                exit
            end if
            if (request == RK54_ACCEPTED) then
                if (state%y(1) < 0.0_dp) then
                    state%y(1) = -state%y(1)
                    state%y(2) = state%y(2) + pi
                end if
                if (state%y(1) >= 1.0_dp) then
                    loss_time = state%t
                    exit
                end if
            else if (request /= RK54_REJECTED) then
                ierr = request
                exit
            end if
        end do
        if (attempt > max_attempts) ierr = RK54_FAILED

        y = state%y
        nfev = state%nfev
    end subroutine gpu_trace_rk54

    subroutine trace_orbits_gpu_rk54_range(si, f, npart, istart, iend, &
            ntimstep, ntau_macro, method, loss_step, loss_time_out, zend, nfev)
        type(symplectic_integrator_t), intent(in) :: si(npart)
        type(field_can_t), intent(in) :: f(npart)
        integer, intent(in) :: npart, istart, iend, ntimstep, method
        integer, intent(in) :: ntau_macro(ntimstep)
        integer, intent(out) :: loss_step(npart), nfev(npart)
        real(dp), intent(out) :: loss_time_out(npart)
        real(dp), intent(out) :: zend(4, npart)

        type(rk54_controls4_t) :: controls
        real(dp) :: duration, loss_time, y(4), elapsed
        integer :: i, it, ierr

        duration = si(istart)%dt*real(sum(ntau_macro(2:ntimstep)), dp)
        controls%rtol = si(istart)%rtol
        controls%atol = si(istart)%rtol
        controls%hmin = 0.0_dp
        controls%hmax = duration
        if (method == CASH_KARP) then
            controls%method = RK54_CASH_KARP
        else
            controls%method = RK54_DORMAND_PRINCE
        end if

        !$acc parallel loop gang vector default(present) &
        !$acc&   copyin(si(istart:iend), f(istart:iend), ntau_macro) &
        !$acc&   copyout(loss_step(istart:iend), zend(:, istart:iend), &
        !$acc&           nfev(istart:iend), loss_time_out(istart:iend)) &
        !$acc&   private(y, loss_time, elapsed, ierr, it) firstprivate(controls, duration)
        do i = istart, iend
            y = si(i)%z
            call gpu_trace_rk54(f(i), y, duration, controls, loss_time, &
                nfev(i), ierr)
            loss_step(i) = ntimstep
            if (loss_time < duration) then
                elapsed = 0.0_dp
                do it = 2, ntimstep
                    elapsed = elapsed + si(i)%dt*real(ntau_macro(it), dp)
                    if (loss_time <= elapsed) then
                        loss_step(i) = it
                        exit
                    end if
                end do
            end if
            if (ierr /= 0) loss_step(i) = -ierr
            loss_time_out(i) = loss_time
            zend(:, i) = y
        end do
    end subroutine trace_orbits_gpu_rk54_range

    subroutine trace_orbits_gpu_midpoint_range(si, f, npart, istart, iend, &
            ntimstep, ntau_macro, loss_step, loss_time_out, zend, nfev)
        type(symplectic_integrator_t), intent(inout) :: si(npart)
        type(field_can_t), intent(inout) :: f(npart)
        integer, intent(in) :: npart, istart, iend, ntimstep
        integer, intent(in) :: ntau_macro(ntimstep)
        integer, intent(out) :: loss_step(npart), nfev(npart)
        real(dp), intent(out) :: loss_time_out(npart)
        real(dp), intent(out) :: zend(4, npart)

        integer :: i, it, ktau, ierr, lstep, step_nfev
        real(dp) :: elapsed
        logical :: warning_mode

        warning_mode = symplectic_newton_warning_mode
        !$acc parallel loop gang vector default(present) &
        !$acc&   copy(si(istart:iend), f(istart:iend)) copyin(ntau_macro) &
        !$acc&   copyout(loss_step(istart:iend), zend(:, istart:iend), &
        !$acc&           nfev(istart:iend), loss_time_out(istart:iend)) &
        !$acc&   private(it, ktau, ierr, lstep, step_nfev, elapsed) &
        !$acc&   firstprivate(warning_mode)
        do i = istart, iend
            ierr = 0
            lstep = ntimstep
            nfev(i) = 0
            elapsed = 0.0_dp
            macro: do it = 2, ntimstep
                do ktau = 1, ntau_macro(it)
                    call gpu_timestep_midpoint(si(i), f(i), warning_mode, ierr, &
                        step_nfev)
                    nfev(i) = nfev(i) + step_nfev
                    elapsed = elapsed + si(i)%dt
                    if (ierr /= 0) exit
                end do
                if (ierr /= 0) then
                    lstep = it
                    exit macro
                end if
            end do macro
            loss_step(i) = lstep
            if (ierr == 0 .and. lstep == ntimstep) &
                elapsed = real(sum(ntau_macro(2:ntimstep)), dp)*si(i)%dt
            loss_time_out(i) = elapsed
            zend(:, i) = si(i)%z
        end do
    end subroutine trace_orbits_gpu_midpoint_range

    subroutine trace_orbits_gpu_range(si, f, npart, istart, iend, ntimstep, &
                                      ntau_macro, loss_step, loss_time, zend, nfev)
        ! Evolve particles istart:iend on the currently selected OpenACC device.
        ! loss_step(i) = first timestep at which the orbit left the plasma (r>1),
        ! or ntimstep if confined for the whole trace. zend(:,i) is the final
        ! integrator state z = (r, th, ph, pphi).
        type(symplectic_integrator_t), intent(inout) :: si(npart)
        type(field_can_t), intent(inout) :: f(npart)
        integer, intent(in) :: npart, istart, iend, ntimstep
        integer, intent(in) :: ntau_macro(ntimstep)
        integer, intent(out) :: loss_step(npart)
        real(dp), intent(out) :: loss_time(npart)
        real(dp), intent(out) :: zend(4, npart)
        integer, intent(out) :: nfev(npart)

        integer :: i, it, ktau, ierr, lstep, step_nfev
        real(dp) :: elapsed
        logical :: warning_mode

        warning_mode = symplectic_newton_warning_mode

        !$acc parallel loop gang vector default(present) &
        !$acc&   copy(si(istart:iend), f(istart:iend)) copyin(ntau_macro) &
        !$acc&   copyout(loss_step(istart:iend), loss_time(istart:iend), &
        !$acc&           zend(:, istart:iend), nfev(istart:iend)) &
        !$acc&   private(it, ktau, ierr, lstep, elapsed, step_nfev) &
        !$acc&   firstprivate(warning_mode)
        do i = istart, iend
            ierr = 0
            lstep = ntimstep
            elapsed = 0.0_dp
            nfev(i) = 0
            macro: do it = 2, ntimstep
                do ktau = 1, ntau_macro(it)
                    call gpu_timestep_euler(si(i), f(i), warning_mode, ierr, &
                        step_nfev)
                    nfev(i) = nfev(i) + step_nfev
                    elapsed = elapsed + si(i)%dt
                    if (ierr /= 0) exit
                end do
                if (ierr /= 0) then
                    lstep = it
                    exit macro
                end if
            end do macro
            loss_step(i) = lstep
            if (ierr == 0 .and. lstep == ntimstep) &
                elapsed = real(sum(ntau_macro(2:ntimstep)), dp)*si(i)%dt
            loss_time(i) = elapsed
            zend(1, i) = si(i)%z(1)
            zend(2, i) = si(i)%z(2)
            zend(3, i) = si(i)%z(3)
            zend(4, i) = si(i)%z(4)
        end do
    end subroutine trace_orbits_gpu_range

    subroutine trace_orbits_gpu(si, f, npart, ntimstep, ntau_macro, loss_step, zend)
        ! Evolve npart pre-initialised particles, splitting the work evenly
        ! across all available NVIDIA GPUs (one host thread per device). Falls
        ! back to a single device / the host when only one is present.
        type(symplectic_integrator_t), intent(inout) :: si(npart)
        type(field_can_t), intent(inout) :: f(npart)
        integer, intent(in) :: npart, ntimstep
        integer, intent(in) :: ntau_macro(ntimstep)
        integer, intent(out) :: loss_step(npart)
        real(dp), intent(out) :: zend(4, npart)

        integer :: ngpu, navail, dev, i0, i1
        real(dp), allocatable :: loss_time(:)
        integer, allocatable :: nfev(:)
        character(16) :: env_val
        integer :: env_len, env_stat, req

        ! Number of devices to use. Default 1: with -gpu=mem:unified the spline
        ! coefficient array is shared and migrates between devices on access, so
        ! splitting across cards thrashes. Real multi-GPU needs a per-device
        ! resident copy of the read-only splines. Opt in with SIMPLE_GPU_NUM_DEVICES.
        req = 1
        call get_environment_variable('SIMPLE_GPU_NUM_DEVICES', env_val, env_len, env_stat)
        if (env_stat == 0 .and. env_len > 0) read (env_val, *, iostat=env_stat) req
        if (env_stat /= 0 .or. req < 1) req = 1

        ngpu = 1
        allocate (loss_time(npart), nfev(npart))
#ifdef _OPENACC
        navail = acc_get_num_devices(acc_device_nvidia)
        if (navail < 1) navail = 1
        ngpu = min(req, navail)
#endif
        if (ngpu <= 1) then
            call trace_orbits_gpu_range(si, f, npart, 1, npart, ntimstep, &
                                        ntau_macro, loss_step, loss_time, zend, nfev)
            return
        end if

        !$omp parallel num_threads(ngpu) private(dev, i0, i1)
        dev = omp_get_thread_num()
#ifdef _OPENACC
        call acc_set_device_num(dev, acc_device_nvidia)
#endif
        i0 = int(int(dev, 8)*npart/ngpu) + 1
        i1 = int(int(dev + 1, 8)*npart/ngpu)
        if (i1 >= i0) then
            call trace_orbits_gpu_range(si, f, npart, i0, i1, ntimstep, &
                                        ntau_macro, loss_step, loss_time, zend, nfev)
        end if
        !$omp end parallel
    end subroutine trace_orbits_gpu

    subroutine trace_orbits_gpu_method(si, f, npart, ntimstep, ntau_macro, &
            method, loss_step, loss_time, zend, nfev)
        ! Compile-time-specialized kernels sit behind this host-only dispatch;
        ! no method branch is executed in the per-stage device hot path.
        type(symplectic_integrator_t), intent(inout) :: si(npart)
        type(field_can_t), intent(inout) :: f(npart)
        integer, intent(in) :: npart, ntimstep, method
        integer, intent(in) :: ntau_macro(ntimstep)
        integer, intent(out) :: loss_step(npart), nfev(npart)
        real(dp), intent(out) :: loss_time(npart)
        real(dp), intent(out) :: zend(4, npart)

        select case (method)
        case (EXPL_IMPL_EULER)
            call trace_orbits_gpu_range(si, f, npart, 1, npart, ntimstep, &
                ntau_macro, loss_step, loss_time, zend, nfev)
        case (MIDPOINT)
            call trace_orbits_gpu_midpoint_range(si, f, npart, 1, npart, &
                ntimstep, ntau_macro, loss_step, loss_time, zend, nfev)
        case (CASH_KARP, DORMAND_PRINCE)
            call trace_orbits_gpu_rk54_range(si, f, npart, 1, npart, &
                ntimstep, ntau_macro, method, loss_step, loss_time, zend, nfev)
        case default
            error stop 'trace_orbits_gpu_method: unsupported GPU integrator'
        end select
    end subroutine trace_orbits_gpu_method

end module simple_gpu
