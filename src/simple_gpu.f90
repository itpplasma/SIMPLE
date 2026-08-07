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
    use util, only: pi, sqrt2, twopi
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
    use boozer_sub, only: boozer_state, splint_boozer_rk_device, &
        BOOZER_SECDERS_RADIAL_MIXED
    use boozer_rk_tables, only: rk_tables_ready, splint_boozer_rk_table_device
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
    public :: trace_orbits_gpu_landreman
    public :: evaluate_rhs_gpu

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

    subroutine gpu_get_derivatives2_radial_mixed(f, pphi)
        ! Euler's two-variable Newton system consumes only derivatives with at
        ! least one radial or canonical-momentum index. Avoid constructing the
        ! angular-only Hessian entries used by the five-variable midpoint map.
        !$acc routine seq
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: pphi

        call get_derivatives(f, pphi)

        f%d2vpar(1:3) = -f%d2Aph(1:3)/f%ro0 - f%d2hph(1:3)*f%vpar
        f%d2vpar(1) = f%d2vpar(1) - 2.0_dp*f%dhph(1)*f%dvpar(1)
        f%d2vpar(2) = f%d2vpar(2) - &
            (f%dhph(1)*f%dvpar(2) + f%dhph(2)*f%dvpar(1))
        f%d2vpar(3) = f%d2vpar(3) - &
            (f%dhph(1)*f%dvpar(3) + f%dhph(3)*f%dvpar(1))
        f%d2vpar(1:3) = f%d2vpar(1:3)/f%hph

        f%d2H(1:3) = f%vpar*f%d2vpar(1:3) + f%mu*f%d2Bmod(1:3)
        f%d2H(1) = f%d2H(1) + f%dvpar(1)**2
        f%d2H(2) = f%d2H(2) + f%dvpar(1)*f%dvpar(2)
        f%d2H(3) = f%d2H(3) + f%dvpar(1)*f%dvpar(3)

        f%d2pth(1:3) = f%d2vpar(1:3)*f%hth + &
            f%vpar*f%d2hth(1:3) + f%d2Ath(1:3)/f%ro0
        f%d2pth(1) = f%d2pth(1) + 2.0_dp*f%dvpar(1)*f%dhth(1)
        f%d2pth(2) = f%d2pth(2) + f%dvpar(1)*f%dhth(2) + &
            f%dvpar(2)*f%dhth(1)
        f%d2pth(3) = f%d2pth(3) + f%dvpar(1)*f%dhth(3) + &
            f%dvpar(3)*f%dhth(1)

        f%d2vpar(7:9) = -f%dhph/f%hph**2
        f%d2H(7:9) = f%dvpar(1:3)/f%hph + f%vpar*f%d2vpar(7:9)
        f%d2pth(7:9) = f%dhth/f%hph + f%hth*f%d2vpar(7:9)
    end subroutine gpu_get_derivatives2_radial_mixed

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

            call eval_field_booz_device(f, x(1), si%z(2), si%z(3), &
                BOOZER_SECDERS_RADIAL_MIXED)
            nfev = nfev + 1
            call gpu_get_derivatives2_radial_mixed(f, x(2))
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
        call eval_field_booz_device(f, x(1), x(2), x(3), 0)
        call get_derivatives(f, x(4))
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

    subroutine gpu_timestep_midpoint(si, f, warning_mode, ierr, nfev, predictor)
        !$acc routine seq
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        logical, intent(in) :: warning_mode
        integer, intent(out) :: ierr, nfev
        real(dp), intent(in), optional :: predictor(5)

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
            if (ktau == 1 .and. present(predictor)) then
                x = predictor
            else
                x(1:4) = si%z
                x(5) = si%z(1)
            end if
            call gpu_newton_midpoint(si, f, x, warning_mode, step_status, step_nfev)
            nfev = nfev + step_nfev
            if (step_status /= SYMPLECTIC_STEP_OK .and. ktau == 1 .and. &
                    present(predictor)) then
                si = accepted_integrator
                f = accepted_field
                x(1:4) = si%z
                x(5) = si%z(1)
                call gpu_newton_midpoint(si, f, x, warning_mode, step_status, &
                    step_nfev)
                nfev = nfev + step_nfev
            end if
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

    subroutine gpu_rhs_canonical(y, mu, ro0, dydt)
        ! Canonical guiding-centre equations used by SIMPLE's CPU RK45 path,
        ! with the Boozer evaluator bound statically for device compilation.
        !$acc routine seq
        real(dp), intent(in) :: y(4)
        real(dp), intent(in) :: mu, ro0
        real(dp), intent(out) :: dydt(4)

        real(dp) :: aphi, daphi, btheta, dbtheta, bphi, dbphi, bmod
        real(dp) :: dbmod(3), hth, hph, vpar, hprime, bmod2
        real(dp) :: dhth_r, dhth_phi
        real(dp) :: dhph_r, dhph_theta, dhph_phi
        real(dp) :: dvpar_r, dvpar_theta, dvpar_phi
        real(dp) :: dh_r, dh_theta, dh_phi, dpth_r, dpth_phi

        if (rk_tables_ready) then
            call splint_boozer_rk_table_device(y(1), y(2), y(3), aphi, &
                daphi, btheta, dbtheta, bphi, dbphi, bmod, dbmod)
        else
            call splint_boozer_rk_device(y(1), y(2), y(3), aphi, daphi, &
                btheta, dbtheta, bphi, dbphi, bmod, dbmod)
        end if
        bmod2 = bmod*bmod
        hth = btheta/bmod
        hph = bphi/bmod
        dhth_r = dbtheta/bmod - btheta*dbmod(1)/bmod2
        dhth_phi = -btheta*dbmod(3)/bmod2
        dhph_r = dbphi/bmod - bphi*dbmod(1)/bmod2
        dhph_theta = -bphi*dbmod(2)/bmod2
        dhph_phi = -bphi*dbmod(3)/bmod2
        vpar = (y(4) - aphi/ro0)/hph
        dvpar_r = -(daphi/ro0 + dhph_r*vpar)/hph
        dvpar_theta = -dhph_theta*vpar/hph
        dvpar_phi = -dhph_phi*vpar/hph
        dh_r = vpar*dvpar_r + mu*dbmod(1)
        dh_theta = vpar*dvpar_theta + mu*dbmod(2)
        dh_phi = vpar*dvpar_phi + mu*dbmod(3)
        dpth_r = dvpar_r*hth + vpar*dhth_r + boozer_state%torflux/ro0
        dpth_phi = dvpar_phi*hth + vpar*dhth_phi
        hprime = dh_r/dpth_r

        dydt(1) = -(dh_theta - hth/hph*dh_phi)/dpth_r
        dydt(2) = hprime
        dydt(3) = (vpar - hprime*hth)/hph
        dydt(4) = -(dh_phi - hprime*dpth_phi)
    end subroutine gpu_rhs_canonical

    subroutine gpu_rhs_canonical_cartesian(y, mu, ro0, dydt)
        ! Integrate the radial Boozer plane as (s*cos(theta), s*sin(theta)).
        ! This removes the theta singularity at the magnetic axis while keeping
        ! SIMPLE's native canonical momentum and canonical Boozer RHS.
        !$acc routine seq
        real(dp), intent(in) :: y(4), mu, ro0
        real(dp), intent(out) :: dydt(4)

        real(dp) :: polar(4), polar_rhs(4), radial, radial_x, radial_y

        radial = sqrt(y(1)*y(1) + y(2)*y(2))
        polar = [radial, atan2(y(2), y(1)), y(3), y(4)]
        call gpu_rhs_canonical(polar, mu, ro0, polar_rhs)
        if (radial > tiny(1.0_dp)) then
            radial_x = y(1)/radial
            radial_y = y(2)/radial
        else
            radial_x = cos(polar(2))
            radial_y = sin(polar(2))
        end if
        dydt(1) = polar_rhs(1)*radial_x - y(2)*polar_rhs(2)
        dydt(2) = polar_rhs(1)*radial_y + y(1)*polar_rhs(2)
        dydt(3:4) = polar_rhs(3:4)
    end subroutine gpu_rhs_canonical_cartesian

    subroutine evaluate_rhs_gpu(y, mu, ro0, npart, dydt)
        ! Batch entry point used to compare the compact device RHS against the
        ! host canonical-field implementation on independently evaluated data.
        integer, intent(in) :: npart
        real(dp), intent(in) :: y(4, npart), mu(npart), ro0(npart)
        real(dp), intent(out) :: dydt(4, npart)

        integer :: i

        !$acc parallel loop gang vector copyin(y, mu, ro0) copyout(dydt)
        do i = 1, npart
            call gpu_rhs_canonical(y(:, i), mu(i), ro0(i), dydt(:, i))
        end do
    end subroutine evaluate_rhs_gpu

    subroutine gpu_trace_rk54(mu, ro0, y, duration, controls, loss_time, &
            nfev, ierr)
        !$acc routine seq
        real(dp), intent(in) :: mu, ro0
        real(dp), intent(inout) :: y(4)
        real(dp), intent(in) :: duration
        type(rk54_controls4_t), intent(in) :: controls
        real(dp), intent(out) :: loss_time
        integer, intent(out) :: nfev, ierr

        integer, parameter :: max_attempts = 100000000
        type(rk54_state4_t) :: state
        real(dp) :: t_eval, y_eval(4), dydt(4), h0
        integer :: request, attempt, rhs_calls

        loss_time = duration
        nfev = 0
        ierr = 0
        if (duration <= 0.0_dp) return

        h0 = min(duration, 1.0e-3_dp*controls%hmax)
        call rk54_initialize4(state, 0.0_dp, y, h0)

        do attempt = 1, max_attempts
            if (state%t >= duration) exit
            if (state%y(1) >= 1.0_dp) then
                loss_time = state%t
                exit
            end if

            state%h = min(state%h, duration - state%t)
            call rk54_request4(state, controls, t_eval, y_eval, request)
            rhs_calls = 0
            do while (request == RK54_NEED_RHS)
                call gpu_rhs_canonical(y_eval, mu, ro0, dydt)
                call rk54_supply4(state, controls, dydt, t_eval, y_eval, request)
                rhs_calls = rhs_calls + 1
                if (rhs_calls > 7) then
                    request = RK54_FAILED
                    exit
                end if
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

    subroutine gpu_trace_rk54_cartesian(mu, ro0, y, duration, controls, &
            loss_time, nfev, ierr)
        !$acc routine seq
        real(dp), intent(in) :: mu, ro0
        real(dp), intent(inout) :: y(4)
        real(dp), intent(in) :: duration
        type(rk54_controls4_t), intent(in) :: controls
        real(dp), intent(out) :: loss_time
        integer, intent(out) :: nfev, ierr

        integer, parameter :: max_attempts = 100000000
        type(rk54_state4_t) :: state
        real(dp) :: t_eval, y_eval(4), dydt(4), h0, radial
        integer :: request, attempt, rhs_calls

        loss_time = duration
        nfev = 0
        ierr = 0
        if (duration <= 0.0_dp) return

        h0 = min(duration, 1.0e-3_dp*controls%hmax)
        call rk54_initialize4(state, 0.0_dp, y, h0)

        do attempt = 1, max_attempts
            if (state%t >= duration) exit
            radial = sqrt(state%y(1)*state%y(1) + state%y(2)*state%y(2))
            if (radial >= 1.0_dp) then
                loss_time = state%t
                exit
            end if

            state%h = min(state%h, duration - state%t)
            call rk54_request4(state, controls, t_eval, y_eval, request)
            rhs_calls = 0
            do while (request == RK54_NEED_RHS)
                call gpu_rhs_canonical_cartesian(y_eval, mu, ro0, dydt)
                call rk54_supply4(state, controls, dydt, t_eval, y_eval, request)
                rhs_calls = rhs_calls + 1
                if (rhs_calls > 7) then
                    request = RK54_FAILED
                    exit
                end if
            end do

            if (request == RK54_FAILED) then
                ierr = RK54_FAILED
                exit
            end if
            if (request == RK54_ACCEPTED) then
                radial = sqrt(state%y(1)*state%y(1) + state%y(2)*state%y(2))
                if (radial >= 1.0_dp) then
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
    end subroutine gpu_trace_rk54_cartesian

    subroutine gpu_rk_segment_hmax(y, duration, hmax)
        ! Match FIRM3D's per-segment quarter-turn cap, evaluated at the
        ! segment's initial state: dt <= (G/B)*pi/(2*v_total). SIMPLE's
        ! canonical independent variable is physical_time*v_total/sqrt(2).
        !$acc routine seq
        real(dp), intent(in) :: y(4), duration
        real(dp), intent(out) :: hmax

        real(dp) :: aphi, daphi, btheta, dbtheta, bphi, dbphi, bmod
        real(dp) :: dbmod(3)

        if (rk_tables_ready) then
            call splint_boozer_rk_table_device(y(1), y(2), y(3), aphi, &
                daphi, btheta, dbtheta, bphi, dbphi, bmod, dbmod)
        else
            call splint_boozer_rk_device(y(1), y(2), y(3), aphi, daphi, &
                btheta, dbtheta, bphi, dbphi, bmod, dbmod)
        end if
        hmax = min(duration, 0.5_dp*pi*abs(bphi/bmod)/sqrt2)
    end subroutine gpu_rk_segment_hmax

    subroutine gpu_rk_segment_hmax_cartesian(y, duration, hmax, &
            momentum_atol_scale)
        !$acc routine seq
        real(dp), intent(in) :: y(4), duration
        real(dp), intent(out) :: hmax, momentum_atol_scale

        real(dp) :: polar(4), aphi, daphi, btheta, dbtheta
        real(dp) :: bphi, dbphi, bmod, dbmod(3)

        polar = [sqrt(y(1)*y(1) + y(2)*y(2)), atan2(y(2), y(1)), &
            y(3), y(4)]
        if (rk_tables_ready) then
            call splint_boozer_rk_table_device(polar(1), polar(2), polar(3), &
                aphi, daphi, btheta, dbtheta, bphi, dbphi, bmod, dbmod)
        else
            call splint_boozer_rk_device(polar(1), polar(2), polar(3), aphi, &
                daphi, btheta, dbtheta, bphi, dbphi, bmod, dbmod)
        end if
        hmax = min(duration, 0.5_dp*pi*abs(bphi/bmod)/sqrt2)
        ! FIRM3D applies atol=rtol*v_total to v_parallel. In SIMPLE's
        ! canonical momentum p_phi=v_parallel*h_phi+A_phi/ro0, the equivalent
        ! absolute scale is sqrt(2)*|h_phi|.
        momentum_atol_scale = sqrt2*abs(bphi/bmod)
    end subroutine gpu_rk_segment_hmax_cartesian

    subroutine gpu_trace_euler_segment(si, f, nsteps, warning_mode, elapsed, &
            nfev, ierr)
        !$acc routine seq
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(in) :: nsteps
        logical, intent(in) :: warning_mode
        real(dp), intent(out) :: elapsed
        integer, intent(out) :: nfev, ierr

        integer :: k, step_nfev

        elapsed = 0.0_dp
        nfev = 0
        ierr = SYMPLECTIC_STEP_OK
        do k = 1, nsteps
            call gpu_timestep_euler(si, f, warning_mode, ierr, step_nfev)
            nfev = nfev + step_nfev
            elapsed = elapsed + si%dt
            if (ierr /= SYMPLECTIC_STEP_OK) exit
        end do
    end subroutine gpu_trace_euler_segment

    subroutine gpu_trace_midpoint_segment(si, f, nsteps, warning_mode, elapsed, &
            nfev, ierr)
        !$acc routine seq
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(in) :: nsteps
        logical, intent(in) :: warning_mode
        real(dp), intent(out) :: elapsed
        integer, intent(out) :: nfev, ierr

        integer :: k, step_nfev
        real(dp) :: z_previous(4), predictor(5)
        logical :: use_predictor

        elapsed = 0.0_dp
        nfev = 0
        ierr = SYMPLECTIC_STEP_OK
        use_predictor = .false.
        do k = 1, nsteps
            if (use_predictor) then
                predictor(1:4) = 2.0_dp*si%z - z_previous
                predictor(5) = 1.5_dp*si%z(1) - 0.5_dp*z_previous(1)
                use_predictor = predictor(1) >= 0.0_dp .and. &
                    predictor(1) <= 1.0_dp .and. predictor(5) >= 0.0_dp .and. &
                    predictor(5) <= 1.0_dp .and. &
                    abs(si%z(2) - z_previous(2)) < 0.5_dp*pi
            end if
            z_previous = si%z
            if (use_predictor) then
                call gpu_timestep_midpoint(si, f, warning_mode, ierr, step_nfev, &
                    predictor)
            else
                call gpu_timestep_midpoint(si, f, warning_mode, ierr, step_nfev)
            end if
            nfev = nfev + step_nfev
            elapsed = elapsed + si%dt
            if (ierr /= SYMPLECTIC_STEP_OK) exit
            use_predictor = .true.
        end do
    end subroutine gpu_trace_midpoint_segment

    subroutine trace_orbits_gpu_rk54_landreman(si, f, npart, method, &
            total_duration, block_duration, time_scale, loss_tau, maxloss, &
            hmin, loss_step, loss_time, zend, nfev, observed_duration, &
            energy_loss_fraction, warp_nfev_slots)
        type(symplectic_integrator_t), intent(inout) :: si(npart)
        type(field_can_t), intent(in) :: f(npart)
        integer, intent(in) :: npart, method
        real(dp), intent(in) :: total_duration, block_duration, time_scale
        real(dp), intent(in) :: loss_tau, maxloss, hmin
        integer, intent(out) :: loss_step(npart), nfev(npart)
        real(dp), intent(out) :: loss_time(npart), zend(4, npart)
        real(dp), intent(out) :: observed_duration, energy_loss_fraction
        integer(int64), intent(out) :: warp_nfev_slots

        type(rk54_controls4_t) :: controls, controls_i
        real(dp), allocatable :: y(:, :), mu(:), ro0(:), loss_time_slot(:)
        real(dp), allocatable :: y_next(:, :), mu_next(:), ro0_next(:)
        real(dp), allocatable :: loss_time_next(:), work_key(:)
        integer, allocatable :: orbit_status(:), orbit_status_next(:)
        integer, allocatable :: nfev_slot(:), nfev_next(:), segment_work(:)
        integer, allocatable :: original_index(:), original_index_next(:)
        integer, allocatable :: order(:), order_scratch(:)
        real(dp) :: current_time, segment_duration, segment_loss_time
        real(dp) :: y_local(4), hmax_local, momentum_atol_scale
        real(dp) :: new_weight, weighted_losses
        real(dp) :: loss_detection_tolerance, pilot_fraction, pilot_duration
        character(16) :: ordering_value, pilot_value
        integer :: i, ierr, segment_nfev, active_count, original
        integer :: ordering_length, ordering_status, pilot_length, pilot_status
        logical :: order_work

        if (block_duration <= 0.0_dp .or. time_scale <= 0.0_dp .or. &
                loss_tau <= 0.0_dp) &
            error stop 'invalid Landreman RK segment controls'

        allocate (y(4, npart), mu(npart), ro0(npart), &
            loss_time_slot(npart), y_next(4, npart), mu_next(npart), &
            ro0_next(npart), loss_time_next(npart), work_key(npart), &
            orbit_status(npart), orbit_status_next(npart), nfev_slot(npart), &
            nfev_next(npart), segment_work(npart), original_index(npart), &
            original_index_next(npart), order(npart), order_scratch(npart))
        orbit_status = 0
        loss_time_slot = total_duration
        weighted_losses = 0.0_dp
        do i = 1, npart
            y(1, i) = si(i)%z(1)*cos(si(i)%z(2))
            y(2, i) = si(i)%z(1)*sin(si(i)%z(2))
            y(3:4, i) = si(i)%z(3:4)
            mu(i) = f(i)%mu
            ro0(i) = f(i)%ro0
            original_index(i) = i
            if (sqrt(y(1, i)*y(1, i) + y(2, i)*y(2, i)) >= &
                    1.0_dp - 64.0_dp*epsilon(1.0_dp)) then
                orbit_status(i) = 1
                loss_time_slot(i) = 0.0_dp
                weighted_losses = weighted_losses + 1.0_dp
            end if
        end do
        nfev_slot = 0
        segment_work = 0
        warp_nfev_slots = 0_int64
        current_time = 0.0_dp
        loss_detection_tolerance = 1.0e-12_dp/time_scale

        order_work = .true.
        call get_environment_variable('SIMPLE_GPU_WORK_ORDER', ordering_value, &
            ordering_length, ordering_status)
        if (ordering_status == 0 .and. ordering_length > 0) &
            order_work = trim(adjustl(ordering_value(:ordering_length))) /= '0' .and. &
                trim(adjustl(ordering_value(:ordering_length))) /= 'none'
#ifndef _OPENACC
        order_work = .false.
#endif
        pilot_fraction = 0.05_dp
        call get_environment_variable('SIMPLE_GPU_PILOT_FRACTION', pilot_value, &
            pilot_length, pilot_status)
        if (pilot_status == 0 .and. pilot_length > 0) then
            read (pilot_value(:pilot_length), *, iostat=pilot_status) pilot_fraction
            if (pilot_status /= 0) &
                error stop 'invalid SIMPLE_GPU_PILOT_FRACTION'
        end if
        if (pilot_fraction < 0.0_dp .or. pilot_fraction > 0.25_dp) &
            error stop 'SIMPLE_GPU_PILOT_FRACTION must be in [0, 0.25]'
        if (.not. order_work .or. npart < 1024) pilot_fraction = 0.0_dp
        if (order_work) then
            do i = 1, npart
                if (orbit_status(i) == 0) then
                    work_key(i) = abs(f(i)%vpar)
                else
                    work_key(i) = -1.0_dp
                end if
            end do
            call stable_order_descending(work_key, order, order_scratch)
            call reorder_rk_slots(order, y, mu, ro0, orbit_status, &
                loss_time_slot, nfev_slot, original_index, y_next, mu_next, &
                ro0_next, orbit_status_next, loss_time_next, nfev_next, &
                original_index_next)
        end if
        active_count = count(orbit_status == 0)

        controls%rtol = si(1)%rtol
        controls%atol = si(1)%rtol
        controls%hmin = hmin
        if (method == CASH_KARP) then
            controls%method = RK54_CASH_KARP
        else
            controls%method = RK54_DORMAND_PRINCE
        end if

        !$acc data copy(y, orbit_status, loss_time_slot, nfev_slot, segment_work) &
        !$acc& copyin(mu, ro0)
        if (pilot_fraction > 0.0_dp .and. active_count > 0) then
            pilot_duration = pilot_fraction*min(block_duration, total_duration)
#ifdef _OPENACC
            !$acc parallel loop gang vector default(present) &
            !$acc&   private(controls_i, y_local, segment_loss_time, &
            !$acc&           hmax_local, momentum_atol_scale, segment_nfev, ierr) &
            !$acc&   firstprivate(controls, pilot_duration)
#else
            !$omp parallel do default(shared) schedule(static) &
            !$omp&   private(i, controls_i, y_local, segment_loss_time, &
            !$omp&           hmax_local, momentum_atol_scale, segment_nfev, ierr)
#endif
            do i = 1, active_count
                y_local = y(:, i)
                controls_i = controls
                call gpu_rk_segment_hmax_cartesian(y_local, pilot_duration, &
                    hmax_local, momentum_atol_scale)
                controls_i%hmax = hmax_local
                controls_i%atol(4) = controls_i%rtol*momentum_atol_scale
                call gpu_trace_rk54_cartesian(mu(i), ro0(i), y_local, &
                    pilot_duration, controls_i, segment_loss_time, &
                    segment_nfev, ierr)
                segment_work(i) = segment_nfev + 1
            end do
#ifndef _OPENACC
            !$omp end parallel do
#endif
            !$acc update self(segment_work)
            nfev_slot(1:active_count) = segment_work(1:active_count)
            do i = 1, active_count, 32
                warp_nfev_slots = warp_nfev_slots + 32_int64*int( &
                    maxval(segment_work(i:min(i + 31, active_count))), int64)
            end do
            do i = 1, npart
                if (orbit_status(i) == 0) then
                    work_key(i) = real(segment_work(i), dp)
                else
                    work_key(i) = -1.0_dp
                end if
            end do
            call stable_order_descending(work_key, order, order_scratch)
            call reorder_rk_slots(order, y, mu, ro0, orbit_status, &
                loss_time_slot, nfev_slot, original_index, y_next, mu_next, &
                ro0_next, orbit_status_next, loss_time_next, nfev_next, &
                original_index_next)
            !$acc update device(y, mu, ro0, orbit_status, loss_time_slot, &
            !$acc& nfev_slot)
        end if
        do while (current_time < total_duration)
            segment_duration = min(block_duration, total_duration - current_time)
            if (weighted_losses/real(npart, dp) > maxloss) then
                current_time = current_time + segment_duration
                exit
            end if
            new_weight = 0.0_dp
#ifdef _OPENACC
            !$acc parallel loop gang vector default(present) &
            !$acc&   private(controls_i, y_local, segment_loss_time, &
            !$acc&           hmax_local, momentum_atol_scale, segment_nfev, ierr) &
            !$acc&   firstprivate(controls, current_time, segment_duration, &
            !$acc&                time_scale, loss_tau, loss_detection_tolerance) &
            !$acc&   copy(new_weight)
#else
            !$omp parallel do default(shared) schedule(static) &
            !$omp&   private(i, controls_i, y_local, segment_loss_time, &
            !$omp&           hmax_local, momentum_atol_scale, segment_nfev, ierr) &
            !$omp&   reduction(+:new_weight)
#endif
            do i = 1, active_count
                segment_work(i) = 0
                if (orbit_status(i) /= 0) cycle
                y_local = y(:, i)
                if (sqrt(y_local(1)*y_local(1) + y_local(2)*y_local(2)) >= &
                        1.0_dp - 64.0_dp*epsilon(1.0_dp)) then
                    orbit_status(i) = 1
                    loss_time_slot(i) = current_time
                    !$acc atomic update
                    new_weight = new_weight + &
                        exp(-loss_time_slot(i)*time_scale/loss_tau)
                    cycle
                end if
                controls_i = controls
                call gpu_rk_segment_hmax_cartesian(y_local, segment_duration, &
                    hmax_local, momentum_atol_scale)
                controls_i%hmax = hmax_local
                controls_i%atol(4) = controls_i%rtol*momentum_atol_scale
                call gpu_trace_rk54_cartesian(mu(i), ro0(i), y_local, &
                    segment_duration, controls_i, segment_loss_time, &
                    segment_nfev, ierr)
                y(:, i) = y_local
                segment_work(i) = segment_nfev + 1
                nfev_slot(i) = nfev_slot(i) + segment_work(i)
                if (ierr /= 0) then
                    orbit_status(i) = -ierr
                else if (segment_loss_time < segment_duration - &
                        loss_detection_tolerance) then
                    orbit_status(i) = 1
                    loss_time_slot(i) = current_time + segment_loss_time
                    !$acc atomic update
                    new_weight = new_weight + &
                        exp(-loss_time_slot(i)*time_scale/loss_tau)
                end if
            end do
#ifndef _OPENACC
            !$omp end parallel do
#endif
            !$acc update self(segment_work)
            do i = 1, active_count, 32
                warp_nfev_slots = warp_nfev_slots + 32_int64*int( &
                    maxval(segment_work(i:min(i + 31, active_count))), int64)
            end do
            weighted_losses = weighted_losses + new_weight
            current_time = current_time + segment_duration
            if (current_time < total_duration .and. &
                    weighted_losses/real(npart, dp) > maxloss) exit
            if (current_time < total_duration .and. order_work) then
                !$acc update self(y, orbit_status, loss_time_slot, nfev_slot)
                do i = 1, npart
                    if (orbit_status(i) == 0) then
                        work_key(i) = real(segment_work(i), dp)
                    else
                        work_key(i) = -1.0_dp
                    end if
                end do
                call stable_order_descending(work_key, order, order_scratch)
                call reorder_rk_slots(order, y, mu, ro0, orbit_status, &
                    loss_time_slot, nfev_slot, original_index, y_next, mu_next, &
                    ro0_next, orbit_status_next, loss_time_next, nfev_next, &
                    original_index_next)
                active_count = count(orbit_status == 0)
                !$acc update device(y, mu, ro0, orbit_status, loss_time_slot, &
                !$acc& nfev_slot)
            end if
        end do
        !$acc end data

        do i = 1, npart
            original = original_index(i)
            si(original)%z(1) = sqrt(y(1, i)*y(1, i) + y(2, i)*y(2, i))
            si(original)%z(2) = atan2(y(2, i), y(1, i))
            si(original)%z(3:4) = y(3:4, i)
            zend(:, original) = si(original)%z
            loss_step(original) = orbit_status(i)
            loss_time(original) = loss_time_slot(i)
            nfev(original) = nfev_slot(i)
        end do
        observed_duration = current_time
        energy_loss_fraction = weighted_losses/real(npart, dp)
    end subroutine trace_orbits_gpu_rk54_landreman

    subroutine stable_order_descending(key, order, scratch)
        real(dp), intent(in) :: key(:)
        integer, intent(out) :: order(size(key))
        integer, intent(inout) :: scratch(size(key))

        integer :: n, width, left, middle, right, i, j, k

        n = size(key)
        do i = 1, n
            order(i) = i
        end do
        width = 1
        do while (width < n)
            do left = 1, n, 2*width
                middle = min(left + width, n + 1)
                right = min(left + 2*width - 1, n)
                i = left
                j = middle
                do k = left, right
                    if (i >= middle) then
                        scratch(k) = order(j)
                        j = j + 1
                    else if (j > right) then
                        scratch(k) = order(i)
                        i = i + 1
                    else if (key(order(i)) >= key(order(j))) then
                        scratch(k) = order(i)
                        i = i + 1
                    else
                        scratch(k) = order(j)
                        j = j + 1
                    end if
                end do
            end do
            order = scratch
            width = 2*width
        end do
    end subroutine stable_order_descending

    subroutine reorder_rk_slots(order, y, mu, ro0, orbit_status, loss_time, &
            nfev, original_index, y_next, mu_next, ro0_next, &
            orbit_status_next, loss_time_next, nfev_next, original_index_next)
        integer, intent(in) :: order(:)
        real(dp), intent(inout) :: y(:, :), mu(:), ro0(:), loss_time(:)
        integer, intent(inout) :: orbit_status(:), nfev(:), original_index(:)
        real(dp), intent(inout) :: y_next(:, :), mu_next(:), ro0_next(:)
        real(dp), intent(inout) :: loss_time_next(:)
        integer, intent(inout) :: orbit_status_next(:), nfev_next(:)
        integer, intent(inout) :: original_index_next(:)

        integer :: i, source

        do i = 1, size(order)
            source = order(i)
            y_next(:, i) = y(:, source)
            mu_next(i) = mu(source)
            ro0_next(i) = ro0(source)
            orbit_status_next(i) = orbit_status(source)
            loss_time_next(i) = loss_time(source)
            nfev_next(i) = nfev(source)
            original_index_next(i) = original_index(source)
        end do
        y = y_next
        mu = mu_next
        ro0 = ro0_next
        orbit_status = orbit_status_next
        loss_time = loss_time_next
        nfev = nfev_next
        original_index = original_index_next
    end subroutine reorder_rk_slots

    subroutine trace_orbits_gpu_symplectic_landreman(si, f, npart, ntimstep, &
            ntau_macro, method, block_duration, time_scale, loss_tau, maxloss, &
            loss_step, loss_time, zend, nfev, observed_duration, &
            energy_loss_fraction, warp_nfev_slots)
        type(symplectic_integrator_t), intent(inout) :: si(npart)
        type(field_can_t), intent(inout) :: f(npart)
        integer, intent(in) :: npart, ntimstep, method
        integer, intent(in) :: ntau_macro(ntimstep)
        real(dp), intent(in) :: block_duration, time_scale, loss_tau, maxloss
        integer, intent(out) :: loss_step(npart), nfev(npart)
        real(dp), intent(out) :: loss_time(npart), zend(4, npart)
        real(dp), intent(out) :: observed_duration, energy_loss_fraction
        integer(int64), intent(out) :: warp_nfev_slots

        integer, allocatable :: orbit_status(:)
        real(dp) :: total_duration, current_time, segment_duration
        real(dp) :: elapsed, new_weight, weighted_losses, expected_duration
        real(dp) :: schedule_tolerance
        integer :: i, it, ierr, segment_nfev
        logical :: warning_mode

        if (block_duration <= 0.0_dp .or. time_scale <= 0.0_dp .or. &
                loss_tau <= 0.0_dp) &
            error stop 'invalid Landreman symplectic segment controls'
        if (method /= EXPL_IMPL_EULER .and. method /= MIDPOINT) &
            error stop 'invalid Landreman symplectic method'

        total_duration = si(1)%dt*real(sum(ntau_macro(2:ntimstep)), dp)
        schedule_tolerance = 64.0_dp*epsilon(1.0_dp)*max(1.0_dp, total_duration)
        current_time = 0.0_dp
        do it = 2, ntimstep
            segment_duration = si(1)%dt*real(ntau_macro(it), dp)
            expected_duration = min(block_duration, total_duration - current_time)
            if (abs(segment_duration - expected_duration) > schedule_tolerance) &
                error stop 'SIMPLE macrostep schedule does not match Landreman block size'
            current_time = current_time + segment_duration
        end do

        allocate (orbit_status(npart))
        orbit_status = 0
        loss_time = total_duration
        nfev = 0
        warp_nfev_slots = 0_int64
        weighted_losses = 0.0_dp
        do i = 1, npart
            if (si(i)%z(1) >= 1.0_dp - 64.0_dp*epsilon(1.0_dp)) then
                orbit_status(i) = 1
                loss_time(i) = 0.0_dp
                weighted_losses = weighted_losses + 1.0_dp
            end if
        end do
        current_time = 0.0_dp
        warning_mode = symplectic_newton_warning_mode

        !$acc data copy(si, f, orbit_status, loss_time, nfev)
        do it = 2, ntimstep
            segment_duration = si(1)%dt*real(ntau_macro(it), dp)
            if (weighted_losses/real(npart, dp) > maxloss) then
                current_time = current_time + segment_duration
                exit
            end if
            new_weight = 0.0_dp
            if (method == EXPL_IMPL_EULER) then
#ifdef _OPENACC
                !$acc parallel loop gang vector default(present) &
                !$acc&   private(elapsed, segment_nfev, ierr) &
                !$acc&   firstprivate(warning_mode, current_time, time_scale, loss_tau) &
                !$acc&   copy(new_weight)
#else
                !$omp parallel do default(shared) schedule(static) &
                !$omp&   private(i, elapsed, segment_nfev, ierr) &
                !$omp&   reduction(+:new_weight)
#endif
                do i = 1, npart
                    if (orbit_status(i) /= 0) cycle
                    call gpu_trace_euler_segment(si(i), f(i), ntau_macro(it), &
                        warning_mode, elapsed, segment_nfev, ierr)
                    nfev(i) = nfev(i) + segment_nfev
                    if (ierr == SYMPLECTIC_STEP_OUTSIDE_DOMAIN) then
                        orbit_status(i) = 1
                        loss_time(i) = current_time + elapsed
                        !$acc atomic update
                        new_weight = new_weight + &
                            exp(-loss_time(i)*time_scale/loss_tau)
                    else if (ierr /= SYMPLECTIC_STEP_OK) then
                        orbit_status(i) = -ierr
                    end if
                end do
#ifndef _OPENACC
                !$omp end parallel do
#endif
            else
#ifdef _OPENACC
                !$acc parallel loop gang vector default(present) &
                !$acc&   private(elapsed, segment_nfev, ierr) &
                !$acc&   firstprivate(warning_mode, current_time, time_scale, loss_tau) &
                !$acc&   copy(new_weight)
#else
                !$omp parallel do default(shared) schedule(static) &
                !$omp&   private(i, elapsed, segment_nfev, ierr) &
                !$omp&   reduction(+:new_weight)
#endif
                do i = 1, npart
                    if (orbit_status(i) /= 0) cycle
                    call gpu_trace_midpoint_segment(si(i), f(i), ntau_macro(it), &
                        warning_mode, elapsed, segment_nfev, ierr)
                    nfev(i) = nfev(i) + segment_nfev
                    if (ierr == SYMPLECTIC_STEP_OUTSIDE_DOMAIN) then
                        orbit_status(i) = 1
                        loss_time(i) = current_time + elapsed
                        !$acc atomic update
                        new_weight = new_weight + &
                            exp(-loss_time(i)*time_scale/loss_tau)
                    else if (ierr /= SYMPLECTIC_STEP_OK) then
                        orbit_status(i) = -ierr
                    end if
                end do
#ifndef _OPENACC
                !$omp end parallel do
#endif
            end if
            weighted_losses = weighted_losses + new_weight
            current_time = current_time + segment_duration
            if (current_time < total_duration .and. &
                    weighted_losses/real(npart, dp) > maxloss) exit
        end do
        !$acc end data

        do i = 1, npart, 32
            warp_nfev_slots = warp_nfev_slots + 32_int64*int( &
                maxval(nfev(i:min(i + 31, npart))), int64)
        end do

        do i = 1, npart
            zend(:, i) = si(i)%z
            loss_step(i) = orbit_status(i)
        end do
        observed_duration = current_time
        energy_loss_fraction = weighted_losses/real(npart, dp)
    end subroutine trace_orbits_gpu_symplectic_landreman

    subroutine trace_orbits_gpu_landreman(si, f, npart, ntimstep, ntau_macro, &
            method, block_duration, time_scale, loss_tau, maxloss, hmin, &
            loss_step, loss_time, zend, nfev, observed_duration, &
            energy_loss_fraction, warp_nfev_slots)
        type(symplectic_integrator_t), intent(inout) :: si(npart)
        type(field_can_t), intent(inout) :: f(npart)
        integer, intent(in) :: npart, ntimstep, method
        integer, intent(in) :: ntau_macro(ntimstep)
        real(dp), intent(in) :: block_duration, time_scale, loss_tau, maxloss, hmin
        integer, intent(out) :: loss_step(npart), nfev(npart)
        real(dp), intent(out) :: loss_time(npart), zend(4, npart)
        real(dp), intent(out) :: observed_duration, energy_loss_fraction
        integer(int64), intent(out) :: warp_nfev_slots

        real(dp) :: total_duration

        total_duration = si(1)%dt*real(sum(ntau_macro(2:ntimstep)), dp)
        select case (method)
        case (CASH_KARP, DORMAND_PRINCE)
            call trace_orbits_gpu_rk54_landreman(si, f, npart, method, &
                total_duration, block_duration, time_scale, loss_tau, maxloss, &
                hmin, loss_step, loss_time, zend, nfev, observed_duration, &
                energy_loss_fraction, warp_nfev_slots)
        case (EXPL_IMPL_EULER, MIDPOINT)
            call trace_orbits_gpu_symplectic_landreman(si, f, npart, ntimstep, &
                ntau_macro, method, block_duration, time_scale, loss_tau, maxloss, &
                loss_step, loss_time, zend, nfev, observed_duration, &
                energy_loss_fraction, warp_nfev_slots)
        case default
            error stop 'trace_orbits_gpu_landreman: unsupported integrator'
        end select
    end subroutine trace_orbits_gpu_landreman

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
        real(dp) :: mu(npart), ro0(npart)
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

        do i = istart, iend
            mu(i) = f(i)%mu
            ro0(i) = f(i)%ro0
        end do

        !$acc parallel loop gang vector default(present) &
        !$acc&   copyin(si(istart:iend), mu(istart:iend), ro0(istart:iend), &
        !$acc&          ntau_macro) &
        !$acc&   copyout(loss_step(istart:iend), zend(:, istart:iend), &
        !$acc&           nfev(istart:iend), loss_time_out(istart:iend)) &
        !$acc&   private(y, loss_time, elapsed, ierr, it) firstprivate(controls, duration)
        do i = istart, iend
            y = si(i)%z
            call gpu_trace_rk54(mu(i), ro0(i), y, duration, controls, loss_time, &
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
        real(dp) :: elapsed, z_previous(4), predictor(5)
        logical :: warning_mode, use_predictor

        warning_mode = symplectic_newton_warning_mode
        !$acc parallel loop gang vector default(present) &
        !$acc&   copy(si(istart:iend), f(istart:iend)) copyin(ntau_macro) &
        !$acc&   copyout(loss_step(istart:iend), zend(:, istart:iend), &
        !$acc&           nfev(istart:iend), loss_time_out(istart:iend)) &
        !$acc&   private(it, ktau, ierr, lstep, step_nfev, elapsed, &
        !$acc&           z_previous, predictor, use_predictor) &
        !$acc&   firstprivate(warning_mode)
        do i = istart, iend
            ierr = 0
            lstep = ntimstep
            nfev(i) = 0
            elapsed = 0.0_dp
            use_predictor = .false.
            macro: do it = 2, ntimstep
                do ktau = 1, ntau_macro(it)
                    if (use_predictor) then
                        predictor(1:4) = 2.0_dp*si(i)%z - z_previous
                        predictor(5) = 1.5_dp*si(i)%z(1) - &
                            0.5_dp*z_previous(1)
                        use_predictor = predictor(1) >= 0.0_dp .and. &
                            predictor(1) <= 1.0_dp .and. &
                            predictor(5) >= 0.0_dp .and. &
                            predictor(5) <= 1.0_dp .and. &
                            abs(si(i)%z(2) - z_previous(2)) < 0.5_dp*pi
                    end if
                    z_previous = si(i)%z
                    if (use_predictor) then
                        call gpu_timestep_midpoint(si(i), f(i), warning_mode, &
                            ierr, step_nfev, predictor)
                    else
                        call gpu_timestep_midpoint(si(i), f(i), warning_mode, &
                            ierr, step_nfev)
                    end if
                    nfev(i) = nfev(i) + step_nfev
                    elapsed = elapsed + si(i)%dt
                    if (ierr /= 0) exit
                    use_predictor = .true.
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
