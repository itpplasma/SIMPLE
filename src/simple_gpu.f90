module simple_gpu
    ! OpenACC GPU offload of the default orbit-tracing hot path:
    ! Boozer field (isw_field_type=2) + explicit-implicit symplectic Euler
    ! (integmode=EXPL_IMPL_EULER). One particle per GPU thread.
    !
    ! The CPU integrator dispatches the field evaluation through a procedure
    ! pointer, which cannot be called inside an OpenACC compute region. This
    ! module keeps the device-side field evaluation local and reuses the shared
    ! symplectic-Euler algebra from orbit_symplectic_euler1. Equivalence is
    ! checked against the CPU path by test_gpu_orbit_bench.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util, only: pi
    use field_can_mod, only: field_can_t, get_derivatives, get_derivatives2
    use field_can_boozer, only: eval_field_booz
    use orbit_symplectic_base, only: symplectic_integrator_t, &
        SYMPLECTIC_STEP_OK, SYMPLECTIC_STEP_OUTSIDE_DOMAIN, &
        SYMPLECTIC_STEP_MAXITER, SYMPLECTIC_STEP_LINEAR_SOLVE, &
        SYMPLECTIC_STEP_BOUNDARY_LIMITED
    use orbit_symplectic_euler1, only: sympl_euler1_newton_iter
    use orbit_symplectic_euler1, only: sympl_euler1_extrapolate_field, sympl_euler1_advance_angles
    use boozer_sub, only: boozer_state
    use omp_lib, only: omp_get_thread_num
#ifdef _OPENACC
    use openacc, only: acc_get_num_devices, acc_set_device_num, acc_device_nvidia
#endif

    implicit none
    private

    integer, parameter :: maxit = 32

    public :: trace_orbits_gpu, trace_orbits_gpu_range

contains

    subroutine gpu_newton1(si, f, x, xlast, status)
        !$acc routine seq
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(inout) :: x(2)
        real(dp), intent(out) :: xlast(2)
        integer, intent(out) :: status

        real(dp) :: tolref(2)
        integer :: kit
        logical :: converged, linear_failed, boundary_limited
        logical :: step_boundary_limited

        status = SYMPLECTIC_STEP_MAXITER
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

            call eval_field_booz(f, x(1), si%z(2), si%z(3), 2)
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
        end if
        ! Non-convergence diagnostics (CPU writes fort.6601) are omitted on device.
    end subroutine gpu_newton1

    subroutine gpu_timestep_euler(si, f, ierr)
        !$acc routine seq
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        real(dp) :: x(2), xlast(2)
        integer :: ktau, newton_status
        logical :: crossed
        type(symplectic_integrator_t) :: accepted_integrator
        type(field_can_t) :: accepted_field

        ierr = 0
        ktau = 0
        do while (ktau < si%ntau)
            accepted_integrator = si
            accepted_field = f
            si%pthold = f%pth

            x(1) = si%z(1)
            x(2) = si%z(4)

            call gpu_newton1(si, f, x, xlast, newton_status)
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
                call eval_field_booz(f, si%z(1), si%z(2), si%z(3), 0)
                call get_derivatives(f, si%z(4))
            else
                call sympl_euler1_extrapolate_field(si, f, x, xlast)
            end if
            call sympl_euler1_advance_angles(si, f)

            ktau = ktau + 1
        end do
    end subroutine gpu_timestep_euler

    subroutine trace_orbits_gpu_range(si, f, npart, istart, iend, ntimstep, &
                                      ntau_macro, loss_step, zend)
        ! Evolve particles istart:iend on the currently selected OpenACC device.
        ! loss_step(i) = first timestep at which the orbit left the plasma (r>1),
        ! or ntimstep if confined for the whole trace. zend(:,i) is the final
        ! integrator state z = (r, th, ph, pphi).
        type(symplectic_integrator_t), intent(inout) :: si(npart)
        type(field_can_t), intent(inout) :: f(npart)
        integer, intent(in) :: npart, istart, iend, ntimstep
        integer, intent(in) :: ntau_macro(ntimstep)
        integer, intent(out) :: loss_step(npart)
        real(dp), intent(out) :: zend(4, npart)

        integer :: i, it, ktau, ierr, lstep

        !$acc parallel loop gang vector default(present) &
        !$acc&   copy(si(istart:iend), f(istart:iend)) copyin(ntau_macro) &
        !$acc&   copyout(loss_step(istart:iend), zend(:, istart:iend)) &
        !$acc&   private(it, ktau, ierr, lstep)
        do i = istart, iend
            ierr = 0
            lstep = ntimstep
            macro: do it = 2, ntimstep
                do ktau = 1, ntau_macro(it)
                    call gpu_timestep_euler(si(i), f(i), ierr)
                    if (ierr /= 0) exit
                end do
                if (ierr /= 0) then
                    lstep = it
                    exit macro
                end if
            end do macro
            loss_step(i) = lstep
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
#ifdef _OPENACC
        navail = acc_get_num_devices(acc_device_nvidia)
        if (navail < 1) navail = 1
        ngpu = min(req, navail)
#endif
        if (ngpu <= 1) then
            call trace_orbits_gpu_range(si, f, npart, 1, npart, ntimstep, &
                                        ntau_macro, loss_step, zend)
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
                                        ntau_macro, loss_step, zend)
        end if
        !$omp end parallel
    end subroutine trace_orbits_gpu

end module simple_gpu
