module simple_gpu
    ! OpenACC GPU offload of the default orbit-tracing hot path:
    ! Boozer field (isw_field_type=2) + explicit-implicit symplectic Euler
    ! (integmode=EXPL_IMPL_EULER). One particle per GPU thread.
    !
    ! The CPU integrator (orbit_symplectic) dispatches the field evaluation
    ! through a procedure pointer, which cannot be called inside an OpenACC
    ! compute region. This module provides device-resident (acc routine seq)
    ! copies of the explicit-implicit Euler step and its 2x2 Newton solver that
    ! call eval_field_booz directly. The arithmetic mirrors
    ! orbit_timestep_sympl_expl_impl_euler / newton1 / f_sympl_euler1 /
    ! jac_sympl_euler1 exactly (extrap_field = .true. branch); equivalence is
    ! checked against the CPU path by test_gpu_orbit.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_mod, only: field_can_t, get_val, get_derivatives, get_derivatives2
    use field_can_boozer, only: eval_field_booz
    use orbit_symplectic_base, only: symplectic_integrator_t
    use boozer_sub, only: torflux_gpu
    use omp_lib, only: omp_get_thread_num
#ifdef _OPENACC
    use openacc, only: acc_get_num_devices, acc_set_device_num, acc_device_nvidia
#endif

    implicit none
    private

    integer, parameter :: maxit = 32

    public :: trace_orbits_gpu, trace_orbits_gpu_range

contains

    subroutine gpu_f_sympl_euler1(si, f, x, fvec)
        !$acc routine seq
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: fvec(2)

        call eval_field_booz(f, x(1), si%z(2), si%z(3), 2)
        call get_derivatives2(f, x(2))

        fvec(1) = f%dpth(1)*(f%pth - si%pthold) &
                  + si%dt*(f%dH(2)*f%dpth(1) - f%dH(1)*f%dpth(2))
        fvec(2) = f%dpth(1)*(x(2) - si%z(4)) &
                  + si%dt*(f%dH(3)*f%dpth(1) - f%dH(1)*f%dpth(3))
    end subroutine gpu_f_sympl_euler1

    subroutine gpu_jac_sympl_euler1(si, f, x, jac)
        !$acc routine seq
        type(symplectic_integrator_t), intent(in) :: si
        type(field_can_t), intent(in) :: f
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: jac(2, 2)

        jac(1, 1) = f%d2pth(1)*(f%pth - si%pthold) + f%dpth(1)**2 &
            + si%dt*(f%d2H(2)*f%dpth(1) + f%dH(2)*f%d2pth(1) &
                     - f%d2H(1)*f%dpth(2) - f%dH(1)*f%d2pth(2))
        jac(1, 2) = f%d2pth(7)*(f%pth - si%pthold) + f%dpth(1)*f%dpth(4) &
            + si%dt*(f%d2H(8)*f%dpth(1) + f%dH(2)*f%d2pth(7) &
                     - f%d2H(7)*f%dpth(2) - f%dH(1)*f%d2pth(8))
        jac(2, 1) = f%d2pth(1)*(x(2) - si%z(4)) &
            + si%dt*(f%d2H(3)*f%dpth(1) + f%dH(3)*f%d2pth(1) &
                     - f%d2H(1)*f%dpth(3) - f%dH(1)*f%d2pth(3))
        jac(2, 2) = f%d2pth(7)*(x(2) - si%z(4)) + f%dpth(1) &
            + si%dt*(f%d2H(9)*f%dpth(1) + f%dH(3)*f%d2pth(7) &
                     - f%d2H(7)*f%dpth(3) - f%dH(1)*f%d2pth(9))
    end subroutine gpu_jac_sympl_euler1

    subroutine gpu_newton1(si, f, x, xlast)
        !$acc routine seq
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(inout) :: x(2)
        real(dp), intent(out) :: xlast(2)

        real(dp) :: fvec(2), fjac(2, 2), ijac(2, 2)
        real(dp) :: tolref(2)
        integer :: kit

        tolref(1) = 1d0
        tolref(2) = dabs(1d1*torflux_gpu/f%ro0)

        do kit = 1, maxit
            if (x(1) > 1d0) return
            if (x(1) < 0d0) x(1) = 0.01d0

            call gpu_f_sympl_euler1(si, f, x, fvec)
            call gpu_jac_sympl_euler1(si, f, x, fjac)

            ijac(1, 1) = 1d0/(fjac(1, 1) - fjac(1, 2)*fjac(2, 1)/fjac(2, 2))
            ijac(1, 2) = -1d0/(fjac(1, 1)*fjac(2, 2)/fjac(1, 2) - fjac(2, 1))
            ijac(2, 1) = -1d0/(fjac(1, 1)*fjac(2, 2)/fjac(2, 1) - fjac(1, 2))
            ijac(2, 2) = 1d0/(fjac(2, 2) - fjac(1, 2)*fjac(2, 1)/fjac(1, 1))

            xlast = x
            x(1) = x(1) - (ijac(1, 1)*fvec(1) + ijac(1, 2)*fvec(2))
            x(2) = x(2) - (ijac(2, 1)*fvec(1) + ijac(2, 2)*fvec(2))

            tolref(2) = max(dabs(x(2)), tolref(2))

            if (dabs(fvec(1)) < si%atol .and. dabs(fvec(2)) < si%atol) return
            if (dabs(x(1) - xlast(1)) < si%rtol*tolref(1) .and. &
                dabs(x(2) - xlast(2)) < si%rtol*tolref(2)) return
        end do
        ! Non-convergence diagnostics (CPU writes fort.6601) are omitted on device.
    end subroutine gpu_newton1

    subroutine gpu_timestep_euler(si, f, ierr)
        !$acc routine seq
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr

        real(dp) :: x(2), xlast(2)
        integer :: ktau

        ierr = 0
        ktau = 0
        do while (ktau < si%ntau)
            si%pthold = f%pth

            x(1) = si%z(1)
            x(2) = si%z(4)

            call gpu_newton1(si, f, x, xlast)

            if (x(1) > 1.0d0) then
                ierr = 1
                return
            end if
            if (x(1) < 0.0d0) x(1) = 0.01d0

            si%z(1) = x(1)
            si%z(4) = x(2)

            ! extrap_field = .true. branch (see orbit_symplectic_base)
            f%pth = f%pth + f%dpth(1)*(x(1) - xlast(1)) + f%dpth(4)*(x(2) - xlast(2))
            f%dH(1) = f%dH(1) + f%d2H(1)*(x(1) - xlast(1)) + f%d2H(7)*(x(2) - xlast(2))
            f%dpth(1) = f%dpth(1) + f%d2pth(1)*(x(1) - xlast(1)) &
                        + f%d2pth(7)*(x(2) - xlast(2))
            f%vpar = f%vpar + f%dvpar(1)*(x(1) - xlast(1)) + f%dvpar(4)*(x(2) - xlast(2))
            f%hth = f%hth + f%dhth(1)*(x(1) - xlast(1))
            f%hph = f%hph + f%dhph(1)*(x(1) - xlast(1))

            si%z(2) = si%z(2) + si%dt*f%dH(1)/f%dpth(1)
            si%z(3) = si%z(3) + si%dt*(f%vpar - f%dH(1)/f%dpth(1)*f%hth)/f%hph

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
