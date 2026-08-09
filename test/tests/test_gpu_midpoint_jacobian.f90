program test_gpu_midpoint_jacobian
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_chartmap, only: load_boozer_from_chartmap
    use field_can_mod, only: field_can_t, field_can_init, get_val
    use field_can_boozer, only: eval_field_booz
    use orbit_symplectic_base, only: symplectic_integrator_t
    use simple_gpu, only: gpu_midpoint_system, gpu_midpoint_newton_update, &
        gpu_midpoint_residual_decreased

    implicit none

    character(len=2048) :: filename
    integer :: status, column
    type(field_can_t) :: field
    type(symplectic_integrator_t) :: integrator
    real(dp) :: z(4), x0(5), h(5), residual(5), residual_plus(5)
    real(dp) :: residual_minus(5), jacobian(5, 5), jacobian_fd(5, 5)
    real(dp) :: max_error, scale, pphi

    call get_command_argument(1, filename, status=status)
    if (status /= 0 .or. len_trim(filename) == 0) &
        error stop 'usage: test_gpu_midpoint_jacobian.x chartmap.nc'

    call load_boozer_from_chartmap(trim(filename))
    call field_can_init(field, mu=0.01_dp, ro0=1.0e8_dp, vpar=0.2_dp)
    call eval_field_booz(field, 0.35_dp, 0.7_dp, 0.4_dp, 0)
    pphi = field%vpar*field%hph + field%Aph/field%ro0
    z = [0.35_dp, 0.7_dp, 0.4_dp, pphi]
    integrator%atol = 1.0e-15_dp
    integrator%rtol = 1.0e-12_dp
    integrator%ntau = 1
    integrator%dt = 1.0e-3_dp
    integrator%z = z
    call get_val(field, z(4))
    integrator%pthold = field%pth

    x0(1:4) = integrator%z + [2.0e-4_dp, -1.0e-3_dp, 1.5e-3_dp, 2.0e-5_dp]
    x0(5) = integrator%z(1) + 1.0e-4_dp
    h = [1.0e-3_dp, 1.0e-3_dp, 1.0e-3_dp, 1.0e-3_dp, 1.0e-3_dp]
    h(4) = 1.0e-3_dp*max(1.0_dp, abs(x0(4)))

    call gpu_midpoint_system(integrator, field, x0, residual, jacobian)
    do column = 1, 5
        x0(column) = x0(column) + 0.5_dp*h(column)
        call gpu_midpoint_system(integrator, field, x0, residual_plus, jacobian_fd)
        x0(column) = x0(column) - h(column)
        call gpu_midpoint_system(integrator, field, x0, residual_minus, jacobian_fd)
        x0(column) = x0(column) + 0.5_dp*h(column)
        jacobian_fd(:, column) = (residual_plus - residual_minus)/h(column)
    end do

    max_error = 0.0_dp
    do column = 1, 5
        scale = max(1.0_dp, maxval(abs(jacobian(:, column))), &
            maxval(abs(jacobian_fd(:, column))))
        max_error = max(max_error, maxval(abs(jacobian(:, column) - &
            jacobian_fd(:, column)))/scale)
    end do
    write (*, '(a,1pE12.4)') 'GPU midpoint analytic/finite-difference Jacobian error: ', &
        max_error
    ! Release builds of gfortran/NVHPC reorder the spline arithmetic enough
    ! to move this finite-difference oracle from 6e-5 to about 3e-4. Keep
    ! the threshold below the scale of an actual algebraic mismatch while
    ! accepting that implementation-dependent roundoff envelope.
    if (max_error > 5.0e-4_dp) &
        error stop 'GPU midpoint Jacobian oracle failed'
    call test_scaled_newton_update()
    call test_residual_decrease_rule()
    print *, 'test_gpu_midpoint_jacobian: PASSED'

contains

    subroutine test_scaled_newton_update()
        real(dp) :: jacobian(5, 5), rhs(5), variable_scale(5), correction(5)
        real(dp) :: residual_error
        integer :: info

        ! J*diag(variable_scale) is a moderate upper-triangular matrix, while
        ! J itself spans sixteen orders of magnitude. Check the returned
        ! correction against the original equation independently of the
        ! implementation's scaled solve.
        jacobian = 0.0_dp
        variable_scale = [1.0e-8_dp, 1.0e8_dp, 1.0e-4_dp, 1.0e4_dp, 1.0_dp]
        jacobian(1, 1) = 2.0e8_dp
        jacobian(1, 2) = 0.1e-8_dp
        jacobian(2, 2) = 1.5e-8_dp
        jacobian(1, 3) = -0.2e4_dp
        jacobian(2, 3) = 0.2e4_dp
        jacobian(3, 3) = 0.75e4_dp
        jacobian(1, 4) = 0.3e-4_dp
        jacobian(2, 4) = -0.1e-4_dp
        jacobian(3, 4) = 0.15e-4_dp
        jacobian(4, 4) = 2.25e-4_dp
        jacobian(1, 5) = -0.1_dp
        jacobian(2, 5) = 0.05_dp
        jacobian(3, 5) = -0.05_dp
        jacobian(4, 5) = 0.2_dp
        jacobian(5, 5) = 1.25_dp
        rhs = [1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp, -1.0_dp]

        call gpu_midpoint_newton_update(jacobian, rhs, variable_scale, .true., &
            .true., correction, info)
        if (info /= 0) error stop 'scaled midpoint Newton solve failed'
        residual_error = maxval(abs(matmul(jacobian, correction) - rhs) / &
            max(1.0_dp, abs(rhs)))
        if (residual_error > 1.0e-9_dp) &
            error stop 'scaled midpoint Newton solve residual failed'
    end subroutine test_scaled_newton_update

    subroutine test_residual_decrease_rule()
        real(dp) :: current_residual(5), trial_residual(5), tolref(5)

        current_residual = [1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp, -1.0_dp]
        tolref = [1.0_dp, 2.0_dp, 1.0_dp, 4.0_dp, 1.0_dp]
        trial_residual = 0.5_dp*current_residual
        if (.not. gpu_midpoint_residual_decreased(current_residual, &
                trial_residual, tolref, 1.0e-12_dp, 1.0e-10_dp)) &
            error stop 'midpoint residual decrease rule rejected a decrease'
        trial_residual(1) = 2.0_dp*current_residual(1)
        if (gpu_midpoint_residual_decreased(current_residual, trial_residual, &
                tolref, 1.0e-12_dp, 1.0e-10_dp)) &
            error stop 'midpoint residual decrease rule accepted an increase'
    end subroutine test_residual_decrease_rule
end program test_gpu_midpoint_jacobian
