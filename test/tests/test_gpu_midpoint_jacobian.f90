program test_gpu_midpoint_jacobian
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_chartmap, only: load_boozer_from_chartmap
    use field_can_mod, only: field_can_t, field_can_init, get_val
    use field_can_boozer, only: eval_field_booz
    use orbit_symplectic_base, only: symplectic_integrator_t
    use simple_gpu, only: gpu_midpoint_system

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
    if (max_error > 2.0e-4_dp) &
        error stop 'GPU midpoint Jacobian oracle failed'
    print *, 'test_gpu_midpoint_jacobian: PASSED'
end program test_gpu_midpoint_jacobian
