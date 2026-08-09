program test_midpoint_jacobian
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_can_mod, only: field_can_t, field_can_from_name, field_can_init, &
        eval_field => evaluate, get_derivatives
    use orbit_symplectic, only: f_midpoint_part1, f_midpoint_part2, &
        jac_midpoint_part1, jac_midpoint_part2
    use orbit_symplectic_base, only: symplectic_integrator_t

    implicit none

    type(field_can_t) :: field
    type(symplectic_integrator_t) :: integrator
    real(dp) :: x(5), residual_plus(5), residual_minus(5), residual(5)
    real(dp) :: analytic(5, 5), finite_difference(5, 5), h
    real(dp) :: scale, relative_error, maximum_error
    integer :: column, point, row

    call field_can_from_name('test')
    call field_can_init(field, mu=0.03_dp, ro0=1.0_dp, vpar=0.1_dp)

    integrator%z = [0.2_dp, 0.7_dp, 1.1_dp, 0.15_dp]
    integrator%dt = 1.0e-3_dp
    integrator%atol = 1.0e-15_dp
    integrator%rtol = 1.0e-10_dp
    call eval_field(field, integrator%z(1), integrator%z(2), &
        integrator%z(3), 0)
    call get_derivatives(field, integrator%z(4))
    integrator%pthold = field%pth

    maximum_error = 0.0_dp
    do point = 1, 4
        x = [0.20_dp + 0.01_dp*real(point, dp), &
            0.70_dp + 0.03_dp*real(point, dp), &
            1.10_dp - 0.02_dp*real(point, dp), &
            0.15_dp + 0.005_dp*real(point, dp), &
            0.205_dp + 0.008_dp*real(point, dp)]

        call midpoint_residual_and_jacobian(integrator, field, x, residual, analytic)

        do column = 1, size(x)
            h = sqrt(epsilon(1.0_dp))*max(1.0_dp, abs(x(column)))
            x(column) = x(column) + h
            call midpoint_residual(integrator, field, x, residual_plus)
            x(column) = x(column) - 2.0_dp*h
            call midpoint_residual(integrator, field, x, residual_minus)
            x(column) = x(column) + h
            finite_difference(:, column) = (residual_plus - residual_minus)/(2.0_dp*h)
        end do

        do column = 1, size(x)
            do row = 1, size(x)
                scale = max(1.0_dp, abs(analytic(row, column)), &
                    abs(finite_difference(row, column)))
                relative_error = abs(analytic(row, column) - &
                    finite_difference(row, column))/scale
                maximum_error = max(maximum_error, relative_error)
            end do
        end do
    end do

    if (maximum_error > 5.0e-6_dp) then
        write(*, '(a, 1pE12.4)') 'midpoint Jacobian finite-difference error: ', &
            maximum_error
        error stop 'midpoint analytic Jacobian failed independent oracle'
    end if

    write(*, '(a, 1pE12.4)') 'midpoint Jacobian oracle passed; max error: ', &
        maximum_error

contains

    subroutine midpoint_residual(integrator, field, x, residual)
        type(symplectic_integrator_t), intent(in) :: integrator
        type(field_can_t), intent(inout) :: field
        real(dp), intent(in) :: x(5)
        real(dp), intent(out) :: residual(5)
        type(symplectic_integrator_t) :: local_integrator

        local_integrator = integrator
        call f_midpoint_part1(local_integrator, field, 5, x, residual)
        call f_midpoint_part2(local_integrator, field, 5, x, residual)
    end subroutine midpoint_residual

    subroutine midpoint_residual_and_jacobian(integrator, field, x, residual, jacobian)
        type(symplectic_integrator_t), intent(in) :: integrator
        type(field_can_t), intent(inout) :: field
        real(dp), intent(in) :: x(5)
        real(dp), intent(out) :: residual(5), jacobian(5, 5)
        type(symplectic_integrator_t) :: local_integrator
        type(field_can_t) :: midpoint_field

        local_integrator = integrator
        call f_midpoint_part1(local_integrator, field, 5, x, residual)
        call jac_midpoint_part1(local_integrator, field, x, jacobian)
        midpoint_field = field
        call f_midpoint_part2(local_integrator, field, 5, x, residual)
        call jac_midpoint_part2(local_integrator, field, midpoint_field, x, jacobian)
    end subroutine midpoint_residual_and_jacobian

end program test_midpoint_jacobian
