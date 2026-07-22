program probe_boozer_chartmap
    !> Print the executable-native signed Boozer field and Cartesian geometry.
    !>
    !> Usage: probe_boozer_chartmap.x <chartmap.nc>
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_boozer_chartmap, only: boozer_chartmap_field_t, create_boozer_chartmap_field

    implicit none

    type(boozer_chartmap_field_t), allocatable :: field
    character(len=1024) :: filename
    real(dp), parameter :: points(3, 4) = reshape([ &
        0.32_dp, 0.31_dp, 0.17_dp, &
        0.43_dp, 1.23_dp, 0.47_dp, &
        0.68_dp, 2.40_dp, 1.13_dp, &
        0.91_dp, 5.20_dp, 2.71_dp], [3, 4])
    real(dp) :: u(3), x_cart(3), Acov(3), hcov(3), h_cart(3), B_cart(3), Bmod
    real(dp) :: e_cov(3, 3), g(3, 3), ginv(3, 3), sqrtg_abs, jacobian_oriented
    integer :: point

    if (command_argument_count() /= 1) then
        error stop 'usage: probe_boozer_chartmap.x <chartmap.nc>'
    end if
    call get_command_argument(1, filename)
    call create_boozer_chartmap_field(trim(filename), field)

    do point = 1, size(points, 2)
        u = points(:, point)
        call field%evaluate(u, Acov, hcov, Bmod)
        call field%coords%evaluate_cart(u, x_cart)
        call field%coords%covariant_basis(u, e_cov)
        call field%coords%metric_tensor(u, g, ginv, sqrtg_abs)
        call field%coords%cov_to_cart(u, hcov, h_cart)
        B_cart = Bmod*h_cart
        jacobian_oriented = determinant3(e_cov)

        write (*, '(a,i0)') 'SIMPLE_PROBE point=', point
        call print_vector('u_rho_thetaB_zetaB=', u)
        call print_vector('position_cartesian_cm=', x_cart)
        call print_vector('Acov_G_cm2=', Acov)
        call print_vector('hcov=', hcov)
        write (*, '(a,es24.16)') 'SIMPLE_PROBE Bmod_G=', Bmod
        call print_vector('B_cartesian_G=', B_cart)
        write (*, '(a,es24.16)') 'SIMPLE_PROBE sqrtg_abs_cm3=', sqrtg_abs
        write (*, '(a,es24.16)') 'SIMPLE_PROBE jacobian_oriented_cm3=', jacobian_oriented
    end do

contains

    subroutine print_vector(label, vector)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: vector(3)

        write (*, '(a,3(1x,es24.16))') 'SIMPLE_PROBE '//label, vector
    end subroutine print_vector

    pure function determinant3(matrix) result(determinant)
        real(dp), intent(in) :: matrix(3, 3)
        real(dp) :: determinant

        determinant = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)) &
            - matrix(1, 2)*(matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) &
            + matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1))
    end function determinant3

end program probe_boozer_chartmap
