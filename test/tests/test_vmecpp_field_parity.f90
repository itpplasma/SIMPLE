program test_vmecpp_field_parity
    use field_vmecpp, only: close_vmecpp_geometry, evaluate_vmecpp_geometry, &
        magfie_vmecpp, open_vmecpp_geometry, vmecpp_geometry_point_t
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: value_tolerance = 1.0e-12_dp
    real(dp), parameter :: derivative_tolerance = 1.0e-9_dp
    character(1024) :: input_file
    real(dp) :: x(3), actual(14), expected(14)
    integer :: i

    call get_command_argument(1, input_file)
    if (len_trim(input_file) == 0) error stop 'expected a VMEC++ input path'

    call open_vmecpp_geometry(trim(input_file))
    do i = 1, 3
        x = [0.27_dp+0.17_dp*i, 0.31_dp*i, 0.17_dp*i]
        call reference_field(x, expected)
        call evaluate(magfie_vmecpp, x, actual)
        call unscale_field(actual)
        call assert_close(actual(1:2), expected(1:2), value_tolerance, &
            'field values')
        call assert_close(actual(6:11), expected(6:11), value_tolerance, &
            'normalized covariant and contravariant fields')
        call assert_close(actual(3:5), expected(3:5), derivative_tolerance, &
            'magnetic-field derivatives')
        call assert_close(actual(12:14), expected(12:14), &
            derivative_tolerance, 'covariant curl')
    end do
    call close_vmecpp_geometry

contains

    subroutine evaluate(field_function, position, values)
        interface
            subroutine field_function(x, bmod, sqrtg, bder, hcovar, hctrvr, &
                    hcurl)
                import :: dp
                real(dp), intent(in) :: x(3)
                real(dp), intent(out) :: bmod, sqrtg
                real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), &
                    hcurl(3)
            end subroutine field_function
        end interface
        real(dp), intent(in) :: position(3)
        real(dp), intent(out) :: values(14)

        call field_function(position, values(1), values(2), values(3:5), &
            values(6:8), values(9:11), values(12:14))
    end subroutine evaluate

    subroutine reference_field(position, values)
        real(dp), intent(in) :: position(3)
        real(dp), intent(out) :: values(14)
        real(dp) :: base(8), plus(8), minus(8), plus2(8), minus2(8)
        real(dp) :: derivative_bmod(3), derivative_hcovar(3, 3)
        real(dp) :: step, cross_value(3)
        integer :: coordinate, component

        call reference_components(position, base)
        values = 0.0_dp
        values(1) = base(1)
        values(2) = base(2)
        values(6:8) = base(3:5)
        values(9:11) = base(6:8)

        do coordinate = 1, 3
            step = 1.0e-4_dp
            call shifted_components(position, coordinate, 2.0_dp*step, plus2)
            call shifted_components(position, coordinate, step, plus)
            call shifted_components(position, coordinate, -step, minus)
            call shifted_components(position, coordinate, -2.0_dp*step, &
                minus2)
            derivative_bmod(coordinate) = (-plus2(1)+8.0_dp*plus(1)- &
                8.0_dp*minus(1)+minus2(1))/(12.0_dp*step)
            do component = 1, 3
                derivative_hcovar(component, coordinate) = &
                    (-plus2(component+2)+8.0_dp*plus(component+2)- &
                    8.0_dp*minus(component+2)+minus2(component+2))/ &
                    (12.0_dp*step)
            end do
        end do

        values(3:5) = derivative_bmod/base(1)
        cross_value = [derivative_hcovar(3, 2)-derivative_hcovar(2, 3), &
            derivative_hcovar(1, 3)-derivative_hcovar(3, 1), &
            derivative_hcovar(2, 1)-derivative_hcovar(1, 2)]
        values(12:14) = cross_value/base(2)
    end subroutine reference_field

    subroutine shifted_components(position, coordinate, shift, values)
        real(dp), intent(in) :: position(3), shift
        integer, intent(in) :: coordinate
        real(dp), intent(out) :: values(8)
        real(dp) :: shifted(3)

        shifted = position
        shifted(coordinate) = shifted(coordinate) + shift
        call reference_components(shifted, values)
    end subroutine shifted_components

    subroutine reference_components(position, values)
        real(dp), intent(in) :: position(3)
        real(dp), intent(out) :: values(8)
        type(vmecpp_geometry_point_t) :: point
        real(dp) :: basis(3, 3), metric(3, 3), bcontrav(3)
        real(dp) :: bcovar(3), bvector(3), cosine, sine, sqrtg
        real(dp) :: poloidal_flux, toroidal_flux, inverse_sqrtg

        call evaluate_vmecpp_geometry(position(1), position(2), position(3), &
            point)
        cosine = cos(position(3))
        sine = sin(position(3))
        basis(:, 1) = [point%r(2)*cosine, point%r(2)*sine, point%z(2)]
        basis(:, 2) = [point%r(3)*cosine, point%r(3)*sine, point%z(3)]
        basis(:, 3) = [point%r(4)*cosine-point%r(1)*sine, &
            point%r(4)*sine+point%r(1)*cosine, point%z(4)]
        sqrtg = dot_product(basis(:, 1), &
            cross_product(basis(:, 2), basis(:, 3)))

        metric = 0.0_dp
        metric(1, 1) = dot_product(basis(:, 1), basis(:, 1))
        metric(1, 2) = dot_product(basis(:, 1), basis(:, 2))
        metric(1, 3) = dot_product(basis(:, 1), basis(:, 3))
        metric(2, 1) = metric(1, 2)
        metric(2, 2) = dot_product(basis(:, 2), basis(:, 2))
        metric(2, 3) = dot_product(basis(:, 2), basis(:, 3))
        metric(3, 1) = metric(1, 3)
        metric(3, 2) = metric(2, 3)
        metric(3, 3) = dot_product(basis(:, 3), basis(:, 3))

        poloidal_flux = point%poloidal_flux(2)
        toroidal_flux = point%toroidal_flux(2)
        inverse_sqrtg = 1.0_dp/(2.0_dp*acos(-1.0_dp)*sqrtg)
        bcontrav = [0.0_dp, &
            -(poloidal_flux-toroidal_flux*point%lambda(4))*inverse_sqrtg, &
            -toroidal_flux*(1.0_dp+point%lambda(3))*inverse_sqrtg]
        bcovar = matmul(metric, bcontrav)
        bvector = matmul(basis, bcontrav)

        values(1) = sqrt(dot_product(bvector, bvector))
        values(2) = sqrtg
        values(3:5) = bcovar/values(1)
        values(6:8) = bcontrav/values(1)
    end subroutine reference_components

    subroutine unscale_field(values)
        real(dp), intent(inout) :: values(14)

        values(1) = values(1)/1.0e4_dp
        values(2) = values(2)/1.0e6_dp
        values(6:8) = values(6:8)/1.0e2_dp
        values(9:11) = values(9:11)*1.0e2_dp
        values(12:14) = values(12:14)*1.0e4_dp
    end subroutine unscale_field

    subroutine assert_close(actual, expected, tolerance, description)
        real(dp), intent(in) :: actual(:), expected(:), tolerance
        character(*), intent(in) :: description
        real(dp) :: error

        error = maxval(abs(actual-expected))
        if (error > tolerance) then
            write(*, '(A,ES12.4)') trim(description)//' absolute error: ', &
                error
            error stop 'VMEC++ direct field failed its analytic oracle'
        end if
    end subroutine assert_close

    pure function cross_product(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)

        c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), &
            a(1)*b(2)-a(2)*b(1)]
    end function cross_product

end program test_vmecpp_field_parity
