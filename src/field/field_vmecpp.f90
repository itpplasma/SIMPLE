module field_vmecpp
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, &
        c_int, c_null_char, c_null_ptr, c_ptr
    implicit none
    private

    integer, parameter :: vmecpp_geometry_jet_size = 10

    ! The entries are value, ds, dtheta, dzeta, dss, ds_dtheta,
    ! ds_dzeta, dtheta2, dtheta_dzeta, and dzeta2.
    type, bind(C), public :: vmecpp_geometry_point_t
        real(c_double) :: r(vmecpp_geometry_jet_size)
        real(c_double) :: z(vmecpp_geometry_jet_size)
        real(c_double) :: lambda(vmecpp_geometry_jet_size)
        real(c_double) :: toroidal_flux(vmecpp_geometry_jet_size)
        real(c_double) :: poloidal_flux(vmecpp_geometry_jet_size)
    end type vmecpp_geometry_point_t

    type, bind(C) :: vmecpp_geometry_metadata_t
        integer(c_int) :: nfp
        real(c_double) :: major_radius
    end type vmecpp_geometry_metadata_t

    type(c_ptr), save :: geometry_handle = c_null_ptr

    public :: close_vmecpp_geometry, evaluate_vmecpp_geometry, &
        get_vmecpp_metadata, magfie_vmecpp, open_vmecpp_geometry

    interface
        function vmecpp_geometry_create(input_path, output) bind(C)
            import :: c_char, c_int, c_ptr
            character(c_char), intent(in) :: input_path(*)
            type(c_ptr), intent(out) :: output
            integer(c_int) :: vmecpp_geometry_create
        end function vmecpp_geometry_create

        subroutine vmecpp_geometry_destroy(handle) bind(C)
            import :: c_ptr
            type(c_ptr), value :: handle
        end subroutine vmecpp_geometry_destroy

        function vmecpp_geometry_get_metadata(handle, output) bind(C)
            import :: c_int, c_ptr, vmecpp_geometry_metadata_t
            type(c_ptr), value :: handle
            type(vmecpp_geometry_metadata_t), intent(out) :: output
            integer(c_int) :: vmecpp_geometry_get_metadata
        end function vmecpp_geometry_get_metadata

        function vmecpp_geometry_evaluate(handle, s, theta, zeta, output) &
                bind(C)
            import :: c_double, c_int, c_ptr, vmecpp_geometry_point_t
            type(c_ptr), value :: handle
            real(c_double), value :: s, theta, zeta
            type(vmecpp_geometry_point_t), intent(out) :: output
            integer(c_int) :: vmecpp_geometry_evaluate
        end function vmecpp_geometry_evaluate
    end interface

contains

    subroutine open_vmecpp_geometry(input_path)
        character(*), intent(in) :: input_path
        character(c_char), allocatable :: c_path(:)
        integer(c_int) :: status
        integer :: i, length

        call close_vmecpp_geometry
        length = len_trim(input_path)
        allocate(c_path(length + 1))
        do i = 1, length
            c_path(i) = input_path(i:i)
        end do
        c_path(length + 1) = c_null_char
        status = vmecpp_geometry_create(c_path, geometry_handle)
        if (status /= 0 .or. .not. c_associated(geometry_handle)) &
            error stop 'VMEC++ could not create geometry from the input file'
    end subroutine open_vmecpp_geometry

    subroutine close_vmecpp_geometry
        if (c_associated(geometry_handle)) then
            call vmecpp_geometry_destroy(geometry_handle)
            geometry_handle = c_null_ptr
        end if
    end subroutine close_vmecpp_geometry

    subroutine get_vmecpp_metadata(nfp, major_radius)
        integer, intent(out) :: nfp
        real(c_double), intent(out) :: major_radius
        type(vmecpp_geometry_metadata_t) :: metadata
        integer(c_int) :: status

        status = vmecpp_geometry_get_metadata(geometry_handle, metadata)
        if (status /= 0) error stop 'VMEC++ geometry metadata query failed'
        nfp = metadata%nfp
        major_radius = metadata%major_radius
    end subroutine get_vmecpp_metadata

    subroutine evaluate_vmecpp_geometry(s, theta, zeta, point)
        real(c_double), intent(in) :: s, theta, zeta
        type(vmecpp_geometry_point_t), intent(out) :: point
        integer(c_int) :: status

        if (.not. c_associated(geometry_handle)) &
            error stop 'VMEC++ geometry is not initialized'
        status = vmecpp_geometry_evaluate( &
            geometry_handle, max(0.0_c_double, min(1.0_c_double, s)), &
            theta, zeta, point)
        if (status /= 0) error stop 'VMEC++ geometry evaluation failed'
    end subroutine evaluate_vmecpp_geometry

    subroutine magfie_vmecpp(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        real(c_double), intent(in) :: x(3)
        real(c_double), intent(out) :: bmod, sqrtg
        real(c_double), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)
        type(vmecpp_geometry_point_t) :: point
        real(c_double) :: basis(3, 3), basis_derivative(3, 3, 3)
        real(c_double) :: metric(3, 3), metric_derivative(3, 3, 3)
        real(c_double) :: bcontrav(3), bcontrav_derivative(3, 3)
        real(c_double) :: bcovar(3), bcovar_derivative(3, 3)
        real(c_double) :: bvector(3), bvector_derivative(3, 3)
        real(c_double) :: sqrtg_derivative(3), hcovar_derivative(3, 3)
        real(c_double) :: d_poloidal_flux(3), d_toroidal_flux(3)
        real(c_double) :: d_lambda_theta(3), d_lambda_zeta(3)
        real(c_double) :: cosine, sine, inverse_sqrtg, derivative_inverse
        real(c_double) :: poloidal_flux, toroidal_flux, lambda_theta, lambda_zeta
        real(c_double) :: numerator, numerator_derivative
        real(c_double) :: cross_value(3)
        integer :: coordinate, i, j

        call evaluate_vmecpp_geometry(x(1), x(2), x(3), point)
        cosine = cos(x(3))
        sine = sin(x(3))

        basis(:, 1) = [point%r(2)*cosine, point%r(2)*sine, point%z(2)]
        basis(:, 2) = [point%r(3)*cosine, point%r(3)*sine, point%z(3)]
        basis(:, 3) = [point%r(4)*cosine-point%r(1)*sine, &
            point%r(4)*sine+point%r(1)*cosine, point%z(4)]

        basis_derivative(:, 1, 1) = [point%r(5)*cosine, point%r(5)*sine, &
            point%z(5)]
        basis_derivative(:, 2, 1) = [point%r(6)*cosine, point%r(6)*sine, &
            point%z(6)]
        basis_derivative(:, 3, 1) = [point%r(7)*cosine-point%r(2)*sine, &
            point%r(7)*sine+point%r(2)*cosine, point%z(7)]

        basis_derivative(:, 1, 2) = [point%r(6)*cosine, point%r(6)*sine, &
            point%z(6)]
        basis_derivative(:, 2, 2) = [point%r(8)*cosine, point%r(8)*sine, &
            point%z(8)]
        basis_derivative(:, 3, 2) = [point%r(9)*cosine-point%r(3)*sine, &
            point%r(9)*sine+point%r(3)*cosine, point%z(9)]

        basis_derivative(:, 1, 3) = [point%r(7)*cosine-point%r(2)*sine, &
            point%r(7)*sine+point%r(2)*cosine, point%z(7)]
        basis_derivative(:, 2, 3) = [point%r(9)*cosine-point%r(3)*sine, &
            point%r(9)*sine+point%r(3)*cosine, point%z(9)]
        basis_derivative(:, 3, 3) = [point%r(10)*cosine- &
            2.0_c_double*point%r(4)*sine-point%r(1)*cosine, &
            point%r(10)*sine+2.0_c_double*point%r(4)*cosine- &
            point%r(1)*sine, point%z(10)]

        sqrtg = dot_product(basis(:, 1), &
            cross_product(basis(:, 2), basis(:, 3)))
        do coordinate = 1, 3
            sqrtg_derivative(coordinate) = dot_product( &
                basis_derivative(:, 1, coordinate), &
                cross_product(basis(:, 2), basis(:, 3))) + &
                dot_product(basis(:, 1), cross_product( &
                basis_derivative(:, 2, coordinate), basis(:, 3))) + &
                dot_product(basis(:, 1), cross_product(basis(:, 2), &
                basis_derivative(:, 3, coordinate)))
        end do

        do i = 1, 3
            do j = 1, 3
                metric(i, j) = dot_product(basis(:, i), basis(:, j))
                do coordinate = 1, 3
                    metric_derivative(i, j, coordinate) = &
                        dot_product(basis_derivative(:, i, coordinate), &
                        basis(:, j)) + dot_product(basis(:, i), &
                        basis_derivative(:, j, coordinate))
                end do
            end do
        end do

        poloidal_flux = point%poloidal_flux(2)
        toroidal_flux = point%toroidal_flux(2)
        d_poloidal_flux = [point%poloidal_flux(5), 0.0_c_double, 0.0_c_double]
        d_toroidal_flux = [point%toroidal_flux(5), 0.0_c_double, 0.0_c_double]
        lambda_theta = point%lambda(3)
        lambda_zeta = point%lambda(4)
        d_lambda_theta = [point%lambda(6), point%lambda(8), &
            point%lambda(9)]
        d_lambda_zeta = [point%lambda(7), point%lambda(9), &
            point%lambda(10)]

        inverse_sqrtg = 1.0_c_double/(2.0_c_double*acos(-1.0_c_double)*sqrtg)
        bcontrav = [0.0_c_double, &
            -(poloidal_flux-toroidal_flux*lambda_zeta)*inverse_sqrtg, &
            -toroidal_flux*(1.0_c_double+lambda_theta)*inverse_sqrtg]
        do coordinate = 1, 3
            derivative_inverse = -inverse_sqrtg*sqrtg_derivative(coordinate)/sqrtg
            numerator = poloidal_flux-toroidal_flux*lambda_zeta
            numerator_derivative = d_poloidal_flux(coordinate)- &
                d_toroidal_flux(coordinate)*lambda_zeta- &
                toroidal_flux*d_lambda_zeta(coordinate)
            bcontrav_derivative(1, coordinate) = 0.0_c_double
            bcontrav_derivative(2, coordinate) = -(numerator_derivative* &
                inverse_sqrtg+numerator*derivative_inverse)
            numerator = toroidal_flux*(1.0_c_double+lambda_theta)
            numerator_derivative = d_toroidal_flux(coordinate)* &
                (1.0_c_double+lambda_theta)+toroidal_flux* &
                d_lambda_theta(coordinate)
            bcontrav_derivative(3, coordinate) = -(numerator_derivative* &
                inverse_sqrtg+numerator*derivative_inverse)
        end do

        bcovar = matmul(metric, bcontrav)
        do coordinate = 1, 3
            bcovar_derivative(:, coordinate) = 0.0_c_double
            do i = 1, 3
                do j = 1, 3
                    bcovar_derivative(i, coordinate) = &
                        bcovar_derivative(i, coordinate)+ &
                        metric_derivative(i, j, coordinate)*bcontrav(j)+ &
                        metric(i, j)*bcontrav_derivative(j, coordinate)
                end do
            end do
        end do

        bvector = matmul(basis, bcontrav)
        do coordinate = 1, 3
            bvector_derivative(:, coordinate) = &
                matmul(basis_derivative(:, :, coordinate), bcontrav)+ &
                matmul(basis, bcontrav_derivative(:, coordinate))
        end do
        bmod = sqrt(dot_product(bvector, bvector))
        do coordinate = 1, 3
            bder(coordinate) = dot_product(bvector, &
                bvector_derivative(:, coordinate))/(bmod*bmod)
        end do
        hcovar = bcovar/bmod
        do coordinate = 1, 3
            hcovar_derivative(:, coordinate) = &
                bcovar_derivative(:, coordinate)/bmod- &
                hcovar*bder(coordinate)
        end do
        hctrvr = bcontrav/bmod
        cross_value = [hcovar_derivative(3, 2)-hcovar_derivative(2, 3), &
            hcovar_derivative(1, 3)-hcovar_derivative(3, 1), &
            hcovar_derivative(2, 1)-hcovar_derivative(1, 2)]
        hcurl = cross_value/sqrtg

        bmod = bmod*1.0e4_c_double
        sqrtg = sqrtg*1.0e6_c_double
        hcovar = hcovar*1.0e2_c_double
        hctrvr = hctrvr/1.0e2_c_double
        hcurl = hcurl/1.0e4_c_double
    end subroutine magfie_vmecpp

    pure function cross_product(a, b) result(c)
        real(c_double), intent(in) :: a(3), b(3)
        real(c_double) :: c(3)

        c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), &
            a(1)*b(2)-a(2)*b(1)]
    end function cross_product

end module field_vmecpp
