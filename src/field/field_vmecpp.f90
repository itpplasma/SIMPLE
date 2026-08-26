module field_vmecpp
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, &
        c_int, c_null_char, c_null_ptr, c_ptr
    implicit none
    private

    type, bind(C), public :: vmecpp_geometry_point_t
        real(c_double) :: r(4)
        real(c_double) :: z(4)
        real(c_double) :: lambda(4)
        real(c_double) :: toroidal_flux(4)
        real(c_double) :: poloidal_flux(4)
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

    subroutine field_components(x, bmod, sqrtg, hcovar, hctrvr)
        real(c_double), intent(in) :: x(3)
        real(c_double), intent(out) :: bmod, sqrtg
        real(c_double), intent(out) :: hcovar(3), hctrvr(3)
        type(vmecpp_geometry_point_t) :: point
        real(c_double) :: e_s(3), e_theta(3), e_zeta(3), metric(3, 3)
        real(c_double) :: bcontrav(3), bcovar(3), bvector(3)
        real(c_double) :: cosine, sine, flux_scale
        call evaluate_vmecpp_geometry(x(1), x(2), x(3), point)
        cosine = cos(x(3))
        sine = sin(x(3))
        e_s = [point%r(2)*cosine, point%r(2)*sine, point%z(2)]
        e_theta = [point%r(3)*cosine, point%r(3)*sine, point%z(3)]
        e_zeta = [point%r(4)*cosine-point%r(1)*sine, &
            point%r(4)*sine+point%r(1)*cosine, point%z(4)]
        sqrtg = dot_product(e_s, cross_product(e_theta, e_zeta))
        flux_scale = 1.0_c_double/(2.0_c_double*acos(-1.0_c_double)*sqrtg)
        bcontrav(1) = 0.0_c_double
        bcontrav(2) = -(point%poloidal_flux(2) - &
            point%toroidal_flux(2)*point%lambda(4))*flux_scale
        bcontrav(3) = -point%toroidal_flux(2)* &
            (1.0_c_double+point%lambda(3))*flux_scale
        metric(:, 1) = [dot_product(e_s, e_s), &
            dot_product(e_theta, e_s), dot_product(e_zeta, e_s)]
        metric(:, 2) = [dot_product(e_s, e_theta), &
            dot_product(e_theta, e_theta), dot_product(e_zeta, e_theta)]
        metric(:, 3) = [dot_product(e_s, e_zeta), &
            dot_product(e_theta, e_zeta), dot_product(e_zeta, e_zeta)]
        bcovar = matmul(metric, bcontrav)
        bvector = bcontrav(2)*e_theta + bcontrav(3)*e_zeta
        bmod = sqrt(dot_product(bvector, bvector))
        hcovar = bcovar/bmod
        hctrvr = bcontrav/bmod
    end subroutine field_components

    subroutine magfie_vmecpp(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        real(c_double), intent(in) :: x(3)
        real(c_double), intent(out) :: bmod, sqrtg
        real(c_double), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(c_double), parameter :: steps(3) = [1.0e-3_c_double, &
            2.0e-3_c_double*acos(-1.0_c_double), &
            4.0e-4_c_double*acos(-1.0_c_double)]
        real(c_double) :: shifted(3), bplus, bminus, unused_jacobian
        real(c_double) :: hcovar_plus(3), hcovar_minus(3), unused(3)
        real(c_double) :: derivatives(3, 3)
        integer :: coordinate

        do coordinate = 1, 3
            shifted = x
            shifted(coordinate) = x(coordinate) + steps(coordinate)
            call field_components(shifted, bplus, unused_jacobian, &
                hcovar_plus, unused)
            shifted(coordinate) = x(coordinate) - steps(coordinate)
            call field_components(shifted, bminus, unused_jacobian, &
                hcovar_minus, unused)
            bder(coordinate) = (bplus-bminus)/(2.0_c_double*steps(coordinate))
            derivatives(:, coordinate) = &
                (hcovar_plus-hcovar_minus)/(2.0_c_double*steps(coordinate))
        end do
        call field_components(x, bmod, sqrtg, hcovar, hctrvr)
        bder = bder/bmod
        hcurl = [derivatives(3, 2)-derivatives(2, 3), &
            derivatives(1, 3)-derivatives(3, 1), &
            derivatives(2, 1)-derivatives(1, 2)]/sqrtg

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
