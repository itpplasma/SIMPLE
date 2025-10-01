module simple_coordinates

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use vmec_coordinates, only: vmec_to_cyl_lib => vmec_to_cyl, &
                                 vmec_to_cart_lib => vmec_to_cart, &
                                 cyl_to_cart_lib  => cyl_to_cart
    use geoflux_coordinates, only: geoflux_to_cyl_lib => geoflux_to_cyl, &
                                     geoflux_to_cart_lib => geoflux_to_cart, &
                                     cyl_to_geoflux_lib => cyl_to_geoflux

    implicit none

    abstract interface
        subroutine transform_i(xfrom, xto, dxto_dxfrom)
            import :: dp
            real(dp), intent(in) :: xfrom(3)
            real(dp), intent(out) :: xto(3)
            real(dp), intent(out), optional :: dxto_dxfrom(3,3)
        end subroutine transform_i
    end interface

contains

subroutine transform_vmec_to_cyl(xfrom, xto, dxto_dxfrom)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out), optional :: dxto_dxfrom(3,3)

    call vmec_to_cyl_lib(xfrom, xto, dxto_dxfrom)
end subroutine transform_vmec_to_cyl


subroutine transform_vmec_to_cart(xfrom, xto, dxto_dxfrom)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out), optional :: dxto_dxfrom(3,3)

    call vmec_to_cart_lib(xfrom, xto, dxto_dxfrom)
end subroutine transform_vmec_to_cart


subroutine transform_cyl_to_cart(xfrom, xto, dxto_dxfrom)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out), optional :: dxto_dxfrom(3,3)

    call cyl_to_cart_lib(xfrom, xto, dxto_dxfrom)
end subroutine transform_cyl_to_cart


subroutine transform_geoflux_to_cyl(xfrom, xto, dxto_dxfrom)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out), optional :: dxto_dxfrom(3,3)

    call geoflux_to_cyl_lib(xfrom, xto, dxto_dxfrom)
end subroutine transform_geoflux_to_cyl


subroutine transform_geoflux_to_cart(xfrom, xto, dxto_dxfrom)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out), optional :: dxto_dxfrom(3,3)

    call geoflux_to_cart_lib(xfrom, xto, dxto_dxfrom)
end subroutine transform_geoflux_to_cart


subroutine transform_cyl_to_geoflux(xfrom, xto, dxto_dxfrom)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out), optional :: dxto_dxfrom(3,3)

    call cyl_to_geoflux_lib(xfrom, xto, dxto_dxfrom)
end subroutine transform_cyl_to_geoflux


function get_transform(from, to)
    procedure(transform_i), pointer :: get_transform
    character(*), intent(in) :: from, to

    get_transform => null()

    select case (trim(from))
    case('cyl')
        select case (trim(to))
        case ('cart')
            get_transform => transform_cyl_to_cart
        case ('geoflux')
            get_transform => transform_cyl_to_geoflux
        case default
            call handle_transform_error(from, to)
        end select
    case ('vmec')
        select case (trim(to))
        case ('cart')
            get_transform => transform_vmec_to_cart
        case ('cyl')
            get_transform => transform_vmec_to_cyl
        case default
            call handle_transform_error(from, to)
        end select
    case ('geoflux')
        select case (trim(to))
        case ('cyl')
            get_transform => transform_geoflux_to_cyl
        case ('cart')
            get_transform => transform_geoflux_to_cart
        case default
            call handle_transform_error(from, to)
        end select
    case default
        print *, "get_transform: Unknown transform from ", from
        error stop
    end select
end function get_transform


subroutine handle_transform_error(from, to)
    character(*), intent(in) :: from
    character(*), intent(in), optional :: to

    if (present(to)) then
        print *, "Unknown transform from ", from, " to ", to
    else
        print *, "Unknown transform from ", from
    end if
    error stop
end subroutine handle_transform_error

end module simple_coordinates
