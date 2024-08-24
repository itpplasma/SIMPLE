module simple_coordinates

use, intrinsic :: iso_fortran_env, only: dp => real64

implicit none

abstract interface
    subroutine transform_i(xfrom, xto, dxto_dxfrom)
        import :: dp
        real(dp), intent(in) :: xfrom(3)
        real(dp), intent(out) :: xto(3)
        real(dp), intent(out) :: dxto_dxfrom(3,3)
    end subroutine transform_i
end interface

contains

subroutine transform_vmec_to_cyl(xfrom, xto, dxto_dxfrom)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out) :: dxto_dxfrom(3,3)

    print *, "VMEC TO CYL"

end subroutine transform_vmec_to_cyl


subroutine transform_vmec_to_cart(xfrom, xto, dxto_dxfrom)
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out) :: dxto_dxfrom(3,3)

    print *, "VMEC TO CART"

end subroutine transform_vmec_to_cart


function get_transform(from, to)
    procedure(transform_i), pointer :: get_transform
    character(*), intent(in) :: from, to

    get_transform => null()

    select case (trim(from))
    case ('vmec')
        select case (trim(to))
        case ('cart')
            get_transform => transform_vmec_to_cart
        case ('cyl')
            get_transform => transform_vmec_to_cyl
        case default
            print *, "get_transform: Unknown transform from ", from, " to ", to
            error stop
        end select
    case default
        print *, "get_transform: Unknown transform from ", from
        error stop
    end select
end function get_transform


end module simple_coordinates
