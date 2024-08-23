module simple_coordinates

use, intrinsic :: iso_fortran_env, only: dp => real64
use magfie_sub, only : TEST, CANFLUX, VMEC, BOOZER, MEISS

implicit none

integer, parameter :: CART = 111, CYL = 222

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



end subroutine transform_vmec_to_cyl


subroutine transform(from, to, xfrom, xto, dxto_dxfrom)
    integer, intent(in) :: from, to
    real(dp), intent(in) :: xfrom(3)
    real(dp), intent(out) :: xto(3)
    real(dp), intent(out) :: dxto_dxfrom(3,3)

    select case (from)
    case (VMEC)
        select case (to)
        case (CYL)
            call transform_vmec_to_cyl(xfrom, xto, dxto_dxfrom)
        end select
    end select
end subroutine transform


end module simple_coordinates
