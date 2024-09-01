module simple_magfie_coils

use, intrinsic :: iso_fortran_env, only: dp => real64
use simple_magfie_base, only: MagneticField
implicit none

type, extends(MagneticField) :: CoilsField
contains
    procedure :: evaluate
end type CoilsField

contains

subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    use simple_coordinates, only: transform_vmec_to_cyl

    class(CoilsField), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out), dimension(3) :: Acov, hcov
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)

    real(dp) :: xcyl(3), dxvmec_dxcyl(3,3)

    call transform_vmec_to_cyl(x, xcyl, dxvmec_dxcyl)

    Acov(1) = 0d0
    Acov(2) = 0d0
    Acov(3) = 0d0

    Bmod = 0d0

    hcov(1) = 0d0
    hcov(2) = 0d0
    hcov(3) = 0d0

    if (present(sqgBctr)) then
        error stop 'sqgBctr not implemented'
    end if
end subroutine evaluate

end module simple_magfie_coils
