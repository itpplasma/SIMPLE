module simple_magfie_coils

use, intrinsic :: iso_fortran_env, only: dp => real64
use neo_biotsavart, only: coils_t, load_coils_from_file, &
    compute_vector_potential, compute_magnetic_field
use simple_magfie_base, only: MagneticField
use simple_coordinates, only: transform_vmec_to_cart

implicit none

type, extends(MagneticField) :: CoilsField
    type(coils_t) :: coils
contains
    procedure :: evaluate
end type CoilsField

contains

subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
    class(CoilsField), intent(in) :: self
    real(dp), intent(in) :: x(3)
    real(dp), intent(out), dimension(3) :: Acov, hcov
    real(dp), intent(out) :: Bmod
    real(dp), intent(out), optional :: sqgBctr(3)

    real(dp) :: xcart(3), dxvmec_dxcart(3,3), dxcart_dxvmec(3,3)
    real(dp), dimension(3) :: Acart, Bcart

    call transform_vmec_to_cart(x, xcart, dxvmec_dxcart)
    dxcart_dxvmec = invert3(dxvmec_dxcart)

    Acart = compute_vector_potential(self%coils, xcart)
    Bcart = compute_magnetic_field(self%coils, xcart)

    Acov = matmul(Acart, dxcart_dxvmec)

    Bmod = sqrt(Bcart(1)**2 + Bcart(2)**2 + Bcart(3)**2)

    hcov = matmul(Bcart, dxcart_dxvmec)/Bmod

    if (present(sqgBctr)) then
        error stop 'sqgBctr not implemented'
    end if
end subroutine evaluate


function create_coils_field(coils_file) result(field)
    class(CoilsField), allocatable :: field
    character(*), intent(in) :: coils_file

    real(dp), parameter :: M_TO_CM = 100.0d0
    real(dp), parameter :: A_TO_STATA = 2997924536.8431d0

    allocate(CoilsField :: field)
    call load_coils_from_file(coils_file, field%coils)

    field%coils%x = field%coils%x * M_TO_CM
    field%coils%y = field%coils%y * M_TO_CM
    field%coils%z = field%coils%z * M_TO_CM
    field%coils%current = field%coils%current * A_TO_STATA

end function create_coils_field


function invert3(M) result(Minv)
    implicit none
    real(8), intent(in) :: M(3,3)
    real(8) :: Minv(3,3)
    real(8) :: det

    det = M(1,1)*(M(2,2)*M(3,3) - M(2,3)*M(3,2)) - &
          M(1,2)*(M(2,1)*M(3,3) - M(2,3)*M(3,1)) + &
          M(1,3)*(M(2,1)*M(3,2) - M(2,2)*M(3,1))

    Minv(1,1) = (M(2,2)*M(3,3) - M(2,3)*M(3,2)) / det
    Minv(1,2) = -(M(1,2)*M(3,3) - M(1,3)*M(3,2)) / det
    Minv(1,3) = (M(1,2)*M(2,3) - M(1,3)*M(2,2)) / det
    Minv(2,1) = -(M(2,1)*M(3,3) - M(2,3)*M(3,1)) / det
    Minv(2,2) = (M(1,1)*M(3,3) - M(1,3)*M(3,1)) / det
    Minv(2,3) = -(M(1,1)*M(2,3) - M(1,3)*M(2,1)) / det
    Minv(3,1) = (M(2,1)*M(3,2) - M(2,2)*M(3,1)) / det
    Minv(3,2) = -(M(1,1)*M(3,2) - M(1,2)*M(3,1)) / det
    Minv(3,3) = (M(1,1)*M(2,2) - M(1,2)*M(2,1)) / det
end function invert3


end module simple_magfie_coils
