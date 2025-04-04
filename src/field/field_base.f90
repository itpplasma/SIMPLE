module field_base

use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none

type, abstract :: MagneticField
contains
    procedure(evaluate), deferred :: evaluate
end type MagneticField

abstract interface
    subroutine evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
        import :: dp, MagneticField
        class(MagneticField), intent(in) :: self
        real(dp), intent(in) :: x(3)  ! r=sqrt(s_vmec), theta_vmec, phi_vmec
        real(dp), intent(out) :: Acov(3)
        real(dp), intent(out) :: hcov(3)
        real(dp), intent(out) :: Bmod
        real(dp), intent(out), optional :: sqgBctr(3)
    end subroutine evaluate
end interface


end module field_base
