module simple_magfie_base

use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none

type, abstract :: MagneticField
contains
    procedure(evaluate), deferred :: evaluate
end type MagneticField

abstract interface
    subroutine evaluate(self, x, Acov, hcov, sqgBctr, Bmod)
        import :: dp, MagneticField
        class(MagneticField), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3)
        real(dp), intent(out) :: hcov(3)
        real(dp), intent(out) :: sqgBctr(3)
        real(dp), intent(out) :: Bmod
    end subroutine evaluate
end interface


end module simple_magfie_base
