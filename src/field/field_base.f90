module field_base

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    type, abstract :: magnetic_field_t
    contains
        procedure(evaluate_interface), deferred :: evaluate
    end type magnetic_field_t

    abstract interface
        subroutine evaluate_interface(self, x, Acov, hcov, Bmod, sqgBctr)
            import :: dp, magnetic_field_t
            class(magnetic_field_t), intent(in) :: self
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: Acov(3)
            real(dp), intent(out) :: hcov(3)
            real(dp), intent(out) :: Bmod
            real(dp), intent(out), optional :: sqgBctr(3)
        end subroutine evaluate_interface
    end interface

end module field_base
