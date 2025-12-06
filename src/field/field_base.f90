module field_base
    !> Abstract base type for magnetic field implementations.
    !> Every field has a coordinate system in which it naturally evaluates.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t

    implicit none

    type, abstract :: magnetic_field_t
        class(coordinate_system_t), allocatable :: coords
    contains
        procedure(evaluate_interface), deferred :: evaluate
    end type magnetic_field_t

    abstract interface
        subroutine evaluate_interface(self, x, Acov, hcov, Bmod, sqgBctr)
            !> Evaluate field at x (in self%coords coordinate system).
            !> Returns covariant components in self%coords.
            import :: dp, magnetic_field_t
            class(magnetic_field_t), intent(in) :: self
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: Acov(3), hcov(3), Bmod
            real(dp), intent(out), optional :: sqgBctr(3)
        end subroutine evaluate_interface
    end interface

end module field_base
