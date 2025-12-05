module field_base

use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none

type, abstract :: magnetic_field_t
contains
    procedure(evaluate_iface), deferred :: evaluate
    procedure(evaluate_spline_input_iface), deferred :: evaluate_spline_input
end type magnetic_field_t

abstract interface
    subroutine evaluate_iface(self, x, Acov, hcov, Bmod, sqgBctr)
        import :: dp, magnetic_field_t
        class(magnetic_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod
        real(dp), intent(out), optional :: sqgBctr(3)
    end subroutine evaluate_iface

    subroutine evaluate_spline_input_iface(self, x, Acov, hcov, Bmod)
        !> Evaluate field for spline construction (no optional args).
        !> Used by build_field_splines to sample the field on a grid.
        import :: dp, magnetic_field_t
        class(magnetic_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod
    end subroutine evaluate_spline_input_iface
end interface

end module field_base
