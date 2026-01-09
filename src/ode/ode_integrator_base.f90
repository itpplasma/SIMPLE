!> Abstract base class for ODE integrators
!>
!> Provides a common interface for different ODE solvers (libneo, VODE, etc.)
!> Uses adapter pattern to allow runtime selection of integrator.
module ode_integrator_base
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: ode_integrator_t, ode_rhs_interface

    !> Abstract interface for ODE right-hand side function
    abstract interface
        subroutine ode_rhs_interface(x, y, dydx, context)
            import :: dp
            real(dp), intent(in) :: x
            real(dp), intent(in) :: y(:)
            real(dp), intent(out) :: dydx(:)
            class(*), intent(in) :: context
        end subroutine ode_rhs_interface
    end interface

    !> Abstract base class for ODE integrators
    type, abstract :: ode_integrator_t
    contains
        procedure(integrate_interface), deferred :: integrate
        procedure(get_name_interface), deferred :: get_name
    end type ode_integrator_t

    abstract interface
        subroutine integrate_interface(self, y, n, context, x1, x2, relerr, rhs)
            import :: dp, ode_integrator_t, ode_rhs_interface
            class(ode_integrator_t), intent(inout) :: self
            real(dp), intent(inout) :: y(:)
            integer, intent(in) :: n
            class(*), intent(in), target :: context
            real(dp), intent(in) :: x1, x2
            real(dp), intent(in) :: relerr
            procedure(ode_rhs_interface) :: rhs
        end subroutine integrate_interface

        function get_name_interface(self) result(name)
            import :: ode_integrator_t
            class(ode_integrator_t), intent(in) :: self
            character(len=32) :: name
        end function get_name_interface
    end interface

end module ode_integrator_base
