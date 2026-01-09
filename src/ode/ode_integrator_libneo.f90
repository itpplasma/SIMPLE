!> Libneo ODE integrator adapter
!>
!> Wraps libneo's odeint_allroutines (Cash-Karp RK45) in the abstract interface.
module ode_integrator_libneo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use ode_integrator_base, only: ode_integrator_t, ode_rhs_interface
    implicit none
    private

    public :: libneo_integrator_t

    !> Libneo integrator implementation
    type, extends(ode_integrator_t) :: libneo_integrator_t
    contains
        procedure :: integrate => libneo_integrate
        procedure :: get_name => libneo_get_name
    end type libneo_integrator_t

    !> Module-level storage for RHS callback (needed for wrapper)
    procedure(ode_rhs_interface), pointer, save :: current_rhs => null()
    class(*), pointer, save :: current_context => null()
    !$omp threadprivate(current_rhs, current_context)

contains

    subroutine libneo_integrate(self, y, n, context, x1, x2, relerr, rhs)
        use odeint_allroutines_sub, only: odeint_allroutines
        class(libneo_integrator_t), intent(inout) :: self
        real(dp), intent(inout) :: y(:)
        integer, intent(in) :: n
        class(*), intent(in), target :: context
        real(dp), intent(in) :: x1, x2
        real(dp), intent(in) :: relerr
        procedure(ode_rhs_interface) :: rhs

        current_rhs => rhs
        current_context => context

        call odeint_allroutines(y, n, context, x1, x2, relerr, rhs_wrapper)

        current_rhs => null()
        current_context => null()
    end subroutine libneo_integrate

    subroutine rhs_wrapper(x, y, dydx, ctx)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydx(:)
        class(*), intent(in) :: ctx

        call current_rhs(x, y, dydx, ctx)
    end subroutine rhs_wrapper

    function libneo_get_name(self) result(name)
        class(libneo_integrator_t), intent(in) :: self
        character(len=32) :: name
        name = "libneo (Cash-Karp RK45)"
    end function libneo_get_name

end module ode_integrator_libneo
