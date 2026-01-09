!> Simple fixed-step RK4 ODE integrator
!>
!> Classic 4th-order Runge-Kutta method with fixed step size.
!> Does not detect stiffness or adapt step size - just integrates robustly.
module ode_integrator_rk4
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use ode_integrator_base, only: ode_integrator_t, ode_rhs_interface
    implicit none
    private

    public :: rk4_integrator_t

    !> RK4 integrator implementation
    type, extends(ode_integrator_t) :: rk4_integrator_t
        integer :: min_steps = 100  ! Minimum steps per integration interval
    contains
        procedure :: integrate => rk4_integrate
        procedure :: get_name => rk4_get_name
    end type rk4_integrator_t

contains

    subroutine rk4_integrate(self, y, n, context, x1, x2, relerr, rhs)
        class(rk4_integrator_t), intent(inout) :: self
        real(dp), intent(inout) :: y(:)
        integer, intent(in) :: n
        class(*), intent(in), target :: context
        real(dp), intent(in) :: x1, x2
        real(dp), intent(in) :: relerr
        procedure(ode_rhs_interface) :: rhs

        real(dp) :: x, h
        real(dp), dimension(size(y)) :: k1, k2, k3, k4, y_temp
        integer :: nsteps, i

        nsteps = max(self%min_steps, ceiling(1.0d0 / relerr))
        h = (x2 - x1) / dble(nsteps)

        x = x1
        do i = 1, nsteps
            call rhs(x, y, k1, context)
            y_temp = y + 0.5d0 * h * k1

            call rhs(x + 0.5d0 * h, y_temp, k2, context)
            y_temp = y + 0.5d0 * h * k2

            call rhs(x + 0.5d0 * h, y_temp, k3, context)
            y_temp = y + h * k3

            call rhs(x + h, y_temp, k4, context)

            y = y + h * (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4) / 6.0d0
            x = x + h
        end do
    end subroutine rk4_integrate

    function rk4_get_name(self) result(name)
        class(rk4_integrator_t), intent(in) :: self
        character(len=32) :: name
        name = "RK4 (fixed-step)"
    end function rk4_get_name

end module ode_integrator_rk4
