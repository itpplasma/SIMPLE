!> Factory module for ODE integrators
!>
!> Creates the appropriate integrator based on configuration.
module ode_integrator_factory
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use ode_integrator_base, only: ode_integrator_t
    use ode_integrator_libneo, only: libneo_integrator_t
    use ode_integrator_vode, only: vode_integrator_t
    use ode_integrator_rk4, only: rk4_integrator_t
    implicit none
    private

    public :: create_integrator, INTEGRATOR_LIBNEO, INTEGRATOR_VODE, INTEGRATOR_RK4

    integer, parameter :: INTEGRATOR_LIBNEO = 1
    integer, parameter :: INTEGRATOR_VODE = 2
    integer, parameter :: INTEGRATOR_RK4 = 3

contains

    function create_integrator(integrator_type) result(integrator)
        integer, intent(in) :: integrator_type
        class(ode_integrator_t), allocatable :: integrator

        select case (integrator_type)
        case (INTEGRATOR_LIBNEO)
            allocate(libneo_integrator_t :: integrator)
        case (INTEGRATOR_VODE)
            allocate(vode_integrator_t :: integrator)
        case (INTEGRATOR_RK4)
            allocate(rk4_integrator_t :: integrator)
        case default
            error stop "Unknown integrator type"
        end select
    end function create_integrator

end module ode_integrator_factory
