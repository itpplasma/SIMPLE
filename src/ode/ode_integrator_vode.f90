!> VODE ODE integrator adapter
!>
!> Wraps DVODE_F90 (stiff/nonstiff solver with BDF/Adams methods) in the abstract interface.
!> VODE is better suited for stiff problems than explicit RK methods.
module ode_integrator_vode
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use ode_integrator_base, only: ode_integrator_t, ode_rhs_interface
    implicit none
    private

    public :: vode_integrator_t

    !> VODE integrator implementation
    type, extends(ode_integrator_t) :: vode_integrator_t
        integer :: mf = 22  ! Method flag: 22 = stiff (BDF) with internal Jacobian
        integer :: max_steps = 10000000  ! Maximum steps (10M, much higher than libneo)
    contains
        procedure :: integrate => vode_integrate
        procedure :: get_name => vode_get_name
    end type vode_integrator_t

    !> Module-level storage for RHS callback (needed for VODE interface)
    procedure(ode_rhs_interface), pointer, save :: current_rhs => null()
    class(*), pointer, save :: current_context => null()
    !$omp threadprivate(current_rhs, current_context)

contains

    subroutine vode_integrate(self, y, n, context, x1, x2, relerr, rhs)
        use DVODE_F90_M, only: DVODE_F90, SET_INTERMEDIATE_OPTS, VODE_OPTS, RELEASE_ARRAYS
        class(vode_integrator_t), intent(inout) :: self
        real(dp), intent(inout) :: y(:)
        integer, intent(in) :: n
        class(*), intent(in), target :: context
        real(dp), intent(in) :: x1, x2
        real(dp), intent(in) :: relerr
        procedure(ode_rhs_interface) :: rhs

        real(dp) :: t, tout
        real(dp) :: atol
        integer :: itask, istate
        type(VODE_OPTS) :: options

        current_rhs => rhs
        current_context => context

        t = x1
        tout = x2
        atol = relerr * 1.0d-3  ! Absolute tolerance (smaller than relative)
        itask = 1  ! Normal computation: integrate to tout
        istate = 1  ! First call

        options = SET_INTERMEDIATE_OPTS( &
            DENSE_J=.true., &
            ABSERR=atol, &
            RELERR=relerr, &
            MXSTEP=self%max_steps &
        )

        call DVODE_F90(vode_rhs_wrapper, n, y, t, tout, itask, istate, options)

        if (istate < 0) then
            write(*, '(A,I0)') "VODE error: ISTATE = ", istate
            if (istate == -1) then
                error stop "VODE: Excessive work done (max steps exceeded)"
            else if (istate == -2) then
                error stop "VODE: Too much accuracy requested"
            else if (istate == -3) then
                error stop "VODE: Illegal input detected"
            else if (istate == -4) then
                error stop "VODE: Repeated error test failures"
            else if (istate == -5) then
                error stop "VODE: Repeated convergence failures"
            else if (istate == -6) then
                error stop "VODE: Error weight became zero"
            else
                error stop "VODE: Unknown error"
            end if
        end if

        call RELEASE_ARRAYS()

        current_rhs => null()
        current_context => null()
    end subroutine vode_integrate

    subroutine vode_rhs_wrapper(neq, t, y, ydot)
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(neq)
        real(dp), intent(out) :: ydot(neq)

        call current_rhs(t, y, ydot, current_context)
    end subroutine vode_rhs_wrapper

    function vode_get_name(self) result(name)
        class(vode_integrator_t), intent(in) :: self
        character(len=32) :: name
        name = "VODE (BDF stiff solver)"
    end function vode_get_name

end module ode_integrator_vode
