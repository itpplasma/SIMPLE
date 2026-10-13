module spectre_refraction
    !> Building blocks of the symplectic interface crossing (crossing_level 2).
    !>
    !> In canonical coordinates (theta, phi, p_theta, p_phi) of a volume chart
    !> the radius is the implicit function r(theta, phi, p_theta, p_phi) defined
    !> by p_theta = P(r, theta, phi, p_phi). The interface r = k is the level
    !> set s = r - k = 0, and its Hamiltonian vector field X_s moves the state
    !> along the interface at fixed r:
    !>
    !>     dtheta/dl = 1/P_r,  dphi/dl = -P_pphi/P_r,  dp_phi/dl = P_phi/P_r.
    !>
    !> A crossing step lands on the interface with a midpoint substep, flows
    !> along X_s for a length l and completes the step with a second midpoint
    !> substep. Choosing l so that the two substeps have equal discrete energy
    !> H(z_mid) makes the composite map exactly symplectic; the condition tends
    !> to H continuous across the jump as the step size vanishes.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use field_can_mod, only: field_can_t, eval_field => evaluate, get_val, &
        get_derivatives
    use odeint_allroutines_sub, only: odeint_allroutines

    implicit none
    private

    public :: midpoint_energy, interface_flow

    integer, parameter, public :: REFRACTION_OK = 0
    integer, parameter, public :: REFRACTION_DEGENERATE = 1
    integer, parameter, public :: REFRACTION_NO_CONVERGENCE = 2

    integer, parameter :: MAX_NEWTON = 40
    real(dp), parameter :: FLOW_RELERR = 1.0e-13_dp
    real(dp), parameter :: RADIUS_TOL = 1.0e-14_dp

    type :: flow_context_t
        real(dp) :: r = 0.0_dp
        real(dp) :: mu = 0.0_dp
        real(dp) :: ro0 = 0.0_dp
    end type flow_context_t

contains

    subroutine midpoint_energy(z0, pth0, z1, pth1, f, energy, status)
        !> H at the midpoint state of the midpoint step z0 -> z1 (same chart):
        !> theta, phi, p_phi and p_theta averaged, r solved from p_theta. f
        !> carries mu and ro0 and is overwritten with the midpoint field.
        real(dp), intent(in) :: z0(4), pth0, z1(4), pth1
        type(field_can_t), intent(inout) :: f
        real(dp), intent(out) :: energy
        integer, intent(out) :: status

        real(dp) :: zm(4), pth_mid, dr
        integer :: it

        zm = 0.5_dp*(z0 + z1)
        pth_mid = 0.5_dp*(pth0 + pth1)
        status = REFRACTION_NO_CONVERGENCE
        do it = 1, MAX_NEWTON
            call eval_field(f, zm(1), zm(2), zm(3), 0)
            call get_derivatives(f, zm(4))
            if (abs(f%dpth(1)) <= tiny(1.0_dp)) then
                status = REFRACTION_DEGENERATE
                return
            end if
            ! p_theta is a small difference of large terms, so the residual
            ! stalls at their round-off; converge on the radial update instead.
            dr = (f%pth - pth_mid)/f%dpth(1)
            zm(1) = zm(1) - dr
            if (abs(dr) <= RADIUS_TOL*max(abs(zm(1)), 1.0_dp)) then
                status = REFRACTION_OK
                exit
            end if
        end do
        call eval_field(f, zm(1), zm(2), zm(3), 0)
        call get_val(f, zm(4))
        energy = f%H
    end subroutine midpoint_energy

    subroutine interface_flow(z, length, f, status)
        !> Advance z = (r, theta, phi, p_phi) along X_s by the given length at
        !> fixed r. f supplies mu and ro0 and holds the end-point field.
        real(dp), intent(inout) :: z(4)
        real(dp), intent(in) :: length
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: status

        type(flow_context_t) :: context
        real(dp) :: y(3)
        integer :: ierr

        status = REFRACTION_OK
        if (length /= 0.0_dp) then
            context%r = z(1)
            context%mu = f%mu
            context%ro0 = f%ro0
            y = z(2:4)
            call odeint_allroutines(y, 3, context, 0.0_dp, length, FLOW_RELERR, &
                flow_rhs, ierr=ierr)
            if (ierr /= 0 .or. .not. all(ieee_is_finite(y))) then
                status = REFRACTION_NO_CONVERGENCE
                return
            end if
            z(2:4) = y
        end if
        call eval_field(f, z(1), z(2), z(3), 0)
        call get_derivatives(f, z(4))
    end subroutine interface_flow

    subroutine flow_rhs(l, y, dydl, context)
        real(dp), intent(in) :: l
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydl(:)
        class(*), intent(in) :: context

        type(field_can_t) :: f

        select type (context)
        type is (flow_context_t)
            f%mu = context%mu
            f%ro0 = context%ro0
            call eval_field(f, context%r, y(1), y(2), 0)
        end select
        call get_derivatives(f, y(3))
        dydl(1) = 1.0_dp/f%dpth(1)
        dydl(2) = -f%dpth(4)/f%dpth(1)
        dydl(3) = f%dpth(3)/f%dpth(1)
    end subroutine flow_rhs

end module spectre_refraction
