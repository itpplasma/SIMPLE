module sheet_setup
    !> Marker initialization on the energy shell H = SPEED^2/2 at fixed p_phi.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use sheet_field, only: field_params_t, psi_value
    implicit none
    private

    public :: initial_state, shell_state

    real(dp), parameter :: TWOPI = 6.283185307179586_dp

contains

    subroutine initial_state(par, speed, im, z)
        !> Deterministic start just inside the sheet, poloidal angle and
        !> velocity direction spread over the markers.
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: speed
        integer, intent(in) :: im
        real(dp), intent(out) :: z(4)

        real(dp) :: theta, angle
        logical :: ok

        theta = TWOPI*(real(im, dp) - 0.5_dp)/8.0_dp
        angle = TWOPI*0.37_dp*real(im, dp)
        call shell_state(par, speed, par%r_sheet - 0.02_dp, theta, angle, z, ok)
        if (.not. ok) error stop 'initial_state: start outside the energy shell'
    end subroutine initial_state

    pure subroutine shell_state(par, speed, r, theta, angle, z, ok)
        !> Canonical state at minor radius r and poloidal angle theta with
        !> |v| = speed, the poloidal velocity direction given by angle and v_phi
        !> fixed by p_phi. ok = .false. where v_phi alone exceeds speed.
        type(field_params_t), intent(in) :: par
        real(dp), intent(in) :: speed, r, theta, angle
        real(dp), intent(out) :: z(4)
        logical, intent(out) :: ok

        real(dp) :: rr, zz, psi, v_phi, v_pol

        rr = par%r_major + r*cos(theta)
        zz = r*sin(theta)
        psi = psi_value(par, merge(1, 2, r < par%r_sheet), [rr, zz])
        v_phi = (par%p_phi - psi)/rr
        z = 0.0_dp
        ok = abs(v_phi) < speed
        if (.not. ok) return
        v_pol = sqrt(speed*speed - v_phi*v_phi)
        z(1) = rr
        z(2) = zz
        z(3) = v_pol*cos(angle)
        z(4) = v_pol*sin(angle) - par%f_tor*log(rr)
    end subroutine shell_state

end module sheet_setup
