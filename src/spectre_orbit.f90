module spectre_orbit
    !> Non-canonical RK45 guiding-center stepping in the SPECTRE stacked-rho
    !> chart with per-volume boundary termination (#438).
    !>
    !> A marker is traced inside one volume at a time: rho_g moves freely until it
    !> reaches an interface (an integer rho_g). Interface crossing physics is
    !> deferred to spectre-08, so the orbit stops there ("boundary-stop"), the
    !> crossing time is located by bisection of the dense RK45 trajectory to
    !> |rho_g - integer| < 1e-10, and the terminal state is reported to the caller.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use odeint_allroutines_sub, only: odeint_allroutines
    use alpha_lifetime_sub, only: velo_can

    implicit none
    private

    public :: spectre_orbit_state_t, spectre_event_t, spectre_state_reset, &
              orbit_timestep_spectre, SPECTRE_OK, SPECTRE_BOUNDARY, SPECTRE_FAULT

    integer, parameter :: SPECTRE_OK = 0
    integer, parameter :: SPECTRE_BOUNDARY = 88
    integer, parameter :: SPECTRE_FAULT = 2

    real(dp), parameter :: rho_tol = 1.0d-10
    integer, parameter :: max_bisect = 80

    type :: spectre_orbit_state_t
        !> Home volume [home_lo, home_hi] is the pair of consecutive integer
        !> interfaces bracketing the marker; it is fixed after the first step so a
        !> start placed exactly on an interface picks the volume it enters.
        logical :: home_set = .false.
        real(dp) :: home_lo = 0.0_dp
        real(dp) :: home_hi = 0.0_dp
        integer :: mvol = 1
    end type spectre_orbit_state_t

    type :: spectre_event_t
        logical :: occurred = .false.
        integer :: iface = 0
        integer :: direction = 0
        real(dp) :: t_frac = 0.0_dp
    end type spectre_event_t

contains

    subroutine spectre_state_reset(state, mvol)
        type(spectre_orbit_state_t), intent(out) :: state
        integer, intent(in) :: mvol

        state%home_set = .false.
        state%home_lo = 0.0_dp
        state%home_hi = 0.0_dp
        state%mvol = mvol
    end subroutine spectre_state_reset

    subroutine orbit_timestep_spectre(state, z, dtaumin, relerr, ierr, event)
        !> Advance z = (rho_g, theta, zeta, p, lambda) by one microstep dtaumin,
        !> then test the home-volume boundaries. On a crossing z is refined to the
        !> interface and the step reports SPECTRE_BOUNDARY.
        type(spectre_orbit_state_t), intent(inout) :: state
        real(dp), intent(inout) :: z(5)
        real(dp), intent(in) :: dtaumin, relerr
        integer, intent(out) :: ierr
        type(spectre_event_t), intent(out) :: event

        real(dp) :: z_start(5), z_end(5)
        real(dp) :: rho_start, rho_end, boundary

        event%occurred = .false.
        ierr = SPECTRE_OK
        rho_start = z(1)

        z_start = z
        call integrate_step(z_start, dtaumin, relerr, z_end, ierr)
        if (ierr /= SPECTRE_OK) return
        rho_end = z_end(1)

        if (.not. state%home_set) then
            call set_home_volume(state, rho_end)
            z = z_end
            return
        end if

        if (rho_end >= state%home_hi) then
            boundary = state%home_hi
            event%direction = 1
        else if (rho_end <= state%home_lo) then
            boundary = state%home_lo
            event%direction = -1
        else
            z = z_end
            return
        end if

        call locate_crossing(z_start, dtaumin, relerr, rho_start, boundary, z, &
                             event%t_frac)
        event%iface = nint(boundary)
        event%occurred = .true.
        ierr = SPECTRE_BOUNDARY
    end subroutine orbit_timestep_spectre

    subroutine set_home_volume(state, rho)
        type(spectre_orbit_state_t), intent(inout) :: state
        real(dp), intent(in) :: rho

        real(dp) :: lo

        lo = real(floor(rho), dp)
        lo = max(0.0_dp, min(lo, real(state%mvol - 1, dp)))
        state%home_lo = lo
        state%home_hi = lo + 1.0_dp
        state%home_set = .true.
    end subroutine set_home_volume

    subroutine locate_crossing(z_start, dtaumin, relerr, rho_start, boundary, &
                               z_hit, t_frac)
        !> Bisect the microstep time so the dense RK45 trajectory reaches rho_g =
        !> boundary to |rho_g - boundary| < 1e-10. rho_start is strictly inside.
        real(dp), intent(in) :: z_start(5), dtaumin, relerr, rho_start, boundary
        real(dp), intent(out) :: z_hit(5)
        real(dp), intent(out) :: t_frac

        real(dp) :: t_lo, t_hi, t_mid, side_start, rho_mid, z_mid(5)
        integer :: it, ierr

        side_start = sign(1.0_dp, rho_start - boundary)
        t_lo = 0.0_dp
        t_hi = dtaumin
        z_hit = z_start
        t_mid = dtaumin

        do it = 1, max_bisect
            t_mid = 0.5_dp*(t_lo + t_hi)
            call integrate_step(z_start, t_mid, relerr, z_mid, ierr)
            if (ierr /= SPECTRE_OK) exit
            rho_mid = z_mid(1)
            z_hit = z_mid
            if (abs(rho_mid - boundary) < rho_tol) exit
            if (sign(1.0_dp, rho_mid - boundary) == side_start) then
                t_lo = t_mid
            else
                t_hi = t_mid
            end if
        end do

        t_frac = t_mid/dtaumin
    end subroutine locate_crossing

    subroutine integrate_step(z_in, tau, relerr, z_out, ierr)
        real(dp), intent(in) :: z_in(5), tau, relerr
        real(dp), intent(out) :: z_out(5)
        integer, intent(out) :: ierr

        integer, parameter :: ndim = 5

        ierr = SPECTRE_OK
        z_out = z_in
        if (tau <= 0.0_dp) return
        call odeint_allroutines(z_out, ndim, 0.0_dp, tau, relerr, velo_can)
    end subroutine integrate_step

end module spectre_orbit
