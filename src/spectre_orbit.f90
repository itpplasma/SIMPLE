module spectre_orbit
    !> Single-volume SPECTRE RK45 stepping, with an absorbing outer boundary
    !> and reflection before the singular coordinate axis.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use odeint_allroutines_sub, only: odeint_allroutines, odeint_has_failed
    use alpha_lifetime_sub, only: velo_can

    implicit none
    private

    public :: spectre_orbit_state_t, spectre_event_t, spectre_state_reset, &
        orbit_timestep_spectre, SPECTRE_OK, SPECTRE_BOUNDARY, SPECTRE_FAULT

    integer, parameter :: SPECTRE_OK = 0
    integer, parameter :: SPECTRE_BOUNDARY = 88
    integer, parameter :: SPECTRE_FAULT = 2

    real(dp), parameter :: axis_offset = 1.0d-3
    real(dp), parameter :: rho_tol = 1.0d-10
    integer, parameter :: max_bisect = 80
    integer, parameter :: max_retry_depth = 8

    !> Boundary limits of the single-volume clamped RHS, private to each thread.
    real(dp) :: ev_lo = 0.0_dp, ev_hi = 0.0_dp
    !$omp threadprivate(ev_lo, ev_hi)

    !> Keep the first accepted state strictly inside the outer boundary.
    real(dp), parameter :: clamp_margin = 1.0d-9
    !> Preserve the existing single-volume boundary regularization.
    real(dp), parameter :: drift_band = 2.0d-2
    !> Backstop for the same Bstar_par -> 0 breakdown when it falls inside the
    !> band (it is a pitch-dependent phase-space surface, not a fixed radius):
    !> cap the RHS magnitude so odeint keeps a finite step through the unphysical
    !> layer instead of collapsing to zero step and exhausting its step budget.
    real(dp), parameter :: drift_cap = 1.0d0

    type :: spectre_orbit_state_t
        logical :: initialized = .false.
    end type spectre_orbit_state_t

    type :: spectre_event_t
        real(dp) :: t_frac = 0.0_dp
    end type spectre_event_t

contains

    subroutine spectre_state_reset(state, mvol)
        type(spectre_orbit_state_t), intent(out) :: state
        integer, intent(in) :: mvol

        if (mvol /= 1) error stop 'SPECTRE orbit stepping requires one volume'
        state%initialized = .false.
    end subroutine spectre_state_reset

    subroutine orbit_timestep_spectre(state, z, dtaumin, relerr, ierr, event)
        !> Advance z = (rho_g, theta, zeta, p, lambda) by one microstep.
        type(spectre_orbit_state_t), intent(inout) :: state
        real(dp), intent(inout) :: z(5)
        real(dp), intent(in) :: dtaumin, relerr
        integer, intent(out) :: ierr
        type(spectre_event_t), intent(out) :: event

        real(dp) :: z_start(5), z_end(5), z_hit(5), boundary
        integer :: direction

        ierr = SPECTRE_OK
        z_start = z

        if (.not. state%initialized) then
            ev_lo = axis_offset
            ev_hi = 1.0_dp
            call integrate_clamped(z_start, dtaumin, relerr, z_end, ierr)
            if (ierr /= SPECTRE_OK) return
            state%initialized = .true.
            z = z_end
            z(1) = max(0.0_dp, min(z(1), 1.0_dp - clamp_margin))
            return
        end if

        ev_lo = 0.0_dp
        ev_hi = 1.0_dp
        call integrate_clamped(z_start, dtaumin, relerr, z_end, ierr)
        if (ierr /= SPECTRE_OK) return

        if (z_end(1) > 1.0_dp .or. &
            (z_end(1) == 1.0_dp .and. z_start(1) < 1.0_dp)) then
            boundary = 1.0_dp
            direction = 1
        else if (z_end(1) < 0.0_dp .or. &
                (z_end(1) == 0.0_dp .and. z_start(1) > 0.0_dp)) then
            boundary = 0.0_dp
            direction = -1
        else
            z = z_end
            return
        end if

        call locate_boundary(z_start, dtaumin, relerr, boundary, z_hit, &
            event%t_frac, ierr)
        if (ierr /= SPECTRE_OK) then
            z = z_start
            return
        end if
        z = z_hit
        if (direction > 0) then
            z(1) = 1.0_dp
            ierr = SPECTRE_BOUNDARY
        else
            z(1) = axis_offset
            z(5) = -z_hit(5)
        end if
    end subroutine orbit_timestep_spectre

    subroutine locate_boundary(z_start, dtaumin, relerr, boundary, z_hit, t_frac, &
            ierr)
        !> Find the microstep time where the clamped-RHS trajectory reaches rho_g =
        !> boundary to |rho_g - boundary| < rho_tol, by the Illinois variant of
        !> false position. z_start is strictly inside the home volume and the full
        !> step already overshot the boundary, so [0, dtaumin] brackets the root.
        real(dp), intent(in) :: z_start(5), dtaumin, relerr, boundary
        real(dp), intent(out) :: z_hit(5)
        real(dp), intent(out) :: t_frac
        integer, intent(out) :: ierr

        real(dp) :: t_lo, t_hi, f_lo, f_hi, t_mid, f_mid, z_mid(5)
        integer :: it, last_side
        logical :: converged

        t_lo = 0.0_dp
        ierr = SPECTRE_OK
        f_lo = z_start(1) - boundary
        t_hi = dtaumin
        call integrate_clamped(z_start, t_hi, relerr, z_mid, ierr)
        if (ierr /= SPECTRE_OK) return
        f_hi = z_mid(1) - boundary
        z_hit = z_mid
        t_mid = t_hi
        last_side = 0
        converged = abs(f_hi) < rho_tol

        do it = 1, max_bisect
            if (converged) exit
            if (abs(f_hi - f_lo) <= tiny(1.0_dp)) exit
            t_mid = t_hi - f_hi*(t_hi - t_lo)/(f_hi - f_lo)
            call integrate_clamped(z_start, t_mid, relerr, z_mid, ierr)
            if (ierr /= SPECTRE_OK) return
            f_mid = z_mid(1) - boundary
            z_hit = z_mid
            if (abs(f_mid) < rho_tol) then
                converged = .true.
                exit
            end if
            if (f_mid*f_lo > 0.0_dp) then
                t_lo = t_mid
                f_lo = f_mid
                if (last_side == 1) f_hi = 0.5_dp*f_hi
                last_side = 1
            else
                t_hi = t_mid
                f_hi = f_mid
                if (last_side == -1) f_lo = 0.5_dp*f_lo
                last_side = -1
            end if
        end do

        if (.not. converged) then
            z_hit = z_start
            t_frac = 0.0_dp
            ierr = SPECTRE_FAULT
            return
        end if
        t_frac = t_mid/dtaumin
    end subroutine locate_boundary

    subroutine velo_can_clamped(tau, z, vz)
        !> velo_can with rho_g clamped to the home volume [ev_lo, ev_hi). Inside
        !> the volume the drift is exact; past a boundary the field is held at
        !> its boundary value. The final state can move beyond the boundary,
        !> allowing the caller to locate the loss or axis-reflection event.
        real(dp), intent(in) :: tau
        real(dp), intent(in) :: z(:)
        real(dp), intent(out) :: vz(:)

        real(dp) :: zc(size(z)), vmag

        zc = z
        zc(1) = max(max(ev_lo + drift_band, axis_offset), &
            min(z(1), ev_hi - drift_band))
        call velo_can(tau, zc, vz)
        vmag = sqrt(sum(vz**2))
        if (vmag > drift_cap) vz = vz*(drift_cap/vmag)
    end subroutine velo_can_clamped

    subroutine integrate_clamped(z_in, tau, relerr, z_out, ierr)
        real(dp), intent(in) :: z_in(5), tau, relerr
        real(dp), intent(out) :: z_out(5)
        integer, intent(out) :: ierr

        call integrate_clamped_retry(z_in, tau, relerr, 0, z_out, ierr)
    end subroutine integrate_clamped

    recursive subroutine integrate_clamped_retry(z_in, tau, relerr, depth, &
            z_out, ierr)
        real(dp), intent(in) :: z_in(5), tau, relerr
        integer, intent(in) :: depth
        real(dp), intent(out) :: z_out(5)
        integer, intent(out) :: ierr

        real(dp) :: z_half(5)
        integer, parameter :: ndim = 5

        ierr = SPECTRE_OK
        z_out = z_in
        if (tau <= 0.0_dp) return
        call odeint_allroutines(z_out, ndim, 0.0_dp, tau, relerr, velo_can_clamped)
        if (.not. odeint_has_failed() .and. all(ieee_is_finite(z_out))) return

        z_out = z_in
        ierr = SPECTRE_FAULT
        if (depth >= max_retry_depth) return
        call integrate_clamped_retry(z_in, 0.5_dp*tau, relerr, depth + 1, &
            z_half, ierr)
        if (ierr /= SPECTRE_OK) return
        call integrate_clamped_retry(z_half, 0.5_dp*tau, relerr, depth + 1, &
            z_out, ierr)
        if (ierr /= SPECTRE_OK) z_out = z_in
    end subroutine integrate_clamped_retry

end module spectre_orbit
