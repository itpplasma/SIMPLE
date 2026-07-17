module spectre_orbit
    !> Non-canonical RK45 guiding-center stepping in the SPECTRE stacked-rho
    !> chart with Level-0 interface crossing (#438, #443).
    !>
    !> A marker is traced inside one volume at a time on a field-clamped RHS whose
    !> |B| is frozen at the home-volume edge, so the RK45 substep never marches
    !> across the interface |B| jump. When the substep leaves [home_lo, home_hi]
    !> the interface time is found by Illinois false position and apply_crossing
    !> switches volume (or reflects) before stepping resumes. The outermost
    !> interface is the loss surface: a crossing outward from volume Mvol reports
    !> SPECTRE_BOUNDARY and terminates the trace.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use odeint_allroutines_sub, only: odeint_allroutines, odeint_has_failed
    use alpha_lifetime_sub, only: velo_can
    use interface_crossing, only: apply_crossing, crossing_info_t, axis_offset, &
        CROSS_LOSS

    implicit none
    private

    public :: spectre_orbit_state_t, spectre_event_t, spectre_state_reset, &
        orbit_timestep_spectre, SPECTRE_OK, SPECTRE_BOUNDARY, SPECTRE_FAULT

    integer, parameter :: SPECTRE_OK = 0
    integer, parameter :: SPECTRE_BOUNDARY = 88
    integer, parameter :: SPECTRE_FAULT = 2

    real(dp), parameter :: rho_tol = 1.0d-10
    integer, parameter :: max_bisect = 80
    integer, parameter :: max_retry_depth = 8

    !> Home-volume boundaries handed to the field-clamped RHS. Each traced marker
    !> owns one thread, so these must be thread-private.
    real(dp) :: ev_lo = 0.0_dp, ev_hi = 0.0_dp
    !$omp threadprivate(ev_lo, ev_hi)

    !> Keep an initially sampled marker strictly inside its inferred home volume.
    real(dp), parameter :: clamp_margin = 1.0d-9
    !> The guiding-center drift diverges (Bstar_par -> 0 from the interface sheet
    !> current) in a thin unphysical layer next to each interface. The RHS field
    !> is frozen at the volume edge minus this band so the substep integrates a
    !> finite, physical drift; the exact interface jump is applied separately by
    !> apply_crossing at rho_g = iface.
    real(dp), parameter :: drift_band = 2.0d-2
    !> Backstop for the same Bstar_par -> 0 breakdown when it falls inside the
    !> band (it is a pitch-dependent phase-space surface, not a fixed radius):
    !> cap the RHS magnitude so odeint keeps a finite step through the unphysical
    !> layer instead of collapsing to zero step and exhausting its step budget.
    real(dp), parameter :: drift_cap = 1.0d0

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
        real(dp) :: t_frac = 0.0_dp
        type(crossing_info_t) :: info
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

    subroutine orbit_timestep_spectre(state, z, dtaumin, relerr, level, ierr, event)
        !> Advance z = (rho_g, theta, zeta, p, lambda) by one microstep dtaumin.
        !> Once the home volume is fixed the substep runs on the field-clamped RHS
        !> (continuous across the interface); if the trajectory leaves the home
        !> volume the crossing time is bisected to the interface and apply_crossing
        !> switches volume or reflects with the Level-`level` map, retaining exact
        !> interface coordinates with explicit side ownership. A crossing outward
        !> through the outermost interface is a loss and reports SPECTRE_BOUNDARY.
        type(spectre_orbit_state_t), intent(inout) :: state
        real(dp), intent(inout) :: z(5)
        real(dp), intent(in) :: dtaumin, relerr
        integer, intent(in) :: level
        integer, intent(out) :: ierr
        type(spectre_event_t), intent(out) :: event

        real(dp) :: z_start(5), z_end(5), z_hit(5), boundary
        integer :: direction, iface

        event%occurred = .false.
        ierr = SPECTRE_OK
        z_start = z

        if (.not. state%home_set) then
            ev_lo = axis_offset
            ev_hi = real(state%mvol, dp)
            call integrate_clamped(z_start, dtaumin, relerr, z_end, ierr)
            if (ierr /= SPECTRE_OK) return
            call set_home_volume(state, z_end(1))
            z = z_end
            z(1) = max(state%home_lo, min(z(1), state%home_hi - clamp_margin))
            return
        end if

        ev_lo = state%home_lo
        ev_hi = state%home_hi
        call integrate_clamped(z_start, dtaumin, relerr, z_end, ierr)
        if (ierr /= SPECTRE_OK) return

        if (z_end(1) > state%home_hi .or. &
            (z_end(1) == state%home_hi .and. z_start(1) < state%home_hi)) then
            boundary = state%home_hi
            direction = 1
        else if (z_end(1) < state%home_lo .or. &
                (z_end(1) == state%home_lo .and. z_start(1) > state%home_lo)) then
            boundary = state%home_lo
            direction = -1
        else
            z = z_end
            return
        end if

        call locate_crossing(z_start, dtaumin, relerr, boundary, z_hit, &
            event%t_frac, ierr)
        if (ierr /= SPECTRE_OK) then
            z = z_start
            event%occurred = .false.
            return
        end if
        iface = nint(boundary)
        call apply_crossing(z_hit, iface, direction, state%mvol, level, &
            z, event%info)
        event%occurred = .true.

        if (event%info%event_type == CROSS_LOSS) then
            ierr = SPECTRE_BOUNDARY
        else
            call set_home_volume_index(state, event%info%vol_to)
            ierr = SPECTRE_OK
        end if
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

    subroutine set_home_volume_index(state, lvol)
        type(spectre_orbit_state_t), intent(inout) :: state
        integer, intent(in) :: lvol

        state%home_lo = real(lvol - 1, dp)
        state%home_hi = real(lvol, dp)
        state%home_set = .true.
    end subroutine set_home_volume_index

    subroutine locate_crossing(z_start, dtaumin, relerr, boundary, z_hit, t_frac, &
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
    end subroutine locate_crossing

    subroutine velo_can_clamped(tau, z, vz)
        !> velo_can with rho_g clamped to the home volume [ev_lo, ev_hi). Inside
        !> the volume the drift is exact; past an interface the field is held at
        !> its boundary value so the RHS stays continuous and the RK45 step never
        !> straddles the |B| jump. z_end itself still moves past the interface, so
        !> the caller detects the crossing and bisects to the exact interface.
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
