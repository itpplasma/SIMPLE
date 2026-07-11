module spectre_sympl_orbit
    !> Symplectic multi-volume microstepping for SPECTRE (#441). When an accepted
    !> implicit step leaves the home volume, the step is rewound and re-taken from
    !> the same pre-step state with a shortened length h* so it lands on the
    !> interface rho_g = k exactly; the map for any fixed h* is symplectic, so the
    !> landed state is on the exact discrete orbit of a symplectic map. The
    !> crossing map (#440/#443) is then applied in physical variables and the
    !> integrator is re-initialized from the mapped state in the target volume's
    !> gauge, the same code path as an orbit start. The remaining microstep budget
    !> is completed with the new volume's data, so a microstep advances exactly
    !> dtaumin and step-halving keeps the scheme's order across crossings.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use util, only: twopi, sqrt2
    use parmot_mod, only: ro0
    use field_can_mod, only: field_can_t, eval_field => evaluate, integ_to_ref, &
                             ref_to_integ
    use field_can_spectre, only: set_spectre_volume_lock
    use orbit_symplectic_base, only: symplectic_integrator_t
    use orbit_symplectic, only: orbit_timestep_sympl, orbit_sympl_init
    use interface_crossing, only: apply_crossing, crossing_info_t, &
                                  crossing_log_record, CROSS_LOSS, CROSS_STOP

    implicit none
    private

    public :: sympl_spectre_state_t, sympl_spectre_reset, recanon_pphi, &
              orbit_microstep_sympl_spectre, sympl_landing_stats_reset, &
              sympl_landing_stats
    public :: SYMPL_SPECTRE_OK, SYMPL_SPECTRE_LOSS, SYMPL_SPECTRE_STOP

    integer, parameter :: SYMPL_SPECTRE_OK = 0
    integer, parameter :: SYMPL_SPECTRE_LOSS = 88
    integer, parameter :: SYMPL_SPECTRE_STOP = 89

    !> Landing target and acceptance for |rho_g - k| at the interface. The
    !> h*-solve aims at landing_tol; the per-volume splines are clamped at the
    !> volume edge, which kinks the landed rho_g just outside the interface, so a
    !> best-effort landing below landing_accept is still committed.
    real(dp), parameter :: landing_tol = 1.0d-10
    real(dp), parameter :: landing_accept = 1.0d-8
    integer, parameter :: max_locate = 64

    !> Events resolved within one microstep budget. A sheet-skimming class can
    !> chatter at an interface (the drift normal reverses across the sheet); past
    !> the cap the rest of the microstep is dropped (the RK45 path drops it after
    !> the first event), never a stop.
    integer, parameter :: max_events = 16
    !> Near-tangent crossings can put the landing root where the implicit solve
    !> of a long substep is multi-branched (Newton basin boundaries), making
    !> g(h) discontinuous. The substep is then halved and re-attempted: committed
    !> half-steps are fixed-length symplectic maps, and the landing solve re-runs
    !> in a smaller-h regime where the step map is single-valued.
    integer, parameter :: max_iters = 200
    real(dp), parameter :: h_min_frac = 1.0d-6
    real(dp), parameter :: budget_eps = 1.0d-12
    !> A committed substep moving rho_g by more than max_step_dr, or changing
    !> the energy 0.5*vpar^2 + mu*B by more than max_step_dh relative, is an
    !> unconverged Newton "teleport" (the steppers return silently after maxit):
    !> reject it and halve the substep. Both bounds are loose sanity checks that
    !> only O(1) garbage states trip; legitimate scheme error stays orders of
    !> magnitude below them at any usable step size.
    real(dp), parameter :: max_step_dr = 0.5d0
    real(dp), parameter :: max_step_dh = 0.1d0

    !> Home-side offset for the per-volume angle-transform evaluation at a landed
    !> point: the volume dispatch keys on int(rho_g), so exactly at rho_g = k it
    !> would pick the outer volume, whose zeta gauge differs by O(1).
    real(dp), parameter :: transform_bias = 1.0d-8

    type :: sympl_spectre_state_t
        real(dp) :: home_lo = 0.0_dp
        real(dp) :: home_hi = 1.0_dp
        integer :: mvol = 1
        integer :: mode = 0
        integer :: level = 1
        real(dp) :: dt_std = 0.0_dp
    end type sympl_spectre_state_t

    integer :: n_landings = 0
    integer :: n_stops = 0
    real(dp) :: max_landing_resid = 0.0_dp

contains

    subroutine sympl_spectre_reset(state, si, mvol, mode, level)
        type(sympl_spectre_state_t), intent(out) :: state
        type(symplectic_integrator_t), intent(in) :: si
        integer, intent(in) :: mvol, mode, level

        if (mode <= 0) then
            error stop 'spectre_sympl_orbit: crossing pipeline requires a '// &
                'symplectic scheme (integmode > 0)'
        end if
        if (si%ntau /= 1) then
            ! Re-canonicalization restarts the scheme from a single point, so no
            ! multistep history may carry over; all schemes here are one-step and
            ! must be driven one step per call.
            error stop 'spectre_sympl_orbit: requires one-step driving (ntau = 1)'
        end if

        state%mvol = mvol
        state%mode = mode
        state%level = level
        state%dt_std = si%dt
        call set_home(state, si%z(1))
    end subroutine sympl_spectre_reset

    subroutine set_home(state, rho)
        !> Fix the home volume [home_lo, home_hi] bracketing rho and pin the
        !> canonical field evaluation to it, so Newton iterates probing past an
        !> interior interface see the home volume's clamped splines instead of
        !> the neighbour's discontinuous field.
        type(sympl_spectre_state_t), intent(inout) :: state
        real(dp), intent(in) :: rho

        real(dp) :: lo

        lo = real(floor(rho), dp)
        lo = max(0.0_dp, min(lo, real(state%mvol - 1, dp)))
        state%home_lo = lo
        state%home_hi = lo + 1.0_dp
        call set_spectre_volume_lock(nint(state%home_hi))
    end subroutine set_home

    subroutine orbit_microstep_sympl_spectre(state, si, f, ipart, t_base_sec, &
                                             dt_sec, ierr, t_frac)
        !> Advance one microstep of length state%dt_std, resolving interface
        !> events by exact-landing substeps. On SYMPL_SPECTRE_LOSS or
        !> SYMPL_SPECTRE_STOP the orbit terminates at fraction t_frac of the
        !> microstep, with (si, f) holding the terminal (landed) state.
        type(sympl_spectre_state_t), intent(inout) :: state
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(in) :: ipart
        real(dp), intent(in) :: t_base_sec, dt_sec
        integer, intent(out) :: ierr
        real(dp), intent(out) :: t_frac

        type(symplectic_integrator_t) :: si0
        type(field_can_t) :: f0
        type(crossing_info_t) :: info
        real(dp) :: budget, used, h_try, h_land, resid, y_iface(5), y_out(5)
        integer :: iters, nev, iface, direction, ierr_step
        logical :: boundary, solved

        ierr = SYMPL_SPECTRE_OK
        t_frac = 1.0_dp
        budget = state%dt_std
        used = 0.0_dp
        h_try = budget
        nev = 0

        do iters = 1, max_iters
            si0 = si
            f0 = f
            si%dt = h_try
            call orbit_timestep_sympl(si, f, ierr_step)

            if (ierr_step == 0) then
                if (step_teleported(si0, f0, si, f)) then
                    si = si0
                    f = f0
                    h_try = 0.5_dp*h_try
                    if (h_try >= h_min_frac*state%dt_std) cycle
                    call record_stop(si, f, ipart, &
                                     t_base_sec + dt_sec*used/state%dt_std, &
                                     nint(state%home_hi), 1)
                    ierr = SYMPL_SPECTRE_STOP
                    t_frac = used/state%dt_std
                    exit
                end if
            end if

            call classify_step(state, si, ierr_step, boundary, iface, direction)

            if (.not. boundary) then
                if (ierr_step /= 0) then
                    ! Newton divergence away from the outer boundary: retry with
                    ! a halved substep, like the teleport and stalled-landing
                    ! branches.
                    si = si0
                    f = f0
                    h_try = 0.5_dp*h_try
                    if (h_try >= h_min_frac*state%dt_std) cycle
                    call record_stop(si, f, ipart, &
                                     t_base_sec + dt_sec*used/state%dt_std, &
                                     nint(state%home_hi), 1)
                    ierr = SYMPL_SPECTRE_STOP
                    t_frac = used/state%dt_std
                    exit
                end if
                used = used + h_try
                budget = state%dt_std - used
                if (budget <= budget_eps*state%dt_std) exit
                ! Grow the working substep gently after halvings: resetting to
                ! the full remaining budget would re-trigger the same failed
                ! solve indefinitely.
                h_try = min(2.0_dp*h_try, budget)
                cycle
            end if

            call locate_landing(si, f, si0, f0, h_try, real(iface, dp), &
                                direction, ierr_step == 0, h_land, resid, solved)
            if (.not. solved) then
                ! The landing root sits where the implicit solve of this substep
                ! length is multi-branched: halve and re-attempt from the same
                ! pre-step state.
                si = si0
                f = f0
                h_try = 0.5_dp*h_try
                if (h_try >= h_min_frac*state%dt_std) cycle
                call record_stop(si, f, ipart, &
                                 t_base_sec + dt_sec*used/state%dt_std, &
                                 iface, direction)
                ierr = SYMPL_SPECTRE_STOP
                t_frac = used/state%dt_std
                exit
            end if
            used = used + h_land
            budget = state%dt_std - used
            call update_landing_stats(resid)

            call landed_state(si, f, real(iface, dp), direction, y_iface)
            call apply_crossing(y_iface, iface, direction, state%mvol, &
                                state%level, y_out, info)
            call crossing_log_record(ipart, t_base_sec + dt_sec*used/state%dt_std, &
                                     info)

            if (info%event_type == CROSS_LOSS) then
                ierr = SYMPL_SPECTRE_LOSS
                t_frac = used/state%dt_std
                exit
            end if

            ! Select the target volume (and its field lock) from the nudged
            ! rho_g before re-canonicalizing in that volume's gauge.
            call set_home(state, y_out(1))
            call recanonicalize(state, si, f, y_out)

            nev = nev + 1
            if (nev >= max_events) exit
            if (budget <= budget_eps*state%dt_std) exit
            h_try = min(2.0_dp*h_try, budget)
        end do

        si%dt = state%dt_std
    end subroutine orbit_microstep_sympl_spectre

    subroutine classify_step(state, si, ierr_step, boundary, iface, direction)
        type(sympl_spectre_state_t), intent(in) :: state
        type(symplectic_integrator_t), intent(in) :: si
        integer, intent(in) :: ierr_step
        logical, intent(out) :: boundary
        integer, intent(out) :: iface, direction

        boundary = .false.
        iface = 0
        direction = 0

        if (ierr_step /= 0) then
            ! The stepper refuses to commit r > sympl_rmax (= Mvol). From the
            ! outermost volume that is the outward-crossing signal, with si%z
            ! still holding the pre-step state.
            if (nint(state%home_hi) == state%mvol) then
                boundary = .true.
                iface = state%mvol
                direction = 1
            end if
            return
        end if

        if (si%z(1) >= state%home_hi) then
            boundary = .true.
            iface = nint(state%home_hi)
            direction = 1
        else if (si%z(1) <= state%home_lo) then
            ! home_lo = 0 is the coordinate axis, not an interface: the chart
            ! flip (r, theta) -> (-r, theta + pi) inside the steppers handles
            ! axis flybys (#370), so only real interfaces raise events here.
            if (nint(state%home_lo) >= 1) then
                boundary = .true.
                iface = nint(state%home_lo)
                direction = -1
            end if
        end if
    end subroutine classify_step

    subroutine locate_landing(si, f, si0, f0, h_full, rho_k, direction, &
                              full_committed, h_land, resid, solved)
        !> Solve g(h) = dir*(z1(1; h) - rho_k) = 0 for the substep length by
        !> Illinois false position; every evaluation is one full implicit step of
        !> length h from the same pre-step state (si0, f0). Spline evaluations
        !> past the volume edge are clamped (field_can_spectre), which kinks g
        !> just outside the interface; the shrinking bracket walks through the
        !> kink, and trials the sympl_rmax guard refuses to commit (outermost
        !> interface only) tighten the upper bound by bisection. On return
        !> (si, f) hold the best landed step of length h_land.
        type(symplectic_integrator_t), intent(inout) :: si
        type(symplectic_integrator_t), intent(in) :: si0
        type(field_can_t), intent(inout) :: f
        type(field_can_t), intent(in) :: f0
        real(dp), intent(in) :: h_full, rho_k
        integer, intent(in) :: direction
        logical, intent(in) :: full_committed
        real(dp), intent(out) :: h_land, resid
        logical, intent(out) :: solved

        type(symplectic_integrator_t) :: si_best
        type(field_can_t) :: f_best
        real(dp) :: dir, t_lo, t_hi, f_lo, f_hi, t_mid, g, best
        integer :: it, ierr_step, last_side
        logical :: hi_valid

        dir = real(direction, dp)
        t_lo = 0.0_dp
        f_lo = dir*(si0%z(1) - rho_k)
        t_hi = h_full
        f_hi = 0.0_dp
        hi_valid = .false.
        last_side = 0

        best = huge(1.0_dp)
        h_land = h_full
        if (full_committed) then
            f_hi = dir*(si%z(1) - rho_k)
            hi_valid = .true.
            best = abs(f_hi)
            si_best = si
            f_best = f
        end if

        do it = 1, max_locate
            if (best < landing_tol) exit
            if (hi_valid) then
                if (abs(f_hi - f_lo) <= tiny(1.0_dp)) exit
                t_mid = t_hi - f_hi*(t_hi - t_lo)/(f_hi - f_lo)
                if (t_mid <= t_lo .or. t_mid >= t_hi) then
                    t_mid = 0.5_dp*(t_lo + t_hi)
                end if
            else
                t_mid = 0.5_dp*(t_lo + t_hi)
            end if

            si = si0
            f = f0
            si%dt = t_mid
            call orbit_timestep_sympl(si, f, ierr_step)
            g = dir*(si%z(1) - rho_k)
            if (ierr_step /= 0 .or. g /= g .or. &
                step_teleported(si0, f0, si, f)) then
                t_hi = t_mid
                hi_valid = .false.
                cycle
            end if

            if (abs(g) < best) then
                best = abs(g)
                h_land = t_mid
                si_best = si
                f_best = f
            end if

            if (g < 0.0_dp) then
                t_lo = t_mid
                f_lo = g
                if (last_side == 1) then
                    if (hi_valid) f_hi = 0.5_dp*f_hi
                end if
                last_side = 1
            else
                t_hi = t_mid
                f_hi = g
                hi_valid = .true.
                if (last_side == -1) f_lo = 0.5_dp*f_lo
                last_side = -1
            end if
        end do

        solved = best < landing_accept
        resid = best
        if (solved) then
            si = si_best
            f = f_best
        end if
    end subroutine locate_landing

    subroutine landed_state(si, f, rho_k, direction, y)
        !> Physical state y = (rho_g, theta, zeta, p, lambda) in SPECTRE reference
        !> coordinates at the landed integrator state. The zeta transform is read
        !> a hair inside the home volume (transform_bias) and unwrapped so the
        !> angle winding survives the modulo in integ_to_ref.
        type(symplectic_integrator_t), intent(in) :: si
        type(field_can_t), intent(in) :: f
        real(dp), intent(in) :: rho_k
        integer, intent(in) :: direction
        real(dp), intent(out) :: y(5)

        real(dp) :: xinteg(3), xref(3), p

        xinteg(1) = rho_k - real(direction, dp)*transform_bias
        xinteg(2:3) = si%z(2:3)
        call integ_to_ref(xinteg, xref)

        y(1) = si%z(1)
        y(2) = si%z(2)
        y(3) = unwrap_near(xref(3), si%z(3))
        p = sqrt(f%mu*f%Bmod + 0.5_dp*f%vpar**2)
        y(4) = p
        y(5) = f%vpar/(p*sqrt2)
    end subroutine landed_state

    subroutine recanonicalize(state, si, f, y)
        !> Re-initialize the integrator from the physical state y in the volume
        !> selected by its (nudged) rho_g, exactly like a fresh orbit start
        !> (init_sympl): per-volume gauges differ, so canonical momenta are
        !> rebuilt from the target volume's own Ath, Aph, hth, hph.
        type(sympl_spectre_state_t), intent(in) :: state
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: y(5)

        real(dp) :: xref(3), xinteg(3), z(4)

        xref = y(1:3)
        call ref_to_integ(xref, xinteg)
        xinteg(2) = y(2)
        xinteg(3) = unwrap_near(xinteg(3), y(3))

        call eval_field(f, xinteg(1), xinteg(2), xinteg(3), 0)
        si%pabs = y(4)
        f%mu = y(4)**2*(1.0_dp - y(5)**2)/f%Bmod
        f%ro0 = ro0/sqrt2
        f%vpar = y(4)*y(5)*sqrt2
        z(1:3) = xinteg
        z(4) = recanon_pphi(f, y(4), y(5))
        call orbit_sympl_init(si, f, z, state%dt_std, 1, si%rtol, state%mode)
    end subroutine recanonicalize

    pure function recanon_pphi(f, p, lambda) result(pphi)
        !> Canonical toroidal momentum of the physical state (p, lambda) in the
        !> gauge of the volume f is evaluated in. Correctness of re-canonicalization
        !> is that f carries the TARGET volume's Aph/hph: decoding pphi via get_val
        !> in the same volume returns vpar = p*lambda*sqrt2. A neighbour volume's
        !> gauge shifts the decoded vpar by [[Aph]]/(f%ro0*f%hph).
        type(field_can_t), intent(in) :: f
        real(dp), intent(in) :: p, lambda
        real(dp) :: pphi

        pphi = p*lambda*sqrt2*f%hph + f%Aph/f%ro0
    end function recanon_pphi

    pure function step_teleported(si0, f0, si, f) result(bad)
        type(symplectic_integrator_t), intent(in) :: si0, si
        type(field_can_t), intent(in) :: f0, f
        logical :: bad

        real(dp) :: h0, h1

        bad = abs(si%z(1) - si0%z(1)) > max_step_dr
        if (bad) return
        ! The steppers commit a negative-radius Newton root as an axis chart
        ! flip (r, theta) -> (|r|, theta + pi). Away from the axis a negative
        ! root is a Newton divergence, not an axis crossing: a committed theta
        ! jump of order pi at large radius is a teleport.
        if (abs(si%z(2) - si0%z(2)) > 1.0_dp) then
            bad = si%z(1) > 0.5_dp
            if (bad) return
        end if
        h0 = 0.5_dp*f0%vpar**2 + f0%mu*f0%Bmod
        h1 = 0.5_dp*f%vpar**2 + f%mu*f%Bmod
        bad = abs(h1 - h0) > max_step_dh*abs(h0)
    end function step_teleported

    pure function unwrap_near(a, near) result(b)
        real(dp), intent(in) :: a, near
        real(dp) :: b

        b = a + twopi*nint((near - a)/twopi)
    end function unwrap_near

    subroutine record_stop(si, f, ipart, t_sec, iface, direction)
        type(symplectic_integrator_t), intent(in) :: si
        type(field_can_t), intent(in) :: f
        integer, intent(in) :: ipart
        real(dp), intent(in) :: t_sec
        integer, intent(in) :: iface, direction

        type(crossing_info_t) :: info
        real(dp) :: vpar

        vpar = f%vpar/sqrt2
        info%event_type = CROSS_STOP
        info%iface = iface
        info%vol_from = merge(iface, iface + 1, direction == 1)
        info%vol_to = info%vol_from
        info%theta = si%z(2)
        info%zeta = si%z(3)
        info%vpar_before = vpar
        info%vpar_after = vpar
        info%mu = 0.5_dp*f%mu
        info%bmod_home = f%Bmod
        info%bmod_target = f%Bmod
        call crossing_log_record(ipart, t_sec, info)

        !$omp critical (spectre_sympl_landing)
        n_stops = n_stops + 1
        !$omp end critical (spectre_sympl_landing)
    end subroutine record_stop

    subroutine update_landing_stats(resid)
        real(dp), intent(in) :: resid

        !$omp critical (spectre_sympl_landing)
        n_landings = n_landings + 1
        max_landing_resid = max(max_landing_resid, resid)
        !$omp end critical (spectre_sympl_landing)
    end subroutine update_landing_stats

    subroutine sympl_landing_stats_reset
        !$omp critical (spectre_sympl_landing)
        n_landings = 0
        n_stops = 0
        max_landing_resid = 0.0_dp
        !$omp end critical (spectre_sympl_landing)
    end subroutine sympl_landing_stats_reset

    subroutine sympl_landing_stats(landings, max_resid, stops)
        integer, intent(out) :: landings, stops
        real(dp), intent(out) :: max_resid

        !$omp critical (spectre_sympl_landing)
        landings = n_landings
        max_resid = max_landing_resid
        stops = n_stops
        !$omp end critical (spectre_sympl_landing)
    end subroutine sympl_landing_stats

end module spectre_sympl_orbit
