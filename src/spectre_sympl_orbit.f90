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
        crossing_log_record, CROSS_LOSS, CROSS_STOP, &
        CROSS_REFLECTION, CROSS_SHEET
    use spectre_sheet_gc, only: sheet_gc_state_t, sheet_gc_initialize, &
        sheet_gc_advance, sheet_gc_to_y, SHEET_GC_OK
    use spectre_fo_hybrid, only: spectre_fo_state_t, spectre_fo_enter, &
        spectre_fo_advance_until_exit, SPECTRE_FO_OK, SPECTRE_FO_LOSS

    implicit none
    private

    public :: sympl_spectre_state_t, sympl_spectre_reset, recanon_pphi, &
        orbit_microstep_sympl_spectre, sympl_landing_stats_reset, &
        sympl_landing_stats, sympl_sheet_stats
    public :: sympl_fo_stats
    public :: SYMPL_SPECTRE_OK, SYMPL_SPECTRE_LOSS, SYMPL_SPECTRE_STOP, &
        SYMPL_SPECTRE_SKIM

    integer, parameter :: SYMPL_SPECTRE_OK = 0
    integer, parameter :: SYMPL_SPECTRE_LOSS = 88
    integer, parameter :: SYMPL_SPECTRE_STOP = 89
    !> A trapped marker whose magnetic-mirror turning point sits on an interface
    !> reflects, remains owned by its home volume, and is returned by its drift within one
    !> microstep -- with v_par ~ 0 the reflection cannot move it away. Left alone
    !> it reflects every microstep forever with no progress in (theta, zeta),
    !> hijacking the run. Such a marker is mirror-confined (it cannot cross the
    !> interface outward), so after SKIM_MAX consecutive same-interface
    !> reflections with no angular progress it is terminated as confined-skimming.
    integer, parameter :: SYMPL_SPECTRE_SKIM = 90
    integer, parameter :: SHEET_STAT_ENTRY = 1, SHEET_STAT_EXIT = 2
    integer, parameter :: SHEET_STAT_INIT_FAIL = 3, SHEET_STAT_ADVANCE_FAIL = 4
    integer, parameter :: STOP_SHEET = 1, STOP_TELEPORT = 2, STOP_STEP = 3
    integer, parameter :: STOP_LANDING = 4, STOP_EVENT_CAP = 5
    integer, parameter :: SKIM_MAX = 256
    real(dp), parameter :: SKIM_PROGRESS = 1.0d-4
    real(dp), parameter :: SKIM_PITCH_MAX = 2.0d-2

    !> Landing target and acceptance for |rho_g - k| at the interface. The
    !> h*-solve aims at landing_tol; the per-volume splines are clamped at the
    !> volume edge, which kinks the landed rho_g just outside the interface, so a
    !> best-effort landing below landing_accept is still committed.
    real(dp), parameter :: landing_tol = 1.0d-10
    real(dp), parameter :: landing_accept = 1.0d-8
    integer, parameter :: max_locate = 64

    !> Events resolved within one microstep budget. Reaching the cap is a reported
    !> stop: silently dropping the remaining time changes the fixed-step map.
    integer, parameter :: max_events = 64
    !> Near-tangent crossings can put the landing root where the implicit solve
    !> of a long substep is multi-branched (Newton basin boundaries), making
    !> g(h) discontinuous. The substep is then halved and re-attempted: committed
    !> half-steps are fixed-length symplectic maps, and the landing solve re-runs
    !> in a smaller-h regime where the step map is single-valued.
    integer, parameter :: max_iters = 200
    real(dp), parameter :: h_min_frac = 1.0d-8
    real(dp), parameter :: budget_eps = 1.0d-12
    !> A committed substep moving rho_g by more than max_step_dr, or changing
    !> the energy 0.5*vpar^2 + mu*B by more than max_step_dh relative, is an
    !> unconverged Newton "teleport" (the steppers return silently after maxit):
    !> reject it and halve the substep. Both bounds are loose sanity checks that
    !> only O(1) garbage states trip; legitimate scheme error stays orders of
    !> magnitude below them at any usable step size.
    real(dp), parameter :: max_step_dr = 0.5d0
    real(dp), parameter :: max_step_dh = 5.0d-5

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
        !> Consecutive same-interface reflections with no angular progress, and
        !> the interface and (theta, zeta) at the start of the current run.
        integer :: skim_count = 0
        integer :: skim_iface = -1
        real(dp) :: skim_theta = 0.0_dp
        real(dp) :: skim_zeta = 0.0_dp
        type(sheet_gc_state_t) :: sheet
        type(spectre_fo_state_t) :: fo
    end type sympl_spectre_state_t

    integer :: n_landings = 0
    integer :: n_stops = 0
    integer :: n_sheet_entries = 0
    integer :: n_sheet_exits = 0
    integer :: n_sheet_init_failures = 0
    integer :: n_sheet_advance_failures = 0
    integer :: n_sheet_failure_status(5) = 0
    integer :: n_stop_reason(5) = 0
    integer :: n_fo_entries = 0
    integer :: n_fo_exits = 0
    integer :: n_fo_losses = 0
    integer :: n_fo_failures = 0
    integer :: n_fo_failure_status(5) = 0
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

    subroutine set_home_volume(state, lvol)
        type(sympl_spectre_state_t), intent(inout) :: state
        integer, intent(in) :: lvol

        state%home_lo = real(lvol - 1, dp)
        state%home_hi = real(lvol, dp)
        call set_spectre_volume_lock(lvol)
    end subroutine set_home_volume

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
        real(dp) :: sheet_dt_used
        integer :: iters, nev, iface, direction, ierr_step, ierr_sheet, ierr_fo
        integer :: fo_owner
        integer :: prev_iface
        logical :: boundary, solved, sheet_exited, fo_exited

        ierr = SYMPL_SPECTRE_OK
        t_frac = 1.0_dp
        budget = state%dt_std
        used = 0.0_dp
        h_try = budget
        nev = 0
        prev_iface = -1

        if (state%fo%active) then
            call continue_fo(state, si, f, budget, used, y_out, fo_owner, &
                fo_exited, ierr_fo)
            if (ierr_fo /= SPECTRE_FO_OK) then
                call finish_fo_error(ierr_fo, si, f, ipart, t_base_sec + &
                    dt_sec*used/state%dt_std, state%fo%iface, 1, STOP_SHEET, &
                    used, state%dt_std, ierr, t_frac)
                si%dt = state%dt_std
                return
            end if
            if (.not. fo_exited .or. budget <= budget_eps*state%dt_std) then
                si%dt = state%dt_std
                return
            end if
            h_try = budget
        end if

        if (state%sheet%active) then
            call sheet_gc_advance(state%sheet, sqrt2*budget, ro0, sheet_exited, &
                sheet_dt_used, ierr_sheet)
            if (sheet_exited) call count_sheet_stat(SHEET_STAT_EXIT)
            if (ierr_sheet /= SHEET_GC_OK) &
                call count_sheet_stat(SHEET_STAT_ADVANCE_FAIL, ierr_sheet)
            used = sheet_dt_used/sqrt2
            budget = state%dt_std - used
            call sheet_gc_to_y(state%sheet, y_out)
            call set_home_volume(state, state%sheet%owner)
            call recanonicalize(state, si, f, y_out)
            if (ierr_sheet /= SHEET_GC_OK) then
                y_iface = y_out
                call start_fo(state, si, f, y_iface, state%sheet%owner, &
                    state%sheet%iface, budget, used, y_out, fo_owner, fo_exited, &
                    ierr_fo)
                if (ierr_fo /= SPECTRE_FO_OK) then
                    direction = merge(1, -1, &
                        state%sheet%owner == state%sheet%iface)
                    call finish_fo_error(ierr_fo, si, f, ipart, &
                        t_base_sec + dt_sec*used/state%dt_std, &
                        state%sheet%iface, direction, STOP_SHEET, used, &
                        state%dt_std, ierr, t_frac)
                    si%dt = state%dt_std
                    return
                end if
                if (.not. fo_exited .or. &
                    budget <= budget_eps*state%dt_std) then
                    si%dt = state%dt_std
                    return
                end if
                sheet_exited = .true.
            end if
            if (.not. sheet_exited .or. budget <= budget_eps*state%dt_std) then
                si%dt = state%dt_std
                return
            end if
            h_try = budget
        end if

        iters = 0
        do
            iters = iters + 1
            if (iters > max_iters) then
                call integrator_state(si, f, y_iface)
                iface = max(1, min(state%mvol - 1, nint(si%z(1))))
                call start_fo(state, si, f, y_iface, nint(state%home_hi), &
                    iface, budget, used, y_out, fo_owner, fo_exited, ierr_fo)
                if (ierr_fo /= SPECTRE_FO_OK) then
                    call finish_fo_error(ierr_fo, si, f, ipart, t_base_sec + &
                        dt_sec*used/state%dt_std, iface, 1, STOP_STEP, used, &
                        state%dt_std, ierr, t_frac)
                    exit
                end if
                if (.not. fo_exited .or. &
                    budget <= budget_eps*state%dt_std) exit
                iters = 0
                h_try = budget
                cycle
            end if
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
                    call integrator_state(si, f, y_iface)
                    iface = max(1, min(state%mvol - 1, nint(si%z(1))))
                    call start_fo(state, si, f, y_iface, nint(state%home_hi), &
                        iface, budget, used, y_out, fo_owner, fo_exited, ierr_fo)
                    if (ierr_fo /= SPECTRE_FO_OK) then
                        call finish_fo_error(ierr_fo, si, f, ipart, &
                            t_base_sec + dt_sec*used/state%dt_std, iface, 1, &
                            STOP_TELEPORT, used, state%dt_std, ierr, t_frac)
                        exit
                    end if
                    if (.not. fo_exited .or. &
                        budget <= budget_eps*state%dt_std) exit
                    h_try = budget
                    cycle
                end if
            end if

            call classify_step(state, si0%z(1), si, ierr_step, boundary, iface, &
                direction)

            if (.not. boundary) then
                if (ierr_step /= 0) then
                    ! Newton divergence away from the outer boundary: retry with
                    ! a halved substep, like the teleport and stalled-landing
                    ! branches.
                    si = si0
                    f = f0
                    h_try = 0.5_dp*h_try
                    if (h_try >= h_min_frac*state%dt_std) cycle
                    call integrator_state(si, f, y_iface)
                    iface = max(1, min(state%mvol - 1, nint(si%z(1))))
                    call start_fo(state, si, f, y_iface, nint(state%home_hi), &
                        iface, budget, used, y_out, fo_owner, fo_exited, ierr_fo)
                    if (ierr_fo /= SPECTRE_FO_OK) then
                        call finish_fo_error(ierr_fo, si, f, ipart, &
                            t_base_sec + dt_sec*used/state%dt_std, iface, 1, &
                            STOP_STEP, used, state%dt_std, ierr, t_frac)
                        exit
                    end if
                    if (.not. fo_exited .or. &
                        budget <= budget_eps*state%dt_std) exit
                    h_try = budget
                    cycle
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
                call integrator_state(si, f, y_iface)
                call start_fo(state, si, f, y_iface, nint(state%home_hi), &
                    iface, budget, used, y_out, fo_owner, fo_exited, ierr_fo)
                if (ierr_fo /= SPECTRE_FO_OK) then
                    call finish_fo_error(ierr_fo, si, f, ipart, t_base_sec + &
                        dt_sec*used/state%dt_std, iface, direction, STOP_LANDING, &
                        used, state%dt_std, ierr, t_frac)
                    exit
                end if
                if (.not. fo_exited .or. &
                    budget <= budget_eps*state%dt_std) exit
                h_try = budget
                cycle
            end if
            used = used + h_land
            budget = state%dt_std - used
            call update_landing_stats(resid)

            call landed_state(si, f, real(iface, dp), direction, y_iface)
            call apply_crossing(y_iface, iface, direction, state%mvol, &
                state%level, y_out, info)
            if (info%event_type == CROSS_LOSS) then
                call crossing_log_record(ipart, &
                    t_base_sec + dt_sec*used/state%dt_std, info)
                ierr = SYMPL_SPECTRE_LOSS
                t_frac = used/state%dt_std
                exit
            end if

            if (iface == prev_iface) then
                call sheet_gc_initialize(iface, direction, y_iface, info%mu, &
                    info%bmod_home, ro0, state%sheet, ierr_sheet)
                if (ierr_sheet == SHEET_GC_OK) then
                    call count_sheet_stat(SHEET_STAT_ENTRY)
                    info%event_type = CROSS_SHEET
                    info%vol_to = info%vol_from
                    call crossing_log_record(ipart, &
                        t_base_sec + dt_sec*used/state%dt_std, info)
                    call sheet_gc_advance(state%sheet, sqrt2*budget, ro0, &
                        sheet_exited, sheet_dt_used, ierr_sheet)
                    if (sheet_exited) call count_sheet_stat(SHEET_STAT_EXIT)
                    if (ierr_sheet /= SHEET_GC_OK) &
                        call count_sheet_stat(SHEET_STAT_ADVANCE_FAIL, ierr_sheet)
                    used = used + sheet_dt_used/sqrt2
                    budget = state%dt_std - used
                    call sheet_gc_to_y(state%sheet, y_out)
                    call set_home_volume(state, state%sheet%owner)
                    call recanonicalize(state, si, f, y_out)
                    if (ierr_sheet /= SHEET_GC_OK) then
                        y_iface = y_out
                        call start_fo(state, si, f, y_iface, state%sheet%owner, &
                            iface, budget, used, y_out, fo_owner, fo_exited, ierr_fo)
                        if (ierr_fo /= SPECTRE_FO_OK) then
                            call finish_fo_error(ierr_fo, si, f, ipart, &
                                t_base_sec + dt_sec*used/state%dt_std, iface, &
                                direction, STOP_SHEET, used, state%dt_std, ierr, &
                                t_frac)
                            exit
                        end if
                        if (.not. fo_exited .or. &
                            budget <= budget_eps*state%dt_std) exit
                        sheet_exited = .true.
                    end if
                    if (.not. sheet_exited .or. &
                        budget <= budget_eps*state%dt_std) exit
                    h_try = budget
                    prev_iface = -1
                    cycle
                end if
                call count_sheet_stat(SHEET_STAT_INIT_FAIL, ierr_sheet)
                call start_fo(state, si, f, y_iface, info%vol_from, iface, &
                    budget, used, y_out, fo_owner, fo_exited, ierr_fo)
                if (ierr_fo /= SPECTRE_FO_OK) then
                    call finish_fo_error(ierr_fo, si, f, ipart, t_base_sec + &
                        dt_sec*used/state%dt_std, iface, direction, STOP_SHEET, &
                        used, state%dt_std, ierr, t_frac)
                    exit
                end if
                if (.not. fo_exited .or. &
                    budget <= budget_eps*state%dt_std) exit
                h_try = budget
                prev_iface = -1
                cycle
            end if

            call crossing_log_record(ipart, t_base_sec + dt_sec*used/state%dt_std, &
                info)

            if (skimming(state, info)) then
                ierr = SYMPL_SPECTRE_SKIM
                t_frac = used/state%dt_std
                exit
            end if

            call set_home_volume(state, info%vol_to)
            call recanonicalize(state, si, f, y_out)

            prev_iface = iface
            nev = nev + 1
            if (nev >= max_events) then
                y_iface = y_out
                call start_fo(state, si, f, y_iface, info%vol_to, iface, &
                    budget, used, y_out, fo_owner, fo_exited, ierr_fo)
                if (ierr_fo /= SPECTRE_FO_OK) then
                    call finish_fo_error(ierr_fo, si, f, ipart, &
                        t_base_sec + dt_sec*used/state%dt_std, iface, direction, &
                        STOP_EVENT_CAP, used, state%dt_std, ierr, t_frac)
                    exit
                end if
                if (.not. fo_exited .or. &
                    budget <= budget_eps*state%dt_std) exit
                nev = 0
                prev_iface = -1
                h_try = budget
                cycle
            end if
            if (budget <= budget_eps*state%dt_std) exit
            h_try = min(2.0_dp*h_try, budget)
        end do

        si%dt = state%dt_std
    end subroutine orbit_microstep_sympl_spectre

    logical function skimming(state, info)
        !> True once a marker has reflected SKIM_MAX times in a row off the same
        !> interface without moving in (theta, zeta): a mirror turning point
        !> pinned on the interface by its drift. Any non-reflection event, a
        !> different interface, or angular progress resets the run counter, so a
        !> physically precessing trapped particle (which advances in theta/zeta
        !> between bounces) is never flagged.
        type(sympl_spectre_state_t), intent(inout) :: state
        type(crossing_info_t), intent(in) :: info

        real(dp) :: dth, dze, speed

        skimming = .false.

        if (info%event_type /= CROSS_REFLECTION) then
            state%skim_count = 0
            state%skim_iface = -1
            return
        end if
        speed = sqrt(info%vpar_before**2 + 2.0_dp*info%mu*info%bmod_home)
        if (abs(info%vpar_before) > SKIM_PITCH_MAX*speed) then
            state%skim_count = 0
            state%skim_iface = -1
            return
        end if

        dth = abs(info%theta - state%skim_theta)
        dze = abs(info%zeta - state%skim_zeta)
        if (info%iface == state%skim_iface .and. dth < SKIM_PROGRESS &
            .and. dze < SKIM_PROGRESS) then
            state%skim_count = state%skim_count + 1
        else
            state%skim_count = 1
            state%skim_iface = info%iface
            state%skim_theta = info%theta
            state%skim_zeta = info%zeta
        end if

        skimming = state%skim_count >= SKIM_MAX
    end function skimming

    subroutine classify_step(state, rho_start, si, ierr_step, boundary, iface, &
            direction)
        type(sympl_spectre_state_t), intent(in) :: state
        real(dp), intent(in) :: rho_start
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

        if (si%z(1) > state%home_hi .or. &
            (si%z(1) == state%home_hi .and. rho_start < state%home_hi)) then
            boundary = .true.
            iface = nint(state%home_hi)
            direction = 1
        else if (si%z(1) < state%home_lo .or. &
                (si%z(1) == state%home_lo .and. rho_start > state%home_lo)) then
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
        p = si%pabs
        y(4) = p
        y(5) = shell_pitch(f, p)
    end subroutine landed_state

    subroutine integrator_state(si, f, y)
        type(symplectic_integrator_t), intent(in) :: si
        type(field_can_t), intent(in) :: f
        real(dp), intent(out) :: y(5)

        real(dp) :: xref(3), p

        call integ_to_ref(si%z(1:3), xref)
        y(1) = si%z(1)
        y(2) = si%z(2)
        y(3) = unwrap_near(xref(3), si%z(3))
        p = si%pabs
        y(4) = p
        y(5) = shell_pitch(f, p)
    end subroutine integrator_state

    pure function shell_pitch(f, p) result(lambda)
        type(field_can_t), intent(in) :: f
        real(dp), intent(in) :: p
        real(dp) :: lambda, lambda2

        lambda2 = max(1.0_dp - f%mu*f%Bmod/p**2, 0.0_dp)
        lambda = sign(sqrt(lambda2), f%vpar)
    end function shell_pitch

    subroutine start_fo(state, si, f, y, owner, iface, budget, used, y_out, &
            owner_out, exited, ierr)
        type(sympl_spectre_state_t), intent(inout) :: state
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: y(5)
        integer, intent(in) :: owner, iface
        real(dp), intent(inout) :: budget, used
        real(dp), intent(out) :: y_out(5)
        integer, intent(out) :: owner_out, ierr
        logical, intent(out) :: exited

        state%sheet%active = .false.
        call spectre_fo_enter(state%fo, y, owner, iface, ro0/sqrt2, ierr)
        if (ierr /= SPECTRE_FO_OK) then
            call count_fo_failure(ierr)
            return
        end if
        call count_fo_entry
        call continue_fo(state, si, f, budget, used, y_out, owner_out, &
            exited, ierr)
    end subroutine start_fo

    subroutine continue_fo(state, si, f, budget, used, y_out, owner_out, &
            exited, ierr)
        type(sympl_spectre_state_t), intent(inout) :: state
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(inout) :: budget, used
        real(dp), intent(out) :: y_out(5)
        integer, intent(out) :: owner_out, ierr
        logical, intent(out) :: exited

        real(dp) :: dt_used

        call spectre_fo_advance_until_exit(state%fo, budget, ro0/sqrt2, &
            dt_used, y_out, owner_out, exited, ierr)
        used = used + dt_used
        budget = state%dt_std - used
        if (ierr /= SPECTRE_FO_OK) then
            if (ierr == SPECTRE_FO_LOSS) then
                !$omp critical (spectre_sympl_landing)
                n_fo_losses = n_fo_losses + 1
                !$omp end critical (spectre_sympl_landing)
            else
                call count_fo_failure(ierr)
            end if
            return
        end if
        if (.not. exited) return
        call count_fo_exit
        call set_home_volume(state, owner_out)
        call recanonicalize(state, si, f, y_out)
    end subroutine continue_fo

    subroutine recanonicalize(state, si, f, y)
        !> Re-initialize the integrator from the physical state y in the explicitly
        !> owned volume, exactly like a fresh orbit start
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

    subroutine record_stop(si, f, ipart, t_sec, iface, direction, reason)
        type(symplectic_integrator_t), intent(in) :: si
        type(field_can_t), intent(in) :: f
        integer, intent(in) :: ipart
        real(dp), intent(in) :: t_sec
        integer, intent(in) :: iface, direction, reason

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
        if (reason >= 1 .and. reason <= size(n_stop_reason)) &
            n_stop_reason(reason) = n_stop_reason(reason) + 1
        !$omp end critical (spectre_sympl_landing)
    end subroutine record_stop

    subroutine finish_fo_error(fo_status, si, f, ipart, t_sec, iface, &
            direction, reason, used, dt, ierr, t_frac)
        integer, intent(in) :: fo_status, ipart, iface, direction, reason
        type(symplectic_integrator_t), intent(in) :: si
        type(field_can_t), intent(in) :: f
        real(dp), intent(in) :: t_sec, used, dt
        integer, intent(out) :: ierr
        real(dp), intent(out) :: t_frac

        if (fo_status == SPECTRE_FO_LOSS) then
            ierr = SYMPL_SPECTRE_LOSS
        else
            call record_stop(si, f, ipart, t_sec, iface, direction, reason)
            ierr = SYMPL_SPECTRE_STOP
        end if
        t_frac = used/dt
    end subroutine finish_fo_error

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
        n_sheet_entries = 0
        n_sheet_exits = 0
        n_sheet_init_failures = 0
        n_sheet_advance_failures = 0
        n_sheet_failure_status = 0
        n_stop_reason = 0
        n_fo_entries = 0
        n_fo_exits = 0
        n_fo_losses = 0
        n_fo_failures = 0
        n_fo_failure_status = 0
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

    subroutine count_sheet_stat(kind, status)
        integer, intent(in) :: kind
        integer, intent(in), optional :: status

        !$omp critical (spectre_sympl_landing)
        select case (kind)
        case (SHEET_STAT_ENTRY)
            n_sheet_entries = n_sheet_entries + 1
        case (SHEET_STAT_EXIT)
            n_sheet_exits = n_sheet_exits + 1
        case (SHEET_STAT_INIT_FAIL)
            n_sheet_init_failures = n_sheet_init_failures + 1
        case (SHEET_STAT_ADVANCE_FAIL)
            n_sheet_advance_failures = n_sheet_advance_failures + 1
        end select
        if (present(status)) then
            if (status >= 1 .and. status <= size(n_sheet_failure_status)) &
                n_sheet_failure_status(status) = n_sheet_failure_status(status) + 1
        end if
        !$omp end critical (spectre_sympl_landing)
    end subroutine count_sheet_stat

    subroutine sympl_sheet_stats(entries, exits, init_failures, advance_failures, &
            failure_status, stop_reason)
        integer, intent(out) :: entries, exits, init_failures, advance_failures
        integer, intent(out) :: failure_status(5)
        integer, intent(out) :: stop_reason(5)

        !$omp critical (spectre_sympl_landing)
        entries = n_sheet_entries
        exits = n_sheet_exits
        init_failures = n_sheet_init_failures
        advance_failures = n_sheet_advance_failures
        failure_status = n_sheet_failure_status
        stop_reason = n_stop_reason
        !$omp end critical (spectre_sympl_landing)
    end subroutine sympl_sheet_stats

    subroutine count_fo_entry
        !$omp critical (spectre_sympl_landing)
        n_fo_entries = n_fo_entries + 1
        !$omp end critical (spectre_sympl_landing)
    end subroutine count_fo_entry

    subroutine count_fo_exit
        !$omp critical (spectre_sympl_landing)
        n_fo_exits = n_fo_exits + 1
        !$omp end critical (spectre_sympl_landing)
    end subroutine count_fo_exit

    subroutine count_fo_failure(status)
        integer, intent(in) :: status

        !$omp critical (spectre_sympl_landing)
        n_fo_failures = n_fo_failures + 1
        if (status >= 1 .and. status <= size(n_fo_failure_status)) &
            n_fo_failure_status(status) = n_fo_failure_status(status) + 1
        !$omp end critical (spectre_sympl_landing)
    end subroutine count_fo_failure

    subroutine sympl_fo_stats(entries, exits, losses, failures, failure_status)
        integer, intent(out) :: entries, exits, losses, failures
        integer, intent(out) :: failure_status(5)

        !$omp critical (spectre_sympl_landing)
        entries = n_fo_entries
        exits = n_fo_exits
        losses = n_fo_losses
        failures = n_fo_failures
        failure_status = n_fo_failure_status
        !$omp end critical (spectre_sympl_landing)
    end subroutine sympl_fo_stats

end module spectre_sympl_orbit
