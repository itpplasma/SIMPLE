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
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    use util, only: twopi, sqrt2
    use parmot_mod, only: ro0
    use field_can_mod, only: field_can_t, eval_field => evaluate, integ_to_ref, &
        ref_to_integ, get_derivatives, get_val
    use field_can_spectre, only: set_spectre_volume_lock
    use orbit_symplectic_base, only: symplectic_integrator_t, MIDPOINT, &
        symplectic_newton_warning_mode, symplectic_spectre_newton_warning_factor, &
        SYMPLECTIC_STEP_OUTSIDE_DOMAIN
    use orbit_symplectic, only: orbit_timestep_sympl, orbit_sympl_init
    use interface_crossing, only: apply_crossing, crossing_info_t, &
        crossing_log_record, CROSS_LOSS, CROSS_STOP, CROSS_RECOVERY, CROSS_INVALID, &
        CROSS_REFLECTION, CROSS_SHEET, CROSS_CROSSING, CROSSING_SYMPLECTIC
    use spectre_refraction, only: midpoint_energy, interface_flow, REFRACTION_OK
    use spectre_sheet_gc, only: sheet_gc_state_t, sheet_gc_initialize, &
        sheet_gc_advance, sheet_gc_to_y, SHEET_GC_OK
    use spectre_fo_hybrid, only: spectre_fo_state_t, spectre_fo_enter, &
        spectre_fo_advance_until_exit, SPECTRE_FO_OK, SPECTRE_FO_LOSS
    use diag_counters, only: count_event, EVT_WARNING_STEP_SKIP

    implicit none
    private

    public :: sympl_spectre_state_t, sympl_spectre_reset, recanon_pphi, &
        orbit_microstep_sympl_spectre, sympl_landing_stats_reset, &
        sympl_landing_stats, sympl_landing_eval_stats, sympl_sheet_stats
    public :: sympl_fo_stats, sympl_symplectic_stats
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
    !> Loose state-sanity bounds for a committed substep. The Hamiltonian is not
    !> required to stay constant from step to step: symplectic discretizations
    !> normally oscillate around a nearby modified Hamiltonian. The shell bound
    !> only rejects a state that has left the tested bounded-error band around
    !> the fixed physical speed pabs. Interface maps re-pin this shell exactly.
    real(dp), parameter :: max_step_dr = 0.5d0
    real(dp), parameter :: max_relative_shell_defect = 5.0d-4

    !> Home-side offset for the per-volume angle-transform evaluation at a landed
    !> point: the volume dispatch keys on int(rho_g), so exactly at rho_g = k it
    !> would pick the outer volume, whose zeta gauge differs by O(1).
    real(dp), parameter :: transform_bias = 1.0d-8

    !> Symplectic crossing (crossing_level 2): landing residual, energy-root
    !> scan along the interface (first point and number of doublings, in
    !> radians of poloidal motion) and the secant tolerance on the discrete
    !> energy, relative to |H|.
    real(dp), parameter :: symplectic_landing_tol = 1.0d-13
    real(dp), parameter :: scan_tau0 = 1.0d-6
    integer, parameter :: n_scan = 24
    real(dp), parameter :: refine_tol = 1.0d-14
    integer, parameter :: max_refine = 40

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
        !> An exhausted warning-mode recovery holds one unresolved interval.
        !> The next microstep retries guiding-centre advancement from the last
        !> valid state. A later identical failure consumes one unresolved
        !> interval directly instead of repeating the full recovery cascade;
        !> default warning mode never turns that retry into a marker stop.
        logical :: warning_hold_latched = .false.
        integer :: warning_hold_count = 0
        integer :: warning_hold_reason = STOP_STEP
        integer :: warning_hold_iface = 0
        integer :: warning_hold_direction = 0
        real(dp) :: warning_hold_z(4) = 0.0_dp
        type(sheet_gc_state_t) :: sheet
        type(spectre_fo_state_t) :: fo
    end type sympl_spectre_state_t

    !> Hot-path statistics are per-thread. Each OpenMP worker writes only its
    !> own slot; totals are reduced after the parallel trace has joined.
    integer, allocatable :: n_landings(:), n_landing_evals(:), n_stops(:)
    integer, allocatable :: n_sheet_entries(:), n_sheet_exits(:)
    integer, allocatable :: n_sheet_init_failures(:), n_sheet_advance_failures(:)
    integer, allocatable :: n_sheet_failure_status(:, :), n_stop_reason(:, :)
    integer, allocatable :: n_fo_entries(:), n_fo_exits(:), n_fo_losses(:)
    integer, allocatable :: n_fo_failures(:), n_fo_failure_status(:, :)
    real(dp), allocatable :: max_landing_resid(:)
    integer, allocatable :: n_sympl_transmissions(:), n_sympl_reflections(:)
    integer, allocatable :: n_sympl_fallbacks(:)

contains

    subroutine sympl_spectre_reset(state, si, mvol, mode, level)
        type(sympl_spectre_state_t), intent(out) :: state
        type(symplectic_integrator_t), intent(in) :: si
        integer, intent(in) :: mvol, mode, level

        if (mode <= 0) then
            error stop 'spectre_sympl_orbit: crossing pipeline requires a '// &
                'symplectic scheme (integmode > 0)'
        end if
        if (level == CROSSING_SYMPLECTIC .and. mode /= MIDPOINT) then
            error stop 'spectre_sympl_orbit: crossing_level = 2 requires '// &
                'integmode = 3 (midpoint)'
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
        logical :: boundary, solved, sheet_exited, fo_exited, resume_gc, crossed

        ierr = SYMPL_SPECTRE_OK
        t_frac = 1.0_dp
        si%warning_factor = symplectic_spectre_newton_warning_factor
        budget = state%dt_std
        used = 0.0_dp
        h_try = budget
        nev = 0
        prev_iface = -1

        if (state%fo%active) then
            call continue_fo(state, si, f, budget, used, y_out, fo_owner, &
                fo_exited, ierr_fo)
            if (ierr_fo /= SPECTRE_FO_OK) then
                if (recover_fo_failure(state, si, f, ierr_fo, fo_owner, &
                    budget, used, ierr, t_frac)) then
                    si%dt = state%dt_std
                    return
                end if
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
                    if (recover_fo_failure(state, si, f, ierr_fo, fo_owner, &
                        budget, used, ierr, t_frac)) then
                        si%dt = state%dt_std
                        return
                    end if
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
                call recover_local_error(state, si, f, ipart, t_base_sec + &
                    dt_sec*used/state%dt_std, STOP_STEP, budget, used, ierr, &
                    t_frac, resume_gc)
                if (resume_gc) then
                    iters = 0
                    h_try = budget
                    cycle
                end if
                exit
            end if
            si0 = si
            f0 = f
            si%dt = h_try
            call orbit_timestep_sympl(si, f, ierr_step)

            if (ierr_step == 0) then
                if (step_teleported(si0, si, f)) then
                    si = si0
                    f = f0
                    if (same_warning_hold(state, si, STOP_TELEPORT, 0, 0)) then
                        if (recover_warning_skip(state, si, budget, used, ierr, &
                            t_frac, STOP_TELEPORT, 0, 0)) exit
                    end if
                    h_try = 0.5_dp*h_try
                    if (h_try >= h_min_frac*state%dt_std) cycle
                    call recover_local_error(state, si, f, ipart, t_base_sec + &
                        dt_sec*used/state%dt_std, STOP_TELEPORT, budget, used, &
                        ierr, t_frac, resume_gc)
                    if (resume_gc) then
                        h_try = budget
                        cycle
                    end if
                    exit
                end if
                state%warning_hold_latched = .false.
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
                    if (same_warning_hold(state, si, STOP_STEP, 0, 0)) then
                        if (recover_warning_skip(state, si, budget, used, ierr, &
                            t_frac, STOP_STEP, 0, 0)) exit
                    end if
                    h_try = 0.5_dp*h_try
                    if (h_try >= h_min_frac*state%dt_std) cycle
                    call recover_local_error(state, si, f, ipart, t_base_sec + &
                        dt_sec*used/state%dt_std, STOP_STEP, budget, used, ierr, &
                        t_frac, resume_gc)
                    if (resume_gc) then
                        h_try = budget
                        cycle
                    end if
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

            if (state%level == CROSSING_SYMPLECTIC .and. iface >= 1 .and. &
                iface < state%mvol) then
                call symplectic_crossing(state, si, f, si0, f0, h_try, &
                    ierr_step == 0, iface, direction, info, crossed)
                if (crossed) then
                    used = used + h_try
                    budget = state%dt_std - used
                    call crossing_log_record(ipart, &
                        t_base_sec + dt_sec*used/state%dt_std, info)
                    if (skimming(state, info)) then
                        ierr = SYMPL_SPECTRE_SKIM
                        t_frac = used/state%dt_std
                        exit
                    end if
                    if (budget <= budget_eps*state%dt_std) exit
                    h_try = budget
                    cycle
                end if
                call count_symplectic_crossing(CROSS_INVALID)
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
                    if (recover_fo_failure(state, si, f, ierr_fo, fo_owner, &
                        budget, used, ierr, t_frac)) exit
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
            call landed_state(si, f, real(iface, dp), direction, y_iface)
            call apply_crossing(y_iface, iface, direction, state%mvol, &
                state%level, y_out, info)
            if (info%event_type == CROSS_INVALID) then
                si = si0
                f = f0
                h_try = 0.5_dp*h_try
                if (h_try >= h_min_frac*state%dt_std) cycle
                call recover_local_error(state, si, f, ipart, t_base_sec + &
                    dt_sec*used/state%dt_std, STOP_LANDING, budget, used, ierr, &
                    t_frac, resume_gc)
                if (resume_gc) then
                    h_try = budget
                    cycle
                end if
                exit
            end if
            used = used + h_land
            budget = state%dt_std - used
            call update_landing_stats(resid)
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
                            if (recover_fo_failure(state, si, f, ierr_fo, &
                                fo_owner, budget, used, ierr, t_frac)) exit
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
                    if (recover_fo_failure(state, si, f, ierr_fo, fo_owner, &
                        budget, used, ierr, t_frac)) exit
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
                    if (recover_fo_failure(state, si, f, ierr_fo, fo_owner, &
                        budget, used, ierr, t_frac)) exit
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

    subroutine symplectic_crossing(state, si, f, si0, f0, h_step, full_committed, &
            iface, direction, info, done)
        !> Exactly symplectic crossing step of the midpoint scheme (crossing_level
        !> 2). The step of length h_step from (si0, f0) hits interface iface. It
        !> lands on the interface with a midpoint substep, flows along the
        !> interface vector field X_s and completes the step with a second
        !> midpoint substep, in the target volume (transmission) or back in the
        !> home volume (reflection). The flow length matches the discrete
        !> energies H(z_mid) of the two substeps (see spectre_refraction). The
        !> state crosses with its physical (position, v_par, mu); p_phi is
        !> rebuilt in the target gauge.
        !>
        !> The exit is the one of the limit of a smoothed jump: inside a layer
        !> H = (1 - chi) H_below + chi H_above the state drifts along X_s with
        !> rate (H_above - H_below) chi'/width, so the flow runs one way only,
        !> and the state leaves at the first point where H_target (transmission)
        !> or H_home (reflection) returns to the entry energy. On done = .false.
        !> (si, f) and the home volume are unchanged and the caller falls back to
        !> the Level-1 map.
        type(sympl_spectre_state_t), intent(inout) :: state
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        type(symplectic_integrator_t), intent(in) :: si0
        type(field_can_t), intent(in) :: f0
        real(dp), intent(in) :: h_step
        logical, intent(in) :: full_committed
        integer, intent(in) :: iface, direction
        type(crossing_info_t), intent(out) :: info
        logical, intent(out) :: done

        type(symplectic_integrator_t) :: si_land, si_out
        type(field_can_t) :: f_land, f_target, f_out, fw
        real(dp) :: h_land, resid, h_rem, e_land, xref(3), xinteg(3), zu(4)
        real(dp) :: tau, roots(2*n_scan), keys(2*n_scan), scale_home, scale_target
        integer :: home, target, status, sigma, n_home, n_target, n, i, j
        integer :: kind(2*n_scan)
        logical :: solved

        done = .false.
        home = nint(state%home_hi)
        target = merge(iface + 1, iface, direction == 1)

        si_land = si
        f_land = f
        call locate_landing(si_land, f_land, si0, f0, h_step, real(iface, dp), &
            direction, full_committed, h_land, resid, solved, &
            tol=symplectic_landing_tol)
        if (.not. solved) return
        h_rem = h_step - h_land

        fw = f_land
        call midpoint_energy(si0%z, f0%pth, si_land%z, f_land%pth, fw, e_land, &
            status)
        if (status /= REFRACTION_OK) return

        call integ_to_ref(si_land%z(1:3), xref)
        xref(3) = unwrap_near(xref(3), si_land%z(3))

        info%iface = iface
        info%vol_from = home
        info%theta = xref(2)
        info%zeta = xref(3)
        info%vpar_before = f_land%vpar/sqrt2
        info%mu = 0.5_dp*f_land%mu
        info%bmod_home = f_land%Bmod

        call set_home_volume(state, target)
        call ref_to_integ(xref, xinteg)
        zu(1) = si_land%z(1)
        zu(2) = unwrap_near(xinteg(2), si_land%z(2))
        zu(3) = unwrap_near(xinteg(3), si_land%z(3))
        f_target = f_land
        call eval_field(f_target, zu(1), zu(2), zu(3), 0)
        zu(4) = f_land%vpar*f_target%hph + f_target%Aph/f_target%ro0
        call get_val(f_target, zu(4))

        ! Flow direction: sign of H_above - H_below at the landing point, with
        ! X_s oriented by s = rho_g - iface.
        sigma = merge(1, -1, real(direction, dp)*(f_target%H - f_land%H) >= 0.0_dp)

        call flow_scale(zu, f_land, scale_target)
        call scan_energy_roots(zu, f_land, e_land, scale_target, sigma, &
            roots(1:n_scan), n_target)
        call set_home_volume(state, home)
        call flow_scale(si_land%z, f_land, scale_home)
        call scan_energy_roots(si_land%z, f_land, e_land, scale_home, sigma, &
            roots(n_target + 1:n_target + n_scan), n_home)
        n = n_target + n_home
        kind(1:n_target) = CROSS_CROSSING
        kind(n_target + 1:n) = CROSS_REFLECTION
        keys(1:n) = abs(roots(1:n))

        j = 0
        do i = 1, n
            j = minloc(keys(1:n), dim=1)
            keys(j) = huge(1.0_dp)
            if (kind(j) == CROSS_CROSSING) then
                call set_home_volume(state, target)
                call refine_discrete_root(state, si0, f_land, zu, h_rem, e_land, &
                    scale_target, roots(j), direction, target, si_out, f_out, tau, &
                    done)
            else
                call set_home_volume(state, home)
                call refine_discrete_root(state, si0, f_land, si_land%z, h_rem, &
                    e_land, scale_home, roots(j), -direction, home, si_out, &
                    f_out, tau, done)
            end if
            if (done) exit
        end do
        if (.not. done) then
            call set_home_volume(state, home)
            return
        end if

        info%event_type = kind(j)
        info%vol_to = merge(target, home, kind(j) == CROSS_CROSSING)
        info%lambda = tau
        info%vpar_after = f_out%vpar/sqrt2
        info%bmod_target = f_out%Bmod
        call integ_to_ref(si_out%z(1:3), xinteg)
        info%dtheta = si_out%z(2) - si_land%z(2)
        info%dzeta = unwrap_near(xinteg(3), xref(3)) - xref(3)
        si = si_out
        f = f_out
        call update_landing_stats(resid)
        call count_symplectic_crossing(info%event_type)
    end subroutine symplectic_crossing

    subroutine flow_scale(z, f_land, scale)
        !> Flow length per radian of poloidal motion along X_s, |P_r|, in the
        !> locked chart.
        real(dp), intent(in) :: z(4)
        type(field_can_t), intent(in) :: f_land
        real(dp), intent(out) :: scale

        type(field_can_t) :: f

        f = f_land
        call eval_field(f, z(1), z(2), z(3), 0)
        call get_derivatives(f, z(4))
        scale = abs(f%dpth(1))
        if (.not. ieee_is_finite(scale)) scale = 0.0_dp
    end subroutine flow_scale

    subroutine scan_energy_roots(z_start, f_land, e_land, scale, sigma, roots, &
            n_roots)
        !> Roots of H(z(tau)) - e_land along X_s in the locked chart, tau in
        !> radians of poloidal motion with sign sigma, on a geometric grid out
        !> to 2 pi, in order of increasing |tau|.
        real(dp), intent(in) :: z_start(4), e_land, scale
        type(field_can_t), intent(in) :: f_land
        integer, intent(in) :: sigma
        real(dp), intent(out) :: roots(:)
        integer, intent(out) :: n_roots

        real(dp) :: z(4), tau_prev, tau, g_prev, g
        integer :: k, status
        type(field_can_t) :: f

        n_roots = 0
        if (scale <= 0.0_dp) return
        f = f_land
        z = z_start
        call interface_flow(z, 0.0_dp, f, status)
        tau_prev = 0.0_dp
        g_prev = f%H - e_land
        do k = 1, n_scan
            tau = real(sigma, dp)*scan_tau0*2.0_dp**(k - 1)
            call interface_flow(z, (tau - tau_prev)*scale, f, status)
            if (status /= REFRACTION_OK) exit
            g = f%H - e_land
            if (g_prev*g <= 0.0_dp) then
                n_roots = n_roots + 1
                roots(n_roots) = tau_prev - g_prev*(tau - tau_prev)/(g - g_prev)
            end if
            tau_prev = tau
            g_prev = g
        end do
    end subroutine scan_energy_roots

    subroutine refine_discrete_root(state, si0, f_land, z_start, h_rem, e_land, &
            scale, tau_guess, side, lvol, si_out, f_out, tau, found)
        !> Secant iteration on G(tau) = H_mid(remainder) - e_land.
        type(sympl_spectre_state_t), intent(in) :: state
        type(symplectic_integrator_t), intent(in) :: si0
        type(field_can_t), intent(in) :: f_land
        real(dp), intent(in) :: z_start(4), h_rem, e_land, scale, tau_guess
        integer, intent(in) :: side, lvol
        type(symplectic_integrator_t), intent(out) :: si_out
        type(field_can_t), intent(out) :: f_out
        real(dp), intent(out) :: tau
        logical, intent(out) :: found

        real(dp) :: t0, t1, g0, g1, t2, tol
        integer :: it
        logical :: ok0, ok1

        found = .false.
        tol = refine_tol*max(abs(e_land), 1.0_dp)
        t0 = tau_guess
        t1 = tau_guess + max(abs(tau_guess), 1.0_dp)*1.0d-7
        call remainder_energy(state, si0, f_land, z_start, h_rem, scale, t0, &
            si_out, f_out, g0, ok0)
        g0 = g0 - e_land
        call remainder_energy(state, si0, f_land, z_start, h_rem, scale, t1, &
            si_out, f_out, g1, ok1)
        g1 = g1 - e_land
        if (.not. (ok0 .and. ok1)) return
        do it = 1, max_refine
            if (abs(g1) <= tol) exit
            if (abs(g1 - g0) <= tiny(1.0_dp)) return
            t2 = t1 - g1*(t1 - t0)/(g1 - g0)
            t0 = t1
            g0 = g1
            t1 = t2
            call remainder_energy(state, si0, f_land, z_start, h_rem, scale, t1, &
                si_out, f_out, g1, ok1)
            if (.not. ok1) return
            g1 = g1 - e_land
        end do
        if (abs(g1) > tol) return
        if (real(side, dp)*(si_out%z(1) - z_start(1)) <= 0.0_dp) return
        if (si_out%z(1) < real(lvol - 1, dp) .or. si_out%z(1) > real(lvol, dp)) return
        tau = t1
        found = .true.
    end subroutine refine_discrete_root

    subroutine remainder_energy(state, si0, f_land, z_start, h_rem, scale, tau, &
            si_out, f_out, energy, ok)
        !> Flow tau radians along X_s from z_start, take the remainder midpoint
        !> substep h_rem and return its discrete energy.
        type(sympl_spectre_state_t), intent(in) :: state
        type(symplectic_integrator_t), intent(in) :: si0
        type(field_can_t), intent(in) :: f_land
        real(dp), intent(in) :: z_start(4), h_rem, scale, tau
        type(symplectic_integrator_t), intent(out) :: si_out
        type(field_can_t), intent(out) :: f_out
        real(dp), intent(out) :: energy
        logical, intent(out) :: ok

        type(symplectic_integrator_t) :: si_start
        type(field_can_t) :: fm
        real(dp) :: zu(4), pth_u
        integer :: status, ierr

        ok = .false.
        energy = 0.0_dp
        f_out = f_land
        zu = z_start
        call interface_flow(zu, tau*scale, f_out, status)
        if (status /= REFRACTION_OK) return
        si_out = si0
        call orbit_sympl_init(si_out, f_out, zu, h_rem, 1, si0%rtol, state%mode)
        si_out%warning_factor = symplectic_spectre_newton_warning_factor
        pth_u = f_out%pth
        si_start = si_out
        if (h_rem > 0.0_dp) then
            call orbit_timestep_sympl(si_out, f_out, ierr)
            if (ierr /= 0) return
            if (step_teleported(si_start, si_out, f_out)) return
        end if
        fm = f_out
        call midpoint_energy(zu, pth_u, si_out%z, f_out%pth, fm, energy, status)
        ok = status == REFRACTION_OK .and. ieee_is_finite(energy)
        si_out%dt = state%dt_std
    end subroutine remainder_energy

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
            if (ierr_step == SYMPLECTIC_STEP_OUTSIDE_DOMAIN .and. &
                    nint(state%home_hi) == state%mvol) then
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
            full_committed, h_land, resid, solved, tol)
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
        !> Target |rho_g - rho_k|; default landing_tol.
        real(dp), intent(in), optional :: tol

        type(symplectic_integrator_t) :: si_best
        type(field_can_t) :: f_best
        real(dp) :: dir, t_lo, t_hi, f_lo, f_hi, t_mid, g, best, target_tol
        integer :: it, ierr_step, last_side, evals
        logical :: hi_valid

        target_tol = landing_tol
        if (present(tol)) target_tol = tol

        dir = real(direction, dp)
        t_lo = 0.0_dp
        f_lo = dir*(si0%z(1) - rho_k)
        t_hi = h_full
        f_hi = 0.0_dp
        hi_valid = .false.
        last_side = 0
        evals = 0

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
            if (best < target_tol) exit
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
            evals = evals + 1
            g = dir*(si%z(1) - rho_k)
            if (ierr_step /= 0 .or. g /= g .or. &
                step_teleported(si0, si, f)) then
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
        if (allocated(n_landing_evals)) then
            n_landing_evals(stats_thread()) = &
                n_landing_evals(stats_thread()) + evals
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

        y_out = y
        owner_out = owner
        exited = .false.
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
                if (allocated(n_fo_losses)) n_fo_losses(stats_thread()) = &
                    n_fo_losses(stats_thread()) + 1
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
        si%warning_factor = symplectic_spectre_newton_warning_factor
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

    pure function step_teleported(si0, si, f) result(bad)
        type(symplectic_integrator_t), intent(in) :: si0, si
        type(field_can_t), intent(in) :: f
        logical :: bad

        bad = .not. all(ieee_is_finite(si%z)) .or. &
            .not. ieee_is_finite(si%pabs) .or. &
            .not. ieee_is_finite(f%Bmod) .or. &
            .not. ieee_is_finite(f%vpar) .or. &
            .not. ieee_is_finite(f%mu) .or. &
            .not. ieee_is_finite(f%H)
        if (bad) return
        bad = abs(si%z(1) - si0%z(1)) > max_step_dr
        if (bad) return
        bad = abs(f%H - si%pabs**2) > &
            max_relative_shell_defect*max(si%pabs**2, tiny(1.0_dp))
        if (bad) return
        ! The steppers commit a negative-radius Newton root as an axis chart
        ! flip (r, theta) -> (|r|, theta + pi). Away from the axis a negative
        ! root is a Newton divergence, not an axis crossing: a committed theta
        ! jump of order pi at large radius is a teleport.
        if (abs(si%z(2) - si0%z(2)) > 1.0_dp) then
            bad = si%z(1) > 0.5_dp
            if (bad) return
        end if
    end function step_teleported

    subroutine recover_local_error(state, si, f, ipart, t_sec, reason, budget, &
            used, ierr, t_frac, resume_gc)
        !> A bulk or axis-local solver failure has no interface metadata. Do not
        !> invent an interior interface from nint(rho): interface FO recovery is
        !> valid only after an actual landed interface event.
        type(sympl_spectre_state_t), intent(inout) :: state
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(in) :: ipart, reason
        real(dp), intent(in) :: t_sec
        real(dp), intent(inout) :: budget, used
        integer, intent(out) :: ierr
        real(dp), intent(out) :: t_frac
        logical, intent(out) :: resume_gc

        real(dp) :: y(5), y_out(5)
        integer :: ierr_fo, owner, owner_out
        logical :: fo_exited

        resume_gc = .false.
        call integrator_state(si, f, y)
        owner = nint(state%home_hi)
        call record_bulk_recovery(si, f, ipart, t_sec)
        call start_fo(state, si, f, y, owner, 0, budget, used, y_out, &
            owner_out, fo_exited, ierr_fo)
        if (ierr_fo == SPECTRE_FO_OK) then
            resume_gc = fo_exited .and. budget > budget_eps*state%dt_std
            return
        end if
        if (ierr_fo == SPECTRE_FO_LOSS) then
            call finish_fo_error(ierr_fo, si, f, ipart, t_sec, 0, 0, reason, &
                used, state%dt_std, ierr, t_frac)
            return
        end if
        ! Bulk recovery has no interface or sheet. Preserve the last valid FO
        ! state if one exists, but retain the original local failure reason and
        ! iface=0 in the marker-local stop record.
        if (state%fo%has_y) then
            call set_home_volume(state, owner_out)
            call recanonicalize(state, si, f, state%fo%last_y)
        end if

        if (recover_warning_skip(state, si, budget, used, ierr, t_frac, reason, &
            0, 0)) return

        call record_stop(si, f, ipart, t_sec, 0, 0, reason)
        ierr = SYMPL_SPECTRE_STOP
        t_frac = used/state%dt_std
    end subroutine recover_local_error

    subroutine record_bulk_recovery(si, f, ipart, t_sec)
        type(symplectic_integrator_t), intent(in) :: si
        type(field_can_t), intent(in) :: f
        integer, intent(in) :: ipart
        real(dp), intent(in) :: t_sec

        type(crossing_info_t) :: info
        real(dp) :: vpar

        vpar = f%vpar/sqrt2
        info%event_type = CROSS_RECOVERY
        info%iface = 0
        info%vol_from = 0
        info%vol_to = 0
        info%theta = si%z(2)
        info%zeta = si%z(3)
        info%vpar_before = vpar
        info%vpar_after = vpar
        info%mu = 0.5_dp*f%mu
        info%bmod_home = f%Bmod
        info%bmod_target = f%Bmod
        call crossing_log_record(ipart, t_sec, info)
    end subroutine record_bulk_recovery

    pure function unwrap_near(a, near) result(b)
        real(dp), intent(in) :: a, near
        real(dp) :: b

        b = a + twopi*nint((near - a)/twopi)
    end function unwrap_near

    logical function recover_fo_failure(state, si, f, fo_status, owner, budget, &
            used, ierr, t_frac)
        !> Full orbit is the final active recovery. If it also fails numerically,
        !> return to its last valid guiding-centre state and use the warning-mode
        !> hold fallback. A physical full-orbit loss is never intercepted.
        type(sympl_spectre_state_t), intent(inout) :: state
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(in) :: fo_status, owner
        real(dp), intent(inout) :: budget, used
        integer, intent(inout) :: ierr
        real(dp), intent(inout) :: t_frac

        recover_fo_failure = .false.
        if (fo_status == SPECTRE_FO_LOSS) return
        if (.not. symplectic_newton_warning_mode) return
        if (state%fo%has_y) then
            call set_home_volume(state, owner)
            call recanonicalize(state, si, f, state%fo%last_y)
        end if
        recover_fo_failure = recover_warning_skip(state, si, budget, used, &
            ierr, t_frac, STOP_SHEET, state%fo%iface, &
            merge(1, -1, owner == state%fo%iface))
    end function recover_fo_failure

    logical function recover_warning_skip(state, si, budget, used, ierr, &
            t_frac, reason, iface, direction)
        !> Last-resort default: after all active recovery has failed, retain the
        !> last valid state for this unresolved interval and resume the marker on
        !> the next microstep. Strict mode returns .false. to expose the stop.
        type(sympl_spectre_state_t), intent(inout) :: state
        type(symplectic_integrator_t), intent(in) :: si
        real(dp), intent(inout) :: budget, used
        integer, intent(inout) :: ierr
        real(dp), intent(inout) :: t_frac
        integer, intent(in) :: reason, iface, direction

        recover_warning_skip = symplectic_newton_warning_mode
        if (.not. recover_warning_skip) return

        used = used + budget
        budget = 0.0_dp
        ierr = SYMPL_SPECTRE_OK
        t_frac = 1.0_dp
        call count_event(EVT_WARNING_STEP_SKIP)

        state%warning_hold_latched = .true.
        state%warning_hold_count = state%warning_hold_count + 1
        state%warning_hold_reason = reason
        state%warning_hold_iface = iface
        state%warning_hold_direction = direction
        state%warning_hold_z = si%z

        state%fo%active = .false.
        state%sheet%active = .false.
    end function recover_warning_skip

    pure logical function same_warning_hold(state, si, reason, iface, direction)
        type(sympl_spectre_state_t), intent(in) :: state
        type(symplectic_integrator_t), intent(in) :: si
        integer, intent(in) :: reason, iface, direction

        same_warning_hold = state%warning_hold_latched .and. &
            reason == state%warning_hold_reason .and. &
            iface == state%warning_hold_iface .and. &
            direction == state%warning_hold_direction .and. &
            maxval(abs(si%z - state%warning_hold_z)) <= &
            100.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
            maxval(abs(state%warning_hold_z)))
    end function same_warning_hold

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
        if (iface == 0) then
            info%vol_from = 0
        else
            info%vol_from = merge(iface, iface + 1, direction == 1)
        end if
        info%vol_to = info%vol_from
        info%theta = si%z(2)
        info%zeta = si%z(3)
        info%vpar_before = vpar
        info%vpar_after = vpar
        info%mu = 0.5_dp*f%mu
        info%bmod_home = f%Bmod
        info%bmod_target = f%Bmod
        call crossing_log_record(ipart, t_sec, info)

        if (allocated(n_stops)) then
            n_stops(stats_thread()) = n_stops(stats_thread()) + 1
            if (reason >= 1 .and. reason <= size(n_stop_reason, 1)) &
                n_stop_reason(reason, stats_thread()) = &
                    n_stop_reason(reason, stats_thread()) + 1
        end if
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

        if (allocated(n_landings)) then
            n_landings(stats_thread()) = n_landings(stats_thread()) + 1
            max_landing_resid(stats_thread()) = &
                max(max_landing_resid(stats_thread()), resid)
        end if
    end subroutine update_landing_stats

    subroutine sympl_landing_stats_reset
        integer :: nthreads

        nthreads = 1
        !$ nthreads = omp_get_max_threads()
        if (allocated(n_landings)) then
            deallocate(n_landings, n_landing_evals, n_stops, n_sheet_entries, &
                n_sheet_exits, n_sheet_init_failures, n_sheet_advance_failures, &
                n_sheet_failure_status, n_stop_reason, n_fo_entries, n_fo_exits, &
                n_fo_losses, n_fo_failures, n_fo_failure_status, &
                max_landing_resid, n_sympl_transmissions, n_sympl_reflections, &
                n_sympl_fallbacks)
        end if
        allocate(n_landings(0:nthreads - 1), n_landing_evals(0:nthreads - 1), &
            n_stops(0:nthreads - 1), n_sheet_entries(0:nthreads - 1), &
            n_sheet_exits(0:nthreads - 1), &
            n_sheet_init_failures(0:nthreads - 1), &
            n_sheet_advance_failures(0:nthreads - 1), &
            n_sheet_failure_status(5, 0:nthreads - 1), &
            n_stop_reason(5, 0:nthreads - 1), n_fo_entries(0:nthreads - 1), &
            n_fo_exits(0:nthreads - 1), n_fo_losses(0:nthreads - 1), &
            n_fo_failures(0:nthreads - 1), &
            n_fo_failure_status(5, 0:nthreads - 1), &
            max_landing_resid(0:nthreads - 1), &
            n_sympl_transmissions(0:nthreads - 1), &
            n_sympl_reflections(0:nthreads - 1), n_sympl_fallbacks(0:nthreads - 1))
        n_landings = 0
        n_landing_evals = 0
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
        n_sympl_transmissions = 0
        n_sympl_reflections = 0
        n_sympl_fallbacks = 0
    end subroutine sympl_landing_stats_reset

    subroutine sympl_landing_stats(landings, max_resid, stops)
        integer, intent(out) :: landings, stops
        real(dp), intent(out) :: max_resid

        landings = 0
        max_resid = 0.0_dp
        stops = 0
        if (allocated(n_landings)) then
            landings = sum(n_landings)
            max_resid = maxval(max_landing_resid)
            stops = sum(n_stops)
        end if
    end subroutine sympl_landing_stats

    subroutine sympl_landing_eval_stats(evals)
        !> Total implicit-step evaluations spent inside locate_landing since
        !> the last reset; evaluations/landing is the root-iteration cost the
        !> seeded bracket shrink is measured by.
        integer, intent(out) :: evals

        evals = 0
        if (allocated(n_landing_evals)) evals = sum(n_landing_evals)
    end subroutine sympl_landing_eval_stats

    subroutine count_sheet_stat(kind, status)
        integer, intent(in) :: kind
        integer, intent(in), optional :: status

        integer :: tid

        if (.not. allocated(n_sheet_entries)) return
        tid = stats_thread()
        select case (kind)
        case (SHEET_STAT_ENTRY)
            n_sheet_entries(tid) = n_sheet_entries(tid) + 1
        case (SHEET_STAT_EXIT)
            n_sheet_exits(tid) = n_sheet_exits(tid) + 1
        case (SHEET_STAT_INIT_FAIL)
            n_sheet_init_failures(tid) = n_sheet_init_failures(tid) + 1
        case (SHEET_STAT_ADVANCE_FAIL)
            n_sheet_advance_failures(tid) = n_sheet_advance_failures(tid) + 1
        end select
        if (present(status)) then
            if (status >= 1 .and. status <= size(n_sheet_failure_status, 1)) &
                n_sheet_failure_status(status, tid) = &
                    n_sheet_failure_status(status, tid) + 1
        end if
    end subroutine count_sheet_stat

    subroutine sympl_sheet_stats(entries, exits, init_failures, advance_failures, &
            failure_status, stop_reason)
        integer, intent(out) :: entries, exits, init_failures, advance_failures
        integer, intent(out) :: failure_status(5)
        integer, intent(out) :: stop_reason(5)

        entries = 0
        exits = 0
        init_failures = 0
        advance_failures = 0
        failure_status = 0
        stop_reason = 0
        if (allocated(n_sheet_entries)) then
            entries = sum(n_sheet_entries)
            exits = sum(n_sheet_exits)
            init_failures = sum(n_sheet_init_failures)
            advance_failures = sum(n_sheet_advance_failures)
            failure_status = sum(n_sheet_failure_status, dim=2)
            stop_reason = sum(n_stop_reason, dim=2)
        end if
    end subroutine sympl_sheet_stats

    subroutine count_fo_entry
        if (allocated(n_fo_entries)) n_fo_entries(stats_thread()) = &
            n_fo_entries(stats_thread()) + 1
    end subroutine count_fo_entry

    subroutine count_fo_exit
        if (allocated(n_fo_exits)) n_fo_exits(stats_thread()) = &
            n_fo_exits(stats_thread()) + 1
    end subroutine count_fo_exit

    subroutine count_fo_failure(status)
        integer, intent(in) :: status
        integer :: tid

        if (.not. allocated(n_fo_failures)) return
        tid = stats_thread()
        n_fo_failures(tid) = n_fo_failures(tid) + 1
        if (status >= 1 .and. status <= size(n_fo_failure_status, 1)) &
            n_fo_failure_status(status, tid) = &
                n_fo_failure_status(status, tid) + 1
    end subroutine count_fo_failure

    subroutine sympl_fo_stats(entries, exits, losses, failures, failure_status)
        integer, intent(out) :: entries, exits, losses, failures
        integer, intent(out) :: failure_status(5)

        entries = 0
        exits = 0
        losses = 0
        failures = 0
        failure_status = 0
        if (allocated(n_fo_entries)) then
            entries = sum(n_fo_entries)
            exits = sum(n_fo_exits)
            losses = sum(n_fo_losses)
            failures = sum(n_fo_failures)
            failure_status = sum(n_fo_failure_status, dim=2)
        end if
    end subroutine sympl_fo_stats

    subroutine sympl_symplectic_stats(transmissions, reflections, fallbacks)
        integer, intent(out) :: transmissions, reflections, fallbacks

        transmissions = 0
        reflections = 0
        fallbacks = 0
        if (allocated(n_sympl_transmissions)) then
            transmissions = sum(n_sympl_transmissions)
            reflections = sum(n_sympl_reflections)
            fallbacks = sum(n_sympl_fallbacks)
        end if
    end subroutine sympl_symplectic_stats

    subroutine count_symplectic_crossing(event_type)
        !> event_type CROSS_CROSSING or CROSS_REFLECTION; any other value counts a
        !> fallback to the Level-1 map.
        integer, intent(in) :: event_type

        integer :: tid

        if (.not. allocated(n_sympl_transmissions)) return
        tid = stats_thread()
        select case (event_type)
        case (CROSS_CROSSING)
            n_sympl_transmissions(tid) = n_sympl_transmissions(tid) + 1
        case (CROSS_REFLECTION)
            n_sympl_reflections(tid) = n_sympl_reflections(tid) + 1
        case default
            n_sympl_fallbacks(tid) = n_sympl_fallbacks(tid) + 1
        end select
    end subroutine count_symplectic_crossing

    integer function stats_thread()
        stats_thread = 0
        !$ stats_thread = omp_get_thread_num()
    end function stats_thread

end module spectre_sympl_orbit
