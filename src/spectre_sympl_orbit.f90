module spectre_sympl_orbit
    !> Per-volume guiding-centre stepping with full-orbit interface excursions.
    !> Interior sheets are crossed only by the full orbit; no fixed-mu map.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    use util, only: twopi, sqrt2
    use parmot_mod, only: ro0
    use field_can_mod, only: field_can_t, eval_field => evaluate, integ_to_ref, &
        ref_to_integ
    use field_can_spectre, only: set_spectre_volume_lock
    use orbit_symplectic_base, only: symplectic_integrator_t, &
        symplectic_newton_warning_mode, symplectic_spectre_newton_warning_factor, &
        SYMPLECTIC_STEP_OUTSIDE_DOMAIN
    use orbit_symplectic, only: orbit_timestep_sympl, orbit_sympl_init
    use spectre_fo_hybrid, only: spectre_fo_state_t, spectre_fo_enter, &
        spectre_fo_advance_until_exit, spectre_fo_radial_larmor, SPECTRE_FO_OK, &
        SPECTRE_FO_LOSS
    use params, only: spectre_fo_entry_larmor
    use diag_counters, only: count_event, EVT_WARNING_STEP_SKIP

    implicit none
    private

    public :: sympl_spectre_state_t, sympl_spectre_reset, recanon_pphi, &
        orbit_microstep_sympl_spectre, sympl_landing_stats_reset, &
        sympl_landing_stats, sympl_landing_eval_stats
    public :: sympl_fo_stats
    public :: SYMPL_SPECTRE_OK, SYMPL_SPECTRE_LOSS, SYMPL_SPECTRE_STOP

    integer, parameter :: SYMPL_SPECTRE_OK = 0
    integer, parameter :: SYMPL_SPECTRE_LOSS = 88
    integer, parameter :: SYMPL_SPECTRE_STOP = 89
    integer, parameter :: STOP_SHEET = 1, STOP_TELEPORT = 2, STOP_STEP = 3
    !> Landing target and acceptance for |rho_g - k| at the interface. The
    !> h*-solve aims at landing_tol; the per-volume splines are clamped at the
    !> volume edge, which kinks the landed rho_g just outside the interface, so a
    !> best-effort landing below landing_accept is still committed.
    real(dp), parameter :: landing_tol = 1.0d-10
    real(dp), parameter :: landing_accept = 1.0d-8
    integer, parameter :: max_locate = 64

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
    !> the fixed physical speed pabs. GC/FO conversion preserves the momentum shell.
    real(dp), parameter :: max_step_dr = 0.5d0
    real(dp), parameter :: max_relative_shell_defect = 5.0d-4

    !> Home-side offset for the per-volume angle-transform evaluation at a landed
    !> point: the volume dispatch keys on int(rho_g), so exactly at rho_g = k it
    !> would pick the outer volume, whose zeta gauge differs by O(1).
    real(dp), parameter :: transform_bias = 1.0d-8

    type :: sympl_spectre_state_t
        real(dp) :: home_lo = 0.0_dp
        real(dp) :: home_hi = 1.0_dp
        integer :: mvol = 1
        integer :: mode = 0
        integer :: level = 3
        real(dp) :: dt_std = 0.0_dp
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
        type(spectre_fo_state_t) :: fo
    end type sympl_spectre_state_t

    !> Hot-path statistics are per-thread. Each OpenMP worker writes only its
    !> own slot; totals are reduced after the parallel trace has joined.
    integer, allocatable :: n_landings(:), n_landing_evals(:), n_stops(:)
    integer, allocatable :: n_fo_entries(:), n_fo_exits(:), n_fo_losses(:)
    integer, allocatable :: n_fo_failures(:), n_fo_failure_status(:, :)
    real(dp), allocatable :: max_landing_resid(:)

contains

    subroutine sympl_spectre_reset(state, si, mvol, mode, level)
        type(sympl_spectre_state_t), intent(out) :: state
        type(symplectic_integrator_t), intent(in) :: si
        integer, intent(in) :: mvol, mode, level

        if (level /= 3) error stop 'crossing_level must be 3; GC maps were removed'
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

    subroutine orbit_microstep_sympl_spectre(state, si, f, ierr, t_frac)
        !> Advance one microstep of length state%dt_std, resolving interface
        !> events by exact-landing substeps. On SYMPL_SPECTRE_LOSS or
        !> SYMPL_SPECTRE_STOP the orbit terminates at fraction t_frac of the
        !> microstep, with (si, f) holding the terminal (landed) state.
        type(sympl_spectre_state_t), intent(inout) :: state
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(out) :: ierr
        real(dp), intent(out) :: t_frac

        type(symplectic_integrator_t) :: si0
        type(field_can_t) :: f0
        real(dp) :: budget, used, h_try, h_land, resid, y_iface(5), y_out(5)
        real(dp) :: exit_rho
        integer :: iters, iface, direction, ierr_step, ierr_fo
        integer :: fo_owner
        logical :: boundary, solved, fo_exited, resume_gc

        ierr = SYMPL_SPECTRE_OK
        t_frac = 1.0_dp
        si%warning_factor = symplectic_spectre_newton_warning_factor
        budget = state%dt_std
        used = 0.0_dp
        h_try = budget

        if (state%fo%active) then
            call continue_fo(state, si, f, budget, used, y_out, fo_owner, &
                fo_exited, ierr_fo)
            if (ierr_fo /= SPECTRE_FO_OK) then
                if (recover_fo_failure(state, si, f, ierr_fo, fo_owner, &
                    budget, used, ierr, t_frac)) then
                    si%dt = state%dt_std
                    return
                end if
                call finish_fo_error(ierr_fo, used, state%dt_std, ierr, t_frac)
                si%dt = state%dt_std
                return
            end if
            if (.not. fo_exited .or. budget <= budget_eps*state%dt_std) then
                si%dt = state%dt_std
                return
            end if
            h_try = budget
        end if

        iters = 0
        do
            iters = iters + 1
            if (iters > max_iters) then
                call recover_local_error(state, si, f, STOP_STEP, budget, used, &
                    ierr, t_frac, resume_gc)
                if (resume_gc) then
                    iters = 0
                    h_try = budget
                    cycle
                end if
                exit
            end if
            if (state%level == 3) then
                call interface_proximity(state, si, f, y_iface, iface, exit_rho)
                if (iface > 0) then
                    call start_fo(state, si, f, y_iface, nint(state%home_hi), &
                        iface, budget, used, y_out, fo_owner, fo_exited, ierr_fo, &
                        exit_rho)
                    if (ierr_fo /= SPECTRE_FO_OK) then
                        if (recover_fo_failure(state, si, f, ierr_fo, fo_owner, &
                            budget, used, ierr, t_frac)) exit
                        call finish_fo_error(ierr_fo, used, state%dt_std, ierr, &
                            t_frac)
                        exit
                    end if
                    if (.not. fo_exited .or. &
                        budget <= budget_eps*state%dt_std) exit
                    h_try = budget
                    cycle
                end if
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
                    call recover_local_error(state, si, f, STOP_TELEPORT, budget, &
                        used, ierr, t_frac, resume_gc)
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
                    call recover_local_error(state, si, f, STOP_STEP, budget, &
                        used, ierr, t_frac, resume_gc)
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
                    call finish_fo_error(ierr_fo, used, state%dt_std, ierr, t_frac)
                    exit
                end if
                if (.not. fo_exited .or. &
                    budget <= budget_eps*state%dt_std) exit
                h_try = budget
                cycle
            end if
            call landed_state(si, f, real(iface, dp), direction, y_iface)
            used = used + h_land
            call update_landing_stats(resid)
            if (state%mvol == 1) then
                ! The separate single-volume GC model retains its GC edge.
                ierr = SYMPL_SPECTRE_LOSS
                t_frac = used/state%dt_std
                exit
            end if
            ! A multi-volume hybrid resolves both sheets and its absorbing
            ! particle boundary in FO if a GC step missed the proximity band.
            budget = state%dt_std - used
            call fo_exit_distance(state, y_iface, exit_rho)
            call start_fo(state, si, f, y_iface, nint(state%home_hi), iface, &
                budget, used, y_out, fo_owner, fo_exited, ierr_fo, exit_rho)
            if (ierr_fo /= SPECTRE_FO_OK) then
                if (recover_fo_failure(state, si, f, ierr_fo, fo_owner, &
                    budget, used, ierr, t_frac)) exit
                call finish_fo_error(ierr_fo, used, state%dt_std, ierr, t_frac)
                exit
            end if
            if (.not. fo_exited .or. budget <= budget_eps*state%dt_std) exit
            h_try = budget
            cycle

        end do

        si%dt = state%dt_std
    end subroutine orbit_microstep_sympl_spectre

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
            ! A trial exceeding the Newton margin (Mvol + 1) from the outer
            ! volume still brackets its physical edge at Mvol, with si%z
            ! holding the pre-step state. Successful trials are checked below.
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
        integer :: it, ierr_step, last_side, evals
        logical :: hi_valid

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

    subroutine interface_proximity(state, si, f, y, iface, exit_rho)
        !> Nearest interface or absorbing outer edge if the GC is closer to
        !> it than spectre_fo_entry_larmor radial gyroradii, else iface = 0. y
        !> is the guiding-centre state and exit_rho the matching exit distance.
        type(sympl_spectre_state_t), intent(in) :: state
        type(symplectic_integrator_t), intent(in) :: si
        type(field_can_t), intent(in) :: f
        real(dp), intent(out) :: y(5), exit_rho
        integer, intent(out) :: iface

        real(dp) :: rho_r
        integer :: k, status

        iface = 0
        exit_rho = 0.0_dp
        call integrator_state(si, f, y)
        if (state%mvol == 1) return
        k = nint(y(1))
        if (k < 1 .or. k > state%mvol) return
        call spectre_fo_radial_larmor(y, nint(state%home_hi), ro0/sqrt2, rho_r, &
            status)
        if (status /= SPECTRE_FO_OK) return
        if (abs(y(1) - real(k, dp)) >= spectre_fo_entry_larmor*rho_r) return
        iface = k
        exit_rho = (spectre_fo_entry_larmor + 1.0_dp)*rho_r
    end subroutine interface_proximity

    subroutine fo_exit_distance(state, y, exit_rho)
        !> Exit distance (rho_g) for a full orbit started on an interface.
        type(sympl_spectre_state_t), intent(in) :: state
        real(dp), intent(in) :: y(5)
        real(dp), intent(out) :: exit_rho

        real(dp) :: rho_r
        integer :: status

        call spectre_fo_radial_larmor(y, nint(state%home_hi), ro0/sqrt2, rho_r, &
            status)
        exit_rho = 0.0_dp
        if (status == SPECTRE_FO_OK) &
            exit_rho = (spectre_fo_entry_larmor + 1.0_dp)*rho_r
    end subroutine fo_exit_distance

    subroutine start_fo(state, si, f, y, owner, iface, budget, used, y_out, &
            owner_out, exited, ierr, exit_rho)
        type(sympl_spectre_state_t), intent(inout) :: state
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        real(dp), intent(in) :: y(5)
        integer, intent(in) :: owner, iface
        real(dp), intent(inout) :: budget, used
        real(dp), intent(out) :: y_out(5)
        integer, intent(out) :: owner_out, ierr
        logical, intent(out) :: exited
        !> Leave the full orbit once the guiding centre is this far (rho_g)
        !> from the interface; default: the gyrophase-based legacy rule.
        real(dp), intent(in), optional :: exit_rho

        y_out = y
        owner_out = owner
        exited = .false.
        call spectre_fo_enter(state%fo, y, owner, iface, ro0/sqrt2, ierr)
        if (ierr == SPECTRE_FO_LOSS) then
            ! A GC initially inside the wall can map to a particle outside it.
            ! This is an immediate physical loss, with no accepted FO time.
            call count_fo_entry
            if (allocated(n_fo_losses)) n_fo_losses(stats_thread()) = &
                n_fo_losses(stats_thread()) + 1
            return
        end if
        if (ierr /= SPECTRE_FO_OK) then
            call count_fo_failure(ierr)
            return
        end if
        if (present(exit_rho)) state%fo%exit_rho = exit_rho
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

    subroutine recover_local_error(state, si, f, reason, budget, used, &
            ierr, t_frac, resume_gc)
        !> A bulk or axis-local solver failure has no interface metadata. Do not
        !> invent an interior interface from nint(rho): interface FO recovery is
        !> valid only after an actual landed interface event.
        type(sympl_spectre_state_t), intent(inout) :: state
        type(symplectic_integrator_t), intent(inout) :: si
        type(field_can_t), intent(inout) :: f
        integer, intent(in) :: reason
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
        call start_fo(state, si, f, y, owner, 0, budget, used, y_out, &
            owner_out, fo_exited, ierr_fo)
        if (ierr_fo == SPECTRE_FO_OK) then
            resume_gc = fo_exited .and. budget > budget_eps*state%dt_std
            return
        end if
        if (ierr_fo == SPECTRE_FO_LOSS) then
            call finish_fo_error(ierr_fo, used, state%dt_std, ierr, t_frac)
            return
        end if
        ! Bulk recovery has no interface or sheet. Preserve the last valid FO
        ! state if one exists, but retain the original local failure reason and
        ! the local stop statistics.
        if (state%fo%has_y) then
            call set_home_volume(state, owner_out)
            call recanonicalize(state, si, f, state%fo%last_y)
        end if

        if (recover_warning_skip(state, si, budget, used, ierr, t_frac, reason, &
            0, 0)) return

        call record_stop
        ierr = SYMPL_SPECTRE_STOP
        t_frac = used/state%dt_std
    end subroutine recover_local_error

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

    subroutine record_stop
        if (allocated(n_stops)) then
            n_stops(stats_thread()) = n_stops(stats_thread()) + 1
        end if
    end subroutine record_stop

    subroutine finish_fo_error(fo_status, used, dt, ierr, t_frac)
        integer, intent(in) :: fo_status
        real(dp), intent(in) :: used, dt
        integer, intent(out) :: ierr
        real(dp), intent(out) :: t_frac

        if (fo_status == SPECTRE_FO_LOSS) then
            ierr = SYMPL_SPECTRE_LOSS
        else
            call record_stop
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
            deallocate(n_landings, n_landing_evals, n_stops, &
                n_fo_entries, n_fo_exits, &
                n_fo_losses, n_fo_failures, n_fo_failure_status, &
                max_landing_resid)
        end if
        allocate(n_landings(0:nthreads - 1), n_landing_evals(0:nthreads - 1), &
            n_stops(0:nthreads - 1), n_fo_entries(0:nthreads - 1), &
            n_fo_exits(0:nthreads - 1), n_fo_losses(0:nthreads - 1), &
            n_fo_failures(0:nthreads - 1), &
            n_fo_failure_status(5, 0:nthreads - 1), &
            max_landing_resid(0:nthreads - 1))
        n_landings = 0
        n_landing_evals = 0
        n_stops = 0
        n_fo_entries = 0
        n_fo_exits = 0
        n_fo_losses = 0
        n_fo_failures = 0
        n_fo_failure_status = 0
        max_landing_resid = 0.0_dp
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

    integer function stats_thread()
        stats_thread = 0
        !$ stats_thread = omp_get_thread_num()
    end function stats_thread

end module spectre_sympl_orbit
