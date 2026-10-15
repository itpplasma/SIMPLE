module spectre_fo_hybrid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie_sub, only: spectre_field
    use spectre_fo_field, only: spectre_fo_field_build, spectre_fo_field_local, &
        spectre_fo_field_locate, spectre_fo_local_t, fo_field_axisymmetric, &
        FO_AXIS_RADIUS
    implicit none
    private
    integer, parameter, public :: SPECTRE_FO_OK = 0
    integer, parameter, public :: SPECTRE_FO_UNINITIALIZED = 1
    integer, parameter, public :: SPECTRE_FO_LOCATE_FAIL = 2
    integer, parameter, public :: SPECTRE_FO_GYROPHASE_FAIL = 3
    integer, parameter, public :: SPECTRE_FO_INVALID_STATE = 4
    integer, parameter, public :: SPECTRE_FO_STEP_LIMIT = 5
    integer, parameter, public :: SPECTRE_FO_LOSS = 6
    real(dp), parameter :: OWNER_EPS = 1.0e-12_dp
    integer, parameter :: STAGE_IN = 0, STAGE_OUT = 1
    !> RK45 relative tolerance, configured before parallel tracing.
    real(dp) :: fo_tolerance = 1.0e-9_dp
    ! At most 0.1 rad (1/63 of a gyro-period). A cross-and-return hidden
    ! between stage samples has circular-orbit depth <= rho_L*(1-cos(0.05)),
    ! about 1.25e-3 rho_L; tolerance alone cannot bound a missed sheet force.
    real(dp), parameter :: RK_GYRO_ANGLE = 0.1_dp
    !> Gyroangle between exit checks. The guiding-centre map costs about as much
    !> as an RK45 substep, so it is not evaluated after every substep; an exit
    !> is detected at most this angle late, and always at the end of a call.
    real(dp), parameter :: EXIT_CHECK_ANGLE = 1.0_dp
    real(dp), parameter :: ENTRY_MOMENTUM_TOL = 1.0e-3_dp
    integer, parameter :: MAX_SUBSTEPS = 100000
    integer, parameter :: NPHASE = 64
    type, public :: spectre_fo_state_t
        real(dp) :: x(3) = 0.0_dp
        real(dp) :: v(3) = 0.0_dp
        real(dp) :: u(3) = 0.0_dp
        real(dp) :: last_y(5) = 0.0_dp
        real(dp) :: entry_rho = 0.0_dp
        real(dp) :: entry_pzeta = 0.0_dp
        ! Actual FO minus target GC momentum, expressed in the entry gauge.
        real(dp) :: entry_pzeta_defect = 0.0_dp
        !> Next RK45 step (adapted between substeps).
        real(dp) :: h_rk = 0.0_dp
        !> Volume whose gauge entry_pzeta is expressed in.
        integer :: entry_owner = 0
        !> Exit when the guiding centre is farther than this from the interface
        !> (rho_g); 0 keeps the gyrophase-based legacy rule.
        real(dp) :: exit_rho = 0.0_dp
        integer :: owner = 0
        integer :: iface = 0
        logical :: has_y = .false.
        logical :: active = .false.
    end type spectre_fo_state_t
    public :: spectre_fo_enter, spectre_fo_to_gc, spectre_fo_configure
    public :: spectre_fo_radial_larmor
    public :: spectre_fo_advance_until_exit
    public :: spectre_fo_canonical_pzeta
contains
    subroutine spectre_fo_configure(tolerance, grid)
        real(dp), intent(in) :: tolerance
        integer, intent(in), optional :: grid(3)
        integer :: n(3)

        if (tolerance <= 0.0_dp) &
            error stop 'spectre_fo_configure: tolerance must be positive'
        fo_tolerance = tolerance
        n = [32, 64, 32]
        if (present(grid)) n = grid
        if (all(spectre_field%data%in == 0)) n(3) = 7
        call spectre_fo_field_build(n(1), n(2), n(3))
    end subroutine spectre_fo_configure
    subroutine spectre_fo_enter(state, y, owner, iface, ro0_bar, status)
        type(spectre_fo_state_t), intent(out) :: state
        real(dp), intent(in) :: y(5), ro0_bar
        integer, intent(in) :: owner, iface
        integer, intent(out) :: status
        real(dp) :: x_gc(3), b(3), e1(3), e2(3), eperp(3), e_zeta(3)
        real(dp) :: acov(3), bmod, vpar, vperp, target, alpha
        integer :: locate_status
        state = spectre_fo_state_t()
        if (ro0_bar <= 0.0_dp .or. y(4) <= 0.0_dp) then
            status = SPECTRE_FO_INVALID_STATE
            return
        end if
        state%u = y(1:3)
        state%last_y = y
        state%entry_rho = y(1)
        state%has_y = .true.
        state%owner = owner
        state%iface = iface
        call logical_sample(state%u, owner, x_gc, b, bmod, acov, e_zeta, status)
        if (status /= SPECTRE_FO_OK) return
        vpar = sqrt(2.0_dp)*y(4)*y(5)
        vperp = sqrt(max(2.0_dp*y(4)**2 - vpar**2, 0.0_dp))
        call perpendicular_frame(b, e1, e2)
        target = vpar*dot_product(b, e_zeta) + acov(3)/ro0_bar
        state%entry_pzeta = target
        state%entry_owner = owner
        alpha = 0.0_dp
        if (axisymmetric()) then
            call solve_gyrophase(x_gc, state%u, owner, b, e1, e2, vpar, &
                vperp, target, ro0_bar, alpha, status)
            if (status /= SPECTRE_FO_OK) return
        end if
        eperp = cos(alpha)*e1 + sin(alpha)*e2
        state%v = vpar*b + vperp*eperp
        state%x = x_gc + (ro0_bar/bmod)*cross(b, vperp*eperp)
        call spectre_fo_field_locate(state%x, state%u, state%owner, locate_status)
        if (locate_status /= 0) then
            status = SPECTRE_FO_LOCATE_FAIL
            if (locate_status == 1 .and. &
                state%u(1) > real(spectre_field%data%Mvol, dp)) then
                state%active = .true.
                status = SPECTRE_FO_LOSS
            end if
            return
        end if
        if (axisymmetric()) then
            call phase_residual(alpha, x_gc, y(1:3), owner, b, e1, e2, vpar, &
                vperp, target, ro0_bar, state%entry_pzeta_defect, status)
            if (status /= SPECTRE_FO_OK) return
        end if
        state%active = .true.
    end subroutine spectre_fo_enter
    subroutine spectre_fo_radial_larmor(y, owner, ro0_bar, rho_r, status)
        !> Gyroradius of the guiding-centre state y in rho_g units.
        real(dp), intent(in) :: y(5), ro0_bar
        integer, intent(in) :: owner
        real(dp), intent(out) :: rho_r
        integer, intent(out) :: status
        real(dp) :: x(3), b(3), bmod, acov(3), e_zeta(3), vpar, vperp

        rho_r = 0.0_dp
        call logical_sample(y(1:3), owner, x, b, bmod, acov, e_zeta, status)
        if (status /= SPECTRE_FO_OK) return
        vpar = sqrt(2.0_dp)*y(4)*y(5)
        vperp = sqrt(max(2.0_dp*y(4)**2 - vpar**2, 0.0_dp))
        rho_r = ro0_bar*vperp/bmod*radial_gradient(y(1:3))
    end subroutine spectre_fo_radial_larmor
    real(dp) function radial_gradient(u) result(grad)
        real(dp), intent(in) :: u(3)
        type(spectre_fo_local_t) :: loc

        call spectre_fo_field_local(u, min(int(u(1)) + 1, &
            spectre_field%data%Mvol), loc)
        grad = norm2(loc%dual(:, 1))
    end function radial_gradient
    subroutine spectre_fo_advance_until_exit(state, dt, ro0_bar, dt_used, y, &
            owner, exited, status, h_scale)
        type(spectre_fo_state_t), intent(inout) :: state
        real(dp), intent(in) :: dt, ro0_bar
        real(dp), intent(out) :: dt_used, y(5)
        integer, intent(out) :: owner, status
        logical, intent(out) :: exited
        !> Gyrostep reduction factor in (0, 1] for convergence studies of the
        !> recorded canonical-momentum defect and the measured return mu; the
        !> production path omits it.
        real(dp), intent(in), optional :: h_scale
        real(dp) :: remaining, h, hfac, angle
        type(spectre_fo_local_t) :: particle
        integer :: nstep, map_status
        dt_used = 0.0_dp
        y = 0.0_dp
        owner = state%owner
        exited = .false.
        if (.not. state%active .or. dt < 0.0_dp .or. ro0_bar <= 0.0_dp) then
            status = SPECTRE_FO_INVALID_STATE
            return
        end if
        hfac = 1.0_dp
        if (present(h_scale)) then
            if (h_scale <= 0.0_dp .or. h_scale > 1.0_dp) then
                status = SPECTRE_FO_INVALID_STATE
                return
            end if
            hfac = h_scale
        end if
        ! Public states may be initialized from Cartesian particle data.
        ! Synchronize once per call; accepted logical steps cache x exactly.
        call synchronize_position(state, status)
        if (status /= SPECTRE_FO_OK) return
        remaining = dt
        angle = 0.0_dp
        call spectre_fo_field_local(state%u, state%owner, particle)
        do nstep = 1, MAX_SUBSTEPS
            if (remaining <= epsilon(1.0_dp)*max(dt, 1.0_dp)) then
                dt_used = dt
                status = SPECTRE_FO_OK
                return
            end if
            h = remaining
            if (state%h_rk > 0.0_dp) h = min(remaining, hfac*state%h_rk)
            h = min(h, hfac*RK_GYRO_ANGLE*ro0_bar/particle%bmod)
            call rk45_substep(state, h, ro0_bar, status, hfac*RK_GYRO_ANGLE, &
                particle)
            if (status /= SPECTRE_FO_OK) then
                if (status == SPECTRE_FO_LOSS) dt_used = dt_used + h
                return
            end if
            dt_used = dt_used + h
            remaining = remaining - h
            angle = angle + h*particle%bmod/ro0_bar
            if (angle < EXIT_CHECK_ANGLE .and. &
                remaining > epsilon(1.0_dp)*max(dt, 1.0_dp)) cycle
            angle = 0.0_dp
            call gc_map(state, particle, ro0_bar, y, owner, map_status, &
                exited=exited)
            if (map_status == SPECTRE_FO_OK) then
                state%last_y = y
                if (exited) then
                    state%active = .false.
                    status = SPECTRE_FO_OK
                    return
                end if
            end if
        end do
        status = SPECTRE_FO_STEP_LIMIT
    end subroutine spectre_fo_advance_until_exit
    subroutine spectre_fo_to_gc(state, ro0_bar, y, owner, status, mu_bar)
        type(spectre_fo_state_t), intent(inout) :: state
        real(dp), intent(in) :: ro0_bar
        real(dp), intent(out) :: y(5)
        integer, intent(out) :: owner, status
        real(dp), intent(out), optional :: mu_bar
        type(spectre_fo_local_t) :: particle

        y = 0.0_dp
        owner = state%owner
        if (present(mu_bar)) mu_bar = 0.0_dp
        if (.not. state%active .or. ro0_bar <= 0.0_dp) then
            status = SPECTRE_FO_INVALID_STATE
            return
        end if
        call synchronize_position(state, status)
        if (status /= SPECTRE_FO_OK) return
        call spectre_fo_field_local(state%u, state%owner, particle)
        call gc_map(state, particle, ro0_bar, y, owner, status, mu_bar)
    end subroutine spectre_fo_to_gc

    subroutine gc_map(state, particle, ro0_bar, y, owner, status, mu_bar, exited)
        type(spectre_fo_state_t), intent(in) :: state
        type(spectre_fo_local_t), intent(in) :: particle
        real(dp), intent(in) :: ro0_bar
        real(dp), intent(out) :: y(5)
        integer, intent(out) :: owner, status
        real(dp), intent(out), optional :: mu_bar
        logical, intent(out), optional :: exited
        type(spectre_fo_local_t) :: gc
        real(dp) :: b(3), b_gc(3), acov(3), e_zeta(3), bmod
        real(dp) :: vpar, vperp(3), rho_l(3), x_gc(3), speed, acov_gc(3)
        real(dp) :: e_zeta_gc(3), pzeta, h_zeta, perpendicular_energy, u_gc(3)
        real(dp) :: radial_larmor
        integer :: owner_gc

        if (present(exited)) exited = .false.
        y = 0.0_dp
        owner = state%owner
        if (present(mu_bar)) mu_bar = 0.0_dp
        if (.not. state%active .or. ro0_bar <= 0.0_dp) then
            status = SPECTRE_FO_INVALID_STATE
            return
        end if
        b = particle%b
        bmod = particle%bmod
        acov = particle%acov
        e_zeta = particle%e(:, 3)
        vpar = dot_product(state%v, b)
        pzeta = dot_product(state%v, e_zeta) + acov(3)/ro0_bar
        ! In an axisymmetric field toroidal canonical momentum is an exact
        ! physical invariant, including across a
        ! discontinuous radial sheet. Use the entry value when mapping back to
        ! GC variables so a short recovery excursion cannot create a hidden
        ! p_zeta kick.
        if (axisymmetric()) pzeta = state%entry_pzeta
        vperp = state%v - vpar*b
        rho_l = (ro0_bar/bmod)*cross(b, vperp)
        x_gc = particle%x - rho_l
        u_gc = state%u - matmul(rho_l, particle%dual)
        owner_gc = state%owner
        call locate_cart(x_gc, u_gc, owner_gc, status)
        if (status /= SPECTRE_FO_OK) return
        call spectre_fo_field_local(u_gc, owner_gc, gc)
        b_gc = gc%b
        bmod = gc%bmod
        acov_gc = gc%acov
        e_zeta_gc = gc%e(:, 3)
        speed = norm2(state%v)
        if (speed <= 0.0_dp) then
            status = SPECTRE_FO_INVALID_STATE
            return
        end if
        h_zeta = dot_product(b_gc, e_zeta_gc)
        if (axisymmetric() .and. abs(h_zeta) <= 100.0_dp*epsilon(1.0_dp)* &
            max(norm2(e_zeta_gc), 1.0_dp)) then
            status = SPECTRE_FO_INVALID_STATE
            return
        end if
        if (axisymmetric()) then
            ! Each volume has its own gauge (A vanishes on its inner boundary);
            ! in axisymmetry A_zeta is constant on an interface, so the entry
            ! momentum moves to the exit volume's gauge by the jumps of A_zeta
            ! across the interfaces in between.
            pzeta = pzeta + gauge_shift(state%entry_owner, owner_gc)/ro0_bar
            vpar = (pzeta - acov_gc(3)/ro0_bar)/h_zeta
        else
            vpar = dot_product(state%v, b_gc)
        end if
        perpendicular_energy = speed**2 - vpar**2
        if (perpendicular_energy < -100.0_dp*epsilon(1.0_dp)*speed**2) then
            status = SPECTRE_FO_INVALID_STATE
            return
        end if
        perpendicular_energy = max(perpendicular_energy, 0.0_dp)
        y(1:3) = u_gc
        y(4) = speed/sqrt(2.0_dp)
        y(5) = vpar/speed
        if (present(mu_bar)) mu_bar = perpendicular_energy/(2.0_dp*bmod)
        owner = owner_gc
        if (.not. present(exited)) return
        radial_larmor = max(abs(state%u(1) - y(1)), 1.0e-8_dp)
        if (state%exit_rho > 0.0_dp) then
            exited = abs(y(1) - real(state%iface, dp)) > state%exit_rho
        else if (state%iface == 0) then
            exited = abs(y(1) - state%entry_rho) > 4.0_dp*radial_larmor
        else
            exited = abs(y(1) - real(state%iface, dp)) > 4.0_dp*radial_larmor
        end if
        exited = exited .and. abs(particle%bmod - gc%bmod)/gc%bmod <= 0.05_dp
        exited = exited .and. norm2(particle%b - gc%b) <= 0.05_dp
    end subroutine gc_map
    real(dp) function gauge_shift(from, to) result(shift)
        !> Sum of the A_zeta jumps (in the field's acov units) across the
        !> interfaces from volume from to volume to, at theta = zeta = 0.
        integer, intent(in) :: from, to
        real(dp) :: x(3), b(3), bmod, a_low(3), a_high(3), e_zeta(3)
        integer :: k, step, status

        shift = 0.0_dp
        if (from == to .or. from < 1 .or. to < 1) return
        step = merge(1, -1, to > from)
        do k = min(from, to), max(from, to) - 1
            call logical_sample([real(k, dp), 0.0_dp, 0.0_dp], k, x, b, bmod, &
                a_low, e_zeta, status)
            call logical_sample([real(k, dp), 0.0_dp, 0.0_dp], k + 1, x, b, &
                bmod, a_high, e_zeta, status)
            shift = shift + real(step, dp)*(a_high(3) - a_low(3))
        end do
    end function gauge_shift
    logical function axisymmetric()
        axisymmetric = allocated(spectre_field)
        if (axisymmetric) axisymmetric = fo_field_axisymmetric
    end function axisymmetric
    subroutine spectre_fo_canonical_pzeta(state, ro0_bar, pzeta, status)
        type(spectre_fo_state_t), intent(inout) :: state
        real(dp), intent(in) :: ro0_bar
        real(dp), intent(out) :: pzeta
        integer, intent(out) :: status
        real(dp) :: b(3), bmod, acov(3), e_zeta(3)
        pzeta = 0.0_dp
        if (ro0_bar <= 0.0_dp) then
            status = SPECTRE_FO_INVALID_STATE
            return
        end if
        call cart_sample(state%x, state%u, state%owner, b, bmod, acov, &
            e_zeta, status)
        if (status /= SPECTRE_FO_OK) return
        pzeta = dot_product(state%v, e_zeta) + acov(3)/ro0_bar
    end subroutine spectre_fo_canonical_pzeta
    subroutine rk45_substep(state, h, ro0_bar, status, angle_limit, particle)
        !> Adaptive Dormand-Prince 5(4) step of x' = v, v' = (|B|/ro0) v x b. On
        !> the sharp field the stages continue the owner's spline polynomial
        !> beyond the volume edge. Stage cross-and-return
        !> trials are rejected; changed endpoints land on the first interface,
        !> where the owner switches. On input h is the trial step, on output the step
        !> taken.
        type(spectre_fo_state_t), intent(inout) :: state
        real(dp), intent(inout) :: h
        real(dp), intent(in) :: ro0_bar, angle_limit
        integer, intent(out) :: status
        type(spectre_fo_local_t), intent(inout) :: particle
        type(spectre_fo_local_t) :: final
        real(dp) :: x1(3), v1(3), u1(3), err, a, b, fa, fb, c, fc, iface_rho
        real(dp) :: trial_angle
        integer :: attempt, owner1, it, side, grazes, stage_sides, opposite_stage
        logical :: root_consistent

        grazes = 0
        do attempt = 1, 60
            call dp_trial(state, h, ro0_bar, x1, v1, u1, owner1, err, status, &
                particle, final, stage_sides, trial_angle)
            if (status /= SPECTRE_FO_OK) then
                ! A failed stage has no embedded error estimate.
                h = 0.5_dp*h
                cycle
            end if
            if (trial_angle > angle_limit*(1.0_dp + 4.0_dp*epsilon(1.0_dp))) then
                ! A zero-time owner change can increase |B|. Recheck the cap
                ! with this trial's initial field, including h_scale.
                h = h*angle_limit/trial_angle
                cycle
            end if
            if (err > 1.0_dp) then
                h = h*max(0.2_dp, min(0.5_dp, 0.9_dp*err**(-0.2_dp)))
                cycle
            end if
            if (stage_sides /= 0) then
                if (owner1 == state%owner .or. &
                    (owner1 > state%owner .and. btest(stage_sides, STAGE_IN)) .or. &
                    (owner1 < state%owner .and. btest(stage_sides, STAGE_OUT))) then
                    ! A stage detects a cross-and-return, or an earlier sheet
                    ! in the opposite direction. Never commit the continuation.
                    h = 0.5_dp*h
                    cycle
                end if
            end if
            state%h_rk = h*min(5.0_dp, 0.9_dp*max(err, 1.0e-10_dp)**(-0.2_dp))
            if (owner1 /= state%owner) then
                ! The first sheet of the initial owner, even if the endpoint
                ! lies beyond a thin intermediate volume.
                iface_rho = real(merge(state%owner, state%owner - 1, &
                    owner1 > state%owner), dp)
                opposite_stage = merge(STAGE_IN, STAGE_OUT, owner1 > state%owner)
                a = 0.0_dp
                fa = state%u(1) - iface_rho
                b = h
                fb = u1(1) - iface_rho
                if (abs(fa) <= OWNER_EPS) then
                    if (real(owner1 - state%owner, dp) &
                        *radial_velocity(state%u, state%v) > 0.0_dp) then
                        ! The outgoing normal velocity leaves the home side.
                        state%owner = merge(state%owner + 1, state%owner - 1, &
                            owner1 > state%owner)
                        if (state%owner > spectre_field%data%Mvol) then
                            state%owner = spectre_field%data%Mvol
                            h = 0.0_dp
                            status = SPECTRE_FO_LOSS
                            return
                        end if
                        grazes = grazes + 1
                        if (grazes >= 2) h = 0.5_dp*h
                        call spectre_fo_field_local(state%u, state%owner, particle)
                        cycle
                    end if
                    ! Starts on the sheet moving into home: exclude the root
                    ! at t=0 and bracket the later return at positive time.
                    a = h/32.0_dp
                    call dp_trial(state, a, ro0_bar, x1, v1, u1, owner1, err, &
                        status, particle, final, stage_sides)
                    if (status /= SPECTRE_FO_OK) then
                        h = 0.5_dp*h
                        cycle
                    end if
                    if (btest(stage_sides, opposite_stage)) then
                        h = 0.5_dp*h
                        cycle
                    end if
                    fa = u1(1) - iface_rho
                    if (fa*fb >= 0.0_dp) then
                        h = 0.5_dp*h
                        cycle
                    end if
                end if
                if (fb == fa) then
                    h = 0.5_dp*h
                    cycle
                end if
                root_consistent = .true.
                side = 0
                do it = 1, 100
                    c = (a*fb - b*fa)/(fb - fa)
                    call dp_trial(state, c, ro0_bar, x1, v1, u1, owner1, err, &
                        status, particle, final, stage_sides)
                    if (status /= SPECTRE_FO_OK) exit
                    if (btest(stage_sides, opposite_stage)) then
                        ! Shortened trials can expose an earlier opposite sheet
                        ! missed by the nominal stage grid. Retry from the start.
                        root_consistent = .false.
                        exit
                    end if
                    fc = u1(1) - iface_rho
                    if (abs(fc) <= 1.0e-12_dp) exit
                    if (fc*fb > 0.0_dp) then
                        b = c
                        fb = fc
                        if (side == -1) fa = 0.5_dp*fa
                        side = -1
                    else
                        a = c
                        fa = fc
                        if (side == 1) fb = 0.5_dp*fb
                        side = 1
                    end if
                end do
                if (status /= SPECTRE_FO_OK .or. .not. root_consistent) then
                    h = 0.5_dp*h
                    cycle
                end if
                if (c <= 1.0e-10_dp*h) then
                    ! The orbit grazes the interface at the step start: it
                    ! belongs to the other volume. Switch and redo the step.
                    ! Repeated grazing must reduce/retry, never bypass landing.
                    state%owner = merge(state%owner + 1, state%owner - 1, &
                        iface_rho >= real(state%owner, dp))
                    state%u(1) = iface_rho
                    if (state%owner > spectre_field%data%Mvol) then
                        state%owner = spectre_field%data%Mvol
                        h = 0.0_dp
                        status = SPECTRE_FO_LOSS
                        return
                    end if
                    grazes = grazes + 1
                    if (grazes >= 2) h = 0.5_dp*h
                    call spectre_fo_field_local(state%u, state%owner, particle)
                    cycle
                end if
                if (it > 100 .or. err > 1.0_dp) then
                    h = 0.5_dp*h
                    cycle
                end if
                h = c
                u1(1) = iface_rho
                owner1 = merge(state%owner + 1, state%owner - 1, &
                    iface_rho >= real(state%owner, dp))
            end if
            if (owner1 == state%owner) then
                ! The seventh RK stage is the accepted endpoint. Reuse its
                ! sample for GC reconstruction and the next first stage.
                particle = final
            else if (owner1 <= spectre_field%data%Mvol) then
                ! A landed sheet changes both the exact coordinate and owner.
                call spectre_fo_field_local(u1, owner1, particle)
            end if
            owner1 = min(owner1, spectre_field%data%Mvol)
            state%x = x1
            state%v = v1
            state%u = u1
            state%owner = owner1
            status = SPECTRE_FO_OK
            if (owner1 == spectre_field%data%Mvol .and. &
                u1(1) >= real(owner1, dp) - OWNER_EPS) status = SPECTRE_FO_LOSS
            return
        end do
        status = SPECTRE_FO_STEP_LIMIT
    end subroutine rk45_substep
    subroutine dp_trial(state, h, ro0_bar, x1, v1, u1, owner1, err, status, &
            initial, final, stage_sides, gyro_angle)
        !> Dormand-Prince on logical position and Cartesian velocity. The
        !> physical Lorentz equation is unchanged; du/dt = E^{-1} v replaces
        !> the iterative Cartesian-to-logical inversion at every RK stage.
        type(spectre_fo_state_t), intent(in) :: state
        real(dp), intent(in) :: h, ro0_bar
        real(dp), intent(out) :: x1(3), v1(3), u1(3), err
        integer, intent(out) :: owner1, status
        type(spectre_fo_local_t), intent(in) :: initial
        type(spectre_fo_local_t), intent(out) :: final
        integer, intent(out), optional :: stage_sides
        real(dp), intent(out), optional :: gyro_angle
        type(spectre_fo_local_t) :: loc
        real(dp) :: ku(3, 7), kv(3, 7), us(3), vs(3), qs(3), q0(3), w(7, 7)
        real(dp) :: rho_l, eu(3), ev(3), ex(3)
        integer :: i, j, owner_geo
        logical :: cartesian

        ! Cartesian position remains a regular chart near the polar axis.
        ! That small recovery region uses the same spline provider, with an
        ! inverse per stage; the regular logical path never performs one.
        cartesian = state%owner == 1 .and. state%u(1) < FO_AXIS_RADIUS
        q0 = state%u
        if (cartesian) q0 = state%x
        u1 = state%u
        if (present(stage_sides)) stage_sides = 0
        status = SPECTRE_FO_OK
        w = 0.0_dp
        w(2, 1) = 1.0_dp/5
        w(3, 1:2) = [3.0_dp/40, 9.0_dp/40]
        w(4, 1:3) = [44.0_dp/45, -56.0_dp/15, 32.0_dp/9]
        w(5, 1:4) = [19372.0_dp/6561, -25360.0_dp/2187, 64448.0_dp/6561, &
            -212.0_dp/729]
        w(6, 1:5) = [9017.0_dp/3168, -355.0_dp/33, 46732.0_dp/5247, &
            49.0_dp/176, -5103.0_dp/18656]
        w(7, 1:6) = [35.0_dp/384, 0.0_dp, 500.0_dp/1113, 125.0_dp/192, &
            -2187.0_dp/6784, 11.0_dp/84]
        do i = 1, 7
            qs = q0
            vs = state%v
            do j = 1, i - 1
                qs = qs + h*w(i, j)*ku(:, j)
                vs = vs + h*w(i, j)*kv(:, j)
            end do
            us = qs
            owner_geo = state%owner
            if (cartesian) then
                call locate_cart(qs, u1, owner_geo, status)
                if (status /= SPECTRE_FO_OK) return
                us = u1
            else
                owner_geo = logical_owner(us(1), state%owner)
            end if
            if (present(stage_sides)) then
                if (owner_geo < state%owner) &
                    stage_sides = ibset(stage_sides, STAGE_IN)
                if (owner_geo > state%owner) &
                    stage_sides = ibset(stage_sides, STAGE_OUT)
            end if
            if (us(1) <= 0.0_dp .and. .not. cartesian) then
                status = SPECTRE_FO_LOCATE_FAIL
                return
            end if
            ! Continue only the initial owner's polynomial for trial stages;
            ! accepted sheet crossings are landed and then change owner.
            if (i == 1 .and. .not. cartesian) then
                loc = initial
            else
                call spectre_fo_field_local(us, state%owner, loc)
            end if
            if (i == 1) then
                rho_l = ro0_bar*norm2(vs)/loc%bmod
                if (present(gyro_angle)) gyro_angle = h*loc%bmod/ro0_bar
            end if
            ku(:, i) = vs
            if (.not. cartesian) ku(:, i) = matmul(vs, loc%dual)
            kv(:, i) = (loc%bmod/ro0_bar)*cross(vs, loc%b)
        end do
        u1 = us
        v1 = vs
        x1 = loc%x
        final = loc
        if (cartesian) x1 = qs
        owner1 = logical_owner(u1(1), state%owner)
        eu = h*(71.0_dp/57600*ku(:, 1) - 71.0_dp/16695*ku(:, 3) + &
            71.0_dp/1920*ku(:, 4) - 17253.0_dp/339200*ku(:, 5) + &
            22.0_dp/525*ku(:, 6) - 1.0_dp/40*ku(:, 7))
        ev = h*(71.0_dp/57600*kv(:, 1) - 71.0_dp/16695*kv(:, 3) + &
            71.0_dp/1920*kv(:, 4) - 17253.0_dp/339200*kv(:, 5) + &
            22.0_dp/525*kv(:, 6) - 1.0_dp/40*kv(:, 7))
        ex = eu
        if (.not. cartesian) ex = matmul(loc%e, eu)
        err = max(norm2(ex)/rho_l, norm2(ev)/norm2(state%v))/fo_tolerance
    end subroutine dp_trial

    integer function logical_owner(rho, home) result(owner)
        real(dp), intent(in) :: rho
        integer, intent(in) :: home

        owner = home
        if (rho < real(home - 1, dp) - OWNER_EPS) &
            owner = max(1, int(rho) + 1)
        if (rho > real(home, dp) + OWNER_EPS) &
            owner = min(spectre_field%data%Mvol + 1, int(rho) + 1)
    end function logical_owner

    subroutine synchronize_position(state, status)
        type(spectre_fo_state_t), intent(inout) :: state
        integer, intent(out) :: status
        type(spectre_fo_local_t) :: loc

        call spectre_fo_field_local(state%u, state%owner, loc)
        status = SPECTRE_FO_OK
        if (norm2(state%x - loc%x) > 1.0e-12_dp*max(norm2(loc%x), 1.0_dp)) &
            call locate_cart(state%x, state%u, state%owner, status)
    end subroutine synchronize_position

    subroutine solve_gyrophase(x_gc, u_gc, owner, b, e1, e2, vpar, vperp, &
            target, ro0_bar, alpha, status)
        real(dp), intent(in) :: x_gc(3), u_gc(3), b(3), e1(3), e2(3)
        real(dp), intent(in) :: vpar, vperp, target, ro0_bar
        integer, intent(in) :: owner
        real(dp), intent(out) :: alpha
        integer, intent(out) :: status
        real(dp) :: a0, a1, f0, f1, amid, fmid, scale, best_residual, best_alpha
        integer :: i, ierr
        logical :: have_previous
        if (vperp <= sqrt(epsilon(1.0_dp))) then
            alpha = 0.0_dp
            call phase_residual(alpha, x_gc, u_gc, owner, b, e1, e2, vpar, &
                vperp, target, ro0_bar, f0, ierr)
            scale = max(abs(target), 1.0_dp)
            if (ierr /= SPECTRE_FO_OK) then
                status = ierr
            else if (abs(f0) <= 1.0e-11_dp*scale) then
                status = SPECTRE_FO_OK
            else
                status = SPECTRE_FO_GYROPHASE_FAIL
            end if
            return
        end if
        have_previous = .false.
        best_residual = huge(1.0_dp)
        best_alpha = 0.0_dp
        do i = 0, NPHASE
            a1 = 8.0_dp*atan(1.0_dp)*real(i, dp)/real(NPHASE, dp)
            call phase_residual(a1, x_gc, u_gc, owner, b, e1, e2, vpar, &
                vperp, target, ro0_bar, f1, ierr)
            if (ierr /= SPECTRE_FO_OK) then
                have_previous = .false.
                cycle
            end if
            if (abs(f1) < best_residual) then
                best_residual = abs(f1)
                best_alpha = a1
            end if
            if (f1 == 0.0_dp) then
                alpha = a1
                status = SPECTRE_FO_OK
                return
            end if
            if (have_previous) then
                if (f0*f1 <= 0.0_dp) exit
            end if
            a0 = a1
            f0 = f1
            have_previous = .true.
        end do
        if (i > NPHASE) then
            ! The gyrophase variation of p_zeta, proportional to v_perp, does
            ! not reach the guiding-centre value (strongly passing markers):
            ! take the closest gyrophase only with a bounded momentum defect.
            ! Entry records the actual signed defect; no order claim is made.
            if (best_residual > ENTRY_MOMENTUM_TOL*max(abs(target), 1.0_dp)) then
                status = SPECTRE_FO_GYROPHASE_FAIL
                return
            end if
            alpha = best_alpha
            status = SPECTRE_FO_OK
            return
        end if
        do i = 1, 60
            amid = 0.5_dp*(a0 + a1)
            call phase_residual(amid, x_gc, u_gc, owner, b, e1, e2, vpar, &
                vperp, target, ro0_bar, fmid, ierr)
            if (ierr /= SPECTRE_FO_OK) then
                status = ierr
                return
            end if
            if (f0*fmid <= 0.0_dp) then
                a1 = amid
            else
                a0 = amid
                f0 = fmid
            end if
        end do
        alpha = 0.5_dp*(a0 + a1)
        status = SPECTRE_FO_OK
    end subroutine solve_gyrophase
    subroutine phase_residual(alpha, x_gc, u_gc, owner, b_gc, e1, e2, vpar, &
            vperp, target, ro0_bar, residual, status)
        real(dp), intent(in) :: alpha, x_gc(3), u_gc(3), b_gc(3), e1(3), e2(3)
        real(dp), intent(in) :: vpar, vperp, target, ro0_bar
        integer, intent(in) :: owner
        real(dp), intent(out) :: residual
        integer, intent(out) :: status
        real(dp) :: eperp(3), x(3), u(3), b(3), bmod, acov(3), e_zeta(3), v(3)
        integer :: trial_owner
        eperp = cos(alpha)*e1 + sin(alpha)*e2
        call logical_sample(u_gc, owner, x, b, bmod, acov, e_zeta, status)
        if (status /= SPECTRE_FO_OK) return
        x = x_gc + (ro0_bar/bmod)*cross(b_gc, vperp*eperp)
        u = u_gc
        trial_owner = owner
        call cart_sample(x, u, trial_owner, b, bmod, acov, e_zeta, status)
        if (status /= SPECTRE_FO_OK) return
        v = vpar*b_gc + vperp*eperp
        residual = dot_product(v, e_zeta) &
            + (acov(3) - gauge_shift(owner, trial_owner))/ro0_bar - target
    end subroutine phase_residual
    subroutine cart_sample(x, u, owner, b, bmod, acov, e_zeta, status)
        real(dp), intent(in) :: x(3)
        real(dp), intent(inout) :: u(3)
        integer, intent(inout) :: owner
        real(dp), intent(out) :: b(3), bmod, acov(3), e_zeta(3)
        integer, intent(out) :: status
        real(dp) :: x_at_u(3)
        call locate_cart(x, u, owner, status)
        if (status /= SPECTRE_FO_OK) return
        call logical_sample(u, owner, x_at_u, b, bmod, acov, e_zeta, status)
    end subroutine cart_sample
    subroutine logical_sample(u, owner, x, b, bmod, acov, e_zeta, status)
        real(dp), intent(in) :: u(3)
        integer, intent(in) :: owner
        real(dp), intent(out) :: x(3), b(3), bmod, acov(3), e_zeta(3)
        integer, intent(out) :: status
        type(spectre_fo_local_t) :: loc

        call spectre_fo_field_local(u, owner, loc)
        x = loc%x
        b = loc%b
        bmod = loc%bmod
        acov = loc%acov
        e_zeta = loc%e(:, 3)
        status = SPECTRE_FO_OK
    end subroutine logical_sample

    subroutine locate_cart(x, u, owner, status)
        real(dp), intent(in) :: x(3)
        real(dp), intent(inout) :: u(3)
        integer, intent(inout) :: owner
        integer, intent(out) :: status
        integer :: ierr

        call spectre_fo_field_locate(x, u, owner, ierr)
        status = SPECTRE_FO_OK
        if (ierr /= 0) status = SPECTRE_FO_LOCATE_FAIL
    end subroutine locate_cart

    real(dp) function radial_velocity(u, v) result(vrho)
        real(dp), intent(in) :: u(3), v(3)
        type(spectre_fo_local_t) :: loc

        call spectre_fo_field_local(u, min(int(u(1)) + 1, &
            spectre_field%data%Mvol), loc)
        vrho = dot_product(v, loc%dual(:, 1))
    end function radial_velocity

    pure subroutine perpendicular_frame(b, e1, e2)
        real(dp), intent(in) :: b(3)
        real(dp), intent(out) :: e1(3), e2(3)
        real(dp) :: axis(3)

        if (abs(b(3)) < 0.9_dp) then
            axis = [0.0_dp, 0.0_dp, 1.0_dp]
        else
            axis = [1.0_dp, 0.0_dp, 0.0_dp]
        end if
        e1 = axis - dot_product(axis, b)*b
        e1 = e1/norm2(e1)
        e2 = cross(b, e1)
    end subroutine perpendicular_frame

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
            a(1)*b(2) - a(2)*b(1)]
    end function cross
end module spectre_fo_hybrid
