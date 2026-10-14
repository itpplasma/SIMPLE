module spectre_fo_hybrid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use spectre_fo_gauss2, only: gauss2_prepare, gauss2_encode, gauss2_decode, &
        gauss2_step, gauss2_gauge_zeta, gauss2_entry_step
    use libneo_coordinates, only: spectre_coordinate_system_t
    use magfie_sub, only: spectre_field, TESLA_TO_GAUSS, M_TO_CM
    implicit none
    private
    integer, parameter, public :: SPECTRE_FO_OK = 0
    integer, parameter, public :: SPECTRE_FO_UNINITIALIZED = 1
    integer, parameter, public :: SPECTRE_FO_LOCATE_FAIL = 2
    integer, parameter, public :: SPECTRE_FO_GYROPHASE_FAIL = 3
    integer, parameter, public :: SPECTRE_FO_INVALID_STATE = 4
    integer, parameter, public :: SPECTRE_FO_STEP_LIMIT = 5
    integer, parameter, public :: SPECTRE_FO_LOSS = 6
    real(dp), parameter :: A_SCALE = TESLA_TO_GAUSS*M_TO_CM**2
    real(dp), parameter :: OWNER_EPS = 1.0e-12_dp
    integer, parameter :: STAGE_IN = 0, STAGE_OUT = 1
    !> RK45 relative tolerance, configured before parallel tracing.
    integer, parameter, public :: FO_RK45 = 1, FO_GAUSS2 = 2
    integer :: fo_integrator = FO_RK45
    real(dp) :: fo_tolerance = 1.0e-9_dp
    ! At most 0.1 rad (1/63 of a gyro-period). A cross-and-return hidden
    ! between stage samples has circular-orbit depth <= rho_L*(1-cos(0.05)),
    ! about 1.25e-3 rho_L; tolerance alone cannot bound a missed sheet force.
    real(dp), parameter :: RK_GYRO_ANGLE = 0.1_dp
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
        ! All Gauss parameters are fixed at FO entry.
        real(dp) :: z_gauss(6) = 0.0_dp
        real(dp) :: h_tau = 0.0_dp, energy = 0.0_dp
        real(dp) :: layer_rho = 0.0_dp, refine_length = 0.0_dp
        real(dp) :: g_min = 1.0_dp, radial_scale = 0.0_dp
        real(dp) :: gc_check_time = 0.0_dp
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
    subroutine spectre_fo_configure(tolerance, integrator)
        real(dp), intent(in) :: tolerance
        integer, intent(in), optional :: integrator

        if (tolerance <= 0.0_dp) &
            error stop 'spectre_fo_configure: tolerance must be positive'
        fo_tolerance = tolerance
        fo_integrator = FO_RK45
        if (present(integrator)) then
            if (integrator /= FO_RK45 .and. integrator /= FO_GAUSS2) &
                error stop 'spectre_fo_configure: unknown integrator'
            fo_integrator = integrator
        end if
        if (fo_integrator == FO_GAUSS2) call gauss2_prepare
    end subroutine spectre_fo_configure
    subroutine spectre_fo_enter(state, y, owner, iface, ro0_bar, status)
        type(spectre_fo_state_t), intent(out) :: state
        real(dp), intent(in) :: y(5), ro0_bar
        integer, intent(in) :: owner, iface
        integer, intent(out) :: status
        real(dp) :: x_gc(3), b(3), e1(3), e2(3), eperp(3), e_zeta(3)
        real(dp) :: acov(3), bmod, vpar, vperp, target, alpha, radial_speed, h_entry
        logical :: ok
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
        call locate_cart(state%x, state%u, state%owner, status)
        if (status /= SPECTRE_FO_OK) return
        if (axisymmetric()) then
            call phase_residual(alpha, x_gc, y(1:3), owner, b, e1, e2, vpar, &
                vperp, target, ro0_bar, state%entry_pzeta_defect, status)
            if (status /= SPECTRE_FO_OK) return
        end if
        if (fo_integrator == FO_GAUSS2) then
            radial_speed = sqrt(2.0_dp)*y(4)/M_TO_CM*radial_gradient(y(1:3))
            state%layer_rho = 0.05_dp*ro0_bar/bmod*radial_speed
            if (state%layer_rho <= 0.0_dp .or. state%layer_rho >= 1.0_dp) then
                status = SPECTRE_FO_INVALID_STATE
                return
            end if
            h_entry = RK_GYRO_ANGLE*ro0_bar/bmod
            state%refine_length = max(16.0_dp*h_entry*radial_speed, &
                16.0_dp*state%layer_rho)
            state%g_min = min(0.1_dp, state%layer_rho/ &
                (256.0_dp*h_entry*radial_speed))
            state%energy = 0.5_dp*dot_product(state%v, state%v)
            state%radial_scale = (state%refine_length/ &
                (16.0_dp*h_entry*sqrt(2.0_dp*state%energy)))**2
            ! Prescribed safety refinement, chosen before any stage solve.
            ! Retain the Hamiltonian's metric reference from the base entry h.
            state%h_tau = gauss2_entry_step(h_entry, state%layer_rho, &
                state%refine_length, state%g_min)
            state%z_gauss(1:3) = state%u
            call gauss2_encode(state%x, state%v, state%layer_rho, ro0_bar, &
                state%z_gauss, ok)
            if (.not. ok) then
                status = SPECTRE_FO_LOCATE_FAIL
                return
            end if
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
        rho_r = ro0_bar*vperp/bmod/M_TO_CM*radial_gradient(y(1:3))
    end subroutine spectre_fo_radial_larmor
    real(dp) function radial_gradient(u) result(grad)
        !> |grad rho_g| per metre at the logical point u.
        real(dp), intent(in) :: u(3)
        real(dp) :: g(3, 3), ginv(3, 3), sqrtg

        call spectre_field%coords%metric_tensor(u, g, ginv, sqrtg)
        grad = sqrt(ginv(1, 1))
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
        real(dp) :: remaining, h, bmod, b(3), acov(3), e_zeta(3), hfac
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
        if (fo_integrator == FO_GAUSS2) then
            call gauss_advance_until_exit(state, dt, ro0_bar, hfac, dt_used, &
                y, owner, exited, status)
            return
        end if
        remaining = dt
        do nstep = 1, MAX_SUBSTEPS
            if (remaining <= epsilon(1.0_dp)*max(dt, 1.0_dp)) then
                dt_used = dt
                status = SPECTRE_FO_OK
                return
            end if
            call cart_sample(state%x, state%u, state%owner, b, bmod, acov, &
                e_zeta, status)
            if (status /= SPECTRE_FO_OK) return
            h = remaining
            if (state%h_rk > 0.0_dp) h = min(remaining, hfac*state%h_rk)
            h = min(h, hfac*RK_GYRO_ANGLE*ro0_bar/bmod)
            call rk45_substep(state, h, ro0_bar, status, hfac*RK_GYRO_ANGLE)
            if (status /= SPECTRE_FO_OK) then
                if (status == SPECTRE_FO_LOSS) dt_used = dt_used + h
                return
            end if
            dt_used = dt_used + h
            remaining = remaining - h
            call spectre_fo_to_gc(state, ro0_bar, y, owner, map_status)
            if (map_status == SPECTRE_FO_OK) then
                state%last_y = y
                call check_exit(state, y, owner, exited, map_status)
                if (map_status /= SPECTRE_FO_OK) then
                    status = map_status
                    return
                end if
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
        real(dp) :: b(3), b_gc(3), acov(3), e_zeta(3), bmod
        real(dp) :: vpar, vperp(3), rho_l(3), x_gc(3), speed, acov_gc(3)
        real(dp) :: e_zeta_gc(3), pzeta, h_zeta, perpendicular_energy, u_gc(3)
        integer :: owner_gc
        y = 0.0_dp
        owner = state%owner
        if (present(mu_bar)) mu_bar = 0.0_dp
        if (.not. state%active .or. ro0_bar <= 0.0_dp) then
            status = SPECTRE_FO_INVALID_STATE
            return
        end if
        call cart_sample(state%x, state%u, state%owner, b, bmod, acov, &
            e_zeta, status)
        if (status /= SPECTRE_FO_OK) return
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
        x_gc = state%x - rho_l
        u_gc = state%u
        owner_gc = state%owner
        call locate_cart(x_gc, u_gc, owner_gc, status)
        if (status /= SPECTRE_FO_OK) return
        call logical_sample(u_gc, owner_gc, x_gc, b_gc, bmod, acov_gc, &
            e_zeta_gc, status)
        if (status /= SPECTRE_FO_OK) return
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
    end subroutine spectre_fo_to_gc
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
    subroutine check_exit(state, y_gc, owner_gc, exited, status)
        type(spectre_fo_state_t), intent(inout) :: state
        real(dp), intent(in) :: y_gc(5)
        integer, intent(in) :: owner_gc
        logical, intent(out) :: exited
        integer, intent(out) :: status
        real(dp) :: bp(3), bgc(3), bmod_p, bmod_gc, acov(3), e_zeta(3), x_gc(3)
        real(dp) :: radial_larmor
        call cart_sample(state%x, state%u, state%owner, bp, bmod_p, acov, &
            e_zeta, status)
        if (status /= SPECTRE_FO_OK) return
        call logical_sample(y_gc(1:3), owner_gc, x_gc, bgc, bmod_gc, acov, &
            e_zeta, status)
        if (status /= SPECTRE_FO_OK) return
        radial_larmor = max(abs(state%u(1) - y_gc(1)), 1.0e-8_dp)
        if (state%exit_rho > 0.0_dp) then
            exited = abs(y_gc(1) - real(state%iface, dp)) > state%exit_rho
        else if (state%iface == 0) then
            exited = abs(y_gc(1) - state%entry_rho) > 4.0_dp*radial_larmor
        else
            exited = abs(y_gc(1) - real(state%iface, dp)) > &
                4.0_dp*radial_larmor
        end if
        exited = exited .and. abs(bmod_p - bmod_gc)/bmod_gc <= 0.05_dp
        exited = exited .and. norm2(bp - bgc) <= 0.05_dp
    end subroutine check_exit
    logical function axisymmetric()
        axisymmetric = allocated(spectre_field)
        if (axisymmetric) axisymmetric = all(spectre_field%data%in == 0)
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
        if (fo_integrator == FO_GAUSS2) then
            pzeta = state%z_gauss(6) - gauss2_gauge_zeta( &
                state%z_gauss(1:3), state%owner, ro0_bar)
            status = SPECTRE_FO_OK
            return
        end if
        call cart_sample(state%x, state%u, state%owner, b, bmod, acov, &
            e_zeta, status)
        if (status /= SPECTRE_FO_OK) return
        pzeta = dot_product(state%v, e_zeta) + acov(3)/ro0_bar
    end subroutine spectre_fo_canonical_pzeta
    subroutine gauss_advance_until_exit(state, dt, ro0_bar, hfac, dt_used, &
            y, owner, exited, status)
        type(spectre_fo_state_t), intent(inout) :: state
        real(dp), intent(in) :: dt, ro0_bar, hfac
        real(dp), intent(out) :: dt_used, y(5)
        integer, intent(out) :: owner, status
        logical, intent(out) :: exited
        real(dp) :: z1(6), clock, x1(3), v1(3)
        integer :: nstep, map_status
        logical :: ok

        dt_used = 0.0_dp
        y = 0.0_dp
        owner = state%owner
        exited = .false.
        status = SPECTRE_FO_INVALID_STATE
        if (state%h_tau <= 0.0_dp .or. state%energy <= 0.0_dp .or. &
            state%layer_rho <= 0.0_dp .or. state%layer_rho >= 1.0_dp) return
        status = SPECTRE_FO_OK
        do nstep = 1, 10000000
            if (dt_used >= dt) return
            call gauss2_step(state%z_gauss, state%layer_rho, ro0_bar, state%energy, &
                state%refine_length, state%g_min, hfac*state%h_tau, z1, clock, ok, &
                radial_scale=state%radial_scale)
            if (.not. ok) then
                ! Retain the last accepted canonical state and refresh its
                ! Cartesian representation before propagating the failure.
                call gauss2_decode(state%z_gauss, state%layer_rho, ro0_bar, &
                    state%x, state%v)
                status = SPECTRE_FO_STEP_LIMIT
                return
            end if
            state%z_gauss = z1
            state%u = z1(1:3)
            state%owner = min(max(int(z1(1)) + 1, 1), spectre_field%data%Mvol)
            owner = state%owner
            dt_used = dt_used + clock
            state%gc_check_time = state%gc_check_time + clock
            if (state%gc_check_time < 0.25_dp*hfac*state%h_tau .and. &
                dt_used < dt .and. z1(1) < real(spectre_field%data%Mvol, dp)) cycle
            call gauss2_decode(z1, state%layer_rho, ro0_bar, x1, v1)
            if (.not. all(ieee_is_finite(v1))) then
                status = SPECTRE_FO_INVALID_STATE
                return
            end if
            state%x = x1
            state%v = v1
            state%gc_check_time = 0.0_dp
            if (z1(1) >= real(spectre_field%data%Mvol, dp)) then
                status = SPECTRE_FO_LOSS
                return
            end if
            call spectre_fo_to_gc(state, ro0_bar, y, owner, map_status)
            if (map_status /= SPECTRE_FO_OK) cycle
            state%last_y = y
            call check_exit(state, y, owner, exited, map_status)
            if (map_status /= SPECTRE_FO_OK) then
                status = map_status
                return
            end if
            if (exited) then
                state%active = .false.
                return
            end if
        end do
        status = SPECTRE_FO_STEP_LIMIT
    end subroutine gauss_advance_until_exit

    subroutine rk45_substep(state, h, ro0_bar, status, angle_limit)
        !> Adaptive Dormand-Prince 5(4) step of x' = v, v' = (|B|/ro0) v x b. On
        !> the sharp field the stages use the owner's field, held at
        !> its interface value beyond the volume edge. Stage cross-and-return
        !> trials are rejected; changed endpoints land on the first interface,
        !> where the owner switches. On input h is the trial step, on output the step
        !> taken.
        type(spectre_fo_state_t), intent(inout) :: state
        real(dp), intent(inout) :: h
        real(dp), intent(in) :: ro0_bar, angle_limit
        integer, intent(out) :: status
        real(dp) :: x1(3), v1(3), u1(3), err, a, b, fa, fb, c, fc, iface_rho
        real(dp) :: trial_angle
        integer :: attempt, owner1, it, side, grazes, stage_sides, opposite_stage
        logical :: root_consistent

        grazes = 0
        do attempt = 1, 60
            call dp_trial(state, h, ro0_bar, x1, v1, u1, owner1, err, status, &
                stage_sides, trial_angle)
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
                        grazes = grazes + 1
                        if (grazes >= 2) h = 0.5_dp*h
                        cycle
                    end if
                    ! Starts on the sheet moving into home: exclude the root
                    ! at t=0 and bracket the later return at positive time.
                    a = h/32.0_dp
                    call dp_trial(state, a, ro0_bar, x1, v1, u1, owner1, err, &
                        status, stage_sides)
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
                        status, stage_sides)
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
                    grazes = grazes + 1
                    if (grazes >= 2) h = 0.5_dp*h
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
            stage_sides, gyro_angle)
        !> One Dormand-Prince 5(4) step from the state; err is the scaled error
        !> (gyroradius for x, |v| for v) over fo_tolerance.
        type(spectre_fo_state_t), intent(in) :: state
        real(dp), intent(in) :: h, ro0_bar
        real(dp), intent(out) :: x1(3), v1(3), u1(3), err
        integer, intent(out) :: owner1, status
        integer, intent(out), optional :: stage_sides
        real(dp), intent(out), optional :: gyro_angle
        real(dp) :: kx(3, 7), kv(3, 7), xs(3), vs(3), w(7, 7), bvec(3), bmod
        real(dp) :: rho_l, ex(3), ev(3)
        integer :: i, j, owner_geo

        if (present(stage_sides)) stage_sides = 0
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
        u1 = state%u
        owner1 = state%owner
        do i = 1, 7
            xs = state%x
            vs = state%v
            do j = 1, i - 1
                xs = xs + h*w(i, j)*kx(:, j)
                vs = vs + h*w(i, j)*kv(:, j)
            end do
            call stage_field(state, xs, u1, bvec, bmod, owner_geo, status)
            if (status /= SPECTRE_FO_OK) return
            if (present(stage_sides)) then
                if (owner_geo < state%owner) &
                    stage_sides = ibset(stage_sides, STAGE_IN)
                if (owner_geo > state%owner) &
                    stage_sides = ibset(stage_sides, STAGE_OUT)
            end if
            if (i == 1) then
                rho_l = ro0_bar*norm2(vs)/bmod
                if (present(gyro_angle)) gyro_angle = h*bmod/ro0_bar
            end if
            kx(:, i) = vs
            kv(:, i) = (bmod/ro0_bar)*cross(vs, bvec)
        end do
        x1 = xs
        v1 = vs
        ex = h*(71.0_dp/57600*kx(:, 1) - 71.0_dp/16695*kx(:, 3) + &
            71.0_dp/1920*kx(:, 4) - 17253.0_dp/339200*kx(:, 5) + &
            22.0_dp/525*kx(:, 6) - 1.0_dp/40*kx(:, 7))
        ev = h*(71.0_dp/57600*kv(:, 1) - 71.0_dp/16695*kv(:, 3) + &
            71.0_dp/1920*kv(:, 4) - 17253.0_dp/339200*kv(:, 5) + &
            22.0_dp/525*kv(:, 6) - 1.0_dp/40*kv(:, 7))
        err = max(norm2(ex)/rho_l, norm2(ev)/norm2(state%v))/fo_tolerance
        call locate_cart(x1, u1, owner1, status)
        if (status /= SPECTRE_FO_OK) return
    end subroutine dp_trial
    subroutine stage_field(state, x, u, bvec, bmod, owner_geo, status)
        !> Field at a stage point: the owner's sharp field, clamped at its edge.
        type(spectre_fo_state_t), intent(in) :: state
        real(dp), intent(in) :: x(3)
        real(dp), intent(inout) :: u(3)
        real(dp), intent(out) :: bvec(3), bmod
        integer, intent(out) :: owner_geo, status
        real(dp) :: acov(3), e_zeta(3), xs(3)

        owner_geo = state%owner
        call locate_cart(x, u, owner_geo, status)
        if (status /= SPECTRE_FO_OK) return
        call logical_sample(u, state%owner, xs, bvec, bmod, acov, e_zeta, status)
        if (status /= SPECTRE_FO_OK) return
    end subroutine stage_field
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
        real(dp) :: ue(3), hcov(3), sqgb(3), sqrtg, g(3, 3), ginv(3, 3)
        real(dp) :: e_cov(3, 3), bcart(3)
        integer :: i

        if (.not. allocated(spectre_field)) then
            status = SPECTRE_FO_UNINITIALIZED
            return
        end if
        ue = u
        ue(1) = max(ue(1), real(owner - 1, dp))
        ue(1) = min(ue(1), real(owner, dp))
        if (owner < spectre_field%data%Mvol) then
            if (ue(1) == real(owner, dp)) ue(1) = ue(1) - OWNER_EPS
        end if
        call spectre_field%coords%evaluate_cart(u, x)
        call spectre_field%evaluate(ue, acov, hcov, bmod, sqgb)
        call spectre_field%coords%covariant_basis(ue, e_cov)
        call spectre_field%coords%metric_tensor(ue, g, ginv, sqrtg)
        bcart = 0.0_dp
        do i = 1, 3
            bcart = bcart + e_cov(:, i)*sqgb(i)/sqrtg
        end do
        b = bcart/max(bmod, tiny(1.0_dp))
        x = x*M_TO_CM
        e_zeta = e_cov(:, 3)*M_TO_CM
        bmod = bmod*TESLA_TO_GAUSS
        acov = acov*A_SCALE
        status = SPECTRE_FO_OK
    end subroutine logical_sample

    subroutine locate_cart(x, u, owner, status)
        real(dp), intent(in) :: x(3)
        real(dp), intent(inout) :: u(3)
        integer, intent(inout) :: owner
        integer, intent(out) :: status
        real(dp) :: xcyl(3), trial(3)
        integer :: candidate, ierr

        if (.not. allocated(spectre_field)) then
            status = SPECTRE_FO_UNINITIALIZED
            return
        end if
        xcyl = [sqrt(x(1)**2 + x(2)**2)/M_TO_CM, atan2(x(2), x(1)), &
            x(3)/M_TO_CM]
        select type (coords => spectre_field%coords)
        type is (spectre_coordinate_system_t)
            trial = u
            call coords%from_cyl_warm(xcyl, trial, owner, ierr)
            if (ierr == 0) then
                u = trial
                call refine_sheet_location(x, u, owner)
                status = SPECTRE_FO_OK
                return
            end if
            do candidate = max(1, owner - 1), min(coords%Mvol, owner + 1)
                if (candidate == owner) cycle
                trial = u
                trial(1) = min(max(trial(1), real(candidate - 1, dp)), &
                    real(candidate, dp))
                call coords%from_cyl_warm(xcyl, trial, candidate, ierr)
                if (ierr == 0) then
                    u = trial
                    owner = candidate
                    call refine_sheet_location(x, u, owner)
                    status = SPECTRE_FO_OK
                    return
                end if
            end do
            trial = u
            call coords%from_cyl(xcyl, trial, ierr)
            if (ierr == 0) then
                candidate = min(int(trial(1)) + 1, coords%Mvol)
                if (abs(trial(1) - real(nint(trial(1)), dp)) > OWNER_EPS) &
                    owner = candidate
                u = trial
                call refine_sheet_location(x, u, owner)
                status = SPECTRE_FO_OK
                return
            end if
        class default
            status = SPECTRE_FO_INVALID_STATE
            return
        end select
        status = SPECTRE_FO_LOCATE_FAIL
    end subroutine locate_cart

    subroutine refine_sheet_location(x, u, owner)
        ! The warm inverse snaps within 5e-10 rho of a sheet. Recover its
        ! signed distance for event bracketing; otherwise fa=fb=0 is possible.
        real(dp), intent(in) :: x(3)
        real(dp), intent(inout) :: u(3)
        integer, intent(inout) :: owner
        real(dp) :: xu(3), e(3, 3), g(3, 3), gi(3, 3), sg
        integer :: k

        k = nint(u(1))
        if (k < 1 .or. k > spectre_field%data%Mvol) return
        if (abs(u(1) - real(k, dp)) > 1.0e-9_dp) return
        call spectre_field%coords%evaluate_cart(u, xu)
        call spectre_field%coords%covariant_basis(u, e)
        call spectre_field%coords%metric_tensor(u, g, gi, sg)
        u(1) = u(1) + dot_product(x/M_TO_CM - xu, matmul(e, gi(:, 1)))
        if (abs(u(1) - real(k, dp)) <= OWNER_EPS) then
            ! Preserve exact terminal coordinates and the landed owner within
            ! the event tolerance, including Cartesian roundoff at an edge.
            u(1) = real(k, dp)
        else
            owner = min(int(u(1)) + 1, spectre_field%data%Mvol)
        end if
    end subroutine refine_sheet_location

    real(dp) function radial_velocity(u, v) result(vrho)
        real(dp), intent(in) :: u(3), v(3)
        real(dp) :: e(3, 3), g(3, 3), gi(3, 3), sg

        call spectre_field%coords%covariant_basis(u, e)
        call spectre_field%coords%metric_tensor(u, g, gi, sg)
        vrho = dot_product(v/M_TO_CM, matmul(e, gi(:, 1)))
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
