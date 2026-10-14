module spectre_fo_hybrid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use libneo_coordinates, only: spectre_coordinate_system_t
    use magfie_sub, only: spectre_field, TESLA_TO_GAUSS, M_TO_CM
    use spectre_fo_gauss2, only: gauss2_prepare, gauss2_encode, gauss2_decode, &
        gauss2_step, gauss2_gauge_zeta
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
    !> Boris step bound h |B| / (2 ro0) <= rotation_parameter, and the width of
    !> the smooth layer that replaces each interior interface for the full
    !> orbit, in gyroradii (0: the sharp per-volume field). Set by
    !> spectre_fo_configure; read-only during tracing.
    real(dp) :: rotation_parameter = 0.05_dp
    real(dp) :: layer_larmor = 0.0_dp
    !> Full-orbit integrator and the relative tolerance of FO_RK45.
    integer, parameter, public :: FO_BORIS = 1, FO_RK45 = 2, FO_GAUSS2 = 3
    integer :: fo_integrator = FO_BORIS
    real(dp) :: fo_tolerance = 1.0e-9_dp
    integer, parameter :: MAX_SUBSTEPS = 100000
    integer, parameter :: NPHASE = 64
    type, public :: spectre_fo_state_t
        real(dp) :: x(3) = 0.0_dp
        real(dp) :: v(3) = 0.0_dp
        real(dp) :: u(3) = 0.0_dp
        real(dp) :: last_y(5) = 0.0_dp
        real(dp) :: entry_rho = 0.0_dp
        real(dp) :: entry_pzeta = 0.0_dp
        !> Next step of FO_RK45 (adapted between substeps).
        real(dp) :: h_rk = 0.0_dp
        !> Canonical Gauss state and frozen parameters of its FO excursion.
        real(dp) :: z_gauss(6) = 0.0_dp
        real(dp) :: h_tau = 0.0_dp, energy = 0.0_dp
        real(dp) :: refine_length = 0.0_dp, g_min = 1.0_dp
        real(dp) :: gc_check_time = 0.0_dp
        !> Volume whose gauge entry_pzeta is expressed in.
        integer :: entry_owner = 0
        !> Smooth-layer width in rho_g, fixed at entry from the gyroradius.
        real(dp) :: layer_rho = 0.0_dp
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
    subroutine spectre_fo_configure(rotation, layer, integrator, tolerance)
        real(dp), intent(in) :: rotation, layer
        integer, intent(in), optional :: integrator
        real(dp), intent(in), optional :: tolerance

        if (rotation <= 0.0_dp .or. layer < 0.0_dp) &
            error stop 'spectre_fo_configure: invalid rotation or layer width'
        rotation_parameter = rotation
        layer_larmor = layer
        if (present(integrator)) then
            if (integrator /= FO_BORIS .and. integrator /= FO_RK45 .and. &
                integrator /= FO_GAUSS2) &
                error stop 'spectre_fo_configure: unknown integrator'
            fo_integrator = integrator
        end if
        if (present(tolerance)) then
            if (tolerance <= 0.0_dp) &
                error stop 'spectre_fo_configure: tolerance must be positive'
            fo_tolerance = tolerance
        end if
        if (fo_integrator == FO_GAUSS2) call gauss2_prepare
    end subroutine spectre_fo_configure
    subroutine spectre_fo_enter(state, y, owner, iface, ro0_bar, status)
        type(spectre_fo_state_t), intent(out) :: state
        real(dp), intent(in) :: y(5), ro0_bar
        integer, intent(in) :: owner, iface
        integer, intent(out) :: status
        real(dp) :: x_gc(3), b(3), e1(3), e2(3), eperp(3), e_zeta(3)
        real(dp) :: acov(3), bmod, vpar, vperp, target, alpha, radial_speed
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
        state%layer_rho = layer_larmor*(ro0_bar*vperp/bmod)/M_TO_CM &
                          *radial_gradient(y(1:3))
        if (fo_integrator == FO_GAUSS2) then
            radial_speed = sqrt(2.0_dp)*y(4)/M_TO_CM*radial_gradient(y(1:3))
            ! Freeze all parameters at entry: state-dependent steps would destroy
            ! the symplectic discrete map. Default layer is rho/20 for Gauss.
            state%layer_rho = merge(layer_larmor, 0.05_dp, layer_larmor > 0.0_dp)* &
                ro0_bar/bmod*radial_speed
            if (state%layer_rho >= 1.0_dp) then
                ! The compact layers around consecutive integer interfaces
                ! must not overlap: nearest-interface dispatch assumes this.
                status = SPECTRE_FO_INVALID_STATE
                return
            end if
            state%h_tau = 2.0_dp*rotation_parameter*ro0_bar/bmod
            state%refine_length = max(16.0_dp*state%h_tau*radial_speed, &
                16.0_dp*state%layer_rho)
            state%g_min = min(0.1_dp, state%layer_rho/ &
                (256.0_dp*state%h_tau*radial_speed))
            state%energy = 0.5_dp*dot_product(state%v, state%v)
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
    subroutine layer_sample(layer_rho, u, owner, b, bmod, status)
        !> Field for the full orbit: inside the layer |rho_g - k| < w/2 around
        !> an interior interface k it blends the two volumes' fields (each held
        !> at its interface value beyond its edge) with a C2 quintic step; B.n
        !> vanishes on both sides, so the blend stays divergence-free to O(w).
        real(dp), intent(in) :: layer_rho, u(3)
        integer, intent(in) :: owner
        real(dp), intent(inout) :: b(3), bmod
        integer, intent(out) :: status
        real(dp) :: x(3), b_other(3), bmod_other, acov(3), e_zeta(3), xi, s
        real(dp) :: b_low(3), b_high(3)
        integer :: k, other

        status = SPECTRE_FO_OK
        if (layer_rho <= 0.0_dp) return
        k = nint(u(1))
        if (k < 1 .or. k >= spectre_field%data%Mvol) return
        if (abs(u(1) - real(k, dp)) >= 0.5_dp*layer_rho) return
        other = merge(k + 1, k, owner == k)
        call logical_sample(u, other, x, b_other, bmod_other, acov, e_zeta, status)
        if (status /= SPECTRE_FO_OK) return
        if (owner == k) then
            b_low = bmod*b
            b_high = bmod_other*b_other
        else
            b_low = bmod_other*b_other
            b_high = bmod*b
        end if
        xi = (u(1) - real(k, dp))/layer_rho + 0.5_dp
        s = xi**3*(10.0_dp - 15.0_dp*xi + 6.0_dp*xi*xi)
        b = (1.0_dp - s)*b_low + s*b_high
        bmod = norm2(b)
        b = b/bmod
    end subroutine layer_sample
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
        integer :: nstep, map_status, decode_status
        logical :: check_gc
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
        remaining = dt
        ! Refinement can require many fictitious steps per physical output
        ! interval. Keep the existing limits for the other two pushers.
        do nstep = 1, merge(10000000, MAX_SUBSTEPS, fo_integrator == FO_GAUSS2)
            if (remaining <= epsilon(1.0_dp)*max(dt, 1.0_dp)) then
                if (fo_integrator /= FO_GAUSS2) dt_used = dt
                status = SPECTRE_FO_OK
                return
            end if
            if (fo_integrator == FO_GAUSS2) then
                call gauss_substep(state, ro0_bar, hfac, h, status)
            else
                call cart_sample(state%x, state%u, state%owner, b, bmod, acov, &
                    e_zeta, status)
                if (status /= SPECTRE_FO_OK) return
                h = min(remaining, hfac*2.0_dp*rotation_parameter*ro0_bar/bmod)
                if (fo_integrator == FO_RK45) then
                    if (state%h_rk > 0.0_dp) h = min(remaining, hfac*state%h_rk)
                    call rk45_substep(state, h, ro0_bar, status)
                else
                    call boris_substep(state, h, ro0_bar, status)
                end if
            end if
            if (status /= SPECTRE_FO_OK) then
                if (fo_integrator == FO_GAUSS2) then
                    ! Refresh the last accepted state, including on a wall hit.
                    call gauss_refresh(state, ro0_bar, decode_status)
                end if
                return
            end if
            dt_used = dt_used + h
            remaining = remaining - h
            if (fo_integrator == FO_GAUSS2) then
                state%gc_check_time = state%gc_check_time + h
                check_gc = state%gc_check_time >= 0.25_dp*hfac*state%h_tau &
                    .or. remaining <= epsilon(1.0_dp)*max(dt, 1.0_dp)
                if (.not. check_gc) cycle
                call gauss_refresh(state, ro0_bar, status)
                if (status /= SPECTRE_FO_OK) return
                state%gc_check_time = 0.0_dp
            end if
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
        ! Standard Boris conserves speed exactly but not the canonical momentum
        ! associated with an ignorable toroidal angle. In an axisymmetric field
        ! that momentum is an exact physical invariant, including across a
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
    subroutine gauss_substep(state, ro0_bar, hfac, dt, status)
        type(spectre_fo_state_t), intent(inout) :: state
        real(dp), intent(in) :: ro0_bar, hfac
        real(dp), intent(out) :: dt
        integer, intent(out) :: status
        real(dp) :: z1(6)
        logical :: ok

        dt = 0.0_dp
        call gauss2_step(state%z_gauss, state%layer_rho, ro0_bar, state%energy, &
            state%refine_length, state%g_min, hfac*state%h_tau, z1, dt, ok)
        if (.not. ok) then
            status = SPECTRE_FO_STEP_LIMIT
            return
        end if
        if (z1(1) >= real(spectre_field%data%Mvol, dp)) then
            status = SPECTRE_FO_LOSS
            return
        end if
        state%z_gauss = z1
        ! Both embeddings agree on the interface itself and retain its side.
        ! The GC mapping below locates Cartesian x in the original chart.
        state%u = z1(1:3)
        state%owner = min(max(int(z1(1)) + 1, 1), spectre_field%data%Mvol)
        status = SPECTRE_FO_OK
    end subroutine gauss_substep

    subroutine gauss_refresh(state, ro0_bar, status)
        type(spectre_fo_state_t), intent(inout) :: state
        real(dp), intent(in) :: ro0_bar
        integer, intent(out) :: status
        real(dp) :: x1(3), v1(3)

        call gauss2_decode(state%z_gauss, state%layer_rho, ro0_bar, x1, v1)
        status = SPECTRE_FO_INVALID_STATE
        if (.not. all(ieee_is_finite(v1))) return
        state%x = x1
        state%v = v1
        status = SPECTRE_FO_OK
    end subroutine gauss_refresh

    subroutine rk45_substep(state, h, ro0_bar, status)
        !> Adaptive Dormand-Prince 5(4) step of x' = v, v' = (|B|/ro0) v x b. On
        !> the sharp field (no layer) the stages use the owner's field, held at
        !> its interface value beyond the volume edge, and a step that ends in
        !> another volume is shortened to land on the interface, where the
        !> owner switches. On input h is the trial step, on output the step
        !> taken.
        type(spectre_fo_state_t), intent(inout) :: state
        real(dp), intent(inout) :: h
        real(dp), intent(in) :: ro0_bar
        integer, intent(out) :: status
        real(dp) :: x1(3), v1(3), u1(3), err, a, b, fa, fb, c, fc, iface_rho
        integer :: attempt, owner1, it, side, grazes

        grazes = 0
        do attempt = 1, 60
            call dp_trial(state, h, ro0_bar, x1, v1, u1, owner1, err, status)
            if (status == SPECTRE_FO_LOSS) then
                h = 0.0_dp
                return
            end if
            if (status /= SPECTRE_FO_OK .or. err > 1.0_dp) then
                h = h*max(0.2_dp, min(0.5_dp, 0.9_dp*max(err, 1.0_dp)**(-0.2_dp)))
                cycle
            end if
            state%h_rk = h*min(5.0_dp, 0.9_dp*max(err, 1.0e-10_dp)**(-0.2_dp))
            if (state%layer_rho <= 0.0_dp .and. owner1 /= state%owner .and. &
                grazes < 2) then
                iface_rho = real(min(owner1, state%owner), dp)
                a = 0.0_dp
                fa = state%u(1) - iface_rho
                b = h
                fb = u1(1) - iface_rho
                side = 0
                do it = 1, 100
                    c = (a*fb - b*fa)/(fb - fa)
                    call dp_trial(state, c, ro0_bar, x1, v1, u1, owner1, err, &
                        status)
                    if (status /= SPECTRE_FO_OK) return
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
                if (c <= 1.0e-10_dp*h) then
                    ! The orbit grazes the interface at the step start: it
                    ! belongs to the other volume. Switch and redo the step; a
                    ! second graze in a row is a tangency, where the step is
                    ! taken with the volume of its end point.
                    state%owner = merge(state%owner + 1, state%owner - 1, &
                        iface_rho >= real(state%owner, dp))
                    state%u(1) = iface_rho
                    grazes = grazes + 1
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
            return
        end do
        status = SPECTRE_FO_STEP_LIMIT
    end subroutine rk45_substep
    subroutine dp_trial(state, h, ro0_bar, x1, v1, u1, owner1, err, status)
        !> One Dormand-Prince 5(4) step from the state; err is the scaled error
        !> (gyroradius for x, |v| for v) over fo_tolerance.
        type(spectre_fo_state_t), intent(in) :: state
        real(dp), intent(in) :: h, ro0_bar
        real(dp), intent(out) :: x1(3), v1(3), u1(3), err
        integer, intent(out) :: owner1, status
        real(dp) :: kx(3, 7), kv(3, 7), xs(3), vs(3), w(7, 7), bvec(3), bmod
        real(dp) :: rho_l, ex(3), ev(3)
        integer :: i, j

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
            call stage_field(state, xs, u1, bvec, bmod, status)
            if (status /= SPECTRE_FO_OK) return
            if (i == 1) rho_l = ro0_bar*norm2(vs)/bmod
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
        if (status /= SPECTRE_FO_OK) then
            if (state%owner == spectre_field%data%Mvol .and. &
                state%u(1) >= real(state%owner, dp) - 0.1_dp) then
                status = SPECTRE_FO_LOSS
            end if
            return
        end if
        if (owner1 == spectre_field%data%Mvol .and. &
            u1(1) >= real(owner1, dp) - OWNER_EPS) status = SPECTRE_FO_LOSS
    end subroutine dp_trial
    subroutine stage_field(state, x, u, bvec, bmod, status)
        !> Field at a stage point: the owner's field on the sharp sheet, the
        !> blended field of the geometric owner with a layer.
        type(spectre_fo_state_t), intent(in) :: state
        real(dp), intent(in) :: x(3)
        real(dp), intent(inout) :: u(3)
        real(dp), intent(out) :: bvec(3), bmod
        integer, intent(out) :: status
        real(dp) :: acov(3), e_zeta(3), xs(3)
        integer :: owner_geo, owner_eval

        owner_geo = state%owner
        call locate_cart(x, u, owner_geo, status)
        if (status /= SPECTRE_FO_OK) return
        owner_eval = state%owner
        if (state%layer_rho > 0.0_dp) owner_eval = owner_geo
        call logical_sample(u, owner_eval, xs, bvec, bmod, acov, e_zeta, status)
        if (status /= SPECTRE_FO_OK) return
        call layer_sample(state%layer_rho, u, owner_eval, bvec, bmod, status)
    end subroutine stage_field
    subroutine boris_substep(state, h, ro0_bar, status)
        type(spectre_fo_state_t), intent(inout) :: state
        real(dp), intent(inout) :: h
        real(dp), intent(in) :: ro0_bar
        integer, intent(out) :: status
        real(dp) :: b(3), bmod, acov(3), e_zeta(3), t(3), s(3), vp(3), vnew(3)
        real(dp) :: xmid(3), xnew(3), h_bound, u_mid(3), u_new(3)
        integer :: owner_mid, owner_new, attempt
        do attempt = 1, 40
            xmid = state%x + 0.5_dp*h*state%v
            u_mid = state%u
            owner_mid = state%owner
            call cart_sample(xmid, u_mid, owner_mid, b, bmod, acov, &
                e_zeta, status)
            if (status == SPECTRE_FO_OK) &
                call layer_sample(state%layer_rho, u_mid, owner_mid, b, bmod, &
                    status)
            if (status /= SPECTRE_FO_OK) then
                if (state%owner == spectre_field%data%Mvol .and. &
                    state%u(1) >= real(state%owner, dp) - OWNER_EPS) then
                    h = 0.0_dp
                    status = SPECTRE_FO_LOSS
                    return
                end if
                h = 0.5_dp*h
                cycle
            end if
            h_bound = 2.0_dp*rotation_parameter*ro0_bar/bmod
            if (h > h_bound*(1.0_dp + 4.0_dp*epsilon(1.0_dp))) then
                h = h_bound
                cycle
            end if
            t = 0.5_dp*h*bmod*b/ro0_bar
            s = 2.0_dp*t/(1.0_dp + dot_product(t, t))
            vp = state%v + cross(state%v, t)
            vnew = state%v + cross(vp, s)
            xnew = xmid + 0.5_dp*h*vnew
            u_new = u_mid
            owner_new = owner_mid
            call locate_cart(xnew, u_new, owner_new, status)
            if (status /= SPECTRE_FO_OK) then
                if (state%owner == spectre_field%data%Mvol .and. &
                    state%u(1) >= real(state%owner, dp) - OWNER_EPS) then
                    h = 0.0_dp
                    status = SPECTRE_FO_LOSS
                    return
                end if
                h = 0.5_dp*h
                cycle
            end if
            if (owner_new == spectre_field%data%Mvol .and. &
                u_new(1) >= real(owner_new, dp) - OWNER_EPS) then
                u_new(1) = real(owner_new, dp)
            end if
            state%x = xnew
            state%v = vnew
            state%u = u_new
            state%owner = owner_new
            status = SPECTRE_FO_OK
            return
        end do
        status = SPECTRE_FO_STEP_LIMIT
    end subroutine boris_substep
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
            ! take the closest gyrophase, an O(gyroradius^2) p_zeta error.
            if (best_residual == huge(1.0_dp)) then
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
        residual = dot_product(v, e_zeta) + acov(3)/ro0_bar - target
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
                status = SPECTRE_FO_OK
                return
            end if
        class default
            status = SPECTRE_FO_INVALID_STATE
            return
        end select
        status = SPECTRE_FO_LOCATE_FAIL
    end subroutine locate_cart

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
