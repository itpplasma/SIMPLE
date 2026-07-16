module spectre_sheet_gc
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use libneo_coordinates, only: spectre_coordinate_system_t
    use magfie_sub, only: spectre_field, TESLA_TO_GAUSS, M_TO_CM, M3_TO_CM3
    use diag_counters, only: count_event, EVT_SPECTRE_INVALID_STATE
    implicit none
    private
    public :: sheet_gc_sample_t, sheet_gc_state_t, evaluate_sheet_profile
    public :: sheet_gc_rhs, sheet_gc_initialize, sheet_gc_advance, sheet_gc_to_y
    public :: SHEET_GC_OK, SHEET_GC_DEGENERATE, SHEET_GC_UNINITIALIZED
    public :: SHEET_GC_STEP_LIMIT, SHEET_GC_ENERGY_BARRIER
    public :: SHEET_GC_VPAR_ZERO
    integer, parameter :: SHEET_GC_OK = 0
    integer, parameter :: SHEET_GC_DEGENERATE = 1
    integer, parameter :: SHEET_GC_UNINITIALIZED = 2
    integer, parameter :: SHEET_GC_STEP_LIMIT = 3
    integer, parameter :: SHEET_GC_ENERGY_BARRIER = 4
    integer, parameter :: SHEET_GC_VPAR_ZERO = 5
    real(dp), parameter :: SIDE_EPS = 1.0e-12_dp
    real(dp), parameter :: DEGENERACY_FACTOR = 100.0_dp
    real(dp), parameter :: PROFILE_EDGE = 4.0_dp
    real(dp), parameter :: MAX_ETA_STEP = 0.05_dp
    real(dp), parameter :: PURE_ETA_STEP = 0.02_dp
    integer, parameter :: MAX_SUBSTEPS = 16384
    type :: sheet_gc_sample_t
        real(dp) :: bmod = 0.0_dp
        real(dp) :: sqrtg = 0.0_dp
        real(dp) :: hcov(3) = 0.0_dp
        real(dp) :: h_eta(3) = 0.0_dp
        real(dp) :: db(3) = 0.0_dp
        real(dp) :: hder(3, 3) = 0.0_dp
        real(dp) :: hctr(3) = 0.0_dp
        real(dp) :: hcurl(3) = 0.0_dp
        real(dp) :: rotation = 0.0_dp
    end type sheet_gc_sample_t
    type :: sheet_gc_state_t
        logical :: active = .false.
        integer :: iface = 0
        integer :: owner = 0
        real(dp) :: eta = 0.0_dp
        real(dp) :: theta = 0.0_dp
        real(dp) :: zeta = 0.0_dp
        real(dp) :: vpar = 0.0_dp
        real(dp) :: p = 0.0_dp
        real(dp) :: mu = 0.0_dp
        real(dp) :: energy = 0.0_dp
    end type sheet_gc_state_t
contains
    elemental logical function ieee_is_finite(value)
        real(dp), intent(in) :: value

        integer(int64), parameter :: exponent_mask = &
            int(z'7FF0000000000000', int64)
        integer(int64) :: bits

        bits = transfer(value, bits)
        ieee_is_finite = iand(bits, exponent_mask) /= exponent_mask
    end function ieee_is_finite

    subroutine evaluate_sheet_profile(iface, eta, theta, zeta, sample, ierr)
        integer, intent(in) :: iface
        real(dp), intent(in) :: eta, theta, zeta
        type(sheet_gc_sample_t), intent(out) :: sample
        integer, intent(out) :: ierr
        real(dp) :: hm(3), hp(3), bm, bp, dhm(3, 3), dhp(3, 3)
        real(dp) :: dbm(3), dbp(3), ginv(3, 3), dginv(3, 3, 3), sqrtg
        real(dp) :: blend, blend_eta
        sample = sheet_gc_sample_t()
        if (.not. all(ieee_is_finite([eta, theta, zeta]))) then
            ierr = SHEET_GC_DEGENERATE
            call count_event(EVT_SPECTRE_INVALID_STATE)
            return
        end if
        if (.not. allocated(spectre_field)) then
            ierr = SHEET_GC_UNINITIALIZED
            return
        end if
        call one_sided_fields(iface, theta, zeta, hm, hp, bm, bp, dhm, dhp, &
            dbm, dbp)
        call interface_metric(iface, theta, zeta, ginv, dginv, sqrtg)
        if (.not. all(ieee_is_finite(hm)) .or. &
            .not. all(ieee_is_finite(hp)) .or. &
            .not. ieee_is_finite(bm) .or. .not. ieee_is_finite(bp) .or. &
            .not. all(ieee_is_finite(dhm)) .or. &
            .not. all(ieee_is_finite(dhp)) .or. &
            .not. all(ieee_is_finite(dbm)) .or. &
            .not. all(ieee_is_finite(dbp)) .or. &
            .not. all(ieee_is_finite(ginv)) .or. &
            .not. all(ieee_is_finite(dginv)) .or. &
            .not. ieee_is_finite(sqrtg) .or. abs(sqrtg) <= tiny(1.0_dp)) then
            ierr = SHEET_GC_DEGENERATE
            call count_event(EVT_SPECTRE_INVALID_STATE)
            return
        end if
        blend = 0.5_dp*(1.0_dp + tanh(eta))
        blend_eta = 2.0_dp*blend*(1.0_dp - blend)
        call normalized_profile(blend, blend_eta, hm, hp, dhm, dhp, ginv, &
            dginv, sample%hcov, sample%h_eta, sample%hder)
        sample%bmod = ((1.0_dp - blend)*bm + blend*bp)*TESLA_TO_GAUSS
        sample%db = ((1.0_dp - blend)*dbm + blend*dbp)*TESLA_TO_GAUSS
        sample%db(1) = blend_eta*(bp - bm)*TESLA_TO_GAUSS
        sample%sqrtg = sqrtg*M3_TO_CM3
        sample%hcov = sample%hcov*M_TO_CM
        sample%h_eta = sample%h_eta*M_TO_CM
        sample%hder = sample%hder*M_TO_CM
        sample%hder(:, 1) = 0.0_dp
        sample%hctr = matmul(ginv/M_TO_CM**2, sample%hcov)
        sample%hcurl(1) = (sample%hder(3, 2) - sample%hder(2, 3))/sample%sqrtg
        sample%hcurl(2) = (sample%hder(1, 3) - sample%hder(3, 1))/sample%sqrtg
        sample%hcurl(3) = (sample%hder(2, 1) - sample%hder(1, 2))/sample%sqrtg
        sample%rotation = sample%hcov(3)*sample%h_eta(2) - &
            sample%hcov(2)*sample%h_eta(3)
        if (.not. ieee_is_finite(sample%bmod) .or. sample%bmod <= 0.0_dp .or. &
            .not. ieee_is_finite(sample%sqrtg) .or. &
            abs(sample%sqrtg) <= tiny(1.0_dp) .or. &
            .not. all(ieee_is_finite(sample%hcov)) .or. &
            .not. all(ieee_is_finite(sample%h_eta)) .or. &
            .not. all(ieee_is_finite(sample%db)) .or. &
            .not. all(ieee_is_finite(sample%hder)) .or. &
            .not. all(ieee_is_finite(sample%hctr)) .or. &
            .not. all(ieee_is_finite(sample%hcurl)) .or. &
            .not. ieee_is_finite(sample%rotation)) then
            sample = sheet_gc_sample_t()
            ierr = SHEET_GC_DEGENERATE
            call count_event(EVT_SPECTRE_INVALID_STATE)
            return
        end if
        ierr = SHEET_GC_OK
    end subroutine evaluate_sheet_profile

    subroutine sheet_gc_rhs(iface, eta, theta, zeta, vpar, mu, rho0, rhs, &
            energy_rate, ierr)
        integer, intent(in) :: iface
        real(dp), intent(in) :: eta, theta, zeta, vpar, mu, rho0
        real(dp), intent(out) :: rhs(4), energy_rate
        integer, intent(out) :: ierr
        type(sheet_gc_sample_t) :: s
        real(dp) :: scale, bstar_s, normal_n
        rhs = 0.0_dp
        energy_rate = 0.0_dp
        if (.not. all(ieee_is_finite([eta, theta, zeta, vpar, mu, rho0]))) then
            ierr = SHEET_GC_DEGENERATE
            call count_event(EVT_SPECTRE_INVALID_STATE)
            return
        end if
        call evaluate_sheet_profile(iface, eta, theta, zeta, s, ierr)
        if (ierr /= SHEET_GC_OK) return
        scale = abs(s%hcov(3)*s%h_eta(2)) + abs(s%hcov(2)*s%h_eta(3))
        if (abs(s%rotation) <= DEGENERACY_FACTOR*epsilon(1.0_dp)*scale &
            .or. rho0 == 0.0_dp) then
            ierr = SHEET_GC_DEGENERATE
            return
        end if
        if (abs(vpar) <= DEGENERACY_FACTOR*sqrt(epsilon(1.0_dp))) then
            ierr = SHEET_GC_VPAR_ZERO
            return
        end if

        bstar_s = s%bmod*s%hctr(1) + rho0*vpar*s%hcurl(1)
        normal_n = vpar*bstar_s + rho0*mu/s%sqrtg* &
            (s%hcov(2)*s%db(3) - s%hcov(3)*s%db(2))
        rhs(1) = s%sqrtg*normal_n/(rho0*vpar*s%rotation)
        rhs(2) = (-vpar**2*s%h_eta(3) + mu*s%hcov(3)*s%db(1))/ &
            (vpar*s%rotation)
        rhs(3) = (vpar**2*s%h_eta(2) - mu*s%hcov(2)*s%db(1))/ &
            (vpar*s%rotation)
        rhs(4) = -mu*(s%sqrtg*bstar_s*s%db(1)/(rho0*vpar*s%rotation) + &
            (-s%h_eta(3)*s%db(2) + s%h_eta(2)*s%db(3))/s%rotation)
        energy_rate = vpar*rhs(4) + mu* &
            (s%db(1)*rhs(1) + s%db(2)*rhs(2) + s%db(3)*rhs(3))
        if (.not. all(ieee_is_finite(rhs)) .or. &
            .not. ieee_is_finite(energy_rate)) then
            rhs = 0.0_dp
            energy_rate = 0.0_dp
            ierr = SHEET_GC_DEGENERATE
            call count_event(EVT_SPECTRE_INVALID_STATE)
            return
        end if
        ierr = SHEET_GC_OK
    end subroutine sheet_gc_rhs

    subroutine sheet_gc_initialize(iface, direction, y, mu, b_home, rho0, state, &
            ierr)
        integer, intent(in) :: iface, direction
        real(dp), intent(in) :: y(5), mu, b_home, rho0
        type(sheet_gc_state_t), intent(out) :: state
        integer, intent(out) :: ierr

        type(sheet_gc_sample_t) :: sample
        real(dp) :: radicand, incoming, rhs(4), energy_rate

        state = sheet_gc_state_t()
        if (.not. all(ieee_is_finite(y)) .or. .not. ieee_is_finite(mu) .or. &
            .not. ieee_is_finite(b_home) .or. .not. ieee_is_finite(rho0) .or. &
            y(4) <= 0.0_dp .or. mu < 0.0_dp .or. b_home <= 0.0_dp) then
            ierr = SHEET_GC_DEGENERATE
            call count_event(EVT_SPECTRE_INVALID_STATE)
            return
        end if
        state%iface = iface
        state%owner = merge(iface, iface + 1, direction > 0)
        state%eta = merge(-PROFILE_EDGE, PROFILE_EDGE, direction > 0)
        state%theta = y(2)
        state%zeta = y(3)
        incoming = y(4)*y(5)
        state%mu = mu
        state%energy = 0.5_dp*incoming**2 + mu*b_home
        state%p = sqrt(2.0_dp*state%energy)
        call evaluate_sheet_profile(iface, state%eta, state%theta, state%zeta, &
            sample, ierr)
        if (ierr /= SHEET_GC_OK) return
        radicand = 2.0_dp*(state%energy - mu*sample%bmod)
        if (.not. ieee_is_finite(radicand) .or. radicand <= 0.0_dp) then
            ierr = SHEET_GC_ENERGY_BARRIER
            return
        end if
        if (rho0 == 0.0_dp) then
            ierr = SHEET_GC_DEGENERATE
            return
        end if
        state%vpar = sign(sqrt(radicand), incoming)
        call sheet_gc_rhs(iface, state%eta, state%theta, state%zeta, &
            state%vpar, state%mu, rho0, rhs, energy_rate, ierr)
        if (ierr == SHEET_GC_DEGENERATE) then
            state%active = .true.
            ierr = SHEET_GC_OK
            return
        end if
        if (ierr /= SHEET_GC_OK) return
        if (.not. all(abs(rhs) <= huge(1.0_dp)) .or. &
            .not. (abs(energy_rate) <= huge(1.0_dp))) then
            ierr = SHEET_GC_DEGENERATE
            return
        end if
        state%active = .true.
    end subroutine sheet_gc_initialize

    subroutine sheet_gc_advance(state, dt, rho0, exited, dt_used, ierr)
        type(sheet_gc_state_t), intent(inout) :: state
        real(dp), intent(in) :: dt, rho0
        logical, intent(out) :: exited
        real(dp), intent(out) :: dt_used
        integer, intent(out) :: ierr

        real(dp) :: remaining, h, h_exit, rhs(4), energy_rate
        type(sheet_gc_state_t) :: state0, valid_state
        integer :: nstep

        exited = .false.
        dt_used = 0.0_dp
        ierr = SHEET_GC_OK
        if (.not. sheet_state_finite(state) .or. .not. ieee_is_finite(dt) .or. &
            .not. ieee_is_finite(rho0)) then
            ierr = SHEET_GC_DEGENERATE
            call count_event(EVT_SPECTRE_INVALID_STATE)
            return
        end if
        if (dt <= 0.0_dp) return
        remaining = dt
        valid_state = state
        do nstep = 1, MAX_SUBSTEPS
            call sheet_gc_rhs(state%iface, state%eta, state%theta, state%zeta, &
                state%vpar, state%mu, rho0, rhs, energy_rate, ierr)
            if (ierr == SHEET_GC_DEGENERATE) then
                state = valid_state
                call pure_magnitude_exit(state, rho0, ierr)
                if (ierr == SHEET_GC_OK) exited = .true.
                return
            end if
            if (ierr /= SHEET_GC_OK) return
            h = min(remaining, MAX_ETA_STEP/max(abs(rhs(1)), tiny(1.0_dp)))
            if (.not. ieee_is_finite(h) .or. h <= 0.0_dp) then
                state = valid_state
                ierr = SHEET_GC_DEGENERATE
                call count_event(EVT_SPECTRE_INVALID_STATE)
                return
            end if
            state0 = state
            call rk4_sheet_step(state, h, rho0, ierr)
            if (ierr == SHEET_GC_DEGENERATE) then
                state = state0
                call pure_magnitude_exit(state, rho0, ierr)
                if (ierr == SHEET_GC_OK) exited = .true.
                return
            end if
            if (ierr /= SHEET_GC_OK) return
            if (.not. sheet_state_finite(state)) then
                state = state0
                ierr = SHEET_GC_DEGENERATE
                call count_event(EVT_SPECTRE_INVALID_STATE)
                return
            end if
            valid_state = state
            remaining = remaining - h
            if (abs(state%eta) >= PROFILE_EDGE) then
                call locate_sheet_exit(state0, h, rho0, &
                    sign(PROFILE_EDGE, state%eta), state, h_exit, ierr)
                if (ierr == SHEET_GC_DEGENERATE) then
                    state = state0
                    call pure_magnitude_exit(state, rho0, ierr)
                    dt_used = dt - remaining - h
                    if (ierr == SHEET_GC_OK) exited = .true.
                    return
                end if
                if (ierr /= SHEET_GC_OK) return
                dt_used = dt - remaining - h + h_exit
                state%owner = merge(state%iface + 1, state%iface, state%eta > 0.0_dp)
                state%active = .false.
                exited = .true.
                return
            end if
            dt_used = dt - remaining
            if (remaining <= epsilon(1.0_dp)*dt) return
        end do
        ierr = SHEET_GC_STEP_LIMIT
    end subroutine sheet_gc_advance

    subroutine sheet_gc_to_y(state, y)
        type(sheet_gc_state_t), intent(in) :: state
        real(dp), intent(out) :: y(5)

        y = [real(state%iface, dp), state%theta, state%zeta, state%p, &
            state%vpar/state%p]
    end subroutine sheet_gc_to_y

    pure logical function sheet_state_finite(state)
        type(sheet_gc_state_t), intent(in) :: state

        sheet_state_finite = all(ieee_is_finite([state%eta, state%theta, &
            state%zeta, state%vpar, state%p, state%mu, state%energy])) .and. &
            state%p > 0.0_dp .and. state%mu >= 0.0_dp .and. &
            state%energy > 0.0_dp
    end function sheet_state_finite

    subroutine rk4_sheet_step(state, h, rho0, ierr)
        type(sheet_gc_state_t), intent(inout) :: state
        real(dp), intent(in) :: h, rho0
        integer, intent(out) :: ierr

        real(dp) :: z(4), k1(4), k2(4), k3(4), k4(4), rate

        z = [state%eta, state%theta, state%zeta, state%vpar]
        call sheet_gc_rhs(state%iface, z(1), z(2), z(3), z(4), state%mu, &
            rho0, k1, rate, ierr)
        if (ierr /= SHEET_GC_OK) return
        call sheet_gc_rhs(state%iface, z(1) + 0.5_dp*h*k1(1), &
            z(2) + 0.5_dp*h*k1(2), z(3) + 0.5_dp*h*k1(3), &
            z(4) + 0.5_dp*h*k1(4), state%mu, rho0, k2, rate, ierr)
        if (ierr /= SHEET_GC_OK) return
        call sheet_gc_rhs(state%iface, z(1) + 0.5_dp*h*k2(1), &
            z(2) + 0.5_dp*h*k2(2), z(3) + 0.5_dp*h*k2(3), &
            z(4) + 0.5_dp*h*k2(4), state%mu, rho0, k3, rate, ierr)
        if (ierr /= SHEET_GC_OK) return
        call sheet_gc_rhs(state%iface, z(1) + h*k3(1), z(2) + h*k3(2), &
            z(3) + h*k3(3), z(4) + h*k3(4), state%mu, rho0, k4, rate, ierr)
        if (ierr /= SHEET_GC_OK) return
        z = z + h*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
        if (.not. all(ieee_is_finite(z))) then
            ierr = SHEET_GC_DEGENERATE
            call count_event(EVT_SPECTRE_INVALID_STATE)
            return
        end if
        call project_sheet_energy(state, z, ierr)
        if (ierr /= SHEET_GC_OK) return
        state%eta = z(1)
        state%theta = z(2)
        state%zeta = z(3)
        state%vpar = z(4)
    end subroutine rk4_sheet_step

    subroutine project_sheet_energy(state, z, ierr)
        type(sheet_gc_state_t), intent(in) :: state
        real(dp), intent(inout) :: z(4)
        integer, intent(out) :: ierr

        type(sheet_gc_sample_t) :: sample
        real(dp) :: radicand

        call evaluate_sheet_profile(state%iface, z(1), z(2), z(3), sample, ierr)
        if (ierr /= SHEET_GC_OK) return
        if (.not. all(ieee_is_finite(z)) .or. .not. sheet_state_finite(state)) then
            ierr = SHEET_GC_DEGENERATE
            call count_event(EVT_SPECTRE_INVALID_STATE)
            return
        end if
        radicand = 2.0_dp*(state%energy - state%mu*sample%bmod)
        if (.not. ieee_is_finite(radicand) .or. radicand <= 0.0_dp) then
            ierr = SHEET_GC_ENERGY_BARRIER
            return
        end if
        z(4) = sign(sqrt(radicand), z(4))
        ierr = SHEET_GC_OK
    end subroutine project_sheet_energy

    subroutine locate_sheet_exit(initial, h_full, rho0, target, state, h_exit, &
            ierr)
        type(sheet_gc_state_t), intent(in) :: initial
        real(dp), intent(in) :: h_full, rho0, target
        type(sheet_gc_state_t), intent(out) :: state
        real(dp), intent(out) :: h_exit
        integer, intent(out) :: ierr

        real(dp) :: lo, hi, mid, z(4)
        integer :: iter

        lo = 0.0_dp
        hi = h_full
        do iter = 1, 48
            mid = 0.5_dp*(lo + hi)
            state = initial
            call rk4_sheet_step(state, mid, rho0, ierr)
            if (ierr /= SHEET_GC_OK) return
            if ((state%eta - target)*(initial%eta - target) > 0.0_dp) then
                lo = mid
            else
                hi = mid
            end if
        end do
        h_exit = hi
        state = initial
        call rk4_sheet_step(state, h_exit, rho0, ierr)
        if (ierr /= SHEET_GC_OK) return
        state%eta = target
        z = [state%eta, state%theta, state%zeta, state%vpar]
        call project_sheet_energy(state, z, ierr)
        if (ierr /= SHEET_GC_OK) return
        state%vpar = z(4)
    end subroutine locate_sheet_exit

    subroutine pure_magnitude_exit(state, rho0, ierr)
        !> Zero-width fast map with t = w*tau: it changes the sheet state while
        !> consuming zero leading-order physical time as w tends to zero.
        type(sheet_gc_state_t), intent(inout) :: state
        real(dp), intent(in) :: rho0
        integer, intent(out) :: ierr
        type(sheet_gc_state_t) :: state0
        real(dp) :: rhs(4), rate, step, edge
        integer :: nstep
        if (.not. sheet_state_finite(state) .or. .not. ieee_is_finite(rho0)) then
            ierr = SHEET_GC_DEGENERATE
            call count_event(EVT_SPECTRE_INVALID_STATE)
            return
        end if
        do nstep = 1, MAX_SUBSTEPS
            call pure_magnitude_rhs(state%iface, state%eta, state%theta, &
                state%zeta, state%vpar, state%mu, rho0, rhs, ierr)
            if (ierr /= SHEET_GC_OK) return
            rate = max(maxval(abs(rhs(1:3))), &
                abs(rhs(4))/max(abs(state%vpar), 1.0e-3_dp))
            step = PURE_ETA_STEP/max(rate, tiny(1.0_dp))
            if (.not. ieee_is_finite(step) .or. step <= 0.0_dp) then
                ierr = SHEET_GC_DEGENERATE
                call count_event(EVT_SPECTRE_INVALID_STATE)
                return
            end if
            state0 = state
            call rk4_pure_magnitude_step(state, step, rho0, ierr)
            if (ierr /= SHEET_GC_OK) then
                state = state0
                return
            end if
            if (.not. sheet_state_finite(state)) then
                state = state0
                ierr = SHEET_GC_DEGENERATE
                call count_event(EVT_SPECTRE_INVALID_STATE)
                return
            end if
            if (abs(state%eta) >= PROFILE_EDGE) then
                edge = sign(PROFILE_EDGE, state%eta)
                call locate_pure_magnitude_exit(state0, step, rho0, edge, state, &
                    ierr)
                if (ierr /= SHEET_GC_OK) return
                state%owner = merge(state%iface + 1, state%iface, edge > 0.0_dp)
                state%active = .false.
                return
            end if
        end do
        ierr = SHEET_GC_STEP_LIMIT
    end subroutine pure_magnitude_exit

    subroutine pure_magnitude_rhs(iface, eta, theta, zeta, vpar, mu, rho0, rhs, &
            ierr)
        integer, intent(in) :: iface
        real(dp), intent(in) :: eta, theta, zeta, vpar, mu, rho0
        real(dp), intent(out) :: rhs(4)
        integer, intent(out) :: ierr
        type(sheet_gc_sample_t) :: s
        real(dp) :: bstar_s, bstar_parallel, normal_n
        rhs = 0.0_dp
        if (.not. all(ieee_is_finite([eta, theta, zeta, vpar, mu, rho0]))) then
            ierr = SHEET_GC_DEGENERATE
            call count_event(EVT_SPECTRE_INVALID_STATE)
            return
        end if
        call evaluate_sheet_profile(iface, eta, theta, zeta, s, ierr)
        if (ierr /= SHEET_GC_OK) return
        bstar_s = s%bmod*s%hctr(1)
        bstar_parallel = s%bmod
        normal_n = vpar*bstar_s + rho0*mu/s%sqrtg* &
            (s%hcov(2)*s%db(3) - s%hcov(3)*s%db(2))
        if (abs(bstar_parallel) <= DEGENERACY_FACTOR*epsilon(1.0_dp)*s%bmod) then
            ierr = SHEET_GC_DEGENERATE
            return
        end if
        rhs(1) = normal_n/bstar_parallel
        rhs(2) = rho0*mu*s%hcov(3)*s%db(1)/(s%sqrtg*bstar_parallel)
        rhs(3) = -rho0*mu*s%hcov(2)*s%db(1)/(s%sqrtg*bstar_parallel)
        rhs(4) = -mu*bstar_s*s%db(1)/bstar_parallel
        if (.not. all(ieee_is_finite(rhs))) then
            rhs = 0.0_dp
            ierr = SHEET_GC_DEGENERATE
            call count_event(EVT_SPECTRE_INVALID_STATE)
            return
        end if
        ierr = SHEET_GC_OK
    end subroutine pure_magnitude_rhs

    subroutine rk4_pure_magnitude_step(state, h, rho0, ierr)
        type(sheet_gc_state_t), intent(inout) :: state
        real(dp), intent(in) :: h, rho0
        integer, intent(out) :: ierr
        real(dp) :: z(4), k1(4), k2(4), k3(4), k4(4)
        z = [state%eta, state%theta, state%zeta, state%vpar]
        call pure_magnitude_rhs(state%iface, z(1), z(2), z(3), z(4), state%mu, &
            rho0, k1, ierr)
        if (ierr /= SHEET_GC_OK) return
        call pure_magnitude_rhs(state%iface, z(1) + 0.5_dp*h*k1(1), &
            z(2) + 0.5_dp*h*k1(2), z(3) + 0.5_dp*h*k1(3), &
            z(4) + 0.5_dp*h*k1(4), state%mu, rho0, k2, ierr)
        if (ierr /= SHEET_GC_OK) return
        call pure_magnitude_rhs(state%iface, z(1) + 0.5_dp*h*k2(1), &
            z(2) + 0.5_dp*h*k2(2), z(3) + 0.5_dp*h*k2(3), &
            z(4) + 0.5_dp*h*k2(4), state%mu, rho0, k3, ierr)
        if (ierr /= SHEET_GC_OK) return
        call pure_magnitude_rhs(state%iface, z(1) + h*k3(1), z(2) + h*k3(2), &
            z(3) + h*k3(3), z(4) + h*k3(4), state%mu, rho0, k4, ierr)
        if (ierr /= SHEET_GC_OK) return
        z = z + h*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
        if (.not. all(ieee_is_finite(z))) then
            ierr = SHEET_GC_DEGENERATE
            call count_event(EVT_SPECTRE_INVALID_STATE)
            return
        end if
        call project_sheet_energy(state, z, ierr)
        if (ierr /= SHEET_GC_OK) return
        state%eta = z(1)
        state%theta = z(2)
        state%zeta = z(3)
        state%vpar = z(4)
    end subroutine rk4_pure_magnitude_step

    subroutine locate_pure_magnitude_exit(initial, h_full, rho0, edge, state, &
            ierr)
        type(sheet_gc_state_t), intent(in) :: initial
        real(dp), intent(in) :: h_full, rho0, edge
        type(sheet_gc_state_t), intent(out) :: state
        integer, intent(out) :: ierr
        real(dp) :: lo, hi, mid, z(4)
        integer :: iter
        lo = 0.0_dp
        hi = h_full
        do iter = 1, 48
            mid = 0.5_dp*(lo + hi)
            state = initial
            call rk4_pure_magnitude_step(state, mid, rho0, ierr)
            if (ierr /= SHEET_GC_OK) return
            if ((state%eta - edge)*(initial%eta - edge) > 0.0_dp) then
                lo = mid
            else
                hi = mid
            end if
        end do
        state = initial
        call rk4_pure_magnitude_step(state, hi, rho0, ierr)
        if (ierr /= SHEET_GC_OK) return
        state%eta = edge
        z = [state%eta, state%theta, state%zeta, state%vpar]
        call project_sheet_energy(state, z, ierr)
        state%vpar = z(4)
    end subroutine locate_pure_magnitude_exit
    subroutine one_sided_fields(iface, theta, zeta, hm, hp, bm, bp, dhm, dhp, &
            dbm, dbp)
        integer, intent(in) :: iface
        real(dp), intent(in) :: theta, zeta
        real(dp), intent(out) :: hm(3), hp(3), bm, bp
        real(dp), intent(out) :: dhm(3, 3), dhp(3, 3), dbm(3), dbp(3)
        real(dp) :: sqgb(3), sqrtg

        call spectre_field%evaluate_der( &
            [real(iface, dp) - SIDE_EPS, theta, zeta], hm, bm, sqgb, sqrtg, &
            dhm, dbm)
        call spectre_field%evaluate_der( &
            [real(iface, dp) + SIDE_EPS, theta, zeta], hp, bp, sqgb, sqrtg, &
            dhp, dbp)
    end subroutine one_sided_fields

    subroutine interface_metric(iface, theta, zeta, ginv, dginv, sqrtg)
        integer, intent(in) :: iface
        real(dp), intent(in) :: theta, zeta
        real(dp), intent(out) :: ginv(3, 3), dginv(3, 3, 3), sqrtg

        real(dp) :: x(3), g(3, 3), dg(3, 3, 3), dsqrtg(3)
        integer :: j

        x = [real(iface, dp), theta, zeta]
        select type (coords => spectre_field%coords)
            type is (spectre_coordinate_system_t)
            call coords%metric_tensor_der(x, g, ginv, sqrtg, dg, dsqrtg)
        class default
            error stop 'spectre_sheet_gc: invalid coordinate system'
        end select
        do j = 1, 3
            dginv(:, :, j) = -matmul(ginv, matmul(dg(:, :, j), ginv))
        end do
    end subroutine interface_metric

    subroutine normalized_profile(p, peta, hm, hp, dhm, dhp, ginv, dginv, &
            h, heta, hder)
        real(dp), intent(in) :: p, peta, hm(3), hp(3)
        real(dp), intent(in) :: dhm(3, 3), dhp(3, 3), ginv(3, 3)
        real(dp), intent(in) :: dginv(3, 3, 3)
        real(dp), intent(out) :: h(3), heta(3), hder(3, 3)

        real(dp) :: u(3), ueta(3), uder(3, 3), gu(3), norm, dnorm
        integer :: j

        u = (1.0_dp - p)*hm + p*hp
        ueta = peta*(hp - hm)
        uder = (1.0_dp - p)*dhm + p*dhp
        gu = matmul(ginv, u)
        norm = sqrt(dot_product(u, gu))
        h = u/norm
        heta = ueta/norm - u*dot_product(ueta, gu)/norm**3
        do j = 1, 3
            dnorm = (2.0_dp*dot_product(uder(:, j), gu) + &
                dot_product(u, matmul(dginv(:, :, j), u)))/(2.0_dp*norm)
            hder(:, j) = uder(:, j)/norm - u*dnorm/norm**2
        end do
    end subroutine normalized_profile
end module spectre_sheet_gc
