module spectre_fo_hybrid
    use, intrinsic :: iso_fortran_env, only: dp => real64
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
    real(dp), parameter :: MAX_ROTATION_PARAMETER = 0.05_dp
    integer, parameter :: MAX_SUBSTEPS = 100000
    integer, parameter :: NPHASE = 64
    type, public :: spectre_fo_state_t
        real(dp) :: x(3) = 0.0_dp
        real(dp) :: v(3) = 0.0_dp
        real(dp) :: u(3) = 0.0_dp
        real(dp) :: last_y(5) = 0.0_dp
        real(dp) :: entry_rho = 0.0_dp
        real(dp) :: entry_pzeta = 0.0_dp
        integer :: owner = 0
        integer :: iface = 0
        logical :: has_y = .false.
        logical :: active = .false.
    end type spectre_fo_state_t
    public :: spectre_fo_enter, spectre_fo_to_gc
    public :: spectre_fo_advance_until_exit
    public :: spectre_fo_canonical_pzeta
contains
    subroutine spectre_fo_enter(state, y, owner, iface, ro0_bar, status)
        type(spectre_fo_state_t), intent(out) :: state
        real(dp), intent(in) :: y(5), ro0_bar
        integer, intent(in) :: owner, iface
        integer, intent(out) :: status
        real(dp) :: x_gc(3), b(3), e1(3), e2(3), eperp(3), e_zeta(3)
        real(dp) :: acov(3), bmod, vpar, vperp, target, alpha
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
        state%active = .true.
    end subroutine spectre_fo_enter
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
            h = min(remaining, hfac*2.0_dp*MAX_ROTATION_PARAMETER*ro0_bar/bmod)
            call boris_substep(state, h, ro0_bar, status)
            if (status /= SPECTRE_FO_OK) return
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
        if (state%iface == 0) then
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
        call cart_sample(state%x, state%u, state%owner, b, bmod, acov, &
            e_zeta, status)
        if (status /= SPECTRE_FO_OK) return
        pzeta = dot_product(state%v, e_zeta) + acov(3)/ro0_bar
    end subroutine spectre_fo_canonical_pzeta
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
            h_bound = 2.0_dp*MAX_ROTATION_PARAMETER*ro0_bar/bmod
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
        real(dp) :: a0, a1, f0, f1, amid, fmid, scale
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
        do i = 0, NPHASE
            a1 = 8.0_dp*atan(1.0_dp)*real(i, dp)/real(NPHASE, dp)
            call phase_residual(a1, x_gc, u_gc, owner, b, e1, e2, vpar, &
                vperp, target, ro0_bar, f1, ierr)
            if (ierr /= SPECTRE_FO_OK) then
                have_previous = .false.
                cycle
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
            status = SPECTRE_FO_GYROPHASE_FAIL
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
