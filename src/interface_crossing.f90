module interface_crossing
    !> Guiding-center interface crossing for SPECTRE stacked-rho charts. The
    !> single entry point apply_crossing selects the map through its level
    !> argument, so callers stay untouched between levels.
    !>
    !> Level 0 (#443) holds (theta, zeta, mu), keeps the sign of v_par, switches
    !> volume, and rescales v_par from exact energy conservation across the
    !> pressure jump [[B]].
    !>
    !> Level 1 (#440) is the thin-current-sheet limit of the guiding-center
    !> dynamics: an impulse Delta z = lambda_k * X along the Hamiltonian vector
    !> field X = {z, rho_g} of the interface function, with the single scalar
    !> lambda_k fixed by exact energy conservation H = v_par^2/2 + mu*B. It adds
    !> the tangential sheet-drift kick and the drift-order v_par term that Level 0
    !> omits. Both levels are energy-exact in the crossing and reflection branch.

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    use magfie_sub, only: magfie
    use parmot_mod, only: ro0

    implicit none
    private

    public :: apply_crossing, crossing_info_t, axis_offset
    public :: crossing_log_reset, crossing_log_record, crossing_log_write, &
        crossing_log_count, crossing_log_count_type
    public :: CROSSING_LEVEL0, CROSSING_LEVEL1, CROSS_CROSSING, CROSS_REFLECTION, &
        CROSS_LOSS, CROSS_STOP, CROSS_SHEET, CROSS_RECOVERY, CROSS_INVALID

    integer, parameter :: CROSSING_LEVEL0 = 0
    integer, parameter :: CROSSING_LEVEL1 = 1
    integer, parameter :: CROSS_CROSSING = 1
    integer, parameter :: CROSS_REFLECTION = 2
    integer, parameter :: CROSS_LOSS = 3
    !> Pathological fallback of the symplectic exact-landing pipeline (SIMPLE#441):
    !> a stop event terminates the orbit when the substep solve or the implicit
    !> step itself fails to converge; regular boundaries cross via apply_crossing.
    integer, parameter :: CROSS_STOP = 4
    integer, parameter :: CROSS_SHEET = 5
    !> Marker-local GC -> full-orbit fallback away from any physical interface.
    !> iface and both volume fields are zero so post-processing cannot confuse
    !> the recovery map with a discontinuity crossing.
    integer, parameter :: CROSS_RECOVERY = 6
    !> A candidate event whose physical state or field values are unusable.
    !> Callers must recover from the pre-step state; this is never a loss.
    integer, parameter :: CROSS_INVALID = 7

    !> Inner cutoff for the innermost volume: rho_g = 0 is the coordinate axis
    !> where sqrt(g) = 0, so the marker reflects trivially before reaching it
    !> rather than crossing a (non-existent) inner interface.
    real(dp), parameter :: axis_offset = 1.0d-3

    !> Both-side |B| is read at rho_g = k -+ iface_eps: int(k -+ 1e-12) selects
    !> the neighbouring volume k-1 resp. k, and the libneo polynomial basis of
    !> each volume is well-defined at that offset from the shared interface point.
    real(dp), parameter :: iface_eps = 1.0d-12
    !> Level-1 refraction solve: a damped (backtracking) Newton on the scalar
    !> lambda_k drives the energy residual |dH|/H below newton_rtol. Backtracking
    !> keeps the step from diverging where the residual is non-monotonic in the
    !> tangential kick (the target |B| varies poloidally). The generator is
    !> refreshed at the midpoint state by sym_iters symmetrization sweeps per
    !> residual, so the discrete map is symplectic. If no refraction root is found
    !> within max_kick radians the map falls back to the energy-exact Level-0
    !> rescale, so every energetically allowed marker still crosses.
    integer, parameter :: max_newton = 40
    integer, parameter :: max_backtrack = 20
    integer, parameter :: sym_iters = 2
    real(dp), parameter :: newton_rtol = 1.0d-14
    real(dp), parameter :: max_kick = 0.3d0
    !> Kicked Level-1 reflection (CA15): a forbidden crossing of a
    !> pure-magnitude sheet is an equal-B relocation along the interface at
    !> frozen (v_par, mu, H): the state follows the drift direction
    !> sign(ro0)*(h_zeta, -h_theta) until the home field returns to the entry
    !> value B* (same-side conjugate exit) or the target field drops to B*
    !> (relocated transmission). reloc_max_dirjump bounds the direction jump
    !> |h_home - h_target| inside which the pure-magnitude limit applies (the
    !> hybrid switch criterion scale); outside it, and on any trace failure,
    !> the Level-0 mirror is kept.
    real(dp), parameter :: reloc_step = 6.283185307179586d-3
    integer, parameter :: reloc_max_steps = 20000
    real(dp), parameter :: reloc_max_dirjump = 0.05d0

    type :: crossing_info_t
        integer :: event_type = 0
        integer :: iface = 0
        integer :: vol_from = 0
        integer :: vol_to = 0
        real(dp) :: theta = 0.0_dp
        real(dp) :: zeta = 0.0_dp
        real(dp) :: vpar_before = 0.0_dp
        real(dp) :: vpar_after = 0.0_dp
        real(dp) :: mu = 0.0_dp
        real(dp) :: bmod_home = 0.0_dp
        real(dp) :: bmod_target = 0.0_dp
        !> Level-1 generator X = {z, rho_g}, scalar lambda_k, and the applied
        !> tangential kicks (theta, zeta receive lambda_k*X); zero at Level 0.
        real(dp) :: xtheta = 0.0_dp
        real(dp) :: xzeta = 0.0_dp
        real(dp) :: xvpar = 0.0_dp
        real(dp) :: lambda = 0.0_dp
        real(dp) :: dtheta = 0.0_dp
        real(dp) :: dzeta = 0.0_dp
    end type crossing_info_t

    type :: crossing_record_t
        integer :: ipart = 0
        real(dp) :: time = 0.0_dp
        type(crossing_info_t) :: info
    end type crossing_record_t

    type :: thread_crossing_log_t
        type(crossing_record_t), allocatable :: record(:)
        integer :: count = 0
    end type thread_crossing_log_t

    type(thread_crossing_log_t), allocatable :: thread_log(:)
    type(crossing_record_t), allocatable :: crossing_log(:)
    integer :: crossing_count = 0

contains

    subroutine apply_crossing(y_iface, iface, direction, mvol, level, y_out, info)
        !> Map the guiding-center state y = (rho_g, theta, zeta, p, lambda) landed
        !> on interface rho_g = iface across to the neighbour volume. direction is
        !> +1 outward (home volume below the interface) or -1 inward (home above).
        real(dp), intent(in) :: y_iface(5)
        integer, intent(in) :: iface, direction, mvol, level
        real(dp), intent(out) :: y_out(5)
        type(crossing_info_t), intent(out) :: info

        real(dp) :: rho_face, rho_home, rho_target
        real(dp) :: bmod_home, bmod_target, perp_inv, mu, vpar, radicand, pitch

        info%event_type = CROSS_INVALID
        y_out = y_iface

        if (level /= CROSSING_LEVEL0 .and. level /= CROSSING_LEVEL1) then
            error stop 'apply_crossing: unknown crossing_level'
        end if
        if ((direction /= -1 .and. direction /= 1) .or. iface < 0 .or. &
                iface > mvol .or. .not. all_finite(y_iface) .or. &
                y_iface(4) <= 0.0_dp .or. abs(y_iface(5)) > 1.0_dp) return

        rho_face = real(iface, dp)
        pitch = max(-1.0_dp, min(1.0_dp, y_iface(5)))
        y_out(5) = pitch
        vpar = y_iface(4)*pitch

        info%iface = iface
        info%theta = y_iface(2)
        info%zeta = y_iface(3)
        info%vpar_before = vpar
        info%vpar_after = vpar

        if (direction == 1) then
            info%vol_from = iface
            info%vol_to = iface + 1
            rho_home = rho_face - iface_eps
            rho_target = rho_face + iface_eps
        else
            info%vol_from = iface + 1
            info%vol_to = iface
            rho_home = rho_face + iface_eps
            rho_target = rho_face - iface_eps
        end if

        bmod_home = bmod_at([rho_home, y_iface(2), y_iface(3)])
        if (.not. finite_value(bmod_home) .or. bmod_home <= 0.0_dp) return
        info%bmod_home = bmod_home
        info%bmod_target = bmod_home

        ! perp_inv = v_perp^2/|B| = z(4)^2 (1 - z(5)^2)/|B| is the perpendicular
        ! invariant the RK45 path stores; mu = perp_inv/2 is the magnetic moment
        ! held fixed across the interface, so the whole map runs on this invariant.
        perp_inv = y_iface(4)**2*max(1.0_dp - pitch**2, 0.0_dp)/bmod_home
        mu = 0.5_dp*perp_inv
        if (.not. finite_value(mu) .or. mu < 0.0_dp) return
        info%mu = mu

        if (direction == 1 .and. iface == mvol) then
            info%event_type = CROSS_LOSS
            info%vol_to = mvol + 1
            return
        end if

        if (direction == -1 .and. iface == 0) then
            info%event_type = CROSS_REFLECTION
            info%vol_to = info%vol_from
            info%vpar_after = -vpar
            y_out(5) = -y_iface(5)
            y_out(1) = axis_offset
            return
        end if

        if (level == CROSSING_LEVEL1) then
            call level1_map(rho_face, rho_home, rho_target, direction, &
                bmod_home, mu, vpar, y_iface, y_out, info)
            if (.not. valid_mapped_crossing(y_out, info)) then
                info%event_type = CROSS_INVALID
                y_out = y_iface
            end if
            return
        end if

        bmod_target = bmod_at([rho_target, y_iface(2), y_iface(3)])
        if (.not. finite_value(bmod_target) .or. bmod_target <= 0.0_dp) return
        info%bmod_target = bmod_target

        radicand = vpar**2 - 2.0_dp*mu*(bmod_target - bmod_home)
        if (.not. finite_value(radicand)) return
        if (radicand >= 0.0_dp) then
            info%event_type = CROSS_CROSSING
            info%vpar_after = sign(sqrt(radicand), vpar)
            y_out(5) = info%vpar_after/y_iface(4)
            y_out(1) = rho_face
        else
            ! Forbidden crossing into the higher-field volume is a magnetic
            ! mirror: the marker stays home with v_par reversed (the same-side
            ! energy root), so |v_par|, mu and |B| are unchanged and H is exactly
            ! conserved while the orbit streams back off the sheet.
            info%event_type = CROSS_REFLECTION
            info%vol_to = info%vol_from
            info%vpar_after = -vpar
            y_out(5) = -y_iface(5)
            y_out(1) = rho_face
        end if
    end subroutine apply_crossing

    pure logical function valid_mapped_crossing(y, info)
        real(dp), intent(in) :: y(5)
        type(crossing_info_t), intent(in) :: info

        valid_mapped_crossing = all_finite(y) .and. y(4) > 0.0_dp .and. &
            abs(y(5)) <= 1.0_dp + 100.0_dp*epsilon(1.0_dp) .and. &
            finite_value(info%mu) .and. info%mu >= 0.0_dp .and. &
            finite_value(info%bmod_home) .and. info%bmod_home > 0.0_dp .and. &
            finite_value(info%bmod_target) .and. info%bmod_target > 0.0_dp .and. &
            (info%event_type == CROSS_CROSSING .or. &
             info%event_type == CROSS_REFLECTION)
    end function valid_mapped_crossing

    pure logical function finite_value(x)
        real(dp), intent(in) :: x
        integer(int64), parameter :: exponent_mask = int(z'7FF0000000000000', int64)
        integer(int64) :: bits

        bits = transfer(x, bits)
        finite_value = iand(bits, exponent_mask) /= exponent_mask
    end function finite_value

    pure logical function all_finite(x)
        real(dp), intent(in) :: x(:)
        integer :: i

        all_finite = .false.
        do i = 1, size(x)
            if (.not. finite_value(x(i))) return
        end do
        all_finite = .true.
    end function all_finite

    function bmod_at(x) result(bmod)
        real(dp), intent(in) :: x(3)
        real(dp) :: bmod

        real(dp) :: sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    end function bmod_at

    subroutine field_at(x, bmod, hcov, hctr)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, hcov(3), hctr(3)

        real(dp) :: sqrtg, bder(3), hcurl(3)

        call magfie(x, bmod, sqrtg, bder, hcov, hctr, hcurl)
    end subroutine field_at

    subroutine relocate_forbidden(rho_home, rho_target, theta0, zeta0, vpar, &
            mu, bmod_home, dtheta, dzeta, vpar_new, info)
        !> Equal-B relocation of a forbidden Level-1 crossing (CA15), with the
        !> Level-0 mirror as fallback. The map conserves H exactly because the
        !> exit field equals the frozen entry value B*; v_par, p, lambda, and
        !> mu are untouched, ownership is assigned explicitly per exit side,
        !> and no physical time is consumed.
        real(dp), intent(in) :: rho_home, rho_target, theta0, zeta0, vpar, mu
        real(dp), intent(in) :: bmod_home
        real(dp), intent(out) :: dtheta, dzeta, vpar_new
        type(crossing_info_t), intent(inout) :: info

        real(dp) :: bstar, sgn, hnorm, cosjump
        real(dp) :: bm, bp, gm, gp, gm_prev, gp_prev
        real(dp) :: th, ze, th_prev, ze_prev, th_hit, ze_hit
        real(dp) :: bm0, bp0, hcov(3), hctr(3), hcov_t(3), hctr_t(3)
        integer :: k
        logical :: reflect_hit, transmit_hit

        ! Level-0 mirror defaults; every early return keeps them.
        dtheta = 0.0_dp
        dzeta = 0.0_dp
        vpar_new = -vpar
        info%event_type = CROSS_REFLECTION
        info%vol_to = info%vol_from
        info%bmod_target = bmod_home

        call field_at([rho_home, theta0, zeta0], bm0, hcov, hctr)
        call field_at([rho_target, theta0, zeta0], bp0, hcov_t, hctr_t)
        ! h^i_home h_{target,i} is the invariant dot product of the two unit
        ! fields; 1 - cos bounds |[[h]]|^2/2.
        cosjump = dot_product(hctr, hcov_t)
        if (1.0_dp - cosjump > 0.5_dp*reloc_max_dirjump**2) return

        bstar = bm0
        sgn = sign(1.0_dp, ro0*mu)
        th = theta0
        ze = zeta0
        bm = bm0
        bp = bp0
        gm = 0.0_dp
        gp = bp0 - bstar
        reflect_hit = .false.
        transmit_hit = .false.

        do k = 1, reloc_max_steps
            th_prev = th
            ze_prev = ze
            gm_prev = gm
            gp_prev = gp
            hnorm = sqrt(hcov(2)**2 + hcov(3)**2)
            if (hnorm <= sqrt(tiny(1.0_dp))) return
            th = th_prev + sgn*reloc_step*hcov(3)/hnorm
            ze = ze_prev - sgn*reloc_step*hcov(2)/hnorm
            call field_at([rho_home, th, ze], bm, hcov, hctr)
            bp = bmod_at([rho_target, th, ze])
            gm = bm - bstar
            gp = bp - bstar
            if (k == 1) then
                ! The drift must carry the state into the sheet; a level line
                ! that exits immediately is the degenerate case the mirror owns.
                if (gm >= 0.0_dp) return
            end if
            if (gm >= 0.0_dp) then
                reflect_hit = .true.
                exit
            end if
            if (gp <= 0.0_dp) then
                transmit_hit = .true.
                exit
            end if
        end do

        if (reflect_hit) then
            call refine_root(rho_home, bstar, th_prev, ze_prev, th, ze, &
                th_hit, ze_hit)
            dtheta = th_hit - theta0
            dzeta = ze_hit - zeta0
            vpar_new = vpar
            info%event_type = CROSS_REFLECTION
            info%vol_to = info%vol_from
            info%bmod_target = bmod_at([rho_target, th_hit, ze_hit])
        else if (transmit_hit) then
            call refine_root(rho_target, bstar, th_prev, ze_prev, th, ze, &
                th_hit, ze_hit)
            dtheta = th_hit - theta0
            dzeta = ze_hit - zeta0
            vpar_new = vpar
            info%event_type = CROSS_CROSSING
            info%bmod_target = bmod_at([rho_target, th_hit, ze_hit])
        end if
    end subroutine relocate_forbidden

    subroutine refine_root(rho_side, bstar, th_a, ze_a, th_b, ze_b, th_hit, ze_hit)
        !> Bisection for B(rho_side, th, ze) = bstar on the segment between the
        !> last two trace points; the sign change is bracketed by construction.
        !> The full 60 halvings put the exit on the level set to round-off, so
        !> the relocation's H error stays at the energy-exactness floor of the
        !> other crossing branches instead of a solver tolerance.
        real(dp), intent(in) :: rho_side, bstar, th_a, ze_a, th_b, ze_b
        real(dp), intent(out) :: th_hit, ze_hit

        real(dp) :: sa, sb, sm, ga, gm_loc
        integer :: k

        sa = 0.0_dp
        sb = 1.0_dp
        ga = bmod_at([rho_side, th_a, ze_a]) - bstar
        do k = 1, 60
            sm = 0.5_dp*(sa + sb)
            gm_loc = bmod_at([rho_side, th_a + sm*(th_b - th_a), &
                              ze_a + sm*(ze_b - ze_a)]) - bstar
            if (gm_loc == 0.0_dp) then
                sa = sm
                sb = sm
                exit
            end if
            if (gm_loc*ga > 0.0_dp) then
                sa = sm
                ga = gm_loc
            else
                sb = sm
            end if
        end do
        sm = 0.5_dp*(sa + sb)
        th_hit = th_a + sm*(th_b - th_a)
        ze_hit = ze_a + sm*(ze_b - ze_a)
    end subroutine refine_root

    subroutine level1_map(rho_face, rho_home, rho_target, direction, &
            bmod_home, mu, vpar, y_iface, y_out, info)
        !> Level-1 refraction map: the impulse Delta z = lambda_k * X along the
        !> interface Hamiltonian vector field. The Level-0 energy criterion at the
        !> landing point (radicand0 = v_par^2 - 2 mu [[B]]) decides crossing vs
        !> mirror, so the event split matches Level 0; Level 1 refines the crossing
        !> with the tangential sheet-drift kick and the drift-order v_par term. The
        !> forbidden crossing reflects as the Level-0 mirror (same-side energy root,
        !> no kick), which keeps the trapped sheet-skimming class identical to
        !> Level 0.
        real(dp), intent(in) :: rho_face, rho_home, rho_target
        integer, intent(in) :: direction
        real(dp), intent(in) :: bmod_home, mu, vpar, y_iface(5)
        real(dp), intent(inout) :: y_out(5)
        type(crossing_info_t), intent(inout) :: info

        real(dp) :: theta0, zeta0, p, hh, radicand0
        real(dp) :: xt, xz, xv, lam, dtheta, dzeta, vpar_new

        theta0 = y_iface(2)
        zeta0 = y_iface(3)
        p = y_iface(4)
        hh = 0.5_dp*vpar**2 + mu*bmod_home
        radicand0 = 2.0_dp*(hh - mu*bmod_at([rho_target, theta0, zeta0]))

        xt = 0.0_dp
        xz = 0.0_dp
        xv = 0.0_dp
        lam = 0.0_dp
        dtheta = 0.0_dp
        dzeta = 0.0_dp

        if (radicand0 >= 0.0_dp) then
            call solve_crossing(rho_home, rho_target, theta0, zeta0, vpar, mu, hh, &
                radicand0, xt, xz, xv, lam, dtheta, dzeta, vpar_new)
            info%event_type = CROSS_CROSSING
            info%bmod_target = bmod_at([rho_target, theta0 + dtheta, zeta0 + dzeta])
            y_out(1) = rho_face
        else
            call relocate_forbidden(rho_home, rho_target, theta0, zeta0, vpar, &
                mu, bmod_home, dtheta, dzeta, vpar_new, info)
            y_out(1) = rho_face
        end if

        info%xtheta = xt
        info%xzeta = xz
        info%xvpar = xv
        info%lambda = lam
        info%dtheta = dtheta
        info%dzeta = dzeta
        info%vpar_after = vpar_new
        y_out(2) = theta0 + dtheta
        y_out(3) = zeta0 + dzeta
        y_out(5) = vpar_new/p
    end subroutine level1_map

    subroutine solve_crossing(rho_home, rho_target, theta0, zeta0, vpar, mu, hh, &
            radicand0, xt, xz, xv, lam, dtheta, dzeta, vpar_new)
        !> Damped Newton for the refraction root of
        !> F(lam) = 0.5*(v_par + lam*X_vpar)^2 + mu*B_target(kick) - H_home. The
        !> residual (via residual_at) refreshes the generator by symmetrization
        !> sweeps at the midpoint state z + 0.5*lam*X. Backtracking halves the step
        !> whenever it fails to reduce |F|, so the walk to the equal-|B| point does
        !> not diverge where F is non-monotonic in the poloidal angle. If |F| is not
        !> driven below tolerance within max_kick the map falls back to the
        !> energy-exact Level-0 rescale (no kick), so every energetically allowed
        !> marker crosses.
        real(dp), intent(in) :: rho_home, rho_target, theta0, zeta0, vpar, mu, hh
        real(dp), intent(in) :: radicand0
        real(dp), intent(out) :: xt, xz, xv, lam, dtheta, dzeta, vpar_new

        integer :: iter, bt
        real(dp) :: fres, dfres, step, lam_try, fres_try, dfres_try
        real(dp) :: xt_try, xz_try, xv_try
        logical :: converged

        lam = 0.0_dp
        call residual_at(rho_home, rho_target, theta0, zeta0, vpar, mu, hh, lam, &
            fres, dfres, xt, xz, xv)
        converged = abs(fres) <= newton_rtol*abs(hh)

        do iter = 1, max_newton
            if (converged) exit
            if (abs(dfres) <= tiny(1.0_dp)) exit
            step = fres/dfres
            do bt = 0, max_backtrack
                lam_try = lam - step
                call residual_at(rho_home, rho_target, theta0, zeta0, vpar, mu, &
                    hh, lam_try, fres_try, dfres_try, &
                    xt_try, xz_try, xv_try)
                if (abs(fres_try) < abs(fres)) exit
                step = 0.5_dp*step
            end do
            if (abs(fres_try) >= abs(fres)) exit
            lam = lam_try
            fres = fres_try
            dfres = dfres_try
            xt = xt_try
            xz = xz_try
            xv = xv_try
            converged = abs(fres) <= newton_rtol*abs(hh)
            if (abs(lam*xt) > max_kick .or. abs(lam*xz) > max_kick) exit
        end do

        if (converged .and. abs(lam*xt) <= max_kick &
            .and. abs(lam*xz) <= max_kick) then
            dtheta = lam*xt
            dzeta = lam*xz
            vpar_new = vpar + lam*xv
        else
            xt = 0.0_dp
            xz = 0.0_dp
            xv = 0.0_dp
            lam = 0.0_dp
            dtheta = 0.0_dp
            dzeta = 0.0_dp
            vpar_new = sign(sqrt(radicand0), vpar)
        end if
    end subroutine solve_crossing

    subroutine residual_at(rho_home, rho_target, theta0, zeta0, vpar, mu, hh, lam, &
            fres, dfres, xt, xz, xv)
        !> Energy residual F(lam) and its lambda-derivative for the refraction
        !> solve, with the generator X evaluated self-consistently at the midpoint
        !> state z + 0.5*lam*X by sym_iters fixed-point sweeps. dfres freezes X, the
        !> dominant term being the poloidal |B| gradient along the tangential kick.
        real(dp), intent(in) :: rho_home, rho_target, theta0, zeta0, vpar, mu, hh
        real(dp), intent(in) :: lam
        real(dp), intent(out) :: fres, dfres, xt, xz, xv

        integer :: k
        real(dp) :: th_mid, ze_mid, vp_mid, th_k, ze_k, vp_k
        real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        xt = 0.0_dp
        xz = 0.0_dp
        xv = 0.0_dp
        do k = 1, sym_iters
            th_mid = theta0 + 0.5_dp*lam*xt
            ze_mid = zeta0 + 0.5_dp*lam*xz
            vp_mid = vpar + 0.5_dp*lam*xv
            call generator_sym(rho_home, rho_target, th_mid, ze_mid, vp_mid, &
                xt, xz, xv)
        end do
        th_k = theta0 + lam*xt
        ze_k = zeta0 + lam*xz
        vp_k = vpar + lam*xv
        call magfie([rho_target, th_k, ze_k], bmod, sqrtg, bder, hcovar, hctrvr, &
            hcurl)
        fres = 0.5_dp*vp_k**2 + mu*bmod - hh
        dfres = vp_k*xv + mu*bmod*(bder(2)*xt + bder(3)*xz)
    end subroutine residual_at

    subroutine generator_sym(rho_home, rho_target, theta, zeta, vpar_mid, &
            xt, xz, xv)
        !> Interface Hamiltonian vector field X = {z, rho_g} in SIMPLE Gaussian
        !> units, from the guiding-center bracket. hcovar/sqrtg/hcurl are the
        !> midpoint average of the two volumes, and the drift constant ro0 and
        !> Bstar_par = bmod + ro0*v_par*(h.curl h) are the same ones velo_can uses
        !> for the per-volume drift (parmot_mod ro0; hpstar = Bstar_par/bmod).
        real(dp), intent(in) :: rho_home, rho_target, theta, zeta, vpar_mid
        real(dp), intent(out) :: xt, xz, xv

        real(dp) :: b_h, sg_h, bder_h(3), hcov_h(3), hctr_h(3), hcurl_h(3)
        real(dp) :: b_t, sg_t, bder_t(3), hcov_t(3), hctr_t(3), hcurl_t(3)
        real(dp) :: sqrtg, hcov(3), hcurl(3), bmod, bstar_par

        call magfie([rho_home, theta, zeta], b_h, sg_h, bder_h, hcov_h, hctr_h, &
            hcurl_h)
        call magfie([rho_target, theta, zeta], b_t, sg_t, bder_t, hcov_t, hctr_t, &
            hcurl_t)
        bmod = 0.5_dp*(b_h + b_t)
        sqrtg = 0.5_dp*(sg_h + sg_t)
        hcov = 0.5_dp*(hcov_h + hcov_t)
        hcurl = 0.5_dp*(hcurl_h + hcurl_t)
        bstar_par = bmod + ro0*vpar_mid*dot_product(hcov, hcurl)
        xt = -ro0*hcov(3)/(sqrtg*bstar_par)
        xz = ro0*hcov(2)/(sqrtg*bstar_par)
        xv = ro0*vpar_mid*hcurl(1)/bstar_par
    end subroutine generator_sym

    subroutine crossing_log_reset(capacity)
        integer, intent(in) :: capacity
        integer :: nthreads, tid

        nthreads = 1
        !$ nthreads = omp_get_max_threads()
        if (allocated(thread_log)) deallocate (thread_log)
        if (allocated(crossing_log)) deallocate (crossing_log)
        allocate (thread_log(0:nthreads - 1))
        do tid = 0, nthreads - 1
            allocate (thread_log(tid)%record(max(capacity, 1)))
            thread_log(tid)%count = 0
        end do
        crossing_count = 0
    end subroutine crossing_log_reset

    subroutine crossing_log_record(ipart, time, info)
        integer, intent(in) :: ipart
        real(dp), intent(in) :: time
        type(crossing_info_t), intent(in) :: info

        type(crossing_record_t), allocatable :: tmp(:)
        integer :: tid, n

        if (.not. allocated(thread_log)) return
        tid = 0
        !$ tid = omp_get_thread_num()
        n = thread_log(tid)%count
        if (n >= size(thread_log(tid)%record)) then
            allocate (tmp(2*size(thread_log(tid)%record)))
            tmp(1:n) = thread_log(tid)%record(1:n)
            call move_alloc(tmp, thread_log(tid)%record)
        end if
        n = n + 1
        thread_log(tid)%count = n
        thread_log(tid)%record(n)%ipart = ipart
        thread_log(tid)%record(n)%time = time
        thread_log(tid)%record(n)%info = info
    end subroutine crossing_log_record

    integer function crossing_log_count()
        integer :: tid

        crossing_log_count = 0
        if (.not. allocated(thread_log)) return
        do tid = lbound(thread_log, 1), ubound(thread_log, 1)
            crossing_log_count = crossing_log_count + thread_log(tid)%count
        end do
    end function crossing_log_count

    integer function crossing_log_count_type(event_type)
        integer, intent(in) :: event_type

        integer :: i, tid

        crossing_log_count_type = 0
        if (.not. allocated(thread_log)) return
        do tid = lbound(thread_log, 1), ubound(thread_log, 1)
            do i = 1, thread_log(tid)%count
                if (thread_log(tid)%record(i)%info%event_type == event_type) &
                    crossing_log_count_type = crossing_log_count_type + 1
            end do
        end do
    end function crossing_log_count_type

    subroutine crossing_log_write(filename)
        character(*), intent(in) :: filename

        integer :: unit, i
        integer, allocatable :: order(:)

        if (.not. allocated(thread_log)) return
        call gather_crossing_log
        call sorted_by_particle(order)
        open (newunit=unit, file=filename, recl=1024)
        write (unit, '(A)') '# particle  time  iface  type  vol_from  vol_to  '// &
            'theta  zeta  vpar_before  vpar_after  mu  bmod_home  bmod_target  '// &
            'xtheta  xzeta  xvpar  lambda  dtheta  dzeta'
        do i = 1, crossing_count
            associate (r => crossing_log(order(i)))
                write (unit, *) r%ipart, r%time, r%info%iface, r%info%event_type, &
                    r%info%vol_from, r%info%vol_to, r%info%theta, r%info%zeta, &
                    r%info%vpar_before, r%info%vpar_after, r%info%mu, &
                    r%info%bmod_home, r%info%bmod_target, &
                    r%info%xtheta, r%info%xzeta, r%info%xvpar, r%info%lambda, &
                    r%info%dtheta, r%info%dzeta
            end associate
        end do
        close (unit)
    end subroutine crossing_log_write

    subroutine gather_crossing_log
        integer :: first, last, tid

        crossing_count = crossing_log_count()
        if (allocated(crossing_log)) deallocate (crossing_log)
        allocate (crossing_log(max(crossing_count, 1)))
        first = 1
        do tid = lbound(thread_log, 1), ubound(thread_log, 1)
            last = first + thread_log(tid)%count - 1
            if (last >= first) crossing_log(first:last) = &
                thread_log(tid)%record(1:thread_log(tid)%count)
            first = last + 1
        end do
    end subroutine gather_crossing_log

    subroutine sorted_by_particle(order)
        !> Stable counting sort of the log by particle index. Each particle is
        !> traced by a single thread, so its events are appended in time order;
        !> grouping by particle therefore yields a thread-order-independent,
        !> reproducible file whose per-particle blocks stay time-ordered.
        integer, allocatable, intent(out) :: order(:)

        integer :: i, p, pmax
        integer, allocatable :: cnt(:), pos(:)

        allocate (order(crossing_count))
        if (crossing_count == 0) return

        pmax = 0
        do i = 1, crossing_count
            pmax = max(pmax, crossing_log(i)%ipart)
        end do

        allocate (cnt(pmax), pos(pmax))
        cnt = 0
        do i = 1, crossing_count
            p = crossing_log(i)%ipart
            cnt(p) = cnt(p) + 1
        end do
        pos(1) = 1
        do p = 2, pmax
            pos(p) = pos(p - 1) + cnt(p - 1)
        end do
        do i = 1, crossing_count
            p = crossing_log(i)%ipart
            order(pos(p)) = i
            pos(p) = pos(p) + 1
        end do
    end subroutine sorted_by_particle

end module interface_crossing
