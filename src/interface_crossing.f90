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

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use magfie_sub, only: magfie
    use parmot_mod, only: ro0

    implicit none
    private

    public :: apply_crossing, crossing_info_t, axis_offset
    public :: crossing_log_reset, crossing_log_record, crossing_log_write, &
              crossing_log_count
    public :: CROSSING_LEVEL0, CROSSING_LEVEL1, CROSS_CROSSING, CROSS_REFLECTION, &
              CROSS_LOSS, CROSS_STOP

    integer, parameter :: CROSSING_LEVEL0 = 0
    integer, parameter :: CROSSING_LEVEL1 = 1
    integer, parameter :: CROSS_CROSSING = 1
    integer, parameter :: CROSS_REFLECTION = 2
    integer, parameter :: CROSS_LOSS = 3
    !> Symplectic orbits terminate at volume boundaries until the exact-landing
    !> crossing is wired (SIMPLE#441); their stop events share this log.
    integer, parameter :: CROSS_STOP = 4

    !> Inner cutoff for the innermost volume: rho_g = 0 is the coordinate axis
    !> where sqrt(g) = 0, so the marker reflects trivially before reaching it
    !> rather than crossing a (non-existent) inner interface.
    real(dp), parameter :: axis_offset = 1.0d-3

    !> Both-side |B| is read at rho_g = k -+ iface_eps: int(k -+ 1e-12) selects
    !> the neighbouring volume k-1 resp. k, and the libneo polynomial basis of
    !> each volume is well-defined at that offset from the shared interface point.
    real(dp), parameter :: iface_eps = 1.0d-12
    !> After the map rho_g is nudged reflect_nudge off the interface so the
    !> restarted RK45 step does not immediately re-detect the same event. It must
    !> clear the odeint event tolerance (~sqrt(epsilon) ~ 1.5e-8) or the next
    !> substep triggers a zero-length event at its start point.
    real(dp), parameter :: reflect_nudge = 1.0d-6

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
        real(dp) :: bmod_home, bmod_target, perp_inv, mu, vpar, radicand

        if (level /= CROSSING_LEVEL0 .and. level /= CROSSING_LEVEL1) then
            error stop 'apply_crossing: unknown crossing_level'
        end if

        rho_face = real(iface, dp)
        y_out = y_iface
        vpar = y_iface(4)*y_iface(5)

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
        info%bmod_home = bmod_home
        info%bmod_target = bmod_home

        ! perp_inv = v_perp^2/|B| = z(4)^2 (1 - z(5)^2)/|B| is the perpendicular
        ! invariant the RK45 path stores; mu = perp_inv/2 is the magnetic moment
        ! held fixed across the interface, so the whole map runs on this invariant.
        perp_inv = y_iface(4)**2*(1.0_dp - y_iface(5)**2)/bmod_home
        mu = 0.5_dp*perp_inv
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
            y_out(1) = axis_offset + reflect_nudge
            return
        end if

        if (level == CROSSING_LEVEL1) then
            call level1_map(rho_face, rho_home, rho_target, direction, &
                            bmod_home, mu, vpar, y_iface, y_out, info)
            return
        end if

        bmod_target = bmod_at([rho_target, y_iface(2), y_iface(3)])
        info%bmod_target = bmod_target

        radicand = vpar**2 - 2.0_dp*mu*(bmod_target - bmod_home)
        if (radicand >= 0.0_dp) then
            info%event_type = CROSS_CROSSING
            info%vpar_after = sign(sqrt(radicand), vpar)
            y_out(5) = info%vpar_after/y_iface(4)
            y_out(1) = rho_face + real(direction, dp)*reflect_nudge
        else
            ! Forbidden crossing into the higher-field volume is a magnetic
            ! mirror: the marker stays home with v_par reversed (the same-side
            ! energy root), so |v_par|, mu and |B| are unchanged and H is exactly
            ! conserved while the orbit streams back off the sheet.
            info%event_type = CROSS_REFLECTION
            info%vol_to = info%vol_from
            info%vpar_after = -vpar
            y_out(5) = -y_iface(5)
            y_out(1) = rho_face - real(direction, dp)*reflect_nudge
        end if
    end subroutine apply_crossing

    function bmod_at(x) result(bmod)
        real(dp), intent(in) :: x(3)
        real(dp) :: bmod

        real(dp) :: sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    end function bmod_at

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
            y_out(1) = rho_face + real(direction, dp)*reflect_nudge
        else
            vpar_new = -vpar
            info%event_type = CROSS_REFLECTION
            info%vol_to = info%vol_from
            info%bmod_target = bmod_home
            y_out(1) = rho_face - real(direction, dp)*reflect_nudge
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

        if (allocated(crossing_log)) deallocate (crossing_log)
        allocate (crossing_log(max(capacity, 1)))
        crossing_count = 0
    end subroutine crossing_log_reset

    subroutine crossing_log_record(ipart, time, info)
        integer, intent(in) :: ipart
        real(dp), intent(in) :: time
        type(crossing_info_t), intent(in) :: info

        type(crossing_record_t), allocatable :: tmp(:)

        !$omp critical (spectre_crossing_log)
        if (.not. allocated(crossing_log)) allocate (crossing_log(64))
        if (crossing_count >= size(crossing_log)) then
            allocate (tmp(2*size(crossing_log)))
            tmp(1:crossing_count) = crossing_log(1:crossing_count)
            call move_alloc(tmp, crossing_log)
        end if
        crossing_count = crossing_count + 1
        crossing_log(crossing_count)%ipart = ipart
        crossing_log(crossing_count)%time = time
        crossing_log(crossing_count)%info = info
        !$omp end critical (spectre_crossing_log)
    end subroutine crossing_log_record

    integer function crossing_log_count()
        crossing_log_count = crossing_count
    end function crossing_log_count

    subroutine crossing_log_write(filename)
        character(*), intent(in) :: filename

        integer :: unit, i
        integer, allocatable :: order(:)

        if (.not. allocated(crossing_log)) return

        ! Checkpoint dumps call this from a worker thread while others append, so
        ! serialize against crossing_log_record to keep the growable array live.
        !$omp critical (spectre_crossing_log)
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
        !$omp end critical (spectre_crossing_log)
    end subroutine crossing_log_write

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
