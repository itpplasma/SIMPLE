program study_crossing
    !> Single sheet crossings: hybrid with SIMPLE's entry lift against the full
    !> orbit everywhere, compared as distributions over markers.
    !>
    !> usage: study_crossing.x SPEED N THETA0 DTHETA LAM0 DLAM [LIFT [EXIT [FILE]]]
    !>   SPEED   particle speed (gyroradius SPEED at |B| = 1)
    !>   N       markers; statistics are printed for the prefixes 256, 1024,
    !>           4096, ... up to N
    !>   THETA0, DTHETA  start poloidal angle of the GC, uniform in
    !>           THETA0 +- DTHETA/2 (rad)
    !>   LAM0, DLAM      start pitch v_par/v of the GC, uniform in LAM0 +- DLAM/2
    !>   LIFT    hybrid entry lift: uniform (default, SIMPLE's rule) |
    !>           uniform_avg (v_par re-solved to the gyroaveraged lifted p_phi) |
    !>           gyrophase (generic lift with the GC gyrophase, for comparison) |
    !>           avgentry (lift iterated so that the gyroaveraged GC mu of the
    !>           lifted particle is the GC mu) | avgentry_p (also its p_phi)
    !>   EXIT    hybrid exit: pinned (default, SIMPLE's rule: GC v_par from the
    !>           entry GC p_phi, mu from the energy) | plain (first-order) |
    !>           avgexit (moments averaged over one gyroperiod) | avgexit_p (and
    !>           v_par from the entry GC p_phi) | avgexit_x (and the GC moved
    !>           along grad psi to the entry GC p_phi, mu kept)
    !>   FILE    optional per-marker output
    !>
    !> Each marker starts as a GC D_START gyroradii inside the sheet
    !> (a_in = 0.25, a_out = 0.75). Hybrid: GC (RK45) until the GC is within
    !> D_IN of the sheet, LIFT to FO (RK45, exact branch switching), FO
    !> until its GC is farther than D_OUT, EXIT to GC. Reference: a GC from the
    !> same start distribution (independent draw, so that the two samples are
    !> independent) lifted at the start with a uniform gyrophase, FO
    !> throughout, with the same entry and exit tests on its GC. After entry
    !> and after the exit both continue FO for one gyroperiod, over which mu
    !> and the GC position are averaged.
    !> Recorded: transmitted (exit outside) or reflected; mu at entry and
    !> after the exit, over mu at the start (the GC value); mu after the exit
    !> over mu of the first FO gyroperiod (after the start for the full orbit,
    !> after entry for the hybrid), the same estimator at both ends; the
    !> poloidal arc of the GC from entry (first within D_IN) to the exit, in
    !> gyroradii; for the hybrid the GC-to-GC ratio mu of the EXIT GC over the
    !> entry GC mu, compared with the full orbit's mu_exit/mu_fo_start.
    !> Printed per prefix: transmitted fraction, mean and standard deviation of
    !> mu_exit/mu_start - 1 and of the arc with standard errors, and the
    !> two-sample Kolmogorov-Smirnov statistic with its asymptotic p-value.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use sheet_field, only: field_params_t, region_of, psi_derivatives, psi_value
    use rk_reference, only: rk_advance
    use gc_hybrid, only: gc_state_t, fo_to_gc, fo_to_gc_pinned, gc_to_fo, &
                         gc_advance, gc_pphi, gc_bmod, gc_from_moments, fo_velocity, &
                         entry_gyrophase, LIFT_GYROPHASE, LIFT_UNIFORM, LIFT_UNIFORM_AVG
    implicit none

    integer, parameter :: HYBRID = 1, FULL = 2
    integer, parameter :: MISSED = 0, TRANSMITTED = 1, REFLECTED = 2, STUCK = 3
    real(dp), parameter :: TWOPI = 6.283185307179586_dp
    real(dp), parameter :: D_START = 4.0_dp, D_IN = 3.0_dp, D_OUT = 4.0_dp
    real(dp), parameter :: RK_TOL = 1.0e-10_dp, GC_TOL = 1.0e-10_dp
    real(dp), parameter :: FO_SAMPLE = 0.1_dp, GC_CHUNK = 1.0_dp
    !> Substeps of the gyroaveraged switches and iteration limits of the entry.
    integer, parameter :: NGYRO = 32, MAX_ENTRY_ITER = 8, SHIFT_ITER = 3
    real(dp), parameter :: ENTRY_TOL = 1.0e-9_dp

    type(field_params_t) :: par
    character(256) :: arg, out_file, lift_name, exit_name
    real(dp) :: speed, theta0, dtheta, lam0, dlam, t_max
    real(dp), allocatable :: mu(:, :), mu_in(:, :), mu_fo(:, :), disp(:, :)
    real(dp), allocatable :: pjump(:), mu_gc(:)
    !> Change of psi at the GC per passage (hybrid: entry GC to exit GC; full:
    !> gyroaveraged first-order X after entry to after the exit).
    real(dp), allocatable :: dpsi(:, :)
    integer, allocatable :: status(:, :), nfail(:), niter(:)
    integer :: n, im, mode, unit, nprefix, lift
    logical :: pinned, avg_entry, pin_entry, avg_exit, shift_exit
    integer(8) :: clock0, clock1, rate

    call get_command_argument(1, arg)
    read (arg, *) speed
    call get_command_argument(2, arg)
    read (arg, *) n
    call get_command_argument(3, arg)
    read (arg, *) theta0
    call get_command_argument(4, arg)
    read (arg, *) dtheta
    call get_command_argument(5, arg)
    read (arg, *) lam0
    call get_command_argument(6, arg)
    read (arg, *) dlam
    lift_name = 'uniform'
    exit_name = 'pinned'
    out_file = ''
    if (command_argument_count() >= 7) call get_command_argument(7, lift_name)
    if (command_argument_count() >= 8) call get_command_argument(8, exit_name)
    if (command_argument_count() >= 9) call get_command_argument(9, out_file)
    pinned = trim(exit_name) == 'pinned' .or. trim(exit_name) == 'avgexit_p'
    avg_exit = index(exit_name, 'avgexit') == 1
    shift_exit = trim(exit_name) == 'avgexit_x'
    avg_entry = index(lift_name, 'avgentry') == 1
    pin_entry = trim(lift_name) == 'avgentry_p'
    select case (trim(lift_name))
    case ('gyrophase', 'avgentry', 'avgentry_p')
        lift = LIFT_GYROPHASE
    case ('uniform_avg')
        lift = LIFT_UNIFORM_AVG
    case default
        lift = LIFT_UNIFORM
    end select
    par%a_in = 0.25_dp
    par%a_out = 0.75_dp
    ! A few poloidal transits at the start pitch.
    t_max = 4.0_dp*TWOPI*par%r_major*2.0_dp/speed

    allocate (mu(n, 2), mu_in(n, 2), mu_fo(n, 2), disp(n, 2), status(n, 2))
    allocate (nfail(n), pjump(n), mu_gc(n), niter(n), dpsi(n, 2))
    dpsi = 0.0_dp
    niter = 0
    pjump = 0.0_dp
    mu_gc = 0.0_dp
    call system_clock(clock0, rate)
    !$omp parallel do default(shared) private(im, mode) schedule(dynamic) collapse(2)
    do im = 1, n
        do mode = HYBRID, FULL
            call run_marker(im, mode, status(im, mode), mu_fo(im, mode), &
                            mu_in(im, mode), &
                            mu(im, mode), disp(im, mode), nfail(im), pjump(im), &
                            mu_gc(im), niter(im))
        end do
    end do
    !$omp end parallel do
    call system_clock(clock1)

    print '(A,F5.3,A,F6.3,A,F5.3,A,F5.3,A,F5.3,4A,A,I0,A,F8.1)', 'v ', speed, &
        '  theta0 ', theta0, ' +- ', 0.5_dp*dtheta, '  lambda ', lam0, ' +- ', &
        0.5_dp*dlam, '  lift ', trim(lift_name), '  exit ', trim(exit_name), &
        '  fails ', sum(nfail), &
        '  wall s ', real(clock1 - clock0, dp)/real(rate, dp)
    print '(A,ES10.3,A,ES10.3,A,F6.2)', '  p_phi jump per round trip at entry: rms ', &
        sqrt(sum(pjump**2)/n), '  mean ', sum(pjump)/n, &
        '  extra gyroperiods per switch pair ', &
        real(sum(niter), dp)/n + merge(1.0_dp, 0.0_dp, avg_exit)
    print '(A)', '      N model   miss stuck  transmitted  mu1/mu0-1: mean'// &
        '              sd       disp/rho: mean                sd'
    nprefix = 256
    do while (nprefix <= n)
        call summarize(nprefix)
        nprefix = 4*nprefix
    end do
    if (len_trim(out_file) > 0) then
        open (newunit=unit, file=trim(out_file), status='replace')
        write (unit, '(A)') '# marker, then for hybrid and full: status '// &
            '(1 transmitted, 2 reflected), mu_entry/mu0, mu_exit/mu0, '// &
            'disp/rho; mu_fo_start/mu0 hybrid, full; hybrid mu_gc_exit/mu_gc_entry'
        do im = 1, n
            write (unit, '(I6,2(I3,3ES16.8),3ES16.8)') im, status(im, HYBRID), &
                mu_in(im, HYBRID), mu(im, HYBRID), disp(im, HYBRID), &
                status(im, FULL), mu_in(im, FULL), mu(im, FULL), disp(im, FULL), &
                mu_fo(im, HYBRID), mu_fo(im, FULL), mu_gc(im)
        end do
        close (unit)
    end if

contains

    real(dp) function uniform(im, k) result(u)
        !> Reproducible uniform number in [0, 1) for marker im, stream k.
        integer, intent(in) :: im, k

        u = entry_gyrophase([real(im, dp), real(k, dp)])/TWOPI
    end function uniform

    function start_gc(im, mode) result(gc)
        !> Independent draws for the two models keep the samples independent.
        integer, intent(in) :: im, mode
        type(gc_state_t) :: gc

        real(dp) :: theta, lam, r, bmod

        theta = theta0 + dtheta*(uniform(im, 3*mode - 2) - 0.5_dp)
        lam = lam0 + dlam*(uniform(im, 3*mode - 1) - 0.5_dp)
        r = par%r_sheet - D_START*speed
        gc%r = par%r_major + r*cos(theta)
        gc%z = r*sin(theta)
        ! |B| = sqrt(F^2 + |grad psi|^2)/R inside the sheet, psi = a_in r^2.
        bmod = hypot(par%f_tor, 2.0_dp*par%a_in*r)/gc%r
        gc%vpar = lam*speed
        gc%mu = 0.5_dp*speed**2*(1.0_dp - lam**2)/bmod
        gc%zeta = TWOPI*uniform(im, 3*mode)
    end function start_gc

    real(dp) function sheet_s(r, z) result(s)
        real(dp), intent(in) :: r, z

        s = hypot(r - par%r_major, z) - par%r_sheet
    end function sheet_s

    real(dp) function poloidal(r, z) result(theta)
        real(dp), intent(in) :: r, z

        theta = atan2(z, r - par%r_major)
    end function poloidal

    subroutine unwrap(theta_acc, theta_last, theta_new)
        real(dp), intent(inout) :: theta_acc, theta_last
        real(dp), intent(in) :: theta_new

        theta_acc = theta_acc + modulo(theta_new - theta_last + 0.5_dp*TWOPI, &
                                       TWOPI) - 0.5_dp*TWOPI
        theta_last = theta_new
    end subroutine unwrap

    subroutine run_marker(im, mode, stat, mu_first, mu_entry, mu_ratio, arc, &
                          fails, dpphi, mu_gc_ratio, iters)
        integer, intent(in) :: im, mode
        integer, intent(out) :: stat
        !> mu_first: mu over the first gyroperiod of the FO (after the start for
        !> the full orbit, after entry for the hybrid), measured as mu_ratio.
        real(dp), intent(out) :: mu_first, mu_entry, mu_ratio, arc
        integer, intent(inout) :: fails
        !> p_phi jump of one further FO -> GC -> FO round trip at entry (hybrid).
        real(dp), intent(inout) :: dpphi
        !> Hybrid: mu of the exit GC over the entry GC mu.
        real(dp), intent(inout) :: mu_gc_ratio
        !> Hybrid: FO gyroperiods spent by the iterated entry.
        integer, intent(inout) :: iters

        type(field_params_t) :: pm
        type(gc_state_t) :: gc, gc0
        real(dp) :: z(4), t, dt, pphi, h_try, s, th_acc, th_last, th_entry
        real(dp) :: pphi_gc
        real(dp) :: mu_avg, th_avg, q_avg(2), psi_in
        integer :: region, ncross
        logical :: approached, entered, ok

        pm = par
        pphi_gc = 0.0_dp
        gc0 = start_gc(im, mode)
        gc = gc0
        t = 0.0_dp
        th_last = poloidal(gc%r, gc%z)
        th_acc = 0.0_dp
        th_entry = 0.0_dp
        stat = MISSED
        mu_entry = 0.0_dp
        mu_first = 0.0_dp
        mu_ratio = 0.0_dp
        arc = 0.0_dp
        if (mode == HYBRID) then
            fails = 0
            entered = .false.
            do while (t < t_max .and. .not. entered)
                call gc_advance(pm, gc, GC_CHUNK, D_IN*speed, GC_TOL, dt, entered)
                t = t + dt
                call unwrap(th_acc, th_last, poloidal(gc%r, gc%z))
            end do
            if (.not. entered) return
            th_entry = th_acc
            pphi_gc = gc_pphi(pm, gc)
            psi_in = psi_at(pm, [gc%r, gc%z])
            if (avg_entry) then
                call averaged_lift(pm, gc, entry_gyrophase([gc%r, gc%z, gc%vpar, &
                                                            gc%mu]), z, pphi, iters, ok)
            else
                call gc_to_fo(pm, gc, lift, pphi_gc, z, pphi, ok)
            end if
            if (.not. ok) fails = fails + 1
            pm%p_phi = pphi
            dpphi = round_trip(pm, z, pphi, pphi_gc, TWOPI*uniform(im, 7))
            approached = .true.
        else
            call gc_to_fo(pm, gc, LIFT_GYROPHASE, 0.0_dp, z, pphi, ok)
            approached = .false.
        end if
        pm%p_phi = pphi
        region = region_of(pm, z(1:2))
        ncross = 0
        h_try = 0.01_dp
        stat = STUCK
        if (approached) then
            call gyro_average(pm, region, z, h_try, ncross, th_acc, th_last, &
                              mu_avg, th_avg, q_avg)
            mu_entry = mu_avg/gc0%mu
            mu_first = mu_entry
        else
            call gyro_average(pm, region, z, h_try, ncross, th_acc, th_last, &
                              mu_avg, th_avg, q_avg)
            mu_first = mu_avg/gc0%mu
        end if
        do while (t < t_max)
            call rk_advance(pm, region, z, FO_SAMPLE, RK_TOL, h_try, ncross)
            t = t + FO_SAMPLE
            gc = fo_to_gc(pm, z)
            call unwrap(th_acc, th_last, poloidal(gc%r, gc%z))
            s = sheet_s(gc%r, gc%z)
            if (.not. approached) then
                if (abs(s) < D_IN*speed) then
                    approached = .true.
                    th_entry = th_acc
                    call gyro_average(pm, region, z, h_try, ncross, th_acc, &
                                      th_last, mu_avg, th_avg, q_avg)
                    mu_entry = mu_avg/gc0%mu
                    psi_in = psi_at(pm, q_avg)
                end if
            else if (abs(s) > D_OUT*speed) then
                stat = merge(TRANSMITTED, REFLECTED, s > 0.0_dp)
                if (mode == HYBRID) then
                    gc = exit_gc(pm, z, pphi_gc, ok)
                    if (.not. ok) fails = fails + 1
                    mu_gc_ratio = gc%mu/gc0%mu
                    dpsi(im, HYBRID) = psi_at(pm, [gc%r, gc%z]) - psi_in
                end if
                exit
            end if
        end do
        if (.not. approached) stat = MISSED
        if (stat == MISSED .or. stat == STUCK) return
        call gyro_average(pm, region, z, h_try, ncross, th_acc, th_last, &
                          mu_avg, th_avg, q_avg)
        mu_ratio = mu_avg/gc0%mu
        if (mode == FULL) dpsi(im, FULL) = psi_at(pm, q_avg) - psi_in
        arc = par%r_sheet*(th_avg - th_entry)/speed
    end subroutine run_marker

    function exit_gc(pm, z, pphi_gc, ok) result(gc)
        !> The selected exit; the pinned one falls back to the plain one.
        type(field_params_t), intent(in) :: pm
        real(dp), intent(in) :: z(4), pphi_gc
        logical, intent(out) :: ok
        type(gc_state_t) :: gc

        ok = .true.
        if (avg_exit) then
            gc = averaged_gc(pm, z, pinned, pphi_gc, ok)
            if (shift_exit) call shift_to_pphi(pm, z, pphi_gc, gc)
            return
        end if
        if (pinned) gc = fo_to_gc_pinned(pm, z, pphi_gc, ok)
        if (.not. pinned .or. .not. ok) gc = fo_to_gc(pm, z)
    end function exit_gc

    function averaged_gc(pm, z, pin, pphi_gc, ok) result(gc)
        !> GC from a copy of the FO continued for one gyroperiod 2 pi/|B(X)|:
        !> first-order X, v_par and v_perp^2 averaged over NGYRO substeps.
        type(field_params_t), intent(in) :: pm
        real(dp), intent(in) :: z(4), pphi_gc
        logical, intent(in) :: pin
        logical, intent(out) :: ok
        type(gc_state_t) :: gc

        type(gc_state_t) :: g
        real(dp) :: zz(4), h, q(2), vperp2, vpar, spd, dt
        integer :: region, ncross, k

        zz = z
        region = region_of(pm, zz(1:2))
        h = 0.01_dp
        ncross = 0
        spd = norm2(fo_velocity(pm, zz))
        g = fo_to_gc(pm, zz)
        dt = TWOPI/gc_bmod(pm, [g%r, g%z])/NGYRO
        q = 0.0_dp
        vperp2 = 0.0_dp
        vpar = 0.0_dp
        do k = 1, NGYRO
            call rk_advance(pm, region, zz, dt, RK_TOL, h, ncross)
            g = fo_to_gc(pm, zz)
            q = q + [g%r, g%z]/NGYRO
            vpar = vpar + g%vpar/NGYRO
            vperp2 = vperp2 + (spd*spd - g%vpar**2)/NGYRO
        end do
        gc = gc_from_moments(pm, q, vperp2, vpar, spd, pin, pphi_gc, ok)
        gc%zeta = g%zeta
    end function averaged_gc

    subroutine shift_to_pphi(pm, z, target, gc)
        !> Move the GC along grad psi so that its p_phi equals target, mu fixed
        !> and v_par from the particle's energy at the new position, sign kept.
        type(field_params_t), intent(in) :: pm
        real(dp), intent(in) :: z(4), target
        type(gc_state_t), intent(inout) :: gc

        real(dp) :: spd, d, psi, gp(2), hp(2, 2)
        integer :: it

        spd = norm2(fo_velocity(pm, z))
        do it = 1, SHIFT_ITER
            d = target - gc_pphi(pm, gc)
            call psi_derivatives(pm, region_of(pm, [gc%r, gc%z]), [gc%r, gc%z], &
                                 psi, gp, hp)
            gc%r = gc%r + d*gp(1)/sum(gp**2)
            gc%z = gc%z + d*gp(2)/sum(gp**2)
            gc%vpar = sign(sqrt(max(spd*spd - 2.0_dp*gc%mu* &
                                    gc_bmod(pm, [gc%r, gc%z]), 0.0_dp)), gc%vpar)
        end do
    end subroutine shift_to_pphi

    subroutine averaged_lift(pm, gc, alpha, z, pphi, iters, ok)
        !> Lift at the phase alpha, iterated so that averaged_gc of the lifted
        !> particle returns the GC mu (and, with pin_entry, the GC p_phi, by
        !> shifting the lift position along grad psi). Each trial costs one
        !> gyroperiod of FO, discarded.
        type(field_params_t), intent(in) :: pm
        type(gc_state_t), intent(in) :: gc
        real(dp), intent(in) :: alpha
        real(dp), intent(out) :: z(4), pphi
        integer, intent(inout) :: iters
        logical, intent(out) :: ok

        type(field_params_t) :: pt
        type(gc_state_t) :: gk, gbar
        real(dp) :: spd, mu_k, shift, target, psi, gp(2), hp(2, 2), bk, dmu, dp_
        integer :: it
        logical :: ok_avg

        spd = sqrt(gc%vpar**2 + 2.0_dp*gc%mu*gc_bmod(pm, [gc%r, gc%z]))
        target = gc_pphi(pm, gc)
        call psi_derivatives(pm, region_of(pm, [gc%r, gc%z]), [gc%r, gc%z], &
                             psi, gp, hp)
        mu_k = gc%mu
        shift = 0.0_dp
        do it = 1, MAX_ENTRY_ITER
            gk = gc
            gk%r = gc%r + shift*gp(1)/sum(gp**2)
            gk%z = gc%z + shift*gp(2)/sum(gp**2)
            bk = gc_bmod(pm, [gk%r, gk%z])
            gk%mu = mu_k
            gk%vpar = sign(sqrt(max(spd*spd - 2.0_dp*mu_k*bk, 0.0_dp)), gc%vpar)
            gk%zeta = alpha
            call gc_to_fo(pm, gk, LIFT_GYROPHASE, 0.0_dp, z, pphi, ok)
            pt = pm
            pt%p_phi = pphi
            iters = iters + 1
            gbar = averaged_gc(pt, z, .false., 0.0_dp, ok_avg)
            dmu = gbar%mu/gc%mu - 1.0_dp
            dp_ = target - gc_pphi(pt, gbar)
            if (abs(dmu) < ENTRY_TOL .and. (.not. pin_entry .or. &
                                  abs(dp_) < ENTRY_TOL*max(abs(target), 1.0_dp))) return
            mu_k = mu_k/(1.0_dp + dmu)
            if (pin_entry) shift = shift + dp_
        end do
        ok = .false.
    end subroutine averaged_lift

    real(dp) function round_trip(pm, z, pphi, pphi_gc, zeta) result(dp_phi)
        !> FO -> GC -> FO with the selected exit and lift: change of the FO
        !> p_phi. The GC gyrophase is set to zeta, as after a GC segment of
        !> arbitrary length.
        type(field_params_t), intent(in) :: pm
        real(dp), intent(in) :: z(4), pphi, pphi_gc, zeta

        type(gc_state_t) :: gc
        real(dp) :: z2(4), pphi2
        integer :: iters
        logical :: ok

        gc = exit_gc(pm, z, pphi_gc, ok)
        gc%zeta = zeta
        if (avg_entry) then
            iters = 0
            call averaged_lift(pm, gc, zeta, z2, pphi2, iters, ok)
        else
            call gc_to_fo(pm, gc, lift, gc_pphi(pm, gc), z2, pphi2, ok)
        end if
        dp_phi = pphi2 - pphi
    end function round_trip

    real(dp) function psi_at(pm, q) result(psi)
        type(field_params_t), intent(in) :: pm
        real(dp), intent(in) :: q(2)

        psi = psi_value(pm, region_of(pm, q), q)
    end function psi_at

    subroutine gyro_average(pm, region, z, h_try, ncross, th_acc, th_last, &
                            mu_avg, th_avg, q_avg)
        !> Advance the FO by one gyroperiod 2 pi R/F and average the GC mu and
        !> poloidal angle over it.
        type(field_params_t), intent(in) :: pm
        integer, intent(inout) :: region, ncross
        real(dp), intent(inout) :: z(4), h_try, th_acc, th_last
        real(dp), intent(out) :: mu_avg, th_avg, q_avg(2)

        type(gc_state_t) :: gc
        integer :: k, nsample

        nsample = max(1, nint(TWOPI*z(1)/pm%f_tor/FO_SAMPLE))
        mu_avg = 0.0_dp
        th_avg = 0.0_dp
        q_avg = 0.0_dp
        do k = 1, nsample
            call rk_advance(pm, region, z, FO_SAMPLE, RK_TOL, h_try, ncross)
            gc = fo_to_gc(pm, z)
            call unwrap(th_acc, th_last, poloidal(gc%r, gc%z))
            mu_avg = mu_avg + gc%mu/nsample
            th_avg = th_avg + th_acc/nsample
            q_avg = q_avg + [gc%r, gc%z]/nsample
        end do
    end subroutine gyro_average

    subroutine summarize(m)
        integer, intent(in) :: m

        character(6), parameter :: label(2) = ['hybrid', 'full  ']
        logical :: done(m, 2)
        real(dp) :: ks_mu(2), ks_arc(2), frac, mean_mu, sd_mu, mean_arc, sd_arc
        integer :: md, nd

        done = status(1:m, :) == TRANSMITTED .or. status(1:m, :) == REFLECTED
        do md = HYBRID, FULL
            nd = count(done(:, md))
            frac = real(count(status(1:m, md) == TRANSMITTED), dp)/max(nd, 1)
            call moments(pack(mu(1:m, md) - 1.0_dp, done(:, md)), mean_mu, sd_mu)
            call moments(pack(disp(1:m, md), done(:, md)), mean_arc, sd_arc)
            print '(I7,1X,A6,2I6,F7.3,"+-",F5.3,'// &
                '2(ES11.3,"+-",ES8.2,ES10.3,"+-",ES8.2))', m, label(md), &
                count(status(1:m, md) == MISSED), &
                count(status(1:m, md) == STUCK), frac, &
                sqrt(frac*(1.0_dp - frac)/max(nd, 1)), &
                mean_mu, sd_mu/sqrt(real(max(nd, 1), dp)), sd_mu, &
                sd_mu/sqrt(2.0_dp*max(nd - 1, 1)), &
                mean_arc, sd_arc/sqrt(real(max(nd, 1), dp)), sd_arc, &
                sd_arc/sqrt(2.0_dp*max(nd - 1, 1))
        end do
        call ks_two_sample(pack(mu(1:m, 1), done(:, 1)), &
                           pack(mu(1:m, 2), done(:, 2)), ks_mu(1), ks_mu(2))
        call ks_two_sample(pack(disp(1:m, 1), done(:, 1)), &
                           pack(disp(1:m, 2), done(:, 2)), ks_arc(1), ks_arc(2))
        print '(I7,A,2F8.4,A,2F8.4)', m, '  KS mu: D, p ', ks_mu, &
            '   KS disp: D, p ', ks_arc
        ! Same estimator at both ends: the lift artefact cancels.
        do md = HYBRID, FULL
            call moments(pack(mu(1:m, md)/mu_fo(1:m, md) - 1.0_dp, done(:, md)), &
                         mean_mu, sd_mu)
            print '(I7,1X,A6,A,ES11.3,"+-",ES8.2,ES10.3,"+-",ES8.2)', m, label(md), &
                '  mu_exit/mu_fo_start-1: mean, sd ', mean_mu, &
                sd_mu/sqrt(real(max(count(done(:, md)), 1), dp)), sd_mu, &
                sd_mu/sqrt(2.0_dp*max(count(done(:, md)) - 1, 1))
        end do
        call ks_two_sample(pack(mu(1:m, 1)/mu_fo(1:m, 1), done(:, 1)), &
                           pack(mu(1:m, 2)/mu_fo(1:m, 2), done(:, 2)), &
                           ks_mu(1), ks_mu(2))
        print '(I7,A,2F8.4)', m, '  KS mu_exit/mu_fo_start: D, p ', ks_mu
        call moments(pack(mu_gc(1:m) - 1.0_dp, done(:, 1)), mean_mu, sd_mu)
        print '(I7,1X,A6,A,ES11.3,"+-",ES8.2,ES10.3,"+-",ES8.2)', m, label(1), &
            '  mu_gc_exit/mu_gc_entry-1: mean, sd ', mean_mu, &
            sd_mu/sqrt(real(max(count(done(:, 1)), 1), dp)), sd_mu, &
            sd_mu/sqrt(2.0_dp*max(count(done(:, 1)) - 1, 1))
        call ks_two_sample(pack(mu_gc(1:m), done(:, 1)), &
                           pack(mu(1:m, 2)/mu_fo(1:m, 2), done(:, 2)), &
                           ks_mu(1), ks_mu(2))
        print '(I7,A,2F8.4)', m, &
            '  KS hybrid GC-to-GC vs full mu_exit/mu_fo_start: D, p ', ks_mu
        do md = HYBRID, FULL
            call moments(pack(dpsi(1:m, md), done(:, md)), mean_mu, sd_mu)
            print '(I7,1X,A6,A,ES11.3,"+-",ES8.2,ES10.3)', m, label(md), &
                '  dpsi per passage: mean, sd ', mean_mu, &
                sd_mu/sqrt(real(max(count(done(:, md)), 1), dp)), sd_mu
        end do
    end subroutine summarize

    subroutine moments(x, mean, sd)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: mean, sd

        mean = sum(x)/max(size(x), 1)
        sd = sqrt(sum((x - mean)**2)/max(size(x) - 1, 1))
    end subroutine moments

    subroutine ks_two_sample(x, y, d, p)
        !> Two-sample Kolmogorov-Smirnov statistic and asymptotic p-value with
        !> the effective-size correction of Stephens (1970).
        real(dp), intent(in) :: x(:), y(:)
        real(dp), intent(out) :: d, p

        real(dp) :: xs(size(x)), ys(size(y)), v, ne, lam, term
        integer :: i, j, nx, ny, k

        nx = size(x)
        ny = size(y)
        d = 0.0_dp
        p = 1.0_dp
        if (nx == 0 .or. ny == 0) return
        xs = x
        ys = y
        call sort(xs)
        call sort(ys)
        i = 1
        j = 1
        do while (i <= nx .and. j <= ny)
            v = min(xs(i), ys(j))
            do while (i <= nx)
                if (xs(i) > v) exit
                i = i + 1
            end do
            do while (j <= ny)
                if (ys(j) > v) exit
                j = j + 1
            end do
            d = max(d, abs(real(i - 1, dp)/nx - real(j - 1, dp)/ny))
        end do
        ne = real(nx, dp)*ny/(nx + ny)
        lam = (sqrt(ne) + 0.12_dp + 0.11_dp/sqrt(ne))*d
        p = 0.0_dp
        do k = 1, 100
            term = 2.0_dp*(-1.0_dp)**(k - 1)*exp(-2.0_dp*k*k*lam*lam)
            p = p + term
            if (abs(term) < 1.0e-12_dp) exit
        end do
        p = min(max(p, 0.0_dp), 1.0_dp)
        if (lam < 0.2_dp) p = 1.0_dp
    end subroutine ks_two_sample

    subroutine sort(a)
        !> Shell sort, ascending.
        real(dp), intent(inout) :: a(:)

        integer :: gap, i, j
        real(dp) :: tmp

        gap = size(a)/2
        do while (gap > 0)
            do i = gap + 1, size(a)
                tmp = a(i)
                j = i
                do while (j > gap)
                    if (a(j - gap) <= tmp) exit
                    a(j) = a(j - gap)
                    j = j - gap
                end do
                a(j) = tmp
            end do
            gap = gap/2
        end do
    end subroutine sort

end program study_crossing
