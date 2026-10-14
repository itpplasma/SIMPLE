module candidate_orbits
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use candidate_field, only: model_t, local_t, geometry, field, cross, frame, &
        distance, side_of, normal, radius, canonical_pphi, connection
    use transforms, only: lift, inverse, gc_energy, order2_rhs, kernel_available
    implicit none
    private
    public :: options_t, orbit_t, initialize, advance, particle, helix, &
        plane_step, gc_rhs, dp_step, reset_old
    real(dp), parameter :: PI = acos(-1.0_dp)

    type :: options_t
        ! 0=reference, 1=old GC, 2=discarded phase, 3=carried phase,
        ! 4=carried phase with frozen local helix patches, 5=Boris layer.
        integer :: method = 0, order = 0, old_level = 0, phase_rule = 0
        real(dp) :: tol = 1.0e-8_dp, din = 3.0_dp, dout = 4.0_dp
        real(dp) :: patch_step = 0.5_dp, boris_step = 0.025_dp
        real(dp) :: sample = 0.5_dp
    end type

    type :: orbit_t
        real(dp) :: y(6) = 0.0_dp, h = 0.05_dp, speed = 0.0_dp
        real(dp) :: p0 = 0.0_dp, energy0 = 0.0_dp, time = 0.0_dp
        real(dp) :: fo_time = 0.0_dp, switch_energy = 0.0_dp
        logical :: fo = .true., failed = .false.
        integer :: region = 1, switches = 0, crossings = 0, roots_failed = 0
        integer :: steps = 0, evaluations = 0
    end type

contains

    subroutine initialize(pm, opt, y, orbit)
        type(model_t), intent(in) :: pm
        type(options_t), intent(in) :: opt
        real(dp), intent(in) :: y(6)
        type(orbit_t), intent(out) :: orbit
        real(dp) :: residual
        logical :: ok
        orbit = orbit_t()
        orbit%y = y
        orbit%speed = norm2(y(4:6))
        orbit%energy0 = 0.5_dp*orbit%speed**2
        orbit%p0 = canonical_pphi(pm, y)
        orbit%region = side_of(pm, y(1:3))
        if (opt%method == 3 .or. opt%method == 4) then
            if (.not. kernel_available(opt%order)) then
                orbit%failed = .true.
                return
            end if
        end if
        if (opt%method == 1) then
            call inverse(pm, y, 0, orbit%y, ok, residual)
            orbit%failed = .not. ok
            orbit%fo = .false.
            orbit%region = side_of(pm, orbit%y(1:3))
        end if
    end subroutine

    subroutine particle(pm, opt, orbit, y, ok)
        type(model_t), intent(in) :: pm
        type(options_t), intent(in) :: opt
        type(orbit_t), intent(in) :: orbit
        real(dp), intent(out) :: y(6)
        logical, intent(out) :: ok
        y = orbit%y
        ok = .true.
        if (.not. orbit%fo) call lift(pm, orbit%y, opt%order, y, ok)
    end subroutine

    function gc_rhs(pm, g, region, order) result(f)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: g(6)
        integer, intent(in) :: region, order
        real(dp) :: f(6), bst(3), bp, rs(3), cr(3)
        type(local_t) :: lf
        logical :: ok
        if (order == 2) then
            if (.not. kernel_available(order)) then
                f = huge(1.0_dp)
                return
            end if
            call order2_rhs(pm, g, f, ok)
            if (.not. ok) f = huge(1.0_dp)
            return
        end if
        call geometry(pm, g(1:3), region, lf)
        rs = lf%gauge
        cr = 0.0_dp
        if (order == 1) call connection(pm, g(1:3), region, rs, cr)
        bst = lf%bvec + g(4)*lf%curlb - g(5)*cr
        bp = dot_product(lf%b, bst)
        f(1:3) = (g(4)*bst + g(5)*cross(lf%b, lf%gradb))/bp
        f(4) = -g(5)*dot_product(bst, lf%gradb)/bp
        f(5) = 0.0_dp
        f(6) = -lf%bmag - dot_product(rs, f(1:3))
    end function

    subroutine dp_step(pm, y, region, fo, order, h, speed, yn, err)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: y(6), h, speed
        integer, intent(in) :: region, order
        logical, intent(in) :: fo
        real(dp), intent(out) :: yn(6), err
        real(dp) :: k(6, 7), ev(6), scale(6)
        k(:, 1) = rhs(y)
        k(:, 2) = rhs(y + h*k(:, 1)/5)
        k(:, 3) = rhs(y + h*(3*k(:, 1) + 9*k(:, 2))/40)
        k(:, 4) = rhs(y + h*(44*k(:, 1)/45 - 56*k(:, 2)/15 &
                          + 32*k(:, 3)/9))
        k(:, 5) = rhs(y + h*(19372*k(:, 1)/6561 - 25360*k(:, 2)/2187 &
                          + 64448*k(:, 3)/6561 - 212*k(:, 4)/729))
        k(:, 6) = rhs(y + h*(9017*k(:, 1)/3168 - 355*k(:, 2)/33 &
                          + 46732*k(:, 3)/5247 + 49*k(:, 4)/176 &
                          - 5103*k(:, 5)/18656))
        yn = y + h*(35*k(:, 1)/384 + 500*k(:, 3)/1113 &
                  + 125*k(:, 4)/192 - 2187*k(:, 5)/6784 + 11*k(:, 6)/84)
        k(:, 7) = rhs(yn)
        ev = h*(71*k(:, 1)/57600 - 71*k(:, 3)/16695 + 71*k(:, 4)/1920 &
               - 17253*k(:, 5)/339200 + 22*k(:, 6)/525 - k(:, 7)/40)
        scale = speed
        if (.not. fo) then
            scale(5) = speed**2
            scale(6) = 1.0_dp
        end if
        err = maxval(abs(ev)/scale)
    contains
        function rhs(yy) result(f)
            real(dp), intent(in) :: yy(6)
            real(dp) :: f(6), b(3)
            if (fo) then
                call field(pm, yy(1:3), region, b)
                f = [yy(4:6), cross(yy(4:6), b)]
            else
                f = gc_rhs(pm, yy, region, order)
            end if
        end function
    end subroutine

    subroutine advance(pm, opt, orbit, span)
        type(model_t), intent(in) :: pm
        type(options_t), intent(in) :: opt
        type(orbit_t), intent(inout) :: orbit
        real(dp), intent(in) :: span
        real(dp) :: t, h, sg, err, yn(6), ym(6), g(6), residual, rho
        real(dp) :: a, b, c, fc, fa, fb, hl, e_before, f(6)
        real(dp) :: va, vb, vc, left, right, peak(6), peak_time
        type(local_t) :: lf
        logical :: ok, event, switch_event
        integer :: it, attempts, ns, target
        if (orbit%failed) return
        t = 0.0_dp
        sg = sign(1.0_dp, span)
        attempts = 0
        do while (t < abs(span) - 1.0e-13_dp)
            attempts = attempts + 1
            if (attempts > 100000 .or. orbit%h < 1.0e-13_dp) then
                orbit%failed = .true.
                return
            end if
            call geometry(pm, orbit%y(1:3), orbit%region, lf)
            rho = orbit%speed/lf%bmag
            h = min(orbit%h, abs(span) - t, 0.8_dp/lf%bmag)
            if (.not. orbit%fo .and. opt%method /= 1) then
                f = gc_rhs(pm, orbit%y, orbit%region, opt%order)
                h = min(h, 0.4_dp*opt%din*rho/max(norm2(f(1:3)), 1.0e-12_dp))
            end if
            if (opt%method == 5) then
                h = min(opt%boris_step, abs(span) - t)
                call boris(pm, orbit%y, sg*h, yn)
                orbit%y = yn
                orbit%steps = orbit%steps + 1
                orbit%evaluations = orbit%evaluations + 1
                t = t + h
                orbit%fo_time = orbit%fo_time + h
                cycle
            end if
            if (opt%method == 4 .and. orbit%fo) then
                h = min(opt%patch_step/lf%bmag, abs(span) - t)
                call local_patch(pm, orbit%y, sg*h, yn, ns, ok)
                orbit%crossings = orbit%crossings + ns
                orbit%failed = .not. ok
                if (.not. ok) return
                err = 0.0_dp
                orbit%evaluations = orbit%evaluations + 2
            else
                call dp_step(pm, orbit%y, orbit%region, orbit%fo, opt%order, &
                             sg*h, orbit%speed, yn, err)
                orbit%evaluations = orbit%evaluations + 7
                err = err/opt%tol
                if (.not. all(ieee_is_finite(yn))) then
                    orbit%failed = .true.
                    return
                end if
                if (err > 1.0_dp) then
                    orbit%h = h*max(0.2_dp, 0.9_dp*err**(-0.2_dp))
                    cycle
                end if
            end if
            event = .false.
            switch_event = .false.
            if (opt%method /= 4 .or. .not. orbit%fo) then
                if (orbit%fo .or. opt%method == 1) then
                    fa = merge(-1.0_dp, 1.0_dp, orbit%region == 1)* &
                         max(abs(distance(pm, orbit%y(1:3))), 1.0e-16_dp)
                    fb = distance(pm, yn(1:3))
                    b = h
                    ! An endpoint sign test alone misses a crossing and return
                    ! in one step. The half-step catches resolvable double roots.
                    if (fa*fb >= 0.0_dp .and. &
                        abs(fa) < 1.1_dp*rho) then
                        call dp_step(pm, orbit%y, orbit%region, orbit%fo, &
                            opt%order, sg*h/2, orbit%speed, ym, fc)
                        orbit%evaluations = orbit%evaluations + 7
                        fc = distance(pm, ym(1:3))
                        if (fa*fc < 0.0_dp) then
                            fb = fc
                            b = h/2
                            yn = ym
                        end if
                    end if
                    if (orbit%fo .and. fa*fb >= 0.0_dp .and. &
                        abs(fa) < 1.1_dp*rho) then
                        va = sg*dot_product(normal(pm, orbit%y(1:3)), &
                            orbit%y(4:6))
                        vb = sg*dot_product(normal(pm, yn(1:3)), yn(4:6))
                        if (va*vb < 0.0_dp .and. &
                            va*merge(-1.0_dp, 1.0_dp, orbit%region == 1) < 0) then
                            ! Home-branch radial extremum: a shallow double
                            ! crossing can occur around it even if the midpoint
                            ! is on the home side. h*Omega<0.8 resolves turns.
                            left = 0.0_dp
                            right = b
                            do it = 1, 35
                                peak_time = left + (right - left)*abs(va) &
                                    /(abs(va) + abs(vb))
                                call dp_step(pm, orbit%y, orbit%region, .true., &
                                    0, sg*peak_time, orbit%speed, peak, vc)
                                orbit%evaluations = orbit%evaluations + 7
                                vc = sg*dot_product(normal(pm, peak(1:3)), peak(4:6))
                                if (abs(vc) < 1.0e-11_dp*orbit%speed) exit
                                if (va*vc > 0) then
                                    left = peak_time
                                    va = vc
                                else
                                    right = peak_time
                                    vb = vc
                                end if
                            end do
                            fc = distance(pm, peak(1:3))
                            if (fa*fc < 0) then
                                fb = fc
                                b = peak_time
                                yn = peak
                            end if
                        end if
                    end if
                    event = fa*fb < 0.0_dp
                else
                    fa = boundary(orbit%y)
                    fb = boundary(yn)
                    b = h
                    switch_event = fb <= 0.0_dp
                    event = switch_event
                end if
                if (event) then
                    a = 0.0_dp
                    do it = 1, 60
                        c = 0.5_dp*(a + b)
                        call dp_step(pm, orbit%y, orbit%region, orbit%fo, &
                            opt%order, sg*c, orbit%speed, ym, err)
                        orbit%evaluations = orbit%evaluations + 7
                        if (switch_event) then
                            fc = boundary(ym)
                        else
                            fc = distance(pm, ym(1:3))
                        end if
                        if (fa*fc > 0.0_dp) then
                            a = c
                            fa = fc
                        else
                            b = c
                        end if
                        if (b - a < 1.0e-11_dp) exit
                    end do
                    hl = 0.5_dp*(a + b)
                    if (hl < 1.0e-9_dp*h .and. .not. switch_event) then
                        orbit%region = 3 - orbit%region
                        cycle
                    end if
                    h = hl
                    call dp_step(pm, orbit%y, orbit%region, orbit%fo, &
                        opt%order, sg*h, orbit%speed, yn, err)
                    orbit%evaluations = orbit%evaluations + 7
                    err = err/opt%tol
                    if (.not. switch_event) then
                        orbit%crossings = orbit%crossings + 1
                        if (opt%method == 1) then
                            target = 3 - orbit%region
                            call reset_old(pm, yn, orbit%region, &
                                opt%old_level, ok, target)
                            if (.not. ok) then
                                orbit%failed = .true.
                                return
                            end if
                            if (target == orbit%region) then
                                f = gc_rhs(pm, yn, orbit%region, 0)
                                fc = sg*dot_product(normal(pm, yn(1:3)), f(1:3))
                                if (fc*merge(-1.0_dp, 1.0_dp, target == 1) <= 0) then
                                    ! The mirrored GC still points out of its
                                    ! assigned volume. Its reset-flow contract
                                    ! fails; do not nudge it through artificially.
                                    orbit%failed = .true.
                                    return
                                end if
                            end if
                            orbit%region = target
                        else
                            orbit%region = 3 - orbit%region
                        end if
                    end if
                end if
            end if
            orbit%y = yn
            t = t + h
            orbit%steps = orbit%steps + 1
            if (orbit%fo) orbit%fo_time = orbit%fo_time + h
            if (switch_event) then
                call lift(pm, orbit%y, opt%order, yn, ok)
                if (opt%method == 2) then
                    orbit%y(6) = 0.0_dp
                    if (opt%phase_rule == 1 .and. abs(pm%delta) < tiny(1.0_dp)) &
                        call match_phase(pm, orbit%y, orbit%p0, ok)
                    if (.not. ok) orbit%roots_failed = orbit%roots_failed + 1
                    call lift(pm, orbit%y, 0, yn, ok)
                end if
                e_before = gc_energy(pm, orbit%y, opt%order)
                orbit%switch_energy = max(orbit%switch_energy, &
                    abs(0.5_dp*dot_product(yn(4:6), yn(4:6))/e_before - 1))
                orbit%y = yn
                orbit%fo = .true.
                orbit%region = side_of(pm, yn(1:3))
                orbit%switches = orbit%switches + 1
            else if (orbit%fo .and. opt%method >= 2 .and. opt%method <= 4) then
                ! Cheap rejection before Newton: a GC differs by about rho.
                if (abs(distance(pm, orbit%y(1:3))) > (opt%dout - 1.4_dp)*rho) then
                    call geometry(pm, orbit%y(1:3), orbit%region, lf)
                    g(1:3) = orbit%y(1:3) &
                        + cross(orbit%y(4:6), lf%b)/lf%bmag
                    if (abs(distance(pm, g(1:3))) < (opt%dout - 0.25_dp)*rho) then
                        orbit%h = max(1.0e-8_dp, h*min(5.0_dp, &
                            0.9_dp*max(err, 1.0e-12_dp)**(-0.2_dp)))
                        cycle
                    end if
                    call inverse(pm, orbit%y, opt%order, g, ok, residual)
                    if (ok) then
                        call geometry(pm, g(1:3), side_of(pm, g(1:3)), lf)
                        rho = orbit%speed/lf%bmag
                        if (abs(distance(pm, g(1:3))) > opt%dout*rho) then
                            e_before = 0.5_dp*sum(orbit%y(4:6)**2)
                            if (opt%method == 2) then
                                call discarded_projection(pm, orbit%y, &
                                    orbit%p0, g, ok)
                            end if
                            if (ok) then
                                orbit%switch_energy = max(orbit%switch_energy, &
                                    abs(gc_energy(pm, g, opt%order)/e_before - 1))
                                orbit%y = g
                                orbit%fo = .false.
                                orbit%region = side_of(pm, g(1:3))
                                orbit%switches = orbit%switches + 1
                            else
                                orbit%roots_failed = orbit%roots_failed + 1
                            end if
                        end if
                    else
                        orbit%roots_failed = orbit%roots_failed + 1
                    end if
                end if
            end if
            orbit%h = max(1.0e-8_dp, h*min(5.0_dp, &
                0.9_dp*max(err, 1.0e-12_dp)**(-0.2_dp)))
        end do
        orbit%time = orbit%time + span
    contains
        function boundary(g) result(bd)
            real(dp), intent(in) :: g(6)
            real(dp) :: bd
            type(local_t) :: lg
            call geometry(pm, g(1:3), orbit%region, lg)
            bd = abs(distance(pm, g(1:3))) - opt%din*orbit%speed/lg%bmag
        end function
    end subroutine

    subroutine discarded_projection(pm, y, p0, g, ok)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: y(6), p0
        real(dp), intent(out) :: g(6)
        logical, intent(out) :: ok
        type(local_t) :: lf
        real(dp) :: rr, bphi, dummy(6), pbase, e
        call geometry(pm, y(1:3), side_of(pm, y(1:3)), lf)
        g(1:3) = y(1:3) + cross(y(4:6), lf%b)/lf%bmag
        call geometry(pm, g(1:3), side_of(pm, g(1:3)), lf)
        g(4) = dot_product(y(4:6), lf%b)
        if (abs(pm%delta) < tiny(1.0_dp)) then
            rr = hypot(g(1), g(2))
            bphi = (-g(2)*lf%b(1) + g(1)*lf%b(2))/rr
            dummy = [g(1:3), 0.0_dp, 0.0_dp, 0.0_dp]
            pbase = canonical_pphi(pm, dummy)
            g(4) = (p0 - pbase)/(rr*bphi)
        end if
        e = 0.5_dp*sum(y(4:6)**2)
        g(5) = (e - 0.5_dp*g(4)**2)/lf%bmag
        g(6) = 0.0_dp
        ok = g(5) > 0.0_dp
    end subroutine

    subroutine match_phase(pm, g, ptarget, ok)
        type(model_t), intent(in) :: pm
        real(dp), intent(inout) :: g(6)
        real(dp), intent(in) :: ptarget
        logical, intent(out) :: ok
        real(dp) :: a, b, c, fa, fb, fc, y(6), gt(6), best, score
        logical :: good
        integer :: k, it
        gt = g
        ok = .false.
        best = huge(1.0_dp)
        do k = 0, 31
            a = 2*PI*k/32
            b = 2*PI*(k + 1)/32
            gt(6) = a
            call lift(pm, gt, 0, y, good)
            fa = canonical_pphi(pm, y) - ptarget
            gt(6) = b
            call lift(pm, gt, 0, y, good)
            fb = canonical_pphi(pm, y) - ptarget
            if (fa*fb > 0) cycle
            do it = 1, 38
                c = 0.5_dp*(a + b)
                gt(6) = c
                call lift(pm, gt, 0, y, good)
                fc = canonical_pphi(pm, y) - ptarget
                if (fa*fc > 0) then
                    a = c
                    fa = fc
                else
                    b = c
                end if
            end do
            score = abs(modulo(c + PI, 2*PI) - PI)
            if (score < best) then
                best = score
                g(6) = c
                ok = .true.
            end if
        end do
    end subroutine

    subroutine reset_old(pm, g, home, level, ok, target)
        type(model_t), intent(in) :: pm
        real(dp), intent(inout) :: g(6)
        integer, intent(in) :: home, level
        logical, intent(out) :: ok
        integer, intent(out), optional :: target
        type(local_t) :: lh, lt
        real(dp) :: e, rad, z0(6), trial(6), mid(6), kick(6), lam
        real(dp) :: n(3), bh(3), ch(3), bp, f, df, ltry, ftry, rr, s, angle
        integer :: it, bt, k
        call geometry(pm, g(1:3), home, lh)
        call geometry(pm, g(1:3), 3 - home, lt)
        e = 0.5_dp*g(4)**2 + g(5)*lh%bmag
        rad = 2.0_dp*(e - g(5)*lt%bmag)
        ok = .true.
        if (present(target)) target = 3 - home
        if (rad < 0.0_dp) then
            g(4) = -g(4)
            if (present(target)) target = home
            return
        end if
        z0 = g
        if (level == 1) then
            lam = 0.0_dp
            do it = 1, 20
                call residual(lam, trial, f, df)
                if (abs(f) < 1.0e-10_dp*e) then
                    g = trial
                    return
                end if
                if (abs(df) < 1.0e-20_dp) exit
                ltry = lam - f/df
                do bt = 1, 10
                    call residual(ltry, mid, ftry, rr)
                    if (abs(ftry) < abs(f)) exit
                    ltry = 0.5_dp*(ltry + lam)
                end do
                if (abs(ftry) >= abs(f)) exit
                lam = ltry
                s = hypot(trial(1), trial(2)) - pm%torus%r_major
                angle = atan2(trial(3), s) &
                    - atan2(z0(3), hypot(z0(1), z0(2)) - pm%torus%r_major)
                if (abs(angle) > 0.5_dp) exit
            end do
        end if
        g = z0
        g(4) = sign(sqrt(rad), z0(4))
    contains
        subroutine residual(l, trial, res, deriv)
            real(dp), intent(in) :: l
            real(dp), intent(out) :: trial(6), res, deriv
            real(dp) :: xp(3), nm(3), r, cp, sp
            mid = z0
            kick = 0.0_dp
            do k = 1, 5
                mid = z0 + 0.5_dp*l*kick
                call geometry(pm, mid(1:3), home, lh)
                call geometry(pm, mid(1:3), 3 - home, lt)
                bh = 0.5_dp*(lh%b + lt%b)
                ch = 0.5_dp*(lh%curlb + lt%curlb)
                bp = 0.5_dp*(lh%bmag + lt%bmag) + mid(4)*dot_product(bh, ch)
                n = normal(pm, mid(1:3))
                kick(1:3) = -cross(bh, n)/bp
                kick(4) = mid(4)*dot_product(ch, n)/bp
            end do
            trial = z0 + l*kick
            ! Re-express the tangential kick on the same toroidal surface.
            r = pm%torus%r_sheet
            cp = trial(1)/hypot(trial(1), trial(2))
            sp = trial(2)/hypot(trial(1), trial(2))
            angle = atan2(trial(3), hypot(trial(1), trial(2)) - pm%torus%r_major)
            xp = [(pm%torus%r_major + r*cos(angle))*cp, &
                  (pm%torus%r_major + r*cos(angle))*sp, r*sin(angle)]
            trial(1:3) = xp
            call geometry(pm, xp, 3 - home, lt)
            nm = normal(pm, xp)
            res = 0.5_dp*trial(4)**2 + z0(5)*lt%bmag - e
            deriv = trial(4)*kick(4) + z0(5)* &
                dot_product(lt%gradb, kick(1:3) - nm*dot_product(nm, kick(1:3)))
        end subroutine
    end subroutine

    subroutine helix(x, v, bvec, t, xn, vn)
        real(dp), intent(in) :: x(3), v(3), bvec(3), t
        real(dp), intent(out) :: xn(3), vn(3)
        real(dp) :: om, b(3), vp(3), parallel(3), vc(3), sn, co
        om = norm2(bvec)
        b = bvec/om
        parallel = dot_product(v, b)*b
        vp = v - parallel
        vc = cross(v, b)
        sn = sin(om*t)
        co = cos(om*t)
        vn = parallel + vp*co + vc*sn
        xn = x + parallel*t + vp*sn/om + vc*(2*sin(om*t/2)**2)/om
    end subroutine

    subroutine plane_step(point, n, bm, bp, x, v, t, xn, vn, crossings, ok)
        real(dp), intent(in) :: point(3), n(3), bm(3), bp(3), x(3), v(3), t
        real(dp), intent(out) :: xn(3), vn(3)
        integer, intent(out) :: crossings
        logical, intent(out) :: ok
        real(dp) :: bvec(3), om, d, aa, bb, rad, alpha, val, roots(2), theta
        real(dp) :: remain, sg, dt, xx(3), vv(3), candidate, root
        real(dp) :: anchor_x(3), anchor_v(3), anchor_remain, cycle_time, cycle_dx(3)
        integer :: region, j, k, explicit_events
        integer(int64) :: cycles
        logical :: on_plane
        xx = x
        vv = v
        sg = sign(1.0_dp, t)
        remain = abs(t)
        d = dot_product(n, xx - point)
        on_plane = abs(d) < 1.0e-13_dp
        region = merge(1, 2, d < 0.0_dp)
        if (abs(d) < 1.0e-13_dp) &
            region = merge(2, 1, sg*dot_product(n, vv) > 0.0_dp)
        crossings = 0
        explicit_events = 0
        anchor_remain = 0.0_dp
        anchor_x = x
        anchor_v = v
        ok = .true.
        do while (remain > 1.0e-13_dp)
            bvec = bm
            if (region == 2) bvec = bp
            om = norm2(bvec)
            d = dot_product(n, xx - point)
            if (on_plane) d = 0.0_dp
            aa = sg*dot_product(n, vv)/om
            bb = dot_product(n, cross(vv, bvec/om))/om
            rad = hypot(aa, bb)
            dt = remain
            root = huge(1.0_dp)
            if (on_plane .and. abs(aa) > tiny(1.0_dp)) then
                ! Stable first return from an oriented sheet crossing:
                ! tan(theta/2)=-aa/bb. Avoid asin/atan cancellation when
                ! an excursion is much shorter than one gyroperiod.
                theta = 2*atan2(abs(aa), -sign(1.0_dp, aa)*bb)
                root = theta/om
            else if (rad > 1.0e-16_dp) then
                val = -(d + bb)/rad
                if (abs(val) <= 1.0_dp) then
                    alpha = atan2(-bb, aa)
                    roots = [asin(val) - alpha, PI - asin(val) - alpha]
                    do j = 1, 2
                        theta = modulo(roots(j), 2*PI)
                        do k = 0, 1
                            candidate = (theta + 2*PI*k)/om
                            if (candidate > 1.0e-10_dp/om) root = min(root, candidate)
                        end do
                    end do
                end if
            end if
            dt = min(dt, root)
            call helix(xx, vv, bvec, sg*dt, xn, vn)
            xx = xn
            vv = vn
            remain = remain - dt
            if (root <= dt + 1.0e-13_dp) then
                crossings = crossings + 1
                explicit_events = explicit_events + 1
                on_plane = .true.
                region = 3 - region
                if (explicit_events == 1) then
                    anchor_x = xx
                    anchor_v = vv
                    anchor_remain = remain
                else if (explicit_events == 3) then
                    ! In a uniform tangent sheet, both tangential canonical
                    ! momenta are conserved. An oriented crossing repeats v
                    ! after two excursions and translates x tangentially.
                    cycle_time = anchor_remain - remain
                    if (cycle_time > 1.0e-12_dp .and. &
                        norm2(vv - anchor_v) < 1.0e-10_dp*norm2(vv)) then
                        cycles = floor(remain/cycle_time, kind=int64)
                        if (cycles > 0) then
                            cycle_dx = xx - anchor_x
                            cycle_dx = cycle_dx - n*dot_product(n, cycle_dx)
                            xx = xx + real(cycles, dp)*cycle_dx
                            remain = max(0.0_dp, remain - real(cycles, dp)*cycle_time)
                            if (cycles <= (int(huge(crossings), int64) &
                                - int(crossings, int64))/2_int64) then
                                crossings = crossings + int(2*cycles)
                            else
                                crossings = huge(crossings)
                            end if
                        end if
                    end if
                end if
            end if
            if (explicit_events > 128) then
                ok = .false.
                return
            end if
        end do
        xn = xx
        vn = vv
    end subroutine

    subroutine local_patch(pm, y, t, yn, crossings, ok)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: y(6), t
        real(dp), intent(out) :: yn(6)
        integer, intent(out) :: crossings
        logical, intent(out) :: ok
        real(dp) :: n(3), point(3), bm(3), bp(3)
        n = normal(pm, y(1:3))
        point = y(1:3) - distance(pm, y(1:3))*n
        ! Taylor expansion about the patch start, not a permanent projection
        ! of the field onto the sheet. Both analytically continued branches
        ! remain tangent to n on nested minor-radius surfaces.
        call field(pm, y(1:3), 1, bm)
        call field(pm, y(1:3), 2, bp)
        call plane_step(point, n, bm, bp, y(1:3), y(4:6), t, &
                         yn(1:3), yn(4:6), crossings, ok)
    end subroutine

    subroutine boris(pm, y, h, yn)
        type(model_t), intent(in) :: pm
        real(dp), intent(in) :: y(6), h
        real(dp), intent(out) :: yn(6)
        real(dp) :: x(3), v(3), b(3), t(3), s(3), vp(3)
        x = y(1:3) + 0.5_dp*h*y(4:6)
        call field(pm, x, side_of(pm, x), b)
        t = 0.5_dp*h*b
        s = 2*t/(1 + dot_product(t, t))
        v = y(4:6)
        vp = v + cross(v, t)
        v = v + cross(vp, s)
        yn = [x + 0.5_dp*h*v, v]
    end subroutine

end module
