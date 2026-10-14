module candidate_orbits
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use candidate_field, only: model_t, local_t, geometry, field, cross, &
        distance, side_of, normal, canonical_pphi, connection
    use transforms, only: lift, inverse, gc_energy, order2_rhs, kernel_available
    implicit none
    private
    public :: options_t, orbit_t, initialize, advance, particle, helix, &
        plane_step, gc_rhs, dp_step
    real(dp), parameter :: PI = acos(-1.0_dp)

    type :: options_t
        ! 0=sharp-field reference, 3=carried-phase hybrid concept.
        integer :: method = 0, order = 0
        real(dp) :: tol = 1.0e-8_dp, din = 3.0_dp, dout = 4.0_dp
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
        orbit = orbit_t()
        orbit%y = y
        orbit%speed = norm2(y(4:6))
        orbit%energy0 = 0.5_dp*orbit%speed**2
        orbit%p0 = canonical_pphi(pm, y)
        orbit%region = side_of(pm, y(1:3))
        if (opt%method == 3) then
            if (.not. kernel_available(opt%order)) then
                orbit%failed = .true.
                return
            end if
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
        integer :: it, attempts
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
            if (.not. orbit%fo) then
                f = gc_rhs(pm, orbit%y, orbit%region, opt%order)
                h = min(h, 0.4_dp*opt%din*rho/max(norm2(f(1:3)), 1.0e-12_dp))
            end if
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
            event = .false.
            switch_event = .false.
            if (orbit%fo) then
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
                    orbit%region = 3 - orbit%region
                end if
            end if
            orbit%y = yn
            t = t + h
            orbit%steps = orbit%steps + 1
            if (orbit%fo) orbit%fo_time = orbit%fo_time + h
            if (switch_event) then
                call lift(pm, orbit%y, opt%order, yn, ok)
                e_before = gc_energy(pm, orbit%y, opt%order)
                orbit%switch_energy = max(orbit%switch_energy, &
                    abs(0.5_dp*dot_product(yn(4:6), yn(4:6))/e_before - 1))
                orbit%y = yn
                orbit%fo = .true.
                orbit%region = side_of(pm, yn(1:3))
                orbit%switches = orbit%switches + 1
            else if (orbit%fo .and. opt%method == 3) then
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
                            orbit%switch_energy = max(orbit%switch_energy, &
                                abs(gc_energy(pm, g, opt%order)/e_before - 1))
                            orbit%y = g
                            orbit%fo = .false.
                            orbit%region = side_of(pm, g(1:3))
                            orbit%switches = orbit%switches + 1
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

end module
