program check_switch_identities
    ! Item D copy of SIMPLE-gcho 598f041 concept/gc_higher_order/gen/derive_gc.f90,
    ! truncated after the selector-0/1 GC->particle maps (value + 6x6 Jacobian).
    ! Physics unchanged; only the emitted module name and provenance strings differ.
    ! Consumer-owned physics. q/m=1/epsilon; theta(paper)=-phi(here).
    ! TB = Tronko & Brizard, PoP 22, 112507 (2015), arXiv:1511.01005v1.
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, func, sqrt, sin, cos, &
                        operator(+), operator(-), operator(*), operator(/), operator(**)
    use fortsym_products, only: jvp
    use fortsym_subs, only: subs_many
    use fortsym, only: simplify
    use fortsym_engine, only: engine_result_t
    use fortsym_string, only: str, chars
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, TARGET_FORTRAN_CPU
    use fortsym_lie, only: lie_series
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_eval, only: binding_t, eval_expr
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: zeros(3)
    type(expr_t) :: p(130), packed_refs(130), field(3), bh(3), e1(3), e2(3), av(3), dv(3)
    type(expr_t) :: dbh(3, 3), gmag(3), rvec(3), rstar(3), cr(3), cb(3)
    type(expr_t) :: kap(3), rho(3), rh(3), vt(3), vel(3), g1(6), coords(6)
    type(expr_t) :: g2x(3), rho1(3), pi1(3), dir, bn, scale, w, tau, al1, al2
    type(expr_t) :: gm, gu, gt, dotrk, h2, avg_g2mu, beta_perp, beta_par
    type(expr_t) :: curl_r(3), alpha_sq, divb, divk, curl_bg(3), gradlog(3)
    type(expr_t) :: work(3), work2(3), neg_g1(6)
    type(expr_t) :: tmap(6), roots(42), rroots(12), bs(3), bp, xd(3), ud, phid
    type(expr_t) :: vdr(3), trial, hgrad(3), hu, hmu
    type(expr_t) :: p2v(3), p3v(3), p4v(3), d1p2(3), d1p4(3), d11p3(3)
    type(expr_t) :: d11p4(3), sigma3, g2(6), g3x(3), g3par, qtheta, fa(3)
    type(expr_t) :: basis(6, 2), fullroots(126), qroots(13), qa, qb
    type(expr_t) :: grouped(18), tangents(18), qtangents(2)
    integer :: a, m
    character(256) :: filename
    character(:), allocatable :: why
    integer :: i, j, k, l, n, unit
    logical :: ok
    type(symengine_engine_t) :: eng
    type(suite_t) :: suite
    type(expr_t) :: wsym, tsym, csym, ksym, gsym, zc(6), bx_field(3), tm(6)
    type(expr_t) :: jac(6, 6), xpart(3), bpart(3), wgc(6, 6), epull, edef, coef
    character(256) :: label

    call arena%init()
    do i = 1, size(p)
        p(i) = sym(arena, 'p'//idx(i))
        packed_refs(i) = func('p', [num(arena, i)])
    end do
    coords = p(1:6)
    zeros = num(arena, 0)
    field = p(8:10)
    bn = sqrt(dot(field, field))
    do i = 1, 3
        bh(i) = field(i)/bn
    end do
    scale = dot(p(47:49), bh)
    do i = 1, 3
        e1(i) = p(46 + i) - scale*bh(i)
    end do
    scale = sqrt(dot(e1, e1))
    do i = 1, 3
        e1(i) = e1(i)/scale
    end do
    e2 = cross(bh, e1)
    do i = 1, 3
        av(i) = cos(p(6))*e1(i) + sin(p(6))*e2(i)
        dv(i) = -sin(p(6))*e1(i) + cos(p(6))*e2(i)
    end do
    wsym = sym(arena, 'wperp')
    w = wsym
    do i = 1, 3
        vt(i) = w*av(i)
        vel(i) = p(4)*bh(i) + vt(i)
    end do
    rho = cross(bh, vt)
    do i = 1, 3
        rho(i) = rho(i)/bn
    end do
    rh = cross(bh, av)
    do j = 1, 3
        gmag(j) = spatial(bn, j)
        rvec(j) = num(arena, 0)
        do i = 1, 3
            dbh(i, j) = spatial(bh(i), j)
            rvec(j) = rvec(j) + e2(i)*spatial(e1(i), j)
        end do
    end do
    kap = matvec(dbh, bh)
    cb = curl(bh)
    tau = dot(bh, cb)
    do i = 1, 3
        gradlog(i) = gmag(i)/bn
    end do
    al1 = num(arena, 0)
    al2 = num(arena, 0)
    do j = 1, 3
        do i = 1, 3
            al1 = al1 - (rh(i)*av(j) + av(i)*rh(j))*dbh(i, j)/2
            al2 = al2 + (av(i)*av(j) - rh(i)*rh(j))*dbh(i, j)/4
        end do
    end do
    ! TB (41)-(44), with theta=-phi. These are not G2 components.
    dotrk = dot(rho, kap)
    gu = -p(4)*dotrk + p(5)*(tau + al1)
    do i = 1, 3
        work(i) = p(5)*gradlog(i) + p(4)**2*kap(i)/bn
        work2(i) = gradlog(i) + p(4)**2*kap(i)/(2*p(5)*bn)
    end do
    gm = dot(rho, work) - p(5)*p(4)*(tau + al1)/bn
    gt = dot(vt, work2)/bn &
         - dot(rho, rvec) + p(4)*al2/bn
    do i = 1, 3
        g1(i) = -rho(i)
    end do
    g1(4) = gu
    g1(5) = gm
    g1(6) = -gt
    do i = 1, 3
        pi1(i) = -p(5)*cb(i)/2
    end do
    ! TB (45),(55). Spatial correction is needed one order beyond velocity.
    work = cross(pi1, bh)
    do i = 1, 3
        g2x(i) = (2*p(4)*dot(vt, kap)/bn**2 + p(5)*al2/bn)*bh(i) &
                 - work(i)/bn &
                 + (p(4)**2*dotrk/bn + p(5)*p(4)*(3*tau - al1)/bn)*rho(i)/(4*p(5)) &
                 + (p(4)*al2/bn + dot(vt, work2)/bn)*vt(i)/(2*bn)
    end do
    do i = 1, 3
        dir = num(arena, 0)
        do j = 1, 6
            dir = dir + g1(j)*derivative(rho(i), j)
        end do
        rho1(i) = -g2x(i) - dir/2
    end do
    ! ---------------- item D exact pullback identities ----------------
    ! The construction above is the prototype's (derive_gc.f90) with the single
    ! change w = wsym (perpendicular speed as independent coordinate), so that
    ! mu = wsym**2/(2*bn) is substituted below and no sqrt of a symbol remains.
    eng = make_symengine_engine(arena)
    call suite_begin(suite, 'item D first-order switch exact identities')
    call cotangent_lift_checks()
    call lie_truncation_checks()
    tsym = sym(arena, 'tq')
    csym = sym(arena, 'cfield')
    ksym = sym(arena, 'kfield')
    gsym = sym(arena, 'gfield')
    zc = [p(1), p(2), p(3), p(4), wsym, p(6)]
    ! Field B = (-c X3, k X1, 1 + g X1): divergence free, grad-B, curvature,
    ! shear and a varying gauge frame (axis = x).
    bx_field = [-csym*p(3), ksym*p(1), 1 + gsym*p(1)]
    do n = 0, 1
        do i = 1, 3
            tmap(i) = p(i) + p(7)*rho(i)
        end do
        tmap(4:6) = vel
        if (n == 1) then
            do i = 1, 3
                tmap(i) = tmap(i) + p(7)**2*rho1(i)
                call lie_series(vel(i), coords, neg_g1, p(7), 1, trial, ok, why)
                if (.not. ok) error stop 'Lie pullback refused'
                tmap(3 + i) = trial
            end do
        end if
        call explicit_field(tmap, tm)
        do a = 1, 6
            do j = 1, 6
                jac(a, j) = at_point(diff(tm(a), zc(j)))
            end do
        end do
        do i = 1, 3
            xpart(i) = at_point(tm(i))
        end do
        bpart = [-csym*xpart(3), ksym*xpart(1), 1 + gsym*xpart(1)]
        call gc_form(n, wgc)
        do i = 1, 6
            do j = i + 1, 6
                epull = num(arena, 0)
                do a = 1, 3
                    do m = 1, 3
                        epull = epull + jac(a, i)*levi_b(a, m, bpart)*jac(m, j)
                    end do
                    epull = epull + p(7)*(jac(3 + a, i)*jac(a, j) - jac(a, i)*jac(3 + a, j))
                end do
                edef = epull - wgc(i, j)
                do k = 0, n + 1
                    coef = subs_many(edef, [p(7)], [num(arena, 0)])
                    write (label, '(a,i0,a,i0,a,i0,a,i0,a)') 'selector ', n, &
                        ': eps^', k - 1, ' coefficient of (T^*omega_p - omega_GC)(', &
                        i, ',', j, ') vanishes (field c,k,g symbolic)'
                    call check_identity(suite, eng, trim(label), clear_t(coef))
                    edef = diff(edef, p(7))
                end do
            end do
        end do
    end do
    call suite_end(suite, 'results/symbolic_identities.json')

contains
    subroutine explicit_field(tin, tout)
        type(expr_t), intent(in) :: tin(6)
        type(expr_t), intent(out) :: tout(6)
        type(expr_t) :: old(123), new(123), bnf
        integer :: i, j, q
        bnf = sqrt(dot(field, field))
        do i = 1, 6
            tout(i) = subs_many(tin(i), [p(5)], [wsym**2/(2*bnf)])
        end do
        q = 0
        do i = 1, 3
            q = q + 1
            old(q) = p(7 + i)
            new(q) = bx_field(i)
        end do
        do j = 1, 3
            do i = 1, 3
                q = q + 1
                old(q) = p(10 + i + 3*(j - 1))
                new(q) = diff(bx_field(i), p(j))
            end do
        end do
        do i = 20, 46
            q = q + 1
            old(q) = p(i)
            new(q) = num(arena, 0)
        end do
        do i = 1, 3
            q = q + 1
            old(q) = p(46 + i)
            new(q) = num(arena, merge(1, 0, i == 1))
        end do
        do i = 50, 130
            q = q + 1
            old(q) = p(i)
            new(q) = num(arena, 0)
        end do
        do i = 1, 6
            tout(i) = subs_many(tout(i), old(1:q), new(1:q))
        end do
    end subroutine explicit_field

    function at_point(e) result(r)
        type(expr_t), intent(in) :: e
        type(expr_t) :: r
        type(engine_result_t) :: s
        r = subs_many(e, [p(1), p(2), p(3), cos(p(6)), sin(p(6))], &
            [num(arena, 0), num(arena, 0), num(arena, 0), &
             (1 - tsym**2)/(1 + tsym**2), 2*tsym/(1 + tsym**2)])
        s = simplify(r)
        if (s%ok) r = s%value
    end function at_point

    function levi_b(a, m, bv) result(r)
        integer, intent(in) :: a, m
        type(expr_t), intent(in) :: bv(3)
        type(expr_t) :: r
        integer :: k
        r = num(arena, 0)
        do k = 1, 3
            if (k == a .or. k == m .or. a == m) cycle
            if (modulo(m - a, 3) == 1) then
                r = r + bv(k)
            else
                r = r - bv(k)
            end if
        end do
    end function levi_b

    !> eps * omega_GC in coordinates (X, u, w, phi), at X = 0.
    !> Gamma_GC = (A/eps + u b - eps mu Rstar).dX - eps mu dphi, mu = w^2/(2|B|),
    !> Rstar = R + curl(b)/2 (selector 1), R (selector 0), R_j = e2.d_j e1.
    subroutine gc_form(n, wout)
        integer, intent(in) :: n
        type(expr_t), intent(out) :: wout(6, 6)
        type(expr_t) :: bnf, bb(3), e1v(3), e2v(3), rv(3), cbv(3), mux, gam(6), sc
        integer :: i, j
        bnf = sqrt(bx_field(1)**2 + bx_field(2)**2 + bx_field(3)**2)
        do i = 1, 3
            bb(i) = bx_field(i)/bnf
        end do
        sc = bb(1)
        do i = 1, 3
            e1v(i) = merge(1, 0, i == 1) - sc*bb(i)
        end do
        sc = sqrt(e1v(1)**2 + e1v(2)**2 + e1v(3)**2)
        do i = 1, 3
            e1v(i) = e1v(i)/sc
        end do
        e2v = cross(bb, e1v)
        do j = 1, 3
            rv(j) = num(arena, 0)
            do i = 1, 3
                rv(j) = rv(j) + e2v(i)*diff(e1v(i), p(j))
            end do
        end do
        cbv = [diff(bb(3), p(2)) - diff(bb(2), p(3)), diff(bb(1), p(3)) - diff(bb(3), p(1)), &
               diff(bb(2), p(1)) - diff(bb(1), p(2))]
        mux = wsym**2/(2*bnf)
        do i = 1, 3
            gam(i) = p(4)*bb(i) - p(7)*mux*rv(i)
            if (n >= 1) gam(i) = gam(i) - p(7)*mux*cbv(i)/2
        end do
        gam(4) = num(arena, 0)
        gam(5) = num(arena, 0)
        gam(6) = -p(7)*mux
        do i = 1, 6
            do j = 1, 6
                wout(i, j) = p(7)*(diff(gam(j), zc(i)) - diff(gam(i), zc(j)))
                if (i <= 3 .and. j <= 3) wout(i, j) = wout(i, j) + levi_b(i, j, bx_field)
                wout(i, j) = at_point(wout(i, j))
            end do
        end do
    end subroutine gc_form

    !> Clear the (1+t^2) denominators of the rational phase substitution.
    function clear_t(e) result(r)
        type(expr_t), intent(in) :: e
        type(expr_t) :: r
        r = e*(1 + tsym**2)**8
    end function clear_t

    !> FO chart: p_j = (v + A/eps).d_j x is an exact cotangent lift.
    subroutine cotangent_lift_checks()
        type(expr_t) :: gm(3, 3), pc(3), av(3), adj(3, 3), det, vv(3), lhs
        character(8) :: nm
        integer :: i, j
        do i = 1, 3
            do j = 1, 3
                write (nm, '(a,i0,i0)') 'G', i, j
                gm(i, j) = sym(arena, trim(nm))
            end do
            write (nm, '(a,i0)') 'Pc', i
            pc(i) = sym(arena, trim(nm))
            write (nm, '(a,i0)') 'Acart', i
            av(i) = sym(arena, trim(nm))
        end do
        ! adj(G)(i,j) = cofactor(j,i); G^{-T} = adj(G)^T/det
        do i = 1, 3
            do j = 1, 3
                adj(i, j) = gm(mod(j, 3) + 1, mod(i, 3) + 1)*gm(mod(j + 1, 3) + 1, mod(i + 1, 3) + 1) &
                    - gm(mod(j, 3) + 1, mod(i + 1, 3) + 1)*gm(mod(j + 1, 3) + 1, mod(i, 3) + 1)
            end do
        end do
        det = gm(1, 1)*adj(1, 1) + gm(1, 2)*adj(2, 1) + gm(1, 3)*adj(3, 1)
        do i = 1, 3
            vv(i) = (adj(1, i)*pc(1) + adj(2, i)*pc(2) + adj(3, i)*pc(3))/det - av(i)/p(7)
        end do
        do j = 1, 3
            lhs = num(arena, 0)
            do i = 1, 3
                lhs = lhs + (vv(i) + av(i)/p(7))*gm(i, j)
            end do
            write (label, '(a,i0,a)') 'FO chart cotangent lift: (v+A/eps).d_', j, &
                'x = p_j (generic 3x3 chart Jacobian)'
            call check_identity(suite, eng, trim(label), (lhs - pc(j))*det)
        end do
    end subroutine cotangent_lift_checks

    !> Generic Lie-series truncation: positions of exp(eps L_g) truncated at
    !> order N+1, velocities at N; the pulled-back 2-form of (A(x)/eps+v).dx
    !> equals d of the Cartan series sum_k eps^k L_g^k Gamma/k! through eps^N.
    subroutine lie_truncation_checks()
        type(expr_t) :: z(4), g(4), xm(4), gam0(4), lg(4), ser(4), cur(4), pb(4)
        type(expr_t) :: ax(2), d2, e
        character(16) :: nm
        integer :: nn, i, j, k, kk, mx, mv
        real :: dummy
        do i = 1, 4
            write (nm, '(a,i0)') 'zz', i
            z(i) = sym(arena, trim(nm))
        end do
        do i = 1, 4
            write (nm, '(a,i0)') 'ga', i
            g(i) = sym(arena, trim(nm))
            do j = 1, 4
                write (nm, '(a,i0,i0)') 'gb', i, j
                g(i) = g(i) + sym(arena, trim(nm))*z(j)
            end do
            write (nm, '(a,i0)') 'gq', i
            g(i) = g(i) + sym(arena, trim(nm))*z(1)*z(mod(i, 4) + 1)
        end do
        ! generic quadratic potential A(x) in the two position variables
        ax(1) = sym(arena, 'aa11')*z(1) + sym(arena, 'aa12')*z(2) + &
            sym(arena, 'aq1')*z(1)*z(2) + sym(arena, 'ar1')*z(1)**2
        ax(2) = sym(arena, 'aa21')*z(1) + sym(arena, 'aa22')*z(2) + &
            sym(arena, 'aq2')*z(1)*z(2) + sym(arena, 'ar2')*z(2)**2
        dummy = 0.0
        do nn = 0, 2
            mx = nn + 1
            mv = nn
            do i = 1, 4
                k = merge(mx, mv, i <= 2)
                call lie_series(z(i), z, g, p(7), k, xm(i), ok, why)
                if (.not. ok) error stop 'lie_series refused'
            end do
            ! pulled-back one-form, scaled by eps
            do j = 1, 4
                pb(j) = num(arena, 0)
                do i = 1, 2
                    pb(j) = pb(j) + (subs_many(ax(i), z, xm) + p(7)*xm(2 + i))*diff(xm(i), z(j))
                end do
            end do
            ! eps-scaled Cartan series through order nn+2
            gam0 = [ax(1) + p(7)*z(3), ax(2) + p(7)*z(4), num(arena, 0), num(arena, 0)]
            ser = gam0
            cur = gam0
            do kk = 1, nn + 2
                do j = 1, 4
                    lg(j) = num(arena, 0)
                    do i = 1, 4
                        lg(j) = lg(j) + g(i)*diff(cur(j), z(i)) + cur(i)*diff(g(i), z(j))
                    end do
                    lg(j) = lg(j)/kk
                end do
                cur = lg
                do j = 1, 4
                    ser(j) = ser(j) + p(7)**kk*cur(j)
                end do
            end do
            do i = 1, 4
                do j = i + 1, 4
                    d2 = diff(pb(j), z(i)) - diff(pb(i), z(j)) &
                        - (diff(ser(j), z(i)) - diff(ser(i), z(j)))
                    do k = 0, nn + 1
                        e = subs_many(d2, [p(7)], [num(arena, 0)])
                        write (label, '(a,i0,a,i0,a,i0,a,i0,a)') 'Lie truncation N=', nn, &
                            ': eps^', k - 1, ' coefficient of d(T_N^*Gamma) - d(Cartan series) (', &
                            i, ',', j, ') vanishes (generic generator, quadratic A)'
                        call check_identity(suite, eng, trim(label), e)
                        d2 = diff(d2, p(7))
                    end do
                end do
            end do
        end do
    end subroutine lie_truncation_checks

    function spatial_many(f, j) result(value)
        type(expr_t), intent(in) :: f(:)
        integer, intent(in) :: j
        type(expr_t) :: value(size(f))
        integer :: a, k, l, slot
        type(expr_t) :: variables(40), directions(40)
        variables(1) = p(j)
        directions(1) = num(arena, 1)
        slot = 1
        do a = 1, 3
            slot = slot + 1
            variables(slot) = p(7 + a)
            directions(slot) = p(10 + a + 3*(j - 1))
            do k = 1, 3
                slot = slot + 1
                variables(slot) = p(10 + a + 3*(k - 1))
                directions(slot) = p(19 + a + 3*(k - 1) + 9*(j - 1))
                do l = 1, 3
                    slot = slot + 1
                    variables(slot) = p(19 + a + 3*(k - 1) + 9*(l - 1))
                    directions(slot) = p(49 + a + 3*(k - 1) &
                        + 9*(l - 1) + 27*(j - 1))
                end do
            end do
        end do
        value = jvp(f, variables, directions)
    end function spatial_many

    function spatial(f, j) result(value)
        type(expr_t), intent(in) :: f
        integer, intent(in) :: j
        type(expr_t) :: value, output(1)
        output = spatial_many([f], j)
        value = output(1)
    end function spatial

    function derivative_many(f, j) result(value)
        type(expr_t), intent(in) :: f(:)
        integer, intent(in) :: j
        type(expr_t) :: value(size(f)), variables(1), directions(1)
        if (j <= 3) then
            value = spatial_many(f, j)
        else
            variables(1) = p(j)
            directions(1) = num(arena, 1)
            value = jvp(f, variables, directions)
        end if
    end function derivative_many

    function along(f) result(value)
        type(expr_t), intent(in) :: f
        type(expr_t) :: value
        integer :: j
        value = num(arena, 0)
        do j = 1, 6
            value = value + g1(j)*derivative(f, j)
        end do
    end function along

    function paper_diff(f, a) result(value)
        type(expr_t), intent(in) :: f
        integer, intent(in) :: a
        type(expr_t) :: value
        value = diff(f, p(3 + a))
        if (a == 3) value = -value
    end function paper_diff

    function d1(v) result(value)
        type(expr_t), intent(in) :: v(3)
        type(expr_t) :: value(3), cv(3)
        integer :: i, a
        cv = cross(rho, curl(v))
        do i = 1, 3
            value(i) = cv(i)
            do a = 4, 6
                value(i) = value(i) + g1(a)*diff(v(i), p(a))
            end do
        end do
    end function d1

    function diff(f, x) result(value)
        type(expr_t), intent(in) :: f, x
        type(expr_t) :: value, root(1), variables(1), directions(1), output(1)
        root(1) = f
        variables(1) = x
        directions(1) = num(arena, 1)
        output = jvp(root, variables, directions)
        value = output(1)
    end function diff

    function derivative(f, j) result(value)
        type(expr_t), intent(in) :: f
        integer, intent(in) :: j
        type(expr_t) :: value
        if (j <= 3) then
            value = spatial(f, j)
        else
            value = diff(f, p(j))
        end if
    end function derivative

    function idx(i) result(s)
        integer, intent(in) :: i
        character(:), allocatable :: s
        character(16) :: buffer
        write (buffer, '(i0)') i
        s = trim(buffer)
    end function idx

    function dot(a, b) result(c)
        type(expr_t), intent(in) :: a(:), b(:)
        type(expr_t) :: c
        integer :: i
        c = num(arena, 0)
        do i = 1, size(a)
            c = c + a(i)*b(i)
        end do
    end function dot

    function cross(a, b) result(c)
        type(expr_t), intent(in) :: a(3), b(3)
        type(expr_t) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
             a(1)*b(2) - a(2)*b(1)]
    end function cross

    function matvec(a, b) result(c)
        type(expr_t), intent(in) :: a(3, 3), b(3)
        type(expr_t) :: c(3)
        integer :: i
        do i = 1, 3
            c(i) = dot(a(i, :), b)
        end do
    end function matvec

    function curl(a) result(c)
        type(expr_t), intent(in) :: a(3)
        type(expr_t) :: c(3)
        c = [spatial(a(3), 2) - spatial(a(2), 3), &
             spatial(a(1), 3) - spatial(a(3), 1), &
             spatial(a(2), 1) - spatial(a(1), 2)]
    end function curl

    subroutine emit(name, expressions)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: expressions(:)
        type(kernel_spec_t) :: spec
        type(expr_t) :: leaves(size(expressions))
        type(engine_result_t) :: simplified
        integer :: i
        do i = 1, size(expressions)
            leaves(i) = subs_many(expressions(i), p(1:3), zeros)
            simplified = simplify(leaves(i))
            if (.not. simplified%ok) error stop 'native simplification refused'
            leaves(i) = subs_many(simplified%value, p, packed_refs)
        end do
        spec%name = str(name//'_leaf')
        spec%generator = str('concept/fo_gc_derivations/D_canonical_switch/gen/gen_switch1.f90')
        spec%regenerate_command = str('bash concept/fo_gc_derivations/D_canonical_switch/scripts/regen.sh')
        spec%target = TARGET_FORTRAN_CPU
        spec%pure_procedure = .true.
        spec%array_temporaries = .true.
        allocate (spec%args(1), spec%arg_shapes(1), spec%outputs(1))
        allocate (spec%output_shapes(1), spec%output_references(size(expressions)))
        spec%args(1) = str('p')
        spec%arg_shapes(1) = str('(130)')
        spec%outputs(1) = str('result')
        spec%output_shapes(1) = str('('//idx(size(expressions))//')')
        do i = 1, size(expressions)
            spec%output_references(i) = str('result('//idx(i)//')')
        end do
        write (unit, '(a)') chars(emit_kernel(leaves, spec))
        write (unit, '(a)') '    pure subroutine '//name//'(p,result)'
        write (unit, '(a)') '        real(dp), intent(in) :: p(130)'
        write (unit, '(a)') '        real(dp), intent(out) :: result('// &
            idx(size(expressions))//')'
        write (unit, '(a)') '        call '//name//'_leaf(p,result)'
        write (unit, '(a)') '    end subroutine '//name
        flush(unit)
        print '(a)', 'Generated '//name
    end subroutine emit
end program check_switch_identities
