program derive_gc
    ! Consumer-owned physics. q/m=1/epsilon; theta(paper)=-phi(here).
    ! TB = Tronko & Brizard, PoP 22, 112507 (2015), arXiv:1511.01005v1.
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, sqrt, sin, cos, &
                        operator(+), operator(-), operator(*), operator(/), operator(**)
    use fortsym_products, only: jvp
    use fortsym_subs, only: subs_many
    use fortsym, only: simplify
    use fortsym_engine, only: engine_result_t
    use fortsym_string, only: str, chars
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, TARGET_FORTRAN_CPU
    use fortsym_lie, only: lie_series
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: zeros(3)
    type(expr_t) :: p(49), field(3), bh(3), e1(3), e2(3), av(3), dv(3)
    type(expr_t) :: dbh(3, 3), gmag(3), rvec(3), rstar(3), cr(3), cb(3)
    type(expr_t) :: kap(3), rho(3), rh(3), vt(3), vel(3), g1(6), coords(6)
    type(expr_t) :: g2x(3), rho1(3), pi1(3), dir, bn, scale, w, tau, al1, al2
    type(expr_t) :: gm, gu, gt, dotrk, h2, avg_g2mu, beta_perp, beta_par
    type(expr_t) :: curl_r(3), alpha_sq, divb, divk, curl_bg(3), gradlog(3)
    type(expr_t) :: work(3), work2(3), neg_g1(6)
    type(expr_t) :: tmap(6), roots(42), rroots(12), bs(3), bp, xd(3), ud, phid
    type(expr_t) :: vdr(3), trial
    character(256) :: filename
    character(:), allocatable :: why
    integer :: i, j, k, l, n, unit
    logical :: ok

    call arena%init()
    do i = 1, size(p)
        p(i) = sym(arena, 'p'//idx(i))
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
    w = sqrt(2*p(5)*bn)
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
    call get_command_argument(1, filename)
    if (len_trim(filename) == 0) error stop 'output filename required'
    open (newunit=unit, file=trim(filename), status='replace')
    write (unit, '(a)') 'module gc_generated'
    write (unit, '(a)') '    use, intrinsic :: iso_fortran_env, only: dp => real64'
    write (unit, '(a)') '    implicit none'
    write (unit, '(a)') 'contains'
    ! The Lie engine sees local B coordinates rather than implicit B(X).
    ! This is the chain-rule image of the physical generator, not new algebra.
    coords(1:3) = p(8:10)
    do i = 1, 3
        neg_g1(i) = num(arena, 0)
        do j = 1, 3
            neg_g1(i) = neg_g1(i) + rho(j)*p(10 + i + 3*(j - 1))
        end do
    end do
    do i = 4, 6
        neg_g1(i) = -g1(i)
    end do
    do n = 0, 1
        do i = 1, 3
            tmap(i) = p(i) + p(7)*rho(i)
        end do
        tmap(4:6) = vel
        if (n == 1) then
            do i = 1, 3
                tmap(i) = tmap(i) + p(7)**2*rho1(i)
            end do
            do i = 1, 3
                ! Generic Lie engine supplies the inverse first-order pullback.
                call lie_series(vel(i), coords, neg_g1, p(7), 1, trial, ok, why)
                if (.not. ok) error stop 'Lie pullback refused'
                tmap(3 + i) = trial
            end do
        end if
        roots(1:6) = tmap
        l = 6
        do j = 1, 6
            do i = 1, 6
                l = l + 1
                roots(l) = derivative(tmap(i), j)
            end do
        end do
        call emit('map'//idx(n), roots)
    end do
    ! Reduced Lagrangian: (A/eps+u*b-eps*mu*Rstar).dX-eps*mu*dphi-Hdt.
    ! This truncates H=u^2/2+mu*B. It is not the full second-order theory.
    do n = 0, 1
        if (n == 0) then
            bs = field
            bp = bn
            do i = 1, 3
                xd(i) = p(4)*bh(i)
            end do
            ud = -p(5)*dot(bh, gmag)
            phid = -bn/p(7) - dot(rvec, xd)
        else
            do i = 1, 3
                rstar(i) = rvec(i) + cb(i)/2
            end do
            cr = curl(rstar)
            do i = 1, 3
                bs(i) = field(i) + p(7)*p(4)*cb(i) - p(7)**2*p(5)*cr(i)
            end do
            bp = dot(bh, bs)
            work = cross(bh, gmag)
            do i = 1, 3
                xd(i) = (p(4)*bs(i) + p(7)*p(5)*work(i))/bp
            end do
            ud = -p(5)*dot(bs, gmag)/bp
            phid = -bn/p(7) - dot(rstar, xd)
        end if
        do i = 1, 3
            vdr(i) = xd(i) - p(4)*bh(i)
        end do
        rroots = [xd, ud, phid, bs, bp, vdr]
        call emit('rhs'//idx(n), rroots)
    end do
    ! Diagnostic TB (33)-(38), Pi1=-mu curl(b)/2, Pi2_parallel=0.
    ! A coefficient transcription plus consistency test, not a derived T2.
    curl_r = curl(rvec)
    divb = num(arena, 0)
    divk = num(arena, 0)
    do i = 1, 3
        divb = divb + dbh(i, i)
        divk = divk + spatial(kap(i), i)
    end do
    alpha_sq = dot(bh, curl_r)/2 + (tau**2 + divb**2)/8
    curl_bg = curl(cross(bh, gradlog))
    beta_perp = -tau**2/2 - dot(bh, curl_r) + alpha_sq &
                - dot(cross(bh, gradlog), cross(bh, gradlog)) + dot(bh, curl_bg)/2
    do i = 1, 3
        work(i) = gradlog(i) - kap(i)
        work2(i) = p(5)*gmag(i) + p(4)**2*kap(i)
    end do
    beta_par = -2*alpha_sq - 3*dot(kap, work) + divk
    vdr = cross(bh, work2)
    do i = 1, 3
        vdr(i) = vdr(i)/bn
    end do
    h2 = p(5)**2*beta_perp/2 + p(5)*p(4)**2*beta_par/(2*bn) &
         - p(4)**4*dot(kap, kap)/(2*bn**2) + dot(pi1, vdr)
    do i = 1, 3
        work(i) = 3*kap(i) - gradlog(i)
    end do
    avg_g2mu = p(5)**2/(2*bn)*(tau**2/2 + dot(bh, curl_r) - alpha_sq &
                               - dot(bh, curl_bg)/2) - p(5)*p(4)**2/(2*bn**2) &
               *(dot(kap, work) + divk - tau**2)
    call emit('diagnostic', [rho, gu, gm, gt, rho1, h2, avg_g2mu, tau, al1, al2])
    write (unit, '(a)') 'end module gc_generated'
    close (unit)
contains
    function spatial(f, j) result(value)
        ! Cartesian total derivative of a scalar built from a local field jet.
        ! The consumer owns the B(X) dependency; fortsym supplies all partials.
        type(expr_t), intent(in) :: f
        integer, intent(in) :: j
        type(expr_t) :: value
        integer :: a, k, slot
        type(expr_t) :: variables(13), directions(13), root(1), output(1)
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
            end do
        end do
        root(1) = f
        output = jvp(root, variables, directions)
        value = output(1)
    end function spatial

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

    recursive subroutine mark_inputs(f, used, visited)
        type(expr_t), intent(in) :: f
        logical, intent(inout) :: used(49), visited(:)
        integer :: k
        if (visited(f%id)) return
        visited(f%id) = .true.
        do k = 1, 49
            if (f%id == p(k)%id) used(k) = .true.
        end do
        do k = 1, f%nargs()
            call mark_inputs(f%arg(k), used, visited)
        end do
    end subroutine mark_inputs

    subroutine emit(name, expressions)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: expressions(:)
        type(kernel_spec_t) :: spec
        type(expr_t) :: leaves(size(expressions))
        type(engine_result_t) :: simplified
        logical :: used(49)
        logical, allocatable :: visited(:)
        integer :: slot
        integer :: i
        do i = 1, size(expressions)
            leaves(i) = subs_many(expressions(i), p(1:3), zeros)
            simplified = simplify(leaves(i))
            if (.not. simplified%ok) error stop 'native simplification refused'
            leaves(i) = simplified%value
        end do
        spec%name = str(name//'_leaf')
        spec%generator = str('concept/gc_higher_order/gen/derive_gc.f90')
        spec%regenerate_command = str('make -f concept/gc_higher_order/Makefile regen')
        spec%target = TARGET_FORTRAN_CPU
        spec%pure_procedure = .true.
        allocate (visited(arena%size()), source=.false.)
        used = .false.
        do i = 1, size(leaves)
            call mark_inputs(leaves(i), used, visited)
        end do
        allocate (spec%args(count(used)), spec%outputs(1), spec%output_shapes(1))
        allocate (spec%output_references(size(expressions)))
        slot = 0
        do i = 1, 49
            if (.not. used(i)) cycle
            slot = slot + 1
            spec%args(slot) = str('p'//idx(i))
        end do
        spec%outputs(1) = str('result')
        spec%output_shapes(1) = str('('//idx(size(expressions))//')')
        do i = 1, size(expressions)
            spec%output_references(i) = str('result('//idx(i)//')')
        end do
        write (unit, '(a)') chars(emit_kernel(leaves, spec))
        write (unit, '(a)') '    pure subroutine '//name//'(p,result)'
        write (unit, '(a)') '        real(dp), intent(in) :: p(49)'
        write (unit, '(a)') '        real(dp), intent(out) :: result('// &
            idx(size(expressions))//')'
        write (unit, '(a)') '        call '//name//'_leaf( &'
        do i = 1, 49
            if (.not. used(i)) cycle
            write (unit, '(a)') '            p('//idx(i)//'), &'
        end do
        write (unit, '(a)') '            result)'
        write (unit, '(a)') '    end subroutine '//name
    end subroutine emit
end program derive_gc
