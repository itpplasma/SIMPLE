module provider_jets_impl
    !> L1 provider-jet identities. Every check uses check_identity (exact zero
    !> test, probe evidence refused). Source expressions are parsed from cited
    !> physical lines of the pinned files given on the command line:
    !>   1: libneo src/interpolate/batch_interpolate_3d.f90 (d30174f)
    !>   2: SIMPLE-fo-gauss2 src/spectre_fo_gauss2.f90 at 3ddc30b
    !>   3: libneo src/spectre/spectre_basis.f90 (d30174f)
    !> run_provider_jets.sh checks their SHA-256 before running this program.
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, func, operator(+), operator(-), &
        operator(*), operator(/), operator(**), exp, cos, sin
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs_many
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use provider_src, only: source_rhs, parse_rhs
    implicit none
    public
    integer, parameter :: NSEC = 16, LN = 48
    type(arena_t), target :: ar
    type(symengine_engine_t) :: eng
    type(suite_t) :: suite
    character(1024) :: f3d, fprov, fbasis
    character(64) :: secname(NSEC)
    integer :: secount(NSEC), isec, nsecs
    ! Global rename table for the libneo tensor routine (array element -> symbol).
    character(LN), allocatable :: tolds(:), tnews(:)

contains
    subroutine run_all
    call get_command_argument(1, f3d)
    call get_command_argument(2, fprov)
    call get_command_argument(3, fbasis)
    if (len_trim(fbasis) == 0) error stop 'usage: gen_provider_jets f3d fprov fbasis'
    call ar%init()
    eng = make_symengine_engine(ar)
    secount = 0
    nsecs = 0
    call suite_begin(suite, 'L1 provider-jet identities')

    call section('spline_horner_symbolic_degree')
    call horner_induction
    call section('spline_der2_nq1_general_all_orders')
    call tensor_der2_all_orders
    call section('spline_der3_basis_all_orders')
    call tensor_der3_all_orders
    call section('spline_affine_chain_and_slots')
    call affine_chain
    call section('gram_cartesian_embedding')
    call gram_cartesian
    call section('gram_cylindrical_3ddc30b')
    call gram_cylindrical
    call section('blend_potential_3ddc30b')
    call blend_potential
    call section('blend_embedding_leibniz_3ddc30b')
    call blend_leibniz
    call section('blend_metric_convexity')
    call blend_metric_pd
    call section('switch_der_3ddc30b')
    call switch_chain
    call section('gauge_projection_3ddc30b')
    call gauge_projection
    call section('spectral_radial_recurrences_libneo')
    call spectral_radial
    call section('spectral_fourier_modes_libneo')
    call spectral_modes
    call section('spectral_volume_chain_3ddc30b')
    call spectral_volume_chain
    call section('spectral_geometry_fourier_3ddc30b')
    call geometry_fourier

    call write_results('results_symbolic.json')
    call suite_end(suite, 'build/fortsym_suite.json')

    end subroutine run_all

    ! ------------------------------------------------------------- utilities

    function s(name) result(e)
        character(*), intent(in) :: name
        type(expr_t) :: e
        e = sym(ar, name)
    end function s

    function n(i) result(e)
        integer, intent(in) :: i
        type(expr_t) :: e
        e = num(ar, i)
    end function n

    function src(file, line, olds, news) result(e)
        character(*), intent(in) :: file
        integer, intent(in) :: line
        character(*), intent(in) :: olds(:), news(:)
        type(expr_t) :: e
        e = parse_rhs(ar, source_rhs(trim(file), line), olds, news)
    end function src

    !> RHS text equality of two cited lines (structural reuse of a formula).
    subroutine same_text(file, l1, l2)
        character(*), intent(in) :: file
        integer, intent(in) :: l1, l2
        if (source_rhs(trim(file), l1) /= source_rhs(trim(file), l2)) then
            print *, 'cited lines differ: ', l1, l2
            error stop 'same_text'
        end if
    end subroutine same_text

    subroutine section(name)
        character(*), intent(in) :: name
        nsecs = nsecs + 1
        isec = nsecs
        secname(isec) = name
        print '(a)', '== '//name
    end subroutine section

    subroutine ident(label, e)
        character(*), intent(in) :: label
        type(expr_t), intent(in) :: e
        integer :: before
        before = suite%failed
        call check_identity(suite, eng, label, e)
        secount(isec) = secount(isec) + 1
        if (suite%failed /= before) print *, 'FAILED: ', label
    end subroutine ident

    subroutine write_results(path)
        character(*), intent(in) :: path
        integer :: u, k
        open (newunit=u, file=path, status='replace')
        write (u, '(a)') '{'
        write (u, '(a)') '  "suite": "L1 provider-jet identities",'
        write (u, '(a,i0,a)') '  "total": ', suite%total, ','
        write (u, '(a,i0,a)') '  "proved_exact": ', suite%proved, ','
        write (u, '(a,i0,a)') '  "by_probe": ', suite%by_probe, ','
        write (u, '(a,i0,a)') '  "failed": ', suite%failed, ','
        write (u, '(a)') '  "sections": {'
        do k = 1, nsecs
            if (k < nsecs) then
                write (u, '(a,a,a,i0,a)') '    "', trim(secname(k)), '": ', secount(k), ','
            else
                write (u, '(a,a,a,i0)') '    "', trim(secname(k)), '": ', secount(k)
            end if
        end do
        write (u, '(a)') '  }'
        write (u, '(a)') '}'
        close (u)
    end subroutine write_results

    !> Replace derivative nodes of f (up to order `order`, all three variables)
    !> by jet symbols prefix_<slot>, highest order first. Mixed partials are
    !> canonicalised by fortsym, so one symbol covers both orders.
    function jet_subs(e, f, q, prefix, order) result(r)
        type(expr_t), intent(in) :: e, f, q(:)
        character(*), intent(in) :: prefix
        integer, intent(in) :: order
        type(expr_t) :: r
        type(expr_t), allocatable :: olds(:), news(:)
        integer :: i, j, k, m
        character(32) :: nm
        r = e
        if (order >= 3) then
            allocate (olds(0), news(0))
            do i = 1, size(q)
                do j = i, size(q)
                    do k = j, size(q)
                        write (nm, '(a,3i0)') prefix, i, j, k
                        olds = [olds, diff(diff(diff(f, q(i)), q(j)), q(k))]
                        news = [news, s(trim(nm))]
                    end do
                end do
            end do
            r = subs_many(r, olds, news)
            deallocate (olds, news)
        end if
        if (order >= 2) then
            allocate (olds(0), news(0))
            do i = 1, size(q)
                do j = i, size(q)
                    write (nm, '(a,2i0)') prefix, i, j
                    olds = [olds, diff(diff(f, q(i)), q(j))]
                    news = [news, s(trim(nm))]
                end do
            end do
            r = subs_many(r, olds, news)
            deallocate (olds, news)
        end if
        allocate (olds(0), news(0))
        do i = 1, size(q)
            write (nm, '(a,i0)') prefix, i
            olds = [olds, diff(f, q(i))]
            news = [news, s(trim(nm))]
        end do
        write (nm, '(a)') prefix
        olds = [olds, f]
        news = [news, s(trim(nm)//'0')]
        r = subs_many(r, olds, news)
        m = 0
    end function jet_subs

    ! ------------------------------------------- libneo spline: 1D Horner

    !> Horner value / first / second derivative accumulators of the libneo
    !> x1 reduction (batch_interpolate_3d.f90 lines 2159-2195), for symbolic
    !> degree. With U_k = sum_{j>=k} c_j x^(j-k), D_k = sum j c_j x^(j-k),
    !> S_k = sum j(j-1) c_j x^(j-k), the invariants
    !>   x U_k' = D_k - k U_k,   x^2 U_k'' = S_k - 2k D_k + k(k+1) U_k
    !> hold at k=N and are preserved by the code's step k+1 -> k. At the end
    !> the code returns U_0, D_1, S_2, which are U_0, U_0', U_0''.
    subroutine horner_induction
        type(expr_t) :: x, k, c, nn, u, du, d2u, dn, sn, uk, dk, sk, e, t1, t2, t3
        type(expr_t) :: hyp1, hyp2, u0, d1, s2
        character(LN) :: o(3), w(3)
        x = s('x1'); k = s('k'); c = s('c'); nn = s('N')
        o = [character(LN) :: 'coeff_23_dx1x1(k2, k3)', 'coeff_23_dx1(k2, k3)', &
            'coeff_23(k2, k3)']
        w = [character(LN) :: 'OLDS', 'OLDD', 'OLDU']
        ! Base k=N (lines 2159-2161): U_N=c, D_N=N c, S_N=N(N-1) c, U_N'=0.
        t1 = subs_many(src(f3d, 2159, o, w), [s('N1')], [nn])
        t2 = subs_many(src(f3d, 2160, o, w), [s('N1')], [nn])
        t3 = subs_many(src(f3d, 2161, o, w), [s('N1')], [nn])
        call ident('horner base first invariant', x*n(0) - (t2 - nn*t1))
        call ident('horner base second invariant', &
            x**2*n(0) - (t3 - 2*nn*t2 + nn*(nn + 1)*t1))
        ! Step (lines 2170-2172) with U_{k+1}=U(x) symbolic and hypotheses.
        u = func('U', [x]); du = diff(u, x); d2u = diff(du, x)
        dn = s('Dn'); sn = s('Sn')
        uk = subs_many(src(f3d, 2170, o, w), [s('OLDU'), s('k1')], [u, k])
        dk = subs_many(src(f3d, 2171, o, w), [s('OLDD'), s('k1')], [dn, k])
        sk = subs_many(src(f3d, 2172, o, w), [s('OLDS'), s('k1')], [sn, k])
        hyp1 = (dn - (k + 1)*u)/x
        hyp2 = (sn - 2*(k + 1)*dn + (k + 1)*(k + 2)*u)/x**2
        e = x*diff(uk, x) - (dk - k*uk)
        e = subs_many(subs_many(e, [d2u, du], [hyp2, hyp1]), [u], [s('Un')])
        call ident('horner step first invariant (symbolic k)', e)
        e = x**2*diff(diff(uk, x), x) - (sk - 2*k*dk + k*(k + 1)*uk)
        e = subs_many(e, [d2u], [hyp2])
        e = subs_many(e, [du], [hyp1])
        e = subs_many(e, [u], [s('Un')])
        call ident('horner step second invariant (symbolic k)', e)
        ! k1=1 lines 2184/2185 are the same formulas; k1=1 skips the second
        ! accumulator, which is exact since k1*(k1-1)*c = 0.
        call same_text(f3d, 2170, 2184)
        call same_text(f3d, 2171, 2185)
        call same_text(f3d, 2170, 2195)
        call ident('horner k=1 second step is x*S_2', &
            subs_many(src(f3d, 2172, o, w), [s('k1')], [n(1)]) - x*s('OLDS'))
        ! Termination: U_0 = c0 + x U_1 (line 2195); invariants at k=1 give
        ! U_0' = D_1 and U_0'' = S_2 (S_1 = x S_2).
        u0 = subs_many(src(f3d, 2195, o, w), [s('OLDU')], [u])
        d1 = s('D1'); s2 = s('S2')
        hyp1 = (d1 - u)/x
        hyp2 = (x*s2 - 2*d1 + 2*u)/x**2
        e = subs_many(subs_many(diff(u0, x) - d1, [du], [hyp1]), [u], [s('U1')])
        call ident('horner termination U0 prime equals D1', e)
        e = diff(diff(u0, x), x) - s2
        e = subs_many(e, [d2u], [hyp2])
        e = subs_many(e, [du], [hyp1])
        e = subs_many(e, [u], [s('U1')])
        call ident('horner termination U0 second equals S2', e)
    end subroutine horner_induction

    ! ------------------------------ libneo spline: tensor der2, nq1 general

    subroutine init_tensor_renames
        tolds = [character(LN) :: 'd2y_batch(1, 1)', 'd2y_batch(2, 1)', &
            'd2y_batch(3, 1)', 'd2y_batch(4, 1)', 'd2y_batch(5, 1)', &
            'd2y_batch(6, 1)', 'dy_batch(1, 1)', 'dy_batch(2, 1)', &
            'dy_batch(3, 1)', 'y_batch(1)', &
            'coeff_23_dx1x1(k2, k3)', 'coeff_23_dx1(k2, k3)', 'coeff_23(k2, k3)', &
            'coeff_3_dx1x1(k3)', 'coeff_3_dx1x2(k3)', 'coeff_3_dx2x2(k3)', &
            'coeff_3_dx1(k3)', 'coeff_3_dx2(k3)', 'coeff_3(k3)', &
            'coeff_3_dx1x1(N3)', 'coeff_3_dx1x2(N3)', 'coeff_3_dx2x2(N3)', &
            'coeff_3_dx1(N3)', 'coeff_3_dx2(N3)', 'coeff_3(N3)']
        tnews = [character(LN) :: 'Y11', 'Y12', 'Y13', 'Y22', 'Y23', 'Y33', &
            'Y1', 'Y2', 'Y3', 'Y', 'A11', 'A1', 'A0', &
            'B11', 'B12', 'B22', 'B1', 'B2', 'B0', &
            'B11', 'B12', 'B22', 'B1', 'B2', 'B0']
    end subroutine init_tensor_renames

    function tline(line) result(e)
        integer, intent(in) :: line
        type(expr_t) :: e
        e = src(f3d, line, tolds, tnews)
    end function tline

    !> Transcribed control flow of evaluate_batch_splines_3d_der2_core_nq1
    !> (lines 2103-2291, general branch) with every assignment RHS parsed from
    !> its cited line. Orders 3..5 are all orders construct accepts.
    subroutine tensor_der2_all_orders
        integer :: n1, n2, n3
        call init_tensor_renames
        do n1 = 3, 5
            do n2 = 3, 5
                do n3 = 3, 5
                    call tensor_der2(n1, n2, n3)
                end do
            end do
        end do
    end subroutine tensor_der2_all_orders

    subroutine tensor_der2(o1, o2, o3)
        integer, intent(in) :: o1, o2, o3
        type(expr_t) :: cf(0:o1, 0:o2, 0:o3), x(3), p, truth(10), out(10)
        type(expr_t) :: a0(0:o2, 0:o3), a1(0:o2, 0:o3), a11(0:o2, 0:o3)
        type(expr_t) :: b0(0:o3), b1(0:o3), b11(0:o3), b2(0:o3), b12(0:o3), b22(0:o3)
        type(expr_t) :: t(2159:2290), names(26), vals(26)
        integer :: k1, k2, k3, i, j, m, line
        integer, parameter :: used(*) = [2159, 2160, 2161, 2170, 2171, 2172, 2184, &
            2185, 2195, 2204, 2205, 2206, 2212, 2213, 2214, 2216, 2217, 2219, 2228, &
            2229, 2230, 2232, 2233, 2241, 2242, 2243, 2247, 2248, 2249, 2250, 2251, &
            2252, 2253, 2254, 2255, 2256, 2259, 2260, 2261, 2262, 2263, 2264, 2265, &
            2266, 2267, 2268, 2273, 2274, 2275, 2276, 2277, 2278, 2279, 2280, 2281, &
            2285, 2286, 2287, 2288, 2289, 2290]
        character(64) :: lab
        do i = 1, size(used)
            line = used(i)
            t(line) = tline(line)
        end do
        call init_names
        x = [s('x1'), s('x2'), s('x3')]
        p = n(0)
        do k3 = 0, o3
            do k2 = 0, o2
                do k1 = 0, o1
                    write (lab, '(a,3(a,i0))') 'c', '_', k1, '_', k2, '_', k3
                    cf(k1, k2, k3) = s(trim(lab))
                    p = p + cf(k1, k2, k3)*x(1)**k1*x(2)**k2*x(3)**k3
                end do
            end do
        end do
        ! First reduction (lines 2156-2197).
        do k3 = 0, o3
            do k2 = 0, o2
                call setc(cf(o1, k2, k3), o1, k2, k3)
                a0(k2, k3) = ev(t(2159)); a1(k2, k3) = ev(t(2160))
                a11(k2, k3) = ev(t(2161))
            end do
        end do
        do k1 = o1 - 1, 2, -1
            do k3 = 0, o3
                do k2 = 0, o2
                    call setc(cf(k1, k2, k3), k1, k2, k3)
                    call seta(a0(k2, k3), a1(k2, k3), a11(k2, k3))
                    a0(k2, k3) = ev(t(2170)); a1(k2, k3) = ev(t(2171))
                    a11(k2, k3) = ev(t(2172))
                end do
            end do
        end do
        k1 = 1
        do k3 = 0, o3
            do k2 = 0, o2
                call setc(cf(k1, k2, k3), k1, k2, k3)
                call seta(a0(k2, k3), a1(k2, k3), a11(k2, k3))
                a0(k2, k3) = ev(t(2184)); a1(k2, k3) = ev(t(2185))
            end do
        end do
        k1 = 0
        do k3 = 0, o3
            do k2 = 0, o2
                call setc(cf(k1, k2, k3), k1, k2, k3)
                call seta(a0(k2, k3), a1(k2, k3), a11(k2, k3))
                a0(k2, k3) = ev(t(2195))
            end do
        end do
        ! Second reduction (lines 2200-2244). Lines 2201-2203 copy row N2.
        do k3 = 0, o3
            b0(k3) = a0(o2, k3); b1(k3) = a1(o2, k3); b11(k3) = a11(o2, k3)
            call setb(b0(k3), b1(k3), b11(k3), b0(k3), b0(k3), b0(k3), o2, k3)
            b2(k3) = ev(t(2204)); b12(k3) = ev(t(2205)); b22(k3) = ev(t(2206))
        end do
        do k2 = o2 - 1, 2, -1
            do k3 = 0, o3
                call seta(a0(k2, k3), a1(k2, k3), a11(k2, k3))
                call setb(b0(k3), b1(k3), b11(k3), b2(k3), b12(k3), b22(k3), k2, k3)
                b0(k3) = ev(t(2212)); b1(k3) = ev(t(2213)); b11(k3) = ev(t(2214))
                b2(k3) = ev(t(2216)); b12(k3) = ev(t(2217)); b22(k3) = ev(t(2219))
            end do
        end do
        k2 = 1
        do k3 = 0, o3
            call seta(a0(k2, k3), a1(k2, k3), a11(k2, k3))
            call setb(b0(k3), b1(k3), b11(k3), b2(k3), b12(k3), b22(k3), k2, k3)
            b0(k3) = ev(t(2228)); b1(k3) = ev(t(2229)); b11(k3) = ev(t(2230))
            b2(k3) = ev(t(2232)); b12(k3) = ev(t(2233))
        end do
        k2 = 0
        do k3 = 0, o3
            call seta(a0(k2, k3), a1(k2, k3), a11(k2, k3))
            call setb(b0(k3), b1(k3), b11(k3), b2(k3), b12(k3), b22(k3), k2, k3)
            b0(k3) = ev(t(2241)); b1(k3) = ev(t(2242)); b11(k3) = ev(t(2243))
        end do
        ! Third reduction (lines 2246-2290).
        call setb(b0(o3), b1(o3), b11(o3), b2(o3), b12(o3), b22(o3), o2, o3)
        call sety([(n(0), m=1, 10)], o3)
        do m = 1, 10
            out(m) = ev(t(2246 + m))
        end do
        do k3 = o3 - 1, 2, -1
            call setb(b0(k3), b1(k3), b11(k3), b2(k3), b12(k3), b22(k3), 0, k3)
            call sety(out, k3)
            do m = 1, 10
                out(m) = ev(t(2258 + m))
            end do
        end do
        k3 = 1
        call setb(b0(k3), b1(k3), b11(k3), b2(k3), b12(k3), b22(k3), 0, k3)
        call sety(out, k3)
        do m = 1, 9
            out(m) = ev(t(2272 + m))
        end do
        k3 = 0
        call setb(b0(k3), b1(k3), b11(k3), b2(k3), b12(k3), b22(k3), 0, k3)
        call sety(out, k3)
        out(1) = ev(t(2285)); out(2) = ev(t(2286)); out(3) = ev(t(2287))
        out(5) = ev(t(2288)); out(6) = ev(t(2289)); out(8) = ev(t(2290))
        ! out ordering follows lines 2247-2256: y, dy1, dy2, dy3, d2y(1..6).
        truth(1) = p
        do i = 1, 3
            truth(1 + i) = diff(p, x(i))
        end do
        m = 4
        do i = 1, 3
            do j = i, 3
                m = m + 1
                truth(m) = diff(diff(p, x(i)), x(j))
            end do
        end do
        do m = 1, 10
            write (lab, '(a,3i1,a,i0)') 'der2 nq1 order ', o1, o2, o3, ' output ', m
            call ident(trim(lab), out(m) - truth(m))
        end do
    contains
        subroutine init_names
            character(4), parameter :: nm(26) = [character(4) :: 'c', 'k1', 'k2', &
                'k3', 'N1', 'N2', 'N3', 'A0', 'A1', 'A11', 'B0', 'B1', 'B11', 'B2', &
                'B12', 'B22', 'Y', 'Y1', 'Y2', 'Y3', 'Y11', 'Y12', 'Y13', 'Y22', &
                'Y23', 'Y33']
            integer :: q
            do q = 1, 26
                names(q) = s(trim(nm(q)))
                vals(q) = n(0)
            end do
            vals(5) = n(o1); vals(6) = n(o2); vals(7) = n(o3)
        end subroutine init_names
        subroutine setc(cv, i1, i2, i3)
            type(expr_t), intent(in) :: cv
            integer, intent(in) :: i1, i2, i3
            vals(1) = cv; vals(2) = n(i1); vals(3) = n(i2); vals(4) = n(i3)
        end subroutine setc
        subroutine seta(v0, v1, v11)
            type(expr_t), intent(in) :: v0, v1, v11
            vals(8) = v0; vals(9) = v1; vals(10) = v11
        end subroutine seta
        subroutine setb(v0, v1, v11, v2, v12, v22, i2, i3)
            type(expr_t), intent(in) :: v0, v1, v11, v2, v12, v22
            integer, intent(in) :: i2, i3
            vals(3) = n(i2); vals(4) = n(i3)
            vals(11) = v0; vals(12) = v1; vals(13) = v11
            vals(14) = v2; vals(15) = v12; vals(16) = v22
        end subroutine setb
        subroutine sety(v, i3)
            type(expr_t), intent(in) :: v(10)
            integer, intent(in) :: i3
            vals(4) = n(i3)
            vals(17:26) = v
        end subroutine sety
        function ev(tt) result(r)
            type(expr_t), intent(in) :: tt
            type(expr_t) :: r
            r = subs_many(tt, names, vals)
        end function ev
    end subroutine tensor_der2

    ! ------------------------------ libneo spline: der3 monomial basis

    subroutine tensor_der3_all_orders
        integer :: n1, n2, n3
        type(expr_t) :: xv, kk, prev
        ! Basis lines 1829/1832/1835/1838 for symbolic kk: with vb0(kk-j) =
        ! xv**(kk-j) the RHS equal the derivatives of xv**kk. Symbolic powers
        ! are avoided (see README: fortsym zero test on symbolic exponents),
        ! so this is checked for every kk = 1..5 the routine can reach.
        integer :: k
        character(LN) :: o(3), w(3)
        character(64) :: lab
        o = [character(LN) :: 'vb0(kk - 1)', 'vb0(kk - 2)', 'vb0(kk - 3)']
        w = [character(LN) :: 'P1', 'P2', 'P3']
        xv = s('xv'); kk = s('kk')
        do k = 1, 5
            prev = xv**(k - 1)
            write (lab, '(a,i0)') 'basis_d3 vb0 kk=', k
            call ident(trim(lab), subs_many(src(f3d, 1829, o, w), [s('P1'), kk], &
                [prev, n(k)]) - xv**k)
            write (lab, '(a,i0)') 'basis_d3 vb1 kk=', k
            call ident(trim(lab), subs_many(src(f3d, 1832, o, w), [s('P1'), kk], &
                [prev, n(k)]) - diff(xv**k, xv))
            if (k >= 2) then
                write (lab, '(a,i0)') 'basis_d3 vb2 kk=', k
                call ident(trim(lab), subs_many(src(f3d, 1835, o, w), [s('P2'), kk], &
                    [xv**(k - 2), n(k)]) - diff(diff(xv**k, xv), xv))
            end if
            if (k >= 3) then
                write (lab, '(a,i0)') 'basis_d3 vb3 kk=', k
                call ident(trim(lab), subs_many(src(f3d, 1838, o, w), [s('P3'), kk], &
                    [xv**(k - 3), n(k)]) - diff(diff(diff(xv**k, xv), xv), xv))
            end if
        end do
        do n1 = 3, 5
            do n2 = 3, 5
                do n3 = 3, 5
                    call tensor_der3(n1, n2, n3)
                end do
            end do
        end do
    end subroutine tensor_der3_all_orders

    !> Accumulation lines 1791-1810 over all cell coefficients, with the
    !> basis values b<d>(k,j) = d-th derivative of x_j**k (proved above).
    subroutine tensor_der3(o1, o2, o3)
        integer, intent(in) :: o1, o2, o3
        type(expr_t) :: x(3), p, t(20), out(20), truth(20), c, names(13), vals(13)
        type(expr_t) :: bas(0:3, 0:5, 3)
        integer :: k1, k2, k3, m, i, j, k, d, ord(3)
        character(LN) :: o(32), w(32)
        character(64) :: lab
        character(4), parameter :: yn(20) = [character(4) :: 'Y', 'Y1', 'Y2', 'Y3', &
            'Y11', 'Y12', 'Y13', 'Y22', 'Y23', 'Y33', 'Z1', 'Z2', 'Z3', 'Z4', 'Z5', &
            'Z6', 'Z7', 'Z8', 'Z9', 'Z10']
        ord = [o1, o2, o3]
        o(1:20) = [character(LN) :: 'y_batch(q)', 'dy_batch(1, q)', 'dy_batch(2, q)', &
            'dy_batch(3, q)', 'd2y_batch(1, q)', 'd2y_batch(2, q)', 'd2y_batch(3, q)', &
            'd2y_batch(4, q)', 'd2y_batch(5, q)', 'd2y_batch(6, q)', &
            'd3y_batch(10, q)', 'd3y_batch(1, q)', 'd3y_batch(2, q)', &
            'd3y_batch(3, q)', 'd3y_batch(4, q)', 'd3y_batch(5, q)', &
            'd3y_batch(6, q)', 'd3y_batch(7, q)', 'd3y_batch(8, q)', 'd3y_batch(9, q)']
        w(1:20) = [character(LN) :: 'OLD', 'OLD', 'OLD', 'OLD', 'OLD', 'OLD', 'OLD', &
            'OLD', 'OLD', 'OLD', 'OLD', 'OLD', 'OLD', 'OLD', 'OLD', 'OLD', 'OLD', &
            'OLD', 'OLD', 'OLD']
        m = 20
        do j = 1, 3
            do d = 0, 3
                m = m + 1
                write (o(m), '(a,i0,a,i0,a,i0,a)') 'b', d, '(k', j, ', ', j, ')'
                write (w(m), '(a,i0,a,i0)') 'BB', d, '_', j
            end do
        end do
        do m = 1, 20
            t(m) = src(f3d, 1790 + m, o, w)
        end do
        x = [s('x1'), s('x2'), s('x3')]
        do j = 1, 3
            do k = 0, ord(j)
                bas(0, k, j) = x(j)**k
                do d = 1, 3
                    bas(d, k, j) = diff(bas(d - 1, k, j), x(j))
                end do
            end do
        end do
        out = n(0)
        p = n(0)
        do k3 = 0, o3
            do k2 = 0, o2
                do k1 = 0, o1
                    write (lab, '(a,3(a,i0))') 'c', '_', k1, '_', k2, '_', k3
                    c = s(trim(lab))
                    p = p + c*x(1)**k1*x(2)**k2*x(3)**k3
                    names(1) = s('c'); vals(1) = c
                    m = 1
                    do j = 1, 3
                        do d = 0, 3
                            m = m + 1
                            write (lab, '(a,i0,a,i0)') 'BB', d, '_', j
                            names(m) = s(trim(lab))
                            select case (j)
                            case (1); vals(m) = bas(d, k1, 1)
                            case (2); vals(m) = bas(d, k2, 2)
                            case (3); vals(m) = bas(d, k3, 3)
                            end select
                        end do
                    end do
                    do i = 1, 20
                        out(i) = subs_many(t(i), [names, s('OLD')], [vals, out(i)])
                    end do
                end do
            end do
        end do
        truth(1) = p
        do i = 1, 3
            truth(1 + i) = diff(p, x(i))
        end do
        m = 4
        do i = 1, 3
            do j = i, 3
                m = m + 1
                truth(m) = diff(diff(p, x(i)), x(j))
            end do
        end do
        do i = 1, 3
            do j = i, 3
                do k = j, 3
                    m = m + 1
                    truth(m) = diff(diff(diff(p, x(i)), x(j)), x(k))
                end do
            end do
        end do
        do m = 1, 20
            write (lab, '(a,3i1,a,a)') 'der3 order ', o1, o2, o3, ' output ', yn(m)
            call ident(trim(lab), out(m) - truth(m))
        end do
    end subroutine tensor_der3

    ! ------------------------------ chain rule for logical coordinates

    !> A spline in s = a*rho + b (or in theta, zeta shifted by a period) gives
    !> rho-jets with factors a, a**2; a periodic shift does not change jets.
    !> Slot order: libneo d2y (11,12,13,22,23,33) -> FO jet slots 5..10
    !> (11,12,22,13,23,33) is the permutation [1,2,4,3,5,6].
    subroutine affine_chain
        type(expr_t) :: rho, th, ze, a, b, per, f, h, fs(3), q(3), e
        integer :: i, j, m, slot, r, t
        integer :: fo_i(10), fo_j(10), lib_i(6), lib_j(6)
        integer, parameter :: perm(6) = [1, 2, 4, 3, 5, 6]
        rho = s('rho'); th = s('theta'); ze = s('zeta'); a = s('a'); b = s('b')
        per = s('per')
        q = [rho, th, ze]
        fs = [a*rho + b, th + per, ze]
        f = func('F', fs)
        h = func('F', [s('u'), s('v'), s('w')])
        ! d/drho F(a rho+b, theta+per, zeta) = a F_1 etc.
        do i = 1, 3
            e = diff(f, q(i))
            call ident('chain first '//achar(48 + i), &
                e - chain_factor(i)*slotjet(f, fs, i, 0))
        end do
        do i = 1, 3
            do j = i, 3
                e = diff(diff(f, q(i)), q(j))
                call ident('chain second '//achar(48 + i)//achar(48 + j), &
                    e - chain_factor(i)*chain_factor(j)*slotjet(f, fs, i, j))
            end do
        end do
        ! Slot bookkeeping, FO side copied from gen/gen_fo_hamiltonian.f90
        ! jet_replacements (slot 5.. by t=1..3, r=1..t) and provider packing.
        slot = 4
        do t = 1, 3
            do r = 1, t
                slot = slot + 1
                fo_i(slot) = r; fo_j(slot) = t
            end do
        end do
        m = 0
        do i = 1, 3
            do j = i, 3
                m = m + 1
                lib_i(m) = i; lib_j(m) = j
            end do
        end do
        do m = 1, 6
            if (fo_i(4 + m) /= lib_i(perm(m)) .or. fo_j(4 + m) /= lib_j(perm(m))) &
                error stop 'slot permutation'
        end do
        print *, 'slot permutation libneo d2y -> FO slots 5..10 = [1,2,4,3,5,6] ok'
    end subroutine affine_chain

    function chain_factor(i) result(e)
        integer, intent(in) :: i
        type(expr_t) :: e
        if (i == 1) then
            e = s('a')
        else
            e = n(1)
        end if
    end function chain_factor

    !> Derivative of F with respect to its own argument slots i (and j),
    !> evaluated at the shifted arguments, obtained by differentiating
    !> F(u,v,w) and substituting the arguments.
    function slotjet(f, fs, i, j) result(e)
        type(expr_t), intent(in) :: f, fs(3)
        integer, intent(in) :: i, j
        type(expr_t) :: e, g, u(3)
        u = [s('u'), s('v'), s('w')]
        g = func('F', u)
        e = diff(g, u(i))
        if (j > 0) e = diff(e, u(j))
        e = subs_many(e, u, fs)
    end function slotjet

    ! ------------------------------ Gram metric of a Cartesian embedding

    !> g_ij = sum_a d_i X_a d_j X_a; dg and d2g are the formulas a spline
    !> provider must assemble from embedding jets to third order.
    subroutine gram_cartesian
        type(expr_t) :: q(3), xf(3), g, e, dgf, d2gf
        integer :: i, j, k, l, a
        character(32) :: lab
        q = [s('q1'), s('q2'), s('q3')]
        do a = 1, 3
            xf(a) = func('X'//achar(48 + a), q)
        end do
        do j = 1, 3
            do i = 1, j
                g = n(0)
                do a = 1, 3
                    g = g + diff(xf(a), q(i))*diff(xf(a), q(j))
                end do
                do k = 1, 3
                    dgf = n(0)
                    do a = 1, 3
                        dgf = dgf + xjet(a, [i, k])*xjet(a, [j]) + &
                            xjet(a, [i])*xjet(a, [j, k])
                    end do
                    e = jets3(diff(g, q(k)) , xf, q) - dgf
                    write (lab, '(a,3i1)') 'cart dg ', i, j, k
                    call ident(trim(lab), e)
                    do l = k, 3
                        d2gf = n(0)
                        do a = 1, 3
                            d2gf = d2gf + xjet(a, [i, k, l])*xjet(a, [j]) + &
                                xjet(a, [i, k])*xjet(a, [j, l]) + &
                                xjet(a, [i, l])*xjet(a, [j, k]) + &
                                xjet(a, [i])*xjet(a, [j, k, l])
                        end do
                        e = jets3(diff(diff(g, q(k)), q(l)), xf, q) - d2gf
                        write (lab, '(a,4i1)') 'cart d2g ', i, j, k, l
                        call ident(trim(lab), e)
                    end do
                end do
            end do
        end do
    end subroutine gram_cartesian

    !> Jet symbol X<a>_<sorted multi-index>.
    function xjet(a, idx) result(e)
        integer, intent(in) :: a, idx(:)
        type(expr_t) :: e
        integer :: srt(size(idx)), i, j, tmp
        character(32) :: nm
        srt = idx
        do i = 1, size(srt)
            do j = i + 1, size(srt)
                if (srt(j) < srt(i)) then
                    tmp = srt(i); srt(i) = srt(j); srt(j) = tmp
                end if
            end do
        end do
        write (nm, '(a,i0,a,*(i0))') 'X', a, '_', srt
        e = s(trim(nm))
    end function xjet

    function jets3(e, xf, q) result(r)
        type(expr_t), intent(in) :: e, xf(3), q(3)
        type(expr_t) :: r
        integer :: a
        character(8) :: pre
        r = e
        do a = 1, 3
            write (pre, '(a,i0,a)') 'X', a, '_'
            r = jet_subs(r, xf(a), q, trim(pre), 3)
        end do
    end function jets3

    ! ------------------------------ cylindrical Gram, 3ddc30b lines 319-345

    !> X = (R cos zeta, R sin zeta, Z) with R,Z functions of q; d(:,i), d2,
    !> d3 are their jets (lines 304-318). Lines 321-345 transcribed with the
    !> two-component dot_product written out; the cited text is checked.
    subroutine gram_cylindrical
        type(expr_t) :: q(3), rf, zf, xc(3), g, e, d(2, 3), d2(2, 3, 3), d3(2, 3, 3, 3)
        type(expr_t) :: rz(2), metric, dg, d2g, ev(3)
        integer :: i, j, k, l, c
        character(40) :: lab
        call expect(fprov, 373, 'x = [rz(1)*cz, rz(1)*sz, rz(2)]')
        call expect(fprov, 375, 'e(:, i) = [d(1, i)*cz, d(1, i)*sz, d(2, i)]')
        call expect(fprov, 377, 'e(:, 3) = e(:, 3) + [-rz(1)*sz, rz(1)*cz, 0.0_dp]')
        call expect(fprov, 378, 'metric = matmul(transpose(d), d)')
        call expect(fprov, 379, 'metric(3, 3) = metric(3, 3) + rz(1)**2')
        call expect(fprov, 384, 'dg(i, j, k) = dot_product(d2(:, i, k), d(:, j)) &')
        call expect(fprov, 385, '+ dot_product(d(:, i), d2(:, j, k))')
        call expect(fprov, 388, &
            'd2g(i, j, k, l) = dot_product(d3(:, i, k, l), d(:, j)) &')
        call expect(fprov, 389, '+ dot_product(d2(:, i, k), d2(:, j, l)) &')
        call expect(fprov, 390, '+ dot_product(d2(:, i, l), d2(:, j, k)) &')
        call expect(fprov, 391, '+ dot_product(d(:, i), d3(:, j, k, l))')
        call expect(fprov, 395, 'dg(3, 3, k) = dg(3, 3, k) + 2.0_dp*rz(1)*d(1, k)')
        call expect(fprov, 398, 'd2g(3, 3, k, l) = d2g(3, 3, k, l) + &')
        call expect(fprov, 399, '2.0_dp*(d(1, k)*d(1, l) + rz(1)*d2(1, k, l))')
        q = [s('rho'), s('theta'), s('zeta')]
        rf = func('R', q); zf = func('Z', q)
        xc = [rf*cos(q(3)), rf*sin(q(3)), zf]
        rz = [s('R_0'), s('Z_0')]
        do i = 1, 3
            d(1, i) = cyl(diff(rf, q(i)), rf, zf, q)
            d(2, i) = cyl(diff(zf, q(i)), rf, zf, q)
            do j = 1, 3
                do c = 1, 2
                    d2(c, i, j) = cyl(diff(diff(pick(c, rf, zf), q(i)), q(j)), rf, zf, q)
                    do k = 1, 3
                        d3(c, i, j, k) = cyl(diff(diff(diff(pick(c, rf, zf), q(i)), &
                            q(j)), q(k)), rf, zf, q)
                    end do
                end do
            end do
        end do
        ! Lines 321-325: x and e(:,i) = d x / d q_i.
        ev = [rz(1)*cos(q(3)), rz(1)*sin(q(3)), rz(2)]
        do c = 1, 3
            call ident('cyl x', cyl(xc(c), rf, zf, q) - ev(c))
        end do
        do i = 1, 3
            ev = [d(1, i)*cos(q(3)), d(1, i)*sin(q(3)), d(2, i)]
            if (i == 3) then
                ev(1) = ev(1) - rz(1)*sin(q(3)); ev(2) = ev(2) + rz(1)*cos(q(3))
            end if
            do c = 1, 3
                write (lab, '(a,2i1)') 'cyl basis e ', c, i
                call ident(trim(lab), cyl(diff(xc(c), q(i)), rf, zf, q) - ev(c))
            end do
        end do
        do j = 1, 3
            do i = 1, 3
                g = n(0)
                do c = 1, 3
                    g = g + diff(xc(c), q(i))*diff(xc(c), q(j))
                end do
                metric = d(1, i)*d(1, j) + d(2, i)*d(2, j)
                if (i == 3 .and. j == 3) metric = metric + rz(1)**2
                write (lab, '(a,2i1)') 'cyl metric ', i, j
                call ident(trim(lab), cyl(g, rf, zf, q) - metric)
                do k = 1, 3
                    dg = d2(1, i, k)*d(1, j) + d2(2, i, k)*d(2, j) + &
                        d(1, i)*d2(1, j, k) + d(2, i)*d2(2, j, k)
                    if (i == 3 .and. j == 3) dg = dg + 2*rz(1)*d(1, k)
                    write (lab, '(a,3i1)') 'cyl dg ', i, j, k
                    call ident(trim(lab), cyl(diff(g, q(k)), rf, zf, q) - dg)
                    do l = 1, 3
                        d2g = d3(1, i, k, l)*d(1, j) + d3(2, i, k, l)*d(2, j) + &
                            d2(1, i, k)*d2(1, j, l) + d2(2, i, k)*d2(2, j, l) + &
                            d2(1, i, l)*d2(1, j, k) + d2(2, i, l)*d2(2, j, k) + &
                            d(1, i)*d3(1, j, k, l) + d(2, i)*d3(2, j, k, l)
                        if (i == 3 .and. j == 3) d2g = d2g + &
                            2*(d(1, k)*d(1, l) + rz(1)*d2(1, k, l))
                        write (lab, '(a,4i1)') 'cyl d2g ', i, j, k, l
                        call ident(trim(lab), &
                            cyl(diff(diff(g, q(k)), q(l)), rf, zf, q) - d2g)
                    end do
                end do
            end do
        end do
    end subroutine gram_cylindrical

    function pick(c, a, b) result(e)
        integer, intent(in) :: c
        type(expr_t), intent(in) :: a, b
        type(expr_t) :: e
        e = a
        if (c == 2) e = b
    end function pick

    function cyl(e, rf, zf, q) result(r)
        type(expr_t), intent(in) :: e, rf, zf, q(3)
        type(expr_t) :: r
        r = jet_subs(e, rf, q, 'R_', 3)
        r = jet_subs(r, zf, q, 'Z_', 3)
    end function cyl

    subroutine expect(file, line, text)
        character(*), intent(in) :: file, text
        integer, intent(in) :: line
        character(:), allocatable :: got
        got = trim(adjustl(source_line_text(file, line)))
        if (got /= text) then
            print *, 'line ', line, ' is: ', got
            print *, 'expected:   ', text
            error stop 'transcription anchor mismatch'
        end if
    end subroutine expect

    function source_line_text(file, line) result(t)
        use provider_src, only: source_line
        character(*), intent(in) :: file
        integer, intent(in) :: line
        character(:), allocatable :: t
        t = source_line(trim(file), line)
    end function source_line_text

    ! ------------------------------ blends, 3ddc30b

    !> Lines 139-151: A = al + s0 (ar - al) with s(1)=w', s(2)=w''.
    subroutine blend_potential
        type(expr_t) :: q(3), al, arr, w, a, e, sw(0:2), dl(3), dr(3), hl(3, 3), hr(3, 3)
        type(expr_t) :: da(3), d2a(3, 3), tv
        integer :: i, j
        character(32) :: lab
        call expect(fprov, 157, 'a = al + s(0)*(ar - al)')
        call expect(fprov, 158, 'da = dl + s(0)*(dr - dl)')
        call expect(fprov, 159, 'da(:, 1) = da(:, 1) + s(1)*(ar - al)')
        call expect(fprov, 166, 'd2a = hl + s(0)*(hr - hl)')
        call expect(fprov, 169, 'if (i == 1) d2a(:, i, j) = d2a(:, i, j) &')
        call expect(fprov, 170, '+ s(1)*(dr(:, j) - dl(:, j))')
        call expect(fprov, 171, 'if (j == 1) d2a(:, i, j) = d2a(:, i, j) &')
        call expect(fprov, 172, '+ s(1)*(dr(:, i) - dl(:, i))')
        call expect(fprov, 175, 'd2a(:, 1, 1) = d2a(:, 1, 1) + s(2)*(ar - al)')
        q = [s('rho'), s('theta'), s('zeta')]
        al = func('AL', q); arr = func('AR', q); w = func('W', [q(1)])
        a = al + w*(arr - al)
        sw(0) = s('W_0'); sw(1) = s('W_1'); sw(2) = s('W_2')
        do i = 1, 3
            dl(i) = s('AL_'//achar(48 + i)); dr(i) = s('AR_'//achar(48 + i))
            do j = 1, 3
                hl(i, j) = s('AL_'//achar(48 + min(i, j))//achar(48 + max(i, j)))
                hr(i, j) = s('AR_'//achar(48 + min(i, j))//achar(48 + max(i, j)))
            end do
        end do
        call ident('blend value', bl(a) - (s('AL_0') + sw(0)*(s('AR_0') - s('AL_0'))))
        do i = 1, 3
            da(i) = dl(i) + sw(0)*(dr(i) - dl(i))
            if (i == 1) da(i) = da(i) + sw(1)*(s('AR_0') - s('AL_0'))
            write (lab, '(a,i1)') 'blend da ', i
            call ident(trim(lab), bl(diff(a, q(i))) - da(i))
        end do
        do j = 1, 3
            do i = 1, 3
                d2a(i, j) = hl(i, j) + sw(0)*(hr(i, j) - hl(i, j))
                if (i == 1) d2a(i, j) = d2a(i, j) + sw(1)*(dr(j) - dl(j))
                if (j == 1) d2a(i, j) = d2a(i, j) + sw(1)*(dr(i) - dl(i))
                if (i == 1 .and. j == 1) d2a(i, j) = d2a(i, j) + sw(2)*(s('AR_0') - s('AL_0'))
                write (lab, '(a,2i1)') 'blend d2a ', i, j
                call ident(trim(lab), bl(diff(diff(a, q(i)), q(j))) - d2a(i, j))
            end do
        end do
    contains
        function bl(e0) result(r)
            type(expr_t), intent(in) :: e0
            type(expr_t) :: r
            type(expr_t) :: wd1, wd2
            wd1 = diff(w, q(1)); wd2 = diff(wd1, q(1))
            r = subs_many(e0, [wd2], [sw(2)])
            r = subs_many(r, [wd1], [sw(1)])
            r = subs_many(r, [w], [sw(0)])
            r = jet_subs(r, al, q, 'AL_', 2)
            r = jet_subs(r, arr, q, 'AR_', 2)
        end function bl
    end subroutine blend_potential

    !> Lines 283-288: coeff(:,r) = cl(:,r) + sum_l binom(l,r) sw(l)(cr-cl)(r-l)
    !> is the r-th rho derivative of cl + w (cr - cl), r = 0..3 (Leibniz);
    !> binom from lines 260-261 (column-major reshape).
    subroutine blend_leibniz
        type(expr_t) :: rho, cl, cr, w, f, e, cf, jl(0:3), jr(0:3), sw(0:3), t
        integer :: r, l
        integer, parameter :: binom(0:3, 0:3) = reshape([ &
            1, 0, 0, 0, 1, 1, 0, 0, 1, 2, 1, 0, 1, 3, 3, 1], [4, 4])
        character(LN) :: o(4), wn(4)
        character(32) :: lab
        call expect(fprov, 298, &
            'integer, parameter :: binom(0:3, 0:3) = reshape([ &')
        call expect(fprov, 299, &
                    '1, 0, 0, 0, 1, 1, 0, 0, 1, 2, 1, 0, 1, 3, 3, 1], [4, 4])')
        o = [character(LN) :: 'coeff(:, r)', 'real(binom(l, r), dp)', 'sw(l)', &
            '(cr(:, r - l) - cl(:, r - l))']
        wn = [character(LN) :: 'OLD', 'BIN', 'SW', 'DIFF']
        rho = s('rho')
        cl = func('CL', [rho]); cr = func('CR', [rho]); w = func('W', [rho])
        f = cl + w*(cr - cl)
        jl(0) = cl; jr(0) = cr; sw(0) = w
        do r = 1, 3
            jl(r) = diff(jl(r - 1), rho); jr(r) = diff(jr(r - 1), rho)
            sw(r) = diff(sw(r - 1), rho)
        end do
        t = src(fprov, 333, o, wn)
        do r = 0, 3
            cf = jl(r)
            do l = 0, r
                cf = subs_many(t, [s('OLD'), s('BIN'), s('SW'), s('DIFF')], &
                    [cf, n(binom(l, r)), sw(l), jr(r - l) - jl(r - l)])
            end do
            e = f
            do l = 1, r
                e = diff(e, rho)
            end do
            write (lab, '(a,i1)') 'embedding blend Leibniz r=', r
            call ident(trim(lab), e - cf)
        end do
    end subroutine blend_leibniz

    !> v^T [(1-w) gL + w gR] v = (1-w) v^T gL v + w v^T gR v: a blend of
    !> positive-definite metrics with 0<=w<=1 is positive definite.
    subroutine blend_metric_pd
        type(expr_t) :: v(3), gl(3, 3), gr(3, 3), w, lhs, ql, qr
        integer :: i, j
        character(8) :: nm
        w = s('w')
        do i = 1, 3
            v(i) = s('v'//achar(48 + i))
            do j = 1, 3
                write (nm, '(a,2i1)') 'gL', min(i, j), max(i, j); gl(i, j) = s(trim(nm))
                write (nm, '(a,2i1)') 'gR', min(i, j), max(i, j); gr(i, j) = s(trim(nm))
            end do
        end do
        lhs = n(0); ql = n(0); qr = n(0)
        do i = 1, 3
            do j = 1, 3
                lhs = lhs + v(i)*(gl(i, j) + w*(gr(i, j) - gl(i, j)))*v(j)
                ql = ql + v(i)*gl(i, j)*v(j)
                qr = qr + v(i)*gr(i, j)*v(j)
            end do
        end do
        call ident('metric blend quadratic form is convex combination', &
            lhs - ((1 - w)*ql + w*qr))
    end subroutine blend_metric_pd

    ! ------------------------------ switch_der, 3ddc30b lines 97-122

    subroutine switch_chain
        type(expr_t) :: rho, k, width, t, phi, small, s0, sp, p1, p2, p3, d(0:3), e
        character(LN) :: o(4), w(4)
        integer :: br, m
        character(48) :: lab
        o = [character(LN) :: 's(0)', 's(1)', 's(2)', 's(3)']
        w = [character(LN) :: 'S0', 'S1', 'S2', 'S3']
        rho = s('rho'); k = s('k'); width = s('width')
        ! Line 97 parsed with real(k, dp) -> (k).
        t = src(fprov, 109, o, w)
        phi = subs_many(src(fprov, 116, o, w), [s('t')], [t])
        call ident('switch phi', phi - (-1/t + 1/(1 - t)))
        do br = 1, 2
            ! br=1: phi<0, small=exp(phi) (line 109 abs), s0 line 111.
            ! br=2: phi>=0, small=exp(-phi), s0 line 113.
            if (br == 1) then
                small = exp(phi)
                s0 = subs_many(src(fprov, 123, o, w), [s('small')], [small])
            else
                small = exp(-phi)
                s0 = subs_many(src(fprov, 125, o, w), [s('small')], [small])
            end if
            sp = subs_many(src(fprov, 127, o, w), [s('small')], [small])
            p1 = subs_many(src(fprov, 128, o, w), [s('t')], [t])
            p2 = subs_many(src(fprov, 129, o, w), [s('t')], [t])
            p3 = subs_many(src(fprov, 130, o, w), [s('t')], [t])
            call ident('switch ph1 = dphi/dt', p1 - dt(phi))
            call ident('switch ph2 = d2phi/dt2', p2 - dt(dt(phi)))
            call ident('switch ph3 = d3phi/dt3', p3 - dt(dt(dt(phi))))
            call ident('switch s0 = logistic(phi)', s0 - 1/(1 + exp(-phi)))
            d(0) = s0
            d(1) = subs_many(src(fprov, 131, o, w), &
                [s('sp'), s('ph1'), s('width')], [sp, p1, width])
            d(2) = subs_many(src(fprov, 132, o, w), &
                [s('sp'), s('ph1'), s('ph2'), s('S0'), s('width')], [sp, p1, p2, s0, width])
            d(3) = subs_many(src(fprov, 133, o, w), &
                [s('sp'), s('ph1'), s('ph2'), s('ph3'), s('S0'), s('width')], &
                [sp, p1, p2, p3, s0, width])
            do m = 1, 3
                write (lab, '(a,i1,a,i1)') 'switch branch ', br, ' s(m)=d s(m-1)/drho m=', m
                e = diff(d(m - 1), rho) - d(m)
                call ident(trim(lab), e)
            end do
        end do
    contains
        function dt(e0) result(r)
            type(expr_t), intent(in) :: e0
            type(expr_t) :: r
            ! d/dt at fixed width: t = (rho-k)/width + 1/2, so d/dt = width d/drho
            r = width*diff(e0, rho)
        end function dt
    end subroutine switch_chain

    ! ------------------------------ gauge projection, 3ddc30b lines 72-82

    subroutine gauge_projection
        type(expr_t) :: m, nn, dth, dze, wave(2), c, proj(2), th, ze, arg, chi
        type(expr_t) :: gc(2), gs(2), val(2), cc, cs
        m = s('m'); nn = s('n'); dth = s('dth'); dze = s('dze')
        call expect(fprov, 84, 'wave = [real(m, dp), -real(n, dp)]')
        call expect(fprov, 85, 'defect = maxval(abs(m*delta(2, :) + n*delta(1, :)))')
        call expect(fprov, 92, &
            'c = dot_product(wave, delta(:, r))/dot_product(wave, wave)')
        call expect(fprov, 93, 'delta(:, r) = c*wave')
        wave = [m, -nn]
        c = (wave(1)*dth + wave(2)*dze)/(wave(1)**2 + wave(2)**2)
        proj = [c*wave(1), c*wave(2)]
        ! Removed part is (defect/|wave|^2) (n, m): the tangential jump left
        ! after projection is bounded by |defect|/sqrt(m^2+n^2).
        call ident('gauge residual theta', (dth - proj(1)) - &
            (m*dze + nn*dth)/(m**2 + nn**2)*nn)
        call ident('gauge residual zeta', (dze - proj(2)) - &
            (m*dze + nn*dth)/(m**2 + nn**2)*m)
        ! A gauge term with coefficients c_r*wave (lines 195-196) is grad chi.
        th = s('theta'); ze = s('zeta')
        arg = m*th - nn*ze
        cc = s('cc'); cs = s('cs')
        gc = [cc*m, -cc*nn]; gs = [cs*m, -cs*nn]
        val = [gc(1)*cos(arg) + gs(1)*sin(arg), gc(2)*cos(arg) + gs(2)*sin(arg)]
        chi = cc*sin(arg) - cs*cos(arg)
        call ident('projected gauge theta is d chi/dtheta', val(1) - diff(chi, th))
        call ident('projected gauge zeta is d chi/dzeta', val(2) - diff(chi, ze))
    end subroutine gauge_projection

    ! ------------------------------ libneo spectral radial bases

    subroutine spectral_radial
        type(expr_t) :: sv, p, q, e, t0, t1, t2, r, mm, f(4), rm, rm1, rm2, z0, z1, z2
        type(expr_t) :: pz, qz
        character(LN) :: o(9), w(9)
        integer :: mi
        character(48) :: lab
        ! Chebyshev step lines 42-47 with T_{l-1}=P(s), T_{l-2}=Q(s).
        o = [character(LN) :: 'cheby(ll - 1, 0)', 'cheby(ll - 1, 1)', 'cheby(ll - 1, 2)', &
            'cheby(ll - 2, 0)', 'cheby(ll - 2, 1)', 'cheby(ll - 2, 2)', 'x', 'x', 'x']
        w = [character(LN) :: 'P0', 'P1', 'P2', 'Q0', 'Q1', 'Q2', 'x', 'x', 'x']
        sv = s('s')
        p = func('P', [sv]); q = func('Q', [sv])
        t0 = jr(src(fbasis, 42, o, w)); t1 = jr(src(fbasis, 43, o, w))
        t2 = jr(src(fbasis, 46, o, w))
        call ident('cheby step d col0 = col1', un(diff(t0, sv) - t1))
        call ident('cheby step d col1 = col2', un(diff(t1, sv) - t2))
        ! Corrections lines 52, 55: constant subtracted from col 0 only, all
        ! columns scaled by 1/(ll+1); base lines 39-40 ([1,0,0], [s,1,0]).
        call ident('cheby correction keeps jet', &
            diff(subs_many(src(fbasis, 52, ['cheby(ll, 0)'], ['T0']), [s('T0')], [p]), sv) &
            - diff(p, sv))
        ! Zernike recurrence lines 104-114, f1..f4 constants.
        o = [character(LN) :: 'zernike(n - 2, m, 0)', 'zernike(n - 2, m, 1)', &
            'zernike(n - 2, m, 2)', 'zernike(n - 4, m, 0)', 'zernike(n - 4, m, 1)', &
            'zernike(n - 4, m, 2)', 'x', 'x', 'x']
        r = s('r')
        pz = func('P', [r]); qz = func('Q', [r])
        z0 = jz(src(fbasis, 104, o, w)); z1 = jz(src(fbasis, 106, o, w))
        z2 = jz(src(fbasis, 110, o, w))
        call ident('zernike step d col0 = col1', unz(diff(z0, r) - z1))
        call ident('zernike step d col1 = col2', unz(diff(z1, r) - z2))
        ! Zernike base lines 83-94 for m = 0..8 with the code's rm,rm1,rm2
        ! (rm=r**m, rm1=r**(m-1) or 0 for m=0, rm2=r**(m-2) or 0 for m<2).
        do mi = 0, 8
            rm = r**mi
            rm1 = n(0); if (mi >= 1) rm1 = r**(mi - 1)
            rm2 = n(0); if (mi >= 2) rm2 = r**(mi - 2)
            f(1) = subs_many(src(fbasis, 83, ['x'], ['x']), [s('rm')], [rm])
            f(2) = subs_many(src(fbasis, 84, ['x'], ['x']), [s('rm1'), s('m')], [rm1, n(mi)])
            f(3) = subs_many(src(fbasis, 85, ['x'], ['x']), [s('rm2'), s('m')], [rm2, n(mi)])
            write (lab, '(a,i0)') 'zernike base (m,m) d0->d1 m=', mi
            call ident(trim(lab), diff(f(1), r) - f(2))
            write (lab, '(a,i0)') 'zernike base (m,m) d1->d2 m=', mi
            call ident(trim(lab), diff(f(2), r) - f(3))
            f(1) = subs_many(src(fbasis, 89, ['x'], ['x']), [s('rm'), s('m')], [rm, n(mi)])
            f(2) = subs_many(src(fbasis, 91, ['x'], ['x']), [s('rm'), s('rm1'), s('m')], &
                [rm, rm1, n(mi)])
            f(3) = subs_many(src(fbasis, 93, ['x'], ['x']), [s('rm'), s('rm2'), s('m')], &
                [rm, rm2, n(mi)])
            write (lab, '(a,i0)') 'zernike base (m+2,m) d0->d1 m=', mi
            call ident(trim(lab), diff(f(1), r) - f(2))
            write (lab, '(a,i0)') 'zernike base (m+2,m) d1->d2 m=', mi
            call ident(trim(lab), diff(f(2), r) - f(3))
        end do
        ! m=1 correction lines 128-133: subtract c*r from col0 and c from col1.
        e = diff(subs_many(src(fbasis, 128, ['zernike(n, 1, 0)'], ['Z0']), [s('Z0')], &
            [pz]), r) - subs_many(src(fbasis, 131, ['zernike(n, 1, 1)'], ['Z1']), &
            [s('Z1')], [diff(pz, r)])
        call ident('zernike m=1 correction keeps jet', e)
        ! Axis chain (eval lines 172, 187-188): sbar=(s+1)/2, tt=[Z,Z'/2,Z''/4].
        mm = func('Zf', [subs_many(src(fbasis, 172, ['x'], ['x']), [s('s')], [sv])])
        call ident('axis chain first', diff(mm, sv) - &
            subs_many(diff(func('Zf', [r]), r), [r], [(sv + 1)/2])/2)
        call ident('axis chain second', diff(diff(mm, sv), sv) - &
            subs_many(diff(diff(func('Zf', [r]), r), r), [r], [(sv + 1)/2])/4)
    contains
        function jr(e0) result(e1)
            type(expr_t), intent(in) :: e0
            type(expr_t) :: e1
            e1 = subs_many(e0, [s('P0'), s('P1'), s('P2'), s('Q0'), s('Q1'), s('Q2')], &
                [p, diff(p, sv), diff(diff(p, sv), sv), q, diff(q, sv), diff(diff(q, sv), sv)])
        end function jr
        function un(e0) result(e1)
            type(expr_t), intent(in) :: e0
            type(expr_t) :: e1
            e1 = jet_subs1(jet_subs1(e0, p, sv, 'P'), q, sv, 'Q')
        end function un
        function jz(e0) result(e1)
            type(expr_t), intent(in) :: e0
            type(expr_t) :: e1
            e1 = subs_many(e0, [s('P0'), s('P1'), s('P2'), s('Q0'), s('Q1'), s('Q2')], &
                [pz, diff(pz, r), diff(diff(pz, r), r), qz, diff(qz, r), &
                diff(diff(qz, r), r)])
        end function jz
        function unz(e0) result(e1)
            type(expr_t), intent(in) :: e0
            type(expr_t) :: e1
            e1 = jet_subs1(jet_subs1(e0, pz, r, 'P'), qz, r, 'Q')
        end function unz
    end subroutine spectral_radial

    function jet_subs1(e, f, x, prefix) result(r)
        type(expr_t), intent(in) :: e, f, x
        character(*), intent(in) :: prefix
        type(expr_t) :: r
        r = subs_many(e, [diff(diff(diff(f, x), x), x)], [s(prefix//'d3')])
        r = subs_many(r, [diff(diff(f, x), x)], [s(prefix//'d2')])
        r = subs_many(r, [diff(f, x)], [s(prefix//'d1')])
        r = subs_many(r, [f], [s(prefix//'d0')])
    end function jet_subs1

    ! ------------------------------ libneo Fourier mode accumulation

    !> Lines 283-306: (e cos(arg) + o sin(arg)) T(s), arg = m theta - n zeta.
    subroutine spectral_modes
        type(expr_t) :: sv, th, ze, m, nn, ce, co, arg, tf, mode, pp, qq, tt(0:2), q(3)
        character(LN) :: o(13), w(13)
        integer :: k, i, j
        integer, parameter :: d2i(6) = [1, 1, 1, 2, 2, 3], d2j(6) = [1, 2, 3, 2, 3, 3]
        character(40) :: lab
        logical :: have_pq
        ! Accumulators (val, d1(k), d2(k)) start from zero; by linearity the
        ! sum over modes and radial degrees is the jet of the sum.
        o = [character(LN) :: 'tt(0)', 'tt(1)', 'tt(2)', 'val', 'd1(1)', 'd1(2)', &
            'd1(3)', 'd2(1)', 'd2(2)', 'd2(3)', 'd2(4)', 'd2(5)', 'd2(6)']
        w = [character(LN) :: 'T0', 'T1', 'T2', '0', '0', '0', '0', '0', '0', '0', &
            '0', '0', '0']
        sv = s('s'); th = s('theta'); ze = s('zeta'); m = s('rmi'); nn = s('rni')
        ce = s('coeff_e'); co = s('coeff_o')
        q = [sv, th, ze]
        arg = m*th - nn*ze
        tf = func('T', [sv])
        tt(0) = tf; tt(1) = diff(tf, sv); tt(2) = diff(tt(1), sv)
        mode = (ce*cos(arg) + co*sin(arg))*tf
        call same_text(fbasis, 283, 300)
        call same_text(fbasis, 284, 301)
        have_pq = .false.
        pp = ins(src(fbasis, 300, o, w))
        qq = ins(src(fbasis, 301, o, w))
        have_pq = .true.
        call ident('mode value', ins(src(fbasis, 303, o, w)) - mode)
        do k = 1, 3
            write (lab, '(a,i1)') 'mode first ', k
            call ident(trim(lab), ins(src(fbasis, 303 + k, o, w)) - diff(mode, q(k)))
        end do
        do k = 1, 6
            i = d2i(k); j = d2j(k)
            write (lab, '(a,i1)') 'mode second slot ', k
            call ident(trim(lab), ins(src(fbasis, 284 + k, o, w)) - &
                diff(diff(mode, q(i)), q(j)))
        end do
    contains
        function ins(e0) result(e1)
            type(expr_t), intent(in) :: e0
            type(expr_t) :: e1
            e1 = subs_many(e0, [s('carg'), s('sarg'), s('T0'), s('T1'), s('T2')], &
                [cos(arg), sin(arg), tt(0), tt(1), tt(2)])
            if (have_pq) e1 = subs_many(e1, [s('p'), s('q')], [pp, qq])
        end function ins
    end subroutine spectral_modes

    ! ------------------------------ volume chain, 3ddc30b lines 173-188

    !> s = 2 (rho - owner + 1) - 1 (line 174): first rho derivative x2 (line
    !> 179), second derivatives x 2**(#rho indices) (line 186), slots from
    !> line 171.
    subroutine spectral_volume_chain
        type(expr_t) :: rho, th, ze, ow, sarg, g, q(3), u(3), gg, e
        integer :: i, j, slot(3, 3), nr
        character(32) :: lab
        call expect(fprov, 196, &
            'integer, parameter :: slot(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])')
        call expect(fprov, 204, 'da(:, 1) = 2.0_dp*da(:, 1)')
     call expect(fprov, 212, 'real(2**(merge(1, 0, i == 1) + merge(1, 0, j == 1)), dp)')
        slot = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
        rho = s('rho'); th = s('theta'); ze = s('zeta'); ow = s('owner')
        sarg = parse_rhs(ar, '2.0_dp*(q(1) - owner + 1) - 1.0_dp', ['q(1)'], ['rho '])
  call expect(fprov, 199, '2.0_dp*(q(1) - owner + 1) - 1.0_dp, q(2), q(3), av, second)')
        q = [rho, th, ze]
        u = [s('u'), s('v'), s('w')]
        g = func('G', [sarg, th, ze])
        gg = func('G', u)
        do i = 1, 3
            e = diff(g, q(i)) - merge(2, 1, i == 1)*subs_many(diff(gg, u(i)), u, [sarg, th, ze])
            write (lab, '(a,i1)') 'volume chain first ', i
            call ident(trim(lab), e)
        end do
        do j = 1, 3
            do i = 1, 3
                nr = merge(1, 0, i == 1) + merge(1, 0, j == 1)
                ! slot(i,j) indexes (ss,st,sz,tt,tz,zz): check it names (i,j).
                if (sl(slot(i, j)) /= 10*min(i, j) + max(i, j)) error stop 'slot table'
                e = diff(diff(g, q(i)), q(j)) - &
                    2**nr*subs_many(diff(diff(gg, u(i)), u(j)), u, [sarg, th, ze])
                write (lab, '(a,2i1)') 'volume chain second ', i, j
                call ident(trim(lab), e)
            end do
        end do
    contains
        integer function sl(k)
            integer, intent(in) :: k
            integer, parameter :: code(6) = [11, 12, 13, 22, 23, 33]
            sl = code(k)
        end function sl
    end subroutine spectral_volume_chain

    ! ------------------------------ geometry Fourier jets, 3ddc30b 273-299

    subroutine geometry_fourier
        type(expr_t) :: rho, th, ze, m, nn, arg, c1, c2, x, cc(0:3), wave(2), e, got
        type(expr_t) :: k1(0:3), k2(0:3)
        integer :: r, t, z, i
        character(32) :: lab
        call expect(fprov, 321, 'wave = [real(m, dp), -real(n, dp)]')
        call expect(fprov, 323, 'c = [cos(arg), -sin(arg), -cos(arg), sin(arg)]')
        call expect(fprov, 44, &
            'angular_factor(t, z, j) = real(spectre_field%data%im(j), dp)**t* &')
        call expect(fprov, 341, 'factor = angular_factor(t, z, ii)')
        call expect(fprov, 342, 'jets(r, t, z, comp) = jets(r, t, z, comp) + factor* &')
        call expect(fprov, 343, '(coeff(1, r)*c(mod(t + z, 4)) &')
        call expect(fprov, 344, '- coeff(2, r)*c(mod(t + z + 1, 4)))')
        rho = s('rho'); th = s('theta'); ze = s('zeta'); m = s('m'); nn = s('n')
        wave = [m, -nn]
        arg = wave(1)*th + wave(2)*ze
        cc = [cos(arg), -sin(arg), -cos(arg), sin(arg)]
        c1 = func('C1', [rho]); c2 = func('C2', [rho])
        k1(0) = c1; k2(0) = c2
        do i = 1, 3
            k1(i) = diff(k1(i - 1), rho); k2(i) = diff(k2(i - 1), rho)
        end do
        x = c1*cos(arg) + c2*sin(arg)
        do z = 0, 3
            do t = 0, 3 - z
                do r = 0, 3 - z - t
                    got = wave(1)**t*wave(2)**z*(k1(r)*cc(mod(t + z, 4)) - &
                        k2(r)*cc(mod(t + z + 1, 4)))
                    e = x
                    do i = 1, r
                        e = diff(e, rho)
                    end do
                    do i = 1, t
                        e = diff(e, th)
                    end do
                    do i = 1, z
                        e = diff(e, ze)
                    end do
                    write (lab, '(a,3i1)') 'geometry Fourier jet ', r, t, z
                    call ident(trim(lab), e - got)
                end do
            end do
        end do
    end subroutine geometry_fourier

end module provider_jets_impl

program gen_provider_jets
    use provider_jets_impl, only: run_all
    implicit none
    call run_all
end program gen_provider_jets
