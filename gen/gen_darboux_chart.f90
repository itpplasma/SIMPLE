!> Pointwise mathematics of the global Darboux chart for guiding-centre
!> orbits through the magnetic axis.
!>
!> Field point u = (x, y) = sqrt(s) (cos theta, sin theta), toroidal angle
!> phi = p and canonical momentum p_phi = pp. With
!>     v_par = (p_phi - A_phi/ro0)/h_phi,   H = v_par**2/2 + mu B,
!>     p_theta = h_theta v_par + A_theta/ro0 = s (q_h v_par + q_A/ro0),
!> where q_A = A_theta/s and q_h = h_theta/s are regular, the chart
!>     kappa = sigma (q_h v_par + q_A/ro0) = J/s > 0,   J = sigma p_theta,
!>     X = sqrt(2 kappa) x,   Y = sigma sqrt(2 kappa) y
!> has X**2 + Y**2 = 2 J and dX ^ dY = dp_theta ^ dtheta, so (X, Y, phi, p_phi)
!> is one Darboux chart with the axis at X = Y = 0. The canonical vector
!> field is
!>     dX/dt = -dK/dY, dY/dt = dK/dX, dphi/dt = dK/dp_phi, dp_phi/dt = -dK/dphi
!> for K(w) = H(u(w), phi, p_phi), where u(w) solves the chart map
!> C(u, phi, p_phi) = (X, Y). Implicit differentiation gives
!>     du/dX = G^-1 e_1, du/dY = G^-1 e_2, du/dphi = -G^-1 C_phi,
!>     du/dp_phi = -G^-1 C_p_phi,   G = dC/du,
!> and every factor is regular at the axis. The generated leaves are
!>   darboux_quotient_d1/d2   q = h/s with jets, in the nested form
!>                            q_x = (h_x - 2 x q)/s that never forms s**2;
!>   darboux_point            kappa, p_theta, H, v_par, C, the vector field
!>                            F(u, phi, p_phi), G and du/dt (first jets);
!>   darboux_stage            C, F and their derivatives with respect to
!>                            u and (phi, p_phi) for the Newton Jacobian of
!>                            the implicit Runge-Kutta stage system
!>                            (second jets);
!>   test_field_regular       the analytic circular test tokamak of
!>                            field_can_test in (x, y) with second jets
!>                            (second output file).
!> Jet order: value, d/dx, d/dy, d/dphi, xx, xy, xphi, yy, yphi, phiphi.
program gen_darboux_chart
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, func, sqrt, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_binding_t
    use simple_gen_support, only: leaf_writer_t, fortsym_revision, binding, &
        output_t, out_scalar, out_array
    implicit none

    character(2), parameter :: heads(5) = [character(2) :: 'QA', 'QH', 'AP', &
        'HP', 'BM']
    character(2), parameter :: names(5) = [character(2) :: 'qa', 'qh', 'ap', &
        'hp', 'bm']
    type(leaf_writer_t) :: w
    type(arena_t), target :: arena
    type(expr_t) :: x, y, p, pp, ro0, mu, sg
    type(expr_t) :: qa, qh, ap, hp, bm, vpar, ham, kap, rt, chart(2)
    type(expr_t) :: g(2, 2), det, ginv(2, 2), cphi(2), cpp(2)
    type(expr_t) :: du(2, 4), dh(4), dk(4), rhs(4), udot(2)
    type(kernel_binding_t), allocatable :: bind1(:), bind2(:)
    character(len=256) :: path, test_path
    integer :: k

    path = 'src/generated/darboux_chart_leaves.inc'
    test_path = 'src/generated/test_field_regular_leaves.inc'
    if (command_argument_count() >= 1) call get_command_argument(1, path)
    if (command_argument_count() >= 2) call get_command_argument(2, test_path)

    call arena%init()
    x = sym(arena, 'x')
    y = sym(arena, 'y')
    p = sym(arena, 'phi')
    pp = sym(arena, 'pphi')
    ro0 = sym(arena, 'ro0')
    mu = sym(arena, 'mu')
    sg = sym(arena, 'sigma')

    allocate (bind1(0), bind2(0))
    do k = 1, size(heads)
        bind1 = [bind1, jet_bindings(heads(k), names(k), 4)]
        bind2 = [bind2, jet_bindings(heads(k), names(k), 10)]
    end do
    qa = func(heads(1), [x, y, p])
    qh = func(heads(2), [x, y, p])
    ap = func(heads(3), [x, y, p])
    hp = func(heads(4), [x, y, p])
    bm = func(heads(5), [x, y, p])

    vpar = (pp - ap/ro0)/hp
    ham = vpar**2/2 + mu*bm
    kap = sg*(qh*vpar + qa/ro0)
    rt = sqrt(2*kap)
    chart = [rt*x, sg*rt*y]

    do k = 1, 2
        g(k, 1) = diff(chart(k), x)
        g(k, 2) = diff(chart(k), y)
        cphi(k) = diff(chart(k), p)
        cpp(k) = diff(chart(k), pp)
    end do
    det = g(1, 1)*g(2, 2) - g(1, 2)*g(2, 1)
    ginv(1, 1) = g(2, 2)/det
    ginv(1, 2) = -g(1, 2)/det
    ginv(2, 1) = -g(2, 1)/det
    ginv(2, 2) = g(1, 1)/det
    do k = 1, 2
        du(k, 1) = ginv(k, 1)
        du(k, 2) = ginv(k, 2)
        du(k, 3) = -(ginv(k, 1)*cphi(1) + ginv(k, 2)*cphi(2))
        du(k, 4) = -(ginv(k, 1)*cpp(1) + ginv(k, 2)*cpp(2))
    end do
    dh = [diff(ham, x), diff(ham, y), diff(ham, p), diff(ham, pp)]
    dk(1) = dh(1)*du(1, 1) + dh(2)*du(2, 1)
    dk(2) = dh(1)*du(1, 2) + dh(2)*du(2, 2)
    dk(3) = dh(3) + dh(1)*du(1, 3) + dh(2)*du(2, 3)
    dk(4) = dh(4) + dh(1)*du(1, 4) + dh(2)*du(2, 4)
    rhs = [-dk(2), dk(1), dk(4), -dk(3)]
    do k = 1, 2
        udot(k) = du(k, 1)*rhs(1) + du(k, 2)*rhs(2) + du(k, 3)*rhs(3) + &
                  du(k, 4)*rhs(4)
    end do

    call w%open(arena, 'gen_darboux_chart', trim(path))
    w%sign = sg
    w%has_sign = .true.
    call w%line('! Generated by gen/gen_darboux_chart.f90; do not edit.')
    call w%line('! FortSym revision: '//fortsym_revision)
    call w%line('! Global Darboux chart leaves, included by')
    call w%line('! src/orbit_symplectic_darboux.f90. Jet order: value, d/dx, d/dy,')
    call w%line('! d/dphi, xx, xy, xphi, yy, yphi, phiphi.')
    call emit_quotients()
    call emit_point()
    call emit_stage()
    call w%close()

    call w%open(arena, 'gen_darboux_chart', trim(test_path))
    w%has_sign = .false.
    call w%line('! Generated by gen/gen_darboux_chart.f90; do not edit.')
    call w%line('! FortSym revision: '//fortsym_revision)
    call w%line('! Analytic circular test tokamak in (x, y) = r (cos, sin) theta,')
    call w%line('! included by src/field/field_can_test.f90. Jet order: value,')
    call w%line('! d/dx, d/dy, d/dphi, xx, xy, xphi, yy, yphi, phiphi.')
    call emit_test_field()
    call w%close()

contains

    !> Bindings of an opaque field function and its partial derivatives in
    !> (x, y, phi) to the elements of a jet array.
    function jet_bindings(head, name, n) result(b)
        character(*), intent(in) :: head, name
        integer, intent(in) :: n
        type(kernel_binding_t), allocatable :: b(:)
        integer, parameter :: idx(2, 6) = reshape([1, 1, 1, 2, 1, 3, 2, 2, &
            2, 3, 3, 3], [2, 6])
        integer :: j

        allocate (b(n))
        b(1) = binding(head, [integer ::], ref(name, 1))
        do j = 1, 3
            b(1 + j) = binding(head, [j], ref(name, 1 + j))
        end do
        if (n < 10) return
        do j = 1, 6
            b(4 + j) = binding(head, idx(:, j), ref(name, 4 + j))
        end do
    end function jet_bindings

    pure function ref(name, j) result(r)
        character(*), intent(in) :: name
        integer, intent(in) :: j
        character(:), allocatable :: r
        character(8) :: buf

        write (buf, '(i0)') j
        r = trim(name)//'('//trim(buf)//')'
    end function ref

    !> Jets of an expression in (x, y, phi): value, gradient, Hessian.
    function jet(e, n) result(r)
        type(expr_t), intent(in) :: e
        integer, intent(in) :: n
        type(expr_t) :: r(n)
        type(expr_t) :: ex, ey, ep

        ex = diff(e, x)
        ey = diff(e, y)
        ep = diff(e, p)
        r(1:4) = [e, ex, ey, ep]
        if (n < 10) return
        r(5:10) = [diff(ex, x), diff(ex, y), diff(ex, p), diff(ey, y), &
                   diff(ey, p), diff(ep, p)]
    end function jet

    subroutine emit_quotients()
        type(expr_t) :: h, s, q(10)
        type(kernel_binding_t), allocatable :: hb(:)

        h = func('HQ', [x, y, p])
        hb = jet_bindings('HQ', 'h', 10)
        s = x**2 + y**2
        q(1) = h/s
        q(2) = (diff(h, x) - 2*x*q(1))/s
        q(3) = (diff(h, y) - 2*y*q(1))/s
        q(4) = diff(h, p)/s
        q(5) = (diff(diff(h, x), x) - 2*q(1) - 4*x*q(2))/s
        q(6) = (diff(diff(h, x), y) - 2*x*q(3) - 2*y*q(2))/s
        q(7) = (diff(diff(h, x), p) - 2*x*q(4))/s
        q(8) = (diff(diff(h, y), y) - 2*q(1) - 4*y*q(3))/s
        q(9) = (diff(diff(h, y), p) - 2*y*q(4))/s
        q(10) = diff(diff(h, p), p)/s
        call w%emit('darboux_quotient_d1', [character(8) :: 'x', 'y', 'h(4)'], &
                    [out_array('q', 4)], q(1:4), hb(1:4))
        call w%emit('darboux_quotient_d2', [character(8) :: 'x', 'y', 'h(10)'], &
                    [out_array('q', 10)], q, hb)
    end subroutine emit_quotients

    subroutine emit_point()
        character(8), parameter :: args(11) = [character(8) :: 'x', 'y', &
            'pphi', 'ro0', 'mu', 'sigma', 'qa(4)', 'qh(4)', 'ap(4)', 'hp(4)', &
            'bm(4)']
        type(expr_t) :: roots(16)

        roots = [kap, (x**2 + y**2)*(qh*vpar + qa/ro0), ham, vpar, chart, rhs, &
                 g(1, 1), g(2, 1), g(1, 2), g(2, 2), udot]
        call w%emit('darboux_point', args, &
                    [out_scalar('kappa'), out_scalar('pth'), out_scalar('ham'), &
                     out_scalar('vpar'), out_array('chart', 2), out_array('rhs', 4), &
                     out_array('gu', 2, 2), out_array('udot', 2)], roots, bind1)
    end subroutine emit_point

    subroutine emit_stage()
        character(8), parameter :: args(11) = [character(8) :: 'x', 'y', &
            'pphi', 'ro0', 'mu', 'sigma', 'qa(10)', 'qh(10)', 'ap(10)', &
            'hp(10)', 'bm(10)']
        type(expr_t) :: roots(30)
        integer :: i

        roots(1:2) = chart
        roots(3:6) = rhs
        do i = 1, 4
            roots(6 + i) = diff(rhs(i), x)
            roots(10 + i) = diff(rhs(i), y)
            roots(14 + i) = diff(rhs(i), p)
            roots(18 + i) = diff(rhs(i), pp)
        end do
        roots(23:26) = [g(1, 1), g(2, 1), g(1, 2), g(2, 2)]
        roots(27:30) = [cphi(1), cphi(2), cpp(1), cpp(2)]
        call w%emit('darboux_stage', args, &
                    [out_array('chart', 2), out_array('rhs', 4), &
                     out_array('fu', 4, 2), out_array('fw', 4, 2), &
                     out_array('gu', 2, 2), out_array('gw', 2, 2)], roots, bind2)
    end subroutine emit_stage

    !> eval_field_test with B0 = iota0 = R0 = 1 and a = 1/2: in
    !> (x, y) = r (cos, sin) theta every component is a polynomial.
    subroutine emit_test_field()
        type(expr_t) :: s, comp(5), roots(50)
        integer :: k

        s = x**2 + y**2
        comp(1) = s/2 - s*x/3
        comp(2) = -(s/2 - s**2)
        comp(3) = (1 - 4*s)*s
        comp(4) = 1 + x
        comp(5) = 1 - x
        do k = 1, 5
            roots(10*k - 9:10*k) = jet(comp(k), 10)
        end do
        call w%emit('test_field_regular', [character(8) :: 'x', 'y'], &
                    [out_array('ath', 10), out_array('aph', 10), &
                     out_array('hth', 10), out_array('hph', 10), &
                     out_array('bmod', 10)], roots)
    end subroutine emit_test_field

end program gen_darboux_chart
