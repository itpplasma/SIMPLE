program generate_spectre_kernels
    ! Pointwise kernels of the SPECTRE full-orbit path: the cylindrical chart
    ! (R, zeta, Z)(rho, theta, zeta), B = curl A with A_rho = 0, and the
    ! guiding-centre <-> particle maps at zeroth order in the Larmor radius,
    ! and grad A_zeta in axisymmetry for the canonical-momentum shift.
    use fortsym_expr, only: expr_t, sym, num, sin, cos, sqrt
    use fortsym_expr, only: operator(+), operator(-), operator(*), &
                            operator(/), operator(**)
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_string, only: str, chars, str_t, words
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(symengine_engine_t) :: engine
    type(expr_t) :: zeta, r, z, dr(3), dz(3), e(3, 3), ec(3, 3), cof(3, 3), jac
    type(expr_t) :: dat(3), daz(3), sqgb(3), bvec(3), bmod
    type(expr_t) :: xg(3), b(3), e1(3), e2(3), bm, vpar, vperp, alpha, ro0
    type(expr_t) :: eperp(3), vp(3), x(3), v(3), rl(3), vpar_out, vperp2
    type(expr_t) :: geometry(22), field(4), lift(6), projection(5)
    type(expr_t) :: d1(3), d2(3), flux(3)
    character(4096) :: output_dir
    integer :: i, j, unit, ios

    output_dir = 'src/generated'
    if (command_argument_count() > 0) call get_command_argument(1, output_dir)
    engine = make_symengine_engine()

    ! Chart x = (R cos zeta, R sin zeta, Z) with R, Z functions of (rho, theta,
    ! zeta); dr(i), dz(i) are their partial derivatives.
    zeta = sym('zeta')
    r = sym('r')
    z = sym('z')
    do i = 1, 3
        dr(i) = sym('dr'//digit(i))
        dz(i) = sym('dz'//digit(i))
    end do
    ! Covariant basis in the orthonormal frame (e_R, e_phi, e_Z), which is
    ! right-handed like (e_x, e_y, e_z), so the Jacobian and the dual basis are
    ! free of zeta; only the final rotation to Cartesian components uses it.
    do i = 1, 3
        ec(1, i) = dr(i)
        ec(2, i) = num(0)
        ec(3, i) = dz(i)
    end do
    ec(2, 3) = r
    cof(:, 1) = cross(ec(:, 2), ec(:, 3))
    cof(:, 2) = cross(ec(:, 3), ec(:, 1))
    cof(:, 3) = cross(ec(:, 1), ec(:, 2))
    jac = dot(ec(:, 1), cof(:, 1))
    geometry(1) = r*cos(zeta)
    geometry(2) = r*sin(zeta)
    geometry(3) = z
    geometry(13) = jac
    do j = 1, 3
        geometry(3 + 3*(j - 1) + 1:3 + 3*j) = to_cartesian(ec(:, j), zeta)
        do i = 1, 3
            cof(i, j) = cof(i, j)/jac
        end do
        geometry(13 + 3*(j - 1) + 1:13 + 3*j) = to_cartesian(cof(:, j), zeta)
    end do

    ! |J| B^i = eps^ijk d_j A_k with A_rho = 0, so B = e_i B^i.
    do i = 1, 3
        dat(i) = sym('dat'//digit(i))
        daz(i) = sym('daz'//digit(i))
    end do
    sqgb(1) = daz(2) - dat(3)
    sqgb(2) = -daz(1)
    sqgb(3) = dat(1)
    jac = sym('jac')
    do i = 1, 3
        do j = 1, 3
            e(i, j) = sym('e'//digit(i)//digit(j))
        end do
    end do
    do i = 1, 3
        bvec(i) = (e(i, 1)*sqgb(1) + e(i, 2)*sqgb(2) + e(i, 3)*sqgb(3))/abs_sym(jac)
    end do
    bmod = sqrt(dot(bvec, bvec))
    do i = 1, 3
        field(i) = bvec(i)/bmod
    end do
    field(4) = bmod

    ! Lift: v = vpar b + vperp (cos alpha e1 + sin alpha e2),
    ! x = X + (ro0/B) b x v_perp.
    do i = 1, 3
        xg(i) = sym('xg'//digit(i))
        b(i) = sym('b'//digit(i))
        e1(i) = sym('e1'//digit(i))
        e2(i) = sym('e2'//digit(i))
        x(i) = sym('x'//digit(i))
        v(i) = sym('v'//digit(i))
    end do
    bm = sym('bmod')
    vpar = sym('vpar')
    vperp = sym('vperp')
    alpha = sym('alpha')
    ro0 = sym('ro0')
    do i = 1, 3
        eperp(i) = cos(alpha)*e1(i) + sin(alpha)*e2(i)
        vp(i) = vperp*eperp(i)
    end do
    rl = cross(b, vp)
    do i = 1, 3
        lift(i) = xg(i) + (ro0/bm)*rl(i)
        lift(3 + i) = vpar*b(i) + vp(i)
    end do

    ! Projection: X = x - (ro0/B) b x v_perp, v_par = v.b, |v_perp|^2.
    vpar_out = dot(v, b)
    do i = 1, 3
        vp(i) = v(i) - vpar_out*b(i)
    end do
    rl = cross(b, vp)
    do i = 1, 3
        projection(i) = x(i) - (ro0/bm)*rl(i)
    end do
    projection(4) = vpar_out
    projection(5) = dot(v, v) - vpar_out**2

    ! Without zeta dependence, d_rho A_zeta = -|J| B^theta and
    ! d_theta A_zeta = |J| B^rho, so grad A_zeta follows from the local B and
    ! the dual basis d_i = grad u^i, with B^i = |B| b.d_i.
    do i = 1, 3
        d1(i) = sym('d'//digit(i)//'1')
        d2(i) = sym('d'//digit(i)//'2')
    end do
    do i = 1, 3
        flux(i) = abs_sym(jac)*bm*(dot(b, d1)*d2(i) - dot(b, d2)*d1(i))
    end do

    open (newunit=unit, file=trim(output_dir)//'/spectre_kernels_generated.f90', &
          status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Cannot open generated SPECTRE kernel module'
    write (unit, '(a)') 'module spectre_kernels_generated'
    write (unit, '(a)') 'implicit none'
    write (unit, '(a)') 'contains'
    call emit('spectre_chart_geometry', geometry, &
              words('zeta r dr1 dr2 dr3 z dz1 dz2 dz3'), &
              words('x1 x2 x3 e11 e21 e31 e12 e22 e32 e13 e23 e33 jac '// &
                    'd11 d21 d31 d12 d22 d32 d13 d23 d33'))
    call emit('spectre_curl_field', field, &
              words('e11 e21 e31 e12 e22 e32 e13 e23 e33 jac '// &
                    'dat1 dat3 daz1 daz2'), words('b1 b2 b3 bmod'))
    call emit('spectre_larmor_lift', lift, &
              words('xg1 xg2 xg3 b1 b2 b3 e11 e12 e13 e21 e22 e23 '// &
                    'bmod vpar vperp alpha ro0'), words('x1 x2 x3 v1 v2 v3'))
    call emit('spectre_gc_projection', projection, &
              words('x1 x2 x3 v1 v2 v3 b1 b2 b3 bmod ro0'), &
              words('xg1 xg2 xg3 vpar vperp2'))
    call emit('spectre_flux_gradient', flux, &
              words('b1 b2 b3 bmod jac d11 d21 d31 d12 d22 d32'), &
              words('g1 g2 g3'))
    write (unit, '(a)') 'end module spectre_kernels_generated'
    close (unit)
contains
    function cross(a, c) result(d)
        type(expr_t), intent(in) :: a(3), c(3)
        type(expr_t) :: d(3)
        d(1) = a(2)*c(3) - a(3)*c(2)
        d(2) = a(3)*c(1) - a(1)*c(3)
        d(3) = a(1)*c(2) - a(2)*c(1)
    end function cross

    function to_cartesian(a, phi) result(d)
        ! (a_R, a_phi, a_Z) -> (a_x, a_y, a_z) at toroidal angle phi.
        type(expr_t), intent(in) :: a(3), phi
        type(expr_t) :: d(3)
        d(1) = a(1)*cos(phi) - a(2)*sin(phi)
        d(2) = a(1)*sin(phi) + a(2)*cos(phi)
        d(3) = a(3)
    end function to_cartesian

    function dot(a, c) result(d)
        type(expr_t), intent(in) :: a(:), c(:)
        type(expr_t) :: d
        integer :: k
        d = a(1)*c(1)
        do k = 2, size(a)
            d = d + a(k)*c(k)
        end do
    end function dot

    function abs_sym(a) result(d)
        ! |J| as sqrt(J**2): the generated code needs no sign branch.
        type(expr_t), intent(in) :: a
        type(expr_t) :: d
        d = sqrt(a*a)
    end function abs_sym

    character(1) function digit(k)
        integer, intent(in) :: k
        digit = achar(iachar('0') + k)
    end function digit

    subroutine emit(name, roots, args, outputs)
        character(*), intent(in) :: name
        type(expr_t), intent(inout) :: roots(:)
        type(str_t), intent(in) :: args(:), outputs(:)
        type(kernel_spec_t) :: spec
        type(engine_result_t) :: result
        type(str_t) :: code
        character(:), allocatable :: message
        logical :: ok
        integer :: k
        do k = 1, size(roots)
            result = engine%simplify(roots(k))
            if (.not. result%ok) error stop 'FortSym simplification failed'
            roots(k) = result%value
        end do
        spec%name = str(name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str('gen/generate_spectre_kernels.f90')
        spec%pure_procedure = .true.
        spec%args = args
        spec%outputs = outputs
        code = emit_kernel(roots, spec, ok=ok, message=message)
        if (.not. ok) then
            write (*, '(a)') name//': '//message
            error stop 'FortSym kernel emission failed'
        end if
        write (unit, '(a)') chars(code)
    end subroutine emit
end program generate_spectre_kernels
