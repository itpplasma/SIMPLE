program gen_fo_hessian
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, func, operator(+), operator(-), &
                            operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs_many
    use fortsym_fparse, only: parse_fortran_expr
    use fortsym_string, only: str_t, str, chars
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE, &
                              cse_result_t, cse_analyse
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    implicit none
    type(arena_t), target :: arena
    type(symengine_engine_t) :: eng
    type(suite_t) :: suite
    type(kernel_spec_t) :: spec
    type(expr_t) :: q(3), p(3), aa(3), g(3, 3), gi(3, 3), w(3), v(3), b(3, 3)
    type(expr_t) :: mass, charge, ham, hh(6, 6), roots(36), term, ref
   type(expr_t) :: proof_hh(6, 6), proof_a(3), proof_g(3, 3), proof_iv(3, 3), proof_v(3)
    type(expr_t) :: before(99), after(99)
    type(expr_t) :: tjet(10), gradient(6), energy, e0, kj(6, 6)
    type(expr_t) :: aj(3, 10), gj(6, 10), iv(3, 3), vv(3), bb(3, 3)
    type(str_t), allocatable :: args(:), outputs(:)
    character(1024) :: out
    character(32) :: name
    integer :: i, j, k, l, ii, u, status, n, c, slot(3, 3), index(3, 10), nsub
    integer, parameter :: si(6) = [1, 1, 1, 2, 2, 3], sj(6) = [1, 2, 3, 2, 3, 3]

    call arena%init()
    eng = make_symengine_engine(arena)
    call get_environment_variable('FO_GENERATOR_OUTPUT', out, status=status)
    if (status /= 0) out = 'concept/fo_symplectic_proof/generated'
    mass = sym(arena, 'mass')
    charge = sym(arena, 'charge')
    do i = 1, 3
        write (name, '(a,i0)') 'q', i
        q(i) = sym(arena, trim(name))
        write (name, '(a,i0)') 'p', i
        p(i) = sym(arena, trim(name))
    end do
    do i = 1, 3
        write (name, '(a,i0)') 'A', i
        aa(i) = func(trim(name), q)
        w(i) = p(i) - charge*aa(i)
    end do
    do l = 1, 6
        write (name, '(a,i0)') 'G', l
        g(si(l), sj(l)) = func(trim(name), q)
        g(sj(l), si(l)) = g(si(l), sj(l))
    end do
    call invert(g, gi)
    do i = 1, 3
        v(i) = num(arena, 0)
        do j = 1, 3
            v(i) = v(i) + gi(i, j)*w(j)
        end do
    end do
    ham = num(arena, 0)
    do i = 1, 3
        ham = ham + w(i)*v(i)/(2*mass)
        do j = 1, 3
            b(i, j) = charge*diff(aa(i), q(j))
            do k = 1, 3
                b(i, j) = b(i, j) + diff(g(i, k), q(j))*v(k)
            end do
        end do
    end do
    call blocks(aa, g, gi, v, b, hh)
    call suite_begin(suite, 'Factored FO Hessian exact identities')
    do j = 1, 6
        do i = 1, 6
            if (i <= 3) then
                ref = diff(ham, q(i))
            else
                ref = diff(ham, p(i - 3))
            end if
            if (j <= 3) then
                ref = diff(ref, q(j))
            else
                ref = diff(ref, p(j - 3))
            end if
            write (name, '(a,i0,a,i0)') 'Hessian ', i, ',', j
            call check_identity(suite, eng, trim(name), ref - hh(i, j))
        end do
    end do
    proof_hh = hh
    proof_a = aa
    proof_g = g
    proof_iv = gi
    proof_v = v
    index = reshape([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, &
                     1, 1, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 2], [3, 10])
    allocate (args(200))
    n = 0
    do i = 1, 3
        write (name, '(a,i0)') 'p', i
        call input(name, p(i))
    end do
    call input('mass', mass)
    call input('charge', charge)
    do j = 1, 10
        do i = 1, 3
            write (name, '(a,i0,a,i0)') 'a', i, '_', j
            call input(name, aj(i, j))
        end do
        do l = 1, 6
            write (name, '(a,i0,a,i0)') 'g', l, '_', j
            call input(name, gj(l, j))
        end do
    end do
    do l = 1, 6
        g(si(l), sj(l)) = gj(l, 1)
        g(sj(l), si(l)) = gj(l, 1)
    end do
    call invert(g, gi)
    do i = 1, 3
        v(i) = num(arena, 0)
        do j = 1, 3
            v(i) = v(i) + gi(i, j)*(p(j) - charge*aj(j, 1))
        end do
    end do
    allocate (outputs(9))
    do l = 1, 6
        write (name, '(a,i0)') 'inverse', l
        outputs(l) = str(trim(name))
        roots(l) = gi(si(l), sj(l))
    end do
    do i = 1, 3
        write (name, '(a,i0)') 'velocity', i
        outputs(6 + i) = str(trim(name))
        roots(6 + i) = v(i)
    end do
    call leaf('fo_factor_velocity', roots(:9), args(:n), outputs)
    do l = 1, 6
        write (name, '(a,i0)') 'inverse', l
        call input(name, iv(si(l), sj(l)))
        iv(sj(l), si(l)) = iv(si(l), sj(l))
    end do
    do i = 1, 3
        write (name, '(a,i0)') 'velocity', i
        call input(name, vv(i))
    end do
    do i = 1, 3
        do j = 1, 3
            bb(i, j) = charge*aj(i, 1 + j)
            do k = 1, 3
                do l = 1, 6
                    if (si(l) == min(i, k) .and. sj(l) == max(i, k)) exit
                end do
                bb(i, j) = bb(i, j) + gj(l, 1 + j)*vv(k)
            end do
        end do
    end do
    slot = reshape([5, 6, 8, 6, 7, 9, 8, 9, 10], [3, 3])
    do j = 1, 3
        do i = 1, 3
            term = num(arena, 0)
            do k = 1, 3
                do l = 1, 3
                    term = term + bb(k, i)*iv(k, l)*bb(l, j)
                end do
                term = term - charge*aj(k, slot(i, j))*vv(k)
            end do
            do l = 1, 6
                c = 1
                if (si(l) /= sj(l)) c = 2
                term = term - c*vv(si(l))*gj(l, slot(i, j))*vv(sj(l))/2
            end do
            hh(i, j) = term/mass
            term = num(arena, 0)
            do k = 1, 3
                term = term - iv(i, k)*bb(k, j)/mass
            end do
            hh(3 + i, j) = term
            hh(j, 3 + i) = term
            hh(3 + i, 3 + j) = iv(i, j)/mass
        end do
    end do
    deallocate (outputs)
    allocate (outputs(36))
    do j = 1, 6
        do i = 1, 6
            k = 6*(j - 1) + i
            roots(k) = hh(i, j)
            write (name, '(a,i0,a,i0)') 'hh', i, '_', j
            outputs(k) = str(trim(name))
        end do
    end do
    nsub = 0
    do j = 1, 10
        do i = 1, 3
            term = proof_a(i)
            do k = 1, 3
                do c = 1, index(k, j)
                    term = diff(term, q(k))
                end do
            end do
            nsub = nsub + 1
            before(nsub) = aj(i, j)
            after(nsub) = term
        end do
        do l = 1, 6
            term = proof_g(si(l), sj(l))
            do k = 1, 3
                do c = 1, index(k, j)
                    term = diff(term, q(k))
                end do
            end do
            nsub = nsub + 1
            before(nsub) = gj(l, j)
            after(nsub) = term
        end do
    end do
    do l = 1, 6
        nsub = nsub + 1
        before(nsub) = iv(si(l), sj(l))
        after(nsub) = proof_iv(si(l), sj(l))
    end do
    do i = 1, 3
        nsub = nsub + 1
        before(nsub) = vv(i)
        after(nsub) = proof_v(i)
    end do
    do j = 1, 6
        do i = 1, 6
            call check_identity(suite, eng, 'jet binding matches derived Hessian', &
                      subs_many(hh(i, j), before(:nsub), after(:nsub)) - proof_hh(i, j))
        end do
    end do
    call leaf('fo_factor_hessian', roots, args(:n), outputs)
    n = 0
    call input('energy', energy)
    call input('e0', e0)
    do i = 1, 10
        write (name, '(a,i0)') 't', i
        call input(name, tjet(i))
    end do
    do i = 1, 6
        write (name, '(a,i0)') 'gh', i
        call input(name, gradient(i))
        do j = 1, 6
            write (name, '(a,i0,a,i0)') 'hh', i, '_', j
            call input(name, hh(i, j))
        end do
    end do
    do j = 1, 6
        do i = 1, 6
            kj(i, j) = tjet(1)*hh(i, j)
            if (i <= 3) kj(i, j) = kj(i, j) + tjet(1 + i)*gradient(j)
            if (j <= 3) kj(i, j) = kj(i, j) + gradient(i)*tjet(1 + j)
           if (i <= 3 .and. j <= 3) kj(i, j) = kj(i, j) + (energy - e0)*tjet(slot(i, j))
        end do
    end do
    do j = 1, 6
        do i = 1, 6
            k = 6*(j - 1) + i
            ii = mod(i + 2, 6) + 1
            roots(k) = kj(ii, j)
            if (i > 3) roots(k) = -roots(k)
            write (name, '(a,i0,a,i0)') 'df', i, '_', j
            outputs(k) = str(trim(name))
        end do
    end do
    call leaf('fo_factor_flow', roots, args(:n), outputs)
    call suite_end(suite, trim(out)//'/fo_hessian_proof.json')
contains
    subroutine invert(a, inverse)
        type(expr_t), intent(in)::a(3, 3)
        type(expr_t), intent(out)::inverse(3, 3)
        integer::i, j, ii, jj
        type(expr_t)::det
        det = a(1, 1)*a(2, 2)*a(3, 3) + 2*a(1, 2)*a(1, 3)*a(2, 3) &
              - a(1, 1)*a(2, 3)**2 - a(2, 2)*a(1, 3)**2 - a(3, 3)*a(1, 2)**2
        do j = 1, 3
            do i = 1, 3
                ii = mod(i, 3) + 1
                jj = mod(j, 3) + 1
                inverse(j, i) = (a(ii, jj)*a(mod(ii, 3) + 1, mod(jj, 3) + 1) &
                                 - a(ii, mod(jj, 3) + 1)*a(mod(ii, 3) + 1, jj))/det
            end do
        end do
    end subroutine invert
    subroutine blocks(a, g, iv, v, b, h)
        type(expr_t), intent(in)::a(3), g(3, 3), iv(3, 3), v(3), b(3, 3)
        type(expr_t), intent(out)::h(6, 6)
        integer::i, j, k, l
        type(expr_t)::t
        do j = 1, 3
            do i = 1, 3
                t = num(arena, 0)
                do k = 1, 3
                    do l = 1, 3
                        t = t + b(k, i)*iv(k, l)*b(l, j) &
                            - v(k)*diff(diff(g(k, l), q(i)), q(j))*v(l)/2
                    end do
                    t = t - charge*diff(diff(a(k), q(i)), q(j))*v(k)
                end do
                h(i, j) = t/mass
                t = num(arena, 0)
                do k = 1, 3
                    t = t - iv(i, k)*b(k, j)/mass
                end do
                h(3 + i, j) = t
                h(j, 3 + i) = t
                h(3 + i, 3 + j) = iv(i, j)/mass
            end do
        end do
    end subroutine blocks
    subroutine input(name, x)
        character(*), intent(in)::name
        type(expr_t), intent(out)::x
        n = n + 1
        args(n) = str(trim(name))
        x = sym(arena, trim(name))
    end subroutine input
    subroutine leaf(name, values, arguments, outputs)
        character(*), intent(in)::name
        type(expr_t), intent(in)::values(:)
        type(str_t), intent(in)::arguments(:), outputs(:)
        type(cse_result_t) :: cse
        type(expr_t), allocatable :: symbols(:), expanded(:)
        type(expr_t) :: parsed
        character(:), allocatable :: message
        character(:), allocatable :: path, code
        logical :: ok
        integer :: i
        spec%name = str(name)
        spec%mode = KERNEL_SUBROUTINE
        spec%temp_prefix = str('work')
        spec%generator = str('gen/gen_fo_hessian.f90')
        spec%regenerate_command = str('bash gen/regenerate_fo.sh')
        spec%args = arguments
        spec%outputs = outputs
        path = trim(out)//'/'//name//'.inc'
        open (newunit=u, file=path, status='replace')
        code = chars(emit_kernel(values, spec))
        if (code(len(code):) == new_line('a')) code = code(:len(code) - 1)
        write (u, '(a)') code
        close (u)
        cse = cse_analyse(values, 'work')
        allocate (symbols(cse%n), expanded(cse%n))
        do i = 1, cse%n
            symbols(i) = sym(arena, chars(cse%names(i)))
            parsed = parse_fortran_expr(arena, path, chars(cse%names(i)), ok, message)
            if (.not. ok) error stop 'factored temporary readback failed'
            if (i > 1) parsed = subs_many(parsed, symbols(:i - 1), expanded(:i - 1))
            expanded(i) = parsed
        end do
        do i = 1, size(outputs)
            parsed = parse_fortran_expr(arena, path, chars(outputs(i)), ok, message)
            if (.not. ok) error stop 'factored output readback failed'
            parsed = subs_many(parsed, symbols, expanded)
      call check_identity(suite, eng, 'emitted '//chars(outputs(i)), parsed - values(i))
        end do
    end subroutine leaf
end program gen_fo_hessian
