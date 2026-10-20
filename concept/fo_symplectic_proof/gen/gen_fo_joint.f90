program gen_fo_joint
    ! One CSE graph from cylindrical embedding jets to the canonical Hessian.
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, func, operator(+), operator(-), &
                            operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs_many
    use fortsym_string, only: str_t, str, chars
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE, &
                              cse_result_t, cse_analyse
    use fortsym_fparse, only: parse_fortran_expr
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    implicit none
    type(arena_t), target :: arena
    type(symengine_engine_t) :: eng
    type(suite_t) :: suite
    type(expr_t) :: q(3), rz(2), g(3, 3), gj(6, 10), iv(3, 3), aj(3, 10)
    type(expr_t) :: p(3), w(3), v(3), b(3, 3), gh(6), hh(6, 6), energy, term, det
    type(expr_t) :: before(40), after(40), roots(103), tj(10), e0, raw
    type(str_t) :: args(100), names(103)
    character(1024) :: out
    character(40) :: name
    integer :: i, j, k, l, c, r, t, z, n, nb, ns, status, mode, ii, jj
    integer, parameter :: si(6) = [1, 1, 1, 2, 2, 3], sj(6) = [1, 2, 3, 2, 3, 3]
    integer, parameter :: slot(3, 3) = reshape([5, 6, 8, 6, 7, 9, 8, 9, 10], [3, 3])
    integer, parameter :: index(3, 10) = reshape([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, &
                         2, 0, 0, 1, 1, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 2], [3, 10])
    call arena%init()
    eng = make_symengine_engine(arena)
    call suite_begin(suite, 'Joint embedding metric and Hamiltonian readback')
    call get_environment_variable('FO_GENERATOR_OUTPUT', out, status=status)
    if (status /= 0) out = 'concept/fo_symplectic_proof/generated'
    do i = 1, 3
        write (name, '(a,i0)') 'q', i
        q(i) = sym(arena, trim(name))
    end do
    rz(1) = func('R', q); rz(2) = func('Z', q)
    n = 0; nb = 0
    do c = 1, 2
        do z = 0, 3
            do t = 0, 3 - z
                do r = 0, 3 - z - t
                    write (name, '(a,4i1)') 'j', r, t, z, c
                    nb = nb + 1
                    before(nb) = derivative(rz(c), [r, t, z])
                    call input(name, after(nb))
                end do
            end do
        end do
    end do
    do j = 1, 10
        do i = 1, 3
            write (name, '(a,i0,a,i0)') 'a', i, '_', j
            call input(name, aj(i, j))
        end do
    end do
    do i = 1, 3
        write (name, '(a,i0)') 'p', i
        call input(name, p(i))
        w(i) = p(i) - aj(i, 1)
    end do
    do j = 1, 3
        do i = 1, 3
            g(i, j) = diff(rz(1), q(i))*diff(rz(1), q(j)) &
                      + diff(rz(2), q(i))*diff(rz(2), q(j))
            if (i == 3 .and. j == 3) g(i, j) = g(i, j) + rz(1)**2
        end do
    end do
    do j = 1, 10
        do l = 1, 6
            term = derivative(g(si(l), sj(l)), index(:, j))
            gj(l, j) = subs_many(term, before, after)
        end do
    end do
    do l = 1, 6
        g(si(l), sj(l)) = gj(l, 1); g(sj(l), si(l)) = gj(l, 1)
    end do
    det = g(1, 1)*g(2, 2)*g(3, 3) + 2*g(1, 2)*g(1, 3)*g(2, 3) &
          - g(1, 1)*g(2, 3)**2 - g(2, 2)*g(1, 3)**2 - g(3, 3)*g(1, 2)**2
    do j = 1, 3
        do i = 1, 3
            ii = mod(i, 3) + 1; jj = mod(j, 3) + 1
            iv(j, i) = (g(ii, jj)*g(mod(ii, 3) + 1, mod(jj, 3) + 1) &
                        - g(ii, mod(jj, 3) + 1)*g(mod(ii, 3) + 1, jj))/det
        end do
    end do
    do i = 1, 3
        v(i) = num(arena, 0)
        do j = 1, 3
            v(i) = v(i) + iv(i, j)*w(j)
        end do
    end do
    energy = num(arena, 0)
    do i = 1, 3
        energy = energy + w(i)*v(i)/2
        gh(i) = num(arena, 0)
        gh(i + 3) = v(i)
        do k = 1, 3
            gh(i) = gh(i) - aj(k, 1 + i)*v(k)
            b(k, i) = aj(k, 1 + i)
            do j = 1, 3
                do l = 1, 6
                    if (si(l) == min(k, j) .and. sj(l) == max(k, j)) exit
                end do
                gh(i) = gh(i) - v(k)*gj(l, 1 + i)*v(j)/2
                b(k, i) = b(k, i) + gj(l, 1 + i)*v(j)
            end do
        end do
    end do
    do j = 1, 3
        do i = 1, 3
            term = num(arena, 0)
            do k = 1, 3
                do l = 1, 3
                    term = term + b(k, i)*iv(k, l)*b(l, j)
                end do
                term = term - aj(k, slot(i, j))*v(k)
            end do
            do l = 1, 6
                c = 1
                if (si(l) /= sj(l)) c = 2
                term = term - c*v(si(l))*gj(l, slot(i, j))*v(sj(l))/2
            end do
            hh(i, j) = term
            term = num(arena, 0)
            do k = 1, 3
                term = term - iv(i, k)*b(k, j)
            end do
            hh(3 + i, j) = term; hh(j, 3 + i) = term
            hh(3 + i, 3 + j) = iv(i, j)
        end do
    end do
    do mode = 1, 2
        ns = 0
        call output('energy', energy)
        do i = 1, 6
            write (name, '(a,i0)') 'gh', i
            call output(name, gh(i))
        end do
        do j = 1, merge(4, 10, mode == 1)
            do l = 1, 6
                write (name, '(a,i0,a,i0)') 'g', l, '_', j
                call output(name, gj(l, j))
            end do
        end do
        if (mode == 2) then
            do j = 1, 6
                do i = 1, 6
                    write (name, '(a,i0,a,i0)') 'hh', i, '_', j
                    call output(name, hh(i, j))
                end do
            end do
        end if
        if (mode == 1) then
            call leaf('fo_joint_first')
        else
            call leaf('fo_joint_full')
        end if
    end do
    ! The monitor follows the existing product rule, after the raw H derivatives.
    n = 0; ns = 0
    call input('energy', raw); call input('e0', e0)
    do i = 1, 10
        write (name, '(a,i0)') 't', i
        call input(name, tj(i))
    end do
    do i = 1, 6
        write (name, '(a,i0)') 'gh', i
        call input(name, gh(i))
        do j = 1, 6
            write (name, '(a,i0,a,i0)') 'hh', i, '_', j
            call input(name, hh(i, j))
        end do
    end do
    call output('kenergy', tj(1)*(raw - e0))
    do i = 1, 6
        ii = mod(i + 2, 6) + 1
        term = tj(1)*gh(ii)
        if (ii <= 3) term = term + tj(1 + ii)*(raw - e0)
        if (i > 3) term = -term
        write (name, '(a,i0)') 'f', i
        call output(name, term)
        do j = 1, 6
            term = tj(1)*hh(ii, j)
            if (ii <= 3) term = term + tj(1 + ii)*gh(j)
            if (j <= 3) term = term + gh(ii)*tj(1 + j)
            if (ii <= 3 .and. j <= 3) term = term + (raw - e0)*tj(slot(ii, j))
            if (i > 3) term = -term
            write (name, '(a,i0,a,i0)') 'df', i, '_', j
            call output(name, term)
        end do
    end do
    call leaf('fo_joint_flow')
    call suite_end(suite, trim(out)//'/fo_joint_proof.json')
contains
    function derivative(x, orders) result(y)
        type(expr_t), intent(in)::x
        integer, intent(in)::orders(3)
        type(expr_t)::y
        integer::axis, order
        y = x
        do axis = 1, 3
            do order = 1, orders(axis)
                y = diff(y, q(axis))
            end do
        end do
    end function
    subroutine input(label, x)
        character(*), intent(in)::label
        type(expr_t), intent(out)::x
        n = n + 1; args(n) = str(trim(label)); x = sym(arena, trim(label))
    end subroutine
    subroutine output(label, x)
        character(*), intent(in)::label
        type(expr_t), intent(in)::x
        ns = ns + 1; names(ns) = str(trim(label)); roots(ns) = x
    end subroutine
    subroutine leaf(label)
        character(*), intent(in)::label
        type(kernel_spec_t)::spec
        type(cse_result_t)::cse
        type(expr_t), allocatable::symbols(:), expanded(:)
        type(expr_t)::parsed
        character(:), allocatable::code, path, message
        logical::ok
        integer::u, i
        spec%name = str(label); spec%mode = KERNEL_SUBROUTINE
        spec%temp_prefix = str('work'); spec%generator = str('gen/gen_fo_joint.f90')
        spec%regenerate_command = str('bash gen/regenerate_fo.sh')
        spec%args = args(:n); spec%outputs = names(:ns)
        path = trim(out)//'/'//label//'.tmp'
        code = chars(emit_kernel(roots(:ns), spec))
        if (code(len(code):) == new_line('a')) code = code(:len(code) - 1)
        open (newunit=u, file=path, status='replace')
        write (u, '(a)') code
        close (u)
        cse = cse_analyse(roots(:ns), 'work')
        allocate (symbols(cse%n), expanded(cse%n))
        do i = 1, cse%n
            symbols(i) = sym(arena, chars(cse%names(i)))
            parsed = parse_fortran_expr(arena, path, chars(cse%names(i)), ok, message)
            if (.not. ok) error stop 'joint temporary readback failed'
            if (i > 1) parsed = subs_many(parsed, symbols(:i - 1), expanded(:i - 1))
            expanded(i) = parsed
        end do
        do i = 1, ns
            parsed = parse_fortran_expr(arena, path, chars(names(i)), ok, message)
            if (.not. ok) error stop 'joint output readback failed'
            parsed = subs_many(parsed, symbols, expanded)
         call check_identity(suite, eng, 'emitted '//chars(names(i)), parsed - roots(i))
        end do
    end subroutine
end program gen_fo_joint
