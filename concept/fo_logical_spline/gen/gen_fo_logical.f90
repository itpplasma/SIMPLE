program gen_fo_logical
    ! Canonical full-orbit equations of motion in SPECTRE logical coordinates
    ! q = (rho, theta, zeta). One scalar Hamiltonian
    !   H(q, p) = (p - A(q)/ro)_i g^ij(q) (p - A(q)/ro)_j / 2
    ! with the metric of the cylindrical embedding X = (R cos zeta, R sin zeta, Z),
    !   g_ij = R_i R_j + Z_i Z_j + R**2 delta_i3 delta_j3,
    ! and SPECTRE's A_rho = 0 gauge. The emitted right-hand side is
    ! (dH/dp, -dH/dq). R, Z enter through their value/gradient/Hessian jets, A
    ! through value/gradient jets, so the kernel takes spline outputs directly.
    use fortsym_string, only: str, chars, str_t
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(/), operator(**), func
    use fortsym_diff, only: diff, partial_derivative
    use fortsym_subs, only: subs_many
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE, &
        cse_result_t, cse_analyse
    use fortsym_fparse, only: parse_fortran_expr
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    implicit none
    integer, parameter :: NARG = 32, NOUT = 7, NJET = 28
    integer, parameter :: si(6) = [1, 1, 1, 2, 2, 3], sj(6) = [1, 2, 3, 2, 3, 3]
    type(arena_t), target :: arena
    type(symengine_engine_t) :: eng
    type(expr_t) :: q(3), p(3), w(3), g(3, 3), gi(3, 3), det, ham, ro
    type(expr_t) :: rfun, zfun, afun(3), dr(3), dz(3), tmp
    type(expr_t) :: old(NJET), replacement(NJET), raw(NOUT), roots(NOUT)
    type(expr_t), allocatable :: temp_symbols(:), temp_values(:)
    type(str_t) :: args(NARG), outputs(NOUT)
    type(engine_result_t) :: res
    type(kernel_spec_t) :: spec
    type(cse_result_t) :: cse
    type(suite_t) :: suite
    character(1024) :: outroot
    character(32) :: name
    character(:), allocatable :: message
    character(*), parameter :: LEAF = 'src/generated/fo_logical_rhs.inc'
    logical :: ok
    integer :: i, j, ii, jj, n, nrep, u, status

    outroot = '.'
    call get_environment_variable('FO_GENERATOR_OUTPUT', outroot, status=status)
    if (status /= 0) outroot = '.'
    call arena%init()
    eng = make_symengine_engine(arena)
    call suite_begin(suite, 'FO logical-coordinate Hamiltonian identities')
    n = 0
    nrep = 0
    do i = 1, 3
        write (name, '(a,i0)') 'q', i
        q(i) = sym(arena, trim(name))
        write (name, '(a,i0)') 'p', i
        p(i) = input(trim(name))
    end do
    ro = input('ro')
    rfun = func('Rfield', q)
    zfun = func('Zfield', q)
    call jet_replacements(rfun, 'r', 10)
    call jet_replacements(zfun, 'z', 10)
    afun(1) = num(arena, 0)
    afun(2) = func('Atheta', q)
    afun(3) = func('Azeta', q)
    call jet_replacements(afun(2), 'at', 4)
    call jet_replacements(afun(3), 'az', 4)
    if (n /= NARG .or. nrep /= NJET) error stop 'input count'
    do i = 1, 3
        dr(i) = diff(rfun, q(i))
        dz(i) = diff(zfun, q(i))
        w(i) = p(i) - afun(i)/ro
    end do
    do i = 1, 6
        tmp = dr(si(i))*dr(sj(i)) + dz(si(i))*dz(sj(i))
        if (si(i) == 3 .and. sj(i) == 3) tmp = tmp + rfun**2
        g(si(i), sj(i)) = tmp
        g(sj(i), si(i)) = tmp
    end do
    ! Adjugate/determinant inverse of the symmetric covariant metric.
    det = g(1, 1)*g(2, 2)*g(3, 3) + 2*g(1, 2)*g(1, 3)*g(2, 3) &
        - g(1, 1)*g(2, 3)**2 - g(2, 2)*g(1, 3)**2 - g(3, 3)*g(1, 2)**2
    do j = 1, 3
        do i = 1, 3
            ii = mod(i, 3) + 1
            jj = mod(j, 3) + 1
            gi(j, i) = (g(ii, jj)*g(mod(ii, 3) + 1, mod(jj, 3) + 1) &
                - g(ii, mod(jj, 3) + 1)*g(mod(ii, 3) + 1, jj))/det
        end do
    end do
    ham = num(arena, 0)
    do j = 1, 3
        do i = 1, 3
            ham = ham + w(i)*gi(i, j)*w(j)/2
        end do
    end do
    raw(1) = ham
    outputs(1) = str('energy')
    do i = 1, 3
        raw(1 + i) = diff(ham, p(i))
        raw(4 + i) = -diff(ham, q(i))
        write (name, '(a,i0)') 'f', i
        outputs(1 + i) = str(trim(name))
        write (name, '(a,i0)') 'f', 3 + i
        outputs(4 + i) = str(trim(name))
    end do
    do i = 1, NOUT
        raw(i) = subs_many(raw(i), old, replacement)
        roots(i) = raw(i)
    end do
    ! Independent identity: dq/dt is the inverse metric applied to the kinetic
    ! momentum, i.e. g_ij dq^j/dt = p_i - A_i/ro.
    do i = 1, 3
        tmp = num(arena, 0)
        do j = 1, 3
            tmp = tmp + subs_many(g(i, j), old, replacement)*roots(1 + j)
        end do
        call check_identity(suite, eng, 'velocity inverts the metric', &
            tmp - subs_many(w(i), old, replacement))
    end do
    ! Energy is the kinetic metric form of the emitted velocity.
    tmp = num(arena, 0)
    do j = 1, 3
        do i = 1, 3
            tmp = tmp + subs_many(g(i, j), old, replacement)*roots(1 + i)* &
                roots(1 + j)/2
        end do
    end do
    call check_identity(suite, eng, 'energy is g_ij qdot^i qdot^j / 2', &
        tmp - roots(1))
    spec%name = str('fo_logical_rhs_scalar')
    spec%mode = KERNEL_SUBROUTINE
    spec%pure_procedure = .true.
    spec%temp_prefix = str('work')
    spec%generator = str('gen/gen_fo_logical.f90')
    spec%regenerate_command = str('bash gen/regenerate_fo_logical.sh')
    spec%args = args
    spec%outputs = outputs
    open (newunit=u, file=trim(outroot)//'/'//LEAF, status='replace')
    write (u, '(a)') chars(emit_kernel(roots, spec))
    close (u)
    ! Read the emitted assignments back through FortSym and compare with the
    ! symbolic roots; no transcription of the generated text is trusted.
    cse = cse_analyse(roots, 'work')
    allocate (temp_symbols(cse%n), temp_values(cse%n))
    do i = 1, cse%n
        name = chars(cse%names(i))
        temp_symbols(i) = sym(arena, trim(name))
        tmp = parse_fortran_expr(arena, trim(outroot)//'/'//LEAF, trim(name), &
            ok, message)
        if (.not. ok) then
            print *, message
            error stop 'emitted temporary parse refused'
        end if
        if (i > 1) tmp = subs_many(tmp, temp_symbols(:i - 1), temp_values(:i - 1))
        temp_values(i) = tmp
    end do
    do i = 1, NOUT
        tmp = parse_fortran_expr(arena, trim(outroot)//'/'//LEAF, &
            chars(outputs(i)), ok, message)
        if (.not. ok) then
            print *, message
            error stop 'emitted output parse refused'
        end if
        tmp = subs_many(tmp, temp_symbols, temp_values)
        call check_identity(suite, eng, 'emitted '//chars(outputs(i)), &
            tmp - roots(i))
    end do
    ! Array facade, generated from the same ordered interface.
    open (newunit=u, file=trim(outroot)//'/'//LEAF, position='append', &
        status='old')
    write (u, '(a)') '', &
        'pure subroutine fo_logical_rhs(p, ro, rj, zj, atj, azj, energy, f)', &
        '    real(dp), intent(in) :: p(3), ro, rj(10), zj(10), atj(4), azj(4)', &
        '    real(dp), intent(out) :: energy, f(6)', &
        '    call fo_logical_rhs_scalar(p(1), p(2), p(3), ro, &'
    call facade_args('rj', 10)
    call facade_args('zj', 10)
    call facade_args('atj', 4)
    call facade_args('azj', 4)
    write (u, '(a)') '        energy, f(1), f(2), f(3), f(4), f(5), f(6))', &
        'end subroutine fo_logical_rhs'
    close (u)
    call suite_end(suite, trim(outroot)//'/build/fo_generator/fo_logical.json')
contains
    subroutine facade_args(label, count)
        character(*), intent(in) :: label
        integer, intent(in) :: count
        integer :: k
        do k = 1, count
            write (u, '(a,a,a,i0,a)') '        ', label, '(', k, '), &'
        end do
    end subroutine facade_args

    subroutine jet_replacements(field, prefix, slots)
        ! Slot order matches libneo's der2 output: value, d1, d2, d3, then the
        ! Hessian d11, d12, d13, d22, d23, d33.
        type(expr_t), intent(in) :: field
        character(*), intent(in) :: prefix
        integer, intent(in) :: slots
        integer, parameter :: ha(6) = [1, 1, 1, 2, 2, 3], hb(6) = [1, 2, 3, 2, 3, 3]
        integer :: slot
        character(32) :: label
        do slot = 1, slots
            nrep = nrep + 1
            if (slot == 1) then
                old(nrep) = field
            else if (slot <= 4) then
                old(nrep) = partial_derivative(field, slot - 1)
            else
                old(nrep) = partial_derivative(partial_derivative(field, &
                    ha(slot - 4)), hb(slot - 4))
            end if
            write (label, '(a,i0)') prefix, slot
            replacement(nrep) = input(trim(label))
        end do
    end subroutine jet_replacements

    function input(label) result(x)
        character(*), intent(in) :: label
        type(expr_t) :: x
        n = n + 1
        args(n) = str(label)
        x = sym(arena, label)
    end function input
end program gen_fo_logical
