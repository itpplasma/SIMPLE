program gen_fo_poly_axis
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), operator(*), operator(**)
    use fortsym_diff, only: diff
    use fortsym_subs, only: subs_many
    use fortsym_fparse, only: parse_fortran_expr
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_string, only: str_t, str, chars
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE, cse_result_t, cse_analyse
    implicit none
    type(arena_t), target::arena
    type(expr_t)::x, c(0:8), p, roots(3)
    type(str_t)::args(10), outputs(3)
    type(kernel_spec_t)::spec
    type(symengine_engine_t)::eng
    type(suite_t)::suite
    character(1024)::out
    character(40)::name
    character(:), allocatable::code
    integer::u, m, n, second, status, top
    call arena%init()
    eng = make_symengine_engine(arena)
    call suite_begin(suite, "Axis polynomial jet readback")
    x = sym(arena, 'x'); args(1) = str('x')
    do n = 0, 8
        write (name, '(a,i0)') 'c', n
        c(n) = sym(arena, trim(name)); args(n + 2) = str(trim(name))
    end do
    outputs = [str('r0'), str('r1'), str('r2')]
    call get_environment_variable('FO_GENERATOR_OUTPUT', out, status=status)
    if (status /= 0) error stop 'output directory required'
    do second = 1, 2
        do m = 0, 8
            top = m + 2*((8 - m)/2)
            p = c(top)
            do n = top - 2, m, -2
                p = c(n) + x**2*p
            end do
            p = x**m*p
            roots(1) = p; roots(2) = diff(p, x); roots(3) = num(arena, 0)
            if (second == 2) roots(3) = diff(roots(2), x)
            write (name, '(a,i0,a,i0)') 'fo_axis_', m, '_', second
            spec%name = str(trim(name)); spec%mode = KERNEL_SUBROUTINE
            spec%temp_prefix = str('work'); spec%generator = str('gen/gen_fo_poly_axis.f90')
            spec%regenerate_command = str('bash gen/regenerate_fo.sh')
            spec%args = args; spec%outputs = outputs
            code = chars(emit_kernel(roots, spec))
            open (newunit=u, file=trim(out)//'/'//trim(name)//'.tmp', status='replace')
            write (u, '(a)') code
            close (u)
            call readback(trim(out)//'/'//trim(name)//'.tmp', roots)
        end do
    end do
    call suite_end(suite, trim(out)//'/fo_poly_axis_proof.json')
contains
    subroutine readback(path, values)
        character(*), intent(in)::path
        type(expr_t), intent(in)::values(3)
        type(cse_result_t)::cse
        type(expr_t), allocatable::symbols(:), expanded(:)
        type(expr_t)::parsed
        character(:), allocatable::message
        logical::ok
        integer::i
        cse = cse_analyse(values, 'work')
        allocate (symbols(cse%n), expanded(cse%n))
        do i = 1, cse%n
            symbols(i) = sym(arena, chars(cse%names(i)))
            parsed = parse_fortran_expr(arena, path, chars(cse%names(i)), ok, message)
            if (.not. ok) error stop 'axis temporary readback failed'
            if (i > 1) parsed = subs_many(parsed, symbols(:i - 1), expanded(:i - 1))
            expanded(i) = parsed
        end do
        do i = 1, 3
            parsed = parse_fortran_expr(arena, path, chars(outputs(i)), ok, message)
            if (.not. ok) error stop 'axis output readback failed'
            parsed = subs_many(parsed, symbols, expanded)
   call check_identity(suite, eng, 'emitted axis radial derivative', parsed - values(i))
        end do
    end subroutine
end program
