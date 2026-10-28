program generate_event_energy
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(*), operator(**)
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_string, only: str, chars
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: energy, p, roots(1)
    type(kernel_spec_t) :: spec
    character(4096) :: output_dir
    integer :: unit

    output_dir = 'src/generated'
    if (command_argument_count() > 0) call get_command_argument(1, output_dir)
    call arena%init()
    energy = sym(arena, 'energy_eV')
    p = sym(arena, 'p')
    roots = [energy*p**2]
    spec%name = str('event_energy')
    spec%mode = KERNEL_SUBROUTINE
    spec%generator = str('gen/generate_event_energy.f90')
    spec%args = [str('energy_eV'), str('p')]
    spec%outputs = [str('kinetic_energy_eV')]
    open (newunit=unit, file=trim(output_dir)//'/generated_event_energy.f90', &
        status='replace', action='write')
    write (unit, '(a)') 'module generated_event_energy'
    write (unit, '(a)') 'implicit none'
    write (unit, '(a)') 'contains'
    write (unit, '(a)') chars(emit_kernel(roots, spec))
    write (unit, '(a)') 'end module generated_event_energy'
    close (unit)
end program generate_event_energy
