program generate_axis_maps
    ! Native q is either normalized flux (R=sqrt(q)) or radius (R=q).
    ! No chart selection or boundary policy belongs in this algebraic kernel.
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, sin, cos, sqrt
    use fortsym_expr, only: operator(+), operator(-), operator(*), &
                            operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_string, only: str, chars, str_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none
    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: q, theta, r2, x, y, qdot, thetadot
    type(expr_t) :: radius, inverse_q, radial_factor
    type(expr_t) :: forward(2), inverse(1), pullback(2)
    character(4096) :: output_dir
    character(:), allocatable :: mode
    integer :: imode, output_unit, ios

    output_dir = 'src/generated'
    if (command_argument_count() > 0) call get_command_argument(1, output_dir)
    call arena%init()
    engine = make_symengine_engine(arena)
    q = sym(arena, 'q')
    theta = sym(arena, 'theta')
    r2 = sym(arena, 'r2')
    x = sym(arena, 'x')
    y = sym(arena, 'y')
    qdot = sym(arena, 'dq')
    thetadot = sym(arena, 'dtheta')
    r2 = x*x + y*y
    open (newunit=output_unit, file=trim(output_dir)//'/axis_map_generated.f90', &
          status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'Cannot open generated axis module'
    write (output_unit, '(a)') 'module axis_map_generated'
    write (output_unit, '(a)') 'implicit none'
    write (output_unit, '(a)') 'contains'
    do imode = 1, 2
        if (imode == 1) then
            mode = 'flux'
            radius = sqrt(q)
            inverse_q = r2
        else
            mode = 'rho'
            radius = q
            inverse_q = sqrt(r2)
        end if
        forward(1) = radius*cos(theta)
        forward(2) = radius*sin(theta)
        call emit('axis_'//mode//'_forward', forward, &
                  [str('q'), str('theta')], [str('x'), str('y')])
        inverse = [inverse_q]
        call emit('axis_'//mode//'_inverse', inverse, &
                  [str('x'), str('y')], [str('q')])
        ! qdot=(dq/dr2)*2*(x*xdot+y*ydot). Angular velocity supplies
        ! the orthogonal component. Derive radial factor from the inverse map.
        radial_factor = 1/(x*diff(inverse_q, x) + y*diff(inverse_q, y))
        pullback(1) = radial_factor*qdot*x - thetadot*y
        pullback(2) = radial_factor*qdot*y + thetadot*x
        call emit('axis_'//mode//'_velocity', pullback, &
                  [str('x'), str('y'), str('dq'), str('dtheta')], &
                  [str('dx'), str('dy')])
    end do
    write (output_unit, '(a)') 'end module axis_map_generated'
    close (output_unit)
contains
    subroutine emit(name, roots, args, outputs)
        character(*), intent(in) :: name
        type(expr_t), intent(inout) :: roots(:)
        type(str_t), intent(in) :: args(:), outputs(:)
        type(kernel_spec_t) :: spec
        type(engine_result_t) :: result
        integer :: i
        do i = 1, size(roots)
            result = engine%simplify(roots(i))
            if (.not. result%ok) error stop 'FortSym simplification failed'
            roots(i) = result%value
        end do
        spec%name = str(name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str('gen/generate_axis_maps.f90')
        spec%args = args
        spec%outputs = outputs
        write (output_unit, '(a)') chars(emit_kernel(roots, spec))
    end subroutine emit
end program generate_axis_maps
