! Shared leaf-writing support for SIMPLE FortSym generators.
! One emit path and one FortSym revision source: SIMPLE_FORTSYM_REVISION in
! CMakeLists.txt, delivered as the configured include file
! simple_fortsym_revision.inc.
module simple_gen_support
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, operator(**)
    use fortsym_string, only: str_t, str, chars
    use fortsym_kernel, only: kernel_spec_t, kernel_binding_t, emit_kernel
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_subs, only: subs_many
    implicit none
    private
    include 'simple_fortsym_revision.inc'
    public :: fortsym_revision, leaf_writer_t, binding, output_t, out_scalar, &
              out_array

    !> One declared output argument: a scalar, or an array whose elements
    !> are filled from consecutive roots in column-major order.
    type :: output_t
        character(:), allocatable :: name
        integer :: rows = 0
        integer :: cols = 0
    end type output_t

    !> Writes generated leaves into one file. The generator owns the arena;
    !> the writer owns the output unit and the leaf conventions.
    type :: leaf_writer_t
        type(arena_t), pointer :: arena => null()
        type(native_engine_t) :: engine
        character(:), allocatable :: generator
        integer :: unit = -1
        !> A symbol known to satisfy sign**2 = 1; its powers are reduced
        !> after simplification.
        type(expr_t) :: sign
        logical :: has_sign = .false.
    contains
        procedure :: open => writer_open
        procedure :: line => writer_line
        procedure :: close => writer_close
        procedure :: emit => writer_emit
        procedure :: simplify => writer_simplify
    end type leaf_writer_t

contains

    subroutine writer_open(w, arena, generator, path)
        class(leaf_writer_t), intent(inout) :: w
        type(arena_t), intent(inout), target :: arena
        character(*), intent(in) :: generator, path

        w%arena => arena
        w%engine = make_native_engine(arena)
        w%generator = generator
        open (newunit=w%unit, file=path, status='replace', action='write')
    end subroutine writer_open

    subroutine writer_line(w, text)
        class(leaf_writer_t), intent(in) :: w
        character(*), intent(in) :: text

        write (w%unit, '(a)') text
    end subroutine writer_line

    subroutine writer_close(w)
        class(leaf_writer_t), intent(inout) :: w

        close (w%unit)
        w%unit = -1
    end subroutine writer_close

    !> Simplify a root with the native engine and, for a declared sign
    !> symbol, reduce its powers with sign**2 = 1. Fails closed.
    function writer_simplify(w, e) result(r)
        class(leaf_writer_t), intent(inout) :: w
        type(expr_t), intent(in) :: e
        type(expr_t) :: r
        type(engine_result_t) :: res

        res = w%engine%simplify(e)
        if (.not. res%ok) error stop chars(res%message)
        r = res%value
        if (.not. w%has_sign) return
        r = subs_many(r, [w%sign**2, w%sign**3, w%sign**4, w%sign**5], &
                      [w%sign**0, w%sign, w%sign**0, w%sign])
        res = w%engine%simplify(r)
        if (.not. res%ok) error stop chars(res%message)
        r = res%value
    end function writer_simplify

    !> A binding of an opaque applied function (or one of its partial
    !> derivatives, by argument index) to a kernel argument expression.
    function binding(head, indices, replacement) result(b)
        character(*), intent(in) :: head, replacement
        integer, intent(in) :: indices(:)
        type(kernel_binding_t) :: b

        b%head = str(head)
        allocate (b%derivative_indices(size(indices)))
        b%derivative_indices = indices
        b%replacement = str(replacement)
    end function binding

    function out_scalar(name) result(o)
        character(*), intent(in) :: name
        type(output_t) :: o

        o%name = name
    end function out_scalar

    function out_array(name, rows, cols) result(o)
        character(*), intent(in) :: name
        integer, intent(in) :: rows
        integer, intent(in), optional :: cols
        type(output_t) :: o

        o%name = name
        o%rows = rows
        o%cols = 1
        if (present(cols)) o%cols = cols
    end function out_array

    !> Emit one pure leaf subroutine without a module wrapper, so the
    !> consumer includes it next to its hot caller. args may carry a shape
    !> suffix such as 'q(10)'; outputs are scalars or column-major arrays.
    subroutine writer_emit(w, name, args, outputs, roots, bindings)
        class(leaf_writer_t), intent(inout) :: w
        character(*), intent(in) :: name, args(:)
        type(output_t), intent(in) :: outputs(:)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_binding_t), intent(in), optional :: bindings(:)
        type(kernel_spec_t) :: spec
        type(str_t) :: source
        type(expr_t) :: clean(size(roots))
        logical :: ok, shaped
        character(:), allocatable :: message, arg
        character(32) :: ref
        integer :: j, k, i, c, n, paren

        spec%name = str(name)
        spec%generator = str(w%generator)
        spec%generator_revision = str(fortsym_revision)
        spec%regenerate_command = str('cmake --build <build> --target '// &
                                      'regenerate_'//w%generator)
        spec%pure_procedure = .true.
        shaped = any(outputs%rows > 0)
        do j = 1, size(args)
            if (index(args(j), '(') > 0) shaped = .true.
        end do
        spec%elemental_procedure = .not. shaped
        allocate (spec%args(size(args)), spec%arg_shapes(size(args)))
        do j = 1, size(args)
            arg = trim(args(j))
            paren = index(arg, '(')
            if (paren > 0) then
                spec%args(j) = str(arg(:paren - 1))
                spec%arg_shapes(j) = str(arg(paren:))
            else
                spec%args(j) = str(arg)
                spec%arg_shapes(j) = str('')
            end if
        end do
        allocate (spec%outputs(size(outputs)), spec%output_shapes(size(outputs)))
        n = 0
        do j = 1, size(outputs)
            spec%outputs(j) = str(outputs(j)%name)
            if (outputs(j)%rows > 0) then
                if (outputs(j)%cols > 1) then
                    write (ref, '(a,i0,a,i0,a)') '(', outputs(j)%rows, ',', &
                        outputs(j)%cols, ')'
                else
                    write (ref, '(a,i0,a)') '(', outputs(j)%rows, ')'
                end if
                spec%output_shapes(j) = str(trim(ref))
                n = n + outputs(j)%rows*outputs(j)%cols
            else
                spec%output_shapes(j) = str('')
                n = n + 1
            end if
        end do
        if (n /= size(roots)) error stop 'writer_emit: root count mismatch'
        allocate (spec%output_references(n))
        k = 0
        do j = 1, size(outputs)
            if (outputs(j)%rows == 0) then
                k = k + 1
                spec%output_references(k) = str(outputs(j)%name)
                cycle
            end if
            do c = 1, outputs(j)%cols
                do i = 1, outputs(j)%rows
                    k = k + 1
                    if (outputs(j)%cols > 1) then
                        write (ref, '(a,i0,a,i0,a)') '(', i, ',', c, ')'
                    else
                        write (ref, '(a,i0,a)') '(', i, ')'
                    end if
                    spec%output_references(k) = str(outputs(j)%name//trim(ref))
                end do
            end do
        end do
        if (present(bindings)) spec%bindings = bindings
        do j = 1, size(roots)
            clean(j) = w%simplify(roots(j))
        end do
        source = emit_kernel(clean, spec, ok, message=message)
        if (.not. ok) error stop message
        write (w%unit, '(a)') chars(source)
    end subroutine writer_emit

end module simple_gen_support
