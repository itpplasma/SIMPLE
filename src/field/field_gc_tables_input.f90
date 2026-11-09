module field_gc_tables_input
    !> Guiding-centre tables from a field file: a solved frozen-current field
    !> (HDF5, libneo frozen_current.x field_file) or a coil file. The
    !> frozen-current types stay out of field_gc_tables: their modules bring
    !> IEEE-module symbols into every module that can see them, and gfortran
    !> then saves and restores the floating-point environment on each call of
    !> every procedure there, which tripled the cost of a tracing step.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_coils, only: coils_field_t, create_coils_field
    use field_gc_tables, only: gc_tables_t, gc_field_source_t, build_gc_tables, &
                               build_gc_tables_coils
    use libneo_coordinates, only: coordinate_system_t
    use neo_frozen_current_field, only: frozen_current_field_t
    use neo_frozen_current_io, only: load_frozen_current_field

    implicit none
    private

    public :: build_gc_tables_from_file, build_gc_tables_frozen

    type, extends(gc_field_source_t) :: frozen_source_t
        type(frozen_current_field_t) :: field
    contains
        procedure :: field_jacobian => frozen_field_jacobian
    end type frozen_source_t

contains

    subroutine build_gc_tables_from_file(path, coords, order, tab)
        use field_file_kind, only: field_file_kind_of, FIELD_FILE_FROZEN_CURRENT, &
                                   FIELD_FILE_COILS
        character(*), intent(in) :: path
        class(coordinate_system_t), intent(in) :: coords
        integer, intent(in) :: order
        type(gc_tables_t), intent(inout) :: tab
        type(coils_field_t) :: raw_coils
        type(frozen_current_field_t) :: frozen
        integer :: status

        select case (field_file_kind_of(path))
        case (FIELD_FILE_FROZEN_CURRENT)
            call load_frozen_current_field(path, frozen, status)
            if (status /= 0) error stop 'build_gc_tables_from_file: cannot load ' &
                //path
            call build_gc_tables_frozen(frozen, coords, order, tab)
        case (FIELD_FILE_COILS)
            call create_coils_field(path, raw_coils)
            call build_gc_tables_coils(raw_coils%coils, coords, order, tab)
        case default
            error stop 'build_gc_tables_from_file: '//trim(path)// &
                ' is neither a frozen-current field nor a coil file'
        end select
    end subroutine build_gc_tables_from_file

    !> Tables of a loaded frozen-current field (libneo, cm and G).
    subroutine build_gc_tables_frozen(field, coords, order, tab, n_r, n_th, n_phi)
        type(frozen_current_field_t), intent(in) :: field
        class(coordinate_system_t), intent(in) :: coords
        integer, intent(in) :: order
        type(gc_tables_t), intent(inout) :: tab
        integer, intent(in), optional :: n_r, n_th, n_phi

        call build_gc_tables(frozen_source_t(field), coords, order, tab, n_r, n_th, &
                             n_phi)
    end subroutine build_gc_tables_frozen

    subroutine frozen_field_jacobian(self, x, b, db)
        class(frozen_source_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3), db(3, 3)

        call self%field%compute_field_jacobian(x, b, db)
    end subroutine frozen_field_jacobian

end module field_gc_tables_input
