module field_gc_tables_input
    !> Guiding-centre tables from a field file: a solved frozen-current field
    !> (.h5, libneo frozen_current.x field_file) or a coil file. The
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

    public :: gc_tables_apply, build_gc_tables_from_file, build_gc_tables_frozen

    type, extends(gc_field_source_t) :: frozen_source_t
        type(frozen_current_field_t) :: field
    contains
        procedure :: field_jacobian => frozen_field_jacobian
    end type frozen_source_t

contains

    !> True when field_input is a coil file or a frozen-current field and the
    !> reference chart is a chartmap with a radius that is not normalised flux,
    !> the cases build_gc_tables accepts.
    logical function gc_tables_apply(path, coords)
        use field_file_kind, only: field_file_kind_of, FIELD_FILE_FROZEN_CURRENT, &
                                   FIELD_FILE_COILS
        use libneo_coordinates, only: chartmap_coordinate_system_t, PSI_TOR_NORM, &
                                      PSI_POL_NORM
        character(*), intent(in) :: path
        class(coordinate_system_t), intent(in) :: coords

        gc_tables_apply = .false.
        select type (coords)
        class is (chartmap_coordinate_system_t)
            if (coords%rho_convention == PSI_TOR_NORM .or. &
                coords%rho_convention == PSI_POL_NORM) return
        class default
            return
        end select
        select case (field_file_kind_of(path))
        case (FIELD_FILE_FROZEN_CURRENT, FIELD_FILE_COILS)
            gc_tables_apply = .true.
        end select
    end function gc_tables_apply

    subroutine build_gc_tables_from_file(path, coords, tab)
        use field_file_kind, only: field_file_kind_of, FIELD_FILE_FROZEN_CURRENT
        character(*), intent(in) :: path
        class(coordinate_system_t), intent(in) :: coords
        type(gc_tables_t), intent(inout) :: tab
        type(coils_field_t) :: raw_coils
        type(frozen_current_field_t) :: frozen
        integer :: status

        if (field_file_kind_of(path) == FIELD_FILE_FROZEN_CURRENT) then
            call load_frozen_current_field(path, frozen, status)
            if (status /= 0) error stop 'build_gc_tables_from_file: cannot load ' &
                //path
            call build_gc_tables_frozen(frozen, coords, tab)
            return
        end if
        call create_coils_field(path, raw_coils)
        call build_gc_tables_coils(raw_coils%coils, coords, tab)
    end subroutine build_gc_tables_from_file

    !> Tables of a loaded frozen-current field (libneo, cm and G).
    subroutine build_gc_tables_frozen(field, coords, tab, n_r, n_th, n_phi)
        type(frozen_current_field_t), intent(in) :: field
        class(coordinate_system_t), intent(in) :: coords
        type(gc_tables_t), intent(inout) :: tab
        integer, intent(in), optional :: n_r, n_th, n_phi

        call build_gc_tables(frozen_source_t(field), coords, tab, n_r, n_th, n_phi)
    end subroutine build_gc_tables_frozen

    subroutine frozen_field_jacobian(self, x, b, db)
        class(frozen_source_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3), db(3, 3)

        call self%field%compute_field_jacobian(x, b, db)
    end subroutine frozen_field_jacobian

end module field_gc_tables_input
