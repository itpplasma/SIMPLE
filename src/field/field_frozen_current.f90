module field_frozen_current
    !> Solved frozen-current field (libneo frozen_current.x field file, HDF5, cm
    !> and G): coils plus the frozen plasma current, with A and B from the axis
    !> to beyond the wall. It is splined onto the reference coordinates like a
    !> coil field, so Meiss and Albert coordinates can be built on it.
    !>
    !> Only the splined field leaves this module. The frozen-current types bring
    !> IEEE-module symbols into every module that can see them, and gfortran then
    !> saves and restores the floating-point environment in each procedure there
    !> (see field_gc_tables_input).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_base, only: magnetic_field_t
    use field_splined, only: splined_field_t, create_splined_field
    use cartesian_coordinates, only: cartesian_coordinate_system_t
    use libneo_coordinates, only: coordinate_system_t
    use neo_frozen_current_field, only: frozen_current_field_t
    use neo_frozen_current_io, only: load_frozen_current_field
    use hdf5_tools, only: hid_t, h5_init, h5_deinit, h5_open, h5_close, h5_exists

    implicit none
    private

    public :: create_frozen_current_splined_field, is_frozen_current_file

    type, extends(magnetic_field_t) :: frozen_field_t
        type(frozen_current_field_t) :: field
    contains
        procedure :: evaluate => frozen_evaluate
    end type frozen_field_t

contains

    subroutine create_frozen_current_splined_field(path, ref_coords, splined)
        character(*), intent(in) :: path
        class(coordinate_system_t), intent(in), target :: ref_coords
        type(splined_field_t), intent(out) :: splined

        type(frozen_field_t) :: raw
        integer :: status

        call load_frozen_current_field(path, raw%field, status)
        if (status /= 0) error stop 'create_frozen_current_splined_field: '// &
            'cannot load '//trim(path)
        allocate (cartesian_coordinate_system_t :: raw%coords)
        call create_splined_field(raw, ref_coords, splined)
    end subroutine create_frozen_current_splined_field

    subroutine frozen_evaluate(self, x, Acov, hcov, Bmod, sqgBctr)
        class(frozen_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod
        real(dp), intent(out), optional :: sqgBctr(3)

        real(dp) :: bcart(3)

        call self%field%compute_abfield(x, Acov, bcart)
        Bmod = norm2(bcart)
        hcov = bcart/Bmod
        if (present(sqgBctr)) sqgBctr = bcart
    end subroutine frozen_evaluate

    !> SPECTRE equilibria are HDF5 as well; a frozen-current file is told
    !> apart by its spline coefficients and mode truncation.
    logical function is_frozen_current_file(path) result(is_frozen)
        character(*), intent(in) :: path

        integer(hid_t) :: h5id
        integer :: unit, ios
        character(len=8) :: magic

        is_frozen = .false.
        open (newunit=unit, file=path, access='stream', form='unformatted', &
              status='old', action='read', iostat=ios)
        if (ios /= 0) return
        read (unit, iostat=ios) magic
        close (unit)
        if (ios /= 0) return
        if (magic /= achar(137)//'HDF'//achar(13)//achar(10)//achar(26)// &
            achar(10)) return

        call h5_init()
        call h5_open(path, h5id)
        is_frozen = h5_exists(h5id, 'u') .and. h5_exists(h5id, 'kmax') .and. &
                    h5_exists(h5id, 'coil_points')
        call h5_close(h5id)
        call h5_deinit()
    end function is_frozen_current_file

end module field_frozen_current
