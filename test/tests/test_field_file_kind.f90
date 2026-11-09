program test_field_file_kind
    ! Field files are told apart by content: a SPECTRE equilibrium and a
    ! frozen-current field (both HDF5) under swapped names, a VMEC wout, a coil
    ! file without the coils name. Last, an HDF5 file of neither kind must stop
    ! field_from_file with its root members named (PASS_REGULAR_EXPRESSION).
    ! usage: test_field_file_kind.x SPECTRE_H5 WOUT_NC
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field, only: field_from_file, is_spectre_file
    use field_base, only: magnetic_field_t
    use field_file_kind, only: field_file_kind_of, FIELD_FILE_FROZEN_CURRENT, &
                               FIELD_FILE_SPECTRE, FIELD_FILE_HDF5, &
                               FIELD_FILE_NETCDF, FIELD_FILE_COILS, FIELD_FILE_UNKNOWN
    use hdf5_tools, only: hid_t, h5_init, h5_deinit, h5_create, h5_close, h5_add, &
                          h5_define_group, h5_close_group
    use neo_biotsavart, only: coils_init
    use neo_frozen_current_field, only: frozen_current_field_t, frozen_current_patch_t
    use fortfem_meridian_bspline_n, only: meridian_bspline_build
    use neo_frozen_current_io, only: save_frozen_current_field
    implicit none

    character(len=4096) :: spectre_file, wout_file
    class(magnetic_field_t), allocatable :: loaded
    integer :: kind_wout

    call get_command_argument(1, spectre_file)
    call get_command_argument(2, wout_file)
    call copy_file(trim(spectre_file), 'spectre_named.field')
    call write_frozen('frozen_named.nc')
    call write_coils('field_points.txt')
    call write_other('other.h5')
    open (newunit=kind_wout, file='empty.txt', status='replace')
    write (kind_wout, '(a)') 'no field here'
    close (kind_wout)

    call expect('spectre_named.field', FIELD_FILE_SPECTRE)
    call expect('frozen_named.nc', FIELD_FILE_FROZEN_CURRENT)
    call expect('field_points.txt', FIELD_FILE_COILS)
    call expect('other.h5', FIELD_FILE_HDF5)
    call expect('empty.txt', FIELD_FILE_UNKNOWN)
    kind_wout = field_file_kind_of(trim(wout_file))
    if (kind_wout /= FIELD_FILE_NETCDF .and. kind_wout /= FIELD_FILE_HDF5) &
        error stop 'wout not recognised as NetCDF'
    if (.not. is_spectre_file('spectre_named.field')) error stop 'SPECTRE by content'
    if (is_spectre_file('frozen_named.nc')) error stop 'frozen-current is not SPECTRE'
    if (is_spectre_file('other.h5')) error stop 'other HDF5 is not SPECTRE'
    print *, 'classification ok'

    call field_from_file('other.h5', loaded)
    error stop 'an unknown HDF5 file must stop field_from_file'

contains

    subroutine copy_file(source, dest)
        character(*), intent(in) :: source, dest
        character(len=:), allocatable :: bytes
        integer :: unit, n

        open (newunit=unit, file=source, access='stream', form='unformatted', &
              status='old', action='read')
        inquire (unit=unit, size=n)
        allocate (character(len=n) :: bytes)
        read (unit) bytes
        close (unit)
        open (newunit=unit, file=dest, access='stream', form='unformatted', &
              status='replace', action='write')
        write (unit) bytes
        close (unit)
    end subroutine copy_file

    subroutine expect(path, kind)
        character(*), intent(in) :: path
        integer, intent(in) :: kind

        if (field_file_kind_of(path) /= kind) then
            print *, path, ': kind ', field_file_kind_of(path), ', expected ', kind
            error stop 'wrong field file kind'
        end if
    end subroutine expect

    subroutine write_frozen(path)
        character(*), intent(in) :: path
        type(frozen_current_field_t) :: f
        type(frozen_current_patch_t) :: patch
        integer :: status

        call coils_init([0.0_dp, 0.0_dp], [0.0_dp, 0.0_dp], [-1.0_dp, 1.0_dp], &
                        [1.0_dp, 0.0_dp], f%coils)
        patch%r0 = 1
        patch%r1 = 2
        patch%z0 = -1
        patch%z1 = 1
        patch%nr = 4
        patch%nz = 4
        patch%kmax = 0
        call meridian_bspline_build(patch%r0, patch%r1, patch%z0, patch%z1, patch%nr, &
                                    patch%nz, patch%degree, f%plasma%sp, status)
        if (status /= 0) error stop 'meridian spline build'
        allocate (f%plasma%u(f%plasma%sp%n1, 0:0))
        f%plasma%u = 0
        call save_frozen_current_field(path, patch, f)
    end subroutine write_frozen

    subroutine write_coils(path)
        character(*), intent(in) :: path
        integer :: unit

        open (newunit=unit, file=path, status='replace')
        write (unit, *) 2
        write (unit, *) 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp
        write (unit, *) 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp
        close (unit)
    end subroutine write_coils

    subroutine write_other(path)
        character(*), intent(in) :: path
        integer(hid_t) :: h5id, gid

        call h5_init()
        call h5_create(path, h5id)
        call h5_add(h5id, 'temperature', 1.0_dp)
        call h5_define_group(h5id, 'mesh', gid)
        call h5_add(gid, 'nodes', 3)
        call h5_close_group(gid)
        call h5_close(h5id)
        call h5_deinit()
    end subroutine write_other

end program test_field_file_kind
