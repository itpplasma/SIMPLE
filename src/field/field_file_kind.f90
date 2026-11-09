module field_file_kind
    !> Kind of a field file from its content, never from its name.
    !>
    !> HDF5 files are a frozen-current field (libneo frozen_current.x), a SPECTRE
    !> equilibrium, NetCDF-4, or none of these; NetCDF classic files are NetCDF;
    !> a text file whose first record is a point count followed by x, y, z,
    !> current rows is a coil file. Use this module inside procedures only: fortio
    !> brings IEEE-module symbols, see field_frozen_current.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortio, only: fortio_file_t, fortio_status_t

    implicit none
    private

    public :: field_file_kind_of, hdf5_root_members
    public :: FIELD_FILE_UNKNOWN, FIELD_FILE_FROZEN_CURRENT, FIELD_FILE_SPECTRE, &
              FIELD_FILE_HDF5, FIELD_FILE_NETCDF, FIELD_FILE_COILS

    integer, parameter :: FIELD_FILE_UNKNOWN = 0
    integer, parameter :: FIELD_FILE_FROZEN_CURRENT = 1
    integer, parameter :: FIELD_FILE_SPECTRE = 2
    !> HDF5 that is neither of the above; NetCDF-4 or unknown.
    integer, parameter :: FIELD_FILE_HDF5 = 3
    integer, parameter :: FIELD_FILE_NETCDF = 4
    integer, parameter :: FIELD_FILE_COILS = 5

    character(len=*), parameter :: frozen_members(3) = &
                                   [character(len=11) :: 'u', 'kmax', 'coil_points']
    character(len=*), parameter :: spectre_members(2) = &
                               [character(len=19) :: 'input/physics/Nvol', 'output/Rbc']

contains

    integer function field_file_kind_of(path) result(kind)
        character(*), intent(in) :: path

        character(len=8) :: magic
        integer :: unit, ios

        kind = FIELD_FILE_UNKNOWN
        open (newunit=unit, file=path, access='stream', form='unformatted', &
              status='old', action='read', iostat=ios)
        if (ios /= 0) return
        read (unit, iostat=ios) magic
        close (unit)
        if (ios == 0) then
            if (magic == achar(137)//'HDF'//achar(13)//achar(10)//achar(26)// &
                achar(10)) then
                kind = hdf5_kind(path)
                return
            end if
            if (magic(1:3) == 'CDF') then
                kind = FIELD_FILE_NETCDF
                return
            end if
        end if
        if (is_coil_text(path)) kind = FIELD_FILE_COILS
    end function field_file_kind_of

    integer function hdf5_kind(path) result(kind)
        character(*), intent(in) :: path

        kind = FIELD_FILE_HDF5
        if (has_all(path, frozen_members)) then
            kind = FIELD_FILE_FROZEN_CURRENT
        else if (has_all(path, spectre_members)) then
            kind = FIELD_FILE_SPECTRE
        end if
    end function hdf5_kind

    logical function has_all(path, members)
        character(*), intent(in) :: path
        character(*), intent(in) :: members(:)

        type(fortio_file_t) :: file
        type(fortio_status_t) :: status
        logical :: found
        integer :: i

        has_all = .false.
        call file%open(path, status)
        if (.not. status%ok()) return
        do i = 1, size(members)
            call file%exists(trim(members(i)), found, status)
            if (.not. status%ok()) exit
            if (.not. found) exit
            if (i == size(members)) has_all = .true.
        end do
        call file%close(status)
    end function has_all

    !> Root members of an HDF5 file, comma separated, for error messages.
    function hdf5_root_members(path) result(list)
        character(*), intent(in) :: path
        character(:), allocatable :: list

        type(fortio_file_t) :: file
        type(fortio_status_t) :: status
        character(len=:), allocatable :: names(:)
        logical, allocatable :: groups(:)
        integer :: i

        list = ''
        call file%open(path, status)
        if (.not. status%ok()) then
            list = '(unreadable: '//status%message//')'
            return
        end if
        call file%list_children('', names, groups, status)
        if (status%ok()) then
            do i = 1, size(names)
                if (i > 1) list = list//', '
                list = list//trim(names(i))
                if (groups(i)) list = list//'/'
            end do
        else
            list = '(unreadable: '//status%message//')'
        end if
        call file%close(status)
    end function hdf5_root_members

    !> First record a positive point count, second record x, y, z, current.
    logical function is_coil_text(path)
        character(*), intent(in) :: path

        integer :: unit, ios, n
        real(dp) :: row(4)

        is_coil_text = .false.
        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        read (unit, *, iostat=ios) n
        if (ios == 0) then
            if (n > 0) then
                read (unit, *, iostat=ios) row
                is_coil_text = ios == 0
            end if
        end if
        close (unit)
    end function is_coil_text

end module field_file_kind
