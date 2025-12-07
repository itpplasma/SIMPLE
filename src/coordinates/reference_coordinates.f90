module reference_coordinates
    !> Reference coordinate system module.
    !>
    !> Provides polymorphic reference coordinates for SIMPLE orbit integration.
    !> Supports VMEC (default for .nc files) and chartmap (for .h5 files).
    !>
    !> Chartmap coordinates use a tabulated (rho, theta, zeta) -> (X, Y, Z)
    !> mapping stored in HDF5/NetCDF format with 3D B-spline interpolation.
    !> This enables using arbitrary reference coordinates derived from VMEC
    !> or other equilibrium codes.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, &
        make_vmec_coordinate_system, make_chartmap_coordinate_system

    implicit none

    class(coordinate_system_t), allocatable, public :: ref_coords

    integer, parameter, public :: REF_COORDS_VMEC = 1
    integer, parameter, public :: REF_COORDS_CHARTMAP = 2

    integer, public :: ref_coords_type = REF_COORDS_VMEC

contains

    subroutine init_reference_coordinates(coord_input)
        !> Initialize reference coordinate system from input file.
        !>
        !> For VMEC (.nc): Uses libneo VMEC splines (requires prior init_vmec)
        !> For chartmap (.h5): Loads tabulated coordinates with spline interpolation
        character(*), intent(in) :: coord_input

        if (len_trim(coord_input) == 0) then
            print *, 'reference_coordinates.init_reference_coordinates: ', &
                'coord_input must be set (see params.apply_config_aliases)'
            error stop
        end if

        if (allocated(ref_coords)) deallocate(ref_coords)

        if (is_chartmap_file(coord_input)) then
            call make_chartmap_coordinate_system(ref_coords, &
                trim(get_chartmap_path(coord_input)))
            ref_coords_type = REF_COORDS_CHARTMAP
            print *, 'Reference coordinates: chartmap from ', &
                trim(get_chartmap_path(coord_input))
        else
            call make_vmec_coordinate_system(ref_coords)
            ref_coords_type = REF_COORDS_VMEC
            print *, 'Reference coordinates: VMEC'
        end if
    end subroutine init_reference_coordinates

    pure function is_chartmap_file(filename) result(is_chartmap)
        !> Check if file is a chartmap file based on prefix or extension.
        !> Chartmap files can be specified as:
        !>   - chartmap:<path> - explicit chartmap prefix
        !>   - *.chartmap.nc - file with .chartmap.nc extension
        !> Otherwise defaults to VMEC interpretation for .nc files.
        character(*), intent(in) :: filename
        logical :: is_chartmap

        integer :: n

        is_chartmap = .false.
        n = len_trim(filename)

        ! Check for explicit chartmap: prefix
        if (n >= 9 .and. filename(1:9) == 'chartmap:') then
            is_chartmap = .true.
            return
        end if

        ! Check for .chartmap.nc extension
        if (n >= 12) then
            if (filename(n-11:n) == '.chartmap.nc') then
                is_chartmap = .true.
                return
            end if
        end if
    end function is_chartmap_file

    pure function get_chartmap_path(filename) result(path)
        !> Extract actual file path from chartmap specification.
        !> Strips chartmap: prefix if present.
        character(*), intent(in) :: filename
        character(len(filename)) :: path

        integer :: n

        n = len_trim(filename)
        if (n >= 9 .and. filename(1:9) == 'chartmap:') then
            path = filename(10:n)
        else
            path = filename
        end if
    end function get_chartmap_path

end module reference_coordinates
