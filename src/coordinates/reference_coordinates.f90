module reference_coordinates

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, make_vmec_coordinate_system, &
        make_chartmap_coordinate_system, detect_refcoords_file_type, &
        refcoords_file_chartmap, refcoords_file_vmec_wout, refcoords_file_unknown

    implicit none

    class(coordinate_system_t), allocatable, public :: ref_coords

contains

    subroutine init_reference_coordinates(coord_input)
        character(*), intent(in) :: coord_input

        integer :: file_type, ierr
        character(len=2048) :: message

        if (len_trim(coord_input) == 0) then
            print *, 'reference_coordinates.init_reference_coordinates: ', &
                'coord_input must be set (see params.apply_config_aliases)'
            error stop
        end if

        if (allocated(ref_coords)) deallocate(ref_coords)

        call detect_refcoords_file_type(coord_input, file_type, ierr, message)
        if (ierr /= 0) then
            print *, 'reference_coordinates.init_reference_coordinates: ', &
                'file detection error for ', trim(coord_input)
            print *, trim(message)
            error stop
        end if

        select case (file_type)
        case (refcoords_file_chartmap)
            call make_chartmap_coordinate_system(ref_coords, coord_input)
        case (refcoords_file_vmec_wout)
            call make_vmec_coordinate_system(ref_coords)
        case (refcoords_file_unknown)
            print *, 'reference_coordinates.init_reference_coordinates: ', &
                'unknown file type for ', trim(coord_input)
            print *, 'Expected VMEC wout (*.nc with rmnc) or chartmap (*.nc with ', &
                'rho/theta/zeta dims and x/y/z vars)'
            error stop
        case default
            print *, 'reference_coordinates.init_reference_coordinates: ', &
                'unexpected file_type ', file_type
            error stop
        end select
    end subroutine init_reference_coordinates

end module reference_coordinates
