module reference_coordinates

  use, intrinsic :: iso_fortran_env, only : dp => real64
  use libneo_coordinates, only : coordinate_system_t, make_vmec_coordinate_system

  implicit none

  class(coordinate_system_t), allocatable, public :: ref_coords

contains

  subroutine init_reference_coordinates(coord_input)
    character(*), intent(in) :: coord_input

    if (len_trim(coord_input) == 0) then
      print *, 'reference_coordinates.init_reference_coordinates: ', &
          'coord_input must be set (see params.apply_config_aliases)'
      error stop
    end if

    if (allocated(ref_coords)) deallocate(ref_coords)

    ! For now we always use VMEC reference coordinates. The params module
    ! is responsible for resolving coord_input versus legacy netcdffile
    ! and field_input; here we only rely on the final coord_input value.
    call make_vmec_coordinate_system(ref_coords)
  end subroutine init_reference_coordinates

end module reference_coordinates
