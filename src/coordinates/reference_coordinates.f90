module reference_coordinates

  use, intrinsic :: iso_fortran_env, only : dp => real64
  use libneo_coordinates, only : coordinate_system_t, make_vmec_coordinate_system, &
      make_chartmap_coordinate_system

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

    if (is_chartmap_file(coord_input)) then
      call make_chartmap_coordinate_system(ref_coords, trim(coord_input))
      print *, 'Reference coordinates: chartmap from ', trim(coord_input)
    else
      call make_vmec_coordinate_system(ref_coords)
      print *, 'Reference coordinates: VMEC'
    end if
  end subroutine init_reference_coordinates


  function is_chartmap_file(filename) result(is_chartmap)
    !> Auto-detect chartmap file by checking for chartmap-specific variables.
    !> Chartmap files have: x, y, z, rho, theta, zeta
    !> VMEC files have: rmnc, zmns, xm, xn
    use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_nowrite, nf90_noerr

    character(*), intent(in) :: filename
    logical :: is_chartmap

    integer :: ncid, varid_x, varid_rho, ierr_x, ierr_rho

    ! Default to false
    is_chartmap = .false.

    ! Try to open file (read-only)
    ierr_x = nf90_open(trim(filename), nf90_nowrite, ncid)
    if (ierr_x /= nf90_noerr) return

    ! Check for chartmap signature variables
    ierr_x = nf90_inq_varid(ncid, 'x', varid_x)
    ierr_rho = nf90_inq_varid(ncid, 'rho', varid_rho)

    ! Chartmap file has both x and rho variables
    if (ierr_x == nf90_noerr .and. ierr_rho == nf90_noerr) then
      is_chartmap = .true.
    end if

    ierr_x = nf90_close(ncid)
  end function is_chartmap_file

end module reference_coordinates
