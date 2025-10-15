module callback

use simple, only : Tracer
use field_can_mod, only : can_to_ref
use netcdf_orbit_output, only : write_orbit_step

implicit none

logical :: output_error = .False.
logical :: output_orbits_macrostep = .False.

contains

subroutine callbacks_macrostep(tracer_, ipart, itime, t, z, ierr_orbit)
    type(Tracer), intent(inout) :: tracer_
    integer, intent(in) :: ipart, itime
    double precision, intent(in) :: t, z(:)
    integer, intent(in) :: ierr_orbit

    if (output_error .and. ierr_orbit /= 0) then
        call write_error(tracer_, ipart, itime, t, z, ierr_orbit)
    end if
    if (output_orbits_macrostep) call write_position(ipart, itime, t, z)
end subroutine callbacks_macrostep


subroutine write_error(tracer_, ipart, itime, t, z, ierr_orbit)
    type(Tracer), intent(inout) :: tracer_
    integer, intent(in) :: ipart, itime
    double precision, intent(in) :: t, z(:)
    integer, intent(in) :: ierr_orbit

    !$omp critical
    if (ierr_orbit == 1) then
        print *, 'Info ', ierr_orbit , 'orbit ', ipart, ' lost in step ', itime
        print *, t, z
    else
        print *, 'Error ', ierr_orbit , 'orbit ', ipart, ' step ', itime
        print *, t, z
    endif
    !$omp end critical
end subroutine write_error


subroutine write_position(ipart, itime, t, z)
    integer, intent(in) :: ipart, itime
    double precision, intent(in) :: t, z(:)

    double precision :: xref(3)

    call can_to_ref(z(1:3), xref)

    ! Write to NetCDF file (thread-safe via critical section in write_orbit_step)
    call write_orbit_step(ipart, itime, t, xref, z(4), z(5))
end subroutine write_position

end module callback
