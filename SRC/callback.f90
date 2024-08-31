module callback

use params, only : Tracer

implicit none

logical :: output_orbits_macro = .True.

contains

subroutine callbacks_macrostep(tracer_, ipart, itime, t, z, ierr_orbit)
    type(Tracer), intent(inout) :: tracer_
    integer, intent(in) :: ipart, itime
    double precision, intent(in) :: t, z(:)
    integer, intent(in) :: ierr_orbit

    if (ierr_orbit /= 0) call write_error(tracer_, ipart, itime, t, z, ierr_orbit)
    if (output_orbits_macro) call write_position(ipart, t, z)
end subroutine callbacks_macrostep


subroutine write_error(tracer_, ipart, itime, t, z, ierr_orbit)
    type(Tracer), intent(inout) :: tracer_
    integer, intent(in) :: ipart, itime
    double precision, intent(in) :: t, z(:)
    integer, intent(in) :: ierr_orbit

    !$omp critical
    print *, 'Error ', ierr_orbit , 'orbit ', ipart, ' step ', itime
    print *, t, z
    !$omp end critical
end subroutine write_error


subroutine write_position(ipart, t, z)
    integer, intent(in) :: ipart
    double precision, intent(in) :: t, z(:)

    !$omp critical
    write(9000+ipart, *) t, z
    !$omp end critical
end subroutine write_position

end module callback
