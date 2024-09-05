module callback

use params, only : Tracer
use field_can_mod, only : can_to_ref

implicit none

logical :: output_err = .True.
logical :: output_orbits_macro = .True.

contains

subroutine callbacks_macrostep(tracer_, ipart, itime, t, z, ierr_orbit)
    type(Tracer), intent(inout) :: tracer_
    integer, intent(in) :: ipart, itime
    double precision, intent(in) :: t, z(:)
    integer, intent(in) :: ierr_orbit

    if (output_err .and. ierr_orbit /= 0) then
        call write_error(tracer_, ipart, itime, t, z, ierr_orbit)
    end if
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

    double precision :: xref(3), xlab(3)

    call can_to_ref(z(1:3), xref)
    xlab(1) = xref(1)**2
    xlab(2) = xref(2)
    xlab(3) = xref(3)

    !$omp critical
    write(9000+ipart, *) t, xlab, z(4), z(5)
    !$omp end critical
end subroutine write_position

end module callback
