module callback

implicit none

logical :: output_error = .False.
logical :: output_orbits_macrostep = .False.

public :: set_output_orbits_macrostep
public :: get_output_orbits_macrostep

contains

subroutine set_output_orbits_macrostep(flag)
    logical, intent(in) :: flag
    output_orbits_macrostep = flag
end subroutine set_output_orbits_macrostep

logical function get_output_orbits_macrostep()
    get_output_orbits_macrostep = output_orbits_macrostep
end function get_output_orbits_macrostep

end module callback
