program test_coordinates

use, intrinsic :: iso_fortran_env, only: dp => real64
use simple_coordinates, only: transform_i, get_transform

implicit none

procedure(transform_i), pointer :: transform

character(32) :: from, to
real(dp) :: xfrom(3), xto(3), dxto_dxfrom(3,3)

xfrom = [0.d0, 0.d0, 0.d0]

call get_command_argument(1, from)
call get_command_argument(2, to)

transform => get_transform(from, to)
call transform(xfrom, xto, dxto_dxfrom)

end program test_coordinates
