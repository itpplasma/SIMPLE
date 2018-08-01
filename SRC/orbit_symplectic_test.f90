program orbit_symplectic_test

use orbit_symplectic, only: f_sympl_euler

implicit none

integer, parameter :: n = 6
double precision, parameter :: tol = 1e-10

double precision, dimension(6) :: x, fvec
integer :: info

x = 1.0

call hybrd1 (f_sympl_euler, n, x, fvec, tol, info)
print *, info
print *, fvec

end program orbit_symplectic_test
