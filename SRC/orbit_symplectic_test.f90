program orbit_symplectic_test

use field_can_mod, only: eval_field

use orbit_symplectic, only: orbit_sympl_init, orbit_timestep_sympl, f_sympl_euler, &
  mu, ro0, f, df, d2f, q, w, pth

implicit none

integer, parameter :: n = 6

double precision :: dt0, vpar
double precision, dimension(6) :: z, fvec
integer :: info

integer :: k

mu = 0.1d0
ro0 = 1d0
dt0 = 0.13*dsqrt(2d0)

w(1) = 0.3d0
q(1) = 1.5d0
q(2) = 0.0d0

vpar = 0.1d0
call eval_field(w(1), q(1), q(2), 0, f, df, d2f)
w(2) = vpar*f%Bph/f%Bmod + f%Aph/ro0
pth = vpar*f%Bth/f%Bmod + f%Ath/ro0

z = 0d0
z(1) = w(1)
z(2:3) = q
z(4) = vpar**2/2d0 + mu*f%Bmod
z(5) = vpar/sqrt(vpar**2/2d0 + mu*f%Bmod)
write(4001,*) z

do k = 1, 1000
  call orbit_timestep_sympl(z, dt0, dt0, info)
  write(4001,*) z
end do

end program orbit_symplectic_test
