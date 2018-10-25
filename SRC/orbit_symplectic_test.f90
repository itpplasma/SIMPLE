program orbit_symplectic_test

use field_can_mod, only: eval_field, neval

use orbit_symplectic, only: orbit_sympl_init, orbit_timestep_sympl, &
  mu, ro0, f, df, d2f, z, pth

implicit none

integer, parameter :: n = 6

double precision :: dt0, vpar
double precision, dimension(6) :: z0, fvec
integer :: info

integer :: k

integer :: nts = 1000
integer :: kwrite = 1

mu = 0.1d0
ro0 = 1d0
dt0 = 0.5*0.13*dsqrt(2d0)

z(1) = 0.3d0
z(2) = 1.5d0
z(3) = 0.0d0

vpar = 0.1d0
call eval_field(z(1), z(2), z(3), 0)
z(4) = vpar*f%hph + f%Aph/ro0
pth = vpar*f%hth + f%Ath/ro0

z0 = 0d0
z0(1:3) = z(1:3)
z0(4) = vpar**2/2d0 + mu*f%Bmod
z0(5) = vpar/sqrt(vpar**2/2d0 + mu*f%Bmod)
write(4001,*) z0

do k = 1, nts
  call orbit_timestep_sympl(z0, dt0, dt0, info)
  if (mod(k, kwrite) == 0) write(4001,*) z0
end do

print *,'done. Evaluations: ', neval

end program orbit_symplectic_test
