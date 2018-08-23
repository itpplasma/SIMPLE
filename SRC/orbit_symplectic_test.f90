program orbit_symplectic_test

use orbit_symplectic, only: orbit_sympl_init, orbit_timestep_sympl, f_sympl_euler, yold,&
  rmumag, ro0, f, df, d2f !, eval_field

implicit none

integer, parameter :: n = 6

double precision :: dt0
double precision, dimension(6) :: z, fvec
integer :: info

integer :: k

rmumag = 0.1d0
ro0 = 1d0
dt0 = 1.0d-1*dsqrt(2d0)

yold = 0d0
yold(1) = 0.3d0
yold(2) = 1.5d0
yold(4) = 0.1d0

call eval_field(yold(1), yold(2), yold(3), 0, f, df, d2f)

print *, f%Ath
print *, f%Aph
print *, f%Bth
print *, f%Bph
print *, f%Bmod
print *, df%dAth
print *, df%dAph
print *, df%dBth
print *, df%dBph
print *, df%dBmod

!stop

yold(5) = yold(4)*f%Bth/f%Bmod + f%Ath/ro0
yold(6) = yold(4)*f%Bph/f%Bmod + f%Aph/ro0

! call f_sympl_euler(n,yold,fvec,info)
! print *, yold
! print *, fvec

z = 0d0
z(1:3) = yold(1:3)
z(4) = yold(4)**2/2d0 + rmumag*f%Bmod
write(4001,*) z

do k = 1, 10000
  call orbit_timestep_sympl(z, dt0, dt0, info)
  write(4001,*) z
end do

end program orbit_symplectic_test
